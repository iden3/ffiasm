#include <thread>
#include <vector>
#include <stdexcept>
#include <cstdint>
#include <cstring>

#include "misc.hpp"

using namespace std;

template <typename Field>
u_int32_t FFT<Field>::log2(u_int64_t n) {
    assert(n!=0);
    u_int32_t res=0;
    while (n!=1) {
        n >>= 1;
        res ++;
    }
    return res;
}

static inline u_int64_t BR(u_int64_t x, u_int64_t domainPow)
{
    x = (x >> 16) | (x << 16);
    x = ((x & 0xFF00FF00) >> 8) | ((x & 0x00FF00FF) << 8);
    x = ((x & 0xF0F0F0F0) >> 4) | ((x & 0x0F0F0F0F) << 4);
    x = ((x & 0xCCCCCCCC) >> 2) | ((x & 0x33333333) << 2);
    return (((x & 0xAAAAAAAA) >> 1) | ((x & 0x55555555) << 1)) >> (32-domainPow);
}

#define ROOT(s,j) (rootsOfUnit[(1<<(s))+(j)])

template <size_t N64>
static inline void shr_words(uint64_t (&out)[N64], const uint64_t (&in)[N64], u_int32_t shift) {
    if (shift == 0) {
        for (size_t i = 0; i < N64; i++) out[i] = in[i];
        return;
    }
    const u_int32_t wordShift = shift / 64;
    const u_int32_t bitShift  = shift % 64;

    for (size_t i = 0; i < N64; i++) out[i] = 0;

    for (size_t i = wordShift; i < N64; i++) {
        uint64_t low = in[i];
        uint64_t hi  = (i + 1 < N64) ? in[i + 1] : 0;

        if (bitShift == 0) {
            out[i - wordShift] = low;
        } else {
            out[i - wordShift] = (low >> bitShift) | (hi << (64 - bitShift));
        }
    }
}

template <size_t N64>
static inline void shr_words_inplace(uint64_t (&a)[N64], u_int32_t shift) {
    uint64_t tmp[N64];
    shr_words(tmp, a, shift);
    for (size_t i = 0; i < N64; i++) a[i] = tmp[i];
}

template <typename Field>
FFT<Field>::FFT(u_int64_t maxDomainSize, uint32_t _nThreads)
    : threadPool(ThreadPool::defaultPool())
{
    f = Field::field;

    u_int32_t domainPow = log2(maxDomainSize);

#if defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__)
    static_assert(__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__,
                  "This FFT exp-cast expects little-endian machine layout.");
#endif

    // qm1_norm = q-1 in NORMAL (non-montgomery) representation
    Element qm1_norm;
    f.fromMontgomery(qm1_norm, f.negOne());

    uint64_t qm1_words[Field::N64];
    std::memcpy(qm1_words, (const void*)qm1_norm.v, sizeof(qm1_words));

    // qm1d2 = (q-1)/2
    uint64_t qm1d2_words[Field::N64];
    shr_words(qm1d2_words, qm1_words, 1);

    // find nqr: cand^((q-1)/2) != 1
    Element cand, res;
    uint64_t cand_ui = 2;
    for (;;) {
        f.fromUI(cand, cand_ui);
        f.exp(res, cand,
              reinterpret_cast<uint8_t*>(qm1d2_words),
              (unsigned)sizeof(qm1d2_words));
        if (!f.eq(res, f.one())) {
            f.copy(nqr, cand);
            break;
        }
        cand_ui++;
    }

    // aux_words starts from (q-1)/2 and we divide by 2 until we reach domainPow
    // so at the end: aux_words = (q-1)/2^domainPow
    uint64_t aux_words[Field::N64];
    for (size_t i = 0; i < (size_t)Field::N64; i++) aux_words[i] = qm1d2_words[i];

    u_int32_t s_tmp = 1;
    while (s_tmp < domainPow) {
        if (aux_words[0] & 1ULL) break;      // odd => stop
        shr_words_inplace(aux_words, 1);
        s_tmp++;
    }

    if (s_tmp < domainPow) {
        throw std::range_error("Domain size too big for the curve");
    }

    s = s_tmp;
    uint64_t nRoots = 1ULL << s;

    roots = new Element[nRoots];
    powTwoInv = new Element[s + 1];

    f.copy(roots[0], f.one());
    f.copy(powTwoInv[0], f.one());

    if (nRoots > 1) {
        // primitive 2^s root of unity: roots[1] = nqr^(aux_words)
        f.exp(roots[1], nqr,
              reinterpret_cast<uint8_t*>(aux_words),
              (unsigned)sizeof(aux_words));

        // powTwoInv[1] = 1/2
        Element two;
        f.fromUI(two, 2);
        f.inv(powTwoInv[1], two);
    }

    threadPool.parallelBlock([&] (uint64_t nThreads, uint64_t idThread) {
        uint64_t increment = nRoots / nThreads;
        uint64_t start = idThread==0 ? 2 : idThread * increment;
        uint64_t end   = idThread==nThreads-1 ? nRoots : (idThread+1) * increment;

        if (end > start) {
            // roots[start] = roots[1]^start
            f.exp(roots[start], roots[1], (uint8_t *)(&start), sizeof(start));
        }
        for (uint64_t i = start + 1; i < end; i++) {
            f.mul(roots[i], roots[i-1], roots[1]);
        }
    });

    Element aux;
    f.mul(aux, roots[nRoots - 1], roots[1]);
    assert(f.eq(aux, f.one()));

    for (uint64_t i = 2; i <= s; i++) {
        f.mul(powTwoInv[i], powTwoInv[i-1], powTwoInv[1]);
    }
}

template <typename Field>
FFT<Field>::~FFT() {
    delete[] roots;
    delete[] powTwoInv;
}

template <typename Field>
void FFT<Field>::reversePermutation(Element *a, u_int64_t n) {
    int domainPow = log2(n);

    threadPool.parallelFor(0, n, [&] (int begin, int end, int numThread) {
        for (u_int64_t i=begin; i<end; i++) {
            Element tmp;
            u_int64_t r = BR(i, domainPow);
            if (i>r) {
                f.copy(tmp, a[i]);
                f.copy(a[i], a[r]);
                f.copy(a[r], tmp);
            }
        }
    });
}

template <typename Field>
void FFT<Field>::fft(Element *a, u_int64_t n) {
    reversePermutation(a, n);
    u_int64_t domainPow = log2(n);
    assert(((u_int64_t)1 << domainPow) == n);

    for (u_int32_t s=1; s<=domainPow; s++) {
        u_int64_t m = 1 << s;
        u_int64_t mdiv2 = m >> 1;

        threadPool.parallelFor(0, (n>>1), [&] (int begin, int end, int numThread) {
            for (u_int64_t i=begin; i< (u_int64_t)end; i++) {
                Element t;
                Element u;
                u_int64_t k=(i/mdiv2)*m;
                u_int64_t j=i%mdiv2;

                f.mul(t, root(s, j), a[k+j+mdiv2]);
                f.copy(u, a[k+j]);
                f.add(a[k+j], t, u);
                f.sub(a[k+j+mdiv2], u, t);
            }
        });
    }
}

template <typename Field>
void FFT<Field>::ifft(Element *a, u_int64_t n ) {
    fft(a, n);
    u_int64_t domainPow = log2(n);
    u_int64_t nDiv2= n >> 1;

    threadPool.parallelFor(1, nDiv2, [&] (int begin, int end, int numThread) {
        for (u_int64_t i=begin; i<(u_int64_t)end; i++) {
            Element tmp;
            u_int64_t r = n-i;
            f.copy(tmp, a[i]);
            f.mul(a[i], a[r], powTwoInv[domainPow]);
            f.mul(a[r], tmp, powTwoInv[domainPow]);
        }
    });

    f.mul(a[0], a[0], powTwoInv[domainPow]);
    f.mul(a[n >> 1], a[n >> 1], powTwoInv[domainPow]);
}

template <typename Field>
void FFT<Field>::printVector(Element *a, u_int64_t n ) {
    cout << "[" << endl;
    for (u_int64_t i=0; i<n; i++) {
        cout << f.toString(a[i]) << endl;
    }
    cout << "]" << endl;
}
