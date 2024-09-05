#include <gmp.h>
#include <iostream>

#include "gtest/gtest.h"
#include "alt_bn128.hpp"
#include "fft.hpp"

using namespace AltBn128;

namespace {

TEST(altBn128, f2_simpleMul) {

    F2Element e1;
    F2.fromString(e1, "(2,2)");

    F2Element e2;
    F2.fromString(e2, "(3,3)");

    F2Element e3;
    F2.mul(e3, e1, e2);

    F2Element e33;
    F2.fromString(e33, "(0,12)");

    // std::cout << F2.toString(e3) << std::endl;

    ASSERT_TRUE(F2.eq(e3, e33));
}

TEST(altBn128, f6_mulDiv) {
    F6Element a;
    F6.fromString(a, "((1,2),(3,4),(5,6))");

    F6Element b;
    F6.fromString(b, "((12,11),(10,9),(8,7))");

    F6Element c,d;

    F6.mul(c,a,b);
    F6.div(d,c,b);

    ASSERT_TRUE(F6.eq(a,d));
}

TEST(altBn128, f6_inv) {
    F6Element a;
    F6.fromString(a, "((239487238491,2356249827341),"
                     "(082659782,182703523765),"
                     "(978236549263,64893242))");

    F6Element inv;
    F6.inv(inv, a);

    F6Element b;
    F6.mul(b, inv, a);

    ASSERT_TRUE(F6.isOne(b));
}

TEST(altBn128, f12_mulDiv) {
    F12Element a;
    F12.fromString(a, "(((1,2),(3,4),(5,6)),((7,8),(9,10),(11,12)))");

    F12Element b;
    F12.fromString(b, "(((12,11),(10,9),(8,7)),((6,5),(4,3),(2,1)))");

    F12Element c,d;

    F12.mul(c,a,b);
    F12.div(d,c,b);

    ASSERT_TRUE(F12.eq(a,d));
}

TEST(altBn128, f12_inv) {
    F12Element a;
    F12.fromString(a,
                   "(((239846234862342323958623,2359862352529835623),"
                     "(928836523,9856234),"
                     "(235635286,5628392833)),"
                    "((252936598265329856238956532167968,23596239865236954178968),"
                     "(95421692834,236548),"
                     "(924523,12954623)))" );

    F12Element inv;
    F12.inv(inv, a);

    F12Element b;
    F12.mul(b, inv, a);

    ASSERT_TRUE(F12.isOne(b));
}

TEST(altBn128, g1_PlusZero) {
    G1Point p1;

    G1.add(p1, G1.one(), G1.zero());

    ASSERT_TRUE(G1.eq(p1, G1.one()));
}

TEST(altBn128, g1_minus_g1) {
    G1Point p1;

    G1.sub(p1, G1.one(), G1.one());

    ASSERT_TRUE(G1.isZero(p1));
}

TEST(altBn128, g1_times_4) {
    G1Point p1;
    G1.add(p1, G1.one(), G1.one());
    G1.add(p1, p1, G1.one());
    G1.add(p1, p1, G1.one());

    G1Point p2;
    G1.dbl(p2, G1.one());
    G1.dbl(p2, p2);

    ASSERT_TRUE(G1.eq(p1,p2));
}


TEST(altBn128, g1_times_3) {
    G1Point p1;
    G1.add(p1, G1.one(), G1.one());
    G1.add(p1, p1, G1.one());

    G1Point p2;
    G1.dbl(p2, G1.one());
    G1.dbl(p2, p2);
    G1.sub(p2, p2, G1.one());

    ASSERT_TRUE(G1.eq(p1,p2));
}

TEST(altBn128, g1_times_3_exp) {
    G1Point p1;
    G1.add(p1, G1.one(), G1.one());
    G1.add(p1, p1, G1.one());

    mpz_t e;
    mpz_init_set_str(e, "3", 10);

    uint8_t scalar[32];
    for (int i=0;i<32;i++) scalar[i] = 0;
    mpz_export((void *)scalar, NULL, -1, 8, -1, 0, e);
    mpz_clear(e);

    G1Point p2;
    G1.mulByScalar(p2, G1.one(), scalar, 32);

    ASSERT_TRUE(G1.eq(p1,p2));
}

TEST(altBn128, g1_times_5) {
    G1Point p1;
    G1.dbl(p1, G1.one());
    G1.dbl(p1, p1);
    G1.add(p1, p1, p1);

    G1Point p2;
    G1Point p3;
    G1Point p4;
    G1Point p5;
    G1Point p6;
    G1.dbl(p2, G1.one());
    G1.dbl(p3, p2);
    G1.dbl(p4, G1.one());
    G1.dbl(p5, p4);
    G1.add(p6, p3, p5);

    ASSERT_TRUE(G1.eq(p1,p6));
}

TEST(altBn128, g1_times_65_exp) {

    G1Point p1;
    G1.dbl(p1, G1.one());
    G1.dbl(p1, p1);
    G1.dbl(p1, p1);
    G1.dbl(p1, p1);
    G1.dbl(p1, p1);
    G1.dbl(p1, p1);
    G1.add(p1, p1, G1.one());

    mpz_t e;
    mpz_init_set_str(e, "65", 10);

    uint8_t scalar[32];
    for (int i=0;i<32;i++) scalar[i] = 0;
    mpz_export((void *)scalar, NULL, -1, 8, -1, 0, e);
    mpz_clear(e);

    G1Point p2;
    G1.mulByScalar(p2, G1.one(), scalar, 32);

    ASSERT_TRUE(G1.eq(p1,p2));
}

TEST(altBn128, g1_expToOrder) {
    mpz_t e;
    mpz_init_set_str(e, "21888242871839275222246405745257275088548364400416034343698204186575808495617", 10);

    uint8_t scalar[32];

    for (int i=0;i<32;i++) scalar[i] = 0;
    mpz_export((void *)scalar, NULL, -1, 8, -1, 0, e);
    mpz_clear(e);

    G1Point p1;

    G1.mulByScalar(p1, G1.one(), scalar, 32);

    ASSERT_TRUE(G1.isZero(p1));
}

TEST(altBn128, g2_expToOrder) {
    mpz_t e;
    mpz_init_set_str(e, "21888242871839275222246405745257275088548364400416034343698204186575808495617", 10);

    uint8_t scalar[32];

    for (int i=0;i<32;i++) scalar[i] = 0;
    mpz_export((void *)scalar, NULL, -1, 8, -1, 0, e);
    mpz_clear(e);

    Curve<F2Field<RawFq>>::Point p1;

    G2.mulByScalar(p1, G2.one(), scalar, 32);

    ASSERT_TRUE(G2.isZero(p1));
}

TEST(altBn128, multiExp) {

    int NMExp = 40000;

    typedef uint8_t Scalar[32];

    Scalar *scalars = new Scalar[NMExp];
    G1PointAffine *bases = new G1PointAffine[NMExp];

    uint64_t acc=0;
    for (int i=0; i<NMExp; i++) {
        if (i==0) {
            G1.copy(bases[0], G1.one());
        } else {
            G1.add(bases[i], bases[i-1], G1.one());
        }
        for (int j=0; j<32; j++) scalars[i][j] = 0;
        *(int *)&scalars[i][0] = i+1;
        acc += (i+1)*(i+1);
    }

    G1Point p1;
    G1.multiMulByScalar(p1, bases, (uint8_t *)scalars, 32, NMExp);

    mpz_t e;
    mpz_init_set_ui(e, acc);

    Scalar sAcc;

    for (int i=0;i<32;i++) sAcc[i] = 0;
    mpz_export((void *)sAcc, NULL, -1, 8, -1, 0, e);
    mpz_clear(e);

    G1Point p2;
    G1.mulByScalar(p2, G1.one(), sAcc, 32);

    ASSERT_TRUE(G1.eq(p1, p2));

    delete[] bases;
    delete[] scalars;
}

TEST(altBn128, multiExpMSM) {

    int NMExp = 40000;

    typedef uint8_t Scalar[32];

    Scalar *scalars = new Scalar[NMExp];
    G1PointAffine *bases = new G1PointAffine[NMExp];

    uint64_t acc=0;
    for (int i=0; i<NMExp; i++) {
        if (i==0) {
            G1.copy(bases[0], G1.one());
        } else {
            G1.add(bases[i], bases[i-1], G1.one());
        }
        for (int j=0; j<32; j++) scalars[i][j] = 0;
        *(int *)&scalars[i][0] = i+1;
        acc += (i+1)*(i+1);
    }

    G1Point p1;
    G1.multiMulByScalarMSM(p1, bases, (uint8_t *)scalars, 32, NMExp);

    mpz_t e;
    mpz_init_set_ui(e, acc);

    Scalar sAcc;

    for (int i=0;i<32;i++) sAcc[i] = 0;
    mpz_export((void *)sAcc, NULL, -1, 8, -1, 0, e);
    mpz_clear(e);

    G1Point p2;
    G1.mulByScalar(p2, G1.one(), sAcc, 32);

    ASSERT_TRUE(G1.eq(p1, p2));

    delete[] bases;
    delete[] scalars;
}

TEST(altBn128, multiExp2) {

    int NMExp = 2;

    AltBn128::FrElement *scalars = new AltBn128::FrElement[NMExp];
    G1PointAffine *bases = new G1PointAffine[NMExp];

    F1.fromString(bases[0].x, "1626275109576878988287730541908027724405348106427831594181487487855202143055");
    F1.fromString(bases[0].y, "18706364085805828895917702468512381358405767972162700276238017959231481018884");
    F1.fromString(bases[1].x, "17245156998235704504461341147511350131061011207199931581281143511105381019978");
    F1.fromString(bases[1].y, "3858908536032228066651712470282632925312300188207189106507111128103204506804");

    Fr.fromString(scalars[0], "1");
    Fr.fromString(scalars[1], "20187316456970436521602619671088988952475789765726813868033071292105413408473");
    Fr.fromMontgomery(scalars[0], scalars[0]);
    Fr.fromMontgomery(scalars[1], scalars[1]);

    G1Point r;
    G1PointAffine ra;
    G1PointAffine ref;

    F1.fromString(ref.x, "9163953212624378696742080269971059027061360176019470242548968584908855004282");
    F1.fromString(ref.y, "20922060990592511838374895951081914567856345629513259026540392951012456141360");

    G1.multiMulByScalar(r, bases, (uint8_t *)scalars, 32, 2);
    G1.copy(ra, r);

    // std::cout << G1.toString(r, 10);

    ASSERT_TRUE(G1.eq(ra, ref));

    delete[] bases;
    delete[] scalars;
}

TEST(altBn128, multiExp2MSM) {

    int NMExp = 2;

    AltBn128::FrElement *scalars = new AltBn128::FrElement[NMExp];
    G1PointAffine *bases = new G1PointAffine[NMExp];

    F1.fromString(bases[0].x, "1626275109576878988287730541908027724405348106427831594181487487855202143055");
    F1.fromString(bases[0].y, "18706364085805828895917702468512381358405767972162700276238017959231481018884");
    F1.fromString(bases[1].x, "17245156998235704504461341147511350131061011207199931581281143511105381019978");
    F1.fromString(bases[1].y, "3858908536032228066651712470282632925312300188207189106507111128103204506804");

    Fr.fromString(scalars[0], "1");
    Fr.fromString(scalars[1], "20187316456970436521602619671088988952475789765726813868033071292105413408473");
    Fr.fromMontgomery(scalars[0], scalars[0]);
    Fr.fromMontgomery(scalars[1], scalars[1]);

    G1Point r;
    G1PointAffine ra;
    G1PointAffine ref;

    F1.fromString(ref.x, "9163953212624378696742080269971059027061360176019470242548968584908855004282");
    F1.fromString(ref.y, "20922060990592511838374895951081914567856345629513259026540392951012456141360");

    G1.multiMulByScalarMSM(r, bases, (uint8_t *)scalars, 32, 2);
    G1.copy(ra, r);

    // std::cout << G1.toString(r, 10);

    ASSERT_TRUE(G1.eq(ra, ref));

    delete[] bases;
    delete[] scalars;
}

TEST(altBn128, multiExp8MSM) {

    int NMExp = 8;

    AltBn128::FrElement *scalars = new AltBn128::FrElement[NMExp];
    G1PointAffine *bases = new G1PointAffine[NMExp];

    F1.fromString(bases[0].x, "1");
    F1.fromString(bases[0].y, "2");
    F1.fromString(bases[1].x, "1368015179489954701390400359078579693043519447331113978918064868415326638035");
    F1.fromString(bases[1].y, "9918110051302171585080402603319702774565515993150576347155970296011118125764");
    F1.fromString(bases[2].x, "3353031288059533942658390886683067124040920775575537747144343083137631628272");
    F1.fromString(bases[2].y, "19321533766552368860946552437480515441416830039777911637913418824951667761761");
    F1.fromString(bases[3].x, "10744596414106452074759370245733544594153395043370666422502510773307029471145");
    F1.fromString(bases[3].y, "848677436511517736191562425154572367705380862894644942948681172815252343932");
    F1.fromString(bases[4].x, "10415861484417082502655338383609494480414113902179649885744799961447382638712");
    F1.fromString(bases[4].y, "10196215078179488638353184030336251401353352596818396260819493263908881608606");
    F1.fromString(bases[5].x, "4444740815889402603535294170722302758225367627362056425101568584910268024244");
    F1.fromString(bases[5].y, "10537263096529483164618820017164668921386457028564663708352735080900270541420");
    F1.fromString(bases[6].x, "12852522211178622728088728121177131998585782282560100422041774753646305409836");
    F1.fromString(bases[6].y, "15918672909255108529698304535345707578139606904951176064731093256171019744261");
    F1.fromString(bases[7].x, "5841054468380737358126901601208759440531393618333695939860021859990434602332");
    F1.fromString(bases[7].y, "14496198492936799798866613472708085382638873795843716312315177042738122955744");

    Fr.fromString(scalars[0], "1");
    Fr.fromString(scalars[1], "4");
    Fr.fromString(scalars[2], "4");
    Fr.fromString(scalars[3], "4");
    Fr.fromString(scalars[4], "8");
    Fr.fromString(scalars[5], "10");
    Fr.fromString(scalars[6], "1");
    Fr.fromString(scalars[7], "5");
    Fr.fromMontgomery(scalars[0], scalars[0]);
    Fr.fromMontgomery(scalars[1], scalars[1]);
    Fr.fromMontgomery(scalars[2], scalars[2]);
    Fr.fromMontgomery(scalars[3], scalars[3]);
    Fr.fromMontgomery(scalars[4], scalars[4]);
    Fr.fromMontgomery(scalars[5], scalars[5]);
    Fr.fromMontgomery(scalars[6], scalars[6]);
    Fr.fromMontgomery(scalars[7], scalars[7]);

    G1Point r;
    G1PointAffine ra;
    G1PointAffine ref;

    F1.fromString(ref.x, "17747920359253913546551417160303297937542312574889904290131615776238588901697");
    F1.fromString(ref.y, "8815119438581789680513912776342567599606944899217792926373871775002956510503");

    G1.multiMulByScalarMSM(r, bases, (uint8_t *)scalars, 32, NMExp);
    G1.copy(ra, r);

    // std::cout << G1.toString(r, 10);

    ASSERT_TRUE(G1.eq(ra, ref));

    delete[] bases;
    delete[] scalars;
}

TEST(altBn128, fft) {
    int NMExp = 1<<10;

    AltBn128::FrElement *a = new AltBn128::FrElement[NMExp];

    for (int i=0; i<NMExp; i++) {
        Fr.fromUI(a[i], i+1);
    }

    FFT<typename Engine::Fr> fft(NMExp);

    fft.fft(a, NMExp);
    fft.ifft(a, NMExp);

    AltBn128::FrElement aux;
    for (int i=0; i<NMExp; i++) {
        Fr.fromUI(aux, i+1);
        ASSERT_TRUE(Fr.eq(a[i], aux));
    }

    delete[] a;
}


}  // namespace

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
