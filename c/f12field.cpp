#include "splitparstr.hpp"
#include "assert.h"
#include <sstream>

template <typename BaseField>
F12Field<BaseField>::F12Field() {
    initField();
    F.F.fromString(xiToPMinus1Over6,
                   "8376118865763821496583973867626364092589906065868298776909617916018768340080,"
                   "16469823323077808223889137241176536799009286646108169935659301613961712198316");

    F.F.F.fromString(xiToPSquaredMinus1Over6,
                  "21888242871839275220042445260109153167277707414472061641714758635765020556617");
}

template <typename BaseField>
void F12Field<BaseField>::initField() {
    F.copy(fZero.x, F.zero());
    F.copy(fZero.y, F.zero());
    F.copy(fOne.x, F.zero());
    F.copy(fOne.y, F.one());
    F.copy(fNegOne.x, F.zero());
    F.copy(fNegOne.y, F.negOne());
}

template <typename BaseField>
void F12Field<BaseField>::fromString(Element &r, std::string s) {

    auto els = splitParStr(s);
    assert(els.size() == 2);

    F.fromString(r.x, els[0]);
    F.fromString(r.y, els[1]);
}

template <typename BaseField>
std::string F12Field<BaseField>::toString(Element &e, uint32_t radix) {
    std::ostringstream stringStream;
    stringStream << "(" << F.toString(e.x, radix) << "," << F.toString(e.y, radix) << ")";
    return stringStream.str();
}

template <typename BaseField>
void F12Field<BaseField>::add(Element &r, Element &a, Element &b) {
    F.add(r.x, a.x, b.x);
    F.add(r.y, a.y, b.y);
}

template <typename BaseField>
void F12Field<BaseField>::sub(Element &r, Element &a, Element &b) {
    F.sub(r.x, a.x, b.x);
    F.sub(r.y, a.y, b.y);
}

template <typename BaseField>
void F12Field<BaseField>::neg(Element &r, Element &a) {
    F.neg(r.x, a.x);
    F.neg(r.y, a.y);
}

template <typename BaseField>
void F12Field<BaseField>::conjugate(Element &r, Element &a) {
    F.neg(r.x, a.x);
    F.copy(r.y, a.y);
}

template <typename BaseField>
void F12Field<BaseField>::copy(Element &r, Element &a) {
    F.copy(r.x, a.x);
    F.copy(r.y, a.y);
}

template <typename BaseField>
void F12Field<BaseField>::mul(Element &r, Element &a, Element &b) {
    typename BaseField::Element tx, t, ty;

    F.mul(tx, a.x, b.y);
    F.mul(t, b.x, a.y);
    F.add(tx, tx, t);

    F.mul(ty, a.y, b.y);
    F.mul(t, a.x, b.x);
    F.mulTau(t, t);
    F.add(r.y, ty, t);
    F.copy(r.x, tx);
}

template <typename BaseField>
void F12Field<BaseField>::mulScalar(Element &r, Element &a, typename BaseField::Element &b) {
    F.mul(r.x, a.x, b);
    F.mul(r.y, a.y, b);
}

#define BIT_IS_SET(s, p) (s[p>>3] & (1 << (p & 0x7)))
template <typename BaseField>
void F12Field<BaseField>::exp(Element &r, Element &a, uint8_t* scalar, unsigned int scalarSize) {

    copy(r, one());
    Element t;

    for (int i=scalarSize*8-1; i>=0; i--) {
        square(t, r);

        if ( BIT_IS_SET(scalar, i) ) {
            mul(r, t, a);
        } else {
            copy(r, t);
        }
    }
}

template <typename BaseField>
void F12Field<BaseField>::square(Element &r, Element &a) {
    typename BaseField::Element v0, t, ty;

    F.mul(v0, a.x, a.y);

    F.mulTau(t, a.x);
    F.add(t, a.y, t);
    F.add(ty, a.x, a.y);
    F.mul(ty, ty, t);
    F.sub(ty, ty, v0);
    F.mulTau(t, v0);
    F.sub(ty, ty, t);

    F.copy(r.y, ty);
    F.dbl(r.x, v0);
}

template <typename BaseField>
void F12Field<BaseField>::Frobenius(Element &r, Element &a) {
   F.Frobenius(r.x, a.x);
   F.Frobenius(r.y, a.y);
   F.mulScalar(r.x, r.x, xiToPMinus1Over6);
}

template <typename BaseField>
void F12Field<BaseField>::FrobeniusP2(Element &r, Element &a) {
   F.FrobeniusP2(r.x, a.x);
   F.mulGFP(r.x, r.x, xiToPSquaredMinus1Over6);
   F.FrobeniusP2(r.y, a.y);
}

template <typename BaseField>
void F12Field<BaseField>::inv(Element &r, Element &a) {
    typename BaseField::Element t1, t2;

    F.square(t1, a.x);
    F.square(t2, a.y);
    F.mulTau(t1, t1);
    F.sub(t1, t2, t1);
    F.inv(t2, t1);

    F.neg(r.x, a.x);
    F.copy(r.y, a.y);
    mulScalar(r, r, t2);
}

template <typename BaseField>
bool F12Field<BaseField>::isZero(Element &a) {
    return F.isZero(a.x) && F.isZero(a.y);
}

template <typename BaseField>
bool F12Field<BaseField>::isOne(Element &a) {
    return F.isZero(a.x) && F.isOne(a.y);
}

template <typename BaseField>
bool F12Field<BaseField>::eq(Element &a, Element &b) {
    return F.eq(a.x, b.x) && F.eq(a.y, b.y);
}
