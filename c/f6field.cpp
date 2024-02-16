#include "splitparstr.hpp"
#include "assert.h"
#include <sstream>

template <typename BaseField>
F6Field<BaseField>::F6Field()
 : F("-1") {
    initField();

    F.fromString(xiTo2PMinus2Over3,
                 "2581911344467009335267311115468803099551665605076196740867805258568234346338,"
                 "19937756971775647987995932169929341994314640652964949448313374472400716661030");

    F.fromString(xiToPMinus1Over3,
                 "21575463638280843010398324269430826099269044274347216827212613867836435027261,"
                 "10307601595873709700152284273816112264069230130616436755625194854815875713954");

    F.F.fromString(xiTo2PSquaredMinus2Over3,
                 "2203960485148121921418603742825762020974279258880205651966");

    F.F.fromString(xiToPSquaredMinus1Over3,
                 "21888242871839275220042445260109153167277707414472061641714758635765020556616");
}

template <typename BaseField>
void F6Field<BaseField>::initField() {
    F.copy(fZero.x, F.zero());
    F.copy(fZero.y, F.zero());
    F.copy(fZero.z, F.zero());
    F.copy(fOne.x, F.zero());
    F.copy(fOne.y, F.zero());
    F.copy(fOne.z, F.one());
    F.copy(fNegOne.x, F.zero());
    F.copy(fNegOne.y, F.zero());
    F.copy(fNegOne.z, F.negOne());
}

template <typename BaseField>
void F6Field<BaseField>::fromString(Element &r, std::string s) {

    auto els = splitParStr(s);
    assert(els.size() == 3);

    F.fromString(r.x, els[0]);
    F.fromString(r.y, els[1]);
    F.fromString(r.z, els[2]);
}

template <typename BaseField>
std::string F6Field<BaseField>::toString(Element &e, uint32_t radix) {
    std::ostringstream stringStream;
    stringStream << "("
                 << F.toString(e.x, radix) << ","
                 << F.toString(e.y, radix) << ","
                 << F.toString(e.z, radix)
                 << ")";
    return stringStream.str();
}

template <typename BaseField>
void F6Field<BaseField>::add(Element &r, Element &a, Element &b) {
    F.add(r.x, a.x, b.x);
    F.add(r.y, a.y, b.y);
    F.add(r.z, a.z, b.z);
}

template <typename BaseField>
void F6Field<BaseField>::sub(Element &r, Element &a, Element &b) {
    F.sub(r.x, a.x, b.x);
    F.sub(r.y, a.y, b.y);
    F.sub(r.z, a.z, b.z);
}

template <typename BaseField>
void F6Field<BaseField>::neg(Element &r, Element &a) {
    F.neg(r.x, a.x);
    F.neg(r.y, a.y);
    F.neg(r.z, a.z);
}

template <typename BaseField>
void F6Field<BaseField>::copy(Element &r, Element &a) {
    F.copy(r.x, a.x);
    F.copy(r.y, a.y);
    F.copy(r.z, a.z);
}

template <typename BaseField>
void F6Field<BaseField>::mul(Element &r, Element &a, Element &b) {

    typename BaseField::Element v0, v1, v2, t0, t1, tz, ty, tx;

    F.mul(v0, a.z, b.z);
    F.mul(v1, a.y, b.y);
    F.mul(v2, a.x, b.x);

    F.add(t0, a.x, a.y);
    F.add(t1, b.x, b.y);
    F.mul(tz, t0, t1);

    F.sub(tz, tz, v1);
    F.sub(tz, tz, v2);
    F.mulXi(tz, tz);
    F.add(tz, tz, v0);

    F.add(t0, a.y, a.z);
    F.add(t1, b.y, b.z);
    F.mul(ty, t0, t1);
    F.sub(ty, ty, v0);
    F.sub(ty, ty, v1);
    F.mulXi(t0, v2);
    F.add(ty, ty, t0);

    F.add(t0, a.x, a.z);
    F.add(t1, b.x, b.z);
    F.mul(tx, t0, t1);
    F.sub(tx, tx, v0);
    F.add(tx, tx, v1);
    F.sub(tx, tx, v2);

    F.copy(r.x, tx);
    F.copy(r.y, ty);
    F.copy(r.z, tz);
}

template <typename BaseField>
void F6Field<BaseField>::mulScalar(Element &r, Element &a, typename BaseField::Element &b) {
    F.mul(r.x, a.x, b);
    F.mul(r.y, a.y, b);
    F.mul(r.z, a.z, b);
}

template <typename BaseField>
void F6Field<BaseField>::mulTau(Element &r, Element &a) {

    typename BaseField::Element tz, ty;

    F.mulXi(tz, a.x);
    F.copy(ty, a.y);
    F.copy(r.y, a.z);
    F.copy(r.x, ty);
    F.copy(r.z, tz);
}

template <typename BaseField>
void F6Field<BaseField>::mulGFP(Element &r, Element &a,  typename BaseField::BaseElement &b) {
    F.mulScalar(r.x, a.x, b);
    F.mulScalar(r.y, a.y, b);
    F.mulScalar(r.z, a.z, b);
}

template <typename BaseField>
void F6Field<BaseField>::square(Element &r, Element &a) {
    typename BaseField::Element v0, v1, v2, c0, c1, c2, xiV2;

    F.square(v0, a.z);
    F.square(v1, a.y);
    F.square(v2, a.x);

    F.add(c0, a.x, a.y);
    F.square(c0, c0);
    F.sub(c0, c0, v1);
    F.sub(c0, c0, v2);
    F.mulXi(c0, c0);
    F.add(c0, c0, v0);

    F.add(c1, a.y, a.z);
    F.square(c1, c1);
    F.sub(c1, c1, v0);
    F.sub(c1, c1, v1);
    F.mulXi(xiV2, v2);
    F.add(c1, c1, xiV2);

    F.add(c2, a.x, a.z);
    F.square(c2, c2);
    F.sub(c2, c2, v0);
    F.add(c2, c2, v1);
    F.sub(c2, c2, v2);

    F.copy(r.x, c2);
    F.copy(r.y, c1);
    F.copy(r.z, c0);
}

template <typename BaseField>
void F6Field<BaseField>::dbl(Element &r, Element &a) {
    F.dbl(r.x, a.x);
    F.dbl(r.y, a.y);
    F.dbl(r.z, a.z);
}

template <typename BaseField>
void F6Field<BaseField>::Frobenius(Element &r, Element &a) {
    F.conjugate(r.x, a.x);
    F.conjugate(r.y, a.y);
    F.conjugate(r.z, a.z);

    F.mul(r.x, r.x, xiTo2PMinus2Over3);
    F.mul(r.y, r.y, xiToPMinus1Over3);
}

template <typename BaseField>
void F6Field<BaseField>::FrobeniusP2(Element &r, Element &a) {
    F.mulScalar(r.x, a.x, xiTo2PSquaredMinus2Over3);
    F.mulScalar(r.y, a.y, xiToPSquaredMinus1Over3);
    F.copy(r.z, a.z);
}

template <typename BaseField>
void F6Field<BaseField>::inv(Element &r, Element &a) {
    typename BaseField::Element t1, A, B, C_, F_;

    F.square(A, a.z);
    F.mul(t1, a.x, a.y);
    F.mulXi(t1, t1);
    F.sub(A, A, t1);

    F.square(B, a.x);
    F.mulXi(B, B);
    F.mul(t1, a.y, a.z);
    F.sub(B, B, t1);

    F.square(C_, a.y);
    F.mul(t1, a.x, a.z);
    F.sub(C_, C_, t1);

    F.mul(F_, C_, a.y);
    F.mulXi(F_, F_);
    F.mul(t1, A, a.z);
    F.add(F_, F_, t1);
    F.mul(t1, B, a.x);
    F.mulXi(t1, t1);
    F.add(F_, F_, t1);

    F.inv(F_, F_);

    F.mul(r.x, C_, F_);
    F.mul(r.y, B, F_);
    F.mul(r.z, A, F_);
}

template <typename BaseField>
bool F6Field<BaseField>::isZero(Element &a) {
    return F.isZero(a.x) && F.isZero(a.y) && F.isZero(a.z);
}

template <typename BaseField>
bool F6Field<BaseField>::isOne(Element &a) {
    return F.isZero(a.x) && F.isZero(a.y) && F.isOne(a.z);
}

template <typename BaseField>
bool F6Field<BaseField>::eq(Element &a, Element &b) {
    return F.eq(a.x, b.x) && F.eq(a.y, b.y) && F.eq(a.z, b.z);
}
