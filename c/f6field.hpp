#include <string>

template <typename BaseField>
class F6Field {

public:
    struct Element {
        typename BaseField::Element x;
        typename BaseField::Element y;
        typename BaseField::Element z;
    };

    BaseField F;

    typedef BaseField Base;
    typedef typename BaseField::Element BaseElement;

private:
    Element fOne;
    Element fZero;
    Element fNegOne;
    typename BaseField::Element xiTo2PMinus2Over3;
    typename BaseField::Element xiToPMinus1Over3;
    typename BaseField::BaseElement xiTo2PSquaredMinus2Over3;
    typename BaseField::BaseElement xiToPSquaredMinus1Over3;

    void initField();
public:

    F6Field();

    Element &zero() { return fZero; }
    Element &one() { return fOne; }
    Element &negOne() { return fNegOne; }

    void copy(Element &r, Element &a);
    void add(Element &r, Element &a, Element &b);
    void sub(Element &r, Element &a, Element &b);
    void neg(Element &r, Element &a);
    void mul(Element &r, Element &a, Element &b);
    void mulScalar(Element &r, Element &a, typename BaseField::Element &b);
    void mulTau(Element &r, Element &a);
    void mulGFP(Element &r,Element &a, typename BaseField::BaseElement& b);
    void square(Element &r, Element &a);
    void dbl(Element &r, Element &a);
    void Frobenius(Element &r, Element &a);
    void FrobeniusP2(Element &r, Element &a);
    void inv(Element &r, Element &a);
    bool isZero(Element &a);
    bool isOne(Element &a);
    bool eq(Element &a, Element &b);

    void fromString(Element &r, std::string s);
    std::string toString(Element &a, uint32_t radix = 10);

};

#include "f6field.cpp"
