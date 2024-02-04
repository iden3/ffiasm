#include <string>

template <typename BaseField>
class F12Field {

public:
    struct Element {
        typename BaseField::Element x;
        typename BaseField::Element y;
    };

    BaseField F;

private:
    Element fOne;
    Element fZero;
    Element fNegOne;
    typename BaseField::BaseElement xiToPMinus1Over6;
    typename BaseField::Base::BaseElement xiToPSquaredMinus1Over6;

    void initField();
public:

    F12Field();

    Element &zero() { return fZero; }
    Element &one() { return fOne; }
    Element &negOne() { return fNegOne; }

    void copy(Element &r, Element &a);
    void add(Element &r, Element &a, Element &b);
    void sub(Element &r, Element &a, Element &b);
    void neg(Element &r, Element &a);
    void conjugate(Element &r, Element &a);
    void mul(Element &r, Element &a, Element &b);
    void mulScalar(Element &r, Element &a, typename BaseField::Element &b);
    void exp(Element &r, Element &a, uint8_t* scalar, unsigned int scalarSize);
    void square(Element &r, Element &a);
    void Frobenius(Element &r, Element &a);
    void FrobeniusP2(Element &r, Element &a);
    void inv(Element &r, Element &a);
    bool isZero(Element &a);
    bool isOne(Element &a);
    bool eq(Element &a, Element &b);

    void fromString(Element &r, std::string s);
    std::string toString(Element &a, uint32_t radix = 10);

};

#include "f12field.cpp"
