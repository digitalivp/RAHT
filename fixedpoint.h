#ifndef FIXEDPOINT_H
#define FIXEDPOINT_H

#include <stdint.h>

#define NUMBER_OF_PRECISION_BITS        9
#define USE_ROUNDING_STEAD_OF_FLOORING  1

class fixedPoint
{
public:
    int64_t val;

    fixedPoint();
    fixedPoint(int64_t val);

    void operator = (double val);
    void operator = (const fixedPoint *that);
    void operator = (const fixedPoint &that);
    void operator = (int64_t val);

    int64_t round();
    double  toDouble();

    void operator += (const fixedPoint &that);
    void operator -= (const fixedPoint &that);
    void operator *= (const fixedPoint &that);
    void operator /= (const fixedPoint &that);

    void operator >>= (const int k);
    void operator <<= (const int k);

    fixedPoint operator + (const fixedPoint &that);
    fixedPoint operator - (const fixedPoint &that);
    fixedPoint operator * (const fixedPoint &that);
    fixedPoint operator / (const fixedPoint &that);

    fixedPoint operator >> (const int k);
    fixedPoint operator << (const int k);
};

fixedPoint sqrt(const fixedPoint &S);

#endif // FIXEDPOINT_H
