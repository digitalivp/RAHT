#ifndef FIXEDPOINT_H
#define FIXEDPOINT_H

#include <stdint.h>
#include <stdlib.h>

/* 32 > _fixedpoint_PRECISION > 0.
 * A integer representing the number of bits after the point */
#define _fixedpoint_PRECISION   8

#define _fixedpoint_MUL         (((int64_t) 0x1)<<_fixedpoint_PRECISION)
#define _fixedpoint_RND         (((int64_t) 0x1)<<(_fixedpoint_PRECISION-1))
#define _fixedpoint_MUL2        (((int64_t) 0x1)<<(_fixedpoint_PRECISION*2))

class fixedPoint
{
public:
    int64_t val;

    fixedPoint() {}
    fixedPoint(const double val);
    fixedPoint(const fixedPoint *that);

    void    randon();

    double  toDouble();
    int64_t round();

    void operator = (const double val);
    void operator = (const int64_t val);

    void operator += (const fixedPoint &that);
    void operator -= (const fixedPoint &that);
    void operator *= (const fixedPoint &that);
    void operator /= (const fixedPoint &that);

    fixedPoint operator + (const fixedPoint &that);
    fixedPoint operator - (const fixedPoint &that);
    fixedPoint operator * (const fixedPoint &that);
    fixedPoint operator / (const fixedPoint &that);
};

int64_t _sqrt(int64_t P);

#endif // FIXEDPOINT_H
