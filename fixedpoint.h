#ifndef FIXEDPOINT_H
#define FIXEDPOINT_H

#include <stdint.h>
#include <stdlib.h>

#define _fixedpoint_PRECISION   8
#define _fixedpoint_MUL         (0x1<<_fixedpoint_PRECISION)
#define _fixedpoint_RND         (0x1<<(_fixedpoint_PRECISION-1))

class fixedPoint
{
public:
    int64_t val;

    fixedPoint() {}
    fixedPoint(const fixedPoint *that);

    void    randon();

    double  toDouble();
    int64_t round();

    void operator = (double val);
    void operator = (int64_t val);

    void operator += (const fixedPoint &that);
    void operator -= (const fixedPoint &that);
    void operator *= (const fixedPoint &that);
    void operator /= (const fixedPoint &that);

    fixedPoint operator + (fixedPoint &that);
    fixedPoint operator - (fixedPoint &that);
    fixedPoint operator * (fixedPoint &that);
    fixedPoint operator / (fixedPoint &that);
};

#endif // FIXEDPOINT_H
