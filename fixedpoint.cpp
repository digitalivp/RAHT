#include "fixedpoint.h"

fixedPoint::fixedPoint(const double val)
{
    *this = val;
}

fixedPoint::fixedPoint(const fixedPoint *that)
{
    this->val = that->val;
}

void fixedPoint::randon()
{
    this->val = rand()%_fixedpoint_MUL;
    if( rand()%2 )
        this->val = -this->val;
}

double fixedPoint::toDouble()
{
    return (double) this->val/_fixedpoint_MUL;
}

int64_t fixedPoint::round()
{
    if( this->val > 0 )
        return (_fixedpoint_RND + this->val)>>_fixedpoint_PRECISION;
    return -((_fixedpoint_RND - this->val)>>_fixedpoint_PRECISION);
}

void fixedPoint::operator = (const double val)
{
    if( val>0 )
        this->val = val*_fixedpoint_MUL + 0.5;
    else
        this->val = val*_fixedpoint_MUL - 0.5;
}

void fixedPoint::operator = (const int64_t val)
{
    if( val > 0 )
        this->val = val<<_fixedpoint_PRECISION;
    else
        this->val = -((-val)<<_fixedpoint_PRECISION);
}

void fixedPoint::operator += (const fixedPoint &that)
{
    this->val += that.val;
}

void fixedPoint::operator -= (const fixedPoint &that)
{
    this->val -= that.val;
}

void fixedPoint::operator *= (const fixedPoint &that)
{
    this->val *= that.val;

    if( this->val < 0 )
        this->val = -( (_fixedpoint_RND - this->val) >> _fixedpoint_PRECISION );
    else
        this->val = +( (_fixedpoint_RND + this->val) >> _fixedpoint_PRECISION );
}

void fixedPoint::operator /= (const fixedPoint &that)
{
    if( this->val < 0 )
    {
        if( that.val < 0 )
            this->val = -(
                    ((-that.val)>>1) +
                    ((-this->val)<<_fixedpoint_PRECISION) )/that.val;
        else
            this->val = -(
                    ((+that.val)>>1) +
                    ((-this->val)<<_fixedpoint_PRECISION) )/that.val;
    }
    else
    {
        if( that.val<0 )
            this->val = +(
                    ((-that.val)>>1) +
                    ((+this->val)<<_fixedpoint_PRECISION) )/that.val;
        else
            this->val = +(
                    ((+that.val)>>1) +
                    ((+this->val)<<_fixedpoint_PRECISION) )/that.val;
    }
}

fixedPoint fixedPoint::operator + (const fixedPoint &that)
{
    fixedPoint  res(this);
    res += that;
    return res;
}

fixedPoint fixedPoint::operator - (const fixedPoint &that)
{
    fixedPoint  res(this);
    res -= that;
    return res;
}

fixedPoint fixedPoint::operator * (const fixedPoint &that)
{
    fixedPoint  res(this);
    res *= that;
    return res;
}

fixedPoint fixedPoint::operator / (const fixedPoint &that)
{
    fixedPoint  res(this);
    res /= that;
    return res;
}

int64_t _sqrt(int64_t P)
{
    int64_t error;

    // Initial guess
    int64_t p = P;
    int64_t A = 1;
    while( p )
    {
        p >>= 2;
        A <<= 1;
    }

    // Newton-Raphson
    while( 1 )
    {
        error = (P/A-A)/2;
        if( !error )
            return A;
        A += error;
        if( !A )
            A = 1;
    }
}
