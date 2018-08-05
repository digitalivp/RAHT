#include "fixedpoint.h"

fixedPoint::fixedPoint()
{
    this->val = 0;
}

fixedPoint::fixedPoint(int64_t val)
{
    this->val = val;
}

void fixedPoint::operator = (double val)
{
    this->val = (val*(0x1<<NUMBER_OF_PRECISION_BITS)) + 0.5;
}

void fixedPoint::operator = (const fixedPoint *that)
{
    this->val = that->val;
}

void fixedPoint::operator = (const fixedPoint &that)
{
    this->val = that.val;
}

void fixedPoint::operator = (int64_t val)
{
    this->val = val;
    *this <<= NUMBER_OF_PRECISION_BITS;
}

int64_t fixedPoint::round()
{
    if( this->val>0 )
        return (this->val + (0x1<<(NUMBER_OF_PRECISION_BITS-1))) >> NUMBER_OF_PRECISION_BITS;
    return -((-this->val + (0x1<<(NUMBER_OF_PRECISION_BITS-1))) >> NUMBER_OF_PRECISION_BITS);
}

double fixedPoint::toDouble()
{
    return (double) this->val / (0x1<<NUMBER_OF_PRECISION_BITS);
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

#if USE_ROUNDING_STEAD_OF_FLOORING
    if( this->val>0 )
        this->val = +(+this->val + (0x1<<(NUMBER_OF_PRECISION_BITS-1)))>>NUMBER_OF_PRECISION_BITS;
    else
        this->val = -(-this->val + (0x1<<(NUMBER_OF_PRECISION_BITS-1)))>>NUMBER_OF_PRECISION_BITS;
#else
    if( this->val>0 )
        this->val >>= NUMBER_OF_PRECISION_BITS;
    else
        this->val = -(-this->val)>>NUMBER_OF_PRECISION_BITS;
#endif
}

void fixedPoint::operator /= (const fixedPoint &that)
{
    uint8_t sign_flag;

    if( this->val>0 )
    {
        this->val <<= NUMBER_OF_PRECISION_BITS;
        sign_flag = 0;
    }
    else
    {
        this->val = (-this->val)<<NUMBER_OF_PRECISION_BITS;
        sign_flag = 1;
    }

    if( that.val>0 )
    {
#if USE_ROUNDING_STEAD_OF_FLOORING
        this->val += (1+that.val)>>1;
#endif
        this->val /= that.val;
    }
    else
    {
#if USE_ROUNDING_STEAD_OF_FLOORING
        this->val += (1-that.val)>>1;
#endif
        this->val /= -that.val;
        sign_flag = sign_flag ? 0 : 1;
    }

    if( sign_flag )
        this->val = -this->val;
}

void fixedPoint::operator >>= (const int k)
{
    if( this->val>0 )
        this->val >>= k;
    else
        this->val = -((-this->val)>>k);
}

void fixedPoint::operator <<= (const int k)
{
    if( this->val>0 )
        this->val <<= k;
    else
        this->val = -((-this->val)<<k);
}

fixedPoint fixedPoint::operator + (const fixedPoint &that)
{
    fixedPoint ans(this->val);
    ans += that;
    return ans;
}

fixedPoint fixedPoint::operator - (const fixedPoint &that)
{
    fixedPoint ans(this->val);
    ans -= that;
    return ans;
}

fixedPoint fixedPoint::operator * (const fixedPoint &that)
{
    fixedPoint ans(this->val);
    ans *= that;
    return ans;
}

fixedPoint fixedPoint::operator / (const fixedPoint &that)
{
    fixedPoint ans(this->val);
    ans /= that;
    return ans;
}

fixedPoint fixedPoint::operator >> (const int k)
{
    fixedPoint ans(this->val);
    ans >>= k;
    return ans;
}

fixedPoint fixedPoint::operator << (const int k)
{
    fixedPoint ans(this->val);
    ans <<= k;
    return ans;
}

fixedPoint sqrt(const fixedPoint &S)
{
    fixedPoint error;
    fixedPoint A;

    // Initial guess
    A.val = 1;
    int64_t p = S.val<<NUMBER_OF_PRECISION_BITS;

    while( p )
    {
        p >>= 2;
        A.val <<= 1;
    }

    // Newton-Raphson
    while( 1 )
    {
        error = S;
        error /= A;
        error -= A;
        error >>= 1;

        if( !error.val )
            return A;

        A += error;
        if( !A.val )
            A.val = 1;
    }
}
