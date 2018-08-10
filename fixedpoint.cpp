#include "fixedpoint.h"

fixedPoint::fixedPoint(const fixedPoint *that)
{
    this->val = that->val;
}

void fixedPoint::randon()
{
    if( rand()%2 )
        this->val = -rand()%65536;
    else
        this->val = rand()%65536;
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

void fixedPoint::operator = (double val)
{
    if( val>0 )
        this->val = val*_fixedpoint_MUL + 0.5;
    else
        this->val = val*_fixedpoint_MUL - 0.5;
}

void fixedPoint::operator = (int64_t val)
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

    if( this->val > 0 )
    {
        this->val += _fixedpoint_RND;
        this->val >>= _fixedpoint_PRECISION;
    }
    else
    {
        this->val = _fixedpoint_RND - this->val;
        this->val >>= _fixedpoint_PRECISION;
        this->val = -this->val;
    }
}

void fixedPoint::operator /= (const fixedPoint &that)
{
    uint8_t signflag = 0;

    if( this->val < 0 )
    {
        this->val = -this->val;
        signflag = 1;
    }

    this->val <<= _fixedpoint_PRECISION;

    if( that.val < 0 )
    {
        this->val += (-that.val)>>1;
        this->val /= -that.val;
        signflag = signflag ? 0 : 1;
    }
    else
    {
        this->val += that.val>>1;
        this->val /= that.val;
    }

    if( signflag )
        this->val = -this->val;
}

fixedPoint fixedPoint::operator + (fixedPoint &that)
{
    fixedPoint  res(this);
    res += that;
    return res;
}

fixedPoint fixedPoint::operator - (fixedPoint &that)
{
    fixedPoint  res(this);
    res -= that;
    return res;
}

fixedPoint fixedPoint::operator * (fixedPoint &that)
{
    fixedPoint  res(this);
    res *= that;
    return res;
}

fixedPoint fixedPoint::operator / (fixedPoint &that)
{
    fixedPoint  res(this);
    res /= that;
    return res;
}
