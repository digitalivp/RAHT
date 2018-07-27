#include "haar3D.h"

int64_t __sqrtIF(int64_t P)
{
#if EMPLOY_STEP_CLEANING
    if( P>=(0x1<<(2*INTEGER_STEPSIZE_PRECISION)) )
        return 0x1<<INTEGER_STEPSIZE_PRECISION;
#endif

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
COMPUTE_ERROR:
    error = (P/A-A)/2;
    if( error )
        goto UPDATE_VALUE;
    return A;

UPDATE_VALUE:
    A += error;
    if( !A )
        A = 1;
    goto COMPUTE_ERROR;
}

int64_t sqrtIF(int64_t A, int64_t W0, int64_t W1)
{
    A *= A;
#if INVERSE_SQUARE_ROOT
    return __sqrtIF((0x1<<(2*INTEGER_STEPSIZE_PRECISION))*(W0*W1)/(W0+W1) * (0x1<<(2*INTEGER_STEPSIZE_PRECISION))/A);
#else
    return __sqrtIF(A/W0 + A/W1);
#endif
}

int64_t sqrtIF(int64_t A, int64_t W)
{
    A *= A;
#if INVERSE_SQUARE_ROOT
    return __sqrtIF((0x1<<(2*INTEGER_STEPSIZE_PRECISION))*W/A*(0x1<<(2*INTEGER_STEPSIZE_PRECISION)));
#else
    return __sqrtIF(A/W);
#endif
}
