#include "haar3D.h"

int64_t __sqrtIF(int64_t P)
{
#if INVERSE_SQUARE_ROOT && EMPLOY_STEP_CLEANING
    if( P>>(2*NUMBER_OF_PRECISION_BITS) )
        return 0x1<<NUMBER_OF_PRECISION_BITS;
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

int64_t sqrtIF(int64_t A, int64_t W0, int64_t W1)
{
#if INVERSE_SQUARE_ROOT
    return __sqrtIF((0x1<<(2*NUMBER_OF_PRECISION_BITS))*(W0*W1)/(W0+W1) * (0x1<<(2*NUMBER_OF_PRECISION_BITS))/(A*A));
#else
    return __sqrtIF((A*A*(W0+W1))/(W0*W1));
#endif
}

int64_t sqrtIF(int64_t A, int64_t W)
{
    A *= A;
#if INVERSE_SQUARE_ROOT
    return __sqrtIF((0x1<<(2*NUMBER_OF_PRECISION_BITS))*W/A*(0x1<<(2*NUMBER_OF_PRECISION_BITS)));
#else
    return __sqrtIF(A/W);
#endif
}
