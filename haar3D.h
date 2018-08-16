#ifndef HAAR3D_H
#define HAAR3D_H

#include <stdint.h>
#include "mex.h"
#include <algorithm>
#include <math.h>

/*
struct vali
{
    double	cx;
    double  cy;
    double  cz;
};
*/

void copyAsort(double *VX, double *CX, size_t N, double *C, uint64_t *W, uint64_t *val, uint64_t *ord);
void haar3D(double *inV, double *inC, size_t K, size_t N, size_t depth, double *outCT, double *outW=NULL, uint64_t *outVAL=NULL);
void inv_haar3D(double *inV, double *inCT, size_t K, size_t N, size_t depth, double *outC);

#endif // EXTRA_H
