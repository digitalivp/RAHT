#ifndef HAAR3D_H
#define HAAR3D_H

#include <stdint.h>
#include "mex.h"
#include <algorithm>
#include <math.h>

struct vali
{
    double	cx;
    double  cy;
    double  cz;
};

void copyAsort(double *VX, double *CX, uint64_t N, vali *C, uint64_t *W, uint64_t *val, uint64_t *ord);
void haar3D(double *inV, double *inC, uint64_t N, uint64_t depth, double *outCT, double *outW=NULL);
void inv_haar3D(double *inV, double *inCT, uint64_t N, uint64_t depth, double *outC);

#endif // EXTRA_H
