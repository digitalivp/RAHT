#ifndef HAAR3D_H
#define HAAR3D_H

#include "inteiros.h"
#include "mex.h"
#include <algorithm>
#include <math.h>

struct vali
{
    double	cx;
    double  cy;
    double  cz;
};

void copyAsort(double *VX, double *CX, uint64 N, vali *C, uint64 *W, uint64 *val, uint64 *ord);
void haar3D(double *inV, double *inC, uint64 N, uint64 depth, double *outCT, double *outW=NULL);
void inv_haar3D(double *inV, double *inCT, uint64 N, uint64 depth, double *outC);

#endif // EXTRA_H
