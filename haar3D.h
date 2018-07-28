#ifndef HAAR3D_H
#define HAAR3D_H

#include <stdint.h>
#include "mex.h"
#include <algorithm>
#include <math.h>

#define INTEGER_TRANSFORM_PRECISION 10
#define INTEGER_STEPSIZE_PRECISION  10
#define INVERSE_SQUARE_ROOT         0
#define EMPLOY_STEP_CLEANING        0

int64_t sqrtIF(int64_t A, int64_t W0, int64_t W1);
int64_t sqrtIF(int64_t A, int64_t W);

void copyAsort(double *VX, double *CX, size_t N, double *C, uint64_t *W, uint64_t *val, uint64_t *ord);
void haar3D(int64_t Qstep, double *inV, double *inC, size_t K, size_t N, size_t depth, double *outCT, double *outW=NULL);
void inv_haar3D(int64_t Qstep, double *inV, double *inCT, size_t K, size_t N, size_t depth, double *outC);

#endif // EXTRA_H
