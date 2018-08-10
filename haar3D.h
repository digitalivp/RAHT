#ifndef HAAR3D_H
#define HAAR3D_H

#include <stdint.h>
#include <algorithm>
#include <math.h>
#include "mex.h"
#include "fixedpoint.h"

#define INVERSE_SQUARE_ROOT         1
#define EMPLOY_STEP_CLEANING        0

int64_t sqrtIF(int64_t A, int64_t W0, int64_t W1);
int64_t sqrtIF(int64_t A, int64_t W);

void copyAsort(double *VX, double *CX, size_t N, double *C, uint64_t *W, uint64_t *val, uint64_t *ord);
void haar3D(fixedPoint Qstep, double *inV, double *inC, size_t K, size_t N, size_t depth, int64_t *outCT);
void inv_haar3D(fixedPoint Qstep, double *inV, int64_t *inCT, size_t K, size_t N, size_t depth, double *outC);

#endif // EXTRA_H
