#ifndef HAAR3D_H
#define HAAR3D_H

#include <stdint.h>
#include <algorithm>
#include <math.h>
#include "fixedpoint.h"

struct _weight
{
    int32_t *val;
    size_t count;
};
    
int64_t sqrtIF(int64_t A, uint64_t W0, uint64_t W1);
int64_t sqrtIF(int64_t A, uint64_t W);

void copyAsort(double *V, size_t N, uint64_t *W, uint64_t *inW, uint64_t *val, uint64_t *ord);
void haar3D(fixedPoint Qstep, double *inV, uint8_t *inC, uint64_t *inW, size_t N, intmax_t *outCT);
void inv_haar3D(fixedPoint Qstep, double *inV, intmax_t *inCT, uint64_t *inW, size_t N, uint8_t *outC);

// f2b -> index_derivate
intmax_t *index_derivate(size_t N, uint8_t *index);

// b2f -> index_integrate
uint8_t *index_integrate(size_t N, intmax_t *derivate);

// f2w -> index2weight
void index2weight(size_t N, uint8_t *index, _weight weight, uint64_t *w);

#endif // EXTRA_H