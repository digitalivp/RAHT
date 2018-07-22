#ifndef VOXEL_H
#define VOXEL_H

#include <stdint.h>
#include <math.h>

#define MIN_NORM    1e-5

class voxel
{
public:
    uint64_t    id;
    uint64_t    weight;
    double      *data;

    voxel(uint64_t id, double *data);
};

uint32_t tmatrix(uint32_t *pos, voxel **v, double *transf_matrix, double *target_matrix);
void direct_haar3d_x3(voxel *vt, voxel **v, double *transf_matrix, double *target_matrix);

#endif // VOXEL_H
