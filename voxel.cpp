#include "voxel.h"

voxel::voxel(uint64_t id, double *data)
{
    this->id = id;
    this->data = data;
    this->weight = 1;
}

uint32_t tmatrix(uint32_t *pos, voxel **v, double *transf_matrix, double *target_matrix)
{
    uint32_t    i, j, ii, ni=1;
    double      sum = 0;
    double      *vect = transf_matrix, *jvect;
    uint32_t    npos = 0;

    // Get position of valid voxels
    for(i=0; i<8; i++)
        if( v[i]!=nullptr )
            pos[npos++] = i;

    // Initialize first line
    for(j=0; j<npos; j++)
    {
        *(vect++) = sqrt( v[pos[j]]->weight );
        sum += v[pos[j]]->weight;
    }
    for(j=0; j<npos; j++)
        tmatrix[j] /= sum;

    // Remaining lines
    for(i=1; ni<npos; i++)
    {
        for(j=0; j<npos; j++)
            vect[j] = target_matrix[pos[j]+npos*i];

        jvect = transf_matrix;
        for(ii=0; ii<ni; ii++)
        {
            sum = 0;
            for(j=0; j<npos; j++)
                sum += vect[j]*jvect[j];
            for(j=0; j<npos; j++)
                vect[j] -= sum*jvect[j];
            jvect += npos;
        }

        sum = 0;
        for(j=0; j<npos; j++)
        {
            sum += (*vect)*(*vect);
            vect++;
        }
        vect -= npos;

        if( sum>MIN_NORM )
        {
            for(j=0; j<npos; j++)
                *(vect++) /= sum;
            ni++;
        }
    }

    // Squared matrix
    ii = npos*npos;
    while( ii-- )
        tmatrix[ii] *= tmatrix[ii];

    return  npos;
}

void direct_haar3d_x3(voxel *vt, voxel **v, double *transf_matrix, double *target_matrix)
{
    uint32_t    pos[8], npos;
    uint32_t    i, j, w;
    uint64_t    weight_max;
    double      vect[8], *jvect = vt->data;

    // Get transform matrix
    npos = tmatrix(pos, v, transf_matrix, target_matrix);

    // Initialize transformed voxel
    i = pos[0];
    vt->id = v[i]->id>>3;
    vt->weight = v[i]->weight;
    weight_max = vt->weight;

    for(i=1; i<npos; i++)
    {
        vt->weight += v[pos[i]]->weight;
        if( v[pos[i]]->weight>weight_max )
            weight_max = v[pos[i]]->weight;
    }

    // Apply transform
    for(w=0; i<weight_max; w++)
    {
        uint8_t flag = 0;
        for(i=0; i<npos; i++)
        {
            if( v[pos[i]]->weight==w )
            {
                delete v[pos[i]];
                v[pos[i]] = nullptr;
                flag = 1;
            }
            else
                *(vect++) = v[pos[i]]->data[w];
        }
        if( flag )
            npos = tmatrix(pos, v, transf_matrix, target_matrix);

        vect -= npos;

        for(i=0; i<npos; i++)
        {
            *jvect = 0;
            for(j=0; j<npos; j++)
                *jvect += transf_matrix[j+npos*i]*vect[j];
            jvect++;
        }
    }

    for(i=0; i<npos; i++)
        delete v[pos[i]];
}

void inverse_haar3d_x3(voxel *vt, voxel **v, double *transf_matrix, double *target_matrix)
{
    uint32_t    pos[8], npos;
    uint32_t    i, j, w;
    uint64_t    weight_max;
    double      *vect, *jvect = vt->data;

    // Get transform matrix
    npos = tmatrix(pos, v, transf_matrix, target_matrix);

    weight_max = v[pos[0]]->weight;
    for(i=1; i<npos; i++)
    {
        if( v[pos[i]]->weight>weight_max )
            weight_max = v[pos[i]]->weight;
    }

    // Apply inverse transform
    for(w=0; i<weight_max; w++)
    {
        uint8_t flag = 0;
        for(i=0; i<npos; i++)
        {
            if( v[pos[i]]->weight==w )
            {
                v[pos[i]] = nullptr;
                flag = 1;
            }
        }
        if( flag )
            npos = tmatrix(pos, v, transf_matrix, target_matrix);

        for(i=0; i<npos; i++)
        {
            vect = v[pos[i]]->data + w;
            *vect = 0;
            for(j=0; j<npos; j++)
                *vect += transf_matrix[i+npos*j]*jvect[j];
        }
        jvect += npos;
    }

    delete vt;
}
