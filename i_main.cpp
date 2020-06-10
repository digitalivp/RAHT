#include <stdint.h>
#include "haar3D.h"
#include <algorithm>
#include <math.h>

#include <sstream>
#include <fstream>
#include <iostream>

#include "file.h"

#include <thread>

using namespace std;

#include "mex.h"

/* 0             0  1               2
 * Cd = RAHT_dec(V, binaryfilepath, weight_roi_vector)
 * Cd = RAHT_dec(V, binaryfilepath)
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
// VERIFICA ENTRADAS E SAÍDAS --------------------------------------------------------------------
    if( !nrhs && !nlhs ) {
        printf("Cd = RAHT_dec(V, binaryfilepath)\n\n"    );
        printf("INPUTS:\n"                                                                      );
        printf("  V ............. (Nx3  double) Voxels vertices\n"                              );
        printf("  binaryfilepath  (char string) file path the the encoded binary file\n\n"      );
        printf("OUTPUTS:\n"                                                                     );
        printf("  Cd ............ (Nx1   uint8) decoded voxels colors (RGB)\n"                  );
        return;
    }
    if( nlhs>1 )
        mexErrMsgTxt("Expected 1 output");
    if( nrhs!=2 )
        mexErrMsgTxt("Expected 2 inputs");

    if( mxGetClassID(prhs[0])!=mxDOUBLE_CLASS )
        mexErrMsgTxt("First input should be DOUBLE");
    if( mxGetClassID(prhs[1])!=mxCHAR_CLASS )
        mexErrMsgTxt("Second input should be CHAR");

    if( mxGetN(prhs[0])!=3 )
        mexErrMsgTxt("First input should have 3 columns");

// ABRE ARQUIVO BINÁRIO ---------------------------------------------------------------------------
    char *binaryfilepath = mxArrayToString(prhs[1]);
    file *fid = new file(binaryfilepath, 0);
    mxFree(binaryfilepath);

    if( fid->openError() ) {
        delete fid;
        mexErrMsgTxt("Unable to open binaryfilepath");
    }

// OBTEM DADOS ------------------------------------------------------------------------------------
    size_t     N;
    fixedPoint Qstep;
    intmax_t   *CT;
    double     *V = mxGetPr(prhs[0]);
    uint8_t    *C;
    _weight    weight;

    N = fid->grRead(20);
    weight.count = fid->grRead(7);
    
    printf("N = %d\n", N);
    printf("weight.count = %d\n", weight.count);
    if( weight.count )
    {
        weight.val = new int32_t[weight.count];
        for(size_t n=0; n<weight.count; n++)
        {
            weight.val[n] = fid->grRead(7)+1;
            printf("  %d\n", weight.val[n]);
        }
            
    }
    uint64_t *W = new uint64_t[N];

    fid->read(&Qstep.val, sizeof(int64_t), 1);
    printf("Qstep = %lf\n", Qstep.toDouble());

    if( weight.count ) {
        intmax_t *derivate = new intmax_t[N];
        fid->rlgrRead(derivate, N);
        uint8_t *index = index_integrate(N, derivate);
        index2weight(N, index, weight, W);

        delete [] derivate;
        delete [] index;
    } else {
        for(size_t n=0; n<N; n++)
            W[n] = 1;
    }

    if( mxGetM(prhs[0])!=N ) {
        delete fid;
        delete [] W;
        mexErrMsgTxt("First input (aka V) and the given file seem to be incompatible");
    }

    plhs[0] = mxCreateNumericMatrix(N, 3, mxUINT8_CLASS, mxREAL);
    C = reinterpret_cast<uint8_t *>(mxGetPr(plhs[0]));

    CT = new intmax_t[N*3];

    fid->rlgrRead(CT, N*3);

    inv_haar3D(Qstep, V, CT, W, N, C);
    delete [] CT;
    delete fid;
    
    if( weight.count )
        delete [] weight.val;
}
