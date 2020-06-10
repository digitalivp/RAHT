#include <stdint.h>
#include "haar3D.h"
#include <algorithm>
#include <math.h>

#include <sstream>
#include <fstream>
#include <iostream>

#include "file.h"
#include "fixedpoint.h"

#define	debug	0

using namespace std;

#include "mex.h"

/*
 *                      0  1  2               3       4          5
 * file_size = RAHT_cod(V, C, binaryfilepath, Qstep, [index_roi, weight_roi])
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
// VERIFICA ENTRADAS E SAÍDAS --------------------------------------------------------------------
    if( (!nrhs && !nlhs) || (nrhs!=4 && nrhs!=6) || (nlhs>1) ) {
        printf("file_size = RAHT_cod(V, C, binaryfilepath, Qstep, [index_roi, weight_roi])\n"    );
        printf("INPUTS:\n"                                                                      );
        printf("  V ............. (Nx3  double) Voxels vertices\n"                              );
        printf("  C ............. (Nx3   uint8) Voxels colors (RGB)\n"                          );
        printf("  binaryfilepath  (char string) file path the the encoded binary file\n"        );
        printf("  Qstep ......... (1x1  double) Quantization step\n\n"                          );
        printf("  index_roi ..... (Nx1   uint8) an index indicating if a voxels is in the ROI\n");
        printf("                  index = 0 -> not in the ROI, weight = 1\n"                    );
        printf("                  index = n (n>0) -> in the ROI, weight = weight_roi(n)\n"      );
        printf("  weight_roi .... (1xK  double) weight of voxels in the ROI\n"                  );
        printf("OUTPUTS:\n"                                                                     );
        printf("  file_size ..... (1x1  uint64) number of bits used to encode the binary file\n");

        if( nrhs || nlhs )
            mexErrMsgTxt("Invalid number of inputs or outputs");

        return;
    }

    if( mxGetClassID(prhs[0])!=mxDOUBLE_CLASS )
        mexErrMsgTxt("First input should be DOUBLE");
    if( mxGetClassID(prhs[1])!=mxUINT8_CLASS )
        mexErrMsgTxt("Second input should be UINT8");
    if( mxGetClassID(prhs[2])!=mxCHAR_CLASS )
        mexErrMsgTxt("Third input should be CHAR");
    if( mxGetClassID(prhs[3])!=mxDOUBLE_CLASS )
        mexErrMsgTxt("Fourth input should be DOUBLE");
    if( mxGetN(prhs[0])!=3 )
        mexErrMsgTxt("First input should have 3 columns");
    if( mxGetN(prhs[1])!=3 )
        mexErrMsgTxt("Second input should have 3 columns");

    size_t N = mxGetM( prhs[0] );
    if( N!=mxGetM(prhs[1]) )
        mexErrMsgTxt("First two inputs should have the same number of lines");

    if( nrhs==6 ) {
        if( mxGetClassID(prhs[4])!=mxUINT8_CLASS && mxGetClassID(prhs[4])!=mxLOGICAL_CLASS )
            mexErrMsgTxt("fifth input should be UINT8");
        if( mxGetClassID(prhs[5])!=mxDOUBLE_CLASS )
            mexErrMsgTxt("Sixth input should be DOUBLE");
        if( N!=mxGetM(prhs[4]) )
            mexErrMsgTxt("First and fifth inputs should have the same number of lines");
        if( mxGetN(prhs[4])!=1 )
            mexErrMsgTxt("Third input should have 1 column");
    }

// OBTEM DADOS ----------------------------------------------------------------------------------
    //uint64_t weight_roi;
    
    _weight weight;
    
    if( nrhs==6 )
    {
        weight.count = mxGetNumberOfElements(prhs[5]);
        weight.val = new int32_t[weight.count];
        
        double *pointer = mxGetPr(prhs[5]);
        
        for(size_t n=0; n<weight.count; n++)
        {
            if( mxGetClassID(prhs[5])==mxDOUBLE_CLASS )
                weight.val[n] = pointer[n];
            else if( mxGetClassID(prhs[5])==mxSINGLE_CLASS )
                weight.val[n] = *(reinterpret_cast<float *>(pointer)+n);
            else if( mxGetClassID(prhs[5])==mxINT64_CLASS )
                weight.val[n] = *(reinterpret_cast<int64_t *>(pointer)+n);
            else if( mxGetClassID(prhs[5])==mxINT32_CLASS )
                weight.val[n] = *(reinterpret_cast<int32_t *>(pointer)+n);
            else if( mxGetClassID(prhs[5])==mxINT16_CLASS )
                weight.val[n] = *(reinterpret_cast<int16_t *>(pointer)+n);
            else if( mxGetClassID(prhs[5])==mxINT8_CLASS )
                weight.val[n] = *(reinterpret_cast<int8_t *>(pointer)+n);
            else if( mxGetClassID(prhs[5])==mxUINT64_CLASS )
                weight.val[n] = *(reinterpret_cast<uint64_t *>(pointer)+n);
            else if( mxGetClassID(prhs[5])==mxUINT32_CLASS )
                weight.val[n] = *(reinterpret_cast<uint32_t *>(pointer)+n);
            else if( mxGetClassID(prhs[5])==mxUINT16_CLASS )
                weight.val[n] = *(reinterpret_cast<uint16_t *>(pointer)+n);
            else if( mxGetClassID(prhs[5])==mxUINT8_CLASS )
                weight.val[n] = *(reinterpret_cast<uint8_t *>(pointer)+n);
            else
                mexErrMsgTxt("Sixth input should be a number type");
            
            if( weight.val[n]<1 )
                mexErrMsgTxt("Invalid weight ROI");
        }
        
        if( weight.count==1 && weight.val[1]==1 )
        {
            delete [] weight.val;
            weight.count = 0;
            weight.val = nullptr;
        }
    }
    else
    {
        weight.count = 0;
        weight.val = nullptr;
    }

    fixedPoint Qstep = *mxGetPr(prhs[3]);

    double   *V = mxGetPr(prhs[0]);
    uint8_t  *C = reinterpret_cast<uint8_t *>(mxGetPr(prhs[1]));
    intmax_t *CT;
    char *binaryfilepath = mxArrayToString(prhs[2]);

// TRANSFORMADA HAAR3D E ESCREVE EM ARQUIVO BINÁRIO ----------------------------------------------
    file *fid = new file(binaryfilepath, 1);
    mxFree(binaryfilepath);

    if( fid->openError() ) {
        delete fid;
        mexErrMsgTxt("Unable to open file binaryfilepath");
    }

    uint64_t *W = new uint64_t[N];
    if( nrhs==6 && weight.count )
        index2weight(N, reinterpret_cast<uint8_t *>(mxGetPr(prhs[4])), weight, W);
    else
        for(size_t n=0; n<N; n++)
            W[n] = 1;

    CT = new intmax_t[3*N];

    haar3D(Qstep, V, C, W, N, CT);

    // write header
    fid->grWrite(N, 20);
    fid->grWrite(weight.count, 7);
    for(size_t n=0; n<weight.count; n++)
        fid->grWrite(weight.val[n]-1, 7);
    fid->write(&Qstep.val, sizeof(int64_t), 1);

    if( weight.count )
    {
        intmax_t *derivate = index_derivate(N, reinterpret_cast<uint8_t *>(mxGetPr(prhs[4])));
        fid->rlgrWrite(derivate, N);
        delete [] derivate;
    }
    fid->rlgrWrite(CT, N*3);
    delete [] CT;

// ----------------------------------------------------------------------------------------------
    if( nlhs ) {
        plhs[0] = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
        reinterpret_cast<uint64_t *>(mxGetPr(plhs[0]))[0] = fid->file_size();
    }
    delete fid;
    
    if( weight.count )
        delete [] weight.val;
}
