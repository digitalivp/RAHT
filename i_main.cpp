#include "mex.h"
#include <stdint.h>
#include "haar3D.h"
#include <algorithm>
#include <math.h>
#include "file.h"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    fixedPoint  Qstep;
    size_t      N, K;
    size_t      depth;
    int64_t    *CT;

	// Test inputs
	if( !nrhs && !nlhs )
	{
        mexPrintf("RAHT_dec    Region Adaptive Hierarchical Transform decoder\n\n");
        mexPrintf("    C = RAHT_dec(V, CT, depth) performs the inverse transform from the given\n");
        mexPrintf("    coefficients\n\n");
        mexPrintf("    C = RAHT_dec(V, filename) decodes the data from the given file\n\n");
        mexPrintf("INPUTS\n");
        mexPrintf("    V:          (Nx3 double matrix) vertices of each voxel in the order\n");
        mexPrintf("                (X, Y, Z). The values should be non negative integers\n\n");
		mexPrintf("    CT:         (NxK double matrix) transformed coefficients coefficients\n\n");
        mexPrintf("    depth:      (1x1 double value) the depth of the octree to be used for\n");
        mexPrintf("                the transform. It needs to be an integer value grater than 0\n\n");
        mexPrintf("    filename:   (string of chars) file name\n\n");
        mexPrintf("OUTPUTS\n");
		mexPrintf("    C:          (NxK) colors associated to each voxel\n");
        return;
	}
    if( nrhs!=2 )
        mexErrMsgTxt("Expected 2 inputs");
	if( nlhs>1 )
		mexErrMsgTxt("Expected 1 output");
    if( mxGetClassID(prhs[0])!=mxDOUBLE_CLASS )
        mexErrMsgTxt("First input should be DOUBLE");
    if( mxGetN(prhs[0])!=3 )
        mexErrMsgTxt("First input (aka V) should have 3 columns");
    if( mxGetNumberOfDimensions(prhs[0])!=2 )
        mexErrMsgTxt("First input (aka V) should be bidimensional");
    if( mxGetClassID(prhs[1])!=mxCHAR_CLASS )
        mexErrMsgTxt("Second input (aka filename) should be a char string");

    char    *filename = mxArrayToString(prhs[1]);
    file    *fid = new file(filename, 0);
    mxFree(filename);

    if( fid->openError() )
    {
        delete fid;
        mexErrMsgTxt("Unable to open file");
    }

    N = fid->grRead(20);
    K = fid->grRead(3);
    depth = fid->grRead(4)+1;
    Qstep.val = fid->read(64);

    if( mxGetM(prhs[0])!=N || depth<=0 )
        mexErrMsgTxt("First input (aka V) and the given file seem to be incompatible");

    CT = (int64_t *) malloc( N*K*sizeof(int64_t) );

    fid->rlgrRead(CT, N*K);
    if( !fid->eof() )
        mexWarnMsgTxt("The given file seens to be incompatible with the given geometry");
    delete fid;

	// Inverse transform
    plhs[0] = mxCreateDoubleMatrix(N, K, mxREAL);
    inv_haar3D(Qstep, mxGetPr(prhs[0]), CT, K, N, depth, mxGetPr(plhs[0]));

    if( nrhs==2 )
        free(CT);
}
