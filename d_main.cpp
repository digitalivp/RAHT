#include "mex.h"
#include <stdint.h>
#include "haar3D.h"
#include <algorithm>
#include <math.h>
#include "file.h"

#define	debug	0

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    fixedPoint  Qstep;
    size_t      N;
    size_t      depth;

    // Test inputs
    if( !nrhs && !nlhs )
    {
        mexPrintf("RAHT_cod    Region Adaptive Hierarchical Transform coder\n\n");
        mexPrintf("    CT = RAHT_cod(V, C, depth) computes the transformed coefficients. The\n");
        mexPrintf("    coefficients are not quantized and not coded to any file\n\n");
        mexPrintf("    CT = RAHT_cod(V, C, depth, filename, Qstep) transforms the coefficients\n");
        mexPrintf("    and write the coded data to the file given by 'filename' with the given\n");
        mexPrintf("    quantization step. The output is optional and is not quantized. Only\n");
        mexPrintf("    the coefficients writen to file are quantized\n\n");
        mexPrintf("INPUTS\n");
        mexPrintf("    V:          (Nx3 double matrix) vertices of each voxel in the order\n");
        mexPrintf("                (X, Y, Z). The values should be non negative integers\n\n");
		mexPrintf("    C:          (NxK double matrix) colors associated to each voxel. This\n");
        mexPrintf("                input is transparente to RGB, YUV, YCbCr or any other color\n");
        mexPrintf("                space to be used as long as it is a triplet. It is also\n");
        mexPrintf("                transparent to the color range: colors can be in the range\n");
        mexPrintf("                0--1, 0--255, 0--65535 or any other range\n\n");
        mexPrintf("    depth:      (1x1 double value) the depth of the octree to be used for\n");
        mexPrintf("                the transform. It needs to be an integer value grater than 0\n\n");
        mexPrintf("    filename:   (string of chars) file name\n\n");
        mexPrintf("    Qstep:      (1x1 double value) quantization step. It is dependent of\n");
        mexPrintf("                the color range\n\n");
        mexPrintf("OUTPUTS\n");
		mexPrintf("    CT:         (NxK) transformed colors coefficients\n");
        return;
    }
    if( nrhs!=5 )
        mexErrMsgTxt("Expected 5 inputs");
    if( nlhs )
        mexErrMsgTxt("Expected no outputs");
    if( mxGetClassID(prhs[0])!=mxDOUBLE_CLASS ||
            mxGetClassID(prhs[1])!=mxDOUBLE_CLASS ||
            mxGetClassID(prhs[2])!=mxDOUBLE_CLASS )
        mexErrMsgTxt("First three inputs should be DOUBLE");
	if( mxGetN(prhs[0])!=3 )
		mexErrMsgTxt("First input (aka V) should have 3 columns");
	if( !mxGetN(prhs[1]) )
		mexErrMsgTxt("Second input (aka C) should have at least 1 column");
    if( mxGetNumberOfElements(prhs[2])!=1 )
		mexErrMsgTxt("Third input (aka depth) shoud have 1 element");
    if( mxGetM(prhs[0])!=mxGetM(prhs[1]) )
		mexErrMsgTxt("First and second (aka V and C) input should have the same number of lines");
    if( mxGetNumberOfDimensions(prhs[0])!=2 || mxGetNumberOfDimensions(prhs[1])!=2 )
		mexErrMsgTxt("First two inputs (aka V and C) should be bidimensional");
	if( ((int64_t) *mxGetPr(prhs[2]))<=0 )
		mexErrMsgTxt("Third input (aka depth) should be greater than 0");
    if( mxGetClassID(prhs[3])!=mxCHAR_CLASS )
        mexErrMsgTxt("Fourth input (aka filename) should be a char string");
    if( mxGetClassID(prhs[4])!=mxDOUBLE_CLASS )
        mexErrMsgTxt("Fifth input (aka Qstep) should be DOUBLE");
    if( mxGetNumberOfElements(prhs[4])!=1 )
        mexErrMsgTxt("Fifth input (aka Qstep) shoud have 1 element");

    Qstep = *mxGetPr(prhs[4]);

    if( Qstep.toDouble()<0 )
        mexErrMsgTxt("Fifth input (aka Qstep) shoud be greater than 0");

    N = mxGetM(prhs[1])*mxGetN(prhs[1]);
    depth = *mxGetPr(prhs[2]);
    int64_t    *data = (int64_t *) malloc( N*sizeof(int64_t) );

	// Haar3d transform
    haar3D(Qstep, mxGetPr(prhs[0]), mxGetPr(prhs[1]), mxGetN(prhs[1]), mxGetM(prhs[1]), depth, data);

	// Encode file
    char    *filename = mxArrayToString(prhs[3]);
    file    *fid = new file(filename, 1);
    mxFree(filename);

    if( fid->openError() )
    {
        delete fid;
        mexErrMsgTxt("Unable to open file");
    }

    // write header
    fid->grWrite(mxGetM(prhs[1]), 20);
    fid->grWrite(mxGetN(prhs[1]), 3);
    fid->grWrite(depth-1, 4);
    fid->write(Qstep.val, 64);

    fid->rlgrWrite(data, mxGetM(prhs[1])*mxGetN(prhs[1]));
    delete fid;
    free(data);
}
