#include "mex.h"
#include "inteiros.h"
#include "haar3D.h"
#include <algorithm>
#include <math.h>
#include "file.h"

#define	debug	0

/* [CT, W] = haar3D_c(V,C,depth)
 * [CT, W] = haar3D_c(V,C,depth,filename,Qstep)
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
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
        mexPrintf("    C:          (Nx3 double matrix) colors associated to each voxel. This\n");
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
        mexPrintf("    CT:         (Nx3) transformed colors coefficients\n");
        return;
    }
    if( nrhs!=3 && nrhs!=5 )
        mexErrMsgTxt("Expected 3 or 5 inputs");
    if( nlhs>2 )
		mexErrMsgTxt("Expected 1 or 2 outputs");
    if( mxGetClassID(prhs[0])!=mxDOUBLE_CLASS ||
            mxGetClassID(prhs[1])!=mxDOUBLE_CLASS ||
            mxGetClassID(prhs[2])!=mxDOUBLE_CLASS )
        mexErrMsgTxt("First three inputs should be DOUBLE");
    if( mxGetN(prhs[0])!=3 || mxGetN(prhs[1])!=3 )
		mexErrMsgTxt("First two inputs (aka V and C) should have 3 columns");
    if( mxGetNumberOfElements(prhs[2])!=1 )
		mexErrMsgTxt("Third input (aka depth) shoud have 1 element");
    if( mxGetM(prhs[0])!=mxGetM(prhs[1]) )
		mexErrMsgTxt("First and second (aka V and C) input should have the same size");
    if( mxGetNumberOfDimensions(prhs[0])!=2 || mxGetNumberOfDimensions(prhs[1])!=2 )
		mexErrMsgTxt("First two inputs (aka V and C) should be bidimensional");
    if( ((int64) *mxGetPr(prhs[2]))<=0 )
		mexErrMsgTxt("Third input (aka depth) should be greater than 0");

	// Haar3d transform
    plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[0]), 3, mxREAL);
    if(nlhs==2)
    {
        plhs[1] = mxCreateDoubleMatrix(mxGetM(prhs[0]), 1, mxREAL);
        haar3D(mxGetPr(prhs[0]), mxGetPr(prhs[1]), mxGetM(prhs[0]), *mxGetPr(prhs[2]), mxGetPr(plhs[0]), mxGetPr(plhs[1]));
    }
    else
        haar3D(mxGetPr(prhs[0]), mxGetPr(prhs[1]), mxGetM(prhs[0]), *mxGetPr(prhs[2]), mxGetPr(plhs[0]));

	// Encode file
    if( nrhs==5 )
    {
        if( mxGetClassID(prhs[3])!=mxCHAR_CLASS )
			mexErrMsgTxt("Fourth input (aka filename) should be a char string");
        if( mxGetClassID(prhs[4])!=mxDOUBLE_CLASS )
            mexErrMsgTxt("Fifth input (aka Qstep) should be DOUBLE");
        if( mxGetNumberOfElements(prhs[4])!=1 )
            mexErrMsgTxt("Fifth input (aka Qstep) shoud have 1 element");

        double  Qstep = *mxGetPr(prhs[4]);

		if( Qstep<=0 )
			mexErrMsgTxt("Fifth input (aka Qstep) shoud be greater than 0");

        char    *filename = mxArrayToString(prhs[3]);
        file    *fid = new file(filename, WRITE_MODE);
        mxFree(filename);

        if( fid->error )
        {
            delete fid;
            mexErrMsgTxt("Unable to open file");
        }

        uint64  N = mxGetM(prhs[0])*3;
        int64   *data = (int64 *) malloc( N*sizeof(int64) );
        double  *CT = mxGetPr(plhs[0]);

        uint64  depth = *mxGetPr(prhs[2]);

		// write header
        fid->grWrite(mxGetM(prhs[0]), 20);
        fid->grWrite(depth-1, 4);
        fid->write(&Qstep, sizeof(double), 1);

        while( N-- )
            data[N] = round(CT[N]/Qstep);

        fid->rlgrStart();
        fid->rlgrWrite(data, mxGetM(prhs[0])*3);
        fid->rlgrStop();

        free(data);

		uint8	error = fid->error;
		delete fid;
		switch (error)
		{
		case 0:
			break;
		case 4:
			mexErrMsgTxt("Overflow: tried to write more than 64 bits at once");
		default:
			mexErrMsgIdAndTxt("RAHT::file", "Unexpected error. Error number %d", error);
		}
    }
}
