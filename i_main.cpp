#include "mex.h"
#include "inteiros.h"
#include "haar3D.h"
#include <algorithm>
#include <math.h>
#include "file.h"

/* C = inv_haar3D_c(V,CT,depth)
 * C = inv_haar3D_c(V,filename)
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
	double	*CT;
	uint64	depth;

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
        mexPrintf("    CT:         (Nx3 double matrix) transformed coefficients coefficients\n\n");
        mexPrintf("    depth:      (1x1 double value) the depth of the octree to be used for\n");
        mexPrintf("                the transform. It needs to be an integer value grater than 0\n\n");
        mexPrintf("    filename:   (string of chars) file name\n\n");
        mexPrintf("OUTPUTS\n");
        mexPrintf("    C:          (Nx3) colors associated to each voxel\n");
        return;
	}
	if( nlhs>1 )
		mexErrMsgTxt("Expected 1 output");
	if( mxGetClassID(prhs[0])!=mxDOUBLE_CLASS )
		mexErrMsgTxt("First input should be DOUBLE");
	if( nrhs==3 )
	{
		if( mxGetClassID(prhs[1])!=mxDOUBLE_CLASS || mxGetClassID(prhs[2])!=mxDOUBLE_CLASS)
			mexErrMsgTxt("Inputs should be DOUBLE");
		if( mxGetN(prhs[0])!=3 || mxGetN(prhs[1])!=3 )
			mexErrMsgTxt("First two inputs (aka V and C) should have 3 columns");
		if( mxGetNumberOfElements(prhs[2])!=1 )
			mexErrMsgTxt("Third input (aka depth) shoud have 1 element");
		if( mxGetM(prhs[0])!=mxGetM(prhs[1]) )
			mexErrMsgTxt("First and second input (aka V and C) should have the same size");
		if( mxGetNumberOfDimensions(prhs[0])!=2 || mxGetNumberOfDimensions(prhs[1])!=2 )
			mexErrMsgTxt("First two inputs (aka V and C) should be bidimensional");

		CT = mxGetPr(prhs[1]);
		depth = *mxGetPr(prhs[2]);

		if( depth<=0 )
			mexErrMsgTxt("Third input (aka depth) should be greater than 0");
	}
	else if( nrhs==2 )
	{
		if( mxGetN(prhs[0])!=3 )
			mexErrMsgTxt("First input (aka V) should have 3 columns");
		if( mxGetNumberOfDimensions(prhs[0])!=2 )
			mexErrMsgTxt("First input (aka V) should be bidimensional");
		if( mxGetClassID(prhs[1])!=mxCHAR_CLASS )
			mexErrMsgTxt("Second input (aka filename) should be a char string");

		double  Qstep;
		char    *filename = mxArrayToString(prhs[1]);
		file    *fid = new file(filename, READ_MODE);
		mxFree(filename);

		if( fid->error )
		{
			delete fid;
			mexErrMsgTxt("Unable to open file");
		}

		uint64	N = fid->grRead(20);
        depth = fid->grRead(4)+1;
		fid->read(&Qstep, sizeof(double), 1);

		if( fid->error )
		{
			delete fid;
			mexErrMsgTxt("Unable to read file");
		}

		if( mxGetM(prhs[0])!=N || depth<=0 )
			mexErrMsgTxt("First input (aka V) and the given file seem to be incompatible");

		N *= 3;
		int64	*data = (int64 *) malloc( N*sizeof(int64) );
		CT = (double *) malloc( N*sizeof(double) );

		fid->rlgrStart();
		fid->rlgrRead(data, N);
		fid->rlgrStop();

		while( N-- )
			CT[N] = Qstep*data[N];

		free(data);

		uint8	error = fid->error;
		delete fid;
		switch (error)
		{
		case 0:
			break;
		case 3:
			mexErrMsgTxt("Reached end of file before reading all necessary bits");
		case 4:
			mexErrMsgTxt("Overflow: tried to read more than 64 bits at once");
		default:
			mexErrMsgIdAndTxt("RAHT::file", "Unexpected error. Error number %d", error);
		}
	}
	else
		mexErrMsgTxt("Expected 3 inputs");

    //
	plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[0]), 3, mxREAL);
	inv_haar3D(mxGetPr(prhs[0]), CT, mxGetM(prhs[0]), depth, mxGetPr(plhs[0]));

	if( nrhs==2 )
		free(CT);
}
