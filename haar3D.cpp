#include "haar3D.h"

struct mem
{
    uint64 val;
    uint64 ord;
};

bool compari(mem &A, mem &B) {return A.val < B.val;}

void copyAsort(double *VX, double *CX, uint64 N, vali *C, uint64 *W, uint64 *val, uint64 *ord)
{
	mem		*MEM = (mem *) malloc( N*sizeof(mem) );
	double	*VY = VX+N, *VZ = VY+N;
	double	*CY = CX+N, *CZ = CY+N;
	uint64	vx, vy, vz;

	for(uint64 i=0; i<N; i++)
	{
		vx = VZ[i];
		vy = VY[i];
		vz = VX[i];

		MEM[i].val =
				((0x000001 & vx)    ) + ((0x000001 & vy)<< 1) + ((0x000001 & vz)<< 2) +
				((0x000002 & vx)<< 2) + ((0x000002 & vy)<< 3) + ((0x000002 & vz)<< 4) +
				((0x000004 & vx)<< 4) + ((0x000004 & vy)<< 5) + ((0x000004 & vz)<< 6) +
				((0x000008 & vx)<< 6) + ((0x000008 & vy)<< 7) + ((0x000008 & vz)<< 8) +
				((0x000010 & vx)<< 8) + ((0x000010 & vy)<< 9) + ((0x000010 & vz)<<10) +
				((0x000020 & vx)<<10) + ((0x000020 & vy)<<11) + ((0x000020 & vz)<<12) +
				((0x000040 & vx)<<12) + ((0x000040 & vy)<<13) + ((0x000040 & vz)<<14) +
				((0x000080 & vx)<<14) + ((0x000080 & vy)<<15) + ((0x000080 & vz)<<16) +
				((0x000100 & vx)<<16) + ((0x000100 & vy)<<17) + ((0x000100 & vz)<<18) +
				((0x000200 & vx)<<18) + ((0x000200 & vy)<<19) + ((0x000200 & vz)<<20) +
				((0x000400 & vx)<<20) + ((0x000400 & vy)<<21) + ((0x000400 & vz)<<22) +
				((0x000800 & vx)<<22) + ((0x000800 & vy)<<23) + ((0x000800 & vz)<<24) +
				((0x001000 & vx)<<24) + ((0x001000 & vy)<<25) + ((0x001000 & vz)<<26) +
				((0x002000 & vx)<<26) + ((0x002000 & vy)<<27) + ((0x002000 & vz)<<28) +
				((0x004000 & vx)<<28) + ((0x004000 & vy)<<29) + ((0x004000 & vz)<<30) +
				((0x008000 & vx)<<30) + ((0x008000 & vy)<<31) + ((0x008000 & vz)<<32) +
				((0x010000 & vx)<<32) + ((0x010000 & vy)<<33) + ((0x010000 & vz)<<34) +
				((0x020000 & vx)<<34) + ((0x020000 & vy)<<35) + ((0x020000 & vz)<<36) +
				((0x040000 & vx)<<36) + ((0x040000 & vy)<<37) + ((0x040000 & vz)<<38) +
				((0x080000 & vx)<<38) + ((0x080000 & vy)<<39) + ((0x080000 & vz)<<40) +
				((0x100000 & vx)<<40) + ((0x100000 & vy)<<41) + ((0x100000 & vz)<<42);
		MEM[i].ord = i;

		C->cx  = CX[i];
		C->cy  = CY[i];
		C->cz  = CZ[i];
		C++;

		W[i] = 1;
	}
	//

	std::sort(MEM, MEM+N, compari);

	for(uint64 i=0; i<N; i++)
	{
		val[i] = MEM[i].val;
		ord[i] = MEM[i].ord;
	}

	free(MEM);
}

void transform(double a0, double a1, vali *C0, vali *C1, vali *CT0, vali *CT1)
{
	CT0->cx = a0*C0->cx + a1*C1->cx;
	CT0->cy = a0*C0->cy + a1*C1->cy;
	CT0->cz = a0*C0->cz + a1*C1->cz;

	CT1->cx = a0*C1->cx - a1*C0->cx;
	CT1->cy = a0*C1->cy - a1*C0->cy;
	CT1->cz = a0*C1->cz - a1*C0->cz;
}

void itransform(double a0, double a1, vali *C0, vali *C1, vali *CT0, vali *CT1)
{
	C0->cx = a0*CT0->cx - a1*CT1->cx;
	C0->cy = a0*CT0->cy - a1*CT1->cy;
	C0->cz = a0*CT0->cz - a1*CT1->cz;

	C1->cx = a0*CT1->cx + a1*CT0->cx;
	C1->cy = a0*CT1->cy + a1*CT0->cy;
	C1->cz = a0*CT1->cz + a1*CT0->cz;
}

void copyFromMEM(uint64 *IN_VAL, uint64 *IN_W, uint64 *OUT_VAL, uint64 *OUT_W, uint64 M)
{
	while(M)
	{
		M--;
		IN_VAL[M] = OUT_VAL[M];
		IN_W[M]   = OUT_W[M];
	}
}

/* uint64	N = mxGetM(prhs[0]);
 * uint64	depth = *mxGetPr(prhs[2])
 * double	*inV = mxGetPr(prhs[0]);
 * double	*inC = mxGetPr(prhs[1]);
 * plhs[0] = mxCreateDoubleMatrix(NN, 3, mxREAL);
 * double	*outCT = mxGetPr(plhs[0])
 */
void haar3D(double *inV, double *inC, uint64 N, uint64 depth, double *outCT, double *outW)
{
	uint64	NN=N;
	uint64	M=N, S, d, i, j;
	double	a;
	depth *= 3;

	vali	*C	  = (vali    *) malloc( N*sizeof(vali) );
	vali	*CT   = (vali    *) malloc( N*sizeof(vali) );
	uint64	*w	  = (uint64  *) malloc( N*sizeof(uint64) );
	uint64	*wT   = (uint64  *) malloc( N*sizeof(uint64) );
	uint64	*val  = (uint64  *) malloc( N*sizeof(uint64) );
	uint64	*valT = (uint64  *) malloc( N*sizeof(uint64) );
	uint64	*TMP  = (uint64  *) malloc( N*sizeof(uint64) );
	vali    *TPV;

	copyAsort(inV, inC, N, CT, w, val, TMP);
	for(i=0; i<N; i++)
		C[i] = CT[TMP[i]];
	free(TMP);

	for(d=0; d<depth; d++)
	{
		i = 0;
		M = 0;
		S = N;

		while( i<S )
		{
			j = i+1;
			valT[M] = val[i] >> 1;

			if( j<S && ((val[i]&0xFFFFFFFFFFFFFFFE)==(val[j]&0xFFFFFFFFFFFFFFFE)) )
			{
				N--;

				wT[M] = w[i]+w[j];
				wT[N] = wT[M];

				a = sqrt(wT[M]);
				transform(sqrt(w[i])/a, sqrt(w[j])/a, C+i, C+j, CT+M, CT+N);

				i += 2;
			}
			else
			{
				wT[M] = w[i];
				CT[M] = C[i];

				i += 1;
			}
			M++;
		}
		for(i=N; i<S; i++)
		{
			C[i] = CT[i];
			w[i] = wT[i];
		}

		TPV  = CT;
		CT   = C;
		C    = TPV;

		TMP  = valT;
		valT = val;
		val  = TMP;

		TMP  = wT;
		wT   = w;
		w    = TMP;
	}

	free(wT);
	free(val);
	free(valT);

	double *CX = outCT, *CY = CX+NN, *CZ = CY+NN;

	if( outW!=NULL )
	{
		double *WX = outW;

		C += NN;
		i = NN;
		while( i )
		{
			C--;
			i--;
			CX[i] = C->cx;
			CY[i] = C->cy;
			CZ[i] = C->cz;
			WX[i] = w[i];
		}
	}
	else
	{
		C += NN;
		i = NN;
		while( i )
		{
			C--;
			i--;
			CX[i] = C->cx;
			CY[i] = C->cy;
			CZ[i] = C->cz;
		}
	}

	free(w);
	free(C);
	free(CT);
}

/* uint64	N = mxGetM(prhs[0])
 * uint64	depth = *mxGetPr(prhs[2]);
 * double	*inV = mxGetPr(prhs[0]);
 * double	*inCT = mxGetPr(prhs[1]);
 * plhs[0] = mxCreateDoubleMatrix(NN, 3, mxREAL);
 * double	*outC = mxGetPr(plhs[0])
 */
void inv_haar3D(double *inV, double *inCT, uint64 N, uint64 depth, double *outC)
{
	uint64	NN=N;
	uint64	M=N, S, d, i, j;
	double	a;
	depth *= 3;

	vali	*C	  = (vali    *) malloc( N*sizeof(vali) );
	vali	*CT   = (vali    *) malloc( N*sizeof(vali) );
	uint64	*w	  = (uint64  *) malloc( N*sizeof(uint64) );
	uint64	*wT   = (uint64  *) malloc( N*sizeof(uint64) );
	uint64	*val  = (uint64  *) malloc( N*sizeof(uint64) );
	uint64	*valT = (uint64  *) malloc( N*sizeof(uint64) );
	uint64	*ord  = (uint64  *) malloc( N*sizeof(uint64) );

	uint64	**VAL = (uint64 **) malloc( depth*sizeof(uint64 *) );
	uint64	**iW  = (uint64 **) malloc( depth*sizeof(uint64 *) );
	uint64	*iM   = (uint64  *) malloc( depth*sizeof(uint64) );

	uint64  *TMP;

	copyAsort(inV, inCT, N, CT, w, val, ord);
	for(i=0; i<N; i++)
		C[i] = CT[i];

	// Transformada direta (partes)
	// executa a mesma ordem de passos da transformada direta, apenas armazenando
	// os valores do pesos e do vetor de valor. Estes valores são armazenados para
	// executar os passos em ordem inversa na transformada inversa
	for(d=0; d<depth; d++)
	{
		VAL[d] = (uint64 *) malloc( M*sizeof(uint64) );
		iW[d]  = (uint64 *) malloc( M*sizeof(uint64) );

		copyFromMEM(VAL[d], iW[d], val, w, M);
		iM[d] = M;

		i = 0;
		S = M;
		M = 0;

		while( i<S )
		{
			valT[M] = val[i] >> 1;
			j = i+1;

			if( j<S && ((val[i]&0xFFFFFFFFFFFFFFFE)==(val[j]&0xFFFFFFFFFFFFFFFE)) )
			{
				wT[M] = w[i]+w[j];
				i += 2;
			}
			else
			{
				wT[M] = w[i];
				i += 1;
			}
			M++;
		}

		TMP  = valT;
		valT = val;
		val  = TMP;

		TMP = wT;
		wT  = w;
		w   = TMP;
	}

	// Transformada inversa
	while( d )
	{
		d--;
		S = iM[d];
		M = d?iM[d-1]:NN;

		copyFromMEM(val, w, VAL[d], iW[d], S);
		for(i=S; i<M; i++)
			C[i] = CT[i];

		free(VAL[d]);
		free(iW[d]);

		M = 0;
		N = S;
		i = 0;

		while( i<S )
		{
			j = i+1;

			if( j<S && ((val[i]&0xFFFFFFFFFFFFFFFE)==(val[j]&0xFFFFFFFFFFFFFFFE)) )
			{
				a  = sqrt(w[i]+w[j]);

				N--;
				itransform(sqrt(w[i])/a, sqrt(w[j])/a, C+i, C+j, CT+M, CT+N);
				i += 2;
			}
			else
			{
				C[i] = CT[M];
				i += 1;
			}
			M++;
		}

		for(i=0; i<S; i++)
			CT[i] = C[i];
	}

	// Copia dados para saida
	free(iM);
	free(iW);
	free(VAL);
	free(valT);
	free(val);
	free(wT);
	free(w);
	free(CT);

	double *CX = outC, *CY = CX+NN, *CZ = CY+NN;

	for(i=0; i<NN; i++)
	{
		CX[ord[i]] = C[i].cx;
		CY[ord[i]] = C[i].cy;
		CZ[ord[i]] = C[i].cz;
	}

	free(C);
	free(ord);
}
