#include "haar3D.h"

struct mem
{
    uint64_t val;
    uint64_t ord;
};

bool compari(const mem &A, const mem &B) {return A.val < B.val;}

void copyAsort(double *V, size_t N, uint64_t *val, uint64_t *ord)
{
    mem			*MEM = new mem[N];
    double      *VX = V;
    double      *VY = VX+N;
    double      *VZ = VY+N;
    uint64_t	vx, vy, vz;

    for(size_t i=0; i<N; i++)
    {
        vz = static_cast<uint64_t>(*(VX++));
        vy = static_cast<uint64_t>(*(VY++));
        vx = static_cast<uint64_t>(*(VZ++));

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
    }

    std::sort(MEM, MEM+N, compari);

    for(size_t i=0; i<N; i++)
    {
        val[i] = MEM[i].val;
        ord[i] = MEM[i].ord;
    }

    delete [] MEM;
}

void transform(fixedPoint Qstep, uint64_t w0, uint64_t w1, fixedPoint *C0, fixedPoint *C1, fixedPoint *CT0, fixedPoint *CT1, size_t K)
{
    fixedPoint  b;
    b.val = static_cast<int64_t>((w1<<_fixedpoint_PRECISION)/(w0+w1));
    Qstep.val = _sqrt( (static_cast<uint64_t>(Qstep.val*Qstep.val)*(w0+w1))/(w0*w1) );

    while( K-- )
    {
        CT1->val  = C1->val;
        CT1->val -= C0->val;
        CT0->val  = CT1->val;
        CT0->val *= b.val;
        if( CT0->val < 0 )
            CT0->val = -( (_fixedpoint_RND - CT0->val) >> _fixedpoint_PRECISION );
        else
            CT0->val = +( (_fixedpoint_RND + CT0->val) >> _fixedpoint_PRECISION );
        CT0->val += C0->val;
        if( CT1->val < 0 )
            CT1->val = -(
                        ((+Qstep.val)>>1) +
                        ((-CT1->val)<<_fixedpoint_PRECISION) )/Qstep.val;
        else
            CT1->val = +(
                        ((+Qstep.val)>>1) +
                        ((+CT1->val)<<_fixedpoint_PRECISION) )/Qstep.val;

        C0++;
        C1++;
        CT0++;
        CT1++;
    }
}

void itransform(fixedPoint Qstep, uint64_t w0, uint64_t w1, fixedPoint *C0, fixedPoint *C1, fixedPoint *CT0, fixedPoint *CT1, size_t K)
{
    fixedPoint  b;
    b.val = static_cast<int64_t>((w1<<_fixedpoint_PRECISION)/(w0+w1));
    Qstep.val = _sqrt( (static_cast<uint64_t>(Qstep.val*Qstep.val)*(w0+w1))/(w0*w1) );

    while( K-- )
    {
        C0->val  = CT1->val;
        C0->val *= Qstep.val;
        if( C0->val < 0 )
            C0->val = -( (_fixedpoint_RND - C0->val) >> _fixedpoint_PRECISION );
        else
            C0->val = +( (_fixedpoint_RND + C0->val) >> _fixedpoint_PRECISION );
        C1->val  = C0->val;
        C0->val *= b.val;
        if( C0->val < 0 )
            C0->val = -( (_fixedpoint_RND - C0->val) >> _fixedpoint_PRECISION );
        else
            C0->val = +( (_fixedpoint_RND + C0->val) >> _fixedpoint_PRECISION );
        C0->val -= CT0->val;
        C0->val  = -C0->val;
        C1->val += C0->val;

        C0++;
        C1++;
        CT0++;
        CT1++;
    }
}

void copyFromMEM(uint64_t *IN_VAL, uint64_t *IN_W, uint64_t *OUT_VAL, uint64_t *OUT_W, uint64_t M)
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
//void haar3D(double *inV, double *inC, size_t K, size_t N, size_t depth, double *outCT, double *outW)
void haar3D(fixedPoint Qstep, double *inV, uint8_t *inC, uint64_t *wT, size_t N, intmax_t *outCT)
{
    size_t	NN=N;
    size_t	M=N, S, d, i, j;

    fixedPoint	*C	  = new fixedPoint[N*3];
    fixedPoint	*CT   = new fixedPoint[N*3];
    uint64_t	*w	  = new uint64_t[N];
    uint64_t	*val  = new uint64_t[N];
    uint64_t	*valT = new uint64_t[N];
    uint64_t	*TMP  = new u_int64_t[N];
    fixedPoint	*TPV;

    copyAsort(inV, N, val, TMP);
    for(i=0; i<N; i++)
    {
        uint8_t *r = inC + TMP[i];
        uint8_t *g = r+N;
        uint8_t *b = g+N;

        double  y;

        w[i] = wT[TMP[i]];

        y = 0.212600*static_cast<double>(*r) + 0.715200*static_cast<double>(*g) + 0.072200*static_cast<double>(*b);
        C[3*i+2] = (static_cast<double>(*r)-y) / 1.57480;
        C[3*i+1] = (static_cast<double>(*b)-y) / 1.85563;
        C[3*i+0] = y;
    }
    delete [] TMP;

    size_t depth = 0;
    {
        uint64_t maxval = val[N-1];
        while( maxval )
        {
            depth += 3;
            maxval >>= 3;
        }
    }

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

                transform(Qstep, w[i], w[j], C+i*3, C+j*3, CT+M*3, CT+N*3, 3);

                i += 2;
            }
            else
            {
                wT[M] = w[i];
                for(size_t k=0; k<3; k++)
                    CT[M*3+k] = C[i*3+k];

                i += 1;
            }
            M++;
        }
        for(i=N; i<S; i++)
        {
            for(size_t k=0; k<3; k++)
                C[i*3+k] = CT[i*3+k];
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

    delete [] wT;
    delete [] val;
    delete [] valT;

    // Quantization of DC coefficients
    Qstep.val = _sqrt( (Qstep.val*Qstep.val)/w[0] );
    for(size_t k=0; k<3; k++)
        C[k] /= Qstep;

    for(i=0; i<NN; i++)
        for(size_t k=0; k<3; k++)
            outCT[i+k*NN] = C[i*3+k].round();

    delete [] w;
    delete [] C;
    delete [] CT;
}

/* uint64	N = mxGetM(prhs[0])
 * uint64	depth = *mxGetPr(prhs[2]);
 * double	*inV = mxGetPr(prhs[0]);
 * double	*inCT = mxGetPr(prhs[1]);
 * plhs[0] = mxCreateDoubleMatrix(NN, 3, mxREAL);
 * double	*outC = mxGetPr(plhs[0])
 */
//void inv_haar3D(double *inV, double *inCT, size_t K, size_t N, size_t depth, double *outC)
void inv_haar3D(fixedPoint Qstep, double *inV, intmax_t *inCT, uint64_t *wT, size_t N, uint8_t *outC)
{
    size_t	NN=N;
    size_t	M=N, S, d, i, j;

    fixedPoint	*C	  = new fixedPoint[N*3];
    fixedPoint	*CT   = new fixedPoint[N*3];
    uint64_t	*w	  = new uint64_t[N];
    uint64_t	*val  = new uint64_t[N];
    uint64_t	*valT = new uint64_t[N];
    uint64_t	*ord  = new uint64_t[N];

    uint64_t    *TMP;

    copyAsort(inV, N, val, ord);
    for(i=0; i<N; i++)
    {
        w[i] = wT[ord[i]];
        for(size_t k=0; k<3; k++)
            CT[i*3+k] = static_cast<double>(inCT[i+k*N]);
    }

    size_t depth = 0;
    {
        uint64_t maxval = val[N-1];
        while( maxval )
        {
            depth += 3;
            maxval >>= 3;
        }
    }

    uint64_t	**VAL = new uint64_t*[depth];
    uint64_t	**iW  = new uint64_t*[depth];
    size_t		*iM   = new size_t[depth];

    // Dequantization of DC coefficients
    {
        fixedPoint Qstep2;
        uint64_t W0 = 0;
        for(size_t n=0; n<N; n++)
            W0 += w[n];
        Qstep2.val = _sqrt( (Qstep.val*Qstep.val)/W0 );
        for(size_t k=0; k<3; k++)
            CT[k] *= Qstep2;
    }

    // Transformada direta (partes)
    // executa a mesma ordem de passos da transformada direta, apenas armazenando
    // os valores do pesos e do vetor de valor. Estes valores são armazenados para
    // executar os passos em ordem inversa na transformada inversa
    for(d=0; d<depth; d++)
    {
        VAL[d] = new uint64_t[M];
        iW[d]  = new uint64_t[M];

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
            for(size_t k=0; k<3; k++)
                C[i*3+k] = CT[i*k+3];

        delete [] VAL[d];
        delete [] iW[d];

        M = 0;
        N = S;
        i = 0;

        while( i<S )
        {
            j = i+1;

            if( j<S && ((val[i]&0xFFFFFFFFFFFFFFFE)==(val[j]&0xFFFFFFFFFFFFFFFE)) )
            {
                N--;
                itransform(Qstep, w[i], w[j], C+i*3, C+j*3, CT+M*3, CT+N*3, 3);
                i += 2;
            }
            else
            {
                for(size_t k=0; k<3; k++)
                    C[i*3+k] = CT[M*3+k];
                i += 1;
            }
            M++;
        }

        for(i=0; i<S; i++)
            for(size_t k=0; k<3; k++)
                CT[i*3+k] = C[i*3+k];
    }

    // Copia dados para saida
    delete [] iM;
    delete [] iW;
    delete [] VAL;
    delete [] valT;
    delete [] val;
    delete [] wT;
    delete [] w;
    delete [] CT;

    //double *CX = outC, *CY = CX+NN, *CZ = CY+NN;

    for(i=0; i<NN; i++)
    {
        uint8_t *color = outC + ord[i];

        double  y = C[3*i+0].toDouble()+0.5;
        double  u = C[3*i+1].toDouble();
        double  v = C[3*i+2].toDouble();

        double  r = y + 1.57480*v;
        double  g = y - 0.18733*u - 0.46813*v;
        double  b = y + 1.85563*u;

        // Clipping
        if( r<0 )
            *color = 0;
        else if( r>255 )
            *color = 255;
        else
            *color = static_cast<uint8_t>(r);
        color += NN;

        if( g<0 )
            *color = 0;
        else if( g>255 )
            *color = 255;
        else
            *color = static_cast<uint8_t>(g);
        color += NN;

        if( b<0 )
            *color = 0;
        else if( b>255 )
            *color = 255;
        else
            *color = static_cast<uint8_t>(b);
    }

    delete [] C;
    delete [] ord;
}

intmax_t *index_derivate(size_t N, uint8_t *index)
{
    intmax_t *derivate = new intmax_t[N];

    derivate[0] = index[0];
    for(size_t n=1; n<N; n++)
        derivate[n] = static_cast<intmax_t>(index[n]) - static_cast<intmax_t>(index[n-1]);

    return derivate;
}

uint8_t *index_integrate(size_t N, intmax_t *derivate)
{
    uint8_t *index = new uint8_t[N];
    
    index[0] = derivate[0];
    for(size_t n=1; n<N; n++)
        index[n] = derivate[n] + index[n-1];
    
    return index;
}

void index2weight(size_t N, uint8_t *index, _weight weight, uint64_t *w)
{
    for(size_t n=0; n<N; n++)
    {
        if( index[n] )
            w[n] = weight.val[index[n]-1];
        else
            w[n] = 1;
    }
}