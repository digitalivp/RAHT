#include "file.h"

#if PROB_BEHAVIOR_RETURN_CDF

/* This function returns the result of
 *	  L*(C1-C0)
 *   -----------
 *    2^64-1-C0
 *
 * It is supposed to be used internally in the code
 */
uint64 _file_multiply(uint64 L, uint64 C1, uint64 C0)
{
	uint64	a1 = C1-C0;
	uint64	a0 = a1 & 0xFFFFFFFF;
	uint64	b1 = L  & 0xFFFFFFFF;
	uint64	b0 = L >> 32;
	a1 >>= 32;

	uint64	ca = a1*b0;
	uint64	cb = a0*b1;
	uint64	d  = a0*b0;

	uint64	low  = (ca&0xFFFFFFFF) + (cb&0xFFFFFFFF) + ((a0*b0)>>32);
	uint64	high = (a1*b1) + (ca>>32) + (cb>>32) + (low>>32);
	low = (low<<32) + (d&0xFFFFFFFF);

	uint64	div = MAXUINT64-C0;
	uint64	q = low<<1;
	uint64	rem = high;
	uint64	carry = low>>63;
	uint64	temp_carry = 0;

	for(int i=0; i<64; i++)
	{
		temp_carry = rem >> 63;
		rem <<= 1;
		rem |= carry;
		carry = temp_carry;

		if( !carry )
		{
			if( rem>=div )
				carry = 1;
			else
			{
				temp_carry = q>>63;
				q <<= 1;
				q |= carry;
				carry = temp_carry;
				continue;
			}
		}

		rem -= div;
		rem -= (1 - carry);
		carry = 1;
		temp_carry = q >> 63;
		q <<= 1;
		q |= carry;
		carry = temp_carry;
	}

	if( (rem<<1)>=b )
		return q+1;
	return q;
}
*/

#else

/* This function returns the result of
 *    a1*(b1+1)
 *    ---------
 *      2^64
 *
 * It is supposed to be used internally in the code
 */
uint64 _file_multiply64(uint64 a1, uint64 b1)
{
	if(b1==MAXUINT64)
		return a1;
	b1++;

	uint64	a0 = a1 & 0xFFFFFFFF;
	uint64	b0 = b1 & 0xFFFFFFFF;
	a1 >>= 32;
	b1 >>= 32;

	uint64	ca = a1*b0;
	uint64	cb = a0*b1;

	return	(a1*b1 ) +
			(ca>>32) +
			(cb>>32) +
			((
				(ca&0xFFFFFFFF) +
				(cb&0xFFFFFFFF) +
				((a0*b0)>>32)
			)>>32);
}

/* This function returns the result of
 *    rem*2^64
 *    --------- - 1
 *       div
 */
uint64 _file_divide128(uint64 rem, uint64 div)
{
	uint64 q = 0;		// quotient
	uint64 carry = 0;
	uint64 temp_carry = 0;

	for(int i=0; i<64; i++)
	{
		temp_carry = rem >> 63;
		rem <<= 1;
		rem |= carry;
		carry = temp_carry;

		if( !carry )
		{
			if( rem>=div )
				carry = 1;
			else
			{
				temp_carry = q>>63;
				q <<= 1;
				q |= carry;
				carry = temp_carry;
				continue;
			}
		}

		rem -= div;
		rem -= (1 - carry);
		carry = 1;
		temp_carry = q >> 63;
		q <<= 1;
		q |= carry;
		carry = temp_carry;
	}

	if( (rem<<1)>=div )
		return q;
	return q-1;
}

#endif

uint64 dist_default(void *par, int64 k)
{
	uint64 *C = ((uint64 *) par) + 1;

	if(k>=(*((uint64 *) par)) || k<0)
		return MAXUINT64;
	return C[k];
}

void file::arithmeticStart()
{
	L = MAXUINT64;

	if( flagWrite )
	{
		B = 0;
		t = 0;
		if( d==NULL )
			d = (uint64 *) malloc( (AC_PAR_MEMORY_SIZE+1)*sizeof(uint64) );
		for(uint64 k=0; k<=AC_PAR_MEMORY_SIZE; k++)
			d[k] = 0x0;
	}
	else
		B = this->read(64);
}

void file::arithmeticStop()
{
	if( flagWrite )
	{
		uint64 k;
		uint64 tt;

		if(t)
		{
			k = t/64;
			tt = t%64;
			if( !tt )
			{
				tt = 64;
				k--;
			}
			this->write(d[k], tt);
			while( k-- )
				this->write(d[k], 64);
		}
		this->write(B, 64);
	}
}

void file::_arithmeticWrite(uint64 X, uint64 Y)
{
	uint64	A = B;
	uint64	k;
	uint8	c;

	B = B+X;
	L = Y-X;
	if( A>B )
	{
		k = 0;
		while( d[k]==MAXUINT64 && k<AC_PAR_MEMORY_SIZE )
		{
			d[k] = 0x0;
			k++;
		}
		d[k]++;

		if( k==AC_PAR_MEMORY_SIZE )
			error = 7;
	}

	c = 0;
	while( L<(0x8000000000000000>>c) && c<64 )
		c++;

	if(c)
	{
		uint8 p = 64-c;
		k = AC_PAR_MEMORY_SIZE;
		t += c;

		while( k )
		{
			d[k] = (d[k]<<c) + (d[k-1]>>p);
			k--;
		}
		d[0] = (d[0]<<c) + (B>>p);

		L <<= c;
		B <<= c;

		if( t>(AC_PAR_MEMORY_SIZE*64) )
		{
			p = t-(AC_PAR_MEMORY_SIZE*64);
			this->write(d[AC_PAR_MEMORY_SIZE], p);
			t -= p;
		}
	}
}

void file::arithmeticWrite(int64 s, void *par, uint64 (*cdf)(void *, int64), uint8 signedFlag)
{
	uint64	k  = 0;
#if PROB_BEHAVIOR_RETURN_CDF
	uint64	C0 = 0;
	uint64	C1 = cdf(par, k);
#else
	uint64	P0 = cdf(par, k);
#endif

	if( s )
	{
#if PROB_BEHAVIOR_RETURN_CDF
		this->_arithmeticWrite(_file_multiply(L,C1,C0), L);
#else
		this->_arithmeticWrite(_file_multiply64(L,P0), L);
#endif

		if( s<0 )
		{
			if( signedFlag )
				this->_arithmeticWrite(L/2, L);
			s = -s-1;
		}
		else
		{
			if( signedFlag )
				this->_arithmeticWrite(0, L/2);
			s =  s-1;
		}

		k  = 1;
#if PROB_BEHAVIOR_RETURN_CDF
		C0 = C1;
		C1 = cdf(par, k);
#else
		P0 = cdf(par, k);
#endif
	}

	while( s )
	{
		s--;

#if PROB_BEHAVIOR_RETURN_CDF
		this->_arithmeticWrite(_file_multiply(L,C1,C0), L);

		C0 = C1;
		C1 = cdf(par, ++k);
#else
		this->_arithmeticWrite(_file_multiply64(L,P0), L);

		P0 = cdf(par, ++k);
#endif
	}
#if PROB_BEHAVIOR_RETURN_CDF
	this->_arithmeticWrite(0, _file_multiply(L,C1,C0));
#else
	this->_arithmeticWrite(0, _file_multiply64(L,P0));
#endif
}

void file::_arithmeticRead(uint64 X, uint64 Y)
{
	B = B-X;
	L = Y-X;

	uint8 c = 0;
	while( L<(0x8000000000000000>>c) && c<64 )
		c++;
	if(c)
	{
		L = (L<<c);
		B = (B<<c) + this->read(c);
	}
}

int64 file::arithmeticRead(void *par, uint64 (*cdf)(void *, int64), uint8 signedFlag)
{
	uint8	signal;
	uint64	k  = 0;
#if PROB_BEHAVIOR_RETURN_CDF
	uint64	C0 = 0;
	uint64	C1 = cdf(par, k);
#else
	uint64	P0 = cdf(par, k);
#endif

	uint64	Z;

#if PROB_BEHAVIOR_RETURN_CDF
	Z = _file_multiply(L,C1,C0);
#else
	Z = _file_multiply64(L,P0);
#endif
	if( Z>B )
	{
		this->_arithmeticRead(0, Z);
		return 0;
	}
	this->_arithmeticRead(Z, L);

	if( signedFlag )
	{
		Z = L>>1;
		if( Z>B )
		{
			this->_arithmeticRead(0, Z);
			signal = 0;
		}
		else
		{
			this->_arithmeticRead(Z, L);
			signal = 1;
		}
	}
	else
		signal = 0;

	while(1)
	{
#if PROB_BEHAVIOR_RETURN_CDF
		C0 = C1;
		C1 = cdf(par, ++k);

		Z = _file_multiply(L,C1,C0);
#else
		P0 = cdf(par, ++k);

		Z = _file_multiply64(L,P0);
#endif
		if( Z>B )
		{
			this->_arithmeticRead(0, Z);
			if( signal )
				return -((int64) k);
			return k;
		}
		this->_arithmeticRead(Z, L);
	}
}
