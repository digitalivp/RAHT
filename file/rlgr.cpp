#include "file.h"

/* Write an integer to file using Golomb-Rice coding. 'bits' specify the number
 * of bits to be employed to code the residual of 'data' when divided by 2^'bits'. The
 * quotient is coded using an unary code.
 *
 * 'data' = p*(2^'bits')+q,     where p is coded using unary code and q with 'bits' bits
 *
 * PARAMETERS:
 *  data: data to be written
 *  bits: number of bits to code the residual
 */
void file::grWrite(uint64 data, uint8 bits)
{
    if( error )
		return;

	uint64 p = data>>bits;

	if( p<64 )
	{
		// Code p with 1s run
		this->write(MAXUINT64-1, p+1);

		// Code residual as integer
		this->write(data, bits);
	}
	else
	{
		// Code 64 1s, and 64 bit unsigned int
		this->write(MAXUINT64, 64);
		this->write(data, 64);
	}
}

/* Read an integer from file using Golomb-Rice coding. 'bits' specify the number
 * of bits to be employed to code the residual of 'data' when divided by 2^'bits'. The
 * quotient is coded using an unary code.
 *
 * 'data' = p*(2^'bits')+q,     where p is coded using unary code and q with 'bits' bits
 *
 * PARAMETERS:
 *  bits: number of bits to code the residual
 *
 * RETURNS:
 *  Read data
 */
uint64 file::grRead(uint8 bits)
{
    if( error )
		return 0;

	uint64 p = 0;

	while( this->read() )
	{
		p++;
		if( p>=64 )
			return this->read(64);
	}

	return (p<<bits) + this->read(bits);
}

/* Convert a signed number to an unsigned represention
 */
uint64 dir(int64 val)
{
	uint64 out;

	if( val<0 )
	{
		out = -val;
		return ((out<<1)-1);
	}
	out = val;
	return (out<<1);
}

/* Convert an unsigned representation of a number to its signed value
 */
int64 inv(uint64 val)
{
	int64 out = val>>1;

	if( val&0x1 )
		return (-out-1);
	return (out);
}

void file::rlgrStart()
{
	this->L = 0;
	this->B = 2*RLGR_PAR_L;
	this->t = 0;

	/*
	k_P = 0;				// k_P virou this->L
	k_RP = 2*RLGR_PAR_L;	// k_RP virou this->B
	m = 0;					// m virou this->t
	*/
}

void file::rlgrStop()
{
	if( this->flagWrite )
		if(this->t)
			this->write(0);
}

/* Write a sequence of 'numel' signed integers given by 'ptr' to file using
 * Adaptive Run-Length / Golomb-Rice Encoding.
 *
 * PARAMETERS:
 *  ptr: pointer to the array containing the data to be written
 *  numel: number of elements to be written/read
 *
 * RETURNS:
 *  Error code
 */
uint8 file::rlgrWrite(int64 *ptr, uint64 numel, uint8 signedFlag)
{
	if( error )
		return error;

	/* K_RP virou this->B
	 * K_P virou this->L
	 * m virou this->t
	 */

	uint64	k_R;
	uint64	k;
	uint64	p;

	uint64	n;
	uint64	u;

	for(n=0; n<numel; n++)
	{
		if( this->B>(32*RLGR_PAR_L) )
			this->B = 32*RLGR_PAR_L;
		if( this->L>(16*RLGR_PAR_L) )
			this->L = 16*RLGR_PAR_L;

		k = this->L/RLGR_PAR_L;
		k_R = this->B/RLGR_PAR_L;

		u = signedFlag?dir(ptr[n]):ptr[n];

		if( k==0 )
		{
			// "No run" mode
			this->grWrite(u, k_R);
			if( error )
				return error;

			// Adapt this->B
			p = u>>k_R;
			if( !p )
			{
				if( this->B<2 )
					this->B = 0;
				else
					this->B -= 2;
			}
			else if( p>1 )
			{
				k_R = this->B+p-1;
				if( k_R>=this->B )
				{
					this->B = k_R;
					if( this->B>(RLGR_PAR_L*32) )
						this->B = RLGR_PAR_L*32;
				}
				else
					this->B = RLGR_PAR_L*32;
			}

			// Adapt this->L
			if( !u )
			{
				k = this->L+RLGR_PAR_U0;
				if( k>=this->L )
					this->L = k;
				else
					this->L = MAXUINT64;
			}
			else
			{
				if( this->L<RLGR_PAR_D0 )
					this->L = 0;
				else
					this->L -= RLGR_PAR_D0;
			}

			this->t = 0;
		}
		else
		{
			// "Run" mode
			if( !u )
			{
				// Continue run of 0s
				this->t++;
				if( this->t==(0x1<<k) )
				{
					this->write(0);
					if( error )
						return error;

					// Adapt this->L
					k = this->L + RLGR_PAR_U1;
					if( k>=this->L )
						this->L = k;
					else
						this->L = MAXUINT64;

					this->t = 0;
				}
			}
			else
			{
				// End run of 0s
				this->write(1);
				this->write(this->t,k);
				this->grWrite(u-1,k_R);
				if( error )
					return error;

				// Adapt this->B
				p = (u-1)>>k_R;
				if( !p )
				{
					if( this->B<2 )
						this->B = 0;
					else
						this->B -= 2;
				}
				else if( p>1 )
				{
					k_R = this->B+p-1;
					if( k_R>=this->B )
					{
						this->B = k_R;
						if( this->B>(RLGR_PAR_L*32) )
							this->B = RLGR_PAR_L*32;
					}
					else
						this->B = RLGR_PAR_L*32;
				}

				// Adapt this->L
				if( this->L>RLGR_PAR_D1 )
					this->L -= RLGR_PAR_D1;
				else
					this->L = 0;

				this->t = 0;
			}
		}
	}

	return 0;
}

/* Read a sequence of 'numel' signed integers given by 'ptr' from file using
 * Adaptive Run-Length / Golomb-Rice Encoding.
 *
 * PARAMETERS:
 *  ptr: pointer to the array containing the data to be read
 *  numel: number of elements to be written/read
 *
 * RETURNS:
 *  Error code
 */
uint8 file::rlgrRead(int64 *ptr, uint64 numel, uint8 signedFlag)
{
	if( error )
		return error;

	/* K_RP virou this->B
	 * K_P virou this->L
	 * m virou this->t
	 */

	uint64	k_R;
	uint64	k;
	uint64	p;

	uint64	n = 0;
	uint64	u;

	while( n<numel )
	{
		if( this->B>(32*RLGR_PAR_L) )
			this->B = 32*RLGR_PAR_L;
		if( this->L>(16*RLGR_PAR_L) )
			this->L = 16*RLGR_PAR_L;

		k = this->L/RLGR_PAR_L;
		k_R = this->B/RLGR_PAR_L;

		if( k==0 )
		{
			// "No run" mode
			u = this->grRead(k_R);
			if( error )
				return error;

			// Write decoded output
			ptr[n++] = signedFlag?inv(u):u;

			// Adapt this->B
			p = u>>k_R;
			if( !p )
			{
				if( this->B<2 )
					this->B = 0;
				else
					this->B -= 2;
			}
			else if( p>1 )
			{
				k_R = this->B+p-1;
				if( k_R>=this->B )
				{
					this->B = k_R;
					if( this->B>(RLGR_PAR_L*32) )
						this->B = RLGR_PAR_L*32;
				}
				else
					this->B = RLGR_PAR_L*32;
			}

			// Adapt this->L
			if( !u )
			{
				k = this->L+RLGR_PAR_U0;
				if( k>=this->L )
					this->L = k;
				else
					this->L = MAXUINT64;
			}
			else
			{
				if( this->L<RLGR_PAR_D0 )
					this->L = 0;
				else
					this->L -= RLGR_PAR_D0;
			}
		}
		else
		{
			// "Run" mode
			if( this->read() )
			{
				this->t = this->read(k);
				u = this->grRead(k_R)+1;
				if( error )
					return error;

				// Write decoded output
				while(this->t--)
					ptr[n++] = 0;
				ptr[n++] = signedFlag?inv(u):u;

				// Adapt this->B
				p = (u-1)>>k_R;
				if( !p )
				{
					if( this->B<2 )
						this->B = 0;
					else
						this->B -= 2;
				}
				else if( p>1 )
				{
					k_R = this->B+p-1;
					if( k_R>=this->B )
					{
						this->B = k_R;
						if( this->B>(RLGR_PAR_L*32) )
							this->B = RLGR_PAR_L*32;
					}
					else
						this->B = RLGR_PAR_L*32;
				}

				// Adapt this->L
				if( this->L>RLGR_PAR_D1 )
					this->L -= RLGR_PAR_D1;
				else
					this->L = 0;
			}
			else
			{
				if( error )
					return error;

				this->t = 0x1<<k;

				// Write 0s to output
				while(this->t-- && n<numel)
					ptr[n++] = 0;

				// Adapt this->L
				k = this->L + RLGR_PAR_U1;
				if( k>=this->L )
					this->L = k;
				else
					this->L = MAXUINT64;
			}
		}
	}

	return 0;
}
