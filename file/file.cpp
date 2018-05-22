#include "file.h"

uintmax_t _s2u(intmax_t val)
{
	uintmax_t dval;

	if( val<0 )
	{
		dval = -val;
		return (dval<<1)-1;
	}
	else
	{
		dval = val;
		return (dval<<1);
	}
}

intmax_t _u2s(uintmax_t val)
{
	intmax_t dval = val>>1;
	if( val&0x1 )
		return -dval-1;
	else
		return dval;
}

file::file(char *filename, uint_least8_t flagWrite)
{
	this->flagWrite = flagWrite;

	if( flagWrite )
		this->fid = fopen(filename, "wb");
	else
		this->fid = fopen(filename, "rb");

	this->bits = 0;
}

file::~file()
{
	if( this->fid==NULL )
		return;

	if( this->flagWrite )
	{
		uint_least8_t r = this->bits%LEAST8_BITS;
		if( r )
			this->write(0, LEAST8_BITS-r);
		else
			this->flush();
	}

	fclose(this->fid);
}

uint_least8_t file::eof()
{
	this->fill();
	return	(feof(this->fid)) &&
			(!(this->bits/LEAST8_BITS)) &&
			(!(this->data&MASK(this->bits)));
}

void file::flush()
{
	uint_least8_t	data;
	while( this->bits>=LEAST8_BITS )
	{
		this->bits -= LEAST8_BITS;
		data = ((this->data)>>(this->bits)) & UINT_LEAST8_MAX;
		fwrite(&data, sizeof(uint_least8_t), 1, this->fid);
	}
}

void file::fill()
{
	uint_least8_t	data;
	while( this->bits<=(UINTMAX_BITS-LEAST8_BITS) )
	{
		if( fread(&data, sizeof(uint_least8_t), 1, this->fid) )
		{
			this->data = (this->data<<LEAST8_BITS) + data;
			this->bits += LEAST8_BITS;
		}
		else
			return;
	}
}

uint_least8_t file::read()
{
	if( !this->bits )
		this->fill();
	this->bits--;
	return (this->data>>this->bits)&0x1;
}

uintmax_t file::read(uint_least8_t bits)
{
	if( bits>(UINTMAX_BITS-LEAST8_BITS) )
	{
		uintmax_t	data = this->read(bits-UINTMAX_BITS/2)<<(UINTMAX_BITS/2);
		return	data + this->read(UINTMAX_BITS/2);
	}

	this->fill();
	this->bits -= bits;
	return	(this->data>>this->bits) & MASK(bits);
}

void file::read(void *ptr, size_t size, size_t count)
{
	if( size==8 )
	{
		while( count-- )
		{
			*((uint64_t *) ptr) = this->read( 64 );
			ptr = ((uint64_t *) ptr) + 1;
		}
		return;
	}
	if( size==4 )
	{
		while( count-- )
		{
			*((uint32_t *) ptr) = this->read( 32 );
			ptr = ((uint32_t *) ptr) + 1;
		}
		return;
	}
	if( size==2 )
	{
		while( count-- )
		{
			*((uint16_t *) ptr) = this->read( 16 );
			ptr = ((uint16_t *) ptr) + 1;
		}
		return;
	}
	count *= size;
	while( count-- )
	{
		*((uint8_t *) ptr) = this->read( 8 );
		ptr = ((uint8_t *) ptr) + 1;
	}
}

void file::write(uint_least8_t data)
{
	this->data <<= 1;
	if( data )
		this->data++;
	this->bits++;

	if( this->bits>=LEAST8_BITS )
		this->flush();
}

void file::write(uintmax_t data, uint_least8_t bits)
{
	if( bits>(UINTMAX_BITS-LEAST8_BITS) )
	{
		this->write(data>>(UINTMAX_BITS/2), bits-UINTMAX_BITS/2);
		this->write(data&MASK(UINTMAX_BITS/2), UINTMAX_BITS/2);
		return;
	}

	this->data = (this->data<<bits) + data;
	this->bits += bits;
	this->flush();
}

void file::write(void *ptr, size_t size, size_t count)
{
	if( size==8 )
	{
		while( count-- )
		{
			this->write( *((uint64_t *) ptr), 64 );
			ptr = ((uint64_t *) ptr) + 1;
		}
		return;
	}
	if( size==4 )
	{
		while( count-- )
		{
			this->write( *((uint32_t *) ptr), 32 );
			ptr = ((uint32_t *) ptr) + 1;
		}
		return;
	}
	if( size==2 )
	{
		while( count-- )
		{
			this->write( *((uint16_t *) ptr), 16 );
			ptr = ((uint16_t *) ptr) + 1;
		}
		return;
	}
	count *= size;
	while( count-- )
	{
		this->write( *((uint8_t *) ptr), 8 );
		ptr = ((uint8_t *) ptr) + 1;
	}
}

uintmax_t file::grRead(uint_least8_t bits)
{
	uintmax_t	p = 0;

	while( this->read() )
	{
		p++;
		if( p>=32 )
			return this->read(32);
	}
	
	return (p<<bits) + this->read(bits);
}

void file::grWrite(uintmax_t data, uint_least8_t bits)
{
	uintmax_t	p = data>>bits;

	if( p<32 )
	{
		this->write(MASK(p+1)-1, p+1);
		this->write(data&MASK(bits), bits);
	}
	else
	{
		this->write(MASK(32), 32);
		this->write(data, 32);
	}
}

void file::rlgrRead(intmax_t *seq, size_t N, uint_least8_t flagSigned)
{
	uintmax_t	u;

	uintmax_t	k_P = 0;
	uintmax_t	k_RP = 2*L;
	uintmax_t	m = 0;

	uintmax_t	k;
	uintmax_t	k_R;
	uintmax_t	p;

	size_t n=0;
	while( n<N )
	{
		k = k_P/L;
		k_R = k_RP/L;

		if( k )
		{
			// "Run" mode
			m = 0;
			while( this->read() )
			{
				m += 0x1<<k;

				// Adapt k_P
				k_P += U1;
				k = k_P/L;
			}
			m += this->read(k);
			while(m--)
				seq[n++] = 0;
			if( n>=N )
				break;
			u = this->grRead(k_R);
			seq[n++] = flagSigned?_u2s(u+1):u+1;

			// Adapt k_RP
			p = u>>k_R;
			if( p )
			{
				k_RP += p-1;
				if( k_RP>32*L )
					k_RP = 32*L;
			}
			else
			{
				if( k_RP<2 )
					k_RP = 0;
				else
					k_RP -= 2;
			}

			// Adapt k_P
			if( k_P<D1 )
				k_P = 0;
			else
				k_P -= D1;
		}
		else
		{
			// "No run" mode
			u = this->grRead(k_R);
			seq[n++] = flagSigned?_u2s(u):u;

			// Adapt k_RP
			p = u>>k_R;
			if( p )
			{
				k_RP = k_RP+p-1;
				if( k_RP>32*L )
					k_RP = 32*L;
			}
			else
			{
				if( k_RP<2 )
					k_RP = 0;
				else
					k_RP -= 2;
			}

			// Adapt k_P
			if( u )
			{
				if( k_P<D0 )
					k_P = 0;
				else
					k_P -= D0;
			}
			else
				k_P += U0;
		}
	}
}

void file::rlgrWrite(intmax_t *seq, size_t N, uint_least8_t flagSigned)
{
	uintmax_t	u;

	uintmax_t	k_P = 0;
	uintmax_t	k_RP = 2*L;
	uintmax_t	m = 0;

	uintmax_t	k;
	uintmax_t	k_R;
	uintmax_t	p;

	for(size_t n=0; n<N; n++)
	{
		u = flagSigned?_s2u(seq[n]):seq[n];

		k = k_P/L;
		k_R = k_RP/L;

		if( k )
		{
			// "Run" mode

			if( u )
			{
				u--;

				// End run of 0s
				this->write(0);
				this->write(m, k);
				this->grWrite(u, k_R);

				// Adapt k_RP
				p = u>>k_R;
				if( p )
				{
					k_RP += p-1;
					if( k_RP>(32*L) )
						k_RP = 32*L;
				}
				else
				{
					if( k_RP<2 )
						k_RP = 0;
					else
						k_RP -= 2;
				}

				// Adapt k_P
				if( k_P<D1 )
					k_P = 0;
				else
					k_P -= D1;

				m = 0;
			}
			else
			{
				// Contine run of 0s
				m++;
				if( m==(0x1<<k) )
				{
					this->write(1);

					// Adapt k_P
					k_P += U1;

					m = 0;
				}
			}
		}
		else
		{
			// "No run" mode
			this->grWrite(u, k_R);

			// Adapt k_RP
			p = u>>k_R;
			if( p )
			{
				k_RP = k_RP+p-1;
				if( k_RP>32*L )
					k_RP = 32*L;
			}
			else
			{
				if( k_RP<2 )
					k_RP = 0;
				else
					k_RP -= 2;
			}

			// Adapt k_P
			if( u )
			{
				if( k_P<D0 )
					k_P = 0;
				else
					k_P -= D0;
			}
			else
				k_P += U0;

			m = 0;
		}
	}

	if( k && !u )
	{
		this->write(0);
		this->write(m, k_P/L);
	}
}

