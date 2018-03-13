#include "file.h"

/* Similar fo fwrite
 *
 * Writes an array of 'numel' elements, each one with 'bytes' bytes, from the block
 * memory pointed by 'ptr'
 *
 * PARAMETERS:
 *  ptr: Pointer to the array of elements to be written, converted to void*
 *  bytes: Size in bytes of each element to be written.
 *  numel: Number of elements
 *
 * RETURNS
 *  Error code
 */
uint8 file::write(void *ptr, size_t bytes, size_t numel)
{
	if( error )
		return error;

#if BYTE_ORDER==LITTLE_ENDIAN
	uint8   *data;
	size_t  count;

	data = (uint8 *) ptr;
	while( numel-- )
	{
		count = bytes/sizeof(uint8);
		while( count-- )
			this->write(data[count], sizeof(uint8)*8);
		data += bytes;
	}

	return error;
#elif BYTE_ORDER==BIG_ENDIAN
	uint8   *data;
	size_t  count;

	data = (uint8 *) ptr;
	while( numel-- )
	{
		count = bytes/sizeof(uint8);
		while( count-- )
			this->write(*(data++), sizeof(uint8)*8);
	}

	return error;
#else
	error = 5;
	return error;
#endif
}

/* Similar fo fread
 *
 * Reads an array of 'numel' elements, each one with 'bytes' bytes, from the block
 * memory pointed by 'ptr'
 *
 * PARAMETERS:
 *  ptr: Pointer to the array of elements to be read, converted to void*
 *  bytes: Size in bytes of each element to be read.
 *  numel: Number of elements
 *
 * RETURNS
 *  Error code
 */
uint8 file::read(void *ptr, size_t bytes, size_t numel)
{
	if( error )
		return error;

#if BYTE_ORDER==LITTLE_ENDIAN
	uint8   *data;
	size_t  count;

	data = (uint8 *) ptr;
	while( numel-- )
	{
		count = bytes/sizeof(uint8);
		while( count-- )
			data[count] = this->read( sizeof(uint8)*8 );
		data += bytes;
	}

	return error;
#elif BYTE_ORDER==BIG_ENDIAN
	uint8   *data;
	size_t  count;

	data = (uint8 *) ptr;
	while( numel-- )
	{
		count = bytes/sizeof(uint8);
		while( count-- )
			*(data++) = this->read( sizeof(uint8)*8 );
	}

	return error;
#else
	error = 5;
	return error;
#endif
}

/* Write 1 bit to file. Only the less significant bit of 'data' is used.
 * The others bits will be ignored.
 *
 * PARAMETERS:
 *  data: data to by written
 */
void file::write(uint8 data)
{
	if( error )
		return;

	if( bufel>=128 )
	{
		this->flush();
		if( bufel>=128 )
		{
			error = 2;
			return;
		}
	}

	bufel++;
	if( bufel>64 )
		buf1 = (buf1<<1) + (buf0>>63);
	buf0 = (buf0<<1) + (data&0x1);
}

/* Write 'bits' bits to file. Only the 'bits' less significant bits of 'data' are used.
 * The others bits will be ignored. If 'bits' is not especified it is set to 1.
 *
 * PARAMETERS:
 *  data: data to by written
 *  bits: number of bits in data
 */
void file::write(uint64 data, uint8 bits)
{
	if( error )
		return;

	if( bits>64 )
	{
		error = 4;
		return;
	}
	if( !bits )
		return;

	if( (bits+bufel)>128 )
	{
		this->flush();
		if( (bits+bufel)>128 )
		{
			error = 2;
			return;
		}
	}

	bufel += bits;
	if(bits==64)
	{
		buf1 = buf0;
		buf0 = data;
	}
	else if(bits)
	{
		if( bufel>64 )
		{
			buf1 = (buf1<<bits) + (buf0>>(64-bits));
			buf0 = (buf0<<bits) + (data & _file_BITMASK(bits));
		}
		else
			buf0 = (buf0<<bits) + (data & _file_BITMASK(bits));
	}
}

/* Read 1 bit from file. Return an integer with the bit given in the position of the
 * less significant bit.
 *
 * RETURNS
 *  Read data
 */
uint8 file::read()
{
	if( error )
		return 0;

	if( !bufel )
	{
		this->fill();
		if( !bufel )
		{
			error = 3;
			return 0;
		}
	}

	if( bufel>64 )
	{
		bufel--;
		return (buf1>>(bufel-64)) & 0x1;
	}
	bufel--;
	return (buf0>>bufel) & 0x1;
}

/* Read 'bits' bits from file. Return an integer with the bits given in the position of the
 * less significant bits. If 'count' is not especified it is set to 1
 *
 * PARAMETERS
 *  bits: number of bits in data
 *
 * RETURNS
 *  Read data
 */
uint64 file::read(uint8 bits)
{
	if( error )
		return 0;

	if( bits>64 )
	{
		error = 4;
		return 0;
	}
	else if( !bits )
		return 0;

	if( bufel<bits )
	{
		this->fill();
		if( bufel<bits )
		{
			error = 3;
			return 0;
		}
	}

	if( bufel>64 )
	{
		bufel -= bits;
		if( bufel>=64 )
			return (buf1>>(bufel-64)) & _file_BITMASK(bits);
		return ((buf1<<(64-bufel))&_file_BITMASK(bits)) + (buf0>>bufel);
	}
	bufel -= bits;
	return (buf0>>bufel) & _file_BITMASK(bits);
}
