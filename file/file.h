#ifndef FILE_H
#define FILE_H

#include "inteiros.h"
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include "math.h"
#include "probability.h"

#define FILE_FAKEWRITE			0

#define WRITE_MODE				1
#define READ_MODE				0

#define MAXUINT64				((uint64) 0xFFFFFFFFFFFFFFFF)
#define _file_BITMASK(k)		(MAXUINT64>>(64-k))

// Parameters of the Arithmetic Coder
#define AC_PAR_MEMORY_SIZE		8	// A positive integer smaller than 256
#define AC_PAR_INTERVAL_JUMP	16	// A power of 2 larger than 2

// Parameters of the Adaptive Run-Length / Golomb-Rice Encoding
#define RLGR_PAR_L				4
#define RLGR_PAR_U0				3
#define RLGR_PAR_D0				1
#define RLGR_PAR_U1				2
#define RLGR_PAR_D1				1

// Compatibility with machines with different endianness and/or smallest addressable memory unit
#ifdef __linux__
#  include <endian.h>
#  if CHAR_BIT==8
#    define htobe(x)	(x)
#    define betoh(x)	(x)
#    define CHAR_MASK	0xFF
#  elif CHAR_BIT==16
#    define htobe(x)	htobe16(x)
#    define betoh(x)	be16toh(x)
#    define CHAR_MASK	0xFFFF
#  elif CHAR_BIT==32
#    define htobe(x)	htobe32(x)
#    define betoh(x)	be32toh(x)
#    define CHAR_MASK	0xFFFFFFFF
#  elif CHAR_BIT==64
#    define htobe(x)	htobe64(x)
#    define betoh(x)	be64toh(x)
#    define CHAR_MASK	0xFFFFFFFFFFFFFFFF
#  endif
#else
#  include <iso646.h>
#  if CHAR_BIT!=8
#    error The smallest addressable memory should be 8 bits
#  endif
#  define htobe(x)		(x)
#  define betoh(x)		(x)
#  define CHAR_MASK		0xFF
#endif


uint64 dist_default(void *ptr, int64 k);

class file
{
public:
#if FILE_FAKEWRITE
	uint64	bytes;
#endif
	uint8	error;	// 0: no error
					// 1: unable to open desired file
                    // 2: unable to write to file. Occurs becouse the file was not correctly
                    //    opened to write
                    // 3: unable to read from file. Occurs when the file was not correctly
                    //    opened to read or has already reached end of file
					// 4: overflow. Occurs when trying to write/read more than 64 bits
                    // 5: unable to write/read. Could occurs when your machine is not LITLE
                    //    nor BIG ENDIAN, or the smallest addressable unit of the machine is
                    //    not 1, 2, 4 or 8 bytes
					// 6: unexpected error
					// 7: occurs in the arithmetic coder when propagating the carry. Means that
					//	  there is no settled bits in the buffer. To correct this, raises the
					//	  value of AC_PAR_MEMORY_SIZE and recompile the code

	/* Similar to "fopen".
	 * Open the file especified by 'filename' to read or write bit a bit (it is not possible
	 * to open a file to update)
	 *
	 * PARAMETERS:
	 *  filename: C string containig the name of the file to be opened
	 *  writeFlag: Number indicating a file access mode. It can be:
	 *	  READ_MODE (aka 0):  Open file for input operations. The file must exist.
	 *    WRITE_MODE (aka 1): Create an empty file for output operations. If a file with
	 *                        the same name exists, its contents are discarded and the file
	 *						  is treated as a new empty file.
	 */
	file(const char *filename, uint8 writeFlag);
	~file();

	/* Similar fo fwrite/fread
	 *
	 * Writes/reads an array of 'numel' elements, each one with 'bytes' bytes, from the block
	 * memory pointed by 'ptr'
	 *
	 * PARAMETERS:
	 *  ptr: Pointer to the array of elements to be written/read, converted to void*
     *  bytes: Size in bytes of each element to be written/read.
	 *  numel: Number of elements
	 *
	 * RETURNS
	 *  Error code
	 */
	uint8	write(void *ptr, size_t bytes, size_t numel);
	uint8	read(void *ptr, size_t bytes, size_t numel);

	/* Checks whether we are at the end-of-file, returning a value different of zero
	 * if it is and zero otherwise. This function always returns zero when in writting mode
	 */
	int		eof();

	/* Write 'bits' bits to file. Only the 'bits' less significant bits of 'data' are used.
	 * The others bits will be ignored. If 'bits' is not especified it is set to 1.
	 *
	 * PARAMETERS:
     *  data: data to be written
	 *  bits: number of bits in data
	 */
	void	write(uint8 data);
	void	write(uint64 data, uint8 bits);

	/* Read 'bits' bits from file. Return an integer with the bits given in the position of the
     * less significant bits. If 'bits' is not especified it is set to 1
	 *
	 * PARAMETERS
	 *  bits: number of bits in data
	 *
	 * RETURNS
	 *  Read data
	 */
	uint8	read();
	uint64	read(uint8 bits);

	/* Write/Read an integer to/from file using Golomb-Rice coding. 'bits' specify the number
	 * of bits to be employed to code the residual of 'data' when divided by 2^'bits'. The
	 * quotient is coded using an unary code.
	 *
	 * 'data' = p*(2^'bits')+q,     where p is coded using unary code and q with 'bits' bits
	 *
	 * PARAMETERS:
	 *  data: data to be written
	 *  bits: number of bits to code the residual
	 *
	 * RETURNS:
	 *  Read data
	 */
	void	grWrite(uint64 data, uint8 bits);
	uint64	grRead(uint8 bits);

	/* Write/Read a sequence of 'numel' signed integers given by 'ptr' to/from file using
	 * Adaptive Run-Length / Golomb-Rice Encoding.
	 *
	 * To use these functions is necessary to first call the function 'rlgrStart' before writing/reading
	 * any simbol. At the end is necessary to call the function 'rlgrStop' to write the remaining bits
	 * on the memory.
	 *
	 * PARAMETERS:
	 *  ptr: pointer to the array containing the data to be written/read
	 *  numel: number of elements to be written/read
     *  signedFlag: 0 when the values pointed by 'ptr' are non-negative, 1 otherwise. When
     *              signedFlag!=0 the funtion will remap posive numbers into even positive
     *              numbers, and negative numbers into odd positive numbers before writing and
     *              the inverse remapping when reading.
	 *
	 * RETURNS:
     *  Error code
	 */
	uint8	rlgrWrite(int64 *ptr, uint64 numel, uint8 signedFlag = 1);
	uint8	rlgrRead(int64 *ptr, uint64 numel, uint8 signedFlag = 1);
	void	rlgrStart();
	void	rlgrStop();

	/* Write/Read the simbols 's' using arithmetic coding. The cumulative count of the occurence of each
	 * simbol is provided by the user by the function 'cdf'.
	 * To use these functions is necessary to first call the function 'arithmeticStart' before writing/reading
	 * any simbol. At the end is necessary to call the function 'arithmeticStop' to write the remaining bits
	 * on the memory.
	 *
	 * PARAMETERS:
	 *	s: simbol to be encoded
	 *  par: a pointer to the parameters to be passed to the 'prob' and 'simbol' function as first element
	 *	cdf(void *par, uint64 s): a function pointer, provided by the user, that returns the cumulative
	 *							count of occurence of the simbol 's', minus 1. The returned value should be
	 *							adjusted so that the last simbol will an associated value of excatly 2^64-1
	 *							This function should satisfy:
	 *								prob(par, s-1) <= prob(par, s) for all s>0
	 *							If prob(par, s-1) = prob(par, s), then simbol 's' should never appear during coding
	 *
	 * RETURNS:
	 *	The decoded simbol
	 */
	void	arithmeticWrite(int64 s, void *par, uint64 (*cdf)(void *, int64) = dist_default, uint8 signedFlag = 1);
	int64	arithmeticRead(void *par, uint64 (*cdf)(void *, int64) = dist_default, uint8 signedFlag = 1);
	void	arithmeticStart();
	void	arithmeticStop();
private:
	FILE	*fid;			// File identifier
	uint64	buf0, buf1;		// Bit's buffer
	uint8	bufel;			// "buffer elements" = number of bits in the buffer
	uint8	byte;
	uint8	flagWrite;		// 1: when file is writable, 0: otherwise

	void	fill();			// Write bits stored in 'buf0' and 'buf1' to file
	void	flush();		// Read bits from file and store in 'buf0' and 'buf1'

	uint64	B, L, t;
	uint64	*d;

	void	_arithmeticWrite(uint64 X, uint64 Y);
	void	_arithmeticRead(uint64 X, uint64 Y);
};

#endif // FILE_H
