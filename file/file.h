#ifndef FILE_H
#define FILE_H

#include <stdio.h>
#include <stdint.h>

#define	UINTMAX_BITS	(sizeof(uintmax_t)*8)
#define LEAST8_BITS		(sizeof(uint_least8_t)*8)
#define MASK(k)			((((uintmax_t) 0x1)<<(k))-1)

#ifndef UINT_LEAST8_MAX
# define UINT_LEAST8_MAX	MASK(LEAST8_BITS)
#endif
#ifndef UINTMAX_MAX
# define UINTMAX_MAX		((uintmax_t) -1)
#endif

// Parameters of the RLGR coder/decoder
#define	L	4
#define	U0	3
#define	D0	1
#define	U1	2
#define	D1	1

class file
{
private:
	uintmax_t		data;
	uint_least8_t	bits;
	uint_least8_t	flagWrite;
	FILE			*fid;

	void	flush();
	void	fill();

public:
	file(char *filename, uint_least8_t flagWrite);
	~file();

	uint_least8_t	openError() { return this->fid==NULL; }

	/* Checks whether the end-of-File has been reached.
	 *
	 * Returns:
	 * 1, when at end-of-file
	 * 0, otherwise
	 */
	uint_least8_t	eof();

	/* Reads "bits" bits from stream. If the number of bits is not
	 * especified it is assumed to be 1
	 */
	uint_least8_t	read();
	uintmax_t		read(uint_least8_t bits);

	/* Writes "data" to stream using "bits" bits. If the number of bits is
	 * not especified it is assumed to be 1
	 */
	void	write(uint_least8_t data);
	void	write(uintmax_t data, uint_least8_t bits);

	/* Functions similar to fread and fwrite
	 *
	 * Writes/Reads an array of count elements, each one with a size of
	 * "size" bytes, from the block of memory pointed by "ptr" to the
	 * current position in the stream.
	 */
	void	read(void *ptr, size_t size, size_t count);
	void	write(void *ptr, size_t size, size_t count);

	/* Reads/Writes an integer using Goulamb-Rice coding with "bits" bits
	 */
	uintmax_t	grRead(uint_least8_t bits);
	void		grWrite(uintmax_t data, uint_least8_t bits);

	/* Reads/Writes a sequence of integers using Run Length Goulamb-Rice
	 * coding
	 *
	 * INPUTS:
	 *	seq:		the sequence to be written/readden
	 *	N:			number of elements on the sequence
	 *	flagSigned:	1, when the sequence may contan positive and negative numbers
	 *				0, when the sequence cotains only non-negative numbers
	 */
	void	rlgrRead(intmax_t *seq, size_t N, uint_least8_t flagSigned=1);
	void	rlgrWrite(intmax_t *seq, size_t N, uint_least8_t flagSigned=1);
};

uintmax_t _s2u(intmax_t val);
intmax_t _u2s(uintmax_t val);

#endif // FILE_H

