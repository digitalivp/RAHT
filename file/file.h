#ifndef FILE_H
#define FILE_H

#include <stdio.h>
#include <stdint.h>

#define LEAST8_BITS         (sizeof(uint_least8_t)*8)
#define MASK(k)             ((((uint64_t) 0x1)<<(k))-1)

#ifndef UINT_LEAST8_MAX
# define UINT_LEAST8_MAX    MASK(LEAST8_BITS)
#endif
#ifndef UINT64_MAX
# define UINT64_MAX        ((uint64_t) -1)
#endif

// Parameters of the RLGR coder/decoder
#define	L   4
#define	U0  3
#define	D0  1
#define	U1  2
#define	D1  1

class file
{
private:
    uint64_t        data;
    uint_least8_t   bits;
    uint_least8_t   flagWrite;
    FILE            *fid;
    uint64_t        filesize;
    void            flush();
    void            fill();

public:

    file(char *filename, uint_least8_t flagWrite);
    ~file();

    uint64_t file_size() {return this->filesize;}

    uint_least8_t   openError()
    {
        return this->fid==NULL;
    }

    /* Checks whether the end-of-File has been reached.
     *
     * Returns:
     * 1, when at end-of-file
     * 0, otherwise
     */
    uint_least8_t   eof();

    /* Reads "bits" bits from stream. If the number of bits is not
     * especified it is assumed to be 1
     */
    uint_least8_t   read();
    uint64_t        read(uint_least8_t bits);

    /* Writes "data" to stream using "bits" bits. If the number of bits is
     * not especified it is assumed to be 1
     */
    void            write(uint_least8_t data);
    void            write(uint64_t data, uint_least8_t bits);

    /* Functions similar to fread and fwrite
     *
     * Writes/Reads an array of count elements, each one with a size of
     * "size" bytes, from the block of memory pointed by "ptr" to the
     * current position in the stream.
     */
    void            read(void *ptr, size_t size, size_t count);
    void            write(void *ptr, size_t size, size_t count);

    /* Reads/Writes an integer using Goulamb-Rice coding with "bits" bits
     */
    uint64_t        grRead(uint_least8_t bits);
    void            grWrite(uint64_t data, uint_least8_t bits);

    /* Reads/Writes a sequence of integers using Run Length Goulamb-Rice
     * coding
     *
     * INPUTS:
     *	seq:		the sequence to be written/readden
     *	N:			number of elements on the sequence
     *	flagSigned:	1, when the sequence may contan positive and negative numbers
     *				0, when the sequence cotains only non-negative numbers
     */
    void            rlgrRead(int64_t *seq, size_t N, uint_least8_t flagSigned=1);
    void            rlgrWrite(int64_t *seq, size_t N, uint_least8_t flagSigned=1);
};

uint64_t _s2u(int64_t val);
int64_t _u2s(uint64_t val);

#endif // FILE_H
