#include "file.h"

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
file::file(const char *filename, uint8 writeFlag)
{
	flagWrite = writeFlag;
	bufel  = 0;
	d = NULL;

#if FILE_FAKEWRITE
	bytes = 0;
	fid = NULL;
	error = 0;
#else
	if(writeFlag)
		fid = fopen(filename, "w");
	else
		fid = fopen(filename, "r");

	if( fid==NULL )
		error = 1;
	else
		error = 0;
#endif
}

file::~file()
{
#if FILE_FAKEWRITE
	if( flagWrite )
	{
		this->flush();

		if( bufel )
			bytes++;
	}
#else
	if( fid!=NULL )
	{
		if( flagWrite )
		{
			this->flush();

			if( bufel )
			{
				byte = htobe( (buf0<<(CHAR_BIT-bufel)) & CHAR_MASK );
				fwrite(&byte, CHAR_BIT/8, 1, fid);
			}
		}

		fclose(fid);
	}
#endif

	if( d!=NULL )
		free(d);
}

/* Checks whether we are at the end-of-File, returning a value different from zero
 * if it is and zero otherwise. This function always returns zero when in writting mode
 */
int file::eof()
{
	if(flagWrite)
		return 0;
	this->fill();
#if FILE_FAKEWRITE
	return 0;
#else
	return feof(fid);
#endif
}

/* Fill 'buf1' and 'buf0' with bits read from the file
 */
void file::fill()
{
	while( bufel<=(128-CHAR_BIT) )
	{
#if FILE_FAKEWRITE
		byte = 0;
#else
		if( !fread(&byte, CHAR_BIT/8, 1, fid) )
			break;
#endif

		bufel += CHAR_BIT;
		if( bufel>64 )
			buf1 = (buf1<<CHAR_BIT) + (buf0>>(64-CHAR_BIT));
		buf0 = (buf0<<CHAR_BIT) + betoh(byte);
	}
}

/* Flush the bits in 'buf1' and 'buf0' to the file
 */
void file::flush()
{
	while( bufel>=CHAR_BIT )
	{
		if( bufel>64 )
		{
			bufel -= CHAR_BIT;
#if !FILE_FAKEWRITE
			if( bufel>=64 )
				byte = htobe( (buf1>>(bufel-64)) & CHAR_MASK );
			else
				byte = htobe( ((buf1<<(64-bufel)) & CHAR_MASK) + (buf0>>bufel) );
#endif
		}
		else
		{
			bufel -= CHAR_BIT;
#if !FILE_FAKEWRITE
			byte = htobe( (buf0>>bufel) & CHAR_MASK );
#endif
		}

#if FILE_FAKEWRITE
		bytes++;
#else
		if( !fwrite(&byte, CHAR_BIT/8, 1, fid) )
			break;
#endif
	}
}
