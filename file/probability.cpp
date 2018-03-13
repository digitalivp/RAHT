#include "probability.h"

/* LAPLACE DISTRIBUTION
 *
 * Parameter type:	double *
 * Parameter value:	{sd}
 *
 *			  1      /   |x| \            __
 * PDF(x) =	----- exp| - --- |, where b*\/2  = sd
 *			 2*b     \    b  /
 */
uint64 dist_laplace(void *par, int64)
{
	//uint64	C = round( exp(-0.7071067811865475244/(*((double *) par)))*18446744073609551616.0 );
	uint64	C = round( exp(-1.4142135623730950488/(*((double *) par)))*18446744073609551616.0 );
	return	MAXUINT64 - C;
}

/* NORMAL DSITRIBUTION
 *
 * Parameter type:	double *
 * Parameter value	{sd}
 *
 *				1        /    x^2   \
 * PDF(x) =	-----____ exp| - ------ |
 *			sd \/2 pi    \   2 sd^2 /
 */
uint64 dist_normal(void *par, int64 s)
{
	double	k = *((double *) par)*2.8284271247461900976;
	uint64	C;

	if( s<0 )
		s = -s;

	if( s )
		C = round( erfc(((double) 2*s+1)/k)/erfc(((double) 2*s-1)/k)*18446744073609551616.0 );
	else
		C = round( erfc(1.0/k)*18446744073609551616.0 );

	return	MAXUINT64 - C;
}

/* CAUCHY DISTRIBUTION
 *
 * Parameter type:	double *
 * Parameter value	{gamma}
 *
 *				1     /    gamma^2    \
 * PDF(x) =	--------- | ------------- |
 *			pi*gamma  \ x^2 + gamma^2 /
 */
uint64 dist_cauchy(void *par, int64 s)
{
	double	k = (*((double *) par))*2;
	uint64	C;

	if( s<0 )
		s = -s;

	if( s )
		C = round( atan(k/(2*s+1))/atan(k/(2*s-1))*18446744073609551616.0 );
	else
		C = round( atan(k)*11743562013064343552.0 );

	return	MAXUINT64 - C;
}

/* UNIFORM DISTRIBUTION
 *
 * Parameter type:	double *
 * Parameter value	{max}
 *
 *			/	  1
 *			| ---------	, if |x| <= max
 * PDF(x) =	{ 2 max + 1
 *			|
 *			\	<< 1	, otherwise (the probability for |x|>max is not exactly 0 to allow this
 *						  distribution to have infinit support)
 */
uint64 dist_uniform(void *par, int64 s)
{
	uint64	C;

	double	k;
	if( s<0 )
		k = (*((double *) par)) + s;
	else
		k = (*((double *) par)) - s;
	if( k<1 )
	{
		if( k<0 )
			return 0x7FFFFFFFFFFFFFFF;
		else
			return 0xFEFFFFFFFFFFFFFF;
	}

	if( s )
		C = round( k/(k+1)*18446744073609551616.0 );
	else
		C = round( k/(k+0.5)*18446744073609551616.0 );
	return	MAXUINT64 - C;
}

/* TRIANGULAR DISTRIBUTION
 *
 * Parameter type:	double *
 * Parameter value	{max}
 *
 *			/ /     |x| \  1
 *			| | 1 - --- | ---	, if |x| < max
 * PDF(x) = { \     max / max
 *			|
 *			\		<< 1		, otherwise (the probability for |x|>=max is not exactly 0 to allow this
 *								  distribution to have infinit support)
 */
uint64 dist_triang(void *par, int64 s)
{
	double	num, den;

	if( s>0 )
		num = *((double *) par) - s;
	else
		num = *((double *) par) + s;

	if( num<1 )
	{
		if( num<0 )
			return 0x7FFFFFFFFFFFFFFF;
		else
			return 0xFEFFFFFFFFFFFFFF;
	}

	if( s )
		den = num+0.5;
	else
		den = num;
	num = (num-0.5)/den;

	uint64	C = round( num*num*18446744073609551616.0 );
	return	MAXUINT64 - C;
}
