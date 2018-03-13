#ifndef PROBABILITY_H
#define PROBABILITY_H

#include "file.h"
#include <math.h>

#define _k_SQRT2	1.4142135623730950488
#define _k_PI		3.141592653589793

/* This file contains functions to be employed in the arithmetic coder to make it encode the symbols with the given distribution
 *
 * All functions receive two parameters:
 *	void *par:	is a pointer to a vector containing the parameters of the given distribution
 *	int64 s:	is the symbol
 *
 * The type and meaning of the first parameter depends on the distribution beeing employed. See table below for more details:
 *
 * | DISTRIBUTION	| PARAMETER TYPE	| PARAMETER VALUE			| PDF(x) (PROBABILITY DENSITY FUNCTION)
 * |----------------|-------------------|---------------------------|-------------------------------------------------------------------------------|
 * |				|					|							|   1      /   |x| \            __												|
 * | LAPLACE		| double *			| {sd}	standard deviation	| ----- exp| - --- |, where b*\/2  = sd											|
 * |				|					|							|  2*b     \    b  /															|
 * |----------------|-------------------|---------------------------|-------------------------------------------------------------------------------|
 * |				|					|							| 	 1         /    x^2   \														|
 * | NORMAL			| double *			| {sd}	standard deviation	| -----____ exp| - ------ |														|
 * |				|					|							| sd \/2 pi    \   2 sd^2 /														|
 * |----------------|-------------------|---------------------------|-------------------------------------------------------------------------------|
 * |				|					|							|	  1      /    gamma^2    \													|
 * | CAUCHY			| double *			| {gamma}					| ---------- | ------------- |													|
 * |				|					|							|  pi*gamma  \ x^2 + gamma^2 /													|
 * |----------------|-------------------|---------------------------|-------------------------------------------------------------------------------|
 * |				|					|							|	  1																			|
 * | UNIFORM		| double *			| {max}	max value			| --------- , if |x| <= max and << 1 otherwise (to allow infinit support)		|
 * |				|					|							| 2 max + 1																		|
 * |----------------|-------------------|---------------------------|-------------------------------------------------------------------------------|
 * |				|					|							| /     |x| \  1																|
 * | TRIANGULAR		| double *			| {max}	max value			| | 1 - --- | --- , if |x| <= max and << 1 otherwise (to allow infinit support) |
 * |				|					|							| \     max / max																|
 * '----------------'-------------------'---------------------------'-------------------------------------------------------------------------------'
 *
 *
 * The value returned by this functions is the normalized probability expressed in fixed point, given by P(s) in the equation:
 *
 *          /     C(s+1)
 *          | 1 - ------	, if s>0
 *          |      C(s)
 *          |
 * P(s)+1   |     C(s-1)
 * ------ = { 1 - ------	, if s<0
 *  2^64    |      C(s)
 *          |
 *          \  C(0) - C(1)	, otherwise
 *
 * where
 *
 *        /  /'infinity
 *        |  |			PDF(x) dx	, if s >= 0
 *        | '/ (s-1/2)
 * C(s) = {
 *        |  /'(s+1/2)
 *        |  |			PDF(x) dx	, otherwise
 *        \ '/-infinity
 */

uint64	dist_laplace	(void *par, int64  );
uint64	dist_normal		(void *par, int64 s);
uint64	dist_cauchy		(void *par, int64 s);
uint64	dist_uniform	(void *par, int64 s);
uint64	dist_triang		(void *par, int64 s);

#endif // PROBABILITY_H
