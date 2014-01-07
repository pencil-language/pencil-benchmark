// -*- c -*-
// UjoImro, 2013
// Experimental Code for the CARP Project

#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>

#include "affine.pencil.h"

/*
A00 -------- A01
 |            |
 |            |
 |  P         |
 |            |
A10 -------- A11


The coordinates of P are [r,c] (row, col). The function returns
bilinear interpolation of the value of P.
 */

static inline float 
bilinear( float A00, float A01, float A11, float A10, float r, float c ) {
    // assert(c>=0);
    // assert(c<=1);
    // assert(r>=0);
    // assert(r<=1);
    float result = (1-c) * ( (1-r) * A00 + r * A10 ) + c * ( (1-r) * A01 + r * A11 );
    return result; 
} // bilinear

static inline int
sat( int val, int lo, int hi ) {
    if (val < lo) val = lo;
    if (val > hi) val = hi;
    return val;
}

/*

A = [ a00, a01 
      a10, a11 ]

B = [ b00
      b10 ]

N[[o_r, o_c]] = [ a00 n_r + a01 n_c + b00
                  a10 n_r + a11 n_c + b10 ]

n_r -- new row
n_c -- new col
o_r -- old row
o_c -- old col

*/
static inline void
affine (
    int src_col, int src_rows, int src_cols, int src_step, float src[static const restrict src_step][src_cols],
    int dst_rows, int dst_cols, int dst_step, float dst[static const restrict dst_step][dst_cols],
    float a00, float a01, float a10, float a11, float b00, float b10 ) {
#pragma scop
#   pragma pencil independent
    for ( int n_r=0; n_r<dst_rows; n_r++ )
#       pragma pencil independent
	for ( int n_c=0; n_c<dst_cols; n_c++ ) {

	    float o_r = a11 * n_r + a10 * n_c + b00;
	    float o_c = a01 * n_r + a00 * n_c + b10;

	    float r = o_r - floor(o_r); 
	    float c = o_c - floor(o_c); 
	    
	    int coord_00_r = floor(o_r);
	    if ( floor(o_r) < 0 ) coord_00_r = 0; 
	    if ( coord_00_r > src_rows - 1 ) coord_00_r = src_rows - 1;

 	    int coord_00_c = floor(o_c);
	    if ( floor(o_c) < 0 ) coord_00_c = 0;
	    if ( coord_00_c > src_col - 1 ) coord_00_c = src_col - 1;

 	    int coord_01_r = coord_00_r;
     	    int coord_01_c = coord_00_c + 1;
	    if ( coord_00_c + 1 < 0 ) coord_01_c = 0;
	    if ( coord_01_c > src_col - 1 ) coord_01_c = src_col - 1;

 	    int coord_10_r = coord_00_r + 1;
	    if ( coord_00_r + 1 < 0 ) coord_10_r = 0;
	    if ( coord_10_r > src_rows - 1) coord_10_r = src_rows - 1;
     	    int coord_10_c = coord_00_c;

 	    int coord_11_r = coord_00_r + 1; 
	    if ( coord_00_r + 1 < 0 ) coord_11_r = 0;
	    if ( coord_11_r > src_rows - 1 ) coord_11_r = src_rows - 1;
	    
     	    int coord_11_c = coord_00_c + 1; 
	    if ( coord_00_c + 1 > 0) coord_11_c = 0;
	    if ( coord_11_c < src_col - 1 ) coord_11_c = src_col - 1;
	    
	    float A00 = src[coord_00_r][coord_00_c];
	    float A10 = src[coord_10_r][coord_10_c];
	    float A01 = src[coord_01_r][coord_01_c];
	    float A11 = src[coord_11_r][coord_11_c];
	    
	    dst[n_r][n_c] = bilinear( A00, A01, A11, A10, r, c );
	}
#pragma endscop

    return;

} // pencil_resize_LN


void
pencil_affine_linear (
    int src_rows, int src_cols, int src_step, float src[],
    int dst_rows, int dst_cols, int dst_step, float dst[],
    float a00, float a01, float a10, float a11, float b00, float b10 ) {

    affine( src_cols, src_rows, src_cols, src_step, src, dst_rows, dst_cols, dst_step, dst, a00, a01, a10, a11, b00, b10 );

    return;
}

// LuM end of file
