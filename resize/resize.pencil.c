// -*- c -*-
// UjoImro, 2013
// Experimental Code for the CARP Project

#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>

#include "resize.pencil.h"

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
    val = (val >= lo) ? val : lo;
    val = (val <= hi) ? val : hi;
    return val;
}

static void
resize (
    int original_rows,
    int original_cols,
    int original_step,
    uint8_t original[static const restrict original_step][original_cols],
    int resampled_rows,
    int resampled_cols,
    int resampled_step,
    uint8_t resampled[static const restrict resampled_step][resampled_cols] ) {

#pragma scop
    // assert(resampled_rows>1);
    // assert(resampled_cols>1);

    float o_h = original_rows;
    float o_w = original_cols;
    float n_h = resampled_rows;
    float n_w = resampled_cols;
#   pragma independent
    for ( int n_r = 0; n_r < resampled_rows; n_r++ )
#       pragma independent
     	for ( int n_c = 0; n_c < resampled_cols; n_c++ ) {
	    float o_r = ( n_r + 0.5 ) * (o_h) / (n_h) - 0.5;
 	    float o_c = ( n_c + 0.5 ) * (o_w) / (n_w) - 0.5;

 	    float r = o_r - floor(o_r); 
	    float c = o_c - floor(o_c); 
	    
	    int coord_00_r = sat( floor(o_r), 0, o_h - 1 );
 	    int coord_00_c = sat( floor(o_c), 0, o_w - 1 );

 	    int coord_01_r = coord_00_r;
     	    int coord_01_c = sat( coord_00_c + 1, 0, o_w - 1 );

 	    int coord_10_r = sat( coord_00_r + 1, 0, o_h - 1 );
     	    int coord_10_c = coord_00_c;

 	    int coord_11_r = sat( coord_00_r + 1, 0, o_h - 1 );
     	    int coord_11_c = sat( coord_00_c + 1, 0, o_w - 1 );
	    
	    uint8_t A00 = original[coord_00_r][coord_00_c];
	    uint8_t A10 = original[coord_10_r][coord_10_c];
	    uint8_t A01 = original[coord_01_r][coord_01_c];
	    uint8_t A11 = original[coord_11_r][coord_11_c];
	    
	    resampled[n_r][n_c] = bilinear( A00, A01, A11, A10, r, c );
     	} // for n_c, n_r

#pragma endscop
    return;
} // pencil_resize_LN

void
pencil_resize_LN (
    int original_rows,  int original_cols,  int original_step,  uint8_t original[],
    int resampled_rows, int resampled_cols, int resampled_step, uint8_t resampled[] ) {

    resize( original_rows, original_cols, original_step, original, 
	    resampled_rows, resampled_cols, resampled_step, resampled );
    return;
} // pencil_resize_LN


// LuM end of file
