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

#define bilinear( A00, A01, A11, A10, r, c ) ((1-c) * ( (1-r) * A00 + r * A10 ) + c * ( (1-r) * A01 + r * A11 ))

#define sat( val, lo, hi ) ((val >= lo) ? ((val <= hi) ? val : hi) : lo)

static void
resize (
    int original_col,
    int original_rows,
    int original_cols,
    int original_step,
    unsigned char original[static const restrict original_step][original_cols],
    int resampled_rows,
    int resampled_cols,
    int resampled_step,
    unsigned char resampled[static const restrict resampled_step][resampled_cols] ) {

#pragma scop
    // assert(resampled_rows>1);
    // assert(resampled_cols>1);

    float o_h = original_rows;
    float o_w = original_col;
    float n_h = resampled_rows;
    float n_w = resampled_cols;
#   pragma pencil independent
    for ( int n_r = 0; n_r < resampled_rows; n_r++ )
#       pragma pencil independent
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
	    
	    unsigned char A00 = original[coord_00_r][coord_00_c];
	    unsigned char A10 = original[coord_10_r][coord_10_c];
	    unsigned char A01 = original[coord_01_r][coord_01_c];
	    unsigned char A11 = original[coord_11_r][coord_11_c];
	    
	    resampled[n_r][n_c] = bilinear( A00, A01, A11, A10, r, c );
     	} // for n_c, n_r

#pragma endscop
    return;
} // pencil_resize_LN

void
pencil_resize_LN (
    int original_rows,  int original_cols,  int original_step,  unsigned char original[],
    int resampled_rows, int resampled_cols, int resampled_step, unsigned char resampled[] ) {

    resize(original_cols, original_rows, original_cols, original_step, original, 
	    resampled_rows, resampled_cols, resampled_step, resampled );
    return;
} // pencil_resize_LN


// LuM end of file
