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

static inline double 
bilinear( double A00, double A01, double A11, double A10, double r, double c ) {
    assert(c>=0);
    assert(c<=1);
    assert(r>=0);
    assert(r<=1);

    return (1-c) * ( (1-r) * A00 + r * A10 ) + c * ( (1-r) * A01 + r * A11 );
} // bilinear

static inline int
sat( int val, int lo, int hi ) {
    val = (val >= lo) ? val : lo;
    val = (val <= hi) ? val : hi;
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
void
pencil_affine_linear (
    int src_rows, int src_cols, int src_step, float src[],
    int dst_rows, int dst_cols, int dst_step, float dst[],
    float a00, float a01, float a10, float a11, float b00, float b10 ) {

    for ( int n_r=0; n_r<dst_rows; n_r++ )
	for ( int n_c=0; n_c<dst_cols; n_c++ ) {

	    double o_r = a11 * n_r + a10 * n_c + b00;
	    double o_c = a01 * n_r + a00 * n_c + b10;
	    // printf( "o_r = %f\n", o_r);
	    // printf( "o_c = %f\n", o_c);

	    double r = o_r - floor(o_r); 
	    double c = o_c - floor(o_c); 
	    
	    int coord_00_r = sat( floor(o_r), 0, src_rows - 1 );
 	    int coord_00_c = sat( floor(o_c), 0, src_cols - 1 );

 	    int coord_01_r = coord_00_r;
     	    int coord_01_c = sat( coord_00_c + 1, 0, src_cols - 1 );

 	    int coord_10_r = sat( coord_00_r + 1, 0, src_rows - 1 );
     	    int coord_10_c = coord_00_c;

 	    int coord_11_r = sat( coord_00_r + 1, 0, src_rows - 1 );
     	    int coord_11_c = sat( coord_00_c + 1, 0, src_cols - 1 );
	    
	    double A00 = src[coord_00_r * src_step + coord_00_c];
	    double A10 = src[coord_10_r * src_step + coord_10_c];
	    double A01 = src[coord_01_r * src_step + coord_01_c];
	    double A11 = src[coord_11_r * src_step + coord_11_c];
	    
	    dst[ n_r * dst_step + n_c ] = bilinear( A00, A01, A11, A10, r, c );
	}

    return;
} // pencil_resize_LN


// LuM end of file
