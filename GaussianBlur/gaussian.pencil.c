// -*- c -*-
// UjoImro, 2013
// Experimental Code for the CARP Project

#include <stdio.h>
#include <stdint.h>
#include <assert.h>

#include "gaussian.pencil.h"

static void 
filter2D(     
    int rows,
    int cols,
    int src_step,
    const float src[static const restrict src_step][cols],
    int kernel_rows,
    int kernel_cols,
    int kernel_step,
    const float kernel_[static const restrict kernel_step][kernel_cols],
    int conv_step, 
    float conv[static const restrict conv_step][cols] ) {

#pragma scop    
    int center_row = kernel_rows / 2;
    int center_col = kernel_cols / 2;
   
#   pragma independent
    for ( int q = 0; q < rows; q++ ) 
#       pragma independent
     	for ( int w = 0; w < cols; w++ ) { 
     	    float prod = 0.;
#           pragma reduction sum on prod
     	    for ( int e = 0; e < kernel_rows; e++ )
     	    	for ( int r = 0; r < kernel_cols; r++ )
	    	{
	    	    int row = q + e - center_row;
     	    	    int col = w + r - center_col;
		    row = row < 0 ? 0 : row;
		    row = row < rows ? row : rows - 1;
		    col = col < 0 ? 0 : col;
		    col = col < cols ? col : cols - 1;
		    if ( (row>=0) && (row<rows) && (col>=0) && (col<cols) )
		      prod += src[row][col] * kernel_[e][r];
     	    	}
	    conv[q][w] = prod;
     	}
#pragma endscop    
    return;

} // filter2D

static void
gaussian (    
    int rows,
    int cols,
    int src_step,
    const float src[static const restrict src_step][cols],
    int kernelX_rows,
    int kernelX_cols,
    int kernelX_step,
    const float kernelX[static const restrict kernelX_step][kernelX_cols],
    int kernelY_rows,
    int kernelY_cols,
    int kernelY_step,
    const float kernelY[static const restrict kernelY_step][kernelY_cols],
    int temp_step,
    float temp[static const restrict temp_step][cols],
    int conv_step, 
    float conv[static const restrict conv_step][cols] ) {
    filter2D( rows, cols, src_step,  src,  kernelX_rows, kernelX_cols, kernelX_step, kernelX, temp_step, temp );    
    filter2D( rows, cols, temp_step, temp, kernelY_rows, kernelY_cols, kernelY_step, kernelY, conv_step, conv );
    return;

} // gaussian

void 
pencil_gaussian( 
    int rows,
    int cols,
    int src_step,
    float src[],
    int kernelX_rows,
    int kernelX_cols,
    int kernelX_step,
    float kernelX[],
    int kernelY_rows,
    int kernelY_cols,
    int kernelY_step,
    float kernelY[],
    int temp_step,
    float temp[],
    int conv_step, 
    float conv[] ) {
    gaussian ( rows, cols, src_step, src, kernelX_rows, kernelX_cols, kernelX_step, kernelX, kernelY_rows, kernelY_cols, kernelY_step, kernelY, temp_step, temp, conv_step, conv );
    return;

} // gaussian


// LuM end of file
