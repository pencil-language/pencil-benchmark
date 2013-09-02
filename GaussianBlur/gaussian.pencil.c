// -*- c -*-
// UjoImro, 2013
// Experimental Code for the CARP Project

#include <stdio.h>
#include <stdint.h>
#include <assert.h>

#include "gaussian.pencil.h"

static void 
gaussian(     
    int rows,
    int cols,
    int src_step,
    const float src[static const restrict src_step][cols],
    int kernel_rows,
    int kernel_cols,
    int kernel_step,
    const float kernel[static const restrict kernel_step][kernel_cols],
    int conv_step, 
    float conv[static const restrict conv_step][cols] ) {
    
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
	    	    // if ( (row>=0) && (row<rows) && (col>=0) && (col<cols) )
		    prod += src[row][col] * kernel[e][r];
     	    	}
	    conv[q][w] = prod;
     	}
    
    return;
} // gaussian 

void 
pencil_gaussian( 
    int rows,
    int cols,
    int src_step,
    float src[],
    int kernel_rows,
    int kernel_cols,
    int kernel_step,
    float kernel[],
    int conv_step, 
    float conv[] ) {
    
    gaussian( rows, cols, src_step, src, kernel_rows, kernel_cols, kernel_step, kernel, conv_step, conv );

    return;
} // gaussian


// LuM end of file
