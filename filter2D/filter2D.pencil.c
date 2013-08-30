// -*- c -*-
// UjoImro, 2013
// Experimental Code for the CARP Project

#include <stdio.h>
#include <stdint.h>
#include <assert.h>

#include "filter2D.pencil.h"

void 
pencil_filter2D( 
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

    int center_row = kernel_rows / 2;
    int center_col = kernel_cols / 2;
   
     for ( int q = 0; q < rows; q++ ) 
     	for ( int w = 0; w < cols; w++ ) { 
     	    float prod = 0.;
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
		    prod += src[ row * src_step + col ] * kernel[ e * kernel_step + r ];
     	    	}
	    conv[ q * conv_step + w ] = prod;
     	} 

    return;
} // filter2D


// LuM end of file
