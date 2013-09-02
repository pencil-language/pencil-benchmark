// -*- c -*-
// UjoImro, 2013
// Experimental Code for the CARP Project

#include <stdio.h>
#include <stdint.h>
#include <assert.h>

#include "dilate.pencil.h"

static void
dilate (
    int rows,
    int cols,
    int cpu_step,
    const uint8_t cpu_gray[static const restrict rows][cpu_step],
    int dilate_step,
    uint8_t dilate[static const restrict rows][dilate_step],
    int se_rows,
    int se_cols,
    int se_step,
    const uint8_t se[static const restrict se_rows][se_step],
    int anchor_row,
    int anchor_col,
    int border_type ) {

#   pragma independent
    for ( int q = 0; q < rows; q++ ) {
#       pragma independent
	for ( int w = 0; w < cols; w++ ) {
	    uint8_t sup = 0;
#           pragma reduce max on sup
	    for ( int e = 0; e < se_rows; e++ ) {
		for ( int r = 0; r < se_cols; r++ ) {
		    int candidate_row = q - anchor_row + e;
		    int candidate_col = w - anchor_col + r;
		    int low_row  = (candidate_row >= 0);
		    int low_col  = (candidate_col >= 0);
		    int high_row = (candidate_row < rows);
		    int high_col = (candidate_col < cols);
		    if ( (low_row) && (low_col) && (high_row) && (high_col) ) {
			uint8_t val = cpu_gray[candidate_row][candidate_col];
			if (se[e][r]!=0)
			    sup = (sup > val) ? sup : val;
		    } // if
		} // r
	    } // e
	    dilate[q][w] = sup;
	} // w
    } // q 
    
} // dilate

void
pencil_dilate (
    int rows,
    int cols,
    int cpu_step,
    uint8_t cpu_gray[],
    int dilate_step,
    uint8_t pdilate[],
    int se_rows,
    int se_cols,
    int se_step,
    uint8_t se[],
    int anchor_row,
    int anchor_col,
    int border_type ) {

    assert(border_type==0);

    dilate( rows, cols, cpu_step, cpu_gray, dilate_step, pdilate, se_rows, se_cols, se_step, se, anchor_row, anchor_col, border_type );

    return;
} // pencil_dilate



// LuM end of file
