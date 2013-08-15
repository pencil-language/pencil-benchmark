// -*- c -*-
// UjoImro, 2013
// Experimental Code for the CARP Project

#include <stdint.h>

#include "cvt_color.pencil.h"

// this function cannot be called from OpenCV
static 
void 
pencil_RGB2Gray_caller(
    int cols,
    int rows,
    int src_step,
    int dst_step,
    int channels,
    int bidx,
    const DATA_TYPE src[static const restrict src_step][cols][channels],
    DATA_TYPE dst[static const restrict dst_step][cols] )
{
#   pragma independent
    for ( int q = 0; q<rows; q++ )
#       pragma independent
	for( int w = 0; w<cols; w++ )
	    dst[q][w] = CV_DESCALE( (src[q][w][bidx] * B2Y + src[q][w][1] * G2Y + src[q][w][(bidx^2)] * R2Y ), yuv_shift );
    return;
} // pencil_RGB2Gray_caller


void 
pencil_RGB2Gray( 
    int cols, 
    int rows, 
    int src_step,
    int dst_step, 
    int channels,
    int bidx, 
    const DATA_TYPE src[], 
    DATA_TYPE dst[] )
{
    pencil_RGB2Gray_caller( cols, rows, src_step, dst_step, channels, bidx, src, dst );

    return;

/* #   pragma independent */
/*     for ( int q = 0; q<cols; q++ ) */
/* #       pragma independent */
/* 	for( int w=0; w<rows; w++ ) { */
/* 	    int src_idx = w * src_step + q * channels;  */
/* 	    int dst_idx = w * dst_step + q;  */

/* 	    dst[dst_idx] = (DATA_TYPE)CV_DESCALE((src[src_idx + bidx] * B2Y + src[src_idx + 1] * G2Y + src[src_idx + (bidx^2)] * R2Y), yuv_shift);  */
/* 	} */
/*     return; */
} // pencil_RGB2Gray



// LuM end of file
