// -*- c -*-
// UjoImro, 2013
// Experimental Code for the CARP Project

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

void
pencil_affine_linear (
    int src_rows, int src_cols, int src_step, float src[],
    int dst_rows, int dst_cols, int dst_step, float dst[],
    float a00, float a01, float a10, float a11, float b00, float b10 );


#ifdef __cplusplus
} // extern "C"
#endif


// LuM end of file
