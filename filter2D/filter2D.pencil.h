// -*- c -*-
// UjoImro, 2013
// Experimental Code for the CARP Project

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

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
        float conv[] );
    
#ifdef __cplusplus
} // extern "C"
#endif


// LuM end of file
