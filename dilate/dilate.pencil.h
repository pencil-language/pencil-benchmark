// -*- c -*-
// UjoImro, 2013
// Experimental Code for the CARP Project

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

void
pencil_dilate (
    uint8_t cpu_gray[],
    int rows,
    int cols,
    int cpu_step,
    uint8_t dilate[],
    int dilate_step,
    uint8_t se[],
    int se_rows,
    int se_cols,
    int se_step,
    int anchor_x,
    int anchor_y,
    int border_type,
    int size_w,
    int size_h );
    

#ifdef __cplusplus
} // extern "C"
#endif


// LuM end of file
