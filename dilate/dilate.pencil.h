// -*- c -*-
// UjoImro, 2013
// Experimental Code for the CARP Project

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

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
    int anchor_x,
    int anchor_y,
    int border_type );
    

#ifdef __cplusplus
} // extern "C"
#endif


// LuM end of file
