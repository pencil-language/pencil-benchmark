// -*- c -*-
// UjoImro, 2013
// Experimental Code for the CARP Project

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

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
	float temp[],
	float conv[] );
    

#ifdef __cplusplus
} // extern "C"
#endif


// LuM end of file
