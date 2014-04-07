// -*- c -*-
// UjoImro, 2013
// Experimental Code for the CARP Project

#ifdef __cplusplus
extern "C" {
#endif

void pencil_affine_linear( const int src_rows, const int src_cols, const int src_step, const float src[]
                         , const int dst_rows, const int dst_cols, const int dst_step,       float dst[]
                         , const float a00, const float a01, const float a10, const float a11, const float b00, const float b10
                         );

#ifdef __cplusplus
} // extern "C"
#endif
