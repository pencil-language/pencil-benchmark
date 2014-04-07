// -*- c -*-
// UjoImro, 2013
// Experimental Code for the CARP Project

#ifndef CVT_COLOR_PENCIL__H__
#define CVT_COLOR_PENCIL__H__

#ifdef __cplusplus
extern "C" {
#endif

#define DEPTH_0

#if defined (DEPTH_0)
#undef  DATA_TYPE
#define DATA_TYPE uint8_t
#define MAX_NUM  255
#define HALF_MAX 128
#define SAT_CAST(num) convert_uchar_sat(num)
#endif

#if defined (DEPTH_2)
#undef  DATA_TYPE
#define DATA_TYPE ushort
#define MAX_NUM  65535
#define HALF_MAX 32768
#define SAT_CAST(num) convert_ushort_sat(num)
#endif

#if defined (DEPTH_5)
#undef  DATA_TYPE
#define DATA_TYPE float
#define MAX_NUM  1.0f
#define HALF_MAX 0.5f
#define SAT_CAST(num) (num)
#endif

#define CV_DESCALE(x,n) (((x) + (1 << ((n)-1))) >> (n))
    enum
    {
        yuv_shift  = 14,
        xyz_shift  = 12,
        R2Y        = 4899,
        G2Y        = 9617,
        B2Y        = 1868,
        BLOCK_SIZE = 256
    };


    void pencil_RGB2Gray( const int rows
                        , const int cols
                        , const int src_step
                        , const int dst_step
                        , const DATA_TYPE src[]
                        , DATA_TYPE dst[]
                        );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* CVT_COLOR_PENCIL__H__ */
// LuM end of file
