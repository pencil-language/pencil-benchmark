#include <stdint.h>

#include "cvt_color.pencil.h"

static void RGB2Gray( const int rows
                    , const int cols
                    , const int src_step
                    , const int dst_step
                    , const DATA_TYPE src[static const restrict rows][src_step][3]
                    , DATA_TYPE dst[static const restrict rows][dst_step]
                    )
{
#pragma scop
#if __PENCIL__
    __pencil_assume(rows     >  0);
    __pencil_assume(cols     >  0);
    __pencil_assume(src_step >= cols);
    __pencil_assume(dst_step >= cols);
#endif
    #pragma pencil independent
    for ( int q = 0; q < rows; q++ )
    {
        #pragma pencil independent
        for( int w = 0; w < cols; w++ )
        {
            dst[q][w] = CV_DESCALE( (src[q][w][2] * B2Y + src[q][w][1] * G2Y + src[q][w][0] * R2Y ), yuv_shift );
        }
    }
#pragma endscop
}

void pencil_RGB2Gray( const int rows
                    , const int cols
                    , const int src_step
                    , const int dst_step
                    , const DATA_TYPE src[]
                    , DATA_TYPE dst[]
                    )
{
    RGB2Gray( rows, cols, src_step, dst_step, src, dst );
}
