#include "resize.h"
#include <pencil.h>


static void resize( const int original_rows
                  , const int original_cols
                  , const int original_step
                  , const unsigned char original[static const restrict original_rows][original_step]
                  , const int resampled_rows
                  , const int resampled_cols
                  , const int resampled_step
                  , unsigned char resampled[static const restrict resampled_rows][resampled_step]
                  )
{
#pragma scop
    __pencil_assume(original_rows  >  0);
    __pencil_assume(original_cols  >  0);
    __pencil_assume(original_step  >= original_cols);
    __pencil_assume(resampled_rows >  0);
    __pencil_assume(resampled_cols >  0);
    __pencil_assume(resampled_step >= resampled_cols);

    __pencil_kill(resampled);
    {
        int o_h = original_rows;
        int o_w = original_cols;
        int n_h = resampled_rows;
        int n_w = resampled_cols;

        #pragma pencil independent
        for ( int n_r = 0; n_r < resampled_rows; n_r++ )
        {
            #pragma pencil independent
            for ( int n_c = 0; n_c < resampled_cols; n_c++ )
            {
                float o_r = ( n_r + 0.5f ) * (o_h) / (n_h) - 0.5f;
                float o_c = ( n_c + 0.5f ) * (o_w) / (n_w) - 0.5f;

                float r = o_r - floor(o_r);
                float c = o_c - floor(o_c);

                int coord_00_r = clamp( (int) floor(o_r), 0, o_h - 1 );
                int coord_00_c = clamp( (int) floor(o_c), 0, o_w - 1 );

                int coord_01_r = coord_00_r;
                int coord_01_c = clamp( coord_00_c + 1, 0, o_w - 1 );

                int coord_10_r = clamp( coord_00_r + 1, 0, o_h - 1 );
                int coord_10_c = coord_00_c;

                int coord_11_r = clamp( coord_00_r + 1, 0, o_h - 1 );
                int coord_11_c = clamp( coord_00_c + 1, 0, o_w - 1 );

                unsigned char A00 = original[coord_00_r][coord_00_c];
                unsigned char A10 = original[coord_10_r][coord_10_c];
                unsigned char A01 = original[coord_01_r][coord_01_c];
                unsigned char A11 = original[coord_11_r][coord_11_c];

                resampled[n_r][n_c] = mixf( mixf(A00, A10, r), mixf(A01, A11, r), c);
            }
        }
    }
    __pencil_kill(original);
#pragma endscop
}

void pencil_resize_LN( const int original_rows
                     , const int original_cols
                     , const int original_step
                     , const unsigned char original[]
                     , const int resampled_rows
                     , const int resampled_cols
                     , const int resampled_step
                     , unsigned char resampled[]
                     )
{
    resize(  original_rows,  original_cols,  original_step, (const unsigned char(*)[ original_step])original
          , resampled_rows, resampled_cols, resampled_step, (      unsigned char(*)[resampled_step])resampled
          );
}
