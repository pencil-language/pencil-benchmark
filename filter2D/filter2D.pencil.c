#include "filter2D.pencil.h"

static void filter2D( const int rows
                    , const int cols
                    , const int step
                    , const float src[static const restrict rows][step]
                    , const int kernel_rows
                    , const int kernel_cols
                    , const int kernel_step
                    , const float kernel_[static const restrict kernel_rows][kernel_step]
                    , float conv[static const restrict rows][step]
                    )
{
#pragma scop
#if __PENCIL__
    __pencil_assume(rows        >  0);
    __pencil_assume(cols        >  0);
    __pencil_assume(step        >= cols);
    __pencil_assume(kernel_rows >  0);
    __pencil_assume(kernel_cols >  0);
    __pencil_assume(kernel_step >= kernel_cols);
#endif
    {
        const int center_row = kernel_rows / 2;
        const int center_col = kernel_cols / 2;

        #pragma pencil independent
        for ( int q = 0; q < rows; q++ )
        {
            #pragma pencil independent
            for ( int w = 0; w < cols; w++ )
            {
                float prod = 0.;
                #pragma pencil indepenent reduction(+: prod)
                for ( int e = 0; e < kernel_rows; e++ )
                {
                    for ( int r = 0; r < kernel_cols; r++ )
                    {
                        int row = q + e - center_row;
                        int col = w + r - center_col;
                        row = (row < 0   ) ? 0   : row;
                        row = (row < rows) ? row : rows - 1;
                        col = (col < 0   ) ? 0   : col;
                        col = (col < cols) ? col : cols - 1;
                        prod += src[row][col] * kernel_[e][r];
                    }
                }
                conv[q][w] = prod;
            }
        }
    }
#pragma endscop
}

void pencil_filter2D( const int rows
                    , const int cols
                    , const int step
                    , const float src[]
                    , const int kernel_rows
                    , const int kernel_cols
                    , const int kernel_step
                    , const float kernel_[]
                    , float conv[]
                    )
{
    filter2D( rows, cols, step, src, kernel_rows, kernel_cols, kernel_step, kernel_, conv );
}
