#include "gaussian.pencil.h"

static void filter2D( int rows
                    , int cols
                    , int step
                    , const float src[static const restrict rows][step]
                    , int kernel_rows
                    , int kernel_cols
                    , int kernel_step
                    , const float kernel_[static const restrict kernel_rows][kernel_step]
                    , float conv[static const restrict step][cols]
                    )
{
#pragma scop
    int center_row = kernel_rows / 2;
    int center_col = kernel_cols / 2;

#   pragma pencil independent
    for ( int q = 0; q < rows; q++ )
#       pragma pencil independent
        for ( int w = 0; w < cols; w++ ) {
            float prod = 0.;
            for ( int e = 0; e < kernel_rows; e++ )
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
            conv[q][w] = prod;
        }
#pragma endscop
}

static void gaussian( int rows
                    , int cols
                    , int step
                    , const float src[static const restrict step][cols]
                    , int kernelX_rows
                    , int kernelX_cols
                    , int kernelX_step
                    , const float kernelX[static const restrict kernelX_step][kernelX_cols]
                    , int kernelY_rows
                    , int kernelY_cols
                    , int kernelY_step
                    , const float kernelY[static const restrict kernelY_step][kernelY_cols]
                    , float temp[static const restrict step][cols]
                    , float conv[static const restrict step][cols]
                    )
{
    filter2D( rows, cols, step,  src, kernelX_rows, kernelX_cols, kernelX_step, kernelX, temp );
    filter2D( rows, cols, step, temp, kernelY_rows, kernelY_cols, kernelY_step, kernelY, conv );
}

void pencil_gaussian( int rows
                    , int cols
                    , int step
                    , float src[]
                    , int kernelX_rows
                    , int kernelX_cols
                    , int kernelX_step
                    , float kernelX[]
                    , int kernelY_rows
                    , int kernelY_cols
                    , int kernelY_step
                    , float kernelY[]
                    , float temp[]
                    , float conv[]
                    )
{
    gaussian ( rows, cols, step, src, kernelX_rows, kernelX_cols, kernelX_step, kernelX, kernelY_rows, kernelY_cols, kernelY_step, kernelY, temp, conv );
}
