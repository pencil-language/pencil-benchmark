#include "gaussian.pencil.h"

static void gaussian( int rows
                    , int cols
                    , int step
                    , const float src[static const restrict rows][step]
                    , int kernelX_rows
                    , int kernelX_cols
                    , int kernelX_step
                    , const float kernelX[static const restrict kernelX_rows][kernelX_step]
                    , int kernelY_rows
                    , int kernelY_cols
                    , int kernelY_step
                    , const float kernelY[static const restrict kernelY_rows][kernelY_step]
                    , float temp[static const restrict rows][step]
                    , float conv[static const restrict rows][step]
                    )
{
#pragma scop
#   pragma pencil independent
    for ( int q = 0; q < rows; q++ )
    {
#       pragma pencil independent
        for ( int w = 0; w < cols; w++ )
        {
            float prod = 0.;
            for ( int e = 0; e < kernelX_rows; e++ )
            {
                for ( int r = 0; r < kernelX_cols; r++ )
                {
                    int row = q + e - kernelX_rows / 2;
                    int col = w + r - kernelX_cols / 2;
                    row = (row < 0   ) ? 0   : row;
                    row = (row < rows) ? row : rows - 1;
                    col = (col < 0   ) ? 0   : col;
                    col = (col < cols) ? col : cols - 1;
                    prod += src[row][col] * kernelX[e][r];
                }
                temp[q][w] = prod;
            }
        }
    }
#   pragma pencil independent
    for ( int q = 0; q < rows; q++ )
    {
#       pragma pencil independent
        for ( int w = 0; w < cols; w++ )
        {
            float prod = 0.;
            for ( int e = 0; e < kernelY_rows; e++ )
            {
                for ( int r = 0; r < kernelY_cols; r++ )
                {
                    int row = q + e - kernelY_rows / 2;
                    int col = w + r - kernelY_cols / 2;
                    row = (row < 0   ) ? 0   : row;
                    row = (row < rows) ? row : rows - 1;
                    col = (col < 0   ) ? 0   : col;
                    col = (col < cols) ? col : cols - 1;
                    prod += temp[row][col] * kernelY[e][r];
                }
                conv[q][w] = prod;
            }
        }
    }
#pragma endscop
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
