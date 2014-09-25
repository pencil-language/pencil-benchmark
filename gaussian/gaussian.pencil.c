#include "gaussian.pencil.h"
#include <pencil.h>

#if !__PENCIL__
#include <stdlib.h>
#endif

static void gaussian( const int rows
                    , const int cols
                    , const int step
                    , const float src[static const restrict rows][step]
                    , const int kernelX_rows
                    , const int kernelX_cols
                    , const int kernelX_step
                    , const float kernelX[static const restrict kernelX_rows][kernelX_step]
                    , const int kernelY_rows
                    , const int kernelY_cols
                    , const int kernelY_step
                    , const float kernelY[static const restrict kernelY_rows][kernelY_step]
                    , float conv[static const restrict rows][step]
                    )
{
#pragma scop
#if __PENCIL__
    __pencil_assume(rows         >  0);
    __pencil_assume(cols         >  0);
    __pencil_assume(step         >= cols);
    __pencil_assume(kernelX_rows >  0);
    __pencil_assume(kernelX_cols >  0);
    __pencil_assume(kernelX_step >= kernelX_cols);
    __pencil_assume(kernelY_rows >  0);
    __pencil_assume(kernelY_cols >  0);
    __pencil_assume(kernelY_step >= kernelY_cols);
    __pencil_assume(kernelX_rows <= 2);
    __pencil_assume(kernelX_cols <= 128);
    __pencil_assume(kernelY_rows <= 128);
    __pencil_assume(kernelY_cols <= 2);
#endif
    {
#if __PENCIL__
        float temp[rows][step];
#else
        float* temp = (float*)malloc(sizeof(float)*rows*step);
#endif
        #pragma pencil independent
        for ( int q = 0; q < rows; q++ )
        {
            #pragma pencil independent
            for ( int w = 0; w < cols; w++ )
            {
                float prod = 0.;
                #pragma pencil independent reduction (+: prod);
                for ( int e = 0; e < kernelX_rows; e++ )
                {
                    for ( int r = 0; r < kernelX_cols; r++ )
                    {
                        int row = clamp(q + e - kernelX_rows / 2, 0, rows-1);
                        int col = clamp(w + r - kernelX_cols / 2, 0, cols-1);
                        prod += src[row][col] * kernelX[e][r];
                    }
                }
#if __PENCIL__
                temp[q][w] = prod;
#else
                temp[q*step+w] = prod;
#endif
            }
        }
        #pragma pencil independent
        for ( int q = 0; q < rows; q++ )
        {
            #pragma pencil independent
            for ( int w = 0; w < cols; w++ )
            {
                float prod = 0.;
                #pragma pencil independent reduction (+: prod);
                for ( int e = 0; e < kernelY_rows; e++ )
                {
                    for ( int r = 0; r < kernelY_cols; r++ )
                    {
                        int row = clamp(q + e - kernelY_rows / 2, 0, rows-1);
                        int col = clamp(w + r - kernelY_cols / 2, 0, cols-1);
#if __PENCIL__
                        prod += temp[row][col] * kernelY[e][r];
#else
                        prod += temp[row*step+col] * kernelY[e][r];
#endif
                    }
                }
                conv[q][w] = prod;
            }
        }
#if !__PENCIL__
        free(temp);
#endif
    }
#pragma endscop
}

void pencil_gaussian( const int rows
                    , const int cols
                    , const int step
                    , const float src[]
                    , const int kernelX_rows
                    , const int kernelX_cols
                    , const int kernelX_step
                    , const float kernelX[]
                    , const int kernelY_rows
                    , const int kernelY_cols
                    , const int kernelY_step
                    , const float kernelY[]
                    , float conv[]
                    )
{
    gaussian ( rows, cols, step, src, kernelX_rows, kernelX_cols, kernelX_step, kernelX, kernelY_rows, kernelY_cols, kernelY_step, kernelY, conv );
}
