#include "gaussian.pencil.h"
#include <pencil.h>
#include <assert.h>

#if !__PENCIL__
#include <stdlib.h>
#endif

static void gaussian( const int rows
                    , const int cols
                    , const int step
                    , const float src[static const restrict rows][step]
                    , const int kernelX_length
                    , const float kernelX[static const restrict kernelX_length]
                    , const int kernelY_length
                    , const float kernelY[static const restrict kernelY_length]
                    , float conv[static const restrict rows][step]
                    )
{
#pragma scop
    __pencil_assume(rows         >  0);
    __pencil_assume(cols         >  0);
    __pencil_assume(step         >= cols);
    __pencil_assume(kernelX_length >  0);
    __pencil_assume(kernelX_length <= 128);
    __pencil_assume(kernelY_length >  0);
    __pencil_assume(kernelY_length <= 128);
    {
#if __PENCIL__
        float temp[rows][step];
#else
        float (*temp)[step] = (float (*)[step])malloc(sizeof(float)*rows*step);
#endif
        #pragma pencil independent
        for ( int q = 0; q < rows; q++ )
        {
            #pragma pencil independent
            for ( int w = 0; w < cols; w++ )
            {
                float prod = 0.;
                #pragma pencil independent reduction (+: prod);
                for ( int r = 0; r < kernelX_length; r++ )
                {
                    int row = q;
                    int col = clamp(w + r - kernelX_length / 2, 0, cols-1);
                    prod += src[row][col] * kernelX[r];
                }
                temp[q][w] = prod;
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
                for ( int e = 0; e < kernelY_length; e++ )
                {
                    int row = clamp(q + e - kernelY_length / 2, 0, rows-1);
                    int col = w;
                    prod += temp[row][col] * kernelY[e];
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
//    printf("%d %d %d / %d %d %d\n", kernelX_rows, kernelX_cols, kernelX_step, kernelY_rows, kernelY_cols, kernelY_step);
    assert( kernelX_rows == 1 );
    assert( kernelX_cols  > 0 );
    assert( kernelX_step == kernelX_cols );
    assert( kernelY_rows  > 0 );
    assert( kernelY_cols == 1 );
    assert( kernelY_step == 1 );
    gaussian( rows, cols, step, (const float(*)[step])src
            , kernelX_cols, kernelX
            , kernelY_rows, kernelY
            , (float(*)[step])conv
            );
}
