#include "filter2D.pencil.h"

#include <pencil.h>

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
    __pencil_assume(kernel_rows <=  8);
    __pencil_assume(kernel_cols <=  8);
    __pencil_assume(kernel_rows >=  3);
    __pencil_assume(kernel_cols >=  3);
    __pencil_assume(   1 <= rows);
    __pencil_assume(   1 <= cols);
    __pencil_assume(cols <= step);
    __pencil_assume(   1 <= kernel_rows);
    __pencil_assume(   1 <= kernel_cols);
    __pencil_assume(cols <= kernel_step);
    {
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
                        int row = iclampi( q + e - kernel_rows / 2, 0, rows-1 );
                        int col = iclampi( w + r - kernel_cols / 2, 0, cols-1 );
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
    filter2D(        rows,        cols,        step, (const float (*)[       step])src
            , kernel_rows, kernel_cols, kernel_step, (const float (*)[kernel_step])kernel_
            ,                                        (      float (*)[       step])conv
            );
}
