#include "dilate.pencil.h"

#include <pencil.h>

static void dilate( const int rows
                  , const int cols
                  , const int step
                  , const unsigned char cpu_gray[static const restrict rows][step]
                  , const int dilate_step
                  , unsigned char dilate[static const restrict rows][dilate_step]
                  , const int se_rows
                  , const int se_cols
                  , const int se_step
                  , const unsigned char se[static const restrict se_rows][se_step]
                  , const int anchor_row
                  , const int anchor_col
                  )
{
#pragma scop
    __pencil_assume(se_rows    <= 9);
    __pencil_assume(se_cols    <= 9);
    __pencil_assume(anchor_row >= 0);
    __pencil_assume(anchor_row < se_rows);
    __pencil_assume(anchor_col >= 0);
    __pencil_assume(anchor_col < se_cols);
    #pragma pencil independent
    for ( int q = 0; q < rows; q++ )
    {
        #pragma pencil independent
        for ( int w = 0; w < cols; w++ )
        {
            unsigned char sup = 0;
            #pragma pencil independent reduction(max: sup)
            for ( int e = 0; e < se_rows; e++ )
            {
                for ( int r = 0; r < se_cols; r++ )
                {
                   int candidate_row = clampi(q - anchor_row + e, 0, rows);
                   int candidate_col = clampi(w - anchor_col + r, 0, cols);

		   sup = (se[e][r]!=0) ? max(sup, cpu_gray[candidate_row][candidate_col]) : sup;
                }
            }
            dilate[q][w] = sup;
        }
    }
#pragma endscop
}

void pencil_dilate( const int rows
                  , const int cols
                  , const int cpu_step
                  , const uint8_t cpu_gray[]
                  , const int dilate_step
                  , uint8_t pdilate[]
                  , const int se_rows
                  , const int se_cols
                  , const int se_step
                  , const uint8_t se[]
                  , const int anchor_row
                  , const int anchor_col
                  )
{
    dilate( rows, cols
          , cpu_step   , (const unsigned char(*)[cpu_step   ])cpu_gray
          , dilate_step, (      unsigned char(*)[dilate_step])pdilate
          , se_rows, se_cols
          , se_step    , (const unsigned char(*)[se_step    ])se
          , anchor_row, anchor_col
          );
}
