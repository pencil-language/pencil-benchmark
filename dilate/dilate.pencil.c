#include <stdio.h>
#include <stdint.h>
#include <assert.h>

#include "dilate.pencil.h"

unsigned char max(unsigned char A, unsigned char B)
{
	return (A>B)?A:B;
}

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
#if __PENCIL__
    __pencil_assume(rows        >  0);
    __pencil_assume(cols        >  0);
    __pencil_assume(step        >= cols);
    __pencil_assume(dilate_step >= cols);
    __pencil_assume(se_rows     >  0);
    __pencil_assume(se_cols     >  0);
    __pencil_assume(se_step     >= se_cols);
    __pencil_assume(anchor_row  >= 0);
    __pencil_assume(anchor_col  >= 0);
#endif
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
                   int candidate_row = q - anchor_row + e;
                   int candidate_col = w - anchor_col + r;

		   int test = (((candidate_row >= 0) && (candidate_row < rows) && (candidate_col >= 0) && (candidate_col < cols)) && ((se[e][r]!=0)));
		   sup = test ? max(sup, cpu_gray[candidate_row][candidate_col]) : sup;
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
    dilate( rows, cols, cpu_step, cpu_gray, dilate_step, pdilate, se_rows, se_cols, se_step, se, anchor_row, anchor_col );
}
