#include "histogram.pencil.h"
#include <pencil.h>

void atomic_inc(int *v);

#ifndef __PENCIL__
void atomic_inc(int *v)
{
    (*v)++;
}
#endif

static void calcHist( const int rows
                    , const int cols
                    , const int step
                    , const unsigned char image[static const restrict rows][step]
                    , int hist[static const restrict HISTOGRAM_BINS]    //out
                    )
{
#pragma scop
    __pencil_assume(rows >  0);
    __pencil_assume(cols >  0);
    __pencil_assume(step >= cols);

    __pencil_kill(hist);

    #pragma pencil independent
    for(int b = 0; b < HISTOGRAM_BINS; ++b)
        hist[b] = 0;

    #pragma pencil independent
    for(int r = 0; r < rows; ++r)
    {
        #pragma pencil independent
        for(int c = 0; c < cols; ++c)
        {
            unsigned char pixel = image[r][c];
            atomic_inc(&hist[pixel]);
        }
    }

    __pencil_kill(image);
#pragma endscop
}

void pencil_calcHist( const int rows, const int cols, const int step, const unsigned char image[], int hist[HISTOGRAM_BINS])
{
    calcHist( rows, cols, step, (const unsigned char(*)[step])image, hist);
}
