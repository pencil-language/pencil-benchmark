#include "histogram.pencil.h"
#include <pencil.h>

static void calcHist( const int rows
                    , const int cols
                    , const int step
                    , const uint8_t image[static const restrict rows][step]
                    , int hist[static const restrict HISTOGRAM_BINS]    //out
                    )
{
#pragma scop
    __pencil_assume(rows >  0);
    __pencil_assume(cols >  0);
    __pencil_assume(step >= cols);

    #pragma pencil independent
    for(int b = 0; b < HISTOGRAM_BINS; ++b)
        hist[b] = 0;

    #pragma pencil independent reduction(+:hist)
    for(int r = 0; r < rows; ++r)
    {
        #pragma pencil independent reduction(+:hist)
        for(int c = 0; c < cols; ++c)
        {
            uint8_t pixel = image[r*step+c];
            ++hist[pixel];
        }
    }
#pragma endscop
}

void pencil_calcHist( const int rows, const int cols, const int step, const uint8_t image[], int hist[HISTOGRAM_BINS])
{
    calcHist( rows, cols, step, (const uint8_t(*)[step])image, hist);
}