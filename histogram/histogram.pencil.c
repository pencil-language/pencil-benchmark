#include "histogram.pencil.h"

void pencil_calcHist( const int rows
                    , const int cols
                    , const int step
                    , const uint8_t image[]
                    , int hist[HISTOGRAM_BINS]    //out
                    )
{
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
}
