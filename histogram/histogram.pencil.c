#include "histogram.pencil.h"

void pencil_calcHist( const int rows
                    , const int cols
                    , const int step
                    , const uint8_t image[]
                    , int hist[HISTOGRAM_BINS]    //out
                    )
{
    for(int b = 0; b < HISTOGRAM_BINS; ++b)
        hist[b] = 0;
    
    for(int r = 0; r < rows; ++r)
    {
        for(int c = 0; c < cols; ++c)
        {
            uint8_t pixel = image[r*step+c];
            ++hist[pixel];
        }
    }
}
