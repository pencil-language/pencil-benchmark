#ifndef HISTOGRAM_PENCIL_H
#define HISTOGRAM_PENCIL_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

    #define HISTOGRAM_BINS 256

    void pencil_calcHist( const int rows
                        , const int cols
                        , const int step
                        , const uint8_t image[]
                        , int hist[HISTOGRAM_BINS]    //out
                        );

#ifdef __cplusplus
} // extern "C"
#endif

#endif