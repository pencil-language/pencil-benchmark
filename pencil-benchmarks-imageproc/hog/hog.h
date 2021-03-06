#ifndef HOG_PENCIL_H
#define HOG_PENCIL_H

#include <pencil.h>
#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

void pencil_hog_static( int NUMBER_OF_CELLS
                      , int NUMBER_OF_BINS
                      , bool GAUSSIAN_WEIGHTS
                      , bool SPARTIAL_WEIGHTS
                      , bool SIGNED_HOG
                      , const int rows
                      , const int cols
                      , const int step
                      , const uint8_t image[]
                      , const int num_locations
                      , const float location[][2]
                      , const float block_size
                      , float hist[]    //out
                      );

void pencil_hog_dynamic( int NUMBER_OF_CELLS
                       , int NUMBER_OF_BINS
                       , bool GAUSSIAN_WEIGHTS
                       , bool SPARTIAL_WEIGHTS
                       , bool SIGNED_HOG
                       , const int rows
                       , const int cols
                       , const int step
                       , const uint8_t image[]
                       , const int num_locations
                       , const float location[][2]
                       , const float block_size[][2]
                       , float hist[]    //out
                       );

#ifdef __cplusplus
} // extern "C"
#endif

#endif //HOG_PENCIL_H
