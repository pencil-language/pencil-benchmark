#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

//parameters
#define NUMBER_OF_CELLS 1
#define NUMBER_OF_BINS 4
#define GAUSSIAN_WEIGHTS 0
#define SPARTIAL_WEIGHTS 0
#define SIGNED_HOG 1
//parameters end

static const int HISTOGRAM_BINS = NUMBER_OF_CELLS * NUMBER_OF_CELLS * NUMBER_OF_BINS;

void pencil_hog( const int rows
               , const int cols
               , const int step
               , const uint8_t image[]
               , const int num_locations
               , const float location_x[]
               , const float location_y[]
               , const float block_size
               , float hist[]    //out
               );

#ifdef __cplusplus
} // extern "C"
#endif
