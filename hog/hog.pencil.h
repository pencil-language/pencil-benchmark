#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

//parameters
#define NUMBER_OF_CELLS 4
#define NUMBER_OF_BINS 8
#define GAUSSIAN_WEIGHTS 1
#define SPARTIAL_WEIGHTS 1
#define SIGNED_HOG 1
//parameters end

static const int HISTOGRAM_BINS = NUMBER_OF_CELLS * NUMBER_OF_CELLS * NUMBER_OF_BINS;

void pencil_hog( const int rows
               , const int cols
               , const int step
               , const uint8_t image[]
               , const int num_locations
               , const double location_x[]
               , const double location_y[]
               , const double block_size
               , double hist[]    //out
               );

#ifdef __cplusplus
} // extern "C"
#endif
