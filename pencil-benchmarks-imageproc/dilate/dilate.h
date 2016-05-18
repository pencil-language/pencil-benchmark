#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

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
                  , const int anchor_x
                  , const int anchor_y
                  );

#ifdef __cplusplus
} // extern "C"
#endif
