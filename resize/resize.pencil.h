#ifndef RESIZE_PENCIL_H
#define RESIZE_PENCIL_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

void pencil_resize_LN( const int original_rows,  const int original_cols,  const int original_step,  const uint8_t original[]
                     , const int resampled_rows, const int resampled_cols, const int resampled_step,       uint8_t resampled[]
                     );

#ifdef __cplusplus
} // extern "C"
#endif

#endif