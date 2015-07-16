#ifndef CVT_COLOR_PENCIL__H__
#define CVT_COLOR_PENCIL__H__

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif
    void pencil_RGB2Gray( const int rows
                        , const int cols
                        , const int src_step
                        , const int dst_step
                        , const uint8_t src[]
                        , uint8_t dst[]
                        );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* CVT_COLOR_PENCIL__H__ */
