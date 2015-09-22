#ifndef NUMBER_OF_CELLS
#error NUMBER_OF_CELLS not defined
#endif
#ifndef NUMBER_OF_BINS
#error NUMBER_OF_BINS not defined
#endif
#ifndef GAUSSIAN_WEIGHTS
#error GAUSSIAN_WEIGHTS not defined
#endif
#ifndef SPARTIAL_WEIGHTS
#error SPARTIAL_WEIGHTS not defined
#endif
#ifndef SIGNED_HOG
#error SIGNED_HOG not defined
#endif
#ifndef STATIC_BLOCK_SIZE
#error STATIC_BLOCK_SIZE not defined
#endif
#if (USE_GLOBAL + USE_LOCAL + USE_PRIVATE != 1) 
#error define one of one of USE_PRIVATE, USE_LOCAL, USE_GLOBAL to 1 to store intermediate histogram data in private, local, or global memory
#endif

#if VECTOR_SIZE == 16
    #define VLOADN     vload16
    #define VSTOREN    vstore16
    #define FLOATN     float16
    #define WG_REDUCEN wg_reduce16
#elif VECTOR_SIZE == 8
    #define VLOADN     vload8
    #define VSTOREN    vstore8
    #define FLOATN     float8
    #define WG_REDUCEN wg_reduce8
#elif VECTOR_SIZE == 4
    #define VLOADN     vload4
    #define VSTOREN    vstore4
    #define FLOATN     float4
    #define WG_REDUCEN wg_reduce4
#elif VECTOR_SIZE == 2
    #define VLOADN     vload2
    #define VSTOREN    vstore2
    #define FLOATN     float2
    #define WG_REDUCEN wg_reduce2
#elif VECTOR_SIZE == 1
    #define VLOADN(x, y)     (y)[x]
    #define VSTOREN(x, y, z) (z)[y] = (x)
    #define FLOATN           float
    #define WG_REDUCEN       wg_reduce
#else
    #error Unsupported vector size!
#endif

#define TOTAL_NUMBER_OF_BINS (NUMBER_OF_CELLS * NUMBER_OF_CELLS * NUMBER_OF_BINS)
#define BIN_INDEX(cellx, celly, bin) ((cellx) * NUMBER_OF_CELLS * NUMBER_OF_BINS + (celly) * NUMBER_OF_BINS + bin)
#define ROUND_DOWN_TO_VECTOR_SIZE(VAL) (VAL & ~(VECTOR_SIZE-1))

//Find the largest 2^n where (2^n) < v
uint lower_power_of_two(uint v)
{
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v++;
    return v >> 1;
}

//Returns the sum of 'value' among all workitems within this workgroup
//Overwrites temp, so add barriers before calling this function when necessary
//Provides functionality similar to OpenCL 2.0's work_group_reduce_add
#define DEFINE_REDUCE(name, datatype)                                            \
datatype name(             datatype                  value                       \
             , local       datatype * const restrict temp                        \
             ,       const size_t                    thread_id                   \
             ,       const size_t                    array_size                  \
             )                                                                   \
{                                                                                \
    /* Store partial sum to temp */                                              \
    if (thread_id < array_size) {                                                \
        temp[thread_id] = value;                                                 \
    }                                                                            \
    barrier(CLK_LOCAL_MEM_FENCE);                                                \
                                                                                 \
    /* Do a reduction on temp */                                                 \
    for ( size_t offset = lower_power_of_two(array_size)                         \
        ; offset > 0                                                             \
        ; offset >>= 1                                                           \
        ) {                                                                      \
        bool check = (thread_id < offset) && (thread_id + offset < array_size ); \
        if (check) {                                                             \
            value += temp[thread_id + offset];                                   \
            temp[thread_id] = value;                                             \
        }                                                                        \
        barrier(CLK_LOCAL_MEM_FENCE);                                            \
    }                                                                            \
    return temp[0];                                                              \
}

DEFINE_REDUCE(wg_reduce  , float  )
DEFINE_REDUCE(wg_reduce2 , float2 )
DEFINE_REDUCE(wg_reduce3 , float3 )
DEFINE_REDUCE(wg_reduce4 , float4 )
DEFINE_REDUCE(wg_reduce8 , float8 )
DEFINE_REDUCE(wg_reduce16, float16)
#undef DEFINE_REDUCE

#if USE_PRIVATE
void atomic_add_float(volatile private float* source, float operand) {
    *source += operand;
}
#else
#if USE_GLOBAL
    #define LOCAL global
#elif USE_LOCAL
    #define LOCAL local
#endif
void atomic_add_float(volatile LOCAL float* source, float operand) {
    volatile LOCAL int* sourceAsInt = (volatile LOCAL int*)source;
    int oldVal;
    int oldVal_2 = 0;
    do {
        oldVal = oldVal_2;
        float newVal = as_float(oldVal) + operand;
        oldVal_2 = atomic_cmpxchg(sourceAsInt, oldVal, as_int(newVal));
    } while (oldVal_2 != oldVal);
}
#undef LOCAL
#endif

kernel void FUNCTIONNAME(        const          int2                    img_size
                        ,        const unsigned int                     step_
                        , global const unsigned char   * const restrict image
                        ,        const unsigned int                     num_locations
                        , global const          float2 * const restrict location_global     //num_locations elements
#if STATIC_BLOCK_SIZE
                        ,        const          float                   blck_size           //Note: static block size: float
#else
                        , global const          float2 * const restrict blocksize_global    //Note: dynamic block size: float2 array, num_locations elements
#endif
                        , global                float  *       restrict hist_global         //num_locations * NUMBER_OF_CELLS * NUMBER_OF_CELLS * NUMBER_OF_BINS elements
#if defined(USE_LOCAL)
                        ,  local                float  *       restrict hist_local          //get_local_size(2) * NUMBER_OF_CELLS * NUMBER_OF_CELLS * NUMBER_OF_BINS elements
#elif defined(USE_PRIVATE)
                        ,  local                float  *       restrict temp                //get_local_size(2) * get_local_size(0) * get_local_size(1) * VECTOR_SIZE elements
#endif
                        )
{
    //get_global_id(0): x coordinate
    //get_global_id(1): y coordinate
    //get_global_id(2): per location
    
    const size_t location_global_idx = get_global_id(2);
    hist_global += TOTAL_NUMBER_OF_BINS * location_global_idx;
#if USE_PRIVATE
    float hist_private[TOTAL_NUMBER_OF_BINS] = { 0 };
    #define hist_temp hist_private
#else
#if USE_GLOBAL
    #define hist_temp hist_global
#elif USE_LOCAL
    hist_local += TOTAL_NUMBER_OF_BINS * get_local_id(2);
    #define hist_temp hist_local
#endif

    //Clear hist_temp
    if (location_global_idx < num_locations) {
        const size_t threadId = get_local_id(0) + get_local_size(0) * get_local_id(1);
        const size_t groupSize = get_local_size(0) * get_local_size(1);
        size_t i = threadId * VECTOR_SIZE;
        for(; i < ROUND_DOWN_TO_VECTOR_SIZE(TOTAL_NUMBER_OF_BINS); i += groupSize * VECTOR_SIZE) {
            VSTOREN(0.0f, 0, &hist_temp[i]);
        }
#if (TOTAL_NUMBER_OF_BINS % VECTOR_SIZE > 0)
        for(; i < TOTAL_NUMBER_OF_BINS; i += groupSize) {
            hist_temp[i] = 0.0f;
        }
#endif
    }
#if USE_GLOBAL
    barrier(CLK_GLOBAL_MEM_FENCE);
#elif USE_LOCAL
    barrier(CLK_LOCAL_MEM_FENCE);
#endif
#endif
    
    if (location_global_idx < num_locations) {
        float2 location = location_global[location_global_idx];
#if !STATIC_BLOCK_SIZE
        float2 blck_size = blocksize_global[location_global_idx];
#endif
        float2 minloc = location - blck_size * 0.5f;
        float2 maxloc = location + blck_size * 0.5f;

        //Not converting to uint before min/max, because without saturation it can underflow, and saturation is expensive (if we could saturate to 1 to img_size-2..)
        uint2 mini = convert_uint2(clamp(convert_int2_rtp(minloc), 1, img_size - 2));
        uint2 maxi = convert_uint2(clamp(convert_int2_rtn(maxloc), 1, img_size - 2));

        for (uint pointy = mini.y + (uint)get_local_id(1); pointy <= maxi.y; pointy += (uint)get_local_size(1)) {
            for (uint pointx = mini.x + (uint)get_local_id(0); pointx <= maxi.x; pointx += (uint)get_local_size(0)) {
                float2 point_f = (float2)(pointx, pointy);
                //Read the image

                uint point_idx = mad24(pointy, step_, pointx);
                uchar3 row_pixels = vload3(0, &image[point_idx-1]);
                uchar  top_pixel  = image[point_idx+step_];
                uchar  bot_pixel  = image[point_idx-step_];
                
                if (top_pixel == bot_pixel && row_pixels.s2 == row_pixels.s0) continue; //If the pixels have the same value, no contribution. Both faster, and correct ( atan2pi(0,0) is undefined ).

                float mdx = row_pixels.s2 - row_pixels.s0;
                float mdy = top_pixel     - bot_pixel;

                //calculate the magnitude
                float magnitude = hypot(mdx, mdy);

                //calculate the orientation
#if SIGNED_HOG
                float orientation = atan2pi(mdy, mdx) * 0.5f;
#else
                float orientation = atan2pi(mdy, mdx) + 0.5f;
#endif

#if GAUSSIAN_WEIGHTS
                {
                    float2 m1p2sigmaSq = -2.0f / (blck_size * blck_size);
                    float2 distance = point_f - location;
                    float2 distanceSq = distance*distance;
                    magnitude *= exp(dot(distanceSq, m1p2sigmaSq));
                }
#endif
                float relative_orientation = orientation * NUMBER_OF_BINS - 0.5f;
            
                int   bin0        = convert_int_rtn(relative_orientation);
                float bin_weight1 = relative_orientation - convert_float(bin0);
            
                int   bin1        = 1    + bin0;
                float bin_weight0 = 1.0f - bin_weight1;
            
                int2   bin        = ((int2)(bin0, bin1) + NUMBER_OF_BINS) % NUMBER_OF_BINS;
                float2 bin_weight = magnitude * (float2)(bin_weight0, bin_weight1);

#if NUMBER_OF_CELLS == 1
                atomic_add_float( &(hist_temp[bin.s0]), bin_weight.s0);
                atomic_add_float( &(hist_temp[bin.s1]), bin_weight.s1);
#elif SPARTIAL_WEIGHTS
                float2 relative_pos = (point_f - minloc) * NUMBER_OF_CELLS / blck_size - 0.5f;
                int2 celli = convert_int2_rtn(relative_pos);

                const float2 scale1 = relative_pos - convert_float2(celli);
                const float2 scale0 = 1.0f - scale1;
                const float4 scale  = (float4)(scale0, scale1);
                const float4 nearly_scaled_bin_weight = scale.s0022 * bin_weight.s0101;

                if (celli.y >= 0) {
                    const float4 scaled_bin_weight = scale.s1 * nearly_scaled_bin_weight;
                    if (celli.x >= 0) {
                        atomic_add_float( &(hist_temp[BIN_INDEX(celli.y  , celli.x  , bin.s0)]), scaled_bin_weight.s0);
                        atomic_add_float( &(hist_temp[BIN_INDEX(celli.y  , celli.x  , bin.s1)]), scaled_bin_weight.s1);
                    }
                    if (celli.x < NUMBER_OF_CELLS - 1) {
                        atomic_add_float( &(hist_temp[BIN_INDEX(celli.y  , celli.x+1, bin.s0)]), scaled_bin_weight.s2);
                        atomic_add_float( &(hist_temp[BIN_INDEX(celli.y  , celli.x+1, bin.s1)]), scaled_bin_weight.s3);
                    }
                }
                if (celli.y < NUMBER_OF_CELLS - 1) {
                    const float4 scaled_bin_weight = scale.s3 * nearly_scaled_bin_weight;
                    if (celli.x >= 0) {
                        atomic_add_float( &(hist_temp[BIN_INDEX(celli.y+1, celli.x  , bin.s0)]), scaled_bin_weight.s0);
                        atomic_add_float( &(hist_temp[BIN_INDEX(celli.y+1, celli.x  , bin.s1)]), scaled_bin_weight.s1);
                    }
                    if (celli.x < NUMBER_OF_CELLS - 1) {
                        atomic_add_float( &(hist_temp[BIN_INDEX(celli.y+1, celli.x+1, bin.s0)]), scaled_bin_weight.s2);
                        atomic_add_float( &(hist_temp[BIN_INDEX(celli.y+1, celli.x+1, bin.s1)]), scaled_bin_weight.s3);
                    }
                }
#else
                int2 celli = convert_int2_rtn((point_f - minloc) * NUMBER_OF_CELLS / blck_size);

                atomic_add_float( &(hist_temp[BIN_INDEX(celli.y, celli.x, bin.s0)]), bin_weight.s0);
                atomic_add_float( &(hist_temp[BIN_INDEX(celli.y, celli.x, bin.s1)]), bin_weight.s1);
#endif
            }
        }
    }
    #undef hist_temp
#if USE_LOCAL
    barrier(CLK_LOCAL_MEM_FENCE);
    
    if (location_global_idx < num_locations) {
        //Add local version to global
        const size_t threadId = get_local_id(0) + get_local_size(0) * get_local_id(1);
        const size_t groupSize = get_local_size(0) * get_local_size(1);
        size_t i = threadId * VECTOR_SIZE;
        for(; i < ROUND_DOWN_TO_VECTOR_SIZE(TOTAL_NUMBER_OF_BINS); i += groupSize * VECTOR_SIZE) {
            VSTOREN(VLOADN(0, &hist_local[i]), 0, &hist_global[i]);
        }
#if (TOTAL_NUMBER_OF_BINS % VECTOR_SIZE > 0)
        for(; i < TOTAL_NUMBER_OF_BINS; i += groupSize) {
            hist_global[i] = hist_local[i];
        }
#endif
    }
#elif defined(USE_PRIVATE)
    //Use lock-step reduction and store result to hist_global
    const size_t threadId = get_local_id(0) + get_local_size(0) * get_local_id(1);
    const size_t groupSize = get_local_size(0) * get_local_size(1);
    temp += get_local_id(2) * groupSize * VECTOR_SIZE;

    int i = 0;
    for(; i < ROUND_DOWN_TO_VECTOR_SIZE(TOTAL_NUMBER_OF_BINS); i += VECTOR_SIZE) {
        barrier(CLK_LOCAL_MEM_FENCE);
        FLOATN value = WG_REDUCEN(VLOADN(0, &hist_private[i]), (local FLOATN *)temp, threadId, groupSize);
        if (0 == threadId && location_global_idx < num_locations) {
            VSTOREN(value, 0, &hist_global[i]);
        }
    }
#if (TOTAL_NUMBER_OF_BINS % VECTOR_SIZE > 0)
    for(; i < TOTAL_NUMBER_OF_BINS; ++i) {
        barrier(CLK_LOCAL_MEM_FENCE);
        float value = wg_reduce(hist_private[i], temp, threadId, groupSize);
        if (0 == threadId && location_global_idx < num_locations) {
            hist_global[i] = value;
        }
    }
#endif
#endif
}
