//parameters
#define NUMBER_OF_CELLS 1
#define NUMBER_OF_BINS 8
#define GAUSSIAN_WEIGHTS 1
#define SPARTIAL_WEIGHTS 0
#define SIGNED_HOG 1
//parameters end

#define TOTAL_NUMBER_OF_BINS (NUMBER_OF_CELLS * NUMBER_OF_CELLS * NUMBER_OF_BINS)

//get_global_id(0): per location
//get_global_id(1): x coordinate
//get_global_id(2): y coordinate

inline void atomic_add_float_global(volatile __global float* source, float operand) {
    volatile __global int* sourceAsInt = (volatile __global int*)source;
    float oldVal;
    float newVal;
    do {
        oldVal = *source;
        newVal = oldVal + operand;
    } while (atomic_cmpxchg(sourceAsInt, as_int(oldVal), as_int(newVal)) != as_int(oldVal));
}

inline void atomic_add_float_local(volatile __local float* source, float operand) {
    volatile __local int* sourceAsInt = (volatile __local int*)source;
    float oldVal;
    float newVal;
    do {
        oldVal = *source;
        newVal = oldVal + operand;
    } while (atomic_cmpxchg(sourceAsInt, as_int(oldVal), as_int(newVal)) != as_int(oldVal));
}

inline size_t get_linear_local_id() {
    return get_local_id(2) * get_local_size(0) * get_local_size(1)
         + get_local_id(1) * get_local_size(0)
         + get_local_id(0);
}

inline size_t get_linear_local_size() {
    return get_local_size(0) * get_local_size(1) * get_local_size(2);
}

#define HIST_INDEX(location, celly, cellx, bin)                         \
        ( hist_local                                                    \
        + (location)*NUMBER_OF_CELLS*NUMBER_OF_CELLS*NUMBER_OF_BINS     \
        + (celly)*NUMBER_OF_CELLS*NUMBER_OF_BINS                        \
        + (cellx)*NUMBER_OF_BINS                                        \
        + (bin)                                                         \
        )

// OpenCL 1.1 does not have clEnqueueFillBuffer...
__kernel void fill_zeros(__global float* arr, int size) {
    int i = get_global_id(0);
    if (i < size) {
        arr[i] = 0.0f;
    }
}

//Calculates the intermediate data
__kernel void calc_histogram( const int rows, const int cols, const int step, __global const unsigned char *image
                            , const int num_locations, __global const float2 *location_global
                            , const float blck_size
                            , __global float *hist_global, __local float *hist_local
                            )
{
    //Clear hist_local
    size_t local_id = get_linear_local_id();
    size_t local_size = get_linear_local_size();
    for(size_t n = local_id; n < TOTAL_NUMBER_OF_BINS*num_locations; n += local_size) {
        *(hist_local + n) = 0.0f;
    }
    
    barrier(CLK_LOCAL_MEM_FENCE);

    int i = get_global_id(0);
    if (i < num_locations) {
        float2 location = location_global[i];
        int2 img_size = (int2)(cols, rows);

        float2 minloc = location - blck_size / 2.0f;
        float2 maxloc = location + blck_size / 2.0f;

        int2 mini = max(convert_int2_rtp(minloc), 1);
        int2 maxi = min(convert_int2_rtn(maxloc), img_size - 2);

        int2 point = mini + (int2)(get_global_id(1),get_global_id(2));

        if (all(point <= maxi)) {
            //Read the image
            float mdx = image[(point.y  )*step + point.x+1] - image[(point.y  )*step + point.x-1];
            float mdy = image[(point.y+1)*step + point.x  ] - image[(point.y-1)*step + point.x  ];

            //calculate the magnitude
            float magnitude = hypot(mdx, mdy);

            //calculate the orientation
#if SIGNED_HOG
            float orientation = atan2pi(mdy, mdx) / 2.0f;
#else
            float orientation = tan2pi(mdy / mdx + DBL_EPSILON) + 0.5f;
#endif

#if GAUSSIAN_WEIGHTS
            {
                float sigma = blck_size / 2.0f;
                float sigmaSq = sigma*sigma;
                float m1p2sigmaSq = -1.0f / (2.0f * sigmaSq);
                float2 distanceSq = pown(convert_float2(point) - location, 2);
                magnitude *= exp(dot(distanceSq, m1p2sigmaSq));
            }
#endif
            float relative_orientation = orientation * NUMBER_OF_BINS - 0.5f;
            int bin1 = ceil(relative_orientation);
            int bin0 = bin1 - 1;
            float bin_weight0 = magnitude * (bin1 - relative_orientation);
            float bin_weight1 = magnitude * (relative_orientation - bin0);
            bin0 = (bin0 + NUMBER_OF_BINS) % NUMBER_OF_BINS;
            bin1 = (bin1 + NUMBER_OF_BINS) % NUMBER_OF_BINS;

#if SPARTIAL_WEIGHTS
            if (cellyi >= 0 && cellxi >= 0) {
                atomic_add_float_local( HIST_INDEX(i, cellyi  , cellxi  , bin0), yscale0 * xscale0 * bin_weight0);
                atomic_add_float_local( HIST_INDEX(i, cellyi  , cellxi  , bin1), yscale0 * xscale0 * bin_weight1);
            }
            if (cellyi >= 0 && cellxi < NUMBER_OF_CELLS - 1) {
                atomic_add_float_local( HIST_INDEX(i, cellyi  , cellxi+1, bin0), yscale0 * xscale1 * bin_weight0);
                atomic_add_float_local( HIST_INDEX(i, cellyi  , cellxi+1, bin1), yscale0 * xscale1 * bin_weight1);
            }
            if (cellyi < NUMBER_OF_CELLS - 1 && cellxi >= 0) {
                atomic_add_float_local( HIST_INDEX(i, cellyi+1, cellxi  , bin0), yscale1 * xscale0 * bin_weight0);
                atomic_add_float_local( HIST_INDEX(i, cellyi+1, cellxi  , bin1), yscale1 * xscale0 * bin_weight1);
            }
            if (cellyi < NUMBER_OF_CELLS - 1 && cellxi < NUMBER_OF_CELLS - 1) {
                atomic_add_float_local( HIST_INDEX(i, cellyi+1, cellxi+1, bin0), yscale1 * xscale1 * bin_weight0);
                atomic_add_float_local( HIST_INDEX(i, cellyi+1, cellxi+1, bin1), yscale1 * xscale1 * bin_weight1);
            }
#elif NUMBER_OF_CELLS == 1
            atomic_add_float_local( HIST_INDEX(i, 0, 0, bin0), bin_weight0);
            atomic_add_float_local( HIST_INDEX(i, 0, 0, bin1), bin_weight1);
#else
            int2 celli = convert_int2_rtn((convert_float2(point) - minloc) / cell_size);

            atomic_add_float_local( HIST_INDEX(i, celli.y, celli.x, bin0), bin_weight0);
            atomic_add_float_local( HIST_INDEX(i, celli.y, celli.x, bin1), bin_weight1);
#endif
        }
    }    
    barrier(CLK_LOCAL_MEM_FENCE);
    
    //Add local version to global
    for(size_t n = local_id; n < TOTAL_NUMBER_OF_BINS*num_locations; n += local_size) {
        atomic_add_float_global( hist_global + n, *(hist_local + n) );
    }
}
