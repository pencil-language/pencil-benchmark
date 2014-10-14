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

// OpenCL 1.1 does not have clEnqueueFillBuffer...
__kernel void fill_zeros(__global float * restrict arr, int size) {
    int i = get_global_id(0);
    if (i < size) {
        arr[i] = 0.0f;
    }
}

//Calculates the intermediate data
__kernel void calc_histogram( const int rows, const int cols, const int step_, __global const unsigned char * restrict image
                            , const int num_locations, __global const float2 * restrict location_global, __global const float2 * restrict blocksize_global
                            , __global float hist_global[][NUMBER_OF_CELLS][NUMBER_OF_CELLS][NUMBER_OF_BINS]
                            ,  __local float  hist_local[][NUMBER_OF_CELLS][NUMBER_OF_CELLS][NUMBER_OF_BINS]
                            )
{
    //Clear hist_local
    for(int a = get_local_id(2); a < num_locations; a+=get_local_size(2)) {
        for(int b = get_local_id(1); b < NUMBER_OF_CELLS; b+=get_local_size(1)) {
            for(int c = get_local_id(0); c < NUMBER_OF_CELLS; c+=get_local_size(0)) {
                for(int n = 0; n < NUMBER_OF_BINS ; ++n) {
                    hist_local[a][b][c][n] = 0.0f;
                }
            }
        }
    }
    
    barrier(CLK_LOCAL_MEM_FENCE);

    int locationIdx = get_global_id(2);
    if (locationIdx < num_locations) {
        float2 location = location_global[locationIdx];
        float2 blck_size = blocksize_global[locationIdx];
        int2 img_size = (int2)(cols, rows);

        float2 minloc = location - blck_size / 2.0f;
        float2 maxloc = location + blck_size / 2.0f;

        int2 mini = max(convert_int2_rtp(minloc), 1);
        int2 maxi = min(convert_int2_rtn(maxloc), img_size - 2);

        int2 point = mini + (int2)(get_global_id(0),get_global_id(1));

        if (all(point <= maxi)) {
            //Read the image
            float mdx = image[(point.y  )*step_ + point.x+1] - image[(point.y  )*step_ + point.x-1];
            float mdy = image[(point.y+1)*step_ + point.x  ] - image[(point.y-1)*step_ + point.x  ];

            //calculate the magnitude
            float magnitude = hypot(mdx, mdy);

            //calculate the orientation
#if SIGNED_HOG
            float orientation = atan2pi(mdy, mdx) / 2.0f;
#else
            float orientation = atan2pi(mdy, mdx) + 0.5f;
#endif

#if GAUSSIAN_WEIGHTS
            {
                float2 sigma = blck_size / 2.0f;
                float2 sigmaSq = sigma*sigma;
                float2 m1p2sigmaSq = -1.0f / (2.0f * sigmaSq);
                float2 distanceSq = pown(convert_float2(point) - location, 2);
                magnitude *= exp(dot(distanceSq, m1p2sigmaSq));
            }
#endif
            float relative_orientation = orientation * NUMBER_OF_BINS - 0.5f;
            int bin1 = convert_int_rtp(relative_orientation);
            int bin0 = bin1 - 1;
            float bin_weight0 = magnitude * (convert_float(bin1) - relative_orientation);
            float bin_weight1 = magnitude * (relative_orientation - convert_float(bin0));
            bin0 = (bin0 + NUMBER_OF_BINS) % NUMBER_OF_BINS;
            bin1 = (bin1 + NUMBER_OF_BINS) % NUMBER_OF_BINS;

#if NUMBER_OF_CELLS == 1
            atomic_add_float_local( &(hist_local[locationIdx][0][0][bin0]), bin_weight0);
            atomic_add_float_local( &(hist_local[locationIdx][0][0][bin1]), bin_weight1);
#elif SPARTIAL_WEIGHTS
            float2 relative_pos = (convert_float2(point) - minloc) * NUMBER_OF_CELLS / blck_size - 0.5f;
            int2 celli = convert_int2_rtn(relative_pos);

            float2 scale1 = relative_pos - convert_float2(celli);
            float2 scale0 = 1.0f - scale1;

            if (celli.y >= 0 && celli.x >= 0) {
                atomic_add_float_local( &(hist_local[locationIdx][celli.y  ][celli.x  ][bin0]), scale0.y * scale0.x * bin_weight0);
                atomic_add_float_local( &(hist_local[locationIdx][celli.y  ][celli.x  ][bin1]), scale0.y * scale0.x * bin_weight1);
            }
            if (celli.y >= 0 && celli.x < NUMBER_OF_CELLS - 1) {
                atomic_add_float_local( &(hist_local[locationIdx][celli.y  ][celli.x+1][bin0]), scale0.y * scale1.x * bin_weight0);
                atomic_add_float_local( &(hist_local[locationIdx][celli.y  ][celli.x+1][bin1]), scale0.y * scale1.x * bin_weight1);
            }
            if (celli.y < NUMBER_OF_CELLS - 1 && celli.x >= 0) {
                atomic_add_float_local( &(hist_local[locationIdx][celli.y+1][celli.x  ][bin0]), scale1.y * scale0.x * bin_weight0);
                atomic_add_float_local( &(hist_local[locationIdx][celli.y+1][celli.x  ][bin1]), scale1.y * scale0.x * bin_weight1);
            }
            if (celli.y < NUMBER_OF_CELLS - 1 && celli.x < NUMBER_OF_CELLS - 1) {
                atomic_add_float_local( &(hist_local[locationIdx][celli.y+1][celli.x+1][bin0]), scale1.y * scale1.x * bin_weight0);
                atomic_add_float_local( &(hist_local[locationIdx][celli.y+1][celli.x+1][bin1]), scale1.y * scale1.x * bin_weight1);
            }
#else
            int2 celli = convert_int2_rtn((convert_float2(point) - minloc) * NUMBER_OF_CELLS / blck_size);

            atomic_add_float_local( &(hist_local[locationIdx][celli.y][celli.x][bin0]), bin_weight0);
            atomic_add_float_local( &(hist_local[locationIdx][celli.y][celli.x][bin1]), bin_weight1);
#endif
        }
    }    
    barrier(CLK_LOCAL_MEM_FENCE);
    
    //Add local version to global
    for(int a = get_local_id(2); a < num_locations; a+=get_local_size(2)) {
        for(int b = get_local_id(1); b < NUMBER_OF_CELLS; b+=get_local_size(1)) {
            for(int c = get_local_id(0); c < NUMBER_OF_CELLS; c+=get_local_size(0)) {
                for(int n = 0; n < NUMBER_OF_BINS ; ++n) {
                    atomic_add_float_global( &(hist_global[a][b][c][n]), hist_local[a][b][c][n] );
                }
            }
        }
    }
}
