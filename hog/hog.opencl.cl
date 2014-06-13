//parameters
#define NUMBER_OF_CELLS 1
#define NUMBER_OF_BINS 8
#define GAUSSIAN_WEIGHTS 1
#define SPARTIAL_WEIGHTS 0
#define SIGNED_HOG 1
//parameters end

#if SIGNED_HOG
#define BINSIZE_IN_DEGREES (360.0f / NUMBER_OF_BINS)
#else
#define BINSIZE_IN_DEGREES (180.0f / NUMBER_OF_BINS)
#endif

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

//get_global_id(0): per location
//get_global_id(1): x coordinate
//get_global_id(2): y coordinate

void atomic_add_float(volatile __global float* source, float operand) {
    volatile __global int* sourceAsInt = (volatile __global int*)source;
    float oldVal;
    float newVal;
    do {
        oldVal = *source;
        newVal = oldVal + operand;
    } while (atomic_cmpxchg(sourceAsInt, as_int(oldVal), as_int(newVal)) != as_int(oldVal));
}

#define HIST_INDEX(location, celly, cellx, bin)                         \
        ( hist                                                          \
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
                            , const int num_locations, __global const float *location_x, __global const float *location_y
                            , const float blck_size
                            , const int num_pixels_in_block //=ceil(blck_size)
                            , __global float *hist
                            )
{
    int i = get_global_id(0);
    if (i >= num_locations)
        return;
    
    float minx = location_x[i] - blck_size / 2.0f;
    float miny = location_y[i] - blck_size / 2.0f;
    float maxx = location_x[i] + blck_size / 2.0f;
    float maxy = location_y[i] + blck_size / 2.0f;

    int minxi = max((int)ceil(minx), 1);
    int minyi = max((int)ceil(miny), 1);
    int maxxi = min((int)floor(maxx), cols - 2);
    int maxyi = min((int)floor(maxy), rows - 2);

    int pointy = minyi + get_global_id(2); //int+size_t !!
    int pointx = minxi + get_global_id(1);

    if (pointy > maxyi || pointx > maxxi)
        return;

    //Read the image
    float mdx = image[(pointy  )*step + pointx+1] - image[(pointy  )*step + pointx-1];
    float mdy = image[(pointy+1)*step + pointx  ] - image[(pointy-1)*step + pointx  ];

    //calculate the magnitude
    float magnitude = hypot(mdx, mdy);

    //calculate the orientation
#if SIGNED_HOG
    float orientation = atan2(mdy, mdx) / (float)M_PI * 180.0f;
#else
    float orientation = tan2(mdy / mdx + DBL_EPSILON) / (float)M_PI * 180.0f + 90.0f;
#endif

#if GAUSSIAN_WEIGHTS
    {
        float sigma = blck_size / 2.0f;
        float sigmaSq = sigma*sigma;
        float m1p2sigmaSq = -1.0f / (2.0f * sigmaSq);
        float dy = pointy - location_y[i];
        float dx = pointx - location_x[i];
        float dySq = dy*dy;
        float dxSq = dx*dx;
        magnitude *= exp((dxSq+dySq) * m1p2sigmaSq);
    }
#endif
    float relative_orientation = (orientation - BINSIZE_IN_DEGREES/2.0f) / BINSIZE_IN_DEGREES;
    int bin1 = ceil(relative_orientation);
    int bin0 = bin1 - 1;
    float bin_weight0 = magnitude * (bin1 - relative_orientation);
    float bin_weight1 = magnitude * (relative_orientation - bin0);
    bin0 = (bin0 + NUMBER_OF_BINS) % NUMBER_OF_BINS;
    bin1 = (bin1 + NUMBER_OF_BINS) % NUMBER_OF_BINS;
#if SPARTIAL_WEIGHTS
    if (cellyi >= 0 && cellxi >= 0) {
        atomic_add_float( HIST_INDEX(i, cellyi  , cellxi  , bin0), yscale0 * xscale0 * bin_weight0);
        atomic_add_float( HIST_INDEX(i, cellyi  , cellxi  , bin1), yscale0 * xscale0 * bin_weight1);
    }
    if (cellyi >= 0 && cellxi < NUMBER_OF_CELLS - 1) {
        atomic_add_float( HIST_INDEX(i, cellyi  , cellxi+1, bin0), yscale0 * xscale1 * bin_weight0);
        atomic_add_float( HIST_INDEX(i, cellyi  , cellxi+1, bin1), yscale0 * xscale1 * bin_weight1);
    }
    if (cellyi < NUMBER_OF_CELLS - 1 && cellxi >= 0) {
        atomic_add_float( HIST_INDEX(i, cellyi+1, cellxi  , bin0), yscale1 * xscale0 * bin_weight0);
        atomic_add_float( HIST_INDEX(i, cellyi+1, cellxi  , bin1), yscale1 * xscale0 * bin_weight1);
    }
    if (cellyi < NUMBER_OF_CELLS - 1 && cellxi < NUMBER_OF_CELLS - 1) {
        atomic_add_float( HIST_INDEX(i, cellyi+1, cellxi+1, bin0), yscale1 * xscale1 * bin_weight0);
        atomic_add_float( HIST_INDEX(i, cellyi+1, cellxi+1, bin1), yscale1 * xscale1 * bin_weight1);
    }
#elif NUMBER_OF_CELLS == 1
    atomic_add_float( HIST_INDEX(i, 0, 0, bin0), bin_weight0);
    atomic_add_float( HIST_INDEX(i, 0, 0, bin1), bin_weight1);
#else
    int cellxi = floor((pointx - minx) / cell_size);
    int cellyi = floor((pointy - miny) / cell_size);

    atomic_add_float( HIST_INDEX(i, cellyi, cellxi, bin0), bin_weight0);
    atomic_add_float( HIST_INDEX(i, cellyi, cellxi, bin1), bin_weight1);
#endif
}
