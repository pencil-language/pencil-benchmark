#include "hog.pencil.h"

#include <math.h>   //sqrt, atan, atan2, ceil, floor
#include <string.h> //memset

#if SIGNED_HOG
static const float BINSIZE_IN_DEGREES = 360.0f / NUMBER_OF_BINS;
#else
static const float BINSIZE_IN_DEGREES = 180.0f / NUMBER_OF_BINS;
#endif

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))

static void normalize(float hist[NUMBER_OF_CELLS][NUMBER_OF_CELLS][NUMBER_OF_BINS]) {
    static const float l2_hys_threshold = 0.15f;
    float sum = 0.0f;
    for (int i = 0; i < NUMBER_OF_CELLS; ++i)
        for (int j = 0; j < NUMBER_OF_CELLS; ++j)
            for (int k = 0; k < NUMBER_OF_BINS; ++k)
                sum += hist[i][j][k] * hist[i][j][k];
    if (sum == 0.0f)
        return;
    float scale = 1.0f/sqrt(sum);
    sum = 0.0f;
    for (int i = 0; i < NUMBER_OF_CELLS; ++i)
        for (int j = 0; j < NUMBER_OF_CELLS; ++j)
            for (int k = 0; k < NUMBER_OF_BINS; ++k) {
                hist[i][j][k] = min(hist[i][j][k] * scale, l2_hys_threshold);
                sum += hist[i][j][k] * hist[i][j][k];
            }
    for (int i = 0; i < NUMBER_OF_CELLS; ++i)
        for (int j = 0; j < NUMBER_OF_CELLS; ++j)
            for (int k = 0; k < NUMBER_OF_BINS; ++k)
                hist[i][j][k] = hist[i][j][k] * scale;
}

static void hog_multi( const int rows
                     , const int cols
                     , const int step
                     , const uint8_t image[static const restrict rows][step]
                     , const int num_locations
                     , const float location_x[static const restrict num_locations]
                     , const float location_y[static const restrict num_locations]
                     , const float block_size
                     , float hist[static const restrict num_locations][NUMBER_OF_CELLS][NUMBER_OF_CELLS][NUMBER_OF_BINS]    //out
                     ) {
#pragma scop

#pragma pencil independent
    for (int i = 0; i < num_locations; ++i) {
	    const float cell_size = block_size / NUMBER_OF_CELLS;
	    const float minx = location_x[i] - block_size / 2.0f;
	    const float miny = location_y[i] - block_size / 2.0f;
	    const float maxx = location_x[i] + block_size / 2.0f;
	    const float maxy = location_y[i] + block_size / 2.0f;

	    const int minxi = max(ceil(minx), 1);
	    const int minyi = max(ceil(miny), 1);
	    const int maxxi = min(ceil(maxx), cols - 2);
	    const int maxyi = min(ceil(maxy), rows - 2);

	    memset(hist[i], 0, HISTOGRAM_BINS);

#if GAUSSIAN_WEIGHTS
	    const float sigma = block_size / 2.0f;
	    const float sigmaSq = sigma*sigma;
	    const float m1p2sigmaSq = -1.0f / (2.0f * sigmaSq);
#endif

	    #pragma pencil independent reduction(+:hist[i])
	    for (int pointy = minyi; pointy <= maxyi; ++pointy) {
#if SPARTIAL_WEIGHTS
            float relative_pos_y = (pointy - miny) / cell_size - 0.5f;
            int cellyi = floor(relative_pos_y);
            float yscale1 = relative_pos_y - cellyi;
            float yscale0 = 1.0f - yscale1;
#endif
#if GAUSSIAN_WEIGHTS
            float dy = pointy - location_y[i];
            float dySq = dy*dy;
#endif

            #pragma pencil independent reduction(+:hist[i])
            for (int pointx = minxi; pointx <= maxxi; ++pointx) {
#if SPARTIAL_WEIGHTS
                float relative_pos_x = (pointx - minx) / cell_size - 0.5f;
                cellxi = floor(relative_pos_x);
                float xscale1 = relative_pos_x - cellxi;
                float xscale0 = 1.0f - xscale1;
#endif

#if GAUSSIAN_WEIGHTS
                float dx = pointx - location_x[i];
                float dxSq = dx*dx;
#endif
                float mdx = image[pointy][pointx+1] - image[pointy][pointx-1];
                float mdy = image[pointy+1][pointx] - image[pointy-1][pointx];

                float magnitude = hypot(mdx, mdy);   //or = sqrt(mdx*mdx + mdy*mdy);
#if SIGNED_HOG
                float orientation = atan2(mdy, mdx) / M_PI * 180.0f;
#else
                float orientation = tan2(mdy / mdx + DBL_EPSILON) / M_PI * 180.0f + 90.0f;
#endif
#if GAUSSIAN_WEIGHTS
                magnitude *= exp((dxSq+dySq) * m1p2sigmaSq);
#endif
                float relative_orientation = (orientation - BINSIZE_IN_DEGREES/2.0) / BINSIZE_IN_DEGREES;
                int bin1 = ceil(relative_orientation);
                int bin0 = bin1 - 1;
                float bin_weight0 = magnitude * (bin1 - relative_orientation);
                float bin_weight1 = magnitude * (relative_orientation - bin0);
                bin0 = (bin0 + NUMBER_OF_BINS) % NUMBER_OF_BINS;
                bin1 = (bin1 + NUMBER_OF_BINS) % NUMBER_OF_BINS;

#if SPARTIAL_WEIGHTS
#if __PENCIL__
                __pencil_assume(cellxi < NUMBER_OF_CELLS);
                __pencil_assume(cellyi < NUMBER_OF_CELLS);
                __pencil_assume(cellxi >= 0);
                __pencil_assume(cellyi >= 0);
#endif
                if (cellyi >= 0 && cellxi >= 0) {
                    hist[i][cellyi][cellxi][bin0] += yscale0 * xscale0 * bin_weight0;
                    hist[i][cellyi][cellxi][bin1] += yscale0 * xscale0 * bin_weight1;
                }
                if (cellyi >= 0 && cellxi < NUMBER_OF_CELLS - 1) {
                    hist[i][cellyi][cellxi+1][bin0] += yscale0 * xscale1 * bin_weight0;
                    hist[i][cellyi][cellxi+1][bin1] += yscale0 * xscale1 * bin_weight1;
                }
                if (cellyi < NUMBER_OF_CELLS - 1 && cellxi >= 0) {
                    hist[i][cellyi+1][cellxi][bin0] += yscale1 * xscale0 * bin_weight0;
                    hist[i][cellyi+1][cellxi][bin1] += yscale1 * xscale0 * bin_weight1;
                }
                if (cellyi < NUMBER_OF_CELLS - 1 && cellxi < NUMBER_OF_CELLS - 1) {
                    hist[i][cellyi+1][cellxi+1][bin0] += yscale1 * xscale1 * bin_weight0;
                    hist[i][cellyi+1][cellxi+1][bin1] += yscale1 * xscale1 * bin_weight1;
                }
#elif NUMBER_OF_CELLS == 1
                hist[i][0][0][bin0] += bin_weight0;
                hist[i][0][0][bin1] += bin_weight1;
#else
                int cellxi = floor((pointx - minx) / cell_size);
                int cellyi = floor((pointy - miny) / cell_size);
#if __PENCIL__
                __pencil_assume(cellxi < NUMBER_OF_CELLS);
                __pencil_assume(cellyi < NUMBER_OF_CELLS);
                __pencil_assume(cellxi >= 0);
                __pencil_assume(cellyi >= 0);
#endif
                hist[i][cellxi][cellyi][bin0] += bin_weight0;
                hist[i][cellxi][cellyi][bin1] += bin_weight1;
#endif
            }
	    }
	    normalize(hist[i]);
    }
#pragma endscop

}

void pencil_hog( const int rows
               , const int cols
               , const int step
               , const uint8_t image[]
               , const int num_locations
               , const float location_x[]
               , const float location_y[]
               , const float block_size
               , float hist[]    //out
               ) {
    hog_multi(rows,cols,step,image,num_locations,location_x,location_y,block_size,hist);
}

