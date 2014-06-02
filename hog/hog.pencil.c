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

static void hog_multi( const int rows
                     , const int cols
                     , const int step
                     , const uint8_t image[static const restrict rows][step]
                     , const int num_locations
                     , const float location_x[static const restrict num_locations]
                     , const float location_y[static const restrict num_locations]
                     , const float blck_size
                     , float hist[static const restrict num_locations][NUMBER_OF_CELLS][NUMBER_OF_CELLS][NUMBER_OF_BINS]    //out
                     ) {
#pragma scop

#pragma pencil independent
    for (int i = 0; i < num_locations; ++i) {
	    float cell_size = blck_size / NUMBER_OF_CELLS;
	    float minx = location_x[i] - blck_size / 2.0f;
	    float miny = location_y[i] - blck_size / 2.0f;
	    float maxx = location_x[i] + blck_size / 2.0f;
	    float maxy = location_y[i] + blck_size / 2.0f;

	    int minxi = max(ceil(minx), 1);
	    int minyi = max(ceil(miny), 1);
	    int maxxi = min(floor(maxx), cols - 2);
	    int maxyi = min(floor(maxy), rows - 2);

	    for (int j1 = 0; j1 < NUMBER_OF_CELLS; j1++)
	      for (int j2 = 0; j2 < NUMBER_OF_CELLS; j2++)
	        for (int j3 = 0; j3 < NUMBER_OF_BINS; j3++)
	          hist[i][j1][j2][j3] = 0;

#if GAUSSIAN_WEIGHTS
	    const float sigma = blck_size / 2.0f;
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
		int temp1 = pointx-1;
		int temp2 = pointy-1;
                float mdx = image[pointy][pointx+1] - image[pointy][temp1];
                float mdy = image[pointy+1][pointx] - image[temp2][pointx];

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
               , const float blck_size
               , float hist[]    //out
               ) {
    hog_multi(rows,cols,step,image,num_locations,location_x,location_y,blck_size,hist);
}

