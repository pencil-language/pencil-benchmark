/* ---------------------------------------------------------------------
* Copyright Realeyes OU 2012-2014
*
* All rights reserved. Realeyes OU PROPRIETARY/CONFIDENTIAL.
* No part of this computer programs(s) may be used, reproduced,
* stored in any retrieval system, or transmitted, in any form or
* by any means, electronic, mechanical, photocopying,
* recording, or otherwise without prior written permission of
* Realeyes OU. Use is subject to license terms.
* ---------------------------------------------------------------------
*/

#ifndef HOGDESCRIPTOR_H
#error This is a private implementation file for HogDescriptor.h, do not include directly
#endif

#include "fast_approximate_math.h"

#include <algorithm>
#include <cfloat>
#include <cassert>
#include <array>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

template<int numberOfCells, int numberOfBins, bool gauss, bool spinterp, bool _signed>
int nel::HOGDescriptor<numberOfCells, numberOfBins, gauss, spinterp, _signed>::getNumberOfBins() {
    return numberOfCells * numberOfCells * numberOfBins;
}

template<int numberOfCells, int numberOfBins, bool gauss, bool spinterp, bool _signed>
float nel::HOGDescriptor<numberOfCells, numberOfBins, gauss, spinterp, _signed>::get_orientation(float mdy, float mdx) {
    if (_signed) {
        return normalized_atan2(mdy, mdx)*90.0f;
    } else {
        return normalized_atan(mdy / (mdx + FLT_EPSILON))*90.0f + 90.0f;
    }
}

template<int numberOfCells, int numberOfBins, bool gauss, bool spinterp, bool _signed>
void nel::HOGDescriptor<numberOfCells, numberOfBins, gauss, spinterp, _signed>::normalize(const cv::Mat_<float> &hist, cv::Mat_<float> &normalizedHist) {
    static const float l2hys_thrs = 0.15f;
    float scale = static_cast<float>(cv::norm(hist, cv::NORM_L2));
    if (scale < std::numeric_limits<float>::epsilon()) {
        normalizedHist = 0.0f;
    } else {
        cv::normalize(cv::min(hist / scale, l2hys_thrs), normalizedHist);
    }
}

template<int numberOfCells, int numberOfBins, bool gauss, bool spinterp, bool _signed>
cv::Mat_<float> nel::HOGDescriptor<numberOfCells, numberOfBins, gauss, spinterp, _signed>::compute( const cv::Mat_<uint8_t>  &image
                                                                                                  , const std::vector<float> &locationsx
                                                                                                  , const std::vector<float> &locationsy
                                                                                                  , const float              &blocksize
                                                                                                  ) {
    assert(locationsx.size() == locationsy.size());
    cv::Mat_<float> descriptors(locationsx.size(), getNumberOfBins(), 0.0f);

    tbb::parallel_for(tbb::blocked_range<size_t>(0, locationsx.size(), 5), [&](const tbb::blocked_range<size_t> range) {
    for (size_t n = range.begin(); n != range.end(); ++n) {
        const float centerx = locationsx[n];
        const float centery = locationsy[n];

        const float cellsize = blocksize / numberOfCells;
        const float halfblocksize = blocksize / 2.0f;

        const float minx = centerx - halfblocksize;
        const float miny = centery - halfblocksize;
        const float maxx = centerx + halfblocksize;
        const float maxy = centery + halfblocksize;

        const int minxi = std::max(fast_ceil(minx) , 1);
        const int minyi = std::max(fast_ceil(miny) , 1);
        const int maxxi = std::min(fast_floor(maxx), image.cols - 2);
        const int maxyi = std::min(fast_floor(maxy), image.rows - 2);

        cv::Matx<float, numberOfCells * numberOfCells, numberOfBins> hist(0.0f);

        float m1p2sigma2;
        if (gauss) {
            float sigma = halfblocksize;
            float sigma2 = sigma*sigma;
            m1p2sigma2 = -1.0f/(2.0f*sigma2);
            //float A = 1.0f / (2 * (float)M_PI*sigma2); // constant multiplier to all bin elements, but we normalize at the end, so we can skip this
        }

        // compute edges, magnitudes, orientations and the histogram
        for (int pointy = minyi; pointy <= maxyi; pointy++) {
            int cellyi;
            float yscale0;
            float yscale1;
            if (spinterp) {
                //Relative position of the pixel compared to the cell centers - y dimension
                float relative_pos_y = (pointy - miny) / cellsize - 0.5f;
                //Calculate the integral part and the fractional part of the relative position - y dimension
                cellyi = fast_floor(relative_pos_y);

                yscale1 = relative_pos_y - cellyi;
                yscale0 = 1.0f - yscale1;
            }

            float dy2;
            if (gauss) {
                float dy = pointy - centery;
                dy2 = dy*dy;
            }

            for (int pointx = minxi; pointx <= maxxi; pointx++) {
                float mdx = static_cast<float>(image(pointy, pointx + 1) - image(pointy, pointx - 1));
                float mdy = static_cast<float>(image(pointy + 1, pointx) - image(pointy - 1, pointx));

                float magnitude = std::sqrt(mdx*mdx + mdy*mdy);
                float orientation = nel::HOGDescriptor<numberOfCells, numberOfBins, gauss, spinterp, _signed>::get_orientation(mdy, mdx);   //GCC 4.6 bug??

                if (gauss) {
                    float dx = pointx - centerx;
                    float dx2 = dx*dx;
                    float B = std::exp((dx2 + dy2) * m1p2sigma2);
                    magnitude *= B;
                }

                // linear/trilinear interpolation of magnitudes
                float bin = (orientation - binSizeInDegrees / 2.0f) / binSizeInDegrees;
                int bin1 = fast_ceil(bin);
                int bin0 = bin1 - 1;
                float magscale0 = magnitude * (bin1 - bin);
                float magscale1 = magnitude * (bin - bin0);
                int binidx0 = (bin0 + numberOfBins) % numberOfBins;
                int binidx1 = (bin1 + numberOfBins) % numberOfBins;

                if (spinterp) {
                    //Relative position of the pixel compared to the cell centers - x dimension
                    float relative_pos_x = (pointx - minx) / cellsize - 0.5f;
                    //Calculate the integral part and the fractional part of the relative position - x dimension
                    int cellxi = fast_floor(relative_pos_x);

                    float xscale1 = relative_pos_x - cellxi;
                    float xscale0 = 1.0f - xscale1;

                    if (cellyi >= 0                && cellxi >= 0) {
                        hist((cellyi + 0) * numberOfCells + cellxi + 0, binidx0) += yscale0 * xscale0 * magscale0;
                        hist((cellyi + 0) * numberOfCells + cellxi + 0, binidx1) += yscale0 * xscale0 * magscale1;
                    }
                    if (cellyi >= 0                && cellxi < numberOfCells - 1) {
                        hist((cellyi + 0) * numberOfCells + cellxi + 1, binidx0) += yscale0 * xscale1 * magscale0;
                        hist((cellyi + 0) * numberOfCells + cellxi + 1, binidx1) += yscale0 * xscale1 * magscale1;
                    }
                    if (cellyi < numberOfCells - 1 && cellxi >= 0) {
                        hist((cellyi + 1) * numberOfCells + cellxi + 0, binidx0) += yscale1 * xscale0 * magscale0;
                        hist((cellyi + 1) * numberOfCells + cellxi + 0, binidx1) += yscale1 * xscale0 * magscale1;
                    }
                    if (cellyi < numberOfCells - 1 && cellxi < numberOfCells - 1) {
                        hist((cellyi + 1) * numberOfCells + cellxi + 1, binidx0) += yscale1 * xscale1 * magscale0;
                        hist((cellyi + 1) * numberOfCells + cellxi + 1, binidx1) += yscale1 * xscale1 * magscale1;
                    }
                } else {
                    if (numberOfCells == 1) {
                        hist(0, binidx0) += magscale0;
                        hist(0, binidx1) += magscale1;
                    } else {
                        int cellxi = fast_floor((pointx - minx) / cellsize);
                        int cellyi = fast_floor((pointy - miny) / cellsize);
                        assert(cellxi < numberOfCells);
                        assert(cellyi < numberOfCells);
                        assert(cellxi >= 0);
                        assert(cellyi >= 0);
                        hist(cellyi * numberOfCells + cellxi, binidx0) += magscale0;
                        hist(cellyi * numberOfCells + cellxi, binidx1) += magscale1;
                    }
                }
            }
        }
        auto row = descriptors.row(n);
        nel::HOGDescriptor<numberOfCells, numberOfBins, gauss, spinterp, _signed>::normalize(cv::Mat(hist, false).reshape(0,1), row); //GCC 4.6 bug??
    }});
    return descriptors;
}
