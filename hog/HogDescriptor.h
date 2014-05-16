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
#define HOGDESCRIPTOR_H

#include <opencv2/core/core.hpp>

namespace nel {
    template<int numberOfCells, int numberOfBins, bool gauss, bool spinterp, bool _signed>
    class HOGDescriptor {
        static_assert(numberOfCells > 1 || !spinterp, "Cannot apply spatial interpolation with only one cell.");
    public:
        static cv::Mat_<float> compute( const cv::Mat_<uint8_t>  &img
                                      , const std::vector<float> &locationsx
                                      , const std::vector<float> &locationsy
                                      , const float              &blocksize
                                      );
        static int getNumberOfBins();

    private:
        static float get_orientation(float mdy, float mdx);
        static void normalize(const cv::Mat_<float> &hist, cv::Mat_<float> &normalizedHist);

    private:
        static const int binSizeInDegrees = (_signed ? 360 : 180) / numberOfBins;
    };
}

#include "HogDescriptor.hpp"

#endif
