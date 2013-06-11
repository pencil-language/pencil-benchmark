/* ---------------------------------------------------------------------
* Copyright Realeyes OU 2012
*
* All rights reserved. Realeyes OU PROPRIETARY/CONFIDENTIAL.
* No part of this computer programs(s) may be used, reproduced,
* stored in any retrieval system, or transmitted, in any form or
* by any means, electronic, mechanical, photocopying,
* recording, or otherwise without prior written permission of
* Realeyes OU. Use is subject to license terms.
* ---------------------------------------------------------------------
*/

#ifndef CLM_HISTO_H
#define CLM_HISTO_H

//#include "Core/defines.h"

#include <opencv2/core/core.hpp>

namespace gel
{

/*!
    \brief windowed mean, min and max using histograms
*/
template <typename T>
class histo
{
public:
    static const int histo_size = 256;

    typedef T value_type;
    typedef histo<T> histo_t;
    typedef cv::Mat_<uint8_t>::value_type sample_type;

    histo();

    void add(const histo_t& other);

    void subtract(const histo_t& other );

    void insert(int value);

    void remove(int value);

    void insert(const cv::Mat_<uint8_t>& sample);

    void remove(const cv::Mat_<uint8_t>& sample);

    int min();

    int max();

    value_type sum();

    value_type mean();

private:
    value_type m_sum; // sum of the elements of the sample
    int m_N; // the number of elements of the sample
    std::vector<int> m_histo; // histogram of the elements of the sample
};

}

#include "histo_impl.h"

#endif /* CLM_HISTO_H */
