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
#error "This file is the implementation part of CLM/histo.h therefore it shouldn't be included directly."
#endif

#ifndef CLM_HISTO_IMPL_H
#define CLM_HISTO_IMPL_H

//#include "CLM/histo.h"

template <typename T>
gel::histo<T>::histo()
    : m_sum(0),
      m_N(0),
      m_histo(histo_size, 0)
{ }

template <typename T>
void gel::histo<T>::add(const histo_t& other)
{
    m_sum += other.m_sum;
    m_N += other.m_N;
    for ( int q=0; q < histo_size; q++ ) {
        m_histo[q] += other.m_histo[q];
    }
}

template <typename T>
void gel::histo<T>::subtract(const histo_t& other)
{
    m_sum -= other.m_sum;
    m_N -= other.m_N;
    for ( int q=0; q < histo_size; q++ ) {
        m_histo[q] -= other.m_histo[q];
    }
}

template <typename T>
void gel::histo<T>::insert(int value)
{
    m_histo[value]++;
    m_sum += value;
    m_N++;
}

template <typename T>
void gel::histo<T>::remove(int value)
{
    m_histo[value]--;
    m_sum -= value;
    m_N--;
}

template <typename T>
void gel::histo<T>::insert(const cv::Mat_<uint8_t>& sample)
{
    int step = sample.step1();
    const sample_type * sample_ptr = sample.ptr<sample_type>();

    for (int q = 0; q < sample.rows; q++) {
        for (int w = 0; w < sample.cols; w++) {
            insert(sample_ptr[ q * step + w ]);
        }
    }
}

template <typename T>
void gel::histo<T>::remove(const cv::Mat_<uint8_t>& sample)
{
    int step = sample.step1();
    const sample_type * sample_ptr = sample.ptr<sample_type>();

    for (int q = 0; q < sample.rows; q++) {
        for (int w = 0; w < sample.cols; w++) {
            remove(sample_ptr[ q * step + w ]);
        }
    }
}

template <typename T>
int gel::histo<T>::min()
{
    int q = 0;
    while ((m_histo[q]<=0) && (q < histo_size))
        q++;
    return q;
}

template <typename T>
int gel::histo<T>::max()
{
    int q = histo_size - 1;
    while ((m_histo[q]<=0) && (q > 0) )
        q--;
    return q;
}

template <typename T>
typename gel::histo<T>::value_type gel::histo<T>::sum()
{
    return m_sum;
}

template <typename T>
typename gel::histo<T>::value_type gel::histo<T>::mean()
{
    return m_sum / m_N;
}


#endif /* CLM_HISTO_IMPL_H */
