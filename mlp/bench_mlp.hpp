// UjoImro, 2013
// Experimental code for the CARP Project
// Copyright (c) RealEyes, 2013
// This version tests the responseMap calculation with input dumps

#ifndef BENCH_MLP__HPP__
#define BENCH_MLP__HPP__

#include <tuple>
#include <chrono>
#include <string>
#include <iomanip>
#include <stdlib.h>
#include <opencv2/core/core.hpp>
#include <boost/preprocessor.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include "mlp.hpp"
#include "errors.hpp"
#include "memory.hpp"
#include "allocator.hpp"
#include "serialization.hpp"

namespace carp {

    clMat /*float*/    convertCVToMatFloat /*float*/( void * self, carp::memory::allocator & pool, const cv::Mat_<double> & input );
    clMat /*uint8_t*/  convertCVToMatChar /*uint8_t*/  ( void * self, carp::memory::allocator & pool, const cv::Mat_<uint8_t> & input );

    // patchsize, m_wIn, m_wOut, m_U, hidden_num, rho2
    template <class allocator>
    std::vector<calcpackage>
    convertHackToMlp ( char * self, std::vector<allocator> & pools, std::vector<int> & memory_segments, const hack_t & hack )
    {
        assert(hack.m_visibleLandmarks_size==hack.m_classifiers.size());
        assert(hack.m_visibleLandmarks_size==hack.responseMaps.size());

        int size = hack.m_visibleLandmarks_size;
        std::vector<calcpackage> result(size);

        // we export each classifier
        for (int q=0; q<size; q++)
        {
            result[q].input.alignedImage = convertCVToMatChar( self + memory_segments[q], pools[q], hack.alignedImage );
            result[q].input.shape        = convertCVToMatFloat( self + memory_segments[q], pools[q], hack.shape );

            result[q].input.m_patchSize  = hack.m_classifiers[q].m_patchSize;

            result[q].input.wIn = convertCVToMatFloat( self + memory_segments[q], pools[q],
                                                       hack.m_classifiers[q].m_wIn( cv::Range(0, hack.m_classifiers[q].m_wIn.rows ), cv::Range(0, hack.m_classifiers[q].m_wIn.cols -1 ) )
                                                       * hack.m_classifiers[q].m_U.t() );

            result[q].input.m_wOut = convertCVToMatFloat( self + memory_segments[q], pools[q], hack.m_classifiers[q].m_wOut );
            result[q].input.bIn = convertCVToMatFloat( self + memory_segments[q], pools[q],
                                                       hack.m_classifiers[q].m_wIn( cv::Range(0, hack.m_classifiers[q].m_wIn.rows), cv::Range(hack.m_classifiers[q].m_wIn.cols - 1, hack.m_classifiers[q].m_wIn.cols )));

            result[q].output.responseMap = carp::CreateMatFloat( pools[q], 2 * hack.m_mapSize + 1, 2 * hack.m_mapSize + 1 );
        } // for q in m_visibleLandmarks_size

        return result;

    }

    cv::Mat_<double> convertMatFloatToCV( void * self, clMat /*float*/input )
    {
        cv::Mat_<double> result( input.rows, input.cols );

        for ( int q=0; q<input.rows; q++)
            for ( int w=0; w<input.cols; w++ )
                result(q,w) = reinterpret_cast<float*>(self)[ q * input.step + w + input.start ];

        return result;
    }

    clMat /*uint8_t*/  convertCVToMatChar /*uint8_t*/  ( void * self, carp::memory::allocator & pool, const cv::Mat_<uint8_t> & input )
    {
        clMat /*uint8_t*/  result = carp::CreateMatChar /*uint8_t*/ ( pool, input.rows, input.cols );

        for ( int q=0; q<input.rows; q++)
            for ( int w=0; w<input.cols; w++ )
                reinterpret_cast<uint8_t*>(self)[ q * result.step + w + result.start ] = input(q,w);

        return result;
    }

    clMat /*float*/convertCVToMatFloat /*float*/( void * self, carp::memory::allocator & pool, const cv::Mat_<double> & input )
    {
        clMat /*float*/result = carp::CreateMatFloat( pool, input.rows, input.cols );

        for ( int q=0; q<input.rows; q++)
            for ( int w=0; w<input.cols; w++ )
                reinterpret_cast<float*>(self)[ q * result.step + w + result.start ] = input(q,w);

        return result;
    }
}

#endif