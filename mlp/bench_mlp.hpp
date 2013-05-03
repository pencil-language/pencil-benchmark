// UjoImro, 2013
// Experimental code for the CARP Project
// Copyright (c) RealEyes, 2013
// This version tests the responseMap calculation with input dumps

#include <tuple>
#include <chrono>
#include <string>
#include <stdlib.h>
#include <opencv2/core/core.hpp>
#include <boost/preprocessor.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include "cast.h"
#include "mlp.hpp"
#include "errors.hpp"
#include "memory.hpp"
#include "mlp_impl.h"
#include "allocator.hpp"

namespace { struct hack_t; }

clMat /*float*/convertCVToMatFloat /*float*/( void * self, carp::memory & pool, const cv::Mat_<double> & input );
clMat /*uint8_t*/  convertCVToMatChar /*uint8_t*/  ( void * self, carp::memory & pool, const cv::Mat_<uint8_t> & input );

namespace {
    
    struct hack_t {
        int m_visibleLandmarks_size;
        int m_mapSize;
        
        cv::Mat_<double> shape;
        cv::Mat_<uint8_t> alignedImage;
        std::vector<gel::MLP<double> > m_classifiers;
        std::vector<cv::Mat_<double> > responseMaps;
        
        template <class MT0>
        void serialize( MT0 & archiver, unsigned int ) {
            GEL_EXPORT_VAR( archiver, m_visibleLandmarks_size );
            GEL_EXPORT_VAR( archiver, m_mapSize );
            GEL_EXPORT_VAR( archiver, shape );
            GEL_EXPORT_VAR( archiver, alignedImage );
            GEL_EXPORT_VAR( archiver, m_classifiers );
            GEL_EXPORT_VAR( archiver, responseMaps );
        } // serialize
    
    }; // struct hack_t
    
    class conductor_t {
    public:
        int id;
        hack_t hack;
        std::ifstream dumpStream;
        boost::archive::binary_iarchive importer;

    public:
        
        conductor_t() : id(0), dumpStream("response_dumps.bin", std::ios::in | std::ios::binary ), importer(dumpStream) {
            
        }; // conductor_t
        
    }; // conductor_t

} // unnamed namespace 

namespace {

    void SetMatToVector( void * self, clVector /*clMat*/ vec, int idx, clMat elem )
    {
        assert(self);
        assert(idx>=0);
        assert(idx<vec.size);
    
        ((clMat*)self)[ vec.start + vec.step * idx ] = elem;
        return;        
    } // GetMatFromVector

    clMat GetMatFromVector( void * self, clVector /*clMat*/ vec, int idx )
    {
        assert(self);
        assert(idx>=0);
        assert(idx<vec.size);
    
        return ((clMat*)self)[ vec.start + vec.step * idx ];
    } // GetMatFromVector

} // unnamed namespace 


// patchsize, m_wIn, m_wOut, m_U, hidden_num, rho2
std::vector<calcpackage>
convertHackToMlp ( void * self, std::vector<carp::memory> & pools, std::vector<int> & memory_segments, const hack_t & hack )
{
    assert(hack.m_visibleLandmarks_size==hack.m_classifiers.size());
    assert(hack.m_visibleLandmarks_size==hack.responseMaps.size());

    int size = hack.m_visibleLandmarks_size;
    std::vector<calcpackage> result(size);
    
    // std::vector<int> m_patchSizes(size);
    // std::vector<clMat /*float*/> m_wIns(size);
    // std::vector<clMat /*float*/> m_wOuts(size);
    // std::vector<clMat /*float*/> m_Us(size);
    // std::vector<int> hidden_nums(size);
    // std::vector<double> rho2s(size);
    
    // we export each classifier
    for (int q=0; q<size; q++)
    {
        result[q].input.alignedImage = convertCVToMatChar( self + memory_segments[q], pools[q], hack.alignedImage );
        result[q].input.shape        = convertCVToMatFloat( self + memory_segments[q], pools[q], hack.shape );
        result[q].input.m_patchSize  = hack.m_classifiers[q].m_patchSize;
        result[q].input.m_wIn        = convertCVToMatFloat( self + memory_segments[q], pools[q], hack.m_classifiers[q].m_wIn );
        result[q].input.m_U          = convertCVToMatFloat( self + memory_segments[q], pools[q], hack.m_classifiers[q].m_U );
        result[q].input.m_wOut       = convertCVToMatFloat( self + memory_segments[q], pools[q], hack.m_classifiers[q].m_wOut );
        result[q].tmp.wIn            = CreateMatFloat ( pools[q], result[q].input.m_wIn.rows, result[q].input.m_U.rows );
        result[q].tmp.patches        = CreateVectorMat( pools[q], size );
        for ( int w = 0; w<size; w++ )
        {
            clMat patch = CreateMatFloat( pools[q], 2 * hack.m_classifiers[q].m_patchSize + 1, 2 * hack.m_classifiers[q].m_patchSize + 1 );
            ::SetMatToVector( self + memory_segments[q], result[q].tmp.patches, w, patch );
        }
        
        result[q].tmp.xOuts          = CreateVectorMat( pools[q], size );
        for ( int w = 0; w<size; w++ )
        {
            clMat xOut = CreateMatFloat( pools[q], result[q].input.m_wIn.rows, 1 );
            ::SetMatToVector( self + memory_segments[q], result[q].tmp.xOuts, w, xOut );
            clMat test = ::GetMatFromVector( self + memory_segments[q], result[q].tmp.xOuts, w );
            assert(test.rows == xOut.rows );
            assert(test.cols == xOut.cols );
            assert(test.step == xOut.step );
            assert(test.start == xOut.start );
        }

        result[q].tmp.es             = CreateVectorMat( pools[q], size );
        for ( int w = 0; w<size; w++ )
        {
            clMat e = CreateMatFloat( pools[q], result[q].input.m_wIn.rows, 1 );
            ::SetMatToVector( self + memory_segments[q], result[q].tmp.es, w, e );
        }
        
        // result[q].input.m_U         = convertCVToMatFloat( self + memory_segments[q] * sizeof(float), pools[q], hack.m_classifiers[q].m_U );
        // hidden_nums[q]  = hack.m_classifiers[q].hidden_num;
        // rho2s[q]        = hack.m_classifiers[q].rho2;
        result[q].output.responseMap = CreateMatFloat( pools[q], 2 * hack.m_mapSize + 1, 2 * hack.m_mapSize + 1 );
    } // for q in m_visibleLandmarks_size
    
    return result;
    
} // convertHackToMlp

// void
// freeClassifiers( void * self, carp::memory & pool, mlp * classifiers[], int size )
// {
//     mlp * result = *classifiers;    
//     for (int q=0; q<size; q++ )
//         freeMLP( self, pool[q], &(result[q]) );

//     free(*classifiers);
//     *classifiers=NULL;

//     return;    
// } // freeClassifiers


clMat /*uint8_t*/  convertCVToMatChar /*uint8_t*/  ( void * self, carp::memory & pool, const cv::Mat_<uint8_t> & input )
{
    clMat /*uint8_t*/  result = CreateMatChar /*uint8_t*/ ( pool, input.rows, input.cols );

    for ( int q=0; q<input.rows; q++)
        for ( int w=0; w<input.cols; w++ )
            reinterpret_cast<uint8_t*>(self)[ q * result.step + w + result.start ] = input(q,w);
    
    return result;    
} // convertCVToclMat /*uint8_t*/ 

clMat /*float*/convertCVToMatFloat /*float*/( void * self, carp::memory & pool, const cv::Mat_<double> & input )
{
    clMat /*float*/result = CreateMatFloat( pool, input.rows, input.cols );
    
    for ( int q=0; q<input.rows; q++)
        for ( int w=0; w<input.cols; w++ )
            reinterpret_cast<float*>(self)[ q * result.step + w + result.start ] = input(q,w);
    
    return result;    
} // convertCVToMatFloat


cv::Mat_<double> convertMatFloatToCV( void * self, clMat /*float*/input )
{
    cv::Mat_<double> result( input.rows, input.cols );
    
    for ( int q=0; q<input.rows; q++)
        for ( int w=0; w<input.cols; w++ )
            result(q,w) = reinterpret_cast<float*>(self)[ q * input.step + w + input.start ];
    
    return result;
} // convertMatFloatToCV


std::vector<clMat>
allocateResponseMaps( void * self, std::vector<carp::memory> & pools, int mapSize, int size )
{
    std::vector<clMat> result(size);

    for ( int q=0; q<size; q++ )
        result[q] = CreateMatFloat( pools[q], 2 * mapSize + 1, 2 * mapSize + 1 );
    
    return result;
}

void freeResponseMaps( void * self, std::vector<carp::memory> pools, std::vector<clMat> & responseMaps, int size )
{
    for ( int q=0; q<size; q++ )
        freeMatFloat( pools[q], &(responseMaps[q]) );

    return;    
}

template <class T0>
auto
microseconds( T0 t0 ) -> decltype(std::chrono::duration_cast<std::chrono::microseconds>(t0).count())
{
    return std::chrono::duration_cast<std::chrono::microseconds>(t0).count();
}

// LuM end of file
