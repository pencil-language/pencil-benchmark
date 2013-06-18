// UjoImro, 2013
// Experimental code for the CARP Project
// Copyright (c) RealEyes, 2013
// This version tests the responseMap calculation with input dumps

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

#include "cast.h"
#include "mlp.hpp"
#include "errors.hpp"
#include "memory.hpp"
#include "allocator.hpp"

namespace { struct hack_t; }

clMat /*float*/convertCVToMatFloat /*float*/( void * self, carp::memory::allocator & pool, const cv::Mat_<double> & input );
clMat /*uint8_t*/  convertCVToMatChar /*uint8_t*/  ( void * self, carp::memory::allocator & pool, const cv::Mat_<uint8_t> & input );

template <class T0>
void
printMatCV( cv::Mat_<T0> & mat, std::string name )
{
    std::cout << std::setprecision(6) << std::fixed;
    
    std::cout << name << " = [\n";

    int q,w;

    for (q=0; q<mat.rows; q++)
    {
        std::cout << "[ ";
        for( w=0; w<mat.cols; w++)
        {
            std::cout << static_cast<float>(mat(q,w)) << ", ";
        }
        std::cout << " ]\n";
    }
    
    std::cout << "]\n";
    
    return;
} // printMatFloat


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
        boost::archive::xml_iarchive importer;

    public:
        
        conductor_t() : id(0), dumpStream("response_dumps.xml", std::ios::in | std::ios::binary ), importer(dumpStream) {
            
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
        
        // result[q].tmp.xOuts          = CreateVectorMat( pools[q], gangsize );
        // for ( int w = 0; w<gangsize; w++ )
        // {
        //     clMat xOut = CreateMatFloat( pools[q], hack.m_classifiers[q].m_wIn.rows, 1 );
        //     ::SetMatToVector( self + memory_segments[q], result[q].tmp.xOuts, w, xOut );
        //     clMat test = ::GetMatFromVector( self + memory_segments[q], result[q].tmp.xOuts, w );
        //     assert(test.rows == xOut.rows );
        //     assert(test.cols == xOut.cols );
        //     assert(test.step == xOut.step );
        //     assert(test.start == xOut.start );
        // }

        
        result[q].output.responseMap = CreateMatFloat( pools[q], 2 * hack.m_mapSize + 1, 2 * hack.m_mapSize + 1 );
    } // for q in m_visibleLandmarks_size
    
    return result;
    
} // convertHackToMlp


clMat /*uint8_t*/  convertCVToMatChar /*uint8_t*/  ( void * self, carp::memory::allocator & pool, const cv::Mat_<uint8_t> & input )
{
    clMat /*uint8_t*/  result = CreateMatChar /*uint8_t*/ ( pool, input.rows, input.cols );

    for ( int q=0; q<input.rows; q++)
        for ( int w=0; w<input.cols; w++ )
            reinterpret_cast<uint8_t*>(self)[ q * result.step + w + result.start ] = input(q,w);
    
    return result;    
} // convertCVToclMat /*uint8_t*/ 

clMat /*float*/convertCVToMatFloat /*float*/( void * self, carp::memory::allocator & pool, const cv::Mat_<double> & input )
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

template <class allocator>
std::vector<clMat>
allocateResponseMaps( void * self, std::vector<allocator> & pools, int mapSize, int size )
{
    std::vector<clMat> result(size);

    for ( int q=0; q<size; q++ )
        result[q] = CreateMatFloat( pools[q], 2 * mapSize + 1, 2 * mapSize + 1 );
    
    return result;
}

template <class allocator>
void freeResponseMaps( void * self, std::vector<allocator> pools, std::vector<clMat> & responseMaps, int size )
{
    for ( int q=0; q<size; q++ )
        freeMatFloat( pools[q], &(responseMaps[q]) );

    return;    
}

template <class T0>
std::chrono::microseconds::rep
microseconds( T0 t0 )
{
    return std::chrono::duration_cast<std::chrono::microseconds>(t0).count();
}

// LuM end of file
