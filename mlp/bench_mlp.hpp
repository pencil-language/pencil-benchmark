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
#include "memory.hpp"
#include "mlp_impl.h"
#include "allocator.hpp"

namespace { struct hack_t; }

cMat /*float*/convertCVToMatFloat /*float*/( void * self, void * allocator, const cv::Mat_<double> & input );

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


// patchsize, m_wIn, m_wOut, m_U, hidden_num, rho2
std::tuple< std::vector<int>, std::vector<cMat /*float*/>, std::vector<cMat /*float*/>, std::vector<cMat /*float*/>, std::vector<int>, std::vector<double> >
convertHackToMlp ( void * self, void * allocator, const hack_t & hack )
{
    assert(hack.m_visibleLandmarks_size==hack.m_classifiers.size());
    assert(hack.m_visibleLandmarks_size==hack.responseMaps.size());

    int size = hack.m_visibleLandmarks_size;
    
    std::vector<int> m_patchSizes(size);
    std::vector<cMat /*float*/> m_wIns(size);
    std::vector<cMat /*float*/> m_wOuts(size);
    std::vector<cMat /*float*/> m_Us(size);
    std::vector<int> hidden_nums(size);
    std::vector<double> rho2s(size);    
    
    // we export each classifier
    for (int q=0; q<size; q++)
    {
        m_patchSizes[q] = hack.m_classifiers[q].m_patchSize;
        m_wIns[q]       = convertCVToMatFloat(self, allocator, hack.m_classifiers[q].m_wIn);
        m_wOuts[q]      = convertCVToMatFloat(self, allocator, hack.m_classifiers[q].m_wOut);
        m_Us[q]         = convertCVToMatFloat(self, allocator, hack.m_classifiers[q].m_U);
        hidden_nums[q]  = hack.m_classifiers[q].hidden_num;
        rho2s[q]        = hack.m_classifiers[q].rho2;
    } // for q in m_visibleLandmarks_size

    return std::make_tuple(m_patchSizes, m_wIns, m_wOuts, m_Us, hidden_nums, rho2s);
} // convertHackToMlp

void
freeClassifiers( void * self, void * allocator, mlp * classifiers[], int size )
{
    mlp * result = *classifiers;    
    for (int q=0; q<size; q++ )
        freeMLP( self, allocator, &(result[q]) );

    free(*classifiers);
    *classifiers=NULL;

    return;    
} // freeClassifiers


cMat /*uint8_t*/  convertCVToMatChar /*uint8_t*/  ( void * self, void * allocator, const cv::Mat_<uint8_t> & input )
{
    cMat /*uint8_t*/  result = CreateMatChar /*uint8_t*/ ( self, allocator, input.rows, input.cols );

    for ( int q=0; q<input.rows; q++)
        for ( int w=0; w<input.cols; w++ )
            reinterpret_cast<uint8_t*>(self)[ q * result.step + w + result.start ] = input(q,w);
    
    return result;    
} // convertCVTocMat /*uint8_t*/ 

cMat /*float*/convertCVToMatFloat /*float*/( void * self, void * allocator, const cv::Mat_<double> & input )
{
    cMat /*float*/result = CreateMatFloat( self, allocator, input.rows, input.cols );
    
    for ( int q=0; q<input.rows; q++)
        for ( int w=0; w<input.cols; w++ )
            reinterpret_cast<float*>(self)[ q * result.step + w + result.start ] = input(q,w);
    
    return result;    
} // convertCVToMatFloat


cv::Mat_<double> convertMatFloatToCV( void * self, cMat /*float*/input )
{
    cv::Mat_<double> result( input.rows, input.cols );
    
    for ( int q=0; q<input.rows; q++)
        for ( int w=0; w<input.cols; w++ )
            result(q,w) = reinterpret_cast<float*>(self)[ q * input.step + w + input.start ];
    
    return result;
} // convertMatFloatToCV


void allocateResponseMaps( void * self, void * allocator, int mapSize, int size, cMat /*float*/* responseMaps[] )
{
    *responseMaps = new cMat[size];
    assert(*responseMaps);    
    cMat /*float*/* result = *responseMaps;    

    for ( int q=0; q<size; q++ )
        result[q] = CreateMatFloat( self, allocator, 2 * mapSize + 1, 2 * mapSize + 1 );
    
    return;
}

void freeResponseMaps( void * self, void * allocator, cMat /*float*/* responseMaps[], int size )
{
    cMat /*float*/* result = *responseMaps;    
    assert(result);
    for ( int q=0; q<size; q++ )
        freeMatFloat( self, allocator, &(result[q]));

    delete [] result;
    *responseMaps=NULL;
    return;    
}

template <class T0>
auto
microseconds( T0 t0 ) -> decltype(std::chrono::duration_cast<std::chrono::microseconds>(t0).count())
{
    return std::chrono::duration_cast<std::chrono::microseconds>(t0).count();
}

// LuM end of file
