// UjoImro, 2013
// Experimental Code for the CARP Projects

#include <stdlib.h>
#include <boost/random.hpp>
#include <opencv2/core/core.hpp>
#include <boost/preprocessor.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/generator_iterator.hpp>

#include "mlp_impl.h"
#include "opencl.hpp"
#include "memory.hpp"
#include "utility.hpp"
#include "bench_mlp.hpp"


const int rows = 5;
const int cols = 7;
const int numel = local_memsize / sizeof(float);

class generateRandomMat {
    
private:
    typedef boost::mt19937 RNGType;
    RNGType rng;
    boost::uniform_int<> allocsize;
    boost::variate_generator< RNGType, boost::uniform_int<> > dice;

    static uint64_t fetchRandomSeed() {
        std::ifstream input("/dev/urandom");
        uint64_t seed = 0;        
        input.read(reinterpret_cast<char*>(&seed), sizeof(seed));

        return seed;        
    }   
    
public:
    
    generateRandomMat() : rng( generateRandomMat::fetchRandomSeed() ), allocsize(1, 10000), dice(rng, allocsize) { }

    clMat operator() ( void * self, carp::memory::allocator & pool, int64_t rows, int64_t cols ) {
        clMat result = CreateMatFloat( pool, rows, cols );

        int q, w;
        for ( q=0; q<result.rows; q++ )
            for ( w=0; w<result.cols; w++ )
                ((float*)self)[ q * result.step + w + result.start ] = 1.0f/dice();        
        
        return result;
    } // operator()
    
}; // class generateRandomMat

    

int main()
{

   

    // opencl initalization
    carp::opencl::device device;
    device.compile( {"operators.cl"}, {"transposeFloat", "expFloat", "expVecFloat"} );
    carp::memory::dense pool({local_memsize, uint8_t()});

    boost::shared_array<uint8_t> buffer( new uint8_t[local_memsize] );    
    void * self = buffer.get();    

    // random variable initialization
    generateRandomMat generator;

    clMat /* float */ exper     = generator( self, pool, rows, cols );
    clMat /* float */ transper  = generator( self, pool, cols, rows );
    clMat /* float */ sample    = generator( self, pool, rows, cols );

    clVector mats = CreateVectorMat( pool, 33 );

    for ( int q=0; q<mats.size; q++) {
        clMat buf = generator( self, pool, rows, cols );
        SetMatToVector( self, mats, q, buf );
    }
    
    printMatFloat(self, sample, "sample" );
    
    carp::opencl::array<uint8_t> clSelf( device, local_memsize, self );

    // expFloat
//    expFloat( self, sample, exper );
    device["expFloat"]( clSelf.cl(), sample, exper ).groupsize({1},{1});

    // transposeFloat
    transposeFloat( self, sample, transper );
    device["transposeFloat"]( clSelf.cl(), sample, transper ).groupsize({1},{1});

    // // repeated expFloat
    // for (int q = 0; q<33; q++ ) {
    //     clMat buf = GetMatFromVector( self, mats, q );
    //     expFloat( self, sample, buf );
    // }

    device["expVecFloat"]( clSelf.cl(), sample, mats ).groupsize({1}, {1});
    
    auto processed = clSelf.get();
    void * results = reinterpret_cast<void*>(processed.data());
    
    cv::Mat_<float> clExperVerify = convertMatFloatToCV( results, exper );
    cv::Mat_<float> clTransposeVerify = convertMatFloatToCV( results, transper );
    
    PRINT( cv::norm( clExperVerify - convertMatFloatToCV( self, exper)));
//    assert( cv::norm( clExperVerify - convertMatFloatToCV( self, exper)) < 0.0001 );

    PRINT( cv::norm( clTransposeVerify - convertMatFloatToCV( self, transper)));
//    assert( cv::norm( clTransposeVerify - convertMatFloatToCV( self, transper)) < 0.0001 );
    
    for (int q=0; q<33; q++) {
        clMat buf = GetMatFromVector( self, mats, q );
        clMat clBuf = GetMatFromVector( results, mats, q );
        assert( buf.rows  == clBuf.rows );
        assert( buf.cols  == clBuf.cols );
        assert( buf.step  == clBuf.step );
        assert( buf.start == clBuf.start );
        cv::Mat_<float> mat = convertMatFloatToCV( self, buf );
        cv::Mat_<float> clmat = convertMatFloatToCV( results, buf );

        assert(cv::norm( mat - clmat ) < 0.0001);
    }
    
    return EXIT_SUCCESS;
} // main



// LuM end of file
