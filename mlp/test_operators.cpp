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
    
    generateRandomMat() : rng( generateRandomMat::fetchRandomSeed() ), allocsize(0, 10), dice(rng, allocsize) { }

    clMat operator() ( void * self, carp::memory::allocator & pool, int64_t rows, int64_t cols ) {
        clMat result = CreateMatFloat( pool, rows, cols );

        int q, w;
        for ( q=0; q<result.rows; q++ )
            for ( w=0; w<result.cols; w++ )
                ((float*)self)[ q * result.step + w + result.start ] = dice();        
        
        return result;
    } // operator()
    
}; // class generateRandomMat

    

int main()
{

   

    // opencl initalization
    carp::opencl::device device;
    device.compile( {"operators.cl"}, {"transposeFloat", "expFloat"} );
    carp::memory::dense pool({local_memsize, uint8_t()});

    boost::shared_array<uint8_t> buffer( new uint8_t[local_memsize] );    
    void * self = buffer.get();    

    // random variable initialization
    generateRandomMat generator;

    
    clMat /* float */ sample = generator( self, pool, rows, cols );
    clMat /* float */ exper  = generator( self, pool, rows, cols );

    carp::opencl::array<uint8_t> clSelf( device, local_memsize, self );
        
    expFloat( clSelf, sample, exper );

    device["expFloat"]( clSelf.cl(), sample, exper ).groupsize({1},{1});
    
    auto processed = clSelf.get();
    void * results = reinterpret_cast<void*>(processed.data);
        
    cv::Mat_<float> clExperVerify = convertMatFloatToCV( results, expert );
    
    PRINT( cv::norm( clExperVerify - convertMatFloatToCV( self, exper)));
    assert( cv::norm( clExperVerify - convertMatFloatToCV( self, exper)) < 0.0001 );
    
        
    return EXIT_SUCCESS;    
} // main



// LuM end of file
