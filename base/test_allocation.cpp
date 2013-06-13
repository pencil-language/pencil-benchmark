// UjoImro, 2013
// Experimental Research Code for the CARP Project

#include <fstream>
#include <iomanip>
#include <iostream>
#include <boost/random.hpp>
#include <boost/preprocessor.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/generator_iterator.hpp>

#include "errors.hpp"
#include "memory.hpp"

const int local_memsize = 48 * KiB;
const int numel = local_memsize / sizeof(float);
const bool debug = false;

struct chunk_t {
    chunk_t( int64_t & size, int64_t & pointer ) : size(size), pointer(pointer) { }
    
    int64_t size;
    int64_t pointer;    
}; // chunk_t

int main()
{
    carp::memory::buddy bpool( carp::memory::allocator::sizer(numel, float()));
    carp::memory::allocator & pool = bpool;
    
    // std::ifstream input("/dev/urandom");
    uint64_t seed = 0;
    // input.read(reinterpret_cast<char*>(&seed), sizeof(seed));

    typedef boost::mt19937 RNGType;
    RNGType rng(seed);
    boost::uniform_int<> allocsize(1, numel/10);
    boost::variate_generator< RNGType, boost::uniform_int<> > dice(rng, allocsize);

    for (int repeat = 0; repeat < 1000; repeat++ ) {
        PRINT(repeat);
        
        std::vector<chunk_t> arrays;

        int64_t real_allocated=0;

        // allocation 
        for ( int q=0; q<14; q++ ) {
            int64_t size = dice();
            int64_t pointer;        
            try {
                pointer = pool.allocate( carp::memory::allocator::sizer(size, float()));
            }
            catch ( carp::memory::exception & exception )
            {
                if (exception.error == carp::memory::INSUFFICIENT_MEMORY ) {
                    if (debug) {                        
                        std::cout << "insufficient memory at " << pool.grossallocated() << "/" << numel
                                  << " (requested size: " << size << ")" << std::endl;

                        std::cout << "real allocated memory: " << real_allocated << "/" << pool.grossallocated() << "=" << 100 * real_allocated / pool.grossallocated() << "%" << std::endl;
                    }                    
                }
                else throw exception;

                continue;            
            }

            arrays.push_back( chunk_t(size, pointer) );
            real_allocated += size;

            if (debug) std::cout << "allocated: " << size << " (size); " << real_allocated << " (net_size); "  << pool.grossallocated() << " (gross_size); " << pointer << " (pointer)" << std::endl;        
        }

        std::cout << "real allocated memory: " << real_allocated << "/" << pool.grossallocated() << "=" << 100 * real_allocated / pool.grossallocated() << "%" << std::endl;
        
        // release
        std::random_shuffle( arrays.begin(), arrays.end() );

        for ( auto & q : arrays )
        {        
            pool.release( carp::memory::allocator::sizer(q.pointer, float()));

            real_allocated -= q.size;

            if (debug)
                if (pool.grossallocated()>0)
                    std::cout << "real allocated memory: " << real_allocated << "/" << pool.grossallocated() << "=" << 100 * real_allocated / pool.grossallocated() << "%" << std::endl;
                else
                    std::cout << "real allocated memory: " << real_allocated << "/" << pool.grossallocated() << std::endl;
        }

        std::cout << "real allocated memory: " << real_allocated << "/" << pool.grossallocated() << std::endl;
    } // repeat
    
    
    return EXIT_SUCCESS;
} // main






























// LuM end of file
