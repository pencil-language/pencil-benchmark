// UjoImro, 2013
// Experimental code for the CARP Project
// Copyright (c) RealEyes, 2013
// This version tests the responseMap calculation with input dumps


/*
  extern int EF_ALIGNMENT = 0;
  extern int EF_PROTECT_BELOW = 0;
  extern int EF_PROTECT_FREE = 0;
  extern int EF_ALLOW_MALLOC_0 = 1;
  extern int EF_FILL = 1922;
*/

#include <vector>
#include <iomanip>
#include <boost/smart_ptr.hpp>

#include "opencl.hpp"
#include "memory.hpp"
#include "bench_mlp.hpp"

// OpenCL includes
#include "mlp_impl.clh"

const int processed_frames = 100;
const int local_memsize = 21 * KiB;

namespace {
    struct opt_t {
        opt_t() : gangsize(-1), fps(0), time(0) {  }
        
        int gangsize;
        float fps;
        float time;    
    }; // opt_t
} // unnamed namespace

    

int main()
{
    std::vector<hack_t> packages;    
    conductor_t conductor; // the class for importing the input from the clm
    carp::opencl::device device;
    device.source_compile( mlp_impl_cl, mlp_impl_cl_len, carp::string_vector("calculateMaps") );
        
    int fail = 0;
    
    for ( conductor.importer >> BOOST_SERIALIZATION_NVP(conductor.id);
          ((conductor.id != -1) && (conductor.id != processed_frames));
          conductor.importer >> BOOST_SERIALIZATION_NVP(conductor.id)
        )
    {
        hack_t hack;
        
        PRINT(conductor.id);
        conductor.importer >> BOOST_SERIALIZATION_NVP(hack);
        packages.push_back(hack);
    } // for conductor

    PRINT("benchmarking");
    opt_t opt;
    
    for ( int gangsize = 32; gangsize<=/*640*/ 512; gangsize+=32 )
    {
        PRINT(gangsize);
        long int elapsed_time = 0;
        int64_t maxnetallocated = 0;
        int64_t maxgrossallocated = 0;        
        for ( auto & package : packages ) {
            int groupsize = package.m_visibleLandmarks_size;
            boost::shared_array<char> buffer( new char[groupsize * local_memsize] );

            std::vector<carp::memory::dense> pools( groupsize, carp::memory::dense(carp::memory::allocator::sizer(local_memsize, uint8_t())));
            carp::memory::local_memory_manager locmm( groupsize * local_memsize, groupsize, local_memsize );
                 
            char * self = buffer.get();
            std::vector<int> segments = locmm.get_segments();

            // here comes the function call
            {           
                auto calcpackages = convertHackToMlp( self, pools, segments, package );
                
                for ( auto & pool : pools ) {
                    maxgrossallocated = std::max( maxgrossallocated, pool.grossallocated() );
                    maxnetallocated = std::max( maxnetallocated, pool.netallocated() );
                }
                
                // preparing the opencl data
                carp::opencl::array<uint8_t> clSelf( device, groupsize * local_memsize, self );
                carp::opencl::array<int> clSegments( device, segments );
                carp::opencl::array<calcpackage> clCalcpackages( device, calcpackages );
                
                auto start = std::chrono::high_resolution_clock::now();
                device["calculateMaps"](
                    clSelf.cl(),
                    clSegments.cl(),
                    package.m_visibleLandmarks_size,
                    package.m_mapSize,
                    clCalcpackages.cl(),
                    local_memsize,
                    carp::opencl::buffer(local_memsize)
                    ).groupsize( carp::make_vector<size_t>(gangsize),carp::make_vector<size_t>( gangsize * package.m_visibleLandmarks_size ) );
                auto end = std::chrono::high_resolution_clock::now();
                elapsed_time += microseconds(end - start);
                
                // copying the data back to the CPU
                auto processed = clSelf.get();
                char * results = reinterpret_cast<char*>(processed.data());
                
                // converting the outputs            
                std::vector< cv::Mat_<double> > calculatedResults;
                for (int q=0; q<package.m_visibleLandmarks_size; q++)
                {
                    cv::Mat_<double> nextResult;
                    nextResult = convertMatFloatToCV( results + segments[q], calcpackages[q].output.responseMap );
                    calculatedResults.push_back(nextResult);                
                }
                
                // testing the output
                for (int q=0; q<package.m_visibleLandmarks_size; q++)
                {
                    if (cv::norm( package.responseMaps[q] - calculatedResults[q] ) > 0.0001) throw std::runtime_error("package.responseMaps[q] - calculatedResults[q] ) < 0.0001 failed");
                }
            }
        } // packages
        float fps = 1000000. * processed_frames / elapsed_time;
        if (opt.fps < fps) {
            opt.fps = fps;
            opt.gangsize = gangsize;
            opt.time = elapsed_time / 1000000.;            
        }
        
        
        std::cout << "total elapsed time = " << elapsed_time / 1000000. << " s." << std::endl;
//        std::cout << std::setprecision(2) << std::fixed;
        std::cout << "processing speed   = " << 1000000. * processed_frames / elapsed_time << "fps" << std::endl;
        PRINT(maxnetallocated);
        PRINT(maxgrossallocated);            
    } // gangsize 

    PRINT(opt.gangsize);
    PRINT(opt.fps);
    PRINT(opt.time);    
    
    return EXIT_SUCCESS;
} // int main


// LuM end of file
