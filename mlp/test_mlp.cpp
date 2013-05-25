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

#include <boost/smart_ptr.hpp>

#include "opencl.hpp"
#include "memory.hpp"
#include "bench_mlp.hpp"


int main()
{
    conductor_t conductor;
    int fail = 0;
    long int elapsed_time = 0;
    int64_t maxnetallocated = 0;
    int64_t maxgrossallocated = 0;
    
    
    for ( conductor.importer >> BOOST_SERIALIZATION_NVP(conductor.id);
          ((conductor.id != -1) and (conductor.id != 1));
          // conductor.id != -1;
          conductor.importer >> BOOST_SERIALIZATION_NVP(conductor.id)
        )
    {
        PRINT(conductor.id);
        conductor.importer >> BOOST_SERIALIZATION_NVP(conductor.hack);

        int groupsize = conductor.hack.m_visibleLandmarks_size;

        boost::shared_array<char> buffer( new char[groupsize * local_memsize] );
        
        // std::vector<carp::memory::buddy> pools( groupsize, carp::memory::buddy({local_memsize, uint8_t()}));
        std::vector<carp::memory::dense> pools( groupsize, carp::memory::dense({local_memsize, uint8_t()}));
        carp::memory::local_memory_manager locmm( groupsize * local_memsize, groupsize, local_memsize );
        
        void * self = buffer.get();
        std::vector<int> segments = locmm.get_segments();
       
        // here comes the function call
        {           
            auto calcpackages = convertHackToMlp( self, pools, segments, conductor.hack );

            for ( auto & pool : pools ) {
//                PRINT(pool.grossallocated());
//                PRINT(pool.netallocated());
                maxgrossallocated = std::max( maxgrossallocated, pool.grossallocated() );
                maxnetallocated = std::max( maxnetallocated, pool.netallocated() );
            }
                                    
            auto start = std::chrono::high_resolution_clock::now();
            calculateMaps(
                self,
                segments.data(),
                conductor.hack.m_visibleLandmarks_size,
                conductor.hack.m_mapSize,
                calcpackages.data() // ,
                );
            auto end = std::chrono::high_resolution_clock::now();
            elapsed_time += microseconds(end - start);

            // converting the outputs
            std::vector< cv::Mat_<double> > calculatedResults;
            for (int q=0; q<conductor.hack.m_visibleLandmarks_size; q++)
            {
                cv::Mat_<double> nextResult;
                nextResult = convertMatFloatToCV( self + segments[q], calcpackages[q].output.responseMap );
                calculatedResults.push_back(nextResult);
            }
            
            // testing the output
            for (int q=9; q<10/*conductor.hack.m_visibleLandmarks_size*/; q++)
            {
                // std::cout << "cv::norm( conductor.hack.responseMaps[" << q << "] - calculatedResults[" << q << "] ) = "
                //             << cv::norm( conductor.hack.responseMaps[q] - calculatedResults[q] ) << std::endl;
                PRINT(cv::norm( conductor.hack.responseMaps[q] - calculatedResults[q] ));
//                assert(cv::norm( conductor.hack.responseMaps[q] - calculatedResults[q] ) < 0.00001);
            }
            
        }
    }
    
    std::cout << "total elapsed time = " << elapsed_time / 1000000. << " s." << std::endl;    
    //conductor.importer >> BOOST_SERIALIZATION_NVP(conductor.hack);

    PRINT(maxnetallocated);
    PRINT(maxgrossallocated);
    
    return EXIT_SUCCESS;
} // int main


// LuM end of file
