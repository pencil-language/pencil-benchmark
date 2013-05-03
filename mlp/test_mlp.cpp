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

const int KiB=1024;
const int MiB=1024*KiB;
// const int memsize = 1.1 * MiB;
const int local_memsize = 8 * 64 * KiB;

int main()
{
    conductor_t conductor;
    int fail = 0;
    long int elapsed_time = 0;
    
    for ( conductor.importer >> BOOST_SERIALIZATION_NVP(conductor.id);
          ((conductor.id != -1) and (conductor.id != 25));
          // conductor.id != -1;
          conductor.importer >> BOOST_SERIALIZATION_NVP(conductor.id)
        )
    {
        PRINT(conductor.id);
        conductor.importer >> BOOST_SERIALIZATION_NVP(conductor.hack);

        int groupsize = conductor.hack.m_visibleLandmarks_size;

        boost::shared_array<char> buffer( new char[groupsize * local_memsize] );
        std::vector<carp::memory> pools( groupsize, {local_memsize, uint8_t()} );
        carp::local_memory_manager locmm( groupsize * local_memsize, groupsize, local_memsize );
        
        void * self = buffer.get();
        std::vector<int> segments = locmm.get_segments();

        // here comes the function call
        {
            // preparing the inputs
//            std::vector<cMat /* uint8_t */> alignedImages(pools.size());
//            std::vector<cMat /* float */> shapes(pools.size());

            // for (int q=0; q<pools.size(); q++)
            // {
            //     alignedImages[q] = convertCVToMatChar( self, pools[q], conductor.hack.alignedImage );
            //     shapes[q]        = convertCVToMatChar( self, pools[q], conductor.hack.shape );
            // }
            
            // cMat /*uint8_t*/  alignedImage = convertCVToMatChar( self, pools, conductor.hack.alignedImage );
            
            // std::vector<cMat /*uint8_t*/> alignedImages = allocateImage( self, pools, alignedImage );
            // std::vector<cMat /*float*/> shapes = allocateImage( self, pools, conductor.hack.shape );
            
            auto calcpackages = convertHackToMlp( self, pools, segments, conductor.hack );
            // cMat /*float*/ * responseMaps;
            
//            std::vector<clMat /*float*/> responseMaps = allocateResponseMaps( self, pools, conductor.hack.m_mapSize, conductor.hack.m_visibleLandmarks_size );
            
            auto start = std::chrono::high_resolution_clock::now();
            calculateMaps(
                self,
                &(segments[0]),
                conductor.hack.m_visibleLandmarks_size,
                conductor.hack.m_mapSize,
                // &(alignedImages[0]),
                // &(shapes[0]),
                &(calcpackages[0]) // ,
                // &(std::get<0>(m_classifiers)[0]), // &patchSizes[0]
                // &(std::get<1>(m_classifiers)[0]), // &m_wIns[0]
                // &(std::get<2>(m_classifiers)[0]), // &m_wOuts[0]
                // &(std::get<3>(m_classifiers)[0]), // &m_Us[0]
                // responseMaps
                );
            auto end = std::chrono::high_resolution_clock::now();
            elapsed_time = microseconds(end - start);

            // inputs will be released automatically, when the pool is destroyed
            // // releasing the inputs
            // freeMatChar( self, allocator, &alignedImage);
            // freeMatFloat( self, allocator, &shape );
            // for ( auto & q : std::get<1>(m_classifiers) ) freeMatFloat( self, allocator, &q );
            // for ( auto & q : std::get<2>(m_classifiers) ) freeMatFloat( self, allocator, &q );
            // for ( auto & q : std::get<3>(m_classifiers) ) freeMatFloat( self, allocator, &q );
            
            // !!!! freeMatFloat(&())
            //freeClassifiers(&m_classifiers, conductor.hack.m_classifiers.size());

            // converting the outputs
            std::vector< cv::Mat_<double> > calculatedResults;
            for (int q=0; q<conductor.hack.m_visibleLandmarks_size; q++)
            {
                cv::Mat_<double> nextResult;
                nextResult = convertMatFloatToCV( self, calcpackages[q].output.responseMap );
                calculatedResults.push_back(nextResult);                
            }
            
            // testing the output
            for (int q=0; q<conductor.hack.m_visibleLandmarks_size; q++)
            {
                // std::cout << "cv::norm( conductor.hack.responseMaps[" << q << "] - calculatedResults[" << q << "] ) = "
                //             << cv::norm( conductor.hack.responseMaps[q] - calculatedResults[q] ) << std::endl;
                PRINT(cv::norm( conductor.hack.responseMaps[q] - calculatedResults[q] ));
//                assert(cv::norm( conductor.hack.responseMaps[q] - calculatedResults[q] ) < 0.00001);
            }
            
            // releasing the outputs
            // freeResponseMaps( self, pools, &responseMaps, conductor.hack.m_visibleLandmarks_size );

        }
        // here comes the test
        // PRINT(cv::norm( ));
    }
    
    std::cout << "total elapsed time = " << elapsed_time / 1000000. << " s." << std::endl;    
    //conductor.importer >> BOOST_SERIALIZATION_NVP(conductor.hack);
    
    return EXIT_SUCCESS;
} // int main


// LuM end of file
