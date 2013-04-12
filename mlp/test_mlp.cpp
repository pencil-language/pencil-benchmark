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
const int memsize = 1.1 * MiB;

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

        boost::shared_array<char> buffer( new char[memsize] );
        carp::memory pool( memsize, uint8_t() );
        void * self = buffer.get();
        void * allocator = &pool;
        // here comes the function call
        {
            // preparing the inputs
            cMat /*uint8_t*/  alignedImage = convertCVToMatChar(self, allocator, conductor.hack.alignedImage);
            cMat /*float*/ shape = convertCVToMatFloat(self, allocator, conductor.hack.shape);
            auto m_classifiers = convertHackToMlp(self, allocator, conductor.hack);
            cMat /*float*/ * responseMaps;
            allocateResponseMaps( self, allocator, conductor.hack.m_mapSize, conductor.hack.m_visibleLandmarks_size, &responseMaps );
            
            auto start = std::chrono::high_resolution_clock::now();
            calculateMaps(
                self,
                allocator,
                conductor.hack.m_visibleLandmarks_size,
                conductor.hack.m_mapSize,
                alignedImage,
                shape,
                &(std::get<0>(m_classifiers)[0]), // &patchSizes[0]
                &(std::get<1>(m_classifiers)[0]), // &m_wIns[0]
                &(std::get<2>(m_classifiers)[0]), // &m_wOuts[0]
                &(std::get<3>(m_classifiers)[0]), // &m_Us[0]
                &responseMaps
                );
            auto end = std::chrono::high_resolution_clock::now();
            elapsed_time = microseconds(end - start);            
            
            // releasing the inputs
            freeMatChar( self, allocator, &alignedImage);
            freeMatFloat( self, allocator, &shape );
            for ( auto & q : std::get<1>(m_classifiers) ) freeMatFloat( self, allocator, &q );
            for ( auto & q : std::get<2>(m_classifiers) ) freeMatFloat( self, allocator, &q );
            for ( auto & q : std::get<3>(m_classifiers) ) freeMatFloat( self, allocator, &q );
            
            // !!!! freeMatFloat(&())
            //freeClassifiers(&m_classifiers, conductor.hack.m_classifiers.size());

            // converting the outputs
            std::vector< cv::Mat_<double> > calculatedResults;
            for (int q=0; q<conductor.hack.m_visibleLandmarks_size; q++)
            {
                cv::Mat_<double> nextResult;
                nextResult = convertMatFloatToCV( self, responseMaps[q] );
                calculatedResults.push_back(nextResult);                
            }
            
            // testing the output
            for (int q=0; q<conductor.hack.m_visibleLandmarks_size; q++)
            {
              // std::cout << "cv::norm( conductor.hack.responseMaps[" << q << "] - calculatedResults[" << q << "] ) = "
              //             << cv::norm( conductor.hack.responseMaps[q] - calculatedResults[q] ) << std::endl;
              assert(cv::norm( conductor.hack.responseMaps[q] - calculatedResults[q] ) < 0.00001);
            }
            
            // releasing the outputs
            freeResponseMaps( self, allocator, &responseMaps, conductor.hack.m_visibleLandmarks_size );

        }
        // here comes the test
        // PRINT(cv::norm( ));
    }
    
    std::cout << "total elapsed time = " << elapsed_time / 1000000. << " s." << std::endl;    
    //conductor.importer >> BOOST_SERIALIZATION_NVP(conductor.hack);
    
    return EXIT_SUCCESS;
} // int main


// LuM end of file
