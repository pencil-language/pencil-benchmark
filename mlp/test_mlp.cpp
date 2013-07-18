// UjoImro, 2013
// Experimental code for the CARP Project
// Copyright (c) RealEyes, 2013
// This version tests the responseMap calculation with input dumps

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
#include "mlp_impl.h"
#include "utility.hpp"
#include "serialization.hpp"

/*
  extern int EF_ALIGNMENT = 0;
  extern int EF_PROTECT_BELOW = 0;
  extern int EF_PROTECT_FREE = 0;
  extern int EF_ALLOW_MALLOC_0 = 1;
  extern int EF_FILL = 1922;
*/

#ifndef PROCESSED_FRAMES
#define PROCESSED_FRAMES 100
#endif
const int processed_frames = PROCESSED_FRAMES;

int main()
{
    carp::conductor_t conductor;
    int fail = 0;
    long int elapsed_time = 0;
    
    
    for ( conductor.importer >> BOOST_SERIALIZATION_NVP(conductor.id);
          ((conductor.id != -1) && (conductor.id != processed_frames));
          // conductor.id != -1;
          conductor.importer >> BOOST_SERIALIZATION_NVP(conductor.id)
        )
    {
        PRINT(conductor.id);
        conductor.importer >> BOOST_SERIALIZATION_NVP(conductor.hack);

        // here comes the function call
        {
            // preparing the inputs
            MatChar alignedImage = carp::convertCVToMatChar(conductor.hack.alignedImage);
            MatFloat shape = carp::convertCVToMatFloat(conductor.hack.shape);
            mlp * m_classifiers = carp::convertHackToMlp(conductor.hack);
            MatFloat * responseMaps;
            carp::allocateResponseMaps( conductor.hack.m_mapSize, conductor.hack.m_visibleLandmarks_size, &responseMaps );

            auto start = std::chrono::high_resolution_clock::now();
            calculateMaps(
                conductor.hack.m_visibleLandmarks_size,
                conductor.hack.m_mapSize,
                alignedImage,
                shape,
                m_classifiers,
                &responseMaps
                );
            auto end = std::chrono::high_resolution_clock::now();
            elapsed_time += carp::microseconds(end - start);

            // releasing the inputs
            freeMatChar(&alignedImage);
            freeMatFloat(&shape);
            carp::freeClassifiers(&m_classifiers, conductor.hack.m_classifiers.size());

            // converting the outputs
            std::vector< cv::Mat_<double> > calculatedResults;
            for (int q=0; q<conductor.hack.m_visibleLandmarks_size; q++)
            {
                cv::Mat_<double> nextResult;
                nextResult = carp::convertMatFloatToCV( responseMaps[q] );
                calculatedResults.push_back(nextResult);
            }
            
            // testing the output
            for (int q=0; q<conductor.hack.m_visibleLandmarks_size; q++)
            {
                if (cv::norm( conductor.hack.responseMaps[q] - calculatedResults[q] ) > 0.0001) throw std::runtime_error("conductor.hack.responseMaps[q] - calculatedResults[q] ) < 0.0001 failed");
            }
            
            // releasing the outputs
            carp::freeResponseMaps( &responseMaps, conductor.hack.m_visibleLandmarks_size );

        }
    }
    
    std::cout << "total elapsed time = " << elapsed_time / 1000000. << " s." << std::endl;
    std::cout << std::setprecision(2) << std::fixed;
    std::cout << "processing speed   = " << 1000000. * processed_frames / elapsed_time << "fps" << std::endl;

    //conductor.importer >> BOOST_SERIALIZATION_NVP(conductor.hack);

    return EXIT_SUCCESS;
} // int main


// LuM end of file
