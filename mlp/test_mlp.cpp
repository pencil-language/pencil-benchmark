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

#include "mlp.hpp"
#include "mlp_impl_pencil.h"
#include "utility.hpp"
#include "serialization.hpp"

const int processed_frames = 100;

int main()
{
    carp::conductor_t conductor;
    std::chrono::duration<double> elapsed_time(0);
    int i = 1;

    for ( conductor.importer >> BOOST_SERIALIZATION_NVP(conductor.id)
        ; (conductor.id != -1) && (conductor.id != processed_frames)
        ; conductor.importer >> BOOST_SERIALIZATION_NVP(conductor.id)
        )
    {
        carp::hack_t package;
        conductor.importer >> BOOST_SERIALIZATION_NVP(package);

        // here comes the function call
        {
	    std::cout << i++ << "/" << processed_frames << std::endl;
            // preparing the inputs
            MatChar alignedImage = carp::convertCVToMatChar(package.alignedImage);
            MatFloat shape = carp::convertCVToMatFloat(package.shape);
            mlp * m_classifiers = carp::convertHackToMlp(package);
            MatFloat * responseMaps;
            carp::allocateResponseMaps( package.m_mapSize, package.m_visibleLandmarks_size, &responseMaps );

            auto start = std::chrono::high_resolution_clock::now();
            calculateMaps( package.m_visibleLandmarks_size
                         , package.m_mapSize
                         , alignedImage
                         , shape
                         , m_classifiers
                         , &responseMaps
                         );
            auto end = std::chrono::high_resolution_clock::now();
            elapsed_time += (end - start);

            // releasing the inputs
            freeMatChar(&alignedImage);
            freeMatFloat(&shape);
            carp::freeClassifiers(&m_classifiers, package.m_classifiers.size());

            // converting-testing the outputs
            for (int q=0; q<package.m_visibleLandmarks_size; q++)
            {
                cv::Mat_<double> nextResult = carp::convertMatFloatToCV( responseMaps[q] );
                if (cv::norm( package.responseMaps[q] - nextResult ) > 0.0001) throw std::runtime_error("package.responseMaps[q] - calculatedResults[q] ) < 0.0001 failed");
            }

            // releasing the outputs
            carp::freeResponseMaps( &responseMaps, package.m_visibleLandmarks_size );

        }
    }

    std::cout << "total elapsed time = " << elapsed_time.count() << " s." << std::endl;
    std::cout << std::setprecision(2) << std::fixed;
    std::cout << "processing speed   = " << processed_frames / elapsed_time.count() << "fps" << std::endl;

    return EXIT_SUCCESS;
}