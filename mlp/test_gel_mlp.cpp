
#include "serialization.hpp"

// UjoImro, 2013
// Experimental code for the CARP Project
// Copyright (c) RealEyes, 2013
// This version tests the responseMap calculation with input dumps

#include "opencl.hpp"
#include "memory.hpp"
#include "bench_mlp.hpp"

#include <tbb/parallel_for.h>
#include <tbb/concurrent_vector.h>
#include <tbb/tbb_thread.h>
#include <tbb/task_scheduler_init.h>


const int processed_frames = 100;

int main()
{
    carp::conductor_t conductor;
    std::chrono::duration<double> elapsed_time(0);

    int numberOfThreads =  tbb::tbb_thread::hardware_concurrency();
    tbb::task_scheduler_init init(numberOfThreads);

    for ( conductor.importer >> BOOST_SERIALIZATION_NVP(conductor.id)
        ; (conductor.id != -1) && (conductor.id != processed_frames)
        ; conductor.importer >> BOOST_SERIALIZATION_NVP(conductor.id)
        )
    {
        carp::hack_t package;
        conductor.importer >> BOOST_SERIALIZATION_NVP(package);

        // here comes the function call
        {
            gel::MLP<double> mlp;
            tbb::concurrent_vector< cv::Mat_<double> > calculatedResults(package.m_visibleLandmarks_size);
            auto start = std::chrono::high_resolution_clock::now();
            tbb::parallel_for(0,package.m_visibleLandmarks_size,[&](int q) {
                const cv::Point2i center(cvRound(package.shape(2*q,0)),cvRound(package.shape(2*q+1,0)));
                calculatedResults[q] = package.m_classifiers[q].generateResponseMap(package.alignedImage, center, package.m_mapSize);
            });
            auto end = std::chrono::high_resolution_clock::now();
            elapsed_time += (end - start);

            // testing the output
            for (int q=0; q<package.m_visibleLandmarks_size; q++)
            {
                if (cv::norm( package.responseMaps[q] - calculatedResults[q] ) > 0.0001) throw std::runtime_error("package.responseMaps[q] - calculatedResults[q] ) < 0.0001 failed");
            }
        }
    }

    std::cout << "total elapsed time = " << elapsed_time.count() << " s." << std::endl;
    std::cout << std::setprecision(2) << std::fixed;
    std::cout << "processing speed   = " << processed_frames / elapsed_time.count() << "fps" << std::endl;

    return EXIT_SUCCESS;
}
