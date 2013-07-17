// UjoImro, 2013
// Experimental code for the CARP Project
// Copyright (c) RealEyes, 2013
// This version tests the responseMap calculation with input dumps

#include <boost/smart_ptr.hpp>

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
    conductor_t conductor;
    int fail = 0;
    long int elapsed_time = 0;
    int64_t maxnetallocated = 0;
    int64_t maxgrossallocated = 0;
    
    int numberOfThreads =  tbb::tbb_thread::hardware_concurrency();
    PRINT(numberOfThreads);    
    tbb::task_scheduler_init init(numberOfThreads);
    //tbb::task_scheduler_init init(1);

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
            gel::MLP<double> mlp;
            tbb::concurrent_vector< cv::Mat_<double> > calculatedResults(conductor.hack.m_visibleLandmarks_size);
            auto start = std::chrono::high_resolution_clock::now();
            tbb::parallel_for(0,conductor.hack.m_visibleLandmarks_size,[&](int q){
//                    for(int q = 0; q < conductor.hack.m_visibleLandmarks_size; ++q) {
                        const cv::Point2i center(cvRound(conductor.hack.shape(2*q,0)),cvRound(conductor.hack.shape(2*q+1,0)));
                        calculatedResults[q] = conductor.hack.m_classifiers[q].generateResponseMap(conductor.hack.alignedImage, center, conductor.hack.m_mapSize);
                }
                    ); // tbb::parallel_for
            auto end = std::chrono::high_resolution_clock::now();
            elapsed_time += carp::microseconds(end - start);

            // testing the output
            for (int q=0; q<conductor.hack.m_visibleLandmarks_size; q++)
            {
                //PRINT(cv::norm( conductor.hack.responseMaps[q] - calculatedResults[q] ));
                if (cv::norm( conductor.hack.responseMaps[q] - calculatedResults[q] ) > 0.0001) throw std::runtime_error("conductor.hack.responseMaps[q] - calculatedResults[q] ) < 0.0001 failed");
            }
        }
    }
    
    std::cout << "total elapsed time = " << elapsed_time / 1000000. << " s." << std::endl;
    std::cout << std::setprecision(2) << std::fixed;
    std::cout << "processing speed   = " << 1000000. * processed_frames / elapsed_time << "fps" << std::endl;

    //conductor.importer >> BOOST_SERIALIZATION_NVP(conductor.hack);

    PRINT(maxnetallocated);
    PRINT(maxgrossallocated);
    
    return EXIT_SUCCESS;
} // int main


// LuM end of file
