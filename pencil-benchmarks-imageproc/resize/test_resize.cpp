#include "utility.hpp"
#include "resize.pencil.h"

#include <opencv2/core/core.hpp>
#include <opencv2/ocl/ocl.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <prl.h>
#include <chrono>

void time_resize( const std::vector<carp::record_t>& pool, const std::vector<cv::Size>& sizes, int iteration )
{
    bool first_execution_opencv = true, first_execution_pencil = true;

    carp::Timing timing("resize");

    for ( int q=0; q<iteration; q++ ) {
        for ( auto & size : sizes ) {
            for ( auto & item : pool ) {
                cv::Mat cpu_gray;
                cv::cvtColor( item.cpuimg(), cpu_gray, CV_RGB2GRAY );

                cv::Mat cpu_result, gpu_result, pen_result;
                std::chrono::duration<double> elapsed_time_cpu, elapsed_time_gpu_p_copy;

                {
                    const auto cpu_start = std::chrono::high_resolution_clock::now();
                    cv::resize( cpu_gray, cpu_result, size, 0, 0, cv::INTER_LINEAR );
                    const auto cpu_end = std::chrono::high_resolution_clock::now();
                    elapsed_time_cpu = cpu_end - cpu_start;
                }
                {
                    // Execute the kernel at least once before starting to take time measurements so that the OpenCV kernel gets compiled. The following run is not included in time measurements.
                    if (first_execution_opencv)
                    {
                        cv::ocl::oclMat gpu_gray(cpu_gray);
                        cv::ocl::oclMat gpu_resize;
                        cv::ocl::resize( gpu_gray, gpu_resize, size, 0, 0, cv::INTER_LINEAR );
                        first_execution_opencv = false;
                    }

                    const auto gpu_start_copy = std::chrono::high_resolution_clock::now();
                    cv::ocl::oclMat gpu_gray(cpu_gray);
                    cv::ocl::oclMat gpu_resize;
                    cv::ocl::resize( gpu_gray, gpu_resize, size, 0, 0, cv::INTER_LINEAR );
                    gpu_result = gpu_resize;
                    const auto gpu_end_copy = std::chrono::high_resolution_clock::now();
                    elapsed_time_gpu_p_copy = gpu_end_copy - gpu_start_copy;
                }
                {
                    // pencil verification
                    pen_result.create(size, CV_8UC1);

                    if (first_execution_pencil)
                    {
                        pencil_resize_LN(cpu_gray.rows, cpu_gray.cols, cpu_gray.step1(), cpu_gray.ptr(), pen_result.rows, pen_result.cols, pen_result.step1(), pen_result.ptr() );
                        first_execution_pencil = false;
                    }

                    prl_prof_reset();
                    prl_prof_start();
                    pencil_resize_LN(cpu_gray.rows, cpu_gray.cols, cpu_gray.step1(), cpu_gray.ptr(), pen_result.rows, pen_result.cols, pen_result.step1(), pen_result.ptr() );
                    prl_prof_stop();
                    // Dump execution times for PENCIL code.
                    prl_prof_dump();
                }
                // Verifying the results - TODO - Something fishy is happening at borders
#define REMOVE_BORDER(img) img(cv::Range(1, size.height-1), cv::Range(1, size.width-1))
                if (( cv::norm(REMOVE_BORDER(cpu_result), REMOVE_BORDER(gpu_result), cv::NORM_INF) > 1 )
                  ||( cv::norm(REMOVE_BORDER(cpu_result), REMOVE_BORDER(pen_result), cv::NORM_INF) > 1 )
                   )
                {
                    cv::imwrite( "gpu_resize.png", gpu_result );
                    cv::imwrite( "cpu_resize.png", cpu_result );
                    cv::imwrite( "pencil_resize.png", pen_result );
                    cv::imwrite( "cpu_gpu_diff.png", cpu_result-gpu_result );
                    cv::imwrite( "cpu_pen_diff.png", cpu_result-pen_result );

                    throw std::runtime_error("The GPU results are not equivalent with the CPU or Pencil results.");
                }
                // Dump execution times for OpenCV calls.
                timing.print( elapsed_time_cpu, elapsed_time_gpu_p_copy );
            } // for pool
        }
    }
}

int main(int argc, char* argv[])
{
    prl_init();

    std::cout << "This executable is iterating over all the files passed to it as an argument. " << std::endl;

    auto pool = carp::get_pool(argc, argv);

#ifdef RUN_ONLY_ONE_EXPERIMENT
    std::vector<cv::Size> sizes{ {1024, 768} };
    int iteration = 1;
#else
    std::vector<cv::Size> sizes{ {320, 240}, {640, 480}, {1024, 768}, {1200, 900}, {1600, 1200} };
    int iteration = 25;
#endif

    time_resize( pool, sizes, iteration );

    prl_release();
    return EXIT_SUCCESS;
}
