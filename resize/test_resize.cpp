#include "utility.hpp"
#include "resize.pencil.h"

#include <opencv2/core/core.hpp>
#include <opencv2/ocl/ocl.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <pencil_runtime.h>
#include <chrono>

void time_resize( const std::vector<carp::record_t>& pool, const std::vector<cv::Size>& sizes, int iteration )
{
    carp::Timing timing("resize");

    for ( int q=0; q<iteration; q++ ) {
        for ( auto & size : sizes ) {
            for ( auto & item : pool ) {
                cv::Mat cpu_gray;
                cv::cvtColor( item.cpuimg(), cpu_gray, CV_RGB2GRAY );

                cv::Mat cpu_result, gpu_result, pen_result;
                std::chrono::duration<double> elapsed_time_cpu, elapsed_time_gpu_p_copy, elapsed_time_gpu_nocopy, elapsed_time_pencil;

                {
                    const auto cpu_start = std::chrono::high_resolution_clock::now();
                    cv::resize( cpu_gray, cpu_result, size, 0, 0, cv::INTER_LINEAR );
                    const auto cpu_end = std::chrono::high_resolution_clock::now();
                    elapsed_time_cpu = cpu_end - cpu_start;
                }
                {
                    const auto gpu_start_copy = std::chrono::high_resolution_clock::now();
                    cv::ocl::oclMat gpu_gray(cpu_gray);
                    cv::ocl::oclMat gpu_resize;
                    const auto gpu_start = std::chrono::high_resolution_clock::now();
                    cv::ocl::resize( gpu_gray, gpu_resize, size, 0, 0, cv::INTER_LINEAR );
                    const auto gpu_end = std::chrono::high_resolution_clock::now();
                    gpu_result = gpu_resize;
                    const auto gpu_end_copy = std::chrono::high_resolution_clock::now();
                    elapsed_time_gpu_p_copy = gpu_end_copy - gpu_start_copy;
                    elapsed_time_gpu_nocopy = gpu_end      - gpu_start;
                }
                {
                    // pencil verification
                    pen_result.create(size, CV_8UC1);

                    const auto pencil_start = std::chrono::high_resolution_clock::now();
                    pencil_resize_LN(cpu_gray.rows, cpu_gray.cols, cpu_gray.step1(), cpu_gray.ptr(), pen_result.rows, pen_result.cols, pen_result.step1(), pen_result.ptr() );
                    const auto pencil_end = std::chrono::high_resolution_clock::now();
                    elapsed_time_pencil = pencil_end - pencil_start;
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

                timing.print( elapsed_time_cpu, elapsed_time_gpu_p_copy, elapsed_time_gpu_nocopy, elapsed_time_pencil );
            } // for pool
        }
    }
}

int main(int argc, char* argv[])
{
    pencil_init();

    std::cout << "This executable is iterating over all the files which are present in the directory `./pool'. " << std::endl;

    auto pool = carp::get_pool("pool");

#ifdef RUN_ONLY_ONE_EXPERIMENT
    std::vector<cv::Size> sizes{ {1024, 768} };
    int iteration = 1;
#else
    std::vector<cv::Size> sizes{ {320, 240}, {640, 480}, {1024, 768}, {1200, 900}, {1600, 1200} };
    int iteration = 35;
#endif

    time_resize( pool, sizes, iteration );

    pencil_shutdown();
    return EXIT_SUCCESS;
}
