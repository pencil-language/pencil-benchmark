// Experimental Research Code for the CARP Project
// UjoImro, 2013

#include <chrono>
#include <opencv2/core/core.hpp>
#include <opencv2/ocl/ocl.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "opencl.hpp"
#include "utility.hpp"
#include "filter2D.pencil.h"

void time_filter2D( const std::vector<carp::record_t>& pool, int iteration )
{
    carp::Timing timing("2D filter");

    for ( int q=0; q<iteration; q++ ) {
        for ( auto & item : pool ) {
            PRINT(item.path());

            cv::Mat cpu_gray;
            cv::cvtColor( item.cpuimg(), cpu_gray, CV_RGB2GRAY );
            cpu_gray.convertTo( cpu_gray, CV_32F, 1.0/255. );

            float kernel_data[] = {-1, -1, -1
                                  , 0,  0,  0
                                  , 1,  1,  1
                                  };
            cv::Mat kernel(3, 3, CV_32F, kernel_data);

            cv::Mat cpu_result, gpu_result, pen_result;
            std::chrono::duration<double> elapsed_time_cpu, elapsed_time_gpu_p_copy, elapsed_time_gpu_nocopy, elapsed_time_pencil;

            {
                const auto cpu_start = std::chrono::high_resolution_clock::now();
                cv::filter2D( cpu_gray, cpu_result, -1, kernel, cv::Point(-1,-1), 0.0, cv::BORDER_REPLICATE );
                const auto cpu_end = std::chrono::high_resolution_clock::now();
                elapsed_time_cpu = cpu_end - cpu_start;
            }
            {
                const auto gpu_start_copy = std::chrono::high_resolution_clock::now();
                cv::ocl::oclMat gpu_gray(cpu_gray);
                cv::ocl::oclMat gpu_convolve;
                const auto gpu_start = std::chrono::high_resolution_clock::now();
                cv::ocl::filter2D( gpu_gray, gpu_convolve, -1, kernel, cv::Point(-1, -1), 0.0, cv::BORDER_REPLICATE );
                const auto gpu_end = std::chrono::high_resolution_clock::now();
                gpu_result = gpu_convolve;
                const auto gpu_end_copy = std::chrono::high_resolution_clock::now();
                elapsed_time_gpu_p_copy = gpu_end_copy - gpu_start_copy;
                elapsed_time_gpu_nocopy = gpu_end      - gpu_start;
            }
            {
                // pencil test:
                pen_result = cv::Mat(cpu_gray.size(), CV_32F);
                const auto pencil_start = std::chrono::high_resolution_clock::now();
                pencil_filter2D( cpu_gray.rows, cpu_gray.cols, cpu_gray.step1(), cpu_gray.ptr<float>(),
                                 kernel.rows, kernel.cols, kernel.step1(), kernel.ptr<float>(),
                                 pen_result.ptr<float>() );
                const auto pencil_end = std::chrono::high_resolution_clock::now();
                elapsed_time_pencil = pencil_end - pencil_start;
            }
            // Verifying the results
            if ( (cv::norm(cpu_result - gpu_result) > 0.01) ||
                 (cv::norm(pen_result - cpu_result) > 0.01) )
            {
                PRINT(cv::norm(gpu_result - cpu_result));
                PRINT(cv::norm(pen_result - cpu_result));

                cv::Mat cpu;
                cv::Mat pencil;
                cv::Mat gpu;
                cpu_result.convertTo( cpu, CV_8U, 255. );
                gpu_result.convertTo( gpu, CV_8U, 255. );
                pen_result.convertTo( pencil, CV_8U, 255. );

                PRINT(cv::norm(gpu-cpu));
                PRINT(cv::norm(gpu-pencil));
                PRINT(cv::norm(cpu-pencil));

                cv::imwrite( "host_convolve.png", cpu );
                cv::imwrite( "gpu_convolve.png", gpu );
                cv::imwrite( "pencil_convolve.png", pencil );
                cv::imwrite( "diff_convolve.png", cv::abs(gpu-pencil) );

                throw std::runtime_error("The GPU results are not equivalent with the CPU results.");
            }

            timing.print( elapsed_time_cpu, elapsed_time_gpu_p_copy, elapsed_time_gpu_nocopy, elapsed_time_pencil );
        }
    }
}

int main(int argc, char* argv[])
{

    std::cout << "This executable is iterating over all the files which are present in the directory `./pool'. " << std::endl;

    auto pool = carp::get_pool("pool");

#ifdef RUN_ONLY_ONE_EXPERIMENT
    time_filter2D( pool, 1 );
#else
    time_filter2D( pool, 35 );
#endif

    return EXIT_SUCCESS;
} // main
