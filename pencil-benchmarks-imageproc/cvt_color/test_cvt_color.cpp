#include "utility.hpp"
#include "cvt_color.pencil.h"

#include <opencv2/core/core.hpp>
#include <opencv2/ocl/ocl.hpp>

#include <prl.h>
#include <chrono>

void time_cvtColor( const std::vector<carp::record_t>& pool, size_t iterations)
{
    bool first_execution_opencv = true, first_execution_pencil = true;

    carp::Timing timing("cvt_color");
    for ( auto & record : pool ) {
        for(size_t i = 0; i < iterations; ++i) {
            cv::Mat cpuimg = record.cpuimg();
            cv::Mat cpu_result, gpu_result, pen_result;

            std::chrono::duration<double> elapsed_time_cpu, elapsed_time_gpu_p_copy;
            {
                const auto start = std::chrono::high_resolution_clock::now();
                cv::cvtColor( cpuimg, cpu_result, CV_RGB2GRAY );
                const auto end = std::chrono::high_resolution_clock::now();
                elapsed_time_cpu = end - start;
            }
            {
                // Execute the kernel at least once before starting to take time measurements so that the OpenCV kernel gets compiled. The following run is not included in time measurements.
                if (first_execution_opencv)
                {
                    cv::ocl::oclMat gpuimg(cpuimg);
                    cv::ocl::oclMat gpu_gray;
                    cv::ocl::cvtColor( gpuimg, gpu_gray, CV_RGB2GRAY );
                    first_execution_opencv = false;
                }

                const auto start_copy = std::chrono::high_resolution_clock::now();
                cv::ocl::oclMat gpuimg(cpuimg);
                cv::ocl::oclMat gpu_gray;
                cv::ocl::cvtColor( gpuimg, gpu_gray, CV_RGB2GRAY );
                gpu_result = gpu_gray;
                const auto end_copy = std::chrono::high_resolution_clock::now();
                elapsed_time_gpu_p_copy = end_copy - start_copy;
            }
            {
                pen_result.create( cpu_result.rows, cpu_result.cols, CV_8U );

                if (first_execution_pencil)
                {
                    pencil_RGB2Gray( cpuimg.rows, cpuimg.cols, cpuimg.step1()/cpuimg.channels(), pen_result.step1(), cpuimg.data, pen_result.data);
                    first_execution_pencil = false;
                }

                prl_timings_reset();
                prl_timings_start();
                pencil_RGB2Gray( cpuimg.rows, cpuimg.cols, cpuimg.step1()/cpuimg.channels(), pen_result.step1()
                               , cpuimg.data, pen_result.data
                               );
                prl_timings_stop();
                // Dump execution times for PENCIL code.
                prl_timings_dump();
            }

            // Verifying the results
            float opencl_err = cv::norm(gpu_result - cpu_result);
            float pencil_err = cv::norm(pen_result - cpu_result);
            if ( opencl_err > 0.01 || pencil_err > 0.01 ) {
                cv::imwrite( "cpu_cvtcolor.png", cpu_result );
                cv::imwrite( "gpu_cvtcolor.png", gpu_result );
                cv::imwrite( "pen_cvtcolor.png", pen_result );
                throw std::runtime_error("The GPU results are not equivalent with the CPU results.");
            }
            // Dump execution times for OpenCV calls.
            timing.print(elapsed_time_cpu, elapsed_time_gpu_p_copy);
        }
    }
}

int main(int argc, char* argv[])
{
    prl_init((prl_init_flags)(PRL_TARGET_DEVICE_DYNAMIC | PRL_PROFILING_ENABLED));

    std::cout << "This executable is iterating over all the files passed to it as an argument. " << std::endl;

    auto pool = carp::get_pool(argc, argv);
    size_t num_iterations;

#ifdef RUN_ONLY_ONE_EXPERIMENT
    num_iterations = 1;
#else
    num_iterations = 65;
#endif

    time_cvtColor( pool, num_iterations );

    prl_shutdown();
    return EXIT_SUCCESS;
}
