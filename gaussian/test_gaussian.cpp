#include "utility.hpp"
#include "gaussian.pencil.h"

#include <opencv2/core/core.hpp>
#include <opencv2/ocl/ocl.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <prl.h>
#include <chrono>

void time_gaussian( const std::vector<carp::record_t>& pool, const std::vector<int>& sizes )
{
    carp::Timing timing("gaussian blur");

    for ( auto & size : sizes ) {
        cv::Size ksize(size, size+4);

        double gaussX = 7.;
        double gaussY = 9.;

        for ( auto & item : pool ) {
            cv::Mat cpu_gray;

            cv::cvtColor( item.cpuimg(), cpu_gray, CV_RGB2GRAY );
            cpu_gray.convertTo( cpu_gray, CV_32F, 1.0/255. );

            cv::Mat cpu_result, gpu_result, pen_result;
            std::chrono::duration<double> elapsed_time_cpu, elapsed_time_gpu_p_copy, elapsed_time_gpu_nocopy;

            {
                const auto cpu_start = std::chrono::high_resolution_clock::now();
                cv::GaussianBlur( cpu_gray, cpu_result, ksize, gaussX, gaussY, cv::BORDER_REPLICATE );
                const auto cpu_end = std::chrono::high_resolution_clock::now();
                elapsed_time_cpu = cpu_end - cpu_start;
                //Free up resources
            }
            {
                const auto gpu_copy_start = std::chrono::high_resolution_clock::now();
                cv::ocl::oclMat src(cpu_gray);
                cv::ocl::oclMat dst;
                const auto gpu_start = std::chrono::high_resolution_clock::now();
                cv::ocl::GaussianBlur( src, dst, ksize, gaussX, gaussY, cv::BORDER_REPLICATE );
                const auto gpu_end = std::chrono::high_resolution_clock::now();
                gpu_result = dst;
                const auto gpu_copy_end = std::chrono::high_resolution_clock::now();
                elapsed_time_gpu_p_copy = gpu_copy_end - gpu_copy_start;
                elapsed_time_gpu_nocopy = gpu_end      - gpu_start;
                //Free up resources
            }
            {
                cv::Mat kernel_x = cv::getGaussianKernel(ksize.width , gaussX, CV_32F);
                cv::Mat kernel_y = cv::getGaussianKernel(ksize.height, gaussY, CV_32F);

                pen_result.create( cpu_gray.size(), CV_32F );

                prl_timings_start();
                pencil_gaussian( cpu_gray.rows, cpu_gray.cols, cpu_gray.step1(), cpu_gray.ptr<float>()
                               , kernel_x.rows, kernel_x.ptr<float>()
                               , kernel_y.rows, kernel_y.ptr<float>()
                               , pen_result.ptr<float>()
                               );
                prl_timings_stop();
                //Free up resources
            }
            // Verifying the results
            if ( (cv::norm(cpu_result - gpu_result) > 0.01) || (cv::norm(cpu_result - pen_result) > 0.01) ) {
                std::cerr << "ERROR: Results don't match. Writing calculated images." << std::endl;
                std::cerr << "CPU norm:" << cv::norm(cpu_result) << std::endl;
                std::cerr << "GPU norm:" << cv::norm(gpu_result) << std::endl;
                std::cerr << "PEN norm:" << cv::norm(pen_result) << std::endl;
                std::cerr << "GPU-CPU norm:" << cv::norm(gpu_result, cpu_result) << std::endl;
                std::cerr << "PEN-CPU norm:" << cv::norm(pen_result, cpu_result) << std::endl;

                cv::imwrite( "gaussian_cpu.png", cpu_result );
                cv::imwrite( "gaussian_gpu.png", gpu_result );
                cv::imwrite( "gaussian_pen.png", pen_result );
                cv::imwrite( "gaussian_cpugpu.png", cv::abs(cpu_result-gpu_result) );
                cv::imwrite( "gaussian_cpupen.png", cv::abs(cpu_result-pen_result) );
                throw std::runtime_error("The OpenCL or PENCIL results are not equivalent with the C++ results.");
            }

            // Dump execution times for OpenCV calls.
            timing.print( elapsed_time_cpu, elapsed_time_gpu_p_copy, elapsed_time_gpu_nocopy );
        }
    }
    // Dump execution times for PENCIL code.
    prl_timings_dump();
}

int main(int argc, char* argv[])
{
    prl_init(PRL_TARGET_DEVICE_DYNAMIC);
    try {
        std::cout << "This executable is iterating over all the files which are present in the directory `./pool'. " << std::endl;

        auto pool = carp::get_pool("pool");
#ifdef RUN_ONLY_ONE_EXPERIMENT
        time_gaussian( pool, {25} );
#else
        time_gaussian( pool, {5, 15, 25, 35, 45} );
#endif
        prl_shutdown();
        return EXIT_SUCCESS;
    } catch(const std::exception& e) {
        std::cout << e.what() << std::endl;

        prl_shutdown();
        return EXIT_FAILURE;
    }
}
