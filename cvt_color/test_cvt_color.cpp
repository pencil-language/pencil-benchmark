// UjoImro, 2013

#include <chrono>
#include <opencv2/opencv.hpp>
#include <opencv2/ocl/ocl.hpp>

#include "opencl.hpp"
#include "utility.hpp"
#include "cvt_color.clh"
#include "cvt_color.pencil.h"

void time_cvtColor( const std::vector<carp::record_t>& pool, size_t iterations)
{
    int bidx = 2;
    carp::Timing timing("cvt_color");
    for ( auto & record : pool ) {
        PRINT(record.path());
        for(size_t i = 0; i < iterations; ++i) {
            PRINT(i);

            cv::Mat cpuimg = record.cpuimg();
            cv::Mat cpu_result, gpu_result, pen_result;

            std::chrono::duration<double> elapsed_time_cpu, elapsed_time_gpu_p_copy, elapsed_time_gpu_nocopy, elapsed_time_pencil;
            {
                const auto start = std::chrono::high_resolution_clock::now();
                cv::cvtColor( cpuimg, cpu_result, CV_RGB2GRAY );
                const auto end = std::chrono::high_resolution_clock::now();
                elapsed_time_cpu = end - start;
            }
            {
                cv::ocl::Context * context = cv::ocl::Context::getContext();
                carp::opencl::device device(context);
                device.source_compile( cvt_color_cl, cvt_color_cl_len, {"RGB2Gray"} );

                const auto start_copy = std::chrono::high_resolution_clock::now();
                cv::ocl::oclMat gpu_gray;
                cv::ocl::oclMat gpuimg(cpuimg);
                gpu_gray.create( gpuimg.rows, gpuimg.cols, CV_8U );

                const auto start = std::chrono::high_resolution_clock::now();
                device["RGB2Gray"]( gpuimg.cols, gpuimg.rows, static_cast<int>(gpuimg.step), static_cast<int>(gpu_gray.step)
                                  , gpuimg.channels()+1, bidx, reinterpret_cast<cl_mem>(gpuimg.data), reinterpret_cast<cl_mem>(gpu_gray.data)
                                  ).groupsize( {16,16}, {cpuimg.cols, cpuimg.rows});
                const auto end = std::chrono::high_resolution_clock::now();
                gpu_result = gpu_gray;
                const auto end_copy = std::chrono::high_resolution_clock::now();
                elapsed_time_gpu_p_copy = end_copy - start_copy;
                elapsed_time_gpu_nocopy = end      - start;
            }
            {
                pen_result.create( cpu_result.rows, cpu_result.cols, CV_8U );

                const auto pencil_start = std::chrono::high_resolution_clock::now();
                pencil_RGB2Gray( cpuimg.rows, cpuimg.cols, cpuimg.step1()/cpuimg.channels(), pen_result.step1()
                               , cpuimg.channels(), bidx, cpuimg.data, pen_result.data
                               );
                const auto pencil_end = std::chrono::high_resolution_clock::now();
                elapsed_time_pencil = pencil_end - pencil_start;
            }

            // Verifying the results
            float pencil_err = cv::norm(pen_result - cpu_result);
            if ( (cv::norm(gpu_result - cpu_result) > 0.01) || (pencil_err>0.01) ) {
                PRINT(cv::norm(gpu_result - cpu_result));
                PRINT(cv::norm(pen_result - cpu_result));
                cv::imwrite( "gpu_img.png", gpu_result );
                cv::imwrite( "cpu_img.png", cpu_result );
		cv::imwrite("pencil_dilate.png", pen_result );
                throw std::runtime_error("The GPU results are not equivalent with the CPU results.");
            }
            timing.print(elapsed_time_cpu, elapsed_time_gpu_p_copy, elapsed_time_gpu_nocopy, elapsed_time_pencil);
        }
    }
}

int main(int argc, char* argv[])
{

    std::cout << "This executable is iterating over all the files which are present in the directory `./pool'. " << std::endl;

    auto pool = carp::get_pool("pool");
    size_t num_iterations = 100;

#ifdef RUN_ONLY_ONE_EXPERIMENT
    num_iterations = 1;
#endif

    time_cvtColor( pool, num_iterations );

    return EXIT_SUCCESS;
}
