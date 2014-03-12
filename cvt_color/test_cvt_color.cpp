// UjoImro, 2013

#include <chrono>
#include <opencv2/opencv.hpp>
#include <opencv2/ocl/ocl.hpp>

#include "opencl.hpp"
#include "utility.hpp"
#include "cvt_color.clh"
#include "cvt_color.pencil.h"

template<class T0>
void
time_cvtColor( T0 & pool, size_t iterations)
{
    int bidx = 2;
    carp::TimingLong timing;
    for ( auto & record : pool ) {
        PRINT(record.path());

        std::chrono::microseconds elapsed_time_gpu(0), elapsed_time_cpu(0), elapsed_time_pencil(0);

        for(size_t i = 0; i < iterations; ++i) {
            PRINT(i);
            cv::Mat cpuimg = record.cpuimg();

            cv::Mat cpu_result, gpu_result, pen_result;

            {
                auto start = std::chrono::high_resolution_clock::now();
                cv::cvtColor( cpuimg, cpu_result, CV_RGB2GRAY );
                auto end = std::chrono::high_resolution_clock::now();
                elapsed_time_cpu += (end - start);
            }
            {
                cv::ocl::Context * context = cv::ocl::Context::getContext();
                carp::opencl::device device(context);
                device.source_compile( cvt_color_cl, cvt_color_cl_len, carp::make_vector<std::string>("RGB2Gray") );

                auto start = std::chrono::high_resolution_clock::now();
                cv::ocl::oclMat gpu_gray;
                cv::ocl::oclMat gpuimg(cpuimg);
                gpu_gray.create( gpuimg.rows, gpuimg.cols, CV_8U );

                device["RGB2Gray"]( gpuimg.cols, gpuimg.rows, static_cast<int>(gpuimg.step), static_cast<int>(gpu_gray.step)
                                  , gpuimg.channels()+1, bidx, reinterpret_cast<cl_mem>(gpuimg.data), reinterpret_cast<cl_mem>(gpu_gray.data)
                                  ).groupsize( carp::make_vector<size_t>(16,16), carp::make_vector<size_t>(cpuimg.cols,cpuimg.rows) );
                gpu_result = gpu_gray;
                auto end = std::chrono::high_resolution_clock::now();
                elapsed_time_gpu += (end - start);
            }
            {
                pen_result.create( cpu_result.rows, cpu_result.cols, CV_8U );

                const auto pencil_start = std::chrono::high_resolution_clock::now();
                pencil_RGB2Gray( cpuimg.rows, cpuimg.cols, cpuimg.step1()/cpuimg.channels(), pen_result.step1()
                               , cpuimg.channels(), bidx, cpuimg.data, pen_result.data
                               );
                const auto pencil_end = std::chrono::high_resolution_clock::now();
                elapsed_time_pencil += (pencil_end - pencil_start);
            }

            // Verifying the results
            float pencil_err = cv::norm(pen_result - cpu_result);
            if ( (cv::norm(gpu_result - cpu_result) > 0.01) || (pencil_err>0.01) ) {
                PRINT(cv::norm(gpu_result - cpu_result));
                PRINT(cv::norm(pen_result - cpu_result));
                cv::imwrite( "gpu_img.png", gpu_result );
                cv::imwrite( "cpu_img.png", cpu_result );
                throw std::runtime_error("The GPU results are not equivalent with the CPU results.");
            }
        }
        timing.print("cvtColor", elapsed_time_cpu, elapsed_time_gpu, elapsed_time_pencil);
    }
    return;
}


int main(int argc, char* argv[])
{

#ifndef BENCHMARK_PRINT_GPU_PENCIL_SPEEDUP_ONLY
    std::cout << "This executable is iterating over all the files which are present in the directory `./pool'. " << std::endl;
#endif

    auto pool = carp::get_pool("pool");
    size_t num_iterations = 10;
    time_cvtColor( pool, num_iterations );

    return EXIT_SUCCESS;
}
