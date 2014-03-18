// Experimental Research Code for the CARP Project
// UjoImro, 2013

#include <chrono>
#include <opencv2/opencv.hpp>
#include <opencv2/ocl/ocl.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "opencl.hpp"
#include "utility.hpp"
#include "imgproc_convolve.clh"
#include "filter2D.pencil.h"

namespace {
    inline int divUp(int total, int grain)
    {
        return (total + grain - 1) / grain;
    }
} // unnamed namespace

namespace carp {

    cv::ocl::oclMat
    filter2D( carp::opencl::device & device, cv::ocl::oclMat gpu_gray, cv::ocl::oclMat kernel_gpu, int border_type ) {
        CV_Assert(gpu_gray.depth() == CV_32F);
        CV_Assert(kernel_gpu.depth() == CV_32F);
        cv::ocl::oclMat gpu_convolve;
        gpu_convolve.create(gpu_gray.size(), gpu_gray.type());
        CV_Assert(gpu_gray.type() == kernel_gpu.type() && gpu_gray.size() == gpu_convolve.size());

        std::string kernelName = "convolve_D5";

        CV_Assert(gpu_gray.depth() == CV_32FC1);
        CV_Assert(kernel_gpu.depth() == CV_32F);
        CV_Assert(kernel_gpu.cols <= 17 && kernel_gpu.rows <= 17);

        gpu_convolve.create(gpu_gray.size(), gpu_gray.type());

        CV_Assert(gpu_gray.cols == gpu_convolve.cols && gpu_gray.rows == gpu_convolve.rows);
        CV_Assert(gpu_gray.type() == gpu_convolve.type());

        int channels = gpu_convolve.oclchannels();

        size_t vector_length = 1;
        int offset_cols = ((gpu_convolve.offset % gpu_convolve.step) / gpu_convolve.elemSize1()) & (vector_length - 1);
        int cols = divUp(gpu_convolve.cols * channels + offset_cols, vector_length);
        int rows = gpu_convolve.rows;

        std::vector<size_t> localThreads  = { 16, 16, 1 };
        std::vector<size_t> globalThreads = { divUp(cols, localThreads[0]) *localThreads[0],
                                              divUp(rows, localThreads[1]) *localThreads[1],
                                              1 };

        device[kernelName] (
            reinterpret_cast<cl_mem>(gpu_gray.data),
            reinterpret_cast<cl_mem>(kernel_gpu.data),
            reinterpret_cast<cl_mem>(gpu_convolve.data),
            gpu_gray.rows,
            cols,
            static_cast<int>(gpu_gray.step),
            static_cast<int>(gpu_convolve.step),
            static_cast<int>(kernel_gpu.step),
            kernel_gpu.rows,
            kernel_gpu.cols
            ).groupsize( localThreads, globalThreads );

        return gpu_convolve;
    } // filter2D

} // namespace carp

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
            cv::Mat kernel_cpu(3, 3, CV_32F, kernel_data);

            cv::Mat cpu_result, gpu_result, pen_result;
            std::chrono::microseconds elapsed_time_cpu, elapsed_time_gpu_p_copy, elapsed_time_gpu_nocopy, elapsed_time_pencil;

            {
                const auto cpu_start = std::chrono::high_resolution_clock::now();
                cv::filter2D( cpu_gray, cpu_result, -1, kernel_cpu, cv::Point(-1,-1), 0.0, cv::BORDER_REPLICATE );
                const auto cpu_end = std::chrono::high_resolution_clock::now();
                elapsed_time_cpu = cpu_end - cpu_start;
            }
            {
                cv::ocl::Context * context = cv::ocl::Context::getContext();
                carp::opencl::device device(context);
                device.source_compile( imgproc_convolve_cl, imgproc_convolve_cl_len, { "convolve_D5" } );
                const auto gpu_start_copy = std::chrono::high_resolution_clock::now();
                cv::ocl::oclMat gpu_gray(cpu_gray);
                cv::ocl::oclMat kernel_gpu(kernel_cpu);
                const auto gpu_start = std::chrono::high_resolution_clock::now();
                cv::ocl::oclMat gpu_convolve = carp::filter2D( device, gpu_gray, kernel_gpu, cv::BORDER_REPLICATE );
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
                                 kernel_cpu.rows, kernel_cpu.cols, kernel_cpu.step1(), kernel_cpu.ptr<float>(),
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

#ifndef BENCHMARK_PRINT_GPU_PENCIL_SPEEDUP_ONLY
    std::cout << "This executable is iterating over all the files which are present in the directory `./pool'. " << std::endl;
#endif

    auto pool = carp::get_pool("pool");
    time_filter2D( pool, 3 );
    return EXIT_SUCCESS;
} // main
