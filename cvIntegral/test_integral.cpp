// UjoImro, 2013

#include <chrono>
#include <opencv2/opencv.hpp>
#include <opencv2/ocl/ocl.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "opencl.hpp"
#include "utility.hpp"
#include "imgproc_integral_sum.clh"


template<class T0>
void
time_integral( T0 & pool, int iteration )
{
    carp::TimingLong timing;

    for ( int q=0; q<iteration; q++ ) {

        for ( auto & item : pool ) {
            PRINT(item.path());
            cv::Mat cpu_gray;
            cv::cvtColor( item.cpuimg(), cpu_gray, CV_RGB2GRAY );
            if(cpu_gray.cols * cpu_gray.rows <= 2901 * 2901)
                cv::resize(cpu_gray, cpu_gray, cv::Size(2901,2901));

            cv::Mat cpu_result, gpu_result;
            std::chrono::microseconds elapsed_time_cpu(0), elapsed_time_gpu(0), elapsed_time_pen(99999);

            {
                auto cpu_start = std::chrono::high_resolution_clock::now();
                cv::integral( cpu_gray, cpu_result, CV_32SC1 );
                auto cpu_end = std::chrono::high_resolution_clock::now();
                elapsed_time_cpu += (cpu_end - cpu_start);
            }
            {
                cv::ocl::Context * context = cv::ocl::Context::getContext();
                carp::opencl::device device(context);
                device.source_compile( imgproc_integral_sum_cl, imgproc_integral_sum_cl_len
                                     , carp::make_vector<std::string>("integral_sum_cols_D4", "integral_sum_rows_D4" )
                                     , " -D DOUBLE_SUPPORT"
                                     );

                auto gpu_start = std::chrono::high_resolution_clock::now();

                cv::ocl::oclMat src(cpu_gray);

                int vlen = 4;
                int offset = src.offset / vlen;
                int pre_invalid = src.offset % vlen;
                int vcols = (pre_invalid + src.cols + vlen - 1) / vlen;

                cv::ocl::oclMat t_sum(src.cols, src.rows, CV_32SC1);
                cv::ocl::oclMat sum(src.rows + 1, src.cols + 1, CV_32SC1);
                int sum_offset = sum.offset / vlen;
                device["integral_sum_cols_D4"]( reinterpret_cast<cl_mem>(src.data), reinterpret_cast<cl_mem>(t_sum.data)
                                              , offset, pre_invalid, src.rows, src.cols, static_cast<int>(src.step), static_cast<int>(t_sum.step)
                                              ).groupsize( {256, 1, 1}, {((vcols + 1) / 2) * 256, 1, 1} );
                device["integral_sum_rows_D4"]( reinterpret_cast<cl_mem>(t_sum.data), reinterpret_cast<cl_mem>(sum.data)
                                              , t_sum.rows, t_sum.cols, static_cast<int>(t_sum.step), static_cast<int>(sum.step), sum_offset
                                              ).groupsize( {256, 1, 1}, {t_sum.cols  * 32, 1, 1} );
                gpu_result = sum;
                auto gpu_end = std::chrono::high_resolution_clock::now();
                elapsed_time_gpu += (gpu_end - gpu_start);
            }

            // Verifying the results
            if ( cv::norm(cpu_result - gpu_result) > 0.01 ) {
                PRINT(cv::norm(gpu_result - cpu_result));
                // no use to write out the results, as they are in float
                throw std::runtime_error("The GPU results are not equivalent with the CPU results.");
            }

            timing.print( "integral image", elapsed_time_cpu, elapsed_time_gpu );
        }
    }
}

int main(int argc, char* argv[])
{

    std::cout << "This executable is iterating over all the files which are present in the directory `./pool'. " << std::endl;

    auto pool = carp::get_pool("pool");

    // Initializing OpenCL
    time_integral( pool, 10 );
    return EXIT_SUCCESS;
}
