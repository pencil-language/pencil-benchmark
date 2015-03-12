#include "utility.hpp"
#include "histogram.pencil.h"

#include <opencv2/core/core.hpp>
#include <opencv2/ocl/ocl.hpp>

#include <prl.h>
#include <chrono>

void time_histogram( const std::vector<carp::record_t>& pool, size_t iterations)
{
    bool first_execution_opencv = true, first_execution_pencil = true;

    carp::Timing timing("histogram");
    for ( auto & item : pool ) {
        for(size_t i = 0; i < iterations; ++i) {
            cv::Mat cpuimg;

            cv::cvtColor( item.cpuimg(), cpuimg, CV_RGB2GRAY );

            cv::Mat cpu_result, gpu_result, pen_result;

            std::chrono::duration<double> elapsed_time_cpu, elapsed_time_gpu_p_copy, elapsed_time_gpu_nocopy;
            {
                cv::Mat tmp;
                const auto start = std::chrono::high_resolution_clock::now();
                const int channels = 0;
                const int histSize = 256;
                const float range[] = {0, 256};
                const float* ranges[] = {range};
                cv::calcHist( &cpuimg, 1, &channels, cv::Mat(), tmp, 1, &histSize, ranges);
                const auto end = std::chrono::high_resolution_clock::now();
                tmp = tmp.t();
                tmp.convertTo(cpu_result, CV_32S);
                elapsed_time_cpu = end - start;
            }
            {
                // Execute the kernel at least once before starting to take time measurements so that the OpenCV kernel gets compiled. The following run is not included in time measurements.
                if (first_execution_opencv)
                {
                    cv::ocl::oclMat gpuimg(cpuimg);
                    cv::ocl::oclMat result;
                    cv::ocl::calcHist(gpuimg, result);
                    first_execution_opencv = false;
                }

                const auto start_copy = std::chrono::high_resolution_clock::now();
                cv::ocl::oclMat gpuimg(cpuimg);
                cv::ocl::oclMat result;

                const auto start = std::chrono::high_resolution_clock::now();
                cv::ocl::calcHist(gpuimg, result);
                const auto end = std::chrono::high_resolution_clock::now();
                gpu_result = result;
                const auto end_copy = std::chrono::high_resolution_clock::now();
                elapsed_time_gpu_p_copy = end_copy - start_copy;
                elapsed_time_gpu_nocopy = end      - start;
            }
            {
                pen_result.create( cpu_result.rows, cpu_result.cols, CV_32S );

                if (first_execution_pencil)
                {
                    pencil_calcHist( cpuimg.rows, cpuimg.cols, cpuimg.step1(), cpuimg.ptr<uint8_t>(), pen_result.ptr<int>());
                    first_execution_pencil = false;
                }

                prl_timings_start();
                pencil_calcHist( cpuimg.rows, cpuimg.cols, cpuimg.step1(), cpuimg.ptr<uint8_t>()
                               , pen_result.ptr<int>()
                               );
                prl_timings_stop();
            }
            // Verifying the results
            float gpu_err = cv::norm(gpu_result - cpu_result);
            float pencil_err = cv::norm(pen_result - cpu_result);
            if ( (gpu_err > 0.01) || (pencil_err>0.01) )
            {
                cv::imwrite( "gpu_img.png", gpu_result );
                cv::imwrite( "cpu_img.png", cpu_result );
                cv::imwrite("pencil_dilate.png", pen_result );
                throw std::runtime_error("The GPU results are not equivalent with the CPU results.");
            }
            // Dump execution times for OpenCV calls.
            timing.print(elapsed_time_cpu, elapsed_time_gpu_p_copy, elapsed_time_gpu_nocopy);
        }
    }
    // Dump execution times for PENCIL code.
    prl_timings_dump();
}

int main(int argc, char* argv[])
{
    prl_init((prl_init_flags)(PRL_TARGET_DEVICE_DYNAMIC | PRL_PROFILING_ENABLED));

    std::cout << "This executable is iterating over all the files which are present in the directory `./pool'. " << std::endl;

    auto pool = carp::get_pool("pool");
    size_t num_iterations = 80;

#ifdef RUN_ONLY_ONE_EXPERIMENT
    num_iterations = 1;
#endif

    time_histogram( pool, num_iterations );

    prl_shutdown();
    return EXIT_SUCCESS;
}
