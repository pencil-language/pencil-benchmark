#include "utility.hpp"
#include "warpAffine.pencil.h"

#include <opencv2/core/core.hpp>
#include <opencv2/ocl/ocl.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <prl.h>
#include <chrono>

namespace
{
    template <class T0>
    void convert_coeffs(T0 * M)
    {
        T0 D = M[0] * M[4] - M[1] * M[3];
        D = D != 0 ? 1. / D : 0;
        T0 A11 = M[4] * D, A22 = M[0] * D;
        M[0] = A11;
        M[1] *= -D;
        M[3] *= -D;
        M[4] = A22;
        T0 b1 = -M[0] * M[2] - M[1] * M[5];
        T0 b2 = -M[3] * M[2] - M[4] * M[5];
        M[2] = b1;
        M[5] = b2;
    }
}

void time_affine( const std::vector<carp::record_t>& pool, int iteration )
{
    carp::Timing timing("affine transform");

    for ( int q=0; q<iteration; q++ ) {
        for ( auto & item : pool ) {
            cv::Mat cpu_gray;
            cv::cvtColor( item.cpuimg(), cpu_gray, CV_RGB2GRAY );
            cpu_gray.convertTo( cpu_gray, CV_32F, 1.0/255. );

            std::vector<float> transform_data = { 2.0f, 0.5f, -500.0f
                                                , 0.333f, 3.0f, -500.0f
                                                };
            cv::Mat transform( 2, 3, CV_32F, transform_data.data() );

            cv::Mat cpu_result, gpu_result, pen_result;
            std::chrono::duration<double> elapsed_time_cpu, elapsed_time_gpu_p_copy, elapsed_time_gpu_nocopy;

            {
                const auto cpu_start = std::chrono::high_resolution_clock::now();
                cv::warpAffine( cpu_gray, cpu_result, transform, cpu_gray.size() );
                const auto cpu_end = std::chrono::high_resolution_clock::now();
                elapsed_time_cpu = cpu_end - cpu_start;
            }
            {
                const auto gpu_start_copy = std::chrono::high_resolution_clock::now();
                cv::ocl::oclMat gpu_gray(cpu_gray);
                cv::ocl::oclMat gpu_affine;
                const auto gpu_start = std::chrono::high_resolution_clock::now();
                cv::ocl::warpAffine( gpu_gray, gpu_affine, transform, gpu_gray.size() );
                const auto gpu_end = std::chrono::high_resolution_clock::now();
                gpu_result = gpu_affine;
                const auto gpu_end_copy = std::chrono::high_resolution_clock::now();
                elapsed_time_gpu_p_copy = gpu_end_copy - gpu_start_copy;
                elapsed_time_gpu_nocopy = gpu_end      - gpu_start;
            }
            {
                // verifying the pencil code
                pen_result.create( cpu_gray.size(), CV_32F );
                convert_coeffs(reinterpret_cast<float*>(transform.data));

                prl_timings_start();
                pencil_affine_linear( cpu_gray.rows, cpu_gray.cols, cpu_gray.step1(), cpu_gray.ptr<float>()
                                    , pen_result.rows, pen_result.cols, pen_result.step1(), pen_result.ptr<float>(),
                        transform.at<float>(0,0), transform.at<float>(0,1), transform.at<float>(1,0), transform.at<float>(1,1),
                        transform.at<float>(1,2), transform.at<float>(0,2) );
                prl_timings_stop();
            }
            // Verifying the results
            if ( (cv::norm(cv::abs(cpu_result - gpu_result), cv::NORM_INF ) > 1 ) || (cv::norm(cv::abs(cpu_result - pen_result), cv::NORM_INF ) > 1 ) )
            {
                cv::Mat gpu_result8;
                cv::Mat cpu_result8;
                cv::Mat pen_result8;

                gpu_result.convertTo( gpu_result8, CV_8UC1, 255. );
                cpu_result.convertTo( cpu_result8, CV_8UC1, 255. );
                pen_result.convertTo( pen_result8, CV_8UC1, 255. );

                cv::imwrite( "gpu_affine.png", gpu_result8 );
                cv::imwrite( "cpu_affine.png", cpu_result8 );
                cv::imwrite( "pencil_affine.png", pen_result8 );

                throw std::runtime_error("The GPU results are not equivalent with the CPU results.");
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
    prl_init((prl_init_flags)(PRL_TARGET_DEVICE_DYNAMIC | PRL_PROFILING_ENABLED));

    std::cout << "This executable is iterating over all the files which are present in the directory `./pool'. " << std::endl;

    auto pool = carp::get_pool("pool");
#ifdef RUN_ONLY_ONE_EXPERIMENT
    time_affine( pool,  1 );
#else
    time_affine( pool, 20 );
#endif

    prl_shutdown();
    return EXIT_SUCCESS;
}
