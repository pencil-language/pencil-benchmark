#include <chrono>
#include <random>

#include <opencv2/opencv.hpp>
//#include <opencv2/ocl/ocl.hpp>
//#include <opencv2/imgproc/imgproc.hpp>

#include "opencl.hpp"
#include "utility.hpp"
#include "hog.pencil.h"
#include "HogDescriptor.h"

void time_hog( const std::vector<carp::record_t>& pool, const std::vector<float>& sizes, int num_positions )
{
    carp::Timing timing("HOG");

    for ( auto & size : sizes ) {
        for ( auto & item : pool ) {
            std::mt19937 rng(0);   //uses same seed, reseed for all iteration
            
            cv::Mat cpu_gray;
            cv::cvtColor( item.cpuimg(), cpu_gray, CV_RGB2GRAY );

            std::vector<float> locations_x, locations_y;
            std::uniform_real_distribution<float> genx(size/2+1, cpu_gray.rows-1-size/2-1);
            std::uniform_real_distribution<float> geny(size/2+1, cpu_gray.cols-1-size/2-1);
            std::generate_n(std::back_inserter(locations_x), num_positions, [&](){ return genx(rng); });
            std::generate_n(std::back_inserter(locations_y), num_positions, [&](){ return geny(rng); });
            
            std::vector<float> cpu_result(num_positions * HISTOGRAM_BINS), gpu_result(num_positions * HISTOGRAM_BINS), pen_result(num_positions * HISTOGRAM_BINS);
            std::chrono::duration<double> elapsed_time_cpu, elapsed_time_gpu_p_copy, elapsed_time_gpu_nocopy, elapsed_time_pencil;

            {
                //CPU implement
                const auto cpu_start = std::chrono::high_resolution_clock::now();

                auto result = nel::HOGDescriptor< NUMBER_OF_CELLS
                                                , NUMBER_OF_BINS
                                                , GAUSSIAN_WEIGHTS
                                                , SPARTIAL_WEIGHTS
                                                , SIGNED_HOG
                                                >::compute(cpu_gray, locations_x, locations_y, size);

                const auto cpu_end = std::chrono::high_resolution_clock::now();

                std::copy(result.begin(), result.end(), cpu_result.begin());
                elapsed_time_cpu = cpu_end - cpu_start;
                //Free up resources
            }
            {
                const auto gpu_compile_start = std::chrono::high_resolution_clock::now();
                const auto gpu_copy_start = std::chrono::high_resolution_clock::now();
                const auto gpu_start = std::chrono::high_resolution_clock::now();

                //!!!TODO!!!

                const auto gpu_end = std::chrono::high_resolution_clock::now();
                const auto gpu_copy_end = std::chrono::high_resolution_clock::now();
                const auto gpu_compile_end = std::chrono::high_resolution_clock::now();

                elapsed_time_gpu_p_copy = gpu_copy_end - gpu_copy_start;
                elapsed_time_gpu_nocopy = gpu_end      - gpu_start;
                auto elapsed_time_gpu_compile = gpu_compile_end - gpu_compile_start;
                //Free up resources
            }
            {
                pen_result.resize(num_positions * HISTOGRAM_BINS, 0.0f);
                const auto pencil_start = std::chrono::high_resolution_clock::now();
                pencil_hog( cpu_gray.rows, cpu_gray.cols, cpu_gray.step1(), cpu_gray.ptr<uint8_t>()
                          , num_positions, locations_x.data(), locations_y.data()
                          , size
                          , pen_result.data()
                          );

                const auto pencil_end = std::chrono::high_resolution_clock::now();
                elapsed_time_pencil = pencil_end - pencil_start;
                //Free up resources
            }
            // Verifying the results
            if ( //cv::norm( cpu_result, gpu_result, cv::NORM_INF) > cv::norm( gpu_result, cv::NORM_INF)*1e-6 ) ||
                 (cv::norm( cpu_result, pen_result, cv::NORM_INF) > cv::norm( cpu_result, cv::NORM_INF)*1e-6 )
               )
            {
                std::cout << "CPU norm:" << cv::norm(cpu_result, cv::NORM_INF) << std::endl;
//                std::cout << "GPU norm:" << cv::norm(gpu_result, cv::NORM_INF) << std::endl;
                std::cout << "PEN norm:" << cv::norm(pen_result, cv::NORM_INF) << std::endl;
//                std::cout << "GPU-CPU norm:" << cv::norm(gpu_result, cpu_result, cv::NORM_INF) << std::endl;
                std::cout << "PEN-CPU norm:" << cv::norm(pen_result, cpu_result, cv::NORM_INF) << std::endl;

#if 0
                cv::Mat gpu_result8;
                cv::Mat cpu_result8;
                cv::Mat pen_result8;
                cv::Mat diff8;

                gpu_result.convertTo( gpu_result8, CV_8UC1, 255. );
                cpu_result.convertTo( cpu_result8, CV_8UC1, 255. );
                pen_result.convertTo( pen_result8, CV_8UC1, 255. );
                cv::Mat absdiff = cv::abs(pen_result - cpu_result);
                absdiff.convertTo( diff8, CV_8UC1, 255. );

                cv::imwrite( "gpu_gaussian.png", gpu_result8 );
                cv::imwrite( "cpu_gaussian.png", cpu_result8 );
                cv::imwrite( "pencil_gaussian.png", pen_result8 );
                cv::imwrite( "diff_gaussian.png", diff8 );
#endif
                throw std::runtime_error("The GPU results are not equivalent with the CPU results.");
            }
            timing.print( elapsed_time_cpu, elapsed_time_gpu_p_copy, elapsed_time_gpu_nocopy, elapsed_time_pencil );
        }
    }
}



int main(int argc, char* argv[])
{
    try {
        std::cout << "This executable is iterating over all the files which are present in the directory `./pool'. " << std::endl;

        auto pool = carp::get_pool("pool");
#ifdef RUN_ONLY_ONE_EXPERIMENT
        time_hog( pool, {64}, 50 );
#else
        time_hog( pool, {16, 32, 64, 128, 192}, 50 );
#endif

        return EXIT_SUCCESS;
    }catch(...) {
        return EXIT_FAILURE;
    }
} // main
