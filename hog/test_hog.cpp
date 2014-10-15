#include <chrono>
#include <random>
#include <array>

#include <opencv2/core/core.hpp>
#include <opencv2/ocl/ocl.hpp>

#include "opencl.hpp"
#include "utility.hpp"
#include "hog.pencil.h"
#include "HogDescriptor.h"
#include "hog.clh"

#define NUMBER_OF_CELLS 1
#define NUMBER_OF_BINS 8
#define GAUSSIAN_WEIGHTS 1
#define SPARTIAL_WEIGHTS 0
#define SIGNED_HOG 1

void time_hog( const std::vector<carp::record_t>& pool, const std::vector<float>& sizes, int num_positions, int repeat )
{
    carp::Timing timing("HOG");

    for (;repeat>0; --repeat) {
        for ( auto & size : sizes ) {
            for ( auto & item : pool ) {
                std::mt19937 rng(0);   //uses same seed, reseed for all iteration

                cv::Mat cpu_gray;
                cv::cvtColor( item.cpuimg(), cpu_gray, CV_RGB2GRAY );

                cv::Mat_<float> locations(num_positions, 2);
                cv::Mat_<float> blocksizes(num_positions, 2);
                //fill locations and blocksizes
                std::uniform_real_distribution<float> genx(size/2+1, cpu_gray.rows-1-size/2-1);
                std::uniform_real_distribution<float> geny(size/2+1, cpu_gray.cols-1-size/2-1);
                for( int i = 0; i < num_positions; ++i) {
                    locations(i, 0) = genx(rng);
                    locations(i, 1) = geny(rng);
                    blocksizes(i, 0) = size;
                    blocksizes(i, 1) = size;
                }

                const int HISTOGRAM_BINS = NUMBER_OF_CELLS * NUMBER_OF_CELLS * NUMBER_OF_BINS;
                std::vector<float> cpu_result(num_positions * HISTOGRAM_BINS), gpu_result(num_positions * HISTOGRAM_BINS), pen_result(num_positions * HISTOGRAM_BINS);
                std::chrono::duration<double> elapsed_time_cpu, elapsed_time_gpu_p_copy, elapsed_time_gpu_nocopy, elapsed_time_pencil;

                {
                    //CPU implement
                    static nel::HOGDescriptorCPP descriptor( NUMBER_OF_CELLS
                                                           , NUMBER_OF_BINS
                                                           , GAUSSIAN_WEIGHTS
                                                           , SPARTIAL_WEIGHTS
                                                           , SIGNED_HOG
                                                           );
                    const auto cpu_start = std::chrono::high_resolution_clock::now();
                    const auto result = descriptor.compute(cpu_gray, locations, blocksizes);
                    const auto cpu_end = std::chrono::high_resolution_clock::now();

                    std::copy(result.begin(), result.end(), cpu_result.begin());
                    elapsed_time_cpu = cpu_end - cpu_start;
                    //Free up resources
                }
                {
                    //GPU implement
                    static nel::HOGDescriptorOCL descriptor( NUMBER_OF_CELLS
                                                           , NUMBER_OF_BINS
                                                           , GAUSSIAN_WEIGHTS
                                                           , SPARTIAL_WEIGHTS
                                                           , SIGNED_HOG
                                                           );
                    const auto gpu_start = std::chrono::high_resolution_clock::now();
                    const auto result = descriptor.compute(cpu_gray, locations, blocksizes, elapsed_time_gpu_nocopy);
                    const auto gpu_end = std::chrono::high_resolution_clock::now();

                    std::copy(result.begin(), result.end(), gpu_result.begin());
                    elapsed_time_gpu_p_copy = gpu_end - gpu_start;
                    //Free up resources
                }
                {
                    pen_result.resize(num_positions * HISTOGRAM_BINS, 0.0f);
                    const auto pencil_start = std::chrono::high_resolution_clock::now();
                    pencil_hog( NUMBER_OF_CELLS, NUMBER_OF_BINS, GAUSSIAN_WEIGHTS, SPARTIAL_WEIGHTS, SIGNED_HOG
                              , cpu_gray.rows, cpu_gray.cols, cpu_gray.step1(), cpu_gray.ptr<uint8_t>()
                              , num_positions
                              , reinterpret_cast<const float (*)[2]>(locations.data)
                              , reinterpret_cast<const float (*)[2]>(blocksizes.data)
                              , pen_result.data()
                              );

                    const auto pencil_end = std::chrono::high_resolution_clock::now();
                    elapsed_time_pencil = pencil_end - pencil_start;
                    //Free up resources
                }
                // Verifying the results
                if ( cv::norm( cpu_result, gpu_result, cv::NORM_INF) > cv::norm( gpu_result, cv::NORM_INF)*1e-5
                  || cv::norm( cpu_result, pen_result, cv::NORM_INF) > cv::norm( cpu_result, cv::NORM_INF)*1e-5
                   )
                {
                    std::vector<float> diff;
                    std::transform(cpu_result.begin(), cpu_result.end(), gpu_result.begin(), std::back_inserter(diff), std::minus<float>());

                    std::cerr << "ERROR: Results don't match." << std::endl;
                    std::cerr << "CPU norm:" << cv::norm(cpu_result, cv::NORM_INF) << std::endl;
                    std::cerr << "GPU norm:" << cv::norm(gpu_result, cv::NORM_INF) << std::endl;
                    std::cerr << "PEN norm:" << cv::norm(pen_result, cv::NORM_INF) << std::endl;
                    std::cerr << "GPU-CPU norm:" << cv::norm(gpu_result, cpu_result, cv::NORM_INF) << std::endl;
                    std::cerr << "PEN-CPU norm:" << cv::norm(pen_result, cpu_result, cv::NORM_INF) << std::endl;

                    throw std::runtime_error("The OpenCL or PENCIL results are not equivalent with the C++ results.");
                }
                timing.print( elapsed_time_cpu, elapsed_time_gpu_p_copy, elapsed_time_gpu_nocopy, elapsed_time_pencil );
            }
        }
    }
}

int main(int argc, char* argv[])
{
    try {
        std::cout << "This executable is iterating over all the files which are present in the directory `./pool'. " << std::endl;

        auto pool = carp::get_pool("pool");
#ifdef RUN_ONLY_ONE_EXPERIMENT
        time_hog( pool, {64}, 50, 1 );
#else
        time_hog( pool, {16, 32, 64, 128, 192}, 50, 10 );
#endif

        return EXIT_SUCCESS;
    }catch(const std::exception& e) {
        std::cout << e.what() << std::endl;
        return EXIT_FAILURE;
    }
} // main
