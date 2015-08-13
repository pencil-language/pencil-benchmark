#include "utility.hpp"
#include "HogDescriptor.h"

#include <opencv2/imgproc/imgproc.hpp>


#ifndef EXCLUDE_PENCIL_TEST
#include <prl.h>
#include "hog.pencil.h"
#else
#define prl_init(x)
#define prl_shutdown(x)
#endif

//Useful set 1:
//#define NUMBER_OF_CELLS 4
//#define NUMBER_OF_BINS 8
//#define GAUSSIAN_WEIGHTS 1
//#define SPARTIAL_WEIGHTS 1
//#define SIGNED_HOG 1
//#define STATIC_HOG 1

//Useful set 2:
//#define NUMBER_OF_CELLS 2
//#define NUMBER_OF_BINS 4
//#define GAUSSIAN_WEIGHTS 1
//#define SPARTIAL_WEIGHTS 1
//#define SIGNED_HOG 0
//#define STATIC_HOG 0

#ifndef BLOCK_SIZE
#define BLOCK_SIZE 32
#endif

#ifndef NUMBER_OF_LOCATIONS
#define NUMBER_OF_LOCATIONS 49
#endif

#ifndef NUMBER_OF_CELLS
#define NUMBER_OF_CELLS 2
#endif

#ifndef NUMBER_OF_BINS
#define NUMBER_OF_BINS 4
#endif

#ifndef GAUSSIAN_WEIGHTS
#define GAUSSIAN_WEIGHTS 1
#endif

#ifndef SPARTIAL_WEIGHTS
#define SPARTIAL_WEIGHTS 1
#endif

#ifndef SIGNED_HOG
#define SIGNED_HOG 0
#endif

#ifndef STATIC_HOG
#define STATIC_HOG 0
#endif

void time_hog( const std::vector<carp::record_t>& pool, const std::vector<float>& sizes, int num_positions, int repeat )
{
    carp::Timing timing("HOG");
    bool first_execution_opencl = true;
#ifndef EXCLUDE_PENCIL_TEST
    bool first_execution_pencil = true;
#endif

    for (;repeat>0; --repeat) {
        for ( auto & size : sizes ) {
            for ( auto & item : pool ) {
                std::mt19937 rng(1);   //uses same seed, reseed for all iteration

                cv::Mat cpu_gray;
                cv::cvtColor( item.cpuimg(), cpu_gray, CV_RGB2GRAY );
                std::cout << "image path: " << item.path()   << std::endl;
                std::cout << "image rows: " << cpu_gray.rows << std::endl;
                std::cout << "image cols: " << cpu_gray.cols << std::endl;

                cv::Mat_<float> locations(num_positions, 2);
#if !STATIC_HOG
                cv::Mat_<float> blocksizes(num_positions, 2);
#else
                float blocksizes = size;
#endif
                //fill locations and blocksizes
                std::uniform_real_distribution<float> genx(size/2+1, cpu_gray.cols-1-size/2-1);
                std::uniform_real_distribution<float> geny(size/2+1, cpu_gray.rows-1-size/2-1);
                for( int i = 0; i < num_positions; ++i) {
                    locations(i, 0) = genx(rng);
                    locations(i, 1) = geny(rng);
#if !STATIC_HOG
                    blocksizes(i, 0) = size;
                    blocksizes(i, 1) = size;
#endif
                }

                cv::Mat_<float> cpu_result, gpu_result, pen_result;
                std::chrono::duration<double> elapsed_time_cpu, elapsed_time_gpu;

                {
                    //CPU implementation
                    static nel::HOGDescriptorCPP<NUMBER_OF_CELLS, NUMBER_OF_BINS, GAUSSIAN_WEIGHTS, SPARTIAL_WEIGHTS, SIGNED_HOG, STATIC_HOG> descriptor;
                    const auto cpu_start = std::chrono::high_resolution_clock::now();
                    cpu_result = descriptor.compute(cpu_gray, locations, blocksizes);
                    const auto cpu_end = std::chrono::high_resolution_clock::now();

                    elapsed_time_cpu = cpu_end - cpu_start;
                    //Free up resources
                }
                {
                    //OpenCL implementation
                    static nel::HOGDescriptorOCL<NUMBER_OF_CELLS, NUMBER_OF_BINS, GAUSSIAN_WEIGHTS, SPARTIAL_WEIGHTS, SIGNED_HOG, STATIC_HOG> descriptor;

                    //First execution includes buffer allocation
                    if (first_execution_opencl)
                    {
                        gpu_result = descriptor.compute(cpu_gray, locations, blocksizes);
                        first_execution_opencl = false;
                    }

                    const auto gpu_start = std::chrono::high_resolution_clock::now();
                    gpu_result = descriptor.compute(cpu_gray, locations, blocksizes);
                    const auto gpu_end = std::chrono::high_resolution_clock::now();

                    elapsed_time_gpu = gpu_end - gpu_start;
                    //Free up resources
                }
#ifndef EXCLUDE_PENCIL_TEST
                {
                    //PENCIL implementation
                    pen_result.create(num_positions, NUMBER_OF_CELLS * NUMBER_OF_CELLS * NUMBER_OF_BINS);

                    if (first_execution_pencil)
                    {
#if STATIC_HOG
                        pencil_hog_static (
#else
                        pencil_hog_dynamic(
#endif
                                NUMBER_OF_CELLS, NUMBER_OF_BINS, GAUSSIAN_WEIGHTS, SPARTIAL_WEIGHTS, SIGNED_HOG
                              , cpu_gray.rows, cpu_gray.cols, cpu_gray.step1(), cpu_gray.ptr<uint8_t>()
                              , num_positions
                              , reinterpret_cast<const float (*)[2]>(locations.data)
#if STATIC_HOG
                              , blocksizes
#else
                              , reinterpret_cast<const float (*)[2]>(blocksizes.data)
#endif
                              , reinterpret_cast<      float  *    >(pen_result.data)
                              );
                        first_execution_pencil = false;
                    }

                    prl_timings_reset();
                    prl_timings_start();

                    pencil_hog( NUMBER_OF_CELLS, NUMBER_OF_BINS, GAUSSIAN_WEIGHTS, SPARTIAL_WEIGHTS, SIGNED_HOG
                              , cpu_gray.rows, cpu_gray.cols, cpu_gray.step1(), cpu_gray.ptr<uint8_t>()
                              , num_positions
                              , reinterpret_cast<const float (*)[2]>(locations.data)
                              , reinterpret_cast<const float (*)[2]>(blocksizes.data)
                              , reinterpret_cast<      float  *    >(pen_result.data)
                              );

                    prl_timings_stop();
                    // Dump execution times for PENCIL code.
                    prl_timings_dump();
                }
#endif
                // Verifying the results
                if ( cv::norm( cpu_result, gpu_result, cv::NORM_INF) > cv::norm( gpu_result, cv::NORM_INF)*1e-5 )
                {
                    std::cerr << "ERROR: Results don't match" << std::endl;
                    std::cerr << "CPU norm:"     << cv::norm(cpu_result,             cv::NORM_INF) << std::endl;
                    std::cerr << "GPU norm:"     << cv::norm(gpu_result,             cv::NORM_INF) << std::endl;
                    std::cerr << "GPU-CPU norm:" << cv::norm(gpu_result, cpu_result, cv::NORM_INF) << std::endl;
                    std::cerr << "CPU result: "         << std::fixed << std::setprecision(6) << cpu_result            << std::endl;
                    std::cerr << "GPU result: "         << std::fixed << std::setprecision(6) << gpu_result            << std::endl;
                    std::cerr << "CPU-GPU difference: " << std::fixed << std::setprecision(6) << cpu_result-gpu_result << std::endl;
                    throw std::runtime_error("The OpenCL results are not equivalent with the C++ results.");
                }
#ifndef EXCLUDE_PENCIL_TEST
                if ( cv::norm( cpu_result, pen_result, cv::NORM_INF) > cv::norm( cpu_result, cv::NORM_INF)*1e-5 )
                {
                    std::cerr << "ERROR: Results don't match." << std::endl;
                    std::cerr << "CPU norm:"     << cv::norm(cpu_result,             cv::NORM_INF) << std::endl;
                    std::cerr << "PEN norm:"     << cv::norm(pen_result,             cv::NORM_INF) << std::endl;
                    std::cerr << "PEN-CPU norm:" << cv::norm(pen_result, cpu_result, cv::NORM_INF) << std::endl;
                    std::cerr << "CPU result: "         << std::fixed << std::setprecision(6) << cpu_result            << std::endl;
                    std::cerr << "PEN result: "         << std::fixed << std::setprecision(6) << pen_result            << std::endl;
                    std::cerr << "CPU-PEN difference: " << std::fixed << std::setprecision(6) << cpu_result-pen_result << std::endl;
                    throw std::runtime_error("The PENCIL results are not equivalent with the C++ results.");
                }
#endif
                timing.print(elapsed_time_cpu, elapsed_time_gpu);
            }
        }
    }
}

int main(int argc, char* argv[])
{
    try
    {
        prl_init((prl_init_flags)(PRL_TARGET_DEVICE_DYNAMIC | PRL_PROFILING_ENABLED));

        std::cout << "This executable is iterating over all the files passed to it as an argument. " << std::endl;
        auto pool = carp::get_pool(argc, argv);

#ifdef RUN_ONLY_ONE_EXPERIMENT
        time_hog( pool, {BLOCK_SIZE}, NUMBER_OF_LOCATIONS, 1 );
#else
        time_hog( pool, {16, 32, 64, 128}, NUMBER_OF_LOCATIONS, 27 );
#endif

        prl_shutdown();

        return EXIT_SUCCESS;
    }
    catch(const std::exception& e)
    {
        std::cout << e.what() << std::endl;

        prl_shutdown();
        return EXIT_FAILURE;
    }
} // main
