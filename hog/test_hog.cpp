#include <chrono>
#include <random>

#include <opencv2/core/core.hpp>
#include <opencv2/ocl/ocl.hpp>

#include "opencl.hpp"
#include "utility.hpp"
#include "hog.pencil.h"
#include "HogDescriptor.h"
#include "hog.clh"

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
                carp::opencl::device device;
                
                size_t max_work_group_size;
                cl_int err;
                err = clGetDeviceInfo( device.get_device_id()
                                     , CL_DEVICE_MAX_WORK_GROUP_SIZE
                                     , sizeof(size_t)
                                     , &max_work_group_size
                                     , nullptr
                                     );
                if (err != CL_SUCCESS) throw std::runtime_error("Cannot query max work group size.");

                const auto gpu_compile_start = std::chrono::high_resolution_clock::now();
                device.source_compile( hog_opencl_cl, hog_opencl_cl_len, {"calc_histogram", "fill_zeros"} );

                const auto gpu_copy_start = std::chrono::high_resolution_clock::now();
                cl_mem gpu_gray = clCreateBuffer( device.get_context()
                                                , CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR
                                                , cpu_gray.elemSize()*cpu_gray.rows*cpu_gray.step1()
                                                , cpu_gray.data
                                                , &err
                                                );
                if (err != CL_SUCCESS) throw std::runtime_error("Cannot copy source image to GPU.");
                cl_mem gpu_locations_x = clCreateBuffer( device.get_context()
                                                       , CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR
                                                       , sizeof(int)*num_positions
                                                       , locations_x.data()
                                                       , &err
                                                       );
                if (err != CL_SUCCESS) throw std::runtime_error("Cannot copy locations_x array to GPU.");
                cl_mem gpu_locations_y = clCreateBuffer( device.get_context()
                                                       , CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR
                                                       , sizeof(int)*num_positions
                                                       , locations_y.data()
                                                       , &err
                                                       );
                if (err != CL_SUCCESS) throw std::runtime_error("Cannot copy locations_y array to GPU.");
                cl_mem gpu_hist = clCreateBuffer(device.get_context(), CL_MEM_WRITE_ONLY, sizeof(cl_float)*HISTOGRAM_BINS*num_positions, nullptr, &err);
                if (err != CL_SUCCESS) throw std::runtime_error("Cannot allocate histogram array.");
                device["fill_zeros"](gpu_hist, HISTOGRAM_BINS*num_positions).groupsize({HISTOGRAM_BINS*num_positions}, {HISTOGRAM_BINS*num_positions});
                
                const auto gpu_start = std::chrono::high_resolution_clock::now();
                device["calc_histogram"]( cpu_gray.rows
                                        , cpu_gray.cols
                                        , (int)cpu_gray.step1()
                                        , gpu_gray
                                        , num_positions
                                        , gpu_locations_x
                                        , gpu_locations_y
                                        , size
                                        , (int)ceil(size)
                                        , gpu_hist
                                        ).groupsize({4,4,4}, {num_positions, (int)ceil(size), (int)ceil(size)});
                const auto gpu_end = std::chrono::high_resolution_clock::now();

                err = clEnqueueReadBuffer(device.get_queue(), gpu_hist, CL_TRUE, 0, sizeof(cl_float)*HISTOGRAM_BINS*num_positions, gpu_result.data(), 0, nullptr, nullptr);
                if (err != CL_SUCCESS) throw std::runtime_error("Cannot copy histogram array to host.");

                err = clReleaseMemObject(gpu_hist);
                if (err != CL_SUCCESS) throw std::runtime_error("Cannot free histogram array.");
                err = clReleaseMemObject(gpu_locations_y);
                if (err != CL_SUCCESS) throw std::runtime_error("Cannot free gpu_locations_y array.");
                err = clReleaseMemObject(gpu_locations_x);
                if (err != CL_SUCCESS) throw std::runtime_error("Cannot free gpu_locations_x array.");
                err = clReleaseMemObject(gpu_gray);
                if (err != CL_SUCCESS) throw std::runtime_error("Cannot free gpu_image array.");
                
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
            if ( cv::norm( cpu_result, gpu_result, cv::NORM_INF) > cv::norm( gpu_result, cv::NORM_INF)*1e-5
              || cv::norm( cpu_result, pen_result, cv::NORM_INF) > cv::norm( cpu_result, cv::NORM_INF)*1e-5
               )
            {
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



int main(int argc, char* argv[])
{
    try {
        std::cout << "This executable is iterating over all the files which are present in the directory `./pool'. " << std::endl;

        auto pool = carp::get_pool("pool");
#ifdef RUN_ONLY_ONE_EXPERIMENT
        time_hog( pool, {64}, 50 );
#else
        time_hog( pool, {16, 32, 64, 128, 192}, 50 );
//        time_hog( pool, {16}, 1 );
#endif

        return EXIT_SUCCESS;
    }catch(const std::exception& e) {
        std::cout << e.what() << std::endl;
        return EXIT_FAILURE;
    }
} // main
