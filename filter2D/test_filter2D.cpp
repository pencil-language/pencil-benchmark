// Experimental Research Code for the CARP Project
// UjoImro, 2013

#include <chrono>
#include <opencv2/opencv.hpp>
#include <opencv2/ocl/ocl.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "opencl.hpp"
#include "utility.hpp"
#include "imgproc_convolve.clh"

namespace {
    inline int divUp(int total, int grain)
    {
        return (total + grain - 1) / grain;
    }
} // unnamed namespace


template<class T0>
void
time_filter2D( carp::opencl::device & device, T0 & pool, int iteration )
{

    double sum_quotient = 0;
    int64_t nums = 0;
        
    for ( int q=0; q<iteration; q++ ) {
        
        for ( auto & item : pool ) {
            PRINT(item.path());

            long int elapsed_time_gpu = 0;
            long int elapsed_time_cpu = 0;
        
            cv::Mat cpu_gray;
            cv::Mat check;    

            cv::cvtColor( item.cpuimg(), cpu_gray, CV_RGB2GRAY );
            cpu_gray.convertTo( cpu_gray, CV_32F, 1.0/255. );
            cv::Mat host_convolve;
            
            float kernel_data[] = {-1, -1, -1
                                   , 0,  0,  0
                                   , 1,  1,  1
            };
            cv::Mat kernel_cpu(3, 3, CV_32F, kernel_data);

            const auto cpu_start = std::chrono::high_resolution_clock::now();
            cv::filter2D( cpu_gray, host_convolve, -1, kernel_cpu, cv::Point(-1,-1), 0.0, cv::BORDER_REPLICATE );
            const auto cpu_end = std::chrono::high_resolution_clock::now();
            elapsed_time_cpu += carp::microseconds(cpu_end - cpu_start);

            {
                cv::ocl::oclMat gpu_gray(cpu_gray);
                cv::ocl::oclMat gpu_convolve;
                cv::ocl::oclMat kernel_gpu(kernel_cpu);                
                
                CV_Assert(gpu_gray.depth() == CV_32F);
                CV_Assert(kernel_gpu.depth() == CV_32F);
                gpu_convolve.create(gpu_gray.size(), gpu_gray.type());
                CV_Assert(gpu_gray.type() == kernel_gpu.type() && gpu_gray.size() == gpu_convolve.size());

                std::string kernelName = "convolve_D5";

                CV_Assert(gpu_gray.depth() == CV_32FC1);
                CV_Assert(kernel_gpu.depth() == CV_32F);
                CV_Assert(kernel_gpu.cols <= 17 && kernel_gpu.rows <= 17);

                gpu_convolve.create(gpu_gray.size(), gpu_gray.type());

                CV_Assert(gpu_gray.cols == gpu_convolve.cols && gpu_gray.rows == gpu_convolve.rows);
                CV_Assert(gpu_gray.type() == gpu_convolve.type());

                cv::ocl::Context  *clCxt = gpu_gray.clCxt;
                int channels = gpu_convolve.oclchannels();
                int depth = gpu_convolve.depth();

                size_t vector_length = 1;
                int offset_cols = ((gpu_convolve.offset % gpu_convolve.step) / gpu_convolve.elemSize1()) & (vector_length - 1);
                int cols = divUp(gpu_convolve.cols * channels + offset_cols, vector_length);
                int rows = gpu_convolve.rows;

                std::vector<size_t> localThreads  = { 16, 16, 1 };
                std::vector<size_t> globalThreads = { divUp(cols, localThreads[0]) *localThreads[0],
                                            divUp(rows, localThreads[1]) *localThreads[1],
                                                      1 };
                
                auto gpu_start = std::chrono::high_resolution_clock::now();
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
                auto gpu_end = std::chrono::high_resolution_clock::now();
                elapsed_time_gpu += carp::microseconds(gpu_end - gpu_start);
                check = gpu_convolve;
            }
            
            // Verifying the results
            if ( cv::norm(host_convolve - check) > 0.01 ) {
                PRINT(cv::norm(check - host_convolve));
                // no use to write out the results, as they are in float
                throw std::runtime_error("The GPU results are not equivalent with the CPU results.");                
            }

            if (elapsed_time_gpu > 1) {
                sum_quotient += static_cast<double>(elapsed_time_cpu) / elapsed_time_gpu;
                nums++;
            }
                        
            carp::Timing::print( "convolve image", elapsed_time_cpu, elapsed_time_gpu );

        } // for pool
            
    } // for q 

    std::cout << "Cumulated Speed Improvement: " << (sum_quotient/nums) << "x" << std::endl;    


    return;
} // text_boxFilter


int main(int argc, char* argv[])
{

    std::cout << "This executable is iterating over all the files which are present in the directory `./pool'. " << std::endl;    

    auto pool = carp::get_pool("pool");

    // Initializing OpenCL
    cv::ocl::Context * context = cv::ocl::Context::getContext();
    carp::Timing::printHeader();
    carp::opencl::device device(context);
    device.source_compile( imgproc_convolve_cl, imgproc_convolve_cl_len,
                           carp::string_vector("convolve_D5" ) );
    time_filter2D( device, pool, 1 );
    return EXIT_SUCCESS;    
} // main


















// LuM end of file
