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
time_cvtColor( carp::opencl::device & device, T0 & pool, size_t iterations)
{
    cv::Mat host_gray;
    cv::Mat cpu_gray;
    cv::Mat check;
    cv::Mat pencil_gray;
    
    long int elapsed_time_gpu = 0;
    long int elapsed_time_cpu = 0;

    for(int i = 0; i < iterations; ++i) {
        PRINT(i);        
        for ( auto & record : pool ) {
            PRINT(record.path());
            cv::Mat cpuimg = record.cpuimg();                
            // CPU Bench
            {
                auto start = std::chrono::high_resolution_clock::now();
                cv::cvtColor( cpuimg, cpu_gray, CV_RGB2GRAY );
                auto end = std::chrono::high_resolution_clock::now();
                elapsed_time_cpu += carp::microseconds(end - start);
            }
            // GPU Bench
            {
                cv::ocl::oclMat gpu_gray;
                cv::ocl::oclMat gpuimg(cpuimg);
                gpu_gray.create( gpuimg.rows, gpuimg.cols, CV_8U );
                                
                // int code = CV_BGR2GRAY;                
                // int bidx = (code == CV_BGR2GRAY || code == CV_BGRA2GRAY) ? 0 : 2;
                int bidx = 2;
                auto start = std::chrono::high_resolution_clock::now();
                device["RGB2Gray"] (
                    gpuimg.cols,
                    gpuimg.rows ,
                    static_cast<int>(gpuimg.step),
                    static_cast<int>(gpu_gray.step),
                    gpuimg.channels()+1,
                    bidx,
                    reinterpret_cast<cl_mem>(gpuimg.data),
                    reinterpret_cast<cl_mem>(gpu_gray.data) 
                    )
                    .groupsize( carp::make_vector<size_t>(16,16), carp::make_vector<size_t>(cpuimg.cols,cpuimg.rows) );               
                auto end = std::chrono::high_resolution_clock::now();
                elapsed_time_gpu += carp::microseconds(end - start);
                check = gpu_gray;
                //check = gpuimg;
                
            }

            // pencil verification
            {
                int bidx = 2;
                pencil_gray.create( cpu_gray.rows, cpu_gray.cols, CV_8U );                
                
                pencil_RGB2Gray(
                    cpuimg.cols,
                    cpuimg.rows,
                    cpuimg.step,
                    pencil_gray.step,
                    cpuimg.channels(),
                    bidx,
                    cpuimg.data,
                    pencil_gray.data
                    );
                
            }
            
            // Verifying the results
            PRINT(cv::norm(pencil_gray - cpu_gray));
            if ( cv::norm(check - cpu_gray) > 0.01 ) {
                cv::imwrite( "gpu_img.png", check );
                cv::imwrite( "cpu_img.png", cpu_gray );
                throw std::runtime_error("The GPU results are not equivalent with the CPU results.");                
            }
        }
    }
    
    carp::Timing::print( "cvtColor", elapsed_time_cpu, elapsed_time_gpu );
    return;
}


int main(int argc, char* argv[])
{

    std::cout << "This executable is iterating over all the files which are present in the directory `./pool'. " << std::endl;    

    auto pool = carp::get_pool("pool");

    // Initializing OpenCL
    cv::ocl::Context * context = cv::ocl::Context::getContext();
    carp::opencl::device device(context);
    device.source_compile( cvt_color_cl, cvt_color_cl_len, carp::string_vector("RGB2Gray") );
    size_t num_iterations = 1;
    carp::Timing::printHeader();
    time_cvtColor( device, pool, num_iterations );

    return EXIT_SUCCESS;    
} // main


















// LuM end of file
