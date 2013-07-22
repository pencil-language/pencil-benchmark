// UjoImro, 2013

#include <chrono>
#include <opencv2/opencv.hpp>
#include <opencv2/ocl/ocl.hpp>

#include "utility.hpp"

template<class T0>
void
time_cvtColor( T0 & pool, size_t iterations)
{
    cv::Mat host_gray;

    // making a copy of the images on the gpu
    
    // PRINT("1");    
    const auto cpu_start = std::chrono::high_resolution_clock::now();
    for(int i = 1; i <= iterations; ++i)
        for ( auto & record : pool )
            cv::cvtColor( record.cpuimg, host_gray, CV_RGB2GRAY);
    const auto cpu_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - cpu_start);

    cv::Mat cpu_gray;
    cv::Mat check;    
    
    long int elapsed_time_gpu = 0;
    long int elapsed_time_cpu = 0;

    for(int i = 1; i <= iterations; ++i)
        for ( auto & record : pool ) {
            PRINT(record.path);            
            // CPU Bench
            {
                auto start = std::chrono::high_resolution_clock::now();
                cv::cvtColor( record.cpuimg, cpu_gray, CV_RGB2GRAY );
                auto end = std::chrono::high_resolution_clock::now();
                elapsed_time_cpu += carp::microseconds(end - start);
            }
            // GPU Bench
            {
                cv::ocl::oclMat gpu_gray;                    
                cv::ocl::oclMat gpuimg(record.cpuimg);                
                auto start = std::chrono::high_resolution_clock::now();
                cv::ocl::cvtColor( gpuimg, gpu_gray, CV_RGB2GRAY );
                auto end = std::chrono::high_resolution_clock::now();
                check = gpu_gray;
            }

            // Verifying the results
            if ( cv::norm(check - cpu_gray) > 0.01 ) {
                cv::imwrite( "gpu_img.png", check );
                cv::imwrite( "cpu_img.png", cpu_gray );
                throw std::runtime_error("The GPU results are not equivalent with the CPU results.");
            }
            
            
        }

    carp::Timing::print( "cvtColor", elapsed_time_cpu, elapsed_time_gpu );
    return;
}


int main(int argc, char* argv[])
{

    std::cout << "This executable is iterating over all the files which are present in the directory `./pool'. " << std::endl;    


    auto pool = carp::get_pool("pool");
    
//     const cv::Mat img_large = cv::imread("pool/test1.jpg");
//     cv::Mat img;
//     cv::resize(img_large, img, cv::Size(5792, 5792));   //Larger: cv::ocl::cvtColor doesn't work (returns memory garbage on some parts of the picture)

    size_t num_iterations = 1;
// //    if (argc > 1)
// //        num_iterations = std::stoi(argv[1]);

    carp::Timing::printHeader();
    
    time_cvtColor(pool, num_iterations);

    int i = 5;
} // main


















// LuM end of file
