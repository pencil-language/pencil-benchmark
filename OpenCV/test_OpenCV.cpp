// Robert DAVID, 2013
// Experimental Research Code for the CARP Project

#include <chrono>
#include <string>
#include <iomanip>
#include <opencv/cv.h>
#include <opencv2/opencv.hpp>
#include <opencv2/ocl/ocl.hpp>
#include <boost/filesystem.hpp>

#include "utility.hpp"


void
time_cvtColor(const cv::Mat& host_img, size_t iterations)
{
    cv::Mat host_gray;

    const auto cpu_start = std::chrono::high_resolution_clock::now();
    for(int i = 1; i <= iterations; ++i)
        cv::cvtColor(host_img, host_gray, CV_RGB2GRAY);
    const auto cpu_time = carp::microseconds( std::chrono::high_resolution_clock::now() - cpu_start);

    const cv::ocl::oclMat gpu_img(host_img);
    cv::ocl::oclMat gpu_gray;

    const auto gpu_start = std::chrono::high_resolution_clock::now();
    for(int i = 1; i <= iterations; ++i)
        cv::ocl::cvtColor(gpu_img, gpu_gray, CV_RGB2GRAY);
    const auto gpu_time = carp::microseconds( std::chrono::high_resolution_clock::now() - gpu_start );

    cv::Mat difference = host_gray - gpu_gray;
    if( cv::norm(difference, cv::NORM_INF) > .01 ) throw std::logic_error("cvtColor: CPU and GPU result differs too much");

    carp::Timing::print("cvtColor", cpu_time, gpu_time );
    return;
}

void
time_boxFilter(const cv::Mat& host_img, size_t iterations)
{
    cv::Mat host_filtered(host_img.size(), host_img.type());

    const auto cpu_start = std::chrono::high_resolution_clock::now();
    for(int i = 1; i <= iterations; ++i)
        cv::boxFilter(host_img, host_filtered, -1, cv::Size(5,5));
    const auto cpu_time = carp::microseconds( std::chrono::high_resolution_clock::now() - cpu_start );

    const cv::ocl::oclMat gpu_img(host_img);
    cv::ocl::oclMat gpu_filtered(gpu_img.size(), gpu_img.type());

    const auto gpu_start = std::chrono::high_resolution_clock::now();
    for(int i = 1; i <= iterations; ++i)
        cv::ocl::boxFilter(gpu_img, gpu_filtered, -1, cv::Size(5,5));
    const auto gpu_time = carp::microseconds( std::chrono::high_resolution_clock::now() - gpu_start );

    cv::Mat difference = host_filtered - gpu_filtered;
    if( cv::norm(difference, cv::NORM_INF) > 1.0 ) throw std::logic_error("boxFilter: CPU and GPU result differs too much");

    carp::Timing::print( "boxFilter", cpu_time, gpu_time );
}

void
time_integral(const cv::Mat& img, size_t iterations)
{
    cv::Mat host_img;
    {
        cv::Mat img_gray;
        cv::cvtColor(img, img_gray, CV_RGB2GRAY);
        cv::resize(img_gray, host_img, cv::Size(500, 500));
    }
    cv::Mat host_integral;

    const auto cpu_start = std::chrono::high_resolution_clock::now();
    for(int i = 1; i <= iterations; ++i)
        cv::integral(host_img, host_integral, CV_32S);
    const auto cpu_time = carp::microseconds( std::chrono::high_resolution_clock::now() - cpu_start );

    const cv::ocl::oclMat gpu_img(host_img);
    cv::ocl::oclMat gpu_integral;

    const auto gpu_start = std::chrono::high_resolution_clock::now();
    for(int i = 1; i <= iterations; ++i)
        cv::ocl::integral(gpu_img, gpu_integral);
    const auto gpu_time = carp::microseconds( std::chrono::high_resolution_clock::now() - gpu_start );

    cv::Mat difference = host_integral - gpu_integral;
    if( cv::norm(difference, cv::NORM_INF) > 1e-5 ) throw std::logic_error("integral: CPU and GPU result differs too much");

    carp::Timing::print("integral", cpu_time, gpu_time);
    return;        
}

void
time_dilate( const cv::Mat& host_img, size_t iterations )
{
    cv::Mat host_dilate;

    cv::Mat element = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(7,7));
    const auto cpu_start = std::chrono::high_resolution_clock::now();
    for(int i = 1; i <= iterations; ++i)
        cv::dilate(host_img, host_dilate, element);
    const auto cpu_time = carp::microseconds( std::chrono::high_resolution_clock::now() - cpu_start );

    const cv::ocl::oclMat gpu_img(host_img);
    cv::ocl::oclMat gpu_dilate;

    const auto gpu_start = std::chrono::high_resolution_clock::now();
    for(int i = 1; i <= iterations; ++i)
        cv::ocl::dilate(gpu_img, gpu_dilate,element);
    const auto gpu_time = carp::microseconds( std::chrono::high_resolution_clock::now() - gpu_start);

    cv::Mat difference = host_dilate - gpu_dilate;
    if( cv::norm(difference, cv::NORM_INF) > 1e-5 ) throw std::logic_error("dilate: CPU and GPU result differs too much");

    carp::Timing::print("dilate", cpu_time, gpu_time);
    return;
}

void
time_convolve(const cv::Mat& img, size_t iterations)
{
    //cv::ocl only support CV_32FC1 images
    cv::Mat host_img;
    {
        cv::Mat img_gray;
        cv::cvtColor(img, img_gray, CV_BGR2GRAY);
        img_gray.convertTo(host_img, CV_32F, 1.0/255.0);
    }
    cv::Mat host_convolve;

    float kernel_data[] = {-1, -1, -1
                          , 0,  0,  0
                          , 1,  1,  1
                          };
    cv::Mat kernel_cpu(3, 3, CV_32F, kernel_data);
    const auto cpu_start = std::chrono::high_resolution_clock::now();
    for(int i = 1; i <= iterations; ++i)
        cv::filter2D(host_img, host_convolve, -1, kernel_cpu, cv::Point(-1,-1), 0.0, cv::BORDER_REPLICATE); //cv::ocl has only a convolve functions, with the extra parameters working like this
    const auto cpu_time = carp::microseconds( std::chrono::high_resolution_clock::now() - cpu_start);

    const cv::ocl::oclMat gpu_img(host_img);
    const cv::ocl::oclMat kernel_gpu(kernel_cpu);
    cv::ocl::oclMat gpu_convolve;

    const auto gpu_start = std::chrono::high_resolution_clock::now();
    for(int i = 1; i <= iterations; ++i)
        cv::ocl::convolve(gpu_img, kernel_gpu, gpu_convolve);
    const auto gpu_time = carp::microseconds( std::chrono::high_resolution_clock::now() - gpu_start);

    cv::Mat difference = host_convolve - gpu_convolve;
    if( cv::norm(difference, cv::NORM_INF) > 1e-5 ) throw std::logic_error("filter2D: CPU and GPU result differs too much");

    carp::Timing::print("filter2D", cpu_time, gpu_time);
    return;
}

void
time_gaussian(const cv::Mat& img, size_t iterations)
{
    cv::Mat host_img;
    cv::resize(img, host_img, cv::Size(500, 500));
    cv::Mat host_gaussian;

    const auto cpu_start = std::chrono::high_resolution_clock::now();
    for(int i = 1; i <= iterations; ++i)
        cv::GaussianBlur(host_img, host_gaussian, cv::Size(5,5), 0, 0, cv::BORDER_REPLICATE);
    const auto cpu_time = carp::microseconds( std::chrono::high_resolution_clock::now() - cpu_start);

    const cv::ocl::oclMat gpu_img(host_img);
    cv::ocl::oclMat gpu_gaussian;

    const auto gpu_start = std::chrono::high_resolution_clock::now();
    for(int i = 1; i <= iterations; ++i)
        cv::ocl::GaussianBlur(gpu_img, gpu_gaussian, cv::Size(5,5), 0, 0, cv::BORDER_REPLICATE);
    const auto gpu_time = carp::microseconds( std::chrono::high_resolution_clock::now() - gpu_start);

    cv::Mat difference = host_gaussian - gpu_gaussian;
    if( cv::norm(difference, cv::NORM_INF) > 1 ) throw std::logic_error("GaussianBlur: CPU and GPU result differs too much");

    carp::Timing::print( "GaussianBlur", cpu_time, gpu_time );
    return;
}

void
time_resize(const cv::Mat& host_img, size_t iterations)
{
    cv::Mat host_resize;

    const auto cpu_start = std::chrono::high_resolution_clock::now();
    for(int i = 1; i <= iterations; ++i)
        cv::resize(host_img, host_resize, cv::Size(1000,1000));
    const auto cpu_time = carp::microseconds( std::chrono::high_resolution_clock::now() - cpu_start);

    const cv::ocl::oclMat gpu_img(host_img);
    cv::ocl::oclMat gpu_resize;

    const auto gpu_start = std::chrono::high_resolution_clock::now();
    for(int i = 1; i <= iterations; ++i)
        cv::ocl::resize(gpu_img, gpu_resize, cv::Size(1000,1000));
    const auto gpu_time = carp::microseconds( std::chrono::high_resolution_clock::now() - gpu_start);

    cv::Mat difference = host_resize - gpu_resize;
    if( cv::norm(difference, cv::NORM_INF) > 1 ) throw std::logic_error("resize: CPU and GPU result differs too much");

    carp::Timing::print("resize", cpu_time, gpu_time);
    return;
}

void
time_warpAffine(const cv::Mat& host_img, size_t iterations)
{
    cv::Mat host_warp;

    float transform_data[] = { 2.0f  , 0.5f, -500.0f
                             , 0.333f, 3.0f, -500.0f
                             };
    cv::Mat transform(2, 3, CV_32F, transform_data);

    const auto cpu_start = std::chrono::high_resolution_clock::now();
    for(int i = 1; i <= iterations; ++i)
        cv::warpAffine(host_img, host_warp, transform, cv::Size(1000,1000));
    const auto cpu_time = carp::microseconds( std::chrono::high_resolution_clock::now() - cpu_start);

    const cv::ocl::oclMat gpu_img(host_img);
    cv::ocl::oclMat gpu_warp;

    const auto gpu_start = std::chrono::high_resolution_clock::now();
    for(int i = 1; i <= iterations; ++i)
        cv::ocl::warpAffine(gpu_img, gpu_warp, transform, cv::Size(1000,1000));
    const auto gpu_time = carp::microseconds( std::chrono::high_resolution_clock::now() - gpu_start);

    cv::Mat difference = host_warp - gpu_warp;
    if( cv::norm(difference, cv::NORM_INF) > 1 ) throw std::logic_error("warpAffine: CPU and GPU result differs too much");

    carp::Timing::print("warpAffine", cpu_time, gpu_time );
    return;
}

int main(int argc, char* argv[])
{


    if ( !boost::filesystem::exists( "test1.jpg" ) )
    {
        throw std::runtime_error("Can't find `test1.jpg' file!");
    }
    
    const cv::Mat img_large = cv::imread("test1.jpg");
    cv::Mat img;
    cv::resize(img_large, img, cv::Size(5792, 5792));   //Larger: cv::ocl::cvtColor doesn't work (returns memory garbage on some parts of the picture)

    size_t num_iterations = 1;
    if (argc > 1)
        num_iterations = std::stoi(argv[1]);

    carp::Timing::printShortHeader();

    try{
        time_boxFilter (img, num_iterations);
        time_integral  (img, num_iterations);
        time_cvtColor  (img, num_iterations);
        time_dilate    (img, num_iterations);
        time_convolve  (img, num_iterations);
        time_gaussian  (img, num_iterations);
        time_resize    (img, num_iterations);
        time_warpAffine(img, num_iterations);
    }
    catch(const std::logic_error& err)
    {
        std::cout << "error: " << err.what();
    }

    int i = 5;
} // main
