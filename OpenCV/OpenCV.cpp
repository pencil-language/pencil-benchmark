#include <opencv/cv.h>
#include <opencv2/opencv.hpp>
#include <opencv2/ocl/ocl.hpp>

#include <chrono>
#include <iomanip>

struct Timing
{
    Timing(const std::chrono::milliseconds& cpu, const std::chrono::milliseconds& gpu)
        : cpu(cpu)
        , gpu(gpu)
    {}
    std::chrono::milliseconds cpu;
    std::chrono::milliseconds gpu;
};

Timing time_cvtColor(const cv::Mat& host_img, size_t iterations)
{
    cv::Mat host_gray;

    const auto cpu_start = std::chrono::high_resolution_clock::now();
    for(int i = 1; i <= iterations; ++i)
        cv::cvtColor(host_img, host_gray, CV_RGB2GRAY);
    const auto cpu_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - cpu_start);

    const cv::ocl::oclMat gpu_img(host_img);
    cv::ocl::oclMat gpu_gray;

    const auto gpu_start = std::chrono::high_resolution_clock::now();
    for(int i = 1; i <= iterations; ++i)
        cv::ocl::cvtColor(gpu_img, gpu_gray, CV_RGB2GRAY);
    const auto gpu_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - gpu_start);

#ifdef _DEBUG
    cv::Mat difference = host_gray - gpu_gray;
    assert( cv::norm(difference, cv::NORM_INF) <.01 );
#endif
    
    return Timing(cpu_time, gpu_time);
}

Timing time_boxFilter(const cv::Mat& host_img, size_t iterations)
{
    cv::Mat host_filtered(host_img.size(), host_img.type());

    const auto cpu_start = std::chrono::high_resolution_clock::now();
    for(int i = 1; i <= iterations; ++i)
        cv::boxFilter(host_img, host_filtered, -1, cv::Size(5,5));
    const auto cpu_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - cpu_start);

    const cv::ocl::oclMat gpu_img(host_img);
    cv::ocl::oclMat gpu_filtered(gpu_img.size(), gpu_img.type());

    const auto gpu_start = std::chrono::high_resolution_clock::now();
    for(int i = 1; i <= iterations; ++i)
        cv::ocl::boxFilter(gpu_img, gpu_filtered, -1, cv::Size(5,5));
    const auto gpu_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - gpu_start);

#ifdef _DEBUG
    cv::Mat difference = host_filtered - gpu_filtered;
    assert( cv::norm(difference, cv::NORM_INF) <= 1.0 );
#endif

    return Timing(cpu_time, gpu_time);
}

Timing time_integral(const cv::Mat& img, size_t iterations)
{
    cv::Mat host_img;
    cv::cvtColor(img, host_img, CV_RGB2GRAY);
    cv::Mat host_integral;//(host_img.size()+cv::Size(1,1), CV_64F);

    const auto cpu_start = std::chrono::high_resolution_clock::now();
    for(int i = 1; i <= iterations; ++i)
        cv::integral(host_img, host_integral, CV_32F);
    const auto cpu_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - cpu_start);

    const cv::ocl::oclMat gpu_img(host_img);
    cv::ocl::oclMat gpu_integral;//(gpu_img.size()+cv::Size(1,1), CV_64F);

    const auto gpu_start = std::chrono::high_resolution_clock::now();
    for(int i = 1; i <= iterations; ++i)
        cv::ocl::integral(gpu_img, gpu_integral);
    const auto gpu_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - gpu_start);

#ifdef _DEBUG
    cv::Mat difference = host_integral - gpu_integral;
    assert( cv::norm(difference, cv::NORM_INF) <= 1e-5 );
#endif

    return Timing(cpu_time, gpu_time);
}

Timing time_dilate(const cv::Mat& host_img, size_t iterations)
{
    cv::Mat host_dilate;

    cv::Mat element = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(7,7));
    const auto cpu_start = std::chrono::high_resolution_clock::now();
    for(int i = 1; i <= iterations; ++i)
        cv::dilate(host_img, host_dilate, element);
    const auto cpu_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - cpu_start);

    const cv::ocl::oclMat gpu_img(host_img);
    cv::ocl::oclMat gpu_dilate;

    const auto gpu_start = std::chrono::high_resolution_clock::now();
    for(int i = 1; i <= iterations; ++i)
        cv::ocl::dilate(gpu_img, gpu_dilate,element);
    const auto gpu_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - gpu_start);

#ifdef _DEBUG
    cv::Mat difference = host_dilate - gpu_dilate;
    assert( cv::norm(difference, cv::NORM_INF) <= 1e-5 );
#endif

    return Timing(cpu_time, gpu_time);
}

Timing time_convolve(const cv::Mat& img, size_t iterations)
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
    const auto cpu_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - cpu_start);

    const cv::ocl::oclMat gpu_img(host_img);
    const cv::ocl::oclMat kernel_gpu(kernel_cpu);
    cv::ocl::oclMat gpu_convolve;

    const auto gpu_start = std::chrono::high_resolution_clock::now();
    for(int i = 1; i <= iterations; ++i)
        cv::ocl::convolve(gpu_img, kernel_gpu, gpu_convolve);
    const auto gpu_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - gpu_start);

#ifdef _DEBUG
    cv::Mat difference = host_convolve - gpu_convolve;
    assert( cv::norm(difference, cv::NORM_INF) <= 1e-5 );
#endif

    return Timing(cpu_time, gpu_time);
}

Timing time_gaussian(const cv::Mat& host_img, size_t iterations)
{
    cv::Mat host_gaussian;

    const auto cpu_start = std::chrono::high_resolution_clock::now();
    for(int i = 1; i <= iterations; ++i)
        cv::GaussianBlur(host_img, host_gaussian, cv::Size(5,5), 0, 0, cv::BORDER_REPLICATE);
    const auto cpu_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - cpu_start);

    const cv::ocl::oclMat gpu_img(host_img);
    cv::ocl::oclMat gpu_gaussian;

    const auto gpu_start = std::chrono::high_resolution_clock::now();
    for(int i = 1; i <= iterations; ++i)
        cv::ocl::GaussianBlur(gpu_img, gpu_gaussian, cv::Size(5,5), 0, 0, cv::BORDER_REPLICATE);
    const auto gpu_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - gpu_start);

#ifdef _DEBUG
    cv::Mat difference = host_gaussian - gpu_gaussian;
    assert( cv::norm(difference, cv::NORM_INF) <= 1 );
#endif

    return Timing(cpu_time, gpu_time);
}

Timing time_resize(const cv::Mat& host_img, size_t iterations)
{
    cv::Mat host_resize;

    const auto cpu_start = std::chrono::high_resolution_clock::now();
    for(int i = 1; i <= iterations; ++i)
        cv::resize(host_img, host_resize, cv::Size(1000,1000));
    const auto cpu_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - cpu_start);

    const cv::ocl::oclMat gpu_img(host_img);
    cv::ocl::oclMat gpu_resize;

    const auto gpu_start = std::chrono::high_resolution_clock::now();
    for(int i = 1; i <= iterations; ++i)
        cv::ocl::resize(gpu_img, gpu_resize, cv::Size(1000,1000));
    const auto gpu_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - gpu_start);

#ifdef _DEBUG
    cv::Mat difference = host_resize - gpu_resize;
    assert( cv::norm(difference, cv::NORM_INF) <= 1 );
#endif

    return Timing(cpu_time, gpu_time);
}

Timing time_warpAffine(const cv::Mat& host_img, size_t iterations)
{
    cv::Mat host_warp;

    float transform_data[] = { 2.0f  , 0.5f, -500.0f
                             , 0.333f, 3.0f, -500.0f
                             };
    cv::Mat transform(2, 3, CV_32F, transform_data);

    const auto cpu_start = std::chrono::high_resolution_clock::now();
    for(int i = 1; i <= iterations; ++i)
        cv::warpAffine(host_img, host_warp, transform, cv::Size(1000,1000));
    const auto cpu_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - cpu_start);

    const cv::ocl::oclMat gpu_img(host_img);
    cv::ocl::oclMat gpu_warp;

    const auto gpu_start = std::chrono::high_resolution_clock::now();
    for(int i = 1; i <= iterations; ++i)
        cv::ocl::warpAffine(gpu_img, gpu_warp, transform, cv::Size(1000,1000));
    const auto gpu_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - gpu_start);

#ifdef _DEBUG
    cv::Mat difference = host_warp - gpu_warp;
    assert( cv::norm(difference, cv::NORM_INF) <= 1 );
#endif

    return Timing(cpu_time, gpu_time);
}

int main()
{
    const cv::Mat img_large = cv::imread("c:/Work/CARP/test1.jpg");
    cv::Mat img;
    cv::resize(img_large, img, cv::Size(5792, 5792));   //Larger: cv::ocl::cvtColor doesn't work (returns memory garbage on some parts of the picture)

#ifdef _DEBUG
    size_t num_iterations = 1;
#else
    size_t num_iterations = 1000;
#endif
    std::cout << "Algorithm - CPU_time - GPU_time" << std::endl;

    const auto boxFilter_times = time_boxFilter (img, num_iterations);
    std::cout << " boxFilter - " << std::setw(6) << boxFilter_times.cpu.count() << "ms - " << std::setw(6) << boxFilter_times.gpu.count() << "ms" << std::endl;
    const auto integral_times  = time_integral  (img, num_iterations);
    std::cout << "  integral - " << std::setw(6) <<  integral_times.cpu.count() << "ms - " << std::setw(6) <<  integral_times.gpu.count() << "ms" << std::endl;
    const auto cvtColor_times  = time_cvtColor  (img, num_iterations);
    std::cout << "  cvtColor - " << std::setw(6) <<  cvtColor_times.cpu.count() << "ms - " << std::setw(6) <<  cvtColor_times.gpu.count() << "ms" << std::endl;
    const auto   dilate_times  = time_dilate    (img, num_iterations);
    std::cout << "    dilate - " << std::setw(6) <<    dilate_times.cpu.count() << "ms - " << std::setw(6) <<    dilate_times.gpu.count() << "ms" << std::endl;
    const auto convolve_times  = time_convolve  (img, num_iterations);
    std::cout << "  filter2D - " << std::setw(6) <<  convolve_times.cpu.count() << "ms - " << std::setw(6) <<  convolve_times.gpu.count() << "ms" << std::endl;
    const auto gaussian_times  = time_gaussian  (img, num_iterations);
    std::cout << "  gaussian - " << std::setw(6) <<  gaussian_times.cpu.count() << "ms - " << std::setw(6) <<  gaussian_times.gpu.count() << "ms" << std::endl;
    const auto   resize_times  = time_resize    (img, num_iterations);
    std::cout << "    resize - " << std::setw(6) <<    resize_times.cpu.count() << "ms - " << std::setw(6) <<    resize_times.gpu.count() << "ms" << std::endl;
    const auto   affine_times  = time_warpAffine(img, num_iterations);
    std::cout << "warpAffine - " << std::setw(6) <<    affine_times.cpu.count() << "ms - " << std::setw(6) <<    affine_times.gpu.count() << "ms" << std::endl;


    int i = 5;
}