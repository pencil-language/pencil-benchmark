// Experimental Research Code for the CARP Project
// UjoImro, 2013

#include <chrono>
#include <opencv2/opencv.hpp>
#include <opencv2/ocl/ocl.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "opencl.hpp"
#include "utility.hpp"
#include "imgproc_resize.clh"
#include "resize.pencil.h"

namespace carp {

void
resize( carp::opencl::device & device,
        const cv::ocl::oclMat & src,
        cv::ocl::oclMat & dst,
        cv::Size dsize,
        int interpolation = cv::INTER_LINEAR )
{
    double fx=0.;
    double fy=0.;
    CV_Assert(src.type() == CV_8UC1 || src.type() == CV_8UC3 || src.type() == CV_8UC4
            || src.type() == CV_32FC1 || src.type() == CV_32FC3 || src.type() == CV_32FC4);
    CV_Assert(interpolation == cv::INTER_LINEAR || interpolation == cv::INTER_NEAREST);
    CV_Assert( src.size().area() > 0 );
    CV_Assert( !(dsize == cv::Size()) || (fx > 0 && fy > 0) );

    if(!(dsize == cv::Size()) && (fx > 0 && fy > 0))
    {
        if(dsize.width != (int)(src.cols * fx) || dsize.height != (int)(src.rows * fy))
        {
            CV_Error(CV_StsUnmatchedSizes, "invalid dsize and fx, fy!");
        }
    }
    if( dsize == cv::Size() )
    {
        dsize = cv::Size(cv::saturate_cast<int>(src.cols * fx), cv::saturate_cast<int>(src.rows * fy));
    }
    else
    {
        fx = (double)dsize.width / src.cols;
        fy = (double)dsize.height / src.rows;
    }

    dst.create(dsize, src.type());

    if(!( interpolation == cv::INTER_NEAREST || interpolation == cv::INTER_LINEAR) )
        CV_Error(CV_StsUnsupportedFormat, "Non-supported interpolation method");

    CV_Assert( (src.channels() == dst.channels()) );
    float ifx = 1. / fx;
    float ify = 1. / fy;
    double ifx_d = 1. / fx;
    double ify_d = 1. / fy;
    int srcStep_in_pixel = src.step1() / src.oclchannels();
    int srcoffset_in_pixel = src.offset / src.elemSize();
    int dstStep_in_pixel = dst.step1() / dst.oclchannels();
    int dstoffset_in_pixel = dst.offset / dst.elemSize();
    //printf("%d %d\n",src.step1() , dst.elemSize());
    std::string kernelName;
    if(interpolation == cv::INTER_LINEAR)
        kernelName = "resizeLN" + std::string("_C1_D0");
    else if(interpolation == cv::INTER_NEAREST)
        kernelName = "resizeNN" + std::string("_C1_D0");

    //TODO: improve this kernel
    size_t blkSizeX = 16, blkSizeY = 16;
    size_t glbSizeX;
    if(src.type() == CV_8UC1)
    {
        size_t cols = (dst.cols + dst.offset % 4 + 3) / 4;
        glbSizeX = cols % blkSizeX == 0 && cols != 0 ? cols : (cols / blkSizeX + 1) * blkSizeX;
    }
    else
    {
        glbSizeX = dst.cols % blkSizeX == 0 && dst.cols != 0 ? dst.cols : (dst.cols / blkSizeX + 1) * blkSizeX;
    }
    size_t glbSizeY = dst.rows % blkSizeY == 0 && dst.rows != 0 ? dst.rows : (dst.rows / blkSizeY + 1) * blkSizeY;

    if(interpolation == cv::INTER_NEAREST)
    {

        throw std::runtime_error("This function is not equivalent!");

        if(src.clCxt->supportsFeature(cv::ocl::Context::CL_DOUBLE))
            device[kernelName]( dst.data, src.data, dstoffset_in_pixel, srcoffset_in_pixel,
                    dstStep_in_pixel, srcStep_in_pixel, src.cols,
                    src.rows, dst.cols, dst.rows, ifx_d, ify_d ).groupsize( {blkSizeX, blkSizeY, 1}, {glbSizeX, glbSizeY, 1} );
        else /* NOT src.clCxt->supportsFeature(cv::ocl::Context::CL_DOUBLE) */
            device[kernelName]( dst.data, src.data, dstoffset_in_pixel, srcoffset_in_pixel,
                    dstStep_in_pixel, srcStep_in_pixel, src.cols,
                    src.rows, dst.cols, dst.rows, ifx, ify ).groupsize( {blkSizeX, blkSizeY, 1}, {glbSizeX, glbSizeY, 1} );
    }
    else /* NOT interpolation == INTER_NEAREST */
        device[kernelName]( reinterpret_cast<cl_mem>(dst.data), reinterpret_cast<cl_mem>(src.data),
                dstoffset_in_pixel, srcoffset_in_pixel, dstStep_in_pixel,
                srcStep_in_pixel, src.cols, src.rows, dst.cols,
                dst.rows, ifx, ify
        ).groupsize( {blkSizeX, blkSizeY, 1}, {glbSizeX, glbSizeY, 1} );
} // resize

} // namespace carp

template<class T0>
void
time_resize( T0 & pool )
{
    // std::vector<cv::Size> sizes = { {10,10}, {503, 786}, {1230, 2341}, {4243, 2324}  };
    std::vector<cv::Size> sizes = { {640, 480 } };
    std::vector<int> methods = { cv::INTER_LINEAR /*, cv::INTER_NEAREST*/ };

    carp::TimingLong timing;

    for ( auto & size : sizes ) {
        for (auto & method : methods ) {
            for ( auto & item : pool ) {
                PRINT(item.path());
                cv::Mat cpu_gray;
                cv::cvtColor( item.cpuimg(), cpu_gray, CV_RGB2GRAY );

                cv::Mat cpu_result, gpu_result, pen_result;
                std::chrono::microseconds elapsed_time_cpu, elapsed_time_gpu, elapsed_time_pencil;

                {
                    const auto cpu_start = std::chrono::high_resolution_clock::now();
                    cv::resize( cpu_gray, cpu_result, size, method );
                    const auto cpu_end = std::chrono::high_resolution_clock::now();
                    elapsed_time_cpu = cpu_end - cpu_start;
                }
                {
                    cv::ocl::Context * context = cv::ocl::Context::getContext();
                    carp::opencl::device device(context);
                    device.source_compile( imgproc_resize_cl, imgproc_resize_cl_len, { "resizeLN_C1_D0", "resizeLN_C4_D0", "resizeLN_C1_D5", "resizeLN_C4_D5", "resizeNN_C1_D0", "resizeNN_C4_D0", "resizeNN_C1_D5", "resizeNN_C4_D5" } );
                    const auto gpu_start = std::chrono::high_resolution_clock::now();
                    cv::ocl::oclMat gpu_gray(cpu_gray);
                    cv::ocl::oclMat gpu_resize;
                    carp::resize( device, gpu_gray, gpu_resize, size, method );
                    gpu_result = gpu_resize;
                    const auto gpu_end = std::chrono::high_resolution_clock::now();
                    elapsed_time_gpu = gpu_end - gpu_start;
                }
                {
                    // pencil verification
                    pen_result.create(size, CV_8UC1);

                    const auto pencil_start = std::chrono::high_resolution_clock::now();
                    pencil_resize_LN(cpu_gray.rows, cpu_gray.cols, cpu_gray.step1(), cpu_gray.ptr(), pen_result.rows, pen_result.cols, pen_result.step1(), pen_result.ptr() );
                    const auto pencil_end = std::chrono::high_resolution_clock::now();
                    elapsed_time_pencil = pencil_end - pencil_start;
                }
                // Verifying the results
                if (( cv::norm(cpu_result - gpu_result, cv::NORM_INF) > 1 ) || ( cv::norm(cpu_result - pen_result, cv::NORM_INF) > 1 ))
                {
                    PRINT(cv::norm(cpu_result - gpu_result, cv::NORM_INF));
                    PRINT(cv::norm(cpu_result - pen_result, cv::NORM_INF));

                    cv::imwrite( "gpu_resize.png", gpu_result );
                    cv::imwrite( "cpu_resize.png", cpu_result );
                    cv::imwrite( "pencil_resize.png", pen_result );

                    throw std::runtime_error("The GPU results are not equivalent with the CPU results.");
                }

                timing.print( "resize", elapsed_time_cpu, elapsed_time_gpu, elapsed_time_pencil );
            } // for pool
        }
    }
} // text_boxFilter

int main(int argc, char* argv[])
{

#ifndef BENCHMARK_PRINT_GPU_PENCIL_SPEEDUP_ONLY
    std::cout << "This executable is iterating over all the files which are present in the directory `./pool'. " << std::endl;
#endif

    auto pool = carp::get_pool("pool");
    time_resize( pool );
    return EXIT_SUCCESS;
} // main
