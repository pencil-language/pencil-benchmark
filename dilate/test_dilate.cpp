// UjoImro, 2013

#include <chrono>
#include <opencv2/opencv.hpp>
#include <opencv2/ocl/ocl.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "opencl.hpp"
#include "utility.hpp"
#include "filtering_morph.clh"
#include "dilate.pencil.h"

namespace {
    inline void normalizeAnchor(int & anchor, int ksize)
    {
        if (anchor < 0)
        {
            anchor = ksize >> 1;
        }

        CV_Assert(0 <= anchor && anchor < ksize);
    } // normalizeAnchor


    inline void normalizeAnchor(cv::Point & anchor, const cv::Size &ksize)
    {
        normalizeAnchor(anchor.x, ksize.width);
        normalizeAnchor(anchor.y, ksize.height);
    } // normalizeAnchor

    inline cv::Mat normalizeKernel(const cv::Mat &kernel, int type = CV_8U, int *nDivisor = 0, bool reverse = false)
    {
        int scale = nDivisor && (kernel.depth() == CV_32F || kernel.depth() == CV_64F) ? 256 : 1;
                if (nDivisor)
        {
            *nDivisor = scale;
        }

        cv::Mat temp(kernel.size(), type);
        kernel.convertTo(temp, type, scale);
        cv::Mat cont_krnl = temp.reshape(1, 1);

        if (reverse)
        {
            int count = cont_krnl.cols >> 1;

            for (int i = 0; i < count; ++i)
            {
                std::swap(cont_krnl.at<int>(0, i), cont_krnl.at<int>(0, cont_krnl.cols - 1 - i));
            }
        }

        return cont_krnl;
    } // normalizeKernel
} // unnamed namespace


namespace carp {
    cv::ocl::oclMat dilate( carp::opencl::device & device, cv::ocl::oclMat src, cv::ocl::oclMat mat_kernel, cv::Point anchor, cv::Size ksize ) {

        cv::ocl::oclMat dst(src.size(), src.type());

        //! data type supported: CV_8UC1, CV_8UC4, CV_32FC1, CV_32FC4
        //Normalize the result by default
        //float alpha = ksize.height * ksize.width;
        CV_Assert(src.clCxt == dst.clCxt);
        CV_Assert((src.cols == dst.cols) &&
                  (src.rows == dst.rows));
        CV_Assert((src.oclchannels() == dst.oclchannels()));

        int srcStep = src.step1() / src.oclchannels();
        int dstStep = dst.step1() / dst.oclchannels();
        int srcOffset = src.offset /  src.elemSize();
        int dstOffset = dst.offset /  dst.elemSize();

        int srcOffset_x = srcOffset % srcStep;
        int srcOffset_y = srcOffset / srcStep;
        std::string kernelName;
        std::vector<size_t> localThreads = {16, 16, 1};
        std::vector<size_t> globalThreads = { (src.cols + localThreads[0] - 1) / localThreads[0] *localThreads[0]
                                            , (src.rows + localThreads[1] - 1) / localThreads[1] *localThreads[1]
                                            , 1
                                            };

        assert(src.type() == CV_8UC1);
        globalThreads[0] = ((src.cols + 3) / 4 + localThreads[0] - 1) / localThreads[0] * localThreads[0];
        CV_Assert(localThreads[0]*localThreads[1] * 8 >= (localThreads[0] * 4 + ksize.width - 1) * (localThreads[1] + ksize.height - 1));

        device["morph_C1_D0"](
            reinterpret_cast<cl_mem>(src.data),
            reinterpret_cast<cl_mem>(dst.data),
            srcOffset_x,
            srcOffset_y,
            src.cols,
            src.rows,
            srcStep,
            dstStep,
            reinterpret_cast<cl_mem>(mat_kernel.data),
            src.wholecols,
            src.wholerows,
            dstOffset
            ).groupsize( localThreads, globalThreads );
        return dst;
    } // dilate

    void dilate( carp::opencl::device & device, cv::Mat cpu_gray, cv::Mat & result, cv::Mat structuring_element, cv::Point anchor, int border_type, cv::Size ksize ) {
    }

} // namespace carp

void time_dilate( const std::vector<carp::record_t>& pool, const std::vector<int>& elemsizes, int iteration )
{
    carp::Timing timing("dilate image");

    for ( int q=0; q<iteration; q++ ) {
        for ( auto & item : pool ) {
            PRINT(item.path());
            for ( auto & elemsize : elemsizes ) {
                PRINT(elemsize);

                // acquiring the image for the test
                cv::Mat cpu_gray;
                cv::cvtColor( item.cpuimg(), cpu_gray, CV_RGB2GRAY );

                cv::Point anchor( elemsize/2, elemsize/2 );
                cv::Size ksize(elemsize, elemsize);
                cv::Mat structuring_element = cv::getStructuringElement( cv::MORPH_ELLIPSE, ksize, anchor );

                cv::Mat cpu_result, gpu_result, pen_result;
                std::chrono::duration<double> elapsed_time_cpu, elapsed_time_gpu_p_copy, elapsed_time_gpu_nocopy, elapsed_time_pencil;

                {
                    const auto cpu_start = std::chrono::high_resolution_clock::now();
                    cv::dilate( cpu_gray, cpu_result, structuring_element, anchor, /*iteration = */1, cv::BORDER_CONSTANT );
                    const auto cpu_end = std::chrono::high_resolution_clock::now();
                    elapsed_time_cpu = cpu_end - cpu_start;
                }
                {
                    cv::ocl::Context * context = cv::ocl::Context::getContext();
                    carp::opencl::device device(context);
                    auto normalizedKernel = normalizeKernel(structuring_element);
                    normalizeAnchor(anchor, ksize);
                    std::string compile_option = " -D RADIUSX=" + std::to_string(anchor.x)
                                               + " -D RADIUSY=" + std::to_string(anchor.y)
                                               + " -D LSIZE0=16 -D LSIZE1=16 -D DILATE -D VAL=0";
                    const auto gpu_start_compile = std::chrono::high_resolution_clock::now();
                    device.source_compile( filtering_morph_cl, filtering_morph_cl_len, {"morph_C1_D0"}, compile_option );

                    const auto gpu_start_copy = std::chrono::high_resolution_clock::now();
                    cv::ocl::oclMat src(cpu_gray);
                    cv::ocl::oclMat mat_kernel(normalizedKernel);
                    const auto gpu_start = std::chrono::high_resolution_clock::now();
                    auto result = carp::dilate( device, src, mat_kernel, anchor, ksize );
                    const auto gpu_end = std::chrono::high_resolution_clock::now();
                    gpu_result = result;
                    const auto gpu_end_copy = std::chrono::high_resolution_clock::now();
                    device.erase("morph_C1_D0");
                    const auto gpu_end_compile = std::chrono::high_resolution_clock::now();
                    elapsed_time_gpu_p_copy = gpu_end_copy - gpu_start_copy;
                    elapsed_time_gpu_nocopy = gpu_end      - gpu_start;
                    auto elapsed_time_gpu_compile = gpu_start_compile - gpu_end_compile;
                }
                {
                    pen_result = cv::Mat(cpu_gray.size(), CV_8U);

                    const auto pencil_start = std::chrono::high_resolution_clock::now();
                    pencil_dilate( cpu_gray.rows, cpu_gray.cols, cpu_gray.step1(), cpu_gray.ptr()
                                 , pen_result.step1(), pen_result.ptr()
                                 , structuring_element.rows, structuring_element.cols, structuring_element.step1(), structuring_element.ptr()
                                 , anchor.x, anchor.y
                                 );
                    const auto pencil_end = std::chrono::high_resolution_clock::now();
                    elapsed_time_pencil = pencil_end - pencil_start;
                }
                // Verifying the results
                if ( (cv::norm(cpu_result - gpu_result) > 0.01) || (cv::norm(cpu_result - pen_result) > 0.01) ) {
                    PRINT(cv::norm(gpu_result - cpu_result));
                    PRINT(cv::norm(cpu_result - pen_result));
                    cv::imwrite("cpu_dilate.png", cpu_result);
                    cv::imwrite("gpu_dilate.png", gpu_result );
                    cv::imwrite("pencil_dilate.png", pen_result );
                    throw std::runtime_error("The GPU results are not equivalent with the CPU results.");
                }
                timing.print( elapsed_time_cpu, elapsed_time_gpu_p_copy, elapsed_time_gpu_nocopy, elapsed_time_pencil );
            }
        }
    }
}

int main(int argc, char* argv[])
{

    std::cout << "This executable is iterating over all the files which are present in the directory `./pool'. " << std::endl;

    auto pool = carp::get_pool("pool");

#ifdef RUN_ONLY_ONE_EXPERIMENT
    time_dilate( pool, { 7 }, 1 );
#else
    time_dilate( pool, { 3, 5, 7, 9 }, 6 );
#endif

    return EXIT_SUCCESS;
}
