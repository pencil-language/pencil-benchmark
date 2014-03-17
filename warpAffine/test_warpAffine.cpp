// Experimental Research Code for the CARP Project
// UjoImro, 2013

#include <chrono>
#include <opencv2/opencv.hpp>
#include <opencv2/ocl/ocl.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "opencl.hpp"
#include "utility.hpp"
#include "imgproc_warpAffine.clh"
#include "warpAffine.pencil.h"

namespace
{

template <class T0>
void convert_coeffs(T0 * M)
{
    T0 D = M[0] * M[4] - M[1] * M[3];
    D = D != 0 ? 1. / D : 0;
    T0 A11 = M[4] * D, A22 = M[0] * D;
    M[0] = A11;
    M[1] *= -D;
    M[3] *= -D;
    M[4] = A22;
    T0 b1 = -M[0] * M[2] - M[1] * M[5];
    T0 b2 = -M[3] * M[2] - M[4] * M[5];
    M[2] = b1;
    M[5] = b2;
} // convert_coeffs

} // unnamed namespace


namespace carp {

void warpAffine( carp::opencl::device & device, const cv::ocl::oclMat & src, cv::ocl::oclMat & dst, cv::Mat & transform, const cv::Size & dsize, int flags = cv::INTER_LINEAR )
{
    int interpolation = flags & cv::INTER_MAX;
    CV_Assert((src.depth() == CV_8U  || src.depth() == CV_32F) && src.oclchannels() != 2 && src.oclchannels() != 3);
    CV_Assert(interpolation == cv::INTER_NEAREST || interpolation == cv::INTER_LINEAR || interpolation == cv::INTER_CUBIC);
    dst.create(dsize, src.type());
    CV_Assert(transform.rows == 2 && transform.cols == 3);
    int warpInd = (flags & cv::WARP_INVERSE_MAP) >> 4;
    double coeffs[2][3];
    double coeffsM[2*3];
    cv::Mat coeffsMat(2, 3, CV_64F, reinterpret_cast<void *>(coeffsM));
    transform.convertTo(coeffsMat, coeffsMat.type());
    if(!warpInd) convert_coeffs(coeffsM);
    for(int i = 0; i < 2; ++i)
        for(int j = 0; j < 3; ++j)
            coeffs[i][j] = coeffsM[i*3+j];
    CV_Assert( (src.oclchannels() == dst.oclchannels()) );
    int srcStep = src.step1();
    int dstStep = dst.step1();
    float float_coeffs[2][3];
    std::string s[3] = {"NN", "Linear", "Cubic"};
    std::string kernelName = "warpAffine" + s[interpolation] + "_C1_D5";

    carp::opencl::array coeffs_cm;

    if(src.clCxt->supportsFeature(cv::ocl::Context::CL_DOUBLE))
        coeffs_cm = carp::opencl::array_<double>( device, 2 * 3, coeffs, CL_MEM_READ_ONLY );
    else
    {
        for(int m = 0; m < 2; m++)
            for(int n = 0; n < 3; n++)
                float_coeffs[m][n] = coeffs[m][n];

        coeffs_cm = carp::opencl::array_<float>( device, 2 * 3, float_coeffs, CL_MEM_READ_ONLY );
    }

    //TODO: improve this kernel
    size_t blkSizeX = 16, blkSizeY = 16;
    size_t glbSizeX;
    int cols;
    //if(src.type() == CV_8UC1 && interpolation != 2)
    if(src.type() == CV_8UC1 && interpolation != 2)
    {
        cols = (dst.cols + dst.offset % 4 + 3) / 4;
        glbSizeX = cols % blkSizeX == 0 ? cols : (cols / blkSizeX + 1) * blkSizeX;
    }
    else
    {
        cols = dst.cols;
        glbSizeX = dst.cols % blkSizeX == 0 ? dst.cols : (dst.cols / blkSizeX + 1) * blkSizeX;
    }
    size_t glbSizeY = dst.rows % blkSizeY == 0 ? dst.rows : (dst.rows / blkSizeY + 1) * blkSizeY;
    std::vector<size_t> globalThreads = {glbSizeX, glbSizeY, 1};
    std::vector<size_t> localThreads  = {blkSizeX, blkSizeY, 1};
    // device[kernelName] ( reinterpret_cast<cl_mem>(src.data), reinterpret_cast<cl_mem>(dst.data),
    //     src.cols, src.rows, dst.cols, dst.rows, srcStep, dstStep, src.offset, dst.offset,
    //                      coeffs_cm.cl(), cols ).groupsize( localThreads, globalThreads );
    device[kernelName] (
            reinterpret_cast<cl_mem>(src.data), reinterpret_cast<cl_mem>(dst.data),
            src.cols, src.rows, dst.cols, dst.rows, srcStep, dstStep, src.offset, dst.offset,
            coeffs_cm.cl(), cols ).groupsize( localThreads, globalThreads );

} // wartAffine

} // namespace carp

void time_affine( const std::vector<carp::record_t>& pool )
{
    carp::Timing timing;

    for ( auto & item : pool ) {
        PRINT(item.path());

        cv::Mat cpu_gray;
        cv::cvtColor( item.cpuimg(), cpu_gray, CV_RGB2GRAY );
        cpu_gray.convertTo( cpu_gray, CV_32F, 1.0/255. );

        std::vector<float> transform_data = { 2.0f, 0.5f, -500.0f
                                            , 0.333f, 3.0f, -500.0f
                                            };
        cv::Mat transform( 2, 3, CV_32F, transform_data.data() );

        cv::Mat cpu_result, gpu_result, pen_result;
        std::chrono::microseconds elapsed_time_cpu, elapsed_time_gpu, elapsed_time_pencil;

        {
            const auto cpu_start = std::chrono::high_resolution_clock::now();
            cv::warpAffine( cpu_gray, cpu_result, transform, cpu_gray.size() );
            const auto cpu_end = std::chrono::high_resolution_clock::now();
            elapsed_time_cpu = cpu_end - cpu_start;
        }
        {
            cv::ocl::Context * context = cv::ocl::Context::getContext();
            carp::opencl::device device(context);
            device.source_compile( imgproc_warpAffine_cl, imgproc_warpAffine_cl_len
                                 , { "warpAffineNN_C1_D0", "warpAffineLinear_C1_D0", "warpAffineCubic_C1_D0", "warpAffineNN_C4_D0"
                                   , "warpAffineLinear_C4_D0", "warpAffineCubic_C4_D0", "warpAffineNN_C1_D5", "warpAffineLinear_C1_D5"
                                   , "warpAffineCubic_C1_D5", "warpAffineNN_C4_D5", "warpAffineLinear_C4_D5", "warpAffineCubic_C4_D5"
                                   }
                                 , "-D DOUBLE_SUPPORT"
                                 );
            const auto gpu_start = std::chrono::high_resolution_clock::now();
            cv::ocl::oclMat gpu_gray(cpu_gray);
            cv::ocl::oclMat gpu_affine;
            carp::warpAffine( device, gpu_gray, gpu_affine, transform, gpu_gray.size() );
            gpu_result = gpu_affine;
            const auto gpu_end = std::chrono::high_resolution_clock::now();
            elapsed_time_gpu = gpu_end - gpu_start;
        }
        {
            // verifying the pencil code
            pen_result.create( cpu_gray.size(), CV_32F );
            convert_coeffs(reinterpret_cast<float*>(transform.data));

            const auto pencil_start = std::chrono::high_resolution_clock::now();
            pencil_affine_linear( cpu_gray.rows, cpu_gray.cols, cpu_gray.step1(), cpu_gray.ptr<float>()
                                , pen_result.rows, pen_result.cols, pen_result.step1(), pen_result.ptr<float>(),
                    transform.at<float>(0,0), transform.at<float>(0,1), transform.at<float>(1,0), transform.at<float>(1,1),
                    transform.at<float>(1,2), transform.at<float>(0,2) );
            const auto pencil_end = std::chrono::high_resolution_clock::now();
            elapsed_time_pencil = pencil_end - pencil_start;
        }
        // Verifying the results
        if ( (cv::norm(cv::abs(cpu_result - gpu_result), cv::NORM_INF ) > 1) || (cv::norm(cv::abs(cpu_result - pen_result), cv::NORM_INF ) > 1) )
        {
            PRINT(cv::norm(cv::abs(cpu_result - gpu_result), cv::NORM_INF ));
            PRINT(cv::norm(cv::abs(cpu_result - pen_result), cv::NORM_INF ));

            cv::Mat gpu_result8;
            cv::Mat cpu_result8;
            cv::Mat pen_result8;

            gpu_result.convertTo( gpu_result8, CV_8UC1, 255. );
            cpu_result.convertTo( cpu_result8, CV_8UC1, 255. );
            pen_result.convertTo( pen_result8, CV_8UC1, 255. );

            cv::imwrite( "gpu_affine.png", gpu_result8 );
            cv::imwrite( "cpu_affine.png", cpu_result8 );
            cv::imwrite( "pencil_affine.png", pen_result8 );

            throw std::runtime_error("The GPU results are not equivalent with the CPU results.");
        }

        timing.print( "affine transform", elapsed_time_cpu, elapsed_time_gpu, elapsed_time_pencil );
    }
}

int main(int argc, char* argv[])
{
#ifndef BENCHMARK_PRINT_GPU_PENCIL_SPEEDUP_ONLY
    std::cout << "This executable is iterating over all the files which are present in the directory `./pool'. " << std::endl;
#endif

    auto pool = carp::get_pool("pool");
    time_affine( pool );
    return EXIT_SUCCESS;
} // main
