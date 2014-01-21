// Experimental Research Code for the CARP Project
// UjoImro, 2013

#include <chrono>
#include <opencv2/opencv.hpp>
#include <opencv2/ocl/ocl.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "opencl.hpp"
#include "utility.hpp"
#include "imgproc_warpAffine.clh"
#include "affine.pencil.h"

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

            carp::opencl::array coeffs_cm = carp::opencl::array_<float>( device, 2 * 3, float_coeffs, CL_MEM_READ_ONLY );
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

template<class T0>
void
time_affine( carp::opencl::device & device, T0 & pool )
{
    carp::TimingLong timing;

    for ( auto & item : pool ) {
        PRINT(item.path());

        cv::Mat cpu_gray;

        cv::cvtColor( item.cpuimg(), cpu_gray, CV_RGB2GRAY );
        cpu_gray.convertTo( cpu_gray, CV_32F, 1.0/255. );
        cv::Mat host_affine;
        
        std::vector<float> transform_data = { 2.0f, 0.5f, -500.0f,
                                              0.333f, 3.0f, -500.0f };
        // std::vector<float> transform_data = { 1.0f, 0.5f, 10.0f,
        //                                       0.0f, 1.0f, -3.0f };
        
        cv::Mat transform( 2, 3, CV_32F, transform_data.data() );
            
        const auto cpu_start = std::chrono::high_resolution_clock::now();
        cv::warpAffine( cpu_gray, host_affine, transform, cpu_gray.size() );
        const auto cpu_end = std::chrono::high_resolution_clock::now();
        auto elapsed_time_cpu = carp::microseconds(cpu_end - cpu_start);

	const auto gpu_start = std::chrono::high_resolution_clock::now();
	cv::ocl::oclMat gpu_gray(cpu_gray);
	cv::ocl::oclMat gpu_affine;
	carp::warpAffine( device, gpu_gray, gpu_affine, transform, gpu_gray.size() );
	cv::Mat check(gpu_affine);
	const auto gpu_end = std::chrono::high_resolution_clock::now();
	auto elapsed_time_gpu = carp::microseconds(gpu_end - gpu_start);

        // verifying the pencil code
        cv::Mat pencil_affine( cpu_gray.size(), CV_32F );
        convert_coeffs(reinterpret_cast<float*>(transform.data));

        const auto pencil_start = std::chrono::high_resolution_clock::now();
        pencil_affine_linear(
            cpu_gray.rows, cpu_gray.cols, cpu_gray.step1(), cpu_gray.ptr<float>(),
            pencil_affine.rows, pencil_affine.cols, pencil_affine.step1(), pencil_affine.ptr<float>(),
            transform.at<float>(0,0), transform.at<float>(0,1), transform.at<float>(1,0), transform.at<float>(1,1),
            transform.at<float>(1,2), transform.at<float>(0,2) );
        const auto pencil_end = std::chrono::high_resolution_clock::now();
        auto elapsed_time_pencil = carp::microseconds(pencil_end - pencil_start);
        
        // Verifying the results
        if ( (cv::norm(cv::abs(host_affine - check), cv::NORM_INF ) > 1) ||
             (cv::norm(cv::abs(pencil_affine - host_affine), cv::NORM_INF ) > 1) ) {
            PRINT(cv::norm(cv::abs(host_affine - check), cv::NORM_INF ));
            PRINT(cv::norm(cv::abs(pencil_affine - host_affine), cv::NORM_INF ));
        
            cv::Mat check8;
            cv::Mat host_affine8;
            cv::Mat pencil_affine8;            
            
            check.convertTo( check8, CV_8UC1, 255. );
            host_affine.convertTo( host_affine8, CV_8UC1, 255. );
            pencil_affine.convertTo( pencil_affine8, CV_8UC1, 255. );
            
            cv::imwrite( "gpu_affine.png", check8 );
            cv::imwrite( "cpu_affine.png", host_affine8 );
            cv::imwrite( "pencil_affine.png", pencil_affine8 );

            throw std::runtime_error("The GPU results are not equivalent with the CPU results.");
        }
                        
        timing.print( "affine transform", elapsed_time_cpu, elapsed_time_gpu, elapsed_time_pencil );
    }
}

int main(int argc, char* argv[])
{
    std::cout << "This executable is iterating over all the files which are present in the directory `./pool'. " << std::endl;    

    auto pool = carp::get_pool("pool");

    // Initializing OpenCL
    cv::ocl::Context * context = cv::ocl::Context::getContext();
    carp::opencl::device device(context);
    device.source_compile( imgproc_warpAffine_cl, imgproc_warpAffine_cl_len
			 , { "warpAffineNN_C1_D0", "warpAffineLinear_C1_D0", "warpAffineCubic_C1_D0", "warpAffineNN_C4_D0", "warpAffineLinear_C4_D0", "warpAffineCubic_C4_D0"
			   , "warpAffineNN_C1_D5", "warpAffineLinear_C1_D5", "warpAffineCubic_C1_D5", "warpAffineNN_C4_D5", "warpAffineLinear_C4_D5", "warpAffineCubic_C4_D5"
			   }, "-D DOUBLE_SUPPORT" );
    
    time_affine( device, pool );
    return EXIT_SUCCESS;    
} // main


















// LuM end of file
