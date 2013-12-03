// Experimental Research Code for the CARP Project
// UjoImro, 2013

#include <chrono>
#include <opencv2/opencv.hpp>
#include <opencv2/ocl/ocl.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "opencl.hpp"
#include "utility.hpp"
#include "filter_sep_col.clh"
#include "filter_sep_row.clh"
#include "gaussian.pencil.h"

namespace
{
    inline void normalizeAnchor(int &anchor, int ksize)
    {
        if (anchor < 0)
        {
            anchor = ksize >> 1;
        }

        CV_Assert(0 <= anchor && anchor < ksize);
    }

    inline void normalizeAnchor(cv::Point &anchor, const cv::Size &ksize)
    {
        normalizeAnchor(anchor.x, ksize.width);
        normalizeAnchor(anchor.y, ksize.height);
    }

    inline void normalizeROI( cv::Rect &roi, const cv::Size &ksize, const cv::Point &anchor, const cv::Size &src_size)
    {
        if (roi == cv::Rect(0, 0, -1, -1))
        {
            roi = cv::Rect(0, 0, src_size.width, src_size.height);
        }

        CV_Assert(ksize.height > 0 && ksize.width > 0 && ((ksize.height & 1) == 1) && ((ksize.width & 1) == 1));
        CV_Assert((anchor.x == -1 && anchor.y == -1) || (anchor.x == ksize.width >> 1 && anchor.y == ksize.height >> 1));
        CV_Assert(roi.x >= 0 && roi.y >= 0 && roi.width <= src_size.width && roi.height <= src_size.height);
    }

} // unnamed namespace 

namespace carp {
    void gaussian( carp::opencl::device & device, const cv::ocl::oclMat & src, cv::ocl::oclMat & dst, const cv::Size & ksize, const double & gaussX, const double & gaussY, int bordertype ) {
        dst.create(src.size(), src.type());
        cv::Point anchor(-1, -1);
        normalizeAnchor(anchor, ksize);
        cv::Size src_size = src.size();
        int cn = src.oclchannels();

        cv::ocl::oclMat dstBuf;
        cv::ocl::oclMat srcROI;
        cv::ocl::oclMat dstROI;
        cv::ocl::oclMat dstBufROI;
        dstBuf.create(src_size.height + ksize.height - 1, src_size.width, CV_MAKETYPE(CV_32F, cn));
        cv::Rect roi = cv::Rect(0, 0, -1, -1);
                
        normalizeROI(roi, ksize, anchor, src_size);

        srcROI = src(roi);
        dstROI = dst(roi);

        int depth = dst.depth();
        CV_Assert(ksize.width > 0 && ksize.width % 2 == 1 && ksize.height > 0 && ksize.height % 2 == 1);

        cv::Mat kx = cv::getGaussianKernel(ksize.width, gaussX, std::max(depth, CV_32F));
        cv::Mat ky;

        if (ksize.height == ksize.width && std::abs(gaussX - gaussY) < DBL_EPSILON)
        {
            ky = kx;
        }
        else
        {
            ky = cv::getGaussianKernel(ksize.height, gaussY, std::max(depth, CV_32F));
        }
                                                
        cv::ocl::oclMat oclKX(kx.reshape(1,1));
        cv::ocl::oclMat oclKY(ky.reshape(1,1));
                
        int sdepth   = srcROI.type();                
        int ddepth   = dstBuf.type();
                
        int bdepth   = std::max(std::max(sdepth, ddepth), CV_32F);
        int bufType  = CV_MAKETYPE(bdepth, cn);
        int channels = srcROI.oclchannels();
                                                                                
        std::vector<size_t> localThreads = {16, 16, 1};

        std::string compile_option;
        compile_option = " -D RADIUSX=" + carp::cast<std::string>(anchor.x)
            + " -D LSIZE0=" + carp::cast<std::string>(localThreads[0])
            + " -D LSIZE1=" + carp::cast<std::string>(localThreads[1])
            + " -D CN=" + carp::cast<std::string>(channels)
            + " -D " + carp::borders[bordertype];

        std::vector<size_t> globalThreads(3);
        globalThreads[1] = (dstBuf.rows + localThreads[1] - 1) / localThreads[1] * localThreads[1];
        globalThreads[2] = (1 + localThreads[2] - 1) / localThreads[2] * localThreads[2];
        globalThreads[0] = (dstBuf.cols + localThreads[0] - 1) / localThreads[0] * localThreads[0];

        //sanity checks
        CV_Assert(srcROI.cols == dstBuf.cols);
        CV_Assert(srcROI.oclchannels() == dstBuf.oclchannels());
        CV_Assert(ksize.width == (anchor.x << 1) + 1);
        int src_pix_per_row, dst_pix_per_row;
        int src_offset_x, src_offset_y;//, dst_offset_in_pixel;
        src_pix_per_row = srcROI.step / srcROI.elemSize();
        src_offset_x = (srcROI.offset % srcROI.step) / srcROI.elemSize();
        src_offset_y = srcROI.offset / srcROI.step;
        dst_pix_per_row = dstBuf.step / dstBuf.elemSize();
        //dst_offset_in_pixel = dstBuf.offset / dstBuf.elemSize();
        int ridusy = (dstBuf.rows - srcROI.rows) >> 1;
        device.source_compile( filter_sep_row_cl, filter_sep_row_cl_len, {"row_filter_C1_D5"}, compile_option );
                
        device["row_filter_C1_D5"]( reinterpret_cast<cl_mem>(srcROI.data), reinterpret_cast<cl_mem>(dstBuf.data), dstBuf.cols, dstBuf.rows,
            srcROI.wholecols, srcROI.wholerows, src_pix_per_row, src_offset_x, src_offset_y, dst_pix_per_row, ridusy, reinterpret_cast<cl_mem>(oclKX.data)
            ).groupsize( localThreads, globalThreads );
                                                                                                                                
        localThreads = {16, 16, 1};
        compile_option = "";
                
        globalThreads.resize(3);
        globalThreads[1] = (dstROI.rows + localThreads[1] - 1) / localThreads[1] * localThreads[1];
        globalThreads[2] = (1 + localThreads[2] - 1) / localThreads[2] * localThreads[2];
        globalThreads[0] = (dstROI.cols + localThreads[0] - 1) / localThreads[0] * localThreads[0];
        compile_option = " -D RADIUSY=" + carp::cast<std::string>(anchor.y)
            + " -D LSIZE0=" + carp::cast<std::string>(localThreads[0])
            + " -D LSIZE1=" + carp::cast<std::string>(localThreads[1]) 
            + " -D CN=" + carp::cast<std::string>(channels)
            + " -D " + carp::borders[bordertype]
            + " -D GENTYPE_SRC=" + "float"
            + " -D GENTYPE_DST=" + "float"
            + " -D convert_to_DST=" + "";                        

        //sanity checks
        CV_Assert(dstBuf.cols == dstROI.cols);
        CV_Assert(dstBuf.oclchannels() == dstROI.oclchannels());
        CV_Assert(ksize.height == (anchor.y << 1) + 1);
        int dst_offset_in_pixel;
        src_pix_per_row = dstBuf.step / dstBuf.elemSize();
        //src_offset_x = (dstBuf.offset % dstBuf.step) / dstBuf.elemSize();
        //src_offset_y = dstBuf.offset / dstBuf.step;
        dst_pix_per_row = dstROI.step / dstROI.elemSize();
        dst_offset_in_pixel = dstROI.offset / dstROI.elemSize();
                
        device.source_compile( filter_sep_col_cl, filter_sep_col_cl_len, { "col_filter" }, compile_option );
        device["col_filter"]( reinterpret_cast<cl_mem>(dstBuf.data), reinterpret_cast<cl_mem>(dstROI.data), dstROI.cols, dstROI.rows,
            dstBuf.wholecols, dstBuf.wholerows, src_pix_per_row, dst_pix_per_row, dst_offset_in_pixel, reinterpret_cast<cl_mem>(oclKY.data)
            ).groupsize(localThreads, globalThreads);

        device.erase("col_filter");
    } // carp::gaussian
    
} // namespace carp 


template<class T0>
void
time_gaussian( carp::opencl::device & device, T0 & pool )
{

    double cpu_gpu_quotient=0;
    double pencil_gpu_quotient=0;
    double pencil_cpu_quotient=0;

    int64_t nums = 0;
    std::vector<int> sizes = {5, 9, 11, 25, 41};

    for ( auto & size : sizes ) {
        PRINT(size);
        cv::Size ksize(size, size+4);
        
        double gaussX = 7.;
        double gaussY = 9.;
        gaussX = std::max(gaussX, 0.0);
        gaussY = std::max(gaussY, 0.0);

        int bordertype = cv::BORDER_REPLICATE;
        
        for ( auto & item : pool ) {

            long int elapsed_time_gpu = 0;
            long int elapsed_time_cpu = 0;
            long int elapsed_time_pencil = 0;
            
            cv::Mat cpu_gray;
            cv::Mat check;    

            cv::cvtColor( item.cpuimg(), cpu_gray, CV_RGB2GRAY );
            cpu_gray.convertTo( cpu_gray, CV_32F, 1.0/255. );
            cv::Mat host_blur;
            
            const auto cpu_start = std::chrono::high_resolution_clock::now();
            cv::GaussianBlur( cpu_gray, host_blur, ksize, gaussX, gaussY, bordertype );
            const auto cpu_end = std::chrono::high_resolution_clock::now();
            elapsed_time_cpu += carp::microseconds(cpu_end - cpu_start);

            cv::ocl::oclMat dst;
            cv::ocl::oclMat src(cpu_gray);                
            const auto gpu_start = std::chrono::high_resolution_clock::now();
            carp::gaussian(device, src, dst, ksize, gaussX, gaussY, bordertype );            
            const auto gpu_end = std::chrono::high_resolution_clock::now();
            elapsed_time_gpu += carp::microseconds(gpu_end - gpu_start);                            
            check = dst;                

            cv::Mat kernel = cv::getGaussianKernel(ksize.width, gaussX, CV_32F).t();
            cv::Mat t_kernel = cv::getGaussianKernel(ksize.height, gaussY, CV_32F);
            
            cv::Mat pencil;            
            {
                cv::Mat intermediate;
                intermediate.create( cpu_gray.size(), CV_32F );

                pencil.create( cpu_gray.size(), CV_32F );
                
                const auto pencil_start = std::chrono::high_resolution_clock::now();
                pencil_gaussian(
                    cpu_gray.rows, cpu_gray.cols, cpu_gray.step1(), cpu_gray.ptr<float>(),
                    kernel.rows, kernel.cols, kernel.step1(), kernel.ptr<float>(),
                    t_kernel.rows, t_kernel.cols, t_kernel.step1(), t_kernel.ptr<float>(),
                    intermediate.step1(), intermediate.ptr<float>(),
                    pencil.step1(), pencil.ptr<float>()
                    );
                const auto pencil_end = std::chrono::high_resolution_clock::now();
                elapsed_time_pencil += carp::microseconds(pencil_end - pencil_start);

            }
            
            
            // Verifying the results
            if ( (cv::norm( host_blur - check ) > 0.01) ||
                 (cv::norm( host_blur - pencil) > 0.01) ) {
                PRINT(cv::norm(check - host_blur ));
                PRINT(cv::norm(pencil - host_blur));
                
                cv::Mat check8;
                cv::Mat host_blur8;
                cv::Mat pencil8;
                cv::Mat diff8;
                                
                check.convertTo( check8, CV_8UC1, 255. );
                host_blur.convertTo( host_blur8, CV_8UC1, 255. );
                pencil.convertTo( pencil8, CV_8UC1, 255. );
                cv::Mat absdiff = cv::abs(pencil - host_blur);
                absdiff.convertTo( diff8, CV_8UC1, 255. );                
                
                cv::imwrite( "gpu_gaussian.png", check8 );
                cv::imwrite( "cpu_gaussian.png", host_blur8 );
                cv::imwrite( "pencil_gaussian.png", pencil8 );
                cv::imwrite( "diff_gaussian.png", diff8 );
                
                throw std::runtime_error("The GPU results are not equivalent with the CPU results.");
            }

            if (elapsed_time_gpu > 1) {
                cpu_gpu_quotient += static_cast<double>(elapsed_time_cpu) / elapsed_time_gpu;
                pencil_gpu_quotient += static_cast<double>(elapsed_time_pencil) / elapsed_time_gpu;
                pencil_cpu_quotient += static_cast<double>(elapsed_time_pencil) / elapsed_time_cpu;
                nums++;
            }
                        
            carp::Timing::print( "blur image", elapsed_time_cpu, elapsed_time_gpu, elapsed_time_pencil );

        } // for pool
    } // sizes
    
    carp::Timing::CSI( cpu_gpu_quotient, pencil_gpu_quotient, pencil_cpu_quotient, nums );    



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
//    device.source_compile( imgproc_convolve_cl, imgproc_convolve_cl_len,
//                           carp::string_vector("convolve_D5" ) );
    time_gaussian( device, pool );
    return EXIT_SUCCESS;    
} // main


















// LuM end of file
