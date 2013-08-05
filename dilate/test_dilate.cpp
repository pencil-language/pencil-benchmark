// UjoImro, 2013

#include <chrono>
#include <opencv2/opencv.hpp>
#include <opencv2/ocl/ocl.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "opencl.hpp"
#include "utility.hpp"
#include "filtering_morph.clh"

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

    inline void normalizeKernel(const cv::Mat &kernel, cv::ocl::oclMat &gpu_krnl, int type = CV_8U, int *nDivisor = 0, bool reverse = false)
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
        
        gpu_krnl.upload(cont_krnl);
    } // normalizeKernel
} // unnamed namespace


template<class T0>
void
time_dilate( carp::opencl::device & device, T0 & pool )
{

    double sum_quotient = 0;
    int64_t nums = 0;
    std::vector<int> elemsizes = { 5, 7, 9 /*, 11, 15, 17*/ };
    
    for ( auto & elemsize : elemsizes ) {
        PRINT(elemsize);        
        long int elapsed_time_gpu = 0;
        long int elapsed_time_cpu = 0;
        
        for ( auto & item : pool ) {
            // acquiring the image for the test
            cv::Mat cpu_gray;
            cv::cvtColor( item.cpuimg(), cpu_gray, CV_RGB2GRAY );
            
            cv::Mat host_dilate;
            cv::Mat check;
            cv::Point anchor = cv::Point( elemsize/2, elemsize/2 );
            cv::Size ksize(elemsize, elemsize);
            cv::Mat structuring_element = cv::getStructuringElement( cv::MORPH_ELLIPSE, ksize, anchor );
                        
            auto cpu_start = std::chrono::high_resolution_clock::now();
            cv::dilate( cpu_gray, host_dilate, structuring_element, anchor, /*iteration = */1, cv::BORDER_CONSTANT );        
            auto cpu_end = std::chrono::high_resolution_clock::now();
            elapsed_time_cpu += carp::microseconds(cpu_end - cpu_start);

            // gpu code
            {
                cv::ocl::oclMat src(cpu_gray);
                cv::ocl::oclMat mat_kernel;
                cv::ocl::oclMat dst;

                dst.create(src.size(), src.type());
                normalizeKernel(structuring_element, mat_kernel);
                normalizeAnchor(anchor, ksize);
                
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
                cv::ocl::Context *clCxt = src.clCxt;
                std::string kernelName;
                size_t localThreads[3] = {16, 16, 1};
                size_t globalThreads[3] = {(src.cols + localThreads[0] - 1) / localThreads[0] *localThreads[0], 
                                           (src.rows + localThreads[1] - 1) / localThreads[1] *localThreads[1], 1};
                
                if (src.type() == CV_8UC1)
                {
                    kernelName = "morph_C1_D0";
                    globalThreads[0] = ((src.cols + 3) / 4 + localThreads[0] - 1) / localThreads[0] * localThreads[0];
                    CV_Assert(localThreads[0]*localThreads[1] * 8 >= (localThreads[0] * 4 + ksize.width - 1) * (localThreads[1] + ksize.height - 1));
                }
                else
                {
                    kernelName = "morph";
                    CV_Assert(localThreads[0]*localThreads[1] * 2 >= (localThreads[0] + ksize.width - 1) * (localThreads[1] + ksize.height - 1));
                }
                
                std::string compile_option;
                
                switch (src.type())
                {
                case CV_8UC1:
                    compile_option = " -D VAL=0 ";
                    break;
                case CV_8UC3:
                case CV_8UC4:
                    compile_option = " -D VAL=0 -D GENTYPE=uchar4 ";
                    break;
                case CV_32FC1:
                    compile_option = " -D VAL=-FLT_MAX -D GENTYPE=float ";
                    break;
                case CV_32FC3:
                case CV_32FC4:
                    compile_option = " -D VAL=-FLT_MAX -D GENTYPE=float4 ";
                    break;
                default:
                    CV_Error(CV_StsUnsupportedFormat, "unsupported type");
                }
                compile_option += " -D RADIUSX=" + carp::cast<std::string>(anchor.x)
                    + " -D RADIUSY=" + carp::cast<std::string>(anchor.y)
                    + " -D LSIZE0=" + carp::cast<std::string>(localThreads[0])
                    + " -D LSIZE1=" + carp::cast<std::string>(localThreads[1])
                    +" -D DILATE ";
                
//                if (rectKernel) compiler_option += " -D RECTKERNEL ";
                auto gpu_start = std::chrono::high_resolution_clock::now();

                bool noZero = true;
                for(int i = 0; i < structuring_element.rows * structuring_element.cols; ++i)
                    if(structuring_element.data[i] != 1)
                        noZero = false;

                if (noZero) compile_option += " -D RECTKERNEL ";
                
                device.source_compile( filtering_morph_cl, filtering_morph_cl_len,  carp::string_vector(kernelName), compile_option );
                device[kernelName](
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
                    ).groupsize( { localThreads[0],  localThreads[1],  localThreads[2]}, { globalThreads[0],  globalThreads[1],  globalThreads[2]} );
                check = dst;
                auto gpu_end = std::chrono::high_resolution_clock::now();
                elapsed_time_gpu += carp::microseconds(gpu_end - gpu_start);
                device.erase(kernelName);
            } // end of gpu code

            // Verifying the results
            if ( cv::norm(host_dilate - check) > 0.01 ) {
                PRINT(cv::norm(check - host_dilate));
                cv::imwrite("cpu_dilate.png", host_dilate);
                cv::imwrite("gpu_dilate.png", check );                
                throw std::runtime_error("The GPU results are not equivalent with the CPU results.");                
            }
            
            if (elapsed_time_gpu > 1) {
                sum_quotient += static_cast<double>(elapsed_time_cpu) / elapsed_time_gpu;
                nums++;
            } // elapsed_time_gpu

        } // for pool
            
        carp::Timing::print( "dilate image", elapsed_time_cpu, elapsed_time_gpu );
            
    } // for elemsizes

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
    time_dilate( device, pool );
    return EXIT_SUCCESS;    
} // main


















// LuM end of file
