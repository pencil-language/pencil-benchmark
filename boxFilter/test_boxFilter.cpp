// UjoImro, 2013

#include <chrono>
#include <opencv2/opencv.hpp>
#include <opencv2/ocl/ocl.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "opencl.hpp"
#include "utility.hpp"
#include "filtering_boxFilter.clh"

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

} // unnamed namespace

    

template<class T0>
void
time_boxFilter( carp::opencl::device & device, T0 & pool )
{
    long int elapsed_time_gpu = 0;
    long int elapsed_time_cpu = 0;
    
    std::vector<int> ksizes = { 11,15,21,25,31,35,41 };
    // std::vector<int> ksizes = {31};

    for ( auto & item : pool ) {
        PRINT(item.path);        
        for ( auto size : ksizes)
        {
            PRINT(size);
            cv::Mat host_filtered;
            cv::Mat cpu_gray;
            cv::Mat check;    

            // the ksize must be 2k+1 form and anchor must be (2k+1)<<1
            cv::Size  ksize(size,size);
            cv::Point anchor((size>>1),(size>>1));
        
            cv::cvtColor( item.cpuimg, cpu_gray, CV_RGB2GRAY );
            cv::boxFilter( cpu_gray, host_filtered, -1, ksize, anchor, true, cv::BORDER_CONSTANT );
                        
            cv::ocl::oclMat gpu_result;
            cv::ocl::oclMat gpu_gray(cpu_gray);
            gpu_result.create( item.cpuimg.rows, item.cpuimg.cols, CV_8U );
            float alpha = ksize.height * ksize.width;
            size_t blockSizeX = 256, blockSizeY = 1;
            size_t gSize = blockSizeX - (ksize.width - 1);
            size_t threads = (gpu_result.offset % gpu_result.step % 4 + gpu_result.cols + 3) / 4;
            size_t globalSizeX = threads % gSize == 0 ? threads / gSize * blockSizeX : (threads / gSize + 1) * blockSizeX;
            size_t globalSizeY = ((gpu_result.rows + 1) / 2) % blockSizeY == 0 ? ((gpu_result.rows + 1) / 2) : (((gpu_result.rows + 1) / 2) / blockSizeY + 1) * blockSizeY;

            PRINT(blockSizeX);
            PRINT(blockSizeY);
            PRINT(globalSizeX);
            PRINT(globalSizeY);
                        
            size_t globalThreads[3] = { globalSizeX, globalSizeY, 1 };
            size_t localThreads[3]  = { blockSizeX, blockSizeY, 1 };
            normalizeAnchor(anchor, ksize);

            device.source_compile(
                filtering_boxFilter_cl,
                filtering_boxFilter_cl_len,
                carp::string_vector("boxFilter_C1_D0"),
                "   -D anX=" + carp::cast<std::string>(anchor.x)
                + " -D anY=" + carp::cast<std::string>(anchor.y)
                + " -D ksX=" + carp::cast<std::string>(ksize.width)
                + " -D ksY=" + carp::cast<std::string>(ksize.height)
                + " -D BORDER_CONSTANT" );
            
            device["boxFilter_C1_D0"] (
                reinterpret_cast<cl_mem>(gpu_gray.data),
                reinterpret_cast<cl_mem>(gpu_result.data),
                static_cast<cl_float>(alpha),
                gpu_gray.offset,
                gpu_gray.wholerows,
                gpu_gray.wholecols,
                static_cast<int>(gpu_gray.step),
                gpu_result.offset,
                gpu_result.rows,
                gpu_result.cols,
                static_cast<int>(gpu_result.step)
                ).groupsize( { blockSizeX, blockSizeY, 1 }, { globalSizeX, globalSizeY, 1 } );
            // //.groupsize( localThreads, globalThreads );

            device.erase("boxFilter_C1_D0");
            check = gpu_result;
            PRINT(gpu_result.type());            
            PRINT(gpu_result.rows);
            PRINT(gpu_result.cols);
            PRINT(gpu_result.step);
            
            // cv::imwrite( "host_filtered.png", host_filtered );
            // cv::imwrite( "gpu_filtered.png", check );
            PRINT(check.type());
            PRINT(check.rows);
            PRINT(check.cols);
            PRINT(check.step);

            PRINT(host_filtered.type());
            PRINT(host_filtered.rows);
            PRINT(host_filtered.cols);
            PRINT(host_filtered.step);
            
            PRINT( cv::norm(check - host_filtered ));
            
        } // for size
    } // for item
    
    
    
    
    // for(int i = 0; i < iterations; ++i) {
    //     PRINT(i);        
    //     for ( auto & record : pool ) {
    //         PRINT(record.path);            
    //         // CPU Bench
    //         {
    //             auto start = std::chrono::high_resolution_clock::now();
    //             cv::cvtColor( record.cpuimg, cpu_gray, CV_RGB2GRAY );
    //             auto end = std::chrono::high_resolution_clock::now();
    //             elapsed_time_cpu += carp::microseconds(end - start);
    //         }
    //         // GPU Bench
    //         {
    //             cv::ocl::oclMat gpu_gray;
    //             cv::ocl::oclMat gpuimg(record.cpuimg);
    //             gpu_gray.create( gpuimg.rows, gpuimg.cols, CV_8U );
                                
    //             // int code = CV_BGR2GRAY;                
    //             // int bidx = (code == CV_BGR2GRAY || code == CV_BGRA2GRAY) ? 0 : 2;
    //             int bidx = 2;
    //             auto start = std::chrono::high_resolution_clock::now();
    //             device["RGB2Gray"] (
    //                 gpuimg.cols,
    //                 gpuimg.rows ,
    //                 static_cast<int>(gpuimg.step),
    //                 static_cast<int>(gpu_gray.step),
    //                 gpuimg.channels()+1,
    //                 bidx,
    //                 reinterpret_cast<cl_mem>(gpuimg.data),
    //                 reinterpret_cast<cl_mem>(gpu_gray.data) 
    //                 )
    //                 .groupsize( carp::make_vector<size_t>(16,16), carp::make_vector<size_t>(record.cpuimg.cols,record.cpuimg.rows) );               
    //             auto end = std::chrono::high_resolution_clock::now();
    //             elapsed_time_gpu += carp::microseconds(end - start);
    //             check = gpu_gray;
    //             //check = gpuimg;
                
    //         }
    //         // Verifying the results
    //         if ( cv::norm(check - cpu_gray) > 0.01 ) {
    //             cv::imwrite( "gpu_img.png", check );
    //             cv::imwrite( "cpu_img.png", cpu_gray );
    //             throw std::runtime_error("The GPU results are not equivalent with the CPU results.");                
    //         }
    //     }
    // }
    
    // carp::Timing::print( "cvtColor", elapsed_time_cpu, elapsed_time_gpu );
    return;
}


int main(int argc, char* argv[])
{

    std::cout << "This executable is iterating over all the files which are present in the directory `./pool'. " << std::endl;    

    auto pool = carp::get_pool("pool");

    // Initializing OpenCL
    cv::ocl::Context * context = cv::ocl::Context::getContext();
    carp::Timing::printHeader();
    carp::opencl::device device(context);
    time_boxFilter( device, pool );

    return EXIT_SUCCESS;    
} // main


















// LuM end of file
