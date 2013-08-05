// UjoImro, 2013

#include <chrono>
#include <opencv2/opencv.hpp>
#include <opencv2/ocl/ocl.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "opencl.hpp"
#include "utility.hpp"
#include "imgproc_integral_sum.clh"
    

template<class T0>
void
time_integral( carp::opencl::device & device, T0 & pool, int iteration )
{

    double sum_quotient = 0;
    int64_t nums = 0;
    
    
    for ( int q=0; q<iteration; q++ ) {
        
        for ( auto & item : pool ) {
            PRINT(item.path);

            long int elapsed_time_gpu = 0;
            long int elapsed_time_cpu = 0;
        
            cv::Mat host_integral;
            cv::Mat cpu_gray;
            cv::Mat check;    

            cv::cvtColor( item.cpuimg, cpu_gray, CV_RGB2GRAY );
            if(cpu_gray.cols * cpu_gray.rows <= 2901 * 2901)
                cv::resize(cpu_gray, cpu_gray, cv::Size(2901,2901));        
            auto cpu_start = std::chrono::high_resolution_clock::now();
            cv::integral( cpu_gray, host_integral, CV_32SC1 );        
            auto cpu_end = std::chrono::high_resolution_clock::now();
            elapsed_time_cpu += carp::microseconds(cpu_end - cpu_start);
        
            {
                cv::ocl::oclMat src(cpu_gray);
                cv::ocl::oclMat sum;            
            
                int vlen = 4;
                int offset = src.offset / vlen;
                int pre_invalid = src.offset % vlen;
                int vcols = (pre_invalid + src.cols + vlen - 1) / vlen;
        
                cv::ocl::oclMat t_sum;
                int w = src.cols + 1, h = src.rows + 1;
                int depth;
                std::string integral_sum_cols;
                std::string integral_sum_rows;            
                // if(src.cols * src.rows <= 2901 * 2901)
                {
                    t_sum.create(src.cols, src.rows, CV_32SC1);
                    sum.create(h, w, CV_32SC1);
                    integral_sum_rows = "integral_sum_rows_D4";
                    integral_sum_cols = "integral_sum_cols_D4";
                }
                // else
                // {
                //     t_sum.create(src.cols, src.rows, CV_32FC1);
                //     sum.create(h, w, CV_32FC1);
                //     integral_sum_rows = "integral_sum_rows_D5";
                //     integral_sum_cols = "integral_sum_cols_D5";                
                // }

                // {
                //     t_sum.create(src.cols, src.rows, CV_64FC1);
                //     sum.create(h, w, CV_64FC1);
                //     integral_sum_rows = "integral_sum_rows_D6";
                //     integral_sum_cols = "integral_sum_cols_D6";                
                // }
                depth = sum.depth();
            
                int sum_offset = sum.offset / vlen;
                auto gpu_start = std::chrono::high_resolution_clock::now();
                device[integral_sum_cols](
                    reinterpret_cast<cl_mem>(src.data),
                    reinterpret_cast<cl_mem>(t_sum.data),
                    offset,
                    pre_invalid,
                    src.rows,
                    src.cols,
                    static_cast<int>(src.step),
                    static_cast<int>(t_sum.step)
                    ).groupsize( {256, 1, 1}, {((vcols + 1) / 2) * 256, 1, 1} );
                device[integral_sum_rows](
                    reinterpret_cast<cl_mem>(t_sum.data),
                    reinterpret_cast<cl_mem>(sum.data),
                    t_sum.rows,
                    t_sum.cols,
                    static_cast<int>(t_sum.step),
                    static_cast<int>(sum.step),
                    sum_offset
                    ).groupsize( {256, 1, 1}, {t_sum.cols  * 32, 1, 1} );

                auto gpu_end = std::chrono::high_resolution_clock::now();
                elapsed_time_gpu += carp::microseconds(gpu_end - gpu_start);

                check = sum;
            }
        
            // Verifying the results
            if ( cv::norm(host_integral - check) > 0.01 ) {
                PRINT(cv::norm(check - host_integral));
                // no use to write out the results, as they are in float
                throw std::runtime_error("The GPU results are not equivalent with the CPU results.");                
            }

            if (elapsed_time_gpu > 1) {                
                sum_quotient += elapsed_time_cpu / elapsed_time_gpu;
                nums++;                
            }
                        
            carp::Timing::print( "integral image", elapsed_time_cpu, elapsed_time_gpu );

        } // for pool
            
    } // for q 

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
    device.source_compile( imgproc_integral_sum_cl, imgproc_integral_sum_cl_len,
                           carp::string_vector("integral_sum_cols_D4", "integral_sum_rows_D4", "integral_sum_cols_D5", "integral_sum_rows_D5", "integral_sum_cols_D6", "integral_sum_rows_D6" ),
                           " -D DOUBLE_SUPPORT" );
    time_integral( device, pool, 10 );
    return EXIT_SUCCESS;    
} // main


















// LuM end of file
