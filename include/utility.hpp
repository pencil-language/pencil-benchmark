// UjoImro, 2013
// Experimental Research Code for the CARP Project

#ifndef __CARP__UTILITY__HPP__
#define __CARP__UTILITY__HPP__

#include <chrono>
#include <string>
#include <iomanip>
#include <iostream>

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/ocl/ocl.hpp>

#include <boost/filesystem.hpp>

//#define PRINT(var)  std::cout << "debug: " << BOOST_PP_STRINGIZE(var) << " = " << var << std::endl
#define PRINT(var)

namespace carp {

class Timing
{
    std::vector<double> cpu_timings;
    std::vector<double> gpu_timings;
    std::vector<double> pen_timings;
    std::vector<double> cpu_div_gpus;
    std::vector<double> cpu_div_pens;
    std::vector<double> gpu_div_pens;
public:
    Timing(const std::string & name) {
#ifndef BENCHMARK_PRINT_GPU_PENCIL_SPEEDUP_ONLY
        std::cout << "Measuring performance of " << name << std::endl;
        std::cout << " CPU Time - GPU Time - PEN Time -  CPU/GPU -  CPU/PEN -  GPU/PEN" << std::endl;
#endif
    }

    void print( const std::chrono::duration<double> &cpu, const std::chrono::duration<double> &gpu, const std::chrono::duration<double> &pen ) {
        auto cpu_div_gpu = cpu / gpu;
        auto cpu_div_pen = cpu / pen;
        auto gpu_div_pen = gpu / pen;
        std::cout << std::fixed << std::setprecision(3);
#ifndef BENCHMARK_PRINT_GPU_PENCIL_SPEEDUP_ONLY
        std::cout << std::setw(8) << cpu.count() << "s -";
        std::cout << std::setw(8) << gpu.count() << "s -";
        std::cout << std::setw(8) << pen.count() << "s -";
        std::cout << std::setw(8) << cpu_div_gpu << "x -";
        std::cout << std::setw(8) << cpu_div_pen << "x -";
        std::cout << std::setw(8) << gpu_div_pen << 'x' << std::endl;
#else
        std::cout << gpu_div_pen;
#endif
        cpu_timings.push_back(cpu.count());
        gpu_timings.push_back(gpu.count());
        pen_timings.push_back(pen.count());
        cpu_div_gpus.push_back(cpu_div_gpu);
        cpu_div_pens.push_back(cpu_div_pen);
        gpu_div_pens.push_back(gpu_div_pen);
    }

    ~Timing() {
#ifndef BENCHMARK_PRINT_GPU_PENCIL_SPEEDUP_ONLY
        std::cout << "    Total CPU time: " << std::accumulate(cpu_timings.begin(),cpu_timings.end(),0.0) << "s\n";
        std::cout << "    Total GPU time: " << std::accumulate(gpu_timings.begin(),gpu_timings.end(),0.0) << "s\n";
        std::cout << "    Total Pen time: " << std::accumulate(pen_timings.begin(),pen_timings.end(),0.0) << "s\n";
        std::cout << "    AVG GPU / CPU speed ratio: " << std::accumulate(cpu_div_gpus.begin(),cpu_div_gpus.end(),0.0) / cpu_div_gpus.size() << "x\n";
        std::cout << "    AVG Pen / CPU speed ratio: " << std::accumulate(cpu_div_pens.begin(),cpu_div_pens.end(),0.0) / cpu_div_pens.size() << "x\n";
        std::cout << "    AVG Pen / GPU speed ratio: " << std::accumulate(gpu_div_pens.begin(),gpu_div_pens.end(),0.0) / gpu_div_pens.size() << "x" << std::endl;
#endif
    }
};

class record_t {
private:
    boost::filesystem::path m_path;

public:
    cv::Mat cpuimg() const {
        return cv::imread(m_path.string());
    }

    std::string path() const {
        return m_path.string();
    }

    record_t( const boost::filesystem::path & path )
    : m_path(path)  { }

}; // record_t


template <class T0>
std::vector<record_t>
get_pool( T0 pathname )
{
    boost::filesystem::path path(pathname);

    if ( (!boost::filesystem::exists(path)) or (!boost::filesystem::is_directory(path)) )
        throw std::runtime_error( std::string("Directory `") + pathname + "' does not exists. The directory should contain the testing images.");

    std::vector<record_t> pool;

    boost::filesystem::directory_iterator end_iter;
    for ( boost::filesystem::directory_iterator iter(path); iter!= end_iter; ++iter )
    {
        if (boost::filesystem::is_regular_file(iter->status()))
        {
            std::string extension = iter->path().extension().string();
            if ( (extension ==".jpg") or (extension==".jpeg") ) {
                PRINT(iter->path().string());
                pool.push_back(record_t(iter->path()));
            }
        }
    }

    return pool;
} // get_pool

static std::map<int, std::string> borders
{{ cv::BORDER_CONSTANT   , "BORDER_CONSTANT"    }
,{ cv::BORDER_REPLICATE  , "BORDER_REPLICATE"   }
,{ cv::BORDER_REFLECT    , "BORDER_REFLECT"     }
,{ cv::BORDER_WRAP       , "BORDER_WRAP"        }
,{ cv::BORDER_REFLECT_101, "BORDER_REFLECT_101" }
};

} // namespace carp

#endif /* __CARP__UTILITY__HPP__ */
