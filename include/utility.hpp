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
    std::vector<double> gpu_p_copy_timings;
    std::vector<double> gpu_nocopy_timings;
    std::vector<double> pen_timings;
public:
    Timing(const std::string & name) {
        std::cout << "Measuring performance of " << name << std::endl;
        std::cout << " CPU Time -  GPU+copy -  GPU Time -  PEN Time" << std::endl;
    }

    void print( const std::chrono::duration<double> &cpu, const std::chrono::duration<double> &gpu_p_copy, const std::chrono::duration<double> &gpu_nocopy, const std::chrono::duration<double> &pen ) {
        std::cout << std::fixed << std::setprecision(6);
        std::cout << std::setw(8) << cpu       .count() << "s - ";
        std::cout << std::setw(8) << gpu_p_copy.count() << "s - ";
        std::cout << std::setw(8) << gpu_nocopy.count() << "s - ";
        std::cout << std::setw(8) << pen       .count() << "s" << std::endl;

	cpu_timings       .push_back(cpu       .count());
        gpu_p_copy_timings.push_back(gpu_p_copy.count());
        gpu_nocopy_timings.push_back(gpu_nocopy.count());
        pen_timings       .push_back(pen       .count());
    }

    ~Timing() {
        std::cout << "    Total CPU time           : " << std::accumulate(cpu_timings       .begin(),cpu_timings       .end(),0.0) << "\n";
        std::cout << "    Total GPU time (inc copy): " << std::accumulate(gpu_p_copy_timings.begin(),gpu_p_copy_timings.end(),0.0) << "\n";
        std::cout << "    Total GPU time (w/o copy): " << std::accumulate(gpu_nocopy_timings.begin(),gpu_nocopy_timings.end(),0.0) << "\n";
        std::cout << "    Total Pen time           : " << std::accumulate(pen_timings       .begin(),pen_timings       .end(),0.0) << "\n";
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
