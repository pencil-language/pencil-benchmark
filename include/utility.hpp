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
    	std::vector<double> gpu_speedups;
    	std::vector<double> pen_speedups;
    	std::vector<double> gpu_div_pens;
    public:
    	Timing() {
#ifndef BENCHMARK_PRINT_GPU_PENCIL_SPEEDUP_ONLY
    		std::cout << "    Operator - CPU Time - GPU Time - Pencil Time - CPU/GPU - CPU/Pencil - GPU/Pencil" << std::endl;
#endif
    	}

        void print( const std::string & name, const std::chrono::duration<double> &cpu, const std::chrono::duration<double> &gpu, const std::chrono::duration<double> &pen ) {
    		auto gpu_speedup = cpu / gpu;
    		auto pen_speedup = cpu / pen;
    		auto gpu_div_pen = gpu / pen;
    		std::cout << std::fixed << std::setprecision(3);
#ifndef BENCHMARK_PRINT_GPU_PENCIL_SPEEDUP_ONLY
    		std::cout << std::setw(12) << name << " - "
    		          << std::setw(7) << cpu.count() << "s - "
    		          << std::setw(7) << gpu.count() << "s - "
    		          << std::setw(7) << pen.count() << "s - "
    		          << std::setw(7) << gpu_speedup << "x - "
    		          << std::setw(7) << pen_speedup << "x - "
    		          << std::setw(7) << gpu_div_pen << 'x' << std::endl;
#else
    		std::cout << gpu_div_pen;
#endif
    		gpu_speedups.push_back(gpu_speedup);
    		pen_speedups.push_back(pen_speedup);
    		gpu_div_pens.push_back(gpu_div_pen);
    	}

        ~Timing() {
#ifndef BENCHMARK_PRINT_GPU_PENCIL_SPEEDUP_ONLY
            std::cout << "Average speed ratios: " << std::endl
                    << "    GPU speed / CPU speed: " << std::accumulate(gpu_speedups.begin(),gpu_speedups.end(),0.0) / gpu_speedups.size() << "x\n"
                    << "    Pen speed / CPU speed: " << std::accumulate(pen_speedups.begin(),pen_speedups.end(),0.0) / pen_speedups.size() << "x\n"
                    << "    Pen speed / GPU speed: " << std::accumulate(gpu_div_pens.begin(),gpu_div_pens.end(),0.0) / gpu_div_pens.size() << "x" << std::endl
                    ;
#endif
        }
    };

    class record_t {
    private:
        boost::filesystem::path m_path;

    public:
        cv::Mat cpuimg() {
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


    template <class RT, class T0>
    RT cast( const T0 & from )
    {
        std::stringstream stringstream;
        RT to;
        stringstream << from;
        stringstream >> to;
        return to;
    } // gel_cast

    static std::map<int, std::string> borders{{ cv::BORDER_CONSTANT   , "BORDER_CONSTANT"    }
                                             ,{ cv::BORDER_REPLICATE  , "BORDER_REPLICATE"   }
                                             ,{ cv::BORDER_REFLECT    , "BORDER_REFLECT"     }
                                             ,{ cv::BORDER_WRAP       , "BORDER_WRAP"        }
                                             ,{ cv::BORDER_REFLECT_101, "BORDER_REFLECT_101" }
                                             };
} // namespace carp


#endif /* __CARP__UTILITY__HPP__ */
