// UjoImro, 2013
// Experimental Research Code for the CARP Project

#ifndef __CARP__UTILITY__HPP__
#define __CARP__UTILITY__HPP__

#include <chrono>
#include <string>
#include <iomanip>
#include <iostream>
#include <boost/assign.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/ocl/ocl.hpp>
#include <boost/filesystem.hpp>
#include <opencv2/core/core.hpp>
#include <boost/preprocessor.hpp>

const int KiB=1024;
const int MiB=1024*KiB;
const int GiB=1024*KiB;

//#define PRINT(var)  std::cout << "debug: " << BOOST_PP_STRINGIZE(var) << " = " << var << std::endl
#define PRINT(var)

namespace carp {

    template<class T0>
    void
    print_image( cv::Mat_<T0> input, std::string name )
    {
        std::cout << name << " = [" << std::endl;

        for (int q=0; q<input.rows; q++)
        {
            std::cout << "[ ";
            for (int w=0; w<input.cols; w++)
            {
                std::cout << input(q,w);
                if (w<input.cols-1)
                    std::cout << ", ";
                else
                    std::cout << " ";
            }

            if (q<input.rows-1)
                std::cout << "], " << std::endl;
            else
                std::cout << "] " << std::endl;
        }

        std::cout << "]" << std::endl;
    } // print_image

    template<class T0, class... Types>
    std::vector<T0>
    make_vector( T0 input0, Types... inputs )
    {
        return std::vector<T0>{ input0, inputs... };
    }

    template <class T0>
    std::chrono::microseconds::rep
    microseconds( T0 t0 )
    {
        return std::chrono::duration_cast<std::chrono::microseconds>(t0).count();
    }


    class TimingLong
    {
    	std::vector<double> gpu_speedups;
    	std::vector<double> pen_speedups;
    public:
    	TimingLong() {
#ifndef BENCHMARK_PRINT_GPU_PENCIL_SPEEDUP_ONLY
    		std::cout << "    Operator - CPU Time - GPU Time - Pencil Time - GPU speedup - Pencil speedup - GPU/Pencil" << std::endl;
#endif
    	}

        void print( const std::string & name, const std::chrono::microseconds &cpu, const std::chrono::microseconds &gpu, const std::chrono::microseconds &pen ) {
            print(name, cpu.count(), gpu.count(), pen.count());
        }

        void print( const std::string & name, const long int & cpu, const long int & gpu, const long int & pen ) {
    		auto gpu_speedup = static_cast<double>(cpu)/static_cast<double>(gpu);
    		auto pen_speedup = static_cast<double>(cpu)/static_cast<double>(pen);
    		double gpu_div_pen = static_cast<double>(gpu)/static_cast<double>(pen);
    		std::cout << std::fixed << std::setprecision(3);
#ifndef BENCHMARK_PRINT_GPU_PENCIL_SPEEDUP_ONLY
    		std::cout << std::setw(12) << name << " - "
    				<< std::setw(7) << (cpu/1000000.) << "s - "
    				<< std::setw(7) << (gpu/1000000.) << "s - "
    				<< std::setw(7) << (pen/1000000.) << "s - "
    				<< std::setw(7) << gpu_speedup << 'x' << " - "
    				<< std::setw(14) << pen_speedup << 'x' << " - "
				<< std::setw(12) << gpu_div_pen  << 'x' << std::endl;
#else
    		std::cout << gpu_div_pen;
#endif
    		gpu_speedups.push_back(gpu_speedup);
    		pen_speedups.push_back(pen_speedup);
    	}

    	~TimingLong() {
#ifndef BENCHMARK_PRINT_GPU_PENCIL_SPEEDUP_ONLY
    		std::cout << "Cumulated Speed Improvement: " << std::endl
    				<< "    GPU speed / CPU speed: " << std::accumulate(gpu_speedups.begin(),gpu_speedups.end(),0.0) / gpu_speedups.size() << "x" << std::endl
    				<< "    Pen speed / CPU speed: " << std::accumulate(pen_speedups.begin(),pen_speedups.end(),0.0) / pen_speedups.size() << "x" << std::endl
    				;
#endif
    	}
    };

    class TimingShort
    {
    	std::vector<double> gpu_speedups;
    public:
    	TimingShort() {
    		std::cout << "    Operator - CPU Time - GPU Time - CPU/GPU" << std::endl;
    	}

    	void print( const std::string & name, const long int & cpu, const long int & gpu ) {
    		auto speedup = static_cast<double>(cpu)/static_cast<double>(gpu);
    		std::cout << std::fixed << std::setprecision(3);
    		std::cout << std::setw(12) << name << " - "
    				<< std::setw(7) << (cpu/1000000.) << "s - "
    				<< std::setw(7) << (gpu/1000000.) << "s - "
    				<< std::setw(7) << speedup << 'x' << std::endl;
    		gpu_speedups.push_back(speedup);
	    }

    	~TimingShort() {
    		std::cout << "Cumulated Speed Improvement: " << std::endl
    				<< "    GPU speed / CPU speed: " << std::accumulate(gpu_speedups.begin(),gpu_speedups.end(),0.0) / gpu_speedups.size() << "x" << std::endl
    				;
    	}
    };

    class record_t {
    private:
        boost::filesystem::path m_path;

    public:
        cv::Mat cpuimg() {
            cv::Mat cpuimg = cv::imread(m_path.string());
            cv::Mat img_resize = cpuimg;
            // cv::resize(cpuimg, img_resize, cv::Size(300,400));

            return img_resize;
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

    static std::map<int, std::string> borders = boost::assign::map_list_of
    ( cv::BORDER_CONSTANT, "BORDER_CONSTANT" )
    ( cv::BORDER_REPLICATE, "BORDER_REPLICATE" )
    ( cv::BORDER_REFLECT, "BORDER_REFLECT" )
    ( cv::BORDER_WRAP, "BORDER_WRAP" )
    ( cv::BORDER_REFLECT_101, "BORDER_REFLECT_101" )
    ; // borders


} // namespace carp


#endif /* __CARP__UTILITY__HPP__ */
