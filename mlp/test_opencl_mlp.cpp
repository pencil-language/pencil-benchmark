// UjoImro, 2013
// Experimental code for the CARP Project
// Copyright (c) RealEyes, 2013
// This version tests the responseMap calculation with input dumps

#include <vector>
#include <iomanip>

#include "opencl.hpp"
#include "memory.hpp"
#include "bench_mlp.hpp"

// OpenCL includes
#include "mlp_impl.clh"

#include <CL/cl.hpp>

const int processed_frames = 100;
const int local_memsize = 21 * 1024;

int main()
{
    std::cout << "MLP benchmark using OpenCL:" << std::endl;

    carp::conductor_t conductor; // the class for importing the input from the clm
    carp::opencl::device device;
    device.source_compile( mlp_impl_cl, mlp_impl_cl_len, {"calculateMaps"} );

    int gangsize = 8*32;
    std::chrono::duration<double> elapsed_time(0), elapsed_time_copy(0);

    for ( conductor.importer >> BOOST_SERIALIZATION_NVP(conductor.id)
        ; (conductor.id != -1) && (conductor.id != processed_frames)
        ; conductor.importer >> BOOST_SERIALIZATION_NVP(conductor.id)
    )
    {
        carp::hack_t package;

        conductor.importer >> BOOST_SERIALIZATION_NVP(package);

        int groupsize = package.m_visibleLandmarks_size;
        std::vector<char> buffer(groupsize * local_memsize);
        std::vector<carp::memory::allocator> pools( groupsize, carp::memory::allocator(local_memsize));
        std::vector<int> segments(groupsize);
        {
            for (int q = 0; q < groupsize; q++ )
                segments[q] = local_memsize * q;
        }
        
        auto calcpackages = convertHackToMlp( buffer.data(), pools, segments, package );

        // preparing the opencl data
        auto start_copy = std::chrono::high_resolution_clock::now();
        carp::opencl::array_<char> clSelf( device, groupsize * local_memsize, buffer.data() );
        carp::opencl::array_<int> clSegments( device, segments.size(), segments.data() );
        carp::opencl::array_<calcpackage> clCalcpackages( device, calcpackages.size(), calcpackages.data() );

        auto start = std::chrono::high_resolution_clock::now();
        device["calculateMaps"]( clSelf.cl()
                , clSegments.cl()
                , package.m_visibleLandmarks_size
                , package.m_mapSize
                , clCalcpackages.cl()
                , local_memsize
                , carp::opencl::buffer(local_memsize)
        ).groupsize( {gangsize}, {gangsize * package.m_visibleLandmarks_size} );
        auto end = std::chrono::high_resolution_clock::now();
        elapsed_time += (end - start);

        // copying the data back to the CPU
        auto processed = clSelf.get();

        auto end_copy = std::chrono::high_resolution_clock::now();
        elapsed_time_copy += (end_copy - start_copy);

        char * results = reinterpret_cast<char*>(processed.data());

        // converting the outputs
        std::vector< cv::Mat_<double> > calculatedResults;
        for (int q=0; q<package.m_visibleLandmarks_size; q++)
        {
            cv::Mat_<double> nextResult;
            nextResult = carp::convertMatFloatToCV( results + segments[q], calcpackages[q].output.responseMap );
            calculatedResults.push_back(nextResult);
        }

        // testing the output
        for (int q=0; q<package.m_visibleLandmarks_size; q++)
        {
            if (cv::norm( package.responseMaps[q] - calculatedResults[q] ) > 0.0001) throw std::runtime_error("package.responseMaps[q] - calculatedResults[q] ) < 0.0001 failed");
        }
    } // packages
    std::cout << "total elapsed time (w/o copy) = " << elapsed_time.count() << " s." << std::endl;
    std::cout << "processing speed   (w/o copy) = " << processed_frames / elapsed_time.count() << "fps" << std::endl;

    std::cout << "total elapsed time (inc copy) = " << elapsed_time_copy.count() << " s." << std::endl;
    std::cout << "processing speed   (inc copy) = " << processed_frames / elapsed_time_copy.count() << "fps" << std::endl;
    return EXIT_SUCCESS;
}
