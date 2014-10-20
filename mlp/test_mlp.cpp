// UjoImro, 2013
// Experimental code for the CARP Project
// Copyright (c) RealEyes, 2013
// This version tests the responseMap calculation with input dumps

#include "mlp.hpp"
#include "mlp_impl_pencil.h"
#include "utility.hpp"
#include "serialization.hpp"
#include "memory.hpp"
#include "bench_mlp.hpp"

#include <tbb/parallel_for.h>
#include <tbb/concurrent_vector.h>

#include <opencv2/core/core.hpp>

#include <boost/preprocessor.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>

#include <chrono>
#include <string>
#include <iomanip>

const int processed_frames = 100;

int main()
{
    carp::conductor_t conductor;
    carp::Timing timing("MLP");
    
    // OpenCL implementation initialization
    //Load source
    std::ifstream source_file{"../mlp/mlp_impl.cl"};
    std::string source{ std::istreambuf_iterator<char>{source_file}, std::istreambuf_iterator<char>{} };

    //Compile source
    std::vector<cl::Device> devices;
    cl::Platform::getDefault().getDevices(CL_DEVICE_TYPE_GPU, &devices);
    if (devices.size()>1)
        std::cout << "warning: more then one GPU device detected, only the first will be used." << std::endl;
    cl::Device device = devices.front();
    cl::Context context = cl::Context(device);
    cl::CommandQueue queue = cl::CommandQueue(context, device);
    cl::Program program(context, source);
    try {
        program.build({device});
    } catch (const cl::Error& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device) << std::endl;
        throw;
    }
    cl::make_kernel<cl::Buffer, cl::Buffer, int, int, cl::Buffer, int, cl::LocalSpaceArg> calculateMaps_cl(program, "calculateMaps");

    static const int gangsize = 8*32;
    static const int local_memsize = 21 * 1024;
    // OpenCL implementation initialization end

    for ( conductor.importer >> BOOST_SERIALIZATION_NVP(conductor.id)
        ; (conductor.id != -1) && (conductor.id != processed_frames)
        ; conductor.importer >> BOOST_SERIALIZATION_NVP(conductor.id)
        )
    {
        //Timings
        std::chrono::duration<double> elapsed_time_cpu, elapsed_time_gpu_p_copy, elapsed_time_gpu_nocopy, elapsed_time_pencil;

        //Load next calculation package
        carp::hack_t package;
        conductor.importer >> BOOST_SERIALIZATION_NVP(package);

        // CPU implementation
        {
            gel::MLP<double> mlp;
            tbb::concurrent_vector< cv::Mat_<double> > calculatedResults(package.m_visibleLandmarks_size);
            auto start = std::chrono::high_resolution_clock::now();
            tbb::parallel_for(0,package.m_visibleLandmarks_size,[&](int q) {
                const cv::Point2i center(cvRound(package.shape(2*q,0)),cvRound(package.shape(2*q+1,0)));
                calculatedResults[q] = package.m_classifiers[q].generateResponseMap(package.alignedImage, center, package.m_mapSize);
            });
            auto end = std::chrono::high_resolution_clock::now();
            elapsed_time_cpu = end - start;

            // testing the output
            for (int q=0; q<package.m_visibleLandmarks_size; q++)
            {
                if (cv::norm( package.responseMaps[q] - calculatedResults[q] ) > 0.0001)
                    throw std::runtime_error("CPU result differs from reference result.");
            }
        }
        // OpenCL implementation
        {
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

            cl::Buffer clSelf( context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, groupsize * local_memsize * sizeof(char), buffer.data() );
            cl::Buffer clSegments( context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, segments.size() * sizeof(int), segments.data() );
            cl::Buffer clCalcpackages( context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, calcpackages.size() * sizeof(calcpackage), calcpackages.data() );
            auto start = std::chrono::high_resolution_clock::now();
            cl::EnqueueArgs enqueArgs( queue, cl::NDRange(gangsize * package.m_visibleLandmarks_size), cl::NDRange(gangsize) );
            calculateMaps_cl( enqueArgs
                            , clSelf
                            , clSegments
                            , package.m_visibleLandmarks_size
                            , package.m_mapSize
                            , clCalcpackages
                            , local_memsize
                            , cl::LocalSpaceArg{ local_memsize }
                            );
            queue.finish();
            auto end = std::chrono::high_resolution_clock::now();
            elapsed_time_gpu_nocopy = end - start;

            // copying the data back to the CPU
            std::vector<char> processed(groupsize * local_memsize);
            queue.enqueueReadBuffer(clSelf, true, 0, groupsize*local_memsize, processed.data());

            auto end_copy = std::chrono::high_resolution_clock::now();
            elapsed_time_gpu_p_copy = end_copy - start_copy;

            char * results = processed.data();

            // converting-checking the outputs
            for (int q=0; q<package.m_visibleLandmarks_size; q++)
            {
                cv::Mat_<double> nextResult = carp::convertMatFloatToCV( results + segments[q], calcpackages[q].output.responseMap );
                if (cv::norm( package.responseMaps[q] - nextResult ) > 0.0001)
                    throw std::runtime_error("OpenCL result differs from reference result.");
            }
        }
        // PENCIL implementation
        {
            // preparing the inputs
            MatChar alignedImage = carp::convertCVToMatChar(package.alignedImage);
            MatFloat shape = carp::convertCVToMatFloat(package.shape);
            mlp * m_classifiers = carp::convertHackToMlp(package);
            MatFloat * responseMaps;
            carp::allocateResponseMaps( package.m_mapSize, package.m_visibleLandmarks_size, &responseMaps );

            auto start = std::chrono::high_resolution_clock::now();
            calculateMaps( package.m_visibleLandmarks_size
                         , package.m_mapSize
                         , alignedImage
                         , shape
                         , m_classifiers
                         , &responseMaps
                         );
            auto end = std::chrono::high_resolution_clock::now();
            elapsed_time_pencil = end - start;

            // releasing the inputs
            freeMatChar(&alignedImage);
            freeMatFloat(&shape);
            carp::freeClassifiers(&m_classifiers, package.m_classifiers.size());

            // converting-testing the outputs
            for (int q=0; q<package.m_visibleLandmarks_size; q++)
            {
                cv::Mat_<double> nextResult = carp::convertMatFloatToCV( responseMaps[q] );
                if (cv::norm( package.responseMaps[q] - nextResult ) > 0.0001)
                    throw std::runtime_error("PENCIL result differs from reference result.");
            }

            // releasing the outputs
            carp::freeResponseMaps( &responseMaps, package.m_visibleLandmarks_size );
        }
        timing.print( elapsed_time_cpu, elapsed_time_gpu_p_copy, elapsed_time_gpu_nocopy, elapsed_time_pencil );
    }
    return EXIT_SUCCESS;
}