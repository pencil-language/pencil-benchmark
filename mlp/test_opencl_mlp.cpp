// UjoImro, 2013
// Experimental code for the CARP Project
// Copyright (c) RealEyes, 2013
// This version tests the responseMap calculation with input dumps

#include <vector>
#include <iomanip>

#include "memory.hpp"
#include "bench_mlp.hpp"

#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>

const int processed_frames = 100;
const int local_memsize = 21 * 1024;

int main()
{
    std::cout << "MLP benchmark using OpenCL:" << std::endl;

    carp::conductor_t conductor; // the class for importing the input from the clm

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
    cl::make_kernel<cl::Buffer, cl::Buffer, int, int, cl::Buffer, int, cl::LocalSpaceArg> calculateMaps(program, "calculateMaps");

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
        
        cl::Buffer clSelf( context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, groupsize * local_memsize * sizeof(char), buffer.data() );
        cl::Buffer clSegments( context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, segments.size() * sizeof(int), segments.data() );
        cl::Buffer clCalcpackages( context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, calcpackages.size() * sizeof(calcpackage), calcpackages.data() );
        auto start = std::chrono::high_resolution_clock::now();
        cl::EnqueueArgs enqueArgs( queue, cl::NDRange(gangsize * package.m_visibleLandmarks_size), cl::NDRange(gangsize) );
        calculateMaps( enqueArgs
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
        elapsed_time += (end - start);

        // copying the data back to the CPU
        std::vector<char> processed(groupsize * local_memsize);
        queue.enqueueReadBuffer(clSelf, true, 0, groupsize*local_memsize, processed.data());

        auto end_copy = std::chrono::high_resolution_clock::now();
        elapsed_time_copy += (end_copy - start_copy);

        char * results = processed.data();

        // converting-checking the outputs
        for (int q=0; q<package.m_visibleLandmarks_size; q++)
        {
            cv::Mat_<double> nextResult = carp::convertMatFloatToCV( results + segments[q], calcpackages[q].output.responseMap );
            if (cv::norm( package.responseMaps[q] - nextResult ) > 0.0001) throw std::runtime_error("package.responseMaps[q] - calculatedResults[q] ) < 0.0001 failed");
        }
    } // packages
    std::cout << "total elapsed time (w/o copy) = " << elapsed_time.count() << " s." << std::endl;
    std::cout << "processing speed   (w/o copy) = " << processed_frames / elapsed_time.count() << "fps" << std::endl;

    std::cout << "total elapsed time (inc copy) = " << elapsed_time_copy.count() << " s." << std::endl;
    std::cout << "processing speed   (inc copy) = " << processed_frames / elapsed_time_copy.count() << "fps" << std::endl;
    return EXIT_SUCCESS;
}
