// UjoImro, 2013

#ifndef __OPENCL_HPP__
#define __OPENCL_HPP__

#include <string>
#include <vector>
#include <CL/cl.h>
#include <fstream>
#include <stdexcept>
#include <opencv2/core/core.hpp>
#include <memory>

#include "cltypes.h"
#include "errors.hpp"
#include "utility.hpp"

#define CL_DEVICE CL_DEVICE_TYPE_GPU

namespace carp {

namespace opencl {

namespace utility {

std::string
readfile( const std::string & filename )
{
    std::ifstream input(filename);
    std::stringstream buffer;
    buffer << input.rdbuf();
    return buffer.str();
} // readfile


void
checkerror(
        const cl_int & error,
        const std::string & filename = __FILE__,
        const int & line = __LINE__,
        const std::string & info = ""
) {
    if ( error != CL_SUCCESS )
        throw std::runtime_error( std::string("error: OpenCL: ") +
                carp::opencl::errors.at(error) +
                ";  " + info + " line: " + std::to_string(line) + " in file " + filename );
} // checkerror

// this function is equivalent to the OpenCL example in the NVidia repository
int roundup( int group_size, int global_size ) {
    if (global_size < group_size) global_size = group_size;

    int r = global_size % group_size;
    if(r == 0)
    {
        return global_size;
    } else
    {
        return global_size + group_size - r;
    }
} // roundup


template <class T0>
std::vector<T0> roundup(
        const std::vector<T0> & group_sizes,
        const std::vector<T0> & global_sizes )
        {
    if (group_sizes.size()!=global_sizes.size())
        throw std::runtime_error("The groupsize and the worksize dimensions are different!" );

    std::vector<T0> result(group_sizes.size());

    auto group = group_sizes.begin();
    auto glob = global_sizes.begin();
    auto rup = result.begin();
    for (; group!= group_sizes.end(); group++, glob++, rup++ )
        *rup = roundup( *group, *glob );

    return result;

        } // roundup for vectors

} // namespace utility



class buffer
{
private:
    size_t m_size;

public:
    buffer(size_t size) : m_size(size) { }

    size_t size() const { return m_size; }

}; // class buffer


class kernel {
private:
    std::shared_ptr<_cl_kernel> cqKernel;
    std::shared_ptr<_cl_command_queue> cqCommandQueue;
    bool m_set;

    template <class MT0>
    bool setparameter( int& pos, const MT0& mt0 ) {
        assert(cqKernel);
        utility::checkerror( clSetKernelArg( cqKernel.get(), pos, sizeof(mt0), reinterpret_cast<const void*>(&mt0) ), __FILE__, __LINE__, "(arg[" + std::to_string(pos) + "]) " );
//        utility::checkerror( clSetKernelArg( cqKernel.get(), pos, sizeof(mt0), reinterpret_cast<const void*>(&mt0) ), __FILE__, __LINE__, "(arg[" + std::to_string(pos) + "]=" + std::to_string(mt0) + ") " );
        pos++; // move the position of the parameter applied
        return true;
    } // setparameter

    bool setparameter( int& pos, const buffer& b ) {
        assert(cqKernel);
        utility::checkerror( clSetKernelArg( cqKernel.get(), pos, b.size(), nullptr ), __FILE__, __LINE__, "(arg[" + std::to_string(pos) + "]) " );
        pos++; // move the position of the parameter applied
        return true;
    } // setparameter

public:
    kernel() : m_set(false) { };

    void
    set( std::shared_ptr<_cl_kernel> cqKernel,
            std::shared_ptr<_cl_command_queue> cqCommandQueue ) {
        this->cqKernel = cqKernel;
        this->cqCommandQueue = cqCommandQueue;
        this->m_set = true;
    } // operator set

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
    template <class ...MT>
    kernel & operator() ( const MT&... mt ) {
        assert(m_set);
        int pos=0; // the position of the parameters
        bool err[] = { setparameter(pos, mt)... };
        return *this;
    } // operator ()
#pragma GCC diagnostic pop

    void
    groupsize( std::vector<size_t> groupsize, std::vector<size_t> worksize ) {
        assert(m_set);
        assert(cqKernel);
        assert(cqCommandQueue);
        std::vector<size_t> kernelsize = utility::roundup(groupsize, worksize);

#ifndef PRINT_OPENCL_PROFILING_KERNEL_EXEC_TIME
        cl_event event;

        utility::checkerror(clEnqueueNDRangeKernel( cqCommandQueue.get(), cqKernel.get(), worksize.size(), NULL, kernelsize.data(), groupsize.data(), 0, NULL, &event ), __FILE__, __LINE__ );
        clWaitForEvents(1, &event);
	clReleaseEvent(event);
#else
        cl_event event;
        long long start, end;
        double total;

        utility::checkerror(clEnqueueNDRangeKernel( cqCommandQueue.get(), cqKernel.get(), worksize.size(), NULL, kernelsize.data(), groupsize.data(), 0, NULL, &event), __FILE__, __LINE__ );
        clWaitForEvents(1, &event);
        clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof start, &start, NULL);
        clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof end, &end, NULL);
        total = (double)(end - start) / 1e9;
        printf("[RealEyes]   Kernel execution time in seconds (without data copy and kernel compilation): %5.9f\n", total);
#endif

        utility::checkerror(clFinish(cqCommandQueue.get()), __FILE__, __LINE__ );
    } // groupsize
}; // class kernel

class device {
private:
    typedef std::map<std::string, opencl::kernel> kernels_t;

    std::shared_ptr<_cl_context> cxGPUContext;         // OpenCL context
    std::shared_ptr<_cl_command_queue> cqCommandQueue; // OpenCL command queue
    std::shared_ptr<cl_platform_id> cpPlatform;       // OpenCL platform
    std::vector<cl_device_id> devices;
    // cl_device_id cdDevice;          // OpenCL device
    std::shared_ptr<_cl_program> cpProgram;           // OpenCL program
    kernels_t kernels; // OpenCL kernel
    cl_uint num_platforms;          // The number of OpenCL platforms
    cl_uint num_devices;            // The number of OpenCL devices

    device( const device & ) = delete; // device is not copyable

public:

    device() {
        cl_int err; // for the error messages

        cpPlatform.reset(new cl_platform_id);
        utility::checkerror(clGetPlatformIDs( 1, cpPlatform.get(), &num_platforms ), __FILE__, __LINE__ );
        assert(num_platforms==1);  // there is only one supported platform at the time

        // getting the number of devices
        utility::checkerror(clGetDeviceIDs( *cpPlatform, CL_DEVICE, 0, NULL, &num_devices), __FILE__, __LINE__ );
        assert(num_devices>0);

        devices.resize(num_devices);
        utility::checkerror(clGetDeviceIDs( *cpPlatform, CL_DEVICE, devices.size(), devices.data(), NULL ), __FILE__, __LINE__ );
        if (num_devices>1)
            std::cout << "warning: more then one GPU device detected, only the first will be used." << std::endl;

        cl_context tmp_cxGPUContext = clCreateContext( NULL, devices.size(), devices.data(), NULL, NULL, &err );
        utility::checkerror( err, __FILE__, __LINE__ );
        cxGPUContext.reset( tmp_cxGPUContext, clReleaseContext );

        cl_command_queue tmp_cqCommandQueue;
#ifndef PRINT_OPENCL_PROFILING_KERNEL_EXEC_TIME
        tmp_cqCommandQueue = clCreateCommandQueue( cxGPUContext.get(), devices[0], 0,  &err ); // kernels are executed in order
#else
        tmp_cqCommandQueue = clCreateCommandQueue( cxGPUContext.get(), devices[0], CL_QUEUE_PROFILING_ENABLE,  &err ); // kernels are executed in order
#endif
        utility::checkerror( err, __FILE__, __LINE__ );
        cqCommandQueue.reset( tmp_cqCommandQueue, clReleaseCommandQueue );

    } // device

    ~device() {
        // if(cpProgram) clReleaseProgram(cpProgram);
        // if (automanage) {
        //     if(cqCommandQueue) clReleaseCommandQueue(cqCommandQueue);
        //     if(cxGPUContext) clReleaseContext(cxGPUContext);
        // } // automanage


    } // ~device

    void compile( const std::vector<std::string> & filenames,
            const std::vector<std::string> & kernelnames,
            const std::string options = std::string("") ) {
        cl_int err;

        // loading the sources
        std::vector<std::string> sources;
        std::vector<const char*> c_strs;
        std::vector<size_t> lengths;
        const char * coptions = NULL;
        if (options!="") coptions = options.c_str();

        for( const std::string & filename : filenames ) {
            sources.push_back(utility::readfile(filename));
            c_strs.push_back(sources.back().c_str());
            lengths.push_back(sources.back().size());
        }

        cl_program tmp_cpProgram = clCreateProgramWithSource( cxGPUContext.get(), sources.size(), &(c_strs[0]), lengths.data(), &err );
        utility::checkerror(err, __FILE__, __LINE__ );
        cpProgram.reset( tmp_cpProgram, clReleaseProgram );

        // building the OpenCL source code
        err = clBuildProgram( cpProgram.get(), devices.size(), devices.data(), coptions, NULL, NULL );
        if ( err!=CL_SUCCESS )
        {
            size_t len;
            std::vector<char> buffer(1<<20);

            utility::checkerror( clGetProgramBuildInfo( cpProgram.get(), devices[0], CL_PROGRAM_BUILD_LOG, buffer.size(), buffer.data(), &len ), __FILE__, __LINE__ );

            std::cerr << buffer.data() << std::endl;
            utility::checkerror(err, __FILE__, __LINE__ );
            throw std::runtime_error("error: OpenCL: The compilation has failed for unknown reason.");
        }

        // extracting the kernel entrances
        for( const std::string & kernel_name : kernelnames )
        {
            cl_int err;
            cl_kernel tmp_kernel = clCreateKernel( cpProgram.get(), kernel_name.c_str(), &err );
            utility::checkerror( err, __FILE__, __LINE__ );
            std::shared_ptr<_cl_kernel> btmp_kernel( tmp_kernel, clReleaseKernel );
            kernels[kernel_name].set( btmp_kernel, cqCommandQueue );
        }
    }

    void source_compile( const unsigned char * source,
            const unsigned int source_len,
            const std::vector<std::string> & kernelnames,
            const std::string options = std::string("") ) {
        cl_int err;
        size_t cl_code_len = source_len;
        const char * csource = reinterpret_cast<const char*>(source);
        const char * coptions = NULL;
        if (options!="") coptions = options.c_str();

        cl_program tmp_cpProgram = clCreateProgramWithSource( cxGPUContext.get(), 1, &csource, &cl_code_len, &err );
        utility::checkerror(err, __FILE__, __LINE__ );
        cpProgram.reset( tmp_cpProgram, clReleaseProgram );

        // building the OpenCL source code
        err = clBuildProgram( cpProgram.get(), 1, devices.data(), coptions, NULL, NULL );
        if ( err!=CL_SUCCESS )
        {
            size_t len;
            std::vector<char> buffer(1<<20);

            utility::checkerror( clGetProgramBuildInfo( cpProgram.get(), devices[0], CL_PROGRAM_BUILD_LOG, buffer.size(), buffer.data(), &len ), __FILE__, __LINE__ );

            std::cerr << buffer.data() << std::endl;
            utility::checkerror(err, __FILE__, __LINE__ );
            throw std::runtime_error("error: OpenCL: The compilation has failed for unknown reason.");
        }

        // extracting the kernel entrances
        for( const std::string & kernel_name : kernelnames )
        {
            cl_int err;
            cl_kernel tmp_kernel;
            tmp_kernel = clCreateKernel( cpProgram.get(), kernel_name.c_str(), &err );
            utility::checkerror( err, __FILE__, __LINE__ );
            std::shared_ptr<_cl_kernel> btmp_kernel( tmp_kernel, clReleaseKernel );
            kernels[kernel_name].set( btmp_kernel, cqCommandQueue );
        }
    }


    kernel & operator [] ( const std::string & kernel_name ) {
        return kernels[kernel_name];
    }

    cl_context get_context() {
        return cxGPUContext.get();
    }

    cl_command_queue get_queue() {
        return cqCommandQueue.get();
    }
    
    cl_device_id get_device_id() const {
        return devices[0];
    }

    void erase( const std::string & kernel_name ) {
        kernels.erase(kernel_name);
    }
}; // class device

template <class T0>
class array_ {
private:
    std::shared_ptr<_cl_mem> cl_ptr;
    size_t m_size;
    cl_context cqContext;
    cl_command_queue cqCommandQueue;

    array_( const array_<T0> & ) = delete; // copy constructor forbidden

public:
    array_( opencl::device & device, size_t size, T0 * ptr ) 
    {
        this->m_size = size;
        this->cqContext = device.get_context();
        this->cqCommandQueue = device.get_queue();

        cl_int err;
        cl_mem tmp;
        tmp = clCreateBuffer( cqContext, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, size * sizeof(T0), reinterpret_cast<void*>(ptr), &err );
        utility::checkerror( err, __FILE__, __LINE__ );
        cl_ptr.reset(tmp, clReleaseMemObject);
    }

    cl_mem cl() {
        return cl_ptr.get();
    }
    
    std::vector<T0> get( ) {
        std::vector<T0> result(m_size);
        utility::checkerror( clEnqueueReadBuffer( cqCommandQueue, cl_ptr.get(), CL_TRUE, 0, sizeof(T0) * m_size, result.data(), 0, NULL, NULL), __FILE__, __LINE__ );

        return result;
    } // extract

}; // class array_

} // namespace opencl

} // namespace carp



#endif /* __OPENCL_HPP__ */
// LuM end of file
