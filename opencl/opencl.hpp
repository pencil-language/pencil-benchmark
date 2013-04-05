// UjoImro, 2013

#ifndef __OPENCL_HPP__
#define __OPENCL_HPP__

#include <string>
#include <vector>
#include <CL/cl.h>
#include <fstream>
#include <stdexcept>
#include <boost/foreach.hpp>
#include <boost/smart_ptr.hpp>

#include "errors.hpp"

namespace carp {
    namespace opencl {
        

        std::string
        readfile( const std::string & filename )
        {
            std::ifstream input( filename, std::ifstream::binary );
    
            if (!input) throw std::runtime_error("error: Couldn't open file `" + filename + "'" );
        
            input.seekg(0, input.end);
            int length = input.tellg();
            input.seekg(0, input.beg);

            boost::scoped_array<char> result_c_str(new char[length+1]);            
            input.read(result_c_str.get(), length);
            if (!input) throw std::runtime_error("error: Could not read file `" + filename + "'");
            
            std::string result(result_c_str.get());

            return result;
        } // readfile


        class kernel {
        private:
            cl_kernel cqKernel;
            
        public:
            kernel() : cqKernel(NULL) { };

            kernel & operator= ( cl_kernel cqKernel ) {
                this->cqKernel = cqKernel;                
                return *this;
            } // operator=

            void release() {
                if (cqKernel) clReleaseKernel(cqKernel);
                cqKernel=NULL;                
            } // release
            
        }; // class kernel
        

        class device {
        private:
            typedef std::map<std::string, opencl::kernel> kernels_t;            
            
            cl_context cxGPUContext;        // OpenCL context
            cl_command_queue cqCommandQueue;// OpenCL command queue
            cl_platform_id cpPlatform;      // OpenCL platform
            cl_device_id cdDevice;          // OpenCL device
            cl_program cpProgram;           // OpenCL program
            kernels_t kernels; // OpenCL kernel            
            cl_uint num_platforms;          // The number of OpenCL platforms
            cl_uint num_devices;            // The number of OpenCL devices
            
        public:

            device() {
                cl_int err; // for the error messages

                device::check(clGetPlatformIDs(1, &cpPlatform, &num_platforms), __FILE__, __LINE__ );
                assert(num_platforms==1);  // there is only one supported platform at the time

                device::check(clGetDeviceIDs(cpPlatform, CL_DEVICE_TYPE_GPU, 1, &cdDevice, &num_devices), __FILE__, __LINE__ );
                assert(num_devices==1);

                cxGPUContext = clCreateContext( NULL, 1, &cdDevice, NULL, NULL, &err );
                device::check( err, __FILE__, __LINE__ );
                
                cqCommandQueue = clCreateCommandQueue( cxGPUContext, cdDevice, CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE, &err );
                device::check( err, __FILE__, __LINE__ );
                
            }
            
            static void check(
                const cl_int & error,
                const std::string & filename = __FILE__,
                const int & line = __LINE__
                ) throw (std::exception&) {
                if ( error != CL_SUCCESS )
                    throw std::runtime_error( std::string("error: OpenCL: ") +
                                              carp::opencl::errors[error] +
                                              "; line: " + std::to_string(__LINE__) + " in file " + __FILE__ );
            }
    
            ~device() {
                if(cpProgram) clReleaseProgram(cpProgram);
                if(cqCommandQueue) clReleaseCommandQueue(cqCommandQueue);
                if(cxGPUContext) clReleaseContext(cxGPUContext);
                std::for_each( kernels.begin(), kernels.end(), [] ( kernels_t::allocator_type::value_type & kernel ) {// kernel->second.release();
                    } );
            } // ~device

            void compile( const std::vector<std::string> & filenames,
                          const std::vector<std::string> & kernelnames ) {
                cl_int err;
                
                // loading the sources
                std::vector<std::string> sources;
                std::vector<const char*> c_strs;
                std::vector<size_t> lengths;                
                
                BOOST_FOREACH( std::string filename, filenames ) {                    
                    sources.push_back(readfile(filename));
                    c_strs.push_back(sources.back().c_str());
                    lengths.push_back(sources.back().size());                    
                }
                                
                cpProgram = clCreateProgramWithSource( cxGPUContext, sources.size(), &(c_strs[0]), &(lengths[0]), &err );
                device::check(err);

                // building the OpenCL source code
                err = clBuildProgram( cpProgram, 0, NULL, NULL, NULL, NULL );
                if ( err!=CL_SUCCESS )
                {
                    size_t len;
                    std::string buffer(1048576, 0);

                    device::check( clGetProgramBuildInfo( cpProgram, cdDevice, CL_PROGRAM_BUILD_LOG, buffer.size(), reinterpret_cast<void*>(&buffer[0]), &len ) );

                    std::cerr << buffer << std::endl;
                    throw std::runtime_error("error: OpenCL: The compilation has failed.");                    
                }

                // extracting the kernel entrances
                BOOST_FOREACH( std::string kernel_name, kernelnames )
                {
                    cl_int err;                    
                    kernels[kernel_name] = clCreateKernel( cpProgram, kernel_name.c_str(), &err);
                    device::check(err);                    
                }                

                return;
            } // compile

            kernel & operator [] ( const std::string & filename ) {
                return kernels[filename];
            } // operator[]
        }; // class device

        
    } // namespace opencl
} // namespace carp



#endif /* __OPENCL_HPP__ */
// LuM end of file
