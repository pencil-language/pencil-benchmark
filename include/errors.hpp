// UjoImro, 2013
// OpenCL code for the CARP Project

#ifndef __ERRORS__HPP__
#define __ERRORS__HPP__

#include <map>
#include <string>
#include <CL/cl.h>

namespace carp {
    namespace opencl {
#define MAKE_ERROR_PAIR(err) { err, #err }
        static const std::map<int, std::string> errors{ MAKE_ERROR_PAIR(CL_SUCCESS)
                                                      , MAKE_ERROR_PAIR(CL_DEVICE_NOT_FOUND)
                                                      , MAKE_ERROR_PAIR(CL_DEVICE_NOT_AVAILABLE)
                                                      , MAKE_ERROR_PAIR(CL_COMPILER_NOT_AVAILABLE)
                                                      , MAKE_ERROR_PAIR(CL_MEM_OBJECT_ALLOCATION_FAILURE)
                                                      , MAKE_ERROR_PAIR(CL_OUT_OF_RESOURCES)
                                                      , MAKE_ERROR_PAIR(CL_OUT_OF_HOST_MEMORY)
                                                      , MAKE_ERROR_PAIR(CL_PROFILING_INFO_NOT_AVAILABLE)
                                                      , MAKE_ERROR_PAIR(CL_MEM_COPY_OVERLAP)
                                                      , MAKE_ERROR_PAIR(CL_IMAGE_FORMAT_MISMATCH)
                                                      , MAKE_ERROR_PAIR(CL_IMAGE_FORMAT_NOT_SUPPORTED)
                                                      , MAKE_ERROR_PAIR(CL_BUILD_PROGRAM_FAILURE)
                                                      , MAKE_ERROR_PAIR(CL_MAP_FAILURE)
                                                      , MAKE_ERROR_PAIR(CL_MISALIGNED_SUB_BUFFER_OFFSET)
                                                      , MAKE_ERROR_PAIR(CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST)
#if __OPENCL_VERSION__ >= 120
                                                      , MAKE_ERROR_PAIR(CL_COMPILE_PROGRAM_FAILURE)
                                                      , MAKE_ERROR_PAIR(CL_LINKER_NOT_AVAILABLE)
                                                      , MAKE_ERROR_PAIR(CL_LINK_PROGRAM_FAILURE)
                                                      , MAKE_ERROR_PAIR(CL_DEVICE_PARTITION_FAILED)
                                                      , MAKE_ERROR_PAIR(CL_KERNEL_ARG_INFO_NOT_AVAILABLE)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_VALUE)
#endif /* __OPENCL_VERSION__ >= 120 */
                                                      , MAKE_ERROR_PAIR(CL_INVALID_DEVICE_TYPE)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_PLATFORM)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_DEVICE)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_CONTEXT)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_QUEUE_PROPERTIES)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_COMMAND_QUEUE)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_HOST_PTR)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_MEM_OBJECT)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_IMAGE_FORMAT_DESCRIPTOR)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_IMAGE_SIZE)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_SAMPLER)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_BINARY)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_BUILD_OPTIONS)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_PROGRAM)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_PROGRAM_EXECUTABLE)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_KERNEL_NAME)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_KERNEL_DEFINITION)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_KERNEL)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_ARG_INDEX)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_ARG_VALUE)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_ARG_SIZE)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_KERNEL_ARGS)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_WORK_DIMENSION)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_WORK_GROUP_SIZE)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_WORK_ITEM_SIZE)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_GLOBAL_OFFSET)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_EVENT_WAIT_LIST)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_EVENT)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_OPERATION)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_GL_OBJECT)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_BUFFER_SIZE)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_MIP_LEVEL)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_GLOBAL_WORK_SIZE)
#if __OPENCL_VERSION__ >= 120
                                                      , MAKE_ERROR_PAIR(CL_INVALID_PROPERTY)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_IMAGE_DESCRIPTOR)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_COMPILER_OPTIONS)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_LINKER_OPTIONS)
                                                      , MAKE_ERROR_PAIR(CL_INVALID_DEVICE_PARTITION_COUNT)
#endif /* __OPENCL_VERSION__ >= 120 */
                                                      };
#undef MAKE_ERROR_PAIR
    }
}
#endif /* __ERRORS__HPP__ */
