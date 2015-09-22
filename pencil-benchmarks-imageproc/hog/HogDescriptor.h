#ifndef HOGDESCRIPTOR_H
#define HOGDESCRIPTOR_H

#include <opencv2/core/core.hpp>

#define NOMINMAX
#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>

namespace nel {
    template<int numberOfCells, int numberOfBins, bool gauss, bool spinterp, bool _signed, bool static_>
    class HOGDescriptorCPP {
        static_assert(numberOfCells > 1 || !spinterp, "Cannot apply spatial interpolation with only one cell.");
    public:
        HOGDescriptorCPP();

        //Wrapper for a single float, to have an operator() to "index" it
        class StaticBlockSizeParameterCPP
        {
        public:
            StaticBlockSizeParameterCPP(float f) : value(f) {}
            float operator()(int, int) const { return value; }
        private:
            float value;
        };
        template<bool is_static>
        using BlockSizeParameter = typename std::conditional<is_static, StaticBlockSizeParameterCPP, cv::Mat_<float>>::type;

        cv::Mat_<float> compute( const cv::Mat_<uint8_t>           &img
                               , const cv::Mat_<float>             &locations
                               , const BlockSizeParameter<static_> &blocksizes
                               ) const;

        static int getNumberOfBins();

    private:
        static float get_orientation(int mdy, int mdx);

    private:
        std::vector<std::pair<float, float>> m_lookupTable;
    };

    template<size_t numberOfCells, size_t numberOfBins, bool gauss, bool spinterp, bool _signed, bool _static>
    class HOGDescriptorOCL {
        static_assert(numberOfCells > 1 || !spinterp, "Cannot apply spatial interpolation with only one cell.");
    public:
        HOGDescriptorOCL();

        //Wrapper for a single cl_float
        class StaticBlockSizeParameterOCL
        {
        public:
            StaticBlockSizeParameterOCL(cl_float f) : value(f) {}
            StaticBlockSizeParameterOCL convert_to_opencl_data(const cl::Context &, const cl::CommandQueue &) const { return *this; }
            cl_float operator()() const { return value; }
            bool compatiblewith(const cv::Mat_<float> &) const { return true; }
        private:
            cl_float value;
        };
        //Wrapper for cv::Mat_<float>, to create a cl::Buffer
        class DynamicBlockSizeParameterOCL
        {
        public:
            DynamicBlockSizeParameterOCL(const cv::Mat_<float> &blocksizes_) : blocksizes(blocksizes_) {}
            cl::Buffer convert_to_opencl_data(const cl::Context &context, const cl::CommandQueue &queue) const {
                size_t blocksizes_bytes = sizeof(cl_float2)*blocksizes.rows;
                cl::Buffer blocksizes_cl = cl::Buffer(context, CL_MEM_READ_ONLY , blocksizes_bytes);
                queue.enqueueWriteBuffer(blocksizes_cl, CL_FALSE, 0, blocksizes_bytes, blocksizes.data);
                return blocksizes_cl;
            }
            bool compatiblewith(const cv::Mat_<float> &locations) const {
                return (2 == blocksizes.cols) && (blocksizes.isContinuous()) && (locations.rows == blocksizes.rows);
            }
        private:
            const cv::Mat_<float> &blocksizes;
        };
        template<bool is_static>
        using BlockSizeParameter = typename std::conditional<is_static, StaticBlockSizeParameterOCL, DynamicBlockSizeParameterOCL>::type;

        cv::Mat_<float> compute( const cv::Mat_<uint8_t>           &img
                               , const cv::Mat_<float>             &locations
                               , const BlockSizeParameter<_static> &blocksizes
                               ) const;

        static size_t getNumberOfBins();

    private:
        enum HOGAlgorithmType { use_private, use_local, use_global };
        HOGAlgorithmType getAlgorithmType() const;

        cl::Device       m_device;
        cl::Context      m_context;
        cl::Program      m_program;
        cl::CommandQueue m_queue;

        mutable cl::Kernel m_calc_hog;
                size_t     m_calc_hog_preferred_multiple;
                size_t     m_calc_hog_group_size;

        bool m_has_local_memory;

        static const size_t ms_calchog_vector_size = 4;
    };
}

#include "HogDescriptor.hpp"

#endif
