/* ---------------------------------------------------------------------
* Copyright Realeyes OU 2012-2015
*
* All rights reserved. Realeyes OU PROPRIETARY/CONFIDENTIAL.
* No part of this computer programs(s) may be used, reproduced,
* stored in any retrieval system, or transmitted, in any form or
* by any means, electronic, mechanical, photocopying,
* recording, or otherwise without prior written permission of
* Realeyes OU. Use is subject to license terms.
* ---------------------------------------------------------------------
*/

#ifndef HOGDESCRIPTOR_H
#error This is a private implementation file for HogDescriptor.h, do not include directly
#endif

#ifdef WITH_TBB
    #include <tbb/parallel_for.h>
    #include <tbb/blocked_range.h>
#endif

namespace {
    template<typename T>
    inline int fast_floor(T f) {
        static_assert(std::is_floating_point<T>::value, "fast_floor: Parameter must be floating point.");
        return static_cast<int>(f)-(f < T(0));
    }

    template<typename T>
    inline int fast_ceil(T f) {
        static_assert(std::is_floating_point<T>::value, "fast_ceil: Parameter must be floating point.");
        return static_cast<int>(f)+(f > T(0));
    }

    inline size_t round_to_multiple(size_t num, size_t factor) {
        return num + factor - 1 - (num - 1) % factor;
    }
}

template<int numberOfCells, int numberOfBins, bool gauss, bool spinterp, bool _signed, bool static_>
nel::HOGDescriptorCPP<numberOfCells, numberOfBins, gauss, spinterp, _signed, static_>::HOGDescriptorCPP()
{
    m_lookupTable.resize(512 * 512);
    for (int mdy = -255; mdy < 256; ++mdy)
        for (int mdx = -255; mdx < 256; ++mdx) {
            m_lookupTable[(mdy + 255) * 512 + mdx + 255].first  = get_orientation(mdy,mdx);
            m_lookupTable[(mdy + 255) * 512 + mdx + 255].second = static_cast<float>(std::hypot(mdx, mdy));
        }
}

template<int numberOfCells, int numberOfBins, bool gauss, bool spinterp, bool _signed, bool static_>
int nel::HOGDescriptorCPP<numberOfCells, numberOfBins, gauss, spinterp, _signed, static_>::getNumberOfBins() {
    return numberOfCells * numberOfCells * numberOfBins;
}

template<int numberOfCells, int numberOfBins, bool gauss, bool spinterp, bool _signed, bool static_>
float nel::HOGDescriptorCPP<numberOfCells, numberOfBins, gauss, spinterp, _signed, static_>::get_orientation(int mdy, int mdx) {
    if (_signed) {
        return static_cast<float>(std::atan2(mdy, mdx) * M_1_PI * 0.5);
    } else {
        return static_cast<float>(std::atan2(mdy, mdx) * M_1_PI + 0.5);
    }
}

template<int numberOfCells, int numberOfBins, bool gauss, bool spinterp, bool _signed, bool static_>
cv::Mat_<float> nel::HOGDescriptorCPP<numberOfCells, numberOfBins, gauss, spinterp, _signed, static_>
    ::compute( const cv::Mat_<uint8_t>           &image
             , const cv::Mat_<float>             &locations
             , const BlockSizeParameter<static_> &blocksizes
             ) const
{
    assert(2 == locations.cols);
    cv::Mat_<float> descriptors(locations.rows, getNumberOfBins(), static_cast<float>(0.0));

#ifdef WITH_TBB
    tbb::parallel_for(tbb::blocked_range<size_t>(0, locations.rows, 5), [&](const tbb::blocked_range<size_t> range) {
    for (size_t n = range.begin(); n != range.end(); ++n) {
#else
    for (size_t n = 0; n < (size_t)locations.rows; ++n) {
#endif

        const float &blocksizeX = blocksizes(n, 0);
        const float &blocksizeY = blocksizes(n, 1);
        const float &centerx    = locations (n, 0);
        const float &centery    = locations (n, 1);

        const float inv_cellsizeX = numberOfCells / blocksizeX;
        const float inv_cellsizeY = numberOfCells / blocksizeY;
        const float halfblocksizeX = blocksizeX * static_cast<float>(0.5);
        const float halfblocksizeY = blocksizeY * static_cast<float>(0.5);

        const float minx = centerx - halfblocksizeX;
        const float miny = centery - halfblocksizeY;
        const float maxx = centerx + halfblocksizeX;
        const float maxy = centery + halfblocksizeY;
        
        const int minxi = std::max(fast_ceil(minx) , 1);
        const int minyi = std::max(fast_ceil(miny) , 1);
        const int maxxi = std::min(fast_floor(maxx), image.cols - 2);
        const int maxyi = std::min(fast_floor(maxy), image.rows - 2);

        cv::Matx<float, numberOfCells * numberOfCells, numberOfBins> hist(0.0f);

        float m1p2sigmaX2;
        float m1p2sigmaY2;
        if (gauss) {
            float sigmaX = halfblocksizeX;
            float sigmaY = halfblocksizeY;
            float sigmaX2 = sigmaX*sigmaX;
            float sigmaY2 = sigmaY*sigmaY;
            m1p2sigmaX2 = static_cast<float>(-0.5) / sigmaX2;
            m1p2sigmaY2 = static_cast<float>(-0.5) / sigmaY2;
        }

        // compute edges, magnitudes, orientations and the histogram
        for (int pointy = minyi; pointy <= maxyi; pointy++) {
            #pragma GCC diagnostic push
            #pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
            int cellyi;
            float yscale0;
            float yscale1;
            if (spinterp) {
                //Relative position of the pixel compared to the cell centers - y dimension
                float relative_pos_y = (pointy - miny) * inv_cellsizeY - static_cast<float>(0.5);
                //Calculate the integral part and the fractional part of the relative position - y dimension
                cellyi = fast_floor(relative_pos_y);

                yscale1 = relative_pos_y - cellyi;
                yscale0 = static_cast<float>(1.0) - yscale1;
            }

            float dy2;
            if (gauss) {
                float dy = pointy - centery;
                dy2 = dy*dy;
            }

            for (int pointx = minxi; pointx <= maxxi; pointx++) {
                int mdxi = image(pointy, pointx + 1) - image(pointy, pointx - 1);
                int mdyi = image(pointy + 1, pointx) - image(pointy - 1, pointx);

                float magnitude   = m_lookupTable[(mdyi + 255) * 512 + mdxi + 255].second;
                float orientation = m_lookupTable[(mdyi + 255) * 512 + mdxi + 255].first;

                if (gauss) {
                    float dx = pointx - centerx;
                    float dx2 = dx*dx;
                    float B = std::exp(dx2 * m1p2sigmaX2 + dy2 * m1p2sigmaY2);
                    magnitude *= B;
                }

                // linear/trilinear interpolation of magnitudes
                float relative_orientation = orientation * numberOfBins - static_cast<float>(0.5);
                int bin1 = fast_ceil(relative_orientation);
                int bin0 = bin1 - 1;
                float magscale0 = magnitude * (bin1 - relative_orientation);
                float magscale1 = magnitude * (relative_orientation - bin0);
                int binidx0 = (bin0 + numberOfBins) % numberOfBins;
                int binidx1 = (bin1 + numberOfBins) % numberOfBins;

                if (spinterp) {
                    //Relative position of the pixel compared to the cell centers - x dimension
                    float relative_pos_x = (pointx - minx) * inv_cellsizeX - static_cast<float>(0.5);
                    //Calculate the integral part and the fractional part of the relative position - x dimension
                    int cellxi = fast_floor(relative_pos_x);

                    float xscale1 = relative_pos_x - cellxi;
                    float xscale0 = static_cast<float>(1.0) - xscale1;

                    if (cellyi >= 0                && cellxi >= 0) {
                        hist((cellyi + 0) * numberOfCells + cellxi + 0, binidx0) += yscale0 * xscale0 * magscale0;
                        hist((cellyi + 0) * numberOfCells + cellxi + 0, binidx1) += yscale0 * xscale0 * magscale1;
                    }
                    if (cellyi >= 0                && cellxi < numberOfCells - 1) {
                        hist((cellyi + 0) * numberOfCells + cellxi + 1, binidx0) += yscale0 * xscale1 * magscale0;
                        hist((cellyi + 0) * numberOfCells + cellxi + 1, binidx1) += yscale0 * xscale1 * magscale1;
                    }
                    if (cellyi < numberOfCells - 1 && cellxi >= 0) {
                        hist((cellyi + 1) * numberOfCells + cellxi + 0, binidx0) += yscale1 * xscale0 * magscale0;
                        hist((cellyi + 1) * numberOfCells + cellxi + 0, binidx1) += yscale1 * xscale0 * magscale1;
                    }
                    if (cellyi < numberOfCells - 1 && cellxi < numberOfCells - 1) {
                        hist((cellyi + 1) * numberOfCells + cellxi + 1, binidx0) += yscale1 * xscale1 * magscale0;
                        hist((cellyi + 1) * numberOfCells + cellxi + 1, binidx1) += yscale1 * xscale1 * magscale1;
                    }
                } else {
                    if (numberOfCells == 1) {
                        hist(0, binidx0) += magscale0;
                        hist(0, binidx1) += magscale1;
                    } else {
                        int cellxi = fast_floor((pointx - minx) * inv_cellsizeX);
                        int cellyi = fast_floor((pointy - miny) * inv_cellsizeY);
                        assert(cellxi < numberOfCells);
                        assert(cellyi < numberOfCells);
                        assert(cellxi >= 0);
                        assert(cellyi >= 0);
                        hist(cellyi * numberOfCells + cellxi, binidx0) += magscale0;
                        hist(cellyi * numberOfCells + cellxi, binidx1) += magscale1;
                    }
                }
            }
            #pragma GCC diagnostic pop
        }
        cv::Mat_<float> destRow = descriptors.row(n);
        cv::Mat(hist).reshape(0,1).copyTo(destRow);
    }
#ifdef WITH_TBB
    });
#endif
    return descriptors;
}

template<size_t numberOfCells, size_t numberOfBins, bool gauss, bool spinterp, bool _signed, bool _static>
typename nel::HOGDescriptorOCL<numberOfCells, numberOfBins, gauss, spinterp, _signed, _static>::HOGAlgorithmType
nel::HOGDescriptorOCL<numberOfCells, numberOfBins, gauss, spinterp, _signed, _static>::getAlgorithmType() const {
    if (getNumberOfBins() <= 36)
        return use_private;
    if (m_has_local_memory)
        return use_local;
    else
        return use_global;
}

template<size_t numberOfCells, size_t numberOfBins, bool gauss, bool spinterp, bool _signed, bool _static>
nel::HOGDescriptorOCL<numberOfCells, numberOfBins, gauss, spinterp, _signed, _static>::HOGDescriptorOCL()
{
    assert(numberOfCells > 1 || !spinterp);

    //Get the used device
    std::vector<cl::Device> devices;
    cl::Platform::getDefault().getDevices(CL_DEVICE_TYPE_GPU, &devices);
    m_device = devices.at(0);

    //Create context
    m_context = cl::Context(m_device);

    m_has_local_memory = (CL_LOCAL == m_device.getInfo<CL_DEVICE_LOCAL_MEM_TYPE>());

    //Load source
    std::ifstream source_file{ "HogDescriptor.cl" };
    std::string source{ std::istreambuf_iterator<char>{source_file}, std::istreambuf_iterator<char>{} };

    //Create program
    cl::Program program(m_context, source);
    std::string functionname = "hog" + std::to_string(numberOfCells) + 'x' + std::to_string(numberOfCells);

    //Build options
    std::stringstream build_opts;
    build_opts << " -cl-single-precision-constant";
    build_opts << " -cl-denorms-are-zero";
    build_opts << " -cl-mad-enable";
    build_opts << " -cl-no-signed-zeros";
    build_opts << " -cl-finite-math-only";
    build_opts << " -Werror";
    build_opts << " -D NUMBER_OF_CELLS="  << numberOfCells;
    build_opts << " -D NUMBER_OF_BINS="   << numberOfBins;
    build_opts << " -D GAUSSIAN_WEIGHTS=" << (gauss    ? "1" : "0");
    build_opts << " -D SPARTIAL_WEIGHTS=" << (spinterp ? "1" : "0");
    build_opts << " -D SIGNED_HOG="       << (_signed  ? "1" : "0");
    build_opts << " -D STATIC_BLOCK_SIZE="<< (_static  ? "1" : "0");
    build_opts << " -D FUNCTIONNAME="     << functionname;
    build_opts << " -D VECTOR_SIZE="      << ms_calchog_vector_size;
    switch (getAlgorithmType()) {
        case use_private: build_opts << " -D USE_PRIVATE=1"; break;
        case use_local  : build_opts << " -D USE_LOCAL=1";   break;
        case use_global: default: break;
    }
#ifndef NDEBUG
    build_opts << " -D DEBUG";
#endif

    //Build program
    try {
        program.build(build_opts.str().c_str());
    } catch (const cl::Error&) {
        auto buildlog = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(m_device);
        std::cerr << "OpenCL compilation error: " << buildlog << std::endl;
        throw;
    }

    //Create queue
    m_queue = cl::CommandQueue(m_context, m_device, CL_QUEUE_PROFILING_ENABLE);

    //Query kernels and work group infos
    m_calc_hog = cl::Kernel(program, functionname.c_str());
    m_calc_hog_preferred_multiple = m_calc_hog.getWorkGroupInfo<CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE>(m_device);
    m_calc_hog_group_size         = m_calc_hog.getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(m_device);
}

template<size_t numberOfCells, size_t numberOfBins, bool gauss, bool spinterp, bool _signed, bool _static>
size_t nel::HOGDescriptorOCL<numberOfCells, numberOfBins, gauss, spinterp, _signed, _static>::getNumberOfBins() {
    return numberOfCells * numberOfCells * numberOfBins;
}

template<size_t numberOfCells, size_t numberOfBins, bool gauss, bool spinterp, bool _signed, bool _static>
cv::Mat_<float> nel::HOGDescriptorOCL<numberOfCells, numberOfBins, gauss, spinterp, _signed, _static>
::compute( const cv::Mat_<uint8_t>           &image
         , const cv::Mat_<float>             &locations
         , const BlockSizeParameter<_static> &blocksizes
         ) const
{
    assert(2 == locations.cols);
    assert(locations.isContinuous());
    assert(blocksizes.compatiblewith(locations));

    int num_locations = locations.rows;
    cv::Mat_<float> descriptors(num_locations, getNumberOfBins());
    assert(descriptors.isContinuous());
    
    //Events to track execution time
    cl::Event event;

    //OPENCL START
    assert(image.elemSize() == sizeof(cl_uchar));
    static_assert(sizeof(cl_float) == sizeof(float), "Error: float and cl_float must be the same type");

    size_t      image_bytes = image.elemSize()*image.rows*image.step1();
    size_t  locations_bytes = sizeof(cl_float2)*num_locations;
    size_t descriptor_bytes = sizeof(cl_float)*getNumberOfBins()*num_locations;

    //Allocate OpenCL buffers
    cl::Buffer      image_cl = cl::Buffer(m_context, CL_MEM_READ_ONLY ,      image_bytes);
    cl::Buffer  locations_cl = cl::Buffer(m_context, CL_MEM_READ_ONLY ,  locations_bytes);
    cl::Buffer descriptor_cl = cl::Buffer(m_context, CL_MEM_READ_WRITE, descriptor_bytes);

    //Write input buffers to device
    m_queue.enqueueWriteBuffer(     image_cl, CL_FALSE, 0,      image_bytes,      image.data);
    m_queue.enqueueWriteBuffer( locations_cl, CL_FALSE, 0,  locations_bytes,  locations.data);
    auto blocksizes_cl = blocksizes.convert_to_opencl_data(m_context, m_queue);   //If dynamic: creates cl::Buffer like the rest. If static: creates a functor that returns the static value as cl_float

#ifndef LWS_X
#define LWS_X 4
#endif

#ifndef LWS_Y
#define LWS_Y 8
#endif

#ifndef LWS_Z
#define LWS_Z std::min( std::max<size_t>(num_locations / m_device.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>(), 1), m_calc_hog_group_size / (lws_x * lws_y))
#endif
    {
        //Figure out the work group sizes
        const size_t lws_x = LWS_X;
        const size_t lws_y = LWS_Y;
        const size_t lws_z = LWS_Z;
        cl::NDRange local_work_size(lws_x, lws_y, lws_z);

        const size_t gws_x = lws_x;
        const size_t gws_y = lws_y;
        const size_t gws_z = round_to_multiple(num_locations, lws_z);
        cl::NDRange global_work_size(gws_x, gws_y, gws_z);

        cl_int2 image_size = { image.cols, image.rows };

        //Execute the kernel
        {
            //ADD MULTITHREAD LOCK HERE (if needed)
            m_calc_hog.setArg(0, (cl_int2)image_size     );
            m_calc_hog.setArg(1, (cl_uint)image.step1()  );
            m_calc_hog.setArg(2, (cl_mem )image_cl()     );
            m_calc_hog.setArg(3, (cl_uint)num_locations  );
            m_calc_hog.setArg(4, (cl_mem )locations_cl() );
            m_calc_hog.setArg(5,          blocksizes_cl()); //Static: cl_float, Dynamic: cl_mem
            m_calc_hog.setArg(6, (cl_mem )descriptor_cl());
            switch (getAlgorithmType()) {
            case use_private:
                m_calc_hog.setArg(7, (size_t)(sizeof(cl_float) * lws_x * lws_y * lws_z * ms_calchog_vector_size), nullptr);
                break;
            case use_local:
                m_calc_hog.setArg(7, (size_t)(sizeof(cl_float) * getNumberOfBins() * lws_z), nullptr);
                break;
            case use_global:
                //nothing to do
                break;
            default:
                assert(false);
            }
        }
        m_queue.enqueueNDRangeKernel(m_calc_hog, cl::NullRange, global_work_size, local_work_size, nullptr, &event);
    }
    
    //Read result buffer from device
    m_queue.enqueueReadBuffer(descriptor_cl, CL_TRUE, 0, descriptor_bytes, descriptors.data);

    //Calculate kernel-only execution time
    double kernel_ns = event.getProfilingInfo<CL_PROFILING_COMMAND_END>() - event.getProfilingInfo<CL_PROFILING_COMMAND_START>();
    double kernel_ms = kernel_ns * 1e-6;
    std::cout << "calc_hog execution time: " << std::fixed << std::setprecision(6) << std::setw(8) << kernel_ms << " ms\n";

    return descriptors;
}
