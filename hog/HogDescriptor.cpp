/* ---------------------------------------------------------------------
* Copyright Realeyes OU 2012-2014
*
* All rights reserved. Realeyes OU PROPRIETARY/CONFIDENTIAL.
* No part of this computer programs(s) may be used, reproduced,
* stored in any retrieval system, or transmitted, in any form or
* by any means, electronic, mechanical, photocopying,
* recording, or otherwise without prior written permission of
* Realeyes OU. Use is subject to license terms.
* ---------------------------------------------------------------------
*/

#include "HogDescriptor.h"

#include <algorithm>
#include <cfloat>
#include <cassert>
#include <array>
#include <fstream>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

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
}

nel::HOGDescriptorCPP::HOGDescriptorCPP(int numberOfCells_, int numberOfBins_, bool gauss_, bool spinterp_, bool _signed_)
    : numberOfCells(numberOfCells_)
    , numberOfBins (numberOfBins_ )
    , gauss        (gauss_        )
    , spinterp     (spinterp_     )
    , _signed      (_signed_      )
{
    assert(numberOfCells > 1 || !spinterp);
    m_lookupTable.resize(512 * 512);
    for (int mdy = -255; mdy < 256; ++mdy)
        for (int mdx = -255; mdx < 256; ++mdx) {
            m_lookupTable[(mdy + 255) * 512 + mdx + 255].first = get_orientation(mdy,mdx);
            m_lookupTable[(mdy + 255) * 512 + mdx + 255].second = std::hypot(mdx, mdy);
        }
}

int nel::HOGDescriptorCPP::getNumberOfBins() const {
    return numberOfCells * numberOfCells * numberOfBins;
}

float nel::HOGDescriptorCPP::get_orientation(float mdy, float mdx) const {
    if (_signed) {
        return std::atan2(mdy, mdx) / static_cast<float>(M_PI) / 2.0f;
    } else {
        return std::atan2(mdy, mdx) / static_cast<float>(M_PI) + 0.5f;
    }
}

cv::Mat_<float> nel::HOGDescriptorCPP::compute( const cv::Mat_<uint8_t>  &image
                                              , const cv::Mat_<float> &locations
                                              , const cv::Mat_<float> &blocksizes
                                              ) const
{
    assert(2 == locations.cols);
    assert(2 == blocksizes.cols);
    cv::Mat_<float> descriptors(locations.rows, getNumberOfBins(), 0.0f);

    tbb::parallel_for(tbb::blocked_range<size_t>(0, locations.rows, 5), [&](const tbb::blocked_range<size_t> range) {
    for (size_t n = range.begin(); n != range.end(); ++n) {

        const float &blocksizeX = blocksizes(n,0);
        const float &blocksizeY = blocksizes(n,1);
        const float centerx = static_cast<float>(locations(n, 0));
        const float centery = static_cast<float>(locations(n, 1));

        const float cellsizeX = blocksizeX / numberOfCells;
        const float cellsizeY = blocksizeY / numberOfCells;
        const float halfblocksizeX = blocksizeX / 2.0f;
        const float halfblocksizeY = blocksizeY / 2.0f;

        const float minx = centerx - halfblocksizeX;
        const float miny = centery - halfblocksizeY;
        const float maxx = centerx + halfblocksizeX;
        const float maxy = centery + halfblocksizeY;
        
        const int minxi = std::max(fast_ceil(minx) , 1);
        const int minyi = std::max(fast_ceil(miny) , 1);
        const int maxxi = std::min(fast_floor(maxx), image.cols - 2);
        const int maxyi = std::min(fast_floor(maxy), image.rows - 2);

        cv::Mat_<float> hist(numberOfCells * numberOfCells, numberOfBins, 0.0f);

        float m1p2sigmaX2;
        float m1p2sigmaY2;
        if (gauss) {
            float sigmaX = halfblocksizeX;
            float sigmaY = halfblocksizeY;
            float sigmaX2 = sigmaX*sigmaX;
            float sigmaY2 = sigmaY*sigmaY;
            m1p2sigmaX2 = -1.0f/(2.0f*sigmaX2);
            m1p2sigmaY2 = -1.0f/(2.0f*sigmaY2);
            //float A = 1.0f / (2 * (float)M_PI*sigma2); // constant multiplier to all bin elements, but we normalize at the end, so we can skip this
        }

        // compute edges, magnitudes, orientations and the histogram
        for (int pointy = minyi; pointy <= maxyi; pointy++) {
            #pragma GCC diagnostic push
            #pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
            int cellyi;
            #pragma GCC diagnostic pop
            float yscale0;
            float yscale1;
            if (spinterp) {
                //Relative position of the pixel compared to the cell centers - y dimension
                float relative_pos_y = (pointy - miny) / cellsizeY - 0.5f;
                //Calculate the integral part and the fractional part of the relative position - y dimension
                cellyi = fast_floor(relative_pos_y);

                yscale1 = relative_pos_y - cellyi;
                yscale0 = 1.0f - yscale1;
            }

            float dy2;
            if (gauss) {
                float dy = pointy - centery;
                dy2 = dy*dy;
            }

            for (int pointx = minxi; pointx <= maxxi; pointx++) {
                int mdxi = image(pointy, pointx + 1) - image(pointy, pointx - 1);
                int mdyi = image(pointy + 1, pointx) - image(pointy - 1, pointx);

                float magnitude = m_lookupTable[(mdyi + 255) * 512 + mdxi + 255].second;
                float orientation = m_lookupTable[(mdyi + 255) * 512 + mdxi + 255].first;

                if (gauss) {
                    float dx = pointx - centerx;
                    float dx2 = dx*dx;
                    #pragma GCC diagnostic push
                    #pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
                    float B = std::exp(dx2 * m1p2sigmaX2 + dy2 * m1p2sigmaY2);
                    #pragma GCC diagnostic pop
                    magnitude *= B;
                }

                // linear/trilinear interpolation of magnitudes
                float relative_orientation = orientation * numberOfBins - 0.5f;
                int bin1 = fast_ceil(relative_orientation);
                int bin0 = bin1 - 1;
                float magscale0 = magnitude * (bin1 - relative_orientation);
                float magscale1 = magnitude * (relative_orientation - bin0);
                int binidx0 = (bin0 + numberOfBins) % numberOfBins;
                int binidx1 = (bin1 + numberOfBins) % numberOfBins;

                if (spinterp) {
                    //Relative position of the pixel compared to the cell centers - x dimension
                    float relative_pos_x = (pointx - minx) / cellsizeX - 0.5f;
                    //Calculate the integral part and the fractional part of the relative position - x dimension
                    int cellxi = fast_floor(relative_pos_x);

                    float xscale1 = relative_pos_x - cellxi;
                    float xscale0 = 1.0f - xscale1;

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
                        int cellxi = fast_floor((pointx - minx) / cellsizeX);
                        int cellyi = fast_floor((pointy - miny) / cellsizeY);
                        assert(cellxi < numberOfCells);
                        assert(cellyi < numberOfCells);
                        assert(cellxi >= 0);
                        assert(cellyi >= 0);
                        hist(cellyi * numberOfCells + cellxi, binidx0) += magscale0;
                        hist(cellyi * numberOfCells + cellxi, binidx1) += magscale1;
                    }
                }
            }
        }
        cv::Mat_<float> destRow = descriptors.row(n);
        hist.reshape(0,1).copyTo(destRow);
    }});
    return descriptors;
}

nel::HOGDescriptorOCL::HOGDescriptorOCL(int numberOfCells_, int numberOfBins_, bool gauss, bool spinterp, bool _signed)
    : numberOfCells(numberOfCells_)
    , numberOfBins(numberOfBins_)
{
    assert(numberOfCells > 1 || !spinterp);

    //Get the used device
    std::vector<cl::Device> devices;
    cl::Platform::getDefault().getDevices(CL_DEVICE_TYPE_GPU, &devices);
    device = devices.at(0);

    //Create context
    context = cl::Context(device);
       
    //Load source
    std::ifstream source_file{"../hog/hog.opencl.cl"};
    std::string source{ std::istreambuf_iterator<char>{source_file}, std::istreambuf_iterator<char>{} };

    //Create program
    program = cl::Program(context, source);

    //Build program
    std::stringstream build_opts;
    //-cl-unsafe-math-optimizations might work, but increases HOG 2x2 error over 1e-6 by a bit.
    //-cl-fast-relaxed-math gives totally wrong results
    build_opts << "-cl-single-precision-constant -cl-denorms-are-zero -cl-no-signed-zeros -cl-finite-math-only";
    build_opts << " -D NUMBER_OF_CELLS=" << std::to_string(numberOfCells);
    build_opts << " -D NUMBER_OF_BINS=" << std::to_string(numberOfBins);
    build_opts << " -D GAUSSIAN_WEIGHTS=" << (gauss ? "1" : "0");
    build_opts << " -D SPARTIAL_WEIGHTS=" << (spinterp ? "1" : "0");
    build_opts << " -D SIGNED_HOG=" << (_signed ? "1" : "0");
    try {
        program.build(build_opts.str().c_str());
    } catch (const cl::Error&) {
        auto buildlog = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device);
        throw;
    }

    //Create queue (apparently it's slow)
    queue = cl::CommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE);

    //Query kernels and work group infos
    calc_hog  = cl::Kernel(program, "calc_histogram");
    calc_hog_preferred_multiple = calc_hog.getWorkGroupInfo<CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE>(device);
    calc_hog_group_size         = calc_hog.getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(device);

#ifndef CL_VERSION_1_2
    fill_zeros = cl::Kernel(program, "fill_zeros");
    fill_zeros_preferred_multiple = fill_zeros.getWorkGroupInfo<CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE>(device);
    fill_zeros_group_size         = fill_zeros.getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(device);
#endif
    is_unified_host_memory = (CL_TRUE == device.getInfo<CL_DEVICE_HOST_UNIFIED_MEMORY>());
}

int nel::HOGDescriptorOCL::getNumberOfBins() const {
    return numberOfCells * numberOfCells * numberOfBins;
}

inline size_t round_to_multiple(size_t num, size_t factor) {
    return num + factor - 1 - (num - 1) % factor;
}

#include <iostream>

cv::Mat_<float> nel::HOGDescriptorOCL::compute( const cv::Mat_<uint8_t> &image
                                              , const cv::Mat_<float>   &locations
                                              , const cv::Mat_<float>   &blocksizes
                                              , const size_t             max_blocksize_x
                                              , const size_t             max_blocksize_y
                                              , std::chrono::duration<double> &elapsed_time_gpu_nocopy
                                              ) const
{
    assert(2 == locations.cols);
    assert(locations.isContinuous());
    assert(2 == blocksizes.cols);
    assert(blocksizes.isContinuous());
    int num_localtions = locations.rows;
    cv::Mat_<float> descriptors(num_localtions, getNumberOfBins());
    
    //Events to track execution time
    cl::Event start, end;

    //OPENCL START
    //Prepare OpenCL buffers
    size_t      image_bytes = image.elemSize()*image.rows*image.step1();
    size_t  locations_bytes = sizeof(cl_float2)*num_localtions;
    size_t blocksizes_bytes = sizeof(cl_float2)*blocksizes.rows;
    size_t descriptor_bytes = sizeof(cl_float)*getNumberOfBins()*num_localtions;
    cl::Buffer      image_cl;
    cl::Buffer  locations_cl;
    cl::Buffer blocksizes_cl;
    cl::Buffer descriptor_cl;
    if (is_unified_host_memory) {
             image_cl = cl::Buffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,      image_bytes,       image.data);
         locations_cl = cl::Buffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,  locations_bytes,   locations.data);
        blocksizes_cl = cl::Buffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, blocksizes_bytes,  blocksizes.data);
        descriptor_cl = cl::Buffer(context, CL_MEM_READ_WRITE| CL_MEM_USE_HOST_PTR, descriptor_bytes, descriptors.data);
    } else {
             image_cl = cl::Buffer(context, CL_MEM_READ_ONLY ,      image_bytes);
         locations_cl = cl::Buffer(context, CL_MEM_READ_ONLY ,  locations_bytes);
        blocksizes_cl = cl::Buffer(context, CL_MEM_READ_ONLY , blocksizes_bytes);
        descriptor_cl = cl::Buffer(context, CL_MEM_READ_WRITE, descriptor_bytes);
        queue.enqueueWriteBuffer(     image_cl, CL_FALSE, 0,      image_bytes,      image.data);
        queue.enqueueWriteBuffer( locations_cl, CL_FALSE, 0,  locations_bytes,  locations.data);
        queue.enqueueWriteBuffer(blocksizes_cl, CL_FALSE, 0, blocksizes_bytes, blocksizes.data);
    }
#ifdef CL_VERSION_1_2
    queue.enqueueFillBuffer(descriptor_cl, 0.0f, 0, descriptor_bytes, nullptr, &start);
#else
    {
        cl::NDRange  local_work_size(fill_zeros_preferred_multiple);
        cl::NDRange global_work_size( round_to_multiple(getNumberOfBins()*num_localtions, fill_zeros_preferred_multiple) );
        {
            //ADD MULTITHREAD LOCK HERE (if needed)
//            std::cout << getNumberOfBins()*num_localtions << '/' << fill_zeros_preferred_multiple << '/' << fill_zeros_group_size << std::endl;
            fill_zeros.setArg(0, descriptor_cl());
            fill_zeros.setArg(1, (cl_int)(getNumberOfBins()*num_localtions));
            queue.enqueueNDRangeKernel(fill_zeros, cl::NullRange, global_work_size, local_work_size, nullptr, &start );
        }
    }
#endif

    {
        //Figure out the work group sizes
        const size_t work_size_x = calc_hog_preferred_multiple;
        const size_t work_size_y = std::max(calc_hog_group_size / num_localtions / work_size_x, (size_t)1u);    //Integer division!
        const size_t work_size_z = std::max(calc_hog_group_size / work_size_y / work_size_x   , (size_t)1u);    //Integer division!
        cl::NDRange local_work_size(work_size_x, work_size_y, work_size_z);

        size_t max_x = round_to_multiple(max_blocksize_x, work_size_x);
        size_t max_y = round_to_multiple(max_blocksize_y, work_size_y);
        size_t max_z = round_to_multiple(num_localtions, work_size_z);
        cl::NDRange global_work_size(max_x, max_y, max_z);

        //Execute the kernel
        {
            //ADD MULTITHREAD LOCK HERE (if needed)
            calc_hog.setArg(0, (cl_int)image.rows);
            calc_hog.setArg(1, (cl_int)image.cols);
            calc_hog.setArg(2, (cl_int)image.step1());
            calc_hog.setArg(3, image_cl());
            calc_hog.setArg(4, (cl_int)num_localtions);
            calc_hog.setArg(5, locations_cl());
            calc_hog.setArg(6, blocksizes_cl());
            calc_hog.setArg(7, descriptor_cl());
            calc_hog.setArg(8, (size_t)(sizeof(cl_float) * getNumberOfBins() * work_size_z), nullptr);
            queue.enqueueNDRangeKernel(calc_hog, cl::NullRange, global_work_size, local_work_size, nullptr, &end );
        }
    }
    
    //Read result
    if (is_unified_host_memory) {
        auto mappedPtr = queue.enqueueMapBuffer(descriptor_cl, CL_TRUE, CL_MAP_READ, 0, descriptor_bytes);
        assert(mappedPtr == descriptors.data);
        queue.enqueueUnmapMemObject(descriptor_cl, mappedPtr);    //Buffer was created with CL_MEM_USE_HOST_PTR -> can unmap immediately
    } else {
        queue.enqueueReadBuffer(descriptor_cl, CL_TRUE, 0, descriptor_bytes, descriptors.data);
    }
    
    //Figure out execution time
    auto time_nanoseconds = end.getProfilingInfo<CL_PROFILING_COMMAND_END>() - start.getProfilingInfo<CL_PROFILING_COMMAND_START>();
    elapsed_time_gpu_nocopy = std::chrono::nanoseconds(time_nanoseconds);

    //Normalize result
//    for (int n = 0; n < num_localtions; ++n)
//        normalize(descriptors.row(n), descriptors.row(n), l2hys_thrs);
    return descriptors;
}
