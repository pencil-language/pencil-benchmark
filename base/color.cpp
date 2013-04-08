// UjoImro, 2013
// This file tests the OpenCL framework for the CARP project

#include <stdlib.h>
#include "opencl.hpp"

const int rows = 10;
const int cols = 15;

int main()
{
    PRINT("ici01");    
    carp::opencl::device device;
    device.compile( {"color.cl"}, {"color"} );
    PRINT("ici02");
    
    carp::opencl::image<int> cl_image(device, 10, 25);
    PRINT("ici03");
    
    // device["color"]( cl_image.cl(), cl_image.rows(), cl_image.cols() )        
    //     .groupsize({16,16}, {cl_image.rows(), cl_image.cols()});
    // PRINT("ici04");

    cv::Mat_<int> image = cl_image.get();
    PRINT("ici05");

    print_image( image, "image");
        
    return EXIT_SUCCESS;
}



// LuM end of file
