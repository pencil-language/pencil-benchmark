// UjoImro, 2013
// This file tests the OpenCL framework for the CARP project

#include <stdlib.h>
#include "opencl.hpp"

const int rows = 10;
const int cols = 15;

int main()
{
    carp::opencl::device device;
    device.compile( {"color.cl"}, {"color"} );

    auto cl_image = device.create_image<float>()
    
    device["color"](a, b);

    return EXIT_SUCCESS;    
}



// LuM end of file
