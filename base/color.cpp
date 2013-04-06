// UjoImro, 2013
// This file tests the OpenCL framework for the CARP project

#include <stdlib.h>
#include "opencl.hpp"

int main()
{
    carp::opencl::device device;
    device.compile( {"color.cl"}, {"color"} );

    int a;
    int b;
    
    device["color"](a, b);

    return EXIT_SUCCESS;    
}



// LuM end of file
