// UjoImro, 2013

#include <stdlib.h>
#include <opencv2/core/core.hpp>

#include "opencl.hpp"
#include "utility.hpp"

const int rows = 10;
const int cols = 15;

int main()
{

    carp::opencl::device device;
    device.compile( carp::string_vector("local.cl"), carp::string_vector("paint") );

    carp::opencl::image<float> cl_image(device, rows, cols);
    
    device["paint"]( cl_image.cl(), cl_image.ptr(), carp::opencl::buffer(1024) )
        .groupsize( carp::make_vector<ulong>(cl_image.cols(), cl_image.rows()), carp::make_vector<ulong>(0,0));

    cv::Mat_<float> image = cl_image.get();

    print_image( image, "image" );
    
    for (int q = 0; q < image.rows; q++ )
        for ( int w = 0; w < image.cols; w++ )
            assert( image(q,w) == q * image.cols + w );
            
    return EXIT_SUCCESS;
    
} // main
    



// LuM end of file
