// -*- c++ -*- 
// UjoImro, 2013
// CARP experimental opencl code


typedef struct {
    int rows;
    int cols;
    int step;
    int start;    
} Mat; // struct Mat

void
transposeFloatGang( Mat input,
                    __global float * input_data,
                    int localid,
                    int gangsize,
                    Mat output,
                    __global float * output_data )
{
    int worksize = input.cols * input.rows;
    int work = 0;
            
    for ( work=0; work < (worksize/gangsize) + 1; work++) {
        int index = work * gangsize + localid;
        if (index>=worksize) continue;
        int q = index / input.cols;
        int w = index % input.cols;
              
        output_data[ w * output.step + q + output.start ]
            = input_data[ q * input.step + w + input.start ];
    }
} // transposeFloatGang



__kernel void
transposeFloat( Mat input,
                __global float * input_data,
                Mat output,
                __global float * output_data )
{
    int localid = get_local_id(0);    
    int gangsize = get_local_size(0);
    
    transposeFloatGang( input, input_data, localid, gangsize, output, output_data );
    
    return;
    
} // transposeFloat


// LuM end of file
