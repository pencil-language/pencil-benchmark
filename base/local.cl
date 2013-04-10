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
paintloc( __local float * buff, int index )    
{
    buff[index] = index;
} // paintloc



__kernel void
paint( Mat output,
       __global float * output_data,
       __local float * buff
    )
{
    int col = get_local_id(0);
    if (col >= output.cols) return;
    
    int row = get_local_id(1);
    if (row >= output.rows) return;

    int groupid = get_group_id(0);
    if (groupid > 0) return;

    int output_index = row * output.step + col + output.start;
    int local_index = row * output.cols + col;

    paintloc( buff, local_index  );

    output_data[output_index] = buff[local_index];

    
    
    return;
    
} // paint


// LuM end of file
