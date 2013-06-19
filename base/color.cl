// -*- c -*-

typedef struct {
    int rows;
    int cols;
    int step;
    int start;    
} clMat; // struct Mat


__kernel void color( clMat image, __global int * data )
{
    // get index into global data array
    int row = get_global_id(0);
    int col = get_global_id(1);

    int id = row * image.step + col + image.start;

    if ( id >= image.rows * image.cols)
    {
	return;
    }
    
    data[id] = id;
    
} // color
