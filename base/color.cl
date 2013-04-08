// -*- c -*-

__kernel void color( __global int* image, int rows, int cols )
{
    // get index into global data array
    int row = get_global_id(0);
    int col = get_global_id(1);

    int id = row*cols + col;

    if ( id >= rows*cols)
    {
	return;
    }
    
    image[id] = id;
    
}
