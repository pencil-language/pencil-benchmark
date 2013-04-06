// -*- c -*-

__kernel void color( __global float* image )
{
    // get index into global data array
    int row = get_global_id(0);
    int col = get_global_id(1);
// int id;
//    row = id/cols;
//    col = id%cols;
    
    int rows = get_global_size(0);
    int cols = get_global_size(1);

    int id = row*cols + col;

    if ( id >= rows*cols)
    {   
        return; 
    }
    
    image[id] = id;
    
}
