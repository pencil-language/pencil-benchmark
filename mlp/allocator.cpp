// UjoImro, 2013
// Experimental code for the CARP Project
// Copyright (c) RealEyes, 2013

#include "memory.hpp"
#include "allocator.hpp"

#include <cassert>

clMat /*float*/
carp::CreateMatFloat( carp::memory::allocator & pool, int rows, int cols )
{
    clMat /*float*/result; 
    assert(rows>0);
    assert(cols>0);
    result.rows  = rows;
    result.cols  = cols;
    int step = cols;    
    result.step  = step;
    result.start = pool.allocate(rows * step * sizeof(float)) / sizeof(float);

    return result;
}

clMat /*uint8_t*/ 
carp::CreateMatChar /*uint8_t*/ ( carp::memory::allocator & pool, int rows, int cols )
{
    clMat /*uint8_t*/result; 
    assert(rows>0);
    assert(cols>0);
    result.rows  = rows;
    result.cols  = cols;
    int step = cols;    
    result.step  = step;
    result.start = pool.allocate( rows * step );

    return result;
}
