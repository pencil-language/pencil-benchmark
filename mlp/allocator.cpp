// UjoImro, 2013
// Experimental code for the CARP Project
// Copyright (c) RealEyes, 2013

#include "memory.hpp"
#include "allocator.hpp"


void carp::freeMLP( carp::memory::allocator & pool, carp::mlp * classifier )
{
    carp::freeMatFloat( pool, &classifier->m_wIn );
    carp::freeMatFloat( pool, &classifier->m_wOut );
    carp::freeMatFloat( pool, &classifier->m_U );
} // freeMLP


clMat /*float*/
carp::CreateMatFloat( carp::memory::allocator & pool, int rows, int cols )
{
    clMat /*float*/result; 
    assert(rows>0);
    assert(cols>0);

    // result.data = NULL;
    // result.data = (float*)malloc( sizeof(float) * rows * cols );
    // assert(result.data);
    result.rows  = rows;
    result.cols  = cols;
    // int r = cols % 32;    
    // int step = cols - r + 32 + 1;
    int step = cols;    
    result.step  = step;
    result.start = pool.allocate(carp::memory::allocator::sizer(rows * step, float()));

    return result;
} // CreateMatFloat

clMat /*uint8_t*/ 
carp::CreateMatChar /*uint8_t*/ ( carp::memory::allocator & pool, int rows, int cols )
{
    clMat /*float*/result; 
    assert(rows>0);
    assert(cols>0);

    // result.data = NULL;
    // result.data = (float*)malloc( sizeof(float) * rows * cols );
    // assert(result.data);
    result.rows  = rows;
    result.cols  = cols;
    // int r = cols % (4 * 32);
    // int step = cols - r + 4 * 32 + 4;
    int step = cols;    
    result.step  = step;
    result.start = pool.allocate( carp::memory::allocator::sizer(rows * step, uint8_t()));

    return result;
}

clVector /* <clMat> */
carp::CreateVectorMat( carp::memory::allocator & pool, int nb_elements )
{
    clVector result;
    assert(nb_elements>0);

    result.size = nb_elements;
    result.step = 1;
    result.start = pool.allocate(carp::memory::allocator::sizer(nb_elements, clMat()));

    return result;        
} // CreateVectorMat
       
    
void
carp::freeMatFloat( carp::memory::allocator & pool, clMat /*float*/* mat )
{
    pool.release(carp::memory::allocator::sizer(mat->start, float()));
    
    mat->rows  = 0;
    mat->cols  = 0;
    mat->step  = 0;
    mat->start = 0;
    return;
} // freeMatFloat

void 
carp::freeMatChar /*uint8_t*/ ( carp::memory::allocator & pool, clMat /*uint8_t*/  * mat )
{
    pool.release(carp::memory::allocator::sizer(mat->start, uint8_t()));
    
    mat->rows  = 0;
    mat->cols  = 0;
    mat->step  = 0;
    mat->start = 0;
    return;
} // freeclMat /*uint8_t*/

void 
carp::freeVectorMat( carp::memory::allocator & pool, clVector * vec )
{
    pool.release(carp::memory::allocator::sizer(vec->start, clMat()));
    vec->size  = 0;        
    vec->step  = 0;        
    vec->start = 0;
} // freeVectorMat <clMat>

// LuM end of file
