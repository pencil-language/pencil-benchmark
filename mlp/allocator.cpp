// UjoImro, 2013
// Experimental code for the CARP Project
// Copyright (c) RealEyes, 2013

#include "memory.hpp"
#include "allocator.hpp"


void freeMLP( carp::memory::allocator & pool, mlp * classifier )
{
    freeMatFloat( pool, &classifier->m_wIn );
    freeMatFloat( pool, &classifier->m_wOut );
    freeMatFloat( pool, &classifier->m_U );
} // freeMLP


clMat /*float*/
CreateMatFloat( carp::memory::allocator & pool, int rows, int cols )
{
    clMat /*float*/result; 
    assert(rows>0);
    assert(cols>0);

    // result.data = NULL;
    // result.data = (float*)malloc( sizeof(float) * rows * cols );
    // assert(result.data);
    result.rows  = rows;
    result.cols  = cols;
    result.step  = cols;
    result.start = pool.allocate({rows * cols, float()});

    return result;
} // CreateMatFloat

clMat /*uint8_t*/ 
CreateMatChar /*uint8_t*/ ( carp::memory::allocator & pool, int rows, int cols )
{
    clMat /*float*/result; 
    assert(rows>0);
    assert(cols>0);

    // result.data = NULL;
    // result.data = (float*)malloc( sizeof(float) * rows * cols );
    // assert(result.data);
    result.rows  = rows;
    result.cols  = cols;
    result.step  = cols;
    result.start = pool.allocate({rows * cols, uint8_t()});

    return result;
}

clVector /* <clMat> */
CreateVectorMat( carp::memory::allocator & pool, int nb_elements )
{
    clVector result;
    assert(nb_elements>0);

    result.size = nb_elements;
    result.step = 1;
    result.start = pool.allocate({nb_elements, clMat()});

    return result;        
} // CreateVectorMat
       
    
void
freeMatFloat( carp::memory::allocator & pool, clMat /*float*/* mat )
{
    pool.release({mat->start, float()});
    
    mat->rows  = 0;
    mat->cols  = 0;
    mat->step  = 0;
    mat->start = 0;
    return;
} // freeMatFloat

void 
freeMatChar /*uint8_t*/ ( carp::memory::allocator & pool, clMat /*uint8_t*/  * mat )
{
    pool.release({mat->start, uint8_t()});
    
    mat->rows  = 0;
    mat->cols  = 0;
    mat->step  = 0;
    mat->start = 0;
    return;
} // freeclMat /*uint8_t*/

void 
freeVectorMat( carp::memory::allocator & pool, clVector * vec )
{
    pool.release({vec->start, clMat()});
    vec->size  = 0;        
    vec->step  = 0;        
    vec->start = 0;
} // freeVectorMat <clMat>

// LuM end of file
