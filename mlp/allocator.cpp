// UjoImro, 2013
// Experimental code for the CARP Project
// Copyright (c) RealEyes, 2013

#include "memory.hpp"
#include "allocator.hpp"

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

void freeMLP( void * self, carp::memory & pool, mlp * classifier )
{
    freeMatFloat( self, pool, &classifier->m_wIn );
    freeMatFloat( self, pool, &classifier->m_wOut );
    freeMatFloat( self, pool, &classifier->m_U );
} // freeMLP


cMat /*float*/
CreateMatFloat( void * self, carp::memory & pool, int rows, int cols )
{
    cMat /*float*/result; 
    assert(rows>0);
    assert(cols>0);

    // result.data = NULL;
    // result.data = (float*)malloc( sizeof(float) * rows * cols );
    // assert(result.data);
    result.rows  = rows;
    result.cols  = cols;
    result.step  = cols;
    result.start = pool.allocate<float>(rows * cols);

    return result;
} // CreateMatFloat

cMat /*uint8_t*/ 
CreateMatChar /*uint8_t*/ ( void * self, carp::memory & pool, int rows, int cols )
{
    cMat /*float*/result; 
    assert(rows>0);
    assert(cols>0);

    // result.data = NULL;
    // result.data = (float*)malloc( sizeof(float) * rows * cols );
    // assert(result.data);
    result.rows  = rows;
    result.cols  = cols;
    result.step  = cols;
    result.start = pool.allocate<uint8_t>(rows * cols);

    return result;
}

void
freeMatFloat( void * self, carp::memory & pool, cMat /*float*/* mat )
{
    pool.release<float>(mat->start);
    
    mat->rows  = 0;
    mat->cols  = 0;
    mat->step  = 0;
    mat->start = 0;
    return;
} // freeMatFloat

void 
freeMatChar /*uint8_t*/ ( void * self, carp::memory & pool, cMat /*uint8_t*/  * mat )
{
    pool.release<uint8_t>(mat->start);
    
    mat->rows  = 0;
    mat->cols  = 0;
    mat->step  = 0;
    mat->start = 0;
    return;
} // freecMat /*uint8_t*/










#ifdef __cplusplus
} // extern C
#endif // __cplusplus


// LuM end of file
