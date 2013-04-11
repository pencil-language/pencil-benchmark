// UjoImro, 2013
// Experimental Research Code for the CARP Project

#include <iostream>
#include <boost/preprocessor.hpp>

#include "memory.hpp"

int main()
{
    carp::memory<float> pool(1000);

    PRINT(pool.allocate(10));
//    PRINT(pool.netallocated());
//    PRINT(pool.grossallocated());
    
    PRINT(pool.allocate(100));

    pool.release(128);

    PRINT(pool.allocate(100));

//    PRINT(pool.netallocated());
//    PRINT(pool.grossallocated());

    PRINT(pool.allocate(500));
//    PRINT(pool.netallocated());
//    PRINT(pool.grossallocated());

    try {
        
        PRINT(pool.allocate(500));
    }
    catch ( carp::exception & exception )
    {
        assert(exception.error == carp::INSUFFICIENT_MEMORY );
        PRINT("insufficient memory exception caught");
    }

    PRINT(pool.allocate(100));
//    PRINT(pool.netallocated());
//    PRINT(pool.grossallocated());

    PRINT(pool.allocate(50));
//    PRINT(pool.netallocated());
//    PRINT(pool.grossallocated());

    PRINT(pool.allocate(100));
    PRINT(pool.netallocated());
    PRINT(pool.grossallocated());
    
    
    PRINT( (1<<7) );
    
    return EXIT_SUCCESS;
} // main






























// LuM end of file
