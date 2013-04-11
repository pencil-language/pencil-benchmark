// UjoImro, 2013
// Experimental Research Code for the CARP Project

#include "errors.hpp"
#include "memory.hpp"

int main()
{
    carp::pool pool(1000);

    PRINT(pool.acquire(10));
    PRINT(pool.acquire(100));
    PRINT(pool.acquire(500));
    PRINT(pool.acquire(500));
    
}






























// LuM end of file
