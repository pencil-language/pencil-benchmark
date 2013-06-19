// UjoImro, 2013

#include <vector>
#include <iostream>

template <class T0>
int printme( int & order, const T0 & t0 )
{
    std::cout << "[ " << order << ", " << t0 << "] " << std::endl;

    order++;

    return true;

}


template <class ...Args>
void
printall( Args... args )
{
    int q = 0;

    int results[] = { printme(q, args)... };

}
// printall


int main()
{
    printall( "1", 34, 1.2, "story" );

}

// LuM end of file
