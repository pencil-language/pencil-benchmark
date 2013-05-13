// UjoImro, 2013
// Experimental Research Code for the CARP Project

#ifndef __CARP__UTILITY__HPP__
#define __CARP__UTILITY__HPP__

#include <iostream>
#include <opencv2/core/core.hpp>
#include <boost/preprocessor.hpp>

#define PRINT(var)  std::cout << "debug: " << BOOST_PP_STRINGIZE(var) << " = " << var << std::endl

template<class T0>
void
print_image( cv::Mat_<T0> input, std::string name )
{
    std::cout << name << " = [" << std::endl;
    
    for (int q=0; q<input.rows; q++)
    {
        std::cout << "[ ";        
        for (int w=0; w<input.cols; w++)
        {
            std::cout << input(q,w);
            if (w<input.cols-1)
                std::cout << ", ";
            else
                std::cout << " ";
        }

        if (q<input.rows-1)
            std::cout << "], " << std::endl;
        else
            std::cout << "] " << std::endl;
    }

    std::cout << "]" << std::endl;    
} // print_image



#endif /* __CARP__UTILITY__HPP__ */
