// UjoImro, 2013
// Experimental Research Code for the CARP Project

#ifndef __CARP__UTILITY__HPP__
#define __CARP__UTILITY__HPP__

#include <iostream>
#include <opencv2/core/core.hpp>
#include <boost/preprocessor.hpp>

const int KiB=1024;
const int MiB=1024*KiB;
const int GiB=1024*KiB;

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

template <class T0>
bool
push( std::vector<T0> & self, T0 t0 )
{
    self.push_back(t0);
    return true;
} // push

namespace carp {

    template<class T0, class... Types>
    std::vector<T0>
    make_vector( T0 t0, Types... inputs )
    {
        std::vector<T0> result;
        result.push_back(t0);        
        bool err[] = { ::push( result, static_cast<T0>(inputs) )... };
        return result;
    } // make_vector


    template <class... Types>
    std::vector<std::string>
    string_vector( Types... inputs)
    {
        return make_vector<std::string>(inputs...);        
    } // string_vector
    
} // namespace carp


#endif /* __CARP__UTILITY__HPP__ */
