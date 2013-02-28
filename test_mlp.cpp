// UjoImro, 2013
// Experimental code for the CARP Project
// Copyright (c) RealEyes, 2013
// This version tests the responseMap calculation with input dumps

#include <string>
#include <stdlib.h>
#include <opencv2/core/core.hpp>
#include <boost/preprocessor.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include "cast.h"
#include "mlp.hpp"


namespace {
     
    struct hack_t {
        int m_visibleLandmarks_size;
        int m_mapSize;
        
        cv::Mat_<double> shape;
        cv::Mat_<uint8_t> alignedImage;
        std::vector<gel::MLP<double> > m_classifiers;
        std::vector<cv::Mat_<double> > responseMaps;
                
        template <class MT0>
        void serialize( MT0 & archiver, unsigned int ) {
            GEL_EXPORT_VAR( archiver, m_visibleLandmarks_size );
            GEL_EXPORT_VAR( archiver, m_mapSize );
            GEL_EXPORT_VAR( archiver, shape );
            GEL_EXPORT_VAR( archiver, alignedImage );
            GEL_EXPORT_VAR( archiver, m_classifiers );
            GEL_EXPORT_VAR( archiver, responseMaps );            
        } // serialize
    
    }; // struct hack_t
    
    class conductor_t {
    public:
        int id;
        hack_t hack;
        std::ifstream dumpStream;
        boost::archive::binary_iarchive importer;

    public:
        
        conductor_t() : id(0), dumpStream("response_dumps.bin", std::ios::in | std::ios::binary ), importer(dumpStream) {
            
        }; // conductor_t
        
    }; // conductor_t

} // unnamed namespace 


int main()
{
    static conductor_t conductor;

    conductor.importer >> BOOST_SERIALIZATION_NVP(conductor.id);
    PRINT(conductor.id);
    
    conductor.importer >> BOOST_SERIALIZATION_NVP(conductor.hack);
    
    return EXIT_SUCCESS;    
}



// LuM end of file
