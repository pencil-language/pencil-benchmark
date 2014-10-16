// UjoImro, 2013

#ifndef SERIALIZATION__HPP__
#define SERIALIZATION__HPP__

#include <boost/filesystem.hpp>
#include <boost/archive/xml_iarchive.hpp>

#include "mlp.hpp"
#include "mlp_impl_pencil.h"

namespace carp {   

    struct hack_t;

    MatChar  convertCVToMatChar ( const cv::Mat_<uint8_t> & input );
    mlp *    convertHackToMlp ( const hack_t & hack );
    void     freeClassifiers( mlp * classifiers[], int size );
    MatFloat convertCVToMatFloat (  const cv::Mat_<double> & input );
    cv::Mat_<double> convertMatFloatToCV( MatFloat input );
    void     allocateResponseMaps( int mapSize, int size, MatFloat * responseMaps[] );
    void     freeResponseMaps( MatFloat * responseMaps[], int size );
     
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

    
    class guard_t {
    public:
        guard_t() {
            if ( !boost::filesystem::exists( "pool/response_dumps.xml" ) )
            {
                throw std::runtime_error("Can't find `pool/response_dumps.xml' file!");                
            } // if
        } // guard_t
    }; // class guard_t

    
    class conductor_t {
    public:
        guard_t guard;        
        int id;
        std::ifstream dumpStream;
        boost::archive::xml_iarchive importer;
    
    public:
    
        conductor_t() : guard(), id(0), dumpStream("pool/response_dumps.xml", std::ios::in | std::ios::binary ), importer(dumpStream)
            { }; // conductor_t
    
    }; // conductor_t

    
        
    
} // namespace carp

    
#endif /* SERIALIZATION__HPP__ */
