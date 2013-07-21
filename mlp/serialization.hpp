// UjoImro, 2013

#ifndef SERIALIZATION__HPP__
#define SERIALIZATION__HPP__

#include <boost/archive/xml_iarchive.hpp>

#include "mlp.hpp"
#include "mlp_impl.h"

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
    
class conductor_t {
public:
    int id;
    hack_t hack;
    std::ifstream dumpStream;
    boost::archive::xml_iarchive importer;
    
public:
    
    conductor_t() : id(0), dumpStream("response_dumps.xml", std::ios::in | std::ios::binary ), importer(dumpStream)
        { }; // conductor_t
    
}; // conductor_t

#endif /* SERIALIZATION__HPP__ */
