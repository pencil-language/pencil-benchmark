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
#include "mlp_impl.h"

extern int EF_ALIGNMENT = 0;
extern int EF_PROTECT_BELOW = 0;
extern int EF_PROTECT_FREE = 0;
extern int EF_ALLOW_MALLOC_0 = 1;
extern int EF_FILL = 1922;


namespace { struct hack_t; }
MatChar  convertCVToMatChar ( const cv::Mat_<uint8_t> & input );
mlp *    convertHackToMlp ( const hack_t & hack );
void     freeClassifiers( mlp * classifiers[], int size );
MatChar  convertCVToMatChar ( const cv::Mat_<uint8_t> & input );
MatFloat convertCVToMatFloat (  const cv::Mat_<double> & input );
cv::Mat_<double> convertMatFloatToCV( MatFloat input );



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



mlp *
convertHackToMlp ( const hack_t & hack )
{
    assert(hack.m_visibleLandmarks_size==hack.m_classifiers.size());
    assert(hack.m_visibleLandmarks_size==hack.responseMaps.size());
        
    mlp * result;
    result = reinterpret_cast<mlp*>( malloc( sizeof(mlp) * hack.m_visibleLandmarks_size ) );

    // we export each classifier
    for (int q=0; q<hack.m_visibleLandmarks_size; q++)
    {
        result[q].m_patchSize = hack.m_classifiers[q].m_patchSize;
        result[q].m_wIn       = convertCVToMatFloat(hack.m_classifiers[q].m_wIn);
        result[q].m_wOut      = convertCVToMatFloat(hack.m_classifiers[q].m_wOut);
        result[q].m_U         = convertCVToMatFloat(hack.m_classifiers[q].m_U);
        result[q].hidden_num  = hack.m_classifiers[q].hidden_num;
        result[q].rho2        = hack.m_classifiers[q].rho2;
    } // for q in m_visibleLandmarks_size

    return result;
} // convertHackToMlp

void
freeClassifiers( mlp * classifiers[], int size )
{
    mlp * result = *classifiers;    
    for (int q=0; q<size; q++ )
        freeMLP( &(result[q]) );

    free(*classifiers);
    *classifiers=NULL;

    return;    
} // freeClassifiers


MatChar convertCVToMatChar ( const cv::Mat_<uint8_t> & input )
{
    MatChar result = CreateMatChar( input.rows, input.cols );

    for ( int q=0; q<input.rows; q++)
        for ( int w=0; w<input.cols; w++ )
            result.data[ q * result.step + w + result.start ] = input(q,w);
    
    return result;    
} // convertCVToMatChar

MatFloat convertCVToMatFloat (  const cv::Mat_<double> & input )
{
    MatFloat result = CreateMatFloat( input.rows, input.cols );
    
    for ( int q=0; q<input.rows; q++)
        for ( int w=0; w<input.cols; w++ )
            result.data[ q * result.step + w + result.start ] = input(q,w);
    
    return result;    
} // convertCVToMatFloat


cv::Mat_<double> convertMatFloatToCV( MatFloat input )
{
    cv::Mat_<double> result( input.rows, input.cols );
    
    for ( int q=0; q<input.rows; q++)
        for ( int w=0; w<input.cols; w++ )
            result(q,w) = input.data[ q * input.step + w + input.start ];
    
    return result;
} // convertMatFloatToCV


void allocateResponseMaps( int mapSize, int size, MatFloat * responseMaps[] )
{
    *responseMaps = (MatFloat*)malloc( sizeof(MatFloat) * size );
    assert(*responseMaps);    
    MatFloat * result = *responseMaps;    

    for ( int q=0; q<size; q++ )
        result[q] = CreateMatFloat( 2 * mapSize + 1, 2 * mapSize + 1 );
    
    return;
}

void freeResponseMaps( MatFloat * responseMaps[], int size )
{
    MatFloat * result = *responseMaps;    
    assert(result);
    for ( int q=0; q<size; q++ )
        freeMatFloat(&(result[q]));

    free(result);
    *responseMaps=NULL;
    return;    
}


int main()
{
    static conductor_t conductor;

    for ( conductor.importer >> BOOST_SERIALIZATION_NVP(conductor.id);
          ((conductor.id != -1) and (conductor.id != 10));
          // conductor.id != -1;
          conductor.importer >> BOOST_SERIALIZATION_NVP(conductor.id)
        )
    {
        PRINT(conductor.id);
        conductor.importer >> BOOST_SERIALIZATION_NVP(conductor.hack);

        // here comes the function call
        {
            // preparing the inputs
            MatChar alignedImage = convertCVToMatChar(conductor.hack.alignedImage);
            MatFloat shape = convertCVToMatFloat(conductor.hack.shape);
            mlp * m_classifiers = convertHackToMlp(conductor.hack);
            MatFloat * responseMaps;
            allocateResponseMaps( conductor.hack.m_mapSize, conductor.hack.m_visibleLandmarks_size, &responseMaps );
                        
            calculateMaps(
                conductor.hack.m_visibleLandmarks_size,
                conductor.hack.m_mapSize,
                alignedImage,
                shape,
                m_classifiers,
                &responseMaps
                );

            // releasing the inputs
            std::cout << "cpp 01" << std::endl;
            freeMatChar(&alignedImage);
            std::cout << "cpp 02" << std::endl;
            freeMatFloat(&shape);
            std::cout << "cpp 03" << std::endl;
            freeClassifiers(&m_classifiers, conductor.hack.m_classifiers.size());
            std::cout << "CHECK ME !!! cpp 04" << std::endl;
            
            // converting the outputs
            std::vector< cv::Mat_<double> > calculatedResults;
            for (int q=0; q<conductor.hack.m_visibleLandmarks_size; q++)
            {
                cv::Mat_<double> nextResult;
                nextResult = convertMatFloatToCV( responseMaps[q] );
                calculatedResults.push_back(nextResult);                
            }
            
            // testing the output
            for (int q=0; q<conductor.hack.m_visibleLandmarks_size; q++)
            {
                std::cout << "cv::norm( conductor.hack.responseMaps[" << q << "] ) = "
                          << cv::norm( conductor.hack.responseMaps[q] ) << std::endl; 
                
                std::cout << "cv::norm( calculatedResults[" << q << "] ) = "
                          << cv::norm( calculatedResults[q] ) << std::endl; 

                std::cout << "cv::norm( conductor.hack.responseMaps[" << q << "] - calculatedResults[" << q << "] ) = "
                          << cv::norm( conductor.hack.responseMaps[q] - calculatedResults[q] ) << std::endl; 
            }
            
            // releasing the outputs
            freeResponseMaps( &responseMaps, conductor.hack.m_visibleLandmarks_size );
            std::cout << "cpp 05" << std::endl;

        }
        // here comes the test
        // PRINT(cv::norm( ));
    }
    
    
    //conductor.importer >> BOOST_SERIALIZATION_NVP(conductor.hack);
    
    return EXIT_SUCCESS;
}



// LuM end of file
