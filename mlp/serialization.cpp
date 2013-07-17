// UjoImro, 2013
// Experimental code for the CARP Project
// Copyright (c) RealEyes, 2013

#include <assert.h>
#include <iostream>
#include <opencv2/core/core.hpp>

#include "serialization.hpp"

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


// LuM end of file

