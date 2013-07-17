// UjoImro, 2013
// Experimental code for the CARP Project
// Copyright (c) RealEyes, 2013
// This is the response-map header exported for testing

#ifndef __MLP_IMPL__H__
#define __MLP_IMPL__H__

#include "cltypes.h"

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus
    typedef enum { none, maxAbs, meanStd } NormalizationMethod;

    typedef unsigned char uint8_t;
        
    typedef struct {
        int rows;
        int cols;
        int step;
        int start;    
        float * data;
    } MatFloat; // struct MatFloat

    typedef struct {
        int rows;
        int cols;
        int step;
        int start;    
        uint8_t * data;    
    } MatChar; // struct MatChar

    typedef struct {
        int m_patchSize;      /*!< \brief Radius like patch size, the true size of the patch is [(2*patchSize+1) x (2*patchSize+1)] */
        MatFloat m_wIn; /*!< \brief */
        MatFloat m_wOut; /*!< \brief  */
        MatFloat m_U; /*!< \brief */
        int hidden_num;
        double rho2;
        // NormalizationMethod preSVDNormalizationMethod /*= none*/;
        // NormalizationMethod postSVDNormalizationMethod;
    } mlp; // struct 

    MatFloat CreateMatFloat( int rows, int cols );

    MatChar CreateMatChar( int rows, int cols );

    void freeMatFloat( MatFloat * mat );

    void freeMatChar( MatChar * mat );

    void freeMLP( mlp * classifier );

    void
    calculateMaps( 
        int m_visibleLandmarks_size, 
        int m_mapSize, 
        MatChar alignedImage, 
        MatFloat shape, 
        mlp m_classifiers[], 
        // results
        MatFloat * responseMaps[] );    



#ifdef __cplusplus
}
#endif // __cplusplus
// LuM end of file

#endif /* __MLP_IMPL__H__ */
