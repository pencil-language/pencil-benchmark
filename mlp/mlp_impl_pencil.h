// UjoImro, 2013
// Experimental code for the CARP Project
// Copyright (c) RealEyes, 2013
// This is the response-map header exported for testing

#ifndef __MLP_IMPL__H__
#define __MLP_IMPL__H__

#include "cltypes.h"

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif
    typedef struct {
        int rows;
        int cols;
        int step;
        int start;
        float * data;
    } MatFloat;

    typedef struct {
        int rows;
        int cols;
        int step;
        int start;
        uint8_t * data;
    } MatChar;

    typedef struct {
        int m_patchSize;      /*!< \brief Radius like patch size, the true size of the patch is [(2*patchSize+1) x (2*patchSize+1)] */
        MatFloat m_wIn;
        MatFloat m_wOut;
        MatFloat m_U;
        int hidden_num;
        double rho2;
    } mlp;

    MatFloat CreateMatFloat( int rows, int cols );

    MatChar CreateMatChar( int rows, int cols );

    void freeMatFloat( MatFloat * mat );

    void freeMatChar( MatChar * mat );

    void freeMLP( mlp * classifier );

    void calculateMaps( int m_visibleLandmarks_size
                      , int m_mapSize
                      , MatChar alignedImage
                      , MatFloat shape
                      , mlp m_classifiers[]
                      , MatFloat * responseMaps[]
                      );



#ifdef __cplusplus
}
#endif // __cplusplus

#endif
