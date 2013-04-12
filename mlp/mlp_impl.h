// UjoImro, 2013
// Experimental code for the CARP Project
// Copyright (c) RealEyes, 2013
// This is the response-map header exported for testing

#ifndef __MLP_IMPL__H__
#define __MLP_IMPL__H__

#include "cltypes.h"
#include "allocator.hpp"

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus


    void
    calculateMaps(
        void * self,
        void * allocator,
        int m_visibleLandmarks_size, 
        int m_mapSize, 
        cMat /*uint8_t*/  alignedImage, 
        cMat /*float*/shape, 
        int m_patchSizes[],      /*!< \brief Radius like patch size, the true size of the patch is [(2*patchSize+1) x (2*patchSize+1)] */
        cMat /*float*/m_wIns[], /*!< \brief */
        cMat /*float*/m_wOuts[], /*!< \brief  */
        cMat /*float*/m_Us[], /*!< \brief */
        
        // results
        cMat /*float*/* responseMaps[] );
    
        
#endif /* __MLP_IMPL__H__ */

#ifdef __cplusplus
} // extern C
#endif // __cplusplus
// LuM end of file
