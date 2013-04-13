// UjoImro, 2013
// Experimental code for the CARP Project
// Copyright (c) RealEyes, 2013
// This is the response-map header exported for testing

#ifndef __ALLOCATOR__H__
#define __ALLOCATOR__H__

#include "cltypes.h"

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus


    typedef struct {
        int rows;
        int cols;
        int step;
        int start;    
    } cMat; // struct cMat

    
    typedef struct {
        int x;
        int y;
    } Point2i; // struct Point2i

    typedef enum { none, maxAbs, meanStd } NormalizationMethod;

    typedef struct {
        int m_patchSize;      /*!< \brief Radius like patch size, the true size of the patch is [(2*patchSize+1) x (2*patchSize+1)] */
        cMat /*float*/m_wIn; /*!< \brief */
        cMat /*float*/m_wOut; /*!< \brief  */
        cMat /*float*/m_U; /*!< \brief */
        int hidden_num;
        double rho2;
        // NormalizationMethod preSVDNormalizationMethod /*= none*/;
        // NormalizationMethod postSVDNormalizationMethod;
    } mlp; // struct 

    typedef struct {
        int m_patchSize;
        cMat alignedImage;
        cMat m_wIn;
        cMat m_wOut;
        cMat m_Us;
    } calcpackage; // struct 
    
    cMat /*float*/CreateMatFloat( void * self, void * allocator, int rows, int cols );

    cMat /*uint8_t*/CreateMatChar /*uint8_t*/ ( void * self, void * allocator, int rows, int cols );

    void freeMatFloat( void * self, void * allocator, cMat /*float*/* mat );

    void freeMatChar( void * self, void * allocator, cMat /*float*/* mat );
    
    void freecMat /*uint8_t*/ ( void * self, void * allocator, cMat /*uint8_t*/  * mat );

    void freeMLP( void * self, void * allocator, mlp * classifier );

    
            

#ifdef __cplusplus
} // extern C
#endif // __cplusplus


#endif /* __ALLOCATOR__H__ */

// LuM end of file
