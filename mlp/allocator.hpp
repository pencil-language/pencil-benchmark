// UjoImro, 2013
// Experimental code for the CARP Project
// Copyright (c) RealEyes, 2013
// This is the response-map header exported for testing

#ifndef __ALLOCATOR__H__
#define __ALLOCATOR__H__

#include "cltypes.h"
#include "memory.hpp"

namespace carp {

    typedef enum { none, maxAbs, meanStd } NormalizationMethod;

    typedef struct {
        int m_patchSize;      /*!< \brief Radius like patch size, the true size of the patch is [(2*patchSize+1) x (2*patchSize+1)] */
        clMat /*float*/m_wIn; /*!< \brief */
        clMat /*float*/m_wOut; /*!< \brief  */
        clMat /*float*/m_U; /*!< \brief */
        int hidden_num;
        double rho2;
        // NormalizationMethod preSVDNormalizationMethod /*= none*/;
        // NormalizationMethod postSVDNormalizationMethod;
    } mlp; // struct 

    clMat /*float*/CreateMatFloat( carp::memory::allocator & pool, int rows, int cols );

    clMat /*uint8_t*/CreateMatChar /*uint8_t*/ ( carp::memory::allocator & pool, int rows, int cols );

    clVector createVectorMatFloat( carp::memory::allocator & pool, int size, int rows, int cols );

    void freeMatFloat( carp::memory::allocator & pool, clMat /*float*/* mat );

    void freeMatChar( carp::memory::allocator & pool, clMat /*float*/* mat );
    
    void freeclMat /*uint8_t*/ ( carp::memory::allocator & pool, clMat /*uint8_t*/  * mat );

    void freeMLP( carp::memory::allocator & pool, mlp * classifier );

    void freeVector( carp::memory::allocator & pool, clVector * vec );

    clVector /* <clMat> */ CreateVectorMat( carp::memory::allocator & pool, int nb_elements );

    void freeVectorMat( carp::memory::allocator & pool, clVector * vec );

}; // namespace carp


#endif /* __ALLOCATOR__H__ */

// LuM end of file
