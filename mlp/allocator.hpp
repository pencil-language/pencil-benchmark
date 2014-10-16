// UjoImro, 2013
// Experimental code for the CARP Project
// Copyright (c) RealEyes, 2013
// This is the response-map header exported for testing

#ifndef __ALLOCATOR__H__
#define __ALLOCATOR__H__

#include "cltypes.h"
#include "memory.hpp"

namespace carp {
    typedef struct {
        int m_patchSize;      /* Radius-like patch size, the true size of the patch is [(2*patchSize+1) x (2*patchSize+1)] */
        clMat /*float*/m_wIn;
        clMat /*float*/m_wOut;
        clMat /*float*/m_U;
        int hidden_num;
        double rho2;
    } mlp;

    clMat /*float*/CreateMatFloat( carp::memory::allocator & pool, int rows, int cols );

    clMat /*uint8_t*/CreateMatChar /*uint8_t*/ ( carp::memory::allocator & pool, int rows, int cols );
};

#endif
