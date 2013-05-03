// UjoImro, 2013
// Experimental code for the CARP Project
// Copyright (c) RealEyes, 2013
// This is the response-map header exported for testing

#ifndef __MLP_IMPL__H__
#define __MLP_IMPL__H__

#include "cltypes.h"
// #include "allocator.hpp"

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

    
void
calculateMaps(
    void * self,
    int memory_segments[], // representing the start of the gang's memory segment
    int m_visibleLandmarks_size, 
    int m_mapSize, 
//    clMat /*float*/shape, 

    calcpackage inputs[] //,

    // // temporaries
    // cMat wIn,
    // cMat patches[],
    // cMat xOuts[],
    // cMat es[],

    // results
    // clMat /*float*/responseMaps[] );
    );
    
    
    void printMatFloat( void * self, clMat /*float*/mat, char * name );

    void printMatChar( void * self, clMat /*uint8_t*/ mat, char * name );
    

        
#endif /* __MLP_IMPL__H__ */

#ifdef __cplusplus
} // extern C
#endif // __cplusplus
// LuM end of file
