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

    
    void calculateMaps( char * self, int memory_segments[], int m_visibleLandmarks_size, int m_mapSize, calcpackage inputs[] );
    
    void printMatFloat( void * self, clMat /*float*/mat, char * name );
    
    void printMatChar( void * self, clMat /*uint8_t*/ mat, char * name );
    
    float GetValueFloat( void * self, clMat /*float*/smat, int row, int col );
    
    uint8_t GetValueChar( void * self, clMat /* uint8_t*/ smat, int row, int col );

    void copyToFloat( void * self, clMat /*float*/input, clMat /*float*/ output );    

    void transposeFloat( void * self, clMat /*float*/input, clMat /*float*/ output );
    
    void transposeFloatGang( void * self, clMat /*float*/input, int localid, clMat /*float*/ output );
    
    float meanChar( void * self, clMat /*uint8_t*/  input );    

    uint8_t minChar( void * self, clMat /*uint8_t*/  input );
    
    uint8_t maxChar( void * self, clMat /*uint8_t*/  input );
    
    clMat /*uint8_t*/ GetBlockChar( void * self, clMat /*uint8_t*/  smat, int row_from, int row_to, int col_from, int col_to );
    
    clMat GetBlockFloat( void * self, clMat /*float*/smat, int row_from, int row_to, int col_from, int col_to );
    
    void convertFromCharToFloat( void * self, clMat /*uint8_t*/  from, float quotient, float shift, clMat /*float*/ to );    

    clMat reshapeFloat( void * self, clMat /*float*/smat, int new_rows );

    void gemmFloatDirDirDir( void * self, clMat /*float*/A, clMat /*float*/B, float alpha, clMat /*float*/C, float beta, clMat /*float*/ result );
    
    void gemmFloatDirDirDirGang( void * self, clMat /*float*/A, clMat /*float*/B, float alpha, clMat /*float*/C, float beta, int localid, clMat /*float*/ result );
    
    void gemmFloatDirTransDirGang( void * self, clMat /*float*/A, clMat /*float*/B, float alpha, clMat /*float*/C, float beta, int localid, clMat /*float*/ result );
    
    void expFloat( void * self, clMat /*float*/input, clMat /*float*/ output );
    
    void addFloat( void * self, clMat /*float*/input, float val, clMat /*float*/ output );
    
    void divideFloat( void * self, float val, clMat /*float*/input, clMat /*float*/ output );
    
    void subtractFloat( void * self, clMat /*float*/input, float val, clMat /*float*/ output );
    
    float GetValueFloat( void * self, clMat /*float*/smat, int row, int col );
    
    uint8_t GetValueChar( void * self, clMat /* uint8_t*/ smat, int row, int col );
    
    void SetValueFloat( void * self, clMat /*float*/ smat, int row, int col, float value );
    
    float dotProductDirDir( void * self, clMat /*float*/A, clMat /*float*/B );
    
    float dotProductTransDir( void * self, clMat /*float*/A, clMat /*float*/B );
        
    void normalizeSample( void * self, clMat /*uint8_t*/  image, clMat /*float*/ * result );
    
    void generateResponseMap( void * self, const clMat /*uint8_t*/  image, const Point2i center, int mapSize, int m_patchSize, clMat /*float*/m_wOut, clMat wIn, clMat bIn, clVector /*cMat float*/ patches, clVector /*cMat float*/ xOuts, clVector /*cMat float*/ es, clMat /*float*/ result );

    int cvRound( float value );
    
        
#ifdef __cplusplus
} // extern C
#endif // __cplusplus

#endif /* __MLP_IMPL__H__ */
// LuM end of file
