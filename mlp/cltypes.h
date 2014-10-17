// UjoImro, 2013

#ifndef __CLTYPES__H__
#define __CLTYPES__H__

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

    typedef struct {
        int rows;
        int cols;
        int step;
        int start;
    } clMat;

    typedef struct {
        int m_patchSize;
        clMat /* uint8_t */ alignedImage;
        clMat /* int32_t */ shape;
        clMat /* float */ m_wOut;
        clMat /* float */ wIn;
        clMat /* float */ bIn;        
    } calcinput;

    typedef struct {
        int x;
        int y;
    } Point2i; // struct Point2i

    typedef struct {
        // results
        clMat /*float*/ responseMap;        
    } calcoutput;
            
    typedef struct {
        calcinput  input;
        calcoutput output;
    } calcpackage;
#ifdef __cplusplus
} // extern C
#endif // __cplusplus

#endif
