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
    } clMat; // struct clMatFloat

    typedef struct {
        int size;
        int step;
        int start;
    } clVector; // struct clVector

    typedef struct {
        int m_patchSize;
        clMat /* uint8_t */ alignedImage;
        clMat /* int32_t */ shape;        
//        clMat /* float */ m_wIn;
        clMat /* float */ m_wOut;
//        clMat /* float */ m_U;
        clMat /* float */ wIn;
        clMat /* float */ bIn;        
    } calcinput; // struct 

    typedef struct {
        int x;
        int y;
    } Point2i; // struct Point2i

    typedef struct {
        // temporaries
//        clVector /* clMat[] */ patches;
//        clVector /* clMat[] */ xOuts;
//        clVector /* clMat[] */ es;        
    } calctemp; // struct
    
    typedef struct {
        // results
        clMat /*float*/ responseMap;        
    } calcoutput;
            
    typedef struct {
        calcinput  input;
        //      calctemp   tmp;
        calcoutput output;
    } calcpackage;

    typedef struct {
        float shift;
        float stride;        
    } normalization;
    
        

#ifdef __cplusplus
} // extern C
#endif // __cplusplus


#endif /* __CLTYPES__H__ */

// LuM end of file
