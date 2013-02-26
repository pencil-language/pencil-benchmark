// UjoImro, 2013
// Experimental code for the CARP Project
// Copyright (c) RealEyes, 2013
// This is the response-map header exported for testing

typedef struct {
    int rows;
    int cols;
    int step;
    float * data;
} MatFloat; // struct MatFloat

typedef struct {
    int rows;
    int cols;
    int step;
    uint8_t * data;    
} MatChar; // struct MatChar

typedef struct {
    int x;
    int y;
} Point2i;

typedef enum { none, maxAbs, meanStd } NormalizationMethod;

MatFloat
generateResponseMap( 
    const MatChar image,
    const Point2i center,
    int mapSize, 
    int m_patchSize,
    int m_patch_line_size,
    int m_number_of_patches_per_line,
    NormalizationMethod postSVDNormalizationMethod,
    MatFloat m_wIn,
    MatFloat m_U,
    MatFloat m_wOut
    );


// LuM end of file
