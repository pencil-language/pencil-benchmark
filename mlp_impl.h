// UjoImro, 2013
// Experimental code for the CARP Project
// Copyright (c) RealEyes, 2013
// This is the response-map header exported for testing

typedef struct {
    int rows;
    int cols;
    int step;
    int start;    
    float * data;
} MatFloat; // struct MatFloat

typedef struct {
    int rows;
    int cols;
    int step;
    int start;    
    uint8_t * data;    
} MatChar; // struct MatChar

typedef struct {
    int x;
    int y;
} Point2i;

typedef enum { none, maxAbs, meanStd } NormalizationMethod;

typedef struct {
    int m_patchSize;      /*!< \brief Radius like patch size, the true size of the patch is [(2*patchSize+1) x (2*patchSize+1)] */
    MatFloat m_wIn; /*!< \brief */
    MatFloat m_wOut; /*!< \brief  */
    MatFloat m_U; /*!< \brief */
    int hidden_num;
    double rho2;
    NormalizationMethod preSVDNormalizationMethod /*= none*/;
    NormalizationMethod postSVDNormalizationMethod;

} mlp; // struct 

// LuM end of file
