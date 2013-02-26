// UjoImro, 2013
// Experimental code for the CARP Project
// Copyright (c) RealEyes, 2013

#include <stdlib.h>
#include <stdint.h>

const int MAX_INT = 1 << sizeof(int) - 1;
const int false = (1!=1);
// const int NULL=0;

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


MatFloat transposeFloat( MatFloat input )
{

}

MatChar GetBlockChar( MatChar self, int row_from, int row_to, int col_from, int col_to )
{
    
}

MatFloat GetBlockFloat( MatFloat self, int row_from, int row_to, int col_from, int col_to )
{
    
}

MatFloat
generatePatch( MatChar sample, int m_patch_line_size )
{

}

MatFloat
reshape( MatFloat self, int ignored, int new_col_size )
{

}

MatFloat evaluateSamples( 
    MatFloat m_wIn_gemm, 
    MatFloat m_all_patches,
    MatFloat m_all_bIns,
    MatFloat m_wOut_gemm,
    MatFloat m_all_bOuts
    )
{
    // PCA
    MatFloat xOut;
    gemmFloat( m_wIn_gemm, m_all_patches, -1., m_all_bIns, -1., xOut ); // cv::gemm

    // input activator
    MatFloat e;
    expFloat( xOut, e ); // cv::exp
    addFloat( e, 1.0, xOut );
    divideFloat( 2.0, xOut, e );
    subtractFloat( e, 1.0, xOut );

    // MLP
    MatFloat dot;
    gemmFloat( m_wOut_gemm, xOut, -1., m_all_bOuts, -1., dot ); // cv::gemm

    // output activator
    MatFloat result;
    expFloat( dot, result ); // expFloat
    addFloat( result, 1.0, dot ); // addFloat
    divideFloat( 1.0, dot, result );

    return result;
} // evaluateSamples

float GetValueFloat( MatFloat self, int x, int y )
{
    
}

void SetValueFloat( MatFloat self, int x, int y, float value )
{
    
}


MatFloat CreateMatFloat( int rows, int cols )
{

}

void update( 
    int mapSize, 
    NormalizationMethod postSVDNormalizationMethod, 
    MatFloat m_wIn,
    MatFloat m_U,
    MatFloat m_wOut,
    int m_patchSize,

    // results 
    MatFloat * m_wIn_gemm,
    MatFloat * m_all_bIns,
    MatFloat * m_all_bOuts,
    MatFloat * m_all_patches,
    MatFloat * m_wOut_gemm
    )
{
    MatFloat target;

    int m_currentMapSize = mapSize;


    // extract to the creation
    switch (postSVDNormalizationMethod)
    {
    case none:
	
    {
	// m_wIn_gemm = m_wIn.colRange( -MAX_INT, MAX_INT, 0, m_wIn.cols - 1 ) *m_U.t();
	MatFloat littleIn;
	littleIn = GetBlockFloat( m_wIn, -MAX_INT, MAX_INT, 0, m_wIn.cols - 1 );
	MatFloat transpU = transposeFloat(m_U);
	gemmFloat( littleIn, transpU, 1., NULL, 1., *m_wIn_gemm );
    }
    break;
    case maxAbs:
    case meanStd:
	// m_wIn_gemm = m_wIn.colRange(0,m_wIn.cols-1);
	*m_wIn_gemm = GetBlockFloat( m_wIn, -MAX_INT, MAX_INT, 0, m_wIn.cols - 1 );
	break;
    default:
	assert(false);
    } // postSVDNormalizationMethod

    // MatFloat bIn = m_wIn.col( m_wIn.cols - 1 );
    MatFloat bIn = GetBlockFloat( m_wIn, -MAX_INT, MAX_INT, m_wIn.cols - 1, m_wIn.cols);

    *m_wOut_gemm = GetBlockFloat( m_wOut, -MAX_INT, MAX_INT, 0, m_wOut.cols - 1 ); ///!!!!!!!!!! NO MORE TRANSPOSE .t();

    float bOut = GetValue(m_wOut, 0, m_wOut.cols - 1);

    int m_number_of_patches_per_line = 2 * m_currentMapSize + 1;
    // here we assume that the patches are always square-like
    int number_of_all_patches  = m_number_of_patches_per_line * m_number_of_patches_per_line;
    int m_patch_line_size = 2 * m_patchSize + 1;
    int patch_size = m_patch_line_size * m_patch_line_size;

    *m_all_bIns = CreateMatFloat( bIn.rows, number_of_all_patches );
    *m_all_bOuts = CreateMatFloat( 1, number_of_all_patches );
    *m_all_patches = CreateMatFloat( m_wIn_gemm->cols, number_of_all_patches );

    int q;
    for ( q=0; q < number_of_all_patches; q++)
    {
	target = GetBlockFloat( *m_all_bIns, -MAX_INT, MAX_INT, q, q+1 );
	copyTo(bIn, target); // bIn.copyTo(target);
	SetValueFloat( *m_all_bOuts, 0, q, bOut);
    }    
} // update

void copyToFloat( MatFloat src, MatFloat dst )
{

}

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
    )
{

    MatFloat m_wIn_gemm;
    MatFloat m_all_bIns;
    MatFloat m_all_bOuts;
    MatFloat m_all_patches;
    MatFloat m_wOut_gemm;

    // make sure that we have the necessary matrices
    // calculated (this is a method; it is NOT thread safe!), but it is called for each classifier once
    update(
	// parameters
	mapSize, 
	postSVDNormalizationMethod, 
	m_wIn,
	m_U,
	m_wOut,
	m_patchSize, 
	// results
	&m_wIn_gemm,
	&m_all_bIns,
	&m_all_bOuts,
	&m_all_patches,
	&m_wOut_gemm
	);

    int ncy, cy;
    for ( ncy = 0, cy = center.y - mapSize; cy < center.y + mapSize + 1; ++ncy, ++cy ) {
	MatChar  sample = GetBlockChar( image, cy - m_patchSize, cy + m_patchSize + 1 , center.x - mapSize - m_patchSize, center.x + mapSize + m_patchSize + 2 );
        MatFloat pack_patches = generatePatch( sample, m_patch_line_size );

        MatFloat target = GetBlockFloat( m_all_patches, -MAX_INT, MAX_INT, ncy * m_number_of_patches_per_line, (ncy + 1) * m_number_of_patches_per_line );
        copyToFloat( pack_patches, target );
    }

    MatFloat result = evaluateSamples(    
	m_wIn_gemm, 
	m_all_patches,
	m_all_bIns,
	m_wOut_gemm,
	m_all_bOuts
	);

    return reshape( result, 0, m_number_of_patches_per_line );

} // generateResponseMap


// LuM end of file
