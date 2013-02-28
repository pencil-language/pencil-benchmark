// UjoImro, 2013
// Experimental code for the CARP Project
// Copyright (c) RealEyes, 2013
// This is a c-implementation of the PCA->MLP response map calculation

#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>

#include "mlp_impl.h"

const int MAX_INT = 1 << sizeof(int) - 1;
const int false = (1!=1);
// const int NULL=0;


MatFloat
CreateMatFloat( int rows, int cols )
{
    MatFloat result; 
    
    result.data  = (float*)malloc( sizeof(float) * rows * cols );
    result.rows  = rows;
    result.cols  = cols;
    result.step  = cols;
    result.start = 0;

    return result;
}

void
copyTo( MatFloat input, MatFloat * output )
{
    assert(input.rows == output->rows);
    assert(input.cols == output->cols);
    
    int q, w;
    for ( q=0; q<input.rows; q++ )
	for ( w=0; w<=input.cols; w++ )
	    output->data[ q * output->cols + w + output->step ] =
		input.data[ q * input.cols + w + input.step ];
    
    return;
} // copyTo

void
transposeFloat( MatFloat input, MatFloat * output )
{
    int q,w;
    assert(output->rows == input.cols);
    assert(output->cols == input.rows);
    for (q=0; q<input.rows; q++)
	for (w=0; w<input.cols; w++)
	    output->data[ w * output->step + q + output->start ] 
		= input.data[ q * input.step + w + input.start ];

    return;
} // transposeFloat

float
meanChar( MatChar input )
{
    int q,w;
    float sum=0;

    for ( q=0; q<input.rows; q++ )
	for ( w=0; w<input.cols; w++ )
	    sum += input.data[ q * input.step + w + input.start ];
    
    return sum / ( input.rows * input.cols );
} // meanFloat

uint8_t
minChar( MatChar input )
{
    int q,w;
    uint8_t minvalue = 255;

    for ( q=0; q<input.rows; q++ )
	for ( w=0; w<input.cols; w++ )
	    minvalue = fmin( minvalue, input.data[ q * input.step + w + input.start ] );
    
    return minvalue;
} // minFloat

uint8_t
maxChar( MatChar input )
{
    int q,w;
    uint8_t maxvalue = 0;

    for ( q=0; q<input.rows; q++ )
	for ( w=0; w<input.cols; w++ )
	    maxvalue = fmax( maxvalue, input.data[ q * input.step + w + input.start ] );
    
    return maxvalue;
} // maxFloat


MatChar
GetBlockChar( MatChar self, int row_from, int row_to, int col_from, int col_to )
{    
    assert(row_from>=0);
    assert(col_from>=0);
    assert(row_from<row_to);
    assert(col_from<col_to);
    assert(row_to<=self.rows);
    assert(col_to<=self.cols);
    
    MatChar result;
    result.rows  = row_to - row_from;
    result.cols  = col_to - col_from;
    result.step  = self.step;
    result.start = self.start + row_from * self.step + col_from;
    result.data  = self.data;

    return result;
}

MatFloat
GetBlockFloat( MatFloat self, int row_from, int row_to, int col_from, int col_to )
{   
    assert(row_from>=0);
    assert(col_from>=0);
    assert(row_from<row_to);
    assert(col_from<col_to);
    assert(row_to<=self.rows);
    assert(col_to<=self.cols);

    MatFloat result;
    result.rows  = row_to - row_from;
    result.cols  = col_to - col_from;
    result.step  = self.step;
    result.start = self.start + row_from * self.step + col_from;
    result.data  = self.data;

    return result;
}

void
convertFromCharToFloat( MatChar from, float quotient, float shift, MatFloat * to )
{
    assert(from.rows == to->rows);
    assert(from.cols == to->cols);
    
    int q, w;
    for ( q=0; q<from.rows; q++ )
	for ( w=0; w<from.cols; w++ )
	    to->data[ q * to->step + w + to->start ] = 
		quotient * from.data[ q * from.step + w + from.start ] + shift;
    return;
}

MatFloat
reshapeFloat( MatFloat self, int new_rows )
{

    assert(self.cols == self.step);
    assert( (self.cols * self.rows) % new_rows == 0 );

    MatFloat result;
    result.rows = new_rows;
    result.cols = self.cols * self.rows / new_rows;
    result.step = result.cols;
    result.start = self.start;
    result.data = self.data;

    return result;
} // reshapeFloat

void
generatePatch( MatChar sample, int patch_size, MatFloat * result )
{   
    // cv::Mat_<pixel_type> result( sample.rows * patch_size, sample.cols - patch_size );
    assert(result->rows == sample.rows * patch_size);
    assert(result->cols == sample.cols - patch_size);
    assert(result->step == result->cols );

    // we are now processing line by line
    int q, col;
    for ( q = patch_size, col = 0;
        q < sample.cols;
        q++, col++ )
    {
        MatChar patch_image = GetBlockChar( sample, -MAX_INT, MAX_INT, col, col + patch_size );

        float sample_mean = meanChar(patch_image);
        float sample_max  = maxChar(patch_image) - sample_mean;
        float sample_min  = minChar(patch_image) - sample_mean;
	
        sample_max = fmax( abs( sample_max ), abs( sample_min ) );
        if (sample_max==0) sample_max = 1;
	
        MatFloat normalized_patch;
	
	convertFromCharToFloat( patch_image, 1.0/sample_max, -(1.0/sample_max)*sample_mean, &normalized_patch );
	
        MatFloat patch_line = reshapeFloat( normalized_patch ,normalized_patch.rows * normalized_patch.cols);
        MatFloat target = GetBlockFloat( *result, -MAX_INT, MAX_INT, col, col + 1 );
	
        copyTo(patch_line, &target);  //patch_line.copyTo( target );
    }

    return;
}


// returns alpha*A*B + beta * C
void 
gemmFloat( MatFloat A, MatFloat B, float alpha, MatFloat C, float beta, MatFloat * result )
{
    assert(A.rows == C.rows);
    assert(A.cols == B.rows); 
    assert(B.cols == C.cols);
    assert(C.rows == result->rows);
    assert(C.cols == result->cols);

    int q, w, e;
    float sum=0;
    for ( q=0; q<C.rows; q++ )
	for ( w=0; w<C.cols; w++ )
	{
	    sum = 0;
	    for ( e=0; e<A.cols; e++ )
		sum += alpha * A.data[ q * A.step + e + A.start ] * B.data[ e * B.step + w + B.start ] + beta * C.data[ q * C.step + w + C.start ];
	    result->data[ q * result->step + w + result->start ] = sum;
	}

    return;
}

void
expFloat( MatFloat input, MatFloat * output )
{
    assert(input.rows == output->rows);
    assert(input.cols == output->cols);

    int q, w;
    for ( q=0; q<input.rows; q++ )
	for ( w=0; w<input.cols; w++ )
	    output->data[ q * output->step + w + output->start ] = 
		exp(input.data[ q * input.step + w + input.start ]);

    return;
}

void
addFloat( MatFloat input, float val, MatFloat * output )
{
    assert(input.rows == output->rows);
    assert(input.cols == output->cols);

    int q, w;
    for ( q=0; q<input.rows; q++ )
	for ( w=0; w<input.cols; w++ )
	    output->data[ q * output->step + w + output->start ] = 
		val + input.data[ q * input.step + w + input.start ];

    return;
}

void 
divideFloat( float val, MatFloat input, MatFloat * output )
{
    assert(input.rows == output->rows);
    assert(input.cols == output->cols);

    int q, w;
    for ( q=0; q<input.rows; q++ )
	for ( w=0; w<input.cols; w++ )
	    output->data[ q * output->step + w + output->start ] = 
		val / input.data[ q * input.step + w + input.start ];

    return;
}

void
subtractFloat( MatFloat input, float val, MatFloat * output )
{
    assert(input.rows == output->rows);
    assert(input.cols == output->cols);

    int q, w;
    for ( q=0; q<input.rows; q++ )
	for ( w=0; w<input.cols; w++ )
	    output->data[ q * output->step + w + output->start ] = 
		input.data[ q * input.step + w + input.start ] - val;
    
    return;
}

MatFloat
evaluateSamples( 
    MatFloat m_wIn_gemm, 
    MatFloat m_all_patches,
    MatFloat m_all_bIns,
    MatFloat m_wOut_gemm,
    MatFloat m_all_bOuts
    )
{
    // PCA
    MatFloat xOut;
    gemmFloat( m_wIn_gemm, m_all_patches, -1., m_all_bIns, -1., &xOut ); // cv::gemm

    // input activator
    MatFloat e;
    expFloat( xOut, &e ); // cv::exp
    addFloat( e, 1.0, &xOut );
    divideFloat( 2.0, xOut, &e );
    subtractFloat( e, 1.0, &xOut );

    // MLP
    MatFloat dot;
    gemmFloat( m_wOut_gemm, xOut, -1., m_all_bOuts, -1., &dot ); // cv::gemm

    // output activator
    MatFloat result;
    expFloat( dot, &result ); // expFloat
    addFloat( result, 1.0, &dot ); // addFloat
    divideFloat( 1.0, dot, &result );

    return result;
} // evaluateSamples

float
GetValueFloat( MatFloat self, int row, int col )
{
    return self.data[ row * self.step + col + self.start ];
}

void
SetValueFloat( MatFloat self, int row, int col, float value )
{
    self.data[ row * self.step + col + self.start ] = value;
}


void
update( 
    int mapSize, 
    NormalizationMethod postSVDNormalizationMethod, 
    MatFloat m_wIn,
    MatFloat m_U,
    MatFloat m_wOut,
    int m_patchSize,

    // results 
    int * m_patch_line_size,
    int * m_number_of_patches_per_line,
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
    /* switch (postSVDNormalizationMethod) // none */
    /* { */
    /* case none: */
	
    /* { */
    /* 	// m_wIn_gemm = m_wIn.colRange( -MAX_INT, MAX_INT, 0, m_wIn.cols - 1 ) *m_U.t(); */
    /* 	MatFloat littleIn; */
    /* 	littleIn = GetBlockFloat( m_wIn, -MAX_INT, MAX_INT, 0, m_wIn.cols - 1 ); */
    /* 	MatFloat transpU = transposeFloat(m_U); */
    /* 	gemmFloat( littleIn, transpU, 1., *m_wIn_gemm, 0., m_wIn_gemm ); */
    /* } */
    /* break; */
    /* case maxAbs: */
    /* case meanStd: */
    /* 	// m_wIn_gemm = m_wIn.colRange(0,m_wIn.cols-1); */
    /* 	*m_wIn_gemm = GetBlockFloat( m_wIn, -MAX_INT, MAX_INT, 0, m_wIn.cols - 1 ); */
    /* 	break; */
    /* default: */
    /* 	assert(false); */
    /* } // postSVDNormalizationMethod */

    MatFloat littleIn;
    littleIn = GetBlockFloat( m_wIn, -MAX_INT, MAX_INT, 0, m_wIn.cols - 1 );
    MatFloat transpU;
    transposeFloat(m_U, &transpU);
    gemmFloat( littleIn, transpU, 1., *m_wIn_gemm, 0., m_wIn_gemm );


    // MatFloat bIn = m_wIn.col( m_wIn.cols - 1 );
    MatFloat bIn = GetBlockFloat( m_wIn, -MAX_INT, MAX_INT, m_wIn.cols - 1, m_wIn.cols);

    *m_wOut_gemm = GetBlockFloat( m_wOut, -MAX_INT, MAX_INT, 0, m_wOut.cols - 1 ); ///!!!!!!!!!! NO MORE TRANSPOSE .t();

    float bOut = GetValueFloat( m_wOut, 0, m_wOut.cols - 1 );

    *m_number_of_patches_per_line = 2 * m_currentMapSize + 1;
    // here we assume that the patches are always square-like
    int number_of_all_patches  = (*m_number_of_patches_per_line) * (*m_number_of_patches_per_line);
    *m_patch_line_size = 2 * m_patchSize + 1;
    int patch_size = (*m_patch_line_size) * (*m_patch_line_size);

    *m_all_bIns = CreateMatFloat( bIn.rows, number_of_all_patches );
    *m_all_bOuts = CreateMatFloat( 1, number_of_all_patches );
    *m_all_patches = CreateMatFloat( m_wIn_gemm->cols, number_of_all_patches );

    int q;
    for ( q=0; q < number_of_all_patches; q++)
    {
	target = GetBlockFloat( *m_all_bIns, -MAX_INT, MAX_INT, q, q+1 );
	copyTo(bIn, &target); // bIn.copyTo(target);
	SetValueFloat( *m_all_bOuts, 0, q, bOut);
    }    
} // update

void
copyToFloat( MatFloat src, MatFloat dst )
{

}

void
generateResponseMap(
    const MatChar image,
    const Point2i center,
    int mapSize, 
    mlp classifier, 
    MatFloat * result
    )
{

    MatFloat m_wIn_gemm;
    MatFloat m_all_bIns;
    MatFloat m_all_bOuts;
    MatFloat m_all_patches;
    MatFloat m_wOut_gemm;
    int m_patch_line_size;
    int m_number_of_patches_per_line;

    // make sure that we have the necessary matrices
    // calculated (this is a method; it is NOT thread safe!), but it is called for each classifier once
    update(
	// parameters
	mapSize,
	classifier.postSVDNormalizationMethod,
	classifier.m_wIn,
	classifier.m_U,
	classifier.m_wOut,
	classifier.m_patchSize,
	// results
	&m_patch_line_size,
	&m_number_of_patches_per_line,
	&m_wIn_gemm,
	&m_all_bIns,
	&m_all_bOuts,
	&m_all_patches,
	&m_wOut_gemm
	);

    int ncy, cy;
    for ( ncy = 0, cy = center.y - mapSize; cy < center.y + mapSize + 1; ++ncy, ++cy ) {
	MatChar  sample = GetBlockChar( image, cy - classifier.m_patchSize, cy + classifier.m_patchSize + 1 , center.x - mapSize - classifier.m_patchSize, center.x + mapSize + classifier.m_patchSize + 2 );
	
        MatFloat pack_patches;
	generatePatch( sample, m_patch_line_size, &pack_patches );

        MatFloat target = GetBlockFloat( m_all_patches, -MAX_INT, MAX_INT, ncy * m_number_of_patches_per_line, (ncy + 1) * m_number_of_patches_per_line );
        copyToFloat( pack_patches, target );
    }

    MatFloat evaluated = evaluateSamples(    
	m_wIn_gemm, 
	m_all_patches,
	m_all_bIns,
	m_wOut_gemm,
	m_all_bOuts
	);

    *result = reshapeFloat( evaluated, m_number_of_patches_per_line );
    return;    
} // generateResponseMap

int 
cvRound( float value )
{

}

void
calculateMaps( 
    int m_visibleLandmarks_size, 
    int m_mapSize, 
    MatChar alignedImage, 
    MatFloat shape, 
    mlp m_classifiers[], 
    MatFloat responseMaps[] )
{
    int q;
    for (q=0; q<m_visibleLandmarks_size; q++ )
    {
        /* const int idx = m_visibleLandmarks[q]; */
        /* assert(idx==q); */
	int idx = q;

        Point2i center;
	float shape_x;
	float shape_y;

	shape_x = GetValueFloat( shape, 2*idx, 0 );
	shape_y = GetValueFloat( shape, 2*idx+1, 0 );
	center.x = cvRound(shape_x);
	center.y = cvRound(shape_y);
        
	// responseMaps[q] = m_classifiers[idx].generateResponseMap( alignedImage, center, m_mapSize );
	generateResponseMap( alignedImage, center, m_mapSize, m_classifiers[idx], (&responseMaps[q]) );
    }

    return;
} // calculateMaps


// LuM end of file
