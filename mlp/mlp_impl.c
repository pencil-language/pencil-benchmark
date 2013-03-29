// -*- c++ -*-
// UjoImro, 2013
// Experimental code for the CARP Project
// Copyright (c) RealEyes, 2013
// This is a c-implementation of the PCA->MLP response map calculation

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>

#include "mlp_impl.h"

const int MAX_INT = ~(1 << (8*sizeof(int) - 1));
const int false = (1!=1);
// const int NULL=0;

void     freeMLP( mlp * classifier );
MatFloat CreateMatFloat( int rows, int cols );
MatChar  CreateMatChar( int rows, int cols );
void     copyToFloat( MatFloat input, MatFloat * output );
void     transposeFloat( MatFloat input, MatFloat * output );
float    meanChar( MatChar input );
uint8_t  minChar( MatChar input );
uint8_t  maxChar( MatChar input );
MatChar  GetBlockChar( MatChar self, int row_from, int row_to, int col_from, int col_to );
MatFloat GetBlockFloat( MatFloat self, int row_from, int row_to, int col_from, int col_to );
void     convertFromCharToFloat( MatChar from, float quotient, float shift, MatFloat * to );
MatFloat reshapeFloat( MatFloat self, int new_rows );
void     generatePatch( MatChar sample, int patch_size, MatFloat * result );
void     freeMatFloat( MatFloat * mat );
void     freeMatChar( MatChar * mat );
void     gemmFloat( MatFloat A, MatFloat B, float alpha, MatFloat C, float beta, MatFloat * result );
void     expFloat( MatFloat input, MatFloat * output );
void     addFloat( MatFloat input, float val, MatFloat * output );
void     divideFloat( float val, MatFloat input, MatFloat * output );
void     subtractFloat( MatFloat input, float val, MatFloat * output );
float    GetValueFloat( MatFloat self, int row, int col );
void     SetValueFloat( MatFloat self, int row, int col, float value );
void     update( int mapSize, MatFloat m_wIn, MatFloat m_U, MatFloat m_wOut, int m_patchSize,  /* results */ int * m_patch_line_size, int * m_number_of_patches_per_line, MatFloat * m_wIn_gemm, MatFloat * m_all_bIns, MatFloat * m_all_bOuts, MatFloat * m_all_patches, MatFloat * m_wOut_gemm  );
void     generateResponseMap( const MatChar image, const Point2i center, int mapSize, mlp classifier, MatFloat * result  );
int      cvRound( float value );
void evaluateSamples( MatFloat m_wIn_gemm, MatFloat m_all_patches, MatFloat m_all_bIns, MatFloat m_wOut_gemm, MatFloat m_all_bOuts, /* result */ MatFloat * result );


void freeMLP( mlp * classifier )
{
    freeMatFloat(&classifier->m_wIn);
    freeMatFloat(&classifier->m_wOut);
    freeMatFloat(&classifier->m_U);
} // freeMLP

MatFloat
CreateMatFloat( int rows, int cols )
{
    MatFloat result; 
    assert(rows>0);
    assert(cols>0);

    result.data = NULL;
    result.data = (float*)malloc( sizeof(float) * rows * cols );
    assert(result.data);
    result.rows  = rows;
    result.cols  = cols;    
    result.step  = cols;
    result.start = 0;

    return result;
}

MatChar
CreateMatChar( int rows, int cols )
{
    MatChar result; 
    assert(rows>0);
    assert(cols>0);

    result.data = NULL;
    result.data = (uint8_t*)malloc( sizeof(uint8_t) * rows * cols );
    assert(result.data);
    result.rows  = rows;
    result.cols  = cols;
    result.step  = cols;
    result.start = 0;

    return result;
}


void
copyToFloat( MatFloat input, MatFloat * output )
{
    assert(input.rows == output->rows);
    assert(input.cols == output->cols);
    
    int q, w;
    for ( q=0; q<input.rows; q++ )
        for ( w=0; w<input.cols; w++ )
            output->data[ q * output->step + w + output->start ] =
        	input.data[ q * input.step + w + input.start ];
    
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
    assert(self.data);
    
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
    assert(self.data);

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
    // printf("generatePatch::00\n");
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
	// printf("generatePatch::01\n");
        MatChar patch_image = GetBlockChar( sample, 0, sample.rows, col, col + patch_size );
	// printf("generatePatch::02\n");
        float sample_mean = meanChar(patch_image);
	// printf("generatePatch::03\n");
        float sample_max  = maxChar(patch_image) - sample_mean;
	// printf("generatePatch::04\n");
        float sample_min  = minChar(patch_image) - sample_mean;
	// printf("generatePatch::05\n");
	
        sample_max = fmax( abs( sample_max ), abs( sample_min ) );
	// printf("generatePatch::06\n");
        if (sample_max==0) sample_max = 1;
	
	{
	    MatFloat normalized_patch = CreateMatFloat( patch_image.rows, patch_image.cols );
	    // printf("generatePatch::07\n");
	
	    convertFromCharToFloat( patch_image, 1.0/sample_max, -(1.0/sample_max)*sample_mean, &normalized_patch );
	
	    // printf("generatePatch::08\n");
	    MatFloat patch_line = reshapeFloat( normalized_patch ,normalized_patch.rows * normalized_patch.cols);
	    // printf("generatePatch::09\n");
	    MatFloat target = GetBlockFloat( *result, 0, result->rows, col, col + 1 );
	    // printf("generatePatch::10\n");
	    copyToFloat(patch_line, &target);  //patch_line.copyTo( target );
	    // printf("generatePatch::11\n");
	    freeMatFloat(&normalized_patch);
	    // printf("generatePatch::12\n");
	}

    } // for q in patch_size..sample.cols

    return;
} // generatePatch

void
freeMatFloat( MatFloat * mat )
{
    assert(mat->data);
    assert(mat->start==0);
    free(mat->data);
    mat->data  = NULL;
    mat->rows  = 0;
    mat->cols  = 0;
    mat->step  = 0;
    mat->start = 0;
    return;
} // freeMatFloat

void 
freeMatChar( MatChar * mat )
{
    assert(mat->data);
    assert(mat->start==0);
    free(mat->data);
    mat->data  = NULL;
    mat->rows  = 0;
    mat->cols  = 0;
    mat->step  = 0;
    mat->start = 0;
    return;
} // freeMatChar


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

    if ( fabs(beta) > 0.000001 ) {
        for ( q=0; q<C.rows; q++ )
            for ( w=0; w<C.cols; w++ )
            {	    
                sum = 0;
                for ( e=0; e<A.cols; e++ )
                {              
                    float y = A.data[ q * A.step + e + A.start ] * B.data[ e * B.step + w + B.start ] - c;
                    float t = sum + y;
                    c = (t - sum) - y;
                    sum = t;              
                }
                
                result->data[ q * result->step + w + result->start ] = alpha * sum  + beta * C.data[ q * C.step + w + C.start ];
            }
    }
    else
    {
        for ( q=0; q<C.rows; q++ )
            for ( w=0; w<C.cols; w++ )
            {	    
                sum = 0;
                for ( e=0; e<A.cols; e++ )
                {              
                    float y = A.data[ q * A.step + e + A.start ] * B.data[ e * B.step + w + B.start ] - c;
                    float t = sum + y;
                    c = (t - sum) - y;
                    sum = t;              
                }
                
                result->data[ q * result->step + w + result->start ] = alpha * sum;
            }        
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

void
evaluateSamples( 
    MatFloat m_wIn_gemm, 
    MatFloat m_all_patches,
    MatFloat m_all_bIns,
    MatFloat m_wOut_gemm,
    MatFloat m_all_bOuts,
    // result
    MatFloat * result
    )
{
    // PCA
    MatFloat xOut = CreateMatFloat( m_wIn_gemm.rows, m_all_patches.cols );
    gemmFloat( m_wIn_gemm, m_all_patches, -1., m_all_bIns, -1., &xOut ); // cv::gemm

    // input activator
    MatFloat e = CreateMatFloat(xOut.rows, xOut.cols);
    expFloat( xOut, &e ); // cv::exp
    addFloat( e, 1.0, &xOut );
    divideFloat( 2.0, xOut, &e );
    subtractFloat( e, 1.0, &xOut );

    // MLP
    MatFloat dot = CreateMatFloat( m_wOut_gemm.rows, xOut.cols );
    gemmFloat( m_wOut_gemm, xOut, -1., m_all_bOuts, -1., &dot ); // cv::gemm

    // output activator
    
    expFloat( dot, result ); // expFloat
    addFloat( *result, 1.0, &dot ); // addFloat
    divideFloat( 1.0, dot, result );

    // printf("evaluateSamples::01\n");    
    freeMatFloat(&dot);
    // printf("evaluateSamples::02\n");
    freeMatFloat(&e);
    // printf("evaluateSamples::03\n");
    freeMatFloat(&xOut);
    // printf("evaluateSamples::04\n");

    return;
} // evaluateSamples

float
GetValueFloat( MatFloat self, int row, int col )
{
    return self.data[ row * self.step + col + self.start ];
    // return -1231.;
}

void
SetValueFloat( MatFloat self, int row, int col, float value )
{
    self.data[ row * self.step + col + self.start ] = value;
}


void
update( 
    int mapSize, 
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

    MatFloat littleIn;
    littleIn = GetBlockFloat( m_wIn, 0, m_wIn.rows, 0, m_wIn.cols - 1 );
    {
	MatFloat transpU = CreateMatFloat( m_U.cols, m_U.rows );
	transposeFloat(m_U, &transpU);
        
	MatFloat tmpMap  = CreateMatFloat( littleIn.rows, transpU.cols );
        
        gemmFloat( littleIn, transpU, 1., tmpMap, 0., m_wIn_gemm );

        freeMatFloat(&tmpMap);
	freeMatFloat(&transpU);
    }
    // MatFloat bIn = m_wIn.col( m_wIn.cols - 1 );
    MatFloat bIn = GetBlockFloat( m_wIn, 0, m_wIn.rows, m_wIn.cols - 1, m_wIn.cols);
     
    *m_wOut_gemm = GetBlockFloat( m_wOut, 0, m_wOut.rows, 0, m_wOut.cols - 1 ); ///!!!!!!!!!! NO MORE TRANSPOSE .t();
    float bOut = GetValueFloat( m_wOut, 0, m_wOut.cols - 1 );
    *m_number_of_patches_per_line = 2 * m_currentMapSize + 1;
    // here we assume that the patches are always square-like
    int number_of_all_patches  = (*m_number_of_patches_per_line) * (*m_number_of_patches_per_line);
    *m_patch_line_size = 2 * m_patchSize + 1;
    int patch_size = (*m_patch_line_size) * (*m_patch_line_size);
    assert( m_all_bIns->rows == bIn.rows );
    assert( m_all_bIns->cols == number_of_all_patches );
    assert( m_all_bOuts->rows == 1 );
    assert( m_all_bOuts->cols == number_of_all_patches );
    assert( m_all_patches->rows == m_wIn_gemm->cols );
    assert( m_all_patches->cols == number_of_all_patches );
    
    int q;
    for ( q=0; q < number_of_all_patches; q++)
    {
        target = GetBlockFloat( *m_all_bIns, 0, m_all_bIns->rows, q, q+1 );
        copyToFloat(bIn, &target); // bIn.copyTo(target);
        SetValueFloat( *m_all_bOuts, 0, q, bOut);
    }    
    return;
} // update

void
generateResponseMap(
    const MatChar image,
    const Point2i center,
    int mapSize, 
    mlp classifier, 
    MatFloat * result
    )
{

    // printf("ici05a\n");
    // printf("classifier.m_wIn.rows = %d\n", classifier.m_wIn.rows );
    // printf("classifier.m_U.rows = %d\n", classifier.m_U.rows );
    MatFloat m_wIn_gemm = CreateMatFloat( classifier.m_wIn.rows, classifier.m_U.rows );
    // printf("ici05b\n");
    MatFloat m_all_bIns = CreateMatFloat( classifier.m_wIn.rows, (2 * mapSize + 1) * (2 * mapSize + 1) );
    // printf("ici05c\n");
    MatFloat m_all_bOuts = CreateMatFloat( 1, (2 * mapSize + 1) * (2 * mapSize + 1) );
    // printf("ici05d\n");
    MatFloat m_all_patches = CreateMatFloat( m_wIn_gemm.cols, (2 * mapSize + 1) * (2 * mapSize + 1) );
    // printf("ici05e\n");
    MatFloat m_wOut_gemm;
    int m_patch_line_size;
    int m_number_of_patches_per_line;
    
    // printf("ici05f\n");
    // make sure that we have the necessary matrices
    // calculated (this is a method; it is NOT thread safe!), but it is called for each classifier once
    update(
    	// parameters
    	mapSize,
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

    // printf("ici05g\n");
    int ncy, cy;
    for ( ncy = 0, cy = center.y - mapSize; cy < center.y + mapSize + 1; ++ncy, ++cy ) {
        // printf("mapSize = %d\n", mapSize);
        // printf("classifier.m_patchSize =%d\n", classifier.m_patchSize );
        // printf("center.y = %d\n", center.y );        
        // printf("cy - classifier.m_patchSize = %d\n", cy - classifier.m_patchSize);
        
    	MatChar  sample = GetBlockChar( image, cy - classifier.m_patchSize, cy + classifier.m_patchSize + 1 , center.x - mapSize - classifier.m_patchSize, center.x + mapSize + classifier.m_patchSize + 2 );	
    	{
    	    // printf("ici05g0\n");
    	    MatFloat pack_patches = CreateMatFloat( sample.rows * m_patch_line_size, sample.cols - m_patch_line_size );
    	    // printf("ici05g1\n");
    	    generatePatch( sample, m_patch_line_size, &pack_patches );
    	    // printf("ici05g2\n");
	    
    	    MatFloat target = GetBlockFloat( m_all_patches, 0, m_all_patches.rows, ncy * m_number_of_patches_per_line, (ncy + 1) * m_number_of_patches_per_line );
    	    // printf("ici05g3\n");
    	    copyToFloat( pack_patches, &target );
    	    // printf("ici05g4\n");
    	    freeMatFloat(&pack_patches);
    	    // printf("ici05g5\n");
    	} // unnamed block
    } // for

    /* printf("ici05h\n"); */
    MatFloat evaluated = CreateMatFloat( m_wOut_gemm.rows, m_all_patches.cols );

    // printf("ici05i\n");
    evaluateSamples(
    	m_wIn_gemm,
    	m_all_patches,
    	m_all_bIns,
    	m_wOut_gemm,
    	m_all_bOuts,
    	// results
    	&evaluated
    	);

    // printf("ici05j\n");
    MatFloat toreturn = reshapeFloat( evaluated, m_number_of_patches_per_line );
    // printf( "toreturn.rows = %d\n", toreturn.rows );
    // printf( "toreturn.cols = %d\n", toreturn.cols );
    assert( toreturn.rows == result->rows );
    assert( toreturn.cols == result->cols );

    // printf("result->cols = %d\n", result->cols );
    // printf("result->rows = %d\n", result->rows );
    copyToFloat(toreturn, result);

    // printf("ici05k\n");
    freeMatFloat(&m_wIn_gemm);
    // printf("ici05l\n");  
    freeMatFloat(&m_all_bIns);
    // printf("ici05m\n");
    freeMatFloat(&m_all_bOuts);
    // printf("ici05n\n");
    freeMatFloat(&m_all_patches);
    // printf("ici05o\n");
    freeMatFloat(&evaluated);
    // printf("ici05p\n");
    return;
} // generateResponseMap

int 
cvRound( float value )
{
    return (int)(value + (value >= 0 ? 0.5 : -0.5));
} // cvRound

void
calculateMaps( 
    int m_visibleLandmarks_size, 
    int m_mapSize, 
    MatChar alignedImage, 
    MatFloat shape, 
    mlp m_classifiers[], 
    // results
    MatFloat * responseMaps[] )
{
    // printf("calculateMaps started\n");    
    int q;
    for (q=0; q<m_visibleLandmarks_size; q++ )
    {
	// printf("processing patch %d/%d\n", q, m_visibleLandmarks_size );
        /* const int idx = m_visibleLandmarks[q]; */
        /* assert(idx==q); */
	int idx = q;

        Point2i center;
	float shape_x;
	float shape_y;

	// printf("ici01\n");
	shape_x = GetValueFloat( shape, 2*idx, 0 );
	// printf("ici02\n");
	shape_y = GetValueFloat( shape, 2*idx+1, 0 );
	// printf("ici03\n");
	center.x = cvRound(shape_x);
	// printf("ici04\n");
	center.y = cvRound(shape_y);
	// printf("ici05\n");
        
	// responseMaps[q] = m_classifiers[idx].generateResponseMap( alignedImage, center, m_mapSize );
	generateResponseMap( alignedImage, center, m_mapSize, m_classifiers[idx], (&(*responseMaps)[q]) );
	// printf("ici06\n");
    }

    // printf("calculateMaps finished\n");
    return;
} // calculateMaps


// LuM end of file
