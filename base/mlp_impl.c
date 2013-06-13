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
const int gangsize = 640;

// const int NULL=0;

void     freeMLP( mlp * classifier );
MatFloat CreateMatFloat( int rows, int cols );
MatChar  CreateMatChar( int rows, int cols );
void     copyToFloat( MatFloat input, MatFloat * output );
void     transposeFloat( MatFloat input, MatFloat * output );
void     transposeFloatGang( MatFloat input, int localid, MatFloat * output );
float    meanChar( MatChar input );
uint8_t  minChar( MatChar input );
uint8_t  maxChar( MatChar input );
MatChar  GetBlockChar( MatChar self, int row_from, int row_to, int col_from, int col_to );
MatFloat GetBlockFloat( MatFloat self, int row_from, int row_to, int col_from, int col_to );
void     convertFromCharToFloat( MatChar from, float quotient, float shift, MatFloat * to );
MatFloat reshapeFloat( MatFloat self, int new_rows );
void     freeMatFloat( MatFloat * mat );
void     freeMatChar( MatChar * mat );
void     gemmFloat( MatFloat A, MatFloat B, float alpha, MatFloat C, float beta, MatFloat * result );
void     expFloat( MatFloat input, MatFloat * output );
void     addFloat( MatFloat input, float val, MatFloat * output );
void     divideFloat( float val, MatFloat input, MatFloat * output );
void     subtractFloat( MatFloat input, float val, MatFloat * output );
float    GetValueFloat( MatFloat self, int row, int col );
void     SetValueFloat( MatFloat * self, int row, int col, float value );
void     generateResponseMap( const MatChar image, const Point2i center, int mapSize, mlp classifier, MatFloat * result  );
int      cvRound( float value );
float    dotProductDirDir( MatFloat A, MatFloat B );
float    dotProductTransDir( MatFloat A, MatFloat B );
void     normalizeSample( MatChar image, MatFloat * result );
void     printMatFloat( MatFloat mat, char * name );

void     printMatFloat( MatFloat mat, char * name )
{
  printf("%s = [\n", name);

  int q,w;

  for (q=0; q<mat.rows; q++)
  {
    printf("[ ");
    for( w=0; w<mat.cols; w++)
    {
      printf( "%f, ", GetValueFloat(mat, q, w) );
    }
    printf(" ]\n");
  }
  
  printf("]\n");
  
  return;
} // printMatFloat

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
    assert(input.data);
    assert(output);
    assert(output->data);    
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

void
transposeFloatGang( MatFloat input, int localid, MatFloat * output )
{
    assert(output->rows == input.cols);
    assert(output->cols == input.rows);
    
    int worksize = input.cols * input.rows;
    int work = 0;
            
    for ( work=0; work < (worksize/gangsize) + 1; work++) {
        int index = work * gangsize + localid;
        if (index>=worksize) continue;
        int q = index / input.cols;
        int w = index % input.cols;
              
        output->data[ w * output->step + q + output->start ] 
            = input.data[ q * input.step + w + input.start ];
    }
} // transposeFloatGang


float
meanChar( MatChar input )
{
    assert(input.data);    
    int q,w;
    float sum=0;
    float c = 0; // kahan summation
    
    for ( q=0; q<input.rows; q++ )
      for ( w=0; w<input.cols; w++ )
      {
        float y = input.data[ q * input.step + w + input.start ] - c;
        float t = sum + y;
        c = (t - sum) - y;        
        sum = t;        
      }
    
    return sum / ( input.rows * input.cols );
} // meanFloat

uint8_t
minChar( MatChar input )
{
    assert(input.data);    
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
    assert(input.data);
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
//    printf("col_to = %d\n", col_to );
//    printf("self.cols = %d\n", self.cols );
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
gemmFloatDirDirDir( MatFloat A, MatFloat B, float alpha, MatFloat C, float beta, MatFloat * result )
{
    assert(A.rows == C.rows);
    assert(A.cols == B.rows); 
    assert(B.cols == C.cols);
    assert(C.rows == result->rows);
    assert(C.cols == result->cols);

    int q, w, e;
    float sum=0.;
    float c;

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
    else // NOT fabs(beta) > 0.000001
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
} // gemmFloatDirDirDir

void 
gemmFloatDirDirDirGang( MatFloat A, MatFloat B, float alpha, MatFloat C, float beta, int localid, MatFloat * result )
{
    assert(A.rows == C.rows);
    assert(A.cols == B.rows); 
    assert(B.cols == C.cols);
    assert(C.rows == result->rows);
    assert(C.cols == result->cols);

    int q, w, e;
    float sum=0.;
    float c;    
    if ( fabs(beta) > 0.000001 ) {
        int worksize = C.rows * C.cols;
        int work = 0;
        
        for ( work=0; work < (worksize/gangsize) + 1; work++) {            
            int index = work * gangsize + localid;
            if (index>=worksize) continue;
            int q = index / C.cols;
            int w = index % C.cols;

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
    else // NOT fabs(beta) > 0.000001
    {
        int worksize = C.rows * C.cols;
        int work = 0;
        
        for ( work=0; work < (worksize/gangsize) + 1; work++) {
            int index = work * gangsize + localid;
            if (index>=worksize) continue;
            int q = index / C.cols;
            int w = index % C.cols;
            
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
} // gemmFloatDirDirDirGang


void 
gemmFloatDirTransDirGang( MatFloat A, MatFloat B, float alpha, MatFloat C, float beta, int localid, MatFloat * result )
{
    assert(A.rows == C.rows);
    assert(A.cols == B.cols); 
    assert(B.rows == C.cols);
    assert(C.rows == result->rows);
    assert(C.cols == result->cols);

    int q, w, e;
    float sum=0.;
    float c;    
    if ( fabs(beta) > 0.000001 ) {
        int worksize = C.rows * C.cols;
        int work = 0;
        
        for ( work=0; work < (worksize/gangsize) + 1; work++) {            
            int index = work * gangsize + localid;
            if (index>=worksize) continue;
            int q = index / C.cols;
            int w = index % C.cols;

            sum = 0;
            for ( e=0; e<A.cols; e++ )
            {              
                float y = A.data[ q * A.step + e + A.start ] * B.data[ w * B.step + e + B.start ] - c;
                float t = sum + y;
                c = (t - sum) - y;
                sum = t;              
            }
                
            result->data[ q * result->step + w + result->start ] = alpha * sum  + beta * C.data[ q * C.step + w + C.start ];
        }
    }
    else // NOT fabs(beta) > 0.000001
    {
        int worksize = C.rows * C.cols;
        int work = 0;
        
        for ( work=0; work < (worksize/gangsize) + 1; work++) {
            int index = work * gangsize + localid;
            if (index>=worksize) continue;
            int q = index / C.cols;
            int w = index % C.cols;
            
            sum = 0;
            for ( e=0; e<A.cols; e++ )
            {              
                float y = A.data[ q * A.step + e + A.start ] * B.data[ w * B.step + e + B.start ] - c;
                float t = sum + y;
                c = (t - sum) - y;
                sum = t;              
            }
            
            result->data[ q * result->step + w + result->start ] = alpha * sum;
        }
    }
    
    return;
} // gemmFloatDirTransDirGang


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

float
GetValueFloat( MatFloat self, int row, int col )
{
    return self.data[ row * self.step + col + self.start ];
    // return -1231.;
}

void
SetValueFloat( MatFloat * self, int row, int col, float value )
{
    self->data[ row * self->step + col + self->start ] = value;
    return;    
}



float dotProductDirDir( MatFloat A, MatFloat B )
{
  assert( A.cols == 1 );
  assert( B.cols == 1 );
  assert( A.rows == B.rows );

  float result = 0.;
  float c = 0.;  

  int q;
  for (q=0; q<A.rows; q++) {
      float y = A.data[ q * A.step + A.start ] * B.data[ q * B.step + B.start ] - c;
      float t = result + y;
      c = (t - result) - y;
      result = t ;
  }
  
  return result;
}

float dotProductTransDir( MatFloat A, MatFloat B )
{
  assert( A.rows == 1 );
  assert( B.cols == 1 );
  assert( A.cols == B.rows );

  float result = 0.;
  float c = 0.;  

  int q;
  for (q=0; q<A.cols; q++) {
      float y = A.data[ q + A.start ] * B.data[ q * B.step + B.start ] - c;
      float t = result + y;
      c = (t - result) - y;
      result = t ;
  }
  
  return result;
}

void normalizeSample( MatChar image, MatFloat * result )
{
  assert(result->cols == image.cols);
  assert(result->rows == image.rows);
    
  float sampleMean = meanChar(image);
  float sampleMin  = minChar(image);
  float sampleMax  = maxChar(image);

  sampleMax -= sampleMean;
  sampleMin -= sampleMean;

  sampleMax = fmax( fabs(sampleMin), fabs(sampleMax));

  if (sampleMax == 0.0) sampleMax = 1.0;

  convertFromCharToFloat( image, 1.0/sampleMax, -(1.0/sampleMax)*sampleMean, result );

  *result = reshapeFloat( *result, image.rows * image.cols );
 
  return;
} // normalizeSample

void
generateResponseMap(
    const MatChar image,
    const Point2i center,
    int mapSize, 
    mlp classifier, 
    MatFloat * result
    )
{

  assert(result->rows == 2 * mapSize + 1);
  assert(result->cols == 2 * mapSize + 1);

  MatFloat wIn_A = GetBlockFloat( classifier.m_wIn, 0, classifier.m_wIn.rows, 0, classifier.m_wIn.cols - 1 );
  MatFloat wIn = CreateMatFloat( wIn_A.rows, classifier.m_U.rows );
  MatFloat bIn = GetBlockFloat( classifier.m_wIn, 0, classifier.m_wIn.rows, classifier.m_wIn.cols - 1, classifier.m_wIn.cols );
  MatFloat wOut_tmp = GetBlockFloat( classifier.m_wOut, 0, classifier.m_wOut.rows, 0, classifier.m_wOut.cols - 1 );

  {
      int localid = 0;
      for ( localid=0; localid<gangsize; localid++ )           
          gemmFloatDirTransDirGang( wIn_A, classifier.m_U, 1.0, wIn, 0.0, localid, &wIn );
  }
  
  {
      int localid = 0;
      for ( localid=0; localid<gangsize; localid++ ) {
          
          int work = 0;      
          int worksize = (2 * mapSize + 1) * (2 * mapSize + 1);
          
          float bOut = GetValueFloat( classifier.m_wOut, 0, classifier.m_wOut.cols - 1 );
          
          int ncy=0; 
          int cy=0;
          int ncx=0;
          int cx=0;
          
          for ( work=0; work < (worksize / gangsize) + 1; work++ )
          {
              int index = work * gangsize + localid;
              if (index>=worksize) continue;
              int ncy = index / (2 * mapSize + 1);
              int cy  = center.y - mapSize + ncy;
              int ncx = index % (2 * mapSize + 1);
              int cx  = center.x - mapSize + ncx;
          
              MatChar  imagePatch = GetBlockChar( image, cy - classifier.m_patchSize, cy + classifier.m_patchSize + 1, cx - classifier.m_patchSize, cx + classifier.m_patchSize + 1 );
              MatFloat patch = CreateMatFloat( imagePatch.rows, imagePatch.cols );
          
              normalizeSample(imagePatch, &patch);
          
              MatFloat xOut = CreateMatFloat( bIn.rows, bIn.cols );
          
              gemmFloatDirDirDir( wIn, patch, -1.0, bIn, -1.0, &xOut );
          
              MatFloat e = CreateMatFloat(xOut.rows, xOut.cols);
          
              expFloat( xOut, &e );
          
              addFloat( e, 1.0, &xOut );
              divideFloat( 2.0, xOut, &e);
              addFloat( e, -1.0, &xOut);
          
              SetValueFloat( result, ncy, ncx, 1./( 1. + exp(- dotProductTransDir(wOut_tmp, xOut) - bOut ) ) );
          
              freeMatFloat(&e);
              freeMatFloat(&xOut);
              freeMatFloat(&patch);
          } // for localid 
      } // for gangsize
  } // localid block
  
//  assert(false);
  
  freeMatFloat(&wIn);
  // end of classic impl
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
    MatChar alignedImages[], 
    MatFloat shapes[], 
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
	shape_x = GetValueFloat( shapes[q], 2*idx, 0 );
	// printf("ici02\n");
	shape_y = GetValueFloat( shapes[q], 2*idx+1, 0 );
	// printf("ici03\n");
	center.x = cvRound(shape_x);
	// printf("ici04\n");
	center.y = cvRound(shape_y);
	// printf("ici05\n");
        
	// responseMaps[q] = m_classifiers[idx].generateResponseMap( alignedImage, center, m_mapSize );
	generateResponseMap( alignedImages[q], center, m_mapSize, m_classifiers[idx], (&(*responseMaps)[q]) );
	// printf("ici06\n");
    }

    // printf("calculateMaps finished\n");
    return;
} // calculateMaps


// LuM end of file
