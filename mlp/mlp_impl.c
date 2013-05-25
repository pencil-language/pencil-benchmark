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
const int gangsize = 32;

// const int NULL=0

clMat GetMatFromVector( void * self, clVector /*clMat*/ vec, int idx )
{
    assert(self);
    assert(idx>=0);
    assert(idx<vec.size);
    
    return ((clMat*)self)[ vec.start + vec.step * idx ];
} // GetMatFromVector

void
printMatFloat( void * self, clMat /*float*/mat, char * name )
{
  printf("%s = [\n", name);

  int q,w;

  for (q=0; q<mat.rows; q++)
  {
    printf("[ ");
    for( w=0; w<mat.cols; w++)
    {
        printf( "%f, ", GetValueFloat(self, mat, q, w) );
    }
    printf(" ]\n");
  }
  
  printf("]\n");
  
  return;
} // printMatFloat

void
printMatChar( void * self, clMat /*uint8_t*/ mat, char * name )
{
  printf("%s = [\n", name);

  int q,w;

  for (q=0; q<mat.rows; q++)
  {
    printf("[ ");
    for( w=0; w<mat.cols; w++)
    {
        printf( "%d, ", (int)GetValueChar(self, mat, q, w) );
    }
    printf(" ]\n");
  }
  
  printf("]\n");
  
  return;
} // printMatChar


void
copyToFloat( void * self, clMat /*float*/input, clMat /*float*/ output )
{
    assert(input.rows == output.rows);
    assert(input.cols == output.cols);
    
    int q, w;
    for ( q=0; q<input.rows; q++ )
        for ( w=0; w<input.cols; w++ )
            ((float*)self)[ q * output.step + w + output.start ] =
                ((float*)self)[ q * input.step + w + input.start ];
    
    return;
} // copyTo

void
transposeFloat( void * self, clMat /*float*/input, clMat /*float*/ output )
{
    int q,w;
    assert(output.rows == input.cols);
    assert(output.cols == input.rows);
    for (q=0; q<input.rows; q++)
        for (w=0; w<input.cols; w++)
            ((float*)self)[ w * output.step + q + output.start ] 
        	= ((float*)self)[ q * input.step + w + input.start ];

    return;
} // transposeFloat

void
transposeFloatGang( void * self, clMat /*float*/input, int localid, clMat /*float*/ output )
{
    assert(output.rows == input.cols);
    assert(output.cols == input.rows);
    
    int worksize = input.cols * input.rows;
    int work = 0;
            
    for ( work=0; work < (worksize/gangsize) + 1; work++) {
        int index = work * gangsize + localid;
        if (index>=worksize) continue;
        int q = index / input.cols;
        int w = index % input.cols;
              
        ((float*)self)[ w * output.step + q + output.start ] 
            = ((float*)self)[ q * input.step + w + input.start ];
    }
} // transposeFloatGang


float
meanChar( void * self, clMat /*uint8_t*/  input )
{    
    int q,w;
    float sum=0;
    float c = 0; // kahan summation
    
    for ( q=0; q<input.rows; q++ )
      for ( w=0; w<input.cols; w++ )
      {
        float y = ((uint8_t*)self)[ q * input.step + w + input.start ] - c;
        float t = sum + y;
        c = (t - sum) - y;        
        sum = t;        
      }
    
    return sum / ( input.rows * input.cols );
} // meanFloat

uint8_t
minChar( void * self, clMat /*uint8_t*/  input )
{    
    int q,w;
    uint8_t minvalue = 255;

    for ( q=0; q<input.rows; q++ )
        for ( w=0; w<input.cols; w++ )
            minvalue = fmin( minvalue, ((uint8_t*)self)[ q * input.step + w + input.start ] );
    
    return minvalue;
} // minFloat

uint8_t
maxChar( void * self, clMat /*uint8_t*/  input )
{
    int q,w;
    uint8_t maxvalue = 0;

    for ( q=0; q<input.rows; q++ )
        for ( w=0; w<input.cols; w++ )
            maxvalue = fmax( maxvalue, ((uint8_t*)self)[ q * input.step + w + input.start ] );
    
    return maxvalue;
} // maxFloat


clMat /*uint8_t*/ 
GetBlockChar( void * self, clMat /*uint8_t*/  smat, int row_from, int row_to, int col_from, int col_to )
{    
    assert(row_from>=0);
    assert(col_from>=0);
    assert(row_from<row_to);
    assert(col_from<col_to);
    assert(row_to<=smat.rows);
    assert(col_to<=smat.cols);
    
    clMat /*uint8_t*/  result;
    result.rows  = row_to - row_from;
    result.cols  = col_to - col_from;
    result.step  = smat.step;
    result.start = smat.start + row_from * smat.step + col_from;

    return result;
}

clMat
GetBlockFloat( void * self, clMat /*float*/smat, int row_from, int row_to, int col_from, int col_to )
{  
    assert(row_from>=0);
    assert(col_from>=0);
    assert(row_from<row_to);
    assert(col_from<col_to);
    assert(row_to<=smat.rows);
    assert(col_to<=smat.cols);

    clMat /*float*/result;
    result.rows  = row_to - row_from;
    result.cols  = col_to - col_from;
    result.step  = smat.step;
    result.start = smat.start + row_from * smat.step + col_from;

    return result;
}

void
convertFromCharToFloat( void * self, clMat /*uint8_t*/  from, float quotient, float shift, clMat /*float*/ to )
{
    assert(from.rows == to.rows);
    assert(from.cols == to.cols);
    
    int q, w;
    for ( q=0; q<from.rows; q++ )
        for ( w=0; w<from.cols; w++ )
            ((float*)self)[ q * to.step + w + to.start ] = 
        	quotient * ((uint8_t*)self)[ q * from.step + w + from.start ] + shift;
    return;
}

clMat
reshapeFloat( void * self, clMat /*float*/smat, int new_rows )
{

    assert(smat.cols == smat.step);
    assert( (smat.cols * smat.rows) % new_rows == 0 );

    clMat /*float*/result;
    result.rows = new_rows;
    result.cols = smat.cols * smat.rows / new_rows;
    result.step = result.cols;
    result.start = smat.start;

    return result;
} // reshapeFloat

// returns alpha*A*B + beta * C
void 
gemmFloatDirDirDir( void * self, clMat /*float*/A, clMat /*float*/B, float alpha, clMat /*float*/C, float beta, clMat /*float*/ result )
{
    assert(A.rows == C.rows);
    assert(A.cols == B.rows); 
    assert(B.cols == C.cols);
    assert(C.rows == result.rows);
    assert(C.cols == result.cols);

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
                    float y = ((float*)self)[ q * A.step + e + A.start ] * ((float*)self)[ e * B.step + w + B.start ] - c;
                    float t = sum + y;
                    c = (t - sum) - y;
                    sum = t;              
                }
                
                ((float*)self)[ q * result.step + w + result.start ] = alpha * sum  + beta * ((float*)self)[ q * C.step + w + C.start ];
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
                    float y = ((float*)self)[ q * A.step + e + A.start ] * ((float*)self)[ e * B.step + w + B.start ] - c;
                    float t = sum + y;
                    c = (t - sum) - y;
                    sum = t;              
                }
                
                ((float*)self)[ q * result.step + w + result.start ] = alpha * sum;
            }        
    }
    
    return;
} // gemmFloatDirDirDir

void 
gemmFloatDirDirDirGang( void * self, clMat /*float*/A, clMat /*float*/B, float alpha, clMat /*float*/C, float beta, int localid, clMat /*float*/ result )
{
    assert(A.rows == C.rows);
    assert(A.cols == B.rows); 
    assert(B.cols == C.cols);
    assert(C.rows == result.rows);
    assert(C.cols == result.cols);

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
                float y = ((float*)self)[ q * A.step + e + A.start ] * ((float*)self)[ e * B.step + w + B.start ] - c;
                float t = sum + y;
                c = (t - sum) - y;
                sum = t;              
            }
                
            ((float*)self)[ q * result.step + w + result.start ] = alpha * sum  + beta * ((float*)self)[ q * C.step + w + C.start ];
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
                float y = ((float*)self)[ q * A.step + e + A.start ] * ((float*)self)[ e * B.step + w + B.start ] - c;
                float t = sum + y;
                c = (t - sum) - y;
                sum = t;              
            }
            
            ((float*)self)[ q * result.step + w + result.start ] = alpha * sum;
        }
    }
    
    return;
} // gemmFloatDirDirDirGang


void 
gemmFloatDirTransDirGang( void * self, clMat /*float*/A, clMat /*float*/B, float alpha, clMat /*float*/C, float beta, int localid, clMat /*float*/ result )
{
    assert(A.rows == C.rows);
    assert(A.cols == B.cols); 
    assert(B.rows == C.cols);
    assert(C.rows == result.rows);
    assert(C.cols == result.cols);

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
                float y = ((float*)self)[ q * A.step + e + A.start ] * ((float*)self)[ w * B.step + e + B.start ] - c;
                float t = sum + y;
                c = (t - sum) - y;
                sum = t;              
            }
                
            ((float*)self)[ q * result.step + w + result.start ] = alpha * sum  + beta * ((float*)self)[ q * C.step + w + C.start ];
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
                float y = ((float*)self)[ q * A.step + e + A.start ] * ((float*)self)[ w * B.step + e + B.start ] - c;
                float t = sum + y;
                c = (t - sum) - y;
                sum = t;              
            }
            
            ((float*)self)[ q * result.step + w + result.start ] = alpha * sum;
        }
    }
    
    return;
} // gemmFloatDirTransDirGang


void
expFloat( void * self, clMat /*float*/input, clMat /*float*/ output )
{
    assert(input.rows == output.rows);
    assert(input.cols == output.cols);

    int q, w;

    for ( q=0; q<input.rows; q++ )
        for ( w=0; w<input.cols; w++ ) {
            ((float*)self)[ q * output.step + w + output.start ] = 
        	exp(((float*)self)[ q * input.step + w + input.start ]);
        }

    return;
}

void
addFloat( void * self, clMat /*float*/input, float val, clMat /*float*/ output )
{
    assert(input.rows == output.rows);
    assert(input.cols == output.cols);

    int q, w;
    for ( q=0; q<input.rows; q++ )
        for ( w=0; w<input.cols; w++ )
            ((float*)self)[ q * output.step + w + output.start ] = 
        	val + ((float*)self)[ q * input.step + w + input.start ];

    return;
}

void 
divideFloat( void * self, float val, clMat /*float*/input, clMat /*float*/ output )
{
    assert(input.rows == output.rows);
    assert(input.cols == output.cols);

    int q, w;
    for ( q=0; q<input.rows; q++ )
        for ( w=0; w<input.cols; w++ )
            ((float*)self)[ q * output.step + w + output.start ] = 
        	val / ((float*)self)[ q * input.step + w + input.start ];

    return;
}

void
subtractFloat( void * self, clMat /*float*/input, float val, clMat /*float*/ output )
{
    assert(input.rows == output.rows);
    assert(input.cols == output.cols);

    int q, w;
    for ( q=0; q<input.rows; q++ )
        for ( w=0; w<input.cols; w++ )
            ((float*)self)[ q * output.step + w + output.start ] = 
        	((float*)self)[ q * input.step + w + input.start ] - val;
    
    return;
}

float
GetValueFloat( void * self, clMat /*float*/smat, int row, int col )
{
    assert(self);
    assert(row >= 0);
    assert(col >= 0);
    assert(row < smat.rows);
    assert(col < smat.cols);
    assert(smat.step >= smat.cols);
    
    return ((float*)self)[ row * smat.step + col + smat.start ];
    // return -1231.;
}

uint8_t
GetValueChar( void * self, clMat /* uint8_t*/ smat, int row, int col )
{
    assert(self);
    assert(row >= 0);
    assert(col >= 0);
    assert(row < smat.rows);
    assert(col < smat.cols);
    assert(smat.step >= smat.cols);
    
    return ((uint8_t*)self)[ row * smat.step + col + smat.start ];
    // return -1231.;
}

void
SetValueFloat( void * self, clMat /*float*/ smat, int row, int col, float value )
{
    assert(self);
    assert(row >= 0);
    assert(col >= 0);
    assert(row < smat.rows);
    assert(col < smat.cols);
    assert(smat.step >= smat.cols);
    assert(smat.step % smat.cols == 0);

    ((float*)self)[ row * smat.step + col + smat.start ] = value;
    return;    
}



float dotProductDirDir( void * self, clMat /*float*/A, clMat /*float*/B )
{
  assert( A.cols == 1 );
  assert( B.cols == 1 );
  assert( A.rows == B.rows );

  float result = 0.;
  float c = 0.;  

  int q;
  for (q=0; q<A.rows; q++) {
      float y = ((float*)self)[ q * A.step + A.start ] * ((float*)self)[ q * B.step + B.start ] - c;
      float t = result + y;
      c = (t - result) - y;
      result = t ;
  }
  
  return result;
}

float dotProductTransDir( void * self, clMat /*float*/A, clMat /*float*/B )
{
  assert( A.rows == 1 );
  assert( B.cols == 1 );
  assert( A.cols == B.rows );

  float result = 0.;
  float c = 0.;  

  int q;
  for (q=0; q<A.cols; q++) {
      float y = ((float*)self)[ q + A.start ] * ((float*)self)[ q * B.step + B.start ] - c;
      float t = result + y;
      c = (t - result) - y;
      result = t ;
  }
  
  return result;
}

void normalizeSample( void * self, clMat /*uint8_t*/  image, clMat /*float*/ * result )
{
  assert(result->cols == image.cols);
  assert(result->rows == image.rows);
    
  float sampleMean = meanChar(self, image);
  float sampleMin  = minChar(self, image);
  float sampleMax  = maxChar(self, image);

  sampleMax -= sampleMean;
  sampleMin -= sampleMean;

  sampleMax = fmax( fabs(sampleMin), fabs(sampleMax));

  if (sampleMax == 0.0) sampleMax = 1.0;

  convertFromCharToFloat( self, image, 1.0/sampleMax, -(1.0/sampleMax)*sampleMean, *result );

  *result = reshapeFloat( self, *result, image.rows * image.cols );
 
  return;
} // normalizeSample

void
generateResponseMap(
    void * self,
    const clMat /*uint8_t*/  image,
    const Point2i center,
    int mapSize, 
    int m_patchSize,
    // clMat /*float*/m_wIn,
    clMat /*float*/m_wOut,
    // clMat /*float*/m_U,
    clMat wIn,
    clMat bIn,
    
    // temporary variables
    clVector /*cMat float*/ patches,
    clVector /*cMat float*/ xOuts,
    clVector /*cMat float*/ es,
    // result
    clMat /*float*/ result
    )
{

  assert( result.rows == 2 * mapSize + 1 );
  assert( result.cols == 2 * mapSize + 1 );

//  clMat /*float*/wIn_A = GetBlockFloat( self, m_wIn, 0, m_wIn.rows, 0, m_wIn.cols - 1 );
  // clMat /*float*/wIn = CreateMatFloat( self, allocator, wIn_A.rows, m_U.rows );
  
//  assert(wIn.rows == wIn_A.rows);
//  assert(wIn.cols == m_U.rows);
  
  // clMat /*float*/bIn = GetBlockFloat( self, m_wIn, 0, m_wIn.rows, m_wIn.cols - 1, m_wIn.cols );
  clMat /*float*/wOut_tmp = GetBlockFloat( self, m_wOut, 0, m_wOut.rows, 0, m_wOut.cols - 1 );

  // {
  //     int localid = 0;
  //     for ( localid=0; localid<gangsize; localid++ )           
  //         gemmFloatDirTransDirGang( self, wIn_A, m_U, 1.0, wIn, 0.0, localid, wIn );
  // }

  {
      int localid = 0;
      for ( localid=0; localid<gangsize; localid++ ) {
          
          int work = 0;      
          int worksize = (2 * mapSize + 1) * (2 * mapSize + 1);
          
          float bOut = GetValueFloat( self, m_wOut, 0, m_wOut.cols - 1 );

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

              clMat /*uint8_t*/   imagePatch = GetBlockChar( self, image, cy - m_patchSize, cy + m_patchSize + 1, cx - m_patchSize, cx + m_patchSize + 1 );

              // printf("cy - m_patchSize = %d\n", cy - m_patchSize );
              // printf("cy + m_patchSize + 1 = %d\n", cy + m_patchSize + 1 );
              // printf("cx - m_patchSize = %d\n", cx - m_patchSize );
              // printf("cx + m_patchSize + 1 = %d\n", cx + m_patchSize + 1);
              
              // printMatChar( self, imagePatch, "imagePatch");
              
              // clMat /*float*/patch = CreateMatFloat( self, allocator, imagePatch.rows, imagePatch.cols );
              clMat patch = GetMatFromVector( self, patches, localid );
                            
              assert( patch.rows == imagePatch.rows );
              assert( patch.cols == imagePatch.cols );

              normalizeSample( self, imagePatch, &patch );
              // printMatFloat(self, patch, "patch");
                            
              // clMat /*float*/xOut = CreateMatFloat( self, allocator, bIn.rows, bIn.cols );
              clMat xOut = GetMatFromVector( self, xOuts, localid );

              assert( xOut.rows == bIn.rows );
              assert( xOut.cols == bIn.cols );
          
              gemmFloatDirDirDir( self, wIn, patch, -1.0, bIn, -1.0, xOut );

              // printMatFloat(self, xOut, "xOut");
              
              // clMat /*float*/e = CreateMatFloat( self, allocator, xOut.rows, xOut.cols);
              clMat e = GetMatFromVector( self, es, localid );
              assert( e.rows == xOut.rows );
              assert( e.cols == xOut.cols );              
          
              expFloat( self, xOut, e );

              // printMatFloat( self, e, "e");
              
              addFloat( self, e, 1.0, xOut );
              divideFloat( self, 2.0, xOut, e);
              addFloat( self, e, -1.0, xOut);

              // printMatFloat( self, xOut, "xOut" );
              
              SetValueFloat( self, result, ncy, ncx, 1./( 1. + exp(- dotProductTransDir( self, wOut_tmp, xOut) - bOut ) ) );

              // freeMatFloat( self, allocator, &e);
              // freeMatFloat( self, allocator, &xOut);
              // freeMatFloat( self, allocator, &patch);
          } // for localid 
      } // for gangsize
  } // localid block
  
//  assert(false);
  
  // freeMatFloat( self, allocator, &wIn);
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
    void * self,
    int memory_segments[], // representing the start of the gang's memory segment
    int m_visibleLandmarks_size,
    int m_mapSize,
    // clMat /*float*/shape,
    calcpackage packages[] // ,

    // // temporaries
    // clMat wIn,
    // clMat patches[],
    // clMat xOuts[],
    // clMat es[],

    // results
    //   clMat /*float*/responseMaps[]
    )
{
    assert(self);    
    assert(memory_segments);
    assert(packages);    
    
    // printf("calculateMaps started\n");
    int q;
    for (q=9; q<10 /*m_visibleLandmarks_size*/; q++ )
    {
        printf("***********************************************************\n");
        printf("***********************************************************\n");
        printf("************* iteration %d *********************************\n", q );
        printf("***********************************************************\n");
        printf("***********************************************************\n");
        
	// printf("processing patch %d/%d\n", q, m_visibleLandmarks_size );
        /* const int idx = m_visibleLandmarks[q]; */
        /* assert(idx==q); */

	int idx = q;

        Point2i center;

        float shape_x;
	float shape_y;
	shape_x = GetValueFloat( self + memory_segments[q], packages[q].input.shape, 2*idx, 0 );
	shape_y = GetValueFloat( self + memory_segments[q], packages[q].input.shape, 2*idx+1, 0 );

        center.x = cvRound(shape_x);
        center.y = cvRound(shape_y);

        // clMat wIn = CreateMatFloat( self, allocator, m_wIns[idx].rows, m_Us[idx].rows );

        // clMat patches[gangsize];
        // clMat xOuts[gangsize];
        // clMat es[gangsize];
        
        // {
        //   int w;
        //   for (w=0; w<gangsize; w++) {
        //     patches[w] = CreateMatFloat( self, allocator, 2 * m_patchSizes[idx] + 1, 2 * m_patchSizes[idx] + 1 );
        //     xOuts[w] = CreateMatFloat( self, allocator, m_wIns[idx].rows, 1 );
        //     es[w] = CreateMatFloat( self, allocator, xOuts[w].rows, xOuts[w].cols );
        //   }
        // }
        
	// responseMaps[q] = m_classifiers[idx].generateResponseMap( alignedImage, center, m_mapSize );
        
	generateResponseMap(
            self + memory_segments[idx],
            packages[idx].input.alignedImage,
            center,
            m_mapSize,
            packages[idx].input.m_patchSize,
//            packages[idx].input.m_wIn,
            packages[idx].input.m_wOut,
//            packages[idx].input.m_U,
            packages[idx].input.wIn,
            packages[idx].input.bIn,

            // temporary
            packages[idx].tmp.patches,
            packages[idx].tmp.xOuts,
            packages[idx].tmp.es,
            // result
            packages[idx].output.responseMap );

        printMatFloat( self + memory_segments[idx], packages[idx].output.responseMap, "packages[idx].output.responseMap");
        

        // freeMatFloat(self, allocator, &wIn);
        // {
        //   int w;
        //   for (w=0; w<gangsize; w++) {
        //     freeMatFloat( self, allocator, &(patches[w]) );
        //     freeMatFloat( self, allocator, &(xOuts[w]) );
        //     freeMatFloat( self, allocator, &(es[w]) );
        //   }
        // }
        
    } // for q in visiblelandmarks_size

    // printf("calculateMaps finished\n");
    return;
} // calculateMaps


// LuM end of file
