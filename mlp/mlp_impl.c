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

// returns alpha*A*B + beta * C
void 
gemmFloatDirDirDir(
    void * self,
    clMat /*float*/A,
    clMat /*uint8_t*/B,
    float alpha,
    clMat /*float*/C,
    float beta,
    normalization norm,
    clMat /*float*/ result
    )
{
    assert(A.rows == C.rows);
    assert(A.cols == B.rows * B.cols ); 
    assert(C.rows == result.rows);
    assert(C.cols == result.cols);

    int q, w, e;
    float sum=0.0f;
    float c=0.0f;

    if ( fabs(beta) > 0.000001 ) {
        for ( q=0; q<C.rows; q++ )
            for ( w=0; w<C.cols; w++ )
            {
                sum = 0;
                for ( e=0; e<A.cols; e++ )
                {
                    float a = ((float*)self)[ q * A.step + e + A.start ];
                    int b_row = e / B.rows;                    
                    int b_col = e % B.cols;                    
                    float b = ((uint8_t*)self)[ b_row * B.step + b_col + B.start ];
                    float y = a * (norm.shift + norm.stride * b) - c;
                    float t = sum + y;
                    c = (t - sum) - y;
                    sum = t;
                }
                
                ((float*)self)[ q * result.step + w + result.start ] = alpha * sum  + beta * ((float*)self)[ q * C.step + w + C.start ];
            }
    }
    else // NOT fabs(beta) > 0.000001
    {
        assert(false);        
    }

//    assert(false);
    
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
    float sum=0.0f;
    float c=0.0f;    
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
    float sum=0.0f;
    float c=0.0f;    
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
activateOutput( void * self, clMat /*float*/input )
{
    int q, w;

    for ( q=0; q<input.rows; q++ )
        for ( w=0; w<input.cols; w++ ) {
            float val = ((float*)self)[ q * input.step + w + input.start ];
            float result = 2.0f/(exp(val) + 1) - 1;
            ((float*)self)[ q * input.step + w + input.start ] = result;
        }

    return;
} // activateOutput

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

  float result = 0.0f;
  float c = 0.0f;  

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

normalization
normalizeSample( void * self, clMat /*uint8_t*/  image ) // , clMat /*float*/ * result )
{
  // assert(result->cols == image.cols);
  // assert(result->rows == image.rows);

    int q,w;
    float sum=0.0f;
    float c = 0.0f; // kahan summation
    uint8_t minvalue = 255;
    uint8_t maxvalue = 0;
    
    for ( q=0; q<image.rows; q++ )
        for ( w=0; w<image.cols; w++ )
        {
	    uint8_t pixel = ((uint8_t*)self)[ q * image.step + w + image.start ];
	    minvalue = fmin( minvalue, pixel );
	    maxvalue = fmax( maxvalue, pixel );
            float y = pixel - c;
            float t = sum + y;
            c = (t - sum) - y;        
            sum = t;        
        }
    
    
    normalization nresult;
    
    float sampleMean = sum / (image.rows * image.cols);
    float sampleMin  = minvalue; 
    float sampleMax  = maxvalue;
    
    sampleMax -= sampleMean;
    sampleMin -= sampleMean;
    
    sampleMax = fmax( fabs(sampleMin), fabs(sampleMax));
    
    if (sampleMax == 0.0) sampleMax = 1.0;
    
    nresult.shift  = -(1.0/sampleMax)*sampleMean;
    nresult.stride = 1.0/sampleMax;
    
    return nresult;
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
//    clVector /*cMat float*/ patches,
//    clVector /*cMat float*/ xOuts,
//    clVector /*cMat float*/ es,
    // result
    clMat /*float*/ result
    )
{

  assert( result.rows == 2 * mapSize + 1 );
  assert( result.cols == 2 * mapSize + 1 );

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

              // clMat patch = GetMatFromVector( self, patches, localid );
                            
              // assert( patch.rows == imagePatch.rows );
              // assert( patch.cols == imagePatch.cols );

              normalization norm = normalizeSample( self, imagePatch );

              // clMat xOut = GetMatFromVector( self, xOuts, localid );

              // assert( xOut.rows == bIn.rows );
              // assert( xOut.cols == bIn.cols );
          
              // gemmFloatDirDirDir( self, wIn, imagePatch /*patch*/, -1.0, bIn, -1.0, norm, xOut );

              // activateOutput( self, xOut );
              
              // SetValueFloat( self, result, ncy, ncx, 1./( 1. + exp(- dotProductTransDir( self, wOut_tmp, xOut) - bOut ) ) );

          } // for localid 
      } // for gangsize
  } // localid block
  
//  assert(false);
  
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
    char * self,
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
    for (q=0; q<m_visibleLandmarks_size; q++ )
    {
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
//            packages[idx].tmp.patches,
            // packages[idx].tmp.xOuts,
//            packages[idx].tmp.es,
            // result
            packages[idx].output.responseMap );
        
    } // for q in visiblelandmarks_size

    // printf("calculateMaps finished\n");
    return;
} // calculateMaps


// LuM end of file
