// -*- c++ -*-
// UjoImro, 2013
// Experimental code for the CARP Project
// Copyright (c) RealEyes, 2013
// This is a c-implementation of the PCA->MLP response map calculation

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
//    clVector /* clMat[] */ patches;
    clVector /* clMat[] */ xOuts;
//    clVector /* clMat[] */ es;        
} calctemp; // struct
    
typedef struct {
    // results
    clMat /*float*/ responseMap;        
} calcoutput;
            
typedef struct {
    calcinput  input;
    calctemp   tmp;
    calcoutput output;
} calcpackage;

typedef struct {
    float shift;
    float stride;        
} normalization;
    

__constant const int MAX_INT = ~(1 << (8*sizeof(int) - 1));
// __constant const int false = (1!=1);
__constant const int gangsize = 32;


clMat GetMatFromVector( __global void * self, clVector /*clMat*/ vec, int idx )
{
//    assert(self);
//    assert(idx>=0);
//    assert(idx<vec.size);
    
    return ((__global clMat*)self)[ vec.start + vec.step * idx ];
} // GetMatFromVector



void SetMatToVector( __global void * self, clVector /*clMat*/ vec, int idx, clMat val )
{
//    assert(self);
//    assert(idx>=0);
//    assert(idx<vec.size);
    
    ((__global clMat*)self)[ vec.start + vec.step * idx ] = val;
    return;
} // GetMatFromVector

clMat /*uint8_t*/ 
GetBlockChar( __global void * self, clMat /*uint8_t*/  smat, int row_from, int row_to, int col_from, int col_to )
{    
    // assert(row_from>=0);
    // assert(col_from>=0);
    // assert(row_from<row_to);
    // assert(col_from<col_to);
    // assert(row_to<=smat.rows);
    // assert(col_to<=smat.cols);
    
    clMat /*uchar*/  result;
    result.rows  = row_to - row_from;
    result.cols  = col_to - col_from;
    result.step  = smat.step;
    result.start = smat.start + row_from * smat.step + col_from;

    return result;
}

clMat
GetBlockFloat( __global void * self, clMat /*float*/smat, int row_from, int row_to, int col_from, int col_to )
{  
    // assert(row_from>=0);
    // assert(col_from>=0);
    // assert(row_from<row_to);
    // assert(col_from<col_to);
    // assert(row_to<=smat.rows);
    // assert(col_to<=smat.cols);

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
    __global void * self,
    clMat /*float*/A,
    clMat /*uint8_t*/B,
    float alpha,
    clMat /*float*/C,
    float beta,
    normalization norm,
    clMat /*float*/ result
    )
{
    // assert(A.rows == C.rows);
    // assert(A.cols == B.rows); 
    // assert(B.cols == C.cols);
    // assert(C.rows == result.rows);
    // assert(C.cols == result.cols);

    int q, w, e;
    float sum=0.0f;
    float c=0.0f;

    for ( q=0; q<C.rows; q++ )
        for ( w=0; w<C.cols; w++ )
        {	    
            sum = 0;
            for ( e=0; e<A.cols; e++ )
            {              
                float a = ((__global float*)self)[ q * A.step + e + A.start ];
                int b_row = e / B.rows;                    
                int b_col = e % B.cols;                    
                float b = ((__global uchar*)self)[ b_row * B.step + b_col + B.start ];
                float y = a * (norm.shift + norm.stride * b) - c;
                float t = sum + y;
                c = (t - sum) - y;
                sum = t;              
            }
                
            ((__global float*)self)[ q * result.step + w + result.start ] = alpha * sum  + beta * ((__global float*)self)[ q * C.step + w + C.start ];
        }
    
    
    return;
} // gemmFloatDirDirDir

void
activateOutput( __global void * self, clMat /*float*/input )
{
    int q, w;

    for ( q=0; q<input.rows; q++ )
        for ( w=0; w<input.cols; w++ ) {
            float val = ((__global float*)self)[ q * input.step + w + input.start ];
            float result = 2.0f/(exp(val) + 1) - 1;
            ((__global float*)self)[ q * input.step + w + input.start ] = result;
        }

    return;
} // activateOutput


float
GetValueFloat( __global void * self, clMat /*float*/smat, int row, int col )
{
    // assert(self);
    // assert(row >= 0);
    // assert(col >= 0);
    // assert(row < smat.rows);
    // assert(col < smat.cols);
    // assert(smat.step >= smat.cols);
    // assert(smat.step % smat.cols == 0);

    return ((__global float*)self)[ row * smat.step + col + smat.start ];
    // return -1231.;
}

uchar
GetValueChar( __global void * self, clMat /* uint8_t*/ smat, int row, int col )
{
    // assert(self);
    // assert(row >= 0);
    // assert(col >= 0);
    // assert(row < smat.rows);
    // assert(col < smat.cols);
    // assert(smat.step >= smat.cols);
    
    return ((__global uchar*)self)[ row * smat.step + col + smat.start ];
    // return -1231.;
}

void
SetValueFloat( __global void * self, clMat /*float*/ smat, int row, int col, float value )
{
    // assert(self);
    // assert(row >= 0);
    // assert(col >= 0);
    // assert(row < smat.rows);
    // assert(col < smat.cols);
    // assert(smat.step >= smat.cols);
    // assert(smat.step % smat.cols == 0);

    ((__global float*)self)[ row * smat.step + col + smat.start ] = value;
    return;    
}

void
SetValueChar( __global void * self, clMat /*float*/ smat, int row, int col, uchar value )
{
    // assert(self);
    // assert(row >= 0);
    // assert(col >= 0);
    // assert(row < smat.rows);
    // assert(col < smat.cols);
    // assert(smat.step >= smat.cols);
    // assert(smat.step % smat.cols == 0);

    ((__global uchar*)self)[ row * smat.step + col + smat.start ] = value;
    return;    
}



float dotProductDirDir( __global void * self, clMat /*float*/A, clMat /*float*/B )
{
  // assert( A.cols == 1 );
  // assert( B.cols == 1 );
  // assert( A.rows == B.rows );

  float result = 0.0f;
  float c = 0.0f;  

  int q;
  for (q=0; q<A.rows; q++) {
      float y = ((__global float*)self)[ q * A.step + A.start ] * ((__global float*)self)[ q * B.step + B.start ] - c;
      float t = result + y;
      c = (t - result) - y;
      result = t ;
  }
  
  return result;
}

float dotProductTransDir( __global void * self, clMat /*float*/A, clMat /*float*/B )
{
  // assert( A.rows == 1 );
  // assert( B.cols == 1 );
  // assert( A.cols == B.rows );

  float result = 0.0f;
  float c = 0.0f;  

  int q;
  for (q=0; q<A.cols; q++) {
      float y = ((__global float*)self)[ q + A.start ] * ((__global float*)self)[ q * B.step + B.start ] - c;
      float t = result + y;
      c = (t - result) - y;
      result = t ;
  }
  
  return result;
}

normalization
normalizeSample( __global void * self, clMat /*uint8_t*/  image ) // , clMat /*float*/ * result )
{
  // assert(result->cols == image.cols);
  // assert(result->rows == image.rows);

    int q,w;
    float sum=0.0f;
    float c = 0.0f; // kahan summation
    uchar minvalue = 255;
    uchar maxvalue = 0;
    
    for ( q=0; q<image.rows; q++ )
        for ( w=0; w<image.cols; w++ )
        {
	    uchar pixel = ((__global uchar*)self)[ q * image.step + w + image.start ];
	    minvalue = min( minvalue, pixel );
	    maxvalue = max( maxvalue, pixel );
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
    __global void * self,
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
    clVector /*cMat float*/ xOuts,
//    clVector /*cMat float*/ es,
    // result
    clMat /*float*/ result
    )
{
    
  // assert( result.rows == 2 * mapSize + 1 );
  // assert( result.cols == 2 * mapSize + 1 );

  //  assert(wIn.rows == wIn_A.rows);
  //  assert(wIn.cols == m_U.rows);

    clMat /*float*/wOut_tmp = GetBlockFloat( self, m_wOut, 0, m_wOut.rows, 0, m_wOut.cols - 1 );

    {
        int localid = get_local_id(0);
//        for ( localid=0; localid<gangsize; localid++ ) {
          
            int work = 0;      
            int worksize = (2 * mapSize + 1) * (2 * mapSize + 1);
          
            float bOut = GetValueFloat( self, m_wOut, 0, m_wOut.cols - 1 );
            
            for ( work=0; work < (worksize / gangsize) + 1; work++ )
            {
              int index = work * gangsize + localid;
              if (index>=worksize) continue;
              int ncy = index / (2 * mapSize + 1);
              int cy  = center.y - mapSize + ncy;
              int ncx = index % (2 * mapSize + 1);
              int cx  = center.x - mapSize + ncx;

              clMat /*uint8_t*/   imagePatch = GetBlockChar( self, image, cy - m_patchSize, cy + m_patchSize + 1, cx - m_patchSize, cx + m_patchSize + 1 );

              normalization norm = normalizeSample( self, imagePatch );

              clMat xOut = GetMatFromVector( self, xOuts, localid );

              gemmFloatDirDirDir( self, wIn, imagePatch, -1.0f, bIn, -1.0f, norm, xOut );

              activateOutput( self, xOut );              
              
              SetValueFloat( self, result, ncy, ncx, 1.0f/( 1.0f + exp(- dotProductTransDir( self, wOut_tmp, xOut) - bOut ) ) );

            } // for gangsize
            //       } // for localid
    } // localid block
  
//  // assert(false);
  
  // freeMatFloat( self, allocator, &wIn);
  // end of classic impl
    return;
} // generateResponseMap

int 
cvRound( float value )
{
    return (int)(value + (value >= 0 ? 0.5f : -0.5f));
} // cvRound


__kernel void
calculateMaps(
    // global variables
    __global void * self,
    __constant int memory_segments[], // representing the start of the gang's memory segment
    int m_visibleLandmarks_size,
    int m_mapSize,

    __constant calcpackage packages[],

    // local buffer
    __local void * buffer
    )
{

    // int q = get_group_id(0);
    // if (q >= m_visibleLandmarks_size) return;
    // int localid = get_local_id(0);
    // if (localid > 0) return;    

    int q = get_group_id(0);
//    if (q > 0) return;
    int localid = get_local_id(0);
//    if (localid > 0) return;    

//    for (q=0; q<m_visibleLandmarks_size; q++) {
            
        int idx = q;

        Point2i center;
        
        float shape_x;
        float shape_y;
        calcpackage package = packages[q];
        shape_x = GetValueFloat( self + memory_segments[q], package.input.shape, 2*idx, 0 );
        shape_y = GetValueFloat( self + memory_segments[q], package.input.shape, 2*idx+1, 0 );
        
        center.x = cvRound(shape_x);
        center.y = cvRound(shape_y);

        generateResponseMap(
            self + memory_segments[idx],
            package.input.alignedImage,
            center,
            m_mapSize,
            package.input.m_patchSize,
            package.input.m_wOut,
            package.input.wIn,
            package.input.bIn,
            
            // temporary
//            package.tmp.patches,
            package.tmp.xOuts,
//            package.tmp.es,
            // result
            package.output.responseMap );

//    } // for q in m_visibleLandmarks
    
    return;
} // calculateMaps


// LuM end of file
