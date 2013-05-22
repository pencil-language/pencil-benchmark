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
    clVector /* clMat[] */ patches;
    clVector /* clMat[] */ xOuts;
    clVector /* clMat[] */ es;        
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


__kernel
void
copyToFloat( __global void * self, clMat /*float*/input, clMat /*float*/ output )
{
    // assert(input.rows == output.rows);
    // assert(input.cols == output.cols);
    
    int q, w;
    for ( q=0; q<input.rows; q++ )
        for ( w=0; w<input.cols; w++ )
            ((__global float*)self)[ q * output.step + w + output.start ] =
                ((__global float*)self)[ q * input.step + w + input.start ];
    
    return;
} // copyTo

__kernel
void
transposeFloat( __global void * self, clMat /*float*/input, clMat /*float*/ output )
{
    int q,w;
    // assert(output.rows == input.cols);
    // assert(output.cols == input.rows);
    for (q=0; q<input.rows; q++)
        for (w=0; w<input.cols; w++)
            ((__global float*)self)[ w * output.step + q + output.start ]
        	= ((__global float*)self)[ q * input.step + w + input.start ];

    return;
} // transposeFloat

__kernel
void
transposeFloatGang( __global void * self, clMat /*float*/input, int localid, clMat /*float*/ output )
{
    // assert(output.rows == input.cols);
    // assert(output.cols == input.rows);
    
    int worksize = input.cols * input.rows;
    int work = 0;
            
    for ( work=0; work < (worksize/gangsize) + 1; work++) {
        int index = work * gangsize + localid;
        if (index>=worksize) continue;
        int q = index / input.cols;
        int w = index % input.cols;
              
        ((__global float*)self)[ w * output.step + q + output.start ] 
            = ((__global float*)self)[ q * input.step + w + input.start ];
    }
} // transposeFloatGang


float
meanChar( __global void * self, clMat /*uint8_t*/  input )
{    
    int q,w;
    float sum=0;
    float c = 0; // kahan summation
    
    for ( q=0; q<input.rows; q++ )
      for ( w=0; w<input.cols; w++ )
      {
        float y = ((__global uchar*)self)[ q * input.step + w + input.start ] - c;
        float t = sum + y;
        c = (t - sum) - y;        
        sum = t;        
      }
    
    return sum / ( input.rows * input.cols );
} // meanFloat

uchar
minChar( __global void * self, clMat /*uint8_t*/  input )
{    
    int q,w;
    uchar minvalue = 255;

    for ( q=0; q<input.rows; q++ )
        for ( w=0; w<input.cols; w++ )
            minvalue = min( minvalue, ((__global uchar*)self)[ q * input.step + w + input.start ] );
    
    return minvalue;
} // minFloat

uchar
maxChar( __global void * self, clMat /*uint8_t*/  input )
{
    int q,w;
    uchar maxvalue = 0;

    for ( q=0; q<input.rows; q++ )
        for ( w=0; w<input.cols; w++ )
            maxvalue = max( maxvalue, ((__global uchar*)self)[ q * input.step + w + input.start ] );
    
    return maxvalue;
} // maxFloat


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

__kernel
void
convertFromCharToFloat( __global void * self, clMat /*uint8_t*/  from, float quotient, float shift, clMat /*float*/ to )
{
    // assert(from.rows == to.rows);
    // assert(from.cols == to.cols);
    
    int q, w;
    for ( q=0; q<from.rows; q++ )
        for ( w=0; w<from.cols; w++ )
            ((__global float*)self)[ q * to.step + w + to.start ] = 
        	quotient * ((__global uchar*)self)[ q * from.step + w + from.start ] + shift;
    return;
}

clMat
reshapeFloat( __global void * self, clMat /*float*/smat, int new_rows )
{

    // assert(smat.cols == smat.step);
    // assert( (smat.cols * smat.rows) % new_rows == 0 );

    clMat /*float*/result;
    result.rows = new_rows;
    result.cols = smat.cols * smat.rows / new_rows;
    result.step = result.cols;
    result.start = smat.start;

    return result;
} // reshapeFloat

// returns alpha*A*B + beta * C
__kernel
void 
gemmFloatDirDirDir( __global void * self, clMat /*float*/A, clMat /*float*/B, float alpha, clMat /*float*/C, float beta, clMat /*float*/ result )
{
    // assert(A.rows == C.rows);
    // assert(A.cols == B.rows); 
    // assert(B.cols == C.cols);
    // assert(C.rows == result.rows);
    // assert(C.cols == result.cols);

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
                    float y = ((__global float*)self)[ q * A.step + e + A.start ] * ((__global float*)self)[ e * B.step + w + B.start ] - c;
                    float t = sum + y;
                    c = (t - sum) - y;
                    sum = t;              
                }
                
                ((__global float*)self)[ q * result.step + w + result.start ] = alpha * sum  + beta * ((__global float*)self)[ q * C.step + w + C.start ];
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
                    float y = ((__global float*)self)[ q * A.step + e + A.start ] * ((__global float*)self)[ e * B.step + w + B.start ] - c;
                    float t = sum + y;
                    c = (t - sum) - y;
                    sum = t;              
                }
                
                ((__global float*)self)[ q * result.step + w + result.start ] = alpha * sum;
            }        
    }
    
    return;
} // gemmFloatDirDirDir

__kernel
void 
gemmFloatDirDirDirGang( __global void * self, clMat /*float*/A, clMat /*float*/B, float alpha, clMat /*float*/C, float beta, int localid, clMat /*float*/ result )
{
    // assert(A.rows == C.rows);
    // assert(A.cols == B.rows); 
    // assert(B.cols == C.cols);
    // assert(C.rows == result.rows);
    // assert(C.cols == result.cols);

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
                float y = ((__global float*)self)[ q * A.step + e + A.start ] * ((__global float*)self)[ e * B.step + w + B.start ] - c;
                float t = sum + y;
                c = (t - sum) - y;
                sum = t;              
            }
                
            ((__global float*)self)[ q * result.step + w + result.start ] = alpha * sum  + beta * ((__global float*)self)[ q * C.step + w + C.start ];
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
                float y = ((__global float*)self)[ q * A.step + e + A.start ] * ((__global float*)self)[ e * B.step + w + B.start ] - c;
                float t = sum + y;
                c = (t - sum) - y;
                sum = t;              
            }
            
            ((__global float*)self)[ q * result.step + w + result.start ] = alpha * sum;
        }
    }
    
    return;
} // gemmFloatDirDirDirGang


__kernel
void 
gemmFloatDirTransDirGang( __global void * self, clMat /*float*/A, clMat /*float*/B, float alpha, clMat /*float*/C, float beta, int localid, clMat /*float*/ result )
{
    // assert(A.rows == C.rows);
    // assert(A.cols == B.cols); 
    // assert(B.rows == C.cols);
    // assert(C.rows == result.rows);
    // assert(C.cols == result.cols);

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
                float y = ((__global float*)self)[ q * A.step + e + A.start ] * ((__global float*)self)[ w * B.step + e + B.start ] - c;
                float t = sum + y;
                c = (t - sum) - y;
                sum = t;              
            }
                
            ((__global float*)self)[ q * result.step + w + result.start ] = alpha * sum  + beta * ((__global float*)self)[ q * C.step + w + C.start ];
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
                float y = ((__global float*)self)[ q * A.step + e + A.start ] * ((__global float*)self)[ w * B.step + e + B.start ] - c;
                float t = sum + y;
                c = (t - sum) - y;
                sum = t;              
            }
            
            ((__global float*)self)[ q * result.step + w + result.start ] = alpha * sum;
        }
    }
    
    return;
} // gemmFloatDirTransDirGang


__kernel
void
expFloat( __global void * self, clMat /*float*/input, clMat /*float*/ output )
{
    // assert(input.rows == output.rows);
    // assert(input.cols == output.cols);
    
    int q, w;
    for ( q=0; q<input.rows; q++ )
        for ( w=0; w<input.cols; w++ )
            ((__global float*)self)[ q * output.step + w + output.start ] =     
                exp(((__global float*)self)[ q * input.step + w + input.start ]);

    return;
}

__kernel
void
expVecFloat( __global void * self, clMat sample, clVector outputs )
{

    for (int q=0; q<33; q++ ) {       
        clMat buf = GetMatFromVector( self, outputs, q );

        expFloat( self, sample, buf );        
    }
        
    return;    
}


__kernel
void
addFloat( __global void * self, clMat /*float*/input, float val, clMat /*float*/ output )
{
    // assert(input.rows == output.rows);
    // assert(input.cols == output.cols);

    int q, w;
    for ( q=0; q<input.rows; q++ )
        for ( w=0; w<input.cols; w++ )
            ((__global float*)self)[ q * output.step + w + output.start ] = 
        	val + ((__global float*)self)[ q * input.step + w + input.start ];

    return;
}

__kernel
void 
divideFloat( __global void * self, float val, clMat /*float*/input, clMat /*float*/ output )
{
    // assert(input.rows == output.rows);
    // assert(input.cols == output.cols);

    int q, w;
    for ( q=0; q<input.rows; q++ )
        for ( w=0; w<input.cols; w++ )
            ((__global float*)self)[ q * output.step + w + output.start ] = 
        	val / ((__global float*)self)[ q * input.step + w + input.start ];

    return;
}

__kernel
void
subtractFloat( __global void * self, clMat /*float*/input, float val, clMat /*float*/ output )
{
    // assert(input.rows == output.rows);
    // assert(input.cols == output.cols);

    int q, w;
    for ( q=0; q<input.rows; q++ )
        for ( w=0; w<input.cols; w++ )
            ((__global float*)self)[ q * output.step + w + output.start ] = 
        	((__global float*)self)[ q * input.step + w + input.start ] - val;
    
    return;
}

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

    ((__global float*)self)[0]=11.;
    


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
    
    return ((uchar*)self)[ row * smat.step + col + smat.start ];
    // return -1231.;
}

__kernel
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



float dotProductDirDir( __global void * self, clMat /*float*/A, clMat /*float*/B )
{
  // assert( A.cols == 1 );
  // assert( B.cols == 1 );
  // assert( A.rows == B.rows );

  float result = 0.;
  float c = 0.;  

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

  float result = 0.;
  float c = 0.;  

  int q;
  for (q=0; q<A.cols; q++) {
      float y = ((__global float*)self)[ q + A.start ] * ((__global float*)self)[ q * B.step + B.start ] - c;
      float t = result + y;
      c = (t - result) - y;
      result = t ;
  }
  
  return result;
}

// __kernel
// void 
// normalizeSample( __global void * self, clMat /*uint8_t*/  image, clMat /*float*/ * result )
// {
//   // assert(result->cols == image.cols);
//   // assert(result->rows == image.rows);
    
//   float sampleMean = meanChar(self, image);
//   float sampleMin  = minChar(self, image);
//   float sampleMax  = maxChar(self, image);

//   sampleMax -= sampleMean;
//   sampleMin -= sampleMean;

//   sampleMax = fmax( fabs(sampleMin), fabs(sampleMax));

//   if (sampleMax == 0.0) sampleMax = 1.0;

//   convertFromCharToFloat( self, image, 1.0/sampleMax, -(1.0/sampleMax)*sampleMean, *result );

//   *result = reshapeFloat( self, *result, image.rows * image.cols );
 
//   return;
// } // normalizeSample

__kernel
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
    clVector /*cMat float*/ patches,
    clVector /*cMat float*/ xOuts,
    clVector /*cMat float*/ es,
    // result
    clMat /*float*/ result
    )
{

  // assert( result.rows == 2 * mapSize + 1 );
  // assert( result.cols == 2 * mapSize + 1 );

  //  assert(wIn.rows == wIn_A.rows);
  //  assert(wIn.cols == m_U.rows);
    ((__global float*)self)[0]=11.;
    


    clMat /*float*/wOut_tmp = GetBlockFloat( self, m_wOut, 0, m_wOut.rows, 0, m_wOut.cols - 1 );

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
            ((__global float*)self)[0]=11.;
            
            for ( work=0; work < (worksize / gangsize) + 1; work++ )
            {
              int index = work * gangsize + localid;
              if (index>=worksize) continue;
              int ncy = index / (2 * mapSize + 1);
              int cy  = center.y - mapSize + ncy;
              int ncx = index % (2 * mapSize + 1);
              int cx  = center.x - mapSize + ncx;

              clMat /*uint8_t*/   imagePatch = GetBlockChar( self, image, cy - m_patchSize, cy + m_patchSize + 1, cx - m_patchSize, cx + m_patchSize + 1 );

              clMat patch = GetMatFromVector( self, patches, localid );
                            
              // normalizeSample( self, imagePatch, &patch );
          
              clMat xOut = GetMatFromVector( self, xOuts, localid );

              // gemmFloatDirDirDir( self, wIn, patch, -1.0, bIn, -1.0, xOut );
              
              clMat e = GetMatFromVector( self, es, localid );
              
              ((__global float*)self)[0]=11.;
              expFloat( self, xOut, e );
          
              // addFloat( self, e, 1.0, xOut );
              // divideFloat( self, 2.0, xOut, e);
              // addFloat( self, e, -1.0, xOut);
              
              // SetValueFloat( self, result, ncy, ncx, 1./( 1. + exp(- dotProductTransDir( self, wOut_tmp, xOut) - bOut ) ) );
          
            } // for gangsize
        } // for localid
    } // localid block
  
//  // assert(false);
  
  // freeMatFloat( self, allocator, &wIn);
  // end of classic impl
    return;
} // generateResponseMap

int 
cvRound( float value )
{
    return (int)(value + (value >= 0 ? 0.5 : -0.5));
} // cvRound


__kernel 
void
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
    if (get_global_id(0)>0) return;

    ((__global float*)self)[0]=11.;
    
    int q;
    for (q=0; q<1 /*m_visibleLandmarks_size*/; q++ )
    {
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
            packages[idx].input.m_wOut,
            packages[idx].input.wIn,
            packages[idx].input.bIn,

            // temporary
            packages[idx].tmp.patches,
            packages[idx].tmp.xOuts,
            packages[idx].tmp.es,
            // result
            packages[idx].output.responseMap );

    } // for q in visiblelandmarks_size

    return;
} // calculateMaps


// LuM end of file
