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
    clMat /* float */ m_wOut;
    clMat /* float */ wIn;
    clMat /* float */ bIn;        
} calcinput; // struct 

typedef struct {
    int x;
    int y;
} Point2i; // struct Point2i

typedef struct {
    // results
    clMat /*float*/ responseMap;        
} calcoutput;
            
typedef struct {
    calcinput  input;
    calcoutput output;
} calcpackage;

typedef struct {
    float shift;
    float stride;
} normalization;
    

__constant const int MAX_INT = ~(1 << (8*sizeof(int) - 1));
// __constant const int false = (1!=1);
__constant const bool precise = false;

clMat /*uint8_t*/ 
GetBlockChar( __local void * self, clMat /*uint8_t*/  smat, int row_from, int row_to, int col_from, int col_to )
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
GetBlockFloat( __local void * self, clMat /*float*/smat, int row_from, int row_to, int col_from, int col_to )
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
float
gemmFloatDirDirDir(
    __local void * self,
    clMat /*float*/A,
    clMat /*uint8_t*/B,
    float alpha,
    clMat /*float*/C,
    float beta,
    normalization norm,
    clMat /* float */ dotmatrix//,
//    clMat /*float*/ result
    )
{
    // assert(A.rows == C.rows);
    // assert(A.cols == B.rows); 
    // assert(B.cols == C.cols);
    // assert(C.rows == result.rows);
    // assert(C.cols == result.cols);
    // assert(C.cols == 1 );
    
    float dot = 0.0f;
    
    for ( int q=0; q<C.rows; q++ ) {
        float sum = 0.0f;

        for ( int e=0; e < A.cols; e++ )
        {              
            float a = ((__local float*)self)[q * A.step + e + A.start];
            int b_row = e / B.rows;
            int b_col = e % B.cols;
            float b = ((__local uchar*)self)[ b_row * B.step + b_col + B.start ];
            sum += a * (norm.shift + norm.stride * b);
        }
        
        float val    = alpha * sum  + beta * ((__local float*)self)[ q * C.step /* + w==0*/ + C.start ];
        float resval = 2.0f/(exp(val) + 1) - 1;
        
        dot += resval * ((__local float*)self)[ q + dotmatrix.start ];     
    } // for q 
    
    return dot;
} // gemmFloatDirDirDir

float
GetValueFloat( __local void * self, clMat /*float*/smat, int row, int col )
{
    // assert(self);
    // assert(row >= 0);
    // assert(col >= 0);
    // assert(row < smat.rows);
    // assert(col < smat.cols);
    // assert(smat.step >= smat.cols);
    // assert(smat.step % smat.cols == 0);

    return ((__local float*)self)[ row * smat.step + col + smat.start ];
    // return -1231.;
}

void
SetValueFloat( __local void * self, clMat /*float*/ smat, int row, int col, float value )
{
    // assert(self);
    // assert(row >= 0);
    // assert(col >= 0);
    // assert(row < smat.rows);
    // assert(col < smat.cols);
    // assert(smat.step >= smat.cols);
    // assert(smat.step % smat.cols == 0);

    ((__local float*)self)[ row * smat.step + col + smat.start ] = value;
    return;    
}

normalization
normalizeSample( __local void * self, clMat /*uint8_t*/  image ) // , clMat /*float*/ * result )
{
  // assert(result->cols == image.cols);
  // assert(result->rows == image.rows);

    float4 sum4      = (float4)(0.0f);
    uchar4 minvalue4 = (uchar4)(255);
    uchar4 maxvalue4 = (uchar4)(0);
    float sum      = 0.0f;
    uchar minvalue = 255;
    uchar maxvalue = 0;

    for ( int q=0; q<image.rows; q++ ) {
        int start_index = image.start + q * image.step;
        int end_index = image.start + q * image.step + image.cols;

        int start_aligned = start_index - start_index % 4 + 4;
        int end_aligned = end_index - end_index % 4;
        
        for ( int index = start_index; index < start_aligned; index++ )
        {
            uchar pixel = ((__local uchar*)self)[index];
            minvalue = min( minvalue, pixel );
            maxvalue = max( maxvalue, pixel );
            sum += pixel;            
        }

        for ( int index = start_aligned / 4; index < end_aligned / 4; index ++ )
        {
            uchar4 pixel4 = ((__local uchar4*)self)[index];
            minvalue4 = min( minvalue4, pixel4 );
            maxvalue4 = max( maxvalue4, pixel4 );
            sum4 += convert_float4(pixel4);
        }
                
        for ( int index = end_aligned; index < end_index; index++ )
        {
            uchar pixel = ((__local uchar*)self)[index];
            minvalue = min( minvalue, pixel );
            maxvalue = max( maxvalue, pixel );
            sum += pixel;
        }
    } // for q

    sum += sum4.s0 + sum4.s1 + sum4.s2 + sum4.s3;    
    
    uchar2 min_a = min( minvalue4.s01, minvalue4.s23 );
    uchar  min_b = min( min_a.s0, min_a.s1 );
    minvalue = min( min_b, minvalue );    
    
    uchar2 max_a = max( maxvalue4.s01, maxvalue4.s23 );
    uchar  max_b = max( max_a.s0, max_a.s1 );
    maxvalue = max( max_b, maxvalue );    
        
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
    __local void * self,
    const clMat /*uint8_t*/  image,
    const Point2i center,
    int mapSize, 
    int m_patchSize,
    // clMat /*float*/m_wIn,
    clMat /*float*/m_wOut,
    // clMat /*float*/m_U,
    clMat wIn,
    clMat bIn,
    
    // result
    clMat /*float*/ result
    )
{

    int gangsize = get_local_size(0);    
  // assert( result.rows == 2 * mapSize + 1 );
  // assert( result.cols == 2 * mapSize + 1 );

  //  assert(wIn.rows == wIn_A.rows);
  //  assert(wIn.cols == m_U.rows);

    clMat /*float*/wOut_tmp = GetBlockFloat( self, m_wOut, 0, m_wOut.rows, 0, m_wOut.cols - 1 );

    {
        int localid = get_local_id(0);
//        for ( localid=0; localid<gangsize; localid++ ) {
          
            int worksize = (2 * mapSize + 1) * (2 * mapSize + 1);
          
            float bOut = GetValueFloat( self, m_wOut, 0, m_wOut.cols - 1 );
            
            for ( int work=0; work < (worksize / gangsize) + 1; work++ )
            {
              int index = work * gangsize + localid;
              if (index>=worksize) continue;
              int ncy = index / (2 * mapSize + 1);
              int cy  = center.y - mapSize + ncy;
              int ncx = index % (2 * mapSize + 1);
              int cx  = center.x - mapSize + ncx;

              clMat /*uint8_t*/   imagePatch = GetBlockChar( self, image, cy - m_patchSize, cy + m_patchSize + 1, cx - m_patchSize, cx + m_patchSize + 1 );

              normalization norm = normalizeSample( self, imagePatch );

              float dot = gemmFloatDirDirDir( self, wIn, imagePatch, -1.0f, bIn, -1.0f, norm, wOut_tmp /*, xOut */);

              SetValueFloat( self, result, ncy, ncx, 1.0f/( 1.0f + exp( - dot - bOut ) ) );
              
            } // for gangsize
            //       } // for localid
    } // localid block
  
    return;
} // generateResponseMap

int 
cvRound( float value )
{
    return (int)(value + (value >= 0 ? 0.5f : -0.5f));
} // cvRound


__kernel void
calculateMaps( __global void * self
             , __constant int memory_segments[]
             , int m_visibleLandmarks_size
             , int m_mapSize
             , __constant calcpackage packages[]
             , int buffersize
             , __local void * buffer
             )
{

    int idx = get_group_id(0);
    int gangsize = get_local_size(0);    
    int localid = get_local_id(0);

    // Copying the calculation into the local memory
    {        
        __global int * glob = ((__global int*)(self + memory_segments[idx]));
        __local  int * loc  = ((__local  int*) buffer);

        for ( int iter = 0; iter < (buffersize/4) / gangsize + 1; iter++ ) {
            int index = gangsize * iter + localid;
            if (index >= buffersize/4) continue;
            loc[index]=glob[index];
        }
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    calcpackage package = packages[idx];
    
    // Calculating
    Point2i center;
    
    float shape_x;
    float shape_y;
    
    shape_x = GetValueFloat( buffer, package.input.shape, 2*idx, 0 );
    shape_y = GetValueFloat( buffer, package.input.shape, 2*idx+1, 0 );
    
    center.x = cvRound(shape_x);
    center.y = cvRound(shape_y);
    
    generateResponseMap(
        buffer,
        package.input.alignedImage,
        center,
        m_mapSize,
        package.input.m_patchSize,
        package.input.m_wOut,
        package.input.wIn,
        package.input.bIn,
        
        // result
        package.output.responseMap );

    barrier(CLK_LOCAL_MEM_FENCE);
    // Copying back the result into the global memory
    {
        int start = package.output.responseMap.start;
        int size  = package.output.responseMap.step * package.output.responseMap.rows;
        
        __global float * glob = ((__global float*)(self + memory_segments[idx]));
        __local  float * loc  = ((__local  float*) buffer);
        
        for ( int iter = 0; iter < size / gangsize + 1; iter++ ) {         
            int index = start + gangsize * iter + localid;
            if (index >= start + size) continue;            
            glob[index]=loc[index];        
        }
    }
    
    return;
} // calculateMaps


// LuM end of file
