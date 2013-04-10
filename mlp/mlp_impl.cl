// -*- c++ -*-
// UjoImro, 2013
// Experimental code for the CARP Project
// Copyright (c) RealEyes, 2013
// This is an opencl-implementation of the PCA->MLP response map calculation

typedef struct {
    int x;
    int y;
} Point2i; // struct Point2i

typedef struct {
    int rows;
    int cols;
    int step;
    int start;    
} clMat; // struct clclMat /*float*/

typedef unsigned char uint8_t;

void
copyToFloat( clMat /*float*/ input, __global float * input_data, clMat /*float*/ output, __global float * output_data )
{
    // assert(input_data);
    // assert(output);
    // assert(output.data);    
    // assert(input.rows == output.rows);
    // assert(input.cols == output.cols);
    
    int q, w;
    for ( q=0; q<input.rows; q++ )
        for ( w=0; w<input.cols; w++ )
            output_data[ q * output.step + w + output.start ] =
        	input_data[ q * input.step + w + input.start ];
    
    return;
} // copyTo

void
transposeFloat( clMat /*float*/ input, __global float * input_data, clMat /*float*/ output, __global float * output_data )
{
    int q,w;
    // assert(output.rows == input.cols);
    // assert(output.cols == input.rows);
    for (q=0; q<input.rows; q++)
        for (w=0; w<input.cols; w++)
            output_data[ w * output.step + q + output.start ] 
        	= input_data[ q * input.step + w + input.start ];

    return;
} // transposeFloat

void
transposeFloatGang( int gangsize, clMat /*float*/ input, __global float * input_data, int localid, clMat /*float*/ output, __global float * output_data )
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
              
        output_data[ w * output.step + q + output.start ] 
            = input_data[ q * input.step + w + input.start ];
    }
} // transposeFloatGang


float
meanChar( clMat /* uint8_t */ input, __global uint8_t * input_data )
{
    // assert(input_data);    
    int q,w;
    float sum=0;
    float c = 0; // kahan summation
    
    for ( q=0; q<input.rows; q++ )
        for ( w=0; w<input.cols; w++ )
        {
            float y = input_data[ q * input.step + w + input.start ] - c;
            float t = sum + y;
            c = (t - sum) - y;        
            sum = t;        
        }
    
    return sum / ( input.rows * input.cols );
} // meanFloat

uint8_t
minChar( clMat /* uint8_t */ input, __global uint8_t * input_data )
{
    // assert(input_data);    
    int q,w;
    uint8_t minvalue = 255;

    for ( q=0; q<input.rows; q++ )
        for ( w=0; w<input.cols; w++ )
            minvalue = min( minvalue, input_data[ q * input.step + w + input.start ] );
    
    return minvalue;
} // minFloat

uint8_t
maxChar( clMat /* uint8_t */ input, __global uint8_t * input_data )
{
    // assert(input_data);
    int q,w;
    uint8_t maxvalue = 0;

    for ( q=0; q<input.rows; q++ )
        for ( w=0; w<input.cols; w++ )
            maxvalue = max( maxvalue, input_data[ q * input.step + w + input.start ] );
    
    return maxvalue;
} // maxFloat


clMat /* uint8_t */
GetBlockChar( clMat /* uint8_t */ self, __global uint8_t * self_data, int row_from, int row_to, int col_from, int col_to )
{    
    // assert(row_from>=0);
    // assert(col_from>=0);
    // assert(row_from<row_to);
    // assert(col_from<col_to);
    // assert(row_to<=self.rows);
//    printf("col_to = %d\n", col_to );
//    printf("self.cols = %d\n", self.cols );
    // assert(col_to<=self.cols);
    // assert(self_data);
    
    clMat /* uint8_t */ result;
    result.rows  = row_to - row_from;
    result.cols  = col_to - col_from;
    result.step  = self.step;
    result.start = self.start + row_from * self.step + col_from;
    // result_data  = self_data;

    return result;
}

clMat /*float*/
GetBlockFloat( clMat /*float*/ self, __global float * self_data, int row_from, int row_to, int col_from, int col_to )
{  
    // assert(row_from>=0);
    // assert(col_from>=0);
    // assert(row_from<row_to);
    // assert(col_from<col_to);
    // assert(row_to<=self.rows);
    // assert(col_to<=self.cols);
    // assert(self_data);

    clMat /*float*/ result;
    result.rows  = row_to - row_from;
    result.cols  = col_to - col_from;
    result.step  = self.step;
    result.start = self.start + row_from * self.step + col_from;
    // result_data  = self_data;

    return result;
}

void
convertFromCharToFloat( clMat /* uint8_t */ from, __global uint8_t * from_data, float quotient, float shift, clMat /*float*/ to, __global float * to_data )
{
    // assert(from.rows == to.rows);
    // assert(from.cols == to.cols);
    
    int q, w;
    for ( q=0; q<from.rows; q++ )
        for ( w=0; w<from.cols; w++ )
            to_data[ q * to.step + w + to.start ] = 
        	quotient * from_data[ q * from.step + w + from.start ] + shift;
    return;
}

clMat /*float*/
reshapeFloat( clMat /*float*/ self, __global float * self_data, int new_rows )
{

    // assert(self.cols == self.step);
    // assert( (self.cols * self.rows) % new_rows == 0 );

    clMat /*float*/ result;
    result.rows = new_rows;
    result.cols = self.cols * self.rows / new_rows;
    result.step = result.cols;
    result.start = self.start;
    // result_data = self_data;

    return result;
} // reshapeFloat


// returns alpha*A*B + beta * C
void 
gemmFloatDirDirDir(
    clMat /*float*/ A, __global float * A_data,
    clMat /*float*/ B, __global float * B_data,
    float alpha,
    clMat /*float*/ C, __global float * C_data,
    float beta,
    clMat /*float*/ result, __global float * result_data )
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
                    float y = A_data[ q * A.step + e + A.start ] * B_data[ e * B.step + w + B.start ] - c;
                    float t = sum + y;
                    c = (t - sum) - y;
                    sum = t;              
                }
                
                result_data[ q * result.step + w + result.start ] = alpha * sum  + beta * C_data[ q * C.step + w + C.start ];
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
                    float y = A_data[ q * A.step + e + A.start ] * B_data[ e * B.step + w + B.start ] - c;
                    float t = sum + y;
                    c = (t - sum) - y;
                    sum = t;              
                }
                
                result_data[ q * result.step + w + result.start ] = alpha * sum;
            }        
    }
    
    return;
} // gemmFloatDirDirDir

void 
gemmFloatDirDirDirGang(
    int gangsize, 
    clMat /*float*/ A, __global float * A_data,
    clMat /*float*/ B, __global float * B_data,
    float alpha,
    clMat /*float*/ C, __global float * C_data,
    float beta,
    int localid,
    clMat /*float*/ result, __global float * result_data )
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
                float y = A_data[ q * A.step + e + A.start ] * B_data[ e * B.step + w + B.start ] - c;
                float t = sum + y;
                c = (t - sum) - y;
                sum = t;              
            }
                
            result_data[ q * result.step + w + result.start ] = alpha * sum  + beta * C_data[ q * C.step + w + C.start ];
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
                float y = A_data[ q * A.step + e + A.start ] * B_data[ e * B.step + w + B.start ] - c;
                float t = sum + y;
                c = (t - sum) - y;
                sum = t;              
            }
            
            result_data[ q * result.step + w + result.start ] = alpha * sum;
        }
    }
    
    return;
} // gemmFloatDirDirDirGang


void 
gemmFloatDirTransDirGang(
    int gangsize, 
    clMat /*float*/ A, __global float * A_data,
    clMat /*float*/ B, __global float * B_data,
    float alpha,
    clMat /*float*/ C, __global float * C_data,
    float beta,
    int localid,
    clMat /*float*/ result, __global float * result_data )
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
                float y = A_data[ q * A.step + e + A.start ] * B_data[ w * B.step + e + B.start ] - c;
                float t = sum + y;
                c = (t - sum) - y;
                sum = t;              
            }
                
            result_data[ q * result.step + w + result.start ] = alpha * sum  + beta * C_data[ q * C.step + w + C.start ];
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
                float y = A_data[ q * A.step + e + A.start ] * B_data[ w * B.step + e + B.start ] - c;
                float t = sum + y;
                c = (t - sum) - y;
                sum = t;              
            }
            
            result_data[ q * result.step + w + result.start ] = alpha * sum;
        }
    }
    
    return;
} // gemmFloatDirTransDirGang


void
expFloat(
    clMat /*float*/ input, __global float * input_data,
    clMat /*float*/ output, __global float * output_data )
{
    // assert(input.rows == output.rows);
    // assert(input.cols == output.cols);

    int q, w;
    for ( q=0; q<input.rows; q++ )
        for ( w=0; w<input.cols; w++ )
            output_data[ q * output.step + w + output.start ] = 
        	exp(input_data[ q * input.step + w + input.start ]);

    return;
}

void
addFloat(
    clMat /*float*/ input, __global float * input_data,
    float val,
    clMat /*float*/ output, __global float * output_data )
{
    // assert(input.rows == output.rows);
    // assert(input.cols == output.cols);

    int q, w;
    for ( q=0; q<input.rows; q++ )
        for ( w=0; w<input.cols; w++ )
            output_data[ q * output.step + w + output.start ] = 
        	val + input_data[ q * input.step + w + input.start ];

    return;
}

void 
divideFloat( float val,
             clMat /*float*/ input, __global float * input_data,
             clMat /*float*/ output, __global float * output_data )
{
    // assert(input.rows == output.rows);
    // assert(input.cols == output.cols);

    int q, w;
    for ( q=0; q<input.rows; q++ )
        for ( w=0; w<input.cols; w++ )
            output_data[ q * output.step + w + output.start ] = 
        	val / input_data[ q * input.step + w + input.start ];

    return;
}

void
subtractFloat( clMat /*float*/ input, __global float * input_data,
               float val,
               clMat /*float*/ output, __global float * output_data )
{
    // assert(input.rows == output.rows);
    // assert(input.cols == output.cols);

    int q, w;
    for ( q=0; q<input.rows; q++ )
        for ( w=0; w<input.cols; w++ )
            output_data[ q * output.step + w + output.start ] = 
        	input_data[ q * input.step + w + input.start ] - val;
    
    return;
}

float
GetValueFloat( clMat /*float*/ self, __global float * self_data, int row, int col )
{
    return self_data[ row * self.step + col + self.start ];
    // return -1231.;
}

void
SetValueFloat( clMat /*float*/ self, __global float * self_data, int row, int col, float value )
{
    self_data[ row * self.step + col + self.start ] = value;
    return;    
}

float dotProductDirDir(
    clMat /*float*/ A, __global float * A_data,
    clMat /*float*/ B, __global float * B_data )
{
    // assert( A.cols == 1 );
    // assert( B.cols == 1 );
    // assert( A.rows == B.rows );

    float result = 0.;
    float c = 0.;  

    int q;
    for (q=0; q<A.rows; q++) {
        float y = A_data[ q * A.step + A.start ] * B_data[ q * B.step + B.start ] - c;
        float t = result + y;
        c = (t - result) - y;
        result = t ;
    }
  
    return result;
}

float dotProductTransDir(
    clMat /*float*/ A, __global float * A_data,
    clMat /*float*/ B, __global float * B_data )
{
    // assert( A.rows == 1 );
    // assert( B.cols == 1 );
    // assert( A.cols == B.rows );

    float result = 0.;
    float c = 0.;  

    int q;
    for (q=0; q<A.cols; q++) {
        float y = A_data[ q + A.start ] * B_data[ q * B.step + B.start ] - c;
        float t = result + y;
        c = (t - result) - y;
        result = t ;
    }
  
    return result;
}

void normalizeSample(
    clMat /* uint8_t */ image, __global uint8_t * image_data,
    clMat /*float*/* result, __global float * result_data )
{
    // assert(result.cols == image.cols);
    // assert(result.rows == image.rows);
    
    float sampleMean = meanChar(image, image_data);
    float sampleMin  = minChar(image, image_data);
    float sampleMax  = maxChar(image, image_data);
    
    sampleMax -= sampleMean;
    sampleMin -= sampleMean;
    
    sampleMax = fmax( fabs(sampleMin), fabs(sampleMax));
    
    if (sampleMax == 0.0) sampleMax = 1.0;

    convertFromCharToFloat( image, image_data, 1.0/sampleMax, -(1.0/sampleMax)*sampleMean, *result, result_data );
    
    *result = reshapeFloat( *result, result_data, image.rows * image.cols );
  
    return;
} // normalizeSample

void
generateResponseMap(
    int gangsize, 
    clMat /* uint8_t */ image,
    __global uint8_t * image_data,
    Point2i center,
    int mapSize, 
    int m_patchSize,
    clMat /*float*/ m_wIn,
    __global float * m_wIn_data,
    clMat /*float*/ m_wOut,
    __global float * m_wOut_data,
    clMat /*float*/ m_U,
    __global float * m_U_data,

    // result
    clMat /*float*/ result,
    __global float * result_data
    )
{

    // assert(result.rows == 2 * mapSize + 1);
    // assert(result.cols == 2 * mapSize + 1);

    clMat /*float*/ wIn_A = GetBlockFloat( m_wIn, m_wIn_data, 0, m_wIn.rows, 0, m_wIn.cols - 1 );
    clMat /*float*/ wIn;   // = CreateMatFloat( wIn_A.rows, m_U.rows );
    clMat /*float*/ bIn   = GetBlockFloat( m_wIn, m_wIn_data, 0, m_wIn.rows, m_wIn.cols - 1, m_wIn.cols );
    clMat /*float*/ wOut_tmp = GetBlockFloat( m_wOut, m_wOut_data, 0, m_wOut.rows, 0, m_wOut.cols - 1 );

    {
        int localid = 0;
///!!!        for ( localid=0; localid<gangsize; localid++ )           
//!!!            gemmFloatDirTransDirGang( wIn_A, m_U, 1.0, wIn, 0.0, localid, &wIn );
    }
  
    {
        int localid = 0;
        for ( localid=0; localid<gangsize; localid++ ) {
          
            int work = 0;      
            int worksize = (2 * mapSize + 1) * (2 * mapSize + 1);
          
            float bOut = GetValueFloat( m_wOut, m_wOut_data, 0, m_wOut.cols - 1 );
          
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
          
                clMat /* uint8_t */  imagePatch = GetBlockChar( image, image_data, cy - m_patchSize, cy + m_patchSize + 1, cx - m_patchSize, cx + m_patchSize + 1 );
                clMat /*float*/ patch; // = CreateMatFloat( imagePatch.rows, imagePatch.cols );
          
//!!!                normalizeSample(imagePatch, &patch);
          
                clMat /*float*/ xOut; // = CreateMatFloat( bIn.rows, bIn.cols );
          
///!!!                gemmFloatDirDirDir( wIn, patch, -1.0, bIn, -1.0, &xOut );
          
                clMat /*float*/ e; // = CreateMatFloat(xOut.rows, xOut.cols);
          
///!!!                expFloat( xOut, &e );
          
///!!!                addFloat( e, 1.0, &xOut );
///!!!                divideFloat( 2.0, xOut, &e);
///!!!                addFloat( e, -1.0, &xOut);
          
///!!!                SetValueFloat( result, ncy, ncx, 1./( 1. + exp(- dotProductTransDir(wOut_tmp, xOut) - bOut ) ) );
          
                // freeMatFloat(&e);
                // freeMatFloat(&xOut);
                // freeMatFloat(&patch);
            } // for localid 
        } // for gangsize
    } // localid block
  
//  // assert(false);
  
    // freeMatFloat(&wIn);
    // end of classic impl
    return;
} // generateResponseMap

int 
cvRound( float value )
{
    return (int)(value + (value >= 0 ? 0.5 : -0.5));
} // cvRound

__kernel void
calculateMaps(
    int m_visibleLandmarks_size, 
    int m_mapSize, 
    clMat /* uint8_t */ alignedImage, __global uint8_t * alignedImage_data,
    clMat /*float*/ shape, __global float * shape_data,
    __global int m_patchSizes[],
    __global clMat /*float*/ m_wIns[], __global float * m_wIns_data[],
    __global clMat /*float*/ m_wOuts[], __global float * m_wOuts_data[],
    __global clMat /*float*/ m_Us[], __global float * m_Us_data[],
    
    // results
    clMat /*float*/ responseMaps[], __global float * responseMaps_data[] )
{

    const int gangsize = get_local_size(0);
    const int localid = get_local_id(0);
    const int groupid = get_group_id(0);

    if (groupid>0) return;
    if (localid>0) return;    
    
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

	shape_x = GetValueFloat( shape, shape_data, 2*idx, 0 );
	shape_y = GetValueFloat( shape, shape_data, 2*idx+1, 0 );
	center.x = cvRound(shape_x);
	center.y = cvRound(shape_y);
        
	// responseMaps[q] = m_classifiers[idx].generateResponseMap( alignedImage, center, m_mapSize );
//        printf("alignedImage.rows = %d\n", alignedImage.rows );
//        printf("alignedImage.cols = %d\n", alignedImage.cols );
        
	generateResponseMap(
            gangsize,
            alignedImage,
            alignedImage_data,
            center,
            m_mapSize,
            m_patchSizes[idx],
            m_wIns[idx],
            m_wIns_data[idx],
            m_wOuts[idx],
            m_wOuts_data[idx],
            m_Us[idx],
            m_Us_data[idx],
            responseMaps[q],
            responseMaps_data[q] );

    } // for cycle

    // printf("calculateMaps finished\n");
    return;
} // calculateMaps


// LuM end of file
