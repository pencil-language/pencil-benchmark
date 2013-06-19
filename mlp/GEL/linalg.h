/* ---------------------------------------------------------------------
* Copyright Realeyes OU 2012
*
* All rights reserved. Realeyes OU PROPRIETARY/CONFIDENTIAL.
* No part of this computer programs(s) may be used, reproduced,
* stored in any retrieval system, or transmitted, in any form or
* by any means, electronic, mechanical, photocopying,
* recording, or otherwise without prior written permission of
* Realeyes OU. Use is subject to license terms.
* ---------------------------------------------------------------------
*/

// UjoImro, 2012 this file is a bridge between gel and other linalg libraries

// IMPORTANT: Make sure, that You are using only one linalg library!!!
// NOTE: if no particular library is selected, the system falls back to opencv!

#ifndef CLM_LINALG_H_
#define CLM_LINALG_H_

#include <cmath>
#include <stdexcept>
//#include "Core/defines.h"
//#include "Core/static_assert.h"

#ifdef WITH_IPP
#  include <ipp.h>
#endif /* WITH_IPP */

#ifdef WITH_MKL
#  include <mkl_vml.h>
#  include <mkl_cblas.h>
#  include <mkl_lapacke.h>
#endif /* WITH_MKL */

#include <opencv2/core/core_c.h> // introduced first for the constants
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/types_c.h>

#include "defines.h"

#include <stdint.h>

namespace gel {

template <typename T0>
void copyTo(const cv::Mat_<T0>& source, cv::Mat_<T0>& dst);

/*!
  Multiplies two matrices and adds a third matrix with
  appropriate weights. The result is
  result = alpha * op(src1) * op(src2) + beta * op(src3)
  where
  op(src1):=src1 if trans1==false
  op(src1):=tr(src1) if trans1==true

  op(src2):=src2 if trans1==false
  op(src2):=tr(src2) if trans1==true


  all input matrices must have the appropriate sizes.
  the result matrix will be resized if necessary.
  result[k,m] := src1[k,l] * src2[l, m] + src3[k, m]

  \param src1 the first multiplier matrix
  \param src2 the second multiplier matrix
  \param alpha the quotient of the multiplication
  \param src3 the addend of the sum
  \param beta the quotient of the addition
  \param result the result of the multiplication
  \param trans1 the transpose operator trigger for src1
  \param trans2 the transpose operator trigger for src2
 */
template <typename T0>
void gemm(const cv::Mat_<T0>& src1, const cv::Mat_<T0>& src2,
          T0 alpha, const cv::Mat_<T0>& src3, T0 beta, cv::Mat_<T0>& result,
          bool trans1=false, bool trans2=false);


template <typename T0>
void exp(const cv::Mat_<T0>& source, cv::Mat_<T0>& result);


template <typename T0>
void solve(cv::Mat_<T0>& A, cv::Mat_<T0>& b);

template <typename src_t, typename dst_t>
void convertTo(const cv::Mat_<src_t>& src,cv::Mat_<dst_t> & dst,
           dst_t alpha,dst_t beta);

template <typename src_t, typename dst_t>
void convertTo(const cv::Mat_<src_t>& src, cv::Mat_<dst_t>& dst);

GEL_EXPORT cv::Mat_<uint8_t> BGRToGRAY(const cv::Mat_<cv::Vec<uint8_t, 3> > & frame );

template <class T0>
cv::Mat_<uint8_t>
convertToUINT8( const cv::Mat_<T0> & image, const T0 & scale, const T0 & shift);

template <class T0>
uint8_t saturate_cast( const T0 & t0 );
      
// *************************************************************************
// HELPER FUNCTIONS
// *************************************************************************

#ifdef WITH_MKL
// gemm helper function for connecting gel with mkl
template <typename T0>
void mkl_gemm( const CBLAS_ORDER Order, const  CBLAS_TRANSPOSE TransA,
          const CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N,
          const MKL_INT K, const T0 alpha, const T0* A,
          const MKL_INT lda, const T0* B, const MKL_INT ldb,
          const T0 beta, T0* C, const MKL_INT ldc )
{
    // STATIC_ASSERT(sizeof(T0)==0); // This value type is not supported or not bridged.
    throw std::runtime_error("error: wrong template specialization.");    
}
    
template <>
void mkl_gemm( const CBLAS_ORDER Order, const  CBLAS_TRANSPOSE TransA,
              const CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N,
              const MKL_INT K, const float alpha, const float* A,
              const MKL_INT lda, const float* B, const MKL_INT ldb,
              const float beta, float* C, const MKL_INT ldc );

    
template <>
void mkl_gemm( const CBLAS_ORDER Order, const  CBLAS_TRANSPOSE TransA,
              const CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N,
              const MKL_INT K, const double alpha, const double* A,
              const MKL_INT lda, const double* B, const MKL_INT ldb,
              const double beta, double* C, const MKL_INT ldc );
    

// exp helper functions
template <typename T0>
void mkl_exp(int n,const T0* source, T0* result)
{
    // STATIC_ASSERT(sizeof(T0)==0); // This value type is not supported or not bridged.
    throw std::runtime_error("error: wrong template specialization.");
}

template <>
void mkl_exp(int n,const float* source, float* result);

template <>
void mkl_exp(int n,const double* source, double* result);

// solve helper functions
template <typename T0>
void mkl_solve( int matrix_order, lapack_int m, lapack_int n,
                           lapack_int nrhs, T0* a,
                           lapack_int lda,  T0* b,
                           lapack_int ldb,  T0* s, T0 rcond,
                           lapack_int* rank  )
{
    // STATIC_ASSERT(sizeof(T0)==0); // This value type is not supported or not bridged.
    throw std::runtime_error("error: wrong template specialization.");
}

template <>
GEL_EXPORT void mkl_solve( int matrix_order, lapack_int m, lapack_int n,
                   lapack_int nrhs, float* a,
                   lapack_int lda,  float* b,
                   lapack_int ldb,  float* s, float rcond,
                   lapack_int* rank  );

template <>
GEL_EXPORT void mkl_solve( int matrix_order, lapack_int m, lapack_int n,
                   lapack_int nrhs, double* a,
                   lapack_int lda,  double* b,
                   lapack_int ldb,  double* s, double rcond,
                   lapack_int* rank  );



#endif /* WITH MKL */


#ifdef WITH_IPP

template <typename T0>
void ipp_copy( const T0 * pSrc, int srcStep, T0 * pDst, int dstStep, IppiSize roiSize )
{
    // STATIC_ASSERT(sizeof(T0)==0); // This value type is not supported or not bridged.
    throw std::runtime_error("error: wrong template specialization.");    
}

template <>
void ipp_copy( const float * pSrc, int srcStep, float * pDst, int dstStep, IppiSize roiSize );

template <>
void ipp_copy( const double * pSrc, int srcStep, double * pDst, int dstStep, IppiSize roiSize );

template <>
void ipp_copy( const int * pSrc, int srcStep, int * pDst, int dstStep, IppiSize roiSize );

template <>
void ipp_copy( const uint8_t * pSrc, int srcStep, uint8_t * pDst, int dstStep, IppiSize roiSize );

template <typename src_t, typename dst_t>
void ipps_convert( src_t* src, dst_t* dst, int32_t length )
{
    // STATIC_ASSERT(sizeof(src_t)==0); // This pixel_type is either not supported or not bridged.
    throw std::runtime_error("error: wrong template specialization.");    
}

template <>
void ipps_convert( uint8_t *, float *, int32_t );
template <>
void ipps_convert( uint8_t *, double *, int32_t );

template <typename T0>
void ipps_AddC( T0 C, T0 * src, int32_t length )
{
    // STATIC_ASSERT(sizeof(T0)==0); // This pixel_type is either not supported or not bridged.
    throw std::runtime_error("error: wrong template specialization.");    
}

template <>
void ipps_AddC( float, float *, int32_t );

template <>
void ipps_AddC( double, double *, int32_t );

template <typename T0>
void ipps_MulC( T0 C, T0 * src, int32_t length )
{
    // STATIC_ASSERT(sizeof(T0)==0); // This pixel_type is either not supported or not bridged.
    throw std::runtime_error("error: wrong template specialization.");    
}

template <>
void ipps_MulC( float, float *, int32_t );
template <>
void ipps_MulC( double, double *, int32_t );

#endif /* WITH_IPP */

}

// *************************************************************************
// Implementation
// *************************************************************************


template <typename T0>
inline void gel::gemm(const cv::Mat_<T0>& src1,
          const cv::Mat_<T0>& src2,
          T0 alpha,
          const cv::Mat_<T0>& src3,
          T0 beta,
          cv::Mat_<T0>& result,
          bool trans1,
          bool trans2)
{
    // PRINT(src1.rows);
    // PRINT(src1.cols);
    // PRINT(src2.rows);
    // PRINT(src2.cols);
    // PRINT(src3.rows);
    // PRINT(src3.cols);

    // check if the sizes are correct
    assert( src1.cols == src2.rows );
    assert( src1.rows == src3.rows );
    assert( src2.cols == src3.cols );

#ifdef WITH_MKL

    // resize the result if necessary
    result.create( src1.rows, src2.cols );
    // these variables specify if we want to transpose the
    // matrix before the multiplication
    CBLAS_TRANSPOSE ctrans1 = trans1 ? CblasTrans : CblasNoTrans;
    CBLAS_TRANSPOSE ctrans2 = trans2 ? CblasTrans : CblasNoTrans;

    //src3.copyTo(result);
    gel::copyTo(src3, result);

    mkl_gemm(
        CblasRowMajor, ctrans1, ctrans2, src1.rows, src2.cols, src1.cols, alpha,
        src1[0], src1.step1(),
        src2[0], result.step1(),
        beta,
        result[0], result.step1()
        );

#else
    cv::gemm( src1, src2, alpha, src3, beta, result );
#endif
}
    
template <typename T0>
inline void gel::copyTo(const cv::Mat_<T0>& source, cv::Mat_<T0>& dst)
{
#ifdef WITH_IPP
    typedef T0 pixel_type;
    assert(source.rows == dst.rows);
    assert(source.cols == dst.cols);

    const T0* source_ptr = source[0];
    T0* dst_ptr = dst[0];

    IppiSize roiSize;
    roiSize.width = source.cols;  // THESE DIMENSIONS ARE TWISTED; CAREFULL !!!
    roiSize.height = source.rows; // THESE DIMENSIONS ARE TWISTED; CAREFULL !!!

    ipp_copy( source_ptr, source.step, dst_ptr, dst.step, roiSize );

#else
    source.copyTo(dst);
#endif
}
    
template <typename T0>
inline void gel::exp(const cv::Mat_<T0>& source, cv::Mat_<T0>& result )
{
#ifdef WITH_MKL
    result.create( source.rows, source.cols );
    const T0* source_ptr = source[0];
    T0* result_ptr = result[0];
    int source_step = source.step1();
    int result_step = result.step1();

    if (source.isContinuous() && result.isContinuous()) {
        // we can calculate the exponential directly
        mkl_exp( source.rows * source.cols, source_ptr, result_ptr );
    } else {
        for ( int q=0; q < source.rows; q++) {
            mkl_exp( source.cols, &source_ptr[q * source_step], &result_ptr[q * result_step] );
        }
    }

#else
    cv::exp( source, result);
#endif
}

template <typename T0>
inline void gel::solve(cv::Mat_<T0> &A,cv::Mat_<T0>& b)
{
    // PRINT(A.rows);
    // PRINT(A.cols);
    // PRINT(b.rows);
    // PRINT(b.cols);
#ifdef WITH_MKL

    int rank; // the number of singular values
    T0* singulars = new T0[ std::max(1, std::min(A.rows, A.cols)) ]; // array for the singular values

    mkl_solve( LAPACK_ROW_MAJOR,
                   A.rows,
                   A.cols,
                   b.cols, // should be 1
                   A[0],
                   A.step1(),
                   b[0],
                   b.step1(),
                   singulars, // the singular values in order
                   static_cast<T0>(-1.), // the tolerance to zero (epsilon), if negative value is given machine precision is used
                   &rank);

    delete [] singulars;

#else
    cv::Mat_<T0> result;
    cv::solve(A, b, result, cv::DECOMP_SVD);
    b = result;
#endif
}

template <typename src_t, typename dst_t>
inline void gel::convertTo(const cv::Mat_<src_t>& src,
           cv::Mat_<dst_t> & dst,
           dst_t alpha,
           dst_t beta )
{
#ifdef WITH_IPP
    dst.create(src.rows, src.cols);

    if (dst.isContinuous()) {
        gel::ipps_convert(
            src[0],
            dst[0],
            src.rows*src.cols);

        ipps_MulC(alpha, dst[0], dst.rows * dst.cols );
        ipps_AddC(beta, dst[0], dst.rows * dst.cols );

    } else {
        for ( int q = 0; q < src.rows; q++ ) {
            gel::ipps_convert(
                src[q],
                dst[q],
                src.rows);

            ipps_MulC(alpha, dst[q], dst.cols );
            ipps_AddC(beta, dst[q], dst.cols );
        }
    }
#else
    src.convertTo(dst, dst.type(), alpha, beta);
#endif
}


template <typename src_t, typename dst_t>
inline void gel::convertTo(const cv::Mat_<src_t>& src, cv::Mat_<dst_t>& dst)
{
#ifdef WITH_IPP
    dst.create(src.rows, src.cols);

    if (dst.isContinuous()) {
        gel::ipps_convert(
                    src[0],
                    dst[0],
                    src.rows*src.cols);
    } else {
        for ( int q = 0; q < src.rows; q++) {
            gel::ipps_convert(
                        src[q],
                        dst[q],
                        src.rows);
        }
    }
#else
    src.convertTo(dst, dst.type(), static_cast<dst_t>(1.0), static_cast<dst_t>(0.0) );
#endif

}

template <typename src_t>
void mulAddC( cv::Mat_<src_t>& src, src_t alpha, src_t beta)
{
#ifdef WITH_IPP
    if (src.isContinuous()) {
        ipps_MulC(alpha, src[0], src.rows * src.cols );
        ipps_AddC(beta, src[0], src.rows * src.cols );
    } else {
        for ( int q = 0; q < src.rows; q++ ) {
            ipps_MulC(alpha, src[q], src.cols );
            ipps_AddC(beta, src[q], src.cols );
        }
    }

#else
    src = alpha * src + beta;
#endif
}

template <class T0>
cv::Mat_<uint8_t>
gel::convertToUINT8( const cv::Mat_<T0> & image, const T0 & scale, const T0 & shift )
{
    cv::Mat_<uint8_t> result;
    result.create(image.rows, image.cols);
    
    for (int q=0; q<image.rows; q++)
        for (int w=0; w<image.cols; w++)
            result(q,w) = gel::saturate_cast<T0>( image(q,w) * scale + shift );

    return result;    
} // convertToUINT8

template <class T0>
uint8_t
gel::saturate_cast( const T0 & t0 )
{
    T0 fractionalPart = 0.0L;
    T0 integralPart = 0.0L;
    
    if (t0>=255.0)
        return 255;
    if (t0<=0.0)
        return 0;

    fractionalPart = std::modf(t0, &integralPart);
    if (fractionalPart >= 0.5)
        return integralPart + 1;

    return integralPart;
} // saturate_cast

#endif /* CLM_LINALG_H_ */
