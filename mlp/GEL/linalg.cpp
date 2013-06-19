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

// UjoImro, 2012; this file is a bridge between gel and other linalg libraries

#define WITH_IPP
#define WITH_MKL

#include "linalg.h"

#ifdef WITH_MKL
#  include <mkl_vml.h>
#  include <mkl_cblas.h>
#  include <mkl_lapacke.h>
#endif /* WITH_MKL */

#ifdef WITH_IPP
#  include <ipp.h>
#endif /* WITH_IPP */

#include <iostream>

#include <opencv2/imgproc/imgproc.hpp>

#ifdef WITH_MKL
    
template <>
void gel::mkl_gemm( const CBLAS_ORDER Order, const  CBLAS_TRANSPOSE TransA,
              const CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N,
              const MKL_INT K, const float alpha, const float* A,
              const MKL_INT lda, const float* B, const MKL_INT ldb,
              const float beta, float* C, const MKL_INT ldc )
{
    cblas_sgemm( Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc );
}

    
template <>
void gel::mkl_gemm( const CBLAS_ORDER Order, const  CBLAS_TRANSPOSE TransA,
              const CBLAS_TRANSPOSE TransB, const MKL_INT M, const MKL_INT N,
              const MKL_INT K, const double alpha, const double* A,
              const MKL_INT lda, const double* B, const MKL_INT ldb,
              const double beta, double* C, const MKL_INT ldc )
{
    cblas_dgemm( Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc );
}

#endif

#ifdef WITH_IPP
template <>
void gel::ipp_copy( const float * pSrc, int srcStep, float * pDst, int dstStep, IppiSize roiSize )
{
    ippiCopy_32f_C1R(pSrc, srcStep, pDst, dstStep, roiSize);
}

template <>
void gel::ipp_copy( const double * pSrc, int srcStep, double * pDst, int dstStep, IppiSize roiSize )
{
    // note: this specialization is a hack there is no ippicopy
    // for doubles, so we copy integers with the double apparent
    // size.

    IppiSize doubleSize;
    doubleSize.height = roiSize.height;
    doubleSize.width = 2*roiSize.width;

    ippiCopy_32s_C1R( reinterpret_cast<const Ipp32s*>(pSrc),
                      srcStep,
                      reinterpret_cast<Ipp32s*>(pDst), dstStep, doubleSize);
}

template <>
void gel::ipp_copy( const int * pSrc, int srcStep, int * pDst, int dstStep, IppiSize roiSize )
{
    ippiCopy_32s_C1R(pSrc, srcStep, pDst, dstStep, roiSize);
}

template <>
void gel::ipp_copy( const uint8_t * pSrc, int srcStep, uint8_t * pDst, int dstStep, IppiSize roiSize )
{
    ippiCopy_8u_C1R(pSrc, srcStep, pDst, dstStep, roiSize);
}

#endif


#ifdef WITH_MKL

template <>
void gel::mkl_exp(int n,const float * source, float * result)
{
    vsExp(n, source, result);
}

template <>
void gel::mkl_exp(int n,const double * source, double * result)
{
    vdExp(n, source, result);
}

#endif




#ifdef WITH_MKL

template <>
void gel::mkl_solve(int matrix_order, lapack_int m, lapack_int n,
                       lapack_int nrhs, float* a,
                       lapack_int lda,  float* b,
                       lapack_int ldb,  float* s, float rcond,
                       lapack_int* rank)
{
    LAPACKE_sgelss(matrix_order, m, n, nrhs, a, lda, b, ldb, s, rcond, rank);
}


template <>
void gel::mkl_solve(int matrix_order, lapack_int m, lapack_int n,
               lapack_int nrhs, double* a,
               lapack_int lda,  double* b,
               lapack_int ldb,  double* s, double rcond,
               lapack_int* rank)
{
    LAPACKE_dgelss(matrix_order, m, n, nrhs, a, lda, b, ldb, s, rcond, rank);
}


#endif


#ifdef WITH_IPP

template <>
void gel::ipps_convert( uint8_t* src, float* dst, int32_t len )
{
    ippsConvert_8u32f(src, dst, len);
}

template <>
void gel::ipps_convert( uint8_t* src, double* dst, int32_t len )
{
    //ippsConvert_8u64f(src, dst, len);
    float * tmp = new float[len];
    ippsConvert_8u32f(src, tmp, len);
    ippsConvert_32f64f(tmp, dst, len);
    delete [] tmp;
}

template <>
void gel::ipps_AddC( float C, float* src, int32_t len )
{
    ippsAddC_32f_I( C, src, len );
}


template <>
void gel::ipps_AddC( double C, double* src, int32_t len )
{
    ippsAddC_64f_I( C, src, len );
}

template <>
void gel::ipps_MulC( float C, float* src, int32_t len )
{
    ippsMulC_32f_I( C, src, len );
}

template <>
void gel::ipps_MulC( double C, double * src, int32_t len )
{
    ippsMulC_64f_I( C, src, len );
}

#endif

cv::Mat_<uint8_t> gel::BGRToGRAY( const cv::Mat_<cv::Vec<uint8_t, 3> > & frame )
{
    assert(frame.channels()==3);
    cv::Mat_<uint8_t> result;
    cv::cvtColor(frame,result,CV_BGR2GRAY);
    return result;
}
