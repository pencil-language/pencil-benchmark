//-------------------------------------------------------------------------//
//                                                                         //
//  This benchmark is an OpenMP C version of the NPB FT code. This OpenMP  //
//  C version is developed by the Center for Manycore Programming at Seoul //
//  National University and derived from the OpenMP Fortran versions in    //
//  "NPB3.3-OMP" developed by NAS.                                         //
//                                                                         //
//  Permission to use, copy, distribute and modify this software for any   //
//  purpose with or without fee is hereby granted. This software is        //
//  provided "as is" without express or implied warranty.                  //
//                                                                         //
//  Information on NPB 3.3, including the technical report, the original   //
//  specifications, source code, results and information on how to submit  //
//  new results, is available at:                                          //
//                                                                         //
//           http://www.nas.nasa.gov/Software/NPB/                         //
//                                                                         //
//  Send comments or suggestions for this OpenMP C version to              //
//  cmp@aces.snu.ac.kr                                                     //
//                                                                         //
//          Center for Manycore Programming                                //
//          School of Computer Science and Engineering                     //
//          Seoul National University                                      //
//          Seoul 151-744, Korea                                           //
//                                                                         //
//          E-mail:  cmp@aces.snu.ac.kr                                    //
//                                                                         //
//-------------------------------------------------------------------------//

//-------------------------------------------------------------------------//
// Authors: Sangmin Seo, Jungwon Kim, Jun Lee, Jeongho Nah, Gangwon Jo,    //
//          and Jaejin Lee                                                 //
//-------------------------------------------------------------------------//

//---------------------------------------------------------------------
// FT benchmark
//---------------------------------------------------------------------

#include "global.h"
#include <math.h>

static dcomplex ty1[MAXDIM][FFTBLOCKPAD_DEFAULT];
static dcomplex ty2[MAXDIM][FFTBLOCKPAD_DEFAULT];
#pragma omp threadprivate(ty1,ty2)

//---------------------------------------------------------------------
// evolve u0 -> u1 (t time steps) in fourier space
//---------------------------------------------------------------------
void evolve(int d1, int d2, int d3,
	    dcomplex u0[static const restrict d3][d2][d1+1],
	    dcomplex u1[static const restrict d3][d2][d1+1],
	    double twiddle[static const restrict d3][d2][d1+1])
{
  int i, j, k;

//  #pragma omp parallel for default(shared) private(i,j,k)
#pragma scop
  for (k = 0; k < d3; k++) {
    for (j = 0; j < d2; j++) {
      for (i = 0; i < d1; i++) {
        u0[k][j][i].real = u0[k][j][i].real*twiddle[k][j][i];
	u0[k][j][i].imag = u0[k][j][i].imag*twiddle[k][j][i];
        u1[k][j][i] = u0[k][j][i];
      }
    }
  }
#pragma endscop
}

//---------------------------------------------------------------------
// compute function from local (i,j,k) to ibar^2+jbar^2+kbar^2 
// for time evolution exponent. 
//---------------------------------------------------------------------
void compute_indexmap(int d1, int d2, int d3,
   		      double twiddle[static const restrict d3][d2][d1+1])
{
  int i, j, k, kk, kk2, jj, kj2, ii;
  double ap;

  //---------------------------------------------------------------------
  // basically we want to convert the fortran indices 
  //   1 2 3 4 5 6 7 8 
  // to 
  //   0 1 2 3 -4 -3 -2 -1
  // The following magic formula does the trick:
  // mod(i-1+n/2, n) - n/2
  //---------------------------------------------------------------------

  ap = -4.0 * ALPHA * PI * PI;

//  #pragma omp parallel for default(shared) private(i,j,k,kk,kk2,jj,kj2,ii)
#pragma scop
  for (k = 0; k < dims[2]; k++) {
    kk = ((k + NZ/2) % NZ) - NZ/2;
    kk2 = kk*kk;
    for (j = 0; j < dims[1]; j++) {
      jj = ((j + NY/2) % NY) - NY/2;
      kj2 = jj*jj + kk2;
      for (i = 0; i < dims[0]; i++) {
        ii = ((i + NX/2) % NX) - NX/2;
        twiddle[k][j][i] = exp(ap * (double)(ii*ii+kj2));
      }
    }
  }
#pragma endscop
}

//---------------------------------------------------------------------
// note: args x1, x2 must be different arrays
// note: args for cfftsx are (direction, layout, xin, xout, scratch)
//       xin/xout may be the same and it can be somewhat faster
//       if they are
//---------------------------------------------------------------------

//void fft(int dir, dcomplex x1[NTOTALP], dcomplex x2[NTOTALP])
void fft(int dir, int d1, int d2, int d3,
	 dcomplex x1[static const restrict d3][d2][d1+1],
	 dcomplex x2[static const restrict d3][d2][d1+1])
{
  if (dir == 1) {
    cffts1_inlined(1, dims[0], dims[1], dims[2], x1, x1);
    cffts2_inlined(1, dims[0], dims[1], dims[2], x1, x1);
    cffts3_inlined(1, dims[0], dims[1], dims[2], x1, x2);
  } else {
    cffts3_inlined(-1, dims[0], dims[1], dims[2], x1, x1);
    cffts2_inlined(-1, dims[0], dims[1], dims[2], x1, x1);
    cffts1_inlined(-1, dims[0], dims[1], dims[2], x1, x2);
  }
}

void cffts1(int is, int d1, int d2, int d3,
	    dcomplex x[static const restrict d3][d2][d1+1],
	    dcomplex xout[static const restrict d3][d2][d1+1])
{
  int logd1;
  int i, j, k, jj;

  logd1 = ilog2(d1);

  /* Requires array expansion.  */
  #pragma omp parallel for default(shared) private(i,j,k,jj)
  for (k = 0; k < d3; k++) {
    for (jj = 0; jj <= d2 - fftblock; jj += fftblock) {
      for (j = 0; j < fftblock; j++) {
        for (i = 0; i < d1; i++) {
          ty1[i][j] = x[k][j+jj][i];
        }
      }

      cfftz(is, logd1, d1, ty1, ty2);

      for (j = 0; j < fftblock; j++) {
        for (i = 0; i < d1; i++) {
          xout[k][j+jj][i] = ty1[i][j];
        }
      }
    }
  }
}

int calc_lk(int l, int shift)
{
  return (1 << (l - shift));
}

int calc_li(int l, int logd1, int shift)
{
  return (1 << (logd1 - l - shift));
}

dcomplex calc_u1(int is, int ku, int i0)
{
  dcomplex u1;

  if (is >= 1) {
    u1 = u[ku+i0];
  } else {
    u1 = dconjg(u[ku+i0]);        
  }
 
  return u1;	  
}

void cffts1_inlined(int is, int d1, int d2, int d3,
	    dcomplex x[static const restrict d3][d2][d1+1],
	    dcomplex xout[static const restrict d3][d2][d1+1])
{
  int logd1;
  int i, j, jj;
  int i0, j0, l;
  int k, kk, n1, li, lj, lk, ku, i11, i12, i21, i22;
  dcomplex u1, x11, x21;

  logd1 = ilog2(d1);

  /* Requires array expansion.  */
//  #pragma omp parallel for default(shared) private(i,j,k,jj)
#pragma scop
  for (k = 0; k < d3; k++) {
    for (jj = 0; jj <= d2 - FFTBLOCK_DEFAULT; jj += FFTBLOCK_DEFAULT) {
      for (j = 0; j < FFTBLOCK_DEFAULT; j++) {
        for (i = 0; i < d1; i++) {
          ty1[i][j] = x[k][j+jj][i];
        }
      }

      //------------------ cfftz_inlined    
      // Perform one variant of the Stockham FFT.
      for (l = 1; l <= logd1; l += 2) {
        // Set initial parameters.
        n1 = d1 / 2;
        lk = calc_lk(l, 1);
        li = calc_li(l, logd1, 0);
        lj = 2 * lk;
        ku = li;
        for (i0 = 0; i0 <= li - 1; i0++) {
          i11 = i0 * lk;
          i12 = i11 + n1;
          i21 = i0 * lj;
          i22 = i21 + lk;
	  u1 = calc_u1(is, ku, i0);

          // This loop is vectorizable.
          for (kk = 0; kk <= lk - 1; kk++) {
            for (j0 = 0; j0 < FFTBLOCK_DEFAULT; j0++) {
              x11 = ty1[i11+kk][j0];
              x21 = ty1[i12+kk][j0];
              ty2[i21+kk][j0].real = x11.real + x21.real;
              ty2[i21+kk][j0].imag = x11.imag + x21.imag;
              ty2[i22+kk][j0].real = (u1.real * (x11.real - x21.real)) -
		      (u1.imag*(x11.imag - x21.imag));
              ty2[i22+kk][j0].imag = (u1.real*(x11.imag - x21.imag)) +
		      (u1.imag*(x11.real - x21.real));
            }
          }
        }  

        // Copy Y to X in the last iteration,
        if (l == logd1) {
          for (j0 = 0; j0 < d1; j0++) {
            for (i0 = 0; i0 < FFTBLOCK_DEFAULT; i0++) {
	      ty1[j0][i0] = ty2[j0][i0];
            }
          }
        } else {
          // Set initial parameters.
          n1 = d1 / 2;
          lk = calc_lk(l, 0);
          li = calc_li(l, logd1, 1); 
          lj = 2 * lk;
          ku = li;
          for (i0 = 0; i0 <= li - 1; i0++) {
            i11 = i0 * lk;
            i12 = i11 + n1;
            i21 = i0 * lj;
            i22 = i21 + lk;
	    u1 = calc_u1(is, ku, i0);

            // This loop is vectorizable.
            for (kk = 0; kk <= lk - 1; kk++) {
              for (j0 = 0; j0 < FFTBLOCK_DEFAULT; j0++) {
	        x11 = ty2[i11+kk][j0];
                x21 = ty2[i12+kk][j0];
                ty1[i21+kk][j0].real = x11.real + x21.real;
                ty1[i21+kk][j0].imag = x11.imag + x21.imag;
                ty1[i22+kk][j0].real = (u1.real * (x11.real - x21.real)) -
		      (u1.imag*(x11.imag - x21.imag));
                ty1[i22+kk][j0].imag = (u1.real*(x11.imag - x21.imag)) +
		      (u1.imag*(x11.real - x21.real));
              }
            }
          }
        }
      }
      // -------------- end of cfftz_inlined

      for (j = 0; j < FFTBLOCK_DEFAULT; j++) {
        for (i = 0; i < d1; i++) {
          xout[k][j+jj][i] = ty1[i][j];
        }
      }
    }
  }
#pragma endscop  
}


void cffts2(int is, int d1, int d2, int d3,
	    dcomplex x[static const restrict d3][d2][d1+1],
	    dcomplex xout[static const restrict d3][d2][d1+1])
{
  int logd2;
  int i, j, k, ii;

  logd2 = ilog2(d2);

  #pragma omp parallel for default(shared) private(i,j,k,ii)
  for (k = 0; k < d3; k++) {
    for (ii = 0; ii <= d1 - fftblock; ii += fftblock) {
      for (j = 0; j < d2; j++) {
        for (i = 0; i < fftblock; i++) {
          ty1[j][i] = x[k][j][i+ii];
        }
      }

      cfftz(is, logd2, d2, ty1, ty2);

      for (j = 0; j < d2; j++) {
        for (i = 0; i < fftblock; i++) {
          xout[k][j][i+ii] = ty1[j][i];
        }
      }
    }
  }
}


void cffts2_inlined (int is, int d1, int d2, int d3,
	             dcomplex x[static const restrict d3][d2][d1+1],
	             dcomplex xout[static const restrict d3][d2][d1+1])
{
  int logd2;
  int i, j, k, ii;

  logd2 = ilog2(d2);

//  #pragma omp parallel for default(shared) private(i,j,k,ii)
#pragma scop
  for (k = 0; k < d3; k++) {
    for (ii = 0; ii <= d1 - FFTBLOCK_DEFAULT; ii += FFTBLOCK_DEFAULT) {
      for (j = 0; j < d2; j++) {
        for (i = 0; i < FFTBLOCK_DEFAULT; i++) {
          ty1[j][i] = x[k][j][i+ii];
        }
      }

      //------------------ cfftz_inlined
      int i0, j0, l;

      // Perform one variant of the Stockham FFT.
      for (l = 1; l <= logd2; l += 2) {
        int k, n1, li, lj, lk, ku, i11, i12, i21, i22;
        dcomplex u1, x11, x21;

        // Set initial parameters.
        n1 = d2 / 2;
        lk = 1 << (l - 1);
        li = 1 << (logd2 - l);
        lj = 2 * lk;
        ku = li;
        for (i0 = 0; i0 <= li - 1; i0++) {
          i11 = i0 * lk;
          i12 = i11 + n1;
          i21 = i0 * lj;
          i22 = i21 + lk;
          if (is >= 1) {
            u1 = u[ku+i0];
          } else {
            u1 = dconjg(u[ku+i0]);
          }

          // This loop is vectorizable.
          for (k = 0; k <= lk - 1; k++) {
            for (j0 = 0; j0 < FFTBLOCK_DEFAULT; j0++) {
              x11 = ty1[i11+k][j0];
              x21 = ty1[i12+k][j0];
              ty2[i21+k][j0] = dcmplx_add(x11, x21);
              ty2[i22+k][j0] = dcmplx_mul(u1, dcmplx_sub(x11, x21));
            }
          }
        }  

        // Copy Y to X in the last iteration,
        if (l == logd2) {
          for (j0 = 0; j0 < d2; j0++) {
            for (i0 = 0; i0 < FFTBLOCK_DEFAULT; i0++) {
	      ty1[j0][i0] = ty2[j0][i0];
            }
          }
        } else {
          // Set initial parameters.
          n1 = d2 / 2;
          lk = 1 << (l);
          li = 1 << (logd2 - l - 1);
          lj = 2 * lk;
          ku = li;
          for (i0 = 0; i0 <= li - 1; i0++) {
            i11 = i0 * lk;
            i12 = i11 + n1;
            i21 = i0 * lj;
            i22 = i21 + lk;
            if (is >= 1) {
              u1 = u[ku+i0];
            } else {
              u1 = dconjg(u[ku+i0]);
            }

            // This loop is vectorizable.
            for (k = 0; k <= lk - 1; k++) {
              for (j0 = 0; j0 < FFTBLOCK_DEFAULT; j0++) {
	        x11 = ty2[i11+k][j0];
                x21 = ty2[i12+k][j0];
                ty1[i21+k][j0] = dcmplx_add(x11, x21);
                ty1[i22+k][j0] = dcmplx_mul(u1, dcmplx_sub(x11, x21));
              }
            }
          }
        }
      }
      // -------------- end of cfftz_inlined

      for (j = 0; j < d2; j++) {
        for (i = 0; i < FFTBLOCK_DEFAULT; i++) {
          xout[k][j][i+ii] = ty1[j][i];
        }
      }
    }
  }
#pragma endscop  
}


void cffts3(int is, int d1, int d2, int d3,
	    dcomplex x[static const restrict d3][d2][d1+1],
	    dcomplex xout[static const restrict d3][d2][d1+1])
{
  int logd3;
  int i, j, k, ii;

  logd3 = ilog2(d3);

  #pragma omp parallel for default(shared) private(i,j,k,ii)
  for (j = 0; j < d2; j++) {
    for (ii = 0; ii <= d1 - fftblock; ii += fftblock) {
      for (k = 0; k < d3; k++) {
        for (i = 0; i < fftblock; i++) {
          ty1[k][i] = x[k][j][i+ii];
        }
      }

      cfftz(is, logd3, d3, ty1, ty2);

      for (k = 0; k < d3; k++) {
        for (i = 0; i < fftblock; i++) {
          xout[k][j][i+ii] = ty1[k][i];
        }
      }
    }
  }
}


void cffts3_inlined(int is, int d1, int d2, int d3,
	            dcomplex x[static const restrict d3][d2][d1+1],
	            dcomplex xout[static const restrict d3][d2][d1+1])
{
  int logd3;
  int i, j, k, ii;

  logd3 = ilog2(d3);

//  #pragma omp parallel for default(shared) private(i,j,k,ii)
#pragma scop
  for (j = 0; j < d2; j++) {
    for (ii = 0; ii <= d1 - FFTBLOCK_DEFAULT; ii += FFTBLOCK_DEFAULT) {
      for (k = 0; k < d3; k++) {
        for (i = 0; i < FFTBLOCK_DEFAULT; i++) {
          ty1[k][i] = x[k][j][i+ii];
        }
      }


      //------------------ cfftz_inlined
      int i0, j0, l;

      // Perform one variant of the Stockham FFT.
      for (l = 1; l <= logd3; l += 2) {
        int k, n1, li, lj, lk, ku, i11, i12, i21, i22;
        dcomplex u1, x11, x21;

        // Set initial parameters.
        n1 = d3 / 2;
        lk = 1 << (l - 1);
        li = 1 << (logd3 - l);
        lj = 2 * lk;
        ku = li;
        for (i0 = 0; i0 <= li - 1; i0++) {
          i11 = i0 * lk;
          i12 = i11 + n1;
          i21 = i0 * lj;
          i22 = i21 + lk;
          if (is >= 1) {
            u1 = u[ku+i0];
          } else {
            u1 = dconjg(u[ku+i0]);
          }

          // This loop is vectorizable.
          for (k = 0; k <= lk - 1; k++) {
            for (j0 = 0; j0 < FFTBLOCK_DEFAULT; j0++) {
              x11 = ty1[i11+k][j0];
              x21 = ty1[i12+k][j0];
              ty2[i21+k][j0] = dcmplx_add(x11, x21);
              ty2[i22+k][j0] = dcmplx_mul(u1, dcmplx_sub(x11, x21));
            }
          }
        }  

        // Copy Y to X in the last iteration,
        if (l == logd3) {
          for (j0 = 0; j0 < d3; j0++) {
            for (i0 = 0; i0 < FFTBLOCK_DEFAULT; i0++) {
	      ty1[j0][i0] = ty2[j0][i0];
            }
          }
        } else {
          // Set initial parameters.
          n1 = d3 / 2;
          lk = 1 << (l);
          li = 1 << (logd3 - l - 1);
          lj = 2 * lk;
          ku = li;
          for (i0 = 0; i0 <= li - 1; i0++) {
            i11 = i0 * lk;
            i12 = i11 + n1;
            i21 = i0 * lj;
            i22 = i21 + lk;
            if (is >= 1) {
              u1 = u[ku+i0];
            } else {
              u1 = dconjg(u[ku+i0]);
            }

            // This loop is vectorizable.
            for (k = 0; k <= lk - 1; k++) {
              for (j0 = 0; j0 < FFTBLOCK_DEFAULT; j0++) {
	        x11 = ty2[i11+k][j0];
                x21 = ty2[i12+k][j0];
                ty1[i21+k][j0] = dcmplx_add(x11, x21);
                ty1[i22+k][j0] = dcmplx_mul(u1, dcmplx_sub(x11, x21));
              }
            }
          }
        }
      }
      // -------------- end of cfftz_inlined


      for (k = 0; k < d3; k++) {
        for (i = 0; i < FFTBLOCK_DEFAULT; i++) {
          xout[k][j][i+ii] = ty1[k][i];
        }
      }
    }
  }
#pragma endscop  
}


//---------------------------------------------------------------------
// Computes NY N-point complex-to-complex FFTs of X using an algorithm due
// to Swarztrauber.  X is both the input and the output array, while Y is a 
// scratch array.  It is assumed that N = 2^M.  Before calling CFFTZ to 
// perform FFTs, the array U must be initialized by calling CFFTZ with IS 
// set to 0 and M set to MX, where MX is the maximum value of M for any 
// subsequent call.
//---------------------------------------------------------------------
void cfftz(int is, int m, int n, 
                  dcomplex x[n][fftblockpad], dcomplex y[n][fftblockpad])
{
  int i, j, l;

  //---------------------------------------------------------------------
  // Perform one variant of the Stockham FFT.
  //---------------------------------------------------------------------
  for (l = 1; l <= m; l += 2) {
    fftz2(is, l, m, n, fftblock, fftblockpad, u, x, y);
    if (l == m) {
      //-----------------------------------------------------------------
      // Copy Y to X.
      //-----------------------------------------------------------------
      for (j = 0; j < n; j++) {
        for (i = 0; i < fftblock; i++) {
          x[j][i] = y[j][i];
        }
      }
      return;
    }
    fftz2(is, l + 1, m, n, fftblock, fftblockpad, u, y, x);
  }
}


void cfftz_inlined(int is, int m, int n, 
                   dcomplex x[n][fftblockpad], dcomplex y[n][fftblockpad])
{
  int i, j, l;

  //---------------------------------------------------------------------
  // Perform one variant of the Stockham FFT.
  //---------------------------------------------------------------------
  for (l = 1; l <= m; l += 2) {
    int k, n1, li, lj, lk, ku, i11, i12, i21, i22;
    dcomplex u1, x11, x21;

    //---------------------------------------------------------------------
    // Set initial parameters.
    //---------------------------------------------------------------------
    n1 = n / 2;
    lk = 1 << (l - 1);
    li = 1 << (m - l);
    lj = 2 * lk;
    ku = li;

    for (i = 0; i <= li - 1; i++) {
      i11 = i * lk;
      i12 = i11 + n1;
      i21 = i * lj;
      i22 = i21 + lk;
      if (is >= 1) {
        u1 = u[ku+i];
      } else {
        u1 = dconjg(u[ku+i]);
      }

      //---------------------------------------------------------------------
      // This loop is vectorizable.
      //---------------------------------------------------------------------
      for (k = 0; k <= lk - 1; k++) {
        for (j = 0; j < FFTBLOCK_DEFAULT; j++) {
          x11 = x[i11+k][j];
          x21 = x[i12+k][j];
          y[i21+k][j] = dcmplx_add(x11, x21);
          y[i22+k][j] = dcmplx_mul(u1, dcmplx_sub(x11, x21));
        }
      }
    }

    if (l == m) {
      //-----------------------------------------------------------------
      // Copy Y to X.
      //-----------------------------------------------------------------
      for (j = 0; j < n; j++) {
        for (i = 0; i < FFTBLOCK_DEFAULT; i++) {
          x[j][i] = y[j][i];
        }
      }
      return;
    }

    //---------------------------------------------------------------------
    // Set initial parameters.
    //---------------------------------------------------------------------
    n1 = n / 2;
    lk = 1 << (l);
    li = 1 << (m - l - 1);
    lj = 2 * lk;
    ku = li;

    for (i = 0; i <= li - 1; i++) {
      i11 = i * lk;
      i12 = i11 + n1;
      i21 = i * lj;
      i22 = i21 + lk;
      if (is >= 1) {
        u1 = u[ku+i];
      } else {
        u1 = dconjg(u[ku+i]);
      }

      //---------------------------------------------------------------------
      // This loop is vectorizable.
      //---------------------------------------------------------------------
      for (k = 0; k <= lk - 1; k++) {
        for (j = 0; j < FFTBLOCK_DEFAULT; j++) {
          x11 = y[i11+k][j];
          x21 = y[i12+k][j];
          x[i21+k][j] = dcmplx_add(x11, x21);
          x[i22+k][j] = dcmplx_mul(u1, dcmplx_sub(x11, x21));
        }
      }
    }
  }
}


//---------------------------------------------------------------------
// Performs the L-th iteration of the second variant of the Stockham FFT.
//---------------------------------------------------------------------
void fftz2(int is, int l, int m, int n, int ny, int ny1, 
                  dcomplex u[n], dcomplex x[n][ny1], dcomplex y[n][ny1])
{
  int k, n1, li, lj, lk, ku, i, j, i11, i12, i21, i22;
  dcomplex u1, x11, x21;

  //---------------------------------------------------------------------
  // Set initial parameters.
  //---------------------------------------------------------------------
  n1 = n / 2;
  lk = 1 << (l - 1);
  li = 1 << (m - l);
  lj = 2 * lk;
  ku = li;

  for (i = 0; i <= li - 1; i++) {
    i11 = i * lk;
    i12 = i11 + n1;
    i21 = i * lj;
    i22 = i21 + lk;
    if (is >= 1) {
      u1 = u[ku+i];
    } else {
      u1 = dconjg(u[ku+i]);
    }

    //---------------------------------------------------------------------
    // This loop is vectorizable.
    //---------------------------------------------------------------------
    for (k = 0; k <= lk - 1; k++) {
      for (j = 0; j < ny; j++) {
        x11 = x[i11+k][j];
        x21 = x[i12+k][j];
        y[i21+k][j] = dcmplx_add(x11, x21);
        y[i22+k][j] = dcmplx_mul(u1, dcmplx_sub(x11, x21));
      }
    }
  }
}

