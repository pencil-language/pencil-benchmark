//-------------------------------------------------------------------------//
//                                                                         //
//  This benchmark is an OpenMP C version of the NPB BT code. This OpenMP  //
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

#include "header.h"
#include "work_lhs.h"
#include "timers.h"

//---------------------------------------------------------------------
// 
// Performs line solves in X direction by first factoring
// the block-tridiagonal matrix into an upper triangular matrix, 
// and then performing back substitution to solve for the unknow
// vectors of each line.  
// 
// Make sure we treat elements zero to cell_size in the direction
// of the sweep.
// 
//---------------------------------------------------------------------
void x_solve()
{
  int i, j, k, m, n, isize;
  isize = PROBLEM_SIZE-1;
  double pivot3, coeff3;
  double pivot2, coeff2;
  double pivot1, coeff1;

  //---------------------------------------------------------------------
  // This function computes the left hand side in the xi-direction
  //---------------------------------------------------------------------

//threadprivate() variables simplified (no duplicates)
//work_lhs.h: #pragma omp threadprivate(fjac,njac,lhs,tmp1,tmp2,tmp3)
//header.h:#pragma omp threadprivate(cuf,q,ue,buf)
// fjac[][][] and njac[][][] need to be privatized with private() to enable
// parallelization.  The first i loop writes in the two arrays and then these
// arrays are read later.
// The arrays u[][][], square[][][], rho_i[][][] and qs[][][] are live-in arrays.  They are
// read in the first i loop.
// lhs[][][] needs also privatization with openmp private(), it is accessed a
// write access in a loop and the read in a following loop (these two loops are
// a part of the same SCC).  lhs[][][] is used all over the program as a
// temporary scalar (do calculations and store them in lhs[][][] then in the
// following loop use these calculations).  A full fusion of i loops can enable
// contraction to reduce the dimensions of lhs[][][][] from 4D to 3D, because
// lhs[][][][] is always used as lhs[i][][][].  Each i loop write in
// lhs[i][*][*][*] and the following loops read from lhs[i][*][*][*], so if we
// fuse the i loops we will not need the i dimension since the fused loop will
// write in lhs[*][*][*] and directly read the lhs[*][*][*] in taht iteration.
// The code uses i loops because originally the loops were a part of different
// functions (modularity), now that inlining is applied, we can fuse the loops
// and contract lhs[][][][] into lhs[][][].  I didn't verify for the rest of
// the arrays and dependences whether they prevent loop fusion or not.
// rhs[k][][][] does not need privatization, it is the array that holds the result
// of the K loop.  Different k iterations write in different places in the
// array.


//#pragma omp parallel for default(shared) shared(isize) private(i,j,k,m,n)
#pragma scop
  for (k = 1; k <= PROBLEM_SIZE-2; k++) {
    for (j = 1; j <= PROBLEM_SIZE-2; j++) {
      for (i = 0; i <= isize; i++) {
        tmp1 = rho_i[k][j][i];
        tmp2 = tmp1 * tmp1;
        tmp3 = tmp1 * tmp2;
        
	fjac[i][0][0] = 0.0;
        fjac[i][1][0] = 1.0;
        fjac[i][2][0] = 0.0;
        fjac[i][3][0] = 0.0;
        fjac[i][4][0] = 0.0;

        fjac[i][0][1] = -(u[k][j][i][1] * tmp2 * u[k][j][i][1])
          + c2 * qs[k][j][i];
        fjac[i][1][1] = ( 2.0 - c2 ) * ( u[k][j][i][1] / u[k][j][i][0] );
        fjac[i][2][1] = - c2 * ( u[k][j][i][2] * tmp1 );
        fjac[i][3][1] = - c2 * ( u[k][j][i][3] * tmp1 );
        fjac[i][4][1] = c2;

        fjac[i][0][2] = - ( u[k][j][i][1]*u[k][j][i][2] ) * tmp2;
        fjac[i][1][2] = u[k][j][i][2] * tmp1;
        fjac[i][2][2] = u[k][j][i][1] * tmp1;
        fjac[i][3][2] = 0.0;
        fjac[i][4][2] = 0.0;

        fjac[i][0][3] = - ( u[k][j][i][1]*u[k][j][i][3] ) * tmp2;
        fjac[i][1][3] = u[k][j][i][3] * tmp1;
        fjac[i][2][3] = 0.0;
        fjac[i][3][3] = u[k][j][i][1] * tmp1;
        fjac[i][4][3] = 0.0;

        fjac[i][0][4] = ( c2 * 2.0 * square[k][j][i] - c1 * u[k][j][i][4] )
          * ( u[k][j][i][1] * tmp2 );
        fjac[i][1][4] = c1 *  u[k][j][i][4] * tmp1 
          - c2 * ( u[k][j][i][1]*u[k][j][i][1] * tmp2 + qs[k][j][i] );
        fjac[i][2][4] = - c2 * ( u[k][j][i][2]*u[k][j][i][1] ) * tmp2;
        fjac[i][3][4] = - c2 * ( u[k][j][i][3]*u[k][j][i][1] ) * tmp2;
        fjac[i][4][4] = c1 * ( u[k][j][i][1] * tmp1 );

        njac[i][0][0] = 0.0;
        njac[i][1][0] = 0.0;
        njac[i][2][0] = 0.0;
        njac[i][3][0] = 0.0;
        njac[i][4][0] = 0.0;

        njac[i][0][1] = - con43 * c3c4 * tmp2 * u[k][j][i][1];
        njac[i][1][1] =   con43 * c3c4 * tmp1;
        njac[i][2][1] =   0.0;
        njac[i][3][1] =   0.0;
        njac[i][4][1] =   0.0;

        njac[i][0][2] = - c3c4 * tmp2 * u[k][j][i][2];
        njac[i][1][2] =   0.0;
        njac[i][2][2] =   c3c4 * tmp1;
        njac[i][3][2] =   0.0;
        njac[i][4][2] =   0.0;

        njac[i][0][3] = - c3c4 * tmp2 * u[k][j][i][3];
        njac[i][1][3] =   0.0;
        njac[i][2][3] =   0.0;
        njac[i][3][3] =   c3c4 * tmp1;
        njac[i][4][3] =   0.0;

        njac[i][0][4] = - ( con43 * c3c4
            - c1345 ) * tmp3 * (u[k][j][i][1]*u[k][j][i][1])
          - ( c3c4 - c1345 ) * tmp3 * (u[k][j][i][2]*u[k][j][i][2])
          - ( c3c4 - c1345 ) * tmp3 * (u[k][j][i][3]*u[k][j][i][3])
          - c1345 * tmp2 * u[k][j][i][4];

        njac[i][1][4] = ( con43 * c3c4
            - c1345 ) * tmp2 * u[k][j][i][1];
        njac[i][2][4] = ( c3c4 - c1345 ) * tmp2 * u[k][j][i][2];
        njac[i][3][4] = ( c3c4 - c1345 ) * tmp2 * u[k][j][i][3];
        njac[i][4][4] = ( c1345 ) * tmp1;
      }
      // now jacobians set, so form left hand side in x direction
      
      //---------------------------------------------------------------------
      // lhsinit(lhs, isize);
      // void lhsinit(double lhs[][3][5][5], int ni)
      //---------------------------------------------------------------------
      for (n = 0; n < 5; n++) {
        for (m = 0; m < 5; m++) {
	      lhs[0][0][n][m] = 0.0;
	      lhs[0][1][n][m] = 0.0;
	      lhs[0][2][n][m] = 0.0;
        }
        lhs[0][1][n][n] = 1.0;
      }
      for (n = 0; n < 5; n++) {
        for (m = 0; m < 5; m++) {
	      lhs[isize][0][n][m] = 0.0;
	      lhs[isize][1][n][m] = 0.0;
	      lhs[isize][2][n][m] = 0.0;
        }
        lhs[isize][1][n][n] = 1.0;
      }

      for (i = 1; i <= isize-1; i++) {
        tmp1 = dt * tx1;
        tmp2 = dt * tx2;

        lhs[i][AA][0][0] = - tmp2 * fjac[i-1][0][0]
          - tmp1 * njac[i-1][0][0]
          - tmp1 * dx1; 
        lhs[i][AA][1][0] = - tmp2 * fjac[i-1][1][0]
          - tmp1 * njac[i-1][1][0];
        lhs[i][AA][2][0] = - tmp2 * fjac[i-1][2][0]
          - tmp1 * njac[i-1][2][0];
        lhs[i][AA][3][0] = - tmp2 * fjac[i-1][3][0]
          - tmp1 * njac[i-1][3][0];
        lhs[i][AA][4][0] = - tmp2 * fjac[i-1][4][0]
          - tmp1 * njac[i-1][4][0];

        lhs[i][AA][0][1] = - tmp2 * fjac[i-1][0][1]
          - tmp1 * njac[i-1][0][1];
        lhs[i][AA][1][1] = - tmp2 * fjac[i-1][1][1]
          - tmp1 * njac[i-1][1][1]
          - tmp1 * dx2;
        lhs[i][AA][2][1] = - tmp2 * fjac[i-1][2][1]
          - tmp1 * njac[i-1][2][1];
        lhs[i][AA][3][1] = - tmp2 * fjac[i-1][3][1]
          - tmp1 * njac[i-1][3][1];
        lhs[i][AA][4][1] = - tmp2 * fjac[i-1][4][1]
          - tmp1 * njac[i-1][4][1];

        lhs[i][AA][0][2] = - tmp2 * fjac[i-1][0][2]
          - tmp1 * njac[i-1][0][2];
        lhs[i][AA][1][2] = - tmp2 * fjac[i-1][1][2]
          - tmp1 * njac[i-1][1][2];
        lhs[i][AA][2][2] = - tmp2 * fjac[i-1][2][2]
          - tmp1 * njac[i-1][2][2]
          - tmp1 * dx3;
        lhs[i][AA][3][2] = - tmp2 * fjac[i-1][3][2]
          - tmp1 * njac[i-1][3][2];
        lhs[i][AA][4][2] = - tmp2 * fjac[i-1][4][2]
          - tmp1 * njac[i-1][4][2];

        lhs[i][AA][0][3] = - tmp2 * fjac[i-1][0][3]
          - tmp1 * njac[i-1][0][3];
        lhs[i][AA][1][3] = - tmp2 * fjac[i-1][1][3]
          - tmp1 * njac[i-1][1][3];
        lhs[i][AA][2][3] = - tmp2 * fjac[i-1][2][3]
          - tmp1 * njac[i-1][2][3];
        lhs[i][AA][3][3] = - tmp2 * fjac[i-1][3][3]
          - tmp1 * njac[i-1][3][3]
          - tmp1 * dx4;
        lhs[i][AA][4][3] = - tmp2 * fjac[i-1][4][3]
          - tmp1 * njac[i-1][4][3];

        lhs[i][AA][0][4] = - tmp2 * fjac[i-1][0][4]
          - tmp1 * njac[i-1][0][4];
        lhs[i][AA][1][4] = - tmp2 * fjac[i-1][1][4]
          - tmp1 * njac[i-1][1][4];
        lhs[i][AA][2][4] = - tmp2 * fjac[i-1][2][4]
          - tmp1 * njac[i-1][2][4];
        lhs[i][AA][3][4] = - tmp2 * fjac[i-1][3][4]
          - tmp1 * njac[i-1][3][4];
        lhs[i][AA][4][4] = - tmp2 * fjac[i-1][4][4]
          - tmp1 * njac[i-1][4][4]
          - tmp1 * dx5;

        lhs[i][BB][0][0] = 1.0
          + tmp1 * 2.0 * njac[i][0][0]
          + tmp1 * 2.0 * dx1;
        lhs[i][BB][1][0] = tmp1 * 2.0 * njac[i][1][0];
        lhs[i][BB][2][0] = tmp1 * 2.0 * njac[i][2][0];
        lhs[i][BB][3][0] = tmp1 * 2.0 * njac[i][3][0];
        lhs[i][BB][4][0] = tmp1 * 2.0 * njac[i][4][0];

        lhs[i][BB][0][1] = tmp1 * 2.0 * njac[i][0][1];
        lhs[i][BB][1][1] = 1.0
          + tmp1 * 2.0 * njac[i][1][1]
          + tmp1 * 2.0 * dx2;
        lhs[i][BB][2][1] = tmp1 * 2.0 * njac[i][2][1];
        lhs[i][BB][3][1] = tmp1 * 2.0 * njac[i][3][1];
        lhs[i][BB][4][1] = tmp1 * 2.0 * njac[i][4][1];

        lhs[i][BB][0][2] = tmp1 * 2.0 * njac[i][0][2];
        lhs[i][BB][1][2] = tmp1 * 2.0 * njac[i][1][2];
        lhs[i][BB][2][2] = 1.0
          + tmp1 * 2.0 * njac[i][2][2]
          + tmp1 * 2.0 * dx3;
        lhs[i][BB][3][2] = tmp1 * 2.0 * njac[i][3][2];
        lhs[i][BB][4][2] = tmp1 * 2.0 * njac[i][4][2];

        lhs[i][BB][0][3] = tmp1 * 2.0 * njac[i][0][3];
        lhs[i][BB][1][3] = tmp1 * 2.0 * njac[i][1][3];
        lhs[i][BB][2][3] = tmp1 * 2.0 * njac[i][2][3];
        lhs[i][BB][3][3] = 1.0
          + tmp1 * 2.0 * njac[i][3][3]
          + tmp1 * 2.0 * dx4;
        lhs[i][BB][4][3] = tmp1 * 2.0 * njac[i][4][3];

        lhs[i][BB][0][4] = tmp1 * 2.0 * njac[i][0][4];
        lhs[i][BB][1][4] = tmp1 * 2.0 * njac[i][1][4];
        lhs[i][BB][2][4] = tmp1 * 2.0 * njac[i][2][4];
        lhs[i][BB][3][4] = tmp1 * 2.0 * njac[i][3][4];
        lhs[i][BB][4][4] = 1.0
          + tmp1 * 2.0 * njac[i][4][4]
          + tmp1 * 2.0 * dx5;

        lhs[i][CC][0][0] =  tmp2 * fjac[i+1][0][0]
          - tmp1 * njac[i+1][0][0]
          - tmp1 * dx1;
        lhs[i][CC][1][0] =  tmp2 * fjac[i+1][1][0]
          - tmp1 * njac[i+1][1][0];
        lhs[i][CC][2][0] =  tmp2 * fjac[i+1][2][0]
          - tmp1 * njac[i+1][2][0];
        lhs[i][CC][3][0] =  tmp2 * fjac[i+1][3][0]
          - tmp1 * njac[i+1][3][0];
        lhs[i][CC][4][0] =  tmp2 * fjac[i+1][4][0]
          - tmp1 * njac[i+1][4][0];

        lhs[i][CC][0][1] =  tmp2 * fjac[i+1][0][1]
          - tmp1 * njac[i+1][0][1];
        lhs[i][CC][1][1] =  tmp2 * fjac[i+1][1][1]
          - tmp1 * njac[i+1][1][1]
          - tmp1 * dx2;
        lhs[i][CC][2][1] =  tmp2 * fjac[i+1][2][1]
          - tmp1 * njac[i+1][2][1];
        lhs[i][CC][3][1] =  tmp2 * fjac[i+1][3][1]
          - tmp1 * njac[i+1][3][1];
        lhs[i][CC][4][1] =  tmp2 * fjac[i+1][4][1]
          - tmp1 * njac[i+1][4][1];

        lhs[i][CC][0][2] =  tmp2 * fjac[i+1][0][2]
          - tmp1 * njac[i+1][0][2];
        lhs[i][CC][1][2] =  tmp2 * fjac[i+1][1][2]
          - tmp1 * njac[i+1][1][2];
        lhs[i][CC][2][2] =  tmp2 * fjac[i+1][2][2]
          - tmp1 * njac[i+1][2][2]
          - tmp1 * dx3;
        lhs[i][CC][3][2] =  tmp2 * fjac[i+1][3][2]
          - tmp1 * njac[i+1][3][2];
        lhs[i][CC][4][2] =  tmp2 * fjac[i+1][4][2]
          - tmp1 * njac[i+1][4][2];

        lhs[i][CC][0][3] =  tmp2 * fjac[i+1][0][3]
          - tmp1 * njac[i+1][0][3];
        lhs[i][CC][1][3] =  tmp2 * fjac[i+1][1][3]
          - tmp1 * njac[i+1][1][3];
        lhs[i][CC][2][3] =  tmp2 * fjac[i+1][2][3]
          - tmp1 * njac[i+1][2][3];
        lhs[i][CC][3][3] =  tmp2 * fjac[i+1][3][3]
          - tmp1 * njac[i+1][3][3]
          - tmp1 * dx4;
        lhs[i][CC][4][3] =  tmp2 * fjac[i+1][4][3]
          - tmp1 * njac[i+1][4][3];

        lhs[i][CC][0][4] =  tmp2 * fjac[i+1][0][4]
          - tmp1 * njac[i+1][0][4];
        lhs[i][CC][1][4] =  tmp2 * fjac[i+1][1][4]
          - tmp1 * njac[i+1][1][4];
        lhs[i][CC][2][4] =  tmp2 * fjac[i+1][2][4]
          - tmp1 * njac[i+1][2][4];
        lhs[i][CC][3][4] =  tmp2 * fjac[i+1][3][4]
          - tmp1 * njac[i+1][3][4];
        lhs[i][CC][4][4] =  tmp2 * fjac[i+1][4][4]
          - tmp1 * njac[i+1][4][4]
          - tmp1 * dx5;
      }

      //---------------------------------------------------------------------
      // performs guaussian elimination on this cell.
      // 
      // assumes that unpacking routines for non-first cells 
      // preload C' and rhs' from previous cell.
      // 
      // assumed send happens outside this routine, but that
      // c'(IMAX) and rhs'(IMAX) will be sent to next cell
      //
      // outer most do loops - sweeping in i direction
      //
      // multiply c[k][j][0] by b_inverse and copy back to c
      // multiply rhs(0) by b_inverse(0) and copy to rhs
      //---------------------------------------------------------------------
      // binvcrhs( lhs[0][BB], lhs[0][CC], rhs[k][j][0] );
      // void binvcrhs(double lhs[5][5], double c[5][5], double r[5])
      {
	  pivot1 = 1.00/lhs[0][BB][0][0];
	  lhs[0][BB][1][0] = lhs[0][BB][1][0]*pivot1;
	  lhs[0][BB][2][0] = lhs[0][BB][2][0]*pivot1;
	  lhs[0][BB][3][0] = lhs[0][BB][3][0]*pivot1;
	  lhs[0][BB][4][0] = lhs[0][BB][4][0]*pivot1;
	  lhs[0][CC][0][0] = lhs[0][CC][0][0]*pivot1;
	  lhs[0][CC][1][0] = lhs[0][CC][1][0]*pivot1;
	  lhs[0][CC][2][0] = lhs[0][CC][2][0]*pivot1;
	  lhs[0][CC][3][0] = lhs[0][CC][3][0]*pivot1;
	  lhs[0][CC][4][0] = lhs[0][CC][4][0]*pivot1;
	  rhs[k][j][0][0]   = rhs[k][j][0][0]  *pivot1;

	  coeff1 = lhs[0][BB][0][1];
	  lhs[0][BB][1][1]= lhs[0][BB][1][1] - coeff1*lhs[0][BB][1][0];
	  lhs[0][BB][2][1]= lhs[0][BB][2][1] - coeff1*lhs[0][BB][2][0];
	  lhs[0][BB][3][1]= lhs[0][BB][3][1] - coeff1*lhs[0][BB][3][0];
	  lhs[0][BB][4][1]= lhs[0][BB][4][1] - coeff1*lhs[0][BB][4][0];
	  lhs[0][CC][0][1] = lhs[0][CC][0][1] - coeff1*lhs[0][CC][0][0];
	  lhs[0][CC][1][1] = lhs[0][CC][1][1] - coeff1*lhs[0][CC][1][0];
	  lhs[0][CC][2][1] = lhs[0][CC][2][1] - coeff1*lhs[0][CC][2][0];
	  lhs[0][CC][3][1] = lhs[0][CC][3][1] - coeff1*lhs[0][CC][3][0];
	  lhs[0][CC][4][1] = lhs[0][CC][4][1] - coeff1*lhs[0][CC][4][0];
	  rhs[k][j][0][1]   = rhs[k][j][0][1]   - coeff1*rhs[k][j][0][0];

	  coeff1 = lhs[0][BB][0][2];
	  lhs[0][BB][1][2]= lhs[0][BB][1][2] - coeff1*lhs[0][BB][1][0];
	  lhs[0][BB][2][2]= lhs[0][BB][2][2] - coeff1*lhs[0][BB][2][0];
	  lhs[0][BB][3][2]= lhs[0][BB][3][2] - coeff1*lhs[0][BB][3][0];
	  lhs[0][BB][4][2]= lhs[0][BB][4][2] - coeff1*lhs[0][BB][4][0];
	  lhs[0][CC][0][2] = lhs[0][CC][0][2] - coeff1*lhs[0][CC][0][0];
	  lhs[0][CC][1][2] = lhs[0][CC][1][2] - coeff1*lhs[0][CC][1][0];
	  lhs[0][CC][2][2] = lhs[0][CC][2][2] - coeff1*lhs[0][CC][2][0];
	  lhs[0][CC][3][2] = lhs[0][CC][3][2] - coeff1*lhs[0][CC][3][0];
	  lhs[0][CC][4][2] = lhs[0][CC][4][2] - coeff1*lhs[0][CC][4][0];
	  rhs[k][j][0][2]   = rhs[k][j][0][2]   - coeff1*rhs[k][j][0][0];

	  coeff1 = lhs[0][BB][0][3];
	  lhs[0][BB][1][3]= lhs[0][BB][1][3] - coeff1*lhs[0][BB][1][0];
	  lhs[0][BB][2][3]= lhs[0][BB][2][3] - coeff1*lhs[0][BB][2][0];
	  lhs[0][BB][3][3]= lhs[0][BB][3][3] - coeff1*lhs[0][BB][3][0];
	  lhs[0][BB][4][3]= lhs[0][BB][4][3] - coeff1*lhs[0][BB][4][0];
	  lhs[0][CC][0][3] = lhs[0][CC][0][3] - coeff1*lhs[0][CC][0][0];
	  lhs[0][CC][1][3] = lhs[0][CC][1][3] - coeff1*lhs[0][CC][1][0];
	  lhs[0][CC][2][3] = lhs[0][CC][2][3] - coeff1*lhs[0][CC][2][0];
	  lhs[0][CC][3][3] = lhs[0][CC][3][3] - coeff1*lhs[0][CC][3][0];
	  lhs[0][CC][4][3] = lhs[0][CC][4][3] - coeff1*lhs[0][CC][4][0];
	  rhs[k][j][0][3]   = rhs[k][j][0][3]   - coeff1*rhs[k][j][0][0];

	  coeff1 = lhs[0][BB][0][4];
	  lhs[0][BB][1][4]= lhs[0][BB][1][4] - coeff1*lhs[0][BB][1][0];
	  lhs[0][BB][2][4]= lhs[0][BB][2][4] - coeff1*lhs[0][BB][2][0];
	  lhs[0][BB][3][4]= lhs[0][BB][3][4] - coeff1*lhs[0][BB][3][0];
	  lhs[0][BB][4][4]= lhs[0][BB][4][4] - coeff1*lhs[0][BB][4][0];
	  lhs[0][CC][0][4] = lhs[0][CC][0][4] - coeff1*lhs[0][CC][0][0];
	  lhs[0][CC][1][4] = lhs[0][CC][1][4] - coeff1*lhs[0][CC][1][0];
	  lhs[0][CC][2][4] = lhs[0][CC][2][4] - coeff1*lhs[0][CC][2][0];
	  lhs[0][CC][3][4] = lhs[0][CC][3][4] - coeff1*lhs[0][CC][3][0];
	  lhs[0][CC][4][4] = lhs[0][CC][4][4] - coeff1*lhs[0][CC][4][0];
	  rhs[k][j][0][4]   = rhs[k][j][0][4]   - coeff1*rhs[k][j][0][0];


	  pivot1 = 1.00/lhs[0][BB][1][1];
	  lhs[0][BB][2][1] = lhs[0][BB][2][1]*pivot1;
	  lhs[0][BB][3][1] = lhs[0][BB][3][1]*pivot1;
	  lhs[0][BB][4][1] = lhs[0][BB][4][1]*pivot1;
	  lhs[0][CC][0][1] = lhs[0][CC][0][1]*pivot1;
	  lhs[0][CC][1][1] = lhs[0][CC][1][1]*pivot1;
	  lhs[0][CC][2][1] = lhs[0][CC][2][1]*pivot1;
	  lhs[0][CC][3][1] = lhs[0][CC][3][1]*pivot1;
	  lhs[0][CC][4][1] = lhs[0][CC][4][1]*pivot1;
	  rhs[k][j][0][1]   = rhs[k][j][0][1]  *pivot1;

	  coeff1 = lhs[0][BB][1][0];
	  lhs[0][BB][2][0]= lhs[0][BB][2][0] - coeff1*lhs[0][BB][2][1];
	  lhs[0][BB][3][0]= lhs[0][BB][3][0] - coeff1*lhs[0][BB][3][1];
	  lhs[0][BB][4][0]= lhs[0][BB][4][0] - coeff1*lhs[0][BB][4][1];
	  lhs[0][CC][0][0] = lhs[0][CC][0][0] - coeff1*lhs[0][CC][0][1];
	  lhs[0][CC][1][0] = lhs[0][CC][1][0] - coeff1*lhs[0][CC][1][1];
	  lhs[0][CC][2][0] = lhs[0][CC][2][0] - coeff1*lhs[0][CC][2][1];
	  lhs[0][CC][3][0] = lhs[0][CC][3][0] - coeff1*lhs[0][CC][3][1];
	  lhs[0][CC][4][0] = lhs[0][CC][4][0] - coeff1*lhs[0][CC][4][1];
	  rhs[k][j][0][0]   = rhs[k][j][0][0]   - coeff1*rhs[k][j][0][1];

	  coeff1 = lhs[0][BB][1][2];
	  lhs[0][BB][2][2]= lhs[0][BB][2][2] - coeff1*lhs[0][BB][2][1];
	  lhs[0][BB][3][2]= lhs[0][BB][3][2] - coeff1*lhs[0][BB][3][1];
	  lhs[0][BB][4][2]= lhs[0][BB][4][2] - coeff1*lhs[0][BB][4][1];
	  lhs[0][CC][0][2] = lhs[0][CC][0][2] - coeff1*lhs[0][CC][0][1];
	  lhs[0][CC][1][2] = lhs[0][CC][1][2] - coeff1*lhs[0][CC][1][1];
	  lhs[0][CC][2][2] = lhs[0][CC][2][2] - coeff1*lhs[0][CC][2][1];
	  lhs[0][CC][3][2] = lhs[0][CC][3][2] - coeff1*lhs[0][CC][3][1];
	  lhs[0][CC][4][2] = lhs[0][CC][4][2] - coeff1*lhs[0][CC][4][1];
	  rhs[k][j][0][2]   = rhs[k][j][0][2]   - coeff1*rhs[k][j][0][1];

	  coeff1 = lhs[0][BB][1][3];
	  lhs[0][BB][2][3]= lhs[0][BB][2][3] - coeff1*lhs[0][BB][2][1];
	  lhs[0][BB][3][3]= lhs[0][BB][3][3] - coeff1*lhs[0][BB][3][1];
	  lhs[0][BB][4][3]= lhs[0][BB][4][3] - coeff1*lhs[0][BB][4][1];
	  lhs[0][CC][0][3] = lhs[0][CC][0][3] - coeff1*lhs[0][CC][0][1];
	  lhs[0][CC][1][3] = lhs[0][CC][1][3] - coeff1*lhs[0][CC][1][1];
	  lhs[0][CC][2][3] = lhs[0][CC][2][3] - coeff1*lhs[0][CC][2][1];
	  lhs[0][CC][3][3] = lhs[0][CC][3][3] - coeff1*lhs[0][CC][3][1];
	  lhs[0][CC][4][3] = lhs[0][CC][4][3] - coeff1*lhs[0][CC][4][1];
	  rhs[k][j][0][3]   = rhs[k][j][0][3]   - coeff1*rhs[k][j][0][1];

	  coeff1 = lhs[0][BB][1][4];
	  lhs[0][BB][2][4]= lhs[0][BB][2][4] - coeff1*lhs[0][BB][2][1];
	  lhs[0][BB][3][4]= lhs[0][BB][3][4] - coeff1*lhs[0][BB][3][1];
	  lhs[0][BB][4][4]= lhs[0][BB][4][4] - coeff1*lhs[0][BB][4][1];
	  lhs[0][CC][0][4] = lhs[0][CC][0][4] - coeff1*lhs[0][CC][0][1];
	  lhs[0][CC][1][4] = lhs[0][CC][1][4] - coeff1*lhs[0][CC][1][1];
	  lhs[0][CC][2][4] = lhs[0][CC][2][4] - coeff1*lhs[0][CC][2][1];
	  lhs[0][CC][3][4] = lhs[0][CC][3][4] - coeff1*lhs[0][CC][3][1];
	  lhs[0][CC][4][4] = lhs[0][CC][4][4] - coeff1*lhs[0][CC][4][1];
	  rhs[k][j][0][4]   = rhs[k][j][0][4]   - coeff1*rhs[k][j][0][1];


	  pivot1 = 1.00/lhs[0][BB][2][2];
	  lhs[0][BB][3][2] = lhs[0][BB][3][2]*pivot1;
	  lhs[0][BB][4][2] = lhs[0][BB][4][2]*pivot1;
	  lhs[0][CC][0][2] = lhs[0][CC][0][2]*pivot1;
	  lhs[0][CC][1][2] = lhs[0][CC][1][2]*pivot1;
	  lhs[0][CC][2][2] = lhs[0][CC][2][2]*pivot1;
	  lhs[0][CC][3][2] = lhs[0][CC][3][2]*pivot1;
	  lhs[0][CC][4][2] = lhs[0][CC][4][2]*pivot1;
	  rhs[k][j][0][2]   = rhs[k][j][0][2]  *pivot1;

	  coeff1 = lhs[0][BB][2][0];
	  lhs[0][BB][3][0]= lhs[0][BB][3][0] - coeff1*lhs[0][BB][3][2];
	  lhs[0][BB][4][0]= lhs[0][BB][4][0] - coeff1*lhs[0][BB][4][2];
	  lhs[0][CC][0][0] = lhs[0][CC][0][0] - coeff1*lhs[0][CC][0][2];
	  lhs[0][CC][1][0] = lhs[0][CC][1][0] - coeff1*lhs[0][CC][1][2];
	  lhs[0][CC][2][0] = lhs[0][CC][2][0] - coeff1*lhs[0][CC][2][2];
	  lhs[0][CC][3][0] = lhs[0][CC][3][0] - coeff1*lhs[0][CC][3][2];
	  lhs[0][CC][4][0] = lhs[0][CC][4][0] - coeff1*lhs[0][CC][4][2];
	  rhs[k][j][0][0]   = rhs[k][j][0][0]   - coeff1*rhs[k][j][0][2];

	  coeff1 = lhs[0][BB][2][1];
	  lhs[0][BB][3][1]= lhs[0][BB][3][1] - coeff1*lhs[0][BB][3][2];
	  lhs[0][BB][4][1]= lhs[0][BB][4][1] - coeff1*lhs[0][BB][4][2];
	  lhs[0][CC][0][1] = lhs[0][CC][0][1] - coeff1*lhs[0][CC][0][2];
	  lhs[0][CC][1][1] = lhs[0][CC][1][1] - coeff1*lhs[0][CC][1][2];
	  lhs[0][CC][2][1] = lhs[0][CC][2][1] - coeff1*lhs[0][CC][2][2];
	  lhs[0][CC][3][1] = lhs[0][CC][3][1] - coeff1*lhs[0][CC][3][2];
	  lhs[0][CC][4][1] = lhs[0][CC][4][1] - coeff1*lhs[0][CC][4][2];
	  rhs[k][j][0][1]   = rhs[k][j][0][1]   - coeff1*rhs[k][j][0][2];

	  coeff1 = lhs[0][BB][2][3];
	  lhs[0][BB][3][3]= lhs[0][BB][3][3] - coeff1*lhs[0][BB][3][2];
	  lhs[0][BB][4][3]= lhs[0][BB][4][3] - coeff1*lhs[0][BB][4][2];
	  lhs[0][CC][0][3] = lhs[0][CC][0][3] - coeff1*lhs[0][CC][0][2];
	  lhs[0][CC][1][3] = lhs[0][CC][1][3] - coeff1*lhs[0][CC][1][2];
	  lhs[0][CC][2][3] = lhs[0][CC][2][3] - coeff1*lhs[0][CC][2][2];
	  lhs[0][CC][3][3] = lhs[0][CC][3][3] - coeff1*lhs[0][CC][3][2];
	  lhs[0][CC][4][3] = lhs[0][CC][4][3] - coeff1*lhs[0][CC][4][2];
	  rhs[k][j][0][3]   = rhs[k][j][0][3]   - coeff1*rhs[k][j][0][2];

	  coeff1 = lhs[0][BB][2][4];
	  lhs[0][BB][3][4]= lhs[0][BB][3][4] - coeff1*lhs[0][BB][3][2];
	  lhs[0][BB][4][4]= lhs[0][BB][4][4] - coeff1*lhs[0][BB][4][2];
	  lhs[0][CC][0][4] = lhs[0][CC][0][4] - coeff1*lhs[0][CC][0][2];
	  lhs[0][CC][1][4] = lhs[0][CC][1][4] - coeff1*lhs[0][CC][1][2];
	  lhs[0][CC][2][4] = lhs[0][CC][2][4] - coeff1*lhs[0][CC][2][2];
	  lhs[0][CC][3][4] = lhs[0][CC][3][4] - coeff1*lhs[0][CC][3][2];
	  lhs[0][CC][4][4] = lhs[0][CC][4][4] - coeff1*lhs[0][CC][4][2];
	  rhs[k][j][0][4]   = rhs[k][j][0][4]   - coeff1*rhs[k][j][0][2];


	  pivot1 = 1.00/lhs[0][BB][3][3];
	  lhs[0][BB][4][3] = lhs[0][BB][4][3]*pivot1;
	  lhs[0][CC][0][3] = lhs[0][CC][0][3]*pivot1;
	  lhs[0][CC][1][3] = lhs[0][CC][1][3]*pivot1;
	  lhs[0][CC][2][3] = lhs[0][CC][2][3]*pivot1;
	  lhs[0][CC][3][3] = lhs[0][CC][3][3]*pivot1;
	  lhs[0][CC][4][3] = lhs[0][CC][4][3]*pivot1;
	  rhs[k][j][0][3]   = rhs[k][j][0][3]  *pivot1;

	  coeff1 = lhs[0][BB][3][0];
	  lhs[0][BB][4][0]= lhs[0][BB][4][0] - coeff1*lhs[0][BB][4][3];
	  lhs[0][CC][0][0] = lhs[0][CC][0][0] - coeff1*lhs[0][CC][0][3];
	  lhs[0][CC][1][0] = lhs[0][CC][1][0] - coeff1*lhs[0][CC][1][3];
	  lhs[0][CC][2][0] = lhs[0][CC][2][0] - coeff1*lhs[0][CC][2][3];
	  lhs[0][CC][3][0] = lhs[0][CC][3][0] - coeff1*lhs[0][CC][3][3];
	  lhs[0][CC][4][0] = lhs[0][CC][4][0] - coeff1*lhs[0][CC][4][3];
	  rhs[k][j][0][0]   = rhs[k][j][0][0]   - coeff1*rhs[k][j][0][3];

	  coeff1 = lhs[0][BB][3][1];
	  lhs[0][BB][4][1]= lhs[0][BB][4][1] - coeff1*lhs[0][BB][4][3];
	  lhs[0][CC][0][1] = lhs[0][CC][0][1] - coeff1*lhs[0][CC][0][3];
	  lhs[0][CC][1][1] = lhs[0][CC][1][1] - coeff1*lhs[0][CC][1][3];
	  lhs[0][CC][2][1] = lhs[0][CC][2][1] - coeff1*lhs[0][CC][2][3];
	  lhs[0][CC][3][1] = lhs[0][CC][3][1] - coeff1*lhs[0][CC][3][3];
	  lhs[0][CC][4][1] = lhs[0][CC][4][1] - coeff1*lhs[0][CC][4][3];
	  rhs[k][j][0][1]   = rhs[k][j][0][1]   - coeff1*rhs[k][j][0][3];

	  coeff1 = lhs[0][BB][3][2];
	  lhs[0][BB][4][2]= lhs[0][BB][4][2] - coeff1*lhs[0][BB][4][3];
	  lhs[0][CC][0][2] = lhs[0][CC][0][2] - coeff1*lhs[0][CC][0][3];
	  lhs[0][CC][1][2] = lhs[0][CC][1][2] - coeff1*lhs[0][CC][1][3];
	  lhs[0][CC][2][2] = lhs[0][CC][2][2] - coeff1*lhs[0][CC][2][3];
	  lhs[0][CC][3][2] = lhs[0][CC][3][2] - coeff1*lhs[0][CC][3][3];
	  lhs[0][CC][4][2] = lhs[0][CC][4][2] - coeff1*lhs[0][CC][4][3];
	  rhs[k][j][0][2]   = rhs[k][j][0][2]   - coeff1*rhs[k][j][0][3];

	  coeff1 = lhs[0][BB][3][4];
	  lhs[0][BB][4][4]= lhs[0][BB][4][4] - coeff1*lhs[0][BB][4][3];
	  lhs[0][CC][0][4] = lhs[0][CC][0][4] - coeff1*lhs[0][CC][0][3];
	  lhs[0][CC][1][4] = lhs[0][CC][1][4] - coeff1*lhs[0][CC][1][3];
	  lhs[0][CC][2][4] = lhs[0][CC][2][4] - coeff1*lhs[0][CC][2][3];
	  lhs[0][CC][3][4] = lhs[0][CC][3][4] - coeff1*lhs[0][CC][3][3];
	  lhs[0][CC][4][4] = lhs[0][CC][4][4] - coeff1*lhs[0][CC][4][3];
	  rhs[k][j][0][4]   = rhs[k][j][0][4]   - coeff1*rhs[k][j][0][3];


	  pivot1 = 1.00/lhs[0][BB][4][4];
	  lhs[0][CC][0][4] = lhs[0][CC][0][4]*pivot1;
	  lhs[0][CC][1][4] = lhs[0][CC][1][4]*pivot1;
	  lhs[0][CC][2][4] = lhs[0][CC][2][4]*pivot1;
	  lhs[0][CC][3][4] = lhs[0][CC][3][4]*pivot1;
	  lhs[0][CC][4][4] = lhs[0][CC][4][4]*pivot1;
	  rhs[k][j][0][4]   = rhs[k][j][0][4]  *pivot1;

	  coeff1 = lhs[0][BB][4][0];
	  lhs[0][CC][0][0] = lhs[0][CC][0][0] - coeff1*lhs[0][CC][0][4];
	  lhs[0][CC][1][0] = lhs[0][CC][1][0] - coeff1*lhs[0][CC][1][4];
	  lhs[0][CC][2][0] = lhs[0][CC][2][0] - coeff1*lhs[0][CC][2][4];
	  lhs[0][CC][3][0] = lhs[0][CC][3][0] - coeff1*lhs[0][CC][3][4];
	  lhs[0][CC][4][0] = lhs[0][CC][4][0] - coeff1*lhs[0][CC][4][4];
	  rhs[k][j][0][0]   = rhs[k][j][0][0]   - coeff1*rhs[k][j][0][4];

	  coeff1 = lhs[0][BB][4][1];
	  lhs[0][CC][0][1] = lhs[0][CC][0][1] - coeff1*lhs[0][CC][0][4];
	  lhs[0][CC][1][1] = lhs[0][CC][1][1] - coeff1*lhs[0][CC][1][4];
	  lhs[0][CC][2][1] = lhs[0][CC][2][1] - coeff1*lhs[0][CC][2][4];
	  lhs[0][CC][3][1] = lhs[0][CC][3][1] - coeff1*lhs[0][CC][3][4];
	  lhs[0][CC][4][1] = lhs[0][CC][4][1] - coeff1*lhs[0][CC][4][4];
	  rhs[k][j][0][1]   = rhs[k][j][0][1]   - coeff1*rhs[k][j][0][4];

	  coeff1 = lhs[0][BB][4][2];
	  lhs[0][CC][0][2] = lhs[0][CC][0][2] - coeff1*lhs[0][CC][0][4];
	  lhs[0][CC][1][2] = lhs[0][CC][1][2] - coeff1*lhs[0][CC][1][4];
	  lhs[0][CC][2][2] = lhs[0][CC][2][2] - coeff1*lhs[0][CC][2][4];
	  lhs[0][CC][3][2] = lhs[0][CC][3][2] - coeff1*lhs[0][CC][3][4];
	  lhs[0][CC][4][2] = lhs[0][CC][4][2] - coeff1*lhs[0][CC][4][4];
	  rhs[k][j][0][2]   = rhs[k][j][0][2]   - coeff1*rhs[k][j][0][4];

	  coeff1 = lhs[0][BB][4][3];
	  lhs[0][CC][0][3] = lhs[0][CC][0][3] - coeff1*lhs[0][CC][0][4];
	  lhs[0][CC][1][3] = lhs[0][CC][1][3] - coeff1*lhs[0][CC][1][4];
	  lhs[0][CC][2][3] = lhs[0][CC][2][3] - coeff1*lhs[0][CC][2][4];
	  lhs[0][CC][3][3] = lhs[0][CC][3][3] - coeff1*lhs[0][CC][3][4];
	  lhs[0][CC][4][3] = lhs[0][CC][4][3] - coeff1*lhs[0][CC][4][4];
	  rhs[k][j][0][3]   = rhs[k][j][0][3]   - coeff1*rhs[k][j][0][4];
      } //END of binvcrhs( lhs[0][BB], lhs[0][CC], rhs[k][j][0] );

      // begin inner most do loop
      // do all the elements of the cell unless last 
      for (i = 1; i <= isize-1; i++) {
        //-------------------------------------------------------------------
        // rhs(i) = rhs(i) - A*rhs(i-1)
        //
        // matvec_sub(       lhs[i][AA],   rhs[k][j][i-1], rhs[k][j][i]);
        // void matvec_sub(double ablock[5][5], double avec[5], double bvec[5])
        //-------------------------------------------------------------------
	{
	  rhs[k][j][i][0] = rhs[k][j][i][0] - lhs[i][AA][0][0]*rhs[k][j][i-1][0]
			    - lhs[i][AA][1][0]*rhs[k][j][i-1][1]
			    - lhs[i][AA][2][0]*rhs[k][j][i-1][2]
			    - lhs[i][AA][3][0]*rhs[k][j][i-1][3]
			    - lhs[i][AA][4][0]*rhs[k][j][i-1][4];
	  rhs[k][j][i][1] = rhs[k][j][i][1] - lhs[i][AA][0][1]*rhs[k][j][i-1][0]
			    - lhs[i][AA][1][1]*rhs[k][j][i-1][1]
			    - lhs[i][AA][2][1]*rhs[k][j][i-1][2]
			    - lhs[i][AA][3][1]*rhs[k][j][i-1][3]
			    - lhs[i][AA][4][1]*rhs[k][j][i-1][4];
	  rhs[k][j][i][2] = rhs[k][j][i][2] - lhs[i][AA][0][2]*rhs[k][j][i-1][0]
			    - lhs[i][AA][1][2]*rhs[k][j][i-1][1]
			    - lhs[i][AA][2][2]*rhs[k][j][i-1][2]
			    - lhs[i][AA][3][2]*rhs[k][j][i-1][3]
			    - lhs[i][AA][4][2]*rhs[k][j][i-1][4];
	  rhs[k][j][i][3] = rhs[k][j][i][3] - lhs[i][AA][0][3]*rhs[k][j][i-1][0]
			    - lhs[i][AA][1][3]*rhs[k][j][i-1][1]
			    - lhs[i][AA][2][3]*rhs[k][j][i-1][2]
			    - lhs[i][AA][3][3]*rhs[k][j][i-1][3]
			    - lhs[i][AA][4][3]*rhs[k][j][i-1][4];
	  rhs[k][j][i][4] = rhs[k][j][i][4] - lhs[i][AA][0][4]*rhs[k][j][i-1][0]
			    - lhs[i][AA][1][4]*rhs[k][j][i-1][1]
			    - lhs[i][AA][2][4]*rhs[k][j][i-1][2]
			    - lhs[i][AA][3][4]*rhs[k][j][i-1][3]
			    - lhs[i][AA][4][4]*rhs[k][j][i-1][4];
	}


        //-------------------------------------------------------------------
        // B(i) = B(i) - C(i-1)*A(i)
        // matmul_sub(lhs[i][AA],                 lhs[i-1][CC],        lhs[i][BB]);
        // void matmul_sub(double ablock[5][5], double bblock[5][5], double cblock[5][5])
        //-------------------------------------------------------------------
	{
	  lhs[i][BB][0][0] = lhs[i][BB][0][0] - lhs[i][AA][0][0]*lhs[i-1][CC][0][0]
				      - lhs[i][AA][1][0]*lhs[i-1][CC][0][1]
				      - lhs[i][AA][2][0]*lhs[i-1][CC][0][2]
				      - lhs[i][AA][3][0]*lhs[i-1][CC][0][3]
				      - lhs[i][AA][4][0]*lhs[i-1][CC][0][4];
	  lhs[i][BB][0][1] = lhs[i][BB][0][1] - lhs[i][AA][0][1]*lhs[i-1][CC][0][0]
				      - lhs[i][AA][1][1]*lhs[i-1][CC][0][1]
				      - lhs[i][AA][2][1]*lhs[i-1][CC][0][2]
				      - lhs[i][AA][3][1]*lhs[i-1][CC][0][3]
				      - lhs[i][AA][4][1]*lhs[i-1][CC][0][4];
	  lhs[i][BB][0][2] = lhs[i][BB][0][2] - lhs[i][AA][0][2]*lhs[i-1][CC][0][0]
				      - lhs[i][AA][1][2]*lhs[i-1][CC][0][1]
				      - lhs[i][AA][2][2]*lhs[i-1][CC][0][2]
				      - lhs[i][AA][3][2]*lhs[i-1][CC][0][3]
				      - lhs[i][AA][4][2]*lhs[i-1][CC][0][4];
	  lhs[i][BB][0][3] = lhs[i][BB][0][3] - lhs[i][AA][0][3]*lhs[i-1][CC][0][0]
				      - lhs[i][AA][1][3]*lhs[i-1][CC][0][1]
				      - lhs[i][AA][2][3]*lhs[i-1][CC][0][2]
				      - lhs[i][AA][3][3]*lhs[i-1][CC][0][3]
				      - lhs[i][AA][4][3]*lhs[i-1][CC][0][4];
	  lhs[i][BB][0][4] = lhs[i][BB][0][4] - lhs[i][AA][0][4]*lhs[i-1][CC][0][0]
				      - lhs[i][AA][1][4]*lhs[i-1][CC][0][1]
				      - lhs[i][AA][2][4]*lhs[i-1][CC][0][2]
				      - lhs[i][AA][3][4]*lhs[i-1][CC][0][3]
				      - lhs[i][AA][4][4]*lhs[i-1][CC][0][4];
	  lhs[i][BB][1][0] = lhs[i][BB][1][0] - lhs[i][AA][0][0]*lhs[i-1][CC][1][0]
				      - lhs[i][AA][1][0]*lhs[i-1][CC][1][1]
				      - lhs[i][AA][2][0]*lhs[i-1][CC][1][2]
				      - lhs[i][AA][3][0]*lhs[i-1][CC][1][3]
				      - lhs[i][AA][4][0]*lhs[i-1][CC][1][4];
	  lhs[i][BB][1][1] = lhs[i][BB][1][1] - lhs[i][AA][0][1]*lhs[i-1][CC][1][0]
				      - lhs[i][AA][1][1]*lhs[i-1][CC][1][1]
				      - lhs[i][AA][2][1]*lhs[i-1][CC][1][2]
				      - lhs[i][AA][3][1]*lhs[i-1][CC][1][3]
				      - lhs[i][AA][4][1]*lhs[i-1][CC][1][4];
	  lhs[i][BB][1][2] = lhs[i][BB][1][2] - lhs[i][AA][0][2]*lhs[i-1][CC][1][0]
				      - lhs[i][AA][1][2]*lhs[i-1][CC][1][1]
				      - lhs[i][AA][2][2]*lhs[i-1][CC][1][2]
				      - lhs[i][AA][3][2]*lhs[i-1][CC][1][3]
				      - lhs[i][AA][4][2]*lhs[i-1][CC][1][4];
	  lhs[i][BB][1][3] = lhs[i][BB][1][3] - lhs[i][AA][0][3]*lhs[i-1][CC][1][0]
				      - lhs[i][AA][1][3]*lhs[i-1][CC][1][1]
				      - lhs[i][AA][2][3]*lhs[i-1][CC][1][2]
				      - lhs[i][AA][3][3]*lhs[i-1][CC][1][3]
				      - lhs[i][AA][4][3]*lhs[i-1][CC][1][4];
	  lhs[i][BB][1][4] = lhs[i][BB][1][4] - lhs[i][AA][0][4]*lhs[i-1][CC][1][0]
				      - lhs[i][AA][1][4]*lhs[i-1][CC][1][1]
				      - lhs[i][AA][2][4]*lhs[i-1][CC][1][2]
				      - lhs[i][AA][3][4]*lhs[i-1][CC][1][3]
				      - lhs[i][AA][4][4]*lhs[i-1][CC][1][4];
	  lhs[i][BB][2][0] = lhs[i][BB][2][0] - lhs[i][AA][0][0]*lhs[i-1][CC][2][0]
				      - lhs[i][AA][1][0]*lhs[i-1][CC][2][1]
				      - lhs[i][AA][2][0]*lhs[i-1][CC][2][2]
				      - lhs[i][AA][3][0]*lhs[i-1][CC][2][3]
				      - lhs[i][AA][4][0]*lhs[i-1][CC][2][4];
	  lhs[i][BB][2][1] = lhs[i][BB][2][1] - lhs[i][AA][0][1]*lhs[i-1][CC][2][0]
				      - lhs[i][AA][1][1]*lhs[i-1][CC][2][1]
				      - lhs[i][AA][2][1]*lhs[i-1][CC][2][2]
				      - lhs[i][AA][3][1]*lhs[i-1][CC][2][3]
				      - lhs[i][AA][4][1]*lhs[i-1][CC][2][4];
	  lhs[i][BB][2][2] = lhs[i][BB][2][2] - lhs[i][AA][0][2]*lhs[i-1][CC][2][0]
				      - lhs[i][AA][1][2]*lhs[i-1][CC][2][1]
				      - lhs[i][AA][2][2]*lhs[i-1][CC][2][2]
				      - lhs[i][AA][3][2]*lhs[i-1][CC][2][3]
				      - lhs[i][AA][4][2]*lhs[i-1][CC][2][4];
	  lhs[i][BB][2][3] = lhs[i][BB][2][3] - lhs[i][AA][0][3]*lhs[i-1][CC][2][0]
				      - lhs[i][AA][1][3]*lhs[i-1][CC][2][1]
				      - lhs[i][AA][2][3]*lhs[i-1][CC][2][2]
				      - lhs[i][AA][3][3]*lhs[i-1][CC][2][3]
				      - lhs[i][AA][4][3]*lhs[i-1][CC][2][4];
	  lhs[i][BB][2][4] = lhs[i][BB][2][4] - lhs[i][AA][0][4]*lhs[i-1][CC][2][0]
				      - lhs[i][AA][1][4]*lhs[i-1][CC][2][1]
				      - lhs[i][AA][2][4]*lhs[i-1][CC][2][2]
				      - lhs[i][AA][3][4]*lhs[i-1][CC][2][3]
				      - lhs[i][AA][4][4]*lhs[i-1][CC][2][4];
	  lhs[i][BB][3][0] = lhs[i][BB][3][0] - lhs[i][AA][0][0]*lhs[i-1][CC][3][0]
				      - lhs[i][AA][1][0]*lhs[i-1][CC][3][1]
				      - lhs[i][AA][2][0]*lhs[i-1][CC][3][2]
				      - lhs[i][AA][3][0]*lhs[i-1][CC][3][3]
				      - lhs[i][AA][4][0]*lhs[i-1][CC][3][4];
	  lhs[i][BB][3][1] = lhs[i][BB][3][1] - lhs[i][AA][0][1]*lhs[i-1][CC][3][0]
				      - lhs[i][AA][1][1]*lhs[i-1][CC][3][1]
				      - lhs[i][AA][2][1]*lhs[i-1][CC][3][2]
				      - lhs[i][AA][3][1]*lhs[i-1][CC][3][3]
				      - lhs[i][AA][4][1]*lhs[i-1][CC][3][4];
	  lhs[i][BB][3][2] = lhs[i][BB][3][2] - lhs[i][AA][0][2]*lhs[i-1][CC][3][0]
				      - lhs[i][AA][1][2]*lhs[i-1][CC][3][1]
				      - lhs[i][AA][2][2]*lhs[i-1][CC][3][2]
				      - lhs[i][AA][3][2]*lhs[i-1][CC][3][3]
				      - lhs[i][AA][4][2]*lhs[i-1][CC][3][4];
	  lhs[i][BB][3][3] = lhs[i][BB][3][3] - lhs[i][AA][0][3]*lhs[i-1][CC][3][0]
				      - lhs[i][AA][1][3]*lhs[i-1][CC][3][1]
				      - lhs[i][AA][2][3]*lhs[i-1][CC][3][2]
				      - lhs[i][AA][3][3]*lhs[i-1][CC][3][3]
				      - lhs[i][AA][4][3]*lhs[i-1][CC][3][4];
	  lhs[i][BB][3][4] = lhs[i][BB][3][4] - lhs[i][AA][0][4]*lhs[i-1][CC][3][0]
				      - lhs[i][AA][1][4]*lhs[i-1][CC][3][1]
				      - lhs[i][AA][2][4]*lhs[i-1][CC][3][2]
				      - lhs[i][AA][3][4]*lhs[i-1][CC][3][3]
				      - lhs[i][AA][4][4]*lhs[i-1][CC][3][4];
	  lhs[i][BB][4][0] = lhs[i][BB][4][0] - lhs[i][AA][0][0]*lhs[i-1][CC][4][0]
				      - lhs[i][AA][1][0]*lhs[i-1][CC][4][1]
				      - lhs[i][AA][2][0]*lhs[i-1][CC][4][2]
				      - lhs[i][AA][3][0]*lhs[i-1][CC][4][3]
				      - lhs[i][AA][4][0]*lhs[i-1][CC][4][4];
	  lhs[i][BB][4][1] = lhs[i][BB][4][1] - lhs[i][AA][0][1]*lhs[i-1][CC][4][0]
				      - lhs[i][AA][1][1]*lhs[i-1][CC][4][1]
				      - lhs[i][AA][2][1]*lhs[i-1][CC][4][2]
				      - lhs[i][AA][3][1]*lhs[i-1][CC][4][3]
				      - lhs[i][AA][4][1]*lhs[i-1][CC][4][4];
	  lhs[i][BB][4][2] = lhs[i][BB][4][2] - lhs[i][AA][0][2]*lhs[i-1][CC][4][0]
				      - lhs[i][AA][1][2]*lhs[i-1][CC][4][1]
				      - lhs[i][AA][2][2]*lhs[i-1][CC][4][2]
				      - lhs[i][AA][3][2]*lhs[i-1][CC][4][3]
				      - lhs[i][AA][4][2]*lhs[i-1][CC][4][4];
	  lhs[i][BB][4][3] = lhs[i][BB][4][3] - lhs[i][AA][0][3]*lhs[i-1][CC][4][0]
				      - lhs[i][AA][1][3]*lhs[i-1][CC][4][1]
				      - lhs[i][AA][2][3]*lhs[i-1][CC][4][2]
				      - lhs[i][AA][3][3]*lhs[i-1][CC][4][3]
				      - lhs[i][AA][4][3]*lhs[i-1][CC][4][4];
	  lhs[i][BB][4][4] = lhs[i][BB][4][4] - lhs[i][AA][0][4]*lhs[i-1][CC][4][0]
				      - lhs[i][AA][1][4]*lhs[i-1][CC][4][1]
				      - lhs[i][AA][2][4]*lhs[i-1][CC][4][2]
				      - lhs[i][AA][3][4]*lhs[i-1][CC][4][3]
				      - lhs[i][AA][4][4]*lhs[i-1][CC][4][4];
	}


        //-------------------------------------------------------------------
        // multiply c[k][j][i] by b_inverse and copy back to c
        // multiply rhs[k][j][0] by b_inverse[k][j][0] and copy to rhs
	//
	// binvcrhs( lhs[i][BB], lhs[i][CC], rhs[k][j][i] );
	// void binvcrhs(double lhs[5][5], double c[5][5], double r[5])
        //-------------------------------------------------------------------
	{

	  pivot2 = 1.00/lhs[i][BB][0][0];
	  lhs[i][BB][1][0] = lhs[i][BB][1][0]*pivot2;
	  lhs[i][BB][2][0] = lhs[i][BB][2][0]*pivot2;
	  lhs[i][BB][3][0] = lhs[i][BB][3][0]*pivot2;
	  lhs[i][BB][4][0] = lhs[i][BB][4][0]*pivot2;
	  lhs[i][CC][0][0] = lhs[i][CC][0][0]*pivot2;
	  lhs[i][CC][1][0] = lhs[i][CC][1][0]*pivot2;
	  lhs[i][CC][2][0] = lhs[i][CC][2][0]*pivot2;
	  lhs[i][CC][3][0] = lhs[i][CC][3][0]*pivot2;
	  lhs[i][CC][4][0] = lhs[i][CC][4][0]*pivot2;
	  rhs[k][j][i][0]   = rhs[k][j][i][0]  *pivot2;

	  coeff2 = lhs[i][BB][0][1];
	  lhs[i][BB][1][1]= lhs[i][BB][1][1] - coeff2*lhs[i][BB][1][0];
	  lhs[i][BB][2][1]= lhs[i][BB][2][1] - coeff2*lhs[i][BB][2][0];
	  lhs[i][BB][3][1]= lhs[i][BB][3][1] - coeff2*lhs[i][BB][3][0];
	  lhs[i][BB][4][1]= lhs[i][BB][4][1] - coeff2*lhs[i][BB][4][0];
	  lhs[i][CC][0][1] = lhs[i][CC][0][1] - coeff2*lhs[i][CC][0][0];
	  lhs[i][CC][1][1] = lhs[i][CC][1][1] - coeff2*lhs[i][CC][1][0];
	  lhs[i][CC][2][1] = lhs[i][CC][2][1] - coeff2*lhs[i][CC][2][0];
	  lhs[i][CC][3][1] = lhs[i][CC][3][1] - coeff2*lhs[i][CC][3][0];
	  lhs[i][CC][4][1] = lhs[i][CC][4][1] - coeff2*lhs[i][CC][4][0];
	  rhs[k][j][i][1]   = rhs[k][j][i][1]   - coeff2*rhs[k][j][i][0];

	  coeff2 = lhs[i][BB][0][2];
	  lhs[i][BB][1][2]= lhs[i][BB][1][2] - coeff2*lhs[i][BB][1][0];
	  lhs[i][BB][2][2]= lhs[i][BB][2][2] - coeff2*lhs[i][BB][2][0];
	  lhs[i][BB][3][2]= lhs[i][BB][3][2] - coeff2*lhs[i][BB][3][0];
	  lhs[i][BB][4][2]= lhs[i][BB][4][2] - coeff2*lhs[i][BB][4][0];
	  lhs[i][CC][0][2] = lhs[i][CC][0][2] - coeff2*lhs[i][CC][0][0];
	  lhs[i][CC][1][2] = lhs[i][CC][1][2] - coeff2*lhs[i][CC][1][0];
	  lhs[i][CC][2][2] = lhs[i][CC][2][2] - coeff2*lhs[i][CC][2][0];
	  lhs[i][CC][3][2] = lhs[i][CC][3][2] - coeff2*lhs[i][CC][3][0];
	  lhs[i][CC][4][2] = lhs[i][CC][4][2] - coeff2*lhs[i][CC][4][0];
	  rhs[k][j][i][2]   = rhs[k][j][i][2]   - coeff2*rhs[k][j][i][0];

	  coeff2 = lhs[i][BB][0][3];
	  lhs[i][BB][1][3]= lhs[i][BB][1][3] - coeff2*lhs[i][BB][1][0];
	  lhs[i][BB][2][3]= lhs[i][BB][2][3] - coeff2*lhs[i][BB][2][0];
	  lhs[i][BB][3][3]= lhs[i][BB][3][3] - coeff2*lhs[i][BB][3][0];
	  lhs[i][BB][4][3]= lhs[i][BB][4][3] - coeff2*lhs[i][BB][4][0];
	  lhs[i][CC][0][3] = lhs[i][CC][0][3] - coeff2*lhs[i][CC][0][0];
	  lhs[i][CC][1][3] = lhs[i][CC][1][3] - coeff2*lhs[i][CC][1][0];
	  lhs[i][CC][2][3] = lhs[i][CC][2][3] - coeff2*lhs[i][CC][2][0];
	  lhs[i][CC][3][3] = lhs[i][CC][3][3] - coeff2*lhs[i][CC][3][0];
	  lhs[i][CC][4][3] = lhs[i][CC][4][3] - coeff2*lhs[i][CC][4][0];
	  rhs[k][j][i][3]   = rhs[k][j][i][3]   - coeff2*rhs[k][j][i][0];

	  coeff2 = lhs[i][BB][0][4];
	  lhs[i][BB][1][4]= lhs[i][BB][1][4] - coeff2*lhs[i][BB][1][0];
	  lhs[i][BB][2][4]= lhs[i][BB][2][4] - coeff2*lhs[i][BB][2][0];
	  lhs[i][BB][3][4]= lhs[i][BB][3][4] - coeff2*lhs[i][BB][3][0];
	  lhs[i][BB][4][4]= lhs[i][BB][4][4] - coeff2*lhs[i][BB][4][0];
	  lhs[i][CC][0][4] = lhs[i][CC][0][4] - coeff2*lhs[i][CC][0][0];
	  lhs[i][CC][1][4] = lhs[i][CC][1][4] - coeff2*lhs[i][CC][1][0];
	  lhs[i][CC][2][4] = lhs[i][CC][2][4] - coeff2*lhs[i][CC][2][0];
	  lhs[i][CC][3][4] = lhs[i][CC][3][4] - coeff2*lhs[i][CC][3][0];
	  lhs[i][CC][4][4] = lhs[i][CC][4][4] - coeff2*lhs[i][CC][4][0];
	  rhs[k][j][i][4]   = rhs[k][j][i][4]   - coeff2*rhs[k][j][i][0];


	  pivot2 = 1.00/lhs[i][BB][1][1];
	  lhs[i][BB][2][1] = lhs[i][BB][2][1]*pivot2;
	  lhs[i][BB][3][1] = lhs[i][BB][3][1]*pivot2;
	  lhs[i][BB][4][1] = lhs[i][BB][4][1]*pivot2;
	  lhs[i][CC][0][1] = lhs[i][CC][0][1]*pivot2;
	  lhs[i][CC][1][1] = lhs[i][CC][1][1]*pivot2;
	  lhs[i][CC][2][1] = lhs[i][CC][2][1]*pivot2;
	  lhs[i][CC][3][1] = lhs[i][CC][3][1]*pivot2;
	  lhs[i][CC][4][1] = lhs[i][CC][4][1]*pivot2;
	  rhs[k][j][i][1]   = rhs[k][j][i][1]  *pivot2;

	  coeff2 = lhs[i][BB][1][0];
	  lhs[i][BB][2][0]= lhs[i][BB][2][0] - coeff2*lhs[i][BB][2][1];
	  lhs[i][BB][3][0]= lhs[i][BB][3][0] - coeff2*lhs[i][BB][3][1];
	  lhs[i][BB][4][0]= lhs[i][BB][4][0] - coeff2*lhs[i][BB][4][1];
	  lhs[i][CC][0][0] = lhs[i][CC][0][0] - coeff2*lhs[i][CC][0][1];
	  lhs[i][CC][1][0] = lhs[i][CC][1][0] - coeff2*lhs[i][CC][1][1];
	  lhs[i][CC][2][0] = lhs[i][CC][2][0] - coeff2*lhs[i][CC][2][1];
	  lhs[i][CC][3][0] = lhs[i][CC][3][0] - coeff2*lhs[i][CC][3][1];
	  lhs[i][CC][4][0] = lhs[i][CC][4][0] - coeff2*lhs[i][CC][4][1];
	  rhs[k][j][i][0]   = rhs[k][j][i][0]   - coeff2*rhs[k][j][i][1];

	  coeff2 = lhs[i][BB][1][2];
	  lhs[i][BB][2][2]= lhs[i][BB][2][2] - coeff2*lhs[i][BB][2][1];
	  lhs[i][BB][3][2]= lhs[i][BB][3][2] - coeff2*lhs[i][BB][3][1];
	  lhs[i][BB][4][2]= lhs[i][BB][4][2] - coeff2*lhs[i][BB][4][1];
	  lhs[i][CC][0][2] = lhs[i][CC][0][2] - coeff2*lhs[i][CC][0][1];
	  lhs[i][CC][1][2] = lhs[i][CC][1][2] - coeff2*lhs[i][CC][1][1];
	  lhs[i][CC][2][2] = lhs[i][CC][2][2] - coeff2*lhs[i][CC][2][1];
	  lhs[i][CC][3][2] = lhs[i][CC][3][2] - coeff2*lhs[i][CC][3][1];
	  lhs[i][CC][4][2] = lhs[i][CC][4][2] - coeff2*lhs[i][CC][4][1];
	  rhs[k][j][i][2]   = rhs[k][j][i][2]   - coeff2*rhs[k][j][i][1];

	  coeff2 = lhs[i][BB][1][3];
	  lhs[i][BB][2][3]= lhs[i][BB][2][3] - coeff2*lhs[i][BB][2][1];
	  lhs[i][BB][3][3]= lhs[i][BB][3][3] - coeff2*lhs[i][BB][3][1];
	  lhs[i][BB][4][3]= lhs[i][BB][4][3] - coeff2*lhs[i][BB][4][1];
	  lhs[i][CC][0][3] = lhs[i][CC][0][3] - coeff2*lhs[i][CC][0][1];
	  lhs[i][CC][1][3] = lhs[i][CC][1][3] - coeff2*lhs[i][CC][1][1];
	  lhs[i][CC][2][3] = lhs[i][CC][2][3] - coeff2*lhs[i][CC][2][1];
	  lhs[i][CC][3][3] = lhs[i][CC][3][3] - coeff2*lhs[i][CC][3][1];
	  lhs[i][CC][4][3] = lhs[i][CC][4][3] - coeff2*lhs[i][CC][4][1];
	  rhs[k][j][i][3]   = rhs[k][j][i][3]   - coeff2*rhs[k][j][i][1];

	  coeff2 = lhs[i][BB][1][4];
	  lhs[i][BB][2][4]= lhs[i][BB][2][4] - coeff2*lhs[i][BB][2][1];
	  lhs[i][BB][3][4]= lhs[i][BB][3][4] - coeff2*lhs[i][BB][3][1];
	  lhs[i][BB][4][4]= lhs[i][BB][4][4] - coeff2*lhs[i][BB][4][1];
	  lhs[i][CC][0][4] = lhs[i][CC][0][4] - coeff2*lhs[i][CC][0][1];
	  lhs[i][CC][1][4] = lhs[i][CC][1][4] - coeff2*lhs[i][CC][1][1];
	  lhs[i][CC][2][4] = lhs[i][CC][2][4] - coeff2*lhs[i][CC][2][1];
	  lhs[i][CC][3][4] = lhs[i][CC][3][4] - coeff2*lhs[i][CC][3][1];
	  lhs[i][CC][4][4] = lhs[i][CC][4][4] - coeff2*lhs[i][CC][4][1];
	  rhs[k][j][i][4]   = rhs[k][j][i][4]   - coeff2*rhs[k][j][i][1];


	  pivot2 = 1.00/lhs[i][BB][2][2];
	  lhs[i][BB][3][2] = lhs[i][BB][3][2]*pivot2;
	  lhs[i][BB][4][2] = lhs[i][BB][4][2]*pivot2;
	  lhs[i][CC][0][2] = lhs[i][CC][0][2]*pivot2;
	  lhs[i][CC][1][2] = lhs[i][CC][1][2]*pivot2;
	  lhs[i][CC][2][2] = lhs[i][CC][2][2]*pivot2;
	  lhs[i][CC][3][2] = lhs[i][CC][3][2]*pivot2;
	  lhs[i][CC][4][2] = lhs[i][CC][4][2]*pivot2;
	  rhs[k][j][i][2]   = rhs[k][j][i][2]  *pivot2;

	  coeff2 = lhs[i][BB][2][0];
	  lhs[i][BB][3][0]= lhs[i][BB][3][0] - coeff2*lhs[i][BB][3][2];
	  lhs[i][BB][4][0]= lhs[i][BB][4][0] - coeff2*lhs[i][BB][4][2];
	  lhs[i][CC][0][0] = lhs[i][CC][0][0] - coeff2*lhs[i][CC][0][2];
	  lhs[i][CC][1][0] = lhs[i][CC][1][0] - coeff2*lhs[i][CC][1][2];
	  lhs[i][CC][2][0] = lhs[i][CC][2][0] - coeff2*lhs[i][CC][2][2];
	  lhs[i][CC][3][0] = lhs[i][CC][3][0] - coeff2*lhs[i][CC][3][2];
	  lhs[i][CC][4][0] = lhs[i][CC][4][0] - coeff2*lhs[i][CC][4][2];
	  rhs[k][j][i][0]   = rhs[k][j][i][0]   - coeff2*rhs[k][j][i][2];

	  coeff2 = lhs[i][BB][2][1];
	  lhs[i][BB][3][1]= lhs[i][BB][3][1] - coeff2*lhs[i][BB][3][2];
	  lhs[i][BB][4][1]= lhs[i][BB][4][1] - coeff2*lhs[i][BB][4][2];
	  lhs[i][CC][0][1] = lhs[i][CC][0][1] - coeff2*lhs[i][CC][0][2];
	  lhs[i][CC][1][1] = lhs[i][CC][1][1] - coeff2*lhs[i][CC][1][2];
	  lhs[i][CC][2][1] = lhs[i][CC][2][1] - coeff2*lhs[i][CC][2][2];
	  lhs[i][CC][3][1] = lhs[i][CC][3][1] - coeff2*lhs[i][CC][3][2];
	  lhs[i][CC][4][1] = lhs[i][CC][4][1] - coeff2*lhs[i][CC][4][2];
	  rhs[k][j][i][1]   = rhs[k][j][i][1]   - coeff2*rhs[k][j][i][2];

	  coeff2 = lhs[i][BB][2][3];
	  lhs[i][BB][3][3]= lhs[i][BB][3][3] - coeff2*lhs[i][BB][3][2];
	  lhs[i][BB][4][3]= lhs[i][BB][4][3] - coeff2*lhs[i][BB][4][2];
	  lhs[i][CC][0][3] = lhs[i][CC][0][3] - coeff2*lhs[i][CC][0][2];
	  lhs[i][CC][1][3] = lhs[i][CC][1][3] - coeff2*lhs[i][CC][1][2];
	  lhs[i][CC][2][3] = lhs[i][CC][2][3] - coeff2*lhs[i][CC][2][2];
	  lhs[i][CC][3][3] = lhs[i][CC][3][3] - coeff2*lhs[i][CC][3][2];
	  lhs[i][CC][4][3] = lhs[i][CC][4][3] - coeff2*lhs[i][CC][4][2];
	  rhs[k][j][i][3]   = rhs[k][j][i][3]   - coeff2*rhs[k][j][i][2];

	  coeff2 = lhs[i][BB][2][4];
	  lhs[i][BB][3][4]= lhs[i][BB][3][4] - coeff2*lhs[i][BB][3][2];
	  lhs[i][BB][4][4]= lhs[i][BB][4][4] - coeff2*lhs[i][BB][4][2];
	  lhs[i][CC][0][4] = lhs[i][CC][0][4] - coeff2*lhs[i][CC][0][2];
	  lhs[i][CC][1][4] = lhs[i][CC][1][4] - coeff2*lhs[i][CC][1][2];
	  lhs[i][CC][2][4] = lhs[i][CC][2][4] - coeff2*lhs[i][CC][2][2];
	  lhs[i][CC][3][4] = lhs[i][CC][3][4] - coeff2*lhs[i][CC][3][2];
	  lhs[i][CC][4][4] = lhs[i][CC][4][4] - coeff2*lhs[i][CC][4][2];
	  rhs[k][j][i][4]   = rhs[k][j][i][4]   - coeff2*rhs[k][j][i][2];


	  pivot2 = 1.00/lhs[i][BB][3][3];
	  lhs[i][BB][4][3] = lhs[i][BB][4][3]*pivot2;
	  lhs[i][CC][0][3] = lhs[i][CC][0][3]*pivot2;
	  lhs[i][CC][1][3] = lhs[i][CC][1][3]*pivot2;
	  lhs[i][CC][2][3] = lhs[i][CC][2][3]*pivot2;
	  lhs[i][CC][3][3] = lhs[i][CC][3][3]*pivot2;
	  lhs[i][CC][4][3] = lhs[i][CC][4][3]*pivot2;
	  rhs[k][j][i][3]   = rhs[k][j][i][3]  *pivot2;

	  coeff2 = lhs[i][BB][3][0];
	  lhs[i][BB][4][0]= lhs[i][BB][4][0] - coeff2*lhs[i][BB][4][3];
	  lhs[i][CC][0][0] = lhs[i][CC][0][0] - coeff2*lhs[i][CC][0][3];
	  lhs[i][CC][1][0] = lhs[i][CC][1][0] - coeff2*lhs[i][CC][1][3];
	  lhs[i][CC][2][0] = lhs[i][CC][2][0] - coeff2*lhs[i][CC][2][3];
	  lhs[i][CC][3][0] = lhs[i][CC][3][0] - coeff2*lhs[i][CC][3][3];
	  lhs[i][CC][4][0] = lhs[i][CC][4][0] - coeff2*lhs[i][CC][4][3];
	  rhs[k][j][i][0]   = rhs[k][j][i][0]   - coeff2*rhs[k][j][i][3];

	  coeff2 = lhs[i][BB][3][1];
	  lhs[i][BB][4][1]= lhs[i][BB][4][1] - coeff2*lhs[i][BB][4][3];
	  lhs[i][CC][0][1] = lhs[i][CC][0][1] - coeff2*lhs[i][CC][0][3];
	  lhs[i][CC][1][1] = lhs[i][CC][1][1] - coeff2*lhs[i][CC][1][3];
	  lhs[i][CC][2][1] = lhs[i][CC][2][1] - coeff2*lhs[i][CC][2][3];
	  lhs[i][CC][3][1] = lhs[i][CC][3][1] - coeff2*lhs[i][CC][3][3];
	  lhs[i][CC][4][1] = lhs[i][CC][4][1] - coeff2*lhs[i][CC][4][3];
	  rhs[k][j][i][1]   = rhs[k][j][i][1]   - coeff2*rhs[k][j][i][3];

	  coeff2 = lhs[i][BB][3][2];
	  lhs[i][BB][4][2]= lhs[i][BB][4][2] - coeff2*lhs[i][BB][4][3];
	  lhs[i][CC][0][2] = lhs[i][CC][0][2] - coeff2*lhs[i][CC][0][3];
	  lhs[i][CC][1][2] = lhs[i][CC][1][2] - coeff2*lhs[i][CC][1][3];
	  lhs[i][CC][2][2] = lhs[i][CC][2][2] - coeff2*lhs[i][CC][2][3];
	  lhs[i][CC][3][2] = lhs[i][CC][3][2] - coeff2*lhs[i][CC][3][3];
	  lhs[i][CC][4][2] = lhs[i][CC][4][2] - coeff2*lhs[i][CC][4][3];
	  rhs[k][j][i][2]   = rhs[k][j][i][2]   - coeff2*rhs[k][j][i][3];

	  coeff2 = lhs[i][BB][3][4];
	  lhs[i][BB][4][4]= lhs[i][BB][4][4] - coeff2*lhs[i][BB][4][3];
	  lhs[i][CC][0][4] = lhs[i][CC][0][4] - coeff2*lhs[i][CC][0][3];
	  lhs[i][CC][1][4] = lhs[i][CC][1][4] - coeff2*lhs[i][CC][1][3];
	  lhs[i][CC][2][4] = lhs[i][CC][2][4] - coeff2*lhs[i][CC][2][3];
	  lhs[i][CC][3][4] = lhs[i][CC][3][4] - coeff2*lhs[i][CC][3][3];
	  lhs[i][CC][4][4] = lhs[i][CC][4][4] - coeff2*lhs[i][CC][4][3];
	  rhs[k][j][i][4]   = rhs[k][j][i][4]   - coeff2*rhs[k][j][i][3];


	  pivot2 = 1.00/lhs[i][BB][4][4];
	  lhs[i][CC][0][4] = lhs[i][CC][0][4]*pivot2;
	  lhs[i][CC][1][4] = lhs[i][CC][1][4]*pivot2;
	  lhs[i][CC][2][4] = lhs[i][CC][2][4]*pivot2;
	  lhs[i][CC][3][4] = lhs[i][CC][3][4]*pivot2;
	  lhs[i][CC][4][4] = lhs[i][CC][4][4]*pivot2;
	  rhs[k][j][i][4]   = rhs[k][j][i][4]  *pivot2;

	  coeff2 = lhs[i][BB][4][0];
	  lhs[i][CC][0][0] = lhs[i][CC][0][0] - coeff2*lhs[i][CC][0][4];
	  lhs[i][CC][1][0] = lhs[i][CC][1][0] - coeff2*lhs[i][CC][1][4];
	  lhs[i][CC][2][0] = lhs[i][CC][2][0] - coeff2*lhs[i][CC][2][4];
	  lhs[i][CC][3][0] = lhs[i][CC][3][0] - coeff2*lhs[i][CC][3][4];
	  lhs[i][CC][4][0] = lhs[i][CC][4][0] - coeff2*lhs[i][CC][4][4];
	  rhs[k][j][i][0]   = rhs[k][j][i][0]   - coeff2*rhs[k][j][i][4];

	  coeff2 = lhs[i][BB][4][1];
	  lhs[i][CC][0][1] = lhs[i][CC][0][1] - coeff2*lhs[i][CC][0][4];
	  lhs[i][CC][1][1] = lhs[i][CC][1][1] - coeff2*lhs[i][CC][1][4];
	  lhs[i][CC][2][1] = lhs[i][CC][2][1] - coeff2*lhs[i][CC][2][4];
	  lhs[i][CC][3][1] = lhs[i][CC][3][1] - coeff2*lhs[i][CC][3][4];
	  lhs[i][CC][4][1] = lhs[i][CC][4][1] - coeff2*lhs[i][CC][4][4];
	  rhs[k][j][i][1]   = rhs[k][j][i][1]   - coeff2*rhs[k][j][i][4];

	  coeff2 = lhs[i][BB][4][2];
	  lhs[i][CC][0][2] = lhs[i][CC][0][2] - coeff2*lhs[i][CC][0][4];
	  lhs[i][CC][1][2] = lhs[i][CC][1][2] - coeff2*lhs[i][CC][1][4];
	  lhs[i][CC][2][2] = lhs[i][CC][2][2] - coeff2*lhs[i][CC][2][4];
	  lhs[i][CC][3][2] = lhs[i][CC][3][2] - coeff2*lhs[i][CC][3][4];
	  lhs[i][CC][4][2] = lhs[i][CC][4][2] - coeff2*lhs[i][CC][4][4];
	  rhs[k][j][i][2]   = rhs[k][j][i][2]   - coeff2*rhs[k][j][i][4];

	  coeff2 = lhs[i][BB][4][3];
	  lhs[i][CC][0][3] = lhs[i][CC][0][3] - coeff2*lhs[i][CC][0][4];
	  lhs[i][CC][1][3] = lhs[i][CC][1][3] - coeff2*lhs[i][CC][1][4];
	  lhs[i][CC][2][3] = lhs[i][CC][2][3] - coeff2*lhs[i][CC][2][4];
	  lhs[i][CC][3][3] = lhs[i][CC][3][3] - coeff2*lhs[i][CC][3][4];
	  lhs[i][CC][4][3] = lhs[i][CC][4][3] - coeff2*lhs[i][CC][4][4];
	  rhs[k][j][i][3]   = rhs[k][j][i][3]   - coeff2*rhs[k][j][i][4];
	} //END of binvcrhs( lhs[i][BB], lhs[i][CC], rhs[k][j][i] );
      } //END of for


      //---------------------------------------------------------------------
      // rhs(isize) = rhs(isize) - A*rhs(isize-1)
      //
      // matvec_sub(lhs[isize][AA], rhs[k][j][isize-1], rhs[k][j][isize]);
      // void matvec_sub(double ablock[5][5], double avec[5], double bvec[5])
      //---------------------------------------------------------------------
      {
	  rhs[k][j][isize][0] = rhs[k][j][isize][0] - lhs[isize][AA][0][0]*rhs[k][j][isize-1][0]
			    - lhs[isize][AA][1][0]*rhs[k][j][isize-1][1]
			    - lhs[isize][AA][2][0]*rhs[k][j][isize-1][2]
			    - lhs[isize][AA][3][0]*rhs[k][j][isize-1][3]
			    - lhs[isize][AA][4][0]*rhs[k][j][isize-1][4];
	  rhs[k][j][isize][1] = rhs[k][j][isize][1] - lhs[isize][AA][0][1]*rhs[k][j][isize-1][0]
			    - lhs[isize][AA][1][1]*rhs[k][j][isize-1][1]
			    - lhs[isize][AA][2][1]*rhs[k][j][isize-1][2]
			    - lhs[isize][AA][3][1]*rhs[k][j][isize-1][3]
			    - lhs[isize][AA][4][1]*rhs[k][j][isize-1][4];
	  rhs[k][j][isize][2] = rhs[k][j][isize][2] - lhs[isize][AA][0][2]*rhs[k][j][isize-1][0]
			    - lhs[isize][AA][1][2]*rhs[k][j][isize-1][1]
			    - lhs[isize][AA][2][2]*rhs[k][j][isize-1][2]
			    - lhs[isize][AA][3][2]*rhs[k][j][isize-1][3]
			    - lhs[isize][AA][4][2]*rhs[k][j][isize-1][4];
	  rhs[k][j][isize][3] = rhs[k][j][isize][3] - lhs[isize][AA][0][3]*rhs[k][j][isize-1][0]
			    - lhs[isize][AA][1][3]*rhs[k][j][isize-1][1]
			    - lhs[isize][AA][2][3]*rhs[k][j][isize-1][2]
			    - lhs[isize][AA][3][3]*rhs[k][j][isize-1][3]
			    - lhs[isize][AA][4][3]*rhs[k][j][isize-1][4];
	  rhs[k][j][isize][4] = rhs[k][j][isize][4] - lhs[isize][AA][0][4]*rhs[k][j][isize-1][0]
			    - lhs[isize][AA][1][4]*rhs[k][j][isize-1][1]
			    - lhs[isize][AA][2][4]*rhs[k][j][isize-1][2]
			    - lhs[isize][AA][3][4]*rhs[k][j][isize-1][3]
			    - lhs[isize][AA][4][4]*rhs[k][j][isize-1][4];
      }


      //---------------------------------------------------------------------
      // B(isize) = B(isize) - C(isize-1)*A(isize)
      //
      // matmul_sub(lhs[isize][AA], lhs[isize-1][CC], lhs[isize][BB]);
      // void matmul_sub(double ablock[5][5], double bblock[5][5], double cblock[5][5])
      //---------------------------------------------------------------------
      {
	  lhs[isize][BB][0][0] = lhs[isize][BB][0][0] - lhs[isize][AA][0][0]*lhs[isize-1][CC][0][0]
				      - lhs[isize][AA][1][0]*lhs[isize-1][CC][0][1]
				      - lhs[isize][AA][2][0]*lhs[isize-1][CC][0][2]
				      - lhs[isize][AA][3][0]*lhs[isize-1][CC][0][3]
				      - lhs[isize][AA][4][0]*lhs[isize-1][CC][0][4];
	  lhs[isize][BB][0][1] = lhs[isize][BB][0][1] - lhs[isize][AA][0][1]*lhs[isize-1][CC][0][0]
				      - lhs[isize][AA][1][1]*lhs[isize-1][CC][0][1]
				      - lhs[isize][AA][2][1]*lhs[isize-1][CC][0][2]
				      - lhs[isize][AA][3][1]*lhs[isize-1][CC][0][3]
				      - lhs[isize][AA][4][1]*lhs[isize-1][CC][0][4];
	  lhs[isize][BB][0][2] = lhs[isize][BB][0][2] - lhs[isize][AA][0][2]*lhs[isize-1][CC][0][0]
				      - lhs[isize][AA][1][2]*lhs[isize-1][CC][0][1]
				      - lhs[isize][AA][2][2]*lhs[isize-1][CC][0][2]
				      - lhs[isize][AA][3][2]*lhs[isize-1][CC][0][3]
				      - lhs[isize][AA][4][2]*lhs[isize-1][CC][0][4];
	  lhs[isize][BB][0][3] = lhs[isize][BB][0][3] - lhs[isize][AA][0][3]*lhs[isize-1][CC][0][0]
				      - lhs[isize][AA][1][3]*lhs[isize-1][CC][0][1]
				      - lhs[isize][AA][2][3]*lhs[isize-1][CC][0][2]
				      - lhs[isize][AA][3][3]*lhs[isize-1][CC][0][3]
				      - lhs[isize][AA][4][3]*lhs[isize-1][CC][0][4];
	  lhs[isize][BB][0][4] = lhs[isize][BB][0][4] - lhs[isize][AA][0][4]*lhs[isize-1][CC][0][0]
				      - lhs[isize][AA][1][4]*lhs[isize-1][CC][0][1]
				      - lhs[isize][AA][2][4]*lhs[isize-1][CC][0][2]
				      - lhs[isize][AA][3][4]*lhs[isize-1][CC][0][3]
				      - lhs[isize][AA][4][4]*lhs[isize-1][CC][0][4];
	  lhs[isize][BB][1][0] = lhs[isize][BB][1][0] - lhs[isize][AA][0][0]*lhs[isize-1][CC][1][0]
				      - lhs[isize][AA][1][0]*lhs[isize-1][CC][1][1]
				      - lhs[isize][AA][2][0]*lhs[isize-1][CC][1][2]
				      - lhs[isize][AA][3][0]*lhs[isize-1][CC][1][3]
				      - lhs[isize][AA][4][0]*lhs[isize-1][CC][1][4];
	  lhs[isize][BB][1][1] = lhs[isize][BB][1][1] - lhs[isize][AA][0][1]*lhs[isize-1][CC][1][0]
				      - lhs[isize][AA][1][1]*lhs[isize-1][CC][1][1]
				      - lhs[isize][AA][2][1]*lhs[isize-1][CC][1][2]
				      - lhs[isize][AA][3][1]*lhs[isize-1][CC][1][3]
				      - lhs[isize][AA][4][1]*lhs[isize-1][CC][1][4];
	  lhs[isize][BB][1][2] = lhs[isize][BB][1][2] - lhs[isize][AA][0][2]*lhs[isize-1][CC][1][0]
				      - lhs[isize][AA][1][2]*lhs[isize-1][CC][1][1]
				      - lhs[isize][AA][2][2]*lhs[isize-1][CC][1][2]
				      - lhs[isize][AA][3][2]*lhs[isize-1][CC][1][3]
				      - lhs[isize][AA][4][2]*lhs[isize-1][CC][1][4];
	  lhs[isize][BB][1][3] = lhs[isize][BB][1][3] - lhs[isize][AA][0][3]*lhs[isize-1][CC][1][0]
				      - lhs[isize][AA][1][3]*lhs[isize-1][CC][1][1]
				      - lhs[isize][AA][2][3]*lhs[isize-1][CC][1][2]
				      - lhs[isize][AA][3][3]*lhs[isize-1][CC][1][3]
				      - lhs[isize][AA][4][3]*lhs[isize-1][CC][1][4];
	  lhs[isize][BB][1][4] = lhs[isize][BB][1][4] - lhs[isize][AA][0][4]*lhs[isize-1][CC][1][0]
				      - lhs[isize][AA][1][4]*lhs[isize-1][CC][1][1]
				      - lhs[isize][AA][2][4]*lhs[isize-1][CC][1][2]
				      - lhs[isize][AA][3][4]*lhs[isize-1][CC][1][3]
				      - lhs[isize][AA][4][4]*lhs[isize-1][CC][1][4];
	  lhs[isize][BB][2][0] = lhs[isize][BB][2][0] - lhs[isize][AA][0][0]*lhs[isize-1][CC][2][0]
				      - lhs[isize][AA][1][0]*lhs[isize-1][CC][2][1]
				      - lhs[isize][AA][2][0]*lhs[isize-1][CC][2][2]
				      - lhs[isize][AA][3][0]*lhs[isize-1][CC][2][3]
				      - lhs[isize][AA][4][0]*lhs[isize-1][CC][2][4];
	  lhs[isize][BB][2][1] = lhs[isize][BB][2][1] - lhs[isize][AA][0][1]*lhs[isize-1][CC][2][0]
				      - lhs[isize][AA][1][1]*lhs[isize-1][CC][2][1]
				      - lhs[isize][AA][2][1]*lhs[isize-1][CC][2][2]
				      - lhs[isize][AA][3][1]*lhs[isize-1][CC][2][3]
				      - lhs[isize][AA][4][1]*lhs[isize-1][CC][2][4];
	  lhs[isize][BB][2][2] = lhs[isize][BB][2][2] - lhs[isize][AA][0][2]*lhs[isize-1][CC][2][0]
				      - lhs[isize][AA][1][2]*lhs[isize-1][CC][2][1]
				      - lhs[isize][AA][2][2]*lhs[isize-1][CC][2][2]
				      - lhs[isize][AA][3][2]*lhs[isize-1][CC][2][3]
				      - lhs[isize][AA][4][2]*lhs[isize-1][CC][2][4];
	  lhs[isize][BB][2][3] = lhs[isize][BB][2][3] - lhs[isize][AA][0][3]*lhs[isize-1][CC][2][0]
				      - lhs[isize][AA][1][3]*lhs[isize-1][CC][2][1]
				      - lhs[isize][AA][2][3]*lhs[isize-1][CC][2][2]
				      - lhs[isize][AA][3][3]*lhs[isize-1][CC][2][3]
				      - lhs[isize][AA][4][3]*lhs[isize-1][CC][2][4];
	  lhs[isize][BB][2][4] = lhs[isize][BB][2][4] - lhs[isize][AA][0][4]*lhs[isize-1][CC][2][0]
				      - lhs[isize][AA][1][4]*lhs[isize-1][CC][2][1]
				      - lhs[isize][AA][2][4]*lhs[isize-1][CC][2][2]
				      - lhs[isize][AA][3][4]*lhs[isize-1][CC][2][3]
				      - lhs[isize][AA][4][4]*lhs[isize-1][CC][2][4];
	  lhs[isize][BB][3][0] = lhs[isize][BB][3][0] - lhs[isize][AA][0][0]*lhs[isize-1][CC][3][0]
				      - lhs[isize][AA][1][0]*lhs[isize-1][CC][3][1]
				      - lhs[isize][AA][2][0]*lhs[isize-1][CC][3][2]
				      - lhs[isize][AA][3][0]*lhs[isize-1][CC][3][3]
				      - lhs[isize][AA][4][0]*lhs[isize-1][CC][3][4];
	  lhs[isize][BB][3][1] = lhs[isize][BB][3][1] - lhs[isize][AA][0][1]*lhs[isize-1][CC][3][0]
				      - lhs[isize][AA][1][1]*lhs[isize-1][CC][3][1]
				      - lhs[isize][AA][2][1]*lhs[isize-1][CC][3][2]
				      - lhs[isize][AA][3][1]*lhs[isize-1][CC][3][3]
				      - lhs[isize][AA][4][1]*lhs[isize-1][CC][3][4];
	  lhs[isize][BB][3][2] = lhs[isize][BB][3][2] - lhs[isize][AA][0][2]*lhs[isize-1][CC][3][0]
				      - lhs[isize][AA][1][2]*lhs[isize-1][CC][3][1]
				      - lhs[isize][AA][2][2]*lhs[isize-1][CC][3][2]
				      - lhs[isize][AA][3][2]*lhs[isize-1][CC][3][3]
				      - lhs[isize][AA][4][2]*lhs[isize-1][CC][3][4];
	  lhs[isize][BB][3][3] = lhs[isize][BB][3][3] - lhs[isize][AA][0][3]*lhs[isize-1][CC][3][0]
				      - lhs[isize][AA][1][3]*lhs[isize-1][CC][3][1]
				      - lhs[isize][AA][2][3]*lhs[isize-1][CC][3][2]
				      - lhs[isize][AA][3][3]*lhs[isize-1][CC][3][3]
				      - lhs[isize][AA][4][3]*lhs[isize-1][CC][3][4];
	  lhs[isize][BB][3][4] = lhs[isize][BB][3][4] - lhs[isize][AA][0][4]*lhs[isize-1][CC][3][0]
				      - lhs[isize][AA][1][4]*lhs[isize-1][CC][3][1]
				      - lhs[isize][AA][2][4]*lhs[isize-1][CC][3][2]
				      - lhs[isize][AA][3][4]*lhs[isize-1][CC][3][3]
				      - lhs[isize][AA][4][4]*lhs[isize-1][CC][3][4];
	  lhs[isize][BB][4][0] = lhs[isize][BB][4][0] - lhs[isize][AA][0][0]*lhs[isize-1][CC][4][0]
				      - lhs[isize][AA][1][0]*lhs[isize-1][CC][4][1]
				      - lhs[isize][AA][2][0]*lhs[isize-1][CC][4][2]
				      - lhs[isize][AA][3][0]*lhs[isize-1][CC][4][3]
				      - lhs[isize][AA][4][0]*lhs[isize-1][CC][4][4];
	  lhs[isize][BB][4][1] = lhs[isize][BB][4][1] - lhs[isize][AA][0][1]*lhs[isize-1][CC][4][0]
				      - lhs[isize][AA][1][1]*lhs[isize-1][CC][4][1]
				      - lhs[isize][AA][2][1]*lhs[isize-1][CC][4][2]
				      - lhs[isize][AA][3][1]*lhs[isize-1][CC][4][3]
				      - lhs[isize][AA][4][1]*lhs[isize-1][CC][4][4];
	  lhs[isize][BB][4][2] = lhs[isize][BB][4][2] - lhs[isize][AA][0][2]*lhs[isize-1][CC][4][0]
				      - lhs[isize][AA][1][2]*lhs[isize-1][CC][4][1]
				      - lhs[isize][AA][2][2]*lhs[isize-1][CC][4][2]
				      - lhs[isize][AA][3][2]*lhs[isize-1][CC][4][3]
				      - lhs[isize][AA][4][2]*lhs[isize-1][CC][4][4];
	  lhs[isize][BB][4][3] = lhs[isize][BB][4][3] - lhs[isize][AA][0][3]*lhs[isize-1][CC][4][0]
				      - lhs[isize][AA][1][3]*lhs[isize-1][CC][4][1]
				      - lhs[isize][AA][2][3]*lhs[isize-1][CC][4][2]
				      - lhs[isize][AA][3][3]*lhs[isize-1][CC][4][3]
				      - lhs[isize][AA][4][3]*lhs[isize-1][CC][4][4];
	  lhs[isize][BB][4][4] = lhs[isize][BB][4][4] - lhs[isize][AA][0][4]*lhs[isize-1][CC][4][0]
				      - lhs[isize][AA][1][4]*lhs[isize-1][CC][4][1]
				      - lhs[isize][AA][2][4]*lhs[isize-1][CC][4][2]
				      - lhs[isize][AA][3][4]*lhs[isize-1][CC][4][3]
				      - lhs[isize][AA][4][4]*lhs[isize-1][CC][4][4];
      } // END of matmul_sub(lhs[isize][AA], lhs[isize-1][CC], lhs[isize][BB]);


      //---------------------------------------------------------------------
      // multiply rhs() by b_inverse() and copy to rhs
      //
      // binvrhs( lhs[isize][BB], rhs[k][j][isize] );
      // void binvrhs(double lhs[5][5], double r[5])
      //---------------------------------------------------------------------
      {
	  pivot3 = 1.00/lhs[isize][BB][0][0];
	  lhs[isize][BB][1][0] = lhs[isize][BB][1][0]*pivot3;
	  lhs[isize][BB][2][0] = lhs[isize][BB][2][0]*pivot3;
	  lhs[isize][BB][3][0] = lhs[isize][BB][3][0]*pivot3;
	  lhs[isize][BB][4][0] = lhs[isize][BB][4][0]*pivot3;
	  rhs[k][j][isize][0]   = rhs[k][j][isize][0]  *pivot3;

	  coeff3 = lhs[isize][BB][0][1];
	  lhs[isize][BB][1][1]= lhs[isize][BB][1][1] - coeff3*lhs[isize][BB][1][0];
	  lhs[isize][BB][2][1]= lhs[isize][BB][2][1] - coeff3*lhs[isize][BB][2][0];
	  lhs[isize][BB][3][1]= lhs[isize][BB][3][1] - coeff3*lhs[isize][BB][3][0];
	  lhs[isize][BB][4][1]= lhs[isize][BB][4][1] - coeff3*lhs[isize][BB][4][0];
	  rhs[k][j][isize][1]   = rhs[k][j][isize][1]   - coeff3*rhs[k][j][isize][0];

	  coeff3 = lhs[isize][BB][0][2];
	  lhs[isize][BB][1][2]= lhs[isize][BB][1][2] - coeff3*lhs[isize][BB][1][0];
	  lhs[isize][BB][2][2]= lhs[isize][BB][2][2] - coeff3*lhs[isize][BB][2][0];
	  lhs[isize][BB][3][2]= lhs[isize][BB][3][2] - coeff3*lhs[isize][BB][3][0];
	  lhs[isize][BB][4][2]= lhs[isize][BB][4][2] - coeff3*lhs[isize][BB][4][0];
	  rhs[k][j][isize][2]   = rhs[k][j][isize][2]   - coeff3*rhs[k][j][isize][0];

	  coeff3 = lhs[isize][BB][0][3];
	  lhs[isize][BB][1][3]= lhs[isize][BB][1][3] - coeff3*lhs[isize][BB][1][0];
	  lhs[isize][BB][2][3]= lhs[isize][BB][2][3] - coeff3*lhs[isize][BB][2][0];
	  lhs[isize][BB][3][3]= lhs[isize][BB][3][3] - coeff3*lhs[isize][BB][3][0];
	  lhs[isize][BB][4][3]= lhs[isize][BB][4][3] - coeff3*lhs[isize][BB][4][0];
	  rhs[k][j][isize][3]   = rhs[k][j][isize][3]   - coeff3*rhs[k][j][isize][0];

	  coeff3 = lhs[isize][BB][0][4];
	  lhs[isize][BB][1][4]= lhs[isize][BB][1][4] - coeff3*lhs[isize][BB][1][0];
	  lhs[isize][BB][2][4]= lhs[isize][BB][2][4] - coeff3*lhs[isize][BB][2][0];
	  lhs[isize][BB][3][4]= lhs[isize][BB][3][4] - coeff3*lhs[isize][BB][3][0];
	  lhs[isize][BB][4][4]= lhs[isize][BB][4][4] - coeff3*lhs[isize][BB][4][0];
	  rhs[k][j][isize][4]   = rhs[k][j][isize][4]   - coeff3*rhs[k][j][isize][0];


	  pivot3 = 1.00/lhs[isize][BB][1][1];
	  lhs[isize][BB][2][1] = lhs[isize][BB][2][1]*pivot3;
	  lhs[isize][BB][3][1] = lhs[isize][BB][3][1]*pivot3;
	  lhs[isize][BB][4][1] = lhs[isize][BB][4][1]*pivot3;
	  rhs[k][j][isize][1]   = rhs[k][j][isize][1]  *pivot3;

	  coeff3 = lhs[isize][BB][1][0];
	  lhs[isize][BB][2][0]= lhs[isize][BB][2][0] - coeff3*lhs[isize][BB][2][1];
	  lhs[isize][BB][3][0]= lhs[isize][BB][3][0] - coeff3*lhs[isize][BB][3][1];
	  lhs[isize][BB][4][0]= lhs[isize][BB][4][0] - coeff3*lhs[isize][BB][4][1];
	  rhs[k][j][isize][0]   = rhs[k][j][isize][0]   - coeff3*rhs[k][j][isize][1];

	  coeff3 = lhs[isize][BB][1][2];
	  lhs[isize][BB][2][2]= lhs[isize][BB][2][2] - coeff3*lhs[isize][BB][2][1];
	  lhs[isize][BB][3][2]= lhs[isize][BB][3][2] - coeff3*lhs[isize][BB][3][1];
	  lhs[isize][BB][4][2]= lhs[isize][BB][4][2] - coeff3*lhs[isize][BB][4][1];
	  rhs[k][j][isize][2]   = rhs[k][j][isize][2]   - coeff3*rhs[k][j][isize][1];

	  coeff3 = lhs[isize][BB][1][3];
	  lhs[isize][BB][2][3]= lhs[isize][BB][2][3] - coeff3*lhs[isize][BB][2][1];
	  lhs[isize][BB][3][3]= lhs[isize][BB][3][3] - coeff3*lhs[isize][BB][3][1];
	  lhs[isize][BB][4][3]= lhs[isize][BB][4][3] - coeff3*lhs[isize][BB][4][1];
	  rhs[k][j][isize][3]   = rhs[k][j][isize][3]   - coeff3*rhs[k][j][isize][1];

	  coeff3 = lhs[isize][BB][1][4];
	  lhs[isize][BB][2][4]= lhs[isize][BB][2][4] - coeff3*lhs[isize][BB][2][1];
	  lhs[isize][BB][3][4]= lhs[isize][BB][3][4] - coeff3*lhs[isize][BB][3][1];
	  lhs[isize][BB][4][4]= lhs[isize][BB][4][4] - coeff3*lhs[isize][BB][4][1];
	  rhs[k][j][isize][4]   = rhs[k][j][isize][4]   - coeff3*rhs[k][j][isize][1];


	  pivot3 = 1.00/lhs[isize][BB][2][2];
	  lhs[isize][BB][3][2] = lhs[isize][BB][3][2]*pivot3;
	  lhs[isize][BB][4][2] = lhs[isize][BB][4][2]*pivot3;
	  rhs[k][j][isize][2]   = rhs[k][j][isize][2]  *pivot3;

	  coeff3 = lhs[isize][BB][2][0];
	  lhs[isize][BB][3][0]= lhs[isize][BB][3][0] - coeff3*lhs[isize][BB][3][2];
	  lhs[isize][BB][4][0]= lhs[isize][BB][4][0] - coeff3*lhs[isize][BB][4][2];
	  rhs[k][j][isize][0]   = rhs[k][j][isize][0]   - coeff3*rhs[k][j][isize][2];

	  coeff3 = lhs[isize][BB][2][1];
	  lhs[isize][BB][3][1]= lhs[isize][BB][3][1] - coeff3*lhs[isize][BB][3][2];
	  lhs[isize][BB][4][1]= lhs[isize][BB][4][1] - coeff3*lhs[isize][BB][4][2];
	  rhs[k][j][isize][1]   = rhs[k][j][isize][1]   - coeff3*rhs[k][j][isize][2];

	  coeff3 = lhs[isize][BB][2][3];
	  lhs[isize][BB][3][3]= lhs[isize][BB][3][3] - coeff3*lhs[isize][BB][3][2];
	  lhs[isize][BB][4][3]= lhs[isize][BB][4][3] - coeff3*lhs[isize][BB][4][2];
	  rhs[k][j][isize][3]   = rhs[k][j][isize][3]   - coeff3*rhs[k][j][isize][2];

	  coeff3 = lhs[isize][BB][2][4];
	  lhs[isize][BB][3][4]= lhs[isize][BB][3][4] - coeff3*lhs[isize][BB][3][2];
	  lhs[isize][BB][4][4]= lhs[isize][BB][4][4] - coeff3*lhs[isize][BB][4][2];
	  rhs[k][j][isize][4]   = rhs[k][j][isize][4]   - coeff3*rhs[k][j][isize][2];


	  pivot3 = 1.00/lhs[isize][BB][3][3];
	  lhs[isize][BB][4][3] = lhs[isize][BB][4][3]*pivot3;
	  rhs[k][j][isize][3]   = rhs[k][j][isize][3]  *pivot3;

	  coeff3 = lhs[isize][BB][3][0];
	  lhs[isize][BB][4][0]= lhs[isize][BB][4][0] - coeff3*lhs[isize][BB][4][3];
	  rhs[k][j][isize][0]   = rhs[k][j][isize][0]   - coeff3*rhs[k][j][isize][3];

	  coeff3 = lhs[isize][BB][3][1];
	  lhs[isize][BB][4][1]= lhs[isize][BB][4][1] - coeff3*lhs[isize][BB][4][3];
	  rhs[k][j][isize][1]   = rhs[k][j][isize][1]   - coeff3*rhs[k][j][isize][3];

	  coeff3 = lhs[isize][BB][3][2];
	  lhs[isize][BB][4][2]= lhs[isize][BB][4][2] - coeff3*lhs[isize][BB][4][3];
	  rhs[k][j][isize][2]   = rhs[k][j][isize][2]   - coeff3*rhs[k][j][isize][3];

	  coeff3 = lhs[isize][BB][3][4];
	  lhs[isize][BB][4][4]= lhs[isize][BB][4][4] - coeff3*lhs[isize][BB][4][3];
	  rhs[k][j][isize][4]   = rhs[k][j][isize][4]   - coeff3*rhs[k][j][isize][3];


	  pivot3 = 1.00/lhs[isize][BB][4][4];
	  rhs[k][j][isize][4]   = rhs[k][j][isize][4]  *pivot3;

	  coeff3 = lhs[isize][BB][4][0];
	  rhs[k][j][isize][0]   = rhs[k][j][isize][0]   - coeff3*rhs[k][j][isize][4];

	  coeff3 = lhs[isize][BB][4][1];
	  rhs[k][j][isize][1]   = rhs[k][j][isize][1]   - coeff3*rhs[k][j][isize][4];

	  coeff3 = lhs[isize][BB][4][2];
	  rhs[k][j][isize][2]   = rhs[k][j][isize][2]   - coeff3*rhs[k][j][isize][4];

	  coeff3 = lhs[isize][BB][4][3];
	  rhs[k][j][isize][3]   = rhs[k][j][isize][3]   - coeff3*rhs[k][j][isize][4];
      }//END of binvrhs( lhs[isize][BB], rhs[k][j][isize] );


      //---------------------------------------------------------------------
      // back solve: if last cell, then generate U(isize)=rhs(isize)
      // else assume U(isize) is loaded in un pack backsub_info
      // so just use it
      // after u(istart) will be sent to next cell
      //---------------------------------------------------------------------
      for (i = isize-1; i >=0; i--) {
        for (m = 0; m < BLOCK_SIZE; m++) {
          for (n = 0; n < BLOCK_SIZE; n++) {
            rhs[k][j][i][m] = rhs[k][j][i][m] 
              - lhs[i][CC][n][m]*rhs[k][j][i+1][n];
          }
        }
      }
    }
  }
#pragma endscop  
}
