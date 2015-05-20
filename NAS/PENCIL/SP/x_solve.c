//-------------------------------------------------------------------------//
//                                                                         //
//  This benchmark is an OpenMP C version of the NPB SP code. This OpenMP  //
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

// Simplified threadprivate extracted from header.h
//#pragma omp threadprivate(cv,rhon,lhs,lhsp,lhsm)

    // The live range on the array cv[] is iteration private for K and for J
    // loops because the first i loop writes in the whole array and the 2nd i
    // loop reads the whole array.  cv[] is not used anywhere later.	
    // The same with rhon[].
    // cv[] and rhon[]should be privatized with ompen private().  If we care
    // only about tiling, then the TACO paper can ignore the fdeps and the code
    // can be tiles.  To parallelize the code we need to privatize the two
    // arrays with private().
    // lhs[][][] is iteration private for the loop k, it can be privatized with
    // private() also to enable parallelism.
    // For finding that the whole array is iteration private you need just to
    // the difference between source and sink for the k dimension and you'll
    // find it equanl to Zero so you conclude that the array is iteration
    // private for that dimension.
    // the fact that lhs[][][] is written in a loop and is read in the
    // following loop and that it is written in that loop and read in the
    // following loop makes the whole program in One SCC for the k dimension
    // (we have an sc graph for each loop dimension).
    // SCCs are described in details in the document
    // Desktop/stmt_clustering.txt
    // Full expansion of lhs[][][] will enable:
    // -- parallelization of the K outermost loop,
    // -- loop distribution, you can distribute the k loop and have one k loop that writes into
    // lhs[][][][] for all k values, then the second k loop will consume all
    // the values written by the first k loop.  Loop distribution in this case
    // does not enhance things, it does not enable any advanced loop nest
    // transformation, so full expansion is not really very interesting.
    // So we just need privatization to the parallelization of the code.
    //
    // The arrays lhsp[][][] and lhsm[][][] are like lhs[][][], the are
    // iteration private to the k loop
    // 
    // The array rhs[][][][] is the array in which this huge x_solve.c function puts the results.
    // The output of the huge K loop is stored in this array.  It does not need
    // any privatization or exapansion.

    // In the original code, we had and irregular array accesses (look at the
    // loops by the end of the SCoP in the original SP):
    // i = grid_points[0];
    // i1 = grid_points[0];
    // lhs[j][i][*] = ....
    // lhs[j][i1][*] = ...
    // Since lhs[][][] is iteration private.  The privatization will guarantee
    // that different iterations access different lines of lhs[][][] which
    // enable the parallelization of the loop although the original loop
    // contained an irregular write access.

    
    
//---------------------------------------------------------------------
// this function performs the solution of the approximate factorization
// step in the x-direction for all five matrix components
// simultaneously. The Thomas algorithm is employed to solve the
// systems for the x-lines. Boundary conditions are non-periodic
//---------------------------------------------------------------------
void x_solve()
{
  int i, j, k, i1, i2, m;
  double ru1, fac1, fac2;

//  #pragma omp parallel for default(shared) private(i,j,k,i1,i2,m, \
                                                   ru1,fac1,fac2)
#pragma scop							   
  for (k = 1; k <=  PROBLEM_SIZE - 2; k++) {
// 	lhsinit(PROBLEM_SIZE - 2+1, ny2);
//	void lhsinit(int ni, int nj)
	{
	  //---------------------------------------------------------------------
	  // zap the whole left hand side for starters
	  // set all diagonal values to 1. This is overkill, but convenient
	  //---------------------------------------------------------------------
	  for (j = 1; j <= PROBLEM_SIZE - 2; j++) {
	    for (m = 0; m < 5; m++) {
	      lhs [j][0][m] = 0.0;
	      lhsp[j][0][m] = 0.0;
	      lhsm[j][0][m] = 0.0;
	      lhs [j][PROBLEM_SIZE - 2+1][m] = 0.0;
	      lhsp[j][PROBLEM_SIZE - 2+1][m] = 0.0;
	      lhsm[j][PROBLEM_SIZE - 2+1][m] = 0.0;
	    }
	    lhs [j][0][2] = 1.0;
	    lhsp[j][0][2] = 1.0;
	    lhsm[j][0][2] = 1.0;
	    lhs [j][PROBLEM_SIZE - 2+1][2] = 1.0;
	    lhsp[j][PROBLEM_SIZE - 2+1][2] = 1.0;
	    lhsm[j][PROBLEM_SIZE - 2+1][2] = 1.0;
	  }
	}

    //---------------------------------------------------------------------
    // Computes the left hand side for the three x-factors  
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    // first fill the lhs for the u-eigenvalue                   
    //---------------------------------------------------------------------
    for (j = 1; j <= PROBLEM_SIZE - 2; j++) {
      for (i = 0; i <= PROBLEM_SIZE-1; i++) {
        ru1 = c3c4*rho_i[k][j][i];
        cv[i] = us[k][j][i];
        rhon[i] = max(max(dx2+con43*ru1,dx5+c1c5*ru1), max(dxmax+ru1,dx1));
      }

      for (i = 1; i <= PROBLEM_SIZE - 2; i++) {
        lhs[j][i][0] =  0.0;
        lhs[j][i][1] = -dttx2 * cv[i-1] - dttx1 * rhon[i-1];
        lhs[j][i][2] =  1.0 + c2dttx1 * rhon[i];
        lhs[j][i][3] =  dttx2 * cv[i+1] - dttx1 * rhon[i+1];
        lhs[j][i][4] =  0.0;
      }
    }

    //---------------------------------------------------------------------
    // add fourth order dissipation                             
    //---------------------------------------------------------------------
    for (j = 1; j <= PROBLEM_SIZE - 2; j++) {
      lhs[j][1][2] = lhs[j][1][2] + comz5;
      lhs[j][1][3] = lhs[j][1][3] - comz4;
      lhs[j][1][4] = lhs[j][1][4] + comz1;

      lhs[j][1+1][1] = lhs[j][1+1][1] - comz4;
      lhs[j][1+1][2] = lhs[j][1+1][2] + comz6;
      lhs[j][1+1][3] = lhs[j][1+1][3] - comz4;
      lhs[j][1+1][4] = lhs[j][1+1][4] + comz1;
    }

    for (j = 1; j <= PROBLEM_SIZE - 2; j++) {
      for (i = 3; i <= PROBLEM_SIZE-4; i++) {
        lhs[j][i][0] = lhs[j][i][0] + comz1;
        lhs[j][i][1] = lhs[j][i][1] - comz4;
        lhs[j][i][2] = lhs[j][i][2] + comz6;
        lhs[j][i][3] = lhs[j][i][3] - comz4;
        lhs[j][i][4] = lhs[j][i][4] + comz1;
      }
    }

    for (j = 1; j <= PROBLEM_SIZE - 2; j++) {
      lhs[j][ PROBLEM_SIZE-3][0] = lhs[j][ PROBLEM_SIZE-3][0] + comz1;
      lhs[j][ PROBLEM_SIZE-3][1] = lhs[j][ PROBLEM_SIZE-3][1] - comz4;
      lhs[j][ PROBLEM_SIZE-3][2] = lhs[j][ PROBLEM_SIZE-3][2] + comz6;
      lhs[j][ PROBLEM_SIZE-3][3] = lhs[j][ PROBLEM_SIZE-3][3] - comz4;

      lhs[j][ PROBLEM_SIZE-3+1][0] = lhs[j][ PROBLEM_SIZE-3+1][0] + comz1;
      lhs[j][ PROBLEM_SIZE-3+1][1] = lhs[j][ PROBLEM_SIZE-3+1][1] - comz4;
      lhs[j][ PROBLEM_SIZE-3+1][2] = lhs[j][ PROBLEM_SIZE-3+1][2] + comz5;
    }

    //---------------------------------------------------------------------
    // subsequently, fill the other factors (u+c), (u-c) by adding to 
    // the first  
    //---------------------------------------------------------------------
    for (j = 1; j <= PROBLEM_SIZE - 2; j++) {
      for (i = 1; i <= PROBLEM_SIZE - 2; i++) {
        lhsp[j][i][0] = lhs[j][i][0];
        lhsp[j][i][1] = lhs[j][i][1] - dttx2 * speed[k][j][i-1];
        lhsp[j][i][2] = lhs[j][i][2];
        lhsp[j][i][3] = lhs[j][i][3] + dttx2 * speed[k][j][i+1];
        lhsp[j][i][4] = lhs[j][i][4];
        lhsm[j][i][0] = lhs[j][i][0];
        lhsm[j][i][1] = lhs[j][i][1] + dttx2 * speed[k][j][i-1];
        lhsm[j][i][2] = lhs[j][i][2];
        lhsm[j][i][3] = lhs[j][i][3] - dttx2 * speed[k][j][i+1];
        lhsm[j][i][4] = lhs[j][i][4];
      }
    }

    //---------------------------------------------------------------------
    // FORWARD ELIMINATION  
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    // perform the Thomas algorithm; first, FORWARD ELIMINATION     
    //---------------------------------------------------------------------
    for (j = 1; j <= PROBLEM_SIZE - 2; j++) {
      for (i = 0; i <= PROBLEM_SIZE-3; i++) {
        fac1 = 1.0/lhs[j][i][2];
        lhs[j][i][3] = fac1*lhs[j][i][3];
        lhs[j][i][4] = fac1*lhs[j][i][4];
        for (m = 0; m < 3; m++) {
          rhs[k][j][i][m] = fac1*rhs[k][j][i][m];
        }
        lhs[j][i+1][2] = lhs[j][i+1][2] - lhs[j][i+1][1]*lhs[j][i][3];
        lhs[j][i+1][3] = lhs[j][i+1][3] - lhs[j][i+1][1]*lhs[j][i][4];
        for (m = 0; m < 3; m++) {
          rhs[k][j][i+1][m] = rhs[k][j][i+1][m] - lhs[j][i+1][1]*rhs[k][j][i][m];
        }
        lhs[j][i+2][1] = lhs[j][i+2][1] - lhs[j][i+2][0]*lhs[j][i][3];
        lhs[j][i+2][2] = lhs[j][i+2][2] - lhs[j][i+2][0]*lhs[j][i][4];
        for (m = 0; m < 3; m++) {
          rhs[k][j][i+2][m] = rhs[k][j][i+2][m] - lhs[j][i+2][0]*rhs[k][j][i][m];
        }
      }
    }

    //---------------------------------------------------------------------
    // The last two rows in this grid block are a bit different,
    // since they for (not have two more rows available for the
    // elimination of off-diagonal entries
    //---------------------------------------------------------------------
    for (j = 1; j <= PROBLEM_SIZE - 2; j++) {
      fac1 = 1.0/lhs[j][PROBLEM_SIZE-2][2];
      lhs[j][PROBLEM_SIZE-2][3] = fac1*lhs[j][PROBLEM_SIZE-2][3];
      lhs[j][PROBLEM_SIZE-2][4] = fac1*lhs[j][PROBLEM_SIZE-2][4];
      for (m = 0; m < 3; m++) {
        rhs[k][j][PROBLEM_SIZE-2][m] = fac1*rhs[k][j][PROBLEM_SIZE-2][m];
      }
      lhs[j][PROBLEM_SIZE-1][2] = lhs[j][PROBLEM_SIZE-1][2] - lhs[j][PROBLEM_SIZE-1][1]*lhs[j][PROBLEM_SIZE-2][3];
      lhs[j][PROBLEM_SIZE-1][3] = lhs[j][PROBLEM_SIZE-1][3] - lhs[j][PROBLEM_SIZE-1][1]*lhs[j][PROBLEM_SIZE-2][4];
      for (m = 0; m < 3; m++) {
        rhs[k][j][PROBLEM_SIZE-1][m] = rhs[k][j][PROBLEM_SIZE-1][m] - lhs[j][PROBLEM_SIZE-1][1]*rhs[k][j][PROBLEM_SIZE-2][m];
      }

      //---------------------------------------------------------------------
      // scale the last row immediately 
      //---------------------------------------------------------------------
      fac2 = 1.0/lhs[j][PROBLEM_SIZE-1][2];
      for (m = 0; m < 3; m++) {
        rhs[k][j][PROBLEM_SIZE-1][m] = fac2*rhs[k][j][PROBLEM_SIZE-1][m];
      }
    }

    //---------------------------------------------------------------------
    // for (the u+c and the u-c factors                 
    //---------------------------------------------------------------------
    for (j = 1; j <= PROBLEM_SIZE - 2; j++) {
      for (i = 0; i <= PROBLEM_SIZE-3; i++) {
        fac1 = 1.0/lhsp[j][i][2];
        lhsp[j][i][3]    = fac1*lhsp[j][i][3];
        lhsp[j][i][4]    = fac1*lhsp[j][i][4];
        rhs[k][j][i][3]  = fac1*rhs[k][j][i][3];
        lhsp[j][i+1][2]   = lhsp[j][i+1][2] - lhsp[j][i+1][1]*lhsp[j][i][3];
        lhsp[j][i+1][3]   = lhsp[j][i+1][3] - lhsp[j][i+1][1]*lhsp[j][i][4];
        rhs[k][j][i+1][3] = rhs[k][j][i+1][3] - lhsp[j][i+1][1]*rhs[k][j][i][3];
        lhsp[j][i+2][1]   = lhsp[j][i+2][1] - lhsp[j][i+2][0]*lhsp[j][i][3];
        lhsp[j][i+2][2]   = lhsp[j][i+2][2] - lhsp[j][i+2][0]*lhsp[j][i][4];
        rhs[k][j][i+2][3] = rhs[k][j][i+2][3] - lhsp[j][i+2][0]*rhs[k][j][i][3];

        fac1 = 1.0/lhsm[j][i][2];
        lhsm[j][i][3]    = fac1*lhsm[j][i][3];
        lhsm[j][i][4]    = fac1*lhsm[j][i][4];
        rhs[k][j][i][4]  = fac1*rhs[k][j][i][4];
        lhsm[j][i+1][2]   = lhsm[j][i+1][2] - lhsm[j][i+1][1]*lhsm[j][i][3];
        lhsm[j][i+1][3]   = lhsm[j][i+1][3] - lhsm[j][i+1][1]*lhsm[j][i][4];
        rhs[k][j][i+1][4] = rhs[k][j][i+1][4] - lhsm[j][i+1][1]*rhs[k][j][i][4];
        lhsm[j][i+2][1]   = lhsm[j][i+2][1] - lhsm[j][i+2][0]*lhsm[j][i][3];
        lhsm[j][i+2][2]   = lhsm[j][i+2][2] - lhsm[j][i+2][0]*lhsm[j][i][4];
        rhs[k][j][i+2][4] = rhs[k][j][i+2][4] - lhsm[j][i+2][0]*rhs[k][j][i][4];
      }
    }

    //---------------------------------------------------------------------
    // And again the last two rows separately
    //---------------------------------------------------------------------
    for (j = 1; j <= PROBLEM_SIZE - 2; j++) {
      fac1 = 1.0/lhsp[j][PROBLEM_SIZE-2][2];
      lhsp[j][PROBLEM_SIZE-2][3]    = fac1*lhsp[j][PROBLEM_SIZE-2][3];
      lhsp[j][PROBLEM_SIZE-2][4]    = fac1*lhsp[j][PROBLEM_SIZE-2][4];
      rhs[k][j][PROBLEM_SIZE-2][3]  = fac1*rhs[k][j][PROBLEM_SIZE-2][3];
      lhsp[j][PROBLEM_SIZE-1][2]   = lhsp[j][PROBLEM_SIZE-1][2] - lhsp[j][PROBLEM_SIZE-1][1]*lhsp[j][PROBLEM_SIZE-2][3];
      lhsp[j][PROBLEM_SIZE-1][3]   = lhsp[j][PROBLEM_SIZE-1][3] - lhsp[j][PROBLEM_SIZE-1][1]*lhsp[j][PROBLEM_SIZE-2][4];
      rhs[k][j][PROBLEM_SIZE-1][3] = rhs[k][j][PROBLEM_SIZE-1][3] - lhsp[j][PROBLEM_SIZE-1][1]*rhs[k][j][PROBLEM_SIZE-2][3];

      fac1 = 1.0/lhsm[j][PROBLEM_SIZE-2][2];
      lhsm[j][PROBLEM_SIZE-2][3]    = fac1*lhsm[j][PROBLEM_SIZE-2][3];
      lhsm[j][PROBLEM_SIZE-2][4]    = fac1*lhsm[j][PROBLEM_SIZE-2][4];
      rhs[k][j][PROBLEM_SIZE-2][4]  = fac1*rhs[k][j][PROBLEM_SIZE-2][4];
      lhsm[j][PROBLEM_SIZE-1][2]   = lhsm[j][PROBLEM_SIZE-1][2] - lhsm[j][PROBLEM_SIZE-1][1]*lhsm[j][PROBLEM_SIZE-2][3];
      lhsm[j][PROBLEM_SIZE-1][3]   = lhsm[j][PROBLEM_SIZE-1][3] - lhsm[j][PROBLEM_SIZE-1][1]*lhsm[j][PROBLEM_SIZE-2][4];
      rhs[k][j][PROBLEM_SIZE-1][4] = rhs[k][j][PROBLEM_SIZE-1][4] - lhsm[j][PROBLEM_SIZE-1][1]*rhs[k][j][PROBLEM_SIZE-2][4];

      //---------------------------------------------------------------------
      // Scale the last row immediately
      //---------------------------------------------------------------------
      rhs[k][j][PROBLEM_SIZE-1][3] = rhs[k][j][PROBLEM_SIZE-1][3]/lhsp[j][PROBLEM_SIZE-1][2];
      rhs[k][j][PROBLEM_SIZE-1][4] = rhs[k][j][PROBLEM_SIZE-1][4]/lhsm[j][PROBLEM_SIZE-1][2];
    }

    //---------------------------------------------------------------------
    // BACKSUBSTITUTION 
    //---------------------------------------------------------------------
    for (j = 1; j <= PROBLEM_SIZE - 2; j++) {
      for (m = 0; m < 3; m++) {
        rhs[k][j][PROBLEM_SIZE-2][m] = rhs[k][j][PROBLEM_SIZE-2][m] - lhs[j][PROBLEM_SIZE-2][3]*rhs[k][j][PROBLEM_SIZE-1][m];
      }

      rhs[k][j][PROBLEM_SIZE-2][3] = rhs[k][j][PROBLEM_SIZE-2][3] - lhsp[j][PROBLEM_SIZE-2][3]*rhs[k][j][PROBLEM_SIZE-1][3];
      rhs[k][j][PROBLEM_SIZE-2][4] = rhs[k][j][PROBLEM_SIZE-2][4] - lhsm[j][PROBLEM_SIZE-2][3]*rhs[k][j][PROBLEM_SIZE-1][4];
    }

    //---------------------------------------------------------------------
    // The first three factors
    //---------------------------------------------------------------------
    for (j = 1; j <= PROBLEM_SIZE - 2; j++) {
      for (i = PROBLEM_SIZE-3; i >= 0; i--) {
        for (m = 0; m < 3; m++) {
          rhs[k][j][i][m] = rhs[k][j][i][m] - 
                            lhs[j][i][3]*rhs[k][j][i+1][m] -
                            lhs[j][i][4]*rhs[k][j][i+2][m];
        }

        //-------------------------------------------------------------------
        // And the remaining two
        //-------------------------------------------------------------------
        rhs[k][j][i][3] = rhs[k][j][i][3] - 
                          lhsp[j][i][3]*rhs[k][j][i+1][3] -
                          lhsp[j][i][4]*rhs[k][j][i+2][3];
        rhs[k][j][i][4] = rhs[k][j][i][4] - 
                          lhsm[j][i][3]*rhs[k][j][i+1][4] -
                          lhsm[j][i][4]*rhs[k][j][i+2][4];
      }
    }
  }
#pragma endscop

  //---------------------------------------------------------------------
  // Do the block-diagonal inversion          
  //---------------------------------------------------------------------
  ninvr();
}

