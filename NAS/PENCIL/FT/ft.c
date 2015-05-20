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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "global.h"
#include "randdp.h"
#include "timers.h"
#include "print_results.h"

//---------------------------------------------------------------------------
static dcomplex ty1[MAXDIM][FFTBLOCKPAD_DEFAULT];
static dcomplex ty2[MAXDIM][FFTBLOCKPAD_DEFAULT];
#pragma omp threadprivate(ty1,ty2)

dcomplex dcmplx_mul_fct(dcomplex u1, dcomplex x11, dcomplex x21)
{
	return dcmplx_mul(u1, dcmplx_sub(x11, x21));
}

//---------------------------------------------------------------------
// u0, u1, u2 are the main arrays in the problem. 
// Depending on the decomposition, these arrays will have different 
// dimensions. To accomodate all possibilities, we allocate them as 
// one-dimensional arrays and pass them to subroutines for different 
// views
//  - u0 contains the initial (transformed) initial condition
//  - u1 and u2 are working arrays
//  - twiddle contains exponents for the time evolution operator. 
//---------------------------------------------------------------------
//---------------------------------------------------------------------
// Large arrays are in common so that they are allocated on the
// heap rather than the stack. This common block is not
// referenced directly anywhere else. Padding is to avoid accidental 
// cache problems, since all array sizes are powers of two.
//---------------------------------------------------------------------
/* common /bigarrays/ */
dcomplex u0[NTOTALP];
dcomplex pad1[3];
dcomplex u1[NTOTALP];
dcomplex pad2[3];
//dcomplex u2[NTOTALP];
double twiddle[NTOTALP];
//---------------------------------------------------------------------------

#include "ft.pencil.c"

int main(int argc, char *argv[])
{
  int i;
  int iter;
  double total_time, mflops;
  logical verified;
  char Class;

  //---------------------------------------------------------------------
  // Run the entire problem once to make sure all data is touched. 
  // This reduces variable startup costs, which is important for such a 
  // short benchmark. The other NPB 2 implementations are similar. 
  //---------------------------------------------------------------------
  for (i = 1; i <= T_max; i++) {
    timer_clear(i);
  }
  setup();
  init_ui(u0, u1, twiddle, dims[0], dims[1], dims[2]);
  compute_indexmap(dims[0], dims[1], dims[2], twiddle);
  compute_initial_conditions(u1, dims[0], dims[1], dims[2]);
  fft_init(dims[0]);
  fft(1, dims[0], dims[1], dims[2], u1, u0);

  //---------------------------------------------------------------------
  // Start over from the beginning. Note that all operations must
  // be timed, in contrast to other benchmarks. 
  //---------------------------------------------------------------------
  for (i = 1; i <= T_max; i++) {
    timer_clear(i);
  }

  timer_start(T_total);
  if (timers_enabled) timer_start(T_setup);

  compute_indexmap(dims[0], dims[1], dims[2], twiddle);

  compute_initial_conditions(u1, dims[0], dims[1], dims[2]);

  fft_init(dims[0]);

  if (timers_enabled) timer_stop(T_setup);
  if (timers_enabled) timer_start(T_fft);
  fft(1, dims[0], dims[1], dims[2], u1, u0);
  if (timers_enabled) timer_stop(T_fft);

  for (iter = 1; iter <= niter; iter++) {
    if (timers_enabled) timer_start(T_evolve);
    evolve(dims[0], dims[1], dims[2], u0, u1, twiddle);
    if (timers_enabled) timer_stop(T_evolve);
    if (timers_enabled) timer_start(T_fft);
    fft(-1, dims[0], dims[1], dims[2], u1, u1);
    if (timers_enabled) timer_stop(T_fft);
    if (timers_enabled) timer_start(T_checksum);
    checksum(iter, u1, dims[0], dims[1], dims[2]);
    if (timers_enabled) timer_stop(T_checksum);
  }

  verify(NX, NY, NZ, niter, &verified, &Class);

  timer_stop(T_total);
  total_time = timer_read(T_total);

  if (total_time != 0.0) {
    mflops = 1.0e-6 * (double)NTOTAL *
            (14.8157 + 7.19641 * log((double)NTOTAL)
            + (5.23518 + 7.21113 * log((double)NTOTAL)) * niter)
            / total_time;
  } else {
    mflops = 0.0;
  }
  print_results("FT", Class, NX, NY, NZ, niter,
                total_time, mflops, "          floating point", verified, 
                NPBVERSION, COMPILETIME, CS1, CS2, CS3, CS4, CS5, CS6, CS7);
  if (timers_enabled) print_timers();

  return 0;
}


//---------------------------------------------------------------------
// touch all the big data
//---------------------------------------------------------------------
static void init_ui(void *ou0, void *ou1, void *ot, int d1, int d2, int d3)
{
  dcomplex (*u0)[d2][d1+1] = (dcomplex (*)[d2][d1+1])ou0;
  dcomplex (*u1)[d2][d1+1] = (dcomplex (*)[d2][d1+1])ou1;
  double (*twiddle)[d2][d1+1] = (double (*)[d2][d1+1])ot;

  int i, j, k;

  #pragma omp parallel for default(shared) private(i,j,k)
  for (k = 0; k < d3; k++) {
    for (j = 0; j < d2; j++) {
      for (i = 0; i < d1; i++) {
        u0[k][j][i] = dcmplx(0.0, 0.0);
        u1[k][j][i] = dcmplx(0.0, 0.0);
        twiddle[k][j][i] = 0.0;
      }
    }
  }
}

//---------------------------------------------------------------------
// Fill in array u0 with initial conditions from 
// random number generator 
//---------------------------------------------------------------------
static void compute_initial_conditions(void *ou0, int d1, int d2, int d3)
{
  dcomplex (*u0)[d2][d1+1] = (dcomplex (*)[d2][d1+1])ou0;

  int k, j;
  double x0, start, an, dummy, starts[NZ];

  start = SEED;
  //---------------------------------------------------------------------
  // Jump to the starting element for our first plane.
  //---------------------------------------------------------------------
  an = ipow46(A, 0);
  dummy = randlc(&start, an);
  an = ipow46(A, 2*NX*NY);

  starts[0] = start;
  for (k = 1; k < dims[2]; k++) {
    dummy = randlc(&start, an);
    starts[k] = start;
  }

  //---------------------------------------------------------------------
  // Go through by z planes filling in one square at a time.
  //---------------------------------------------------------------------
  #pragma omp parallel for default(shared) private(k,j,x0)
  for (k = 0; k < dims[2]; k++) {
    x0 = starts[k];
    for (j = 0; j < dims[1]; j++) {
      vranlc(2*NX, &x0, A, (double *)&u0[k][j][0]);
    }
  }
}


//---------------------------------------------------------------------
// compute a^exponent mod 2^46
//---------------------------------------------------------------------
static double ipow46(double a, int exponent)
{
  double result, dummy, q, r;
  int n, n2;

  //---------------------------------------------------------------------
  // Use
  //   a^n = a^(n/2)*a^(n/2) if n even else
  //   a^n = a*a^(n-1)       if n odd
  //---------------------------------------------------------------------
  result = 1;
  if (exponent == 0) return result;
  q = a;
  r = 1;
  n = exponent;

  while (n > 1) {
    n2 = n / 2;
    if (n2 * 2 == n) {
      dummy = randlc(&q, q);
      n = n2;
    } else {
      dummy = randlc(&r, q);
      n = n-1;
    }
  }
  dummy = randlc(&r, q);
  result = r;
  return result;
}


static void setup()
{
  FILE *fp;
  debug = false;

  if ((fp = fopen("timer.flag", "r")) != NULL) {
    timers_enabled = true;
    fclose(fp);
  } else {
    timers_enabled = false;
  }

  niter = NITER_DEFAULT;

  printf("\n\n NAS Parallel Benchmarks (NPB3.3-OMP-C) - FT Benchmark\n\n");
  printf(" Size                : %4dx%4dx%4d\n", NX, NY, NZ);
  printf(" Iterations                  :%7d\n", niter);
  printf(" Number of available threads :%7d\n", omp_get_max_threads());
  printf("\n");

  dims[0] = NX;
  dims[1] = NY;
  dims[2] = NZ;

  //---------------------------------------------------------------------
  // Set up info for blocking of ffts and transposes.  This improves
  // performance on cache-based systems. Blocking involves
  // working on a chunk of the problem at a time, taking chunks
  // along the first, second, or third dimension. 
  //
  // - In cffts1 blocking is on 2nd dimension (with fft on 1st dim)
  // - In cffts2/3 blocking is on 1st dimension (with fft on 2nd and 3rd dims)

  // Since 1st dim is always in processor, we'll assume it's long enough 
  // (default blocking factor is 16 so min size for 1st dim is 16)
  // The only case we have to worry about is cffts1 in a 2d decomposition. 
  // so the blocking factor should not be larger than the 2nd dimension. 
  //---------------------------------------------------------------------

  fftblock = FFTBLOCK_DEFAULT;
  fftblockpad = FFTBLOCKPAD_DEFAULT;

  if (fftblock != FFTBLOCK_DEFAULT) fftblockpad = fftblock+3;
}

static void print_timers()
{
  int i;
  double t, t_m;
  char *tstrings[T_max+1];
  tstrings[1] = "          total "; 
  tstrings[2] = "          setup "; 
  tstrings[3] = "            fft "; 
  tstrings[4] = "         evolve "; 
  tstrings[5] = "       checksum "; 
  tstrings[6] = "           fftx "; 
  tstrings[7] = "           ffty "; 
  tstrings[8] = "           fftz ";

  t_m = timer_read(T_total);
  if (t_m <= 0.0) t_m = 1.00;
  for (i = 1; i <= T_max; i++) {
    t = timer_read(i);
    printf(" timer %2d(%16s) :%9.4f (%6.2f%%)\n", 
        i, tstrings[i], t, t*100.0/t_m);
  }
}


//---------------------------------------------------------------------
// compute the roots-of-unity array that will be used for subsequent FFTs. 
//---------------------------------------------------------------------
static void fft_init(int n)
{
  int m, nu, ku, i, j, ln;
  double t, ti;

  //---------------------------------------------------------------------
  // Initialize the U array with sines and cosines in a manner that permits
  // stride one access at each FFT iteration.
  //---------------------------------------------------------------------
  nu = n;
  m = ilog2(n);
  u[0] = dcmplx(m, 0.0);
  ku = 2;
  ln = 1;

  for (j = 1; j <= m; j++) {
    t = PI / ln;

    for (i = 0; i <= ln - 1; i++) {
      ti = i * t;
      u[i+ku-1] = dcmplx(cos(ti), sin(ti));
    }

    ku = ku + ln;
    ln = 2 * ln;
  }
}

static void checksum(int i, void *ou1, int d1, int d2, int d3)
{
  dcomplex (*u1)[d2][d1+1] = (dcomplex (*)[d2][d1+1])ou1;

  int j, q, r, s;
  dcomplex chk = dcmplx(0.0, 0.0);

  #pragma omp parallel default(shared) private(i,q,r,s)
  {
    dcomplex my_chk = dcmplx(0.0, 0.0);

    #pragma omp for nowait
    for (j = 1; j <= 1024; j++) {
      q = j % NX;
      r = 3*j % NY;
      s = 5*j % NZ;
      my_chk = dcmplx_add(my_chk, u1[s][r][q]);
    }

    #pragma omp critical
    {
      chk = dcmplx_add(chk, my_chk);
    }
  }

  chk = dcmplx_div2(chk, (double)(NTOTAL));

  printf(" T =%5d     Checksum =%22.12E%22.12E\n", i, chk.real, chk.imag);
  sums[i] = chk;
}


static void verify(int d1, int d2, int d3, int nt, 
                   logical *verified, char *Class)
{
  int i;
  double err, epsilon;

  //---------------------------------------------------------------------
  // Reference checksums
  //---------------------------------------------------------------------
  dcomplex csum_ref[25+1];

  *Class = 'U';

  epsilon = 1.0e-12;
  *verified = false;

  if (d1 == 64 && d2 == 64 && d3 == 64 && nt == 6) {
    //---------------------------------------------------------------------
    //   Sample size reference checksums
    //---------------------------------------------------------------------
    *Class = 'S';
    csum_ref[1] = dcmplx(5.546087004964E+02, 4.845363331978E+02);
    csum_ref[2] = dcmplx(5.546385409189E+02, 4.865304269511E+02);
    csum_ref[3] = dcmplx(5.546148406171E+02, 4.883910722336E+02);
    csum_ref[4] = dcmplx(5.545423607415E+02, 4.901273169046E+02);
    csum_ref[5] = dcmplx(5.544255039624E+02, 4.917475857993E+02);
    csum_ref[6] = dcmplx(5.542683411902E+02, 4.932597244941E+02);

  } else if (d1 == 128 && d2 == 128 && d3 == 32 && nt == 6) {
    //---------------------------------------------------------------------
    //   Class W size reference checksums
    //---------------------------------------------------------------------
    *Class = 'W';
    csum_ref[1] = dcmplx(5.673612178944E+02, 5.293246849175E+02);
    csum_ref[2] = dcmplx(5.631436885271E+02, 5.282149986629E+02);
    csum_ref[3] = dcmplx(5.594024089970E+02, 5.270996558037E+02);
    csum_ref[4] = dcmplx(5.560698047020E+02, 5.260027904925E+02);
    csum_ref[5] = dcmplx(5.530898991250E+02, 5.249400845633E+02);
    csum_ref[6] = dcmplx(5.504159734538E+02, 5.239212247086E+02);

  } else if (d1 == 256 && d2 == 256 && d3 == 128 && nt == 6) {
    //---------------------------------------------------------------------
    //   Class A size reference checksums
    //---------------------------------------------------------------------
    *Class = 'A';
    csum_ref[1] = dcmplx(5.046735008193E+02, 5.114047905510E+02);
    csum_ref[2] = dcmplx(5.059412319734E+02, 5.098809666433E+02);
    csum_ref[3] = dcmplx(5.069376896287E+02, 5.098144042213E+02);
    csum_ref[4] = dcmplx(5.077892868474E+02, 5.101336130759E+02);
    csum_ref[5] = dcmplx(5.085233095391E+02, 5.104914655194E+02);
    csum_ref[6] = dcmplx(5.091487099959E+02, 5.107917842803E+02);

  } else if (d1 == 512 && d2 == 256 && d3 == 256 && nt == 20) {
    //---------------------------------------------------------------------
    //   Class B size reference checksums
    //---------------------------------------------------------------------
    *Class = 'B';
    csum_ref[1]  = dcmplx(5.177643571579E+02, 5.077803458597E+02);
    csum_ref[2]  = dcmplx(5.154521291263E+02, 5.088249431599E+02);
    csum_ref[3]  = dcmplx(5.146409228649E+02, 5.096208912659E+02);
    csum_ref[4]  = dcmplx(5.142378756213E+02, 5.101023387619E+02);
    csum_ref[5]  = dcmplx(5.139626667737E+02, 5.103976610617E+02);
    csum_ref[6]  = dcmplx(5.137423460082E+02, 5.105948019802E+02);
    csum_ref[7]  = dcmplx(5.135547056878E+02, 5.107404165783E+02);
    csum_ref[8]  = dcmplx(5.133910925466E+02, 5.108576573661E+02);
    csum_ref[9]  = dcmplx(5.132470705390E+02, 5.109577278523E+02);
    csum_ref[10] = dcmplx(5.131197729984E+02, 5.110460304483E+02);
    csum_ref[11] = dcmplx(5.130070319283E+02, 5.111252433800E+02);
    csum_ref[12] = dcmplx(5.129070537032E+02, 5.111968077718E+02);
    csum_ref[13] = dcmplx(5.128182883502E+02, 5.112616233064E+02);
    csum_ref[14] = dcmplx(5.127393733383E+02, 5.113203605551E+02);
    csum_ref[15] = dcmplx(5.126691062020E+02, 5.113735928093E+02);
    csum_ref[16] = dcmplx(5.126064276004E+02, 5.114218460548E+02);
    csum_ref[17] = dcmplx(5.125504076570E+02, 5.114656139760E+02);
    csum_ref[18] = dcmplx(5.125002331720E+02, 5.115053595966E+02);
    csum_ref[19] = dcmplx(5.124551951846E+02, 5.115415130407E+02);
    csum_ref[20] = dcmplx(5.124146770029E+02, 5.115744692211E+02);

  } else if (d1 == 512 && d2 == 512 && d3 == 512 && nt == 20) {
    //---------------------------------------------------------------------
    //   Class C size reference checksums
    //---------------------------------------------------------------------
    *Class = 'C';
    csum_ref[1]  = dcmplx(5.195078707457E+02, 5.149019699238E+02);
    csum_ref[2]  = dcmplx(5.155422171134E+02, 5.127578201997E+02);
    csum_ref[3]  = dcmplx(5.144678022222E+02, 5.122251847514E+02);
    csum_ref[4]  = dcmplx(5.140150594328E+02, 5.121090289018E+02);
    csum_ref[5]  = dcmplx(5.137550426810E+02, 5.121143685824E+02);
    csum_ref[6]  = dcmplx(5.135811056728E+02, 5.121496764568E+02);
    csum_ref[7]  = dcmplx(5.134569343165E+02, 5.121870921893E+02);
    csum_ref[8]  = dcmplx(5.133651975661E+02, 5.122193250322E+02);
    csum_ref[9]  = dcmplx(5.132955192805E+02, 5.122454735794E+02);
    csum_ref[10] = dcmplx(5.132410471738E+02, 5.122663649603E+02);
    csum_ref[11] = dcmplx(5.131971141679E+02, 5.122830879827E+02);
    csum_ref[12] = dcmplx(5.131605205716E+02, 5.122965869718E+02);
    csum_ref[13] = dcmplx(5.131290734194E+02, 5.123075927445E+02);
    csum_ref[14] = dcmplx(5.131012720314E+02, 5.123166486553E+02);
    csum_ref[15] = dcmplx(5.130760908195E+02, 5.123241541685E+02);
    csum_ref[16] = dcmplx(5.130528295923E+02, 5.123304037599E+02);
    csum_ref[17] = dcmplx(5.130310107773E+02, 5.123356167976E+02);
    csum_ref[18] = dcmplx(5.130103090133E+02, 5.123399592211E+02);
    csum_ref[19] = dcmplx(5.129905029333E+02, 5.123435588985E+02);
    csum_ref[20] = dcmplx(5.129714421109E+02, 5.123465164008E+02);

  } else if (d1 == 2048 && d2 == 1024 && d3 == 1024 && nt == 25) {
    //---------------------------------------------------------------------
    //   Class D size reference checksums
    //---------------------------------------------------------------------
    *Class = 'D';
    csum_ref[1]  = dcmplx(5.122230065252E+02, 5.118534037109E+02);
    csum_ref[2]  = dcmplx(5.120463975765E+02, 5.117061181082E+02);
    csum_ref[3]  = dcmplx(5.119865766760E+02, 5.117096364601E+02);
    csum_ref[4]  = dcmplx(5.119518799488E+02, 5.117373863950E+02);
    csum_ref[5]  = dcmplx(5.119269088223E+02, 5.117680347632E+02);
    csum_ref[6]  = dcmplx(5.119082416858E+02, 5.117967875532E+02);
    csum_ref[7]  = dcmplx(5.118943814638E+02, 5.118225281841E+02);
    csum_ref[8]  = dcmplx(5.118842385057E+02, 5.118451629348E+02);
    csum_ref[9]  = dcmplx(5.118769435632E+02, 5.118649119387E+02);
    csum_ref[10] = dcmplx(5.118718203448E+02, 5.118820803844E+02);
    csum_ref[11] = dcmplx(5.118683569061E+02, 5.118969781011E+02);
    csum_ref[12] = dcmplx(5.118661708593E+02, 5.119098918835E+02);
    csum_ref[13] = dcmplx(5.118649768950E+02, 5.119210777066E+02);
    csum_ref[14] = dcmplx(5.118645605626E+02, 5.119307604484E+02);
    csum_ref[15] = dcmplx(5.118647586618E+02, 5.119391362671E+02);
    csum_ref[16] = dcmplx(5.118654451572E+02, 5.119463757241E+02);
    csum_ref[17] = dcmplx(5.118665212451E+02, 5.119526269238E+02);
    csum_ref[18] = dcmplx(5.118679083821E+02, 5.119580184108E+02);
    csum_ref[19] = dcmplx(5.118695433664E+02, 5.119626617538E+02);
    csum_ref[20] = dcmplx(5.118713748264E+02, 5.119666538138E+02);
    csum_ref[21] = dcmplx(5.118733606701E+02, 5.119700787219E+02);
    csum_ref[22] = dcmplx(5.118754661974E+02, 5.119730095953E+02);
    csum_ref[23] = dcmplx(5.118776626738E+02, 5.119755100241E+02);
    csum_ref[24] = dcmplx(5.118799262314E+02, 5.119776353561E+02);
    csum_ref[25] = dcmplx(5.118822370068E+02, 5.119794338060E+02);

  } else if (d1 == 4096 && d2 == 2048 && d3 == 2048 && nt == 25) {
    //---------------------------------------------------------------------
    //   Class E size reference checksums
    //---------------------------------------------------------------------
    *Class = 'E';
    csum_ref[1]  = dcmplx(5.121601045346E+02, 5.117395998266E+02);
    csum_ref[2]  = dcmplx(5.120905403678E+02, 5.118614716182E+02);
    csum_ref[3]  = dcmplx(5.120623229306E+02, 5.119074203747E+02);
    csum_ref[4]  = dcmplx(5.120438418997E+02, 5.119345900733E+02);
    csum_ref[5]  = dcmplx(5.120311521872E+02, 5.119551325550E+02);
    csum_ref[6]  = dcmplx(5.120226088809E+02, 5.119720179919E+02);
    csum_ref[7]  = dcmplx(5.120169296534E+02, 5.119861371665E+02);
    csum_ref[8]  = dcmplx(5.120131225172E+02, 5.119979364402E+02);
    csum_ref[9]  = dcmplx(5.120104767108E+02, 5.120077674092E+02);
    csum_ref[10] = dcmplx(5.120085127969E+02, 5.120159443121E+02);
    csum_ref[11] = dcmplx(5.120069224127E+02, 5.120227453670E+02);
    csum_ref[12] = dcmplx(5.120055158164E+02, 5.120284096041E+02);
    csum_ref[13] = dcmplx(5.120041820159E+02, 5.120331373793E+02);
    csum_ref[14] = dcmplx(5.120028605402E+02, 5.120370938679E+02);
    csum_ref[15] = dcmplx(5.120015223011E+02, 5.120404138831E+02);
    csum_ref[16] = dcmplx(5.120001570022E+02, 5.120432068837E+02);
    csum_ref[17] = dcmplx(5.119987650555E+02, 5.120455615860E+02);
    csum_ref[18] = dcmplx(5.119973525091E+02, 5.120475499442E+02);
    csum_ref[19] = dcmplx(5.119959279472E+02, 5.120492304629E+02);
    csum_ref[20] = dcmplx(5.119945006558E+02, 5.120506508902E+02);
    csum_ref[21] = dcmplx(5.119930795911E+02, 5.120518503782E+02);
    csum_ref[22] = dcmplx(5.119916728462E+02, 5.120528612016E+02);
    csum_ref[23] = dcmplx(5.119902874185E+02, 5.120537101195E+02);
    csum_ref[24] = dcmplx(5.119889291565E+02, 5.120544194514E+02);
    csum_ref[25] = dcmplx(5.119876028049E+02, 5.120550079284E+02);
  }

  if (*Class != 'U') {
    *verified = true;
    for (i = 1; i <= nt; i++) {
      err = dcmplx_abs(dcmplx_div(dcmplx_sub(sums[i], csum_ref[i]),
                                  csum_ref[i]));
      if (!(err <= epsilon)) {
        *verified = false;
        break;
      }
    }
  }

  if (*Class != 'U') {
    if (*verified) {
      printf(" Result verification successful\n");
    } else {
      printf(" Result verification failed\n");
    }
  }
  printf(" class = %c\n", *Class);
}

