/*This is a short description of the parallel-radix sort algorithm.
  The list of numbers is divided into (parallel) groups.
  Each group has a separate bucket. Each bucket has 16 bins(0xF).
  The sorting is performed over chunks of each 4-bit (hex) digits.

  The radix-sort first counts (function: reduce) the numbers in its
  group that fall into different bins of its bucket.

  Next, the counts per-bin across all buckets are added. The addition
  is performed using parallel prefix-scan (top-scan & bottom-scan)
  algorithm. The scans leave the information about local-offsets.

  The sorting takes place (function: spread) by moving the numbers from
  old position to new position (global-offset + local-offset).
  The movement is performed for each hex digit (LSB first) and at the end
  of all the movements the numbers are in sorted order.
 */
#include "pencil.h"

static void reduce( int n, int group_size, int shift,
             unsigned int in[const restrict static n],
             unsigned int bucket[const restrict static no_of_groups][16])
{
#pragma scop
  {
  __pencil_assume(n > 0);
  __pencil_assume(n > no_of_groups);
  __pencil_assume(n > group_size);
  __pencil_assume(group_size >= 2);
  __pencil_assume(shift > 0);

  {
#pragma pencil independent
    for(int g = 0; g < no_of_groups; g++)
    {
      for(int k = 0; k < 16; k++) bucket[g][k] = 0;

      int start = g*group_size;
      int stop;
      if( g == no_of_groups-1) stop = n; else stop = (g+1)*group_size;

      for( int k = start; k < stop; k++)
      {
        unsigned int index;
        index = (in[k] >> shift) & 0x000F;
        bucket[g][index]++;
      }
    }
  }
  }
#pragma endscop
}

static void top_scan( int n, int group_size, int hdigit,
    unsigned int in[const restrict static n],
    unsigned int bucket[const restrict static no_of_groups][16])
{
  __pencil_assume((hdigit >=0) && (hdigit < 16));
  __pencil_assume( (n>0) && (n > no_of_groups) && (n>group_size));
  __pencil_assume( group_size >= 2);

  for(int d = 0; d < NGRP_LOG2; d++){
    int step = 1 << (d+1);
    int dual = 1 << (d);

#pragma scop
    {
#pragma pencil independent
      for(int k = 0; k < ((no_of_groups + step - 1) / step); k++){
          int step_index = k * step +step-1;
          int dual_index = k * step +dual-1;
          bucket[step_index][hdigit] += bucket[dual_index][hdigit];
      }
    }
#pragma endscop

  }
}

static void bottom_scan( int n, int group_size, int hdigit,
    unsigned int in[const restrict static n],
    unsigned int bucket[const restrict static no_of_groups][16])
{
  __pencil_assume(hdigit >=0);
  __pencil_assume(hdigit < 16);
  __pencil_assume(n>0);
  __pencil_assume(n > no_of_groups);
  __pencil_assume(n>group_size);
  __pencil_assume(group_size >= 2);

  bucket[no_of_groups-1][hdigit] = 0;

  for( int d = NGRP_LOG2-1; d >= 0 ; d--){
    int step = 1 << (d+1);
    int dual = 1 << (d);

#pragma scop
    {
#pragma pencil independent
      for(int k = 0; k < ((no_of_groups + step - 1) / step); k++){
          int dual_index = k * step +dual-1;
          int step_index = k * step +step-1;

          int t = bucket[dual_index][hdigit];
          bucket[dual_index][hdigit] = bucket[step_index][hdigit];
          bucket[step_index][hdigit] += t;
      }
    }
#pragma endscop
  }
}

static unsigned int compute_global_offset( unsigned int bucket[const restrict static no_of_groups][16],
    unsigned int global_offset[const restrict static 16])
{
#pragma scop
  {
#pragma pencil independent
    for(int d = 0; d < 16; d++){
      unsigned int sum = 0;
      for(int g = 0; g < no_of_groups; g++ ){
        sum += bucket[g][d];
      }
      global_offset[d] = sum;
    }
  }
#pragma endscop

  unsigned int sum = 0;
  for(int d = 0; d < 16; d++){
    unsigned int t = global_offset[d];
    global_offset[d] = sum;
    sum += t;
  }
  return sum;
}

static void spread( int n, int group_size, unsigned int shift,
    unsigned int in[const restrict static n], unsigned int out[const restrict static n],
    unsigned int bucket[const restrict static no_of_groups][16], unsigned int global_offset[const restrict static 16]);// __attribute__ (( __pencil_access(spread_summary)));

static void spread_summary( int n, int group_size, unsigned int shift,
    unsigned int in[const restrict static n], unsigned int out[const restrict static n],
    unsigned int bucket[const restrict static no_of_groups][16], unsigned int global_offset[const restrict static 16])
{
  for(int i = 0; i < 16; i++){
    __pencil_assume((global_offset[i] < n));
  }

  for(int i = 0; i < no_of_groups; i++){
    for(int j = 0; j < 16; j++){
      __pencil_assume((bucket[i][j] < n));
    }
  }
}

static void spread( int n, int group_size, unsigned int shift,
    unsigned int in[const restrict static n], unsigned int out[const restrict static n],
    unsigned int bucket[const restrict static no_of_groups][16], unsigned int global_offset[const restrict static 16])
{
  __pencil_assume( (n>0) && (n > no_of_groups) && (n>group_size));
  __pencil_assume( group_size >= 2);

#pragma scop
  {
#pragma pencil independent
    for(int g = 0; g < no_of_groups; g++)
    {
      int start = g*group_size;
      int stop;
      if( g == no_of_groups-1) stop = n; else stop = (g+1)*group_size;

      unsigned int local_bucket[16];
      for(int i = 0 ; i < 16; i++) local_bucket[i] = bucket[g][i];

      for(int k = start; k < stop; k++)
      {
        unsigned int digit = (in[k] >> shift) & 0x000F;
        unsigned int dst = global_offset[digit] + local_bucket[digit];

        __pencil_assume(digit <= 15);
        __pencil_assume((dst < n));
        __pencil_assume((k >= 0) && (k < n));

        local_bucket[digit]++;
        out[dst] = in[k];
      }
    }
  }
#pragma endscop
}

void parallel_sort_int(int n, unsigned int in[const restrict static 2][n])
{
  const int num_digits  = 16;//number of digits

  int t = n/no_of_groups;
  const int group_size = (t*no_of_groups == n) ? t : (t+1);

  unsigned int global_offset[16];
  unsigned int bucket[no_of_groups][16]; //ERROR --convert to malloc later

  //sort for each digit i.e. radix_width
  for(int shift = 0; shift < sizeof(unsigned int)*8; shift +=4){
    /*#pragma scop*/
    {
    int in_idx = (shift / 4) % 2;
    int out_idx = 1 - in_idx;
    //walk over the list and increment the counter for each bucket[grp][digit]
    reduce(n, group_size, shift, in[in_idx], bucket);

    unsigned int K = compute_global_offset(bucket, global_offset);

    //sum the elements falling into each bin (repeat for each digit)
    for(int digit = 0; digit < num_digits; digit++)
    {
      top_scan(n, group_size, digit, in[in_idx], bucket);
      bottom_scan(n, group_size, digit, in[in_idx], bucket);
    }

    //actual sorting based on the digits (put in the right place using scan infor)
    spread(n, group_size, shift, in[in_idx], in[out_idx], bucket, global_offset);
    }
    /*#pragma endscop*/

  }
}


