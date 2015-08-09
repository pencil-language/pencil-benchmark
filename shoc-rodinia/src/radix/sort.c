#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "pencil_runtime.h"

#include "measure-time.h"

typedef unsigned int uint;
#define no_of_groups 32

void parallel_sort_int(int n, unsigned int in[const restrict static 2][n]);

//debug-printing functions
void print_list( int n, uint in[n], char *msg){
  printf("%s:\n", msg);
  for( int i = 0; i < n; i++ ) printf("%X ",in[i]); printf("\n");
}

void print_bucket(uint bucket[no_of_groups][16], char *msg, int x){
  printf("%s %d:\n", msg,x);

  printf("bucket[grp=X]: ");
  for( int d = 0; d < 16; d++ ) printf("%2X ", d);
  printf("\n");

  for( int g = 0; g < no_of_groups; g++ ){
    printf("bucket[grp=%d]: ",g);
    for(int k = 0; k < 16; k++ ){
      printf("%2d ", bucket[g][k]);
    }
    printf("\n");
  }
}

void print_global_offset(uint global_offset[16]){
  printf("global offset: ");
  for( int d = 0; d < 16; d++ ) printf("%2d ", global_offset[d]);
  printf("\n");
}

//simple serial sort - keeps input data intact
void serial_sort( int N, uint keys[N], uint sorted_keys[N])
{
  for(int i = 0; i < N; i++) sorted_keys[i] = keys[i];
  int swapped = 1;
  while(swapped){
    swapped = 0;
    for(int i = 0; i < N-1; i++){
      if( sorted_keys[i] > sorted_keys[i+1] ){
        uint t = sorted_keys[i];
        sorted_keys[i] = sorted_keys[i+1];
        sorted_keys[i+1] = t;
        swapped = 1;
      }
    }
  }
}

int verify( int N, uint keys[N], uint sorted_keys[N])
{
  for(int i = 0; i < N-1; i++) assert( keys[i] <= keys[i+1]);
  for(int i = 0; i < N; i++)   assert( keys[i] == sorted_keys[i]);
  printf("SUCCESS \n");
  return 0;
}


void parallel_sort( int n, unsigned int in[const restrict static n])
{
  unsigned int *buff = malloc(sizeof(unsigned int) * n * 2);
  memcpy(buff, in, n * sizeof(unsigned int));
  parallel_sort_int(n, buff);
  memcpy(in, buff, n * sizeof(unsigned int));
  free(buff);
}


int main(int argc, char *argv[]){
  int check = 0;
  if( argc < 2 ){
    printf("Usage: %s data_file \n", argv[0]); return 1;
  }

  if( argc > 2 && !strcmp(argv[2], "--check"))
  {
    check = 1;
  }

  FILE *file = fopen(argv[1],"r");
  if( file == NULL){
    printf("%s: unable to open file %s \n", argv[0], argv[1]); return 1;
  }

  int N;
  fscanf(file,"%d", &N);
  assert(N > 0);

  pencil_init(PENCIL_TARGET_DEVICE_DYNAMIC);

  uint *keys         = (uint *)malloc(sizeof(uint)*N);
  uint *ssorted_keys = (uint *)malloc(sizeof(uint)*N);
  uint *psorted_keys = (uint *)pencil_alloc(sizeof(uint)*N);

  for(int i = 0; i < N; i++ ) fscanf(file, "%d", &keys[i]);
  for(int i = 0; i < N; i++ ) ssorted_keys[i] = psorted_keys[i] = keys[i];
  START_MEASURE_TIME
  parallel_sort(N, psorted_keys);
  STOP_MEASURE_TIME

  if (check)
  {
    serial_sort(N, keys, ssorted_keys);
    verify(N, ssorted_keys, psorted_keys);
  }

  free(keys);
  free(ssorted_keys);
  pencil_free(psorted_keys);
  pencil_shutdown();
  PRINT_TIME
  return 0;
}
