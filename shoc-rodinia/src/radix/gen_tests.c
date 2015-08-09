//Generates a list of random numbers for sorting
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <stdlib.h>
#include <limits.h>

int main(int argc, char *argv[]){
  if( argc < 2 ){
    printf("Usage: %s number \n", argv[0]); return 1;
  }

  int N = atoi(argv[1]);

  char file_name[100];
  sprintf(file_name,"n%d",N);

  FILE *file = fopen(file_name,"w");
  if( file == NULL){
    printf("%s: unable to open file %s \n", argv[0], file_name); return 1;
  }

  srand( time(NULL));
  fprintf(file,"%d\n",N);
  for(int i = 0; i < N; i++ ) fprintf(file, "%d\n", rand() % INT_MAX);
  fclose(file);
}
