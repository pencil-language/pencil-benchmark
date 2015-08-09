#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifndef HOST_DEBUG
#include "pencil_runtime.h"
#else
#define pencil_alloc malloc
#define pencil_free free
#define pencil_init()
#define pencil_shutdown()
#endif

#include "measure-time.h"

void diffusion(int Nr, int Nc, int Ne, float q0sqr,
    int iN[const restrict static Nr], int iS[const restrict static Nr],
    int jE[const restrict static Nc], int jW[const restrict static Nc],
    float dN[const restrict static Ne], float dS[const restrict static Ne],
    float dE[const restrict static Ne], float dW[const restrict static Ne],
    float image[const restrict static Ne], float c[const restrict static Ne]);

void divergence(int Nr, int Nc, int Ne, float lambda,
    int iS[const restrict static Nr], int jE[const restrict static Nc],
    float dN[const restrict static Ne], float dS[const restrict static Ne],
    float dE[const restrict static Ne], float dW[const restrict static Ne],
    float image[const restrict static Ne], float c[const restrict static Ne]);

void extract(int Ne, float image[const restrict static Ne]);
void compress(int Ne, float image[const restrict static Ne]);

void srad(int niter, int NeROI, int Nr, int Nc, int Ne, float lambda,
    int iN[const restrict static Nr], int iS[const restrict static Nr],
    int jE[const restrict static Nc], int jW[const restrict static Nc],
    float dN[const restrict static Ne], float dS[const restrict static Ne],
    float dE[const restrict static Ne], float dW[const restrict static Ne],
    float image[const restrict static Ne], float c[const restrict static Ne]);

/////////////////////////////////////////////////////////////////////////////////////////////////

void resize(float* input_data, int input_rows, int input_cols, float* output_data, int output_rows, int output_cols, int format){

  int i, j, i2, j2;

  if(format == 0){/* do if data is saved row major*/

    for(i=0, i2=0; i<output_rows; i++, i2++){
      if(i2>=input_rows){
        i2 = i2 - input_rows;
      }
      for(j=0, j2=0; j<output_cols; j++, j2++){
        if(j2>=input_cols){
          j2 = j2 - input_cols;
        }
        output_data[i*output_cols+j] = input_data[i2*input_cols+j2];
      }
    }
  }
  else{/* do if data is saved column major*/
    for(j=0, j2=0; j<output_cols; j++, j2++){
      if(j2>=input_cols){
        j2 = j2 - input_cols;
      }
      for(i=0, i2=0; i<output_rows; i++, i2++){
        if(i2>=input_rows){
          i2 = i2 - input_rows;
        }
        output_data[j*output_rows+i] = input_data[j2*input_rows+i2];
      }
    }
  }
}

void write_pgm_file(char* filename, float* input_data, int data_rows, int data_cols, int format, int data_range){
  FILE* fid;
  int i, j;

  fid = fopen(filename, "w");
  if( fid == NULL ){
    printf( "The file %s could not be created for writing\n", filename);
    return;
  }

  fprintf(fid, "P2\n");
  fprintf(fid, "%d %d\n", data_cols, data_rows);
  fprintf(fid, "%d\n", data_range);

  // if matrix is saved row major in memory (C)
  if(format==0){
    for(i=0; i<data_rows; i++){
      for(j=0; j<data_cols; j++){
        fprintf(fid, "%d ", (int)input_data[i*data_cols+j]);
      }
      fprintf(fid, "\n");
    }
  }
  // if matrix is saved column major in memory (MATLAB)
  else{
    for(i=0; i<data_rows; i++){
      for(j=0; j<data_cols; j++){
        fprintf(fid, "%d ", (int)input_data[j*data_rows+i]);
      }
      fprintf(fid, "\n");
    }
  }

  fclose(fid);
}

void read_pgm_file(const char* filename, float* input_data, int data_rows, int data_cols, int format){

  FILE* fid;
  int i, j, temp;
  char c;

  fid = fopen(filename, "r");
  if( fid == NULL ){
    printf( "Unable to open file %s for reading\n", filename);
    return;
  }

  i = 0;
  while(i<3){
    c = fgetc(fid);
    if(c == '\n'){
      i = i+1;
    }
  }

  if(format==0){// if matrix is saved row major in memory (C)
    for(i=0; i<data_rows; i++){
      for(j=0; j<data_cols; j++){
        fscanf(fid, "%d", &temp);
        input_data[i*data_cols+j] = (float)temp;
      }
    }
  }
  else{// if matrix is saved column major in memory (MATLAB)
    for(i=0; i<data_rows; i++){
      for(j=0; j<data_cols; j++){
        fscanf(fid, "%d", &temp);
        input_data[j*data_rows+i] = (float)temp;
      }
    }
  }
  fclose(fid);
}

void fix_surrounding_index(int Nr, int Nc,
    int iN[const restrict static Nr], int iS[const restrict static Nr],
    int jE[const restrict static Nc], int jW[const restrict static Nc]){

  for (int i=0; i<Nr; i++){
    iN[i] = i-1;/*index of pixel above*/
    iS[i] = i+1;/*index of pixel below*/
  }
  for (int j=0; j<Nc; j++){
    jW[j] = j-1;/*pixel to left*/
    jE[j] = j+1;/*pixel to right*/
  }

  /*boundary conditions*/
  iN[0]    = 0;
  iS[Nr-1] = Nr-1;
  jW[0]    = 0;
  jE[Nc-1] = Nc-1;
}


int main(int argc, char *argv []){

  if(argc != 6){
    printf("Usage: %s filename.pgm nIterations lambda nRows nColumns\n", argv[0]);
    return EXIT_FAILURE;
  }

  const char *file_name = argv[1];
  int   niter = atoi(argv[2]);
  float lambda = atof(argv[3]);
  int   Nr = atoi(argv[4]);/* it is 502 in the original image*/
  int   Nc = atoi(argv[5]);/* it is 458 in the original image*/

  int   Ne = Nr*Nc;
  int   NeROI = Nr*Nc;// number of elements in ROI, ROI size

  pencil_init(PENCIL_TARGET_DEVICE_DYNAMIC);

  /*surrounding pixel indices*/
  int *iN = pencil_alloc(sizeof(int*)*Nr);
  int *iS = pencil_alloc(sizeof(int*)*Nr);
  int *jW = pencil_alloc(sizeof(int*)*Nc);
  int *jE = pencil_alloc(sizeof(int*)*Nc);

  float *image = (float*)pencil_alloc(sizeof(float) * Ne);

  int original_Nr = 502;
  int original_Nc = 458;
  float* original_image = (float*)malloc(sizeof(float) * original_Nr * original_Nc);

  read_pgm_file(file_name, original_image, original_Nr, original_Nc, 1);/*assume col-major for original*/

  resize(original_image, original_Nr, original_Nc, image, Nr, Nc, 1);

  fix_surrounding_index(Nr, Nc, iN, iS, jE, jW);

  // Actual work starts here
  START_MEASURE_TIME

  /*directional derivatives*/
  float *dN = pencil_alloc(sizeof(float)*Ne);
  float *dS = pencil_alloc(sizeof(float)*Ne);
  float *dW = pencil_alloc(sizeof(float)*Ne);
  float *dE = pencil_alloc(sizeof(float)*Ne);

  float *c  = pencil_alloc(sizeof(float)*Ne) ;  /* diffusion coefficient*/

  srad(niter, NeROI, Nr, Nc, Ne, lambda, iN, iS, jE, jW, dN, dS, dE, dW, image, c);

  STOP_MEASURE_TIME

  write_pgm_file(	"image_out.pgm", image, Nr, Nc, 1, 255);
  printf("generated output file: image_out.pgm\n");

  pencil_free(iN); pencil_free(iS); pencil_free(jW); pencil_free(jE); pencil_free(dN); pencil_free(dS); pencil_free(dW); pencil_free(dE); pencil_free(c);
  pencil_free(image);
  free(original_image);
  pencil_shutdown();

  PRINT_TIME

  return EXIT_SUCCESS;
}
