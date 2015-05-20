#include "../include/vobla_headers.h"

#define N 1024
#define M 1024
#define DATA_TYPE double
#define DATA_PRINTF_FORMAT "%f "

void print_array(int n, int m, DATA_TYPE B[n][m])
{
	int i, j;
	for (i=0; i<n; i++)
  	  for (j=0; j<m; j++)
		fprintf(DATA_OUTPUT_FILE, DATA_PRINTF_FORMAT, B[i][j]);
}

void init_data(int n, int m, DATA_TYPE B[n][m])
{
	int i, j;
	for (i=0; i<n; i++)
   	  for (j=0; j<m; j++)
		B[i][j] = 2.0;
}

int main()
{
	double start_time, end_time;

	DATA_TYPE (*A)[M][N] = (DATA_TYPE (*)[M][N]) malloc(M*N * sizeof(DATA_TYPE));
	DATA_TYPE (*B)[N][1] = (DATA_TYPE (*)[N][1]) malloc(N*1 * sizeof(DATA_TYPE));
	DATA_TYPE (*C)[M][1] = (DATA_TYPE (*)[M][1]) malloc(M*1 * sizeof(DATA_TYPE));

        init_data(M, N, *A);
	init_data(N, 1, *B);
	init_data(M, 1, *C);


	start_time = get_time();
	pencil_dgemv_t_nn(M, N, 10.0, 0, *A, 1, *B, 10.0, 1, *C);
	end_time = get_time();

#ifdef BENCHMARK_TIME
	fprintf(TIME_OUTPUT_FILE, "%0.6lf", end_time - start_time);
	fprintf(TIME_OUTPUT_FILE, "\n");
#endif

#ifdef BENCHMARK_DUMP_ARRAYS
	print_array(M, 1, *C);
#endif
	return 0;
}

void pencil_dgemv_t_nn (int m_10065_10129,int n_10066_10130,double alpha_10067_10131,int ldA_10068_10132,double A_10069_10133[restrict const static m_10065_10129][n_10066_10130],int incX_10070_10134,double X_10071_10135[restrict const static n_10066_10130][incX_10070_10134],double beta_10072_10136,int incY_10073_10137,double Y_10074_10138[restrict const static m_10065_10129][incY_10073_10137]){

#pragma scop
{
#pragma pencil independent
for(int i_4193_4429_8620_12181 = 0;i_4193_4429_8620_12181 <= (m_10065_10129-1);i_4193_4429_8620_12181 += 1){
Y_10074_10138[i_4193_4429_8620_12181][0] = (Y_10074_10138[i_4193_4429_8620_12181][0]*beta_10072_10136);
}
for(int i_4195_4430_8621_12183 = 0;i_4195_4430_8621_12183 <= (n_10066_10130-1);i_4195_4430_8621_12183 += 1){
for(int i_4194_4431_8622_12182 = 0;i_4194_4431_8622_12182 <= (m_10065_10129-1);i_4194_4431_8622_12182 += 1){
Y_10074_10138[i_4194_4431_8622_12182][0] = (Y_10074_10138[i_4194_4431_8622_12182][0]+((alpha_10067_10131*A_10069_10133[i_4195_4430_8621_12183][i_4194_4431_8622_12182])*X_10071_10135[i_4195_4430_8621_12183][0]));
}
}
}
#pragma endscop

}
