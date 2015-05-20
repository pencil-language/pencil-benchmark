#include "../include/vobla_headers.h"

#define N 1024
#define M 1024
#define K 1024
#define DATA_TYPE float
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

	DATA_TYPE (*A)[M][K] = (DATA_TYPE (*)[M][K]) malloc(M*K * sizeof(DATA_TYPE));
	DATA_TYPE (*B)[K][N] = (DATA_TYPE (*)[K][N]) malloc(K*N * sizeof(DATA_TYPE));
	DATA_TYPE (*C)[M][N] = (DATA_TYPE (*)[M][N]) malloc(M*N * sizeof(DATA_TYPE));

        init_data(M, K, *A);
	init_data(K, N, *B);
	init_data(M, N, *C);


	/* pencil_saxpy_nn */
	start_time = get_time();
	pencil_sgemm_ttt(M, N, K, 10, 0, *A, 0, *B, 10, 0, *C);
	end_time = get_time();

#ifdef BENCHMARK_TIME
	fprintf(TIME_OUTPUT_FILE, "%0.6lf", end_time - start_time);
	fprintf(TIME_OUTPUT_FILE, "\n");
#endif

#ifdef BENCHMARK_DUMP_ARRAYS
	print_array(M, N, *C);
#endif
	return 0;
}

/*---- VOBLA function.  ----*/
void pencil_sgemm_ttt (int m_3219_3413,int n_3220_3414,int k_3221_3415,float alpha_3222_3416,int ldA_3223_3417,float A_3224_3418[restrict const static m_3219_3413][k_3221_3415],int ldB_3225_3419,float B_3226_3420[restrict const static k_3221_3415][n_3220_3414],float beta_3227_3421,int ldC_3228_3422,float C_3229_3423[restrict const static m_3219_3413][n_3220_3414]){

#pragma scop
{
#pragma pencil independent
for(int i_689_1002_2260_4927 = 0;i_689_1002_2260_4927 <= (n_3220_3414-1);i_689_1002_2260_4927 += 1){
#pragma pencil independent
for(int i_688_1003_2261_4928 = 0;i_688_1003_2261_4928 <= (m_3219_3413-1);i_688_1003_2261_4928 += 1){
C_3229_3423[i_689_1002_2260_4927][i_688_1003_2261_4928] = (C_3229_3423[i_689_1002_2260_4927][i_688_1003_2261_4928]*beta_3227_3421);
}
}
for(int i_691_1004_2262_4931 = 0;i_691_1004_2262_4931 <= (k_3221_3415-1);i_691_1004_2262_4931 += 1){
for(int i_690_1005_2263_4930 = 0;i_690_1005_2263_4930 <= (m_3219_3413-1);i_690_1005_2263_4930 += 1){
for(int i_692_1006_2264_4929 = 0;i_692_1006_2264_4929 <= (n_3220_3414-1);i_692_1006_2264_4929 += 1){
C_3229_3423[i_692_1006_2264_4929][i_690_1005_2263_4930] = (C_3229_3423[i_692_1006_2264_4929][i_690_1005_2263_4930]+((alpha_3222_3416*A_3224_3418[i_691_1004_2262_4931][i_690_1005_2263_4930])*B_3226_3420[i_692_1006_2264_4929][i_691_1004_2262_4931]));
}
}
}
}
#pragma endscop

}
