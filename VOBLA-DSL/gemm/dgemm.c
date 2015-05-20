#include "../include/vobla_headers.h"

#define N 1024
#define M 1024
#define K 1024
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

	DATA_TYPE (*A)[M][K] = (DATA_TYPE (*)[M][K]) malloc(M*K * sizeof(DATA_TYPE));
	DATA_TYPE (*B)[K][N] = (DATA_TYPE (*)[K][N]) malloc(K*N * sizeof(DATA_TYPE));
	DATA_TYPE (*C)[M][N] = (DATA_TYPE (*)[M][N]) malloc(M*N * sizeof(DATA_TYPE));

        init_data(M, K, *A);
	init_data(K, N, *B);
	init_data(M, N, *C);


	/* pencil_saxpy_nn */
	start_time = get_time();
	pencil_dgemm_ttt(M, N, K, 10, 0, *A, 0, *B, 10, 0, *C);
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
void pencil_dgemm_ttt (int m_3367_3469,int n_3368_3470,int k_3369_3471,double alpha_3370_3472,int ldA_3371_3473,double A_3372_3474[restrict const static m_3367_3469][k_3369_3471],int ldB_3373_3475,double B_3374_3476[restrict const static k_3369_3471][n_3368_3470],double beta_3375_3477,int ldC_3376_3478,double C_3377_3479[restrict const static m_3367_3469][n_3368_3470]){
#pragma scop
{
#pragma pencil independent
for(int i_1003_972_2194_5003 = 0;i_1003_972_2194_5003 <= (n_3368_3470-1);i_1003_972_2194_5003 += 1){
#pragma pencil independent
for(int i_1002_973_2195_5004 = 0;i_1002_973_2195_5004 <= (m_3367_3469-1);i_1002_973_2195_5004 += 1){
C_3377_3479[i_1003_972_2194_5003][i_1002_973_2195_5004] = (C_3377_3479[i_1003_972_2194_5003][i_1002_973_2195_5004]*beta_3375_3477);
}
}
for(int i_1005_974_2196_5007 = 0;i_1005_974_2196_5007 <= (k_3369_3471-1);i_1005_974_2196_5007 += 1){
for(int i_1004_975_2197_5006 = 0;i_1004_975_2197_5006 <= (m_3367_3469-1);i_1004_975_2197_5006 += 1){
for(int i_1006_976_2198_5005 = 0;i_1006_976_2198_5005 <= (n_3368_3470-1);i_1006_976_2198_5005 += 1){
C_3377_3479[i_1006_976_2198_5005][i_1004_975_2197_5006] = (C_3377_3479[i_1006_976_2198_5005][i_1004_975_2197_5006]+((alpha_3370_3472*A_3372_3474[i_1005_974_2196_5007][i_1004_975_2197_5006])*B_3374_3476[i_1006_976_2198_5005][i_1005_974_2196_5007]));
}
}
}
}
#pragma endscop
}
