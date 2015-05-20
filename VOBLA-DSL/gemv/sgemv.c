#include "../include/vobla_headers.h"

#define N 1024
#define M 1024
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

	DATA_TYPE (*A)[M][N] = (DATA_TYPE (*)[M][N]) malloc(M*N * sizeof(DATA_TYPE));
	DATA_TYPE (*B)[N][1] = (DATA_TYPE (*)[N][1]) malloc(N*1 * sizeof(DATA_TYPE));
	DATA_TYPE (*C)[M][1] = (DATA_TYPE (*)[M][1]) malloc(M*1 * sizeof(DATA_TYPE));

        init_data(M, N, *A);
	init_data(N, 1, *B);
	init_data(M, 1, *C);


	start_time = get_time();
	pencil_sgemv_t_nn(M, N, 10.0, 0, *A, 1, *B, 10.0, 1, *C);
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

/*---- VOBLA function.  ----*/
void pencil_sgemv_t_nn (int m_9793_10025,int n_9794_10026,float alpha_9795_10027,int ldA_9796_10028,float A_9797_10029[restrict const static m_9793_10025][n_9794_10026],int incX_9798_10030,float X_9799_10031[restrict const static n_9794_10026][incX_9798_10030],float beta_9800_10032,int incY_9801_10033,float Y_9802_10034[restrict const static m_9793_10025][incY_9801_10033]){
#pragma scop
{
#pragma pencil independent
for(int i_2515_2299_5268_12053 = 0;i_2515_2299_5268_12053 <= (m_9793_10025-1);i_2515_2299_5268_12053 += 1){
Y_9802_10034[i_2515_2299_5268_12053][0] = (Y_9802_10034[i_2515_2299_5268_12053][0]*beta_9800_10032);
}
for(int i_2517_2300_5269_12055 = 0;i_2517_2300_5269_12055 <= (n_9794_10026-1);i_2517_2300_5269_12055 += 1){
for(int i_2516_2301_5270_12054 = 0;i_2516_2301_5270_12054 <= (m_9793_10025-1);i_2516_2301_5270_12054 += 1){
Y_9802_10034[i_2516_2301_5270_12054][0] = (Y_9802_10034[i_2516_2301_5270_12054][0]+((alpha_9795_10027*A_9797_10029[i_2517_2300_5269_12055][i_2516_2301_5270_12054])*X_9799_10031[i_2517_2300_5269_12055][0]));
}
}
}
#pragma endscop

}
