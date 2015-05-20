#include "../include/vobla_headers.h"

#define N 1024
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

	DATA_TYPE (*A)[N][K] = (DATA_TYPE (*)[N][K]) malloc(N*K * sizeof(DATA_TYPE));
	DATA_TYPE (*B)[K][N] = (DATA_TYPE (*)[K][N]) malloc(K*N * sizeof(DATA_TYPE));
	DATA_TYPE (*C)[N][N] = (DATA_TYPE (*)[N][N]) malloc(N*N * sizeof(DATA_TYPE));

        init_data(N, K, *A);
        init_data(K, N, *A);
	init_data(N, N, *C);


	start_time = get_time();
	pencil_ssyr2k_nn_u(N, K, 10.0, 1, *A, 1, *B, 10.0, 1, *C);
	end_time = get_time();

#ifdef BENCHMARK_TIME
	fprintf(TIME_OUTPUT_FILE, "%0.6lf", end_time - start_time);
	fprintf(TIME_OUTPUT_FILE, "\n");
#endif

#ifdef BENCHMARK_DUMP_ARRAYS
	print_array(N, N, *C);
#endif
	return 0;
}

/*---- VOBLA function.  ----*/
void pencil_ssyr2k_nn_u (int n_115_555,int k_116_556,float alpha_117_557,int ldA_118_558,float A_119_559[restrict const static n_115_555][k_116_556],int ldB_120_560,float B_121_561[restrict const static k_116_556][n_115_555],float beta_122_562,int ldCT_123_563,float CT_124_564[restrict const static n_115_555][n_115_555]){

#pragma scop
{
#pragma pencil independent
for(int i_236_382_1467_2001 = 0;i_236_382_1467_2001 <= (n_115_555-1);i_236_382_1467_2001 += 1){
#pragma pencil independent
for(int i_238_384_1468_2003 = 0;i_238_384_1468_2003 <= ((((n_115_555-1)-i_236_382_1467_2001)+1)-1);i_238_384_1468_2003 += 1){
float accum_239_385_1469_2000;
accum_239_385_1469_2000 = 0.0f;
for(int i_240_386_1470_2002 = 0;i_240_386_1470_2002 <= (k_116_556-1);i_240_386_1470_2002 += 1){
accum_239_385_1469_2000 = (accum_239_385_1469_2000+((A_119_559[i_236_382_1467_2001][i_240_386_1470_2002]*B_121_561[(i_238_384_1468_2003+i_236_382_1467_2001)][i_240_386_1470_2002])+(B_121_561[i_236_382_1467_2001][i_240_386_1470_2002]*A_119_559[(i_238_384_1468_2003+i_236_382_1467_2001)][i_240_386_1470_2002])));
}
CT_124_564[i_236_382_1467_2001][(i_236_382_1467_2001+i_238_384_1468_2003)] = ((beta_122_562*CT_124_564[i_236_382_1467_2001][(i_236_382_1467_2001+i_238_384_1468_2003)])+(alpha_117_557*accum_239_385_1469_2000));
}
}
}
#pragma endscop

}
