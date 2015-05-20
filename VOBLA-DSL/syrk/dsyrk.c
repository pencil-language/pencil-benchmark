#include "../include/vobla_headers.h"

#define N 1024
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

	DATA_TYPE (*A)[N][K] = (DATA_TYPE (*)[N][K]) malloc(N*K * sizeof(DATA_TYPE));
	DATA_TYPE (*C)[N][N] = (DATA_TYPE (*)[N][N]) malloc(N*N * sizeof(DATA_TYPE));

        init_data(N, K, *A);
	init_data(N, N, *C);


	start_time = get_time();
	pencil_dsyrk_n_u(N, K, 10.0, 1, *A, 10.0, 1, *C);
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
void pencil_dsyrk_n_u (int n_181_445,int k_182_446,double alpha_183_447,int ldA_184_448,double A_185_449[restrict const static n_181_445][k_182_446],double beta_186_450,int ldCT_187_451,double CT_188_452[restrict const static n_181_445][n_181_445]){

#pragma scop
{
#pragma pencil independent
for(int i_190_142_758_1529 = 0;i_190_142_758_1529 <= (n_181_445-1);i_190_142_758_1529 += 1){
#pragma pencil independent
for(int i_192_144_759_1531 = 0;i_192_144_759_1531 <= ((((n_181_445-1)-i_190_142_758_1529)+1)-1);i_192_144_759_1531 += 1){
double accum_193_145_760_1528;
accum_193_145_760_1528 = 0.0;
for(int i_194_146_761_1530 = 0;i_194_146_761_1530 <= (k_182_446-1);i_194_146_761_1530 += 1){
accum_193_145_760_1528 = (accum_193_145_760_1528+(A_185_449[i_190_142_758_1529][i_194_146_761_1530]*A_185_449[(i_192_144_759_1531+i_190_142_758_1529)][i_194_146_761_1530]));
}
CT_188_452[i_190_142_758_1529][(i_190_142_758_1529+i_192_144_759_1531)] = ((beta_186_450*CT_188_452[i_190_142_758_1529][(i_190_142_758_1529+i_192_144_759_1531)])+(alpha_183_447*accum_193_145_760_1528));
}
}
}
#pragma endscop

}
