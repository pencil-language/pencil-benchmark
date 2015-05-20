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
	DATA_TYPE (*B)[K][N] = (DATA_TYPE (*)[K][N]) malloc(K*N * sizeof(DATA_TYPE));
	DATA_TYPE (*C)[N][N] = (DATA_TYPE (*)[N][N]) malloc(N*N * sizeof(DATA_TYPE));

        init_data(N, K, *A);
        init_data(K, N, *A);
	init_data(N, N, *C);


	start_time = get_time();
	pencil_dsyr2k_nn_u(N, K, 10.0, 1, *A, 1, *B, 10.0, 1, *C);
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
void pencil_dsyr2k_nn_u (int n_255_607,int k_256_608,double alpha_257_609,int ldA_258_610,double A_259_611[restrict const static n_255_607][k_256_608],int ldB_260_612,double B_261_613[restrict const static k_256_608][n_255_607],double beta_262_614,int ldCT_263_615,double CT_264_616[restrict const static n_255_607][n_255_607]){

#pragma scop
{
#pragma pencil independent
for(int i_323_273_1243_2069 = 0;i_323_273_1243_2069 <= (n_255_607-1);i_323_273_1243_2069 += 1){
#pragma pencil independent
for(int i_325_275_1244_2071 = 0;i_325_275_1244_2071 <= ((((n_255_607-1)-i_323_273_1243_2069)+1)-1);i_325_275_1244_2071 += 1){
double accum_326_276_1245_2068;
accum_326_276_1245_2068 = 0.0;
for(int i_327_277_1246_2070 = 0;i_327_277_1246_2070 <= (k_256_608-1);i_327_277_1246_2070 += 1){
accum_326_276_1245_2068 = (accum_326_276_1245_2068+((A_259_611[i_323_273_1243_2069][i_327_277_1246_2070]*B_261_613[(i_325_275_1244_2071+i_323_273_1243_2069)][i_327_277_1246_2070])+(B_261_613[i_323_273_1243_2069][i_327_277_1246_2070]*A_259_611[(i_325_275_1244_2071+i_323_273_1243_2069)][i_327_277_1246_2070])));
}
CT_264_616[i_323_273_1243_2069][(i_323_273_1243_2069+i_325_275_1244_2071)] = ((beta_262_614*CT_264_616[i_323_273_1243_2069][(i_323_273_1243_2069+i_325_275_1244_2071)])+(alpha_257_609*accum_326_276_1245_2068));
}
}
}
#pragma endscop

}
