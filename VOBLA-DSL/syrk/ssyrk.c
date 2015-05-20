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
	DATA_TYPE (*C)[N][N] = (DATA_TYPE (*)[N][N]) malloc(N*N * sizeof(DATA_TYPE));

        init_data(N, K, *A);
	init_data(N, N, *C);


	start_time = get_time();
	pencil_ssyrk_n_u(N, K, 10.0, 1, *A, 10.0, 1, *C);
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
void pencil_ssyrk_n_u (int n_81_405,int k_82_406,float alpha_83_407,int ldA_84_408,float A_85_409[restrict const static n_81_405][k_82_406],float beta_86_410,int ldCT_87_411,float CT_88_412[restrict const static n_81_405][n_81_405]){

#pragma scop
{
#pragma pencil independent
for(int i_323_231_924_1477 = 0;i_323_231_924_1477 <= (n_81_405-1);i_323_231_924_1477 += 1){
#pragma pencil independent
for(int i_325_233_925_1479 = 0;i_325_233_925_1479 <= ((((n_81_405-1)-i_323_231_924_1477)+1)-1);i_325_233_925_1479 += 1){
float accum_326_234_926_1476;
accum_326_234_926_1476 = 0.0f;
for(int i_327_235_927_1478 = 0;i_327_235_927_1478 <= (k_82_406-1);i_327_235_927_1478 += 1){
accum_326_234_926_1476 = (accum_326_234_926_1476+(A_85_409[i_323_231_924_1477][i_327_235_927_1478]*A_85_409[(i_325_233_925_1479+i_323_231_924_1477)][i_327_235_927_1478]));
}
CT_88_412[i_323_231_924_1477][(i_323_231_924_1477+i_325_233_925_1479)] = ((beta_86_410*CT_88_412[i_323_231_924_1477][(i_323_231_924_1477+i_325_233_925_1479)])+(alpha_83_407*accum_326_234_926_1476));
}
}
}
#pragma endscop

}
