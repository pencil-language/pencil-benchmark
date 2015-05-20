#include "../include/vobla_headers.h"

#define N 1024
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

	DATA_TYPE (*A)[N][N] = (DATA_TYPE (*)[N][N]) malloc(N*N * sizeof(DATA_TYPE));
	DATA_TYPE (*C)[N][1] = (DATA_TYPE (*)[N][1]) malloc(N*1 * sizeof(DATA_TYPE));

        init_data(N, N, *A);
	init_data(N, 1, *C);


	start_time = get_time();
	pencil_dtrmv_n_nun(N, 10.0, *A, 1, *C);
	end_time = get_time();

#ifdef BENCHMARK_TIME
	fprintf(TIME_OUTPUT_FILE, "%0.6lf", end_time - start_time);
	fprintf(TIME_OUTPUT_FILE, "\n");
#endif

#ifdef BENCHMARK_DUMP_ARRAYS
	print_array(N, 1, *C);
#endif
	return 0;
}

/*---- VOBLA function.  ----*/
void pencil_dtrmv_n_nun (int n_8292_8308,int ldAT_8293_8309,double AT_8294_8310[restrict const static n_8292_8308][n_8292_8308],int incX_8295_8311,double X_8296_8312[restrict const static n_8292_8308][incX_8295_8311]){

#pragma scop
{
for(int i_2447_1705_4072_10751 = 0;i_2447_1705_4072_10751 <= (n_8292_8308-1);i_2447_1705_4072_10751 += 1){
double accum_2449_1707_4073_10750[n_8292_8308];
accum_2449_1707_4073_10750[i_2447_1705_4072_10751] = 0.0;
for(int i_2450_1708_4074_10752 = 0;i_2450_1708_4074_10752 <= ((((n_8292_8308-1)-i_2447_1705_4072_10751)+1)-1);i_2450_1708_4074_10752 += 1){
accum_2449_1707_4073_10750[i_2447_1705_4072_10751] = (accum_2449_1707_4073_10750[i_2447_1705_4072_10751]+(AT_8294_8310[i_2447_1705_4072_10751][(i_2447_1705_4072_10751+i_2450_1708_4074_10752)]*X_8296_8312[(i_2447_1705_4072_10751+i_2450_1708_4074_10752)][0]));
}
X_8296_8312[i_2447_1705_4072_10751][0] = accum_2449_1707_4073_10750[i_2447_1705_4072_10751];
}
}
#pragma endscop

}
