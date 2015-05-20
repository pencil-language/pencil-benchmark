#include "../include/vobla_headers.h"

#define N 1024*1024*10
#define DATA_TYPE double
#define DATA_PRINTF_FORMAT "%lf "

void print_array(DATA_TYPE B[N][1], int n)
{
	int i;
	for (i=0; i<n; i++)
		fprintf(DATA_OUTPUT_FILE, DATA_PRINTF_FORMAT, B[i][0]);
}

void init_data(DATA_TYPE B[N][1], int n)
{
	int i;
	for (i=0; i<n; i++)
		B[i][0] = 2.0;
}

int main()
{
	double start_time, end_time;

	DATA_TYPE (*A)[N][1] = (DATA_TYPE (*)[N][1]) malloc(N * sizeof(DATA_TYPE));
	DATA_TYPE (*B)[N][1] = (DATA_TYPE (*)[N][1]) malloc(N * sizeof(DATA_TYPE));
	
        init_data(*A, N);
	init_data(*B, N);

	/* pencil_saxpy_nn */
	start_time = get_time();
	pencil_daxpy_nn(N, 10, 1, *A, 1, *B); 
	end_time = get_time();

#ifdef BENCHMARK_TIME
	fprintf(TIME_OUTPUT_FILE, "%0.6lf", end_time - start_time);
	fprintf(TIME_OUTPUT_FILE, "\n");
#endif

#ifdef BENCHMARK_DUMP_ARRAYS
	print_array(*B, N);
#endif
	return 0;
}

/*---- VOBLA function.  ----*/
void pencil_daxpy_nn (int n_137_345,double alpha_138_346,int incX_139_347,double X_140_348[restrict const static n_137_345][incX_139_347],int incY_141_349,double Y_142_350[restrict const static n_137_345][incY_141_349]){

#pragma scop
{
#pragma pencil independent
for(int i_123_110_631_1076 = 0;i_123_110_631_1076 <= (n_137_345-1);i_123_110_631_1076 += 1){
Y_142_350[i_123_110_631_1076][0] = (Y_142_350[i_123_110_631_1076][0]+(alpha_138_346*X_140_348[i_123_110_631_1076][0]));
}
}
#pragma endscop

}
