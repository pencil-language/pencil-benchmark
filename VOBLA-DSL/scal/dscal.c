#include "../include/vobla_headers.h"

#define N 1024*1024*10
#define DATA_TYPE double
#define DATA_PRINTF_FORMAT "%f "

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
	
        init_data(*A, N);

	/* pencil_saxpy_nn */
	start_time = get_time();
	pencil_dscal_n(N, 10.0, 1, *A); 
	end_time = get_time();

#ifdef BENCHMARK_TIME
	fprintf(TIME_OUTPUT_FILE, "%0.6lf", end_time - start_time);
	fprintf(TIME_OUTPUT_FILE, "\n");
#endif

#ifdef BENCHMARK_DUMP_ARRAYS
	print_array(*A, N);
#endif
	return 0;
}

/*---- VOBLA function.  ----*/
void pencil_dscal_n (int n_243_247,double alpha_244_248,int incX_245_249,double X_246_250[restrict const static n_243_247][incX_245_249])
{
#pragma scop
{
#pragma pencil independent
for(int i_64_64_145_388 = 0;i_64_64_145_388 <= (n_243_247-1);i_64_64_145_388 += 1){
X_246_250[i_64_64_145_388][0] = (alpha_244_248*X_246_250[i_64_64_145_388][0]);
}
}
#pragma endscop
}
