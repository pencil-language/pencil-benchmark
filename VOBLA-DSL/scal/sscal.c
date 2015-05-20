#include "../include/vobla_headers.h"

#define N 1024*1024*10
#define DATA_TYPE float
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
	pencil_sscal_n(N, 10.0, 1, *A); 
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
void pencil_sscal_n (int n_219_237,float alpha_220_238,int incX_221_239,float X_222_240[restrict const static n_219_237][incX_221_239])
{
#pragma scop
{
#pragma pencil independent
for(int i_59_73_163_377 = 0;i_59_73_163_377 <= (n_219_237-1);i_59_73_163_377 += 1){
X_222_240[i_59_73_163_377][0] = (alpha_220_238*X_222_240[i_59_73_163_377][0]);
}
}
#pragma endscop
}
