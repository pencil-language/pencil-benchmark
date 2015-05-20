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
	DATA_TYPE (*B)[N][1] = (DATA_TYPE (*)[N][1]) malloc(N * sizeof(DATA_TYPE));
	
        init_data(*A, N);
	init_data(*B, N);

	/* pencil_saxpy_nn */
	start_time = get_time();
	pencil_saxpy_nn(N, 10, 1, *A, 1, *B); 
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

void pencil_saxpy_nn (int n_53_313,float alpha_54_314,int incX_55_315,float X_56_316[restrict const static n_53_313][incX_55_315],int incY_57_317,float Y_58_318[restrict const static n_53_313][incY_57_317]){

#pragma scop
{
#pragma pencil independent
for(int i_138_158_744_1040 = 0;i_138_158_744_1040 <= (n_53_313-1);i_138_158_744_1040 += 1){
Y_58_318[i_138_158_744_1040][0] = (Y_58_318[i_138_158_744_1040][0]+(alpha_54_314*X_56_316[i_138_158_744_1040][0]));
}
}
#pragma endscop

}


