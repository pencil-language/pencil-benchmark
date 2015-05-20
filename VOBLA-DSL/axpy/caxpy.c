#include "../include/vobla_headers.h"

#define N 1024*1024*10
#define DATA_TYPE struct ComplexFloat
#define DATA_PRINTF_FORMAT "%f %f"

void print_array(DATA_TYPE B[N][1], int n)
{
	int i;
	for (i=0; i<n; i++)
		fprintf(DATA_OUTPUT_FILE, DATA_PRINTF_FORMAT, B[i][0].Re, B[i][0].Im);
}

void init_data(DATA_TYPE B[N][1], int n)
{
	int i;
	for (i=0; i<n; i++)
	{
		B[i][0].Re = 2.0;
		B[i][0].Im = 2.0;
	}
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
	pencil_caxpy_nn(N, 10.0, 10.0, 1, *A, 1, *B); 
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

void pencil_caxpy_nn (int n_225_377,float alpha_226__Re_705,float alpha_226__Im_706,int incX_227_379,struct ComplexFloat X_228_380[restrict const static n_225_377][incX_227_379],int incY_229_381,struct ComplexFloat Y_230_382[restrict const static n_225_377][incY_229_381])
{
#pragma scop
{
#pragma pencil independent
for(int i_132_86_575_1114 = 0;i_132_86_575_1114 <= (n_225_377-1);i_132_86_575_1114 += 1){
float var_133__field1_370_576_1115;
float var_134__field0_371_578_1116;
float var_134__field1_372_579_1117;
float var_133__field0_369_577_1113;
var_133__field0_369_577_1113 = ((alpha_226__Re_705*X_228_380[i_132_86_575_1114][0].Re)-(alpha_226__Im_706*X_228_380[i_132_86_575_1114][0].Im));
var_133__field1_370_576_1115 = ((alpha_226__Im_706*X_228_380[i_132_86_575_1114][0].Re)+(alpha_226__Re_705*X_228_380[i_132_86_575_1114][0].Im));
var_134__field0_371_578_1116 = (Y_230_382[i_132_86_575_1114][0].Re+var_133__field0_369_577_1113);
var_134__field1_372_579_1117 = (Y_230_382[i_132_86_575_1114][0].Im+var_133__field1_370_576_1115);
Y_230_382[i_132_86_575_1114][0].Re = var_134__field0_371_578_1116;
Y_230_382[i_132_86_575_1114][0].Im = var_134__field1_372_579_1117;
}
}
#pragma endscop
}
