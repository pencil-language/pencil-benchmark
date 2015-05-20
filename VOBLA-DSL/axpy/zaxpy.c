#include "../include/vobla_headers.h"

#define N 1024*1024*10
#define DATA_TYPE struct ComplexDouble
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
	pencil_zaxpy_nn(N, 10.0, 10.0, 1, *A, 1, *B); 
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

void pencil_zaxpy_nn (int n_313_409,double alpha_314__Re_713,double alpha_314__Im_714,int incX_315_411,struct ComplexDouble X_316_412[restrict const static n_313_409][incX_315_411],int incY_317_413,struct ComplexDouble Y_318_414[restrict const static n_313_409][incY_317_413]){
#pragma scop
{
#pragma pencil independent
for(int i_113_102_613_1170 = 0;i_113_102_613_1170 <= (n_313_409-1);i_113_102_613_1170 += 1){
double var_115__field0_379_615_1172;
double var_115__field1_380_617_1173;
double var_114__field0_377_616_1169;
double var_114__field1_378_614_1171;
var_114__field0_377_616_1169 = ((alpha_314__Re_713*X_316_412[i_113_102_613_1170][0].Re)-(alpha_314__Im_714*X_316_412[i_113_102_613_1170][0].Im));
var_114__field1_378_614_1171 = ((alpha_314__Im_714*X_316_412[i_113_102_613_1170][0].Re)+(alpha_314__Re_713*X_316_412[i_113_102_613_1170][0].Im));
var_115__field0_379_615_1172 = (Y_318_414[i_113_102_613_1170][0].Re+var_114__field0_377_616_1169);
var_115__field1_380_617_1173 = (Y_318_414[i_113_102_613_1170][0].Im+var_114__field1_378_614_1171);
Y_318_414[i_113_102_613_1170][0].Re = var_115__field0_379_615_1172;
Y_318_414[i_113_102_613_1170][0].Im = var_115__field1_380_617_1173;
}
}
#pragma endscop

}
