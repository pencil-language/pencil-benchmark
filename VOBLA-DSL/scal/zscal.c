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
	
        init_data(*A, N);

	/* pencil_saxpy_nn */
	start_time = get_time();
	pencil_zscal_n(N, 10.0, 10.0, 1, *A); 
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
void pencil_zscal_n (int n_295_267,double alpha_296__Re_281,double alpha_296__Im_282,int incX_297_269,struct ComplexDouble X_298_270[restrict const static n_295_267][incX_297_269]){

#pragma scop
{
#pragma pencil independent
for(int i_56_49_115_418 = 0;i_56_49_115_418 <= (n_295_267-1);i_56_49_115_418 += 1){
double var_57__field1_164_117_419;
double var_57__field0_163_116_417;
var_57__field0_163_116_417 = ((alpha_296__Re_281*X_298_270[i_56_49_115_418][0].Re)-(alpha_296__Im_282*X_298_270[i_56_49_115_418][0].Im));
var_57__field1_164_117_419 = ((alpha_296__Im_282*X_298_270[i_56_49_115_418][0].Re)+(alpha_296__Re_281*X_298_270[i_56_49_115_418][0].Im));
X_298_270[i_56_49_115_418][0].Re = var_57__field0_163_116_417;
X_298_270[i_56_49_115_418][0].Im = var_57__field1_164_117_419;
}
}
#pragma endscop

}
