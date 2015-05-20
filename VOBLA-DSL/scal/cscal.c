#include "../include/vobla_headers.h"

#define N 1024*1024*10
#define DATA_TYPE struct ComplexFloat
#define DATA_PRINTF_FORMAT "%f %f "

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
	pencil_cscal_n(N, 10.0, 10.0, 1, *A); 
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
void pencil_cscal_n (int n_269_257,float alpha_270__Re_277,float alpha_270__Im_278,int incX_271_259,struct ComplexFloat X_272_260[restrict const static n_269_257][incX_271_259]){

#pragma scop
{
#pragma pencil independent
for(int i_54_86_191_401 = 0;i_54_86_191_401 <= (n_269_257-1);i_54_86_191_401 += 1){
float var_55__field1_174_192_402;
float var_55__field0_173_193_400;
var_55__field0_173_193_400 = ((alpha_270__Re_277*X_272_260[i_54_86_191_401][0].Re)-(alpha_270__Im_278*X_272_260[i_54_86_191_401][0].Im));
var_55__field1_174_192_402 = ((alpha_270__Im_278*X_272_260[i_54_86_191_401][0].Re)+(alpha_270__Re_277*X_272_260[i_54_86_191_401][0].Im));
X_272_260[i_54_86_191_401][0].Re = var_55__field0_173_193_400;
X_272_260[i_54_86_191_401][0].Im = var_55__field1_174_192_402;
}
}
#pragma endscop

}
