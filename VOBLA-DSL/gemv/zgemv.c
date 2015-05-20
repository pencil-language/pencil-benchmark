#include "../include/vobla_headers.h"

#define N 1024
#define M 1024
#define DATA_TYPE struct ComplexDouble
#define DATA_PRINTF_FORMAT "%f %f"

void print_array(int n, int m, DATA_TYPE B[n][m])
{
	int i, j;
	for (i=0; i<n; i++)
  	  for (j=0; j<m; j++)
		fprintf(DATA_OUTPUT_FILE, DATA_PRINTF_FORMAT, B[i][j].Re, B[i][j].Im);

	fprintf(DATA_OUTPUT_FILE, "\n");
}

void init_data(int n, int m, DATA_TYPE B[n][m])
{
	int i, j;
	for (i=0; i<n; i++)
   	  for (j=0; j<m; j++)
	  {
		B[i][j].Re = 2.0;
		B[i][j].Im = 2.0;
	  }
}

int main()
{
	double start_time, end_time;

	DATA_TYPE (*A)[M][N] = (DATA_TYPE (*)[M][N]) malloc(M*N * sizeof(DATA_TYPE));
	DATA_TYPE (*B)[N][1] = (DATA_TYPE (*)[N][1]) malloc(N*1 * sizeof(DATA_TYPE));
	DATA_TYPE (*C)[M][1] = (DATA_TYPE (*)[M][1]) malloc(M*1 * sizeof(DATA_TYPE));

        init_data(M, N, *A);
	init_data(N, 1, *B);
	init_data(M, 1, *C);


	start_time = get_time();
	pencil_zgemv_t_nn(M, N, 10.0, 10.0, 0, *A, 1, *B, 10.0, 10.0, 1, *C);
	end_time = get_time();

#ifdef BENCHMARK_TIME
	fprintf(TIME_OUTPUT_FILE, "%0.6lf", end_time - start_time);
	fprintf(TIME_OUTPUT_FILE, "\n");
#endif

#ifdef BENCHMARK_DUMP_ARRAYS
	print_array(M, 1, *C);
#endif
	return 0;
}

void pencil_zgemv_t_nn (int m_10877_10389,int n_10878_10390,double alpha_10879__Re_10593,double alpha_10879__Im_10594,int ldA_10880_10392,struct ComplexDouble A_10881_10393[restrict const static m_10877_10389][n_10878_10390],int incX_10882_10394,struct ComplexDouble X_10883_10395[restrict const static n_10878_10390][incX_10882_10394],double beta_10884__Re_10595,double beta_10884__Im_10596,int incY_10885_10397,struct ComplexDouble Y_10886_10398[restrict const static m_10877_10389][incY_10885_10397]){

#pragma scop
{
#pragma pencil independent
for(int i_2873_3382_7009_12624 = 0;i_2873_3382_7009_12624 <= (m_10877_10389-1);i_2873_3382_7009_12624 += 1){
double var_2874__field1_8894_7010_12625;
double var_2874__field0_8893_7011_12623;
var_2874__field0_8893_7011_12623 = ((Y_10886_10398[i_2873_3382_7009_12624][0].Re*beta_10884__Re_10595)-(Y_10886_10398[i_2873_3382_7009_12624][0].Im*beta_10884__Im_10596));
var_2874__field1_8894_7010_12625 = ((Y_10886_10398[i_2873_3382_7009_12624][0].Im*beta_10884__Re_10595)+(Y_10886_10398[i_2873_3382_7009_12624][0].Re*beta_10884__Im_10596));
Y_10886_10398[i_2873_3382_7009_12624][0].Re = var_2874__field0_8893_7011_12623;
Y_10886_10398[i_2873_3382_7009_12624][0].Im = var_2874__field1_8894_7010_12625;
}
for(int i_2876_3384_7012_12627 = 0;i_2876_3384_7012_12627 <= (n_10878_10390-1);i_2876_3384_7012_12627 += 1){
for(int i_2875_3385_7013_12628 = 0;i_2875_3385_7013_12628 <= (m_10877_10389-1);i_2875_3385_7013_12628 += 1){
double var_2878__field0_8897_7016_12630;
double var_2877__field0_8895_7019_12626;
double var_2878__field1_8898_7018_12631;
double var_2877__field1_8896_7017_12629;
double var_2879__field1_8900_7014_12633;
double var_2879__field0_8899_7015_12632;
var_2877__field0_8895_7019_12626 = ((alpha_10879__Re_10593*A_10881_10393[i_2876_3384_7012_12627][i_2875_3385_7013_12628].Re)-(alpha_10879__Im_10594*A_10881_10393[i_2876_3384_7012_12627][i_2875_3385_7013_12628].Im));
var_2877__field1_8896_7017_12629 = ((alpha_10879__Im_10594*A_10881_10393[i_2876_3384_7012_12627][i_2875_3385_7013_12628].Re)+(alpha_10879__Re_10593*A_10881_10393[i_2876_3384_7012_12627][i_2875_3385_7013_12628].Im));
var_2878__field0_8897_7016_12630 = ((var_2877__field0_8895_7019_12626*X_10883_10395[i_2876_3384_7012_12627][0].Re)-(var_2877__field1_8896_7017_12629*X_10883_10395[i_2876_3384_7012_12627][0].Im));
var_2878__field1_8898_7018_12631 = ((var_2877__field1_8896_7017_12629*X_10883_10395[i_2876_3384_7012_12627][0].Re)+(var_2877__field0_8895_7019_12626*X_10883_10395[i_2876_3384_7012_12627][0].Im));
var_2879__field0_8899_7015_12632 = (Y_10886_10398[i_2875_3385_7013_12628][0].Re+var_2878__field0_8897_7016_12630);
var_2879__field1_8900_7014_12633 = (Y_10886_10398[i_2875_3385_7013_12628][0].Im+var_2878__field1_8898_7018_12631);
Y_10886_10398[i_2875_3385_7013_12628][0].Re = var_2879__field0_8899_7015_12632;
Y_10886_10398[i_2875_3385_7013_12628][0].Im = var_2879__field1_8900_7014_12633;
}
}
}
#pragma endscop

}
