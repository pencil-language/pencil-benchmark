#include "../include/vobla_headers.h"

#define N 1024
#define K 1024
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

	DATA_TYPE (*A)[N][K] = (DATA_TYPE (*)[N][K]) malloc(N*K * sizeof(DATA_TYPE));
	DATA_TYPE (*B)[K][N] = (DATA_TYPE (*)[K][N]) malloc(K*N * sizeof(DATA_TYPE));
	DATA_TYPE (*C)[N][N] = (DATA_TYPE (*)[N][N]) malloc(N*N * sizeof(DATA_TYPE));

        init_data(N, K, *A);
        init_data(K, N, *A);
	init_data(N, N, *C);


	start_time = get_time();
	pencil_zsyr2k_nn_u(N, K, 10.0, 10.0, 1, *A, 1, *B, 10.0, 10.0, 1, *C);
	end_time = get_time();

#ifdef BENCHMARK_TIME
	fprintf(TIME_OUTPUT_FILE, "%0.6lf", end_time - start_time);
	fprintf(TIME_OUTPUT_FILE, "\n");
#endif

#ifdef BENCHMARK_DUMP_ARRAYS
	print_array(N, N, *C);
#endif
	return 0;
}

/*---- VOBLA function.  ----*/
void pencil_zsyr2k_nn_u (int n_551_711,int k_552_712,double alpha_553__Re_1313,double alpha_553__Im_1314,int ldA_554_714,struct ComplexDouble A_555_715[restrict const static n_551_711][k_552_712],int ldB_556_716,struct ComplexDouble B_557_717[restrict const static k_552_712][n_551_711],double beta_558__Re_1315,double beta_558__Im_1316,int ldCT_559_719,struct ComplexDouble CT_560_720[restrict const static n_551_711][n_551_711]){

#pragma scop
{
#pragma pencil independent
for(int i_356_185_1064_2309 = 0;i_356_185_1064_2309 <= (n_551_711-1);i_356_185_1064_2309 += 1){
#pragma pencil independent
for(int i_358_187_1065_2310 = 0;i_358_187_1065_2310 <= ((((n_551_711-1)-i_356_185_1064_2309)+1)-1);i_358_187_1065_2310 += 1){
double var_367__field0_735_1069_2325;
double accum_360__field1_724_1072_2313;
double var_359__field1_722_1067_2311;
double var_367__field1_736_1070_2326;
double accum_360__field0_723_1071_2312;
double var_359__field0_721_1073_2308;
double var_366__field0_733_1066_2323;
double var_366__field1_734_1068_2324;
var_359__field0_721_1073_2308 = ((beta_558__Re_1315*CT_560_720[i_356_185_1064_2309][(i_356_185_1064_2309+i_358_187_1065_2310)].Re)-(beta_558__Im_1316*CT_560_720[i_356_185_1064_2309][(i_356_185_1064_2309+i_358_187_1065_2310)].Im));
var_359__field1_722_1067_2311 = ((beta_558__Im_1316*CT_560_720[i_356_185_1064_2309][(i_356_185_1064_2309+i_358_187_1065_2310)].Re)+(beta_558__Re_1315*CT_560_720[i_356_185_1064_2309][(i_356_185_1064_2309+i_358_187_1065_2310)].Im));
accum_360__field0_723_1071_2312 = 0.0;
accum_360__field1_724_1072_2313 = 0.0;
for(int i_361_192_1074_2315 = 0;i_361_192_1074_2315 <= (k_552_712-1);i_361_192_1074_2315 += 1){
double var_363__field0_727_1081_2317;
double var_364__field0_729_1076_2319;
double var_364__field1_730_1077_2320;
double var_365__field1_732_1079_2322;
double var_362__field1_726_1075_2316;
double var_363__field1_728_1082_2318;
double var_365__field0_731_1080_2321;
double var_362__field0_725_1078_2314;
var_362__field0_725_1078_2314 = ((A_555_715[i_356_185_1064_2309][i_361_192_1074_2315].Re*B_557_717[(i_358_187_1065_2310+i_356_185_1064_2309)][i_361_192_1074_2315].Re)-(A_555_715[i_356_185_1064_2309][i_361_192_1074_2315].Im*B_557_717[(i_358_187_1065_2310+i_356_185_1064_2309)][i_361_192_1074_2315].Im));
var_362__field1_726_1075_2316 = ((A_555_715[i_356_185_1064_2309][i_361_192_1074_2315].Im*B_557_717[(i_358_187_1065_2310+i_356_185_1064_2309)][i_361_192_1074_2315].Re)+(A_555_715[i_356_185_1064_2309][i_361_192_1074_2315].Re*B_557_717[(i_358_187_1065_2310+i_356_185_1064_2309)][i_361_192_1074_2315].Im));
var_363__field0_727_1081_2317 = ((B_557_717[i_356_185_1064_2309][i_361_192_1074_2315].Re*A_555_715[(i_358_187_1065_2310+i_356_185_1064_2309)][i_361_192_1074_2315].Re)-(B_557_717[i_356_185_1064_2309][i_361_192_1074_2315].Im*A_555_715[(i_358_187_1065_2310+i_356_185_1064_2309)][i_361_192_1074_2315].Im));
var_363__field1_728_1082_2318 = ((B_557_717[i_356_185_1064_2309][i_361_192_1074_2315].Im*A_555_715[(i_358_187_1065_2310+i_356_185_1064_2309)][i_361_192_1074_2315].Re)+(B_557_717[i_356_185_1064_2309][i_361_192_1074_2315].Re*A_555_715[(i_358_187_1065_2310+i_356_185_1064_2309)][i_361_192_1074_2315].Im));
var_364__field0_729_1076_2319 = (var_362__field0_725_1078_2314+var_363__field0_727_1081_2317);
var_364__field1_730_1077_2320 = (var_362__field1_726_1075_2316+var_363__field1_728_1082_2318);
var_365__field0_731_1080_2321 = (accum_360__field0_723_1071_2312+var_364__field0_729_1076_2319);
var_365__field1_732_1079_2322 = (accum_360__field1_724_1072_2313+var_364__field1_730_1077_2320);
accum_360__field0_723_1071_2312 = var_365__field0_731_1080_2321;
accum_360__field1_724_1072_2313 = var_365__field1_732_1079_2322;
}
var_366__field0_733_1066_2323 = ((alpha_553__Re_1313*accum_360__field0_723_1071_2312)-(alpha_553__Im_1314*accum_360__field1_724_1072_2313));
var_366__field1_734_1068_2324 = ((alpha_553__Im_1314*accum_360__field0_723_1071_2312)+(alpha_553__Re_1313*accum_360__field1_724_1072_2313));
var_367__field0_735_1069_2325 = (var_359__field0_721_1073_2308+var_366__field0_733_1066_2323);
var_367__field1_736_1070_2326 = (var_359__field1_722_1067_2311+var_366__field1_734_1068_2324);
CT_560_720[i_356_185_1064_2309][(i_356_185_1064_2309+i_358_187_1065_2310)].Re = var_367__field0_735_1069_2325;
CT_560_720[i_356_185_1064_2309][(i_356_185_1064_2309+i_358_187_1065_2310)].Im = var_367__field1_736_1070_2326;
}
}
}
#pragma endscop

}
