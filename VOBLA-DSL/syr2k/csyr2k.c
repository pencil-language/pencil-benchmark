#include "../include/vobla_headers.h"

#define N 1024
#define K 1024
#define DATA_TYPE struct ComplexFloat
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
	pencil_csyr2k_nn_u(N, K, 10.0, 10.0, 1, *A, 1, *B, 10.0, 10.0, 1, *C);
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
void pencil_csyr2k_nn_u (int n_403_659,int k_404_660,float alpha_405__Re_1297,float alpha_405__Im_1298,int ldA_406_662,struct ComplexFloat A_407_663[restrict const static n_403_659][k_404_660],int ldB_408_664,struct ComplexFloat B_409_665[restrict const static k_404_660][n_403_659],float beta_410__Re_1299,float beta_410__Im_1300,int ldCT_411_667,struct ComplexFloat CT_412_668[restrict const static n_403_659][n_403_659]){

#pragma scop
{
#pragma pencil independent
for(int i_244_253_1202_2173 = 0;i_244_253_1202_2173 <= (n_403_659-1);i_244_253_1202_2173 += 1){
#pragma pencil independent
for(int i_246_255_1203_2174 = 0;i_246_255_1203_2174 <= ((((n_403_659-1)-i_244_253_1202_2173)+1)-1);i_246_255_1203_2174 += 1){
float accum_248__field1_756_1211_2177;
float var_255__field0_767_1208_2189;
float accum_248__field0_755_1206_2176;
float var_254__field0_765_1209_2187;
float var_254__field1_766_1210_2188;
float var_247__field1_754_1207_2175;
float var_247__field0_753_1204_2172;
float var_255__field1_768_1205_2190;
var_247__field0_753_1204_2172 = ((beta_410__Re_1299*CT_412_668[i_244_253_1202_2173][(i_244_253_1202_2173+i_246_255_1203_2174)].Re)-(beta_410__Im_1300*CT_412_668[i_244_253_1202_2173][(i_244_253_1202_2173+i_246_255_1203_2174)].Im));
var_247__field1_754_1207_2175 = ((beta_410__Im_1300*CT_412_668[i_244_253_1202_2173][(i_244_253_1202_2173+i_246_255_1203_2174)].Re)+(beta_410__Re_1299*CT_412_668[i_244_253_1202_2173][(i_244_253_1202_2173+i_246_255_1203_2174)].Im));
accum_248__field0_755_1206_2176 = 0.0f;
accum_248__field1_756_1211_2177 = 0.0f;
for(int i_249_260_1212_2179 = 0;i_249_260_1212_2179 <= (k_404_660-1);i_249_260_1212_2179 += 1){
float var_250__field1_758_1218_2180;
float var_253__field1_764_1217_2186;
float var_251__field1_760_1214_2182;
float var_250__field0_757_1219_2178;
float var_251__field0_759_1213_2181;
float var_253__field0_763_1216_2185;
float var_252__field0_761_1220_2183;
float var_252__field1_762_1215_2184;
var_250__field0_757_1219_2178 = ((A_407_663[i_244_253_1202_2173][i_249_260_1212_2179].Re*B_409_665[(i_246_255_1203_2174+i_244_253_1202_2173)][i_249_260_1212_2179].Re)-(A_407_663[i_244_253_1202_2173][i_249_260_1212_2179].Im*B_409_665[(i_246_255_1203_2174+i_244_253_1202_2173)][i_249_260_1212_2179].Im));
var_250__field1_758_1218_2180 = ((A_407_663[i_244_253_1202_2173][i_249_260_1212_2179].Im*B_409_665[(i_246_255_1203_2174+i_244_253_1202_2173)][i_249_260_1212_2179].Re)+(A_407_663[i_244_253_1202_2173][i_249_260_1212_2179].Re*B_409_665[(i_246_255_1203_2174+i_244_253_1202_2173)][i_249_260_1212_2179].Im));
var_251__field0_759_1213_2181 = ((B_409_665[i_244_253_1202_2173][i_249_260_1212_2179].Re*A_407_663[(i_246_255_1203_2174+i_244_253_1202_2173)][i_249_260_1212_2179].Re)-(B_409_665[i_244_253_1202_2173][i_249_260_1212_2179].Im*A_407_663[(i_246_255_1203_2174+i_244_253_1202_2173)][i_249_260_1212_2179].Im));
var_251__field1_760_1214_2182 = ((B_409_665[i_244_253_1202_2173][i_249_260_1212_2179].Im*A_407_663[(i_246_255_1203_2174+i_244_253_1202_2173)][i_249_260_1212_2179].Re)+(B_409_665[i_244_253_1202_2173][i_249_260_1212_2179].Re*A_407_663[(i_246_255_1203_2174+i_244_253_1202_2173)][i_249_260_1212_2179].Im));
var_252__field0_761_1220_2183 = (var_250__field0_757_1219_2178+var_251__field0_759_1213_2181);
var_252__field1_762_1215_2184 = (var_250__field1_758_1218_2180+var_251__field1_760_1214_2182);
var_253__field0_763_1216_2185 = (accum_248__field0_755_1206_2176+var_252__field0_761_1220_2183);
var_253__field1_764_1217_2186 = (accum_248__field1_756_1211_2177+var_252__field1_762_1215_2184);
accum_248__field0_755_1206_2176 = var_253__field0_763_1216_2185;
accum_248__field1_756_1211_2177 = var_253__field1_764_1217_2186;
}
var_254__field0_765_1209_2187 = ((alpha_405__Re_1297*accum_248__field0_755_1206_2176)-(alpha_405__Im_1298*accum_248__field1_756_1211_2177));
var_254__field1_766_1210_2188 = ((alpha_405__Im_1298*accum_248__field0_755_1206_2176)+(alpha_405__Re_1297*accum_248__field1_756_1211_2177));
var_255__field0_767_1208_2189 = (var_247__field0_753_1204_2172+var_254__field0_765_1209_2187);
var_255__field1_768_1205_2190 = (var_247__field1_754_1207_2175+var_254__field1_766_1210_2188);
CT_412_668[i_244_253_1202_2173][(i_244_253_1202_2173+i_246_255_1203_2174)].Re = var_255__field0_767_1208_2189;
CT_412_668[i_244_253_1202_2173][(i_244_253_1202_2173+i_246_255_1203_2174)].Im = var_255__field1_768_1205_2190;
}
}
}
#pragma endscop

}
