#include "../include/vobla_headers.h"

#define N 1024
#define M 1024
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

	DATA_TYPE (*A)[M][N] = (DATA_TYPE (*)[M][N]) malloc(M*N * sizeof(DATA_TYPE));
	DATA_TYPE (*B)[N][1] = (DATA_TYPE (*)[N][1]) malloc(N*1 * sizeof(DATA_TYPE));
	DATA_TYPE (*C)[M][1] = (DATA_TYPE (*)[M][1]) malloc(M*1 * sizeof(DATA_TYPE));

        init_data(M, N, *A);
	init_data(N, 1, *B);
	init_data(M, 1, *C);


	start_time = get_time();
	pencil_cgemv_t_nn(M, N, 10.0, 10.0, 0, *A, 1, *B, 10.0, 10.0, 1, *C);
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

void pencil_cgemv_t_nn (int m_10445_10233,int n_10446_10234,float alpha_10447__Re_10545,float alpha_10447__Im_10546,int ldA_10448_10236,struct ComplexFloat A_10449_10237[restrict const static m_10445_10233][n_10446_10234],int incX_10450_10238,struct ComplexFloat X_10451_10239[restrict const static n_10446_10234][incX_10450_10238],float beta_10452__Re_10547,float beta_10452__Im_10548,int incY_10453_10241,struct ComplexFloat Y_10454_10242[restrict const static m_10445_10233][incY_10453_10241]){

#pragma scop
{
#pragma pencil independent
for(int i_4819_4461_8672_12312 = 0;i_4819_4461_8672_12312 <= (m_10445_10233-1);i_4819_4461_8672_12312 += 1){
float var_4820__field0_9265_8673_12311;
float var_4820__field1_9266_8674_12313;
var_4820__field0_9265_8673_12311 = ((Y_10454_10242[i_4819_4461_8672_12312][0].Re*beta_10452__Re_10547)-(Y_10454_10242[i_4819_4461_8672_12312][0].Im*beta_10452__Im_10548));
var_4820__field1_9266_8674_12313 = ((Y_10454_10242[i_4819_4461_8672_12312][0].Im*beta_10452__Re_10547)+(Y_10454_10242[i_4819_4461_8672_12312][0].Re*beta_10452__Im_10548));
Y_10454_10242[i_4819_4461_8672_12312][0].Re = var_4820__field0_9265_8673_12311;
Y_10454_10242[i_4819_4461_8672_12312][0].Im = var_4820__field1_9266_8674_12313;
}
for(int i_4822_4463_8675_12315 = 0;i_4822_4463_8675_12315 <= (n_10446_10234-1);i_4822_4463_8675_12315 += 1){
for(int i_4821_4464_8676_12316 = 0;i_4821_4464_8676_12316 <= (m_10445_10233-1);i_4821_4464_8676_12316 += 1){
float var_4825__field0_9271_8681_12320;
float var_4824__field0_9269_8677_12318;
float var_4825__field1_9272_8682_12321;
float var_4824__field1_9270_8678_12319;
float var_4823__field1_9268_8680_12317;
float var_4823__field0_9267_8679_12314;
var_4823__field0_9267_8679_12314 = ((alpha_10447__Re_10545*A_10449_10237[i_4822_4463_8675_12315][i_4821_4464_8676_12316].Re)-(alpha_10447__Im_10546*A_10449_10237[i_4822_4463_8675_12315][i_4821_4464_8676_12316].Im));
var_4823__field1_9268_8680_12317 = ((alpha_10447__Im_10546*A_10449_10237[i_4822_4463_8675_12315][i_4821_4464_8676_12316].Re)+(alpha_10447__Re_10545*A_10449_10237[i_4822_4463_8675_12315][i_4821_4464_8676_12316].Im));
var_4824__field0_9269_8677_12318 = ((var_4823__field0_9267_8679_12314*X_10451_10239[i_4822_4463_8675_12315][0].Re)-(var_4823__field1_9268_8680_12317*X_10451_10239[i_4822_4463_8675_12315][0].Im));
var_4824__field1_9270_8678_12319 = ((var_4823__field1_9268_8680_12317*X_10451_10239[i_4822_4463_8675_12315][0].Re)+(var_4823__field0_9267_8679_12314*X_10451_10239[i_4822_4463_8675_12315][0].Im));
var_4825__field0_9271_8681_12320 = (Y_10454_10242[i_4821_4464_8676_12316][0].Re+var_4824__field0_9269_8677_12318);
var_4825__field1_9272_8682_12321 = (Y_10454_10242[i_4821_4464_8676_12316][0].Im+var_4824__field1_9270_8678_12319);
Y_10454_10242[i_4821_4464_8676_12316][0].Re = var_4825__field0_9271_8681_12320;
Y_10454_10242[i_4821_4464_8676_12316][0].Im = var_4825__field1_9272_8682_12321;
}
}
}
#pragma endscop

}
