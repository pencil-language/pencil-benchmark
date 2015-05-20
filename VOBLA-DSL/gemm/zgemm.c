#include "../include/vobla_headers.h"

#define N 1024
#define M 1024
#define K 1024
#define DATA_TYPE struct ComplexFloat
#define DATA_PRINTF_FORMAT "%f %f"

void print_array(int n, int m, DATA_TYPE B[n][m])
{
	int i, j;
	for (i=0; i<n; i++)
  	  for (j=0; j<m; j++)
		fprintf(DATA_OUTPUT_FILE, DATA_PRINTF_FORMAT, B[i][j].Re, B[i][j].Im);
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

	DATA_TYPE (*A)[M][K] = (DATA_TYPE (*)[M][K]) malloc(M*K * sizeof(DATA_TYPE));
	DATA_TYPE (*B)[K][N] = (DATA_TYPE (*)[K][N]) malloc(K*N * sizeof(DATA_TYPE));
	DATA_TYPE (*C)[M][N] = (DATA_TYPE (*)[M][N]) malloc(M*N * sizeof(DATA_TYPE));

        init_data(M, K, *A);
	init_data(K, N, *B);
	init_data(M, N, *C);


	/* pencil_saxpy_nn */
	start_time = get_time();
	pencil_cgemm_ttt(M, N, K, 10.0, 10.0, 0, *A, 0, *B, 10.0, 10.0, 0, *C);
	end_time = get_time();

#ifdef BENCHMARK_TIME
	fprintf(TIME_OUTPUT_FILE, "%0.6lf", end_time - start_time);
	fprintf(TIME_OUTPUT_FILE, "\n");
#endif

#ifdef BENCHMARK_DUMP_ARRAYS
	print_array(M, N, *C);
#endif
	return 0;
}

/*---- VOBLA function.  ----*/
void pencil_cgemm_ttt (int m_3648_3525,int n_3649_3526,int k_3650_3527,float alpha_3651__Re_3777,float alpha_3651__Im_3778,int ldA_3652_3529,struct ComplexFloat A_3653_3530[restrict const static m_3648_3525][k_3650_3527],int ldB_3654_3531,struct ComplexFloat B_3655_3532[restrict const static k_3650_3527][n_3649_3526],float beta_3656__Re_3779,float beta_3656__Im_3780,int ldC_3657_3534,struct ComplexFloat C_3658_3535[restrict const static m_3648_3525][n_3649_3526]){

#pragma scop
{
#pragma pencil independent
for(int i_1345_809_1880_5082 = 0;i_1345_809_1880_5082 <= (n_3649_3526-1);i_1345_809_1880_5082 += 1){
#pragma pencil independent
for(int i_1344_810_1881_5083 = 0;i_1344_810_1881_5083 <= (m_3648_3525-1);i_1344_810_1881_5083 += 1){
float var_1346__field1_2728_1883_5084;
float var_1346__field0_2727_1882_5081;
var_1346__field0_2727_1882_5081 = ((C_3658_3535[i_1345_809_1880_5082][i_1344_810_1881_5083].Re*beta_3656__Re_3779)-(C_3658_3535[i_1345_809_1880_5082][i_1344_810_1881_5083].Im*beta_3656__Im_3780));
var_1346__field1_2728_1883_5084 = ((C_3658_3535[i_1345_809_1880_5082][i_1344_810_1881_5083].Im*beta_3656__Re_3779)+(C_3658_3535[i_1345_809_1880_5082][i_1344_810_1881_5083].Re*beta_3656__Im_3780));
C_3658_3535[i_1345_809_1880_5082][i_1344_810_1881_5083].Re = var_1346__field0_2727_1882_5081;
C_3658_3535[i_1345_809_1880_5082][i_1344_810_1881_5083].Im = var_1346__field1_2728_1883_5084;
}
}
for(int i_1348_812_1884_5086 = 0;i_1348_812_1884_5086 <= (k_3650_3527-1);i_1348_812_1884_5086 += 1){
for(int i_1347_813_1885_5087 = 0;i_1347_813_1885_5087 <= (m_3648_3525-1);i_1347_813_1885_5087 += 1){
for(int i_1349_814_1886_5090 = 0;i_1349_814_1886_5090 <= (n_3649_3526-1);i_1349_814_1886_5090 += 1){
float var_1351__field0_2731_1887_5089;
float var_1351__field1_2732_1888_5091;
float var_1350__field1_2730_1889_5088;
float var_1352__field1_2734_1892_5093;
float var_1352__field0_2733_1891_5092;
float var_1350__field0_2729_1890_5085;
var_1350__field0_2729_1890_5085 = ((alpha_3651__Re_3777*A_3653_3530[i_1348_812_1884_5086][i_1347_813_1885_5087].Re)-(alpha_3651__Im_3778*A_3653_3530[i_1348_812_1884_5086][i_1347_813_1885_5087].Im));
var_1350__field1_2730_1889_5088 = ((alpha_3651__Im_3778*A_3653_3530[i_1348_812_1884_5086][i_1347_813_1885_5087].Re)+(alpha_3651__Re_3777*A_3653_3530[i_1348_812_1884_5086][i_1347_813_1885_5087].Im));
var_1351__field0_2731_1887_5089 = ((var_1350__field0_2729_1890_5085*B_3655_3532[i_1349_814_1886_5090][i_1348_812_1884_5086].Re)-(var_1350__field1_2730_1889_5088*B_3655_3532[i_1349_814_1886_5090][i_1348_812_1884_5086].Im));
var_1351__field1_2732_1888_5091 = ((var_1350__field1_2730_1889_5088*B_3655_3532[i_1349_814_1886_5090][i_1348_812_1884_5086].Re)+(var_1350__field0_2729_1890_5085*B_3655_3532[i_1349_814_1886_5090][i_1348_812_1884_5086].Im));
var_1352__field0_2733_1891_5092 = (C_3658_3535[i_1349_814_1886_5090][i_1347_813_1885_5087].Re+var_1351__field0_2731_1887_5089);
var_1352__field1_2734_1892_5093 = (C_3658_3535[i_1349_814_1886_5090][i_1347_813_1885_5087].Im+var_1351__field1_2732_1888_5091);
C_3658_3535[i_1349_814_1886_5090][i_1347_813_1885_5087].Re = var_1352__field0_2733_1891_5092;
C_3658_3535[i_1349_814_1886_5090][i_1347_813_1885_5087].Im = var_1352__field1_2734_1892_5093;
}
}
}
}
#pragma endscop

}
