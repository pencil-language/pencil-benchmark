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
	DATA_TYPE (*C)[N][N] = (DATA_TYPE (*)[N][N]) malloc(N*N * sizeof(DATA_TYPE));

        init_data(N, K, *A);
	init_data(N, N, *C);


	start_time = get_time();
	pencil_zsyrk_n_u(N, K, 10.0, 10.0, 1, *A, 10.0, 10.0, 1, *C);
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
void pencil_zsyrk_n_u (int n_397_525,int k_398_526,double alpha_399__Re_977,double alpha_399__Im_978,int ldA_400_528,struct ComplexDouble A_401_529[restrict const static n_397_525][k_398_526],double beta_402__Re_979,double beta_402__Im_980,int ldCT_403_531,struct ComplexDouble CT_404_532[restrict const static n_397_525][n_397_525]){

#pragma scop
{
#pragma pencil independent
for(int i_248_199_862_1713 = 0;i_248_199_862_1713 <= (n_397_525-1);i_248_199_862_1713 += 1){
#pragma pencil independent
for(int i_250_201_863_1714 = 0;i_250_201_863_1714 <= ((((n_397_525-1)-i_248_199_862_1713)+1)-1);i_250_201_863_1714 += 1){
double var_256__field1_574_869_1724;
double var_256__field0_573_871_1723;
double var_257__field0_575_870_1725;
double var_251__field0_565_868_1712;
double var_257__field1_576_865_1726;
double var_251__field1_566_866_1715;
double accum_252__field0_567_867_1716;
double accum_252__field1_568_864_1717;
var_251__field0_565_868_1712 = ((beta_402__Re_979*CT_404_532[i_248_199_862_1713][(i_248_199_862_1713+i_250_201_863_1714)].Re)-(beta_402__Im_980*CT_404_532[i_248_199_862_1713][(i_248_199_862_1713+i_250_201_863_1714)].Im));
var_251__field1_566_866_1715 = ((beta_402__Im_980*CT_404_532[i_248_199_862_1713][(i_248_199_862_1713+i_250_201_863_1714)].Re)+(beta_402__Re_979*CT_404_532[i_248_199_862_1713][(i_248_199_862_1713+i_250_201_863_1714)].Im));
accum_252__field0_567_867_1716 = 0.0;
accum_252__field1_568_864_1717 = 0.0;
for(int i_253_206_872_1719 = 0;i_253_206_872_1719 <= (k_398_526-1);i_253_206_872_1719 += 1){
double var_254__field1_570_874_1720;
double var_254__field0_569_876_1718;
double var_255__field0_571_873_1721;
double var_255__field1_572_875_1722;
var_254__field0_569_876_1718 = ((A_401_529[i_248_199_862_1713][i_253_206_872_1719].Re*A_401_529[(i_250_201_863_1714+i_248_199_862_1713)][i_253_206_872_1719].Re)-(A_401_529[i_248_199_862_1713][i_253_206_872_1719].Im*A_401_529[(i_250_201_863_1714+i_248_199_862_1713)][i_253_206_872_1719].Im));
var_254__field1_570_874_1720 = ((A_401_529[i_248_199_862_1713][i_253_206_872_1719].Im*A_401_529[(i_250_201_863_1714+i_248_199_862_1713)][i_253_206_872_1719].Re)+(A_401_529[i_248_199_862_1713][i_253_206_872_1719].Re*A_401_529[(i_250_201_863_1714+i_248_199_862_1713)][i_253_206_872_1719].Im));
var_255__field0_571_873_1721 = (accum_252__field0_567_867_1716+var_254__field0_569_876_1718);
var_255__field1_572_875_1722 = (accum_252__field1_568_864_1717+var_254__field1_570_874_1720);
accum_252__field0_567_867_1716 = var_255__field0_571_873_1721;
accum_252__field1_568_864_1717 = var_255__field1_572_875_1722;
}
var_256__field0_573_871_1723 = ((alpha_399__Re_977*accum_252__field0_567_867_1716)-(alpha_399__Im_978*accum_252__field1_568_864_1717));
var_256__field1_574_869_1724 = ((alpha_399__Im_978*accum_252__field0_567_867_1716)+(alpha_399__Re_977*accum_252__field1_568_864_1717));
var_257__field0_575_870_1725 = (var_251__field0_565_868_1712+var_256__field0_573_871_1723);
var_257__field1_576_865_1726 = (var_251__field1_566_866_1715+var_256__field1_574_869_1724);
CT_404_532[i_248_199_862_1713][(i_248_199_862_1713+i_250_201_863_1714)].Re = var_257__field0_575_870_1725;
CT_404_532[i_248_199_862_1713][(i_248_199_862_1713+i_250_201_863_1714)].Im = var_257__field1_576_865_1726;
}
}
}
#pragma endscop

}
