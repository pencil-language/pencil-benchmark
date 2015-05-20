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
	DATA_TYPE (*C)[N][N] = (DATA_TYPE (*)[N][N]) malloc(N*N * sizeof(DATA_TYPE));

        init_data(N, K, *A);
	init_data(N, N, *C);


	start_time = get_time();
	pencil_csyrk_n_u(N, K, 10.0, 10.0, 1, *A, 10.0, 10.0, 1, *C);
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
void pencil_csyrk_n_u (int n_289_485,int k_290_486,float alpha_291__Re_961,float alpha_291__Im_962,int ldA_292_488,struct ComplexFloat A_293_489[restrict const static n_289_485][k_290_486],float beta_294__Re_963,float beta_294__Im_964,int ldCT_295_491,struct ComplexFloat CT_296_492[restrict const static n_289_485][n_289_485]){

#pragma scop
{
#pragma pencil independent
for(int i_261_288_1028_1609 = 0;i_261_288_1028_1609 <= (n_289_485-1);i_261_288_1028_1609 += 1){
#pragma pencil independent
for(int i_263_290_1029_1610 = 0;i_263_290_1029_1610 <= ((((n_289_485-1)-i_261_288_1028_1609)+1)-1);i_263_290_1029_1610 += 1){
float var_264__field0_613_1030_1608;
float var_269__field1_622_1037_1620;
float var_270__field0_623_1034_1621;
float accum_265__field0_615_1031_1612;
float var_270__field1_624_1032_1622;
float var_269__field0_621_1035_1619;
float accum_265__field1_616_1033_1613;
float var_264__field1_614_1036_1611;
var_264__field0_613_1030_1608 = ((beta_294__Re_963*CT_296_492[i_261_288_1028_1609][(i_261_288_1028_1609+i_263_290_1029_1610)].Re)-(beta_294__Im_964*CT_296_492[i_261_288_1028_1609][(i_261_288_1028_1609+i_263_290_1029_1610)].Im));
var_264__field1_614_1036_1611 = ((beta_294__Im_964*CT_296_492[i_261_288_1028_1609][(i_261_288_1028_1609+i_263_290_1029_1610)].Re)+(beta_294__Re_963*CT_296_492[i_261_288_1028_1609][(i_261_288_1028_1609+i_263_290_1029_1610)].Im));
accum_265__field0_615_1031_1612 = 0.0f;
accum_265__field1_616_1033_1613 = 0.0f;
for(int i_266_295_1038_1615 = 0;i_266_295_1038_1615 <= (k_290_486-1);i_266_295_1038_1615 += 1){
float var_268__field1_620_1041_1618;
float var_267__field0_617_1039_1614;
float var_268__field0_619_1042_1617;
float var_267__field1_618_1040_1616;
var_267__field0_617_1039_1614 = ((A_293_489[i_261_288_1028_1609][i_266_295_1038_1615].Re*A_293_489[(i_263_290_1029_1610+i_261_288_1028_1609)][i_266_295_1038_1615].Re)-(A_293_489[i_261_288_1028_1609][i_266_295_1038_1615].Im*A_293_489[(i_263_290_1029_1610+i_261_288_1028_1609)][i_266_295_1038_1615].Im));
var_267__field1_618_1040_1616 = ((A_293_489[i_261_288_1028_1609][i_266_295_1038_1615].Im*A_293_489[(i_263_290_1029_1610+i_261_288_1028_1609)][i_266_295_1038_1615].Re)+(A_293_489[i_261_288_1028_1609][i_266_295_1038_1615].Re*A_293_489[(i_263_290_1029_1610+i_261_288_1028_1609)][i_266_295_1038_1615].Im));
var_268__field0_619_1042_1617 = (accum_265__field0_615_1031_1612+var_267__field0_617_1039_1614);
var_268__field1_620_1041_1618 = (accum_265__field1_616_1033_1613+var_267__field1_618_1040_1616);
accum_265__field0_615_1031_1612 = var_268__field0_619_1042_1617;
accum_265__field1_616_1033_1613 = var_268__field1_620_1041_1618;
}
var_269__field0_621_1035_1619 = ((alpha_291__Re_961*accum_265__field0_615_1031_1612)-(alpha_291__Im_962*accum_265__field1_616_1033_1613));
var_269__field1_622_1037_1620 = ((alpha_291__Im_962*accum_265__field0_615_1031_1612)+(alpha_291__Re_961*accum_265__field1_616_1033_1613));
var_270__field0_623_1034_1621 = (var_264__field0_613_1030_1608+var_269__field0_621_1035_1619);
var_270__field1_624_1032_1622 = (var_264__field1_614_1036_1611+var_269__field1_622_1037_1620);
CT_296_492[i_261_288_1028_1609][(i_261_288_1028_1609+i_263_290_1029_1610)].Re = var_270__field0_623_1034_1621;
CT_296_492[i_261_288_1028_1609][(i_261_288_1028_1609+i_263_290_1029_1610)].Im = var_270__field1_624_1032_1622;
}
}
}
#pragma endscop

}
