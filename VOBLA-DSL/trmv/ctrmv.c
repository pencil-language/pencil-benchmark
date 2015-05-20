#include "../include/vobla_headers.h"

#define N 1024
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

	DATA_TYPE (*A)[N][N] = (DATA_TYPE (*)[N][N]) malloc(N*N * sizeof(DATA_TYPE));
	DATA_TYPE (*C)[N][1] = (DATA_TYPE (*)[N][1]) malloc(N*1 * sizeof(DATA_TYPE));

        init_data(N, N, *A);
	init_data(N, 1, *C);


	start_time = get_time();
	pencil_ctrmv_n_nun(N, 10.0, *A, 1, *C);
	end_time = get_time();

#ifdef BENCHMARK_TIME
	fprintf(TIME_OUTPUT_FILE, "%0.6lf", end_time - start_time);
	fprintf(TIME_OUTPUT_FILE, "\n");
#endif

#ifdef BENCHMARK_DUMP_ARRAYS
	print_array(N, 1, *C);
#endif
	return 0;
}

/*---- VOBLA function.  ----*/
void pencil_ctrmv_n_nun (int n_8692_8420,int ldAT_8693_8421,struct ComplexFloat AT_8694_8422[restrict const static n_8692_8420][n_8692_8420],int incX_8695_8423,struct ComplexFloat X_8696_8424[restrict const static n_8692_8420][incX_8695_8423]){

#pragma scop
{
for(int i_2151_1302_3429_10982 = 0;i_2151_1302_3429_10982 <= (n_8692_8420-1);i_2151_1302_3429_10982 += 1){
float accum_2153__field1_6578_3431_10980[n_8692_8420];
float accum_2153__field0_6577_3430_10979[n_8692_8420];
accum_2153__field0_6577_3430_10979[i_2151_1302_3429_10982] = 0.0f;
accum_2153__field1_6578_3431_10980[i_2151_1302_3429_10982] = 0.0f;
for(int i_2154_1305_3432_10983 = 0;i_2154_1305_3432_10983 <= ((((n_8692_8420-1)-i_2151_1302_3429_10982)+1)-1);i_2154_1305_3432_10983 += 1){
float var_2156__field0_6581_3435_10985;
float var_2155__field0_6579_3434_10981;
float var_2156__field1_6582_3436_10986;
float var_2155__field1_6580_3433_10984;
var_2155__field0_6579_3434_10981 = ((AT_8694_8422[i_2151_1302_3429_10982][(i_2151_1302_3429_10982+i_2154_1305_3432_10983)].Re*X_8696_8424[(i_2151_1302_3429_10982+i_2154_1305_3432_10983)][0].Re)-(AT_8694_8422[i_2151_1302_3429_10982][(i_2151_1302_3429_10982+i_2154_1305_3432_10983)].Im*X_8696_8424[(i_2151_1302_3429_10982+i_2154_1305_3432_10983)][0].Im));
var_2155__field1_6580_3433_10984 = ((AT_8694_8422[i_2151_1302_3429_10982][(i_2151_1302_3429_10982+i_2154_1305_3432_10983)].Im*X_8696_8424[(i_2151_1302_3429_10982+i_2154_1305_3432_10983)][0].Re)+(AT_8694_8422[i_2151_1302_3429_10982][(i_2151_1302_3429_10982+i_2154_1305_3432_10983)].Re*X_8696_8424[(i_2151_1302_3429_10982+i_2154_1305_3432_10983)][0].Im));
var_2156__field0_6581_3435_10985 = (accum_2153__field0_6577_3430_10979[i_2151_1302_3429_10982]+var_2155__field0_6579_3434_10981);
var_2156__field1_6582_3436_10986 = (accum_2153__field1_6578_3431_10980[i_2151_1302_3429_10982]+var_2155__field1_6580_3433_10984);
accum_2153__field0_6577_3430_10979[i_2151_1302_3429_10982] = var_2156__field0_6581_3435_10985;
accum_2153__field1_6578_3431_10980[i_2151_1302_3429_10982] = var_2156__field1_6582_3436_10986;
}
X_8696_8424[i_2151_1302_3429_10982][0].Re = accum_2153__field0_6577_3430_10979[i_2151_1302_3429_10982];
X_8696_8424[i_2151_1302_3429_10982][0].Im = accum_2153__field1_6578_3431_10980[i_2151_1302_3429_10982];
}
}
#pragma endscop

}
