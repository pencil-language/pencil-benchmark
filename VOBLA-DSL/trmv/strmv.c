#include "../include/vobla_headers.h"

#define N 1024
#define DATA_TYPE float
#define DATA_PRINTF_FORMAT "%f "

void print_array(int n, int m, DATA_TYPE B[n][m])
{
	int i, j;
	for (i=0; i<n; i++)
  	  for (j=0; j<m; j++)
		fprintf(DATA_OUTPUT_FILE, DATA_PRINTF_FORMAT, B[i][j]);
}

void init_data(int n, int m, DATA_TYPE B[n][m])
{
	int i, j;
	for (i=0; i<n; i++)
   	  for (j=0; j<m; j++)
		B[i][j] = 2.0;
}

int main()
{
	double start_time, end_time;

	DATA_TYPE (*A)[N][N] = (DATA_TYPE (*)[N][N]) malloc(N*N * sizeof(DATA_TYPE));
	DATA_TYPE (*C)[N][1] = (DATA_TYPE (*)[N][1]) malloc(N*1 * sizeof(DATA_TYPE));

        init_data(N, N, *A);
	init_data(N, 1, *C);


	start_time = get_time();
	pencil_strmv_n_nun(N, 10.0, *A, 1, *C);
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
void pencil_strmv_n_nun (int n_7988_8196,int ldAT_7989_8197,float AT_7990_8198[restrict const static n_7988_8196][n_7988_8196],int incX_7991_8199,float X_7992_8200[restrict const static n_7988_8196][incX_7991_8199]){

#pragma scop
{
for(int i_1992_4008_7662_10587 = 0;i_1992_4008_7662_10587 <= (n_7988_8196-1);i_1992_4008_7662_10587 += 1){
float accum_1994_4010_7663_10586[n_7988_8196];
accum_1994_4010_7663_10586[i_1992_4008_7662_10587] = 0.0f;
for(int i_1995_4011_7664_10588 = 0;i_1995_4011_7664_10588 <= ((((n_7988_8196-1)-i_1992_4008_7662_10587)+1)-1);i_1995_4011_7664_10588 += 1){
accum_1994_4010_7663_10586[i_1992_4008_7662_10587] = (accum_1994_4010_7663_10586[i_1992_4008_7662_10587]+(AT_7990_8198[i_1992_4008_7662_10587][(i_1992_4008_7662_10587+i_1995_4011_7664_10588)]*X_7992_8200[(i_1992_4008_7662_10587+i_1995_4011_7664_10588)][0]));
}
X_7992_8200[i_1992_4008_7662_10587][0] = accum_1994_4010_7663_10586[i_1992_4008_7662_10587];
}
}
#pragma endscop

}
