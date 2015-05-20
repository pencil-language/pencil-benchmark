#include "../include/vobla_headers.h"

#define N 1024
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

	DATA_TYPE (*A)[N][N] = (DATA_TYPE (*)[N][N]) malloc(N*N * sizeof(DATA_TYPE));
	DATA_TYPE (*C)[N][1] = (DATA_TYPE (*)[N][1]) malloc(N*1 * sizeof(DATA_TYPE));

        init_data(N, N, *A);
	init_data(N, 1, *C);


	start_time = get_time();
	pencil_ztrmv_n_nun(N, 10.0, *A, 1, *C);
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
void pencil_ztrmv_n_nun (int n_9148_8588,int ldAT_9149_8589,struct ComplexDouble AT_9150_8590[restrict const static n_9148_8588][n_9148_8588],int incX_9151_8591,struct ComplexDouble X_9152_8592[restrict const static n_9148_8588][incX_9151_8591]){

#pragma scop
{
for(int i_2733_3682_7173_11396 = 0;i_2733_3682_7173_11396 <= (n_9148_8588-1);i_2733_3682_7173_11396 += 1){
double accum_2735__field1_7458_7174_11394[n_9148_8588];
double accum_2735__field0_7457_7175_11393[n_9148_8588];
accum_2735__field0_7457_7175_11393[i_2733_3682_7173_11396] = 0.0;
accum_2735__field1_7458_7174_11394[i_2733_3682_7173_11396] = 0.0;
for(int i_2736_3685_7176_11397 = 0;i_2736_3685_7176_11397 <= ((((n_9148_8588-1)-i_2733_3682_7173_11396)+1)-1);i_2736_3685_7176_11397 += 1){
double var_2738__field0_7461_7178_11399;
double var_2737__field1_7460_7180_11398;
double var_2737__field0_7459_7179_11395;
double var_2738__field1_7462_7177_11400;
var_2737__field0_7459_7179_11395 = ((AT_9150_8590[i_2733_3682_7173_11396][(i_2733_3682_7173_11396+i_2736_3685_7176_11397)].Re*X_9152_8592[(i_2733_3682_7173_11396+i_2736_3685_7176_11397)][0].Re)-(AT_9150_8590[i_2733_3682_7173_11396][(i_2733_3682_7173_11396+i_2736_3685_7176_11397)].Im*X_9152_8592[(i_2733_3682_7173_11396+i_2736_3685_7176_11397)][0].Im));
var_2737__field1_7460_7180_11398 = ((AT_9150_8590[i_2733_3682_7173_11396][(i_2733_3682_7173_11396+i_2736_3685_7176_11397)].Im*X_9152_8592[(i_2733_3682_7173_11396+i_2736_3685_7176_11397)][0].Re)+(AT_9150_8590[i_2733_3682_7173_11396][(i_2733_3682_7173_11396+i_2736_3685_7176_11397)].Re*X_9152_8592[(i_2733_3682_7173_11396+i_2736_3685_7176_11397)][0].Im));
var_2738__field0_7461_7178_11399 = (accum_2735__field0_7457_7175_11393[i_2733_3682_7173_11396]+var_2737__field0_7459_7179_11395);
var_2738__field1_7462_7177_11400 = (accum_2735__field1_7458_7174_11394[i_2733_3682_7173_11396]+var_2737__field1_7460_7180_11398);
accum_2735__field0_7457_7175_11393[i_2733_3682_7173_11396] = var_2738__field0_7461_7178_11399;
accum_2735__field1_7458_7174_11394[i_2733_3682_7173_11396] = var_2738__field1_7462_7177_11400;
}
X_9152_8592[i_2733_3682_7173_11396][0].Re = accum_2735__field0_7457_7175_11393[i_2733_3682_7173_11396];
X_9152_8592[i_2733_3682_7173_11396][0].Im = accum_2735__field1_7458_7174_11394[i_2733_3682_7173_11396];
}
}
#pragma endscop

}
