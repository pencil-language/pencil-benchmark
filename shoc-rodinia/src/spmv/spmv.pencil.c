#include <pencil.h>

void pencil_csr_scalar_dp(int nRows, int nNonzero,
        int rowDelims[restrict const static nRows+1],
        int columns[restrict const static nNonzero],
        double mat[restrict const static nNonzero],
        double vec[restrict const static nRows],
        double out[restrict const static nRows])
{
#pragma scop
    __pencil_assume(nRows > 0);
    __pencil_assume(nNonzero > 0);

#pragma pencil independent
    for (int i = 0; i < nRows; i++) {
        int start = rowDelims[i];
        int end = rowDelims[i+1];
        double accum = 0.0;

        for (int j = start; j < end; j++) {
            int col = columns[j];

            accum += mat[j] * vec[col];
        }
        out[i] = accum;
    }
#pragma endscop
}

void pencil_ellpackr_dp(int nRows, int maxRowLength,
        int rowLengths[restrict const static nRows],
        int columns[restrict const static maxRowLength][nRows],
        double mat[restrict const static maxRowLength][nRows],
        double vec[restrict const static nRows],
        double out[restrict const static nRows])
{
#pragma scop
    __pencil_assume(nRows > 0);
    __pencil_assume(maxRowLength > 0);

#pragma pencil independent
    for (int i = 0; i < nRows; i++) {
        int rowLength = rowLengths[i];
        double accum = 0.0;

        for (int j = 0; j < rowLength; j++) {
            int col = columns[j][i];

            accum += mat[j][i] * vec[col];
        }
        out[i] = accum;
    }
#pragma endscop
}
