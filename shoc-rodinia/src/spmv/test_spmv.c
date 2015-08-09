#include <float.h>
#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>

#include "prl.h"
#include "measure-time.h"

#define PAD_FACTOR 16
#define DATATYPE double
#define IS_EQUAL(a,b) (fabs((a) - (b)) <= DBL_EPSILON * fabs(a))

void initRandomMatrix(int *cols, int *rowDelimiters, const int n, const int dim);
void convertToColMajor(DATATYPE *A, int *cols, int dim, int *rowDelimiters,
                       DATATYPE *newA, int *newcols, int *rl, int maxrl,
                       int padded);

void pencil_csr_scalar_dp(int nRows, int nNonzero,
        int rowDelims[restrict const static nRows+1],
        int columns[restrict const static nNonzero],
        double mat[restrict const static nNonzero],
        double vec[restrict const static nRows],
        double out[restrict const static nRows]);

void pencil_ellpackr_dp(int nRows, int maxRowLength,
        int rowLengths[restrict const static nRows],
        int columns[restrict const static maxRowLength][nRows],
        double mat[restrict const static maxRowLength][nRows],
        double vec[restrict const static nRows],
        double out[restrict const static nRows]);

void ref_csr_scalar_dp(int nRows, int nNonzero,
        int rowDelims[restrict const static nRows+1],
        int columns[restrict const static nNonzero],
        double mat[restrict const static nNonzero],
        double vec[restrict const static nRows],
        double out[restrict const static nRows])
{
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
}

void ref_ellpackr_dp(int nRows, int maxRowLength,
        int rowLengths[restrict const static nRows],
        int columns[restrict const static maxRowLength][nRows],
        double mat[restrict const static maxRowLength][nRows],
        double vec[restrict const static nRows],
        double out[restrict const static nRows])
{
    for (int i = 0; i < nRows; i++) {
        int rowLength = rowLengths[i];
        double accum = 0.0;

        for (int j = 0; j < rowLength; j++) {
            int col = columns[j][i];

            accum += mat[j][i] * vec[col];
        }
        out[i] = accum;
    }
}

// Original version by: Kyle Spafford
void initRandomMatrix(int *cols, int *rowDelimiters, const int n, const int dim)
{
    int nnzAssigned = 0;

    // Figure out the probability that a nonzero should be assigned to a given
    // spot in the matrix
    double prob = (double)n / ((double)dim * (double)dim);

    // Seed random number generator
    srand48(8675309L);

    // Randomly decide whether entry i,j gets a value, but ensure n values
    // are assigned
    int fillRemaining = 0;
    int i,j;
    for (i = 0; i < dim; i++)
    {
        rowDelimiters[i] = nnzAssigned;
        for (j = 0; j < dim; j++)
        {
            int numEntriesLeft = (dim * dim) - ((i * dim) + j);
            int needToAssign   = n - nnzAssigned;
            if (numEntriesLeft <= needToAssign) {
                fillRemaining = 1;
            }
            if ((nnzAssigned < n && drand48() <= prob) || fillRemaining)
            {
                // Assign (i,j) a value
                cols[nnzAssigned] = j;
                nnzAssigned++;
            }
        }
    }
    // Observe the convention to put the number of non zeroes at the end of the
    // row delimiters array
    rowDelimiters[dim] = n;
    assert(nnzAssigned == n);
}

// Original version by: Lukasz Wesolowski
void convertToColMajor(DATATYPE *A, int *cols, int dim, int *rowDelimiters,
                       DATATYPE *newA, int *newcols, int *rl, int maxrl,
                       int padded)
{
    int pad = 0;
    if (padded && dim % PAD_FACTOR != 0)
    {
        pad = PAD_FACTOR - dim % PAD_FACTOR;
    }

    int newIndex = 0;
    int i,j,p;
    for (j=0; j<maxrl; j++)
    {
        for (i=0; i<dim; i++)
        {
            if (rowDelimiters[i] + j < rowDelimiters[i+1])
            {
                newA[newIndex] = A[rowDelimiters[i]+j];
                newcols[newIndex] = cols[rowDelimiters[i]+j];
            }
            else
            {
                newA[newIndex] = 0;
            }
            newIndex++;
        }
        if (padded)
        {
            for (p=0; p<pad; p++)
            {
                newA[newIndex] = 0;
                newIndex++;
            }
        }
    }
}



enum Benchmark {
    bm_csr,
    bm_ellpackr
};

struct Options {
    /* Number of passes. */
    int nPasses;

    /* Matrix size. */
    int size;

    /* Benchmark version to run. */
    enum Benchmark benchmark;
};

static struct Options options;

/* Print vector to stdout.
 */
void printVector(int n, DATATYPE *array) {
    for (int i = 0; i < n; i++) {
        printf("%f\n", array[i]);
    }
}

/* Compare two vectors. */
void validateVector(const int dim, DATATYPE expected[restrict const static dim], DATATYPE actual[restrict const static dim]) {
    for (int i = 0; i < dim; i++) {
        if (!IS_EQUAL(expected[i], actual[i])) {
            fprintf(stderr, "Validation error at element %d: ", i);
            fprintf(stderr, "expected %f", expected[i]);
            fprintf(stderr, "   found %f", actual[i]);
            fprintf(stderr, "\n");
        }
    }
}

/* Fill an array with random values.
 */
void fillArray(int n, DATATYPE *array) {
    for (int i = 0; i < n; i++) {
        array[i] = (DATATYPE)(10.0 * rand() / (RAND_MAX + 1.0f));
    }
}

/* Generate a square matrix in CSR format with the given number of rows and 1% nonzero entries.
 */
void generateCsrMatrix(int nRows, int *nNonzero, int **rowDelims, int **columns, DATATYPE **values) {
    *nNonzero = nRows * nRows / 100;    // As with SHOC, make 1% of matrix entries nonzero
    *values = prl_alloc(*nNonzero * sizeof(DATATYPE));
    *rowDelims = prl_alloc((nRows+1) * sizeof(int));
    *columns = prl_alloc(*nNonzero * sizeof(int));
    fillArray(*nNonzero, *values);
    initRandomMatrix(*columns, *rowDelims, *nNonzero, nRows);
}

/* Compute row lengths for a CSR matrix and return the length of the widest row.
 */
int getRowLengths(int nRows, int rowDelims[restrict const static nRows+1], int rowLengths[restrict const static nRows]) {
    int maxLength = 0;
    for (int i = 0; i < nRows; i++) {
        rowLengths[i] = rowDelims[i+1] - rowDelims[i];
        if (rowLengths[i] > maxLength) {
            maxLength = rowLengths[i];
        }
    }
    return maxLength;
}

/* Generate an ELLPACK-R matrix in column-major order.
 *
 * First generates a CSR matrix and then transforms it to a column-major
 * ELLPACK-R matrix.
 */
void generateEllpackrColMajMatrix(int nRows, int *nNonzero, int *maxRowLength, int **rowLengths, int **columns, DATATYPE **values) {
    int *csrRowDelims;
    int *csrColumns;
    DATATYPE *csrValues;
    generateCsrMatrix(nRows, nNonzero, &csrRowDelims, &csrColumns, &csrValues);

    *rowLengths = prl_alloc(nRows * sizeof(int));
    *maxRowLength = getRowLengths(nRows, csrRowDelims, *rowLengths);

    *columns = prl_alloc(nRows * *maxRowLength * sizeof(int));
    *values = prl_alloc(nRows * *maxRowLength * sizeof(DATATYPE));
    convertToColMajor(csrValues, csrColumns, nRows, csrRowDelims, *values, *columns, *rowLengths, *maxRowLength, 0);

    prl_free(csrRowDelims);
    prl_free(csrColumns);
    prl_free(csrValues);
}

/* Test SpMV on a CSR matrix.
 */
void testCsr(int nRows) {
    int *rowDelims;
    int *columns;
    DATATYPE *values;
    int nNonzero;
    generateCsrMatrix(nRows, &nNonzero, &rowDelims, &columns, &values);

    // Get an input vector with random values
    DATATYPE *vector = prl_alloc(nRows * sizeof(DATATYPE));
    fillArray(nRows, vector);

    // Fill reference output with random values
    DATATYPE *refOut = malloc(nRows * sizeof(DATATYPE));
    fillArray(nRows, refOut);

    // Fill PPCG output with random values
    DATATYPE *ppcgOut = prl_alloc(nRows * sizeof(DATATYPE));
    fillArray(nRows, ppcgOut);

/*    START_MEASURE_TIME_HOST
    ref_csr_scalar_dp(nRows, nNonzero, rowDelims, columns, values, vector, refOut);
    STOP_MEASURE_TIME_HOST
    PRINT_TIME_HOST*/

    START_MEASURE_TIME
    for (int i = 0; i < options.nPasses; i++) {
        pencil_csr_scalar_dp(nRows, nNonzero, rowDelims, columns, values, vector, ppcgOut);
    }
    STOP_MEASURE_TIME

    validateVector(nRows, refOut, ppcgOut);
    PRINT_TIME

    prl_free(values);
    prl_free(rowDelims);
    prl_free(columns);
    prl_free(vector);
    free(refOut);
    prl_free(ppcgOut);
}

/* Test SpMV on an ELLPACK-R matrix.
 */
void testEllpackR(int nRows) {
    int nNonzero;
    int maxRowLength;
    int *rowLengths;
    int *columns;
    DATATYPE *values;
    generateEllpackrColMajMatrix(nRows, &nNonzero, &maxRowLength, &rowLengths, &columns, &values);

    // Get an input vector with random values
    DATATYPE *vector = prl_alloc(nRows * sizeof(DATATYPE));
    fillArray(nRows, vector);

    // Fill reference output with random values
    DATATYPE *refOut = malloc(nRows * sizeof(DATATYPE));
    fillArray(nRows, refOut);

    // Fill PPCG output with random values
    DATATYPE *ppcgOut = prl_alloc(nRows * sizeof(DATATYPE));
    fillArray(nRows, ppcgOut);

    START_MEASURE_TIME_HOST
    ref_ellpackr_dp(nRows, maxRowLength, rowLengths, columns, values, vector, refOut);
    STOP_MEASURE_TIME_HOST
    PRINT_TIME_HOST

    START_MEASURE_TIME
    for (int i = 0; i < options.nPasses; i++) {
        pencil_ellpackr_dp(nRows, maxRowLength, rowLengths, columns, values, vector, ppcgOut);
    }
    STOP_MEASURE_TIME

    validateVector(nRows, refOut, ppcgOut);
    PRINT_TIME

    prl_free(values);
    prl_free(rowLengths);
    prl_free(columns);
    prl_free(vector);
    free(refOut);
    prl_free(ppcgOut);
}

void printUsage(const char *name) {
    fprintf(stderr, "Usage: %s [-i iters] [-n passes] [-s size] csr|ellpackr\n", name);
}

void parseCommandLine(int argc, char *argv[]) {
    options.nPasses = 1;
    options.size = 16384;

    int opt;
    while ((opt = getopt(argc, argv, "i:n:s:")) != -1) {
        switch (opt) {
        case 'i':
            fprintf(stderr, "Varying the number of kernel invocations is not supported.\n");
            exit(EXIT_FAILURE);
            break;
        case 'n':
            options.nPasses = atoi(optarg);
            break;
        case 's':
            options.size = atoi(optarg);
            break;
        default:
            printUsage(argv[0]);
            exit(EXIT_FAILURE);
        }
    }

    if (optind >= argc) {
        fprintf(stderr, "Insufficient number of arguments.\n");
        printUsage(argv[0]);
        exit(EXIT_FAILURE);
    }

    if (strcmp(argv[optind], "csr") == 0) {
        options.benchmark = bm_csr;
    }
    else if (strcmp(argv[optind], "ellpackr") == 0) {
        options.benchmark = bm_ellpackr;
    }
    else {
        fprintf(stderr, "Invalid benchmark '%s'\n", argv[optind]);
        exit(EXIT_FAILURE);
    }
}

int main(int argc, char *argv[]) {
    parseCommandLine(argc, argv);

    if (options.benchmark == bm_csr) {
        testCsr(options.size);
    }
    else if (options.benchmark == bm_ellpackr) {
        testEllpackR(options.size);
    }

    return EXIT_SUCCESS;
}
