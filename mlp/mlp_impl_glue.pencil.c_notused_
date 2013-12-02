// -*- c++ -*-
// UjoImro, 2013
// Experimental code for the CARP Project
// Copyright (c) RealEyes, 2013
// This is a c-implementation of the PCA->MLP response map calculation

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include "pencil.h"

#include "mlp_impl_pencil.h"

static void copyMatFloatToArray(MatFloat In, int n, int m, float Out[][m]) {
    // assert(In.rows == n);
    // assert(In.cols == m);
#pragma scop
    for (int i = 0; i < In.rows; i++)
        for (int j = 0; j < In.cols; j++)
            Out[i][j] = In.data[In.start + i * In.step + j];
#pragma endscop
    return;
}

static void copyMatCharToArray(MatChar In, int n, int m, uint8_t Out[][m]) {
    // assert(In.rows == n);
    // assert(In.cols == m);
#pragma scop
    for (int i = 0; i < In.rows; i++)
        for (int j = 0; j < In.cols; j++)
            Out[i][j] = In.data[In.start + i * In.step + j];
#pragma endscop
    return;    
}

static void copyArrayToMatFloat(int n, int m, float In[][m], MatFloat Out) {
    // assert(Out.rows == n);
    // assert(Out.cols == m);
    // assert(Out.step == m);
#pragma scop
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            Out.data[Out.start + i * Out.step + j] = In[i][j];
#pragma endscop
    return;    
}

static void copyArrayToMatChar(int n, int m, uint8_t In[][m], MatChar Out) {
    // assert(Out.rows == n);
    // assert(Out.cols == m);
    // assert(Out.step == m);
#pragma scop
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            Out.data[Out.start + i * Out.step + j] = In[i][j];
#pragma endscop
    return;    
}

void printArray(int n, int m, float Array[n][m]) {
    printf("%s = [\n", "NAME");


    for (int q = 0; q < n; q++) {
        printf("[ ");
        for (int w = 0; w < m; w++) {
            printf("%.40f, ", Array[q][w]);
        }
        printf(" ]\n");
    }

    printf("]\n");
    return;
}

void freeMLP(mlp *classifier) {
#pragma scop
    freeMatFloat(&classifier->m_wIn);
    freeMatFloat(&classifier->m_wOut);
    freeMatFloat(&classifier->m_U);
#pragma endscop
    return;    
}

MatFloat CreateMatFloat(int rows, int cols) {
#pragma scop
    MatFloat result;
    // assert(rows > 0);
    // assert(cols > 0);

    result.data = NULL;
    result.data = (float *)malloc(sizeof(float) * rows * cols);
    // assert(result.data);
    result.rows = rows;
    result.cols = cols;
    result.step = cols;
    result.start = 0;

#pragma endscop
    return result;
}

MatChar CreateMatChar(int rows, int cols) {
#pragma scop
    MatChar result;
    // assert(rows > 0);
    // assert(cols > 0);

    result.data = NULL;
    result.data = (uint8_t *)malloc(sizeof(uint8_t) * rows * cols);
    // assert(result.data);
    result.rows = rows;
    result.cols = cols;
    result.step = cols;
    result.start = 0;

#pragma endscop
    return result;
}

static void convertFromCharToFloatArray(int imageRows, int imageCols,
                                        uint8_t In[static const restrict imageRows][imageCols],
					int imageOffsetRow, int imageOffsetCol,
                                        float quotient, float shift,
                                        int OutRows, int OutCols,
                                        float Out[static const restrict OutRows][OutCols]) {

#pragma scop
    for (int i = 0; i < OutRows; i++)
        for (int j = 0; j < OutCols; j++)
            Out[i][j] = quotient * (float)In[i + imageOffsetRow][j + imageOffsetCol] + shift;
#pragma endscop
    return;
}

void freeMatFloat(MatFloat *mat) {
#pragma scop
    // assert(mat->data);
    // assert(mat->start == 0);
    free(mat->data);
    mat->data = NULL;
    mat->rows = 0;
    mat->cols = 0;
    mat->step = 0;
    mat->start = 0;
#pragma endscop
    return;
}

void freeMatChar(MatChar *mat) {
#pragma scop
    // assert(mat->data);
    // assert(mat->start == 0);
    free(mat->data);
    mat->data = NULL;
    mat->rows = 0;
    mat->cols = 0;
    mat->step = 0;
    mat->start = 0;
#pragma endscop
    return;

}

// returns alpha*A*B + beta * C
    static void gemmFloatArray(int ARows, int ACols,
        float A[static const restrict ARows][ACols],
        int BRows, int BCols,
        float B[static const restrict BRows][BCols],
        float alpha, int CRows, int CCols,
        float C[static const restrict CRows][CCols],
        float beta, int ResRows,
        int ResCols,
        float Res[static const restrict ResRows][ResCols]) {
    // assert(ARows == CRows);
    // assert(ACols == BRows);
    // assert(BCols == CCols);
    // assert(CRows == ResRows);
    // assert(CCols == ResCols)
#pragma scop
    for (int i = 0; i < CRows; i++)
        for (int j = 0; j < CCols; j++) {
    Res[i][j] = beta * C[i][j];
    for (int k = 0; k < ACols; k++) {
    Res[i][j] += alpha * A[i][k] * B[k][j];
}
}
#pragma endscop
    return;
}

static void expFloat(int rows, int cols,
                     float Mat[static const restrict rows][cols]) {
#pragma scop
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            Mat[i][j] = expf(Mat[i][j]);
#pragma endscop
    return;
}

    static void addFloat(int rows, int cols,
        float Mat[static const restrict rows][cols], float Val) {
#pragma scop
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            Mat[i][j] = Val + Mat[i][j];
#pragma endscop
    return;
}

static void divideFloat(float Val, int rows, int cols,
                        float Mat[static const restrict rows][cols]) {

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            Mat[i][j] = Val / Mat[i][j];
#pragma endscop
    return;
}

    static float dotProduct(int LeftRows, int LeftCols, int RightRows,
        int RightCols, float Left[const restrict][LeftCols],
        float Right[const restrict][RightCols]) {
#pragma scop
    // assert(LeftCols == 1);
    // assert(RightCols == 1);
    // assert(LeftRows == RightRows);

    float Result = 0.;

    for (int i = 0; i < LeftRows; i++) {
    Result += Left[i][0] * Right[i][0];
}

#pragma endscop
    return Result;

}

    void calculateRespondMaps(
    int m_visibleLandmarks_size, int MapSize, int ImageRows, int ImageCols,
        uint8_t Image[static const restrict ImageRows][ImageCols], int shape_rows,
        int shape_cols, float shape_data[static const restrict shape_rows][shape_cols],
        int shape_start, mlp m_classifiers[static const restrict m_visibleLandmarks_size],
        float ResponseMaps[static const restrict m_visibleLandmarks_size][MapSize + MapSize + 1][MapSize + MapSize + 1]);

/// @brief Calculate a set of response maps.
///
/// This function implements the external interface. In its implementation
/// we transform the input parameters in a way, such that the internal
/// implementation is as close to PENCIL as possible.
    void calculateMaps(int NumLandMarks, int MapSize, MatChar Image, MatFloat Shape,
        mlp Classifiers[], MatFloat *MatResponseMaps[]) {
#pragma scop

    int Width = MapSize * 2 + 1;
    float (*ResponseMaps)[Width][Width] =
        malloc(sizeof(float) * NumLandMarks * Width * Width + 1);

    int ImageRows = Image.rows;
    int ImageCols = Image.cols;
    uint8_t (*ImageArray)[ImageCols] =
        malloc(sizeof(uint8_t) * ImageRows * ImageCols);
    copyMatCharToArray(Image, ImageRows, ImageCols, ImageArray);

    calculateRespondMaps(NumLandMarks, MapSize, ImageRows, ImageCols, ImageArray,
        Shape.step, Shape.cols, (float (*) [Shape.cols]) Shape.data,
        Shape.start, Classifiers, ResponseMaps);

    for (int i = 0; i < NumLandMarks; ++i)
        copyArrayToMatFloat(Width, Width, ResponseMaps[i],
        *(&(*MatResponseMaps)[i]));

    free(ImageArray);
    free(ResponseMaps);
#pragma endscop
    return;
    
}
