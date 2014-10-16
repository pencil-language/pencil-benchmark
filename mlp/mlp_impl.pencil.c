// -*- c++ -*-
// UjoImro, 2013
// Experimental code for the CARP Project
// Copyright (c) RealEyes, 2013
// This is a c-implementation of the PCA->MLP response map calculation

#include "mlp_impl_pencil.h"

#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <pencil.h>

static void copyMatCharToArray(MatChar In, int n, int m, uint8_t Out[][m]) {
    assert(In.rows == n);
    assert(In.cols == m);

    for (int i = 0; i < In.rows; i++)
        for (int j = 0; j < In.cols; j++)
            Out[i][j] = In.data[In.start + i * In.step + j];
}

static void copyArrayToMatFloat(int n, int m, float In[][m], MatFloat Out) {
    assert(Out.rows == n);
    assert(Out.cols == m);
    assert(Out.step == m);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            Out.data[Out.start + i * Out.step + j] = In[i][j];
}

void freeMLP(mlp *classifier) {
    freeMatFloat(&classifier->m_wIn);
    freeMatFloat(&classifier->m_wOut);
    freeMatFloat(&classifier->m_U);
}

MatFloat CreateMatFloat(int rows, int cols) {
    MatFloat result;
    assert(rows > 0);
    assert(cols > 0);

    result.data = (float *)malloc(sizeof(float) * rows * cols);
    assert(result.data);
    result.rows = rows;
    result.cols = cols;
    result.step = cols;
    result.start = 0;

    return result;
}

MatChar CreateMatChar(int rows, int cols) {
    MatChar result;
    assert(rows > 0);
    assert(cols > 0);

    result.data = (uint8_t *)malloc(sizeof(uint8_t) * rows * cols);
    assert(result.data);
    result.rows = rows;
    result.cols = cols;
    result.step = cols;
    result.start = 0;

    return result;
}

static void transposeFloat( int InRows, int InCols, float In[InRows][InCols]
                          , int OutRows, int OutCols, float Out[OutRows][OutCols])
{
#pragma scop
    __pencil_assume(InRows == OutCols);
    __pencil_assume(OutCols == InRows);

    for (int i = 0; i < InRows; i++)
        for (int j = 0; j < InCols; j++)
            Out[j][i] = In[i][j];
#pragma endscop
}

static float meanChar(int subImageRows, int subImageCols, int imageRows,
		      int imageCols, uint8_t image[imageRows][imageCols],
		      int imageOffsetRow, int imageOffsetCol) {
#pragma scop
    float sum = 0;

    for (int i = 0; i < subImageRows; i++) {
        for (int j = 0; j < subImageCols; j++) {
            sum += image[i + imageOffsetRow][j + imageOffsetCol];
        }
    }
#pragma endscop
    return sum / (subImageRows * subImageCols);
}

static uint8_t minChar(int subImageRows, int subImageCols, int imageRows,
		       int imageCols, uint8_t image[imageRows][imageCols],
		       int imageOffsetRow, int imageOffsetCol) {
#pragma scop
    uint8_t minvalue = 255;

    for (int i = 0; i < subImageRows; i++)
        for (int j = 0; j < subImageCols; j++)
            minvalue = min(minvalue, image[i + imageOffsetRow][j+imageOffsetCol]);
#pragma endscop
    return minvalue;
}

static uint8_t maxChar(int subImageRows, int subImageCols, int imageRows,
		       int imageCols, uint8_t image[imageRows][imageCols],
		       int imageOffsetRow, int imageOffsetCol) {
#pragma scop
    uint8_t maxvalue = 0;

    for (int i = 0; i < subImageRows; i++)
        for (int j = 0; j < subImageCols; j++)
            maxvalue = max(maxvalue, image[i + imageOffsetRow][j+imageOffsetCol]);
#pragma endscop
    return maxvalue;
}

void freeMatFloat(MatFloat *mat) {
    assert(mat->data);
    assert(mat->start == 0);
    free(mat->data);
    mat->data = NULL;
    mat->rows = 0;
    mat->cols = 0;
    mat->step = 0;
    mat->start = 0;
}

void freeMatChar(MatChar *mat) {
    assert(mat->data);
    assert(mat->start == 0);
    free(mat->data);
    mat->data = NULL;
    mat->rows = 0;
    mat->cols = 0;
    mat->step = 0;
    mat->start = 0;
}

// returns alpha*A*B + beta * C
static void gemmFloatArray(int ARows, int ACols, float A[ARows][ACols],
                           int BRows, int BCols, float B[BRows][BCols],
                           float alpha, int CRows, int CCols,
                           float C[CRows][CCols], float beta, int ResRows,
                           int ResCols, float Res[ResRows][ResCols]) {
#pragma scop
    __pencil_assume(ARows == CRows);
    __pencil_assume(ACols == BRows);
    __pencil_assume(BCols == CCols);
    __pencil_assume(CRows == ResRows);
    __pencil_assume(CCols == ResCols);

    for (int i = 0; i < CRows; i++) {
        for (int j = 0; j < CCols; j++) {
            Res[i][j] = beta * C[i][j];
            for (int k = 0; k < ACols; k++) {
                Res[i][j] += alpha * A[i][k] * B[k][j];
            }
        }
    }
#pragma endscop
}

float GetValueFloat(MatFloat self, int row, int col) {
    return self.data[row * self.step + col + self.start];
}

static void copySubArrayFloat(int arrayRows, int arrayCols,
                              float Array[arrayRows][arrayCols],
                              int subArrayRows, int subArrayCols,
                              float subArray[subArrayRows][subArrayCols],
                              int offsetRow, int offsetCol) {
#pragma scop
    for (int i = 0; i < subArrayRows; i++)
        for (int j = 0; j < subArrayCols; j++)
            subArray[i][j] = Array[i + offsetRow][j + offsetCol];
#pragma endscop
}

static float generateResponseMapPatch(
    int mapSize, int ncx, int ncy, int bInRows, int bInCols,
    float bInArray[bInRows][bInCols], mlp classifier, int ImageRows,
    int ImageCols, uint8_t Image[ImageRows][ImageCols], Point2i center,
    int wInRows, int wInCols, float wInArray[wInRows][wInCols], int wOutRows,
    int wOutCols, float wOutArray[wOutRows][wOutCols], float bOut) {
#pragma scop
    int cy = ncy + center.y - mapSize;
    int cx = ncx + center.x - mapSize;

    int imagePatchRows = 2 * classifier.m_patchSize + 1; // m_patchSize is always 5
    int imagePatchCols = 2 * classifier.m_patchSize + 1;

    int imageOffsetRow = cy - classifier.m_patchSize;
    int imageOffsetCol = cx - classifier.m_patchSize;
    float sampleMean = meanChar(imagePatchRows, imagePatchCols, ImageRows,
                                ImageCols, Image, imageOffsetRow, imageOffsetCol);
    float sampleMin = minChar(imagePatchRows, imagePatchCols, ImageRows,
                              ImageCols, Image, imageOffsetRow, imageOffsetCol);
    float sampleMax = maxChar(imagePatchRows, imagePatchCols, ImageRows,
                              ImageCols, Image, imageOffsetRow, imageOffsetCol);

    sampleMax -= sampleMean;
    sampleMin -= sampleMean;

    sampleMax = fmaxf(fabsf(sampleMin), fabsf(sampleMax));

    if (sampleMax == 0.0)
        sampleMax = 1.0;

    float quotient = 1.0f / sampleMax;
    float shift = -(1.0f / sampleMax) * sampleMean;

    float alpha = -1.0;
    float beta = -1.0;
    float result = 0;

    for (int i = 0; i < bInRows; i++) {
        for (int j = 0; j < bInCols; j++) { // This loop seems to have a single
            // iteration? Is this always true?
            float xOutArray;
            xOutArray = beta * bInArray[i][j];
            for (int k = 0; k < wInCols; k++) {
		int index1 = k / imagePatchCols + imageOffsetRow;
		int index2 = k % imagePatchRows + j + imageOffsetCol;
                xOutArray += alpha * wInArray[i][k] *
                    (quotient * Image[index1][index2] +
                     shift);
            }
            xOutArray = exp(xOutArray);	//This should be expf in C and exp in OpenCL
            xOutArray = xOutArray + 1.0f;
            xOutArray = 2.0f / xOutArray;
            xOutArray = xOutArray + -1.0f;
            result += wOutArray[i][0] * xOutArray;
        }
    }

    result = - result;
    result -= bOut;
    result = 1.0f / (1.0f + exp(result)); //This should be expf in C and exp in OpenCL
#pragma endscop
    return result;
}

/// @brief Calculate a single response map for an image.
///
/// @param Image The image to process.
static void generateResponseMap(
    int ImageRows, int ImageCols, uint8_t Image[ImageRows][ImageCols],
    const Point2i center, int mapSize, mlp classifier,
    float ResponseMap[mapSize + mapSize + 1][mapSize + mapSize + 1]) {

    // Translate input arrays into C99 Arrays.
    //
    // The problem here is that the size of the arrays is not know to
    // be identical for each respond map calculation. Hence, we can
    // not easily move those allocations out of the core computation.
    __pencil_assume(classifier.m_U.start == 0);
    __pencil_assume(classifier.m_U.cols == classifier.m_U.step);
    int m_URows = classifier.m_U.rows; // Always 121
    int m_UCols = classifier.m_U.cols; // Varies between 17 and 30
    float (*m_UArray)[m_UCols] = (void*)classifier.m_U.data;


    // Translate input arrays into C99 Arrays
    __pencil_assume(classifier.m_wIn.start == 0);
    __pencil_assume(classifier.m_wIn.cols == classifier.m_wIn.step);
    int m_wInRows = classifier.m_wIn.rows; // Always 25
    int m_wInCols = classifier.m_wIn.cols; // Varies between 17 and 31
    float (*m_wInArray)[m_wInCols] = (void*)classifier.m_wIn.data;

    // Translate input arrays into C99 Arrays
    __pencil_assume(classifier.m_wOut.start == 0);
    __pencil_assume(classifier.m_wOut.cols == classifier.m_wOut.step);
    int m_wOutRows = classifier.m_wOut.rows; // Always 1
    int m_wOutCols = classifier.m_wOut.cols; // Always 26
    float (*m_wOutArray)[m_wOutCols] = (void*)classifier.m_wOut.data;

    // This is a temporary array.
    //
    // Also, I wonder if this computation can not be precalculated? Or that
    // we could keep track of the index dimensions being transposed and we
    // just access the elements at their previous locations?
    //
    // Also, the size of this array varies / is data-dependent.
    int m_U_transposeRows = m_UCols;
    int m_U_transposeCols = m_URows;
    float (*m_U_transposeArray)[m_U_transposeCols] =
        malloc(sizeof(float) * m_U_transposeRows * m_U_transposeCols);

    transposeFloat(m_URows, m_UCols, m_UArray, m_U_transposeRows,
                   m_U_transposeCols, m_U_transposeArray);

    // A sub array.
    //
    // Instead of explicitly calculating this, it would be nice if ppcg
    // could just be told that we only access a subset of this array.
    //
    // Also the size of this array is data-dependent
    int wIn_ARows = m_wInRows;
    int wIn_ACols = m_wInCols - 1;
    float (*wIn_AArray)[wIn_ACols] = malloc(sizeof(float) * wIn_ARows * wIn_ACols);
    copySubArrayFloat(m_wInRows, m_wInCols, m_wInArray, m_wInRows,
                      m_wInCols - 1, wIn_AArray, 0, 0);

    // A temporary array.
    //
    // Why again can this not be precomputed?
    //
    // The size of this array seems to be constant for all test cases.
    int wInRows = wIn_ARows;
    int wInCols = m_U_transposeCols;
    float (*wInArray)[wInCols] = malloc(sizeof(float) * wInRows * wInCols);
    gemmFloatArray(wIn_ARows, wIn_ACols, wIn_AArray, m_U_transposeRows,
                   m_U_transposeCols, m_U_transposeArray, 1.0, wInRows, wInCols,
                   wInArray, 0.0, wInRows, wInCols, wInArray);

    // A sub array.
    //
    // Instead of explicitly calculating this, it would be nice if ppcg
    // could just be told that we only access a subset of this array.
    int bInRows = m_wInRows;
    int bInCols = 1;
    float (*bInArray)[bInCols] = malloc(sizeof(float) * bInRows * bInCols);
    copySubArrayFloat(m_wInRows, m_wInCols, m_wInArray, bInRows,
                      bInCols, bInArray, 0, m_wInCols - 1);


    // A sub array.
    //
    // Instead of explicitly calculating this, it would be nice if ppcg
    // could just be told that we only access a subset of this array.
    int wOut_tmpRows = m_wOutRows;
    int wOut_tmpCols = m_wOutCols - 1;
    float (*wOut_tmpArray)[wOut_tmpCols] = malloc(sizeof(float) * wOut_tmpRows * wOut_tmpCols);
    copySubArrayFloat(m_wOutRows, m_wOutCols, m_wOutArray, wOut_tmpRows,
                      wOut_tmpCols, wOut_tmpArray, 0, 0);

    int wOutRows = wOut_tmpCols;
    int wOutCols = wOut_tmpRows;
    float (*wOutArray)[wOutCols] = malloc(sizeof(float) * wOutRows * wOutCols);
    transposeFloat(wOut_tmpRows, wOut_tmpCols, wOut_tmpArray, wOutRows, wOutCols,
                   wOutArray);

    float bOut = m_wOutArray[0][m_wOutCols - 1];

    for (int ncy = 0; ncy <= 2 * mapSize; ++ncy) {
        for (int ncx = 0; ncx <= 2 * mapSize; ++ncx) {
            ResponseMap[ncy][ncx] = generateResponseMapPatch(
                mapSize, ncx, ncy, bInRows, bInCols, bInArray, classifier,
                ImageRows, ImageCols, Image, center, wInRows, wInCols, wInArray,
                wOutRows, wOutCols, wOutArray, bOut);
        }
    }

    free(wIn_AArray);
    free(m_U_transposeArray);
    free(wInArray);
    free(wOut_tmpArray);
    free(wOutArray);
}

static int cvRound(float value) {
    return (int)(value + (value >= 0 ? 0.5 : -0.5));
}

// Calculate a set of respone maps for a given image.
//
// @param Image The image to process
// @param responseMaps A pointer at which results of the calculation,
//                     the response maps, are stored.
void calculateRespondMaps(
    int m_visibleLandmarks_size, int MapSize, int ImageRows, int ImageCols,
    uint8_t Image[ImageRows][ImageCols], MatFloat shape, mlp m_classifiers[],
    float ResponseMaps[][MapSize + MapSize + 1][MapSize + MapSize + 1]) {
#pragma scop

    // The response maps calculated in this loop are in general calculated
    // from non-overlapping parts of the image. However, even if those parts
    // would overlap, this loop is still parallel, as memory is either only
    // read from or, if it is written to, each iteration writes into a distinct
    // subarray. When translating this to GPU code, we could map this loop to
    // distinct thread groups.
    for (int i = 0; i < m_visibleLandmarks_size; i++) {

        // This array is interesting as it gives the coordinates of the different
        // subimages that we need to process. It is not trivially polyhedral, as
        // when inlined, it is an indirection array, that makes all loop bounds
        // data-dependent. However, it is read only and basically used as a
        // two-dimensional set of parameters.
        //
        // TODO: I do not see why 'shape' is defined using floats. Also, this
        //       is a trivial case where we could use a C99 array to represent
        //       the matrix.
        //       Finally, Point2i is non-polyhedral but a plain structure. We
        //       could represent it as a 2-element array, but it may actually
        //       be worth supporting such simple structures in ppcg.
        //
        Point2i center;
        float shape_x;
        float shape_y;

        shape_x = GetValueFloat(shape, 2 * i, 0);
        shape_y = GetValueFloat(shape, 2 * i + 1, 0);
        center.x = cvRound(shape_x);
        center.y = cvRound(shape_y);

        generateResponseMap(ImageRows, ImageCols, Image, center, MapSize,
                            m_classifiers[i], ResponseMaps[i]);
    }

#pragma endscop
}

/// @brief Calculate a set of response maps.
///
/// This function implements the external interface. In its implementation
/// we transform the input parameters in a way, such that the internal
/// implementation is as close to PENCIL as possible.
void calculateMaps( int NumLandMarks
                  , int MapSize
                  , MatChar Image
                  , MatFloat Shape
                  , mlp Classifiers[static const restrict NumLandMarks]
                  , MatFloat *MatResponseMaps[static const restrict NumLandMarks]
                  )
{
    int Width = MapSize * 2 + 1;
    float (*ResponseMaps)[Width][Width] =
        malloc(sizeof(float) * NumLandMarks * Width * Width + 1);

    int ImageRows = Image.rows;
    int ImageCols = Image.cols;
    uint8_t (*ImageArray)[ImageCols] =
        malloc(sizeof(uint8_t) * ImageRows * ImageCols);
    copyMatCharToArray(Image, ImageRows, ImageCols, ImageArray);

    calculateRespondMaps(NumLandMarks, MapSize, ImageRows, ImageCols, ImageArray,
                         Shape, Classifiers, ResponseMaps);

    for (int i = 0; i < NumLandMarks; ++i)
        copyArrayToMatFloat(Width, Width, ResponseMaps[i],
                            *(&(*MatResponseMaps)[i]));

    free(ImageArray);
    free(ResponseMaps);
}
