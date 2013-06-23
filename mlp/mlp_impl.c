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

#include "mlp_impl.h"

static void copyMatFloatToArray(MatFloat In, int n, int m, float Out[][m]) {
  assert(In.rows == n);
  assert(In.cols == m);

  for (int i = 0; i < In.rows; i++)
    for (int j = 0; j < In.cols; j++)
      Out[i][j] = In.data[In.start + i * In.step + j];
}

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

static void copyArrayToMatChar(int n, int m, uint8_t In[][m], MatChar Out) {
  assert(Out.rows == n);
  assert(Out.cols == m);
  assert(Out.step == m);

  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      Out.data[Out.start + i * Out.step + j] = In[i][j];
}

void printArray(int n, int m, float Array[n][m]) {
  printf("%s = [\n", "NAME");

  int q, w;

  for (q = 0; q < n; q++) {
    printf("[ ");
    for (w = 0; w < m; w++) {
      printf("%f, ", Array[n][m]);
    }
    printf(" ]\n");
  }

  printf("]\n");

  return;
}

float GetValueFloat(MatFloat self, int row, int col);

void printMatFloat(MatFloat mat, char *name) {
  printf("%s = [\n", name);

  int q, w;

  for (q = 0; q < mat.rows; q++) {
    printf("[ ");
    for (w = 0; w < mat.cols; w++) {
      printf("%f, ", GetValueFloat(mat, q, w));
    }
    printf(" ]\n");
  }

  printf("]\n");

  return;
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

  result.data = NULL;
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

  result.data = NULL;
  result.data = (uint8_t *)malloc(sizeof(uint8_t) * rows * cols);
  assert(result.data);
  result.rows = rows;
  result.cols = cols;
  result.step = cols;
  result.start = 0;

  return result;
}

static void transposeFloat(MatFloat input, MatFloat *output) {
  int q, w;
  assert(output->rows == input.cols);
  assert(output->cols == input.rows);
  for (q = 0; q < input.rows; q++)
    for (w = 0; w < input.cols; w++)
      output->data[w * output->step + q + output->start] =
          input.data[q * input.step + w + input.start];

  return;
}

static float meanChar(int Rows, int Cols, uint8_t M[Rows][Cols]) {
  float sum = 0;

  for (int i = 0; i < Rows; i++)
    for (int j = 0; j < Cols; j++) {
      sum += M[i][j];
    }

  return sum / (Rows * Cols);
}

static uint8_t min(uint8_t a, uint8_t b) { return a < b ? a : b; }

static uint8_t minChar(int Rows, int Cols, uint8_t M[Rows][Cols]) {
  uint8_t minvalue = 255;

  for (int i = 0; i < Rows; i++)
    for (int j = 0; j < Cols; j++)
      minvalue = min(minvalue, M[i][j]);

  return minvalue;
}

static uint8_t max(uint8_t a, uint8_t b) { return a > b ? a : b; }

static uint8_t maxChar(int Rows, int Cols, uint8_t M[Rows][Cols]) {
  uint8_t maxvalue = 0;

  for (int i = 0; i < Rows; i++)
    for (int j = 0; j < Cols; j++)
      maxvalue = max(maxvalue, M[i][j]);

  return maxvalue;
}

static MatChar GetBlockChar(MatChar self, int row_from, int row_to,
                            int col_from, int col_to) {
  assert(row_from >= 0);
  assert(col_from >= 0);
  assert(row_from < row_to);
  assert(col_from < col_to);
  assert(row_to <= self.rows);
  assert(col_to <= self.cols);
  assert(self.data);

  MatChar result;
  result.rows = row_to - row_from;
  result.cols = col_to - col_from;
  result.step = self.step;
  result.start = self.start + row_from * self.step + col_from;
  result.data = self.data;

  return result;
}

static MatFloat GetBlockFloat(MatFloat self, int row_from, int row_to,
                              int col_from, int col_to) {
  assert(row_from >= 0);
  assert(col_from >= 0);
  assert(row_from < row_to);
  assert(col_from < col_to);
  assert(row_to <= self.rows);
  assert(col_to <= self.cols);
  assert(self.data);

  MatFloat result;
  result.rows = row_to - row_from;
  result.cols = col_to - col_from;
  result.step = self.step;
  result.start = self.start + row_from * self.step + col_from;
  result.data = self.data;

  return result;
}

static void convertFromCharToFloatArray(int InRows, int InCols,
                                        uint8_t In[InRows][InCols],
                                        float quotient, float shift,
                                        int OutRows, int OutCols,
                                        float Out[OutRows][OutCols]) {
  assert(InRows == OutRows);
  assert(InCols == OutCols);

  for (int i = 0; i < InRows; i++)
    for (int j = 0; j < InCols; j++)
      Out[i][j] = quotient * (float)In[i][j] + shift;
  return;
}

static MatFloat reshapeFloat(MatFloat self, int new_rows) {

  assert(self.cols == self.step);
  assert((self.cols * self.rows) % new_rows == 0);

  MatFloat result;
  result.rows = new_rows;
  result.cols = self.cols * self.rows / new_rows;
  assert(result.cols == 1);
  result.step = result.cols;
  result.start = self.start;
  result.data = self.data;

  return result;
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
  return;
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
  return;
}

// returns alpha*A*B + beta * C
static void gemmFloat(MatFloat A, MatFloat B, float alpha, MatFloat C,
                      float beta, MatFloat *result) {
  assert(A.rows == C.rows);
  assert(A.cols == B.rows);
  assert(B.cols == C.cols);
  assert(C.rows == result->rows);
  assert(C.cols == result->cols);

  float sum = 0;

  for (int q = 0; q < C.rows; q++)
    for (int w = 0; w < C.cols; w++) {
      sum = 0;
      for (int e = 0; e < A.cols; e++) {
        sum += A.data[q * A.step + e + A.start] *
                      B.data[e * B.step + w + B.start];
      }

      result->data[q * result->step + w + result->start] =
          alpha * sum + beta * C.data[q * C.step + w + C.start];
    }

  return;
}

// returns alpha*A*B + beta * C
static void gemmFloatArray(int ARows, int ACols, float A[ARows][ACols],
                           int BRows, int BCols, float B[ARows][BCols],
                           float alpha, int CRows, int CCols,
                           float C[CRows][CCols], float beta, int ResRows,
                           int ResCols, float Res[ResRows][ResCols]) {
  assert(ARows == CRows);
  assert(ACols == BRows);
  assert(BCols == CCols);
  assert(CRows == ResRows);
  assert(CCols == ResCols);

  for (int i = 0; i < CRows; i++)
    for (int j = 0; j < CCols; j++) {
      Res[i][j] = beta * C[i][j];
      for (int k = 0; k < ACols; k++) {
        Res[i][j] += alpha * A[i][k] * B[k][j];
      }
    }

  return;
}

static void expFloat(int rows, int cols, float Mat[rows][cols]) {

  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      Mat[i][j] = expf(Mat[i][j]);

  return;
}

static void addFloat(int rows, int cols, float Mat[rows][cols], float Val) {

  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      Mat[i][j] = Val + Mat[i][j];

  return;
}

static void divideFloat(float Val, int rows, int cols, float Mat[rows][cols]) {

  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      Mat[i][j] = Val / Mat[i][j];

  return;
}

float GetValueFloat(MatFloat self, int row, int col) {
  return self.data[row * self.step + col + self.start];
}

static void SetValueFloat(MatFloat *self, int row, int col, float value) {
  self->data[row * self->step + col + self->start] = value;
  return;
}

static float dotProduct(int LeftRows, int LeftCols, int RightRows,
                        int RightCols, float Left[][LeftCols],
                        float Right[][RightCols]) {
  assert(LeftCols == 1);
  assert(RightCols == 1);
  assert(LeftRows == RightRows);

  float Result = 0.;

  for (int i = 0; i < LeftRows; i++) {
    Result += Left[i][0] * Right[i][0];
  }

  return Result;
}

static void normalizeSample(int imageRows, int imageCols,
                            uint8_t imageArray[imageRows][imageCols],
                            int resultRows, int resultCols,
                            float resultArray[resultRows][resultCols]) {

  float sampleMean = meanChar(imageRows, imageCols, imageArray);
  float sampleMin = minChar(imageRows, imageCols, imageArray);
  float sampleMax = maxChar(imageRows, imageCols, imageArray);

  sampleMax -= sampleMean;
  sampleMin -= sampleMean;

  sampleMax = fmaxf(fabsf(sampleMin), fabsf(sampleMax));

  if (sampleMax == 0.0)
    sampleMax = 1.0;

  convertFromCharToFloatArray(imageRows, imageCols, imageArray, 1.0 / sampleMax,
                              -(1.0 / sampleMax) * sampleMean, resultRows,
                              resultCols, resultArray);
  return;
}

static void copySubArrayChar(int arrayRows, int arrayCols,
                             uint8_t Array[arrayRows][arrayCols],
                             int subArrayRows, int subArrayCols,
                             uint8_t subArray[subArrayRows][subArrayCols],
                             int offsetRow, int offsetCol) {
  for (int i = 0; i < subArrayRows; i++)
    for (int j = 0; j < subArrayCols; j++)
      subArray[i][j] = Array[i + offsetRow][j + offsetCol];
}

/// @brief Calculate a single response map for an image.
///
/// @param Image The image to process.
static void generateResponseMap(
    int ImageRows, int ImageCols, uint8_t Image[ImageRows][ImageCols],
    const Point2i center, int mapSize, mlp classifier,
    float ResponseMap[mapSize + mapSize + 1][mapSize + mapSize + 1]) {

  MatFloat m_U_transpose =
      CreateMatFloat(classifier.m_U.cols, classifier.m_U.rows);
  transposeFloat(classifier.m_U, &m_U_transpose);

  // Subarray
  //
  // This function defines a sub array of a given 2D array. We need to
  // represent this in some way.
  MatFloat wIn_A = GetBlockFloat(classifier.m_wIn, 0, classifier.m_wIn.rows, 0,
                                 classifier.m_wIn.cols - 1);
  MatFloat wIn = CreateMatFloat(wIn_A.rows, m_U_transpose.cols);
  gemmFloat(wIn_A, m_U_transpose, 1.0, wIn, 0.0, &wIn);

  int wInRows = wIn.rows;
  int wInCols = wIn.cols;
  float (*wInArray)[wInCols] = malloc(sizeof(float) * wInRows * wInCols);
  copyMatFloatToArray(wIn, wInRows, wInCols, wInArray);

  // Subarray
  MatFloat bIn =
      GetBlockFloat(classifier.m_wIn, 0, classifier.m_wIn.rows,
                    classifier.m_wIn.cols - 1, classifier.m_wIn.cols);

  // Subarray
  MatFloat wOut_tmp =
      GetBlockFloat(classifier.m_wOut, 0, classifier.m_wOut.rows, 0,
                    classifier.m_wOut.cols - 1);

  MatFloat wOut = CreateMatFloat(wOut_tmp.cols, wOut_tmp.rows);
  transposeFloat(wOut_tmp, &wOut);

  int bInRows = bIn.rows;
  int bInCols = bIn.cols;
  float (*bInArray)[bInCols] = malloc(sizeof(float) * bInRows * bInCols);
  copyMatFloatToArray(bIn, bInRows, bInCols, bInArray);

  int wOutRows = wOut.rows;
  int wOutCols = wOut.cols;
  float (*wOutArray)[wOutCols] = malloc(sizeof(float) * wOutRows * wOutCols);
  copyMatFloatToArray(wOut, wOutRows, wOutCols, wOutArray);

  float bOut = GetValueFloat(classifier.m_wOut, 0, classifier.m_wOut.cols - 1);

  for (int ncy = 0; ncy <= 2 * mapSize; ++ncy) {
    int cy = ncy + center.y - mapSize;
    for (int ncx = 0; ncx <= 2 * mapSize; ++ncx) {
      int cx = ncx + center.x - mapSize;

      // We allocate memory within the program.
      //
      // This is a temporary array, that is just used within this loop.
      // Requering the user to declare this outside cause problems:
      //
      //   a) We bloat the code, as this array needs to be propagated
      //      up to this location.
      //   b) We introduce memory dependences, that would not be here
      //      if we would create different temporary arrays for
      //      different parallel execution streams.
      int xOutRows = bIn.rows;
      int xOutCols = bIn.cols;
      float (*xOutArray)[xOutCols] = malloc(sizeof(float) * xOutRows * xOutCols);

      int imagePatchRows = 2 * classifier.m_patchSize + 1;
      int imagePatchCols = 2 * classifier.m_patchSize + 1;
      uint8_t (*imagePatchArray)[imagePatchCols] =
          malloc(sizeof(uint8_t) * imagePatchRows * imagePatchCols);

      int patchRows = imagePatchRows;
      int patchCols = imagePatchCols;
      float (*patchArray)[patchCols] =
          malloc(sizeof(float) * patchRows * patchCols);

      // Copy a subarray
      //
      // Instead of creating an explicit copy, we could read from the
      // original array and just offset the reads by the offsets the
      // subarray we copy out is set off.
      copySubArrayChar(ImageRows, ImageCols, Image, imagePatchRows,
                       imagePatchCols, imagePatchArray,
                       cy - classifier.m_patchSize,
                       cx - classifier.m_patchSize);

      normalizeSample(imagePatchRows, imagePatchCols, imagePatchArray,
                      patchRows, patchCols, patchArray);

      // Reshape the C99 Array.
      //
      // Is this necessary or can we write the later functions such that they
      // can work on 2D arrays.
      int patchReshapedRows = imagePatchRows * imagePatchCols;
      int patchReshapedCols = 1;
      float (*patchReshapedArray)[patchReshapedCols] = patchArray;

      gemmFloatArray(wInRows, wInCols, wInArray, patchReshapedRows, patchReshapedCols,
                     patchReshapedArray, -1.0, bInRows, bInCols, bInArray, -1.0,
                     xOutRows, xOutCols, xOutArray);

      expFloat(xOutRows, xOutCols, xOutArray);
      addFloat(xOutRows, xOutCols, xOutArray, 1.0f);
      divideFloat(2.0f, xOutRows, xOutCols, xOutArray);
      addFloat(xOutRows, xOutCols, xOutArray, -1.0f);

      ResponseMap[ncy][ncx] =
          (1.0f / (1.0f + expf(-dotProduct(wOutRows, wOutCols, xOutRows,
                                           xOutCols, wOutArray, xOutArray) -
                               bOut)));

      free(patchArray);
      free(imagePatchArray);
      free(xOutArray);
    }
  }

  free(wOutArray);
  freeMatFloat(&wOut);
  freeMatFloat(&wIn);
  freeMatFloat(&m_U_transpose);

  return;
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

  return;
}

/// @brief Calculate a set of response maps.
///
/// This function implements the external interface. In its implementation
/// we transform the input parameters in a way, such that the internal
/// implementation is as close to PENCIL as possible.
void calculateMaps(int NumLandMarks, int MapSize, MatChar Image, MatFloat Shape,
                   mlp Classifiers[], MatFloat *MatResponseMaps[]) {
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
