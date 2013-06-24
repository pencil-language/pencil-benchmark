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

static void transposeFloat(int InRows, int InCols,
                                float In[InRows][InCols], int OutRows,
                                int OutCols, float Out[OutRows][OutCols]) {
  assert(InRows == OutCols);
  assert(OutCols == InRows);

  for (int i = 0; i < InRows; i++)
    for (int j = 0; j < InCols; j++)
	Out[j][i] = In[i][j];

  return;
}

static float meanChar(int subImageRows, int subImageCols, int imageRows,
		      int imageCols, uint8_t image[imageRows][imageCols],
		      int imageOffsetRow, int imageOffsetCol) {
  float sum = 0;

  for (int i = 0; i < subImageRows; i++)
    for (int j = 0; j < subImageCols; j++) {
      sum += image[i + imageOffsetRow][j + imageOffsetCol];
    }

  return sum / (subImageRows * subImageCols);
}

static uint8_t min(uint8_t a, uint8_t b) { return a < b ? a : b; }

static uint8_t minChar(int subImageRows, int subImageCols, int imageRows,
		       int imageCols, uint8_t image[imageRows][imageCols],
		       int imageOffsetRow, int imageOffsetCol) {
  uint8_t minvalue = 255;

  for (int i = 0; i < subImageRows; i++)
    for (int j = 0; j < subImageCols; j++)
      minvalue = min(minvalue, image[i + imageOffsetRow][j+imageOffsetCol]);

  return minvalue;
}

static uint8_t max(uint8_t a, uint8_t b) { return a > b ? a : b; }

static uint8_t maxChar(int subImageRows, int subImageCols, int imageRows,
		       int imageCols, uint8_t image[imageRows][imageCols],
		       int imageOffsetRow, int imageOffsetCol) {
  uint8_t maxvalue = 0;

  for (int i = 0; i < subImageRows; i++)
    for (int j = 0; j < subImageCols; j++)
      maxvalue = max(maxvalue, image[i + imageOffsetRow][j+imageOffsetCol]);

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

static void convertFromCharToFloatArray(int imageRows, int imageCols,
                                        uint8_t In[imageRows][imageCols],
					int imageOffsetRow, int imageOffsetCol,
                                        float quotient, float shift,
                                        int OutRows, int OutCols,
                                        float Out[OutRows][OutCols]) {

  for (int i = 0; i < OutRows; i++)
    for (int j = 0; j < OutCols; j++)
      Out[i][j] = quotient * (float)In[i + imageOffsetRow][j + imageOffsetCol] + shift;
  return;
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

/* This function takes an image[imageRows][imageCols] as input and then works on a subset of that
 * image. imageOffsetRow and imageOffsetCol are used as offsets for image[][] while
 * subImageRows ad subImageCols define the size of the subimage.
 */
static void normalizeSample(int subImageRows, int subImageCols, int imageRows,
		            int imageCols,
			    uint8_t imageArray[imageRows][imageCols],
			    int imageOffsetRow, int imageOffsetCol,
                            int resultRows, int resultCols,
                            float resultArray[resultRows][resultCols]) {

  float sampleMean = meanChar(subImageRows, subImageCols, imageRows,
		              imageCols, imageArray, imageOffsetRow,
                              imageOffsetCol);
  float sampleMin = minChar(subImageRows, subImageCols, imageRows,
		            imageCols, imageArray, imageOffsetRow,
                            imageOffsetCol);
  float sampleMax = maxChar(subImageRows, subImageCols, imageRows,
		            imageCols, imageArray, imageOffsetRow,
                            imageOffsetCol);

  sampleMax -= sampleMean;
  sampleMin -= sampleMean;

  sampleMax = fmaxf(fabsf(sampleMin), fabsf(sampleMax));

  if (sampleMax == 0.0)
    sampleMax = 1.0;

  convertFromCharToFloatArray(imageRows, imageCols, imageArray, imageOffsetRow,
                              imageOffsetCol, 1.0 / sampleMax,
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

static void copySubArrayFloat(int arrayRows, int arrayCols,
                              float Array[arrayRows][arrayCols],
                              int subArrayRows, int subArrayCols,
                              float subArray[subArrayRows][subArrayCols],
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

  // Translate input arrays into C99 Arrays.
  //
  // The problem here is that the size of the arrays is not know to
  // be identical for each respond map calculation. Hence, we can
  // not easily move those allocations out of the core computation.
  int m_URows = classifier.m_U.rows; // Always 121
  int m_UCols = classifier.m_U.cols; // Varies between 17 and 30
  float (*m_UArray)[m_UCols] =
  malloc(sizeof(float) * m_URows * m_UCols);
  copyMatFloatToArray(classifier.m_U, m_URows, m_UCols, m_UArray);

  // Translate input arrays into C99 Arrays
  int m_wInRows = classifier.m_wIn.rows; // Always 25
  int m_wInCols = classifier.m_wIn.cols; // Varies between 17 and 31
  float (*m_wInArray)[m_wInCols] =
  malloc(sizeof(float) * m_wInRows * m_wInCols);
  copyMatFloatToArray(classifier.m_wIn, m_wInRows, m_wInCols, m_wInArray);

  // Translate input arrays into C99 Arrays
  int m_wOutRows = classifier.m_wOut.rows; // Always 1
  int m_wOutCols = classifier.m_wOut.cols; // Always 26
  float (*m_wOutArray)[m_wOutCols] =
      malloc(sizeof(float) * m_wOutRows * m_wOutCols);
  copyMatFloatToArray(classifier.m_wOut, m_wOutRows, m_wOutCols, m_wOutArray);

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
      int xOutRows = bInRows;
      int xOutCols = bInCols;
      float (*xOutArray)[xOutCols] = malloc(sizeof(float) * xOutRows * xOutCols);

      int imagePatchRows = 2 * classifier.m_patchSize + 1; // m_patchSize is always 5
      int imagePatchCols = 2 * classifier.m_patchSize + 1;

      int patchRows = imagePatchRows;
      int patchCols = imagePatchCols;
      float (*patchArray)[patchCols] =
          malloc(sizeof(float) * patchRows * patchCols);

      normalizeSample(imagePatchRows, imagePatchCols, ImageRows, ImageCols, Image,
		      cy - classifier.m_patchSize, cx - classifier.m_patchSize,
                      patchRows, patchCols, patchArray);

      // Reshape the C99 Array.
      //
      // This reshaping linearizes the array. It would be interesting to see if/why this is
      // necessary.
      //
      int patchReshapedRows = imagePatchRows * imagePatchCols; // is always 121
      int patchReshapedCols = 1;

      gemmFloatArray(wInRows, wInCols, wInArray, patchReshapedRows, patchReshapedCols,
                     patchArray, -1.0, bInRows, bInCols, bInArray, -1.0,
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
      free(xOutArray);
    }
  }

  free(m_wInArray);
  free(wIn_AArray);
  free(m_U_transposeArray);
  free(m_UArray);
  free(wInArray);
  free(m_wOutArray);
  free(wOut_tmpArray);
  free(wOutArray);

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
