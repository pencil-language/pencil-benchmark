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
  assert(In.start == 0);
  assert(In.rows == n);
  assert(In.cols == m);
  assert(In.step == m);

  for (int i = 0; i < In.rows; i++)
    for (int j = 0; j < In.cols; j++)
      Out[i][j] = In.data[i * In.step + j];
}

static void copyArrayToMatFloat(int n, int m, float In[][m], MatFloat Out) {
  assert(Out.start == 0);
  assert(Out.rows == n);
  assert(Out.cols == m);
  assert(Out.step == m);

  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      Out.data[i * Out.step + j] = In[i][j];
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

static float meanChar(MatChar input) {
  assert(input.data);
  int q, w;
  float sum = 0;
  float c = 0; // kahan summation

  for (q = 0; q < input.rows; q++)
    for (w = 0; w < input.cols; w++) {
      float y = input.data[q * input.step + w + input.start] - c;
      float t = sum + y;
      c = (t - sum) - y;
      sum = t;
    }

  return sum / (input.rows * input.cols);
}

static uint8_t min(uint8_t a, uint8_t b) { return a < b ? a : b; }

static uint8_t minChar(MatChar input) {
  assert(input.data);
  int q, w;
  uint8_t minvalue = 255;

  for (q = 0; q < input.rows; q++)
    for (w = 0; w < input.cols; w++)
      minvalue = min(minvalue, input.data[q * input.step + w + input.start]);

  return minvalue;
}

static uint8_t max(uint8_t a, uint8_t b) { return a > b ? a : b; }

static uint8_t maxChar(MatChar input) {
  assert(input.data);
  int q, w;
  uint8_t maxvalue = 0;

  for (q = 0; q < input.rows; q++)
    for (w = 0; w < input.cols; w++)
      maxvalue = max(maxvalue, input.data[q * input.step + w + input.start]);

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

static void convertFromCharToFloat(MatChar from, float quotient, float shift,
                                   MatFloat *to) {
  assert(from.rows == to->rows);
  assert(from.cols == to->cols);

  int q, w;
  for (q = 0; q < from.rows; q++)
    for (w = 0; w < from.cols; w++)
      to->data[q * to->step + w + to->start] =
          quotient * (float)from.data[q * from.step + w + from.start] + shift;
  return;
}

static MatFloat reshapeFloat(MatFloat self, int new_rows) {

  assert(self.cols == self.step);
  assert((self.cols * self.rows) % new_rows == 0);

  MatFloat result;
  result.rows = new_rows;
  result.cols = self.cols * self.rows / new_rows;
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

  int q, w, e;
  float sum = 0;
  float c;

  for (q = 0; q < C.rows; q++)
    for (w = 0; w < C.cols; w++) {
      sum = 0;
      for (e = 0; e < A.cols; e++) {
        float y = A.data[q * A.step + e + A.start] *
                      B.data[e * B.step + w + B.start] -
                  c;
        float t = sum + y;
        c = (t - sum) - y;
        sum = t;
      }

      result->data[q * result->step + w + result->start] =
          alpha * sum + beta * C.data[q * C.step + w + C.start];
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

static float dotProduct(MatFloat A, MatFloat B) {
  assert(A.cols == 1);
  assert(B.cols == 1);
  assert(A.rows == B.rows);

  float result = 0.;
  float c = 0.;

  int q;
  for (q = 0; q < A.rows; q++) {
    float y = A.data[q * A.step + A.start] * B.data[q * B.step + B.start] - c;
    float t = result + y;
    c = (t - result) - y;
    result = t;
  }

  return result;
}

static void normalizeSample(MatChar image, MatFloat *result) {
  assert(result->cols == image.cols);
  assert(result->rows == image.rows);

  float sampleMean = meanChar(image);
  float sampleMin = minChar(image);
  float sampleMax = maxChar(image);

  sampleMax -= sampleMean;
  sampleMin -= sampleMean;

  sampleMax = fmaxf(fabsf(sampleMin), fabsf(sampleMax));

  if (sampleMax == 0.0)
    sampleMax = 1.0;

  convertFromCharToFloat(image, 1.0 / sampleMax,
                         -(1.0 / sampleMax) * sampleMean, result);

  *result = reshapeFloat(*result, image.rows * image.cols);

  return;
}

/// @brief Calculate a single response map for an image.
///
/// @param Image The image to process.
static void generateResponseMap(
    const MatChar Image, const Point2i center, int mapSize, mlp classifier,
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
      int bRows = bIn.rows;
      int bCols = bIn.cols;
      float (*OutArray)[bCols] = malloc(sizeof(float) * bRows * bCols);

      MatChar imagePatch = GetBlockChar(
          Image, cy - classifier.m_patchSize, cy + classifier.m_patchSize + 1,
          cx - classifier.m_patchSize, cx + classifier.m_patchSize + 1);
      MatFloat patch = CreateMatFloat(imagePatch.rows, imagePatch.cols);

      // Array reshaping
      //
      // This function is changing the interpretation of the array from a
      // 2D array to a 1D array, such that the dot-product can be applied.
      // We need to figure out how to interpret this. Either we support
      // array reshaping or we try to rewrite the code to not include array
      // reshaping. One solution to get around the need for linearization is to
      // define some kind of 2D dot product.
      normalizeSample(imagePatch, &patch);

      MatFloat xOut = CreateMatFloat(bIn.rows, bIn.cols);

      gemmFloat(wIn, patch, -1.0, bIn, -1.0, &xOut);

      copyMatFloatToArray(xOut, bRows, bCols, OutArray);

      expFloat(bRows, bCols, OutArray);
      addFloat(bRows, bCols, OutArray, 1.0f);
      divideFloat(2.0f, bRows, bCols, OutArray);
      addFloat(bRows, bCols, OutArray, -1.0f);

      copyArrayToMatFloat(bRows, bCols, OutArray, xOut);

      ResponseMap[ncy][ncx] =
          (1.0f / (1.0f + expf(-dotProduct(wOut, xOut) - bOut)));

      free(OutArray);
      freeMatFloat(&xOut);
      freeMatFloat(&patch);
    }
  }

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
    int m_visibleLandmarks_size, int MapSize, MatChar Image, MatFloat shape,
    mlp m_classifiers[],
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

    generateResponseMap(Image, center, MapSize, m_classifiers[i],
                        ResponseMaps[i]);
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

  calculateRespondMaps(NumLandMarks, MapSize, Image, Shape, Classifiers,
                       ResponseMaps);

  for (int i = 0; i < NumLandMarks; ++i)
    copyArrayToMatFloat(Width, Width, ResponseMaps[i],
                        *(&(*MatResponseMaps)[i]));

  free(ResponseMaps);
}
