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

float GetValueFloat(MatFloat self, int row, int col) {
  return self.data[row * self.step + col + self.start];
}

static void SetValueFloat(MatFloat *self, int row, int col, float value) {
  self->data[row * self->step + col + self->start] = value;
  return;
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

static int cvRound(float value) {
  return (int)(value + (value >= 0 ? 0.5 : -0.5));
}

// Calculate a set of respone maps for a given image.
//
// @param Image The image to process
// @param responseMaps A pointer at which results of the calculation,
//                     the response maps, are stored.
void calculateRespondMaps(
    int m_visibleLandmarks_size, int mapSize, int ImageRows, int ImageCols,
    uint8_t Image[ImageRows][ImageCols], MatFloat shape, mlp m_classifiers[],
    float ResponseMaps[][mapSize + mapSize + 1][mapSize + mapSize + 1]) {

    // Parallel loop
#pragma scop	
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
    int m_URows = m_classifiers[i].m_U.rows;
    int m_UCols = m_classifiers[i].m_U.cols;

    // Reshape the flat array m_classifiers[i].m_U.data into a 2D array
    float (*m_UArray)[m_UCols] = (void*) m_classifiers[i].m_U.data;

    shape_x = GetValueFloat(shape, 2 * i, 0);
    shape_y = GetValueFloat(shape, 2 * i + 1, 0);
    center.x = cvRound(shape_x);
    center.y = cvRound(shape_y);

  // Translate input arrays into C99 Arrays
  assert(m_classifiers[i].m_wIn.start == 0);
  assert(m_classifiers[i].m_wIn.cols == m_classifiers[i].m_wIn.step);
  int m_wInRows = m_classifiers[i].m_wIn.rows;
  int m_wInCols = m_classifiers[i].m_wIn.cols;
  float (*m_wInArray)[m_wInCols] = (void*)m_classifiers[i].m_wIn.data;

  // Translate input arrays into C99 Arrays
  assert(m_classifiers[i].m_wOut.start == 0);
  assert(m_classifiers[i].m_wOut.cols == m_classifiers[i].m_wOut.step);
  int m_wOutRows = m_classifiers[i].m_wOut.rows;
  int m_wOutCols = m_classifiers[i].m_wOut.cols;
  float (*m_wOutArray)[m_wOutCols] = (void*)m_classifiers[i].m_wOut.data;

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

  for (int i = 0; i < m_URows; i++)
    for (int j = 0; j < m_UCols; j++)
	m_U_transposeArray[j][i] = m_UArray[i][j];

  // A temporary array.
  //
  // Why again can this not be precomputed?
  //
  // The size of this array seems to be constant for all test cases.
  int wInRows = m_wInRows;
  int wInCols = m_U_transposeCols;
  float (*wInArray)[wInCols] = malloc(sizeof(float) * wInRows * wInCols);

  // The following gemm uses m_wInArray[][] directly instead of using a sub
  // array of m_wInArray[][], it has been modified accordingly.
  for (int i = 0; i < wInRows; i++)
    for (int j = 0; j < wInCols; j++) {
      wInArray[i][j] = 0;
      for (int k = 0; k < m_wInCols-1; k++) {
        wInArray[i][j] += m_wInArray[i][k] * m_U_transposeArray[k][j];
      }
    }

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

  for (int i = 0; i < wOut_tmpRows; i++)
    for (int j = 0; j < wOut_tmpCols; j++)
	wOutArray[j][i] = wOut_tmpArray[i][j];

  float bOut = m_wOutArray[0][m_wOutCols - 1];

  for (int ncy = 0; ncy <= 2 * mapSize; ++ncy) {
    for (int ncx = 0; ncx <= 2 * mapSize; ++ncx) {

	  int cy = ncy + center.y - mapSize;
	  int cx = ncx + center.x - mapSize;

	  int imagePatchRows =
	      2 * m_classifiers[i].m_patchSize + 1;
	  int imagePatchCols = 2 * m_classifiers[i].m_patchSize + 1;

	  int imageOffsetRow = cy - m_classifiers[i].m_patchSize;
	  int imageOffsetCol = cx - m_classifiers[i].m_patchSize;

	  float sum = 0;
	  for (int l = 0; l < imagePatchRows; l++)
	    for (int j = 0; j < imagePatchCols; j++) {
	      sum += Image[l + imageOffsetRow][j + imageOffsetCol];
	    }
	  float sampleMean = sum / (imagePatchRows * imagePatchCols);

	  uint8_t minvalue = 255;
	  for (int l = 0; l < imagePatchRows; l++)
	    for (int j = 0; j < imagePatchCols; j++)
	      minvalue = min(minvalue, Image[l + imageOffsetRow][j+imageOffsetCol]);
	  float sampleMin = minvalue;

	  uint8_t maxvalue = 0;
	  for (int l = 0; l < imagePatchRows; l++)
	    for (int j = 0; j < imagePatchCols; j++)
	      maxvalue = max(maxvalue, Image[l + imageOffsetRow][j+imageOffsetCol]);
	  float sampleMax = maxvalue;

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

	  for (int l = 0; l < m_wInRows; l++) {
	    for (int j = 0; j < 1; j++) {
					
	      float xOutArray;
	      xOutArray = beta *  m_wInArray[l][j + m_wInCols - 1];
	      for (int k = 0; k < wInCols; k++) {
		xOutArray += alpha * wInArray[l][k] *
			     (quotient * Image[k / imagePatchCols + imageOffsetRow][
					     k % imagePatchRows + j + imageOffsetCol] +
			      shift);
	      }
	      xOutArray = expf(xOutArray);
	      xOutArray = xOutArray + 1.0f;
	      xOutArray = 2.0f / xOutArray;
	      xOutArray = xOutArray + -1.0f;
	      result += wOutArray[l][0] * xOutArray;
	    }
	  }

	  result = - result;
	  result -= bOut;
	  result = 1.0f / (1.0f + expf(result));

          ResponseMaps[i][ncy][ncx] = result; 
    }
  }

  free(m_U_transposeArray);
  free(wInArray);
  free(wOut_tmpArray);
  free(wOutArray);
  }
#pragma endscop

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
