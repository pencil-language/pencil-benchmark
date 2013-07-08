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
#include "mlp_impl.h"

static void transposeFloat(int InRows, int InCols,
                           float In[static const restrict InRows][InCols],
			   int OutRows, int OutCols,
			   float Out[static const restrict OutRows][OutCols])
			   PENCIL {
  __pencil_assert(InRows == OutCols);
  __pencil_assert(OutCols == InRows);

  for (int i = 0; i < InRows; i++)
    for (int j = 0; j < InCols; j++)
	Out[j][i] = In[i][j];

  return;
}

static float meanChar(int subImageRows, int subImageCols, int imageRows,
		      int imageCols,
		      uint8_t image[static const restrict imageRows][imageCols],
		      int imageOffsetRow, int imageOffsetCol) PENCIL {
  float sum = 0;

  for (int i = 0; i < subImageRows; i++)
    for (int j = 0; j < subImageCols; j++) {
      sum += image[i + imageOffsetRow][j + imageOffsetCol];
    }

  return sum / (subImageRows * subImageCols);
}

static uint8_t min(uint8_t a, uint8_t b) { return a < b ? a : b; }

static uint8_t minChar(int subImageRows, int subImageCols, int imageRows,
		       int imageCols,
		       uint8_t image[static const restrict imageRows][imageCols],
		       int imageOffsetRow, int imageOffsetCol) PENCIL {
  uint8_t minvalue = 255;

  for (int i = 0; i < subImageRows; i++)
    for (int j = 0; j < subImageCols; j++)
      minvalue = min(minvalue, image[i + imageOffsetRow][j+imageOffsetCol]);

  return minvalue;
}

static uint8_t max(uint8_t a, uint8_t b) { return a > b ? a : b; }

static uint8_t maxChar(int subImageRows, int subImageCols, int imageRows,
		       int imageCols,
		       uint8_t image[static const restrict imageRows][imageCols],
		       int imageOffsetRow, int imageOffsetCol) PENCIL {
  uint8_t maxvalue = 0;

  for (int i = 0; i < subImageRows; i++)
    for (int j = 0; j < subImageCols; j++)
      maxvalue = max(maxvalue, image[i + imageOffsetRow][j+imageOffsetCol]);

  return maxvalue;
}

// returns alpha*A*B + beta * C
static void gemmFloatArray_subArray(int ARows, int ACols,
		           float A[static const restrict ARows][ACols],
                           int BRows, int BCols,
			   float B[static const restrict BRows][BCols],
                           float alpha, int CRows, int CCols,
                           float C[static const restrict CRows][CCols],
			   float beta, int ResRows, int ResCols,
			   float Res[static const restrict ResRows][ResCols]) PENCIL {
  __pencil_assert(BCols == CCols);
  __pencil_assert(CRows == ResRows);
  __pencil_assert(CCols == ResCols);

  for (int i = 0; i < CRows; i++)
    for (int j = 0; j < CCols; j++) {
      Res[i][j] = beta * C[i][j];
      for (int k = 0; k < ACols-1; k++) {
        Res[i][j] += alpha * A[i][k] * B[k][j];
      }
    }

  return;
}

static void copySubArrayFloat(int arrayRows, int arrayCols,
                              float Array[static const restrict arrayRows][arrayCols],
                              int subArrayRows, int subArrayCols,
                              float subArray[static const restrict subArrayRows][subArrayCols],
                              int offsetRow, int offsetCol) PENCIL {
  for (int i = 0; i < subArrayRows; i++)
    for (int j = 0; j < subArrayCols; j++)
      subArray[i][j] = Array[i + offsetRow][j + offsetCol];
}

static float generateResponseMapPatchNoMemory(
    int mapSize, int ncx, int ncy, int bInRows, int bInCols,
    float bInArray[static const restrict bInRows][bInCols],
    mlp classifier, int ImageRows,
    int ImageCols, uint8_t Image[static const restrict ImageRows][ImageCols],
    Point2i center, int wInRows, int wInCols,
    float wInArray[static const restrict wInRows][wInCols], int wOutRows,
    int wOutCols, float wOutArray[static const restrict wOutRows][wOutCols],
    float bOut) PENCIL {
  int cy = ncy + center.y - mapSize;
  int cx = ncx + center.x - mapSize;

  int imagePatchRows = 2 * classifier.m_patchSize + 1;
  int imagePatchCols = 2 * classifier.m_patchSize + 1;

  int patchRows = imagePatchRows;
  int patchCols = imagePatchCols;

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
        xOutArray += alpha * wInArray[i][k] *
                     (quotient * Image[k / imagePatchCols + imageOffsetRow][
                                     k % imagePatchRows + j + imageOffsetCol] +
                      shift);
      }
      xOutArray = expf(xOutArray);
      xOutArray = xOutArray + 1.0f;
      xOutArray = 2.0f / xOutArray;
      xOutArray = xOutArray + -1.0f;
      result += wOutArray[i][0] * xOutArray;
    }
  }

  result = - result;
  result -= bOut;
  result = 1.0f / (1.0f + expf(result));

  return result;
}

/// @brief Calculate a single response map for an image.
///
/// @param Image The image to process.
static void generateResponseMap(
    int ImageRows, int ImageCols,
    uint8_t Image[static const restrict ImageRows][ImageCols],
    const Point2i center, int mapSize, mlp classifier,
    float ResponseMap[static const restrict mapSize + mapSize + 1][mapSize + mapSize + 1],
    int m_URows, int m_UCols,
    float m_UArray[static const restrict m_URows][m_UCols],
    int m_wInRows, int m_wInCols,
    float m_wInArray[static const restrict m_wInRows][m_wInCols],
    int m_wOutRows, int m_wOutCols,
    float m_wOutArray[static const restrict m_wOutRows][m_wOutCols] 
    )
    PENCIL {

  // This is a temporary array.
  //
  // Also, I wonder if this computation can not be precalculated? Or that
  // we could keep track of the index dimensions being transposed and we
  // just access the elements at their previous locations?
  //
  // Also, the size of this array varies / is data-dependent.
  int m_U_transposeRows = m_UCols;
  int m_U_transposeCols = m_URows;
  float m_U_transposeArray[m_U_transposeRows][m_U_transposeCols];

  transposeFloat(m_URows, m_UCols, m_UArray, m_U_transposeRows,
                 m_U_transposeCols, m_U_transposeArray);

  // A temporary array.
  //
  // Why again can this not be precomputed?
  //
  // The size of this array seems to be constant for all test cases.
  int wInRows = m_wInRows;
  int wInCols = m_U_transposeCols;
  float wInArray[wInRows][wInCols];

  gemmFloatArray_subArray(m_wInRows, m_wInCols, m_wInArray, m_U_transposeRows,
                 m_U_transposeCols, m_U_transposeArray, 1.0, wInRows, wInCols,
                 wInArray, 0.0, wInRows, wInCols, wInArray);

  // A sub array.
  //
  // Instead of explicitly calculating this, it would be nice if ppcg
  // could just be told that we only access a subset of this array.
  int bInRows = m_wInRows;
  int bInCols = 1;
  float bInArray[bInRows][bInCols];

  copySubArrayFloat(m_wInRows, m_wInCols, m_wInArray, bInRows,
                    bInCols, bInArray, 0, m_wInCols - 1);


  // A sub array.
  //
  // Instead of explicitly calculating this, it would be nice if ppcg
  // could just be told that we only access a subset of this array.
  int wOut_tmpRows = m_wOutRows;
  int wOut_tmpCols = m_wOutCols - 1;
  float wOut_tmpArray[wOut_tmpRows][wOut_tmpCols];

  copySubArrayFloat(m_wOutRows, m_wOutCols, m_wOutArray, wOut_tmpRows,
                    wOut_tmpCols, wOut_tmpArray, 0, 0);

  int wOutRows = wOut_tmpCols;
  int wOutCols = wOut_tmpRows;
  float wOutArray[wOutRows][wOutCols];

  transposeFloat(wOut_tmpRows, wOut_tmpCols, wOut_tmpArray, wOutRows, wOutCols,
                 wOutArray);

  float bOut = m_wOutArray[0][m_wOutCols - 1];

  for (int ncy = 0; ncy <= 2 * mapSize; ++ncy) {
    for (int ncx = 0; ncx <= 2 * mapSize; ++ncx) {
      ResponseMap[ncy][ncx] = generateResponseMapPatchNoMemory(
          mapSize, ncx, ncy, bInRows, bInCols, bInArray, classifier,
          ImageRows, ImageCols, Image, center, wInRows, wInCols, wInArray,
          wOutRows, wOutCols, wOutArray, bOut);
    }
  }

  return;
}

static int cvRound(float value) {
  return (int)(value + (value >= 0 ? 0.5 : -0.5));
}

float GetValueFloat(int N, int M, float A[static const restrict N][M],
		    int row, int col, int offset) PENCIL {
  return A[row][col + offset];
}

// Calculate a set of respone maps for a given image.
//
// @param Image The image to process
// @param responseMaps A pointer at which results of the calculation,
//                     the response maps, are stored.
void calculateRespondMaps(
    int m_visibleLandmarks_size, int MapSize, int ImageRows, int ImageCols,
    uint8_t Image[static const restrict ImageRows][ImageCols], int shape_rows,
    int shape_cols, float shape_data[static const restrict shape_rows][shape_cols],
    int shape_start, mlp m_classifiers[const restrict],
    float ResponseMaps[const restrict][MapSize + MapSize + 1][MapSize + MapSize + 1]) PENCIL {

// This loop is parallel, ppcg can find parallelism without the need for
// the independent directive.
#pragma indepdent  
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

    shape_x = GetValueFloat(shape_rows, shape_cols, shape_data, 2 * i, 0, shape_start);
    shape_y = GetValueFloat(shape_rows, shape_cols, shape_data, 2 * i + 1, 0, shape_start);
    center.x = cvRound(shape_x);
    center.y = cvRound(shape_y);

    generateResponseMap(ImageRows, ImageCols, Image, center, MapSize,
                        m_classifiers[i], ResponseMaps[i], m_classifiers[i].m_U.rows,
			m_classifiers[i].m_U.cols, m_classifiers[i].m_U.data,
			m_classifiers[i].m_wIn.rows, m_classifiers[i].m_wIn.cols,
			m_classifiers[i].m_wIn.data, m_classifiers[i].m_wOut.rows,
			m_classifiers[i].m_wOut.cols, m_classifiers[i].m_wOut.data);
  }

  return;
}
