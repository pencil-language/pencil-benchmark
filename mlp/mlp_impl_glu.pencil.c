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

float GetValueFloat(int N, int M, float A[static const restrict N][M],
		    int row, int col, int offset) PENCIL {
  return A[row][col + offset];
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
                           float In[static const restrict InRows][InCols],
			   int OutRows, int OutCols,
			   float Out[static const restrict OutRows][OutCols])
			   PENCIL {
  assert(InRows == OutCols);
  assert(OutCols == InRows);

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

static void convertFromCharToFloatArray(int imageRows, int imageCols,
                                        uint8_t In[static const restrict imageRows][imageCols],
					int imageOffsetRow, int imageOffsetCol,
                                        float quotient, float shift,
                                        int OutRows, int OutCols,
                                        float Out[static const restrict OutRows][OutCols]) {

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
static void gemmFloatArray_subArray(int ARows, int ACols,
		           float A[static const restrict ARows][ACols],
                           int BRows, int BCols,
			   float B[static const restrict BRows][BCols],
                           float alpha, int CRows, int CCols,
                           float C[static const restrict CRows][CCols],
			   float beta, int ResRows, int ResCols,
			   float Res[static const restrict ResRows][ResCols]) PENCIL {
  assert(BCols == CCols);
  assert(CRows == ResRows);
  assert(CCols == ResCols);

  for (int i = 0; i < CRows; i++)
    for (int j = 0; j < CCols; j++) {
      Res[i][j] = beta * C[i][j];
      for (int k = 0; k < ACols-1; k++) {
        Res[i][j] += alpha * A[i][k] * B[k][j];
      }
    }

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

static void expFloat(int rows, int cols,
		     float Mat[static const restrict rows][cols]) {

  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      Mat[i][j] = expf(Mat[i][j]);

  return;
}

static void addFloat(int rows, int cols,
		     float Mat[static const restrict rows][cols], float Val) {

  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      Mat[i][j] = Val + Mat[i][j];

  return;
}

static void divideFloat(float Val, int rows, int cols,
		        float Mat[static const restrict rows][cols]) {

  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      Mat[i][j] = Val / Mat[i][j];

  return;
}

static float dotProduct(int LeftRows, int LeftCols, int RightRows,
                        int RightCols, float Left[const restrict][LeftCols],
                        float Right[const restrict][RightCols]) {
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

static void copySubArrayFloat(int arrayRows, int arrayCols,
                              float Array[static const restrict arrayRows][arrayCols],
                              int subArrayRows, int subArrayCols,
                              float subArray[static const restrict subArrayRows][subArrayCols],
                              int offsetRow, int offsetCol) PENCIL {
  for (int i = 0; i < subArrayRows; i++)
    for (int j = 0; j < subArrayCols; j++)
      subArray[i][j] = Array[i + offsetRow][j + offsetCol];
}

static float generateResponseMapPatch(
    int mapSize, int ncx, int ncy, int bInRows, int bInCols,
    float bInArray[bInRows][bInCols], mlp classifier, int ImageRows,
    int ImageCols, uint8_t Image[ImageRows][ImageCols], Point2i center,
    int wInRows, int wInCols, float wInArray[wInCols][wInRows], int wOutRows,
    int wOutCols, float wOutArray[wOutRows][wOutCols], float bOut) {
  int cy = ncy + center.y - mapSize;
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

  int imagePatchRows =
      2 * classifier.m_patchSize + 1; // m_patchSize is always 5
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
  // This reshaping linearizes the array. It would be interesting to see if/why
  // this is
  // necessary.
  //
  int patchReshapedRows = imagePatchRows * imagePatchCols; // is always 121
  int patchReshapedCols = 1;

  gemmFloatArray(wInRows, wInCols, wInArray, patchReshapedRows,
                 patchReshapedCols, patchArray, -1.0, bInRows, bInCols,
                 bInArray, -1.0, xOutRows, xOutCols, xOutArray);

  expFloat(xOutRows, xOutCols, xOutArray);
  addFloat(xOutRows, xOutCols, xOutArray, 1.0f);
  divideFloat(2.0f, xOutRows, xOutCols, xOutArray);
  addFloat(xOutRows, xOutCols, xOutArray, -1.0f);

  float result;
  result =
      -dotProduct(wOutRows, wOutCols, xOutRows, xOutCols, wOutArray, xOutArray);

  result -= bOut;
  result = 1.0f / (1.0f + expf(result));

  free(patchArray);
  free(xOutArray);
  return result;
}

// This function has identical behavior as the generateResponseMapPatch
// function, but it does not need to allocate any memory.
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
    float m_U_transposeArray[static const restrict m_UCols][m_URows], //temp array
    int m_wInRows, int m_wInCols,
    float m_wInArray[static const restrict m_wInRows][m_wInCols],
    float wInArray[static const restrict m_wInRows][m_URows], //temp array
    float bInArray[static const restrict m_wInRows][1], //temp array
    int m_wOutRows, int m_wOutCols,
    float m_wOutArray[static const restrict m_wOutRows][m_wOutCols],
    float wOut_tmpArray[static const restrict m_wOutRows][m_wOutCols - 1], //temp array
    float wOutArray[static const restrict m_wOutRows][m_wOutCols - 1] //temp array
    )
    PENCIL {

  int m_U_transposeRows = m_UCols;
  int m_U_transposeCols = m_URows;

  transposeFloat(m_URows, m_UCols, m_UArray, m_U_transposeRows,
                 m_U_transposeCols, m_U_transposeArray);

  int wInRows = m_wInRows;
  int wInCols = m_U_transposeCols;
  
  gemmFloatArray_subArray(m_wInRows, m_wInCols, m_wInArray, m_U_transposeRows,
                 m_U_transposeCols, m_U_transposeArray, 1.0, wInRows, wInCols,
                 wInArray, 0.0, wInRows, wInCols, wInArray);

  // A sub array.
  //
  // Instead of explicitly calculating this, we agreed during the Paris
  // F2F meeting to modify this part of the program so that it works on
  // the original array and to use offsets to determine the subarray,
  int bInRows = m_wInRows;
  int bInCols = 1;
  
  copySubArrayFloat(m_wInRows, m_wInCols, m_wInArray, bInRows,
                    bInCols, bInArray, 0, m_wInCols - 1);


  // A sub array.
  //
  // Instead of explicitly calculating this, we agreed during the Paris
  // F2F meeting to modify this part of the program so that it works on
  // the original array and to use offsets to determine the subarray,
  int wOut_tmpRows = m_wOutRows;
  int wOut_tmpCols = m_wOutCols - 1;
  
  copySubArrayFloat(m_wOutRows, m_wOutCols, m_wOutArray, wOut_tmpRows,
                    wOut_tmpCols, wOut_tmpArray, 0, 0);

  int wOutRows = wOut_tmpCols;
  int wOutCols = wOut_tmpRows;

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

// Calculate a set of respone maps for a given image.
//
// @param Image The image to process
// @param responseMaps A pointer at which results of the calculation,
//                     the response maps, are stored.
void calculateRespondMaps(
    int m_visibleLandmarks_size,
    int MapSize,
    int ImageRows, int ImageCols,
    uint8_t Image[static const restrict ImageRows][ImageCols],
    int shape_rows, int shape_cols,
    float shape_data[static const restrict shape_rows][shape_cols],
    int shape_start,
    mlp m_classifiers[const restrict],
    float ResponseMaps[const restrict][MapSize + MapSize + 1][MapSize + MapSize + 1],
    int max_m_U_transposeRows, int max_m_U_transposeCols,
    float m_U_transposeArray[static const restrict m_visibleLandmarks_size][max_m_U_transposeRows][max_m_U_transposeCols],
    int max_wInRows, int max_wInCols,
    float wInArray[static const restrict m_visibleLandmarks_size][max_wInRows][max_wInCols],
    int max_bInRows, int max_bInCols,
    float bInArray[static const restrict m_visibleLandmarks_size][max_bInRows][max_bInCols],
    int max_wOut_tmpRows, int max_wOut_tmpCols,
    float wOut_tmpArray[static const restrict m_visibleLandmarks_size][max_wOut_tmpRows][max_wOut_tmpCols],
    int max_wOutRows, int max_wOutCols,
    float wOutArray[static const restrict m_visibleLandmarks_size][max_wOutRows][max_wOutCols])
    PENCIL {

#pragma indepdent  
    for (int i = 0; i < m_visibleLandmarks_size; i++)
    {
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
			  m_U_transposeArray[i],
			  m_classifiers[i].m_wIn.rows, m_classifiers[i].m_wIn.cols,	
			  m_classifiers[i].m_wIn.data,
			  wInArray[i],
			  bInArray[i],
			  m_classifiers[i].m_wOut.rows,
			  m_classifiers[i].m_wOut.cols, m_classifiers[i].m_wOut.data,
			  wOut_tmpArray[i], wOutArray[i]);
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

  int max_m_U_transposeRows = 0;
  int max_m_U_transposeCols = 0;
  int max_wInRows = 0;
  int max_wInCols = 0;
  int max_bInRows = 0;
  int max_wOut_tmpRows = 0;
  int max_wOut_tmpCols = 0;
  int max_wOutRows = 0;
  int max_wOutCols = 0;

  // Allocate all the temporary arrays used in the PENCIL code here,
  // but we need first to calculate the maximal possible size of these arrays.
  for (int i = 0; i < NumLandMarks; i++) {
	  int m_UCols = Classifiers[i].m_U.cols;
	  int m_URows = Classifiers[i].m_U.rows;
	  int m_wInCols = Classifiers[i].m_wIn.cols;
	  int m_wInRows = Classifiers[i].m_wIn.cols;
	  int m_wOutRows = Classifiers[i].m_wOut.rows;
	  int m_wOutCols = Classifiers[i].m_wOut.cols;

	  max_m_U_transposeRows = max(m_UCols, max_m_U_transposeRows);
	  max_m_U_transposeCols = max(m_URows, max_m_U_transposeCols);

	  max_wInRows = max(m_wInRows, max_wInRows);
	  max_wInCols = max(max_m_U_transposeCols, max_wInCols);

	  max_bInRows = max(m_wInRows, max_bInRows);

	  max_wOut_tmpRows = max(m_wOutRows, max_wOut_tmpRows);
	  max_wOut_tmpCols = max(m_wOutCols - 1, max_wOut_tmpCols);

	  max_wOutRows = max(max_wOut_tmpCols, max_wOutRows);
	  max_wOutCols = max(max_wOut_tmpRows, max_wOutCols);
  }

  float (* restrict m_U_transposeArray)[max_m_U_transposeRows][max_m_U_transposeCols] =
	     malloc(sizeof(float) * NumLandMarks * max_m_U_transposeRows * max_m_U_transposeCols);
  float (* restrict wInArray)[max_wInRows][max_wInCols] =
	  malloc(sizeof(float) * NumLandMarks * max_wInRows * max_wInCols);
  float (* restrict bInArray)[max_bInRows][1] =
	  malloc(sizeof(float) * NumLandMarks * max_bInRows * 1);
  float (* restrict wOut_tmpArray)[max_wOut_tmpRows][max_wOut_tmpCols] =
	  malloc(sizeof(float) * NumLandMarks * max_wOut_tmpRows * max_wOut_tmpCols);
  float (* restrict wOutArray)[max_wOutRows][max_wOutCols] =
	  malloc(sizeof(float) * NumLandMarks * max_wOutRows * max_wOutCols);

  calculateRespondMaps(NumLandMarks, MapSize, ImageRows, ImageCols, ImageArray,
		  Shape.step, Shape.cols, (float (*) [Shape.cols]) Shape.data,
		  Shape.start, Classifiers, ResponseMaps,
		  max_m_U_transposeRows, max_m_U_transposeCols, m_U_transposeArray,
		  max_wInRows, max_wInCols, wInArray,
		  max_bInRows, 1, bInArray,
		  max_wOut_tmpRows, max_wOut_tmpCols, wOut_tmpArray,
		  max_wOutRows, max_wOutCols, wOutArray);

  for (int i = 0; i < NumLandMarks; ++i)
    copyArrayToMatFloat(Width, Width, ResponseMaps[i],
                        *(&(*MatResponseMaps)[i]));

  free(ImageArray);
  free(ResponseMaps);
}
