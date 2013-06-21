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

static void expFloat(MatFloat input, MatFloat *output) {
  assert(input.rows == output->rows);
  assert(input.cols == output->cols);

  int q, w;
  for (q = 0; q < input.rows; q++)
    for (w = 0; w < input.cols; w++)
      output->data[q * output->step + w + output->start] =
          exp(input.data[q * input.step + w + input.start]);

  return;
}

static void addFloat(MatFloat input, float val, MatFloat *output) {
  assert(input.rows == output->rows);
  assert(input.cols == output->cols);

  int q, w;
  for (q = 0; q < input.rows; q++)
    for (w = 0; w < input.cols; w++)
      output->data[q * output->step + w + output->start] =
          val + input.data[q * input.step + w + input.start];

  return;
}

static void divideFloat(float val, MatFloat input, MatFloat *output) {
  assert(input.rows == output->rows);
  assert(input.cols == output->cols);

  int q, w;
  for (q = 0; q < input.rows; q++)
    for (w = 0; w < input.cols; w++)
      output->data[q * output->step + w + output->start] =
          val / input.data[q * input.step + w + input.start];

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

/*  convertFromCharToFloat(image, 1.0 / sampleMax,
                         -(1.0 / sampleMax) * sampleMean, result);*/
  int q, w;
  for (q = 0; q < image.rows; q++)
    for (w = 0; w < image.cols; w++)
      result->data[q * result->step + w + result->start] =
          (1.0 / sampleMax) * (float)image.data[q * image.step + w + image.start] - (1.0 / sampleMax) * sampleMean;

//  assert(result.cols == result.step);
//  assert((result.cols * result.rows) % new_rows == 0);

  result->cols = result->cols * result->rows / (image.rows * image.cols);
  result->rows = image.rows * image.cols;
  result->step = result->cols;


//  *result = reshapeFloat(*result, image.rows * image.cols);

  return;
}

static void generateResponseMap(const MatChar image, const Point2i center,
                                int mapSize, mlp classifier, MatFloat *result) {

  int q, w, e;
  float sum = 0;
  float c;

  assert(result->rows == 2 * mapSize + 1);
  assert(result->cols == 2 * mapSize + 1);

  MatFloat m_U_transpose =
      CreateMatFloat(classifier.m_U.cols, classifier.m_U.rows);
 


  /* transposeFloat(classifier.m_U, &m_U_transpose); */
  assert(m_U_transpose.rows == classifier.m_U.cols);
  assert(m_U_transpose.cols == classifier.m_U.rows);
  for (q = 0; q < classifier.m_U.rows; q++)
    for (w = 0; w < classifier.m_U.cols; w++)
      m_U_transpose.data[w * m_U_transpose.step + q + m_U_transpose.start] =
          classifier.m_U.data[q * classifier.m_U.step + w + classifier.m_U.start];
  /* --------------------- END ------------------- */



  MatFloat wIn_A = GetBlockFloat(classifier.m_wIn, 0, classifier.m_wIn.rows, 0,
                                 classifier.m_wIn.cols - 1);
  MatFloat wIn = CreateMatFloat(wIn_A.rows, m_U_transpose.cols);



  /* gemmFloat(wIn_A, m_U_transpose, 1.0, wIn, 0.0, &wIn); */
  for (q = 0; q < wIn.rows; q++)
    for (w = 0; w < wIn.cols; w++) {
      sum = 0;
      for (e = 0; e < wIn_A.cols; e++) {
        float y = wIn_A.data[q * wIn_A.step + e + wIn_A.start] *
                      m_U_transpose.data[e * m_U_transpose.step + w + m_U_transpose.start] - c;
        float t = sum + y;
        c = (t - sum) - y;
        sum = t;
      }

      wIn.data[q * wIn.step + w + wIn.start] = sum;
    }
  /* --------------------- END ------------------- */



  MatFloat bIn =
      GetBlockFloat(classifier.m_wIn, 0, classifier.m_wIn.rows,
                    classifier.m_wIn.cols - 1, classifier.m_wIn.cols);

  MatFloat wOut_tmp =
      GetBlockFloat(classifier.m_wOut, 0, classifier.m_wOut.rows, 0,
                    classifier.m_wOut.cols - 1);
  MatFloat wOut = CreateMatFloat(wOut_tmp.cols, wOut_tmp.rows);



  /* transposeFloat(wOut_tmp, &wOut); */
  assert(wOut.rows == wOut_tmp.cols);
  assert(wOut.cols == wOut_tmp.rows);
  for (q = 0; q < wOut_tmp.rows; q++)
    for (w = 0; w < wOut_tmp.cols; w++)
      wOut.data[w * wOut.step + q + wOut.start] =
          wOut_tmp.data[q * wOut_tmp.step + w + wOut_tmp.start];
  /* --------------------- END ------------------- */



  float bOut = GetValueFloat(classifier.m_wOut, 0, classifier.m_wOut.cols - 1);

  int ncy = 0;
  int cy = 0;
  int ncx = 0;
  int cx = 0;

  for (ncy = 0, cy = center.y - mapSize; cy <= center.y + mapSize;
       ++ncy, ++cy) {
    for (ncx = 0, cx = center.x - mapSize; cx <= center.x + mapSize;
         ++ncx, ++cx) {

      MatChar imagePatch = GetBlockChar(
          image, cy - classifier.m_patchSize, cy + classifier.m_patchSize + 1,
          cx - classifier.m_patchSize, cx + classifier.m_patchSize + 1);
      MatFloat patch = CreateMatFloat(imagePatch.rows, imagePatch.cols);



      /* normalizeSample(imagePatch, &patch); */
      float sampleMean = meanChar(imagePatch);
      float sampleMin = minChar(imagePatch);
      float sampleMax = maxChar(imagePatch);

      sampleMax -= sampleMean;
      sampleMin -= sampleMean;

      sampleMax = fmaxf(fabsf(sampleMin), fabsf(sampleMax));

      if (sampleMax == 0.0)
        sampleMax = 1.0;

      int q, w;
      for (q = 0; q < imagePatch.rows; q++)
        for (w = 0; w < imagePatch.cols; w++)
          patch.data[q * patch.step + w + patch.start] =
             (1.0 / sampleMax) * (float)imagePatch.data[q * imagePatch.step + w + imagePatch.start] - (1.0 / sampleMax) * sampleMean;

      patch.cols = patch.cols * patch.rows / (imagePatch.rows * imagePatch.cols);
      patch.rows = imagePatch.rows * imagePatch.cols;
      patch.step = patch.cols;
      /* --------------------- END ------------------- */



      MatFloat xOut = CreateMatFloat(bIn.rows, bIn.cols);



      /* gemmFloat(wIn, patch, -1.0, bIn, -1.0, &xOut);*/
      assert(wIn.rows == bIn.rows);
      assert(wIn.cols == patch.rows);
      assert(patch.cols == bIn.cols);
      assert(bIn.rows == xOut.rows);
      assert(bIn.cols == xOut.cols);

      for (q = 0; q < bIn.rows; q++)
        for (w = 0; w < bIn.cols; w++) {
          sum = 0;
          for (e = 0; e < wIn.cols; e++) {
            float y = wIn.data[q * wIn.step + e + wIn.start] *
                      patch.data[e * patch.step + w + patch.start] - c;
            float t = sum + y;
            c = (t - sum) - y;
            sum = t;
          }

        xOut.data[q * xOut.step + w + xOut.start] =
            (-1.0) * sum + (-1.0) * bIn.data[q * bIn.step + w + bIn.start];
      }
      /* --------------------- END ------------------- */



      MatFloat e = CreateMatFloat(xOut.rows, xOut.cols);

      assert(xOut.rows == e.rows);
      assert(xOut.cols == e.cols);

      

      /* expFloat(xOut, &e); */
      for (q = 0; q < xOut.rows; q++)
        for (w = 0; w < xOut.cols; w++)
          e.data[q * e.step + w + e.start] =
            exp(xOut.data[q * xOut.step + w + xOut.start]);
      /* --------------------- END ------------------- */



      assert(xOut.rows == e.rows);
      assert(xOut.cols == e.cols);



      /* addFloat(e, 1.0, &xOut); */
      for (q = 0; q < xOut.rows; q++)
        for (w = 0; w < xOut.cols; w++)
          xOut.data[q * xOut.step + w + xOut.start] = 1.0 +
            e.data[q * e.step + w + e.start];      
      /* --------------------- END ------------------- */



      assert(xOut.rows == e.rows);
      assert(xOut.cols == e.cols);



      /* divideFloat(2.0, xOut, &e); */
      for (q = 0; q < xOut.rows; q++)
        for (w = 0; w < xOut.cols; w++)
          e.data[q * e.step + w + e.start] = 2.0 /
            xOut.data[q * xOut.step + w + xOut.start];
      /* --------------------- END ------------------- */



      assert(xOut.rows == e.rows);
      assert(xOut.cols == e.cols);



      /* addFloat(e, -1.0, &xOut); */
      for (q = 0; q < xOut.rows; q++)
        for (w = 0; w < xOut.cols; w++)
          xOut.data[q * xOut.step + w + xOut.start] = -1.0 +
            e.data[q * e.step + w + e.start];
      /* --------------------- END ------------------- */



      result->data[ncy * (result->step) + ncx + (result->start)] = (1. / (1. + expf(-dotProduct(wOut, xOut)) - bOut));
      result->data[ncy * (result->step) + ncx + (result->start)] = (1. / (1. + expf(-dotProduct(wOut, xOut) - bOut)));

      freeMatFloat(&e);
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

void calculateMaps(int m_visibleLandmarks_size, int m_mapSize,
                   MatChar alignedImage, MatFloat shape, mlp m_classifiers[],
                   // results
                   MatFloat *responseMaps[]) {
  int q;
  for (q = 0; q < m_visibleLandmarks_size; q++) {
    int idx = q;

    Point2i center;
    float shape_x;
    float shape_y;

    shape_x = GetValueFloat(shape, 2 * idx, 0);
    shape_y = GetValueFloat(shape, 2 * idx + 1, 0);
    center.x = cvRound(shape_x);
    center.y = cvRound(shape_y);

    generateResponseMap(alignedImage, center, m_mapSize, m_classifiers[idx],
                        (&(*responseMaps)[q]));
  }

  return;
}
