#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <sys/types.h>
#include <SDL2/SDL_render.h>
#include <forward_list>
#include <fftw3.h>
#include <x86intrin.h>
#define SIZE 1024
#define SEARCH_RADIUS 64
#define AP_DURATION 8
#define REST_DURATION 4
#define AP_THRESHOLD 16

enum CellType {
  // A heart cell here is represented either as a pacemaker cell, or a normal tissue cell
  Pacemaker,
  Tissue,
  RestingTissue
};

// A cell is declared as a struct since it has no methods
struct Cell {
  CellType type;
  // The state is a positive integer
  uint state;
  uint orientationIndex;
};

struct Orientation {
  float xDir;
  float yDir;
  uint cellCount;
  std::forward_list<Cell*> cells;
};

struct Cells {
  uint width;
  uint height;
  // 2D array is represented as one contiguous block of memory, for performance reasons
  // (indexed [i][j] would be [i * height + j])
  Cell* cells;
  uint numOrientations;
  Orientation* orientations;
};

class NeighbourCounter {
  public:
    Cells* cells;
    double* stateArray;
    double** distanceCoefficients;
    double** distanceCoefficientsPadded;
    fftw_complex** distanceCoefficientsTransformed;
    fftw_plan* distanceCoefficientsFFT;
    fftw_complex** neighbourArraysTransformed;
    fftw_plan* stateArrayFFT;
    fftw_plan* stateArrayIFFT;
    // 2D array here
    double* neighbourArray;
    double** neighbourArrays;
    uint numOrientations;
    NeighbourCounter(Cells* cells, double* stateArray) {
      this->stateArray = stateArray;
      this->cells = cells;
      distanceCoefficients = new double*[cells->numOrientations];
      distanceCoefficientsPadded = new double*[cells->numOrientations];
      distanceCoefficientsTransformed = new fftw_complex*[cells->numOrientations];
      neighbourArraysTransformed = new fftw_complex*[cells->numOrientations];
      neighbourArrays = new double*[cells->numOrientations];
      distanceCoefficientsFFT = new fftw_plan[cells->numOrientations];
      stateArrayFFT = new fftw_plan[cells->numOrientations];
      stateArrayIFFT = new fftw_plan[cells->numOrientations];
      for (int i = 0; i < cells->numOrientations; i++) {
        distanceCoefficients[i] = fftw_alloc_real(cells->height * cells->width);
        distanceCoefficientsPadded[i] = fftw_alloc_real(cells->height * cells->width);
        neighbourArrays[i] = fftw_alloc_real(cells->height * cells->width);
        neighbourArraysTransformed[i] = fftw_alloc_complex(cells->height * (cells->width / 2 + 1));
        distanceCoefficientsTransformed[i] = fftw_alloc_complex(cells->height * (cells->width / 2 + 1));
        distanceCoefficientsFFT[i] = fftw_plan_dft_r2c_2d(cells->height, cells->width, distanceCoefficientsPadded[i], distanceCoefficientsTransformed[i], 0);
        stateArrayFFT[i] = fftw_plan_dft_r2c_2d(cells->height, cells->width, stateArray, neighbourArraysTransformed[i], 0);
        stateArrayIFFT[i] = fftw_plan_dft_c2r_2d(cells->height, cells->width, neighbourArraysTransformed[i], neighbourArrays[i], 0);
      }
      numOrientations = cells->numOrientations;
      neighbourArray = fftw_alloc_real(cells->height * cells->width * numOrientations);
      initialize();
    }
    ~NeighbourCounter() {
      delete[] stateArrayFFT;
      delete[] distanceCoefficientsFFT;
      delete[] stateArrayIFFT;
      for (int i = 0; i < cells->numOrientations; i++) {
        fftw_free(distanceCoefficients[i]);
        fftw_free(distanceCoefficientsPadded[i]);
        fftw_free(distanceCoefficientsTransformed[i]);
        fftw_free(neighbourArraysTransformed[i]);
        fftw_free(neighbourArrays[i]);
      }
      delete[] distanceCoefficients;
      delete[] distanceCoefficientsPadded;
      delete[] neighbourArrays;
      delete[] neighbourArraysTransformed;
      delete[] distanceCoefficientsTransformed;
    }
    void reinitialize() {
      // If the number of orientations has changed, then all the arrays must be reinitialized
      if (cells->numOrientations != numOrientations) {
        numOrientations = cells->numOrientations;
        delete[] distanceCoefficientsFFT;
        delete[] stateArrayFFT;
        delete[] stateArrayIFFT;
        for (int i = 0; i < numOrientations; i++) {
          fftw_free(distanceCoefficients[i]);
          fftw_free(distanceCoefficientsPadded[i]);
          fftw_free(distanceCoefficientsTransformed[i]);
          fftw_free(neighbourArraysTransformed[i]);
          fftw_free(neighbourArrays[i]);
        }
        delete[] distanceCoefficients;
        delete[] distanceCoefficientsPadded;
        delete[] distanceCoefficientsTransformed;
        delete[] neighbourArraysTransformed;
        delete[] neighbourArrays;
        distanceCoefficients = new double*[cells->numOrientations];
        distanceCoefficientsPadded = new double*[cells->numOrientations];
        distanceCoefficientsTransformed = new fftw_complex*[cells->numOrientations];
        neighbourArraysTransformed = new fftw_complex*[cells->numOrientations];
        neighbourArrays = new double*[cells->numOrientations];
        distanceCoefficientsFFT = new fftw_plan[cells->numOrientations];
        stateArrayFFT = new fftw_plan[cells->numOrientations];
        stateArrayIFFT = new fftw_plan[cells->numOrientations];
        for (int i = 0; i < cells->numOrientations; i++) {
          distanceCoefficients[i] = fftw_alloc_real(cells->height * cells->width);
          distanceCoefficientsPadded[i] = fftw_alloc_real(cells->height * cells->width);
          neighbourArrays[i] = fftw_alloc_real(cells->height * cells->width);
          neighbourArraysTransformed[i] = fftw_alloc_complex(cells->height * (cells->width / 2 + 1));
          distanceCoefficientsTransformed[i] = fftw_alloc_complex(cells->height * (cells->width / 2 + 1));
          distanceCoefficientsFFT[i] = fftw_plan_dft_r2c_2d(cells->height, cells->width, distanceCoefficientsPadded[i], distanceCoefficientsTransformed[i], 0);
          stateArrayFFT[i] = fftw_plan_dft_r2c_2d(cells->height, cells->width, stateArray, neighbourArraysTransformed[i], 0);
          stateArrayIFFT[i] = fftw_plan_dft_c2r_2d(cells->height, cells->width, neighbourArraysTransformed[i], neighbourArrays[i], 0);
        }
      }
      initialize();
    }

    void calculateNeighbourCounts() {
      // TODO: switch to standard convolution if the number of cells is low enough
      for (int i = 0; i < numOrientations; i++) {
        fftw_execute(stateArrayFFT[i]);
        multiply(neighbourArraysTransformed[i], distanceCoefficientsTransformed[i]);
        fftw_execute(stateArrayIFFT[i]);
      }
      for (int i = 0; i < cells->height * cells->width; i++) {
        for (int j = 0; j < numOrientations; j++) {
          neighbourArray[(i * numOrientations) + j] = neighbourArrays[j][i];
        }
      }
    }
  private:
    void calculateDistanceCoefficients(Orientation orientation, double* coefficients) {
      for (int i = 0; i < SEARCH_RADIUS; i++) {
        for (int j = 0; j < SEARCH_RADIUS; j++) {
          if (i == SEARCH_RADIUS / 2 && j == SEARCH_RADIUS / 2) continue;
          double xCoord = (j - SEARCH_RADIUS / 2.0);
          double yCoord = (i - SEARCH_RADIUS / 2.0);
          double distance = xCoord * xCoord + yCoord * yCoord;
          double dotProduct = xCoord * orientation.xDir + yCoord * orientation.yDir;
          double cosTheta = dotProduct / (std::sqrt(distance) * std::sqrt(orientation.xDir * orientation.xDir + orientation.yDir * orientation.yDir));
          coefficients[(i * SEARCH_RADIUS) + j] = (1.0 / distance) * (cosTheta + 1) * 0.5;
        }
      }
    }
    // Cyclically shifts a convolution to expand it, i.e. for padding reasons
    void shiftConvolution(double* originalConvolution, double* shiftedConvolution, int convWidth, int dataHeight, int dataLength) {
      //Top left corner (shifted so that the middle element of the convolution is now at (0, 0))
      for (int i = convWidth / 2; i < convWidth; i++) {
        for (int j = convWidth / 2; j < convWidth; j++) {
          shiftedConvolution[(i - (convWidth / 2)) * dataLength + (j - (convWidth / 2))] = originalConvolution[i * convWidth + j];
        }
      }
      // Top right corner
      for (int i = convWidth / 2; i < convWidth; i++) {
        for (int j = 0; j < convWidth / 2; j++) {
          shiftedConvolution[(i - (convWidth / 2)) * dataLength + (dataLength + j - convWidth / 2)] = originalConvolution[i * convWidth + j];
        }
      }
      // Bottom left corner
      for (int i = 0; i < convWidth / 2; i++) {
        for (int j = convWidth / 2; j < convWidth; j++) {
          shiftedConvolution[(dataHeight + i - convWidth / 2) * dataLength + (j - convWidth / 2)] = originalConvolution[i * convWidth + j];
        }
      }
      // Bottom right corner
      for (int i = 0; i < (convWidth / 2); i++) {
        for (int j = 0; j < (convWidth / 2); j++) {
          shiftedConvolution[(dataHeight + i - convWidth / 2) * dataLength + (dataLength + j - convWidth / 2)] = originalConvolution[i * convWidth + j];
        }
      }
    }
    // Multiplies two complex arrays, storing the result in the first operand
    void multiply(fftw_complex* array1, fftw_complex* array2) {
      __m256d normalizationFactor;
      __m256d array1Real;
      __m256d array1Imag;
      __m256d array2Real;
      __m256d array2Imag;
      __m256d realResult;
      __m256d imagResult;
      double real[4];
      double imag[4];
      for (int i = 0; i < cells->height * (cells->width / 2 + 1); i += 4) {
        normalizationFactor = _mm256_set1_pd(cells->height * cells->width);
        array1Real = _mm256_set_pd(array1[i + 3][0], array1[i + 2][0], array1[i + 1][0], array1[i][0]);
        array1Imag = _mm256_set_pd(array1[i + 3][1], array1[i + 2][1], array1[i + 1][1], array1[i][1]);
        array2Real = _mm256_set_pd(array2[i + 3][0], array2[i + 2][0], array2[i + 1][0], array2[i][0]);
        array2Imag = _mm256_set_pd(array2[i + 3][1], array2[i + 2][1], array2[i + 1][1], array2[i][1]);
        realResult = _mm256_sub_pd(_mm256_mul_pd(array1Real, array2Real), _mm256_mul_pd(array1Imag, array2Imag));
        imagResult = _mm256_add_pd(_mm256_mul_pd(array1Real, array2Imag), _mm256_mul_pd(array1Imag, array2Real));
        realResult = _mm256_div_pd(realResult, normalizationFactor);
        imagResult = _mm256_div_pd(imagResult, normalizationFactor);
        _mm256_storeu_pd(real, realResult);
        _mm256_storeu_pd(imag, imagResult);
        for (int j = 0; j < 4; j++) {
          array1[i + j][0] = real[j];
          array1[i + j][1] = imag[j];
        }
      }
    }
    // Calculates all the convolutions, shifts them, and transforms them
    void initialize() {
      for (int i = 0; i < numOrientations; i++) {
        calculateDistanceCoefficients(cells->orientations[i], distanceCoefficients[i]);
        shiftConvolution(distanceCoefficients[i], distanceCoefficientsPadded[i], SEARCH_RADIUS, cells->height, cells->width);
        fftw_execute(distanceCoefficientsFFT[i]);
      }
    }
};

const char* cellTypeToString(CellType type);

void advanceCells(Cells cells, int* searchOffsets, int offsetLength);

// Turn a 2D array of cells into a 1D array of bytes (i.e. for dumping to a file)
unsigned char* serializeCells(Cells cells);

// Inverse of serializeCells
Cells readCells(unsigned char* serializedCells);

void renderCells(Cells cells, SDL_Renderer* renderer, int xOffset, int yOffset, float zoomFactor);
