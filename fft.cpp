#include <complex>
#include "fft.h"

// Calculates the FFT of a 1D array (in) in-place
bool FFT(std::complex<double>* in, int length, int direction) {
  // If length is not a power of two, return false
  if ((length & (length - 1)) != 0) return false;
  // Base case of the recursion
  if (length <= 1) return true;
  std::complex<double>* even = new std::complex<double>[length / 2];
  std::complex<double>* odd= new std::complex<double>[length / 2];
  for (int i = 0; i < length / 2; i++) {
      even[i] = in[i * 2];
      odd[i] = in[i * 2 + 1];
  }
  // Recurse
  FFT(even, length / 2, direction);
  FFT(odd, length / 2, direction);
  // Combine the results of the recursion
  for (int i = 0; i < length / 2; i++) {
    std::complex<double> t = std::polar(-1.0, -2 * PI * i * direction / length) * odd[i];
    in[i] = even[i] + t;
    in[i + length / 2] = even[i] - t;
    if (direction == -1) {
      in[i] /= 2;
      in[i + length / 2] /= 2;
    }
  }
  return true;
}

// Transforms the complex 2D array inplace
// This assumes that width and height are both powers of 2
bool FFT2D(std::complex<double>* in, int width, int height, int direction) {
  // Transform the rows
  std::complex<double>* toTransform = new std::complex<double>[width];
  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      toTransform[j] = in[i * height + j];
    }
    if (!FFT(toTransform, width, direction)) return false;
    for (int j = 0; j < width; j++) {
      in[i * height + j] = toTransform[j];
    }
  }
  // Transform the columns
  delete[] toTransform;
  toTransform = new std::complex<double>[height];
  for (int j = 0; j < width; j++) {
    for (int i = 0; i < height; i++) {
      toTransform[i] = in[i * height + j];
    }
    if (!FFT(toTransform, height, direction)) return false;
    for (int i = 0; i < height; i++) {
      in[i * height + j] = toTransform[i];
    }
  }
  delete[] toTransform;
  return true;
}
