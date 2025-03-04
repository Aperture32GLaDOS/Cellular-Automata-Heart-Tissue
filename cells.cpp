#include <SDL2/SDL.h>
#include <SDL2/SDL_pixels.h>
#include <SDL2/SDL_rect.h>
#include <SDL2/SDL_render.h>
#include <SDL2/SDL_surface.h>
#include <SDL2/SDL_ttf.h>
#include <SDL2/SDL_error.h>
#include <chrono>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <immintrin.h>
#include <iostream>
#include <thread>
#include <x86intrin.h>
#include <fftw3.h>
#include "cells.h"

uint getSizeOfData(Cells data) {
  return sizeof(uint) * 2 + sizeof(Cell) * data.height * data.width;
}

const char* cellTypeToString(CellType type) {
  switch (type) {
    case CellType::Tissue:
      return "Normal Cell";
    case CellType::Pacemaker:
      return "Pacemaker Cell";
    case CellType::RestingTissue:
      return "Resting Cell";
    default:
      std::cout << "ERROR!!! Unkown Cell Type Encountered";
      return NULL;
  }
}

unsigned char* serializeCells(Cells currentState) {
  unsigned char* serializedData = new uint8_t[getSizeOfData(currentState)];
  // Global index to keep track of where in the serializedData array we are (to avoid annoying index calculations)
  uint index = 0;
  // Serialize the width and height
  for (int i = 0; i < sizeof(uint); i++) {
    serializedData[index] = (unsigned char) (currentState.width >> (i * 8));
    index++;
  }
  for (int i = 0; i < sizeof(uint); i++) {
    serializedData[index] = (unsigned char) (currentState.height >> (i * 8));
    index++;
  }
  // Serialize all the actual data
  memcpy(&serializedData[index], currentState.cells, sizeof(Cell) * currentState.width * currentState.height);
  return serializedData;
}

Cells readCells(unsigned char* serializedData) {
  Cells cells;
  uint index = 0;
  cells.width = 0;
  for (int i = 0; i < sizeof(uint); i++) {
    cells.width = cells.width | (serializedData[index] << i * 8);
    index++;
  }
  cells.height = 0;
  for (int i = 0; i < sizeof(uint); i++) {
    cells.height = cells.height | (serializedData[index] << i * 8);
    index++;
  }
  cells.cells = new Cell[cells.height * cells.width];
  memcpy(cells.cells, &serializedData[index], sizeof(Cell) * cells.width * cells.height);
  return cells;
}

void saveCellsToFile(Cells cells, const char* fileName) {
  unsigned char* data = serializeCells(cells);
  std::ofstream outputStream;
  outputStream.open(fileName);
  outputStream.write((char*) data, getSizeOfData(cells));
  outputStream.close();
  delete[] data;
}

Cells readCellsFromFile(const char* fileName) {
  std::ifstream inputStream;
  inputStream.open(fileName);
  inputStream.seekg(0, std::ios::end);
  size_t length = inputStream.tellg();
  inputStream.seekg(0, std::ios::beg);
  unsigned char* data = new unsigned char[length];
  inputStream.read((char*) data, length);
  Cells output = readCells(data);
  // Free the memory space
  delete[] data;
  return output;
}

void updateCellsArea(Cells* currentState, double* distanceArray, double* stateArray, int start, int end) {
  int pacemaker = 0;
  int tissue = 1;
  int restingTissue = 2;
  __m256i pacemakerAVX = _mm256_set1_epi32(pacemaker);
  __m256i tissueAVX = _mm256_set1_epi32(tissue);
  __m256i restingTissueAVX = _mm256_set1_epi32(restingTissue);
  __m256d firstHalfNeighbours;
  __m256d secondHalfNeighbours;
  __m128 firstHalfNeighboursSingle;
  __m128 secondHalfNeighboursSingle;
  __m256 neighbours;
  __m256i neighboursRounded;
  __m256i cellStates;
  __m256i cellTypes;
  __m256i stateArrayAVX;
  __m256i isPacemaker;
  __m256i isTissue;
  __m256i isResting;
  __m256i wasActive;
  __m256i isZeroState;
  __m256i isAboveThreshold;
  __m256i negativeBitAVX = _mm256_set1_epi32(1 << 31);
  __m256 floatNegativeBitAVX = *((__m256*) &negativeBitAVX);
  __m256i zeroAVX = _mm256_set1_epi32(0);
  __m256i oneAVX = _mm256_set1_epi32(1);
  __m256i allOneBitsAVX = _mm256_sub_epi32(zeroAVX, oneAVX);
  __m256i restingDurationAVX = _mm256_set1_epi32(REST_DURATION);
  __m256i maxStateAVX = _mm256_set1_epi32(AP_DURATION);
  __m256 thresholdAVX = _mm256_set1_ps(AP_THRESHOLD);
  __m256i restingToNormal = _mm256_sub_epi32(tissueAVX, restingTissueAVX);
  __m256i normalToResting = _mm256_sub_epi32(restingTissueAVX, tissueAVX);
  int cellStatesUpdated[8];
  int cellTypesUpdatedInt[8];
  int stateArrayUpdated[8];
  CellType* cellTypesUpdated = (CellType*) cellTypesUpdatedInt;
  for (int i = start; i < end; i+=8) {
    // Convert neighbourhood counts to single-precision, so 8 of them can be stored at once
    firstHalfNeighbours = _mm256_set_pd(distanceArray[(i + 3)], distanceArray[(i + 2)], distanceArray[(i + 1)], distanceArray[i]);
    secondHalfNeighbours = _mm256_set_pd(distanceArray[(i + 7)], distanceArray[(i + 6)], distanceArray[(i + 5)], distanceArray[(i + 4)]);
    firstHalfNeighboursSingle = _mm256_cvtpd_ps(firstHalfNeighbours);
    secondHalfNeighboursSingle = _mm256_cvtpd_ps(secondHalfNeighbours);
    neighbours = _mm256_castps128_ps256(firstHalfNeighboursSingle);
    neighbours = _mm256_insertf128_ps(neighbours, secondHalfNeighboursSingle, 1);
    neighbours = _mm256_sub_ps(neighbours, thresholdAVX);
    neighbours = _mm256_and_ps(neighbours, floatNegativeBitAVX);
    // We now have 0 if the count is above the threshold, and 1 otherwise, but this 
    // is on the first bit
    neighboursRounded = *((__m256i*) &neighbours);
    // Shift 31 bits to the right
    neighboursRounded = _mm256_srli_epi32(neighboursRounded, 31);
    neighboursRounded = _mm256_sub_epi32(neighboursRounded, oneAVX);

    // Store the next eight cell states
    cellStates = _mm256_set_epi32(currentState->cells[(i + 7)].state, currentState->cells[(i + 6)].state, currentState->cells[(i + 5)].state, currentState->cells[(i + 4)].state, currentState->cells[(i + 3)].state, currentState->cells[(i + 2)].state, currentState->cells[(i + 1)].state, currentState->cells[i].state);
    // And the cell types
    cellTypes = _mm256_set_epi32(currentState->cells[(i + 7)].type, currentState->cells[(i + 6)].type, currentState->cells[(i + 5)].type, currentState->cells[(i + 4)].type, currentState->cells[(i + 3)].type, currentState->cells[(i + 2)].type, currentState->cells[(i + 1)].type, currentState->cells[i].type);
    // This value is 0 if the cell is a pacemaker, and some other value otherwise
    isPacemaker = _mm256_xor_si256(pacemakerAVX, cellTypes);
    // Clamp value to 0 or 1
    isPacemaker = _mm256_min_epu32(isPacemaker, oneAVX);
    // Subtract one so that if the cell is a pacemaker, the bits are all one, and otherwise the bits are all zero
    isPacemaker = _mm256_sub_epi32(isPacemaker, oneAVX);
    isTissue = _mm256_xor_si256(tissueAVX, cellTypes);
    isTissue = _mm256_min_epu32(isTissue, oneAVX);
    isTissue = _mm256_sub_epi32(isTissue, oneAVX);
    isResting = _mm256_xor_si256(restingTissueAVX, cellTypes);
    isResting = _mm256_min_epu32(isResting, oneAVX);
    isResting = _mm256_sub_epi32(isResting, oneAVX);

    wasActive = _mm256_xor_si256(cellStates, zeroAVX);
    wasActive = _mm256_min_epu32(wasActive, oneAVX);
    wasActive = _mm256_sub_epi32(wasActive, oneAVX);
    wasActive = _mm256_xor_si256(wasActive, allOneBitsAVX);

    // If the cell state is not initially 0, then reduce it by one
    cellStates = _mm256_sub_epi32(cellStates, _mm256_and_si256(wasActive, oneAVX));
    isZeroState = _mm256_xor_si256(cellStates, zeroAVX);
    isZeroState = _mm256_min_epu32(isZeroState, oneAVX);
    isZeroState = _mm256_sub_epi32(isZeroState, oneAVX);
    // If the cell state is 0, and the cell is a pacemaker then set the cell state to AP_DURATION
    cellStates = _mm256_add_epi32(cellStates, _mm256_and_si256(maxStateAVX, _mm256_and_si256(isZeroState, isPacemaker)));
    // If the cell state is 0 and the cell is resting, then the cell is now set to normal tissue
    cellTypes = _mm256_add_epi32(cellTypes, _mm256_and_si256(restingToNormal, _mm256_and_si256(isZeroState, isResting)));
    // If the cell state is 0 and the cell was active, then the cell is now resting tissue with a state of REST_DURATION
    cellStates = _mm256_add_epi32(cellStates, _mm256_and_si256(restingDurationAVX, _mm256_and_si256(isZeroState, _mm256_and_si256(wasActive, isTissue))));
    cellTypes = _mm256_add_epi32(cellTypes, _mm256_and_si256(normalToResting, _mm256_and_si256(isZeroState, _mm256_and_si256(wasActive, isTissue))));
    isTissue = _mm256_xor_si256(tissueAVX, cellTypes);
    isTissue = _mm256_min_epu32(isTissue, oneAVX);
    isTissue = _mm256_sub_epi32(isTissue, oneAVX);
    // Alternatively, then if the cell is normal tissue and not already active, then set the cell's state to AP_DURATION only if the neighbor count is greater than AP_THRESHOLD
    cellStates = _mm256_add_epi32(cellStates, _mm256_and_si256(maxStateAVX, _mm256_and_si256(neighboursRounded, _mm256_and_si256(isZeroState, isTissue))));

    // Only count the cells as part of the state array if they are a pacemaker or normal tissue cell
    stateArrayAVX = _mm256_and_si256(cellStates, _mm256_or_si256(isPacemaker, isTissue));
    _mm256_storeu_si256((__m256i*) cellStatesUpdated, cellStates);
    _mm256_storeu_si256((__m256i*) cellTypesUpdatedInt, cellTypes);
    _mm256_storeu_si256((__m256i*) stateArrayUpdated, stateArrayAVX);

    for (int j = 0; j < 8; j++) {
      currentState->cells[(i + j)].state = cellStatesUpdated[j];
      currentState->cells[(i + j)].type = cellTypesUpdated[j];
      stateArray[i + j] = (double) stateArrayUpdated[j];
    }
  }
}

void advanceCells(Cells* currentState, fftw_complex* distanceCoefficients, double* stateArray, fftw_complex* stateArrayTransformed, double* distanceArray, fftw_plan stateArrayFFT, fftw_plan stateArrayIFFT) {
  auto start = std::chrono::high_resolution_clock::now();
  // Calculate the FFT of the stateArray
  fftw_execute(stateArrayFFT);
  // Multiply the respective elements of stateArray and distanceCoefficients (distanceCoefficients is already transformed)
  // We use AVX intrinsics here for a ~4x speedup
  __m256d normalizationFactor;
  __m256d realStateArray;
  __m256d imagStateArray;
  __m256d realDistanceCoefficient;
  __m256d imagDistanceCoefficient;
  __m256d realResult;
  __m256d imagResult;
  double real[4];
  double imag[4];
  for (int i = 0; i < currentState->height * (currentState->width / 2 + 1); i += 4) {
    normalizationFactor = _mm256_set1_pd(currentState->height * currentState->width);
    realStateArray = _mm256_set_pd(stateArrayTransformed[i + 3][0], stateArrayTransformed[i + 2][0], stateArrayTransformed[i + 1][0], stateArrayTransformed[i][0]);
    imagStateArray = _mm256_set_pd(stateArrayTransformed[i + 3][1], stateArrayTransformed[i + 2][1], stateArrayTransformed[i + 1][1], stateArrayTransformed[i][1]);
    realDistanceCoefficient = _mm256_set_pd(distanceCoefficients[i + 3][0], distanceCoefficients[i + 2][0], distanceCoefficients[i + 1][0], distanceCoefficients[i][0]);
    imagDistanceCoefficient = _mm256_set_pd(distanceCoefficients[i + 3][1], distanceCoefficients[i + 2][1], distanceCoefficients[i + 1][1], distanceCoefficients[i][1]);
    realResult = _mm256_sub_pd(_mm256_mul_pd(realStateArray, realDistanceCoefficient), _mm256_mul_pd(imagStateArray, imagDistanceCoefficient));
    imagResult = _mm256_add_pd(_mm256_mul_pd(realStateArray, imagDistanceCoefficient), _mm256_mul_pd(imagStateArray, realDistanceCoefficient));
    realResult = _mm256_div_pd(realResult, normalizationFactor);
    imagResult = _mm256_div_pd(imagResult, normalizationFactor);
    _mm256_storeu_pd(real, realResult);
    _mm256_storeu_pd(imag, imagResult);
    for (int j = 0; j < 4; j++) {
      stateArrayTransformed[i + j][0] = real[j];
      stateArrayTransformed[i + j][1] = imag[j];
    }
  }
  fftw_execute(stateArrayIFFT);
  // Safe to thread here as mutex is locked when this function is called
  constexpr int NUM_THREADS = 8;
  std::thread threads[NUM_THREADS];
  int delta = (currentState->width * currentState->height) / NUM_THREADS;
  for (int i = 0; i < NUM_THREADS; i++) {
    threads[i] = std::thread(updateCellsArea, currentState, distanceArray, stateArray, delta * i, delta * (i + 1));
  }
  for (int i = 0; i < NUM_THREADS; i++) {
    threads[i].join();
  }
  auto end = std::chrono::high_resolution_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  std::cout << "Time taken to calculate cells: " << elapsed.count() << "ms" << std::endl;
}

void renderCells(Cells cells, SDL_Renderer* render, TTF_Font* font, float xOffset, float yOffset, float zoomFactor, int selectedCellI, int selectedCellJ,
    int firstCornerX, int secondCornerX, int firstCornerY, int secondCornerY) {
  auto start = std::chrono::high_resolution_clock::now();
  SDL_SetRenderDrawColor(render, 0, 0, 0, 0);
  SDL_RenderClear(render);
  SDL_SetRenderDrawColor(render, 255, 0, 0, 255);
  if (firstCornerX > secondCornerX) {
    std::swap(firstCornerX, secondCornerX);
  }
  if (firstCornerY > secondCornerY) {
    std::swap(firstCornerY, secondCornerY);
  }
  int firstCornerYScreenSpace = (firstCornerX + yOffset) * zoomFactor;
  int firstCornerXScreenSpace = (firstCornerY + xOffset) * zoomFactor;
  int secondCornerYScreenSpace = (secondCornerX + yOffset) * zoomFactor;
  int secondCornerXScreenSpace = (secondCornerY + xOffset) * zoomFactor;
  SDL_FRect cell;
  bool hasSelectedCell = false;
  Cell selectedCell;
  for (int i = std::max((int) -yOffset - 1, 0); i < std::min(cells.height, (uint) (cells.height / zoomFactor - yOffset) + 1); i++) {
    for (int j = std::max((int) -xOffset - 1, 0); j < std::min(cells.width, (uint) (cells.width / zoomFactor - xOffset) + 1); j++) {
      SDL_SetRenderDrawColor(render, 255, 0, 0, 255);
      Cell currentCell = cells.cells[i * cells.width + j];
      if (i == selectedCellI && j == selectedCellJ) {
        selectedCell = currentCell;
        hasSelectedCell = true;
        SDL_SetRenderDrawColor(render, 100, 100, 100, 255);
        cell.x = (j + xOffset) * zoomFactor;
        cell.y = (i + yOffset) * zoomFactor;
        cell.w = zoomFactor;
        cell.h = zoomFactor;
        SDL_RenderFillRectF(render, &cell);
      }
      else {
        if (currentCell.state > 0 && currentCell.type != CellType::RestingTissue) {
          if (currentCell.type == Pacemaker) {
            SDL_SetRenderDrawColor(render, 255, 0, 255, 255);
          }
          cell.x = (j + xOffset) * zoomFactor;
          cell.y = (i + yOffset) * zoomFactor;
          cell.w = zoomFactor;
          cell.h = zoomFactor;
          SDL_RenderFillRectF(render, &cell);
        }
      }
    }
  }
  SDL_Rect selectedRect;
  selectedRect.x = firstCornerXScreenSpace;
  selectedRect.y = firstCornerYScreenSpace;
  selectedRect.w = secondCornerXScreenSpace - firstCornerXScreenSpace;
  selectedRect.h = secondCornerYScreenSpace - firstCornerYScreenSpace;
  SDL_SetRenderDrawColor(render, 0, 255, 0, 255);
  if (firstCornerX != secondCornerX && firstCornerY != secondCornerY) {
    SDL_RenderDrawRect(render, &selectedRect);
  }
  if (hasSelectedCell) {
    // TODO: automatic file location OR have a font folder in the project
    SDL_Color textColor = {255, 255, 255, 255};
    char* message = new char[100];
    snprintf(message, 100, "Cell type: %s  Cell state: %d", cellTypeToString(selectedCell.type), selectedCell.state);
    SDL_Surface* textSurface = TTF_RenderText_Solid(font, message, textColor);
    if (textSurface == NULL) {
      std::cout << SDL_GetError() << std::endl;
    }
    else {
      SDL_Texture* textTexture = SDL_CreateTextureFromSurface(render, textSurface);
      SDL_Rect textRect = {SIZE - textSurface->w, 0, textSurface->w, textSurface->h};
      SDL_RenderCopy(render, textTexture, NULL, &textRect);
    }
    delete[] message;
  }
  SDL_RenderPresent(render);
  auto end = std::chrono::high_resolution_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  std::cout << "Time taken to render cells: " << elapsed.count() << "ms" << std::endl;
}

