#include <SDL2/SDL.h>
#include <SDL2/SDL_pixels.h>
#include <SDL2/SDL_rect.h>
#include <SDL2/SDL_render.h>
#include <SDL2/SDL_surface.h>
#include <SDL2/SDL_ttf.h>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <immintrin.h>
#include <iostream>
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
  for (int i = 0; i < currentState.height; i++) {
    for (int j = 0; j < currentState.width; j++) {
      Cell currentCell = currentState.cells[i * currentState.width + j];
      // Serialize the cell type, followed by the state
      for (int k = 0; k < sizeof(CellType); k++) {
        serializedData[index] = (unsigned char) (currentCell.type >> (k * 8));
        index++;
      }
      for (int k = 0; k < sizeof(uint); k++) {
        serializedData[index] = (unsigned char) (currentCell.state >> (k * 8));
        index++;
      }
    }
  }
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
  for (int i = 0; i < cells.height; i++) {
    for (int j = 0; j < cells.width; j++) {
      Cell currentCell;
      currentCell.type = (CellType) 0;
      for (int k = 0; k < sizeof(CellType); k++) {
        currentCell.type = (CellType) (currentCell.type | (serializedData[index] << k * 8));
        index++;
      }
      currentCell.state = 0;
      for (int k = 0; k < sizeof(CellType); k++) {
        currentCell.state = currentCell.state | (serializedData[index] << k * 8);
        index++;
      }
      cells.cells[i * cells.width + j] = currentCell;
    }
  }
  return cells;
}

void saveCellsToFile(Cells cells, char* fileName) {
  unsigned char* data = serializeCells(cells);
  std::ofstream outputStream;
  outputStream.open(fileName);
  outputStream.write((char*) data, getSizeOfData(cells));
  outputStream.close();
  delete[] data;
}

Cells readCellsFromFile(char* fileName) {
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

void advanceCells(Cells* currentState, fftw_complex* distanceCoefficients, double* stateArray, fftw_complex* stateArrayTransformed, double* distanceArray, fftw_plan stateArrayFFT, fftw_plan stateArrayIFFT) {
  // Calculate the FFT of the stateArray
  fftw_execute(stateArrayFFT);
  // Multiply the respective elements of stateArray and distanceCoefficients (distanceCoefficients is already transformed)
  // We use AVX intrinsics here for a ~4x speedup
  for (int i = 0; i < currentState->height * (currentState->width / 2 + 1); i += 4) {
    __m256d normalizationFactor = _mm256_set1_pd(currentState->height * currentState->width);
    __m256d realStateArray = _mm256_set_pd(stateArrayTransformed[i + 3][0], stateArrayTransformed[i + 2][0], stateArrayTransformed[i + 1][0], stateArrayTransformed[i][0]);
    __m256d imagStateArray = _mm256_set_pd(stateArrayTransformed[i + 3][1], stateArrayTransformed[i + 2][1], stateArrayTransformed[i + 1][1], stateArrayTransformed[i][1]);
    __m256d realDistanceCoefficient = _mm256_set_pd(distanceCoefficients[i + 3][0], distanceCoefficients[i + 2][0], distanceCoefficients[i + 1][0], distanceCoefficients[i][0]);
    __m256d imagDistanceCoefficient = _mm256_set_pd(distanceCoefficients[i + 3][1], distanceCoefficients[i + 2][1], distanceCoefficients[i + 1][1], distanceCoefficients[i][1]);
    __m256d realResult = _mm256_sub_pd(_mm256_mul_pd(realStateArray, realDistanceCoefficient), _mm256_mul_pd(imagStateArray, imagDistanceCoefficient));
    __m256d imagResult = _mm256_add_pd(_mm256_mul_pd(realStateArray, imagDistanceCoefficient), _mm256_mul_pd(imagStateArray, realDistanceCoefficient));
    realResult = _mm256_div_pd(realResult, normalizationFactor);
    imagResult = _mm256_div_pd(imagResult, normalizationFactor);
    double real[4];
    double imag[4];
    _mm256_storeu_pd(real, realResult);
    _mm256_storeu_pd(imag, imagResult);
    for (int j = 0; j < 4; j++) {
      stateArrayTransformed[i + j][0] = real[j];
      stateArrayTransformed[i + j][1] = imag[j];
    }
  }
  fftw_execute(stateArrayIFFT);
  // TODO: multithread this loop
  for (int i = 0; i < currentState->height; i++) {
    for (int j = 0; j < currentState->width; j++) {
      double neighborCount = distanceArray[i * currentState->width + j];
      Cell currentCell = currentState->cells[i * currentState->width + j];
      Cell newCell = currentCell;
      // Pacemaker cells are constantly in their action potential
      if (currentCell.type == CellType::Pacemaker) {
        newCell.type = CellType::Pacemaker;
        if (currentCell.state == 0) {
          newCell.state = AP_DURATION;
        }
        else {
          newCell.state = currentCell.state - 1;
        }
        stateArray[i * currentState->width + j] = newCell.state;
      }
      // Resting tissue reactivates into normal tissue eventually
      else if (currentCell.type == CellType::RestingTissue) {
        if (currentCell.state == 1) {
          newCell.type = CellType::Tissue;
          newCell.state = 0;
        }
        else {
          newCell.type = CellType::RestingTissue;
          newCell.state = currentCell.state - 1;
        }
      }
      // If condition here so that more enum variants may be added
      else if (currentCell.type == CellType::Tissue) {
        bool doChange = true;
        // Only search if the cell is not excited
        if (currentCell.state == 0) {
          doChange = false;
          // Search in a radius around the current cell
          doChange = neighborCount >= AP_THRESHOLD;
        }
        if (doChange) {
          if (currentCell.state == 0) {
            newCell.type = CellType::Tissue;
            newCell.state = AP_DURATION;
          }
          else if (currentCell.state == 1) {
            newCell.type = CellType::RestingTissue;
            newCell.state = REST_DURATION;
          }
          else {
            newCell.type = CellType::Tissue;
            newCell.state = currentCell.state - 1;
          }
          if (newCell.type == CellType::RestingTissue) {
            stateArray[i * currentState->width + j] = 0;
          }
          else {
            stateArray[i * currentState->width + j] = newCell.state;
          }
        }
      }
      currentState->cells[i * currentState->width + j] = newCell;
    }
  }
}

// Renders the cells at a (currently) 1-1 ratio of cells to pixels
void renderCells(Cells cells, SDL_Renderer* render, float xOffset, float yOffset, float zoomFactor, int selectedCellI, int selectedCellJ) {
  // TODO: render different colours depending on cell state
  SDL_SetRenderDrawColor(render, 0, 0, 0, 0);
  SDL_RenderClear(render);
  SDL_SetRenderDrawColor(render, 255, 0, 0, 255);
  SDL_FRect cell;
  bool hasSelectedCell = false;
  Cell selectedCell;
  for (int i = 0; i < cells.height; i++) {
    for (int j = 0; j < cells.width; j++) {
      SDL_SetRenderDrawColor(render, 255, 0, 0, 255);
      Cell currentCell = cells.cells[i * cells.width + j];
      if (i == selectedCellI && j == selectedCellJ) {
        selectedCell = cells.cells[selectedCellI * cells.width + selectedCellJ];
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
  if (hasSelectedCell) {
    // TODO: automatic file location OR have a font folder in the project
    TTF_Font* font = TTF_OpenFont("/usr/share/fonts/Ubuntu.ttf", 32);
    SDL_Color textColor = {255, 255, 255, 255};
    char* message = new char[100];
    snprintf(message, 100, "Cell type: %s  Cell state: %d", cellTypeToString(selectedCell.type), selectedCell.state);
    SDL_Surface* textSurface = TTF_RenderText_Solid(font, message, textColor);
    delete[] message;
    SDL_Texture* textTexture = SDL_CreateTextureFromSurface(render, textSurface);
    SDL_Rect textRect = {SIZE - textSurface->w, 0, textSurface->w, textSurface->h};
    SDL_RenderCopy(render, textTexture, NULL, &textRect);
  }
  SDL_RenderPresent(render);
}

