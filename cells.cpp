#include <SDL2/SDL_rect.h>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <tuple>
#include <vector>
#include <fftw3.h>
#include "cells.h"

uint getSizeOfData(Cells data) {
  return sizeof(uint) * 2 + sizeof(Cell) * data.height * data.width;
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

void advanceCells(Cells currentState, fftw_complex* distanceCoefficients, double* stateArray, fftw_complex* stateArrayTransformed, double* distanceArray, fftw_plan stateArrayFFT, fftw_plan stateArrayIFFT) {
  std::vector<std::tuple<int, int>> cellsToChange;
  std::vector<Cell> newCells;
  // Calculate the FFT of the stateArray
  fftw_execute(stateArrayFFT);
  // Multiply the respective elements of stateArray and distanceCoefficients (distanceCoefficients is already transformed)
  for (int i = 0; i < currentState.height; i++) {
    // Only go up to about half the width, since beyond that is just complex conjugates
    for (int j = 0; j < (currentState.width / 2 + 1); j++) {
      fftw_complex stateArrayElement = {0.0, 0.0};
      fftw_complex distanceCoefficient = {0.0, 0.0};
      stateArrayElement[0] = stateArrayTransformed[i * (currentState.width / 2 + 1) + j][0];
      stateArrayElement[1] = stateArrayTransformed[i * (currentState.width / 2 + 1) + j][1];
      distanceCoefficient[0] = distanceCoefficients[i * (currentState.width / 2 + 1) + j][0];
      distanceCoefficient[1] = distanceCoefficients[i * (currentState.width / 2 + 1) + j][1];
      // Multiply the two complex numbers stateArrayElement and distanceCoefficient
      // Note the constant division by the size of the array, for normalization purposes
      stateArrayTransformed[i * (currentState.width / 2 + 1) + j][0] = (stateArrayElement[0] * distanceCoefficient[0] - stateArrayElement[1] * distanceCoefficient[1])
        / (currentState.height * currentState.width);
      stateArrayTransformed[i * (currentState.width / 2 + 1) + j][1] = (stateArrayElement[0] * distanceCoefficient[1] + stateArrayElement[1] * distanceCoefficient[0])
        / (currentState.height * currentState.width);
    }
  }
  // Take the inverse FFT to get the neighbor counts
  fftw_execute(stateArrayIFFT);
  for (int i = 0; i < currentState.height; i++) {
    for (int j = 0; j < currentState.width; j++) {
      double neighborCount = distanceArray[i * currentState.width + j];
      Cell currentCell = currentState.cells[i * currentState.width + j];
      Cell newCell;
      // Pacemaker cells are constantly in their action potential
      if (currentCell.type == CellType::Pacemaker) {
        cellsToChange.push_back(std::make_tuple(i, j));
        newCell.type = CellType::Pacemaker;
        if (currentCell.state == 0) {
          newCell.state = AP_DURATION;
        }
        else {
          newCell.state = currentCell.state - 1;
        }
        stateArray[i * currentState.width + j] = newCell.state;
        newCells.push_back(newCell);
      }
      // Resting tissue reactivates into normal tissue eventually
      else if (currentCell.type == CellType::RestingTissue) {
        cellsToChange.push_back(std::make_tuple(i, j));
        if (currentCell.state == 1) {
          newCell.type = CellType::Tissue;
          newCell.state = 0;
        }
        else {
          newCell.type = CellType::RestingTissue;
          newCell.state = currentCell.state - 1;
        }
        newCells.push_back(newCell);
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
          Cell newCell;
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
            stateArray[i * currentState.width + j] = 0;
          }
          else {
            stateArray[i * currentState.width + j] = newCell.state;
          }
          cellsToChange.push_back(std::make_tuple(i, j));
          newCells.push_back(newCell);
        }
      }
      if (false) {
        // Push the location of the cell and its new state
        cellsToChange.push_back(std::make_tuple(i, j));
        Cell newCell;
        newCells.push_back(newCell);
      }
    }
  }

  // Set all changed states
  for (int i = 0; i < cellsToChange.size(); i++) {
    std::tuple<int, int> cellLocation = cellsToChange[i];
    currentState.cells[std::get<0>(cellLocation) * currentState.width + std::get<1>(cellLocation)] = newCells[i];
  }
}

// Renders the cells at a (currently) 1-1 ratio of cells to pixels
void renderCells(Cells cells, SDL_Renderer* render, float xOffset, float yOffset, float zoomFactor, int selectedCellI, int selectedCellJ) {
  // TODO: render different colours depending on cell state
  SDL_SetRenderDrawColor(render, 0, 0, 0, 0);
  SDL_RenderClear(render);
  SDL_SetRenderDrawColor(render, 255, 0, 0, 255);
  SDL_FRect cell;
  for (int i = 0; i < cells.height; i++) {
    for (int j = 0; j < cells.width; j++) {
      SDL_SetRenderDrawColor(render, 255, 0, 0, 255);
      Cell currentCell = cells.cells[i * cells.width + j];
      if (i == selectedCellI && j == selectedCellJ) {
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
  SDL_RenderPresent(render);
}

