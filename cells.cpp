#include <SDL2/SDL_rect.h>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <tuple>
#include <vector>
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
      Cell currentCell = currentState.cells[i * currentState.height + j];
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
      cells.cells[i * cells.height + j] = currentCell;
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

void advanceCells(Cells currentState, int* searchOffsets, int offsetLength, double* searchCoefficients) {
  std::vector<std::tuple<int, int>> cellsToChange;
  std::vector<Cell> newCells;
  for (int i = 0; i < currentState.height; i++) {
    for (int j = 0; j < currentState.width; j++) {
      // TODO: optimize the hell out of this!!!
      // Potential ideas:
      // 1.
      //  GPU utilization - would definitely make it go like the clappers, but will be VERY hard to implement (complete rewrite of rendering engine,
      //  and working with Vulkan is notoriously difficult)
      // 2.
      //  Efficient convolution algorithms - the search algorithm can be described as a convolution problem, which have very clever algorithms
      //  To solve them in very quick time
      double neighborCount = 0;
      Cell currentCell = currentState.cells[i * currentState.height + j];
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
          for (int k = 0; k < offsetLength; k+=2) {
            double coefficient = searchCoefficients[k / 2];
            int searchI = searchOffsets[k];
            int searchJ = searchOffsets[k + 1];
            if (i + searchI >= 0 && i + searchI < currentState.height && j + searchJ >= 0 && j + searchJ < currentState.width) {
              Cell neighbor = currentState.cells[(i + searchI) * currentState.height + (j + searchJ)];
              // If a cell is resting, then it cannot excite surrounding tissue
              if (neighbor.type != CellType::RestingTissue) {
                neighborCount += ((double) neighbor.state) * coefficient;
              }
            }
          }
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
    currentState.cells[std::get<0>(cellLocation) * currentState.height + std::get<1>(cellLocation)] = newCells[i];
  }
}

// Renders the cells at a (currently) 1-1 ratio of cells to pixels
void renderCells(Cells cells, SDL_Renderer* render, float xOffset, float yOffset, float zoomFactor) {
  // TODO: render different colours depending on cell state
  SDL_SetRenderDrawColor(render, 0, 0, 0, 0);
  SDL_RenderClear(render);
  SDL_SetRenderDrawColor(render, 255, 0, 0, 255);
  SDL_FRect cell;
  for (int i = 0; i < cells.height; i++) {
    for (int j = 0; j < cells.width; j++) {
      if (cells.cells[i * cells.height + j].state > 0 && cells.cells[i * cells.height + j].type != CellType::RestingTissue) {
        cell.x = (j + xOffset) * zoomFactor;
        cell.y = (i + yOffset) * zoomFactor;
        cell.w = zoomFactor;
        cell.h = zoomFactor;
        SDL_RenderFillRectF(render, &cell);
      }
    }
  }
  SDL_RenderPresent(render);
}

