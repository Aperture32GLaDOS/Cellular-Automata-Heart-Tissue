#include <cstddef>
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
      Cell currentCell = currentState.cells[i][j];
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

Cells readCells(uint8_t* serializedData) {
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
  cells.cells = new Cell*[cells.height];
  for (int i = 0; i < cells.height; i++) {
    cells.cells[i] = new Cell[cells.width];
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
      cells.cells[i][j] = currentCell;
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
  return readCells(data);
}

void advanceCells(Cells currentState) {
  std::vector<std::tuple<int, int>> cellsToChange;
  std::vector<uint> newStates;
  for (int i = 0; i < currentState.height; i++) {
    for (int j = 0; j < currentState.width; j++) {
      // TODO: implement rules for changing states
      // If the cell must be changed,
      if (false) {
        // Push the location of the cell and its new state
        cellsToChange.push_back(std::make_tuple(i, j));
        newStates.push_back(0);
      }
    }
  }

  // Set all changed states
  for (int i = 0; i < cellsToChange.size(); i++) {
    std::tuple<int, int> cellLocation = cellsToChange[i];
    currentState.cells[std::get<0>(cellLocation)][std::get<1>(cellLocation)].state = newStates[i];
  }
}

int main (int argc, char *argv[]) {
  // Declare the 2D plane of cells
  Cells cells;
  cells.width = SIZE;
  cells.height = SIZE;
  cells.cells = new Cell*[cells.height];
  for (int i = 0; i < cells.height; i++) {
    cells.cells[i] = new Cell[cells.width];
    for (int j = 0; j < cells.width; j++) {
      Cell newCell;
      newCell.type = CellType::Tissue;
      newCell.state = 0;
      cells.cells[i][j] = newCell;
    }
  }
  return 0;
}
