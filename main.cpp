#include "cells.cpp"

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
