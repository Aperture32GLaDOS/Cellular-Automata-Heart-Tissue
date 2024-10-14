#include <iostream>
#include "cells.h"

// TODO: this is quite slow currently, and deep-copies a 1000*1000 array every iteration
Cell** advanceCells(Cell** currentState) {
  Cell** newState = new Cell*[SIZE];
  for (int i = 0; i < SIZE; i++) {
    newState[i] = new Cell[SIZE];
    for (int j = 0; j < SIZE; j++) {
      // TODO: apply advancing rules for cells

    }
  }
  // Free memory used for currentState
  // TODO: remember to remove this when the array is no longer deep-copied
  for (int i = 0; i < SIZE; i++) {
    delete [] currentState[i];
  }
  delete currentState;
  return newState;
}

int main (int argc, char *argv[]) {
  // Declare the 2D plane of cells
  // TODO: make a way of serializing this data, so certain situations can be saved for
  // further research
  Cell** cells = new Cell*[SIZE];
  for (int i = 0; i < SIZE; i++) {
    cells[i] = new Cell[SIZE];
  }
  std::cout << "Hello, World!" << std::endl;
  return 0;
}
