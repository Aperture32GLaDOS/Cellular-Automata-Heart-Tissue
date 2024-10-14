#include <sys/types.h>
#define SIZE 1000

enum CellType{
  // A heart cell here is represented either as a pacemaker cell, or a normal tissue cell
  Pacemaker,
  Tissue
};

// A cell is declared as a struct since it has no methods
struct Cell {
  CellType type;
  // The state is a positive integer
  uint state;
};

Cell** advanceCells(Cell** currentState);
