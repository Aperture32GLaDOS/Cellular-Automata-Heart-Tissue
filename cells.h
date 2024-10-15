#include <cstdint>
#include <sys/types.h>
#include <SDL2/SDL_render.h>
#define SIZE 1000

enum CellType{
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
};

struct Cells {
  uint width;
  uint height;
  Cell** cells;
};

void advanceCells(Cells currentState);

// Turn a 2D array of cells into a 1D array of bytes (i.e. for dumping to a file)
unsigned char* serializeCells(Cells currentState);

// Inverse of serializeCells
Cells readCells(unsigned char* serializedCells);

void renderCells(Cells cells, SDL_Renderer* renderer, int xOffset, int yOffset, float zoomFactor);
