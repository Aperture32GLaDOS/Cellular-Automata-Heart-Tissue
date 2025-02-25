#include <cstdint>
#include <sys/types.h>
#include <SDL2/SDL_render.h>
#define SIZE 1024
#define SEARCH_RADIUS 35
#define AP_DURATION 32
#define REST_DURATION 8
#define AP_THRESHOLD 16

enum CellType {
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
  // 2D array is represented as one contiguous block of memory, for performance reasons
  // (indexed [i][j] would be [i * height + j])
  Cell* cells;
};

const char* cellTypeToString(CellType type);

void advanceCells(Cells currentState, int* searchOffsets, int offsetLength);

// Turn a 2D array of cells into a 1D array of bytes (i.e. for dumping to a file)
unsigned char* serializeCells(Cells currentState);

// Inverse of serializeCells
Cells readCells(unsigned char* serializedCells);

void renderCells(Cells cells, SDL_Renderer* renderer, int xOffset, int yOffset, float zoomFactor);
