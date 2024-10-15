#include "cells.cpp"
#include <SDL2/SDL.h>
#include <SDL2/SDL_events.h>
#include <SDL2/SDL_syswm.h>
#include <SDL2/SDL_video.h>

SDL_Window* window = NULL;
SDL_Renderer* renderer = NULL;

int main (int argc, char *argv[]) {
  if (SDL_Init(SDL_INIT_VIDEO) < 0) {
    return 0;
  }
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
      if (std::rand() % 100 <= 50) {
        newCell.state = 0;
      }
      else {
        newCell.state = 1;
      }
      cells.cells[i][j] = newCell;
    }
  }
  window = SDL_CreateWindow("Heart Tissue", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
      cells.width, cells.height, SDL_WINDOW_SHOWN);
  renderer = SDL_CreateRenderer(window, -1, 0);
  SDL_RenderPresent(renderer);
  bool quit = false;
  SDL_Event currentEvent;
  while (!quit) {
    while (SDL_PollEvent(&currentEvent) != 0) {
      renderCells(cells, renderer);
      advanceCells(cells);
      if (currentEvent.type == SDL_QUIT) {
        quit = true;
      }
    }
  }
  SDL_DestroyWindow(window);
  SDL_Quit();
  return 0;
}
