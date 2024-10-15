#include "cells.cpp"
#include <SDL2/SDL.h>
#include <SDL2/SDL_events.h>
#include <SDL2/SDL_syswm.h>
#include <SDL2/SDL_video.h>
#include <chrono>
#include <thread>

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
  long int startTime;
  long int elapsedTime;
  while (!quit) {
    // TODO: poll events while sleeping
    // TODO: zoom in and pan around
    while (SDL_PollEvent(&currentEvent) != 0) {
      if (currentEvent.type == SDL_QUIT) {
        quit = true;
      }
    }
    startTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    renderCells(cells, renderer, 0, 0, 1);
    advanceCells(cells);
    elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count() - startTime;
    if (elapsedTime <= 500) {
      std::this_thread::sleep_for(std::chrono::milliseconds(500 - elapsedTime));
    }
  }
  SDL_DestroyWindow(window);
  SDL_Quit();
  return 0;
}
