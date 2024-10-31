#include "cells.cpp"
#include "fft.cpp"
#include <SDL2/SDL.h>
#include <SDL2/SDL_events.h>
#include <SDL2/SDL_keycode.h>
#include <SDL2/SDL_mouse.h>
#include <SDL2/SDL_syswm.h>
#include <SDL2/SDL_video.h>
#include <chrono>
#include <cstdlib>
#include <thread>
#include <mutex>

SDL_Window* window = NULL;
SDL_Renderer* renderer = NULL;


std::mutex mu;

// Updates the cells in a seperate thread, so as to keep the render updates fast
void updateCells(Cells* cells, bool* quit, bool* paused, int* frameTime) {
  long int startTime;
  long int elapsedTime;
  int* searchOffsets = new int[2 * SEARCH_RADIUS * SEARCH_RADIUS];
  int offsetLength = 0;
  int counter = -1;
  for (int i = -(SEARCH_RADIUS / 2); i <= (SEARCH_RADIUS / 2); i++) {
    for (int j = -(SEARCH_RADIUS / 2); j <= (SEARCH_RADIUS / 2); j++) {
      if (i == 0 && j == 0) continue;
      if ((i * i) + (j * j) <= (SEARCH_RADIUS / 2) * (SEARCH_RADIUS / 2)) {
        searchOffsets[offsetLength] = i;
        searchOffsets[offsetLength + 1] = j;
        offsetLength += 2;
      }
    }
  }
  offsetLength -= 2;
  while (!(*quit)) {
    startTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    // Lock the mutex, as data is being written
    std::unique_lock<std::mutex> lock(mu);
    advanceCells(*cells, searchOffsets, offsetLength);
    lock.unlock();
    elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count() - startTime;
    if (elapsedTime <= *frameTime) {
      std::this_thread::sleep_for(std::chrono::milliseconds(*frameTime - elapsedTime));
    }
    while (*paused && !(*quit)) {
      std::this_thread::sleep_for(std::chrono::milliseconds(250));
    }
  }
}


int main (int argc, char *argv[]) {
  if (SDL_Init(SDL_INIT_VIDEO) < 0) {
    return 0;
  }
  // Test 2D fft
  std::complex<double>* array = new std::complex<double>[32];
  for (int i = 0; i < 32; i++) {
    if (i % 2 == 0) {
      array[i * 32] = std::complex(420.0, 0.0);
    }
    else {
      array[i * 32] = std::complex(0.0, 0.0);
    }
  }
  std::cout << array[16] << std::endl;
  FFT(array, 32, 1);
  std::cout << array[16] << std::endl;
  FFT(array, 32, -1);
  std::cout << array[16] << std::endl;
  // Declare the 2D plane of cells
  Cells cells;
  cells.width = SIZE;
  cells.height = SIZE;
  cells.cells = new Cell[cells.height * cells.width];
  for (int i = 0; i < cells.height; i++) {
    for (int j = 0; j < cells.width; j++) {
      Cell newCell;
      newCell.type = CellType::Tissue;
      newCell.state = 0;
      cells.cells[i * cells.height + j] = newCell;
    }
  }
  for (int i = -3; i <= 3; i++) {
    for (int j = -3; j <= 3; j++) {
      cells.cells[(cells.height / 2 + i) * cells.height + (cells.width / 2 + j)].type = CellType::Pacemaker;
    }
  }
  window = SDL_CreateWindow("Heart Tissue", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
      cells.width, cells.height, SDL_WINDOW_SHOWN);
  renderer = SDL_CreateRenderer(window, -1, 0);
  SDL_RenderPresent(renderer);
  bool quit = false;
  bool paused = false;
  SDL_Event currentEvent;
  float xOffset = 0;
  float yOffset = 0;
  float zoomFactor = 1;
  int mousePosX;
  int mousePosY;
  int frameTime = 500;
  std::thread updateThread(updateCells, &cells, &quit, &paused, &frameTime);
  while (!quit) {
    while (SDL_PollEvent(&currentEvent) != 0) {
      if (currentEvent.type == SDL_QUIT) {
        quit = true;
      }
      else if (currentEvent.type == SDL_KEYDOWN) {
        if (currentEvent.key.keysym.sym == SDLK_w) {
          yOffset += 10;
        }
        else if (currentEvent.key.keysym.sym == SDLK_s) {
          yOffset -= 10;
        }
        else if (currentEvent.key.keysym.sym == SDLK_a) {
          xOffset += 10;
        }
        else if (currentEvent.key.keysym.sym == SDLK_d) {
          xOffset -= 10;
        }
        else if (currentEvent.key.keysym.sym == SDLK_SPACE) {
          paused = !paused;
        }
        else if (currentEvent.key.keysym.sym == SDLK_EQUALS) {
          frameTime -= 50;
          if (frameTime == 0) {
            frameTime = 50;
          }
        }
        // Saves the current state to a file
        else if (currentEvent.key.keysym.sym == SDLK_F1) {
          std::unique_lock<std::mutex> lock(mu);
          saveCellsToFile(cells, "cells.dmp");
          lock.unlock();
        }
        else if (currentEvent.key.keysym.sym == SDLK_F2) {
          std::unique_lock<std::mutex> lock(mu);
          // Delete the old array so as to avoid a memory leak
          delete[] cells.cells;
          cells = readCellsFromFile("cells.dmp");
          SDL_SetWindowSize(window, cells.width, cells.height);
          lock.unlock();
        }
        else if (currentEvent.key.keysym.sym == SDLK_MINUS) {
          frameTime += 50;
        }
      }
      else if (currentEvent.type == SDL_MOUSEWHEEL) {
        // Zooms into where the user's cursor is
        SDL_GetMouseState(&mousePosX, &mousePosY);
        // Calculate the selected X and Y pixel
        float selectedX = (mousePosX / zoomFactor) - xOffset;
        float selectedY = (mousePosY / zoomFactor) - yOffset;
        zoomFactor = zoomFactor + currentEvent.wheel.y;
        if (zoomFactor < 1) {
          zoomFactor = 1;
        }
        else {
          // Rearrange the equation to calculate the new x and y offsets
          xOffset = (mousePosX / zoomFactor) - selectedX;
          yOffset = (mousePosY / zoomFactor) - selectedY;
        }
      }
      // When the user presses the mouse button, change the state of the cellular automata
      else if (currentEvent.type == SDL_MOUSEBUTTONDOWN) {
        SDL_GetMouseState(&mousePosX, &mousePosY);
        std::unique_lock<std::mutex> lock(mu);
        Cell* selectedCell = &cells.cells[((int) ((mousePosY / zoomFactor) - yOffset)) * cells.height + (int) ((mousePosX / zoomFactor) - xOffset)];
        if (currentEvent.button.button == SDL_BUTTON_LEFT) {
          selectedCell->state = AP_DURATION;
        }
        else if (currentEvent.button.button == SDL_BUTTON_RIGHT) {
          if (selectedCell-> state != 0) selectedCell->state = selectedCell->state - 1;
        }
        // TODO: change tissue type on shift-right click (or similar)
        lock.unlock();
      }
    }
    std::unique_lock<std::mutex> lock(mu);
    renderCells(cells, renderer, xOffset, yOffset, zoomFactor);
    lock.unlock();
    // Use fewer CPU cycles if paused
    if (paused) {
      SDL_Delay(25);
    }
    else {
      SDL_Delay(10);
    }
  }
  updateThread.join();
  SDL_DestroyWindow(window);
  SDL_Quit();
  return 0;
}
