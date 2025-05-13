/*
This program is a cellular automata which models heart tissue - all source files are to be
licensed under the conditions defined in LICENSE.md
Copyright (C) 2025 Eshe Hinchliffe
*/

#include "cells.cpp"
#include <SDL2/SDL.h>
#include <SDL2/SDL_events.h>
#include <SDL2/SDL_keycode.h>
#include <SDL2/SDL_mouse.h>
#include <SDL2/SDL_render.h>
#include <SDL2/SDL_syswm.h>
#include <SDL2/SDL_video.h>
#include <SDL2/SDL_ttf.h>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fftw3.h>
#include <thread>
#include <mutex>
#include <utility>

SDL_Window* window = NULL;
SDL_Renderer* renderer = NULL;


std::mutex mu;

// Cyclically shifts a convolution to expand it, i.e. for padding reasons
void shiftConvolution(double* originalConvolution, double* shiftedConvolution, int convWidth, int dataHeight, int dataLength) {
  //Top left corner (shifted so that the middle element of the convolution is now at (0, 0))
  for (int i = convWidth / 2; i < convWidth; i++) {
    for (int j = convWidth / 2; j < convWidth; j++) {
      shiftedConvolution[(i - (convWidth / 2)) * dataLength + (j - (convWidth / 2))] = originalConvolution[i * convWidth + j];
    }
  }
  // Top right corner
  for (int i = convWidth / 2; i < convWidth; i++) {
    for (int j = 0; j < convWidth / 2; j++) {
      shiftedConvolution[(i - (convWidth / 2)) * dataLength + (dataLength + j - convWidth / 2)] = originalConvolution[i * convWidth + j];
    }
  }
  // Bottom left corner
  for (int i = 0; i < convWidth / 2; i++) {
    for (int j = convWidth / 2; j < convWidth; j++) {
      shiftedConvolution[(dataHeight + i - convWidth / 2) * dataLength + (j - convWidth / 2)] = originalConvolution[i * convWidth + j];
    }
  }
  // Bottom right corner
  for (int i = 0; i < (convWidth / 2); i++) {
    for (int j = 0; j < (convWidth / 2); j++) {
      shiftedConvolution[(dataHeight + i - convWidth / 2) * dataLength + (dataLength + j - convWidth / 2)] = originalConvolution[i * convWidth + j];
    }
  }
}

// Updates the cells in a seperate thread, so as to keep the render updates fast
void updateCells(Cells* cells, bool* quit, bool* paused, bool* step, int* frameTime, double* stateArray, NeighbourCounter* neighbourCounter) {
  long int startTime;
  long int elapsedTime;
  while (!(*quit)) {
    startTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    // Lock the mutex, as data is being written
    std::unique_lock<std::mutex> lock(mu);
    advanceCells(cells, neighbourCounter);
    lock.unlock();
    elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count() - startTime;
    if (elapsedTime <= *frameTime) {
      std::this_thread::sleep_for(std::chrono::milliseconds(*frameTime - elapsedTime));
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(25));
    while (*paused && !(*quit)) {
      std::this_thread::sleep_for(std::chrono::milliseconds(250));
      if (*step) {
        *step = false;
        std::unique_lock<std::mutex> lock(mu);
        advanceCells(cells, neighbourCounter);
        lock.unlock();
      }
    }
  }
}


int main (int argc, char *argv[]) {
  if (SDL_Init(SDL_INIT_VIDEO) < 0 || TTF_Init() < 0) {
    return 0;
  }
  // Declare the 2D plane of cells
  Cells cells;
  cells.width = SIZE;
  cells.height = SIZE;
  cells.cells = new Cell[cells.height * cells.width];
  cells.numOrientations = 1;
  cells.orientations = new Orientation[1];
  cells.orientations[0].xDir = 1.0;
  cells.orientations[0].yDir = 0.0;
  cells.orientations[0].cellCount = cells.height * cells.width * 0.5;
  // Initialize all cells to be inactive normal tissue
  for (int i = 0; i < cells.height; i++) {
    for (int j = 0; j < cells.width; j++) {
      Cell newCell;
      newCell.type = CellType::Tissue;
      newCell.state = 0;
      if (j < cells.width / 2) {
        newCell.orientationIndex = 0;
        cells.orientations[0].cells.push_front(&newCell);
      }
      else {
        newCell.orientationIndex = 0;
        cells.orientations[0].cells.push_front(&newCell);
      }
      cells.cells[i * cells.width + j] = newCell;
    }
  }
  // for (int i = 0; i < 7; i++) {
  //   for (int j = 0; j < 7; j++) {
  //     cells.cells[((i - 3 + cells.height / 2) * cells.width) + (j - 3) + cells.width / 2].type = CellType::Pacemaker;
  //   }
  // }
  window = SDL_CreateWindow("Heart Tissue", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
      cells.width, cells.height, SDL_WINDOW_SHOWN);
  renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
  // TODO: automatic file location OR have a font folder in the project
  TTF_Font* font = TTF_OpenFont("/usr/share/fonts/TTF/FiraCode-Regular.ttf", 32);
  SDL_RenderPresent(renderer);
  bool quit = false;
  bool paused = true;
  bool step = false;
  SDL_Event currentEvent;
  float xOffset = 0;
  float yOffset = 0;
  float zoomFactor = 1;
  int mousePosX;
  int mousePosY;
  int selectedCellX;
  int selectedCellY;
  int frameTime = 500;
  bool isSelectingRect = false;
  bool isUsingRect = false;
  int firstCornerY;
  int firstCornerX;
  int secondCornerY;
  int secondCornerX;
  int highlightedX = -1;
  int highlightedY = -1;
  double* stateArray = fftw_alloc_real(cells.height * cells.width);
  for (int i = 0; i < cells.height * cells.width; i++) {
    stateArray[i] = 0.0;
  }
  NeighbourCounter neighbourCounter(&cells, stateArray);
  std::thread updateThread(updateCells, &cells, &quit, &paused, &step, &frameTime, stateArray, &neighbourCounter);
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
        else if (currentEvent.key.keysym.sym == SDLK_r) {
          if (!isSelectingRect) {
            firstCornerY = selectedCellY;
            firstCornerX = selectedCellX;
          }
          else {
            secondCornerY = selectedCellY;
            secondCornerX = selectedCellX;
          }
          isSelectingRect = !isSelectingRect;
        }
        else if (currentEvent.key.keysym.sym == SDLK_LSHIFT) {
          isUsingRect = true;
        }
        else if (currentEvent.key.keysym.sym == SDLK_h) {
          if (highlightedX == -1 && highlightedY == -1) {
            highlightedX = selectedCellX;
            highlightedY = selectedCellY;
          }
          else {
            float distance = std::sqrt(std::pow((highlightedX - selectedCellX), 2.0f) + std::pow(highlightedY - selectedCellY, 2.0f));
            std::cout << "Distance between cells: " << distance << std::endl;
            highlightedX = -1;
            highlightedY = -1;
          }
        }
        else if (currentEvent.key.keysym.sym == SDLK_SPACE) {
          paused = !paused;
        }
        else if (currentEvent.key.keysym.sym == SDLK_PERIOD) {
          step = true;
        }
        else if (currentEvent.key.keysym.sym == SDLK_EQUALS) {
          frameTime -= 50;
          if (frameTime < 0) {
            frameTime = 0;
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
          delete[] cells.orientations;
          cells = readCellsFromFile("cells.dmp");
          SDL_SetWindowSize(window, cells.width, cells.height);
          for (int i = 0; i < cells.height; i++) {
            for (int j = 0; j < cells.width; j++) {
              Cell currentCell = cells.cells[i * cells.width + j];
              if (currentCell.type == CellType::RestingTissue) {
                stateArray[i * cells.width + j] = 0;
              }
              else {
                stateArray[i * cells.width + j] = currentCell.state;
              }
            }
          }
          neighbourCounter.reinitialize();
          lock.unlock();
        }
        else if (currentEvent.key.keysym.sym == SDLK_MINUS) {
          frameTime += 50;
        }
        // Equivalent to giving a shock to the whole heart
        else if (currentEvent.key.keysym.sym == SDLK_g) {
          std::unique_lock<std::mutex> lock(mu);
          for (int i = 0; i < cells.height * cells.width; i++) {
            Cell selectedCell = cells.cells[i];
            if (selectedCell.type == CellType::RestingTissue) {
              continue;
            }
            selectedCell.state = AP_DURATION;
            cells.cells[i] = selectedCell;
            stateArray[i] = 0.0;
          }
          lock.unlock();
        }
      }
      else if (currentEvent.type == SDL_KEYUP) {
        if (currentEvent.key.keysym.sym == SDLK_LSHIFT) {
          isUsingRect = false;
        }
      }
      // Keep the mousePosX and mousePosY updated
      else if (currentEvent.type == SDL_MOUSEMOTION) {
        SDL_GetMouseState(&mousePosX, &mousePosY);
        selectedCellY = (int) ((mousePosY / zoomFactor) - yOffset) % cells.height;
        selectedCellX = (int) ((mousePosX / zoomFactor) - xOffset) % cells.width;
        if ((mousePosX / zoomFactor) - xOffset < 0) {
          selectedCellX--;
        }
      }
      else if (currentEvent.type == SDL_MOUSEWHEEL) {
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
        if (!isUsingRect) {
          Cell* selectedCell = &cells.cells[selectedCellY * cells.width + selectedCellX];
          if (currentEvent.button.button == SDL_BUTTON_LEFT) {
            selectedCell->state = AP_DURATION;
            if (selectedCell->type != CellType::RestingTissue) {
              stateArray[selectedCellY * cells.width + selectedCellX] = selectedCell->state;
            }
          }
          else if (currentEvent.button.button == SDL_BUTTON_RIGHT) {
            if (selectedCell-> state != 0 && selectedCell->type != CellType::RestingTissue) selectedCell->state = 0;
            if (selectedCell->type != CellType::RestingTissue) {
              stateArray[selectedCellY * cells.width + selectedCellX] = selectedCell->state;
            }
          }
          else if (currentEvent.button.button == SDL_BUTTON_MIDDLE) {
            if (selectedCell->type == CellType::RestingTissue) {
              selectedCell->type = CellType::Tissue;
              stateArray[selectedCellY * cells.width + selectedCellX] = selectedCell->state;
            }
            else if (selectedCell->type == CellType::Tissue) {
              selectedCell->type = CellType::RestingTissue;
              stateArray[selectedCellY * cells.width + selectedCellX] = 0;
            }
          }
        }
        else {
          if (firstCornerY > secondCornerY) {
            std::swap(firstCornerY, secondCornerY);
          }
          if (firstCornerX > secondCornerX) {
            std::swap(firstCornerX, secondCornerX);
          }
          for (int i = firstCornerY; i < secondCornerY; i++) {
            for (int j = firstCornerX; j < secondCornerX; j++) {
              Cell* selectedCell = &cells.cells[i * cells.width + j];
              if (currentEvent.button.button == SDL_BUTTON_LEFT) {
                selectedCell->state = AP_DURATION;
                if (selectedCell->type != CellType::RestingTissue) {
                  stateArray[i * cells.width + j] = selectedCell->state;
                }
              }
              else if (currentEvent.button.button == SDL_BUTTON_RIGHT) {
                if (selectedCell-> state != 0 && selectedCell->type != CellType::RestingTissue) selectedCell->state = 0;
                if (selectedCell->type != CellType::RestingTissue) {
                  stateArray[i * cells.width + j] = selectedCell->state;
                }
              }
              else if (currentEvent.button.button == SDL_BUTTON_MIDDLE) {
                if (selectedCell->type == CellType::RestingTissue) {
                  selectedCell->type = CellType::Tissue;
                  stateArray[i * cells.width + j] = selectedCell->state;
                }
                else if (selectedCell->type == CellType::Tissue) {
                  selectedCell->type = CellType::RestingTissue;
                  stateArray[i * cells.width + j] = 0;
                }
              }
            }
          }
        }
        // TODO: change tissue type on shift-right click (or similar)
        lock.unlock();
      }
    }
    std::unique_lock<std::mutex> lock(mu);
    renderCells(cells, renderer, font, xOffset, yOffset, zoomFactor, selectedCellY, selectedCellX, firstCornerY, secondCornerY, firstCornerX, secondCornerX);
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
  delete[] cells.cells;
  fftw_free(stateArray);
  TTF_CloseFont(font);
  SDL_DestroyWindow(window);
  SDL_Quit();
  return 0;
}
