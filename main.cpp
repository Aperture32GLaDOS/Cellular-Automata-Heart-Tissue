#include "cells.cpp"
#include <SDL2/SDL.h>
#include <SDL2/SDL_events.h>
#include <SDL2/SDL_keycode.h>
#include <SDL2/SDL_mouse.h>
#include <SDL2/SDL_syswm.h>
#include <SDL2/SDL_video.h>
#include <SDL2/SDL_ttf.h>
#include <chrono>
#include <cstdlib>
#include <fftw3.h>
#include <thread>
#include <mutex>

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
void updateCells(Cells* cells, bool* quit, bool* paused, bool* step, int* frameTime, double* stateArray) {
  long int startTime;
  long int elapsedTime;
  double* distanceCoefficientsPadded = fftw_alloc_real(cells->height * cells->width);
  // Complex arrays roughly half the size since FFTW takes advantage of complex conjugates
  fftw_complex* distanceCoefficientsTransformed = fftw_alloc_complex(cells->height * (cells->width / 2 + 1));
  fftw_complex* stateArrayTransformed = fftw_alloc_complex(cells->height * (cells->width / 2 + 1));
  double* distanceArray = fftw_alloc_real(cells->height * cells->width);
  // Plan FFTs before array initialization
  fftw_plan distanceCoefficientsFFT = fftw_plan_dft_r2c_2d(cells->height, cells->width, distanceCoefficientsPadded, distanceCoefficientsTransformed, 0);
  fftw_plan stateArrayFFT = fftw_plan_dft_r2c_2d(cells->height, cells->width, stateArray, stateArrayTransformed, 0);
  fftw_plan stateArrayIFFT = fftw_plan_dft_c2r_2d(cells->height, cells->width, stateArrayTransformed, distanceArray, 0);
  double* distanceCoefficients = new double[SEARCH_RADIUS * SEARCH_RADIUS];
  int offsetLength = 0;
  for (int i = 0; i < SEARCH_RADIUS; i++) {
    for (int j = 0; j < SEARCH_RADIUS; j++) {
      if (i == SEARCH_RADIUS / 2 && j == SEARCH_RADIUS / 2) continue;
      double distance = (i - SEARCH_RADIUS / 2) * (i - SEARCH_RADIUS / 2) + (j - SEARCH_RADIUS / 2) * (j - SEARCH_RADIUS / 2);
      distanceCoefficients[(i * SEARCH_RADIUS) + j] = 1.0 / distance;
    }
  }
  shiftConvolution(distanceCoefficients, distanceCoefficientsPadded, SEARCH_RADIUS, cells->height, cells->width);
  delete[] distanceCoefficients;
  fftw_execute(distanceCoefficientsFFT);
  // The coeffients never change, and therefore only need to be transformed once
  fftw_destroy_plan(distanceCoefficientsFFT);
  while (!(*quit)) {
    startTime = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    // Lock the mutex, as data is being written
    std::unique_lock<std::mutex> lock(mu);
    advanceCells(cells, distanceCoefficientsTransformed, stateArray, stateArrayTransformed, distanceArray, stateArrayFFT, stateArrayIFFT);
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
        advanceCells(cells, distanceCoefficientsTransformed, stateArray, stateArrayTransformed, distanceArray, stateArrayFFT, stateArrayIFFT);
        lock.unlock();
      }
    }
  }
  // Cleanup
  fftw_destroy_plan(stateArrayFFT);
  fftw_destroy_plan(stateArrayIFFT);
  fftw_free(distanceCoefficientsPadded);
  fftw_free(distanceCoefficientsTransformed);
  fftw_free(stateArrayTransformed);
  fftw_free(distanceArray);
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
  // Initialize all cells to be inactive normal tissue
  for (int i = 0; i < cells.height; i++) {
    for (int j = 0; j < cells.width; j++) {
      Cell newCell;
      newCell.type = CellType::Tissue;
      newCell.state = 0;
      cells.cells[i * cells.width + j] = newCell;
    }
  }
  for (int i = -3; i <= 3; i++) {
    for (int j = -3; j <= 3; j++) {
      cells.cells[(cells.height / 2 + i) * cells.width + (cells.width / 2 + j)].type = CellType::Pacemaker;
    }
  }
  cells.cells[(cells.height / 2) * cells.width + (cells.width / 2)].state = AP_DURATION;
  window = SDL_CreateWindow("Heart Tissue", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
      cells.width, cells.height, SDL_WINDOW_SHOWN);
  renderer = SDL_CreateRenderer(window, -1, 0);
  // TODO: automatic file location OR have a font folder in the project
  TTF_Font* font = TTF_OpenFont("/usr/share/fonts/Ubuntu.ttf", 32);
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
  int frameTime = 500;
  double* stateArray = fftw_alloc_real(cells.height * cells.width);
  std::thread updateThread(updateCells, &cells, &quit, &paused, &step, &frameTime, stateArray);
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
      // Keep the mousePosX and mousePosY updated
      else if (currentEvent.type == SDL_MOUSEMOTION) {
        SDL_GetMouseState(&mousePosX, &mousePosY);
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
        Cell* selectedCell = &cells.cells[((int) ((mousePosY / zoomFactor) - yOffset)) * cells.width + (int) ((mousePosX / zoomFactor) - xOffset)];
        if (currentEvent.button.button == SDL_BUTTON_LEFT) {
          selectedCell->state = AP_DURATION;
          if (selectedCell->type != CellType::RestingTissue) {
            stateArray[((int) ((mousePosY / zoomFactor) - yOffset)) * cells.width + (int) ((mousePosX / zoomFactor) - xOffset)] = selectedCell->state;
          }
        }
        else if (currentEvent.button.button == SDL_BUTTON_RIGHT) {
          if (selectedCell-> state != 0 && selectedCell->type != CellType::RestingTissue) selectedCell->state = 0;
          if (selectedCell->type != CellType::RestingTissue) {
            stateArray[((int) ((mousePosY / zoomFactor) - yOffset)) * cells.width + (int) ((mousePosX / zoomFactor) - xOffset)] = selectedCell->state;
          }
        }
        // TODO: change tissue type on shift-right click (or similar)
        lock.unlock();
      }
    }
    std::unique_lock<std::mutex> lock(mu);
    renderCells(cells, renderer, font, xOffset, yOffset, zoomFactor, ((int) ((mousePosY / zoomFactor) - yOffset)), (int) ((mousePosX / zoomFactor) - xOffset));
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
