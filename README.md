# Cellular automata model of heart tissue
This program models heart tissue using a cellular automata. Some features include
- Read/write capability
- Multithreading
- SIMD speedup
## Compiling and using
To compile simply create a build directory (i.e. mkdir build), and from there, run cmake .., followed by make. This project only requires SDL2, and FFTW.
To use the program, you can use the scroll wheel to zoom in and out of the simulation, and WASD to pan around. Left and right click changes the states of cells, and pressing "R" defines a rectangle. Pressing "G" stimulates the entire simulation at once, and clicking with the middle mouse button will toggle cells between an active and resting state.
Hold shift to apply the operators on the area defined in the rectangle.
