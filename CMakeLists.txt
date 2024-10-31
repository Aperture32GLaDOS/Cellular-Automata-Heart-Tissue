cmake_minimum_required(VERSION 3.5)
set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra")


project(Cellular_Automata_Heart_Tissue
  VERSION 1.0
  DESCRIPTION "A cellular automata which simulates heart tissue"
  LANGUAGES CXX)

add_executable(main main.cpp)
find_package(SDL2 REQUIRED)
target_include_directories(main PRIVATE ${SDL2_INCLUDE_DIRS})
target_link_libraries(main ${SDL2_LIBRARIES})
