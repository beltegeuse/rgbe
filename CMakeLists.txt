# Tested only on Windows and Linux (Ubuntu 10.10)
cmake_minimum_required(VERSION 2.8.3 FATAL_ERROR)
project(rgbe)

# Add openmp for the //
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_subdirectory(pybind11)
pybind11_add_module(rgbe sources/rgbe.cpp)

  
  
