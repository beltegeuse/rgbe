# Tested only on Windows and Linux (Ubuntu 10.10)
cmake_minimum_required(VERSION 2.8.3 FATAL_ERROR)
project(rgbe)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/data/cmake")

# Add openmp for the parallel computation
# This package is optional but it is highly advised to have it!
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Add OpenEXR support
# As this support depends on external library, this support is optional
find_package (OpenEXR 2.0 REQUIRED)
#OpenEXR 2.2 still has problems with importing ImathInt64.h unqualified
#thus need for ilmbase/OpenEXR
include_directories ("${OPENEXR_INCLUDE_DIR}"
        "${ILMBASE_INCLUDE_DIR}"
        "${ILMBASE_INCLUDE_DIR}/OpenEXR")
if (${OPENEXR_VERSION} VERSION_LESS 2.0.0)
    message (FATAL_ERROR "OpenEXR/Ilmbase is too old")
endif ()
if (NOT OpenEXR_FIND_QUIETLY)
    message (STATUS "OPENEXR_INCLUDE_DIR = ${OPENEXR_INCLUDE_DIR}")
    message (STATUS "OPENEXR_LIBRARIES = ${OPENEXR_LIBRARIES}")
endif ()

add_subdirectory(pybind11)
pybind11_add_module(rgbe sources/rgbe.cpp)
target_link_libraries(rgbe PRIVATE ${OPENEXR_LIBRARIES})
  
  
