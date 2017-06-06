#Set minimum version requered
cmake_minimum_required(VERSION 2.4.6)

#set project name
project(ROPTLIB)

find_package(OpenCV REQUIRED)

#set the include directories
#include_directories (${CMAKE_CURRENT_SOURCE_DIR})

# Use C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11 -O3 -march=native")

include_directories(${OpenCV_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/cwrapper/blas ${CMAKE_CURRENT_SOURCE_DIR}/cwrapper/lapack )

# Find all the cpp files
FILE(GLOB_RECURSE source_files ${CMAKE_CURRENT_SOURCE_DIR}/Manifolds/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/Others/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/Problems/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/Solvers/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/cwrapper/*.cpp)

# Libraries
add_library(ropt STATIC ${source_files})

# Mandatory for running the test example file
add_definitions(-DTESTSIMPLEEXAMPLE)

# Executable
add_executable(TestSimpleExample ${CMAKE_SOURCE_DIR}/test/test1.cpp )

# Link
target_link_libraries(TestSimpleExample ropt lapack blas ${OpenCV_LIBS} )

