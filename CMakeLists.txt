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


FILE(GLOB_RECURSE additional_source_files ${CMAKE_CURRENT_SOURCE_DIR}/cwrapper/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/Others/SparseBLAS/*.cpp)
FILE(GLOB source_files ${CMAKE_CURRENT_SOURCE_DIR}/Manifolds/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/Manifolds/Stiefel/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/Manifolds/Sphere/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/Problems/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/Problems/estimateRotation/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/Problems/estimateRotation_Temporal/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/Problems/SphereConvexHull/*.cpp  )
SET(solvers ${CMAKE_CURRENT_SOURCE_DIR}/Solvers/RTRNewton.cpp ${CMAKE_CURRENT_SOURCE_DIR}/Solvers/SolversTR.cpp ${CMAKE_CURRENT_SOURCE_DIR}/Solvers/QuasiNewton.cpp ${CMAKE_CURRENT_SOURCE_DIR}/Solvers/Solvers.cpp )
SET(additional_debug_files ${CMAKE_CURRENT_SOURCE_DIR}/Others/ForDebug.cpp ${CMAKE_CURRENT_SOURCE_DIR}/Others/randgen.cpp ${CMAKE_CURRENT_SOURCE_DIR}/Others/Timer.cpp ${CMAKE_CURRENT_SOURCE_DIR}/Others/MyMatrix.cpp ${CMAKE_CURRENT_SOURCE_DIR}/Others/MinPNormConHull.cpp)
# Libraries
add_library(ropt STATIC ${source_files} ${solvers} ${additional_source_files} ${additional_debug_files})

# Mandatory for running the test example file
add_definitions(-DTESTSIMPLEEXAMPLE)

# Executable
add_executable(TestSimpleExample ${CMAKE_SOURCE_DIR}/test/test1.cpp )

# Link
target_link_libraries(TestSimpleExample ropt lapack blas ${OpenCV_LIBS} )

