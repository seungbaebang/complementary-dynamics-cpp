cmake_minimum_required(VERSION 3.1)
project(complementary_dynamics)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
#list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
set(CMAKE_CXX_STANDARD 14)

set(bartels_ROOT ../Bartels)
set(bartels_INCLUDE_DIR ${bartels_ROOT}/include/)
set(bartels_LIB_DIR ${bartels_ROOT}/build/)

message("bartels_INCLUDE_DIR: ${bartels_INCLUDE_DIR}")
# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_MATLAB            "Use Matlab"         ON)
#option(LIBIGL_WITH_CGAL              "Use CGAL"           ON)
option(LIBIGL_USE_STATIC_LIBRARY     "Use libIGL as static librarie" ON)

find_package(LIBIGL REQUIRED QUIET)

#include(${bartels_SOURCE_DIR}/cmake/bartels.cmake)

# Add your project files
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include/")
#include_directories(${bartels_INCLUDE_DIR})
include_directories("${bartels_INCLUDE_DIR}")
link_directories("${bartels_LIB_DIR}")

file(GLOB SRCFILES *.cpp src/*.cpp)

add_library(core ${SRCFILES})
#add_library(bartels ${SRCFILES})

target_link_libraries(core igl::core igl::opengl_glfw igl::matlab)

add_executable(linear_demo linear_demo.cpp)
target_link_libraries(linear_demo core bartels)

add_executable(nonlinear_demo nonlinear_demo.cpp)
target_link_libraries(nonlinear_demo core bartels)

