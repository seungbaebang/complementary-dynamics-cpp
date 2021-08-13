cmake_minimum_required(VERSION 3.1)
project(complementary_dynamics)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
#list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)
set(CMAKE_CXX_STANDARD 17)

set(bartels_ROOT Bartels)
set(bartels_INCLUDE_DIR ${bartels_ROOT}/include/)
set(bartels_LIB_DIR ${bartels_ROOT}/build/)

set(CMAKE_BUILD_TYPE Release)

# libigl
option(LIBIGL_USE_PREBUILT_LIBRARIES "Use prebuilt libraries"       ON)
option(LIBIGL_USE_STATIC_LIBRARY     "Use libigl as static library" ON)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
# option(LIBIGL_WITH_MATLAB            "Use Matlab"         ON)

find_package(LIBIGL REQUIRED QUIET)

# include(${bartels_ROOT}/cmake/bartels.cmake)

# Add your project files
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include/")
#include_directories(${bartels_INCLUDE_DIR})
include_directories("${bartels_INCLUDE_DIR}")
link_directories("${bartels_LIB_DIR}")

file(GLOB SRCFILES src/*.cpp include/*.h)

# add_library(compdyna ${SRCFILES})
#add_library(bartels ${SRCFILES})

# target_link_libraries(core igl::core igl::opengl igl::opengl_glfw)

# add_executable(lbs ${SRCFILES} demo/lbs.cpp)
# target_link_libraries(lbs igl::core igl::opengl igl::opengl_glfw bartels)

add_executable(complementary_dynamics ${SRCFILES} demo/complementary_dynamics.cpp)
target_link_libraries(complementary_dynamics igl::core igl::opengl igl::opengl_glfw bartels)



