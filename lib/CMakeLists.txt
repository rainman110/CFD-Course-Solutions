cmake_minimum_required(VERSION 3.0)
project(cfdcourselib)


add_library(cfd
  cfdlib/fields.hpp
  cfdlib/functions.hpp
  cfdlib/blas.hpp
  cfdlib/vtk-io.hpp
  cfdlib/vtk-io.cpp
  cfdlib/poisson.cpp
  cfdlib/poisson.hpp
  cfdlib/real.hpp
  cfdlib/solver.hpp
  cfdlib/timer.cpp
  cfdlib/timer.hpp
)

target_include_directories(cfd PUBLIC ${PROJECT_SOURCE_DIR})
