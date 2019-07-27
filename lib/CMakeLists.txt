cmake_minimum_required(VERSION 3.0)
project(cfdcourselib)


add_library(cfd
  cfdlib/fields.hpp
  cfdlib/functions.hpp
  cfdlib/blas.hpp
  cfdlib/vtk-io.hpp
  cfdlib/vtk-io.cpp
)

target_include_directories(cfd PUBLIC ${PROJECT_SOURCE_DIR})
