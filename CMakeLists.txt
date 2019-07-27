cmake_minimum_required(VERSION 3.1)
project(Exercise1)

set (CMAKE_CXX_STANDARD 11)

add_subdirectory(lib)

file(COPY ${PROJECT_SOURCE_DIR}/data DESTINATION ${PROJECT_BINARY_DIR})
add_subdirectory(tests)
