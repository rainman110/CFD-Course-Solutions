cmake_minimum_required(VERSION 3.0)
project(Exercise1)

add_subdirectory(lib)

file(COPY ${PROJECT_SOURCE_DIR}/data DESTINATION ${PROJECT_BINARY_DIR})
add_subdirectory(tests)
