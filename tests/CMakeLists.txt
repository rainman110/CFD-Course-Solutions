cmake_minimum_required(VERSION 3.4)
project(cfdlib-tests)

if (NOT TARGET "cfd")
  find_package(cfdlib REQUIRED)
endif()

enable_testing()
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include(AddGoogleTest)

set(TESTSOURCES
  testmain.cpp
  ex1-field1d.cpp
  ex1-field2d.cpp
  ex1-integration.cpp
  ex2_integration.cpp
)

add_executable(runUnitTests ${TESTSOURCES})
target_link_libraries(runUnitTests PRIVATE cfd gtest_main)
add_gtest(runUnitTests)
