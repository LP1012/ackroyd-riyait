# CMakeLists for RCANA software. Parts of this file are inspired by Cyclus
# (https://github.com/cyclus/cyclus).
#

cmake_minimum_required(VERSION 3.24)

project(
  ackroyd-riyait
  VERSION 0.1.0
  LANGUAGES CXX)

# ----------------------------------------------------------------------------#
# CMake load modules
# ----------------------------------------------------------------------------#
list(APPEND CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
     "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
message("-- CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")

include(Policies)
require_out_of_source_build()

include(GlobalVariables)


# ----------------------------------------------------------------------------#
# CMake configuration
# ----------------------------------------------------------------------------#
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

add_subdirectory(src)

# ----------------------------------------------------------------------------#
# Testing
# ----------------------------------------------------------------------------#
# must be at top-level cmake file OPTION(USE_TESTING "Build testing" ON)
option(ENABLE_TESTS "Enable tests" ON)
if(ENABLE_TESTS)
  include(CTest)
  enable_testing()

  add_subdirectory(tests)
endif()

