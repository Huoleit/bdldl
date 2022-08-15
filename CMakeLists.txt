cmake_minimum_required(VERSION 3.16.0)
project(bdldl VERSION 0.1.0 LANGUAGES C ASM)

# ###########################
# Fetch Dependencies
# ###########################
include(FetchContent)

# BLASFEO Settings
set(BUILD_SHARED_LIBS ON CACHE STRING "Build shared libraries" FORCE)
set(TARGET X64_AUTOMATIC CACHE STRING "Target architecture" FORCE)
set(BLASFEO_EXAMPLES OFF CACHE BOOL "Examples enabled")

FetchContent_Declare(blasfeoPackage
  GIT_REPOSITORY https://github.com/giaf/blasfeo
  GIT_TAG 466812522c8b3840b845f74a5da2e11b35bdcc22 # release-0.1.3
)
FetchContent_MakeAvailable(blasfeoPackage)

# ###########################
# bdldl configuration
# ###########################
option(BDLDL_USE_SINGLE_PRECISION "Use single precision floating point number" OFF)
message(STATUS "Using single precision floats - ${BDLDL_USE_SINGLE_PRECISION}")

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/configuration/bdldl_types.h.in
  ${CMAKE_CURRENT_SOURCE_DIR}/include/${PROJECT_NAME}/bdldl_types.h
)

# ############################
# Build bdldl
# ############################
add_library(bdldl
  src/bdldl.c
)
target_include_directories(bdldl
  PUBLIC
  include
)
target_link_libraries(bdldl
  blasfeo
)

# #############################
# Documentation
# #############################
find_package(Doxygen)
set(DOXYGEN_OUTPUT_DIRECTORY docs)

doxygen_add_docs(doxygen
  include/ src/
)
