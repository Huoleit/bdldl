cmake_minimum_required(VERSION 3.16.0)
project(sparseLDL VERSION 0.1.0 LANGUAGES C ASM)

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
# sparseLDL configuration
# ###########################
option(SPARSE_LDL_USE_SINGLE_PRECISION "Use single precision floating point number" TRUE)

if(${SPARSE_LDL_USE_SINGLE_PRECISION})
  message(STATUS "Using single precision floats")
else()
  message(STATUS "Using double precision floats")
endif()

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/configuration/sparseLDL_types.h.in
  ${CMAKE_CURRENT_SOURCE_DIR}/include/${PROJECT_NAME}/sparseLDL_types.h
)

# ############################
# Build sparseLDL
# ############################
add_library(sparseLDL
  src/sparseLDL.c
)
target_include_directories(sparseLDL
  PUBLIC
  include
)
target_link_libraries(sparseLDL
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
