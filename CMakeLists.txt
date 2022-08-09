cmake_minimum_required(VERSION 3.16.0)
project(sparseLDL VERSION 0.1.0 LANGUAGES C ASM)
add_subdirectory(blasfeo)

option(SPARSE_LDL_USE_SINGLE_PRECISION "Use single precision floating point number" TRUE)

if(${SPARSE_LDL_USE_SINGLE_PRECISION})
  message(STATUS "Using single precision floats")
else()
  message(STATUS "Using double precision floats")
endif()

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/configuration/ldl_types.h.in
  ${CMAKE_CURRENT_SOURCE_DIR}/include/${PROJECT_NAME}/ldl_types.h)

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