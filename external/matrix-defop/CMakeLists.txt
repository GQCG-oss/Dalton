cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(matrix-defop Fortran C)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/compilers
    )

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)
include(ConfigParentSettings)

set(SOURCES
    ${PROJECT_SOURCE_DIR}/src/matrix_backend.F90
    ${PROJECT_SOURCE_DIR}/src/matrix_defop.F90
    ${PROJECT_SOURCE_DIR}/src/matrix_lowlevel.F90
    )

include(ConfigCompilerFlags)

add_library(
    matrix-defop
    ${SOURCES}
    )

install(TARGETS matrix-defop ARCHIVE DESTINATION lib)
