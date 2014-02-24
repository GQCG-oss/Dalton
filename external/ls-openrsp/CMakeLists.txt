cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(openrsp Fortran C)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/compilers
    )

option(ENABLE_64BIT_INTEGERS "Enable 64-bit integers" OFF)

include(ConfigSafeGuards)

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)
include(ConfigParentSettings)

set(SOURCES
    ${PROJECT_SOURCE_DIR}/src/rsp_contribs.F90
    ${PROJECT_SOURCE_DIR}/src/rsp_equations.F90
    )

include(ConfigCompilerFlags)

add_library(
    openrsp
    ${SOURCES}
    )

install(TARGETS openrsp ARCHIVE DESTINATION lib)
