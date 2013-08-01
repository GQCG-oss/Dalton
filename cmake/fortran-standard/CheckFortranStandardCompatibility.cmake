include(CheckFortranSourceCompiles)

# test whether we are able to compile a simple Fortran 2003 program
file(READ "${CMAKE_SOURCE_DIR}/cmake/fortran-standard/test-fortran-2003-compatibility.F90" _source)
check_fortran_source_compiles(
    ${_source}
    COMPILER_UNDERSTANDS_FORTRAN_2003
    )

if(COMPILER_UNDERSTANDS_FORTRAN_2003)
    add_definitions(-DCOMPILER_UNDERSTANDS_FORTRAN_2003)
else()
    message(FATAL_ERROR "Your compiler is not compatible with some Fortran 2003 features that we need.")
endif()
