include(CheckFortranSourceCompiles)

# test whether we are able to compile a simple Fortran 2003 program
file(READ "${CMAKE_SOURCE_DIR}/cmake/fortran-standard/test-fortran-2003-compatibility.F90" _source)
check_fortran_source_compiles(
    ${_source}
    COMPILER_UNDERSTANDS_FORTRAN03
    )

if(COMPILER_UNDERSTANDS_FORTRAN03)
    add_definitions(-DCOMPILER_UNDERSTANDS_FORTRAN_2003)
    # radovan: it is better if code stops at runtim iff
    #          user asks for the functionality that requires -DCOMPILER_UNDERSTANDS_FORTRAN_2003
#else()
#    message(FATAL_ERROR "Your compiler is not compatible with some Fortran 2003 features that we need.")
endif()

# test whether we are able to compile a simple Fortran 2003 program with ptr reshapes
file(READ "${CMAKE_SOURCE_DIR}/cmake/fortran-standard/test-fortran-pointer-reshape.F90" _source)
check_fortran_source_compiles(
    ${_source}
    PTR_RESHAPE_WORKS
    )

if(PTR_RESHAPE_WORKS)
    add_definitions(-DVAR_PTR_RESHAPE)
endif()

