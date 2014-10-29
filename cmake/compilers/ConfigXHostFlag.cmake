
# Description
#     We use a Fortran test program
#     to find out whether we can compile and run with -xHost.

# Authors
#     Ivan Hrasko
#     Miroslav Ilias
#     Radovan Bast

# Variables used
#     CMAKE_Fortran_COMPILER_WORKS
#     CMAKE_Fortran_COMPILER_ID
#     CMAKE_SYSTEM_NAME

# Variables set
#     XHOST_FLAG_AVAILABLE

#-------------------------------------------------------------------------------

# default to OFF
set(XHOST_FLAG_AVAILABLE OFF)

if(CMAKE_Fortran_COMPILER_WORKS AND CMAKE_Fortran_COMPILER_ID MATCHES "Intel")

    # compile test program
    if(CMAKE_SYSTEM_NAME MATCHES "Windows")
        set(_xhost_binary ${CMAKE_BINARY_DIR}/test_xhost_flag.exe)
        execute_process(COMMAND ifort /QxHost ${CMAKE_SOURCE_DIR}/cmake/compilers/test_xhost_flag.f90 /exe:${_xhost_binary})
    else()
        set(_xhost_binary ${CMAKE_BINARY_DIR}/test_xhost_flag.x)
        execute_process(COMMAND ifort -xHost ${CMAKE_SOURCE_DIR}/cmake/compilers/test_xhost_flag.f90 -o ${_xhost_binary})
    endif()

    # run test program
    execute_process(
        COMMAND ${_xhost_binary}
        OUTPUT_VARIABLE _xhost_binary_output
        ERROR_VARIABLE _xhost_binary_output
    )

    # depending on result set XHOST_FLAG_AVAILABLE
    if(_xhost_binary_output MATCHES "Intel xHost flag is available on this computer.")
        set(XHOST_FLAG_AVAILABLE ON)
    endif()

endif()
