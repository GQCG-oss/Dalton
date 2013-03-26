add_library(
    xc_interface_host
    xc_host.F90
    )

set(PARENT_DEFINITIONS "-DBUILD_XCINT -DVAR_LSDALTON")

set(ExternalProjectCMakeArgs
    -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DPARENT_INCLUDE_DIR=${PROJECT_SOURCE_DIR}/include
    -DPARENT_MODULE_DIR=${PROJECT_BINARY_DIR}/modules
    -DPARENT_DEFINITIONS=${PARENT_DEFINITIONS}
    )

add_dependencies(xc_interface_host xcint)
