add_library(
    cgto_diff_eri_interface_host
    cgto_diff_eri_host.F90
    )

set(PARENT_DEFINITIONS "-DBUILD_CGTODIFF -DVAR_LSDALTON")

set(ExternalProjectCMakeArgs
    -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DPARENT_INCLUDE_DIR=${PROJECT_SOURCE_DIR}/include
    -DPARENT_MODULE_DIR=${PROJECT_BINARY_DIR}/modules
    -DPARENT_DEFINITIONS=${PARENT_DEFINITIONS}
    )

add_dependencies(cgto_diff_eri_interface_host cgto-diff-eri)
