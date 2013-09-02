add_library(
    gen1int_interface
    gen1int_shell.F90
    gen1int_api.F90
    gen1int_host.F90
    )

set(PARENT_DEFINITIONS "-DBUILD_GEN1INT_LSDALTON -DVAR_LSDALTON")

set(ExternalProjectCMakeArgs
    -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DPARENT_MODULE_DIR=${PROJECT_BINARY_DIR}/modules
    -DPARENT_DEFINITIONS=${PARENT_DEFINITIONS}
    )

add_external(gen1int)

#add_dependencies(gen1int_interface matrix-defop)
add_dependencies(gen1int_interface gen1int)
