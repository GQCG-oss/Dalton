# radovan: this is so far not used but
#          soon it will be

add_definitions(-DBUILD_GEN1INT)
add_definitions(-DPRG_DIRAC)

add_library(
    gen1int_interface
    dirac_gen1int_interface.F90
    gen1int_api.F90
    gen1int_cube.F90
    gen1int_host.F90
    gen1int_matrix.F90
    gen1int_shell.F90
    )

# set(PARENT_DEFINITIONS "-DBUILD_GEN1INT -DBUILD_OPENRSP -DVAR_DALTON")
set(ExternalProjectCMakeArgs
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DPARENT_INCLUDE_DIR=${PROJECT_SOURCE_DIR}/src/include
    -DPARENT_MODULE_DIR=${PROJECT_BINARY_DIR}/modules
    -DPARENT_DEFINITIONS=${PARENT_DEFINITIONS}
    -DENABLE_64BIT_INTEGERS=${ENABLE_64BIT_INTEGERS}
    )

add_external(gen1int)

add_dependencies(gen1int_interface gen1int)
