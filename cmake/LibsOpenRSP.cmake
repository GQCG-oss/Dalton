set(PARENT_DEFINITIONS "-DBUILD_GEN1INT -DBUILD_OPENRSP -DPRG_DALTON")
if(MPI_FOUND)
    set(PARENT_DEFINITIONS "${PARENT_DEFINITIONS} -DVAR_MPI")
endif()
set(ExternalProjectCMakeArgs
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DPARENT_INCLUDE_DIR=${CMAKE_SOURCE_DIR}/DALTON/include
    -DPARENT_MODULE_DIR=${PROJECT_BINARY_DIR}/modules
    -DPARENT_DEFINITIONS=${PARENT_DEFINITIONS}
    )

add_external(xcint)
add_external(matrix-defop)
add_external(openrsp)
add_external(cgto-diff-eri)

add_dependencies(cgto-diff-eri     matrix-defop)
add_dependencies(dalton            cgto-diff-eri)
add_dependencies(dalton            gen1int_interface)
add_dependencies(dalton            matrix-defop)
add_dependencies(dalton            openrsp)
add_dependencies(dalton            xcint)
add_dependencies(gen1int           matrix-defop)
add_dependencies(gen1int_interface matrix-defop)
add_dependencies(openrsp           cgto-diff-eri)
add_dependencies(openrsp           gen1int_interface)
add_dependencies(openrsp           matrix-defop)
add_dependencies(openrsp           xcint)

set(EXTERNAL_LIBS
    ${PROJECT_BINARY_DIR}/external/lib/libopenrsp.a
    dalton
    gen1int_interface
    ${PROJECT_BINARY_DIR}/external/lib/libcgto-diff-eri.a
    ${PROJECT_BINARY_DIR}/external/lib/libgen1int.a
    ${PROJECT_BINARY_DIR}/external/lib/libxcint.a
    ${PROJECT_BINARY_DIR}/external/lib/librhoType.a
    ${PROJECT_BINARY_DIR}/external/lib/libmatrix-defop.a
    ${PROJECT_BINARY_DIR}/external/xcint-build/external/lib/libxcfun_f90_bindings.a
    ${PROJECT_BINARY_DIR}/external/xcint-build/external/lib/libxcfun.a
    stdc++
    ${EXTERNAL_LIBS}
    )
