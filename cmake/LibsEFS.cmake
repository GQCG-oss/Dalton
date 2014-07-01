set(_efs_definitions "-DPRG_DALTON")
if(MPI_FOUND)
    set(_efs_definitions "${_efs_definitions} -DVAR_MPI")
endif()

set(ExternalProjectCMakeArgs
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DPARENT_DEFINITIONS=${_efs_definitions}
    )

unset(_efs_definitions)

add_external(efs)

set(OMP_LIB)
if(CMAKE_COMPILER_IS_GNUCXX)
    set(OMP_LIB gomp)
elseif(CMAKE_CXX_COMPILER_ID MATCHES Intel)
    set(OMP_LIB iomp5)
endif()

set(DALTON_LIBS
    ${PROJECT_BINARY_DIR}/external/lib/libechidna.a
    ${PROJECT_BINARY_DIR}/external/lib/libquimera.a
    ${PROJECT_BINARY_DIR}/external/lib/libnewcontractions.a
    ${PROJECT_BINARY_DIR}/external/efs-build/libk2c_library.a
    stdc++
    ${OMP_LIB}
    ${DALTON_LIBS}
    )

if(MPI_FOUND)
    set(DALTON_LIBS
        ${DALTON_LIBS}
        mpi_cxx
        )
endif()

add_definitions(-DENABLE_EFS)

unset(OMP_LIB)
