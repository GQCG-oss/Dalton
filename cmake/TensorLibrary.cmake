option(ENABLE_TENSORS "Enable parallel distributed tensors" OFF)

if(ENABLE_TENSORS)
    add_definitions(-DVAR_ENABLE_TENSORS)

    include(ExternalProject)
    set(_tensor_lib_args
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
        -DENABLE_MPI=${ENABLE_MPI}
        -DENABLE_OMP=${ENABLE_OMP}
        -DENABLE_GPU=${ENABLE_GPU}
        -DENABLE_64BIT_INTEGERS=${ENABLE_64BIT_INTEGERS}
        -DENABLE_MPI_32BIT_INT=${ENABLE_MPI_32BIT_INT}
        -DENABLE_REAL_SP=${ENABLE_REAL_SP}
        -DENABLE_TITANBUILD=${ENABLE_TITANBUILD}
        -DENABLE_CRAY_WRAPPERS=${ENABLE_CRAY_WRAPPERS}
        -DMPI_ITYPE_MATCHES=${MPI_ITYPE_MATCHES}
        -DMPI_COMPILER_MATCHES=${MPI_COMPILER_MATCHES}
        )

    ExternalProject_Add(tensor_lib
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/tensor_lib
        BINARY_DIR ${PROJECT_BINARY_DIR}/external/tensor_lib-build
        STAMP_DIR ${PROJECT_BINARY_DIR}/external/tensor_lib-stamp
        TMP_DIR ${PROJECT_BINARY_DIR}/external/tensor_lib-tmp
        INSTALL_DIR ${PROJECT_BINARY_DIR}/external
        CMAKE_ARGS ${_tensor_lib_args}
        )

    unset(_tensor_lib_args)

    include_directories(${PROJECT_BINARY_DIR}/external/tensor_lib-build/modules)

    set(LSDALTON_EXTERNAL_LIBS
        ${PROJECT_BINARY_DIR}/external/lib/libTENSOR_LIB.a
        ${LSDALTON_EXTERNAL_LIBS}
        )

    add_dependencies(lsutillib_common1 tensor_lib)
endif()
