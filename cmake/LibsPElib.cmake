add_library(pelib_interface DALTON/pelib/pelib_interface.F90)
add_dependencies(dalton pelib_interface)
set(DALTON_LIBS pelib_interface ${DALTON_LIBS})
if(ENABLE_PELIB)
    add_definitions(-DBUILD_PELIB)
    set(PE_HOST_PROGRAM "DALTON")
    if(ENABLE_GEN1INT)
        set(PE_INTEGRAL_LIBRARY "GEN1INT")
    else()
        message(FATAL_ERROR "-- PElib requires Gen1Int, use -DENABLE_GEN1INT=ON to enable Gen1Int or -DENABLE_PELIB=OFF to disable PElib")
    endif()
    set(PE_INCLUDE_DIR)
    set(PE_MPIF OFF)
    if(MPI_FOUND)
        if(NOT MPI_COMPILER_MATCHES)
            set(PE_MPIF ON)
        endif()
        set(PE_INCLUDE_DIR ${MPI_INCLUDE_PATH})
    endif()
    set(ExternalProjectCMakeArgs
        -DPARENT_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DPARENT_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        -DPARENT_INCLUDE_DIR=${PE_INCLUDE_DIR}
        -DPARENT_MODULE_DIR=${PROJECT_BINARY_DIR}/modules
        -DINTEGRAL_LIBRARY=${PE_INTEGRAL_LIBRARY}
        -DENABLE_64BIT_INTEGERS=${ENABLE_64BIT_INTEGERS}
        -DENABLE_STATIC_LINKING=${ENABLE_STATIC_LINKING}
        -DENABLE_MPI=${ENABLE_MPI}
        -DENABLE_MPIF=${PE_MPIF}
        -DHOST_PROGRAM=${PE_HOST_PROGRAM}
        )
    add_external(pelib)
    add_dependencies(pelib_interface pelib)
    if(ENABLE_GEN1INT)
        add_dependencies(pelib gen1int_interface)
    endif()
    set(EXTERNAL_LIBS ${PROJECT_BINARY_DIR}/external/lib/libpelib.a ${EXTERNAL_LIBS})
endif()
