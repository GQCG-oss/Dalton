add_library(qfitlib_interface DALTON/qfitlib/qfitlib_interface.F90)
add_dependencies(dalton qfitlib_interface)
set(DALTON_LIBS qfitlib_interface ${DALTON_LIBS})
if(ENABLE_QFITLIB)
    add_definitions(-DBUILD_QFITLIB)
    if(ENABLE_GEN1INT)
        set(QFIT_INTEGRAL_LIBRARY "GEN1INT")
    else()
        message(FATAL_ERROR "-- QFITlib requires Gen1Int, use -DENABLE_GEN1INT=ON to enable Gen1Int or -DENABLE_QFITLIB=OFF to disable QFITlib")
    endif()
    set(QFIT_INCLUDE_DIR)
    set(QFIT_MPIF OFF)
    if(MPI_FOUND)
        if(NOT MPI_COMPILER_MATCHES)
            set(QFIT_MPIF ON)
        endif()
        set(QFIT_INCLUDE_DIR ${MPI_INCLUDE_PATH})
    endif()
    set(ExternalProjectCMakeArgs
        -DPARENT_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DPARENT_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        -DPARENT_INCLUDE_DIR=${QFIT_INCLUDE_DIR}
        -DPARENT_MODULE_DIR=${PROJECT_BINARY_DIR}/modules
        -DINTEGRAL_LIBRARY=${QFIT_INTEGRAL_LIBRARY}
        -DENABLE_64BIT_INTEGERS=${ENABLE_64BIT_INTEGERS}
        -DENABLE_STATIC_LINKING=${ENABLE_STATIC_LINKING}
        -DENABLE_MPI=${ENABLE_MPI}
        -DENABLE_MPIF=${QFIT_MPIF}
        -DHOST_PROGRAM=DALTON
        )
    add_external(qfitlib)
    add_dependencies(qfitlib_interface qfitlib)
    if(ENABLE_GEN1INT)
        add_dependencies(qfitlib gen1int_interface)
    endif()
    set(EXTERNAL_LIBS ${PROJECT_BINARY_DIR}/external/lib/libqfitlib.a ${EXTERNAL_LIBS})
endif()
