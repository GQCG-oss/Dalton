add_library(pelib_interface DALTON/pelib/pelib_interface.F90)
add_dependencies(dalton pelib_interface)
set(DALTON_LIBS pelib_interface ${DALTON_LIBS})
if(ENABLE_64BIT_INTEGERS)
    message(STATUS "PElib is not compatible with 64-bit integer builds")
    message(WARNING "PElib DISABLED")
    set(ENABLE_PELIB OFF)
endif()
if(ENABLE_PELIB)
    include(GNUInstallDirs)
    add_definitions(-DBUILD_PELIB)
    set(PE_HOST_PROGRAM "DALTON")
    if(ENABLE_GEN1INT)
        set(PE_INTEGRAL_LIBRARY "GEN1INT")
    else()
        message(FATAL_ERROR "PElib requires Gen1Int, use -DENABLE_GEN1INT=ON to enable Gen1Int or -DENABLE_PELIB=OFF to disable PElib")
    endif()
    if(ENABLE_64BIT_INTEGERS)
        set(PE_INTEGER_PRECISION INT64)
    else()
        set(PE_INTEGER_PRECISION INT32)
    endif()
    set(PE_INCLUDE_DIRS)
    set(PE_MPIF OFF)
    if(MPI_FOUND)
        if(NOT MPI_COMPILER_MATCHES)
            set(PE_MPIF ON)
        endif()
        set(PE_INCLUDE_DIRS ${PE_INCLUDE_DIRS} ${MPI_INCLUDE_DIRS})
    endif()
    if(ENABLE_PDE)
        find_package(HDF5 REQUIRED COMPONENTS Fortran)
        set(PE_INCLUDE_DIRS ${PE_INCLUDE_DIRS} ${HDF5_Fortran_INCLUDE_DIRS})
        add_definitions(${HDF5_Fortran_DEFINITIONS})
    endif()
    set(ExternalProjectCMakeArgs
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        -DPARENT_INCLUDE_DIRS=${PE_INCLUDE_DIRS}
        -DPARENT_MODULE_DIRS=${PROJECT_BINARY_DIR}/modules
        -DINTEGRAL_LIBRARY=${PE_INTEGRAL_LIBRARY}
        -DINTEGER_PRECISION=${PE_INTEGER_PRECISION}
        -DENABLE_MPI=${ENABLE_MPI}
        -DUSE_MPIF=${PE_MPIF}
        -DENABLE_PDE=${ENABLE_PDE}
        -DHOST_PROGRAM=${PE_HOST_PROGRAM}
        )
    add_external(pelib)
    add_dependencies(pelib_interface pelib)
    if(ENABLE_GEN1INT)
        add_dependencies(pelib gen1int_interface)
    endif()
    set(EXTERNAL_LIBS ${EXTERNAL_LIBS} ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/libPElib.a)
    if(ENABLE_PDE)
        set(EXTERNAL_LIBS ${EXTERNAL_LIBS} ${HDF5_Fortran_LIBRARIES})
    endif()
    include_directories(${PROJECT_BINARY_DIR}/include)
endif()
