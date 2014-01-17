include(CheckFortranSourceCompiles)

if(ENABLE_SGI_MPT)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -lmpi")
    set(CMAKE_C_FLAGS       "${CMAKE_C_FLAGS}       -lmpi")
    set(MPI_FOUND TRUE)
endif()

if(ENABLE_MPI)
    if(ENABLE_CRAY_WRAPPERS)
        message("-- Use CRAY wrappers; this disables MPI detection")
        set(MPI_FOUND TRUE)
    else()
        find_package(MPI)
        if(MPI_FOUND)
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_COMPILE_FLAGS}")
            include_directories(${MPI_INCLUDE_PATH})
        else()
            message(FATAL_ERROR "-- You asked for MPI, but CMake could not find any MPI installation, check $PATH")
        endif()
    endif()
endif()

if(MPI_FOUND)

    # test whether we are able to compile a simple MPI program
    file(READ "${CMAKE_SOURCE_DIR}/cmake/mpi/test-MPI-compatibility.F90" _source)
    check_fortran_source_compiles(
        ${_source}
        MPI_COMPATIBLE
        )
  # if(NOT MPI_COMPATIBLE)
  #     message(FATAL_ERROR "Your compiler is not MPI compatible")
  # endif()

    add_definitions(-DVAR_MPI)

    # test whether MPI module is compatible with compiler
    file(READ "${CMAKE_SOURCE_DIR}/cmake/mpi/test-MPI-compiler-compatibility.F90" _source)
    check_fortran_source_compiles(
        ${_source}
        MPI_COMPILER_MATCHES
        )
    if(MPI_COMPILER_MATCHES)
        message("-- mpi.mod matches current compiler, setting -DUSE_MPI_MOD_F90")
        add_definitions(-DUSE_MPI_MOD_F90)
    else()
        message("-- WARNING: mpi.mod compiled with different compiler")
    endif()

    # test whether MPI integer type matches
    file(READ "${CMAKE_SOURCE_DIR}/cmake/mpi/test-MPI-itype-compatibility.F90" _source)
    check_fortran_source_compiles(
        ${_source}
        MPI_ITYPE_MATCHES
        )
    if(NOT MPI_ITYPE_MATCHES)
        if(ENABLE_64BIT_INTEGERS)
            message("-- No 64-bit integer MPI interface found, will use 32-bit integer MPI interface")
            add_definitions(-DVAR_MPI_32BIT_INT)
            set(USE_32BIT_MPI_INTERFACE TRUE)
        else()
            message("-- WARNING: Cannot determine whether MPI is built for 32bit integers")
        endif()
    endif()

    if(ENABLE_64BIT_INTEGERS AND FORCE_32BIT_MPI_INTERFACE)
        message("-- 32-bit integer MPI interface activated by the user")
        add_definitions(-DVAR_MPI_32BIT_INT)
        set(USE_32BIT_MPI_INTERFACE TRUE)
    endif()

    if(ENABLE_MPI2_DETECTION)
        # test whether MPI-2 is available
        file(READ "${CMAKE_SOURCE_DIR}/cmake/mpi/test-MPI-2-compatibility.F90" _source)
        check_fortran_source_compiles(
            ${_source}
            MPI_2_COMPATIBLE
            )
        if(MPI_2_COMPATIBLE)
            add_definitions(-DVAR_MPI2)
            message("-- MPI-2 support found")
        else()
            message("-- no MPI-2 support found, will try with MPI-1")
        endif()
    else()
        message("-- MPI-2 support disabled, will try with MPI-1")
    endif()

endif()
