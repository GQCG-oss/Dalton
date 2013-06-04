include(CheckFortranSourceCompiles)

if(ENABLE_MPI)
    if(ENABLE_CRAY_WRAPPERS)
        message("-- Use CRAY wrappers; this disables MPI detection")
        add_definitions(-DVAR_MPI)
        add_definitions(-DVAR_LSMPI)
    else()
        find_package(MPI)
        if(MPI_FOUND)
            add_definitions(-DVAR_MPI)
            add_definitions(-DVAR_LSMPI)
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_COMPILE_FLAGS}")
            include_directories(${MPI_INCLUDE_PATH})

            # FIXME: probably we should run this test also on the CRAY
            # test whether MPI integer type matches
            file(READ "${CMAKE_SOURCE_DIR}/../cmake/mpi/test-MPI-itype-compatibility.F90" _source)
            check_fortran_source_compiles(
                ${_source}
                MPI_ITYPE_MATCHES
                )
            if(NOT MPI_ITYPE_MATCHES)
                if(ENABLE_64BIT_INTEGERS)
                    message("-- No 64-bit integer MPI interface found, will use 32-bit integer MPI interface")
                    set(ENABLE_MPI32 TRUE)
                    add_definitions(-DVAR_LSMPI_32)
                else()
                    message("-- WARNING: Cannot determine whether MPI is built for 32bit integers")
                endif()
            endif()
        endif()
    endif()
endif()
