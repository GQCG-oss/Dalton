include(CheckFortranSourceCompiles)

set(MPI_FOUND FALSE)

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

   if(NOT MPI_COMPATIBLE)
      message("WARNING: Your compiler does not seem to be MPI compatible")
   endif()

   add_definitions(-DVAR_MPI)

   # test whether we can use an mpi module ("use mpi")
   file(READ "${CMAKE_SOURCE_DIR}/cmake/mpi/test-MPI-f90-mod-i4.F90" _source)
   check_fortran_source_compiles(${_source} MPI_F90_I4)
   file(READ "${CMAKE_SOURCE_DIR}/cmake/mpi/test-MPI-f90-mod-i8.F90" _source)
   check_fortran_source_compiles(${_source} MPI_F90_I8)

   if(MPI_F90_I8)

      if(ENABLE_64BIT_INTEGERS)

         message("-- found 64bit integer mpi module, setting -DUSE_MPI_MOD_F90")
         add_definitions(-DUSE_MPI_MOD_F90)

         if(FORCE_32BIT_MPI_INTERFACE)
            message(FATAL_ERROR "32-bit integer MPI interface forced by the user, but check showed MPI is with 64-bit integers")
         endif()

      else()

         message("INFO -- found 64bit integer mpi module for 32bit integer code, will not set -DUSE_MPI_MOD_F90")
         message("   (i.e. use mpif.h, because 'use mpi' then gives problems for some mpi calls)")

      endif()

   elseif(MPI_F90_I4)


      if(ENABLE_64BIT_INTEGERS)

         message("INFO -- found 32bit integer mpi module, setting -DVAR_MPI_32BIT_INT and not -DUSE_MPI_MOD_F90")
         message("   (i.e. use mpif.h, because 'use mpi' then gives problems for some mpi calls)")

         # hjaaj March 2020
         add_definitions(-DVAR_MPI_32BIT_INT)
         set(USE_32BIT_MPI_INTERFACE TRUE)

      else()

         message("-- found 32bit integer mpi module, setting -DUSE_MPI_MOD_F90")
         add_definitions(-DUSE_MPI_MOD_F90)

      endif()

   else() # any found mpi.mod is not compatible with current fortran compiler

      message("-- WARNING: compatible Fortran mpi module not found, will use mpif.h instead")
      message("   (mpif.h should work also if mpi and code are not using same default integer - 32bit vs. 64bit)")

      if(ENABLE_64BIT_INTEGERS AND FORCE_32BIT_MPI_INTERFACE)

         message("-- 32-bit integer MPI interface activated by the user")
         add_definitions(-DVAR_MPI_32BIT_INT)
         set(USE_32BIT_MPI_INTERFACE TRUE)

      endif()

   endif()

   # test current setup and check mpi 3 features
   if(FORCE_DISABLE_MPI3)
      message("-- Disabled MPI3 features")
   else()
      file(READ "${CMAKE_SOURCE_DIR}/cmake/mpi/test-MPI-3-features-simple.F90" _source)
      check_fortran_source_compiles( ${_source} ENABLE_MPI3_FEATURES)

      if(ENABLE_MPI3_FEATURES)
         message("-- found an MPI 3 compatible MPI lib, setting -DVAR_HAVE_MPI3")
         add_definitions(-DVAR_HAVE_MPI3)
         set(MPI_DEFS "${MPI_DEFS} -DVAR_HAVE_MPI3")
      endif()
   endif()
endif()
