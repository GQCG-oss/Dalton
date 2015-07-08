if(CMAKE_Fortran_COMPILER_ID MATCHES GNU) # this is gfortran
    add_definitions(-DVAR_GFORTRAN)
    set(CMAKE_Fortran_FLAGS         "-DVAR_GFORTRAN -DGFORTRAN=445 -ffloat-store -fcray-pointer")
    if(${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "i386")
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -m32"
            )
    endif()
    if(${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "x86_64")
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -m64"
            )
    endif()
    if(NOT DEVELOPMENT_CODE)
        # suppress warnings in exported code
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w")
    endif()
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fbacktrace -fcray-pointer -Wuninitialized")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ffast-math -funroll-loops -ftree-vectorize")
    set(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_RELEASE} -g -pg")
    if(ENABLE_STATIC_LINKING)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -static"
            )
    endif()
    if(ENABLE_OMP)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -fopenmp"
            )
    endif()
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -fdefault-integer-8"
            )
    endif()
    if(ENABLE_BOUNDS_CHECK)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -fbounds-check -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wunderflow"
            )
    endif()
    if(ENABLE_CODE_COVERAGE)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -fprofile-arcs -ftest-coverage"
            )
    endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    add_definitions(-DVAR_IFORT)
    set(CMAKE_Fortran_FLAGS         "-fpp -assume byterecl -DVAR_IFORT")
    if(NOT DEVELOPMENT_CODE)
        # suppress warnings in exported code
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w")
    endif()
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -traceback")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ip")
    set(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_RELEASE} -g -pg")

    if(DEFINED MKL_FLAG)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MKL_FLAG}")
    endif()

    if(ENABLE_STATIC_LINKING)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -static-libgcc -static-intel"
            )
    endif()
    if(ENABLE_OMP)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -openmp"
            )
    endif()
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -i8"
            )
    endif()
    if(ENABLE_BOUNDS_CHECK)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -check all -traceback -debug all -fpstkchk -ftrapuv"
            )
    endif()
    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        message("--Switch off warnings due to incompatibility XCode 4 and Intel 11 on OsX 10.6")
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -Qoption,ld,-w"
            )
    endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)

    add_definitions(-DVAR_PGI)

# Patrick: mcmodel=medium is not available on PGI Free for MacOS X
    set(CMAKE_Fortran_FLAGS         "-DVAR_PGF90")
    if(NOT ${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
       set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mcmodel=medium")
    endif()

# Simen: added to include c++ libraries needed for the final linking 
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -pgcpplibs")

    set(CMAKE_Fortran_FLAGS_DEBUG   "-g -O0 -Mframe")
# I would like to add -fast but this makes certain dec tests fails
    #set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -Mipa=fast")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
    set(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_RELEASE} -g -pg")
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -i8 -i8storage"
            )
    endif()
    if(ENABLE_OMP) 
        set(CMAKE_Fortran_FLAGS
            #"${CMAKE_Fortran_FLAGS} -mp -Mconcur"
            "${CMAKE_Fortran_FLAGS} -mp "
            )
    endif()
# WARNING you may need to add -Mcuda=5.5 
# For now use --extra-fc-flags="-Mcuda=5.5"
# ./setup --fc=pgf90 --cc=pgcc --cxx=pgcpp --openacc --extra-fc-flags="-Mcuda=5.5" build
    if(ENABLE_OPENACC) 
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -acc"
            )
    endif()
    if(ENABLE_BOUNDS_CHECK)
        set(CMAKE_Fortran_FLAGS
#add -Mbounds at some point
            "${CMAKE_Fortran_FLAGS} "
            )
    endif()
    if(ENABLE_CODE_COVERAGE)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} "
            )
    endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES XL)
    add_definitions(-DVAR_XLF)
    set(CMAKE_Fortran_FLAGS         "-qzerosize -qextname")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-g")
    set(CMAKE_Fortran_FLAGS_RELEASE "-qstrict -O3")
    set(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_RELEASE} -g -pg")
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -q64"
            )
    endif()
    set_source_files_properties(${DALTON_FREE_FORTRAN_SOURCES}    PROPERTIES COMPILE_FLAGS "-qfree -qlanglvl=extended -qinit=f90ptr")
    set_source_files_properties(${LSDALTON_FREE_FORTRAN_SOURCES}  PROPERTIES COMPILE_FLAGS "-qfree -qlanglvl=extended -qinit=f90ptr")
    set_source_files_properties(${DALTON_FIXED_FORTRAN_SOURCES}   PROPERTIES COMPILE_FLAGS "-qfixed")
    set_source_files_properties(${LSDALTON_FIXED_FORTRAN_SOURCES} PROPERTIES COMPILE_FLAGS "-qfixed")
    set_source_files_properties(${DALTON_OWN_BLAS_SOURCES}        PROPERTIES COMPILE_FLAGS "-qfixed")
    set_source_files_properties(${DALTON_OWN_LAPACK_SOURCES}      PROPERTIES COMPILE_FLAGS "-qfixed")
    set_source_files_properties(${LSDALTON_OWN_BLAS_SOURCES}      PROPERTIES COMPILE_FLAGS "-qfixed")
    set_source_files_properties(${LSDALTON_OWN_LAPACK_SOURCES}    PROPERTIES COMPILE_FLAGS "-qfixed")
    if(ENABLE_BOUNDS_CHECK)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -C"
            )
    endif()

endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES Cray) 
    add_definitions(-DVAR_CRAY)

    set(CMAKE_Fortran_FLAGS         "-DVAR_CRAY -eZ")
    # Patrick: For cray we want to use the system allocator since it is faster and has less memory requirements than the cray allocator
    #if(ENABLE_TITANBUILD)
    #   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -hsystem_alloc")
    #endif()

    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
    set(CMAKE_Fortran_FLAGS_RELEASE " ")
    set(CMAKE_Fortran_FLAGS_PROFILE "-g")
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -s integer64"
            )
    endif()
# WARNING OpenMP (OMP) is activated as default with cray 
    if(ENABLE_OMP) 
      #do nothing OpenMP activated as default with cray 
    else()
      #can be deactivated using -x omp or -h noomp
      set(CMAKE_Fortran_FLAGS
        "${CMAKE_Fortran_FLAGS} -h noomp"
        )
    endif()
# WARNING OpenACC is activated as default on cray 
    if(ENABLE_OPENACC) 
      #do nothing OpenACC activated as default with cray 
    else()
      #can be deactivated using -x acc or -h noacc
      set(CMAKE_Fortran_FLAGS
        "${CMAKE_Fortran_FLAGS} -h noacc"
        )	  
    endif()
    if(ENABLE_BOUNDS_CHECK)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -R bps"
            )
    endif()
endif()

if(DEFINED EXTRA_Fortran_FLAGS)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${EXTRA_Fortran_FLAGS}")
endif()

save_compiler_flags(Fortran)
