if(NOT DEFINED DEFAULT_Fortran_FLAGS_SET)

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU) # this is gfortran
    add_definitions(-DVAR_GFORTRAN)
    set(CMAKE_Fortran_FLAGS         "-DVAR_GFORTRAN -DGFORTRAN=445 -ffloat-store")
    if(${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "i386")
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -m64"
            )
    endif()
    if(${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "x86_64")
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -m64"
            )
    endif()
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ffast-math -funroll-loops -ftree-vectorize")
    set(CMAKE_Fortran_FLAGS_PROFILE "-O3 -ffast-math -funroll-loops -ftree-vectorize -g -pg")
    if(ENABLE_CODE_COVERAGE)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -fprofile-arcs -ftest-coverage"
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
            "${CMAKE_Fortran_FLAGS} -ftest-coverage"
            )
    endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES G95)
    add_definitions(-DVAR_G95)
    set(CMAKE_Fortran_FLAGS         "-fno-second-underscore -ftrace=full -DVAR_G95")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fsloppy-char")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fsloppy-char -g -pg")
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -i8"
            )
    endif()
    if(ENABLE_BOUNDS_CHECK)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -Wall -fbounds-check"
            )
    endif()
    if(ENABLE_CODE_COVERAGE)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS}"
            )
    endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)

    # example code to determine ifort version
    exec_program(ifort
        ARGS -v
        OUTPUT_VARIABLE IFORT_OUTPUT
        RETURN_VALUE IFORT_INFO
        )
    STRING(REGEX MATCH "[0-9]+" IFORT_VERSION "${IFORT_OUTPUT}")
  # message("-- hey i found the following ifort version: ${IFORT_VERSION}")
  # if(${IFORT_VERSION} STREQUAL "11")
  #     ...
  # else()
  #     ...
  # endif()
    add_definitions(-DVAR_IFORT)
    set(CMAKE_Fortran_FLAGS         "-w -assume byterecl -DVAR_IFORT")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -traceback")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ip")
    set(CMAKE_Fortran_FLAGS_PROFILE "-O3 -ip -g -pg")

    if(ENABLE_OMP)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -openmp -parallel"
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
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
    add_definitions(-DVAR_PGI)
    set(CMAKE_Fortran_FLAGS         "-DVAR_PGI -Mpreprocess")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-g -O0 -Mlarge_arrays -mcmodel=medium")
# I would like to add -fast but this makes certain dec tests fails
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -Mlarge_arrays -mcmodel=medium -Mipa=fast")
    set(CMAKE_Fortran_FLAGS_PROFILE "-O3 -g -pg -Mlarge_arrays -mcmodel=medium")
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -i8 -i8storage"
            )
    endif()
    if(ENABLE_OMP) 
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -mp -Mconcur"
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
    set(CMAKE_Fortran_FLAGS         "-qextname -qlanglvl=extended -qinit=f90ptr")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-g")
    set(CMAKE_Fortran_FLAGS_RELEASE "-qstrict -O3")
    set(CMAKE_Fortran_FLAGS_PROFILE "-qstrict -O3 -p -pg")
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -q64"
            )
    endif()
    set_source_files_properties(${FREE_FORTRAN_SOURCES}
        PROPERTIES COMPILE_FLAGS
        "-qfree"
        )
# -qsuffix=f=f90:cpp=f90 -d
    set_source_files_properties(${FIXED_FORTRAN_SOURCES}
        PROPERTIES COMPILE_FLAGS
        "-qfixed"
        )
    set_source_files_properties(${OWN_BLAS_SOURCES}
        PROPERTIES COMPILE_FLAGS
        "-qfixed"
        )
    set_source_files_properties(${OWN_LAPACK_SOURCES}
        PROPERTIES COMPILE_FLAGS
        "-qfixed"
        )
    if(ENABLE_BOUNDS_CHECK)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -C"
            )
    endif()

endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES Cray) 
    add_definitions(-DVAR_CRAY)
    set(CMAKE_Fortran_FLAGS         "-DVAR_CRAY -eZ")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
    set(CMAKE_Fortran_FLAGS_RELEASE " ")
    set(CMAKE_Fortran_FLAGS_PROFILE "-g")
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -s integer64"
            )
    endif()
    if(ENABLE_BOUNDS_CHECK)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -R bps"
            )
    endif()
endif()

save_compiler_flags(Fortran)
endif()
