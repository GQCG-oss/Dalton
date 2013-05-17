message("-- Fortran compiler is  \"${CMAKE_Fortran_COMPILER}\" (\"${CMAKE_Fortran_COMPILER_ID}\") ")
message("--       C compiler is  \"${CMAKE_C_COMPILER}\" (\"${CMAKE_C_COMPILER_ID}\") ")
message("my CMAKE_HOST_SYSTEM_PROCESSOR is ${CMAKE_HOST_SYSTEM_PROCESSOR}")

# Fortran compilers

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
    set_source_files_properties(${LSDALTONMAIN_FORTRAN_SOURCES} ${LINEARS_SOURCES} ${RSPSOLVER_SOURCES} ${SOLVERUTIL_SOURCES} ${GEOOPT_SOURCES} ${LSUTIL_PRECISION_SOURCES} ${LSUTIL_COMMON_SOURCES} ${LSUTIL_MATRIXM_SOURCES} ${LSUTIL_MATRIXO_SOURCES} ${LSUTIL_MATRIXU_SOURCES} ${LSUTIL_TYPE_SOURCES} ${LSUTILLIB_SOURCES} ${INTERESTLIB_SOURCES} ${FMM_SOURCES} ${DFTFUNC_F_SOURCES}  ${LSINT_SOURCES} ${PBC_FORTRAN_SOURCES} ${DDYNAM_SOURCES} ${DEC_SOURCES} ${RSP_PROPERTIES_SOURCES} ${LSLIB_SOURCES}
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
#CRAY
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

# C compilers

if(CMAKE_C_COMPILER_ID MATCHES GNU)
    set(CMAKE_C_FLAGS         "-std=c99 -DRESTRICT=restrict -DFUNDERSCORE=1 -DHAVE_NO_LSEEK64 -DUSE_UNDERSCORES -ffloat-store")
    if(${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "i386")
        set(CMAKE_C_FLAGS
            "${CMAKE_C_FLAGS} -m64"
            )
    endif()
    if(${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "x86_64")
        set(CMAKE_C_FLAGS
            "${CMAKE_C_FLAGS} -m64"
            )
    endif()
    set(CMAKE_C_FLAGS_DEBUG   "-O0 -g3")
    set(CMAKE_C_FLAGS_RELEASE "-O3 -ffast-math -funroll-loops -ftree-vectorize")
    set(CMAKE_C_FLAGS_PROFILE "-O3 -ffast-math -funroll-loops -ftree-vectorize -g -pg")
#  PROBLEM WITH FMM C CODE
#    if(ENABLE_CODE_COVERAGE)
#        set(CMAKE_C_FLAGS
#            "${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage"
#            )
#    endif()
    if(ENABLE_OMP)
        set(CMAKE_C_FLAGS
            "${CMAKE_C_FLAGS} -fopenmp"
            )
    endif()
endif()

if(CMAKE_C_COMPILER_ID MATCHES Intel)
    set(CMAKE_C_FLAGS         "-wd981 -wd279 -wd383 -vec-report0 -wd1572 -wd177")
    set(CMAKE_C_FLAGS_DEBUG   "-g -O0")
    set(CMAKE_C_FLAGS_RELEASE "-g -O2")
    set(CMAKE_C_FLAGS_PROFILE "-g -O2 -g -pg")
    set(CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} -shared-intel")
    if(ENABLE_OMP)
        set(CMAKE_C_FLAGS
            "${CMAKE_C_FLAGS} -openmp"
            )
    endif()
endif()

if(CMAKE_C_COMPILER_ID MATCHES PGI)
    set(CMAKE_C_FLAGS         "-Mpreprocess")
    set(CMAKE_C_FLAGS_DEBUG   "-g -O0 -c9x")
    set(CMAKE_C_FLAGS_RELEASE "-O3 -c9x")
    set(CMAKE_C_FLAGS_PROFILE "-O3 -g -pg -c9x")
    if(ENABLE_OMP)
        set(CMAKE_C_FLAGS
            "${CMAKE_C_FLAGS} -mp"
            )
    endif()
endif()

if(CMAKE_C_COMPILER_ID MATCHES XL)
    set(CMAKE_C_FLAGS         " ")
    set(CMAKE_C_FLAGS_DEBUG   "-DVAR_DEBUG ")
    set(CMAKE_C_FLAGS_RELEASE " ")
    set(CMAKE_C_FLAGS_PROFILE " ")
endif()

if(CMAKE_C_COMPILER_ID MATCHES Cray)
    set(CMAKE_C_FLAGS         "-DVAR_CRAY -eZ")
    set(CMAKE_C_FLAGS_DEBUG   "-g -O0")
    set(CMAKE_C_FLAGS_RELEASE " ")
    set(CMAKE_C_FLAGS_PROFILE "-g")
endif()
