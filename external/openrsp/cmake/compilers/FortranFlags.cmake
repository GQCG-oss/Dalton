if(NOT DEFINED DEFAULT_Fortran_FLAGS_SET)

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU) # this is gfortran
    add_definitions(-DVAR_GFORTRAN)
    set(CMAKE_Fortran_FLAGS         "-DVAR_GFORTRAN -DGFORTRAN=445 -ffloat-store -fcray-pointer")
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
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fbacktrace -fcray-pointer -Wuninitialized")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ffast-math -funroll-loops -ftree-vectorize -w")
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

if(CMAKE_Fortran_COMPILER_ID MATCHES G95)
    add_definitions(-DVAR_G95)
    set(CMAKE_Fortran_FLAGS         "-fno-second-underscore -ftrace=full -DVAR_G95")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fsloppy-char")
    set(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_RELEASE} -g -pg")
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
    add_definitions(-DVAR_IFORT)
    set(CMAKE_Fortran_FLAGS         "-w -fpp -assume byterecl -DVAR_IFORT")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -traceback")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -xW -ip")
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
    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        message("--Switch off warnings due to incompatibility XCode 4 and Intel 11 on OsX 10.6")
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -Qoption,ld,-w"
            )
    endif()
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
    add_definitions(-DVAR_PGI)
    set(CMAKE_Fortran_FLAGS         "-DVAR_PGF90 -mcmodel=medium")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-g -O0 -Mframe")
# I would like to add -fast but this makes certain dec tests fails
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -Mipa=fast")
    set(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_RELEASE} -g -pg")
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
    set(CMAKE_Fortran_FLAGS         "-qzerosize -qextname")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-g")
    set(CMAKE_Fortran_FLAGS_RELEASE "-qstrict -O3")
    set(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_RELEASE} -g -pg")
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -q64"
            )
    endif()
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

if(DEFINED EXTRA_Fortran_FLAGS)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${EXTRA_Fortran_FLAGS}")
endif()

save_compiler_flags(Fortran)
endif()
