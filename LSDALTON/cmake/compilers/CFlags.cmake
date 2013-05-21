if(NOT DEFINED DEFAULT_C_FLAGS_SET)

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

save_compiler_flags(C)
endif()
