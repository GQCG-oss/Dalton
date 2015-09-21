if (CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS         "-g -Wall -fno-rtti -fno-exceptions")
    if(NOT DEVELOPMENT_CODE)
        # suppress warnings in exported code
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
    endif()
    set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g3")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffast-math")
    if(NOT ${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
        # radovan: vpotdamp code needs this
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
    endif()
    if (ENABLE_CODE_COVERAGE)
        set (CMAKE_CXX_FLAGS
            "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
        set (CMAKE_CXX_LINK_FLAGS "-fprofile-arcs -ftest-coverage")
    endif()
    if(ENABLE_OMP)
        set(CMAKE_CXX_FLAGS
            "${CMAKE_CXX_FLAGS} -fopenmp"
            )
    endif()
    if(ENABLE_STATIC_LINKING)
        set(CMAKE_CXX_FLAGS
            "${CMAKE_CXX_FLAGS} -static -fpic"
            )
    endif()
elseif (CMAKE_CXX_COMPILER_ID MATCHES Intel)
    set(CMAKE_CXX_FLAGS         "-g -wd981 -wd279 -wd383 -vec-report0 -wd1572 -wd177 -fno-rtti")
    if(DEVELOPMENT_CODE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
    else()
        # suppress warnings in exported code
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
    endif()
    set(CMAKE_CXX_FLAGS_DEBUG   "-O0")
    set(CMAKE_CXX_FLAGS_RELEASE "-O2")
    set (CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -shared-intel")
    if(ENABLE_OMP)
        set(CMAKE_CXX_FLAGS
            "${CMAKE_CXX_FLAGS} -openmp"
            )
    endif()
endif ()

if(DEFINED EXTRA_CXX_FLAGS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_CXX_FLAGS}")
endif()

save_compiler_flags(CXX)
