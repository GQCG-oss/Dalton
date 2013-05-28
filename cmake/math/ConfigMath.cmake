include(MathLibsFunctions)

set(MATH_LANG "Fortran")

set(EXPLICIT_LIBS
    "${EXPLICIT_LIBS}"
    CACHE STRING
    "User set math libraries"
    FORCE
    )

set(USE_OWN_BLAS   FALSE)
set(USE_OWN_LAPACK FALSE)

# this should move outside of cmake/math
if(EXPLICIT_LIBS)
    set(LIBS
        ${LIBS}
        ${EXPLICIT_LIBS}
        )
    message("-- User set explicit libraries: ${EXPLICIT_LIBS}")
endif()

if(DEFINED MKL_FLAG)
    set(LIBS
        ${LIBS}
        ${MKL_FLAG}
        )
    message("-- User set explicit MKL flag which is passed to the compiler and linker: ${MKL_FLAG}")
else()
    config_math_service(LAPACK)
    config_math_service(BLAS)
endif()
