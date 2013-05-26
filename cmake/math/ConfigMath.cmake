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

if(EXPLICIT_LIBS)
    set(LIBS
        ${LIBS}
        ${EXPLICIT_LIBS}
        )
    message("-- User set math libraries: ${EXPLICIT_LIBS}")
else()
    config_math_service(LAPACK)
    config_math_service(BLAS)
endif()
