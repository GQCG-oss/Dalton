include(MathLibsFunctions)

set(MATH_LANG "Fortran")

set(USE_OWN_BLAS   FALSE)
set(USE_OWN_LAPACK FALSE)

set(EXPLICIT_LIBS
    "${EXPLICIT_LIBS}"
    CACHE STRING
    "User set math libraries"
    FORCE
    )

if(ENABLE_CRAY_WRAPPERS)
    message("-- Use CRAY wrappers; this disables math detection and explicit math libraries")
else()
    if(EXPLICIT_LIBS)
        set(LIBS
            ${LIBS}
            ${EXPLICIT_LIBS}
            )
        message("-- User set explicit libraries (skipping BLAS/LAPACK detection): ${EXPLICIT_LIBS}")
    else()
        if(DEFINED MKL_FLAG)
            set(LIBS
                ${LIBS}
                ${MKL_FLAG}
                )
            message("-- User set explicit MKL flag which is passed to the compiler and linker (skipping BLAS/LAPACK detection): ${MKL_FLAG}")
        else()
             config_math_service(LAPACK)
             config_math_service(BLAS)
        endif()
    endif()
endif()
