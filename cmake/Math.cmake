set(MATH_LANG "Fortran")

set(BLAS_FOUND   FALSE)
set(LAPACK_FOUND FALSE)

option(FORCE_OWN_BLAS   "Force use of own BLAS"   OFF)
option(FORCE_OWN_LAPACK "Force use of own LAPACK" OFF)
mark_as_advanced(FORCE_OWN_BLAS)
mark_as_advanced(FORCE_OWN_LAPACK)

if(NOT ENABLE_INTERNAL_MATH)

    # user sets
    set(USERDEFINED_MATH
        "${USERDEFINED_MATH}"
        CACHE STRING
        "User set math libraries"
        FORCE
        )
    if(NOT "${USERDEFINED_MATH}" STREQUAL "")
        set(MATH_LIBS
            "${USERDEFINED_MATH}"
            CACHE STRING
            "User set math libraries"
            FORCE
            )
        message("-- User set math libraries: ${MATH_LIBS}")
        set(BLAS_FOUND   TRUE)
        set(LAPACK_FOUND TRUE)
        set(LIBS
            ${LIBS}
            ${MATH_LIBS}
            )
    endif()

    # try to find the best libraries using environment variables
    if(NOT FORCE_OWN_BLAS)
        if(NOT BLAS_FOUND)
            find_package(BLAS)
            if(BLAS_FOUND)
                set(LIBS
                    ${LIBS}
                    ${BLAS_LIBRARIES}
                    )
            endif()
        endif()
        if(NOT BLAS_FOUND)
            message("-- No external BLAS library found")
        endif()
    endif()
    if(NOT FORCE_OWN_LAPACK)
        if(NOT LAPACK_FOUND)
            find_package(LAPACK)
            if(LAPACK_FOUND)
                set(LIBS
                    ${LAPACK_LIBRARIES}
                    ${LIBS}
                    )
            endif()
        endif()
        if(NOT LAPACK_FOUND)
            message("-- No external LAPACK library found")
        endif()
    endif()
endif()

if(NOT BLAS_FOUND OR FORCE_OWN_BLAS)
    message("-- Using own BLAS implementation (slow)")
    set(FIXED_FORTRAN_SOURCES
        ${FIXED_FORTRAN_SOURCES}
        ${OWN_BLAS_SOURCES}
        )
endif()
if(NOT LAPACK_FOUND OR FORCE_OWN_LAPACK)
    message("-- Using own LAPACK implementation (slow)")
    set(FIXED_FORTRAN_SOURCES
        ${FIXED_FORTRAN_SOURCES}
        ${OWN_LAPACK_SOURCES}
        )
endif()
