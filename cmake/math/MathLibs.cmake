if(ENABLE_64BIT_INTEGERS)
    set(MATH_LIB_SEARCH_ORDER MKL)
else()
    set(MATH_LIB_SEARCH_ORDER MKL ESSL ATLAS ACML GENERIC)
endif()

#-------------------------------------------------------------------------------
# GENERIC

set(GENERIC_BLAS_INCLUDE_PATH_SUFFIXES)
set(GENERIC_LAPACK_INCLUDE_PATH_SUFFIXES)

set(GENERIC_BLAS_HEADERS   cblas.h)
set(GENERIC_LAPACK_HEADERS clapack.h)

set(GENERIC_BLAS_LIBRARY_PATH_SUFFIXES)
set(GENERIC_LAPACK_LIBRARY_PATH_SUFFIXES)

set(GENERIC_BLAS_LIBS   blas)
set(GENERIC_LAPACK_LIBS lapack)

#-------------------------------------------------------------------------------
# ESSL

set(ESSL_BLAS_INCLUDE_PATH_SUFFIXES)
set(ESSL_LAPACK_INCLUDE_PATH_SUFFIXES)

set(ESSL_BLAS_HEADERS)
set(ESSL_LAPACK_HEADERS)

set(ESSL_BLAS_LIBRARY_PATH_SUFFIXES)
set(ESSL_LAPACK_LIBRARY_PATH_SUFFIXES)

set(ESSL_BLAS_LIBS   essl)
set(ESSL_LAPACK_LIBS essl)

#-------------------------------------------------------------------------------
# ACML

set(ACML_BLAS_INCLUDE_PATH_SUFFIXES)
set(ACML_LAPACK_INCLUDE_PATH_SUFFIXES)

set(ACML_BLAS_HEADERS   cblas.h)
set(ACML_LAPACK_HEADERS clapack.h)

set(ACML_BLAS_LIBRARY_PATH_SUFFIXES   libso)
set(ACML_LAPACK_LIBRARY_PATH_SUFFIXES libso)

set(ACML_BLAS_LIBS   acml)
set(ACML_LAPACK_LIBS acml)

#-------------------------------------------------------------------------------
# ATLAS

set(ATLAS_BLAS_INCLUDE_PATH_SUFFIXES   atlas)
set(ATLAS_LAPACK_INCLUDE_PATH_SUFFIXES atlas)

set(ATLAS_BLAS_HEADERS   cblas.h)
set(ATLAS_LAPACK_HEADERS clapack.h)

set(ATLAS_BLAS_LIBRARY_PATH_SUFFIXES   atlas atlas-base atlas-base/atlas atlas-sse3)
set(ATLAS_LAPACK_LIBRARY_PATH_SUFFIXES atlas atlas-base atlas-base/atlas atlas-sse3)

set(ATLAS_BLAS_LIBS   f77blas cblas atlas)
set(ATLAS_LAPACK_LIBS atlas lapack)

#-------------------------------------------------------------------------------
# MKL

set(MKL_BLAS_INCLUDE_PATH_SUFFIXES)
set(MKL_LAPACK_INCLUDE_PATH_SUFFIXES)

set(MKL_BLAS_HEADERS   mkl_cblas.h)
set(MKL_LAPACK_HEADERS mkl_clapack.h)

if(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")
    set(MKL_BLAS_LIBRARY_PATH_SUFFIXES   intel64 em64t)
    set(MKL_LAPACK_LIBRARY_PATH_SUFFIXES intel64 em64t)
else()
    set(MKL_BLAS_LIBRARY_PATH_SUFFIXES   ia32 32)
    set(MKL_LAPACK_LIBRARY_PATH_SUFFIXES ia32 32)
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    set(COMPILER_SPECIFIC_THREAD_LIB mkl_intel_thread)
    set(COMPILER_SPECIFIC_LIB_PREFIX mkl_intel)
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
    set(COMPILER_SPECIFIC_THREAD_LIB mkl_pgi_thread)
    set(COMPILER_SPECIFIC_LIB_PREFIX mkl_intel)
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    set(COMPILER_SPECIFIC_THREAD_LIB mkl_gnu_thread)
    set(COMPILER_SPECIFIC_LIB_PREFIX mkl_gf)
endif()

set(_mysuffix)
if(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")
    if(ENABLE_64BIT_INTEGERS)
        set(_mysuffix _ilp64)
    else()
        set(_mysuffix _lp64)
    endif()
endif()

set(MKL_BLAS_LIBS mkl_core ${COMPILER_SPECIFIC_THREAD_LIB} guide pthread m)
set(MKL_BLAS_LIBS ${MKL_BLAS_LIBS} ${COMPILER_SPECIFIC_LIB_PREFIX}${_mysuffix})

# newer MKL versions don't have libguide
set(MKL_BLAS_LIBS2 mkl_core ${COMPILER_SPECIFIC_THREAD_LIB} pthread m)
set(MKL_BLAS_LIBS2 ${MKL_BLAS_LIBS2} ${COMPILER_SPECIFIC_LIB_PREFIX}${_mysuffix})

# ancient MKL BLAS
set(MKL_BLAS_LIBS3 mkl guide m)

set(MKL_LAPACK_LIBS mkl_lapack95${_mysuffix} ${COMPILER_SPECIFIC_LIB_PREFIX}${_mysuffix})

# older MKL LAPACK
set(MKL_LAPACK_LIBS2 mkl_lapack)

unset(COMPILER_SPECIFIC_THREAD_LIB)
unset(_mysuffix)
