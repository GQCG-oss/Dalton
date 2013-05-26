# - Find a BLAS library
#
# This module will first look in BLAS_ROOT before considering the default
# system pahts.
# The linker language can be defined by setting the varable BLAS_LANG
#
# This module defines:
#
#  BLAS_INCLUDE_DIRS Where to find blas.h (or equivalent)
#  BLAS_LIBRARIES Libraries to link against to use BLAS
#  BLAS_FOUND Defined if BLAS is available
#  HAVE_BLAS To be used in #ifdefs
#
# None of the above will be defined unless BLAS can be found.
#
#=============================================================================
# Copyright 2011-2013 Jonas Juselius <jonas.juselius@uit.no>
#                     Radovan Bast   <radovan.bast@gmail.com>
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distributed this file outside of CMake, substitute the full
#  License text for the above reference.)

if(EXISTS $ENV{MATH_ROOT})
    if(NOT DEFINED BLAS_ROOT})
        set(BLAS_ROOT $ENV{MATH_ROOT})
    endif()
endif()

if(EXISTS $ENV{BLAS_ROOT})
    if(NOT DEFINED BLAS_ROOT})
        set(BLAS_ROOT $ENV{BLAS_ROOT})
    endif()
endif()

if(BLAS_INCLUDE_DIRS AND BLAS_LIBRARIES)
    set(BLAS_FIND_QUIETLY TRUE)
endif()

if(NOT BLAS_FIND_COMPONENTS)
    if(DEFINED BLAS_TYPE)
        set(BLAS_FIND_COMPONENTS ${BLAS_TYPE})
    elseif(ENABLE_64BIT_INTEGERS)
        set(BLAS_FIND_COMPONENTS MKL)
    else()
        set(BLAS_FIND_COMPONENTS MKL ESSL Atlas ACML GENERIC)
    endif()
endif()

find_service(BLAS)
