# - Find a LAPACK library
#
# This module will first look in LAPACK_ROOT before considering the default
# system pahts.
# The linker language can be defined by setting the varable LAPACK_LANG
#
# This module defines:
#
#  LAPACK_INCLUDE_DIRS Where to find blas.h (or equivalent)
#  LAPACK_LIBRARIES Libraries to link against to use LAPACK
#  LAPACK_FOUND Defined if LAPACK is available
#  HAVE_LAPACK To be used in #ifdefs
#
# None of the above will be defined unless LAPACK can be found.
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
    if(NOT DEFINED LAPACK_ROOT})
        set(LAPACK_ROOT $ENV{MATH_ROOT})
    endif()
endif()

if(EXISTS $ENV{LAPACK_ROOT})
    if(NOT DEFINED LAPACK_ROOT})
        set(LAPACK_ROOT $ENV{LAPACK_ROOT})
    endif()
endif()

if(LAPACK_INCLUDE_DIRS AND LAPACK_LIBRARIES)
    set(LAPACK_FIND_QUIETLY TRUE)
endif()

if(NOT LAPACK_FIND_COMPONENTS)
    if(DEFINED LAPACK_TYPE)
        set(LAPACK_FIND_COMPONENTS ${LAPACK_TYPE})
    elseif(ENABLE_64BIT_INTEGERS)
        set(LAPACK_FIND_COMPONENTS MKL)
    else()
        set(LAPACK_FIND_COMPONENTS MKL ESSL Atlas ACML GENERIC)
    endif()
endif()

find_service(LAPACK)
