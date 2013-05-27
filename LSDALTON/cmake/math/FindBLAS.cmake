

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
#  BLAS_H Name of BLAS header file
#
# None of the above will be defined unless BLAS can be found.
#
#=============================================================================
# Copyright 2011 Jonas Juselius <jonas.juselius@uit.no>
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

include(MathLibFunctions)

if (EXISTS $ENV{MATH_ROOT})
	if (NOT DEFINED BLAS_ROOT})
		set(BLAS_ROOT $ENV{MATH_ROOT})
	endif()
endif()

if (EXISTS $ENV{BLAS_ROOT})
	if (NOT DEFINED BLAS_ROOT})
		set(BLAS_ROOT $ENV{BLAS_ROOT})
	endif()
endif()

if (MATH_LANG STREQUAL "C")
	set(blas_h cblas.h)
	set(blas_libs cblas)
else()
	set(blas_libs blas)
endif()

if (BLAS_INCLUDE_DIRS AND BLAS_LIBRARIES)
  set(BLAS_FIND_QUIETLY TRUE)
endif ()

if (NOT BLAS_FIND_COMPONENTS)
	set(BLAS_FIND_COMPONENTS MKL Atlas ACML default)
endif()

function(find_blas)
	foreach (blas ${BLAS_FIND_COMPONENTS})
		if (${blas} MATCHES "MKL")
			find_mkl()
		elseif (${blas} MATCHES "Atlas")
			find_atlas()
		elseif (${blas} MATCHES "ACML")
			find_acml()
		else()
			find_default()
		endif()
		if (BLAS_FOUND)
			break()
		endif()
	endforeach()
endfunction()

macro(find_default)
	if(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")
	        set(path_suffixes lib64)
	else()
	    set(path_suffixes lib)
	endif()
	set(blas_libs blas)

	find_math_header(blas)
	find_math_libs(blas)
	cache_math_result(default blas)
endmacro()

macro(find_acml)
    set(path_suffixes lib libso)
    set(blas_libs acml)

	find_math_header(blas)
	find_math_libs(blas)
	cache_math_result(ACML blas)
endmacro()

macro(find_atlas)
	set(path_suffixes lib lib/atlas)
	if (MATH_LANG STREQUAL "C")
		set(blas_libs cblas atlas f77blas)
	else()
		set(blas_libs atlas f77blas)
	endif()

	find_math_header(blas)
	find_math_libs(blas)
	cache_math_result(Atlas blas)
endmacro()

macro(find_mkl)
	if (MATH_LANG STREQUAL "C")
		set(blas_h mkl_cblas.h)
	endif()

        set(_mysuffix)
	if(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")
            set(path_suffixes lib/intel64 lib/em64t)
            if(ENABLE_64BIT_INTEGERS)
                set(_mysuffix _ilp64)
            else()
                set(_mysuffix _lp64)
            endif()
	else()
            set(path_suffixes lib/ia32 lib/32)
	endif()

        set(_mythreadtype)
	if(ENABLE_OMP)
            set(_mythreadtype mkl_intel_thread)
            set(_mythreadtype mkl_sequential)
	else()
	endif()

       	if(ENABLE_SCALAPACK)
            set(blas_libs mkl_scalapack${_mysuffix} mkl_core mkl_intel${_mysuffix} ${_mythreadtype} mkl_blacs_intelmpi${_mysuffix} guide pthread m)
        else()
            set(blas_libs                           mkl_core mkl_intel${_mysuffix} ${_mythreadtype}                                guide pthread m)
	endif()
        unset(_mysuffix)
        unset(_mythreadtype)

	find_math_header(blas)
	find_math_libs(blas)
	if(blas_libraries)
		set(blas_libraries -Wl,--start-group ${blas_libraries} -Wl,--end-group )
	endif()
	cache_math_result(MKL blas)
endmacro()

find_blas()

if(BLAS_LIBRARIES)
   set(BLAS_FOUND TRUE)
endif()

unset(blas_h)
unset(blas_libs)
