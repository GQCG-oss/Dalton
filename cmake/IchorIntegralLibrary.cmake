option(ENABLE_ICHOR "Enable ichor integrals" OFF)

if(ENABLE_ICHOR)
    add_definitions(-DVAR_ICHOR)

    include(ExternalProject)
#Even if Dalton is built with debug option Ichor should still be compiled 
#with release options withour code coverage and bounds checking.  
    set(_ichor_args
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
#        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
#        -DENABLE_CODE_COVERAGE=${ENABLE_CODE_COVERAGE}
#	 -DENABLE_BOUNDS_CHECK=${DENABLE_BOUNDS_CHECK}
        -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
	-DENABLE_64BIT_INTEGER=${ENABLE_64BIT_INTEGER}
	-DENABLE_GPU=${DENABLE_GPU}
        -DENABLE_OPENACC=${ENABLE_OPENACC}
	-DENABLE_CUDA=${DENABLE_CUDA}
	-DENABLE_CUBLAS=${DENABLE_CUBLAS}
        -DENABLE_MPI=${ENABLE_MPI}
        -DENABLE_OMP=${ENABLE_OMP}
        -DENABLE_TIMINGS=${ENABLE_TIMINGS}	
        )

    ExternalProject_Add(IchorIntegralLibrary
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/IchorIntegralLibrary
        BINARY_DIR ${PROJECT_BINARY_DIR}/external/IchorIntegralLibrary-build
        STAMP_DIR ${PROJECT_BINARY_DIR}/external/IchorIntegralLibrary-stamp
        TMP_DIR ${PROJECT_BINARY_DIR}/external/IchorIntegralLibrary-tmp
        INSTALL_DIR ${PROJECT_BINARY_DIR}/external
        CMAKE_ARGS ${_ichor_args}
        )

    unset(_ichor_args)

    include_directories(${PROJECT_BINARY_DIR}/external/IchorIntegralLibrary-build/modules)

    set(LSDALTON_EXTERNAL_LIBS
        ${PROJECT_BINARY_DIR}/external/lib/libichorlib.a
        ${LSDALTON_EXTERNAL_LIBS}
        )

    add_dependencies(lsintlib IchorIntegralLibrary)
endif()
