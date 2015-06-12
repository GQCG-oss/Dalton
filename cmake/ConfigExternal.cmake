include(ExternalProject)

macro(add_external _project)

    add_custom_target(
        check_external_timestamp_${_project}
        COMMAND python ${PROJECT_SOURCE_DIR}/cmake/check_external_timestamp.py
                       ${PROJECT_BINARY_DIR}/external/${_project}-stamp/${_project}-configure
                       ${PROJECT_BINARY_DIR}/external/${_project}-stamp
                       ${PROJECT_SOURCE_DIR}/external/${_project}
    )
    
    set(_testing_command "${ARGV1}")
    
    if("${_testing_command}" STREQUAL "")
       ExternalProject_Add(${_project}                                 
           DOWNLOAD_COMMAND echo
           DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}
           SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/${_project}
           BINARY_DIR ${PROJECT_BINARY_DIR}/external/${_project}-build
           STAMP_DIR ${PROJECT_BINARY_DIR}/external/${_project}-stamp
           TMP_DIR ${PROJECT_BINARY_DIR}/external/${_project}-tmp
           INSTALL_DIR ${PROJECT_BINARY_DIR}/external
           CMAKE_ARGS ${ExternalProjectCMakeArgs}
           )
    else()
       # For unfathomable reasons, CMake expects the TEST_COMMAND to be ;-separated list...
       separate_arguments(_testing_command)
       ExternalProject_Add(${_project}                                 
           DOWNLOAD_COMMAND echo
           DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}
           SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/${_project}
           BINARY_DIR ${PROJECT_BINARY_DIR}/external/${_project}-build
           STAMP_DIR ${PROJECT_BINARY_DIR}/external/${_project}-stamp
           TMP_DIR ${PROJECT_BINARY_DIR}/external/${_project}-tmp
           INSTALL_DIR ${PROJECT_BINARY_DIR}/external
           CMAKE_ARGS ${ExternalProjectCMakeArgs}
	   TEST_BEFORE_INSTALL 1
	   TEST_COMMAND "${_testing_command}"
	   LOG_TEST 1
           )
    endif()	    

    link_directories(${PROJECT_BINARY_DIR}/external/lib)
    link_directories(${PROJECT_BINARY_DIR}/external/${_project}-build/external/lib)

    add_dependencies(${_project} check_external_timestamp_${_project})

endmacro()

macro(add_PCMSOLVER)
    set(PCMSOLVER_TESTS ON)
    set(ExternalProjectCMakeArgs
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
	-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
	-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DENABLE_64BIT_INTEGERS=${ENABLE_64BIT_INTEGERS}
        -DENABLE_BOUNDS_CHECK=${ENABLE_BOUNDS_CHECK}
        -DENABLE_CODE_COVERAGE=${ENABLE_CODE_COVERAGE}
        -DENABLE_STATIC_LINKING=${ENABLE_STATIC_LINKING}
        -DPARENT_MODULE_DIR=${PROJECT_BINARY_DIR}/modules
        -DPARENT_DEFINITIONS=${PARENT_DEFINITIONS}
	-DEIGEN3_ROOT=${EIGEN3_ROOT}
	-DDISABLE_EIGEN_OWN=${DISABLE_EIGEN_OWN}
	-DENABLE_EIGEN_MKL=${ENABLE_EIGEN_MKL}
	-DBOOST_INCLUDEDIR=${BOOST_INCLUDEDIR}
	-DBOOST_LIBRARYDIR=${BOOST_LIBRARYDIR}
	-DENABLE_64BIT_INTEGERS=${ENABLE_64BIT_INTEGERS}
	-DENABLE_PYTHON_EMBEDDING=OFF
	-DENABLE_TESTS=${PCMSOLVER_TESTS}
        )

    if(PCMSOLVER_TESTS)
	    set(testing_command "make test")
	    add_external(pcmsolver ${testing_command})
    else()
	    add_external(pcmsolver)
    endif()
endmacro()
