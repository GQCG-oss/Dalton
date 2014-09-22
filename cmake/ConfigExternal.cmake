
if(DEVELOPMENT_CODE)
    include(FindGit)
    if (GIT_FOUND)
        add_custom_target(
            git_update
            COMMAND ${GIT_EXECUTABLE} submodule init
            COMMAND ${GIT_EXECUTABLE} submodule update
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        )
    else()
        message("-- Git not found. You need Git for the Git submodule mechanism to work.")
    endif()
endif()

include(ExternalProject)

macro(add_external _project)

    if(DEVELOPMENT_CODE AND GIT_FOUND)
        set(UPDATE_COMMAND ${GIT_EXECUTABLE} submodule update)
    else()
        set(UPDATE_COMMAND echo)
    endif()

    # Added by Bin Gao, Jun. 18, 2014
    # If we want to compile the same external project separately for DALTON and LSDALTON,
    # extra arguments should be passed to this macro.
    #
    # Checks if there are extra arguments (note: ARGN can not be used directly with list() command)
    set(extra_macro_args ${ARGN})
    list(LENGTH extra_macro_args num_extra_args)
    if(${num_extra_args} EQUAL 0)
        set(_project_src_dir ${_project})
        # Automatically adds the external project build directory to the include directories
        include_directories(${PROJECT_BINARY_DIR}/external/${_project}-build)
    elseif(${num_extra_args} EQUAL 1)
        # Gets the name of the project source directory
        list(GET extra_macro_args 0 _project_src_dir)
        include_directories(${PROJECT_BINARY_DIR}/external/${_project}-build)
    # The external project build directory needs to be included manually
    elseif(${num_extra_args} EQUAL 2)
        list(GET extra_macro_args 0 _project_src_dir)
        #FIXME: how to save the build directory into a variable named as _project_INCLUDE_DIRS
        # and could be used later??
        #set(_project_INCLUDE_DIRS ${PROJECT_BINARY_DIR}/external/${_project}-build)
    else()
        message(FATAL_ERROR "-- Invalid number of arguments!")
    endif()
    message("-- External project '${_project}' source directory: '${_project_src_dir}'")

    ExternalProject_Add(${_project}
        DOWNLOAD_COMMAND ${UPDATE_COMMAND}
        DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/${_project_src_dir}
        BINARY_DIR ${PROJECT_BINARY_DIR}/external/${_project}-build
        STAMP_DIR ${PROJECT_BINARY_DIR}/external/${_project}-stamp
        TMP_DIR ${PROJECT_BINARY_DIR}/external/${_project}-tmp
        INSTALL_DIR ${PROJECT_BINARY_DIR}/external
        CMAKE_ARGS ${ExternalProjectCMakeArgs}
        )
    link_directories(${PROJECT_BINARY_DIR}/external/lib)
    link_directories(${PROJECT_BINARY_DIR}/external/${_project}-build/external/lib)
    if(DEVELOPMENT_CODE)
        add_dependencies(${_project} git_update)
    endif()

    if(ALWAYS_RESET_EXTERNAL)
        # remove stamps for external builds so that they are rebuilt every time
        add_custom_command(
            TARGET ${_project}
            PRE_BUILD
            COMMAND rm -rf ${PROJECT_BINARY_DIR}/external/${_project}-stamp
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
            )
    endif()
endmacro()

macro(add_PCMSOLVER)
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
        )
    add_external(pcmsolver)
endmacro()
