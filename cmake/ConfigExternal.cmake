
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

    add_custom_target(
        check_external_timestamp_${_project}
        COMMAND python ${PROJECT_SOURCE_DIR}/cmake/check_external_timestamp.py
                       ${PROJECT_BINARY_DIR}/external/${_project}-stamp/${_project}-configure
                       ${PROJECT_BINARY_DIR}/external/${_project}-stamp
                       ${PROJECT_SOURCE_DIR}/external/${_project}
    )

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
        add_dependencies(check_external_timestamp_${_project} git_update)
    endif()
    add_dependencies(${_project} check_external_timestamp_${_project})

endmacro()
