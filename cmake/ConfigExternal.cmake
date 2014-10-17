
# variables used:
#     DEVELOPMENT_CODE - True if within Dalton Git repository
#                      - False if within exported tarball

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
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/${_project}
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
