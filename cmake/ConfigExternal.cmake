include(ExternalProject)

macro(add_external _project)

    add_custom_target(
        check_external_timestamp_${_project}
        COMMAND python ${PROJECT_SOURCE_DIR}/cmake/check_external_timestamp.py
                       ${PROJECT_BINARY_DIR}/external/${_project}-stamp/${_project}-configure
                       ${PROJECT_BINARY_DIR}/external/${_project}-stamp
                       ${PROJECT_SOURCE_DIR}/external/${_project}
    )

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

    link_directories(${PROJECT_BINARY_DIR}/external/lib)
    link_directories(${PROJECT_BINARY_DIR}/external/${_project}-build/external/lib)

    add_dependencies(${_project} check_external_timestamp_${_project})

endmacro()
