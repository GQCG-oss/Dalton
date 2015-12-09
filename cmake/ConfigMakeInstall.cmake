set(INSTALL_DIRECTORY "dalton")

# create install directory
install(
    DIRECTORY
    DESTINATION
    ${INSTALL_DIRECTORY}
    )

if(NOT ENABLE_CHEMSHELL)
    foreach(_executable dalton.x)
        install(
            TARGETS ${_executable}
            DESTINATION ${INSTALL_DIRECTORY}
            PERMISSIONS
            OWNER_READ OWNER_WRITE OWNER_EXECUTE
            GROUP_READ             GROUP_EXECUTE
            WORLD_READ             WORLD_EXECUTE
            )
    endforeach()
endif()

foreach(_script dalton)
    install(
        FILES ${CMAKE_BINARY_DIR}/${_script}
        DESTINATION ${INSTALL_DIRECTORY}
        PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ             GROUP_EXECUTE
        WORLD_READ             WORLD_EXECUTE
        )
endforeach()

foreach(_directory ${CMAKE_SOURCE_DIR}/basis ${CMAKE_BINARY_DIR}/tools)
    install(
        DIRECTORY ${_directory}
        DESTINATION ${INSTALL_DIRECTORY}
        )
endforeach()

# write git hash to build dir
file(WRITE ${CMAKE_BINARY_DIR}/GIT_HASH "${GIT_REVISION}")

# copy version info to install dir
install(
    FILES ${CMAKE_BINARY_DIR}/GIT_HASH ${CMAKE_SOURCE_DIR}/VERSION
    DESTINATION ${INSTALL_DIRECTORY}
    PERMISSIONS
    OWNER_READ OWNER_WRITE
    GROUP_READ
    WORLD_READ
    )
