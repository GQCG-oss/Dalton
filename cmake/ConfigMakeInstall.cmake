set(INSTALL_DIRECTORY "dalton")

# create install directory
install(
    DIRECTORY
    DESTINATION
    ${INSTALL_DIRECTORY}
    )

foreach(
    EXECUTABLE
    dalton.x
    )
    install(
        TARGETS ${EXECUTABLE}
        DESTINATION ${INSTALL_DIRECTORY}
        PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ             GROUP_EXECUTE
        WORLD_READ             WORLD_EXECUTE
        )
endforeach()

install(
    FILES ${CMAKE_BINARY_DIR}/dalton
    DESTINATION ${INSTALL_DIRECTORY}
    PERMISSIONS
    OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ             GROUP_EXECUTE
    WORLD_READ             WORLD_EXECUTE
    )

install(
    DIRECTORY ${PROJECT_SOURCE_DIR}/../basis
    DESTINATION ${INSTALL_DIRECTORY}
    PATTERN .git EXCLUDE
    )

# write git hash to build dir
file(WRITE ${PROJECT_BINARY_DIR}/GIT_HASH "${GIT_REVISION}")

# copy version info to install dir
install(
    FILES ${PROJECT_BINARY_DIR}/GIT_HASH ${PROJECT_SOURCE_DIR}/../VERSION
    DESTINATION ${INSTALL_DIRECTORY}
    PERMISSIONS
    OWNER_READ OWNER_WRITE
    GROUP_READ
    WORLD_READ
    )
