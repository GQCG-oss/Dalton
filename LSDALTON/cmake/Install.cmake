set(INSTALL_DIRECTORY "lsdalton")

install(
    DIRECTORY
    DESTINATION
    ${INSTALL_DIRECTORY}
    )

foreach(
    EXECUTABLE
    lsdalton.x
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
    FILES
    ${CMAKE_BINARY_DIR}/lsdalton
    DESTINATION
    ${INSTALL_DIRECTORY}
    PERMISSIONS
    OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ             GROUP_EXECUTE
    WORLD_READ             WORLD_EXECUTE
    )

install(
    DIRECTORY
    ${PROJECT_SOURCE_DIR}/basis
    DESTINATION
    ${INSTALL_DIRECTORY}
    PATTERN .git EXCLUDE
    )
