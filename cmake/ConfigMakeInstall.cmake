set(INSTALL_DIRECTORY "dalton")

# create install directory
install(
    DIRECTORY
    DESTINATION
    ${INSTALL_DIRECTORY}
    )

if(ENABLE_CHEMSHELL)
    set(LIST_OF_EXECUTABLES          lsdalton.x lslib_tester.x)
else()
    set(LIST_OF_EXECUTABLES dalton.x lsdalton.x lslib_tester.x)
endif()

foreach(_executable ${LIST_OF_EXECUTABLES})
    install(
        TARGETS ${_executable}
        DESTINATION ${INSTALL_DIRECTORY}
        PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ             GROUP_EXECUTE
        WORLD_READ             WORLD_EXECUTE
        )
endforeach()

foreach(_script dalton lsdalton)
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
