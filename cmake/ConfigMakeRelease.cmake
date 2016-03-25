if(DEVELOPMENT_CODE)

include(InstallRequiredSystemLibraries)

set(EXPORT_DIR ${CMAKE_BINARY_DIR}/export)

set(CPACK_SOURCE_GENERATOR             "TGZ")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY  "DALTON")
set(CPACK_PACKAGE_VENDOR               "DALTON developers")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER    "DALTON developers")
set(CPACK_RESOURCE_FILE_LICENSE        "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
set(CPACK_PACKAGE_FILE_NAME            "DALTON")
set(CPACK_SOURCE_PACKAGE_FILE_NAME     "DALTON-Source")
set(CPACK_SOURCE_INSTALLED_DIRECTORIES "${EXPORT_DIR};/")

include(CPack)

set(DIRS_TO_RELEASE basis cmake CMakeLists.txt COPYING CTestConfig.cmake DALTON external setup VERSION .gitignore CHANGELOG.md dalton.in)

add_custom_target(
    release
    COMMAND ${CMAKE_SOURCE_DIR}/maintenance/release/fetch_external_sources ${CMAKE_SOURCE_DIR}
    COMMAND mkdir -p ${EXPORT_DIR}
    COMMAND tar cf ${EXPORT_DIR}/export.tar -C ${CMAKE_SOURCE_DIR} ${DIRS_TO_RELEASE}
    COMMAND tar xf ${EXPORT_DIR}/export.tar -C ${EXPORT_DIR}
    COMMAND rm     ${EXPORT_DIR}/export.tar
    COMMAND ${CMAKE_SOURCE_DIR}/maintenance/release/remove_unreleased_code ${EXPORT_DIR}
    COMMAND echo "${GIT_REVISION}" > ${EXPORT_DIR}/cmake/GIT_HASH
    COMMAND make package_source
    COMMENT "Packaging source files"
    VERBATIM
    )

endif(DEVELOPMENT_CODE)
