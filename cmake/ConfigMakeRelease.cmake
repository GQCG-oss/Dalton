if(DEVELOPMENT_CODE)

include(InstallRequiredSystemLibraries)

set(CPACK_SOURCE_GENERATOR "TGZ")

set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "DALTON")
set(CPACK_PACKAGE_VENDOR "DALTON developers")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "DALTON developers")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/../INSTALL.rst")
set(CPACK_RESOURCE_FILE_LICENSE    "${CMAKE_CURRENT_SOURCE_DIR}/../COPYING")
set(CPACK_PACKAGE_FILE_NAME        "DALTON-${DALTON_VERSION}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "DALTON-${DALTON_VERSION}-Source")
set(CPACK_PACKAGE_INSTALL_DIRECTORY
    "DALTON ${DALTON_VERSION}"
    )

SET(EXPORT_DIR ${CMAKE_BINARY_DIR}/export)
SET(CPACK_SOURCE_INSTALLED_DIRECTORIES
    "${EXPORT_DIR};/"
    )

include(CPack)

# create list of include files
execute_process(
    COMMAND
    ${CMAKE_SOURCE_DIR}/../maintenance/release/get_list_of_include_files.py ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR}
    )
include(IncludeFiles)

add_custom_target(
    release
    COMMAND ${CMAKE_SOURCE_DIR}/../maintenance/release/fetch_external_sources ${CMAKE_SOURCE_DIR}/..
    COMMAND mkdir -p ${EXPORT_DIR}/DALTON
    COMMAND tar cf ${EXPORT_DIR}/DALTON/sources.tar -C ${CMAKE_SOURCE_DIR} ${DALTON_SOURCES} ${INCLUDE_FILES}
    COMMAND tar xf ${EXPORT_DIR}/DALTON/sources.tar -C ${EXPORT_DIR}/DALTON
    COMMAND rm     ${EXPORT_DIR}/DALTON/sources.tar
    COMMAND cp -r  ${CMAKE_SOURCE_DIR}/../basis          ${EXPORT_DIR}
    COMMAND cp -r  ${CMAKE_SOURCE_DIR}/../cmake          ${EXPORT_DIR}
    COMMAND cp     ${CMAKE_SOURCE_DIR}/CMakeLists.txt    ${EXPORT_DIR}/DALTON
    COMMAND cp     COPYING                               ${EXPORT_DIR}
    COMMAND cp     INSTALL.rst                           ${EXPORT_DIR}
    COMMAND cp     ${CMAKE_SOURCE_DIR}/CTestConfig.cmake ${EXPORT_DIR}/DALTON
    COMMAND cp     ${CMAKE_SOURCE_DIR}/dalton.in         ${EXPORT_DIR}/DALTON
    COMMAND cp     ${CMAKE_SOURCE_DIR}/setup             ${EXPORT_DIR}/DALTON
    COMMAND cp     ${CMAKE_SOURCE_DIR}/../VERSION        ${EXPORT_DIR}
    COMMAND cp -r  ${CMAKE_SOURCE_DIR}/abacus/dalton.F   ${EXPORT_DIR}/DALTON/abacus
    COMMAND cp -r  ${CMAKE_SOURCE_DIR}/cmake             ${EXPORT_DIR}/DALTON
    COMMAND cp -r  ${CMAKE_SOURCE_DIR}/Doc               ${EXPORT_DIR}/DALTON
    COMMAND cp -r  ${CMAKE_SOURCE_DIR}/Doc_cc            ${EXPORT_DIR}/DALTON
    COMMAND cp -r  ${CMAKE_SOURCE_DIR}/test              ${EXPORT_DIR}/DALTON
    COMMAND cp -r  ${CMAKE_SOURCE_DIR}/test_cc           ${EXPORT_DIR}/DALTON
    COMMAND cp -r  ${CMAKE_SOURCE_DIR}/tools             ${EXPORT_DIR}/DALTON
    COMMAND cp -r  ${CMAKE_SOURCE_DIR}/gen1int           ${EXPORT_DIR}/DALTON
    COMMAND ${CMAKE_SOURCE_DIR}/../maintenance/release/remove_unreleased_code ${EXPORT_DIR}
    COMMAND cp -r  ${CMAKE_SOURCE_DIR}/../external       ${EXPORT_DIR}
    COMMAND rm -rf ${EXPORT_DIR}/external/gen1int/.git
    COMMAND rm -rf ${EXPORT_DIR}/external/xcfun/.git
    COMMAND echo "${GIT_REVISION}"                     > ${EXPORT_DIR}/DALTON/cmake/GIT_HASH
    COMMAND make package_source
    COMMENT "Packaging source files"
    VERBATIM
    )

endif(DEVELOPMENT_CODE)
