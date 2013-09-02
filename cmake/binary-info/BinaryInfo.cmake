include(ConfigGitRevision)

execute_process(
    COMMAND whoami
    TIMEOUT 1
    OUTPUT_VARIABLE USER_NAME
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

execute_process(
    COMMAND hostname
    TIMEOUT 1
    OUTPUT_VARIABLE HOST_NAME
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

execute_process(
    COMMAND python ${CMAKE_SOURCE_DIR}/cmake/binary-info/get_compiler_version.py ${CMAKE_Fortran_COMPILER}
    TIMEOUT 1
    OUTPUT_VARIABLE FORTRAN_COMPILER_VERSION
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

execute_process(
    COMMAND python ${CMAKE_SOURCE_DIR}/cmake/binary-info/get_compiler_version.py ${CMAKE_C_COMPILER}
    TIMEOUT 1
    OUTPUT_VARIABLE C_COMPILER_VERSION
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

execute_process(
    COMMAND python ${CMAKE_SOURCE_DIR}/cmake/binary-info/get_compiler_version.py ${CMAKE_CXX_COMPILER}
    TIMEOUT 1
    OUTPUT_VARIABLE CXX_COMPILER_VERSION
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/binary-info/binary_info.py.in
    ${CMAKE_BINARY_DIR}/binary_info.py
    )

add_custom_target(
    generate_binary_info
    COMMAND python ${CMAKE_BINARY_DIR}/binary_info.py > ${CMAKE_BINARY_DIR}/binary_info.F90
#   COMMAND rm     ${CMAKE_BINARY_DIR}/binary_info.py
    )

set_source_files_properties(${CMAKE_BINARY_DIR}/binary_info.F90 PROPERTIES GENERATED 1)
