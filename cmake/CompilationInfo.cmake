# create compilation info
configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/compilation_info.py.in
    ${CMAKE_BINARY_DIR}/compilation_info.py
    )
execute_process(
    COMMAND python ${CMAKE_BINARY_DIR}/compilation_info.py
    TIMEOUT 1
    OUTPUT_FILE ${CMAKE_SOURCE_DIR}/gp/compilation_info.F90
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
