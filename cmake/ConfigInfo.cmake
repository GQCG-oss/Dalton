message("-- System                : ${CMAKE_SYSTEM_NAME}")
message("-- Processor type        : ${CMAKE_HOST_SYSTEM_PROCESSOR}")
message("-- Fortran compiler flags: ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${cmake_build_type_toupper}}")
message("-- C compiler flags      : ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${cmake_build_type_toupper}}")
message("-- Libraries             : ${EXTERNAL_LIBS}")

get_directory_property(LIST_OF_DEFINITIONS DIRECTORY ${CMAKE_SOURCE_DIR} COMPILE_DEFINITIONS)
message("-- Definitions           : ${LIST_OF_DEFINITIONS}")
unset(LIST_OF_DEFINITIONS)

# get size of static allocations
foreach(_binary ${STATIC_MEM_INFO_BINARIES})
    add_custom_target(
        static_mem_${_binary}
            COMMAND ${CMAKE_SOURCE_DIR}/cmake/binary-info/get_static_size.py ${_binary}.x
        )
endforeach()
