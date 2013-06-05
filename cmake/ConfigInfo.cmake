message("-- System: ${CMAKE_SYSTEM_NAME}")
message("-- Fortran Compiler Flags: ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${cmake_build_type_toupper}}")
message("-- C Compiler Flags: ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${cmake_build_type_toupper}}")
message("-- Libraries: ${LIBS}")

# get size of static allocations
foreach(_binary ${STATIC_MEM_INFO_BINARIES})
    add_custom_target(
        static_mem_${_binary}
            COMMAND ${CMAKE_SOURCE_DIR}/../cmake/binary-info/get_static_size.py ${_binary}.x
        )
endforeach()
