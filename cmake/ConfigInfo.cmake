message("-- System: ${CMAKE_SYSTEM_NAME}")

if(cmake_build_type_tolower STREQUAL "debug")
    message("-- Fortran Compiler Flags: ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG}")
    message("-- C Compiler Flags:       ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG}")
else()
    message("-- Fortran Compiler Flags: ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE}")
    message("-- C Compiler Flags:       ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}")
endif()

# get size of static allocations
add_custom_target(
    get_static_dalton
    COMMAND ${CMAKE_SOURCE_DIR}/../cmake/binary-info/get_static_size.py dalton.x
    )
add_custom_target(
    get_static_lsdalton
    COMMAND ${CMAKE_SOURCE_DIR}/../cmake/binary-info/get_static_size.py lsdalton.x
    )
