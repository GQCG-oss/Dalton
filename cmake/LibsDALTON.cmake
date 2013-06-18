add_library(
    dalton
    ${DALTON_C_SOURCES}
    ${DALTON_FREE_FORTRAN_SOURCES}
    ${DALTON_FIXED_FORTRAN_SOURCES}
    ${GENERATED_FILES}
    )

if(ENABLE_GEN1INT)
    add_subdirectory(DALTON/gen1int ${CMAKE_BINARY_DIR}/gen1int)
    add_dependencies(dalton gen1int_interface)
    set(LIBS
        gen1int_interface
        ${PROJECT_BINARY_DIR}/external/lib/libgen1int.a
        ${LIBS}
        )
endif()

add_executable(
    dalton.x
    ${CMAKE_SOURCE_DIR}/DALTON/abacus/dalton.F
    )

set_property(TARGET dalton.x PROPERTY LINKER_LANGUAGE Fortran)

target_link_libraries(
    dalton.x
    dalton
    ${LIBS}
    )

# compile Peter's utilities
add_executable(aces2dalton ${CMAKE_SOURCE_DIR}/DALTON/tools/aces2dalton.f90 ${CMAKE_SOURCE_DIR}/DALTON/tools/blocks.f90)
add_executable(xyz2dalton  ${CMAKE_SOURCE_DIR}/DALTON/tools/xyz2dalton.f90  ${CMAKE_SOURCE_DIR}/DALTON/tools/blocks.f90)
add_executable(distances   ${CMAKE_SOURCE_DIR}/DALTON/tools/distances.f90   ${CMAKE_SOURCE_DIR}/DALTON/tools/blocks.f90)
