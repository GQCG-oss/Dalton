set(DALTON_LIBS)

add_library(
    dalton
    ${DALTON_C_SOURCES}
    ${DALTON_FREE_FORTRAN_SOURCES}
    ${DALTON_FIXED_FORTRAN_SOURCES}
    ${CMAKE_BINARY_DIR}/binary_info.F90
    )

add_dependencies(dalton generate_binary_info)

if(ENABLE_GEN1INT)
    add_subdirectory(DALTON/gen1int ${CMAKE_BINARY_DIR}/gen1int)
    add_dependencies(dalton gen1int_interface)
    set(DALTON_LIBS
        gen1int_interface
        ${PROJECT_BINARY_DIR}/external/lib/libgen1int.a
        ${DALTON_LIBS}
        )
endif()

if(ENABLE_PELIB)
    include(LibsPElib)
    add_dependencies(dalton pelib)
endif()

add_executable(
    dalton.x
    ${CMAKE_SOURCE_DIR}/DALTON/abacus/dalton.F
    )

set_property(TARGET dalton.x PROPERTY LINKER_LANGUAGE Fortran)

target_link_libraries(
    dalton.x
    dalton
    ${DALTON_LIBS}
    ${EXTERNAL_LIBS}
    )

# compile utilities

file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/tools) 

add_library(peter_utils_blocks ${CMAKE_SOURCE_DIR}/DALTON/tools/blocks.f90)

add_executable(tools/aces2dalton ${CMAKE_SOURCE_DIR}/DALTON/tools/aces2dalton.f90)
add_executable(tools/xyz2dalton  ${CMAKE_SOURCE_DIR}/DALTON/tools/xyz2dalton.f90)
add_executable(tools/distances   ${CMAKE_SOURCE_DIR}/DALTON/tools/distances.f90)

target_link_libraries(tools/aces2dalton peter_utils_blocks)
target_link_libraries(tools/xyz2dalton  peter_utils_blocks)
target_link_libraries(tools/distances   peter_utils_blocks)

add_executable(tools/FChk2HES ${CMAKE_SOURCE_DIR}/DALTON/tools/FChk2HES.f)
add_executable(tools/labread  ${CMAKE_SOURCE_DIR}/DALTON/tools/labread.f)
# radovan: compilation broken
#add_executable(tools/ODCPRG   ${CMAKE_SOURCE_DIR}/DALTON/tools/ODCPRG.f)
