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
    set(LIBS
        gen1int_interface
        ${PROJECT_BINARY_DIR}/external/lib/libgen1int.a
        ${LIBS}
        )
endif()

if(ENABLE_PELIB)
    set(PARENT_DEFINITIONS "-DPRG_DALTON -DDALTON_MASTER")
    if(ENABLE_GEN1INT)
        set(PARENT_DEFINITIONS "${PARENT_DEFINITIONS} -DBUILD_GEN1INT")
    else()
        message(FATAL_ERROR "-- Gen1Int not enabled. The PE library requires the Gen1Int library.")
    endif()
    if(MPI_FOUND)
        set(PARENT_DEFINITIONS "${PARENT_DEFINITIONS} -DVAR_MPI")
    endif()
    set(ExternalProjectCMakeArgs
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        -DPARENT_INCLUDE_DIR=${PROJECT_SOURCE_DIR}/DALTON/include
        -DPARENT_MODULE_DIR=${PROJECT_BINARY_DIR}/modules
        -DPARENT_DEFINITIONS=${PARENT_DEFINITIONS}
        )
    add_external(pelib)
    add_dependencies(dalton pelib)
    add_dependencies(pelib gen1int_interface)
    add_definitions(-DBUILD_PELIB)
    set(LIBS
        ${PROJECT_BINARY_DIR}/external/lib/libpelib.a
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
add_library(peter_utils_blocks ${CMAKE_SOURCE_DIR}/DALTON/tools/blocks.f90)
add_executable(aces2dalton ${CMAKE_SOURCE_DIR}/DALTON/tools/aces2dalton.f90)
add_executable(xyz2dalton  ${CMAKE_SOURCE_DIR}/DALTON/tools/xyz2dalton.f90)
add_executable(distances   ${CMAKE_SOURCE_DIR}/DALTON/tools/distances.f90)
target_link_libraries(aces2dalton peter_utils_blocks)
target_link_libraries(xyz2dalton peter_utils_blocks)
target_link_libraries(distances peter_utils_blocks)
