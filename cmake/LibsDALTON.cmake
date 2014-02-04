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
    if(ENABLE_OPENRSP)
        set(PARENT_DEFINITIONS "-DPRG_DALTON -DDALTON_MASTER -DBUILD_OPENRSP")
    else()
        set(PARENT_DEFINITIONS "-DPRG_DALTON -DDALTON_MASTER")
    endif()
    if(ENABLE_GEN1INT)
        set(PARENT_DEFINITIONS "${PARENT_DEFINITIONS} -DBUILD_GEN1INT")
    else()
        message(FATAL_ERROR "-- Gen1Int not enabled. The PE library requires the Gen1Int library.")
    endif()
    if(MPI_FOUND)
        set(PARENT_DEFINITIONS "${PARENT_DEFINITIONS} -DVAR_MPI")
        if(MPI_COMPILER_MATCHES)
            set(PARENT_DEFINITIONS "${PARENT_DEFINITIONS} -DUSE_MPI_MOD_F90")
        endif()
    endif()
    set(ExternalProjectCMakeArgs
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        -DENABLE_64BIT_INTEGERS=${ENABLE_64BIT_INTEGERS}
        -DENABLE_BOUNDS_CHECK=${ENABLE_BOUNDS_CHECK}
        -DENABLE_CODE_COVERAGE=${ENABLE_CODE_COVERAGE}
        -DENABLE_STATIC_LINKING=${ENABLE_STATIC_LINKING}
        -DPARENT_MODULE_DIR=${PROJECT_BINARY_DIR}/modules
        -DPARENT_DEFINITIONS=${PARENT_DEFINITIONS}
        )
    add_external(pelib)
    add_dependencies(dalton pelib)
    add_dependencies(pelib gen1int_interface)
    add_definitions(-DBUILD_PELIB)
    set(DALTON_LIBS
        ${PROJECT_BINARY_DIR}/external/lib/libpelib.a
        ${DALTON_LIBS}
        )
endif()

if(ENABLE_OPENRSP)
    include(LibsOpenRSP)
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
