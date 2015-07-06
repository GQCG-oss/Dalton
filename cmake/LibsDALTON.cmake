set(DALTON_LIBS)
if(ENABLE_VPOTDAMP)
    add_definitions(-DENABLE_VPOTDAMP)
    add_subdirectory(DALTON/1e_cpp ${CMAKE_BINARY_DIR}/vpotdamp)
    set(DALTON_LIBS
        vpotdamp
        ${LIBS}
        )
endif()

if(ENABLE_EFS)
    include(LibsEFS)
    set(DALTON_FIXED_FORTRAN_SOURCES
        DALTON/abacus/efs_interface.F90
        ${DALTON_FIXED_FORTRAN_SOURCES}
        )
    add_subdirectory(DALTON/efs ${CMAKE_BINARY_DIR}/efs_interface)
    set(DALTON_LIBS
        efs_interface
        ${DALTON_LIBS}
        )
endif()

if(ENABLE_CHEMSHELL)
    set(DALTON_FIXED_FORTRAN_SOURCES
        ${DALTON_FIXED_FORTRAN_SOURCES}
        ${CMAKE_SOURCE_DIR}/DALTON/abacus/dalton.F
        )
endif()

add_library(
    dalton
    ${DALTON_C_SOURCES}
    ${DALTON_FREE_FORTRAN_SOURCES}
    ${DALTON_FIXED_FORTRAN_SOURCES}
    ${CMAKE_BINARY_DIR}/binary_info.F90
    )

add_dependencies(dalton generate_binary_info)


# XCint interface https://github.com/rbast/xcint
option(ENABLE_XCINT "Enable XCint interface" OFF)
if(ENABLE_XCINT)
    add_definitions(-DENABLE_XCINT)

    # has to be deactivated for multiple frequencies
    add_definitions(-DENABLE_XCINT_RESPONSE)

    add_library(dalton_xcint_interface DALTON/xcint/dalton_xcint_interface.F90)

    set(ExternalProjectCMakeArgs
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DPARENT_INCLUDE_DIR=${CMAKE_SOURCE_DIR}/DALTON/include
        -DPARENT_MODULE_DIR=${PROJECT_BINARY_DIR}/modules
        -DENABLE_FORTRAN_INTERFACE=ON
        -DENABLE_MPI=${ENABLE_MPI}
        )
    add_external(xcint)

    add_dependencies(dalton_xcint_interface xcint)
    add_dependencies(dalton dalton_xcint_interface)

    set(EXTERNAL_LIBS
        ${EXTERNAL_LIBS}
        dalton_xcint_interface
        ${PROJECT_BINARY_DIR}/external/lib/libxcint.a
        ${PROJECT_BINARY_DIR}/external/xcint-build/external/lib/libxcfun.a
        ${PROJECT_BINARY_DIR}/external/xcint-build/external/lib/libnumgrid.a
        stdc++
        )
endif()

if(ENABLE_EFS)
    add_dependencies(dalton efs)
endif()

if(ENABLE_GEN1INT)
    add_subdirectory(DALTON/gen1int ${CMAKE_BINARY_DIR}/gen1int)
    add_dependencies(dalton gen1int_interface)
    set(DALTON_LIBS
        gen1int_interface
        ${PROJECT_BINARY_DIR}/external/lib/libgen1int.a
        ${DALTON_LIBS}
        )
endif()

if(ENABLE_LSLIB)
    add_definitions( -DBUILD_LSLIB )
    add_dependencies(dalton lsdalton)
    set(DALTON_LIBS
        lsdalton
        ${DALTON_LIBS}
        )
endif()

if(ENABLE_PELIB)
    include(LibsPElib)
    add_dependencies(dalton pelib)
endif()

if(ENABLE_QFITLIB)
    include(LibsQFITlib)
    add_dependencies(dalton qfitlib)
endif()

if(ENABLE_OPENRSP)
    include(LibsOpenRSP)
endif()

if(ENABLE_PCMSOLVER)
    set(PARENT_DEFINITIONS "-DPRG_DALTON -DDALTON_MASTER")
    if(MPI_FOUND)
        set(PARENT_DEFINITIONS "${PARENT_DEFINITIONS} -DVAR_MPI")
    endif()
    add_dependencies(dalton pcmsolver)
    set(DALTON_LIBS
        ${PCMSOLVER_LIBS}
        ${DALTON_LIBS}
        )
endif()

if(ENABLE_QMMM_CUDA)
    add_subdirectory(external/qmmm_cuda)
    add_dependencies(dalton qmmm_cuda)
    set(DALTON_LIBS
        ${PROJECT_BINARY_DIR}/lib/libqmmm_cuda.a
        ${DALTON_LIBS}
        )
endif()

if(NOT ENABLE_CHEMSHELL)
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
endif()

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
