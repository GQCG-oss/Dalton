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
        ${CMAKE_SOURCE_DIR}/DALTON/main/dalton.F
        )
endif()

add_library(
    dalton
    ${DALTON_C_SOURCES}
    ${DALTON_FREE_FORTRAN_SOURCES}
    ${DALTON_FIXED_FORTRAN_SOURCES}
    ${CMAKE_BINARY_DIR}/binary_info.F90
    )

if(ENABLE_PCMSOLVER)
  add_dependencies(dalton pcmsolver)
  get_target_property(_incdirs dalton INCLUDE_DIRECTORIES)
  set(_incdirs ${_incdirs} ${PROJECT_BINARY_DIR}/external/pcmsolver/src/pcmsolver-build/modules)
  set_target_properties(dalton PROPERTIES INCLUDE_DIRECTORIES "${_incdirs}")
  set(DALTON_LIBS
    ${PCMSOLVER_LIBS}
    ${DALTON_LIBS}
    )
endif()

add_dependencies(dalton generate_binary_info)

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

include(LibsPElib)

include(LibsQFITlib)

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
        ${CMAKE_SOURCE_DIR}/DALTON/main/dalton.F
        )

    set_property(TARGET dalton.x PROPERTY LINKER_LANGUAGE Fortran)

    target_link_libraries(
        dalton.x
        dalton
        ${DALTON_LIBS}
        ${EXTERNAL_LIBS}
        )
endif()

add_subdirectory(DALTON/tools ${CMAKE_BINARY_DIR}/tools)
