if(ENABLE_SCALASCA)
    set(SCALASCA_INSTRUMENT ${CMAKE_Fortran_COMPILER})
    configure_script(
        ${CMAKE_SOURCE_DIR}/LSDALTON/scalasca.in
        ${PROJECT_BINARY_DIR}/scalascaf90.sh
        )
    set(SCALASCA_INSTRUMENT ${CMAKE_C_COMPILER})
    configure_script(
        ${CMAKE_SOURCE_DIR}/LSDALTON/scalasca.in
        ${PROJECT_BINARY_DIR}/scalascaCC.sh
        )
    unset(SCALASCA_INSTRUMENT)
    SET(CMAKE_Fortran_COMPILER "${PROJECT_BINARY_DIR}/scalascaf90.sh")
    SET(CMAKE_C_COMPILER "${PROJECT_BINARY_DIR}/scalascaCC.sh")
endif()

add_library(
    lsutillib_precision
    ${LSUTIL_PRECISION_SOURCES}
    )

add_library(
    matrixmlib
    ${LSUTIL_MATRIXM_SOURCES}
    )

target_link_libraries(matrixmlib lsutillib_precision)

# automatially generate the manual_reorderdings.F90
set(MANUAL_REORDERING_SOURCES
    ${CMAKE_BINARY_DIR}/manual_reordering/reorder_frontend.F90
    ${CMAKE_BINARY_DIR}/manual_reordering/reord2d_2_reord.F90
    ${CMAKE_BINARY_DIR}/manual_reordering/reord3d_1_reord.F90
    ${CMAKE_BINARY_DIR}/manual_reordering/reord3d_2_reord.F90
    ${CMAKE_BINARY_DIR}/manual_reordering/reord3d_3_reord.F90
    ${CMAKE_BINARY_DIR}/manual_reordering/reord4d_1_reord.F90
    ${CMAKE_BINARY_DIR}/manual_reordering/reord4d_2_reord.F90
    ${CMAKE_BINARY_DIR}/manual_reordering/reord4d_3_reord.F90
    ${CMAKE_BINARY_DIR}/manual_reordering/reord4d_4_reord.F90
    ${CMAKE_BINARY_DIR}/manual_reordering/reord4d_1_utils_f2t.F90
    ${CMAKE_BINARY_DIR}/manual_reordering/reord4d_2_utils_f2t.F90
    ${CMAKE_BINARY_DIR}/manual_reordering/reord4d_3_utils_f2t.F90
    ${CMAKE_BINARY_DIR}/manual_reordering/reord4d_4_utils_f2t.F90
    ${CMAKE_BINARY_DIR}/manual_reordering/reord4d_1_utils_t2f.F90
    ${CMAKE_BINARY_DIR}/manual_reordering/reord4d_2_utils_t2f.F90
    ${CMAKE_BINARY_DIR}/manual_reordering/reord4d_3_utils_t2f.F90
    ${CMAKE_BINARY_DIR}/manual_reordering/reord4d_4_utils_t2f.F90
    )

get_directory_property(LIST_OF_DEFINITIONS DIRECTORY ${CMAKE_SOURCE_DIR} COMPILE_DEFINITIONS)
add_custom_command(
    OUTPUT
    ${MANUAL_REORDERING_SOURCES}
    COMMAND
    python ${CMAKE_SOURCE_DIR}/LSDALTON/lsutil/autogen/generate_man_reord.py nocollapse CMAKE_BUILD=${CMAKE_BINARY_DIR}/manual_reordering ${LIST_OF_DEFINITIONS}
    DEPENDS
    ${CMAKE_SOURCE_DIR}/LSDALTON/lsutil/autogen/generate_man_reord.py
    )
unset(LIST_OF_DEFINITIONS)

add_library(
    lsutillib_common
    ${MANUAL_REORDERING_SOURCES}
    ${LSUTIL_COMMON_SOURCES}
    )

target_link_libraries(lsutillib_common matrixmlib)

add_library(
    matrixolib
    ${LSUTIL_MATRIXO_SOURCES}
    ${LSUTIL_MATRIXO_C_SOURCES}
    )

target_link_libraries(matrixolib lsutillib_common)

add_library(
    matrixulib
    ${LSUTIL_MATRIXU_SOURCES}
    )

target_link_libraries(matrixulib matrixolib)

add_library(
    pdpacklib
    ${LSDALTON_FIXED_FORTRAN_SOURCES}
    )

target_link_libraries(pdpacklib matrixulib)



add_library(
    lsutiltypelib_common
    ${LSUTIL_TYPE_SOURCES}
    )


target_link_libraries(lsutiltypelib_common pdpacklib)

add_library(
    lsutillib
    ${LSUTILLIB_SOURCES}
    ${CMAKE_BINARY_DIR}/binary_info.F90
    )

add_dependencies(lsutillib generate_binary_info)

target_link_libraries(lsutillib lsutiltypelib_common)

add_library(
    xcfun_interface
    LSDALTON/xcfun_host/xcfun_host.F90
    )

add_dependencies(xcfun_interface lsutillib_precision)

if(ENABLE_XCFUN)
    set(ExternalProjectCMakeArgs
        -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DPARENT_INCLUDE_DIR=${PROJECT_SOURCE_DIR}/include
        -DPARENT_MODULE_DIR=${PROJECT_BINARY_DIR}/modules
        -DPARENT_DEFINITIONS="-DVAR_LSDALTON"
        )
    add_external(xcfun)
    add_dependencies(xcfun_interface xcfun)
    add_definitions(-DVAR_XCFUN)
    set(EXTERNAL_LIBS
        ${PROJECT_BINARY_DIR}/external/lib/libxcfun_f90_bindings.a
        ${PROJECT_BINARY_DIR}/external/lib/libxcfun.a
        ${EXTERNAL_LIBS}
        )
endif()

if(ENABLE_INTEREST)
    add_definitions(-DVAR_INTEREST)
    add_library(
        interestlib
        ${INTERESTLIB_SOURCES}
        )
    if(ENABLE_XCFUN)
        target_link_libraries(interestlib xcfun_interface)
    endif()
endif()

add_library(
    fmmlib
    ${FMM_SOURCES}
    ${FMM_C_SOURCES}
    )

if(ENABLE_XCFUN)
    target_link_libraries(fmmlib xcfun_interface)
endif()

add_dependencies(fmmlib lsutillib_precision)
add_dependencies(fmmlib lsutillib_common)
add_dependencies(fmmlib lsutiltypelib_common)

if(ENABLE_INTEREST)
    target_link_libraries(fmmlib interestlib)
endif()

add_library(
    dftfunclib
    ${DFTFUNC_SOURCES}
    ${DFTFUNC_F_SOURCES}
    )

target_link_libraries(dftfunclib fmmlib)

add_library(
    lsintlib
    ${LSINT_SOURCES}
    )

target_link_libraries(lsintlib dftfunclib)
add_dependencies(lsintlib pdpacklib)
add_dependencies(lsintlib lsutillib)
add_dependencies(lsintlib xcfun_interface)

add_library(
    pbclib
    ${PBC_FORTRAN_SOURCES}
    )

target_link_libraries(pbclib lsintlib)

add_library(
    ddynamlib
    ${DDYNAM_SOURCES}
    )

target_link_libraries(ddynamlib lsintlib)

add_library(
    declib
    ${DEC_SOURCES}
    ${DEC_C_SOURCES}
    )

target_link_libraries(declib lsintlib)

add_library(
    solverutillib
    ${SOLVERUTIL_SOURCES}
    )

target_link_libraries(solverutillib lsintlib)

add_library(
    rspsolverlib
    ${RSPSOLVER_SOURCES}
    )

target_link_libraries(rspsolverlib solverutillib)

add_library(
    linearslib
    ${LINEARS_SOURCES}
    )

target_link_libraries(linearslib rspsolverlib)

if(DEVELOPMENT_CODE)
    add_library(
        rsp_propertieslib
        ${RSP_PROPERTIES_SOURCES}
        )
    target_link_libraries(rsp_propertieslib lsintlib)
endif()

add_library(
    geooptlib
    ${GEOOPT_SOURCES}
    )

target_link_libraries(geooptlib lsintlib)

add_library(
    lsdaltonmain 
    ${LSDALTONMAIN_FORTRAN_SOURCES}
    )

target_link_libraries(lsdaltonmain pbclib)
target_link_libraries(lsdaltonmain geooptlib)
target_link_libraries(lsdaltonmain linearslib)
target_link_libraries(lsdaltonmain declib)
target_link_libraries(lsdaltonmain ddynamlib)
if(DEVELOPMENT_CODE)
    target_link_libraries(lsdaltonmain rsp_propertieslib)
endif()
target_link_libraries(lsdaltonmain rspsolverlib)
target_link_libraries(lsdaltonmain xcfun_interface)

add_executable(
    lsdalton.x
    ${CMAKE_SOURCE_DIR}/LSDALTON/lsdaltonsrc/lsdalton_wrapper.f90
    ${LINK_FLAGS}
    )

if(NOT ENABLE_CHEMSHELL)
    add_executable(
        lslib_tester.x
        ${LSLIB_SOURCES}
        ${LINK_FLAGS}
        )

    # we always want to compile lslib_tester.x along with lsdalton.x
    add_dependencies(lsdalton.x lslib_tester.x)
endif()

if(MPI_FOUND)
    # Simen's magic fix for Mac/GNU/OpenMPI
    if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
        if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
            SET_TARGET_PROPERTIES(lsdalton.x     PROPERTIES LINK_FLAGS "-Wl,-commons,use_dylibs")
            if(NOT ENABLE_CHEMSHELL)
                SET_TARGET_PROPERTIES(lslib_tester.x PROPERTIES LINK_FLAGS "-Wl,-commons,use_dylibs")
            endif()
        endif()
    endif()
endif()

if(ENABLE_INTEREST)
    MERGE_STATIC_LIBS(
        lsint
        lsintlib
        interestlib
        )
else()
    MERGE_STATIC_LIBS(
        lsint
        lsintlib
        )
endif()

set(LIBS_TO_MERGE
    lsutillib_precision
    matrixmlib
    lsutillib_common
    matrixolib
    matrixulib
    pdpacklib
    lsutiltypelib_common
    lsutillib
    fmmlib
    dftfunclib
    lsint
    pbclib
    ddynamlib
    declib
    solverutillib
    rspsolverlib
    linearslib
    geooptlib
    xcfun_interface
    lsdaltonmain 
    )

if(DEVELOPMENT_CODE)
    MERGE_STATIC_LIBS(
        rsp_prop
        rsp_propertieslib
        )
    set(LIBS_TO_MERGE ${LIBS_TO_MERGE} rsp_prop)
endif()

MERGE_STATIC_LIBS(
    lsdalton
    ${LIBS_TO_MERGE}
    )

target_link_libraries(
    lsdalton
    lsdaltonmain
    ${EXTERNAL_LIBS}
    )

target_link_libraries(
    lsdalton.x
    lsdalton
    ) 

if(NOT ENABLE_CHEMSHELL)
    target_link_libraries(
        lslib_tester.x
        lsdalton
        )
endif()
