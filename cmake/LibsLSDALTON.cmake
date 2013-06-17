if(ENABLE_SCALASCA)
    set(SCALASCA_INSTRUMENT ${CMAKE_Fortran_COMPILER})
    configure_script(
        ${CMAKE_SOURCE_DIR}/scalasca.in
        ${CMAKE_SOURCE_DIR}/scalascaf90.sh
        )
    set(SCALASCA_INSTRUMENT ${CMAKE_C_COMPILER})
    configure_script(
        ${CMAKE_SOURCE_DIR}/scalasca.in
        ${CMAKE_SOURCE_DIR}/scalascaCC.sh
        )
    unset(SCALASCA_INSTRUMENT)
    SET(CMAKE_Fortran_COMPILER "../scalascaf90.sh")
    SET(CMAKE_C_COMPILER "../scalascaCC.sh")
    MESSAGE(STATUS "Fortran Compiler " ${CMAKE_Fortran_COMPILER})
    MESSAGE(STATUS "C Compiler       " ${CMAKE_C_COMPILER})
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

add_library(
    lsutillib_common
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

# automatially generate the manual_reorderdings.F90
SET_SOURCE_FILES_PROPERTIES(${CMAKE_SOURCE_DIR}/LSDALTON/lsutil/manual_reorderings.F90 PROPERTIES GENERATED 1)
execute_process(COMMAND python ${CMAKE_SOURCE_DIR}/LSDALTON/lsutil/autogen/generate_man_reord.py nocollapse ${LIST_OF_DEFINITIONS})

add_library(
    lsutiltypelib_common
    ${LSUTIL_TYPE_SOURCES}
    )

target_link_libraries(lsutiltypelib_common pdpacklib)

add_library(
    lsutillib
    ${LSUTILLIB_SOURCES}
    ${GENERATED_FILES}
    )

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
    set(LIBS
        ${PROJECT_BINARY_DIR}/external/lib/libxcfun_f90_bindings.a
        ${PROJECT_BINARY_DIR}/external/lib/libxcfun.a
        ${LIBS}
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

add_library(
    rsp_propertieslib
    ${RSP_PROPERTIES_SOURCES}
    )

target_link_libraries(rsp_propertieslib lsintlib)

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
target_link_libraries(lsdaltonmain rsp_propertieslib)
target_link_libraries(lsdaltonmain rspsolverlib)
target_link_libraries(lsdaltonmain xcfun_interface)

add_executable(
    lsdalton.x
    ${CMAKE_SOURCE_DIR}/LSDALTON/lsdaltonsrc/lsdalton_wrapper.f90
    ${LINK_FLAGS}
    )

add_executable(
    lslib_tester.x
    ${LSLIB_SOURCES}
    ${LINK_FLAGS}
    )

if(MPI_FOUND)
    # Simen's magic fix for Mac/GNU/OpenMPI
    if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
        if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
            SET_TARGET_PROPERTIES(lsdalton.x     PROPERTIES LINK_FLAGS "-Wl,-commons,use_dylibs")
            SET_TARGET_PROPERTIES(lslib_tester.x PROPERTIES LINK_FLAGS "-Wl,-commons,use_dylibs")
        endif()
    endif()
endif()

MERGE_STATIC_LIBS(
    rsp_prop
    rsp_propertieslib
    )

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

MERGE_STATIC_LIBS(
    lsdalton
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
    rsp_prop
    geooptlib
    xcfun_interface
    lsdaltonmain 
    )

target_link_libraries(
    lsdalton
    lsdaltonmain
    ${LIBS}
    )

target_link_libraries(
    lsdalton.x
    lsdalton
    ) 

target_link_libraries(
    lslib_tester.x
    lsdalton
    ) 
