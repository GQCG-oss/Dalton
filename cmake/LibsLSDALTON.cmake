set(LSDALTON_EXTERNAL_LIBS)

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
    set(SCALASCA_INSTRUMENT ${CMAKE_CXX_COMPILER})
    configure_script(
        ${CMAKE_SOURCE_DIR}/LSDALTON/scalasca.in
        ${PROJECT_BINARY_DIR}/scalascaCXX.sh
        )
    unset(SCALASCA_INSTRUMENT)
    SET(CMAKE_Fortran_COMPILER "${PROJECT_BINARY_DIR}/scalascaf90.sh")
    SET(CMAKE_C_COMPILER "${PROJECT_BINARY_DIR}/scalascaCC.sh")
    SET(CMAKE_CXX_COMPILER "${PROJECT_BINARY_DIR}/scalascaCXX.sh")
endif()
if(ENABLE_VAMPIRTRACE)
#    set(VAMPIRTRACE_INSTRUMENT ${CMAKE_Fortran_COMPILER})
#    configure_script(
#        ${CMAKE_SOURCE_DIR}/LSDALTON/vampirtrace.in
#        ${PROJECT_BINARY_DIR}/vampirtracef90.sh
#        )
#    set(VAMPIRTRACE_INSTRUMENT ${CMAKE_C_COMPILER})
#    configure_script(
#        ${CMAKE_SOURCE_DIR}/LSDALTON/vampirtrace.in
#        ${PROJECT_BINARY_DIR}/vampirtraceCC.sh
#        )
#    set(VAMPIRTRACE_INSTRUMENT ${CMAKE_CXX_COMPILER})
#    configure_script(
#        ${CMAKE_SOURCE_DIR}/LSDALTON/vampirtrace.in
#        ${PROJECT_BINARY_DIR}/vampirtraceCXX.sh
#        )
#    unset(VAMPIRTRACE_INSTRUMENT)
    SET(CMAKE_Fortran_COMPILER "vtfort")
    SET(CMAKE_C_COMPILER "vtcc")
    SET(CMAKE_CXX_COMPILER "vtc++")
endif()


add_library(
    lsutillib_precision
    ${LSUTIL_PRECISION_SOURCES}
    )

add_library(
    cuda_gpu_interfaces
    ${CUDA_GPU_INTERFACE_SOURCES}
    )

target_link_libraries(cuda_gpu_interfaces lsutillib_precision)

add_library(
    matrixmlib
    ${LSUTIL_MATRIXM_SOURCES}
    )

target_link_libraries(matrixmlib cuda_gpu_interfaces)

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
    ${CMAKE_BINARY_DIR}/manual_reordering/reord2d_1_utils_f2t.F90
    ${CMAKE_BINARY_DIR}/manual_reordering/reord2d_2_utils_f2t.F90
    ${CMAKE_BINARY_DIR}/manual_reordering/reord3d_1_utils_f2t.F90
    ${CMAKE_BINARY_DIR}/manual_reordering/reord3d_2_utils_f2t.F90
    ${CMAKE_BINARY_DIR}/manual_reordering/reord3d_3_utils_f2t.F90
    ${CMAKE_BINARY_DIR}/manual_reordering/reord4d_1_utils_f2t.F90
    ${CMAKE_BINARY_DIR}/manual_reordering/reord4d_2_utils_f2t.F90
    ${CMAKE_BINARY_DIR}/manual_reordering/reord4d_3_utils_f2t.F90
    ${CMAKE_BINARY_DIR}/manual_reordering/reord4d_4_utils_f2t.F90
    ${CMAKE_BINARY_DIR}/manual_reordering/reord2d_1_utils_t2f.F90
    ${CMAKE_BINARY_DIR}/manual_reordering/reord2d_2_utils_t2f.F90
    ${CMAKE_BINARY_DIR}/manual_reordering/reord3d_1_utils_t2f.F90
    ${CMAKE_BINARY_DIR}/manual_reordering/reord3d_2_utils_t2f.F90
    ${CMAKE_BINARY_DIR}/manual_reordering/reord3d_3_utils_t2f.F90
    ${CMAKE_BINARY_DIR}/manual_reordering/reord4d_1_utils_t2f.F90
    ${CMAKE_BINARY_DIR}/manual_reordering/reord4d_2_utils_t2f.F90
    ${CMAKE_BINARY_DIR}/manual_reordering/reord4d_3_utils_t2f.F90
    ${CMAKE_BINARY_DIR}/manual_reordering/reord4d_4_utils_t2f.F90
    )
if(ENABLE_GPU)
    set(MANUAL_REORDERING_SOURCES ${MANUAL_REORDERING_SOURCES}
        ${CMAKE_BINARY_DIR}/manual_reordering/reord2d_acc_reord.F90
        ${CMAKE_BINARY_DIR}/manual_reordering/reord3d_acc_reord.F90
        ${CMAKE_BINARY_DIR}/manual_reordering/reord4d_acc_reord.F90
       )
endif()

get_directory_property(LIST_OF_DEFINITIONS DIRECTORY ${CMAKE_SOURCE_DIR} COMPILE_DEFINITIONS)
if(ENABLE_GPU)
add_custom_command(
    OUTPUT
    ${MANUAL_REORDERING_SOURCES}
    COMMAND
    python ${CMAKE_SOURCE_DIR}/LSDALTON/lsutil/autogen/generate_man_reord.py CMAKE_BUILD=${CMAKE_BINARY_DIR}/manual_reordering acc ${LIST_OF_DEFINITIONS}
    DEPENDS
    ${CMAKE_SOURCE_DIR}/LSDALTON/lsutil/autogen/generate_man_reord.py
    )
elseif(ENABLE_COLLAPSE)
add_custom_command(
    OUTPUT
    ${MANUAL_REORDERING_SOURCES}
    COMMAND
    python ${CMAKE_SOURCE_DIR}/LSDALTON/lsutil/autogen/generate_man_reord.py CMAKE_BUILD=${CMAKE_BINARY_DIR}/manual_reordering ${LIST_OF_DEFINITIONS}
    DEPENDS
    ${CMAKE_SOURCE_DIR}/LSDALTON/lsutil/autogen/generate_man_reord.py
    )
else()
add_custom_command(
    OUTPUT
    ${MANUAL_REORDERING_SOURCES}
    COMMAND
    python ${CMAKE_SOURCE_DIR}/LSDALTON/lsutil/autogen/generate_man_reord.py CMAKE_BUILD=${CMAKE_BINARY_DIR}/manual_reordering nocollapse ${LIST_OF_DEFINITIONS}
    DEPENDS
    ${CMAKE_SOURCE_DIR}/LSDALTON/lsutil/autogen/generate_man_reord.py
    )
endif()
unset(LIST_OF_DEFINITIONS)

add_library(
    lsutillib_common
    ${MANUAL_REORDERING_SOURCES}
    ${LSUTIL_COMMON_C_SOURCES}
    ${LSUTIL_COMMON_SOURCES}
    )

target_link_libraries(lsutillib_common matrixmlib)

add_library(
    lsutil_tensor_lib
    ${LSUTIL_TENSOR_SOURCES}
    )

target_link_libraries(lsutil_tensor_lib lsutillib_common)

add_library(
    matrixolib
    ${LSUTIL_MATRIXO_SOURCES}
    ${LSUTIL_MATRIXO_C_SOURCES}
    )

target_link_libraries(matrixolib lsutil_tensor_lib)

add_library(
    matrixulib
    ${LSUTIL_MATRIXU_SOURCES}
    )

target_link_libraries(matrixulib matrixolib)

if(ENABLE_PCMSOLVER)
    set(EXTERNAL_LIBS ${PCMSOLVER_LIBS} ${EXTERNAL_LIBS})
    add_library(
        lspcm
        ${LSDALTON_PCM_SOURCES})
endif()

set(ExternalProjectCMakeArgs
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    -DENABLE_64BIT_INTEGERS=${ENABLE_64BIT_INTEGERS}
    -DPARENT_MODULE_DIR=${PROJECT_BINARY_DIR}/modules
    )
if(ENABLE_RSP)
add_external(ls-matrix-defop)
set(LSDALTON_EXTERNAL_LIBS
    ${PROJECT_BINARY_DIR}/external/lib/libmatrix-defop.a
    ${LSDALTON_EXTERNAL_LIBS}
    )

add_dependencies(ls-matrix-defop matrixmlib)
add_dependencies(ls-matrix-defop matrixolib)
add_dependencies(ls-matrix-defop matrixulib)
endif()

add_library(
    pdpacklib
    ${LSDALTON_FIXED_FORTRAN_SOURCES}
    )
# Bin Gao: matrixulib needs subroutines in pdpacklib
target_link_libraries(matrixulib pdpacklib)


add_library(
    lsutiltypelib_common
    ${LSUTIL_TYPE_SOURCES}
    )
add_dependencies(lsutiltypelib_common lsutillib_common)
add_dependencies(lsutiltypelib_common matrixulib)

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
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DPARENT_INCLUDE_DIR=${PROJECT_SOURCE_DIR}/include
        -DPARENT_MODULE_DIR=${PROJECT_BINARY_DIR}/modules
        -DENABLE_64BIT_INTEGERS=${ENABLE_64BIT_INTEGERS}
        -DPARENT_DEFINITIONS="-DVAR_LSDALTON"
        )
    add_external(xcfun)
    include_directories(${PROJECT_BINARY_DIR}/external/xcfun-build)
    add_dependencies(xcfun_interface xcfun)
    add_definitions(-DVAR_XCFUN)
    set(LSDALTON_EXTERNAL_LIBS
        ${PROJECT_BINARY_DIR}/external/lib/libxcfun_f90_bindings.a
        ${PROJECT_BINARY_DIR}/external/lib/libxcfun.a
        ${LSDALTON_EXTERNAL_LIBS}
        )
endif()

if(ENABLE_INTEREST)
    add_definitions(-DVAR_INTEREST)
    add_library(
        interestlib
        ${INTERESTLIB_SOURCES}
        )
    target_link_libraries(interestlib xcfun_interface)
endif()

add_library(
    fmmlib
    ${FMM_SOURCES}
    ${FMM_C_SOURCES}
    )


add_dependencies(fmmlib lsutillib_precision)
add_dependencies(fmmlib lsutillib_common)
add_dependencies(fmmlib lsutiltypelib_common)

if(ENABLE_INTEREST)
    target_link_libraries(fmmlib interestlib)
endif()

if(ENABLE_ICHOR)
add_library(
    ichorintlib
    ${ICHORINT_SOURCES}
    )
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
add_dependencies(lsintlib xcfun_interface)
add_dependencies(lsintlib pdpacklib)
add_dependencies(lsintlib lsutillib)
add_dependencies(lsintlib xcfun_interface)
if(ENABLE_ICHOR)
     add_dependencies(lsintlib ichorintlib)
endif()

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
    solverutillib
    ${SOLVERUTIL_SOURCES}
    )

target_link_libraries(solverutillib lsintlib)

add_library(
    rspsolverlib
    ${RSPSOLVER_SOURCES}
    )

target_link_libraries(rspsolverlib solverutillib)

set(ExternalProjectCMakeArgs
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    -DENABLE_64BIT_INTEGERS=${ENABLE_64BIT_INTEGERS}
    -DPARENT_MODULE_DIR=${PROJECT_BINARY_DIR}/modules
    )
if(ENABLE_RSP)
add_external(ls-openrsp)
set(LSDALTON_EXTERNAL_LIBS
    ${PROJECT_BINARY_DIR}/external/lib/libopenrsp.a
    ${LSDALTON_EXTERNAL_LIBS}
    )

add_dependencies(ls-openrsp ls-matrix-defop)
add_dependencies(ls-openrsp solverutillib)
add_dependencies(ls-openrsp rspsolverlib)
endif()

add_library(
    linearslib
    ${LINEARS_SOURCES}
    )

target_link_libraries(linearslib rspsolverlib)
if(ENABLE_RSP)
add_dependencies(linearslib ls-openrsp)
add_dependencies(linearslib ls-matrix-defop)
endif()

add_library(
    declib
    ${DEC_SOURCES}
    )

target_link_libraries(declib lsutiltypelib_common)
target_link_libraries(declib lsutillib_common)
target_link_libraries(declib lsintlib)
target_link_libraries(declib linearslib)

add_library(
    rsp_propertieslib
    ${RSP_PROPERTIES_SOURCES}
    )
add_dependencies(rsp_propertieslib linearslib)
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

if(ENABLE_PCMSOLVER)
    target_link_libraries(lsdaltonmain lspcm)
    add_dependencies(lsdaltonmain  pcmsolver lspcm)
    add_dependencies(linearslib    pcmsolver lspcm)
    add_dependencies(solverutillib pcmsolver lspcm)
    add_dependencies(lspcm lsutillib lsintlib)
endif()

if(ENABLE_CUDA)
    find_package(CUDA)
endif()
if(CUDA_FOUND)
    # this below is a bit convoluted but here we make
    # sure that the CUDA sources are compiled with GNU always
    # this makes life easier if LSDalton is compiled with Intel
    add_definitions(-DENABLE_CUDA)
    set(ExternalProjectCMakeArgs
        -DCMAKE_C_COMPILER=gcc
        -DCMAKE_CXX_COMPILER=g++
        )
    ExternalProject_Add(cuda_interface
        SOURCE_DIR  ${PROJECT_SOURCE_DIR}/LSDALTON/cuda
        BINARY_DIR  ${PROJECT_BINARY_DIR}/cuda/build
        STAMP_DIR   ${PROJECT_BINARY_DIR}/cuda/stamp
        TMP_DIR     ${PROJECT_BINARY_DIR}/cuda/tmp
        DOWNLOAD_COMMAND ""
        INSTALL_COMMAND ""
        )
    include_directories(${PROJECT_SOURCE_DIR}/src/cuda)
    set(LSDALTON_EXTERNAL_LIBS
        ${PROJECT_BINARY_DIR}/cuda/build/libcuda_interface.a
        ${CUDA_LIBRARIES}
        ${LSDALTON_EXTERNAL_LIBS}
        )
    add_dependencies(lsdaltonmain cuda_interface)
endif()

if(NOT ENABLE_CHEMSHELL)
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

    # we always want to compile lslib_tester.x along with lsdalton.x
    add_dependencies(lsdalton.x lslib_tester.x)

    if(MPI_FOUND)
        # Simen's magic fix for Mac/GNU/OpenMPI
        if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
            if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
                SET_TARGET_PROPERTIES(lsdalton.x     PROPERTIES LINK_FLAGS "-Wl,-commons,use_dylibs")
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
  if(ENABLE_ICHOR)
    MERGE_STATIC_LIBS(
        lsint
	ichorintlib
        lsintlib
        )
  else()
    MERGE_STATIC_LIBS(
        lsint
        lsintlib
        )
  endif()
endif()

set(LIBS_TO_MERGE
    lsutillib_precision
    cuda_gpu_interfaces
    matrixmlib
    lsutillib_common
    lsutil_tensor_lib
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
if(ENABLE_PCMSOLVER)
	set(LIBS_TO_MERGE ${LIBS_TO_MERGE} lspcm)
endif()	

MERGE_STATIC_LIBS(
    rsp_prop
    rsp_propertieslib
    )
set(LIBS_TO_MERGE ${LIBS_TO_MERGE} rsp_prop)

MERGE_STATIC_LIBS(
    lsdalton
    ${LIBS_TO_MERGE}
    )

#DO NOT ALWAYS USE stdc++ SINCE THIS IS ONLY!!!! THE GNU STDC++ LIB
if(CMAKE_Fortran_COMPILER_ID MATCHES Cray)
   set(USE_GNU_STDCXX_LIB "")
else()
   set(USE_GNU_STDCXX_LIB "stdc++")
endif()

target_link_libraries(
    lsdalton
    ${EXTERNAL_LIBS}
    ${LSDALTON_EXTERNAL_LIBS}
    ${USE_GNU_STDCXX_LIB}
    )

if(NOT ENABLE_CHEMSHELL)
    target_link_libraries(
        lsdalton.x
        lsdalton
	${PCMSOLVER_LIBS}
        )

    target_link_libraries(
        lslib_tester.x
        lsdalton
	${PCMSOLVER_LIBS}
        )
endif()

# check the LSDALTON source with a python script
add_custom_command(
   OUTPUT ${PROJECT_BINARY_DIR}/check-source # this is just a dummy
   COMMAND ${CMAKE_COMMAND} -P ${PROJECT_SOURCE_DIR}/cmake/CheckLSDALTON.cmake
   WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)

add_custom_target(
   CheckLSDALTON
   ALL DEPENDS ${PROJECT_BINARY_DIR}/check-source
)

add_dependencies(lsdalton CheckLSDALTON)
