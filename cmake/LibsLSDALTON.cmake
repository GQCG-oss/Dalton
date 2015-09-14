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

add_library(
    lsutillib_common1
    ${LSUTIL_COMMON_C_SOURCES}
    ${LSUTIL_TYPE_SOURCES}
    )

if(NOT ENABLE_TENSORS)
   target_link_libraries(lsutillib_common1 matrixmlib)
endif()

add_library(
    lsutillib_common2
    ${LSUTIL_COMMON_SOURCES2}
    )

target_link_libraries(lsutillib_common2 lsutillib_common1)

add_library(
    lsutillib_common3
    ${LSUTIL_COMMON_SOURCES3}
    )

target_link_libraries(lsutillib_common3 lsutillib_common2)

add_library(
    lsutillib_common4
    ${LSUTIL_COMMON_SOURCES4}
    )

target_link_libraries(lsutillib_common4 lsutillib_common3)

add_library(
    lsutillib_common5
    ${LSUTIL_COMMON_SOURCES5}
    )

target_link_libraries(lsutillib_common5 lsutillib_common4)

add_library(
    lsutillib_common6
    ${LSUTIL_COMMON_SOURCES6}
    )

target_link_libraries(lsutillib_common6 lsutillib_common5)

add_library(
    lsutillib_common7
    ${LSUTIL_COMMON_SOURCES7}
    )

target_link_libraries(lsutillib_common7 lsutillib_common6)

add_library(
    lsutillib_common8
    ${LSUTIL_COMMON_SOURCES8}
    )

target_link_libraries(lsutillib_common8 lsutillib_common7)

add_library(
    matrixolib
    ${LSUTIL_MATRIXO_SOURCES}
    ${LSUTIL_MATRIXO_C_SOURCES}
    )

target_link_libraries(matrixolib lsutillib_common8)

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
    ${LSUTIL_TYPEOP_SOURCES}
    )
add_dependencies(lsutiltypelib_common lsutillib_common8)
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
add_dependencies(fmmlib lsutillib_common1)
add_dependencies(fmmlib lsutillib_common2)
add_dependencies(fmmlib lsutillib_common3)
add_dependencies(fmmlib lsutillib_common4)
add_dependencies(fmmlib lsutillib_common5)
add_dependencies(fmmlib lsutillib_common6)
add_dependencies(fmmlib lsutillib_common7)
add_dependencies(fmmlib lsutillib_common8)
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
add_dependencies(lsintlib xcfun_interface)
add_dependencies(lsintlib pdpacklib)
add_dependencies(lsintlib lsutillib)
add_dependencies(lsintlib xcfun_interface)

# https://gitlab.com/dalton/IchorIntegralLibrary
include(IchorIntegralLibrary)

# https://gitlab.com/pett/tensor_lib
include(TensorLibrary)

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

if(ENABLE_REAL_SP)
   set(CCSDPT_SINGLE_PREC_SOURCE
       ${CMAKE_SOURCE_DIR}/LSDALTON/deccc/ccsdpt_kernels_sp.F90
       ${CMAKE_SOURCE_DIR}/LSDALTON/deccc/ccsdpt_full_sp.F90
       ${CMAKE_SOURCE_DIR}/LSDALTON/deccc/ccsdpt_dec_sp.F90
       )
endif()

if(ENABLE_REAL_SP)
   get_directory_property(LIST_OF_DEFINITIONS DIRECTORY ${CMAKE_SOURCE_DIR} COMPILE_DEFINITIONS)
   if(${CMAKE_SOURCE_DIR}/LSDALTON/deccc/ccsdpt_kernels.F90 IS_NEWER_THAN ${CMAKE_SOURCE_DIR}/LSDALTON/deccc/ccsdpt_kernels_sp.F90 OR
     ${CMAKE_SOURCE_DIR}/LSDALTON/deccc/ccsdpt_full.F90 IS_NEWER_THAN ${CMAKE_SOURCE_DIR}/LSDALTON/deccc/ccsdpt_full_sp.F90 OR
     ${CMAKE_SOURCE_DIR}/LSDALTON/deccc/ccsdpt_dec.F90 IS_NEWER_THAN ${CMAKE_SOURCE_DIR}/LSDALTON/deccc/ccsdpt_dec_sp.F90)
     add_custom_command(
     OUTPUT
     ${CCSDPT_SINGLE_PREC_SOURCE}
     COMMAND
     bash ${CMAKE_SOURCE_DIR}/LSDALTON/deccc/ccsdpt_sp.sh CMAKE_BUILD=${CMAKE_SOURCE_DIR}/LSDALTON/deccc ${LIST_OF_DEFINITIONS}
     DEPENDS
     ${CMAKE_SOURCE_DIR}/LSDALTON/deccc/ccsdpt_kernels.F90
     ${CMAKE_SOURCE_DIR}/LSDALTON/deccc/ccsdpt_full.F90
     ${CMAKE_SOURCE_DIR}/LSDALTON/deccc/ccsdpt_dec.F90
     ${CMAKE_SOURCE_DIR}/LSDALTON/deccc/ccsdpt_sp.sh
     )
   endif()
   unset(LIST_OF_DEFINITIONS)
endif()

if(ENABLE_DEC)
  if(ENABLE_REAL_SP)
    add_library(
      declib
      ${CCSDPT_SINGLE_PREC_SOURCE}
      ${DEC_SOURCES}
      )
  else()
    add_library(
      declib
      ${DEC_SOURCES}
      )
  endif()

  target_link_libraries(declib lsutiltypelib_common)
  target_link_libraries(declib lsutillib_common1)
  target_link_libraries(declib lsutillib_common2)
  target_link_libraries(declib lsutillib_common3)
  target_link_libraries(declib lsutillib_common4)
  target_link_libraries(declib lsutillib_common5)
  target_link_libraries(declib lsutillib_common6)
  target_link_libraries(declib lsutillib_common7)
  target_link_libraries(declib lsutillib_common8)
  target_link_libraries(declib lsintlib)
  target_link_libraries(declib linearslib)
endif()
  
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
if(ENABLE_DEC)
  target_link_libraries(lsdaltonmain declib)
endif()
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
    MERGE_STATIC_LIBS(
        lsint
        lsintlib
        )
endif()

set(LIBS_TO_MERGE
    lsutillib_precision
    cuda_gpu_interfaces
    matrixmlib
    lsutillib_common1
    lsutillib_common2
    lsutillib_common3
    lsutillib_common4
    lsutillib_common5
    lsutillib_common6
    lsutillib_common7
    lsutillib_common8
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
)

if(ENABLE_DEC)
  set(LIBS_TO_MERGE ${LIBS_TO_MERGE} declib)
endif()

set(LIBS_TO_MERGE ${LIBS_TO_MERGE} 
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
