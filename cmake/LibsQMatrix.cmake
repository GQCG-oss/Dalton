# However, both DALTON and LSDALTON do not use HDF5 library.
# Therefore, I (Bin Gao) disable it for the time being.
#
# QMatrix uses HDF5 library to save the structures of matrices,
# modified from the CMakeLists.txt in the QMatrix library
#-INCLUDE(FindHDF5)
#-IF(HDF5_FOUND)
#-    # Location of the HDF5 includes
#-    INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})
#-    # Required compiler definitions for HDF5, the latest HDF5 APs are used, so it is necessary
#-    # to define H5_NO_DEPRECATED_SYMBOLS on some systems
#-    ADD_DEFINITIONS(${HDF5_DEFINITIONS})
#-    # The HDF5-1.8 APIs are used in the QMatrix library, so it is necessary to define
#-    # H5_NO_DEPRECATED_SYMBOLS on some systems
#-    ADD_DEFINITIONS(-DH5_NO_DEPRECATED_SYMBOLS)
#-    # Required libraries for all requested bindings are defined in ${HDF5_LIBRARIES_RELEASE}
#-    FOREACH(_LIB ${HDF5_LIBRARIES_RELEASE})
#-        STRING(FIND ${_LIB} "libpthread.a" _LIB_PTHREAD)
#-        IF(${_LIB_PTHREAD} EQUAL -1)
#-            SET(LIB_HDF5_RELEASE ${LIB_HDF5_RELEASE} ${_LIB})
#-        ELSE()
#-            # To prevent, for instance undefined reference to `__syscall_error'
#-            SET(LIB_HDF5_RELEASE ${LIB_HDF5_RELEASE} "-lpthread")
#-        ENDIF()
#-    ENDFOREACH()
#-    #ADD_LIBRARY(LIB_HDF5_RELEASE UNKNOWN IMPORTED)
#-    #SET_PROPERTY(TARGET LIB_HDF5_RELEASE PROPERTY IMPORTED_LOCATION "${HDF5_LIBRARIES_RELEASE}")
#-ELSE()
#-    MESSAGE(FATAL_ERROR "HDF5 library is needed for QMatrix!")
#-ENDIF()
# Sets the external QMatrix library, and a library libqmatrix.a will be generated
set(LIB_QMATRIX qmatrix)
set(ExternalProjectCMakeArgs
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    -DQMATRIX_ZERO_BASED=OFF
    -DQMATRIX_ROW_MAJOR=OFF
    -DQMATRIX_SINGLE_PRECISION=OFF
    -DQMATRIX_STORAGE_MODE=OFF
    -DQMATRIX_ENABLE_VIEW=ON
    -DQMATRIX_ENABLE_HDF5=OFF
    -DQMATRIX_3M_METHOD=ON
    -DQMATRIX_STRASSEN_METHOD=ON
    -DQMATRIX_AUTO_ERROR_EXIT=ON
    -DQMATRIX_TEST_EXECUTABLE=OFF
    -DQMATRIX_TEST_3M_METHOD=OFF
    -DLIB_QMATRIX_NAME=${LIB_QMATRIX}
    -DQMATRIX_BUILD_ADAPTER=ON
    -DEXTERNAL_BLOCK_MATRIX=OFF
    -DEXTERNAL_COMPLEX_MATRIX=OFF
    -DQMATRIX_ADAPTER_TYPE=F03
    -DQMATRIX_EXTERNAL_PATH=${CMAKE_Fortran_MODULE_DIRECTORY}
    -DQMATRIX_EXTERNAL_LIBRARIES=None
    -DLANG_F_MODULE=qmatrix_backend
    -DLANG_F_MATRIX=real_mat_t
    -DQMATRIX_Fortran_API=F03
    -DBLAS_LIBRARIES=None
    -DLAPACK_LIBRARIES=None
    -DPARENT_MODULE_DIR=${CMAKE_Fortran_MODULE_DIRECTORY}
   )
#-    -DLIB_HDF5_RELEASE=${LIB_HDF5_RELEASE}
add_external(${LIB_QMATRIX} qmatrix)
# LSDALTON/qmatrix/qmatrix_backend.F90 needs the header files of the QMatrix library
include_directories(${PROJECT_SOURCE_DIR}/external/qmatrix/include)
# Since QMatrix will use the matrix module of LSDALTON, we can not put it after
# matrixulib as in the external libraries.
#
# Instead, we put QMatrix library into LIB_LS_QMATRIX, which will be used
# for linking together with Gen1Int, QMatrix, OpenRSP and TDRSP libraries
set(LIB_LS_QMATRIX
    ${PROJECT_BINARY_DIR}/external/lib/lib${LIB_QMATRIX}.a)
#-    ${LIB_HDF5_RELEASE})
# QMatrix will use the matrix module of LSDALTON
#add_dependencies(${LIB_QMATRIX} matrixmlib)
#add_dependencies(${LIB_QMATRIX} matrixolib)
add_dependencies(${LIB_QMATRIX} matrixulib)
