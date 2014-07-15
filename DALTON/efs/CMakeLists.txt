set (EFS_HEADERS_DIR ${PROJECT_SOURCE_DIR}/external/efs)

include_directories(${EFS_HEADERS_DIR})

set(EFS_INTER_HEADERS ${EFS_HEADERS_DIR}/libquimera/libquimera.hpp ${EFS_HEADERS_DIR}/libechidna/libechidna.hpp ${EFS_HEADERS_DIR}/low/MPIwrap.hpp ${EFS_HEADERS_DIR}/math/tensors.hpp)

set_source_files_properties(${EFS_HEADERS_DIR}/libquimera/libquimera.hpp PROPERTIES GENERATED TRUE)
set_source_files_properties(${EFS_HEADERS_DIR}/libechidna/libechidna.hpp PROPERTIES GENERATED TRUE)
set_source_files_properties(${EFS_HEADERS_DIR}/low/MPIwrap.hpp           PROPERTIES GENERATED TRUE)
set_source_files_properties(${EFS_HEADERS_DIR}/math/tensors.hpp          PROPERTIES GENERATED TRUE)

if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    set(CMAKE_CXX_FLAGS "-fopenmp")
else()
    set(CMAKE_CXX_FLAGS "-march=native -fopenmp")
endif()

add_library(
    efs_interface
    efs_interface.cpp
    ${EFS_INTER_HEADERS}
    )

add_dependencies(efs_interface efs)
