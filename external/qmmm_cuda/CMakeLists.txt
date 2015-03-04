project(qmmm_cuda)

cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/cmake
    ${CMAKE_BINARY_DIR}
    )

include(mergestaticlibs)

find_package(CUDA)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-arch=sm_20)
cuda_add_library(integrals_1el_potential_cuda src/integrals_1el_potential_cuda_lib.cu)

add_library(integrals_1el_potential src/integrals_1el_potential_lib.cc)
add_library(qmmm_utility            src/qmmm_utility.cc)

MERGE_STATIC_LIBS(
    qmmm_cuda
    integrals_1el_potential_cuda
    integrals_1el_potential
    qmmm_utility
    )
