cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(cuda_interface)

enable_language(CXX)

find_package(CUDA)

if(CUDA_FOUND)
    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-gencode arch=compute_35,code=sm_35; -ccbin gcc)

    cuda_compile(
        cuda_objects
        cuda_api.cu
        )

    cuda_add_library(cuda_interface ${cuda_objects})
endif()
