option(ENABLE_OPENMP "Enable OpenMP parallelization" OFF)

if(ENABLE_OPENMP)
    add_definitions(-DVAR_OMP)
    set(ENABLE_THREADED_MKL TRUE)
else()
    set(ENABLE_THREADED_MKL FALSE)
endif()
