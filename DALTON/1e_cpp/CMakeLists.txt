set(CMAKE_CXX_FLAGS "-march=native -fopenmp")

add_library(
    vpotdamp
    vpotdamp.cpp
    ababps.cpp
    )
