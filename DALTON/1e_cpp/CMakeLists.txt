if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    set(CMAKE_CXX_FLAGS "-fopenmp")
else()
    set(CMAKE_CXX_FLAGS "-march=native -fopenmp")
endif()

add_library(
    vpotdamp
    vpotdamp.cpp
    ababps.cpp
    )
