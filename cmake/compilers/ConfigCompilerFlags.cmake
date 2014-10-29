include(SaveCompilerFlags)

# check whether we can use -xHost on this machine
include(ConfigXHostFlag)

if(CMAKE_C_COMPILER_WORKS)
    include(CFlags)
endif()

if(CMAKE_CXX_COMPILER_WORKS)
    include(CXXFlags)
endif()

if(CMAKE_Fortran_COMPILER_WORKS)
    include(FortranFlags)
endif()
