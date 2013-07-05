# set cdash buildname
set(BUILDNAME
    "${CMAKE_SYSTEM_NAME}-${CMAKE_HOST_SYSTEM_PROCESSOR}-${CMAKE_Fortran_COMPILER_ID}-${BLAS_TYPE}-${CMAKE_BUILD_TYPE}"
    CACHE STRING
    "Name of build on the dashboard"
    )

# set ctest own timeout
set(DART_TESTING_TIMEOUT
    "1200"
    CACHE STRING
    "Set timeout in seconds for every single test"
    )

include(TestsCommon)
include(TestsDALTON)
include(TestsLSDALTON)

include(CTest)
enable_testing()
