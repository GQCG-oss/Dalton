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

configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/CTestCustom.cmake.in
    ${CMAKE_BINARY_DIR}/CTestCustom.cmake
    @ONLY
    )
configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/ConfigPreTest.cmake.in
    ${CMAKE_BINARY_DIR}/ConfigPreTest.cmake
    @ONLY
    )

include(CTest)
enable_testing()
