# set cdash buildname
set(BUILDNAME
    "BUILDNAME-not-set"
    CACHE STRING
    "Name of build on the dashboard"
    )
# set ctest own timeout
set(DART_TESTING_TIMEOUT
    "1200"
    CACHE STRING
    "Set timeout in seconds for every single test"
    )
include(Tests)
include(CTest)
enable_testing()
