
set(CTEST_NIGHTLY_START_TIME "00:00:00 CEST")
set(CTEST_DROP_METHOD        "http")
set(CTEST_DROP_SITE          "repo.ctcc.no")
set(CTEST_DROP_SITE_CDASH    TRUE)

if(DEFINED ENV{CTEST_PROJECT_NAME})
    set(CTEST_PROJECT_NAME  "$ENV{CTEST_PROJECT_NAME}")
    set(CTEST_DROP_LOCATION "/CDash/submit.php?project=$ENV{CTEST_PROJECT_NAME}")
else()
    set(CTEST_PROJECT_NAME  "DALTON")
    set(CTEST_DROP_LOCATION "/CDash/submit.php?project=DALTON")
endif()

if(DEFINED ENV{CTEST_MAKE_NUM_PROCS})
    set(MAKECOMMAND "make -j$ENV{CTEST_MAKE_NUM_PROCS}" CACHE STRING "Custom make command")
endif()

# total allowed time for all tests in seconds
# see http://www.itk.org/Wiki/CTest:Using_CTEST_and_CDASH_without_CMAKE#General_CTest_settings
set(CTEST_TIMEOUT           "21600")
