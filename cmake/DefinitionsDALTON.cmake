if(DEVELOPMENT_CODE AND NOT ENABLE_RELEASE)
    add_definitions(-DMOD_UNRELEASED)
endif()

if(ENABLE_GEN1INT)
    add_definitions(-DBUILD_GEN1INT)
endif()

add_definitions(-DVAR_MFDS)
add_definitions(-D_FILE_OFFSET_BITS=64)
add_definitions(-DIMPLICIT_NONE)

if(ENABLE_64BIT_INTEGERS)
    add_definitions(-DVAR_INT64)
endif()

# forward CPP directly to the code
set(CPP)
if(NOT "${CPP}" STREQUAL "")
    add_definitions(${CPP})
endif()
