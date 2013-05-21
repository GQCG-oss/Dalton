if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    add_definitions(-DSYS_LINUX)
    if(${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "i686")
        add_definitions(-DARCH32BIT)
    endif()
endif()

if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    add_definitions(-DSYS_LINUX)
    # fixme: HAVE_NO_LSEEK64 should be tested by cmake
    #        not hardcoded for Mac OSX
    add_definitions(-DHAVE_NO_LSEEK64)
    # amt: This definition is a workaround for some missing strdup() used by DFT
    #      c code on OSX Lion and newer. The workaround also works on older OSX but
    #      not needed there. For now live with always on.
    add_definitions(-DSYS_OSXLION)
endif()

if(${CMAKE_SYSTEM_NAME} STREQUAL "AIX")
    add_definitions(-DSYS_AIX)
endif()

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    add_definitions(-DSYS_WINDOWS)
endif()
