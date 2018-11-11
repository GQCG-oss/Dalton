if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    add_definitions(-DSYS_LINUX)
    add_definitions(-DSYS_UNIX)
    if(${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "i686")
        add_definitions(-DARCH32BIT)
    endif()
endif()

if(${CMAKE_SYSTEM_NAME} STREQUAL "FreeBSD")
    add_definitions(-DSYS_FREEBSD)
    add_definitions(-DSYS_UNIX)
    if(${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "i386")
        add_definitions(-DARCH32BIT)
    endif()
endif()

if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    add_definitions(-DSYS_DARWIN)
    add_definitions(-DSYS_UNIX)
    # fixme: HAVE_NO_LSEEK64 should be tested by cmake
    #        not hardcoded for Mac OSX
    add_definitions(-DHAVE_NO_LSEEK64)
    # amt: This definition is a workaround for some missing strdup() used by DFT
    #      c code on OSX Lion and newer. The workaround also works on older OSX but
    #      not needed there. For now live with always on.

    # work-around for error in Macports cmake on OSX
    if(${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "i386")
        set(CMAKE_HOST_SYSTEM_PROCESSOR x86_64)
    endif()
endif()

if(${CMAKE_SYSTEM_NAME} STREQUAL "CYGWIN")
    add_definitions(-DSYS_CYGWIN)
    add_definitions(-DSYS_LINUX)
    add_definitions(-DSYS_UNIX)
    # fixme: HAVE_NO_LSEEK64 should be tested by cmake
    #        now just trying to get it compiled
    # add_definitions(-DHAVE_NO_LSEEK64)
endif()

if(${CMAKE_SYSTEM_NAME} STREQUAL "AIX")
    add_definitions(-DSYS_AIX)
    add_definitions(-DSYS_UNIX)
endif()
