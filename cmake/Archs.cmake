message("-- Detected CMAKE_SYSTEM_NAME is \"${CMAKE_SYSTEM_NAME}\"")
set(LINUX_UBUNTU FALSE)
set(LINUX_FEDORA FALSE)

if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    add_definitions(-DSYS_LINUX)
    if(${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "i686")
        add_definitions(-DARCH32BIT)
    elseif(${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "x86_64")
    else()
    endif()

    execute_process(
        COMMAND lsb_release -d
        TIMEOUT 1
        OUTPUT_VARIABLE LINUX_DISTRIBUTION
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/cmake
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(${LINUX_DISTRIBUTION} MATCHES "Ubuntu")
        set(LINUX_UBUNTU TRUE)
    else()
    endif()
    if(${LINUX_DISTRIBUTION} MATCHES "Fedora")
        set(LINUX_FEDORA TRUE)
    else()
    endif()
endif()

if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    add_definitions(-DSYS_LINUX)
    # fixme: HAVE_NO_LSEEK64 should be tested by cmake, not hardcoded for Mac OSX
    #        now just trying to get it compiled
    add_definitions(-DHAVE_NO_LSEEK64)
endif()

if(${CMAKE_SYSTEM_NAME} STREQUAL "AIX")
    add_definitions(-DSYS_AIX)
endif()

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
     add_definitions(-DSYS_WINDOWS)
endif()

