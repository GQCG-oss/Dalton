# only in the development code we get the hash from git
# in the released code we read it from file, in this case
# it is already set at this stage
if(DEVELOPMENT_CODE)
    find_package(Git)
    if(GIT_FOUND)
        execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-list --max-count=1 HEAD
            OUTPUT_VARIABLE GIT_REVISION
            ERROR_QUIET
            )
        if(NOT ${GIT_REVISION} STREQUAL "")
            string(STRIP ${GIT_REVISION} GIT_REVISION)
        endif()

        execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
            OUTPUT_VARIABLE GIT_BRANCH
            ERROR_QUIET
            )
        if(NOT ${GIT_BRANCH} STREQUAL "")
            string(STRIP ${GIT_BRANCH} GIT_BRANCH)
        endif()
    endif()
endif()
