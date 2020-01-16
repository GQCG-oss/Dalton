    find_package(Git)
    if(GIT_FOUND)
        execute_process(
            COMMAND           ${GIT_EXECUTABLE} rev-list --max-count=1 HEAD
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            OUTPUT_VARIABLE   GIT_REVISION
            ERROR_QUIET
            )
        if(NOT ${GIT_REVISION} STREQUAL "")
            string(STRIP ${GIT_REVISION} GIT_REVISION)
        endif()

      # radovan: deactivated 2014-03-26
      #          not 100% portable and fails on some machines
      # execute_process(
      #     COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
      #     OUTPUT_VARIABLE GIT_BRANCH
      #     ERROR_QUIET
      #     )
      execute_process(
          COMMAND           ${GIT_EXECUTABLE} branch
          COMMAND           grep "*"
          WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
          OUTPUT_VARIABLE   GIT_BRANCH
          ERROR_QUIET
          )
      if(NOT ${GIT_BRANCH} STREQUAL "")
          string(REPLACE "*" " " GIT_BRANCH ${GIT_BRANCH})
          string(STRIP ${GIT_BRANCH} GIT_BRANCH)
           message("-- GIT_BRANCH            : ${GIT_BRANCH}")
      endif()

      execute_process(
          COMMAND           ${GIT_EXECUTABLE} status -uno
          WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
          OUTPUT_FILE       ${PROJECT_BINARY_DIR}/GIT_STATUS_AT_BUILD
          ERROR_QUIET
          )

    endif()
