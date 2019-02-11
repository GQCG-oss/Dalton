file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/tools)

add_library(peter_utils_blocks "")

target_sources(
  peter_utils_blocks
  PRIVATE
    blocks.f90
  )

foreach(_util IN ITEMS aces2dalton xyz2dalton distances)
  add_executable(
    ${_util}
    ${_util}.f90
    )

  target_link_libraries(
    ${_util}
    PRIVATE
      peter_utils_blocks
    )

  set_target_properties(
    ${_util}
    PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tools"
  )
endforeach()

add_executable(FChk2HES FChk2HES.f)

set_target_properties(
  FChk2HES
  PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tools"
)

add_executable(labread labread.f)

set_target_properties(
  labread
  PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tools"
)
