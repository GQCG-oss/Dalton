add_custom_target(python_fcomp
    COMMAND python ${CMAKE_SOURCE_DIR}/cmake/get_compiler_version.py ${CMAKE_Fortran_COMPILER}
    OUTPUT_VARIABLE FORTRAN_COMPILER_VERSION
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Fetching current fortran compiler version"
)

add_custom_target(python_ccomp
    COMMAND python ${CMAKE_SOURCE_DIR}/cmake/get_compiler_version.py ${CMAKE_C_COMPILER}
    OUTPUT_VARIABLE C_COMPILER_VERSION
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Fetching current C compiler version"
)

configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/compilation_info.py.in
    ${CMAKE_BINARY_DIR}/compilation_info.py
    )

add_custom_target(generate_compinfo
    COMMAND python "${CMAKE_BINARY_DIR}/compilation_info.py" > "${CMAKE_SOURCE_DIR}/lsutil/compilation_info.f90" 
    COMMENT "Generating lsutil/compilation_info.f90"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
#    DEPENDS python_fcomp python_ccomp
)
add_dependencies(generate_compinfo python_fcomp python_ccomp)

add_custom_target(python_rev
    COMMAND python ${CMAKE_SOURCE_DIR}/cmake/get_last_rev.py ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE LAST_GIT_REVISION
    COMMENT "Fetching latest git revision number"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/git_revision_info.py.in
    ${CMAKE_BINARY_DIR}/git_revision_info.py
    )

add_custom_target(generate_revinfo
    COMMAND python "${CMAKE_BINARY_DIR}/git_revision_info.py" > "${CMAKE_SOURCE_DIR}/lsutil/git_revision_info.f90" 
    COMMENT "Generating lsutil/git_revision_info.f90"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
add_dependencies(generate_revinfo python_rev)
