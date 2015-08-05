option(ENABLE_ICHOR "Enable ichor integrals" OFF)

if(ENABLE_ICHOR)
    add_definitions(-DVAR_ICHOR)

    include(ExternalProject)

    set(_ichor_args
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
        )

    ExternalProject_Add(IchorIntegralLibrary
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/IchorIntegralLibrary
        BINARY_DIR ${PROJECT_BINARY_DIR}/external/IchorIntegralLibrary-build
        STAMP_DIR ${PROJECT_BINARY_DIR}/external/IchorIntegralLibrary-stamp
        TMP_DIR ${PROJECT_BINARY_DIR}/external/IchorIntegralLibrary-tmp
        INSTALL_DIR ${PROJECT_BINARY_DIR}/external
        CMAKE_ARGS ${_ichor_args}
        )

    unset(_ichor_args)

    include_directories(${PROJECT_BINARY_DIR}/external/IchorIntegralLibrary-build/modules)

    set(LSDALTON_EXTERNAL_LIBS
        ${PROJECT_BINARY_DIR}/external/lib/libichorlib.a
        ${LSDALTON_EXTERNAL_LIBS}
        )

    add_dependencies(lsintlib IchorIntegralLibrary)
endif()
