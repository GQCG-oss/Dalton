set(CUSTOM_PCMSolver_LOCATION ${PROJECT_BINARY_DIR}/external/pcmsolver)
get_filename_component(ZLIB_ROOT ${ZLIB_LIBRARIES} PATH)

# PCMSolver does not know profile
if(CMAKE_BUILD_TYPE MATCHES "profile")
  set(PCM_BUILD_TYPE "release")
else()
  set(PCM_BUILD_TYPE ${CMAKE_BUILD_TYPE})
endif()

list(APPEND PCMSolverCMakeArgs
  -DCMAKE_BUILD_TYPE=${PCM_BUILD_TYPE}
  -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
  -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
  -DEXTRA_Fortran_FLAGS=${PCM_EXTRA_Fortran_FLAGS}
  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
  -DEXTRA_C_FLAGS=${PCM_EXTRA_C_FLAGS}
  -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
  -DEXTRA_CXX_FLAGS=${PCM_EXTRA_CXX_FLAGS}
  -DENABLE_CXX11_SUPPORT=OFF
  -DBOOST_INCLUDEDIR=${BOOST_INCLUDE_DIRS}
  -DBOOST_LIBRARYDIR=${BOOST_LIBRARIES}
  -DENABLE_64BIT_INTEGERS=${ENABLE_64BIT_INTEGERS}
  -DENABLE_TESTS=OFF
  -DENABLE_LOGGER=OFF
  -DENABLE_TIMER=OFF
  -DBUILD_STANDALONE=OFF
  -DENABLE_FORTRAN_API=ON
  -DSTATIC_LIBRARY_ONLY=ON
  -DENABLE_GENERIC=${ENABLE_STATIC_LINKING}
  -DZLIB_ROOT=${ZLIB_ROOT}
  -DPYTHON_INTERPRETER=${PYTHON_EXECUTABLE}
  -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
  -DCMAKE_INSTALL_LIBDIR=lib
  )

ExternalProject_Add(pcmsolver
  PREFIX ${CUSTOM_PCMSolver_LOCATION}
  GIT_REPOSITORY https://github.com/PCMSolver/pcmsolver
  GIT_TAG v1.1.0
  CMAKE_ARGS "${PCMSolverCMakeArgs}"
  INSTALL_DIR "${CUSTOM_PCMSolver_LOCATION}/install"
  )

# Set also variables usually set by find_package
ExternalProject_Get_Property(pcmsolver INSTALL_DIR)
set(PCMSolver_LIBRARY "${INSTALL_DIR}/lib/libpcm.a")
file(MAKE_DIRECTORY ${INSTALL_DIR}/include/pcmsolver)  # note [1] below
set(PCMSolver_INCLUDE_DIRS "${INSTALL_DIR}/include" ${ZLIB_INCLUDE_DIRS})
set(PCMSolver_LIBRARIES ${PCMSolver_LIBRARY} ${ZLIB_LIBRARIES})

# Set target for DALTON to depend upon as set by find_package
add_library(PCMSolver::PCMSolver STATIC IMPORTED GLOBAL)
add_dependencies(PCMSolver::PCMSolver pcmsolver)
set_target_properties(PCMSolver::PCMSolver PROPERTIES
  IMPORTED_LOCATION "${PCMSolver_LIBRARY}"
  INTERFACE_LINK_LIBRARIES "${PCMSolver_LIBRARIES}"
  INTERFACE_INCLUDE_DIRECTORIES "${PCMSolver_INCLUDE_DIRS}"
  )

include_directories(SYSTEM "${PCMSolver_INCLUDE_DIRS}")

link_directories(${INSTALL_DIR}/lib)
message("${PCMSolver_LIBRARIES}")
set(PCMSOLVER_LIBS
  ${CMAKE_CXX_IMPLICIT_LINK_LIBRARIES}
  ${PCMSolver_LIBRARIES}
  ${CMAKE_CXX_IMPLICIT_LINK_LIBRARIES}
  ${ZLIB_LIBRARIES}
  PARENT_SCOPE)
set(PCMSolver_PARSE_DIR ${INSTALL_DIR}/bin PARENT_SCOPE)
add_definitions(-DHAS_PCMSOLVER)

# Set up sources DALTON-side
set(pcm_sources
  ${PROJECT_SOURCE_DIR}/DALTON/pcm/pcm_config.F90
  ${PROJECT_SOURCE_DIR}/DALTON/pcm/pcm_write.F90
  ${PROJECT_SOURCE_DIR}/DALTON/pcm/pcm_integrals.F90
  ${PROJECT_SOURCE_DIR}/DALTON/pcm/pcm_scf.F90
  ${PROJECT_SOURCE_DIR}/DALTON/pcm/pcm_linrsp.F90
  )
if(ENABLE_MPI)
  set(pcm_sources
    ${pcm_sources}
    ${PROJECT_SOURCE_DIR}/DALTON/pcm/pcm_parallel.F90
    )
endif()
set(DALTON_FREE_FORTRAN_SOURCES ${DALTON_FREE_FORTRAN_SOURCES} ${pcm_sources} PARENT_SCOPE)

# [1] It's nice to have a full PCMSolver::PCMSolver target that has embedded
# the library, the linking library paths with dependencies, and the include
# paths with dependencies, just like FindPCMSolver supplies.
# Problem is that conventional target
# derived from ExternalProject_Add is a highly sought but not quite
# certified cmake pattern. Hence INTERFACE_INCLUDE_DIRECTORIES complains
# that the directories don't exist at configure time. Hence the hack to
# create an empty directory.
