# Set variables to help find Python library that is compatible with interpreter
# This is needed since the FindPythonLibs standard CMake module tends to find
# a version of the libraries different from the interpreter.
# Here we are just setting up "hint" variables, based on executing simple scripts
# with the Python interpreter, before actually running the
# CMake builtin module. This makes sure that we find the EXACT version of
# the libraries matching that of the interpreter.
# Copied from https://bitbucket.org/fenics-project/dolfin
function(set_python_hint_variables pyLibsFound pyLibs pyInclude pyInclude2) 
	find_package(PythonInterp REQUIRED)
	if (PYTHONINTERP_FOUND)
		# Get Python include path from Python interpretter
  		execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c
                  "import distutils.sysconfig, sys; sys.stdout.write(distutils.sysconfig.get_python_inc())"
                  OUTPUT_VARIABLE _PYTHON_INCLUDE_PATH
                  RESULT_VARIABLE _PYTHON_INCLUDE_RESULT)
                # Get Python library path from interpreter                                                                                      
                execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c
                                "import os, sys, inspect; sys.stdout.write(os.path.split(os.path.split(inspect.getfile(inspect))[0])[0])"
                                OUTPUT_VARIABLE _PYTHON_LIB_PATH
                                RESULT_VARIABLE _PYTHON_LIB_RESULT)
                # Set include path, if returned by interpreter
                if ("${_PYTHON_INCLUDE_RESULT}" STREQUAL "0")
                  set(PYTHON_INCLUDE_DIR ${_PYTHON_INCLUDE_PATH})
                endif()
                # Add a search path for Python library based on output from
                # interpreter
                set(CMAKE_LIBRARY_PATH_SAVE ${CMAKE_LIBRARY_PATH})
                if ("${_PYTHON_LIB_RESULT}" STREQUAL "0")
                  set(CMAKE_LIBRARY_PATH ${_PYTHON_LIB_PATH})
                endif()
                # Find Pythons libs
                find_package(PythonLibs ${PYTHON_VERSION_STRING} EXACT REQUIRED)
		
		set(pyLibsFound "${PYTHONLIBS_FOUND}" PARENT_SCOPE)
		
		set(pyLibs      "${PYTHON_LIBRARIES}" PARENT_SCOPE)
		
		list(GET PYTHON_INCLUDE_DIRS -1 pyInclude)
		set(pyInclude "${pyInclude}" PARENT_SCOPE)
		list(GET PYTHON_INCLUDE_DIRS  0 pyInclude2)
		set(pyInclude2 "${pyInclude2}" PARENT_SCOPE)

		set(PYTHON_LIBRARIES "${PYTHON_LIBRARIES}" PARENT_SCOPE)
        endif()
endfunction()
