find_package(Sphinx)
if(SPHINX_FOUND)
    add_custom_target(
        html
        COMMAND sphinx-build -b html -d _build/doctrees ${CMAKE_SOURCE_DIR}/doc _build/html
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        )
    add_custom_target(
        latex
        COMMAND sphinx-build -b latex -d _build/doctrees ${CMAKE_SOURCE_DIR}/doc _build/latex
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        )
else()
    add_custom_target(
        html
        COMMAND echo error: please install python-sphinx and python-matplotlib first
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        )
    add_custom_target(
        latex
        COMMAND echo error: please install python-sphinx and python-matplotlib first
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        )
endif()
