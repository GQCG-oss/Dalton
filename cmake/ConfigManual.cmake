add_custom_target(
    pdfmanual
    COMMAND DALTON_VERSION=${CMAKE_BINARY_DIR}/../VERSION DOC_DIRECTORY=${CMAKE_BINARY_DIR}/Doc/ ${CMAKE_BINARY_DIR}/Doc/makePDFmanual.sh
    COMMAND mv ${CMAKE_BINARY_DIR}/Doc/*.pdf ${CMAKE_BINARY_DIR}
    )

add_custom_target(
    htmlmanual
    COMMAND DALTON_VERSION=${CMAKE_BINARY_DIR}/../VERSION DOC_DIRECTORY=${CMAKE_BINARY_DIR}/Doc/ ${CMAKE_BINARY_DIR}/Doc/makeHTMLmanual.sh
    )
