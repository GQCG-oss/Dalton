execute_process(
   COMMAND python ${CMAKE_SOURCE_DIR}/LSDALTON/tools/check_decbcast.py
   OUTPUT_VARIABLE _output
   OUTPUT_STRIP_TRAILING_WHITESPACE
   )

if(NOT ${_output} MATCHES "TESTSTATUS: GOOD")
   message( "LSDALTON BCAST TYPE TEST RETURNED:" )
   message( ${_output} )
   message(FATAL_ERROR "LSDALTON BCAST TYPE TEST NOT PASSED" )
endif()
