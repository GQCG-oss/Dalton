
message( "TESTING TESTING TESTING" )

execute_process(
   COMMAND python ${PROJECT_SOURCE_DIR}/LSDALTON/tools/check_decbcast.py
   OUTPUT_VARIABLE _output
   OUTPUT_STRIP_TRAILING_WHITESPACE
   )

if( NOT ${_output} MATCHES "TESTSTATUS: GOOD")
   message(FATAL_ERROR "LSDALTON BCAST TYPE TEST NOT PASSED, PLEASE RERUN LSDALTON/tools/check_decbcast.py TO GET MORE INFORMATION")
else()
   message( "FOO" )
endif()

