!#define LSLIB_RESTART
PROGRAM lslib_test
use lslib_tester_mod
implicit none
logical :: OnMaster,meminfo_slaves
Integer :: lupri,luerr
luerr          = 0
lupri          = 0
! setup the calculation 
call lslib_init(OnMaster,lupri,luerr)

IF(OnMaster) call LSlib_test_driver(OnMaster,lupri,luerr,meminfo_slaves)

! free everything take time and close the files
call lslib_free(OnMaster,lupri,luerr,meminfo_slaves)

END PROGRAM LSLIB_TEST
