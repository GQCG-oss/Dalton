!> @file
!> Contains variables for general library routines for integral evaluation
MODULE lslib_state
use configurationType, only: configitem
use typedeftype, only: lsitem
use precision, only: realk
type(lsitem),save     :: ls
Integer,save          :: nbasis
type(configItem),save :: LSlibconfig
logical,save          :: state_set = .FALSE.
real(realk),save      :: tstart,tend
END MODULE lslib_state

