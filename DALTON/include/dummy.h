! --- dummy.h ---

!     Make DUMMY, VDUMMY, IDUMMY  write protected (on most computers)
!     Used in parameter lists for parameters which are not supposed
!     to be used by the subroutine/function in this call.
!     VDUMMY, IVDUMMY is for arrays, so we can quench compiler warnings.

      REAL*8,  parameter ::  DUMMY = 1.0D20,    VDUMMY(1:2) = 1.0D19
      INTEGER, parameter :: IDUMMY = -9999999, IVDUMMY(1:2) = -9999998

! --- end of dummy.h ---
