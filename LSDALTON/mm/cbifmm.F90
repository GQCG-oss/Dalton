Module cbifmm_module
! A COMMON BLOCK 
!
! variables required in a Multipole Method (FMM) run
! external to the main MM code
! (e.g. multipole moment generation order via HERMIT)
!
LOGICAL :: MMSKIP, MMCUT,DOSAVEMM
LOGICAL :: USEBUFMM
INTEGER :: IPRINT, LMAX
INTEGER :: MMSAVE, MMBUFLEN, MAXBUFI, MAXBUFR, MAXBUFN
!
!     defaults:
!
LOGICAL,PARAMETER :: MMSKIPDF=.FALSE.,MMCUTDF=.TRUE.
INTEGER,PARAMETER :: MMPRDF = 0, LMAXDF = 8
contains 
   subroutine cbifmm_dummy
      implicit none
   end subroutine cbifmm_dummy
end module cbifmm_module

