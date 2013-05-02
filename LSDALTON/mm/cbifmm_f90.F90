Module cbifmm
! A FUCKING COMMON BLOCK 
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
end module
