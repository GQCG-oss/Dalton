!
! variables required in a Multipole Method (FMM) run
! external to the main MM code
! (e.g. multipole moment generation order via HERMIT)
!
      LOGICAL         MMSKIP, MMCUT,MMSKIPDF, MMCUTDF, DOSAVEMM
      LOGICAL         USEBUFMM
      INTEGER         IPRINT, LMAX, MMPRDF, LMAXDF,                                & 
     &                MMSAVE, MMBUFLEN, MAXBUFI, MAXBUFR, MAXBUFN
!
      COMMON /CBIFMM/ IPRINT, LMAX, MMSKIP, MMCUT, MMSAVE, DOSAVEMM,               &
     &                USEBUFMM, MMBUFLEN, MAXBUFI, MAXBUFR, MAXBUFN
!
!     defaults:
!
      PARAMETER ( MMSKIPDF=.FALSE.,MMCUTDF=.TRUE.)
      PARAMETER ( MMPRDF = 0, LMAXDF = 8)

