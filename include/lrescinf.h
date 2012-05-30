!
!     File: lrescinf.h
!     Purpose: Control of what to do in LRESC module
!
!     NOTE: 
!
cx jim-gesc : RNLRSC login included in abainf.h, jimprt for debugging prints
      INTEGER NATOM, JIMPRT
!     PARAMETER (NSYMLasc = 8)
      REAL*8 LRGAUG(3)
      LOGICAL SIGMAP1S, SIGMAP1T, SIGMAD1S, SIGMAD0S,                   &
     &        SIGMAP3S, SIGMAP3T, LRESCALL,GAUCHANG
      COMMON /LRESCINF/ SIGMAP1S, SIGMAP1T, SIGMAD1S, SIGMAD0S,         &
     &        GAUCHANG, SIGMAP3S, SIGMAP3T, JIMPRT, NATOM, LRGAUG,      &
     &        LRESCALL

! -- end of abainf.h --
