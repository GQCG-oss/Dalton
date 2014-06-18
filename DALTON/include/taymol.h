! FILE: taymol.h
      REAL*8  ERGMOL
      COMMON /CB_ERGMOL/ ERGMOL
! May 2014 hjaaj: GRDMOL and HESMOL removed to reduce static memory,
!     they are now saved and read from ABACUS.RESTART using ABAWRIT_TAYMOL
!     and ABAREAD_TAYMOL.
!old  REAL*8  ERGMOL, GRDMOL, HESMOL
!old  COMMON /TAYMOL/ ERGMOL, GRDMOL(MXCOOR), HESMOL(MXCOOR,MXCOOR)
! --- end of taymol.h ---
