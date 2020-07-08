!
!     cbiexc.h - Control common block for abacus/abaexc.F
!
      INTEGER, PARAMETER :: MAXPP = 200
      LOGICAL     :: SKIP, CUT, DIPSTR, ROTSTR, ROTVEL, FNAC, EXCTRP,
     &               SUMRUL, OOTV, MAGPRP, STOPPW
      INTEGER     :: LQ, LVEL, NEXCIT(8), MAXITE, MXNEXI, MXRM, 
     &               MXPHP, NABAPP, IPREXC, IPR1IN
      REAL*8      :: QMIN, QMAX, QSTEP, VMIN, VMAX, VSTEP, QINP,
     &               THREXC
      CHARACTER*8 :: LABAPP(MAXPP)
      INTEGER     :: LABSYM(MAXPP)
!
!
      COMMON /PPLBL / LABAPP, LABSYM
! LOGICAL
      COMMON /CBIEXC/ SKIP, CUT, DIPSTR, ROTSTR, ROTVEL, FNAC,
     &                EXCTRP, SUMRUL, OOTV, MAGPRP, STOPPW,
! INTEGER             
     &                LQ, LVEL, NEXCIT, MAXITE, MXNEXI, MXRM,
     &                MXPHP, NABAPP, IPREXC, IPR1IN,
! REAL
     &                QINP, THREXC, QMIN, QMAX, QSTEP, VMIN, VMAX,
     &                VSTEP
! -- end of abaexc.h --
