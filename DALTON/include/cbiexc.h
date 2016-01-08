!
!     cbiexc.h - Control common block for abacus/abaexc.F
!
      LOGICAL         SKIP, CUT, DIPSTR, ROTSTR, ROTVEL, FNAC, EXCTRP,
CPFP
C     &        SUMRUL, OOTV
     &        SUMRUL, OOTV, MAGPRP,
Cend-PFP
CClark:7/1/2016
     &        STOPPW
      REAL*8          QMIN,QMAX,QSTEP
CClark:end
      PARAMETER       (MAXPP = 200)
      CHARACTER*8     LABAPP
      COMMON /PPLBL / LABAPP(MAXPP), LABSYM(MAXPP)
      COMMON /CBIEXC/ THREXC,
CClark:7/1/2016
     &                QMIN,QMAX,QSTEP,
CClark:end
     &                NEXCIT(8), MAXITE, MXNEXI, MXRM,
     &                MXPHP, NABAPP, IPREXC, IPR1IN,
     &                SKIP, CUT, DIPSTR, ROTSTR, ROTVEL, FNAC, EXCTRP,
CPFP
C     &                SUMRUL, OOTV
     &                SUMRUL, OOTV, MAGPRP,
Cend-PFP
CClark:7/1/2016
     &                STOPPW
CClark:end
! -- end of abaexc.h --
