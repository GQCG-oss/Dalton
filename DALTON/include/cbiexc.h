!
!     cbiexc.h - Control common block for abacus/abaexc.F
!
      LOGICAL         SKIP, CUT, DIPSTR, ROTSTR, ROTVEL, FNAC, EXCTRP,
CPFP
C     &        SUMRUL, OOTV
     &        SUMRUL, OOTV, MAGPRP
Cend-PFP
      PARAMETER       (MAXPP = 200)
      CHARACTER*8     LABAPP
      COMMON /PPLBL / LABAPP(MAXPP), LABSYM(MAXPP)
      COMMON /CBIEXC/ THREXC,
     &                NEXCIT(8), MAXITE, MXNEXI, MXRM,
     &                MXPHP, NABAPP, IPREXC, IPR1IN,
     &                SKIP, CUT, DIPSTR, ROTSTR, ROTVEL, FNAC, EXCTRP,
CPFP
C     &                SUMRUL, OOTV
     &                SUMRUL, OOTV, MAGPRP
Cend-PFP
! -- end of abaexc.h --
