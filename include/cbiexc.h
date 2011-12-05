!
!     cbiexc.h - Control common block for abacus/abaexc.F
!
      LOGICAL SKIP, CUT, DIPSTR, ROTSTR, ROTVEL, FNAC, EXCTRP,
!PFP
     &        SUMRUL, OOTV
!     &        SUMRUL, OOTV, DIRECT, MAGPRP
!     &        SUMRUL, OOTV, MAGPRP
!end-PFP
      PARAMETER       (MAXPP = 200)
      CHARACTER*8     LABAPP
      COMMON /PPLBL / LABAPP(MAXPP), LABSYM(MAXPP)
      COMMON /CBIEXC/ THREXC,
     &                NEXCIT(8), MAXITE, MXNEXI, MXRM,
     &                MXPHP, NABAPP, IPREXC, IPR1IN,
     &                SKIP, CUT, DIPSTR, ROTSTR, ROTVEL, FNAC, EXCTRP,
!PFP
     &                SUMRUL, OOTV
!     &                SUMRUL, OOTV, DIRECT, MAGPRP
!     &                SUMRUL, OOTV, MAGPRP
!end-PFP
! -- end of abaexc.h --
