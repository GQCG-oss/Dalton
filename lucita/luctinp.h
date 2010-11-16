!     Common block for LUCITA input information
!     read in dirrdn.F/PSIINP/LUCTINP
!     Information is transferred to LUCITA routine dirluct.F
!     for further processing.
!
      INTEGER NTABLE
      PARAMETER (NTABLE = 30)

      INTEGER MXNGAS
      PARAMETER (MXNGAS = 16)

      CHARACTER*72 WAFFCD, CALCTP, SZCALD, TITLUC, CRDINA,              &
     &             CRDGAS(MXNGAS), CRDGOC(MXNGAS), CRDFRO, CRDRS1,      &
     &             CRDRS2, CRDRS3

      COMMON /LUCTINFC/ CRDINA, CRDGAS, CRDGOC, CRDFRO, CRDRS1, CRDRS2, &
     &                  CRDRS3, WAFFCD, CALCTP, SZCALD, TITLUC

      INTEGER IMOKW, NROOTD, ISSYMD, NACTED, IMULTD, IPRNGD, IPRNLD,    &
     &        IDENSD, MXHL1D, MXEL3D, INGASD, NSEQCD, IRSTLT, MXCIVE,   &
     &        IPARMODEL, IPARIOMOD, ICIMAXITER, IMAXBLKSIZE, IFILESYS,  &
     &        I_USE_DIST_ROUTE, IFRACTFILE
      COMMON /LUCTINFI/ IMOKW, NROOTD, ISSYMD, NACTED, IMULTD, IPRNGD,  &
     &                  IPRNLD, IDENSD, MXHL1D, MXEL3D, INGASD, NSEQCD, &
     &                  IRSTLT, MXCIVE,IPARMODEL,IPARIOMOD,ICIMAXITER,  &
     &                  IMAXBLKSIZE, IFILESYS, I_USE_DIST_ROUTE,        &
     &                  IFRACTFILE
      REAL*8 CTRUNC_FAC
	COMMON /LUCTINFR/ CTRUNC_FAC
