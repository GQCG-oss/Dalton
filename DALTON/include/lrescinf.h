!
!     File: lrescinf.h
!     Purpose: Control of what to do in LRESC module
!
!     NOTE:
!
cx jim-gesc : RNLRSC login included in abainf.h, jimprt for debugging prints
      REAL*8 calfa, CFCZK, CSDZK, CFCBS, CSDBS, CPSOK, CLKIN, CDIAM,
     &       CDIAD, CDIAK, CANGP,CFCAV

       PARAMETER (calfa=1.0/137.036, CFCZK=1.0/3.0, CSDZK=-0.25D0,
     &            CFCBS=-0.25D0, CSDBS=-0.25D0, CPSOK=0.50D0,
     &            CLKIN=1.0D0, CDIAM=-1.0D0, CDIAD=-1.0D0,
     &            CDIAK=10.0/137.036, CANGP=-0.5D0,
     &            CFCAV=-1.0*7.0/16.0)

      LOGICAL SIGMAP1S, SIGMAP1T, SIGMAD1S, SIGMAD0S,
     &        SIGMAP3S, SIGMAP3T, LRESCALL, GAUCHANG,
     &        PRTALL

      INTEGER JIMPRT, LRATOM

      DOUBLE PRECISION LRFCAV(3,3), LRDIAK(3,3), LRANGP(3,3),
     &    LRDIAM(3,3), LRDIAD(3,3), LRLKIN(3,3), LRPSOK(3,3),
     &    LRPSKI(3,3), LRFCZK(3,3), LRSDZK(3,3), LRFCBS(3,3),
     &    LRSDBS(3,3), LRGAUG(3)

      COMMON /LRESCINF/ SIGMAP1S, SIGMAP1T, SIGMAD1S, SIGMAD0S,
     &    SIGMAP3S, SIGMAP3T, LRESCALL, GAUCHANG, JIMPRT, LRATOM,
     &    LRFCAV, LRDIAK, LRANGP, LRDIAM, LRDIAD, LRLKIN, LRPSOK,
     &    LRPSKI, LRFCZK, LRSDZK, LRFCBS, LRSDBS, LRGAUG, PRTALL
! -- end of lrescinf.h --
