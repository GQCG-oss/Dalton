!
!     File: lrescinf.h
!     Purpose: Control of what to do in LRESC module
!
!     NOTE: 
!
cx jim-gesc : RNLRSC login included in abainf.h, jimprt for debugging prints
!     PARAMETER (NSYMLasc = 8)

      LOGICAL SIGMAP1S, SIGMAP1T, SIGMAD1S, SIGMAD0S,
     &        SIGMAP3S, SIGMAP3T, LRESCALL, GAUCHANG

      INTEGER JIMPRT, LRATOM

      DOUBLE PRECISION LRFCAV(3,3), LRDIAK(3,3), LRANGP(3,3),
     &    LRDIAM(3,3), LRDIAD(3,3), LRLKIN(3,3), LRPSOK(3,3),
     &    LRPSKI(3,3), LRFCZK(3,3), LRSDZK(3,3), LRFCBS(3,3), 
     &    LRSDBS(3,3), LRGAUG(3)

      COMMON /LRESCINF/ SIGMAP1S, SIGMAP1T, SIGMAD1S, SIGMAD0S,
     &    SIGMAP3S, SIGMAP3T, LRESCALL, GAUCHANG, JIMPRT, LRATOM, 
     &    LRFCAV, LRDIAK, LRANGP, LRDIAM, LRDIAD, LRLKIN, LRPSOK,
     &    LRPSKI, LRFCZK, LRSDZK, LRFCBS, LRSDBS, LRGAUG
! -- end of lrescinf.h --
