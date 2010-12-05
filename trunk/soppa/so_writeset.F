C
C  /* Deck so_writeset */
      SUBROUTINE SO_WRITESET(ARRAY,LVEC,NVEC,LUFILE,FILENM,LREC,JOFF)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, May 1997
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Write a set of arrays to a file.
C
#include <implicit.h>
#include <soppinf.h>
C
      DIMENSION   ARRAY(LVEC,NVEC)
      CHARACTER*7 FILENM
C
      DTIME      = SECOND()
C
      IF (LVEC .EQ. 0) RETURN
C
      IOFF = JOFF - LREC
C
      DO IVEC = 1,NVEC
C
         IOFF = IOFF + LREC
C
         CALL PUTWA2(LUFILE,FILENM,ARRAY(1,IVEC),IOFF,LVEC)
C
      END DO
C
      DTIME     = SECOND()  - DTIME
      SOORWC(2) = SOORWC(2) + DTIME
C
      RETURN
C
      END
