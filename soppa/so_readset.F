C
C  /* Deck so_readset */
      SUBROUTINE SO_READSET(ARRAY,LVEC,NVEC,LUFILE,FILENM,LREC,JOFF)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, May 1997
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Read a set of arrays from a file.
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
         IERR = 0
         IOFF = IOFF + LREC
C
         CALL GETWA2(LUFILE,FILENM,ARRAY(1,IVEC),IOFF,LVEC)
C
      END DO
C
      DTIME     = SECOND()  - DTIME
      SOORWC(1) = SOORWC(1) + DTIME
C
      RETURN
C
      END
