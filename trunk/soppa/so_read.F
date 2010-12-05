C
C  /* Deck so_read */
      SUBROUTINE SO_READ(ARRAY,LARRAY,LUFILE,FILENM,IREC)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, May 1996
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Read array from LUFILE. Notice that it is assumed
C              that the file is open and positioned correctly.
C
#include <implicit.h>
#include <soppinf.h>
C
      DIMENSION   ARRAY(LARRAY)
      CHARACTER*7 FILENM
C
      DTIME      = SECOND()
C
      IF (LARRAY .EQ. 0) RETURN
C
      IOFF = LARRAY * (IREC - 1) + 1
C
      CALL GETWA2(LUFILE,FILENM,ARRAY,IOFF,LARRAY)
C
      DTIME     = SECOND()  - DTIME
      SOORWC(1) = SOORWC(1) + DTIME
C
      RETURN
C
      END
