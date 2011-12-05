C
C  /* Deck so_write */
      SUBROUTINE SO_WRITE(ARRAY,LARRAY,LUFILE,FILENM,IREC)
C
C     Keld Bak, May 1996
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Write array from LUFILE. Notice that it is assumed
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
      CALL PUTWA2(LUFILE,FILENM,ARRAY,IOFF,LARRAY)
C
      DTIME     = SECOND()  - DTIME
      SOORWC(2) = SOORWC(2) + DTIME
C
      RETURN
C
      END
