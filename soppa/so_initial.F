C  /* Deck so_initial */
      SUBROUTINE SO_INITIAL(IREC1,IREC2,ARRAY,LARRAY,LUFILE,FILENM)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, May 1996
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Initialize the NR arrays ARRAY of length LARRAY
C              kept in file LUFILE to zero.
C
#include <implicit.h>
#include <soppinf.h>
C
      DIMENSION   ARRAY(LARRAY)
      CHARACTER*7 FILENM
C
      CALL DZERO(ARRAY,LARRAY)
C
      CALL SO_OPEN(LUFILE,FILENM,LARRAY)
C
      DO 100 I = IREC1,IREC2
C
         CALL SO_WRITE(ARRAY,LARRAY,LUFILE,FILENM,I)
C
  100 CONTINUE
C
      CALL SO_CLOSE(LUFILE,FILENM,'KEEP')
C
      RETURN
      END
