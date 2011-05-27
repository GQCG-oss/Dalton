C
C  /* Deck so_collect_tm */
      SUBROUTINE SO_COLLECT_TM(ISYMTR,LABEL,IEXCI,TRMOM,TRLEN,TRVEL,
CSPAS:23/5-11: second and third moment sum rules
     &                         TQLEN,TQVEL,TTLEN,
CKeinSPASmehr
     &                         TRLON,TRMAG,BSRLON,EXENG)
C
C     This routine is part of the atomic integral direct SOPPA program.
C   
C     Keld Bak, August 1997
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Collect transition moments in the right arrays.
C
#include <implicit.h>
#include <priunit.h>
C
CSPAS:15/03-2006: merge with Dalton-2.0  
C#include <cbilrs.h>
#include <cbiexc.h>
CKeinSPASmehr
#include <ccorb.h>
#include <mxcent.h>
#include <nuclei.h>
C
      CHARACTER*8 LABEL
C
      DIMENSION TRLEN(3,NSYM,MXNEXI),        TRVEL(3,NSYM,MXNEXI)
CSPAS:23/5-11: second and third moment sum rules
      DIMENSION TQLEN(3,3,NSYM,MXNEXI), TQVEL(3,3,NSYM,MXNEXI)
      DIMENSION TTLEN(10,NSYM,MXNEXI)
CKeinSPASmehr
      DIMENSION TRLON(3,NSYM,MXNEXI),        TRMAG(3,NSYM,MXNEXI)
      DIMENSION BSRLON(3,NSYM,MXNEXI),       EXENG(NSYM,MXNEXI)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_COLLECT_TM')
C
C-----------------------------------------
C     Check that IEXCI is less that MXNEXI.
C-----------------------------------------
C
      IF ( IEXCI .GT. MXNEXI ) THEN
         WRITE(LUPRI,*) 'WARNING: IEXCI greater than MXNEXI'
         WRITE(LUPRI,*) 'Reset MXNEXI in the file cbiexc.h and resubmit'
         RETURN
      END IF
C 
C------------------------------------------------------------
C     The calculated transition moments are <f|Operator|0>. 
C     We want to keep them as <O|Operator|f>. Therefore, for 
C     anti-hermitian operators we multiply with a minus-sign.
C------------------------------------------------------------
C
      IF (LABEL(2:7).EQ.'DIPLEN') THEN
C
         IF (LABEL(1:1).EQ.'X') TRLEN(1,ISYMTR,IEXCI) =   TRMOM
         IF (LABEL(1:1).EQ.'Y') TRLEN(2,ISYMTR,IEXCI) =   TRMOM
         IF (LABEL(1:1).EQ.'Z') TRLEN(3,ISYMTR,IEXCI) =   TRMOM
C
      ELSE IF (LABEL(2:7).EQ.'DIPVEL') THEN
C
         IF (LABEL(1:1).EQ.'X') TRVEL(1,ISYMTR,IEXCI) = - TRMOM
         IF (LABEL(1:1).EQ.'Y') TRVEL(2,ISYMTR,IEXCI) = - TRMOM
         IF (LABEL(1:1).EQ.'Z') TRVEL(3,ISYMTR,IEXCI) = - TRMOM
C
CSPAS:23/5-11: second and third moment sum rules
      ELSE IF (LABEL(3:8).EQ.'SECMOM') THEN
         IF (LABEL(1:2) .EQ. 'XX') THEN
                                TQLEN(1,1,ISYMTR,IEXCI) = TRMOM
         ELSE IF (LABEL(1:2) .EQ. 'XY') THEN
                                TQLEN(1,2,ISYMTR,IEXCI) = TRMOM
                                TQLEN(2,1,ISYMTR,IEXCI) = TRMOM
         ELSE IF (LABEL(1:2) .EQ. 'XZ') THEN
                                TQLEN(1,3,ISYMTR,IEXCI) = TRMOM
                                TQLEN(3,1,ISYMTR,IEXCI) = TRMOM
         ELSE IF (LABEL(1:2) .EQ. 'YY') THEN
                                TQLEN(2,2,ISYMTR,IEXCI) = TRMOM
         ELSE IF (LABEL(1:2) .EQ. 'YZ') THEN
                                TQLEN(2,3,ISYMTR,IEXCI) = TRMOM
                                TQLEN(3,2,ISYMTR,IEXCI) = TRMOM
         ELSE IF (LABEL(1:2) .EQ. 'ZZ') THEN
                                TQLEN(3,3,ISYMTR,IEXCI) = TRMOM
         END IF
      ELSE IF (LABEL(5:8).EQ.'3MOM') THEN
         IF (LABEL(1:3) .EQ. 'XXX') THEN
                                TTLEN(1,ISYMTR,IEXCI) = TRMOM
         ELSE IF (LABEL(1:3) .EQ. 'YYY') THEN
                                TTLEN(2,ISYMTR,IEXCI) = TRMOM
         ELSE IF (LABEL(1:3) .EQ. 'ZZZ') THEN
                                TTLEN(3,ISYMTR,IEXCI) = TRMOM
         ELSE IF (LABEL(1:3) .EQ. 'XXY') THEN
                                TTLEN(4,ISYMTR,IEXCI) = TRMOM
         ELSE IF (LABEL(1:3) .EQ. 'XXZ') THEN
                                TTLEN(5,ISYMTR,IEXCI) = TRMOM
         ELSE IF (LABEL(1:3) .EQ. 'XYY') THEN
                                TTLEN(6,ISYMTR,IEXCI) = TRMOM
         ELSE IF (LABEL(1:3) .EQ. 'XYZ') THEN
                                TTLEN(7,ISYMTR,IEXCI) = TRMOM
         ELSE IF (LABEL(1:3) .EQ. 'XZZ') THEN
                                TTLEN(8,ISYMTR,IEXCI) = TRMOM
         ELSE IF (LABEL(1:3) .EQ. 'YYZ') THEN
                                TTLEN(9,ISYMTR,IEXCI) = TRMOM
         ELSE IF (LABEL(1:3) .EQ. 'YZZ') THEN
                                TTLEN(10,ISYMTR,IEXCI) = TRMOM
         END IF
CKeinSPASmehr
C
      ELSE IF (LABEL(2:7).EQ.'ANGMOM') THEN
C
         IF (LABEL(1:1).EQ.'X') TRMAG(1,ISYMTR,IEXCI) = - TRMOM
         IF (LABEL(1:1).EQ.'Y') TRMAG(2,ISYMTR,IEXCI) = - TRMOM
         IF (LABEL(1:1).EQ.'Z') TRMAG(3,ISYMTR,IEXCI) = - TRMOM
C
      ELSE IF (LABEL(2:7).EQ.'LONMAG') THEN
C
         IF (LABEL(1:1).EQ.'X') TRLON(1,ISYMTR,IEXCI) = - TRMOM
         IF (LABEL(1:1).EQ.'Y') TRLON(2,ISYMTR,IEXCI) = - TRMOM
         IF (LABEL(1:1).EQ.'Z') TRLON(3,ISYMTR,IEXCI) = - TRMOM
C
      ELSE IF (LABEL(2:6).EQ.'HBDO ') THEN
C
         IF (LABEL(7:7).EQ.'X') BSRLON(1,ISYMTR,IEXCI) =
     &                          -TRMOM * EXENG(ISYMTR,IEXCI)
         IF (LABEL(7:7).EQ.'Y') BSRLON(2,ISYMTR,IEXCI) =
     &                          -TRMOM * EXENG(ISYMTR,IEXCI)
         IF (LABEL(7:7).EQ.'Z') BSRLON(3,ISYMTR,IEXCI) =
     &                          -TRMOM * EXENG(ISYMTR,IEXCI)
C
      ENDIF
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_COLLECT_TM')
C
      RETURN
      END
