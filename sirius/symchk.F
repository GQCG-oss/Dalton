!
!...   Copyright (c) 2011 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2011 (2011), see http://daltonprogram.org"
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dalton,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For further information, including how to get a licence, see:
!...      http://daltonprogram.org
!
!
#include <single.h>
C
C  /* Deck cichck */
      SUBROUTINE CICHCK (WRK,LFREE,CCHECK)
C 890302-hjaaj -- interface to ci check routine
#include <implicit.h>
      DIMENSION WRK(LFREE)
      INTEGER   CCHECK
C
C Used from common blocks:
C   INFINP: ISTATE,LROOTS,NROOTS
C
#include <maxorb.h>
#include <infinp.h>
      KTMAT = 1
      KVEC  = KTMAT + LROOTS*LROOTS
      LVEC  = LFREE - KVEC
      CALL CICHC2 (WRK(KTMAT),WRK(KVEC),LVEC,CCHECK)
      RETURN
      END
      SUBROUTINE CICHC2 (TMAT,VEC,LVEC,CCHECK)
C
C  7-AUG. 1986 Hans Agren
C  Revisions:
C
C Purpose:
C        Called after SIROPT.OPTST.CICTL to check
C        symmetry of LROOTS start CI-vectors, and to remove those
C        which have wrong symmetry (e.g. delta vectors in a sigma
C        symmetry calculation). The check is performed on the matrix
C        element <VEC1 | oper | VEC2>, where "oper" is the CI-diagonal.
C        There is two options for this check:
C
C        CCHECK = 1 : Remove those vectors which do not have the same
C                     symmetry as the ISTATE vector, reassign ISTATE
C        CCHECK = 2 : Remove those vectors which do not have the same
C                     symmetry as the lowest state vector
C        other values: check symmetry, do not remove any.
C
C        In both cases NROOTS and LROOTS are updated.
C
C Input:
C        LVEC  : length of VEC()
C
C Scratch:
C        TMAT  : matrix containing operator product of VEC1 and VEC2
C               (operator = CI diagonal)
C        VEC   : used for CI vectors and CI diagonal from unit LUIT2
C
C
#include <implicit.h>
#include <dummy.h>
#include <infvar.h>
      DIMENSION VEC(NCONF,*),TMAT(LROOTS,*)
      INTEGER   CCHECK
C
      PARAMETER ( D0=0.0D0, D1=1.0D0 )
      PARAMETER ( THRSYM = 1.D-7)
C
C Used from common blocks:
C   INFINP: ISTATE,LROOTS,NROOTS
C   INFVAR: NCONF
C   INFTAP: LUIT1, LUTEMP
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inftap.h>
#include <infpri.h>
C
      LOGICAL WRTEMP
      CHARACTER*8 TABLE(3),LABIT1(2)
C
      DATA TABLE/'********','CIDIAG2 ','STARTVEC'/
C
      CALL QENTER('CICHCK')
      NCNF4 = MAX(4,NCONF)
      IF (LROOTS .LE. 1) GO TO 9999
      JCIDIA = LVEC / NCONF
      MCVEC  = MIN(LROOTS,JCIDIA - 1)
      IF (MCVEC .LT. 2) CALL ERRWRK('CICHCK',3*NCONF,LVEC)
C
C     Read CI - diagonal from LUIT2
C
      REWIND LUIT2
      CALL MOLLAB(TABLE(2),LUIT2,lupri)
      CALL READT (LUIT2,NCONF,VEC(1,JCIDIA))
C
C     Read CI - vectors from LUIT1
C
      LUTEMP = -1
      IF ((JCIDIA-1) .LT. LROOTS) THEN
         CALL GPOPEN(LUTEMP,' ','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         WRTEMP = .TRUE.
         REWIND LUTEMP
      ELSE
         WRTEMP = .FALSE.
      END IF
      REWIND LUIT1
      CALL MOLLAB(TABLE(3),LUIT1,lupri)
      DO 50 I = 1,MCVEC
         CALL READT (LUIT1,NCONF,VEC(1,I))
         IF (WRTEMP) CALL WRITT (LUTEMP,NCONF,VEC(1,I))
   50 CONTINUE
C
      DO 100 I = 1,MCVEC
         DO 200 J = 1,I
            TMAT(I,J) = DV3DOT(NCONF,VEC(1,I),VEC(1,JCIDIA),VEC(1,J))
            TMAT(J,I) = TMAT(I,J)
  200    CONTINUE
  100 CONTINUE
C
  400 CONTINUE
      IF (MCVEC .LT. LROOTS) THEN
         NCVEC = MIN(JCIDIA - 2, LROOTS - MCVEC)
         ICVEC = JCIDIA - 1
         DO 450 I = 1,NCVEC
            CALL READT (LUIT1,NCONF,VEC(1,I))
            IF (WRTEMP) CALL WRITT (LUTEMP,NCONF,VEC(1,I))
  450    CONTINUE
         REWIND LUTEMP
         DO 560 I = 1,MCVEC
            CALL READT (LUTEMP,NCONF,VEC(1,ICVEC))
            DO 540 J = 1,NCVEC
               TMAT(I,MCVEC+J) =
     *         DV3DOT(NCONF,VEC(1,ICVEC),VEC(1,JCIDIA),VEC(1,J))
               TMAT(MCVEC+J,I) = TMAT(I,MCVEC+J)
  540       CONTINUE
  560    CONTINUE
         DO 640 I = 1,NCVEC
            READ (LUTEMP)
            DO 620 J = 1,I
               TMAT(MCVEC+I,MCVEC+J) =
     *         DV3DOT(NCONF,VEC(1,I),VEC(1,JCIDIA),VEC(1,J))
               TMAT(MCVEC+J,MCVEC+I) = TMAT(MCVEC+I,MCVEC+J)
  620       CONTINUE
  640    CONTINUE
         MCVEC = MCVEC + NCVEC
         GO TO 400
      END IF
C
C  Check vectors and write the correct ones out on LUIT1
C
      REWIND LUIT1
      CALL MOLLAB(TABLE(3),LUIT1,lupri)
      LROOTX = 0
      IREMOV = 0
      JSTATE = ISTATE
      IF (CCHECK .EQ. 1) THEN
C
C        Option 1. remove vectors with symmetry different from ISTATE:s
C
         DO 700 J = 1,LROOTS
            IF ( ABS( TMAT(ISTATE,J)) .GT.  THRSYM) THEN
               LROOTX = LROOTX + 1
               IF (J .EQ. ISTATE) JSTATE = LROOTX
            ELSE
               IREMOV = IREMOV + 1
               WRITE (LUPRI,'(/A,I3)')
     *         ' >>> CICHCK: removed CI start vector no.',J
            END IF
  700    CONTINUE
         IF (LROOTX .NE. LROOTS) THEN
            BACKSPACE LUIT1
            WRITE (LABIT1(1),'(2I4)') LROOTX,JSTATE
            WRITE (LABIT1(2),'(A8)') '(CICHCK)'
            WRITE (LUIT1) TABLE(1),LABIT1,TABLE(3)
            DO 710 J = 1,LROOTS
               IF ( ABS( TMAT(ISTATE,J)) .GT.  THRSYM) THEN
                  IF (WRTEMP) THEN
                     CALL READT(LUTEMP,NCONF,VEC)
                     JJ = 1
                  ELSE
                     JJ = J
                  END IF
                  CALL WRITT (LUIT1,NCNF4,VEC(1,JJ))
               ELSE
                  IF (WRTEMP) READ (LUTEMP)
               END IF
  710       CONTINUE
         END IF
C
      ELSE IF (CCHECK .EQ. 2) THEN
C
C        Option 2. remove vectors with symmetry different from
C                  lowest state
C
         DO 800 J = 1,LROOTS
            IF ( ABS( TMAT(1,J)) .GT.  THRSYM) THEN
               LROOTX = LROOTX + 1
            ELSE
               IREMOV = IREMOV + 1
               WRITE (LUPRI,'(/A,I3)')
     *         ' >>> CICHCK: removed CI start vector no.',J
            END IF
  800    CONTINUE
         JSTATE = ISTATE
         IF (LROOTX .NE. LROOTS) THEN
            BACKSPACE LUIT1
            WRITE (LABIT1(1),'(2I4)') LROOTX,JSTATE
            WRITE (LABIT1(2),'(A8)') '(CICHCK)'
            WRITE (LUIT1) TABLE(1),LABIT1,TABLE(3)
            DO 810 J = 1,LROOTS
               IF ( ABS( TMAT(1,J)) .GT.  THRSYM) THEN
                  IF (WRTEMP) THEN
                     CALL READT(LUTEMP,NCONF,VEC)
                     JJ = 1
                  ELSE
                     JJ = J
                  END IF
                  CALL WRITT (LUIT1,NCNF4,VEC(1,JJ))
               ELSE
                  IF (WRTEMP) READ (LUTEMP)
               END IF
  810       CONTINUE
         END IF
      ELSE
         LROOTX = LROOTS
      END IF
      NROOTX = MIN(NROOTS+JSTATE-ISTATE,LROOTX)
C
      IF (WRTEMP) CALL GPCLOSE(LUTEMP,'DELETE')
C
      IF (CCHECK .GE. 1 .AND. CCHECK .LE. 2) THEN
         IF (IREMOV .EQ. 0 .AND. IPRI6 .GE. 10) IREMOV = -1
      ELSE
C        do not remove CI vectors, just make a check
         IREMOV = -1
      END IF
C
      IF (IREMOV .NE. 0) THEN
         IREMOV = MAX(0,IREMOV)
         WRITE (LUPRI,'(///A/A,I3,A,I3,3(/A,2I4)//A)')
     *   ' >>> CICHCK; REMOVE CI VECTORS OF NON-DESIRED SYMMETRY',
     *   '             OPTION =',CCHECK,
     *   ', number of CI trial vectors removed :',IREMOV,
     *   '             ISTATE, orignal and after reduce:',ISTATE,JSTATE,
     *   '             LROOTS, orignal and after reduce:',LROOTS,LROOTX,
     *   '             NROOTS, orignal and after reduce:',NROOTS,NROOTX,
     *   ' >>> CICHCK; THE SYMMETRY TEST MATRIX:'
         CALL OUTPUT(TMAT,1,LROOTS,1,LROOTS,LROOTS,LROOTS,1,LUPRI)
         IF (NROOTX .LT. JSTATE) THEN
            WRITE (LUPRI,'(//A/A,I5/A,I3)')
     *      ' >>> CICHCK option 2, cannot proceed:',
     *      '     number of trial vectors of right symmetry is',NROOTX,
     *      '     which is less than current reference state no.',JSTATE
            CALL QUIT(' >>> CICHCK: too few trial CI '//
     &                'vecs. of right symmetry')
         END IF
      END IF
C
      IF (CCHECK .GE. 1 .AND. CCHECK .LE. 2) THEN
         ISTATE = JSTATE
         LROOTS = LROOTX
         NROOTS = NROOTX
         IROOT(1) = ISTATE
         I = 1
         DO 900 J = 1,LROOTS
            IF (J .NE. ISTATE) THEN
               I = I + 1
               IROOT(I) = J
            END IF
  900    CONTINUE
      END IF
C
 9999 CALL QEXIT('CICHCK')
      RETURN
C
C     END OF CICHCK.
C
      END
