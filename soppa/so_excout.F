C
C  /* Deck so_excout */
CSPAS:23/5-11: second and third moment sum rules
C     SUBROUTINE SO_EXCOUT(TRLEN,TRVEL,TRMAG,TRLON,BSRLON,EXENG,
C    &                     FONAC,FONA2,SECMAT,WORK,LWORK)
      SUBROUTINE SO_EXCOUT(TRLEN,TRVEL,TQLEN,TQVEL,TRMAG,TRLON,BSRLON,
     &                     TTLEN,EXENG,FONAC,FONA2,SECMAT,WORK,LWORK)
CKeinSPASmehr
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, November 1997
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Driver for printing final output for the RPA, RPA(D), 
C              SOPPA, and SOPPA(CCSD) excitation calculations.
C
#include <implicit.h>
#include <priunit.h>
C
#include <ccsdsym.h>
#include <ccorb.h>
#include <cbiexc.h>
CSPAS:15/03-2006: merge with Dalton-2.0  
C#include <cbilrs.h>
CKeinSPASmehr
#include <soppinf.h>
#include <maxorb.h>
#include <secmom.h>
#include <codata.h>
C
      DIMENSION TRLEN(3,NSYM,MXNEXI,4), TRVEL(3,NSYM,MXNEXI,4)
CSPAS:23/5-11: second and third moment sum rules
      DIMENSION TQLEN(3,3,NSYM,MXNEXI,4), TQVEL(3,3,NSYM,MXNEXI,4)
CKeinSPASmehr
      DIMENSION TRLON(3,NSYM,MXNEXI,4), TRMAG(3,NSYM,MXNEXI,4)
CSPAS:23/5-11: second and third moment sum rules
      DIMENSION TTLEN(10,NSYM,MXNEXI,4)
CKeinSPASmehr
      DIMENSION BSRLON(3,NSYM,MXNEXI,4),EXENG(NSYM,MXNEXI,4)
      DIMENSION FONAC(*)
      DIMENSION FONA2(*)
      DIMENSION SECMAT(3,MXNEXI,NSYM)
      DIMENSION WORK(LWORK)
      DIMENSION NEXCTX(8)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_EXCOUT')
C
C===============================================
C     Write calculated RPA properties to output.
C===============================================
C
      IF (AORPA) THEN
C
         DO 100 ISYM = 1, NSYM
C
            NEXCTX(ISYM) = NEXCIT(ISYM)
            NEXCIT(ISYM) = MIN(NEXCIT(ISYM),NT1AM(ISYM))
C
  100    CONTINUE
C
         WRITE(LUPRI,9000)
         WRITE(LUPRI,'(30X,A)') ' RPA results:'
         WRITE(LUPRI,9001)
C
CSPAS:23/5-11: second and third moment sum rules
C        CALL EXCOUT(TRLEN(1,1,1,1),TRVEL(1,1,1,1),DUMMY,DUMMY,
C    &               TRMAG(1,1,1,1),TRLON(1,1,1,1),BSRLON(1,1,1,1),
C    &               EXENG(1,1,1),FONAC,FONA2,DUMMY,DUMMY,DUMMY,DUMMY,
         CALL EXCOUT(TRLEN(1,1,1,1),TRVEL(1,1,1,1),TQLEN(1,1,1,1,1),
     &               TQVEL(1,1,1,1,1),TRMAG(1,1,1,1),TRLON(1,1,1,1),
     &               BSRLON(1,1,1,1),TTLEN(1,1,1,1),EXENG(1,1,1),
     &               FONAC,FONA2,DUMMY,DUMMY,DUMMY,DUMMY,
CKeinSPASmehr
     &               DUMMY,DUMMY,DUMMY,DUMMY,WORK,LWORK)
C
         DO 200 ISYM = 1, NSYM
C
            NEXCIT(ISYM) = NEXCTX(ISYM)
C
  200    CONTINUE
C
      END IF
C
C--------------------------------------------------
C     Write <R**2> of molecular orbitals to output.
C--------------------------------------------------
C
      WRITE(LUPRI,'(//,1X,A)')
     &   '         Extent of molecular orbitals (au)'
      WRITE(LUPRI,'(1X,A/)')
     &   '         ---------------------------------'
      WRITE(LUPRI,'(1X,A)')
     &   ' Sym.  MO.   <X**2>    <Y**2>    <Z**2>    <R**2>'
      WRITE(LUPRI,'(1X,A)')
     &   ' ------------------------------------------------'
C
      DO 220 ISYM = 1, NSYM 
C
         DO 210 IMO = 1, NRHF(ISYM)
C
            SECSUM = SECOMO(1,IMO,ISYM) + SECOMO(2,IMO,ISYM) 
     &             + SECOMO(3,IMO,ISYM)
            WRITE(LUPRI,'(2I5,4F10.3)') ISYM, IMO,
     &         (SECOMO(I,IMO,ISYM), I = 1,3), SECSUM
C
  210    CONTINUE
C
  220 CONTINUE
C
      WRITE(LUPRI,'(1X)')
C
      DO 221 ISYMA = 1, NSYM
C
         DO 211 IMOA = 1, NVIR(ISYMA)
C
            SECSUM = SECVMO(1,IMOA,ISYMA) + SECVMO(2,IMOA,ISYMA)
     &             + SECVMO(3,IMOA,ISYMA)
            WRITE(LUPRI,'(2I5,4F10.3)') ISYMA, IMOA,
     &         (SECVMO(I,IMOA,ISYMA), I = 1,3), SECSUM
C
  211    CONTINUE
C
  221 CONTINUE
C--------------------------------------------
C     Write <R**2> of ground state to output.
C--------------------------------------------
C
      WRITE(LUPRI,'(//,35X,A)') 'Ground state (au)'
      WRITE(LUPRI,'(25X,A)') 
     &   '<X**2>    <Y**2>    <Z**2>    <R**2>'
      WRITE(LUPRI,'(25X,A)') 
     &   '------------------------------------'
      SECSUM = SECGR(1) + SECGR(2) + SECGR(3)
      WRITE(LUPRI,'(21X,4F10.3)') (SECGR(I), I = 1,3), SECSUM
C
C-------------------------------------------------------------------
C     Write change of <R**2> from ground to excited state to output.
C-------------------------------------------------------------------
C
      WRITE(LUPRI,'(//,7X,A)') 
     &   'Change in <R**2> from ground to excited state (au)'
      WRITE(LUPRI,'(7X,A/)') 
     &   '--------------------------------------------------'
      WRITE(LUPRI,'(1X,A)') 
     &   ' Sym. State Freq. (eV)  <X**2>    <Y**2>    <Z**2>    <R**2>'
      WRITE(LUPRI,'(1X,A)') 
     &   ' -----------------------------------------------------------'
C
      DO 240 ISYM = 1, NSYM
C
CSPAS:24/5-11: remove a lot of garbage output
C        DO 230 IEXCI = 1, NEXCIT(ISYM)
         DO 230 IEXCI = 1, MIN(NEXCIT(ISYM),NT1AM(ISYM))
CKeinSPASmehr
C
            SECSUM = SECMAT(1,IEXCI,ISYM) 
     &             + SECMAT(2,IEXCI,ISYM) + SECMAT(3,IEXCI,ISYM)
C
            WRITE(LUPRI,'(2I5,F11.4,4F10.3)') ISYM, IEXCI, 
     &         XTEV * EXENG(ISYM,IEXCI,1),
     &         ( SECMAT(I,IEXCI,ISYM), I = 1,3), SECSUM
C
  230    CONTINUE
C
  240 CONTINUE
C
C==================================================
C     Write calculated RPA(D) properties to output.
C==================================================
C
      IF (DCRPA) THEN
C
         DO 300 ISYM = 1, NSYM
C
            NEXCTX(ISYM) = NEXCIT(ISYM)
            NEXCIT(ISYM) = MIN(NEXCIT(ISYM),(NT1AM(ISYM)+NT2AM(ISYM)))
C
  300    CONTINUE
C
         WRITE(LUPRI,9000)
         WRITE(LUPRI,'(30X,A)') 'RPA(D) results:'
         WRITE(LUPRI,9001)
C
CSPAS:23/5-11: second and third moment sum rules
C        CALL EXCOUT(TRLEN(1,1,1,2),TRVEL(1,1,1,2),DUMMY,DUMMY,
C    &               TRMAG(1,1,1,2),TRLON(1,1,1,2),BSRLON(1,1,1,2),
C    &               EXENG(1,1,2),FONAC,FONA2,DUMMY,DUMMY,DUMMY,DUMMY,
         CALL EXCOUT(TRLEN(1,1,1,2),TRVEL(1,1,1,2),TQLEN(1,1,1,1,2),
     &               TQVEL(1,1,1,1,2),TRMAG(1,1,1,2),TRLON(1,1,1,2),
     &               BSRLON(1,1,1,2),TTLEN(1,1,1,2),EXENG(1,1,2),
     &               FONAC,FONA2,DUMMY,DUMMY,DUMMY,DUMMY,
CKeinSPASmehr
     &               DUMMY,DUMMY,DUMMY,DUMMY,WORK,LWORK)
C
         DO 400 ISYM = 1, NSYM
C
            NEXCIT(ISYM) = NEXCTX(ISYM)
C
  400    CONTINUE
C
      END IF
C
C=================================================
C     Write calculated SOPPA properties to output.
C=================================================
C
      IF (AOSOP) THEN
C
         DO 500 ISYM = 1, NSYM
C
            NEXCTX(ISYM) = NEXCIT(ISYM)
            NEXCIT(ISYM) = MIN(NEXCIT(ISYM),(NT1AM(ISYM)+NT2AM(ISYM)))
C
  500    CONTINUE
C
         WRITE(LUPRI,9000)
         WRITE(LUPRI,'(30X,A)') 'SOPPA results:'
         WRITE(LUPRI,9001)
C
CSPAS:23/5-11: second and third moment sum rules
C        CALL EXCOUT(TRLEN(1,1,1,3),TRVEL(1,1,1,3),DUMMY,DUMMY,
C    &               TRMAG(1,1,1,3),TRLON(1,1,1,3),BSRLON(1,1,1,3),
C    &               EXENG(1,1,3),FONAC,FONA2,DUMMY,DUMMY,DUMMY,DUMMY,
         CALL EXCOUT(TRLEN(1,1,1,3),TRVEL(1,1,1,3),TQLEN(1,1,1,1,3),
     &               TQVEL(1,1,1,1,3),TRMAG(1,1,1,3),TRLON(1,1,1,3),
     &               BSRLON(1,1,1,3),TTLEN(1,1,1,3),EXENG(1,1,3),
     &               FONAC,FONA2,DUMMY,DUMMY,DUMMY,DUMMY,
CKeinSPASmehr
     &               DUMMY,DUMMY,DUMMY,DUMMY,WORK,LWORK)
C
         DO 600 ISYM = 1, NSYM
C
            NEXCIT(ISYM) = NEXCTX(ISYM)
C
  600    CONTINUE
C
      END IF
C
C=======================================================
C     Write calculated SOPPA(CCSD) properties to output.
C=======================================================
C
      IF (AOSOC) THEN
C
         DO 700 ISYM = 1, NSYM
C
            NEXCTX(ISYM) = NEXCIT(ISYM)
            NEXCIT(ISYM) = MIN(NEXCIT(ISYM),(NT1AM(ISYM)+NT2AM(ISYM)))
C
  700    CONTINUE
C
         WRITE(LUPRI,9000)
         WRITE(LUPRI,'(30X,A)') 'SOPPA(CCSD) results:'
         WRITE(LUPRI,9001)
C
CSPAS:23/5-11: second and third moment sum rules
C        CALL EXCOUT(TRLEN(1,1,1,4),TRVEL(1,1,1,4),DUMMY,DUMMY,
C    &               TRMAG(1,1,1,4),TRLON(1,1,1,4),BSRLON(1,1,1,4),
C    &               EXENG(1,1,4),FONAC,FONA2,DUMMY,DUMMY,DUMMY,DUMMY,
         CALL EXCOUT(TRLEN(1,1,1,4),TRVEL(1,1,1,4),TQLEN(1,1,1,1,4),
     &               TQVEL(1,1,1,1,4),TRMAG(1,1,1,4),TRLON(1,1,1,4),
     &               BSRLON(1,1,1,4),TTLEN(1,1,1,4),EXENG(1,1,4),
     &               FONAC,FONA2,DUMMY,DUMMY,DUMMY,DUMMY,
CKeinSPASmehr
     &               DUMMY,DUMMY,DUMMY,DUMMY,WORK,LWORK)
C
         DO 800 ISYM = 1, NSYM
C
            NEXCIT(ISYM) = NEXCTX(ISYM)
C
  800    CONTINUE
C
      END IF
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_EXCOUT')
C
      RETURN
C
 9000 FORMAT(//' =========================================',
     &       '=============================')
 9001 FORMAT(' -----------------------------------------',
     &       '-----------------------------')
      END
