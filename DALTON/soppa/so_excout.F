C
C  /* Deck so_excout */
      SUBROUTINE SO_EXCOUT(TRLEN,TRVEL,TQLEN,TQVEL,TRMAG,TRLON,BSRLON,
     &                     TTLEN,EXENG,
CClark:11/01/2016
     &                     BETHE,STOPP,
CClark:end
     &                     FONAC,FONA2,SECMAT,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, November 1997
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C     PFP & SPAS, November 2013: HRPA and triplet excitation energies
C
C     PURPOSE: Driver for printing final output for the RPA, HRPA, RPA(D), 
C              SOPPA, and SOPPA(CCSD) excitation calculations.
C
#include "implicit.h"
#include "priunit.h"
C
#include "ccsdsym.h"
#include "ccorb.h"
#include "cbiexc.h"
#include "soppinf.h"
#include "maxorb.h"
#include "secmom.h"
#include "codata.h"
C
CSPAS: HRPA i.e. 5 models now
C     DIMENSION TRLEN(3,NSYM,MXNEXI,4), TRVEL(3,NSYM,MXNEXI,4)
C     DIMENSION TQLEN(3,3,NSYM,MXNEXI,4), TQVEL(3,3,NSYM,MXNEXI,4)
C     DIMENSION TRLON(3,NSYM,MXNEXI,4), TRMAG(3,NSYM,MXNEXI,4)
C     DIMENSION TTLEN(10,NSYM,MXNEXI,4)
C     DIMENSION BSRLON(3,NSYM,MXNEXI,4),EXENG(NSYM,MXNEXI,4)
      DIMENSION TRLEN(3,NSYM,MXNEXI,5), TRVEL(3,NSYM,MXNEXI,5)
      DIMENSION TQLEN(3,3,NSYM,MXNEXI,5), TQVEL(3,3,NSYM,MXNEXI,5)
      DIMENSION TRLON(3,NSYM,MXNEXI,5), TRMAG(3,NSYM,MXNEXI,5)
      DIMENSION TTLEN(10,NSYM,MXNEXI,5)
      DIMENSION BSRLON(3,NSYM,MXNEXI,5),EXENG(NSYM,MXNEXI,5)
CClark:11/01/2016
      DIMENSION BETHE(3,LQ,5),STOPP(3,LVEL,5)
CClark:end
CKeinSPASmehr
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
         DO ISYM = 1, NSYM
C
            NEXCTX(ISYM) = NEXCIT(ISYM)
            NEXCIT(ISYM) = MIN(NEXCIT(ISYM),NT1AM(ISYM))
C
         END DO  
C
         WRITE(LUPRI,9000)
         WRITE(LUPRI,'(30X,A)') ' RPA results:'
         WRITE(LUPRI,9001)
C
         CALL EXCOUT(TRLEN(1,1,1,1),TRVEL(1,1,1,1),TQLEN(1,1,1,1,1),
     &               TQVEL(1,1,1,1,1),TRMAG(1,1,1,1),TRLON(1,1,1,1),
     &               BSRLON(1,1,1,1),TTLEN(1,1,1,1),EXENG(1,1,1),
     &               FONAC,FONA2,DUMMY,DUMMY,DUMMY,DUMMY,
     &               DUMMY,DUMMY,DUMMY,DUMMY,WORK,LWORK)
C
CSPAS:26/01-08: 
C     this gives a lot of non-sense output for the orbital extend
C
C        DO 200 ISYM = 1, NSYM
C
C           NEXCIT(ISYM) = NEXCTX(ISYM)
C
C 200    END DO
C
C     END IF
CKeinSPASmehr
C
CClark:11/01/2016
C
C---------------------------------------------
C     Write Bethe sum rule and Stopping power.
C---------------------------------------------

      WRITE(LUPRI,'(//,1X,A)')
     &   '         Bethe sum rule                   '
      WRITE(LUPRI,'(1X,A/)')
     &   '         --------------'
      WRITE(LUPRI,'(1X,A)')
     &   ' Q value      X Component    Y Component    Z Component'
      WRITE(LUPRI,'(1X,A)')
     &   ' ------------------------------------------------------'
      WRITE(LUPRI,'(1X,A,I4,3(1X,A,G10.4))')
     &   'LQ=',LQ,'QMAX=',QMAX,'QMIN=',QMIN,'QSTEP=',QSTEP
      DO IQ = 1, LQ
C
         WRITE(LUPRI,'(2X,G9.4,3(4X,G11.5))')
     &                QMIN+(IQ-1)*QSTEP,(BETHE(K,IQ,1),K=1,3)
C
      END DO
C
      WRITE(LUPRI,'(//,1X,A)')
     &   '         Stopping Power                 '
      WRITE(LUPRI,'(1X,A/)')
     &   '         --------------'
      WRITE(LUPRI,'(1X,A)')
     &   ' Velocity        Component       Value           '
      WRITE(LUPRI,'(1X,A)')
     &   ' ------------------------------------------------'
C
      DO IVEL = 1, LVEL
C
         VELOC = VMIN+(IVEL-1)*VSTEP
C
         WRITE(LUPRI,'(2X,G9.4,7X,A,10X,G17.6)')
     &        VELOC,'X',STOPP(1,IVEL,1)
         WRITE(LUPRI,'(2X,G9.4,7X,A,10X,G17.6)')
     &        VELOC,'Y',STOPP(2,IVEL,1)
         WRITE(LUPRI,'(2X,G9.4,7X,A,10X,G17.6)')
     &        VELOC,'Z',STOPP(3,IVEL,1)
C
      END DO
C
CClark:end
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
      DO ISYM = 1, NSYM 
C
         DO IMO = 1, NRHF(ISYM)
C
            SECSUM = SECOMO(1,IMO,ISYM) + SECOMO(2,IMO,ISYM) 
     &             + SECOMO(3,IMO,ISYM)
            WRITE(LUPRI,'(2I5,4F10.3)') ISYM, IMO,
     &         (SECOMO(I,IMO,ISYM), I = 1,3), SECSUM
C
         END DO  
C
      END DO  
C
      WRITE(LUPRI,'(1X)')
C
      DO ISYMA = 1, NSYM
C
         DO IMOA = 1, NVIR(ISYMA)
C
            SECSUM = SECVMO(1,IMOA,ISYMA) + SECVMO(2,IMOA,ISYMA)
     &             + SECVMO(3,IMOA,ISYMA)
            WRITE(LUPRI,'(2I5,4F10.3)') ISYMA, IMOA,
     &         (SECVMO(I,IMOA,ISYMA), I = 1,3), SECSUM
C
         END DO  
C
      END DO  
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
      DO ISYM = 1, NSYM
C
         DO IEXCI = 1, NEXCIT(ISYM)
C
            SECSUM = SECMAT(1,IEXCI,ISYM) 
     &             + SECMAT(2,IEXCI,ISYM) + SECMAT(3,IEXCI,ISYM)
C
            WRITE(LUPRI,'(2I5,F11.4,4F10.3)') ISYM, IEXCI, 
     &         XTEV * EXENG(ISYM,IEXCI,1),
     &         ( SECMAT(I,IEXCI,ISYM), I = 1,3), SECSUM
C
         END DO
C
      END DO
C
CSPAS:26/01-08: 
C     moved because otherwise we get a lot of non-sense output 
C     for the orbital extend
         DO ISYM = 1, NSYM
C
            NEXCIT(ISYM) = NEXCTX(ISYM)
C
         END DO
C
      END IF
CKeinSPASmehr
CPFP HRPA included
C
C================================================
C     Write calculated HRPA properties to output.
C================================================
C
      IF (AOHRP) THEN
C
         DO ISYM = 1, NSYM
C
            NEXCTX(ISYM) = NEXCIT(ISYM)
            NEXCIT(ISYM) = MIN(NEXCIT(ISYM),NT1AM(ISYM))
C
         END DO
C
         WRITE(LUPRI,9000)
         WRITE(LUPRI,'(30X,A)') 'HRPA results:'
         WRITE(LUPRI,9001)
C
         CALL EXCOUT(TRLEN(1,1,1,2),TRVEL(1,1,1,2),TQLEN(1,1,1,1,2),
     &               TQVEL(1,1,1,1,2),TRMAG(1,1,1,2),TRLON(1,1,1,2),
     &               BSRLON(1,1,1,2),TTLEN(1,1,1,2),EXENG(1,1,2),
     &               FONAC,FONA2,DUMMY,DUMMY,DUMMY,DUMMY,
     &               DUMMY,DUMMY,DUMMY,DUMMY,WORK,LWORK)
C
         DO ISYM = 1, NSYM
C
            NEXCIT(ISYM) = NEXCTX(ISYM)
C
         END DO
C
      END IF
C
Cend-PFP
C
C==================================================
C     Write calculated RPA(D) properties to output.
C==================================================
C
      IF (DCRPA) THEN
C
         DO ISYM = 1, NSYM
C
            NEXCTX(ISYM) = NEXCIT(ISYM)
CSPAS:26/01-08: in RPA(D) we have only single excitations !!
C           NEXCIT(ISYM) = MIN(NEXCIT(ISYM),(NT1AM(ISYM)+NT2AM(ISYM)))
            NEXCIT(ISYM) = MIN(NEXCIT(ISYM),NT1AM(ISYM))
CKeinSPASmehr
C
         END DO
C
         WRITE(LUPRI,9000)
         WRITE(LUPRI,'(30X,A)') 'RPA(D) results:'
         WRITE(LUPRI,9001)
C
CSPAS:15/11-13: RPA(D) is now model 3 due to HRPA
         CALL EXCOUT(TRLEN(1,1,1,3),TRVEL(1,1,1,3),TQLEN(1,1,1,1,3),
     &               TQVEL(1,1,1,1,3),TRMAG(1,1,1,3),TRLON(1,1,1,3),
     &               BSRLON(1,1,1,3),TTLEN(1,1,1,3),EXENG(1,1,3),
CKeinSPASmehr
     &               FONAC,FONA2,DUMMY,DUMMY,DUMMY,DUMMY,
     &               DUMMY,DUMMY,DUMMY,DUMMY,WORK,LWORK)
C
         DO ISYM = 1, NSYM
C
            NEXCIT(ISYM) = NEXCTX(ISYM)
C
         END DO
C
      END IF
C
C=================================================
C     Write calculated SOPPA properties to output.
C=================================================
C
      IF (AOSOP) THEN
C
         DO ISYM = 1, NSYM
C
            NEXCTX(ISYM) = NEXCIT(ISYM)
CPFP   2009-Jun: new triplet variables
            NEXCIT(ISYM) = MIN(NEXCIT(ISYM),(NT1AM(ISYM)+NT2AM(ISYM)))
            !NEXCIT(ISYM) = MIN(NEXCIT(ISYM),(NT1AM(ISYM)+N2P2HOP(ISYM))) ! todo triple excitations are not fully implemented yet
Cend-PFP
C
         END DO
C
         WRITE(LUPRI,9000)
         WRITE(LUPRI,'(30X,A)') 'SOPPA results:'
         WRITE(LUPRI,9001)
C
CSPAS:15/11-13: SOPPA is now model 4 due to HRPA
         CALL EXCOUT(TRLEN(1,1,1,4),TRVEL(1,1,1,4),TQLEN(1,1,1,1,4),
     &               TQVEL(1,1,1,1,4),TRMAG(1,1,1,4),TRLON(1,1,1,4),
     &               BSRLON(1,1,1,4),TTLEN(1,1,1,4),EXENG(1,1,4),
CKeinSPASmehr
     &               FONAC,FONA2,DUMMY,DUMMY,DUMMY,DUMMY,
     &               DUMMY,DUMMY,DUMMY,DUMMY,WORK,LWORK)
C
         DO ISYM = 1, NSYM
C
            NEXCIT(ISYM) = NEXCTX(ISYM)
C
         END DO
C
      END IF
C
C=======================================================
C     Write calculated SOPPA(CCSD) properties to output.
C=======================================================
C
      IF (AOSOC) THEN
C
         DO ISYM = 1, NSYM
C
            NEXCTX(ISYM) = NEXCIT(ISYM)
CPFP   2009-Jun: new triplet variables
            NEXCIT(ISYM) = MIN(NEXCIT(ISYM),(NT1AM(ISYM)+NT2AM(ISYM)))
            !NEXCIT(ISYM) = MIN(NEXCIT(ISYM),(NT1AM(ISYM)+N2P2HOP(ISYM))) ! todo triple excitations are not fully implemented yet
Cend-PFP
C
         END DO
C
         WRITE(LUPRI,9000)
         WRITE(LUPRI,'(30X,A)') 'SOPPA(CCSD) results:'
         WRITE(LUPRI,9001)
C
CSPAS:15/11-13: SOPPA(CCSD) is now model 5 due to HRPA
         CALL EXCOUT(TRLEN(1,1,1,5),TRVEL(1,1,1,5),TQLEN(1,1,1,1,5),
     &               TQVEL(1,1,1,1,5),TRMAG(1,1,1,5),TRLON(1,1,1,5),
     &               BSRLON(1,1,1,5),TTLEN(1,1,1,5),EXENG(1,1,5),
CKeinSPASmehr
     &               FONAC,FONA2,DUMMY,DUMMY,DUMMY,DUMMY,
     &               DUMMY,DUMMY,DUMMY,DUMMY,WORK,LWORK)
C
         DO ISYM = 1, NSYM
C
            NEXCIT(ISYM) = NEXCTX(ISYM)
C
         END DO
C
      END IF
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL FLSHFO(LUPRI)
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
