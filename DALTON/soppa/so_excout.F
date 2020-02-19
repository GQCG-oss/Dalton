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
      use so_info, only: sop_num_models, sop_mod_fullname,
     &                   sop_model_rpa, sop_model_rpad,
     &                   so_get_active_models, so_has_doubles,
     &                   sop_model_hrpad, sop_model_shrpad
#include "implicit.h"
#include "priunit.h"
C
#include "ccsdsym.h"
#include "ccorb.h"
#include "cbiexc.h"
#include "soppinf.h"
#include "maxorb.h"
C#include "secmom.h"
#include "codata.h"
C
CRF: Defer the last dimension (useful, so we need only allocate the
C    actual amount of space, that we use)
      DIMENSION TRLEN(3,NSYM,MXNEXI,*),
     &          TRVEL(3,NSYM,MXNEXI,*),
     &          TQLEN(3,3,NSYM,MXNEXI,*),
     &          TQVEL(3,3,NSYM,MXNEXI,*)
      DIMENSION TRLON(3,NSYM,MXNEXI,*),
     &          TRMAG(3,NSYM,MXNEXI,*)
      DIMENSION TTLEN(10,NSYM,MXNEXI,*)
      DIMENSION BSRLON(3,NSYM,MXNEXI,*),
     &          EXENG(NSYM,MXNEXI,*)
CClark:11/01/2016
      DIMENSION BETHE(3,LQ,*),STOPP(3,LVEL,2,*)
CClark:end
CKeinSPASmehr
      DIMENSION FONAC(*)
      DIMENSION FONA2(*)
      DIMENSION SECMAT(3,MXNEXI,NSYM)
      DIMENSION WORK(LWORK)
      DIMENSION NEXCTX(8)
C
      LOGICAL   ACTIVE_MODELS(sop_num_models),
     &          TREAT_DOUBLES
      INTEGER   IMODEL
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_EXCOUT')
C
C     Get list of active models
      CALL SO_GET_ACTIVE_MODELS( ACTIVE_MODELS)

C     Result pointer
      IOUT = 0
C
C===============================================
C     LOOP OVER METHODS
C===============================================
C
      DO IMODEL = 1, SOP_NUM_MODELS

C        Skip models not used in this calculation
         IF ( .NOT. ACTIVE_MODELS(IMODEL) ) CYCLE
         IOUT = IOUT + 1
C
C        Check if we there is a treatment of doubly excited states
C       - RPA(D) has cannot do this!
         TREAT_DOUBLES = SO_HAS_DOUBLES(IMODEL)
     &                   .AND. (.NOT.IMODEL.EQ.SOP_MODEL_RPAD)
     &                   .AND. (.NOT.IMODEL.EQ.SOP_MODEL_HRPAD)
     &                   .AND. (.NOT.IMODEL.EQ.SOP_MODEL_SHRPAD)
C
         DO ISYM = 1, NSYM
C
C  Adjust number of excitations down, to the maximum allowed by the
C  method/basis set
C
            IF (TREAT_DOUBLES) THEN
               MAX_EXCITATION = NT1AM(ISYM) + N2P2HOP(ISYM)
            ELSE
               MAX_EXCITATION = NT1AM(ISYM)
            ENDIF
            NEXCTX(ISYM) = NEXCIT(ISYM)
            NEXCIT(ISYM) = MIN(NEXCIT(ISYM),MAX_EXCITATION)
C
         END DO
C
C----------------------------------------------
C        Write calculated properties to output.
C----------------------------------------------
C
         WRITE(LUPRI,9000)
         WRITE(LUPRI,'(31X,A,A)')  TRIM(SOP_MOD_FULLNAME(IMODEL)),
     &                              ' results:'
         WRITE(LUPRI,9001)
C
         CALL EXCOUT(TRLEN(1,1,1,IOUT),TRVEL(1,1,1,IOUT),
     &               TQLEN(1,1,1,1,IOUT),TQVEL(1,1,1,1,IOUT),
     &               TRMAG(1,1,1,IOUT),TRLON(1,1,1,IOUT),
     &               BSRLON(1,1,1,IOUT),TTLEN(1,1,1,IOUT),
     &               EXENG(1,1,IOUT),
     &               FONAC,FONA2,DUMMY,DUMMY,DUMMY,DUMMY,
     &               DUMMY,DUMMY,DUMMY,DUMMY,WORK,LWORK)
C
         IF (IMODEL .EQ. SOP_MODEL_RPA)
     &            CALL RP_WRITE_EXTEND(EXENG(1,1,IOUT),SECMAT)
C
         IF (STOPPW) THEN
C
            CALL STOPP_WRITE(BETHE(1,1,IOUT),STOPP(1,1,1,IOUT) )
C
         ENDIF
CSPAS:26/01-08:
C     moved because otherwise we get a lot of non-sense output
C     for the orbital extend
         DO ISYM = 1, NSYM
C
            NEXCIT(ISYM) = NEXCTX(ISYM)
C
         END DO
C
      END DO
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
C
CClark:01/03/2016
C
      SUBROUTINE STOPP_WRITE(BETHE,STOPP)
C
#include "implicit.h"
#include "priunit.h"
#include "cbiexc.h"
#include "ccorb.h"
#include "pi.h"
C
      DIMENSION BETHE(3,LQ),STOPP(3,LVEL,2)
      REAL*8    LTWO
      PARAMETER ( D2 = 2.0D00 )
      PARAMETER ( D4 = 4.0D00 )
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
      WRITE(LUPRI,'(1X,A)')
     &   ' ------------------------------------------------------'
C
      DO IQ = 1, LQ
C
         WRITE(LUPRI,'(2X,G8.3,3(4X,G11.5))')
     &                QMIN+(IQ-1)*QSTEP,(BETHE(K,IQ),K=1,3)
C
      END DO
C
      WRITE(LUPRI,'(//,1X,A)')
     &   '       Stopping Power ( X Component )  [a.u.]           '
      WRITE(LUPRI,'(1X,A/)')
     &   '    ---------------------------------------------'
      WRITE(LUPRI,'(1X,A)')
     &   ' V[a.u.] Ekinet[KeV] Full  L1  + L2   = split [a.u.]'
      WRITE(LUPRI,'(1X,A)')
     &   ' -----------------------------------------------------------'
C
      DO IVEL = 1, LVEL
C
         VELOC = VMIN+(IVEL-1)*VSTEP
C
C 2.0*NRFHT is the number of electrons.
C
         IF (VELOC*2.0 .GE. QINP) THEN
C
            LTWO = PI*D4*LOG(VELOC*D2/QINP)*D2*NRHFT/(VELOC*VELOC)
C
         ELSE
C
            LTWO = 0
C
         ENDIF
C
         WRITE(LUPRI,
     &         '(1X,F5.2,1X,F7.1,1X,F8.4,5X,F7.3,2X,F7.3,2X,F8.4)')
     &          VELOC,VELOC*VELOC*24.98,
     &          STOPP(1,IVEL,1),STOPP(1,IVEL,2),LTWO,
     &          STOPP(1,IVEL,2)+LTWO
C
      END DO
C
      WRITE(LUPRI,'(//,1X,A)')
     &   '       Stopping Power ( Y Component )  [a.u.]           '
      WRITE(LUPRI,'(1X,A/)')
     &   '    ---------------------------------------------'
      WRITE(LUPRI,'(1X,A)')
     &   ' V[a.u.] Ekinet[KeV] Full  L1  + L2   = split [a.u.]'
      WRITE(LUPRI,'(1X,A)')
     &   ' -----------------------------------------------------------'
C
      DO IVEL = 1, LVEL
C
         VELOC = VMIN+(IVEL-1)*VSTEP
C
C 2.0*NRFHT is the number of electrons.
C
         IF (VELOC*2.0 .GE. QINP) THEN
C
            LTWO = PI*D4*LOG(VELOC*D2/QINP)*D2*NRHFT/(VELOC*VELOC)
C
         ELSE
C
            LTWO = 0
C
         ENDIF
C
         WRITE(LUPRI,
     &         '(1X,F5.2,1X,F7.1,1X,F8.4,5X,F7.3,2X,F7.3,2X,F8.4)')
     &          VELOC,VELOC*VELOC*24.98,
     &          STOPP(2,IVEL,1),STOPP(2,IVEL,2),LTWO,
     &          STOPP(2,IVEL,2)+LTWO
C
      END DO
C
      WRITE(LUPRI,'(//,1X,A)')
     &   '       Stopping Power ( Z Component )  [a.u.]           '
      WRITE(LUPRI,'(1X,A/)')
     &   '    ---------------------------------------------'
      WRITE(LUPRI,'(1X,A)')
     &   ' V[a.u.] Ekinet[KeV] Full  L1  + L2   = split [a.u.]'
      WRITE(LUPRI,'(1X,A)')
     &   ' -----------------------------------------------------------'
C
      DO IVEL = 1, LVEL
C
         VELOC = VMIN+(IVEL-1)*VSTEP
C
C 2.0*NRFHT is the number of electrons.
C
         IF (VELOC*2.0 .GE. QINP) THEN
C
            LTWO = PI*D4*LOG(VELOC*D2/QINP)*D2*NRHFT/(VELOC*VELOC)
C
         ELSE
C
            LTWO = 0
C
         ENDIF
C
         WRITE(LUPRI,
     &         '(1X,F5.2,1X,F7.1,1X,F8.4,5X,F7.3,2X,F7.3,2X,F8.4)')
     &          VELOC,VELOC*VELOC*24.98,
     &          STOPP(3,IVEL,1),STOPP(3,IVEL,2),LTWO,
     &          STOPP(3,IVEL,2)+LTWO
C
      END DO
C
      END
C
CClark:end
C
      SUBROUTINE RP_WRITE_EXTEND(EXENG,SECMAT)
C
C     The output of RPA orbital extend is put here to clean up the above
C     routine
C
#include "implicit.h"
#include "priunit.h"

C secmom needs maxorb
#include "maxorb.h"
C Need the SECOMO, SECVMO, SECGR array.
C Secmom is a common block that lives ONLY to
C pass info from rp_charge to here... find some other way to move this
C info?
#include "secmom.h"
C Need symmetry info (NSYM, NVIR, NRHF...)
#include "ccorb.h"
C Need NEXCIT
#include "cbiexc.h"
C Need XTEV
#include "codata.h"

      DIMENSION SECMAT(3,MXNEXI,NSYM),
     &          EXENG(NSYM,MXNEXI)
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
     &         XTEV * EXENG(ISYM,IEXCI),
     &         ( SECMAT(I,IEXCI,ISYM), I = 1,3), SECSUM
C
         END DO
C
      END DO


      RETURN
      END SUBROUTINE


