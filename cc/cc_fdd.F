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
*=====================================================================*
      SUBROUTINE CC_FDD(NC1VEC,NC2VEC,LISTA,ITAMPA,LISTB,ITAMPB,
     &                  LISTC,ITAMPC,TYAM,RESULT,WORK,LWORK)
C
C---------------------------------------------------------------------
C Test routine for calculating the CC D matrix by finite difference on
C the C matrix transformation.
C Ch. Haettig, maj 1997, based on Oves CCLR_FDF routine
C---------------------------------------------------------------------
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxorb.h>
#include <iratdef.h>
#include <ccorb.h>
#include <aovec.h>
#include <ccsdinp.h>
#include <cclr.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <leinf.h>
C
      DIMENSION WORK(LWORK),ITADR(2),RESULT(*)
      PARAMETER (XHALF = 0.5D00,XMTWO = -2.0D00, DELTA = 1.0D-07)
      PARAMETER (ONE = 1.0d0, ZERO = 0.0d0, TWO = 2.0d0)
      CHARACTER MODEL*10
      CHARACTER FILCMA*8
      LOGICAL L1TST,L2TST, LETST
      INTEGER ICTRAN(4)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      MODEL = 'CCSD      '
      IF (CCS) MODEL = 'CCS       '
      IF (CC2) MODEL = 'CC2       '
C
      IF (CCR12) CALL QUIT('Finite-difference D-matrix for CCR12 '//
     &                     'not adapted')
C
      IF (IPRINT.GT.5) THEN
         CALL AROUND( 'IN CC_FDD  : MAKING FINITE DIFF. CC D Matrix')
      ENDIF
C
C----------------------------
C     Work space allocations.
C----------------------------
C
      ISYMTR     = 1
      ISYMOP     = 1
C
      NTAMP      = NT1AM(ISYMTR) + NT2AM(ISYMTR)
      NTAMP2     = NTAMP*(NC1VEC + NC2VEC )
      KF         = 1
      KRHO1      = KF       + NTAMP2
      KRHO2      = KRHO1    + NT1AMX
      KC1AM      = KRHO2    + MAX(NT2AMX,NT2AM(ISYMTR))
      KC2AM      = KC1AM    + NT1AM(ISYMTR)
      KEND1      = KC2AM 
     *           + MAX(NT2AMX,NT2AM(ISYMTR),NT2AO(ISYMTR),
     *                 2*NT2ORT(ISYMTR))
      LWRK1      = LWORK    - KEND1
C
      KRHO1D     = KEND1
      KRHO2D     = KRHO1D   + NT1AMX
      KEND2      = KRHO2D     
     *           + MAX(NT2AMX,NT2AM(ISYMTR),NT2AO(ISYMTR),
     *                 2*NT2ORT(ISYMTR))
      LWRK2      = LWORK      - KEND1
C
      IF (IPRINT .GT. 100 ) THEN
         WRITE(LUPRI,*) ' IN CC_FDD: KF      =  ',KF     
         WRITE(LUPRI,*) ' IN CC_FDD: KRHO1   =  ',KRHO1
         WRITE(LUPRI,*) ' IN CC_FDD: KRHO2   =  ',KRHO2
         WRITE(LUPRI,*) ' IN CC_FDD: KC1AM   =  ',KC1AM
         WRITE(LUPRI,*) ' IN CC_FDD: KC2AM   =  ',KC2AM
         WRITE(LUPRI,*) ' IN CC_FDD: KRHO1D  =  ',KRHO1D
         WRITE(LUPRI,*) ' IN CC_FDD: KRHO2D  =  ',KRHO2D
         WRITE(LUPRI,*) ' IN CC_FDD: KEND2   =  ',KEND2
         WRITE(LUPRI,*) ' IN CC_FDD: LWRK2   =  ',LWRK2
      ENDIF
      IF (LWRK2.LT.0 ) THEN
         WRITE(LUPRI,*) 'Too little work space in CC_FDD '
         WRITE(LUPRI,*) 'AVAILABLE: LWORK   =  ',LWORK
         WRITE(LUPRI,*) 'NEEDED (AT LEAST)  =  ',KEND2
         CALL QUIT('TOO LITTLE WORKSPACE IN CC_FDD ')
      ENDIF
      KF2   = KF      + NC1VEC*NTAMP
C
C---------------------
C     Initializations.
C---------------------
C
      CALL DZERO(WORK(KC1AM),NT1AMX)
      CALL DZERO(WORK(KC2AM),NT2AMX)
      CALL DZERO(WORK(KF),NTAMP2)
      IF (ABS(DELTA) .GT. 1.0D-15 ) THEN 
         DELTAI = 1.0D00/DELTA
      ELSE
         DELTAI = 1
      ENDIF
      X11 = 0.0D00
      X12 = 0.0D00
      X21 = 0.0D00
      X22 = 0.0D00
      XNJ = 0.0D00
C
C------------------------------------------------
C     Read the CC reference amplitudes From disk.
C------------------------------------------------
C
      IOPT = 3
      CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KC1AM),WORK(KC2AM))
C
C----------------------------------------------
C     Save the CC reference amplitudes on disk.
C----------------------------------------------
C
      LUTAM = -1
      CALL GPOPEN(LUTAM,'TAM_SAV','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     *            .FALSE.)
      REWIND(LUTAM)
      WRITE(LUTAM) (WORK(KC1AM + I -1 ), I = 1, NT1AMX)
      WRITE(LUTAM) (WORK(KC2AM + I -1 ), I = 1, NT2AMX)
      CALL GPCLOSE(LUTAM,'KEEP')
C
      IF (IPRINT.GT.125) THEN
         RHO1N = DDOT(NT1AMX,WORK(KC1AM),1,WORK(KC1AM),1)
         RHO2N = DDOT(NT2AMX,WORK(KC2AM),1,WORK(KC2AM),1)
         WRITE(LUPRI,*) 'Norm of T1AM: ',RHO1N
         WRITE(LUPRI,*) 'Norm of T2AM: ',RHO2N
         CALL CC_PRP(WORK(KC1AM),WORK(KC2AM),1,1,1)
      ENDIF
      RSPIM = .TRUE.
C
C------------------------------------------
C     Calculate reference A*T vector.
C------------------------------------------
C
      ICTRAN(1) = ITAMPA
      ICTRAN(2) = ITAMPB
      ICTRAN(3) = ITAMPC
      ICTRAN(4) = 0

      NCTRAN = 1
      IOPT   = 1
      FILCMA = 'CCCMAT'
      
      CALL CC_CMAT(ICTRAN,NCTRAN,LISTA,LISTB,LISTC,IOPT,
     &             FILCMA,IDUM,RDUM,0,WORK(KRHO1D),LWORK-KRHO1D)

C
C-------------------------
C     Zero out components.
C-------------------------
C
      IF (LCOR .OR. LSEC) THEN
C
         CALL CC_CORE(WORK(KRHO1D),WORK(KRHO2D),ISYMTR)
C
      ENDIF
C
      IF (IPRINT.GT.2) THEN
         RHO1N = DDOT(NT1AMX,WORK(KRHO1D),1,WORK(KRHO1D),1)
         RHO2N = DDOT(NT2AMX,WORK(KRHO2D),1,WORK(KRHO2D),1)
         WRITE(LUPRI,*) 'Norm of RHO1: ',RHO1N,'ref'
         WRITE(LUPRI,*) 'Norm of RHO2: ',RHO2N,'ref'
      ENDIF
      IF (IPRINT.GT.125) THEN
         CALL CC_PRP(WORK(KRHO1D),WORK(KRHO2D),1,1,1)
      ENDIF

      CALL DCOPY(NT1AMX,WORK(KRHO1D),1,WORK(KRHO1),1)
      CALL DCOPY(NT2AMX,WORK(KRHO2D),1,WORK(KRHO2),1)
C
C=============================================
C     Calculate C-matrix by finite difference.
C=============================================
C
      DO 100 I = 1, NC1VEC
         WRITE (LUPRI,*) 'singles index:',I
C
C----------------------------------------
C        Add finite displadement to t and 
C        calculate new intermediates.
C----------------------------------------
C
         LUTAM = -1
         CALL GPOPEN(LUTAM,'TAM_SAV','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     *               .FALSE.)
         READ(LUTAM) (WORK(KC1AM + J -1 ) , J = 1, NT1AMX)
         READ(LUTAM) (WORK(KC2AM + J -1 ) , J = 1, NT2AMX)
         CALL GPCLOSE(LUTAM,'KEEP')
C
         TI   = SECOND()
         WORK(KC1AM +I -1) = WORK(KC1AM +I -1 ) + DELTA
         IF (LCOR .OR. LSEC) THEN
            CALL CC_CORE(WORK(KC1AM),WORK(KC2AM),ISYMTR)
         ENDIF
C
         IOPT = 3
         CALL CC_WRRSP('R0',0,1,IOPT,MODEL,DUMMY,WORK(KC1AM),
     *                 WORK(KC2AM),WORK(KEND2),LWRK2)
C
         RSPIM = .TRUE.
         CALL CCRHSN(WORK(KRHO1D),WORK(KRHO2D),WORK(KC1AM),
     *               WORK(KC2AM),WORK(KEND2),LWRK2,'XXX') 
C
C------------------
C        Transform.
C------------------
C
         ICTRAN(1) = ITAMPA
         ICTRAN(2) = ITAMPB
         ICTRAN(3) = ITAMPC
         ICTRAN(4) = 0

         NCTRAN = 1
         IOPT   = 1
         FILCMA = 'CCCMAT'
      
         CALL CC_CMAT(ICTRAN,NCTRAN,LISTA,LISTB,LISTC,IOPT,
     &                FILCMA,IDUM,RDUM,0,WORK(KRHO1D),LWORK-KRHO1D)


         IF (LCOR .OR. LSEC) THEN
            CALL CC_CORE(WORK(KRHO1D),WORK(KRHO2D),ISYMTR)
         ENDIF
C
         IF (IPRINT.GT.2) THEN
            RHO1N = DDOT(NT1AMX,WORK(KRHO1D),1,WORK(KRHO1D),1)
            RHO2N = DDOT(NT2AMX,WORK(KRHO2D),1,WORK(KRHO2D),1)
            WRITE(LUPRI,*) 'Norm of RHO1: ',RHO1N,'ai=',I
            WRITE(LUPRI,*) 'Norm of RHO2: ',RHO2N,'ai=',I
         ENDIF
         IF (IPRINT.GT.125) THEN
            CALL CC_PRP(WORK(KRHO1D),WORK(KRHO2D),1,1,1)
         ENDIF
         CALL DAXPY(NT1AMX,-1.0D00,WORK(KRHO1),1,WORK(KRHO1D),1)
         CALL DAXPY(NT2AMX,-1.0D00,WORK(KRHO2),1,WORK(KRHO2D),1)
         CALL DSCAL(NT1AMX,DELTAI,WORK(KRHO1D),1)
         CALL DSCAL(NT2AMX,DELTAI,WORK(KRHO2D),1)
         CALL DCOPY(NT1AMX,WORK(KRHO1D),1,
     *              WORK(KF+NTAMP*(I-1)),1)
         CALL DCOPY(NT2AMX,WORK(KRHO2D),1,
     *             WORK(KF+NTAMP*(I-1)+NT1AMX),1)
         X11 = X11 + DDOT(NT1AMX,WORK(KRHO1D),1,WORK(KRHO1D),1)
         X21 = X21 + DDOT(NT2AMX,WORK(KRHO2D),1,WORK(KRHO2D),1)
C
         TI   = SECOND() - TI
         IF (IPRINT.GT.5 ) THEN
            WRITE(LUPRI,*) '  '
            WRITE(LUPRI,*) 'FDB ROW NR. ',I,' DONE IN ',TI,' SEC.'
         ENDIF
C
 100  CONTINUE
C
C----------------------------------------------------------------
C     Loop over T2 amplitudes. Take care of diagonal t2 elements
C     is in a different convention in the energy code.
C     Factor 1/2 from right , and factor 2 from left.
C----------------------------------------------------------------
C
      DO 200 NAI = 1, NT1AMX
        DO 300 NBJ = 1, NAI
         I = INDEX(NAI,NBJ)
C
         IF (I.LE.NC2VEC) THEN
           WRITE (LUPRI,*) 'doubles index:',I
C
C--------------------------------------------
C          Add finite displacement to t and
C          calculate new intermediates.
C-------------------------------------------
C
           LUTAM = -1
           CALL GPOPEN(LUTAM,'TAM_SAV','UNKNOWN',' ','UNFORMATTED',
     *                 IDUMMY,.FALSE.)
           READ(LUTAM) (WORK(KC1AM + J -1 ) , J = 1, NT1AMX)
           READ(LUTAM) (WORK(KC2AM + J -1 ) , J = 1, NT2AMX)
           CALL GPCLOSE(LUTAM,'KEEP')
C
           TI   = SECOND()
           DELT = DELTA
           IF (NAI.EQ.NBJ) DELT = 2*DELTA
           WORK(KC2AM + I -1) = WORK(KC2AM+I -1) + DELT
           IF (LCOR .OR. LSEC) THEN
             CALL CC_CORE(WORK(KC1AM),WORK(KC2AM),ISYMTR)
           ENDIF
C
           IOPT = 3
           CALL CC_WRRSP('R0',0,1,IOPT,MODEL,DUMMY,WORK(KC1AM),
     *                   WORK(KC2AM),WORK(KEND2),LWRK2)
C
           RSPIM = .TRUE.
           CALL CCRHSN(WORK(KRHO1D),WORK(KRHO2D),WORK(KC1AM),
     *                 WORK(KC2AM),WORK(KEND2),LWRK2,'XXX') 
C
C--------------------
C          Transform.
C--------------------
C
           ICTRAN(1) = ITAMPA
           ICTRAN(2) = ITAMPB
           ICTRAN(3) = ITAMPC
           ICTRAN(4) = 0

           NCTRAN = 1
           IOPT   = 1
           FILCMA = 'CCCMAT'
      
           CALL CC_CMAT(ICTRAN,NCTRAN,LISTA,LISTB,LISTC,IOPT,
     &                  FILCMA,IDUM,RDUM,0,WORK(KRHO1D),LWORK-KRHO1D)


           IF (LCOR .OR. LSEC) THEN
              CALL CC_CORE(WORK(KRHO1D),WORK(KRHO2D),ISYMTR)
           ENDIF
C
           IF (IPRINT.GT.2) THEN
             RHO1N = DDOT(NT1AMX,WORK(KRHO1D),1,WORK(KRHO1D),1)
             RHO2N = DDOT(NT2AMX,WORK(KRHO2D),1,WORK(KRHO2D),1)
             WRITE(LUPRI,*) 'Norm of RHO1: ',RHO1N,'aibj=',I
             WRITE(LUPRI,*) 'Norm of RHO2: ',RHO2N,'aibj=',I
           ENDIF
           IF (IPRINT.GT.125) THEN
            CALL CC_PRP(WORK(KRHO1D),WORK(KRHO2D),1,1,1)
           ENDIF
C
           CALL DAXPY(NT1AMX,-1.0D00,WORK(KRHO1),1,WORK(KRHO1D),1)
           CALL DAXPY(NT2AMX,-1.0D00,WORK(KRHO2),1,WORK(KRHO2D),1)
           CALL DSCAL(NT1AMX,DELTAI,WORK(KRHO1D),1)
           CALL DSCAL(NT2AMX,DELTAI,WORK(KRHO2D),1)
           CALL DCOPY(NT1AMX,WORK(KRHO1D),1,
     *              WORK(KF2+NTAMP*(I-1)),1)
           CALL DCOPY(NT2AMX,WORK(KRHO2D),1,
     *              WORK(KF2+NTAMP*(I-1)+NT1AMX),1)
C
           X12 = X12 + DDOT(NT1AMX,WORK(KRHO1D),1,WORK(KRHO1D),1)
           X22 = X22 + DDOT(NT2AMX,WORK(KRHO2D),1,WORK(KRHO2D),1)
           TI   = SECOND() - TI
           IF (IPRINT.GT.5 ) THEN
              WRITE(LUPRI,*) '  '
              WRITE(LUPRI,*) 'FDB ROW NR. ',I+NT1AMX,
     *                  ' DONE IN ',TI,' SEC.'
           ENDIF
C
         ENDIF
C
 300    CONTINUE
 200  CONTINUE
C
      WRITE(LUPRI,*) '    '
      WRITE(LUPRI,*) '**  FINITE DIFF WITH DELTA ',DELTA, '**'
      WRITE(LUPRI,*) '    '
      IF ((IPRINT .GT. 4).AND.(.TRUE.)) THEN
       CALL AROUND('FINITE DIFF. CC D*TA*TB*TC-Matrix - 11 & 21 PART' )
       CALL OUTPUT(WORK(KF),1,NTAMP,1,NC1VEC,NTAMP,NC1VEC,1,LUPRI)
       CALL AROUND('FINITE DIFF. CC D*TA*TB*TC-Matrix - 12 & 22 PART' )
       CALL OUTPUT(WORK(KF+NTAMP*NC1VEC),1,NTAMP,1,NC2VEC,
     *               NTAMP,NC2VEC,1,LUPRI)
      ENDIF
      IF (.TRUE.) THEN
       XNJ = X11 + X12 + X21 + X22
       WRITE(LUPRI,*)'  '
       WRITE(LUPRI,*)' NORM OF FIN. DIFF. D*TA*TB*TC-Matrix.', SQRT(XNJ)
       WRITE(LUPRI,*)'  '
       WRITE(LUPRI,*)' NORM OF 11 PART OF FD. D*TA*TB*TC-mat.: ',
     &      SQRT(X11)
       WRITE(LUPRI,*)' NORM OF 21 PART OF FD. D*TA*TB*TC-mat.: ', 
     &      SQRT(X21)
       WRITE(LUPRI,*)' NORM OF 12 PART OF FD. D*TA*TB*TC-mat.: ', 
     &      SQRT(X12)
       WRITE(LUPRI,*)' NORM OF 22 PART OF FD. D*TA*TB*TC-mat.: ', 
     &      SQRT(X22)
      ENDIF
C
C--------------------------------------
C     Calculate Matrix times Ty vector.
C--------------------------------------
C
      CALL DGEMV('N',NTAMP,NTAMP,ONE,WORK(KF),NTAMP,TYAM,1,
     *           ZERO,RESULT,1)

C--------------------------------------
C     scale diagonal with 1/2:
C--------------------------------------
C     CALL CCLR_DIASCL(RESULT(NT1AM(1)+1),TWO,1)

      WRITE (LUPRI,*) 'NTAMP:',NTAMP
      WRITE (LUPRI,*) 'NORM^2 OF TXAM VECTOR:',
     *   DDOT(NT1AM(1)+NT2AM(1),TXAM,1,TXAM,1)
      WRITE (LUPRI,*) 'NORM^2 OF TYAM VECTOR:',
     *   DDOT(NT1AM(1)+NT2AM(1),TYAM,1,TYAM,1)
      WRITE (LUPRI,*) 'NORM^2 OF RESULT VECTOR:',
     *   DDOT(NTAMP,RESULT,1,RESULT,1)

C
C-------------------------------------------------
C     Restore the CC reference amplitudes on disk.
C-------------------------------------------------
C
      LUTAM = -1
      CALL GPOPEN(LUTAM,'TAM_SAV','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     *            .FALSE.)
      REWIND(LUTAM)
      READ(LUTAM) (WORK(KC1AM + I -1 ) , I = 1, NT1AMX)
      READ(LUTAM) (WORK(KC2AM + I -1 ) , I = 1, NT2AMX)
      CALL GPCLOSE(LUTAM,'DELETE')
C
      IOPT = 3
      CALL CC_WRRSP('R0',0,1,IOPT,MODEL,DUMMY,WORK(KC1AM),
     *              WORK(KC2AM),WORK(KEND2),LWRK2)
C
      RSPIM = .TRUE.
      CALL CCRHSN(WORK(KRHO1D),WORK(KRHO2D),WORK(KC1AM),
     &            WORK(KC2AM),
     *            WORK(KEND2),LWRK2,'XXX') 
C
      IF (IPRINT .GT. 10) THEN
         CALL AROUND(' END OF CC_FDD ')
      ENDIF
C
      RETURN
      END
*=====================================================================*
