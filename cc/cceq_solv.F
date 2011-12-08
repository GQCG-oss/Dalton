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
C  /* Deck cceq_sol */
      SUBROUTINE CCEQ_SOL(LIST,LPROJECT,ISTATPRJ,ECURR,
     *                    FRHO1,LUFR1,FRHO2,LUFR2,FRHO12,LUFR12,
     *                    FC1AM,LUFC1,FC2AM,LUFC2,FC12AM,LUFC12,
     *                    LREDS,REDS,FS12AM,LUFS12,FS2AM,LUFS2,
     *                    LINQCC,TRIPLET,ISIDE,ISTRVE,NVEC,NUPVEC,
     *                    NREDH,REDH,EIVAL,SOLEQ,
     *                    WRK,LWRK,APROXR12)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C 
C input
C
C   FRHO1,FRHO2,FRHO12,FC1AM,FC2AM,FC12AM are file names for files where
C   transformed vectors (rho1,rho2,rhor12) and trial vectors 
C   (c1am,c2am,cr12am) are stored.
C   LUFR1,LUFR2,LUFR12,LUFC1,LUFC2,LUFC12 are the corresponding unit 
C   numbers.
C
C   TRIPLET.EQ.T solve triplet equations
C          .EQ.F solve singlet equations
C   
C   ISIDE.EQ.1  SOLVE EQUATIONS FROM RIGHT
C   ISIDE.EQ.-1 SOLVE EQUATIONS FROM LEFT
C
C   LPROJECT = LOGICAL FOR PROJECTION, (sonia)
C   ISTATPRJ = STATE INDEX FOR PROJECTION (sonia)
C
C   LINQCC.EQ.T ,SOLVE SET OF NVEC LINEAR EQUATIONS
C   ISTRVE = START NUMBER IN THE LIST FOR SOLVING LINEAR EQUATIONS
C   NVEC   = NUMBER OF EQUATIONS THAT MUST BE SOLVED
C   NUPVEC = NUMBER OF CURRENT LINEAR INDEPENDENT NEW VECTORS.
C
C   EIVAL CONTAINS NVEC FREQUENCIES 
C   EIVAL( ISTRVE ... (ISTRVE+NVEC-1) )
C
C   LINQCC.EQ.F ,SOLVE EIGENVALUE EQUATIONS FOR NVEC LOWEST ROOTS
C
C IF (LINQCC) THEN
C
C    THIS SUBROUTINE DIRECT THE SOLUTIONS OF NVECS SIMULTANEOUS
C    SET OF NEWTON-RAPHSON EQUATIONS
C
C    IF (ISIDE.EQ.1) THEN
C
C       (A-EIVAL*I)*X(J) + Fright(J) = 0
C
C    ELSE IF ( ISIDE.EQ.-1) THEN
C
C       [ (A-EIVAL*I) ] T *X(J) + Fleft(J) = 0
C       ( the transposed reduced matrix is constructed, this is
C         done without change in code because the linear 
C         transformations are from the left, the reduced space is
C         set up as
C         REDH(I,J) = B(I) * S(J) )
C
C    END IF
C ELSE
C
C    SOLVE for NVECS lowest eigenvalues
C
C    IF (ISIDE.EQ.1) THEN
C
C       (A-EIVAL*I)*X(J)  = 0
C
C    ELSE IF ( ISIDE.EQ.-1) THEN
C
C       [ (A-EIVAL*I) ] T *X(J)  = 0
C       ( the transposed reduced matrix is constructed, this is
C         done without change in code because the linear 
C         transformations are from the left)
C
C    END IF
C
C END IF
C
C
C----------------------------------------------------------------
C CC_TRDRV carry out linear transformations on new trial vectors 
C          on LUFC1,LUFC2,LUFC12 and return linear transformed 
C          trial vectors on LUFR1,LUFR2,LUFR12
C CCRED finds the solution in reduced space.
C CCNEX finds residual and the next trial vectors (saved on disk)
C
C Work space flow:
C  The whole work space is released after each routine is called
C  all communication takes place via file with trial vectors
C  and linear transformed vectors on LUFC1,LUFC2,LUFC12 and 
C  LUFR1,LUFR2,LUFR12 respectively
C
C  MAXLE:  MAXIMUM NUMBER OF MICROITERATIONS
C          from common LEINF
C  MAXRED  MAXIMUM DIMENSION OF REDUCED SPACE from common CCLR   
C----------------------------------------------------------------
C
C----------------------------------------------------------------
C
C JK&OC 080802:
C Dirty hack to allow small number of iterations for CCSLV/MM calc.
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include <implicit.h>
#include <priunit.h>
#include <ccsdinp.h>
#include <ccsections.h>
#include <cclr.h>
#include <leinf.h>
#include <ccslvinf.h>
Cholesky
#include <maxorb.h>
#include <ccdeco.h>
C
      PARAMETER (TOL = 1.0D-14)
Cholesky
      PARAMETER (D0 = 0.0D0 )
C
      CHARACTER*(*) FRHO1,FRHO2,FC1AM,FC2AM,FRHO12,FC12AM,FS12AM,FS2AM
      CHARACTER*(*) LIST
      CHARACTER*3 APROXR12
      DIMENSION REDH(*), SOLEQ(*), REDS(*)
      DIMENSION EIVAL(*)
      DIMENSION WRK(*)
      LOGICAL LINQCC, LPROJECT, TRIPLET, LREDS
C
      CALL QENTER('CCEQ_SOL')
      IF (IPRLE .GT. 3 ) THEN
         WRITE(LUPRI,*) 'CCEQ_SOL: IPRLE   =    ',IPRLE
         WRITE(LUPRI,*) 'CCEQ_SOL: TRIPLET =    ',TRIPLET
         WRITE(LUPRI,*) 'CCEQ_SOL: LINQCC  =    ',LINQCC
         WRITE(LUPRI,*) 'CCEQ_SOL: NREDH   =    ',NREDH
         WRITE(LUPRI,*) 'CCEQ_SOL: MAXRED  =    ',MAXRED
         WRITE(LUPRI,*) 'CCEQ_SOL: MAXLE   =    ',MAXLE
         WRITE(LUPRI,*) 'CCEQ_SOL: ISIDE   =    ',ISIDE
         WRITE(LUPRI,*) 'CCEQ_SOL: ISTRVE  =    ',ISTRVE
         WRITE(LUPRI,*) 'CCEQ_SOL: NVEC    =    ',NVEC
         WRITE(LUPRI,*) 'CCEQ_SOL: NUPVEC  =    ',NUPVEC
         WRITE(LUPRI,*) 'CCEQ_SOL: NCCVAR  =    ',NCCVAR
      ENDIF
C
      CALL GETTIM(CSTR,WSTR)
C
      TIMMIC = SECOND()
      TIMLIN = D0
      TIMRED = D0
      TIMNEX = D0
      ITLE   = 0
C
Cholesky
C
C     Check that this is N-R for Cholesky CC2: Eigenvalue equations
C     are non-linear and cannot be determined here...
C
      IF ((CHOINT .AND. CC2) .AND. (.NOT. LINQCC)) THEN
         WRITE(LUPRI,'(//,A)')
     &   'CCEQ_SOL: Cholesky CC2 eigenvalues not available:'
         WRITE(LUPRI,'(A,A,/,A,/)')
     &  'The eigenvalue equations for Cholesky CC2 are ',
     &  'non-linear and thus outside the scope of CCEQ_SOL'
        CALL QUIT('CCEQ_SOL: no eigenvalues for Cholesky CC2')
      ENDIF
C
C     Check that all frequencies in EIVAL are identical for Cholesky CC2.
C
      IF (CHOINT .AND. CC2.AND. (NVEC.GT.0)) THEN
         ICOUN = 0
         REFRQ = EIVAL(1)
         DO I = 2,NVEC
            DIFF = DABS(EIVAL(I) - REFRQ)
            IF (DIFF .GT.  TOL) ICOUN = ICOUN + 1
         ENDDO
         IF (ICOUN .GT. 0) THEN
            WRITE(LUPRI,'(//,A)')
     &      'CCEQ_SOL: frequencies MUST be identical for all',
     &      ' perturbations for Cholesky CCC2.'
            WRITE(LUPRI,'(A)')
     &      'Frequencies passed to CCEQ_SOL in EIVAL array:'
            WRITE(LUPRI,'(4D20.12)') (EIVAL(I), I = 1,NVEC)
            CALL QUIT('CCEQ_SOL: frequency error for Cholesky CC2')
         ENDIF
      ENDIF
C
Cholesky
C
C     print banner for solver:
C
      IF (IPRLE.GT.0) THEN
         WRITE (LUPRI,'(///A/)') 
     &        ' >>>>> COUPLED CLUSTER RESPONSE SOLVER <<<<<'
Chol
         IF (CHOINT .AND. CC2.AND. (NVEC.GT.0))
     &       WRITE(LUPRI,'(A,F12.5,/)')
     &            ' Frequnecy :',EIVAL(1)
Chol
         WRITE (LUPRI,'(3X,A)')
     &        ' Iter  #Vectors  time (min)   residual'
         WRITE (LUPRI,'(3X,A)')
     &        ' --------------------------------------'
      END IF 
C
C
C     Loop over micro iterations. 
C
C
 100  CONTINUE 
        ITLE = ITLE + 1
C
        IF (IPRLE .GE. 2) TIMIT = SECOND()
C
         TIM    = SECOND()
C
C
C--------------------------------
C LINEAR TRANSFORMATIONS ON NUPVEC TRIAL VECTORS ON LUB
C THE LINEAR TRANSFORMED VECTORS ARE RETURNED ON LUS
C--------------------------------
C
         IF (IPRLE .GT. 5) THEN
            WRITE (LUPRI,*)
     *      ' >>> Call CC_TRDRV with TRIPLET and ECURR set to:',
     *      TRIPLET, ECURR
         ENDIF                      
C
Cholesky     Pass EIVAL (before it was DUMMY) to CC_TRDRV
Cholesky     to be used by CC_CHOATR
C
         IST = NREDH - NUPVEC + 1
         CALL CC_TRDRV(ECURR,FRHO1,LUFR1,FRHO2,LUFR2,FRHO12,LUFR12,
     *                       FC1AM,LUFC1,FC2AM,LUFC2,FC12AM,LUFC12,
     *                 TRIPLET,.FALSE.,IDUMMY,EIVAL,FS12AM,LUFS12,
     *                 FS2AM,LUFS2,
     *                 ISIDE,IST,NUPVEC,WRK,LWRK,APROXR12)
C
         CALL FLSHFO(LUPRI)
         TIM    = SECOND() - TIM
         TIMLIN =  TIMLIN + TIM
         NUPVECSAVE = NUPVEC
C
C
C--------------------------------
C        CALL CCRED(), INCREASE DIMENSION OF REDUCED SPACE,
C        AND FIND SOLUTIONS
C--------------------------------
C
         TIMRED = TIMRED - SECOND()

         CALL CCRED(FRHO1,LUFR1,FRHO2,LUFR2,FRHO12,LUFR12,
     *              FC1AM,LUFC1,FC2AM,LUFC2,FC12AM,LUFC12,
     *              CCR12,FS12AM,LUFS12,FS2AM,LUFS2,
     *              LINQCC,TRIPLET,ISIDE,LIST,
     *              LPROJECT,ISTATPRJ,LREDS,REDS,
     *              REDH,NREDH,ISTRVE,NUPVEC,NVEC,
     *              EIVAL,SOLEQ,WRK,LWRK,DEBUG)
         CALL FLSHFO(LUPRI)
         TIMRED = TIMRED + SECOND()
C
C
C        CREATE NVEC NEW LINEAR INDEPENDENT TRIAL VECTORS IN CCNEX()
C        FROM THE SOLUTIONS OF THE REDUCED  EQUATIONS
C        REDUCE THE NUMBER FOR CONVERGED VECTORS AND IF LINEAR 
C        DEPENDENCE IS ENCOUNTERED. 
C        NUPVEC IS NVEC REDUCED FOR CONVERGENCE AND
C        LINEAR DEPENDENCE
C
         IF (ITLE .LT. MAXLE) THEN
            JCONV = 0
         ELSE
            JCONV = -1
         END IF
C
C--------------------------------
C     Find the next trial vector.
C--------------------------------
C
         TIMNEX = TIMNEX - SECOND()
         CALL CCNEX(LIST,ITLE,LPROJECT,ISTATPRJ,
     *              FRHO1,LUFR1,FRHO2,LUFR2,FRHO12,LUFR12,
     *              FC1AM,LUFC1,FC2AM,LUFC2,FC12AM,LUFC12,
     *              FS12AM,LUFS12,FS2AM,LUFS2,
     *              LINQCC,TRIPLET,ISIDE,NREDH,ISTRVE,
     *              NVEC,NUPVEC,JCONV,EIVAL,SOLEQ,WRK,LWRK,RMXNORM)
C
C-----------------------------------------------------
C     Print timing & convergence statistic
C-----------------------------------------------------
C
         IF (IPRLE.GT.0) THEN
            WRITE (LUPRI,'(2X,I5,3X,I5,1X,F12.2,1X,E12.2)')
     *         ITLE,NUPVECSAVE,TIM/60.0D0,RMXNORM
         END IF
         CALL FLSHFO(LUPRI)
C
C-----------------------------------------------------
C     Print out the lowest eigenvalues, Ove 24-9-1995.
C-----------------------------------------------------
C
         IF (IPRLE .GT. 2 ) THEN
C
            IF (.NOT. LINQCC) THEN
               WRITE (LUPRI,'(/A)')
     *         ' Lowest eigenvalues in this iteration: '
               CALL OUTPUT(EIVAL,1,NVEC,1,1,NVEC,1,1,LUPRI) 
            ENDIF
C
            WRITE (LUPRI,'(A,I5)')
     *         ' Dimension of the reduced space: ',NREDH
C
         ENDIF
C
         CALL FLSHFO(LUPRI)
         TIMNEX = TIMNEX + SECOND()
C
         IF (IPRLE .GE. 2) THEN
            TIMIT = SECOND() - TIMIT
            WRITE (LUPRI,'(A,F12.2,/)')
     *         ' >>> TIME USED IN THIS MICRO ITERATION',TIMIT
         END IF
C
         IF (JCONV.LT.0) THEN
C           ( LINEAR DEPENDENCE BETWEEN ALL NEW TRIAL VECTORS )
            WRITE (LUPRI,'(/A/A)')
     *         ' ***  CCEQ_SOL -MICROITERATIONS STOPPED  ',
     *         '     LINEAR DEPENDENCE BETWEEN ALL NEW TRIAL VECTORS'
CCN           IF (RMXNORM.GT.THRLE) THEN 
CCN             WRITE(LUPRI,*) 
CCN    *           '    BUT NOT ALL SOLUTION VECTORS ARE'//
CCN    *           ' CONVERGED!'
CCN             CALL QUIT(' SOLUTION VECTORS NOT CONVERGED ')
CCN           END IF
         ELSE IF (JCONV.GT.0) THEN
C           (CONVERGED)
            IF (IPRLE .GT. 10) THEN
               WRITE(LUPRI,'(/,(A,I4,A,A,I4,A))')
     *         ' MICROITERATIONS CONVERGED FOR',NVEC,
     *         ' SOLUTION VECTORS',
     *         ' IN',ITLE,' ITERATIONS.'
            END IF
            IF (CCSLV) NSLVINIT = NSLVINIT + NREDH
            IF (CCSLV) LRSPFUL = .TRUE.
         ELSE IF (NREDH + NUPVEC .GT. MAXRED) THEN
C           (NOT CONVERGED)
               WRITE(LUPRI,'(/A,I5,A)')
     *          ' *** CCEQ_SOL-MAXIMUM DIMENSION OF REDUCED EQUATIONS',
     *            MAXRED,' REACHED.'
               WRITE(LUPRI,*)' ALL SOLUTION VECTORS NOT CONVERGED '
               CALL QUIT(' ALL SOLUTION VECTORS NOT CONVERGED ')
         ELSE IF (ITLE .LT. MAXLE) THEN
            CALL FLSHFO(LUPRI)
              IF (LINQCC.AND.CCSLV.AND.(ITLE.GE.MXLINIT).AND.
     *           (LIST(1:1) .EQ. 'L')) THEN
C             additions to reduced space should be ignored.
C
               NREDH = NREDH-NUPVEC
               NSLVINIT = NSLVINIT + NREDH
               LRSPFUL = .FALSE.
               WRITE(LUPRI,*)
     *        ' CCSLV: We stop for now though not fully converged yet'
                WRITE(LUPRI,*)
     *        ' Accumulated inner iterations are:', NSLVINIT
C               WRITE(LUPRI,*) 'last - NREDH,NUPVEC:',NREDH,NUPVEC
              ELSE
                GO TO 100
              END IF
         ELSE
            WRITE(LUPRI,'(/A,I5,A)')
     *         ' *** CCEQ_SOL-MAXIMUM NUMBER OF MICROITERATIONS',
     *         ITLE,' REACHED.'
               CALL QUIT(' *** CCEQ_SOL-MAX. MICROITERATIONS REACHED')
         END IF
C
C
C=====================
C     End of 100 Loop.
C=====================
C
      TIMMIC = SECOND() - TIMMIC
      CALL GETTIM(CEND,WEND)
      CTOT = CEND - CSTR
      WTOT = WEND - WSTR
C
      IF (IPRLE .GT. 0) THEN
         WRITE (LUPRI,'(3X,A)')
     &        ' --------------------------------------'
         WRITE (LUPRI,'(3X,A,I3,A)') ' converged in',ITLE,' iterations'
         WRITE (LUPRI,'(3X,A,E12.2)') ' threshold:',THRLE
         WRITE (LUPRI,'(//T10,A)') 'Routine          Time (min)'
         WRITE (LUPRI,'(T10,A)')  '---------------------------'
         WRITE (LUPRI,'(T10,A,F15.2)') 'CC_TRDRV ',TIMLIN/60.0D0
         WRITE (LUPRI,'(T10,A,F15.2)') 'CCRED    ',TIMRED/60.0D0
         WRITE (LUPRI,'(T10,A,F15.2)') 'CCNEX    ',TIMNEX/60.0D0
         WRITE (LUPRI,'(T10,A)') '---------------------------'
         WRITE (LUPRI,'(T10,A,F14.2,//)') 'Total time',TIMMIC/60.0D0
         CALL TIMTXT('>>>> Total CPU  time used in CCEQ_SOLV:',CTOT,
     &        LUPRI)
         CALL TIMTXT('>>>> Total wall time used in CCEQ_SOLV:',WTOT,
     &        LUPRI)
         WRITE (LUPRI,'(//A/)')
     &        ' >>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<'
      END IF
C
      IF (IPRLE .GT. 0 .AND. .NOT. LINQCC) THEN
         WRITE (LUPRI,'(//A/A/A,1P,D10.2,A/)')
     *      '  Final eigenvalues :',
     *      ' =====================',
     *      ' (convergence threshold:',THRLE,')'
         DO 920 I = 1,NVEC
            WRITE (LUPRI,'(I5,F15.10)') I,EIVAL(I)
  920    CONTINUE
      END IF
      IF (IPRLE .GE. 10) THEN
         IF (LINQCC) THEN
            WRITE (LUPRI,'(//A)')
     *         ' FINAL LINEAR SOLUTION VECTORS IN REDUCED BASIS FOR :'
         ELSE
            WRITE (LUPRI,'(//A)')
     *         ' FINAL EIGENVECTORS IN REDUCED BASIS FOR :'
         END IF
         IOFF = 0
         DO 930 I = 1,NVEC
            IF (LINQCC) THEN
               WRITE (LUPRI,'(/A,I4/)') ' - GRADIENT VECTOR NO.',I
            ELSE
               WRITE (LUPRI,'(/A,I4,A,1P,D15.6/A)')
     *            ' - EIGENVALUE NO.',I, ' - EIGENVALUE',EIVAL(I),
     *            ' - EIGENVECTOR :'
            END IF
            WRITE (LUPRI,'(5F15.8)') (SOLEQ(IOFF+J),J=1,NREDH)
            IOFF = IOFF + MAXRED
  930    CONTINUE
      END IF
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('CCEQ_SOL')
      RETURN
      END
C  /* Deck cceq_str */
      SUBROUTINE CCEQ_STR(FC1AM,LUFC1,FC2AM,LUFC2,FC12AM,LUFC12,
     *                    FS12AM,LUFS12,FS2AM,LUFS2,LPROJECT,ISTATPRJ,
     *                    TRIPLET,ISIDE,ISTRVE,NVEC,NUPVEC,
     *                    NREDH,EIVAL,IPLACE,WORK,LWORK,LIST)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C Based on original EVALST routine.
C
C Modified by Ove Christiansen to calculate diagonal when needed.
C Included projection, Sonia Coriani 1997
C 
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include <implicit.h>
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0 ,
     *            THRLDP = 1.0D-20 , THRDIA = 1.0D-6 )
#include <priunit.h>
#include <ibndxdef.h>
#include <ccorb.h>
#include <cclr.h>
#include <ccinf.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <leinf.h>
#include <ccinf.h>
#include <ccexci.h>
#include <dummy.h>
#include <r12int.h>
C
      LOGICAL LOCDBG, LNOREAD
      PARAMETER (LOCDBG = .FALSE.) 
      CHARACTER*(*)  FC1AM,FC2AM,FC12AM,FS12AM,FS2AM
      CHARACTER*(*) LIST
      CHARACTER*10  MODEL
      DIMENSION IPLACE(*) , WORK(*),EIVAL(*)
      LOGICAL NEWSTR, LPROJECT, TRIPLET
C
      CALL QENTER('CCEQ_STR')
      IF (IPRLE. GT. 10 ) THEN
          CALL AROUND( ' START OF CCEQ_STR ' )
          WRITE(LUPRI,*) 'CCEQ_STR: LIST     =',LIST(1:3)
          WRITE(LUPRI,*) 'CCEQ_STR: TRIPLET  =',TRIPLET
          WRITE(LUPRI,*) 'CCEQ_STR: LPROJECT =',LPROJECT
          WRITE(LUPRI,*) 'NCCVAR             =',NCCVAR
          CALL FLSHFO(LUPRI)
      ENDIF
C
      MODEL = 'CCSD      '
      IF (CCS) MODEL = 'CCS       '
      IF (CIS) MODEL = 'CIS       '
      IF (CC2) MODEL = 'CC2       '
      IF (CC3) MODEL = 'CC3       '
      IF (CC1A) MODEL = 'CC1A      '
C
      IF (TRIPLET) THEN
        IMULT = 3
      ELSE
        IMULT = 1
      ENDIF
C
C     ----------------------------
C     Initialize NUPVEC and NREDH:
C     ----------------------------
      NUPVEC = 0
      NREDH  = 0
C
      KV    = 1
      KEND1 = KV    + NCCVAR
      LEND1 = LWORK - KEND1
C
      KDIA  = KEND1
      KEND2 = KDIA   + NCCVAR
      LEND2 = LWORK  - KEND2
      IF (LEND2 .LT. 0) CALL QUIT('Too little workspace in CCEQ_STR-1')
C
C     -------------------------------------
C     Get an approximate diagonal Jacobian:
C     -------------------------------------
      IF (LIST(2:2) .EQ.'E') THEN
        CALL CCLR_DIA(WORK(KDIA),ISYMTR,TRIPLET,WORK(KEND2),LEND2)
        IF (LOCDBG) THEN
           CALL AROUND(' CCEQ_STR:  approximate diagonal Jacobian')
           CALL OUTPUT(WORK(KDIA),1,NCCVAR,1,1,NCCVAR,1,1,LUPRI)
        ENDIF
      ENDIF
C
C     ---------------------------------------------------------------
C     For eigenvalue problems select elements for unit start vectors:
C     ---------------------------------------------------------------
      IF (STVEC.AND.(LIST(2:2) .EQ.'E')) THEN
         WRITE (LUPRI,'(/A/A/A)')
     *      ' CC elements selected for start vectors:',
     *      ' =======================================',
     *      ' Start vector no.      CC element no.'
         DO I = 1,NVEC       
            IPLACE(I) = ISTVEC(I,ISYMTR)
            WRITE (LUPRI,'(I10,I20)') I,IPLACE(I)
         END DO
      ELSE IF (LIST(2:2) .EQ.'E') THEN
         MXELMN = MIN(NVEC+3,NCCVAR,MAXRED)
         CALL FNDMN3(WORK(KDIA),NCCVAR,MXELMN,IPLACE,
     *               NELMN,IPRLE,THRDIA)
         DO I = 1, NVEC
            IF (WORK(KDIA-1+IPLACE(I)) .GT. 1.0D6) THEN
              WRITE (LUPRI,*) 'Problem in CCEQ_STR:'
              WRITE (LUPRI,*) 'not enough allowed elements in vector:'
              WRITE (LUPRI,*) 'requested:',NVEC
              WRITE (LUPRI,*) 'found    :',I-1
              WRITE (LUPRI,*) 'reset NVEC variable to ',I-1
              NELMN = MIN(NELMN,I-1)
              NVEC  = MIN(NVEC,I-1)
            END IF
         END DO                  
      END IF
C
C--------------------------------------------------------
C     Add one vector at the time to list of startvectors.
C--------------------------------------------------------
C
      DO 100 I = 1,NVEC       
         NEWSTR = .TRUE.
C
C        ----------------------------------------------
C        try to restart from vectors available on file: 
C        ----------------------------------------------
         IF (CCRSTR) THEN
            IOPT   = 33
            ILSTNR = ISTRVE + I - 1
            ISYM   = ILSTSYM(LIST,ILSTNR)
            CALL DZERO(WORK(KV),NCCVAR)
            CALL CC_RDRSP(LIST,ILSTNR,ISYM,IOPT,MODEL,
     &                    WORK(KV),WORK(KV+NT1AM(ISYM)))
            IF (CCR12 .AND. IOPT.NE.33) THEN
              IOPTR12 = -32
              CALL CC_RDRSP(LIST,ILSTNR,ISYM,IOPTR12,MODEL,
     &                      DUMMY,WORK(KV+NT1AM(ISYM)+NT2AM(ISYM)))
            END IF
            ! check if readin of vector was succesfull - else calc. new.
            LEN_LIST = LEN(LIST)
            IF (IOPT .EQ. 33) THEN
               NEWSTR = .TRUE.
               WRITE(LUPRI,'(A,I3,3A,I3,A)') ' Vector nr.',ILSTNR,
     &         ' of type ',LIST(1:LEN_LIST),' and symmetry',ISYM,
     &         ' was not found on file.'
            ELSE
               NEWSTR = .FALSE.
               WRITE(LUPRI,'(A,I3,A,I3,A)') ' Vector nr.',ILSTNR,
     &         ' of symmetry',ISYM,' found on file - RESTART SUCCESS'
               WRITE(LUPRI,'(4A)') ' Start vector is a ',
     &         MODEL,LIST(1:3),' vector'
Casm &         MODEL,LIST(1:LEN_LIST),' vector'
            ENDIF
         ENDIF
C
C        -----------------------------------------------------------
C        for left eigenvectors try to start from right eigenvectors:
C        -----------------------------------------------------------
         IF (NEWSTR .AND. (LIST .EQ. 'LE')) THEN
            IOPT   = 33
            ILSTNR = ISTRVE + I - 1
            ISYM   = ILSTSYM(LIST,ILSTNR)
            WRITE(LUPRI,'(2A)') 'Try to use right eigenvector as start',
     &         ' for a left eigenvector.'
            CALL DZERO(WORK(KV),NCCVAR)
            CALL CC_RDRSP('RE ',ILSTNR,ISYM,IOPT,MODEL,
     &                    WORK(KV),WORK(KV+NT1AM(ISYM)))
            IF (CCR12 .AND. IOPT.NE.33) THEN
              IOPTR12 = -32
              CALL CC_RDRSP('RE ',ILSTNR,ISYM,IOPTR12,MODEL,
     &                      DUMMY,WORK(KV+NT1AM(ISYM)+NT2AM(ISYM)))
            END IF
            IF (IOPT .EQ. 33) THEN
               NEWSTR = .TRUE.
               WRITE(LUPRI,'(1X,A,I3,A,I3,A)') 'Vector nr.',ILSTNR,
     *         ' of type RE and symmetry',ISYM,' was ot found on file.'
            ELSE
               NEWSTR = .FALSE.
               WRITE(LUPRI,'(1X,A,I3,A,I3,A)') 'Vector nr.',ILSTNR,
     *         ' of symmetry',ISYM,' found on file - RESTART SUCCESS'
               WRITE(LUPRI,'(1X,A,A13,A)') 'Start vector is a ',
     *         MODEL//'RE',' vector'
            ENDIF
         ENDIF
C
C        ------------------------------------------------------------
C        if a restart was not possible generate a trial vector from
C        the gradient vector (for lin.eq.) or by choosing a unit vec.
C        ------------------------------------------------------------
         IF (NEWSTR) THEN

           IF (LIST(2:2) .EQ.'E') THEN
C             ---------------------------------------------------------
C             for excited state pick unit vector according to diagonal:
C             ---------------------------------------------------------
              WRITE(LUPRI,'(1X,A)')'Start vector guessed from diagonal'
              WRITE(LUPRI,'(1X,A,I3)')'... selected element no.',
     &              IPLACE(I)
              CALL DZERO(WORK(KV),NCCVAR)
              WORK(KV-1+IPLACE(I)) = D1
           ELSE 
C              --------------------
C              get gradient vector:
C              --------------------
               ILSTNR = ISTRVE + I - 1
               CALL CC_GETGD(WORK(KV),ILSTNR,ISIDE,LIST,
     *                       WORK(KEND1),LEND1)
               IF (IPRLE.GT.0) THEN
                 IF (LIST(1:2).EQ.'L0') THEN
                   WRITE(LUPRI,'(5X,A)') 
     *               'Start vector generated from gradient'
                 ELSE
                   WRITE(LUPRI,'(5X,2A,I3,A,I3,A)') LIST,
     *               ' start vector nr.',ILSTNR,
     *               ' of symmetry',ISYMTR,' generated from gradient'
                 END IF
               END IF

C              ------------------------------------------
C              Scale with diagonal.(including frequency).
C              (not for Cauchy vectors. CH.H. 20-11-97)
C              ------------------------------------------
               IF ( LIST(1:2).NE.'RC' .AND. LIST(1:2).NE.'LC' .AND.
     *              LIST(1:3).NE.'CRC'.AND. LIST(1:3).NE.'CLC' ) THEN
                 KOMEG1 = KV
                 KOMEG2 = KV + NT1AM(ISYMTR)
                 CALL CC_VSCAL(WORK(KOMEG1),WORK(KOMEG2),EIVAL(I),
     *                         WORK(KEND1),LEND1,ISYMTR)
               ENDIF
           ENDIF

         ENDIF
C
         IF ( DEBUG .OR. LOCDBG ) THEN
            RHO = DDOT(NCCVAR,WORK(KV),1,WORK(KV),1)
            WRITE(LUPRI,*) 'Norm of Start vector before CCORT: ',RHO
         ENDIF
C
         IF (LOCDBG) THEN
            CALL AROUND(' CCEQ_STR:VECTORS BEFORE CCORT/PRJDRV: ')
            WRITE(LUPRI,*) 'this is vector number:',I
            CALL OUTPUT(WORK(KV),1,NCCVAR,1,1,NCCVAR,1,1,LUPRI)
         ENDIF

         NFIN = 1

C
C-------------------------------------------------------------
C        if requested, project excited reference state out:
C -------------------------------------------------------------  
C
         IF (LPROJECT) THEN
            IF (TRIPLET) CALL QUIT('LPROJECT & TRIPLET not tested!')
            IF (DEBUG) THEN
               WRITE(LUPRI,*)' ---- CCEQ_STR: PROJECTION REQUIRED ----'
            END IF
            KPRJC = KEND2
            ISYMPR = ILSTSYM(LIST,ILSTNR)
            CALL PRJDRV(ISIDE,ISTATPRJ,ISYMPR,NFIN,NCCVAR,
     *                             WORK(KV),WORK(KPRJC)) 
         ENDIF
CCN
C
C        calculate metric for new b-vectors:
C
         IF (CCR12) THEN
           LNOREAD = .TRUE.
           IF (ISIDE .EQ. 1) THEN
             CALL CC_R12METRIC(ISYMTR,BRASCL,KETSCL,
     *                         WORK(KEND1),LEND1,DUMMY,DUMMY,
     *                         DUMMY,DUMMY,FS12AM,LUFS12,
     *                         FS2AM,LUFS2,NREDH+1,LNOREAD,WORK(KV))
           ELSE IF (ISIDE .EQ. -1) THEN
             CALL CC_R12METRIC(ISYMTR,0.5D0*KETSCL,2.0D0*BRASCL,
     *                         WORK(KEND1),LEND1,DUMMY,DUMMY,
     *                         DUMMY,DUMMY,FS12AM,LUFS12,
     *                         FS2AM,LUFS2,NREDH+1,LNOREAD,WORK(KV))
           END IF
         END IF
CCN
C
         CALL CCORT(LUFC1,FC1AM,LUFC2,FC2AM,LUFC12,FC12AM,
     *              LUFS2,FS2AM,LUFS12,FS12AM,
     *              TRIPLET,NFIN,NREDH,NCCVAR, 
     *              WORK(KV),ISYMTR,THRLDP,WORK(KEND1),LEND1,IPRLE)
         NUPVEC = NUPVEC + NFIN
         IF ( DEBUG .OR. LOCDBG ) THEN
            RHO = DDOT(NCCVAR,WORK(KV),1,WORK(KV),1)
            WRITE(LUPRI,*) 'Norm of Start vectors after CCORT:',RHO
         ENDIF
C
 100  CONTINUE
      IF (IPRLE .GT. 10)
     *   WRITE (LUPRI,'(/A,I5)') ' CCEQ_STR: NREDH =',NREDH
C
      IF (IPRLE. GT. 10 ) THEN
          CALL AROUND( ' END OF CCEQ_STR ' )
      ENDIF
C
      CALL QEXIT('CCEQ_STR')
      RETURN
      END
C  /* Deck ccconv */
      SUBROUTINE CCCONV(LUFC1,C1FIL,LUFC2,C2FIL,LUFC12,C12FIL,
     *                  TRIPLET,CCR12,NREDH,LST,LED,SOLEQ,RD,WRK)
C
C PURPOSE:
C  Calculate converged solutions to the LED sets of linear
C  equations in RD starting with equation number LEDST
C
#include <implicit.h>
#include <priunit.h>
#include <cclr.h>
#include <leinf.h>
C
      CHARACTER*(*) C1FIL,C2FIL,C12FIL
      LOGICAL TRIPLET,CCR12
C
      DIMENSION SOLEQ(MAXRED,*),RD(NCCVAR,*),WRK(*)
C
      CALL QENTER('CCCONV')
      CALL DZERO(RD,LED*NCCVAR)
      DO I = 1, NREDH
         CALL CC_GETVEC(LUFC1,C1FIL,LUFC2,C2FIL,LUFC12,C12FIL,
     *                  TRIPLET,I,WRK)
         DO K = 1, LED
            KVENU = LST -1 + K
            CALL DAXPY(NCCVAR,SOLEQ(I,KVENU),WRK,1,RD(1,K),1)
         END DO
      END DO
      CALL QEXIT('CCCONV')
      RETURN
      END
C  /* Deck ccnex */
      SUBROUTINE CCNEX(LIST,ITLE,LPROJECT,ISTATPRJ,
     *                 FRHO1,LUFR1,FRHO2,LUFR2,FRHO12,LUFR12,
     *                 FC1AM,LUFC1,FC2AM,LUFC2,FC12AM,LUFC12,
     *                 FS12AM,LUFS12,FS2AM,LUFS2,
     *                 LINQCC,TRIPLET,ISIDE,NREDH,ISTRVE,
     *                 NVEC,NUPVEC,JCONV,EIVAL,SOLEQ,WRK,LWRK,RMXNORM)
C
C PURPOSE: 1) Construct residual (A)*X(I) - EIVAL(I)*X(I) + F
C             for solution X(I) of reduced  equations
C             ( F is zero for eigenvalue equation )
C          2) Test for convergence of solutions
C             convergence criterium:
C             ||((A)*X(I) - EIVAL(I)*X(I) + F|| / ||X|| .LE. THRLE
C          3) Use generalized conjugate gradient algorithm
C             or davidson (for eigenvalue equation)
C             to determine next guess of trial vectors
C
C JCONV  input: if JCONV .lt. 0 do not calculate new trial vectors.
C       output: =  1 converged
C               =  0 not converged
C               = -1 not converged, linear dependency among all
C                                   trial vectors.
C
C RMXNORM : maximum norm of residuals
C
#include <implicit.h>
#include <priunit.h>
#include <ccsdinp.h>
#include <ccexci.h>
#include <ccsdsym.h>
#include <cclr.h>
#include <ccfro.h>
#include <leinf.h>
#include <r12int.h>
Cholesky
#include <maxorb.h>
#include <ccdeco.h>
Cholesky
      LOGICAL LOCDBG,LNOREAD
      PARAMETER ( LOCDBG = .FALSE. )
      PARAMETER ( THRLDP = 1.D-20 )
      PARAMETER ( DTEST = 1.0D-04 )
      PARAMETER ( DM1=-1.0D0,D1 =1.0D0, D0=0.0D0 )
C
      CHARACTER*(*) FRHO1,FRHO2,FRHO12,FC1AM,FC2AM,FC12AM,FS12AM
      CHARACTER*(*) LIST
      LOGICAL LINQCC,CCRSTRS,LPROJECT,TRIPLET
      DIMENSION EIVAL(*), SOLEQ(MAXRED,*),WRK(*)
C
C
C Space for CCNEX:
C
C MAXNEX: Maximum number of simultaneous vectors in CCNEX
C
      IF (DEBUG)  CALL AROUND(' Start of CCNEX ')

      MAXNEX = (LWRK-3*NCCVAR)/NCCVAR
      NUPVEC = 0
      NOTCON = 0
      DO 4000 ISIMC = 1,NVEC,MAXNEX
         NBX   = MIN(MAXNEX,(NVEC+1-ISIMC))
c        NBX   = MIN(MAXNEX,(NVEC+1-ISIMC),(NREDH+1-ISIMC))
C
C        CONSTRUCT RESIDUAL IN WRK(KRES)
C        RESIDUAL: (A)*X(I)-EIVAL(I)*X(I) +F
C
         KRES  = 1
         KTST  = KRES + NBX*NCCVAR 
         KWRK1 = KTST + NBX
         LWRK1 = LWRK - KWRK1

         IF (LINQCC) THEN
            LGD = KRES
            DO JR = 1,NBX 

C              ---------------------
C              Get gradient vectors:
C              ---------------------
               IVENU  = ISIMC-1+JR
               ILSTNR = IVENU + ISTRVE - 1

               CALL CC_GETGD(WRK(LGD),ILSTNR,ISIDE,LIST,
     *                       WRK(KWRK1),LWRK1)
C
C              ----------------------
C              Project from gradient:
C              ----------------------
               IF (LPROJECT) THEN
                  ISYPR = ILSTSYM(LIST,ILSTNR)
                  IF (DEBUG) WRITE(LUPRI,*)
     &                    '--- CCNEX projection from gradient --- '
                  CALL PRJDRV(ISIDE,ISTATPRJ,ISYPR,1,NCCVAR,
     *                            WRK(LGD),WRK(KWRK1))
               END IF
               LGD = LGD + NCCVAR
            END DO
         ELSE
            CALL DZERO(WRK(KRES),NBX*NCCVAR)
         END IF

         IF (IPRLE.GT.110 .OR. LOCDBG) THEN
             WRITE(LUPRI,*) 'CCNEX> RESIDUAL AFTER GD CONTRIBUTION'
             CALL OUTPUT(WRK(KRES),1,NCCVAR,1,NBX,NCCVAR,NBX,1,LUPRI)
         END IF
C
C        -----------------------------------------------------------
C        add -EIVAL(I)*X(I), where X(I) is the i'th solution vector:
C        -----------------------------------------------------------
         ISYM = ILSTSYM(LIST,ISTRVE)
         NVARPT = NCCVAR + 2*NALLAI(ISYM)
         MAXBVE = ( LWRK1 - NVARPT )/NCCVAR
         DO 60 JR = 1,NBX,MAXBVE
            NBVEC = MIN(MAXBVE,NBX+1-JR)
            KB = KWRK1
            KWRK2 = KB   + NCCVAR*NBVEC
            LWRK2 = LWRK - KWRK2
            IF (LWRK2 .LT. 0 ) THEN 
               CALL QUIT('Too little work in CCNEX xx')
            END IF

C           ----------------------------------------------------------
C           Construct the solution vectors in full space:
C             Save norm of solution vector in WRK(KTST) and the
C             vectors itself for restart on the files CC//LIST//{ivec}
C           ----------------------------------------------------------
            IVECNU = ISIMC -1 + JR
            CALL CCCONV(LUFC1,FC1AM,LUFC2,FC2AM,LUFC12,FC12AM,TRIPLET,
     *               CCR12,NREDH,IVECNU,NBVEC,SOLEQ,WRK(KB),WRK(KWRK2))

            DO IBVEC = 1, NBVEC
               WRK(KTST+JR-1+IBVEC-1) = 
     *         DNORM2(NCCVAR,WRK(KB+(IBVEC-1)*NCCVAR),1)
               IVEC = ISIMC - 1 + JR -1 + IBVEC   
               CALL CC_SAVE(WRK(KB+(IBVEC-1)*NCCVAR),IVEC+ISTRVE-1,
     *                      LIST,WRK(KWRK2),LWRK2)
            END DO

            IF (IPRLE.GT.105 .OR. LOCDBG) THEN
             WRITE (LUPRI,'(/A)')' CCNEX: solution vectors'
             CALL OUTPUT(WRK(KB),1,NCCVAR,1,NBVEC,NCCVAR,NBVEC,1,LUPRI)
             CALL PROVLP(WRK(KB),WRK(KB),NCCVAR,NBVEC,WRK(KWRK2),LUPRI)
            END IF

C           ----------------------------------------------------------
C           Construct S * X in full space and add - Eval * S * X to
C           the residual vectors:
C             Note for conventional CC the metric S is a unit matrix
C             and thus S * X are the solution vectors already
C             constructed above. Only for CCR12 we need to construct
C             S * X from the transformations with the metric matrix
C           ----------------------------------------------------------
            IF (CCR12) THEN
             IVECNU = ISIMC -1 + JR
             IF (IANR12.EQ.1 .OR. IANR12.EQ.3) THEN
               CALL CCCONV(LUFC1,FC1AM,LUFC2,FC2AM,LUFS12,FS12AM,
     *                     TRIPLET,CCR12,NREDH,IVECNU,NBVEC,SOLEQ,
     *                     WRK(KB),WRK(KWRK2))
             ELSE IF (IANR12.EQ.2) THEN
               CALL CCCONV(LUFC1,FC1AM,LUFS2,FS2AM,LUFS12,FS12AM,
     *                     TRIPLET,CCR12,NREDH,IVECNU,NBVEC,SOLEQ,
     *                     WRK(KB),WRK(KWRK2))
             ELSE
               CALL QUIT('This R12 ansatz is not implemented yet.')
             END IF
            END IF

            DO IBVEC = 1,NBVEC
               IVEC = ISIMC - 1 + JR -1 + IBVEC   
               IVEOFF = JR -1 + IBVEC -1
               CALL DAXPY(NCCVAR,-EIVAL(IVEC),WRK(KB+(IBVEC-1)*NCCVAR),
     *                    1,WRK(KRES+IVEOFF*NCCVAR),1)
            END DO
C
 60      CONTINUE

         IF (IPRLE.GT.110 .OR. LOCDBG) THEN
           WRITE(LUPRI,*)' RESIDUAL AFTER - w * S * X CONTRIBUTION'
           CALL OUTPUT(WRK(KRES),1,NCCVAR,1,NBX,NCCVAR,NBX,1,LUPRI)
         END IF

C        -----------------------------------------------
C        Add A * X(I) where X(I) is the solution to the
C        I'th set of Newton-Raphson equations
C        -----------------------------------------------
         DO K = 1,NREDH
            CALL CC_GETVEC(LUFR1,FRHO1,LUFR2,FRHO2,LUFR12,FRHO12,
     *                     TRIPLET,K,WRK(KWRK1))
            DO JR = 1,NBX
               JROOTJ = ISIMC - 1 + JR
               EVAL1  = SOLEQ(K,JROOTJ)
               CALL DAXPY(NCCVAR,EVAL1,WRK(KWRK1),1,
     *                       WRK(KRES+(JR-1)*NCCVAR),1)
            END DO
         END DO

C        -------------------------------
C        Residual is now in WRK(KRES) 
C        -------------------------------
         IF (IPRLE.GT.45 .OR. LOCDBG) THEN
            WRITE (LUPRI,'(/A)') ' CCNEX: residual vectors '
            CALL OUTPUT(WRK(KRES),1,NCCVAR,1,NBX,NCCVAR,NBX,1,LUPRI)
            CALL PROVLP(WRK(KRES),WRK(KRES),NCCVAR,NBX,WRK(KWRK1),LUPRI)
         END IF

C        -----------------------------------------
C        Projection of eigenvectors from residual:
C        -----------------------------------------
         IF (DEBUG) WRITE (LUPRI,*) 'CCNEX> LPROJECT = ', LPROJECT
         IF (LPROJECT) THEN
           KPRJ   = KWRK2
           ISYMPR = ILSTSYM(LIST,ISTRVE)
           IF (DEBUG) THEN
              WRITE(LUPRI,*)'-- CCNEX projection from  residual -- '
              WRITE(LUPRI,*)' NBX IN CCNEX',NBX
           END IF
           CALL PRJDRV(ISIDE,ISTATPRJ,ISYMPR,NBX,NCCVAR,
     *                             WRK(KRES),WRK(KWRK2))
         ENDIF
C
C        ------------------------------------------
C        Test for convergence of the N-R solutions:
C        ------------------------------------------
         RMXNORM = 0.0D0
         NFIN  = 0
         IF (IPRLE .GT. 1) WRITE(LUPRI,*) ' '
         DO 2000 JR = 1,NBX
            JROOTJ = ISIMC - 1 + JR
            QNORM = DNORM2(NCCVAR,WRK(KRES+(JR-1)*NCCVAR),1)
     *             /MAX(WRK(KTST+JR-1),1.0D-14)
            RMXNORM = MAX(RMXNORM,QNORM)
C
            IF (IPRLE .GT. 1 .OR. DEBUG) THEN
               WRITE(LUPRI,'(A,I3,1P,1(A,D12.5))')
     *         ' ROOT:',JROOTJ,' RESIDUAL TOT:',QNORM
            END IF

            IF (QNORM.GT.THRLE) THEN
               NFIN  = NFIN  + 1
               IF ( NFIN.NE.JR) THEN
                  CALL DCOPY(NCCVAR,WRK(KRES+(JR-1)*NCCVAR),
     *                       1,WRK(KRES+(NFIN-1)*NCCVAR),1)
               END IF
            ELSE
               IF (IPRLE .GT. 1) THEN
                  WRITE(LUPRI,'(A,I3,1P,1(A,D12.5))')
     *            ' ROOT:',JROOTJ,' HAS CONVERGED. RESIDUAL TOT:',QNORM
               END IF
            END IF
 2000    CONTINUE
         NOTCON = NOTCON + NFIN
C
         IF (NFIN.EQ.0 .OR. JCONV.LT.0) GO TO 3999
C
C        --------------------------------------------
C        Use generalized conjugate gradient algorithm
C        to form new trial vectors
C        --------------------------------------------
Cholesky
C
C    Cholesky CC2 did not converge using precond, so
C    uncomment next lines and go back to old safe method
C
       IF (CHOINT) THEN
          KDIA   = KWRK1
          KWRK2  = KDIA  + NCCVAR
          LWRK2  = LWRK  - KWRK2
          IF (LWRK2 .LT. 0 ) THEN
             CALL QUIT('Too little work in CCNEX')
          END IF
          CALL CCLR_DIA(WRK(KDIA),ISYMTR,TRIPLET,WRK(KWRK2),LWRK2)
          CALL LEDIAG(NCCVAR,WRK(KDIA),DTEST)
C
          DO JR = 1,NFIN
             IOFF = (JR-1)*NCCVAR
             DO I=1,NCCVAR
                WRK(IOFF+I) = WRK(IOFF+I)*WRK(KDIA-1+I)
             END DO
          END DO
C
       ELSE
CCH
C         WRITE(LUPRI,'(/I4,1X,A)') NFIN,
C     *         'new trial vectors before PRECOND:'
C         WRITE (LUPRI,*) 'NORM: ',DDOT(NCCVAR*NFIN,WRK(KRES),1,
C     *                                   WRK(KRES),1)
C         CALL OUTPUT(WRK(KRES),1,NCCVAR,1,NFIN,NCCVAR,NFIN,1,LUPRI)
CCH
         CALL CC_PRECOND(WRK(KRES),NFIN,NCCVAR,ISYMTR,TRIPLET,
     &                   WRK(KWRK1),LWRK1)
       END IF
C
Cholesky
C
         IF ((IPRLE.GT.105 .OR. LOCDBG).AND. NFIN .GT. 0) THEN
           WRITE(LUPRI,'(/I4,1X,A)') NFIN,
     *         'new trial vectors before CCORT/PRJDRV:'
           CALL OUTPUT(WRK(KRES),1,NCCVAR,1,NFIN,NCCVAR,NFIN,1,LUPRI)
           CALL PROVLP(WRK(KRES),WRK(KRES),NCCVAR,NFIN,WRK(KWRK1),LUPRI)
           CALL FLSHFO(LUPRI)
         END IF

C        ------------------------------
C        Projection from updated trials
C        ------------------------------
         IF (LPROJECT) THEN
           KPRJ   = KWRK2
           ISYMPR = ILSTSYM(LIST,ISTRVE)
           IF (DEBUG) THEN
             WRITE(LUPRI,*)' NFIN IN CCNEX',NFIN
             WRITE(LUPRI,*)
     &            ' ---- CCNEX Before projection on trials ----'
           END IF
           CALL PRJDRV(ISIDE,ISTATPRJ,ISYMPR,NFIN,NCCVAR,
     *                             WRK(KRES),WRK(KPRJ))
         ENDIF
CCN
C        -----------------------------------
C        calculate metric for trial vectors:
C        -----------------------------------
         IF (CCR12) THEN
           DO I = 1, NFIN
             LNOREAD = .TRUE.
             IF (ISIDE .EQ. 1) THEN
               CALL CC_R12METRIC(ISYMTR,BRASCL,KETSCL,
     *                           WRK(KWRK1),LWRK1,DUMMY,DUMMY,
     *                           DUMMY,DUMMY,FS12AM,LUFS12,
     *                           FS2AM,LUFS2,NREDH+I,LNOREAD,
     *                           WRK(KRES+NCCVAR*(I-1)))
             ELSE IF (ISIDE .EQ. -1) THEN
               CALL CC_R12METRIC(ISYMTR,0.5D0*KETSCL,2.0D0*BRASCL,
     *                           WRK(KWRK1),LWRK1,DUMMY,DUMMY,
     *                           DUMMY,DUMMY,FS12AM,LUFS12,
     *                           FS2AM,LUFS2,NREDH+I,LNOREAD,
     *                           WRK(KRES+NCCVAR*(I-1)))
             END IF
           END DO
         END IF
CCN
C        --------------------------------------------------------------
C        Orthogonalize trial vectors and examine for linear dependence:
C        --------------------------------------------------------------
         CALL CCORT(LUFC1,FC1AM,LUFC2,FC2AM,LUFC12,FC12AM,
     *              LUFS2,FS2AM,LUFS12,FS12AM,
     *              TRIPLET,NFIN,NREDH,NCCVAR,
     *              WRK(KRES),ISYMTR,THRLDP,WRK(KWRK1),LWRK1,IPRLE)
C
         IF ((IPRLE.GT.105.OR.LOCDBG)  .AND. NFIN.GT.0) THEN
           WRITE (LUPRI,'(/I4,1X,A)') NFIN,
     *       'new trial vectors after CCORT/PRJDRV:'
           CALL OUTPUT(WRK(KRES),1,NCCVAR,1,NFIN,NCCVAR,NFIN,1,LUPRI)
           CALL PROVLP(WRK(KRES),WRK(KRES),NCCVAR,NFIN,WRK(KWRK1),LUPRI)
         END IF

 3999    CONTINUE
         NUPVEC = NUPVEC + NFIN
 4000 CONTINUE

C we have now:
C NUPVEC : number of new trial vectors
C NOTCON : number of equations not yet converged
C
      IF (NOTCON .EQ. 0 .AND. NUPVEC.EQ.0) THEN 
         ! ALL EQUATIONS HAVE CONVERGED
         IF (IPRLE.GT.10) WRITE(LUPRI,'(A)') ' *** EQUATIONS CONVERGED'
         JCONV = 1
         RETURN
      END IF
C
C     Do NOT calculate new trial vectors
C
      IF (JCONV .LT. 0) THEN
         JCONV = 0
         RETURN
      END IF
C
C     Linear dependence between trial vectors
C
      IF (NUPVEC .EQ. 0 .AND. NOTCON.GT.0) THEN
         WRITE(LUPRI,5010)
         WRITE(LUPRI,*)(NVEC-NOTCON),' trial vectors converged'
         WRITE(LUPRI,*)' LINEAR DEPENDENCE BETWEEN',NOTCON,
     *         ' NON CONVERGED TRIAL VECTORS'
         JCONV = -1
         RETURN
      ENDIF
 5010 FORMAT(/' *** MICROITERATIONS STOPPED DUE TO LINEAR ',
     *        'DEPENDENCE BETWEEN ALL NEW TRIAL VECTORS')
C
C-----------------------------------------------------------
C     Ove Christiansen 2-4-1996:
C     Quick and dirty way of skipping all vectors up to this
C     point and restart from the optimal solution so far.
C-----------------------------------------------------------
C
      IF ( NREDH .GE. MXLRV ) THEN
C        IF ( IPRLE .GT. 20 ) THEN
            WRITE(LUPRI,*) ' NREDH : ',NREDH
            WRITE(LUPRI,*) ' MXLRV : ',MXLRV
C        ENDIF
C
         WRITE(LUPRI,'(/,1X,A)') 
     *      'ATTENTION!!!! '
         WRITE(LUPRI,'(/,1X,A,/,A)') 
     *      'All trial and transformed vectors are skipped '//
     *      'and the iterative procedure ',
     *      ' is continued from the current optimal solution '
         WRITE(LUPRI,'(/1X,A,/)') 
     *      'Notice: You asked for it by setting a low MXLRV.'
C
         IF(ITLE.LE.1) THEN
          WRITE(LUPRI,*)' MXLRV:',MXLRV,
     &           ' ALLOWS ONE ITERATION ,ITLE:',ITLE
          WRITE(LUPRI,*)' NEW INFORMATION NOT CONTAINED IN'
          WRITE(LUPRI,*)' NEW REDUCED SPACE, SPECIFY NEW MXLRV'
          CALL QUIT('CCNEX SPECIFY NEW MXLRV')
         END IF
         CCRSTRS = CCRSTR
         CCRSTR  = .TRUE.
         KIPLAC = 1
         KWRK1  = KIPLAC + MAXRED
         LWRK1  = LWRK   - KWRK1
         CALL CCEQ_STR(FC1AM,LUFC1,FC2AM,LUFC2,FC12AM,LUFC12,
     *                 FS12AM,LUFS12,FS2AM,LUFS2,LPROJECT,ISTATPRJ,
     *                 TRIPLET,ISIDE,ISTRVE,NVEC,NUPVEC,
     *                 NREDH,EIVAL,WRK(KIPLAC),
     *                 WRK(KWRK1),LWRK1,LIST)
         CCRSTR  = CCRSTRS
         IF ( IPRLE .GT. 20 ) THEN
            WRITE(LUPRI,*) ' NREDH : ', NREDH
            WRITE(LUPRI,*) ' NUPVEC: ', NUPVEC
            WRITE(LUPRI,*) ' NVEC:   ', NVEC
            WRITE(LUPRI,*) ' MAXRED: ', MAXRED
            WRITE(LUPRI,*) ' MXLRV : ', MXLRV
         ENDIF
         CALL FLSHFO(LUPRI)
C
      ENDIF
C
      IF (DEBUG)  THEN 
         CALL AROUND(' End of CCNEX ')
         CALL FLSHFO(LUPRI)
      ENDIF
C
C End of CCNEX
C
      RETURN
      END
C  /* Deck ccred */
      SUBROUTINE CCRED(FRHO1,LUFR1,FRHO2,LUFR2,FRHO12,LUFR12,
     *                 FC1AM,LUFC1,FC2AM,LUFC2,FC12AM,LUFC12,
     *                 CCR12,FS12AM,LUFS12,FS2AM,LUFS2,
     *                 LINQCC,TRIPLET,ISIDE,LIST,
     *                 LPROJECT,ISTATPRJ,LREDS,REDS,
     *                 REDH,NREDH,ISTRVE,NUPVEC,NVEC,
     *                 EIVAL,SOLEQ,WRK,LWRK,DEBUG)
C
C
C Input:
C  ISIDE, 1 for right transformation ,-1 for left
C  NREDH,  dimension of the (updated) Hessian/Jacobian matrix in the
C          reduced space (NREDH is not changed in CCRED)
C  ISTRVE  is startnumber on list of equations to be solved.
C  NUPVEC, number of new b-vectors for which REDH is extended
C  NVEC,   number of linear equations or roots to be solved for 
C          in the reduced space
C  
C  LINQCC.EQ.T : solve nvec set of linear equations 
C                EIVAL contains NVEC frequencies 
C                EIVAL( ISTRVE ... (ISTRVE+NVEC-1) )
C
C  LINQCC.EQ.F : solve for NVEC lowest eigenvalues
C  
C  (The reduced PROJECTED HESSIAN matrix is the projection
C   on the basis of b-vectors.
C   REDH(I,J) = B(I) * S(J)   )
C
C Output:
C  REDH,   Jacobian in reduced space (dimension: NREDH)
C  REDS,   overlap in reduced space (dimension: NREDH)
C  SOLEQ,  solutions to the NVEC set of NEWTON equations
C          or eigenvalue equations
C  EIVAL,  eigenvalues or frequencies
C
C Flow:
C 1. if NUPVEC gt 0 update reduced jacobian (and metric) with nupvec 
C    vectors
C 2. determine NVEC solution vectors, returned in SOLEQ for eigenvalue 
C    equation the reduced eigenvalues are returned in EIVAL
C
#include <implicit.h>
#include <priunit.h>
#include <cclr.h>
#include <ccfro.h>
#include <leinf.h>
#include <ccrc1rsp.h>
#include <r12int.h>
#include <ccsdsym.h>
C
      PARAMETER (D1 = 1.0D0 , D0 = 0.0D0, THRZER = 1.0D-99 )
C
      CHARACTER*(*) FRHO1,FRHO2,FC1AM,FC2AM,FRHO12,FC12AM,FS12AM
      CHARACTER*(*) LIST,FS2AM
      LOGICAL   LINQCC, LPROJECT, TRIPLET, DEBUG, CCR12, LREDS, LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      DIMENSION REDH(*),SOLEQ(MAXRED,*),REDS(*)
      DIMENSION EIVAL(*),WRK(*)
C
C **************************************************************
C Section 1:
C extend reduced PROJECTED HESSIAN matrix with NUPVEC new b-vectors
C **************************************************************
C
      IF (IPRLE.GT.5) THEN
         WRITE(LUPRI,*)' CCRED  '
         WRITE(LUPRI,*)' ISIDE  ',ISIDE
         WRITE(LUPRI,*)' LREDS  ',LREDS
         WRITE(LUPRI,*)' NREDH  ',NREDH
         WRITE(LUPRI,*)' ISTRVE ',ISTRVE
         WRITE(LUPRI,*)' NUPVEC ',NUPVEC
         WRITE(LUPRI,*)' NVEC   ',NVEC
         WRITE(LUPRI,*)' LWRK   ',LWRK
      END IF
C
      IF (NUPVEC.GT.0) THEN 
         IF (LREDS) THEN
           MAXVEC = (LWRK-NCCVAR)/(2*NCCVAR)
         ELSE
           MAXVEC = (LWRK-NCCVAR)/NCCVAR
         END IF
         
         IREDH  = NREDH - NUPVEC 
         DO 100 IVEC = 1,NUPVEC,MAXVEC 
            NSIM = MIN(MAXVEC,NUPVEC+1-IVEC)
C
C           work space allocation
C
            KSVEC  = 1
            KBVEC  = KSVEC + NSIM*NCCVAR
            KSBVEC = KBVEC + NCCVAR
C
C           read in sig vectors in batches
C   
            LSVEC = KSVEC
            DO INUM = 1,NSIM
               IVENU = IREDH + IVEC-1 + INUM
               CALL CC_GETVEC(LUFR1,FRHO1,LUFR2,FRHO2,LUFR12,FRHO12,
     *                        TRIPLET,IVENU,WRK(LSVEC))
               LSVEC = LSVEC + NCCVAR
            END DO
            IF (LOCDBG) THEN
               WRITE(LUPRI,'(a,i3,a,i3)')
     *           'A * C for ',IVEC,' - ',IVEC+NSIM-1
               CALL OUTPUT(WRK(KSVEC),1,NCCVAR,1,NSIM,
     *                     NCCVAR,NSIM,1,LUPRI)
            END IF

            IF (LREDS) THEN
              ! read in S * b vectors in batches, assumes that the
              ! conventional CC amplitudes are in the orthogonal basis
              LSVEC = KSBVEC
              DO INUM = 1,NSIM
                 IVENU = IREDH + IVEC-1 + INUM
                 IF (IANR12.EQ.1 .OR. IANR12.EQ.3) THEN
                   CALL CC_GETVEC(LUFC1,FC1AM,LUFC2,FC2AM,LUFS12,FS12AM,
     *                            TRIPLET,IVENU,WRK(LSVEC))
                 ELSE IF (IANR12.EQ.2) THEN
                   CALL CC_GETVEC(LUFC1,FC1AM,LUFS2,FS2AM,LUFS12,FS12AM,
     *                            TRIPLET,IVENU,WRK(LSVEC))
                 ELSE
                   CALL QUIT('This R12 ansatz is not implemented yet')
                 END IF
                 LSVEC = LSVEC + NCCVAR
              END DO
              IF (LOCDBG) THEN
                 WRITE(LUPRI,'(a,i3,a,i3)') 
     *             'S * C for ',IVEC,' - ',IVEC+NSIM-1
                 CALL OUTPUT(WRK(KSBVEC),1,NCCVAR,1,NSIM,
     *                       NCCVAR,NSIM,1,LUPRI)
              END IF
            END IF
C
C           read in b vectors and extend colomns of the Jacobian 
C           (and the metric) in the reduced space
C
            DO I = 1,NREDH
               CALL CC_GETVEC(LUFC1,FC1AM,LUFC2,FC2AM,LUFC12,FC12AM,
     *                        TRIPLET,I,WRK(KBVEC))
               IF (LOCDBG) THEN
                 WRITE(LUPRI,'(a,i3)') 'C for ',I
                 CALL OUTPUT(WRK(KBVEC),1,NCCVAR,1,1,NCCVAR,1,1,LUPRI)
               END IF
               DO JNUM = 1,NSIM
                  J = IREDH + IVEC -1 + JNUM
                  REDH( I + (J-1)*MAXRED ) =
     *            DDOT(NCCVAR,WRK(KBVEC),1,WRK(KSVEC+(JNUM-1)*NCCVAR),1)
                  IF (LREDS) THEN
                    REDS( I + (J-1)*MAXRED ) = DDOT(NCCVAR,
     *                   WRK(KBVEC),1,WRK(KSBVEC+(JNUM-1)*NCCVAR),1)
                  END IF
               END DO
            END DO

 100     CONTINUE
C
C        extend row vectors in reduced spaces matrices if IREDH gt 0
C
         IF (IREDH.GT.0) THEN
            IF (LREDS) THEN
              MAXVEC = (LWRK-2*NCCVAR)/NCCVAR
            ELSE
              MAXVEC = (LWRK-NCCVAR)/NCCVAR
            END IF

            DO 600 IVEC = 1,NUPVEC,MAXVEC 
               NSIM = MIN(MAXVEC,NUPVEC+1-IVEC)
C
               KBVEC  = 1
               KSVEC  = KBVEC + NSIM*NCCVAR
               KSBVEC = KSVEC + NCCVAR
C
               LBVEC = KBVEC
               DO INUM = 1,NSIM
                  IVENU = IREDH + IVEC -1 + INUM
                  CALL CC_GETVEC(LUFC1,FC1AM,LUFC2,FC2AM,LUFC12,FC12AM,
     *                           TRIPLET,IVENU,WRK(LBVEC))
                  LBVEC = LBVEC + NCCVAR
               END DO
               IF (LOCDBG) THEN
                 WRITE(LUPRI,'(a,i3,a,i3)')'C for ',IVEC,'-',IVEC+NSIM-1
                 CALL OUTPUT(WRK(KBVEC),1,NCCVAR,1,NSIM,
     *                       NCCVAR,NSIM,1,LUPRI)
               END IF
C
C              read in s vectors and extend rows of projected matrices
C
               DO J = 1,IREDH
                  CALL CC_GETVEC(LUFR1,FRHO1,LUFR2,FRHO2,LUFR12,FRHO12,
     *                           TRIPLET,J,WRK(KSVEC))
                  IF (LREDS) THEN
                   IF (IANR12.EQ.1 .OR. IANR12.EQ.3) THEN
                     CALL CC_GETVEC(LUFC1,FC1AM,LUFC2,FC2AM,LUFS12,
     *                              FS12AM,TRIPLET,J,WRK(KSBVEC))
                   ELSE IF (IANR12.EQ.2) THEN
                     CALL CC_GETVEC(LUFC1,FC1AM,LUFS2,FS2AM,LUFS12,
     *                              FS12AM,TRIPLET,J,WRK(KSBVEC))
                   ELSE
                     CALL QUIT('This R12 ansatz is not implemented yet')
                   END IF
                  END IF
                  DO INUM = 1,NSIM
                     I = IREDH + IVEC -1 + INUM
                     REDH( I + (J-1)*MAXRED ) = DDOT(NCCVAR,
     *                   WRK(KSVEC),1, WRK(KBVEC+(INUM-1)*NCCVAR),1)
                     IF (LREDS) THEN
                       REDS( I + (J-1)*MAXRED ) = DDOT(NCCVAR,
     *                   WRK(KSBVEC),1,WRK(KBVEC+(INUM-1)*NCCVAR),1)
                     END IF
                  END DO
               END DO

 600        CONTINUE
         END IF
C     
         IF (LINQCC) THEN
           NCAU = 0
           IF (NEWCAU) THEN
             IF (LIST(1:2).EQ.'RC') THEN
               DO IVEC = 1, NVEC
                 ILSTNR = ISTRVE-1 + IVEC
                 NCAU = MAX(NCAU,ILRCAU(ILSTNR))
               END DO
               IF (IPRLE.GE.1) THEN
                WRITE(LUPRI,'(3X,A,I5)') 
     &           'Nb. of Cauchy iterations in reduced space:', NCAU
               END IF
             END IF
           END IF
         END IF
C
      END IF ! (NUPVEC.GT.0) THEN 
C

 1111    CONTINUE

C
      IF (NUPVEC.GT.0) THEN 
         IF (LINQCC) THEN
C
C set up reduced gradient
C
C determine maximum number of simultaneous gradients 
C
            MAXVEC = (LWRK-3*NCCVAR)/NCCVAR 
            DO 1000 IVEC = 1,NVEC,MAXVEC
               NSIM = MIN(MAXVEC,NVEC+1-IVEC)
C
               KGDVEC = 1
               KWRK1  = KGDVEC + NSIM*NCCVAR
               LWRK1  = LWRK - KWRK1
               IF (LWRK1 .LT. 0 ) CALL QUIT('Too little work in CCRED')
C
               LGDVEC = KGDVEC
               DO 1100 INUM = 1,NSIM
C
C-------------------------------
C           Get gradient vector.
C-------------------------------
C
C GETGD bruger GD + 3 exxtra. NCCVAR - skal vaere foerst i memory.
                  ILSTNR = ISTRVE + IVEC -1 + INUM - 1 
                  CALL CC_GETGD(WRK(LGDVEC),ILSTNR,ISIDE,LIST,
     *                          WRK(KWRK1),LWRK1)
C
cs ------------------
C            Project from gradient
cs ------------------
               IF (LPROJECT) THEN
                  ISYPR = ILSTSYM(LIST,ILSTNR)
                  IF (DEBUG) THEN
                     WRITE(LUPRI,*)'CCRED projection from gradient'
                  END IF
                  CALL PRJDRV(ISIDE,ISTATPRJ,ISYPR,1,NCCVAR,
     *                            WRK(LGDVEC),WRK(KWRK1))
               END IF
cs ------------------
C
                  LGDVEC = LGDVEC + NCCVAR
 1100          CONTINUE
C
               KBVEC = KWRK1
               KWRK2 = KBVEC  + NCCVAR
               LWRK2 = LWRK   - KWRK2
               IF (LWRK2 .LT. 0 ) CALL QUIT('Too little work in CCRED')
C
C
C read in b vectors and set up reduced gradient in SOLEQ
C
               DO 1200 I = 1,NREDH
                  CALL CC_GETVEC(LUFC1,FC1AM,LUFC2,FC2AM,LUFC12,FC12AM,
     *                           TRIPLET,I,WRK(KBVEC))
                  DO 1300 JNUM = 1,NSIM
                     J =  IVEC -1 + JNUM
                     SOLEQ(I,J ) =
     *              -DDOT(NCCVAR,WRK(KBVEC),1,
     *                    WRK(KGDVEC+(JNUM-1)*NCCVAR),1)
 1300             CONTINUE
 1200          CONTINUE
               IF (IPRLE.GT.105) THEN
                  WRITE(LUPRI,'(/A)')' ** REDUCED GRADIENT i constr **'
                  CALL OUTPUT(SOLEQ(1,1),1,NREDH,1,NSIM,
     *                        MAXRED,NSIM,1,LUPRI)
               END IF
 1000       CONTINUE
C
         END IF

         IF (IPRLE.GE.15 .OR. LOCDBG) THEN
            IF (LREDS) THEN 
              WRITE(LUPRI,'(/A)')' ** REDUCED OVERLAP MATRIX **'
              CALL OUTPUT(REDS,1,NREDH,1,NREDH,MAXRED,MAXRED,1,LUPRI)
            END IF
            WRITE(LUPRI,'(/A)')  ' ** REDUCED HESSIAN MATRIX **'
            CALL OUTPUT(REDH,1,NREDH,1,NREDH,MAXRED,MAXRED,1,LUPRI)
         END IF
C
      END IF ! (NUPVEC.GT.0) THEN 
C
C **************************************************************
C Section 2
C   find solution vectors in reduced space
C **************************************************************
C
      IF ( LINQCC ) THEN
C
C
C **************************************************************
C Section 2: FIND SOLUTIONS TO THE NVEC SET OF NEWTON-RAPHSON EQ.
C **************************************************************
C
C     SOLVE NVEC SET OF N-R EQ.
C     Copy REDH to SLEVEC 
C
C     SOLEQ contain reduced gradients
C
         KMA = MAXRED*MAXRED
         DO 2100 IVEC = 1,NVEC
            CALL DCOPY(KMA,REDH, 1,WRK,1)
            IF (CCR12 .AND. DABS(EIVAL(IVEC)).GT.1.0D-15) THEN
              CALL DAXPY(NREDH*NREDH,-EIVAL(IVEC),REDS,1,WRK,1)
            ELSE
              DO 2200 I = 1,NREDH
                 IADR = (I-1)*MAXRED + I
                 WRK(IADR) = WRK(IADR) - EIVAL(IVEC)
 2200         CONTINUE
            END IF
            IF (IPRLE.GE.15 .OR. LOCDBG) THEN
               IF (LREDS) THEN 
                 WRITE(LUPRI,'(/A)')' ** REDUCED OVERLAP MATRIX **'
                 CALL OUTPUT(REDS,1,NREDH,1,NREDH,MAXRED,MAXRED,1,LUPRI)
               END IF
               WRITE(LUPRI,'(/A)')' ** REDUCED HESSIAN - w1 MATRIX **'
               WRITE(LUPRI,'(A,F10.6)') '  FREQUENCY: ',EIVAL(IVEC)
               CALL OUTPUT(WRK,1,NREDH,1,NREDH,MAXRED,MAXRED,1,LUPRI)
               WRITE(LUPRI,'(/A)')' ** REDUCED GRADIENT  **'
               CALL OUTPUT(SOLEQ(1,IVEC),1,NREDH,1,1,MAXRED,1,1,LUPRI)
            END IF
            NGD = 1
            CALL DGESOL(NGD,NREDH,MAXRED,MAXRED,WRK,
     *                  SOLEQ(1,IVEC),WRK(1+KMA),INFO)
C           CALL DGESOL(NSIM,N,LDA,LDB,A,B,KPVT,INFO)
            IF (INFO.NE.0) WRITE(LUPRI,7000)IVEC,INFO
 7000 FORMAT(/' ***CCRED*** SOLUTION NOT OBTAINED TO LINEAR EQUATIONS'
     *      ,I5,/' CHECK IF HESSIAN MATRIX IS SINGULAR,',
     *        ' INFO from DGESOL/DSPSOL =',I7)
            IF (IPRLE.GE.15 .OR. LOCDBG) THEN
               WRITE(LUPRI,'(/A)')' ** REDUCED SOLUTION VECTOR  **'
               CALL OUTPUT(SOLEQ(1,IVEC),1,NREDH,1,1,MAXRED,1,1,LUPRI)
            END IF
 2100    CONTINUE
         IF (IPRLE.GT.15 .OR. LOCDBG) THEN
            WRITE(LUPRI,*)'  * SOLUTION TO',NVEC,' NEWTON EQUATIONS *'
            CALL OUTPUT(SOLEQ(1,1),1,NREDH,1,NVEC,MAXRED,NVEC,
     *                     1,LUPRI)
         END IF
C
      END IF
C
      IF (LINQCC.AND.NEWCAU.AND.NCAU.GT.1) THEN
        NCAU = NCAU - 1
        ISYM = ILSTSYM(LIST,ISTRVE)
        NVARPT = NCCVAR + 2*NALLAI(ISYM)
        MAXVEC = (LWRK-NVARPT)/NCCVAR
        DO IVEC = 1, NVEC, MAXVEC
          NBVEC = MIN(MAXVEC,NVEC+1-IVEC)
          KWRK1 = 1 + NCCVAR*NBVEC
          LWRK1 = LWRK - KWRK1
          IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Out of memory in CCRED')
          END IF 
            CALL CCCONV(LUFC1,FC1AM,LUFC2,FC2AM,LUFC12,FC12AM,TRIPLET,
     *                  CCR12,NREDH,IVEC,NBVEC,SOLEQ,WRK(1),WRK(KWRK1))
          DO I = 1, NBVEC
            IDX = IVEC-1+I
            CALL CC_SAVE(WRK(1+(I-1)*NCCVAR),ISTRVE-1+IDX,
     *                   LIST,WRK(KWRK1),LWRK1)
          END DO
        END DO
        GOTO 1111
      END IF
C
      IF (.NOT. LINQCC ) THEN
C
C
C **************************************************************
C Section 2: Solve reduced generalized eigenvalue problem in subspace
C **************************************************************
C
         MATZ = 1
         KWI  = 1
         KAMAT= KWI  + MAXRED
         KIV1 = KAMAT+ MAXRED*MAXRED
         KFV1 = KIV1 + MAXRED
         KEND = KFV1 + MAXRED
         IF (LREDS) THEN
           KSMAT  = KEND
           KEND   = KSMAT  + MAXRED*MAXRED
           KDENOM = KFV1
         END IF

         LEND = LWRK - KEND
         IF (KEND .GT. LWRK) CALL ERRWRK('LERED.RG',KEND,LWRK)

         IF (LREDS) THEN
           ! solve generalized eigenvalue problem for a real general
           ! jacobian and a nontrivial metric 
           CALL DCOPY(MAXRED*MAXRED,REDH,1,WRK(KAMAT),1)
           CALL DCOPY(MAXRED*MAXRED,REDS,1,WRK(KSMAT),1)
c
           CALL RGG(MAXRED,NREDH,WRK(KAMAT),WRK(KSMAT),EIVAL,
     *              WRK(KWI),WRK(KDENOM),MATZ,SOLEQ,IERR)
           DO I = 1, NREDH
             IF (ABS(WRK(KDENOM-1+I)).GT.THRZER) THEN
               EIVAL(I) = EIVAL(I)/WRK(KDENOM-1+I)
               WRK(KWI-1+i) = WRK(KWI-1+I)/WRK(KDENOM-1+I)
             ELSE
               EIVAL(I) = D1/THRZER
               WRK(KWI-1+i) = WRK(KWI-1+I)/THRZER
             END IF
           END DO
         ELSE
           ! solve eigenvalue problem for a real general jacobian
           CALL DCOPY(MAXRED*MAXRED,REDH,1,WRK(KAMAT),1)
           CALL RG(MAXRED,NREDH,WRK(KAMAT),EIVAL,WRK(KWI),MATZ,SOLEQ,
     *             WRK(KIV1),WRK(KFV1),IERR)
         END IF

         IF (IPRLE .GE. 70 .OR. IERR .NE. 0) THEN
            WRITE (LUPRI,'(/A)') ' REDUCED EIGENVALUES real part:'
            WRITE (LUPRI,'(A)') ' before sort of eigenvalues'
            CALL OUTPUT(EIVAL,1,NREDH,1,1,NREDH,MAXRED,1,LUPRI)
            WRITE (LUPRI,'(/A)')
     *           ' REDUCED EIGENVALUES imaginary part:'
            WRITE (LUPRI,'(A)') ' before sort of eigenvalues'
            CALL OUTPUT(WRK(KWI),1,NREDH,1,1,NREDH,MAXRED,1,LUPRI)
            WRITE (LUPRI,'(/A)') ' REDUCED EIGENVECTORS :'
            WRITE (LUPRI,'(A)') ' before sort of eigenvalues'
            CALL OUTPUT(SOLEQ,1,NREDH,1,NREDH,MAXRED,MAXRED,1,LUPRI)
         END IF
         IF ( IERR.NE.0 ) THEN
            WRITE(LUPRI,'(/A,I5)')
     *      ' EIGENVALUE PROBLEM NOT CONVERGED IN RG, IERR =',IERR
               CALL QUIT(' CCRED: EIGENVALUE EQUATION NOT CONVERGED ')
         END IF
C
         CALL RGORD(MAXRED,NREDH,EIVAL,WRK(KWI),SOLEQ,.FALSE.)
C
         ICPLX = 0
         DO I=1,NVEC
            IF (WRK(I) .NE. D0) THEN
               ICPLX = ICPLX + 1
               WRITE(LUPRI,'(I10,1P,2D15.8,A/)') I,EIVAL(I),WRK(I),
     *            ' *** CCRED  WARNING **** COMPLEX VALUE.'
            END IF
         END DO
C
         IF (IPRLE .GE.1 .OR. ICPLX .GT. 0) THEN
            WRITE (LUPRI,'(/A)') ' REDUCED EIGENVALUES real part:'
            CALL OUTPUT(EIVAL,1,NVEC,1,1,MAXRED,1,1,LUPRI)
         ENDIF
C
         IF (IPRLE .GE.11 .OR. ICPLX .GT. 0) THEN
            WRITE (LUPRI,'(/A)') ' REDUCED EIGENVALUES imaginary part:'
            CALL OUTPUT(WRK(KWI),1,NVEC,1,1,MAXRED,1,1,LUPRI)
         END IF
C
         IF (IPRLE .GE. 15) THEN
            WRITE (LUPRI,'(/A)') ' REDUCED EIGENVECTORS :'
            CALL OUTPUT(SOLEQ,1,NREDH,1,NVEC,MAXRED,MAXRED,1,LUPRI)
         END IF
C
         IF (ICPLX .GT. 0)
     *      WRITE(LUPRI,*) '**WARNING CCRED: COMPLEX EIGENVALUES.'
C
      END IF
C
      IF (IPRLE.GE.15 ) THEN
         WRITE(LUPRI,'(/A)')' REDUCED HESSIAN MATRIX:'
         CALL OUTPUT(REDH,1,NREDH,1,NREDH,MAXRED,MAXRED,1,LUPRI)
         IF (LREDS) THEN
           WRITE(LUPRI,'(/A)')' REDUCED METRIC MATRIX:'
           CALL OUTPUT(REDS,1,NREDH,1,NREDH,MAXRED,MAXRED,1,LUPRI)
         END IF 
      END IF
C
C *** End of subroutine CCRED
C
      RETURN
      END
C  /* Deck rgord */
      SUBROUTINE RGORD(NM,N,WR,WI,Z,LSORT)
C
C  15-Aug-1989 Hans Joergen Aa. Jensen
C  modified by Ove Christiansen 20-dec 1994
C  to normalize correct!
C
C  Reorder and normalize eigenpairs from RG.
C
C  LSORT = .FALSE. --> ascending order
C  LSORT = .TRUE.  --> descending order
#include <implicit.h>
      LOGICAL LSORT
      DIMENSION WR(N), WI(N), Z(NM,N)
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0 )
      DO 100 I = 1,N-1
         ILOW = I
         XLOW = WR(I)
         DO 50 J = I+1,N
          IF (.NOT.LSORT) THEN
            IF (WR(J) .LT. XLOW) THEN
               ILOW = J
               XLOW = WR(J)
            END IF
          ELSE
            IF (WR(J) .GT. XLOW) THEN
               ILOW = J
               XLOW = WR(J)
            END IF
          END IF
   50    CONTINUE
         IF (ILOW .NE. I) THEN
            WR(ILOW) = WR(I)
            WR(I)    = XLOW
            XLOW     = WI(ILOW)
            WI(ILOW) = WI(I)
            WI(I)    = XLOW
            CALL DSWAP(N,Z(1,I),1,Z(1,ILOW),1)
         END IF
  100 CONTINUE
C
C     Normalize eigenvectors
C
      DO 200 I = 1,N
         ZNRM = DDOT(N,Z(1,I),1,Z(1,I),1)
         ZNRM = D1 / SQRT(ZNRM)
         IMAX = IDAMAX(N,Z(1,I),1)
         IF (Z(IMAX,I) .LT. D0) ZNRM = -ZNRM
         CALL DSCAL(N,ZNRM,Z(1,I),1)
  200 CONTINUE
      RETURN
      END
C  /* Deck ccort */
      SUBROUTINE CCORT(LUFC1,C1FIL,LUFC2,C2FIL,LUFC12,C12FIL,
     *                 LUFS2,S2FIL,LUFS12,S12FIL,
     *                 TRIPLET,NBADD,NBPREV,NCCVAR,
     *                 BVECS,ISYMTR,THRLDP,WRK,LWRK,IPRLE)
C 
C Modified version of ORTBVLE of Hans Jorgen Aa Jensen
C specialized for CC 
C
C Purpose:
C  Orthogonalize new b-vectors against all previous b-vectors
C  and among themselves, and renormalize.
C  (Orthogonalization is performed twice if round-off is large,
C   if larger than THRRND).
C
C Input:
C  LUFC1  - file with singles parts of trial vectors with name C1FIL
C  LUFC2  - file with doubles parts of trial vectors with name C2FIL
C  LUFC12 - file with R12 doubles parts of trial vec. with name C12FIL
C  BVECS  - non-orthogonal  b-vectors
C  NBADD  - number of BVECS
C  NBPREV - number of trial vectors on LUB
C  THRLDP - threshold for linear dependence
C  
C
C Output:
C  NBADD trial vectors added to previous
C  NBPREV updated to contain NBADD new trial vectors
C
C Scratch:
C  require WRK of length NCCVAR
C
#include <implicit.h>
#include <priunit.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <r12int.h>
C
      CHARACTER C1FIL*(*),C2FIL*(*),C12FIL*(*),S2FIL*(*),S12FIL*(*)
      DIMENSION BVECS(NCCVAR,*),WRK(NCCVAR)
      LOGICAL TRIPLET, LOCDBG
      PARAMETER (THRRND=1.D-4, THRTT=1.D-6, THRZER=1.0D-35)
      PARAMETER (THRTST=1.0D-12)
      PARAMETER (D1 = 1.0D0)
      PARAMETER (LOCDBG = .FALSE.)
      INTEGER SVEC, RVEC, LVEC
C
C
C
      NSAVE = NBADD
C
      IF ((IPRLE.GT.120) .OR. LOCDBG) THEN
         WRITE (LUPRI,'(//A/)') ' >>> Debug output from CCORT <<<'
      END IF
C
C     work space allocation
C
      KEND1  = 1
      IF (CCR12) THEN
        KSVECS = KEND1
        KEND1  = KSVECS+ NCCVAR*NBADD
      END IF
      LWRK1  = LWRK - KEND1
      IF (LWRK1.LT.0) THEN
        WRITE(LUPRI,*) 'LWRK, KEND1: ',LWRK, KEND1
        CALL QUIT('Insufficient memory in CCORT')
      END IF
C
C     Normalize input BVECS and remove
C     vectors with norm .le. THRZER
C
      NNUM = 0
      DO I = 1, NBADD

         IF (CCR12) THEN
           CALL DCOPY(NCCVAR,BVECS(1,I),1,WRK(KSVECS+NCCVAR*(I-1)),1)
           CALL CC_RVEC(LUFS12,S12FIL,NTR12AM(ISYMTR),NTR12AM(ISYMTR),
     *                  NBPREV+I,WRK(KSVECS+NCCVAR*(I-1)+NT1AM(ISYMTR)+
     *                  NT2AM(ISYMTR)))
           IF (IANR12.EQ.2) THEN
             CALL CC_RVEC(LUFS2,S2FIL,NT2AM(ISYMTR),NT2AM(ISYMTR),
     *                    NBPREV+I,WRK(KSVECS+NCCVAR*(I-1)+
     *                    NT1AM(ISYMTR)))
           END IF
         END IF   

         IF (CCR12) THEN
C          Check norm of r12 part:
           TST = DDOT(NTR12AM(ISYMTR),BVECS((1+NT1AM(ISYMTR)+
     *                NT2AM(ISYMTR)),I),1,WRK(KSVECS+NCCVAR*(I-1)+
     *                NT1AM(ISYMTR)+NT2AM(ISYMTR)),1)           
           IF (TST .LT. -THRTST) THEN
             WRITE(LUPRI,*) '(CCORT) Norm^2 (R12-part) = ',TST
             WRITE(LUPRI,'(A,I4,A)') 'WARNING: R12 part of input '//
     *                      'vector no. ',I,' in CCORT has norm < 0'
C            WRITE(LUPRI,*) '>>>>>>> Check your '//
C    *                      'auxiliary basis and/or approximations!'
C            CALL QUIT('R12 input vector in CCORT has norm < 0')
           END IF
           TT = DDOT(NCCVAR,BVECS(1,I),1,WRK(KSVECS+NCCVAR*(I-1)),1)
         ELSE
           TT = DDOT(NCCVAR,BVECS(1,I),1,BVECS(1,I),1)
         END IF
         IF ((IPRLE.GT.20) .OR. LOCDBG) THEN
            WRITE (LUPRI,*) ' I, initial norm  of BVECS(1,I) :',I,TT
         END IF
         IF (TT.LT.THRZER) THEN
            IF (IPRLE.GT.2) WRITE (LUPRI,3240) I,TT
 3240       FORMAT(/' (CCORT) new  b-vector no.',I3,
     *       /' is removed at input because of too small norm;'
     *       /' norm of input vector is ',1P,D12.5)
         ELSE
            NNUM = NNUM + 1
            TT = SQRT(TT)
            IF (TT.LT.THRTT) THEN
               CALL DSCAL(NCCVAR,(D1/THRTT),BVECS(1,I),1)
               IF (CCR12) CALL 
     &           DSCAL(NCCVAR,(D1/THRTT),WRK(KSVECS+NCCVAR*(I-1)),1)
               IF (CCR12) THEN
                 TT=DDOT(NCCVAR,BVECS(1,I),1,WRK(KSVECS+NCCVAR*(I-1)),1)
               ELSE
                 TT = DDOT(NCCVAR,BVECS(1,I),1,BVECS(1,I),1)
               END IF
               TT = SQRT(TT)
            END IF
            CALL DSCAL (NCCVAR,(D1/TT),BVECS(1,I),1)
            IF (CCR12) CALL 
     &        DSCAL(NCCVAR,(D1/TT),WRK(KSVECS+NCCVAR*(I-1)),1)
            IF ( I.GT.NNUM) THEN 
               CALL DCOPY(NCCVAR,BVECS(1,I),1,BVECS(1,NNUM),1)
               IF (CCR12) CALL DCOPY(NCCVAR,WRK(KSVECS+NCCVAR*(I-1)),1,
     *                               WRK(KSVECS+NCCVAR*(NNUM-1)),1)
            END IF
         END IF
      END DO
C
      IF (NBADD.NE.NNUM) THEN
         WRITE(LUPRI,*)NNUM-NBADD,' vectors with norm le:',
     *                 THRZER,' REMOVED'
         NBADD = NNUM
      END IF
      IF (NBADD .EQ. 0) THEN
         WRITE (LUPRI,'(/1X,2A,I5)')
     *      'CCORT,**WARNING**  number vectors added is zero',
     *      ' after check for vector norm , NBADD :',NBADD
      END IF
C
C work space allocation
C
      KBPRE  = KEND1
      KEND1  = KBPRE + NCCVAR
      IF (CCR12) THEN
        KSPRE  = KEND1
        KEND1  = KSPRE + NCCVAR
      END IF
      IROUND = 0
      ITURN  = 0

 1500 CONTINUE
      ITURN  = ITURN + 1
C
C     Orthogonalize new b-vectors against previous b-vectors
C
      DO K = 1, NBPREV
         IF (NBADD .GT. 0) THEN
            CALL CC_GETVEC(LUFC1,C1FIL,LUFC2,C2FIL,LUFC12,C12FIL,
     *                     TRIPLET,K,WRK(KBPRE))
            IF (CCR12) THEN
C             put S * b into WRK(KSPRE)
              CALL DCOPY(NCCVAR,WRK(KBPRE),1,WRK(KSPRE),1) 
              IF (IANR12.EQ.1 .OR. IANR12.EQ.3) THEN 
                CALL CC_GETVEC(LUFC1,C1FIL,LUFC2,C2FIL,LUFS12,S12FIL,
     *                      TRIPLET,K,WRK(KSPRE))
              ELSE IF (IANR12.EQ.2) THEN
                CALL CC_GETVEC(LUFC1,C1FIL,LUFS2,S2FIL,LUFS12,S12FIL,
     *                      TRIPLET,K,WRK(KSPRE))
              ELSE 
                CALL QUIT('Unknown IANR12 in CCORT')
              END IF
            END IF
            DO J = 1, NBADD
               IF (CCR12) THEN
                TT=-DDOT(NCCVAR,WRK(KBPRE),1,WRK(KSVECS+NCCVAR*(J-1)),1)
               ELSE
                 TT = -DDOT(NCCVAR,WRK(KBPRE),1,BVECS(1,J),1)
               END IF
               IF ((IPRLE.GT.20) .OR. LOCDBG) THEN
                  WRITE (LUPRI,*)
     *            ' K (prev.vec.), J (new type a), overlap:',
     *            K,J,TT
               END IF
               CALL DAXPY(NCCVAR,TT,WRK(KBPRE),1,BVECS(1,J),1)
               IF (CCR12) THEN
                 CALL DAXPY(NCCVAR,TT,WRK(KSPRE),1,WRK(KSVECS+
     *                      NCCVAR*(J-1)),1)
               END IF
            END DO
         ELSE
            WRITE (LUPRI,'(/1X,A,2I5)')
     *        'CCORT, number vectors added is zero , NBADD :',NBADD
         END IF
      END DO
CCN
C     IF (LOCDBG) THEN
C       IF (CCR12) THEN
C         KOLP  = KEND1
C         KEND2 = KOLP + NBADD*NBADD
C         IF (KEND2.LT.0) 
C    *      CALL QUIT('Insufficient memory in CCORT')
C         DO I = 1, NBADD
C           DO J = 1, NBADD
C             WRK(KOLP-1+NBADD*(J-1)+I) = DDOT(NTR12AM(ISYMTR),
C    *                   BVECS((1+NT1AM(ISYMTR)+NT2AM(ISYMTR)),I),1,
C    *                   WRK(KSVECS+NCCVAR*(J-1)+
C    *                   NT1AM(ISYMTR)+NT2AM(ISYMTR)),1)  
C           END DO
C         END DO
C         WRITE(LUPRI,*)' R12 PART OF OVERLAP:'
C         CALL OUTPUT(WRK(KOLP),1,NBADD,1,NBADD,NBADD,NBADD,1,LUPRI)
C         WRITE(LUPRI,*)' R12 PART OF B vectors:'
C         CALL OUTPUT(BVECS,1+NT1AM(ISYMTR)+NT2AM(ISYMTR),NCCVAR,
C    *                      1,NBADD,NCCVAR,NBADD,1,LUPRI)
C         WRITE(LUPRI,*)' R12 PART OF S x B vectors:'
C         CALL OUTPUT(WRK(KSVECS),1+NT1AM(ISYMTR)+NT2AM(ISYMTR),NCCVAR,
C    *                      1,NBADD,NCCVAR,NBADD,1,LUPRI)
C       END IF
C     END IF
CCN
C
C     Normalize and orthogonalize new vectors against each other
C     and remove vectors that are linear dependent 
C
      NNUM = 0
      DO I = 1, NBADD
         DO J = 1, NNUM
            IF (CCR12) THEN
              T1 =  DDOT(NCCVAR,BVECS(1,J),1,WRK(KSVECS+NCCVAR*(J-1)),1)
              TT = -DDOT(NCCVAR,BVECS(1,J),1,WRK(KSVECS+NCCVAR*(I-1)),1)
            ELSE
              T1 =  DDOT(NCCVAR,BVECS(1,J),1,BVECS(1,J),1)
              TT = -DDOT(NCCVAR,BVECS(1,J),1,BVECS(1,I),1)
            END IF
            IF ((IPRLE.GT.20) .OR. LOCDBG) THEN
               WRITE (LUPRI,*)
     *         ' I (new type a), J (new type a), norm bJ, overlap:',
     *         I,J,T1,TT
            END IF
            TT = TT / T1
            CALL DAXPY(NCCVAR,TT,BVECS(1,J),1,BVECS(1,I),1)
            IF (CCR12) THEN
              CALL DAXPY(NCCVAR,TT,WRK(KSVECS+NCCVAR*(J-1)),1,
     *                             WRK(KSVECS+NCCVAR*(I-1)),1)
            END IF
         END DO
         IF (CCR12) THEN
C          Check norm of r12 part:
           TST = DDOT(NTR12AM(ISYMTR),BVECS((1+NT1AM(ISYMTR)+
     *                NT2AM(ISYMTR)),I),1,WRK(KSVECS+NCCVAR*(I-1)+
     *                NT1AM(ISYMTR)+NT2AM(ISYMTR)),1)           
           IF (TST .LT. -THRTST) THEN
             WRITE(LUPRI,*) '(CCORT) Norm^2 (R12-part) = ',TST
             WRITE(LUPRI,'(A,I4,A)') 'WARNING: R12 part of result '//
     *                      'vector no. ',I,' in CCORT has norm < 0'
C            WRITE(LUPRI,*) '>>>>>>> Check your '//
C    *                      'auxiliary basis and/or approximations!'
C            CALL QUIT('R12 result vector in CCORT has norm < 0')
           END IF
           TT = DDOT(NCCVAR,BVECS(1,I),1,WRK(KSVECS+NCCVAR*(I-1)),1)
         ELSE
           TT = DDOT(NCCVAR,BVECS(1,I),1,BVECS(1,I),1)
         END IF
         IF ((IPRLE.GT.20) .OR. LOCDBG) THEN
            WRITE (LUPRI,*)
     *      ' I (new type a), norm before renormalization:',
     *      I,SQRT(TT)
         END IF
         IF (TT .LE. THRLDP) THEN
            IF (IPRLE.GT.2) WRITE (LUPRI,3250) I,TT
 3250        FORMAT(/' (CCORT) new  b-vector no.',I3,
     *       /' is removed because of linear dependence;'
     *       /' norm of vector after Gram-Schmidt''s orthogonalization '
     *        ,1P,D12.5)
         ELSE
            NNUM = NNUM + 1
            IF (TT .LT. THRRND) IROUND = IROUND+1
            TT = SQRT(TT)
            IF (TT.LT.THRTT) THEN
               CALL DSCAL(NCCVAR,(D1/THRTT),BVECS(1,I),1)
               IF (CCR12) CALL 
     &           DSCAL(NCCVAR,(D1/THRTT),WRK(KSVECS+NCCVAR*(I-1)),1)
               IF (CCR12) THEN
                 TT=DDOT(NCCVAR,BVECS(1,I),1,WRK(KSVECS+NCCVAR*(I-1)),1)
               ELSE
                 TT = DDOT(NCCVAR,BVECS(1,I),1,BVECS(1,I),1)
               END IF
               TT = SQRT(TT)
            END IF
            CALL DSCAL(NCCVAR,(D1/TT),BVECS(1,I),1)
            IF (CCR12) CALL 
     &        DSCAL(NCCVAR,(D1/TT),WRK(KSVECS+NCCVAR*(I-1)),1)
            IF ( I.NE.NNUM) THEN
               CALL DCOPY(NCCVAR,BVECS(1,I),1,BVECS(1,NNUM),1)
               IF (CCR12) CALL DCOPY(NCCVAR,WRK(KSVECS+NCCVAR*(I-1)),1,
     *                    WRK(KSVECS+NCCVAR*(NNUM-1)),1)
            END IF
         END IF
      END DO
C
      IF (NBADD.NE.NNUM) THEN
         IF (IPRLE.GT.2) THEN
           WRITE(LUPRI,*)' After Gram Schmidt orthonormalization'
           WRITE(LUPRI,*) NNUM-NBADD,' vectors with norm le:'
     *                  ,THRZER,' REMOVED'
         END IF
         NBADD = NNUM
      END IF
      IF (IROUND.GT.0  .AND. ITURN.EQ.1 ) GO TO 1500
C
C Save new vector together with old ones on LUS
C
      NLAST = NBPREV + NBADD 
      DO I = 1, NBADD
         IVENU = NBPREV + I
         CALL CC_PUTVEC(LUFC1,C1FIL,LUFC2,C2FIL,LUFC12,C12FIL,
     *                  TRIPLET,IVENU,BVECS(1,I))
         IF (CCR12) THEN
           CALL CC_WVEC(LUFS12,S12FIL,NTR12AM(ISYMTR),NTR12AM(ISYMTR),
     *                  IVENU,WRK(KSVECS+NT1AM(ISYMTR)+NT2AM(ISYMTR)+
     *                  NCCVAR*(I-1)))
           IF (IANR12.EQ.2) THEN
             CALL CC_WVEC(LUFS2,S2FIL,NT2AM(ISYMTR),NT2AM(ISYMTR),
     *                    IVENU,WRK(KSVECS+NT1AM(ISYMTR)+NCCVAR*(I-1)))
           END IF
         END IF
      END DO
C
C     update NBPREV
C
      NBPREV = NBPREV + NBADD
C
C
      IF ( NSAVE .NE. NBADD ) THEN
         IF (IPRLE.GT.2) WRITE(LUPRI,4310)NSAVE,NBADD
 4310    FORMAT(/' (CCORT)  trial vectors reduced from',I3,' to',I3)
      END IF
C
C test if trial vectors on disk are orthonormal
C
      IF (( DEBUG .AND. IPRLE .GT. 20) .OR. LOCDBG) THEN
         KVEC  = 1
         KOLP  = KVEC + NBPREV*NCCVAR
         KWRK1 = KOLP + NBPREV*NBPREV
         IF (CCR12) THEN
           SVEC  = KWRK1
           KWRK1 = SVEC + NBPREV*NCCVAR
         END IF 
         LWRK1 = LWRK - KWRK1
         IF ( LWRK1 .LE.0) CALL QUIT('CCORT TOO LITTLE WRK')
         LVEC = KVEC
         RVEC = SVEC
         DO I = 1, NBPREV
            CALL CC_GETVEC(LUFC1,C1FIL,LUFC2,C2FIL,LUFC12,C12FIL,
     *                     TRIPLET,I,WRK(LVEC))
            IF (CCR12) THEN
              IF (IANR12.EQ.1 .OR. IANR12.EQ.3) THEN
                CALL CC_GETVEC(LUFC1,C1FIL,LUFC2,C2FIL,LUFS12,S12FIL,
     *                         TRIPLET,I,WRK(RVEC))
              ELSE IF (IANR12.EQ.2) THEN
                CALL CC_GETVEC(LUFC1,C1FIL,LUFS2,S2FIL,LUFS12,S12FIL,
     *                         TRIPLET,I,WRK(RVEC))
              ELSE
                CALL QUIT('Unknown IANR12 in CCORT')
              END IF 
            END IF
            LVEC = LVEC + NCCVAR
            RVEC = RVEC + NCCVAR
         END DO
         WRITE(LUPRI,*)'*** CCORT TEST ****'
         WRITE(LUPRI,*)' OVERLAP BETWEEN',NBPREV,' TRIAL VECTORS'
         IF (CCR12) THEN
          CALL PROVLP(WRK(KVEC),WRK(SVEC),NCCVAR,NBPREV,WRK(KOLP),LUPRI)
         ELSE
          CALL PROVLP(WRK(KVEC),WRK(KVEC),NCCVAR,NBPREV,WRK(KOLP),LUPRI)
         END IF
         IF (IPRLE.GT.50) THEN
           CALL OUTPUT(WRK(KVEC),1,NCCVAR,1,NBPREV,NCCVAR,NBPREV,
     *                     1,LUPRI)
         END IF
         CALL FLSHFO(LUPRI)
      END IF
C
C *** End of subroutine CCORT
C
      IF (IPRLE .GT. 20) THEN
        WRITE (LUPRI,'(/A//)')' >>> End of debug output from CCORT <<<'
      END IF

      RETURN
      END
C  /* Deck provlp */
      SUBROUTINE PROVLP(VECS,SVECS,LVEC,NVECS,OVLP,LUPRI)
C
C 22-Aug-1989 Hans Joergen Aa. Jensen
C
C Print overlap matrix between the NVECS input vectors.
C
#include <implicit.h>
      DIMENSION VECS(LVEC,NVECS), OVLP(*), SVECS(LVEC,NVECS)
C
      IJ = 0
      DO 200 I = 1,NVECS
         DO 100 J = 1,NVECS
            IJ = IJ + 1
            OVLP(IJ) = DDOT(LVEC,VECS(1,I),1,SVECS(1,J),1)
  100    CONTINUE
  200 CONTINUE
      WRITE (LUPRI,'(/A)') ' Overlap matrix :'
      CALL OUTPUT(OVLP,1,NVECS,1,NVECS,NVECS,NVECS,1,LUPRI)
      RETURN
      END
C  /* Deck lediag */
      SUBROUTINE LEDIAG(NDIAG,DIAG,DTEST)
C
C  5-Nov-1988 Hans Joergen Aa. Jensen
C
C Calculate inverse diagonal,
C for generalized conjugate gradient algorithm.
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION DIAG(NDIAG)
      PARAMETER ( D1 = 1.0D0 )
C
      DO 100 I = 1,NDIAG
         IF (ABS(DIAG(I)).LE.DTEST) THEN
            DIAG(I) = D1 / SIGN(DTEST,DIAG(I))
         ELSE
            DIAG(I) = D1 / DIAG(I)
         END IF
  100 CONTINUE
      RETURN
      END
C  /* Deck cc_save */
      SUBROUTINE CC_SAVE(VEC,ILSTNR,LIST,WORK,LWORK)
C
C     Ove Christiansen 14-11-1996
C     Write vec to file based on list and list number.
C
#include <implicit.h>
#include <priunit.h>
#include <cclr.h>
#include <ccorb.h>
#include <leinf.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
Cholesky
#include <maxorb.h>
#include <ccdeco.h>
Cholesky
C
      CHARACTER*(*) LIST
      CHARACTER*14  MODEL, LABR12
C
      DIMENSION VEC(*), WORK(LWORK)
C
      INTEGER ILSTSYM
C
      MODEL = '              '
      IF (CCR12 .AND. (.NOT.(CCS.OR.CIS))) THEN
        ILEN   = 4
        LABR12 = '-R12'
      ELSE
        ILEN = 1
        LABR12 = ' '
      END IF 

      MODEL = 'CCSD'//LABR12(1:ILEN)
      IF (CCS)  MODEL = 'CCS'
      IF (CC2)  MODEL = 'CC2'//LABR12(1:ILEN)
      IF (CC3)  MODEL = 'CC3'//LABR12(1:ILEN)
      IF (CC1A) MODEL = 'CCSDT-1A'//LABR12(1:ILEN)
      IF (CC1B) MODEL = 'CCSDT-1B'//LABR12(1:ILEN)
C
      ISYM = ILSTSYM(LIST,ILSTNR)
C
      KT1   = 1
      KT2   = 1 + NT1AM(ISYM)
C
      IOPT   = 3
Cholesky
      IF (CHOINT .AND. CC2) IOPT = 1
Cholesky
      CALL CC_WRRSP(LIST,ILSTNR,ISYM,IOPT,MODEL,DUMMY,VEC(KT1),
     *              VEC(KT2),WORK,LWORK)
      IF (CCR12 .AND. (.NOT.(CCS.OR.CIS))) THEN
        IOPT = 32
        KT12 = KT2 + NT2AM(ISYM)
        CALL CC_WRRSP(LIST,ILSTNR,ISYM,IOPT,MODEL,DUMMY,DUMMY,
     *                VEC(KT12),WORK,LWORK)
      END IF
C
      END 
*=====================================================================*
      SUBROUTINE PRJDRV(ISIDE,ISTATPRJ,ISYMPRJ,NVECS,NDIM,VECT,WORK)
*
*     driver subroutine for projection.
*     ISIDE = 1 : solve right response equations (A-w)R = X
*     ISIDE =-1 : solve left response equations  L(A+w) = X
*     ISTATPRJ  : index on list of the excited state used for proj.
*     ISYMPRJ   : symmetry of BOTH the vector and the eigevector.
*     Please note: The vector we project from and the eigenvector
*                  we project away MUST have the same symmetry
*
*     Sonia Coriani sept 1997. Revised march 2000
*=====================================================================*
#include <implicit.h>
#include <priunit.h>
#include <ccsdinp.h>
#include <leinf.h>
#include <ccexci.h>
#include <ccsdsym.h>
#include <cbirea.h>
*
      DIMENSION VECT(*), WORK(*)
      CHARACTER*10 MODFIL
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
*
      IF (LOCDBG) THEN
         WRITE(LUPRI,*) 'PRJDRV: IPRLE  =    ',IPRLE
         WRITE(LUPRI,*) 'PRJDRV: ISIDE  =    ',ISIDE
         WRITE(LUPRI,*) 'PRJDRV: ISTATPRJ =  ',ISTATPRJ
         WRITE(LUPRI,*) 'PRJDRV: ISYMPRJ  =  ',ISYMPRJ
         WRITE(LUPRI,*) 'PRJDRV: NCCVAR =    ',NDIM
         WRITE(LUPRI,*) 'PRJDRV: NVECS  =    ',NVECS
      ENDIF
*
      IF (LMULBS) CALL QUIT('PRJDRV not yet adapted for CC-R12')

      IF (CCS) THEN
        IOPT = 1
      ELSE
        IOPT = 3
      END IF
*
* set dimension of work space
*
      NCCVAR1 = NT1AM(ISYMPRJ)
      NCCVAR2 = NT2AM(ISYMPRJ)
*
* allocate WORK space for RE and LE
*
      KREVT1   = 1                
      KREVT2   = KREVT1 + NCCVAR1 
      KLEVT1   = KREVT2 + NCCVAR2 
      KLEVT2   = KLEVT1 + NCCVAR1 
*
*  puntatore (offset) 
*
      IOFFSTP = ISTATPRJ 
*
* get right eigenvector RE
*
      CALL CC_RDRSP('RE',IOFFSTP,ISYMPRJ,IOPT,MODFIL,
     *                                   WORK(KREVT1),WORK(KREVT2))
*
* get left eigenvector LE
*
      CALL CC_RDRSP('LE',IOFFSTP,ISYMPRJ,IOPT,MODFIL,
     *                                   WORK(KLEVT1),WORK(KLEVT2))
*     
* check dimension matching
*
      IF (CCS) THEN 
         NDIMV = NCCVAR1
      ELSE
         NDIMV = NCCVAR1+NCCVAR2
      ENDIF
      IF (NDIMV.NE.NDIM) THEN
         WRITE(LUPRI,*)' Warning PRJDRV: vects dimension inconsistent'
         CALL QUIT('PRJDRV: NDIMV.ne.NDIM')
      ENDIF
c
      DO IVC = 1, NVECS
        IF (ISIDE.EQ.1) THEN
*
* use LE for constant, RE for subtraction
*
          KOFFV = (IVC-1)*NDIMV+1   !offset vector to be projected
          CALL PRJECT(NDIMV,VECT(KOFFV),WORK(KLEVT1),WORK(KREVT1))
*
          CHECKL = DDOT(NDIMV,VECT(KOFFV),1,WORK(KLEVT1),1)
          CHECKR = DDOT(NDIMV,VECT(KOFFV),1,WORK(KREVT1),1)
          IF (LOCDBG) THEN
            WRITE(LUPRI,*) 'After PRJECT from right: <LE|VP> = ', CHECKL
            WRITE(LUPRI,*) 'After PRJECT from right: <RET|VP> = ',CHECKR
          END IF
*
        ELSE 
*
* use RE for the constant, LE for subtraction
*
          KOFFV = (IVC-1)*NDIMV+1   !offset vector to be projected
          CALL PRJECT(NDIMV,VECT(KOFFV),WORK(KREVT1),WORK(KLEVT1))
C
          CHECKR = DDOT(NDIMV,VECT(KOFFV),1,WORK(KREVT1),1)
          CHECKL = DDOT(NDIMV,VECT(KOFFV),1,WORK(KLEVT1),1)
ccs          CALL PRJECT(NDIMV,VECT,WORK(KREVT1),WORK(KLEVT1))
ccs          CHECKR = DDOT(NDIMV,VECT,1,WORK(KREVT1),1)
ccs          CHECKL = DDOT(NDIMV,VECT,1,WORK(KLEVT1),1)
          IF (LOCDBG) THEN
            WRITE(LUPRI,*) 'After PRJECT from left: <VP|RE> = ', CHECKR
            WRITE(LUPRI,*) 'After PRJECT from left: <VP|LET> = ',CHECKL
          END IF
        ENDIF
      END DO

      RETURN
      END
*================================================================*
      SUBROUTINE PRJECT(NDIMVPR,VECT,CVECT,SVECT)
*================================================================*
* Sonia Coriani sept. 1999
* subroutine for projecting out eigenstates
* NDIMVPR = dimension of the vectors
* VECT  = input vector to project/output projected vector
* CVECT = left/right vector  for the constant
* SVECT = right/left vector to be subtracted
*================================================================*
#include <implicit.h>
#include <priunit.h>
#include <ccsdinp.h>
*
      DIMENSION VECT(NDIMVPR),CVECT(NDIMVPR),SVECT(NDIMVPR)
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
*
* compute the constant for projection
* const = <left|vect> or const = <vect|right>
*
      CONST = - DDOT(NDIMVPR,VECT,1,CVECT,1)
      IF (LOCDBG) THEN
        WRITE(LUPRI,*) '=== We are inside PRJECT === '
        WRITE(LUPRI,*) 'constant for projection = ', CONST
      END IF
*  
* get projected vector
* |vectp> = |vect> - |right>*<left|vect>
* or 
* <vectp| = <vect| - <vect|right>*<left|
*
      CALL DAXPY(NDIMVPR,CONST,SVECT,1,VECT,1) 
*
      RETURN
      END
*======================================================================*
      SUBROUTINE CC_PRECOND(VECTORS,NVEC,NCCVAR,ISYMTR,TRIPLET,
     &                      WORK,LWORK)
      IMPLICIT NONE
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <dummy.h>
#include <r12int.h>

      LOGICAL TRIPLET, LCCEQ
      INTEGER NVEC, NCCVAR, ISYMTR, LWORK

#if defined (SYS_CRAY)
      REAL VECTORS(NCCVAR,NVEC), WORK(*), DTEST
#else
      DOUBLE PRECISION VECTORS(NCCVAR,NVEC), WORK(*), DTEST
#endif
      PARAMETER ( DTEST = 1.0D-04 )

      INTEGER KDIA, KEND1, KSCR, LWRK1, NCCVAR0, JVEC

      KDIA  = 1
      KEND1 = KDIA  + NCCVAR

      IF (CCR12) THEN
        KSCR = KEND1
        KEND1= KSCR + NTR12AM(ISYMTR) 
      END IF
  
      LWRK1 = LWORK - KEND1
      IF (LWRK1.LT.0) THEN
        CALL QUIT('Too little work space in CC_PRECOND')
      END IF

      CALL CCLR_DIA(WORK(KDIA),ISYMTR,TRIPLET,WORK(KEND1),LWRK1)
      CALL LEDIAG(NCCVAR,WORK(KDIA),DTEST)

      IF (CCS) THEN
        NCCVAR0 = NT1AM(ISYMTR)
      ELSE
        NCCVAR0 = NT1AM(ISYMTR) + NT2AM(ISYMTR)
        IF (TRIPLET) NCCVAR0 = NCCVAR0 + NT2AMA(ISYMTR)
      END IF

      DO JVEC = 1,NVEC
         ! for conventional part divide by orbital energy differences
         DO I=1,NCCVAR0
            VECTORS(I,JVEC) = VECTORS(I,JVEC)*WORK(KDIA-1+I)
         END DO

         IF (CCR12) THEN
            ! for R12 part divide by B matrix
            LCCEQ = .FALSE.
            CALL CC_R12NXTAM(VECTORS(NCCVAR0+1,JVEC),ISYMTR,
     &                       WORK(KSCR),LCCEQ,
     &                       DUMMY,WORK(KEND1),LWRK1)
            CALL DSCAL(NTR12AM(ISYMTR),-1.0D0,VECTORS(NCCVAR0+1,JVEC),1)
         END IF 
      END DO

      RETURN
      END
*======================================================================*
