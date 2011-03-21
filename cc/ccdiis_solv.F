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
C  /* Deck ccdiis_sol */
      SUBROUTINE CCDIIS_SOL(FRHO1,LUFR1,FRHO2,LUFR2,
     *                    FC1AM,LUFC1,FC2AM,LUFC2,LIST,ISTART,
     *                    TRIPLET,ISIDE,NVEC,NUPVEC,
     *                    NREDH,REDH,REDS,EIVAL,SOLEQ,
     *                    AMAT,ITRAN,CONVERGED,RNORM,WRK,LWRK)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C 
C input
C
C   FRHO1,FRHO2,FC1AM,FC2AM is file names for files where
C   transformed vectors (rho1,rho2) and trial vectors (c1am,c2am) are
C   stored.
C   LUFR1,LUFR2,LUFC1,LUFC2 is the corresponding unit numbers.
C
C   TRIPLET.EQ.T solve triplet equations
C          .EQ.F solve singlet equations
C   
C   ISIDE.EQ.1  SOLVE EQUATIONS FROM RIGHT
C   ISIDE.EQ.-1 SOLVE EQUATIONS FROM LEFT
C
C   NVEC   = NUMBER OF EQUATIONS THAT MUST BE SOLVED
C   NUPVEC = NUMBER OF CURRENT LINEAR INDEPENDENT NEW VECTORS.
C
C   EIVAL CONTAINS NVEC FREQUENCIES 
C
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include <implicit.h>
#include <priunit.h>
#include <ccsdinp.h>
#include <cclr.h>
#include <leinf.h>
      PARAMETER (D0 = 0.0D0 )
C
      CHARACTER*(*)  FRHO1,FRHO2,FC1AM,FC2AM
      CHARACTER*8 FS12AM, FRHO12, FC12AM, FS2AM
      CHARACTER*(*) LIST
      CHARACTER*3 APROXR12
      DIMENSION REDH(*), REDS(*), SOLEQ(*)
      DIMENSION EIVAL(*), AMAT((MXDIIS+1)*(MXDIIS+1),*)
      DIMENSION WRK(LWRK), RNORM(*)
      LOGICAL TRIPLET, CONVERGED(*)
      INTEGER ITRAN(*)
C
      CALL QENTER('CCDIIS_SOL')
      IF (IPRLE .GT. 3 ) THEN
         WRITE(LUPRI,*) 'CCDIIS_SOL: IPRLE   =    ',IPRLE
         WRITE(LUPRI,*) 'CCDIIS_SOL: TRIPLET =    ',TRIPLET
         WRITE(LUPRI,*) 'CCDIIS_SOL: NREDH   =    ',NREDH
         WRITE(LUPRI,*) 'CCDIIS_SOL: MAXRED  =    ',MAXRED
         WRITE(LUPRI,*) 'CCDIIS_SOL: MAXLE   =    ',MAXLE
         WRITE(LUPRI,*) 'CCDIIS_SOL: ISIDE   =    ',ISIDE
         WRITE(LUPRI,*) 'CCDIIS_SOL: NVEC    =    ',NVEC
         WRITE(LUPRI,*) 'CCDIIS_SOL: NUPVEC  =    ',NUPVEC
         WRITE(LUPRI,*) 'CCDIIS_SOL: NCCVAR  =    ',NCCVAR
      ENDIF
C
      IF (CCR12) CALL QUIT('No R12 yet in CCDIIS_SOL.')
C
      IF (NUPVEC.NE.NVEC) 
     *  CALL QUIT('Not enough start vectors in CCDIIS_SOLV.')
      IF (NREDH.NE.NVEC) 
     *  CALL QUIT('Wrong dimension of reduced space in CCDIIS_SOLV.')
C
      CALL GETTIM(CSTR,WSTR)
C
      TIMMIC = SECOND()
      TIMLIN = D0
      TIMRED = D0
      TIMNEX = D0
      ITLE   = 0
C
C     print banner for solver:
C
      IF (IPRLE.GT.0) THEN
         WRITE (LUPRI,'(///A/)') 
     &        ' >>>>> COUPLED CLUSTER DIIS SOLVER <<<<<'
         WRITE (LUPRI,'(3X,A)')
     &        ' Iter  #Vectors  time (min)   residual'
         WRITE (LUPRI,'(3X,A)')
     &        ' --------------------------------------'
      END IF 
C
C---------------------------------------------------------
C     mark all vectors as not converged, 
C     and start loop over iterations
C---------------------------------------------------------
C
      DO I = 1, NVEC
         CONVERGED(I) = .FALSE.
      ENDDO
C
C     Loop over iterations. 
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
C-----------------------------------------------------
C        Transform vectors :
C        Nredh  : Dim. of reduced space
C        Nupvec : # new trial vectors in this ite.
C-----------------------------------------------------
C
          NUPVEC = 0
          DO I = 1, NVEC
            IF (.NOT.CONVERGED(I)) THEN
               NUPVEC = NUPVEC + 1
               ITRAN(NUPVEC) = I
            END IF
          ENDDO

         IF (IPRLE .GT. 5) THEN
            WRITE (LUPRI,*)
     *      ' >>> Call CC_TRDRV with TRIPLET flag set to:',
     *      TRIPLET
         ENDIF                      
C
         CALL CC_TRDRV(ECURR,FRHO1,LUFR1,FRHO2,LUFR2,FRHO12,LUFR12,
     *                 FC1AM,LUFC1,FC2AM,LUFC2,FC12AM,LUFC12,
     *                 TRIPLET,.TRUE.,ITRAN,EIVAL,FS12AM,IDUMMY,
     *                 FS2AM,IDUMMY,
     *                 ISIDE,0,NUPVEC,WRK,LWRK,APROXR12)
C
         CALL FLSHFO(LUPRI)
         TIM    = SECOND() - TIM
         TIMLIN =  TIMLIN + TIM
C
C--------------------------------
C        CALL CCRED(), INCREASE DIMENSION OF REDUCED SPACE,
C        AND FIND SOLUTIONS
C--------------------------------
C
         TIMRED = TIMRED - SECOND()
C
         CALL CCRED2(FRHO1,LUFR1,FRHO2,LUFR2,
     *               FC1AM,LUFC1,FC2AM,LUFC2,
     *               TRIPLET,REDH,REDS,NREDH,
     *               EIVAL,SOLEQ,WRK,LWRK,DEBUG)
C
         CALL FLSHFO(LUPRI)
         TIMRED = TIMRED + SECOND()
C
C
C        CONSTRUCT SOLUTION VECTORS IN FULL SPACE, SAVE THEM ON FILE,
C        CHECK CONVERGENCE AND GET NEW TRIAL VECTORS FOR NON-CONVERGED
C        EIGENVECTORS
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
         CALL CCNEX_DIIS(FRHO1,LUFR1,FRHO2,LUFR2,
     *                   FC1AM,LUFC1,FC2AM,LUFC2,LIST,ISTART,
     *                   ITLE,TRIPLET,NREDH,
     *                   JCONV,EIVAL,SOLEQ,AMAT,
     *                   CONVERGED,WRK,LWRK,RNORM,RMXNORM)
C
C-----------------------------------------------------
C        Print timing & convergence statistic
C-----------------------------------------------------
C
         IF (IPRLE.GT.0) THEN
            WRITE (LUPRI,'(2X,I5,3X,I5,1X,F12.2,1X,E12.2)')
     *         ITLE,NUPVEC,TIM/60.0D0,RMXNORM
         END IF
         CALL FLSHFO(LUPRI)
C
C-----------------------------------------------------
C        Print eigenvalues and residuals:
C-----------------------------------------------------
C
         IF (IPRLE .GT. 1 ) THEN
           WRITE(LUPRI,'(/5X,A)')' eigenvalues    residual   converged'
           WRITE(LUPRI,'(5X,A)')'.....................................'
           DO I = 1, NVEC
              WRITE(LUPRI,'(5X,F12.9,E13.4,5X,L2)')
     *          EIVAL(I), RNORM(I), CONVERGED(I)
           END DO
           WRITE(LUPRI,'(5X,A)')'.....................................'
         ENDIF
C
         CALL FLSHFO(LUPRI)
         TIMNEX = TIMNEX + SECOND()
C
         IF (IPRLE .GE. 2) THEN
            TIMIT = SECOND() - TIMIT
            WRITE (LUPRI,'(/A,F12.2,/)')
     *         ' >>> TIME USED IN THIS ITERATION',TIMIT
         END IF
C
         IF (JCONV.LT.0) THEN
C           ( LINEAR DEPENDENCE BETWEEN ALL NEW TRIAL VECTORS )
            WRITE (LUPRI,'(/A/A)')
     *         ' ***  CCDIIS_SOL - ITERATIONS STOPPED  ',
     *         '     LINEAR DEPENDENCE BETWEEN ALL NEW TRIAL VECTORS'
         ELSE IF (JCONV.GT.0) THEN
C           (CONVERGED)
            IF (IPRLE .GT. 10) THEN
               WRITE(LUPRI,'(/,(A,I4,A,A,I4,A))')
     *         ' ITERATIONS CONVERGED FOR',NVEC,
     *         ' SOLUTION VECTORS',
     *         ' IN',ITLE,' ITERATIONS.'
            END IF
         ELSE IF (ITLE .LT. MAXLE) THEN
            CALL FLSHFO(LUPRI)
            GO TO 100
         ELSE
            WRITE(LUPRI,'(/A,I5,A)')
     *         ' *** CCDIIS_SOL - MAXIMUM NUMBER OF ITERATIONS',
     *         ITLE,' REACHED.'
               CALL QUIT(' *** CCDIIS_SOL-MAX. ITERATIONS REACHED')
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
      IF (IPRLE .GT. 0) THEN
         WRITE (LUPRI,'(//A/A/A,1P,D10.2,A/)')
     *      '  Final eigenvalues :',
     *      ' =====================',
     *      ' (convergence threshold:',THRLE,')'
         DO 920 I = 1,NVEC
            WRITE (LUPRI,'(I5,F15.10)') I,EIVAL(I)
  920    CONTINUE
      END IF
      IF (IPRLE .GE. 10) THEN
         WRITE (LUPRI,'(//A)')
     *      ' FINAL EIGENVECTORS IN REDUCED BASIS FOR :'
         IOFF = 0
         DO 930 I = 1,NVEC
            WRITE (LUPRI,'(/A,I4,A,1P,D15.6/A)')
     *         ' - EIGENVALUE NO.',I, ' - EIGENVALUE',EIVAL(I),
     *         ' - EIGENVECTOR :'
            WRITE (LUPRI,'(5F15.8)') (SOLEQ(IOFF+J),J=1,NREDH)
            IOFF = IOFF + MAXRED
  930    CONTINUE
      END IF
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('CCDIIS_SOL')
      RETURN
      END
C  /* Deck ccred2 */
      SUBROUTINE CCRED2(FRHO1,LUFR1,FRHO2,LUFR2,
     *                  FC1AM,LUFC1,FC2AM,LUFC2,
     *                  TRIPLET,REDH,REDS,NREDH,
     *                  EIVAL,SOLEQ,WRK,LWRK,DEBUG)
C
C
C Input:
C  NREDH,  dimension of new reduced PROJECTED HESSIAN matrix, NREDH is
C          not changed in CCRED
C  
C  EIVAL CONTAINS NVEC FREQUENCIES 
C
C Output:
C  REDH,   the new, extended reduced PROJECTED HESSIAN matrix
C          (dimension: NREDH)
C  SOLEQ,  solutions to the NVEC set of NEWTON equations
C          or eigenvalue equations
C  EIVAL,  eigenvalues or frequencies
C
C Flow:
C 1. set up reduced jacobian and overlap matrix
C 2. determine NREDH solution vectors, returned in SOLEQ,
C    and eigenvalues, returned in EIVAL.
C
#include <implicit.h>
#include <priunit.h>
#include <cclr.h>
#include <ccfro.h>
#include <leinf.h>
#include <ccrc1rsp.h>
C
      PARAMETER (D1 = 1.0D0 , D0 = 0.0D0 , THRZER = 1.0D-99)
C
      CHARACTER*(*)  FRHO1,FRHO2,FC1AM,FC2AM
      CHARACTER*(8)  FC12AM, FRHO12
      LOGICAL   TRIPLET, DEBUG
      DIMENSION REDH(*), REDS(*), SOLEQ(MAXRED,*)
      DIMENSION EIVAL(*),WRK(*)
C
      CALL QENTER('CCRED2')
C
C **************************************************************
C Section 1: set up reduced JACOBIAN and OVERLAP matrix 
C **************************************************************
C
      IF (IPRLE.GT.5) THEN
         WRITE(LUPRI,*)' CCRED  '
         WRITE(LUPRI,*)' NREDH  ',NREDH
         WRITE(LUPRI,*)' LWRK   ',LWRK
      END IF
C
      MAXVEC =  (LWRK-NCCVAR)/NCCVAR
      IF (MAXVEC.LT.1) CALL QUIT('Insufficient memory in CCRED2.')
C
      DO IVEC = 1,NREDH,MAXVEC 
         NSIM = MIN(MAXVEC,NREDH+1-IVEC)

         KBVEC = 1
         KSVEC = KBVEC + NSIM*NCCVAR
C
C        read in b vectors 
C
         LBVEC = KBVEC
         DO INUM = 1,NSIM
            I = IVEC - 1 + INUM
            CALL CC_GETVEC(LUFC1,FC1AM,LUFC2,FC2AM,LUC12AM,FC12AM,
     *                     TRIPLET,I,WRK(KBVEC+(INUM-1)*NCCVAR))
         END DO
C
C        read in s vectors and extend projected jacobian matrix
C
         DO J = 1,NREDH
            CALL CC_GETVEC(LUFR1,FRHO1,LUFR2,FRHO2,LUFR12,FRHO12,
     *                           TRIPLET,J,WRK(KSVEC))
            DO INUM = 1,NSIM
               I = IVEC - 1 + INUM
               REDH( I + (J-1)*MAXRED ) =
     *           DDOT(NCCVAR,WRK(KSVEC),1,WRK(KBVEC+(INUM-1)*NCCVAR),1)
            END DO
         END DO
C
C        read in b vectors and extend projected overlap matrix
C
         DO J = 1,NREDH
            CALL CC_GETVEC(LUFC1,FC1AM,LUFC2,FC2AM,LUC12,FC12AM,
     *                     TRIPLET,J,WRK(KSVEC))
            DO INUM = 1,NSIM
               I = IVEC -1 + INUM
               REDS( I + (J-1)*MAXRED ) =
     *           DDOT(NCCVAR,WRK(KSVEC),1,WRK(KBVEC+(INUM-1)*NCCVAR),1)
            END DO
         END DO

      END DO
C
C **************************************************************
C Section 2: Solve reduced generalized eigenvalue problem in subspace
C **************************************************************
C
C     Use EISPACK routine for real general matrices in
C     generalized eigenvalue problem
C
      MATZ = 1
      KWI    = 1
      KDENOM = KWI    + MAXRED
      KAMAT  = KDENOM + MAXRED
      KSMAT  = KAMAT  + MAXRED*MAXRED
      KEND   = KSMAT  + MAXRED*MAXRED
      LEND    = LWRK - KEND
      IF (KEND .GT. LWRK) CALL ERRWRK('LERED.RG',KEND,LWRK)
      CALL DCOPY(MAXRED*MAXRED,REDH,1,WRK(KAMAT),1)
      CALL DCOPY(MAXRED*MAXRED,REDS,1,WRK(KSMAT),1)
      CALL RGG(MAXRED,NREDH,WRK(KAMAT),WRK(KSMAT),EIVAL,
     *         WRK(KWI),WRK(KDENOM),MATZ,SOLEQ,IERR)
      DO I = 1, NREDH
        IF (ABS(WRK(KDENOM-1+I)).GT.THRZER) THEN
          EIVAL(I)  = EIVAL(I)/WRK(KDENOM-1+I)
          WRK(KWI-1+i)  = WRK(KWI-1+I)/WRK(KDENOM-1+I)
        ELSE
          EIVAL(I)  = 1.0D0/THRZER
          WRK(KWI-1+i)  = WRK(KWI-1+I)/THRZER
        END IF
      END DO
      IF (IPRLE .GE. 70 .OR. IERR .NE. 0) THEN
         WRITE (LUPRI,'(/A)') ' REDUCED EIGENVALUES real part:'
         WRITE (LUPRI,'(A)') ' before sort of eigenvalues'
         CALL OUTPUT(EIVAL,1,NREDH,1,1,NREDH,MAXRED,1,LUPRI)
         WRITE (LUPRI,'(/A)')
     *        ' REDUCED EIGENVALUES imaginary part:'
         WRITE (LUPRI,'(A)') ' before sort of eigenvalues'
         CALL OUTPUT(WRK(KWI),1,NREDH,1,1,NREDH,MAXRED,1,LUPRI)
         WRITE (LUPRI,'(/A)') ' REDUCED EIGENVECTORS :'
         WRITE (LUPRI,'(A)') ' before sort of eigenvalues'
         CALL OUTPUT(SOLEQ,1,NREDH,1,NREDH,MAXRED,MAXRED,1,LUPRI)
      END IF
      IF ( IERR.NE.0 ) THEN
         WRITE(LUPRI,'(/A,I5)')
     *   ' EIGENVALUE PROBLEM NOT CONVERGED IN RG, IERR =',IERR
            CALL QUIT(' CCRED: EIGENVALUE EQUATION NOT CONVERGED ')
      END IF
C
      CALL RGORD(MAXRED,NREDH,EIVAL,WRK,SOLEQ,.FALSE.)
C
      ICPLX = 0
      DO I=1,NREDH
         IF (WRK(I) .NE. D0) THEN
            ICPLX = ICPLX + 1
            WRITE(LUPRI,'(I10,1P,2D15.8,A/)') I,EIVAL(I),WRK(I),
     *         ' *** CCRED  WARNING **** COMPLEX VALUE.'
         END IF
      END DO
C
      IF (IPRLE .GE.11) THEN
         WRITE (LUPRI,'(/A)') ' REDUCED EIGENVALUES real part:'
         CALL OUTPUT(EIVAL,1,NREDH,1,1,NREDH,MAXRED,1,LUPRI)
         WRITE (LUPRI,'(/A)') ' REDUCED EIGENVALUES imaginary part:'
         CALL OUTPUT(WRK(KWI),1,NREDH,1,1,NREDH,MAXRED,1,LUPRI)
      END IF
C
      IF (IPRLE.GE.15 ) THEN
         WRITE(LUPRI,'(/A)')' *** REDUCED HESSIAN MATRIX ***'
         CALL OUTPUT(REDH,1,NREDH,1,NREDH,MAXRED,MAXRED,1,LUPRI)
      END IF
C
C *** End of subroutine CCRED2
C
      CALL QEXIT('CCRED2')
C
      RETURN
      END
C  /* Deck ccnex2 */
      SUBROUTINE CCNEX_DIIS(FRHO1,LUFR1,FRHO2,LUFR2,
     *                      FC1AM,LUFC1,FC2AM,LUFC2,LIST,ISTART,
     *                      ITLE,TRIPLET,NREDH,
     *                      JCONV,EIVAL,SOLEQ,AMAT,
     *                      CONVERGED,WRK,LWRK,RNORM,RMXNORM)
C
C PURPOSE: 1) Construct residual (A)*X(I) - EIVAL(I)*X(I) 
C             for solution X(I) of reduced  equations
C          2) Test for convergence of solutions
C             convergence criterium:
C             ||((A)*X(I) - EIVAL(I)*X(I)|| / ||X|| .LE. THRLE
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
#include <ccorb.h>
#include <ccfro.h>
#include <leinf.h>
#include <inftap.h>
      PARAMETER (  THRLDP = 1.D-20 )
      PARAMETER ( DTEST = 1.0D-04 )
      PARAMETER ( DM1=-1.0D0,D1 =1.0D0, D0=0.0D0 )
C
      CHARACTER*(*)  FRHO1,FRHO2,FC1AM,FC2AM
      CHARACTER*(*) LIST
      CHARACTER*(8) FC12AM
      CHARACTER*(16) FILDIISS, FILDIIST
      LOGICAL   CCRSTRS,TRIPLET
      DIMENSION EIVAL(*), SOLEQ(MAXRED,*),WRK(*),RNORM(*)
      DIMENSION AMAT((MXDIIS+1)*(MXDIIS+1),NREDH)
      LOGICAL CONVERGED(NREDH)
C
C
C Space for CCNEX_DIIS:
C
C MAXNEX: Maximum number of simultaneous vectors in CCNEX_DIIS
C
      IF (DEBUG)  THEN 
         CALL AROUND(' Start of CCNEX_DIIS ')
      ENDIF
      MAXNEX = (LWRK-3*NCCVAR)/(2*NCCVAR)
      NOTCON = 0
      RMXNORM = 0.0D0
C
C     read orbital energies from file (needed for preconditioning)
C
      KFOCKD= 1
      KWRK1 = KFOCKD + NORBTS
      LWRK1 = LWRK   - KWRK1
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WRK(KFOCKD+I-1), I=1,NORBTS)
      CALL GPCLOSE(LUSIFC,'KEEP')
      IF (FROIMP.OR.FROEXP) 
     &  CALL CCSD_DELFRO(WRK(KFOCKD),WRK(KWRK1),LWRK1) 
      CALL FOCK_REORDER(WRK(KFOCKD),WRK(KWRK1),LWRK1) 
C
      DO ISIMC = 1,NREDH,MAXNEX
         NBX   = MIN(MAXNEX,(NREDH+1-ISIMC))
C
C        CONSTRUCT RESIDUAL IN WRK(KRES)
C        AND SOLUTION VECTORS IN WRK(KSOL)
C        RESIDUAL: R = (A-EIVAL(I))*X(I)
C
         KRES  = KFOCKD + NORBTS
         KSOL  = KRES + NBX*NCCVAR 
         KWRK1 = KSOL + NBX*NCCVAR 
         LWRK1 = LWRK - KWRK1
         IF (LWRK1 .LT. NCCVAR ) THEN 
            CALL QUIT('Too little work in CCNEX_DIIS xx')
         END IF
C
         CALL DZERO(WRK(KRES),NBX*NCCVAR)
C
         CALL CCCONV(LUFC1,FC1AM,LUFC2,FC2AM,LFC12,FC12AM,
     &               TRIPLET,.FALSE.,NREDH,
     *               ISIMC,NBX,SOLEQ,WRK(KSOL),WRK(KWRK1))
C
C        Set norm of solution vector to 1 and save on file
C
         DO IBVEC = 1, NBX
           XNRM = DNORM2(NCCVAR,WRK(KSOL+(IBVEC-1)*NCCVAR),1)
           CALL DSCAL(NCCVAR,1.0D0/XNRM,WRK(KSOL+(IBVEC-1)*NCCVAR),1)

           IVEC = ISIMC - 1 + IBVEC   
           IF (.NOT.CONVERGED(IVEC)) THEN
             CALL CC_SAVE(WRK(KSOL+(IBVEC-1)*NCCVAR),ISTART-1+IVEC,
     *                    LIST,WRK(KWRK1),LWRK1)
           END IF
         END DO
C
         IF (IPRLE.GT.105) THEN
            WRITE (LUPRI,'(/A)')
     *      ' CCNEX: solution vectors'
            CALL PROVLP(WRK(KSOL),WRK(KSOL),NCCVAR,NBX,WRK(KWRK1),LUPRI)
            CALL OUTPUT(WRK(KSOL),1,NCCVAR,1,NBX,NCCVAR,NBX,
     *                  1,LUPRI)
         END IF
C
C        -------------------------------
C        add -EIVAL(I)*X(I) to residual:
C        -------------------------------
C 
         DO IBVEC = 1,NBX
            IVEC = ISIMC - 1 + IBVEC   
            CALL DAXPY(NCCVAR,-EIVAL(IVEC),WRK(KSOL+(IBVEC-1)*NCCVAR),
     *                 1,WRK(KRES+(IBVEC-1)*NCCVAR),1)
         END DO

         IF (IPRLE.GT.110) THEN
           WRITE(LUPRI,*)' RESIDUAL AFTER eival CONTRIBUTION'
               CALL OUTPUT(WRK(KRES),1,NCCVAR,1,NBX,NCCVAR,NBX,
     *                     1,LUPRI)
         END IF
C
C------------------------------------------------
C Add  (A)*X(I) where X(I) is the solution to the
C I'th set of Newton-Raphson equations
C------------------------------------------------
C
         DO 900 K = 1,NREDH
            CALL CC_GETVEC(LUFR1,FRHO1,LUFR2,FRHO2,LUFR12,FRHO12,
     &                     TRIPLET,K,WRK(KWRK1))
            DO 700 JR = 1,NBX
               JROOTJ = ISIMC - 1 + JR
               EVAL1  = SOLEQ(K,JROOTJ)
               CALL DAXPY(NCCVAR,EVAL1,WRK(KWRK1),1,
     *                       WRK(KRES+(JR-1)*NCCVAR),1)
 700        CONTINUE
 900     CONTINUE
C
C-------------------------------
C Residual is now in WRK(KRES) 
C-------------------------------
C
         IF (IPRLE.GT.45) THEN
            WRITE (LUPRI,'(/A)')
     *         ' CCNEX_DIIS: residual vectors '
               CALL PROVLP(WRK(KRES),WRK(KRES),NCCVAR,NBX,WRK(KWRK1),
     *                     LUPRI)
               CALL OUTPUT(WRK(KRES),1,NCCVAR,1,NBX,NCCVAR,NBX,
     *                     1,LUPRI)
         END IF
C
C------------------------------------------
C        Test for convergence
C------------------------------------------
C
         IF (IPRLE .GT. 1) WRITE(LUPRI,*) ' '

         DO JR = 1,NBX
            JROOTJ  = ISIMC - 1 + JR
            QNORM   = DNORM2(NCCVAR,WRK(KRES+(JR-1)*NCCVAR),1)
            RNORM(JROOTJ) = QNORM
            RMXNORM = MAX(RMXNORM,QNORM)
C
            IF (QNORM.GT.THRLE .AND. .NOT.CONVERGED(JROOTJ)) THEN
             NOTCON = NOTCON + 1
             CONVERGED(JROOTJ) = .FALSE.
             IF (.NOT.(JCONV.LT.0)) THEN
               KSOL1 = KSOL  + (JR-1)*NCCVAR
               KSOL2 = KSOL1 + NT1AM(ISYMTR)
               KSOL3 = KSOL2 + NT2AM(ISYMTR)
               KRES1 = KRES  + (JR-1)*NCCVAR
               KRES2 = KRES1 + NT1AM(ISYMTR)
               KRES3 = KRES2 + NT2AM(ISYMTR)
               CALL CCSD_NXTAM(WRK(KSOL1),WRK(KSOL2),WRK(KSOL3),
     *                         WRK(KRES1),WRK(KRES2),WRK(KRES3),
     *                         WRK(KFOCKD),TRIPLET,ISYMTR,EIVAL(JROOTJ))
            
               WRITE(FILDIISS,'(a7,i3)') 'CCDIISS',JROOTJ
               WRITE(FILDIIST,'(a7,i3)') 'CCDIIST',JROOTJ
               DO I = 8, 10
                 IF (FILDIISS(I:I) .EQ. ' ') FILDIISS(I:I) = '_'
                 IF (FILDIIST(I:I) .EQ. ' ') FILDIIST(I:I) = '_'
               ENDDO

               CALL CCEX_DIIS(WRK(KSOL1),WRK(KRES1),FILDIISS,
     *                        FILDIIST,AMAT(1,JROOTJ),
     *                        NCCVAR,MXDIIS,ITLE,
     *                        WRK(KWRK1),LWRK1)


               XNRM = DNORM2(NCCVAR,WRK(KSOL1),1)
               CALL DSCAL(NCCVAR,1.0D0/XNRM,WRK(KSOL1),1)

               CALL CC_PUTVEC(LUFC1,FILC1,LUFC2,FILC2,LUFC12,FILC12,
     &                        TRIPLET,JROOTJ,WRK(KSOL1))
             END IF
            ELSE
             CONVERGED(JROOTJ) = .TRUE.
            END IF
         END DO
      END DO
C
C     NOTCON : Number of vectors not converged.
C 
      IF (NOTCON.EQ.0) THEN 
C        ALL EQUATIONS HAVE CONVERGED
         IF (IPRLE .GT. 10) WRITE(LUPRI,'(A)')' *** EQUATIONS CONVERGED'
         JCONV = 1
         RETURN
      ELSE
         JCONV = 0
      END IF
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
C  /* Deck ccex_diis */
      SUBROUTINE CCEX_DIIS(CAMP,CPERT,FILDIISS,FILDIIST,AMAT,
     *                     NCCVAR,MXDIIS,NITER,WORK,LWORK)
C
C     CAMP  : Input  : R amplitudes from last iteration
C             Output : New R amplitudes
C
C     CPERT : Perturbative estimate for the R amplitudes.
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
C
      INTEGER NCCVAR, NITER, MXDIIS, LWORK
      INTEGER KCVEC, MXSIZE, KBVEC, KAMAT, KSCR1, KSCR2, KEND
      INTEGER NVEC, IVEC, LUSS, LUST, NSIZE, IFAIL
      CHARACTER FILDIISS*(*), FILDIIST*(*)
C
#if defined (SYS_CRAY)
      REAL CAMP(*), CPERT(*), AMAT(MXDIIS+1,MXDIIS+1)
      REAL WORK(LWORK), XMONE, ERROR, ONE, ZERO, DDOT, DNORM2
#else
      DOUBLE PRECISION CAMP(*), CPERT(*), AMAT(MXDIIS+1,MXDIIS+1)
      DOUBLE PRECISION WORK(LWORK), XMONE, ERROR, ONE, ZERO, DDOT, 
     *                 DNORM2
#endif
C
      PARAMETER (XMONE = -1.0D0, ONE = 1.0D0, ZERO = 0.0D0)
C
      CALL QENTER('CCEX_DIIS')
C
C------------------------------------
C     Allocation
C------------------------------------
C
      MXSIZE = MXDIIS + 1
C
      KCVEC  = 1
      KBVEC  = KCVEC + MXSIZE
      KAMAT  = KBVEC + MXSIZE
      KSCR1  = KAMAT + MXSIZE*MXSIZE
      KSCR2  = KSCR1 + MXSIZE
      KEND   = KSCR2 + MXSIZE
C
      IF (KEND .GT. LWORK) THEN
         CALL QUIT('Not enough mem in CCEX_DIIS')
      ENDIF
C
C--------------------------------------------------
C     Set CAMP = CPERT - CAMP
C--------------------------------------------------
C
      CALL DAXPY(NCCVAR,XMONE,CPERT,1,CAMP,1)
      CALL DSCAL(NCCVAR,XMONE,CAMP,1)
C
      NVEC = NITER
      IVEC = NITER
C
      if (niter .gt. mxdiis) then
         nvec = mxdiis
         ivec = niter - ((niter-1)/mxdiis)*mxdiis
      end if
C
C----------------------------------------------------
C     Open files and dump CPERT and CAMP:
C----------------------------------------------------
C
      LUSS = -1
      LUST = -1
      CALL WOPEN2(LUSS,FILDIISS,64,0)
      CALL WOPEN2(LUST,FILDIIST,64,0)
C
      CALL PUTWA2(LUST,FILDIIST,CPERT,1+(IVEC-1)*NCCVAR,NCCVAR)
      CALL PUTWA2(LUSS,FILDIISS,CAMP, 1+(IVEC-1)*NCCVAR,NCCVAR)
C
C----------------------------------------------------
C     Set up DIIS equations:
C----------------------------------------------------
C
      DO I = 1, NVEC
 
         CALL GETWA2(LUSS,FILDIISS,CPERT,1+(I-1)*NCCVAR,NCCVAR)
 
         error = ddot(nccvar,camp,1,cpert,1)
 
         amat(i,ivec)   = error
         amat(ivec,i)   = error
         amat(i,nvec+1) = -one
         amat(nvec+1,i) = -one
         work(kbvec-1+i)= zero
 
      end do
 
      amat(nvec+1,nvec+1) = zero
      work(kbvec+nvec)    = -one
 
      nsize = nvec + 1

*----------------------------------------------------------------------*
* solve the DIIS equations:
*----------------------------------------------------------------------*
      call f04atf(amat,mxsize,work(kbvec),nsize,work(kcvec),
     &            work(kamat),mxsize,work(kscr1),work(kscr2),ifail)
 
      if (ifail.ne.0) call quit('Error in ccex_diis: F04ATF failed.')
*----------------------------------------------------------------------*
* build the DIIS estimate for the new t amplitudes:
*----------------------------------------------------------------------*
      call dzero(camp,nccvar)
      do i = 1, nvec
         call getwa2(lust,fildiist,cpert,1+(I-1)*NCCVAR,NCCVAR)
         call daxpy(nccvar,work(kcvec-1+i),cpert,1,camp,1)
      end do
C
C----------------------------------------------------
C     Close files
C----------------------------------------------------
C
      CALL WCLOSE2(LUSS,FILDIISS,'KEEP')
      CALL WCLOSE2(LUST,FILDIIST,'KEEP')
C
C-----------------------------------------------
C     Finish.
C-----------------------------------------------
      CALL QEXIT('CCEX_DIIS')
C
      RETURN
      END
