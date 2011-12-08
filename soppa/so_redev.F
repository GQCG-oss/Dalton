C  /* Deck so_redev */
      SUBROUTINE SO_REDEV(NEXCI,NOLDTR,NNEWTR,ISYMTR,REDE,LREDE,REDS,
     &                    LREDS,LREDOL,EXVAL,LEXVAL,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Stephan Sauer, May 1996
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C     Andrea Ligabue, December 2003: linear response functions 
C                                    implemented
C
C        Set up and solve the reduced eigenvalue problem. Save the
C        new reduced matrices and the optimized trial and solution
C        vectors which may be written to disk.
C
C NEXCI           # of excitation energies to be calculated in the
C                 given symmetry
C NOLDTR          # of old trial vectors
C NNEWTR          # of new trial vectors
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C
#include <implicit.h>
#include <priunit.h>
#include <soppinf.h>
#include <ccsdsym.h>
C
      PARAMETER (ZERO = 0.0D+00, ONE = 1.0D0)
C
C--------------------------------
C     Dimensions of the arguments
C--------------------------------
C
      DIMENSION REDE(LREDE,LREDE),REDS(LREDS,LREDS)
      DIMENSION EXVAL(LEXVAL)
      DIMENSION WORK(LWORK)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_REDEV')
C
C-------------------------------------
C     Check dimensions of the matrices
C-------------------------------------
C
      IF (LREDE .NE. LREDS) THEN
         WRITE(LUPRI,*) 
     &        'SO_REDEV : Dimensions of reduced E[2] and S[2] matrices',
     &        ' LREDE : ',LREDE,' and LREDS : ',LREDS,' are different'
         CALL QUIT('Incompatible dimensions of reduced matrices in //
     &             SO_REDEV.1')
      ENDIF
C
C---------------------------------------
C     Initialize size of reduced problem
C---------------------------------------
C
      NTRIAL = 2 * (NOLDTR + NNEWTR)
C
      IF (NTRIAL .GT. LREDE) THEN
         WRITE(LUPRI,*) 'SO_REDEV : Number of trial vectors ',
     &              NTRIAL,' exceeds dimensions of reduced ',
     &              'E[2] and S[2] matrices LREDE/LREDS : ',LREDE
         CALL QUIT('Dimensions of reduced matrices exceeded in //
     &             SO_REDEV.2')
      ENDIF
C
C==========================================
C     Set up the reduced eigenvalue problem
C==========================================
C
C--------------------------------------------------------------------
C     Work space allocation no. 1.
C     Notice that the E[2] linear transformed trial vector and the 
C     S[2] linear transformed trial vector are of equal length and
C     that they use the same work space.
C--------------------------------------------------------------------
C
      LTR1E    = NT1AM(ISYMTR)
      LTR2E    = NT2AM(ISYMTR)
      LI       = NOLDTR + NNEWTR
      LJ       = NOLDTR + NNEWTR
C
      KIS     = 1
      KIE     = KIS + LI
      KIN     = KIE + LI
      KJS     = KIN + LI
      KJE     = KJS + LJ
      KJN     = KJE + LJ
      KEND1   = KJN + LJ
      LWORK1  = LWORK  - KEND1
C
      CALL SO_MEMMAX ('SO_SRES',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('SO_REDEV.1',' ',KEND1,LWORK)
C
C---------------------------------------------------
C     Calculate new elements of the reduced matrices
C---------------------------------------------------
C
      DTIME      = SECOND()
      CALL SO_INCRED(NOLDTR,NNEWTR,ISYMTR,REDE,LREDE,REDS,LREDS,
     &               LREDOL,LTR1E,LTR2E,WORK(KIS),WORK(KIE),
     &               WORK(KIN),LI,WORK(KJS),WORK(KJE),WORK(KJN),LJ,
     &               WORK(KEND1),LWORK1)
      DTIME      = SECOND()   - DTIME
      SOTIME(28) = SOTIME(28) + DTIME
C
C
      IF ( IPRSOP .GE. 5 ) THEN
C
C---------------------------------------------
C        Print reduced E[2] and S[2] matrices.
C---------------------------------------------
C
         CALL AROUND('Reduced E[2] Matrix')
C
         CALL OUTPUT(REDE,1,NTRIAL,1,NTRIAL,LREDE,LREDE,1,LUPRI)
C
C
         CALL AROUND('Reduced S[2] Matrix')
C
         CALL OUTPUT(REDS,1,NTRIAL,1,NTRIAL,LREDS,LREDS,1,LUPRI)
C
      END IF
C
C==========================================
C     Solve the reduced eigenvalue problem.
C==========================================
C
C----------------------------------
C     Work space allocation no. 2.
C----------------------------------
C
      LEIVEC  = LREDE
      LALFAR  = NTRIAL
      LALFAI  = NTRIAL
      LBETA   = NTRIAL
      LREDES  = LREDE
      LREDSS  = LREDS
      LREDST  = LREDS
      LSCR1   = NTRIAL      
C     
      KEIVEC   = 1
      KEND2A   = KEIVEC + LEIVEC*LEIVEC
      LWORK2A  = LWORK  - KEND2A
C
      KALFAR  = KEND2A
      KALFAI  = KALFAR + LALFAR
      KBETA   = KALFAI + LALFAI
      KREDES  = KBETA  + LBETA
      KREDSS  = KREDES + LREDES*LREDES
      KREDST  = KREDSS
      KEND2B  = KREDSS + LREDSS*LREDSS   
C
      KSCR1   = KEND2A
      KEND2C  = KEND2A + LSCR1*LSCR1
C
      KEND2   = MAX(KEND2B,KEND2C)
      LWORK2  = LWORK  - KEND2
C
      CALL SO_MEMMAX ('SO_REDEV.2',LWORK2)
      IF (LWORK2 .LT. 0) CALL STOPIT('SO_REDEV.2',' ',KEND2,LWORK)
C
C-------------------------------------------------------
C     Use EISPACK routine for real general matrices 
C     in generalized eigenvalue problem.
C     The reduced E[2] and S[2] matrices are backuped to 
C     WORK(KREDES) and WORK(KREDSS).
C-------------------------------------------------------
C
      CALL DCOPY(LREDE*LREDE,REDE,1,WORK(KREDES),1)
      CALL DCOPY(LREDS*LREDS,REDS,1,WORK(KREDSS),1)
C
      MATZ  = 1
      CALL RGG(LREDE,NTRIAL,REDE,REDS,WORK(KALFAR),WORK(KALFAI),
     &         WORK(KBETA),MATZ,WORK(KEIVEC),IERR)
C
      IF (IERR.NE.0) THEN
         WRITE(LUPRI,*) ' SO_REDEV : Reduced eigenvalue problem not ',
     &              'converged.'
         CALL QUIT(' Reduced eigenvalue equation not converged in//
     &               SO_REDEV')
      END IF
C
      CALL DCOPY(LREDE*LREDE,WORK(KREDES),1,REDE,1)
      CALL DCOPY(LREDS*LREDS,WORK(KREDSS),1,REDS,1)
C
C-----------------------------------------------
C     Order reduced eigenvalues and eigenvectors 
C-----------------------------------------------
C
      DTIME      = SECOND()
      CALL SO_ORDEIG (NTRIAL,rede,REDS,LREDS,WORK(KALFAR),LALFAR,
     &                WORK(KALFAI),LALFAI,WORK(KBETA),LBETA,
     &                WORK(KEIVEC),LEIVEC,WORK(KREDST),LREDST,
     &                WORK(KEND2),LWORK2)
      DTIME      = SECOND()   - DTIME
      SOTIME(36) = SOTIME(36) + DTIME
C
C----------------------------------------------
C     Save the lowest NEXCI excitation energies
C----------------------------------------------
C
      DO I = 1, NEXCI
         EXVAL(I) = WORK(KALFAR-1+I)
      END DO
C
C-------------------------------------------------------------------
C     Calculate orthonormalized reduced vectors spanning the optimum
C     space.
C-------------------------------------------------------------------
C
      MXDIM = NEXCI * NSAVMX
C
C---------------------------------
C     Work space allocation no. 3.
C---------------------------------
C
      KOVLM   = KEND2A
      KPEIV   = KOVLM  + LEIVEC*LEIVEC
      KEND3   = KPEIV  + LEIVEC
      LWORK3  = LWORK  - KEND3
C
      CALL SO_MEMMAX ('SO_REDEV.3',LWORK3)
      IF (LWORK3 .LT. 0) CALL STOPIT('SO_REDEV.3',' ',KEND3,LWORK)
C
      DTIME      = SECOND()
      CALL SO_ROPT(NOLDTR,NNEWTR,NEXCI,MXDIM,WORK(KEIVEC),LEIVEC,
     &             EXVAL,LEXVAL,WORK(KOVLM),WORK(KPEIV),REDE,LREDE,
     &             REDS,LREDS,.FALSE.,DUMMY,WORK(KEND3),LWORK3)
      DTIME      = SECOND()   - DTIME
      SOTIME(38) = SOTIME(38) + DTIME
C
C-----------------------------------------------------------------
C     Calculate the new optimal trial and linear transformed trial
C     (solution) vectors
C-----------------------------------------------------------------
C
      DTIME      = SECOND()   - DTIME
      CALL SO_OPTVEC(NOLDTR,NNEWTR,NEXCI,WORK(KEIVEC),LEIVEC,ISYMTR,
     &               WORK(KEND2A),LWORK2A)
      DTIME      = SECOND()   - DTIME
      SOTIME(37) = SOTIME(37) + DTIME
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_REDEV')
C
      RETURN
      END
