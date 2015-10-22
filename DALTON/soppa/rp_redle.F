Cend-PFP
C
C  /* Deck rp_redle */
      SUBROUTINE RP_REDLE(NEXCI,NOLDTR,NNEWTR,LABEL,ISYMTR,REDE,LREDE,
     &                    REDS,LREDS,REDC,LREDC,LREDOL,FRVAL,NFRVAL,
     &                    IFREQ,ENORM,PROP,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     A.L. December 2002
C
C        Set up and solve the reduced linear resposne problem. Save the
C        new reduced matrices and the optimized trial and solution
C        vectors which may be written to disk.
C
C NEXCI         # of frequencies for that property (for us is alwyas 1)
C NOLDTR        # of old trial vectors
C NNEWTR        # of new trial vectors
C LABEL         name of the property
C ISYMTR        symmetry of the property
C REDE          E[2] reduced matrix
C LREDE         dimension of the E[2] reduced matrix
C REDS          S[2] reduced matrix
C LREDS         dimension of the S[2] reduced matrix
C LREDOL        dimension of the old E[2] and S[2] matrix
C FRVAL         array with the frequencies fot that property
C NFRVAL        dimension of the array FRVAL
C IFREQ        for which components of the FRVAL array we are computing
C ENORM         is the norm used to generate the new Z vector; i need that 
C               to compute the residual 
C WORK
C LWORK
C
#include "implicit.h"
#include "soppinf.h"
#include "ccsdsym.h"
#include "priunit.h"
C
      PARAMETER (ZERO = 0.0D+00, ONE = 1.0D0)
C
C--------------------------------
C     Dimensions of the arguments
C--------------------------------
C
      DIMENSION REDE(LREDE,LREDE),REDS(LREDS,LREDS),REDC(LREDC)
      DIMENSION FRVAL(NFRVAL)
cLig  <> I have to check if I really need to dimension these here!
cLig  and if it is true, I have to use the work array and check the 
cLig  memory requirment in teh calling subroutine (RP_LRSOEQ)
      DIMENSION INERT(3),DET(2),KPVT(LREDE)
      DIMENSION WORK(LWORK)
      CHARACTER*8 LABEL
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('RP_REDLE')
C
cLig  <> used in DSPLI 
      NSIM = 1
C
C-------------------------------------
C     Check dimensions of the matrices
C-------------------------------------
C
cLig  we can remove that check since the LREDE, LREDS and LREDC are
cLig  defined just few lines before the call to rp_redle
      IF (LREDE .NE. LREDS) THEN
       WRITE(LUPRO,*) 'RP_REDLE : Dimensions of reduced E[2] and S[2] ',
     &              'matrices LREDE : ',LREDE,' and LREDS : ',LREDS,
     &              ' are different'
         CALL QUIT('Incompatible dimensions of reduced matrices in //
     &             RP_REDLE.1')
      ENDIF
C
C---------------------------------------
C     Initialize size of reduced problem
C---------------------------------------
C
      NTRIAL = 2 * (NOLDTR + NNEWTR)
C
      IF (NTRIAL .GT. LREDE) THEN
       WRITE(LUPRI,*) 'RP_REDLE : Number of trial vectors ',
     &              NTRIAL,' exceeds dimensions of reduced ',
     &              'E[2] and S[2] matrices LREDE/LREDS : ',LREDE
         CALL QUIT('Dimensions of reduced matrices exceeded in //
     &             RP_REDLE.2')
      ENDIF
C
C===============================================
C     Set up the reduced linear response problem
C===============================================
C
C--------------------------------------------------------------------
C     Work space allocation no. 1.
C     Notice that the E[2] linear transformed trial vector and the 
C     S[2] linear transformed trial vector are of equal length and
C     that they use the same work space.
C--------------------------------------------------------------------
C
      LTR1E    = NT1AM(ISYMTR)
      LI       = NOLDTR + NNEWTR
      LJ       = NOLDTR + NNEWTR
C
      KIS      = 1
      KIE      = KIS + LI
      KIN      = KIE + LI
      KJS      = KIN + LI
      KJE      = KJS + LJ
      KJN      = KJE + LJ
      KEND1    = KJN + LJ
      LWORK1   = LWORK  - KEND1
C
      CALL SO_MEMMAX ('RP_REDLE.1',LWORK1)
      IF (LWORK1 .LT.0) CALL STOPIT('RP_REDLE.1',' ',KEND1,LWORK)
C
      IF ( IPRSOP .GE. 10 ) THEN
C
C---------------------------------------------
C        Print reduced E[2] and S[2] matrices.
C---------------------------------------------
C
         CALL AROUND('Reduced E[2] Matrix before RP_INCRED')
C
         CALL OUTPUT(REDE,1,LREDOL,1,LREDOL,LREDOL,LREDOL,1,LUPRI)
C
C
         CALL AROUND('Reduced S[2] Matrix before RP_INCRED')
C
         CALL OUTPUT(REDS,1,LREDOL,1,LREDOL,LREDOL,LREDOL,1,LUPRI)
C
      END IF
C
C---------------------------------------------------
C     Calculate new elements of the reduced matrices
C---------------------------------------------------
C
      DTIME      = SECOND()
      CALL RP_INCRED(NOLDTR,NNEWTR,ISYMTR,REDE,LREDE,REDS,LREDS,
     &               LREDOL,LTR1E,WORK(KIS),WORK(KIE),
     &               WORK(KIN),LI,WORK(KJS),WORK(KJE),WORK(KJN),LJ,
     &               WORK(KEND1),LWORK1)
      DTIME      = SECOND()   - DTIME
      SOTIME(28) = SOTIME(28) + DTIME
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
C-----------------------------------------------------
C     Calculate new elements of the reduced GP vectors
C-----------------------------------------------------
C
cLig  
      DTIME      = SECOND()
      CALL RP_REDGP(NOLDTR,NNEWTR,ISYMTR,REDC,LREDC,WORK(KEND1),LWORK1)
      DTIME      = SECOND()   - DTIME
      SOTIME(28) = SOTIME(28) + DTIME
C
      IF ( IPRSOP .GE. 5 ) THEN
C
C------------------------------
C        Print reduced C vector
C------------------------------
C
         CALL AROUND('Reduced C Vector')
C
         CALL OUTPUT(REDC,1,NTRIAL,1,1,LREDC,LREDC,1,LUPRI)
C
      END IF
C
C
C===============================================
C     Solve the reduced linear response problem.
C===============================================
C
C----------------------------------
C     Work space allocation no. 2.
C----------------------------------
C
      LEIVEC  = LREDE
C     
      KEIVEC  = 1
      KEMOS   = KEIVEC + LEIVEC*LEIVEC
      KZETA   = KEMOS  + (LREDE*(LREDE+1))/2
      KEND2   = KZETA  + LEIVEC*LEIVEC
      LWORK2  = LWORK  - KEND2
C
      CALL SO_MEMMAX ('RP_REDLE.2',LWORK2)
      IF (LWORK2 .LT.0) CALL STOPIT('RP_REDLE.2',' ',KEND2,LWORK)
C
C-------------------------------------------------------
C     Use EISPACK routine for real matrices 
C     in packed form 
C-------------------------------------------------------
C
      IJ = 0
      DO 100 JCOMP = 1,LREDE
C
        DO 200 ICOMP = 1,JCOMP
C
          IJ = IJ +1
          WORK(KEMOS+IJ-1) = REDE(ICOMP,JCOMP) - 
     &                       FRVAL(IFREQ) * REDS(ICOMP,JCOMP)
C
  200   CONTINUE
C
  100 CONTINUE
C
C---------------------------------------------------
C     Solve the LR problem Z[1]R = REDP-1 * C[1]R 
C---------------------------------------------------
C
cLig  I'll copy REDC on WORK(KZETA) since I need to compute GP * SV
      CALL DZERO(WORK(KZETA),LEIVEC*LEIVEC)
C
      CALL DCOPY(LEIVEC,REDC,1,WORK(KZETA),1)
C
      CALL DSPSLI(LREDE,NSIM,WORK(KEMOS),REDC,KPVT,INFO,DET,INERT)
C
      IF(IPRSOP.GE.7) THEN
C
              CALL AROUND("Reduced Solution Vector after DSPSLI")
C
              CALL OUTPUT(REDC,1,LREDE,1,1,1,LREDE,1,LUPRI)
C
              POL = DDOT(LREDE,REDC,1,WORK(KZETA),1)
C
              WRITE(LUPRI,'(/,1X,A,F15.8,A,F12.5,/)') "REDUCED GP * SV",
     &         POL," AT FREQUENCY ", FRVAL(IFREQ)
      ENDIF
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
      KOVLM   = KEND2
      KPEIV   = KOVLM  + LEIVEC*LEIVEC
      KEND3   = KPEIV  + LEIVEC
      LWORK3  = LWORK  - KEND3
C
      CALL SO_MEMMAX ('RP_RELDE.3',LWORK3)
      IF (LWORK3 .LT.0) CALL STOPIT('RP_REDLE.3',' ',KEND3,LWORK)
C
cLig  I'll copy the th solution vector REDC in WORK(KZETA)
      CALL DZERO(WORK(KZETA),LEIVEC*LEIVEC)
C
      CALL DCOPY(LEIVEC,REDC,1,WORK(KZETA),1)
C
      DTIME      = SECOND()
      CALL SO_ROPT(NOLDTR,NNEWTR,NEXCI,MXDIM,WORK(KZETA),LEIVEC,
     &             FRVAL,NFRVAL,WORK(KOVLM),WORK(KPEIV),REDE,LREDE,
     &             REDS,LREDS,.TRUE.,ENORM,WORK(KEND3),LWORK3)
      DTIME      = SECOND()   - DTIME
      SOTIME(38) = SOTIME(38) + DTIME
C
C-----------------------------------------------------------------
C     Calculate the value of the property with the trial vectors used for the
C     linear problem before they will be changed by the RP_OPTVEC
C     and write it on the file      
C-----------------------------------------------------------------
C
      DTIME      = SECOND()
      CALL RP_PROP(NOLDTR,NNEWTR,ISYMTR,REDC,LREDC,LABEL,FRVAL(IFREQ),
     &             PROP,WORK(KEND2),LWORK2)
      DTIME      = SECOND()   - DTIME
cLig  I'm quite sure to have to change this SOTIME 
      SOTIME(38) = SOTIME(38) + DTIME
cLig >>
C
C
C-----------------------------------------------------------------
C     Calculate the new optimal trial and linear transformed trial
C     (solution) vectors
C-----------------------------------------------------------------
C
      DTIME      = SECOND()   - DTIME
      CALL RP_OPTVEC(NOLDTR,NNEWTR,NEXCI,WORK(KZETA),LREDC,ISYMTR,
     &               WORK(KEND2),LWORK2)
      DTIME      = SECOND()   - DTIME
      SOTIME(37) = SOTIME(37) + DTIME
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('RP_REDLE')
C
      RETURN
      END
