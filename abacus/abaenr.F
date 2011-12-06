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
c===========================================================================
c/* old abaenr_log */
c920325-hjaaj: ANRCTL: increment NWARN if not converged
c   ANRNEX: fixed for JCONV .lt. 0 as input parameter
c920302-hjaaj: corrected NCONST checks (was .gt. 1, should be .gt. 0)
c910108-prt: changing some DCOPY calls in ANRNEX, since with overlapping
c            vectors the local (not system) routine misbehaves when
c            compiled with vectorization and optimization on the CRAY
c900216-hjaaj: changed names from NR... to ANR... (A for Abacus).
c- save perturbation symmetry on LUNR1 for restart, new restart label.
c===========================================================================
C  /* Deck anrctl */
      SUBROUTINE ANRCTL(LMAXVE,REDH,REDGD,NREDH,NGD,NLOAD,DOXKAP,
     &                  SOLEQ,CMO,CREF,GORB,DV,PV,FC,FV,FCAC,H2AC,
     &                  DIACI,DIAOR,INDXCI,IBNDX,KNR,GD,
     &                  WRK,LWRK,MAXRED)
C
C Written Jan-1985 by Poul Joergensen
C Revisions: Jul-1988 HJJ (TUH) - new CI interface
C
C This subroutine directs the solutions of NGR simultaneous
C set of Newton equations
C
C (A-B)*X(J)+F(J)=0
C
C KNRSIM: Desired number of simultaneous roots in a
C         microiteration.  When convergence is
C         obtained for one root KNRSIM is decreased
C         by one.
C JNRSIM: is KNRSIM where linear dependence is removed
C         from the new trial vectors
C MAXNR:  maximum number of microiterations
C
#include <implicit.h>
#include <priunit.h>
#include <ibndxdef.h>
#include <iratdef.h>
#include <maxash.h>
#include <maxorb.h>
      PARAMETER (D0 = 0.0D0)
C
      DIMENSION REDH(*),REDGD(*),SOLEQ(*),CMO(*),CREF(*),GORB(*)
      DIMENSION DV(*),PV(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION DIACI(*),DIAOR(*),INDXCI(*)
      DIMENSION IBNDX(*),KNR(*),GD(NVARPT,*),WRK(*)
C
C Used from common blocks:
C   ANRINF : MAXNR, IPRNR, ?
C   INFDIM : LACIMX, LBCIMX,?
C   INFLIN : NCONST,...
C   INFPRI : NWARN, IPRKAP
C
#include <infvar.h>
#include <anrinf.h>
#include <inftap.h>
#include <infdim.h>
#include <inforb.h>
#include <infind.h>
#include <infpri.h>
#include <inflin.h>
#include <linaba.h>
C
      LOGICAL DOXKAP
      LOGICAL BINMEM
C
      CALL QENTER('ANRCTL')
C
C
      IF (IPRNR.GT.100) THEN
         WRITE(LUPRI,'(/A)')' ***** ENTERING ANRCTL **** '
         IF (NASHT.GT.0) THEN
            WRITE(LUPRI,'(/A)')
     *      ' ONE ELECTRON ACTIVE DENSITY MATRIX  PACKED'
            CALL OUTPAK(DV,NASHT,1,LUPRI)
            WRITE(LUPRI,'(/A)') ' FOCK VALENCE MATRIX'
            CALL OUTPKB(FV,NORB,NSYM,1,LUPRI)
         END IF
         WRITE(LUPRI,'(/A)') ' FOCK CORE MATRIX'
         CALL OUTPKB(FC,  NORB,NSYM,1,LUPRI)
      END IF
C
      IF (IPRNR .GT. 0) CALL TIMER ('START ',TIMEIN,TIMOUT)
      IPRKAP = IPRNR
      IF (L2PRI) THEN
         CALL ABALMT(IPRNR,REDH,REDGD,NREDH,NGD,SOLEQ,CMO,
     *               CREF,GORB,DV,PV,FC,FV,FCAC,H2AC,
     *               INDXCI,GD,WRK,LWRK)
      END IF
C
C LMAXVE: maximum number of simultaneous vectors in ABALIN
C
      NCONSX=MAX(4,NCONST)
      NWOPPX=MAX(4,NWOPPT)
C
C Initialize pointer to the NGD gradient vectors
C
      KNRSIM = NGD
      DO 95 I = 1,KNRSIM
95       KNR(I) = I
C
      IF (IPRNR.GT.30) THEN
         CALL AROUND('Test output of GD vectors from ANRCTL')
         DO 950 J = 1,NGD
            WRITE(LUPRI,'(/I3,A,I3)') J,'. GD gradient vector of',NGD
            WRITE(LUPRI,'(2X,10F7.4)') (GD(I,J),I=1,NVARPT)
  950    CONTINUE
      END IF
C
C CONSTRUCT NREDH START VECTORS
C
      IF (IPRNR.GT.15) THEN
         WRITE(LUPRI,'(/A,I8)')' NLOAD:',NLOAD
         WRITE(LUPRI,'(/A,I8,A,I8)')
     *   '    NREDH:',NREDH,'     LMAXVE:',LMAXVE
      END IF
      REWIND LUNR3
      IF (LOFFRD.EQ.1) READ (LUNR3)
      IF (NLOAD.EQ.1) THEN
         IF (NREDH.LE.LMAXVE) THEN
            BINMEM=.TRUE.
            NCSIM = 0
            NOSIM = 0
            DO 94 I = 1,NREDH
               IF (IBNDX(LOFFRD+I).EQ.JBCNDX) THEN
                  NCSIM = NCSIM + 1
               ELSE
                  NOSIM = NOSIM + 1
               ENDIF
   94       CONTINUE
            IF (IPRNR.GT.15) WRITE(LUPRI,'(/A,I5,A,I5)')
     *         ' NCSIM:',NCSIM,'    NOSIM:',NOSIM
            IF (IPRNR.GT.20) THEN
               WRITE(LUPRI,'(/A)')' IBNDX(LOFFRD+I),I=1,NREDH'
               WRITE(LUPRI,'(5I5,5X,5I5)')
     &            (IBNDX(LOFFRD+I),I=1,NREDH)
            END IF
            ISTBC = 1
            ISTBO = ISTBC + NCSIM*NCONST
            DO 96 I=1,NREDH
               IF (IBNDX(LOFFRD+I).EQ.JBCNDX) THEN
                  CALL READT(LUNR3,NCONST,WRK(ISTBC))
                  ISTBC=ISTBC+NCONST
               ELSE
                  CALL READT(LUNR3,NWOPPT,WRK(ISTBO))
                  ISTBO=ISTBO+NWOPPT
               ENDIF
 96         CONTINUE
         ELSE
            BINMEM=.FALSE.
         ENDIF
         IF (IPRNR.GT.30 .AND. BINMEM) THEN
            IOFF = 0
            DO 960 J = 1,NCSIM
               WRITE(LUPRI,'(/I3,A,I3)') J,'. start CI vector of',NCSIM
               WRITE(LUPRI,'(2X,10F7.4)') (WRK(IOFF+I),I=1,NCONST)
               IOFF = IOFF + NCONST
  960       CONTINUE
            DO 965 J = 1,NOSIM
               WRITE(LUPRI,'(/I3,A,I3)') 
     *             J,'. start orbital vector of',NOSIM
               WRITE(LUPRI,'(2X,10F7.4)') (WRK(IOFF+I),I=1,NWOPPT)
               IOFF = IOFF + NWOPPT
  965       CONTINUE
         END IF
      ELSE
         DO 97 I=1,NREDH
            IF (IBNDX(LOFFRD+I).EQ.JBONDX) THEN
               NTOT=NWOPPT
               IOFF=1+NCONST
            ELSE
               NTOT=NCONST
               IOFF=1
            ENDIF
            CALL READT(LUNR3,NTOT,WRK(1))
            DO 98 J=1,NGD
               REDGD((J-1)*NREDH+I)=-DDOT(NTOT,GD(IOFF,J),1,WRK(1),1)
 98         CONTINUE
 97      CONTINUE
      ENDIF
C
C
C
      TIMMIC = SECOND()
      TIMLIN = D0
      TIMRED = D0
      TIMNEX = D0
      MLIN  = 0
      MCLIN = 0
      MOLIN = 0
      MRED  = 0
      MNEX  = 0
C
C
C
      IF (NLOAD .EQ. 1) THEN
         JNRSIM = NREDH
         NREDH = 0
         ITNR = 0
      ELSE
         ITNR = 1
         ICTL = 2
         KBCNRV = 1
         KBONRV = 1
         KWRK1  = 1
         KIWRK  = KWRK1 + IROW(NREDH+1)
         LWRKNR = KIWRK + (NREDH+1)
         IF (LWRKNR.GT.LWRK) CALL ERRWRK('ANRCTL.ANRRED',LWRKNR,LWRK)
         IF (IPRNR .GT. 0) THEN
            WRITE(LUPRI,105) ITNR,0,0
            WRITE(LUPRI,'(A/)') ' (Restart with old trial vectors.)'
         END IF
         GO TO 810
C        ------------v
      END IF
 100  CONTINUE
      IF ( IPRNR .GT. 5) THEN
         WRITE(LUPRI,'(A,3I5)')' JNRSIM,ITNR,NLOAD',JNRSIM,ITNR,NLOAD
         WRITE(LUPRI,'(A,3I5)')' LMAXVE,NREDH,LOFFRD'
     *                         ,LMAXVE,NREDH,LOFFRD
         WRITE(LUPRI,'(A,L8)')' BINMEM',BINMEM
      END IF
      IF (IPRNR.GT.20) THEN
         WRITE(LUPRI,'(/A)')' IBNDX(I),I=1,(LOFFRD+NREDH+JNRSIM)'
         WRITE(LUPRI,'(5I5,5X,5I5)')(IBNDX(I),I=1,LOFFRD+NREDH+JNRSIM)
      END IF
      ITNR  = ITNR + 1
      IF (IPRNR .GE. 2) TIMIT = SECOND()
C
C     LOOP OVER BATCHES OF SIMULTANEOUS TRIAL VECTORS
C
 200  CONTINUE
C
         NSIM =MIN(LMAXVE,JNRSIM)
         NCSIM=0
         NOSIM=0
         DO 380 ISIM=1,NSIM
            IF (IBNDX(NREDH+ISIM+LOFFRD).EQ.JBCNDX) THEN
               NCSIM=NCSIM+1
            ELSE
               NOSIM=NOSIM+1
            ENDIF
  380    CONTINUE
         IF (IPRNR .GT. 0) WRITE(LUPRI,105) ITNR,NCSIM,NOSIM
  105    FORMAT(//' **MICROITERATION NUMBER',I5,
     *           /' ----------------------------------------',
     *           /'   Number of CI trial vectors      :',I5,
     *           /'   Number of orbital trial vectors :',I5,/)
C
C        WORKSPACE FOR ABALIN AND ANRRED
C
C        ANRRED USE MAX(NCONST,NWOPPT) AS WORKING SPACE FOR BVECTORS
C
         KBCNRV=1
         KBONRV=KBCNRV+NCSIM*NCONST
         KWRK1 =KBCNRV+MAX((NCSIM*NCONST+NOSIM*NWOPPT),NWOPPT,NCONST)
         LWRK1 =LWRK-KWRK1
C
C        WORK SPACE IN ANRRED
C
         KIWRK  = KWRK1 + IROW(NREDH+1+NSIM)
         LWRKNR = KIWRK + (NREDH+1+NSIM)
         IF (LWRKNR.GT.LWRK) CALL ERRWRK('ANRCTL.ANRRED',LWRKNR,LWRK)
C
         IF (.NOT.BINMEM) THEN
            REWIND LUNR3
            IF (LOFFRD.EQ.1) READ (LUNR3)
            DO 210 I=1,NREDH
               READ (LUNR3)
 210        CONTINUE
            ISTBC=KBCNRV
            ISTBO=KBONRV
            DO 400 ISIM=1,NSIM
               IF(IBNDX(NREDH+ISIM+LOFFRD).EQ.JBCNDX) THEN
                  CALL READT(LUNR3,NCONST,WRK(ISTBC))
                  ISTBC=ISTBC+NCONST
               ELSE
                  CALL READT(LUNR3,NWOPPT,WRK(ISTBO))
                  ISTBO=ISTBO+NWOPPT
               ENDIF
 400        CONTINUE
         ENDIF
C
         MLIN  = MLIN + 1
         MCLIN = MCLIN + NCSIM
         MOLIN = MOLIN + NOSIM
         TIM  = SECOND()
         CALL ABALIN (NCSIM,NOSIM,WRK(KBCNRV),WRK(KBONRV),
     *                CMO,CREF,GORB,DV,PV,
     *                FC,FV,FCAC,H2AC,INDXCI,WRK(KWRK1),1,LWRK1)
         TIM = SECOND() - TIM
         TIMLIN = TIMLIN + TIM
         IF (IPRNR .GE. 2) THEN
            WRITE (LUPRI,'(/A,F12.2,/A,I4,A,I4,A/)')
     *      ' >>> TIME USED IN ABALIN ',TIM,
     *      '     FOR',NCSIM,' CI AND',NOSIM,' ORBITAL TRIAL VECTORS.'
         END IF
         IF (IPRNR .GE. 10) THEN
            WRITE(LUPRI,'(/A,3I8)')' NSIM,NOSIM,NCSIM',NSIM,NOSIM,NCSIM
            WRITE(LUPRI,'(/A,I8)')' NREDH',NREDH
         END IF
C        Write (A-B)*X on LUNR5
C
         REWIND LUNR5
         IF (LOFFRD.EQ.1) THEN
            IF (NCONST.GT.0) READ (LUNR5)
            IF (NWOPPT.GT.0) READ (LUNR5)
         END IF
         DO 215 I=1,NREDH
            IF (NCONST.GT.0) READ (LUNR5)
            IF (NWOPPT.GT.0) READ (LUNR5)
 215     CONTINUE
         ISTSC=KWRK1
         ISTSO=KWRK1+NVARPT*NCSIM
         DO 216 I=1,NSIM
            IF (IBNDX(NREDH+LOFFRD+I).EQ.JBCNDX) THEN
               CALL WRITT(LUNR5,NCONSX,WRK(ISTSC))
               IF (NWOPPT.GT.0)
     *            CALL WRITT(LUNR5,NWOPPX,WRK(ISTSC+NCONST))
               ISTSC=ISTSC+NVARPT
            ELSE
               IF (NCONST.GT.0) CALL WRITT(LUNR5,NCONSX,WRK(ISTSO))
               CALL WRITT(LUNR5,NWOPPX,WRK(ISTSO+NCONST))
               ISTSO=ISTSO+NVARPT
            ENDIF
 216     CONTINUE
         NREDH =NREDH +NSIM
         IF (NREDH .GT. MAXRED) THEN
            WRITE (LUPRI,'(3X,A/3X,A,I5/3X,A,I5/3X,A/3X,A)')
     &      'ANRCTL: Allocation for reduced space too small.',
     &      'Dimension allocated :',MAXRED,
     &      'Dimension needed now:',NREDH,
     &      'Increase dimension with e.g. 50% with the .MAXRED option'//
     &          ' under *RESPON',
     &      'in the ABACUS module'//
     &          ' (i.e. under **START, **EACH STEP, and **PROPERTIES)'
            CALL QUIT('MAXRED too small in ANRCTL.')
         END IF
         JNRSIM=JNRSIM-NSIM
         IF (JNRSIM.LE.0) GO TO 800
C                         ------------v
C
C        CALL ANRRED(1,..) Increase dimension of reduced linear equation
C
         TIMRED = TIMRED - SECOND()
         MRED = MRED + 1
         CALL ANRRED(1,REDH,REDGD,NREDH,NSIM,NGD,SOLEQ,GD,IBNDX,
     *               WRK(KBCNRV),WRK(KBONRV),WRK(KWRK1),WRK(KIWRK))
C        CALL ANRRED(ICTL,REDH,REDGD,NREDH,N,NGD,SOLEQ,GD,IBNDX,
C    *               BCNRVE,BONRVE,SNRVEC,IWRK)
         TIMRED = TIMRED + SECOND()
         GO TO 200
C     ^-----------
C
C     CALL ANRRED(,3,..) Increase dimension of reduced linear equation
C     and find solutions
C
  800 CONTINUE
      ICTL=3
C
  810 CONTINUE
      MRED = MRED + 1
      TIMRED = TIMRED - SECOND()
      CALL ANRRED(ICTL,REDH,REDGD,NREDH,NSIM,NGD,SOLEQ,GD,IBNDX,
     *            WRK(KBCNRV),WRK(KBONRV),WRK(KWRK1),WRK(KIWRK))
      TIMRED = TIMRED + SECOND()
C
C     Save restart information
C
      CALL ANRSAV (NREDH, IBNDX, REDH)
C
C     Create KNRSIM new linear independent trial vectors in ANRNEX()
C     from the solutions of the reduced Newton set of equations
C
      IF (ITNR .LT. MAXNR) THEN
         JCONV = 0
      ELSE
         JCONV = -1
      END IF
      MNEX = MNEX + 1
      TIMNEX = TIMNEX - SECOND()
      CALL ANRNEX(NREDH,KNRSIM,JNRSIM,LMAXVE,DOXKAP,JCONV,IBNDX,
     *     KNR,GD,DIAOR,DIACI,SOLEQ,CMO,GORB,DV,PV,FC,FV,
     *     WRK,LWRK)
      TIMNEX = TIMNEX + SECOND()
      IF (JNRSIM.LE.LMAXVE) THEN
         BINMEM=.TRUE.
      ELSE
         BINMEM=.FALSE.
      ENDIF
C     
      IF (IPRNR .GE. 2) THEN
         TIMIT = SECOND() - TIMIT
         WRITE (LUPRI,'(/A,F12.2,/)')
     *      ' >>> Time used in this microiteration',TIMIT
      END IF
C
      IF (JCONV.LT.0) THEN
C        (LINEAR DEPENDENCE BETWEEN NEW TRIAL VECTOR )
         NWARN = NWARN + 1
         WRITE (LUPRI,'(/A/A/A,I5)')
     *      ' *** ANRCTL WARNING : Microiterations stopped - ',
     *      ' Linear dependence between all new trial vectors.',
     *      ' Perturbation symmetry:',LSYMPT
         WRITE(LUPRI,'(/A,I3,A)')
     *      ' Microiterations converged for',NGD-KNRSIM,
     *      ' gradient vectors,',
     *      ' Microiterations not converged for',KNRSIM,' vectors.'
      ELSE IF (JCONV.GT.0) THEN
C        (CONVERGED)
         IF (IPRNR .GT. 0) WRITE(LUPRI,'(/A,I4,A)')
     *     ' Microiterations converged for all',NGD,' gradient vectors.'
      ELSE
C        (NOT CONVERGED)
         IF (ITNR .LT. MAXNR) THEN
            GO TO 100
         ELSE
            NWARN = NWARN + 1
            WRITE(LUPRI,'(/A,I5,A/A,I5)')
     *         ' *** ANRCTL WARNING: Maximum number of microiterations',
     *         ITNR,' reached.',
     *         ' Perturbation symmetry:',LSYMPT
            WRITE(LUPRI,'(/A,I3,A)')
     *         ' Microiterations converged for',NGD-KNRSIM,
     *         ' gradient vectors,',
     *         ' Microiterations not converged for',KNRSIM,' vectors.'
         END IF
      END IF
C
      TIMMIC = SECOND() - TIMMIC
      IF (IPRNR .GE. 3) THEN
         WRITE (LUPRI,'(//A)') ' Final NR vectors in reduced basis for:'
         IOFF = 0
         DO 930 I = 1,NGD
            WRITE (LUPRI,'(/A,I4/)') ' - for gradient vector no.',I
            WRITE (LUPRI,'(5F15.8)') (SOLEQ(IOFF+J),J=1,NREDH)
            IOFF = IOFF + NREDH
  930    CONTINUE
      END IF
      IF (IPRNR .GE. 1) THEN
         WRITE (LUPRI,'(//A,F10.2,A//A,2I4/)')
     *      ' >>> Total time in microiterations',TIMMIC,' s.',
     *      '   Total number of CI and orbital trial vectors :',
     *      MCLIN,MOLIN
        IF (IPRNR .GE. 2) THEN
         WRITE (LUPRI,'(T14,A)') 'ROUTINE   CALLS    TIME'
         WRITE (LUPRI,'(T14,A,I8,F10.2)')   'ABALIN',MLIN,TIMLIN
         WRITE (LUPRI,'(T14,A,I8,F10.2)')   'ANRRED',MRED,TIMRED
         WRITE (LUPRI,'(T14,A,I8,F10.2)')   'ANRNEX',MNEX,TIMNEX
        END IF
      END IF
      IF (IPRNR .GT. 0) CALL TIMER ('ANRCTL',TIMEIN,TIMOUT)
      CALL QEXIT('ANRCTL')
      RETURN
      END
C  /* Deck anrcon */
      SUBROUTINE ANRCON(NREDH,NRD,SOLEQ,RD,IBNDX,WRK)
C
C PURPOSE:
C  CALCULATE CONVERGED SOLUTIONS TO THE NRD SETS OF LINEAR
C  EQUATIONS IN RD
C
#include <implicit.h>
#include <priunit.h>
#include <ibndxdef.h>
#include <maxash.h>
#include <maxorb.h>
C
      DIMENSION SOLEQ(NREDH,*),RD(NVARPT,*),IBNDX(*),WRK(*)
C
#include <inftap.h>
#include <anrinf.h>
#include <infind.h>
#include <infvar.h>
#include <inflin.h>
#include <linaba.h>
#include <infpri.h>
C
      REWIND LUNR3
      IF (LOFFRD.EQ.1) READ (LUNR3)
      CALL DZERO(RD,NRD*NVARPT)
      DO 100 I=1,NREDH
         IF (IBNDX(LOFFRD+I).EQ.JBCNDX) THEN
            CALL READT(LUNR3,NCONST,WRK)
            DO 200 K=1,NRD
               CALL DAXPY(NCONST,SOLEQ(I,K),WRK,1,RD(1,K),1)
 200        CONTINUE
         ELSE
            CALL READT(LUNR3,NWOPPT,WRK)
            DO 210 K=1,NRD
               CALL DAXPY(NWOPPT,SOLEQ(I,K),WRK,1,RD(1+NCONST,K),1)
 210        CONTINUE
         ENDIF
 100  CONTINUE
      IF (IPRNR.GT.100) THEN
         CALL AROUND('Test output of solution vectors from ANRCON')
         CALL OUTPUT(RD,1,NVARPT,1,NRD,NVARPT,NRD,1,LUPRI)
      END IF
C
      RETURN
      END
C  /* Deck anrnex */
      SUBROUTINE ANRNEX(NREDH,KNRSIM,JNRSIM,MAXSIM,DOXKAP,JCONV,
     *                 IBNDX,KNR,GD,DIAOR,DIACI,SOLEQ,CMO,GORB,
     *                 DV,PV,FC,FV,WRK,LWRK)
C
C PURPOSE: 1)CONSTRUCT RESIDUAL (A-B)*X(I)+F
C            FOR SOLUTION X(I) OF REDUCED N-R EQUATIONS
C          2)TEST FOR CONVERGENCE OF N-R SOLUTIONS
C            CONVERGENCE CRITERIUM:
C            //((A-B)*X(I)+F// .LE. THRNR
C          3)USE GENERALIZED CONJUGATE GRADIENT ALGORITHM
C            TO DETERMINE NEXT GUESS OF TRIAL VECTORS
C            IN CSF SPACE AND ORBITAL OPTIMIZATION
C            IN NEXKAP TO DETERMINE TRIAL ORBITAL
C            VECTORS.EITHER THE ORBITAL OR THE CSF PARAMETERS
C            ARE IMPROVED IN ONE ITERATION NEVER BOTH.
C
C JCONV  input: if JCONV .lt. 0 do not calculate new trial vectors.
C       output: =  1 converged
C               =  0 not converged
C               = -1 not converged, linear dependency among all
C                                   trial vectors.
C
C MAXSIM is maximum simultaneous bnrvec in nrctl.
C
C PJ JAN 1985
C
#include <implicit.h>
#include <priunit.h>
#include <ibndxdef.h>
#include <dummy.h>
#include <maxash.h>
#include <maxorb.h>
      PARAMETER ( FKAPC = 0.4D0, FKAPO = 0.01D0, FKAPNR = 0.25D0)
      PARAMETER ( MAXTST=30, THRLDP = 1.D-20, DTEST=1.0D-4  )
      PARAMETER ( MAXGDM=15 ,D10=1.0D1 ,D1000=1.0D3 )
      PARAMETER ( DM1 =-1.0D0,D1=1.0D0 , D0=0.0D0 )
C
      LOGICAL DOXKAP
      DIMENSION IBNDX(*),KNR(*),GD(NVARPT,*),DIAOR(*),DIACI(*)
      DIMENSION SOLEQ(NREDH,*)
      DIMENSION CMO(*),GORB(*),DV(*),PV(*),FC(*),FV(*),WRK(*)
      DIMENSION TEST(MAXTST),THRKAP(MAXGDM),SHIFNR(MAXGDM)
C
      SAVE TEST
C
C Used from common blocks:
C   ANRINF : THRNR
C
#include <infvar.h>
#include <anrinf.h>
#include <inftap.h>
#include <infind.h>
#include <inflin.h>
#include <linaba.h>
#include <infpri.h>
C
C
C Space for ANRNEX:
C
C MAXNEX: MAXIMUM NUMBER OF SIMULTANEOUS VECTORS IN ANRNEX
C
      MAXNEX=MIN((LWRK-NVARPT)/(NVARPT+NWOPPT),KNRSIM)
C
C SPACE ALLOCATION
C
      KB    = 1
      KS    = KB + MAXNEX*NVARPT
      KWKAP = KS + NVARPT
C
      NEXEC = 0
      NTEST = 0
      IBOFF = 0
      NBNR3 = NREDH + LOFFRD
      NSNR5 = NBNR3
C     ... NSNR5 = no. of sigma vectors on LUNR5
      NLOAD = 0
      DO 4000 ISIMC = 1,KNRSIM,MAXNEX
         NLOAD = NLOAD + 1
         NBX   = MIN(MAXNEX,(KNRSIM+1-ISIMC))
C
C CONSTRUCT RESIDUAL IN BNRVEC
C RESIDUAL: (A-B)*X(I)+F
C
C PUT F IN BNRVEC
C
         DO 50 JR=1,NBX
            JROOTJ = KNR(IBOFF+JR)
            CALL DCOPY(NVARPT,GD(1,JROOTJ),1,WRK(KB+(JR-1)*NVARPT),1)
  50     CONTINUE
         NTOT = NBX*NVARPT
         CALL DZERO(WRK(KWKAP),NWOPPT*NBX)
C
C CONSTRUCT  (A-B)*X(I) WHERE X(I) IS THE SOLUTION TO THE
C I'TH SET OF NEWTON-RAPHSON EQUATIONS
C
         REWIND LUNR5
C
C FIRST VECTOR CORRESPOND TO THE CI VECTOR AT THE EXPANSION
C POINT
C
         IF (LOFFRD.EQ.1) THEN
            IF (NCONST.GT.0) READ (LUNR5)
            IF (NWOPPT.GT.0) READ (LUNR5)
         END IF
         DO 900 K=1,NREDH
            IF(NCONST.GT.0) THEN
               CALL READT(LUNR5,NCONST,WRK(KS))
            ELSE
               WRK(KS) = D0
            ENDIF
            IF (NWOPPT.GT.0) CALL READT(LUNR5,NWOPPT,WRK(KS+NCONST))
            DO 700 JR=1,NBX
               JROOTJ=KNR(IBOFF+JR)
               EVAL1=SOLEQ(K,JROOTJ)
               IF (IBNDX(K+LOFFRD).EQ.JBONDX)THEN
                  CALL DAXPY(NVARPT,EVAL1,WRK(KS),1,
     *                       WRK(KB+(JR-1)*NVARPT),1)
               ELSE
                  CALL DAXPY(NCONST,EVAL1,WRK(KS),1,
     *                       WRK(KB+(JR-1)*NVARPT),1)
                  CALL DAXPY(NWOPPT,EVAL1,WRK(KS+NCONST),1,
     *                       WRK(KWKAP+(JR-1)*NWOPPT),1)
               ENDIF
 700        CONTINUE
 900     CONTINUE
         DO 950 JR=1,NBX
            CALL DAXPY(NWOPPT,D1,WRK(KWKAP+(JR-1)*NWOPPT),1,
     *                 WRK(KB+(JR-1)*NVARPT+NCONST),1)
            JROOTJ = KNR(IBOFF+JR)
            CALL DAXPY(NWOPPT,D1,GD(1+NCONST,JROOTJ),1,
     *                 WRK(KWKAP+(JR-1)*NWOPPT),1)
 950     CONTINUE
C
C RESIDUAL IS NOW IN BNRVEC AND MODIFIED GRADIENT IN WRK
C
C TEST FOR CONVERGENCE OF THE N-R SOLUTIONS
C
         IF (IPRNR .GT. 1) WRITE (LUPRI,'(//A,T20,A,T40,A,T60,A/)')
     *      ' GRADIENT NO.','CSF RESIDUAL','ORBITAL RESIDUAL',
     *      'TOTAL RESIDUAL'
         NFIN =0
         NCSIM=0
         NOSIM=0
         DO 2000 JR=1,NBX
            JROOTJ=KNR(IBOFF+JR)
            QCNORM=DDOT(NCONST,WRK(KB+NFIN*NVARPT),1,
     *                         WRK(KB+NFIN*NVARPT),1)
            QONORM=DDOT(NWOPPT,WRK(KB+NCONST+NFIN*NVARPT),1,
     *                        WRK(KB+NCONST+NFIN*NVARPT),1)
            QBTEST=SQRT(QCNORM+QONORM)
            QCNORM=SQRT(QCNORM)
            QONORM=SQRT(QONORM)
            IF (IPRNR .GT. 1) WRITE(LUPRI,'(I10,T15,1P,3D20.8)')
     *         JROOTJ,QCNORM,QONORM,QBTEST
            IF (JROOTJ.LE.MAXTST) TEST(JROOTJ)=QBTEST
            IF (QBTEST.GT.THRNR) THEN
               NEXEC=NEXEC+1
               NFIN=NFIN+1
               THRKAP(NFIN)=MAX(FKAPC*QCNORM,FKAPO*QONORM,FKAPNR*THRNR)
               KNR(NEXEC)  =KNR(IBOFF+JR)
               IF (QONORM.GT.QCNORM) THEN
                  NOSIM=NOSIM+1
                  IBNDX(NBNR3+NFIN)=JBONDX
               ELSE
                  NCSIM=NCSIM+1
                  IBNDX(NBNR3+NFIN)=JBCNDX
               ENDIF
            ELSE
               IF (JR.LT.NBX) THEN
#if defined (VAR_OLDCODE)
                  CALL DCOPY(((NBX-JR)*NVARPT),WRK(KB+(NFIN+1)*NVARPT),
     *                       1,WRK(KB+NFIN*NVARPT),1)
                  CALL DCOPY(((NBX-JR)*NWOPPT),
     *                       WRK(KWKAP+(NFIN+1)*NWOPPT),
     *                       1,WRK(KWKAP+NFIN*NWOPPT),1)
#else
C
C      Use loops instead of DCOPY to avoid overlapping vector errors
C
                  DO 2200 IQ = 0, (NBX-JR)*NVARPT-1
                     WRK(KB+NFIN*NVARPT+IQ) =
     *                      WRK(KB+(NFIN+1)*NVARPT+IQ)
2200              CONTINUE
                  DO 2300 IQ = 0, (NBX-JR)*NWOPPT-1
                     WRK(KWKAP+NFIN*NWOPPT+IQ) =
     *                      WRK(KWKAP+(NFIN+1)*NWOPPT+IQ)
2300              CONTINUE
#endif
               END IF
            ENDIF
 2000    CONTINUE
C
C     Skip generation of new trial vectors if either all vectors
C     in this load are converged or requested by input value of JCONV
C
      IF (NFIN.EQ.0 .OR. JCONV.LT.0) GO TO 3999
C
C USE CONJUGATE GRAD ALGORITHM TO FORM NEW CSF TRIAL VECTORS
C OR SOLVE LINEAR EQUATIONS TO GET IMPROVED ORBITAL VECTORS
C
         DO 1300 JR=1,NFIN
            JROOTJ=KNR(NEXEC-NFIN+JR)
            IF (IBNDX(NBNR3+JR).EQ.JBONDX) THEN
               IF (DOXKAP) THEN
                  CALL DCOPY(NWOPPT,WRK(KWKAP+(JR-1)*NWOPPT),1,
     *                       WRK(KB+(JR-1)*NVARPT+NCONST),1)
               ELSE
                  DO 1350 I=1,NWOPPT
                     IOFFI=(KB-1+(JR-1)*NVARPT+NCONST) + I
                     IF (ABS(DIAOR(I)).LE.DTEST) THEN
                        WRK(IOFFI)=WRK(IOFFI)/SIGN(DTEST,DIAOR(I))
                     ELSE
                        WRK(IOFFI)=WRK(IOFFI)/DIAOR(I)
                     ENDIF
 1350             CONTINUE
               ENDIF
            ELSE
               DO 1400 I=1,NCONST
                  IF (ABS(DIACI(I)).LE.DTEST) THEN
                     WRK(KB-1+I+(JR-1)*NVARPT)=
     *               WRK(KB-1+I+(JR-1)*NVARPT)/SIGN(DTEST,DIACI(I))
                  ELSE
                     WRK(KB-1+I+(JR-1)*NVARPT)=
     *               WRK(KB-1+I+(JR-1)*NVARPT)/DIACI(I)
                  ENDIF
 1400          CONTINUE
            ENDIF
 1300    CONTINUE
C
C CHANGE COLUMNS TO GET THE FIRST NCIVEC VECTORS TO
C CONTAIN IMPROVED CSF VECTORS AND THE NEXT
C (NFIN-NCIVEC) VECTORS IMPROVED ORBITAL VECTORS
C
         DO 1600 JR=1,NCSIM
            IF (IBNDX(NBNR3+JR).EQ.JBONDX) THEN
               DO 1700 KR=(NCSIM+1),NFIN
                  IF (IBNDX(NBNR3+KR).EQ.JBCNDX) THEN
                     CALL DSWAP(NVARPT,WRK(KB+(JR-1)*NVARPT),1,
     *                               WRK(KB+(KR-1)*NVARPT),1)
                     IBNDX(NBNR3+JR)=JBCNDX
                     IBNDX(NBNR3+KR)=JBONDX
                     THRKAP(KR)=THRKAP(JR)
                     GO TO 1600
                  ENDIF
 1700          CONTINUE
            ENDIF
 1600    CONTINUE
         ISTBC=NCONST+KB
         DO 1800 I=2,NCSIM
#if defined (VAR_OLDCODE)
            CALL DCOPY(NCONST,WRK(KB+(I-1)*NVARPT),1,WRK(ISTBC),1)
#else
C
C      Use loops instead of DCOPY to avoid overlapping vector errors
C
            DO 2400 IQ = 0,NCONST-1
               WRK(ISTBC+IQ) = WRK(KB+(I-1)*NVARPT+IQ)
2400        CONTINUE
#endif
            ISTBC=ISTBC+NCONST
 1800    CONTINUE
         ISTBO=NCSIM*NCONST+KB
         DO 1850 I=(NCSIM+1),NFIN
#if defined (VAR_OLDCODE)
            CALL DCOPY(NWOPPT,WRK(KB+NCONST+(I-1)*NVARPT),1,
     *                        WRK(ISTBO),1)
#else
            DO 2500 IQ = 0,NWOPPT-1
               WRK(ISTBO+IQ) = WRK(KB+NCONST+(I-1)*NVARPT+IQ)
2500        CONTINUE
#endif
            ISTBO=ISTBO+NWOPPT
 1850    CONTINUE
C
C
         KBCVEC=KB
         KBOVEC=KBCVEC+NCSIM*NCONST
         KWRK1 =KBOVEC+NOSIM*NWOPPT
         LWRK1 = LWRK - KWRK1
C
C        Find optimal orbital trial vectors, if requested
C
         IF (DOXKAP .AND. NOSIM .GT. 0) THEN
            DO 1545 I=1,NOSIM
 1545          SHIFNR(I) = D0
            CALL NEXKAP(NOSIM,WRK(KBOVEC),NSNR5,LUNR3,LUNR5,IBNDX,1,
     *                  DUMMY,THRKAP(NCSIM+1),SHIFNR,DAMP,DIAOR,CMO,
     *                  GORB,DV,PV,FC,FV,WRK(KWRK1),LWRK1)
         END IF
C
C ORTHONORMALIZE TRIAL VECTORS AND
C REMOVE NUMERICALLY LINEAR DEPENDENT VECTORS
C
         KWRK2 = KWRK1 + NCSIM + NOSIM
         NCONSM = MAX(NCONST,1)
         NWOPPM = MAX(NWOPPT,1)
         CALL ORTBVC(NCSIM,NOSIM,NCONSM,NWOPPM,NBNR3,IBNDX,LUNR3,
     *               WRK(KBCVEC),WRK(KBOVEC),THRLDP,
     *               WRK(KWRK1),WRK(KWRK2))
         NTEST  = NTEST + NCSIM + NOSIM
C
 3999 CONTINUE
         IBOFF=IBOFF+NBX
 4000 CONTINUE
      KNRSIM=NEXEC
      JNRSIM=NTEST
C
C KNRSIM: THE NUMBER OF DESIRED NEW TRIAL VECTORS TO BE
C         CREATED IN THE NEXT MICROITERATION
C JNRSIM: THE ACTUAL NUMBER OF NEW TRIAL VECTORS IN THE
C         NEXT MICROITERATION
C
C ---
C
C Output:
C
      IF (IPRNR.GT.0) THEN
         IF (KNRSIM.GT.MAXTST ) WRITE (LUPRI,95) KNRSIM,MAXTST
         K = MIN(KNRSIM,MAXTST)
         IF (K .GT. 0) THEN
            WRITE (LUPRI,96) (KNR(I),TEST(KNR(I)),I=1,K)
         END IF
      ENDIF
 95   FORMAT(/' ** ANRNEX ** NOTE: NUMBER OF SIMULTANEOUS TRIAL',
     *        ' VECTORS',I5,/,' EXCEEDS THE PRINT MAXIMUM OF',I5/)
 96   FORMAT(//' Residual of Newton-Raphson vectors not converged:'
     *       //' Newton-Raphson vector no.       norm of residual'
     *        /,(I15,15X,1P,D15.8))
C
C
C All Newton-Raphson linear equations have converged
C
      IF (KNRSIM .EQ. 0) THEN
         IF (IPRNR .GT. 0) WRITE(LUPRI,5030)
         JCONV = 1
         GO TO 9999
      ENDIF
 5030 FORMAT(/' *** Response equations converged for this symmetry.')
C
C MICROITERATIONS NOT CONVERGED
C =============================
C
C Case A) DO NOT CALCULATE NEW TRIAL VECTORS
C
      IF (JCONV .LT. 0) THEN
         JCONV = 0
         GO TO 9999
      END IF
C
C Case B) LINEAR DEPENDENCE BETWEEN ALL TRIAL VECTORS
C
      IF (JNRSIM .EQ. 0) THEN
         WRITE(LUPRI,5010)
         JCONV=-1
         GO TO 9999
      ENDIF
 5010 FORMAT(/' *** Microiterations stopped due to linear ',
     *        'dependence between all new trial vectors')
C
C Case C) NEW TRIAL VECTORS AVAILABLE
C IF JNRSIM .LE. MAXSIM : PLACE BNRVEC CORRECTLY IN MEMORY
C (FIRST ALL CI THEN ALL ORBITAL).
C
      JCONV = 0
      IF (JNRSIM .LE. MAXSIM) THEN
         IF (NLOAD.EQ.1) THEN
            INEW = KBCVEC + NCSIM*NCONST
            IOLD = KBOVEC
            IF (INEW .NE. IOLD) THEN
#if defined (VAR_OLDCODE)
               CALL DCOPY((NOSIM*NWOPPT),WRK(IOLD),1,WRK(INEW),1)
#else
C
C      Use loops instead of DCOPY to avoid overlapping vector errors
C
               DO 2600 IQ = 0,NOSIM*NWOPPT-1
                  WRK(INEW+IQ) = WRK(IOLD+IQ)
2600           CONTINUE
#endif
            ENDIF
         ELSE
            NCSIM = 0
            DO 800 I = (NREDH+1+LOFFRD),(NREDH+LOFFRD+JNRSIM)
  800          IF (IBNDX(I).EQ.JBCNDX) NCSIM = NCSIM + 1
            IBC = KBCVEC
            IBO = IBC + NCSIM*NCONST
            REWIND LUNR3
            IF (LOFFRD.EQ.1) READ (LUNR3)
            DO 810 I = 1,NREDH
  810          READ (LUNR3)
            DO 820 I = (NREDH+1+LOFFRD),(NREDH+LOFFRD+JNRSIM)
               IF (IBNDX(I).EQ.JBCNDX) THEN
                  CALL READT(LUNR3,NCONST,WRK(IBC))
                  IBC = IBC + NCONST
               ELSE
                  CALL READT(LUNR3,NWOPPT,WRK(IBO))
                  IBO = IBO + NWOPPT
               END IF
  820       CONTINUE
         END IF
      END IF
C
C END OF ANRNEX
C
 9999 RETURN
      END
C  /* Deck anrred */
      SUBROUTINE ANRRED (ICTL,REDH,REDGD,NREDH,N,NGD,SOLEQ,GD,
     *                  IBNDX,BCNRVE,BONRVE,SNRVEC,IWRK)
C
C Written JANUAR 85 PJ based on redmat
C
C INPUT:
C  ICTL, flow control
C         =1 ,extend the reduced PROJECTED HESSIAN and gradients
C         =2 ,find new solutions to reduced problem
C         =3 ,both
C  NREDH, dimension of new reduced PROJECTED HESSIAN matrix
C  N, number of new b-vectors
C  NGD, number of grad vectors
C  BCNRVE, the new csf b-vectors (is destroyed in NRRED)
C  BONRVE, the new orbital b-vectors (is destroyed in NRRED)
C  SNRVEC, the sigma vectors of the N new b-vectors
C
C  (The reduced PROJECTED HESSIAN matrix is the projection
C   on the basis of b-vectors.)
C
C Output:
C  REDH, the new, extended reduced PROJECTED HESSIAN matrix
c        (dimension: NREDH)
C  REDGD, the new, extended set of gradient vectors (dim:NREDH*NGD)
C  SOLEQ, solutions to the NGD set of NEWTON RAPHSON equations
C
C Scratch:
C  BCNRVE are used for B vectors from unit LUNR3 (section 1)
C  SNRVEC is used for reduced PROJECTED HESSIAN matrix (section 2)
C  IWRK, integer scratch array(dimension nredh)
C
#include <implicit.h>
#include <priunit.h>
#include <maxash.h>
#include <maxorb.h>
#include <ibndxdef.h>
      PARAMETER (D0 = 0.0D0)
C
      DIMENSION REDH(*),REDGD(*),SOLEQ(NREDH,*),GD(NVARPT,*)
      DIMENSION IBNDX(*),BCNRVE(*),BONRVE(*),SNRVEC(NVARPT,*)
      DIMENSION IWRK(*)
C
#include <anrinf.h>
#include <inftap.h>
#include <infinp.h>
#include <infind.h>
#include <infvar.h>
#include <inflin.h>
#include <linaba.h>
#include <ibtdef.h>
C
C Section 1: extend reduced PROJECTED HESSIAN matrix
C            with N new b-vectors
C
      IF (IBTAND(ICTL,1).EQ.0) GO TO 1000
      IF (N.LT.1) GO TO 1000
      IREDH = NREDH - N
C
C New b-vectors (are in BNRVEC)
C SET UP REDUCED GRAD VECTOR
C
C REORDER PREVIOUS REDGD ELEMENTS
C
      NTOT=NGD*NREDH
      CALL DZERO(SOLEQ,NTOT)
      DO 14 I=1,NGD
         CALL DCOPY(IREDH,REDGD((I-1)*IREDH+1),1,SOLEQ(1,I),1)
 14   CONTINUE
C
C EXTEND SIZE OF REDGD VECTRORS
C
      DO 25 I=1,NGD
         ISTBC=1
         ISTBO=1
         DO 30 K=1,N
            IF(IBNDX(IREDH+K+LOFFRD).EQ.JBCNDX) THEN
               SOLEQ(IREDH+K,I) =
     *            -DDOT(NCONST,BCNRVE(ISTBC),1,GD(1,I),1)
               ISTBC=ISTBC+NCONST
            ELSE
               SOLEQ(IREDH+K,I) =
     *            -DDOT(NWOPPT,BONRVE(ISTBO),1,GD(NCONST+1,I),1)
               ISTBO=ISTBO+NWOPPT
            END IF
 30      CONTINUE
 25   CONTINUE
      CALL DCOPY(NTOT,SOLEQ,1,REDGD,1)
C
C EXTEND THE REDUCED PROJECTED HESSIAN MATRIX
C
      ISTBC = 1
      ISTBO = 1
      DO 35 K = 1,N
         KL = IROW(IREDH+K) + IREDH
         IF (IBNDX(IREDH+K+LOFFRD).EQ.JBCNDX) THEN
            DO 40 L=1,K
               REDH(KL+L) = DDOT(NCONST,BCNRVE(ISTBC),1,SNRVEC(1,L),1)
 40         CONTINUE
            ISTBC=ISTBC+NCONST
         ELSE
            DO 41 L=1,K
               REDH(KL+L) = DDOT(NWOPPT,BONRVE(ISTBO),1,
     *            SNRVEC(NCONST+1,L),1)
 41         CONTINUE
            ISTBO=ISTBO+NWOPPT
         END IF
   35 CONTINUE
C
C Old b-vectors.
C Rewind LUNR3 and skip CREF, the first b-vector
C (THE PROJECTED HESSIAN MATRIX HAS NO CO COMPONENT)
C
      REWIND LUNR3
      IF (LOFFRD.EQ.1) READ (LUNR3)
      LL = 0
      DO 15 J = 1,IREDH
         IF (IBNDX(J+LOFFRD).EQ.JBCNDX) THEN
            CALL READT(LUNR3,NCONST,BCNRVE)
            LL=LL+1
            DO 27 K=1,N
               REDH(LL+IROW(IREDH+K)) =
     *            DDOT(NCONST,BCNRVE,1,SNRVEC(1,K),1)
   27        CONTINUE
          ELSE
             CALL READT(LUNR3,NWOPPT,BCNRVE)
             LL=LL+1
             DO 26 K=1,N
                REDH(LL+IROW(IREDH+K)) =
     *             DDOT(NWOPPT,BCNRVE,1,SNRVEC(NCONST+1,K),1)
   26        CONTINUE
          ENDIF
   15  CONTINUE
C
C *********************************************************
C Section 2: FIND SOLUTIONS TO THE NGR SET OF NEWTON-RAPHSON EQ.
C
 1000 CONTINUE
      IF (IBTAND(ICTL,2).EQ.0) GO TO 2000
C
C SOLVE NGD SET OF N-R EQ.
C Copy REDH to SNRVEC and REDGD to soleq
C
      KMA = IROW(NREDH+1)
      CALL DCOPY(KMA,REDH(1),1,SNRVEC,1)
      KMA = NREDH*NGD
      CALL DCOPY(KMA,REDGD(1),1,SOLEQ(1,1),1)
      CALL DSPSOL(NREDH,NGD,SNRVEC,SOLEQ,IWRK,INFO)
      IF (IPRNR.GT.5 .OR. INFO.NE.0) THEN
         WRITE(LUPRI,'(/A)')' *** REDUCED HESSIAN MATRIX ***'
         CALL OUTPAK(REDH(1),NREDH,1,LUPRI)
         WRITE(LUPRI,'(/A)')' *** REDUCED GRADIENT VECTOR ***'
         CALL OUTPUT(REDGD,1,NREDH,1,NGD,NREDH,NGD,1,LUPRI)
         WRITE(LUPRI,'(/A)')' * SOLUTION TO NEWTON RAPHSON EQUATIONS *'
         CALL OUTPUT(SOLEQ(1,1),1,NREDH,1,NGD,NREDH,NGD,1,LUPRI)
      END IF
      IF (INFO.NE.0)THEN
         WRITE(LUPRI,7000)
         CALL QUIT('ANRRED: NO SOLUTION FOUND TO LINEAR EQUATIONS')
      ENDIF
 7000 FORMAT(/' ***ANRRED*** SOLUTION NOT OBTAINED TO LINEAR EQUATIONS',
     *       /' CHECK IF HESSIAN MATRIX IS SINGULAR')
C
C *** End of subroutine ANRRED
C
 2000 CONTINUE
      RETURN
      END
C  /* Deck anrrst */
      LOGICAL FUNCTION ANRRST (NREDH, IBNDX, REDH)
C
C 13-Aug-1985 hjaaj
C
C Purpose: recover restart information in for ANRCTL module.
C
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
#include <maxorb.h>
      DIMENSION IBNDX(*), REDH(*)
C
C Used from common blocks:
C  INFVAR: NCONST
C  INFTAP: LUNR1, LUNR3, LUNR5
C
#include <infvar.h>
#include <inftap.h>
#include <inflin.h>
#include <linaba.h>
C
      LOGICAL LTEST, FNDLAB, AROUND
C
C DECLARE EXTERNAL FUNCTIONS
C
      LTEST = .FALSE.
      CALL GPINQ('ABANR1','EXIST',AROUND)
      IF (.NOT.AROUND) GO TO 1000
      CALL GPINQ('ABANR3','EXIST',AROUND)
      IF (.NOT.AROUND) GO TO 1000
      CALL GPINQ('ABANR5','EXIST',AROUND)
      IF (.NOT.AROUND) GO TO 1000
C
      CALL GPOPEN(LUNR1,ABANR1,'OLD',' ',' ',IDUMMY,.FALSE.)
      CALL GPOPEN(LUNR3,ABANR3,'OLD',' ',' ',IDUMMY,.FALSE.)
      CALL GPOPEN(LUNR5,ABANR5,'OLD',' ',' ',IDUMMY,.FALSE.)
C
      REWIND LUNR1
      LTEST = FNDLAB('ANRRESTR', LUNR1)
      IF (.NOT. LTEST) THEN
         WRITE (LUPRI,'(/A)') ' ANRRST: RESTART LABEL NOT FOUND'
         GO TO 1000
      END IF
      LTEST = .FALSE.
C
      READ (LUNR1) JSYMPT,NREDH, (IBNDX(I), I = 1,NREDH+1)
      IF (JSYMPT .NE. LSYMPT) THEN
         WRITE (LUPRI,'(/A,2(/A,I4))')
     &      ' ANRRST: Restart is not possible',
     &      ' Perturbation symmetry on restart file:',JSYMPT,
     &      ' Current perturbation symmetry        :',LSYMPT
         GO TO 1000
      END IF
C
C TEST IF RECORDS EXIST ON LUNR3 AND LUNR5
C
      REWIND LUNR3
      NREC3 = NREDH + LOFFRD
      DO 130 I = 1,NREC3
         READ (LUNR3,END=830,ERR=830) DUM
  130 CONTINUE
      REWIND LUNR3
      REWIND LUNR5
      NREC5 = NREC3
      IF (NCONST .GT. 0) NREC5 = 2*NREC5
      DO 150 I = 1,NREC5
         READ (LUNR5,END=850,ERR=850) DUM
  150 CONTINUE
      REWIND LUNR5
      LTEST = .TRUE.
C
      NNREDH = NREDH * (NREDH + 1) / 2
      READ (LUNR1) (REDH(I), I = 1,NNREDH)
      GO TO 1000
C
  830 WRITE (LUPRI,'(/A,I3)') ' ANRRST: error testing LUNR3 =',LUNR3
      GO TO 1000
  850 WRITE (LUPRI,'(/A,I3)') ' ANRRST: error testing LUNR5 =',LUNR5
C
 1000 ANRRST = LTEST
      CALL GPCLOSE(LUNR1,'KEEP')
      CALL GPCLOSE(LUNR3,'KEEP')
      CALL GPCLOSE(LUNR5,'KEEP')
      RETURN
      END
C  /* Deck anrsav */
      SUBROUTINE ANRSAV(NREDH, IBNDX, REDH)
C
C 13-Aug-1985 hjaaj
C
C Purpose: save restart information in for ANRCTL module.
C
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
      DIMENSION IBNDX(*), REDH(*)
C
C Used from common blocks:
C   INFTAP : LUNR1
C
#include <inflin.h>
#include <linaba.h>
#include <inftap.h>
C
      CHARACTER*8 TABLE(3)
      DATA TABLE /'********', 'ANRRESTR', 'EODATA  '/
C
      CALL GPOPEN(LUNR1,ABANR1,'UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      REWIND LUNR1
      WRITE (LUNR1) TABLE(1),TABLE(1),TABLE(1),TABLE(2)
      NBNDX = MAX(7,(NREDH+LOFFRD))
      WRITE (LUNR1) LSYMPT,NREDH,(IBNDX(I), I = 1,NBNDX)
      NNREDH = NREDH * (NREDH + 1) / 2
      WRITE (LUNR1) (REDH(I), I = 1,NNREDH)
      WRITE (LUNR1) TABLE(1),TABLE(1),TABLE(1),TABLE(3)
      CALL GPCLOSE(LUNR1,'KEEP')
      RETURN
      END
C  /* Deck abalmt */
      SUBROUTINE ABALMT(IPRNR,REDH,REDGD,NREDH,NGD,SOLEQ,CMO,
     *                  CREF,GORB,DV,PV,FC,FV,FCAC,H2AC,
     *                  INDXCI,GD,WRK,LWRK)
C
#include <implicit.h>
#include <priunit.h>
C
C     This is a test routine that calculates the MC electronic Hessian
C     matrix explicitly by carrying out linear transformations of unit
C     vectors
C
#include <inflin.h>
#include <linaba.h>
#include <infpri.h>
C
      PARAMETER ( D0 = 0.0D0 , D1 = 1.0D0 , DTOL = 1.0D-8 )
      DIMENSION REDH(*),REDGD(*),SOLEQ(*),CMO(*),CREF(*),GORB(*)
      DIMENSION DV(*),PV(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION INDXCI(*),GD(*),WRK(*)
C
      CALL AROUND('Test output from ABALMT')
C
C     Allocate work space
C
      NCOMAX = MAX(NCONST,NWOPPT)
      KL2    = 1
      KL2END = KL2 + NVARPT*NVARPT
      KBVEC  = KL2END
      KWRK1  = KBVEC  + NCOMAX
      LWRK1  = LWRK   - KWRK1
      IF (KWRK1 .GT. LWRK) CALL STOPIT('ABALMT','MC Hess',-KWRK1,LWRK)
C
      CALL DZERO(WRK(KBVEC),NCOMAX)
      DO 100 I = 1,NVARPT
         IF (I.LE.NCONST) THEN
            NCSIM = 1
            NOSIM = 0
            IOFF  = I
         ELSE
            NCSIM = 0
            NOSIM = 1
            IOFF  = I - NCONST
         END IF
         WRK(KBVEC-1+IOFF) = D1
         IF (( NCSIM.GT.0 ).AND.( LSYMST.EQ.LSYMRF)) THEN
            CALL DAXPY(NCONST,-CREF(IOFF),CREF,1,WRK(KBVEC),1)
         END IF
         IF (IPRNR .GT. 90) THEN
            IF(NOSIM.GT.0)  THEN
               NDIM = NWOPPT
               CALL HEADER('Orbital trial vector in ABALMT',-1)
            END IF
            IF(NCSIM.GT.0) THEN
               NDIM = NCONST
               CALL HEADER('Configuration trial vector in ABALMT',-1)
            END IF
            WRITE (LUPRI,'(1X,A,I5)') ' This is vector No.',I
            CALL OUTPUT(WRK(KBVEC),1,NDIM,1,1,NDIM,1,1,LUPRI)
         ELSE IF (IPRNR .GT. 10) THEN
            IF(NOSIM.GT.0)  THEN
              WRITE (LUPRI,'(//A,2I5)')
     &        ' Orbital trial vector in ABALMT, vector no.',I,IOFF
            ELSE
              WRITE (LUPRI,'(//A,2I5)')
     &        ' Configuration trial vector in ABALMT, vector no.',I,IOFF
            END IF
         END IF
C
C        Carry out linear transformation
C
         CALL ABALIN(NCSIM,NOSIM,WRK(KBVEC),WRK(KBVEC),
     &               CMO,CREF,GORB,DV,PV,FC,FV,FCAC,H2AC,INDXCI,
     &               WRK(KWRK1),1,LWRK1)
C
C        Project out reference state components from linearly
C        transformed vectors
C
         IF ((.NOT.ABAHF ).AND.( LSYMST.EQ.LSYMRF)) THEN
            XL2OVL = DDOT(NCONST,CREF,1,WRK(KWRK1),1)
            CALL DAXPY(NCONST,-XL2OVL,CREF,1,WRK(KWRK1),1)
         END IF
C
         CALL DCOPY(NVARPT,WRK(KWRK1),1,WRK(KL2+(I-1)*NVARPT),1)
         IF (( NCSIM.GT.0 ).AND.( LSYMST.EQ.LSYMRF)) THEN
            CALL DZERO(WRK(KBVEC),NCONST)
         ELSE
            WRK(KBVEC-1+IOFF) = D0
         END IF
 100  CONTINUE
C
C     Investigate deviation from symmetry:
C
      CALL HEADER('Deviation from symmetry in MC electronic Hessian',-1)
      WRITE(LUPRI,'(A,I8)') ' Perturbation symm.:',LSYMPT
      WRITE(LUPRI,'(A,I8)') ' Reference symmetry:',LSYMRF
      WRITE(LUPRI,'(A,I8)') ' CI symmetry       :',LSYMST
      WRITE(LUPRI,'(A,I8)') ' CI       dimension:',NCONST
      WRITE(LUPRI,'(A,I8)') ' Orb.rot. dimension:',NWOPPT
      WRITE(LUPRI,'(A,I8)') ' Total    dimension:',NVARPT
C
C     Orbital - orbital block
C
      ZAMAX = D0
      IA    = 0
      JA    = 0
      DO 150 IZ = NCONST + 1,NVARPT
         DO 160 JZ = NCONST + 1,IZ-1
            ZAIJ  = WRK(KL2-1+(IZ-1)*NVARPT+JZ)
            ZAJI  = WRK(KL2-1+(JZ-1)*NVARPT+IZ)
            ZDEV  = ABS(ZAIJ-ZAJI)
            IF (ZDEV.GT.ZAMAX) THEN
               ZAMAX = ZDEV
               IA    = IZ
               JA    = JZ
            END IF
 160     CONTINUE
 150  CONTINUE
      IF (NWOPPT .GT. 1) THEN
         WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *   ' Maximum deviation in orbital-orbital block I=',
     *   IA,' J=',JA,' DEV=',ZAMAX
         IF ( ZAMAX.GT.DTOL) WRITE(LUPRI,'(/2(/A,I5,A,I5,A,1P,G16.8))')
     *      ' I=',IA,' J=',JA,' element',WRK(KL2-1+(IA-1)*NVARPT+JA),
     *      ' I=',JA,' J=',IA,' element',WRK(KL2-1+(JA-1)*NVARPT+IA)
      END IF
C
C     CI - CI block
C
      ZBMAX = D0
      IA    = 0
      JA    = 0
      DO 151 IZ = 1,NCONST
         DO 161 JZ = 1,IZ-1
            ZAIJ  = WRK(KL2-1+(IZ-1)*NVARPT+JZ)
            ZAJI  = WRK(KL2-1+(JZ-1)*NVARPT+IZ)
            ZDEV  = ABS(ZAIJ-ZAJI)
            IF (ZDEV.GT.ZBMAX) THEN
               ZBMAX = ZDEV
               IA    = IZ
               JA    = JZ
            END IF
 161     CONTINUE
 151  CONTINUE
      IF (NCONST .GT. 1) THEN
         WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *   ' Maximum deviation in CI-CI block           I=',
     *   IA,' J=',JA,' DEV=',ZBMAX
         IF ( ZBMAX.GT.DTOL) WRITE(LUPRI,'(/2(/A,I5,A,I5,A,1P,G16.8))')
     *      ' I=',IA,' J=',JA,' element',WRK(KL2-1+(IA-1)*NVARPT+JA),
     *      ' I=',JA,' J=',IA,' element',WRK(KL2-1+(JA-1)*NVARPT+IA)
      END IF
C
C     CI - orbital block
C
      ZCMAX = D0
      IA    = 0
      JA    = 0
      DO 152 IZ = 1 + NCONST,NVARPT
         DO 162 JZ = 1, NCONST
            ZAIJ  = WRK(KL2-1+(IZ-1)*NVARPT+JZ)
            ZAJI  = WRK(KL2-1+(JZ-1)*NVARPT+IZ)
            ZDEV  = ABS(ZAIJ-ZAJI)
            IF (ZDEV.GT.ZCMAX) THEN
               ZCMAX = ZDEV
               IA    = IZ
               JA    = JZ
            END IF
 162     CONTINUE
 152  CONTINUE
      IF (NWOPPT .GT. 0 .AND. NCONST .GT. 0) THEN
         WRITE(LUPRI,'(/A,I5,A,I5,A,1P,G16.8)')
     *   ' Maximum deviation in orbital-CI block      I=',
     *   IA,' J=',JA,' DEV=',ZCMAX
         IF ( ZCMAX.GT.DTOL) WRITE(LUPRI,'(/2(/A,I5,A,I5,A,1P,G16.8))')
     *      ' I=',IA,' J=',JA,' element',WRK(KL2-1+(IA-1)*NVARPT+JA),
     *      ' I=',JA,' J=',IA,' element',WRK(KL2-1+(JA-1)*NVARPT+IA)
      END IF
      ZMAX = MAX(ZAMAX,ZBMAX,ZCMAX)
      IF (IPRNR .GE. 5 .OR. ZMAX .GT. DTOL) THEN
         CALL AROUND('MC electronic Hessian in ABALMT')
         CALL OUTPUT(WRK(KL2),1,NVARPT,1,NVARPT,NVARPT,NVARPT,1,LUPRI)
      END IF
C
C     **************************************************
C     ***** Testing simultaneous vectors in ABALIN *****
C     **************************************************
C
      CALL AROUND('Simultaneous linear transformations - ABALMT')
C
C     Allocate work space
C
      IF (ABAHF) THEN
         NCSIM = 0
         NOSIM = MIN(2,NWOPPT)
      ELSE IF (ABACI) THEN
         NCSIM = MIN(2,NCONST)
         NOSIM = 0
      ELSE
         NCSIM = MIN(2,NCONST)
         NOSIM = MIN(2,NWOPPT)
      END IF
      KBCVEC = KL2END
      KBOVEC = KBCVEC + NCSIM*NCONST
      KWRK1  = KBOVEC + NOSIM*NWOPPT
      LWRK1  = LWRK   - KWRK1
      IF (KWRK1 .GT. LWRK) THEN
         CALL STOPIT('ABALMT','MC sim vec Hess',-KWRK1,LWRK)
      END IF
C
C     Orbital vectors
C
      CALL DZERO(WRK(KBOVEC),NOSIM*NWOPPT)
      IOFF = KBOVEC - 1
      DO 200 I = 1, NOSIM
         WRK(IOFF + I) = D1
         IOFF = IOFF + NWOPPT
  200 CONTINUE
C
C     Configuration vectors
C
      CALL DZERO(WRK(KBCVEC),NCSIM*NCONST)
      IOFF = KBCVEC
      DO 300 I = 1, NCSIM
         WRK(IOFF + I - 1) = D1
         IF (LSYMST .EQ. LSYMRF) THEN
            CALL DAXPY(NCONST,-CREF(I),CREF,1,WRK(IOFF),1)
         END IF
         IOFF = IOFF + NCONST
  300 CONTINUE
      IF (IPRNR .GT. 90) THEN
         CALL HEADER('Orbital trial vectors in ABALMT',-1)
         CALL OUTPUT(WRK(KBOVEC),1,NWOPPT,1,NOSIM,NWOPPT,NOSIM,1,LUPRI)
         IF (.NOT.ABAHF) THEN
            CALL HEADER('Configuration trial vectors in ABALMT',-1)
            CALL OUTPUT(WRK(KBCVEC),1,NCONST,1,NCSIM,
     &                  NCONST,NCSIM,1,LUPRI)
         END IF
      END IF
C
C     Carry out linear transformation
C
      CALL ABALIN(NCSIM,NOSIM,WRK(KBCVEC),WRK(KBOVEC),
     &            CMO,CREF,GORB,DV,PV,FC,FV,FCAC,H2AC,INDXCI,
     &            WRK(KWRK1),1,LWRK1)
C
C     Project out reference state components from linearly
C     transformed vectors and copy to final matrix
C
      IF ((.NOT.ABAHF ).AND.( LSYMST.EQ.LSYMRF)) THEN
         IOFF = 0
         DO 400 I = 1, NCSIM + NOSIM
            XL2OVL = DDOT(NCONST,CREF,1,WRK(KWRK1 + IOFF),1)
            CALL DAXPY(NCONST,-XL2OVL,CREF,1,WRK(KWRK1 + IOFF),1)
            IOFF = IOFF + NVARPT
  400    CONTINUE
      END IF
C
C     Print transformed vectors
C
      IF (IPRNR .GE. 10) THEN
         CALL AROUND('MC electronic Hessian  (sim. vectors) in ABALMT')
         CALL OUTPUT(WRK(KWRK1),1,NVARPT,1,NCSIM+NOSIM,NVARPT,
     &               NCSIM+NOSIM,1,LUPRI)
      END IF
C
C     Compare with results using single vectors
C
      IOFF = KWRK1
      JOFF = KL2
      DO 500 I = 1, NCSIM + NOSIM
         IF (I .EQ. NCSIM + 1) JOFF = KL2 + NCONST*NVARPT
         CALL DAXPY(NVARPT,-D1,WRK(JOFF),1,WRK(IOFF),1)
         IOFF = IOFF + NVARPT
         JOFF = JOFF + NVARPT
  500 CONTINUE
      CALL AROUND('MC Hessian (sim. - single vectors) in ABALMT')
C
      CALL HEADER(
     &  'Vector    Largest element        CI and orbital norms      ',
     &   LUPRI)
      DNRMC = D0
      DNRMO = D0
      ZMAX  = D0
      IOFF  = KWRK1
      DO 600 I = 1, NCSIM + NOSIM
         IMAXEL = IDAMAX(NVARPT,WRK(IOFF),1)
         IF (NCONST .GT. 0) DNRMC = DNRM2 (NCONST,WRK(IOFF),1)
         IF (NWOPPT .GT. 0) DNRMO = DNRM2 (NWOPPT,WRK(IOFF+NCONST),1)
         IOFF = IOFF + NVARPT
         WRITE (LUPRI,'(5X,I6,5X,I8,9X,1P,2D16.8)') I,IMAXEL,DNRMC,DNRMO
         ZMAX = MAX(ZMAX,DNRMC,DNRMO)
  600 CONTINUE
      IF (IPRNR .GE. 5 .OR. ZMAX .GT. DTOL) THEN
         CALL OUTPUT(WRK(KWRK1),1,NVARPT,1,NCSIM+NOSIM,NVARPT,
     &               NCSIM+NOSIM,1,LUPRI)
      END IF
      RETURN
      END
