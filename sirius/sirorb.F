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
C  /* Deck sirorb */
      SUBROUTINE SIRORB(MAXITO,NITORB,CMO,GORB,DV,PV,FC,FV,FQ,
     &                  WRK,LFREE)
C
C Nov-1985 HJAaJ
C DFT modification tuh
C
C Purpose:
C  Carry out intermediate optimization of orbitals for fixed
C  CSF coefficients ("absorption")
C
C MOTECC-90: The purpose of this module, SIRORB, and the algorithms used
C            are described in Chapter 8 Section E.5 of MOTECC-90
C            "Intermediate Optimization of Orbitals for Fixed
C            Configuration Coefficients"
C
C
#include <implicit.h>
C
      DIMENSION CMO(*),GORB(*),DV(*),PV(*),FC(*),FV(*),FQ(*),
     &          WRK(*)
C
      PARAMETER ( TKFAC=0.1D0, LUDMMY=-1)
      PARAMETER ( MAXBCK=5, CVOFAC=0.8D0 )
      PARAMETER ( D0=0.0D0, D1=1.0D0, DP5=0.5D0 )
#include <dummy.h>
C
C  INFINP : FLAG(*),ISTATE,ITRLVL,?
C  INFTAP : ?
C  INFOPT : ITMAC, RTRUST
C  INFDIM : NWOPMA
C  INFORB : N2BAST,..
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inftap.h>
#include <infopt.h>
#include <infvar.h>
#include <inforb.h>
#include <infdim.h>
#include <infpri.h>
#include <inftra.h>
#include <dfterg.h>
C
      LOGICAL GETWOP, ACTROT, RHFCLC, NEWLVL
      CHARACTER*8 TABLE(4)
      DATA TABLE/'ABS1SYM1','ABS2SYM1','ABS3SYM1', 'EXOPSYM1'/
C
C
      CALL QENTER('SIRORB')
      WRITE (LUW4,'(//A/)')
     *    ' <<<< OUTPUT FROM ORBITAL ABSORPTION (SIRORB) >>>>'
      ACTROT   = FLAG(23)
      FLAG(23) = .TRUE.
C
C ALLOCATE WORK SPACE
C
      KDIAOR = 1
      KCMOOL = KDIAOR + NWOPMA
      KXKAP  = KCMOOL + NCMOT
      KWRK0  = KXKAP  + NWOPMA
      LWRK0  = LFREE  - KWRK0
C
      IF (LWRK0 .LT. 0) CALL ERRWRK('SIRORB',-KWRK0,LFREE)
C
C     Determine level of "absorption":
C     LVLABS = 1 : inac-act + act-act rotations
C            = 2 : + inac-sec rotations
C            = 3 : + act-sec rotations
C
      IF (FLAG(51)) THEN
         LVLABS = 1
!        JTRLVL = 0 ! Mar2011: TR1H2M sets NEEDMU(5) = 1 and we
!                     must therefore include secondary orbitals in the
!                     integral transformation although we will not need
!                     the integrals. There is no transparent way of
!                     transferring this information to TR1H2M, and as
!                     the CPU time saved will be a small fraction of
!                     the total CPU time, it is not worth doing more
!                     about. With this change, all levels of absorption
!                     should work again. /hjaaj
         JTRLVL = ITRLVL
      ELSE IF (FLAG(52)) THEN
         LVLABS = 2
!        JTRLVL = 0 ! see comment under FLAG(51) above.
         JTRLVL = ITRLVL
      ELSE IF (FLAG(53)) THEN
         LVLABS = 3
         JTRLVL = ITRLVL
      ELSE
         WRITE (LUW4,'(//A/)') '*** ERROR, SIRORB called but no'//
     *      ' absorption level specified.'
         CALL QTRACE(LUW4)
         CALL QUIT('*** ERROR (SIRORB) - SOFTWARE ERROR, PLEASE REPORT')
      END IF
#if defined (VAR_NOITRABS)
C --- 851216 -- EXPLICIT ORBITAL HESSIAN DIAGONAL IS WRONG
C               IF KTRLVL = 0 !
      JTRLVL   = ITRLVL
#endif
      WRITE (LUW4,'(/A,I3)') ' *** ORBITAL ABSORPTION LEVEL =',LVLABS
C
C     START-UP SECTION:
C     =================
C     GET ORBITAL ROTATION SPECIFICATION CORRESPONDING TO ABSORPTION
C     LEVEL USING LOGICAL FUNCTION GETWOP(XLABEL,LUNIT):
C
      IF ( .NOT. GETWOP(TABLE(LVLABS)) ) THEN
C        NWARN = NWARN + 1
C        WRITE (LUW4,'(//A,A8,A/A/)')
C    *      ' *** WARNING (SIRORB), label "',
C    *      TABLE(LVLABS),'" not found on file LUINDF',
C    *      '     Continuation without absorption will be attempted.'
         WRITE (LUW4,'(//A,A8,A/A/)') ' *** ERROR (SIRORB), label "',
     *      TABLE(LVLABS),'" not found on file LUINDF'
         CALL QTRACE(LUW4)
         CALL QUIT('*** ERROR (SIRORB) label not found on LUINDF')
      END IF
      WRITE (LUW4,'(/T6,A,I5)') 'Number of orbital rotations =',NWOPT
      NEWLVL = .FALSE.
C
      CALL UPDGRD (3,DUMMY,CMO,DV,PV,GORB,EMCMY,EMCACT,FC,FV,FQ,
     *             WRK(KDIAOR),DUMMY,WRK(KWRK0),LWRK0)
C     CALL UPDGRD (ICTL,CREF,CMO,DV,PV,G,EMCMY,EMCACT,FC,FV,FQ,
C    *             DIAOR,INDXCI,WRK,LFREE)
      MCONF = NCONF
      NCONF = 0
      CALL GRDINF(GNORM,GORB)
      NCONF = MCONF
C
C     Write orbital gradient which probably is greater than before
C     because of active-active rotations
C
      WRITE (LUW4,1020) GNORM(3)
      IF (LUPRI.NE.LUW4) WRITE (LUPRI,1020) GNORM(3)
C
      ITORB  = 0
      ITBCK  = 0
      ROTRST = RTRUST
      OSTP0  = D0
C
C
C **********************************************************************
C START OF NON-STANDARD LOOPS
C REPEAT OUTER LOOP UNTIL (converged or ITORB .gt. MAXITO or error)
C
 1000 CONTINUE
         ITORB = ITORB + 1
         WRITE (LUW4,21000) ITMAC,ITORB
         IF (LUPRI .NE. LUW4) WRITE (LUPRI,21000) ITMAC,ITORB
21000    FORMAT(//,' <<< ORBITAL ABSORPTION ITERATION NO. (',I3,',',
     *          I3,') >>>'/)
C
         CALL FLSHFO(LUW4)
         IF (LUPRI .NE. LUW4) CALL FLSHFO(LUPRI)
C
         EMCOLD = EMCSCF
         CALL DCOPY(NCMOT,CMO,1,WRK(KCMOOL),1)
         CALL DCOPY(NWOPT,GORB,1,WRK(KXKAP),1)
C
C        DAMP = initial beta value
C        Convergence in nexkap is tkfac of ORBITAL GRADIENT NORM.
C
         DAMP   = D1
         THRKAP = TKFAC * GNORM(3)
         EVAL   = D1
C        EVAL defined to avoid ftnchek error message, not needed here
         CALL NEXKAP(1,WRK(KXKAP),0,LUDMMY,LUDMMY,DUMMY,
     *               0,ROTRST,THRKAP,EVAL,DAMP,WRK(KDIAOR),
     *               CMO,GORB,DV,PV,FC,FV,WRK(KWRK0),LWRK0)
C        CALL NEXKAP(NGDM,GDM,NBFIL3,IFIL3,IFIL5,IBNDX,
C    *               NRSEQ,ROTRST,THRKAP,EVAL,DAMP,DIAOR,
C    *               CMO,GORB,DV,PV,FC,FV,WRK,LFREE)
         OSTP0  = DNORM2(NWOPT,WRK(KXKAP),1)
         XKAPN2 = OSTP0**2
         DEPRED = DP5*EVAL*(XKAPN2 + D1/(DAMP*DAMP))
         WRITE (LUW4,'(/A,F15.10)')
     *      ' Orbital rotation step length =',OSTP0
C
 1100    CONTINUE
            CALL MO_ROTATE(WRK(KXKAP),CMO,WRK(KWRK0),LWRK0)
C           CALL MO_ROTATE(XKAP,CMO,WRK,LFREE)
            CALL NEWORB('(ABSORP)',CMO,.FALSE.)
            FLAG(14) = FLAG(34)
C           ... new orbitals, integral transformation may be needed.
            IF ( ITORB .GE. MAXITO ) THEN
               WRITE (LUW4,'(//A/A/)')
     *           ' *** Orbital absorption discontinued,',
     *           '     maximum number of absorption iterations reached.'
               GO TO 9000
            ELSE IF ( OSTP0 .LE. CVOFAC*ROTRST ) THEN
               IF (FLAG(53) .AND. LVLABS .LT. 3) THEN
C                 Low level absorption is only useful until once
C                 converged (861016 - hjaaj)
                  FLAG(50+LVLABS) = .FALSE.
               END IF
 1200          LVLABS = LVLABS + 1
               IF (LVLABS .GT. 3) THEN
                  WRITE (LUW4,'(//A)')
     *               ' *** Orbital absorption converged.'
C                 Low level absorption is only useful in first
C                 macro iteration
                  FLAG(51) = .FALSE.
                  FLAG(52) = .FALSE.
C
                  GO TO 9000
               ELSE IF (.NOT. FLAG(50+LVLABS) ) THEN
                  GO TO 1200
C              ^------------
               END IF
               NEWLVL = .TRUE.
            END IF
            IF ( ITORB .EQ. (MAXITO - 1) .AND. LVLABS .LT. 3 ) THEN
C              always last absorption at highest level
               IF (FLAG(53)) THEN
                  LVLABS = 3
                  NEWLVL = .TRUE.
               ELSE IF (LVLABS .EQ. 1 .AND. FLAG(52)) THEN
                  LVLABS = 2
                  NEWLVL = .TRUE.
               END IF
            END IF
            IF (NEWLVL) THEN
               WRITE (LUW4,'(/A,I3)')
     *            ' *** NEW ORBITAL ABSORPTION LEVEL =',LVLABS
               IF ( .NOT. GETWOP(TABLE(LVLABS)) ) THEN
                  WRITE (LUW4,'(//A,A8,A/)')
     *               ' *** ERROR (SIRORB), label "',TABLE(LVLABS),
     *               '" not found on file LUINDF'
                  CALL QTRACE(LUW4)
                  CALL QUIT('*ERROR (SIRORB) label not found on LUINDF')
               END IF
               WRITE (LUW4,'(/T6,A,I5)')
     *            'Number of orbital rotations =',NWOPT
               IF (LVLABS .EQ. 3) JTRLVL = ITRLVL
               ROTRST = RTRUST
               NEWLVL = .FALSE.
            END IF
C
            CALL FLSHFO(LUW4)
            IF (LUPRI .NE. LUW4) CALL FLSHFO(LUPRI)
C
C
C           Transform integrals
C
            IF (.NOT. FLAG(14)) THEN
               CALL TRACTL(JTRLVL,CMO,WRK(KWRK0),LWRK0)
C              CALL TRACTL(ITRLVL,CMO,WRK,LFREE)
C
            END IF
C
            CALL UPDGRD (3,DUMMY,CMO,DV,PV,GORB,EMCMY,EMCACT,FC,FV,FQ,
     *                   WRK(KDIAOR),DUMMY,WRK(KWRK0),LWRK0)
C           CALL UPDGRD (ICTL,CREF,CMO,DV,PV,G,EMCMY,EMCACT,FC,FV,FQ,
C    *                   DIAOR,INDXCI,WRK,LFREE)
            MCONF  = NCONF
            NCONF  = 0
            CALL GRDINF(GNORM,GORB)
            NCONF  = MCONF
            EMCSCF = POTNUC + EMCMY + EMCACT + EDFTY
            WRITE (LUW4,1010) EMCSCF,ITORB
            IF ( P4FLAG(2) .OR. (LUPRI.EQ.LUW4 .AND. P6FLAG(1)) )
     *         WRITE (LUW4,1015) POTNUC,EMCMY,EMCACT
            WRITE (LUW4,1020) GNORM(3)
            IF (LUPRI.NE.LUW4) THEN
               WRITE (LUPRI,1010) EMCSCF,ITORB
               IF (P6FLAG(1))WRITE (LUPRI,1015) POTNUC,EMCMY,EMCACT
               WRITE (LUPRI,1020) GNORM(3)
            END IF
 1010       FORMAT(//' Total MCSCF energy       :',F25.15,T54,
     *               '(after absorp ',I3,')' )
 1015       FORMAT(//' - Nuclear repulsion      :',F25.15,
     *              /' - Inactive energy        :',F25.15,
     *              /' - Active energy          :',F25.15)
 1020       FORMAT(//' Norm of orbital gradient :',F22.12)
C
C           PRINT GRADIENT, IF REQUESTED
C
            IF (P4FLAG(14)) THEN
               WRITE (LUW4,3020) ITORB
               IF (MPRI4.GT.100) THEN
                  PRFAC = 0.0D0
               ELSE IF (MPRI4.GT.10) THEN
                  PRFAC = 0.1D0
               ELSE
                  PRFAC = 0.2D0
               END IF
               CALL PRKAP (NWOPT,GORB,PRFAC,LUW4)
            END IF
 3020       FORMAT(/' Orbital gradient after absorption it.no.',I3,
     *             /' -------------------------------------------')
C
            CALL FLSHFO(LUW4)
            IF (LUPRI.NE.LUW4) CALL FLSHFO(LUPRI)
C
C
C           BACKSTEP ?
C
            CALL ORBSTP(ISTEP,DAMP,EVAL,OSTP0,CMO,WRK(KCMOOL),
     *                  WRK(KXKAP),ROTRST)
C           CALL ORBSTP(ISTEP,DAMP,EVAL,OSTP0,CMO,CMOOL,XKAP,ROTRST)
            IF (ISTEP.NE.0) THEN
               ITBCK = ITBCK + 1
               WRITE (LUW4,21100) ITBCK
               IF (LUPRI .NE. LUW4) WRITE (LUPRI,21100) ITBCK
               IF (ITBCK .GT. MAXBCK) THEN
                  WRITE (LUW4,'(//A/A)')
     *               ' *** Orbital absorption discontinued,'//
     *               ' maximum number of absorption backsteps reached.',
     *               ' Old orbitals restored.'
                  CALL DCOPY(NCMOT,WRK(KCMOOL),1,CMO,1)
                  CALL NEWORB('ABANDABS',CMO,.FALSE.)
                  FLAG(14) = FLAG(34)
C                 ... integral transformation may be needed.
                  GO TO 9000
               END IF
               GO TO 1100
C        ^--------------- backstep loop
            ELSE
               ITBCK = 0
               GO TO 1000
C     ^------------------ orbital iteration loop
            END IF
21100       FORMAT(//' *** BACK STEP -- THIS ORBITAL ABSORPTION STEP ',
     *         'IS REJECTED ***',
     *        /'  -  BACKUP NUMBER',I3,' IN THIS ORBITAL ABSORPTION.')
C END OF NON-STANDARD LOOPS
C **********************************************************************
C Recover original orbital rotation specifications and return
C
 9000 CONTINUE
      FLAG(23) = ACTROT
      IF ( .NOT. GETWOP(TABLE(4)) ) THEN
         WRITE (LUW4,'(//A,A8,A/)')' *** ERROR (SIRORB), label "',
     *      TABLE(4),'" not found on file LUINDF'
         CALL QTRACE(LUW4)
         CALL QUIT('*** ERROR (SIRORB) label not found on LUINDF')
      END IF
      NITORB = ITORB
C
      CALL FLSHFO(LUW4)
      IF (LUPRI .NE. LUW4) CALL FLSHFO(LUPRI)
C
      CALL QEXIT('SIRORB')
      RETURN
C
C END OF SIRORB.
C
      END
C  /* Deck orbstp */
      SUBROUTINE ORBSTP(ISTEP,DAMP,EVAL,OSTP0,CMO,CMOOLD,XKAP,ROTRST)
C
C
C Purpose:
C     STEP CONTROL -- check if step in intermediate orbital
C     optimization is too large, that is if actual energy
C     change deviates too much from the second order prediction in
C     last absorption iteration.
C     When step is rejected the new predicted second order energy
C     is returned in DEPRED , scaled step is returned in
C     kappa and molecular orbital expansion coefficients at previous
C     expansion point is in CMO. OSTP0 is orbital  step length which is
C     returned from nexkap
C     Signal to calling program with ISTEP = 0 if step is OK,
C     with ISTEP = 1 if step is to large.
C
C Output:
C  ISTEP, =0 if step is OK
C         =1 if step is too large
C
#include <implicit.h>
C
      DIMENSION CMO(*),CMOOLD(*),XKAP(*)
C
      PARAMETER (DP1=0.1D0, DP5=0.5D0, D1=1.0D0, D2=2.0D0)
      PARAMETER (THDE = 1.D-10)
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <infpri.h>
#include <inforb.h>
#include <infvar.h>
#include <infopt.h>
C
C
      CALL QENTER('ORBSTP')
C
C
C *****************************************************************
C *** Trust region control
C
      RTSAVE= ROTRST
      ISTEP = 0
      DEACT = EMCSCF - EMCOLD
      IF (ABS(DEPRED).GT.THDE) THEN
         RATIO = DEACT / DEPRED
         WRITE (LUW4,1100) DEACT,DEPRED,RATIO
      ELSE
         RATIO = D1
         WRITE (LUW4,1110) DEACT,DEPRED
      END IF
 1100 FORMAT (/' (ORBSTP) Energy difference;',
     *         ' actual, predicted and ratio:',/5X,3F15.10)
 1110 FORMAT (/' (ORBSTP) Close to convergence, ratio set to one.',
     *        /' Energy difference; actual and predicted:',
     *        /5X,1P,2D15.5)
C
C
      IF (ISTATE .EQ. 1) THEN
C     (case 1: ground state optimization)
         IF (RATIO.LT.RATMIN) THEN
            ROTRST = STPRED*ROTRST
            IF (RATIO.LT.RATREJ) ISTEP = 1
         ELSE IF (RATIO.GT.RATGOD) THEN
            ROTRST = MIN(STPMAX,STPINC*ROTRST)
         END IF
      ELSE
C     (case 2: excited state optimization)
         IF (RATIO.LT.RATMIN .OR. RATIO.GT.(D2-RATMIN)) THEN
            ROTRST = STPRED*ROTRST
            IF (RATIO.LT.RATREJ .OR. RATIO.GT.(D2-RATREJ)) ISTEP = 1
         ELSE IF (RATIO.GT.RATGOD .AND. RATIO.LT.(D2-RATGOD)) THEN
            ROTRST = MIN(STPMAX,STPINC*ROTRST)
         END IF
      END IF
C
C *****************************************************************
C *** Step is acceptable:
C
      IF (ISTEP.EQ.0) THEN
         GO TO 9000
      END IF
C
C *****************************************************************
C *** Step is too large:
C
      WRITE (LUW4,'(/A,F10.5)')
     *   ' (ORBSTP) step is too large -- STEP IS SCALED TO',ROTRST
C
C     PUT OLD MOLECULAR ORBITALS IN CMO
C
      CALL DCOPY(NCMOT,CMOOLD,1,CMO,1)
C
C     SCALE STEP TO ROTRST
C
      OSTP1  = DNORM2(NWOPT,XKAP,1)
      FAC    = ROTRST/OSTP1
      CALL DSCAL(NWOPT,FAC,XKAP,1)
C
C     FIND NEW PREDICTED ENERGY
C
      FAC    = ROTRST / OSTP0
      DEPRED = DP5*EVAL*FAC*FAC*( (D2/FAC-D1)/(DAMP*DAMP) + OSTP0*OSTP0)
      OSTP0  = OSTP1
C
C *****************************************************************
C *** End of subroutine ORBSTP
C
 9000 CONTINUE
      WRITE (LUW4,'(/A,2F15.10)')
     *   ' (ORBSTP) Old and new trust radius:',RTSAVE,ROTRST
      CALL QEXIT('ORBSTP')
      RETURN
      END
C  /* Deck orbabs */
      LOGICAL FUNCTION ORBABS(ISTATE,IBNDX,WRK)
C
C  CHECK IF THE CSF PART OF THE HESSIAN CONTAINS THE DESIRED NUMBER
C  OF NEGATIVE EIGENVALUES. ORBABS IS THEN .TRUE. ELSE .FALSE.
C
#include <implicit.h>
      DIMENSION IBNDX(*),WRK(*)
C
#include <priunit.h>
      PARAMETER ( D0=0.0D0 )
C
#include <inftap.h>
C
      LOGICAL FNDLAB
C
C Recover IBNDX and REDL from previous macro iteration
C
      REWIND LUIT1
      IF (.NOT.FNDLAB('RESTART ',LUIT1)) GO TO 8000
      READ (LUIT1) DUM,DUM,DUM,DUM,DUM,IDUM,
     *             MREDL,(IBNDX(I),I=1,MREDL)
      IF (.NOT.FNDLAB('LREDUCED',LUIT1)) GO TO 8000
      NMREDL = MREDL*(MREDL+1)/2
      CALL READT(LUIT1,NMREDL,WRK)
C
C ALLOCATE WORK SPACE
C
      KREDL = 1
      KREDCI= KREDL + NMREDL
      KEVEC = KREDCI+ NMREDL
      KWRK1 = KEVEC + MREDL*MREDL
      KIWRK1= KWRK1 + MREDL
C
C RECOVER CSF PART OF REDUCED HESSIAN
C
      KCIDIM=0
      KIJ   = KREDCI -1
      IJ    = KREDL  -1
      DO 100 I=1,MREDL
         IF (I.NE.1) THEN
            IF (IBNDX(I).EQ.-1) KCIDIM = KCIDIM + 1
         END IF
         DO 110 J=1,I
            IJ = IJ + 1
            IF (IBNDX(I).NE.-1) GO TO 110
            IF ((J.NE.1) .AND. (IBNDX(J).EQ.-1))THEN
               KIJ = KIJ + 1
               WRK(KIJ) = WRK(IJ)
            END IF
 110     CONTINUE
 100  CONTINUE
C
C DIAGONALIZE CSF PART OF REDUCED HESSIAN
C
      CALL DUNIT(WRK(KEVEC),KCIDIM)
      CALL JACO(WRK(KREDCI),WRK(KEVEC),KCIDIM,KCIDIM,KCIDIM,
     *          WRK(KWRK1),WRK(KIWRK1))
C
C COUNT NUMBER OF NEGATIVE EIGENVALUES IN CSF PART OF REDUCED HESSIAN
C
      NNEG = 0
      II   = KREDCI - 1
      DO 150 I=1,KCIDIM
         II = II + I
         IF (WRK(II).LT.D0) NNEG = NNEG + 1
 150  CONTINUE
      IF ( ISTATE .EQ. NNEG ) THEN
         ORBABS = .TRUE.
      ELSE
         ORBABS = .FALSE.
      END IF
 1000 CONTINUE
      RETURN
C
 8000 CONTINUE
         ORBABS = .FALSE.
         WRITE (LUPRI,'(/A,/A)')
     *         ' Failed to read in reduced Hessian in ORBABS.',
     *         ' Orbital absorption will not be performed.'
         GO TO 1000
C     END OF ORBABS.
      END
C  /* Deck rotflp */
      FUNCTION ROTFLP(ISTATE,IBNDX,NREDL,REDL,WRK)
C
C  17-May-1986 PJ
C
C  IF DESIRED NUMBER OF NEGATIVE EIGENVALUES IN THE
C  HESSIAN IS IN THE CSF SPACE, THEN ROTFLP IS FALSE,
C  ELSE ROTFLP IS TRUE
C
#include <implicit.h>
      LOGICAL ROTFLP
      DIMENSION IBNDX(*),REDL(*),WRK(*)
      PARAMETER ( TOLFLP = 0.75D0 )
#include <ibndxdef.h>
C
C INFPRI: IPRI6
C
#include <priunit.h>
#include <infpri.h>
C
C ALLOCATE WORK SPACE
C
      KREDC = 1
      KEVEC = KREDC + NREDL*(NREDL+1)/2
      KBIG  = KEVEC + NREDL*NREDL
      KJBIG = KBIG  + NREDL
C
      K = 0
      DO 100 I=1,NREDL
         IF (IBNDX(I).EQ.JBCNDX) THEN
            IOFF = I * (I-1) / 2
            K = K + 1
            KOFF = K * (K-1) /2
            DO 200 J=1,I
               IF (IBNDX(J).EQ.JBCNDX) THEN
                  KOFF = KOFF + 1
                  WRK(KOFF) = REDL(IOFF+J)
               END IF
 200        CONTINUE
         END IF
 100  CONTINUE
      CALL DUNIT(WRK(KEVEC),K)
C
      CALL JACO(WRK,WRK(KEVEC),K,K,K,WRK(KBIG),WRK(KJBIG))
      II = 0
      DO 300 I=1,K
         II = II + I
         WRK(I) = WRK( II )
 300  CONTINUE
      CALL ORDER(WRK(KEVEC),WRK,K,K)
      IF (WRK(KEVEC+(ISTATE-1)*K ) .GT. TOLFLP) THEN
         ROTFLP=.FALSE.
      ELSE
         ROTFLP=.TRUE.
      END IF
C
      IF (IPRI6 .GT. 5) THEN
         WRITE (LUPRI,'(//A,L8)')
     *      ' *** CHECK OF ROOT FLIPPING; ROTFLP =',ROTFLP
      END IF
      IF (IPRI6 .GT. 15) THEN
         WRITE (LUPRI,'(/A,I5//A/,(5F15.8))')
     *      '     Desired state is no.',ISTATE,
     *      '     Eigenvalues of reduced CI matrix :',
     *      (WRK(I), I = 1,K)
         WRITE (LUPRI,'(/A,/A,F12.8,A)')
     *      '     Eigenvectors of reduced CI matrix :',
     *      '     (tolerance for root flip =',TOLFLP,
     *      ' on basis vector no. 1)'
         CALL OUTPUT(WRK(KEVEC),1,K,1,K,K,K,1,LUPRI)
      END IF
      RETURN
C     END OF ROTFLP
      END
