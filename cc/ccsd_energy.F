C  /* Deck ccsd_energy */
      SUBROUTINE CCSD_ENERGY(LSORT,WORK,LWORK)
C
C     Written by Henrik Koch 27-Mar-1990.
C     DIIS and Brueckner bit by Rika Kobayashi 1992.
C
C     Ove juli-sept. 1995: RSP intermediates
C                          noccit
C     Ove februar    1997: CCS, FD gradient hacks and restart.
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxorb.h>
      PARAMETER (XMONE = -1.0D0, IZERO = 0, TWO = 2.0D0, ZERO = 0.0D00)
      LOGICAL   CCSAV, CC1BSV, CC1ASV, CCPTSV, CCP3SV, MP2DIF, LSORT
      DIMENSION WORK(LWORK)
      COMMON /LUDIIS/ LUTDIS, LUSDIS
#include <ccorb.h>
#include <iratdef.h>
#include <ccsdinp.h>
#include <ccsections.h>
#include <ccsdsym.h>
#include <ccfro.h>
#include <ccsdio.h>
#include <ccinftap.h>
#include <inftap.h>
#include <cclr.h>
#include <ccslvinf.h>
#include <ccfdgeo.h>
      LOGICAL LCONVG,RSPIM2,EX
      CHARACTER*5 ETY0, ETY1, ETY2, MODEL*10
C
#if defined (SYS_CRAY)
C     Open file for diis extrapolation
C
      CALL WOPEN('CC_DIIS',64,0,IERR)
C
      IF (IERR .NE. 0) CALL QUIT('Error opening CC_DIIS')
#endif
C
      ETY0 = 'SCF  '
C
C     Call the CCSD initialization routine.
C
      ISYMOP = 1
C
      RSPIM2 = .FALSE.
      OMEGSQ = .FALSE.
      OMEGOR = .TRUE.
      DUMPCD = .TRUE.
      CC3LR  = .FALSE.
      NEWGAM = .TRUE.
      CCPTSV = .FALSE.
      CCP3SV = .FALSE.
      ECURR  = ZERO
      EX     = .FALSE.
C
      IALPR  = IPRINT
C
      IF (INTTR) CALL CCSD_INIT1(WORK,LWORK)
C
C----------------------------------------------------------------------
C     Save RSPIM flag to calculate response global intermediates later.
C     If CCS or MP2 no intermediates is calculated.
C----------------------------------------------------------------------
C
      IF (RSPIM .AND. ( .NOT. (CCS .OR.(MP2.AND.(.NOT.CCP2))))) THEN
         RSPIM2 = RSPIM
         RSPIM  = .FALSE.
      ENDIF
C
C------------------------------
C     Print information header.
C------------------------------
C
      WRITE (LUPRI,'(1x,A,/)') '  '
      WRITE (LUPRI,'(1x,A)')
     *'*********************************************************'//
     *'**********'
      WRITE (LUPRI,'(1x,A)')
     *'*                                                        '//
     *'         *'
      WRITE (LUPRI,'(1x,A)')
     *'*<<<<<<<<<<                                             >'//
     *'>>>>>>>>>*'
      WRITE (LUPRI,'(1x,A)')
     *'*<<<<<<<<<< OUTPUT FROM COUPLED CLUSTER ENERGY PROGRAM  >'//
     *'>>>>>>>>>*'
      WRITE (LUPRI,'(1x,A)')
     *'*<<<<<<<<<<                                             >'//
     *'>>>>>>>>>*'
      WRITE (LUPRI,'(1x,A)')
     *'*                                                        '//
     *'         *'
      WRITE (LUPRI,'(1x,A,/)')
     *'*********************************************************'//
     *'**********'
C
      WRITE(LUPRI,*)
      WRITE(LUPRI,'(13X,A)')
     *     'The Direct Coupled Cluster Energy Program'
      WRITE(LUPRI,'(13X,A)')
     *     '-----------------------------------------'
      WRITE(LUPRI,'(//10X,A,I8)')
     *     'Number of t1 amplitudes                 :  ',NT1AMX
      WRITE(LUPRI,'(10X,A,I8)')
     *     'Number of t2 amplitudes                 :  ',NT2AMX
      WRITE(LUPRI,'(10X,A,I8)')
     *     'Total number of amplitudes in ccsd      :  ',
     *     (NT1AMX+NT2AMX)
      WRITE(LUPRI,*)
C
C----------------------------------------------------
C     IF iterative triples, open files for integrals.
C----------------------------------------------------
C
      IF (CCSDT.OR.CCPT.OR.CCP3.OR.(CCR3.OR.CCRT.OR.CCR1A.OR.CCR1B))
     *    CALL CC3_FILOP()
C
C------------------------------------------
C     Sort AO-integrals into distributions.
C------------------------------------------
C
      IF ((.NOT. DIRECT).AND. LSORT) THEN
         CALL CCSD_SORTAO(WORK,LWORK,IALPR)
      ENDIF
C
C----------------------------------------------------------------
C     If CCS then no the wavefunction optimization.
C     CCS energy is equal to HF energy -> find and put in ECCGRS.
C     For polarizabilities and oscillatro strengths,
C     we need the (ia|jb) integrals.
C----------------------------------------------------------------
C
      IF (CCS ) THEN
         WRITE(LUPRI,'(//10X,A,I8)')
     *                'CCS CALC. - NO WAVEFUNCTION OPTIMIZATION'
         IF (INTTR.AND.(.NOT.LISKIP)) THEN
            KT1AM = 1
            KT2AM = KT1AM + NT1AMX
            KEND1 = KT2AM + NT2AMX
            LWRK1 = LWORK - KEND1
            IF ( KEND1 .GT. LWORK ) THEN
               CALL QUIT('Insufficient spaces in CCSD_ENERGY')
            ENDIF
            CALL DZERO(WORK(KT1AM),NT1AMX)
            CALL CCSD_IAJB(WORK(KT2AM),WORK(KT1AM),WORK(KEND1),LWRK1)
            REWIND(LUIAJB)
            CALL WRITI(LUIAJB,IRAT*NT2AM(ISYMOP),WORK(KT2AM))
         ENDIF
C
         CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND LUSIFC
         CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
         READ (LUSIFC) POTNUC,EMY,EACTIV,EMCSCF,ISTATE,ISPIN,NACTEL,
     *              LSYM
         ESCF   = EMCSCF
         ECCGRS = EMCSCF
         CALL GPCLOSE(LUSIFC,'KEEP')
C
C        write SCF energy to summary file:
         WRITE(LURES,'(/12X,A,A,A,F32.10)')
     *               'Total ',ETY0,' energy: ',ESCF
C
         RETURN
      ENDIF
C
C------------------------
C     Dynamic allocation.
C------------------------
C
      KFOCKD  = 1
      KT1AM   = KFOCKD  + NORBTS
      KOMEG1  = KT1AM   + NT1AMX
      KOMEG2  = KOMEG1  + NT1AM(ISYMOP)
      IF (OMEGSQ) THEN
         KT2AM = KOMEG2
     *         + MAX(NT2AMX,NT2AM(ISYMOP),NT2AO(ISYMOP),NT2AOS(ISYMOP))
      ELSE
         KT2AM = KOMEG2
     *        + MAX(NT2AMX,NT2AM(ISYMOP),NT2AO(ISYMOP),2*NT2ORT(ISYMOP))
      ENDIF
      KSCR2   = KT2AM   + NT2AMX
      KEND1   = KSCR2   + NT2AMX + NT1AMX
      LWRK1   = LWORK   - KEND1
C
      IF ( KEND1 .GT. LWORK ) THEN
         CALL QUIT('Insufficient spaces in CCSD_ENERGY')
      ENDIF
C
C-------------------------------------
C     Read canonical orbital energies.
C-------------------------------------
C
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
C
      CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
      READ (LUSIFC) POTNUC,EMY,EACTIV,EMCSCF,ISTATE,ISPIN,NACTEL,
     *              LSYM
C
      ESCF = EMCSCF
C
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KFOCKD+I-1), I=1,NORBTS)
C
      CALL GPCLOSE(LUSIFC,'KEEP')
C
C----------------------------------------------------------
C     Change symmetry-ordering of the Fock-matrix diagonal.
C----------------------------------------------------------
C
      IF (FROIMP .OR. FROEXP)
     *    CALL CCSD_DELFRO(WORK(KFOCKD),WORK(KEND1),LWRK1)
C
      CALL FOCK_REORDER(WORK(KFOCKD),WORK(KEND1),LWRK1)
C
C-----------------------------------------------------------
C     Calculate the ( ia | jb ) integrals and write to disk.
C-----------------------------------------------------------
C
      IF (INTTR.AND.(.NOT.LISKIP)) THEN
         CALL DZERO(WORK(KT1AM),NT1AMX)
         CALL CCSD_IAJB(WORK(KT2AM),WORK(KT1AM),WORK(KEND1),LWRK1)
         REWIND(LUIAJB)
         CALL WRITI(LUIAJB,IRAT*NT2AM(ISYMOP),WORK(KT2AM))
      ELSE
         CALL DCOPY(NT2AM(ISYMOP),99.99D0,0,WORK(KT2AM),1)
         REWIND(LUIAJB)
         CALL READI(LUIAJB,IRAT*NT2AM(ISYMOP),WORK(KT2AM))
      ENDIF
C
C------------------------------------
C     Setup the initial guess vector.
C------------------------------------
C
      IF (CCRSTR.AND.(.NOT.MP2)) THEN
         ILSTNR = 1
         ISYM   = 1
         IOPT   = 33
         CALL CC_RDRSP('R0',ILSTNR,ISYM,IOPT,MODEL,WORK(KT1AM),
     *                 WORK(KT2AM))
         ETY1 = 'RSTAR'
      ENDIF
   
      IF  ( (.NOT.CCRSTR) .OR. MP2 .OR. (IOPT.EQ.33) ) THEN
         CALL CCSD_GUESS(WORK(KT1AM),WORK(KT2AM),WORK(KFOCKD),
     *                   IALPR)
C
C-----------------------------------------------------------------
C        Print largest amplitudes in the zero order wave function.
C-----------------------------------------------------------------
C
         IF (IPRINT .GT. 4) THEN
C
            CALL AROUND('Largest amplitudes in MP2 guess')
C
            CALL DCOPY(NT1AMX,WORK(KT1AM),1,WORK(KOMEG1),1)
            CALL DCOPY(NT2AMX,WORK(KT2AM),1,WORK(KOMEG1+NT1AMX),1)
C
            CALL CC_PRAM(WORK(KOMEG1),PT1,1)
C
         ENDIF
C
         ETY1 = 'MP2  '
C
      ENDIF
C
C
      EN1=0D0
      EN2=99D0
      LCONVG=.FALSE.
      ITER=1
C
C     START OF ITERATIVE LOOP
C
C
      IF (LCOR .OR. LSEC ) THEN
C
         CALL CC_CORE(WORK(KT1AM),WORK(KT2AM),1)
C
      ENDIF
C
      IF (CCD) THEN
         CALL DZERO(WORK(KT1AM),NT1AMX)
      ENDIF
      IF (CCSTST) THEN
         CALL DZERO(WORK(KT2AM),NT2AMX)
      ENDIF
C
      IT1 = 0
      IF ( CCRSTR ) IT1 = 1
      CALL CCSD_ECCSD(WORK(KT1AM),WORK(KT2AM),WORK(KFOCKD),
     *                WORK(KEND1),LWRK1,EN2,POTNUC,ESCF,
     *                ETY1,IALPR,IT1,ITER)
C
      EINI = EN2
C
C----------------------------------------
C     If MP2 or NOCCIT do not enter loop.
C----------------------------------------
C
      IF (MP2 .OR. NOCCIT) GOTO 500
C
      IF (CCPT .OR. CCP3) THEN
         CCSAV = CCSDT
         CCSDT = .FALSE.
      ENDIF
C
  200 CONTINUE
C
C---------------------------------
C        Write amplitudes to disk.
C---------------------------------
C
C        needed in loop over models if CCDERI=.true. !?
C        IF (LUTAM.GT.0) CALL GPCLOSE(LUTAM,'KEEP')
C
         LUTAM = 0
         CALL GPOPEN(LUTAM,'CCSD_TAM','UNKNOWN',' ','UNFORMATTED',
     *               IDUMMY,.FALSE.)
         REWIND (LUTAM)
         WRITE(LUTAM) (WORK(KT1AM+I-1), I = 1,NT1AMX)
         WRITE(LUTAM) (WORK(KT2AM+I-1), I = 1,NT2AMX)
         CALL GPCLOSE(LUTAM,'KEEP')
C
         EN1 = EN2
C
         IF ( IPRINT .GT. 2 ) THEN
            WRITE(LUPRI,249) ITER
  249       FORMAT(/,3X,'    Iteration no.:',I3)
            WRITE(LUPRI,*)'      -----------------'
            WRITE(LUPRI,*)
         ENDIF
         CALL CCRHSN(WORK(KOMEG1),WORK(KOMEG2),WORK(KT1AM),WORK(KT2AM),
     *               WORK(KEND1),LWRK1)
C
  240    CONTINUE
C
         IF (LCOR .OR. LSEC ) THEN
            CALL CC_CORE(WORK(KOMEG1),WORK(KOMEG2),1)
         ENDIF
         IF (CCSTST) THEN
            CALL DZERO(WORK(KOMEG2),NT2AMX)
         ENDIF
         IF (CCD) THEN
            CALL DZERO(WORK(KOMEG1),NT1AMX)
         ENDIF
C
         IF (IPRINT .GE. 5) THEN
            WRITE(LUPRI,529) 'Norm of t1am   after ccvec:',
     *               DDOT(NT1AMX,WORK(KT1AM),1,WORK(KT1AM),1)
            WRITE(LUPRI,529) 'Norm of t2am   after ccvec:',
     *               DDOT(NT2AMX,WORK(KT2AM),1,WORK(KT2AM),1)
         ENDIF
         IF (IPRINT .GE. 3) THEN
            WRITE(LUPRI,529) 'Norm of omega1 after ccvec:',
     *               DDOT(NT1AM(ISYMOP),WORK(KOMEG1),1,WORK(KOMEG1),1)
            WRITE(LUPRI,529) 'Norm of omega2 after ccvec:',
     *               DDOT(NT2AM(ISYMOP),WORK(KOMEG2),1,WORK(KOMEG2),1)
         END IF
C
         IF (CCSLV) THEN
           IF (IPRINT .GE. 3) THEN
              WRITE(LUPRI,529) 'Norm of omega1 in sol. part.:',
     *              DDOT(NT1AM(ISYMOP),WORK(KOMEG1),1,WORK(KOMEG1),1)
              WRITE(LUPRI,529) 'Norm of omega2 in sol. part.:',
     *              DDOT(NT2AM(ISYMOP),WORK(KOMEG2),1,WORK(KOMEG2),1)
           END IF
           LUSLV = 0
           CALL GPOPEN(LUSLV,'CC_OME','UNKNOWN',' ','UNFORMATTED',
     *                 IDUMMY,.FALSE.)
           REWIND (LUSLV)
           WRITE(LUSLV) (WORK(KOMEG1+I-1), I = 1,NT1AMX)
           WRITE(LUSLV) (WORK(KOMEG2+I-1), I = 1,NT2AMX)
           CALL GPCLOSE(LUSLV,'KEEP')
         ENDIF
C
  529    FORMAT(7X,A,D25.10)
C
         CALL CCSD_NXTAM(WORK(KT1AM),WORK(KT2AM),WORK(KOMEG1),
     *                   WORK(KOMEG2),WORK(KFOCKD))
C
         IF (CCSTST) THEN
            CALL DZERO(WORK(KT2AM),NT2AMX)
         ENDIF
         IF (CCD) THEN
            CALL CCSD_DIIS(WORK(KOMEG2),WORK(KT2AM),NT2AMX,ITER)
         ELSE
	    CALL DCOPY(NT2AMX,WORK(KT2AM),1,WORK(KSCR2+NT1AMX),1)
            CALL DCOPY(NT1AMX,WORK(KT1AM),1,WORK(KT2AM),1)
            CALL DCOPY(NT2AMX,WORK(KSCR2+NT1AMX),1,WORK(KT2AM+NT1AMX),1)
            CALL CCSD_DIIS(WORK(KOMEG1),WORK(KT2AM),NT2AMX+NT1AMX,ITER)
            CALL DCOPY(NT1AMX,WORK(KT2AM),1,WORK(KT1AM),1)
            CALL DCOPY(NT2AMX,WORK(KT2AM+NT1AMX),1,WORK(KSCR2+NT1AMX),1)
            CALL DCOPY(NT2AMX,WORK(KSCR2+NT1AMX),1,WORK(KT2AM),1)
         ENDIF
C
C---------------------------------
C        The order is important !!
C---------------------------------
C
         ETY2 = 'CCSD '
         IF (CCD)   ETY2 = 'CCD  '
         IF (CCSDT) ETY2 = 'CC3  '
         IF (CC1B)  ETY2 = 'CC-1b '
         IF (CC1A)  ETY2 = 'CC-1a '
         IF (CC2)   ETY2 = 'CC2  '
         IF (MCC2)  ETY2 = 'MCC2'
C
         IT1 = 1
         CALL CCSD_ECCSD(WORK(KT1AM),WORK(KT2AM),WORK(KFOCKD),
     *                   WORK(KEND1),LWRK1,EN2,POTNUC,ESCF,
     *                   ETY2,IPRINT,IT1,ITER)
C
         CALL FLSHFO(LUPRI)
C
         IF (DABS(EN2-EN1) .LT. THRENR) LCONVG=.TRUE.
	 ITER = ITER+1
	 IF (ITER .GT. MAXITE) THEN
	    WRITE(LUPRI,*) 'Energy not converged in ',MAXITE,
     &           ' iterations'
	    CALL QUIT(' ')
         ENDIF
      IF (.NOT. LCONVG) GOTO 200
C
      WRITE(LUPRI,250) ETY2,THRENR,EN2
 250  FORMAT(/,1x,A5,' energy converged to within ',D10.2,' is ',F16.9)
C
C-----------------
C     end of loop.
C-----------------
C
#if defined (SYS_CRAY)
      CALL WCLOSE('CC_DIIS',IERR)
      INFO = ISHELL('rm CC_DIIS')
#endif
#if !defined (SYS_CRAY)
      CALL GPCLOSE(LUTDIS,'DELETE')
      CALL GPCLOSE(LUSDIS,'DELETE')
#endif
C
C--------------------------------------------------------------
C     Print largest amplitudes in the zero order wave function.
C--------------------------------------------------------------
C
      IF (IPRINT .GT. 2) THEN
C
         CALL AROUND('Largest amplitudes in converged solution')
C
         CALL DCOPY(NT1AMX,WORK(KT1AM),1,WORK(KOMEG1),1)
         CALL DCOPY(NT2AMX,WORK(KT2AM),1,WORK(KOMEG1+NT1AMX),1)
C
         CALL CC_PRAM(WORK(KOMEG1),PT1,1)
C
      ENDIF
C
 500  CONTINUE
C
C------------------------------------
C     Write final amplitudes to disk.
C------------------------------------
C
      IF (CCSTST) THEN
         CALL DZERO(WORK(KT2AM),NT2AMX)
      ENDIF
      LUTAM = 0
      CALL GPOPEN(LUTAM,'CCSD_TAM','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     *            .FALSE.)
      REWIND (LUTAM)
      WRITE(LUTAM) (WORK(KT1AM+I-1), I = 1,NT1AMX)
      WRITE(LUTAM) (WORK(KT2AM+I-1), I = 1,NT2AMX)
      CALL GPCLOSE(LUTAM,'KEEP')
C
C----------------------------------
C     save a copy on file CCR0___0
C----------------------------------
C
      KT0AM = KT2AM + NT2AMX
      KEND0 = KT0AM + NALLAI(1)
      LRWK0 = LWORK - KEND0

      CALL DZERO(WORK(KT0AM),2*NALLAI(1))
C
C     -------------------------------------------------------------
C     set model for which the current t-amplitudes were calculated:
C     (this should be more or less the same as ETY2)
C     -------------------------------------------------------------
C
      MODEL = 'UNKOWN'
      IF (CIS)  MODEL = 'CIS'
      IF (CCS)  MODEL = 'CCS'
      IF (MP2)  MODEL = 'MP2'
      IF (CC2)  MODEL = 'CC2'
      IF (CCD)  MODEL = 'CCD'
      IF (CCSD) MODEL = 'CCSD'
      IF (CC3)  MODEL = 'CC3'
      IF (CC1A) MODEL = 'CCSDT-1a'
      IF (CC1B) MODEL = 'CCSDT-1b'
      IF (CCPT) MODEL = 'CCSD'
      IF (CCP3) MODEL = 'CCSD'
      IF (CCRT) MODEL = 'CCSD'
      IF (CCR3) MODEL = 'CCSD'
      IF (CCR1A) MODEL = 'CCSD'
      IF (CCR1B) MODEL = 'CCSD'
C
      IOPT = 7
      CALL CC_WRRSP('R0',0,1,IOPT,MODEL,WORK(KT0AM),WORK(KT1AM),
     *              WORK(KT2AM),WORK(KEND0),LWRK0)
C
C---------------------------------
C     Comparison with MP2 vectors.
C---------------------------------
C
      MP2DIF = .FALSE.
      IF (MP2DIF) THEN
C
         CALL DCOPY(NT1AMX,WORK(KT1AM),1,WORK(KOMEG1),1)
         CALL DCOPY(NT2AMX,WORK(KT2AM),1,WORK(KOMEG1+NT1AMX),1)
C
         REWIND(LUIAJB)
         CALL READI(LUIAJB,IRAT*NT2AM(ISYMOP),WORK(KT2AM))
C
         CALL CCSD_GUESS(WORK(KT1AM),WORK(KT2AM),WORK(KFOCKD),
     *                   IALPR)
C
         CALL GPOPEN(LUTAM,'CCSD_TAM','UNKNOWN',' ','UNFORMATTED',
     &               IDUMMY,.FALSE.)
         REWIND(LUTAM)
         READ(LUTAM) (WORK(KT1AM+I-1), I = 1,NT1AMX)
         CALL GPCLOSE(LUTAM,'KEEP')
C        
         CALL DAXPY(NT1AMX,-1.0D0,WORK(KT1AM),1,WORK(KOMEG1),1)
C
C        CALL DAXPY(NT2AMX,-1.0D0,WORK(KT2AM),1,WORK(KOMEG1+NT1AMX),1)
C
C        CALL CC_PRAM(WORK(KOMEG1),PT1,1)
C
         LUTMP = 0
         CALL GPOPEN(LUTMP,' ','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND(LUTMP)
         WRITE(LUTMP) (WORK(KT2AM+I-1), I = 1,NT2AMX) 
         XTH = 1.0D0
         DO 600 I = 1,15
            REWIND(LUTMP)
            READ(LUTMP) (WORK(KT2AM+K-1), K = 1,NT2AMX) 
C
C        CALL DZERO(WORK(KT2AM),NT2AMX)
C
            ICOUNT = 0
            DO 610 J = 1,NT2AMX
C              IF ( ABS(WORK(KOMEG1+NT1AMX+J-1)) .GE. XTH) THEN
               IF ( ABS(WORK(KT2AM+J-1)) .LE. XTH) THEN
                  WORK(KT2AM+J-1) = ZERO                     
                  ICOUNT = ICOUNT + 1
               ELSE
                  WORK(KT2AM+J-1) = WORK(KOMEG1+NT1AMX+J-1) 
C
C                 WORK(KT2AM+J-1) = WORK(KT2AM+J-1) 
C    *                            + WORK(KOMEG1+NT1AMX+J-1)
C
C                 ICOUNT = ICOUNT + 1
               ENDIF
  610       CONTINUE
C
            WRITE(LUPRI,*)'ICOUNT=',ICOUNT
            IT1 = 1
            CALL CCSD_ECCSD(WORK(KT1AM),WORK(KT2AM),WORK(KFOCKD),
     *                      WORK(KEND1),LWRK1,EN2,POTNUC,ESCF,
     *                      ETY2,IPRINT,IT1,I)
 
            XTH = XTH/10.0D0
C
  600    CONTINUE 
         CALL GPCLOSE(LUTMP,'DELETE')
C
      ENDIF
C
C SLV98,OC
C
      IF (CCSLV) THEN
        XTNCCCU = DDOT(NT2AMX,WORK(KT2AM),1,WORK(KT2AM),1)
     *          + DDOT(NT1AMX,WORK(KT1AM),1,WORK(KT1AM),1)
        IF (ABS(XTNCCPR-XTNCCCU).LT.CVGTSOL) LSLTCVG = .TRUE.
        IF (IPRINT.GT.2) THEN
          WRITE(LUPRI,*)'Norm of T-amplitudes in this solvent it.:',XTNCCCU
          WRITE(LUPRI,*)'Norm of T-amplitudes in prev solvent it.:',XTNCCPR
          WRITE(LUPRI,*)'LSLTCVG: ',LSLTCVG
        ENDIF
        WRITE(LUPRI,*)
     *  ' Change in norm^2 of T-amplitudes in this solvent it.:',
     *  XTNCCCU-XTNCCPR
        XTNCCPR = XTNCCCU
      ENDIF
C
C-----------------------------------------------------
C     IF Triples corrections open files for integrals.
C-----------------------------------------------------
C
      IF (CCPT .OR. CCP3) THEN
C
C--------------------------------------------------
C        Calculate energy EN2 for CCSD(T) or CC(3).
C--------------------------------------------------
C
         CCSDT = .TRUE.
C
         IF (CCPT) THEN
            CC1BSV = CC1B
            CC1B   = .TRUE.
            CC1ASV = CC1A
            CC1A   = .TRUE.
         ENDIF
C
         CALL CCRHSN(WORK(KOMEG1),WORK(KOMEG2),WORK(KT1AM),
     *               WORK(KT2AM),WORK(KEND1),LWRK1)
C
         IF (CCPT) THEN
            CC1B   = CC1BSV
            CC1A   = CC1ASV
         ENDIF
C
         CALL CCSD_TCMEPK(WORK(KT2AM),0.5d0,1)
         ECCP1 = TWO*DDOT(NT1AMX,WORK(KT1AM),1,WORK(KOMEG1),1)
         ECCP2 = TWO*DDOT(NT2AMX,WORK(KT2AM),1,WORK(KOMEG2),1)
C
         CALL GPOPEN(LUTAM,'CCSD_TAM','UNKNOWN',' ','UNFORMATTED',
     &               IDUMMY,.FALSE.)
         REWIND (LUTAM)
         READ(LUTAM) (WORK(KT1AM+I-1), I = 1,NT1AMX)
         READ(LUTAM) (WORK(KT2AM+I-1), I = 1,NT2AMX)
         CALL GPCLOSE(LUTAM,'KEEP')
C
         IF ( CCR3 ) THEN
C
C-------------------------------------------
C           for perturbative correction CCT:
C           scale vector and add to t.
C-------------------------------------------
C
            WRITE(LUPRI,'(/,1X,A,/)')
     *       'Perturbational corrected amplitudes calculated'
            CALL CC_VSCAL(WORK(KOMEG1),WORK(KOMEG2),ZERO,
     *                    WORK(KEND1),LWRK1,1)
C
            CALL DAXPY(NT1AM(ISYMOP),XMONE,WORK(KOMEG1),1,
     *                 WORK(KT1AM),1)
            CALL DAXPY(NT2AM(ISYMOP),XMONE,WORK(KOMEG2),1,
     *                 WORK(KT2AM),1)
C
            IF ( IPRINT .GT. 10 ) THEN
               CALL AROUND('CCSD_ENERGY: third order (T1,T2)')
               RHO1N = DDOT(NT1AMX,WORK(KT1AM),1,WORK(KT1AM),1)
               RHO2N = DDOT(NT2AMX,WORK(KT2AM),1,WORK(KT2AM),1)
               WRITE(LUPRI,*) 'Norm of T1AM: ',RHO1N
               WRITE(LUPRI,*) 'Norm of T2AM: ',RHO2N
            ENDIF
C
            IF (IPRINT .GT. 45) THEN
               CALL CC_PRP(WORK(KOMEG1),WORK(KOMEG2),1,1,1)
            ENDIF
C
            CALL GPOPEN(LUTAM,'CCSD_TAM','UNKNOWN',' ','UNFORMATTED',
     *                  IDUMMY,.FALSE.)
            REWIND (LUTAM)
            WRITE(LUTAM) (WORK(KT1AM+I-1), I = 1,NT1AMX)
            WRITE(LUTAM) (WORK(KT2AM+I-1), I = 1,NT2AMX)
            CALL GPCLOSE(LUTAM,'KEEP')
C
            IF (IPRINT .GT. 4) THEN
C
               CALL AROUND('Largest amplitudes in pert. corr. ampl.')
C
               CALL DCOPY(NT1AMX,WORK(KT1AM),1,WORK(KOMEG1),1)
               CALL DCOPY(NT2AMX,WORK(KT2AM),1,WORK(KOMEG1+NT1AMX),1)
C
               CALL CC_PRAM(WORK(KOMEG1),PT1,1)
C
            ENDIF
C
         ENDIF
C
         ETOT = EN2 + ECCP1 + ECCP2
c        IF (CCPT) THEN
c           WRITE(LUPRI,'(20X,A,F20.10)') ' Total energy CCSD(T):',ETOT
c        ELSE
c           WRITE(LUPRI,'(20X,A,F20.10)') ' Total energy CC(3):',ETOT
c        ENDIF
c        WRITE(LUPRI,'(A,F12.10)')' T1 contribution:', ECCP1
c        WRITE(LUPRI,'(A,F12.10)')' T2 contribution:', ECCP2
C
         CCSDT  = CCSAV
         CCPTSV = CCPT
         CCP3SV = CCP3
         CCPT   = .FALSE.
         CCP3   = .FALSE.
C
         IF (CCSDT) THEN
            ITER = 1
            GOTO 240
         ENDIF
C
      ENDIF
C
C------------------------------------------------------------
C     Print and save (in ECCGRS) final ground state energies.
C------------------------------------------------------------
C
      WRITE(LUPRI,'(//)')
      CALL AROUND
     * ('Final results from the Coupled Cluster energy program')
      WRITE(LUPRI,'(/)')
      WRITE(LUPRI,'(12X,A,A,A,F32.10,/)') 'Total ',ETY0,' energy: ',ESCF
      WRITE(LUPRI,'(12X,A,A,A,F32.10,/)') 'Total ',ETY1,' energy: ',EINI
      WRITE(LURES,'(/)')
      WRITE(LURES,'(12X,A,A,A,F32.10)') 'Total ',ETY0,' energy: ',ESCF
      WRITE(LURES,'(12X,A,A,A,F32.10)') 'Total ',ETY1,' energy: ',EINI
      IF (.NOT. (MP2 .OR. NOCCIT))
     *WRITE(LUPRI,'(12X,A,A,A,F32.10)')   'Total ',ETY2,' energy: ',EN2
      IF (.NOT. (MP2 .OR. NOCCIT))
     *WRITE(LURES,'(12X,A,A,A,F32.10)')   'Total ',ETY2,' energy: ',EN2
      ECCGRS = EN2
      IF (CCPTSV .OR. CCP3SV) THEN
         WRITE(LUPRI,'(//,21X,A)') 'Perturbative triples corrections'
         WRITE(LUPRI,'(21X,A,/)')  '--------------------------------'
         WRITE(LUPRI,'(12X,A,F25.10)')
     *        'The E4 doubles and triples:', ECCP2
         WRITE(LUPRI,'(12X,A,F25.10)')
     *        'The E5 singles and triples:', ECCP1
         WRITE(LURES,'(//,21X,A)') 'Perturbative triples corrections'
         WRITE(LURES,'(21X,A,/)')  '--------------------------------'
         WRITE(LURES,'(12X,A,F25.10)')
     *        'The E4 doubles and triples:', ECCP2
         WRITE(LURES,'(12X,A,F25.10)')
     *        'The E5 singles and triples:', ECCP1
         IF (CCPTSV) THEN
            WRITE(LUPRI,'(/,12X,A,F31.10)') 'Total energy CCSD(T):',ETOT
            WRITE(LURES,'(/,12X,A,F31.10)') 'Total energy CCSD(T):',ETOT
         ELSE
            WRITE(LUPRI,'(/,12X,A,F31.10)') 'Total energy CC(3):  ',ETOT
            WRITE(LURES,'(/,12X,A,F31.10)') 'Total energy CC(3):  ',ETOT
         ENDIF
         ECCGRS = ETOT
      END IF
C
      CCP3 = CCP3SV
      CCPT = CCPTSV
C
C=======================================================
C     Calculate Intermediates for response calculations:
C
C      for cc2: E-intermediates
C
C      for ccsd also:
C
C         BF intermediate in ao.,
C         C & D intermediates,
C         Gamma intermediates.
C
C      OC 26-7-1995
C=======================================================
C
      IF (RSPIM2.AND.(.NOT.IMSKIP)) THEN
C
         RSPIM  = RSPIM2
C
         WRITE(LUPRI,'(/)')
         CALL AROUND( 'Calculating singlet intermediates for CCLR ')
         WRITE(LUPRI,'(/)')
C
         CCSAV = CCSDT
         CCSDT = .FALSE.
         CALL CCRHSN(WORK(KOMEG1),WORK(KOMEG2),WORK(KT1AM),WORK(KT2AM),
     *               WORK(KEND1),LWRK1)
         CCSDT = CCSAV
C
         IF (IPRINT .GT. 1) WRITE(LUPRI,'(/)')
         WRITE(LUPRI,'(12X,A)') 'E-intermediates calculated '
         WRITE(LUPRI,'(12X,A)') 'Fock-intermediate calculated '
C
         IF (.NOT. CC2 ) THEN
C
            WRITE(LUPRI,'(12X,A)') 'Gamma-intermediate calculated '
            WRITE(LUPRI,'(12X,A)') 'BF-intermediate calculated '
            WRITE(LUPRI,'(12X,A)') 'C-intermediate calculated '
            WRITE(LUPRI,'(12X,A)') 'D-intermediate calculated '
C
         ENDIF
         WRITE(LUPRI,'(/)')
C
      ELSE IF (RSPIM2.AND.IMSKIP) THEN
C
         RSPIM  = RSPIM2
         WRITE(LUPRI,'(12X,A)')
     &        'Intermediates assumed to be restart IM. '
C
      ENDIF
!
!----------------------------------------------------
!
!     Calculate the triplet global intermediates.
!
!----------------------------------------------------
!
      IF (TRIPIM .AND. (.NOT.IMSKIP) .AND. (.NOT.(CCS.OR.CC2))) THEN
!
         RSPIM = RSPIM2
!
         WRITE(LUPRI,'(/)')
         CALL AROUND( 'Calculating triplet intermediates for CCLR ')
         WRITE(LUPRI,'(/)')
!
         CALL CCRHSN3(WORK,LWORK)
!
         WRITE(LUPRI,'(12X,A)')
     &        'Triplet D and CD intermediate calculated '
!
      ENDIF
                                        
      RETURN
      END
C  /* Deck ccsd_guess */
      SUBROUTINE CCSD_GUESS(T1AM,T2AM,FCDIAG,IPRINT)
C
C     Written by Henrik Koch 27-Mar-1990.
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, TWO = 2.0D0, THREE = 3.0D0)
      DIMENSION T1AM(*),T2AM(*)
      DIMENSION FCDIAG(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
C-------------------------------------
C     Initial guess for t1 amplitudes.
C-------------------------------------
C
      CALL DZERO(T1AM,NT1AMX)
C
C-------------------------------------
C     Initial guess for t2 amplitudes.
C-------------------------------------
C
      DO 100 ISYMBJ = 1,NSYM
         ISYMAI = ISYMBJ
         DO 110 ISYMJ = 1,NSYM
            ISYMB = MULD2H(ISYMJ,ISYMBJ)
            DO 120 ISYMI = 1,NSYM
               ISYMA = MULD2H(ISYMI,ISYMAI)
               DO 130 J = 1,NRHF(ISYMJ)
                  KOFFJ = IRHF(ISYMJ) + J
                  DO 140 B = 1,NVIR(ISYMB)
                     KOFFB = IVIR(ISYMB) + B
                     NBJ = IT1AM(ISYMB,ISYMJ)+NVIR(ISYMB)*(J-1)+B
                     DO 150 I = 1,NRHF(ISYMI)
                        KOFFI = IRHF(ISYMI) + I
                        DO 160 A = 1,NVIR(ISYMA)
                           KOFFA = IVIR(ISYMA) + A
                           NAI = IT1AM(ISYMA,ISYMI)+NVIR(ISYMA)*(I-1)+A
C
                           IF (NAI .GT. NBJ) GOTO 160
C
                           NAIBJ = IT2AM(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
C
                           T2AM(NAIBJ) = T2AM(NAIBJ)/
     *                                 (FCDIAG(KOFFI) + FCDIAG(KOFFJ)
     *                                - FCDIAG(KOFFA) - FCDIAG(KOFFB))
C
  160                   CONTINUE
  150                CONTINUE
  140             CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      IF (IPRINT .GT. 15) THEN
         CALL AROUND('T1 Guess vector')
         DO 200 ISYMI = 1,NSYM
            ISYMA = ISYMI
            KOFF  = IT1AM(ISYMA,ISYMI) + 1
            NVIRA = NVIR(ISYMA)
            NRHFI = NRHF(ISYMI)
            CALL OUTPUT(T1AM(KOFF),1,NVIRA,1,NRHFI,NVIRA,NRHFI,1,LUPRI)
  200    CONTINUE
C
         CALL AROUND('T2 Guess vector')
         DO 250 ISYMBJ = 1,NSYM
            ISYMAI = ISYMBJ
            KOFF   = IT2AM(ISYMAI,ISYMBJ) + 1
            NTOTAI = NT1AM(ISYMAI)
            CALL OUTPAK(T2AM(KOFF),NTOTAI,1,LUPRI)
  250    CONTINUE
      ENDIF
C
      RETURN
      END
C  /* Deck ccsd_nxtam */
      SUBROUTINE CCSD_NXTAM(T1AM,T2AM,OMEGA1,OMEGA2,FCDIAG)
C
C     Written by Henrik Koch 27-Mar-1990.
C     Brueckner bit by Rika Kobayashi 1992.
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, TWO = 2.0D0, THREE = 3.0D0)
      DIMENSION T1AM(*),T2AM(*)
      DIMENSION OMEGA1(*),OMEGA2(*)
      DIMENSION FCDIAG(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
c     IF (.NOT. (CCD.OR.LBRUK)) THEN
C
         DO 100 ISYMI = 1,NSYM
            ISYMA = ISYMI
            DO 110 I = 1,NRHF(ISYMI)
               KOFFI = IRHF(ISYMI) + I
               DO 120 A = 1,NVIR(ISYMA)
C
                  KOFFA = IVIR(ISYMA) + A
                  NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A
C
                  OMEGA1(NAI) = T1AM(NAI) + OMEGA1(NAI)/
     *                       (FCDIAG(KOFFI) - FCDIAG(KOFFA))
C
  120          CONTINUE
  110       CONTINUE
  100    CONTINUE
C
c     ENDIF
c     IF (LBRUK) CALL DCOPY(NVIRT*NRHFT,OMEGA1,1,T1AM,1)
C
      DO 200 ISYMBJ = 1,NSYM
         ISYMAI = ISYMBJ
         DO 210 ISYMJ = 1,NSYM
            ISYMB = MULD2H(ISYMJ,ISYMBJ)
            DO 220 ISYMI = 1,NSYM
               ISYMA = MULD2H(ISYMI,ISYMAI)
               DO 230 J = 1,NRHF(ISYMJ)
                  KOFFJ = IRHF(ISYMJ) + J
                  DO 240 B = 1,NVIR(ISYMB)
                     KOFFB = IVIR(ISYMB) + B
                     NBJ = IT1AM(ISYMB,ISYMJ)+NVIR(ISYMB)*(J-1)+B
                     DO 250 I = 1,NRHF(ISYMI)
                        KOFFI = IRHF(ISYMI) + I
                        DO 260 A = 1,NVIR(ISYMA)
                           KOFFA = IVIR(ISYMA) + A
                           NAI = IT1AM(ISYMA,ISYMI)+NVIR(ISYMA)*(I-1)+A
C
                           IF (NAI .GT. NBJ) GOTO 260
C
                           NAIBJ = IT2AM(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
C
                           OMEGA2(NAIBJ) = T2AM(NAIBJ) + OMEGA2(NAIBJ)/
     *                                   (FCDIAG(KOFFI) + FCDIAG(KOFFJ)
     *                                  - FCDIAG(KOFFA) - FCDIAG(KOFFB))
C
  260                   CONTINUE
  250                CONTINUE
  240             CONTINUE
  230          CONTINUE
  220       CONTINUE
  210    CONTINUE
  200 CONTINUE
C
      RETURN
      END
C  /* Deck ccsd_eccsd */
      SUBROUTINE CCSD_ECCSD(T1AM,T2AM,FCDIAG,WORK,LWORK,XECCSD,POTNUC,
     *                      ESCF,ETY,IALPR,IT1,ITER)
C
C     Written by Henrik Koch 27-Mar-1990.
C
C     Ove Christiansen 23-1-1996: Introduction of Finite field contribution.
C                                 IT1 = 0 : no amplitudes on disk
C                                 IT1 = 1 : t1 amplitudes read from disk
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
#include <iratdef.h>
      DIMENSION FCDIAG(*)
      DIMENSION T1AM(*),T2AM(*),WORK(*)
      CHARACTER*5 ETY
      LOGICAL LEXIST
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <ccfield.h>
#include <ccinftap.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      XECCSD = ESCF
C
C---------------------------------
C     Dynamic allocation of space.
C---------------------------------
C
      KIAJB  = 1
      KEND1  = KIAJB + NT2AMX
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient spaces in ECCSD')
      ENDIF
C
      REWIND(LUIAJB)
      CALL READI(LUIAJB,IRAT*NT2AMX,WORK)
C
      DO 100 ISYMJ = 1,NSYM
         DO 110 ISYMB = 1,NSYM
            ISYMBJ = MULD2H(ISYMB,ISYMJ)
            ISYMAI = ISYMBJ
            DO 120 ISYMI = 1,NSYM
               ISYMBI = MULD2H(ISYMB,ISYMI)
               ISYMA  = MULD2H(ISYMI,ISYMAI)
               ISYMAJ = ISYMBI
C
               DO 130 J = 1,NRHF(ISYMJ)
                  DO 140 B = 1,NVIR(ISYMB)
C
                     KBJ = IT1AM(ISYMB,ISYMJ)
                     NBJ = KBJ + NVIR(ISYMB)*(J-1) + B
C
                     DO 150 I = 1,NRHF(ISYMI)
C
                        KBI = IT1AM(ISYMB,ISYMI)
                        NBI = KBI + NVIR(ISYMB)*(I-1) + B
C
                        DO 160 A = 1,NVIR(ISYMA)
C
                           KAI = IT1AM(ISYMA,ISYMI)
                           NAI = KAI + NVIR(ISYMA)*(I-1) + A
                           KAJ = IT1AM(ISYMA,ISYMJ)
                           NAJ = KAJ + NVIR(ISYMA)*(J-1) + A
C
                           NAIBJ = IT2AM(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
                           NAJBI = IT2AM(ISYMAJ,ISYMBI) + INDEX(NAJ,NBI)
C
                           IF (ISYMB .EQ. ISYMJ) THEN
                              XECCSD = XECCSD
     *                            + (TWO*WORK(NAIBJ) - WORK(NAJBI))*
     *                              (T2AM(NAIBJ) + T1AM(NAI)*T1AM(NBJ))
                           ELSE
                              XECCSD = XECCSD
     *                    + (TWO*WORK(NAIBJ) - WORK(NAJBI))*T2AM(NAIBJ)
                           ENDIF
C
  160                   CONTINUE
  150                CONTINUE
  140             CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
C-------------------------------------------------------------------
C     Add field dependent energy in case of finite field ONEelectron
C     Perturbation. The AO integral from ONEP is already scaled with
C     the fieldstrengths!!!
C-------------------------------------------------------------------
C
      DO 13 IF = 1, NFIELD
        IF (NONHF) THEN
C
         KONEP  = 1
         KT1AM  = KONEP  + N2BST(ISYMOP)
         KLAMDP = KT1AM  + NT1AMX
         KLAMDH = KLAMDP + NLAMDT
         KEND1  = KLAMDH + NLAMDT
         LWRK1  = LWORK  - KEND1
         IF ( LWRK1 .LT. 0 )
     *     CALL QUIT(' Too little workspace in ccsd_eccsd-2')
C
         CALL DZERO(WORK(KONEP),N2BST(ISYMOP))
         FF = EFIELD(IF)
         CALL CC_ONEP(WORK(KONEP),WORK(KEND1),LWRK1,FF,1,LFIELD(IF))
C
         IF (.NOT.(CCS.OR.CCP2)) THEN
C
            IF ( IT1 .EQ. 1 ) THEN
               LUTAM = 0
               CALL GPOPEN(LUTAM,'CCSD_TAM','UNKNOWN',' ',
     *                     'UNFORMATTED',IDUMMY,.FALSE.)
               REWIND(LUTAM)
               READ(LUTAM) (WORK(KT1AM - 1 + J), J = 1,NT1AMX)
               CALL GPCLOSE(LUTAM,'KEEP')
            ELSE IF (IT1 .EQ. 0) THEN
               CALL DZERO(WORK(KT1AM),NT1AMX)
            ELSE
               CALL QUIT('IT1 should be 0 or 1 in ccsd_eccsd')
            ENDIF
         ENDIF
C
         CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),
     *               WORK(KEND1),LWRK1)
C
         ISYFAO = 1
         ISYMPA = 1
         ISYMHO = 1
C
         CALL CC_FCKMO(WORK(KONEP),WORK(KLAMDP),WORK(KLAMDH),
     *                 WORK(KEND1),LWRK1,ISYFAO,ISYMPA,ISYMHO)
C
         IF (IALPR .GT. 50) THEN
            CALL AROUND( 'In CCSD_ECCSD: Pert. matrix' )
            CALL CC_PRFCKMO(WORK(KONEP),ISYMOP)
         ENDIF
C
         DO 200 ISYMI = 1, NSYM
C
            DO 210 I = 1, NRHF(ISYMI)
C
               KII  = KONEP - 1 + IFCRHF(ISYMI,ISYMI)
     *              + NORB(ISYMI)*(I - 1) + I
C
               XECCSD = XECCSD + 2.0D00*WORK(KII)
C
  210       CONTINUE
  200    CONTINUE
C
        ENDIF
  13  CONTINUE
C
      WRITE(LUPRI,'(1X,A,I3,A,A,A,F23.16)')
     *  'Iter.',ITER,': Coupled cluster ',ETY,' energy :  ',XECCSD
C
      CALL FLSHFO(LUPRI)
C
      RETURN
      END
C  /* Deck ccsd_iajb */
      SUBROUTINE CCSD_IAJB(XAIBJ,T1AM,WORK,LWORK)
C
C     Written by Henrik Koch 27-Mar-1990.
C
C     Small modifications by Asger Halkier 22/5 - 1998 for extra
C     MO integrals needed for gradients and frozen core FOP.
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
#include <maxorb.h>
#include <maxash.h>
#include <mxorb.h>
#include <mxcent.h>
#include <aovec.h>
#include <iratdef.h>
      DIMENSION XAIBJ(*),T1AM(*),WORK(*),INDEXA(MXCORB)
#include <ccorb.h>
#include <ccisao.h>
#include <blocks.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <cbieri.h>
#include <distcl.h>
#include <eribuf.h>
#include <ccfro.h>
#include <ccfop.h>
#include <ccsections.h>
#include <ccfield.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
C-----------------------------------------
C     Initialize the XAIBJ integral array.
C-----------------------------------------
C
      CALL DZERO(XAIBJ,NT2AM(ISYMOP))
C
C---------------------------------
C     Dynamic allocation of space.
C---------------------------------
C
      KLAMDP = 1
      KLAMDH = KLAMDP + NLAMDT
      KEND1  = KLAMDH + NLAMDT
C
      KCMO   = KEND1
      KDNSHF = KCMO   + NLAMDS
      KFCKHF = KDNSHF + N2BAST
      KEND1  = KFCKHF + N2BAST
C
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSD_IAJB')
      ENDIF
C
C-----------------------------------------------------
C     Calculate the lamda matrices and get CMO vector:
C-----------------------------------------------------
C
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),T1AM,WORK(KEND1),LWRK1)
C
C---------------------------------------------------------------------
C     initialize CMO vector, SCF density and SCF AO-Fock matrix:
C       we include in the SCF AO-Fock matrix ONLY fields added 
C       already at the SCF level (i.e. the ``relaxed'' fields)
C       this matrix is needed for relaxed CC2 response, the
C       numerical Xksi and Eta vectors (CC_FDXI, CC_FDETA)
C---------------------------------------------------------------------
C
      LUSIFC = 0
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',
     *            IDUMMY,.FALSE.)
      REWIND(LUSIFC)
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ(LUSIFC)
      READ(LUSIFC)
      READ(LUSIFC) (WORK(KCMO+I-1),I=1,NLAMDS)
      CALL GPCLOSE(LUSIFC,'KEEP')
C
      CALL CMO_REORDER(WORK(KCMO),WORK(KEND1),LWRK1) 
C

      CALL CC_AODENS(WORK(KCMO),WORK(KCMO),WORK(KDNSHF),1,1,
     *               WORK(KEND1),LWRK1)
C
      CALL CCRHS_ONEAO(WORK(KFCKHF),WORK(KEND1),LWRK1)
      DO IF = 1, NFIELD
        IF ( .NOT. NHFFIELD(IF) ) THEN
          CALL CC_ONEP(WORK(KFCKHF),WORK(KEND1),LWRK1,EFIELD(IF),
     *                 1,LFIELD(IF))
        END IF
      END DO
C
C--------------------------------------
C     Additional work space allocation.
C--------------------------------------
C
      IF ((FROIMP) .AND. ((RELORB) .OR. ((CCFOP) .AND. (MP2)))) THEN
C
         KCMO  = KEND1
         KFRIN = KCMO  + NLAMDS
         KEND1 = KFRIN + NT2FRO(1)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient space in CCSD_IAJB')
         ENDIF
C
         CALL DZERO(WORK(KCMO),NLAMDS)
         CALL DZERO(WORK(KFRIN),NT2FRO(1))
C
C----------------------------------------------
C     Calculate the FULL MO coefficient matrix.
C----------------------------------------------
C
         CALL CMO_ALL(WORK(KCMO),WORK(KEND1),LWRK1)
C
      ENDIF
C
C====================================================
C     Start the loop over distributions of integrals.
C====================================================
C
      IF (DEBUG) THEN
C        IPRERI = 5
         WRITE(LUPRI,'(1X,A,I10)') 'LWORK = ',LWORK
      END IF
C
      IF (DIRECT) THEN
         DTIME  = SECOND()
         IF (HERDIR) THEN
            CALL HERDI1(WORK(KEND1),LWRK1,IPRERI)
         ELSE
            KCCFB1 = KEND1
            KINDXB = KCCFB1 + MXPRIM*MXCONT
            KEND1  = KINDXB + (8*MXSHEL*MXCONT + 1)/IRAT
            LWRK1  = LWORK  - KEND1
            CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     &                  KODPP1,KODPP2,KRDPP1,KRDPP2,
     &                  KFREE,LFREE,KEND1,WORK(KCCFB1),WORK(KINDXB),
     &                  WORK(KEND1),LWRK1,IPRERI)
            KEND1 = KFREE
            LWRK1 = LFREE
         ENDIF
         NTOSYM = 1
      ELSE
         NTOSYM = NSYM
      ENDIF
C
      THRDIS = 1.0D-8
      ICOUNT1 = 0
      ICOUNT2 = 0
C
      KENDSV = KEND1
      LWRKSV = LWRK1
C
      DO 100 ISYMD1 = 1,NTOSYM
C
         IF (DIRECT) THEN
            IF (HERDIR) THEN
               NTOT = MAXSHL
            ELSE
               NTOT = MXCALL
            ENDIF
         ELSE
            NTOT = NBAS(ISYMD1)
         ENDIF
C
         DO 110 ILLL = 1,NTOT
C
C---------------------------------------------
C           If direct calculate the integrals.
C---------------------------------------------
C
            IF (DIRECT) THEN
C
               KEND1 = KENDSV
               LWRK1 = LWRKSV
C
               IF (HERDIR) THEN
                  CALL HERDI2(WORK(KEND1),LWRK1,INDEXA,ILLL,NUMDIS,
     &                        IPRERI)
               ELSE
                  CALL ERIDI2(ILLL,INDEXA,NUMDIS,0,0,
     &                        WORK(KODCL1),WORK(KODCL2),WORK(KODBC1),
     &                        WORK(KODBC2),WORK(KRDBC1),WORK(KRDBC2),
     &                        WORK(KODPP1),WORK(KODPP2),WORK(KRDPP1),
     &                        WORK(KRDPP2),WORK(KCCFB1),WORK(KINDXB),
     &                        WORK(KEND1), LWRK1,IPRERI)
               ENDIF
C
               KRECNR = KEND1
               KEND1  = KRECNR + (NBUFX(0) - 1)/IRAT + 1
               LWRK1  = LWORK  - KEND1
               IF (LWRK1 .LT. 0) THEN
                  CALL QUIT('Insufficient core in CCRHSN')
               END IF
C
            ELSE
               NUMDIS = 1
            ENDIF
C
C-----------------------------------------------------
C           Loop over number of distributions in disk.
C-----------------------------------------------------
C
            DO 120 IDEL2 = 1,NUMDIS
C
               IF (DIRECT) THEN
                  IDEL  = INDEXA(IDEL2)
                  ISYMD = ISAO(IDEL)
               ELSE
                  IDEL  = IBAS(ISYMD1) + ILLL
                  ISYMD = ISYMD1
               ENDIF
C
               ISYMB  = ISYMD
               ISYDIS = MULD2H(ISYMD,ISYMOP)
C
C-----------------------------------------------
C              Dynamic allocation of work space.
C-----------------------------------------------
C
               KXINT = KEND1
               KSCR1 = KXINT + NDISAO(ISYDIS)
               KSCR2 = KSCR1 + NBAST*NBAST
               KEND2 = KSCR2 + NBAST*NRHFT
               LWRK2 = LWORK - KEND2
C
               IF (LWRK2 .LT. 0) THEN
                  CALL QUIT('Insufficient space in CCSD_IAJB')
               ENDIF
C
C-----------------------------------------
C              Read in batch of integrals.
C-----------------------------------------
C
               CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND2),LWRK2,
     *                     WORK(KRECNR),DIRECT)
C
C-----------------------------------------
C              compute the AO-Fock matrix:
C-----------------------------------------
C
               CALL CC_AOFOCK(WORK(KXINT),WORK(KDNSHF),WORK(KFCKHF),
     *                        WORK(KEND2),LWRK2,IDEL,ISYMD,1)
C
C-----------------------------------------------
C              Calculate integrals (cJ|dk)
C              needed for frozen core gradients.
C-----------------------------------------------
C
               IF ((FROIMP) .AND.
     *             ((RELORB) .OR. ((CCFOP) .AND. (MP2)))) THEN
C
                  CALL CC_FRCOIN(WORK(KFRIN),WORK(KXINT),WORK(KCMO),
     *                           WORK(KEND2),LWRK2,IDEL,ISYMD)
C
               ENDIF
C
C---------------------------------------------------
C              Transform one index in the integrals.
C---------------------------------------------------
C
               DO 130 ISYMG = 1,NSYM
C
                  ISYMAB = MULD2H(ISYMG,ISYDIS)
                  ISYMJ  = ISYMG
                  ISYMBJ = MULD2H(ISYMB,ISYMJ)
                  ISYMAI = MULD2H(ISYMBJ,ISYMOP)
C
                  IF (ISYMAI .GT. ISYMBJ) GOTO 130
C
                  KOFF1 = KXINT  + IDSAOG(ISYMG,ISYDIS)
                  KOFF2 = KLAMDP + ILMRHF(ISYMJ)
C
                  IF (LWRK2 .LT. NNBST(ISYMAB)*NRHF(ISYMJ)) THEN
                     CALL QUIT('Insufficient core in CCSD_IAJB')
                  ENDIF
C
C--------------------------------------------------------
C                 Analyse size of integral distributions.
C--------------------------------------------------------
C
                  DO 140 G = 1,NBAS(ISYMG)
C
                     KOFFG = KXINT + IDSAOG(ISYMG,ISYDIS)
     *                             + NNBST(ISYMAB)*(G - 1)
                     NAB   = NNBST(ISYMAB)
C
                     DO 150 IAB = 1,NAB
                        IF (ABS(WORK(KOFFG+IAB)) .GT. THRDIS) GOTO 158
  150                CONTINUE
C
C                    WRITE(LUPRI,*) 'ISYMD,IDEL,ISYMG,G : ',
C    *                   ISYMD,IDEL,ISYMG,G
                     ICOUNT1 = ICOUNT1 + 1
C
  158                CONTINUE
C
                     ICOUNT2 = ICOUNT2 + 1
C
  140             CONTINUE
C
C-------------------------------------------------------------------
C                 Transform the gamma index in the integral (AB|GD).
C-------------------------------------------------------------------
C
                  NNBSAB = MAX(NNBST(ISYMAB),1)
                  NBASG  = MAX(NBAS(ISYMG),1)
C
                  CALL DGEMM('N','N',NNBST(ISYMAB),NRHF(ISYMJ),
     *                       NBAS(ISYMG),ONE,WORK(KOFF1),NNBSAB,
     *                       WORK(KOFF2),NBASG,ZERO,WORK(KEND2),
     *                       NNBSAB)
C
C------------------------------------------------------------------
C                 Transform integrals and add to the result vector.
C------------------------------------------------------------------
C
                  KOFF4  = IT2AM(ISYMAI,ISYMBJ) + 1
C
                  CALL CCSD_AIBJ2(WORK(KEND2),XAIBJ(KOFF4),WORK(KLAMDP),
     *                            WORK(KLAMDH),WORK(KSCR1),WORK(KSCR2),
     *                            IDEL,ISYMD,ISYMJ,ISYMAB)
C
                  IF (CCSDT.OR.(CCPT.OR.CCP3).OR.
     *                (CCRT.OR.CCR3.OR.CCR1A.OR.CCR1B )) THEN
C
C------------------------------------------------------
C                    Construct I(kd,c) for fixed alpha.
C------------------------------------------------------
C
                     KINT3 = KEND2
                     KINT4 = KINT3 + NT1AM(ISYMAB)*NVIR(ISYMG)
                     KSCR3 = KINT4 + NT1AM(ISYMAB)*NVIR(ISYMG)
                     KEND3 = KSCR3 + NT1AM(ISYMAB)*NBAS(ISYMG)
                     LWRK3 = LWORK - KEND3
C
                     IF (LWRK3 .LT. 0) THEN
                        CALL QUIT('Insufficient space in CCSD_IAJB')
                     END IF

                     KOFF5 = KXINT + IDSAOG(ISYMG,ISYDIS)
C
                     CALL CCSD_AIBJ3(WORK(KOFF5),WORK(KINT3),
     *                               WORK(KINT4),WORK(KLAMDP),
     *                               WORK(KLAMDH),WORK(KSCR1),
     *                               WORK(KSCR2),WORK(KSCR3),
     *                               IDEL,ISYMD,ISYMG,ISYMAB)
C
                  END IF
C
  130          CONTINUE
C
  120       CONTINUE
C
  110    CONTINUE
C
  100 CONTINUE
C
      KENDSV = KEND1
      LWRKSV = LWRK1
C
C-------------------------------------
C     write AO-Fock matrix to file:
C-------------------------------------
C
      LUFCK = 0
      CALL GPOPEN(LUFCK,'CC_FCKREF','UNKNOWN',' ','UNFORMATTED',
     *            IDUMMY,.FALSE.)
      REWIND(LUFCK)
      WRITE(LUFCK)(WORK(KFCKHF + I-1),I = 1,N2BST(ISYMOP))
      CALL GPCLOSE(LUFCK,'KEEP' )
C
      IF (IPRINT .GT.150) THEN
         CALL AROUND( 'Fock AO matrix for reference state:' )
         CALL CC_PRFCKAO(WORK(KFCKHF),1)
      ENDIF                  
C
      IF (ANAAOD) THEN
        CALL AROUND('Analysis of integral distributions')
C
        WRITE(LUPRI,'(10X,/,A,D10.5)') 'Thereshold in analysis : ',
     &       THRDIS
        WRITE(LUPRI,'(10X,A,I7)')'Total number of dist.           : ',
     *                        ICOUNT2
        WRITE(LUPRI,'(10X,A,I7)')'Total number larger than thr.   : ',
     *                        ICOUNT2 - ICOUNT1
        WRITE(LUPRI,'(10X,A,I7)')'Total number smaller than thr.  : ',
     *                        ICOUNT1
C
        IF (IPRINT .GT. 45) THEN
           CALL AROUND('(ia|jb) integral vector')
           DO 250 ISYMBJ = 1,NSYM
              ISYMAI = ISYMBJ
              KOFF   = IT2AM(ISYMAI,ISYMBJ) + 1
              NTOTAI = NT1AM(ISYMAI)
              CALL OUTPAK(XAIBJ(KOFF),NTOTAI,1,LUPRI)
  250      CONTINUE
        ENDIF
      END IF
C
C-------------------------------------
C     Write integrals (cJ|dk) to disk.
C-------------------------------------
C
      IF ((FROIMP) .AND. ((RELORB) .OR. ((CCFOP) .AND. (MP2)))) THEN
         LUCJDK = 0
         CALL GPOPEN(LUCJDK,'INCJDK','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND(LUCJDK)
         WRITE(LUCJDK) (WORK(KFRIN+I-1), I = 1,NT2FRO(1))
         CALL GPCLOSE(LUCJDK,'KEEP')
      ENDIF
C
      IF (CCSDT.OR.(CCPT.OR.CCP3).OR.
     *    (CCRT.OR.CCR3.OR.CCR1A.OR.CCR1B)) THEN
C
C------------------------------------
C        Sort integrals (kc,d alpha).
C------------------------------------
C
         ISYINT = ISYMOP
         CALL CC3_SORT1(WORK,LWORK,1,ISYINT)
C
      END IF
C
      RETURN
      END
C  /* Deck ccsd_aibj2 */
      SUBROUTINE CCSD_AIBJ2(XINT,XAIBJ,XLAMDP,XLAMDH,
     *                      SCR1,SCR2,IDEL,ISYMD,ISYMJ,ISYMAB)
C
C     Written by Henrik Koch 27-Mar-1990.
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION XINT(*),XAIBJ(*), SCR1(*),SCR2(*)
      DIMENSION XLAMDP(*),XLAMDH(*)
#include <priunit.h>
#include <ccinftap.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      DO 100 J = 1,NRHF(ISYMJ)
C
         KOFF1 = NNBST(ISYMAB)*(J-1) + 1
C
         CALL CCSD_SYMSQ(XINT(KOFF1),ISYMAB,SCR1)
C
C--------------------------------------------------
C        Transformation of the A-index to occupied.
C--------------------------------------------------
C
         KOFF3 = 1
         DO 110 ISYMI = 1,NSYM
C
            ISYMA = ISYMI
            ISYMB = MULD2H(ISYMA,ISYMAB)
C
            KOFF1 = IAODIS(ISYMA,ISYMB) + 1
            KOFF2 = ILMRHF(ISYMI) + 1
C
            NBASA = MAX(NBAS(ISYMA),1)
            NBASB = MAX(NBAS(ISYMB),1)
C
            CALL DGEMM('T','N',NBAS(ISYMB),NRHF(ISYMI),NBAS(ISYMA),
     *                 ONE,SCR1(KOFF1),NBASA,XLAMDP(KOFF2),
     *                 NBASA,ZERO,SCR2(KOFF3),NBASB)
C
            KOFF3 = KOFF3 + NBAS(ISYMB)*NRHF(ISYMI)
C
  110    CONTINUE
C
C-------------------------------------------------
C        Transformation of the B-index to virtual.
C-------------------------------------------------
C
         KOFF2 = 1
         DO 120 ISYMI = 1,NSYM
C
            ISYMB = MULD2H(ISYMI,ISYMAB)
            ISYMA = ISYMB
C
            KOFF1 = ILMVIR(ISYMA) + 1
            KOFF3 = IT1AM(ISYMA,ISYMI) + 1
C
            NBASB = MAX(NBAS(ISYMB),1)
            NVIRA = MAX(NVIR(ISYMA),1)
C
            CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMB),
     *                 ONE,XLAMDH(KOFF1),NBASB,SCR2(KOFF2),
     *                 NBASB,ZERO,SCR1(KOFF3),NVIRA)
C
            KOFF2 = KOFF2 + NBAS(ISYMB)*NRHF(ISYMI)
C
  120    CONTINUE
C
C------------------------------------------
C        Write out integrals used in CCSDT.
C------------------------------------------
C
         IF (CCSDT.OR.(CCPT.OR.CCP3).OR.
     *      (CCRT.OR.CCR3.OR.CCR1A.OR.CCR1B)) THEN
C
            ISYMI  = ISYMJ
            ISYMCK = ISYMAB
            ISYCKI = MULD2H(ISYMCK,ISYMI)
C
            I  = J
            ID = IDEL - IBAS(ISYMD)
C
            IOFF = ICKID(ISYCKI,ISYMD) + NCKI(ISYCKI)*(ID - 1)
     *           + ICKI(ISYMCK,ISYMI) + NT1AM(ISYMCK)*(I - 1) + 1
C
            IF (NT1AM(ISYMCK) .GT. 0) THEN
               CALL PUTWA2(LUTOC,FNTOC,SCR1,IOFF,NT1AM(ISYMCK))
            ENDIF
         ENDIF
C
C--------------------------------------------------
C        Add the contribution to the result vector.
C--------------------------------------------------
C
         ISYMB  = ISYMD
         ISYMBJ = MULD2H(ISYMB,ISYMJ)
         ISYMAI = ISYMAB
C
         DO 130 B = 1,NVIR(ISYMB)
C
            NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J-1) + B
C
            IF (ISYMAI .EQ. ISYMBJ) THEN
               NTOTAI = NBJ
               KOFF1 = NBJ*(NBJ - 1)/2 + 1
            ELSE
               NTOTAI = NT1AM(ISYMAI)
               KOFF1  = NT1AM(ISYMAI)*(NBJ - 1) + 1
            ENDIF
C
            KOFF2 = ILMVIR(ISYMB) + NBAS(ISYMD)*(B-1) + IDEL
     *              - IBAS(ISYMD)
C
            CALL DAXPY(NTOTAI,XLAMDH(KOFF2),SCR1,1,XAIBJ(KOFF1),1)
C
  130    CONTINUE
C
  100 CONTINUE
C
      RETURN
      END
C  /* Deck ccsd_aibj3 */
      SUBROUTINE CCSD_AIBJ3(XINT,XINT3,XINT4,XLAMDP,XLAMDH,SCR1,SCR2,
     *                      SCR3,IDEL,ISYDEL,ISYMG,ISYMAB)
C
C     Written by Henrik Koch 27-Mar-1990.
C     Modified asm
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION XINT(*),XINT3(*),XINT4(*),SCR1(*),SCR2(*),SCR3(*)
      DIMENSION XLAMDP(*),XLAMDH(*)
#include <priunit.h>
#include <ccinftap.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      ISYMKD = ISYMAB
C
      DO 100 G = 1,NBAS(ISYMG)
C
         KOFF1 = NNBST(ISYMAB)*(G-1) + 1
C
         CALL CCSD_SYMSQ(XINT(KOFF1),ISYMAB,SCR1)
C
C--------------------------------------------------
C        Transformation of the A-index to occupied.
C--------------------------------------------------
C
         KOFF3 = 1
         DO 110 ISYMK = 1,NSYM
C
            ISYMA = ISYMK
            ISYMB = MULD2H(ISYMA,ISYMAB)
C
            KOFF1 = IAODIS(ISYMA,ISYMB) + 1
            KOFF2 = ILMRHF(ISYMK) + 1
C
            NBASA = MAX(NBAS(ISYMA),1)
            NBASB = MAX(NBAS(ISYMB),1)
C
            CALL DGEMM('T','N',NBAS(ISYMB),NRHF(ISYMK),NBAS(ISYMA),
     *                 ONE,SCR1(KOFF1),NBASA,XLAMDP(KOFF2),
     *                 NBASA,ZERO,SCR2(KOFF3),NBASB)
C
            KOFF3 = KOFF3 + NBAS(ISYMB)*NRHF(ISYMK)
C
  110    CONTINUE
C
C-------------------------------------------------
C        Transformation of the B-index to virtual.
C-------------------------------------------------
C
         KOFF2 = 1
         DO 120 ISYMK = 1,NSYM
C
            ISYMB  = MULD2H(ISYMK,ISYMAB)
            ISYMC  = ISYMB
            ISYMCK = MULD2H(ISYMC,ISYMK)
C
            KOFF1 = ILMVIR(ISYMC) + 1
            KOFF3 = NT1AM(ISYMCK)*(G - 1) + IT1AM(ISYMC,ISYMK) + 1
C
            NBASB = MAX(NBAS(ISYMB),1)
            NVIRC = MAX(NVIR(ISYMC),1)
C
            CALL DGEMM('T','N',NVIR(ISYMB),NRHF(ISYMK),NBAS(ISYMB),
     *                 ONE,XLAMDH(KOFF1),NBASB,SCR2(KOFF2),
     *                 NBASB,ZERO,SCR3(KOFF3),NVIRC)
C
            KOFF2 = KOFF2 + NBAS(ISYMB)*NRHF(ISYMK)
C
  120    CONTINUE
C
  100 CONTINUE
C
C--------------------------------
C     Transform gamma index to d.
C--------------------------------
C
      ISYMCK = ISYMAB
      ISYMD  = ISYMG
C
      NBASG  = MAX(NBAS(ISYMG),1)
      NTOTCK = MAX(NT1AM(ISYMCK),1)
C
      KOFF = ILMVIR(ISYMG) + 1
C
      CALL DGEMM('N','N',NT1AM(ISYMCK),NVIR(ISYMD),NBAS(ISYMG),ONE,
     *           SCR3,NTOTCK,XLAMDH(KOFF),NBASG,ZERO,XINT3,NTOTCK)
C
C-------------------------------
C     Dump to disk (kc|d alpha).
C-------------------------------
C
      IA     = IDEL - IBAS(ISYDEL)
      ISYMA  = ISYDEL
      ISYCKD = MULD2H(ISYMCK,ISYMD)
C
      LENGTH = NT1AM(ISYMCK)*NVIR(ISYMD)
C
      IOFF = ICKDAO(ISYCKD,ISYMA) + NCKATR(ISYCKD)*(IA - 1)
     *     + ICKATR(ISYMCK,ISYMD) + 1
C
      IF (LENGTH .GT. 0) THEN
         CALL PUTWA2(LU3SRT,FN3SRT,XINT3,IOFF,LENGTH)
      ENDIF
C
      RETURN
      END
C  /* Deck inidat */
      BLOCK DATA INIDAT
C
C     Initialize MULD2H, FRORHF and FROVIR in common block /CCORB/
C
#include <ccorb.h>
C
      PARAMETER (NLIST = 8*MAXFRO)
      DATA FRORHF /NLIST * .FALSE./
      DATA FROVIR /NLIST * .FALSE./
C
      DATA MULD2H/1,2,3,4,5,6,7,8,
     *            2,1,4,3,6,5,8,7,
     *            3,4,1,2,7,8,5,6,
     *            4,3,2,1,8,7,6,5,
     *            5,6,7,8,1,2,3,4,
     *            6,5,8,7,2,1,4,3,
     *            7,8,5,6,3,4,1,2,
     *            8,7,6,5,4,3,2,1/
C
      END
C  /* Deck ccsd_init1 */
      SUBROUTINE CCSD_INIT1(WORK,LWORK)
C
C     Henrik Koch and Alfredo Sanchez.       29-Jun-1994
C
C     Set up indexing arrays
C
C     FREEZE OC230899
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
      DIMENSION WORK(LWORK)
C
      EXTERNAL INIDAT
C
#include <maxorb.h>
#include <mxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <inftap.h>
#include <symsq.h>
#include <ccisao.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
C-------------------------------------
C     Read in information from sirius.
C-------------------------------------
C
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
C
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC) NSYM, NORBTS, NBAST, NLAMDS, (NRHFS(I),I=1,NSYM),
     &              (NORBS(I),I=1,NSYM), (NBAS(I),I=1,NSYM), DUM, DUM
C
C
      IF (FREEZE) THEN
         WRITE(LUPRI,*) ' I am freezing!'
C
         KFOCKD = 1
         KFCS   = KFOCKD + NORBTS
         KFVS   = KFCS   + NSYM
         KEND1  = KFVS   + NSYM
         LEND1  = LWORK  - KEND1
C
c        CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
c        READ (LUSIFC)
         READ (LUSIFC) (WORK(KFOCKD+I-1), I=1,NORBTS)
C
         CALL CC_FREEZER(WORK(KFOCKD),NORBTS,WORK(KFCS),WORK(KFVS),
     *                   WORK(KEND1),LEND1)
C
      ENDIF
C
      CALL GPCLOSE(LUSIFC,'KEEP')
C
C-----------------------------
C     Construct rest of CCORB.
C-----------------------------
C
      NNBASX = (NBAST*(NBAST+1))/2
      N2BASX = NBAST*NBAST
C
      NORBT  = 0
      NRHFT  = 0
      NRHFTS = 0
      N2BAST = 0
C
      ICOUN1 = 0
      ICOUN2 = 0
      IOFF   = 0
C
      DO 5 ISYM = 1,NSYM
C
         NVIRS(ISYM) = NORBS(ISYM) - NRHFS(ISYM)
C
         NRHF(ISYM)  = NRHFS(ISYM) - NRHFFR(ISYM)
         NVIR(ISYM)  = NVIRS(ISYM) - NVIRFR(ISYM)
         NORB(ISYM)  = NRHF(ISYM)  + NVIR(ISYM)
C
         NORBT  = NORBT  + NORB(ISYM)
         NRHFT  = NRHFT  + NRHF(ISYM)
         NRHFTS = NRHFTS + NRHFS(ISYM)
         N2BAST = N2BAST + NBAS(ISYM)*NBAS(ISYM)
C
        IORB(ISYM) = ICOUN1
        IBAS(ISYM) = ICOUN2
C
        ICOUN1 = ICOUN1 + NORB(ISYM)
        ICOUN2 = ICOUN2 + NBAS(ISYM)
C
        DO 6 I = 1,NBAS(ISYM)
C
             IOFF = IOFF + 1
             ISAO(IOFF) = ISYM
C
    6   CONTINUE
C
    5 CONTINUE
C
      NVIRT  = NORBT  - NRHFT
      NVIRTS = NORBTS - NRHFTS
C
c     IF (IPRINT .GT. 20) THEN
         CALL AROUND('Information from CCORB')
         WRITE(LUPRI,1) 'NBAS   :',(NBAS(I),   I=1,NSYM)
         WRITE(LUPRI,1) 'IBAS   :',(IBAS(I),   I=1,NSYM)
         WRITE(LUPRI,1) 'NRHF   :',(NRHF(I),   I=1,NSYM)
         WRITE(LUPRI,1) 'NVIR   :',(NVIR(I),   I=1,NSYM)
         WRITE(LUPRI,1) 'NRHFS  :',(NRHFS(I),  I=1,NSYM)
         WRITE(LUPRI,1) 'NVIRS  :',(NVIRS(I),  I=1,NSYM)
         WRITE(LUPRI,1) 'NRHFFR :',(NRHFFR(I), I=1,NSYM)
         WRITE(LUPRI,1) 'NVIRFR :',(NVIRFR(I), I=1,NSYM)
	 WRITE(LUPRI,1) 'NORBTS :',NORBTS
	 WRITE(LUPRI,1) 'NORBT  :',NORBT
	 WRITE(LUPRI,1) 'N2BAST :',N2BAST
	 WRITE(LUPRI,1) 'N2BASX :',N2BASX
	 WRITE(LUPRI,1) 'NNBASX :',NNBASX
c     END IF
C
C-------------------------------
C     Construct frozen matrices.
C-------------------------------
C
      IF (FROIMP) THEN
C
         DO 50 ISYM = 1,NSYM
C
            IF (NRHFFR(ISYM) .GT. MAXFRO) THEN
               WRITE(LUPRI,'(//,1X,2A,I3)') 'ERROR: Maximum number of ',
     &              'frozen orbitals per symmetry is:',MAXFRO
               CALL QUIT('Too many frozen orbitals')
            END IF
C
            DO 51 I = 1,NRHFFR(ISYM)
                  FRORHF(I,ISYM) = .TRUE.
                  KFRRHF(I,ISYM) = I
   51       CONTINUE
C
            IF (NVIRFR(ISYM) .GT. MAXFRO) THEN
               WRITE(LUPRI,'(//,1X,2A,I3)') 'ERROR: Maximum number of ',
     &              'frozen orbitals per symmetry is:',MAXFRO
               CALL QUIT('Too many frozen orbitals')
            END IF
C
            DO 52 I = 1,NVIRFR(ISYM)
               JORB = NVIRS(ISYM) - I + 1
               FROVIR(JORB,ISYM) = .TRUE.
               KFRVIR(I,ISYM)    = JORB
   52       CONTINUE
C
   50    CONTINUE
C
      ELSE IF (FROEXP) THEN
C
         DO 60 ISYM = 1,NSYM
C
            DO 61 I = 1,NRHFFR(ISYM)
               II = KFRRHF(I,ISYM)
               FRORHF(II,ISYM) = .TRUE.
   61       CONTINUE
C
            DO 62 I = 1,NVIRFR(ISYM)
               II = KFRVIR(I,ISYM)
               FROVIR(II,ISYM) = .TRUE.
   62       CONTINUE
C
   60    CONTINUE
C
      END IF
C
C------------------------------------------
C     Calculate the number of t-amplitudes.
C------------------------------------------
C
      DO 100 ISYMAI = 1,NSYM
         NT1AM(ISYMAI) = 0
         NT1AO(ISYMAI) = 0
	 DO 200 ISYMI = 1,NSYM
            ISYMA = MULD2H(ISYMAI,ISYMI)
            NT1AM(ISYMAI) = NT1AM(ISYMAI) + NVIR(ISYMA) * NRHF(ISYMI)
            NT1AO(ISYMAI) = NT1AO(ISYMAI) + NBAS(ISYMA) * NRHF(ISYMI)
  200    CONTINUE
  100 CONTINUE
C
      DO 300 ISAIBJ = 1,NSYM
         NT2AM(ISAIBJ)  = 0
         NT2AO(ISAIBJ)  = 0
         NT2AMA(ISAIBJ) = 0
         NT2AMT(ISAIBJ) = 0
         DO 400 ISYMBJ = 1,NSYM
            ISYMAI = MULD2H(ISYMBJ,ISAIBJ)
            IF (ISYMBJ .GT. ISYMAI) THEN
               NT2AM(ISAIBJ) = NT2AM(ISAIBJ) +
     &                         NT1AM(ISYMAI) * NT1AM(ISYMBJ)
               NT2AO(ISAIBJ) = NT2AO(ISAIBJ) +
     &                         NT1AO(ISYMAI) * NT1AO(ISYMBJ)
               NT2AMA(ISAIBJ)= NT2AM(ISAIBJ)
               NT2AMT(ISAIBJ)= NT2AM(ISAIBJ) + NT2AMA(ISAIBJ)     
            ELSE IF (ISYMBJ .EQ. ISYMAI) THEN
               NT2AM(ISAIBJ) = NT2AM(ISAIBJ) +
     &                         NT1AM(ISYMAI) * (NT1AM(ISYMBJ)+1)/2
               NT2AO(ISAIBJ) = NT2AO(ISAIBJ) +
     &                         NT1AO(ISYMAI) * (NT1AO(ISYMBJ)+1)/2
               NT2AMA(ISAIBJ)= NT2AM(ISAIBJ)
               NT2AMT(ISAIBJ)= NT2AM(ISAIBJ) + NT2AMA(ISAIBJ)  
            END IF
  400    CONTINUE
  300 CONTINUE
C
      NT1AMX = NT1AM(1)
      NT1AOX = NT1AO(1)
      NT2AMX = NT2AM(1)
      NT2AOX = NT2AO(1)
C
      ICOUN1 = 0
      DO 450 ISYM = 1,NSYM
C
         NNBST(ISYM) = 0
         N2BST(ISYM) = 0
C
         DO 460 ISYMB = 1,NSYM
C
            ISYMA = MULD2H(ISYMB,ISYM)
C
            N2BST(ISYM) = N2BST(ISYM) + NBAS(ISYMA)*NBAS(ISYMB)
C
            IF (ISYMB .GT. ISYMA) THEN
               NNBST(ISYM) = NNBST(ISYM) + NBAS(ISYMA)*NBAS(ISYMB)
            ELSE IF (ISYMB .EQ. ISYMA) THEN
               NNBST(ISYM) = NNBST(ISYM) + NBAS(ISYMA)*(NBAS(ISYMA)+1)/2
            ENDIF
C
  460    CONTINUE
C
         I2BST(ISYM) = ICOUN1
C
         ICOUN1 = ICOUN1 + N2BST(ISYM)
C
  450 CONTINUE
C
      DO 500 ISYMD = 1,NSYM
         NDISAO(ISYMD)   = 0
         NDSRHF(ISYMD)   = 0
         NDISAOSQ(ISYMD) = 0
         NDSRHFSQ(ISYMD) = 0
         NT2BCD(ISYMD)   = 0
         NT2BGD(ISYMD)   = 0
         DO 510 ISYMG = 1,NSYM
            ISYMAB = MULD2H(ISYMG,ISYMD)
            NDISAO(ISYMD) = NDISAO(ISYMD) + NNBST(ISYMAB)*NBAS(ISYMG)
            NDSRHF(ISYMD) = NDSRHF(ISYMD) + NNBST(ISYMAB)*NRHF(ISYMG)
            NDISAOSQ(ISYMD)=NDISAOSQ(ISYMD)+N2BST(ISYMAB)*NBAS(ISYMG)
            NDSRHFSQ(ISYMD)=NDSRHFSQ(ISYMD)+N2BST(ISYMAB)*NRHF(ISYMG)
            NT2BCD(ISYMD) = NT2BCD(ISYMD) + NT1AM(ISYMAB)*NRHF(ISYMG)
            NT2BGD(ISYMD) = NT2BGD(ISYMD) + NT1AO(ISYMAB)*NRHF(ISYMG)
  510    CONTINUE
  500 CONTINUE
C
      ICOUN1 = 0
      ICOUN2 = 0
      ICOUN3 = 0
      ICOUN4 = NRHFT
      ICOUN5 = 0
      ICOUN6 = 0
      DO 600 ISYMP = 1,NSYM
         ICOUN1 = ICOUN1 + NBAS(ISYMP)*NORB(ISYMP)
         ICOUN2 = ICOUN2 + NBAS(ISYMP)*NRHF(ISYMP)
         ICOUN5 = ICOUN5 + NBAS(ISYMP)*NRHFS(ISYMP)
C
         IRHF(ISYMP) = ICOUN3
         IVIR(ISYMP) = ICOUN4
         ICOUN3 = ICOUN3 + NRHF(ISYMP)
         ICOUN4 = ICOUN4 + NVIR(ISYMP)
C
  600 CONTINUE
      NLAMDT = ICOUN1
      NLMRHF = ICOUN2
      NLRHSI = ICOUN5
C
      DO 610 ISYMK = 1,NSYM
         ICOUN1 = 0
         ICOUN2 = 0
         ICOUN3 = 0
         ICOUN4 = 0
         ICOUN5 = 0
         ICOUN6 = 0
         ICOUN7 = 0
         ICOUN8 = 0
         ICOUN9 = 0
         ICOU10 = 0
         ICOU11 = 0
         ICOU12 = 0
         ICOU13 = 0
         ICOU14 = 0
         ICOU15 = 0
         ICOU16 = 0
         ICOU17 = 0
         ICOU18 = 0
         ICOU19 = 0
         DO 620 ISYMJ = 1,NSYM
C
            ISYMI  = MULD2H(ISYMJ,ISYMK)
C
            IT1AM(ISYMI,ISYMJ)  = ICOUN1
            IT1AO(ISYMI,ISYMJ)  = ICOUN5
            IT1AMT(ISYMI,ISYMJ) = ICOU11
            IT1AOT(ISYMI,ISYMJ) = ICOU12
            IEMAT1(ISYMI,ISYMJ) = ICOU15
            IMATAV(ISYMI,ISYMJ) = ICOU18
C
            ICOUN1 = ICOUN1 + NRHF(ISYMJ)*NVIR(ISYMI)
            ICOUN5 = ICOUN5 + NRHF(ISYMJ)*NBAS(ISYMI)
            ICOU11 = ICOU11 + NRHF(ISYMI)*NVIR(ISYMJ)
            ICOU12 = ICOU12 + NRHF(ISYMI)*NBAS(ISYMJ)
            ICOU15 = ICOU15 + NVIR(ISYMI)*NBAS(ISYMJ)
            ICOU18 = ICOU18 + NBAS(ISYMI)*NVIR(ISYMJ)
C
            IF (ISYMJ .GT. ISYMI) THEN
               IT2AM(ISYMI,ISYMJ) = ICOUN2
               IT2AM(ISYMJ,ISYMI) = ICOUN2
               ICOUN2 = ICOUN2 + NT1AM(ISYMI)*NT1AM(ISYMJ)
               IT2AO(ISYMI,ISYMJ) = ICOUN6
               IT2AO(ISYMJ,ISYMI) = ICOUN6
               ICOUN6 = ICOUN6 + NT1AO(ISYMI)*NT1AO(ISYMJ)
            ELSE IF (ISYMK .EQ. 1) THEN
               IT2AM(ISYMI,ISYMJ) = ICOUN2
               ICOUN2 = ICOUN2 + NT1AM(ISYMI)*(NT1AM(ISYMI)+1)/2
               IT2AO(ISYMI,ISYMJ) = ICOUN6
               ICOUN6 = ICOUN6 + NT1AO(ISYMI)*(NT1AO(ISYMI)+1)/2
            ENDIF
C
            IT2BGD(ISYMI,ISYMJ)   = ICOUN8
            IT2BCD(ISYMI,ISYMJ)   = ICOUN9
            IDSRHF(ISYMI,ISYMJ)   = ICOU10
            IT2BGT(ISYMI,ISYMJ)   = ICOU13
            IT2BCT(ISYMI,ISYMJ)   = ICOU14
            ICKALP(ISYMI,ISYMJ)   = ICOU16
            ICKATR(ISYMI,ISYMJ)   = ICOU17
            IDSRHFSQ(ISYMI,ISYMJ) = ICOU19
C
            ICOUN3 = ICOUN3 + NVIR(ISYMI)*NBAS(ISYMJ)
            ICOUN4 = ICOUN4 + NRHF(ISYMI)*NRHF(ISYMJ)
C
            IT2SQ(ISYMI,ISYMJ) = ICOUN7
C
            ICOUN7 = ICOUN7 + NT1AM(ISYMI)*NT1AM(ISYMJ)
            ICOUN8 = ICOUN8 + NT1AO(ISYMI)*NRHF(ISYMJ)
            ICOUN9 = ICOUN9 + NT1AM(ISYMI)*NRHF(ISYMJ)
            ICOU10 = ICOU10 + NNBST(ISYMI)*NRHF(ISYMJ)
            ICOU13 = ICOU13 + NT1AO(ISYMJ)*NRHF(ISYMI)
            ICOU14 = ICOU14 + NT1AM(ISYMJ)*NRHF(ISYMI)
            ICOU16 = ICOU16 + NT1AM(ISYMI)*NBAS(ISYMJ)
            ICOU17 = ICOU17 + NT1AM(ISYMI)*NVIR(ISYMJ)
            ICOU19 = ICOU19 + N2BST(ISYMI)*NRHF(ISYMJ)
C
  620    CONTINUE
C
         NEMAT1(ISYMK) = ICOUN3
         NMATIJ(ISYMK) = ICOUN4
         NMATAV(ISYMK) = ICOU18
C
  610 CONTINUE
C
      DO 630 ISYMK = 1,NSYM
         ICOUN1 = 0
         ICOUN2 = 0
         DO 640 ISYMJ = 1,NSYM
            ISYMI = MULD2H(ISYMJ,ISYMK)
C
            IF (ISYMJ .GT. ISYMI) THEN
               ICOUN1 = ICOUN1 + NMATIJ(ISYMI)*NMATIJ(ISYMJ)
            ELSE IF (ISYMK .EQ. 1) THEN
               ICOUN1 = ICOUN1 + NMATIJ(ISYMI)*(NMATIJ(ISYMI)+1)/2
            ENDIF
C
            ICOUN2 = ICOUN2 + NVIR(ISYMI)*NVIR(ISYMJ)
C
  640    CONTINUE
C
         NGAMMA(ISYMK) = ICOUN1
         NMATAB(ISYMK) = ICOUN2
C
  630 CONTINUE
C
C
C--------------------------------------------------------
C     Section for calculating index arrays needed in left
C     hand side transformation. Asger Halkier 30/10-1995!
C     Revised 7/3-1996 for index arrays for densities!
C--------------------------------------------------------
C
      DO 550 ISYIJK = 1,NSYM
         ICOUN1 = 0
         ICOUN2 = 0
         ICOUN3 = 0
         ICOUN4 = 0
         ICOUN5 = 0
         ICOUN6 = 0
         DO 560 ISYMK = 1,NSYM
            ISYMIJ = MULD2H(ISYMK,ISYIJK)
            IMAIJK(ISYMIJ,ISYMK) = ICOUN1
            IT2AIJ(ISYMIJ,ISYMK) = ICOUN2
            IMAIJA(ISYMIJ,ISYMK) = ICOUN3
            ID2IJG(ISYMIJ,ISYMK) = ICOUN4
            ID2AIG(ISYMIJ,ISYMK) = ICOUN5
            ID2ABG(ISYMIJ,ISYMK) = ICOUN6
            ICOUN1 = ICOUN1 + NMATIJ(ISYMIJ)*NRHF(ISYMK)
            ICOUN2 = ICOUN2 + NVIR(ISYMIJ)*NMATIJ(ISYMK)
            ICOUN3 = ICOUN3 + NMATIJ(ISYMIJ)*NVIR(ISYMK)
            ICOUN4 = ICOUN4 + NMATIJ(ISYMIJ)*NBAS(ISYMK)
            ICOUN5 = ICOUN5 + NT1AM(ISYMIJ)*NBAS(ISYMK)
            ICOUN6 = ICOUN6 + NMATAB(ISYMIJ)*NBAS(ISYMK)
  560    CONTINUE
         NMAIJK(ISYIJK) = ICOUN1
         NT2AIJ(ISYIJK) = ICOUN2
         NMAIJA(ISYIJK) = ICOUN3
         ND2IJG(ISYIJK) = ICOUN4
         ND2AIG(ISYIJK) = ICOUN5
         ND2ABG(ISYIJK) = ICOUN6
  550 CONTINUE
C
      DO 570 ISIJKD = 1,NSYM
         ICOUN1 = 0
         ICOUN2 = 0
         ICOUN3 = 0
         DO 580 ISYMD = 1,NSYM
            ISYIJK = MULD2H(ISYMD,ISIJKD)
            I3ODEL(ISYIJK,ISYMD) = ICOUN1
            I3ORHF(ISYIJK,ISYMD) = ICOUN2
            I3OVIR(ISYIJK,ISYMD) = ICOUN3
            ICOUN1 = ICOUN1 + NMAIJK(ISYIJK)*NBAS(ISYMD)
            ICOUN2 = ICOUN2 + NMAIJK(ISYIJK)*NRHF(ISYMD)
            ICOUN3 = ICOUN3 + NMAIJK(ISYIJK)*NVIR(ISYMD)
  580    CONTINUE
         N3ODEL(ISIJKD) = ICOUN1
         N3ORHF(ISIJKD) = ICOUN2
         N3OVIR(ISIJKD) = ICOUN3
  570 CONTINUE
C
      ICOUN = 0
C
      DO 590 ISYM = 1,NSYM
C
         IFCKDO(ISYM) = ICOUN
         ICOUN = ICOUN + NORB(ISYM)*NRHF(ISYM)
         IFCKDV(ISYM) = ICOUN
         ICOUN = ICOUN + NORB(ISYM)*NVIR(ISYM)
C
  590 CONTINUE
C
      ICOUN1 = 0
      ICOUN2 = NLMRHF
      ICOUN7 = 0
      ICOUN8 = NLRHSI
      DO 700 ISYMI = 1,NSYM
C
         ILMRHF(ISYMI) = ICOUN1
         ILMVIR(ISYMI) = ICOUN2
         ICOUN1 = ICOUN1 + NBAS(ISYMI)*NRHF(ISYMI)
         ICOUN2 = ICOUN2 + NBAS(ISYMI)*NVIR(ISYMI)
C
         ILRHSI(ISYMI) = ICOUN7
         ILVISI(ISYMI) = ICOUN8
         ICOUN7 = ICOUN7 + NBAS(ISYMI)*NRHFS(ISYMI)
         ICOUN8 = ICOUN8 + NBAS(ISYMI)*NVIRS(ISYMI)
C
         ICOUN3 = 0
         ICOUN4 = 0
         ICOUN5 = 0
         ICOUN6 = 0
         ICOUN9 = 0
C
         DO 710 ISYMJ = 1,NSYM
C
            ISYMK = MULD2H(ISYMJ,ISYMI)
C
            IDSAOG(ISYMJ,ISYMI)   = ICOUN3
            IMATIJ(ISYMK,ISYMJ)   = ICOUN4
            IGAMMA(ISYMK,ISYMJ)   = ICOUN5
            IGAMMA(ISYMJ,ISYMK)   = ICOUN5
            IMATAB(ISYMK,ISYMJ)   = ICOUN6
            IDSAOGSQ(ISYMJ,ISYMI) = ICOUN9
C
            ICOUN3 = ICOUN3 + NNBST(ISYMK)*NBAS(ISYMJ)
            ICOUN4 = ICOUN4 + NRHF(ISYMK)*NRHF(ISYMJ)
            ICOUN6 = ICOUN6 + NVIR(ISYMK)*NVIR(ISYMJ)
            ICOUN9 = ICOUN9 + N2BST(ISYMK)*NBAS(ISYMJ)
C
            IF (ISYMJ .GT. ISYMK) THEN
               ICOUN5 = ICOUN5 + NMATIJ(ISYMK)*NMATIJ(ISYMJ)
            ELSE IF (ISYMI .EQ. 1) THEN
               ICOUN5 = ICOUN5 + NMATIJ(ISYMJ)*(NMATIJ(ISYMJ)+1)/2
            ENDIF
C
  710    CONTINUE
  700 CONTINUE
C
      DO 720 ISYMAB = 1,NSYM
         ICOUN1 = 0
         ICOUN2 = 0
         DO 730 ISYMB = 1,NSYM
C
            ISYMA = MULD2H(ISYMB,ISYMAB)
C
            IAODIS(ISYMA,ISYMB) = ICOUN1
            IAODPK(ISYMA,ISYMB) = ICOUN2
            IAODPK(ISYMB,ISYMA) = ICOUN2
C
            ICOUN1 = ICOUN1 + NBAS(ISYMA)*NBAS(ISYMB)
            IF (ISYMB .GT. ISYMA) THEN
               ICOUN2 = ICOUN2 + NBAS(ISYMA)*NBAS(ISYMB)
            ELSE IF (ISYMAB .EQ. 1) THEN
               ICOUN2 = ICOUN2 + NBAS(ISYMB)*(NBAS(ISYMB)+1)/2
            ENDIF
C
  730    CONTINUE
  720 CONTINUE
C
      DO 800 ISYM = 1,NSYM
C
         ICOUNT = 0
         DO 810 ISYMK = 1,NSYM
C
            ISYMP = MULD2H(ISYMK,ISYM)
C
            IFCRHF(ISYMP,ISYMK) = ICOUNT
C
            ICOUNT = ICOUNT + NORB(ISYMP)*NRHF(ISYMK)
C
  810    CONTINUE
C
         DO 820 ISYMC = 1,NSYM
C
            ISYMP = MULD2H(ISYMC,ISYM)
C
            IFCVIR(ISYMP,ISYMC) = ICOUNT
C
            ICOUNT = ICOUNT + NORB(ISYMP)*NVIR(ISYMC)
C
  820    CONTINUE
C
  800 CONTINUE
C
C
      DO 900 ISYM = 1,NSYM
C
         ICOUNT = 0
         ICOUN1 = 0
         DO 910 ISYMJ = 1,NSYM
C
            ISYMI  = MULD2H(ISYMJ,ISYM)
            IT2AOS(ISYMI,ISYMJ) = ICOUNT
C
            ICOUNT = ICOUNT + NT1AO(ISYMI)*NT1AO(ISYMJ)
            ICOUN1 = ICOUN1 + NT1AM(ISYMI)*NT1AM(ISYMJ)
C
  910    CONTINUE
C
         NT2AOS(ISYM) = ICOUNT
         NT2SQ(ISYM)  = ICOUN1
C
  900 CONTINUE
C
      DO 920 ISYMAB = 1,NSYM
C
         DO 930 ISYMB = 1,NSYM
C
            ISYMA = MULD2H(ISYMB,ISYMAB)
C
            DO 940 B = 1,NBAS(ISYMB)
C
               DO 950 A = 1,NBAS(ISYMA)
C
                  NABSQ = I2BST(ISYMAB) + IAODIS(ISYMA,ISYMB)
     *                  + NBAS(ISYMA)*(B - 1) + A
C
                  IF (ISYMA .EQ. ISYMB) THEN
                     NAB = INDEX(A,B)
                  ELSE IF (ISYMB .GT.ISYMA) THEN
                     NAB = NBAS(ISYMA)*(B - 1) + A
                  ELSE
                     NAB = NBAS(ISYMB)*(A - 1) + B
                  END IF
C
                  NABPK = IAODPK(ISYMA,ISYMB) + NAB
C
                  IADRPK(NABSQ) = NABPK
C
  950          CONTINUE
  940       CONTINUE
  930    CONTINUE
  920 CONTINUE
C
      DO 1000 ISYM = 1,NSYM
C
         ICOUN1 = 0
         DO 1010 ISYMJ = 1,NSYM
C
            ISYMI = MULD2H(ISYMJ,ISYM)
C
            IF (ISYMI .GT. ISYMJ) GOTO 1010
C
            IMIJP(ISYMI,ISYMJ) = ICOUN1
            IMIJP(ISYMJ,ISYMI) = ICOUN1
C
            IF (ISYMI .EQ. ISYMJ) THEN
               ICOUN1 = ICOUN1 + NRHF(ISYMI)*(NRHF(ISYMI) + 1)/2
            ELSE
               ICOUN1 = ICOUN1 + NRHF(ISYMI)*NRHF(ISYMJ)
            ENDIF
C
 1010    CONTINUE
C
         NMIJP(ISYM) = ICOUN1
C
 1000 CONTINUE
C
C
      DO 1020 ISYM = 1,NSYM
C
         ICOUNT = 0
         ICOUN1 = 0
         DO 1030 ISYMJ = 1,NSYM
C
            ISYMI = MULD2H(ISYMJ,ISYM)
C
            IT2ORT(ISYMI,ISYMJ)  = ICOUNT
            IT2AOIJ(ISYMI,ISYMJ) = ICOUN1
C
            ICOUNT = ICOUNT + NNBST(ISYMI)*NMIJP(ISYMJ)
            ICOUN1 = ICOUN1 + NT1AO(ISYMI)*NMATIJ(ISYMJ)
C
 1030    CONTINUE
C
         NT2ORT(ISYM)  = ICOUNT
         NT2AOIJ(ISYM) = ICOUN1
C
 1020 CONTINUE
C
      DO 1040 ISYCKA = 1,NSYM
C
         ICOUN1 = 0
         ICOUN2 = 0
         ICOUN3 = 0
         ICOUN4 = 0
         DO 1050 ISYMA = 1,NSYM
C
            ISYMCK = MULD2H(ISYMA,ISYCKA)
C
            ICKA(ISYMCK,ISYMA)   = ICOUN1
            ICKI(ISYMCK,ISYMA)   = ICOUN2
            ISAIK(ISYMCK,ISYMA)  = ICOUN2
            ICKATR(ISYMCK,ISYMA) = ICOUN3
            ICKASR(ISYMCK,ISYMA) = ICOUN4
C
            ICOUN1 = ICOUN1 + NT1AM(ISYMCK)*NBAS(ISYMA)
            ICOUN2 = ICOUN2 + NT1AM(ISYMCK)*NRHF(ISYMA)
            ICOUN3 = ICOUN3 + NT1AM(ISYMCK)*NVIR(ISYMA)
            ICOUN4 = ICOUN4 + NMATAB(ISYMCK)*NRHF(ISYMA)
C
 1050    CONTINUE
C
         NCKA(ISYCKA)   = ICOUN1
         NCKI(ISYCKA)   = ICOUN2
         NCKATR(ISYCKA) = ICOUN3
         NCKASR(ISYCKA) = ICOUN4
C
 1040 CONTINUE
C
      DO 1060 ISYMJ = 1,NSYM
C
         ICOUN2 = 0
         ICOUN3 = 0
C
         DO 1065 ISYMD = 1,NSYM
C
            ISYCKA = MULD2H(ISYMD,ISYMJ)
C
            ICOUN2 = ICOUN2 + NCKI(ISYCKA)*NBAS(ISYMD)
            ICOUN3 = ICOUN3 + NCKI(ISYCKA)*NRHF(ISYMD)
C
 1065    CONTINUE
C
         NTOTOC(ISYMJ) = ICOUN2
         NTRAOC(ISYMJ) = ICOUN3
C
 1060 CONTINUE
C
      DO 1070 JSAIKJ = 1,NSYM
C
         ICOUN1 = 0
         ICOUN2 = 0
         ICOUN3 = 0
         ICOUN4 = 0
         ICOUN5 = 0
         ICOUN6 = 0
         DO 1080 ISYMJ = 1, NSYM
C
            ISYAIK = MULD2H(JSAIKJ,ISYMJ)
C
            ISAIKJ(ISYAIK,ISYMJ) = ICOUN1
            ICKITR(ISYAIK,ISYMJ) = ICOUN1
            ICKID(ISYAIK,ISYMJ)  = ICOUN2
            ICKAD(ISYAIK,ISYMJ)  = ICOUN3
            ICKDAO(ISYAIK,ISYMJ) = ICOUN4
            ICKBD(ISYAIK,ISYMJ)  = ICOUN5
            IT2SP(ISYAIK,ISYMJ)  = ICOUN6
C
            ICOUN1 = ICOUN1 + NCKI(ISYAIK)*NRHF(ISYMJ)
            ICOUN2 = ICOUN2 + NCKI(ISYAIK)*NBAS(ISYMJ)
            ICOUN3 = ICOUN3 + NCKA(ISYAIK)*NVIR(ISYMJ)
            ICOUN4 = ICOUN4 + NCKATR(ISYAIK)*NBAS(ISYMJ)
            ICOUN5 = ICOUN5 + NCKATR(ISYAIK)*NVIR(ISYMJ)
            ICOUN6 = ICOUN6 + NCKI(ISYAIK)*NVIR(ISYMJ)
C
 1080    CONTINUE
C
         NCKIJ(JSAIKJ) = ICOUN1
C
 1070 CONTINUE
C
      DO 1090 ISYJIK = 1,NSYM
C
         ICOUN1 = 0
         DO 1100 ISYMK = 1,NSYM
C
            ISYMJI = MULD2H(ISYJIK,ISYMK)
C
            ICOUN1 = ICOUN1 + NMATIJ(ISYMJI)*NRHF(ISYMK)
 1100    CONTINUE
C
         NMAJIK(ISYJIK) = ICOUN1
C
 1090 CONTINUE
C
      DO 1110 JSJIKA = 1,NSYM
C
         ICOUN1 = 0
         ICOUN2 = 0
         ICOUN3 = 0
         DO 1120 ISYMA = 1,NSYM
C
            ISYJIK = MULD2H(JSJIKA,ISYMA)
C
            ISJIKA(ISYJIK,ISYMA) = ICOUN1
            ISJIK(ISYJIK,ISYMA)  = ICOUN2
            ISAIKL(ISYJIK,ISYMA) = ICOUN3
C
            ICOUN1 = ICOUN1 + NMAJIK(ISYJIK)*NVIR(ISYMA)
            ICOUN2 = ICOUN2 + NMATIJ(ISYJIK)*NRHF(ISYMA)
            ICOUN3 = ICOUN3 + NT1AM(ISYJIK)*NMATIJ(ISYMA)
C
 1120    CONTINUE
 1110 CONTINUE
C
C------------------------------------------------------------------
C     Section for making index matrices for general Lamda matrices.
C     Needed for linear transformation. OC 10-2-1995
C------------------------------------------------------------------
C
      DO 1200 ISYM = 1,NSYM
C
         ICOUN1 = 0
         ICOUN2 = 0
         ICOUN3 = 0
C
         DO 1210 ISYM2 = 1,NSYM
C
            ISYM1  = MULD2H(ISYM,ISYM2)
            ICOUN1 = ICOUN1 + NBAS(ISYM1)*NORB(ISYM2)
            ICOUN2 = ICOUN2 + NBAS(ISYM1)*NRHF(ISYM2)
            ICOUN3 = ICOUN3 + NORB(ISYM1)*NRHF(ISYM2)
C
 1210    CONTINUE
C
         NGLMDT(ISYM) = ICOUN1
         NGLMRH(ISYM) = ICOUN2
         NLRHFR(ISYM) = ICOUN3
         ICOUN1 = 0
C
         DO 1220 ISYM2 = 1,NSYM
C
            ISYM1  = MULD2H(ISYM,ISYM2)
            IGLMRH(ISYM1,ISYM2) = ICOUN1
            IGLMVI(ISYM1,ISYM2) = ICOUN2
C
            ICOUN1 = ICOUN1 + NBAS(ISYM1)*NRHF(ISYM2)
            ICOUN2 = ICOUN2 + NBAS(ISYM1)*NVIR(ISYM2)
C
 1220    CONTINUE
C
 1200 CONTINUE
C
      DO 1230 ISYMD = 1,NSYM
         DO 1240 ISYMTR  = 1,NSYM
            NT2MMO(ISYMD,ISYMTR) = 0
            NT2MAO(ISYMD,ISYMTR) = 0
            ISYCIJ = MULD2H(ISYMD,ISYMTR)
            DO 1250 ISYMJ = 1,NSYM
               ISYMCI = MULD2H(ISYMJ,ISYCIJ)
               NT2MMO(ISYMD,ISYMTR) = NT2MMO(ISYMD,ISYMTR) +
     *                                NT1AM(ISYMCI)*NRHF(ISYMJ)
               NT2MAO(ISYMD,ISYMTR) = NT2MAO(ISYMD,ISYMTR) +
     *                                NT1AO(ISYMCI)*NRHF(ISYMJ)
 1250       CONTINUE
 1240    CONTINUE
 1230 CONTINUE
C
C----------------------------------------------------
C     Section for extra frozen core gradient indices.
C     Asger Halkier 22/5 - 1998.
C----------------------------------------------------
C
      CALL CC_INIFRO(WORK,LWORK)
C
C----------------------------------------------------------
C     Extra index array needed for F-matrix transformation.
C     Ove Christiansen 17-6-1996
C----------------------------------------------------------
C
      DO 1490 ISYMT = 1,NSYM
         DO 1500 ISYMD = 1,NSYM
            NDSGRH(ISYMD,ISYMT) = 0
            ISYABL = MULD2H(ISYMD,ISYMT)
            DO 1510 ISYMG = 1,NSYM
               ISYMAB = MULD2H(ISYMG,ISYMD)
               ISYML  = MULD2H(ISYABL,ISYMAB)
               NDSGRH(ISYMD,ISYMT) = NDSGRH(ISYMD,ISYMT) 
     *                       + NNBST(ISYMAB)*NRHF(ISYML)
 1510       CONTINUE
 1500    CONTINUE
 1490 CONTINUE
C
      IF (IPRINT .GT. 5) THEN
         CALL AROUND('Information from CCSDSYM')
         WRITE(LUPRI,1) 'NNBST  :',(NNBST(I),  I=1,NSYM)
         WRITE(LUPRI,1) 'NT1AM  :',(NT1AM(I),  I=1,NSYM)
         WRITE(LUPRI,1) 'NT2AM  :',(NT2AM(I),  I=1,NSYM)
         WRITE(LUPRI,1) 'NDISAO :',(NDISAO(I), I=1,NSYM)
         WRITE(LUPRI,1) 'NDSRHF :',(NDSRHF(I), I=1,NSYM)
         WRITE(LUPRI,1) 'ILMRHF :',(ILMRHF(I), I=1,NSYM)
         WRITE(LUPRI,1) 'ILMVIR :',(ILMVIR(I), I=1,NSYM)
         WRITE(LUPRI,1) 'NT1AO  :',(NT1AO(I),  I=1,NSYM)
         WRITE(LUPRI,1) 'NT2AO  :',(NT2AO(I),  I=1,NSYM)
         WRITE(LUPRI,1) 'N2BST  :',(N2BST(I),  I=1,NSYM)
         WRITE(LUPRI,1) 'NT2BCD :',(NT2BCD(I), I=1,NSYM)
         WRITE(LUPRI,1) 'NT2BGD :',(NT2BGD(I), I=1,NSYM)
         WRITE(LUPRI,1) 'NMATIJ :',(NMATIJ(I), I=1,NSYM)
         WRITE(LUPRI,1) 'NGAMMA :',(NGAMMA(I), I=1,NSYM)
         WRITE(LUPRI,1) 'NEMAT1 :',(NEMAT1(I), I=1,NSYM)
         WRITE(LUPRI,1) 'NMATAB :',(NMATAB(I), I=1,NSYM)
         WRITE(LUPRI,1) 'NT2AOS :',(NT2AOS(I), I=1,NSYM)
         WRITE(LUPRI,1) 'NT2SQ  :',(NT2SQ(I) , I=1,NSYM)
         WRITE(LUPRI,1) 'NMIJP  :',(NMIJP(I) , I=1,NSYM)
         WRITE(LUPRI,1) 'NT2ORT :',(NT2ORT(I), I=1,NSYM)
         WRITE(LUPRI,1) 'NGLMDT :',(NGLMDT(I), I=1,NSYM)
         WRITE(LUPRI,1) 'NGLMRH :',(NGLMRH(I), I=1,NSYM)
         WRITE(LUPRI,1) 'NLRHFR :',(NLRHFR(I), I=1,NSYM)
         WRITE(LUPRI,*)
         DO 9901 I = 1,NSYM
            WRITE(LUPRI,1) 'IDSAOG :',(IDSAOG(I,J), J=1,NSYM)
 9901    CONTINUE
         WRITE(LUPRI,*)
         DO 9902 I = 1,NSYM
            WRITE(LUPRI,1) 'IT1AM  :',(IT1AM(I,J), J=1,NSYM)
 9902    CONTINUE
         WRITE(LUPRI,*)
         DO 9903 I = 1,NSYM
            WRITE(LUPRI,1) 'IT2AM  :',(IT2AM(I,J), J=1,NSYM)
 9903    CONTINUE
         WRITE(LUPRI,*)
         DO 9904 I = 1,NSYM
            WRITE(LUPRI,1) 'IT1AO  :',(IT1AO(I,J), J=1,NSYM)
 9904    CONTINUE
         WRITE(LUPRI,*)
         DO 9905 I = 1,NSYM
            WRITE(LUPRI,1) 'IT2AO  :',(IT2AO(I,J), J=1,NSYM)
 9905    CONTINUE
         WRITE(LUPRI,*)
         DO 9906 I = 1,NSYM
            WRITE(LUPRI,1) 'IT2SQ  :',(IT2SQ(I,J), J=1,NSYM)
 9906    CONTINUE
         WRITE(LUPRI,*)
         DO 9907 I = 1,NSYM
            WRITE(LUPRI,1) 'IAODIS :',(IAODIS(I,J), J=1,NSYM)
 9907    CONTINUE
         WRITE(LUPRI,*)
         DO 9908 I = 1,NSYM
            WRITE(LUPRI,1) 'IT2BCD :',(IT2BCD(I,J), J=1,NSYM)
 9908    CONTINUE
         WRITE(LUPRI,*)
         DO 9909 I = 1,NSYM
            WRITE(LUPRI,1) 'IT2BGD :',(IT2BGD(I,J), J=1,NSYM)
 9909    CONTINUE
         WRITE(LUPRI,*)
         DO 9910 I = 1,NSYM
            WRITE(LUPRI,1) 'IMATIJ :',(IMATIJ(I,J), J=1,NSYM)
 9910    CONTINUE
         WRITE(LUPRI,*)
         DO 9911 I = 1,NSYM
            WRITE(LUPRI,1) 'IGAMMA :',(IGAMMA(I,J), J=1,NSYM)
 9911    CONTINUE
         WRITE(LUPRI,*)
         DO 9912 I = 1,NSYM
            WRITE(LUPRI,1) 'IDSRHF :',(IDSRHF(I,J), J=1,NSYM)
 9912    CONTINUE
         WRITE(LUPRI,*)
         DO 9913 I = 1,NSYM
            WRITE(LUPRI,1) 'IT1AMT :',(IT1AMT(I,J), J=1,NSYM)
 9913    CONTINUE
         WRITE(LUPRI,*)
         DO 9914 I = 1,NSYM
            WRITE(LUPRI,1) 'IT1AOT :',(IT1AOT(I,J), J=1,NSYM)
 9914    CONTINUE
         WRITE(LUPRI,*)
         DO 9915 I = 1,NSYM
            WRITE(LUPRI,1) 'IT2BCT :',(IT2BCT(I,J), J=1,NSYM)
 9915    CONTINUE
         WRITE(LUPRI,*)
         DO 9916 I = 1,NSYM
            WRITE(LUPRI,1) 'IT2BGT :',(IT2BGT(I,J), J=1,NSYM)
 9916    CONTINUE
         WRITE(LUPRI,*)
         DO 9917 I = 1,NSYM
            WRITE(LUPRI,1) 'IFCRHF :',(IFCRHF(I,J), J=1,NSYM)
 9917    CONTINUE
         WRITE(LUPRI,*)
         DO 9918 I = 1,NSYM
            WRITE(LUPRI,1) 'IFCVIR :',(IFCVIR(I,J), J=1,NSYM)
 9918    CONTINUE
         WRITE(LUPRI,*)
         DO 9919 I = 1,NSYM
            WRITE(LUPRI,1) 'IEMAT1 :',(IEMAT1(I,J), J=1,NSYM)
 9919    CONTINUE
         WRITE(LUPRI,*)
         DO 9920 I = 1,NSYM
            WRITE(LUPRI,1) 'IMATAB :',(IMATAB(I,J), J=1,NSYM)
 9920    CONTINUE
         WRITE(LUPRI,*)
         DO 9921 I = 1,NSYM
            WRITE(LUPRI,1) 'IT2AOS :',(IT2AOS(I,J), J=1,NSYM)
 9921    CONTINUE
         WRITE(LUPRI,*)
         DO 9922 I = 1,NSYM
            WRITE(LUPRI,1) 'IMIJP  :',(IMIJP(I,J), J=1,NSYM)
 9922    CONTINUE
         WRITE(LUPRI,*)
         DO 9923 I = 1,NSYM
            WRITE(LUPRI,1) 'IT2ORT :',(IT2ORT(I,J), J=1,NSYM)
 9923    CONTINUE
         WRITE(LUPRI,*)
         DO 9924 I = 1,NSYM
            WRITE(LUPRI,1) 'IGLMRH :',(IGLMRH(I,J), J=1,NSYM)
 9924    CONTINUE
         WRITE(LUPRI,*)
         DO 9925 I = 1,NSYM
            WRITE(LUPRI,1) 'IGLMVI :',(IGLMVI(I,J), J=1,NSYM)
 9925    CONTINUE
         WRITE(LUPRI,*)
         DO 9926 I = 1,NSYM
            WRITE(LUPRI,1) 'NT2MMO :',(NT2MMO(I,J), J=1,NSYM)
 9926    CONTINUE
         WRITE(LUPRI,*)
         DO 9927 I = 1,NSYM
            WRITE(LUPRI,1) 'NT2MAO :',(NT2MAO(I,J), J=1,NSYM)
 9927    CONTINUE
         DO 9928 I = 1,NSYM
            WRITE(LUPRI,1) 'NDSGRH :',(NDSGRH(I,J), J=1,NSYM)
 9928    CONTINUE
         WRITE(LUPRI,*)
         WRITE(LUPRI,1) 'NLAMDS :',NLAMDS
         WRITE(LUPRI,1) 'NLRHSI :',NLRHSI
         WRITE(LUPRI,1) 'NLAMDT :',NLAMDT
         WRITE(LUPRI,1) 'NLMRHF :',NLMRHF
      END IF
C
      RETURN
C
    1 FORMAT(3X,A8,8I8)
C
      END
C  /* Deck fock_reorder */
      SUBROUTINE FOCK_REORDER(FOCK,WORK,LWORK)
C
C     Henrik Koch and Alfredo Sanchez.       29-Jun-1994
C
C     Reorder the symmetry ordering of the fock matrix.
C     First occupied orbitals in different symmetries and then
C     the virtuals in different symmetries.
C
#include <implicit.h>
      DIMENSION FOCK(NORBT),WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      IF (LWORK .LT. NORBT) THEN
         CALL QUIT('Insufficient space in FOCK_REORDER')
      ENDIF
C
      ICRHF  = 0
      ICVIR  = NRHFT
      ICOUNT = 0
      DO 100 ISYM = 1,NSYM
C
         DO 110 I = 1,NRHF(ISYM)
            ICRHF  = ICRHF  + 1
            ICOUNT = ICOUNT + 1
            WORK(ICRHF) = FOCK(ICOUNT)
  110    CONTINUE
C
         DO 120 A = 1,NVIR(ISYM)
            ICVIR  = ICVIR  + 1
            ICOUNT = ICOUNT + 1
            WORK(ICVIR) = FOCK(ICOUNT)
  120    CONTINUE
C
  100 CONTINUE
C
      IF (IPRINT .GT. 20) THEN
         CALL AROUND('Fock matrix diagonal in FOCK_REORDER')
         WRITE(LUPRI,1)
         DO 200 I = 1,NORBT
            WRITE(LUPRI,2) FOCK(I),WORK(I)
  200    CONTINUE
      END IF
C
      CALL DCOPY(NORBT,WORK,1,FOCK,1)
C
      RETURN
C
    1 FORMAT(7X,'Sirius order',5X,'CCSD order')
    2 FORMAT(6X,F14.10,3X,F14.10)
C
      END
C  /* Deck cmo_reorder */
      SUBROUTINE CMO_REORDER(CMO,WORK,LWORK)
C
C     Henrik Koch and Alfredo Sanchez.       30-Jun-1994
C
C     Reorder the symmetry ordering of the MO coefficient matrix.
C     First occupied orbitals in different symmetries and then
C     the virtuals in different symmetries.
C
#include <implicit.h>
      DIMENSION CMO(*),WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
C-----------------------
C     Memory allocation.
C-----------------------
C
      KSCR1 = 1
      KSCR2 = KSCR1 + NLAMDS
      KEND  = KSCR2 + NLAMDT
      LWRK1 = LWORK - KEND
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CMO_REORDER')
      ENDIF
C
C----------------------------------
C     Reorder all orbitals in work.
C----------------------------------
C
      ICRHF  = KSCR1
      ICVIR  = KSCR1 + NLRHSI
      ICOUNT = 1
      DO 100 ISYM = 1,NSYM
C
         CALL DCOPY(NBAS(ISYM)*NRHFS(ISYM),CMO(ICOUNT),1,WORK(ICRHF),1)
         ICRHF  = ICRHF  + NBAS(ISYM)*NRHFS(ISYM)
         ICOUNT = ICOUNT + NBAS(ISYM)*NRHFS(ISYM)
C
         CALL DCOPY(NBAS(ISYM)*NVIRS(ISYM),CMO(ICOUNT),1,WORK(ICVIR),1)
         ICVIR  = ICVIR  + NBAS(ISYM)*NVIRS(ISYM)
         ICOUNT = ICOUNT + NBAS(ISYM)*NVIRS(ISYM)
C
  100 CONTINUE
C
C----------------------------
C     Delete frozen orbitals.
C----------------------------
C
      IF ((.NOT. FROIMP) .AND. (.NOT. FROEXP)) THEN
C
         CALL DCOPY(NLAMDT,WORK(KSCR1),1,WORK(KSCR2),1)
C
      ELSE IF (FROIMP) THEN
C
         DO 110 ISYM = 1, NSYM
C
            KOFF1 = KSCR1 + ILRHSI(ISYM) + NBAS(ISYM)*NRHFFR(ISYM)
            KOFF2 = KSCR2 + ILMRHF(ISYM)
C
            LENGTH = NBAS(ISYM)*NRHF(ISYM)
            CALL DCOPY(LENGTH,WORK(KOFF1),1,WORK(KOFF2),1)
C
            KOFF1 = KSCR1 + ILVISI(ISYM)
            KOFF2 = KSCR2 + ILMVIR(ISYM)
C
            LENGTH = NBAS(ISYM)*NVIR(ISYM)
            CALL DCOPY(LENGTH,WORK(KOFF1),1,WORK(KOFF2),1)
C
  110    CONTINUE
C
      ELSE
C
         DO 120 ISYM = 1,NSYM
C
             KOFF1 = KSCR1 + ILRHSI(ISYM)
             KOFF2 = KSCR2 + ILMRHF(ISYM)
C
             DO 130 IOCC = 1,NRHFS(ISYM)
C
                IF (.NOT. FRORHF(IOCC,ISYM)) THEN
                   CALL DCOPY(NBAS(ISYM),WORK(KOFF1),1,WORK(KOFF2),1)
                   KOFF2 = KOFF2 + NBAS(ISYM)
                END IF
C
                KOFF1 = KOFF1 + NBAS(ISYM)
C
  130        CONTINUE
C
             KOFF1 = KSCR1 + ILVISI(ISYM)
             KOFF2 = KSCR2 + ILMVIR(ISYM)
C
             DO 140 IVIR1 = 1,NVIRS(ISYM)
C
                IF (.NOT. FROVIR(IVIR1,ISYM)) THEN
                   CALL DCOPY(NBAS(ISYM),WORK(KOFF1),1,WORK(KOFF2),1)
                   KOFF2 = KOFF2 + NBAS(ISYM)
                END IF
C
                KOFF1 = KOFF1 + NBAS(ISYM)
C
  140        CONTINUE
C
  120    CONTINUE
C
      END IF
C
C----------------------
C     Print if desired.
C----------------------
C
      IF (IPRINT .GT. 20) THEN
         CALL AROUND('MO-coefficient matrix in CMO_REORDER')
         KOFF1 = 1
         KOFF2 = KSCR2
         KOFF3 = KSCR2 + NLMRHF
         DO 200 ISYM = 1,NSYM
            WRITE(LUPRI,1) ISYM
            IF (NORB(ISYM) .EQ. 0) THEN
               WRITE(LUPRI,8)
               GOTO 200
            ENDIF
            WRITE(LUPRI,2)
            WRITE(LUPRI,3)
            CALL OUTPUT(CMO(KOFF1),1,NBAS(ISYM),1,NORBS(ISYM),
     *                  NBAS(ISYM),NORBS(ISYM),1,LUPRI)
            WRITE(LUPRI,4)
            WRITE(LUPRI,5)
            CALL OUTPUT(WORK(KOFF2),1,NBAS(ISYM),1,NRHF(ISYM),
     *                  NBAS(ISYM),NRHF(ISYM),1,LUPRI)
            WRITE(LUPRI,6)
            WRITE(LUPRI,7)
            CALL OUTPUT(WORK(KOFF3),1,NBAS(ISYM),1,NVIR(ISYM),
     *                  NBAS(ISYM),NVIR(ISYM),1,LUPRI)
            KOFF1 = KOFF1 + NBAS(ISYM)*NORBS(ISYM)
            KOFF2 = KOFF2 + NBAS(ISYM)*NRHF(ISYM)
            KOFF3 = KOFF3 + NBAS(ISYM)*NVIR(ISYM)
  200    CONTINUE
      END IF
C
      CALL DCOPY(NLAMDT,WORK(KSCR2),1,CMO,1)
C
      RETURN
C
    1 FORMAT(/,/,7X,'Symmetry number :',I5)
    2 FORMAT(/,/,7X,'Sirius ordering')
    3 FORMAT(7X,'---------------')
    4 FORMAT(/,/,7X,'CCSD ordering occupied part')
    5 FORMAT(7X,'---------------------------')
    6 FORMAT(/,/,7X,'CCSD ordering virtual part')
    7 FORMAT(7X,'--------------------------')
    8 FORMAT(/,/,7X,'This symmetry is empty')
C
      END
C  /* Deck ccsd_symsqo */
      SUBROUTINE CCSD_SYMSQO(DISTAB,ISYMAB,SCR)
C
C     Henrik Koch and Alfredo Sanchez.       1-July-1994
C
C     Squareup the integral distribution.
C
#include <implicit.h>
      DIMENSION DISTAB(*), SCR(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      IF (ISYMAB .EQ. 1) THEN
C
         KOFF1 = 1
         KOFF2 = 1
         DO 100 ISYMB = 1,NSYM
            CALL SQMATR(NBAS(ISYMB),DISTAB(KOFF1),SCR(KOFF2))
            KOFF1 = KOFF1 + NBAS(ISYMB)*(NBAS(ISYMB)+1)/2
            KOFF2 = KOFF2 + NBAS(ISYMB)*NBAS(ISYMB)
  100    CONTINUE
C
      ELSE
         KOFF1 = 1
         KOFF2 = 1
         DO 200 ISYMB = 1,NSYM
C
            ISYMA = MULD2H(ISYMB,ISYMAB)
            IF (ISYMB .GT. ISYMA) THEN
C
               NTOT  = NBAS(ISYMA)*NBAS(ISYMB)
C
               KOFF2 = KOFF1
               KOFF3 = IAODIS(ISYMB,ISYMA) + 1
               DO 210 B = 1,NBAS(ISYMB)
                  CALL DCOPY(NBAS(ISYMA),DISTAB(KOFF2),1,SCR(KOFF3),
     *                       NBAS(ISYMB))
                  KOFF2 = KOFF2 + NBAS(ISYMA)
                  KOFF3 = KOFF3 + 1
  210          CONTINUE
C
               KOFF4 = IAODIS(ISYMA,ISYMB) + 1
               CALL DCOPY(NTOT,DISTAB(KOFF1),1,SCR(KOFF4),1)
C
               KOFF1 = KOFF1 + NTOT
C
            ENDIF
C
  200    CONTINUE
C
      ENDIF
C
      RETURN
      END
      SUBROUTINE CCSD_SYMSQ(DISTAB,ISYMAB,SCR)
C
C     Henrik Koch and Alfredo Sanchez.       1-July-1994
C
C     Squareup the integral distribution.
C
#include <implicit.h>
      DIMENSION DISTAB(*), SCR(*)
#include <priunit.h>
#include <maxorb.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <symsq.h>
C
C
#if !defined (SYS_CRAY)
      DO 100 IJSQ = 1,N2BST(ISYMAB)
C
         KOFF = I2BST(ISYMAB) + IJSQ
         IJPK = IADRPK(KOFF)
C
         SCR(IJSQ) = DISTAB(IJPK)
C
  100 CONTINUE
#endif
C
#if defined (SYS_CRAY)
      KOFF = I2BST(ISYMAB) + 1
      CALL GATHER(N2BST(ISYMAB),SCR,DISTAB,IADRPK(KOFF))
#endif
C
      RETURN
      END
C  /* Deck cc3_sort1 */
      SUBROUTINE CC3_SORT1(WORK,LWORK,IOPT,ISYINT)
C
C     Henrik Koch and Alfredo Sanchez.       28-May-1995
C
C     Sort virtual integrals for perturbative triples.
C
#include <implicit.h>
      DIMENSION WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccinftap.h>
#include <ccsdsym.h>
C
C
      IF ((IOPT .NE. 1) .AND. (IOPT .NE. 2)) THEN
         CALL QUIT('IOPT error in CC3_SORT1')
      END IF
C
C-----------------------------------------
C     Start loop over symmetries of delta.
C-----------------------------------------
C
      MAXCK = 0
      DO 50 ISYMCK = 1,NSYM
         IF (NT1AM(ISYMCK) .GT. MAXCK) MAXCK = NT1AM(ISYMCK)
   50 CONTINUE
C
      DO 100 ISYMD = 1,NSYM
C
         IF (NBAS(ISYMD) .EQ. 0) GOTO 100
C
C--------------------------
C        Memory allocation.
C--------------------------
C
         ISYCKB = MULD2H(ISYMD,ISYINT)
C
         LENMIN = NCKATR(ISYCKB) + MAXCK
         NDISTR = MIN(LWORK/LENMIN,NBAS(ISYMD))
C
         IF (NDISTR .EQ. 0) THEN
            CALL QUIT('Insufficient work space in CC3_SORT1')
         ENDIF
C
         NBATCH = (NBAS(ISYMD) - 1)/NDISTR + 1
C
         KSCR1 = 1
         KSCR2 = KSCR1 + NCKATR(ISYCKB)*NDISTR
         KEND1 = KSCR2 + MAXCK*NDISTR
C
         DO 110 IBATCH = 1,NBATCH
C
            NUMD = NDISTR
            IF (IBATCH .EQ. NBATCH) THEN
               NUMD = NBAS(ISYMD) - NDISTR*(NBATCH - 1)
            ENDIF
C
            ID1 = NDISTR*(IBATCH - 1) + 1
C
C--------------------------
C           Read integrals.
C--------------------------
C
            LENGTH = NCKATR(ISYCKB)*NUMD
C
            IOFF = ICKDAO(ISYCKB,ISYMD) + NCKATR(ISYCKB)*(ID1 - 1) + 1
C
            IF (LENGTH .GT. 0) THEN
               CALL GETWA2(LU3SRT,FN3SRT,WORK(KSCR1),IOFF,LENGTH)
            ENDIF
C
C-----------------------------------------------------
C           Sort integrals (ck,del,b) from (ck,b,del).
C-----------------------------------------------------
C
C
            DO 120 ISYMB = 1,NSYM
C
               ISYMCK = MULD2H(ISYCKB,ISYMB)
               ISYCKD = MULD2H(ISYMCK,ISYMD)
C
               DO 130 B = 1,NVIR(ISYMB)
C
                  DO 140 I = 1,NUMD
C
                     ID = ID1 + I - 1
C
                     KOFF1 = KSCR1
     *                     + NCKATR(ISYCKB)*(I - 1)
     *                     + ICKATR(ISYMCK,ISYMB)
     *                     + NT1AM(ISYMCK)*(B - 1)
                     KOFF2 = KSCR2
     *                     + NT1AM(ISYMCK)*(I - 1)

                     CALL DCOPY(NT1AM(ISYMCK),WORK(KOFF1),1,
     *                          WORK(KOFF2),1)
C
  140             CONTINUE
C
C----------------------------------------
C                 Write sorted integrals.
C----------------------------------------
C
                  LENGTH = NT1AM(ISYMCK)*NUMD
C
                  IF (LENGTH .GT. 0) THEN
C
                     IOFF = ICKAD(ISYCKD,ISYMB)
     *                    + NCKA(ISYCKD)*(B - 1)
     *                    + ICKA(ISYMCK,ISYMD)
     *                    + NT1AM(ISYMCK)*(ID1 - 1) + 1
C
                     IF (IOPT .EQ. 1) THEN
                        CALL PUTWA2(LU3VI,FN3VI,WORK(KSCR2),IOFF,LENGTH)
                     ELSE
                        CALL PUTWA2(LUDELD,FNDELD,WORK(KSCR2),IOFF,
     &                              LENGTH)
                     END IF
                  ENDIF
C
  130          CONTINUE
  120       CONTINUE
C
            IF (IOPT .EQ. 2) GOTO 110
C
C-----------------------------------------------------
C           Sort integrals (bk,del,c) from (ck,b,del).
C-----------------------------------------------------
C
            DO 150 ISYMC = 1,NSYM
C
               ISYMBK = MULD2H(ISYCKB,ISYMC)
               ISYBKD = MULD2H(ISYMBK,ISYMD)
C
               DO 160 C = 1,NVIR(ISYMC)
C
                  DO 170 I = 1,NUMD
C
                     ID = ID1 + I - 1
C
                     DO 180 ISYMK = 1,NSYM
C
                        ISYMB  = MULD2H(ISYMBK,ISYMK)
                        ISYMCK = MULD2H(ISYMC,ISYMK)
C
                        NTOTCK = MAX(NT1AM(ISYMCK),1)
C
                        DO 190 K = 1,NRHF(ISYMK)

C
                           KOFF1 = KSCR1
     *                           + NCKATR(ISYCKB)*(I - 1)
     *                           + ICKATR(ISYMCK,ISYMB)
     *                           + IT1AM(ISYMC,ISYMK)
     *                           + NVIR(ISYMC)*(K - 1) + C - 1
C
                           KOFF2 = KSCR2
     *                           + NT1AM(ISYMBK)*(I - 1)
     *                           + IT1AM(ISYMB,ISYMK)
     *                           + NVIR(ISYMB)*(K - 1)
C
                           CALL DCOPY(NVIR(ISYMB),WORK(KOFF1),NTOTCK,
     *                                WORK(KOFF2),1)
C
  190                   CONTINUE
  180                CONTINUE
  170             CONTINUE
C
C----------------------------------------
C                 Write sorted integrals.
C----------------------------------------
C
                  LENGTH = NT1AM(ISYMBK)*NUMD
C
                  IF (LENGTH .GT. 0) THEN
C
                     IOFF = ICKAD(ISYBKD,ISYMC)
     *                    + NCKA(ISYBKD)*(C - 1)
     *                    + ICKA(ISYMBK,ISYMD)
     *                    + NT1AM(ISYMBK)*(ID1 - 1) + 1
C
                     CALL PUTWA2(LU3VI2,FN3VI2,WORK(KSCR2),IOFF,LENGTH)
                  ENDIF
C
  160          CONTINUE
  150       CONTINUE
C
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck ccsd_delfro */
      SUBROUTINE CCSD_DELFRO(FOCDIA,WORK,LWORK)
C
#include <implicit.h>
C
      DIMENSION FOCDIA(*),WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      IF (LWORK .LT. NORBT) THEN
         WRITE(LUPRI,*) 'Insufficient space in CCSD_DELFRO'
         CALL QUIT(' ')
      END IF
C
      KOFF1 = 0
      KOFF2 = 0
C
      DO 100 ISYM = 1,NSYM
C
         DO 110 I = 1,NRHFS(ISYM)
            KOFF1 = KOFF1 + 1
            IF (.NOT. FRORHF(I,ISYM)) THEN
               KOFF2 = KOFF2 + 1
               WORK(KOFF2) = FOCDIA(KOFF1)
            END IF
  110    CONTINUE
C
         DO 120 A = 1,NVIRS(ISYM)
            KOFF1 = KOFF1 + 1
            IF (.NOT. FROVIR(A,ISYM)) THEN
               KOFF2 = KOFF2 + 1
               WORK(KOFF2) = FOCDIA(KOFF1)
            END IF
  120    CONTINUE
C
  100 CONTINUE
C
      CALL DCOPY(NORBT,WORK,1,FOCDIA,1)
C
      RETURN
      END
C  /* Deck CC_freezer*/
      SUBROUTINE CC_FREEZER(FOCDIA,NF,NFCS,NFVS,WORK,LWORK)
C
C     Ove Christiansen 230899, Find and freeze NFC lowest-lying/NFV highest lying
C     canonical orbitals in CC calculation.
C
#include <implicit.h>
#include <mxorb.h>
C
      DIMENSION FOCDIA(NF),WORK(LWORK),NFCS(8),NFVS(8),IPLACE(MXCORB)
#include <ccorb.h>
#include <priunit.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
C
      IF (IPRINT.GT.5) WRITE(LUPRI,*) ' In CC_FREEZER: ' 
      IF (IPRINT.GT.5) WRITE(LUPRI,*) ' Freezing occupied, virtual:',NFC,NFV
      IF (LWORK .LT. NORBT) THEN
         WRITE(LUPRI,*) 'Insufficient space in CCSD_DELFRO'
         CALL QUIT( 'Insufficient space in CCSD_DELFRO')
      END IF
      CALL FLSHFO(LUPRI)
C
      DO ISYM=1,NSYM
         NRHFFR(ISYM) = 0
         NVIRFR(ISYM) = 0
      ENDDO
C
C-----------------------------------------------------------------------
C     Find NFC lowest orbital energies
C-----------------------------------------------------------------------
C
      IF (NFC.GT.0) THEN
       MXELMN = NFC
       NELMN  = NFC
       THRDIA = 1.0D-06
       CALL FNDMN3(FOCDIA,NF,MXELMN,IPLACE,
     *            NELMN,IPRINT,THRDIA)
C
C-------------------------------------------------------------
C        Find # frozen orbitals in each symmetryclass: NRHFFR 
C-------------------------------------------------------------
C
       DO I = 1,NFC 
         IHFO = IPLACE(I)
         CALL CC_SYMHFO(IHFO,ISYMHFO)
         WRITE(LUPRI,'(A,I3,A,I3,A,F10.4)') 
     *    ' Freezing HF-orbital ',IHFO,' of symmetry '
     *    ,ISYMHFO,' and with orbital energy',FOCDIA(IHFO)
         NRHFFR(ISYMHFO) = NRHFFR(ISYMHFO)+1
       ENDDO
       WRITE(LUPRI,'(A,8I3)')  
     *      ' In total frozen-core per symmetry-class:', 
     *                     (NRHFFR(ISYM),ISYM=1,NSYM)
       WRITE(LUPRI,'(A)')  ' ' 
      ENDIF
C
C-----------------------------------------------------------------------
C     Find NFV highest orbital energies
C-----------------------------------------------------------------------
C
      IF (NFV.GT.0) THEN
       MXELMN = NFV
       NELMN  = NFV
       THRDIA = 1.0D-06
       ONEM = -1.0D0
       CALL DSCAL(NF,ONEM,FOCDIA,1)
       CALL FNDMN3(FOCDIA,NF,MXELMN,IPLACE,
     *             NELMN,IPRINT,THRDIA)
       CALL DSCAL(NF,ONEM,FOCDIA,1)
C
C-------------------------------------------------------------
C        Find # frozen orbitals in each symmetryclass: NVIRFR 
C-------------------------------------------------------------
C
C
       DO I = 1,NFV
         IHFO = IPLACE(I)
         CALL CC_SYMHFO(IHFO,ISYMHFO)
         WRITE(LUPRI,'(A,I3,A,I3,A,F10.4)')
     *    ' Freezing HF-orbital ',IHFO,' of symmetry '
     *    ,ISYMHFO,' and with orbital energy',FOCDIA(IHFO)
         NVIRFR(ISYMHFO) = NVIRFR(ISYMHFO)+1
       ENDDO
       WRITE(LUPRI,'(A,8I3)')
     *         ' In total frozen-virtual per symmetry-class:', 
     *                     (NVIRFR(ISYM),ISYM=1,NSYM)
       WRITE(LUPRI,'(A)')  ' ' 
      ENDIF
C
C-----------------------------------------------------------------------
c     Put orbitals lowest and highest obital energies on the list of 
C     orbitals to be deleted.
C-----------------------------------------------------------------------
C
      DO ISYM = 1,NSYM

         IF (NRHFFR(ISYM) .GT. MAXFRO) THEN
            WRITE(LUPRI,'(//,1X,2A,I3)') 'ERROR: Maximum number of ',
     *           'frozen orbitals per symmetry is:',MAXFRO
            CALL QUIT('Too many frozen orbitals')
         END IF
         DO I = 1,NRHFFR(ISYM)
            FRORHF(I,ISYM) = .TRUE.
            KFRRHF(I,ISYM) = I
         ENDDO

         IF (NVIRFR(ISYM) .GT. MAXFRO) THEN
            WRITE(LUPRI,'(//,1X,2A,I3)') 'ERROR: Maximum number of ',
     *           'frozen orbitals per symmetry is:',MAXFRO
            CALL QUIT('Too many frozen orbitals')
         END IF
         DO I = 1,NVIRFR(ISYM)
            JORB = NVIRS(ISYM) - I + 1
            FROVIR(JORB,ISYM) = .TRUE.
            KFRVIR(I,ISYM)    = JORB
         ENDDO

      ENDDO    
      CALL FLSHFO(LUPRI)
      RETURN
      END
C  /* Deck CC_SYMHFO*/
      SUBROUTINE CC_SYMHFO(IHFO,ISYMHFO)
C
C     OC 230899, find symmetry ISYMHFO of HF orbital nr. IHFO
C
#include <implicit.h>
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
C
      ISYMHFO = 0
      ICOUNT = 0
      DO ISYM = 1, NSYM
        IF ((IHFO.GT.ICOUNT).AND.(IHFO.LE.(ICOUNT+NORBS(ISYM))))
     *    ISYMHFO = ISYM
          ICOUNT = ICOUNT + NORBS(ISYM)
      ENDDO
      IF (ISYMHFO.EQ.0) WRITE(LUPRI,*) 'Something is wrong in CC_SYMHFO'
      RETURN
      END
