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
C  /* Deck ccsd_energy */
      SUBROUTINE CCSD_ENERGY(WORK,LWORK,APROXR12,CCR12RSP,CCR12LIM)
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
      LOGICAL   CCSAV, CC1BSV, CC1ASV, CCPTSV, CCP3SV, LCCEQ
      LOGICAL   MP2R12TST,LRES
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
#include <gnrinf.h>
#include <ccfdgeo.h>
#include <cbirea.h>
#include <r12int.h>
#include <ccr12int.h>
!Sonia
#include <ccfop.h>
#include <ccnoddy.h>
Cholesky
#include <ccdeco.h>
#include <chodbg.h>
#include <cc_cho.h>
#include <chocc2.h>
C
      LOGICAL CPTDBG
Cholesky
      LOGICAL LCONVG,RSPIM2,EX,LEXIST,LHTF,MKVABKL
      LOGICAL CCR12RSP, CCR12LIM
      LOGICAL LCONV1,LCONV2
      CHARACTER*5 ETY0, ETY1, ETY2
      CHARACTER MODEL*10, MODELR*10, ETYPE*24, MODELR12*24, MOPRPC*10
      CHARACTER MODREF*10
      CHARACTER*3 APROXR12
      CHARACTER*24 BLANKS
      DATA BLANKS /'                        '/
      INTEGER LENMOD
      CHARACTER*8 LABEL1
C
      CALL QENTER('CCSD_ENERGY')
celena
      IF (R12PRP) INTTR = .TRUE.
celena
C
C     -------------------------------------------------------------
C     set model for which the current t-amplitudes were calculated:
C     -------------------------------------------------------------
C
      MODEL = 'UNKNOWN   '
      IF (CIS)   MODEL = 'CIS       '
      IF (CCS)   MODEL = 'CCS       '
      IF (MP2)   MODEL = 'MP2       '
      IF (CC2)   MODEL = 'CC2       '
      IF (CCD)   MODEL = 'CCD       '
      IF (CCSD)  MODEL = 'CCSD      '
      IF (CC3)   MODEL = 'CC3       '
      IF (CC1A)  MODEL = 'CCSDT-1a  '
      IF (CC1B)  MODEL = 'CCSDT-1b  '
      IF (CCPT)  MODEL = 'CCSD(T)   '
      IF (CCP3)  MODEL = 'CC(3)     '
      IF (CCRT)  MODEL = 'CCSDR(T)  '
      IF (CCR3)  MODEL = 'CCSDR(3)  '
      IF (CCR1A) MODEL = 'CCSDR(1A) '
      IF (CCR1B) MODEL = 'CCSDR(1B) '
      MOPRPC = MODEL
      ! set model for CCR12
      CALL CCSD_MODEL(MODELR12,LENMOD,24,MODEL,10,APROXR12)
      MODEL = MODELR12(1:10)
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
      EX     = .FALSE.
C
C-------------------------------------------------
C     Employ MP2-R12 method (WK/UniKA/04-11-2002).
C-------------------------------------------------
C
      R12NOP = R12NOP .OR. .NOT. R12XXL
      IF (R12CAL.AND..NOT.LISKIP) THEN
        IPRSAVE = IPRINT
        IPRINT  = IPRINT / 10
        CALL GETTIM(T0,W0)
        CALL CCSD_R12(WORK,LWORK,WORK,LWORK,CCR12RSP)
        CALL GETTIM(T1,W1)
        WRITE(LUPRI,*)'Time for MP2-R12 part cpu :', T1-T0
        WRITE(LUPRI,*)'Time for MP2-R12 part wall:', W1-W0
        CALL FLSHFO(LUPRI)
        ! restore print level
        IPRINT = IPRSAVE
C
C       Use LABEL (WK/UniKA/04-11-2002).
        LABEL = 'TRCCINT '
        IF (LMULBS) THEN
          NOAUXB = .TRUE.
          IF (HERDIR) THEN
            WRITE (LUPRI,'(/A/)') 'NOAUXB with HERDIR not implemented'
            GOTO 9999
          ENDIF
C         IF (.NOT. DIRECT) 
C    &               CALL QUIT('NOAUXB without DIRECT not implemented')
        END IF
C
C       reset nbas, etc. to original values:
        CALL CCSD_INIT1(WORK,LWORK)
      END IF
C
C     switch off R12-MP12 for future calls
      R12CAL = .FALSE.
      CC2R12INT = .FALSE.
      CCSDR12INT= .FALSE.
C
C     use V^(alpha beta)_(kl)?
      USEVABKL = CCR12 .AND. (USEVABKL .OR. .NOT.CC2) 
      MKVABKL  = USEVABKL
C
      IF (CCR12.AND.MP2 .AND. .NOT. R12PRP) THEN
         CALL QEXIT('CCSD_ENERGY')
         RETURN
      ELSE IF (CCR12.AND.(.NOT.LISKIP)) THEN
        IF (MKVABKL .AND. .NOT. MP2) THEN
          WRITE(LUPRI,*)'Preparing R12 V-interm. ... ONEAUX=',ONEAUX
          CALL GETTIM(T0,W0)
          CALL CC_R12PREPCCSD(WORK,LWORK)
          CALL GETTIM(T1,W1)
          WRITE(LUPRI,*)'Time used for V^albe_kl cpu:', T1-T0
          WRITE(LUPRI,*)'Time used for V^albe_kl wall:',W1-W0
          WRITE(LUPRI,*)
        END IF

        IF (MP2 .OR. (IANR12.EQ.1)) THEN
          CONTINUE
        ELSE IF (IANR12.EQ.2 .OR. IANR12.EQ.3) THEN
          IF (.NOT.CC2 .AND. .NOT.R12CBS)
     *      CALL QUIT('This CC-R12 model is not implemented w/o CABS')
          WRITE(LUPRI,*)'Preparing R12 Ansatz 2/3 ... ONEAUX=',ONEAUX
          CALL CCR12PREP2(WORK,LWORK)
        ELSE
          WRITE(LUPRI,*) 'IANR12 = ',IANR12
          CALL QUIT('This CC-R12 Ansatz is currently not implemented')
        END IF
      END IF

C     ----------------------------------------------------------------
C     Read packed r12 amplitudes from file CCR12_D for present
C     Ansatz and approximation and put on CCR12_C and CCR0_1___1
C     ----------------------------------------------------------------
      IOPT = 0
      IF (CCR12.AND..NOT.(CIS.OR.CCS.OR.MP2)) THEN
           KTAMP12 = 1
           KEND1   = KTAMP12 + NTR12AM(1)
           LWRK1   = LWORK - KEND1
           IF (LWRK1 .LT. 0) THEN
             CALL QUIT('Not enough work space for R12')
           END IF
           LU43  = -43
           CALL GPOPEN(LU43,FCCR12D,'UNKNOWN',' ','UNFORMATTED',
     &                 IDUM,LDUM)
 1816      READ(LU43,end=1817) IAN,IAP,APROXR12
           READ(LU43) (WORK(KTAMP12-1+I),I=1,NTR12AM(1))
           IF ((IAN.NE.IANR12).OR.(IAP.NE.IAPR12)) GOTO 1816
           CALL GPCLOSE(LU43,'KEEP') 
           CALL GPOPEN(LU43,FCCR12C,'UNKNOWN',' ','UNFORMATTED',
     &                 IDUM,LDUM)
           WRITE(LU43) (WORK(KTAMP12-1+I),I=1,NTR12AM(1))
           CALL GPCLOSE(LU43,'KEEP') 
C          
           IF (.NOT.CCRSTR) THEN
C            WRITE(LUPRI,*) 'Writing R12 amplitudes to disk, MODEL=',MODEL
             IOPT = 32
             CALL CC_WRRSP('R0 ',0,1,IOPT,MODEL,DUMMY,DUMMY,
     &                     WORK(KTAMP12),WORK(KEND1),LWRK1) 
           END IF
      END IF
C
C----------------------------------------------------------------------------
C     Calculate X^V intermediates needed for CCR12 response and finite fields
C----------------------------------------------------------------------------
C
C     IF (CCR12RSP) THEN
C       call cc_r12vxint(work,lwork,.false.)
C     END IF
C
C-------------------
C     Cholesky debug
C-------------------
C
      IF (CHODBG) CALL CC_CHODBG(WORK,LWORK)
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
      WRITE(LUPRI,*)
      WRITE(LUPRI,'(13X,A)')
     *     'The Direct Coupled Cluster Energy Program'
      WRITE(LUPRI,'(13X,A)')
     *     '-----------------------------------------'
      WRITE(LUPRI,'(//10X,A,I8)')
     *     'Number of t1 amplitudes                 :  ',NT1AMX
      WRITE(LUPRI,'(10X,A,I8)')
     *     'Number of t2 amplitudes                 :  ',NT2AMX
      NCCVAR = NT1AMX + NT2AMX
      IF (CCR12) THEN
        WRITE(LUPRI,'(10X,A,I8)')
     *     "Number of t2' amplitudes for R12 part   :  ",NTR12AM(1)
        NCCVAR = NCCVAR + NTR12AM(1)
      END IF
      WRITE(LUPRI,'(10X,A,I8)')
     *     'Total number of amplitudes in ccsd      :  ',NCCVAR
      WRITE(LUPRI,*)
      CALL FLSHFO(LUPRI)
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
         LABEL1 = 'ENERGY  '
         MODREF = 'CCS/SCF   '
         CALL CC_PRPC(ESCF,MODREF,0,
     *                LABEL1,LABEL1,LABEL1,LABEL1,
     *                ZERO,ZERO,ZERO,1,0,0,0)
         CALL WRIPRO(ESCF,MODREF,0,
     *               LABEL1,LABEL1,LABEL1,LABEL1,
     *               ZERO,ZERO,ZERO,1,0,0,0)

         GO TO 9999
C        ... exit this routine
      ENDIF
C
C     
C--------------------
C     Cholesky stuff.
C--------------------
C        
C     Cholesky MP2 section.
C     ---------------------
         
      IF (MP2 .AND. CHOINT) THEN
      
         CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',
     &               IDUMMY,.FALSE.) 
         REWIND LUSIFC
         CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
         READ (LUSIFC) POTNUC,EMY,EACTIV,EMCSCF,ISTATE,ISPIN,NACTEL,
     &              LSYM
         ESCF = EMCSCF
         CALL GPCLOSE(LUSIFC,'KEEP')
         
C        Calculate MP2 energy correction.
C        --------------------------------

         CALL CC_CHOMP2(WORK,LWORK,EMP2)
C        
C        Write SCF and MP2 energies to output and summary.
C        -------------------------------------------------

         ETY1   = 'MP2  '
         ECCGRS = ESCF + EMP2
         
         CALL AROUND
     &   ('Final results from the Coupled Cluster energy program')
     
         WRITE(LUPRI,'(/12X,A,A,A,F32.10)')
     &               'Total ',ETY0,' energy: ',ESCF
         WRITE(LUPRI,'(12X,A,A,A,F32.10)')
     &               'Total ',ETY1,' energy: ',ECCGRS
         WRITE(LURES,'(/12X,A,A,A,F32.10)')
     &               'Total ',ETY0,' energy: ',ESCF
         WRITE(LURES,'(12X,A,A,A,F32.10)')
     &               'Total ',ETY1,' energy: ',ECCGRS

         LABEL1 = 'ENERGY  '
         MODREF = 'MP2/CHOLES'
         CALL CC_PRPC(ECCGRS,MODREF,0,
     *                LABEL1,LABEL1,LABEL1,LABEL1,
     *                ZERO,ZERO,ZERO,1,0,0,0)
         CALL WRIPRO(ECCGRS,MODREF,0,
     *               LABEL1,LABEL1,LABEL1,LABEL1,
     *               ZERO,ZERO,ZERO,1,0,0,0)

         GOTO 9999

      ENDIF
C
C     Cholesky CC2 section.
C     ---------------------

      IF (CC2 .AND. CHOINT) THEN

         CALL CC_CHOECC2(WORK,LWORK,ESCF,ECC2,RSPIM2)

         ETY1   = 'CC2  '
         ECCGRS = ECC2

         CALL AROUND
     &   ('Final results from the Coupled Cluster energy program')

         WRITE(LUPRI,'(/12X,A,A,A,F32.10)')
     &               'Total ',ETY0,' energy: ',ESCF
         WRITE(LUPRI,'(12X,A,A,A,F32.10)')
     &               'Total ',ETY1,' energy: ',ECCGRS
         WRITE(LURES,'(/12X,A,A,A,F32.10)')
     &               'Total ',ETY0,' energy: ',ESCF
         WRITE(LURES,'(12X,A,A,A,F32.10)')
     &               'Total ',ETY1,' energy: ',ECCGRS


         LABEL1 = 'ENERGY  '
         MODREF = 'CC2/CHOLES'
         CALL CC_PRPC(ECCGRS,MODREF,0,
     *                LABEL1,LABEL1,LABEL1,LABEL1,
     *                ZERO,ZERO,ZERO,1,0,0,0)
         CALL WRIPRO(ECCGRS,MODREF,0,
     *               LABEL1,LABEL1,LABEL1,LABEL1,
     *               ZERO,ZERO,ZERO,1,0,0,0)


         GOTO 9999

      ENDIF
C
C------------------------
C     Dynamic allocation.
C------------------------
C
      NTAMR12 = 0 
      IF (CCD) THEN
        NTAMP = NT2AMX
      ELSE
        NTAMP = NT1AMX  + NT2AMX
      END IF
      IF (LMULBS) THEN
C       add length of R12 part
        NTAMR12 = NTR12AM(1)
        NTAMP = NTAMP + NTAMR12
      ENDIF 
C
C     CCRHSN assumes that T2AM can hold the cluster amplitudes 
C     or the vector function in different storage schemes 
C     (triangular, squared, half transformed)
C     --> R12 doubles cannot be stored directly after conv. doubles
      KFOCKD  = 1
      KT1AM   = KFOCKD  + NORBTS
      KOMEG1  = KT1AM   + NT1AMX
      KOMEG2  = KOMEG1  + NT1AM(ISYMOP)
      KTAMP12 = KOMEG2  + NT2AMX
      KT2AM   = KOMEG2  +
     *          MAX(NTAMP,NT2AO(ISYMOP),2*NT2ORT(ISYMOP))
      IF ( (KTAMP12 + NTAMR12) .GT. KT2AM )
     *     CALL QUIT('Allocation error for KTAMP12 in CCSD_ENERGY!')
      KEND1   = KT2AM   + 
     *          MAX(NT2SQ(ISYMOP),(NT2AMX+NTAMR12),NT2R12(1),NTG2SQ(1)) 
                ! CCRHSN uses T2AM for a squared array.
                ! This implies also, that we cannot store 
                ! the R12 doubles right after the doubles
                ! before calling ccrhs.
      IF (CCPAIR) THEN
C        Work space for printing of pair energies (WK/UniKA/21-11-2002).
         KES     = KEND1
         KET     = KES     + NRHFT * (NRHFT + 1)/2
         KQS     = KET     + NRHFT * (NRHFT + 1)/2
         KQT     = KQS     + NRHFT * (NRHFT + 1)/2
         KT1S    = KQT     + NRHFT * (NRHFT + 1)/2
         KT1T    = KT1S    + NRHFT * (NRHFT + 1)/2
         KT2S    = KT1T    + NRHFT * (NRHFT + 1)/2
         KT2T    = KT2S    + NRHFT * (NRHFT + 1)/2
         KEND1   = KT2T    + NRHFT * (NRHFT + 1)/2
      END IF 
      
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
      IF (INTTR) THEN
         CALL DZERO(WORK(KT1AM),NT1AMX)
         LHTF = .FALSE.
         CALL CCSD_IAJB(WORK(KT2AM),WORK(KT1AM),LHTF,
     *                       .FALSE.,.FALSE.,WORK(KEND1),LWRK1)
         REWIND(LUIAJB)
         CALL WRITI(LUIAJB,IRAT*NT2AM(ISYMOP),WORK(KT2AM))
      ELSE
         CALL DCOPY(NT2AM(ISYMOP),99.99D0,0,WORK(KT2AM),1)
         REWIND(LUIAJB)
         CALL READI(LUIAJB,IRAT*NT2AM(ISYMOP),WORK(KT2AM))
      ENDIF
C
C----------------------------------------------------------------------
C     Setup the initial guess vector:
C       1) if CCRSTR flag set try to restart from old amplitude vector
C          (ignored for MP2 calculations)
C       2) if CCRSTR flag not set or if restart failed or if we do
C          a MP2 calculation, construct MP2 amplitude vector from
C          the integrals, which we have in memory
C----------------------------------------------------------------------
C
      IF (CCRSTR.AND.(.NOT.MP2)) THEN
         ETY1   = 'RSTAR'
         IOPT   = 33
         CALL CC_RDRSP('R0',0,1,IOPT,MODELR,WORK(KT1AM),WORK(KT2AM))
         IF (IOPT.EQ.33) THEN
           INQUIRE(FILE='CCSD_TAM',EXIST=LEXIST,IOSTAT=IOS,ERR=990)
           IF (LEXIST) THEN ! read old CCSD_TAM file
             LUTAM = -1
             CALL GPOPEN(LUTAM,'CCSD_TAM','UNKNOWN',' ','UNFORMATTED',
     *                   IDUMMY,.FALSE.)
             REWIND (LUTAM)
             WRITE(LUTAM) (WORK(KT1AM+I-1), I = 1,NT1AMX)
             IF (.NOT.CCS) WRITE(LUTAM) (WORK(KT2AM+I-1), I = 1,NT2AMX)
             CALL GPCLOSE(LUTAM,'KEEP')
             IOPT = 3
           END IF
990        CONTINUE ! nothing to restart from ...
         END IF
      ENDIF
   
      IF  ( (.NOT.CCRSTR) .OR. MP2 .OR. (IOPT.EQ.33) ) THEN
         IF (CCPAIR) THEN
C           Print MP2 pair energies (WK/UniKA/21-11-2002).
            CALL CCSD_CBS1(WORK(KT2AM),WORK(KFOCKD),
     *                     WORK(KES),WORK(KET),
     *                     WORK(KQS),WORK(KQT))
         END IF
         IF (CCR12.AND.CC2.AND.(IANR12.EQ.2.OR.IANR12.EQ.3)) THEN
           LRES = .FALSE.
           CALL CCRHS_EPPP(WORK(KT2AM),WORK(KEND1),LWRK1,APROXR12,
     &                     LRES,IDUMMY,CDUMMY,IDUMMY,1)
         END IF
         CALL CCSD_GUESS(WORK(KT1AM),WORK(KT2AM),WORK(KFOCKD),IPRINT)
         IF (IPRINT .GT. 4) THEN
            CALL AROUND('Largest amplitudes in MP2 guess')
            CALL DCOPY(NT1AMX,WORK(KT1AM),1,WORK(KOMEG1),1)
            CALL DCOPY(NT2AMX,WORK(KT2AM),1,WORK(KOMEG1+NT1AMX),1)
            CALL CC_PRAM(WORK(KOMEG1),PT1,1,.FALSE.)
         ENDIF
         ETY1 = 'MP2  '
 
      ENDIF
C
C-----------------------------------------------------------------------
C     START OF ITERATIVE LOOP
C-----------------------------------------------------------------------
C
      EN1=0D0
      EN2=99D0
      EN1R12 = 0.0d0
      EN2R12 = 0.0d0
      LCONVG=.FALSE.

!radovan: otherwise er12 and en1r12 and en2r12 become undefined if lr12 = .false.
      er12 = 0.0d0

      ITER=1
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
      IF (CCR12.AND..NOT.(CCS.OR.CIS)) THEN
C       read R12 amplitudes
        IF (R12PRP) THEN
          CALL DZERO(WORK(KTAMP12),NTR12AM(1))
        ELSE
          IOPT = 32
          CALL CC_RDRSP('R0 ',0,1,IOPT,MODELR,DUMMY,WORK(KTAMP12))
        END IF
      END IF
C
      IT1 = 0
      IF ( ETY1.EQ.'RSTAR' ) IT1 = 1
      CALL CCSD_ECCSD(WORK(KT1AM),WORK(KT2AM),WORK(KFOCKD),
     *                WORK(KTAMP12),
     *                WORK(KEND1),LWRK1,EN2,POTNUC,ESCF,
     *                ETY1,ER12,LMULBS,IT1,ITER,APROXR12)
C
      EINI = EN2
C
CSPAS: 15.11.2009 adding AO-SOPPA
C-----------------------------------------------
C     For AO-SOPPA Write MP2 amplitudes to disk.
C-----------------------------------------------
C
      IF (AOSOPPA) THEN
         LUTAM = -1
         CALL GPOPEN(LUTAM,'MP2__TAM',' ',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND LUTAM
         WRITE(LUTAM) (WORK(KT1AM+I-1), I = 1,NT1AMX)
         WRITE(LUTAM) (WORK(KT2AM+I-1), I = 1,NT2AMX)
         CALL GPCLOSE(LUTAM,'KEEP')
      ENDIF
CKeinSPASmehr
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
         IOPT = 3
         CALL CC_WRRSP('R0',0,1,IOPT,MODEL,DUMMY,WORK(KT1AM),
     *                 WORK(KT2AM),WORK(KEND1),LWRK1)
         IF (CCR12.AND..NOT.(CCS.OR.CIS)) THEN
           IOPT =32
           CALL CC_WRRSP('R0',0,1,IOPT,MODEL,DUMMY,DUMMY,
     *                 WORK(KTAMP12),WORK(KEND1),LWRK1)
         END IF
C
         EN1 = EN2
         EN1R12 = ER12
C
         IF ( IPRINT .GT. 2 ) THEN
            WRITE(LUPRI,249) ITER
  249       FORMAT(/,3X,'    Iteration no.:',I3)
            WRITE(LUPRI,*)'      -----------------'
            WRITE(LUPRI,*)
         ENDIF

         CALL CCRHSN(WORK(KOMEG1),WORK(KOMEG2),WORK(KT1AM),WORK(KT2AM),
     *               WORK(KEND1),LWRK1,APROXR12)

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
            WRITE(LUPRI,529) 'Norm^2 of t1am   after ccvec:',
     *               DDOT(NT1AMX,WORK(KT1AM),1,WORK(KT1AM),1)
            WRITE(LUPRI,529) 'Norm^2 of t2am   after ccvec:',
     *               DDOT(NT2AMX,WORK(KT2AM),1,WORK(KT2AM),1)
         ENDIF
         OMNM1 = DDOT(NT1AM(ISYMOP),WORK(KOMEG1),1,WORK(KOMEG1),1)
         OMNM2 = DDOT(NT2AM(ISYMOP),WORK(KOMEG2),1,WORK(KOMEG2),1)
         OMNM  = DSQRT(OMNM1+OMNM2)
         IF (IPRINT .GE. 3) THEN
            WRITE(LUPRI,529) 'Norm^2 of omega1 after ccvec:',OMNM1
            WRITE(LUPRI,529) 'Norm^2 of omega2 after ccvec:',OMNM2
         END IF

C
         IF (CCSLV) THEN
           IF (IPRINT .GE. 3) THEN
              WRITE(LUPRI,529) 'Norm^2 of omega1 in sol. part.:',
     *              DDOT(NT1AM(ISYMOP),WORK(KOMEG1),1,WORK(KOMEG1),1)
              WRITE(LUPRI,529) 'Norm^2 of omega2 in sol. part.:',
     *              DDOT(NT2AM(ISYMOP),WORK(KOMEG2),1,WORK(KOMEG2),1)
           END IF
           LUSLV = -1
           CALL GPOPEN(LUSLV,'CC_OME','UNKNOWN',' ','UNFORMATTED',
     *                 IDUMMY,.FALSE.)
           REWIND (LUSLV)
           WRITE(LUSLV) (WORK(KOMEG1+I-1), I = 1,NT1AMX)
           WRITE(LUSLV) (WORK(KOMEG2+I-1), I = 1,NT2AMX)
           CALL GPCLOSE(LUSLV,'KEEP')
         ENDIF
C
  529    FORMAT(7X,A,D24.10)
C
         CALL CCSD_NXTAM(WORK(KT1AM),WORK(KT2AM),DUMMY,WORK(KOMEG1),
     *                   WORK(KOMEG2),DUMMY,WORK(KFOCKD),.FALSE.,
     *                   ISYMOP,0.0D0)
C        IF (IPRINT .GE. 5) THEN
C           WRITE(LUPRI,529) 'Norm^2 of t1am   after NXTAM:',
C    *               DDOT(NT1AMX,WORK(KT1AM),1,WORK(KT1AM),1)
C           WRITE(LUPRI,529) 'Norm^2 of t2am   after NXTAM:',
C    *               DDOT(NT2AMX,WORK(KT2AM),1,WORK(KT2AM),1)
C        ENDIF
C        IF (IPRINT .GE. 3) THEN
C           WRITE(LUPRI,529) 'Norm^2 of omega1 after NXTAM:',
C    *             DDOT(NT1AM(ISYMOP),WORK(KOMEG1),1,WORK(KOMEG1),1)
C           WRITE(LUPRI,529) 'Norm^2 of omega2 after NXTAM:',
C    *             DDOT(NT2AM(ISYMOP),WORK(KOMEG2),1,WORK(KOMEG2),1)
C        END IF
C         ------------------------------------------------
C         set address for singles amplitudes in DIIS 
C         (have to be just before the doubles amplitudes)
C         ------------------------------------------------
          KTAMP1  = KT2AM   - NT1AMX
         IF ( KTAMP1 .LT. (KOMEG2+NT2AMX+NTAMR12) )
     *     CALL QUIT('Allocation error in CCSD_ENERGY!')

         IF (LMULBS) THEN
C          -----------------------------------------------------------
C          for DIIS we put the R12 part of the vector function behind
C          the conventional doubles part 
C          -----------------------------------------------------------
           KOMEG12 = KOMEG2   + NT2AMX
           IF ( (KOMEG12+NTR12AM(1)) .GT. KTAMP1 )
     *       CALL QUIT('Allocation error in CCSD_ENERGY!')
C          ---------------------------------------
C          and similar for the cluster amplitudes
C          ---------------------------------------
           KTAMP12 = KT2AM   + NT2AMX
C          --------------------------------------------------------
C          read r12 part of the vector function and the amplitudes
C          and apply the perturbational update for the amplitudes:
C          --------------------------------------------------------
             LCCEQ = .TRUE.
             CALL CC_R12NXTAM(WORK(KOMEG12),1,
     &                      WORK(KTAMP12),LCCEQ,
     &                      ER12,WORK(KEND1),LWRK1)
C            IF (IPRINT .GE. 5) THEN
C              WRITE(LUPRI,529) 'Norm^2 of R12amp. after NXTAM:',
C    *               DDOT(NTR12AM(1),WORK(KTAMP12),1,WORK(KTAMP12),1)
C            ENDIF
         END IF 

         IF (CCSTST) THEN
            CALL DZERO(WORK(KT2AM),NT2AMX)
         ENDIF
         IF (CCD) THEN
            CALL CCSD_DIIS(WORK(KOMEG2),WORK(KT2AM),NTAMP,ITER)
         ELSE
            CALL DCOPY(NT1AMX,WORK(KT1AM),1,WORK(KTAMP1),1)
            CALL CCSD_DIIS(WORK(KOMEG1),WORK(KTAMP1),NTAMP,ITER)
            CALL DCOPY(NT1AMX,WORK(KTAMP1),1,WORK(KT1AM),1)
         ENDIF
C
C---------------------------------
C        The order is important !!
C---------------------------------
C
         ETY2 = MODEL(1:5)//'     '
         IF (CCSD.OR.CCPT.OR.CCP3)  ETY2 = 'CCSD '
         IF (CCD)   ETY2 = 'CCD  '
         IF (CCSDT) ETY2 = 'CC3  '
         IF (CC1B)  ETY2 = 'CC-1b'
         IF (CC1A)  ETY2 = 'CC-1a'
         IF (CC2)   ETY2 = 'CC2  '
C
         IT1 = 1
         CALL CCSD_ECCSD(WORK(KT1AM),WORK(KT2AM),WORK(KFOCKD),
     *                   WORK(KTAMP12),
     *                   WORK(KEND1),LWRK1,EN2,POTNUC,ESCF,
     *                   ETY2,ER12,LMULBS,IT1,ITER,
     *                   APROXR12)
C
         EN2R12 = ER12
C
         IF (CCR12.AND..NOT.(CCS.OR.CIS)) THEN
C          Save new R12 amplitudes on disk:
           LUNIT = -1
           CALL GPOPEN(LUNIT,FCCR12C,'UNKNOWN',' ','UNFORMATTED',
     &                       IDUM,LDUM)
           WRITE(LUNIT) (WORK(KTAMP12-1+I),I=1,NTR12AM(1))
           CALL GPCLOSE(LUNIT,'KEEP')
         END IF
C
         CALL FLSHFO(LUPRI)
C
         LCONV1 = DABS(EN2-EN1)      .LT. THRENR .AND.
     &            DABS(EN2R12-EN1R12).LT. THRENR
         LCONV2 = OMNM    .LT. THRVEC
         LCONVG = LCONV1 .AND. LCONV2

         ITER = ITER+1
         IF (ITER .GT. MAXITE) THEN
           WRITE(LUPRI,*) 'Energy not converged in ',MAXITE,
     &           ' iterations'
           CALL QUIT('CC equations not converged.')
         ENDIF
C
         NSLVINIT = NSLVINIT + 1

         IF (CCSLV.AND.((ITER-1).GE.MXTINIT)) THEN
           WRITE(LUPRI,241) ETY2
 241       FORMAT(/,1x,A5,
     *    'energy will not be converged further'
     *     //' right now in CCSLV calc.')
           WRITE(LUPRI,242) NSLVINIT
 242       FORMAT(' Accumulated inner iterations at this point are ',I5)
C
         ELSE
           IF (.NOT. LCONVG) THEN
              IF (IPRINT.GE.3 .AND. LCONV1) THEN
                 WRITE(LUPRI,'(3X,A,D15.6,A,D15.6)')
     &           'Energy difference ',DABS(EN2-EN1),
     &           'is less then threshold ',THRENR
                 WRITE(LUPRI,'(3X,A,A,D15.6,/,3X,A,D15.6)')
     &           'Iterations continue, as the 2-norm of the vector ',
     &           ' function: ',OMNM,
     &           'is larger than the threshold: ',THRVEC
                 CALL FLSHFO(LUPRI)
              END IF
              GOTO 200 ! go to next iteration
           END IF
         ENDIF
C
      CALL CCSD_MODEL(ETYPE,LENET,24,ETY2,5,APROXR12)

      WRITE(LUPRI,250) ETYPE(1:LENET),THRENR,EN2
      WRITE(LUPRI,'(1X,A,1P,D15.8)')
     & 'Final 2-norm of the CC vector function: ',OMNM
 250  FORMAT(/,1x,A,' energy converged to within ',D10.2,' is ',F16.9)
      IF (CCR12) THEN
        WRITE(LUPRI,'(A,D10.2,a,F16.9)') 
     &  ' R12       energy converged to within ',THRENR,' is ',ER12
      END IF
c
      IF (CCPAIR) THEN 
C        Print CC pair energies (WK/UniKA/21-11-2002).
         CALL CCSD_CBS2(WORK(KT1AM),WORK(KT2AM),WORK(KEND1),LWRK1,
     *                  WORK(KT1S),WORK(KT1T),WORK(KT2S),WORK(KT2T),
     *                  ETY2)
      END IF
C
C-----------------
C     end of loop.
C-----------------
C

#if defined (SYS_CRAY)
      CALL WCLOSE('CC_DIIS',IERR)
      INFO = ISHELL('rm CC_DIIS')
#else
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
         CALL CC_PRAM(WORK(KOMEG1),PT1,1,.TRUE.)
C
      ENDIF
C
 500  CONTINUE
C
CSPAS 22.10.2003 implementing SOPPA(CCSD)
C
C------------------------------------------------------
C     Write interface to SIRIUS SOPPA response program.
C------------------------------------------------------
C
      IF (SIRSOP) THEN
C
         CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &              .FALSE.)
         REWIND LUSIFC
C
C
         IERR = -1
         CALL MOLLAB('MP2INFO ',LUSIFC,IERR)
         IF (IERR.EQ.0) THEN
            READ(LUSIFC)
            READ(LUSIFC)
         ELSE IF (IERR.EQ.-1) THEN
            REWIND LUSIFC
            CALL MOLLAB('EODATA  ',LUSIFC,LUPRI)
            BACKSPACE(LUSIFC)
         ENDIF
         CALL GETDAT(CDATE,CTIME)
         IF (MP2) THEN
            WRITE(LUSIFC) '********',CDATE,CTIME,'MP2INFO '
         ELSE
            WRITE(LUSIFC) '********',CDATE,CTIME,'CCSDINFO'
         ENDIF
C
         KSCR1 = KEND1
         KSCR2 = KSCR1 + NT2SQ(1)
         KSCR3 = KSCR2 + NT2SQ(1)
         KEND2 = KSCR3 + NORBTS*NORBTS
         LWRK2 = LWORK - KEND2
C
        IF (LWRK2 .LT. 0) THEN
            STOP 'Insufficient memory in CCSD_ENERGY'
         ENDIF
C
         CALL CC_T2SQ(WORK(KT2AM),WORK(KSCR1),1)
         CALL CCRHS_T2TR(WORK(KSCR1),WORK(KEND2),LWRK2,1)
         CALL T2AM_REORDER(WORK(KSCR1),WORK(KSCR2),IPRINT)
C
         WRITE(LUSIFC) (WORK(KSCR2+I-1), I = 1,NT2SQ(1))
C
C----------------------------------------------------------
C        Calculate density matrices D(ij), D(ab) and D(ia).
C----------------------------------------------------------
C
         CALL CC_T2SQ(WORK(KT2AM),WORK(KSCR2),1)
C
         CALL SOPPA_DENSITY(WORK(KSCR3),WORK(KT1AM),WORK(KSCR2),
     *                       WORK(KSCR1),IPRINT)
C
         WRITE(LUSIFC) (WORK(KSCR3+I-1), I = 1,NORBTS*NORBTS)
C
         WRITE(LUSIFC) '********',CDATE,CTIME,'EODATA  '
C
         CALL GPCLOSE(LUSIFC,'KEEP')
C
C
      END IF
CKeinSPASmehr
C
C------------------------------------
C     Write final amplitudes to disk.
C------------------------------------
C
      IF (CCSTST) THEN
         CALL DZERO(WORK(KT2AM),NT2AMX)
      ENDIF
C
CSPAS: 15.11.09 adding AO-SOPPA
      IF (AOSOPPA .AND. .NOT. MP2) THEN
         LUTAM = -1
         CALL GPOPEN(LUTAM,'CCSD_TAM',' ',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND LUTAM
         WRITE(LUTAM) (WORK(KT1AM+I-1), I = 1,NT1AMX)
         WRITE(LUTAM) (WORK(KT2AM+I-1), I = 1,NT2AMX)
         CALL GPCLOSE(LUTAM,'KEEP')
      ENDIF
CKeinSPASmehr
C
C----------------------------------
C     save a copy on file CCR0___0
C----------------------------------
C
      KT0AM   = KEND1
      KEND2   = KT0AM   + 2*NALLAI(1)
      LWRK2   = LWORK - KEND2

      IF ( LWRK2 .LT. 0 ) THEN
         write(lupri,*) 'LWORK, LWRK2: ',LWORK, LWRK2
         CALL QUIT('Insufficient spaces in CCSD_ENERGY (2)')
      ENDIF

      CALL DZERO(WORK(KT0AM),2*NALLAI(1))
C
      IOPT = 7
      CALL CC_WRRSP('R0',0,1,IOPT,MODEL,WORK(KT0AM),WORK(KT1AM),
     *              WORK(KT2AM),WORK(KEND2),LWRK2)
C
C     --------------------------------------------------------
C     for CC-R12 save also the R12 amplitudes on CCR0... file
C     --------------------------------------------------------
C
      !R12 amps. are still on KTAMP12 = KT2AM + NT2AMX!
      !do not overwrite them before!
      IF (CCR12) THEN
        IOPT=32
        CALL CC_WRRSP('R0 ',0,1,IOPT,MODEL,WORK(KT0AM),DUMMY,
     &                WORK(KTAMP12),WORK(KEND2),LWRK2)
C
      END IF
 
C
C SLV98,OC
C
      IF (CCSLV) THEN
        XTNCCCU = DDOT(NT2AMX,WORK(KT2AM),1,WORK(KT2AM),1)
     *          + DDOT(NT1AMX,WORK(KT1AM),1,WORK(KT1AM),1)
        IF (ABS(XTNCCPR-XTNCCCU).LT.CVGTSOL) LSLTCVG = .TRUE.
        IF (IPRINT.GT.2) THEN
          WRITE(LUPRI,*)'Norm^2 of T-amplitudes in this solvent it.:',
     &                  XTNCCCU
          WRITE(LUPRI,*)'Norm^2 of T-amplitudes in prev solvent it.:',
     &                  XTNCCPR
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
         IF (.NOT. CHOPT) THEN
C
            IF (CCPT) THEN
               CC1BSV = CC1B
               CC1B   = .TRUE.
               CC1ASV = CC1A
               CC1A   = .TRUE.
            ENDIF
C
            CALL CCRHSN(WORK(KOMEG1),WORK(KOMEG2),WORK(KT1AM),
     *               WORK(KT2AM),WORK(KEND1),LWRK1,APROXR12)
C
            IF (CCPT) THEN
               CC1B   = CC1BSV
               CC1A   = CC1ASV
            ENDIF
C
            IOPTTCME = 1
            CALL CCSD_TCMEPK(WORK(KT2AM),0.5d0,1,IOPTTCME)
            ECCP1 = TWO*DDOT(NT1AMX,WORK(KT1AM),1,WORK(KOMEG1),1)
            ECCP2 = TWO*DDOT(NT2AMX,WORK(KT2AM),1,WORK(KOMEG2),1)
C
            IOPT = 3
            CALL CC_RDRSP('R0',0,1,IOPT,MODELR,WORK(KT1AM),WORK(KT2AM))
C
         ELSE
C
C----------------------------
C           Cholesky CCSD(T)
C----------------------------
C
            INQUIRE(FILE='CC_CHOPT_DBG',EXIST=CPTDBG)
            IF (CPTDBG) THEN
               WRITE(LUPRI,'(//,15X,A,/,15X,A)')
     &               '*** NOTICE ***',
     &               'File CC_CHOPT_DBG found. Calling CC_CHOPT_DBG.'
               CALL CC_CHOPT_DBG(WORK,LWORK)
               GOTO 9999
            ENDIF
C
C-----------------------------------------------------------------------
C           NB! This assumes that the orbital energies and
C           T1 amplitudes are in fact in these positions:
C-----------------------------------------------------------------------
C
            KFOCKD  = 1
            KT1AM   = KFOCKD + NORBTS
            KEND1   = KT1AM  + NT1AMX
            LWRK1   = LWORK  - KEND1
C
            CHOTIM = SECOND()
            CALL CC_CHOPT(WORK(KFOCKD),WORK(KT1AM),WORK(KEND1),LWRK1)
            CHOTIM = SECOND() - CHOTIM
            WRITE(LUPRI,9998) 'CC_CHOPT',CHOTIM
 9998       FORMAT(7X,'Time used in',2X,A12,2X,': ',F10.2,' seconds')
C
            ECCP1 = XEN5
            ECCP2 = XEN4
C
         END IF
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
               WRITE(LUPRI,*) 'Norm^2 of T1AM: ',RHO1N
               WRITE(LUPRI,*) 'Norm^2 of T2AM: ',RHO2N
            ENDIF
C
            IF (IPRINT .GT. 45) THEN
               CALL CC_PRP(WORK(KOMEG1),WORK(KOMEG2),1,1,1)
            ENDIF
C
            IOPT = 3
            CALL CC_WRRSP('R0',0,1,IOPT,MODEL,DUMMY,WORK(KT1AM),
     *                    WORK(KT2AM),WORK(KEND1),LWRK1)
C
            IF (IPRINT .GT. 4) THEN
C
               CALL AROUND('Largest amplitudes in pert. corr. ampl.')
C
               CALL DCOPY(NT1AMX,WORK(KT1AM),1,WORK(KOMEG1),1)
               CALL DCOPY(NT2AMX,WORK(KT2AM),1,WORK(KOMEG1+NT1AMX),1)
C
               CALL CC_PRAM(WORK(KOMEG1),PT1,1,.FALSE.)
C
            ENDIF
C
         ENDIF ! end CCR3
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
      IF (.NOT. QM3) THEN
      WRITE(LUPRI,'(//)')
      CALL AROUND
     * ('Final results from the Coupled Cluster energy program')
      WRITE(LUPRI,'(//12X,A,A,A,F32.10,/)')
     &   'Total ',ETY0,' energy: ',ESCF
      WRITE(LURES,'(//12X,A,A,A,F32.10)') 'Total ',ETY0,' energy: ',ESCF
      CALL CCSD_MODEL(ETYPE,LENET,24,ETY1,5,APROXR12)
      WRITE(LUPRI,'(12X,A,A,A,A,F25.10,/)') 'Total ',ETYPE(1:LENET),
     &    ' energy: ',BLANKS(1:12-LENET),EINI
      WRITE(LURES,'(12X,A,A,A,A,F25.10)')   'Total ',ETYPE(1:LENET),
     &    ' energy: ',BLANKS(1:12-LENET),EINI
      IF (.NOT. (MP2 .OR. NOCCIT)) THEN
        CALL CCSD_MODEL(ETYPE,LENET,24,ETY2,5,APROXR12)
        WRITE(LUPRI,'(12X,A,A,A,A,F25.10)') 'Total ',ETYPE(1:LENET),
     &     ' energy: ',BLANKS(1:12-LENET),EN2
        WRITE(LURES,'(12X,A,A,A,A,F25.10)') 'Total ',ETYPE(1:LENET),
     &     ' energy: ',BLANKS(1:12-LENET),EN2
      END IF
      END IF
      ECCGRS = EN2
      IF ( (.NOT.CCSLV) .OR. (ICCSLIT.EQ.0)) THEN
         LABEL1 = 'ENERGY  '
         MODREF = 'SCF       '
         CALL CC_PRPC(ESCF,MODREF,0,
     *                LABEL1,LABEL1,LABEL1,LABEL1,
     *                ZERO,ZERO,ZERO,1,0,0,0)
      ENDIF
      IF (.NOT.CCSLV) THEN
         LABEL1 = 'ENERGY  '
         CALL CC_PRPC(EN2,MOPRPC,0,
     *                LABEL1,LABEL1,LABEL1,LABEL1,
     *                ZERO,ZERO,ZERO,1,0,0,0)
         IF (CCPTSV .OR. CCP3SV) THEN
            CALL WRIPRO(ETOT,MOPRPC,0,
     *                  LABEL1,LABEL1,LABEL1,LABEL1,
     *                  ZERO,ZERO,ZERO,1,0,0,0)
         ELSE
            CALL WRIPRO(EN2,MOPRPC,0,
     *                  LABEL1,LABEL1,LABEL1,LABEL1,
     *                  ZERO,ZERO,ZERO,1,0,0,0)
         ENDIF
      ENDIF
C
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
         IF (CCPTSV.AND..NOT.CCR12) THEN
            WRITE(LUPRI,'(/,12X,A,F31.10)') 'Total energy CCSD(T):',ETOT
            WRITE(LURES,'(/,12X,A,F31.10)') 'Total energy CCSD(T):',ETOT
         ELSE IF ((CCPTSV.AND.CCR12).AND.(IANR12.EQ.1)) THEN
            WRITE(LUPRI,'(/,12X,A,A,A,F23.10)') 
     &            'Total CCSD(R12)(T)/',APROXR12,'energy:',ETOT
            WRITE(LURES,'(/,12X,A,A,A,F23.10)') 
     &            'Total CCSD(R12)(T)/',APROXR12,'energy:',ETOT
         ELSE IF ((CCP3SV.AND.CCR12).AND.(IANR12.EQ.1)) THEN
            WRITE(LUPRI,'(/,12X,3A,F25.10)') 'Total CC(3)(R12)/',
     &            APROXR12,'energy:',ETOT
            WRITE(LURES,'(/,12X,3A,F25.10)') 'Total CC(3)(R12)/',
     &            APROXR12,'energy:',ETOT
         ELSE IF (CCP3SV.AND..NOT.CCR12) THEN
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
     *               WORK(KEND1),LWRK1,APROXR12)
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

C------------------------------------------------------------
C     Precompute intermediates needed for CC-R12 left transf.
C------------------------------------------------------------
      IF (CCR12LIM .AND. .NOT.(CCS.OR.CC2)) THEN
C
         KT1AM  = 1
         KLAMDP = KT1AM + NT1AMX
         KLAMDH = KLAMDP + NLAMDT
         KVABKL = KLAMDH + NLAMDT
         KEND1  = KVABKL + NVABKL(1)
         LWRK1  = LWORK - KEND1 
         IF (LWRK1.LT.0) 
     &     CALL QUIT('Insufficient memory for VABKL in CCSD_ENERGY')
C
         IOPT = 1
         CALL CC_RDRSP('R0 ',0,1,IOPT,MODELR,WORK(KT1AM),DUMMY)
         CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),
     &               WORK(KEND1),LWRK1)
C  
         LUNIT = -1
         CALL GPOPEN(LUNIT,FVABKL,'OLD',' ','UNFORMATTED',IDUM,.FALSE.)
         READ(LUNIT)(WORK(KVABKL+I-1),I=1,NVABKL(1))
         CALL GPCLOSE(LUNIT,'KEEP')
C
         ! calculate V_(\tilde{a} \tilde{b})^(kl) and save on disk:
         IOPT = 1
         CALL CC_R12MKVIRT(WORK(KVABKL),WORK(KLAMDP),1,WORK(KLAMDP),1,
     &                     'R12VCTDTKL',IOPT,WORK(KEND1),LWRK1)
      END IF
C--------------------------------------------------
C     Precompute some integrals and amplitudes for
C     the CC3 noddy response code:
C--------------------------------------------------
      IF (NODDY_INIT) THEN
        CALL CCSDT_INIT_NODDY(WORK,LWORK,.FALSE.) 
      END IF

 9999 CALL QEXIT('CCSD_ENERGY')
C
      RETURN
 1817 CALL QUIT('R12 amplitudes not found on disk')
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
      CALL QENTER('CCSD_GUESS')
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
      CALL QEXIT('CCSD_GUESS')
C
      RETURN
      END
C  /* Deck ccsd_nxtam */
      SUBROUTINE CCSD_NXTAM(T1AM,T2AM,T2AM2,OMEGA1,OMEGA2,OMEGA22,
     *                      FCDIAG,TRIPLET,ISYMT,FREQ)
C
C     Written by Henrik Koch 27-Mar-1990.
C     Brueckner bit by Rika Kobayashi 1992.
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, TWO = 2.0D0, THREE = 3.0D0)
      DIMENSION T1AM(*),T2AM(*), T2AM2(*)
      DIMENSION OMEGA1(*),OMEGA2(*), OMEGA22(*)
      DIMENSION FCDIAG(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
Cholesky
#include <ccsdinp.h>
#include <maxorb.h>
#include <ccdeco.h>
Cholesky
C
      LOGICAL TRIPLET
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CCSD_NXTAM')
C
c     IF (.NOT. (CCD.OR.LBRUK)) THEN
C
         DO 100 ISYMI = 1,NSYM
            ISYMA = MULD2H(ISYMT,ISYMI)
            DO 110 I = 1,NRHF(ISYMI)
               KOFFI = IRHF(ISYMI) + I
               DO 120 A = 1,NVIR(ISYMA)
C
                  KOFFA = IVIR(ISYMA) + A
                  NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A
C
                  OMEGA1(NAI) = T1AM(NAI) + OMEGA1(NAI)/
     *                       (FREQ + FCDIAG(KOFFI) - FCDIAG(KOFFA))
C
  120          CONTINUE
  110       CONTINUE
  100    CONTINUE
C
c     ENDIF
c     IF (LBRUK) CALL DCOPY(NVIRT*NRHFT,OMEGA1,1,T1AM,1)
C
C
      IF (CC2 .AND. CHOINT) GOTO 1000  ! Skip doubles part for Cholesky CC2.
C
      DO 200 ISYMBJ = 1,NSYM
         ISYMAI = MULD2H(ISYMBJ,ISYMT)
         IF (ISYMAI .LE. ISYMBJ) THEN
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
                           IF (ISYMAI.EQ.ISYMBJ .AND. NAI.LE.NBJ) THEN
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                            + INDEX(NAI,NBJ)
                              OMEGA2(NAIBJ) = T2AM(NAIBJ)+OMEGA2(NAIBJ)/
     *                                   (FREQ + 
     *                                    FCDIAG(KOFFI) + FCDIAG(KOFFJ)
     *                                  - FCDIAG(KOFFA) - FCDIAG(KOFFB))
                           ELSE IF (ISYMAI.LT.ISYMBJ) THEN
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                            + NT1AM(ISYMAI)*(NBJ-1) + NAI
                              OMEGA2(NAIBJ) = T2AM(NAIBJ)+OMEGA2(NAIBJ)/
     *                                   (FREQ + 
     *                                    FCDIAG(KOFFI) + FCDIAG(KOFFJ)
     *                                  - FCDIAG(KOFFA) - FCDIAG(KOFFB))
                           ENDIF
C
C
  260                   CONTINUE
  250                CONTINUE
  240             CONTINUE
  230          CONTINUE
  220       CONTINUE
  210    CONTINUE
         ENDIF
  200 CONTINUE
C
C     Do second double block if triplet.
C
C
      IF (TRIPLET) THEN
C
      DO ISYMBJ = 1,NSYM
         ISYMAI = MULD2H(ISYMBJ,ISYMT)
         IF (ISYMAI .LE. ISYMBJ) THEN
         DO ISYMJ = 1,NSYM
            ISYMB = MULD2H(ISYMJ,ISYMBJ)
            DO ISYMI = 1,NSYM
               ISYMA = MULD2H(ISYMI,ISYMAI)
               DO J = 1,NRHF(ISYMJ)
                  KOFFJ = IRHF(ISYMJ) + J
                  DO B = 1,NVIR(ISYMB)
                     KOFFB = IVIR(ISYMB) + B
                     NBJ = IT1AM(ISYMB,ISYMJ)+NVIR(ISYMB)*(J-1)+B
                     DO I = 1,NRHF(ISYMI)
                        KOFFI = IRHF(ISYMI) + I
                        DO A = 1,NVIR(ISYMA)
                           KOFFA = IVIR(ISYMA) + A
                           NAI = IT1AM(ISYMA,ISYMI)+NVIR(ISYMA)*(I-1)+A
C
                           IF (ISYMAI.EQ.ISYMBJ .AND. NAI.LE.NBJ) THEN
                             NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                           + INDEX(NAI,NBJ)
                             OMEGA22(NAIBJ)=T2AM2(NAIBJ)+OMEGA22(NAIBJ)/
     *                                  (FREQ + 
     *                                   FCDIAG(KOFFI) + FCDIAG(KOFFJ)
     *                                 - FCDIAG(KOFFA) - FCDIAG(KOFFB))
                           ELSE IF (ISYMAI.LT.ISYMBJ) THEN
                             NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                           + NT1AM(ISYMAI)*(NBJ-1) + NAI
                             OMEGA22(NAIBJ)=T2AM2(NAIBJ)+OMEGA22(NAIBJ)/
     *                                  (FREQ + 
     *                                   FCDIAG(KOFFI) + FCDIAG(KOFFJ)
     *                                 - FCDIAG(KOFFA) - FCDIAG(KOFFB))
                           ENDIF
C
C
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
         ENDIF
      END DO
      ENDIF
C
 1000 CONTINUE
      CALL QEXIT('CCSD_NXTAM')
C
      RETURN
      END
C  /* Deck ccsd_eccsd */
      SUBROUTINE CCSD_ECCSD(T1AM,T2AM,FCDIAG,TAMR12,
     *                      WORK,LWORK,XECCSD,POTNUC,
     *                      ESCF,ETY,ER12,LR12,IT1,ITER,
     *                      APROXR12)
C
C     Written by Henrik Koch 27-Mar-1990.
C
C     Ove Christiansen 23-1-1996: Introduction of Finite field contribution.
C                                 IT1 = 0 : no amplitudes on disk
C                                 IT1 = 1 : t1 amplitudes read from disk
C
C  Bug fix for frozen core in finite field calculations, 
C  C.Haettig, 23.3.05
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
#include <iratdef.h>
      DIMENSION FCDIAG(*)
      DIMENSION T1AM(*),T2AM(*),TAMR12(*),WORK(*)
      CHARACTER ETY*5, ETYPE*24, MODEL*10
      CHARACTER*(*) APROXR12
      LOGICAL LEXIST, LR12, LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      INTEGER ICMO(8,8), NCMO(8), IGLMRHS(8,8), IGLMVIS(8,8), NGLMDS(8)
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <ccfield.h>
#include <ccinftap.h>
#include <r12int.h>
#include <ccr12int.h>

C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CCSD_ECCSD')
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
         DO ISYM = 1, NSYM
            ICOUNT = 0
            ICOUNT3 = 0
            DO ISYM2 = 1, NSYM
               ISYM1 = MULD2H(ISYM,ISYM2)
               ICMO(ISYM1,ISYM2)    = ICOUNT
               ICOUNT  = ICOUNT  + NBAS(ISYM1)*NORBS(ISYM2)
               ICOUNT3 = ICOUNT3 + NBAS(ISYM1)*NRHFS(ISYM2)
            END DO
            NCMO(ISYM)   = ICOUNT
            NGLMDS(ISYM) = ICOUNT

            ICOUNT2 = 0
            DO ISYM2 = 1, NSYM
               ISYM1 = MULD2H(ISYM,ISYM2)
               IGLMRHS(ISYM1,ISYM2) = ICOUNT2
               IGLMVIS(ISYM1,ISYM2) = ICOUNT3
               ICOUNT2 = ICOUNT2 + NBAS(ISYM1)*NRHFS(ISYM2)
               ICOUNT3 = ICOUNT3 + NBAS(ISYM1)*NVIRS(ISYM2)
            END DO
         END DO
C
         KONEP  = 1
         KT1AM  = KONEP  + N2BST(ISYMOP)
         KLAMDPS= KT1AM  + NT1AMX
         KLAMDHS= KLAMDPS+ NGLMDS(1)
         KEND1  = KLAMDHS+ NGLMDS(1)
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
               IOPT = 1
               CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KT1AM),DUMMY)
            ELSE IF (IT1 .EQ. 0) THEN
               CALL DZERO(WORK(KT1AM),NT1AMX)
            ELSE
               CALL QUIT('IT1 should be 0 or 1 in ccsd_eccsd')
            ENDIF
         ENDIF
         CALL LAMMATS(WORK(KLAMDPS),WORK(KLAMDHS),WORK(KT1AM),
     &                1,.FALSE.,.FALSE.,
     &                NGLMDS,IGLMRHS,IGLMVIS,ICMO,WORK(KEND1),LWRK1)

         DO ISYM = 1, NSYM

           KSCR1 = KEND1
           KEND2 = KSCR1 + NBAS(ISYM) * NRHFS(ISYM)
           LWRK2 = LWORK  - KEND2
           IF ( LWRK2 .LT. 0 )
     *       CALL QUIT(' Too little workspace in ccsd_eccsd-3')

           NBAS1 = MAX(NBAS(ISYM),1)
           KOFF1 = KONEP   + IAODIS(ISYM,ISYM) 
           KOFF2 = KLAMDHS + IGLMRHS(ISYM,ISYM)

           CALL DGEMM('N','N',NBAS(ISYM),NRHFS(ISYM),NBAS(ISYM),
     *                ONE,WORK(KOFF1),NBAS1,WORK(KOFF2),NBAS1,
     *                ZERO,WORK(KSCR1),NBAS1)

           KOFF2 = KLAMDPS + IGLMRHS(ISYM,ISYM)

           TRACE = DDOT(NBAS(ISYM)*NRHFS(ISYM),
     &                    WORK(KOFF2),1,WORK(KSCR1),1)
           XECCSD = XECCSD + TWO * TRACE
         END DO
 
        ENDIF
  13  CONTINUE
C
      XCORR = XECCSD - ESCF

      ETYPE(1:5) = ETY(1:5)
      LENET = 5

      IF (LR12) THEN
C       NRHFTRIA= NRHFT * (NRHFT+1) / 2
C       N2 = NRHFTRIA * NRHFTRIA
C
C       KVR12S = 1
C       KVR12T = KVR12S + N2
C       KEND1  = KVR12T + N2
C       LWRK1  = LWORK  - KEND1
C       IF ( LWRK1 .LT. 0 )
C    *   CALL QUIT(' Too little workspace in ccsd_eccsd-3')
C
C        read V matrices
C        LUNIT = -1
C        CALL GPOPEN(LUNIT,FCCR12V,'UNKNOWN',' ','FORMATTED',
C    &                    IDUM,LDUM)
C6666     READ(LUNIT,'(I3)') IAN
C        READ(LUNIT,'(4E30.20)') (WORK(KVR12S+IJ), IJ = 0, N2-1)
C        READ(LUNIT,'(4E30.20)') (WORK(KVR12T+IJ), IJ = 0, N2-1)
C        IF (IAN.NE.IANR12) GOTO 6666
C        CALL GPCLOSE(LUNIT,'KEEP')
C
C        ER12S = DDOT(N2,WORK(KVR12S),1,TAMR12S,1)
C        ER12T = 3.0D0*DDOT(N2,WORK(KVR12T),1,TAMR12T,1)
C
C        XECCSD = XECCSD + ER12S + ER12T

        KVR12 = 1
        KEND1  = KVR12 + NTR12AM(1)
        LWRK1  = LWORK - KEND1
        IF ( LWRK1 .LT. 0 )
     *    CALL QUIT(' Too little workspace in ccsd_eccsd-3')
C
C       read V matrices
        LUNIT = -1
        CALL GPOPEN(LUNIT,FCCR12V,'UNKNOWN',' ','UNFORMATTED',
     &              IDUM,LDUM)
6666    READ(LUNIT) IAN
        READ(LUNIT) (WORK(KVR12-1+I), I=1, NTR12AM(1))
        IF (IAN.NE.IANR12) GOTO 6666
        CALL GPCLOSE(LUNIT,'KEEP')
        CALL CC_R12TCMEPK(WORK(KVR12),1,.FALSE.)
        CALL CCLR_DIASCLR12(WORK(KVR12),0.5D0,1)

        ER12 = 2.0D0*DDOT(NTR12AM(1),TAMR12,1,WORK(KVR12),1)

        XECCSD = XECCSD + ER12
       
        CALL CCSD_MODEL(ETYPE,LENET,24,ETY,5,APROXR12)
      END IF
C
      WRITE(LUPRI,'(1X,A,I3,A,A,A,F23.16)')
     *  'Iter.',ITER,': Coupled cluster ',ETYPE(1:LENET),
     *  ' energy :  ',XECCSD
C
      IF (IPRINT .GE. 2) THEN
        WRITE(LUPRI,'(5X,A,F23.16)')
     &    'Conventional correlation energy:',XCORR
        IF (LR12) THEN 
          WRITE(LUPRI,'(3(5X,A,F23.16,/))')
C    &    'Singlet R12 correlation energy :',ER12S,
C    &    'Triplet R12 correlation energy :',ER12T,
     &    'R12 correlation energy         :',ER12,
     &    'Total correlation energy       :',XCORR+ER12
        END IF
      END IF
      IF (LOCDBG) THEN
        CALL AROUND('Amplitudes at this iteration:')
        CALL CC_PRP(T1AM,T2AM,1,1,1)
        IF (CCR12) CALL CC_PRPR12(TAMR12,1,1,.TRUE.)
      END IF
C
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('CCSD_ECCSD')
C
      RETURN
      END
      SUBROUTINE CCSD_MODEL(MODELR12,LENMR12,LMAX,MODEL,LENM,APROXR12)
      IMPLICIT NONE
#include <r12int.h>
#include <ccsdinp.h>
     
      INTEGER LENM,LENMR12,LMAX,I
      CHARACTER*(*) MODELR12, MODEL, APROXR12

      IF (LMAX.LT.LENM) CALL QUIT('LMAX too small in CCSD_MODEL')
      
      IF (CCR12) THEN
        MODELR12(1:LENM) = MODEL(1:LENM)
        LENMR12 = LENM 
        DO WHILE (LENMR12.GT.0 .AND. MODELR12(LENMR12:LENMR12).EQ.' ')
          LENMR12 = LENMR12 -1
        END DO

        IF (LMAX.LT.LENMR12+5) CALL QUIT('LMAX too small in CCSD_MODEL')
        IF (MP2 .OR. CC2) THEN
          MODELR12(LENMR12+1:LENMR12+5) = '-R12/'
          LENMR12 = LENMR12 + 5
        ELSE IF (MODELR12(1:LENMR12).EQ.'MP2') THEN
            MODELR12(LENMR12+1:LENMR12+5) = '-R12/'
            LENMR12 = LENMR12 + 5
        ELSE 
            MODELR12(LENMR12+1:LENMR12+6) = '(R12)/'
            LENMR12 = LENMR12 + 6
        END IF

        I = 1 
        DO WHILE(I.LE.LEN(APROXR12) .AND. APROXR12(I:I).NE.' ')
          IF (LMAX.LT.LENMR12+1) 
     &       CALL QUIT('LMAX too small in CCSD_MODEL')
          LENMR12 = LENMR12 + 1
          MODELR12(LENMR12:LENMR12) = APROXR12(I:I)
          I = I + 1
        END DO
     
      ELSE
        MODELR12(1:5) = MODEL(1:5)
        LENMR12 = 5
      END IF

      RETURN
      END
C  /* Deck ccsd_iajb */
      SUBROUTINE CCSD_IAJB(XAIBJ,T1AM,LHTF,CCR12RSP,MKVAJKL,WORK,LWORK)
C
C     Written by Henrik Koch 27-Mar-1990.
C
C     Small modifications by Asger Halkier 22/5 - 1998 for extra
C     MO integrals needed for gradients and frozen core FOP.
C
C     Added calculation of V(alpha j,kl) for CC2-R12, H.Fliegl, C. Haettig
C
C     Added flag for computation of additional half transformed 
C     integrals needed for R12 response (e.g. r12 integrals with 
C     two auxiliary basis functions): CCR12RSP
C     C. Neiss, 10.12.2004
C     
#include <implicit.h>
#include <ccr12int.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
#include <priunit.h>
#include <dummy.h>
#include <maxorb.h>
#include <maxash.h>
#include <mxcent.h>
#include <aovec.h>
#include <iratdef.h>
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
#include <r12int.h>
      DIMENSION XAIBJ(*),T1AM(*),WORK(*),INDEXA(MXCORB)
      INTEGER IDUM,LUNITR12,LUNITR12_2
      LOGICAL LDUM,LHTF,MKVAJKL,CCR12RSP
      INTEGER KO2AM,YS2AM,KOFFH,KOFFD
      integer ilmorb(8)
      INTEGER IGABJ(8),IBASX(8),ICMO(8,8),IGLMRHS(8,8),NCMO(8),
     &        NGLMDS(8),KLAMDHS,KLAMDPS,KEND0,LWRK0,IGLMVIS(8,8)
      INTEGER IMAIJM(8,8),NMAIJM(8),IMATIJM(8,8),NMATIJM(8),
     &        IGAMSM(8,8),NGAMSM(8),IRGIJS(8,8),NRGIJS(8),
     &        IR1BASM(8,8),NR1BASM(8),IR2BASM(8,8),NR2BASM,
     &        IR1XBASM(8,8),NR1XBASM(8),IR2XBASM(8,8),IMATF(8,8),
     &        NMATF(8)
      INTEGER IMAKLM(8,8),NMAKLM(8)
C
      CHARACTER*5 FN3FOP
      CHARACTER*6 FN3VI, FN3FOP2
      CHARACTER*8 FN3SRT, FN3VI2, FNTOC
      CHARACTER*8 FILER12, FILER12_2
      CHARACTER*8 FILBACK
      CHARACTER*10 FILE_BACK 
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CCSD_IAJB')
C-----------------------------------------------------------
C     calculate some offsets and dimensions needed for R12 
C-----------------------------------------------------------
      KEND0 = 1
      IF (MKVAJKL) THEN
         TIMVAJKL = 0.0D0
   
celena
         IF (R12PRP) THEN
           DO ISYMAI = 1,NSYM
             ICOUN2 = 0
             DO  ISYMI = 1,NSYM
                ISYMA = MULD2H(ISYMAI,ISYMI)            
                ILMORB(ISYMI) = ICOUN2
                ICOUN2 = ICOUN2 + NBAS(ISYMI)*
     &                   (NORB1(ISYMI)-NRHFFR(ISYMI))
             ENDDO
           ENDDO
         ENDIF
celena
  
C        CALL CC_R12OFFSET(NR1ORB,NR1XORB,NR1BAS,NR1XBAS,NR2BAS,
C    &        NRGKL,NRXGKL,N2BST1,IR1ORB,IR1XORB,IR1BAS,IR1XBAS,IR2BAS,
C    &        IR2XBAS,IRGKL,IRXGKL,IAODIS1,NALPHAJ,IALPHAJ)
c
         IF (IANR12.EQ.2 .OR. IANR12.EQ.3) then
c          calculate some offsets and dimensions needed for Lambda
c          including active and inactive occupied molecular orbitals  

           CALL CC_R12OFFS23(IGLMRHS,IGLMVIS,NGLMDS,ICMO,NCMO,
     &                       IMAIJM,NMAIJM,IMAKLM,NMAKLM,
     &                       IMATIJM,NMATIJM,
     &                       IGAMSM,NGAMSM,IRGIJS,NRGIJS,
     &                       IR1BASM,NR1BASM,IR2BASM,NR2BASM,IR1XBASM,
     &                       NR1XBASM,IR2XBASM,IMATF,NMATF) 

           KLAMDHS = KEND0
           KLAMDPS = KLAMDHS + NGLMDS(1)
           KT1AM   = KLAMDPS + NGLMDS(1)
           KEND0   = KT1AM + NT1AMX
           LWRK0   = LWORK - KEND0
           IF (LWRK0.LT.0) THEN
             CALL QUIT('Insufficient work space in ccsd_iajb')
           END IF 
           CALL DZERO(WORK(KT1AM),NT1AMX)
           CALL LAMMATS(WORK(KLAMDPS),WORK(KLAMDHS),WORK(KT1AM),
     &                  1,.TRUE.,.FALSE.,
     &                  NGLMDS,IGLMRHS,IGLMVIS,ICMO,WORK(KEND0),LWRK0)
         END IF
      END IF
C-----------------------------------------
C     Initialize the XAIBJ integral array.
C-----------------------------------------
C
      IF (ONEAUX) THEN
         CALL DZERO(XAIBJ,NH2AM(ISYMOP))
      ELSE IF (U12INT .OR. R12SQR) THEN
C        Zero space for non-Hermitean integrals (WK/UniKA/04-11-2002).
         CALL DZERO(XAIBJ,NU2AM(ISYMOP))
      ELSE
         CALL DZERO(XAIBJ,NT2AM(ISYMOP))
      END IF

C----------------------------------------
C     Open files needed for CC-R12:
C----------------------------------------
      IF (LHTF) THEN
        LUNITR12 = -1
        IF (R12EOR.AND.CCR12SM) THEN
c       case for new correlation factor: need (ialpha|f12/r12|jbeta) on file
c                                        for initialisation if V-interm.
          FILER12  = FR12F12HTF
        ELSE
          FILER12  = FRHTF
        END IF
        CALL WOPEN2(LUNITR12,FILER12,64,0)
        IF (CCR12RSP.AND..NOT.CCR12SM) THEN
          LUNITR12_2 = -1
          FILER12_2  = FRHTF2
          CALL WOPEN2(LUNITR12_2,FILER12_2,64,0)
          CALL FLSHFO(LUPRI)
        END IF
      END IF

      IF (CCR12.AND.V12INT.AND.LHTF.AND.
     &           (IANR12.EQ.2.OR.IANR12.EQ.3)) THEN
        LUNITR12 = -1
        FILER12  = FGHTF
        CALL WOPEN2(LUNITR12,FILER12,64,0)
        LU44 = -1
        CALL WOPEN2(LU44,FCCGMNAB,64,0)
      END IF
C
C---------------------------------
C     Dynamic allocation of space.
C---------------------------------
C
      KLAMDP = KEND0
      IF (R12SQR) THEN
C        Read MO coefficients from GUMAT.n for n=1,2 (WK/UniKA/04-11-2002).
         KLAMDQ = KLAMDP + NLAMDT
         LU43 = -43 
         IF (COMBSS) THEN
            CALL GPOPEN(LU43,'GUMAT.2','UNKNOWN',' ','UNFORMATTED',
     &                  IDUM,LDUM)
         ELSE
            CALL GPOPEN(LU43,'GUMAT.1','UNKNOWN',' ','UNFORMATTED',
     &                  IDUM,LDUM)
         END IF
         REWIND(LU43)
         READ(LU43) NTOTGU
         READ(LU43) (WORK(KLAMDQ+I-1), I = 1, NTOTGU)
         IF ((R12EIN .AND. INTGAC .EQ. 4) .OR. (R12PRP) 
     *        ) THEN
            CALL GPCLOSE(LU43,'KEEP')
         ELSE
            CALL GPCLOSE(LU43,'KEEP')
         END IF
         KLAMDH = KLAMDQ + NTOTGU
      ELSE
         KLAMDQ = KLAMDP
         KLAMDH = KLAMDQ + NLAMDT
      END IF
      KEND1  = KLAMDH + NLAMDT
C
      KCMO   = KEND1
      KDNSHF = KCMO   + NLAMDS
      KFCKHF = KDNSHF + N2BAST
      KEND1  = KFCKHF + N2BAST

C
      IF (MKVAJKL) THEN
         KVAJKL = KEND1
         KEND1  = KVAJKL + NVAJKL(1)
         IF (R12PRP) THEN
            KXAJKL = KEND1
            KEND1  = KXAJKL+ NVAJKL(1)
         ENDIF
      END IF

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
      LUSIFC = -1
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',
     *            IDUMMY,.FALSE.)
      REWIND(LUSIFC)
C     Use LABEL (WK/UniKA/04-11-2002).
      CALL MOLLAB(LABEL,LUSIFC,LUPRI)
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
      IF ((FROIMP .OR. FROEXP) .AND. (.NOT. R12INT
     *    .AND. .NOT. R12EIN .AND. .NOT. U12INT) .AND.
     *    (R12TRA .OR. RELORB .OR. MP2) .OR. (FROIMP .AND. 
     *     R12PRP .AND. MKVAJKL)) THEN
C        Not needed for R12 integrals (WK/UniKA/04-11-2002).
!     *    (RELORB .OR. (CCFOP .AND. MP2))) THEN
!     Sonia: remove "FOP" condition to be able to do gradients MP2

         KCMO  = KEND1
         KFRIN = KCMO  + NLAMDS
         KFRGR = KFRIN + NT2FRO(1)
         KFRGR1= KFRGR + NFROVR(1)
         KEND1 = KFRGR1+ NFROVF(1)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient space in CCSD_IAJB')
         ENDIF
C
         CALL DZERO(WORK(KCMO),NLAMDS)
         IF (R12TRA .AND. .NOT. R12PRP) THEN
            CALL DZERO(WORK(KFRIN),NF2FRO(1))
         ELSE
            CALL DZERO(WORK(KFRIN),NT2FRO(1))
            CALL DZERO(WORK(KFRGR),NFROVR(1))
            CALL DZERO(WORK(KFRGR1),NFROVF(1))
         END IF
C
C----------------------------------------------
C     Calculate the FULL MO coefficient matrix.
C----------------------------------------------
C
         CALL CMO_ALL(WORK(KCMO),WORK(KEND1),LWRK1)
C
      ENDIF
C----------------------------------------------------
C     initialize V(alpha j,kl)
c----------------------------------------------------
      IF (MKVAJKL .AND. (.NOT. FNVAJKL .EQ. 'CCR12XAJKL'
     &    .AND. ( .NOT. FNVAJKL .EQ. 'CCR12QAJKL')
     &    .AND. ( .NOT. FNVAJKL .EQ. 'CCR12QIJAL')
     &    .AND. ( .NOT. FNVAJKL .EQ. 'CCR12UAJKL')
     &    .AND. ( .NOT. FNVAJKL .EQ. 'CCR12UIJAL') )) THEN
         DTIME = SECOND()
         IOPT = 1 
         CALL CC_R12MKVAMKL0(WORK(KVAJKL),NVAJKL(1),IOPT,WORK(KLAMDH),1,
     &        WORK(KEND1),LWRK1) 

         TIMVAJKL = TIMVAJKL + ( SECOND() - DTIME )
      ELSEIF ((MKVAJKL .AND. FNVAJKL .EQ. 'CCR12UIJAL') 
     &        .OR. (MKVAJKL .AND. FNVAJKL .EQ. 'CCR12QAJKL')
     &        .OR. (MKVAJKL .AND. FNVAJKL .EQ. 'CCR12QIJAL')
     &        .OR. (MKVAJKL .AND. FNVAJKL .EQ. 'CCR12UAJKL')
     &        .OR. (MKVAJKL .AND. FNVAJKL .EQ. 'CCR12XAJKL')) THEN
         CALL DZERO(WORK(KVAJKL),NVAJKL(1))
         CALL DZERO(WORK(KXAJKL),NVAJKL(1))
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
               KRECNR = KEND1
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
                  IF (NOAUXB) THEN
                     IDUM = 1
                     CALL IJKAUX(IDEL,IDUM,IDUM,IDUM)
                  END IF
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
               IF (U21INT) KSCR1 = KSCR1 + NDISAO(ISYDIS)
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
               IOFFU21 = NDISAO(ISYDIS)
               CALL DZERO(WORK(KXINT),2*NDISAO(ISYDIS))
               CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND2),LWRK2,
     *                     WORK(KRECNR),DIRECT)
C
C-----------------------------------------
C              compute the AO-Fock matrix:
C-----------------------------------------
C
C              Not needed for R12 part (WK/UniKA/28-04-2003).
               IF (.NOT. R12TRA) 
     *         CALL CC_AOFOCK(WORK(KXINT),WORK(KDNSHF),WORK(KFCKHF),
     *                        WORK(KEND2),LWRK2,IDEL,ISYMD,.FALSE.,
     *                        DUMMY,1)
C
C-----------------------------------------------
C              Calculate integrals (cJ|dk)
C              needed for frozen core gradients.
C-----------------------------------------------
C
C             Modified for R12 method (WK/UniKA/04-11-2002).
C             Modified (RELORB .OR. (CCFOP .AND. MP2)) for MP2 frozen-core gradients
C             Sonia
              IF ((FROIMP .OR. FROEXP) .AND. (.NOT. R12INT
     *            .AND. .NOT. R12EIN .AND. .NOT. U12INT) .AND.
     *            (R12TRA .OR. RELORB .OR. MP2)) THEN
C
                  IF (ONEAUX) THEN
                     CALL CC_FRCR12(WORK(KFRIN),WORK(KXINT),WORK(KCMO),
     *                              WORK(KEND2),LWRK2,IDEL,ISYMD)
                  ELSE
                     CALL CC_FRCOIN(WORK(KFRIN),WORK(KXINT),WORK(KCMO),
     *                              WORK(KEND2),LWRK2,IDEL,ISYMD)
                     IF (R12PRP .AND. R12TRA) THEN
                         CALL CC_FRCOGR(WORK(KFRGR),WORK(KXINT),
     *                        WORK(KCMO),WORK(KEND2),LWRK2,IDEL,ISYMD)
                         CALL CC_FRCOGR1(WORK(KFRGR1),WORK(KXINT),
     *                        WORK(KCMO),WORK(KEND2),LWRK2,IDEL,ISYMD)
                     ENDIF 

                  END IF
               END IF
C
C-----------------------------------------------------------------------
C              For CC-R12 with Ansatz 2 calculate two-index transformed
C              coulomb integrals (M alpha | N beta), where M,N are 
C              frozen and active occupied orbitals; 
C              integrals stored on file FCCGMNAB
C-----------------------------------------------------------------------
C
chf
               IF (CCR12.AND.V12INT.AND.LHTF.AND.
     &             MKVAJKL.AND.(IANR12.EQ.2.OR.IANR12.EQ.3)) THEN
                 CALL CC_R12MKGMNAB(WORK(KXINT),WORK(KLAMDHS),1,IDEL,
     &                              ISYMD,IGLMRHS,NGLMDS,LU44,
     &                              FCCGMNAB,WORK(KEND2),LWRK2)
               END IF
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
                  IF (U21INT) KOFFT = KOFF1  + NDISAO(ISYDIS)
C                 Use KLAMDQ instead of KLAMDP (WK/UniKA/04-11-2002).
                  IF (FNVAJKL .EQ. 'CCR12QIJAL') THEN
                      KOFF2 = KLAMDQ + ILMORB(ISYMJ)
                     
                  ELSE
                      KOFF2 = KLAMDQ + ILMRHF(ISYMJ)
                  ENDIF 
                  KOFF6 = KLAMDP + ILMRHF(ISYMJ)
C
                  IF (LWRK2 .LT. 2*NNBST(ISYMAB)*NRHF(ISYMJ)) THEN
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
                  CALL DGEMM('N','N',NNBST(ISYMAB),NRHF(ISYMJ),
     *                       NBAS(ISYMG),ONE,WORK(KOFF1),NNBSAB,
     *                       WORK(KOFF2),NBASG,ZERO,WORK(KEND2),
     *                       NNBSAB)
                  IF ((ONEAUX .OR. R12PRP) .AND. U21INT) THEN
                    KENDT = KEND2 + NNBST(ISYMAB)*NRHF(ISYMJ)
                    CALL DGEMM('N','N',NNBST(ISYMAB),NRHF(ISYMJ),
     *                          NBAS(ISYMG),ONE,WORK(KOFFT),NNBSAB,
     *                          WORK(KOFF2),NBASG,ZERO,WORK(KENDT),
     *                          NNBSAB)
                  ELSE IF (ONEAUX .AND. R12SQR) THEN
                     KENDT = KEND2 + NNBST(ISYMAB)*NRHF(ISYMJ)
                     CALL DGEMM('N','N',NNBST(ISYMAB),NRHF(ISYMJ),
     *                          NBAS(ISYMG),ONE,WORK(KOFF1),NNBSAB,
     *                          WORK(KOFF6),NBASG,ZERO,WORK(KENDT),
     *                          NNBSAB)
                  END IF
C------------------------------------------------------------------
C                 Transform integrals and add to the result vector.
C------------------------------------------------------------------
C
                  IF (CCSDT .OR. CCPT .OR. CHOPT .OR.
     *                CCP3 .OR. (CCRT.OR.CCR3.OR.CCR1A.OR.CCR1B)) THEN
                      LUTOC = -1
                      FNTOC = 'CCSDT_OC'
                      CALL WOPEN2(LUTOC,FNTOC,64,0)
                  ENDIF
C
                  IF (ONEAUX) THEN
                     KOFF4  = IH2AM(ISYMAI,ISYMBJ) + 1
                  ELSE IF (U12INT .OR. R12SQR ) THEN
C                    KOFF4 for non-Hermitean integrals (WK/UniKA/04-11-2002).
                     KOFF4  = IU2AM(ISYMAI,ISYMBJ) + 1
                  ELSE
                     KOFF4  = IT2AM(ISYMAI,ISYMBJ) + 1
                  END IF
C
                  CALL CCSD_AIBJ2(WORK(KEND2),XAIBJ(KOFF4),WORK(KLAMDP),
     *                            WORK(KLAMDH),WORK(KSCR1),WORK(KSCR2),
     *                            IDEL,ISYMD,ISYMJ,ISYMAB,LUTOC,FNTOC,
     *                            .FALSE.,LUNITR12,FILER12,LUNITR12_2,
     *                            FILER12_2,LHTF,CCR12RSP)
                  IF (ONEAUX .AND. U21INT) THEN
                  CALL CCSD_AIBJ2(WORK(KENDT),XAIBJ(KOFF4),WORK(KLAMDP),
     *                            WORK(KLAMDH),WORK(KSCR1),WORK(KSCR2),
     *                            IDEL,ISYMD,ISYMJ,ISYMAB,LUTOC,FNTOC,
     *                            .TRUE.,LUNITR12,FILER12,LUNITR12_2,
     *                            FILER12_2,LHTF,CCR12RSP)
                  ELSE IF (ONEAUX .AND. R12SQR) THEN
                  CALL CCSD_AIBJ2(WORK(KENDT),XAIBJ(KOFF4),WORK(KLAMDQ),
     *                            WORK(KLAMDH),WORK(KSCR1),WORK(KSCR2),
     *                            IDEL,ISYMD,ISYMJ,ISYMAB,LUTOC,FNTOC,
     *                            .FALSE.,LUNITR12,FILER12,LUNITR12_2,
     *                            FILER12_2,LHTF,CCR12RSP)
                  END IF
C---------------------------------------------------------
C                 compute contributions to V(alpha j,kl)
C---------------------------------------------------------
                  IF (MKVAJKL) THEN
                   DTIME = SECOND()
                   IF (MBAS1(ISYMG).GT.0 .OR. NRHF(ISYMJ).GT.0) THEN
                     IBASX(1) = 0
                     DO ISYM = 2, NSYM
                       IBASX(ISYM) = IBASX(ISYM-1)+MBAS2(ISYM-1)
                     END DO
                     KGABJD = KEND2
                     KEND3 = KGABJD + NNBST(ISYMAB)*NRHF(ISYMJ)
                     IF (U21INT) THEN
                       KTABJD = KENDT
                       KEND3 = KTABJD + NNBST(ISYMAB)*NRHF(ISYMJ)
                     END IF
                     LWRK3 = LWORK - KEND3
 
                     IF (LWRK3 .LT. 0) THEN
                        CALL QUIT('Insufficient space in CCSD_IAJB')
                     END IF

                     KOFF5 = KXINT + IDSAOG(ISYMG,ISYDIS)
                     IF(U21INT) KOFF6 = KOFF5  + NDISAO(ISYDIS)

                     IF      (IANR12.EQ.1 .AND. (.NOT.R12PRP)) THEN
                        FILBACK = FNBACK
celena
                     ELSEIF (R12PRP) THEN
                           IF (FNVAJKL .EQ. 'CCR12VIJAL') THEN
                               FILE_BACK = FV12BACK
                           ELSE IF (FNVAJKL .EQ. 'CCR12VAJKL') THEN 
                               FILBACK = FNBACK
                           ELSE IF (FNVAJKL .EQ. 'CCR12BIJAL') THEN
                               FILE_BACK = FT12BACK
                           ELSE IF (FNVAJKL .EQ. 'CCR12BAJKL') THEN
                               FILBACK = FNBACK
                           ELSE IF (FNVAJKL .EQ. 'CCR12QAJKL') THEN
                               FILBACK = FNBACK
                           ELSE IF (FNVAJKL .EQ. 'CCR12QIJAL') THEN
                               FILBACK = FNBACK
                           ELSE IF (FNVAJKL .EQ. 'CCR12UAJKL') THEN
                               FILE_BACK = FU12BACK
                           ELSE IF (FNVAJKL .EQ. 'CCR12UIJAL') THEN
                               FILE_BACK = FQ12BACK
                           ELSE IF (FNVAJKL .EQ. 'CCR12XAJKL') THEN
                               FILBACK = FNBACK
                           ENDIF
celena
                     ELSE IF (IANR12.EQ.2) THEN
c                       FILBACK = FRHTF
                        FILBACK = FNBACK2
                     ELSE IF (IANR12.EQ.3) THEN
                        IDELTA = IDEL - IBAS(ISYMD)
                        IF (IDELTA.LE.MBAS1(ISYMD)) THEN
                          FILBACK = FNBACK
                        ELSE 
c                         FILBACK = FRHTF
                          FILBACK = FNBACK2
                        END IF
                     ELSE
                        WRITE(LUPRI,*) 'IANR12 = ',IANR12
                        CALL QUIT('Illegal IANR12.')
                     END IF
                     IF (FNVAJKL .EQ. 'CCR12QIJAL') THEN
                        CALL R12MKVAMKL(FILBACK,WORK(KGABJD),
     &                       WORK(KTABJD),WORK(KVAJKL),
     &                       WORK(KLAMDQ),1,WORK(KLAMDHS),WORK(KLAMDPS),
     &                       WORK(KOFF5),WORK(KOFF6),
     &                       IDEL,ISYMD,ISYMJ,ISYMAB,ISYMG,
     &                       WORK(KSCR1),IBASX,IGLMRHS,NGLMDS,
     &                       WORK(KEND3),LWRK3)
                     ELSEIF (FNVAJKL .EQ. 'CCR12BIJAL' .OR.
     &                       FNVAJKL .EQ. 'CCR12UAJKL' .OR.
     &                       FNVAJKL .EQ. 'CCR12VIJAL' .OR.
     &                       FNVAJKL .EQ. 'CCR12UIJAL') THEN
                        CALL R12MKVAMKL(FILE_BACK,WORK(KGABJD),
     &                       WORK(KTABJD),WORK(KVAJKL),
     &                       WORK(KLAMDQ),1,WORK(KLAMDHS),WORK(KLAMDPS),
     &                       WORK(KOFF5),WORK(KOFF6),
     &                       IDEL,ISYMD,ISYMJ,ISYMAB,ISYMG,
     &                       WORK(KSCR1),IBASX,IGLMRHS,NGLMDS,
     &                       WORK(KEND3),LWRK3)
                     ELSE 
                        CALL R12MKVAMKL(FILBACK,WORK(KGABJD),
     &                       WORK(KTABJD),WORK(KVAJKL),
     &                       WORK(KLAMDH),1,WORK(KLAMDHS),WORK(KLAMDPS),
     &                       WORK(KOFF5),WORK(KOFF6),
     &                       IDEL,ISYMD,ISYMJ,ISYMAB,ISYMG,
     &                       WORK(KSCR1),IBASX,IGLMRHS,NGLMDS,
     &                       WORK(KEND3),LWRK3)
                        IF (IANR12.EQ.3 .AND. R12CBS) THEN
                          !once more in this case...
                          IDELTA = IDEL - IBAS(ISYMD)
                          IF (IDELTA.LE.MBAS1(ISYMD)) THEN
                            IANR12 = 2
                            FILBACK = FNBACK2
                            CALL R12MKVAMKL(FILBACK,WORK(KGABJD),
     &                           WORK(KTABJD),WORK(KVAJKL),
     &                           WORK(KLAMDH),1,WORK(KLAMDHS),
     &                           WORK(KLAMDPS),WORK(KOFF5),WORK(KOFF6),
     &                           IDEL,ISYMD,ISYMJ,ISYMAB,ISYMG,
     &                           WORK(KSCR1),IBASX,IGLMRHS,NGLMDS,
     &                           WORK(KEND3),LWRK3)
                            IANR12 = 3
                          END IF
                        END IF
                     ENDIF 
                   END IF
                   TIMVAJKL = TIMVAJKL + ( SECOND() - DTIME )
                  END IF
C
C----------------------------------------------------------------------
C                  Construct I(kd,c) for fixed alpha.
C                  Not needed for R12 integrals (WK/UniKA/04-11-2002).
C----------------------------------------------------------------------
C
                  IF (CCSDT .OR. CCPT .OR. CHOPT .OR.
     *                CCP3 .OR. (CCRT.OR.CCR3.OR.CCR1A.OR.CCR1B)) THEN
                     CALL WCLOSE2(LUTOC,FNTOC,'KEEP')
                  ENDIF
C
                  IF (.NOT. R12TRA .AND. (CCSDT.OR.(CCPT.OR.CCP3).OR.
     *                (CCRT.OR.CCR3.OR.CCR1A.OR.CCR1B.OR.CHOPT))) THEN

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
                     LU3SRT = -1
                     FN3SRT = 'CC3_SORT'
                     CALL WOPEN2(LU3SRT,FN3SRT,64,0)
C
                     CALL CCSD_AIBJ3(WORK(KOFF5),WORK(KINT3),
     *                               WORK(KINT4),WORK(KLAMDP),
     *                               WORK(KLAMDH),WORK(KSCR1),
     *                               WORK(KSCR2),WORK(KSCR3),
     *                               IDEL,ISYMD,ISYMG,ISYMAB,
     *                               LU3SRT,FN3SRT)
C
                     CALL WCLOSE2(LU3SRT,FN3SRT,'KEEP')
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
      KEND1 = KENDSV
      LWRK1 = LWRKSV
C
c-------------------------------------
C     write AO-Fock matrix to file:
C-------------------------------------
C
      IF (.NOT. R12TRA) THEN
         LUFCK = -1
         CALL GPOPEN(LUFCK,'CC_FCKREF','UNKNOWN',' ','UNFORMATTED',
     *               IDUMMY,.FALSE.)
         REWIND(LUFCK)
         WRITE(LUFCK)(WORK(KFCKHF + I-1),I = 1,N2BST(ISYMOP))
         CALL GPCLOSE(LUFCK,'KEEP' )
C
         IF (IPRINT .GT.150) THEN
            CALL AROUND( 'Fock AO matrix for reference state:' )
            CALL CC_PRFCKAO(WORK(KFCKHF),1)
         ENDIF                  
      ENDIF                  
C
      IF (ANAAOD) THEN
        CALL AROUND('Analysis of integral distributions')
C
        WRITE(LUPRI,'(10X,/,A,D10.5)') 'Threshold in analysis : ',
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
          IF (ONEAUX) THEN
            DO 250 ISYMBJ = 1,NSYM
              ISYMAI = ISYMBJ
              KOFF   = IH2AM(ISYMAI,ISYMBJ) + 1
              NTOTAI = NH1AM(ISYMAI)
              CALL OUTPAK(XAIBJ(KOFF),NTOTAI,1,LUPRI)
              KOFF   = KOFF + NTOTAI * (NTOTAI + 1) / 2
              NTOTBJ = NG1AM(ISYMAI)
              CALL OUTPUT(XAIBJ(KOFF),1,NTOTAI,1,NTOTBJ,
     &                    NTOTAI,NTOTBJ,1,LUPRI)
  250       CONTINUE
          ELSE IF (U12INT .OR. R12SQR) THEN
C           Output of non-Hermitean integrals (WK/UniKA/04-11-2002).
            DO 251 ISYMBJ = 1,NSYM
              ISYMAI = ISYMBJ
              KOFF   = IU2AM(ISYMAI,ISYMBJ) + 1
              NTOTAI = NT1AM(ISYMAI)
              CALL OUTPUT(XAIBJ(KOFF),1,NTOTAI,1,NTOTAI,
     &                    NTOTAI,NTOTAI,1,LUPRI)
  251       CONTINUE
          ELSE
            DO 252 ISYMBJ = 1,NSYM
              ISYMAI = ISYMBJ
              KOFF   = IT2AM(ISYMAI,ISYMBJ) + 1
              NTOTAI = NT1AM(ISYMAI)
              CALL OUTPAK(XAIBJ(KOFF),NTOTAI,1,LUPRI)
  252       CONTINUE
          END IF
        ENDIF
      END IF

C     -----------------------------------------
C     write V(alpha j,kl) to disk
C     -----------------------------------------
      IF (MKVAJKL) THEN
         DTIME = SECOND() 

         IF (DEBUG) THEN
           CALL CC_R12MKVIJKL(WORK(KVAJKL),1,WORK(KLAMDH),1,
     &                        WORK(KEND1),LWRK1,.FALSE.,DUMMY,DUMMY)
         END IF

         LUVAJKL = -1
         CALL GPOPEN(LUVAJKL,FVAJKL,'UNKNOWN',' ','UNFORMATTED',
     &               IDUMMY,.FALSE.)
         REWIND(LUVAJKL)
         WRITE(LUVAJKL) (WORK(KVAJKL+I-1), I = 1,NVAJKL(1))
         CALL GPCLOSE(LUVAJKL,'KEEP')

C     Compute Y(a,j,k,l) for MP2-R12 first order properties (Y=B,V,X)       
         IF (R12PRP) THEN
            IF (FNVAJKL .EQ. 'CCR12QAJKL') THEN 
               CALL CC_R12MKXAJKL(WORK(KVAJKL),WORK(KCMO),WORK(KEND1),
     &                            LWRK1,.true.)
               LUVAJKL = -1
               CALL GPOPEN(LUVAJKl,FNVAJKL,'UNKNOWN',' ','UNFORMATTED',
     &                     IDUMMY,.FALSE.)
               REWIND(LUVAJKL)
               WRITE(LUVAJKL) (WORK(KVAJKL+I-1), I = 1,NVAJKL(1))
               CALL GPCLOSE(LUVAJKL,'KEEP')

             ELSE 
               CALL CC_R12MKXAJKL(WORK(KVAJKL),WORK(KCMO),WORK(KEND1),
     &                            LWRK1,.false.)
             ENDIF

             IF (FROIMP) THEN
               IF (FNVAJKL .EQ. 'CCR12QAJKL') THEN
                  CALL CC_R12MKXIJKL(WORK(KVAJKL),WORK(KCMO),WORK(KEND1)
     &                               ,LWRK1,.true.)
                ELSE
                 CALL CC_R12MKXIJKL(WORK(KVAJKL),WORK(KCMO),WORK(KEND1)
     &                               ,LWRK1,.false.)
                ENDIF
            END IF
         END IF
         

         TIMVAJKL = TIMVAJKL + ( SECOND() - DTIME )
          
         WRITE(LUPRI,'(1X,A)')
     &     'Computation of V^aj_kl intermediate done'
         WRITE(LUPRI,'(/1X,A,F7.2,A)')
     &     '>>> Time used for V^aj_kl is ',TIMVAJKL,' seconds'
         WRITE(LUPRI,*)
      END IF
C
C-------------------------------------
C     Write integrals (cJ|dk) to disk.
C-------------------------------------
C
C     Modified for R12 method (WK/UniKA/04-11-2002).
C     Modified for MP2 frozen-core geometry opt. Sonia 2002
      IF ((FROIMP .OR. FROEXP) .AND. (.NOT. R12INT
     *    .AND. .NOT. R12EIN .AND. .NOT. U12INT) .AND.
     *    (R12TRA .OR. RELORB .OR. MP2)) THEN
C
         LUCJDK = -1
         CALL GPOPEN(LUCJDK,'INCJDK','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND(LUCJDK)
         WRITE(LUCJDK) (WORK(KFRIN+I-1), I = 1,NT2FRO(1))
         CALL GPCLOSE(LUCJDK,'KEEP')
        IF (R12PRP) THEN
            LUCJDK = -1
            CALL GPOPEN(LUCJDK,'INCJDA','UNKNOWN',' ','UNFORMATTED',
     &               IDUMMY,.FALSE.)
            REWIND(LUCJDK)
            WRITE(LUCJDK) (WORK(KFRGR+I-1), I = 1,NFROVR(1))
            CALL GPCLOSE(LUCJDK,'KEEP')

            LUCJDK = -1
            CALL GPOPEN(LUCJDK,'INCJDI','UNKNOWN',' ','UNFORMATTED',
     &               IDUMMY,.FALSE.)
            REWIND(LUCJDK)
            WRITE(LUCJDK) (WORK(KFRGR1+I-1), I = 1,NFROVF(1))
            CALL GPCLOSE(LUCJDK,'KEEP')
         ENDIF 
      ENDIF
C
      IF (.NOT.R12TRA.AND.(CCSDT.OR.(CCPT.OR.CCP3).OR.
     *    (CCRT.OR.CCR3.OR.CCR1A.OR.CCR1B.OR.CHOPT))) THEN
C
C------------------------------------
C        Sort integrals (kc,d alpha).
C------------------------------------
C
         LU3SRT  = -1
         LU3VI   = -1
         LU3VI2  = -1
         LU3FOP  = -1
         LU3FOP2 = -1
         FN3SRT  = 'CC3_SORT'
         FN3VI   = 'CC3_VI'
         FN3VI2  = 'CC3_VI12'
         FN3FOP  = 'PTFOP'
         FN3FOP2 = 'PTFOP2'
         CALL WOPEN2(LU3SRT,FN3SRT,64,0)
         CALL WOPEN2(LU3VI,FN3VI,64,0)
         CALL WOPEN2(LU3VI2,FN3VI2,64,0)
         CALL WOPEN2(LU3FOP,FN3FOP,64,0)
         CALL WOPEN2(LU3FOP2,FN3FOP2,64,0)
C
         ISYINT = ISYMOP
         CALL CC3_SORT1(WORK,LWORK,1,ISYINT,LU3SRT,FN3SRT,
     *                  LU3VI,FN3VI,LU3VI2,FN3VI2,LU3FOP,FN3FOP,
     *                  LU3FOP2,FN3FOP2)
C
         CALL WCLOSE2(LU3SRT,FN3SRT,'KEEP')
         CALL WCLOSE2(LU3VI,FN3VI,'KEEP')
         CALL WCLOSE2(LU3VI2,FN3VI2,'KEEP')
         CALL WCLOSE2(LU3FOP,FN3FOP,'KEEP')
         CALL WCLOSE2(LU3FOP2,FN3FOP2,'KEEP')
C
      ENDIF
C  
      IF (LHTF) THEN
        CALL WCLOSE2(LUNITR12,FILER12,'KEEP')
        IF (CCR12RSP) THEN
          CALL WCLOSE2(LUNITR12_2,FRHTF2,'KEEP')
        END IF
      END IF
C
      IF (CCR12.AND.V12INT.AND.LHTF.AND.
     &        (IANR12.EQ.2.OR.IANR12.EQ.3)) THEN
        CALL WCLOSE2(LU44,FCCGMNAB,'KEEP')
      END IF
C  
      CALL QEXIT('CCSD_IAJB')

      RETURN
      END
C  /* Deck ccsd_aibj2 */
      SUBROUTINE CCSD_AIBJ2(XINT,XAIBJ,XLAMDP,XLAMDH,
     *                      SCR1,SCR2,IDEL,ISYMD,ISYMJ,ISYMAB,
     *                      LUFILE,FNFILE,ANTISYM,
     *                      LUNITR12,FILER12,LUNITR12_2,FILER12_2,
     *                      LHTF,CCR12RSP)
C
C     Written by Henrik Koch 27-Mar-1990.
C
#include <implicit.h>
      INTEGER LU43,LUNITR12,LUNITR12_2
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION XINT(*),XAIBJ(*), SCR1(*),SCR2(*)
      DIMENSION XLAMDP(*),XLAMDH(*)
#include <priunit.h>
#include <ccinftap.h>
#include <ccorb.h>
#include <r12int.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <ccfop.h>
C
      LOGICAL ANTISYM,LHTF,CCR12RSP
      CHARACTER*(*) FNFILE,FILER12,FILER12_2
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CCSD_AIBJ2')
C
      IF (ANTISYM) THEN
         FACDG = -ONE
      ELSE
         FACDG = ONE
      END IF
C
      DO 100 J = 1,NRHF(ISYMJ)
C
         KOFF1 = NNBST(ISYMAB)*(J-1) + 1
C
         IF (ANTISYM) THEN
            CALL CCSD_ASYMSQ(XINT(KOFF1),ISYMAB,SCR1,0,0)
         ELSE
            CALL CCSD_SYMSQ(XINT(KOFF1),ISYMAB,SCR1)
         END IF
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
            CALL DGEMM('T','N',NBAS(ISYMB),NRHF(ISYMI),NBAS(ISYMA),
     *                 FACDG,SCR1(KOFF1),NBASA,XLAMDP(KOFF2),
     *                 NBASA,ZERO,SCR2(KOFF3),NBASB)
C
            KOFF3 = KOFF3 + NBAS(ISYMB)*NRHF(ISYMI)
C
  110    CONTINUE
C
         IF (LHTF) THEN
           NSCR1 = NBAST*NBAST
           CALL CC_R12WHTF(SCR2,IDEL,ISYMD,J,ISYMJ,ISYMAB,CCR12RSP,
     &                     LUNITR12,FILER12,LUNITR12_2,FILER12_2,
     &                     SCR1,NSCR1)
         END IF
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
            NBASB = MAX(NBAS(ISYMB),1)
C       
            IF (ONEAUX) THEN
              KOFF3 = IH1AM(ISYMA,ISYMI) + 1
              NVIRA = MAX(NORB1(ISYMA),1)
              CALL DGEMM('T','N',NORB1(ISYMA),NRHF(ISYMI),NBAS(ISYMB),
     *                   ONE,XLAMDH(KOFF1),NBASB,SCR2(KOFF2),
     *                   NBASB,ZERO,SCR1(KOFF3),NVIRA)
            ELSE
              KOFF3 = IT1AM(ISYMA,ISYMI) + 1
              NVIRA = MAX(NVIR(ISYMA),1)
              CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMB),
     *                   ONE,XLAMDH(KOFF1),NBASB,SCR2(KOFF2),
     *                   NBASB,ZERO,SCR1(KOFF3),NVIRA)
            END IF
            KOFF2 = KOFF2 + NBAS(ISYMB)*NRHF(ISYMI)
C       
  120    CONTINUE
c------------------------------------------------------------
CHF write and grep out here occupied g_ijkdelta integrals 
c------------------------------------------------------------
C
C------------------------------------------
C        Write out integrals used in CCSDT.
C------------------------------------------
C
         IF (CCSDT.OR.CCPT.OR.CCP3.OR.CCRT
     *       .OR.CCR3.OR.CCR1A.OR.CCR1B .OR. CHOPT) THEN
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
               CALL PUTWA2(LUFILE,FNFILE,SCR1,IOFF,NT1AM(ISYMCK))
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
         IF (ONEAUX) THEN
          DO 131 B = 1, NORB1(ISYMB)
            NBJ = IH1AM(ISYMB,ISYMJ) + NORB1(ISYMB)*(J-1) + B
            NTOTAI = NBJ
            KOFF1 = NBJ*(NBJ - 1)/2 + 1
            KOFF2 = ILMVIR(ISYMB) + NBAS(ISYMD)*(B-1) + IDEL
     *              - IBAS(ISYMD)
            CALL DAXPY(NTOTAI,XLAMDH(KOFF2),SCR1,1,XAIBJ(KOFF1),1)
  131     CONTINUE
          NTOTAI = NH1AM(ISYMAI)
          KOFF0 = NTOTAI * (NTOTAI + 1) / 2 + 1
          DO 132 B = 1, NORB2(ISYMB)
            KKB = B + NORB1(ISYMB)
            NBJ = IG1AM(ISYMB,ISYMJ) + NORB2(ISYMB)*(J-1) + B
            KOFF1  = NTOTAI*(NBJ - 1) + KOFF0
            KOFF2 = ILMVIR(ISYMB) + NBAS(ISYMD)*(KKB-1) + IDEL
     *              - IBAS(ISYMD)
            CALL DAXPY(NTOTAI,XLAMDH(KOFF2),SCR1,1,XAIBJ(KOFF1),1)
  132     CONTINUE
         ELSE
          DO 130 B = 1, NVIR(ISYMB)
            NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J-1) + B
            IF (ISYMAI .EQ. ISYMBJ .AND. .NOT.
     &                        (U12INT .OR. R12SQR)) THEN
               NTOTAI = NBJ
               KOFF1 = NBJ*(NBJ - 1)/2 + 1
            ELSE
               NTOTAI = NT1AM(ISYMAI)
               KOFF1  = NTOTAI*(NBJ - 1) + 1
            ENDIF
            KOFF2 = ILMVIR(ISYMB) + NBAS(ISYMD)*(B-1) + IDEL
     *              - IBAS(ISYMD)
            CALL DAXPY(NTOTAI,XLAMDH(KOFF2),SCR1,1,XAIBJ(KOFF1),1)
  130     CONTINUE
         END IF
C
  100 CONTINUE
C
      CALL QEXIT('CCSD_AIBJ2')
C
      RETURN
      END
C  /* Deck ccsd_aibj3 */
      SUBROUTINE CCSD_AIBJ3(XINT,XINT3,XINT4,XLAMDP,XLAMDH,SCR1,SCR2,
     *                      SCR3,IDEL,ISYDEL,ISYMG,ISYMAB,LUFILE,FNFILE)
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
      CHARACTER*(*) FNFILE
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CCSD_AIBJ3')
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
         CALL PUTWA2(LUFILE,FNFILE,XINT3,IOFF,LENGTH)
      ENDIF
C
      CALL QEXIT('CCSD_AIBJ3')

      RETURN
      END
C  /* Deck inidat */
      BLOCK DATA INIDAT
C
C     Initialize MULD2H in common block /CCORB/
C
#include <ccorb.h>
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
C     Frozen orbital bug-fix, tbp July 2003.
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
      DIMENSION WORK(LWORK)
C
      EXTERNAL INIDAT
C
#include <maxorb.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <inftap.h>
#include <symsq.h>
#include <ccisao.h>
#include <r12int.h>
#include <cc3t3d.h>
Cholesky
#include <dccorb.h>
#include <dccsdsym.h>
#include <ccisvi.h>
Cholesky
C
      INTEGER NMATAK(8)
      LOGICAL FIRST
      DATA FIRST /.TRUE./
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CCSD_INIT1')
C
C-------------------------------------
C     Read in information from sirius.
C-------------------------------------
C
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
C
C     LABEL is used (WK/UniKA/04-11-2002).
      CALL MOLLAB(LABEL,LUSIFC,LUPRI)
      READ (LUSIFC) NSYM, NORBTS, NBAST, NLAMDS, (NRHFS(I),I=1,NSYM),
     &              (NORBS(I),I=1,NSYM), (NBAS(I),I=1,NSYM), PDUM, EDUM
cccms IF (FIRST .AND. LGLO) THEN
c         DO ISYM = 1 , NSYM
c            NORB1(ISYM) = NORB1(ISYM) + NRHFS(ISYM)
c         ENDDO
c         FIRST = .FALSE.
c     END IF
C
      IF (FREEZE) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) ' I am freezing!'
C
         KFOCKD = 1
         KFCS   = KFOCKD + NORBTS
         KFVS   = KFCS   + NSYM
         KEND1  = KFVS   + NSYM
         LEND1  = LWORK  - KEND1
C
         READ (LUSIFC) (WORK(KFOCKD+I-1), I=1,NORBTS)
C
         CALL CC_FREEZER(WORK(KFOCKD),NORBTS,WORK(KFCS),WORK(KFVS),
     *                   WORK(KEND1),LEND1,LABEL)
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
      NRHFTB = 0
C
      ICOUN1 = 0
      ICOUN2 = 0
      IOFF   = 0
      JOFF   = 0  ! ALFREDO OK? JOFF was not initialized in your ccsd_energy.F
C
      DO 5 ISYM = 1,NSYM
C
         NVIRS(ISYM) = NORBS(ISYM) - NRHFS(ISYM)
C
         NRHF(ISYM)  = NRHFS(ISYM) - NRHFFR(ISYM)
         NVIR(ISYM)  = NVIRS(ISYM) - NVIRFR(ISYM)
         NORB(ISYM)  = NRHF(ISYM)  + NVIR(ISYM)
C
         XRHF(ISYM) = 1.0D0 * NRHF(ISYM)
         XVIR(ISYM) = 1.0D0 * NVIR(ISYM)
C
         IF (LABEL.EQ.'TRCCINT ') THEN
           NRHFA(ISYM)  = NRHF(ISYM) 
           NRHFSA(ISYM) = NRHFS(ISYM) 
           NRHFB(ISYM)  = NRHF(ISYM) + NRXR12(ISYM)
           NRHFSB(ISYM) = NRHFS(ISYM) + NRXR12(ISYM)
           NRHFTB = NRHFTB + NRHFB(ISYM)
         ELSE
           NRHFA(ISYM)  = NRHF(ISYM) - NRXR12(ISYM)
           NRHFSA(ISYM) = NRHFS(ISYM) - NRXR12(ISYM)
           NRHFB(ISYM)  = NRHF(ISYM)
           NRHFSB(ISYM) = NRHFS(ISYM)
           NRHFTB = NRHFTB + NRHFB(ISYM)
         END IF
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
Cholesky
        DO I = 1,NVIR(ISYM)

           JOFF = JOFF + 1
           ISVI(JOFF) = ISYM

        ENDDO
Cholesky
    5 CONTINUE
C
      NVIRT  = NORBT  - NRHFT
      NVIRTS = NORBTS - NRHFTS
C
      IF (IPRINT .GT. 20) THEN
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
      END IF
C
C--------------------------------------------------------
C     Construct implicitly frozen matrices.
C     (Matrices for FROEXP constructed in input routine.)
C--------------------------------------------------------
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
               KFRVIR(I,ISYM) = JORB
   52       CONTINUE
C
   50    CONTINUE
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
Chol
         XT1AM(ISYMAI) = 0.0D0
         XT1AO(ISYMAI) = 0.0D0
Chol
         NH1AM(ISYMAI) = 0
         NG1AM(ISYMAI) = 0
celena
         NT1VM(ISYMAI) = 0
celena
         DO 200 ISYMI = 1,NSYM
            ISYMA = MULD2H(ISYMAI,ISYMI)
            NT1AM(ISYMAI) = NT1AM(ISYMAI) + NVIR(ISYMA) * NRHF(ISYMI)
            NT1AO(ISYMAI) = NT1AO(ISYMAI) + NBAS(ISYMA) * NRHF(ISYMI)
Chol
            XT1AM(ISYMAI) = XT1AM(ISYMAI) + XVIR(ISYMA) * XRHF(ISYMI)
            XT1AO(ISYMAI) = XT1AO(ISYMAI) + NBAS(ISYMA) * XRHF(ISYMI)
Chol
            NH1AM(ISYMAI) = NH1AM(ISYMAI) + NORB1(ISYMA) * NRHF(ISYMI)
            NG1AM(ISYMAI) = NG1AM(ISYMAI) + NORB2(ISYMA) * NRHF(ISYMI)
celena
            NT1VM(ISYMAI) = NT1VM(ISYMAI) + NVIR(ISYMA) * 
     &                      (NORB1(ISYMI)-NRHFFR(ISYMI))
celena
  200    CONTINUE
  100 CONTINUE
C
      DO 300 ISAIBJ = 1,NSYM
         NT2AM(ISAIBJ)  = 0
         NT2AO(ISAIBJ)  = 0
         NT2AMA(ISAIBJ) = 0
         NT2AMT(ISAIBJ) = 0
         NH2AM(ISAIBJ)  = 0
         NU2AM(ISAIBJ)  = 0
Chol
         XT2AM(ISAIBJ) = 0.0D0
Chol
         DO 400 ISYMBJ = 1,NSYM
            ISYMAI = MULD2H(ISYMBJ,ISAIBJ)
            IF (ISYMBJ .GT. ISYMAI) THEN
               NT2AM(ISAIBJ) = NT2AM(ISAIBJ) +
     &                         NT1AM(ISYMAI) * NT1AM(ISYMBJ)
               NT2AO(ISAIBJ) = NT2AO(ISAIBJ) +
     &                         NT1AO(ISYMAI) * NT1AO(ISYMBJ)
               NT2AMA(ISAIBJ)= NT2AM(ISAIBJ)
               NT2AMT(ISAIBJ)= NT2AM(ISAIBJ) + NT2AMA(ISAIBJ)     
Chol
               XT2AM(ISAIBJ) = XT2AM(ISAIBJ) +
     &                         XT1AM(ISYMAI) * XT1AM(ISYMBJ)
Chol
               NH2AM(ISAIBJ) = NH2AM(ISAIBJ) +
     &                         NH1AM(ISYMAI) * NT1AM(ISYMBJ)
            ELSE IF (ISYMBJ .EQ. ISYMAI) THEN
               NT2AM(ISAIBJ) = NT2AM(ISAIBJ) +
     &                         NT1AM(ISYMAI) * (NT1AM(ISYMBJ)+1)/2
               NT2AO(ISAIBJ) = NT2AO(ISAIBJ) +
     &                         NT1AO(ISYMAI) * (NT1AO(ISYMBJ)+1)/2
               NT2AMA(ISAIBJ)= NT2AM(ISAIBJ)
               NT2AMT(ISAIBJ)= NT2AM(ISAIBJ) + NT2AMA(ISAIBJ)  
Chol
               XT2AM(ISAIBJ) = XT2AM(ISAIBJ) +
     &               XT1AM(ISYMAI) * (XT1AM(ISYMBJ)+1.0D0) / 2.0D0
Chol
               NH2AM(ISAIBJ) = NH2AM(ISAIBJ) +
     &                         NH1AM(ISYMAI) * (NH1AM(ISYMBJ)+1)/2 +
     &                         NH1AM(ISYMAI) * NG1AM(ISYMBJ)
            END IF
C           For [T1+T2,r12] integrals (WK/UniKA/04-11-2002).
            NU2AM(ISAIBJ) = NU2AM(ISAIBJ) +
     &                      NT1AM(ISYMAI) * NT1AM(ISYMBJ)
  400    CONTINUE
  300 CONTINUE
C
      NT1AMX = NT1AM(1)
      NT1AOX = NT1AO(1)
      NH1AMX = NH1AM(1)
      NT2AMX = NT2AM(1)
      NT2AOX = NT2AO(1)
      NU2AMX = NU2AM(1)
      NH2AMX = NH2AM(1)
      NT1VMX = NT1VM(1)
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
      ICOUN7 = 0
      ICOUN8 = 0
      DO 600 ISYMP = 1,NSYM
         ICOUN1 = ICOUN1 + NBAS(ISYMP)*NORB(ISYMP)
         ICOUN2 = ICOUN2 + NBAS(ISYMP)*NRHF(ISYMP)
         ICOUN5 = ICOUN5 + NBAS(ISYMP)*NRHFS(ISYMP)
C
         IRHF(ISYMP) = ICOUN3
         IRHFA(ISYMP) = ICOUN7
         IRHFB(ISYMP) = ICOUN8
         IVIR(ISYMP) = ICOUN4
         ICOUN3 = ICOUN3 + NRHF(ISYMP)
         ICOUN4 = ICOUN4 + NVIR(ISYMP)
         ICOUN7 = ICOUN7 + NRHFA(ISYMP)
         ICOUN8 = ICOUN8 + NRHFB(ISYMP)
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
Chol
         XCOUN4 = 0.0D0
Chol
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
C        For [T1+T2,r12] integrals (WK/UniKA/04-11-2002).
         ICOU20 = 0
         ICOU21 = 0
         ICOU22 = 0
         ICOU23 = 0
C        For R12-index pairs (C. Neiss):
         ICOU24 = 0
         ICOU25 = 0
         ICOU26 = 0
         DO 620 ISYMJ = 1,NSYM
C
            ISYMI  = MULD2H(ISYMJ,ISYMK)
C
            IT1AM(ISYMI,ISYMJ)  = ICOUN1
            IH1AM(ISYMI,ISYMJ)  = ICOU21
            IG1AM(ISYMI,ISYMJ)  = ICOU23
            IT1AO(ISYMI,ISYMJ)  = ICOUN5
            IT1AMT(ISYMI,ISYMJ) = ICOU11
            IT1AOT(ISYMI,ISYMJ) = ICOU12
            IEMAT1(ISYMI,ISYMJ) = ICOU15
            IMATAV(ISYMI,ISYMJ) = ICOU18
C
            IF (ISYMJ .GE. ISYMI) THEN
C              For [T1+T2,r12] integrals (WK/UniKA/04-11-2002).
               IU2AM(ISYMI,ISYMJ) = ICOU20
               IU2AM(ISYMJ,ISYMI) = ICOU20
               IH2AM(ISYMI,ISYMJ) = ICOU22
               IH2AM(ISYMJ,ISYMI) = ICOU22
               ICOU20 = ICOU20 + NT1AM(ISYMI)*NT1AM(ISYMJ)
               IF (ISYMJ .EQ. ISYMI) THEN
                  ICOU22 = ICOU22 + NH1AM(ISYMI)*(NH1AM(ISYMJ)+1)/2 +
     &                              NH1AM(ISYMI)*NG1AM(ISYMJ)
               ELSE
                  ICOU22 = ICOU22 + NH1AM(ISYMI)*NT1AM(ISYMJ)
               END IF
            END IF
C
            ICOUN1 = ICOUN1 + NRHF(ISYMJ)*NVIR(ISYMI)
            ICOU21 = ICOU21 + NRHF(ISYMJ)*NORB1(ISYMI)
            ICOU23 = ICOU23 + NRHF(ISYMJ)*NORB2(ISYMI)
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
Chol
            XCOUN4 = XCOUN4 + XRHF(ISYMI)*XRHF(ISYMJ)
Chol
            ICOU24 = ICOU24 + NRHFB(ISYMI)*NRHFB(ISYMJ)
            ICOU25 = ICOU25 + NRHFB(ISYMI)*NRHFA(ISYMJ)
            ICOU26 = ICOU26 + NVIR(ISYMI)*NRHFB(ISYMJ)
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
Chol
         XMATIJ(ISYMK) = XCOUN4
Chol
         NMATKL(ISYMK) = ICOU24
         NMATKI(ISYMK) = ICOU25
         NMATAK(ISYMK) = ICOU26
C
  610 CONTINUE
C
      DO 630 ISYMK = 1,NSYM
         ICOUN1 = 0
         ICOUN2 = 0
         ICOUN3 = 0
C        For R12 (C. Neiss):
         ICOUN4 = 0
         ICOUN5 = 0
         ICOUN6 = 0
         ICOUN7 = 0
         ICOUN8 = 0
         DO 640 ISYMJ = 1,NSYM
            ISYMI = MULD2H(ISYMJ,ISYMK)
C
            IF (ISYMJ .GT. ISYMI) THEN
               ICOUN1 = ICOUN1 + NMATIJ(ISYMI)*NMATIJ(ISYMJ)
               ICOUN4 = ICOUN4 + NMATKI(ISYMI)*NMATKI(ISYMJ)
               ICOUN5 = ICOUN5 + NMATKL(ISYMI)*NMATKL(ISYMJ)
               IT2R12(ISYMI,ISYMJ) = ICOUN8
               IT2R12(ISYMJ,ISYMI) = ICOUN8
               ICOUN8 = ICOUN8 + NMATAK(ISYMI)* NMATAK(ISYMJ)
            ELSE IF (ISYMK .EQ. 1) THEN
               ICOUN1 = ICOUN1 + NMATIJ(ISYMI)*(NMATIJ(ISYMI)+1)/2
               ICOUN4 = ICOUN4 + NMATKI(ISYMI)*(NMATKI(ISYMI)+1)/2
               ICOUN5 = ICOUN5 + NMATKL(ISYMI)*(NMATKL(ISYMI)+1)/2
               IT2R12(ISYMI,ISYMJ) = ICOUN8
               ICOUN8 = ICOUN8 + NMATAK(ISYMI)*(NMATAK(ISYMJ)+1)/2
            ENDIF
C
            ICOUN2 = ICOUN2 + NVIR(ISYMI)*NVIR(ISYMJ)
            ICOUN3 = ICOUN3 + NMATIJ(ISYMI)*NMATIJ(ISYMJ)
            ICOUN6 = ICOUN6 + NMATIJ(ISYMI)*NMATKL(ISYMJ)
            ICOUN7 = ICOUN7 + NMATKL(ISYMI)*NMATKL(ISYMJ)
C
  640    CONTINUE
C
         NGAMMA(ISYMK) = ICOUN1
         NMATAB(ISYMK) = ICOUN2
         NGAMSQ(ISYMK) = ICOUN3
C        For R12 (C. Neiss):
         NTR12AM(ISYMK)  = ICOUN4
         NR12R12P(ISYMK) = ICOUN5
         NTR12SQ(ISYMK)  = ICOUN6
         NR12R12SQ(ISYMK)= ICOUN7
         NT2R12(ISYMK)  = ICOUN8
C
  630 CONTINUE
C
      IF ((.NOT. ONEAUX) .AND. (.NOT.LABEL.EQ.'TRCCINT ')) THEN
c     IF (.NOT. ONEAUX) THEN
         NH1AMX = NT1AMX
         NH2AMX = NT2AMX
         DO ISYMI = 1,NSYM
            NH1AM(ISYMI) = NT1AM(ISYMI)
            NH2AM(ISYMI) = NT2AM(ISYMI)
            DO ISYMJ = 1,NSYM
               IH1AM(ISYMI,ISYMJ) = IT1AM(ISYMI,ISYMJ)
               IH2AM(ISYMI,ISYMJ) = IT2AM(ISYMI,ISYMJ)
            ENDDO
         ENDDO
      END IF
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
         ICOUN7 = 0
         ICOUN8 = 0
         ICOUN9 = 0
         ICOUN10 = 0
         ICOUN11 = 0
         DO 560 ISYMK = 1,NSYM
            ISYMIJ = MULD2H(ISYMK,ISYIJK)
            IMAIJK(ISYMIJ,ISYMK) = ICOUN1
            IT2AIJ(ISYMIJ,ISYMK) = ICOUN2
            IMAIJA(ISYMIJ,ISYMK) = ICOUN3
            ID2IJG(ISYMIJ,ISYMK) = ICOUN4
            ID2AIG(ISYMIJ,ISYMK) = ICOUN5
            ID2ABG(ISYMIJ,ISYMK) = ICOUN6
            IMAABC(ISYMIJ,ISYMK) = ICOUN7
            IMAABI(ISYMIJ,ISYMK) = ICOUN8
            IMAIAB(ISYMIJ,ISYMK) = ICOUN9
            IMAIAJ(ISYMIJ,ISYMK) = ICOUN10
Cholesky
            IT2VO(ISYMIJ,ISYMK)  = ICOUN11
Cholesky
            ICOUN1 = ICOUN1 + NMATIJ(ISYMIJ)*NRHF(ISYMK)
            ICOUN2 = ICOUN2 + NVIR(ISYMIJ)*NMATIJ(ISYMK)
            ICOUN3 = ICOUN3 + NMATIJ(ISYMIJ)*NVIR(ISYMK)
            ICOUN4 = ICOUN4 + NMATIJ(ISYMIJ)*NBAS(ISYMK)
            ICOUN5 = ICOUN5 + NT1AM(ISYMIJ)*NBAS(ISYMK)
            ICOUN6 = ICOUN6 + NMATAB(ISYMIJ)*NBAS(ISYMK)
            ICOUN7 = ICOUN7 + NMATAB(ISYMIJ)*NVIR(ISYMK)
            ICOUN8 = ICOUN8 + NMATAB(ISYMIJ)*NRHF(ISYMK)
            ICOUN9 = ICOUN9 + NT1AM(ISYMIJ)*NVIR(ISYMK)
            ICOUN10 = ICOUN10 + NRHF(ISYMIJ)*NT1AM(ISYMK)
Cholesky
            ICOUN11 = ICOUN11 + NMATAB(ISYMIJ)*NMATIJ(ISYMK)
Cholesky
  560    CONTINUE
         NMAIJK(ISYIJK) = ICOUN1
         NT2AIJ(ISYIJK) = ICOUN2
         NMAIJA(ISYIJK) = ICOUN3
         ND2IJG(ISYIJK) = ICOUN4
         ND2AIG(ISYIJK) = ICOUN5
         ND2ABG(ISYIJK) = ICOUN6
         NMAABC(ISYIJK) = ICOUN7
         NMAABI(ISYIJK) = ICOUN8
         NMAIAB(ISYIJK) = ICOUN9
         NMAIAJ(ISYIJK) = ICOUN10
  550 CONTINUE
C
      DO 570 ISIJKD = 1,NSYM
         ICOUN1 = 0
         ICOUN2 = 0
         ICOUN3 = 0
         ICOUN4 = 0
         ICOUN5 = 0
         ICOUN6 = 0
         ICOUN7 = 0
         ICOUN8 = 0
         ICOUN9 = 0!added by FP 16-03-04, needed for new CC3 LHTR
         ICOUN10 = 0!added by FP 16-03-04, needed for new CC3 LHTR
         ICOUN11 = 0!added by FP 16-03-04, needed for new CC3 LHTR

C
         DO 580 ISYMD = 1,NSYM
            ISYIJK = MULD2H(ISYMD,ISIJKD)
            I3ODEL(ISYIJK,ISYMD) = ICOUN1
            I3ORHF(ISYIJK,ISYMD) = ICOUN2
            I3OVIR(ISYIJK,ISYMD) = ICOUN3
            I3VDEL(ISYIJK,ISYMD) = ICOUN4
            I3VVIR(ISYIJK,ISYMD) = ICOUN5
            I3VOOO(ISYIJK,ISYMD) = ICOUN6
            IMAABCI(ISYIJK,ISYMD) = ICOUN7
            IMAAB_CI(ISYIJK,ISYMD) = ICOUN8
            I3AORHF(ISYIJK,ISYMD) = ICOUN9!added by FP 16-03-04 (new CC3 LHTR)
            I3AO(ISYIJK,ISYMD) = ICOUN10!added by FP 16-03-04 (new CC3 LHTR)
            IRHF3O(ISYIJK,ISYMD) = ICOUN11!added by FP 29-03-04 (new CC3 LHTR)

            ICOUN1 = ICOUN1 + NMAIJK(ISYIJK)*NBAS(ISYMD)
            ICOUN2 = ICOUN2 + NMAIJK(ISYIJK)*NRHF(ISYMD)
            ICOUN3 = ICOUN3 + NMAIJK(ISYIJK)*NVIR(ISYMD)
            ICOUN4 = ICOUN4 + NMAABC(ISYIJK)*NBAS(ISYMD)
            ICOUN5 = ICOUN5 + NMAABC(ISYIJK)*NVIR(ISYMD)
            ICOUN6 = ICOUN6 + NVIR(ISYIJK)*NMAIJK(ISYMD)
            ICOUN7 = ICOUN7 + NMAABC(ISYIJK)*NRHF(ISYMD)
            ICOUN8 = ICOUN8 + NMATAB(ISYIJK)*NT1AM(ISYMD)
            ICOUN9 = ICOUN9 + NDISAOSQ(ISYIJK)*NRHF(ISYMD)!FP 16-03-04(CC3 LHTR)
            ICOUN10 = ICOUN10 + N2BST(ISYIJK)*NBAS(ISYMD)!FP (CC3 LHTR)
            ICOUN11 = ICOUN11 + NRHF(ISYIJK)*NMAIJK(ISYMD)!FP (CC3 LHTR)
C
  580    CONTINUE
         N3ODEL(ISIJKD) = ICOUN1
         N3ORHF(ISIJKD) = ICOUN2
         N3OVIR(ISIJKD) = ICOUN3
         N3VDEL(ISIJKD) = ICOUN4
         N3VVIR(ISIJKD) = ICOUN5
         N3VOOO(ISIJKD) = ICOUN6
         NMAABCI(ISIJKD) = ICOUN7
         NMAAB_CI(ISIJKD) = ICOUN8
         N3AORHF(ISIJKD) = ICOUN9!FP 16-03-04(CC3 LHTR)
         N3AO(ISIJKD) = ICOUN10!FP 16-03-04(CC3 LHTR)
         NRHF3O(ISIJKD) = ICOUN11!FP 29-03-04(CC3 LHTR)

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
         ICOU10 = 0
C        For R12 (C. Neiss):
         ICOU11 = 0
         ICOU12 = 0
         ICOU13 = 0
         ICOU14 = 0
         ICOU15 = 0
         ICOU16 = 0
         ICOU17 = 0
C
         DO 710 ISYMJ = 1,NSYM
C
            ISYMK = MULD2H(ISYMJ,ISYMI)
C
            IDSAOG(ISYMJ,ISYMI)   = ICOUN3
            IMATIJ(ISYMK,ISYMJ)   = ICOUN4
            IGAMMA(ISYMK,ISYMJ)   = ICOUN5
            IGAMMA(ISYMJ,ISYMK)   = ICOUN5
            IGAMSQ(ISYMJ,ISYMK)   = ICOU10
            IMATAB(ISYMK,ISYMJ)   = ICOUN6
            IDSAOGSQ(ISYMJ,ISYMI) = ICOUN9
            IMATKL(ISYMK,ISYMJ)   = ICOU11
            IMATKI(ISYMK,ISYMJ)   = ICOU12
            ITR12AM(ISYMK,ISYMJ)  = ICOU13
            ITR12AM(ISYMJ,ISYMK)  = ICOU13
            ITR12SQ(ISYMJ,ISYMK)  = ICOU14
            ITR12SQT(ISYMJ,ISYMK) = ICOU17
C            ITR12SQT(ISYMJ,ISYMK) = ITR12SQ(ISYMJ,ISYMK)
            IR12R12P(ISYMK,ISYMJ) = ICOU15
            IR12R12P(ISYMJ,ISYMK) = ICOU15
            IR12R12SQ(ISYMJ,ISYMK)= ICOU16
C
            ICOUN3 = ICOUN3 + NNBST(ISYMK)*NBAS(ISYMJ)
            ICOUN4 = ICOUN4 + NRHF(ISYMK)*NRHF(ISYMJ)
            ICOUN6 = ICOUN6 + NVIR(ISYMK)*NVIR(ISYMJ)
            ICOUN9 = ICOUN9 + N2BST(ISYMK)*NBAS(ISYMJ)
            ICOU10 = ICOU10 + NMATIJ(ISYMK)*NMATIJ(ISYMJ)
            ICOU11 = ICOU11 + NRHFB(ISYMK)*NRHFB(ISYMJ)
            ICOU12 = ICOU12 + NRHFB(ISYMK)*NRHFA(ISYMJ)
            ICOU14 = ICOU14 + NMATIJ(ISYMK)*NMATKL(ISYMJ)
            ICOU17 = ICOU17 + NMATKL(ISYMK)*NMATIJ(ISYMJ)
            ICOU16 = ICOU16 + NMATKL(ISYMK)*NMATKL(ISYMJ)
C
            IF (ISYMJ .GT. ISYMK) THEN
               ICOUN5 = ICOUN5 + NMATIJ(ISYMK)*NMATIJ(ISYMJ)
               ICOU13 = ICOU13 + NMATKI(ISYMK)*NMATKI(ISYMJ)
               ICOU15 = ICOU15 + NMATKL(ISYMK)*NMATKL(ISYMJ)
            ELSE IF (ISYMI .EQ. 1) THEN
               ICOUN5 = ICOUN5 + NMATIJ(ISYMJ)*(NMATIJ(ISYMJ)+1)/2
               ICOU13 = ICOU13 + NMATKI(ISYMJ)*(NMATKI(ISYMJ)+1)/2
               ICOU15 = ICOU15 + NMATKL(ISYMJ)*(NMATKL(ISYMJ)+1)/2
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
         ICOUN2 = 0
C           
         XCOUN1 = 0.0D0
C
         DO 910 ISYMJ = 1,NSYM
C
            ISYMI  = MULD2H(ISYMJ,ISYM)
            IT2AOS(ISYMI,ISYMJ) = ICOUNT
            ITG2SQ(ISYMI,ISYMJ) = ICOUN2
C
            ICOUNT = ICOUNT + NT1AO(ISYMI)*NT1AO(ISYMJ)
            ICOUN1 = ICOUN1 + NT1AM(ISYMI)*NT1AM(ISYMJ)
            ICOUN2 = ICOUN2 + NT1AM(ISYMI)*NG1AM(ISYMJ)
C           
            XCOUN1 = XCOUN1 + XT1AM(ISYMI)*XT1AM(ISYMJ)
C
  910    CONTINUE
C
         NT2AOS(ISYM) = ICOUNT
         NT2SQ(ISYM)  = ICOUN1
         NTG2SQ(ISYM) = ICOUN2
C
         XT2SQ(ISYM)  = XCOUN1
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
Chol
         XCOUN2 = 0.0D0
Chol
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
            XCOUN2 = XCOUN2 + XT1AM(ISYMCK)*XRHF(ISYMA)
C
 1050    CONTINUE
C
         NCKA(ISYCKA)   = ICOUN1
         NCKI(ISYCKA)   = ICOUN2
         NCKATR(ISYCKA) = ICOUN3
         NCKASR(ISYCKA) = ICOUN4
C
         XCKI(ISYCKA)   = XCOUN2
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
         ICOUN7 = 0
         ICOUN8 = 0
         ICOUN9 = 0
         ICOUN10 = 0
         ICOUN11 = 0
C
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
            ICDKAO(ISYAIK,ISYMJ) = ICOUN7
            ICDKVI(ISYAIK,ISYMJ) = ICOUN8
            IMAJBAI(ISYAIK,ISYMJ)  = ICOUN9
            IMAAOBCI(ISYAIK,ISYMJ)  = ICOUN10
            IMAJBAIT(ISYAIK,ISYMJ)  = ICOUN11
C
            ICOUN1 = ICOUN1 + NCKI(ISYAIK)*NRHF(ISYMJ)
            ICOUN2 = ICOUN2 + NCKI(ISYAIK)*NBAS(ISYMJ)
            ICOUN3 = ICOUN3 + NCKA(ISYAIK)*NVIR(ISYMJ)
            ICOUN4 = ICOUN4 + NCKATR(ISYAIK)*NBAS(ISYMJ)
            ICOUN5 = ICOUN5 + NCKATR(ISYAIK)*NVIR(ISYMJ)
            ICOUN6 = ICOUN6 + NCKI(ISYAIK)*NVIR(ISYMJ)
            ICOUN7 = ICOUN7 + NCKASR(ISYAIK)*NBAS(ISYMJ)
            ICOUN8 = ICOUN8 + NCKASR(ISYAIK)*NVIR(ISYMJ)
            ICOUN9 = ICOUN9 + NRHF(ISYAIK)*NCKATR(ISYMJ)
            ICOUN10 = ICOUN10 + NVIR(ISYAIK)*NCKATR(ISYMJ)
            ICOUN11 = ICOUN11 + NCKATR(ISYAIK)*NRHF(ISYMJ)
C
            IF (CCSDT.OR.CCPT.OR.CCP3.OR.CCRT.OR.
     *                CHOPT .OR. CCR3.OR.CCR1A.OR.CCR1B) THEN
               IF (ICOUN1 .LT. 0) WRITE(LUPRI,*)
     &                             'Negative ICKITR in CCSD_INIT1'
               IF (ICOUN2 .LT. 0) WRITE(LUPRI,*)
     &                             'Negative ICKID in CCSD_INIT1'
               IF (ICOUN6 .LT. 0) WRITE(LUPRI,*)
     &                             'Negative IT2SP in CCSD_INIT1'
               IF (ICOUN9 .LT. 0) WRITE(LUPRI,*)
     &                             'Negative ICKDAO in CCSD_INIT1'
               IF (ICOUN11 .LT. 0) WRITE(LUPRI,*)
     &                             'Negative IMAJBAIT in CCSD_INIT1'
               IF ((ICOUN1 .LT. 0) .OR. (ICOUN2 .LT. 0) .OR.
     &             (ICOUN6 .LT. 0) .OR. (ICOUN9 .LT. 0)) THEN
                   WRITE(LUPRI,'(A,A)') 
     &                  'Calculation too large for 32-bit integers',
     &                  'Try rebuilding Dalton using 64-bit integers'
                   STOP 'Negative index in CCSD_INIT1'
               END IF
            END IF
C
 1080    CONTINUE
C
         NCKIJ(JSAIKJ) = ICOUN1
         NMAAOBCI(JSAIKJ) = ICOUN10
C
 1070 CONTINUE
C
C     FIND MAX LENGTH OF NCKIJ(JSAIKJ)
C
      NCKIJMAX = 0
      DO I = 1,NSYM
         NCKIJMAX = MAX(NCKIJMAX,NCKIJ(I))
      ENDDO
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
C------------------------------------------------------------
C     set offsets and dimensions for CCR12
C------------------------------------------------------------
      DO ISYMAK = 1, NSYM
         NVAJKL(ISYMAK) = 0
         NVABKL(ISYMAK) = 0
         ICOUNT1 = 0
         ICOUNT2 = 0
         DO ISYMK = 1, NSYM
            ISYMA = MULD2H(ISYMAK,ISYMK)
            IVAJKL(ISYMA,ISYMK) = ICOUNT1
            IVABKL(ISYMA,ISYMK) = ICOUNT2
            NVAJKL(ISYMAK) = NVAJKL(ISYMAK) + NT1AO(ISYMA)*NMATKL(ISYMK)
            NVABKL(ISYMAK) = NVABKL(ISYMAK) + N2BST(ISYMA)*NMATKL(ISYMK)
            ICOUNT1 = ICOUNT1 + NT1AO(ISYMA)*NMATKL(ISYMK)
            ICOUNT2 = ICOUNT2 + N2BST(ISYMA)*NMATKL(ISYMK)
         END DO
      END DO

C------------------------------------------------------------
C     set offset arrays ISWTL and ISTLN and dimensions NIMFN:
C------------------------------------------------------------
C
      DO ISYMDL = 1, NSYM
        IOFF = 0
        DO ISYML = 1, NSYM
          ISWMAT = MULD2H(ISYMDL,ISYML)
          ISWTL(ISWMAT,ISYML) = IOFF
          IOFF = IOFF + NT2SQ(ISWMAT)*NRHF(ISYML)
        END DO
      END DO

      DO ISAIBJ = 1, NSYM
        IOFF = 0
        DO ISYMJ = 1, NSYM
          ISAIB = MULD2H(ISAIBJ,ISYMJ)
          ISTLN(ISAIB,ISYMJ) = IOFF
          IOFF = IOFF + NCKATR(ISAIB)*NRHF(ISYMJ)
        END DO
      END DO

      DO ISYM = 1, NSYM
        ILEN = 0
        DO ISYMFN = 1, NSYM
          ISYMIM = MULD2H(ISYM,ISYMFN)
          ILEN   = ILEN + NMATIJ(ISYMIM)*NT1AM(ISYMFN)
        END DO
        NIMFN(ISYM) = ILEN
      END DO
C
      IF (IPRINT .GT. 9) THEN
         CALL AROUND('Information from CCSDSYM')
         WRITE(LUPRI,1) 'NRHF   :',(NRHF(I),   I=1,NSYM)
         WRITE(LUPRI,1) 'NRHFS  :',(NRHFS(I),  I=1,NSYM)
         WRITE(LUPRI,1) 'NRHFA  :',(NRHFA(I),  I=1,NSYM)
         WRITE(LUPRI,1) 'NRHFSA :',(NRHFSA(I), I=1,NSYM)
         WRITE(LUPRI,1) 'NRHFB  :',(NRHFB(I),  I=1,NSYM)
         WRITE(LUPRI,1) 'NORBS  :',(NORBS(I),  I=1,NSYM)
         WRITE(LUPRI,1) 'NNBST  :',(NNBST(I),  I=1,NSYM)
         WRITE(LUPRI,1) 'NT1AM  :',(NT1AM(I),  I=1,NSYM)
         WRITE(LUPRI,1) 'NT2AM  :',(NT2AM(I),  I=1,NSYM)
         WRITE(LUPRI,1) 'NG1AM  :',(NG1AM(I),  I=1,NSYM)
         WRITE(LUPRI,1) 'NH1AM  :',(NH1AM(I),  I=1,NSYM)
         WRITE(LUPRI,1) 'NH2AM  :',(NH2AM(I),  I=1,NSYM)
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
         WRITE(LUPRI,1) 'NMATKI :',(NMATKI(I), I=1,NSYM)
         WRITE(LUPRI,1) 'NMATKL :',(NMATKL(I), I=1,NSYM)
         WRITE(LUPRI,1) 'NGAMMA :',(NGAMMA(I), I=1,NSYM)
         WRITE(LUPRI,1) 'NTR12AM:',(NTR12AM(I),I=1,NSYM)
         WRITE(LUPRI,1) 'NGAMSQ :',(NGAMSQ(I), I=1,NSYM)
         WRITE(LUPRI,1) 'NTR12SQ:',(NTR12SQ(I),I=1,NSYM)
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
         DO I = 1,NSYM
            WRITE(LUPRI,1) 'IH1AM  :',(IH1AM(I,J), J=1,NSYM)
         END DO
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
         DO I = 1,NSYM
            WRITE(LUPRI,1) 'IMATKI :',(IMATKI(I,J), J=1,NSYM)
         END DO
         WRITE(LUPRI,*)
         DO I = 1,NSYM
            WRITE(LUPRI,1) 'IMATKL :',(IMATKL(I,J), J=1,NSYM)
         END DO
         WRITE(LUPRI,*)
         DO 9911 I = 1,NSYM
            WRITE(LUPRI,1) 'IGAMMA :',(IGAMMA(I,J), J=1,NSYM)
 9911    CONTINUE
         WRITE(LUPRI,*)
         DO I = 1, NSYM
            WRITE(LUPRI,1) 'ITR12AM:',(ITR12AM(I,J),J=1,NSYM)
         END DO
         WRITE(LUPRI,*)
         DO I = 1, NSYM
            WRITE(LUPRI,1) 'IGAMSQ :',(IGAMSQ(I,J), J=1,NSYM)
         END DO
         WRITE(LUPRI,*)
         DO I = 1, NSYM
            WRITE(LUPRI,1) 'ITR12SQ:',(ITR12SQ(I,J),J=1,NSYM)
         END DO
         WRITE(LUPRI,*)
         DO I = 1, NSYM
            WRITE(LUPRI,1) 'ITR12SQT:',(ITR12SQT(I,J),J=1,NSYM)
         END DO
         WRITE(LUPRI,*)
         DO I = 1, NSYM
            WRITE(LUPRI,1) 'IR12R12SQ:',(IR12R12SQ(I,J),J=1,NSYM)
         END DO
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
C
         CALL AROUND('Information from DCCSDSYM')
         WRITE(LUPRI,2) 'XT1AM  :',(XT1AM(I), I=1,NSYM)
         WRITE(LUPRI,2) 'XT2AM  :',(XT2AM(I), I=1,NSYM)
         WRITE(LUPRI,2) 'XT2SQ  :',(XT2SQ(I), I=1,NSYM)
         WRITE(LUPRI,2) 'XCKI   :',(XCKI(I),  I=1,NSYM)
         WRITE(LUPRI,2) 'XMATIJ :',(XMATIJ(I), I=1,NSYM)
         WRITE(LUPRI,2) 'XT1AO  :',(XT1AO(I),  I=1,NSYM)
C
      END IF
C
      CALL QEXIT('CCSD_INIT1')
C
      RETURN
C
    1 FORMAT(3X,A8,8I8)
    2 FORMAT(3X,A8,8D8.3)
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
#include <priunit.h>
#include <ccorb.h>
      DIMENSION FOCK(NORBT),WORK(LWORK)
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      CALL QENTER('FOCK_REORDER')
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
            WRITE(55,'(4e30.20)') FOCK(I)
  200    CONTINUE
      END IF
C
      CALL DCOPY(NORBT,WORK,1,FOCK,1)
C
      CALL QEXIT('FOCK_REORDER')
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
      LOGICAL FRORHF, FROVIR
C
      CALL QENTER('CMO_REORDER')
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
      IF (IPRINT .GT. 200) THEN
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
      CALL QEXIT('CMO_REORDER')
C
      RETURN
C
    1 FORMAT(//,7X,'Symmetry number :',I5)
    2 FORMAT(//,7X,'Sirius ordering')
    3 FORMAT(7X,'---------------')
    4 FORMAT(//,7X,'CCSD ordering occupied part')
    5 FORMAT(7X,'---------------------------')
    6 FORMAT(//,7X,'CCSD ordering virtual part')
    7 FORMAT(7X,'--------------------------')
    8 FORMAT(//,7X,'This symmetry is empty')
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
      CALL QENTER('CCSD_SYMSQO')
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
      CALL QEXIT('CCSD_SYMSQO')
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
      CALL QENTER('CCSD_SYMSQ')
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
#else

C     SYS_CRAY code

      KOFF = I2BST(ISYMAB) + 1
      CALL GATHER(N2BST(ISYMAB),SCR,DISTAB,IADRPK(KOFF))
#endif
C
      CALL QEXIT('CCSD_SYMSQ')
C
      RETURN
      END
      SUBROUTINE CCSD_SYMSQT(DISTAB1,ISYMAB,SCR)
C
#include <implicit.h>
      DIMENSION DISTAB1(*), SCR(*)
#include <priunit.h>
#include <maxorb.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <symsq.h>
C
      CALL QENTER('CCSD_SYMSQT')
C
#if defined (SYS_CRAY)
       WRITE(LUPRI,*) 'computation of Bijal on SYS_CRAY not implemented'
       CALL QUIT('computation of Bijal on SYS_CRAY not implemented')
#endif
C
       DO IJSQ = 1,N2BST(ISYMAB)
         KOFF = I2BST(ISYMAB) + IJSQ
         IJPK = IADRPK(KOFF)
         SCR(IJSQ) = -DISTAB1(IJPK) 
       END DO  
C
C
      CALL QEXIT('CCSD_SYMSQT')
C
      RETURN
      END
C  /* Deck cc3_sort1 */
      SUBROUTINE CC3_SORT1(WORK,LWORK,IOPT,ISYINT,LU3SRT,FN3SRT,
     *                     LU3VI,FN3VI,LU3VI2,FN3VI2,
     *                     LU3FOP,FN3FOP,LU3FOP2,FN3FOP2)
C
C     Henrik Koch and Alfredo Sanchez.       28-May-1995
C
C     Kasper Hald fall 2001 - Added 2*C-E for ccsd(t) f.o.p.
C
C     FN3VI can be FN3VI or FNDELD
C
C     Sort virtual integrals for perturbative triples.
C
#include <implicit.h>
      DIMENSION WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccinftap.h>
#include <ccsdsym.h>
#include <ccfop.h>
#include <ccsdinp.h>
#include <ccsections.h>
C
      PARAMETER (TWO = 2.0D0)
C
      CHARACTER*(*) FN3SRT, FN3VI, FN3VI2, FN3FOP, FN3FOP2
C
      CALL QENTER('CC3_SORT1')
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
Casm     Apparently, it is not possible to read more than 2 Gb (268435456 dw)
C        
         MXDALF = 268435455 / NCKATR(ISYCKB)
         NDISTR = MIN(NDISTR,MXDALF)
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
                     CALL PUTWA2(LU3VI,FN3VI,WORK(KSCR2),IOFF,LENGTH)
C
                  ENDIF
C
  130          CONTINUE
  120       CONTINUE
C
C----------------------------------------------------------------------
C           Sort integrals (ck,del,b) from (ck,b,del).
C           for (ccpt and ccfop) = true  and iopt = 1 
C           and construct 2*C-E.
C----------------------------------------------------------------------
C
            IF ((IOPT .EQ. 1) .AND.
     &          (CC3 .OR. (CCPT .AND. (CCFOP.OR.ETACCPT)))) THEN
C
               DO ISYMB = 1,NSYM
C
                  ISYMCK = MULD2H(ISYCKB,ISYMB)
                  ISYCKD = MULD2H(ISYMCK,ISYMD)
C
                  DO B = 1,NVIR(ISYMB)
C
                     DO I = 1,NUMD
C
                        ID = ID1 + I - 1
C
                        DO ISYMK = 1, NSYM
C
                           ISYMC  = MULD2H(ISYMCK,ISYMK)
                           ISYMBK = MULD2H(ISYMB,ISYMK)
C
                           DO K = 1, NRHF(ISYMK)
                           DO C = 1, NVIR(ISYMC)
C
                              KOFF1 = KSCR1
     *                              + NCKATR(ISYCKB)*(I - 1)
     *                              + ICKATR(ISYMCK,ISYMB)
     *                              + NT1AM(ISYMCK)*(B - 1)
     *                              + IT1AM(ISYMC,ISYMK)
     *                              + NVIR(ISYMC)*(K-1) + C - 1
                              KOFF2 = KSCR1
     *                              + NCKATR(ISYCKB)*(I - 1)
     *                              + ICKATR(ISYMBK,ISYMC)
     *                              + NT1AM(ISYMBK)*(C - 1)
     *                              + IT1AM(ISYMB,ISYMK)
     *                              + NVIR(ISYMB)*(K-1) + B - 1
                              KOFF3 = KSCR2
     *                              + NT1AM(ISYMCK)*(I - 1)
     *                              + IT1AM(ISYMC,ISYMK)
     *                              + NVIR(ISYMC)*(K-1) + C - 1
C
                              WORK(KOFF3) = TWO*WORK(KOFF1) 
     *                                    - WORK(KOFF2)
C
                           ENDDO  ! B
                           ENDDO  ! K
                        ENDDO     ! ISYMK
                     ENDDO        ! I
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
     *                       + NCKA(ISYCKD)*(B - 1)
     *                       + ICKA(ISYMCK,ISYMD)
     *                       + NT1AM(ISYMCK)*(ID1 - 1) + 1
C
                        CALL PUTWA2(LU3FOP,FN3FOP,WORK(KSCR2),
     *                              IOFF,LENGTH)
C
                     ENDIF
C
                  ENDDO   ! B
               ENDDO      ! ISYMB
C
            ENDIF
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
                        CALL PUTWA2(LU3VI2,FN3VI2,WORK(KSCR2),IOFF,
     *                              LENGTH)
                  ENDIF
C
  160          CONTINUE
  150       CONTINUE
C
C------------------------------------------------------------------------
C           For iopt = 1 and (ccpt.and.ccfop) = .true. construct 2*C - E
C           Sort integrals 2*(bk,c,del) - (ck,b,del)
C           from (ck,b,del).
C------------------------------------------------------------------------
C
            IF ((IOPT .EQ. 1 .AND. 
     *           (CC3 .OR. (CCPT .AND. (CCFOP.OR.ETACCPT))) )
     *           .OR. (IOPT .EQ. 3)) THEN

C
               DO ISYMC = 1,NSYM
C
               ISYMBK = MULD2H(ISYCKB,ISYMC)
               ISYBKD = MULD2H(ISYMBK,ISYMD)
C
               DO C = 1,NVIR(ISYMC)
                  DO I = 1,NUMD
C
                     ID = ID1 + I - 1
C
                     DO ISYMK = 1,NSYM
C
                        ISYMB  = MULD2H(ISYMBK,ISYMK)
                        ISYMCK = MULD2H(ISYMC,ISYMK)
C
                        NTOTCK = MAX(NT1AM(ISYMCK),1)
C
                        DO K = 1,NRHF(ISYMK)
                        DO B = 1,NVIR(ISYMB)

C
                           KOFF1 = KSCR1
     *                           + NCKATR(ISYCKB)*(I - 1)
     *                           + ICKATR(ISYMCK,ISYMB)
     *                           + NT1AM(ISYMCK)*(B-1)
     *                           + IT1AM(ISYMC,ISYMK)
     *                           + NVIR(ISYMC)*(K - 1) + C - 1
                           KOFF2 = KSCR1
     *                           + NCKATR(ISYCKB)*(I - 1)
     *                           + ICKATR(ISYMBK,ISYMC)
     *                           + NT1AM(ISYMBK)*(C-1)
     *                           + IT1AM(ISYMB,ISYMK)
     *                           + NVIR(ISYMB)*(K - 1) + B - 1
C
                           KOFF3 = KSCR2
     *                           + NT1AM(ISYMBK)*(I - 1)
     *                           + IT1AM(ISYMB,ISYMK)
     *                           + NVIR(ISYMB)*(K - 1) + B - 1
C
                           WORK(KOFF3) = TWO*WORK(KOFF1)-WORK(KOFF2)
C
                        ENDDO ! B
                        ENDDO ! K
                     ENDDO    ! ISYMK
                  ENDDO       ! I
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
                     CALL PUTWA2(LU3FOP2,FN3FOP2,WORK(KSCR2),
     *                           IOFF,LENGTH)
                  ENDIF
C
C
               ENDDO          ! C
               ENDDO          ! ISYMC
            ENDIF
C
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CC3_SORT1')
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
      LOGICAL FRORHF, FROVIR
C
      CALL QENTER('CCSD_DELFRO')
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
      CALL QEXIT('CCSD_DELFRO')
C
      RETURN
      END
C  /* Deck CC_freezer*/
      SUBROUTINE CC_FREEZER(FOCDIA,NF,NFCS,NFVS,WORK,LWORK,LABEL)
C
C     Ove Christiansen 230899, Find and freeze NFC lowest-lying/NFV highest lying
C     canonical orbitals in CC calculation.
C
#include <implicit.h>
#include <maxorb.h>
C
      DIMENSION FOCDIA(NF),WORK(LWORK),NFCS(8),NFVS(8),IPLACE(MXCORB)
      CHARACTER*8 LABEL
#include <ccorb.h>
#include <priunit.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
C
      CALL QENTER('CC_FREEZER')
C
      IF (LABEL .EQ. 'FULLBAS ') THEN
         NFC0 = NFC*2
      ELSE
         NFC0 = NFC
      END IF
C
      IF (IPRINT.GT.5) WRITE(LUPRI,*) ' In CC_FREEZER: ' 
      IF (IPRINT.GT.5) WRITE(LUPRI,*) 
     *' Freezing occupied, virtual:',NFC,NFV
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
       MXELMN = NFC0
       NELMN  = NFC0
       THRDIA = 1.0D-06
       CALL FNDMN3(FOCDIA,NF,MXELMN,IPLACE,
     *            NELMN,IPRINT,THRDIA)
C
C-------------------------------------------------------------
C        Find # frozen orbitals in each symmetryclass: NRHFFR 
C-------------------------------------------------------------
C
       DO I = 1,NFC0
         IF (LABEL .NE. 'FULLBAS ' .OR.
     *      (LABEL .EQ. 'FULLBAS ' .AND. MOD(I,2) .NE. 0)) THEN
           IHFO = IPLACE(I)
           CALL CC_SYMHFO(IHFO,ISYMHFO)
           WRITE(LUPRI,'(A,I3,A,I3,A,F10.4)') 
     *      ' Freezing HF-orbital ',IHFO,' of symmetry '
     *      ,ISYMHFO,' and with orbital energy',FOCDIA(IHFO)
           NRHFFR(ISYMHFO) = NRHFFR(ISYMHFO)+1
         END IF
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
            KFRRHF(I,ISYM) = I
         ENDDO

         IF (NVIRFR(ISYM) .GT. MAXFRO) THEN
            WRITE(LUPRI,'(//,1X,2A,I3)') 'ERROR: Maximum number of ',
     *           'frozen orbitals per symmetry is:',MAXFRO
            CALL QUIT('Too many frozen orbitals')
         END IF
         DO I = 1,NVIRFR(ISYM)
            JORB = NVIRS(ISYM) - I + 1
            KFRVIR(I,ISYM) = JORB
         ENDDO

      ENDDO    
C
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('CC_FREEZER')
C
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
      CALL QENTER('CC_SYMHFO(IHFO,ISYMHFO')
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
C
      CALL QEXIT('CC_SYMHFO(IHFO,ISYMHFO')
C
      RETURN
      END
C  /* Deck ccsd_cbs1 */
      SUBROUTINE CCSD_CBS1(T2AM,FCDIAG,ES,ET,QS,QT)
C
C     Written by Wim Klopper (University of Karlsruhe, 21 November 2002).
C
C     MP2 pair energies and CBS scaling factors
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
      PARAMETER (DP25 = 0.25D0, DP75 = 0.75D0)
#include <priunit.h>
      DIMENSION T2AM(*),FCDIAG(*),ES(*),ET(*),QS(*),QT(*)
#include <ccorb.h>
#include <ccsdsym.h>
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
      N12 = NRHFT * (NRHFT + 1)/2
      CALL DZERO(ES,N12)
      CALL DZERO(ET,N12)
      CALL DZERO(QS,N12)
      CALL DZERO(QT,N12)
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
                           NAIBJ = IT2AM(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
                           DENOM = ONE/(FCDIAG(KOFFI) + FCDIAG(KOFFJ)
     *                              -   FCDIAG(KOFFA) - FCDIAG(KOFFB))
                           IJ = INDEX(KOFFI,KOFFJ)
                           ISYMJA = MULD2H(ISYMJ,ISYMA) 
                           ISYMIB = MULD2H(ISYMI,ISYMB)
                           NAJ = IT1AM(ISYMA,ISYMJ) +
     *                           NVIR(ISYMA)*(J-1) + A
                           NBI = IT1AM(ISYMB,ISYMI) +
     *                           NVIR(ISYMB)*(I-1) + B
                           NAJBI = IT2AM(ISYMJA,ISYMIB) +
     *                             INDEX(NAJ,NBI)  
                           VAIBJ = T2AM(NAIBJ)
                           VAJBI = T2AM(NAJBI)
                           CS = ABS(VAIBJ + VAJBI)
                           CT = ABS(VAIBJ - VAJBI)
                           VS = CS**2
                           VT = CT**2
                           ES(IJ) = ES(IJ) + VS * DENOM
                           ET(IJ) = ET(IJ) + VT * DENOM 
                           QS(IJ) = QS(IJ) + CS * DENOM
                           QT(IJ) = QT(IJ) + CT * DENOM 
  160                   CONTINUE
  150                CONTINUE
  140             CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
      E2S = ZERO
      E2T = ZERO
      CALL AROUND('SECOND-ORDER PAIR ENERGIES')
      WRITE(LUPRI,'(4X,A8,4(7X,A8))') 
     *  '   I   J','T_2(s=0)','T_2(s=1)'
CQST *             'Q_2(s=0)','Q_2(s=1)'
      I = 0
      DO 230 KI=1,NRHFT
       DO 230 KJ=1,KI
        I = I + 1 
        ES(I) = ES(I) * DP25
        ET(I) = ET(I) * DP75
        E2S = E2S + ES(I)
        E2T = E2T + ET(I)
        WRITE(LUPRI,'(4X,2I4,4F15.9)') KI,KJ,ES(I),ET(I)
CQST *        (ONE+QS(I))**2,(ONE+QT(I))**2
  230 CONTINUE
      E2 = E2S + E2T
      WRITE(LUPRI,'(/A5,7X,2F15.9 )') ' SUM ',E2S,E2T 
      WRITE(LUPRI,'(/A5,7X, F15.9/)') ' TOT.',E2
      CALL FLSHFO(LUPRI)
      RETURN
      END
C  /* Deck ccsd_cbs2 */
      SUBROUTINE CCSD_CBS2(T1AM,T2AM,WORK,LWORK,
     *                     ET1S,ET1T,ET2S,ET2T,ETY)
C
C     Written by Wim Klopper (University of Karlsruhe, 21 November 2002).
C
C     Coupled-cluster pair energies and CBS scaling factors
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (DP5 = 0.5D0, D1P5 = 1.5D0)
#include <priunit.h>
#include <iratdef.h>
      DIMENSION T1AM(*),T2AM(*),WORK(*),
     *          ET1S(*),ET1T(*),ET2S(*),ET2T(*)
      CHARACTER*5 ETY
      LOGICAL LEXIST
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <ccfield.h>
#include <ccinftap.h>
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
      N12 = NRHFT * (NRHFT + 1)/2 
      CALL DZERO(ET1S,N12)
      CALL DZERO(ET1T,N12)
      CALL DZERO(ET2S,N12)
      CALL DZERO(ET2T,N12)
      KIAJB  = 1
      KEND1  = KIAJB + NT2AMX
      LWRK1  = LWORK - KEND1
      IF (LWRK1 .LT. 0) THEN
         STOP 'Insufficient spaces in CCSD_CBS2'
      ENDIF
      REWIND(LUIAJB)
      CALL READI(LUIAJB,IRAT*NT2AMX,WORK)
      DO 100 ISYMJ = 1,NSYM
         DO 110 ISYMB = 1,NSYM
            ISYMBJ = MULD2H(ISYMB,ISYMJ)
            ISYMAI = ISYMBJ
            DO 120 ISYMI = 1,NSYM
               ISYMBI = MULD2H(ISYMB,ISYMI)
               ISYMA  = MULD2H(ISYMI,ISYMAI)
               ISYMAJ = ISYMBI
               DO 130 J = 1,NRHF(ISYMJ)
                  KOFFJ = IRHF(ISYMJ) + J
                  DO 140 B = 1,NVIR(ISYMB)
                     KBJ = IT1AM(ISYMB,ISYMJ)
                     NBJ = KBJ + NVIR(ISYMB)*(J-1) + B
                     DO 150 I = 1,NRHF(ISYMI)
                        KOFFI = IRHF(ISYMI) + I
                        IJ = INDEX(KOFFI,KOFFJ)
                        KBI = IT1AM(ISYMB,ISYMI)
                        NBI = KBI + NVIR(ISYMB)*(I-1) + B
                        DO 160 A = 1,NVIR(ISYMA)
                           KAI = IT1AM(ISYMA,ISYMI)
                           NAI = KAI + NVIR(ISYMA)*(I-1) + A
                           KAJ = IT1AM(ISYMA,ISYMJ)
                           NAJ = KAJ + NVIR(ISYMA)*(J-1) + A
                           NAIBJ = IT2AM(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
                           NAJBI = IT2AM(ISYMAJ,ISYMBI) + INDEX(NAJ,NBI)
                           IF (ISYMB .EQ. ISYMJ) THEN
                              ET1S(IJ) = ET1S(IJ) + 
     *                                   (WORK(NAIBJ) + WORK(NAJBI)) *
     *                                   T1AM(NAI)*T1AM(NBJ)
                              ET1T(IJ) = ET1T(IJ) + 
     *                                   (WORK(NAIBJ) - WORK(NAJBI)) *
     *                                   T1AM(NAI)*T1AM(NBJ)
                           ENDIF
                           ET2S(IJ) = ET2S(IJ) + 
     *                                (WORK(NAIBJ) + WORK(NAJBI)) *
     *                                T2AM(NAIBJ)        
                           ET2T(IJ) = ET2T(IJ) + 
     *                                (WORK(NAIBJ) - WORK(NAJBI))*
     *                                T2AM(NAIBJ)        
  160                   CONTINUE
  150                CONTINUE
  140             CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
      EET1S = ZERO
      EET1T = ZERO
      EET2S = ZERO
      EET2T = ZERO
      CALL AROUND(ETY//' PAIR ENERGIES')
      WRITE(LUPRI,'(4X,A8,4(7X,A8))')
     *  '   I   J','T_1(s=0)','T_1(s=1)',
     *             'T_2(s=0)','T_2(s=1)'  
      I = 0
      DO 230 KI=1,NRHFT
       DO 230 KJ=1,KI        
        I = I + 1
        ET1S(I) = ET1S(I) * DP5
        ET2S(I) = ET2S(I) * DP5
        ET1T(I) = ET1T(I) * D1P5
        ET2T(I) = ET2T(I) * D1P5
        WRITE(LUPRI,'(4X,2I4,4F15.9)') KI,KJ,
     *  ET1S(I), ET1T(I), ET2S(I), ET2T(I)
        EET1S = EET1S + ET1S(I)
        EET1T = EET1T + ET1T(I)
        EET2S = EET2S + ET2S(I)
        EET2T = EET2T + ET2T(I)
  230 CONTINUE
      EE = EET1S + EET1T + EET2S + EET2T
      WRITE(LUPRI,'(/A5,7X,4F15.9 )') ' SUM ',EET1S,EET1T,EET2S,EET2T
      WRITE(LUPRI,'(/A5,7X, F15.9/)') ' TOT.',EE
      CALL FLSHFO(LUPRI)
      RETURN
      END
C  /* Deck frorhf */
      LOGICAL FUNCTION FRORHF(I,ISYM)
C
C     Thomas Bondo Pedersen, July 2003.
C
C     Returns .TRUE. if occupied orbital I of symmetry ISYM is frozen.
C
#include <implicit.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      FRORHF = .FALSE.

      IF (FROIMP) THEN

         IF (I .LE. NRHFFR(ISYM)) FRORHF = .TRUE.

      ELSE IF (FROEXP) THEN

         DO II = 1,NRHFFR(ISYM)
            IF (I .EQ. KFRRHF(II,ISYM)) THEN
               FRORHF = .TRUE.
               GO TO 100
            END IF
         END DO
  100    CONTINUE

      END IF

      IF (LOCDBG) THEN
         IF (FRORHF) THEN
            WRITE(LUPRI,'(A,I6,A,I2,A)')
     &      'Occupied orbital',I,' of sym.' ,ISYM,' is frozen'
         ELSE
            WRITE(LUPRI,'(A,I6,A,I2,A)')
     &      'Occupied orbital',I,' of sym.' ,ISYM,' is NOT frozen'
         END IF
      END IF

      RETURN
      END
C  /* Deck frovir */
      LOGICAL FUNCTION FROVIR(A,ISYM)
C
C     Thomas Bondo Pedersen, July 2003.
C
C     Returns .TRUE. if virtual orbital A of symmetry ISYM is frozen.
C
#include <implicit.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>

      INTEGER A
      INTEGER AA

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      FROVIR = .FALSE.

      IF (FROIMP) THEN

         IF (A .GT. NVIR(ISYM)) FROVIR = .TRUE.

      ELSE IF (FROEXP) THEN

         DO AA = 1,NVIRFR(ISYM)
            IF (A .EQ. KFRVIR(AA,ISYM)) THEN
               FROVIR = .TRUE.
               GO TO 100
            END IF
         END DO
  100    CONTINUE

      END IF

      IF (LOCDBG) THEN
         IF (FROVIR) THEN
            WRITE(LUPRI,'(A,I6,A,I2,A)')
     &      'Virtual  orbital',A,' of sym.' ,ISYM,' is frozen'
         ELSE
            WRITE(LUPRI,'(A,I6,A,I2,A)')
     &      'Virtual  orbital',A,' of sym.' ,ISYM,' is NOT frozen'
         END IF
      END IF

      RETURN
      END
