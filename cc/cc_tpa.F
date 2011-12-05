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
c /* deck CC_TPA */
*=====================================================================*
       SUBROUTINE CC_TPA(WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: direct calculation of second-order transition properties
*             (two-photon absorption strengths)
*             for the Coupled Cluster models
*
*                        CCS, CC2, CCSD, CC3
*
*     Written by Christof Haettig October 2003.
*
*=====================================================================*
      IMPLICIT NONE  
#include <priunit.h>
#include <cclists.h>
#include <cctpainf.h>
#include <ccsdinp.h>
#include <dummy.h>
#include <second.h>
#include <ccr1rsp.h>
#include <cclrmrsp.h>
#include <ccexcinf.h>
#include <ccorb.h>

* local parameters:
      CHARACTER*(16) MSGDBG
      PARAMETER (MSGDBG = '[debug] CC_TPA> ')

#if defined (SYS_CRAY)
      REAL ZERO
#else
      DOUBLE PRECISION ZERO
#endif
      PARAMETER (ZERO = 0.0d0)

      CHARACTER*10 MODEL
      INTEGER LWORK

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL TIM0, TIM1, TIMG, TIMFA, TIMF, TIMAA, 
     &                 TIMX2, TIMO1, TIMO2
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION TIM0, TIM1, TIMG, TIMFA, TIMF, TIMAA, 
     &                 TIMX2, TIMO1, TIMO2
#endif

      LOGICAL LADD
      INTEGER NSTATES, NBTPA, MXTRAN, MXVEC, KRESULT
      INTEGER KGTRAN,   KGDOTS,   KGCONS,   NGTRAN,
     &        KFATRAN,  KFADOTS,  KFACONS,  NFATRAN,
     &        KF1TRAN,  KF1DOTS,  KF1CONS,  NF1TRAN,
     &        KF2TRAN,  KF2DOTS,  KF2CONS,  NF2TRAN,
     &        KF3TRAN,  KF3DOTS,  KF3CONS,  NF3TRAN,
     &        KAA1TRAN, KAA1DOTS, KAA1CONS, NAA1TRAN,
     &        KAA2TRAN, KAA2DOTS, KAA2CONS, NAA2TRAN,
     &        KAA3TRAN, KAA3DOTS, KAA3CONS, NAA3TRAN,
     &        KXTRAN,   KXDOTS,   KXCONS,   NXTRAN,
     &        KO1TRAN,  KO1DOTS,  KO1CONS,  NO1TRAN,
     &        KO2TRAN,  KO2DOTS,  KO2CONS,  NO2TRAN
      INTEGER ISYM, KEND0, LEND0, IOPT, IORDER

* external functions:

*---------------------------------------------------------------------*
* print header for second-order property section:
*---------------------------------------------------------------------*
      WRITE (LUPRI,'(7(/1X,2A),/)')
     & '************************************',
     &                               '******************************',
     & '*                                   ',
     &                               '                             *',
     & '*<<<<<<    OUTPUT FROM COUPLED CLUST',
     &                               'ER LINEAR RESPONSE    >>>>>>>*',
     & '*                                   ',
     &                               '                             *',
     & '*<<<<<<  CALCULATION OF TWO-PHOTON A',
     &                               'ABSORPTION STRENGHTS  >>>>>>>*',
     & '*                                   ',
     &                               '                             *',
     & '************************************',
     &                               '******************************' 

*---------------------------------------------------------------------*
      IF (.NOT. (CCS .OR. CC2 .OR. CCSD .OR. CC3) ) THEN
         CALL QUIT('CC_TPA called for unknown Coupled Cluster.')
      END IF

* print some debug/info output
      IF (IPRINT .GT. 10) WRITE(LUPRI,*) 'CC_TPA Workspace:',LWORK
  
      TIM0  = SECOND()

*---------------------------------------------------------------------*
* allocate & initialize work space for property contributions:
*---------------------------------------------------------------------*
      ! number of excited states
      NSTATES = 0
      DO ISYM = 1, NSYM
        NSTATES = NSTATES + NCCEXCI(ISYM,1)
      END DO

      ! maximum number of transition moments to compute
      NBTPA    = 2 * NSMOPER * NSMSEL

      ! maximum number of transformations or vector calculations
      MXTRAN   = MAX(NSTATES,NLRM,NLRTLBL) * NLRTLBL

      ! maximum number of vectors to dot on 
      MXVEC    = MAX(NSTATES,NLRM,NLRTLBL)


      KRESULT  = 1
      KEND0    = KRESULT  + NBTPA
               
      KGTRAN   = KEND0
      KGDOTS   = KGTRAN   + MXTRAN * MXDIM_GTRAN
      KGCONS   = KGDOTS   + MXVEC  * MXTRAN
      KEND0    = KGCONS   + MXVEC  * MXTRAN

      KFATRAN  = KEND0
      KFADOTS  = KFATRAN  + MXTRAN * MXDIM_FATRAN
      KFACONS  = KFADOTS  + MXVEC  * MXTRAN
      KEND0    = KFACONS  + MXVEC  * MXTRAN

      KF1TRAN  = KEND0
      KF1DOTS  = KF1TRAN  + MXTRAN * MXDIM_FTRAN
      KF1CONS  = KF1DOTS  + MXVEC  * MXTRAN
      KEND0    = KF1CONS  + MXVEC  * MXTRAN

      KF2TRAN  = KEND0
      KF2DOTS  = KF2TRAN  + MXTRAN * MXDIM_FTRAN
      KF2CONS  = KF2DOTS  + MXVEC  * MXTRAN
      KEND0    = KF2CONS  + MXVEC  * MXTRAN

      KF3TRAN  = KEND0
      KF3DOTS  = KF3TRAN  + MXTRAN * MXDIM_FTRAN
      KF3CONS  = KF3DOTS  + MXVEC  * MXTRAN
      KEND0    = KF3CONS  + MXVEC  * MXTRAN

      KAA1TRAN = KEND0
      KAA1DOTS = KAA1TRAN + MXTRAN * MXDIM_XEVEC
      KAA1CONS = KAA1DOTS + MXVEC  * MXTRAN
      KEND0    = KAA1CONS + MXVEC  * MXTRAN

      KAA2TRAN = KEND0
      KAA2DOTS = KAA2TRAN + MXTRAN * MXDIM_XEVEC
      KAA2CONS = KAA2DOTS + MXVEC  * MXTRAN
      KEND0    = KAA2CONS + MXVEC  * MXTRAN

      KAA3TRAN = KEND0
      KAA3DOTS = KAA3TRAN + MXTRAN * MXDIM_XEVEC
      KAA3CONS = KAA3DOTS + MXVEC  * MXTRAN
      KEND0    = KAA3CONS + MXVEC  * MXTRAN

      KXTRAN   = KEND0
      KXDOTS   = KXTRAN   + MXTRAN
      KXCONS   = KXDOTS   + MXVEC  * MXTRAN
      KEND0    = KXCONS   + MXVEC  * MXTRAN

      KO1TRAN  = KEND0
      KO1DOTS  = KO1TRAN  + MXTRAN
      KO1CONS  = KO1DOTS  + MXVEC  * MXTRAN
      KEND0    = KO1CONS  + MXVEC  * MXTRAN

      KO2TRAN  = KEND0
      KO2DOTS  = KO2TRAN  + MXTRAN
      KO2CONS  = KO2DOTS  + MXVEC  * MXTRAN
      KEND0    = KO2CONS  + MXVEC  * MXTRAN

      LEND0 = LWORK - KEND0
      IF (LEND0 .LT. 0) THEN
        CALL QUIT('Insufficient memory in CC_TPA. (1)')
      END IF

      CALL DZERO(WORK(KRESULT),NBTPA)

      CALL DZERO(WORK(KGTRAN),  MXTRAN*(2*MXVEC+MXDIM_GTRAN))
      CALL DZERO(WORK(KFATRAN), MXTRAN*(2*MXVEC+MXDIM_FATRAN))
      CALL DZERO(WORK(KF1TRAN), MXTRAN*(2*MXVEC+MXDIM_FTRAN))
      CALL DZERO(WORK(KF2TRAN), MXTRAN*(2*MXVEC+MXDIM_FTRAN))
      CALL DZERO(WORK(KF3TRAN), MXTRAN*(2*MXVEC+MXDIM_FTRAN))
      CALL DZERO(WORK(KAA1TRAN),MXTRAN*(2*MXVEC+MXDIM_XEVEC))
      CALL DZERO(WORK(KAA2TRAN),MXTRAN*(2*MXVEC+MXDIM_XEVEC))
      CALL DZERO(WORK(KAA3TRAN),MXTRAN*(2*MXVEC+MXDIM_XEVEC))
      CALL DZERO(WORK(KXTRAN),  MXTRAN*(2*MXVEC+1))
      CALL DZERO(WORK(KO1TRAN), MXTRAN*(2*MXVEC+1))
      CALL DZERO(WORK(KO2TRAN), MXTRAN*(2*MXVEC+1))

*---------------------------------------------------------------------*
* set up lists for F transformations, ETA{O} and Xi{O} vectors:
*---------------------------------------------------------------------*
      LADD = .FALSE.

      CALL CCTPA_SETUP(MXTRAN,  MXVEC,
     &       WORK(KGTRAN),   WORK(KGDOTS),   WORK(KGCONS),   NGTRAN,
     &       WORK(KFATRAN),  WORK(KFADOTS),  WORK(KFACONS),  NFATRAN,
     &       WORK(KF1TRAN),  WORK(KF1DOTS),  WORK(KF1CONS),  NF1TRAN,
     &       WORK(KF2TRAN),  WORK(KF2DOTS),  WORK(KF2CONS),  NF2TRAN,
     &       WORK(KF3TRAN),  WORK(KF3DOTS),  WORK(KF3CONS),  NF3TRAN,
     &       WORK(KAA1TRAN), WORK(KAA1DOTS), WORK(KAA1CONS), NAA1TRAN,
     &       WORK(KAA2TRAN), WORK(KAA2DOTS), WORK(KAA2CONS), NAA2TRAN,
     &       WORK(KAA3TRAN), WORK(KAA3DOTS), WORK(KAA3CONS), NAA3TRAN,
     &       WORK(KXTRAN),   WORK(KXDOTS),   WORK(KXCONS),   NXTRAN,
     &       WORK(KO1TRAN),  WORK(KO1DOTS),  WORK(KO1CONS),  NO1TRAN,
     &       WORK(KO2TRAN),  WORK(KO2DOTS),  WORK(KO2CONS),  NO2TRAN,
     &       WORK(KRESULT),  LADD,           WORK(KEND0),    LEND0)


*---------------------------------------------------------------------*
* calculate G matrix contributions:
*---------------------------------------------------------------------*
      TIM1 = SECOND()
 
      IF (.NOT. LTPA_USE_X2) THEN
       IOPT = 5
       CALL CC_GMATRIX('L0 ','R1 ','R1 ','RE ',NGTRAN, MXVEC,
     &               WORK(KGTRAN),WORK(KGDOTS),WORK(KGCONS),
     &               WORK(KEND0), LEND0, IOPT )
      END IF

      TIMG = SECOND() - TIM1
      IF (NGTRAN.GT.0) WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")')
     &  '>>> Time used for',NGTRAN,' G matrix transformations:   ',TIMG
      CALL FLSHFO(LUPRI)

*---------------------------------------------------------------------*
* calculate F{O} matrix contributions:
*---------------------------------------------------------------------*
      TIM1 = SECOND()
 
      IF (.NOT.LTPA_USE_X2) THEN
       CALL CCQR_FADRV('L0 ','o1 ','RE ','R1 ',NFATRAN, MXVEC,
     &                  WORK(KFATRAN),WORK(KFADOTS),WORK(KFACONS),
     &                  WORK(KEND0), LEND0, 'DOTP' )
      END IF
 
      TIMFA = SECOND() - TIM1
      IF (NFATRAN.GT.0) WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")')
     &'>>> Time used for',NFATRAN,' F{O} matrix transformations:',TIMFA
      CALL FLSHFO(LUPRI)

*---------------------------------------------------------------------*
* calculate F matrix contributions:
*---------------------------------------------------------------------*
      TIM1 = SECOND()

      IF (.NOT. LTPA_USE_O2) THEN
        ! M1(w_f) x B x R1A(-w_A) x R1B(-w_B)
        IOPT = 5
        CALL CC_FMATRIX(WORK(KF1TRAN),NF1TRAN,'M1 ','R1 ',IOPT,'R1 ',
     &                  WORK(KF1DOTS),WORK(KF1CONS),MXVEC,
     &                  WORK(KEND0), LEND0)
      END IF

      IF (.NOT. LTPA_USE_O2) THEN
        ! LE(-w_f) x B x R1A(w_A) x R1B(w_B)
        IOPT = 5
        CALL CC_FMATRIX(WORK(KF2TRAN),NF2TRAN,'LE ','R1 ',IOPT,'R1 ',
     &                  WORK(KF2DOTS),WORK(KF2CONS),MXVEC,
     &                  WORK(KEND0), LEND0)
      END IF

      IF (.NOT. LTPA_USE_X2) THEN
        IOPT = 5
        CALL CC_FMATRIX(WORK(KF3TRAN),NF3TRAN,'L1 ','RE ',IOPT,'R1 ',
     &                  WORK(KF3DOTS),WORK(KF3CONS),MXVEC,
     &                  WORK(KEND0), LEND0)
      END IF

      TIMF = SECOND() - TIM1

      IF ((NF1TRAN+NF2TRAN+NF3TRAN).GT.0) 
     &  WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")')
     &   '>>> Time used for',NF1TRAN+NF2TRAN+NF3TRAN,
     &     ' F matrix transformations:',TIMF
      CALL FLSHFO(LUPRI)

*---------------------------------------------------------------------*
* calculate ETA{O} vector contributions:
*---------------------------------------------------------------------*
      TIM1 = SECOND()

      IF (.NOT. LTPA_USE_X2) THEN
        ! L1A(-w_A) x A{B} x RE(w_f)
        IOPT   = 5
        IORDER = 1
        CALL CC_XIETA(WORK(KAA3TRAN), NAA3TRAN, IOPT, IORDER, 'L1 ',
     &                '---',IDUMMY,        DUMMY,
     &                'RE ',WORK(KAA3DOTS),WORK(KAA3CONS),
     &                .FALSE.,MXVEC, WORK(KEND0), LEND0 )
      END IF

      IF (.NOT. LTPA_USE_O2) THEN
        ! LE(-w_f) x A{A} x R1B(w_B)
        IOPT   = 5
        IORDER = 1
        CALL CC_XIETA(WORK(KAA2TRAN), NAA2TRAN, IOPT, IORDER, 'LE ',
     &                '---',IDUMMY,        DUMMY,
     &                'R1 ',WORK(KAA2DOTS),WORK(KAA2CONS),
     &                .FALSE.,MXVEC, WORK(KEND0), LEND0 )
      END IF

      IF (.NOT. LTPA_USE_O2) THEN
        ! M1(w_f) x A{A} x R1B(-w_B)
        IOPT   = 5
        IORDER = 1
        CALL CC_XIETA(WORK(KAA1TRAN), NAA1TRAN, IOPT, IORDER, 'M1 ',
     &                '---',IDUMMY,        DUMMY,
     &                'R1 ',WORK(KAA1DOTS),WORK(KAA1CONS),
     &                .FALSE.,MXVEC, WORK(KEND0), LEND0 )
      END IF

      TIMAA = SECOND() - TIM1
      IF ((NAA1TRAN+NAA2TRAN+NAA3TRAN).GT.0) 
     &  WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")') 
     &   '>>> Time used for',NAA1TRAN+NAA2TRAN+NAA3TRAN,
     &       ' Eta{O} vector calculation:',TIMAA
      CALL FLSHFO(LUPRI)

*---------------------------------------------------------------------*
* CHI2 x RE vector contributions:
*---------------------------------------------------------------------*
      IF (LTPA_USE_X2) THEN
        TIM1 = SECOND()
 
        ! CHI2(-w_A,-w_B) x RE(w_f) 
        CALL CC_DOTDRV('X2 ','RE ',NXTRAN,MXVEC,
     &                 WORK(KXTRAN), WORK(KXDOTS), WORK(KXCONS),
     &                 WORK(KEND0), LEND0 )
 
        TIMX2 = SECOND() - TIM1
        WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")')
     &    '>>> Time used for',NXTRAN, ' X2 x RE dot products:',TIMX2
        CALL FLSHFO(LUPRI)
 
      END IF 

*---------------------------------------------------------------------*
* M1 x O2 and LE x O2 dot products:
*---------------------------------------------------------------------*
      IF (LTPA_USE_O2) THEN
        TIM1 = SECOND()
 
        CALL CC_DOTDRV('M1 ','O2 ',NO1TRAN,MXVEC,
     &                 WORK(KO1TRAN), WORK(KO1DOTS), WORK(KO1CONS),
     &                 WORK(KEND0), LEND0 )
 
        TIMO1 = SECOND() - TIM1
        WRITE (LUPRI,'(/A,I5,A,5X,F12.2," seconds.")')
     &    '>>> Time used for',NO1TRAN, ' M1 x O2 dot products:  ',TIMO1
        CALL FLSHFO(LUPRI)
 

        TIM1 = SECOND()
 
        CALL CC_DOTDRV('LE ','O2 ',NO2TRAN,MXVEC,
     &                 WORK(KO2TRAN), WORK(KO2DOTS), WORK(KO2CONS),
     &                 WORK(KEND0), LEND0 )
 
        TIMO2 = SECOND() - TIM1
        WRITE (LUPRI,'(/A,I5,A,5X,F12.2," seconds.")')
     &    '>>> Time used for',NO2TRAN, ' LE x O2 dot products:  ',TIMO2
        CALL FLSHFO(LUPRI)
 
      END IF

*---------------------------------------------------------------------*
* collect contributions and sum them up to the final results:
*---------------------------------------------------------------------*
      LADD = .TRUE.

cch
c     call dcopy(mxvec*mxtran,1.234d0,0,work(kaa1cons),1)
c     call dcopy(mxvec*mxtran,1.234d0,0,work(kaa2cons),1)
c     call dcopy(mxvec*mxtran,1.234d0,0,work(kaa3cons),1)
cch

      CALL CCTPA_SETUP(MXTRAN,  MXVEC,
     &       WORK(KGTRAN),   WORK(KGDOTS),   WORK(KGCONS),   NGTRAN,
     &       WORK(KFATRAN),  WORK(KFADOTS),  WORK(KFACONS),  NFATRAN,
     &       WORK(KF1TRAN),  WORK(KF1DOTS),  WORK(KF1CONS),  NF1TRAN,
     &       WORK(KF2TRAN),  WORK(KF2DOTS),  WORK(KF2CONS),  NF2TRAN,
     &       WORK(KF3TRAN),  WORK(KF3DOTS),  WORK(KF3CONS),  NF3TRAN,
     &       WORK(KAA1TRAN), WORK(KAA1DOTS), WORK(KAA1CONS), NAA1TRAN,
     &       WORK(KAA2TRAN), WORK(KAA2DOTS), WORK(KAA2CONS), NAA2TRAN,
     &       WORK(KAA3TRAN), WORK(KAA3DOTS), WORK(KAA3CONS), NAA3TRAN,
     &       WORK(KXTRAN),   WORK(KXDOTS),   WORK(KXCONS),   NXTRAN,
     &       WORK(KO1TRAN),  WORK(KO1DOTS),  WORK(KO1CONS),  NO1TRAN,
     &       WORK(KO2TRAN),  WORK(KO2DOTS),  WORK(KO2CONS),  NO2TRAN,
     &       WORK(KRESULT),  LADD,           WORK(KEND0),    LEND0)

*---------------------------------------------------------------------*
* print timing:
*---------------------------------------------------------------------*
      WRITE (LUPRI,'(/A,I4,A,F12.2," seconds.")') '>>> Total time for',
     &  NBTPA,' quadratic response func.:', SECOND() - TIM0

*---------------------------------------------------------------------*
* print one-photon absorption properties and return:
*---------------------------------------------------------------------*
      CALL  CCTPAPRT(WORK(KRESULT),.FALSE.,NSMOPER,NSMSEL)

      CALL FLSHFO(LUPRI)

      RETURN
      END

*=====================================================================*
*              END OF SUBROUTINE CC_TPA                               *
*=====================================================================*
c /* deck cctpaprt */
*=====================================================================*
      SUBROUTINE CCTPAPRT(TRANMOM,XST,NOPPAIRS,NTRANSIT)
*---------------------------------------------------------------------*
*
*    Purpose: print second-order transition moments and two-photon 
*             absorption properties calculated from these
*
*    Written by Christof Haettig in winter 2003/2004.
*
*=====================================================================*
      IMPLICIT NONE  
#include <priunit.h>
#include <pgroup.h>
#include <ccroper.h>
#include <cctpainf.h>
#include <ccsdinp.h>
#include <ccexci.h>
#include <codata.h>
#include <ccorb.h>


      LOGICAL XST
      INTEGER NOPPAIRS, NTRANSIT

#if defined (SYS_CRAY)
      REAL TRANMOM(2,NOPPAIRS,NTRANSIT)
      REAL ZERO, THIRTY, DELTAF, DELTAG, DELTAH
      REAL DELTA_TPlpa, DELTA_TPlpe, DELTA_TPc
      REAL DELTA_TPlpaE, DELTA_TPlpeE, DELTA_TPcE
      REAL AUTIME, K_1, K_2, K_3, K_lpa, K_lpe, K_c
      REAL SMFREQ, EIGV, SIGN, STRN, TMF0, TM0F
      REAL DIPLEN(3,3,2)
#else
      DOUBLE PRECISION TRANMOM(2,NOPPAIRS,NTRANSIT)
      DOUBLE PRECISION ZERO, THIRTY, DELTAF, DELTAG, DELTAH 
      DOUBLE PRECISION DELTA_TPlpa, DELTA_TPlpe, DELTA_TPc
      DOUBLE PRECISION DELTA_TPlpaE, DELTA_TPlpeE, DELTA_TPcE
      DOUBLE PRECISION AUTIME, K_1, K_2, K_3, K_lpa, K_lpe, K_c
      DOUBLE PRECISION SMFREQ, EIGV, SIGN, STRN, TMF0, TM0F
      DOUBLE PRECISION DIPLEN(3,3,2)
#endif
      PARAMETER (ZERO=0.0d0, THIRTY=30.0d0)

  
      CHARACTER*5  BLANKS
      CHARACTER*8  LABELA, LABELB
      CHARACTER*10 MODEL
      CHARACTER*8 LAB
      CHARACTER*80 STRING

      LOGICAL LDIPLEN, LDIPL(3,3)
      INTEGER IRSD, ISYME, IEXCI, ISTATE, IMULE, IOPA, IOPB, ISYMA,
     &        ISYMB, ISAMA, ISAMB, IOPPAIR, IDX, JDX, ISYMAB

#if defined (SYS_CRAY)
      REAL DDOT
#else
      DOUBLE PRECISION DDOT
#endif

*---------------------------------------------------------------------*
* print header for second-order properties: 
*---------------------------------------------------------------------*
      BLANKS = '     '
      STRING = ' RESULTS FOR TWO-PHOTON ABSORPTION STRENGTHS '

      IF (CCS) THEN
         CALL AROUND( BLANKS//'FINAL CCS'//STRING(1:45)//BLANKS ) 
         MODEL = 'CCS'
      ELSE IF (CC2) THEN
         CALL AROUND( BLANKS//'FINAL CC2'//STRING(1:45)//BLANKS )
         MODEL = 'CC2'
      ELSE IF (CCSD) THEN
         CALL AROUND( BLANKS//'FINAL CCSD'//STRING(1:45)//BLANKS )
         MODEL = 'CCSD'
      ELSE IF (CC3) THEN
         CALL AROUND( BLANKS//'FINAL CC3'//STRING(1:45)//BLANKS )
         MODEL = 'CC3'
      ELSE
         CALL QUIT('CCTPAPRT called for an unknown Coupled '//
     &             'Cluster model.')
      END IF


C     WRITE(LUPRI,*) 'XTANGM10:',XTANGM10
C     WRITE(LUPRI,*) 'XFSEC:',XFSEC
C     WRITE(LUPRI,*) 'ALPHA2:',ALPHA2
C     WRITE(LUPRI,*) 'a_0^4 t_0 alpha^2:',XTANGM10**4 * XFSEC * ALPHA2
C     WRITE(LUPRI,*) '(2*PI)**2:',(2.0d0 * PI)**2

      DO IRSD = 1, NTRANSIT
       
       IF (XST) THEN
        CALL QUIT('unfinished cas in CCTPAPRT.')
       ELSE

        ! ground to excited state transition:
        ISYME  = ISMSEL(IRSD,1)   ! irrep
        IEXCI  = ISMSEL(IRSD,2)   ! state number within irrep
        ISTATE = ISYOFE(ISYME) + IEXCI  ! number over all irreps
        SMFREQ = BSMFR(IRSD) 
        IMULE  = 1
        EIGV   = EIGVAL(ISTATE)
        
        WRITE(LUPRI,'(//3X,A,9X,A,A3,A,I5,"^",I1,A3)') 
     &  'transition',':   X ^1',REP(0),'  <--',IEXCI,IMULE,REP(ISYME-1)
        WRITE(LUPRI,'(3X,A,F15.10,A,F8.3,A,F9.0,A)')
     &    'excitation energy  :',EIGV, ' a.u.  ',
     &           EIGV*XTEV,' e.V.  ',EIGV*XTKAYS,' cm^-1'
        WRITE(LUPRI,'((3X,A,F15.10,A,F8.3,A,F9.0,A))')
     &    'photon energy for A:', SMFREQ,' a.u.  ',
     &       SMFREQ*XTEV,' e.V.  ',SMFREQ*XTKAYS,' cm^-1',
     &    'photon energy for B:', EIGV-SMFREQ,' a.u.  ',
     &      (EIGV-SMFREQ)*XTEV,' e.V.  ',(EIGV-SMFREQ)*XTKAYS,' cm^-1'

       END IF

       WRITE(LUPRI,'(/3x,a,a,/3x,a,3(f8.4,a),/3x,a,a)')
     &    '+-------------------+----------------+----------------+',
     &                                        '----------------+',
     &    '|    A        B     | M_0f(',-SMFREQ,') | M_f0(',SMFREQ,
     &           ') | S^0f(',SMFREQ,') |',
     &    '+-------------------+----------------+----------------+',
     &                                        '----------------+'

       DO IDX = 1, 3
        DO JDX = 1, 3
         LDIPL(IDX,JDX) = .FALSE.
        END DO 
       END DO 

       DO IOPPAIR = 1, NOPPAIRS
         IF (XST) THEN
           CONTINUE
         ELSE
           IOPA = IASMOP(IOPPAIR)
           IOPB = IBSMOP(IOPPAIR)
         END IF
         ISYMA  = ISYOPR(IOPA)
         ISYMB  = ISYOPR(IOPB)
         ISAMA  = ISYMAT(IOPA)
         ISAMB  = ISYMAT(IOPB)
         LABELA = LBLOPR(IOPA)
         LABELB = LBLOPR(IOPB)

         ISYMAB = MULD2H(ISYMA,ISYMB)

         IF (ISYMAB.EQ.ISYME) THEN
C          ! compute the effect of complex conjugation upon
C          ! the sign of M_0f and M_f0 (not used at the moment)
C          SIGN  = DBLE(ISAMA*ISAMB)
C          IF (ISAMA*ISAMB.EQ.0) SIGN = +ONE
          
           TMF0  = TRANMOM(1,IOPPAIR,IRSD)
           TM0F  = TRANMOM(2,IOPPAIR,IRSD)
           STRN  = TM0F * TMF0
          
           WRITE(LUPRI,
     &       '("   | ",2(a8,1x),"|",2(f15.8," |"),f15.8," |")')
     &        LABELA,LABELB,TM0F,TMF0,STRN

         ELSE
           TMF0  = ZERO
           TM0F  = ZERO
           STRN  = ZERO
          
           WRITE(LUPRI,
     &       '("   | ",2(a8,1x),"|",2(6x,a,7x,"|"),6x,a,7x"|")')
     &        LABELA,LABELB,'---','---','---'
         END IF

         
         ! collect transition moments for special operators

         ! -> here some dipole vectors:
         IF (LABELA(1:1).EQ.'X') IDX = 1
         IF (LABELA(1:1).EQ.'Y') IDX = 2
         IF (LABELA(1:1).EQ.'Z') IDX = 3
         IF (LABELB(1:1).EQ.'X') JDX = 1
         IF (LABELB(1:1).EQ.'Y') JDX = 2
         IF (LABELB(1:1).EQ.'Z') JDX = 3

         IF (LABELA(2:7).EQ.'DIPLEN' .AND. LABELB(2:7).EQ.'DIPLEN') THEN
            DIPLEN(IDX,JDX,1) = TM0F
            DIPLEN(IDX,JDX,2) = TMF0
            LDIPL( IDX,JDX)   = .TRUE.
         END IF

       END DO

       WRITE(LUPRI,'(3x,a,a)')
     &    '+-------------------+----------------+----------------+',
     &                                        '----------------+'

       LDIPLEN = .TRUE.
       DO IDX = 1, 3
         DO JDX = 1, 3
           LDIPLEN = LDIPLEN .AND. LDIPL(IDX,JDX)
         END DO
       END DO

       ! dipole strength strengths
       IF (LDIPLEN) THEN
         DELTAF = ZERO
         DELTAG = ZERO
         DELTAH = ZERO
         DO IDX = 1, 3
           DO JDX = 1, 3
             DELTAF = DELTAF + DIPLEN(IDX,IDX,1)*DIPLEN(JDX,JDX,2)
             DELTAG = DELTAG + DIPLEN(IDX,JDX,1)*DIPLEN(IDX,JDX,2)
             DELTAH = DELTAH + DIPLEN(IDX,JDX,1)*DIPLEN(JDX,IDX,2)
           END DO
         END DO
         DELTAF = DELTAF / THIRTY
         DELTAG = DELTAG / THIRTY
         DELTAH = DELTAH / THIRTY
       
         WRITE(LUPRI,'(6x,a,f15.8)')
     &     ' delta_F = S^0f_iijj/30   (length gauge)   : ',DELTAF,
     &     ' delta_G = S^0f_ijij/30   (length gauge)   : ',DELTAG,
     &     ' delta_H = S^0f_ijji/30   (length gauge)   : ',DELTAH
       END IF

*********************************************************************
*                                                                   *
*                 Compute rotationally averaged                     *
*                 two-photon transition strengths                   *
*                 for different light beam polarizations            *
*                                                                   *
*	Martin J. Paterson	22/02/05                            *
*                                                                   *
*********************************************************************
         WRITE(LUPRI,'(a)')
     &    ' ' 
         WRITE(LUPRI,'(6x,a)')
     &    ' Rotationally averaged two-photon transition strengths',
     &    '                 and rate constants'
 
         DELTA_TPlpa = ZERO 
         DELTA_TPlpa = 2 * (DELTAF + DELTAG + DELTAH)
         DELTA_TPlpe = ZERO
         DELTA_TPlpe = ((-1 * DELTAF ) + (4 * DELTAG) + DELTAH)
         DELTA_TPc = ZERO
         DELTA_TPc = ((-2 * DELTAF ) + (3 * DELTAG) + (3 * DELTAH))
         AUTIME = HBAR/XTJ 
         K_1 = ZERO
         K_1 = (XTANG * XTANG * XTANG * XTANG) * AUTIME  
         K_2 = ZERO
         K_2 = (8 * PI * PI *PI) * (ALPHAC * ALPHAC)
         K_3 = ZERO
         K_3 = (SMFREQ * (EIGV-SMFREQ))
         K_lpa = ZERO
         K_lpa = K_1 * K_2 * K_3 * DELTA_TPlpa
         K_lpe = ZERO
         K_lpe = K_1 * K_2 * K_3 * DELTA_TPlpe
         K_c = ZERO
         K_c = K_1 * K_2 * K_3 * DELTA_TPc 

         WRITE(LUPRI,'(a)')
     &    ' '
         WRITE(LUPRI,'(6x,a)')
     &    ' +-------------------+----------------+----------------+'
         WRITE(LUPRI,'(6x,a)')
     &    ' |   Polarization    |    DELTA_TP    |   K (0 -> f)   |'
         WRITE(LUPRI,'(6x,a)')
     &    ' +-------------------+----------------+----------------+'
         WRITE(LUPRI,'(6x,a,f13.8,a,E14.6,a)')
     &    ' |  linear (para)    |', DELTA_TPlpa, '   |', K_lpa,'  |'
         WRITE(LUPRI,'(6x,a)')
     &    ' +-------------------+----------------+----------------+'
         WRITE(LUPRI,'(6x,a,f13.8,a,E14.6,a)')
     &    ' |  linear (perp)    |', DELTA_TPlpe, '   |', K_lpe,'  |'
         WRITE(LUPRI,'(6x,a)')
     &    ' +-------------------+----------------+----------------+'
         WRITE(LUPRI,'(6x,a,f13.8,a,E14.6,a)')
     &    ' |    circular       |', DELTA_TPc,   '   |', K_c,'  |'
         WRITE(LUPRI,'(6x,a)')
     &    ' +-------------------+----------------+----------------+'

C
C        Output to property file, both Delta and w_1*w_2*Delta 
C
         LAB = 'D_pa    ' 
         CALL WRIPRO(DELTA_TPlpa,MODEL,-30,LAB,LAB,LAB,LAB,
     *               SMFREQ,EIGV-SMFREQ,EIGV,1,ISYME,1,ISTATE)
         LAB = 'D_pe    ' 
         CALL WRIPRO(DELTA_TPlpe,MODEL,-31,LAB,LAB,LAB,LAB,
     *               SMFREQ,EIGV-SMFREQ,EIGV,1,ISYME,1,ISTATE)
         LAB = 'D_pc    ' 
         CALL WRIPRO(DELTA_TPc,MODEL,-32,LAB,LAB,LAB,LAB,
     *               SMFREQ,EIGV-SMFREQ,EIGV,1,ISYME,1,ISTATE)
         DELTA_TPlpaE=SMFREQ*(EIGV-SMFREQ)*DELTA_TPlpa
         DELTA_TPlpeE=SMFREQ*(EIGV-SMFREQ)*DELTA_TPlpe 
         DELTA_TPcE=SMFREQ*(EIGV-SMFREQ)*DELTA_TPc
         LAB = 'w1w2D_pa' 
         CALL WRIPRO(DELTA_TPlpaE,MODEL,-33,LAB,LAB,LAB,LAB,
     *               SMFREQ,EIGV-SMFREQ,EIGV,1,ISYME,1,ISTATE)
         LAB = 'w1w2D_pe' 
         CALL WRIPRO(DELTA_TPlpeE,MODEL,-34,LAB,LAB,LAB,LAB,
     *               SMFREQ,EIGV-SMFREQ,EIGV,1,ISYME,1,ISTATE)
         LAB = 'w1w2D_pc' 
         CALL WRIPRO(DELTA_TPcE,MODEL,-35,LAB,LAB,LAB,LAB,
     *               SMFREQ,EIGV-SMFREQ,EIGV,1,ISYME,1,ISTATE)

*********************************************************************


      END DO ! IRSD

      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCTPAPRT                             *
*---------------------------------------------------------------------*
c /* deck cctpa_setup */
*=====================================================================*
      SUBROUTINE CCTPA_SETUP(MXTRAN,  MXVEC,
     &                       IGTRAN,   IGDOTS,   GCONS,   NGTRAN,
     &                       IFATRAN,  IFADOTS,  FACONS,  NFATRAN,
     &                       IF1TRAN,  IF1DOTS,  F1CONS,  NF1TRAN,
     &                       IF2TRAN,  IF2DOTS,  F2CONS,  NF2TRAN,
     &                       IF3TRAN,  IF3DOTS,  F3CONS,  NF3TRAN,
     &                       IAA1TRAN, IAA1DOTS, AA1CONS, NAA1TRAN,
     &                       IAA2TRAN, IAA2DOTS, AA2CONS, NAA2TRAN,
     &                       IAA3TRAN, IAA3DOTS, AA3CONS, NAA3TRAN,
     &                       IXTRAN,   IXDOTS,   XCONS,   NXTRAN,
     &                       IO1TRAN,  IO1DOTS,  O1CONS,  NO1TRAN,
     &                       IO2TRAN,  IO2DOTS,  O2CONS,  NO2TRAN,
     &                       RESULT,   LADD,     WORK,    LWORK )
*---------------------------------------------------------------------*
*
*    Purpose: set up for CC second-order transition moments
*
*     Written by Christof Haettig, Oct 2003
*
*=====================================================================*
      IMPLICIT NONE  
#include <priunit.h>
#include <cclists.h>
#include <cctpainf.h>
#include <ccroper.h>
#include <ccexci.h>
#include <ccsdinp.h>
#include <ccorb.h>

* local parameters:
      CHARACTER*(21) MSGDBG
      PARAMETER (MSGDBG = '[debug] CCTPA_SETUP> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      LOGICAL LADD
      INTEGER MXVEC, MXTRAN, LWORK

      INTEGER IGTRAN(MXDIM_GTRAN,MXTRAN)
      INTEGER IFATRAN(MXDIM_FATRAN,MXTRAN)
      INTEGER IF1TRAN(MXDIM_FTRAN,MXTRAN)
      INTEGER IF2TRAN(MXDIM_FTRAN,MXTRAN)
      INTEGER IF3TRAN(MXDIM_FTRAN,MXTRAN)
      INTEGER IAA1TRAN(MXDIM_XEVEC,MXTRAN)
      INTEGER IAA2TRAN(MXDIM_XEVEC,MXTRAN)
      INTEGER IAA3TRAN(MXDIM_XEVEC,MXTRAN)
      INTEGER IXTRAN(MXTRAN)
      INTEGER IO1TRAN(MXTRAN)
      INTEGER IO2TRAN(MXTRAN)

      INTEGER IGDOTS(MXVEC,MXTRAN)
      INTEGER IFADOTS(MXVEC,MXTRAN)
      INTEGER IF1DOTS(MXVEC,MXTRAN)
      INTEGER IF2DOTS(MXVEC,MXTRAN)
      INTEGER IF3DOTS(MXVEC,MXTRAN)
      INTEGER IAA1DOTS(MXVEC,MXTRAN)
      INTEGER IAA2DOTS(MXVEC,MXTRAN)
      INTEGER IAA3DOTS(MXVEC,MXTRAN)
      INTEGER IXDOTS(MXVEC,MXTRAN)
      INTEGER IO1DOTS(MXVEC,MXTRAN)
      INTEGER IO2DOTS(MXVEC,MXTRAN)

      INTEGER NGTRAN, NFATRAN, NF1TRAN, NF2TRAN, NF3TRAN, NAA1TRAN,
     &        NAA2TRAN, NAA3TRAN, NXTRAN, NO1TRAN, NO2TRAN

#if defined (SYS_CRAY)
      REAL RESULT(2,NSMOPER,NSMSEL)
      REAL GCONS(MXVEC,MXTRAN)
      REAL FACONS(MXVEC,MXTRAN)
      REAL F1CONS(MXVEC,MXTRAN)
      REAL F2CONS(MXVEC,MXTRAN)
      REAL F3CONS(MXVEC,MXTRAN)
      REAL AA1CONS(MXVEC,MXTRAN)
      REAL AA2CONS(MXVEC,MXTRAN)
      REAL AA3CONS(MXVEC,MXTRAN)
      REAL XCONS(MXVEC,MXTRAN)
      REAL O1CONS(MXVEC,MXTRAN)
      REAL O2CONS(MXVEC,MXTRAN)
      REAL WORK(LWORK)
      REAL ZERO, EIGV, FREQB, FREQA
      REAL GCON, FACON1, FACON2, F3CON1, F3CON2, F1CON,
     &     AA1CON1, AA1CON2, XCON, O1CON, TMAB0F,
     &     F2CON, AA2CON1, AA2CON2, O2CON, TMABF0
#else
      DOUBLE PRECISION RESULT(2,NSMOPER,NSMSEL)
      DOUBLE PRECISION GCONS(MXVEC,MXTRAN)
      DOUBLE PRECISION FACONS(MXVEC,MXTRAN)
      DOUBLE PRECISION F1CONS(MXVEC,MXTRAN)
      DOUBLE PRECISION F2CONS(MXVEC,MXTRAN)
      DOUBLE PRECISION F3CONS(MXVEC,MXTRAN)
      DOUBLE PRECISION AA1CONS(MXVEC,MXTRAN)
      DOUBLE PRECISION AA2CONS(MXVEC,MXTRAN)
      DOUBLE PRECISION AA3CONS(MXVEC,MXTRAN)
      DOUBLE PRECISION XCONS(MXVEC,MXTRAN)
      DOUBLE PRECISION O1CONS(MXVEC,MXTRAN)
      DOUBLE PRECISION O2CONS(MXVEC,MXTRAN)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION ZERO, EIGV, FREQB, FREQA
      DOUBLE PRECISION GCON, FACON1, FACON2, F3CON1, F3CON2, F1CON,
     &                 AA1CON1, AA1CON2, XCON, O1CON, TMAB0F,
     &                 F2CON, AA2CON1, AA2CON2, O2CON, TMABF0,
     &                 AA3CON1, AA3CON2
#endif
      PARAMETER (ZERO = 0.0D0)

      CHARACTER*(8) LABELA, LABELB
      LOGICAL LORXA, LORXB, LPDBSA, LPDBSB
      INTEGER IOPPAIR, IOPA, IOPB, ISYMA, ISYMB, ISYMAB, IRSD, ISTATE,
     &        NBTPA, IRE0, ILE0, IMBAR, IR1AP, IR1BP, IR1AM, IR1BM,
     &        IL1AM, IL1BM, IO2P, IO2M, IX2M, IDX, ITRAN, IVEC, I,
     &        MXGVEC, MXFAVEC, MXF1VEC, MXF2VEC, MXF3VEC, MXAA1VEC,
     &        MXAA2VEC, MXAA3VEC, MXXVEC, MXO1VEC, MXO2VEC

* external functions:
      INTEGER IR1TAMP
      INTEGER ILRMAMP
      INTEGER IL1ZETA
      INTEGER IRHSR2
      INTEGER ICHI2

*---------------------------------------------------------------------*
* initialize lists of transformations and dots products:
*---------------------------------------------------------------------*
      IF (.NOT. LADD) THEN

        NGTRAN   = 0
        NFATRAN  = 0
        NF1TRAN  = 0
        NF2TRAN  = 0
        NF3TRAN  = 0
        NAA1TRAN = 0
        NAA2TRAN = 0
        NAA3TRAN = 0
        NXTRAN   = 0
        NO1TRAN  = 0
        NO2TRAN  = 0

        DO ITRAN = 1, MXTRAN
         DO IDX = 1, MXDIM_GTRAN
           IGTRAN(IDX,ITRAN) = 0
         END DO
        
         DO IDX = 1, MXDIM_FATRAN
           IFATRAN(IDX,ITRAN) = 0
         END DO
        
         DO IDX = 1, MXDIM_FTRAN
           IF1TRAN(IDX,ITRAN) = 0
           IF2TRAN(IDX,ITRAN) = 0
           IF3TRAN(IDX,ITRAN) = 0
         END DO
        
         DO IDX = 1, MXDIM_XEVEC
           IAA1TRAN(1,ITRAN) = 0
           IAA2TRAN(1,ITRAN) = 0
           IAA3TRAN(1,ITRAN) = 0
         END DO
        
         IAA1TRAN(3,ITRAN)  = -1
         IAA2TRAN(3,ITRAN)  = -1
         IAA3TRAN(3,ITRAN)  = -1
         IAA1TRAN(4,ITRAN)  = -1
         IAA2TRAN(4,ITRAN)  = -1
         IAA3TRAN(4,ITRAN)  = -1
        
         IXTRAN(ITRAN)  = 0
         IO1TRAN(ITRAN) = 0
         IO2TRAN(ITRAN) = 0
        
        
         DO IVEC  = 1, MXVEC
          IGDOTS(IVEC,ITRAN)   = 0
          IFADOTS(IVEC,ITRAN)  = 0
          IF1DOTS(IVEC,ITRAN)  = 0
          IF2DOTS(IVEC,ITRAN)  = 0
          IF3DOTS(IVEC,ITRAN)  = 0
          IAA1DOTS(IVEC,ITRAN) = 0
          IAA2DOTS(IVEC,ITRAN) = 0
          IAA3DOTS(IVEC,ITRAN) = 0
          IXDOTS(IVEC,ITRAN)   = 0
          IO1DOTS(IVEC,ITRAN)  = 0
          IO2DOTS(IVEC,ITRAN)  = 0
         END DO
        END DO
        
      END IF
        
      MXGVEC   = 0
      MXFAVEC  = 0
      MXF1VEC  = 0
      MXF2VEC  = 0
      MXF3VEC  = 0
      MXAA1VEC = 0
      MXAA2VEC = 0
      MXAA3VEC = 0
      MXXVEC   = 0
      MXO1VEC  = 0
      MXO2VEC  = 0

*---------------------------------------------------------------------*
* initialize array with final results:
*---------------------------------------------------------------------*
      IF (LADD) THEN
        CALL DZERO(RESULT,NSMSEL*NSMOPER*2)
      END IF

*---------------------------------------------------------------------*
* print some information:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
        WRITE(LUPRI,*) 'LTPA_USE_O2:',LTPA_USE_O2
        WRITE(LUPRI,*) 'LTPA_USE_X2:',LTPA_USE_X2
      END IF

*---------------------------------------------------------------------*
* start loop over all requested transition moments: 
* note that for the transition strength S^0f_AB,AB(w_B) the two 
* moments M^AB_0f(-w_B) and M^AB_f0(w_B) are needed. it is assumed 
* that the operators A and B are either real or pure imaginary.
*---------------------------------------------------------------------*
      NBTPA = 0

      DO IOPPAIR = 1, NSMOPER
        IOPA   = IASMOP(IOPPAIR)
        IOPB   = IBSMOP(IOPPAIR)
        LORXA  = .FALSE.
        LORXB  = .FALSE.
        ISYMA  = ISYOPR(IOPA)
        ISYMB  = ISYOPR(IOPB)
        ISYMAB = MULD2H(ISYMA,ISYMB)
        LABELA = LBLOPR(IOPA)
        LABELB = LBLOPR(IOPB)
        LPDBSA = LPDBSOP(IOPA)
        LPDBSB = LPDBSOP(IOPB)

        IF (LPDBSA .OR. LPDBSB) CALL QUIT('perturbation-dependent '//
     &    'basis sets not implemented in CCTPA_SETUP.')


        DO IRSD = 1, NSMSEL
         IF (ISMSEL(IRSD,1) .EQ. ISYMAB) THEN
          ISTATE = ISYOFE(ISYMAB) + ISMSEL(IRSD,2)
          EIGV   = EIGVAL(ISTATE)
          FREQB  = BSMFR(IRSD)
          FREQA  = EIGV - FREQB

          NBTPA = NBTPA + 1

          IRE0  = ISYOFE(ISYMAB) + ISMSEL(IRSD,2)
          ILE0  = ISYOFE(ISYMAB) + ISMSEL(IRSD,2)
          IMBAR = ILRMAMP(ISTATE,EIGV,ISYMAB)

          IR1BP  = IR1TAMP(LABELB,.FALSE.,+FREQB,ISYMB)
          IR1AP  = IR1TAMP(LABELA,.FALSE.,+FREQA,ISYMA)
          IR1BM  = IR1TAMP(LABELB,.FALSE.,-FREQB,ISYMB)
          IR1AM  = IR1TAMP(LABELA,.FALSE.,-FREQA,ISYMA)

          IL1BM  = IL1ZETA(LABELB,.FALSE.,-FREQB,ISYMB)
          IL1AM  = IL1ZETA(LABELA,.FALSE.,-FREQA,ISYMA)

          IF (LTPA_USE_O2) THEN
           IO2P =IRHSR2(LABELA,.FALSE.,+FREQA,ISYMA,
     &                  LABELB,.FALSE.,+FREQB,ISYMB) 
           IO2M =IRHSR2(LABELA,.FALSE.,-FREQA,ISYMA,
     &                  LABELB,.FALSE.,-FREQB,ISYMB) 
          END IF

          IF (LTPA_USE_X2) THEN
           IX2M = ICHI2(LABELA,.FALSE.,-FREQA,ISYMA,
     &                  LABELB,.FALSE.,-FREQB,ISYMB)
          END IF

          IF (LOCDBG) THEN
            WRITE(LUPRI,*) 'ISYME,IEXCI,EIGV:',ISMSEL(IRSD,1),
     &                                    ISMSEL(IRSD,2),EIGV
            WRITE(LUPRI,*) 'LABELA,ISYMA,FREQA:',LABELA,ISYMA,FREQA
            WRITE(LUPRI,*) 'R1^A(+w_A):',IR1AP
            WRITE(LUPRI,*) 'R1^A(-w_A):',IR1AM
            WRITE(LUPRI,*) 'L1^A(-w_A):',IL1AM
            WRITE(LUPRI,*) 'LABELB,ISYMB,FREQB:',LABELB,ISYMB,FREQB
            WRITE(LUPRI,*) 'R1^B(+w_B):',IR1BP
            WRITE(LUPRI,*) 'R1^B(-w_B):',IR1BM
            WRITE(LUPRI,*) 'L1^B(-w_B):',IL1BM
            IF (LTPA_USE_O2) THEN
              WRITE(LUPRI,*) 'O2^AB(+w_A,+w_B):',IO2P
              WRITE(LUPRI,*) 'O2^AB(-w_A,-w_B):',IO2M
            END IF
            IF (LTPA_USE_X2) WRITE(LUPRI,*) 'X2^AB(-w_A,-w_B):',IX2M
          END IF
 
*---------------------------------------------------------------------*
*         G matrix transformation for M^AB_0f(-w_B): 
*            G x R1A(-w_A) x R1B(-w_B) x RE(w_f)
*---------------------------------------------------------------------*
          GCON = ZERO
          IF (.NOT. LTPA_USE_X2) THEN
            CALL CC_SETG112(IGTRAN,IGDOTS,MXTRAN,MXVEC,   
     &                      0,IR1AM,IR1BM,IRE0,ITRAN,IVEC)
            NGTRAN = MAX(NGTRAN,ITRAN)
            MXGVEC = MAX(MXGVEC,IVEC) 
            GCON   = GCONS(IVEC,ITRAN)
          END IF

*---------------------------------------------------------------------*
*         F{O} matrix transformations for M^AB_0f(-w_B): 
*            F{A} x RE(w_f) x R1B(-w_B)
*            F{B} x RE(w_f) x R1A(-w_A)
*---------------------------------------------------------------------*
          FACON1 = ZERO
          FACON2 = ZERO
          IF (.NOT. LTPA_USE_X2) THEN
            CALL CC_SETFA12(IFATRAN,IFADOTS,MXTRAN,MXVEC, 
     &                      0,IOPA,IRE0,IR1BM,ITRAN,IVEC)
            NFATRAN = MAX(NFATRAN,ITRAN)
            MXFAVEC = MAX(MXFAVEC,IVEC) 
            FACON1  = FACONS(IVEC,ITRAN)
            
            CALL CC_SETFA12(IFATRAN,IFADOTS,MXTRAN,MXVEC, 
     &                      0,IOPB,IRE0,IR1AM,ITRAN,IVEC)
            NFATRAN = MAX(NFATRAN,ITRAN)
            MXFAVEC = MAX(MXFAVEC,IVEC) 
            FACON2  = FACONS(IVEC,ITRAN)
          END IF

*---------------------------------------------------------------------*
*         B matrix transformations for M^AB_0f(-w_B): 
*            L1A(-w_A) x B x RE(w_f) x R1B(-w_B)
*            L1B(-w_B) x B x RE(w_f) x R1A(-w_A)
*         (computed as generalized F matrix contraction)
*---------------------------------------------------------------------*
          F3CON1 = ZERO
          F3CON2 = ZERO
          IF (.NOT. LTPA_USE_X2) THEN
            CALL CC_SETF12(IF3TRAN,IF3DOTS,MXTRAN,MXVEC, 
     &                     IL1AM,IRE0,IR1BM,ITRAN,IVEC)
            NF3TRAN = MAX(NF3TRAN,ITRAN)
            MXF3VEC = MAX(MXF3VEC,IVEC) 
            F3CON1  = F3CONS(IVEC,ITRAN)
            
            CALL CC_SETF12(IF3TRAN,IF3DOTS,MXTRAN,MXVEC, 
     &                     IL1BM,IRE0,IR1AM,ITRAN,IVEC)
            NF3TRAN = MAX(NF3TRAN,ITRAN)
            MXF3VEC = MAX(MXF3VEC,IVEC) 
            F3CON2  = F3CONS(IVEC,ITRAN)
          END IF

*---------------------------------------------------------------------*
*         A{O} matrix transformations for M^AB_0f(-w_B): 
*            L1A(-w_A) x A{B} x RE(w_f)
*            L1B(-w_B) x A{B} x RE(w_f)
*         (computed as generalized Eta vector dot product)
*---------------------------------------------------------------------*
          AA3CON1 = ZERO
          AA3CON2 = ZERO
          IF (.NOT. LTPA_USE_X2) THEN
            CALL CC_SETXE('Eta',IAA3TRAN,IAA3DOTS,MXTRAN,MXVEC,
     &                    IL1AM,IOPB,0,0,0,0,IRE0,ITRAN,IVEC)
            NAA3TRAN = MAX(NAA3TRAN,ITRAN)
            MXAA3VEC = MAX(MXAA3VEC,IVEC) 
            AA3CON1  = AA3CONS(IVEC,ITRAN)
            
            CALL CC_SETXE('Eta',IAA3TRAN,IAA3DOTS,MXTRAN,MXVEC,
     &                    IL1BM,IOPA,0,0,0,0,IRE0,ITRAN,IVEC)
            NAA3TRAN = MAX(NAA3TRAN,ITRAN)
            MXAA3VEC = MAX(MXAA3VEC,IVEC) 
            AA3CON2  = AA3CONS(IVEC,ITRAN)
          END IF

*---------------------------------------------------------------------*
*         CHI2(-w_A,-w_B) x RE(w_f) dot products for M^AB_0f(-w_B):
*         (replaces the above G, F{O}, B, and A{O} transformations 
*          if LTPA_USE_X2 is set to true)
*---------------------------------------------------------------------*
          XCON = ZERO
          IF (LTPA_USE_X2) THEN
            CALL CC_SETDOT(IXTRAN,IXDOTS,MXTRAN,MXVEC,
     &                     IX2M,IRE0,ITRAN,IVEC)
            NXTRAN = MAX(NXTRAN,ITRAN)
            MXXVEC = MAX(MXXVEC,IVEC) 
            XCON   = XCONS(IVEC,ITRAN)
          END IF

*---------------------------------------------------------------------*
*         B matrix transformation for M^AB_0f(-w_B): 
*            M1(w_f) x B x R1A(-w_A) x R1B(-w_B)
*         (computed as generalized F matrix contraction)
*---------------------------------------------------------------------*
          F1CON = ZERO
          IF (.NOT. LTPA_USE_O2) THEN
            CALL CC_SETF12(IF1TRAN,IF1DOTS,MXTRAN,MXVEC, 
     &                     IMBAR,IR1AM,IR1BM,ITRAN,IVEC)
            NF1TRAN = MAX(NF1TRAN,ITRAN)
            MXF1VEC = MAX(MXF1VEC,IVEC) 
            F1CON   = F1CONS(IVEC,ITRAN)
          END IF

*---------------------------------------------------------------------*
*         A{O} matrix transformations for M^AB_0f(-w_B): 
*            M1(w_f) x A{A} x R1B(-w_B)
*            M1(w_f) x A{B} x R1A(-w_A)
*         (computed as generalized Eta vector dot product)
*---------------------------------------------------------------------*
          AA1CON1 = ZERO
          AA1CON2 = ZERO

          IF (.NOT. LTPA_USE_O2) THEN
            CALL CC_SETXE('Eta',IAA1TRAN,IAA1DOTS,MXTRAN,MXVEC,
     &                    IMBAR,IOPA,0,0,0,0,IR1BM,ITRAN,IVEC)
            NAA1TRAN = MAX(NAA1TRAN,ITRAN)
            MXAA1VEC = MAX(MXAA1VEC,IVEC) 
            AA1CON1  = AA1CONS(IVEC,ITRAN)

            CALL CC_SETXE('Eta',IAA1TRAN,IAA1DOTS,MXTRAN,MXVEC,
     &                    IMBAR,IOPB,0,0,0,0,IR1AM,ITRAN,IVEC)
            NAA1TRAN = MAX(NAA1TRAN,ITRAN)
            MXAA1VEC = MAX(MXAA1VEC,IVEC) 
            AA1CON2  = AA1CONS(IVEC,ITRAN)
          END IF

*---------------------------------------------------------------------*
*         Mbar(w_f) x O2(-w_B,-w_A) dot products for M^AB_0f(-w_B):
*         (replaces above B and A{O} transformations for LTPA_USE_O2)
*---------------------------------------------------------------------*
          O1CON = ZERO
          IF (LTPA_USE_O2) THEN
            CALL CC_SETDOT(IO1TRAN,IO1DOTS,MXTRAN,MXVEC,
     &                     IMBAR,IO2M,ITRAN,IVEC)
            NO1TRAN = MAX(NO1TRAN,ITRAN)
            MXO1VEC = MAX(MXO1VEC,IVEC) 
            O1CON   = O1CONS(IVEC,ITRAN)
          END IF

*---------------------------------------------------------------------*
*         B matrix transformations for M^AB_f0(w_B): 
*            LE(-w_f) x B x R1A(w_A) x R1B(w_B)
*         (computed as generalized F matrix contraction)
*---------------------------------------------------------------------*
          F2CON = ZERO
          IF (.NOT. LTPA_USE_O2) THEN
            CALL CC_SETF12(IF2TRAN,IF2DOTS,MXTRAN,MXVEC, 
     &                     ILE0,IR1AP,IR1BP,ITRAN,IVEC)
            NF2TRAN = MAX(NF2TRAN,ITRAN)
            MXF2VEC = MAX(MXF2VEC,IVEC) 
            F2CON   = F2CONS(IVEC,ITRAN)
          END IF

*---------------------------------------------------------------------*
*         A{O} matrix transformations for M^AB_f0(w_B): 
*            LE(-w_f) x A{A} x R1B(w_B)
*            LE(-w_f) x A{B} x R1A(w_A)
*         (computed as generalized Eta vector dot product)
*---------------------------------------------------------------------*
          AA2CON1 = ZERO
          AA2CON2 = ZERO
          IF (.NOT. LTPA_USE_O2) THEN
            CALL CC_SETXE('Eta',IAA2TRAN,IAA2DOTS,MXTRAN,MXVEC,
     &                    ILE0,IOPA,0,0,0,0,IR1BP,ITRAN,IVEC)
            NAA2TRAN = MAX(NAA2TRAN,ITRAN)
            MXAA2VEC = MAX(MXAA2VEC,IVEC) 
            AA2CON1  = AA2CONS(IVEC,ITRAN)
            
            CALL CC_SETXE('Eta',IAA2TRAN,IAA2DOTS,MXTRAN,MXVEC,
     &                    ILE0,IOPB,0,0,0,0,IR1AP,ITRAN,IVEC)
            NAA2TRAN = MAX(NAA2TRAN,ITRAN)
            MXAA2VEC = MAX(MXAA2VEC,IVEC) 
            AA2CON2  = AA2CONS(IVEC,ITRAN)
          END IF

*---------------------------------------------------------------------*
*         LE(-w_f) x O2(w_A,w_B) dot products for M^AB_f0(w_B):
*         (replaces above B and A{O} transf. for LTPA_USE_O2=.true.)
*---------------------------------------------------------------------*
          O2CON = ZERO
          IF (LTPA_USE_O2) THEN
            CALL CC_SETDOT(IO2TRAN,IO2DOTS,MXTRAN,MXVEC,
     &                     ILE0,IO2P,ITRAN,IVEC)
            NO2TRAN = MAX(NO2TRAN,ITRAN)
            MXO2VEC = MAX(MXO2VEC,IVEC) 
            O2CON   = O2CONS(IVEC,ITRAN)
          END IF

*---------------------------------------------------------------------*
*       range checking:
*---------------------------------------------------------------------*
        IF (NGTRAN.GT.MXTRAN .OR. MXGVEC.GT.MXVEC) 
     &   CALL QUIT('Out of range for G trans. in CCTPA_SETUP')
        IF (NFATRAN.GT.MXTRAN .OR. MXFAVEC.GT.MXVEC) 
     &   CALL QUIT('Out of range for F{O} trans. in CCTPA_SETUP')
        IF (NF1TRAN.GT.MXTRAN .OR. MXF1VEC.GT.MXVEC) 
     &   CALL QUIT('Out of range for gen. F trans. (1) in CCTPA_SETUP')
        IF (NF2TRAN.GT.MXTRAN .OR. MXF2VEC.GT.MXVEC) 
     &   CALL QUIT('Out of range for gen. F trans. (2) in CCTPA_SETUP')
        IF (NF3TRAN.GT.MXTRAN .OR. MXF3VEC.GT.MXVEC) 
     &   CALL QUIT('Out of range for gen. F trans. (3) in CCTPA_SETUP')
        IF (NAA1TRAN.GT.MXTRAN .OR. MXAA1VEC.GT.MXVEC) 
     &   CALL QUIT('Out of range for gen. A{O} tr. (1) in CCTPA_SETUP')
        IF (NAA2TRAN.GT.MXTRAN .OR. MXAA2VEC.GT.MXVEC) 
     &   CALL QUIT('Out of range for gen. A{O} tr. (2) in CCTPA_SETUP')
        IF (NAA3TRAN.GT.MXTRAN .OR. MXAA3VEC.GT.MXVEC) 
     &   CALL QUIT('Out of range for gen. A{O} tr. (3) in CCTPA_SETUP')
        IF (NXTRAN.GT.MXTRAN .OR. MXXVEC.GT.MXVEC) 
     &   CALL QUIT('Out of range for Chi2 vec in CCTPA_SETUP')
        IF (NO1TRAN.GT.MXTRAN .OR. MXO1VEC.GT.MXVEC) 
     &   CALL QUIT('Out of range for O2 (1) vec in CCTPA_SETUP')
        IF (NO2TRAN.GT.MXTRAN .OR. MXO2VEC.GT.MXVEC) 
     &   CALL QUIT('Out of range for O2 (2) vec in CCTPA_SETUP')

*---------------------------------------------------------------------*
*       accumulate M^AB_0f(-w_B) and M^AB_f0(w_B):
*---------------------------------------------------------------------*
        IF (LADD) THEN

           TMABF0 = F2CON + AA2CON1 + AA2CON2 + O2CON

           TMAB0F = XCON + O1CON + GCON + FACON1 + FACON2 + 
     &              F3CON1 + F3CON2 + AA3CON1 + AA3CON2 + 
     &              F1CON + AA1CON1 + AA1CON2

           RESULT(1,IOPPAIR,IRSD) = TMABF0
           RESULT(2,IOPPAIR,IRSD) = TMAB0F

           IF (LOCDBG) THEN
             WRITE(LUPRI,'(3A)') 'OPERATORS:',LABELA,LABELB
             WRITE(LUPRI,'(A,I5,3F18.6)') 'ISTATE,EIGV,FREQA,FREQB:',
     &                                     ISTATE,EIGV,FREQA,FREQB
             WRITE(LUPRI,*) 'GCON:',GCON
             WRITE(LUPRI,*) 'FACON:',FACON1,FACON2
             WRITE(LUPRI,*) 'F3CON:',F3CON1,F3CON2
             WRITE(LUPRI,*) 'AA3CON:',AA3CON1,AA3CON2
             WRITE(LUPRI,*) 'F1CON:',F1CON
             WRITE(LUPRI,*) 'AA1CON:',AA1CON1,AA1CON2
             WRITE(LUPRI,*) 'XCON,O1CON:',XCON,O1CON
             WRITE(LUPRI,*) 'TMAB0F:',TMAB0F
             WRITE(LUPRI,*) 'F2CON:',F2CON
             WRITE(LUPRI,*) 'AA2CON:',AA2CON1,AA2CON2
             WRITE(LUPRI,*) 'O2CON:',O2CON
             WRITE(LUPRI,*) 'TMABF0:',TMABF0
           END IF

        END IF

*---------------------------------------------------------------------*
*       end loop over transition moments
*---------------------------------------------------------------------*
         END IF
        END DO
      END DO

*---------------------------------------------------------------------*
* print the lists: 
*---------------------------------------------------------------------*
* general statistics:
      IF ((.NOT.LADD) .OR. LOCDBG) THEN
       WRITE(LUPRI,'(/,/3X,A,I3,A)') 'For the requested',NBTPA,
     &      ' two-photon transition moments'
       WRITE(LUPRI,'((8X,A,I3,A))') 
     & ' - ',NF1TRAN, ' generalized F matrix transf. M1 x B x R1',
     & ' - ',NAA1TRAN,' generalized Eta vector calc. M1 x A{O}',
     & ' - ',NF2TRAN, ' generalized F matrix transf. LE x B x R1',
     & ' - ',NAA2TRAN,' generalized Eta vector calc. LE x A{O}',
     & ' - ',NGTRAN,  ' G matrix transformations     G  x R1 x R1',
     & ' - ',NFATRAN, ' F{O} matrix transformations  F{O} x RE',
     & ' - ',NF3TRAN, ' generalized F matrix transf. L1 x B x RE',
     & ' - ',NAA3TRAN,' generalized Eta vector calc. L1 x A{O}',
     & ' - ',NO1TRAN, ' M1 x O2 dot products',
     & ' - ',NO2TRAN, ' LE x O2 dot products',
     & ' - ',NXTRAN,  ' RE x X2 dot products'
       WRITE(LUPRI,'(3X,A,/,/)') 'will be performed.'
      END IF

      IF (LOCDBG) THEN

         ! G matrix transformations:
         WRITE(LUPRI,*)'List G matrix transf. G x R1 x R1 x RE:'
         DO ITRAN = 1, NGTRAN
           WRITE(LUPRI,'(A,3I5,5X,(25I3,20X))') MSGDBG,
     &      (IGTRAN(I,ITRAN),I=1,3),(IGDOTS(I,ITRAN),I=1,MXGVEC)
         END DO
         WRITE(LUPRI,*)

         ! F{O} matrix transformations:
         WRITE(LUPRI,*)'List F{O} matrix transf. F{O} x RE x R1:'
         DO ITRAN = 1, NFATRAN
           WRITE(LUPRI,'(A,3I5,5X,(25I3,20X))') MSGDBG,
     &      (IFATRAN(I,ITRAN),I=1,3),(IFADOTS(I,ITRAN),I=1,MXFAVEC)
         END DO
         WRITE(LUPRI,*)

         ! F matrix transformations:
         WRITE(LUPRI,*)'List gen. F matrix transf. M1 x B x R1 x R1:'
         DO ITRAN = 1, NF1TRAN
           WRITE(LUPRI,'(A,2I5,5X,(25I3,20X))') MSGDBG,
     &      (IF1TRAN(I,ITRAN),I=1,2),(IF1DOTS(I,ITRAN),I=1,MXF1VEC)
         END DO
         WRITE(LUPRI,*)

         ! F matrix transformations:
         WRITE(LUPRI,*)'List gen. F matrix transf. RE x B x R1 x R1:'
         DO ITRAN = 1, NF2TRAN
           WRITE(LUPRI,'(A,2I5,5X,(25I3,20X))') MSGDBG,
     &      (IF2TRAN(I,ITRAN),I=1,2),(IF2DOTS(I,ITRAN),I=1,MXF2VEC)
         END DO
         WRITE(LUPRI,*)

         ! F matrix transformations:
         WRITE(LUPRI,*)'List gen. F matrix transf. L1 x B x RE x R1:'
         DO ITRAN = 1, NF3TRAN
           WRITE(LUPRI,'(A,2I5,5X,(25I3,20X))') MSGDBG,
     &      (IF3TRAN(I,ITRAN),I=1,2),(IF3DOTS(I,ITRAN),I=1,MXF3VEC)
         END DO
         WRITE(LUPRI,*)

         ! ETA{O} vector calculations:
         WRITE(LUPRI,*) 'List of gen. ETA{O} vec. calc. M1 x A{O} x R1:'
         DO ITRAN = 1, NAA1TRAN
           WRITE(LUPRI,'(A,5I5,5X,(25I3,20X))') MSGDBG,
     &      (IAA1TRAN(I,ITRAN),I=1,5),(IAA1DOTS(I,ITRAN),I=1,MXAA1VEC)
         END DO
         WRITE(LUPRI,*)

         ! ETA{O} vector calculations:
         WRITE(LUPRI,*) 'List of gen. ETA{O} vec. calc. LE x A{O} x R1:'
         DO ITRAN = 1, NAA2TRAN
           WRITE(LUPRI,'(A,5I5,5X,(25I3,20X))') MSGDBG,
     &      (IAA2TRAN(I,ITRAN),I=1,5),(IAA2DOTS(I,ITRAN),I=1,MXAA2VEC)
         END DO
         WRITE(LUPRI,*)

         ! ETA{O} vector calculations:
         WRITE(LUPRI,*) 'List of gen. ETA{O} vec. calc. L1 x A{O} x RE:'
         DO ITRAN = 1, NAA3TRAN
           WRITE(LUPRI,'(A,5I5,5X,(25I3,20X))') MSGDBG,
     &      (IAA3TRAN(I,ITRAN),I=1,5),(IAA3DOTS(I,ITRAN),I=1,MXAA3VEC)
         END DO
         WRITE(LUPRI,*)

         ! extra O2 vector calculations:
         WRITE(LUPRI,*) 'List of extra O2 vec. calc. for dotting on M1:'
         DO ITRAN = 1, NO1TRAN
           WRITE(LUPRI,'(A,I5,5X,(25I3,20X))') MSGDBG,
     &      IO1TRAN(ITRAN),(IO1DOTS(I,ITRAN),I=1,MXO1VEC)
         END DO
         WRITE(LUPRI,*)

         ! extra O2 vector calculations:
         WRITE(LUPRI,*) 'List of extra O2 vec. calc. for dotting on LE:'
         DO ITRAN = 1, NO2TRAN
           WRITE(LUPRI,'(A,I5,5X,(25I3,20X))') MSGDBG,
     &      IO2TRAN(ITRAN),(IO2DOTS(I,ITRAN),I=1,MXO2VEC)
         END DO
         WRITE(LUPRI,*)

         ! extra X2 vector calculations:
         WRITE(LUPRI,*) 'List of extra X2 vec. calc. for dotting on RE:'
         DO ITRAN = 1, NXTRAN
           WRITE(LUPRI,'(A,I5,5X,(25I3,20X))') MSGDBG,
     &      IXTRAN(ITRAN),(IXDOTS(I,ITRAN),I=1,MXXVEC)
         END DO
         WRITE(LUPRI,*)

      END IF

      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCTPA_SETUP                          *
*---------------------------------------------------------------------*
