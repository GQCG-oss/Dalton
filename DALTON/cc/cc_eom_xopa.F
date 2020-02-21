!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2015 (2015), see http://daltonprogram.org"
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
C
c /* deck CC_EOM_XOPA */
*=====================================================================*
       SUBROUTINE CC_EOM_XOPA(WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: direct calculation of first-order transition properties
*             (transition moments and oscillator strengths)
*             for transitions between two excited states with the 
*             Coupled Cluster models
*
*                        CCS, CC2, CCSD, CC3
*
*             and partially with SCF and CIS
*
*     Written by Christof Haettig winter 2002/2003.
*     Modified version by Sonia, to include EOM models, 2015
*=====================================================================*
      IMPLICIT NONE  
#include "priunit.h"
#include "cclists.h"
#include "ccxopainf.h"
#include "ccsdinp.h"
#include "dummy.h"
#include "second.h"
#include "ccexcinf.h"
#include "ccorb.h"

* local parameters:
      CHARACTER*(16) MSGDBG
      PARAMETER (MSGDBG = '[debug] CC_EOM_XOPA> ')

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
      REAL TIM0, TIM1, TIMF, TIMXE1, TIMXE2
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION TIM0, TIM1, TIMF, TIMXE1, TIMXE2
#endif

      LOGICAL LADD
      INTEGER NBOPA, MXFTRAN, MXATRAN, MXXTRAN, MXFVEC, MXAVEC, MXXVEC,
     &        NFTRAN, NXE1TRAN, NXE2TRAN, NSTATES,
     &        KRESULT, KFTRAN, KFDOTS, KFCONS, KEND0, LEND0,
     &        KE1TRAN, KE1DOTS, KE1CONS,
     &        KX2TRAN, KX2DOTS, KX2CONS,
     &        IOPT, IORDER, ISYM
      INTEGER KEOMX2TRAN, KEOMX2DOTS, KEOMX2CONS 
      INTEGER KEOML0TRAN, KEOML0DOTS, KEOML0CONS
      INTEGER NEOMXE2TRAN, MXEOMXTRAN, MXEOMXVEC
      INTEGER NEOML0TRAN,MXEOML0TRAN,MXEOML0VEC

* external functions: none

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
     & '*<<<<<<  CALCULATION OF ONE-PHOTON A',
     &                               'BSORPTION STRENGTHS  >>>>>>>*',
     & '*<<<<<<     FOR EXCITED TO EXCITED S',
     &                               'TATE TRANSITIONS      >>>>>>>*',
     & '*                                   ',
     &                               '                             *',
     & '************************************',
     &                               '******************************' 

*---------------------------------------------------------------------*
      IF (.NOT. (CCS .OR. CC2 .OR. CCSD .OR. CC3) ) THEN
         CALL QUIT('CC_EOM_XOPA called for unknown Coupled Cluster.')
      END IF

* print some debug/info output
      IF(IPRINT .GT. 10) WRITE(LUPRI,*) 'CC_EOM_XOPA Workspace:',LWORK
  
      TIM0  = SECOND()

*---------------------------------------------------------------------*
* allocate & initialize work space for property contributions:
*---------------------------------------------------------------------*
      ! maximum number of transition moments to compute
      NBOPA   = 2 * NQR2OP * NXQR2ST

      ! number of excited states
      NSTATES = 0
      DO ISYM = 1, NSYM
        NSTATES = NSTATES + NCCEXCI(ISYM,1) + NCCEXCI(ISYM,3)
      END DO

      ! maximum number of transformations or vector calculations
      ! NSTATES * NQR2OP   LE x Eta{X} transformations
      ! NQR2OP             Xi{X} vectors
      ! 2*NXQR2ST          LE x B x RE transformations   
      MXATRAN    = NSTATES * NQR2OP
      MXXTRAN    = NQR2OP
      MXEOMXTRAN = NQR2OP
      MXEOML0TRAN  = 1
      MXFTRAN    = 2*NXQR2ST

      ! maximum number of vectors to dot on 
      ! NSTATES    RE vectors dotted on a LE x Eta{X} transformation
      ! 2*NXQR2ST  N2 vectors dotted on a Xi{X} vector
      ! NQR2OP     R1 vectors dotted on a LE x B x RE transformation
      MXAVEC      = NSTATES
      MXXVEC      = 2*NXQR2ST
      MXEOMXVEC   = 2*NXQR2ST !I am not sure I understand why this number...
      MXFVEC      = NQR2OP
      MXEOML0VEC  = NSTATES

      KRESULT  = 1
      KEND0    = KRESULT  + NBOPA
               
      KFTRAN   = KEND0
      KFDOTS   = KFTRAN   + MXFTRAN * MXDIM_FTRAN
      KFCONS   = KFDOTS   + MXFVEC  * MXFTRAN
      KEND0    = KFCONS   + MXFVEC  * MXFTRAN

      KE1TRAN  = KEND0
      KE1DOTS  = KE1TRAN  + MXATRAN * MXDIM_XEVEC
      KE1CONS  = KE1DOTS  + MXAVEC  * MXATRAN
      KEND0    = KE1CONS  + MXAVEC  * MXATRAN

      KX2TRAN  = KEND0
      KX2DOTS  = KX2TRAN  + MXXTRAN * MXDIM_XEVEC
      KX2CONS  = KX2DOTS  + MXXVEC  * MXXTRAN
      KEND0    = KX2CONS  + MXXVEC  * MXXTRAN

      KEOMX2TRAN= KEND0
      KEOMX2DOTS= KEOMX2TRAN  + MXEOMXTRAN * MXDIM_XEVEC
      KEOMX2CONS= KEOMX2DOTS  + MXEOMXVEC * MXEOMXTRAN
      KEND0     = KEOMX2CONS  + MXEOMXVEC * MXEOMXTRAN    

      KEOML0TRAN= KEND0
      KEOML0DOTS= KEOML0TRAN  + MXEOML0TRAN
      KEOML0CONS= KEOML0DOTS  + MXEOML0VEC * MXEOML0TRAN
      KEND0     = KEOML0CONS  + MXEOML0VEC * MXEOML0TRAN

      LEND0 = LWORK - KEND0
      IF (LEND0 .LT. 0) THEN
        CALL QUIT('Insufficient memory in CC_EOM_XOPA. (1)')
      END IF

      CALL DZERO(WORK(KRESULT),NBOPA)

      !sonia
      CALL DZERO(WORK(KEOMX2CONS),MXEOMXVEC * MXEOMXTRAN)

*---------------------------------------------------------------------*
* set up lists for F transformations, ETA{O} and Xi{O} vectors:
*---------------------------------------------------------------------*
      LADD = .FALSE.

      CALL CCXOPA_EOMSETUP(WORK(KFTRAN),WORK(KFDOTS),WORK(KFCONS),
     &                  NFTRAN,MXFTRAN,MXFVEC,
     &                  WORK(KE1TRAN),WORK(KE1DOTS),WORK(KE1CONS),
     &                  NXE1TRAN,MXATRAN,MXAVEC,
     &                  WORK(KX2TRAN),WORK(KX2DOTS),WORK(KX2CONS),
     &                  NXE2TRAN,MXXTRAN,MXXVEC,
     &                  WORK(KEOMX2TRAN),WORK(KEOMX2DOTS),
     &                  WORK(KEOMX2CONS),
     &                  NEOMXE2TRAN,MXEOMXTRAN,MXEOMXVEC,
     &                  WORK(KEOML0TRAN),WORK(KEOML0DOTS),
     &                  WORK(KEOML0CONS),
     &                  NEOML0TRAN,MXEOML0TRAN,MXEOML0VEC,
     &                  WORK(KRESULT),NBOPA,LADD,WORK(KEND0),LEND0)

*---------------------------------------------------------------------*
* calculate F matrix contributions:
*---------------------------------------------------------------------*
      TIM1 = SECOND()

      CALL DZERO(WORK(KFCONS),MXFVEC*NFTRAN)

      IOPT = 5
      CALL CC_FMATRIX(WORK(KFTRAN),NFTRAN,'LE ','RE ',IOPT,'R1 ',
     &                WORK(KFDOTS),WORK(KFCONS),MXFVEC,
     &                WORK(KEND0), LEND0)

      TIMF = SECOND() - TIM1

      IF (NFTRAN.GT.0) WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")')
     & '>>> Time used for',NFTRAN,' F matrix transformations:',TIMF
      CALL FLSHFO(LUPRI)

*---------------------------------------------------------------------*
* calculate LE x A{O} x RE contributions:
*---------------------------------------------------------------------*
      TIM1 = SECOND()

      CALL DZERO(WORK(KE1CONS),MXAVEC*NXE1TRAN)



      IOPT   = 5
      IORDER = 1
      if (leomxopa) then
         CALL CCEOM_XIETA( WORK(KE1TRAN),NXE1TRAN,IOPT,IORDER,'LE ',
     &                 '---',IDUMMY,       DUMMY,
     &                 'RE ',WORK(KE1DOTS),WORK(KE1CONS),
     &                       MXAVEC, WORK(KEND0), LEND0 )
      else
       CALL CC_XIETA( WORK(KE1TRAN), NXE1TRAN, IOPT, IORDER, 'LE ',
     &               '---',DUMMY,DUMMY,
     &               'RE ',WORK(KE1DOTS),WORK(KE1CONS),
     &               .FALSE.,MXAVEC, WORK(KEND0), LEND0 )
      end if

      TIMXE1 = SECOND() - TIM1
      IF (NXE1TRAN.GT.0) WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")') 
     & '>>> Time used for',NXE1TRAN,' A{X} matrix transformations:',
     & TIMXE1
      CALL FLSHFO(LUPRI)

*---------------------------------------------------------------------*
* calculate N2 x Xksi{O} vector contributions:
*---------------------------------------------------------------------*
      TIM1 = SECOND()

      CALL DZERO(WORK(KX2CONS),MXXVEC*NXE2TRAN)

      IOPT   = 5
      IORDER = 1
      CALL CC_XIETA( WORK(KX2TRAN), NXE2TRAN, IOPT, IORDER, '---',
     &               'N2 ',WORK(KX2DOTS),WORK(KX2CONS),
     &               '---',IDUMMY,DUMMY,
     &               .FALSE.,MXXVEC, WORK(KEND0), LEND0 )

      TIMXE2 = SECOND() - TIM1
      IF (NXE2TRAN.GT.0) WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")') 
     & '>>> Time used for',NXE2TRAN,' O1/X1 vector calculation:',TIMXE2
      CALL FLSHFO(LUPRI)

*---------------------------------------------------------------------*
* calculate LE x Xksi{O} vector contributions:
* calculate RE x Tbar0   vector contributions:
* multiply them together to get the EOM contrib
*---------------------------------------------------------------------*
      TIM1 = SECOND()

      if (leomxopa) then
         CALL DZERO(WORK(KEOMX2CONS),MXEOMXVEC*NEOMXE2TRAN)
         CALL DZERO(WORK(KEOML0CONS),MXEOML0VEC*NEOML0TRAN)

         IOPT   = 5
         IORDER = 1
         CALL CC_XIETA( WORK(KEOMX2TRAN), NEOMXE2TRAN, 
     &                IOPT, IORDER, 'L0 ',
     &               'LE ',WORK(KEOMX2DOTS),WORK(KEOMX2CONS),
     &               '---',IDUMMY,DUMMY,
     &               .FALSE.,MXEOMXVEC, WORK(KEND0), LEND0 )

         CALL CC_DOTDRV('L0 ','RE ',NEOML0TRAN,MXEOML0VEC,
     &                 WORK(KEOML0TRAN), WORK(KEOML0DOTS), 
     &                 WORK(KEOML0CONS),
     &                 WORK(KEND0), LEND0 )

         TIMXE2 = SECOND() - TIM1
      IF (NEOMXE2TRAN.GT.0) WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")') 
     & '>>> Time used for',NEOMXE2TRAN,
     &  ' O1/X1 vector calculation:',TIMXE2
         CALL FLSHFO(LUPRI)
      end if

*---------------------------------------------------------------------*
* collect contributions and sum them up to the final results:
*---------------------------------------------------------------------*
      LADD = .TRUE.

      CALL CCXOPA_EOMSETUP(WORK(KFTRAN),WORK(KFDOTS),WORK(KFCONS),
     &                  NFTRAN,MXFTRAN,MXFVEC,
     &                  WORK(KE1TRAN),WORK(KE1DOTS),WORK(KE1CONS),
     &                  NXE1TRAN,MXATRAN,MXAVEC,
     &                  WORK(KX2TRAN),WORK(KX2DOTS),WORK(KX2CONS),
     &                  NXE2TRAN,MXXTRAN,MXXVEC,
     &                  WORK(KEOMX2TRAN),WORK(KEOMX2DOTS),
     &                  WORK(KEOMX2CONS),
     &                  NEOMXE2TRAN,MXEOMXTRAN,MXEOMXVEC,
     &                  WORK(KEOML0TRAN),WORK(KEOML0DOTS),
     &                  WORK(KEOML0CONS),
     &                  NEOML0TRAN,MXEOML0TRAN,MXEOML0VEC,
     &                  WORK(KRESULT),NBOPA,LADD,WORK(KEND0),LEND0)

*---------------------------------------------------------------------*
* print timing:
*---------------------------------------------------------------------*
      WRITE (LUPRI,'(/A,I4,A,F12.2," seconds.")') '>>> Total time for',
     &  NBOPA,' quadratic response func.:', SECOND() - TIM0

*---------------------------------------------------------------------*
* print one-photon absorption properties and return:
*---------------------------------------------------------------------*
      CALL  CCOPAPRT(WORK(KRESULT),.TRUE.,NQR2OP,NXQR2ST)

      CALL FLSHFO(LUPRI)

      RETURN
      END

*=====================================================================*
*              END OF SUBROUTINE CC_EOM_XOPA                              *
*=====================================================================*
