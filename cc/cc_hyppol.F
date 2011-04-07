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
c /* deck cc_hyppol */
*=====================================================================*
       SUBROUTINE CC_HYPPOL(WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: direct calculation of frequency-dependent 
*             first hyperpolarizabilities for the Couple Cluster models
*
*                        CCS, CC2, CCSD
*
*     Written by Christof Haettig summer/fall 1996.
*     Dispersion coefficients, october 1997, Christof Haettig.
*     Orbital relax. for one operator, june 1999, Christof Haettig.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccqrinf.h>
#include <ccr1rsp.h>
#include <ccrc1rsp.h>
#include <ccroper.h>
#include <cclists.h>
#include <second.h>
#include <dummy.h>
#include <ccsections.h>
#include <ccslvinf.h>

* local parameters:
      CHARACTER*(19) MSGDBG
      PARAMETER (MSGDBG = '[debug] CC_HYPPOL> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

#if defined (SYS_CRAY)
      REAL ZERO
#else
      DOUBLE PRECISION ZERO
#endif
      PARAMETER (ZERO = 0.0d0)

      CHARACTER*10 MODEL
      CHARACTER*8  LABSOP
      LOGICAL LADD
      LOGICAL LORXA, LORXB, LORXC, LPDBSA, LPDBSB, LPDBSC
      INTEGER LWORK
      INTEGER ISYMA, ITAMPA, IZETVA, IOPERA, IKAPA, IR2AB
      INTEGER ISYMB, ITAMPB, IZETVB, IOPERB, IKAPB, IR2AC
      INTEGER ISYMC, ITAMPC, IZETVC, IOPERC, IKAPC, IR2BC, IOPSOP
      INTEGER KHYPPOL, NBHYPPOL, KEND0, LEND0, KEND1, LEND1
      INTEGER KEXPCOF, NBEXPCOF, NBINOM, KBINOM, KDSPCF
      INTEGER KDSPCFA, KSHGCFA, KORCOFA, KEOPECA, KABCDE
      INTEGER KSHGCF, KORCOF, KEOPEC, ISGNSOP
      INTEGER MXTRAN, MXVEC, ITRAN, IVEC, ISIGN, IORDER, ISYSOP, INUM
      INTEGER IOPT, IDX, IFREQ, ISYM, ISYML, ISYML1, ISYML2, IOPER
      INTEGER KGTRAN, KGDOTS, KGCONS, NGTRAN, MXGTRAN, MXGDOTS
      INTEGER KFTRAN, KFDOTS, KFCONS, NFTRAN, MXFTRAN, MXFDOTS
      INTEGER KBTRAN, KBDOTS, KBCONS, NBTRAN, MXBTRAN, MXBDOTS
      INTEGER KKTRAN, KKDOTS, KKCONS, NKTRAN, MXKTRAN, MXKDOTS
      INTEGER KFATRAN, KFADOTS, KFACONS, NFATRAN, MXFATRAN, MXFADOTS
      INTEGER KEATRAN, KEADOTS, KEACONS, NEATRAN, MXEATRAN, MXEADOTS
      INTEGER KR2TRAN, KR2DOTS, KR2CONS, NR2TRAN, MXR2TRAN, MXR2DOTS
      INTEGER KXETRAN, KEDOTS,  KECONS,  NXETRAN, MXXETRAN, MXXEDOTS
      INTEGER KAATRAN, KAADOTS, KAACONS, NAATRAN, MXAATRAN, MXAADOTS
      INTEGER KXDOTS, KXCONS

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL GCON, FACON1,FACON2,FACON3, FCON1,FCON2,FCON3
      REAL EACON1, EACON2, EACON3, EACON4, EACON5, EACON6
      REAL BCON1, BCON2, BCON3, R2CON1, R2CON2, R2CON3
      REAL SLVKCON1, SLVKCON2, SLVKCON3
      REAL X2CON1, X2CON2, X2CON3, E2CON1, E2CON2, E2CON3
      REAL SIGN
      REAL TIMTOT, TIM0, TIM1, TIMG, TIMF, TIMB, TIMFA, TIMEA, TIM2
      REAL TIMR2, TIMXE, TIMAA, TIMK
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION GCON, FACON1,FACON2,FACON3, FCON1,FCON2,FCON3
      DOUBLE PRECISION EACON1, EACON2, EACON3, EACON4, EACON5, EACON6
      DOUBLE PRECISION BCON1, BCON2, BCON3, R2CON1, R2CON2, R2CON3
      DOUBLE PRECISION SLVKCON1, SLVKCON2, SLVKCON3
      DOUBLE PRECISION X2CON1, X2CON2, X2CON3, E2CON1, E2CON2, E2CON3
      DOUBLE PRECISION SIGN
      DOUBLE PRECISION TIMTOT,TIM0,TIM1,TIMG,TIMF,TIMB,TIMFA,TIMEA,TIM2
      DOUBLE PRECISION TIMR2, TIMXE, TIMAA, TIMK
#endif

* external functions
      INTEGER IR1TAMP
      INTEGER IR1KAPPA
      INTEGER IL1ZETA
      INTEGER IR2TAMP
      INTEGER IROPER

*---------------------------------------------------------------------*
* print header for hyperpolarizability section
*---------------------------------------------------------------------*
      WRITE (LUPRI,'(7(/1X,2A),/)')
     & '************************************',
     &                               '*******************************',
     & '*                                   ',
     &                               '                              *',
     & '*<<<<<<<< OUTPUT FROM COUPLED CLUSTER Q',
     &                                  'UADRATIC RESPONSE >>>>>>>>>*',
     & '*                                   ',
     &                               '                              *',
     & '*<<<<<<<<     CALCULATION OF CC HYPER',
     &                                'POLARIZABILITIES    >>>>>>>>>*',
     & '*                                   ',
     &                               '                              *',
     & '************************************',
     &                               '*******************************' 

*---------------------------------------------------------------------*
      IF (.NOT. (CCS .OR. CC2 .OR. CCSD .OR. CC3) ) THEN
         CALL QUIT('CC_HYPPOL called for unknown Coupled Cluster.')
      END IF

* print some debug/info output
      IF (IPRINT .GT. 10) WRITE(LUPRI,*) 'CC_HYPPOL Workspace:',LWORK
  
      TIM0  = SECOND()

      call FLSHFO(LUPRI)

*---------------------------------------------------------------------*
* allocate & initialize work space for hyperpolarizabilities
*---------------------------------------------------------------------*
      NBHYPPOL = NQROPER * NQRFREQ

      MXTRAN  = 2 * NLRTLBL * NLRTLBL
      MXVEC   = NLRTLBL

      MXGTRAN  = MXDIM_GTRAN  * MXTRAN
      MXFTRAN  = MXDIM_FTRAN  * MXTRAN
      MXBTRAN  = MXDIM_BTRAN  * MXTRAN
      MXFATRAN = MXDIM_FATRAN * MXTRAN
      MXAATRAN = MXDIM_AATRAN * MXTRAN
      MXEATRAN = MXDIM_XEVEC  * MXTRAN
      MXR2TRAN = 1 * MXTRAN
      MXXETRAN = MXDIM_XEVEC  * MXTRAN
      MXKTRAN  = MXDIM_KTRAN  * MXTRAN

      MXGDOTS  = MXVEC * MXTRAN
      MXFDOTS  = MXVEC * MXTRAN
      MXBDOTS  = MXVEC * MXTRAN
      MXFADOTS = MXVEC * MXTRAN
      MXEADOTS = MXVEC * MXTRAN
      MXR2DOTS = MXVEC * MXTRAN
      MXXEDOTS = MXVEC * MXTRAN
      MXAADOTS = MXVEC * MXTRAN
      MXKDOTS  = MXVEC * MXTRAN

      KHYPPOL = 1
      KGTRAN  = KHYPPOL + 2 * NBHYPPOL
      KGDOTS  = KGTRAN  + MXGTRAN
      KGCONS  = KGDOTS  + MXGDOTS
      KFTRAN  = KGCONS  + MXGDOTS
      KFDOTS  = KFTRAN  + MXFTRAN
      KFCONS  = KFDOTS  + MXFDOTS
      KBTRAN  = KFCONS  + MXFDOTS
      KBDOTS  = KBTRAN  + MXBTRAN
      KBCONS  = KBDOTS  + MXBDOTS
      KKTRAN  = KBCONS  + MXBDOTS
      KKDOTS  = KKTRAN  + MXKTRAN
      KKCONS  = KKDOTS  + MXKDOTS
      KFATRAN = KKCONS  + MXKDOTS
      KFADOTS = KFATRAN + MXFATRAN
      KFACONS = KFADOTS + MXFADOTS
      KEATRAN = KFACONS + MXFADOTS
      KEADOTS = KEATRAN + MXEATRAN
      KEACONS = KEADOTS + MXEADOTS
      KR2TRAN = KEACONS + MXEADOTS
      KR2DOTS = KR2TRAN + MXR2TRAN
      KR2CONS = KR2DOTS + MXR2DOTS
      KXETRAN = KR2CONS + MXR2DOTS
      KXDOTS  = KXETRAN + MXXETRAN
      KXCONS  = KXDOTS  + MXXEDOTS
      KEDOTS  = KXCONS  + MXXEDOTS
      KECONS  = KEDOTS  + MXXEDOTS
      KAATRAN = KECONS  + MXXEDOTS
      KAADOTS = KAATRAN + MXAATRAN
      KAACONS = KAADOTS + MXAADOTS
      KEND0   = KAACONS + MXAADOTS
      LEND0   = LWORK   - KEND0

      IF (LEND0 .LT. 0) THEN
        CALL QUIT('Insufficient memory in CC_HYPPOL. (1)')
      END IF

      CALL DZERO(WORK(KHYPPOL),2*NBHYPPOL)

*---------------------------------------------------------------------*
* set up lists for G, F and F{A} transformations and ETA{O} vectors:
*---------------------------------------------------------------------*
      CALL CCQR_SETUP(MXTRAN, MXVEC,
     &                WORK(KGTRAN), WORK(KGDOTS), NGTRAN,
     &                WORK(KFTRAN), WORK(KFDOTS), NFTRAN,
     &                WORK(KBTRAN), WORK(KBDOTS), NBTRAN,
     &                WORK(KKTRAN), WORK(KKDOTS), NKTRAN,
     &                WORK(KFATRAN),WORK(KFADOTS),NFATRAN,
     &                WORK(KAATRAN),WORK(KAADOTS),NAATRAN,
     &                WORK(KEATRAN),WORK(KEADOTS),NEATRAN,
     &                WORK(KR2TRAN),WORK(KR2DOTS),NR2TRAN,
     &                WORK(KXETRAN),WORK(KXDOTS),WORK(KEDOTS),NXETRAN,
     &                WORK(KEND0), LEND0)

*---------------------------------------------------------------------*
* calculate G matrix contributions:
*---------------------------------------------------------------------*

      TIM1 = SECOND()

      IOPT = 5
      CALL CC_GMATRIX('L0 ','R1 ','R1 ','R1 ',NGTRAN, MXVEC,
     &              WORK(KGTRAN),WORK(KGDOTS),WORK(KGCONS),
     &              WORK(KEND0), LEND0, IOPT )


      TIMG = SECOND() - TIM1

      WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")')
     &  '>>> Time used for',NGTRAN,' G matrix transformations:',TIMG
      CALL FLSHFO(LUPRI)
      
      
*---------------------------------------------------------------------*
* calculate B matrix contributions:
*---------------------------------------------------------------------*
      IF ( USEBTRAN .AND. (.NOT.USE_R2) ) THEN
        TIM1 = SECOND()

        IOPT = 5
        CALL CC_BMATRIX(WORK(KBTRAN),NBTRAN,'R1 ','R1 ',IOPT,'L1 ',
     &               WORK(KBDOTS),WORK(KBCONS),MXVEC,
     &               .FALSE.,WORK(KEND0),LEND0)

        TIMB = SECOND() - TIM1

        WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")')
     &   '>>> Time used for',NBTRAN,' B matrix transformations:',TIMB
        CALL FLSHFO(LUPRI)

      ELSE
        TIMB = 0.0d0
      END IF
*---------------------------------------------------------------------*
* calculate K matrix contributions:
*---------------------------------------------------------------------*
      IF (CCSLV) THEN
        TIM1 = SECOND()

        IOPT = 5
        CALL CC_KMATRIX(WORK(KKTRAN),NKTRAN,'L1 ','L1 ',IOPT,'R1 ',
     &               WORK(KKDOTS),WORK(KKCONS),MXVEC,WORK(KEND0),LEND0)

        TIMK = SECOND() - TIM1

        WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")')
     &   '>>> Time used for',NKTRAN,' K matrix transformations:',TIMK
        CALL FLSHFO(LUPRI)

      ELSE
        TIMK = 0.0d0
      END IF

*---------------------------------------------------------------------*
* calculate F matrix contributions:
*---------------------------------------------------------------------*
      IF ( (.NOT.USEBTRAN) .AND. (.NOT.USE_R2) ) THEN
        TIM1 = SECOND()

        IOPT = 5
        CALL CC_FMATRIX(WORK(KFTRAN),NFTRAN,'L1 ','R1 ',IOPT,'R1 ',
     &                  WORK(KFDOTS),WORK(KFCONS),MXVEC,
     &                  WORK(KEND0), LEND0)

        TIMF = SECOND() - TIM1

        WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")')
     &   '>>> Time used for',NFTRAN,' F matrix transformations:',TIMF
        CALL FLSHFO(LUPRI)

      ELSE
        TIMF = 0.0d0
      END IF

*---------------------------------------------------------------------*
* calculate F{O} matrix contributions:
*---------------------------------------------------------------------*
      TIM1 = SECOND()

      CALL CCQR_FADRV('L0 ','o1 ','R1 ','R1 ',NFATRAN, MXVEC,
     &                 WORK(KFATRAN),WORK(KFADOTS),WORK(KFACONS),
     &                 WORK(KEND0), LEND0, 'DOTP' )
  
      TIMFA = SECOND() - TIM1

      WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")')
     & '>>> Time used for',NFATRAN,' F{O} matrix transform.:  ',TIMFA
      CALL FLSHFO(LUPRI)

*---------------------------------------------------------------------*
* calculate ETA{O} vector contributions:
*---------------------------------------------------------------------*
      CALL DZERO(WORK(KEACONS),MXEADOTS)

      IF ((.NOT. USE_R2) .AND. (.NOT. USE_AAMAT)) THEN
        TIM1 = SECOND()

        IOPT   = 5
        IORDER = 1
        CALL CC_XIETA( WORK(KEATRAN), NEATRAN, IOPT, IORDER, 'L1 ',
     &                 '---',IDUMMY,       DUMMY,
     &                 'R1 ',WORK(KEADOTS),WORK(KEACONS),
     &                 .FALSE.,MXVEC, WORK(KEND0), LEND0 )

        TIMEA = SECOND() - TIM1
        WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")')
     &       '>>> Time used for',
     &    NEATRAN,' ETA{O} vector calculat.: ', SECOND() - TIM1
        CALL FLSHFO(LUPRI)

      ELSE
        TIMEA = 0.0d0
      END IF

*---------------------------------------------------------------------*
* calculate A{O} matrix contributions:
*---------------------------------------------------------------------*
      CALL DZERO(WORK(KAACONS),MXAADOTS)

      IF ((.NOT. USE_R2) .AND. USE_AAMAT) THEN
        TIM1 = SECOND()

        IOPT   = 5
        CALL CC_AAMAT(WORK(KAATRAN),NAATRAN,'o1 ','R1 ',IOPT,'L1 ',
     &                WORK(KAADOTS),WORK(KAACONS),
     &                MXVEC, WORK(KEND0), LEND0 )        

        TIMAA = SECOND() - TIM1
        WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")')
     &       '>>> Time used for',
     &    NAATRAN,' A{O} matrix transform. : ', SECOND() - TIM1
        CALL FLSHFO(LUPRI)

      ELSE
        TIMAA = 0.0d0
      END IF

*---------------------------------------------------------------------*
* calculate X1 x R2 dots products:
*---------------------------------------------------------------------*
      IF ( USE_R2 ) THEN
        TIM1 = SECOND()

C       CALL CC_DOTDRV('O2 ','L1 ',NR2TRAN,MXVEC,
C    &                 WORK(KR2TRAN), WORK(KR2DOTS), WORK(KR2CONS),
C    &                 WORK(KEND0), LEND0 )
        CALL CC_DOTDRV('R2 ','X1B',NR2TRAN,MXVEC,
     &                 WORK(KR2TRAN), WORK(KR2DOTS), WORK(KR2CONS),
     &                 WORK(KEND0), LEND0 )

        TIMR2 = SECOND() - TIM1
        WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")')
     &       '>>> Time used for',
     &    NR2TRAN,' X1 x R2 dot products: ', SECOND() - TIM1
        CALL FLSHFO(LUPRI)

      ELSE
        TIMR2 = 0.0d0
      END IF

*---------------------------------------------------------------------*
* calculate ETA{O2} x R1 and Xksi{O2} x L1 dots products:
*---------------------------------------------------------------------*
      TIM1 = SECOND()

      CALL DZERO(WORK(KXCONS),MXXEDOTS)
      CALL DZERO(WORK(KECONS),MXXEDOTS)

      IOPT   = 5
      IORDER = 2
      CALL CC_XIETA( WORK(KXETRAN), NXETRAN, IOPT, IORDER, 'L0 ',
     &               'L1 ',WORK(KXDOTS),WORK(KXCONS),
     &               'R1 ',WORK(KEDOTS),WORK(KECONS),
     &               .FALSE.,MXVEC, WORK(KEND0), LEND0 )

      TIMXE = SECOND() - TIM1
      WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")') '>>> Time used for',
     & NXETRAN,' Xksi{O} and Eta{O2} vector calculat.: ', SECOND()-TIM1
      CALL FLSHFO(LUPRI)

*=====================================================================*
* loop over triples of operator labels and over frequencies and 
* collect all contributions to the final hyperpolarizabilities:
*=====================================================================*
      DO IOPER = 1, NQROPER

        IOPERA = IAQROP(IOPER)
        IOPERB = IBQROP(IOPER)
        IOPERC = ICQROP(IOPER)

        ISYMB  = ISYOPR(IOPERB)
        ISYMC  = ISYOPR(IOPERC)
        ISYMA  = ISYOPR(IOPERA)


      IF (MULD2H(ISYMB,ISYMC) .EQ. ISYMA) THEN

        LORXA  = LAQLRX(IOPER)
        LORXB  = LBQLRX(IOPER)
        LORXC  = LCQLRX(IOPER)

        LPDBSA = LPDBSOP(IOPERA)
        LPDBSB = LPDBSOP(IOPERB)
        LPDBSC = LPDBSOP(IOPERC) 

      DO IFREQ = 1, NQRFREQ

      DO ISIGN = +1, -1, -2
        SIGN = DBLE(ISIGN)
    
        ITAMPA = IR1TAMP(LBLOPR(IOPERA),LORXA,SIGN*AQRFR(IFREQ),ISYML)
        ITAMPB = IR1TAMP(LBLOPR(IOPERB),LORXB,SIGN*BQRFR(IFREQ),ISYML)
        ITAMPC = IR1TAMP(LBLOPR(IOPERC),LORXC,SIGN*CQRFR(IFREQ),ISYML)

        IZETVA = IL1ZETA(LBLOPR(IOPERA),LORXA,SIGN*AQRFR(IFREQ),ISYML)
        IZETVB = IL1ZETA(LBLOPR(IOPERB),LORXB,SIGN*BQRFR(IFREQ),ISYML)
        IZETVC = IL1ZETA(LBLOPR(IOPERC),LORXC,SIGN*CQRFR(IFREQ),ISYML)

        IKAPA  = 0
        IKAPB  = 0
        IKAPC  = 0
        IF(LORXA)IKAPA=IR1KAPPA(LBLOPR(IOPERA),SIGN*AQRFR(IFREQ),ISYML)
        IF(LORXB)IKAPB=IR1KAPPA(LBLOPR(IOPERB),SIGN*BQRFR(IFREQ),ISYML)
        IF(LORXC)IKAPC=IR1KAPPA(LBLOPR(IOPERC),SIGN*CQRFR(IFREQ),ISYML)

        IF ( USE_R2 ) THEN
          IR2AB=IR2TAMP(LBLOPR(IOPERA),.FALSE.,SIGN*AQRFR(IFREQ),ISYML1,
     &                  LBLOPR(IOPERB),.FALSE.,SIGN*BQRFR(IFREQ),ISYML2)
          IR2AC=IR2TAMP(LBLOPR(IOPERA),.FALSE.,SIGN*AQRFR(IFREQ),ISYML1,
     &                  LBLOPR(IOPERC),.FALSE.,SIGN*CQRFR(IFREQ),ISYML2)
          IR2BC=IR2TAMP(LBLOPR(IOPERB),.FALSE.,SIGN*BQRFR(IFREQ),ISYML1,
     &                  LBLOPR(IOPERC),.FALSE.,SIGN*CQRFR(IFREQ),ISYML2)
        END IF

*---------------------------------------------------------------------*

        CALL CCQR_SETG(WORK(KGTRAN),WORK(KGDOTS),NGTRAN,MXVEC,
     &                 0,ITAMPA,ITAMPB,ITAMPC,ITRAN,IVEC)
        GCON = WORK(KGCONS-1 + (ITRAN-1)*MXVEC + IVEC)

*---------------------------------------------------------------------*

        CALL CC_SETFA(WORK(KFATRAN),WORK(KFADOTS),NFATRAN,MXVEC,
     &                0,IOPERA,IKAPA,ITAMPB,ITAMPC,ITRAN,IVEC)
        FACON1 = WORK(KFACONS-1 + (ITRAN-1)*MXVEC + IVEC)

        CALL CC_SETFA(WORK(KFATRAN),WORK(KFADOTS),NFATRAN,MXVEC,
     &                0,IOPERB,IKAPB,ITAMPC,ITAMPA,ITRAN,IVEC)
        FACON2 = WORK(KFACONS-1 + (ITRAN-1)*MXVEC + IVEC)

        CALL CC_SETFA(WORK(KFATRAN),WORK(KFADOTS),NFATRAN,MXVEC,
     &                0,IOPERC,IKAPC,ITAMPA,ITAMPB,ITRAN,IVEC)
        FACON3 = WORK(KFACONS-1 + (ITRAN-1)*MXVEC + IVEC)

C--------------------------------------------------------------------*
C     K-Matrix contribution:
C--------------------------------------------------------------------*
C
        IF (CCSLV) THEN
          CALL CC_SETB11(WORK(KKTRAN),WORK(KKDOTS),NKTRAN,MXVEC,
     &                   ITAMPA,IZETVB,IZETVC,ITRAN,IVEC)

          SLVKCON1 = WORK(KKCONS-1 + (ITRAN-1)*MXVEC + IVEC)

          CALL CC_SETB11(WORK(KKTRAN),WORK(KKDOTS),NKTRAN,MXVEC,
     &                   ITAMPB,IZETVC,IZETVA,ITRAN,IVEC)
          SLVKCON2 = WORK(KKCONS-1 + (ITRAN-1)*MXVEC + IVEC)

          CALL CC_SETB11(WORK(KKTRAN),WORK(KKDOTS),NKTRAN,MXVEC,
     &                   ITAMPC,IZETVA,IZETVB,ITRAN,IVEC)
          SLVKCON3 = WORK(KKCONS-1 + (ITRAN-1)*MXVEC + IVEC)
        ELSE
           SLVKCON1 = 0.0d0
           SLVKCON2 = 0.0d0
           SLVKCON3 = 0.0d0
        END IF
C
C---------------------------------------------------------------------*
C
        IF ( USEBTRAN .AND. (.NOT.USE_R2) ) THEN
          CALL CC_SETB11(WORK(KBTRAN),WORK(KBDOTS),NBTRAN,MXVEC,
     &                   IZETVA,ITAMPB,ITAMPC,ITRAN,IVEC)
          BCON1 = WORK(KBCONS-1 + (ITRAN-1)*MXVEC + IVEC)

          CALL CC_SETB11(WORK(KBTRAN),WORK(KBDOTS),NBTRAN,MXVEC,
     &                   IZETVB,ITAMPC,ITAMPA,ITRAN,IVEC)
          BCON2 = WORK(KBCONS-1 + (ITRAN-1)*MXVEC + IVEC)

          CALL CC_SETB11(WORK(KBTRAN),WORK(KBDOTS),NBTRAN,MXVEC,
     &                   IZETVC,ITAMPA,ITAMPB,ITRAN,IVEC)
          BCON3 = WORK(KBCONS-1 + (ITRAN-1)*MXVEC + IVEC)

          FCON1  = 0.0d0 
          FCON2  = 0.0d0
          FCON3  = 0.0d0
          R2CON1 = 0.0d0 
          R2CON2 = 0.0d0
          R2CON3 = 0.0d0

        ELSE IF ( (.NOT.USEBTRAN) .AND. (.NOT.USE_R2) ) THEN

          CALL CCQR_SETF(WORK(KFTRAN),WORK(KFDOTS),NFTRAN,MXVEC,
     &                   IZETVA,ITAMPB,ITAMPC,ITRAN,IVEC)
          FCON1 = WORK(KFCONS-1 + (ITRAN-1)*MXVEC + IVEC)

          CALL CCQR_SETF(WORK(KFTRAN),WORK(KFDOTS),NFTRAN,MXVEC,
     &                   IZETVB,ITAMPC,ITAMPA,ITRAN,IVEC)
          FCON2 = WORK(KFCONS-1 + (ITRAN-1)*MXVEC + IVEC)

          CALL CCQR_SETF(WORK(KFTRAN),WORK(KFDOTS),NFTRAN,MXVEC,
     &                   IZETVC,ITAMPA,ITAMPB,ITRAN,IVEC)
          FCON3 = WORK(KFCONS-1 + (ITRAN-1)*MXVEC + IVEC)

          BCON1  = 0.0d0 
          BCON2  = 0.0d0
          BCON3  = 0.0d0
          R2CON1 = 0.0d0 
          R2CON2 = 0.0d0
          R2CON3 = 0.0d0

        ELSE IF ( USE_R2 ) THEN
 
          CALL CC_SETDOT(WORK(KR2TRAN),WORK(KR2DOTS),NR2TRAN,MXVEC,
     &                   IR2AB,IZETVC,ITRAN,IVEC)
          R2CON1 = WORK(KR2CONS-1 + (ITRAN-1)*MXVEC + IVEC)

          CALL CC_SETDOT(WORK(KR2TRAN),WORK(KR2DOTS),NR2TRAN,MXVEC,
     &                   IR2AC,IZETVB,ITRAN,IVEC)
          R2CON2 = WORK(KR2CONS-1 + (ITRAN-1)*MXVEC + IVEC)
          
          CALL CC_SETDOT(WORK(KR2TRAN),WORK(KR2DOTS),NR2TRAN,MXVEC,
     &                   IR2BC,IZETVA,ITRAN,IVEC)
          R2CON3 = WORK(KR2CONS-1 + (ITRAN-1)*MXVEC + IVEC)

          BCON1 = 0.0d0 
          BCON2 = 0.0d0
          BCON3 = 0.0d0
          FCON1 = 0.0d0 
          FCON2 = 0.0d0
          FCON3 = 0.0d0
        ELSE 
          CALL QUIT('Error in CC_HYPPOL.')
        END IF

*---------------------------------------------------------------------*

        IF ( .NOT. USE_R2 ) THEN
         IF (USE_AAMAT) THEN
          IF (LOCDBG) WRITE(LUPRI,*) 'Use A{O} matrix transformation.'
          IF (IKAPA.NE.0 .OR. IKAPB.NE.0 .OR. IKAPC.NE.0) THEN
            CALL QUIT('A{O} matrix transformation not yet relaxed.')
          END IF
 
          CALL CC_SETAA(WORK(KAATRAN),WORK(KAADOTS),MXTRAN,MXVEC,
     &                  IZETVA,IOPERB,ITAMPC,ITRAN,IVEC)
          EACON1 = WORK(KAACONS-1 + (ITRAN-1)*MXVEC + IVEC)

          CALL CC_SETAA(WORK(KAATRAN),WORK(KAADOTS),MXTRAN,MXVEC,
     &                  IZETVA,IOPERC,ITAMPB,ITRAN,IVEC)
          EACON2 = WORK(KAACONS-1 + (ITRAN-1)*MXVEC + IVEC)

          CALL CC_SETAA(WORK(KAATRAN),WORK(KAADOTS),MXTRAN,MXVEC,
     &                  IZETVB,IOPERC,ITAMPA,ITRAN,IVEC)
          EACON3 = WORK(KAACONS-1 + (ITRAN-1)*MXVEC + IVEC)

          CALL CC_SETAA(WORK(KAATRAN),WORK(KAADOTS),MXTRAN,MXVEC,
     &                  IZETVB,IOPERA,ITAMPC,ITRAN,IVEC)
          EACON4 = WORK(KAACONS-1 + (ITRAN-1)*MXVEC + IVEC)

          CALL CC_SETAA(WORK(KAATRAN),WORK(KAADOTS),MXTRAN,MXVEC,
     &                  IZETVC,IOPERA,ITAMPB,ITRAN,IVEC)
          EACON5 = WORK(KAACONS-1 + (ITRAN-1)*MXVEC + IVEC)

          CALL CC_SETAA(WORK(KAATRAN),WORK(KAADOTS),MXTRAN,MXVEC,
     &                  IZETVC,IOPERB,ITAMPA,ITRAN,IVEC)
          EACON6 = WORK(KAACONS-1 + (ITRAN-1)*MXVEC + IVEC)
         ELSE
          CALL CC_SETXE('Eta',WORK(KEATRAN),WORK(KEADOTS),NEATRAN,MXVEC,
     &                  IZETVA,IOPERB,IKAPB,0,0,0,ITAMPC,ITRAN,IVEC)
          EACON1 = WORK(KEACONS-1 + (ITRAN-1)*MXVEC + IVEC)
         
          CALL CC_SETXE('Eta',WORK(KEATRAN),WORK(KEADOTS),NEATRAN,MXVEC,
     &                  IZETVA,IOPERC,IKAPC,0,0,0,ITAMPB,ITRAN,IVEC)
          EACON2 = WORK(KEACONS-1 + (ITRAN-1)*MXVEC + IVEC)
         
          CALL CC_SETXE('Eta',WORK(KEATRAN),WORK(KEADOTS),NEATRAN,MXVEC,
     &                  IZETVB,IOPERC,IKAPC,0,0,0,ITAMPA,ITRAN,IVEC)
          EACON3 = WORK(KEACONS-1 + (ITRAN-1)*MXVEC + IVEC)
         
          CALL CC_SETXE('Eta',WORK(KEATRAN),WORK(KEADOTS),NEATRAN,MXVEC,
     &                  IZETVB,IOPERA,IKAPA,0,0,0,ITAMPC,ITRAN,IVEC)
          EACON4 = WORK(KEACONS-1 + (ITRAN-1)*MXVEC + IVEC)
         
          CALL CC_SETXE('Eta',WORK(KEATRAN),WORK(KEADOTS),NEATRAN,MXVEC,
     &                  IZETVC,IOPERA,IKAPA,0,0,0,ITAMPB,ITRAN,IVEC)
          EACON5 = WORK(KEACONS-1 + (ITRAN-1)*MXVEC + IVEC)
         
          CALL CC_SETXE('Eta',WORK(KEATRAN),WORK(KEADOTS),NEATRAN,MXVEC,
     &                  IZETVC,IOPERB,IKAPB,0,0,0,ITAMPA,ITRAN,IVEC)
          EACON6 = WORK(KEACONS-1 + (ITRAN-1)*MXVEC + IVEC)
         END IF
        ELSE

         EACON1 = 0.0d0
         EACON2 = 0.0d0
         EACON3 = 0.0d0
         EACON4 = 0.0d0
         EACON5 = 0.0d0
         EACON6 = 0.0d0

        END IF


*---------------------------------------------------------------------*
        IF (LPDBSA.OR.LORXA .OR. LPDBSB.OR.LORXB) THEN
          CALL CC_FIND_SO_OP(LBLOPR(IOPERA),LBLOPR(IOPERB),
     &                       LABSOP,ISYSOP,ISGNSOP,INUM,WORK,LWORK)
          IOPSOP = IROPER(LABSOP,ISYSOP)

          CALL CC_SETXE('Xi ',WORK(KXETRAN),WORK(KXDOTS),MXTRAN,MXVEC,
     &                  0,IOPSOP,IKAPA,IKAPB,0,0,IZETVC,ITRAN,IVEC)
          X2CON1 = WORK(KXCONS-1 + (ITRAN-1)*MXVEC + IVEC)

          CALL CC_SETXE('Eta',WORK(KXETRAN),WORK(KEDOTS),MXTRAN,MXVEC,
     &                  0,IOPSOP,IKAPA,IKAPB,0,0,ITAMPC,ITRAN,IVEC)
          E2CON1 = WORK(KECONS-1 + (ITRAN-1)*MXVEC + IVEC)
        ELSE
          X2CON1 = 0.0d0
          E2CON1 = 0.0d0
        END IF                     

        IF (LPDBSA.OR.LORXA .OR. LPDBSC.OR.LORXC) THEN
          CALL CC_FIND_SO_OP(LBLOPR(IOPERA),LBLOPR(IOPERC),
     &                       LABSOP,ISYSOP,ISGNSOP,INUM,WORK,LWORK)
          IOPSOP = IROPER(LABSOP,ISYSOP)

          CALL CC_SETXE('Xi ',WORK(KXETRAN),WORK(KXDOTS),MXTRAN,MXVEC,
     &                  0,IOPSOP,IKAPA,IKAPC,0,0,IZETVB,ITRAN,IVEC)
          X2CON2 = WORK(KXCONS-1 + (ITRAN-1)*MXVEC + IVEC)

          CALL CC_SETXE('Eta',WORK(KXETRAN),WORK(KEDOTS),MXTRAN,MXVEC,
     &                  0,IOPSOP,IKAPA,IKAPC,0,0,ITAMPB,ITRAN,IVEC)
          E2CON2 = WORK(KECONS-1 + (ITRAN-1)*MXVEC + IVEC)
        ELSE
          X2CON2 = 0.0d0
          E2CON2 = 0.0d0
        END IF                         

        IF (LPDBSC.OR.LORXC .OR. LPDBSB.OR.LORXB) THEN
          CALL CC_FIND_SO_OP(LBLOPR(IOPERC),LBLOPR(IOPERB),
     &                       LABSOP,ISYSOP,ISGNSOP,INUM,WORK,LWORK)
          IOPSOP = IROPER(LABSOP,ISYSOP)

          CALL CC_SETXE('Xi ',WORK(KXETRAN),WORK(KXDOTS),MXTRAN,MXVEC,
     &                  0,IOPSOP,IKAPC,IKAPB,0,0,IZETVA,ITRAN,IVEC)
          X2CON3 = WORK(KXCONS-1 + (ITRAN-1)*MXVEC + IVEC)

          CALL CC_SETXE('Eta',WORK(KXETRAN),WORK(KEDOTS),MXTRAN,MXVEC,
     &                  0,IOPSOP,IKAPC,IKAPB,0,0,ITAMPA,ITRAN,IVEC)
          E2CON3 = WORK(KECONS-1 + (ITRAN-1)*MXVEC + IVEC)
        ELSE
          X2CON3 = 0.0d0
          E2CON3 = 0.0d0
        END IF                
*---------------------------------------------------------------------*

        IDX = NQRFREQ*(IOPER-1)+IFREQ
        
        IF (ISIGN.EQ.-1) IDX = IDX + NBHYPPOL

        WORK(KHYPPOL-1+IDX) = WORK(KHYPPOL-1+IDX) + 
     &         GCON   + 
     &         FACON1 + FACON2 + FACON3 +
     &         FCON1  + FCON2  + FCON3  +
     &         BCON1  + BCON2  + BCON3  +
     &         SLVKCON1  + SLVKCON2  + SLVKCON3  +
     &         R2CON1 + R2CON2 + R2CON3 +
     &         EACON1 + EACON2 + EACON3 + EACON4 + EACON5 + EACON6 +
     &         E2CON1 + E2CON2 + E2CON3 + X2CON1 + X2CON2 + X2CON3

*---------------------------------------------------------------------*

        IF (LOCDBG .OR. IPRQHYP.GE.15) THEN
         WRITE(LUPRI,*) 'ITAMPA,ITAMPB,ITAMPC:',ITAMPA,ITAMPB,ITAMPC
         WRITE(LUPRI,*) 'IZETVA,IZETVB,IZETVC:',IZETVA,IZETVB,IZETVC
         WRITE(LUPRI,*) 'IOPERA,IOPERB,IOPERC:',IOPERA,IOPERB,IOPERC
         WRITE(LUPRI,*) 'IKAPA, IKAPB, IKAPC :',IKAPA, IKAPB, IKAPC
         WRITE(LUPRI,*) 'GCON :',GCON
         WRITE(LUPRI,*) 'FACON:',FACON1,FACON2,FACON3
         WRITE(LUPRI,*) 'FCON :',FCON1,FCON2,FCON3
         WRITE(LUPRI,*) 'BCON :',BCON1,BCON2,BCON3
         WRITE(LUPRI,*) 'SLVKCON :',SLVKCON1,SLVKCON2,SLVKCON3
         WRITE(LUPRI,*) 'EACON:',EACON1,EACON2,EACON3,
     &                           EACON4,EACON5,EACON6
         WRITE(LUPRI,*) 'E2CON:',E2CON1,E2CON2,E2CON3
         WRITE(LUPRI,*) 'X2CON:',X2CON1,X2CON2,X2CON3
         WRITE(LUPRI,*) 'R2CON:',R2CON1,R2CON2,R2CON3
         WRITE(LUPRI,*) 'HYPPOL:',WORK(KHYPPOL-1+IDX)
        END IF

      END DO 
      END DO 
      END IF
      END DO 

*---------------------------------------------------------------------*
* print timing:
*---------------------------------------------------------------------*
      WRITE (LUPRI,'(/A,I4,A,F12.2," seconds.")') '>>> Total time for',
     &  NBHYPPOL,' quadratic response func.:', SECOND() - TIM0

*---------------------------------------------------------------------*
* print frequency-dependent hyperpolarizabilities:
*---------------------------------------------------------------------*

      CALL  CCQHYPPRT(WORK(KHYPPOL))

      CALL FLSHFO(LUPRI)

      IF (NQRDISP.EQ.0) THEN
        RETURN
      END IF
  
*---------------------------------------------------------------------*
* allocate & initialize work space for expansion coefficients
*---------------------------------------------------------------------*
      NBEXPCOF = NQROPER * NQRDISP

      MXTRAN  = 2 * NLRCLBL * NLRCLBL
      MXVEC   = NLRCLBL

      MXGTRAN  = 4 * MXTRAN
      MXFTRAN  = 3 * MXTRAN
      MXBTRAN  = 3 * MXTRAN
      MXFATRAN = 5 * MXTRAN
      MXEATRAN = 3 * MXTRAN
      MXR2TRAN = 1 * MXTRAN

      MXGDOTS  = MXVEC * MXTRAN
      MXFDOTS  = MXVEC * MXTRAN
      MXBDOTS  = MXVEC * MXTRAN
      MXFADOTS = MXVEC * MXTRAN
      MXEADOTS = MXVEC * MXTRAN
      MXR2DOTS = MXVEC * MXTRAN

      KEXPCOF = 1
      KGTRAN  = KEXPCOF + NBEXPCOF
      KGDOTS  = KGTRAN  + MXGTRAN
      KGCONS  = KGDOTS  + MXGDOTS
      KFTRAN  = KGCONS  + MXGDOTS
      KFDOTS  = KFTRAN  + MXFTRAN
      KFCONS  = KFDOTS  + MXFDOTS
      KBTRAN  = KFCONS  + MXFDOTS
      KBDOTS  = KBTRAN  + MXBTRAN
      KBCONS  = KBDOTS  + MXBDOTS
      KFATRAN = KBCONS  + MXBDOTS
      KFADOTS = KFATRAN + MXFATRAN
      KFACONS = KFADOTS + MXFADOTS
      KEATRAN = KFACONS + MXFADOTS
      KEADOTS = KEATRAN + MXEATRAN
      KEACONS = KEADOTS + MXEADOTS
      KR2TRAN = KEACONS + MXEADOTS
      KR2DOTS = KR2TRAN + MXR2TRAN
      KR2CONS = KR2DOTS + MXR2DOTS
      KEND0   = KR2CONS + MXR2DOTS
      LEND0   = LWORK - KEND0

      IF (LEND0 .LT. 0) THEN
        CALL QUIT('Insufficient memory in CC_HYPPOL. (2)')
      END IF

* initialze hyperpolarizabilities and the lists of contributions:
      CALL DZERO(WORK(KEXPCOF),NBEXPCOF)
      CALL DZERO(WORK(KGTRAN), MXGTRAN +2*MXGDOTS)
      CALL DZERO(WORK(KFTRAN), MXFTRAN +2*MXFDOTS)
      CALL DZERO(WORK(KBTRAN), MXBTRAN +2*MXBDOTS)
      CALL DZERO(WORK(KFATRAN),MXFATRAN+2*MXFADOTS)
      CALL DZERO(WORK(KEATRAN),MXEATRAN+2*MXEADOTS)
      CALL DZERO(WORK(KR2TRAN),MXR2TRAN+2*MXR2DOTS)

* initialize timing:
      TIM2 = SECOND()

*---------------------------------------------------------------------*
* set up lists for G, F/B and F{A} transformations and ETA{O} vectors:
*---------------------------------------------------------------------*
      LADD = .FALSE.

      CALL CCQR_DISP_SETUP(MXTRAN, MXVEC,
     &       WORK(KGTRAN), WORK(KGDOTS), WORK(KGCONS), NGTRAN,
     &       WORK(KFTRAN), WORK(KFDOTS), WORK(KFCONS), NFTRAN,
     &       WORK(KBTRAN), WORK(KBDOTS), WORK(KBCONS), NBTRAN,
     &       WORK(KFATRAN),WORK(KFADOTS),WORK(KFACONS),NFATRAN,
     &       WORK(KEATRAN),WORK(KEADOTS),WORK(KEACONS),NEATRAN,
     &       WORK(KR2TRAN),WORK(KR2DOTS),WORK(KR2CONS),NR2TRAN,
     &       WORK(KEXPCOF),NBEXPCOF, LADD                       )

*---------------------------------------------------------------------*
* calculate G matrix contributions:
*---------------------------------------------------------------------*
      TIM1 = SECOND()

      IOPT = 5
      CALL CC_GMATRIX('L0 ','RC ','RC ','RC ',NGTRAN, MXVEC,
     &              WORK(KGTRAN),WORK(KGDOTS),WORK(KGCONS),
     &              WORK(KEND0), LEND0, IOPT )

      TIMG = TIMG + SECOND() - TIM1

      WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")') '>>> Time used for',
     &  NGTRAN,' G matrix transformations:', SECOND() - TIM1
      CALL FLSHFO(LUPRI)

*---------------------------------------------------------------------*
* calculate B matrix contributions:
*---------------------------------------------------------------------*
      IF ( USEBTRAN .AND. (.NOT.USE_R2) ) THEN

        TIM1 = SECOND()

        IOPT = 5
        CALL CC_BMATRIX(WORK(KBTRAN),NBTRAN,'RC ','RC ',IOPT,'LC ',
     &               WORK(KBDOTS),WORK(KBCONS),MXVEC,
     &               .FALSE.,WORK(KEND0),LEND0)

        TIMB = TIMB + SECOND() - TIM1

        WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")') '>>> Time used for',
     &    NBTRAN,' B matrix transformations:', SECOND() - TIM1
        CALL FLSHFO(LUPRI)

      END IF

*---------------------------------------------------------------------*
* calculate F matrix contributions:
*---------------------------------------------------------------------*
      IF ( (.NOT.USEBTRAN) .AND. (.NOT.USE_R2) ) THEN

        TIM1 = SECOND()

        IOPT = 5
        CALL CC_FMATRIX(WORK(KFTRAN),NFTRAN,'LC ','RC ',IOPT,'RC ',
     &                  WORK(KFDOTS),WORK(KFCONS),MXVEC,
     &                  WORK(KEND0), LEND0)

        TIMF = TIMF + SECOND() - TIM1

        WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")') '>>> Time used for',
     &    NFTRAN,' F matrix transformations:', SECOND() - TIM1
        CALL FLSHFO(LUPRI)

      END IF

*---------------------------------------------------------------------*
* calculate F{O} matrix contributions:
*---------------------------------------------------------------------*
      TIM1 = SECOND()

      CALL CCQR_FADRV('L0 ','o1 ','RC ','RC ',NFATRAN, MXVEC,
     &                 WORK(KFATRAN),WORK(KFADOTS),WORK(KFACONS),
     &                 WORK(KEND0), LEND0, 'DOTP' )

      TIMFA = TIMFA + SECOND() - TIM1

      WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")') '>>> Time used for',
     &   NFATRAN,' F{O} matrix transform.:  ', SECOND() - TIM1
      CALL FLSHFO(LUPRI)

*---------------------------------------------------------------------*
* calculate ETA{O} vector contributions:
*---------------------------------------------------------------------*
      IF ( .NOT. USE_R2 ) THEN
        TIM1 = SECOND()

        CALL CCQR_EADRV('LC ','o1 ','RC ',NEATRAN, MXVEC,
     &                   WORK(KEATRAN),WORK(KEADOTS),WORK(KEACONS),
     &                   WORK(KEND0), LEND0, 'DOTP' )

        TIMEA = TIMEA + SECOND() - TIM1

        WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")') '>>> Time used for',
     &      NEATRAN,' ETA{O} vector calculat.: ', SECOND() - TIM1
        CALL FLSHFO(LUPRI)

      ELSE
        TIMEA = 0.0d0
      END IF

*---------------------------------------------------------------------*
* calculate X1 x R2 dots products:
*---------------------------------------------------------------------*
      IF ( USE_R2 ) THEN
        TIM1 = SECOND()

        CALL CC_DOTDRV('CR2','XC ',NR2TRAN,MXVEC,
     &                 WORK(KR2TRAN), WORK(KR2DOTS), WORK(KR2CONS),
     &                 WORK(KEND0), LEND0 )

        TIMR2 = SECOND() - TIM1
        WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")') '>>> Time used for',
     &    NR2TRAN,' XC x CR2 dot products: ', SECOND() - TIM1
        CALL FLSHFO(LUPRI)

      ELSE
        TIMR2 = 0.0d0
      END IF

*---------------------------------------------------------------------*
* collect contributions and add them to hyperpolarizabilities
*---------------------------------------------------------------------*
      LADD = .TRUE.

      CALL CCQR_DISP_SETUP(MXTRAN, MXVEC,
     &       WORK(KGTRAN), WORK(KGDOTS), WORK(KGCONS), NGTRAN,
     &       WORK(KFTRAN), WORK(KFDOTS), WORK(KFCONS), NFTRAN,
     &       WORK(KBTRAN), WORK(KBDOTS), WORK(KBCONS), NBTRAN,
     &       WORK(KFATRAN),WORK(KFADOTS),WORK(KFACONS),NFATRAN,
     &       WORK(KEATRAN),WORK(KEADOTS),WORK(KEACONS),NEATRAN,
     &       WORK(KR2TRAN),WORK(KR2DOTS),WORK(KR2CONS),NR2TRAN,
     &       WORK(KEXPCOF),NBEXPCOF, LADD                       )

*---------------------------------------------------------------------*
* print timing:
*---------------------------------------------------------------------*
      WRITE (LUPRI,'(/A,I4,A,F12.2," seconds.")') '>>> Total time for',
     &  NBEXPCOF,' expansion coeffients:    ', SECOND() - TIM2

*---------------------------------------------------------------------*
* print expansion coefficients for first hyperpolarizabilities:
*---------------------------------------------------------------------*

      IF (IPRQHYP .GE. 15) THEN
        CALL  CCQEXPPRT(WORK(KEXPCOF))
      END IF

      CALL FLSHFO(LUPRI)

*---------------------------------------------------------------------*
* calculate from the d_ABC expansion coefficients the D_ABC disp. cof.
*---------------------------------------------------------------------*
      NBINOM = (NQRDSPE+2)*(NQRDSPE+1)/2

      KBINOM = KEND0 
      KDSPCF = KBINOM + NBINOM
      KSHGCF = KDSPCF + NBINOM*NQROPER
      KORCOF = KSHGCF + (NQRDSPE+1)*NQROPER
      KEOPEC = KORCOF + (NQRDSPE+1)*NQROPER
      KEND1  = KEOPEC + (NQRDSPE+1)*NQROPER
      LEND1  = LWORK - KEND1

      IF (BETA_AVERAGE) THEN
        KDSPCFA = KEND1 
        KSHGCFA = KDSPCFA + 4*NBINOM
        KORCOFA = KSHGCFA + 4*(NQRDSPE+1)
        KEOPECA = KORCOFA + 4*(NQRDSPE+1)
        KABCDE  = KEOPECA + 4*(NQRDSPE+1)
        KEND1   = KABCDE  + 12*2
        LEND1   = LWORK - KEND1
      END IF

      IF (LEND1 .LT. 0) THEN
        CALL QUIT('Insufficient memory in CC_HYPPOL. (3)')
      END IF

      CALL DZERO(WORK(KBINOM),KEND1-KBINOM)

      CALL CCQRDISP(WORK(KDSPCF), WORK(KDSPCFA),WORK(KABCDE),
     &              WORK(KSHGCF), WORK(KORCOF), WORK(KEOPEC),
     &              WORK(KSHGCFA),WORK(KORCOFA),WORK(KEOPECA),
     &              WORK(KEXPCOF),WORK(KBINOM))

*---------------------------------------------------------------------*
* print dispersion coefficients for first hyperpolarizabilities:
*---------------------------------------------------------------------*

      CALL CCQDISPRT(WORK(KDSPCF), WORK(KDSPCFA),WORK(KABCDE),
     &               WORK(KSHGCF), WORK(KORCOF), WORK(KEOPEC),
     &               WORK(KSHGCFA),WORK(KORCOFA),WORK(KEOPECA))

*---------------------------------------------------------------------*
* print timing & return:
*---------------------------------------------------------------------*
      WRITE (LUPRI,'(/A,A,F12.2," seconds.")') '>>> Total time used ',
     &  'in quadratic response section:', SECOND() - TIM0

      CALL FLSHFO(LUPRI)

      RETURN
      END

*=====================================================================*
*              END OF SUBROUTINE CC_HYPPOL                            *
*=====================================================================*
c /* deck ccqrdisp */
*=====================================================================*
      SUBROUTINE CCQRDISP(DISPCF,AVEDISP,ABCDE,
     &                    SHGCF, ORCOF,EOPECF,
     &                    AVESHG,AVEOR,AVEEOPE,
     &                    EXPCOF,BINOM)
*---------------------------------------------------------------------*
*
*   Purpose: calculate from the expansion coefficients defined by
*            beta_ABC = sum_{klm} w_A^k w_B^l w_C^m d_ABC(klm) the
*            physical more relevant dispersion coefficients defined
*            by beta_ABC = sum_{mn} w_B^m w_C^n D_ABC(mn)
*
*
*   Written by Christof Haettig in October 1997.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccqrinf.h>
#include <ccroper.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER KA, KB, KC, KCP, KD, KDP, KE, KEP
      PARAMETER (KA=1,KB=2,KC=3,KCP=4,KD=5,KDP=6,KE=7,KEP=8)
      INTEGER KBP, KDP2, KEP2, KEP3
      PARAMETER (KBP=9,KDP2=10,KEP2=11,KEP3=12)

      INTEGER IDISP, ICOMP, ISAMA, ISAMB, ISAMC, ISAPROP, MN0, MNS


#if defined (SYS_CRAY)
      REAL DISPCF((NQRDSPE+2)*(NQRDSPE+1)/2,NQROPER)
      REAL AVEDISP((NQRDSPE+2)*(NQRDSPE+1)/2,4)
      REAL ABCDE(12,2)
      REAL SHGCF(NQRDSPE+1,NQROPER)
      REAL ORCOF(NQRDSPE+1,NQROPER)
      REAL EOPECF(NQRDSPE+1,NQROPER)
      REAL AVESHG(NQRDSPE+1,4)
      REAL AVEOR(NQRDSPE+1,4)
      REAL AVEEOPE(NQRDSPE+1,4)
      REAL EXPCOF(NQRDISP,NQROPER)
      REAL BINOM((NQRDSPE+2)*(NQRDSPE+1)/2)
      REAL SIGNP, SIGNPQ, SIGNN, ERROR
      REAL BETA0,A0,A1,A2,B0,B1, B2, B3, B4, C0, C1, C2, C3
      REAL CP0,CP1,CP2,D0,D1,D2,D3, DP0, DP1, DP2, DP3, DP4
      REAL E0,E1,E2,E3, EP0, EP1, EP2, EP3, EP4, EP5, EP6
      REAL BP0,BP1,CP3,DPP0, DPP1, DPP2, DPP3
      REAL EPP0, EPP1, EPP2, EPPP0, EPPP1, EPPP2, EPPP3
      REAL BORTFAC(7), BM23FAC(7)
#else
      DOUBLE PRECISION DISPCF((NQRDSPE+2)*(NQRDSPE+1)/2,NQROPER)
      DOUBLE PRECISION AVEDISP((NQRDSPE+2)*(NQRDSPE+1)/2,4)
      DOUBLE PRECISION ABCDE(12,2)
      DOUBLE PRECISION SHGCF(NQRDSPE+1,NQROPER)
      DOUBLE PRECISION ORCOF(NQRDSPE+1,NQROPER)
      DOUBLE PRECISION EOPECF(NQRDSPE+1,NQROPER)
      DOUBLE PRECISION AVESHG(NQRDSPE+1,4)
      DOUBLE PRECISION AVEOR(NQRDSPE+1,4)
      DOUBLE PRECISION AVEEOPE(NQRDSPE+1,4)
      DOUBLE PRECISION EXPCOF(NQRDISP,NQROPER)
      DOUBLE PRECISION BINOM((NQRDSPE+2)*(NQRDSPE+1)/2)
      DOUBLE PRECISION SIGNP, SIGNPQ, SIGNN, ERROR
      DOUBLE PRECISION BETA0,A0,A1,A2,B0,B1, B2, B3, B4, C0, C1, C2, C3
      DOUBLE PRECISION CP0,CP1,CP2,D0,D1,D2,D3, DP0, DP1, DP2, DP3, DP4
      DOUBLE PRECISION E0,E1,E2,E3, EP0, EP1, EP2, EP3, EP4, EP5, EP6
      DOUBLE PRECISION BP0,BP1,CP3,DPP0, DPP1, DPP2, DPP3
      DOUBLE PRECISION EPP0, EPP1, EPP2, EPPP0, EPPP1, EPPP2, EPPP3
      DOUBLE PRECISION BORTFAC(7), BM23FAC(7), BKERFAC(7)
#endif
C
      DATA BORTFAC / 0.2d0, 0.4d0,-0.6d0,0.4d0, 0.4d0,-0.6d0,0.4d0 /
      DATA BM23FAC / 0.0d0, 0.5d0,-1.0d0,0.5d0, 0.5d0,-1.0d0,0.5d0 /
      DATA BKERFAC / 0.6d0,-0.3d0,1.2d0,-0.3d0,-0.3d0,1.2d0,-0.3d0 /
C
      INTEGER N, M, MN, P, Q, I, J, ITRI, IOPER, OFFEX, IEXPCF
      ITRI(I,J) = (MAX(I,J)+1)*(MAX(I,J)-2)/2 + I + J + 2
C

      DO M = 0, NQRDSPE
        BINOM(ITRI(M,0)) = 1.0d0
        DO N = 1, M-1
          BINOM(ITRI(M,N)) = 
     &       BINOM(ITRI(M-1,N-1)) + BINOM(ITRI(M-1,N))
        END DO
        BINOM(ITRI(M,M)) = 1.0d0
      END DO

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'binomial coefficients:'
        DO M = 0, NQRDSPE
          WRITE (LUPRI,'(10F8.2)') (BINOM(ITRI(M,N)),N=0,M)
        END DO
      END IF

      DO IOPER = 1, NQROPER
        ISAMA   = ISYMAT(IAQROP(IOPER))
        ISAMB   = ISYMAT(IBQROP(IOPER))
        ISAMC   = ISYMAT(ICQROP(IOPER))
        ISAPROP = ISAMA * ISAMB * ISAMC
        IF (ALLDSPCF) ISAPROP = 0

        MN0 = 0
        MNS = 2

        IF (ISAPROP.EQ.-1) MN0 = 1
        IF (ISAPROP.EQ. 0) MNS = 1

        IF (LOCDBG) THEN
          WRITE (LUPRI,*) 'LABELS:',LBLOPR(IAQROP(IOPER)),
     &             LBLOPR(IBQROP(IOPER)), LBLOPR(ICQROP(IOPER))
          WRITE (LUPRI,*) 'ISAPROP,MN0,MNS:',ISAPROP,MN0,MNS
        END IF

        OFFEX = 0
        IF (MN0.EQ.1) OFFEX = 1
        DO MN = MN0, NQRDSPE, MNS
          IF (LOCDBG) THEN
            WRITE (LUPRI,'(//A,62X,A,10X,A,10X,A)')
     &      '  M  N    ', '   SHGCF ','ORCOF','EOPECF,'
          END IF
          SHGCF(MN+1,IOPER)  = 0.0d0
          ORCOF(MN+1,IOPER)  = 0.0d0
          EOPECF(MN+1,IOPER) = 0.0d0
          SIGNN = 1.0d0
          DO N  = 0, MN
            IF (LOCDBG) THEN
              WRITE (LUPRI,'(//A,A)')
     &        '  M  N  P  Q    SIGN   BINOM',
     &        '    I J K      EXPCOF        DISPCF' 
            END IF
            M = MN - N
            DISPCF(ITRI(M+N,N),IOPER) = 0.0d0
            SIGNP = +1.0d0
            DO P = 0, M
              SIGNPQ = SIGNP
              DO Q = 0, N
                IEXPCF = OFFEX + ITRI(M+N-P-Q,N-Q)
                DISPCF(ITRI(M+N,N),IOPER) = 
     &             DISPCF(ITRI(M+N,N),IOPER) +
     &              SIGNPQ * BINOM(ITRI(P+Q,P)) * EXPCOF(IEXPCF,IOPER)
                IF (LOCDBG) THEN
                  WRITE (LUPRI,'(4I3,2F8.2,3X,3I2,4E16.8)')
     &             M,N,P,Q,SIGNPQ,BINOM(ITRI(P+Q,P)),P+Q,M-P,N-Q,
     &             EXPCOF(IEXPCF,IOPER),DISPCF(ITRI(M+N,N),IOPER)
                END IF
                SIGNPQ = -SIGNPQ
              END DO 
              SIGNP = -SIGNP
            END DO 
            SHGCF(MN+1,IOPER) = SHGCF(MN+1,IOPER) + 
     &                           DISPCF(ITRI(M+N,N),IOPER)
            ORCOF(MN+1,IOPER) = ORCOF(MN+1,IOPER) +
     &                   SIGNN * DISPCF(ITRI(M+N,N),IOPER)
            IF (N.EQ.0) EOPECF(MN+1,IOPER)=DISPCF(ITRI(M+N,N),IOPER)
            SIGNN = -SIGNN
            IF (LOCDBG) THEN
              WRITE (LUPRI,'(2I3,66X,3E16.8)') M,N,SHGCF(MN+1,IOPER),
     &          ORCOF(MN+1,IOPER),EOPECF(MN+1,IOPER)
            END IF
          END DO
          OFFEX = OFFEX + (MN+2)*(MN+1)/2
          IF (MNS.EQ.2) OFFEX = OFFEX + (MN+3)*(MN+2)/2
        END DO
      END DO

      IF (BETA_AVERAGE) THEN

* calculate averaged vector components parallel/orthogonal to z axis:
        DO IDISP = 1, (NQRDSPE+2)*(NQRDSPE+1)/2
          AVEDISP(IDISP,1) = 0.6d0*DISPCF(IDISP,1)
          AVEDISP(IDISP,2) = BORTFAC(1)*DISPCF(IDISP,1)
          AVEDISP(IDISP,3) = BM23FAC(1)*DISPCF(IDISP,1)
          AVEDISP(IDISP,4) = BKERFAC(1)*DISPCF(IDISP,1)
        END DO
        DO IDISP = 1, NQRDSPE+1
          AVESHG(IDISP,1)  = 0.6d0*SHGCF(IDISP,1)
          AVEOR(IDISP,1)   = 0.6d0*ORCOF(IDISP,1)
          AVEEOPE(IDISP,1) = 0.6d0*EOPECF(IDISP,1)
          AVESHG(IDISP,2)  = BORTFAC(1)*SHGCF(IDISP,1)
          AVEOR(IDISP,2)   = BORTFAC(1)*ORCOF(IDISP,1)
          AVEEOPE(IDISP,2) = BORTFAC(1)*EOPECF(IDISP,1)
          AVESHG(IDISP,3)  = BM23FAC(1)*SHGCF(IDISP,1)
          AVEOR(IDISP,3)   = BM23FAC(1)*ORCOF(IDISP,1)
          AVEEOPE(IDISP,3) = BM23FAC(1)*EOPECF(IDISP,1)
          AVESHG(IDISP,4)  = BKERFAC(1)*SHGCF(IDISP,1)
          AVEOR(IDISP,4)   = BKERFAC(1)*ORCOF(IDISP,1)
          AVEEOPE(IDISP,4) = BKERFAC(1)*EOPECF(IDISP,1)
        END DO
        IF (XY_DEGENERAT) THEN
         DO ICOMP = 2, 4
          DO IDISP = 1, (NQRDSPE+2)*(NQRDSPE+1)/2
           AVEDISP(IDISP,1)= AVEDISP(IDISP,1)+0.4d0*DISPCF(IDISP,ICOMP)
           AVEDISP(IDISP,2)= AVEDISP(IDISP,2) +
     &                        2.0d0*BORTFAC(ICOMP)*DISPCF(IDISP,ICOMP)
           AVEDISP(IDISP,3)= AVEDISP(IDISP,3) +
     &                        2.0d0*BM23FAC(ICOMP)*DISPCF(IDISP,ICOMP)
           AVEDISP(IDISP,4)= AVEDISP(IDISP,4) +
     &                        2.0d0*BKERFAC(ICOMP)*DISPCF(IDISP,ICOMP)
          END DO
          DO IDISP = 1, NQRDSPE+1
           AVESHG(IDISP,1) = AVESHG(IDISP,1) +0.4d0*SHGCF(IDISP,ICOMP)
           AVEOR(IDISP,1)  = AVEOR(IDISP,1)  +0.4d0*ORCOF(IDISP,ICOMP)
           AVEEOPE(IDISP,1)= AVEEOPE(IDISP,1)+0.4d0*EOPECF(IDISP,ICOMP)
           AVESHG(IDISP,2) = AVESHG(IDISP,2) +
     &                         2.0d0*BORTFAC(ICOMP)*SHGCF(IDISP,ICOMP)
           AVEOR(IDISP,2)  = AVEOR(IDISP,2)  +
     &                         2.0d0*BORTFAC(ICOMP)*ORCOF(IDISP,ICOMP)
           AVEEOPE(IDISP,2)= AVEEOPE(IDISP,2)+
     &                         2.0d0*BORTFAC(ICOMP)*EOPECF(IDISP,ICOMP)
           AVESHG(IDISP,3) = AVESHG(IDISP,3) +
     &                         2.0d0*BM23FAC(ICOMP)*SHGCF(IDISP,ICOMP)
           AVEOR(IDISP,3)  = AVEOR(IDISP,3)  +
     &                         2.0d0*BM23FAC(ICOMP)*ORCOF(IDISP,ICOMP)
           AVEEOPE(IDISP,3)= AVEEOPE(IDISP,3)+
     &                         2.0d0*BM23FAC(ICOMP)*EOPECF(IDISP,ICOMP)
           AVESHG(IDISP,4) = AVESHG(IDISP,4) +
     &                         2.0d0*BKERFAC(ICOMP)*SHGCF(IDISP,ICOMP)
           AVEOR(IDISP,4)  = AVEOR(IDISP,4)  +
     &                         2.0d0*BKERFAC(ICOMP)*ORCOF(IDISP,ICOMP)
           AVEEOPE(IDISP,4)= AVEEOPE(IDISP,4)+
     &                         2.0d0*BKERFAC(ICOMP)*EOPECF(IDISP,ICOMP)
          END DO
         END DO
        ELSE
         DO ICOMP = 2, 7
          DO IDISP = 1, (NQRDSPE+2)*(NQRDSPE+1)/2
           AVEDISP(IDISP,1)= AVEDISP(IDISP,1)+0.2d0*DISPCF(IDISP,ICOMP)
           AVEDISP(IDISP,2)= AVEDISP(IDISP,2) +
     &                              BORTFAC(ICOMP)*DISPCF(IDISP,ICOMP)
           AVEDISP(IDISP,3)= AVEDISP(IDISP,3) +
     &                              BM23FAC(ICOMP)*DISPCF(IDISP,ICOMP)
           AVEDISP(IDISP,4)= AVEDISP(IDISP,4) +
     &                              BKERFAC(ICOMP)*DISPCF(IDISP,ICOMP)
          END DO
          DO IDISP = 1, NQRDSPE+1
           AVESHG(IDISP,1) = AVESHG(IDISP,1) +0.2d0*SHGCF(IDISP,ICOMP)
           AVEOR(IDISP,1)  = AVEOR(IDISP,1)  +0.2d0*ORCOF(IDISP,ICOMP)
           AVEEOPE(IDISP,1)= AVEEOPE(IDISP,1)+0.2d0*EOPECF(IDISP,ICOMP)
           AVESHG(IDISP,2) = AVESHG(IDISP,2) +
     &                               BORTFAC(ICOMP)*SHGCF(IDISP,ICOMP)
           AVEOR(IDISP,2)  = AVEOR(IDISP,2)  +
     &                               BORTFAC(ICOMP)*ORCOF(IDISP,ICOMP)
           AVEEOPE(IDISP,2)= AVEEOPE(IDISP,2)+
     &                               BORTFAC(ICOMP)*EOPECF(IDISP,ICOMP)
           AVESHG(IDISP,3) = AVESHG(IDISP,3) +
     &                               BM23FAC(ICOMP)*SHGCF(IDISP,ICOMP)
           AVEOR(IDISP,3)  = AVEOR(IDISP,3)  +
     &                               BM23FAC(ICOMP)*ORCOF(IDISP,ICOMP)
           AVEEOPE(IDISP,3)= AVEEOPE(IDISP,3)+
     &                               BM23FAC(ICOMP)*EOPECF(IDISP,ICOMP)
           AVESHG(IDISP,4) = AVESHG(IDISP,4) +
     &                               BKERFAC(ICOMP)*SHGCF(IDISP,ICOMP)
           AVEOR(IDISP,4)  = AVEOR(IDISP,4)  +
     &                               BKERFAC(ICOMP)*ORCOF(IDISP,ICOMP)
           AVEEOPE(IDISP,4)= AVEEOPE(IDISP,4)+
     &                               BKERFAC(ICOMP)*EOPECF(IDISP,ICOMP)
          END DO
         END DO
        END IF


        IF (LOCDBG) THEN
          WRITE (LUPRI,'(//,A)')
     &    'DISPERSION COEFF. FOR PARALLEL/ORTHOGONAL/23/Kerr AVERAGE:'
          DO IDISP = 1, (NQRDSPE+2)*(NQRDSPE+1)/2
          END DO
          DO MN = 0, NQRDSPE
          DO M = 0, MN
            N = MN - M
            WRITE (LUPRI,'(2I5,4G16.8)') N, M,
     &           (AVEDISP(ITRI(MN,M),I),I=1,4)
          END DO
          END DO

          WRITE (LUPRI,'(//,A)')
     &       'SHG COEFFICIENTS FOR PARALLEL/ORTOGONAL/23/Kerr AVERAGE:'
          DO IDISP = 1, NQRDSPE+1
            WRITE (LUPRI,'(I5,4G16.8)') IDISP-1, (AVESHG(IDISP,I),I=1,4)
          END DO

          WRITE (LUPRI,'(//,A)')
     &      'OR COEFFICIENTS FOR PARALLEL/ORTHOGONAL/23/Kerr AVERAGE:'
          DO IDISP = 1, NQRDSPE+1
            WRITE (LUPRI,'(I5,4G16.8)') IDISP-1, (AVEOR(IDISP,I),I=1,4)
          END DO

          WRITE (LUPRI,'(//,A)')
     &     'EOPE COEFFICIENTS FOR PARALLEL/ORTHOGONAL/23/Kerr AVERAGE:'
          DO IDISP = 1, NQRDSPE+1
            WRITE (LUPRI,'(I5,4G16.8)') IDISP-1,(AVEEOPE(IDISP,I),I=1,4)
          END DO
        END IF

* calculate A, B, C, etc. coefficients:
        BETA0 = AVEDISP(ITRI(0,0),1)
        IF (NQRDSPE.GE.2) THEN
          A0 = AVEDISP(ITRI(2,0),1) / (2.0d0 * BETA0)
          A1 = AVEDISP(ITRI(2,1),1) / (2.0d0 * BETA0)
          A2 = AVEDISP(ITRI(2,2),1) / (2.0d0 * BETA0)
          ERROR = ( DABS(A1-A0) + DABS(A2-A0) ) / DABS(A0)
          IF ( ERROR .GT. 1.0d-10 ) THEN
            WRITE (LUPRI,*) 'ERROR DURING CALCULATION OF A COEFFICIENT:'
            WRITE (LUPRI,*) 'A0:',A0
            WRITE (LUPRI,*) 'A1:',A1
            WRITE (LUPRI,*) 'A2:',A2
            CALL QUIT('ERROR DURING CALCULATION OF A COEFFICIENT.')
          END IF
          ABCDE(KA,1) = A0
          IF (LOCDBG) THEN
           WRITE (LUPRI,'(//,A)') 
     &            'A,B,C,D,E COEFF. FOR PARALLEL/23 AVERAGES:'
           WRITE (LUPRI,*) ' A     coefficient: ',A0
          END IF
          
          A0 = AVEDISP(ITRI(2,0),3) / ( 1.0d0 * BETA0)
          A1 = AVEDISP(ITRI(2,1),3) / (-2.0d0 * BETA0)
          A2 = AVEDISP(ITRI(2,2),3) / (-2.0d0 * BETA0)
          ERROR = ( DABS(A1-A0) + DABS(A2-A0) ) / DABS(A0)
          IF ( ERROR .GT. 1.0d-10 ) THEN
            WRITE (LUPRI,*) 
     &            'ERROR DURING CALCULATION OF A_23 COEFFICIENT:'
            WRITE (LUPRI,*) 'A0:',A0
            WRITE (LUPRI,*) 'A1:',A1
            WRITE (LUPRI,*) 'A2:',A2
            CALL QUIT('ERROR DURING CALCULATION OF A_23 COEFFICIENT.')
          END IF
          ABCDE(KA,2) = A0
          IF (LOCDBG) THEN
           WRITE (LUPRI,*) ' A_23  coefficient: ',A0
          END IF
        END IF
        IF (NQRDSPE.GE.4) THEN
          B0 = AVEDISP(ITRI(4,0),1) / ( 4.0d0 * BETA0)
          B1 = AVEDISP(ITRI(4,1),1) / ( 8.0d0 * BETA0)
          B2 = AVEDISP(ITRI(4,2),1) / (12.0d0 * BETA0)
          B3 = AVEDISP(ITRI(4,3),1) / ( 8.0d0 * BETA0)
          B4 = AVEDISP(ITRI(4,4),1) / ( 4.0d0 * BETA0)
          ERROR = ( DABS(B1-B0)+DABS(B2-B0)+DABS(B3-B0)+DABS(B4-B0) )/
     &                    DABS(B0)
          IF ( ERROR .GT. 1.0d-10 ) THEN
            WRITE (LUPRI,*) 'ERROR DURING CALCULATION OF B COEFFICIENT:'
            WRITE (LUPRI,*) 'B0:',B0
            WRITE (LUPRI,*) 'B1:',B1
            WRITE (LUPRI,*) 'B2:',B2
            WRITE (LUPRI,*) 'B3:',B3
            WRITE (LUPRI,*) 'B4:',B4
            CALL QUIT('ERROR DURING CALCULATION OF B COEFFICIENT.')
          END IF
          ABCDE(KB,1) = B0
          IF (LOCDBG) THEN
            WRITE (LUPRI,*) ' B     coefficient: ',B0
          END IF

          B0  = AVEDISP(ITRI(4,0),3) / ( 2.0d0 * BETA0)
          B1  = AVEDISP(ITRI(4,3),3) / (-8.0d0 * BETA0)
          B2  = AVEDISP(ITRI(4,4),3) / (-4.0d0 * BETA0)
          BP0 = AVEDISP(ITRI(4,1),3) + 1.0d0 * AVEDISP(ITRI(4,0),3)
          BP1 = AVEDISP(ITRI(4,2),3) + 3.0d0 * AVEDISP(ITRI(4,0),3)
          BP0 = BP0 / (-9.0d0 * BETA0)
          BP1 = BP1 / (-9.0d0 * BETA0)
          ERROR = ( DABS(B1-B0)+DABS(B2-B0)+DABS(BP1-BP0))/ DABS(B0)
          IF ( ERROR .GT. 1.0d-10 ) THEN
            WRITE (LUPRI,*) 
     &            'ERROR DURING CALCULATION OF B_23 COEFFICIENT:'
            WRITE (LUPRI,*) 'B0:',B0
            WRITE (LUPRI,*) 'B1:',B1
            WRITE (LUPRI,*) 'B2:',B2
            WRITE (LUPRI,*) 'BP0:',BP0
            WRITE (LUPRI,*) 'BP1:',BP1
            CALL QUIT('ERROR DURING CALCULATION OF B_23 COEFFICIENT.')
          END IF
          ABCDE(KB,2)  = B0
          ABCDE(KBP,2) = BP0
          IF (LOCDBG) THEN
            WRITE (LUPRI,*) " B_23, B_23'  coefficients: ",B0, BP0
          END IF
        END IF
        IF (NQRDSPE.GE.6) THEN
          C0  = AVEDISP(ITRI(6,0),1) / ( 8.0d0 * BETA0)
          C1  = AVEDISP(ITRI(6,1),1) / (24.0d0 * BETA0)
          C2  = AVEDISP(ITRI(6,5),1) / (24.0d0 * BETA0)
          C3  = AVEDISP(ITRI(6,6),1) / ( 8.0d0 * BETA0)
          CP0 = AVEDISP(ITRI(6,2),1) - 6.0d0 * AVEDISP(ITRI(6,0),1) 
          CP1 = AVEDISP(ITRI(6,3),1) - 7.0d0 * AVEDISP(ITRI(6,0),1) 
          CP2 = AVEDISP(ITRI(6,4),1) - 6.0d0 * AVEDISP(ITRI(6,0),1) 
          CP0 = CP0 / ( 9.0d0 * BETA0)
          CP1 = CP1 / (18.0d0 * BETA0)
          CP2 = CP2 / ( 9.0d0 * BETA0)
          ERROR = ( DABS(C1-C0)   + DABS(C2-C0) + DABS(C3-C0) +
     &              DABS(CP1-CP0) + DABS(CP2-CP0) ) / DABS(C0)
          IF ( ERROR .GT. 1.0d-10 ) THEN
            WRITE (LUPRI,*)
     &            'ERROR DURING CALCULATION OF C COEFFICIENTS:'
            WRITE (LUPRI,*) 'C0:',C0
            WRITE (LUPRI,*) 'C1:',C1
            WRITE (LUPRI,*) 'C2:',C2
            WRITE (LUPRI,*) 'C3:',C3
            WRITE (LUPRI,*) 'CP0:',CP0
            WRITE (LUPRI,*) 'CP1:',CP1
            WRITE (LUPRI,*) 'CP2:',CP2
            CALL QUIT('ERROR DURING CALCULATION OF C COEFFICIENTS.')
          END IF
          ABCDE(KC,1)  = C0
          ABCDE(KCP,1) = CP0
          IF (LOCDBG) THEN
            WRITE (LUPRI,*) " C, C' coefficients:",C0,CP0
          END IF

          C0  = AVEDISP(ITRI(6,0),3) / (  4.0d0 * BETA0)
          C1  = AVEDISP(ITRI(6,5),3) / (-24.0d0 * BETA0)
          C2  = AVEDISP(ITRI(6,6),3) / ( -8.0d0 * BETA0)
          CP0 = AVEDISP(ITRI(6,1),3) 
          CP1 = AVEDISP(ITRI(6,2),3) + 3.0d0 * AVEDISP(ITRI(6,0),3) 
          CP2 = AVEDISP(ITRI(6,3),3) + 8.0d0 * AVEDISP(ITRI(6,0),3) 
          CP3 = AVEDISP(ITRI(6,4),3) + 9.0d0 * AVEDISP(ITRI(6,0),3) 
          CP0 = CP0 / (-18.0d0 * BETA0)
          CP1 = CP1 / (-36.0d0 * BETA0)
          CP2 = CP2 / (-36.0d0 * BETA0)
          CP3 = CP3 / (-18.0d0 * BETA0)
          ERROR = ( DABS(C1-C0)   + DABS(C2-C0) + DABS(CP3-CP0) +
     &              DABS(CP1-CP0) + DABS(CP2-CP0) ) / DABS(C0)
          IF ( ERROR .GT. 1.0d-10 ) THEN
            WRITE (LUPRI,*) 
     &            'ERROR DURING CALCULATION OF C_23 COEFFICIENTS:'
            WRITE (LUPRI,*) 'C0:',C0
            WRITE (LUPRI,*) 'C1:',C1
            WRITE (LUPRI,*) 'C2:',C2
            WRITE (LUPRI,*) 'CP0:',CP0
            WRITE (LUPRI,*) 'CP1:',CP1
            WRITE (LUPRI,*) 'CP2:',CP2
            WRITE (LUPRI,*) 'CP3:',CP3
            CALL QUIT('ERROR DURING CALCULATION OF C_23 COEFFICIENTS.')
          END IF
          ABCDE(KC,2)  = C0
          ABCDE(KCP,2) = CP0
          IF (LOCDBG) THEN
            WRITE (LUPRI,*) " C_23, C_23' coefficients:",C0,CP0
          END IF
        END IF
        IF (NQRDSPE.GE.8) THEN
          D0 = AVEDISP(ITRI(8,0),1) / (16.0d0 * BETA0)
          D1 = AVEDISP(ITRI(8,1),1) / (64.0d0 * BETA0)
          D2 = AVEDISP(ITRI(8,7),1) / (64.0d0 * BETA0)
          D3 = AVEDISP(ITRI(8,8),1) / (16.0d0 * BETA0)
          DP0 = AVEDISP(ITRI(8,2),1) - 10.0d0 * AVEDISP(ITRI(8,0),1)
          DP1 = AVEDISP(ITRI(8,3),1) - 16.0d0 * AVEDISP(ITRI(8,0),1)
          DP2 = AVEDISP(ITRI(8,4),1) - 19.0d0 * AVEDISP(ITRI(8,0),1)
          DP3 = AVEDISP(ITRI(8,5),1) - 16.0d0 * AVEDISP(ITRI(8,0),1)
          DP4 = AVEDISP(ITRI(8,6),1) - 10.0d0 * AVEDISP(ITRI(8,0),1)
          DP0 = DP0 / (18.0d0 * BETA0)
          DP1 = DP1 / (54.0d0 * BETA0)
          DP2 = DP2 / (72.0d0 * BETA0)
          DP3 = DP3 / (54.0d0 * BETA0)
          DP4 = DP4 / (18.0d0 * BETA0)
          ERROR = ( DABS(D1-D0)   + DABS(D2-D0)   + DABS(D3-D0) +
     &              DABS(DP1-DP0) + DABS(DP2-DP0) + 
     &              DABS(DP3-DP0) + DABS(DP4-DP0)  ) / DABS(D0)
          IF ( ERROR .GT. 1.0d-10 ) THEN
            WRITE (LUPRI,*) 
     &           'ERROR DURING CALCULATION OF D COEFFICIENTS:'
            WRITE (LUPRI,*) 'D0:',D0
            WRITE (LUPRI,*) 'D1:',D1
            WRITE (LUPRI,*) 'D2:',D2
            WRITE (LUPRI,*) 'D3:',D3
            WRITE (LUPRI,*) 'DP0:',DP0
            WRITE (LUPRI,*) 'DP1:',DP1
            WRITE (LUPRI,*) 'DP2:',DP2
            WRITE (LUPRI,*) 'DP3:',DP3
            WRITE (LUPRI,*) 'DP4:',DP4
            CALL QUIT('ERROR DURING CALCULATION OF D COEFFICIENTS.')
          END IF
          ABCDE(KD,1)  = D0
          ABCDE(KDP,1) = DP0
          IF (LOCDBG) THEN
            WRITE (LUPRI,*) " D, D' coefficients:",D0,DP0
          END IF

          D0   = AVEDISP(ITRI(8,0),3) / (  8.0d0 * BETA0)
          D1   = AVEDISP(ITRI(8,7),3) / (-64.0d0 * BETA0)
          D2   = AVEDISP(ITRI(8,8),3) / (-16.0d0 * BETA0)
          DP0  = AVEDISP(ITRI(8,1),3) -  1.0d0 * AVEDISP(ITRI(8,0),3)
          DP1  = AVEDISP(ITRI(8,3),3) + 11.0d0 * AVEDISP(ITRI(8,0),3)
          DPP0 = AVEDISP(ITRI(8,2),3) -  3.0d0 * AVEDISP(ITRI(8,1),3) +
     &                                   5.0d0 * AVEDISP(ITRI(8,0),3)
          DPP1 = AVEDISP(ITRI(8,4),3) -  5.0d0 * AVEDISP(ITRI(8,1),3) +
     &                                  25.0d0 * AVEDISP(ITRI(8,0),3)
          DPP2 = AVEDISP(ITRI(8,5),3) -  3.0d0 * AVEDISP(ITRI(8,1),3) +
     &                                  26.0d0 * AVEDISP(ITRI(8,0),3)
          DPP3 = AVEDISP(ITRI(8,6),3) -  1.0d0 * AVEDISP(ITRI(8,1),3) +
     &                                  18.0d0 * AVEDISP(ITRI(8,0),3)
          DP0  = DP0  / (  -36.0d0 * BETA0)
          DP1  = DP1  / ( -180.0d0 * BETA0)
          DPP0 = DPP0 / (   9.0d0 * BETA0)
          DPP1 = DPP1 / ( -45.0d0 * BETA0)
          DPP2 = DPP2 / ( -54.0d0 * BETA0)
          DPP3 = DPP3 / ( -18.0d0 * BETA0)
          ERROR = ( DABS(D1-D0)     + DABS(D2-D0)     + 
     &              DABS(DP1-DP0)   + DABS(DPP1-DPP0) + 
     &              DABS(DPP2-DPP0) + DABS(DPP3-DPP0)  ) / DABS(D0)
          IF ( ERROR .GT. 1.0d-10 ) THEN
            WRITE (LUPRI,*)
     &            'ERROR DURING CALCULATION OF D COEFFICIENTS:'
            WRITE (LUPRI,*) 'D0:',D0
            WRITE (LUPRI,*) 'D1:',D1
            WRITE (LUPRI,*) 'D2:',D2
            WRITE (LUPRI,*) 'DP0:',DP0
            WRITE (LUPRI,*) 'DP1:',DP1
            WRITE (LUPRI,*) 'DPP0:',DPP0
            WRITE (LUPRI,*) 'DPP1:',DPP1
            WRITE (LUPRI,*) 'DPP2:',DPP2
            WRITE (LUPRI,*) 'DPP3:',DPP3
            CALL QUIT('ERROR DURING CALCULATION OF D COEFFICIENTS.')
          END IF
          ABCDE(KD,2)   = D0
          ABCDE(KDP,2)  = DP0
          ABCDE(KDP2,2) = DPP0
          IF (LOCDBG) THEN
            WRITE (LUPRI,'(A,3G16.8)')
     &          " D_23, D_23', D_23'' coefficients:",D0,DP0,DPP0
          END IF
        END IF
        IF (NQRDSPE.GE.10) THEN
          E0 = AVEDISP(ITRI(10,0),1)  / ( 32.0d0 * BETA0)
          E1 = AVEDISP(ITRI(10,1),1)  / (160.0d0 * BETA0)
          E2 = AVEDISP(ITRI(10,9),1)  / (160.0d0 * BETA0)
          E3 = AVEDISP(ITRI(10,10),1) / ( 32.0d0 * BETA0)
          EP0 = AVEDISP(ITRI(10,2),1) - 15.0d0 * AVEDISP(ITRI(10,0),1)
          EP1 = AVEDISP(ITRI(10,3),1) - 30.0d0 * AVEDISP(ITRI(10,0),1)
          EP2 = AVEDISP(ITRI(10,4),1) - 45.0d0 * AVEDISP(ITRI(10,0),1)
          EP3 = AVEDISP(ITRI(10,5),1) - 51.0d0 * AVEDISP(ITRI(10,0),1)
          EP4 = AVEDISP(ITRI(10,6),1) - 45.0d0 * AVEDISP(ITRI(10,0),1)
          EP5 = AVEDISP(ITRI(10,7),1) - 30.0d0 * AVEDISP(ITRI(10,0),1)
          EP6 = AVEDISP(ITRI(10,8),1) - 15.0d0 * AVEDISP(ITRI(10,0),1)
          EP0 = EP0 / ( 36.0d0 * BETA0)
          EP1 = EP1 / (144.0d0 * BETA0)
          EP2 = EP2 / (288.0d0 * BETA0)
          EP3 = EP3 / (360.0d0 * BETA0)
          EP4 = EP4 / (288.0d0 * BETA0)
          EP5 = EP5 / (144.0d0 * BETA0)
          EP6 = EP6 / ( 36.0d0 * BETA0)
          ERROR = ( DABS(E1-E0)   + DABS(E2-E0)   + DABS(E3-E0) +
     &              DABS(EP1-EP0) + DABS(EP2-EP0) + 
     &              DABS(EP3-EP0) + DABS(EP4-EP0) +
     &              DABS(EP5-EP0) + DABS(EP6-EP0)  ) / DABS(E0)
          IF ( ERROR .GT. 1.0d-10 ) THEN
            WRITE (LUPRI,*)
     &            'ERROR DURING CALCULATION OF E COEFFICIENTS:'
            WRITE (LUPRI,*) 'E0: ',E0
            WRITE (LUPRI,*) 'E1: ',E1
            WRITE (LUPRI,*) 'E2: ',E2
            WRITE (LUPRI,*) 'E3: ',E3
            WRITE (LUPRI,*) 'EP0:',EP0
            WRITE (LUPRI,*) 'EP1:',EP1
            WRITE (LUPRI,*) 'EP2:',EP2
            WRITE (LUPRI,*) 'EP3:',EP3
            WRITE (LUPRI,*) 'EP4:',EP4
            WRITE (LUPRI,*) 'EP5:',EP5
            WRITE (LUPRI,*) 'EP6:',EP6
            CALL QUIT('ERROR DURING CALCULATION OF E COEFFICIENTS.')
          END IF
          ABCDE(KE,1)  = E0
          ABCDE(KEP,1) = EP0
          IF (LOCDBG) THEN
            WRITE (LUPRI,*) " E, E' coefficients:",E0,EP0
          END IF

          E0 = AVEDISP(ITRI(10,0),3)  / (  16.0d0 * BETA0)
          E1 = AVEDISP(ITRI(10,9),3)  / (-160.0d0 * BETA0)
          E2 = AVEDISP(ITRI(10,10),3) / ( -32.0d0 * BETA0)
          EP0 = AVEDISP(ITRI(10,1),3) - 2.0d0 * AVEDISP(ITRI(10,0),3)
          EP1 = AVEDISP(ITRI(10,8),3) +27.0d0 * AVEDISP(ITRI(10,0),3)+
     &                                  2.0d0 * AVEDISP(ITRI(10,2),3)
          EP0 = EP0 / (-72.0d0 * BETA0)
          EP1 = EP1 / (-648.0d0 * BETA0)
      
          EPP0 = AVEDISP(ITRI(10,2),3) - 4.0d0*AVEDISP(ITRI(10,1),3) +
     &                                   8.0d0*AVEDISP(ITRI(10,0),3)
          EPP1 = AVEDISP(ITRI(10,7),3) - 4.0d0*AVEDISP(ITRI(10,1),3) +
     &                                  56.0d0*AVEDISP(ITRI(10,0),3)
          EPP2 = AVEDISP(ITRI(10,8),3) - 1.0d0*AVEDISP(ITRI(10,1),3) +
     &                                  29.0d0*AVEDISP(ITRI(10,0),3)
          EPP0 = EPP0 / (  18.0d0 * BETA0)
          EPP1 = EPP1 / (-144.0d0 * BETA0)
          EPP2 = EPP2 / ( -36.0d0 * BETA0)

          EPPP0= AVEDISP(ITRI(10,3),3) - 1.0d0 * AVEDISP(ITRI(10,2),3) 
     &     - 5.0d0*AVEDISP(ITRI(10,1),3) + 22.0d0*AVEDISP(ITRI(10,0),3)
          EPPP1= AVEDISP(ITRI(10,4),3) + 4.0d0 * AVEDISP(ITRI(10,2),3)
     &     -29.0d0*AVEDISP(ITRI(10,1),3) + 91.0d0*AVEDISP(ITRI(10,0),3)
          EPPP2= AVEDISP(ITRI(10,5),3) + 11.0d0 * AVEDISP(ITRI(10,2),3)
     &     -57.0d0*AVEDISP(ITRI(10,1),3) +168.0d0*AVEDISP(ITRI(10,0),3)
          EPPP3= AVEDISP(ITRI(10,6),3) + 13.0d0 * AVEDISP(ITRI(10,2),3)
     &     -61.0d0*AVEDISP(ITRI(10,1),3) +182.0d0*AVEDISP(ITRI(10,0),3)
          EPPP0 = EPPP0 / ( -81.0d0 * BETA0)
          EPPP1 = EPPP1 / (-243.0d0 * BETA0)
          EPPP2 = EPPP2 / (-243.0d0 * BETA0)
          EPPP3 = EPPP3 / ( -81.0d0 * BETA0)
          ERROR = ( DABS(E1-E0)       + DABS(E2-E0)   + DABS(EP1-EP0) +
     &              DABS(EPP1-EPP0)   + DABS(EPP2-EPP0) + 
     &              DABS(EPPP1-EPPP0) + DABS(EPPP2-EPPP0) + 
     &              DABS(EPPP3-EPPP0) ) / DABS(E0)
          IF ( ERROR .GT. 1.0d-10 ) THEN
            WRITE (LUPRI,*)
     &            'ERROR DURING CALCULATION OF E_23 COEFFICIENTS:'
            WRITE (LUPRI,*) 'E0: ',E0
            WRITE (LUPRI,*) 'E1: ',E1
            WRITE (LUPRI,*) 'E2: ',E2
            WRITE (LUPRI,*) 'EP0:',EP0
            WRITE (LUPRI,*) 'EP1:',EP1
            WRITE (LUPRI,*) 'EPP0:',EPP0
            WRITE (LUPRI,*) 'EPP1:',EPP1
            WRITE (LUPRI,*) 'EPP2:',EPP2
            WRITE (LUPRI,*) 'EPPP0:',EPPP0
            WRITE (LUPRI,*) 'EPPP1:',EPPP1
            WRITE (LUPRI,*) 'EPPP2:',EPPP2
            WRITE (LUPRI,*) 'EPPP3:',EPPP3
            CALL QUIT('ERROR DURING CALCULATION OF E_23 COEFFICIENTS.')
          END IF
          ABCDE(KE,2)   = E0
          ABCDE(KEP,2)  = EP0
          ABCDE(KEP2,2) = EPP0
          ABCDE(KEP3,2) = EPPP0
          IF (LOCDBG) THEN
            WRITE (LUPRI,'(A,4G16.8)')
     &        " E_23, E_23', E_23'', E_23''' coefficients:",
     &        E0, EP0, EPP0, EPPP0
          END IF
        END IF

      END IF

      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCQRDISP                             *
*---------------------------------------------------------------------*

c /* deck ccqhypprt */
*=====================================================================*
       SUBROUTINE CCQHYPPRT(HYPERPOL)
*---------------------------------------------------------------------*
*
*    Purpose: print frequency-dependent first hyperpolarizabilities 
*
*
*     Written by Christof Haettig in October 1996.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccqrinf.h>
#include <ccroper.h>


      CHARACTER*15  BLANKS
      CHARACTER*10 RLXA, RLXB, RLXC
      CHARACTER*80 STRING
      CHARACTER*8 BETALAB
      INTEGER ISYMA, ISYMB, ISYMC, ISAMA, ISAMB, ISAMC, ISAPROP
      INTEGER IFREQ, IOPER, ICMP, I, J, K, L


#if defined (SYS_CRAY)
      REAL HYPERPOL(NQRFREQ,NQROPER,2)
      REAL HALF, ONE, SIGN, BETA, PROP, ERROR, THIRD, ZERO
      REAL FACTOR
      REAL BETANEW(NQRFREQ,3,3,3)
      REAL MEANBETA, MUBETANORM, MU_DOT_BETA
      REAL BETAVEC(3), MUBETA(3)
#else
      DOUBLE PRECISION HYPERPOL(NQRFREQ,NQROPER,2)
      DOUBLE PRECISION HALF, ONE, SIGN, BETA, PROP, ERROR, THIRD, ZERO
      DOUBLE PRECISION FACTOR
      DOUBLE PRECISION BETANEW(NQRFREQ,3,3,3)
      DOUBLE PRECISION MEANBETA, MUBETANORM, MU_DOT_BETA
      DOUBLE PRECISION BETAVEC(3), MUBETA(3)
#endif
      PARAMETER (HALF = 0.5d0, ONE = 1.0d0, THIRD = 1.0d0/3.0d0)
      PARAMETER (ZERO = 0.0d0, FACTOR = 1.0d0/6.0d0)

      CHARACTER*10  MODEL,MOPRPC
      INTEGER ISYMAB,ISYMABC
*---------------------------------------------------------------------*
* print header for first hyperpolarizabilities: 
*---------------------------------------------------------------------*
      BLANKS = '               '
      STRING = ' RESULTS FOR THE FIRST HYPERPOLARIZABILITIES '

      IF (CCS) THEN
         CALL AROUND120( BLANKS//'FINAL CCS'//STRING(1:45)//BLANKS ) 
      ELSE IF (CC2) THEN
         CALL AROUND120( BLANKS//'FINAL CC2'//STRING(1:45)//BLANKS )
      ELSE IF (CCSD) THEN
         CALL AROUND120( BLANKS//'FINAL CCSD'//STRING(1:45)//BLANKS )
      ELSE IF (CC3) THEN
         CALL AROUND120( BLANKS//'FINAL CC3'//STRING(1:45)//BLANKS )
      ELSE
         CALL QUIT('CCQHYPPRT called for an unknown '//
     &             'Coupled Cluster model.')
      END IF

      IF (CCSD) THEN
         MODEL = 'CCSD'
         MOPRPC = 'CCSD      '
      ELSE IF (CIS) THEN
         MODEL  = 'CIS'
         MOPRPC = 'CIS       '
      ELSE IF (CCS) THEN
         MODEL  = 'CCS'
         MOPRPC = 'CCS       '
      ELSE IF (CC2) THEN
         MODEL  = 'CC2'
         MOPRPC = 'CC2       '
      ELSE IF (CC3) THEN
         MODEL  = 'CC3'
         MOPRPC = 'CC3       '
      END IF

      IF (IPRQHYP.GT.5) THEN
       WRITE(LUPRI,'(/1X,3(1X,A," operator",17X),4X,A,6X,A,/,130("-"))')
     &      'A','B','C','property ','(asy./sym. Resp.)'
      ELSE
       WRITE(LUPRI,'(/1X,3(1X,A," operator",17X),4X,A,/,102("-"))')
     &      'A','B','C','property '
      END IF


      DO IOPER = 1, NQROPER
         ISYMA = ISYOPR(IAQROP(IOPER))
         ISYMB = ISYOPR(IBQROP(IOPER))
         ISYMC = ISYOPR(ICQROP(IOPER))

         ISAMA = ISYMAT(IAQROP(IOPER))
         ISAMB = ISYMAT(IBQROP(IOPER))
         ISAMC = ISYMAT(ICQROP(IOPER))

         IF (       LAQLRX(IOPER) ) RLXA = ' (relax.) '
         IF ( .NOT. LAQLRX(IOPER) ) RLXA = ' (unrel.) '
         IF (       LBQLRX(IOPER) ) RLXB = ' (relax.) '
         IF ( .NOT. LBQLRX(IOPER) ) RLXB = ' (unrel.) '
         IF (       LCQLRX(IOPER) ) RLXC = ' (relax.) '
         IF ( .NOT. LCQLRX(IOPER) ) RLXC = ' (unrel.) '

         ISAPROP = ISAMA * ISAMB * ISAMC
         SIGN    = DBLE(ISAPROP)
         IF (ISAPROP.EQ.0) SIGN = +ONE

         IFREQ = 1
cmbh
         PROP=0.0d0
cmbh end
         IF (MULD2H(ISYMA,ISYMB).EQ.ISYMC) THEN

           PROP  = -HALF*(HYPERPOL(IFREQ,IOPER,1) +
     &                     SIGN * HYPERPOL(IFREQ,IOPER,2))
           ERROR = -HALF*(HYPERPOL(IFREQ,IOPER,1) -
     &                     SIGN * HYPERPOL(IFREQ,IOPER,2))

           IF (IPRQHYP.GT.5 .OR. ISAPROP.EQ.0) THEN
            IF (ISAPROP.NE.0) THEN
              WRITE(LUPRI,
     &               '(/1X,3(A8,A10,F7.4,3X),G18.10," (",G18.10,")")')
     &          LBLOPR(IAQROP(IOPER)),RLXA,AQRFR(IFREQ),
     &          LBLOPR(IBQROP(IOPER)),RLXB,BQRFR(IFREQ),
     &          LBLOPR(ICQROP(IOPER)),RLXC,CQRFR(IFREQ), PROP, ERROR
            ELSE
              WRITE(LUPRI,'(/1X,2A)')
     &          'Cannot determine if real or imaginary property...   ',
     &          'print: symmetric (antisymmetric) parts in +/- w'
              WRITE(LUPRI,
     &               '(1X,3(A8,A10,F7.4,3X),G18.10," (",G18.10,")")')
     &          LBLOPR(IAQROP(IOPER)),RLXA,AQRFR(IFREQ),
     &          LBLOPR(IBQROP(IOPER)),RLXB,BQRFR(IFREQ),
     &          LBLOPR(ICQROP(IOPER)),RLXC,CQRFR(IFREQ), PROP, ERROR
            END IF
           ELSE
              WRITE(LUPRI,'(/1X,3(A8,A10,F7.4,3X),G16.8)')
CCN              WRITE(LUPRI,'(/1X,3(A8,A10,F7.4,3X),F24.20)')
     &        LBLOPR(IAQROP(IOPER)),RLXA,AQRFR(IFREQ),
     &        LBLOPR(IBQROP(IOPER)),RLXB,BQRFR(IFREQ),
     &        LBLOPR(ICQROP(IOPER)),RLXC,CQRFR(IFREQ), PROP
           ENDIF
         ELSE
           IF (IPRQHYP.GT.5 .OR. ISAPROP.EQ.0) THEN
            WRITE(LUPRI,
     &              '(/1X,3(A8,A10,A7,3X),7X,A,8X," (",9X,A,10X,")")')
     &        LBLOPR(IAQROP(IOPER)),RLXA,'   -.-  ',
     &        LBLOPR(IBQROP(IOPER)),RLXB,'   -.-  ',
     &        LBLOPR(ICQROP(IOPER)),RLXC,'   -.-  ', 
     &        '---',
     &        '---'
           ELSE IF (IPRQHYP.GT.0) THEN
            WRITE(LUPRI,'(/1X,3(A8,A10,A7,3X),6X,A,7X)')
     &        LBLOPR(IAQROP(IOPER)),RLXA,'   -.-  ',
     &        LBLOPR(IBQROP(IOPER)),RLXB,'   -.-  ',
     &        LBLOPR(ICQROP(IOPER)),RLXC,'   -.-  ', 
     &        '---'
           END IF
         END IF
         ISYMAB  = MULD2H(ISYMA,ISYMB)
         ISYMABC = MULD2H(ISYMAB,ISYMC)
         CALL WRIPRO(PROP,MOPRPC,3,
     &               LBLOPR(IAQROP(IOPER)),LBLOPR(IBQROP(IOPER)),
     &               LBLOPR(ICQROP(IOPER)),LBLOPR(ICQROP(IOPER)),
     &               BQRFR(IFREQ),CQRFR(IFREQ),CQRFR(IFREQ),ISYMABC,
     &               0,0,0)
         DO IFREQ = 2, NQRFREQ
           IF (MULD2H(ISYMA,ISYMB).EQ.ISYMC) THEN

             PROP  = -HALF*(HYPERPOL(IFREQ,IOPER,1) +
     &                      SIGN * HYPERPOL(IFREQ,IOPER,2))
             ERROR = -HALF*(HYPERPOL(IFREQ,IOPER,1) -
     &                      SIGN * HYPERPOL(IFREQ,IOPER,2))

             IF (IPRQHYP.GT.5 .OR. ISAPROP.EQ.0) THEN
               WRITE(LUPRI,'(1X,3(18X,F7.4,3X),G18.10," (",G18.10,")")')
     &           AQRFR(IFREQ), BQRFR(IFREQ), CQRFR(IFREQ), PROP, ERROR
             ELSE
               WRITE(LUPRI,'(1X,3(18X,F7.4,3X),G16.8)')
     &           AQRFR(IFREQ), BQRFR(IFREQ), CQRFR(IFREQ), PROP
             END IF
           END IF
          CALL WRIPRO(PROP,MOPRPC,3,
     &                LBLOPR(IAQROP(IOPER)),LBLOPR(IBQROP(IOPER)),
     &                LBLOPR(ICQROP(IOPER)),LBLOPR(ICQROP(IOPER)),
     &                BQRFR(IFREQ),CQRFR(IFREQ),CQRFR(IFREQ),ISYMABC,
     &                0,0,0)
         END DO

      END DO

      IF (BETA_AVERAGE) THEN
       WRITE(LUPRI,'(/////A,/,55("-"))')'  average         frequencies'

* calculate & print beta_{||}:
       IFREQ = 1
       BETA = -0.3d0*(HYPERPOL(IFREQ,1,1)+HYPERPOL(IFREQ,1,2))
       IF (XY_DEGENERAT) THEN
         DO ICMP = 2, 4
           BETA = BETA-
     &            0.2d0*(HYPERPOL(IFREQ,ICMP,1)+HYPERPOL(IFREQ,ICMP,2))
         END DO
       ELSE
         DO ICMP = 2, 7
           BETA = BETA-
     &            0.1d0*(HYPERPOL(IFREQ,ICMP,1)+HYPERPOL(IFREQ,ICMP,2))
         END DO
       END IF

       WRITE(LUPRI,'(/1X,A8,3X,3(F7.4,2X),G16.8)')
     &   ' beta_||', AQRFR(IFREQ), BQRFR(IFREQ), CQRFR(IFREQ), BETA
C
       BETALAB='beta_|| '
       CALL WRIPRO(BETA,MOPRPC,3,
     *             BETALAB,BETALAB,BETALAB,BETALAB,
     *             AQRFR(IFREQ),BQRFR(IFREQ),CQRFR(IFREQ),ISYMABC,
     *             0,0,0)
C
       DO IFREQ = 2, NQRFREQ
        BETA = -0.3d0*(HYPERPOL(IFREQ,1,1)+HYPERPOL(IFREQ,1,2))
        IF (XY_DEGENERAT) THEN
          DO ICMP = 2, 4
            BETA = BETA-
     &           0.2d0*(HYPERPOL(IFREQ,ICMP,1)+HYPERPOL(IFREQ,ICMP,2))
          END DO
        ELSE
          DO ICMP = 2, 7
            BETA = BETA-
     &           0.1d0*(HYPERPOL(IFREQ,ICMP,1)+HYPERPOL(IFREQ,ICMP,2))
          END DO
        END IF

        WRITE(LUPRI,'(1X,A8,3X,3(F7.4,2X),G16.8)')
     &    '        ', AQRFR(IFREQ), BQRFR(IFREQ), CQRFR(IFREQ), BETA
C
       BETALAB='beta_|| '
       CALL WRIPRO(BETA,MOPRPC,3,
     *             BETALAB,BETALAB,BETALAB,BETALAB,
     *             AQRFR(IFREQ),BQRFR(IFREQ),CQRFR(IFREQ),ISYMABC,
     *             0,0,0)
C
       END DO
   
* calculate & print beta^K:
       IFREQ = 1
       BETA = -0.3d0*(HYPERPOL(IFREQ,1,1)+HYPERPOL(IFREQ,1,2))
       BETA = BETA + 0.15d0*(HYPERPOL(IFREQ,2,1)+HYPERPOL(IFREQ,2,2))
       BETA = BETA - 0.60d0*(HYPERPOL(IFREQ,3,1)+HYPERPOL(IFREQ,3,2))
       BETA = BETA + 0.15d0*(HYPERPOL(IFREQ,4,1)+HYPERPOL(IFREQ,4,2))
       IF (XY_DEGENERAT) THEN
         BETA = BETA + 0.15d0*(HYPERPOL(IFREQ,2,1)+HYPERPOL(IFREQ,2,2))
         BETA = BETA - 0.60d0*(HYPERPOL(IFREQ,3,1)+HYPERPOL(IFREQ,3,2))
         BETA = BETA + 0.15d0*(HYPERPOL(IFREQ,4,1)+HYPERPOL(IFREQ,4,2))
       ELSE
         BETA = BETA + 0.15d0*(HYPERPOL(IFREQ,5,1)+HYPERPOL(IFREQ,5,2))
         BETA = BETA - 0.60d0*(HYPERPOL(IFREQ,6,1)+HYPERPOL(IFREQ,6,2))
         BETA = BETA + 0.15d0*(HYPERPOL(IFREQ,7,1)+HYPERPOL(IFREQ,7,2))
       END IF

       WRITE(LUPRI,'(/1X,A8,3X,3(F7.4,2X),G16.8)')
     &   ' beta^K ', AQRFR(IFREQ), BQRFR(IFREQ), CQRFR(IFREQ), BETA
C
       BETALAB='beta^K  '
       CALL WRIPRO(BETA,MOPRPC,3,
     *             BETALAB,BETALAB,BETALAB,BETALAB,
     *             AQRFR(IFREQ),BQRFR(IFREQ),CQRFR(IFREQ),ISYMABC,
     *             0,0,0)
C
       DO IFREQ = 2, NQRFREQ
        BETA = -0.3d0*(HYPERPOL(IFREQ,1,1)+HYPERPOL(IFREQ,1,2))
        BETA = BETA + 0.15d0*(HYPERPOL(IFREQ,2,1)+HYPERPOL(IFREQ,2,2))
        BETA = BETA - 0.60d0*(HYPERPOL(IFREQ,3,1)+HYPERPOL(IFREQ,3,2))
        BETA = BETA + 0.15d0*(HYPERPOL(IFREQ,4,1)+HYPERPOL(IFREQ,4,2))
        IF (XY_DEGENERAT) THEN
         BETA = BETA + 0.15d0*(HYPERPOL(IFREQ,2,1)+HYPERPOL(IFREQ,2,2))
         BETA = BETA - 0.60d0*(HYPERPOL(IFREQ,3,1)+HYPERPOL(IFREQ,3,2))
         BETA = BETA + 0.15d0*(HYPERPOL(IFREQ,4,1)+HYPERPOL(IFREQ,4,2))
        ELSE
         BETA = BETA + 0.15d0*(HYPERPOL(IFREQ,5,1)+HYPERPOL(IFREQ,5,2))
         BETA = BETA - 0.60d0*(HYPERPOL(IFREQ,6,1)+HYPERPOL(IFREQ,6,2))
         BETA = BETA + 0.15d0*(HYPERPOL(IFREQ,7,1)+HYPERPOL(IFREQ,7,2))
        END IF

        WRITE(LUPRI,'(1X,A8,3X,3(F7.4,2X),G16.8)')
     &    '        ', AQRFR(IFREQ), BQRFR(IFREQ), CQRFR(IFREQ), BETA
C
        BETALAB='beta^K  '
        CALL WRIPRO(BETA,MOPRPC,3,
     *              BETALAB,BETALAB,BETALAB,BETALAB,
     *              AQRFR(IFREQ),BQRFR(IFREQ),CQRFR(IFREQ),ISYMABC,
     *              0,0,0)
C
       END DO
   

* calculate & print beta_{_|_}:
       IFREQ = 1
       BETA = -0.1d0*(HYPERPOL(IFREQ,1,1)+HYPERPOL(IFREQ,1,2))
       BETA = BETA - 0.2d0*(HYPERPOL(IFREQ,2,1)+HYPERPOL(IFREQ,2,2))
       BETA = BETA + 0.3d0*(HYPERPOL(IFREQ,3,1)+HYPERPOL(IFREQ,3,2))
       BETA = BETA - 0.2d0*(HYPERPOL(IFREQ,4,1)+HYPERPOL(IFREQ,4,2))
       IF (XY_DEGENERAT) THEN
         BETA = BETA - 0.2d0*(HYPERPOL(IFREQ,2,1)+HYPERPOL(IFREQ,2,2))
         BETA = BETA + 0.3d0*(HYPERPOL(IFREQ,3,1)+HYPERPOL(IFREQ,3,2))
         BETA = BETA - 0.2d0*(HYPERPOL(IFREQ,4,1)+HYPERPOL(IFREQ,4,2))
       ELSE
         BETA = BETA - 0.2d0*(HYPERPOL(IFREQ,5,1)+HYPERPOL(IFREQ,5,2))
         BETA = BETA + 0.3d0*(HYPERPOL(IFREQ,6,1)+HYPERPOL(IFREQ,6,2))
         BETA = BETA - 0.2d0*(HYPERPOL(IFREQ,7,1)+HYPERPOL(IFREQ,7,2))
       END IF

       WRITE(LUPRI,'(/1X,A8,3X,3(F7.4,2X),G16.8)')
     &   ' beta_|_', AQRFR(IFREQ), BQRFR(IFREQ), CQRFR(IFREQ), BETA
C
       BETALAB='beta_|_ '
       CALL WRIPRO(BETA,MOPRPC,3,
     *             BETALAB,BETALAB,BETALAB,BETALAB,
     *             AQRFR(IFREQ),BQRFR(IFREQ),CQRFR(IFREQ),ISYMABC,
     *             0,0,0)
C
       DO IFREQ = 2, NQRFREQ
        BETA = -0.1d0*(HYPERPOL(IFREQ,1,1)+HYPERPOL(IFREQ,1,2))
        BETA = BETA - 0.2d0*(HYPERPOL(IFREQ,2,1)+HYPERPOL(IFREQ,2,2))
        BETA = BETA + 0.3d0*(HYPERPOL(IFREQ,3,1)+HYPERPOL(IFREQ,3,2))
        BETA = BETA - 0.2d0*(HYPERPOL(IFREQ,4,1)+HYPERPOL(IFREQ,4,2))
        IF (XY_DEGENERAT) THEN
          BETA = 2.0d0 * BETA
        ELSE
          BETA = BETA - 0.2d0*(HYPERPOL(IFREQ,5,1)+HYPERPOL(IFREQ,5,2))
          BETA = BETA + 0.3d0*(HYPERPOL(IFREQ,6,1)+HYPERPOL(IFREQ,6,2))
          BETA = BETA - 0.2d0*(HYPERPOL(IFREQ,7,1)+HYPERPOL(IFREQ,7,2))
        END IF

        WRITE(LUPRI,'(1X,A8,3X,3(F7.4,2X),G16.8)')
     &    '        ', AQRFR(IFREQ), BQRFR(IFREQ), CQRFR(IFREQ), BETA
C
        BETALAB='beta_|_ '
        CALL WRIPRO(BETA,MOPRPC,3,
     *              BETALAB,BETALAB,BETALAB,BETALAB,
     *              AQRFR(IFREQ),BQRFR(IFREQ),CQRFR(IFREQ),ISYMABC,
     *              0,0,0)
C
       END DO
   
* calculate & print beta_ms:
       IFREQ = 1
       BETA = 0.0d0
       BETA = BETA - 0.25d0*(HYPERPOL(IFREQ,2,1)+HYPERPOL(IFREQ,2,2))
       BETA = BETA + 0.50d0*(HYPERPOL(IFREQ,3,1)+HYPERPOL(IFREQ,3,2))
       BETA = BETA - 0.25d0*(HYPERPOL(IFREQ,4,1)+HYPERPOL(IFREQ,4,2))
       IF (XY_DEGENERAT) THEN
         BETA = BETA - 0.25d0*(HYPERPOL(IFREQ,2,1)+HYPERPOL(IFREQ,2,2))
         BETA = BETA + 0.50d0*(HYPERPOL(IFREQ,3,1)+HYPERPOL(IFREQ,3,2))
         BETA = BETA - 0.25d0*(HYPERPOL(IFREQ,4,1)+HYPERPOL(IFREQ,4,2))
       ELSE
         BETA = BETA - 0.25d0*(HYPERPOL(IFREQ,5,1)+HYPERPOL(IFREQ,5,2))
         BETA = BETA + 0.50d0*(HYPERPOL(IFREQ,6,1)+HYPERPOL(IFREQ,6,2))
         BETA = BETA - 0.25d0*(HYPERPOL(IFREQ,7,1)+HYPERPOL(IFREQ,7,2))
       END IF

       WRITE(LUPRI,'(/1X,A11,3(F7.4,2X),G16.8)')
     &   ' beta_ms   ', AQRFR(IFREQ), BQRFR(IFREQ), CQRFR(IFREQ), BETA
C
        BETALAB='beta_ms '
        CALL WRIPRO(BETA,MOPRPC,3,
     *              BETALAB,BETALAB,BETALAB,BETALAB,
     *              AQRFR(IFREQ),BQRFR(IFREQ),CQRFR(IFREQ),ISYMABC,
     *              0,0,0)
C
    

       DO IFREQ = 2, NQRFREQ
        BETA = 0.0d0
        BETA = BETA - 0.25d0*(HYPERPOL(IFREQ,2,1)+HYPERPOL(IFREQ,2,2))
        BETA = BETA + 0.50d0*(HYPERPOL(IFREQ,3,1)+HYPERPOL(IFREQ,3,2))
        BETA = BETA - 0.25d0*(HYPERPOL(IFREQ,4,1)+HYPERPOL(IFREQ,4,2))
        IF (XY_DEGENERAT) THEN
         BETA = 2.0d0 * BETA
        ELSE
         BETA = BETA - 0.25d0*(HYPERPOL(IFREQ,5,1)+HYPERPOL(IFREQ,5,2))
         BETA = BETA + 0.50d0*(HYPERPOL(IFREQ,6,1)+HYPERPOL(IFREQ,6,2))
         BETA = BETA - 0.25d0*(HYPERPOL(IFREQ,7,1)+HYPERPOL(IFREQ,7,2))
        END IF

        WRITE(LUPRI,'(1X,A8,3X,3(F7.4,2X),G16.8)')
     &    '        ', AQRFR(IFREQ), BQRFR(IFREQ), CQRFR(IFREQ), BETA
C
         BETALAB='beta_ms '
         CALL WRIPRO(BETA,MOPRPC,3,
     *               BETALAB,BETALAB,BETALAB,BETALAB,
     *               AQRFR(IFREQ),BQRFR(IFREQ),CQRFR(IFREQ),ISYMABC,
     *               0,0,0)
C
       END DO
      END IF
C-----------------------------------------------------------------------
C     June 04: AO+JK
C     .AVANEW in input block *CCQR -> LAVANEW =.TRUE.
C     Calculates:
C     i)   Beta_i; i=x,y,z
C     ii)  Mu_i*Beta_i -> |Mu*Beta|
C     iii) <Beta>
C-----------------------------------------------------------------------
C
      IF (LAVANEW) THEN
C
        BLANKS = '               '
        STRING = '   AVERAGES FOR SECOND HARMONIC GENERATION   '
C
        IF (CCS) THEN
           CALL AROUND120( BLANKS//'FINAL CCS'//STRING(1:45)//BLANKS )
        ELSE IF (CC2) THEN
           CALL AROUND120( BLANKS//'FINAL CC2'//STRING(1:45)//BLANKS )
        ELSE IF (CCSD) THEN
           CALL AROUND120( BLANKS//'FINAL CCSD'//STRING(1:45)//BLANKS )
        ELSE
           CALL QUIT('CCQHYPPRT called for an unknown '//
     &             'Coupled Cluster model.')
        END IF
C
C----------------------------------
C       Reset BETANEW(IFREQ,I,J,K):
C----------------------------------
C
        DO IFREQ=1,NQRFREQ
          DO I=1,3
            DO J=1,3
              DO K=1,3
                BETANEW(IFREQ,I,J,K) = ZERO
              END DO
            END DO
          END DO
        END DO
C
C---------------------------------------
C       Rearrange the components of beta
C---------------------------------------
C
        WRITE(LUPRI,'(/////20X,A)')
     *  '+----------+----------+----------+----------+----------+'
        WRITE(LUPRI,'(20X,A)')
     *  '| Property:| Freq_A:  | Freq_B:  | Freq_C:  | Value:   |'
        WRITE(LUPRI,'(20X,A)')
     *  '|----------+----------+----------+----------+----------|'
        WRITE(LUPRI,'(20X,A,A,A,F8.4,A,F8.4,A,F8.4,A,F8.4,A)')
     *  '| ','Dipole_x',' | ',ZERO,' | ',ZERO,' | ',ZERO,' | ',
     *  DIPSAVE(1),' |'
        WRITE(LUPRI,'(20X,A,A,A,F8.4,A,F8.4,A,F8.4,A,F8.4,A)')
     *  '| ','Dipole_y',' | ',ZERO,' | ',ZERO,' | ',ZERO,' | ',
     *  DIPSAVE(2),' |'
        WRITE(LUPRI,'(20X,A,A,A,F8.4,A,F8.4,A,F8.4,A,F8.4,A)')
     *  '| ','Dipole_z',' | ',ZERO,' | ',ZERO,' | ',ZERO,' | ',
     *  DIPSAVE(3),' |'
        WRITE(LUPRI,'(20X,A)')
     *  '|----------+----------+----------+----------+----------|'
C
        L = 0
        DO I = 1, 3
          DO J = 1, 3
            DO K = 1, 3
              L = L + 1
              ISAMA = ISYMAT(IAQROP(L))
              ISAMB = ISYMAT(IBQROP(L))
              ISAMC = ISYMAT(ICQROP(L))
              ISYMA = ISYOPR(IAQROP(L))
              ISYMB = ISYOPR(IBQROP(L))
              ISYMC = ISYOPR(ICQROP(L))
              ISYMAB  = MULD2H(ISYMA,ISYMB)
              ISYMABC = MULD2H(ISYMAB,ISYMC)
              ISAPROP = ISAMA * ISAMB * ISAMC
              SIGN = DBLE(ISAPROP)
              IF (ISAPROP .EQ. 0) SIGN = +ONE
              DO IFREQ = 1, NQRFREQ
                BETANEW(IFREQ,I,J,K) = -HALF*(HYPERPOL(IFREQ,L,1)
     *                               +  SIGN * HYPERPOL(IFREQ,L,2))
              END DO
            END DO
          END DO
        END DO
C
C
C----------------------------------------------------------------
C         Calculate the components of the beta_i; i=x,y,z vector:
C----------------------------------------------------------------
C
        DO IFREQ = 1, NQRFREQ
          MUBETANORM = ZERO
          MU_DOT_BETA= ZERO
          MEANBETA   = ZERO
          DO I = 1, 3
            BETAVEC(I) = ZERO
            MUBETA(I)  = ZERO
            DO J = 1, 3
              BETAVEC(I) = BETAVEC(I) + THIRD * (BETANEW(IFREQ,I,J,J) + 
     *                     BETANEW(IFREQ,J,J,I) + BETANEW(IFREQ,J,I,J))
            END DO
              MUBETA(I)  = DIPSAVE(I) * BETAVEC(I) 
C
            IF (I .EQ. 1) BETALAB='beta_x  '
            IF (I .EQ. 2) BETALAB='beta_y  '
            IF (I .EQ. 3) BETALAB='beta_z  '
C
            WRITE(LUPRI,'(20X,A,A,A,F8.4,A,F8.4,A,F8.4,A,F8.4,A)')
     *                  '| ',BETALAB,' | ',AQRFR(IFREQ),' | ',
     *                  BQRFR(IFREQ),' | ',CQRFR(IFREQ),' | ',
     *                  BETAVEC(I),' |'
C
            CALL WRIPRO(BETAVEC(I),MOPRPC,3,
     *                  BETALAB,BETALAB,BETALAB,BETALAB,
     *                  AQRFR(IFREQ),BQRFR(IFREQ),CQRFR(IFREQ),ISYMABC,
     *                  0,0,0)
C 
          END DO
C
          WRITE(LUPRI,'(20X,A)')
     *    '+----------+----------+----------+----------+----------+'
C
          DO I = 1, 3
            MUBETANORM = MUBETANORM + MUBETA(I)*MUBETA(I)
            MU_DOT_BETA= MU_DOT_BETA+ MUBETA(I)
            IF (I .EQ. 1) BETALAB='mu*bet_x'
            IF (I .EQ. 2) BETALAB='mu*bet_y'
            IF (I .EQ. 3) BETALAB='mu*bet_z'
            WRITE(LUPRI,'(20X,A,A,A,F8.4,A,F8.4,A,F8.4,A,F8.4,A)')
     *                  '| ',BETALAB,' | ',AQRFR(IFREQ),' | ',
     *                  BQRFR(IFREQ),' | ',CQRFR(IFREQ),' | ',
     *                  MUBETA(I),' |' 
C
            CALL WRIPRO(MUBETA(I),MOPRPC,3,
     *                  BETALAB,BETALAB,BETALAB,BETALAB,
     *                  AQRFR(IFREQ),BQRFR(IFREQ),CQRFR(IFREQ),ISYMABC,
     *                  0,0,0)
          END DO
C
          WRITE(LUPRI,'(20X,A)')
     *    '+----------+----------+----------+----------+----------+'
C
          BETALAB = ' mu*bet '
          WRITE(LUPRI,'(20X,A,A,A,F8.4,A,F8.4,A,F8.4,A,F8.4,A)')
     *                '| ',BETALAB,' | ',AQRFR(IFREQ),' | ',
     *                BQRFR(IFREQ),' | ',CQRFR(IFREQ),' | ',
     *                MU_DOT_BETA,' |'
C
            CALL WRIPRO(MU_DOT_BETA,MOPRPC,3,
     *                  BETALAB,BETALAB,BETALAB,BETALAB,
     *                  AQRFR(IFREQ),BQRFR(IFREQ),CQRFR(IFREQ),ISYMABC,
     *                  0,0,0)
C
          WRITE(LUPRI,'(20X,A)')
     *    '+----------+----------+----------+----------+----------+'
          MUBETANORM = SQRT(MUBETANORM)
          BETALAB = '|mu*bet|'
          WRITE(LUPRI,'(20X,A,A,A,F8.4,A,F8.4,A,F8.4,A,F8.4,A)')
     *                '| ',BETALAB,' | ',AQRFR(IFREQ),' | ',
     *                BQRFR(IFREQ),' | ',CQRFR(IFREQ),' | ',
     *                MUBETANORM,' |'
C
            CALL WRIPRO(MUBETANORM,MOPRPC,3,
     *                  BETALAB,BETALAB,BETALAB,BETALAB,
     *                  AQRFR(IFREQ),BQRFR(IFREQ),CQRFR(IFREQ),ISYMABC,
     *                  0,0,0)
C
          WRITE(LUPRI,'(20X,A)')
     *    '+----------+----------+----------+----------+----------+'
          MEANBETA   = FACTOR * ( BETANEW(IFREQ,1,2,3) -
     *                 BETANEW(IFREQ,1,3,2) + BETANEW(IFREQ,2,3,1) -
     *                 BETANEW(IFREQ,2,1,3) + BETANEW(IFREQ,3,1,2) -
     *                 BETANEW(IFREQ,3,2,1) )
          BETALAB = '<beta>  '
          WRITE(LUPRI,'(20X,A,A,A,F8.4,A,F8.4,A,F8.4,A,F8.4,A)')
     *                '| ',BETALAB,' | ',AQRFR(IFREQ),' | ',
     *                BQRFR(IFREQ),' | ',CQRFR(IFREQ),' | ',
     *                MEANBETA,' |'
C
            CALL WRIPRO(MEANBETA,MOPRPC,3,
     *                  BETALAB,BETALAB,BETALAB,BETALAB,
     *                  AQRFR(IFREQ),BQRFR(IFREQ),CQRFR(IFREQ),ISYMABC,
     *                  0,0,0)
C
          WRITE(LUPRI,'(20X,A)')
     *    '+----------+----------+----------+----------+----------+'
        END DO
      END IF

      WRITE(LUPRI,'(//,102("-"),///)')

      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCQHYPPRT                            *
*---------------------------------------------------------------------*
c /* deck ccqexpprt */
*=====================================================================*
       SUBROUTINE CCQEXPPRT(EXPCOF)
*---------------------------------------------------------------------*
*
*   Purpose: print expansion coefficients for
*            first hyperpolarizabilities 
*
*
*   Written by Christof Haettig in October 1997.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccqrinf.h>
#include <ccroper.h>


      LOGICAL TSTISA
      CHARACTER*5  BLANKS
      CHARACTER*80 STRING
      INTEGER ISYMA, ISYMB, ISYMC, ISAMA, ISAMB, ISAMC, ISAPROP
      INTEGER IDISP, IOPER, ICTOT, ISACAU, IDISP2


#if defined (SYS_CRAY)
      REAL EXPCOF(NQRDISP,NQROPER)
#else
      DOUBLE PRECISION EXPCOF(NQRDISP,NQROPER)
#endif

*---------------------------------------------------------------------*
* print header for expansion coefficients:
*---------------------------------------------------------------------*
      BLANKS = '     '
      STRING = ' RESULTS FOR EXPANSION COEFFICIENTS'

      IF (CCS) THEN
         CALL AROUND( BLANKS//'FINAL CCS'//STRING(1:45)//BLANKS ) 
      ELSE IF (CC2) THEN
         CALL AROUND( BLANKS//'FINAL CC2'//STRING(1:45)//BLANKS )
      ELSE IF (CCSD) THEN
         CALL AROUND( BLANKS//'FINAL CCSD'//STRING(1:45)//BLANKS )
      ELSE
         CALL QUIT('CCQEXPPRT called for an unknown '//
     &        'Coupled Cluster model.')
      END IF

      WRITE(LUPRI,'(/1X,3(1X,A," operator",3X),4X,A,/,68("-"))')
     &      'A','B','C','d_ABC'


      DO IOPER = 1, NQROPER
         ISYMA = ISYOPR(IAQROP(IOPER))
         ISYMB = ISYOPR(IBQROP(IOPER))
         ISYMC = ISYOPR(ICQROP(IOPER))

         ISAMA = ISYMAT(IAQROP(IOPER))
         ISAMB = ISYMAT(IBQROP(IOPER))
         ISAMC = ISYMAT(ICQROP(IOPER))

         ISAPROP = ISAMA * ISAMB * ISAMC

         IDISP   = 1

         ICTOT   = IQCAUA(IDISP) + IQCAUB(IDISP) + IQCAUC(IDISP)
         ISACAU  = 2 * ( (ICTOT/2)*2 - ICTOT ) + 1
         TSTISA  = ISACAU.EQ.ISAPROP .OR. ISAPROP.EQ.0 .OR. ALLDSPCF

         IF (.NOT. TSTISA) IDISP = IDISP + 1

         IF (MULD2H(ISYMA,ISYMB).EQ.ISYMC) THEN
            WRITE(LUPRI,'(/1X,3(A8,I3,3X),G16.8)')
     &        LBLOPR(IAQROP(IOPER)),IQCAUA(IDISP),
     &        LBLOPR(IBQROP(IOPER)),IQCAUB(IDISP),
     &        LBLOPR(ICQROP(IOPER)),IQCAUC(IDISP),
     &        -EXPCOF(IDISP,IOPER)
         ELSE
           IF (IPRQHYP.GT.0) THEN
            WRITE(LUPRI,'(/1X,3(A8,A3,3X),6X,A,7X)')
     &        LBLOPR(IAQROP(IOPER)),' - ',
     &        LBLOPR(IBQROP(IOPER)),' - ',
     &        LBLOPR(ICQROP(IOPER)),' - ', 
     &        '---'
           END IF
         END IF

         IDISP2 = IDISP + 1

         DO IDISP = IDISP2, NQRDISP
          ICTOT   = IQCAUA(IDISP) + IQCAUB(IDISP) + IQCAUC(IDISP)
          ISACAU  = 2 * ( (ICTOT/2)*2 - ICTOT ) + 1
          TSTISA  = ISACAU.EQ.ISAPROP .OR. ISAPROP.EQ.0 .OR. ALLDSPCF

          IF (MULD2H(ISYMA,ISYMB).EQ.ISYMC .AND. TSTISA) THEN
            WRITE(LUPRI,'(1X,3(8X,I3,3X),G16.8)')
     &        IQCAUA(IDISP), IQCAUB(IDISP), IQCAUC(IDISP),
     &        -EXPCOF(IDISP,IOPER)
          END IF
         END DO

      END DO

      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCQEXPPRT                            *
*---------------------------------------------------------------------*
c /* deck ccqdisprt */
*=====================================================================*
       SUBROUTINE CCQDISPRT (DISPCF, AVEDISP,ABCDE,
     &                       SHGCF, ORCOF, EOPECF,
     &                       AVESHG,AVEOR,AVEEOPE)
*---------------------------------------------------------------------*
*
*   Purpose: print dispersion coefficients first hyperpolarizabilities 
*
*
*   Written by Christof Haettig in November 1997.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccqrinf.h>
#include <ccroper.h>

      INTEGER KA, KB, KC, KCP, KD, KDP, KE, KEP
      PARAMETER (KA=1,KB=2,KC=3,KCP=4,KD=5,KDP=6,KE=7,KEP=8)
      INTEGER KBP, KDP2, KEP2, KEP3
      PARAMETER (KBP=9,KDP2=10,KEP2=11,KEP3=12)

      CHARACTER*5  BLANKS
      CHARACTER*80 STRING
      INTEGER ISYMA, ISYMB, ISYMC, ISAMA, ISAMB, ISAMC, ISAPROP
      INTEGER IOPER, MN, N, M, MN0, MNS


#if defined (SYS_CRAY)
      REAL DISPCF((NQRDSPE+2)*(NQRDSPE+1)/2,NQROPER)
      REAL AVEDISP((NQRDSPE+2)*(NQRDSPE+1)/2,4)
      REAL ABCDE(12,2)
      REAL SHGCF(NQRDSPE+1,NQROPER)
      REAL ORCOF(NQRDSPE+1,NQROPER)
      REAL EOPECF(NQRDSPE+1,NQROPER)
      REAL AVESHG(NQRDSPE+1,4)
      REAL AVEOR(NQRDSPE+1,4)
      REAL AVEEOPE(NQRDSPE+1,4)
#else
      DOUBLE PRECISION DISPCF((NQRDSPE+2)*(NQRDSPE+1)/2,NQROPER)
      DOUBLE PRECISION AVEDISP((NQRDSPE+2)*(NQRDSPE+1)/2,4)
      DOUBLE PRECISION ABCDE(12,2)
      DOUBLE PRECISION SHGCF(NQRDSPE+1,NQROPER)
      DOUBLE PRECISION ORCOF(NQRDSPE+1,NQROPER)
      DOUBLE PRECISION EOPECF(NQRDSPE+1,NQROPER)
      DOUBLE PRECISION AVESHG(NQRDSPE+1,4)
      DOUBLE PRECISION AVEOR(NQRDSPE+1,4)
      DOUBLE PRECISION AVEEOPE(NQRDSPE+1,4)
#endif

C
      INTEGER ITRI, I, J
      ITRI(I,J) = (MAX(I,J)+1)*(MAX(I,J)-2)/2 + I + J + 2
C

*---------------------------------------------------------------------*
* print header for hyperpolarizability dispersion coefficient section
*---------------------------------------------------------------------*
      BLANKS = '     '
      STRING = ' RESULTS FOR DISPERSION COEFFICIENTS'

      IF (CCS) THEN
         CALL AROUND( BLANKS//'FINAL CCS'//STRING(1:45)//BLANKS ) 
      ELSE IF (CC2) THEN
         CALL AROUND( BLANKS//'FINAL CC2'//STRING(1:45)//BLANKS )
      ELSE IF (CCSD) THEN
         CALL AROUND( BLANKS//'FINAL CCSD'//STRING(1:45)//BLANKS )
      ELSE
         CALL QUIT('CCQDISPRT called for an unknown '//
     &        'Coupled Cluster model.')
      END IF

      WRITE(LUPRI,'(/1X,3(A," operator",3X),3(4X,A),/,68("-"))')
     &      'A','B','C','N','M','D_ABC'


      DO IOPER = 1, NQROPER
         ISYMA = ISYOPR(IAQROP(IOPER))
         ISYMB = ISYOPR(IBQROP(IOPER))
         ISYMC = ISYOPR(ICQROP(IOPER))

         ISAMA = ISYMAT(IAQROP(IOPER))
         ISAMB = ISYMAT(IBQROP(IOPER))
         ISAMC = ISYMAT(ICQROP(IOPER))

         ISAPROP = ISAMA * ISAMB * ISAMC
         IF (ALLDSPCF) ISAPROP = 0

         MN0 = 0
         MNS = 2

         IF (ISAPROP.EQ.-1) MN0 = 1
         IF (ISAPROP.EQ. 0) MNS = 1

         MN = MN0
         M  = 0
         N  = MN - M
         IF (MULD2H(ISYMA,ISYMB).EQ.ISYMC) THEN
            WRITE(LUPRI,'(/1X,3(A8,5X),2I5,1X,G16.8)')
     &        LBLOPR(IAQROP(IOPER)),
     &        LBLOPR(IBQROP(IOPER)),
     &        LBLOPR(ICQROP(IOPER)),N,M,
     &        -DISPCF(ITRI(M+N,N),IOPER)
         ELSE
           IF (IPRQHYP.GT.0) THEN
            WRITE(LUPRI,'(/1X,3(A8,5X),3A)')
     &        LBLOPR(IAQROP(IOPER)),
     &        LBLOPR(IBQROP(IOPER)),
     &        LBLOPR(ICQROP(IOPER)),'  -  ', '  -  ', '---'
           END IF
         END IF

         IF (MN0.EQ.1 .AND. MULD2H(ISYMA,ISYMB).EQ.ISYMC) THEN
           M = 1
           N = MN - M
           WRITE(LUPRI,'(1X,3(8X,5X),2I5,1X,G16.8)')
     &          N,M, -DISPCF(ITRI(M+N,N),IOPER)
         END IF

         IF (MULD2H(ISYMA,ISYMB).EQ.ISYMC) THEN
           DO MN = MN0+MNS, NQRDSPE, MNS
           DO M = 0, MN
             N = MN - M
             WRITE(LUPRI,'(1X,3(8X,5X),2I5,1X,G16.8)')
     &          N,M, -DISPCF(ITRI(M+N,N),IOPER)
           END DO
           END DO
         END IF

      END DO

      IF (BETA_AVERAGE) THEN
         MN = 0
         N  = 0
         M  = 0
         WRITE(LUPRI,'(/1X,A10,29X,2I5,1X,G16.8)')
     &     'beta_{||} ',N,M, -AVEDISP(ITRI(M+N,M),1)
         DO MN = 2, NQRDSPE, 2
         DO M = 0, MN
           N = MN - M
           WRITE(LUPRI,'(1X,3(8X,5X),2I5,1X,G16.8)')
     &        N,M, -AVEDISP(ITRI(M+N,M),1)
         END DO
         END DO

         MN = 0
         N  = 0
         M  = 0
         WRITE(LUPRI,'(/1X,A10,29X,2I5,1X,G16.8)')
     &     'beta^K    ',N,M, -AVEDISP(ITRI(M+N,M),4)
         DO MN = 2, NQRDSPE, 2
         DO M = 0, MN
           N = MN - M
           WRITE(LUPRI,'(1X,3(8X,5X),2I5,1X,G16.8)')
     &        N,M, -AVEDISP(ITRI(M+N,M),4)
         END DO
         END DO

         MN = 0
         N  = 0
         M  = 0
         WRITE(LUPRI,'(/1X,A10,29X,2I5,1X,G16.8)')
     &     'beta_{_|_}',N,M, -AVEDISP(ITRI(M+N,M),2)
         DO MN = 2, NQRDSPE, 2
         DO M = 0, MN
           N = MN - M
           WRITE(LUPRI,'(1X,3(8X,5X),2I5,1X,G16.8)')
     &        N,M, -AVEDISP(ITRI(M+N,M),2)
         END DO
         END DO

         MN = 0
         N  = 0
         M  = 0
         WRITE(LUPRI,'(/1X,A10,29X,2I5,1X,G16.8)')
     &     'beta_ms   ',N,M, -AVEDISP(ITRI(M+N,M),3)
         DO MN = 2, NQRDSPE, 2
         DO M = 0, MN
           N = MN - M
           WRITE(LUPRI,'(1X,3(8X,5X),2I5,1X,G16.8)')
     &        N,M, -AVEDISP(ITRI(M+N,M),3)
         END DO
         END DO
      END IF

      WRITE(LUPRI,'(68("-"),//)')

      WRITE(LUPRI,'(////1X,3(A," operator",3X),2(4X,A),/,91("-"))')
     &    'A','B','C','N','D^SHG           D^OR            D^EOPE'

      DO IOPER = 1, NQROPER
         ISYMA = ISYOPR(IAQROP(IOPER))
         ISYMB = ISYOPR(IBQROP(IOPER))
         ISYMC = ISYOPR(ICQROP(IOPER))

         ISAMA = ISYMAT(IAQROP(IOPER))
         ISAMB = ISYMAT(IBQROP(IOPER))
         ISAMC = ISYMAT(ICQROP(IOPER))

         ISAPROP = ISAMA * ISAMB * ISAMC
         IF (ALLDSPCF) ISAPROP = 0

         MN0 = 0
         MNS = 2

         IF (ISAPROP.EQ.-1) MN0 = 1
         IF (ISAPROP.EQ. 0) MNS = 1

         MN = MN0
         IF (MULD2H(ISYMA,ISYMB).EQ.ISYMC) THEN
            WRITE(LUPRI,'(/1X,3(A8,5X),I5,1X,3G16.8)')
     &        LBLOPR(IAQROP(IOPER)),
     &        LBLOPR(IBQROP(IOPER)),
     &        LBLOPR(ICQROP(IOPER)),MN,
     &        -SHGCF(MN+1,IOPER),
     &        -ORCOF(MN+1,IOPER),
     &        -EOPECF(MN+1,IOPER)
         ELSE
           IF (IPRQHYP.GT.0) THEN
            WRITE(LUPRI,'(/1X,3(A8,5X),A,3(A,8X))')
     &        LBLOPR(IAQROP(IOPER)),
     &        LBLOPR(IBQROP(IOPER)),
     &        LBLOPR(ICQROP(IOPER)),'  -  ',
     &        '---','---','---'
           END IF
         END IF

         IF (MULD2H(ISYMA,ISYMB).EQ.ISYMC) THEN
           DO MN = MN0+MNS, NQRDSPE, MNS
             WRITE(LUPRI,'(1X,3(8X,5X),I5,1X,3G16.8)') MN,
     &        -SHGCF(MN+1,IOPER),-ORCOF(MN+1,IOPER),-EOPECF(MN+1,IOPER)
           END DO
         END IF

      END DO

      IF (BETA_AVERAGE) THEN
        MN  = 0
        WRITE(LUPRI,'(/1X,A10,29X,I5,1X,3G16.8)') 'beta_{||} ',MN,
     &    -AVESHG(MN+1,1), -AVEOR(MN+1,1), -AVEEOPE(MN+1,1)
        DO MN = 2, NQRDSPE, 2
          WRITE(LUPRI,'(1X,39X,I5,1X,3G16.8)') MN,
     &      -AVESHG(MN+1,1), -AVEOR(MN+1,1), -AVEEOPE(MN+1,1)
        END DO
        MN  = 0
        WRITE(LUPRI,'(/1X,A10,29X,I5,1X,3G16.8)') 'beta^K    ',MN,
     &    -AVESHG(MN+1,4), -AVEOR(MN+1,4), -AVEEOPE(MN+1,4)
        DO MN = 2, NQRDSPE, 2
          WRITE(LUPRI,'(1X,39X,I5,1X,3G16.8)') MN,
     &      -AVESHG(MN+1,4), -AVEOR(MN+1,4), -AVEEOPE(MN+1,4)
        END DO
        MN  = 0
        WRITE(LUPRI,'(/1X,A10,29X,I5,1X,3G16.8)') 'beta_{_|_}',MN,
     &    -AVESHG(MN+1,2), -AVEOR(MN+1,2), -AVEEOPE(MN+1,2)
        DO MN = 2, NQRDSPE, 2
          WRITE(LUPRI,'(1X,39X,I5,1X,3G16.8)') MN,
     &      -AVESHG(MN+1,2), -AVEOR(MN+1,2), -AVEEOPE(MN+1,2)
        END DO
        MN  = 0
        WRITE(LUPRI,'(/1X,A10,29X,I5,1X,3G16.8)') 'beta_ms   ',MN,
     &    -AVESHG(MN+1,3), -AVEOR(MN+1,3), -AVEEOPE(MN+1,3)
        DO MN = 2, NQRDSPE, 2
          WRITE(LUPRI,'(1X,39X,I5,1X,3G16.8)') MN,
     &      -AVESHG(MN+1,3), -AVEOR(MN+1,3), -AVEEOPE(MN+1,3)
        END DO
      END IF

      WRITE(LUPRI,'(91("-"),//)')

      IF (BETA_AVERAGE) THEN
        WRITE(LUPRI,'(//////1X,A,/,49("-"))')
     &   'PROCESS INDEPENDENT COEFFICIENTS FOR beta_{||}:'
        WRITE(LUPRI,'(1X,A,G16.8)')'beta(0)',-AVEDISP(1,1)
        IF (NQRDSPE.GE.2)WRITE(LUPRI,'(1X,A,5X,G16.8)')"A ",ABCDE(KA,1)
        IF (NQRDSPE.GE.4)WRITE(LUPRI,'(1X,A,5X,G16.8)')"B ",ABCDE(KB,1)
        IF (NQRDSPE.GE.6) THEN
           WRITE(LUPRI,'(1X,A,5X,G16.8,5X,A,5X,G16.8)')
     &       "C ",ABCDE(KC,1), "C'",ABCDE(KCP,1)
        END IF
        IF (NQRDSPE.GE.8) THEN
           WRITE(LUPRI,'(1X,A,5X,G16.8,5X,A,5X,G16.8)')
     &       "D ",ABCDE(KD,1), "D'",ABCDE(KDP,1)
        END IF
        IF (NQRDSPE.GE.10) THEN
           WRITE(LUPRI,'(1X,A,5X,G16.8,5X,A,5X,G16.8)')
     &       "E ",ABCDE(KE,1), "E'",ABCDE(KEP,1)
        END IF
        WRITE(LUPRI,'(49("-"),//)')
      END IF

      IF (BETA_AVERAGE) THEN
        WRITE(LUPRI,'(///1X,A,/,101("-"))')
     &   'PROCESS INDEPENDENT COEFFICIENTS FOR beta_ms:'
        IF (NQRDSPE.GE.2)
     &     WRITE(LUPRI,'(1X,A,2X,G16.8)') "A   ",ABCDE(KA,2)
        IF (NQRDSPE.GE.4) THEN
           WRITE(LUPRI,'(1X,A,2X,G16.8,5X,A,2X,G16.8)')
     &       "B   ",ABCDE(KB,2), "B'  ",ABCDE(KBP,2)
        END IF
        IF (NQRDSPE.GE.6) THEN
           WRITE(LUPRI,'(1X,A,2X,G16.8,5X,A,2X,G16.8)')
     &       "C   ",ABCDE(KC,2), "C'  ",ABCDE(KCP,2)
        END IF
        IF (NQRDSPE.GE.8) THEN
           WRITE(LUPRI,'(1X,3(A,2X,G16.8,5X))')
     &      "D   ",ABCDE(KD,2),"D'  ",ABCDE(KDP,2),"D'' ",ABCDE(KDP2,2)
        END IF
        IF (NQRDSPE.GE.10) THEN
           WRITE(LUPRI,'(1X,3(A,2X,G16.8,5X),A,2X,G16.8)')
     &       "E   ",ABCDE(KE,2),   "E'  ",ABCDE(KEP,2),
     &       "E'' ",ABCDE(KEP2,2), "E'''",ABCDE(KEP3,2)
        END IF
        WRITE(LUPRI,'(101("-"),////)')
      END IF


      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCQDISPRT                            *
*---------------------------------------------------------------------*
C  /* Deck around120 */
      SUBROUTINE AROUND120(HEAD)
      CHARACTER HEAD*(*)
#include <priunit.h>
      LNG = LEN(HEAD) + 2
      IND = (120 - LNG)/2 + 1
      WRITE (LUPRI, '(//,120A)') (' ',I=1,IND), '+', ('-',I=1,LNG), '+'
      WRITE (LUPRI, '(120A)')    (' ',I=1,IND), '! ', HEAD, ' !'
      WRITE (LUPRI, '(120A)')    (' ',I=1,IND), '+', ('-',I=1,LNG), '+'
      WRITE (LUPRI, '()')
      CALL FLSHFO(LUPRI)
      RETURN
      END
*---------------------------------------------------------------------*
