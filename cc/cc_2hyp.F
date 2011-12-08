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
c /* deck CC_2HYP */
*=====================================================================*
       SUBROUTINE CC_2HYP(WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: direct calculation of frequency-dependent second
*             hyperpolarizabilities for the Couple Cluster models
*
*                        CCS, CC2, CCSD
*
*     Written by Christof Haettig, february 1997.
*     Dispersion coefficients, march 1998, Christof Haettig
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
#include <cccrinf.h>
#include <ccroper.h>
#include <ccr1rsp.h>
#include <ccrc1rsp.h>
#include <cccperm.h>
#include <second.h>
#include <cclists.h>
#include <dummy.h>

* local parameters:
      CHARACTER*(19) MSGDBG
      PARAMETER (MSGDBG = '[debug] CC_2HYP> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

#if defined (SYS_CRAY)
      REAL ZERO
#else
      DOUBLE PRECISION ZERO
#endif
      PARAMETER (ZERO = 0.0d0)

      LOGICAL LCHI2_SAVE, LXKS3_SAVE, LADD, LERROR
      CHARACTER*10 MODEL, MOPRPC
      INTEGER LWORK
      INTEGER IZT(4),IAM(4),IOP(4),IR2(6),IX2(6),IO3(4),IL2(6),IO2(6)
      INTEGER ISYML, ISYM1, ISYM2, ISYM3, ISYMA, ISYMB, ISYMC, ISYMD
      INTEGER KHYPPOL, NBHYPPOL, IKAP(4)
      INTEGER ITRAN, IVEC, IDX, I, ISIGN, IFREQ, IOPER, IOPT
      INTEGER K0HTRAN, K0HDOTS, K0HCONS, KXCONS, KXDOTS, KXTRAN
      INTEGER K0GTRAN, K0GDOTS, K0GCONS, KAGTRAN, KAGDOTS, KAGCONS
      INTEGER K0FTRAN, K0FDOTS, K0FCONS, KAFTRAN, KAFDOTS, KAFCONS
      INTEGER K0FATRAN,K0FADOTS,K0FACONS,KAFATRAN,KAFADOTS,KAFACONS
      INTEGER K0EATRAN,K0EADOTS,K0EACONS,KAEATRAN,KAEADOTS,KAEACONS
      INTEGER KOCONS, KODOTS, KOTRAN, KLCONS, KLDOTS, KLTRAN, KEXPCOF
      INTEGER N0HTRAN, N0GTRAN, NAGTRAN, N0FTRAN, NAFTRAN, NBEXPCOF
      INTEGER N0FATRAN, NAFATRAN, NAEATRAN, NXTRAN, NOTRAN, NLTRAN
      INTEGER MX0HTRAN, MX0GTRAN, MXAGTRAN, MX0FTRAN, MXAFTRAN
      INTEGER MX0FATRAN, MXAFATRAN, MXAEATRAN, MXXTRAN, MXOTRAN,MXLTRAN
      INTEGER MX0HDOTS, MX0GDOTS, MXAGDOTS, MX0FDOTS, MXAFDOTS
      INTEGER MX0FADOTS, MXAFADOTS, MXAEADOTS, MXXDOTS, MXODOTS,MXLDOTS
      INTEGER MXTRAN3, MXTRAN2, MXVEC1, MXVEC2, IORDER
      INTEGER KEND0, LEND0
      INTEGER P
      INTEGER NTRINOM, KTRINOM, KEND0A, KEND1, LEND1
      INTEGER KDSPCF, KTHGCF, KESHGCF, KKERRCF, KDFWMCF
      INTEGER KDSPCFA, KTHGCFA, KESHGCA, KKERRCA, KDFWMCA, KABCDE


#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL TIM0, TIM1, TIMH, TIMG, TIMF, TIMFA, TIMEA, TIMX2, TIMO2
      REAL TIM2, TIML2
      REAL H0CON, G0CON(6), GACON(4), F0CON(3), FACON(12)
      REAL F0ACON(12), FAACON(12), EAACON(12), XCON(6)
      REAL OCON(4), LCON(3)
      REAL SIGN
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION TIM0,TIM1,TIMH,TIMG,TIMF,TIMFA,TIMEA,TIMX2,TIMO2
      DOUBLE PRECISION TIM2,TIML2
      DOUBLE PRECISION H0CON, G0CON(6), GACON(4), F0CON(3), FACON(12)
      DOUBLE PRECISION F0ACON(12), FAACON(12), EAACON(12), XCON(6)
      DOUBLE PRECISION OCON(4), LCON(3)
      DOUBLE PRECISION SIGN
#endif

* external functions
      INTEGER IR3TAMP
      INTEGER IR2TAMP
      INTEGER IR1TAMP
      INTEGER IL1ZETA
      INTEGER IL2ZETA
      INTEGER ICHI2
      INTEGER IRHSR2
      INTEGER IRHSR3


*---------------------------------------------------------------------*
* print header for hyperpolarizability section
*---------------------------------------------------------------------*
      WRITE (LUPRI,'(7(/1X,2A),/)')
     & '************************************',
     &                               '*******************************',
     & '*                                   ',
     &                               '                              *',
     & '*<<<<<<<<<  OUTPUT FROM COUPLED CLUST',
     &                                'ER CUBIC RESPONSE   >>>>>>>>>*',
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
         CALL QUIT('CC_2HYP called for unknown Coupled Cluster.')
      END IF

* print some debug/info output
      IF (IPRCHYP .GT. 10) WRITE(LUPRI,*) 'CC_2HYP Workspace:',LWORK
  
      TIM0 = SECOND()

      IF (NCRFREQ.EQ.0 .AND. NCRDISP.EQ.0) THEN
        WRITE(LUPRI,*) 'Nothing to be done in this section...'
        RETURN
      END IF

      IF (NCRFREQ.GT.0) THEN
*---------------------------------------------------------------------*
* allocate & initialize work space for hyperpolarizabilities
*---------------------------------------------------------------------*
      NBHYPPOL = NCROPER * NCRFREQ

      KHYPPOL = 1
      KEND0   = KHYPPOL + 2*NBHYPPOL


      MXTRAN2 = 2 * NLRTLBL * NLRTLBL
      MXVEC2  = (NLRTLBL+1) * NLRTLBL /2
      MXTRAN3 = 2 * NLRTLBL * NLRTLBL * NLRTLBL
      MXVEC1  = NLRTLBL
      MXTRAN2 = MIN(MXTRAN2,2*12*NBHYPPOL) ! 12 = max. # of permut.
      MXTRAN3 = MIN(MXTRAN3,2*12*NBHYPPOL) ! 12 = max. # of permut.
      MXVEC2  = MIN(MXVEC2,2*12*NBHYPPOL)  ! 12 = max. # of permut.
      MXVEC1  = MIN(MXVEC1,2*12*NBHYPPOL)  ! 12 = max. # of permut.

      MX0HTRAN  = MXDIM_HTRAN  * MXTRAN3
      MX0GTRAN  = MXDIM_GTRAN  * MXTRAN2
      MXAGTRAN  = MXDIM_GTRAN  * MXTRAN3
      MX0FTRAN  = MXDIM_FTRAN  * MXTRAN2
      MXAFTRAN  = MXDIM_FTRAN  * MXTRAN2
      MX0FATRAN = MXDIM_FATRAN * MXTRAN2
      MXAFATRAN = MXDIM_FATRAN * MXTRAN3
      MXAEATRAN = MXDIM_XEVEC  * MXTRAN2
      MXXTRAN   = 1 * MXTRAN2
      MXOTRAN   = 1 * MXTRAN3
      MXLTRAN   = 1 * MXTRAN2

      MX0HDOTS  = MXTRAN3 * MXVEC1
      MX0GDOTS  = MXTRAN2 * MXVEC2
      MXAGDOTS  = MXTRAN3 * MXVEC1
      MX0FDOTS  = MXTRAN2 * MXVEC2
      MXAFDOTS  = MXTRAN2 * MXVEC2
      MX0FADOTS = MXTRAN2 * MXVEC2
      MXAFADOTS = MXTRAN3 * MXVEC1
      MXAEADOTS = MXTRAN2 * MXVEC2
      MXXDOTS   = MXTRAN2 * MXVEC2
      MXODOTS   = MXTRAN3 * MXVEC1
      MXLDOTS   = MXTRAN2 * MXVEC2

      K0HTRAN  = KEND0
      K0HDOTS  = K0HTRAN  + MX0HTRAN
      K0HCONS  = K0HDOTS  + MX0HDOTS
      KEND0    = K0HCONS  + MX0HDOTS

      K0GTRAN  = KEND0
      K0GDOTS  = K0GTRAN  + MX0GTRAN
      K0GCONS  = K0GDOTS  + MX0GDOTS
      KEND0    = K0GCONS  + MX0GDOTS

      KAGTRAN  = KEND0
      KAGDOTS  = KAGTRAN  + MXAGTRAN
      KAGCONS  = KAGDOTS  + MXAGDOTS
      KEND0    = KAGCONS  + MXAGDOTS

      K0FTRAN  = KEND0
      K0FDOTS  = K0FTRAN  + MX0FTRAN
      K0FCONS  = K0FDOTS  + MX0FDOTS
      KEND0    = K0FCONS  + MX0FDOTS

      KAFTRAN  = KEND0
      KAFDOTS  = KAFTRAN  + MXAFTRAN
      KAFCONS  = KAFDOTS  + MXAFDOTS
      KEND0    = KAFCONS  + MXAFDOTS

      K0FATRAN = KEND0
      K0FADOTS = K0FATRAN + MX0FATRAN
      K0FACONS = K0FADOTS + MX0FADOTS
      KEND0    = K0FACONS + MX0FADOTS

      KAFATRAN = KEND0
      KAFADOTS = KAFATRAN + MXAFATRAN
      KAFACONS = KAFADOTS + MXAFADOTS
      KEND0    = KAFACONS + MXAFADOTS

      KAEATRAN = KEND0
      KAEADOTS = KAEATRAN + MXAEATRAN
      KAEACONS = KAEADOTS + MXAEADOTS
      KEND0    = KAEACONS + MXAEADOTS

      KXTRAN   = KEND0
      KXDOTS   = KXTRAN   + MXXTRAN
      KXCONS   = KXDOTS   + MXXDOTS
      KEND0    = KXCONS   + MXXDOTS

      KOTRAN   = KEND0
      KODOTS   = KOTRAN   + MXOTRAN
      KOCONS   = KODOTS   + MXODOTS
      KEND0    = KOCONS   + MXODOTS

      KLTRAN   = KEND0
      KLDOTS   = KLTRAN   + MXLTRAN
      KLCONS   = KLDOTS   + MXLDOTS
      KEND0    = KLCONS   + MXLDOTS

      LEND0    = LWORK - KEND0

      IF (LEND0 .LT. 0) THEN
       CALL QUIT('Insufficient memory in CC_2HYP. (1)')
      END IF

* initialize hyperpolarizabilities and the lists of contributions:
      CALL DZERO(WORK(KHYPPOL),2*NBHYPPOL)
      CALL DZERO(WORK(K0HTRAN),MX0HTRAN+2*MX0HDOTS)
      CALL DZERO(WORK(K0GTRAN),MX0GTRAN+2*MX0GDOTS)
      CALL DZERO(WORK(KAGTRAN),MXAGTRAN+2*MXAGDOTS)
      CALL DZERO(WORK(K0FTRAN),MX0FTRAN+2*MX0FDOTS)
      CALL DZERO(WORK(KAFTRAN),MXAFTRAN+2*MXAFDOTS)
      CALL DZERO(WORK(K0FATRAN),MX0FATRAN+2*MX0FADOTS)
      CALL DZERO(WORK(KAFATRAN),MXAFATRAN+2*MXAFADOTS)
      CALL DZERO(WORK(KAEATRAN),MXAEATRAN+2*MXAEADOTS)
      CALL DZERO(WORK(KXTRAN),MXXTRAN+2*MXXDOTS)
      CALL DZERO(WORK(KOTRAN),MXOTRAN+2*MXODOTS)
      CALL DZERO(WORK(KLTRAN),MXLTRAN+2*MXLDOTS)

*---------------------------------------------------------------------*
* set up lists for H, G, F and F{O} transformations and ETA{O} vectors:
*---------------------------------------------------------------------*
      CALL CCCR_SETUP(MXTRAN2, MXVEC2, MXTRAN3, MXVEC1,
     &                WORK(K0HTRAN), WORK(K0HDOTS), N0HTRAN,
     &                WORK(K0GTRAN), WORK(K0GDOTS), N0GTRAN,
     &                WORK(KAGTRAN), WORK(KAGDOTS), NAGTRAN,
     &                WORK(K0FTRAN), WORK(K0FDOTS), N0FTRAN,
     &                WORK(KAFTRAN), WORK(KAFDOTS), NAFTRAN,
     &                WORK(K0FATRAN),WORK(K0FADOTS),N0FATRAN,
     &                WORK(KAFATRAN),WORK(KAFADOTS),NAFATRAN,
     &                WORK(KAEATRAN),WORK(KAEADOTS),NAEATRAN,
     &                WORK(KXTRAN),  WORK(KXDOTS),  NXTRAN,    
     &                WORK(KOTRAN),  WORK(KODOTS),  NOTRAN,
     &                WORK(KLTRAN),  WORK(KLDOTS),  NLTRAN    )

*---------------------------------------------------------------------*
* calculate H matrix contributions:
*---------------------------------------------------------------------*
      TIM1 = SECOND()

      IOPT = 5
      CALL CC_HMAT('L0 ','R1 ','R1 ','R1 ','R1 ',N0HTRAN, MXVEC1,
     &             WORK(K0HTRAN),WORK(K0HDOTS),WORK(K0HCONS),
     &             WORK(KEND0), LEND0, IOPT )

      TIMH = SECOND() - TIM1
      WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")')
     & '>>> Time used for',N0HTRAN,' H matrix transformations:   ',TIMH
      CALL FLSHFO(LUPRI)

*---------------------------------------------------------------------*
* calculate G matrix contributions:
*---------------------------------------------------------------------*
      TIM1 = SECOND()

      IF (.NOT.L_USE_CHI2) THEN
       IOPT = 5
       CALL CC_GMATRIX('L0 ','R1 ','R1 ','R2 ',N0GTRAN, MXVEC2,
     &               WORK(K0GTRAN),WORK(K0GDOTS),WORK(K0GCONS),
     &               WORK(KEND0), LEND0, IOPT )
      END IF

      IF (.NOT.L_USE_XKS3) THEN
        IOPT = 5
        CALL CC_GMATRIX('L1 ','R1 ','R1 ','R1 ',NAGTRAN, MXVEC1,
     &                WORK(KAGTRAN),WORK(KAGDOTS),WORK(KAGCONS),
     &                WORK(KEND0), LEND0, IOPT )
      END IF
     
      TIMG = SECOND() - TIM1
      WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")')
     &  '>>> Time used for',N0GTRAN+NAGTRAN,
     &           ' G matrix transformations:   ',TIMG
      CALL FLSHFO(LUPRI)
      
*---------------------------------------------------------------------*
* calculate F matrix contributions:
*---------------------------------------------------------------------*
      TIM1 = SECOND()

      IOPT = 5
      CALL CC_FMATRIX(WORK(K0FTRAN),N0FTRAN,'L0 ','R2 ',IOPT,'R2 ',
     &                WORK(K0FDOTS),WORK(K0FCONS),MXVEC2,
     &                WORK(KEND0),LEND0)
      IF (L_USE_CHI2) THEN
        CALL DSCAL(MX0FDOTS,-1.0d0,WORK(K0FCONS),1)
      END IF

      IF (.NOT. (L_USE_CHI2 .OR. L_USE_XKS3) ) THEN
         IOPT = 5
         CALL CC_FMATRIX(WORK(KAFTRAN),NAFTRAN,'L1 ','R1 ',IOPT,'R2 ',
     &                   WORK(KAFDOTS),WORK(KAFCONS),MXVEC2,
     &                   WORK(KEND0),LEND0)
      END IF

      TIMF = SECOND() - TIM1
      WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")')
     &  '>>> Time used for',N0FTRAN+NAFTRAN,
     &          ' F matrix transformations:   ',TIMF
      CALL FLSHFO(LUPRI)
      
*---------------------------------------------------------------------*
* calculate F{O} matrix contributions:
*---------------------------------------------------------------------*
      TIM1 = SECOND()

      IF (.NOT.L_USE_CHI2) THEN
       CALL CCQR_FADRV('L0 ','o1 ','R1 ','R2 ',N0FATRAN, MXVEC2,
     &                  WORK(K0FATRAN),WORK(K0FADOTS),
     &                  WORK(K0FACONS),
     &                  WORK(KEND0), LEND0, 'DOTP' )
      END IF

      IF (.NOT.L_USE_XKS3) THEN
        CALL CCQR_FADRV('L1 ','o1 ','R1 ','R1 ',NAFATRAN, MXVEC1,
     &                   WORK(KAFATRAN),WORK(KAFADOTS),
     &                   WORK(KAFACONS),
     &                   WORK(KEND0), LEND0, 'DOTP' )
      END IF

      TIMFA = SECOND() - TIM1
      WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")')
     &  '>>> Time used for',N0FATRAN+NAFATRAN,
     &          ' F{O} matrix transformations:',TIMFA
      CALL FLSHFO(LUPRI)
      
*---------------------------------------------------------------------*
* calculate ETA{O} vector contributions:
*---------------------------------------------------------------------*
      TIM1 = SECOND()

      IF (.NOT. (L_USE_CHI2 .OR. L_USE_XKS3) ) THEN

       IOPT   = 5
       IORDER = 1
       CALL CC_XIETA( WORK(KAEATRAN), NAEATRAN, IOPT, IORDER, 'L1 ',
     &                '---',IDUMMY,       DUMMY,
     &                'R2 ',WORK(KAEADOTS),WORK(KAEACONS),
     &                .FALSE.,MXVEC2, WORK(KEND0), LEND0 )

      END IF

      TIMEA = SECOND() - TIM1
      WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")') '>>> Time used for',
     &  NAEATRAN,' ETA{O} vector calculations :',TIMEA
      CALL FLSHFO(LUPRI)
      
*---------------------------------------------------------------------*
* chi vector contributions:
*---------------------------------------------------------------------*
      IF (L_USE_CHI2) THEN
        TIM1 = SECOND()

        CALL CC_DOTDRV('X2B','R2 ',NXTRAN,MXVEC2,
     &                 WORK(KXTRAN), WORK(KXDOTS), WORK(KXCONS),
     &                 WORK(KEND0), LEND0 )

        TIMX2 = SECOND() - TIM1
        WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")')
     &    '>>> Time used for',NXTRAN, ' X2 x R2 dot products:',TIMX2
        CALL FLSHFO(LUPRI)
      
      END IF

*---------------------------------------------------------------------*
* L2 x O2 dot products:
*---------------------------------------------------------------------*
      IF (USE_L2BC .OR. USE_LBCD) THEN
        TIM1 = SECOND()

        CALL CC_DOTDRV('L2 ','O2 ',NLTRAN,MXVEC2,
     &                 WORK(KLTRAN), WORK(KLDOTS), WORK(KLCONS),
     &                 WORK(KEND0), LEND0 )

        TIML2 = SECOND() - TIM1
        WRITE (LUPRI,'(/A,I5,A,5X,F12.2," seconds.")')
     &    '>>> Time used for',NLTRAN, ' L2 x R2 dot products:  ',TIML2
        CALL FLSHFO(LUPRI)
      
      END IF
*---------------------------------------------------------------------*
* xksi vector contributions:
*---------------------------------------------------------------------*
      IF (L_USE_XKS3) THEN

        TIM1 = SECOND()

        CALL CC_DOTDRV('R3 ','X1B',NOTRAN,MXVEC1,
     &                 WORK(KOTRAN), WORK(KODOTS), WORK(KOCONS),
     &                 WORK(KEND0), LEND0 )

        TIMO2 = SECOND() - TIM1
        WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")')
     &    '>>> Time used for',NOTRAN, ' R3 x X1 dot products:',TIMO2
        CALL FLSHFO(LUPRI)

      END IF
*=====================================================================*
* loop over quadruples of operator labels and over frequencies and 
* collect all contributions to the final hyperpolarizabilities:
*=====================================================================*
      DO IOPER = 1, NCROPER

        IOP(A) = IACROP(IOPER)
        IOP(B) = IBCROP(IOPER)
        IOP(C) = ICCROP(IOPER)
        IOP(D) = IDCROP(IOPER)

        IKAP(A)= 0
        IKAP(B)= 0
        IKAP(C)= 0
        IKAP(D)= 0

        ISYMB  = ISYOPR(IOP(B))
        ISYMC  = ISYOPR(IOP(C))
        ISYMA  = ISYOPR(IOP(A))
        ISYMD  = ISYOPR(IOP(D))


      IF (MULD2H(ISYMA,ISYMB) .EQ. MULD2H(ISYMC,ISYMD)) THEN

      DO IFREQ = 1, NCRFREQ

      DO ISIGN = +1, -1, -2
        SIGN = DBLE(ISIGN)

        IZT(A) =IL1ZETA(LBLOPR(IOP(A)),.FALSE.,SIGN*ACRFR(IFREQ),ISYML)
        IZT(B) =IL1ZETA(LBLOPR(IOP(B)),.FALSE.,SIGN*BCRFR(IFREQ),ISYML)
        IZT(C) =IL1ZETA(LBLOPR(IOP(C)),.FALSE.,SIGN*CCRFR(IFREQ),ISYML)
        IZT(D) =IL1ZETA(LBLOPR(IOP(D)),.FALSE.,SIGN*DCRFR(IFREQ),ISYML)

        IAM(A) =IR1TAMP(LBLOPR(IOP(A)),.FALSE.,SIGN*ACRFR(IFREQ),ISYML)
        IAM(B) =IR1TAMP(LBLOPR(IOP(B)),.FALSE.,SIGN*BCRFR(IFREQ),ISYML)
        IAM(C) =IR1TAMP(LBLOPR(IOP(C)),.FALSE.,SIGN*CCRFR(IFREQ),ISYML)
        IAM(D) =IR1TAMP(LBLOPR(IOP(D)),.FALSE.,SIGN*DCRFR(IFREQ),ISYML)

        IF ( USE_LBCD ) THEN
         IL2(BC)=IL2ZETA(LBLOPR(IOP(B)),SIGN*BCRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(C)),SIGN*CCRFR(IFREQ),ISYM2)
         IL2(BD)=IL2ZETA(LBLOPR(IOP(B)),SIGN*BCRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(D)),SIGN*DCRFR(IFREQ),ISYM2)
         IL2(CD)=IL2ZETA(LBLOPR(IOP(C)),SIGN*CCRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(D)),SIGN*DCRFR(IFREQ),ISYM2)
         IO2(AB)=IRHSR2( LBLOPR(IOP(A)),.FALSE.,SIGN*ACRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(B)),.FALSE.,SIGN*BCRFR(IFREQ),ISYM2)
         IO2(AC)=IRHSR2( LBLOPR(IOP(A)),.FALSE.,SIGN*ACRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(C)),.FALSE.,SIGN*CCRFR(IFREQ),ISYM2)
         IO2(AD)=IRHSR2( LBLOPR(IOP(A)),.FALSE.,SIGN*ACRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(D)),.FALSE.,SIGN*DCRFR(IFREQ),ISYM2)
        ELSE IF ( USE_L2BC) THEN
         IR2(AB)=IR2TAMP(LBLOPR(IOP(A)),.FALSE.,SIGN*ACRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(B)),.FALSE.,SIGN*BCRFR(IFREQ),ISYM2)
         IR2(AC)=IR2TAMP(LBLOPR(IOP(A)),.FALSE.,SIGN*ACRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(C)),.FALSE.,SIGN*CCRFR(IFREQ),ISYM2)
         IL2(BC)=IL2ZETA(LBLOPR(IOP(B)),SIGN*BCRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(C)),SIGN*CCRFR(IFREQ),ISYM2)
         IO2(AD)=IRHSR2( LBLOPR(IOP(A)),.FALSE.,SIGN*ACRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(D)),.FALSE.,SIGN*DCRFR(IFREQ),ISYM2)
        ELSE
         IR2(AB)=IR2TAMP(LBLOPR(IOP(A)),.FALSE.,SIGN*ACRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(B)),.FALSE.,SIGN*BCRFR(IFREQ),ISYM2)
         IR2(AC)=IR2TAMP(LBLOPR(IOP(A)),.FALSE.,SIGN*ACRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(C)),.FALSE.,SIGN*CCRFR(IFREQ),ISYM2)
         IR2(AD)=IR2TAMP(LBLOPR(IOP(A)),.FALSE.,SIGN*ACRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(D)),.FALSE.,SIGN*DCRFR(IFREQ),ISYM2)
        END IF

        IR2(BC)=IR2TAMP(LBLOPR(IOP(B)),.FALSE.,SIGN*BCRFR(IFREQ),ISYM1,
     &                  LBLOPR(IOP(C)),.FALSE.,SIGN*CCRFR(IFREQ),ISYM2)
        IR2(BD)=IR2TAMP(LBLOPR(IOP(B)),.FALSE.,SIGN*BCRFR(IFREQ),ISYM1,
     &                  LBLOPR(IOP(D)),.FALSE.,SIGN*DCRFR(IFREQ),ISYM2)
        IR2(CD)=IR2TAMP(LBLOPR(IOP(C)),.FALSE.,SIGN*CCRFR(IFREQ),ISYM1,
     &                  LBLOPR(IOP(D)),.FALSE.,SIGN*DCRFR(IFREQ),ISYM2)

        IF (L_USE_CHI2) THEN
c         IX2(AB)=ICHI2(LBLOPR(IOP(A)),.FALSE.,SIGN*ACRFR(IFREQ),ISYM1,
c    &                  LBLOPR(IOP(B)),.FALSE.,SIGN*BCRFR(IFREQ),ISYM2)
c         IX2(AC)=ICHI2(LBLOPR(IOP(A)),.FALSE.,SIGN*ACRFR(IFREQ),ISYM1,
c    &                  LBLOPR(IOP(C)),.FALSE.,SIGN*CCRFR(IFREQ),ISYM2)
c         IX2(AD)=ICHI2(LBLOPR(IOP(A)),.FALSE.,SIGN*ACRFR(IFREQ),ISYM1,
c    &                  LBLOPR(IOP(D)),.FALSE.,SIGN*DCRFR(IFREQ),ISYM2)
c         IX2(BC)=ICHI2(LBLOPR(IOP(B)),.FALSE.,SIGN*BCRFR(IFREQ),ISYM1,
c    &                  LBLOPR(IOP(C)),.FALSE.,SIGN*CCRFR(IFREQ),ISYM2)
c         IX2(BD)=ICHI2(LBLOPR(IOP(B)),.FALSE.,SIGN*BCRFR(IFREQ),ISYM1,
c    &                  LBLOPR(IOP(D)),.FALSE.,SIGN*DCRFR(IFREQ),ISYM2)
c         IX2(CD)=ICHI2(LBLOPR(IOP(C)),.FALSE.,SIGN*CCRFR(IFREQ),ISYM1,
c    &                  LBLOPR(IOP(D)),.FALSE.,SIGN*DCRFR(IFREQ),ISYM2)

         IX2(AB)=IL2ZETA(LBLOPR(IOP(A)),        SIGN*ACRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(B)),        SIGN*BCRFR(IFREQ),ISYM2)
         IX2(AC)=IL2ZETA(LBLOPR(IOP(A)),        SIGN*ACRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(C)),        SIGN*CCRFR(IFREQ),ISYM2)
         IX2(AD)=IL2ZETA(LBLOPR(IOP(A)),        SIGN*ACRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(D)),        SIGN*DCRFR(IFREQ),ISYM2)
         IX2(BC)=IL2ZETA(LBLOPR(IOP(B)),        SIGN*BCRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(C)),        SIGN*CCRFR(IFREQ),ISYM2)
         IX2(BD)=IL2ZETA(LBLOPR(IOP(B)),        SIGN*BCRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(D)),        SIGN*DCRFR(IFREQ),ISYM2)
         IX2(CD)=IL2ZETA(LBLOPR(IOP(C)),        SIGN*CCRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(D)),        SIGN*DCRFR(IFREQ),ISYM2)
        END IF

        IF (L_USE_XKS3) THEN
          IO3(ABC) = IR3TAMP(LBLOPR(IOP(A)),SIGN*ACRFR(IFREQ),ISYM1,
     &                      LBLOPR(IOP(B)),SIGN*BCRFR(IFREQ),ISYM2,
     &                      LBLOPR(IOP(C)),SIGN*CCRFR(IFREQ),ISYM3)
          IO3(ABD) = IR3TAMP(LBLOPR(IOP(A)),SIGN*ACRFR(IFREQ),ISYM1,
     &                      LBLOPR(IOP(B)),SIGN*BCRFR(IFREQ),ISYM2,
     &                      LBLOPR(IOP(D)),SIGN*DCRFR(IFREQ),ISYM3)
          IO3(ACD) = IR3TAMP(LBLOPR(IOP(A)),SIGN*ACRFR(IFREQ),ISYM1,
     &                      LBLOPR(IOP(C)),SIGN*CCRFR(IFREQ),ISYM2,
     &                      LBLOPR(IOP(D)),SIGN*DCRFR(IFREQ),ISYM3)
          IO3(BCD) = IR3TAMP(LBLOPR(IOP(B)),SIGN*BCRFR(IFREQ),ISYM1,
     &                      LBLOPR(IOP(C)),SIGN*CCRFR(IFREQ),ISYM2,
     &                      LBLOPR(IOP(D)),SIGN*DCRFR(IFREQ),ISYM3)
        END IF


*---------------------------------------------------------------------*
* get H^0 matrix transformations, 1 permutation
*---------------------------------------------------------------------*
        CALL CC_SETH1111(WORK(K0HTRAN),WORK(K0HDOTS),MXTRAN3,MXVEC1,
     &                   0,IAM(A),IAM(B),IAM(C),IAM(D),ITRAN,IVEC)
        H0CON = WORK(K0HCONS-1 + (ITRAN-1)*MXVEC1 + IVEC)

*---------------------------------------------------------------------*
* get G^0 matrix transformations, 6 permutations
*---------------------------------------------------------------------*
      DO P = 1, 6
        IF ((USE_L2BC.AND.P.EQ.4) .OR. (USE_LBCD.AND.P.GT.3) ) THEN
          G0CON(P) = ZERO
        ELSE
          CALL CC_SETG112(WORK(K0GTRAN),WORK(K0GDOTS),MXTRAN2,MXVEC2,   
     &                0,IAM(I1(P)),IAM(I2(P)),IR2(IP(P)),ITRAN,IVEC)
          G0CON(P) = WORK(K0GCONS-1 + (ITRAN-1)*MXVEC2 + IVEC)
        END IF
      END DO

*---------------------------------------------------------------------*
* get G^A matrix transformations, 4 permutations
*---------------------------------------------------------------------*
        CALL CCQR_SETG(WORK(KAGTRAN),WORK(KAGDOTS),MXTRAN3,MXVEC1,  
     &                 IZT(A),IAM(B),IAM(C),IAM(D),ITRAN,IVEC)
        GACON(1) = WORK(KAGCONS-1 + (ITRAN-1)*MXVEC1 + IVEC)

        CALL CCQR_SETG(WORK(KAGTRAN),WORK(KAGDOTS),MXTRAN3,MXVEC1,
     &                 IZT(B),IAM(A),IAM(C),IAM(D),ITRAN,IVEC)
        GACON(2) = WORK(KAGCONS-1 + (ITRAN-1)*MXVEC1 + IVEC)

        CALL CCQR_SETG(WORK(KAGTRAN),WORK(KAGDOTS),MXTRAN3,MXVEC1,
     &                 IZT(C),IAM(B),IAM(A),IAM(D),ITRAN,IVEC)
        GACON(3) = WORK(KAGCONS-1 + (ITRAN-1)*MXVEC1 + IVEC)

        CALL CCQR_SETG(WORK(KAGTRAN),WORK(KAGDOTS),MXTRAN3,MXVEC1,
     &                 IZT(D),IAM(B),IAM(C),IAM(A),ITRAN,IVEC)
        GACON(4) = WORK(KAGCONS-1 + (ITRAN-1)*MXVEC1 + IVEC)

*---------------------------------------------------------------------*
* get F^0 matrix transformations, 3 permutations
*---------------------------------------------------------------------*
       IF (.NOT. USE_LBCD) THEN
          CALL CCQR_SETF(WORK(K0FTRAN),WORK(K0FDOTS),MXTRAN2,MXVEC2,
     &                   0,IR2(AB),IR2(CD),ITRAN,IVEC)
          F0CON(1) = WORK(K0FCONS-1 + (ITRAN-1)*MXVEC2 + IVEC)
       ELSE
          F0CON(1) = ZERO
       END IF

       IF (.NOT. USE_LBCD) THEN
          CALL CCQR_SETF(WORK(K0FTRAN),WORK(K0FDOTS),MXTRAN2,MXVEC2,
     &                   0,IR2(AC),IR2(BD),ITRAN,IVEC)
          F0CON(2) = WORK(K0FCONS-1 + (ITRAN-1)*MXVEC2 + IVEC)
       ELSE
          F0CON(2) = ZERO
       END IF

       IF (.NOT. (USE_LBCD .OR. USE_L2BC)) THEN
          CALL CCQR_SETF(WORK(K0FTRAN),WORK(K0FDOTS),MXTRAN2,MXVEC2,
     &                   0,IR2(AD),IR2(BC),ITRAN,IVEC)
          F0CON(3) = WORK(K0FCONS-1 + (ITRAN-1)*MXVEC2 + IVEC)
       ELSE
          F0CON(3) = ZERO
       END IF

*---------------------------------------------------------------------*
* get F^A matrix transformations, 12 permutations
*---------------------------------------------------------------------*
      DO P = 1, 6
        IF ((USE_L2BC.AND.P.EQ.4) .OR. (USE_LBCD.AND.P.GT.3) ) THEN
          FACON(P)   = ZERO
          FACON(6+P) = ZERO
        ELSE
          CALL CC_SETF12(WORK(KAFTRAN),WORK(KAFDOTS),MXTRAN2,MXVEC2, 
     &                   IZT(I1(P)),IAM(I2(P)),IR2(IP(P)),ITRAN,IVEC)
          FACON(P) = WORK(KAFCONS-1 + (ITRAN-1)*MXVEC2 + IVEC)

          CALL CC_SETF12(WORK(KAFTRAN),WORK(KAFDOTS),MXTRAN2,MXVEC2, 
     &                   IZT(I2(P)),IAM(I1(P)),IR2(IP(P)),ITRAN,IVEC)
          FACON(6+P) = WORK(KAFCONS-1 + (ITRAN-1)*MXVEC2 + IVEC)
        END IF
      END DO

*---------------------------------------------------------------------*
* get F^0{O} matrix transformations, 12 permutations
*---------------------------------------------------------------------*
      DO P = 1, 6
        IF ((USE_L2BC.AND.P.EQ.4) .OR. (USE_LBCD.AND.P.GT.3) ) THEN
          F0ACON(P)   = ZERO
          F0ACON(6+P) = ZERO
        ELSE
          CALL CC_SETFA12(WORK(K0FATRAN),WORK(K0FADOTS),MXTRAN2,MXVEC2, 
     &                 0,IOP(I1(P)),IAM(I2(P)),IR2(IP(P)),ITRAN,IVEC)
          F0ACON(P) = WORK(K0FACONS-1 + (ITRAN-1)*MXVEC2 + IVEC)

          CALL CC_SETFA12(WORK(K0FATRAN),WORK(K0FADOTS),MXTRAN2,MXVEC2, 
     &                 0,IOP(I2(P)),IAM(I1(P)),IR2(IP(P)),ITRAN,IVEC)
          F0ACON(6+P) = WORK(K0FACONS-1 + (ITRAN-1)*MXVEC2 + IVEC)
        END IF
      END DO

*---------------------------------------------------------------------*
* get F^A{O} matrix transformations, 12 permutations
*---------------------------------------------------------------------*
      DO P = 1, 6
        CALL CCQR_SETFA(WORK(KAFATRAN),WORK(KAFADOTS),MXTRAN3,MXVEC1, 
     &      IZT(I1(P)),IOP(I2(P)),IAM(I3(P)),IAM(I4(P)),ITRAN,IVEC)
        FAACON(P) = WORK(KAFACONS-1 + (ITRAN-1)*MXVEC1 + IVEC)

        CALL CCQR_SETFA(WORK(KAFATRAN),WORK(KAFADOTS),MXTRAN3,MXVEC1, 
     &      IZT(I2(P)),IOP(I1(P)),IAM(I3(P)),IAM(I4(P)),ITRAN,IVEC)
        FAACON(6+P) = WORK(KAFACONS-1 + (ITRAN-1)*MXVEC1 + IVEC)
      END DO

*---------------------------------------------------------------------*
* get ETA{O} vector calculations, 12 permutations
*---------------------------------------------------------------------*
      DO P = 1, 6
        IF ((USE_L2BC.AND.P.EQ.4) .OR. (USE_LBCD.AND.P.GT.3) ) THEN
          EAACON(P)   = ZERO
          EAACON(6+P) = ZERO
        ELSE
          CALL CC_SETXE('Eta',WORK(KAEATRAN),WORK(KAEADOTS), !1x2x3,4
     &                  MXTRAN2,MXVEC2,
     &                  IZT(I1(P)),IOP(I2(P)),IKAP(I2(P)),0,0,0,
     &                  IR2(IP(P)),ITRAN,IVEC)
          EAACON(P) = WORK(KAEACONS-1 + (ITRAN-1)*MXVEC2 + IVEC)

          CALL CC_SETXE('Eta',WORK(KAEATRAN),WORK(KAEADOTS), !2x1x3,4
     &                  MXTRAN2,MXVEC2,
     &                  IZT(I2(P)),IOP(I1(P)),IKAP(I1(P)),0,0,0,
     &                  IR2(IP(P)),ITRAN,IVEC)
          EAACON(6+P) = WORK(KAEACONS-1 + (ITRAN-1)*MXVEC2 + IVEC)
        END IF
      END DO

*---------------------------------------------------------------------*
* get L2 x O2 vector dot products, max. 3 permutations
*---------------------------------------------------------------------*
      IF (USE_LBCD .OR. USE_L2BC) THEN
        CALL CC_SETDOT(WORK(KLTRAN),WORK(KLDOTS),MXTRAN2,MXVEC2,
     &                 IL2(BC),IO2(AD),ITRAN,IVEC)
        LCON(1) = WORK(KLCONS-1 + (ITRAN-1)*MXVEC2 + IVEC)
      ELSE 
        LCON(1) = ZERO
      END IF

      IF (USE_LBCD) THEN
        CALL CC_SETDOT(WORK(KLTRAN),WORK(KLDOTS),MXTRAN2,MXVEC2,
     &                 IL2(BD),IO2(AC),ITRAN,IVEC)
        LCON(2) = WORK(KLCONS-1 + (ITRAN-1)*MXVEC2 + IVEC)

        CALL CC_SETDOT(WORK(KLTRAN),WORK(KLDOTS),MXTRAN2,MXVEC2,
     &                 IL2(CD),IO2(AB),ITRAN,IVEC)
        LCON(3) = WORK(KLCONS-1 + (ITRAN-1)*MXVEC2 + IVEC)
      ELSE 
        LCON(2) = ZERO
        LCON(3) = ZERO
      END IF
*---------------------------------------------------------------------*
* get chi vector dot products, 6 permutations
*---------------------------------------------------------------------*
      IF (L_USE_CHI2) THEN
        CALL CC_SETDOT(WORK(KXTRAN),WORK(KXDOTS),MXTRAN2,MXVEC2,
     &                 IX2(AB),IR2(CD),ITRAN,IVEC)
        XCON(1) = WORK(KXCONS-1 + (ITRAN-1)*MXVEC2 + IVEC)

        CALL CC_SETDOT(WORK(KXTRAN),WORK(KXDOTS),MXTRAN2,MXVEC2,
     &                 IX2(AC),IR2(BD),ITRAN,IVEC)
        XCON(2) = WORK(KXCONS-1 + (ITRAN-1)*MXVEC2 + IVEC)

        CALL CC_SETDOT(WORK(KXTRAN),WORK(KXDOTS),MXTRAN2,MXVEC2,
     &                 IX2(AD),IR2(BC),ITRAN,IVEC)
        XCON(3) = WORK(KXCONS-1 + (ITRAN-1)*MXVEC2 + IVEC)

        CALL CC_SETDOT(WORK(KXTRAN),WORK(KXDOTS),MXTRAN2,MXVEC2,
     &                 IX2(CD),IR2(AB),ITRAN,IVEC)
        XCON(4) = WORK(KXCONS-1 + (ITRAN-1)*MXVEC2 + IVEC)

        CALL CC_SETDOT(WORK(KXTRAN),WORK(KXDOTS),MXTRAN2,MXVEC2,
     &                 IX2(BD),IR2(AC),ITRAN,IVEC)
        XCON(5) = WORK(KXCONS-1 + (ITRAN-1)*MXVEC2 + IVEC)

        CALL CC_SETDOT(WORK(KXTRAN),WORK(KXDOTS),MXTRAN2,MXVEC2,
     &                 IX2(BC),IR2(AD),ITRAN,IVEC)
        XCON(6) = WORK(KXCONS-1 + (ITRAN-1)*MXVEC2 + IVEC)
      ELSE 
        XCON(1) = ZERO
        XCON(2) = ZERO
        XCON(3) = ZERO
        XCON(4) = ZERO
        XCON(5) = ZERO
        XCON(6) = ZERO
      END IF


*---------------------------------------------------------------------*
* get xksi vector dot products, 4 permutations
*---------------------------------------------------------------------*
      IF (L_USE_XKS3) THEN
        CALL CC_SETDOT(WORK(KOTRAN),WORK(KODOTS),MXTRAN3,MXVEC1,
     &                 IO3(ABC),IZT(D),ITRAN,IVEC)
        OCON(1) = WORK(KOCONS-1 + (ITRAN-1)*MXVEC1 + IVEC)

        CALL CC_SETDOT(WORK(KOTRAN),WORK(KODOTS),MXTRAN3,MXVEC1,
     &                 IO3(ABD),IZT(C),ITRAN,IVEC)
        OCON(2) = WORK(KOCONS-1 + (ITRAN-1)*MXVEC1 + IVEC)

        CALL CC_SETDOT(WORK(KOTRAN),WORK(KODOTS),MXTRAN3,MXVEC1,
     &                 IO3(ACD),IZT(B),ITRAN,IVEC)
        OCON(3) = WORK(KOCONS-1 + (ITRAN-1)*MXVEC1 + IVEC)

        CALL CC_SETDOT(WORK(KOTRAN),WORK(KODOTS),MXTRAN3,MXVEC1,
     &                 IO3(BCD),IZT(A),ITRAN,IVEC)
        OCON(4) = WORK(KOCONS-1 + (ITRAN-1)*MXVEC1 + IVEC)
      ELSE 
        OCON(1) = ZERO
        OCON(2) = ZERO
        OCON(3) = ZERO
        OCON(4) = ZERO
      END IF

*---------------------------------------------------------------------*
* add contributions to hyperpolarizabilities:
*---------------------------------------------------------------------*

      IDX = NCRFREQ*(IOPER-1)+ IFREQ

      IF (ISIGN.EQ.-1) IDX = IDX + NCRFREQ*NCROPER

      WORK(KHYPPOL-1+IDX) = WORK(KHYPPOL-1+IDX) + 
     &  H0CON    +
     &  G0CON(1) +G0CON(2) +G0CON(3) +G0CON(4)  +G0CON(5)  +G0CON(6)  + 
     &  GACON(1) +GACON(2) +GACON(3) +GACON(4)  +
     &  F0CON(1) +F0CON(2) +F0CON(3) +
     &  FACON(1) +FACON(2) +FACON(3) +FACON(4)  +FACON(5)  +FACON(6)  + 
     &  FACON(7) +FACON(8) +FACON(9) +FACON(10) +FACON(11) +FACON(12) + 
     &  F0ACON(1)+F0ACON(2)+F0ACON(3)+F0ACON(4) +F0ACON(5) +F0ACON(6) + 
     &  F0ACON(7)+F0ACON(8)+F0ACON(9)+F0ACON(10)+F0ACON(11)+F0ACON(12)+
     &  FAACON(1)+FAACON(2)+FAACON(3)+FAACON(4) +FAACON(5) +FAACON(6) + 
     &  FAACON(7)+FAACON(8)+FAACON(9)+FAACON(10)+FAACON(11)+FAACON(12)+
     &  EAACON(1)+EAACON(2)+EAACON(3)+EAACON(4) +EAACON(5) +EAACON(6) + 
     &  EAACON(7)+EAACON(8)+EAACON(9)+EAACON(10)+EAACON(11)+EAACON(12)+
     &  XCON(1)  +XCON(2)  +XCON(3)  +XCON(4)   +XCON(5)   +XCON(6)   +
     &  OCON(1)  +OCON(2)  +OCON(3)  +OCON(4)   +
     &  LCON(1)  +LCON(2)  +LCON(3)

      IF (LOCDBG) THEN
        WRITE(LUPRI,'(A,3I5)') 'IOPER, IFREQ, ISIGN:',IOPER,IFREQ,ISIGN
        WRITE(LUPRI,'(A,4I5)') 'IZT:',(IZT(I),I=1,4)
        WRITE(LUPRI,'(A,4I5)') 'IAM:',(IAM(I),I=1,4)
        WRITE(LUPRI,'(A,6I5)') 'IR2:',(IR2(I),I=1,6)
        IF (L_USE_CHI2) WRITE(LUPRI,'(A,6I5)') 'IX2:',(IX2(I),I=1,6)
        IF (L_USE_XKS3) WRITE(LUPRI,'(A,6I5)') 'IO3:',(IO3(I),I=1,4)
        WRITE(LUPRI,*) 'H0CON: ',H0CON
        WRITE(LUPRI,*) 'G0CON: ',(G0CON(I),I=1,6)
        WRITE(LUPRI,*) 'GACON: ',(GACON(I),I=1,4)
        WRITE(LUPRI,*) 'F0CON: ',(F0CON(I),I=1,3)
        WRITE(LUPRI,*) 'FACON: ',(FACON(I),I=1,12)
        WRITE(LUPRI,*) 'F0ACON:',(F0ACON(I),I=1,12)
        WRITE(LUPRI,*) 'FAACON:',(FAACON(I),I=1,12)
        WRITE(LUPRI,*) 'EAACON:',(EAACON(I),I=1,12)
        WRITE(LUPRI,*) 'XCON:',(XCON(I),I=1,6)
        WRITE(LUPRI,*) 'OCON:',(OCON(I),I=1,4)
        WRITE(LUPRI,*) 'LCON:',(LCON(I),I=1,3)
      END IF

      END DO 
      END DO 
      END IF
      END DO 

*---------------------------------------------------------------------*
* print frequency-dependent hyperpolarizabilities
*---------------------------------------------------------------------*
      WRITE (LUPRI,'(/A,I4,A,F12.2," seconds.")') '>>> Total time for',
     &  NBHYPPOL,' cubic response functions:   ', SECOND() - TIM0

      CALL CCCHYPPRT(WORK(KHYPPOL))

*---------------------------------------------------------------------*
* now we are finished with frequency-dependent hyperpolarizabilities...
*---------------------------------------------------------------------*
      END IF ! (NCRFREQ.GT.0) 

      IF (NCRDISP.EQ.0) THEN
        RETURN
      END IF

*---------------------------------------------------------------------*
* allocate & initialize work space for expansion coefficients
*---------------------------------------------------------------------*
      NBEXPCOF = NCROPER * NCRDISP

      KEXPCOF = 1
      KEND0   = KEXPCOF + NBEXPCOF
      KEND0A  = KEND0


      MXTRAN2 = 2 * NLRCLBL * NLRCLBL
      MXVEC2  = (NLRCLBL+1) * NLRCLBL /2
      MXTRAN3 = 2 * NLRCLBL * NLRCLBL * NLRCLBL
      MXVEC1  = NLRCLBL
      MXTRAN2 = MIN(MXTRAN2,2*12*NBEXPCOF) ! 12 = max. # of permut.
      MXTRAN3 = MIN(MXTRAN3,2*12*NBEXPCOF) ! 12 = max. # of permut.
      MXVEC2  = MIN(MXVEC2,2*12*NBEXPCOF)  ! 12 = max. # of permut.
      MXVEC1  = MIN(MXVEC1,2*12*NBEXPCOF)  ! 12 = max. # of permut.

      MX0HTRAN  = 5 * MXTRAN3
      MX0GTRAN  = 4 * MXTRAN2
      MXAGTRAN  = 4 * MXTRAN3
      MX0FTRAN  = 3 * MXTRAN2
      MXAFTRAN  = 3 * MXTRAN2
      MX0FATRAN = 5 * MXTRAN2
      MXAFATRAN = 5 * MXTRAN3
      MXAEATRAN = 3 * MXTRAN2
      MXXTRAN   = 1 * MXTRAN2
      MXOTRAN   = 1 * MXTRAN2
      MXLTRAN   = 1 * MXTRAN2

      MX0HDOTS  = MXTRAN3 * MXVEC1
      MX0GDOTS  = MXTRAN2 * MXVEC2
      MXAGDOTS  = MXTRAN3 * MXVEC1
      MX0FDOTS  = MXTRAN2 * MXVEC2
      MXAFDOTS  = MXTRAN2 * MXVEC2
      MX0FADOTS = MXTRAN2 * MXVEC2
      MXAFADOTS = MXTRAN3 * MXVEC1
      MXAEADOTS = MXTRAN2 * MXVEC2
      MXXDOTS   = MXTRAN2 * MXVEC2
      MXODOTS   = MXTRAN2 * MXVEC2
      MXLDOTS   = MXTRAN2 * MXVEC2

      K0HTRAN  = KEND0
      K0HDOTS  = K0HTRAN  + MX0HTRAN
      K0HCONS  = K0HDOTS  + MX0HDOTS
      KEND0    = K0HCONS  + MX0HDOTS

      K0GTRAN  = KEND0
      K0GDOTS  = K0GTRAN  + MX0GTRAN
      K0GCONS  = K0GDOTS  + MX0GDOTS
      KEND0    = K0GCONS  + MX0GDOTS

      KAGTRAN  = KEND0
      KAGDOTS  = KAGTRAN  + MXAGTRAN
      KAGCONS  = KAGDOTS  + MXAGDOTS
      KEND0    = KAGCONS  + MXAGDOTS

      K0FTRAN  = KEND0
      K0FDOTS  = K0FTRAN  + MX0FTRAN
      K0FCONS  = K0FDOTS  + MX0FDOTS
      KEND0    = K0FCONS  + MX0FDOTS

      KAFTRAN  = KEND0
      KAFDOTS  = KAFTRAN  + MXAFTRAN
      KAFCONS  = KAFDOTS  + MXAFDOTS
      KEND0    = KAFCONS  + MXAFDOTS

      K0FATRAN = KEND0
      K0FADOTS = K0FATRAN + MX0FATRAN
      K0FACONS = K0FADOTS + MX0FADOTS
      KEND0    = K0FACONS + MX0FADOTS

      KAFATRAN = KEND0
      KAFADOTS = KAFATRAN + MXAFATRAN
      KAFACONS = KAFADOTS + MXAFADOTS
      KEND0    = KAFACONS + MXAFADOTS

      KAEATRAN = KEND0
      KAEADOTS = KAEATRAN + MXAEATRAN
      KAEACONS = KAEADOTS + MXAEADOTS
      KEND0    = KAEACONS + MXAEADOTS

      KXTRAN   = KEND0
      KXDOTS   = KXTRAN   + MXXTRAN
      KXCONS   = KXDOTS   + MXXDOTS
      KEND0    = KXCONS   + MXXDOTS

      KOTRAN   = KEND0
      KODOTS   = KOTRAN   + MXOTRAN
      KOCONS   = KODOTS   + MXODOTS
      KEND0    = KOCONS   + MXODOTS

      KLTRAN   = KEND0
      KLDOTS   = KLTRAN   + MXLTRAN
      KLCONS   = KLDOTS   + MXLDOTS
      KEND0    = KLCONS   + MXLDOTS

      LEND0    = LWORK - KEND0

      IF (LEND0 .LT. 0) THEN
       WRITE (LUPRI,*) 'Insufficient memory in CC_2HYP. (2)'
       WRITE (LUPRI,*) 'available:',LWORK
       WRITE (LUPRI,*) '   needed:',KEND0
       WRITE (LUPRI,*) '  MXTRAN2:',MXTRAN2
       WRITE (LUPRI,*) '  MXTRAN3:',MXTRAN3
       WRITE (LUPRI,*) '  MXVEC1 :',MXVEC1
       WRITE (LUPRI,*) '  MXVEC2 :',MXVEC2
       CALL QUIT('Insufficient memory in CC_2HYP. (2)')
      END IF


* initialize expansion coefficients and the lists of contributions:
      CALL DZERO(WORK(KEXPCOF),NBEXPCOF)
      CALL DZERO(WORK(K0HTRAN),MX0HTRAN+2*MX0HDOTS)
      CALL DZERO(WORK(K0GTRAN),MX0GTRAN+2*MX0GDOTS)
      CALL DZERO(WORK(KAGTRAN),MXAGTRAN+2*MXAGDOTS)
      CALL DZERO(WORK(K0FTRAN),MX0FTRAN+2*MX0FDOTS)
      CALL DZERO(WORK(KAFTRAN),MXAFTRAN+2*MXAFDOTS)
      CALL DZERO(WORK(K0FATRAN),MX0FATRAN+2*MX0FADOTS)
      CALL DZERO(WORK(KAFATRAN),MXAFATRAN+2*MXAFADOTS)
      CALL DZERO(WORK(KAEATRAN),MXAEATRAN+2*MXAEADOTS)
      CALL DZERO(WORK(KXTRAN),MXXTRAN+2*MXXDOTS)
      CALL DZERO(WORK(KOTRAN),MXOTRAN+2*MXODOTS)
      CALL DZERO(WORK(KLTRAN),MXLTRAN+2*MXLDOTS)

* initialize timing:
      TIM2 = SECOND()

* save options L_USE_CHI2 and L_USE_XKS3, which are not implemented
* for the dispersion coefficients:
      LCHI2_SAVE = L_USE_CHI2 
      LXKS3_SAVE = L_USE_XKS3
      L_USE_CHI2 = .FALSE.
      L_USE_XKS3 = .FALSE.
     
*---------------------------------------------------------------------*
* set up lists for H, G, F and F{O} transformations and ETA{O} vectors:
*---------------------------------------------------------------------*
      LADD = .FALSE.

      CALL CCCR_DISP_SETUP(MXTRAN2, MXVEC2, MXTRAN3, MXVEC1,
     &          WORK(K0HTRAN), WORK(K0HDOTS), WORK(K0HCONS), N0HTRAN,
     &          WORK(K0GTRAN), WORK(K0GDOTS), WORK(K0GCONS), N0GTRAN,
     &          WORK(KAGTRAN), WORK(KAGDOTS), WORK(KAGCONS), NAGTRAN,
     &          WORK(K0FTRAN), WORK(K0FDOTS), WORK(K0FCONS), N0FTRAN,
     &          WORK(KAFTRAN), WORK(KAFDOTS), WORK(KAFCONS), NAFTRAN,
     &          WORK(K0FATRAN),WORK(K0FADOTS),WORK(K0FACONS),N0FATRAN,
     &          WORK(KAFATRAN),WORK(KAFADOTS),WORK(KAFACONS),NAFATRAN,
     &          WORK(KAEATRAN),WORK(KAEADOTS),WORK(KAEACONS),NAEATRAN,
     &          WORK(KXTRAN),  WORK(KXDOTS),  WORK(KXCONS),  NXTRAN,    
     &          WORK(KOTRAN),  WORK(KODOTS),  WORK(KOCONS),  NOTRAN,
     &          WORK(KLTRAN),  WORK(KLDOTS),  WORK(KLCONS),  NLTRAN,
     &          WORK(KEXPCOF), NBEXPCOF, LADD    )

*---------------------------------------------------------------------*
* calculate H matrix contributions:
*---------------------------------------------------------------------*
      TIM1 = SECOND()

      IOPT = 5
      CALL CC_HMAT('L0 ','RC ','RC ','RC ','RC ',N0HTRAN, MXVEC1,
     &             WORK(K0HTRAN),WORK(K0HDOTS),WORK(K0HCONS),
     &             WORK(KEND0), LEND0, IOPT )

      TIMH = SECOND() - TIM1
      WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")')
     & '>>> Time used for',N0HTRAN,' H matrix transformations:   ',TIMH
      CALL FLSHFO(LUPRI)

*---------------------------------------------------------------------*
* calculate G matrix contributions:
*---------------------------------------------------------------------*
      TIM1 = SECOND()

      IF ((.NOT.L_USE_CHI2)  .AND. NO_2NP1_RULE) THEN
       IOPT = 5
       CALL CC_GMATRIX('L0 ','RC ','RC ','CR2',N0GTRAN, MXVEC2,
     &               WORK(K0GTRAN),WORK(K0GDOTS),WORK(K0GCONS),
     &               WORK(KEND0), LEND0, IOPT )
      END IF

      IF (.NOT.L_USE_XKS3) THEN
        IOPT = 5
        CALL CC_GMATRIX('LC ','RC ','RC ','RC ',NAGTRAN, MXVEC1,
     &                WORK(KAGTRAN),WORK(KAGDOTS),WORK(KAGCONS),
     &                WORK(KEND0), LEND0, IOPT )
      END IF
     
      TIMG = SECOND() - TIM1
      WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")')
     &  '>>> Time used for',N0GTRAN+NAGTRAN,
     &           ' G matrix transformations:   ',TIMG
      CALL FLSHFO(LUPRI)
      
*---------------------------------------------------------------------*
* calculate F matrix contributions:
*---------------------------------------------------------------------*
      TIM1 = SECOND()

      IOPT = 5
      CALL CC_FMATRIX(WORK(K0FTRAN),N0FTRAN,'L0 ','CR2',IOPT,'CR2',
     &                WORK(K0FDOTS),WORK(K0FCONS),MXVEC2,
     &                WORK(KEND0),LEND0)

      IF (L_USE_CHI2) THEN
        CALL DSCAL(MX0FDOTS,-1.0d0,WORK(K0FCONS),1)
      END IF

      IF ((.NOT. (L_USE_CHI2 .OR. L_USE_XKS3)) .AND. NO_2NP1_RULE) THEN
         IOPT = 5
         CALL CC_FMATRIX(WORK(KAFTRAN),NAFTRAN,'LC ','RC ',IOPT,'CR2',
     &                   WORK(KAFDOTS),WORK(KAFCONS),MXVEC2,
     &                   WORK(KEND0),LEND0)
      END IF

      TIMF = SECOND() - TIM1
      WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")')
     &  '>>> Time used for',N0FTRAN+NAFTRAN,
     &          ' F matrix transformations:   ',TIMF
      CALL FLSHFO(LUPRI)
      
*---------------------------------------------------------------------*
* calculate F{O} matrix contributions:
*---------------------------------------------------------------------*
      TIM1 = SECOND()

      IF ((.NOT.L_USE_CHI2) .AND. NO_2NP1_RULE) THEN
       CALL CCQR_FADRV('L0 ','o1 ','RC ','CR2',N0FATRAN, MXVEC2,
     &                  WORK(K0FATRAN),WORK(K0FADOTS),
     &                  WORK(K0FACONS),
     &                  WORK(KEND0), LEND0, 'DOTP' )
      END IF

      IF (.NOT.L_USE_XKS3) THEN
        CALL CCQR_FADRV('LC ','o1 ','RC ','RC ',NAFATRAN, MXVEC1,
     &                   WORK(KAFATRAN),WORK(KAFADOTS),
     &                   WORK(KAFACONS),
     &                   WORK(KEND0), LEND0, 'DOTP' )
      END IF

      TIMFA = SECOND() - TIM1
      WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")')
     &  '>>> Time used for',N0FATRAN+NAFATRAN,
     &          ' F{O} matrix transformations:',TIMFA
      CALL FLSHFO(LUPRI)
      
*---------------------------------------------------------------------*
* calculate ETA{O} vector contributions:
*---------------------------------------------------------------------*
      TIM1 = SECOND()

      IF ((.NOT.(L_USE_CHI2 .OR. L_USE_XKS3)) .AND. NO_2NP1_RULE ) THEN
       CALL CCQR_EADRV('LC ','o1 ','CR2',NAEATRAN, MXVEC2,
     &                  WORK(KAEATRAN),WORK(KAEADOTS),WORK(KAEACONS),
     &                  WORK(KEND0), LEND0, 'DOTP' )
      END IF

      TIMEA = SECOND() - TIM1
      WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")') '>>> Time used for',
     &  NAEATRAN,' ETA{O} vector calculations :',TIMEA
      CALL FLSHFO(LUPRI)
      
*---------------------------------------------------------------------*
* chi vector contributions:
*---------------------------------------------------------------------*
      IF (.NOT. NO_2NP1_RULE) THEN
        TIM1 = SECOND()

        CALL CC_DOTDRV('CX2','CR2',NXTRAN,MXVEC2,
     &                 WORK(KXTRAN), WORK(KXDOTS), WORK(KXCONS),
     &                 WORK(KEND0), LEND0 )

        TIMX2 = SECOND() - TIM1
        WRITE (LUPRI,'(/A,I5,A,5X,F12.2," seconds.")')
     &    '>>> Time used for',NXTRAN, ' CX2 x CR2 dot products:',TIMX2
        CALL FLSHFO(LUPRI)
      
      END IF
*---------------------------------------------------------------------*
* CL2 x CR2 dot products:
*---------------------------------------------------------------------*
      IF (.NOT. NO_2NP1_RULE) THEN
        TIM1 = SECOND()

        CALL CC_DOTDRV('CL2','CR2',NLTRAN,MXVEC2,
     &                 WORK(KLTRAN), WORK(KLDOTS), WORK(KLCONS),
     &                 WORK(KEND0), LEND0 )

        TIML2 = SECOND() - TIM1
        WRITE (LUPRI,'(/A,I5,A,5X,F12.2," seconds.")')
     &    '>>> Time used for',NLTRAN, ' CL2 x CR2 dot products:',TIML2
        CALL FLSHFO(LUPRI)
      
      END IF
*---------------------------------------------------------------------*
* xksi vector contributions:
*---------------------------------------------------------------------*
      IF (.NOT. NO_2NP1_RULE) THEN

        TIM1 = SECOND()

        CALL CC_DOTDRV('CO2','CL2',NOTRAN,MXVEC2,
     &                 WORK(KOTRAN), WORK(KODOTS), WORK(KOCONS),
     &                 WORK(KEND0), LEND0 )

        TIMO2 = SECOND() - TIM1
        WRITE (LUPRI,'(/A,I5,A,5X,F12.2," seconds.")')
     &    '>>> Time used for',NOTRAN, ' CL2 x CO2 dot products:',TIMO2
        CALL FLSHFO(LUPRI)

      END IF

*---------------------------------------------------------------------*
* collect contributions and add them to hyperpolarizabilities
*---------------------------------------------------------------------*
      LADD = .TRUE.

      CALL CCCR_DISP_SETUP(MXTRAN2, MXVEC2, MXTRAN3, MXVEC1,
     &          WORK(K0HTRAN), WORK(K0HDOTS), WORK(K0HCONS), N0HTRAN,
     &          WORK(K0GTRAN), WORK(K0GDOTS), WORK(K0GCONS), N0GTRAN,
     &          WORK(KAGTRAN), WORK(KAGDOTS), WORK(KAGCONS), NAGTRAN,
     &          WORK(K0FTRAN), WORK(K0FDOTS), WORK(K0FCONS), N0FTRAN,
     &          WORK(KAFTRAN), WORK(KAFDOTS), WORK(KAFCONS), NAFTRAN,
     &          WORK(K0FATRAN),WORK(K0FADOTS),WORK(K0FACONS),N0FATRAN,
     &          WORK(KAFATRAN),WORK(KAFADOTS),WORK(KAFACONS),NAFATRAN,
     &          WORK(KAEATRAN),WORK(KAEADOTS),WORK(KAEACONS),NAEATRAN,
     &          WORK(KXTRAN),  WORK(KXDOTS),  WORK(KXCONS),  NXTRAN,    
     &          WORK(KOTRAN),  WORK(KODOTS),  WORK(KOCONS),  NOTRAN,
     &          WORK(KLTRAN),  WORK(KLDOTS),  WORK(KLCONS),  NLTRAN,
     &          WORK(KEXPCOF), NBEXPCOF, LADD    )


* recover options L_USE_CHI2 and L_USE_XKS3, which are not implemented
* for the dispersion coefficients:
      L_USE_CHI2 = LCHI2_SAVE
      L_USE_XKS3 = LXKS3_SAVE

*---------------------------------------------------------------------*
* print timing:
*---------------------------------------------------------------------*
      WRITE (LUPRI,'(/A,I5,A,F12.2," seconds.")') '>>> Total time for',
     &  NBEXPCOF,' expansion coefficients:    ', SECOND() - TIM2

*---------------------------------------------------------------------*
* print expansion coefficients for second hyperpolarizabilities:
*---------------------------------------------------------------------*
      IF (IPRCHYP .GE. 15) THEN
        CALL CCCEXPPRT(WORK(KEXPCOF))
      END IF

      CALL FLSHFO(LUPRI)

*---------------------------------------------------------------------*
* calculate from the d_ABCD expansion coeff. the D_ABCD disp coeff.:
*---------------------------------------------------------------------*
      IF (NCRDSPE.GE.0) THEN
        NTRINOM = (NCRDSPE+3)*(NCRDSPE+2)*(NCRDSPE+1)/6

        KTRINOM = KEND0A
        KDSPCF  = KTRINOM + NTRINOM
        KTHGCF  = KDSPCF  + NTRINOM * NCROPER
        KESHGCF = KTHGCF  + (NCRDSPE+1) * NCROPER
        KKERRCF = KESHGCF + (NCRDSPE+1) * NCROPER
        KDFWMCF = KKERRCF + (NCRDSPE+1) * NCROPER
        KEND1   = KDFWMCF + (NCRDSPE+1) * NCROPER
        LEND1   = LWORK - KEND1

        IF (GAMMA_PAR .OR. GAMMA_ORT) THEN
          KDSPCFA = KEND1
          KTHGCFA = KDSPCFA + NTRINOM * 4
          KESHGCA = KTHGCFA + (NCRDSPE+1) * 4
          KKERRCA = KESHGCA + (NCRDSPE+1) * 4
          KDFWMCA = KKERRCA + (NCRDSPE+1) * 4
          KABCDE  = KDFWMCA + (NCRDSPE+1) * 4
          KEND1   = KABCDE  + 15 * 2
          LEND1   = LWORK - KEND1
        END IF

        IF (LEND1 .LT. 0) THEN
          CALL QUIT('Insufficient memory in CC_2HYP. (3)')
        END IF

        CALL DZERO(WORK(KTRINOM),KEND1-KTRINOM)

        CALL CCCRDISP(WORK(KDSPCF),WORK(KDSPCFA),WORK(KABCDE),
     &                WORK(KTHGCF),WORK(KESHGCF),
     &                WORK(KKERRCF),WORK(KDFWMCF),
     &                WORK(KTHGCFA),WORK(KESHGCA),
     &                WORK(KKERRCA),WORK(KDFWMCA),
     &                WORK(KEXPCOF),WORK(KTRINOM), NTRINOM,LERROR)

        CALL CCCDISPRT(WORK(KDSPCF),WORK(KDSPCFA),WORK(KABCDE),
     &                 WORK(KTHGCF),WORK(KESHGCF),
     &                 WORK(KKERRCF),WORK(KDFWMCF),
     &                 WORK(KTHGCFA),WORK(KESHGCA),
     &                 WORK(KKERRCA),WORK(KDFWMCA),LERROR)
      
      END IF

      RETURN
      END

*=====================================================================*
*              END OF SUBROUTINE CC_2HYP                              *
*=====================================================================*
       SUBROUTINE CCCHYPPRT(HYPERPOL)
*---------------------------------------------------------------------*
*
*    Purpose: print second hyperpolarizabilities 
*
*
*     Written by Christof Haettig in February 1997.
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
#include <cccrinf.h>
#include <ccroper.h>


      CHARACTER*10 BLANKS
      CHARACTER*80 STRING
      INTEGER ISYMA, ISYMB, ISYMC, ISYMD
      INTEGER IFREQ, IOPER, ISYOLD, ICOMP, IADD


#if defined (SYS_CRAY)
      REAL HYPERPOL(NCRFREQ,NCROPER,2)
      REAL HALF, GAMMA, PROP
#else
      DOUBLE PRECISION HYPERPOL(NCRFREQ,NCROPER,2)
      DOUBLE PRECISION HALF, GAMMA, PROP
#endif
      PARAMETER (HALF = 0.5D0)

      CHARACTER MODEL*10,MOPRPC
      INTEGER ISYMABCD
*---------------------------------------------------------------------*
* print header for hyperpolarizability section
*---------------------------------------------------------------------*
      BLANKS = '          '
      STRING = ' RESULTS FOR THE SECOND HYPERPOLARIZABILITIES '

      IF (CCS) THEN
         CALL AROUND( BLANKS//'FINAL CCS'//STRING(1:45)//BLANKS ) 
      ELSE IF (CC2) THEN
         CALL AROUND( BLANKS//'FINAL CC2'//STRING(1:45)//BLANKS )
      ELSE IF (CCSD) THEN
         CALL AROUND( BLANKS//'FINAL CCSD'//STRING(1:45)//BLANKS )
      ELSE IF (CC3) THEN
         CALL AROUND( BLANKS//'FINAL CC3'//STRING(1:45)//BLANKS )
      ELSE
         CALL QUIT('CCCHYPPRT called for an unknown Coupled '//
     &             'Cluster model.')
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

      IF (IPRCHYP.GT.5) THEN
       WRITE(LUPRI,'(/1X,4(1X,A," operator",7X),4X,A,8X,A,/,107("-"))')
     &     'A','B','C','D','gamma','(asy. Resp.)'
      ELSE
       WRITE(LUPRI,'(/1X,4(1X,A," operator",7X),4X,A,/,86("-"))')
     &     'A','B','C','D','gamma'
      END IF

      ISYOLD = 1

      DO IOPER = 1, NCROPER
         ISYMA = ISYOPR(IACROP(IOPER))
         ISYMB = ISYOPR(IBCROP(IOPER))
         ISYMC = ISYOPR(ICCROP(IOPER))
         ISYMD = ISYOPR(IDCROP(IOPER))
         ISYMABCD = MULD2H(MULD2H(ISYMA,ISYMB),MULD2H(ISYMC,ISYMD))


         IFREQ = 1
cmbh initialize
         PROP=0.0d0
cmbh end
         IF (MULD2H(ISYMA,ISYMB).EQ.MULD2H(ISYMC,ISYMD)) THEN
          PROP = HALF*(HYPERPOL(IFREQ,IOPER,1)+HYPERPOL(IFREQ,IOPER,2))
          IF (IPRCHYP.GT.5) THEN
            WRITE(LUPRI,'(/1X,4(A8,F7.4,3X),G18.10," (",G18.10,")")')
     &        LBLOPR(IACROP(IOPER)),ACRFR(IFREQ),
     &        LBLOPR(IBCROP(IOPER)),BCRFR(IFREQ),
     &        LBLOPR(ICCROP(IOPER)),CCRFR(IFREQ),
     &        LBLOPR(IDCROP(IOPER)),DCRFR(IFREQ),
     &        HALF*(HYPERPOL(IFREQ,IOPER,1)+HYPERPOL(IFREQ,IOPER,2)),
     &        HALF*(HYPERPOL(IFREQ,IOPER,1)-HYPERPOL(IFREQ,IOPER,2))
          ELSE
            WRITE(LUPRI,'(/1X,4(A8,F7.4,3X),G16.8)')
     &        LBLOPR(IACROP(IOPER)),ACRFR(IFREQ),
     &        LBLOPR(IBCROP(IOPER)),BCRFR(IFREQ),
     &        LBLOPR(ICCROP(IOPER)),CCRFR(IFREQ),
     &        LBLOPR(IDCROP(IOPER)),DCRFR(IFREQ),
     &        HALF*(HYPERPOL(IFREQ,IOPER,1)+HYPERPOL(IFREQ,IOPER,2))
          END IF
          ISYOLD = 1
         ELSE
          IF (IPRCHYP.GT.5) THEN
            IF (ISYOLD.EQ.1) WRITE(LUPRI,*)
            WRITE(LUPRI,'(1X,4(A8,A7,3X),7X,A,9X," (",5X,A,6X,")")')
     &        LBLOPR(IACROP(IOPER)),'   -.-  ',
     &        LBLOPR(IBCROP(IOPER)),'   -.-  ',
     &        LBLOPR(ICCROP(IOPER)),'   -.-  ', 
     &        LBLOPR(IDCROP(IOPER)),'   -.-  ', 
     &        '---',
     &        '---'
          ELSE IF (IPRCHYP.GT.0) THEN
            IF (ISYOLD.EQ.1) WRITE(LUPRI,*)
            WRITE(LUPRI,'(1X,4(A8,A7,3X),6X,A,7X)')
     &        LBLOPR(IACROP(IOPER)),'   -.-  ',
     &        LBLOPR(IBCROP(IOPER)),'   -.-  ',
     &        LBLOPR(ICCROP(IOPER)),'   -.-  ', 
     &        LBLOPR(IDCROP(IOPER)),'   -.-  ', 
     &        '---'
          END IF
          ISYOLD = 0
         END IF
         CALL WRIPRO(PROP,model,4,
     &               LBLOPR(IACROP(IOPER)),LBLOPR(IBCROP(IOPER)),
     &               LBLOPR(ICCROP(IOPER)),LBLOPR(IDCROP(IOPER)),
     &               BCRFR(IFREQ),CCRFR(IFREQ),DCRFR(IFREQ),ISYMABCD,
     &               0,0,0)


         DO IFREQ = 2, NCRFREQ
          PROP = 0.0D0
          IF (MULD2H(ISYMA,ISYMB).EQ.MULD2H(ISYMC,ISYMD)) THEN
           PROP = HALF*(HYPERPOL(IFREQ,IOPER,1)+HYPERPOL(IFREQ,IOPER,2))
           IF (IPRCHYP.GT.5) THEN
            WRITE(LUPRI,'(1X,4(8X,F7.4,3X),G18.10," (",G18.10,")")')
     &        ACRFR(IFREQ), BCRFR(IFREQ), CCRFR(IFREQ), DCRFR(IFREQ),
     &        HALF*(HYPERPOL(IFREQ,IOPER,1)+HYPERPOL(IFREQ,IOPER,2)),
     &        HALF*(HYPERPOL(IFREQ,IOPER,1)-HYPERPOL(IFREQ,IOPER,2))
           ELSE
            WRITE(LUPRI,'(1X,4(8X,F7.4,3X),G16.8)')
     &        ACRFR(IFREQ), BCRFR(IFREQ), CCRFR(IFREQ), DCRFR(IFREQ),
     &        HALF*(HYPERPOL(IFREQ,IOPER,1)+HYPERPOL(IFREQ,IOPER,2))
           END IF
          END IF
cmbh
          CALL WRIPRO(PROP,model,4,
     &               LBLOPR(IACROP(IOPER)),LBLOPR(IBCROP(IOPER)),
     &               LBLOPR(ICCROP(IOPER)),LBLOPR(IDCROP(IOPER)),
     &               BCRFR(IFREQ),CCRFR(IFREQ),DCRFR(IFREQ),ISYMABCD,
     &               0,0,0)
cmbh end
         END DO

      END DO

      WRITE(LUPRI,'(86("-"))')


      IF (GAMMA_PAR .OR. GAMMA_ORT) THEN
       WRITE(LUPRI,'(/////A,40X,A,/,75("-"))')
     &   ' average      frequencies','value'

       IF (GAMMA_PAR) THEN

* calculate & print gamma_{||} for all requested frequencies:
       DO IFREQ = 1, NCRFREQ

         IF (CSYM(1:6).EQ.'ATOMIC') THEN
           GAMMA = HALF*(HYPERPOL(IFREQ,1,1)+HYPERPOL(IFREQ,1,2))
         ELSE IF (CSYM(1:6).EQ.'SPHTOP') THEN
           GAMMA = 0.3d0 * (HYPERPOL(IFREQ,1,1)+HYPERPOL(IFREQ,1,2))
           DO ICOMP = 2, 4
             GAMMA = GAMMA + 0.2d0 * 
     &          (HYPERPOL(IFREQ,ICOMP,1)+HYPERPOL(IFREQ,ICOMP,2))
           END DO
         ELSE IF (CSYM(1:6).EQ.'LINEAR') THEN
           GAMMA=      1.5d0*(HYPERPOL(IFREQ,1,1)+HYPERPOL(IFREQ,1,2))
           GAMMA=GAMMA+4.0d0*(HYPERPOL(IFREQ,8,1)+HYPERPOL(IFREQ,8,2))
           DO ICOMP = 2, 7
             GAMMA = GAMMA + 
     &          (HYPERPOL(IFREQ,ICOMP,1)+HYPERPOL(IFREQ,ICOMP,2))
           END DO
           GAMMA = GAMMA / 15.0d0
         ELSE 
           GAMMA =         HYPERPOL(IFREQ,1,1)  + HYPERPOL(IFREQ,1,2)
           GAMMA = GAMMA + HYPERPOL(IFREQ,8,1)  + HYPERPOL(IFREQ,8,2)
           GAMMA = GAMMA + HYPERPOL(IFREQ,15,1) + HYPERPOL(IFREQ,15,2)
           DO ICOMP = 1, 21
             GAMMA = GAMMA + 
     &        HALF*(HYPERPOL(IFREQ,ICOMP,1)+HYPERPOL(IFREQ,ICOMP,2))
           END DO
           GAMMA = GAMMA / 15.0d0
         END IF

         IF (IFREQ.EQ.1) THEN
           WRITE(LUPRI,'(1X,A8,4X,F7.4,2X,3(3X,F7.4,2X),G16.8)')
     &        'gamma_||', ACRFR(IFREQ), BCRFR(IFREQ), 
     &                    CCRFR(IFREQ), DCRFR(IFREQ), GAMMA
         ELSE
           WRITE(LUPRI,'(1X,8X,4X,F7.4,2X,3(3X,F7.4,2X),G16.8)')
     &                    ACRFR(IFREQ), BCRFR(IFREQ), 
     &                    CCRFR(IFREQ), DCRFR(IFREQ), GAMMA
         END IF
       END DO 

       END IF

       IF (GAMMA_ORT) THEN

* calculate & print gamma_{_|_} for all requested frequencies:
       DO IFREQ = 1, NCRFREQ

         IF (CSYM(1:6).EQ.'ATOMIC') THEN
           GAMMA =      (HYPERPOL(IFREQ,2,1)+HYPERPOL(IFREQ,2,2))
     &           - HALF*(HYPERPOL(IFREQ,1,1)+HYPERPOL(IFREQ,1,2))
     &           +      (HYPERPOL(IFREQ,3,1)+HYPERPOL(IFREQ,3,2))
         ELSE IF (CSYM(1:6).EQ.'SPHTOP') THEN
           GAMMA = 0.1d0 * (HYPERPOL(IFREQ,1,1)+HYPERPOL(IFREQ,1,2))
     &           + 0.4d0 * (HYPERPOL(IFREQ,2,1)+HYPERPOL(IFREQ,2,2))
     &           - 0.1d0 * (HYPERPOL(IFREQ,3,1)+HYPERPOL(IFREQ,3,2))
     &           - 0.1d0 * (HYPERPOL(IFREQ,4,1)+HYPERPOL(IFREQ,4,2))
         ELSE IF (CSYM(1:6).EQ.'LINEAR') THEN
           GAMMA=  1.0d0 * (HYPERPOL(IFREQ,1,1)+HYPERPOL(IFREQ,1,2))
     &           + 4.0d0 * (HYPERPOL(IFREQ,2,1)+HYPERPOL(IFREQ,2,2))
     &           - 1.0d0 * (HYPERPOL(IFREQ,3,1)+HYPERPOL(IFREQ,3,2))
     &           - 1.0d0 * (HYPERPOL(IFREQ,4,1)+HYPERPOL(IFREQ,4,2))
     &           + 4.0d0 * (HYPERPOL(IFREQ,5,1)+HYPERPOL(IFREQ,5,2))
     &           - 1.0d0 * (HYPERPOL(IFREQ,6,1)+HYPERPOL(IFREQ,6,2))
     &           - 1.0d0 * (HYPERPOL(IFREQ,7,1)+HYPERPOL(IFREQ,7,2))
     &           + 1.0d0 * (HYPERPOL(IFREQ,8,1)+HYPERPOL(IFREQ,8,2))
     &           + 5.0d0 * (HYPERPOL(IFREQ,9,1)+HYPERPOL(IFREQ,9,2))
           GAMMA = GAMMA / 30.0d0
         ELSE 
           GAMMA = 0.0d0
           DO IADD = 0, 14, 7
             GAMMA = GAMMA + 
     &        1.0d0*(HYPERPOL(IFREQ,1+IADD,1)+HYPERPOL(IFREQ,1+IADD,2))
     &       +2.0d0*(HYPERPOL(IFREQ,2+IADD,1)+HYPERPOL(IFREQ,2+IADD,2))
     &       -0.5d0*(HYPERPOL(IFREQ,3+IADD,1)+HYPERPOL(IFREQ,3+IADD,2))
     &       -0.5d0*(HYPERPOL(IFREQ,4+IADD,1)+HYPERPOL(IFREQ,4+IADD,2))
     &       +2.0d0*(HYPERPOL(IFREQ,5+IADD,1)+HYPERPOL(IFREQ,5+IADD,2))
     &       -0.5d0*(HYPERPOL(IFREQ,6+IADD,1)+HYPERPOL(IFREQ,6+IADD,2))
     &       -0.5d0*(HYPERPOL(IFREQ,7+IADD,1)+HYPERPOL(IFREQ,7+IADD,2))
           END DO
           GAMMA = GAMMA / 30.0d0
         END IF

         IF (IFREQ.EQ.1) THEN
           WRITE(LUPRI,'(/1X,A8,4X,F7.4,2X,3(3X,F7.4,2X),G16.8)')
     &        'gamma_|_', ACRFR(IFREQ), BCRFR(IFREQ), 
     &                    CCRFR(IFREQ), DCRFR(IFREQ), GAMMA
         ELSE
           WRITE(LUPRI,'(1X,8X,4X,F7.4,2X,3(3X,F7.4,2X),G16.8)')
     &                    ACRFR(IFREQ), BCRFR(IFREQ), 
     &                    CCRFR(IFREQ), DCRFR(IFREQ), GAMMA
         END IF
       END DO 

       END IF

       IF (GAMMA_ORT) THEN

* calculate & print gamma_{ms} for all requested frequencies:
       DO IFREQ = 1, NCRFREQ

         IF (CSYM(1:6).EQ.'ATOMIC') THEN
           GAMMA = 3.0d0*HALF*(HYPERPOL(IFREQ,2,1)+HYPERPOL(IFREQ,2,2))
     &           - 2.0d0*HALF*(HYPERPOL(IFREQ,1,1)+HYPERPOL(IFREQ,1,2))
     &           + 3.0d0*HALF*(HYPERPOL(IFREQ,3,1)+HYPERPOL(IFREQ,3,2))
         ELSE IF (CSYM(1:6).EQ.'SPHTOP') THEN
           GAMMA =  0.5d0 * (HYPERPOL(IFREQ,2,1)+HYPERPOL(IFREQ,2,2))
     &            + 0.5d0 * (HYPERPOL(IFREQ,3,1)+HYPERPOL(IFREQ,3,2))
     &            - 1.0d0 * (HYPERPOL(IFREQ,4,1)+HYPERPOL(IFREQ,4,2))
         ELSE IF (CSYM(1:6).EQ.'LINEAR') THEN
           GAMMA =  0.5d0 * (HYPERPOL(IFREQ,2,1)+HYPERPOL(IFREQ,2,2))
     &            + 0.5d0 * (HYPERPOL(IFREQ,3,1)+HYPERPOL(IFREQ,3,2))
     &            - 1.0d0 * (HYPERPOL(IFREQ,4,1)+HYPERPOL(IFREQ,4,2))
     &            + 0.5d0 * (HYPERPOL(IFREQ,5,1)+HYPERPOL(IFREQ,5,2))
     &            + 0.5d0 * (HYPERPOL(IFREQ,6,1)+HYPERPOL(IFREQ,6,2))
     &            - 1.0d0 * (HYPERPOL(IFREQ,7,1)+HYPERPOL(IFREQ,7,2))
     &            + 1.5d0 * (HYPERPOL(IFREQ,9,1)+HYPERPOL(IFREQ,9,2))
     &            + 1.5d0 * (HYPERPOL(IFREQ,10,1)+HYPERPOL(IFREQ,10,2))
     &            - 1.0d0 * (HYPERPOL(IFREQ,8,1)+HYPERPOL(IFREQ,8,2))
           GAMMA = GAMMA / 3.0d0
         ELSE 
           GAMMA = 0.0d0
           DO IADD = 0, 14, 7
             GAMMA = GAMMA + 
     &        0.5d0*(HYPERPOL(IFREQ,2+IADD,1)+HYPERPOL(IFREQ,2+IADD,2))
     &       +0.5d0*(HYPERPOL(IFREQ,3+IADD,1)+HYPERPOL(IFREQ,3+IADD,2))
     &       -1.0d0*(HYPERPOL(IFREQ,4+IADD,1)+HYPERPOL(IFREQ,4+IADD,2))
     &       +0.5d0*(HYPERPOL(IFREQ,5+IADD,1)+HYPERPOL(IFREQ,5+IADD,2))
     &       +0.5d0*(HYPERPOL(IFREQ,6+IADD,1)+HYPERPOL(IFREQ,6+IADD,2))
     &       -1.0d0*(HYPERPOL(IFREQ,7+IADD,1)+HYPERPOL(IFREQ,7+IADD,2))
           END DO
           GAMMA = GAMMA / 6.0d0
         END IF

         IF (IFREQ.EQ.1) THEN
           WRITE(LUPRI,'(/1X,A8,4X,F7.4,2X,3(3X,F7.4,2X),F16.7)')
     &        'gamma_ms', ACRFR(IFREQ), BCRFR(IFREQ), 
     &                    CCRFR(IFREQ), DCRFR(IFREQ), GAMMA
         ELSE
           WRITE(LUPRI,'(1X,8X,4X,F7.4,2X,3(3X,F7.4,2X),G16.8)')
     &                    ACRFR(IFREQ), BCRFR(IFREQ), 
     &                    CCRFR(IFREQ), DCRFR(IFREQ), GAMMA
         END IF
       END DO 

       END IF

      END IF

      WRITE(LUPRI,'(75("-"))')

      CALL FLSHFO(LUPRI)

      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCCHYPPRT                            *
*---------------------------------------------------------------------*
       SUBROUTINE CCCEXPPRT(EXPCOF)
*---------------------------------------------------------------------*
*
*    Purpose: print expansion coefficients for
*             second hyperpolarizabilities 
*
*
*     Written by Christof Haettig in March 1998.
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
#include <cccrinf.h>
#include <ccroper.h>


      CHARACTER*10 BLANKS
      CHARACTER*80 STRING
      INTEGER ISYMA, ISYMB, ISYMC, ISYMD
      INTEGER IDISP, IOPER, ISYOLD, ICOMP, IADD


#if defined (SYS_CRAY)
      REAL EXPCOF(NCRDISP,NCROPER)
      REAL GAMMA
#else
      DOUBLE PRECISION EXPCOF(NCRDISP,NCROPER)
      DOUBLE PRECISION GAMMA
#endif

*---------------------------------------------------------------------*
* print header for hyperpolarizability section
*---------------------------------------------------------------------*
      BLANKS = '          '
      STRING = ' RESULTS FOR EXPANSION COEFFICIENTS'

      IF (CCS) THEN
         CALL AROUND( BLANKS//'FINAL CCS'//STRING(1:45)//BLANKS ) 
      ELSE IF (CC2) THEN
         CALL AROUND( BLANKS//'FINAL CC2'//STRING(1:45)//BLANKS )
      ELSE IF (CCSD) THEN
         CALL AROUND( BLANKS//'FINAL CCSD'//STRING(1:45)//BLANKS )
      ELSE IF (CC3) THEN
         CALL AROUND( BLANKS//'FINAL CC3'//STRING(1:45)//BLANKS )
      ELSE
         CALL QUIT('CCCEXPPRT called for an unknown Coupled '//
     &             'Cluster model.')
      END IF

      WRITE(LUPRI,'(/1X,4(1X,A," operator",3X),4X,A,/,72("-"))')
     &     'A','B','C','D','d_ABCD'

      ISYOLD = 1

      DO IOPER = 1, NCROPER
         ISYMA = ISYOPR(IACROP(IOPER))
         ISYMB = ISYOPR(IBCROP(IOPER))
         ISYMC = ISYOPR(ICCROP(IOPER))
         ISYMD = ISYOPR(IDCROP(IOPER))


         IDISP = 1
         IF (MULD2H(ISYMA,ISYMB).EQ.MULD2H(ISYMC,ISYMD)) THEN
            WRITE(LUPRI,'(/1X,4(A8,I3,3X),G16.8)')
     &        LBLOPR(IACROP(IOPER)),ICCAUA(IDISP),
     &        LBLOPR(IBCROP(IOPER)),ICCAUB(IDISP),
     &        LBLOPR(ICCROP(IOPER)),ICCAUC(IDISP),
     &        LBLOPR(IDCROP(IOPER)),ICCAUD(IDISP),
     &        -EXPCOF(IDISP,IOPER)
         ELSE
           IF (IPRCHYP.GT.0) THEN
             WRITE(LUPRI,'(1X,4(A8,A3,3X),6X,A,7X)')
     &         LBLOPR(IACROP(IOPER)),' - ',
     &         LBLOPR(IBCROP(IOPER)),' - ',
     &         LBLOPR(ICCROP(IOPER)),' - ', 
     &         LBLOPR(IDCROP(IOPER)),' - ', 
     &         '---'
           END IF
         END IF

         DO IDISP = 2, NCRDISP
          IF (MULD2H(ISYMA,ISYMB).EQ.MULD2H(ISYMC,ISYMD)) THEN
            WRITE(LUPRI,'(1X,4(8X,I3,3X),G16.8)')
     &        ICCAUA(IDISP),ICCAUB(IDISP),ICCAUC(IDISP),ICCAUD(IDISP),
     &        -EXPCOF(IDISP,IOPER)
          END IF
         END DO

      END DO

      IF (GAMMA_PAR .OR. GAMMA_ORT) THEN
        WRITE(LUPRI,'(/////A,/,86("-"))')'  average                    '
      END IF

      IF (GAMMA_PAR) THEN
* calculate & print gamma_{||} for all requested orders:
       DO IDISP = 1, NCRDISP

         IF (CSYM(1:6).EQ.'ATOMIC') THEN
           GAMMA = EXPCOF(IDISP,1)
         ELSE IF (CSYM(1:6).EQ.'SPHTOP') THEN
           GAMMA = 0.6d0 * EXPCOF(IDISP,1)
           DO ICOMP = 2, 4
             GAMMA = GAMMA + 0.4d0 * EXPCOF(IDISP,ICOMP)
           END DO
         ELSE IF (CSYM(1:6).EQ.'LINEAR') THEN
           GAMMA =         3.0d0 * EXPCOF(IDISP,1)
           GAMMA = GAMMA + 8.0d0 * EXPCOF(IDISP,8)
           DO ICOMP = 2, 7
             GAMMA = GAMMA + 2.0d0 * EXPCOF(IDISP,ICOMP)
           END DO
           GAMMA = GAMMA / 15.0d0
         ELSE 
           GAMMA =         2.0d0 * EXPCOF(IDISP,1)
           GAMMA = GAMMA + 2.0d0 * EXPCOF(IDISP,8)
           GAMMA = GAMMA + 2.0d0 * EXPCOF(IDISP,15)
           DO ICOMP = 1, 21
             GAMMA = GAMMA + EXPCOF(IDISP,ICOMP)
           END DO
           GAMMA = GAMMA / 15.0d0
         END IF

         IF (IDISP.EQ.1) THEN
           WRITE(LUPRI,'(1X,A8,4X,I3,2X,3(8X,I3,2X),G16.8)')
     &        'gamma_||', ICCAUA(IDISP), ICCAUB(IDISP), 
     &                    ICCAUC(IDISP), ICCAUD(IDISP), GAMMA
         ELSE
           WRITE(LUPRI,'(1X,8X,4X,I3,2X,3(8X,I3,2X),G16.8)')
     &                    ICCAUA(IDISP), ICCAUB(IDISP), 
     &                    ICCAUC(IDISP), ICCAUD(IDISP), GAMMA
         END IF
       END DO 
      END IF ! (GAMMA_PAR) 

      IF (GAMMA_ORT) THEN
* calculate & print gamma_{_|_} for all requested orders:
       DO IDISP = 1, NCRDISP

         IF (CSYM(1:6).EQ.'ATOMIC') THEN
           GAMMA = 2.0d0*EXPCOF(IDISP,2)
     &           - 1.0d0*EXPCOF(IDISP,1)
     &           + 2.0d0*EXPCOF(IDISP,3)
         ELSE IF (CSYM(1:6).EQ.'SPHTOP') THEN
           GAMMA = 0.2d0 * EXPCOF(IDISP,1)
     &           + 0.8d0 * EXPCOF(IDISP,2)
     &           - 0.2d0 * EXPCOF(IDISP,3)
     &           - 0.2d0 * EXPCOF(IDISP,4)
         ELSE IF (CSYM(1:6).EQ.'LINEAR') THEN
           GAMMA=  1.0d0 * EXPCOF(IDISP,1)
     &           + 4.0d0 * EXPCOF(IDISP,2)
     &           - 1.0d0 * EXPCOF(IDISP,3)
     &           - 1.0d0 * EXPCOF(IDISP,4)
     &           + 4.0d0 * EXPCOF(IDISP,5)
     &           - 1.0d0 * EXPCOF(IDISP,6)
     &           - 1.0d0 * EXPCOF(IDISP,7)
     &           + 1.0d0 * EXPCOF(IDISP,8)
     &           + 5.0d0 * EXPCOF(IDISP,9)
           GAMMA = GAMMA / 15.0d0
         ELSE 
           GAMMA = 0.0d0
           DO IADD = 0, 14, 7
             GAMMA = GAMMA + 
     &                       1.0d0 * EXPCOF(IDISP,1+IADD)
     &                      +2.0d0 * EXPCOF(IDISP,2+IADD)
     &                      -0.5d0 * EXPCOF(IDISP,3+IADD)
     &                      -0.5d0 * EXPCOF(IDISP,4+IADD)
     &                      +2.0d0 * EXPCOF(IDISP,5+IADD)
     &                      -0.5d0 * EXPCOF(IDISP,6+IADD)
     &                      -0.5d0 * EXPCOF(IDISP,7+IADD)
           END DO
           GAMMA = GAMMA / 15.0d0
         END IF

         IF (IDISP.EQ.1) THEN
           WRITE(LUPRI,'(/1X,A8,4X,I3,2X,3(8X,I3,2X),G16.8)')
     &        'gamma_|_', ICCAUA(IDISP), ICCAUB(IDISP), 
     &                    ICCAUC(IDISP), ICCAUD(IDISP), GAMMA
         ELSE
           WRITE(LUPRI,'(1X,8X,4X,I3,2X,3(8X,I3,2X),G16.8)')
     &                    ICCAUA(IDISP), ICCAUB(IDISP), 
     &                    ICCAUC(IDISP), ICCAUD(IDISP), GAMMA
         END IF
       END DO 
      END IF ! (GAMMA_ORT) 

      IF (GAMMA_ORT) THEN
* calculate & print gamma_{ms} for all requested orders:
       DO IDISP = 1, NCRDISP

         IF (CSYM(1:6).EQ.'ATOMIC') THEN
           GAMMA = 3.0d0*EXPCOF(IDISP,2)
     &           - 2.0d0*EXPCOF(IDISP,1)
     &           + 3.0d0*EXPCOF(IDISP,3)
         ELSE IF (CSYM(1:6).EQ.'SPHTOP') THEN
           GAMMA =  1.0d0 * EXPCOF(IDISP,2)
     &            + 1.0d0 * EXPCOF(IDISP,3)
     &            - 2.0d0 * EXPCOF(IDISP,4)
         ELSE IF (CSYM(1:6).EQ.'LINEAR') THEN
           GAMMA =  0.5d0 * EXPCOF(IDISP,2)
     &            + 0.5d0 * EXPCOF(IDISP,3)
     &            - 1.0d0 * EXPCOF(IDISP,4)
     &            + 0.5d0 * EXPCOF(IDISP,5)
     &            + 0.5d0 * EXPCOF(IDISP,6)
     &            - 1.0d0 * EXPCOF(IDISP,7)
     &            + 1.5d0 * EXPCOF(IDISP,9)
     &            + 1.5d0 * EXPCOF(IDISP,10)
     &            - 1.0d0 * EXPCOF(IDISP,8)
           GAMMA = GAMMA / 1.5d0
         ELSE 
           GAMMA = 0.0d0
           DO IADD = 0, 14, 7
             GAMMA = GAMMA + 
     &                          0.5d0 * EXPCOF(IDISP,2+IADD)
     &                        + 0.5d0 * EXPCOF(IDISP,3+IADD)
     &                        - 1.0d0 * EXPCOF(IDISP,4+IADD)
     &                        + 0.5d0 * EXPCOF(IDISP,5+IADD)
     &                        + 0.5d0 * EXPCOF(IDISP,6+IADD)
     &                        - 1.0d0 * EXPCOF(IDISP,7+IADD)
           END DO
           GAMMA = GAMMA / 3.0d0
         END IF

         IF (IDISP.EQ.1) THEN
           WRITE(LUPRI,'(/1X,A8,4X,I3,2X,3(8X,I3,2X),G16.8)')
     &        'gamma_ms', ICCAUA(IDISP), ICCAUB(IDISP), 
     &                    ICCAUC(IDISP), ICCAUD(IDISP), GAMMA
         ELSE
           WRITE(LUPRI,'(1X,8X,4X,I3,2X,3(8X,I3,2X),G16.8)')
     &                    ICCAUA(IDISP), ICCAUB(IDISP), 
     &                    ICCAUC(IDISP), ICCAUD(IDISP), GAMMA
         END IF
       END DO 
      END IF ! (GAMMA_ORT) 


      CALL FLSHFO(LUPRI)

      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCCEXPPRT                            *
*---------------------------------------------------------------------*
c /* deck CCCRDISP */
*=====================================================================*
      SUBROUTINE CCCRDISP(DISPCF,AVEDISP,ABCDE,
     &                    THGCF, ESHGCF, KERRCF,DFWMCF,
     &                    AVETHG,AVESHG, AVEKERR, AVEDFWM,
     &                    EXPCOF,TRINOM, NDSPCOF,LERROR)
*---------------------------------------------------------------------*
*
*   Purpose: calculate from the expansion coefficients defined by
*            gamma_ABCD=sum_{klmn} w_A^k w_B^l w_C^m w_D^n d_ABCD(klmn)
*            the physical more relevant dispersion coefficients defined
*            by gamma_ABCD = sum_{lmn} w_B^l w_C^m w_D^l D_ABCD(lmn),
*            special versions thereof for THG, ESHG, Kerr & DFWM,
*            as well as coefficients for the isotropic averages and
*            the `universal' dispersion coefficients A, B, ...
*
*
*   Written by Christof Haettig in March 1998.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <cccrinf.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER KA, KB, KBP, KC, KCP, KCPP, KD, KDP, KDPP, KDPPP
      INTEGER KE, KEP, KEPP, KEPPP, KEPPPP
      PARAMETER (KA=1, KB=2, KBP=3, KC=4, KCP=5, KCPP=6)
      PARAMETER (KD=7, KDP=8, KDPP=9, KDPPP=10)
      PARAMETER (KE=11, KEP=12, KEPP=13, KEPPP=14, KEPPPP=15)

      INTEGER NDSPCOF ! leading dimension of DISPCF, AVEDISP, TRINOM
      LOGICAL LERROR

#if defined (SYS_CRAY)
      REAL DISPCF(NDSPCOF,NCROPER)
      REAL AVEDISP(NDSPCOF,4)
      REAL ABCDE(15,2)
      REAL THGCF(NCRDSPE+1,NCROPER)
      REAL ESHGCF(NCRDSPE+1,NCROPER)
      REAL KERRCF(NCRDSPE+1,NCROPER)
      REAL DFWMCF(NCRDSPE+1,NCROPER)
      REAL AVESHG(NCRDSPE+1,4)
      REAL AVETHG(NCRDSPE+1,4)
      REAL AVEKERR(NCRDSPE+1,4)
      REAL AVEDFWM(NCRDSPE+1,4)
      REAL EXPCOF(NCRDISP,NCROPER)
      REAL TRINOM(NDSPCOF)
      REAL GAMMA0,A,B,BP,C,CP,CPP
      REAL D, DP, DPP, DPPP, E, EP, EPP, EPPP, EPPPP
      REAL ERROR, ERR1, ERR2, ERR3, ERR4, ERR5, ERR6, ERR7
      REAL ERR0, ERR8, ERR9, ERR10, ERR11, ERR12, ERR13, ERR14
      REAL ERR15, ERR16
      REAL Y, X1, X2, X3, X4, Z, W
      REAL THRSDISP
#else
      DOUBLE PRECISION DISPCF(NDSPCOF,NCROPER)
      DOUBLE PRECISION AVEDISP(NDSPCOF,4)
      DOUBLE PRECISION ABCDE(15,2)
      DOUBLE PRECISION THGCF(NCRDSPE+1,NCROPER)
      DOUBLE PRECISION ESHGCF(NCRDSPE+1,NCROPER)
      DOUBLE PRECISION KERRCF(NCRDSPE+1,NCROPER)
      DOUBLE PRECISION DFWMCF(NCRDSPE+1,NCROPER)
      DOUBLE PRECISION AVESHG(NCRDSPE+1,4)
      DOUBLE PRECISION AVETHG(NCRDSPE+1,4)
      DOUBLE PRECISION AVEKERR(NCRDSPE+1,4)
      DOUBLE PRECISION AVEDFWM(NCRDSPE+1,4)
      DOUBLE PRECISION EXPCOF(NCRDISP,NCROPER)
      DOUBLE PRECISION TRINOM(NDSPCOF)
      DOUBLE PRECISION GAMMA0,A,B,BP,C,CP,CPP
      DOUBLE PRECISION D, DP, DPP, DPPP, E, EP, EPP, EPPP, EPPPP
      DOUBLE PRECISION ERROR, ERR1, ERR2, ERR3, ERR4, ERR5, ERR6, ERR7
      DOUBLE PRECISION ERR0,ERR8,ERR9,ERR10,ERR11,ERR12,ERR13,ERR14
      DOUBLE PRECISION ERR15, ERR16
      DOUBLE PRECISION Y, X1, X2, X3, X4, Z, W
      DOUBLE PRECISION THRSDISP
#endif
      PARAMETER (THRSDISP = 1.0d-5)
 
      INTEGER IDISP, ICOMP, IOPER,OFFEX, IEXPCF
      INTEGER L, N, M, LMN, MN, P, Q, R, I, J, IADD

      INTEGER ITRI
#if defined (SYS_CRAY)
      REAL SIGNEO
#else
      DOUBLE PRECISION SIGNEO
#endif

* index for 3-dim array with L >= M >= N and L,M,N >= 0 :
      ITRI(L,M,N) = L*(L+1)*(L+2)/6 + M*(M+1)/2 + N + 1

* this gives (-1)^I:
      SIGNEO(I)   = DBLE(1-MOD(I,2)*2)
      
C
*---------------------------------------------------------------------*
* initialze flag for numerical errors in A,B... coefficients:
*---------------------------------------------------------------------*
      LERROR = .FALSE.

*---------------------------------------------------------------------*
* precompute trinomial coefficients: L! / ( (L-M)! (M-N)! N! )
*---------------------------------------------------------------------*
      DO L = 0, NCRDSPE
        TRINOM(ITRI(L,0,0)) = 1.0d0
        TRINOM(ITRI(L,L,0)) = 1.0d0
        TRINOM(ITRI(L,L,L)) = 1.0d0
        DO M = 1, L-1
          TRINOM(ITRI(L,M,0)) = 
     &       TRINOM(ITRI(L-1,M,0)) + TRINOM(ITRI(L-1,M-1,0))
          TRINOM(ITRI(L,M,M)) = 
     &       TRINOM(ITRI(L-1,M,M)) + TRINOM(ITRI(L-1,M-1,M-1))
          DO N = 1, M-1
            TRINOM(ITRI(L,M,N)) = TRINOM(ITRI(L-1,M,N)) +
     &         TRINOM(ITRI(L-1,M-1,N)) + TRINOM(ITRI(L-1,M-1,N-1))
          END DO
        END DO
        DO N = 1, L-1
          TRINOM(ITRI(L,L,N)) = 
     &       TRINOM(ITRI(L-1,L-1,N)) + TRINOM(ITRI(L-1,L-1,N-1))
        END DO
      END DO

      IF (LOCDBG) THEN
        WRITE (LUPRI,'(///)')
        WRITE (LUPRI,*) 'DEBUG_CCCRDISP> trinomial coefficients:'
        DO L = 0, NCRDSPE
          DO M = 0, L
             WRITE (LUPRI,'(2I3,4X,15F8.2)') 
     &            L,M,(TRINOM(ITRI(L,M,N)),N=0,M)
          END DO
        END DO
        WRITE (LUPRI,*)
      END IF

      DO IOPER = 1, NCROPER
        OFFEX = 0
        DO LMN = 0, NCRDSPE, 2
        DO MN  = 0, LMN
        DO N   = 0, MN
          L = LMN - MN
          M = MN - N

          IF (LOCDBG) THEN
            WRITE (LUPRI,'(//A,A)')
     &      '  L  M  N  P  Q  R    SIGN   TRINOM',
     &      '   I J K L      EXPCOF        DISPCF' 
          END IF

          DISPCF(ITRI(LMN,MN,N),IOPER) = 0.0d0

          DO P = 0, L
          DO Q = 0, M
          DO R = 0, N
            IEXPCF = OFFEX + ITRI(LMN-P-Q-R,MN-Q-R,N-R)

            DISPCF(ITRI(LMN,MN,N),IOPER) = 
     &        DISPCF(ITRI(LMN,MN,N),IOPER) + SIGNEO(P+Q+R) * 
     &          TRINOM(ITRI(P+Q+R,Q+R,R)) * EXPCOF(IEXPCF,IOPER)

            IF (LOCDBG) THEN
              WRITE (LUPRI,'(6I3,2F8.2,3X,4I2,2E16.8,5I5)')
     &         L,M,N,P,Q,R,SIGNEO(P+Q+R),TRINOM(ITRI(P+Q+R,Q+R,R)),
     &          P+Q+R,L-P,M-Q,N-R, EXPCOF(IEXPCF,IOPER),
     &           DISPCF(ITRI(LMN,MN,N),IOPER),
     &            ITRI(P+Q+R,Q+R,R),ITRI(LMN,MN,N),IEXPCF,IOPER
            END IF

          END DO ! R
          END DO ! Q
          END DO ! P

        END DO ! N
        END DO ! MN

        OFFEX = OFFEX + (LMN+3)*(LMN+2)*(LMN+1)/6

        END DO ! LMN
      END DO ! IOPER


      DO IOPER = 1, NCROPER
        DO LMN = 0, NCRDSPE, 2
          IF (LOCDBG) THEN
            WRITE (LUPRI,'(//A,5(5X,A,5X))') '  L  M  N ', 
     &         'DISPCF', 'THGCF ', 'ESHGCF', 'KERRCF', 'DFWMCF'
          END IF
          THGCF(LMN+1,IOPER)  = 0.0d0
          ESHGCF(LMN+1,IOPER) = 0.0d0
          KERRCF(LMN+1,IOPER) = 0.0d0
          DFWMCF(LMN+1,IOPER) = 0.0d0

          KERRCF(LMN+1,IOPER) = KERRCF(LMN+1,IOPER) +
     &                              DISPCF(ITRI(LMN,LMN,LMN),IOPER)
          DO MN = 0, LMN
            ESHGCF(LMN+1,IOPER) = ESHGCF(LMN+1,IOPER) +
     &                                DISPCF(ITRI(LMN,MN,0),IOPER)
            DO N  = 0, MN
              THGCF(LMN+1,IOPER)  = THGCF(LMN+1,IOPER) + 
     &                                 DISPCF(ITRI(LMN,MN,N),IOPER)
       
              DFWMCF(LMN+1,IOPER) = DFWMCF(LMN+1,IOPER) +
     &                     SIGNEO(N) * DISPCF(ITRI(LMN,MN,N),IOPER)

              IF (LOCDBG) THEN
                WRITE (LUPRI,'(3I3,5E16.8)') LMN-MN,MN-N,N,
     &                 DISPCF(ITRI(LMN,MN,N),IOPER), 
     &                 THGCF(LMN+1,IOPER), ESHGCF(LMN+1,IOPER), 
     &                 KERRCF(LMN+1,IOPER),DFWMCF(LMN+1,IOPER)
              END IF

            END DO
          END DO

        END DO
      END DO

      IF (GAMMA_PAR .OR. GAMMA_ORT) THEN

* calculate averaged vector components parallel/orthogonal to z axis:
        DO IDISP = 1, (NCRDSPE+3)*(NCRDSPE+2)*(NCRDSPE+1)/6
          IF (CSYM(1:6).EQ.'ATOMIC') THEN
            IF (GAMMA_PAR) THEN
              AVEDISP(IDISP,1) = 1.0d0*DISPCF(IDISP,1)
            END IF
            IF (GAMMA_ORT) THEN
              AVEDISP(IDISP,2) = 2.0d0*DISPCF(IDISP,2)
     &                         - 1.0d0*DISPCF(IDISP,1)
     &                         + 2.0d0*DISPCF(IDISP,3)
              AVEDISP(IDISP,3) = 3.0d0*DISPCF(IDISP,2)
     &                         - 2.0d0*DISPCF(IDISP,1)
     &                         + 3.0d0*DISPCF(IDISP,3)
            END IF
          ELSE IF (CSYM(1:6).EQ.'SPHTOP') THEN
            IF (GAMMA_PAR) THEN
              AVEDISP(IDISP,1) = 0.6d0 * DISPCF(IDISP,1)
     &                         + 0.4d0 * DISPCF(IDISP,2)
     &                         + 0.4d0 * DISPCF(IDISP,3)
     &                         + 0.4d0 * DISPCF(IDISP,4)
            END IF
            IF (GAMMA_ORT) THEN
              AVEDISP(IDISP,2) = 0.2d0 * DISPCF(IDISP,1)
     &                         + 0.8d0 * DISPCF(IDISP,2)
     &                         + 0.8d0 * DISPCF(IDISP,3)
     &                         - 1.2d0 * DISPCF(IDISP,4)
              AVEDISP(IDISP,3) = 1.0d0 * DISPCF(IDISP,2)
     &                         + 1.0d0 * DISPCF(IDISP,3)
     &                         - 2.0d0 * DISPCF(IDISP,4)
            END IF
          ELSE IF (CSYM(1:6).EQ.'LINEAR') THEN
            IF (GAMMA_PAR) THEN
              AVEDISP(IDISP,1) = 3.0d0 * DISPCF(IDISP,1)
              AVEDISP(IDISP,1) = AVEDISP(IDISP,1) + 
     &                                    8.0d0 * DISPCF(IDISP,8)
              DO ICOMP = 2, 7
                AVEDISP(IDISP,1) = AVEDISP(IDISP,1) + 
     &                                  2.0d0 * DISPCF(IDISP,ICOMP)
              END DO
              AVEDISP(IDISP,1) = AVEDISP(IDISP,1) / 15.0d0
            END IF
            IF (GAMMA_ORT) THEN
              AVEDISP(IDISP,2) =  1.0d0 * DISPCF(IDISP,1)
     &                          + 4.0d0 * DISPCF(IDISP,2)
     &                          + 4.0d0 * DISPCF(IDISP,3)
     &                          - 6.0d0 * DISPCF(IDISP,4)
     &                          + 4.0d0 * DISPCF(IDISP,5)
     &                          + 4.0d0 * DISPCF(IDISP,6)
     &                          - 6.0d0 * DISPCF(IDISP,7)
     &                          + 10.0d0 * DISPCF(IDISP,9)
     &                          + 10.0d0 * DISPCF(IDISP,10)
     &                          - 4.0d0 * DISPCF(IDISP,8)
              AVEDISP(IDISP,2) = AVEDISP(IDISP,2) / 15.0d0
              AVEDISP(IDISP,3) =  0.5d0 * DISPCF(IDISP,2)
     &                          + 0.5d0 * DISPCF(IDISP,3)
     &                          - 1.0d0 * DISPCF(IDISP,4)
     &                          + 0.5d0 * DISPCF(IDISP,5)
     &                          + 0.5d0 * DISPCF(IDISP,6)
     &                          - 1.0d0 * DISPCF(IDISP,7)
     &                          + 1.5d0 * DISPCF(IDISP,9)
     &                          + 1.5d0 * DISPCF(IDISP,10)
     &                          - 1.0d0 * DISPCF(IDISP,8)
              AVEDISP(IDISP,3) = AVEDISP(IDISP,3) / 1.5d0
            END IF
          ELSE
            IF (GAMMA_PAR) THEN
              AVEDISP(IDISP,1) =                  2.0d0*DISPCF(IDISP,1)
              AVEDISP(IDISP,1) = AVEDISP(IDISP,1)+2.0d0*DISPCF(IDISP,8)
              AVEDISP(IDISP,1) =AVEDISP(IDISP,1)+2.0d0*DISPCF(IDISP,15)
              DO ICOMP = 1, 21
                AVEDISP(IDISP,1)=AVEDISP(IDISP,1) + DISPCF(IDISP,ICOMP)
              END DO
              AVEDISP(IDISP,1) = AVEDISP(IDISP,1) / 15.0d0
            END IF
            IF (GAMMA_ORT) THEN
              AVEDISP(IDISP,2) = 0.0d0
              DO IADD = 0, 14, 7
                AVEDISP(IDISP,2) = AVEDISP(IDISP,2) +
     &                          1.0d0 * DISPCF(IDISP,1+IADD)
     &                         +2.0d0 * DISPCF(IDISP,2+IADD)
     &                         +2.0d0 * DISPCF(IDISP,3+IADD)
     &                         -3.0d0 * DISPCF(IDISP,4+IADD)
     &                         +2.0d0 * DISPCF(IDISP,5+IADD)
     &                         -2.0d0 * DISPCF(IDISP,6+IADD)
     &                         -3.0d0 * DISPCF(IDISP,7+IADD)
              END DO
              AVEDISP(IDISP,2) = AVEDISP(IDISP,2) / 15.0d0
              AVEDISP(IDISP,3) = 0.0d0
              DO IADD = 0, 14, 7
                AVEDISP(IDISP,3) = AVEDISP(IDISP,3) +
     &                             0.5d0 * DISPCF(IDISP,2+IADD)
     &                           + 0.5d0 * DISPCF(IDISP,3+IADD)
     &                           - 1.0d0 * DISPCF(IDISP,4+IADD)
     &                           + 0.5d0 * DISPCF(IDISP,5+IADD)
     &                           + 0.5d0 * DISPCF(IDISP,6+IADD)
     &                           - 1.0d0 * DISPCF(IDISP,7+IADD)
              END DO
              AVEDISP(IDISP,3) = AVEDISP(IDISP,3) / 3.0d0
            END IF
          END IF
        END DO

        DO IDISP = 1, NCRDSPE+1
          IF (CSYM(1:6).EQ.'ATOMIC') THEN
            IF (GAMMA_PAR) THEN
              AVETHG(IDISP,1)  = THGCF(IDISP,1)
              AVESHG(IDISP,1)  = ESHGCF(IDISP,1)
              AVEKERR(IDISP,1) = KERRCF(IDISP,1)
              AVEDFWM(IDISP,1) = DFWMCF(IDISP,1)
            END IF
            IF (GAMMA_ORT) THEN
              AVETHG(IDISP,2)  = 2.0d0*THGCF(IDISP,2)
     &                         - 1.0d0*THGCF(IDISP,1)
     &                         + 2.0d0*THGCF(IDISP,3)
              AVETHG(IDISP,3)  = 3.0d0*THGCF(IDISP,2)
     &                         - 2.0d0*THGCF(IDISP,1)
     &                         + 3.0d0*THGCF(IDISP,3)
              AVESHG(IDISP,2)  = 2.0d0*ESHGCF(IDISP,2)
     &                         - 1.0d0*ESHGCF(IDISP,1)
     &                         + 2.0d0*ESHGCF(IDISP,3)
              AVESHG(IDISP,3)  = 3.0d0*ESHGCF(IDISP,2)
     &                         - 2.0d0*ESHGCF(IDISP,1)
     &                         + 3.0d0*ESHGCF(IDISP,3)
              AVEKERR(IDISP,2) = 2.0d0*KERRCF(IDISP,2)
     &                         - 1.0d0*KERRCF(IDISP,1)
     &                         + 2.0d0*KERRCF(IDISP,3)
              AVEKERR(IDISP,3) = 3.0d0*KERRCF(IDISP,2)
     &                         - 2.0d0*KERRCF(IDISP,1)
     &                         + 3.0d0*KERRCF(IDISP,3)
              AVEDFWM(IDISP,2) = 2.0d0*DFWMCF(IDISP,2)
     &                         - 1.0d0*DFWMCF(IDISP,1)
     &                         + 2.0d0*DFWMCF(IDISP,3)
              AVEDFWM(IDISP,3) = 3.0d0*DFWMCF(IDISP,2)
     &                         - 2.0d0*DFWMCF(IDISP,1)
     &                         + 3.0d0*DFWMCF(IDISP,3)
            END IF
          ELSE IF (CSYM(1:6).EQ.'SPHTOP') THEN
            IF (GAMMA_PAR) THEN
              AVETHG(IDISP,1) = 0.6d0 * THGCF(IDISP,1)
              AVESHG(IDISP,1) = 0.6d0 * ESHGCF(IDISP,1)
              AVEKERR(IDISP,1) = 0.6d0 * KERRCF(IDISP,1)
              AVEDFWM(IDISP,1) = 0.6d0 * DFWMCF(IDISP,1)
              DO ICOMP = 2, 4
                AVETHG(IDISP,1) = AVETHG(IDISP,1) + 
     &                                    0.4d0 * THGCF(IDISP,ICOMP)
                AVESHG(IDISP,1) = AVESHG(IDISP,1) + 
     &                                    0.4d0 * ESHGCF(IDISP,ICOMP)
                AVEKERR(IDISP,1) = AVEKERR(IDISP,1) + 
     &                                    0.4d0 * KERRCF(IDISP,ICOMP)
                AVEDFWM(IDISP,1) = AVEDFWM(IDISP,1) + 
     &                                    0.4d0 * DFWMCF(IDISP,ICOMP)
              END DO
            END IF
            IF (GAMMA_ORT) THEN
              AVETHG(IDISP,2)  = 0.2d0 * THGCF(IDISP,1)
     &                         + 0.8d0 * THGCF(IDISP,2)
     &                         - 0.2d0 * THGCF(IDISP,3)
     &                         - 0.2d0 * THGCF(IDISP,4)
              AVETHG(IDISP,3)  = 1.0d0 * THGCF(IDISP,2)
     &                         + 1.0d0 * THGCF(IDISP,3)
     &                         - 2.0d0 * THGCF(IDISP,4)
              AVESHG(IDISP,2)  = 0.2d0 * ESHGCF(IDISP,1)
     &                         + 0.8d0 * ESHGCF(IDISP,2)
     &                         - 0.2d0 * ESHGCF(IDISP,3)
     &                         - 0.2d0 * ESHGCF(IDISP,4)
              AVESHG(IDISP,3)  = 1.0d0 * ESHGCF(IDISP,2)
     &                         + 1.0d0 * ESHGCF(IDISP,3)
     &                         - 2.0d0 * ESHGCF(IDISP,4)
              AVEKERR(IDISP,2) = 0.2d0 * KERRCF(IDISP,1)
     &                         + 0.8d0 * KERRCF(IDISP,2)
     &                         - 0.2d0 * KERRCF(IDISP,3)
     &                         - 0.2d0 * KERRCF(IDISP,4)
              AVEKERR(IDISP,3) = 1.0d0 * KERRCF(IDISP,2)
     &                         + 1.0d0 * KERRCF(IDISP,3)
     &                         - 2.0d0 * KERRCF(IDISP,4)
              AVEDFWM(IDISP,2) = 0.2d0 * DFWMCF(IDISP,1)
     &                         + 0.8d0 * DFWMCF(IDISP,2)
     &                         - 0.2d0 * DFWMCF(IDISP,3)
     &                         - 0.2d0 * DFWMCF(IDISP,4)
              AVEDFWM(IDISP,3) = 1.0d0 * DFWMCF(IDISP,2)
     &                         + 1.0d0 * DFWMCF(IDISP,3)
     &                         - 2.0d0 * DFWMCF(IDISP,4)
            END IF
          ELSE IF (CSYM(1:6).EQ.'LINEAR') THEN
            IF (GAMMA_PAR) THEN
              AVETHG(IDISP,1) = 3.0d0 * THGCF(IDISP,1)
              AVESHG(IDISP,1) = 3.0d0 * ESHGCF(IDISP,1)
              AVEKERR(IDISP,1) = 3.0d0 * KERRCF(IDISP,1)
              AVEDFWM(IDISP,1) = 3.0d0 * DFWMCF(IDISP,1)

              AVETHG(IDISP,1) = AVETHG(IDISP,1) + 
     &                                    8.0d0 * THGCF(IDISP,8)
              AVESHG(IDISP,1) = AVESHG(IDISP,1) + 
     &                                    8.0d0 * ESHGCF(IDISP,8)
              AVEKERR(IDISP,1) = AVEKERR(IDISP,1) + 
     &                                    8.0d0 * KERRCF(IDISP,8)
              AVEDFWM(IDISP,1) = AVEDFWM(IDISP,1) + 
     &                                    8.0d0 * DFWMCF(IDISP,8)

              DO ICOMP = 2, 7
                AVETHG(IDISP,1) = AVETHG(IDISP,1) + 
     &                                    2.0d0 * THGCF(IDISP,ICOMP)
                AVESHG(IDISP,1) = AVESHG(IDISP,1) + 
     &                                    2.0d0 * ESHGCF(IDISP,ICOMP)
                AVEKERR(IDISP,1) = AVEKERR(IDISP,1) + 
     &                                    2.0d0 * KERRCF(IDISP,ICOMP)
                AVEDFWM(IDISP,1) = AVEDFWM(IDISP,1) + 
     &                                    2.0d0 * DFWMCF(IDISP,ICOMP)
              END DO

              AVETHG(IDISP,1) = AVETHG(IDISP,1) / 15.0d0
              AVESHG(IDISP,1) = AVESHG(IDISP,1) / 15.0d0
              AVEKERR(IDISP,1) = AVEKERR(IDISP,1) / 15.0d0
              AVEDFWM(IDISP,1) = AVEDFWM(IDISP,1) / 15.0d0
            END IF
            IF (GAMMA_ORT) THEN
              AVETHG(IDISP,2)  =  1.0d0 * THGCF(IDISP,1)
     &                          + 4.0d0 * THGCF(IDISP,2)
     &                          - 1.0d0 * THGCF(IDISP,3)
     &                          - 1.0d0 * THGCF(IDISP,4)
     &                          + 4.0d0 * THGCF(IDISP,5)
     &                          - 1.0d0 * THGCF(IDISP,6)
     &                          - 1.0d0 * THGCF(IDISP,7)
     &                          + 1.0d0 * THGCF(IDISP,8)
     &                          + 5.0d0 * THGCF(IDISP,9)
              AVETHG(IDISP,2)  = AVETHG(IDISP,2) / 15.0d0
              AVETHG(IDISP,3)  =  0.5d0 * THGCF(IDISP,2)
     &                          + 0.5d0 * THGCF(IDISP,3)
     &                          - 1.0d0 * THGCF(IDISP,4)
     &                          + 0.5d0 * THGCF(IDISP,5)
     &                          + 0.5d0 * THGCF(IDISP,6)
     &                          - 1.0d0 * THGCF(IDISP,7)
     &                          + 1.5d0 * THGCF(IDISP,9)
     &                          + 1.5d0 * THGCF(IDISP,10)
     &                          - 1.0d0 * THGCF(IDISP,8)
              AVETHG(IDISP,3)  = AVETHG(IDISP,3) / 1.5d0
              AVESHG(IDISP,2)  =  1.0d0 * ESHGCF(IDISP,1)
     &                          + 4.0d0 * ESHGCF(IDISP,2)
     &                          - 1.0d0 * ESHGCF(IDISP,3)
     &                          - 1.0d0 * ESHGCF(IDISP,4)
     &                          + 4.0d0 * ESHGCF(IDISP,5)
     &                          - 1.0d0 * ESHGCF(IDISP,6)
     &                          - 1.0d0 * ESHGCF(IDISP,7)
     &                          + 1.0d0 * ESHGCF(IDISP,8)
     &                          + 5.0d0 * ESHGCF(IDISP,9)
              AVESHG(IDISP,2)  = AVESHG(IDISP,2) / 15.0d0
              AVESHG(IDISP,3)  =  0.5d0 * ESHGCF(IDISP,2)
     &                          + 0.5d0 * ESHGCF(IDISP,3)
     &                          - 1.0d0 * ESHGCF(IDISP,4)
     &                          + 0.5d0 * ESHGCF(IDISP,5)
     &                          + 0.5d0 * ESHGCF(IDISP,6)
     &                          - 1.0d0 * ESHGCF(IDISP,7)
     &                          + 1.5d0 * ESHGCF(IDISP,9)
     &                          + 1.5d0 * ESHGCF(IDISP,10)
     &                          - 1.0d0 * ESHGCF(IDISP,8)
              AVESHG(IDISP,3)  = AVESHG(IDISP,3) / 1.5d0
              AVEKERR(IDISP,2) =  1.0d0 * KERRCF(IDISP,1)
     &                          + 4.0d0 * KERRCF(IDISP,2)
     &                          - 1.0d0 * KERRCF(IDISP,3)
     &                          - 1.0d0 * KERRCF(IDISP,4)
     &                          + 4.0d0 * KERRCF(IDISP,5)
     &                          - 1.0d0 * KERRCF(IDISP,6)
     &                          - 1.0d0 * KERRCF(IDISP,7)
     &                          + 1.0d0 * KERRCF(IDISP,8)
     &                          + 5.0d0 * KERRCF(IDISP,9)
              AVEKERR(IDISP,2) = AVEKERR(IDISP,2) / 15.0d0
              AVEKERR(IDISP,3) =  0.5d0 * KERRCF(IDISP,2)
     &                          + 0.5d0 * KERRCF(IDISP,3)
     &                          - 1.0d0 * KERRCF(IDISP,4)
     &                          + 0.5d0 * KERRCF(IDISP,5)
     &                          + 0.5d0 * KERRCF(IDISP,6)
     &                          - 1.0d0 * KERRCF(IDISP,7)
     &                          + 1.5d0 * KERRCF(IDISP,9)
     &                          + 1.5d0 * KERRCF(IDISP,10)
     &                          - 1.0d0 * KERRCF(IDISP,8)
              AVEKERR(IDISP,3) = AVEKERR(IDISP,3) / 1.5d0
              AVEDFWM(IDISP,2) =  1.0d0 * DFWMCF(IDISP,1)
     &                          + 4.0d0 * DFWMCF(IDISP,2)
     &                          - 1.0d0 * DFWMCF(IDISP,3)
     &                          - 1.0d0 * DFWMCF(IDISP,4)
     &                          + 4.0d0 * DFWMCF(IDISP,5)
     &                          - 1.0d0 * DFWMCF(IDISP,6)
     &                          - 1.0d0 * DFWMCF(IDISP,7)
     &                          + 1.0d0 * DFWMCF(IDISP,8)
     &                          + 5.0d0 * DFWMCF(IDISP,9)
              AVEDFWM(IDISP,2) = AVEDFWM(IDISP,2) / 15.0d0
              AVEDFWM(IDISP,3) =  0.5d0 * DFWMCF(IDISP,2)
     &                          + 0.5d0 * DFWMCF(IDISP,3)
     &                          - 1.0d0 * DFWMCF(IDISP,4)
     &                          + 0.5d0 * DFWMCF(IDISP,5)
     &                          + 0.5d0 * DFWMCF(IDISP,6)
     &                          - 1.0d0 * DFWMCF(IDISP,7)
     &                          + 1.5d0 * DFWMCF(IDISP,9)
     &                          + 1.5d0 * DFWMCF(IDISP,10)
     &                          - 1.0d0 * DFWMCF(IDISP,8)
              AVEDFWM(IDISP,3) = AVEDFWM(IDISP,3) / 1.5d0
            END IF
          ELSE
            IF (GAMMA_PAR) THEN
              AVETHG(IDISP,1) =                 2.0d0*THGCF(IDISP,1)
              AVESHG(IDISP,1) =                 2.0d0*ESHGCF(IDISP,1)
              AVEKERR(IDISP,1)=                 2.0d0*KERRCF(IDISP,1)
              AVEDFWM(IDISP,1)=                 2.0d0*DFWMCF(IDISP,1)

              AVETHG(IDISP,1) =AVETHG(IDISP,1) +2.0d0*THGCF(IDISP,8)
              AVESHG(IDISP,1) =AVESHG(IDISP,1) +2.0d0*ESHGCF(IDISP,8)
              AVEKERR(IDISP,1)=AVEKERR(IDISP,1)+2.0d0*KERRCF(IDISP,8)
              AVEDFWM(IDISP,1)=AVEDFWM(IDISP,1)+2.0d0*DFWMCF(IDISP,8)

              AVETHG(IDISP,1) =AVETHG(IDISP,1) +2.0d0*THGCF(IDISP,15)
              AVESHG(IDISP,1) =AVESHG(IDISP,1) +2.0d0*ESHGCF(IDISP,15)
              AVEKERR(IDISP,1)=AVEKERR(IDISP,1)+2.0d0*KERRCF(IDISP,15)
              AVEDFWM(IDISP,1)=AVEDFWM(IDISP,1)+2.0d0*DFWMCF(IDISP,15)

              DO ICOMP = 1, 21
                AVETHG(IDISP,1) =AVETHG(IDISP,1)  + THGCF(IDISP,ICOMP)
                AVESHG(IDISP,1) =AVESHG(IDISP,1)  + ESHGCF(IDISP,ICOMP)
                AVEKERR(IDISP,1)=AVEKERR(IDISP,1) + KERRCF(IDISP,ICOMP)
                AVEDFWM(IDISP,1)=AVEDFWM(IDISP,1) + DFWMCF(IDISP,ICOMP)
              END DO

              AVETHG(IDISP,1)  = AVETHG(IDISP,1) / 15.0d0
              AVESHG(IDISP,1)  = AVESHG(IDISP,1) / 15.0d0
              AVEKERR(IDISP,1) = AVEKERR(IDISP,1) / 15.0d0
              AVEDFWM(IDISP,1) = AVEDFWM(IDISP,1) / 15.0d0
            END IF
            IF (GAMMA_ORT) THEN
              AVETHG(IDISP,2) = 0.0d0
              DO IADD = 0, 14, 7
                AVETHG(IDISP,2) = AVETHG(IDISP,2) +
     &                          1.0d0 * THGCF(IDISP,1+IADD)
     &                         +2.0d0 * THGCF(IDISP,2+IADD)
     &                         -0.5d0 * THGCF(IDISP,3+IADD)
     &                         -0.5d0 * THGCF(IDISP,4+IADD)
     &                         +2.0d0 * THGCF(IDISP,5+IADD)
     &                         -0.5d0 * THGCF(IDISP,6+IADD)
     &                         -0.5d0 * THGCF(IDISP,7+IADD)
              END DO
              AVETHG(IDISP,2) = AVETHG(IDISP,2) / 15.0d0
              AVETHG(IDISP,3) = 0.0d0
              DO IADD = 0, 14, 7
                AVETHG(IDISP,3) = AVETHG(IDISP,3) +
     &                             0.5d0 * THGCF(IDISP,2+IADD)
     &                           + 0.5d0 * THGCF(IDISP,3+IADD)
     &                           - 1.0d0 * THGCF(IDISP,4+IADD)
     &                           + 0.5d0 * THGCF(IDISP,5+IADD)
     &                           + 0.5d0 * THGCF(IDISP,6+IADD)
     &                           - 1.0d0 * THGCF(IDISP,7+IADD)
              END DO
              AVETHG(IDISP,3) = AVETHG(IDISP,3) / 3.0d0

              AVESHG(IDISP,2) = 0.0d0
              DO IADD = 0, 14, 7
                AVESHG(IDISP,2) = AVESHG(IDISP,2) +
     &                          1.0d0 * ESHGCF(IDISP,1+IADD)
     &                         +2.0d0 * ESHGCF(IDISP,2+IADD)
     &                         -0.5d0 * ESHGCF(IDISP,3+IADD)
     &                         -0.5d0 * ESHGCF(IDISP,4+IADD)
     &                         +2.0d0 * ESHGCF(IDISP,5+IADD)
     &                         -0.5d0 * ESHGCF(IDISP,6+IADD)
     &                         -0.5d0 * ESHGCF(IDISP,7+IADD)
              END DO
              AVESHG(IDISP,2) = AVESHG(IDISP,2) / 15.0d0
              AVESHG(IDISP,3) = 0.0d0
              DO IADD = 0, 14, 7
                AVESHG(IDISP,3) = AVESHG(IDISP,3) +
     &                             0.5d0 * ESHGCF(IDISP,2+IADD)
     &                           + 0.5d0 * ESHGCF(IDISP,3+IADD)
     &                           - 1.0d0 * ESHGCF(IDISP,4+IADD)
     &                           + 0.5d0 * ESHGCF(IDISP,5+IADD)
     &                           + 0.5d0 * ESHGCF(IDISP,6+IADD)
     &                           - 1.0d0 * ESHGCF(IDISP,7+IADD)
              END DO
              AVESHG(IDISP,3) = AVESHG(IDISP,3) / 3.0d0

              AVEKERR(IDISP,2) = 0.0d0
              DO IADD = 0, 14, 7
                AVEKERR(IDISP,2) = AVEKERR(IDISP,2) +
     &                          1.0d0 * KERRCF(IDISP,1+IADD)
     &                         +2.0d0 * KERRCF(IDISP,2+IADD)
     &                         -0.5d0 * KERRCF(IDISP,3+IADD)
     &                         -0.5d0 * KERRCF(IDISP,4+IADD)
     &                         +2.0d0 * KERRCF(IDISP,5+IADD)
     &                         -0.5d0 * KERRCF(IDISP,6+IADD)
     &                         -0.5d0 * KERRCF(IDISP,7+IADD)
              END DO
              AVEKERR(IDISP,2) = AVEKERR(IDISP,2) / 15.0d0
              AVEKERR(IDISP,3) = 0.0d0
              DO IADD = 0, 14, 7
                AVEKERR(IDISP,3) = AVEKERR(IDISP,3) +
     &                             0.5d0 * KERRCF(IDISP,2+IADD)
     &                           + 0.5d0 * KERRCF(IDISP,3+IADD)
     &                           - 1.0d0 * KERRCF(IDISP,4+IADD)
     &                           + 0.5d0 * KERRCF(IDISP,5+IADD)
     &                           + 0.5d0 * KERRCF(IDISP,6+IADD)
     &                           - 1.0d0 * KERRCF(IDISP,7+IADD)
              END DO
              AVEKERR(IDISP,3) = AVEKERR(IDISP,3) / 3.0d0

              AVEDFWM(IDISP,2) = 0.0d0
              DO IADD = 0, 14, 7
                AVEDFWM(IDISP,2) = AVEDFWM(IDISP,2) +
     &                          1.0d0 * DFWMCF(IDISP,1+IADD)
     &                         +2.0d0 * DFWMCF(IDISP,2+IADD)
     &                         -0.5d0 * DFWMCF(IDISP,3+IADD)
     &                         -0.5d0 * DFWMCF(IDISP,4+IADD)
     &                         +2.0d0 * DFWMCF(IDISP,5+IADD)
     &                         -0.5d0 * DFWMCF(IDISP,6+IADD)
     &                         -0.5d0 * DFWMCF(IDISP,7+IADD)
              END DO
              AVEDFWM(IDISP,2) = AVEDFWM(IDISP,2) / 15.0d0
              AVEDFWM(IDISP,3) = 0.0d0
              DO IADD = 0, 14, 7
                AVEDFWM(IDISP,3) = AVEDFWM(IDISP,3) +
     &                             0.5d0 * DFWMCF(IDISP,2+IADD)
     &                           + 0.5d0 * DFWMCF(IDISP,3+IADD)
     &                           - 1.0d0 * DFWMCF(IDISP,4+IADD)
     &                           + 0.5d0 * DFWMCF(IDISP,5+IADD)
     &                           + 0.5d0 * DFWMCF(IDISP,6+IADD)
     &                           - 1.0d0 * DFWMCF(IDISP,7+IADD)
              END DO
              AVEDFWM(IDISP,3) = AVEDFWM(IDISP,3) / 3.0d0
            END IF
          END IF
        END DO

        IF (LOCDBG) THEN
          WRITE (LUPRI,'(//,A)')
     &    'DISPERSION COEFF. FOR PARALLEL/ORTHOGONAL/MS AVERAGE:'
          DO LMN = 0, NCRDSPE
          DO MN  = 0, LMN
          DO N   = 0, MN
            L = LMN - MN
            M = MN  - N
            WRITE (LUPRI,'(3I5,4G16.8)')
     &           L,M,N,(AVEDISP(ITRI(LMN,MN,N),I),I=1,3)
          END DO
          END DO
          END DO

          WRITE (LUPRI,'(//,A)')
     &       'THG COEFFICIENTS FOR PARALLEL/ORTHOGONAL/MS AVERAGE:'
          DO IDISP = 1, NCRDSPE+1
            WRITE (LUPRI,'(I5,4G16.8)') IDISP-1,(AVETHG(IDISP,I),I=1,3)
          END DO

          WRITE (LUPRI,'(//,A)')
     &       'ESHG COEFFICIENTS FOR PARALLEL/ORTHOGONAL/MS AVERAGE:'
          DO IDISP = 1, NCRDSPE+1
            WRITE (LUPRI,'(I5,4G16.8)') IDISP-1,(AVESHG(IDISP,I),I=1,3)
          END DO

          WRITE (LUPRI,'(//,A)')
     &      'KERR COEFFICIENTS FOR PARALLEL/ORTHOGONAL/MS AVERAGE:'
          DO IDISP = 1, NCRDSPE+1
            WRITE (LUPRI,'(I5,4G16.8)') IDISP-1,(AVEKERR(IDISP,I),I=1,3)
          END DO

          WRITE (LUPRI,'(//,A)')
     &     'DFWM COEFFICIENTS FOR PARALLEL/ORTHOGONAL/MS AVERAGE:'
          DO IDISP = 1, NCRDSPE+1
            WRITE (LUPRI,'(I5,4G16.8)') IDISP-1,(AVEDFWM(IDISP,I),I=1,3)
          END DO
        END IF

* calculate A, B, C, etc. coefficients:
        GAMMA0 = AVEDISP(ITRI(0,0,0),1)

        IF (NCRDSPE.GE.2) THEN

          A = AVEDISP(ITRI(2,0,0),1) / (2.0d0 * GAMMA0)
          ERR1 = AVEDISP(ITRI(2,1,0),1) / (2.0d0 * GAMMA0) - A
          ERROR = DABS(ERR1) / DABS(A)
          IF (ERROR .GT. THRSDISP) THEN
            WRITE (LUPRI,*) 
     &                'WARNING: LARGE NUMERICAL (?) ERRORS DURING',
     &                ' CALCULATION OF A COEFFICIENT:'
            WRITE (LUPRI,*) 'TEST #1 GAVE:',ERR1
            LERROR = .TRUE.
          END IF
          ABCDE(KA,1) = A
          IF (LOCDBG) THEN
            WRITE (LUPRI,'(//,A)') 
     &            'A,B,C,D,E COEFF. FOR PARALLEL AVERAGE:'
            WRITE (LUPRI,*) " A                       coefficient :",A
          END IF

          IF (GAMMA_ORT) THEN
            A = AVEDISP(ITRI(2,0,0),3) / GAMMA0
            ERR1 = AVEDISP(ITRI(2,1,0),3) / (-2.0d0 * GAMMA0) - A
            ERROR = DABS(ERR1) / DABS(A)
            IF (ERROR .GT. THRSDISP) THEN
              WRITE (LUPRI,*) 
     &              'WARNING: LARGE NUMERICAL (?) ERRORS DURING',
     &                        ' CALCULATION OF A_ms COEFFICIENT:'
              WRITE (LUPRI,*) 'TEST #1 GAVE:',ERR1
              LERROR = .TRUE.
            END IF
            ABCDE(KA,2) = A
            IF (LOCDBG) THEN
              WRITE (LUPRI,*) " A_ms                    coefficient :",A
            END IF
          END IF
        END IF

        IF (NCRDSPE.GE.4) THEN
          B  = 0.25d0*AVEDISP(ITRI(4,0,0),1) 
          BP=0.25d0*AVEDISP(ITRI(4,2,1),1)-0.5d0*AVEDISP(ITRI(4,1,0),1)
          ERR1=AVEDISP(ITRI(4,1,0),1) -  8.0d0 * B
          ERR2=AVEDISP(ITRI(4,2,0),1) - 12.0d0 * B
          ERROR = (DABS(ERR1)+DABS(ERR2)) / DABS(2.0d0*B)
          B  = B  / GAMMA0
          BP = BP / GAMMA0
          ABCDE(KB,1)  = B 
          ABCDE(KBP,1) = BP
          IF (ERROR .GT. THRSDISP) THEN
            WRITE (LUPRI,*) 
     &                'WARNING: LARGE NUMERICAL (?) ERRORS DURING',
     &                ' CALCULATION OF B COEFFICIENTS:'
            WRITE (LUPRI,*) 'TEST #1 GAVE:',ERR1 
            WRITE (LUPRI,*) 'TEST #2 GAVE:',ERR2  
            LERROR = .TRUE.
          END IF
          IF (LOCDBG) THEN
           WRITE(LUPRI,*)" B, B'                   coefficients:",B,BP
          END IF


          IF (GAMMA_ORT) THEN
C           B  = AVEDISP(ITRI(4,2,2),3) / 24.0d0 
C           BP = 2.0d0 * B - 0.5d0 * AVEDISP(ITRI(4,0,0),3)
C           ERR1 = AVEDISP(ITRI(4,0,0),3) - (  4.0d0*B - 2.0d0*BP)
C           ERR2 = AVEDISP(ITRI(4,1,0),3) - ( -4.0d0*B - 4.0d0*BP)
C           ERR3 = AVEDISP(ITRI(4,1,1),3) - ( 20.0d0*B - 4.0d0*BP)
C           ERR4 = AVEDISP(ITRI(4,2,0),3) - (-12.0d0*B - 0.0d0*BP)
C           ERR5 = AVEDISP(ITRI(4,2,1),3) - (  4.0d0*B - 4.0d0*BP)
C           ERR6 = AVEDISP(ITRI(4,2,2),3) - ( 24.0d0*B - 0.0d0*BP)
C           ERR7 = AVEDISP(ITRI(4,3,0),3) - (-16.0d0*B + 8.0d0*BP)
C           ERR8 = AVEDISP(ITRI(4,3,1),3) - ( -8.0d0*B + 8.0d0*BP)
C           ERR9 = AVEDISP(ITRI(4,4,0),3) - ( -8.0d0*B + 4.0d0*BP)
            B  = AVEDISP(ITRI(4,0,0),3) / 2.0d0 
            BP =-(AVEDISP(ITRI(4,1,0),3)+AVEDISP(ITRI(4,0,0),3))/9.0d0 
            ERR1 = AVEDISP(ITRI(4,0,0),3) - (  2.0d0*B + 0.0d0*BP)
            ERR2 = AVEDISP(ITRI(4,1,0),3) - ( -2.0d0*B - 9.0d0*BP)
            ERR3 = AVEDISP(ITRI(4,1,1),3) - ( 10.0d0*B + 9.0d0*BP)
            ERR4 = AVEDISP(ITRI(4,2,0),3) - ( -6.0d0*B - 9.0d0*BP)
            ERR5 = AVEDISP(ITRI(4,2,1),3) - (  2.0d0*B - 3.0d0*BP)
            ERR6 = AVEDISP(ITRI(4,2,2),3) - ( 12.0d0*B +18.0d0*BP)
            ERR7 = AVEDISP(ITRI(4,3,0),3) - ( -8.0d0*B + 0.0d0*BP)
            ERR8 = AVEDISP(ITRI(4,3,1),3) - ( -4.0d0*B + 6.0d0*BP)
            ERR9 = AVEDISP(ITRI(4,4,0),3) - ( -4.0d0*B + 0.0d0*BP)
            ERROR = (DABS(ERR1) + DABS(ERR2) + DABS(ERR3) + DABS(ERR4)+
     &               DABS(ERR5) + DABS(ERR6) + DABS(ERR7) + DABS(ERR8)+
     &               DABS(ERR9) ) / DABS(9.0d0*B)
            IF (ERROR .GT. THRSDISP) THEN
              WRITE(LUPRI,*) 
     &              'WARNING: LARGE NUMERICAL (?) ERRORS DURING',
     &              ' CALCULATION OF B_ms COEFFICIENTS:'
              WRITE (LUPRI,*) 'TEST #1 GAVE:',ERR1 
              WRITE (LUPRI,*) 'TEST #2 GAVE:',ERR2  
              WRITE (LUPRI,*) 'TEST #3 GAVE:',ERR3   
              WRITE (LUPRI,*) 'TEST #4 GAVE:',ERR4   
              WRITE (LUPRI,*) 'TEST #5 GAVE:',ERR5   
              WRITE (LUPRI,*) 'TEST #6 GAVE:',ERR6   
              WRITE (LUPRI,*) 'TEST #7 GAVE:',ERR7   
              WRITE (LUPRI,*) 'TEST #8 GAVE:',ERR8   
              WRITE (LUPRI,*) 'TEST #9 GAVE:',ERR9   
              LERROR = .TRUE.
            END IF
            B  = B  / GAMMA0
            BP = BP / GAMMA0
            ABCDE(KB,2)  = B
            ABCDE(KBP,2) = BP
            IF (LOCDBG) THEN
              WRITE (LUPRI,*) 
     &              " B_ms, B_ms'             coefficients:",B,BP
            END IF
          END IF
        END IF

        IF (NCRDSPE.GE.6) THEN
          C   = AVEDISP(ITRI(6,0,0),1) / 8.0d0
          CP  =(AVEDISP(ITRI(6,2,0),1) 
     &                - 6.0d0*AVEDISP(ITRI(6,0,0),1) ) /  9.0d0
          CPP =(AVEDISP(ITRI(6,2,1),1) -2.0d0*AVEDISP(ITRI(6,2,0),1)
     &                + AVEDISP(ITRI(6,1,0),1) ) /  8.0d0
          ERR1 = AVEDISP(ITRI(6,1,0),1)-24.0d0*C
          ERR2 = AVEDISP(ITRI(6,3,0),1)-56.0d0*C-18.0d0*CP
          ERR3 = AVEDISP(ITRI(6,3,1),1)-120.0d0*C-16.0d0*CPP-54.0d0*CP
          ERR4 = AVEDISP(ITRI(6,4,2),1)-168.0d0*C-24.0d0*CPP-90.0d0*CP
          C   = C   /  GAMMA0
          CP  = CP  /  GAMMA0
          CPP = CPP /  GAMMA0
          ERROR = ( DABS(ERR1) + DABS(ERR2) + DABS(ERR3) 
     &                         + DABS(ERR4) ) / DABS(4.0d0*C)
          IF (ERROR .GT. THRSDISP) THEN
              WRITE (LUPRI,*) 
     &            'WARNING: LARGE NUMERICAL (?) ERRORS DURING',
     &            ' CALCULATION OF C COEFFICIENTS:'
              WRITE (LUPRI,*) 'TEST #1 GAVE:',ERR1 
              WRITE (LUPRI,*) 'TEST #2 GAVE:',ERR2  
              WRITE (LUPRI,*) 'TEST #3 GAVE:',ERR3   
              WRITE (LUPRI,*) 'TEST #4 GAVE:',ERR4   
              LERROR = .TRUE.
          END IF
          ABCDE(KC,1)   = C
          ABCDE(KCP,1)  = CP
          ABCDE(KCPP,1) = CPP
          IF (LOCDBG) THEN
            WRITE (LUPRI,*) " C, C', C''              coefficients:",
     &                 C, CP, CPP
          END IF


          IF (GAMMA_ORT) THEN
C           CPP = - AVEDISP(ITRI(6,1,0),3) / 12.0d0
C           CP  = AVEDISP(ITRI(6,3,1),3)/16.0d0 + 0.50d0 * CPP
C           C   = (AVEDISP(ITRI(6,0,0),3)-AVEDISP(ITRI(6,1,0),3)/3.0d0
C    &              -4.0d0 * CP ) / 8.0d0 
C           ERR1 = (  8.d0*C+ 4.d0*CP- 4.d0*CPP)-AVEDISP(ITRI(6,0,0),3)
C           ERR2 = (  0.d0*C+ 0.d0*CP-12.d0*CPP)-AVEDISP(ITRI(6,1,0),3)
C           ERR3 = ( 48.d0*C+24.d0*CP-12.d0*CPP)-AVEDISP(ITRI(6,1,1),3)
C           ERR4 = (-24.d0*C-12.d0*CP-12.d0*CPP)-AVEDISP(ITRI(6,2,0),3)
C           ERR5 = ( 48.d0*C+ 8.d0*CP-28.d0*CPP)-AVEDISP(ITRI(6,2,1),3)
C           ERR6 = ( 96.d0*C+48.d0*CP-12.d0*CPP)-AVEDISP(ITRI(6,2,2),3)
C           ERR7 = (-64.d0*C-32.d0*CP+ 8.d0*CPP)-AVEDISP(ITRI(6,3,0),3)
C           ERR8 = (  0.d0*C+16.d0*CP- 8.d0*CPP)-AVEDISP(ITRI(6,3,1),3)
C           ERR9 = ( 96.d0*C-32.d0*CP-32.d0*CPP)-AVEDISP(ITRI(6,3,2),3)
C           ERR10 =(128.d0*C+64.d0*CP-16.d0*CPP)-AVEDISP(ITRI(6,3,3),3)
C           ERR11 =(-72.d0*C-36.d0*CP+24.d0*CPP)-AVEDISP(ITRI(6,4,0),3)
C           ERR12 =(-96.d0*C+16.d0*CP+40.d0*CPP)-AVEDISP(ITRI(6,4,1),3)
C           ERR13 =(  0.d0*C+ 0.d0*CP+ 0.d0*CPP)-AVEDISP(ITRI(6,4,2),3)
C           ERR14 =(-48.d0*C-24.d0*CP+24.d0*CPP)-AVEDISP(ITRI(6,5,0),3)
C           ERR15 =(-96.d0*C-16.d0*CP+56.d0*CPP)-AVEDISP(ITRI(6,5,1),3)
C           ERR16 =(-16.d0*C- 8.d0*CP+ 8.d0*CPP)-AVEDISP(ITRI(6,6,0),3)
            C   = AVEDISP(ITRI(6,0,0),3) / 4.0d0 
            CP  = -AVEDISP(ITRI(6,1,0),3) / 18.0d0
            CPP = ( 3.0d0*AVEDISP(ITRI(6,3,1),3) - 
     &                 2.0d0*AVEDISP(ITRI(6,1,0),3) ) / 24.0d0
            ERR1 = (  4.d0*C+ 0.d0*CPP+ 0.d0*CP)-AVEDISP(ITRI(6,0,0),3)
            ERR2 = (  0.d0*C+ 0.d0*CPP-18.d0*CP)-AVEDISP(ITRI(6,1,0),3)
            ERR3 = ( 24.d0*C+ 0.d0*CPP+18.d0*CP)-AVEDISP(ITRI(6,1,1),3)
            ERR4 = (-12.d0*C+ 0.d0*CPP-36.d0*CP)-AVEDISP(ITRI(6,2,0),3)
            ERR5 = ( 24.d0*C- 8.d0*CPP- 6.d0*CP)-AVEDISP(ITRI(6,2,1),3)
            ERR6 = ( 48.d0*C+ 0.d0*CPP+54.d0*CP)-AVEDISP(ITRI(6,2,2),3)
            ERR7 = (-32.d0*C+ 0.d0*CPP-36.d0*CP)-AVEDISP(ITRI(6,3,0),3)
            ERR8 = (  0.d0*C+ 8.d0*CPP-12.d0*CP)-AVEDISP(ITRI(6,3,1),3)
            ERR9 = ( 48.d0*C-40.d0*CPP+24.d0*CP)-AVEDISP(ITRI(6,3,2),3)
            ERR10 =( 64.d0*C+ 0.d0*CPP+72.d0*CP)-AVEDISP(ITRI(6,3,3),3)
            ERR11 =(-36.d0*C+ 0.d0*CPP-18.d0*CP)-AVEDISP(ITRI(6,4,0),3)
            ERR12 =(-48.d0*C+32.d0*CPP-12.d0*CP)-AVEDISP(ITRI(6,4,1),3)
            ERR13 =(  0.d0*C+ 0.d0*CPP+ 0.d0*CP)-AVEDISP(ITRI(6,4,2),3)
            ERR14 =(-24.d0*C+ 0.d0*CPP+ 0.d0*CP)-AVEDISP(ITRI(6,5,0),3)
            ERR15 =(-48.d0*C+16.d0*CPP+12.d0*CP)-AVEDISP(ITRI(6,5,1),3)
            ERR16 =( -8.d0*C+ 0.d0*CPP+ 0.d0*CP)-AVEDISP(ITRI(6,6,0),3)
            ERROR = ( DABS(ERR1) + DABS(ERR2) + DABS(ERR3) +
     &                DABS(ERR4) + DABS(ERR5) + DABS(ERR6) +
     &                DABS(ERR7) + DABS(ERR8) + DABS(ERR9) +
     &                DABS(ERR10)+ DABS(ERR11)+ DABS(ERR12)+
     &                DABS(ERR13)+ DABS(ERR14)+ DABS(ERR15)+
     &                DABS(ERR16)     ) / DABS(16.0d0*C)
            IF ( ERROR .GT. THRSDISP) THEN
              WRITE (LUPRI,*) 
     &              'WARNING: LARGE NUMERICAL (?) ERRORS DURING',
     &              ' CALCULATION OF C_ms COEFFICIENTS:'
              WRITE (LUPRI,*) 'TEST #1  GAVE:',ERR1 
              WRITE (LUPRI,*) 'TEST #2  GAVE:',ERR2  
              WRITE (LUPRI,*) 'TEST #3  GAVE:',ERR3   
              WRITE (LUPRI,*) 'TEST #4  GAVE:',ERR4   
              WRITE (LUPRI,*) 'TEST #5  GAVE:',ERR5   
              WRITE (LUPRI,*) 'TEST #6  GAVE:',ERR6   
              WRITE (LUPRI,*) 'TEST #7  GAVE:',ERR7   
              WRITE (LUPRI,*) 'TEST #8  GAVE:',ERR8   
              WRITE (LUPRI,*) 'TEST #9  GAVE:',ERR9   
              WRITE (LUPRI,*) 'TEST #10 GAVE:',ERR10   
              WRITE (LUPRI,*) 'TEST #11 GAVE:',ERR11   
              WRITE (LUPRI,*) 'TEST #12 GAVE:',ERR12   
              WRITE (LUPRI,*) 'TEST #13 GAVE:',ERR13   
              WRITE (LUPRI,*) 'TEST #14 GAVE:',ERR14   
              WRITE (LUPRI,*) 'TEST #15 GAVE:',ERR15   
              WRITE (LUPRI,*) 'TEST #16 GAVE:',ERR16   
              LERROR = .TRUE.
            END IF
            C   = C   / GAMMA0
            CP  = CP  / GAMMA0
            CPP = CPP / GAMMA0
            ABCDE(KC,2)   = C
            ABCDE(KCP,2)  = CP
            ABCDE(KCPP,2) = CPP
            IF (LOCDBG) THEN
               WRITE (LUPRI,*) 
     &              " C_ms, C_ms', Cms''      coefficients:",C,CP,CPP
            END IF
          END IF

        END IF

        IF (NCRDSPE.GE.8) THEN
          D    = AVEDISP(ITRI(8,0,0),1) / 16.0d0
          DP   = (AVEDISP(ITRI(8,2,1),1) -2.d0*AVEDISP(ITRI(8,2,0),1)
     &               + AVEDISP(ITRI(8,1,0),1) ) / 16.0d0
          DPP  = (AVEDISP(ITRI(8,4,2),1) - AVEDISP(ITRI(8,4,1),1)
     &      -AVEDISP(ITRI(8,3,1),1)+4.d0*AVEDISP(ITRI(8,1,0),1))/16.d0
          DPPP = (AVEDISP(ITRI(8,2,0),1) 
     &              - 10.0d0 * AVEDISP(ITRI(8,0,0),1) ) / 18.0d0
          ERR1 =  AVEDISP(ITRI(8,1,0),1) -   64.0d0*D    
          ERR2 =  AVEDISP(ITRI(8,3,0),1) - (256.0d0*D +  54.0d0*DPPP)
          ERR3 =  AVEDISP(ITRI(8,3,1),1) -
     &           ( 576.0d0*D +  48.0d0*DP + 162.0d0*DPPP)
          ERR4 =  AVEDISP(ITRI(8,4,0),1) - (304.0d0*D +  72.0d0*DPPP)
          ERR5 =  AVEDISP(ITRI(8,4,1),1) -
     &           ( 832.0d0*D +  80.0d0*DP + 306.0d0*DPPP)
          ERR6 =  AVEDISP(ITRI(8,5,2),1) -
     &           (1408.0d0*D + 176.0d0*DP +  32.0d0*DPP + 648.0d0*DPPP)

          DPPP = DPPP / GAMMA0
          DPP  = DPP  / GAMMA0
          DP   = DP   / GAMMA0
          D    = D    / GAMMA0

          ERROR = ( DABS(ERR1) + DABS(ERR2) + DABS(ERR3) + DABS(ERR4) +
     &              DABS(ERR5) + DABS(ERR6) ) / DABS(6.0d0*D)
          IF ( ERROR .GT. THRSDISP) THEN
            WRITE (LUPRI,*) 
     &                'WARNING: LARGE NUMERICAL (?) ERRORS DURING',
     &                ' CALCULATION OF D COEFFICIENTS:'
            WRITE (LUPRI,*) 'TEST #1  GAVE:',ERR1 
            WRITE (LUPRI,*) 'TEST #2  GAVE:',ERR2  
            WRITE (LUPRI,*) 'TEST #3  GAVE:',ERR3   
            WRITE (LUPRI,*) 'TEST #4  GAVE:',ERR4   
            WRITE (LUPRI,*) 'TEST #5  GAVE:',ERR5   
            WRITE (LUPRI,*) 'TEST #6  GAVE:',ERR6   
            LERROR = .TRUE.
          END IF
          ABCDE(KD,1)    = D
          ABCDE(KDP,1)   = DP
          ABCDE(KDPP,1)  = DPP
          ABCDE(KDPPP,1) = DPPP
          IF (LOCDBG) THEN
            WRITE (LUPRI,*) " D, D', D'', D'''        coefficients:",
     &                 D, DP, DPP, DPPP
          END IF
        END IF

        IF (NCRDSPE.GE.10) THEN
          E     = AVEDISP(ITRI(10,0,0),1) / 32.0d0
          EP    =(AVEDISP(ITRI(10,2,1),1)-2.d0*AVEDISP(ITRI(10,2,0),1)
     &               + AVEDISP(ITRI(10,1,0),1) ) / 32.0d0
          EPP = (AVEDISP(ITRI(10,4,2),1)-2.0d0*AVEDISP(ITRI(10,4,1),1)
     &   +3.0d0*AVEDISP(ITRI(10,2,1),1)+10.0d0*AVEDISP(ITRI(10,2,0),1)
     &                -29.0d0*AVEDISP(ITRI(10,1,0),1) ) / 32.0d0
          EPPP  = (AVEDISP(ITRI(10,2,0),1) 
     &                 -15.0d0 * AVEDISP(ITRI(10,0,0),1) ) / 36.0d0
          EPPPP=(AVEDISP(ITRI(10,4,1),1)-9.0d0*AVEDISP(ITRI(10,2,1),1)
     &                -14.0d0*AVEDISP(ITRI(10,2,0),1)
     &                +61.0d0*AVEDISP(ITRI(10,1,0),1) ) / 36.0d0

          ERR1 =   32.d0*E - AVEDISP(ITRI(10,0,0),1)
          ERR2 =  160.d0*E - AVEDISP(ITRI(10,1,0),1)
          ERR3 =  480.d0*E+ 36.d0*EPPP - AVEDISP(ITRI(10,2,0),1)
          ERR4 =  800.d0*E+ 32.d0*EP+ 72.d0*EPPP 
     &                                - AVEDISP(ITRI(10,2,1),1)
          ERR5 =  960.d0*E+ 144.d0*EPPP - AVEDISP(ITRI(10,3,0),1)
          ERR6 = 2240.d0*E+128.d0*EP+ 432.d0*EPPP
     &                   - AVEDISP(ITRI(10,3,1),1)
          ERR7 = 1440.d0*E+288.d0*EPPP - AVEDISP(ITRI(10,4,0),1)
          ERR8 = 4160.d0*E+288.d0*EP+1152.d0*EPPP+ 36.d0*EPPPP 
     &                      - AVEDISP(ITRI(10,4,1),1)
          ERR9 = 5760.d0*E+480.d0*EP+32.d0*EPP+1728.d0*EPPP+ 
     &                   72.d0*EPPPP - AVEDISP(ITRI(10,4,2),1)
          ERR10 = 1632.d0*E+ 360.d0*EPPP - AVEDISP(ITRI(10,5,0),1)
          ERR11= 5600.d0*E+416.d0*EP+1800.d0*EPPP+ 108.d0*EPPPP 
     &                 - AVEDISP(ITRI(10,5,1),1)
          ERR12= 9600.d0*E+960.d0*EP+96.d0*EPP+3600.d0*EPPP+
     &                  324.d0*EPPPP - AVEDISP(ITRI(10,5,2),1)
          ERR13=11360.d0*E+1184.d0*EP+128.d0*EPP+4536.d0*EPPP+
     &                  504.d0*EPPPP - AVEDISP(ITRI(10,6,2),1)
          ERR14=14080.d0*E+1632.d0*EP+224.d0*EPP+6048.d0*EPPP+
     &                  792.d0*EPPPP - AVEDISP(ITRI(10,6,3),1)

          E     = E     / GAMMA0
          EP    = EP    / GAMMA0
          EPP   = EPP   / GAMMA0
          EPPP  = EPPP  / GAMMA0
          EPPPP = EPPPP / GAMMA0

          ERROR = ( DABS(ERR1) + DABS(ERR2) + DABS(ERR3) + DABS(ERR4) +
     &              DABS(ERR5) + DABS(ERR6) + DABS(ERR7) + DABS(ERR8) +
     &              DABS(ERR9) + DABS(ERR10)+ DABS(ERR11)+ DABS(ERR12)+
     &              DABS(ERR13)+ DABS(ERR14)  ) / DABS(14.0d0*E)
          IF ( ERROR .GT. THRSDISP ) THEN
            WRITE (LUPRI,*) 
     &                 'WARNING: LARGE NUMERICAL (?) ERRORS DURING',
     &                 ' CALCULATION OF E COEFFICIENTS:'
            WRITE (LUPRI,*) 'TEST #1  GAVE:',ERR1 
            WRITE (LUPRI,*) 'TEST #2  GAVE:',ERR2  
            WRITE (LUPRI,*) 'TEST #3  GAVE:',ERR3   
            WRITE (LUPRI,*) 'TEST #4  GAVE:',ERR4   
            WRITE (LUPRI,*) 'TEST #5  GAVE:',ERR5   
            WRITE (LUPRI,*) 'TEST #6  GAVE:',ERR6   
            WRITE (LUPRI,*) 'TEST #7  GAVE:',ERR7   
            WRITE (LUPRI,*) 'TEST #8  GAVE:',ERR8   
            WRITE (LUPRI,*) 'TEST #9  GAVE:',ERR9   
            WRITE (LUPRI,*) 'TEST #10 GAVE:',ERR10   
            WRITE (LUPRI,*) 'TEST #11 GAVE:',ERR11   
            WRITE (LUPRI,*) 'TEST #12 GAVE:',ERR12   
            WRITE (LUPRI,*) 'TEST #13 GAVE:',ERR13   
            WRITE (LUPRI,*) 'TEST #14 GAVE:',ERR14   
            LERROR = .TRUE.
          END IF
          ABCDE(KE,1)     = E
          ABCDE(KEP,1)    = EP
          ABCDE(KEPP,1)   = EPP
          ABCDE(KEPPP,1)  = EPPP
          ABCDE(KEPPPP,1) = EPPPP
          IF (LOCDBG) THEN
            WRITE (LUPRI,*) " E, E', E'', E''', E'''' coefficients:",
     &                 E, EP, EPP, EPPP, EPPPP
          END IF
        END IF
      END IF

      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCCRDISP                             *
*---------------------------------------------------------------------*
c /* deck CCCDISPRT */
*=====================================================================*
       SUBROUTINE CCCDISPRT (DISPCF, AVEDISP,ABCDE,
     &                       THGCF,  ESHGCF, KERRCF,  DFWMCF,
     &                       AVETHG, AVESHG, AVEKERR, AVEDFWM, LERROR)
*---------------------------------------------------------------------*
*
*   Purpose: print dispersion coefficients second hyperpolarizabilities 
*
*
*   Written by Christof Haettig in March 1998.
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
#include <cccrinf.h>
#include <ccroper.h>

      INTEGER KA, 
     &        KB, KBP, 
     &        KC, KCP, KCPP, 
     &        KD, KDP, KDPP, KDPPP, 
     &        KE, KEP, KEPP, KEPPP, KEPPPP
      PARAMETER (KA=1,
     &           KB=2, KBP=3,
     &           KC=4, KCP=5, KCPP=6,
     &           KD=7, KDP=8, KDPP=9, KDPPP=10,
     &           KE=11,KEP=12,KEPP=13,KEPPP=14,KEPPPP=15)


#if defined (SYS_CRAY)
      REAL DISPCF( (NCRDSPE+3)*(NCRDSPE+2)*(NCRDSPE+1)/6,*)
      REAL AVEDISP((NCRDSPE+3)*(NCRDSPE+2)*(NCRDSPE+1)/6,4)
      REAL ABCDE(15,2)
      REAL THGCF( NCRDSPE+1,NCROPER)
      REAL ESHGCF(NCRDSPE+1,NCROPER)
      REAL KERRCF(NCRDSPE+1,NCROPER)
      REAL DFWMCF(NCRDSPE+1,NCROPER)
      REAL AVETHG( NCRDSPE+1,4)
      REAL AVESHG( NCRDSPE+1,4)
      REAL AVEKERR(NCRDSPE+1,4)
      REAL AVEDFWM(NCRDSPE+1,4)
#else
      DOUBLE PRECISION DISPCF( (NCRDSPE+3)*(NCRDSPE+2)*(NCRDSPE+1)/6,*)
      DOUBLE PRECISION AVEDISP((NCRDSPE+3)*(NCRDSPE+2)*(NCRDSPE+1)/6,4)
      DOUBLE PRECISION ABCDE(15,2)
      DOUBLE PRECISION THGCF( NCRDSPE+1,NCROPER)
      DOUBLE PRECISION ESHGCF(NCRDSPE+1,NCROPER)
      DOUBLE PRECISION KERRCF(NCRDSPE+1,NCROPER)
      DOUBLE PRECISION DFWMCF(NCRDSPE+1,NCROPER)
      DOUBLE PRECISION AVETHG( NCRDSPE+1,4)
      DOUBLE PRECISION AVESHG( NCRDSPE+1,4)
      DOUBLE PRECISION AVEKERR(NCRDSPE+1,4)
      DOUBLE PRECISION AVEDFWM(NCRDSPE+1,4)
#endif

      CHARACTER*5  BLANKS
      CHARACTER*80 STRING
      LOGICAL LERROR
      INTEGER ISYMA, ISYMB, ISYMC, ISYMD
      INTEGER IOPER, MN, LMN

C
      INTEGER ITRI, L, M, N
      ITRI(L,M,N) = L*(L+1)*(L+2)/6 + M*(M+1)/2 + N + 1
C

*---------------------------------------------------------------------*
* print header for hyperpolarizability section
*---------------------------------------------------------------------*
      BLANKS = '     '
      STRING = ' RESULTS FOR DISPERSION COEFFICIENTS'

      IF (CCS) THEN
         CALL AROUND( BLANKS//'FINAL CCS'//STRING(1:45)//BLANKS ) 
      ELSE IF (CC2) THEN
         CALL AROUND( BLANKS//'FINAL CC2'//STRING(1:45)//BLANKS )
      ELSE IF (CCSD) THEN
         CALL AROUND( BLANKS//'FINAL CCSD'//STRING(1:45)//BLANKS )
      ELSE IF (CC3) THEN
         CALL AROUND( BLANKS//'FINAL CC3'//STRING(1:45)//BLANKS )
      ELSE
         CALL QUIT('CCCDISPRT called for an unknown Coupled '//
     &             'Cluster model.')
      END IF

      WRITE(LUPRI,'(/1X,4(A," op.",6X),4(4X,A),/,74("-"))')
     &      'A','B','C','D','L','M','N','D_ABCD'


      DO IOPER = 1, NCROPER
         ISYMA = ISYOPR(IACROP(IOPER))
         ISYMB = ISYOPR(IBCROP(IOPER))
         ISYMC = ISYOPR(ICCROP(IOPER))
         ISYMD = ISYOPR(IDCROP(IOPER))

         IF (MULD2H(ISYMA,ISYMB).EQ.MULD2H(ISYMC,ISYMD)) THEN
            WRITE(LUPRI,'(/1X,4(A8,3X),3I5,1X,G16.8)')
     &        LBLOPR(IACROP(IOPER)), LBLOPR(IBCROP(IOPER)),
     &        LBLOPR(ICCROP(IOPER)), LBLOPR(IDCROP(IOPER)),0,0,0,
     &        -DISPCF(ITRI(0,0,0),IOPER)
         ELSE
           IF (IPRCHYP.GT.0) THEN
            WRITE(LUPRI,'(/1X,4(A8,3X),4A)') LBLOPR(IACROP(IOPER)),
     &        LBLOPR(IBCROP(IOPER)), LBLOPR(ICCROP(IOPER)),
     &        LBLOPR(IDCROP(IOPER)),'  -  ','  -  ', '  -  ', '---'
           END IF
         END IF

         IF (MULD2H(ISYMA,ISYMB).EQ.MULD2H(ISYMC,ISYMD)) THEN
           DO LMN = 2, NCRDSPE, 2
           DO MN = 0, LMN
           DO N = 0, MN
             L = LMN - MN
             M = MN  - N
             WRITE(LUPRI,'(1X,4(8X,3X),3I5,1X,G16.8)')
     &          L,M,N, -DISPCF(ITRI(LMN,MN,N),IOPER)
           END DO
           END DO
           END DO
         END IF

      END DO

      IF (GAMMA_PAR) THEN
         WRITE(LUPRI,'(/1X,A10,34X,3I5,1X,G16.8)')
     &     'gamma_{||} ',0,0,0, -AVEDISP(ITRI(0,0,0),1)
         DO LMN = 2, NCRDSPE, 2
         DO MN = 0, LMN
         DO N = 0, MN
           L = LMN - MN
           M = MN  - N
           WRITE(LUPRI,'(6X,3(8X,5X),3I5,1X,G16.8)')
     &        L,M,N, -AVEDISP(ITRI(LMN,MN,N),1)
         END DO
         END DO
         END DO
      END IF

      IF (GAMMA_PAR .AND. GAMMA_ORT) THEN
         WRITE(LUPRI,'(/1X,A10,34X,3I5,1X,G16.8)') 'gamma_K   ',0,0,0, 
     &      -1.5d0*(AVEDISP(ITRI(0,0,0),1)-AVEDISP(ITRI(0,0,0),2))
         DO LMN = 2, NCRDSPE, 2
         DO MN = 0, LMN
         DO N = 0, MN
           L = LMN - MN
           M = MN  - N
           WRITE(LUPRI,'(1X,10X,34X,3I5,1X,G16.8)') L,M,N, 
     &     -1.5d0*(AVEDISP(ITRI(LMN,MN,N),1)-AVEDISP(ITRI(LMN,MN,N),2))
         END DO
         END DO
         END DO
      END IF

      IF (GAMMA_ORT) THEN
         WRITE(LUPRI,'(/1X,A11,33X,3I5,1X,G16.8)')
     &     'gamma_{_|_}',0,0,0, -AVEDISP(ITRI(0,0,0),2)
         DO LMN = 2, NCRDSPE, 2
         DO MN = 0, LMN
         DO N = 0, MN
           L = LMN - MN
           M = MN  - N
           WRITE(LUPRI,'(6X,3(8X,5X),3I5,1X,G16.8)')
     &        L,N,M, -AVEDISP(ITRI(LMN,MN,N),2)
         END DO
         END DO
         END DO

         WRITE(LUPRI,'(/1X,A10,34X,3I5,1X,G16.8)')
     &     'gamma_{ms} ',0,0,0, -AVEDISP(ITRI(0,0,0),3)
         DO LMN = 2, NCRDSPE, 2
         DO MN = 0, LMN
         DO N = 0, MN
           L = LMN - MN
           M = MN  - N
           WRITE(LUPRI,'(6X,3(8X,5X),3I5,1X,G16.8)')
     &        L,N,M, -AVEDISP(ITRI(LMN,MN,N),3)
         END DO
         END DO
         END DO
      END IF

      WRITE(LUPRI,'(74("-"),//)')

      WRITE(LUPRI,'(////1X,4(A," op.",6X),2(4X,A),/,112("-"))')
     & 'A','B','C','D','N',
     & 'D^THG           D^ESHG          D^DFWM          D^KERR'

      DO IOPER = 1, NCROPER
         ISYMA = ISYOPR(IACROP(IOPER))
         ISYMB = ISYOPR(IBCROP(IOPER))
         ISYMC = ISYOPR(ICCROP(IOPER))
         ISYMD = ISYOPR(IDCROP(IOPER))

         IF (MULD2H(ISYMA,ISYMB).EQ.MULD2H(ISYMC,ISYMD)) THEN
            WRITE(LUPRI,'(/1X,4(A8,3X),I5,1X,4G16.8)')
     &        LBLOPR(IACROP(IOPER)),
     &        LBLOPR(IBCROP(IOPER)),
     &        LBLOPR(ICCROP(IOPER)),
     &        LBLOPR(IDCROP(IOPER)),0,
     &        -THGCF(1,IOPER),  -ESHGCF(1,IOPER),
     &        -DFWMCF(1,IOPER), -KERRCF(1,IOPER)
         ELSE
           IF (IPRCHYP.GT.0) THEN
            WRITE(LUPRI,'(/1X,4(A8,3X),A,4(A,8X))')
     &        LBLOPR(IACROP(IOPER)),
     &        LBLOPR(IBCROP(IOPER)),
     &        LBLOPR(ICCROP(IOPER)),
     &        LBLOPR(IDCROP(IOPER)),'  -  ',
     &        '---','---','---','---'
           END IF
         END IF

         IF (MULD2H(ISYMA,ISYMB).EQ.MULD2H(ISYMC,ISYMD)) THEN
           DO N = 2, NCRDSPE, 2
             WRITE(LUPRI,'(1X,4(8X,3X),I5,1X,4G16.8)') N,
     &        -THGCF(N+1,IOPER), -ESHGCF(N+1,IOPER),
     &        -DFWMCF(N+1,IOPER),-KERRCF(N+1,IOPER)
           END DO
         END IF

      END DO

      IF (GAMMA_PAR) THEN
        WRITE(LUPRI,'(/1X,A10,34X,I5,1X,4G16.8)') 'gamma_{||} ',0,
     &    -AVETHG(1,1),-AVESHG(1,1),-AVEDFWM(1,1),-AVEKERR(1,1)
        DO N = 2, NCRDSPE, 2
          WRITE(LUPRI,'(1X,44X,I5,1X,4G16.8)') N,
     &    -AVETHG(N+1,1),-AVESHG(N+1,1),-AVEDFWM(N+1,1),-AVEKERR(N+1,1)
        END DO
      END IF
      IF (GAMMA_PAR .AND. GAMMA_ORT) THEN
        WRITE(LUPRI,'(/1X,A10,34X,I5,1X,4G16.8)') 'gamma^K    ',0,
     &    -1.5d0*(AVETHG(1,1) -AVETHG(1,2)),
     &    -1.5d0*(AVESHG(1,1) -AVESHG(1,2)),
     &    -1.5d0*(AVEDFWM(1,1)-AVEDFWM(1,2)),
     &    -1.5d0*(AVEKERR(1,1)-AVEKERR(1,2))
        DO N = 2, NCRDSPE, 2
          WRITE(LUPRI,'(1X,44X,I5,1X,4G16.8)') N,
     &    -1.5d0*(AVETHG(N+1,1) -AVETHG(N+1,2)),
     &    -1.5d0*(AVESHG(N+1,1) -AVESHG(N+1,2)),
     &    -1.5d0*(AVEDFWM(N+1,1)-AVEDFWM(N+1,2)),
     &    -1.5d0*(AVEKERR(N+1,1)-AVEKERR(N+1,2))
        END DO
      END IF
      IF (GAMMA_ORT) THEN
        WRITE(LUPRI,'(/1X,A11,33X,I5,1X,4G16.8)') 'gamma_{_|_}',0,
     &    -AVETHG(1,2),-AVESHG(1,2),-AVEDFWM(1,2),-AVEKERR(1,2)
        DO N = 2, NCRDSPE, 2
          WRITE(LUPRI,'(1X,44X,I5,1X,4G16.8)') N,
     &    -AVETHG(N+1,2),-AVESHG(N+1,2),-AVEDFWM(N+1,2),-AVEKERR(N+1,2)
        END DO
        WRITE(LUPRI,'(/1X,A10,34X,I5,1X,4G16.7)') 'gamma_{ms} ',0,
     &    -AVETHG(1,3),-AVESHG(1,3),-AVEDFWM(1,3),-AVEKERR(1,3)
        DO N = 2, NCRDSPE, 2
          WRITE(LUPRI,'(1X,44X,I5,1X,4G16.8)') N,
     &    -AVETHG(N+1,3),-AVESHG(N+1,3),-AVEDFWM(N+1,3),-AVEKERR(N+1,3)
        END DO
      END IF

      WRITE(LUPRI,'(112("-"),//)')

      IF (GAMMA_PAR) THEN
        WRITE(LUPRI,'(////1X,A,/,78("-"))')
     &   'Process independent dispersion coefficients for gamma_{||}:'
        WRITE(LUPRI,'(1X,A,G16.8)')         'gamma_0',-AVEDISP(1,1)
        IF (NCRDSPE.GE.2) THEN
          WRITE(LUPRI,'(1X,A,2X,G16.8,5X)') "A    ",ABCDE(KA,1)
        END IF
        IF (NCRDSPE.GE.4) THEN
          WRITE(LUPRI,'(1X,2(A,2X,G16.8,5X))')
     &       "B    ",ABCDE(KB,1),   "B'   ",ABCDE(KBP,1)
        END IF
        IF (NCRDSPE.GE.6) THEN
          WRITE(LUPRI,'(1X,3(A,2X,G16.8,5X))')
     &       "C    ",ABCDE(KC,1),   "C'   ",ABCDE(KCP,1), 
     &       "C''  ",ABCDE(KCPP,1)
        END IF
        IF (NCRDSPE.GE.8) THEN
          WRITE(LUPRI,'(1X,4(A,2X,G16.8,5X))')
     &       "D    ",ABCDE(KD,1),   "D'   ",ABCDE(KDP,1),
     &       "D''  ",ABCDE(KDPP,1), "D''' ",ABCDE(KDPPP,1)
        END IF
        IF (NCRDSPE.GE.10) THEN
          WRITE(LUPRI,'(1X,5(A,2X,G16.8,5X))')
     &       "E    ",ABCDE(KE,1),   "E'   ",ABCDE(KEP,1),
     &       "E''  ",ABCDE(KEPP,1), "E''' ",ABCDE(KEPPP,1),
     &       "E''''",ABCDE(KEPPPP,1)
        END IF
        IF (NCRDSPE.GE.2 .AND. LERROR) THEN
          WRITE(LUPRI,'(/1X,A,/)') 'large numerical (?) errors ',
     &     'encountered ... check output for warnings.'
        END IF
        WRITE(LUPRI,'(78("-"))')
      END IF

      IF (GAMMA_ORT) THEN
        WRITE(LUPRI,'(///1X,A,/,78("-"))')
     &   'Process independent dispersion coefficients for gamma_{ms}:'
        IF (NCRDSPE.GE.2)
     &    WRITE(LUPRI,'(1X,A,2X,G16.8)') "A    ",ABCDE(KA,2)
        IF (NCRDSPE.GE.4) THEN
          WRITE(LUPRI,'(1X,A,2X,G16.8,5X,A,2X,G16.8)')
     &      "B    ",ABCDE(KB,2), "B'   ",ABCDE(KBP,2)
        END IF
        IF (NCRDSPE.GE.6) THEN
          WRITE(LUPRI,'(1X,3(A,2X,G16.8,5X))') "C    ",ABCDE(KC,2),
     &      "C'   ",ABCDE(KCP,2),"C''  ",ABCDE(KCPP,2)
        END IF
        IF (NCRDSPE.GE.2 .AND. LERROR) THEN
          WRITE(LUPRI,'(/1X,A,/)') 'large numerical (?) errors ',
     &     'encountered ... check output for warnings.'
        END IF
        WRITE(LUPRI,'(78("-"))')
      END IF

      WRITE(LUPRI,'(///)')

      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCCDISPRT                            *
*---------------------------------------------------------------------*
