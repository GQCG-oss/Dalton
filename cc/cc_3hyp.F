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
c /* deck CC_3HYP */
*=====================================================================*
       SUBROUTINE CC_3HYP(WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: direct calculation of frequency-dependent third
*             hyperpolarizabilities for the Couple Cluster models
*
*                        CCS, CC2, CCSD
*
*     Written by Christof Haettig, april 1997.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdinp.h>
#include <cc4rinf.h>
#include <ccr1rsp.h>

* local parameters:
      CHARACTER*(19) MSGDBG
      PARAMETER (MSGDBG = '[debug] CC_3HYP> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE. )

      CHARACTER*10 MODEL
      LOGICAL LADD
      INTEGER LWORK
      INTEGER KHYPPOL, NBHYPPOL
      INTEGER K0HTRAN, K0HDOTS, K0HCONS, N0HTRAN, MX0HTRAN, MX0HDOTS
      INTEGER K1HTRAN, K1HDOTS, K1HCONS, N1HTRAN, MX1HTRAN, MX1HDOTS
      INTEGER K0GTRAN, K0GDOTS, K0GCONS, N0GTRAN, MX0GTRAN, MX0GDOTS
      INTEGER K1GTRAN, K1GDOTS, K1GCONS, N1GTRAN, MX1GTRAN, MX1GDOTS
      INTEGER K2GTRAN, K2GDOTS, K2GCONS, N2GTRAN, MX2GTRAN, MX2GDOTS
      INTEGER K1FTRAN, K1FDOTS, K1FCONS, N1FTRAN, MX1FTRAN, MX1FDOTS
      INTEGER K2FTRAN, K2FDOTS, K2FCONS, N2FTRAN, MX2FTRAN, MX2FDOTS
      INTEGER K0FATRAN,K0FADOTS,K0FACONS,N0FATRAN,MX0FATRAN,MX0FADOTS
      INTEGER K1FATRAN,K1FADOTS,K1FACONS,N1FATRAN,MX1FATRAN,MX1FADOTS
      INTEGER K2FATRAN,K2FADOTS,K2FACONS,N2FATRAN,MX2FATRAN,MX2FADOTS
      INTEGER K2EATRAN,K2EADOTS,K2EACONS,N2EATRAN,MX2EATRAN,MX2EADOTS
      INTEGER KXTRAN,  KXDOTS,  KXCONS,  NXTRAN,  MXXTRAN,  MXXDOTS
      INTEGER MXTRAN3, MXTRAN4, MXVEC1, MXVEC2
      INTEGER KEND0, LEND0, IOPT

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
#else
      DOUBLE PRECISION WORK(LWORK)
#endif

*---------------------------------------------------------------------*
* print header for hyperpolarizability section
*---------------------------------------------------------------------*
      WRITE (LUPRI,'(7(/1X,2A),/)')
     & '************************************',
     &                               '*******************************',
     & '*                                   ',
     &                               '                              *',
     & '*<<<<<<<<< OUTPUT FROM COUPLED CLUSTE',
     &                                'R QUARTIC RESPONSE  >>>>>>>>>*',
     & '*                                   ',
     &                               '                              *',
     & '*<<<<<<<<     CALCULATION OF CC HYPER',
     &                                'POLARIZABILITIES    >>>>>>>>>*',
     & '*                                   ',
     &                               '                              *',
     & '************************************',
     &                               '*******************************' 

*---------------------------------------------------------------------*
      IF (.NOT. (CCS .OR. CC2 .OR. CCSD) ) THEN
         CALL QUIT('CC_3HYP called for unknown Coupled Cluster.')
      END IF

* print some debug/info output
      IF (IPR4HYP .GT. 10) WRITE(LUPRI,*) 'CC_3HYP Workspace:',LWORK
  
*---------------------------------------------------------------------*
* allocate & initialize work space for hyperpolarizabilities
*---------------------------------------------------------------------*
      NBHYPPOL = N4ROPER * N4RFREQ

      KHYPPOL = 1
      KEND0   = KHYPPOL + 2*NBHYPPOL


      MXTRAN3 = 2 * NLRTLBL * NLRTLBL * NLRTLBL
      MXVEC2  = (NLRTLBL+1) * NLRTLBL / 2
      MXTRAN4 = 2 * NLRTLBL * NLRTLBL * NLRTLBL * NLRTLBL
      MXVEC1  = NLRTLBL
      MXTRAN3 = MIN(MXTRAN3,2*60*NBHYPPOL) ! 60 = # perm. for F^A{O}
      MXTRAN4 = MIN(MXTRAN4,2*30*NBHYPPOL) ! 30 = # perm. for F^AB{O}
      MXVEC2  = MIN(MXVEC2,2*60*NBHYPPOL)  ! 60 = # perm. for F^A{O}
      MXVEC1  = MIN(MXVEC1,2*30*NBHYPPOL)  ! 30 = # perm. for F^AB{O}

      MX0HTRAN  = 5 * MXTRAN3
      MX1HTRAN  = 5 * MXTRAN4
      MX0GTRAN  = 4 * MXTRAN3
      MX1GTRAN  = 4 * MXTRAN3
      MX2GTRAN  = 4 * MXTRAN4
      MX1FTRAN  = 3 * MXTRAN3
      MX2FTRAN  = 3 * MXTRAN3
      MX0FATRAN = 5 * MXTRAN3
      MX1FATRAN = 5 * MXTRAN3
      MX2FATRAN = 5 * MXTRAN4
      MX2EATRAN = 3 * MXTRAN3
      MXXTRAN   = 1 * MXTRAN3

      MX0HDOTS  = MXTRAN3 * MXVEC2
      MX1HDOTS  = MXTRAN4 * MXVEC1
      MX0GDOTS  = MXTRAN3 * MXVEC2
      MX1GDOTS  = MXTRAN3 * MXVEC2
      MX2GDOTS  = MXTRAN4 * MXVEC1
      MX1FDOTS  = MXTRAN3 * MXVEC2
      MX2FDOTS  = MXTRAN3 * MXVEC2
      MX0FADOTS = MXTRAN3 * MXVEC2
      MX1FADOTS = MXTRAN3 * MXVEC2
      MX2FADOTS = MXTRAN4 * MXVEC1
      MX2EADOTS = MXTRAN3 * MXVEC2
      MXXDOTS   = MXTRAN3 * MXVEC2

      K0HTRAN  = KEND0
      K0HDOTS  = K0HTRAN  + MX0HTRAN
      K0HCONS  = K0HDOTS  + MX0HDOTS
      KEND0    = K0HCONS  + MX0HDOTS

      K1HTRAN  = KEND0
      K1HDOTS  = K1HTRAN  + MX1HTRAN
      K1HCONS  = K1HDOTS  + MX1HDOTS
      KEND0    = K1HCONS  + MX1HDOTS

      K0GTRAN  = KEND0
      K0GDOTS  = K0GTRAN  + MX0GTRAN
      K0GCONS  = K0GDOTS  + MX0GDOTS
      KEND0    = K0GCONS  + MX0GDOTS

      K1GTRAN  = KEND0
      K1GDOTS  = K1GTRAN  + MX1GTRAN
      K1GCONS  = K1GDOTS  + MX1GDOTS
      KEND0    = K1GCONS  + MX1GDOTS

      K2GTRAN  = KEND0
      K2GDOTS  = K2GTRAN  + MX2GTRAN
      K2GCONS  = K2GDOTS  + MX2GDOTS
      KEND0    = K2GCONS  + MX2GDOTS

      K1FTRAN  = KEND0
      K1FDOTS  = K1FTRAN  + MX1FTRAN
      K1FCONS  = K1FDOTS  + MX1FDOTS
      KEND0    = K1FCONS  + MX1FDOTS

      K2FTRAN  = KEND0
      K2FDOTS  = K2FTRAN  + MX2FTRAN
      K2FCONS  = K2FDOTS  + MX2FDOTS
      KEND0    = K2FCONS  + MX2FDOTS

      K0FATRAN = KEND0
      K0FADOTS = K0FATRAN + MX0FATRAN
      K0FACONS = K0FADOTS + MX0FADOTS
      KEND0    = K0FACONS + MX0FADOTS

      K1FATRAN = KEND0
      K1FADOTS = K1FATRAN + MX1FATRAN
      K1FACONS = K1FADOTS + MX1FADOTS
      KEND0    = K1FACONS + MX1FADOTS

      K2FATRAN = KEND0
      K2FADOTS = K2FATRAN + MX2FATRAN
      K2FACONS = K2FADOTS + MX2FADOTS
      KEND0    = K2FACONS + MX2FADOTS

      K2EATRAN = KEND0
      K2EADOTS = K2EATRAN + MX2EATRAN
      K2EACONS = K2EADOTS + MX2EADOTS
      KEND0    = K2EACONS + MX2EADOTS

      KXTRAN   = KEND0
      KXDOTS   = KXTRAN   + MXXTRAN
      KXCONS   = KXDOTS   + MXXDOTS
      KEND0    = KXCONS   + MXXDOTS

      LEND0    = LWORK - KEND0

      IF (LEND0.LT.0) THEN
        WRITE (LUPRI,*) 'KEND0,LEND0:',KEND0,LEND0
        CALL QUIT('Insufficient memory in CC_3HYP.')
      END IF

* initialize hyperpolarizabilities and the lists of contributions:
      CALL DZERO(WORK(KHYPPOL),2*NBHYPPOL)
      CALL DZERO(WORK(K0HTRAN),MX0HTRAN+2*MX0HDOTS)
      CALL DZERO(WORK(K1HTRAN),MX1HTRAN+2*MX1HDOTS)
      CALL DZERO(WORK(K0GTRAN),MX0GTRAN+2*MX0GDOTS)
      CALL DZERO(WORK(K1GTRAN),MX1GTRAN+2*MX1GDOTS)
      CALL DZERO(WORK(K2GTRAN),MX2GTRAN+2*MX2GDOTS)
      CALL DZERO(WORK(K1FTRAN),MX1FTRAN+2*MX1FDOTS)
      CALL DZERO(WORK(K2FTRAN),MX2FTRAN+2*MX2FDOTS)
      CALL DZERO(WORK(K0FATRAN),MX0FATRAN+2*MX0FADOTS)
      CALL DZERO(WORK(K1FATRAN),MX1FATRAN+2*MX1FADOTS)
      CALL DZERO(WORK(K2FATRAN),MX2FATRAN+2*MX2FADOTS)
      CALL DZERO(WORK(K2EATRAN),MX2EATRAN+2*MX2EADOTS)
      CALL DZERO(WORK(KXTRAN),MXXTRAN+2*MXXDOTS)

*---------------------------------------------------------------------*
* set up lists for H, G and F transformations and ETA vectors:
*---------------------------------------------------------------------*

      LADD = .FALSE.

      CALL CC4R_SETUP(MXTRAN3, MXVEC2, MXTRAN4, MXVEC1,
     &     WORK(K0HTRAN), WORK(K0HDOTS), WORK(K0HCONS), N0HTRAN,
     &     WORK(K1HTRAN), WORK(K1HDOTS), WORK(K1HCONS), N1HTRAN,
     &     WORK(K0GTRAN), WORK(K0GDOTS), WORK(K0GCONS), N0GTRAN,
     &     WORK(K1GTRAN), WORK(K1GDOTS), WORK(K1GCONS), N1GTRAN,
     &     WORK(K2GTRAN), WORK(K2GDOTS), WORK(K2GCONS), N2GTRAN,
     &     WORK(K1FTRAN), WORK(K1FDOTS), WORK(K1FCONS), N1FTRAN,
     &     WORK(K2FTRAN), WORK(K2FDOTS), WORK(K2FCONS), N2FTRAN,
     &     WORK(K0FATRAN),WORK(K0FADOTS),WORK(K0FACONS),N0FATRAN,
     &     WORK(K1FATRAN),WORK(K1FADOTS),WORK(K1FACONS),N1FATRAN,
     &     WORK(K2FATRAN),WORK(K2FADOTS),WORK(K2FACONS),N2FATRAN,
     &     WORK(K2EATRAN),WORK(K2EADOTS),WORK(K2EACONS),N2EATRAN,
     &     WORK(KXTRAN),  WORK(KXDOTS),  WORK(KXCONS),  NXTRAN,
     &     WORK(KHYPPOL), NBHYPPOL, LADD )

*---------------------------------------------------------------------*
* calculate H matrix contributions:
*---------------------------------------------------------------------*
      IF (.NOT.L_USE_CHI3) THEN
        IOPT = 5
        CALL CC_HMAT('L0','R1','R1','R1','R2',N0HTRAN, MXVEC2,
     &                WORK(K0HTRAN),WORK(K0HDOTS),WORK(K0HCONS),
     &                WORK(KEND0), LEND0, IOPT )
      END IF

      IOPT = 5
      CALL CC_HMAT('L1 ','R1 ','R1 ','R1 ','R1 ',N1HTRAN, MXVEC1,
     &              WORK(K1HTRAN),WORK(K1HDOTS),WORK(K1HCONS),
     &              WORK(KEND0), LEND0, IOPT )

*---------------------------------------------------------------------*
* calculate G matrix contributions:
*---------------------------------------------------------------------*
      IOPT = 5
      CALL CC_GMATRIX('L0 ','R2 ','R1 ','R2 ',N0GTRAN, MXVEC2,
     &              WORK(K0GTRAN),WORK(K0GDOTS),WORK(K0GCONS),
     &              WORK(KEND0), LEND0, IOPT )

      IF (L_USE_CHI3) CALL DSCAL(MX0GDOTS,-1.0d0,WORK(K0GCONS),1)

      IF (.NOT.L_USE_CHI3) THEN
        IOPT = 5
        CALL CC_GMATRIX('L1 ','R1 ','R1 ','R2 ',N1GTRAN, MXVEC2,
     &                WORK(K1GTRAN),WORK(K1GDOTS),WORK(K1GCONS),
     &                WORK(KEND0), LEND0, IOPT )
      END IF

      IOPT = 5
      CALL CC_GMATRIX('L2 ','R1 ','R1 ','R1 ',N2GTRAN, MXVEC1,
     &              WORK(K2GTRAN),WORK(K2GDOTS),WORK(K2GCONS),
     &              WORK(KEND0), LEND0, IOPT )

      
*---------------------------------------------------------------------*
* calculate F matrix contributions:
*---------------------------------------------------------------------*
      IOPT = 5
      CALL CC_FMATRIX(WORK(K1FTRAN),N1FTRAN,'L1 ','R2 ',IOPT,'R2 ',
     &                WORK(K1FDOTS),WORK(K1FCONS),MXVEC2,
     &                WORK(KEND0), LEND0)

      IF (L_USE_CHI3) CALL DSCAL(MX1FDOTS,-1.0d0,WORK(K1FCONS),1)

      IF (.NOT.L_USE_CHI3) THEN
        IOPT = 5
        CALL CC_FMATRIX(WORK(K2FTRAN),N2FTRAN,'L2 ','R1 ',IOPT,'R2 ',
     &                  WORK(K2FDOTS),WORK(K2FCONS),MXVEC2,
     &                  WORK(KEND0), LEND0)
      END IF

*---------------------------------------------------------------------*
* calculate F{O} matrix contributions:
*---------------------------------------------------------------------*
      CALL CCQR_FADRV('L0 ','o1 ','R2 ','R2 ',N0FATRAN, MXVEC2,
     &                 WORK(K0FATRAN),WORK(K0FADOTS),
     &                 WORK(K0FACONS),
     &                 WORK(KEND0), LEND0, 'DOTP' )

      IF (L_USE_CHI3) CALL DSCAL(MX0FADOTS,-1.0d0,WORK(K0FACONS),1)

      IF (.NOT.L_USE_CHI3) THEN
        CALL CCQR_FADRV('L1 ','o1 ','R1 ','R2 ',N1FATRAN, MXVEC2,
     &                   WORK(K1FATRAN),WORK(K1FADOTS),
     &                   WORK(K1FACONS),
     &                   WORK(KEND0), LEND0, 'DOTP' )
      END IF

      CALL CCQR_FADRV('L2 ','o1 ','R1 ','R1 ',N2FATRAN, MXVEC1,
     &                 WORK(K2FATRAN),WORK(K2FADOTS),
     &                 WORK(K2FACONS),
     &                 WORK(KEND0), LEND0, 'DOTP' )

*---------------------------------------------------------------------*
* calculate ETA{O} vector contributions:
*---------------------------------------------------------------------*
      IF (.NOT.L_USE_CHI3) THEN
        CALL CCQR_EADRV('L2 ','o1 ','R2 ',N2EATRAN, MXVEC2,
     &                   WORK(K2EATRAN),WORK(K2EADOTS),WORK(K2EACONS),
     &                   WORK(KEND0), LEND0, 'DOTP' )
      END IF

*---------------------------------------------------------------------*
* chi vector contributions:
*---------------------------------------------------------------------*
      IF (L_USE_CHI3) THEN
        CALL CC_DOTDRV('X3 ','R2 ',NXTRAN,MXVEC2,
     &                 WORK(KXTRAN), WORK(KXDOTS), WORK(KXCONS),
     &                 WORK(KEND0), LEND0 )
C       CALL CC_DOTDRV('L3 ','O2 ',NXTRAN,MXVEC2,
C    &                 WORK(KXTRAN), WORK(KXDOTS), WORK(KXCONS),
C    &                 WORK(KEND0), LEND0 )
      END IF

*---------------------------------------------------------------------*
* collect contributions and add them to hyperpolarizabilities
*---------------------------------------------------------------------*

      LADD = .TRUE.

      CALL CC4R_SETUP(MXTRAN3, MXVEC2, MXTRAN4, MXVEC1,
     &     WORK(K0HTRAN), WORK(K0HDOTS), WORK(K0HCONS), N0HTRAN,
     &     WORK(K1HTRAN), WORK(K1HDOTS), WORK(K1HCONS), N1HTRAN,
     &     WORK(K0GTRAN), WORK(K0GDOTS), WORK(K0GCONS), N0GTRAN,
     &     WORK(K1GTRAN), WORK(K1GDOTS), WORK(K1GCONS), N1GTRAN,
     &     WORK(K2GTRAN), WORK(K2GDOTS), WORK(K2GCONS), N2GTRAN,
     &     WORK(K1FTRAN), WORK(K1FDOTS), WORK(K1FCONS), N1FTRAN,
     &     WORK(K2FTRAN), WORK(K2FDOTS), WORK(K2FCONS), N2FTRAN,
     &     WORK(K0FATRAN),WORK(K0FADOTS),WORK(K0FACONS),N0FATRAN,
     &     WORK(K1FATRAN),WORK(K1FADOTS),WORK(K1FACONS),N1FATRAN,
     &     WORK(K2FATRAN),WORK(K2FADOTS),WORK(K2FACONS),N2FATRAN,
     &     WORK(K2EATRAN),WORK(K2EADOTS),WORK(K2EACONS),N2EATRAN,
     &     WORK(KXTRAN),  WORK(KXDOTS),  WORK(KXCONS),  NXTRAN,
     &     WORK(KHYPPOL), NBHYPPOL, LADD )


*---------------------------------------------------------------------*
* print output & return:
*---------------------------------------------------------------------*

      CALL CC4HYPPRT(WORK(KHYPPOL))

      RETURN
      END

*=====================================================================*
*              END OF SUBROUTINE CC_3HYP                              *
*=====================================================================*
c /* deck CC4HYPPRT */
*=====================================================================*
       SUBROUTINE CC4HYPPRT(HYPERPOL)
*---------------------------------------------------------------------*
*
*    Purpose: print third hyperpolarizabilities 
*
*    Written by Christof Haettig in april 1997.
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
#include <cc4rinf.h>
#include <ccroper.h>


      CHARACTER*10 BLANKS
      CHARACTER*80 STRING
      INTEGER ISYMA, ISYMB, ISYMC, ISYMD, ISYME
      INTEGER IFREQ, IOPER, ISYOLD, ISYTOT


#if defined (SYS_CRAY)
      REAL HYPERPOL(N4RFREQ,N4ROPER,2)
      REAL HALF
#else
      DOUBLE PRECISION HYPERPOL(N4RFREQ,N4ROPER,2)
      DOUBLE PRECISION HALF
#endif
      PARAMETER (HALF = 0.5D0)

*---------------------------------------------------------------------*
* print header for hyperpolarizability section
*---------------------------------------------------------------------*
      BLANKS = '          '
      STRING = ' RESULTS FOR THE THIRD HYPERPOLARIZABILITIES '

      IF (CCS) THEN
         CALL AROUND( BLANKS//'FINAL CCS'//STRING(1:45)//BLANKS ) 
      ELSE IF (CC2) THEN
         CALL AROUND( BLANKS//'FINAL CC2'//STRING(1:45)//BLANKS )
      ELSE IF (CCSD) THEN
         CALL AROUND( BLANKS//'FINAL CCSD'//STRING(1:45)//BLANKS )
      ELSE
         CALL QUIT('CC4HYPPRT called for an unknown Coupled '//
     &             'Cluster model.')
      END IF

      IF (IPR4HYP.GT.5) THEN
       WRITE(LUPRI,'(/1X,5(1X,A," operator",5X),4X,A,8X,A,/,121("-"))')
     &     'A','B','C','D','E','delta','(asy. Resp.)'
      ELSE
       WRITE(LUPRI,'(/1X,5(1X,A," operator",5X),4X,A,/,86("-"))')
     &     'A','B','C','D','E','delta'
      END IF

      ISYOLD = 1

      DO IOPER = 1, N4ROPER
         ISYTOT = 1

         ISYMA  = ISYOPR(IA4ROP(IOPER))
         ISYTOT = MULD2H(ISYTOT,ISYMA)

         ISYMB  = ISYOPR(IB4ROP(IOPER))
         ISYTOT = MULD2H(ISYTOT,ISYMB)

         ISYMC  = ISYOPR(IC4ROP(IOPER))
         ISYTOT = MULD2H(ISYTOT,ISYMC)

         ISYMD  = ISYOPR(ID4ROP(IOPER))
         ISYTOT = MULD2H(ISYTOT,ISYMD)

         ISYME  = ISYOPR(IE4ROP(IOPER))
         ISYTOT = MULD2H(ISYTOT,ISYME)


         IFREQ = 1
         IF (ISYTOT.EQ.1) THEN
          IF (IPR4HYP.GT.5) THEN
            WRITE(LUPRI,'(/1X,5(A7,F7.4,2X),G16.8," (",G10.3,")")')
     &        LBLOPR(IA4ROP(IOPER))(1:7),A4RFR(IFREQ),
     &        LBLOPR(IB4ROP(IOPER))(1:7),B4RFR(IFREQ),
     &        LBLOPR(IC4ROP(IOPER))(1:7),C4RFR(IFREQ),
     &        LBLOPR(ID4ROP(IOPER))(1:7),D4RFR(IFREQ),
     &        LBLOPR(IE4ROP(IOPER))(1:7),E4RFR(IFREQ),
     &        HALF*(HYPERPOL(IFREQ,IOPER,1)+HYPERPOL(IFREQ,IOPER,2)),
     &        HALF*(HYPERPOL(IFREQ,IOPER,1)-HYPERPOL(IFREQ,IOPER,2))
          ELSE
            WRITE(LUPRI,'(/1X,5(A7,F7.4,2X),G16.8)')
     &        LBLOPR(IA4ROP(IOPER))(1:7),A4RFR(IFREQ),
     &        LBLOPR(IB4ROP(IOPER))(1:7),B4RFR(IFREQ),
     &        LBLOPR(IC4ROP(IOPER))(1:7),C4RFR(IFREQ),
     &        LBLOPR(ID4ROP(IOPER))(1:7),D4RFR(IFREQ),
     &        LBLOPR(IE4ROP(IOPER))(1:7),E4RFR(IFREQ),
     &        HALF*(HYPERPOL(IFREQ,IOPER,1)+HYPERPOL(IFREQ,IOPER,2))
          END IF
          ISYOLD = 1
         ELSE
          IF (IPR4HYP.GT.5) THEN
            IF (ISYOLD.EQ.1) WRITE(LUPRI,*)
            WRITE(LUPRI,'(1X,5(A7,A7,2X),6X,A,7X," (",5X,A,6X,")")')
     &        LBLOPR(IA4ROP(IOPER))(1:7),'   -.-  ',
     &        LBLOPR(IB4ROP(IOPER))(1:7),'   -.-  ',
     &        LBLOPR(IC4ROP(IOPER))(1:7),'   -.-  ', 
     &        LBLOPR(ID4ROP(IOPER))(1:7),'   -.-  ', 
     &        LBLOPR(IE4ROP(IOPER))(1:7),'   -.-  ', 
     &        '---',
     &        '---'
          ELSE IF (IPR4HYP.GT.0) THEN
            IF (ISYOLD.EQ.1) WRITE(LUPRI,*)
            WRITE(LUPRI,'(1X,5(A7,A7,2X),6X,A,7X)')
     &        LBLOPR(IA4ROP(IOPER))(1:7),'   -.-  ',
     &        LBLOPR(IB4ROP(IOPER))(1:7),'   -.-  ',
     &        LBLOPR(IC4ROP(IOPER))(1:7),'   -.-  ', 
     &        LBLOPR(ID4ROP(IOPER))(1:7),'   -.-  ', 
     &        LBLOPR(IE4ROP(IOPER))(1:7),'   -.-  ', 
     &        '---'
          END IF
          ISYOLD = 0
         END IF

         DO IFREQ = 2, N4RFREQ
          IF (ISYTOT.EQ.1) THEN
           IF (IPR4HYP.GT.5) THEN
            WRITE(LUPRI,'(1X,5(7X,F7.4,2X),G16.8," (",G10.3,")")')
     &        A4RFR(IFREQ), B4RFR(IFREQ), C4RFR(IFREQ), D4RFR(IFREQ),
     &        E4RFR(IFREQ),
     &        HALF*(HYPERPOL(IFREQ,IOPER,1)+HYPERPOL(IFREQ,IOPER,2)),
     &        HALF*(HYPERPOL(IFREQ,IOPER,1)-HYPERPOL(IFREQ,IOPER,2))
           ELSE
            WRITE(LUPRI,'(1X,5(7X,F7.4,2X),G16.8)')
     &        A4RFR(IFREQ), B4RFR(IFREQ), C4RFR(IFREQ), D4RFR(IFREQ),
     &        E4RFR(IFREQ),
     &        HALF*(HYPERPOL(IFREQ,IOPER,1)+HYPERPOL(IFREQ,IOPER,2))
           END IF
          END IF
         END DO

      END DO

      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC4HYPPRT                            *
*---------------------------------------------------------------------*
c /* deck CC4R_SETUP */
*=====================================================================*
      SUBROUTINE CC4R_SETUP(MXTRAN3, MXVEC2, MXTRAN4, MXVEC1,
     &                      I0HTRAN, I0HDOTS, W0H, N0HTRAN,
     &                      I1HTRAN, I1HDOTS, W1H, N1HTRAN,
     &                      I0GTRAN, I0GDOTS, W0G, N0GTRAN,
     &                      I1GTRAN, I1GDOTS, W1G, N1GTRAN,
     &                      I2GTRAN, I2GDOTS, W2G, N2GTRAN,
     &                      I1FTRAN, I1FDOTS, W1F, N1FTRAN,
     &                      I2FTRAN, I2FDOTS, W2F, N2FTRAN,
     &                      I0FATRAN,I0FADOTS,W0FA,N0FATRAN,
     &                      I1FATRAN,I1FADOTS,W1FA,N1FATRAN,
     &                      I2FATRAN,I2FADOTS,W2FA,N2FATRAN,
     &                      I2EATRAN,I2EADOTS,W2EA,N2EATRAN,
     &                      IXTRAN,  IXDOTS,  WX,  NXTRAN,
     &                      HYPPOL,  MXHYPPOL,LADD     )
*---------------------------------------------------------------------*
*
*    Purpose: set up for CC4R section
*                - list of H^0  matrix transformations 
*                - list of H^A  matrix transformations 
*                - list of G^0  matrix transformations 
*                - list of G^A  matrix transformations 
*                - list of G^AB matrix transformations 
*                - list of F^A  matrix transformations 
*                - list of F^AB matrix transformations 
*                - list of F^0{O}  matrix transformations 
*                - list of F^A{O}  matrix transformations 
*                - list of F^AB{O} matrix transformations 
*                - list of ETA^AB{O} vector calculations 
*                - list of chi vector dot products
*
*     LADD = .FALSE.  --> set lists
*     LADD = .TRUE.   --> add contributions to hyperpolarizabilities
*
*     Written by Christof Haettig, april 1997.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <cc4rinf.h>
#include <cc4perm.h>
#include <ccroper.h>

* local parameters:
      CHARACTER*(20) MSGDBG
      PARAMETER (MSGDBG = '[debug] CC4R_SETUP> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      LOGICAL LADD

      INTEGER MXVEC2, MXTRAN3, MXVEC1, MXTRAN4, MXHYPPOL

      INTEGER I0HTRAN(5,MXTRAN3)
      INTEGER I0HDOTS(MXVEC2,MXTRAN3)
      INTEGER I1HTRAN(5,MXTRAN4)
      INTEGER I1HDOTS(MXVEC1,MXTRAN4)

      INTEGER I0GTRAN(4,MXTRAN3)
      INTEGER I0GDOTS(MXVEC2,MXTRAN3)
      INTEGER I1GTRAN(4,MXTRAN3)
      INTEGER I1GDOTS(MXVEC2,MXTRAN3)
      INTEGER I2GTRAN(4,MXTRAN4)
      INTEGER I2GDOTS(MXVEC1,MXTRAN4)

      INTEGER I1FTRAN(3,MXTRAN3)
      INTEGER I1FDOTS(MXVEC2,MXTRAN3)
      INTEGER I2FTRAN(3,MXTRAN3)
      INTEGER I2FDOTS(MXVEC2,MXTRAN3)
  
      INTEGER I0FATRAN(5,MXTRAN3)
      INTEGER I0FADOTS(MXVEC2,MXTRAN3)
      INTEGER I1FATRAN(5,MXTRAN3)
      INTEGER I1FADOTS(MXVEC2,MXTRAN3)
      INTEGER I2FATRAN(5,MXTRAN4)
      INTEGER I2FADOTS(MXVEC1,MXTRAN4)

      INTEGER I2EATRAN(3,MXTRAN3)
      INTEGER I2EADOTS(MXVEC2,MXTRAN3)

      INTEGER IXTRAN(MXTRAN3)
      INTEGER IXDOTS(MXVEC2,MXTRAN3)

      INTEGER N0HTRAN, N0GTRAN,          N0FATRAN
      INTEGER N1HTRAN, N1GTRAN, N1FTRAN, N1FATRAN
      INTEGER          N2GTRAN, N2FTRAN, N2FATRAN, N2EATRAN, NXTRAN
      INTEGER N4RRESF
      INTEGER MXVH0, MXVH1, MXVG0, MXVG1, MXVG2, MXVF1, MXVF2
      INTEGER MXVFA0, MXVFA1, MXVFA2, MXVEA2, MXX

      INTEGER IVEC, ITRAN, I, IDX
      INTEGER ISYML, ISYM1, ISYM2, ISYTOT, ISYM3
      INTEGER IFREQ, IOPER
      INTEGER P, ISIGN

#if defined (SYS_CRAY)
      REAL HYPPOL(MXHYPPOL)
      REAL SIGN, SUM, SSUM
      REAL H0CON(10),  H1CON(5)
      REAL G0CON(15),  G1CON(30),  G2CON(10)
      REAL             F1CON(15),  F2CON(30)
      REAL F0ACON(15), F1ACON(60), F2ACON(30)
      REAL                         E2ACON(30)
      REAL XCON(10)
      REAL W0H(MXVEC2,MXTRAN3)
      REAL W1H(MXVEC1,MXTRAN4)
      REAL W0G(MXVEC2,MXTRAN3)
      REAL W1G(MXVEC2,MXTRAN3)
      REAL W2G(MXVEC1,MXTRAN4)
      REAL W1F(MXVEC2,MXTRAN3)
      REAL W2F(MXVEC2,MXTRAN3)
      REAL W0FA(MXVEC2,MXTRAN3)
      REAL W1FA(MXVEC2,MXTRAN3)
      REAL W2FA(MXVEC1,MXTRAN4)
      REAL W2EA(MXVEC2,MXTRAN3)
      REAL WX(MXVEC2,MXTRAN3)
#else
      DOUBLE PRECISION SIGN, SUM, SSUM
      DOUBLE PRECISION HYPPOL(MXHYPPOL)
      DOUBLE PRECISION H0CON(10),  H1CON(5)
      DOUBLE PRECISION G0CON(15),  G1CON(30),  G2CON(10)
      DOUBLE PRECISION             F1CON(15),  F2CON(30)
      DOUBLE PRECISION F0ACON(15), F1ACON(60), F2ACON(30)
      DOUBLE PRECISION                         E2ACON(30)
      DOUBLE PRECISION XCON(10)
      DOUBLE PRECISION W0H(MXVEC2,MXTRAN3)
      DOUBLE PRECISION W1H(MXVEC1,MXTRAN4)
      DOUBLE PRECISION W0G(MXVEC2,MXTRAN3)
      DOUBLE PRECISION W1G(MXVEC2,MXTRAN3)
      DOUBLE PRECISION W2G(MXVEC1,MXTRAN4)
      DOUBLE PRECISION W1F(MXVEC2,MXTRAN3)
      DOUBLE PRECISION W2F(MXVEC2,MXTRAN3)
      DOUBLE PRECISION W0FA(MXVEC2,MXTRAN3)
      DOUBLE PRECISION W1FA(MXVEC2,MXTRAN3)
      DOUBLE PRECISION W2FA(MXVEC1,MXTRAN4)
      DOUBLE PRECISION W2EA(MXVEC2,MXTRAN3)
      DOUBLE PRECISION WX(MXVEC2,MXTRAN3)
#endif


      INTEGER IOP(5), ISY(5), IL1(5), IR1(5), IR2(10), IL2(10)
      INTEGER IX3(10)


* external functions:
      INTEGER IR2TAMP
      INTEGER IR1TAMP
      INTEGER IL1ZETA
      INTEGER IL2ZETA
      INTEGER ICHI3


*---------------------------------------------------------------------*
* initializations:
*---------------------------------------------------------------------*
      IF (.NOT. LADD) THEN
        N0HTRAN  = 0
        N1HTRAN  = 0
        N0GTRAN  = 0
        N1GTRAN  = 0
        N2GTRAN  = 0
        N1FTRAN  = 0
        N2FTRAN  = 0
        N0FATRAN = 0
        N1FATRAN = 0
        N2FATRAN = 0
        N2EATRAN = 0
        NXTRAN   = 0
      END IF

      MXVH0  = 0
      MXVH1  = 0
      MXVG0  = 0
      MXVG1  = 0
      MXVG2  = 0
      MXVF1  = 0
      MXVF2  = 0
      MXVFA0 = 0
      MXVFA1 = 0
      MXVFA2 = 0
      MXVEA2 = 0
      MXX    = 0

      N4RRESF  = 0

      CALL DZERO(HYPPOL,MXHYPPOL)
 
*---------------------------------------------------------------------*
* start loop over all requested third hyperpolarizabilities
*---------------------------------------------------------------------*
 
      DO IOPER = 1, N4ROPER
        IOP(A) = IA4ROP(IOPER)
        IOP(B) = IB4ROP(IOPER)
        IOP(C) = IC4ROP(IOPER)
        IOP(D) = ID4ROP(IOPER)
        IOP(E) = IE4ROP(IOPER)

        ISY(A) = ISYOPR(IOP(A))
        ISY(B) = ISYOPR(IOP(B))
        ISY(C) = ISYOPR(IOP(C))
        ISY(D) = ISYOPR(IOP(D))
        ISY(E) = ISYOPR(IOP(E))
 
        ISYTOT = 1
        DO I = A, E
          ISYTOT = MULD2H(ISYTOT,ISY(I))
        END DO

      IF (ISYTOT .EQ. 1) THEN

      DO IFREQ = 1, N4RFREQ

        N4RRESF = N4RRESF + 1
      
      DO ISIGN = 1, -1, -2
        SIGN = DBLE(ISIGN)

        IL1(A) =IL1ZETA(LBLOPR(IOP(A)),.FALSE.,SIGN*A4RFR(IFREQ),ISYML)
        IL1(B) =IL1ZETA(LBLOPR(IOP(B)),.FALSE.,SIGN*B4RFR(IFREQ),ISYML)
        IL1(C) =IL1ZETA(LBLOPR(IOP(C)),.FALSE.,SIGN*C4RFR(IFREQ),ISYML)
        IL1(D) =IL1ZETA(LBLOPR(IOP(D)),.FALSE.,SIGN*D4RFR(IFREQ),ISYML)
        IL1(E) =IL1ZETA(LBLOPR(IOP(E)),.FALSE.,SIGN*E4RFR(IFREQ),ISYML)

        IR1(A) =IR1TAMP(LBLOPR(IOP(A)),.FALSE.,SIGN*A4RFR(IFREQ),ISYML)
        IR1(B) =IR1TAMP(LBLOPR(IOP(B)),.FALSE.,SIGN*B4RFR(IFREQ),ISYML)
        IR1(C) =IR1TAMP(LBLOPR(IOP(C)),.FALSE.,SIGN*C4RFR(IFREQ),ISYML)
        IR1(D) =IR1TAMP(LBLOPR(IOP(D)),.FALSE.,SIGN*D4RFR(IFREQ),ISYML)
        IR1(E) =IR1TAMP(LBLOPR(IOP(E)),.FALSE.,SIGN*E4RFR(IFREQ),ISYML)

        IR2(AB)=IR2TAMP(LBLOPR(IOP(A)),.FALSE.,SIGN*A4RFR(IFREQ),ISYM1,
     &                  LBLOPR(IOP(B)),.FALSE.,SIGN*B4RFR(IFREQ),ISYM2)
        IR2(AC)=IR2TAMP(LBLOPR(IOP(A)),.FALSE.,SIGN*A4RFR(IFREQ),ISYM1,
     &                  LBLOPR(IOP(C)),.FALSE.,SIGN*C4RFR(IFREQ),ISYM2)
        IR2(AD)=IR2TAMP(LBLOPR(IOP(A)),.FALSE.,SIGN*A4RFR(IFREQ),ISYM1,
     &                  LBLOPR(IOP(D)),.FALSE.,SIGN*D4RFR(IFREQ),ISYM2)
        IR2(AE)=IR2TAMP(LBLOPR(IOP(A)),.FALSE.,SIGN*A4RFR(IFREQ),ISYM1,
     &                  LBLOPR(IOP(E)),.FALSE.,SIGN*E4RFR(IFREQ),ISYM2)
        IR2(BC)=IR2TAMP(LBLOPR(IOP(B)),.FALSE.,SIGN*B4RFR(IFREQ),ISYM1,
     &                  LBLOPR(IOP(C)),.FALSE.,SIGN*C4RFR(IFREQ),ISYM2)
        IR2(BD)=IR2TAMP(LBLOPR(IOP(B)),.FALSE.,SIGN*B4RFR(IFREQ),ISYM1,
     &                  LBLOPR(IOP(D)),.FALSE.,SIGN*D4RFR(IFREQ),ISYM2)
        IR2(BE)=IR2TAMP(LBLOPR(IOP(B)),.FALSE.,SIGN*B4RFR(IFREQ),ISYM1,
     &                  LBLOPR(IOP(E)),.FALSE.,SIGN*E4RFR(IFREQ),ISYM2)
        IR2(CD)=IR2TAMP(LBLOPR(IOP(C)),.FALSE.,SIGN*C4RFR(IFREQ),ISYM1,
     &                  LBLOPR(IOP(D)),.FALSE.,SIGN*D4RFR(IFREQ),ISYM2)
        IR2(CE)=IR2TAMP(LBLOPR(IOP(C)),.FALSE.,SIGN*C4RFR(IFREQ),ISYM1,
     &                  LBLOPR(IOP(E)),.FALSE.,SIGN*E4RFR(IFREQ),ISYM2)
        IR2(DE)=IR2TAMP(LBLOPR(IOP(D)),.FALSE.,SIGN*D4RFR(IFREQ),ISYM1,
     &                  LBLOPR(IOP(E)),.FALSE.,SIGN*E4RFR(IFREQ),ISYM2)

        IL2(AB) = IL2ZETA(LBLOPR(IOP(A)),SIGN*A4RFR(IFREQ),ISYM1,
     &                    LBLOPR(IOP(B)),SIGN*B4RFR(IFREQ),ISYM2)
        IL2(AC) = IL2ZETA(LBLOPR(IOP(A)),SIGN*A4RFR(IFREQ),ISYM1,
     &                    LBLOPR(IOP(C)),SIGN*C4RFR(IFREQ),ISYM2)
        IL2(AD) = IL2ZETA(LBLOPR(IOP(A)),SIGN*A4RFR(IFREQ),ISYM1,
     &                    LBLOPR(IOP(D)),SIGN*D4RFR(IFREQ),ISYM2)
        IL2(AE) = IL2ZETA(LBLOPR(IOP(A)),SIGN*A4RFR(IFREQ),ISYM1,
     &                    LBLOPR(IOP(E)),SIGN*E4RFR(IFREQ),ISYM2)
        IL2(BC) = IL2ZETA(LBLOPR(IOP(B)),SIGN*B4RFR(IFREQ),ISYM1,
     &                    LBLOPR(IOP(C)),SIGN*C4RFR(IFREQ),ISYM2)
        IL2(BD) = IL2ZETA(LBLOPR(IOP(B)),SIGN*B4RFR(IFREQ),ISYM1,
     &                    LBLOPR(IOP(D)),SIGN*D4RFR(IFREQ),ISYM2)
        IL2(BE) = IL2ZETA(LBLOPR(IOP(B)),SIGN*B4RFR(IFREQ),ISYM1,
     &                    LBLOPR(IOP(E)),SIGN*E4RFR(IFREQ),ISYM2)
        IL2(CD) = IL2ZETA(LBLOPR(IOP(C)),SIGN*C4RFR(IFREQ),ISYM1,
     &                    LBLOPR(IOP(D)),SIGN*D4RFR(IFREQ),ISYM2)
        IL2(CE) = IL2ZETA(LBLOPR(IOP(C)),SIGN*C4RFR(IFREQ),ISYM1,
     &                    LBLOPR(IOP(E)),SIGN*E4RFR(IFREQ),ISYM2)
        IL2(DE) = IL2ZETA(LBLOPR(IOP(D)),SIGN*D4RFR(IFREQ),ISYM1,
     &                    LBLOPR(IOP(E)),SIGN*E4RFR(IFREQ),ISYM2)

      IF (L_USE_CHI3) THEN 
        IX3(CDE) = ICHI3(LBLOPR(IOP(C)),SIGN*C4RFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(D)),SIGN*D4RFR(IFREQ),ISYM2,
     &                   LBLOPR(IOP(E)),SIGN*E4RFR(IFREQ),ISYM3)
        IX3(BDE) = ICHI3(LBLOPR(IOP(B)),SIGN*B4RFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(D)),SIGN*D4RFR(IFREQ),ISYM2,
     &                   LBLOPR(IOP(E)),SIGN*E4RFR(IFREQ),ISYM3)
        IX3(BCE) = ICHI3(LBLOPR(IOP(B)),SIGN*B4RFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(C)),SIGN*C4RFR(IFREQ),ISYM2,
     &                   LBLOPR(IOP(E)),SIGN*E4RFR(IFREQ),ISYM3)
        IX3(BCD) = ICHI3(LBLOPR(IOP(B)),SIGN*B4RFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(C)),SIGN*C4RFR(IFREQ),ISYM2,
     &                   LBLOPR(IOP(D)),SIGN*D4RFR(IFREQ),ISYM3)
        IX3(ADE) = ICHI3(LBLOPR(IOP(A)),SIGN*A4RFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(D)),SIGN*D4RFR(IFREQ),ISYM2,
     &                   LBLOPR(IOP(E)),SIGN*E4RFR(IFREQ),ISYM3)
        IX3(ACE) = ICHI3(LBLOPR(IOP(A)),SIGN*A4RFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(C)),SIGN*C4RFR(IFREQ),ISYM2,
     &                   LBLOPR(IOP(E)),SIGN*E4RFR(IFREQ),ISYM3)
        IX3(ACD) = ICHI3(LBLOPR(IOP(A)),SIGN*A4RFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(C)),SIGN*C4RFR(IFREQ),ISYM2,
     &                   LBLOPR(IOP(D)),SIGN*D4RFR(IFREQ),ISYM3)
        IX3(ABE) = ICHI3(LBLOPR(IOP(A)),SIGN*A4RFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(B)),SIGN*B4RFR(IFREQ),ISYM2,
     &                   LBLOPR(IOP(E)),SIGN*E4RFR(IFREQ),ISYM3)
        IX3(ABD) = ICHI3(LBLOPR(IOP(A)),SIGN*A4RFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(B)),SIGN*B4RFR(IFREQ),ISYM2,
     &                   LBLOPR(IOP(D)),SIGN*D4RFR(IFREQ),ISYM3)
        IX3(ABC) = ICHI3(LBLOPR(IOP(A)),SIGN*A4RFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(B)),SIGN*B4RFR(IFREQ),ISYM2,
     &                   LBLOPR(IOP(C)),SIGN*C4RFR(IFREQ),ISYM3)
      END IF
*---------------------------------------------------------------------*
* set up list of H^0 matrix transformations, 10 permutations
*---------------------------------------------------------------------*
      DO P = 1, 5
        CALL CC_SETH1112(I0HTRAN,I0HDOTS,MXTRAN3,MXVEC2,
     &    0,IR1(J3(P)),IR1(J4(P)),IR1(J5(P)),IR2(JP1(P)),ITRAN,IVEC)
        N0HTRAN  = MAX(N0HTRAN,ITRAN)
        MXVH0    = MAX(MXVH0,IVEC)
        H0CON(P) = W0H(IVEC,ITRAN)

        CALL CC_SETH1112(I0HTRAN,I0HDOTS,MXTRAN3,MXVEC2,
     &    0,IR1(J2(P)),IR1(J4(P)),IR1(J5(P)),IR2(JP2(P)),ITRAN,IVEC)
        N0HTRAN    = MAX(N0HTRAN,ITRAN)
        MXVH0      = MAX(MXVH0,IVEC)
        H0CON(P+5) = W0H(IVEC,ITRAN)
      END DO

*---------------------------------------------------------------------*
* set up list of H^A matrix transformations, 5 permutations
*---------------------------------------------------------------------*
      DO P = 1, 5
        CALL CC_SETH1111(I1HTRAN,I1HDOTS,MXTRAN4,MXVEC1,
     &    IL1(J1(P)),IR1(J2(P)),IR1(J3(P)),IR1(J4(P)),IR1(J5(P)),
     &    ITRAN,IVEC)
        N1HTRAN  = MAX(N1HTRAN,ITRAN)
        MXVH1    = MAX(MXVH1,IVEC)
        H1CON(P) = W1H(IVEC,ITRAN)
      END DO

*---------------------------------------------------------------------*
* set up list of G^0 matrix transformations, 15 permutations
*---------------------------------------------------------------------*
      DO P = 1, 15
        CALL CC_SETG212(I0GTRAN,I0GDOTS,MXTRAN3,MXVEC2,   
     &              0,IR2(IP1(P)),IR1(I1(P)),IR2(IP2(P)),ITRAN,IVEC)
        N0GTRAN  = MAX(N0GTRAN,ITRAN)
        MXVG0    = MAX(MXVG0,IVEC)
        G0CON(P) = W0G(IVEC,ITRAN)
      END DO

*---------------------------------------------------------------------*
* set up list of G^A matrix transformations, 30 permutations
*---------------------------------------------------------------------*
      DO P = 1, 15
        CALL CC_SETG112(I1GTRAN,I1GDOTS,MXTRAN3,MXVEC2,   
     &     IL1(I1(P)),IR1(I2(P)),IR1(I3(P)),IR2(IP2(P)),ITRAN,IVEC)
        N1GTRAN  = MAX(N1GTRAN,ITRAN)
        MXVG1    = MAX(MXVG1,IVEC)
        G1CON(P) = W1G(IVEC,ITRAN)

        CALL CC_SETG112(I1GTRAN,I1GDOTS,MXTRAN3,MXVEC2,   
     &     IL1(I1(P)),IR1(I4(P)),IR1(I5(P)),IR2(IP1(P)),ITRAN,IVEC)
        N1GTRAN     = MAX(N1GTRAN,ITRAN)
        MXVG1       = MAX(MXVG1,IVEC)
        G1CON(P+15) = W1G(IVEC,ITRAN)
      END DO

*---------------------------------------------------------------------*
* set up list of G^AB matrix transformations, 10 permutations
*---------------------------------------------------------------------*
      DO P = 1, 5
        CALL CCQR_SETG(I2GTRAN,I2GDOTS,MXTRAN4,MXVEC1,   
     &     IL2(JP1(P)),IR1(J3(P)),IR1(J4(P)),IR1(J5(P)),ITRAN,IVEC)
        N2GTRAN  = MAX(N2GTRAN,ITRAN)
        MXVG2    = MAX(MXVG2,IVEC)
        G2CON(P) = W2G(IVEC,ITRAN)

        CALL CCQR_SETG(I2GTRAN,I2GDOTS,MXTRAN4,MXVEC1,   
     &     IL2(JP2(P)),IR1(J2(P)),IR1(J4(P)),IR1(J5(P)),ITRAN,IVEC)
        N2GTRAN    = MAX(N2GTRAN,ITRAN)
        MXVG2      = MAX(MXVG2,IVEC)
        G2CON(P+5) = W2G(IVEC,ITRAN)
      END DO

*---------------------------------------------------------------------*
* set up list of F^A matrix transformations, 15 permutations
*---------------------------------------------------------------------*
      DO P = 1, 15
        CALL CCQR_SETF(I1FTRAN,I1FDOTS,MXTRAN3,MXVEC2, 
     &               IL1(I1(P)),IR2(IP1(P)),IR2(IP2(P)),ITRAN,IVEC)
        N1FTRAN  = MAX(N1FTRAN,ITRAN)
        MXVF1    = MAX(MXVF1,IVEC)
        F1CON(P) = W1F(IVEC,ITRAN)
      END DO

*---------------------------------------------------------------------*
* set up list of F^AB matrix transformations, 30 permutations
*---------------------------------------------------------------------*
      DO P = 1, 15
        CALL CC_SETF12(I2FTRAN,I2FDOTS,MXTRAN3,MXVEC2, 
     &               IL2(IP1(P)),IR1(I1(P)),IR2(IP2(P)),ITRAN,IVEC)
        N2FTRAN  = MAX(N2FTRAN,ITRAN)
        MXVF2    = MAX(MXVF2,IVEC)
        F2CON(P) = W2F(IVEC,ITRAN)

        CALL CC_SETF12(I2FTRAN,I2FDOTS,MXTRAN3,MXVEC2, 
     &               IL2(IP2(P)),IR1(I1(P)),IR2(IP1(P)),ITRAN,IVEC)
        N2FTRAN     = MAX(N2FTRAN,ITRAN)
        MXVF2       = MAX(MXVF2,IVEC)
        F2CON(P+15) = W2F(IVEC,ITRAN)
      END DO

*---------------------------------------------------------------------*
* set up list of F^0{O} matrix transformations, 15 permutations
*---------------------------------------------------------------------*
      DO P = 1, 15
        CALL CCQR_SETFA(I0FATRAN,I0FADOTS,MXTRAN3,MXVEC2, 
     &               0,IOP(I1(P)),IR2(IP1(P)),IR2(IP2(P)),ITRAN,IVEC)
        N0FATRAN  = MAX(N0FATRAN,ITRAN)
        MXVFA0    = MAX(MXVFA0,IVEC)
        F0ACON(P) = W0FA(IVEC,ITRAN)
      END DO

*---------------------------------------------------------------------*
* set up list of F^A{O} matrix transformations, 60 permutations
*---------------------------------------------------------------------*
      DO P = 1, 15
        CALL CC_SETFA12(I1FATRAN,I1FADOTS,MXTRAN3,MXVEC2,
     &      IL1(I1(P)),IOP(I2(P)),IR1(I3(P)),IR2(IP2(P)),ITRAN,IVEC)
        N1FATRAN  = MAX(N1FATRAN,ITRAN)
        MXVFA1    = MAX(MXVFA1,IVEC)
        F1ACON(P) = W1FA(IVEC,ITRAN)

        CALL CC_SETFA12(I1FATRAN,I1FADOTS,MXTRAN3,MXVEC2,
     &      IL1(I1(P)),IOP(I3(P)),IR1(I2(P)),IR2(IP2(P)),ITRAN,IVEC)
        N1FATRAN     = MAX(N1FATRAN,ITRAN)
        MXVFA1       = MAX(MXVFA1,IVEC)
        F1ACON(P+15) = W1FA(IVEC,ITRAN)

        CALL CC_SETFA12(I1FATRAN,I1FADOTS,MXTRAN3,MXVEC2,
     &      IL1(I1(P)),IOP(I4(P)),IR1(I5(P)),IR2(IP1(P)),ITRAN,IVEC)
        N1FATRAN     = MAX(N1FATRAN,ITRAN)
        MXVFA1       = MAX(MXVFA1,IVEC)
        F1ACON(P+30) = W1FA(IVEC,ITRAN)

        CALL CC_SETFA12(I1FATRAN,I1FADOTS,MXTRAN3,MXVEC2,
     &      IL1(I1(P)),IOP(I5(P)),IR1(I4(P)),IR2(IP1(P)),ITRAN,IVEC)
        N1FATRAN     = MAX(N1FATRAN,ITRAN)
        MXVFA1       = MAX(MXVFA1,IVEC)
        F1ACON(P+45) = W1FA(IVEC,ITRAN)
      END DO

*---------------------------------------------------------------------*
* set up list of F^AB{O} matrix transformations, 30 permutations
*---------------------------------------------------------------------*
      DO P = 1, 15
        CALL CCQR_SETFA(I2FATRAN,I2FADOTS,MXTRAN4,MXVEC1, 
     &        IL2(IP1(P)),IOP(I1(P)),IR1(I4(P)),IR1(I5(P)),ITRAN,IVEC)
        N2FATRAN  = MAX(N2FATRAN,ITRAN)
        MXVFA2    = MAX(MXVFA2,IVEC)
        F2ACON(P) = W2FA(IVEC,ITRAN)

        CALL CCQR_SETFA(I2FATRAN,I2FADOTS,MXTRAN4,MXVEC1, 
     &        IL2(IP2(P)),IOP(I1(P)),IR1(I2(P)),IR1(I3(P)),ITRAN,IVEC)
        N2FATRAN     = MAX(N2FATRAN,ITRAN)
        MXVFA2       = MAX(MXVFA2,IVEC)
        F2ACON(P+15) = W2FA(IVEC,ITRAN)
      END DO

*---------------------------------------------------------------------*
* set up list of ETA{O} vector calculations, 30 permutations
*---------------------------------------------------------------------*
      DO P = 1, 15
        CALL CCQR_SETEA(I2EATRAN,I2EADOTS,MXTRAN3,MXVEC2, 
     &                  IL2(IP1(P)),IOP(I1(P)),IR2(IP2(P)),ITRAN,IVEC)
        N2EATRAN  = MAX(N2EATRAN,ITRAN)
        MXVEA2    = MAX(MXVEA2,IVEC)
        E2ACON(P) = W2EA(IVEC,ITRAN)

        CALL CCQR_SETEA(I2EATRAN,I2EADOTS,MXTRAN3,MXVEC2,
     &                  IL2(IP2(P)),IOP(I1(P)),IR2(IP1(P)),ITRAN,IVEC)
        N2EATRAN     = MAX(N2EATRAN,ITRAN)
        MXVEA2       = MAX(MXVEA2,IVEC)
        E2ACON(P+15) = W2EA(IVEC,ITRAN)
      END DO

*---------------------------------------------------------------------*
* set up list of chi vector dot products, 6 permutations
*---------------------------------------------------------------------*
      IF (L_USE_CHI3) THEN 
        DO P = 1, 10
          CALL CC_SETDOT(IXTRAN,IXDOTS,MXTRAN3,MXVEC2,
     &                   IX3(11-P),IR2(P),ITRAN,IVEC)
          NXTRAN  = MAX(NXTRAN,ITRAN)
          MXX     = MAX(MXX,IVEC)
          XCON(P) = WX(IVEC,ITRAN)
        END DO
      ELSE
        DO P = 1, 10
          XCON(P) = 0.0d0
        END DO
      END IF


*---------------------------------------------------------------------*
* add all 250 contributions to the hyperpolarizabilities:
*---------------------------------------------------------------------*
      IF (LADD) THEN
        IDX = N4RFREQ*(IOPER-1) + IFREQ
        IF (ISIGN.EQ.-1) IDX = IDX + N4RFREQ*N4ROPER

        DO P = 1, 5
         HYPPOL(IDX) = HYPPOL(IDX) + H0CON(P) + H0CON(P+5) + H1CON(P) +
     &                               G2CON(P) + G2CON(P+5) +
     &                               XCON(P) + XCON(P+5)
        END DO

        DO P = 1, 15
         HYPPOL(IDX) = HYPPOL(IDX) + 
     &    G0CON(P)  + G1CON(P)     + G1CON(P+15) +
     &    F1CON(P)  + F2CON(P)     + F2CON(P+15)  + F0ACON(P) + 
     &    F1ACON(P) + F1ACON(P+15) + F1ACON(P+30) + F1ACON(P+45) + 
     &    F2ACON(P) + F2ACON(P+15) + E2ACON(P)    + E2ACON(P+15)
        END DO

        IF (LOCDBG) THEN
          WRITE (LUPRI,*) 'IOP:',IOP
          WRITE (LUPRI,*) 'IDX:',IDX

          SSUM= 0.0d0
          SUM = 0.0d0
          DO I = 1, 10
            SUM = SUM + H0CON(I)
          END DO
          SSUM = SSUM + SUM
          WRITE (LUPRI,*) 'H0CON:',H0CON
          WRITE (LUPRI,*) 'SUM:', SUM, SSUM

          SUM = 0.0d0
          DO I = 1, 5
            SUM = SUM + H1CON(I)
          END DO
          SSUM = SSUM + SUM
          WRITE (LUPRI,*) 'H1CON:',H1CON
          WRITE (LUPRI,*) 'SUM:', SUM, SSUM


          SUM = 0.0d0
          DO I = 1, 15
            SUM = SUM + G0CON(I)
          END DO
          SSUM = SSUM + SUM
          WRITE (LUPRI,*) 'G0CON:',G0CON
          WRITE (LUPRI,*) 'SUM:', SUM, SSUM


          SUM = 0.0d0
          DO I = 1, 30
            SUM = SUM + G1CON(I)
          END DO
          SSUM = SSUM + SUM
          WRITE (LUPRI,*) 'G1CON:',G1CON
          WRITE (LUPRI,*) 'SUM:', SUM, SSUM

          SUM = 0.0d0
          DO I = 1, 10
            SUM = SUM + G2CON(I)
          END DO
          SSUM = SSUM + SUM
          WRITE (LUPRI,*) 'G2CON:',G2CON
          WRITE (LUPRI,*) 'SUM:', SUM, SSUM

          SUM = 0.0d0
          DO I = 1, 15
            SUM = SUM + F1CON(I)
          END DO
          SSUM = SSUM + SUM
          WRITE (LUPRI,*) 'F1CON:',F1CON
          WRITE (LUPRI,*) 'SUM:', SUM, SSUM

          SUM = 0.0d0
          DO I = 1, 30
            SUM = SUM + F2CON(I)
          END DO
          SSUM = SSUM + SUM
          WRITE (LUPRI,*) 'F2CON:',F2CON
          WRITE (LUPRI,*) 'SUM:', SUM, SSUM


          SUM = 0.0d0
          DO I = 1, 15
            SUM = SUM + F0ACON(I)
          END DO
          SSUM = SSUM + SUM
          WRITE (LUPRI,*) 'F0ACON:',F0ACON
          WRITE (LUPRI,*) 'SUM:', SUM, SSUM


          SUM = 0.0d0
          DO I = 1, 60
            SUM = SUM + F1ACON(I)
          END DO
          SSUM = SSUM + SUM
          WRITE (LUPRI,*) 'F1ACON:',F1ACON
          WRITE (LUPRI,*) 'SUM:', SUM, SSUM


          SUM = 0.0d0
          DO I = 1, 30
            SUM = SUM + F2ACON(I)
          END DO
          SSUM = SSUM + SUM
          WRITE (LUPRI,*) 'F2ACON:',F2ACON
          WRITE (LUPRI,*) 'SUM:', SUM, SSUM


          SUM = 0.0d0
          DO I = 1, 30
            SUM = SUM + E2ACON(I)
          END DO
          SSUM = SSUM + SUM
          WRITE (LUPRI,*) 'E2ACON:',E2ACON
          WRITE (LUPRI,*) 'SUM:', SUM, SSUM

          SUM = 0.0d0
          DO I = 1, 30
            SUM = SUM + XCON(I)
          END DO
          SSUM = SSUM + SUM
          WRITE (LUPRI,*) 'XCON:',XCON
          WRITE (LUPRI,*) 'SUM:', SUM, SSUM

          WRITE (LUPRI,*) 'HYPPOL:',HYPPOL(IDX)
        END IF

      END IF

*---------------------------------------------------------------------*
* end loop over all requested hyperpolarizabilities
*---------------------------------------------------------------------*
      END DO
      END DO
      END IF
      END DO

*---------------------------------------------------------------------*
* print statistics and lists: 
*---------------------------------------------------------------------*
      IF (.NOT. LADD) THEN

* general statistics:
      WRITE(LUPRI,'(/,/3X,A,I4,A)') 'For the requested',N4RRESF,
     &      ' quartic response functions '
      WRITE(LUPRI,'((8X,A,I3,A))') 
     &   ' - ',N0HTRAN,  ' H^0 matrix transformations ',
     &   ' - ',N1HTRAN,  ' H^A matrix transformations ',
     &   ' - ',N0GTRAN,  ' G^0 matrix transformations ',
     &   ' - ',N1GTRAN,  ' G^A matrix transformations ',
     &   ' - ',N2GTRAN,  ' G^AB matrix transformations ',
     &   ' - ',N1FTRAN,  ' F^A matrix transformations ',
     &   ' - ',N2FTRAN,  ' F^AB matrix transformations ',
     &   ' - ',N0FATRAN, ' F^0{O} matrix transformations ',
     &   ' - ',N1FATRAN, ' F^A{O} matrix transformations ',
     &   ' - ',N2FATRAN, ' F^AB{O} matrix transformations ',
     &   ' - ',N2EATRAN, ' ETA^AB{O} vector calculations ' 
      WRITE(LUPRI,'(3X,A,/,/)') 'will be performed.'


* H^0 matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of H^0 matrix transformations:'
        DO ITRAN = 1, N0HTRAN
          WRITE(LUPRI,'(A,4I5,5X,(12I5,20X))') MSGDBG,
     &     (I0HTRAN(I,ITRAN),I=1,4),(I0HDOTS(I,ITRAN),I=1,MXVH0)
        END DO
        WRITE (LUPRI,*)
      END IF

* H^A matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of H^A matrix transformations:'
        DO ITRAN = 1, N1HTRAN
          WRITE(LUPRI,'(A,4I5,5X,(12I5,20X))') MSGDBG,
     &     (I1HTRAN(I,ITRAN),I=1,4),(I1HDOTS(I,ITRAN),I=1,MXVH1)
        END DO
        WRITE (LUPRI,*)
      END IF

* G^0 matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of G^0 matrix transformations:'
        DO ITRAN = 1, N0GTRAN
          WRITE(LUPRI,'(A,3I5,5X,(12I5,20X))') MSGDBG,
     &     (I0GTRAN(I,ITRAN),I=1,3),(I0GDOTS(I,ITRAN),I=1,MXVG0)
        END DO
        WRITE (LUPRI,*)
      END IF

* G^A matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of G^A matrix transformations:'
        DO ITRAN = 1, N1GTRAN
          WRITE(LUPRI,'(A,3I5,5X,(12I5,20X))') MSGDBG,
     &     (I1GTRAN(I,ITRAN),I=1,3),(I1GDOTS(I,ITRAN),I=1,MXVG1)
        END DO
        WRITE (LUPRI,*)
      END IF

* G^AB matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of G^AB matrix transformations:'
        DO ITRAN = 1, N2GTRAN
          WRITE(LUPRI,'(A,3I5,5X,(12I5,20X))') MSGDBG,
     &     (I2GTRAN(I,ITRAN),I=1,3),(I2GDOTS(I,ITRAN),I=1,MXVG2)
        END DO
        WRITE (LUPRI,*)
      END IF

* F^A matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of F^A matrix transformations:'
        DO ITRAN = 1, N1FTRAN
          WRITE(LUPRI,'(A,2I5,5X,(12I5,20X))') MSGDBG,
     &     (I1FTRAN(I,ITRAN),I=1,2),(I1FDOTS(I,ITRAN),I=1,MXVF1)
        END DO
        WRITE (LUPRI,*)
      END IF

* F^AB matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of F^AB matrix transformations:'
        DO ITRAN = 1, N2FTRAN
          WRITE(LUPRI,'(A,2I5,5X,(12I5,20X))') MSGDBG,
     &     (I2FTRAN(I,ITRAN),I=1,2),(I2FDOTS(I,ITRAN),I=1,MXVF2)
        END DO
        WRITE (LUPRI,*)
      END IF

* F^0{O} matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of F{O} matrix transformations:'
        DO ITRAN = 1, N0FATRAN
          WRITE(LUPRI,'(A,5I5,5X,(12I5,20X))') MSGDBG,
     &     (I0FATRAN(I,ITRAN),I=1,5),(I0FADOTS(I,ITRAN),I=1,MXVFA0)
        END DO
        WRITE (LUPRI,*)
      END IF

* F^A{O} matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of F^A{O} matrix transformations:'
        DO ITRAN = 1, N1FATRAN
          WRITE(LUPRI,'(A,5I5,5X,(12I5,20X))') MSGDBG,
     &     (I1FATRAN(I,ITRAN),I=1,5),(I1FADOTS(I,ITRAN),I=1,MXVFA1)
        END DO
        WRITE (LUPRI,*)
      END IF

* F^AB{O} matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of F^AB{O} matrix transformations:'
        DO ITRAN = 1, N2FATRAN
          WRITE(LUPRI,'(A,5I5,5X,(12I5,20X))') MSGDBG,
     &     (I2FATRAN(I,ITRAN),I=1,5),(I2FADOTS(I,ITRAN),I=1,MXVFA2)
        END DO
        WRITE (LUPRI,*)
      END IF

* ETA^AB{O} vector calculations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of ETA^AB{O} vector calculations:'
        DO ITRAN = 1, N2EATRAN
          WRITE(LUPRI,'(A,2I5,5X,(12I5,20X))') MSGDBG,
     &     (I2EATRAN(I,ITRAN),I=1,2),(I2EADOTS(I,ITRAN),I=1,MXVEA2)
        END DO
        WRITE (LUPRI,*)
      END IF

      CALL FLSHFO(LUPRI)

      END IF ! (.NOT. LADD)
*---------------------------------------------------------------------*
* return
*---------------------------------------------------------------------*

      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC4R_SETUP                           *
*---------------------------------------------------------------------*

c /* deck CC_SETH1112 */
*=====================================================================*
      SUBROUTINE CC_SETH1112(IHTRAN,IHDOTS,MXTRAN,MXVEC,IZETAV,
     &                       ITAMPA,ITAMPB,ITAMPC,ITAMPD,ITRAN,IVEC)
*---------------------------------------------------------------------*
*
*    Purpose: maintain a list of dot products of H matrix 
*             transformations with 3 right amplitude vectors:  
*                       (Z*K*T^A*T^B*T^C) * T^D
*             assumes that T^A, T^B and T^C belong to the same list
*             and T^D belongs to a second list
*
*             IHTRAN - list of H matrix transformations
*             IHDOTS - list of vectors it should be dottet on
*             
*             MXTRAN - maximum list dimension
*             MXVEC  - maximum second dimesion for IHDOTS
*      
*             IZETAV - index of lagrangian multiplier vector
*             ITAMPA - index of amplitude vector A
*             ITAMPB - index of amplitude vector B
*             ITAMPC - index of amplitude vector C
*             ITAMPD - index of amplitude vector D
*
*             ITRAN - index in IHTRAN list
*             IVEC  - second index in IHDOTS list
*
*    Written by Christof Haettig, april 1997.
*
*=====================================================================*
      IMPLICIT NONE  
      INTEGER MXVEC, MXTRAN
      INTEGER IHTRAN(5,MXTRAN)
      INTEGER IHDOTS(MXVEC,MXTRAN)

      LOGICAL LFNDA, LFNDB, LFNDC, LFNDD
      INTEGER IZETAV, ITAMPA, ITAMPB, ITAMPC, ITAMPD
      INTEGER ITRAN, IVEC
      INTEGER ITAMP, I, IDX

* statement  functions:
#include <priunit.h>      
      LOGICAL LHTST, LHEND
      INTEGER IL, IA, IB,IC
      LHTST(ITRAN,IL,IA,IB,IC) = IHTRAN(1,ITRAN).EQ.IL .AND. (
     &     (IHTRAN(2,ITRAN).EQ.IA .AND. IHTRAN(3,ITRAN).EQ.IB
     &                            .AND. IHTRAN(4,ITRAN).EQ.IC) .OR. 
     &     (IHTRAN(2,ITRAN).EQ.IB .AND. IHTRAN(3,ITRAN).EQ.IA
     &                            .AND. IHTRAN(4,ITRAN).EQ.IC) .OR.
     &     (IHTRAN(2,ITRAN).EQ.IC .AND. IHTRAN(3,ITRAN).EQ.IA
     &                            .AND. IHTRAN(4,ITRAN).EQ.IB) .OR. 
     &     (IHTRAN(2,ITRAN).EQ.IA .AND. IHTRAN(3,ITRAN).EQ.IC
     &                            .AND. IHTRAN(4,ITRAN).EQ.IB) .OR.
     &     (IHTRAN(2,ITRAN).EQ.IB .AND. IHTRAN(3,ITRAN).EQ.IC
     &                            .AND. IHTRAN(4,ITRAN).EQ.IA) .OR. 
     &     (IHTRAN(2,ITRAN).EQ.IC .AND. IHTRAN(3,ITRAN).EQ.IB
     &                            .AND. IHTRAN(4,ITRAN).EQ.IA)      )
      LHEND(ITRAN) = ITRAN.GT.MXTRAN .OR.
     &               (IHTRAN(1,ITRAN)+IHTRAN(2,ITRAN)+
     &                IHTRAN(3,ITRAN)+IHTRAN(4,ITRAN) ).LE.0

*---------------------------------------------------------------------*
* set up list of K matrix transformations
*---------------------------------------------------------------------*
      ITRAN = 1
      LFNDD = LHTST(ITRAN,IZETAV,ITAMPA,ITAMPB,ITAMPC)

      DO WHILE ( .NOT. (LFNDD.OR.LHEND(ITRAN)) )
        ITRAN = ITRAN + 1
        LFNDD = LHTST(ITRAN,IZETAV,ITAMPA,ITAMPB,ITAMPC)
      END DO

      IF (.NOT.LFNDD) THEN
        IHTRAN(1,ITRAN) = IZETAV
        IHTRAN(2,ITRAN) = ITAMPA
        IHTRAN(3,ITRAN) = ITAMPB
        IHTRAN(4,ITRAN) = ITAMPC
        ITAMP = ITAMPD
      ELSE 
        ITAMP = ITAMPD
      END IF

      IVEC = 1
      DO WHILE (IHDOTS(IVEC,ITRAN).NE.ITAMP .AND.
     &           IHDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
        IVEC = IVEC + 1
      END DO

      IHDOTS(IVEC,ITRAN) = ITAMP

*---------------------------------------------------------------------*
      IF (IVEC.GT.MXVEC .OR. ITRAN.GT.MXTRAN) THEN
        WRITE (LUPRI,*) 'IVEC :',IVEC
        WRITE (LUPRI,*) 'ITRAN:',ITRAN
        WRITE (LUPRI,*) 'ITAMPA,ITAMPB,ITAMPC,ITAMPD:',
     &              ITAMPA,ITAMPB,ITAMPC,ITAMPD
        IDX = 1
        DO WHILE (.NOT. LHEND(IDX))
          WRITE(LUPRI,'(A,4I5,5X,(12I5,20X))') 'CC_SETH1112>',
     &       (IHTRAN(I,IDX),I=1,4),(IHDOTS(I,IDX),I=1,MXVEC)
          IDX = IDX + 1
        END DO
        CALL FLSHFO(LUPRI)
        CALL QUIT('Overflow error in CC_SETH1112.')
      END IF
      
      RETURN
      END 
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC_SETH1112                          *
*---------------------------------------------------------------------*
c /* deck CC_SETG212 */
*=====================================================================*
      SUBROUTINE CC_SETG212(IGTRAN,IGDOTS,MXTRAN,MXVEC,
     &                      IZETAV,ITAMPA,ITAMPB,ITAMPC,ITRAN,IVEC)
*---------------------------------------------------------------------*
*
*    Purpose: maintain a list of dot products of G matrix 
*             transformations with right amplitude vectors:  
*                       (Z*G*T^A*T^B) * T^C
*             assumes that T^A and T^C belong to one list,
*             and T^B belongs to a second list
*
*             IGTRAN - list of G matrix transformations
*             IGDOTS - list of vectors it should be dottet on
*             
*             MXTRAN - maximum list dimension
*             MXVEC  - maximum second dimesion for IGDOTS
*      
*             IZETAV - index of lagrangian multiplier vector
*             ITAMPA - index of amplitude vector A
*             ITAMPB - index of amplitude vector B
*             ITAMPC - index of amplitude vector C
*
*             ITRAN - index in IGTRAN list
*             IVEC  - second index in IGDOTS list
*
*    Written by Christof Haettig, april 1997.
*
*=====================================================================*
      IMPLICIT NONE  

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER MXVEC, MXTRAN
      INTEGER IGTRAN(4,MXTRAN)
      INTEGER IGDOTS(MXVEC,MXTRAN)

      LOGICAL LFNDC, LFNDA
      INTEGER IZETAV, ITAMPA, ITAMPB, ITAMPC
      INTEGER ITRAN, IVEC
      INTEGER ITAMP, I, IDX

#include <priunit.h>
* statement  functions:
      LOGICAL LGTST, LGEND
      INTEGER IL, IA, IB
      LGTST(ITRAN,IL,IA,IB) = IGTRAN(1,ITRAN).EQ.IL .AND.
     &     (IGTRAN(2,ITRAN).EQ.IA .AND. IGTRAN(3,ITRAN).EQ.IB)
      LGEND(ITRAN) = ITRAN.GT.MXTRAN .OR.
     &   (IGTRAN(1,ITRAN)+IGTRAN(2,ITRAN)+IGTRAN(3,ITRAN)).LE.0



       IF (LOCDBG) THEN
          WRITE (LUPRI,*) 'CC_SETG212> on entry:'
          WRITE (LUPRI,*) 'IVEC :',IVEC
          WRITE (LUPRI,*) 'ITRAN:',ITRAN
          WRITE (LUPRI,*) 'ITAMPA,ITAMPB,ITAMPC:',ITAMPA,ITAMPB,ITAMPC
          IDX = 1
          DO WHILE (.NOT. LGEND(IDX))
            WRITE(LUPRI,'(A,3I5,5X,(12I5,20X))') 'CC_SETG212>',
     &         (IGTRAN(I,IDX),I=1,3),(IGDOTS(I,IDX),I=1,MXVEC)
            IDX = IDX + 1
          END DO
          CALL FLSHFO(LUPRI)
       END IF

*---------------------------------------------------------------------*
* set up list of G matrix transformations
*---------------------------------------------------------------------*
        ITRAN = 1
        LFNDA = LGTST(ITRAN,IZETAV,ITAMPC,ITAMPB)
        LFNDC = LGTST(ITRAN,IZETAV,ITAMPA,ITAMPB)

        DO WHILE ( .NOT. (LFNDA.OR.LFNDC.OR.LGEND(ITRAN)) )
         ITRAN = ITRAN + 1
         LFNDA = LGTST(ITRAN,IZETAV,ITAMPC,ITAMPB)
         LFNDC = LGTST(ITRAN,IZETAV,ITAMPA,ITAMPB)
        END DO

        IF (.NOT.(LFNDC.OR.LFNDA)) THEN
          IGTRAN(1,ITRAN) = IZETAV
          IGTRAN(2,ITRAN) = ITAMPA
          IGTRAN(3,ITRAN) = ITAMPB
          ITAMP = ITAMPC
        ELSE 
          IF (LFNDC) ITAMP = ITAMPC
          IF (LFNDA) ITAMP = ITAMPA
        END IF

        IVEC = 1
        DO WHILE (IGDOTS(IVEC,ITRAN).NE.ITAMP .AND.
     &             IGDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
          IVEC = IVEC + 1
        END DO

        IGDOTS(IVEC,ITRAN) = ITAMP

*---------------------------------------------------------------------*
        IF (LOCDBG .OR. IVEC.GT.MXVEC .OR. ITRAN.GT.MXTRAN) THEN
          WRITE (LUPRI,*) 'CC_SETG212> on exit:'
          WRITE (LUPRI,*) 'IVEC :',IVEC
          WRITE (LUPRI,*) 'ITRAN:',ITRAN
          WRITE (LUPRI,*) 'ITAMPA,ITAMPB,ITAMPC:',ITAMPA,ITAMPB,ITAMPC
          IDX = 1
          DO WHILE (.NOT. LGEND(IDX))
            WRITE(LUPRI,'(A,3I5,5X,(12I5,20X))') 'CC_SETG212>',
     &         (IGTRAN(I,IDX),I=1,3),(IGDOTS(I,IDX),I=1,MXVEC)
            IDX = IDX + 1
          END DO
          CALL FLSHFO(LUPRI)
        END IF

        IF (IVEC.GT.MXVEC .OR. ITRAN.GT.MXTRAN) THEN
          CALL QUIT('Overflow error in CC_SETG212.')
        END IF
      
        RETURN
        END 
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC_SETG212                           *
*---------------------------------------------------------------------*
