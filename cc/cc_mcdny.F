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
*=====================================================================*
c /* deck cc_mcdcal */
*=====================================================================*
       SUBROUTINE CC_MCDCAL(WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: drive final calculation of B terms for MCD
*
*    Sonia Coriani Fall 1997
*    Restructured January 2000. Sonia
*    Calculation of LAO-B term introduced in Spring 2000. Sonia
*    New output style March 2001. Sonia
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
#include <ccroper.h>
C#include <ccropr2.h>
#include <ccmcdinf.h>
#include <ccexci.h>
#include <cclists.h>
#include <ccr1rsp.h>
#include <ccpl1rsp.h>
#include <ccel1rsp.h>
#include <ccer1rsp.h>
#include <cclrmrsp.h>
#include <cco1rsp.h>
#include <ccx1rsp.h>
#include <cco2rsp.h>
#include <dummy.h>
* arguments
      INTEGER LWORK
#if defined (SYS_CRAY)
      REAL WORK(LWORK)
#else 
      DOUBLE PRECISION WORK(LWORK)
#endif
* local parameters:
      CHARACTER*(19) MSGDBG
      PARAMETER (MSGDBG = '[debug] CC_MCDCAL> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
#if defined (SYS_CRAY)
      REAL DDOT,DNRM2
      REAL ZERO, D05, TWO
      REAL FREQEX, XLEFTC, XRIGHTC, XRIGHTAB,XLEFTAB
      REAL BTERM_A,BTERM_S, XMAB_MC, XMC_MAB
      REAL FORZADIP, SQRFORZA
#else
      DOUBLE PRECISION DDOT,DNRM2
      DOUBLE PRECISION ZERO, D05, TWO, ONE
      DOUBLE PRECISION FREQEX,  XLEFTC, XRIGHTC, XRIGHTAB,XLEFTAB
      DOUBLE PRECISION BTERM_A,BTERM_S,XMAB_MC,XMC_MAB
      DOUBLE PRECISION FORZADIP, SQRFORZA
#endif
      PARAMETER (ZERO=0.0D00, D05=0.5D00, TWO=2.0D00, ONE=1.0D00)
C      PARAMETER (DUMMY = -12345.67890)

      CHARACTER MODEL*10, MODPRI*5
      CHARACTER*8  LABELA, LABELB, LABELC, LABSOP
      CHARACTER*8  CRLXA,  CRLXB,  CRLXC
      LOGICAL LORXA, LORXB, LORXC, LPDBSA, LPDBSB, LPDBSC, LPROJ 
      INTEGER IOPERA,IOPERB,IOPERC,IOPSOP,ISGNSOP,INUM
      INTEGER ISYMA,ISYMB,ISYMC,ISYSOP
      INTEGER ITAMPA, ITAMPB, ITAMPC
      INTEGER IKAPA, IKAPB, IKAPC, ISYMAB
      INTEGER KEND0, LEND0, KEND1, LEND1
      INTEGER MXTRAN, MXVEC, ITRAN, IVEC, IORDER
      INTEGER IOPT, IOPTRD, IFREQ, ISYM, ISYML, IOPER
      INTEGER KBTERMS
      INTEGER KGTRAN, KGDOTS, KGCONS, NGTRAN, MXGTRAN, MXGDOTS
      INTEGER KFATRAN, KFADOTS, KFACONS, NFATRAN, MXFATRAN, MXFADOTS
      INTEGER KEATRAN, KEADOTS, KEACONS, NEATRAN, MXEATRAN, MXEADOTS
      INTEGER KEATRA1, KEADOT1, KEACON1, NEATRA1
      INTEGER KEATRA2, KEADOT2, KEACON2, NEATRA2
      INTEGER KFTRAN, KFDOTS, KFCONS, NFTRAN, MXFTRAN, MXFDOTS
      INTEGER KFTRA1, KFDOT1, KFCON1, NFTRA1
      INTEGER KXE2TRA, KX2DOTS, KE2DOTS, KX2CONS,KE2CONS
      INTEGER NXE2TRA, MXXETRAN,MXXEDOTS
      INTEGER KX2TRAN, KX2DOT1, KX2CON1, NX2TRAN
      INTEGER IEXCIF, ISTATE, ISYMSF, ISTATF
      INTEGER NTAMPA, NTAMPB, NTAMPC, NBCONTR

      INTEGER IZETA0,IO2AB,IM1F,IER1A,IER1B,IEO1B,IEX1B,IEL1B,IPL1A
      INTEGER IETAA,IETAB,IXKSIA,IXKSIB,IFTAMA,IFTAMB
      INTEGER IETAC,IXKSIC,IFTAMC

      INTEGER KRE_1,KRE_2,KM1F_1, KM1F_2, KTAMA_1,KTAMA_2
      INTEGER KEX1B_1,KEX1B_2,KER1B_1,KER1B_2,KER1A_1, KER1A_2
      INTEGER KPL1A_1,KPL1A_2,KEL1B_1,KEL1B_2,KFTA__1,KFTA__2
      INTEGER KEO1B_1,KEO1B_2,KFTB__1,KFTB__2
      INTEGER KLE_1,KLE_2,KO2AB_1,KO2AB_2
      INTEGER KETAA_1,KETAA_2,KETAB_1,KETAB_2
      INTEGER KXKSIA_1,KXKSIA_2
      INTEGER KXKSIC_1,KXKSIC_2,KETAC_1,KETAC_2

#if defined (SYS_CRAY)
      REAL GCON, FACON1, FACON2, BCON1, BCON2
      REAL EACON1, EACON2, EACON3, EACON4
      REAL DOTCON1, DOTCON2, DOTCON3, DOTCON4, DOTCON5
      REAL DWFB0, XE2CON1, XE2CON2, XE2CON3
#else
      DOUBLE PRECISION GCON, FACON1, FACON2, BCON1, BCON2
      DOUBLE PRECISION EACON1, EACON2, EACON3, EACON4
      DOUBLE PRECISION DOTCON1, DOTCON2, DOTCON3, DOTCON4, DOTCON5 
      DOUBLE PRECISION DWFB0, XE2CON1, XE2CON2, XE2CON3
#endif
* external functions
      INTEGER IROPER
      INTEGER IR1TAMP
      INTEGER IR1KAPPA
      INTEGER IRHSR2
      INTEGER IER1AMP
      INTEGER IEL1AMP
      INTEGER ILRMAMP
      INTEGER IETA1 
      INTEGER IRHSR1
      INTEGER IPL1ZETA

*---------------------------------------------------------------------*
* print header for magnetic circular dichroism section
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
     & '*<<<<<<<< CALCULATION OF B TERMS IN MAG',
     &                                  'NETIC CIRCULAR D. >>>>>>>>>*',
     & '*                                   ',
     &                               '                              *',
     & '************************************',
     &                               '*******************************' 

*---------------------------------------------------------------------*
* print some debug/info output
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
         WRITE(LUPRI,*) MSGDBG, ' Workspace = ',LWORK
         WRITE(LUPRI,*) MSGDBG, ' LUSEPL1 = ', LUSEPL1
      END IF
*---------------------------------------------------------------------*
* set IOPTRD for calls to CC_RDRSP 
*---------------------------------------------------------------------*
      IF (CCS) THEN
         MODEL = 'CCS       '
         IOPTRD = 1
         MODPRI = 'CCS  '
      ELSE IF (CC2) THEN 
         MODEL = 'CC2       '
         IOPTRD = 3
         MODPRI = 'CC2  '
      ELSE IF (CCSD) THEN
         MODEL = 'CCSD      '
         IOPTRD = 3
         MODPRI = 'CCSD '
      ELSE
         CALL QUIT('Unknown coupled cluster model in CC_MCD')
      END IF
*---------------------------------------------------------------------*
* allocate & initialize work space for B terms
*---------------------------------------------------------------------*
      LPROJ = .FALSE.

      NBCONTR = NMCDOPER*NMCDST 

      !max number of transformations: conservative estimate
      MXTRAN  = NLRTLBL * MAX(NLRM,NPL1LBL,NLRTLBL)
      !max number of vector products....
      MXVEC  = MAX(NLRTLBL,NER1LBL,NEL1LBL,NLRM,NO2LBL,NPL1LBL,NX1LBL)                

      MXGTRAN  = MXDIM_GTRAN  * MXTRAN
      MXFATRAN = MXDIM_FATRAN * MXTRAN
      MXEATRAN = MXDIM_XEVEC  * MXTRAN
      !for the use-left case
      MXFTRAN  = MXDIM_FTRAN  * MXTRAN
      !for the additional contrib from relax
      MXXETRAN = MXDIM_XEVEC  * MXTRAN    

      MXGDOTS  = MXVEC * MXTRAN
      MXFADOTS = MXVEC * MXTRAN
      MXEADOTS = MXVEC * MXTRAN
      !for the use-left case
      MXFDOTS  = MXVEC * MXTRAN
      !for the additional contrib from relax
      MXXEDOTS = MXVEC * MXTRAN            
   
      KBTERMS = 1
      KGTRAN  = KBTERMS + NBCONTR
      KGDOTS  = KGTRAN  + MXGTRAN
      KGCONS  = KGDOTS  + MXGDOTS
      KFATRAN = KGCONS  + MXGDOTS
      KFADOTS = KFATRAN + MXFATRAN
      KFACONS = KFADOTS + MXFADOTS
      KEATRAN = KFACONS + MXFADOTS
      KEADOTS = KEATRAN + MXEATRAN
      KEACONS = KEADOTS + MXEADOTS
      KEND0   = KEACONS + MXEADOTS
      LEND0   = LWORK   - KEND0

      !additional allocations for integer arrays ETA{A} and BMAT
      KEATRA1 = KEND0
      KEADOT1 = KEATRA1 + MXEATRAN
      KEACON1 = KEADOT1 + MXEADOTS
      KEATRA2 = KEACON1 + MXEADOTS
      KEADOT2 = KEATRA2 + MXEATRAN
      KEACON2 = KEADOT2 + MXEADOTS
      KFTRAN  = KEACON2 + MXEADOTS
      KFDOTS  = KFTRAN  + MXFTRAN
      KFCONS  = KFDOTS  + MXFDOTS
      KFTRA1  = KFCONS  + MXFDOTS
      KFDOT1  = KFTRA1  + MXFTRAN
      KFCON1  = KFDOT1  + MXFDOTS
      KEND0   = KFCON1  + MXFDOTS

      !extra contributions from PDBS
      KXE2TRA = KEND0
      KX2DOTS = KXE2TRA + MXXETRAN    !Xksi^{A^B} contrib.
      KX2CONS = KX2DOTS + MXXEDOTS
      KE2DOTS = KX2CONS + MXXEDOTS    !Eta^{A^B} contrib. 
      KE2CONS = KE2DOTS + MXXEDOTS
      KX2TRAN = KE2CONS + MXXEDOTS
      KX2DOT1 = KX2TRAN + MXXETRAN
      KX2CON1 = KX2DOT1 + MXXEDOTS
      KEND0   = KX2CON1 + MXXEDOTS
      LEND0   = LWORK   - KEND0

      IF (LEND0 .LT. 0 ) THEN
        CALL QUIT('Insufficient work space in CC_MCDCAL (1)')
      END IF

      CALL DZERO(WORK(KBTERMS),NBCONTR)

*---------------------------------------------------------------------*
* set up lists for G and F{A} transformations and ETA{A} vectors:
* as well as B (generalized F)
*---------------------------------------------------------------------*
      CALL CCMCD_SETUP(MXTRAN, MXVEC,
     &                 WORK(KGTRAN), WORK(KGDOTS), NGTRAN,
     &                 WORK(KFATRAN),WORK(KFADOTS),NFATRAN,
     &                 WORK(KEATRAN),WORK(KEADOTS),NEATRAN,
     &                 WORK(KEATRA1),WORK(KEADOT1),NEATRA1,
     &                 WORK(KEATRA2),WORK(KEADOT2),NEATRA2,
     &                 WORK(KFTRAN), WORK(KFDOTS), NFTRAN,
     &                 WORK(KFTRA1), WORK(KFDOT1), NFTRA1,
     &                 WORK(KXE2TRA),WORK(KX2DOTS),
     &                               WORK(KE2DOTS),NXE2TRA,
     &                 WORK(KX2TRAN),WORK(KX2DOT1),NX2TRAN,
     &                 WORK(KEND0), LEND0 )
*---------------------------------------------------------------------*
* calculate G matrix contributions: (G*T^A*T^B)*E^f     GCON
*---------------------------------------------------------------------*
      IOPT = 5                                  !calculate ddot product
      CALL CC_GMATRIX('L0 ','R1 ','R1 ','RE ',NGTRAN, MXVEC,
     &              WORK(KGTRAN),WORK(KGDOTS),WORK(KGCONS),
     &              WORK(KEND0), LEND0, IOPT )
*---------------------------------------------------------------------*
* calculate F{A} matrix contributions: (F^A*T^B)*RE and (F^B*T^A)*RE
*                                      FACON1            FACON2
*---------------------------------------------------------------------*
      CALL CCQR_FADRV('L0 ','o1 ','R1 ','RE ',NFATRAN, MXVEC,
     &                 WORK(KFATRAN),WORK(KFADOTS),WORK(KFACONS),
     &                 WORK(KEND0), LEND0, 'DOTP' )
*---------------------------------------------------------------------*
* calculate ETA{A} vector contributions:
*    ETAA(LE)*TB  (always)      EACON1
*    ETAA(M1)*TB,ETAB(M1)*TA    EACON2, EACON3
*    ETAB(ZA)*RE                EACON4
*---------------------------------------------------------------------*
      CALL DZERO(WORK(KEACONS),MXEADOTS)
      IOPT   = 5
      IORDER = 1
      CALL CC_XIETA( WORK(KEATRAN), NEATRAN, IOPT, IORDER, 'LE ',
     &               '---',IDUMMY,       DUMMY,
     &               'R1 ',WORK(KEADOTS),WORK(KEACONS),
     &                .FALSE.,MXVEC, WORK(KEND0), LEND0 )

      IF (LUSEPL1) THEN
         CALL DZERO(WORK(KEACON1),MXEADOTS)
         CALL DZERO(WORK(KEACON2),MXEADOTS)
         IOPT   = 5
         IORDER = 1
         CALL CC_XIETA( WORK(KEATRA1), NEATRA1, IOPT, IORDER, 'M1 ',
     &                '---',IDUMMY,       DUMMY,
     &                'R1 ',WORK(KEADOT1),WORK(KEACON1),
     &                 .FALSE.,MXVEC, WORK(KEND0), LEND0 )
         CALL FLSHFO(LUPRI)
         IOPT   = 5
         IORDER = 1
         CALL CC_XIETA( WORK(KEATRA2), NEATRA2, IOPT, IORDER, 'PL1',
     &                '---',IDUMMY,       DUMMY,
     &                'RE ',WORK(KEADOT2),WORK(KEACON2),
     &                 .FALSE.,MXVEC, WORK(KEND0), LEND0 )

      END IF

*---------------------------------------------------------------------*
* calculate B matrix contributions (use generalized F matrix): 
*           (M1*B*TA)*TB    BCON1
*           (ZA*B*TB)*RE    BCON2
*---------------------------------------------------------------------*
      IF (LUSEPL1) THEN  
        IOPT = 5
        CALL CC_FMATRIX(WORK(KFTRAN),NFTRAN,'M1 ','R1 ',IOPT,
     &                 'R1 ',WORK(KFDOTS),WORK(KFCONS),MXVEC,
     &                 WORK(KEND0), LEND0)

        IOPT = 5
        CALL CC_FMATRIX(WORK(KFTRA1),NFTRA1,'PL1','R1 ',IOPT,
     &                 'RE ',WORK(KFDOT1),WORK(KFCON1),MXVEC,
     &                 WORK(KEND0), LEND0)
      END IF
*---------------------------------------------------------------------*
* calculate PDBS extra dot product contributions:
*    M1 * Xksi{O2} and ETA{O2} * RE   (for left  moment) XECON1,XECON2
*    LE * Xksi{O2}                    (for right moment) XECON3
* It would be better to have Xksi{O2} of file and reuse it... 
*---------------------------------------------------------------------*
      CALL DZERO(WORK(KX2CONS),MXXEDOTS)
      CALL DZERO(WORK(KE2CONS),MXXEDOTS)
      CALL DZERO(WORK(KX2CON1),MXXEDOTS)

      IF (LUSEPL1) THEN    
        IOPT   = 5
        IORDER = 2
        CALL CC_XIETA( WORK(KXE2TRA), NXE2TRA, IOPT, IORDER, 'L0 ',
     &               'M1 ',WORK(KX2DOTS),WORK(KX2CONS),
     &               'RE ',WORK(KE2DOTS),WORK(KE2CONS),
     &               .FALSE.,MXVEC, WORK(KEND0), LEND0 )
        CALL CC_XIETA( WORK(KX2TRAN), NX2TRAN, IOPT, IORDER, 'L0 ',
     &               'LE ',WORK(KX2DOT1),WORK(KX2CON1),
     &               '---',IDUMMY,DUMMY,
     &               .FALSE.,MXVEC, WORK(KEND0), LEND0 )
      END IF

*---------------------------------------------------------------------*
C  Gather the different contributions and calculate remaining 
C  dot products:
*  M1  * R1      DOTCON1    M^f*T^A
*  EX1 * RE      DWFB0      DeltaW^B(0)
*  PL1 * RE      DOTCON2    
*  LE  * EO1     DWFB0      DeltaW^B(0)
*  M1  * O2      DOTCON3    M^f*Xksi^AB
C---------------------------------------------------------------------*

      DO 100  IOPER = 1, NMCDOPER

         LPROJ = .FALSE.

         IOPERA = IAMCDOP(IOPER)
         IOPERB = IBMCDOP(IOPER)
         IOPERC = ICMCDOP(IOPER)
 
         LABELA = LBLOPR(IAMCDOP(IOPER))
         LABELB = LBLOPR(IBMCDOP(IOPER))
         LABELC = LBLOPR(ICMCDOP(IOPER))
 
         ISYMA  = ISYOPR(IAMCDOP(IOPER))
         ISYMB  = ISYOPR(IBMCDOP(IOPER))
         ISYMC  = ISYOPR(ICMCDOP(IOPER))
         ISYMAB = MULD2H(ISYMA,ISYMB)
 
         IF (LOCDBG) THEN
            WRITE(LUPRI,*) MSGDBG, 'A oper: ',LABELA, ' symm: ', ISYMA
            WRITE(LUPRI,*) MSGDBG, 'B oper: ',LABELB, ' symm: ', ISYMB
            WRITE(LUPRI,*) MSGDBG, 'C oper: ',LABELC, ' symm: ', ISYMC
         END IF
 
         CALL FLSHFO(LUPRI)
  
         IF (ISYMAB.EQ.ISYMC) THEN

            LORXA  = LAMCDRX(IOPER)
            LORXB  = LBMCDRX(IOPER)
            LORXC  = LCMCDRX(IOPER)

            LPDBSA = LPDBSOP(IOPERA)
            LPDBSB = LPDBSOP(IOPERB)
            LPDBSC = LPDBSOP(IOPERC)

            IF (       LORXA) CRLXA = '(relax.)'
            IF ( .NOT. LORXA) CRLXA = '(unrel.)'
            IF (       LORXB) CRLXB = '(relax.)'
            IF ( .NOT. LORXB) CRLXB = '(unrel.)'
            IF (       LORXC) CRLXC = '(relax.)'
            IF ( .NOT. LORXC) CRLXC = '(unrel.)'


            DO ISTATE = 1, NMCDST
               ISYMSF = IMCDSTSY(ISTATE)  
               ISTATF = IMCDSTNR(ISTATE)          !RE,LE index within sym
               IEXCIF = ISYOFE(ISYMSF) + ISTATF   !RE,LE absolute index 
               FREQEX = EIGVAL(IEXCIF)
               
               IF (ISYMSF.EQ.ISYMC) THEN
 
C------------------------------------------------------------------
C           Calculate residues.
C           - request indices of vectors in the vector list
C------------------------------------------------------------------

               IZETA0 = 0
               ITAMPA = IR1TAMP(LABELA,LORXA,-FREQEX,ISYMA)   
               ITAMPB = IR1TAMP(LABELB,LORXB, ZERO,  ISYMB)     
               IKAPA  = 0
               IKAPB  = 0
               IF(LORXA)IKAPA=IR1KAPPA(LABELA,-FREQEX,ISYMA)
               IF(LORXB)IKAPB=IR1KAPPA(LABELB, ZERO,  ISYMB)

               IF (.NOT.LUSE2N1) THEN
                 ITAMPC = IR1TAMP(LABELC,LORXC, FREQEX,ISYMC)     
                 IKAPC  = 0
                 IF (LORXC) IKAPC = IR1KAPPA(LABELC, FREQEX,ISYMC)
                 IFTAMC = ITAMPC
               END IF

               IM1F   = ILRMAMP(IEXCIF,FREQEX,ISYMC) 
               IER1A  = IER1AMP(IEXCIF,FREQEX,ISYMC, 
     &                          LABELA,-FREQEX,ISYMA,.FALSE.) 
               IETAB  = IETA1(LABELB,LORXB,ZERO,ISYMB)
               IETAC  = IETA1(LABELC,LORXC,FREQEX,ISYMC)
               IFTAMB = ITAMPB
               IF (ISYMB .EQ. 1) LPROJ = .TRUE.
               IEL1B  = IEL1AMP(IEXCIF,FREQEX,ISYMC,
     &                          LABELB,ZERO,ISYMB,LORXB,LPROJ) 
               IEX1B  = IEL1B
               IXKSIA = IRHSR1(LABELA,LORXA,-FREQEX,ISYMA)
               IXKSIC = IRHSR1(LABELC,LORXC,FREQEX,ISYMC)

               IF (LUSEPL1) THEN
                  IF (ISYMB .EQ. 1) LPROJ = .TRUE.
                  IPL1A = IPL1ZETA(LABELA,LORXA,-FREQEX,ISYMA,
     &                             LPROJ,IEXCIF, FREQEX,ISYMC)
               ELSE
                  IF (ISYMB .EQ. 1) LPROJ = .TRUE.
                  IER1B  = IER1AMP(IEXCIF,FREQEX,ISYMC,
     &                             LABELB,ZERO,ISYMB,LPROJ)
                  IEO1B  = IER1B
                  IO2AB  = IRHSR2(LABELA,LORXA,-FREQEX,ISYMA,
     &                           LABELB,LORXB,   ZERO,ISYMB) 
                  IETAA  = IETA1(LABELA,LORXA,-FREQEX,ISYMA)
                  IFTAMA = ITAMPA

               END IF
*-----------------------------------------------------------------*
*           Calculate linear response residues.                   *
*-----------------------------------------------------------------*
               XLEFTC  = ZERO
               XRIGHTC = ZERO
               FORZADIP = ZERO
               SQRFORZA = ZERO

               NTAMPC = NT1AM(ISYMC) + NT2AM(ISYMC)
               IF (CCS) NTAMPC = NT1AM(ISYMC)

               KXKSIC_1 = KEND0
               KXKSIC_2 = KXKSIC_1 + NT1AM(ISYMC)
               KLE_1    = KXKSIC_2 + NT2AM(ISYMC)
               KLE_2    = KLE_1 + NT1AM(ISYMC)
               KEND1    = KLE_2 + NT2AM(ISYMC)
               LEND1    = LWORK - KEND1

               IF (LEND1.LE.0) THEN
                 WRITE(LUPRI,*)'Insuff. work in CC_MCDCAL (LR1)'
                 WRITE(LUPRI,*)'Need: ', KEND1, 'Available: ',LWORK
                 CALL QUIT('Insufficient work space in CC_MCDCAL')
               END IF

               CALL CC_RDRSP('O1 ',IXKSIC,ISYMC,IOPTRD,MODEL,
     &                       WORK(KXKSIC_1), WORK(KXKSIC_2))
               CALL CC_RDRSP('LE ',IEXCIF,ISYMC,IOPTRD,MODEL,
     &                       WORK(KLE_1), WORK(KLE_2))

               XRIGHTC = DDOT(NTAMPC,WORK(KLE_1),1,WORK(KXKSIC_1),1)
        
               IF (LUSE2N1) THEN

                 KM1F_1  = KLE_1
                 KM1F_2  = KM1F_1 + NT1AM(ISYMC)
                 KEND1   = KM1F_2 + NT2AM(ISYMC)
                 LEND1   = LWORK - KEND1

                 CALL CC_RDRSP('M1 ',IM1F,ISYMC,IOPTRD,MODEL,
     &                          WORK(KM1F_1), WORK(KM1F_2))

                 XLEFTC = DDOT(NTAMPC,WORK(KM1F_1),1,WORK(KXKSIC_1),1)

                 KETAC_1 = KM1F_1
                 KETAC_2 = KETAC_1 + NT1AM(ISYMC)
                 KRE_1   = KETAC_2 + NT2AM(ISYMC)
                 KRE_2   = KRE_1 + NT1AM(ISYMC)
                 KEND1   = KRE_2 + NT2AM(ISYMC)
                 LEND1   = LWORK - KEND1

                 CALL CC_RDRSP('RE ',IEXCIF,ISYMC,IOPTRD,MODEL,
     &                         WORK(KRE_1), WORK(KRE_2))
                 CALL CC_RDRSP('X1 ',IETAC,ISYMC,IOPTRD,MODEL,
     &                         WORK(KETAC_1), WORK(KETAC_2))
                 XLEFTC = XLEFTC + 
     &                    DDOT(NTAMPC,WORK(KETAC_1),1,WORK(KRE_1),1)
               ELSE
                 !Unfinished 
                 WRITE(LUPRI,*)'MCD: LUSE2N1=FALSE not yet implemented'
                 CALL QUIT('MCD: LUSE2N1=FALSE not yet implemented')
               END IF
c
c Dipole strength: S^of_C,C = 0.5(M^C_of M^C_fo+[M^C_of M^C_fo]^*)
c
               FORZADIP = XLEFTC*XRIGHTC
c
               IF (FORZADIP.GE.ZERO) THEN
                  SQRFORZA = SQRT(FORZADIP)
               ELSE 
                  SQRFORZA = -SQRT(-FORZADIP)
               ENDIF
C------------------------------------------------------------------
C           Calculate quadratic response residues.
C------------------------------------------------------------------
C             - LEFT PART 
C             - RIGHT PART 
C the call to *SET* routines recover ITRAN,IVEC
*----------------------------------------------------------------------*
* First contribution: {G*T^A(-w_f)*T^B(0)}*E^f (see at the beginning)
*----------------------------------------------------------------------*
              CALL CC_SETG112(WORK(KGTRAN),WORK(KGDOTS),NGTRAN,MXVEC,
     &                        IZETA0,ITAMPA,ITAMPB,IEXCIF,ITRAN,IVEC)
              GCON = WORK(KGCONS-1 + (ITRAN-1)*MXVEC + IVEC)
*----------------------------------------------------------------------*
* Second contribution: {F^A * T^B} * E^f
*----------------------------------------------------------------------*
c              CALL CC_SETFA12(WORK(KFATRAN),WORK(KFADOTS),NFATRAN,MXVEC, 
c     &                        IZETA0,IOPERA,ITAMPB,IEXCIF,ITRAN,IVEC)
              CALL CC_SETFB12(WORK(KFATRAN),WORK(KFADOTS),NFATRAN,MXVEC, 
     &                     IZETA0,IOPERA,IKAPA,ITAMPB,IEXCIF,ITRAN,IVEC)
              FACON1 = WORK(KFACONS-1 + (ITRAN-1)*MXVEC + IVEC)
*----------------------------------------------------------------------*
* Third contribution: {F^B*T^A} * E^f
*----------------------------------------------------------------------*
c              CALL CC_SETFA12(WORK(KFATRAN),WORK(KFADOTS),NFATRAN,MXVEC, 
c     &                        IZETA0,IOPERB,ITAMPA,IEXCIF,ITRAN,IVEC)
              CALL CC_SETFB12(WORK(KFATRAN),WORK(KFADOTS),NFATRAN,MXVEC, 
     &                     IZETA0,IOPERB,IKAPB,ITAMPA,IEXCIF,ITRAN,IVEC)
              FACON2 = WORK(KFACONS-1 + (ITRAN-1)*MXVEC + IVEC)
*----------------------------------------------------------------------*
* Fourth contribution M^f*{various terms} (excluding relaxation and DW)*
*----------------------------------------------------------------------*
              IF (LUSEPL1) THEN
                 CALL CC_SETXE('Eta',WORK(KEATRA1),WORK(KEADOT1),
     &                          NEATRA1,MXVEC,IM1F,IOPERA,IKAPA,0,0,0,
     &                          ITAMPB,ITRAN,IVEC)
                 EACON2 = WORK(KEACON1-1 + (ITRAN-1)*MXVEC + IVEC)
*
                 CALL CC_SETXE('Eta',WORK(KEATRA1),WORK(KEADOT1),
     &                          NEATRA1,MXVEC,IM1F,IOPERB,IKAPB,0,0,0,
     &                          ITAMPA,ITRAN,IVEC)
                 EACON3 = WORK(KEACON1-1 + (ITRAN-1)*MXVEC + IVEC)
*
                 CALL CCQR_SETF(WORK(KFTRAN),WORK(KFDOTS),NFTRAN,
     &                          MXVEC,IM1F,ITAMPA,ITAMPB,ITRAN,IVEC)
                 BCON1  = WORK(KFCONS-1 + (ITRAN-1)*MXVEC + IVEC)
*
              ELSE
                 NTAMPC = NT1AM(ISYMC) + NT2AM(ISYMC)
                 IF (CCS) NTAMPC = NT1AM(ISYMC)
                 KM1F_1  = KEND0
                 KM1F_2  = KM1F_1 + NT1AM(ISYMC)
                 KO2AB_1 = KM1F_2 + NT2AM(ISYMC)
                 KO2AB_2 = KO2AB_1 + NT1AM(ISYMC)
                 KEND1   = KO2AB_2 + NT2AM(ISYMC)
                 LEND1   = LWORK - KEND1

                 IF (LEND1.LE.0) THEN
                    WRITE(LUPRI,*)'Insuff. work sp. in  CC_MCDCAL (Mf)'
                    WRITE(LUPRI,*)'Need: ', KEND1,'Available: ', LWORK
                    CALL QUIT('Insufficient work space in CC_MCDCAL')
                 END IF

                 CALL CC_RDRSP('M1 ',IM1F,ISYMC,IOPTRD,MODEL,
     &                          WORK(KM1F_1), WORK(KM1F_2))

                 CALL CC_RDRSP('O2 ',IO2AB,ISYMC,IOPTRD,MODEL,
     &                          WORK(KO2AB_1), WORK(KO2AB_2))
                 CALL CCLR_DIASCL(WORK(KO2AB_2),D05,ISYMC)
                 DOTCON3 = DDOT(NTAMPC,WORK(KM1F_1),1,WORK(KO2AB_1),1)

              END IF 
*----------------------------------------------------------------------*
* Fifth contribution: Xksibar^A(-wf) * E^fB  (if not LUSEPL1: DOTCON4)
*----------------------------------------------------------------------*
              IF (LUSEPL1) THEN

                 CALL CC_SETF12(WORK(KFTRA1),WORK(KFDOT1),NFTRA1,MXVEC,
     &                          IPL1A,ITAMPB,IEXCIF,ITRAN,IVEC)
                 BCON2 = WORK(KFCON1-1 + (ITRAN-1)*MXVEC + IVEC)

                 CALL CC_SETXE('Eta',WORK(KEATRA2),WORK(KEADOT2),
     &                         NEATRA2,MXVEC,IPL1A,IOPERB,IKAPB,0,0,0,
     &                         IEXCIF,ITRAN,IVEC)
                 EACON4 = WORK(KEACON2-1 + (ITRAN-1)*MXVEC + IVEC)
              ELSE
                 NTAMPA = NT1AM(ISYMA) + NT2AM(ISYMA)
                 IF (CCS) NTAMPA = NT1AM(ISYMA)
                 KER1B_1 = KEND0
                 KER1B_2 = KER1B_1 + NT1AM(ISYMA)
                 KETAA_1 = KER1B_2 + NT2AM(ISYMA)
                 KETAA_2 = KETAA_1 + NT1AM(ISYMA)
                 KFTA__1 = KETAA_2 + NT2AM(ISYMA)
                 KFTA__2 = KFTA__1 + NT1AM(ISYMA)
                 KEND1   = KFTA__2 + NT2AM(ISYMA)
                 LEND1   = LWORK - KEND1

                 IF (LEND1.LE.0) THEN
                    WRITE(LUPRI,*)'Insuff. work  in  CC_MCDCAL (XibA)'
                    WRITE(LUPRI,*)'Need: ', KEND1, 'Available: ', LWORK
                    CALL QUIT('Insufficient work space in CC_MCDCAL')
                 END IF

                 CALL CC_RDRSP('ER1',IER1B,ISYMA,IOPTRD,MODEL,
     &                          WORK(KER1B_1), WORK(KER1B_2))
                 CALL CC_RDRSP('X1 ',IETAA,ISYMA,IOPTRD,MODEL,
     &                          WORK(KETAA_1), WORK(KETAA_2))
                 CALL CC_RDRSP('F1 ',IFTAMA,ISYMA,IOPTRD,MODEL,
     &                          WORK(KFTA__1), WORK(KFTA__2))
                 CALL DAXPY(NTAMPA,ONE,WORK(KFTA__1),1,WORK(KETAA_1),1)

                 DOTCON4 = DDOT(NTAMPA,WORK(KETAA_1),1,WORK(KER1B_1),1)

              END IF 
*----------------------------------------------------------------------*
* Sixth contribution Xksibar^B(-wf) * E^fA  DOTCON5
*----------------------------------------------------------------------*
              NTAMPB = NT1AM(ISYMB) + NT2AM(ISYMB)
              IF (CCS) NTAMPB = NT1AM(ISYMB)
              KER1A_1 = KEND0
              KER1A_2 = KER1A_1 + NT1AM(ISYMB)
              KETAB_1 = KER1A_2 + NT2AM(ISYMB)
              KETAB_2 = KETAB_1 + NT1AM(ISYMB)
              KFTB__1 = KETAB_2 + NT2AM(ISYMB)
              KFTB__2 = KFTB__1 + NT1AM(ISYMB)
              KEND1   = KFTB__2 + NT2AM(ISYMB)
              LEND1   = LWORK - KEND1

              IF (LEND1.LT.NT2AM(ISYMA)) THEN
                 WRITE(LUPRI,*)'Insuff. work sp. in  CC_MCDCAL (XibB*)'
                 WRITE(LUPRI,*)'Need: ', KEND1,'Available: ', LWORK
                 CALL QUIT('Insufficient work space in CC_MCDCAL')
              END IF

              CALL CC_RDRSP('ER1',IER1A,ISYMB,IOPTRD,MODEL,
     &                       WORK(KER1A_1), WORK(KER1A_2))

              CALL CC_RDRSP('X1 ',IETAB,ISYMB,IOPTRD,MODEL,
     &                       WORK(KETAB_1), WORK(KETAB_2))
              CALL CC_RDRSP('F1 ',IFTAMB,ISYMB,IOPTRD,MODEL,
     &                       WORK(KFTB__1), WORK(KFTB__2))
              CALL DAXPY(NTAMPB,ONE,WORK(KFTB__1),1,WORK(KETAB_1),1)

              DOTCON5 = DDOT(NTAMPB,WORK(KETAB_1),1,WORK(KER1A_1),1)

*----------------------------------------------------------------------* 
*   Extra contribution: {DeltaW}_of^B(0) * Mbar^f(w_f) * T^A(-w_f)
*                 and : {DeltaW}_of^B(0) * PL1A * RE
*
*   where    DeltaW^B_of(0) = Ebar^f(-wf) * Csi^fB(w_f) = 
*                           = Ebar^f(-wf) * (B*T^B + A^B)
* ---------------------------------------------------------------------* 
               IF (ISYMB.EQ.1) THEN

                 NTAMPA = NT1AM(ISYMA) + NT2AM(ISYMA)
                 IF (CCS) NTAMPA = NT1AM(ISYMA)
                 KM1F_1  = KEND0
                 KM1F_2  = KM1F_1 + NT1AM(ISYMA)
                 KTAMA_1 = KM1F_2 + NT2AM(ISYMA)
                 KTAMA_2 = KTAMA_1 + NT1AM(ISYMA)
                 KEND1   = KTAMA_2 + NT2AM(ISYMA)
                 LEND1   = LWORK - KEND1

                 IF (LEND1.LE.0) THEN
                   WRITE(LUPRI,*)'Insuff. work sp.in CC_MCDCAL (DeltaW)'
                   WRITE(LUPRI,*)'Need: ', KEND1,'Available: ', LWORK
                   CALL QUIT('Insufficient work space in CC_MCDCAL')
                 END IF

                 CALL CC_RDRSP('M1 ',IM1F,ISYMA,IOPTRD,MODEL,
     &                          WORK(KM1F_1), WORK(KM1F_2))

                 CALL CC_RDRSP('R1 ',ITAMPA,ISYMA,IOPTRD,MODEL,
     &                          WORK(KTAMA_1), WORK(KTAMA_2))
                 DOTCON1 = DDOT(NTAMPA,WORK(KM1F_1),1,WORK(KTAMA_1),1)

                 IF (LUSEPL1) THEN
                   NTAMPA = NT1AM(ISYMA) + NT2AM(ISYMA)
                   IF (CCS) NTAMPA = NT1AM(ISYMA)
                   KRE_1   = KEND0
                   KRE_2   = KRE_1 + NT1AM(ISYMA)
                   KEX1B_1 = KRE_2 + NT2AM(ISYMA)
                   KEX1B_2 = KEX1B_1 + NT1AM(ISYMA)
                   KPL1A_1 = KEX1B_2 + NT2AM(ISYMA)
                   KPL1A_2 = KPL1A_1 + NT1AM(ISYMA)
                   KEND1   = KPL1A_2 + NT2AM(ISYMA)
                   LEND1   = LWORK - KEND1
                   IF (LEND1.LE.0) THEN
                     WRITE(LUPRI,*)'Insuff. work sp.in CC_MCDCAL (LPL1)'
                     WRITE(LUPRI,*)'Need: ', KEND1,' Available: ', LWORK
                     CALL QUIT('Insufficient work space in CC_MCDCAL')
                   END IF

                   CALL CC_RDRSP('RE ',IEXCIF,ISYMA,IOPTRD,MODEL,
     &                            WORK(KRE_1), WORK(KRE_2))

                   CALL CC_RDRSP('EX1',IEX1B,ISYMA,IOPTRD,MODEL,
     &                            WORK(KEX1B_1), WORK(KEX1B_2))

                   CALL CC_RDRSP('PL1',IPL1A,ISYMA,IOPTRD,MODEL,
     &                            WORK(KPL1A_1), WORK(KPL1A_2))

                   DWFB0 = DDOT(NTAMPA,WORK(KEX1B_1),1,WORK(KRE_1),1)
                   DOTCON2 = DDOT(NTAMPA,WORK(KPL1A_1),1,WORK(KRE_1),1)
                 ELSE
                   NTAMPA = NT1AM(ISYMA) + NT2AM(ISYMA)
                   IF (CCS) NTAMPA = NT1AM(ISYMA)
                   KLE_1   = KEND0
                   KLE_2   = KRE_1 + NT1AM(ISYMA)
                   KEO1B_1 = KRE_2 + NT2AM(ISYMA)
                   KEO1B_2 = KEO1B_1 + NT1AM(ISYMA)
                   KEND1   = KEO1B_2 + NT2AM(ISYMA)
                   LEND1   = LWORK - KEND1
                   IF (LEND1.LE.0) THEN
                     WRITE(LUPRI,*)'Insuff. work sp.in CC_MCDCAL (NPL1)'
                     WRITE(LUPRI,*)'Need: ', KEND1,'Available: ', LWORK
                     CALL QUIT('Insufficient work space in CC_MCDCAL')
                   END IF

                   CALL CC_RDRSP('LE ',IEXCIF,ISYMA,IOPTRD,MODEL,
     &                            WORK(KLE_1), WORK(KLE_2))

                   CALL CC_RDRSP('EO1',IEO1B,ISYMA,IOPTRD,MODEL,
     &                            WORK(KEO1B_1), WORK(KEO1B_2))

                   DWFB0 = DDOT(NTAMPA,WORK(KLE_1),1,WORK(KEO1B_1),1)
                   DOTCON2 = ZERO
                 END IF
               ELSE
                 DWFB0   = ZERO
                 DOTCON1 = ZERO
                 DOTCON2 = ZERO
               END IF
*----------------------------------------------------------------------*
* PDBS/Relax contribution to Left 2nd order moment
*----------------------------------------------------------------------*
               IF ((LUSEPL1).AND.(LORXB.OR.LPDBSB)) THEN
                 CALL CC_FIND_SO_OP(LBLOPR(IOPERA),LBLOPR(IOPERB),
     &                      LABSOP,ISYSOP,ISGNSOP,INUM,
     &                      WORK(KEND0),LEND0)
                 IOPSOP = IROPER(LABSOP,ISYSOP)
                 CALL CC_SETXE('Xi ',WORK(KXE2TRA),WORK(KX2DOTS),
     &                          MXTRAN,MXVEC,
     &                          IZETA0,IOPSOP,IKAPA,IKAPB,0,0,IM1F,
     &                          ITRAN,IVEC)
                 XE2CON1 = WORK(KX2CONS-1 + (ITRAN-1)*MXVEC + IVEC)
                 CALL CC_SETXE('Eta',WORK(KXE2TRA),WORK(KE2DOTS),
     &                          MXTRAN,MXVEC,
     &                          IZETA0,IOPSOP,IKAPA,IKAPB,0,0,IEXCIF,
     &                          ITRAN,IVEC)
                 XE2CON2 = WORK(KE2CONS-1 + (ITRAN-1)*MXVEC + IVEC)
                 IF (LOCDBG) THEN
                       WRITE(LUPRI,*) MSGDBG,
     &                'XE2CON1: ',XE2CON1,' XE2CON2: ',XE2CON2
                 END IF
               ELSE
                 XE2CON1 = ZERO
                 XE2CON2 = ZERO
               END IF

*----------------------------------------------------------------------*
* TOTAL LEFT TRANSITION M^AB{n<-f} * M^C{f<-n} = XMAB_MC
*----------------------------------------------------------------------*
              XLEFTAB = ZERO
              XMAB_MC = ZERO
c
              XLEFTAB = GCON + FACON1 + FACON2 + DOTCON5 +
     &                  DWFB0*(DOTCON1-DOTCON2) +
     &                  XE2CON1 + XE2CON2
              IF (LUSEPL1) THEN
                 XLEFTAB = XLEFTAB + EACON2 + EACON3 + BCON1 +
     &                     EACON4 + BCON2
              ELSE
                 XLEFTAB = XLEFTAB + DOTCON3 + DOTCON4
              END IF

              XMAB_MC  = XLEFTAB * XRIGHTC 


************************************************************************
*                     RIGHT PART                                       *
************************************************************************
*----------------------------------------------------------------------*
* First contribution: Ebar^fB(-w_f,0) * Xksi^A                         *
*----------------------------------------------------------------------*
               
              NTAMPA   = NT1AM(ISYMA) + NT2AM(ISYMA)
              KEL1B_1  = KEND0
              KEL1B_2  = KEL1B_1 + NT1AM(ISYMA)
              KXKSIA_1 = KEL1B_2 + NT2AM(ISYMA)
              KXKSIA_2 = KXKSIA_1 + NT1AM(ISYMA) 
              KEND1    = KXKSIA_2 + NT2AM(ISYMA)
              LEND1    = LWORK - KEND1

              CALL CC_RDRSP('EL1',IEL1B,ISYMA,IOPTRD,MODEL,
     &                       WORK(KEL1B_1),WORK(KEL1B_2))
              CALL CC_RDRSP('O1 ',IXKSIA,ISYMA,IOPTRD,MODEL,
     &                       WORK(KXKSIA_1),WORK(KXKSIA_2))

C             CALL CC_XKSI(WORK(KXKSIA_1),LABELA,ISYMA,DUMMY,
C    *                     WORK(KEND1),LEND1)

              DOTCON5 = DDOT(NTAMPA,WORK(KXKSIA_1),1,WORK(KEL1B_1),1)
* ---------------------------------------------------------------------*
*  Second contribution: Ebar^f * A^A * T^B(0) 
*                       computed as Eta^A(EL1) * T^B(0)
* ---------------------------------------------------------------------*
             CALL CC_SETXE('Eta',WORK(KEATRAN),WORK(KEADOTS),NEATRAN,
     &                      MXVEC,IEXCIF,IOPERA,IKAPA,0,0,0,ITAMPB,
     &                      ITRAN,IVEC)
             EACON1 = WORK(KEACONS-1 + (ITRAN-1)*MXVEC + IVEC)
*----------------------------------------------------------------------*
* PDBS/Relax contribution to Right 2nd order moment
*----------------------------------------------------------------------*
               IF ((LUSEPL1).AND.(LORXB.OR.LPDBSB)) THEN
c                 CALL CC_FIND_SO_OP(LBLOPR(IOPERA),LBLOPR(IOPERB),
c     &                      LABSOP,ISYSOP,ISGNSOP,INUM,
c     &                      WORK(KEND0),LEND0)
                 CALL CC_SETXE('Xi ',WORK(KX2TRAN),WORK(KX2DOT1),
     &                          MXTRAN,MXVEC,
     &                          IZETA0,IOPSOP,IKAPA,IKAPB,0,0,IEXCIF,
     &                          ITRAN,IVEC)
                 XE2CON3 = WORK(KX2CON1-1 + (ITRAN-1)*MXVEC + IVEC)
                 IF (LOCDBG) THEN
                    WRITE(LUPRI,*) MSGDBG,'XE2CON3: ',XE2CON3
                 END IF
               ELSE
                 XE2CON3 = ZERO
               END IF
*
*----------------------------------------------------------------------*
* TOTAL RIGHT TRANSITION M^C{o<-f} * M^AB{f<-o} = XMC_MAB
*----------------------------------------------------------------------*
*
              XRIGHTAB = ZERO
              XMC_MAB  = ZERO

              XRIGHTAB = DOTCON5 + EACON1 + XE2CON3
              XMC_MAB  = XLEFTC*XRIGHTAB  
*----------------------------------------------------------------------*
* Debug information
*----------------------------------------------------------------------*
              IF (LOCDBG) THEN
                 WRITE(LUPRI,'(/1X,A,A,2(/1X,A,F16.8))')
     *           'For C operator ', LABELC,
     *           ' Left  1st single moment = ', XLEFTC,
     *           ' Right 1st single moment = ', XRIGHTC
                 WRITE(LUPRI,'(1X,A,F16.8,A1,F16.8,A1)') 
     *           ' Dipole Strength (a.u.)  = ',FORZADIP,'(',SQRFORZA,')'
                 WRITE(LUPRI,'(/1X,A,A,A,A,3(/1X,A,F16.8))')
     *           'For A, B operators ', LABELA,', ', LABELB,
     *           ' Left  2nd residue AB    = ', XLEFTAB,
     *           ' Right 2nd residue AB    = ', XRIGHTAB,
     *           ' S^of_AB,AB              = ', XLEFTAB*XRIGHTAB
                 CALL FLSHFO(LUPRI)
              END IF
*
* ----------------------------------------------------------------------*
*  combine various terms to get B term contribution
* ----------------------------------------------------------------------*
*
              BTERM_S = D05 * (XMAB_MC + XMC_MAB)
              BTERM_A = D05 * (XMAB_MC - XMC_MAB)
*
* ----------------------------------------------------------------------*
*  Write output
*  Note: final result divided by 2 because of (mu = -0.5 L)
* ----------------------------------------------------------------------*
*
         WRITE(LUPRI,'(1X,58("-"))')
         WRITE(LUPRI,'(/1x,a,f9.5,a,i1)') 
     &  'For transition |o> -> |f(',FREQEX,')>, of symm. ', ISYMC 
         WRITE(LUPRI,'(3(/1x,a,a,a,a1,f9.5,a,i1))') 
     &     ' A oper.: ',LABELA,CRLXA,'(',-FREQEX,'), symm. ',ISYMA,
     &     ' B oper.: ',LABELB,CRLXB,'(',ZERO   ,'), symm. ',ISYMB,
     &     ' C oper.: ',LABELC,CRLXC,'(',FREQEX ,'), symm. ',ISYMC
C         WRITE(LUPRI,'(3(/1X,A,F16.8))') 
C     &      ' S^of_C,C  = (Dipole strength (au))      = ', FORZADIP,
C     &      ' S^of_AB,C = M^AB_of(0.0) x M^C_fo   = ', XMAB_MC,
C     &      ' S^of_C,AB = M^C*_of x M^AB*_fo(0.0) = ', XMC_MAB
         WRITE(LUPRI,'(3(/1X,A,F16.8))') 
     &     ' Dipole strength (au) (C oper) = ',FORZADIP,
     &     ' M^AB_of(0.0) x M^C_fo         = ',XMAB_MC,
     &     ' M^C*_of x M^AB*_fo(0.0)       = ',XMC_MAB
         IF (IPRINT.GT.5) THEN 
            WRITE(LUPRI,'(/1X,a5,A,F16.6,A,F16.8,A)')
     &     MODPRI,
     &     ' B term contribution (au): ', D05*BTERM_A, ' (antisym) ', 
     &                                    D05*BTERM_S, ' (symm)    '
            WRITE(LUPRI,'(1X,69("-"))')
         ELSE
            WRITE(LUPRI,'(/1X,A5,A,F16.8,A)')
     &     MODPRI,
     &     ' B term contribution (au): ', D05*BTERM_A, ' (antisym) ' 
            WRITE(LUPRI,'(1X,58("-"))')
         END IF
* ---------------------------------------------------------------------*
           END IF
           END DO 

         ENDIF
 100  CONTINUE

      RETURN
      END
*=====================================================================*
*              END OF SUBROUTINE CC_MCDCAL                            *
*=====================================================================*
