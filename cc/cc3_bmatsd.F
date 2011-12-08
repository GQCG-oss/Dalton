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
C    /* Deck cc3_bmatsd */ 
      SUBROUTINE CC3_BMATSD(OMEGA1,OMEGA2,OMEGA1EFF,OMEGA2EFF,ISYRES,
     *                       LISTB,IDLSTB,LISTC,IDLSTC,
     *                       WORK,LWORK)
*
*******************************************************************
*
* This routine calculates the part of the right hand side vector
* for the paritioned second-order amplitude equations.
*
* Calculate the B matrix contributions (using W intermmediate): 
*
* W^BD = 1/2<mu3|[H^BC,T2^0]|HF> + <mu3|[H^B,T2^C]|HF>
*
* projected up into the Singles-and-Doubles space:
*
* omega1eff = omega1eff + <mu1|[H,W^BD]|HF>
* omega2eff = omega2eff + <mu2|[H,W^BD]|HF>
*
*
* At the end put omega to omegaeff:
*
* omega1eff = omega1eff + omega1
* omega2eff = omega2eff + omega2
*
*
* F. Pawlowski, P. Jorgensen, Aarhus Spring 2003.
*
*******************************************************************
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccr1rsp.h>
#include <ccinftap.h>
#include <iratdef.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <dummy.h>
#include <ccorb.h>
#include <ccexci.h>
C
      CHARACTER*6 FN3VI
      CHARACTER*8 FNTOC,FN3VI2
C        
      PARAMETER (FNTOC   = 'CCSDT_OC' )
      PARAMETER (FN3VI   = 'CC3_VI'  , FN3VI2  = 'CC3_VI12')

      CHARACTER*12 FN3SRTR, FNCKJDR, FNDELDR, FNDKBCR
      CHARACTER*13 FNCKJDR1, FNDKBCR1
      CHARACTER*13 FNCKJDR2, FNDKBCR2
C
      PARAMETER(FN3SRTR  = 'CCSDT_FBMAT1', FNCKJDR  = 'CCSDT_FBMAT2',
     *          FNDELDR  = 'CCSDT_FBMAT3', FNDKBCR  = 'CCSDT_FBMAT4')
C
      PARAMETER(FNCKJDR1  = 'CCSDT_FBMAT21',FNDKBCR1  = 'CCSDT_FBMAT41')
      PARAMETER(FNCKJDR2  = 'CCSDT_FBMAT22',FNDKBCR2  = 'CCSDT_FBMAT42')
C
      INTEGER LUTOC,LU3VI,LU3VI2
      INTEGER LU3SRTR, LUCKJDR, LUDELDR, LUDKBCR
      INTEGER LUCKJDR1, LUDKBCR1
      INTEGER LUCKJDR2, LUDKBCR2
C
      CHARACTER*3 LISTB, LISTC
      CHARACTER*8 LABELB, LABELC
      INTEGER IDLSTB,IDLSTC
      INTEGER ISYMRB,ISYMRC
C
      CHARACTER CDUMMY*1
C
      INTEGER ISYRES
      INTEGER LWORK
C
      INTEGER ISYM0,ISINT1,ISINT2,ISYMT2
      INTEGER KT2TP,KFOCKD,KLAMDP,KLAMDH,KXIAJB,KFOCK0,KEND00,LWRK00
      INTEGER KT1AM,KEND01,LWRK01
      INTEGER IOPTTCME,IOPT
      INTEGER KT2C,KT1B,KEND1,LWRK1
      INTEGER ISINTR1B,ISINTR2B
      INTEGER ISINTBC
      INTEGER KRBJIA,KTROC,KTROC1,KTROC0Y,KTROC3,KINTOC
      INTEGER MAXOCC,KEND2,LWRK2
      INTEGER IOFF
      INTEGER ISYMD,ISCKB1,ISCKB2Y,ISCKBBC
      INTEGER KTRVI,KTRVI1,KTRVI0Y,KEND3,LWRK3,KTRVI7,KINTVI
      INTEGER ISYMB,ISYMBC,ISYMBD,ISWMAT,ISCKD2Y,ISCKDBC
      INTEGER KDIAGW,KWMAT,KTMAT,KEND4,KTRVI8Y,LWRK4,KTRVI9
      INTEGER LENSQW,KINDSQW
      INTEGER KEND0,LWRK0
      INTEGER ISINTR1C,ISINTR2C
      INTEGER KT2B,KT1C,KTROC0X
      INTEGER ISCKB2X,KTRVI0X
      INTEGER ISCKD2X,KTRVI8X
      INTEGER MMAXOCC
      INTEGER ILSTSYM
c
      integer kx3am
c
C
#if defined (SYS_CRAY)
      REAL OMEGA1(*),OMEGA1EFF(*),OMEGA2(*),OMEGA2EFF(*)
      REAL WORK(LWORK)
      REAL FREQB,FREQC,FREQBC
      REAL XNORMVAL,DDOT
#else
      DOUBLE PRECISION OMEGA1(*),OMEGA1EFF(*),OMEGA2(*),OMEGA2EFF(*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION FREQB,FREQC,FREQBC
      DOUBLE PRECISION XNORMVAL,DDOT
#endif
C

      CALL QENTER('CC3_BMATSD')
C
C----------------------
C     Lists check
C----------------------
C
      IF ( (LISTB(1:3).EQ.'R1 ') .AND. (LISTC(1:3).EQ.'R1 ') ) THEN
         CONTINUE
      ELSE IF ( (LISTB(1:3).EQ.'R1 ') .AND. (LISTC(1:3).EQ.'RE ') ) THEN
         CONTINUE
      ELSE
         WRITE(LUPRI,*)'LISTB = ', LISTB
         WRITE(LUPRI,*)'LISTC = ', LISTC
         WRITE(LUPRI,*)'Only implemented when '
     *   //'LISTB = R1 and LISTC = R1 or RE '
         CALL QUIT('Case not implemented in CC3_BMATSD')
      END IF
C
C--------------------
C     Lists handling
C--------------------
C
      IF (LISTB(1:3).EQ.'R1 ') THEN 
         ISYMRB = ISYLRT(IDLSTB)
         FREQB  = FRQLRT(IDLSTB)
         LABELB = LRTLBL(IDLSTB)
      ELSE
         CALL QUIT('Unknown list for LISTB in CC3_BMATSD')
      END IF
C
      IF (LISTC(1:3).EQ.'R1 ') THEN
         ISYMRC = ISYLRT(IDLSTC)
         FREQC  = FRQLRT(IDLSTC)
         LABELC = LRTLBL(IDLSTC)
      ELSE IF (LISTC(1:3).EQ.'RE ') THEN
         ISYMRC = ILSTSYM(LISTC,IDLSTC)
         FREQC  = EIGVAL(IDLSTC)
         LABELC = '- none -'
      ELSE
         CALL QUIT('Unknown list for LISTC in CC3_BMATSD')
      END IF
C
C-------------------------------------------
C     Construct FREQBC = (omega_B + omega_C)
C-------------------------------------------
C
      FREQBC = FREQB + FREQC
C
C--------------------------
C     Save and set flags.
C--------------------------
C
      ISYM0   = 1
      ISINT1  = 1
      ISINT2  = 1
      ISYMT2  = 1
C
C------------------------------------------
C     Open files (integrals in contraction)
C------------------------------------------
C
      LUTOC    = -1
      LU3VI    = -1
      LU3VI2   = -1
C
      CALL WOPEN2(LUTOC,FNTOC,64,0)
      CALL WOPEN2(LU3VI,FN3VI,64,0)
      CALL WOPEN2(LU3VI2,FN3VI2,64,0)
C
C-----------------------------------------
C     Open temporary files (H^B integrals)
C-----------------------------------------
C
      LU3SRTR  = -1
      LUCKJDR  = -1
      LUDELDR  = -1
      LUDKBCR  = -1
C
      CALL WOPEN2(LU3SRTR,FN3SRTR,64,0)
      CALL WOPEN2(LUCKJDR,FNCKJDR,64,0)
      CALL WOPEN2(LUDELDR,FNDELDR,64,0)
      CALL WOPEN2(LUDKBCR,FNDKBCR,64,0)
C
C------------------------------------------
C     Open temporary files (H^BC integrals)
C------------------------------------------
C
      LUCKJDR1  = -1
      LUDKBCR1  = -1
C
      CALL WOPEN2(LUCKJDR1,FNCKJDR1,64,0)
      CALL WOPEN2(LUDKBCR1,FNDKBCR1,64,0)
C
C------------------------------------------
C     Open temporary files (H^C integrals)
C------------------------------------------
C
      LUCKJDR2  = -1
      LUDKBCR2  = -1
C
      CALL WOPEN2(LUCKJDR2,FNCKJDR2,64,0)
      CALL WOPEN2(LUDKBCR2,FNDKBCR2,64,0)
C
C----------------------------------------------
C     Calculate the zeroth order stuff once
C----------------------------------------------
C
      KT2TP  = 1
      KFOCKD = KT2TP  + NT2SQ(ISYMT2)
      KLAMDP = KFOCKD + NORBTS
      KLAMDH = KLAMDP + NLAMDT
      KXIAJB = KLAMDH + NLAMDT
      KFOCK0 = KXIAJB + NT2AM(ISINT1)
      KEND00 = KFOCK0 + N2BST(ISYM0)
      LWRK00 = LWORK  - KEND00
C
      KT1AM  = KEND00
      KEND01 = KT1AM  + NT1AM(ISYMT2)
      LWRK01 = LWORK  - KEND01
C
      IF (LWRK01 .LT. 0) THEN
         WRITE(LUPRI,*)'Memory available: ',LWORK
         WRITE(LUPRI,*)'Memory needed: ',KEND01
         CALL QUIT('Out of memory in CC3_BMATSD (zeroth allo.)')
      ENDIF
C
C------------------------
C     Construct L(ia,jb).
C------------------------
C
      REWIND(LUIAJB)
      CALL READI(LUIAJB,IRAT*NT2AM(ISINT1),WORK(KXIAJB))
      IOPTTCME = 1
      CALL CCSD_TCMEPK(WORK(KXIAJB),1.0D0,ISINT1,IOPTTCME)
C
      IF ( IPRINT .GT. 55) THEN
         XNORMVAL = DDOT(NT2AM(ISINT1),WORK(KXIAJB),1,
     *                WORK(KXIAJB),1)
         WRITE(LUPRI,*) 'Norm of IAJB ',XNORMVAL
      ENDIF
C
C----------------------------------------------------------------
C     Read t1 and calculate the zero'th order Lambda matrices
C----------------------------------------------------------------
C
      CALL GET_LAMBDA0(WORK(KLAMDP),WORK(KLAMDH),WORK(KEND01),LWRK01)
C
C-------------------------------------------
C     Read in t2 , square it and reorder 
C-------------------------------------------
C
      IOPT = 2
      CALL GET_T1_T2(IOPT,.FALSE.,DUMMY,WORK(KT2TP),'R0',0,ISYMT2,
     *               WORK(KEND01),LWRK01)

      IF (IPRINT .GT. 55) THEN
         XNORMVAL = DDOT(NT2SQ(ISYMT2),WORK(KT2TP),1,WORK(KT2TP),1)
         WRITE(LUPRI,*) 'Norm of T2TP ',XNORMVAL
      ENDIF
C
C------------------------------------------------------------------
C     Read in Fock matrix in AO basis (from file) and transform to    
C     Lambda_0 basis.             
C------------------------------------------------------------------
C
      CALL GET_FOCK0(WORK(KFOCK0),WORK(KLAMDP),WORK(KLAMDH),
     *               WORK(KEND01),LWRK01)
C
C--------------------------------------
C     Read in orbital energies
C--------------------------------------
C
      CALL GET_ORBEN(WORK(KFOCKD),WORK(KEND01),LWRK01)
c
c If we want to sum the T3 amplitudes
c
      if (.false.) then
         kx3am  = kend00
         kend00 = kx3am + nt1amx*nt1amx*nt1amx
         call dzero(work(kx3am),nt1amx*nt1amx*nt1amx)
         lwrk00 = lwork - kend00
         if (lwrk00 .lt. 0) then
            write(lupri,*) 'Memory available : ',lwork
            write(lupri,*) 'Memory needed    : ',kend00
            call quit('Insufficient space (T3) in CC3_BMATSD (1)')
         END IF
      endif
C
C-------------------------------------------------
C     Read T1^B and T2^C
C     Calculate (ck|de)-Btilde  and (ck|lm)-Btilde
C     Calculate (ck|de)-BCtilde and (ck|lm)-BCtilde
C     Used to construct WBD intermmediate
C-------------------------------------------------
C
      KT2C   = KEND00
      KEND0  = KT2C   + NT2SQ(ISYMRC)
      LWRK0  = LWORK  - KEND0
C
      KT2B   = KEND0
      KEND0  = KT2B   + NT2SQ(ISYMRB)
      LWRK0  = LWORK  - KEND0
C
      KT1C   = KEND0
      KEND0  = KT1C   + NT1AM(ISYMRC)
      LWRK0  = LWORK  - KEND0
C
      KT1B   = KEND0
      KEND1  = KT1B   + NT1AM(ISYMRB)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. NT2SQ(ISYMRB)) THEN
         CALL QUIT('Out of memory in CC3_BMATSD (TOT_T3Y) ')
      ENDIF
C
C--------------------------
C     Read in T1^B and T2^B
C--------------------------
C
      IOPT = 3
      CALL GET_T1_T2(IOPT,.TRUE.,WORK(KT1B),WORK(KT2B),LISTB,IDLSTB,
     *               ISYMRB,WORK(KEND1),LWRK1)
C
      IF (IPRINT .GT. 55) THEN
         XNORMVAL = DDOT(NT1AM(ISYMRB),WORK(KT1B),1,WORK(KT1B),1)
         WRITE(LUPRI,*) 'Norm of T1B  ',XNORMVAL
         XNORMVAL = DDOT(NT2SQ(ISYMRB),WORK(KT2B),1,WORK(KT2B),1)
         WRITE(LUPRI,*) 'Norm of T2B  ',XNORMVAL
      ENDIF
C
C--------------------------
C     Read in T1^C and T2^C
C--------------------------
C
      IOPT = 3
      CALL GET_T1_T2(IOPT,.TRUE.,WORK(KT1C),WORK(KT2C),LISTC,IDLSTC,
     *               ISYMRC,WORK(KEND1),LWRK1)
C
      IF (IPRINT .GT. 55) THEN
         XNORMVAL = DDOT(NT1AM(ISYMRC),WORK(KT1C),1,WORK(KT1C),1)
         WRITE(LUPRI,*) 'Norm of T1C  ',XNORMVAL
         XNORMVAL = DDOT(NT2SQ(ISYMRC),WORK(KT2C),1,WORK(KT2C),1)
         WRITE(LUPRI,*) 'Norm of T2C  ',XNORMVAL
      ENDIF
C
C----------------------------------------------------
C     Integrals (ck|de)-tilde(B) and (ck|lm)-tilde(B)
C----------------------------------------------------
C
      ISINTR1B = MULD2H(ISINT1,ISYMRB)
      ISINTR2B = MULD2H(ISINT2,ISYMRB)
C
      CALL CC3_BARINT(WORK(KT1B),ISYMRB,WORK(KLAMDP),
     *                WORK(KLAMDH),WORK(KEND1),LWRK1,
     *                LU3SRTR,FN3SRTR,LUCKJDR,FNCKJDR)
C
      CALL CC3_SORT1(WORK(KEND1),LWRK1,2,ISINTR1B,LU3SRTR,FN3SRTR,
     *               LUDELDR,FNDELDR,IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *               IDUMMY,CDUMMY)
C
      CALL CC3_SINT(WORK(KLAMDH),WORK(KEND1),LWRK1,ISINTR1B,
     *              LUDELDR,FNDELDR,LUDKBCR,FNDKBCR)
C
C--------------------------------
C    Re-use some temporary files
C--------------------------------
C
      CALL WCLOSE2(LU3SRTR,FN3SRTR,'DELETE')
      CALL WCLOSE2(LUDELDR,FNDELDR,'DELETE')
C
      CALL WOPEN2(LU3SRTR,FN3SRTR,64,0)
      CALL WOPEN2(LUDELDR,FNDELDR,64,0)
C
C------------------------------------------------------
C     Calculate the (ck|de)-BCtilde and (ck|lm)-BCtilde
C------------------------------------------------------
C
      ISYMBC   = MULD2H(ISYMRB,ISYMRC)
C
      CALL CC3_3BARINT(ISYMRB,LISTB,IDLSTB,ISYMRC,LISTC,IDLSTC,
     *                 IDUMMY,CDUMMY,IDUMMY,.FALSE.,
     *                 WORK(KLAMDP),WORK(KLAMDH),WORK(KEND1),LWRK1,
     *                 LU3SRTR,FN3SRTR,LUCKJDR1,FNCKJDR1)
C
      CALL CC3_SORT1(WORK(KEND1),LWRK1,2,ISYMBC,LU3SRTR,FN3SRTR,
     *               LUDELDR,FNDELDR,IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *               IDUMMY,CDUMMY)
C
      CALL CC3_SINT(WORK(KLAMDH),WORK(KEND1),LWRK1,ISYMBC,
     *              LUDELDR,FNDELDR,LUDKBCR1,FNDKBCR1)
C
C--------------------------------
C    Re-use some temporary files
C--------------------------------
C
      CALL WCLOSE2(LU3SRTR,FN3SRTR,'DELETE')
      CALL WCLOSE2(LUDELDR,FNDELDR,'DELETE')
C
      CALL WOPEN2(LU3SRTR,FN3SRTR,64,0)
      CALL WOPEN2(LUDELDR,FNDELDR,64,0)
C
C----------------------------------------------------
C     Integrals (ck|de)-tilde(C) and (ck|lm)-tilde(C)
C----------------------------------------------------
C
      ISINTR1C = MULD2H(ISINT1,ISYMRC)
      ISINTR2C = MULD2H(ISINT2,ISYMRC)
C
      CALL CC3_BARINT(WORK(KT1C),ISYMRC,WORK(KLAMDP),
     *                WORK(KLAMDH),WORK(KEND1),LWRK1,
     *                LU3SRTR,FN3SRTR,LUCKJDR2,FNCKJDR2)
C
      CALL CC3_SORT1(WORK(KEND1),LWRK1,2,ISINTR1C,LU3SRTR,FN3SRTR,
     *               LUDELDR,FNDELDR,IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *               IDUMMY,CDUMMY)
C
      CALL CC3_SINT(WORK(KLAMDH),WORK(KEND1),LWRK1,ISINTR1C,
     *              LUDELDR,FNDELDR,LUDKBCR2,FNDKBCR2)
C
C--------------------------------
C        Read occupied integrals
C-------------------------------
C
      ISINTBC = MULD2H(ISINT2,ISYMBC)
C
      !Use KEND0, because KT1B is not needed any more
      KRBJIA = KEND0
      KTROC  = KRBJIA + NT2SQ(ISYRES)
      KTROC1 = KTROC  + NTRAOC(ISINT1)   !KTROC - int. in contraction
      KEND1  = KTROC1 + NTRAOC(ISINT1)   !KTROC1 - int. in contraction
      LWRK1  = LWORK  - KEND1
C
      KTROC0Y = KEND1                   
      KEND1   = KTROC0Y + NTRAOC(ISINTR2B)!KTROC0Y - int. in <mu3|[H^B,T2^C]|HF>
C
      KTROC0X = KEND1                   
      KEND1   = KTROC0X + NTRAOC(ISINTR2C)!KTROC0X - int. in <mu3|[H^C,T2^B]|HF>
C                                                                 ===
      KTROC3 = KEND1
      KEND1  = KTROC3   + NTRAOC(ISINTBC)!KTROC3  - int. in <mu3|[H^BC,T2]|HF>
C                                                                 ===
      KINTOC  = KEND1
      MAXOCC  = MAX(NTOTOC(ISINTR2B),NTOTOC(ISINTBC))
      MMAXOCC = MAX(NTOTOC(ISINTR2C),MAXOCC)
      KEND2   = KINTOC  + MAX(NTOTOC(ISINT1),MMAXOCC)!KINTOC - temporary storage
      LWRK2   = LWORK   - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in CC3_BMATSD (2)')
      END IF
C
      CALL DZERO(WORK(KRBJIA),NT2SQ(ISYRES))
C
C
C-------------------------------------------------------------
C     B-transformed occupied integrals for <mu3|[H^B,T2^C]|HF>
C-------------------------------------------------------------
C
         IOFF = 1
         IF (NTOTOC(ISINTR2B) .GT. 0) THEN
            CALL GETWA2(LUCKJDR,FNCKJDR,WORK(KINTOC),IOFF,
     *                  NTOTOC(ISINTR2B))
         ENDIF
C
         IF (IPRINT .GT. 55) THEN
            XNORMVAL = DDOT(NTOTOC(ISINTR2B),WORK(KINTOC),1,
     *                   WORK(KINTOC),1)
            WRITE(LUPRI,*) 'Norm of CCSDT_OC-INT (B transformed) ',
     *                      XNORMVAL
         ENDIF
C
         CALL CC3_TROCC(WORK(KINTOC),WORK(KTROC0Y),WORK(KLAMDP),
     *                  WORK(KEND2),LWRK2,ISINTR2B)
C
C-------------------------------------------------------------
C     C-transformed occupied integrals for <mu3|[H^C,T2^B]|HF>
C-------------------------------------------------------------
C
         IOFF = 1
         IF (NTOTOC(ISINTR2C) .GT. 0) THEN
            CALL GETWA2(LUCKJDR2,FNCKJDR2,WORK(KINTOC),IOFF,
     *                  NTOTOC(ISINTR2C))
         ENDIF
C
         IF (IPRINT .GT. 55) THEN
            XNORMVAL = DDOT(NTOTOC(ISINTR2C),WORK(KINTOC),1,
     *                   WORK(KINTOC),1)
            WRITE(LUPRI,*) 'Norm of CCSDT_OC-INT (C transformed) ',
     *                      XNORMVAL
         ENDIF
C
         CALL CC3_TROCC(WORK(KINTOC),WORK(KTROC0X),WORK(KLAMDP),
     *                  WORK(KEND2),LWRK2,ISINTR2C)
C
C--------------------------------------------------------------------------
C    BC-transformed occupied integrals for <mu3|[H^BC,T2]|HF>
C--------------------------------------------------------------------------
C
      IOFF = 1 
      IF (NTOTOC(ISINTBC) .GT. 0) THEN
         CALL GETWA2(LUCKJDR1,FNCKJDR1,WORK(KINTOC),IOFF,
     *               NTOTOC(ISINTBC))
      ENDIF
C
      CALL CC3_TROCC(WORK(KINTOC),WORK(KTROC3),WORK(KLAMDP),
     *                    WORK(KEND2),LWRK2,ISINTBC)
C
C-----------------------------------
C   Occupied integrals in contraction
C-----------------------------------
C
      IOFF = 1
      IF (NTOTOC(ISINT1) .GT. 0) THEN
         CALL GETWA2(LUTOC,FNTOC,WORK(KINTOC),IOFF,NTOTOC(ISINT1))
      ENDIF
C
      !Write out norms of integrals.
      IF (IPRINT .GT. 55) THEN
         XNORMVAL  = DDOT(NTOTOC(ISINT1),WORK(KINTOC),1,
     *                WORK(KINTOC),1)
         WRITE(LUPRI,*) 'Norm of CCSDT_OC-INT ',XNORMVAL
      ENDIF
C
      !Transform (ia|j delta) integrals to (ia|j k) and sort as (i,j,k,a)
      CALL CCSDT_TROCC(WORK(KINTOC),WORK(KTROC),WORK(KLAMDH),
     *                  WORK(KEND2),LWRK2)
C
      !sort (i,j,k,a) as (a,i,j,k)
      CALL CCSDT_SRTOC2(WORK(KTROC),WORK(KTROC1),ISINT1,
     *                  WORK(KEND2),LWRK2)
C
C---------------------
C     Start ISYMD-loop
C---------------------
C
      DO ISYMD = 1,NSYM
C
         ISCKB1  = MULD2H(ISINT1,ISYMD)      
         ISCKB2Y = MULD2H(ISINTR2B,ISYMD)
         ISCKB2X = MULD2H(ISINTR2C,ISYMD)
         ISCKBBC = MULD2H(ISINTBC,ISYMD)
C
C----------------------------------------
C        Read virtual integrals (fixed D)
C----------------------------------------
C
         !Use KEND1, because KINTOC is not needed any more
         KTRVI  = KEND1 
         KTRVI1 = KTRVI  + NCKATR(ISCKB1)   !KTRVI - int. in contraction
         KEND2 = KTRVI1 + NCKATR(ISCKB1)    !KTRVI1 - int. in contraction
         LWRK2  = LWORK  - KEND2
C
         KTRVI0Y  = KEND2
         KEND3 = KTRVI0Y + NCKATR(ISCKB2Y)!KTRVI0Y - int. in <mu3|[H^B,T2^C]|HF>
         LWRK3    = LWORK    - KEND3      !                        ===
C
         KTRVI0X  = KEND3
         KEND3 = KTRVI0X + NCKATR(ISCKB2X)!KTRVI0X - int. in <mu3|[H^C,T2^B]|HF>
         LWRK3    = LWORK    - KEND3      !                        ===
C
         KTRVI7 = KEND3
         KEND3 = KTRVI7 + NCKATR(ISCKBBC)!KTRVI7 - int. in <mu3|[H^BC,T2^0]|HF>
         LWRK3  = LWORK  - KEND3         !                       ====
C
         KINTVI = KEND3
         KEND4  = KINTVI + NCKA(ISCKB1)!KINTVI - temporary storage
         LWRK4  = LWORK  - KEND4
C
         IF (LWRK4 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND4
            CALL QUIT('Insufficient space in CC3_BMATSD (3)')
         END IF
C
C--------------------
C        Start D-loop
C--------------------
C
         DO D = 1,NVIR(ISYMD)
C
C----------------------------------------------------------------------------
C           B-transformed virtual integrals for <mu3|[H^B,T2^C]|HF> (fixed D)
C----------------------------------------------------------------------------
C
            IOFF = ICKBD(ISCKB2Y,ISYMD) + NCKATR(ISCKB2Y)*(D - 1) + 1
            IF (NCKATR(ISCKB2Y) .GT. 0) THEN
               CALL GETWA2(LUDKBCR,FNDKBCR,WORK(KTRVI0Y),IOFF,
     &                     NCKATR(ISCKB2Y))
            ENDIF
C
C----------------------------------------------------------------------------
C           C-transformed virtual integrals for <mu3|[H^C,T2^B]|HF> (fixed D)
C----------------------------------------------------------------------------
C
            IOFF = ICKBD(ISCKB2X,ISYMD) + NCKATR(ISCKB2X)*(D - 1) + 1
            IF (NCKATR(ISCKB2X) .GT. 0) THEN
               CALL GETWA2(LUDKBCR2,FNDKBCR2,WORK(KTRVI0X),IOFF,
     &                     NCKATR(ISCKB2X))
            ENDIF
C
C-----------------------------------------------------------------------------
C           B-transformed virtual integrals for <mu3|[H^BC,T2^0]|HF> (fixed D)
C-----------------------------------------------------------------------------
C
            IOFF = ICKBD(ISCKBBC,ISYMD) + NCKATR(ISCKBBC)*(D - 1) + 1
            IF (NCKATR(ISCKBBC) .GT. 0) THEN
               CALL GETWA2(LUDKBCR1,FNDKBCR1,WORK(KTRVI7),IOFF,
     &                     NCKATR(ISCKBBC))
            ENDIF
C
C-----------------------------------------------------
C           Virtual integrals in contraction (fixed D)
C-----------------------------------------------------
C
            IOFF = ICKAD(ISCKB1,ISYMD) + NCKA(ISCKB1)*(D - 1) + 1
            IF (NCKA(ISCKB1) .GT. 0) THEN
               CALL GETWA2(LU3VI2,FN3VI2,WORK(KINTVI),IOFF,
     *                     NCKA(ISCKB1))
            ENDIF
            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI),WORK(KLAMDP),
     *                       ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
C
            IOFF = ICKAD(ISCKB1,ISYMD) + NCKA(ISCKB1)*(D - 1) + 1
            IF (NCKA(ISCKB1) .GT. 0) THEN
               CALL GETWA2(LU3VI,FN3VI,WORK(KINTVI),IOFF,
     *                     NCKA(ISCKB1))
            ENDIF
            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI1),WORK(KLAMDP),
     *                          ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
C
C---------------------------
C           Start ISYMB-loop
C---------------------------
C
            DO ISYMB = 1,NSYM
C
               ISYMBC  = MULD2H(ISYMRB,ISYMRC)
               ISYMBD  = MULD2H(ISYMB,ISYMD)
               ISWMAT  = MULD2H(ISYMBC,ISYMBD)
               ISCKD2Y = MULD2H(ISINTR2B,ISYMB)
               ISCKD2X = MULD2H(ISINTR2C,ISYMB)
               ISCKDBC = MULD2H(ISINTBC,ISYMB)
C
               !Use KEND3, because KINTVI is not needed any more
               KDIAGW  = KEND3
               KWMAT   = KDIAGW  + NCKIJ(ISWMAT)
               KINDSQW = KWMAT   + NCKIJ(ISWMAT)
               KTMAT   = KINDSQW + (6*NCKIJ(ISWMAT) - 1)/IRAT + 1
               KEND4   = KTMAT   + NCKIJ(ISWMAT)
C
               KTRVI8Y = KEND4
               KEND4   = KTRVI8Y + NCKATR(ISCKD2Y)!KTRVI8Y: <mu3|[H^B,T2^C]|HF>
               LWRK4   = LWORK   - KEND4          !               ===
C
               KTRVI8X = KEND4
               KEND4   = KTRVI8X + NCKATR(ISCKD2X)!KTRVI8X: <mu3|[H^C,T2^B]|HF>
               LWRK4   = LWORK   - KEND4          !               ===
C
               KTRVI9 = KEND4
               KEND4  = KTRVI9 + NCKATR(ISCKDBC)!KTRVI9: <mu3|[H^BC,T2^0]|HF>
               LWRK4   = LWORK   - KEND4        !              ====
C
               IF (LWRK4 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Memory available : ',LWORK
                  WRITE(LUPRI,*) 'Memory needed    : ',KEND4
                  CALL QUIT('Insufficient space in CC3_BMATSD (4)')
               END IF
C
C---------------------------------------------
C              Construct part of the diagonal.
C---------------------------------------------
C
               CALL CC3_DIAG(WORK(KDIAGW),WORK(KFOCKD),ISWMAT)
C
               IF (IPRINT .GT. 55) THEN
                  XNORMVAL = DDOT(NCKIJ(ISWMAT),WORK(KDIAGW),1,
     *                    WORK(KDIAGW),1)
                  WRITE(LUPRI,*) 'Norm of DIA  ',XNORMVAL
               ENDIF
C
C-------------------------------------
C              Construct index arrays.
C-------------------------------------
C
               LENSQW = NCKIJ(ISWMAT)
               CALL CC3_INDSQ(WORK(KINDSQW),LENSQW,ISWMAT)
C
C--------------------------
C              Start B-loop
C--------------------------
C
               DO B = 1,NVIR(ISYMB)
C
C---------------------------------------
C                 Initialise
C---------------------------------------
C
                  CALL DZERO(WORK(KWMAT),NCKIJ(ISWMAT))
C
C----------------------------------------------------------------------------
C                 B-transformed virtual integrals for <mu3|[H^B,T2^C]|HF> 
C                 (fixed B)
C----------------------------------------------------------------------------
C
                  IOFF = ICKBD(ISCKD2Y,ISYMB) + NCKATR(ISCKD2Y)*(B-1) +1
                  IF (NCKATR(ISCKD2Y) .GT. 0) THEN
                     CALL GETWA2(LUDKBCR,FNDKBCR,WORK(KTRVI8Y),IOFF,
     *                           NCKATR(ISCKD2Y))
                  ENDIF
C
C----------------------------------------------------------------------------
C                 C-transformed virtual integrals for <mu3|[H^C,T2^B]|HF> 
C                 (fixed B)
C----------------------------------------------------------------------------
C
                  IOFF = ICKBD(ISCKD2X,ISYMB) + NCKATR(ISCKD2X)*(B-1) +1
                  IF (NCKATR(ISCKD2X) .GT. 0) THEN
                     CALL GETWA2(LUDKBCR2,FNDKBCR2,WORK(KTRVI8X),IOFF,
     *                           NCKATR(ISCKD2X))
                  ENDIF
C
C-----------------------------------------------------------------------------
C                 B-transformed virtual integrals for <mu3|[H^BC,T2^0]|HF> 
C                 (fixed B)
C-----------------------------------------------------------------------------
C
                  IOFF = ICKBD(ISCKDBC,ISYMB) + NCKATR(ISCKDBC)*(B-1) +1
                  IF (NCKATR(ISCKDBC) .GT. 0) THEN
                     CALL GETWA2(LUDKBCR1,FNDKBCR1,WORK(KTRVI9),IOFF,
     *                           NCKATR(ISCKDBC))
                  ENDIF

C
C----------------------------------------------
C                 Calculate <mu3|[H^C,T2^B]|HF>
C----------------------------------------------
C
                  CALL WBD_GROUND(WORK(KT2B),ISYMRB,WORK(KTMAT),
     *                            WORK(KTRVI0X),WORK(KTRVI8X),
     *                            WORK(KTROC0X),ISINTR2C,WORK(KWMAT),
     *                            WORK(KEND4),LWRK4,
     *                            WORK(KINDSQW),LENSQW,
     *                            ISYMB,B,ISYMD,D)
C
C----------------------------------------------
C                 Calculate <mu3|[H^B,T2^C]|HF>
C----------------------------------------------
C
                  CALL WBD_GROUND(WORK(KT2C),ISYMRC,WORK(KTMAT),
     *                            WORK(KTRVI0Y),WORK(KTRVI8Y),
     *                            WORK(KTROC0Y),ISINTR2B,WORK(KWMAT),
     *                            WORK(KEND4),LWRK4,
     *                            WORK(KINDSQW),LENSQW,
     *                            ISYMB,B,ISYMD,D)
C
C----------------------------------------------
C                 Calculate <mu3|[H^BC,T2^0]|HF>
C----------------------------------------------
C
                  CALL WBD_GROUND(WORK(KT2TP),ISYMT2,WORK(KTMAT),
     *                            WORK(KTRVI7),WORK(KTRVI9),
     *                            WORK(KTROC3),ISINTBC,WORK(KWMAT),
     *                            WORK(KEND4),LWRK4,
     *                            WORK(KINDSQW),LENSQW,
     *                            ISYMB,B,ISYMD,D)
C
C----------------------------------------------------
C                 Divide by orbital energy difference
C                 and remove the forbidden elements
C----------------------------------------------------
C
C

                  CALL T3_FORBIDDEN(WORK(KWMAT),ISYMBC,
     *                              ISYMB,B,ISYMD,D)

c                  call sum_pt3(work(kwmat),isymb,b,isymd,d,
c    *                          iswmat,work(kx3am),4)
c                  write(lupri,*) 'Total norm of WBD : ',
c    *            ddot(nckij(iswmat),work(kwmat),1,work(kwmat),1)

                  CALL WBD_DIA(B,ISYMB,D,ISYMD,FREQBC,
     *                         ISWMAT,WORK(KWMAT),
     *                         WORK(KDIAGW),WORK(KFOCKD))

C
C-----------------------------------------------
C                 Carry out the contractions...
C-----------------------------------------------
C
C
C--------------------------------------------------------
C                 Calculate the  term <mu1|[H,W^BD(3)]|HF> 
C                 added in OMEGA1EFF 
C--------------------------------------------------------
C
                  CALL CC3_CY1(OMEGA1EFF,ISYRES,WORK(KWMAT),ISWMAT,
     *                         WORK(KTMAT),WORK(KXIAJB),
     *                         ISINT1,WORK(KINDSQW),LENSQW,
     *                         WORK(KEND4),LWRK4,
     *                         ISYMB,B,ISYMD,D)
C
C------------------------------------------------------
C                 Calculate the  term <mu2|[H,W^BD(3)]|HF> 
C                 ( Fock matrix cont ) 
C                 added in OMEGA2EFF 
C------------------------------------------------------
C
                  CALL CC3_CY2F(OMEGA2EFF,ISYRES,WORK(KWMAT),ISWMAT,
     *                          WORK(KTMAT),WORK(KFOCK0),ISYM0,
     *                          WORK(KINDSQW),
     *                          LENSQW,WORK(KEND4),LWRK4,
     *                          ISYMB,B,ISYMD,D)
c                 
C------------------------------------------------------
C                 Calculate the  term <mu2|[H,W^BD(3)]|HF> 
C                 ( Occupied  cont ) 
C                 added in OMEGA2EFF 
C------------------------------------------------------
C
                 CALL CC3_CY2O(OMEGA2EFF,ISYRES,WORK(KWMAT),ISWMAT,
     *                          WORK(KTMAT),WORK(KTROC),
     *                          WORK(KTROC1),ISINT1,WORK(KEND4),LWRK4,
     *                          WORK(KINDSQW),LENSQW,
     *                          ISYMB,B,ISYMD,D)
C
C
C------------------------------------------------------
C                 Calculate the  term <mu2|[H,W^BD(3)]|HF> 
C                 ( Virtual  cont ) 
C                 added in OMEGA2EFF 
C------------------------------------------------------
C
                  CALL CC3_CY2V(OMEGA2EFF,ISYRES,WORK(KRBJIA),
     *                          WORK(KWMAT),
     *                          ISWMAT,WORK(KTMAT),WORK(KTRVI),
     *                          WORK(KTRVI1),ISINT1,WORK(KEND4),LWRK4,
     *                          WORK(KINDSQW),LENSQW,
     *                          ISYMB,B,ISYMD,D)
C
               END DO !B
            END DO    !ISYMB
C
         END DO       !D
      END DO          !ISYMD
c
c      
c      write(lupri,*) 'T3XY in CC3_BMATSD  : '
c      call print_pt3(work(kx3am),1,4)
c

C
C
C------------------------------------------------------
C     Accumulate RBJIA from <mu2|[H,W^BD(3)]|HF> ( Virtual  cont ) 
C     in OMEGA2EFF 
C------------------------------------------------------
C
      CALL CC3_RBJIA(OMEGA2EFF,ISYRES,WORK(KRBJIA))
c
c     write(lupri,*)'OMEGA1EFF after CC3_BMATSD'
c     call output(OMEGA1EFF,1,nvir(1),1,nrhf(1),nvir(1),nrhf(1),1,lupri)
c
c     write(lupri,*)'OMEGA2EFF after CC3_BMATSD'
c     call outpak(OMEGA2EFF,NT1AM(1),1,lupri)

C
      IF (IPRINT .GT. 55) THEN
         XNORMVAL = DDOT(NT2AM(ISYRES),OMEGA2EFF,1,OMEGA2EFF,1)
         WRITE(LUPRI,*)'Norm of OMEGA2EFF final before added  ',XNORMVAL
      ENDIF
C
      DO I = 1,NT2AM(ISYRES)
         OMEGA2EFF(I) = OMEGA2EFF(I) + OMEGA2(I)
      END DO
C
      IF (IPRINT .GT. 55) THEN
         XNORMVAL = DDOT(NT2AM(ISYRES),OMEGA2EFF,1,OMEGA2EFF,1)
         WRITE(LUPRI,*)'Norm of OMEGA2EFF final, OMEGA2EFF + OMEGA2F  ',
     *                  XNORMVAL
      ENDIF
C
      DO I = 1,NT1AM(ISYRES)
         OMEGA1EFF(I) = OMEGA1EFF(I) + OMEGA1(I)
      END DO
C
      IF (IPRINT .GT. 55) THEN
         XNORMVAL = DDOT(NT1AM(ISYRES),OMEGA1EFF,1,OMEGA1EFF,1)
         WRITE(LUPRI,*) 'Norm of OMEGA1EFF final, OMEGA1EFF + OMEGA1  ',
     *                   XNORMVAL
      ENDIF
C
C-------------------------------------------
C     Close files (integrals in contraction)
C-------------------------------------------
C
      CALL WCLOSE2(LUTOC,FNTOC,'KEEP')
      CALL WCLOSE2(LU3VI,FN3VI,'KEEP')
      CALL WCLOSE2(LU3VI2,FN3VI2,'KEEP')
C
C------------------------------------------
C     Close temporary files (H^B integrals)
C------------------------------------------
C
      CALL WCLOSE2(LU3SRTR,FN3SRTR,'DELETE')
      CALL WCLOSE2(LUCKJDR,FNCKJDR,'DELETE')
      CALL WCLOSE2(LUDELDR,FNDELDR,'DELETE')
      CALL WCLOSE2(LUDKBCR,FNDKBCR,'DELETE')
C
C------------------------------------------
C     Close temporary files (H^BC integrals)
C------------------------------------------
C
      CALL WCLOSE2(LUCKJDR1,FNCKJDR1,'DELETE')
      CALL WCLOSE2(LUDKBCR1,FNDKBCR1,'DELETE')
C
C------------------------------------------
C     Close temporary files (H^C integrals)
C------------------------------------------
C
      CALL WCLOSE2(LUCKJDR2,FNCKJDR2,'DELETE')
      CALL WCLOSE2(LUDKBCR2,FNDKBCR2,'DELETE')
C
C-------------
C     End
C-------------
C
      CALL QEXIT('CC3_BMATSD ')
C
      RETURN
      END
C
