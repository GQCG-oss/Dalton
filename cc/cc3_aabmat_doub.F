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
C    /* Deck cc3_aabmat_doub */ 
      SUBROUTINE CC3_AABMAT_DOUB(OMEGA2,LISTB,IDLSTB,LISTC,
     *                     IDLSTC,WORK,LWORK)

******************************************************************
*
* Calculate the triples correction to the doubles part of the 
* right-hand side of second-order amplitudes equation:
*
* omega2 = omega2 + <mu2|[X,T3Y]|HF> + <mu2|[[H0,T1X],T3Y]|HF> 
*
* T3Y calculated with W intermmediate.
*
*
* 
* F. Pawlowski, P. Jorgensen, Aarhus, Spring 2003.
*
* (based on CCSDT_FBMAT routine)
*
******************************************************************
*
      IMPLICIT NONE
C
      INTEGER LWORK
      INTEGER ISINT1, ISINT2, ISYMT2, ISYMT3, ISYMWB
      INTEGER LU3SRTR, LUCKJDR, LUDELDR, LUDKBCR
      INTEGER LUCKJDR2, LUDKBCR2
      INTEGER KT2TP, KRBJIA, KFOCKD, KEND0, LWRK0, IOPT
      INTEGER IDLSTB, ISYMB, KT1B, KT2B
      INTEGER ISINTR1B, ISINTR2B, KLAMDP, KLAMDH, KEND1, LWRK1
      INTEGER KTROCY, KTROC1Y, KTROC0, KTROC02, KTROC0X
      INTEGER KINTOC, KEND2, LWRK2
      INTEGER IOFF, ISYMD, ISCKB1, ISCKB2, ISCKB2X
      INTEGER KTRVIY, KTRVI1Y, KTRVI2, KTRVI3, KTRVI0, KTRVI0X
      INTEGER KEND3, LWRK3, KEND4, LWRK4, KINTVI, LENSQ, LENSQWB
      INTEGER ISYALJ, ISYALJ2, ISYMBD, ISCKIJ, ISCKD2, ISCKD2X
      INTEGER ISWMATB, KSMAT, KSMAT3, KUMAT, KUMAT3, KDIAG, KDIAGWB
      INTEGER KWMAT, KINDSQ, KINDSQWB, KINDEX, KINDEX2, KTMAT
      INTEGER KTRVI8, KTRVI8X, KTRVI9, KTRVI10
      INTEGER LUCKJD, LUTOC, LU3VI, LU3VI2, LUDKBC, LUDELD
      INTEGER KSAVE, KEND00, LWRK00
      INTEGER IDLSTC, ISYMRB, ISYMRC, ISINT1C, ISYRES
      INTEGER KLAMPB, KLAMHB, KFCKB, LUFCK
      INTEGER KLAMPC, KLAMHC, KT1C, ISYFCKB, ISYFCKC, KFCKC
      INTEGER IRREP, IERR, KXIAJB, IOPTTCME, ISYMK, KOFF1, KOFF2
      INTEGER ISYMC, IRREP2, ISCKB1C, ISYCKB
      INTEGER ISINT1B
      INTEGER KTROCX,KTROC1X
      INTEGER ISCKB1B
      INTEGER KTRVIX,KTRVI1X
      INTEGER ISINTR2C,ISINTR1C,ISCKB2Y,ISCKD2Y
      INTEGER KTRVI0Y,KTRVI8Y,KTROC0Y
      INTEGER ISYMWC,ISWMATC
      INTEGER KT2C,KFCKBB,KFCKCC,KOFF3,MAXBC,KT3MAT,KDIAGWC,KINDSQWC
      INTEGER LENSQWC
      INTEGER ILSTSYM
c
      INTEGER kx3am
c
#if defined (SYS_CRAY)
      REAL OMEGA2(*)
      REAL WORK(LWORK)
      REAL FREQB, FREQC
      REAL XNORMVAL, DDOT
#else
      DOUBLE PRECISION OMEGA2(*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION FREQB, FREQC
      DOUBLE PRECISION XNORMVAL, DDOT
#endif
      CHARACTER*6 FNCKJD, FNDELD, FN3VI
      CHARACTER*8 FNTOC,FN3VI2
      CHARACTER*4 FNDKBC
C        
      PARAMETER (FNCKJD ='CKJDEL' , FNTOC   = 'CCSDT_OC' )
      PARAMETER (FNDELD = 'CKDELD', FNDKBC  = 'DKBC' )
      PARAMETER (FN3VI  = 'CC3_VI', FN3VI2  = 'CC3_VI12')
C
      CHARACTER*12 FN3SRTR, FNCKJDR, FNDELDR, FNDKBCR, FNDKBCR4
      CHARACTER*13 FNCKJDR2, FNDKBCR2
      CHARACTER*10 MODEL
      CHARACTER*8 LABELB, LABELC
      CHARACTER*3 LISTB, LISTC
      CHARACTER CDUMMY*1
C
#include <priunit.h>
#include <dummy.h>
#include <iratdef.h>
#include <ccsdsym.h>
#include <inftap.h>
#include <ccinftap.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccr1rsp.h>
#include <ccexci.h>
C
      PARAMETER(FN3SRTR  = 'CCSDT_FBMAT1', FNCKJDR  = 'CCSDT_FBMAT2', 
     *          FNDELDR  = 'CCSDT_FBMAT3', FNDKBCR  = 'CCSDT_FBMAT4',
     *          FNDKBCR4 = 'CCSDT_FBMAT5')
C
      PARAMETER(FNCKJDR2  = 'CCSDT_FBMAT22',FNDKBCR2  = 'CCSDT_FBMAT42')
C
      CALL QENTER('AABDB')
C
C     Lists check
C
      IF ( (LISTB(1:3).EQ.'R1 ') .AND. (LISTC(1:3).EQ.'R1 ') ) THEN
         CONTINUE
      ELSE IF ( (LISTB(1:3).EQ.'R1 ') .AND. (LISTC(1:3).EQ.'RE ') ) THEN
         CONTINUE
      ELSE
         WRITE(LUPRI,*)'LISTB = ', LISTB
         WRITE(LUPRI,*)'LISTC = ', LISTC
         WRITE(LUPRI,*)'Only implemented when  '
     *   //'LISTB = R1 and LISTC = R1 or RE '
         CALL QUIT('Case not implemented in CC3_AABMAT_DOUB')
      END IF
      
c
c     write(lupri,*)'OMEGA2 at the beginning of CC3_AABMAT_DOUB '
c     call outpak(OMEGA2,NT1AM(1),1,lupri)

C
C--------------------------
C     Save and set flags.
C--------------------------
C
      IPRCC   = IPRINT
      ISINT1  = 1
      ISINT2  = 1
      ISYMT2  = 1
      ISYMT3  = MULD2H(ISINT2,ISYMT2)
C
C--------------------------------
C     Open files
C--------------------------------
C
      LUCKJD   = -1
      LUTOC    = -1
      LU3VI    = -1
      LU3VI2   = -1
      LUDKBC   = -1
      LUDELD   = -1
C
      CALL WOPEN2(LUCKJD,FNCKJD,64,0)
      CALL WOPEN2(LUTOC,FNTOC,64,0)
      CALL WOPEN2(LU3VI,FN3VI,64,0)
      CALL WOPEN2(LU3VI2,FN3VI2,64,0)
      CALL WOPEN2(LUDKBC,FNDKBC,64,0)
      CALL WOPEN2(LUDELD,FNDELD,64,0)
C
C--------------------------------
C     Open temporary files
C--------------------------------
C
      LU3SRTR   = -1
      LUCKJDR   = -1
      LUDELDR   = -1
      LUDKBCR   = -1
      LUCKJDR2  = -1
      LUDKBCR2  = -1
C
      CALL WOPEN2(LU3SRTR,FN3SRTR,64,0)
      CALL WOPEN2(LUCKJDR,FNCKJDR,64,0)
      CALL WOPEN2(LUDELDR,FNDELDR,64,0)
      CALL WOPEN2(LUDKBCR,FNDKBCR,64,0)
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
      KEND00 = KXIAJB + NT2AM(ISINT1)
      LWRK00 = LWORK  - KEND00
C
      IF (LWRK00 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_AABMAT_DOUB (zeroth allo.')
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
      CALL GET_LAMBDA0(WORK(KLAMDP),WORK(KLAMDH),WORK(KEND00),LWRK00)
C
C-------------------------------------------
C     Read in t2 , square it and reorder 
C-------------------------------------------
C
      IOPT = 2
      CALL GET_T1_T2(IOPT,.FALSE.,DUMMY,WORK(KT2TP),'R0',0,ISYMT2,
     *               WORK(KEND00),LWRK00)
C
      IF (IPRINT .GT. 55) THEN
         XNORMVAL = DDOT(NT2SQ(ISYMT2),WORK(KT2TP),1,WORK(KT2TP),1)
         WRITE(LUPRI,*) 'Norm of T2TP ',XNORMVAL
      ENDIF
C
C--------------------------------------
C     Read in orbital energies
C--------------------------------------
C
      CALL GET_ORBEN(WORK(KFOCKD),WORK(KEND00),LWRK00)
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT If we want to sum the T3 amplitudes
COMMENT COMMENT COMMENT
      if (.false.) then
         kx3am  = kend00
         kend00 = kx3am + nt1amx*nt1amx*nt1amx
         call dzero(work(kx3am),nt1amx*nt1amx*nt1amx)
         lwrk00 = lwork - kend00
         if (lwrk00 .lt. 0) then
            write(lupri,*) 'Memory available : ',lwork
            write(lupri,*) 'Memory needed    : ',kend00
            call quit('Insufficient space (T3) in CC3_AABMAT_DOUB')
         END IF
      endif
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
C
C-------------------
C     Lists handling
C-------------------
C


      IF (LISTB(1:3).EQ.'R1 ') THEN
         ISYMRB = ISYLRT(IDLSTB)
         FREQB  = FRQLRT(IDLSTB)
         LABELB = LRTLBL(IDLSTB)
      ELSE 
         CALL QUIT('Unknown list for LISTB in CC3_AABMAT_DOUB.')
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
         CALL QUIT('Unknown list for LISTC in CC3_AABMAT_DOUB.')
      END IF
C
      ISYMWB   = MULD2H(ISYMT3,ISYMRB)
      ISYMWC   = MULD2H(ISYMT3,ISYMRC)
      ISINT1C = MULD2H(ISINT1,ISYMRC)
      ISINT1B = MULD2H(ISINT1,ISYMRB)
      ISYRES  = MULD2H(ISYMWB,ISINT1C)
      !Symmetry check
      IF ( ISYRES .NE. MULD2H(ISYMWC,ISINT1B) ) THEN
         CALL QUIT('Symmetry mismatch in CC3_AABMAT_DOUB (ISYRES)')
      END IF
C
      ISYFCKB = MULD2H(ISYMOP,ISYMRB)
      ISYFCKC = MULD2H(ISYMOP,ISYMRC)
C
C----------------------------------------------
C     Read T1^B and T2^B
C     Read T1^C and T2^C
C     Calculate (ck|de)-tilde(B) and (ck|lm)-tilde(B)
C     Used to construct T3^B
C
C     Calculate (ck|de)-tilde(C) and (ck|lm)-tilde(C)
C     Used to construct T3^C
C----------------------------------------------
C
      KT2B   = KEND00
      KEND0  = KT2B   + NT2SQ(ISYMRB)
      LWRK0  = LWORK  - KEND0
C
      KT2C   = KEND0
      KEND0  = KT2C   + NT2SQ(ISYMRC)
C
      KT1B   = KEND0
      KEND0  = KT1B   + NT1AM(ISYMRB)
      LWRK0  = LWORK  - KEND0
C
      KT1C   = KEND0
      KEND0  = KT1C   + NT1AM(ISYMRC)
      LWRK1  = LWORK  - KEND0
C
      IF (LWRK0 .LT. NT2SQ(ISYMRB)) THEN
         CALL QUIT('Out of memory in CC3_AABMAT_DOUB (TOT_T3Y) ')
      ENDIF
C
C     Readin T1B and T2B
C
      IOPT = 3
      CALL GET_T1_T2(IOPT,.TRUE.,WORK(KT1B),WORK(KT2B),LISTB,IDLSTB,
     *            ISYMRB,WORK(KEND0),LWRK0)
C
      IF (IPRINT .GT. 55) THEN
         XNORMVAL = DDOT(NT1AM(ISYMRB),WORK(KT1B),1,WORK(KT1B),1)
         WRITE(LUPRI,*) 'Norm of T1B  ',XNORMVAL
         XNORMVAL = DDOT(NT2SQ(ISYMRB),WORK(KT2B),1,WORK(KT2B),1)
         WRITE(LUPRI,*) 'Norm of T2B  ',XNORMVAL
      ENDIF
C
C     Readin T1C and T2C
C
      IOPT = 3
      CALL GET_T1_T2(IOPT,.TRUE.,WORK(KT1C),WORK(KT2C),LISTC,IDLSTC,
     *            ISYMRC,WORK(KEND0),LWRK0)
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
     *                WORK(KLAMDH),WORK(KEND0),LWRK0,
     *                LU3SRTR,FN3SRTR,LUCKJDR,FNCKJDR)
C
      CALL CC3_SORT1(WORK(KEND0),LWRK0,2,ISINTR1B,LU3SRTR,FN3SRTR,
     *               LUDELDR,FNDELDR,IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *               IDUMMY,CDUMMY)
C
      CALL CC3_SINT(WORK(KLAMDH),WORK(KEND0),LWRK0,ISINTR1B,
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
C----------------------------------------------------
C     Integrals (ck|de)-tilde(C) and (ck|lm)-tilde(C)
C----------------------------------------------------
C
      ISINTR1C = MULD2H(ISINT1,ISYMRC)
      ISINTR2C = MULD2H(ISINT2,ISYMRC)
C
      CALL CC3_BARINT(WORK(KT1C),ISYMRC,WORK(KLAMDP),
     *                WORK(KLAMDH),WORK(KEND0),LWRK0,
     *                LU3SRTR,FN3SRTR,LUCKJDR2,FNCKJDR2)
C
      CALL CC3_SORT1(WORK(KEND0),LWRK0,2,ISINTR1C,LU3SRTR,FN3SRTR,
     *               LUDELDR,FNDELDR,IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *               IDUMMY,CDUMMY)
C
      CALL CC3_SINT(WORK(KLAMDH),WORK(KEND0),LWRK0,ISINTR1C,
     *              LUDELDR,FNDELDR,LUDKBCR2,FNDKBCR2)
C
C---------------------------------------------------
C     Regain work space and create lambdaC and lambdaB
C     lambda-transformed (B and C) matrices
C     and C and B property matices.
C---------------------------------------------------
C
      ! Use KEND0 again, because neither KT1B nor KT1C are needed any more
      KFCKB  = KEND0
      KFCKC  = KFCKB  + N2BST(ISYFCKB)
      KEND0  = KFCKC  + N2BST(ISYFCKC)
C
      KFCKBB = KEND0
      KFCKCC = KFCKBB + N2BST(ISYFCKB)
      KEND0  = KFCKCC + N2BST(ISYFCKC)
C
      KLAMPC = KEND0
      KLAMHC = KLAMPC + NLAMDT
      KEND0  = KLAMHC + NLAMDT
      LWRK0  = LWORK  - KEND0
C
      KLAMPB = KEND0
      KLAMHB = KLAMPB + NLAMDT
      KEND0  = KLAMHB + NLAMDT
      LWRK0  = LWORK  - KEND0
C
      IF (LWRK0 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND0
         CALL QUIT('Insufficient space in CC3_AABMAT_DOUB')
      END IF
C
      CALL GET_LAMBDAX(WORK(KLAMPC),WORK(KLAMHC),LISTC,IDLSTC,ISYMRC,
     *                 WORK(KLAMDP),WORK(KLAMDH),WORK(KEND0),LWRK0)
C
      CALL GET_LAMBDAX(WORK(KLAMPB),WORK(KLAMHB),LISTB,IDLSTB,ISYMRB,
     *                 WORK(KLAMDP),WORK(KLAMDH),WORK(KEND0),LWRK0)
C
C---------------------------------------
C     Calculate the F^B matrix
C---------------------------------------
C
      CALL DZERO(WORK(KFCKB),N2BST(ISYFCKB))
      IF (LISTB(1:3).EQ.'R1 ') THEN
C
        CALL GET_FOCKX(WORK(KFCKB),LABELB,IDLSTB,ISYFCKB,
     *                 WORK(KLAMDP),1,WORK(KLAMDH),1,WORK(KEND0),LWRK0)
C
        IF (IPRINT .GT. 50) THEN
           CALL AROUND( 'In CC3_AABMAT_DOUB: Fock^B MO matrix' )
           CALL CC_PRFCKMO(WORK(KFCKB),ISYFCKB)
        ENDIF
C
      END IF
C
C------------------------------------------------
C     Calculate the F^C property-operator matrix 
C------------------------------------------------
C
      CALL DZERO(WORK(KFCKC),N2BST(ISYFCKC))
C
      IF (LISTC(1:3).EQ.'R1 ') THEN
C
        CALL GET_FOCKX(WORK(KFCKC),LABELC,IDLSTC,ISYFCKC,
     *                 WORK(KLAMDP),1,WORK(KLAMDH),1,WORK(KEND0),LWRK0)
C
        IF (IPRINT .GT. 50) THEN
           CALL AROUND( 'In CC3_AABMAT_DOUB: F^C property matrix' )
           CALL CC_PRFCKMO(WORK(KFCKC),ISYFCKC)
        ENDIF
C
      END IF
C
C-------------------------------------------------------------------------
C        Calculate the F^C Fock matrix and add C property matrix
C        (to be able to calculate <mu2|[C,T3B]|HF> contribution also here)
C-------------------------------------------------------------------------
C
      CALL CC3LR_MFOCK(WORK(KEND0),WORK(KT1C),WORK(KXIAJB),ISYFCKC)
C
C     Put the F_{kc} part into correct F_{pq} and add the property matrix
C
      DO ISYMC = 1, NSYM
         ISYMK = MULD2H(ISYFCKC,ISYMC)
         DO K = 1, NRHF(ISYMK)
            DO C = 1, NVIR(ISYMC)
               KOFF1 = KFCKC - 1
     *               + IFCVIR(ISYMK,ISYMC)
     *               + NORB(ISYMK)*(C - 1)
     *               + K
               KOFF2 = KEND0 - 1
     *               + IT1AM(ISYMC,ISYMK)
     *               + NVIR(ISYMC)*(K - 1)
     *               + C
               KOFF3 = KFCKCC - 1
     *               + IFCVIR(ISYMK,ISYMC)
     *               + NORB(ISYMK)*(C - 1)
     *               + K
C
               WORK(KOFF3) = WORK(KOFF1) + WORK(KOFF2)
C
            ENDDO
         ENDDO
      ENDDO
C
C-------------------------------------------------------------------------
C        Calculate the F^B Fock matrix and add B property matrix
C        (to be able to calculate <mu2|[B,T3C]|HF> contribution also here)
C-------------------------------------------------------------------------
C
      CALL CC3LR_MFOCK(WORK(KEND0),WORK(KT1B),WORK(KXIAJB),ISYFCKB)
C
C     Put the F_{kc} part into correct F_{pq} and add the property matrix
C
      DO ISYMC = 1, NSYM
         ISYMK = MULD2H(ISYFCKB,ISYMC)
         DO K = 1, NRHF(ISYMK)
            DO C = 1, NVIR(ISYMC)
               KOFF1 = KFCKB - 1
     *               + IFCVIR(ISYMK,ISYMC)
     *               + NORB(ISYMK)*(C - 1)
     *               + K
               KOFF2 = KEND0 - 1
     *               + IT1AM(ISYMC,ISYMK)
     *               + NVIR(ISYMC)*(K - 1)
     *               + C
               KOFF3 = KFCKBB - 1
     *               + IFCVIR(ISYMK,ISYMC)
     *               + NORB(ISYMK)*(C - 1)
     *               + K
C
               WORK(KOFF3) = WORK(KOFF1) + WORK(KOFF2)
C
            ENDDO
         ENDDO
      ENDDO
C
      IF (IPRINT .GT. 50) THEN
         CALL AROUND( 'In CC3_AABMAT_DOUB: Fock^C MO matrix' )
         CALL CC_PRFCKMO(WORK(KFCKC),ISYFCKC)
      ENDIF
C
C--------------------------
C     Read occupied integrals.
C--------------------------
C
C     Memory allocation.
C
C------------------------------------------------------
C     Work allocation 
C------------------------------------------------------
C
      KRBJIA = KEND0
      KTROCY  = KRBJIA + NT2SQ(ISYRES)
      KTROC1Y = KTROCY  + NTRAOC(ISINT1C)
      KTROC0 = KTROC1Y + NTRAOC(ISINT1C)
      KTROC02= KTROC0 + NTRAOC(ISINT2)
      KEND1  = KTROC02+ NTRAOC(ISINT2)
      LWRK1  = LWORK  - KEND1
C
      KTROCX  = KEND1
      KTROC1X = KTROCX  + NTRAOC(ISINT1B)
      KEND1   = KTROC1X + NTRAOC(ISINT1B)
C
      KTROC0X = KEND1
      KEND1   = KTROC0X + NTRAOC(ISINTR2B)
C
      KTROC0Y = KEND1
      KEND1   = KTROC0Y + NTRAOC(ISINTR2C)
C
      KINTOC = KEND1
      KEND2  = KINTOC  + NTOTOC(ISINT1)
      LWRK2  = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in CC3_XI')
      END IF
C
C-------------------------------------
C     Initialise the result vectors
C-------------------------------------
C
      CALL DZERO(WORK(KRBJIA),NT2SQ(ISYRES))
C
C------------------------
C     Occupied integrals.
C
C     Read in integrals for SMAT etc.
C-----------------------
C
      CALL INTOCC_T30(LUCKJD,FNCKJD,WORK(KLAMDP),ISINT2,WORK(KTROC0),
     *                WORK(KTROC02),WORK(KEND2),LWRK2)
C
C------------------------------------------
C     B transformed Occupied integrals.
C-----------------------------------------
C
      CALL INTOCC_T3X(LUCKJDR,FNCKJDR,WORK(KLAMDP),ISINTR2B,
     *                WORK(KTROC0X),WORK(KEND2),LWRK2)
C
C------------------------------------------
C     C transformed Occupied integrals.
C-----------------------------------------
C
      CALL INTOCC_T3X(LUCKJDR2,FNCKJDR2,WORK(KLAMDP),ISINTR2C,
     *                WORK(KTROC0Y),WORK(KEND2),LWRK2)
C
C------------------------
C     Occupied integrals.
C
C     Read in integrals for WMAT etc.
C-----------------------
C
      IOFF = 1
      IF (NTOTOC(ISYMOP) .GT. 0) THEN
         CALL GETWA2(LUTOC,FNTOC,WORK(KINTOC),IOFF,NTOTOC(ISYMOP))
      ENDIF
C
C----------------------------------
C     Write out norms of Integrals.
C----------------------------------
C
      IF (IPRINT .GT. 55) THEN
         XNORMVAL = DDOT(NTOTOC(ISYMOP),WORK(KINTOC),1,
     *                WORK(KINTOC),1)
         WRITE(LUPRI,*) 'Norm of CCSDT_OC-INT ',XNORMVAL
      ENDIF
C
C----------------------------------------------------------------------
C     Transform (ia|j delta) integrals to (ia|j k)-tildeC and sort as (i,j,k,a)
C----------------------------------------------------------------------
C
      CALL CCLR_TROCC(WORK(KINTOC),WORK(KTROCY),
     *                 WORK(KLAMHC),ISYMRC,
     *                 WORK(KEND2),LWRK2)
C
C----------------------------------------------------------------------
C     sort (i,j,k,a) as (a,i,j,k)
C----------------------------------------------------------------------
C
         CALL CCSDT_SRTOC2(WORK(KTROCY),WORK(KTROC1Y),ISINT1C,
     *                     WORK(KEND2),LWRK2)
C
C----------------------------------------------------------------------
C     Transform (ia|j delta) integrals to (ia|j k)-tildeB and sort as (i,j,k,a)
C----------------------------------------------------------------------
C
      CALL CCLR_TROCC(WORK(KINTOC),WORK(KTROCX),
     *                 WORK(KLAMHB),ISYMRB,
     *                 WORK(KEND2),LWRK2)
C
C----------------------------------------------------------------------
C     sort (i,j,k,a) as (a,i,j,k)
C----------------------------------------------------------------------
C
         CALL CCSDT_SRTOC2(WORK(KTROCX),WORK(KTROC1X),ISINT1B,
     *                     WORK(KEND2),LWRK2)
C
C----------------------------
C     General loop structure.
C----------------------------
C
      DO ISYMD = 1,NSYM
C
         ISYCKB  = MULD2H(ISYMOP,ISYMD)
         ISCKB2  = MULD2H(ISINT2,ISYMD)
         ISCKB2X = MULD2H(ISINTR2B,ISYMD)
         ISCKB2Y = MULD2H(ISINTR2C,ISYMD)
         ISCKB1  = MULD2H(ISINT1,ISYMD)
         ISCKB1C = MULD2H(ISINT1C,ISYMD)
         ISCKB1B = MULD2H(ISINT1B,ISYMD)
C
         IF (IPRINT .GT. 55) THEN
C
            WRITE(LUPRI,*) 'In CC3_AABMAT_DOUB: ISYCKB :',ISYCKB
            WRITE(LUPRI,*) 'In CC3_AABMAT_DOUB: ISCKB2 :',ISCKB2
            WRITE(LUPRI,*) 'In CC3_AABMAT_DOUB: ISCKB2X:',ISCKB2X
            WRITE(LUPRI,*) 'In CC3_AABMAT_DOUB: ISCKB1 :',ISCKB1
            WRITE(LUPRI,*) 'In CC3_AABMAT_DOUB: ISCKB1C:',ISCKB1C
C
         ENDIF

C
C--------------------------
C        Memory allocation.
C--------------------------
C
         KTRVIY  = KEND1
         KTRVI1Y = KTRVIY  + NCKATR(ISCKB1C)
         KTRVI3 = KTRVI1Y + NCKATR(ISCKB1C)
         KTRVI2 = KTRVI3 + NCKATR(ISCKB2)
         KEND2  = KTRVI2 + NCKATR(ISCKB2)
         LWRK2  = LWORK  - KEND2
C
         KTRVIX  = KEND2
         KTRVI1X = KTRVIX  + NCKATR(ISCKB1B)
         KEND2 = KTRVI1X + NCKATR(ISCKB1B)
C
         KTRVI0  = KEND2
         KEND3   = KTRVI0  + NCKATR(ISCKB2)
         LWRK3   = LWORK  - KEND3
C
         KTRVI0X  = KEND3
         KEND3    = KTRVI0X  + NCKATR(ISCKB2X)
         LWRK3    = LWORK    - KEND3
C
         KTRVI0Y  = KEND3
         KEND3    = KTRVI0Y  + NCKATR(ISCKB2Y)
         LWRK3    = LWORK    - KEND3
C
         KINTVI = KEND3
         KEND4  = KINTVI + NCKA(ISYCKB)
         LWRK4  = LWORK  - KEND4
C
         IF (LWRK4 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND4
            CALL QUIT('Insufficient space in CC3_XI')
         END IF
C
C---------------------
C        Sum over D
C---------------------
C
         DO D = 1,NVIR(ISYMD)
C
C------------------------------------------------------------
C           Read and transform integrals used in contraction
C           with W intermediate.
C------------------------------------------------------------
C
            IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D - 1) + 1
            IF (NCKA(ISYCKB) .GT. 0) THEN
               CALL GETWA2(LU3VI2,FN3VI2,WORK(KINTVI),IOFF,
     &                     NCKA(ISYCKB))
            ENDIF
C
            ! Transform for C operator
            CALL CCLR_TRVIR(WORK(KINTVI),WORK(KTRVIY),
     *                       WORK(KLAMPC),ISYMRC,
     *                       ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
C
            ! Transform for B operator
            CALL CCLR_TRVIR(WORK(KINTVI),WORK(KTRVIX),
     *                       WORK(KLAMPB),ISYMRB,
     *                       ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
C
            IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D - 1) + 1
            IF (NCKA(ISYCKB) .GT. 0) THEN
               CALL GETWA2(LU3VI,FN3VI,WORK(KINTVI),IOFF,
     &                     NCKA(ISYCKB))
            ENDIF
C
            ! Transform for C operator
            CALL CCLR_TRVIR(WORK(KINTVI),WORK(KTRVI1Y),
     *                       WORK(KLAMPC),ISYMRC,
     *                       ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
C
            ! Transform for B operator
            CALL CCLR_TRVIR(WORK(KINTVI),WORK(KTRVI1X),
     *                       WORK(KLAMPB),ISYMRB,
     *                       ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
            

C
            CALL INTVIR_T30_D(LUDKBC,FNDKBC,LUDELD,FNDELD,ISINT2,
     *                        WORK(KTRVI0),WORK(KTRVI2),WORK(KTRVI3),
     *                        WORK(KLAMDH),ISYMD,D,WORK(KEND4),LWRK4)
C
C-----------------------------------------------------------------------
C           Read B transformed virtual integrals used for W for TOT_T3Y
C-----------------------------------------------------------------------
C
            IOFF = ICKBD(ISCKB2X,ISYMD) + NCKATR(ISCKB2X)*(D - 1) + 1
            IF (NCKATR(ISCKB2X) .GT. 0) THEN
               CALL GETWA2(LUDKBCR,FNDKBCR,WORK(KTRVI0X),IOFF,
     &                     NCKATR(ISCKB2X))
            ENDIF
C
C-----------------------------------------------------------------------
C           Read C transformed virtual integrals used for W for TOT_T3Y
C-----------------------------------------------------------------------
C
            IOFF = ICKBD(ISCKB2Y,ISYMD) + NCKATR(ISCKB2Y)*(D - 1) + 1
            IF (NCKATR(ISCKB2Y) .GT. 0) THEN
               CALL GETWA2(LUDKBCR2,FNDKBCR2,WORK(KTRVI0Y),IOFF,
     &                     NCKATR(ISCKB2Y))
            ENDIF
C
C-----------------------------------------------------
C           Sum over ISYMB
C-----------------------------------------------------
C
            DO ISYMB = 1,NSYM
C
               ISYALJ  = MULD2H(ISYMB,ISYMT2)
               ISYALJ2 = MULD2H(ISYMD,ISYMT2)
               ISYMBD  = MULD2H(ISYMB,ISYMD)
               ISCKIJ  = MULD2H(ISYMBD,ISYMT3)
               ISCKD2  = MULD2H(ISINT2,ISYMB)
               ISCKD2X = MULD2H(ISINTR2B,ISYMB)
               ISCKD2Y = MULD2H(ISINTR2C,ISYMB)
               ISWMATB  = MULD2H(ISYMRB,ISCKIJ)
               ISWMATC  = MULD2H(ISYMRC,ISCKIJ)
C
               IF (IPRINT .GT. 55) THEN
C 
                  WRITE(LUPRI,*) 'In CC3_AABMAT_DOUB: ISYMD  :',ISYMD
                  WRITE(LUPRI,*) 'In CC3_AABMAT_DOUB: ISYMB  :',ISYMB
                  WRITE(LUPRI,*) 'In CC3_AABMAT_DOUB: ISYALJ :',ISYALJ
                  WRITE(LUPRI,*) 'In CC3_AABMAT_DOUB: ISYALJ2:',ISYALJ2
                  WRITE(LUPRI,*) 'In CC3_AABMAT_DOUB: ISYMBD :',ISYMBD
                  WRITE(LUPRI,*) 'In CC3_AABMAT_DOUB: ISCKIJ :',ISCKIJ
                  WRITE(LUPRI,*) 'In CC3_AABMAT_DOUB: ISWMATB :',ISWMATB
C
               ENDIF
C
C
               MAXBC   = MAX(NCKIJ(ISWMATB),NCKIJ(ISWMATC))
C
C              Can use kend3 since we do not need the integrals anymore.
               KSMAT   = KEND3
               KSMAT3  = KSMAT   + NCKIJ(ISCKIJ)
               KUMAT   = KSMAT3  + NCKIJ(ISCKIJ)
               KUMAT3  = KUMAT   + NCKIJ(ISCKIJ)
               KDIAG   = KUMAT3  + NCKIJ(ISCKIJ)
               KDIAGWB  = KDIAG   + NCKIJ(ISCKIJ)
               KWMAT   = KDIAGWB  + NCKIJ(ISWMATB)
               KINDSQWB = KWMAT   + MAXBC
               KINDSQ  = KINDSQWB + (6*NCKIJ(ISWMATB) - 1)/IRAT + 1
               KINDEX  = KINDSQ  + (6*NCKIJ(ISCKIJ) - 1)/IRAT + 1
               KINDEX2 = KINDEX  + (NCKI(ISYALJ) - 1)/IRAT + 1
               KTMAT   = KINDEX2 + (NCKI(ISYALJ2) - 1)/IRAT + 1
               KT3MAT  = KTMAT   + MAX(NCKIJ(ISCKIJ),NCKIJ(ISWMATB))
               KTRVI8  = KT3MAT   + MAX(NCKIJ(ISCKIJ),NCKIJ(ISWMATC))
               KTRVI9  = KTRVI8  + NCKATR(ISCKD2)
               KTRVI10 = KTRVI9  + NCKATR(ISCKD2)
               KEND4   = KTRVI10 + NCKATR(ISCKD2)
               LWRK4   = LWORK   - KEND4
C
               KDIAGWC  = KEND4
               KINDSQWC = KDIAGWC + NCKIJ(ISWMATC)
               KEND4    = KINDSQWC + (6*NCKIJ(ISWMATC) - 1)/IRAT + 1
C
               KTRVI8X = KEND4
               KEND4   = KTRVI8X + NCKATR(ISCKD2X)
               LWRK4   = LWORK   - KEND4
C
               KTRVI8Y = KEND4
               KEND4   = KTRVI8Y + NCKATR(ISCKD2Y)
               LWRK4   = LWORK   - KEND4
C
               IF (LWRK4 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Memory available : ',LWORK
                  WRITE(LUPRI,*) 'Memory needed    : ',KEND4
                  CALL QUIT('Insufficient space in CC3_XI')
               END IF
C
C---------------------------------------------
C              Construct part of the diagonal.
C---------------------------------------------
C
               CALL CC3_DIAG(WORK(KDIAG),WORK(KFOCKD),ISCKIJ)
               CALL CC3_DIAG(WORK(KDIAGWB),WORK(KFOCKD),ISWMATB)
               CALL CC3_DIAG(WORK(KDIAGWC),WORK(KFOCKD),ISWMATC)
C
               IF (IPRINT .GT. 55) THEN
                  XNORMVAL = DDOT(NCKIJ(ISCKIJ),WORK(KDIAG),1,
     *                    WORK(KDIAG),1)
                  WRITE(LUPRI,*) 'Norm of DIA  ',XNORMVAL
               ENDIF
C
C-------------------------------------
C              Construct index arrays.
C-------------------------------------
C
               LENSQ = NCKIJ(ISCKIJ)
               CALL CC3_INDSQ(WORK(KINDSQ),LENSQ,ISCKIJ)
               CALL CC3_INDEX(WORK(KINDEX),ISYALJ)
               CALL CC3_INDEX(WORK(KINDEX2),ISYALJ2)
               LENSQWB = NCKIJ(ISWMATB)
               CALL CC3_INDSQ(WORK(KINDSQWB),LENSQWB,ISWMATB) 
               LENSQWC = NCKIJ(ISWMATC)
               CALL CC3_INDSQ(WORK(KINDSQWC),LENSQWC,ISWMATC) 
C
               DO B = 1,NVIR(ISYMB)
C
C------------------------------------
C                 Integrals for t30
C------------------------------------
C
                  CALL INTVIR_T30_D(LUDKBC,FNDKBC,LUDELD,FNDELD,
     *                              ISINT2,WORK(KTRVI8),WORK(KTRVI9),
     *                              WORK(KTRVI10),WORK(KLAMDH),
     *                              ISYMB,B,WORK(KEND4),LWRK4)
C
C--------------------------------------------------------------------
C                 Read B transformed virtual integrals used in W
C--------------------------------------------------------------------
C
                  IOFF = ICKBD(ISCKD2X,ISYMB) + 
     &                   NCKATR(ISCKD2X)*(B - 1) + 1
                  IF (NCKATR(ISCKD2X) .GT. 0) THEN
                     CALL GETWA2(LUDKBCR,FNDKBCR,WORK(KTRVI8X),IOFF,
     &                           NCKATR(ISCKD2X))
                  ENDIF
C
C--------------------------------------------------------------------
C                 Read C transformed virtual integrals used in W
C--------------------------------------------------------------------
C
                  IOFF = ICKBD(ISCKD2Y,ISYMB) + 
     &                   NCKATR(ISCKD2Y)*(B - 1) + 1
                  IF (NCKATR(ISCKD2Y) .GT. 0) THEN
                     CALL GETWA2(LUDKBCR2,FNDKBCR2,WORK(KTRVI8Y),IOFF,
     &                           NCKATR(ISCKD2Y))
                  ENDIF
C
                  CALL GET_T30_BD(ISYMT3,ISINT2,WORK(KT2TP),ISYMT2,
     *                            WORK(KTMAT),WORK(KFOCKD),
     *                            WORK(KDIAG),WORK(KINDSQ),LENSQ,
     *                            WORK(KSMAT),WORK(KTRVI0),
     *                            WORK(KTRVI2),WORK(KTROC0),
     *                            WORK(KINDEX),WORK(KSMAT3),
     *                            WORK(KTRVI8),WORK(KTRVI9),
     *                            WORK(KINDEX2),WORK(KUMAT),
     *                            WORK(KTRVI3),WORK(KTROC02),
     *                            WORK(KUMAT3),WORK(KTRVI10),
     *                            ISYMB,B,ISYMD,D,
     *                            ISCKIJ,WORK(KEND4),LWRK4)

c                 call sum_pt3(work(KTMAT),isymb,b,isymd,d,
c    *                        isckij,work(kx3am),4)

C
                  
                  ! Copy T3^0 from KTMAT to KT3MAT
                  ! KTMAT is used as help array and overwritten.
                  ! But we still need T3^0  for WC (and there we use KT3MAT)
                  CALL DCOPY(NCKIJ(ISCKIJ),WORK(KTMAT),1,
     *                       WORK(KT3MAT),1)
C

C
C-----------------------------------------------------------------------
C                 Based on T30 BD amplitudes calculate WB BD intermediates
C-----------------------------------------------------------------------
C
C
                  CALL DZERO(WORK(KWMAT),MAXBC)
C
                  CALL GET_T3X_BD(ISYMRB,WORK(KWMAT),ISWMATB,
     *                           WORK(KTMAT),ISCKIJ,WORK(KFCKB),
     *                           ISYMRB,WORK(KT2TP),ISYMT2,
     *                           WORK(KT2B),ISYMRB,
     *                           WORK(KTRVI0),WORK(KTRVI8),
     *                           WORK(KTROC0),ISINT2,
     *                           WORK(KTRVI0X),WORK(KTRVI8X),
     *                           WORK(KTROC0X),ISINTR2B,
     *                           FREQB,WORK(KDIAGWB),WORK(KFOCKD),
     *                           WORK(KINDSQWB),LENSQWB,
     *                           B,ISYMB,D,ISYMD,WORK(KEND4),LWRK4)
c
 
c                 call sum_pt3(work(kwmat),isymb,b,isymd,d,
c    *                        isckij,work(kx3am),4)
c
  
C
C---------------------------------------------------------
C                 Calculate the  term <mu2|[H,W^BD(3)]|HF> 
C                 ( Fock matrix cont ) 
C                 added in OMEGA2
C---------------------------------------------------------
C
                  CALL CC3_CY2F(OMEGA2,ISYRES,WORK(KWMAT),
     *                          ISWMATB,WORK(KTMAT),WORK(KFCKCC),
     *                          ISYFCKC,WORK(KINDSQWB),LENSQWB,
     *                          WORK(KEND4),LWRK4,
     *                          ISYMB,B,ISYMD,D)
C
C---------------------------------------------------------
C                 Calculate the  term <mu2|[H,W^BD(3)]|HF> 
C                 ( Occupied  cont ) 
C                 added in OMEGA2 
C---------------------------------------------------------
C
                  CALL CC3_CY2O(OMEGA2,ISYRES,WORK(KWMAT),
     *                          ISWMATB,WORK(KTMAT),WORK(KTROCY),
     *                          WORK(KTROC1Y),ISINT1C,
     *                          WORK(KEND4),LWRK4,
     *                          WORK(KINDSQWB),LENSQWB,
     *                          ISYMB,B,ISYMD,D)
C
C---------------------------------------------------------
C                 Calculate the  term <mu2|[H,W^BD(3)]|HF> 
C                 ( Virtual cont ) 
C                 added in OMEGA2 
C---------------------------------------------------------
C
                  CALL CC3_CY2V(OMEGA2,ISYRES,WORK(KRBJIA),
     *                          WORK(KWMAT),ISWMATB,WORK(KTMAT),
     *                          WORK(KTRVIY),WORK(KTRVI1Y),ISINT1C,
     *                          WORK(KEND4),LWRK4,WORK(KINDSQWB),
     *                          LENSQWB,ISYMB,B,ISYMD,D)
C
C-----------------------------------------------------------------------
C                 Based on T30 BD amplitudes calculate WC BD intermediates
C-----------------------------------------------------------------------
C
C
                  CALL DZERO(WORK(KWMAT),MAXBC)
C
                  IF (LISTC(1:3).EQ.'R1 ') THEN
                    CALL GET_T3X_BD(ISYMRC,WORK(KWMAT),ISWMATC,
     *                             WORK(KT3MAT),ISCKIJ,WORK(KFCKC),
     *                             ISYMRC,WORK(KT2TP),ISYMT2,
     *                             WORK(KT2C),ISYMRC,
     *                             WORK(KTRVI0),WORK(KTRVI8),
     *                             WORK(KTROC0),ISINT2,
     *                             WORK(KTRVI0Y),WORK(KTRVI8Y),
     *                             WORK(KTROC0Y),ISINTR2C,
     *                             FREQC,WORK(KDIAGWC),WORK(KFOCKD),
     *                             WORK(KINDSQWC),LENSQWC,
     *                             B,ISYMB,D,ISYMD,WORK(KEND4),LWRK4)
c
c
c                   call sum_pt3(work(kwmat),isymb,b,isymd,d,
c    *                          isckij,work(kx3am),4)
c
                  END IF
C
                  IF (LISTC(1:3).EQ.'RE ') THEN
                    CALL GET_E3_BD(ISYMRC,WORK(KWMAT),ISWMATC,
     *                             WORK(KT3MAT),
     *                             WORK(KT2TP),ISYMT2,
     *                             WORK(KT2C),ISYMRC,
     *                             WORK(KTRVI0),WORK(KTRVI8),
     *                             WORK(KTROC0),ISINT2,
     *                             WORK(KTRVI0Y),WORK(KTRVI8Y),
     *                             WORK(KTROC0Y),ISINTR2C,
     *                             FREQC,WORK(KDIAGWC),WORK(KFOCKD),
     *                             WORK(KINDSQWC),LENSQWC,
     *                             B,ISYMB,D,ISYMD,WORK(KEND4),LWRK4)
c
c
c                   call sum_pt3(work(kwmat),isymb,b,isymd,d,
c    *                          isckij,work(kx3am),4)
c
                  END IF
  
C
C---------------------------------------------------------
C                 Calculate the  term <mu2|[H,W^BD(3)]|HF> 
C                 ( Fock matrix cont ) 
C                 added in OMEGA2
C---------------------------------------------------------
C
                  CALL CC3_CY2F(OMEGA2,ISYRES,WORK(KWMAT),
     *                          ISWMATC,WORK(KTMAT),WORK(KFCKBB),
     *                          ISYFCKB,WORK(KINDSQWC),LENSQWC,
     *                          WORK(KEND4),LWRK4,
     *                          ISYMB,B,ISYMD,D)
c
c---------------------------------------------------------
c                 Calculate the  term <mu2|[H,W^BD(3)]|HF> 
c                 ( Occupied  cont ) 
c                 added in OMEGA2 
c---------------------------------------------------------
c
                  CALL CC3_CY2O(OMEGA2,ISYRES,WORK(KWMAT),
     *                          ISWMATC,WORK(KTMAT),WORK(KTROCX),
     *                          WORK(KTROC1X),ISINT1B,
     *                          WORK(KEND4),LWRK4,
     *                          WORK(KINDSQWC),LENSQWC,
     *                          ISYMB,B,ISYMD,D)
c
c---------------------------------------------------------
c                 Calculate the  term <mu2|[H,W^BD(3)]|HF> 
c                 ( Virtual cont ) 
c                 added in OMEGA2 
c---------------------------------------------------------
 
                  CALL CC3_CY2V(OMEGA2,ISYRES,WORK(KRBJIA),
     *                          WORK(KWMAT),ISWMATC,WORK(KTMAT),
     *                          WORK(KTRVIX),WORK(KTRVI1X),ISINT1B,
     *                          WORK(KEND4),LWRK4,WORK(KINDSQWC),
     *                          LENSQWC,ISYMB,B,ISYMD,D)
C 
               END DO  ! B
            END DO     ! ISYMB
C
         END DO        ! D
      END DO           ! ISYMD
C

C
C------------------------------------------------------
C     Accumulate RBJIA from <mu2|[H,W^BD(3)]|HF> ( Vccupied  cont ) 
C     in OMEGA2 
C------------------------------------------------------
C
      CALL CC3_RBJIA(OMEGA2,ISYRES,WORK(KRBJIA))
C
      IF (IPRINT .GT. 55) THEN
         XNORMVAL = DDOT(NT2AM(ISYRES),OMEGA2,1,OMEGA2,1)
         WRITE(LUPRI,*) 'Norm of final OMEGA2 ',XNORMVAL
      ENDIF
C
C--------------------------------
C     Close files
C--------------------------------
C
      CALL WCLOSE2(LUCKJD,FNCKJD,'KEEP')
      CALL WCLOSE2(LUTOC,FNTOC,'KEEP')
      CALL WCLOSE2(LU3VI,FN3VI,'KEEP')
      CALL WCLOSE2(LU3VI2,FN3VI2,'KEEP')
      CALL WCLOSE2(LUDKBC,FNDKBC,'KEEP')
      CALL WCLOSE2(LUDELD,FNDELD,'KEEP')
C
C--------------------------------
C     Close files for "response"
C--------------------------------
C
      CALL WCLOSE2(LU3SRTR,FN3SRTR,'DELETE')
      CALL WCLOSE2(LUCKJDR,FNCKJDR,'DELETE')
      CALL WCLOSE2(LUDELDR,FNDELDR,'DELETE')
      CALL WCLOSE2(LUDKBCR,FNDKBCR,'DELETE')
      CALL WCLOSE2(LUCKJDR2,FNCKJDR2,'DELETE')
      CALL WCLOSE2(LUDKBCR2,FNDKBCR2,'DELETE')
C
C-------------
C     End
C-------------
C
      CALL QEXIT('AABDB')
C
      RETURN
      END
