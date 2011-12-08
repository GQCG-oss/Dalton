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
C    /* Deck cc3_ft3b */ 
      SUBROUTINE CC3_FT3B(LISTL,IDLSTL,LISTB,IDLSTB,IOPTRES,
     *                    OMEGA1,ISYRES,WORK,LWORK,
     *                    FNCKJD,FNDKBC,FNDELD,FNTOC)

C
C     Calculate <L2|[[H,E_AI],R_3]|HF>
C
C     Initially construct 
C
C     t3(ai,Bj,Ck) = S^BC(aikj) + U^BC(aikj) + S^CB(aijk) + U^CB(aijk)
C
C     W^BC(aikj) = sum_d X(ad)*t3(di,Bj,Ck) - sum_l X(li)*t3(al,Bj,Ck)
C                 - sum_ld X(ld)*t2(ai,Cl)*t2(Bj,dk) 
C                 - sum_ld X(ld)*t2(ai,Bl)*t2(Ck,dj) 
C                 + ( -<mu3 | [[H,T1^Y],T2] | HF> -<mu3| [H,T2^Y] |HF> )
C
C     Written by Poul Jorgensen and Filip Pawlowski, Aarhus, Fall 2002
C
      IMPLICIT NONE
C
      INTEGER LWORK, IOPTRES
      INTEGER ISINT1, ISINT2, ISYMT2, ISYMT3, ISYMW
      INTEGER LU3SRTR, LUCKJDR, LUDELDR, LUDKBCR, LUDKBCR4
      INTEGER KT1AM, KT2TP, KFOCKD, KEND0, LWRK0, IOPT
      INTEGER IDLSTB, ISYMB, KT1B, KT2B
      INTEGER ISINTR1, ISINTR2, KLAMDP, KLAMDH, KEND1, LWRK1
      INTEGER KTROC, KTROC1, KW3XOG1, KT3OG2, KW3XOGX1
      INTEGER KINTOC, KEND2, LWRK2
      INTEGER IOFF, ISYMD, ISCKB1, ISCKB2, ISCKB2Y
      INTEGER KTRVI, KTRVI1, KT3VDG2, KT3VDG3, KW3XVDG1, KW3XVDGX1
      INTEGER KEND3, LWRK3, KEND4, LWRK4, KINTVI, LENSQ, LENSQW
      INTEGER ISYALJ, ISYALJ2, ISYMBD, ISCKIJ, ISCKD2, ISCKD2Y
      INTEGER ISWMAT, KSMAT, KSMAT3, KUMAT, KUMAT3, KDIAG, KDIAGW
      INTEGER KWMAT, KINDSQ, KINDSQW, KINDEX, KINDEX2, KTMAT
      INTEGER KW3XVBG1, KW3XVDGX2, KT3VBG2, KT3VBG3, KEND5, LWRK5
      INTEGER LUCKJD, LUTOC, LUDKBC, LUDELD
      INTEGER KSAVE, KEND00, LWRK00, KEND01, LWRK01
      INTEGER ISYMRB, ISYRES
      INTEGER KLAMPB, KLAMHB, KFCKB, LUFCK, IVEC, IDLSTL, ISYML
      INTEGER ISYFCKB
      INTEGER IRREP, IERR, KXIAJB, KYIAJB, IOPTTCME, ISYMK, KOFF1, KOFF2
      INTEGER IRREP2, ISYCKB
      INTEGER KXAIJB, KYAIJB, KYAIBJ
      INTEGER ISINT3, KL2TP, ISYML2
C
C     Functions :
C
      INTEGER ILSTSYM
C
      integer kx3am
C
#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL OMEGA1(*)
      REAL FREQB,  TCON
      REAL XNORMVAL, DDOT, HALF, TWO, XT2TP, XUMAT
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION OMEGA1(*)
      DOUBLE PRECISION FREQB,  TCON
      DOUBLE PRECISION XNORMVAL, DDOT, HALF, TWO, XT2TP, XUMAT
#endif
      CHARACTER*(*) FNCKJD, FNDKBC, FNDELD, FNTOC
      CHARACTER*12 FN3SRTR, FNCKJDR, FNDELDR, FNDKBCR, FNDKBCR4
      CHARACTER*10 MODEL
      CHARACTER*8 LABELB
      CHARACTER*3 LISTL, LISTB
      CHARACTER CDUMMY*1
C
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
      PARAMETER(HALF = 0.5D0, TWO = 2.0D0)
      PARAMETER(FN3SRTR  = 'CCSDT_FBMAT1', FNCKJDR  = 'CCSDT_FBMAT2', 
     *          FNDELDR  = 'CCSDT_FBMAT3', FNDKBCR  = 'CCSDT_FBMAT4',
     *          FNDKBCR4 = 'CCSDT_FBMAT5')
C
      CALL QENTER('CC3_FT3B')
C
C--------------------------
C     Save and set flags.
C--------------------------
C
C     Set symmetry flags.
C
C
C     ISYMT2 is symmetry of T2TP
C     ISINT2 is symmetry of integrals in triples equation (ISINT2=1)
C     ISINT1 is symmetry of integrals in contraction (ISINT1=1)
C     ISYMT3  is symmetry of S and U intermediate
C     ISYMW   is symmetry of W intermediate
C     ISINT3  is symmetry of integrals used in contraction with W
C     ISYRES  is symmetry of result vector 
C
C-------------------------------------------------------------
C
COMMENT COMMENT COMMENT 
c
c    summing up contributions in omega1eff
c    must be changed when introduced in real implementation
c
c       CALL DZERO(OMEGA1EFF,NT1AM(1))
c
COMMENT COMMENT COMMENT 
      IPRCC   = IPRINT
      ISINT1  = 1
      ISINT2  = 1
      ISINT3  = 1
      ISYMT2  = 1
      ISYMT3  = MULD2H(ISINT2,ISYMT2)
C
C--------------------------------
C     Open files
C--------------------------------
C
      LUCKJD   = -1
      LUTOC    = -1
      LUDKBC   = -1
      LUDELD   = -1
C
      CALL WOPEN2(LUCKJD,FNCKJD,64,0)
      CALL WOPEN2(LUTOC,FNTOC,64,0)
      CALL WOPEN2(LUDKBC,FNDKBC,64,0)
      CALL WOPEN2(LUDELD,FNDELD,64,0)
C
C--------------------------------
C     Open temporary files
C--------------------------------
C
      LU3SRTR  = -1
      LUCKJDR  = -1
      LUDELDR  = -1
      LUDKBCR  = -1
      LUDKBCR4 = -1
C
      CALL WOPEN2(LU3SRTR,FN3SRTR,64,0)
      CALL WOPEN2(LUCKJDR,FNCKJDR,64,0)
      CALL WOPEN2(LUDELDR,FNDELDR,64,0)
      CALL WOPEN2(LUDKBCR,FNDKBCR,64,0)
      CALL WOPEN2(LUDKBCR4,FNDKBCR4,64,0)
C
C----------------------------------------------
C     Calculate the zeroth order stuff once
C----------------------------------------------
C
      KT2TP  = 1 
      KL2TP  = KT2TP  + NT2SQ(ISYMT2)
      KFOCKD = KL2TP  + NT2SQ(ISYMT2)
      KLAMDP = KFOCKD + NORBTS
      KLAMDH = KLAMDP + NLAMDT
      KXIAJB = KLAMDH + NLAMDT
      KYIAJB = KXIAJB + NT2AM(ISINT1)
      KXAIJB = KYIAJB + NT2AM(ISINT1)
      KYAIJB = KXAIJB + NT2SQ(ISINT1)
      KYAIBJ = KYAIJB + NT2SQ(ISINT1) 
      KEND00 = KYAIBJ + NT2SQ(ISINT1)
      LWRK00 = LWORK  - KEND00
C
      KT1AM  = KEND00
      KEND01 = KT1AM  + NT1AM(ISYMT2)
      LWRK01 = LWORK  - KEND01
C
      IF (LWRK01 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_FT3B (zeroth allo.')
      ENDIF
C
C-----------------------
C read in doubles component of L0  
C-----------------------
C
      ISYML2  = ILSTSYM(LISTL,IDLSTL)
      IF (LWRK01 .LT. NT2SQ(ISYML2)) THEN
        CALL QUIT('Not enough memory to construct L2TP in CC3_FT3B')
      ENDIF

      IOPT = 2
      CALL GET_T1_T2(IOPT,.FALSE.,DUMMY,WORK(KL2TP),LISTL,IDLSTL,ISYML2,
     *               WORK(KEND01),LWRK01)
C
C------------------------
C     Construct L(ia,jb).
C------------------------
C
      REWIND(LUIAJB)
      CALL READI(LUIAJB,IRAT*NT2AM(ISINT1),WORK(KXIAJB))
      CALL DCOPY(NT2AM(ISINT1),WORK(KXIAJB),1,WORK(KYIAJB),1)
C
      IOPTTCME = 1
      CALL CCSD_TCMEPK(WORK(KYIAJB),1.0D0,ISINT1,IOPTTCME)
C
      IF ( IPRINT .GT. 55) THEN
         XNORMVAL = DDOT(NT2AM(ISINT1),WORK(KYIAJB),1,
     *                WORK(KYIAJB),1)
         WRITE(LUPRI,*) 'Norm of L(IAJB) ',XNORMVAL
      ENDIF
      CALL CC_T2SQ(WORK(KXIAJB),WORK(KEND01),ISINT1)
      CALL CC_T2SQ(WORK(KYIAJB),WORK(KYAIBJ),ISINT1)
C
      CALL CC3_T2TP(WORK(KXAIJB),WORK(KEND01),ISINT1)
      CALL CC3_T2TP(WORK(KYAIJB),WORK(KYAIBJ),ISINT1)
C
      IF (IPRINT .GT. 55) THEN
         XT2TP = DDOT(NT2SQ(ISINT1),WORK(KYAIJB),1,WORK(KYAIJB),1)
         WRITE(LUPRI,*) 'Norm of L(AIJB) ', XT2TP
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
C--------------------------------------
C     Read in orbital energies
C--------------------------------------
C
C
      CALL GET_ORBEN(WORK(KFOCKD),WORK(KEND01),LWRK01)
C If we want to sum the T3 amplitudes
      if (.false.) then
         kx3am  = kend00
         kend00 = kx3am + nt1amx*nt1amx*nt1amx
         call dzero(work(kx3am),nt1amx*nt1amx*nt1amx)
         lwrk00 = lwork - kend00
         if (lwrk00 .lt. 0) then
            write(lupri,*) 'Memory available : ',lwork
            write(lupri,*) 'Memory needed    : ',kend00
            call quit('Insufficient space (T3) in CC3_FT3B')
         END IF
      endif
C
C determining R1 labels
C
      IF (LISTB(1:3).EQ.'R1 ') THEN
         ISYMRB = ISYLRT(IDLSTB)
         FREQB  = FRQLRT(IDLSTB)
         LABELB = LRTLBL(IDLSTB)
      ELSE IF (LISTB(1:3).EQ.'RE ') THEN
         ISYMRB  = ILSTSYM(LISTB,IDLSTB)
         FREQB  = EIGVAL(IDLSTB)
         LABELB = '- none -'
C
      ELSE
         CALL QUIT('Unknown list in CC3_FT3B.') 
      END IF
C
      !symmetry check
c     ISYMW   = MULD2H(ISYMT3,ISYMRB)
c     IF (ISYRES.NE.MULD2H(ISYMW,ISINT3)) THEN
c        WRITE(LUPRI,*) 'INCONSISTENT SYMMETRY SPECIFICATION'
c        WRITE(LUPRI,*) 
c    *      'ISYRES =',ISYRES,'ISYMW =',ISYMW,'ISINT3 =',ISINT3
c        STOP' CC3_FT3B inconsistent symmetry specification' 
c     END IF
C
      ISYFCKB = MULD2H(ISYMOP,ISYMRB)
C
C-------------------------------------------------
C     Read T1^B and T2^B
C     Calculate (ck|de)-tilde(B) and (ck|lm)-tilde(B)
C     Used to construct T3^B
C-------------------------------------------------
C
      KT2B   = KEND00
      KFCKB  = KT2B   + NT2SQ(ISYMRB)
      KEND0  = KFCKB  + N2BST(ISYFCKB)
      LWRK0  = LWORK  - KEND0
C
      KT1B   = KEND0
      KLAMPB = KT1B   + NT1AM(ISYMRB)
      KLAMHB = KLAMPB + NLAMDT
      KEND1  = KLAMHB + NLAMDT
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. NT2SQ(ISYMRB)) THEN
         CALL QUIT('Out of memory in CC3_XI (TOT_T3Y) ')
      ENDIF
C
C     Readin
C
      IOPT = 3
      CALL GET_T1_T2(IOPT,.TRUE.,WORK(KT1B),WORK(KT2B),LISTB,IDLSTB,
     *               ISYMRB,WORK(KEND1),LWRK1)
C
      IF (IPRINT .GT. 55) THEN
         XNORMVAL = DDOT(NT2SQ(ISYMRB),WORK(KT2B),1,WORK(KT2B),1)
         WRITE(LUPRI,*) 'Norm of T2B  ',XNORMVAL
         XNORMVAL = DDOT(NT1AM(ISYMRB),WORK(KT1B),1,WORK(KT1B),1)
         WRITE(LUPRI,*) 'Norm of T1B  ',XNORMVAL
      ENDIF
C
C     Integrals

      ISINTR1 = MULD2H(ISINT1,ISYMRB)
      ISINTR2 = MULD2H(ISINT2,ISYMRB)
C
      CALL CC3_BARINT(WORK(KT1B),ISYMRB,WORK(KLAMDP),
     *                WORK(KLAMDH),WORK(KEND1),LWRK1,
     *                LU3SRTR,FN3SRTR,LUCKJDR,FNCKJDR)
C
      CALL CC3_SORT1(WORK(KEND1),LWRK1,2,ISINTR1,LU3SRTR,FN3SRTR,
     *               LUDELDR,FNDELDR,IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *               IDUMMY,CDUMMY)
C
      CALL CC3_SINT(WORK(KLAMDH),WORK(KEND1),LWRK1,ISINTR1,
     *              LUDELDR,FNDELDR,LUDKBCR,FNDKBCR)

      CALL CC3_TCME(WORK(KLAMDP),ISINTR1,WORK(KEND1),LWRK1,
     *              IDUMMY,CDUMMY,LUDKBCR,FNDKBCR,
     *              IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *              IDUMMY,CDUMMY,LUDKBCR4,FNDKBCR4,2)

C---------------------------------------
C        Calculate the Lambda^B
C---------------------------------------
C
      CALL CCLR_LAMTRA(WORK(KLAMDP),WORK(KLAMPB),WORK(KLAMDH),
     *                 WORK(KLAMHB),WORK(KT1B),ISYMRB)
C
C------------------------------------------
C     Calculate the F^B matrix
C------------------------------------------
C
      IF ( .NOT.(LISTB(1:3).EQ.'RE ') ) THEN
         CALL GET_FOCKX(WORK(KFCKB),LABELB,IDLSTB,ISYMRB,WORK(KLAMDP),1,
     *                  WORK(KLAMDH),1,WORK(KEND1),LWRK1)
      END IF
C
C-----------------------------
C     Read occupied integrals.
C-----------------------------
C
C     Memory allocation.
C
C---------------------------------------------------------
C     Work allocation 
C---------------------------------------------------------
C
      KW3XOG1 = KEND0
      KT3OG2= KW3XOG1 + NTRAOC(ISINT2)
      KEND1  = KT3OG2+ NTRAOC(ISINT2)
      LWRK1  = LWORK  - KEND1
C
      KW3XOGX1 = KEND1
      KEND1   = KW3XOGX1 + NTRAOC(ISINTR2)
C
      KINTOC = KEND1
      KEND2  = KINTOC 
     *       + MAX(NTOTOC(ISINT2),NTOTOC(ISINT1),NTOTOC(ISYMOP))
      LWRK2  = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in CC3_XI')
      END IF
C
C     Occupied integrals.
C
C-----------------------
C     Read in integrals for SMAT etc.
C-----------------------
C
      CALL INTOCC_T30(LUCKJD,FNCKJD,WORK(KLAMDP),ISINT2,WORK(KW3XOG1),
     *                WORK(KT3OG2),WORK(KEND2),LWRK2)
C
C------------------------------------------
C     B transformed Occupied integrals.
C-----------------------------------------
C
      CALL INTOCC_T3X(LUCKJDR,FNCKJDR,WORK(KLAMDP),ISINTR2,
     *                WORK(KW3XOGX1),WORK(KEND2),LWRK2)
C
C----------------------------
C     General loop structure.
C----------------------------
C
      DO ISYMD = 1,NSYM

         ISYCKB  = MULD2H(ISYMOP,ISYMD)
         ISCKB2  = MULD2H(ISINT2,ISYMD)
         ISCKB2Y = MULD2H(ISINTR2,ISYMD)
         ISCKB1  = MULD2H(ISINT1,ISYMD)
C
         IF (IPRINT .GT. 55) THEN
C
            WRITE(LUPRI,*) 'In CC3_FT3B: ISYCKB :',ISYCKB
            WRITE(LUPRI,*) 'In CC3_FT3B: ISCKB2 :',ISCKB2
            WRITE(LUPRI,*) 'In CC3_FT3B: ISCKB2Y:',ISCKB2Y
            WRITE(LUPRI,*) 'In CC3_FT3B: ISCKB1 :',ISCKB1

         ENDIF

C
C--------------------------
C        Memory allocation.
C--------------------------
C
         KT3VDG3 = KEND1 
         KT3VDG2 = KT3VDG3 + NCKATR(ISCKB2)
         KEND2  = KT3VDG2 + NCKATR(ISCKB2)
         LWRK2  = LWORK  - KEND2

         KW3XVDG1  = KEND2
         KEND3   = KW3XVDG1  + NCKATR(ISCKB2)
         LWRK3   = LWORK  - KEND3
C
         KW3XVDGX1  = KEND3
         KEND3    = KW3XVDGX1  + NCKATR(ISCKB2Y)
         LWRK3    = LWORK    - KEND3

         KINTVI = KEND3
         KEND4  = KINTVI 
     *          + MAX(NCKA(ISCKB2),NCKA(ISCKB1),NCKA(ISYCKB))
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
C-----------------------------------------------
C        Read virtual integrals used in first s3am.
C-----------------------------------------------
C
            CALL INTVIR_T30_D(LUDKBC,FNDKBC,LUDELD,FNDELD,ISINT2,
     *                        WORK(KW3XVDG1),WORK(KT3VDG2),
     *                        WORK(KT3VDG3),WORK(KLAMDH),ISYMD,D,
     *                        WORK(KEND4),LWRK4)
C
C-----------------------------------------------------------------------
C        Read B transformed virtual integrals used for W for TOT_T3Y
C-----------------------------------------------------------------------
C
            IOFF = ICKBD(ISCKB2Y,ISYMD) + NCKATR(ISCKB2Y)*(D - 1) + 1
            IF (NCKATR(ISCKB2Y) .GT. 0) THEN
               CALL GETWA2(LUDKBCR,FNDKBCR,WORK(KW3XVDGX1),IOFF,
     &                     NCKATR(ISCKB2Y))
            ENDIF
C
C--------------------------------------------------------
C           Sum over ISYMB
C--------------------------------------------------------
C
            DO ISYMB = 1,NSYM

               ISYALJ  = MULD2H(ISYMB,ISYMT2)
               ISYALJ2 = MULD2H(ISYMD,ISYMT2)
               ISYMBD  = MULD2H(ISYMB,ISYMD)
               ISCKIJ  = MULD2H(ISYMBD,ISYMT3)
               ISCKD2  = MULD2H(ISINT2,ISYMB)
               ISCKD2Y = MULD2H(ISINTR2,ISYMB)
               ISWMAT  = MULD2H(ISYMRB,ISCKIJ)

               IF (IPRINT .GT. 55) THEN

                  WRITE(LUPRI,*) 'In CC3_FT3B: ISYMD  :',ISYMD
                  WRITE(LUPRI,*) 'In CC3_FT3B: ISYMB  :',ISYMB
                  WRITE(LUPRI,*) 'In CC3_FT3B: ISYALJ :',ISYALJ
                  WRITE(LUPRI,*) 'In CC3_FT3B: ISYALJ2:',ISYALJ2
                  WRITE(LUPRI,*) 'In CC3_FT3B: ISYMBD :',ISYMBD
                  WRITE(LUPRI,*) 'In CC3_FT3B: ISCKIJ :',ISCKIJ
                  WRITE(LUPRI,*) 'In CC3_FT3B: ISWMAT :',ISWMAT

               ENDIF

C           Can use kend3 since we do not need the integrals anymore.
               KSMAT   = KEND3
               KSMAT3  = KSMAT   + NCKIJ(ISCKIJ)
               KUMAT   = KSMAT3  + NCKIJ(ISCKIJ)
               KUMAT3  = KUMAT   + NCKIJ(ISCKIJ)
               KDIAG   = KUMAT3  + NCKIJ(ISCKIJ)
               KDIAGW  = KDIAG   + NCKIJ(ISCKIJ)
               KWMAT   = KDIAGW  + NCKIJ(ISWMAT)
               KINDSQW = KWMAT   + NCKIJ(ISWMAT)
               KINDSQ  = KINDSQW + (6*NCKIJ(ISWMAT) - 1)/IRAT + 1
               KINDEX  = KINDSQ  + (6*NCKIJ(ISCKIJ) - 1)/IRAT + 1
               KINDEX2 = KINDEX  + (NCKI(ISYALJ) - 1)/IRAT + 1
               KTMAT   = KINDEX2 + (NCKI(ISYALJ2) - 1)/IRAT + 1
               KW3XVBG1  = KTMAT   + MAX(NCKIJ(ISCKIJ),NCKIJ(ISWMAT))
               KT3VBG2  = KW3XVBG1  + NCKATR(ISCKD2)
               KT3VBG3 = KT3VBG2  + NCKATR(ISCKD2)
               KEND4   = KT3VBG3 + NCKATR(ISCKD2)
               LWRK4   = LWORK   - KEND4
C
               KW3XVDGX2 = KEND4
               KEND4   = KW3XVDGX2 + NCKATR(ISCKD2Y)

C
               KINTVI  = KEND4
               KEND5   = KINTVI  + NCKA(ISCKD2)
               LWRK5   = LWORK   - KEND5
C
               IF (LWRK5 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Memory available : ',LWORK
                  WRITE(LUPRI,*) 'Memory needed    : ',KEND5
                  CALL QUIT('Insufficient space in CC3_XI')
               END IF
C
C---------------------------------------------
C           Construct part of the diagonal.
C---------------------------------------------
C
               CALL CC3_DIAG(WORK(KDIAG),WORK(KFOCKD),ISCKIJ)
               CALL CC3_DIAG(WORK(KDIAGW),WORK(KFOCKD),ISWMAT)
C
               IF (IPRINT .GT. 55) THEN
                  XNORMVAL = DDOT(NCKIJ(ISCKIJ),WORK(KDIAG),1,
     *                    WORK(KDIAG),1)
                  WRITE(LUPRI,*) 'Norm of DIA  ',XNORMVAL
               ENDIF
C
C-------------------------------------
C           Construct index arrays.
C----------------------------------
C
               LENSQ = NCKIJ(ISCKIJ)
               CALL CC3_INDSQ(WORK(KINDSQ),LENSQ,ISCKIJ)
               CALL CC3_INDEX(WORK(KINDEX),ISYALJ)
               CALL CC3_INDEX(WORK(KINDEX2),ISYALJ2)
               LENSQW = NCKIJ(ISWMAT)
               CALL CC3_INDSQ(WORK(KINDSQW),LENSQW,ISWMAT) 

               DO B = 1,NVIR(ISYMB)
C
C------------------------------------
C           Initialise
C---------------------------------------
C
               CALL DZERO(WORK(KWMAT),NCKIJ(ISWMAT))
C
C-----------------------------------------------
C           Read virtual integrals used in second s3am.
C-----------------------------------------------
C
                  CALL INTVIR_T30_D(LUDKBC,FNDKBC,LUDELD,FNDELD,ISINT2,
     *                              WORK(KW3XVBG1),WORK(KT3VBG2),
     *                              WORK(KT3VBG3),WORK(KLAMDH),ISYMB,B,
     *                              WORK(KEND5),LWRK5)
C
C--------------------------------------------------------------------
C           Read B transformed virtual integrals used in W
C--------------------------------------------------------------------
C
                  IOFF = ICKBD(ISCKD2Y,ISYMB) + 
     &                   NCKATR(ISCKD2Y)*(B - 1) + 1
                  IF (NCKATR(ISCKD2Y) .GT. 0) THEN
                     CALL GETWA2(LUDKBCR,FNDKBCR,WORK(KW3XVDGX2),IOFF,
     &                           NCKATR(ISCKD2Y))
                  ENDIF
C
C------------------------------------------------------------------------
C           Calculate the S(ci,bk,dj) matrix for T3 for B,D.
C-------------------------------------------------------------------
C
                  CALL GET_T30_BD(ISYMT3,ISINT2,WORK(KT2TP),ISYMT2,
     *                            WORK(KTMAT),WORK(KFOCKD),WORK(KDIAG),
     *                            WORK(KINDSQ),LENSQ,WORK(KSMAT),
     *                            WORK(KW3XVDG1),WORK(KT3VDG2),
     *                            WORK(KW3XOG1),WORK(KINDEX),
     *                            WORK(KSMAT3),WORK(KW3XVBG1),
     *                            WORK(KT3VBG2),WORK(KINDEX2),
     *                            WORK(KUMAT),WORK(KT3VDG3),
     *                            WORK(KT3OG2),WORK(KUMAT3),
     *                            WORK(KT3VBG3),ISYMB,B,ISYMD,D,
     *                            ISCKIJ,WORK(KEND4),LWRK4)
C
C------------------------------------------------------
C Based on T3 BD amplitudes calculate W BD intermediates
C------------------------------------------------------
C

                  IF (LISTB(1:3).EQ.'R1 ') THEN

                     CALL GET_T3X_BD(ISYMRB,WORK(KWMAT),ISWMAT,
     *                               WORK(KTMAT),ISCKIJ,WORK(KFCKB),
     *                               ISYMRB,WORK(KT2TP),ISYMT2,
     *                               WORK(KT2B),ISYMRB,WORK(KW3XVDG1),
     *                               WORK(KW3XVBG1),WORK(KW3XOG1),
     *                               ISINT2,WORK(KW3XVDGX1),
     *                               WORK(KW3XVDGX2),WORK(KW3XOGX1),
     *                               ISINTR2,FREQB,WORK(KDIAGW),
     *                               WORK(KFOCKD),WORK(KINDSQW),LENSQW,
     *                               B,ISYMB,D,ISYMD,WORK(KEND4),LWRK4)

                  END IF
C
                  IF (LISTB(1:3).EQ.'RE ') THEN
                     CALL WBD_GROUND(WORK(KT2B),ISYMRB,WORK(KTMAT),
     *                               WORK(KW3XVDG1),WORK(KW3XVBG1),
     *                               WORK(KW3XOG1),ISINT2,WORK(KWMAT),
     *                               WORK(KEND4),LWRK4,
     *                               WORK(KINDSQW),LENSQW,
     *                               ISYMB,B,ISYMD,D)
C
                     CALL WBD_GROUND(WORK(KT2TP),ISYMT2,WORK(KTMAT),
     *                               WORK(KW3XVDGX1),WORK(KW3XVDGX2),
     *                               WORK(KW3XOGX1),ISINTR2,WORK(KWMAT),
     *                               WORK(KEND4),LWRK4,
     *                               WORK(KINDSQW),LENSQW,
     *                               ISYMB,B,ISYMD,D)

C
C------------------------------------------------
C     Divide by the energy difference and
C     remove the forbidden elements
C------------------------------------------------
C
                     CALL WBD_DIA(B,ISYMB,D,ISYMD,FREQB,
     *                            ISWMAT,WORK(KWMAT),
     *                            WORK(KDIAGW),WORK(KFOCKD))

                     CALL T3_FORBIDDEN(WORK(KWMAT),ISYMRB,
     *                                 ISYMB,B,ISYMD,D)
                  END IF
C
C----------------------------------------------------------------------
C     Calculate the  term <L2|[H,E_ai],W^BD(3)]]|HF> 
C----------------------------------------------------------------------
C 
*      call sum_pt3(work(kwmat),isymb,b,isymd,d,
*    *             iswmat,work(kx3am),4)
C
                  CALL T3_ONEL1W(OMEGA1,WORK(KWMAT),WORK(KTMAT),
     *                           ISWMAT,WORK(KXIAJB),WORK(KXAIJB),
     *                           WORK(KYIAJB),WORK(KYAIJB),ISINT3,
     *                           WORK(KL2TP),ISYML2,WORK(KINDSQW),
     *                           LENSQW,WORK(KEND4),
     *                           LWRK4,ISYMB,B,ISYMD,D)
C
                  CALL T3_ONEL2W(OMEGA1,WORK(KWMAT),WORK(KTMAT),
     *                           ISWMAT,WORK(KXIAJB),WORK(KYIAJB),
     *                           ISINT3,
     *                           WORK(KL2TP),ISYML2,WORK(KINDSQW),
     *                           LENSQW,WORK(KEND4),
     *                           LWRK4,ISYMB,B,ISYMD,D)
C 
                  CALL T3_ONEL3W(OMEGA1,WORK(KWMAT),WORK(KTMAT),
     *                           ISWMAT,WORK(KYIAJB),WORK(KYAIBJ),
     *                           ISINT3,
     *                           WORK(KL2TP),ISYML2,WORK(KINDSQW),
     *                           LENSQW,WORK(KEND4),
     *                           LWRK4,ISYMB,B,ISYMD,D)
C
               END DO  ! B
            END DO     ! ISYMB
         END DO        ! D
      END DO           ! ISYMD
C
*     write(lupri,*)'omega1 '
*     do i = 1,nt1am(isyres)
*        write(lupri,*) i , OMEGA1(i)
*     end do
*
*      write(lupri,*) 'cont to t3x in CC3_FT3B'
*      call print_pt3(work(kx3am),1,4)
*      write(lupri,*) 'The summed S terms : '
*      call print_pt3(work(kx3am),1,1)
C
C--------------------------------
C     Close files
C--------------------------------
C
      CALL WCLOSE2(LUCKJD,FNCKJD,'KEEP')
      CALL WCLOSE2(LUTOC,FNTOC,'KEEP')
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
      CALL WCLOSE2(LUDKBCR4,FNDKBCR4,'DELETE')
C
C-------------
C     End
C-------------
C
      CALL QEXIT('CC3_FT3B')
C
      RETURN
      END
C  /* Deck t3_onel1w */
      SUBROUTINE T3_ONEL1W(OMEGA1,WMAT,TMAT,ISYMIM,XIAJB,XAIJB,YIAJB, 
     *                    YAIJB,ISYINT,C2TP,ISYMC2,INDSQ,LENSQ,WORK,
     *                    LWORK,ISYMB,B,ISYMD,D)
C
C     Written by Filip Pawlowski and Poul Jorgensen based on
C     SUBROUTINE T3_ONEL1  by K. Hald, Spring 2002.
C
C     Calculate the term t^{def}_{lmn} T^{fd}_{mi} g_{nela} 
C                      - t^{def}_{lnm} T^{fd}_{mi} L_{nela}
C     using W intermediates
C
C     g contributions
C
C
C 1)
C     W^{ed}(flnm) T^{df}_{mi} g_{nela}
C 2)
C     W^{ed}(fnlm) T^{fd}_{mi} g_{nela}
C 3)
C     W^{fd}(emln) T^{fd}_{mi} g_{nela}
C 
C     L contributions
C
C 4)
C     W^{ed}(flmn) T^{df}_{mi} L_{nela}
C 5)
C     W^{ed}(fmln) T^{fd}_{mi} L_{nela}
C 6)
C     W^{fd}(enlm) T^{fd}_{mi} L_{nela}
C
C
C     XIAJB contains g and YIAJB contains L
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYMIM, ISYINT, ISYMC2, LENSQ, LWORK, ISYMB, ISYMD
      INTEGER INDSQ(LENSQ,6), INDEX
      INTEGER ISYRE1, ISYRES, ISYMBD, ISFLMN, ISYANL, LENGTH
      INTEGER ISYFIM, KTMAT, KC2TEMP, KINT, KEND1, LWRK1
      INTEGER ISYMM, ISYMFI, ISYMF, KOFF1, KOFF2
      INTEGER ISYML, ISYMAN, ISYMA, ISYMN, ISYMLN, NBN, NAN, NAL
      INTEGER ISYMI, ISYMAB, ISYMFM, KOFF3, NUMBFM, NUMBLN, NUMBA
      INTEGER ISYMBN, ISYMAL, KINT2
      INTEGER ISYMMI, KMIMAT, KEND2, LWRK2, ISYBMI, ISYMBM
      INTEGER ISENLM, ISYENL, NUMENL, NUMM, NUMA
      INTEGER ISYRE2
C
#if defined (SYS_CRAY)
      REAL OMEGA1(*), WMAT(*),  TMAT(*), XIAJB(*), XAIJB(*), YAIJB(*)
      REAL YIAJB(*), C2TP(*), WORK(LWORK), ZERO, ONE
#else
      DOUBLE PRECISION OMEGA1(*), WMAT(*),  TMAT(*), XIAJB(*), XAIJB(*)
      DOUBLE PRECISION YAIJB(*)
      DOUBLE PRECISION YIAJB(*), C2TP(*), WORK(LWORK), ZERO, ONE
#endif
C
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('T3_ONEL1W')
C
      ISYMBD = MULD2H(ISYMB,ISYMD)
C
      ISYRE1 = MULD2H(ISYMIM,ISYMC2)
      ISYRE2 = MULD2H(ISYRE1,ISYMBD)
      ISYRES = MULD2H(ISYRE2,ISYINT)
C
      ISFLMN = ISYMIM
      ISYANL = MULD2H(ISYMB,ISYINT)
C
      LENGTH = NCKIJ(ISFLMN)
C
C calculate the terms
C
C 1)
C     W^{ed}(flnm) T^{df}_{mi} g_{nela}
C 2)
C     W^{ed}(fnlm) T^{fd}_{mi} g_{nela}
C
C-----------------------------
C     Sort C2
C-----------------------------
C
      ISYFIM = MULD2H(ISYMC2,ISYMD)
C
      KTMAT   = 1
      KC2TEMP = KTMAT   + NCKIJ(ISFLMN)
      KINT    = KC2TEMP + NMAIJA(ISYFIM)
      KINT2   = KINT    + NCKI(ISYANL)
      KEND1   = KINT2   + NCKI(ISYANL)
      LWRK1   = LWORK   - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory in T3_ONEL1W (sort)')
      ENDIF
C
      DO ISYMM = 1, NSYM
         ISYMFI = MULD2H(ISYFIM,ISYMM)
         DO ISYMF = 1, NSYM
            ISYMI = MULD2H(ISYMFI,ISYMF)
            ISYMFM = MULD2H(ISYMF,ISYMM)
C
            DO M = 1, NRHF(ISYMM)
               DO I = 1, NRHF(ISYMI)
C
                  KOFF1 = IT2SP(ISYFIM,ISYMD)
     *                  + NCKI(ISYFIM)*(D-1)
     *                  + ICKI(ISYMFI,ISYMM)
     *                  + NT1AM(ISYMFI)*(M-1)
     *                  + IT1AM(ISYMF,ISYMI)
     *                  + NVIR(ISYMF)*(I-1) 
     *                  + 1
C
                  KOFF2 = KC2TEMP
     *                  + ICKI(ISYMFM,ISYMI)
     *                  + NT1AM(ISYMFM)*(I-1)
     *                  + IT1AM(ISYMF,ISYMM)
     *                  + NVIR(ISYMF)*(M-1) 

C
                  CALL DCOPY(NVIR(ISYMF),C2TP(KOFF1),1,WORK(KOFF2),1)
C
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C---------------------------
C     Sort integrals.
C---------------------------
C
      DO ISYML = 1, NSYM
         ISYMAN = MULD2H(ISYANL,ISYML)
         DO ISYMA = 1, NSYM
            ISYMN  = MULD2H(ISYMAN,ISYMA)
            ISYMLN = MULD2H(ISYMN,ISYML)
            ISYMBN = MULD2H(ISYMB,ISYMN)
            ISYMAL = MULD2H(ISYMA,ISYML)
C
            DO N = 1, NRHF(ISYMN)
               NBN = IT1AM(ISYMB,ISYMN) + NVIR(ISYMB)*(N-1) + B
               DO A = 1, NVIR(ISYMA)
                  NAN = IT1AM(ISYMA,ISYMN) + NVIR(ISYMA)*(N-1) + A
                  DO L = 1, NRHF(ISYML)
                     NAL = IT1AM(ISYMA,ISYML) + NVIR(ISYMA)*(L-1) + A
C
                     KOFF1 = IT2AM(ISYMBN,ISYMAL) + INDEX(NBN,NAL)
                     KOFF2 = KINT - 1
     *                     + IMAIJA(ISYMLN,ISYMA)
     *                     + NMATIJ(ISYMLN)*(A-1)
     *                     + IMATIJ(ISYML,ISYMN)
     *                     + NRHF(ISYML)*(N-1)
     *                     + L
                     KOFF3 = KINT2 - 1
     *                     + IMAIJA(ISYMLN,ISYMA)
     *                     + NMATIJ(ISYMLN)*(A-1)
     *                     + IMATIJ(ISYML,ISYMN)
     *                     + NRHF(ISYML)*(N-1)
     *                     + L
C
                     WORK(KOFF2) = XIAJB(KOFF1)
                     WORK(KOFF3) = YIAJB(KOFF1)
C
                  ENDDO
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C
C----------------------
C     Construct TMAT for the g term
C----------------------
C
      DO I = 1, LENGTH
         TMAT(I) = WMAT(INDSQ(I,2))
         WORK(KTMAT-1+I) = WMAT(INDSQ(I,5))
      ENDDO
C
C---------------------------------------------
C     Symmetry sorting if symmetry
C---------------------------------------------
C
      IF (NSYM .GT. 1) THEN
         CALL CC_GATHER(LENGTH,WORK(KEND1),TMAT,INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK(KEND1),1,TMAT,1)
C
         CALL CC_GATHER(LENGTH,WORK(KEND1),WORK(KTMAT),INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK(KEND1),1,WORK(KTMAT),1)
      ENDIF
C
C-------------------------------------
C     Contract
C-------------------------------------
C
      DO ISYMA = 1, NSYM
         ISYMI = MULD2H(ISYRES,ISYMA)
         ISYMAB = MULD2H(ISYMA,ISYMB)
         ISYMLN = MULD2H(ISYINT,ISYMAB)
         ISYMFM = MULD2H(ISYFIM,ISYMI)
C
         CALL DZERO(WORK(KEND1),NMATIJ(ISYMLN)*NRHF(ISYMI))
C
         KOFF1 = ISAIKL(ISYMFM,ISYMLN) + 1
         KOFF2 = KC2TEMP
     *         + ICKI(ISYMFM,ISYMI)
         KOFF3 = KEND1
C
         NUMBFM = MAX(1,NT1AM(ISYMFM))
         NUMBLN = MAX(1,NMATIJ(ISYMLN))
C
         CALL DGEMM('T','N',NMATIJ(ISYMLN),NRHF(ISYMI),NT1AM(ISYMFM),
     *              ONE,TMAT(KOFF1),NUMBFM,WORK(KOFF2),NUMBFM,
     *              ONE,WORK(KOFF3),NUMBLN)
C
         KOFF1 = KTMAT
     *         + ISAIKL(ISYMFM,ISYMLN)
         KOFF2 = IT2SP(ISYFIM,ISYMD)
     *         + NCKI(ISYFIM)*(D-1)
     *         + ICKI(ISYMFM,ISYMI) + 1
         KOFF3 = KEND1
C
         NUMBFM = MAX(1,NT1AM(ISYMFM))
         NUMBLN = MAX(1,NMATIJ(ISYMLN))
C
         CALL DGEMM('T','N',NMATIJ(ISYMLN),NRHF(ISYMI),NT1AM(ISYMFM),
     *              ONE,WORK(KOFF1),NUMBFM,C2TP(KOFF2),NUMBFM,
     *              ONE,WORK(KOFF3),NUMBLN)
C
         KOFF1 = KINT
     *         + IMAIJA(ISYMLN,ISYMA) 
         KOFF2 = KEND1
         KOFF3 = IT1AM(ISYMA,ISYMI) + 1
C
         NUMBLN = MAX(1,NMATIJ(ISYMLN))
         NUMBA  = MAX(1,NVIR(ISYMA))
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NMATIJ(ISYMLN),
     *              -ONE,WORK(KOFF1),NUMBLN,WORK(KOFF2),NUMBLN,
     *              ONE,OMEGA1(KOFF3),NUMBA)
C
      ENDDO
C
C-------------------------------------
C 3)
C     W^{fd}(emln) T^{fd}_{mi} g_{nela}
C--------------------------------------
C
C -------------------------------
C     Sort T^{BD}_{mi} as T_{mi}
C -------------------------------
C
      ISYMMI = MULD2H(ISYMBD,ISYMC2)
      ISYBMI = MULD2H(ISYMMI,ISYMB)
C    
      KMIMAT  = KEND1
      KEND2   = KMIMAT + NMATIJ(ISYMMI)
      LWRK2   = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         CALL QUIT('Out of memory in T3_ONEL1W (sort)')
      ENDIF
C
      DO ISYMI = 1,NSYM
         ISYMM = MULD2H(ISYMMI,ISYMI)
         ISYMBM = MULD2H(ISYBMI,ISYMI)
         DO I = 1,NRHF(ISYMI)
            DO M = 1,NRHF(ISYMM)
               KOFF1 = IT2SP(ISYBMI,ISYMD)
     *                  + NCKI(ISYBMI)*(D-1)
     *                  + ICKI(ISYMBM,ISYMI)
     *                  + NT1AM(ISYMBM)*(I-1)
     *                  + IT1AM(ISYMB,ISYMM)
     *                  + NVIR(ISYMB)*(M-1) 
     *                  + B
               KOFF2 = IMATIJ(ISYMM,ISYMI)
     *                  + NRHF(ISYMM)*(I-1)
     *                  + M 
               WORK(KMIMAT-1+KOFF2) = C2TP(KOFF1)
            ENDDO 
         ENDDO 
      ENDDO 
C
C
C----------------------------------
C     Construct TMAT 
C----------------------------------
C
      DO I = 1, LENGTH
         TMAT(I) =  WMAT(INDSQ(I,5))
      ENDDO
C
C-------------------------------------
C     Omega(ai) = Omega(ai) +
C     W^{BD}(emln) T^{BD}_{mi} g_{nela} 
C
C     Tmat(enl,m) * T_{mi} * Xaijb(enla)
C--------------------------------------
C
      ISENLM = ISYMIM
      DO ISYMA = 1,NSYM
         ISYMI = MULD2H(ISYRES,ISYMA)
         ISYMMI = MULD2H(ISYMBD,ISYMC2)
         ISYMM  = MULD2H(ISYMMI,ISYMI)
         ISYENL = MULD2H(ISENLM,ISYMM)
C
         KOFF1  = ISAIKJ(ISYENL,ISYMM) + 1
         KOFF2  = KMIMAT + IMATIJ(ISYMM,ISYMI)  
         KOFF3  = KEND2  + ISAIKJ(ISYENL,ISYMI) 
         NUMENL = MAX(1,NCKI(ISYENL))
         NUMM   = MAX(1,NRHF(ISYMM))
C
         CALL DGEMM('N','N',NCKI(ISYENL),NRHF(ISYMI),NRHF(ISYMM),
     *              ONE,TMAT(KOFF1),NUMENL,WORK(KOFF2),NUMM,
     *              ZERO,WORK(KOFF3),NUMENL) 
C
         KOFF1  = IT2SP(ISYENL,ISYMA)  + 1
         KOFF2  = KEND2 + ISAIKJ(ISYENL,ISYMI) 
         KOFF3  = IT1AM(ISYMA,ISYMI)   + 1 
         NUMENL = MAX(1,NCKI(ISYENL))
         NUMA   = MAX(1,NVIR(ISYMA)) 
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NCKI(ISYENL),
     *              -ONE,XAIJB(KOFF1),NUMENL,WORK(KOFF2),NUMENL,
     *              ONE,OMEGA1(KOFF3),NUMA)
      END DO
C
C-------------------------------------
C 4)
C     W^{ed}(flmn) T^{df}_{mi} L_{nela}
C 5)
C     W^{ed}(fmln) T^{fd}_{mi} L_{nela}
C-------------------------------------
C
C----------------------------------
C     Construct TMAT for L term
C----------------------------------
C
      DO I = 1, LENGTH
         TMAT(I) =  WMAT(INDSQ(I,1))
C
         WORK(KTMAT-1+I) =  WMAT(I)
      ENDDO
C
C---------------------------------------------
C     Symmetry sorting if symmetry
C---------------------------------------------
C
      IF (NSYM .GT. 1) THEN
         CALL CC_GATHER(LENGTH,WORK(KEND2),TMAT,INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK(KEND2),1,TMAT,1)
C
         CALL CC_GATHER(LENGTH,WORK(KEND2),WORK(KTMAT),INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK(KEND2),1,WORK(KTMAT),1)
      ENDIF
C
C-------------------------------------
C     Contract
C-------------------------------------
C
      DO ISYMA = 1, NSYM
         ISYMI = MULD2H(ISYRES,ISYMA)
         ISYMAB = MULD2H(ISYMA,ISYMB)
         ISYMLN = MULD2H(ISYINT,ISYMAB)
         ISYMFM = MULD2H(ISYFIM,ISYMI)
C
         CALL DZERO(WORK(KEND2),NMATIJ(ISYMLN)*NRHF(ISYMI))
C
         KOFF1 = ISAIKL(ISYMFM,ISYMLN) + 1
         KOFF2 = KC2TEMP
     *         + ICKI(ISYMFM,ISYMI)
         KOFF3 = KEND2
C
         NUMBFM = MAX(1,NT1AM(ISYMFM))
         NUMBLN = MAX(1,NMATIJ(ISYMLN))
C
         CALL DGEMM('T','N',NMATIJ(ISYMLN),NRHF(ISYMI),NT1AM(ISYMFM),
     *              ONE,TMAT(KOFF1),NUMBFM,WORK(KOFF2),NUMBFM,
     *              ONE,WORK(KOFF3),NUMBLN)
C
         KOFF1 = KTMAT
     *         + ISAIKL(ISYMFM,ISYMLN)
         KOFF2 = IT2SP(ISYFIM,ISYMD)
     *         + NCKI(ISYFIM)*(D-1)
     *         + ICKI(ISYMFM,ISYMI) + 1
         KOFF3 = KEND2
C
         NUMBFM = MAX(1,NT1AM(ISYMFM))
         NUMBLN = MAX(1,NMATIJ(ISYMLN))
C
         CALL DGEMM('T','N',NMATIJ(ISYMLN),NRHF(ISYMI),NT1AM(ISYMFM),
     *              ONE,WORK(KOFF1),NUMBFM,C2TP(KOFF2),NUMBFM,
     *              ONE,WORK(KOFF3),NUMBLN)
C
         KOFF1 = KINT2
     *         + IMAIJA(ISYMLN,ISYMA)
         KOFF2 = KEND2
         KOFF3 = IT1AM(ISYMA,ISYMI) + 1
C
         NUMBLN = MAX(1,NMATIJ(ISYMLN))
         NUMBA  = MAX(1,NVIR(ISYMA))
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NMATIJ(ISYMLN),
     *              ONE,WORK(KOFF1),NUMBLN,WORK(KOFF2),NUMBLN,
     *              ONE,OMEGA1(KOFF3),NUMBA)
      ENDDO
C
C--------------------------------------
C 6)
C     W^{fd}(enlm) T^{fd}_{mi} L_{nela}
C--------------------------------------
C
C----------------------------------
C     Construct TMAT
C----------------------------------
C
      DO I = 1, LENGTH
         TMAT(I) =  WMAT(I)
      ENDDO
C
C-------------------------------------
C     Omega(ai) = Omega(ai) +
C     W^{BD}(emln) T^{BD}_{mi} L_{nela}
C
C     Tmat(enl,m) * T_{mi} * Yaijb(enla)
C--------------------------------------
C

      ISENLM = ISYMIM
      DO ISYMA = 1,NSYM
         ISYMI = MULD2H(ISYRES,ISYMA)
         ISYMMI = MULD2H(ISYMBD,ISYMC2)
         ISYMM  = MULD2H(ISYMMI,ISYMI)
         ISYENL = MULD2H(ISENLM,ISYMM)
C
         KOFF1  = ISAIKJ(ISYENL,ISYMM) + 1
         KOFF2  = KMIMAT + IMATIJ(ISYMM,ISYMI)
         KOFF3  = KEND2  + ISAIKJ(ISYENL,ISYMI)
         NUMENL = MAX(1,NCKI(ISYENL))
         NUMM   = MAX(1,NRHF(ISYMM))
C
         CALL DGEMM('N','N',NCKI(ISYENL),NRHF(ISYMI),NRHF(ISYMM),
     *              ONE,TMAT(KOFF1),NUMENL,WORK(KOFF2),NUMM,
     *              ZERO,WORK(KOFF3),NUMENL)
C
         KOFF1  = IT2SP(ISYENL,ISYMA)  + 1
         KOFF2  = KEND2 + ISAIKJ(ISYENL,ISYMI)
         KOFF3  = IT1AM(ISYMA,ISYMI)   + 1
         NUMENL = MAX(1,NCKI(ISYENL))
         NUMA   = MAX(1,NVIR(ISYMA))
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NCKI(ISYENL),
     *              ONE,YAIJB(KOFF1),NUMENL,WORK(KOFF2),NUMENL,
     *              ONE,OMEGA1(KOFF3),NUMA)
C
      END DO

C

C
C----------------------------
C     End.
C----------------------------
C
      CALL QEXIT('T3_ONEL1W')
C
      RETURN
      END
C  /* Deck t3_onel2w */
      SUBROUTINE T3_ONEL2W(OMEGA1,WMAT,TMAT,ISYMIM,XIAJB,YIAJB,
     *                    ISYINT,C2TP,ISYMC2,INDSQ,LENSQ,WORK,LWORK,
     *                    ISYMB,B,ISYMD,D)
C
C     Written Filip Pawlowski and Poul Jorgensen based on
C     SUBROUTINE T3_ONEL2 by K. Hald, Spring 2002.
C
C     Calculate the term t^{def}_{lmn} L^{ad}_{mn} g_{ielf}
C                       -t^{def}_{nml} L^{ad}_{mn} L_{ielf}
C     using W intermediates
C
C     g contributions
C 1)
C     W^{df}(enml) T^{ad}_{mn} g_{ifle}
C 2)
C     W^{df}(emnl) T^{ad}_{mn} g_{ielf}
C 3)
C     W^{ef}(dlnm) T^{ad}_{mn} g_{ielf}
C
C     L contributions
C 4)
C     W^{df}(elmn) T^{ad}_{mn} L_{ifle}
C 5)
C     W^{df}(emln) T^{ad}_{mn} L_{ielf}
C 6)
C     W^{ef}(dnlm) T^{ad}_{mn} L_{ielf}


C     XIAJB contains g and YIAJB contains L
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYMIM, ISYINT, ISYMC2, LENSQ, LWORK, ISYMB, ISYMD
      INTEGER INDSQ(LENSQ,6), INDEX
      INTEGER ISYRE1, ISYRES, ISYMBD, ISELMN, ISYAMN, ISYELI
      INTEGER LENGTH, KTMAT, KINT1, KINT2, KC2TEMP, KEND1, LWRK1
      INTEGER ISYMN, ISYMAM, ISYMA, ISYMM, ISYMMN, KOFF1, KOFF2, KOFF3
      INTEGER ISYML, ISYMEI, ISYMDL, ISYME, ISYMI, ISYMEL, NDL
      INTEGER NEI, NUMBEL, NUMBMN, NUMBA
      INTEGER KINTIL, ISYMIL, ISYMBI, NBI, ISYENM, NUMENM, NUML, NUMA
      INTEGER ISYRE2
C
#if defined (SYS_CRAY)
      REAL OMEGA1(*), WMAT(*),  TMAT(*), XIAJB(*)
      REAL YIAJB(*), C2TP(*), WORK(LWORK), ZERO, ONE
#else
      DOUBLE PRECISION OMEGA1(*), WMAT(*),  TMAT(*), XIAJB(*)
      DOUBLE PRECISION YIAJB(*), C2TP(*), WORK(LWORK), ZERO, ONE
#endif
C
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('T3_ONEL2W')
C
      ISYMBD = MULD2H(ISYMB,ISYMD)
C
      ISYRE1 = MULD2H(ISYMIM,ISYINT)
      ISYRE2 = MULD2H(ISYMBD,ISYRE1)
      ISYRES = MULD2H(ISYRE2,ISYMC2)
C
      ISELMN = ISYMIM
      ISYMIL = MULD2H(ISYMBD,ISYINT)
C
      ISYAMN = MULD2H(ISYMB,ISYMC2)
      ISYELI = MULD2H(ISYMD,ISYINT)
C
      LENGTH = NCKIJ(ISELMN)
C
      KTMAT   = 1
      KINT1   = KTMAT   + NCKIJ(ISELMN)
      KINT2   = KINT1   + NCKI(ISYELI)
      KC2TEMP = KINT2   + NCKI(ISYELI)
      KINTIL  = KC2TEMP + NMAIJA(ISYAMN)
      KEND1   = KINTIL  + NMATIJ(ISYMIL)
      LWRK1   = LWORK   - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory in T3_ONEL2W (sort)')
      ENDIF
C
C     Calculate the terms
C 1)
C     W^{df}(enml) T^{ad}_{mn} g_{ifle}
C 2)
C     W^{df}(emnl) T^{ad}_{mn} g_{ielf}
C
C-----------------------------
C     Sort C2
C-----------------------------
C
      DO ISYMN = 1, NSYM
         ISYMAM = MULD2H(ISYAMN,ISYMN)
         DO ISYMA = 1, NSYM
            ISYMM  = MULD2H(ISYMAM,ISYMA)
            ISYMMN = MULD2H(ISYMM,ISYMN)
C
            DO M = 1, NRHF(ISYMM)
               DO N = 1, NRHF(ISYMN)
C
                  KOFF1 = IT2SP(ISYAMN,ISYMB)
     *                  + NCKI(ISYAMN)*(B-1)
     *                  + ICKI(ISYMAM,ISYMN)
     *                  + NT1AM(ISYMAM)*(N-1)
     *                  + IT1AM(ISYMA,ISYMM)
     *                  + NVIR(ISYMA)*(M-1) 
     *                  + 1
C
                  KOFF2 = KC2TEMP - 1
     *                  + IMAIJA(ISYMMN,ISYMA)
     *                  + IMATIJ(ISYMM,ISYMN)
     *                  + NRHF(ISYMM)*(N-1)
     *                  + M

C
                  CALL DCOPY(NVIR(ISYMA),C2TP(KOFF1),1,
     *                       WORK(KOFF2),NMATIJ(ISYMMN))
C
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C---------------------------
C     Sort g integrals.
C---------------------------
C
      DO ISYML = 1, NSYM
         ISYMEI = MULD2H(ISYELI,ISYML)
         ISYMDL = MULD2H(ISYML,ISYMD)
         DO ISYME = 1, NSYM
            ISYMI  = MULD2H(ISYMEI,ISYME)
            ISYMEL = MULD2H(ISYME,ISYML)
C
            DO L = 1, NRHF(ISYML)
               NDL = IT1AM(ISYMD,ISYML) + NVIR(ISYMD)*(L-1) + D
               DO E = 1, NVIR(ISYME)
                  DO I = 1, NRHF(ISYMI)
                     NEI = IT1AM(ISYME,ISYMI) + NVIR(ISYME)*(I-1) + E
C
                     KOFF1 = IT2AM(ISYMDL,ISYMEI) + INDEX(NDL,NEI)
                     KOFF2 = KINT1 - 1
     *                     + ICKI(ISYMEL,ISYMI)
     *                     + NT1AM(ISYMEL)*(I-1)
     *                     + IT1AM(ISYME,ISYML)
     *                     + NVIR(ISYME)*(L-1)
     *                     + E
                     KOFF3 = KINT2 - 1
     *                     + ICKI(ISYMEI,ISYML)
     *                     + NT1AM(ISYMEI)*(L-1)
     *                     + IT1AM(ISYME,ISYMI)
     *                     + NVIR(ISYME)*(I-1)
     *                     + E
C
                     WORK(KOFF2) = XIAJB(KOFF1)
                     WORK(KOFF3) = XIAJB(KOFF1)
C
                  ENDDO
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C
C----------------------
C     Construct TMAT
C----------------------
C
      DO I = 1, LENGTH
         TMAT(I) = WMAT(INDSQ(I,2))
C
      WORK(KTMAT-1+I) = WMAT(INDSQ(I,5))
      ENDDO
C
C---------------------------------------------
C     Symmetry sorting if symmetry
C---------------------------------------------
C
      IF (NSYM .GT. 1) THEN
         CALL CC_GATHER(LENGTH,WORK(KEND1),TMAT,INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK(KEND1),1,TMAT,1)
C
         CALL CC_GATHER(LENGTH,WORK(KEND1),WORK(KTMAT),INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK(KEND1),1,WORK(KTMAT),1)
      ENDIF
C
C-------------------------------------
C     Contract
C-------------------------------------
C
      DO ISYMA = 1, NSYM
         ISYMI  = MULD2H(ISYRES,ISYMA)
         ISYMEL = MULD2H(ISYELI,ISYMI)
         ISYMMN = MULD2H(ISYMEL,ISELMN)
C
         KOFF1 = ISAIKL(ISYMEL,ISYMMN) + 1
         KOFF2 = KINT1
     *         + ICKI(ISYMEL,ISYMI)
         KOFF3 = KEND1
C
         NUMBEL = MAX(1,NT1AM(ISYMEL))
         NUMBMN = MAX(1,NMATIJ(ISYMMN))
C
         CALL DGEMM('T','N',NMATIJ(ISYMMN),NRHF(ISYMI),NT1AM(ISYMEL),
     *              ONE,TMAT(KOFF1),NUMBEL,WORK(KOFF2),NUMBEL,
     *              ZERO,WORK(KOFF3),NUMBMN)
C
         KOFF1 = KTMAT
     *         + ISAIKL(ISYMEL,ISYMMN)
         KOFF2 = KINT2
     *         + ICKI(ISYMEL,ISYMI)
         KOFF3 = KEND1
C
         NUMBEL = MAX(1,NT1AM(ISYMEL))
         NUMBMN = MAX(1,NMATIJ(ISYMMN))
C
         CALL DGEMM('T','N',NMATIJ(ISYMMN),NRHF(ISYMI),NT1AM(ISYMEL),
     *              ONE,WORK(KOFF1),NUMBEL,WORK(KOFF2),NUMBEL,
     *              ONE,WORK(KOFF3),NUMBMN)
C
         KOFF1 = KC2TEMP
     *         + IMAIJA(ISYMMN,ISYMA)
         KOFF2 = KEND1
         KOFF3 = IT1AM(ISYMA,ISYMI) + 1
C
         NUMBMN = MAX(1,NMATIJ(ISYMMN))
         NUMBA  = MAX(1,NVIR(ISYMA))
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NMATIJ(ISYMMN),
     *              -ONE,WORK(KOFF1),NUMBMN,WORK(KOFF2),NUMBMN,
     *              ONE,OMEGA1(KOFF3),NUMBA)
C
      ENDDO
C
C--------------------------------------
C 3)
C     W^{ef}(dlnm) T^{ad}_{mn} g_{ielf}
C--------------------------------------
C
C-------------------------------------
C     Sort integrals
C     g(BiDl) sorted as g^{BD}_{li}
C--------------------------------------
C
      ISYMIL = MULD2H(ISYMBD,ISYINT)
      DO  ISYMI = 1,NSYM
         ISYML  = MULD2H(ISYMIL,ISYMI)
         ISYMBI = MULD2H(ISYMB,ISYMI)
         ISYMDL = MULD2H(ISYMD,ISYML)
         DO I = 1,NRHF(ISYMI)
            NBI = IT1AM(ISYMB,ISYMI) + NVIR(ISYMB)*(I-1) + B
            DO L = 1,NRHF(ISYML)
               NDL = IT1AM(ISYMD,ISYML) + NVIR(ISYMD)*(L-1) + D
               KOFF1 = IT2AM(ISYMBI,ISYMDL) + INDEX(NBI,NDL) 
               KOFF2 = KINTIL + IMATIJ(ISYML,ISYMI) 
     *                        + NRHF(ISYML)*(I-1) + L - 1
               WORK(KOFF2) = XIAJB(KOFF1)
            ENDDO
         ENDDO
      ENDDO
C
C----------------------
C     Construct TMAT
C----------------------
C
      DO I = 1, LENGTH
         TMAT(I) = WMAT(INDSQ(I,4))
      ENDDO
C
C-------------------------------------
C     Omega(ai) = Omega(ai) +
C     W^{BD}(elnm) T^{ea}_{nm} g_{iBlD}
C
C     Tmat(enml) * T_{enma} * Xiajb^{BD}_{li}
C--------------------------------------
C
      DO ISYMA = 1,NSYM
         ISYMI  = MULD2H(ISYRES,ISYMA)
         ISYENM = MULD2H(ISYMC2,ISYMA)
         ISYML  = MULD2H(ISELMN,ISYENM)
C
         KOFF1  = ISAIKJ(ISYENM,ISYML) + 1
         KOFF2  = KINTIL + IMATIJ(ISYML,ISYMI)
         KOFF3  = KEND1  + ISAIKJ(ISYENM,ISYMI)
         NUMENM = MAX(1,NCKI(ISYENM))
         NUML   = MAX(1,NRHF(ISYML))
         CALL DGEMM('N','N',NCKI(ISYENM),NRHF(ISYMI),NRHF(ISYML),
     *              ONE,TMAT(KOFF1),NUMENM,WORK(KOFF2),NUML,
     *              ZERO,WORK(KOFF3),NUMENM)
C
         KOFF1  = IT2SP(ISYENM,ISYMA)  + 1
         KOFF2  = KEND1 + ISAIKJ(ISYENM,ISYMI)
         KOFF3  = IT1AM(ISYMA,ISYMI)   + 1
         NUMENM = MAX(1,NCKI(ISYENM))
         NUMA   = MAX(1,NVIR(ISYMA))
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NCKI(ISYENM),
     *              -ONE,C2TP(KOFF1),NUMENM,WORK(KOFF2),NUMENM,
     *              ONE,OMEGA1(KOFF3),NUMA)
      ENDDO
C
C-------------------------------------
C 4)
C     W^{df}(elmn) T^{ad}_{mn} L_{ifle}
C 5)
C     W^{df}(emln) T^{ad}_{mn} L_{ielf}
C--------------------------------------
C
C---------------------------
C     Sort L integrals.
C---------------------------
C
      DO ISYML = 1, NSYM
         ISYMEI = MULD2H(ISYELI,ISYML)
         ISYMDL = MULD2H(ISYML,ISYMD)
         DO ISYME = 1, NSYM
            ISYMI  = MULD2H(ISYMEI,ISYME)
            ISYMEL = MULD2H(ISYME,ISYML)
C
            DO L = 1, NRHF(ISYML)
               NDL = IT1AM(ISYMD,ISYML) + NVIR(ISYMD)*(L-1) + D
               DO E = 1, NVIR(ISYME)
                  DO I = 1, NRHF(ISYMI)
                     NEI = IT1AM(ISYME,ISYMI) + NVIR(ISYME)*(I-1) + E
C
                     KOFF1 = IT2AM(ISYMDL,ISYMEI) + INDEX(NDL,NEI)
                     KOFF2 = KINT1 - 1
     *                     + ICKI(ISYMEL,ISYMI)
     *                     + NT1AM(ISYMEL)*(I-1)
     *                     + IT1AM(ISYME,ISYML)
     *                     + NVIR(ISYME)*(L-1)
     *                     + E
                     KOFF3 = KINT2 - 1
     *                     + ICKI(ISYMEI,ISYML)
     *                     + NT1AM(ISYMEI)*(L-1)
     *                     + IT1AM(ISYME,ISYMI)
     *                     + NVIR(ISYME)*(I-1)
     *                     + E
C
                     WORK(KOFF2) = YIAJB(KOFF1)
                     WORK(KOFF3) = YIAJB(KOFF1)
C
                  ENDDO
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C
C----------------------
C     Construct TMAT
C----------------------
C
      DO I = 1, LENGTH
         TMAT(I) =  WMAT(INDSQ(I,1))
C
      WORK(KTMAT-1+I) =  WMAT(I)
      ENDDO
C
C---------------------------------------------
C     Symmetry sorting if symmetry
C---------------------------------------------
C
      IF (NSYM .GT. 1) THEN
         CALL CC_GATHER(LENGTH,WORK(KEND1),TMAT,INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK(KEND1),1,TMAT,1)
C
         CALL CC_GATHER(LENGTH,WORK(KEND1),WORK(KTMAT),INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK(KEND1),1,WORK(KTMAT),1)
      ENDIF
C
C-------------------------------------
C     Contract
C-------------------------------------
C
      DO ISYMA = 1, NSYM
         ISYMI  = MULD2H(ISYRES,ISYMA)
         ISYMEL = MULD2H(ISYELI,ISYMI)
         ISYMMN = MULD2H(ISYMEL,ISELMN)
C
         KOFF1 = ISAIKL(ISYMEL,ISYMMN) + 1
         KOFF2 = KINT1
     *         + ICKI(ISYMEL,ISYMI)
         KOFF3 = KEND1
C
         NUMBEL = MAX(1,NT1AM(ISYMEL))
         NUMBMN = MAX(1,NMATIJ(ISYMMN))
C
         CALL DGEMM('T','N',NMATIJ(ISYMMN),NRHF(ISYMI),NT1AM(ISYMEL),
     *              ONE,TMAT(KOFF1),NUMBEL,WORK(KOFF2),NUMBEL,
     *              ZERO,WORK(KOFF3),NUMBMN)
C
         KOFF1 = KTMAT
     *         + ISAIKL(ISYMEL,ISYMMN)
         KOFF2 = KINT2
     *         + ICKI(ISYMEL,ISYMI)
         KOFF3 = KEND1
C
         NUMBEL = MAX(1,NT1AM(ISYMEL))
         NUMBMN = MAX(1,NMATIJ(ISYMMN))
C
         CALL DGEMM('T','N',NMATIJ(ISYMMN),NRHF(ISYMI),NT1AM(ISYMEL),
     *              ONE,WORK(KOFF1),NUMBEL,WORK(KOFF2),NUMBEL,
     *              ONE,WORK(KOFF3),NUMBMN)
C
         KOFF1 = KC2TEMP
     *         + IMAIJA(ISYMMN,ISYMA)
         KOFF2 = KEND1
         KOFF3 = IT1AM(ISYMA,ISYMI) + 1
C
         NUMBMN = MAX(1,NMATIJ(ISYMMN))
         NUMBA  = MAX(1,NVIR(ISYMA))
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NMATIJ(ISYMMN),
     *              ONE,WORK(KOFF1),NUMBMN,WORK(KOFF2),NUMBMN,
     *              ONE,OMEGA1(KOFF3),NUMBA)
C
      ENDDO
C--------------------------------------
C 6)
C     W^{ef}(dnlm) T^{ad}_{mn} L_{ielf}
C--------------------------------------
C
C-------------------------------------
C     Sort integrals
C     L(BiDl) sorted as L^{BD}_{li}
C--------------------------------------
C
      ISYMIL = MULD2H(ISYMBD,ISYINT)
      DO  ISYMI = 1,NSYM
         ISYML  = MULD2H(ISYMIL,ISYMI)
         ISYMBI = MULD2H(ISYMB,ISYMI)
         ISYMDL = MULD2H(ISYMD,ISYML)
         DO I = 1,NRHF(ISYMI)
            NBI = IT1AM(ISYMB,ISYMI) + NVIR(ISYMB)*(I-1) + B
            DO L = 1,NRHF(ISYML)
               NDL = IT1AM(ISYMD,ISYML) + NVIR(ISYMD)*(L-1) + D
               KOFF1 = IT2AM(ISYMBI,ISYMDL) + INDEX(NBI,NDL)  
               KOFF2 = KINTIL + IMATIJ(ISYML,ISYMI)
     *                        + NRHF(ISYML)*(I-1) + L - 1
               WORK(KOFF2) = YIAJB(KOFF1)
            ENDDO
         ENDDO
      ENDDO
C

C
C----------------------
C     Construct TMAT
C----------------------
C
      DO I = 1, LENGTH
         TMAT(I) = WMAT(INDSQ(I,3))
      ENDDO
C
C-------------------------------------
C     Omega(ai) = Omega(ai) +
C     W^{BD}(enlm) T^{ea}_{nm} L_{iBlD}
C
C     Tmat(enml) * T_{enma} * Yiajb^{BD}_{li}
C--------------------------------------
C
      DO ISYMA = 1,NSYM
         ISYMI  = MULD2H(ISYRES,ISYMA)
         ISYENM = MULD2H(ISYMC2,ISYMA)
         ISYML  = MULD2H(ISELMN,ISYENM)
C
         KOFF1  = ISAIKJ(ISYENM,ISYML) + 1
         KOFF2  = KINTIL + IMATIJ(ISYML,ISYMI)
         KOFF3  = KEND1  
         NUMENM = MAX(1,NCKI(ISYENM))
         NUML   = MAX(1,NRHF(ISYML))
C
         CALL DGEMM('N','N',NCKI(ISYENM),NRHF(ISYMI),NRHF(ISYML),
     *              ONE,TMAT(KOFF1),NUMENM,WORK(KOFF2),NUML,
     *              ZERO,WORK(KOFF3),NUMENM)
C
         KOFF1  = IT2SP(ISYENM,ISYMA)  + 1
         KOFF2  = KEND1 
         KOFF3  = IT1AM(ISYMA,ISYMI)   + 1
         NUMENM = MAX(1,NCKI(ISYENM))
         NUMA   = MAX(1,NVIR(ISYMA))
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NCKI(ISYENM),
     *              ONE,C2TP(KOFF1),NUMENM,WORK(KOFF2),NUMENM,
     *              ONE,OMEGA1(KOFF3),NUMA)
C
C
      ENDDO
C
C----------------------------
C     End.
C----------------------------
C
      CALL QEXIT('T3_ONEL2W')
C
      RETURN
      END
C  /* Deck t3_onel3w */
      SUBROUTINE T3_ONEL3W(OMEGA1,WMAT,TMAT,ISYMIM,XIAJB,XAIBJ,ISYINT,
     *                    C2TP,ISYMC2,INDSQ,LENSQ,WORK,LWORK,
     *                    ISYMB,B,ISYMD,D)
C
C     Written by Filip Pawlowski and Poul Jorgensen based on
C     SUBROUTINE T3_ONEL3 by K. Hald, Spring 2002.
C
C     Calculate the term (t^{def}_{lmn} - t^{def}_{lnm}) T^{de}_{lm} L_{ianf}
C     based on W intermediates
C
C     1)
C        ( W^{fe}(dlmn) - W^{fe}(dlnm) ) T^{de}_{lm} L_{ianf}
C
C     +  ( W^{fd}(emln) - W^{fd}(enlm) ) T^{de}_{lm} L_{ianf}
C
C      2)
C        ( W^{de}(fnml) - W^{de}(fmnl) ) T^{de}_{lm} L_{ianf}
C
C
C
C     Note : XIAJB is coming in as L and not g.
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYMIM, ISYINT, ISYMC2, LENSQ, LWORK, ISYMB, ISYMD
      INTEGER INDSQ(LENSQ,6), INDEX
      INTEGER ISYRE1, ISYRES, ISYMBD, ISFLMN, ISYAIN, ISYFLM, LENGTH
      INTEGER KTMAT, KC2TEMP, KINT, KEND1, LWRK1, ISYMM, ISYMFL
      INTEGER ISYMF, ISYML, ISYMLM, ISYMFM, KOFF1, KOFF2, KOFF3
      INTEGER ISYMI, ISYMAN, ISYMA, ISYMN, ISYMAI, ISYMBN, NBN
      INTEGER NAI, NUMFLM, NUMBAI
      INTEGER KTLM, ISYBLM, ISYMBL, ISYMFN
      INTEGER NUMFN, NUMAI
      INTEGER ISYRE2, IS1FLM
C
#if defined (SYS_CRAY)
      REAL OMEGA1(*), WMAT(*),  TMAT(*), XIAJB(*)
      REAL C2TP(*), WORK(LWORK), ZERO, ONE, XAIBJ(*)
      REAL DDOT
#else
      DOUBLE PRECISION OMEGA1(*), WMAT(*), TMAT(*), XIAJB(*)
      DOUBLE PRECISION C2TP(*), WORK(LWORK), ZERO, ONE, XAIBJ(*)
      DOUBLE PRECISION DDOT
#endif
C
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('T3_ONEL3W')
C
      ISYMBD = MULD2H(ISYMB,ISYMD)
C
      ISYRE1 = MULD2H(ISYMIM,ISYMC2)
      ISYRE2 = MULD2H(ISYRE1,ISYMBD)
      ISYRES = MULD2H(ISYRE2,ISYINT)
C
      ISFLMN = ISYMIM
      ISYAIN = MULD2H(ISYMB,ISYINT)
      ISYFLM = MULD2H(ISYMC2,ISYMD)
      ISYMLM = MULD2H(ISYMBD,ISYMC2)
C
      LENGTH = NCKIJ(ISFLMN)
C
      KTMAT   = 1
      KC2TEMP = KTMAT   + NCKIJ(ISFLMN)
      KINT    = KC2TEMP + NCKI(ISYFLM)
      KTLM    = KINT    + NCKI(ISYAIN) 
      KEND1    = KTLM    + NMATIJ(ISYMLM)
      LWRK1   = LWORK   - KEND1
C
C--------------------------------------------
C Calculate the terms:
C     1)
C        ( W^{fe}(dlmn) - W^{fe}(dlnm) ) T^{de}_{lm} L_{ianf}
C
C     +  ( W^{fd}(emln) - W^{fd}(enlm) ) T^{de}_{lm} L_{ianf}
C---------------------------------------------
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory in T3_ONEL3W (sort)')
      ENDIF
C
C-----------------------------
C     Sort C2
C-----------------------------
C
      DO ISYMM = 1, NSYM
         ISYMFL = MULD2H(ISYFLM,ISYMM)
         DO ISYMF = 1, NSYM
            ISYML = MULD2H(ISYMFL,ISYMF)
            ISYMLM = MULD2H(ISYMM,ISYML)
            ISYMFM = MULD2H(ISYMF,ISYMM)
C
            DO M = 1, NRHF(ISYMM)
               DO L = 1, NRHF(ISYML)
C
                  KOFF1 = IT2SP(ISYFLM,ISYMD)
     *                  + NCKI(ISYFLM)*(D-1)
     *                  + ICKI(ISYMFM,ISYML)
     *                  + NT1AM(ISYMFM)*(L-1)
     *                  + IT1AM(ISYMF,ISYMM)
     *                  + NVIR(ISYMF)*(M-1) 
     *                  + 1
C
                  KOFF2 = KC2TEMP - 1
     *                  + ICKI(ISYMFL,ISYMM)
     *                  + NT1AM(ISYMFL)*(M-1)
     *                  + IT1AM(ISYMF,ISYML)
     *                  + NVIR(ISYMF)*(L-1) 
     *                  + 1

C
                  CALL DCOPY(NVIR(ISYMF),C2TP(KOFF1),1,WORK(KOFF2),1)
C
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C---------------------------
C     Sort integrals.
C---------------------------
C
      DO ISYMI = 1, NSYM
         ISYMAN = MULD2H(ISYAIN,ISYMI)
         DO ISYMA = 1, NSYM
            ISYMN  = MULD2H(ISYMAN,ISYMA)
            ISYMAI = MULD2H(ISYMA,ISYMI)
            ISYMBN = MULD2H(ISYMB,ISYMN)
C
            DO N = 1, NRHF(ISYMN)
               NBN = IT1AM(ISYMB,ISYMN) + NVIR(ISYMB)*(N-1) + B
               DO A = 1, NVIR(ISYMA)
                  DO I = 1, NRHF(ISYMI)
                     NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
C
                     KOFF1 = IT2AM(ISYMBN,ISYMAI) + INDEX(NBN,NAI)
                     KOFF2 = KINT - 1
     *                     + ICKI(ISYMAI,ISYMN)
     *                     + NT1AM(ISYMAI)*(N-1)
     *                     + IT1AM(ISYMA,ISYMI)
     *                     + NVIR(ISYMA)*(I-1)
     *                     + A
C
                     WORK(KOFF2) = XIAJB(KOFF1)
C
                  ENDDO
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C
C----------------------
C     Construct TMAT
C----------------------
C
      DO I = 1, LENGTH
         TMAT(I) = WMAT(I)
     *           - WMAT(INDSQ(I,3))
C
      WORK(KTMAT-1+I) = WMAT(INDSQ(I,1))
     *                - WMAT(INDSQ(I,4))
      ENDDO
C
C
C-------------------------------------
C     Contract
C-------------------------------------
C
      ISYMAI = ISYRES
C
      ISYMN = MULD2H(ISYRES,ISYAIN)
C
      IS1FLM = MULD2H(ISYMIM,ISYMN)
      KOFF1 = ISAIKJ(IS1FLM,ISYMN) + 1
      KOFF2 = IT2SP(ISYFLM,ISYMD)
     *      + NCKI(ISYFLM)*(D-1)
     *      + 1
      KOFF3 = KEND1
C
      CALL DZERO(WORK(KOFF3),NRHF(ISYMN))
C
      NUMFLM = MAX(1,NCKI(ISYFLM))
C
C
      CALL DGEMV('T',NCKI(ISYFLM),NRHF(ISYMN),ONE,
     *           TMAT(KOFF1),NUMFLM,C2TP(KOFF2),1,
     *           ZERO,WORK(KOFF3),1)
C
C
      KOFF1 = KTMAT
     *      + ISAIKJ(IS1FLM,ISYMN)
      KOFF2 = KC2TEMP
      KOFF3 = KEND1
C
      NUMFLM = MAX(1,NCKI(ISYFLM))
C
C
      CALL DGEMV('T',NCKI(ISYFLM),NRHF(ISYMN),ONE,
     *           WORK(KOFF1),NUMFLM,WORK(KOFF2),1,
     *           ONE,WORK(KOFF3),1)
C
C
      KOFF1 = KINT
     *      + ICKI(ISYMAI,ISYMN)
      KOFF2 = KEND1
      KOFF3 = 1
C
      NUMBAI = MAX(1,NT1AM(ISYMAI))
C
C
      CALL DGEMV('N',NT1AM(ISYMAI),NRHF(ISYMN),-ONE,
     *           WORK(KOFF1),NUMBAI,WORK(KOFF2),1,
     *           ONE,OMEGA1(KOFF3),1)
C
C------------------------------------------------------------
C Calculate the terms:
C      2)
C        ( W^{de}(fnml) - W^{de}(fmnl) ) T^{de}_{lm} L_{ianf}
C------------------------------------------------------------
C
C-------------------------------------
C     Sort multipliers
C     T^{BD}_{lm}  sorted as T_{lm} 
C--------------------------------------
C
      ISYMLM = MULD2H(ISYMBD,ISYMC2)
      DO  ISYMM = 1,NSYM
         ISYML  = MULD2H(ISYMLM,ISYMM)
         ISYBLM = MULD2H(ISYMB,ISYMLM)
         ISYMBL = MULD2H(ISYMB,ISYML)
         DO L = 1,NRHF(ISYML)
            DO M = 1,NRHF(ISYMM)
C
               KOFF1 = IT2SP(ISYBLM,ISYMD)
     *               + NCKI(ISYBLM)*(D-1)
     *               + ICKI(ISYMBL,ISYMM)
     *               + NT1AM(ISYMBL)*(M-1)
     *               + IT1AM(ISYMB,ISYML)
     *               + NVIR(ISYMB)*(L-1)
     *               + B
C
               KOFF2 = KTLM + IMATIJ(ISYML,ISYMM)
     *                        + NRHF(ISYML)*(M-1) + L - 1
C

               WORK(KOFF2) = C2TP(KOFF1)
            ENDDO
         ENDDO
      ENDDO
C     write(lupri,*)'Norm T_{lm} = ',ddot(NMATIJ(ISYMLM),
C    *  WORK(KTLM),1,WORK(KTLM),1)
C
C----------------------
C     Construct TMAT
C----------------------
C
      DO I = 1, LENGTH
         TMAT(I) = WMAT(INDSQ(I,3)) - WMAT(INDSQ(I,4))
      ENDDO
C     write(lupri,*)'Norm TMAT = ',ddot(LENGTH,
C    *  TMAT,1,TMAT,1)
C
C---------------------------------------------
C     Symmetry sorting if symmetry
C---------------------------------------------
C
      IF (NSYM .GT. 1) THEN
         CALL CC_GATHER(LENGTH,WORK(KEND1),TMAT,INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK(KEND1),1,TMAT,1)
      ENDIF
C
C------------------------------------------------------
C     Omega(ai) = Omega(ai) +
C   ( W^{BD}(fnml) - W^{BD}(fmnl) ) T^{BD}_{lm} L_{ianf}
C
C     Tmat(fnlm) * T_{lm} * Xaibj^{aifn}
C------------------------------------------------------
C
      ISYMFN = MULD2H(ISFLMN,ISYMLM)
      ISYMAI = MULD2H(ISYMFN,ISYINT)
C
      KOFF1 = ISAIKL(ISYMFN,ISYMLM) + 1
      KOFF2  = KTLM 
      KOFF3  = KEND1 
      NUMFN  = MAX(1,NT1AM(ISYMFN))
C
C
      CALL DZERO(WORK(KOFF3),NT1AM(ISYMFN))
C
      CALL DGEMV('N',NT1AM(ISYMFN),NMATIJ(ISYMLM),ONE,
     *           TMAT(KOFF1),NUMFN,WORK(KOFF2),1,
     *           ZERO,WORK(KOFF3),1)
C
      KOFF1  = IT2SQ(ISYMAI,ISYMFN) + 1
      KOFF2  = KEND1 
      KOFF3  =  1
      NUMAI  = MAX(1,NT1AM(ISYMAI))
C
      CALL DGEMV('N',NT1AM(ISYMAI),NT1AM(ISYMFN),-ONE,
     *           XAIBJ(KOFF1),NUMAI,WORK(KOFF2),1,
     *           ONE,OMEGA1(KOFF3),1)
C
C----------------------------
C     End.
C----------------------------
C
      CALL QEXIT('T3_ONEL3W')
      RETURN
      END
