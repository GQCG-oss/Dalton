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
C Deck /* cc3_abmatt3zu */
      SUBROUTINE CC3_ABMATT3ZU(LISTX,IDLSTX,LISTZU,IDLSTZU,
     &                                OMEGA2,
     &                                OMEGA2EFF,
     *                                ISYRES,WORK,LWORK)
C
      IMPLICIT NONE
#include <priunit.h>
C
      CHARACTER*3 LISTX,LISTZU
      INTEGER IDLSTX,IDLSTZU
C
      INTEGER ISYRES, LWORK
C
#if defined (SYS_CRAY)
      REAL OMEGA2(*),OMEGA2EFF(*)
      REAL WORK(LWORK)
#else
      DOUBLE PRECISION OMEGA2(*),OMEGA2EFF(*)
      DOUBLE PRECISION WORK(LWORK)
#endif
C
      CALL QENTER('ABT3ZU')
C
      CALL CC3_AMATT3ZU(LISTX,IDLSTX,LISTZU,IDLSTZU,
     &                 OMEGA2,OMEGA2EFF,
     *                 ISYRES,WORK,LWORK)
C
      CALL CC3_BMATT3ZU(LISTX,IDLSTX,LISTZU,IDLSTZU,
     &                 OMEGA2,OMEGA2EFF,
     *                 ISYRES,WORK,LWORK)
C
C----------
C     End .
C----------
C
      CALL QEXIT('ABT3ZU')
C
      RETURN
      END
C /* Deck cc3_amatt3zu */
      SUBROUTINE CC3_AMATT3ZU(LISTX,IDLSTX,LISTZU,IDLSTZU,
     &                                OMEGA2,
     &                                OMEGA2EFF,
     *                                ISYRES,WORK,LWORK)
C
      IMPLICIT NONE
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccsdinp.h>
C
      INTEGER LUCKJD, LUDELD, LU3VI, LUTOC,LU3VI2, LUDKBC, LU3FOP
C
      CHARACTER*6 FNCKJD, FNDELD, FN3VI
      CHARACTER*8 FNTOC,FN3VI2 
      CHARACTER*4 FNDKBC
      CHARACTER*5 FN3FOP
C        
      PARAMETER (FNCKJD='CKJDEL', FNTOC   = 'CCSDT_OC' )
      PARAMETER (FNDELD  = 'CKDELD'  , FNDKBC  = 'DKBC' )
      PARAMETER (FN3VI   = 'CC3_VI'  , FN3VI2  = 'CC3_VI12')
      PARAMETER (FN3FOP='PTFOP' )
C
      CHARACTER*3 LISTX,LISTZU,
      INTEGER IDLSTX,IDLSTZU
C
      INTEGER ISYRES, LWORK
      integer isym
C
#if defined (SYS_CRAY)
      REAL OMEGA2(*),OMEGA2EFF(*)
      REAL WORK(LWORK)
#else
      DOUBLE PRECISION OMEGA2(*),OMEGA2EFF(*)
      DOUBLE PRECISION WORK(LWORK)
#endif
C
      CALL QENTER('AT3ZU')
C
C--------------------------------
C     Open files
C--------------------------------
C
      LUCKJD   = -1
      LUTOC    = -1
      LUDELD   = -1
      LUDKBC   = -1
      LU3VI2   = -1
      LU3VI    = -1
      LU3FOP   = -1
C
      CALL WOPEN2(LUCKJD,FNCKJD,64,0)
      CALL WOPEN2(LUTOC,FNTOC,64,0)
      CALL WOPEN2(LUDELD,FNDELD,64,0)
      CALL WOPEN2(LUDKBC,FNDKBC,64,0)
      CALL WOPEN2(LU3VI2,FN3VI2,64,0)
      CALL WOPEN2(LU3VI,FN3VI,64,0)
      CALL WOPEN2(LU3FOP,FN3FOP,64,0)
C
      CALL CC3_AMATT3ZUVIR(LISTX,IDLSTX,LISTZU,IDLSTZU,
     &                                OMEGA2,
     &                                OMEGA2EFF,
     *                                ISYRES,WORK,LWORK,
     *                                LUCKJD,FNCKJD,LUTOC,FNTOC,
     *                                LUDELD,FNDELD,LUDKBC,FNDKBC,
     *                                LU3VI2,FN3VI2,LU3VI,FN3VI)
C
c
      CALL CC3_AMATT3ZUOCC(LISTX,IDLSTX,LISTZU,IDLSTZU,
     &                                OMEGA2,
     &                                OMEGA2EFF,
     *                                ISYRES,WORK,LWORK,
     *                                LUCKJD,FNCKJD,LUTOC,FNTOC,
     *                                LUDELD,FNDELD,LU3FOP,FN3FOP,
     *                                LU3VI,FN3VI)
c
      CALL WCLOSE2(LUCKJD,FNCKJD,'KEEP')
      CALL WCLOSE2(LUTOC,FNTOC,'KEEP')
      CALL WCLOSE2(LUDELD,FNDELD,'KEEP')
      CALL WCLOSE2(LUDKBC,FNDKBC,'KEEP')
      CALL WCLOSE2(LU3VI2,FN3VI2,'KEEP')
      CALL WCLOSE2(LU3VI,FN3VI,'KEEP')
      CALL WCLOSE2(LU3FOP,FN3FOP,'KEEP')
C
C----------
C     End .
C----------
C
      CALL QEXIT('AT3ZU')
C
      RETURN
      END
C    /* Deck cc3_amatt3zuvir */ 
      SUBROUTINE CC3_AMATT3ZUVIR(LISTX,IDLSTX,LISTZU,IDLSTZU,
     &                                OMEGA2,
     &                                OMEGA2EFF,
     *                                ISYRES,WORK,LWORK,
     *                                LUCKJD,FNCKJD,LUTOC,FNTOC,
     *                                LUDELD,FNDELD,LUDKBC,FNDKBC,
     *                                LU3VI2,FN3VI2,LU3VI,FN3VI)
C
*
*******************************************************************
*
* This routine calculates the contribution to cubic response 
* (originating from both F matrix and B matrix):
*
* omega2eff = <mu2|[H^X,T3^ZU]|HF>
*
* which is contracted outside with the multipliers (zeroth-order for 
* F matrix or first-order for B matrix).
*
*
*
* T3^ZU is calculated using W^BD intermmediate.
*
* !!! ACTUALLY HERE WE CALCULATE ONLY THE PART OF THE TERM, WHICH ORIGINATES
*     FROM THE AY-MATRIX CONTRIBUTIONS TO T3^ZU 
*     (to get the rest of contribution call cc3_bmatt3zu routine) !!!
*
* Thus here:
*
* W^BD =  <mu3|[[Z,T1^U],T3^0]|HF> + <mu3|[[Z,T2^U],T2^0]|HF>
*      +  <mu3|[Z,T3^U]|HF>
*
* (!!! TO GET THE TOTAL CONTRIBUTION FROM THE LAST TERM CALL CC3_AMATT3ZUOCC !!!)
*
* F. Pawlowski, P. Jorgensen, Aarhus Spring 2003.
*
**************************************************************

C
      IMPLICIT NONE
#include <priunit.h>
#include <ccsdsym.h>
#include <dummy.h>
#include <iratdef.h>
#include <inftap.h>
#include <ccinftap.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccr1rsp.h>
#include <ccexci.h>
C
      CHARACTER*(*) FNCKJD,FNTOC,FNDELD,FNDKBC,FN3VI2,FN3VI
      INTEGER LUCKJD,LUTOC,LUDELD,LUDKBC,LU3VI2,LU3VI
C
      CHARACTER*12 FN3SRTR, FNCKJDRZ, FNDELDRZ, FNDKBCRZ
      CHARACTER*13 FNCKJDRU,FNDELDRU,FNDKBCRU
      INTEGER LU3SRTR, LUCKJDRZ, LUDELDRZ, LUDKBCRZ
      INTEGER LUCKJDRU,LUDELDRU,LUDKBCRU
C
      CHARACTER CDUMMY*1
C
      PARAMETER(FN3SRTR  = 'CCSDT_FBMAT1', FNCKJDRZ  = 'CCSDT_FBMAT2', 
     *          FNDELDRZ  = 'CCSDT_FBMAT3', FNDKBCRZ  = 'CCSDT_FBMAT4')
C
      PARAMETER(FNCKJDRU  = 'CCSDT_FBMAT22',FNDELDRU  = 'CCSDT_FBMAT33',
     *          FNDKBCRU  = 'CCSDT_FBMAT44')
C
      CHARACTER*3 LISTRZ, LISTRU, LISTX, LISTZU
      CHARACTER*8 LABELRZ,LABELRU
      INTEGER     IDLSTRZ,IDLSTRU,IDLSTX,IDLSTZU
      INTEGER     ISYMRZ, ISYMRU
C
      LOGICAL T2XNET2Z
      LOGICAL LORXRZ,LORXRU
C
      INTEGER ISYRES, LWORK
      INTEGER ISYM0,ISINT1,ISINT2,ISYMT2,ISYMT3
      INTEGER KT2TP,KFOCKD,KLAMDP,KLAMDH,KEND00,LWRK00,KXIAJB
      INTEGER IOPTTCME,IOPT,IOFF
      INTEGER ISYMWZ,ISYMWZU,ISYFCKZ,ISYFCKU,ISYMZU
      INTEGER KT2Z,KT2U,KFCKZ,KFCKU,KEND0,LWRK0
      INTEGER KFCKZUV,KFCKZUO,KFCKUZV,KFCKUZO,KT1Z,KT1U,KEND1,LWRK1
      INTEGER KLAMDPU,KLAMDHU,KLAMDPZ,KLAMDHZ
      INTEGER ISINTRZ1,ISINTRZ2,ISINTRU1,ISINTRU2
      INTEGER KRBJIAZU,KT3OG1,KT3OG2,KW3ZOGZ1,KW3UOGU1,KTROC,KTROC1
      INTEGER KINTOC
      INTEGER ISYMD,ISCKB1,ISCKB2,ISCKB2Z,ISCKB2U
      INTEGER KT3VDG3,KT3VDG2,KEND2,LWRK2
      INTEGER KT3VDG1,KEND3,LWRK3,KW3ZVDGZ1,KW3UVDGU1,KT3VDGZ3,KT3VDGU3
      INTEGER KTRVI,KTRVI1,MAXZU
      INTEGER ISYMB,ISYALJ,ISYALJ2,ISYMBD,ISCKIJ,ISCKD2,ISCKD2Z,ISCKD2U
      INTEGER ISWMATZ,ISWMATU,ISWMATZU
      INTEGER KSMAT,KSMAT3,KUMAT,KUMAT3,KDIAG,KT3MAT
      INTEGER KDIAGWZ,KWMATZ,KINDSQWZ,KDIAGWU,KWMATU,KINDSQWU
      INTEGER KINDSQ,KINDEX,KINDEX2
      INTEGER MAXZZUU,KTMATW,KT3VBG1,KT3VBG2,KT3VBG3,KEND4,LWRK4
      INTEGER KW3ZVDGZ2,KW3UVDGU2,KT3VBGZ3,KT3VBGU3,KINTVI
      INTEGER KINDSQWZU,KDIAGWZU,KWMATZU
      INTEGER LENSQ,LENSQWZ,LENSQWU,LENSQWZU
      INTEGER AIBJCK_PERM
c
      INTEGER ISYMRX,ISYFCKX
      INTEGER KLAMDPX,KLAMDHX,KFOCKX,KT1X,ISYMC,ISYMK,KOFF1,KOFF2
      INTEGER ISINT1X,ISCKB1X
c
      integer kx3am
c
#if defined (SYS_CRAY)
      REAL OMEGA2(*),OMEGA2EFF(*)
      REAL WORK(LWORK)
      REAL FREQRZ,FREQRU,FREQZU
      REAL XNORMVAL,DDOT
#else
      DOUBLE PRECISION OMEGA2(*),OMEGA2EFF(*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION FREQRZ,FREQRU,FREQZU
      DOUBLE PRECISION XNORMVAL,DDOT
#endif
C
      CALL QENTER('AT3ZUV')
c
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
C     ISYRES  is symmetry of result vector 
C
C-------------------------------------------------------------
C
      IPRCC   = IPRINT
      ISYM0   = 1
      ISINT1  = 1
      ISINT2  = 1
      ISYMT2  = 1
      ISYMT3  = MULD2H(ISINT2,ISYMT2)
C
C--------------------------------
C     Open temporary files
C--------------------------------
C
      LU3SRTR   = -1
      LUCKJDRZ  = -1
      LUCKJDRU  = -1
      LUDELDRZ  = -1
      LUDELDRU  = -1
      LUDKBCRZ  = -1
      LUDKBCRU  = -1
C
      CALL WOPEN2(LU3SRTR,FN3SRTR,64,0)
      CALL WOPEN2(LUCKJDRZ,FNCKJDRZ,64,0)
      CALL WOPEN2(LUCKJDRU,FNCKJDRU,64,0)
      CALL WOPEN2(LUDELDRZ,FNDELDRZ,64,0)
      CALL WOPEN2(LUDELDRU,FNDELDRU,64,0)
      CALL WOPEN2(LUDKBCRZ,FNDKBCRZ,64,0)
      CALL WOPEN2(LUDKBCRU,FNDKBCRU,64,0)
C
C----------------------------------------------
C     Calculate the zeroth order stuff once
C----------------------------------------------
C
      KT2TP  = 1 
      KFOCKD  = KT2TP  + NT2SQ(ISYMT2)
      KLAMDP = KFOCKD + NORBTS
      KLAMDH = KLAMDP + NLAMDT
      KEND00 = KLAMDH + NLAMDT
      LWRK00 = LWORK  - KEND00
C
      KXIAJB = KEND00
      KEND00 = KXIAJB + NT2AM(ISINT1)
      LWRK00 = LWORK  - KEND00
C
      IF (LWRK00 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_AMATT3ZUVIR (zeroth allo.')
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
            call quit('Insufficient space (T3) in CC3_AMATT3ZUVIR')
         END IF
      endif
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
C
C determining R1 labels
C
      ISYMRX = ISYLRT(IDLSTX)
      !FREQuency not needed for LISTX
      !LABELX = LRTLBL(IDLSTX) not needed for LISTX
C
      IF (LISTZU(1:3).EQ.'R2 ') THEN
        LISTRZ  = 'R1 '
        LABELRZ = LBLR2T(IDLSTZU,1)
        ISYMRZ  = ISYR2T(IDLSTZU,1)
        FREQRZ  = FRQR2T(IDLSTZU,1)
        LORXRZ  = LORXR2T(IDLSTZU,1)
        IDLSTRZ = IR1TAMP(LABELRZ,LORXRZ,FREQRZ,ISYMRZ)
C
        IF (LORXRZ) CALL QUIT('NO ORBITAL RELAX. IN CC3_AMATT3ZUVIR')

        LISTRU  = 'R1 '
        LABELRU = LBLR2T(IDLSTZU,2)
        ISYMRU  = ISYR2T(IDLSTZU,2)
        FREQRU  = FRQR2T(IDLSTZU,2)
        LORXRU  = LORXR2T(IDLSTZU,2)
        IDLSTRU = IR1TAMP(LABELRU,LORXRU,FREQRU,ISYMRU)
C
        IF (LORXRU) CALL QUIT('NO ORBITAL RELAX. IN CC3_AMATT3ZUVIR')
C
      ELSE 
         CALL QUIT('Unknown list for LISTZU in CC3_AMATT3ZUVIR.')
      END IF
C
      FREQZU = FREQRZ + FREQRU
C
      ISYMWZ   = MULD2H(ISYMT3,ISYMRZ)
C
      ISYMWZU   = MULD2H(ISYMWZ,ISYMRU)
C
      IF (ISYRES.NE.ISYMWZU) THEN
        WRITE(LUPRI,*) 'INCONSISTENT SYMMETRY SPECIFICATION'
        WRITE(LUPRI,*) 
     *      'ISYRES =',ISYRES,'ISYMWZU =',ISYMWZU
        CALL QUIT('CC3_AMATT3ZUVIR inconsistent symmetry specification')
      END IF
C
      ISYFCKZ = ISYMRZ
      ISYFCKU = ISYMRU
      ISYMZU = MULD2H(ISYMRZ,ISYMRU)
C
C-------------------------------------------------
C     Read T1^Z and T2^Z
C     Calculate (ck|de)-tilde(Z) and (ck|lm)-tilde(Z)
C     Used to construct T3^Z
C
C     Read T1^U and T2^U
C     Calculate (ck|de)-tilde(U) and (ck|lm)-tilde(U)
C     Used to construct T3^U
C-------------------------------------------------
C
      ISYFCKX = MULD2H(ISYMOP,ISYMRX)
C
      KT2Z   = KEND00
      KT2U   = KT2Z   + NT2SQ(ISYMRZ)
      KFCKZ  = KT2U   + NT2SQ(ISYMRU)
      KFCKU  = KFCKZ  + N2BST(ISYFCKZ)
      KEND0  = KFCKU  + N2BST(ISYFCKU)
      LWRK0  = LWORK  - KEND0
C
      KFCKZUV = KEND0
      KFCKZUO = KFCKZUV + N2BST(ISYMZU)
      KFCKUZV = KFCKZUO + N2BST(ISYMZU)
      KFCKUZO = KFCKUZV + N2BST(ISYMZU)
      KEND0   = KFCKUZO + N2BST(ISYMZU)
C
      KLAMDPX = KEND0
      KLAMDHX = KLAMDPX + NLAMDT
      KEND0   = KLAMDHX + NLAMDT
      LWRK0   = LWORK   - KEND0
C
      KFOCKX  = KEND0
      KEND0   = KFOCKX  + N2BST(ISYFCKX)
      LWRK0   = LWORK  - KEND0
C
      KT1Z   = KEND0
      KT1U   = KT1Z   + NT1AM(ISYMRZ)
      KEND1  = KT1U   + NT1AM(ISYMRU)
C
      KT1X   = KEND1
      KEND1  = KT1X   + NT1AM(ISYMRX)
      LWRK1  = LWORK  - KEND1
C
      KLAMDPU = KEND1
      KLAMDHU = KLAMDPU + NLAMDT
      KLAMDPZ = KLAMDHU + NLAMDT
      KLAMDHZ = KLAMDPZ + NLAMDT
      KEND1   = KLAMDHZ + NLAMDT
      LWRK1   = LWORK   - KEND1
C
      IF (LWRK1 .LT. NT2SQ(ISYMRZ)) THEN
         CALL QUIT('Out of memory in CC3_AMATT3ZUVIR (TOT_T3Y) ')
      ENDIF
C        
C---------------------------
C     Read Lambda_X matrices
C---------------------------
C     
      CALL GET_LAMBDAX(WORK(KLAMDPX),WORK(KLAMDHX),LISTX,IDLSTX,
     *                 ISYMRX,WORK(KLAMDP),WORK(KLAMDH),WORK(KEND1),
     *                 LWRK1)
C
C     T1^Z and T2^Z amplitudes
C     
      IOPT = 3
      CALL GET_T1_T2(IOPT,.TRUE.,WORK(KT1Z),WORK(KT2Z),LISTRZ,IDLSTRZ,
     *               ISYMRZ,WORK(KEND1),LWRK1)
C
      IF (IPRINT .GT. 55) THEN
         XNORMVAL = DDOT(NT2SQ(ISYMRZ),WORK(KT2Z),1,WORK(KT2Z),1)
         WRITE(LUPRI,*) 'Norm of T2X  ',XNORMVAL
         XNORMVAL = DDOT(NT1AM(ISYMRZ),WORK(KT1Z),1,WORK(KT1Z),1)
         WRITE(LUPRI,*) 'Norm of T1X  ',XNORMVAL
      ENDIF
C
C     T1^U and T2^U amplitudes
C     
      IOPT = 3
      CALL GET_T1_T2(IOPT,.TRUE.,WORK(KT1U),WORK(KT2U),LISTRU,IDLSTRU,
     *               ISYMRU,WORK(KEND1),LWRK1)
C
      IF (IPRINT .GT. 55) THEN
         XNORMVAL = DDOT(NT2SQ(ISYMRU),WORK(KT2U),1,WORK(KT2U),1)
         WRITE(LUPRI,*) 'Norm of T2Y  ',XNORMVAL
         XNORMVAL = DDOT(NT1AM(ISYMRU),WORK(KT1U),1,WORK(KT1U),1)
         WRITE(LUPRI,*) 'Norm of T1Y  ',XNORMVAL
      ENDIF
C
C------------------
C     Read in T1^X 
C------------------
C
      IOPT = 1
      CALL GET_T1_T2(IOPT,.FALSE.,WORK(KT1X),DUMMY,LISTX,IDLSTX,
     *               ISYMRX,WORK(KEND1),LWRK1)
C
      IF (IPRINT .GT. 55) THEN
         XNORMVAL = DDOT(NT1AM(ISYMRX),WORK(KT1X),1,WORK(KT1X),1)
         WRITE(LUPRI,*) 'Norm of T1A  ',XNORMVAL
      ENDIF
C
C------------------------------------------
C        Calculate the F^X matrix
C------------------------------------------
C
      CALL CC3LR_MFOCK(WORK(KEND1),WORK(KT1X),WORK(KXIAJB),ISYFCKX)
C
C     Put the F_{kc} part into correct F_{pq}
C
         CALL DZERO(WORK(KFOCKX),N2BST(ISYFCKX))
C
         DO ISYMC = 1, NSYM
            ISYMK = MULD2H(ISYFCKX,ISYMC)
            DO K = 1, NRHF(ISYMK)
               DO C = 1, NVIR(ISYMC)
                  KOFF1 = KFOCKX - 1
     *                  + IFCVIR(ISYMK,ISYMC)
     *                  + NORB(ISYMK)*(C - 1)
     *                  + K
                  KOFF2 = KEND1 - 1
     *                  + IT1AM(ISYMC,ISYMK)
     *                  + NVIR(ISYMC)*(K - 1)
     *                  + C
C
                  WORK(KOFF1) = WORK(KOFF2)
C
               ENDDO
            ENDDO
         ENDDO
C
         IF (IPRINT .GT. 55) THEN
            CALL AROUND( 'In CC3_BMATT3ZU: Fock^A MO matrix' )
            CALL CC_PRFCKMO(WORK(KFOCKX),ISYFCKX)
         ENDIF
C

C
C     Integrals H^Z in T3^Z
C

      ISINTRZ1 = MULD2H(ISINT1,ISYMRZ)
      ISINTRZ2 = MULD2H(ISINT2,ISYMRZ)
C
      CALL CC3_BARINT(WORK(KT1Z),ISYMRZ,WORK(KLAMDP),
     *                WORK(KLAMDH),WORK(KEND1),LWRK1,
     *                LU3SRTR,FN3SRTR,LUCKJDRZ,FNCKJDRZ)
C
      CALL CC3_SORT1(WORK(KEND1),LWRK1,2,ISINTRZ1,LU3SRTR,FN3SRTR,
     *               LUDELDRZ,FNDELDRZ,IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *               IDUMMY,CDUMMY)
C
      CALL CC3_SINT(WORK(KLAMDH),WORK(KEND1),LWRK1,ISINTRZ1,
     *              LUDELDRZ,FNDELDRZ,LUDKBCRZ,FNDKBCRZ)
C
C     Integrals H^U in T3^U
C

      ISINTRU1 = MULD2H(ISINT1,ISYMRU)
      ISINTRU2 = MULD2H(ISINT2,ISYMRU)
C
      CALL CC3_BARINT(WORK(KT1U),ISYMRU,WORK(KLAMDP),
     *                WORK(KLAMDH),WORK(KEND1),LWRK1,
     *                LU3SRTR,FN3SRTR,LUCKJDRU,FNCKJDRU)
C
      CALL CC3_SORT1(WORK(KEND1),LWRK1,2,ISINTRU1,LU3SRTR,FN3SRTR,
     *               LUDELDRU,FNDELDRU,IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *               IDUMMY,CDUMMY)
C
      CALL CC3_SINT(WORK(KLAMDH),WORK(KEND1),LWRK1,ISINTRU1,
     *              LUDELDRU,FNDELDRU,LUDKBCRU,FNDKBCRU)
C
C------------------------------------------
C     Calculate the F^Z matrix
C------------------------------------------
C
      CALL GET_FOCKX(WORK(KFCKZ),LABELRZ,IDLSTRZ,ISYMRZ,WORK(KLAMDP),
     *                  ISYM0,WORK(KLAMDH),ISYM0,WORK(KEND1),LWRK1)
C
C------------------------------------------
C     Calculate the F^U matrix
C------------------------------------------
C
      CALL GET_FOCKX(WORK(KFCKU),LABELRU,IDLSTRU,ISYMRU,WORK(KLAMDP),
     *                  ISYM0,WORK(KLAMDH),ISYM0,WORK(KEND1),LWRK1)
C
C------------------------------------------
C     Calculate the [Z,T1^U] matrix
C     Recall that we only need the occ-occ and vir-vir block.
C------------------------------------------
C
      CALL GET_LAMBDAX(WORK(KLAMDPU),WORK(KLAMDHU),LISTRU,IDLSTRU,
     *                 ISYMRU,WORK(KLAMDP),WORK(KLAMDH),WORK(KEND1),
     *                 LWRK1)
      ! get vir-vir block Z_(c-,d)
      CALL GET_FOCKX(WORK(KFCKZUV),LABELRZ,IDLSTRZ,ISYMRZ,WORK(KLAMDPU),
     *                  ISYMRU,WORK(KLAMDH),ISYM0,WORK(KEND1),LWRK1)
      ! get occ-occ block Z_(l,k-)
      CALL GET_FOCKX(WORK(KFCKZUO),LABELRZ,IDLSTRZ,ISYMRZ,WORK(KLAMDP),
     *                  ISYM0,WORK(KLAMDHU),ISYMRU,WORK(KEND1),LWRK1)
C
C------------------------------------------
C     Calculate the [U,T1^Z] matrix
C     Recall that we only need the occ-occ and vir-vir block.
C------------------------------------------
C
      CALL GET_LAMBDAX(WORK(KLAMDPZ),WORK(KLAMDHZ),LISTRZ,IDLSTRZ,
     *                 ISYMRZ,WORK(KLAMDP),WORK(KLAMDH),WORK(KEND1),
     *                 LWRK1)
      ! get vir-vir block U_(c-,d)
      CALL GET_FOCKX(WORK(KFCKUZV),LABELRU,IDLSTRU,ISYMRU,WORK(KLAMDPZ),
     *                  ISYMRZ,WORK(KLAMDH),ISYM0,WORK(KEND1),LWRK1)
      ! get occ-occ block U_(l,k-)
      CALL GET_FOCKX(WORK(KFCKUZO),LABELRU,IDLSTRU,ISYMRU,WORK(KLAMDP),
     *                  ISYM0,WORK(KLAMDHZ),ISYMRZ,WORK(KEND1),LWRK1)
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
      ISINT1X  = MULD2H(ISINT1,ISYMRX)
C
      KRBJIAZU = KEND0
      KT3OG1   = KRBJIAZU + NT2SQ(ISYRES)
      KT3OG2   = KT3OG1 + NTRAOC(ISINT2)
      KEND1    = KT3OG2+ NTRAOC(ISINT2)
      LWRK1    = LWORK  - KEND1
C
      KW3ZOGZ1 = KEND1
      KEND1   = KW3ZOGZ1 + NTRAOC(ISINTRZ2)
C
      KW3UOGU1 = KEND1
      KEND1   = KW3UOGU1 + NTRAOC(ISINTRU2)
      LWRK1  = LWORK  - KEND1
C
      KTROC  = KEND1
      KTROC1 = KTROC  + NTRAOC(ISINT1X)   !KTROC - int. in contraction
      KEND1  = KTROC1 + NTRAOC(ISINT1X)   !KTROC1 - int. in contraction
      LWRK1  = LWORK  - KEND1
C
      KINTOC = KEND1
      KEND1  = KINTOC  + NTOTOC(ISINT1) !KINTOC - temporary storage
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in CC3_AMATT3ZUVIR')
      END IF
C
      CALL DZERO(WORK(KRBJIAZU),NT2SQ(ISYRES))
C
C     Occupied integrals.
C
C-----------------------
C     Read in integrals for SMAT etc.
C-----------------------
C
      CALL INTOCC_T30(LUCKJD,FNCKJD,WORK(KLAMDP),ISINT2,WORK(KT3OG1),
     *                WORK(KT3OG2),WORK(KEND1),LWRK1)
C
C------------------------------------------
C     Z transformed Occupied integrals.
C-----------------------------------------
C
      CALL INTOCC_T3X(LUCKJDRZ,FNCKJDRZ,WORK(KLAMDP),ISINTRZ2,
     *                WORK(KW3ZOGZ1),WORK(KEND1),LWRK1)
C
C------------------------------------------
C     U transformed Occupied integrals.
C-----------------------------------------
C
      CALL INTOCC_T3X(LUCKJDRU,FNCKJDRU,WORK(KLAMDP),ISINTRU2,
     *                WORK(KW3UOGU1),WORK(KEND1),LWRK1)
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
         XNORMVAL  = DDOT(NTOTOC(ISINT1),WORK(KINTOC),1,WORK(KINTOC),1)
         WRITE(LUPRI,*) 'Norm of CCSDT_OC-INT ',XNORMVAL
      ENDIF
C
      !Transform (ia|j delta) integrals to (ia|j k) and sort as (i,j,k,a)
      CALL CCLR_TROCC(WORK(KINTOC),WORK(KTROC),WORK(KLAMDHX),ISYMRX,
     *                  WORK(KEND1),LWRK1)
C
      !sort (i,j,k,a) as (a,i,j,k)
      CALL CCSDT_SRTOC2(WORK(KTROC),WORK(KTROC1),ISINT1X,
     *                  WORK(KEND1),LWRK1)
c
C
C----------------------------
C     General loop structure.
C----------------------------
C
      DO ISYMD = 1,NSYM

         ISCKB1  = MULD2H(ISINT1,ISYMD)
         ISCKB1X  = MULD2H(ISINT1X,ISYMD)
         ISCKB2  = MULD2H(ISINT2,ISYMD)
         ISCKB2Z = MULD2H(ISINTRZ2,ISYMD)
         ISCKB2U = MULD2H(ISINTRU2,ISYMD)
C
         IF (IPRINT .GT. 55) THEN
C
            WRITE(LUPRI,*) 'In CC3_AMATT3ZUVIR: ISCKB2 :',ISCKB2
            WRITE(LUPRI,*) 'In CC3_AMATT3ZUVIR: ISCKB2Z:',ISCKB2Z
            WRITE(LUPRI,*) 'In CC3_AMATT3ZUVIR: ISCKB2U:',ISCKB2U

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

         KT3VDG1  = KEND2
         KEND3   = KT3VDG1  + NCKATR(ISCKB2)
         LWRK3   = LWORK  - KEND3
C
         KW3ZVDGZ1  = KEND3
         KEND3    = KW3ZVDGZ1  + NCKATR(ISCKB2Z)
         LWRK3    = LWORK    - KEND3
C
         KW3UVDGU1  = KEND3
         KEND3    = KW3UVDGU1  + NCKATR(ISCKB2U)
         LWRK3    = LWORK    - KEND3
C
         KT3VDGZ3 = KEND3
         KEND3    = KT3VDGZ3 + NCKATR(ISCKB2Z)
C
         KT3VDGU3 = KEND3
         KEND3    = KT3VDGU3 + NCKATR(ISCKB2U)
C
         KTRVI  = KEND3
         KTRVI1 = KTRVI  + NCKATR(ISCKB1X)   !KTRVI - int. in contraction
         KEND3 = KTRVI1 + NCKATR(ISCKB1X)    !KTRVI1 - int. in contraction
         LWRK3  = LWORK  - KEND3
C
         KINTVI  = KEND3
         MAXZU   = MAX(NCKA(ISCKB2Z),NCKA(ISCKB2U))
         KEND3   = KINTVI  + MAX(MAXZU,NCKA(ISCKB1))
         LWRK3   = LWORK   - KEND3
C
         IF (LWRK3 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND3
            CALL QUIT('Insufficient space in CC3_AMATT3ZUVIR')
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
     *                        WORK(KT3VDG1),WORK(KT3VDG2),
     *                        WORK(KT3VDG3),WORK(KLAMDH),ISYMD,D,
     *                        WORK(KEND3),LWRK3)
C
C-----------------------------------------------------------------------
C        Read Z transformed virtual integrals (H^Z) used for W^Z 
C-----------------------------------------------------------------------
C
            IOFF = ICKBD(ISCKB2Z,ISYMD) + NCKATR(ISCKB2Z)*(D - 1) + 1
            IF (NCKATR(ISCKB2Z) .GT. 0) THEN
               CALL GETWA2(LUDKBCRZ,FNDKBCRZ,WORK(KW3ZVDGZ1),IOFF,
     &                     NCKATR(ISCKB2Z))
            ENDIF
C
C-----------------------------------------------------------------------
C        Read U transformed virtual integrals (H^U) used for W^U 
C-----------------------------------------------------------------------
C
            IOFF = ICKBD(ISCKB2U,ISYMD) + NCKATR(ISCKB2U)*(D - 1) + 1
            IF (NCKATR(ISCKB2U) .GT. 0) THEN
               CALL GETWA2(LUDKBCRU,FNDKBCRU,WORK(KW3UVDGU1),IOFF,
     &                     NCKATR(ISCKB2U))
            ENDIF
C
C--------------------------------------------------------------------
C           Read virtual integrals [H,T1Z] where Z is LISTRZ (used in W^Z)
C--------------------------------------------------------------------
C
            IF (NCKA(ISCKB2Z) .GT. 0) THEN
               IOFF = ICKAD(ISCKB2Z,ISYMD) +
     &                NCKA(ISCKB2Z)*(D - 1) + 1
               CALL GETWA2(LUDELDRZ,FNDELDRZ,WORK(KINTVI),IOFF,
     *              NCKA(ISCKB2Z))
            ENDIF
C
            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KT3VDGZ3),
     *                       WORK(KLAMDH),ISYMD,D,ISINTRZ2,
     *                       WORK(KEND3),LWRK3)
C
C--------------------------------------------------------------------
C           Read virtual integrals [H,T1U] where U is LISTRU (used in W^U)
C--------------------------------------------------------------------
C
            IF (NCKA(ISCKB2U) .GT. 0) THEN
               IOFF = ICKAD(ISCKB2U,ISYMD) +
     &                NCKA(ISCKB2U)*(D - 1) + 1
               CALL GETWA2(LUDELDRU,FNDELDRU,WORK(KINTVI),IOFF,
     *              NCKA(ISCKB2U))
            ENDIF
C
            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KT3VDGU3),
     *                       WORK(KLAMDH),ISYMD,D,ISINTRU2,
     *                       WORK(KEND3),LWRK3)
C
C-----------------------------------------------------
C           Virtual integrals in pojection (fixed D)
C-----------------------------------------------------
C
            IOFF = ICKAD(ISCKB1,ISYMD) + NCKA(ISCKB1)*(D - 1) + 1
            IF (NCKA(ISCKB1) .GT. 0) THEN
               CALL GETWA2(LU3VI2,FN3VI2,WORK(KINTVI),IOFF,
     *                     NCKA(ISCKB1))
            ENDIF
            CALL CCLR_TRVIR(WORK(KINTVI),WORK(KTRVI),WORK(KLAMDPX),
     *                      ISYMRX,
     *                      ISYMD,D,ISYMOP,WORK(KEND3),LWRK3)
C
            IOFF = ICKAD(ISCKB1,ISYMD) + NCKA(ISCKB1)*(D - 1) + 1
            IF (NCKA(ISCKB1) .GT. 0) THEN
               CALL GETWA2(LU3VI,FN3VI,WORK(KINTVI),IOFF,
     *                     NCKA(ISCKB1))
            ENDIF
            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI1),WORK(KLAMDPX),
     *                      ISYMRX,
     *                      ISYMD,D,ISYMOP,WORK(KEND3),LWRK3)
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
               ISCKD2Z = MULD2H(ISINTRZ2,ISYMB)
               ISCKD2U = MULD2H(ISINTRU2,ISYMB)
               ISWMATZ  = MULD2H(ISYMRZ,ISCKIJ)
               ISWMATU  = MULD2H(ISYMRU,ISCKIJ)
               ISWMATZU  = MULD2H(ISWMATZ,ISYMRU)

               IF (IPRINT .GT. 55) THEN

                  WRITE(LUPRI,*) 'In CC3_AMATT3ZUVIR: ISYMD  :',ISYMD
                  WRITE(LUPRI,*) 'In CC3_AMATT3ZUVIR: ISYMB  :',ISYMB
                  WRITE(LUPRI,*) 'In CC3_AMATT3ZUVIR: ISYALJ :',ISYALJ
                  WRITE(LUPRI,*) 'In CC3_AMATT3ZUVIR: ISYALJ2:',ISYALJ2
                  WRITE(LUPRI,*) 'In CC3_AMATT3ZUVIR: ISYMBD :',ISYMBD
                  WRITE(LUPRI,*) 'In CC3_AMATT3ZUVIR: ISCKIJ :',ISCKIJ
                  WRITE(LUPRI,*) 'In CC3_AMATT3ZUVIR: ISWMATZ :',ISWMATZ
                  WRITE(LUPRI,*) 'In CC3_AMATT3ZUVIR: ISWMATU :',ISWMATU

               ENDIF
C
               KSMAT   = KEND3
               KSMAT3  = KSMAT   + NCKIJ(ISCKIJ)
               KUMAT   = KSMAT3  + NCKIJ(ISCKIJ)
               KUMAT3  = KUMAT   + NCKIJ(ISCKIJ)
               KDIAG   = KUMAT3  + NCKIJ(ISCKIJ)
               KT3MAT   = KDIAG   + NCKIJ(ISCKIJ)
               KEND3   = KT3MAT   + NCKIJ(ISCKIJ)
C
               KDIAGWZ  = KEND3
               KWMATZ   = KDIAGWZ  + NCKIJ(ISWMATZ)
               KINDSQWZ = KWMATZ   + NCKIJ(ISWMATZ)
               KEND3  = KINDSQWZ + (6*NCKIJ(ISWMATZ) - 1)/IRAT + 1
C
               KDIAGWU  = KEND3
               KWMATU   = KDIAGWU  + NCKIJ(ISWMATU)
               KINDSQWU = KWMATU   + NCKIJ(ISWMATU)
               KEND3  = KINDSQWU + (6*NCKIJ(ISWMATU) - 1)/IRAT + 1
C
               KINDSQ  = KEND3
               KINDEX  = KINDSQ  + (6*NCKIJ(ISCKIJ) - 1)/IRAT + 1
               KINDEX2 = KINDEX  + (NCKI(ISYALJ) - 1)/IRAT + 1
               KEND3   = KINDEX2 + (NCKI(ISYALJ2) - 1)/IRAT + 1
C
               MAXZU    = MAX(NCKIJ(ISWMATZ),NCKIJ(ISWMATU))
               MAXZZUU    = MAX(MAXZU,NCKIJ(ISWMATZU))
C
               KTMATW   = KEND3
               KT3VBG1  = KTMATW   + MAX(MAXZZUU,NCKIJ(ISCKIJ))
               KT3VBG2  = KT3VBG1  + NCKATR(ISCKD2)
               KT3VBG3 = KT3VBG2  + NCKATR(ISCKD2)
               KEND4   = KT3VBG3 + NCKATR(ISCKD2)
               LWRK4   = LWORK   - KEND4
C
               KW3ZVDGZ2 = KEND4
               KEND4   = KW3ZVDGZ2 + NCKATR(ISCKD2Z)
C
               KW3UVDGU2 = KEND4
               KEND4   = KW3UVDGU2 + NCKATR(ISCKD2U)
C
               KT3VBGZ3 = KEND4
               KEND4    = KT3VBGZ3 + NCKATR(ISCKD2Z)
C
               KT3VBGU3 = KEND4
               KEND4    = KT3VBGU3 + NCKATR(ISCKD2U)
C
               KINTVI  = KEND4
               KEND4   = KINTVI  + MAX(NCKA(ISCKD2Z),NCKA(ISCKD2U))
               LWRK4   = LWORK   - KEND4
C
               KINDSQWZU = KEND4
               KDIAGWZU  = KINDSQWZU + (6*NCKIJ(ISWMATZU) - 1)/IRAT + 1
               KWMATZU   = KDIAGWZU  + NCKIJ(ISWMATZU)
               KEND4     = KWMATZU   + NCKIJ(ISWMATZU)
               LWRK4     = LWORK     - KEND4
C
               IF (LWRK4 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Memory available : ',LWORK
                  WRITE(LUPRI,*) 'Memory needed    : ',KEND4
                  CALL QUIT('Insufficient space in CC3_AMATT3ZUVIR')
               END IF
C
C---------------------------------------------
C           Construct part of the diagonal.
C---------------------------------------------
C
               CALL CC3_DIAG(WORK(KDIAG),WORK(KFOCKD),ISCKIJ)
               CALL CC3_DIAG(WORK(KDIAGWZ),WORK(KFOCKD),ISWMATZ)
               CALL CC3_DIAG(WORK(KDIAGWU),WORK(KFOCKD),ISWMATU)
               CALL CC3_DIAG(WORK(KDIAGWZU),WORK(KFOCKD),ISWMATZU)
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
               LENSQWZ = NCKIJ(ISWMATZ)
               CALL CC3_INDSQ(WORK(KINDSQWZ),LENSQWZ,ISWMATZ) 
               LENSQWU = NCKIJ(ISWMATU)
               CALL CC3_INDSQ(WORK(KINDSQWU),LENSQWU,ISWMATU) 
               LENSQWZU = NCKIJ(ISWMATZU)
               CALL CC3_INDSQ(WORK(KINDSQWZU),LENSQWZU,ISWMATZU) 


               DO B = 1,NVIR(ISYMB)
C
C------------------------------------
C           Initialise
C---------------------------------------
C
               CALL DZERO(WORK(KWMATZ),NCKIJ(ISWMATZ))
               CALL DZERO(WORK(KWMATU),NCKIJ(ISWMATU))
               CALL DZERO(WORK(KWMATZU),NCKIJ(ISWMATZU))
C
C-----------------------------------------------
C           Read virtual integrals used in second s3am.
C-----------------------------------------------
C
                  CALL INTVIR_T30_D(LUDKBC,FNDKBC,LUDELD,FNDELD,ISINT2,
     *                              WORK(KT3VBG1),WORK(KT3VBG2),
     *                              WORK(KT3VBG3),WORK(KLAMDH),ISYMB,B,
     *                              WORK(KEND4),LWRK4)
C
C--------------------------------------------------------------------
C           Read Z transformed virtual integrals (H^Z) used in W^Z
C--------------------------------------------------------------------
C
                  IOFF = ICKBD(ISCKD2Z,ISYMB) + 
     &                   NCKATR(ISCKD2Z)*(B - 1) + 1
                  IF (NCKATR(ISCKD2Z) .GT. 0) THEN
                     CALL GETWA2(LUDKBCRZ,FNDKBCRZ,WORK(KW3ZVDGZ2),IOFF,
     &                           NCKATR(ISCKD2Z))
                  ENDIF
C
C--------------------------------------------------------------------
C           Read U transformed virtual integrals (H^U) used in W^U
C--------------------------------------------------------------------
C
                  IOFF = ICKBD(ISCKD2U,ISYMB) + 
     &                   NCKATR(ISCKD2U)*(B - 1) + 1
                  IF (NCKATR(ISCKD2U) .GT. 0) THEN
                     CALL GETWA2(LUDKBCRU,FNDKBCRU,WORK(KW3UVDGU2),IOFF,
     &                           NCKATR(ISCKD2U))
                  ENDIF
C
C--------------------------------------------------------------------
C           Read virtual integrals [H,T1Z] where Z is LISTRZ (used in W^Z)
C--------------------------------------------------------------------
C
                  IF (NCKA(ISCKD2Z) .GT. 0) THEN
                     IOFF = ICKAD(ISCKD2Z,ISYMB) +
     &                      NCKA(ISCKD2Z)*(B - 1) + 1
                     CALL GETWA2(LUDELDRZ,FNDELDRZ,WORK(KINTVI),IOFF,
     *                    NCKA(ISCKD2Z))
                  ENDIF
C
                  CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KT3VBGZ3),
     *                             WORK(KLAMDH),ISYMB,B,ISINTRZ2,
     *                             WORK(KEND4),LWRK4)
C
C--------------------------------------------------------------------
C           Read virtual integrals [H,T1U] where U is LISTRU (used in W^U)
C--------------------------------------------------------------------
C
                  IF (NCKA(ISCKD2U) .GT. 0) THEN
                     IOFF = ICKAD(ISCKD2U,ISYMB) +
     &                      NCKA(ISCKD2U)*(B - 1) + 1
                     CALL GETWA2(LUDELDRU,FNDELDRU,WORK(KINTVI),IOFF,
     *                    NCKA(ISCKD2U))
                  ENDIF
C
                  CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KT3VBGU3),
     *                             WORK(KLAMDH),ISYMB,B,ISINTRU2,
     *                             WORK(KEND4),LWRK4)
C
C------------------------------------------------------------------------
C           Calculate the S(ci,bk,dj) matrix for T3 for B,D.
C-------------------------------------------------------------------
C
                  CALL GET_T30_BD(ISYMT3,ISINT2,WORK(KT2TP),ISYMT2,
     *                            WORK(KT3MAT),WORK(KFOCKD),WORK(KDIAG),
     *                            WORK(KINDSQ),LENSQ,WORK(KSMAT),
     *                            WORK(KT3VDG1),WORK(KT3VDG2),
     *                            WORK(KT3OG1),WORK(KINDEX),
     *                            WORK(KSMAT3),WORK(KT3VBG1),
     *                            WORK(KT3VBG2),WORK(KINDEX2),
     *                            WORK(KUMAT),WORK(KT3VDG3),
     *                            WORK(KT3OG2),WORK(KUMAT3),
     *                            WORK(KT3VBG3),ISYMB,B,ISYMD,D,
     *                            ISCKIJ,WORK(KEND4),LWRK4)

c
c       call sum_pt3(work(KT3MAT),isymb,b,isymd,d,
c    *             1,work(kx3am),3)

C
C------------------------------------------------------
C Calculate WBD^Z
C------------------------------------------------------
C

                  ! Copy T3^0 from KT3MAT to KTMATW. KTMATW used later
                  ! as help array and KTMATW overwritten.
c                 CALL DCOPY(NCKIJ(ISCKIJ),WORK(KT3MAT),1,
c    *                       WORK(KTMATW),1)
C

C
C------------------------------------------------------
C     Calculate the  term <mu3|[Z,T30]|HF> occupied contribution 
C     added in W^BD(a,i-,k-,j-), where i- means one-index transformation of i
C     (the three contributions are added: fixed BD, fixed aB, fixed Da 
C      added in W^BD).
C------------------------------------------------------
C
                  AIBJCK_PERM = 4 
C     write(lupri,*)'entering WXBD_O... '
                  CALL WXBD_O(AIBJCK_PERM,WORK(KT3MAT),ISCKIJ,
     *                        WORK(KFCKZ),ISYMRZ,
     *                        WORK(KWMATZ),ISWMATZ,WORK(KEND4),LWRK4)
C
C------------------------------------------------------
C     Calculate the  term <mu3|[[Z,T2],T2]|HF> 
C     added in W^BD(a,i-,k-,j-).
C------------------------------------------------------
C
C     write(lupri,*)'entering WXBD_T2... '
                  CALL WXBD_T2(AIBJCK_PERM,B,ISYMB,D,ISYMD,WORK(KT2TP),
     *                         ISYMT2,WORK(KFCKZ),
     *                         ISYMRZ,WORK(KINDSQWZ),LENSQWZ,
     *                         WORK(KWMATZ),ISWMATZ,WORK(KEND4),LWRK4)

C
C----------------------------------------------------------------
C     Calculate the  terms <mu3|[H^Z,T2]|HF> + <mu3|[H,T2^Z]|HF>
C     added in W^BD(a,i-,k-,j-).
C----------------------------------------------------------------
C
C     write(lupri,*)'entering WXBD_GROUND(1)... '
                  CALL WXBD_GROUND(AIBJCK_PERM,
     *                       WORK(KT2Z),ISYMRZ,WORK(KTMATW),
     *                       WORK(KT3VDG1),WORK(KT3VBG1),WORK(KT3VBG3),
     *                       WORK(KT3VDG3),
     *                       WORK(KT3OG1),ISINT2,
     *                       WORK(KWMATZ),WORK(KEND4),LWRK4,
     *                       WORK(KINDSQWZ),LENSQWZ,
     *                       ISYMB,B,ISYMD,D)
C
C     write(lupri,*)'entering WXBD_GROUND(2)... '
                  CALL WXBD_GROUND(AIBJCK_PERM,
     *                       WORK(KT2TP),ISYMT2,WORK(KTMATW),
     *                       WORK(KW3ZVDGZ1),WORK(KW3ZVDGZ2),
     *                       WORK(KT3VBGZ3),WORK(KT3VDGZ3),
     *                       WORK(KW3ZOGZ1),ISINTRZ2,
     *                       WORK(KWMATZ),WORK(KEND4),LWRK4,
     *                       WORK(KINDSQWZ),LENSQWZ,
     *                       ISYMB,B,ISYMD,D)
C
C------------------------------------------------
C     Divide by the energy difference and
C     remove the forbidden elements
C------------------------------------------------
C
                  CALL WBD_DIA(B,ISYMB,D,ISYMD,FREQRZ,ISWMATZ,
     *                         WORK(KWMATZ),WORK(KDIAGWZ),WORK(KFOCKD))
                  CALL T3_FORBIDDEN(WORK(KWMATZ),ISYMRZ,ISYMB,B,
     *                              ISYMD,D)

c
c       call sum_pt3(work(KWMATZ),isymb,b,isymd,d,
c    *             ISWMATZ,work(kx3am),5)


C
c                 CALL GET_T3X_BD(ISYMRZ,WORK(KWMATZ),ISWMATZ,
c    *                            WORK(KTMATW),ISCKIJ,WORK(KFCKZ),
c    *                            ISYMRZ,WORK(KT2TP),ISYMT2,
c    *                            WORK(KT2Z),ISYMRZ,WORK(KT3VDG1),
c    *                            WORK(KT3VBG1),WORK(KT3OG1),
c    *                            ISINT2,WORK(KW3ZVDGZ1),
c    *                            WORK(KW3ZVDGZ2),WORK(KW3ZOGZ1),
c    *                            ISINTRZ2,FREQRZ,WORK(KDIAGWZ),
c    *                            WORK(KFOCKD),WORK(KINDSQWZ),LENSQWZ,
c    *                            B,ISYMB,D,ISYMD,WORK(KEND4),LWRK4)

C                 Calculate WUZ intermmediate for <mu3|[U,wZ^{aBC}_{ijk}]|HF>
C                 where wZ^{aBC}_{ijk} is part of T3^Z with the virtual
C                 contribution (WXBD_V routine) taken out.

c
                  CALL WBD_O(WORK(KWMATZ),ISWMATZ,WORK(KFCKU),ISYMRU,
     *                 WORK(KWMATZU),ISWMATZU,WORK(KEND4),LWRK4)
C
                  CALL WBD_V(WORK(KWMATZ),ISWMATZ,WORK(KFCKU),ISYMRU,
     *                 WORK(KWMATZU),ISWMATZU,WORK(KEND4),LWRK4)
C
C                 Calculate <mu3|[[U,T1Z],T30]|HF>
C
                  CALL WBD_O(WORK(KT3MAT),ISCKIJ,WORK(KFCKUZO),ISYMZU,
     *                 WORK(KWMATZU),ISWMATZU,WORK(KEND4),LWRK4)
C
                  CALL WBD_V(WORK(KT3MAT),ISCKIJ,WORK(KFCKUZV),ISYMZU,
     *                 WORK(KWMATZU),ISWMATZU,WORK(KEND4),LWRK4)
C
C                 Calculate <mu3|[[U,T2Z],T20]|HF>
C
                  T2XNET2Z = .TRUE.
                  CALL WBD_T2(T2XNET2Z,B,ISYMB,D,ISYMD,
     *                        WORK(KT2Z),ISYMRZ,WORK(KT2TP),ISYMT2,
     *                        WORK(KFCKU),ISYMRU,
     *                        WORK(KINDSQWZU),LENSQWZU,WORK(KWMATZU),
     *                        ISWMATZU,WORK(KEND4),LWRK4)
C
C------------------------------------------------
C     Divide by the energy difference and
C     remove the forbidden elements (here only for debugging)
C------------------------------------------------
C
c                 CALL WBD_DIA(B,ISYMB,D,ISYMD,FREQZU,ISWMATZU,
c    *                        WORK(KWMATZU),WORK(KDIAGWZU),WORK(KFOCKD))
c                 CALL T3_FORBIDDEN(WORK(KWMATZU),ISYMZU,ISYMB,B,
c    *                              ISYMD,D)

 
c       call sum_pt3(work(KWMATZU),isymb,b,isymd,d,
c    *             ISWMATZU,work(kx3am),5)


C
C------------------------------------------------------
C Calculate WBD^U
C------------------------------------------------------
C

                  ! Copy T3^0 from KT3MAT to KTMATW. KTMATW used later
                  ! as help array and KTMATW overwritten.
c                 CALL DCOPY(NCKIJ(ISCKIJ),WORK(KT3MAT),1,
c    *                       WORK(KTMATW),1)
C
C
C------------------------------------------------------
C     Calculate the  term <mu3|[U,T30]|HF> occupied contribution 
C     added in W^BD(a,i-,k-,j-), where i- means one-index transformation of i
C     (the three contributions are added: fixed BD, fixed aB, fixed Da 
C      added in W^BD).
C------------------------------------------------------
C
                  AIBJCK_PERM = 4 ! transform all occupied indecies
c     write(lupri,*)'entering WXBD_O... U'
                  CALL WXBD_O(AIBJCK_PERM,WORK(KT3MAT),ISCKIJ,
     *                        WORK(KFCKU),ISYMRU,
     *                        WORK(KWMATU),ISWMATU,WORK(KEND4),LWRK4)
C
C------------------------------------------------------
C     Calculate the  term <mu3|[[U,T2],T2]|HF> 
C     added in W^BD(a,i-,k-,j-).
C------------------------------------------------------
C
c     write(lupri,*)'entering WXBD_T2... U'
                  CALL WXBD_T2(AIBJCK_PERM,B,ISYMB,D,ISYMD,WORK(KT2TP),
     *                         ISYMT2,WORK(KFCKU),
     *                         ISYMRU,WORK(KINDSQWU),LENSQWU,
     *                         WORK(KWMATU),ISWMATU,WORK(KEND4),LWRK4)

C
C----------------------------------------------------------------
C     Calculate the  terms <mu3|[H^U,T2]|HF> + <mu3|[H,T2^U]|HF>
C     added in W^BD(a,i-,k-,j-).
C----------------------------------------------------------------
C
c     write(lupri,*)'entering WXBD_GROUND(1)... U'
                  CALL WXBD_GROUND(AIBJCK_PERM,
     *                       WORK(KT2U),ISYMRU,WORK(KTMATW),
     *                       WORK(KT3VDG1),WORK(KT3VBG1),WORK(KT3VBG3),
     *                       WORK(KT3VDG3),
     *                       WORK(KT3OG1),ISINT2,
     *                       WORK(KWMATU),WORK(KEND4),LWRK4,
     *                       WORK(KINDSQWU),LENSQWU,
     *                       ISYMB,B,ISYMD,D)
C
c     write(lupri,*)'entering WXBD_GROUND(2)... U'
                  CALL WXBD_GROUND(AIBJCK_PERM,
     *                       WORK(KT2TP),ISYMT2,WORK(KTMATW),
     *                       WORK(KW3UVDGU1),WORK(KW3UVDGU2),
     *                       WORK(KT3VBGU3),WORK(KT3VDGU3),
     *                       WORK(KW3UOGU1),ISINTRU2,
     *                       WORK(KWMATU),WORK(KEND4),LWRK4,
     *                       WORK(KINDSQWU),LENSQWU,
     *                       ISYMB,B,ISYMD,D)
C
C------------------------------------------------
C     Divide by the energy difference and
C     remove the forbidden elements
C------------------------------------------------
C
                  CALL WBD_DIA(B,ISYMB,D,ISYMD,FREQRU,ISWMATU,
     *                         WORK(KWMATU),WORK(KDIAGWU),WORK(KFOCKD))
                  CALL T3_FORBIDDEN(WORK(KWMATU),ISYMRU,ISYMB,B,
     *                              ISYMD,D)

c
c       call sum_pt3(work(KWMATU),isymb,b,isymd,d,
c    *             ISWMATU,work(kx3am),5)


c                 CALL GET_T3X_BD(ISYMRU,WORK(KWMATU),ISWMATU,
c    *                            WORK(KTMATW),ISCKIJ,WORK(KFCKU),
c    *                            ISYMRU,WORK(KT2TP),ISYMT2,
c    *                            WORK(KT2U),ISYMRU,WORK(KT3VDG1),
c    *                            WORK(KT3VBG1),WORK(KT3OG1),
c    *                            ISINT2,WORK(KW3UVDGU1),
c    *                            WORK(KW3UVDGU2),WORK(KW3UOGU1),
c    *                            ISINTRU2,FREQRU,WORK(KDIAGWU),
c    *                            WORK(KFOCKD),WORK(KINDSQWU),LENSQWU,
c    *                            B,ISYMB,D,ISYMD,WORK(KEND4),LWRK4)

C                 Calculate WZU intermmediate for <mu3|[Z,wY^{aBC}_{ijk}]|HF>
C                 where wY^{aBC}_{ijk} is part of T3^U with the virtual
C                 contribution (WXBD_V routine) taken out.

                  CALL WBD_O(WORK(KWMATU),ISWMATU,WORK(KFCKZ),ISYMRZ,
     *                 WORK(KWMATZU),ISWMATZU,WORK(KEND4),LWRK4)
C
                  CALL WBD_V(WORK(KWMATU),ISWMATU,WORK(KFCKZ),ISYMRZ,
     *                 WORK(KWMATZU),ISWMATZU,WORK(KEND4),LWRK4)
C
C                 Calculate <mu3|[[Z,T1U],T30]|HF>
C
                  CALL WBD_O(WORK(KT3MAT),ISCKIJ,WORK(KFCKZUO),ISYMZU,
     *                 WORK(KWMATZU),ISWMATZU,WORK(KEND4),LWRK4)
C
                  CALL WBD_V(WORK(KT3MAT),ISCKIJ,WORK(KFCKZUV),ISYMZU,
     *                 WORK(KWMATZU),ISWMATZU,WORK(KEND4),LWRK4)
C
C                 Calculate <mu3|[[Z,T2U],T20]|HF>
C
                  T2XNET2Z = .TRUE.
                  CALL WBD_T2(T2XNET2Z,B,ISYMB,D,ISYMD,
     *                        WORK(KT2U),ISYMRU,WORK(KT2TP),ISYMT2,
     *                        WORK(KFCKZ),ISYMRZ,
     *                        WORK(KINDSQWZU),LENSQWZU,WORK(KWMATZU),
     *                        ISWMATZU,WORK(KEND4),LWRK4)
C
C------------------------------------------------
C     Divide by the energy difference and
C     remove the forbidden elements 
C------------------------------------------------
C
                  CALL T3_FORBIDDEN(WORK(KWMATZU),ISYMZU,ISYMB,B,
     *                              ISYMD,D)

c                 call sum_pt3(work(KWMATZU),isymb,b,isymd,d,
c    *                       ISWMATZU,work(kx3am),4)

                  CALL WBD_DIA(B,ISYMB,D,ISYMD,FREQZU,ISWMATZU,
     *                        WORK(KWMATZU),WORK(KDIAGWZU),WORK(KFOCKD))
 

C
C-----------------------------------------------
C                 Carry out the contractions...
C-----------------------------------------------
C
C
C------------------------------------------------------
C                 Calculate the  term <mu2|[H,WBD^ZU(3)]|HF> 
C                 ( Fock matrix cont ) 
C                 added in OMEGA2EFF 
C------------------------------------------------------
C
                  CALL CC3_CY2F(OMEGA2EFF,ISYRES,WORK(KWMATZU),ISWMATZU,
     *                          WORK(KTMATW),WORK(KFOCKX),ISYFCKX,
     *                          WORK(KINDSQWZU),
     *                          LENSQWZU,WORK(KEND4),LWRK4,
     *                          ISYMB,B,ISYMD,D)
c                 
C------------------------------------------------------
C                 Calculate the  term <mu2|[H,WBD^ZU(3)]|HF> 
C                 ( Occupied  cont ) 
C                 added in OMEGA2EFF 
C------------------------------------------------------
C
                 CALL CC3_CY2O(OMEGA2EFF,ISYRES,WORK(KWMATZU),ISWMATZU,
     *                          WORK(KTMATW),WORK(KTROC),
     *                          WORK(KTROC1),ISINT1,WORK(KEND4),LWRK4,
     *                          WORK(KINDSQWZU),LENSQWZU,
     *                          ISYMB,B,ISYMD,D)
C
C------------------------------------------------------
C                 Calculate the  term <mu2|[H,WBD^ZU(3)]|HF> 
C                 ( Virtual  cont ) 
C                 added in OMEGA2EFF 
C------------------------------------------------------
C
                  CALL CC3_CY2V(OMEGA2EFF,ISYRES,WORK(KRBJIAZU),
     *                          WORK(KWMATZU),
     *                          ISWMATZU,WORK(KTMATW),WORK(KTRVI),
     *                          WORK(KTRVI1),ISINT1,WORK(KEND4),LWRK4,
     *                          WORK(KINDSQWZU),LENSQWZU,
     *                          ISYMB,B,ISYMD,D)
C
               END DO !B
            END DO    !ISYMB
C
         END DO       !D
      END DO          !ISYMD
c
c      
c      write(lupri,*) 'KWMATZ in CC3_AAMAT  : '
c      call print_pt3(work(kx3am),1,4)
c

C
C
C------------------------------------------------------
C     Accumulate RBJIA from <mu2|[H,WBD^ZU(3)]|HF> ( Virtual  cont ) 
C     in OMEGA2EFF 
C------------------------------------------------------
C
      CALL CC3_RBJIA(OMEGA2EFF,ISYRES,WORK(KRBJIAZU))
C
      IF (IPRINT .GT. 55) THEN
         XNORMVAL = DDOT(NT2AM(ISYRES),OMEGA2EFF,1,OMEGA2EFF,1)
         WRITE(LUPRI,*)'Norm of OMEGA2EFF final ',XNORMVAL
      ENDIF
C
C
C
c      write(lupri,*) 'cont to t3x in CC3_AMATT3ZUVIR'
c      call print_pt3(work(kx3am),1,4)
C      write(lupri,*) 'The summed S terms : '
C      call print_pt3(work(kx3am),1,1)
C
C--------------------------------
C     Close files for "response"
C--------------------------------
C
      CALL WCLOSE2(LU3SRTR,FN3SRTR,'DELETE')
      CALL WCLOSE2(LUCKJDRZ,FNCKJDRZ,'DELETE')
      CALL WCLOSE2(LUCKJDRU,FNCKJDRU,'DELETE')
      CALL WCLOSE2(LUDELDRZ,FNDELDRZ,'DELETE')
      CALL WCLOSE2(LUDELDRU,FNDELDRU,'DELETE')
      CALL WCLOSE2(LUDKBCRZ,FNDKBCRZ,'DELETE')
      CALL WCLOSE2(LUDKBCRU,FNDKBCRU,'DELETE')
C
C-------------
C     End
C-------------
C
      CALL QEXIT('AT3ZUV')
C
      RETURN
      END
C  /* Deck cc3_amatt3zuocc */
      SUBROUTINE CC3_AMATT3ZUOCC(LISTX,IDLSTX,LISTZU,IDLSTZU,
     &                          OMEGA2,
     &                          OMEGA2EFF,
     *                          ISYRES,WORK,LWORK,
     *                          LUCKJD,FNCKJD,LUTOC,FNTOC,
     *                          LUDELD,FNDELD,LU3FOP,FN3FOP,
     *                          LU3VI,FN3VI)
*
*******************************************************************
*
* This routine calculates the contribution to cubic response 
* (originating form both F matrix and B matrix):
*
* omega2eff = <mu2|[H^X,T3^ZU]|HF>
*
* which is contracted outside with the multipliers (zeroth-order for 
* F matrix or first-order for B matrix).
*
*
*
* T3^ZU is calculated using W^JK intermmediate.
*
* !!! ACTUALLY HERE WE CALCULATE ONLY THE PART OF THE TERM, WHICH ORIGINATES
*     FROM THE A^U-MATRIX CONTRIBUTION TO T3^ZU 
*     (to get the rest of contribution call cc3_bmatt3zu routine) !!!
*
* Thus here:
*
* W^BD =    <mu3|[Z,T3^U]|HF>
*
* (!!! TO GET THE TOTAL CONTRIBUTION FROM THIS TERM AND THE OTHER 
*       A^U-MATRIX TERMS CALL CC3_AMATT3ZUVIR !!!)
*
*
* F. Pawlowski, P. Jorgensen, Aarhus Spring 2003.
*
**************************************************************
C
      IMPLICIT NONE
#include <priunit.h>
#include <ccsdsym.h>
#include <ccl1rsp.h>
#include <ccr1rsp.h>
#include <ccorb.h>
#include <dummy.h>
#include <iratdef.h>
#include <ccinftap.h>
#include <ccsdinp.h>
C
      CHARACTER*(*) FNCKJD,FNTOC,FNDELD,FN3FOP,FN3VI
      INTEGER LUCKJD,LUTOC,LUDELD,LU3FOP,LU3VI
C
      CHARACTER CDUMMY*1
C
      CHARACTER*3 LISTX, LISTRZU, LISTRZ, LISTRU
      CHARACTER*8 LABELRZ,LABELRU
      INTEGER     IDLSTX,IDLSTRZU,IDLSTRZ,IDLSTRU
      INTEGER     ISYMRZ, ISYMRU
C
      INTEGER ISYRES, LWORK
      INTEGER ISYM0,ISYMT2,ISYMT3
      INTEGER KT2TP,KFOCKD,KLAMDP,KLAMDH,KEND00,LWRK00,KXIAJB,KFOCK0
      INTEGER IOPTTCME,IOPT
      INTEGER ISYMWZ,ISYMWU,ISYMWZU
      INTEGER ISYFCKZ,ISYFCKU,ISYMZU,KFCKZ,KFCKU,KEND1,LWRK1
      INTEGER ISINT1,ISINT2,KT3BOG2,KT3OG1
      INTEGER KT3OG2,KT3VIJG1,KXGADCK,KXLADCK
      INTEGER ISYMTETAZ,ISYMTETAU,ISYMTETAZU
      INTEGER ISYMK,ISYML,ISYMKL,ISYT30KL,ISTETAZKL,ISTETAUKL,ISTETAZUKL
      INTEGER KT30KL,KTETAZKL,KTETAUKL,KTETAZUKL,KTMAT
c
      INTEGER ISYMRX,ISYFCKX,KLAMDPX,KLAMDHX,KFOCKX,KT1X,ISYMC
      INTEGER KOFF1,KOFF2,ISINT2X,KT3BOG1,KT3BOL2,KINTOC,IOFF
c
      integer kx3am
c
#if defined (SYS_CRAY)
      REAL OMEGA2(*),OMEGA2EFF(*)
      REAL WORK(LWORK)
      REAL FREQRZ,FREQRU,FREQZU
      REAL XNORMVAL,DDOT
#else
      DOUBLE PRECISION OMEGA2(*),OMEGA2EFF(*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION FREQRZ,FREQRU,FREQZU
      DOUBLE PRECISION XNORMVAL,DDOT
#endif
C
      CALL QENTER('AT3ZUOCC')
C
C--------------------------
C     Save and set flags.
C--------------------------
C
      IPRCC   = IPRINT
      ISYM0   = 1
      ISINT1  = 1
      ISINT2  = 1
      ISYMT2  = 1
      ISYMT3  = MULD2H(ISINT2,ISYMT2)
C
C----------------------------------------------
C     Calculate the zeroth order stuff once
C----------------------------------------------
C
      KT2TP  = 1 
      KFOCKD  = KT2TP  + NT2SQ(ISYMT2)
      KLAMDP = KFOCKD + NORBTS
      KLAMDH = KLAMDP + NLAMDT
      KEND00 = KLAMDH + NLAMDT
      LWRK00 = LWORK  - KEND00
C
      KXIAJB = KEND00
      KEND00 = KXIAJB + NT2AM(ISINT1)
      LWRK00 = LWORK  - KEND00
C
      IF (LWRK00 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_AMATT3ZUOCC (zeroth allo.')
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
            call quit('Insufficient space (T3) in CC3_AMATT3ZUOCC')
         END IF
      endif
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
C
C determining R1 labels
C
      ISYMRX = ISYLRT(IDLSTX)
      !FREQuency not needed for LISTX
      !LABELX = LRTLBL(IDLSTX) not needed for LISTX
C
      IF (LISTZU(1:3).EQ.'R2 ') THEN
        LISTRZ  = 'R1 '
        LABELRZ = LBLR2T(IDLSTZU,1)
        ISYMRZ  = ISYR2T(IDLSTZU,1)
        FREQRZ  = FRQR2T(IDLSTZU,1)
        LORXRZ  = LORXR2T(IDLSTZU,1)
        IDLSTRZ = IR1TAMP(LABELRZ,LORXRZ,FREQRZ,ISYMRZ)
C
        IF (LORXRZ) CALL QUIT('NO ORBITAL RELAX. IN CC3_AMATT3ZUOCC')

        LISTRU  = 'R1 '
        LABELRU = LBLR2T(IDLSTZU,2)
        ISYMRU  = ISYR2T(IDLSTZU,2)
        FREQRU  = FRQR2T(IDLSTZU,2)
        LORXRU  = LORXR2T(IDLSTZU,2)
        IDLSTRU = IR1TAMP(LABELRU,LORXRU,FREQRU,ISYMRU)
C
        IF (LORXRU) CALL QUIT('NO ORBITAL RELAX. IN CC3_AMATT3ZUOCC')
C
      ELSE
         CALL QUIT('Unknown list for LISTZU in CC3_AMATT3ZUOCC.')
      END IF
C
      FREQZU = FREQRZ + FREQRU
C
      ISYMWZ   = MULD2H(ISYMT3,ISYMRZ)
      ISYMWU   = MULD2H(ISYMT3,ISYMRU)
C
      ISYMWZU   = MULD2H(ISYMWZ,ISYMRU)
C
      IF (ISYRES.NE.ISYMWZU) THEN
       WRITE(LUPRI,*) 'INCONSISTENT SYMMETRY SPECIFICATION'
       WRITE(LUPRI,*) 
     *      'ISYRES =',ISYRES,'ISYMWZU =',ISYMWZU
       CALL QUIT('CC3_AMATT3ZUOCC inconsistent symmetry specification') 
      END IF
C
      ISYFCKZ = ISYMRZ
      ISYFCKU = ISYMRU
      ISYMZU = MULD2H(ISYMRZ,ISYMRU)
      ISYFCKX = MULD2H(ISYMOP,ISYMRX)
C
C-------------------------------------------------
C     Read property matrices and Lambda matrices
C-------------------------------------------------
C
      KFCKZ  = KEND00
      KFCKU  = KFCKZ  + N2BST(ISYFCKZ)
      KEND1  = KFCKU  + N2BST(ISYFCKU)
      LWRK1  = LWORK  - KEND1
C
      KLAMDPX = KEND1
      KLAMDHX = KLAMDPX + NLAMDT
      KEND1   = KLAMDHX + NLAMDT
      LWRK1   = LWORK   - KEND1
C
      KFOCKX  = KEND1
      KEND1   = KFOCKX  + N2BST(ISYFCKX)
      LWRK1   = LWORK  - KEND1
C
      KT1X   = KEND1
      KEND1  = KT1X   + NT1AM(ISYMRX)
      LWRK1  = LWORK  - KEND1
C        
      IF (LWRK1 .LT. NT2SQ(ISYMRZ)) THEN
         CALL QUIT('Out of memory in CC3_AMATT3ZUOCC (TOT_T3Y) ')
      ENDIF
C
C---------------------------
C     Read Lambda_X matrices
C---------------------------
C
      CALL GET_LAMBDAX(WORK(KLAMDPX),WORK(KLAMDHX),LISTX,IDLSTX,
     *                 ISYMRX,WORK(KLAMDP),WORK(KLAMDH),WORK(KEND1),
     *                 LWRK1)
C
C------------------
C     Read in T1^X 
C------------------
C
      IOPT = 1
      CALL GET_T1_T2(IOPT,.FALSE.,WORK(KT1X),DUMMY,LISTX,IDLSTX,
     *               ISYMRX,WORK(KEND1),LWRK1)
C
C------------------------------------------
C        Calculate the F^X matrix
C------------------------------------------
C
      CALL CC3LR_MFOCK(WORK(KEND1),WORK(KT1X),WORK(KXIAJB),ISYFCKX)
C
C     Put the F_{kc} part into correct F_{pq}
C
         CALL DZERO(WORK(KFOCKX),N2BST(ISYFCKX))
C
         DO ISYMC = 1, NSYM
            ISYMK = MULD2H(ISYFCKX,ISYMC)
            DO K = 1, NRHF(ISYMK)
               DO C = 1, NVIR(ISYMC)
                  KOFF1 = KFOCKX - 1
     *                  + IFCVIR(ISYMK,ISYMC)
     *                  + NORB(ISYMK)*(C - 1)
     *                  + K
                  KOFF2 = KEND1 - 1
     *                  + IT1AM(ISYMC,ISYMK)
     *                  + NVIR(ISYMC)*(K - 1)
     *                  + C
C
                  WORK(KOFF1) = WORK(KOFF2)
C
               ENDDO
            ENDDO
         ENDDO

C
C------------------------------------------
C     Calculate the F^Z matrix
C------------------------------------------
C
      CALL GET_FOCKX(WORK(KFCKZ),LABELRZ,IDLSTRZ,ISYMRZ,WORK(KLAMDP),
     *                  ISYM0,WORK(KLAMDH),ISYM0,WORK(KEND1),LWRK1)
C
C------------------------------------------
C     Calculate the F^U matrix
C------------------------------------------
C
      CALL GET_FOCKX(WORK(KFCKU),LABELRU,IDLSTRU,ISYMRU,WORK(KLAMDP),
     *                  ISYM0,WORK(KLAMDH),ISYM0,WORK(KEND1),LWRK1)
C
C-----------------------------
C     Read integrals.
C-----------------------------
C
C
C-----------------------------
C     Memory allocation.
C-----------------------------
C
C        isint1, isint2  - symmetry of integrals in standard H, transformed
C                  with LambdaH_0
C        isintr1 - symmetry of integrals in standard H, transformed
C                  with LambdaH_R1

      ISINT1    = 1
      ISINT2    = 1
      ISINT2X   = MULD2H(ISINT2,ISYMRX)
C
      KT3BOG1   = KEND1
      KT3BOG2   = KT3BOG1   + NTRAOC(ISINT2X)
      KT3BOL2   = KT3BOG2   + NTRAOC(ISINT2X)
      KT3OG1    = KT3BOL2   + NTRAOC(ISINT2X)
      KT3OG2    = KT3OG1    + NTRAOC(ISINT2)
      KEND1     = KT3OG2    + NTRAOC(ISINT2)
C
      KT3VIJG1  = KEND1
      KEND1     = KT3VIJG1  + NMAABCI(ISYM0)
      LWRK1     = LWORK     - KEND1
C
      KXGADCK   = KEND1 
      KXLADCK   = KXGADCK + NMAABCI(ISYMRX)
      KEND1     = KXLADCK + NMAABCI(ISYMRX)
      LWRK1     = LWORK     - KEND1
C
      KINTOC = KEND1
      KEND1  = KINTOC + NTOTOC(ISINT2)
      LWRK1     = LWORK     - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in CC3_AMATT3ZUOCC (3)')
      END IF
C
C-----------------------------------------------------------------
C     Construct occupied integrals which are required for occupied
C     part of contraction (KT3BOG2)
C-----------------------------------------------------------------
C
c     CALL INTOCC_T3BAR0(LUTOC,FNTOC,WORK(KLAMDH),ISYM0,WORK(KT3BOG1),
c    *                   WORK(KT3BOL1),WORK(KT3BOG2),WORK(KT3BOL2),
c    *                   WORK(KEND1),LWRK1)
      IOFF = 1
      IF (NTOTOC(ISINT2) .GT. 0) THEN
         CALL GETWA2(LUTOC,FNTOC,WORK(KINTOC),IOFF,NTOTOC(ISINT2))
      ENDIF
C
      !Transform (ia|j delta) integrals to (ia|j k-) and sort  G(j,i,k-,a)
      CALL CCLR_TROCC(WORK(KINTOC),WORK(KT3BOG1),WORK(KLAMDHX),ISYMRX,
     *               WORK(KEND1),LWRK1)
C
      !integrals g(ia|j k-) sorted as G(j,i,k-,a) resorted as G(k-,i,j,a)
      !KT3BOL2 are not used 
      CALL CC3_SRTOCC(WORK(KT3BOG1),WORK(KT3BOL2),ISINT2X)
C
      ! (ia|jk) sorted as T3BOG1(kija) and now sorted as T3BOG2(ajik)
      CALL CCFOP_SORT(WORK(KT3BOG1),WORK(KT3BOG2),ISINT2X,1)
C
C-----------------------------------------------------------------
C     Construct occupied integrals which are required to calculate    
C     t3_0 amplitudes
C-----------------------------------------------------------------
C
      CALL INTOCC_T30(LUCKJD,FNCKJD,WORK(KLAMDP),ISINT2,WORK(KT3OG1),
     *                WORK(KT3OG2),WORK(KEND1),LWRK1)
C
C----------------------------------------------
C     Get virtual integrals for t30 amplitudes
C     KT3VIJG1 : (ck|da) sorted as I(ad|ck)
C----------------------------------------------
C
      CALL INTVIR_T30_IJ(WORK(KT3VIJG1),ISYM0,WORK(KLAMDH),LUDELD,
     *                   FNDELD,WORK(KEND1),LWRK1)
C
C-------------------------------------------------------------------
C     Get virtual integrals for virtual part of contraction (KXGADCK)
C     KXGADCK g(kcad) = (kc ! ad) sorted as I(adck)
C     KXLADCK L(kcad) sorted as I(adck)
C-------------------------------------------------------------------
C
      CALL INTVIR_T3B0_JK(WORK(KXGADCK),WORK(KXLADCK),ISYM0,
     *                    WORK(KLAMDPX),ISYMRX,LU3VI,FN3VI,
     *                    LU3FOP,FN3FOP,
     *                    WORK(KEND1),LWRK1)
C
C----------------------------
C     Loop over K
C----------------------------
C
      ISYMTETAZ = MULD2H(ISYM0,ISYMRZ)
      ISYMTETAU = MULD2H(ISYM0,ISYMRU)
      ISYMTETAZU = MULD2H(ISYM0,ISYMZU)
      DO ISYMK = 1,NSYM

         DO K = 1,NRHF(ISYMK)
C
C----------------------------
C           Loop over L
C----------------------------
C
            DO ISYML = 1,NSYM
C
               ISYMKL = MULD2H(ISYMK,ISYML)
               ISYT30KL = MULD2H(ISYMKL,ISYMT3)
               ISTETAZKL  = MULD2H(ISYMKL,ISYMTETAZ)
               ISTETAUKL  = MULD2H(ISYMKL,ISYMTETAU)
               ISTETAZUKL  = MULD2H(ISYMKL,ISYMTETAZU)
C
               KT30KL = KEND1
               KTETAZKL  = KT30KL + NMAABCI(ISYT30KL)
               KTETAUKL   = KTETAZKL + NMAABCI(ISTETAZKL)
               KTETAZUKL   = KTETAUKL + NMAABCI(ISTETAUKL)
               KEND1   = KTETAZUKL + NMAABCI(ISTETAZUKL)
               LWRK1  = LWORK  - KEND1
C
               KTMAT = KEND1
               KEND1 = KTMAT + NMAABCI(ISTETAZUKL)
               LWRK1  = LWORK  - KEND1
C
               IF (LWRK1 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Memory available : ',LWORK
                  WRITE(LUPRI,*) 'Memory needed    : ',KEND1
                  CALL QUIT('Insufficient space in CC3_AMATT3ZUOCC (4)')
               END IF
C
               DO L = 1,NRHF(ISYML)
C
C
C-------------------------------------------
C                 Get T30^KL amplitudes
C-------------------------------------------
C
                  CALL DZERO(WORK(KT30KL),NMAABCI(ISYT30KL))
C
                  CALL GET_T30_IJ_O(WORK(KT30KL),ISYT30KL,WORK(KT2TP),
     *                              ISYM0,
     *                              WORK(KT3OG2),ISYM0,ISYML,L,ISYMK,K,
     *                              WORK(KEND1),LWRK1)
C
                  CALL GET_T30_IJ_V(WORK(KT30KL),ISYT30KL,WORK(KT2TP),
     *                              ISYM0,WORK(KT3VIJG1),
     *                              ISYM0,ISYML,L,ISYMK,K,
     *                              WORK(KEND1),LWRK1)
C
C--------------------------------------------------------------
C                 Divide by orbital energy difference and remove 
C                 forbidden elements
C--------------------------------------------------------------
C
                  CALL T3JK_DIA(WORK(KT30KL),ISYT30KL,ISYML,L,ISYMK,K,
     *                         WORK(KFOCKD))
                  CALL T3_FORBIDDEN_JK(WORK(KT30KL),ISYMT3,ISYML,L,
     *                                ISYMK,K)
C
c                 call sum_pt3_jk(work(kt30kl),isyml,l,isymk,k,isyt30kl,
c    *                           work(kx3am),3)
C
                  IF (IPRINT .GT. 55) THEN
                    WRITE(LUPRI,*)'ISYML,L,ISYMK,K ', ISYML,L,ISYMK,K
                    XNORMVAL = DDOT(NMAABCI(ISYT30KL),WORK(KT30KL),1,
     *                              WORK(KT30KL),1)
                    WRITE(LUPRI,*)'NORM OF KT30KL IN CC3_AMATT3ZUOCC ', 
     *                             XNORMVAL
                  END IF
C
C---------------------------------
C                 Calculate KTETAZKL
C---------------------------------
C
                  CALL DZERO(WORK(KTETAZKL),NMAABCI(ISTETAZKL))
C
                  CALL TETAX_JK_BC(WORK(KT30KL),ISYT30KL,WORK(KFCKZ),
     *                             ISYFCKZ,WORK(KTETAZKL),ISTETAZKL,
     *                             WORK(KEND1),LWRK1)
C
c                 call w3jk_dia(work(KTETAZKL),ISTETAZKL,isyml,l,isymk,
c    *                          k,work(kfockd),freqr1)
c                 call t3_forbidden_jk(work(KTETAZKL),ISYMTETAZ,isyml,l,
c    *                                  isymk,k)
c
c                 call sum_pt3_jk(work(KTETAZKL),isyml,l,isymk,k,
c    *                            ISTETAZKL,
c    *                            work(kx3am),1)
c
c                 xnormval = ddot(nmaabci(ISTETAZKL),work(KTETAZKL),1,
c    *                              work(KTETAZKL),1)
c                 write(lupri,*)'norm of KTETAZKL in CC3_AMATT3ZUOCC ',
c    *                             xnormval
C
                  CALL TETAX_JK_A(WORK(KT30KL),ISYT30KL,WORK(KFCKZ),
     *                             ISYFCKZ,WORK(KTETAZKL),ISTETAZKL,
     *                             WORK(KEND1),LWRK1)
C
                  CALL W3JK_DIA(WORK(KTETAZKL),ISTETAZKL,ISYML,L,ISYMK,
     *                           K,WORK(KFOCKD),FREQRZ)
                  CALL T3_FORBIDDEN_JK(WORK(KTETAZKL),ISYMTETAZ,ISYML,L,
     *                                  ISYMK,K)

c                 call sum_pt3_jk(work(KTETAZKL),isyml,l,isymk,k,
c    *                            ISTETAZKL,
c    *                            work(kx3am),4)

C
                  IF (IPRINT .GT. 55) THEN
                    XNORMVAL = DDOT(NMAABCI(ISTETAZKL),WORK(KTETAZKL),1,
     *                              WORK(KTETAZKL),1)
                    WRITE(LUPRI,*)'(FINAL) NORM OF KTETAZKL ', XNORMVAL
                  END IF
C
C---------------------------------
C                 Calculate KTETAUKL
C---------------------------------
C
                  CALL DZERO(WORK(KTETAUKL),NMAABCI(ISTETAUKL))
C
                  CALL TETAX_JK_BC(WORK(KT30KL),ISYT30KL,WORK(KFCKU),
     *                             ISYFCKU,WORK(KTETAUKL),ISTETAUKL,
     *                             WORK(KEND1),LWRK1)
C
C
                  CALL TETAX_JK_A(WORK(KT30KL),ISYT30KL,WORK(KFCKU),
     *                             ISYFCKU,WORK(KTETAUKL),ISTETAUKL,
     *                             WORK(KEND1),LWRK1)
C
                  CALL W3JK_DIA(WORK(KTETAUKL),ISTETAUKL,ISYML,L,ISYMK,
     *                           K,WORK(KFOCKD),FREQRU)
                  CALL T3_FORBIDDEN_JK(WORK(KTETAUKL),ISYMTETAU,ISYML,L,
     *                                  ISYMK,K)

c                 call sum_pt3_jk(work(KTETAUKL),isyml,l,isymk,k,
c    *                            ISTETAUKL,
c    *                            work(kx3am),1)

C
                  IF (IPRINT .GT. 55) THEN
                    XNORMVAL = DDOT(NMAABCI(ISTETAUKL),WORK(KTETAUKL),1,
     *                              WORK(KTETAUKL),1)
                    WRITE(LUPRI,*)'(FINAL) NORM OF KTETAUKL ', XNORMVAL
C
                  END IF
C
C-------------------------------------------------------------
C                 Calculate KTETAZUKL = <mu3|[U,KTETAZKL]|HF>
C-------------------------------------------------------------
C
                  CALL DZERO(WORK(KTETAZUKL),NMAABCI(ISTETAZUKL))
C
                  ! virtual contribution:
                  ! KTETAZUKL = sum_d U_ad KTETAZKL^{d-,b-,c-}_{iKL}
                  CALL TETAX_JK_A(WORK(KTETAZKL),ISTETAZKL,WORK(KFCKU),
     *                             ISYFCKU,WORK(KTETAZUKL),ISTETAZUKL,
     *                             WORK(KEND1),LWRK1)
C
                  ! occupied contribution:
                  ! KTETAZUKL = sum_m U_mi KTETAZKL^{a-,b-,c-}_{mKL}
                  CALL TETAX_JK_I(WORK(KTETAZKL),ISTETAZKL,WORK(KFCKU),
     *                             ISYFCKU,WORK(KTETAZUKL),ISTETAZUKL,
     *                             WORK(KEND1),LWRK1)
C
C-------------------------------------------------------------
C                 Calculate KTETAZUKL = <mu3|[Z,KTETAUKL]|HF>
C-------------------------------------------------------------
C
                  ! virtual contribution:
                  ! KTETAZUKL = sum_d Z_ad KTETAUKL^{d-,b-,c-}_{iKL}
                  CALL TETAX_JK_A(WORK(KTETAUKL),ISTETAUKL,WORK(KFCKZ),
     *                             ISYFCKZ,WORK(KTETAZUKL),ISTETAZUKL,
     *                             WORK(KEND1),LWRK1)
C
                  ! occupied contribution:
                  ! KTETAZUKL = sum_m Z_mi KTETAUKL^{a-,b-,c-}_{mKL}
                  CALL TETAX_JK_I(WORK(KTETAUKL),ISTETAUKL,WORK(KFCKZ),
     *                             ISYFCKZ,WORK(KTETAZUKL),ISTETAZUKL,
     *                             WORK(KEND1),LWRK1)
C
C---------------------------------------------------------------
C                 Divide by orbital energy difference and remove
C                 forbidden elements
C---------------------------------------------------------------
C
                  CALL T3_FORBIDDEN_JK(WORK(KTETAZUKL),ISYMTETAZU,ISYML,
     *                                 L,ISYMK,K)

c                 call sum_pt3_jk(work(KTETAZUKL),isyml,l,isymk,k,
c    *                            ISTETAZUKL,
c    *                            work(kx3am),4)

                  CALL W3JK_DIA(WORK(KTETAZUKL),ISTETAZUKL,ISYML,L,
     *                          ISYMK,K,WORK(KFOCKD),FREQZU)


C
                  IF (IPRINT .GT. 55) THEN
                    XNORMVAL = DDOT(NMAABCI(ISTETAZUKL),WORK(KTETAZUKL),
     *                              1,WORK(KTETAZUKL),1)
                    WRITE(LUPRI,*)'(FINAL) NORM OF KTETAZUKL ', XNORMVAL
C
                  END IF
C
C-------------------------------------------------------
C                 Project up to doubles space:
C                 omega2eff = - <mu2|[H,THETAZU^LK]|HF>
C-------------------------------------------------------
C
                  ! Fock matrix contribution
                  CALL CC3_CY2F_JK(OMEGA2EFF,ISYRES,WORK(KTETAZUKL),
     *                             ISTETAZUKL,WORK(KTMAT),
     *                             WORK(KFOCKX),ISYFCKX,
     *                             WORK(KEND1),LWRK1,
     *                             ISYML,L,ISYMK,K)
C
                  ! occupied integrals contribution
                  CALL CC3_CY2O_JK(OMEGA2EFF,ISYRES,WORK(KTETAZUKL),
     *                             ISTETAZUKL,WORK(KTMAT),WORK(KT3BOG2),
     *                             ISYM0,WORK(KEND1),LWRK1,
     *                             ISYML,L,ISYMK,K)
                  ! virtual integrals contribution
                  CALL CC3_CY2V_JK(OMEGA2EFF,ISYRES,WORK(KTETAZUKL),
     *                             ISTETAZUKL,WORK(KTMAT),WORK(KXGADCK),
     *                             ISYM0,WORK(KEND1),LWRK1,
     *                             ISYML,L,ISYMK,K)
C
               ENDDO   ! L
            ENDDO      ! ISYML
         ENDDO       ! K
      ENDDO          ! ISYMK 
C
c      write(lupri,*) 'T30KL in CC3_AMATT3ZUOCC'
c      call print_pt3(work(kx3am),isym0,4)
C
C
C------------------------------------------------------
C add: omega contributions to omegaeff     
C------------------------------------------------------
C
c     write(lupri,*)'OMEGA2EFF in CC3_AASDOCC'
c     call outpak(OMEGA2EFF,NT1AM(isyres),1,lupri)
 
      IF (IPRINT .GT. 55) THEN
         XNORMVAL = DDOT(NT2AM(ISYRES),OMEGA2EFF,1,OMEGA2EFF,1)
         WRITE(LUPRI,*)'Norm of OMEGA2EFF final ',XNORMVAL
      ENDIF
C 
C-------------
C     End
C-------------
C
      CALL QEXIT('AT3ZUOCC')
C
      RETURN
      END
C    /* Deck cc3_bmatt3zu */ 
      SUBROUTINE CC3_BMATT3ZU(LISTX,IDLSTX,LISTZU,IDLSTZU,
     *                        OMEGA2,OMEGA2EFF,
     *                        ISYRES,WORK,LWORK)
*
*******************************************************************
*
* This routine calculates the contribution to cubic response 
* (originating form both F matrix and B matrix):
*
* omega2eff = <mu2|[H^X,T3^ZU]|HF>
*
* which is contracted outside with the multipliers (zeroth-order for 
* F matrix or first-order for B matrix).
*
*
*
* T3^ZU is calculated using W^BD intermmediate.
*
* !!! ACTUALLY HERE WE CALCULATE ONLY THE PART OF THE TERM, WHICH ORIGINATES
*     FROM THE B-MATRIX CONTRIBUTIONS TO T3^ZU + 2 EXTRA CONTRIBUTIONS FROM
*     JACOBIAN (to get the rest of contribution call cc3_amatt3zu routine) !!!
*
* Thus here:
*
* W^BD = 1/2<mu3|[[[H,T1^Z],T1^U],T2^0]|HF> + <mu3|[H^Z,T2^U]|HF> !b-mat 
*      + <mu3|[[H,T1^ZU],T2^0]|HF> + <mu3|[H,T2^ZU]|HF>           !jacobian
*
* Note (from the expression for T3^ZU) that the permutation operator P(Z,U) 
* does NOT act on the last two terms.
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
C
      CHARACTER*6 FN3VI
      CHARACTER*8 FNTOC,FN3VI2
      CHARACTER*6 FNCKJD
      CHARACTER*4 FNDKBC
C        
      PARAMETER (FNTOC   = 'CCSDT_OC' )
      PARAMETER (FN3VI   = 'CC3_VI'  , FN3VI2  = 'CC3_VI12')
      PARAMETER (FNCKJD  = 'CKJDEL' )
      PARAMETER (FNDKBC  = 'DKBC' )

      CHARACTER*12 FN3SRTR, FNCKJDR, FNDELDR, FNDKBCR
      CHARACTER*13 FNCKJDR1, FNDKBCR1
      CHARACTER*13 FNCKJDR2, FNDKBCR2
      CHARACTER*13 FNCKJDR3, FNDKBCR3
C
      PARAMETER(FN3SRTR  = 'CCSDT_FBMAT1', FNCKJDR  = 'CCSDT_FBMAT2',
     *          FNDELDR  = 'CCSDT_FBMAT3', FNDKBCR  = 'CCSDT_FBMAT4')
C
      PARAMETER(FNCKJDR1  = 'CCSDT_FBMAT21',FNDKBCR1  = 'CCSDT_FBMAT41')
      PARAMETER(FNCKJDR2  = 'CCSDT_FBMAT22',FNDKBCR2  = 'CCSDT_FBMAT42')
      PARAMETER(FNCKJDR3  = 'CCSDT_FBMAT23',FNDKBCR3  = 'CCSDT_FBMAT43')
C
      INTEGER LUTOC,LU3VI,LU3VI2,LUCKJD,LUDKBC
      INTEGER LU3SRTR, LUCKJDR, LUDELDR, LUDKBCR
      INTEGER LUCKJDR1, LUDKBCR1
      INTEGER LUCKJDR2, LUDKBCR2
      INTEGER LUCKJDR3, LUDKBCR3
C
      CHARACTER*3 LISTZ, LISTU, LISTZU, LISTX
      CHARACTER*8 LABELZ, LABELU, LABELX, LABELZU
      INTEGER IDLSTZ,IDLSTU,IDLSTZU,IDLSTX
      INTEGER ISYMRZ,ISYMRU
C
      CHARACTER CDUMMY*1
C
      INTEGER ISYRES
      INTEGER LWORK
C
      INTEGER ISYM0,ISINT1,ISINT2,ISYMT2
      INTEGER KT2TP,KFOCKD,KLAMDP,KLAMDH,KXIAJB,KEND00,LWRK00
      INTEGER IOPTTCME,IOPT
      INTEGER KT2U,KT1Z,KEND1,LWRK1
      INTEGER ISINTR1Z,ISINTR2Z
      INTEGER ISINTZU
      INTEGER KRBJIA,KTROC,KTROC1,KTROC0Z,KTROC3,KINTOC
      INTEGER MAXOCC,KEND2,LWRK2
      INTEGER IOFF
      INTEGER ISYMD,ISCKB1,ISCKB2Z,ISCKBZU
      INTEGER KTRVI,KTRVI1,KTRVI0Z,KEND3,LWRK3,KTRVI7,KINTVI
      INTEGER ISYMB,ISYMZU,ISYMBD,ISWMAT,ISCKD2Z,ISCKDZU
      INTEGER KDIAGW,KWMAT,KTMAT,KEND4,KTRVI8Z,LWRK4,KTRVI9
      INTEGER LENSQW,KINDSQW
      INTEGER KEND0,LWRK0
      INTEGER ISINTR1U,ISINTR2U
      INTEGER KT2Z,KT1U,KTROC0U
      INTEGER ISCKB2U,KTRVI0U
      INTEGER ISCKD2U,KTRVI8U
      INTEGER MMAXOCC
c
      INTEGER ISYMRX,ISYMRZU,ISYFCKX,KT2ZU,KLAMDPX,KLAMDHX
      INTEGER KFOCKX,KT1ZU,KT1X,ISYMC,ISYMK,KOFF1,KOFF2,ISINTR1ZU
      INTEGER ISINTR2ZU,ISINT1X,KTROC0ZU,KT3OG1,MMMAXOCC,ISCKB1X
      INTEGER ISCKB2ZU,KTRVI0ZU,KT3VDG1,ISCKD2ZU,ISCKD1,KTRVI8ZU
      INTEGER KT3VBG1
c
      integer kx3am
c
C
#if defined (SYS_CRAY)
      REAL OMEGA2(*),OMEGA2EFF(*)
      REAL WORK(LWORK)
      REAL FREQZ,FREQU,FREQZU
      REAL XNORMVAL,DDOT
#else
      DOUBLE PRECISION OMEGA2(*),OMEGA2EFF(*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION FREQZ,FREQU,FREQZU
      DOUBLE PRECISION XNORMVAL,DDOT
#endif
C

      CALL QENTER('CC3_BMATT3ZU')
C
C-----------------------------
C     Lists handling and check
C-----------------------------
C
      ISYMRX = ISYLRT(IDLSTX)
      !FREQuency not needed for LISTX
      !LABELX = LRTLBL(IDLSTX) not needed for LISTX
C
C
      IF (LISTZU(1:3).EQ.'R2 ') THEN
        LISTZ  = 'R1 '
        LABELZ = LBLR2T(IDLSTZU,1)
        ISYMRZ  = ISYR2T(IDLSTZU,1)
        FREQZ  = FRQR2T(IDLSTZU,1)
        LORXRZ  = LORXR2T(IDLSTZU,1)
        IDLSTZ = IR1TAMP(LABELRZ,LORXRZ,FREQRZ,ISYMRZ)
C
        IF (LORXRZ) CALL QUIT('NO ORBITAL RELAX. IN CC3_BMATT3ZU')

        LISTU  = 'R1 '
        LABELU = LBLR2T(IDLSTZU,2)
        ISYMRU  = ISYR2T(IDLSTZU,2)
        FREQU  = FRQR2T(IDLSTZU,2)
        LORXRU  = LORXR2T(IDLSTZU,2)
        IDLSTU = IR1TAMP(LABELRU,LORXRU,FREQRU,ISYMRU)
C
        IF (LORXRU) CALL QUIT('NO ORBITAL RELAX. IN CC3_BMATT3ZU')
C
      ELSE
         CALL QUIT('Unknown list for LISTZU in CC3_BMATT3ZU.')
      END IF
C
      ISYMZU   = MULD2H(ISYMRZ,ISYMRU)
C
      IF (LISTZU(1:3).EQ.'R2 ') THEN 
         ISYMRZU = ISYLRT(IDLSTZU)
         !FREQuency for second-order amplitudes not needed - we get it
         ! below as the sum of "first-order" frequencies
         !LABELZU = LRTLBL(IDLSTZU) not needed for LISTZU
      ELSE
         CALL QUIT('Unknown list for LISTZU in CC3_BMATT3ZU')
      END IF
C
      !Symmetry check
      IF (ISYMZU .NE. ISYMRZU) THEN
         WRITE(LUPRI,*)'ISYMZU = ', ISYMZU
         WRITE(LUPRI,*)'ISYMRZU = ', ISYMRZU
         WRITE(LUPRI,*)'ISYMZU should be the same as ISYMRZU !!! '
         CALL QUIT('Symmetry mismatch in CC3_BMATT3ZU')
      END IF
C
C-------------------------------------------
C     Construct FREQZU = (omega_Z + omega_U)
C-------------------------------------------
C
      FREQZU = FREQZ + FREQU
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
C------------------------------------------
C     Open files (H0 integrals for WBD)
C------------------------------------------
C
      LUCKJD = -1
      LUDKBC = -1
C
      CALL WOPEN2(LUCKJD,FNCKJD,64,0)
      CALL WOPEN2(LUDKBC,FNDKBC,64,0)
C
C-----------------------------------------
C     Open temporary files (H^Z integrals)
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
C     Open temporary files (H^{Z,U} integrals) !double one-index transformed
C                                              !with first-order amplitudes
C------------------------------------------
C
      LUCKJDR1  = -1
      LUDKBCR1  = -1
C
      CALL WOPEN2(LUCKJDR1,FNCKJDR1,64,0)
      CALL WOPEN2(LUDKBCR1,FNDKBCR1,64,0)
C
C------------------------------------------
C     Open temporary files (H^U integrals)
C------------------------------------------
C
      LUCKJDR2  = -1
      LUDKBCR2  = -1
C
      CALL WOPEN2(LUCKJDR2,FNCKJDR2,64,0)
      CALL WOPEN2(LUDKBCR2,FNDKBCR2,64,0)
C
C------------------------------------------
C     Open temporary files (H^ZU integrals) !one-index transformed 
C                                           !with second-order amplitudes
C------------------------------------------
C
      LUCKJDR3  = -1
      LUDKBCR3  = -1
C
      CALL WOPEN2(LUCKJDR3,FNCKJDR3,64,0)
      CALL WOPEN2(LUDKBCR3,FNDKBCR3,64,0)
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
         WRITE(LUPRI,*)'Memory available: ',LWORK
         WRITE(LUPRI,*)'Memory needed: ',KEND00
         CALL QUIT('Out of memory in CC3_BMATT3ZU (zeroth allo.)')
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
            call quit('Insufficient space (T3) in CC3_BMATT3ZU (1)')
         END IF
      endif
C
C-------------------------------------------------
C     Read T1^Z and T2^Z
C     Read T1^U and T2^U
C     Read T1^ZU and T2^ZU   !second-order amplitudes
C
C     Calculate (ck|de)-Utilde  and (ck|lm)-Utilde
C     Calculate (ck|de)-{Z,U}tilde and (ck|lm)-{Z,U}tilde
C     Calculate (ck|de)-Ztilde  and (ck|lm)-Ztilde
C
C     Used to construct WBD intermmediate
C-------------------------------------------------
C
      ISYFCKX = MULD2H(ISYMOP,ISYMRX)
C
      KT2U   = KEND00
      KEND0  = KT2U   + NT2SQ(ISYMRU)
      LWRK0  = LWORK  - KEND0
C
      KT2Z   = KEND0
      KEND0  = KT2Z   + NT2SQ(ISYMRZ)
      LWRK0  = LWORK  - KEND0
C
      KT2ZU  = KEND0
      KEND0  = KT2ZU  + NT2SQ(ISYMRZU)
      LWRK0  = LWORK  - KEND0
C
      KLAMDPX = KEND0
      KLAMDHX = KLAMDPX + NLAMDT
      KEND0   = KLAMDHX + NLAMDT
      LWRK0   = LWORK   - KEND0
C
      KFOCKX  = KEND0
      KEND0   = KFOCKX  + N2BST(ISYFCKX)
      LWRK0   = LWORK  - KEND0
C    
      KT1U   = KEND0
      KEND1  = KT1U   + NT1AM(ISYMRU)
      LWRK1  = LWORK  - KEND1
C
      KT1Z   = KEND1
      KEND1  = KT1Z   + NT1AM(ISYMRZ)
      LWRK1  = LWORK  - KEND1
C
      KT1ZU  = KEND1
      KEND1  = KT1ZU  + NT1AM(ISYMRZU)
      LWRK1  = LWORK  - KEND1
C
      KT1X   = KEND1
      KEND1  = KT1X   + NT1AM(ISYMRX)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. NT2SQ(ISYMRZ)) THEN
         CALL QUIT('Out of memory in CC3_BMATT3ZU (TOT_T3Y) ')
      ENDIF
C
C---------------------------
C     Read Lambda_X matrices
C---------------------------
C
      CALL GET_LAMBDAX(WORK(KLAMDPX),WORK(KLAMDHX),LISTX,IDLSTX,
     *                 ISYMRX,WORK(KLAMDP),WORK(KLAMDH),WORK(KEND1),
     *                 LWRK1)
C
C--------------------------
C     Read in T1^Z and T2^Z
C--------------------------
C
      IOPT = 3
      CALL GET_T1_T2(IOPT,.TRUE.,WORK(KT1Z),WORK(KT2Z),LISTZ,IDLSTZ,
     *               ISYMRZ,WORK(KEND1),LWRK1)
C
      IF (IPRINT .GT. 55) THEN
         XNORMVAL = DDOT(NT1AM(ISYMRZ),WORK(KT1Z),1,WORK(KT1Z),1)
         WRITE(LUPRI,*) 'Norm of T1B  ',XNORMVAL
         XNORMVAL = DDOT(NT2SQ(ISYMRZ),WORK(KT2Z),1,WORK(KT2Z),1)
         WRITE(LUPRI,*) 'Norm of T2B  ',XNORMVAL
      ENDIF
C
C--------------------------
C     Read in T1^U and T2^U
C--------------------------
C
      IOPT = 3
      CALL GET_T1_T2(IOPT,.TRUE.,WORK(KT1U),WORK(KT2U),LISTU,IDLSTU,
     *               ISYMRU,WORK(KEND1),LWRK1)
C
      IF (IPRINT .GT. 55) THEN
         XNORMVAL = DDOT(NT1AM(ISYMRU),WORK(KT1U),1,WORK(KT1U),1)
         WRITE(LUPRI,*) 'Norm of T1C  ',XNORMVAL
         XNORMVAL = DDOT(NT2SQ(ISYMRU),WORK(KT2U),1,WORK(KT2U),1)
         WRITE(LUPRI,*) 'Norm of T2C  ',XNORMVAL
      ENDIF
C
C-------------------------------------------------------
C     Read in T1^ZU and T2^ZU   !second-order amplitudes
C-------------------------------------------------------
C
      IOPT = 3
      CALL GET_T1_T2(IOPT,.TRUE.,WORK(KT1ZU),WORK(KT2ZU),LISTZU,IDLSTZU,
     *               ISYMRZU,WORK(KEND1),LWRK1)
C
      IF (IPRINT .GT. 55) THEN
         XNORMVAL = DDOT(NT1AM(ISYMRZU),WORK(KT1ZU),1,WORK(KT1ZU),1)
         WRITE(LUPRI,*) 'Norm of T1BC  ',XNORMVAL
         XNORMVAL = DDOT(NT2SQ(ISYMRZU),WORK(KT2ZU),1,WORK(KT2ZU),1)
         WRITE(LUPRI,*) 'Norm of T2BC  ',XNORMVAL
      ENDIF
C
C------------------
C     Read in T1^X 
C------------------
C
      IOPT = 1
      CALL GET_T1_T2(IOPT,.FALSE.,WORK(KT1X),DUMMY,LISTX,IDLSTX,
     *               ISYMRX,WORK(KEND1),LWRK1)
C
      IF (IPRINT .GT. 55) THEN
         XNORMVAL = DDOT(NT1AM(ISYMRX),WORK(KT1X),1,WORK(KT1X),1)
         WRITE(LUPRI,*) 'Norm of T1A  ',XNORMVAL
      ENDIF
C
C------------------------------------------
C        Calculate the F^X matrix
C------------------------------------------
C
      CALL CC3LR_MFOCK(WORK(KEND1),WORK(KT1X),WORK(KXIAJB),ISYFCKX)
C
C     Put the F_{kc} part into correct F_{pq}
C
         CALL DZERO(WORK(KFOCKX),N2BST(ISYFCKX))
C
         DO ISYMC = 1, NSYM
            ISYMK = MULD2H(ISYFCKX,ISYMC)
            DO K = 1, NRHF(ISYMK)
               DO C = 1, NVIR(ISYMC)
                  KOFF1 = KFOCKX - 1
     *                  + IFCVIR(ISYMK,ISYMC)
     *                  + NORB(ISYMK)*(C - 1)
     *                  + K
                  KOFF2 = KEND1 - 1
     *                  + IT1AM(ISYMC,ISYMK)
     *                  + NVIR(ISYMC)*(K - 1)
     *                  + C
C
                  WORK(KOFF1) = WORK(KOFF2)
C
               ENDDO
            ENDDO
         ENDDO
C
         IF (IPRINT .GT. 55) THEN
            CALL AROUND( 'In CC3_BMATT3ZU: Fock^A MO matrix' )
            CALL CC_PRFCKMO(WORK(KFOCKX),ISYFCKX)
         ENDIF
C
C----------------------------------------------------
C     Integrals (ck|de)-tilde(Z) and (ck|lm)-tilde(Z)
C----------------------------------------------------
C
      ISINTR1Z = MULD2H(ISINT1,ISYMRZ)
      ISINTR2Z = MULD2H(ISINT2,ISYMRZ)
C
      CALL CC3_BARINT(WORK(KT1Z),ISYMRZ,WORK(KLAMDP),
     *                WORK(KLAMDH),WORK(KEND1),LWRK1,
     *                LU3SRTR,FN3SRTR,LUCKJDR,FNCKJDR)
C
      CALL CC3_SORT1(WORK(KEND1),LWRK1,2,ISINTR1Z,LU3SRTR,FN3SRTR,
     *               LUDELDR,FNDELDR,IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *               IDUMMY,CDUMMY)
C
      CALL CC3_SINT(WORK(KLAMDH),WORK(KEND1),LWRK1,ISINTR1Z,
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
C     Calculate the (ck|de)-{Z,U}tilde and (ck|lm)-{Z,U}tilde
C     (double one-index transformed with first-order amplitudes)
C------------------------------------------------------
C
      CALL CC3_3BARINT(ISYMRZ,LISTZ,IDLSTZ,ISYMRU,LISTU,IDLSTU,
     *                 IDUMMY,CDUMMY,IDUMMY,.FALSE.,
     *                 WORK(KLAMDP),WORK(KLAMDH),WORK(KEND1),LWRK1,
     *                 LU3SRTR,FN3SRTR,LUCKJDR1,FNCKJDR1)
C
      CALL CC3_SORT1(WORK(KEND1),LWRK1,2,ISYMZU,LU3SRTR,FN3SRTR,
     *               LUDELDR,FNDELDR,IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *               IDUMMY,CDUMMY)
C
      CALL CC3_SINT(WORK(KLAMDH),WORK(KEND1),LWRK1,ISYMZU,
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
C     Integrals (ck|de)-tilde(U) and (ck|lm)-tilde(U)
C----------------------------------------------------
C
      ISINTR1U = MULD2H(ISINT1,ISYMRU)
      ISINTR2U = MULD2H(ISINT2,ISYMRU)
C
      CALL CC3_BARINT(WORK(KT1U),ISYMRU,WORK(KLAMDP),
     *                WORK(KLAMDH),WORK(KEND1),LWRK1,
     *                LU3SRTR,FN3SRTR,LUCKJDR2,FNCKJDR2)
C
      CALL CC3_SORT1(WORK(KEND1),LWRK1,2,ISINTR1U,LU3SRTR,FN3SRTR,
     *               LUDELDR,FNDELDR,IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *               IDUMMY,CDUMMY)
C
      CALL CC3_SINT(WORK(KLAMDH),WORK(KEND1),LWRK1,ISINTR1U,
     *              LUDELDR,FNDELDR,LUDKBCR2,FNDKBCR2)
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
C----------------------------------------------------------
C     Integrals (ck|de)-tilde(ZU) and (ck|lm)-tilde(ZU)
C     (one-index transformed with second-order amplitudes)
C----------------------------------------------------------
C
      ISINTR1ZU = MULD2H(ISINT1,ISYMRZU)
      ISINTR2ZU = MULD2H(ISINT2,ISYMRZU)
C
      CALL CC3_BARINT(WORK(KT1ZU),ISYMRZU,WORK(KLAMDP),
     *                WORK(KLAMDH),WORK(KEND1),LWRK1,
     *                LU3SRTR,FN3SRTR,LUCKJDR3,FNCKJDR3)
C
      CALL CC3_SORT1(WORK(KEND1),LWRK1,2,ISINTR1ZU,LU3SRTR,FN3SRTR,
     *               LUDELDR,FNDELDR,IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *               IDUMMY,CDUMMY)
C
      CALL CC3_SINT(WORK(KLAMDH),WORK(KEND1),LWRK1,ISINTR1ZU,
     *              LUDELDR,FNDELDR,LUDKBCR3,FNDKBCR3)
C
C--------------------------------
C        Read occupied integrals
C-------------------------------
C
      ISINTZU = MULD2H(ISINT2,ISYMZU)
      ISINT1X = MULD2H(ISINT1,ISYMRX)
C
      !Use KEND0, because KT1Z, KT1U, KT1X, and KT1ZU are not needed any more
      KRBJIA = KEND0
      KTROC  = KRBJIA + NT2SQ(ISYRES)
      KTROC1 = KTROC  + NTRAOC(ISINT1X)   !KTROC - int. in contraction
      KEND1  = KTROC1 + NTRAOC(ISINT1X)   !KTROC1 - int. in contraction
      LWRK1  = LWORK  - KEND1
C
      KTROC0Z = KEND1                   
      KEND1   = KTROC0Z + NTRAOC(ISINTR2Z)!KTROC0Z - int. in <mu3|[H^Z,T2^U]|HF>
C
      KTROC0U = KEND1                   
      KEND1   = KTROC0U + NTRAOC(ISINTR2U)!KTROC0U - int. in <mu3|[H^U,T2^Z]|HF>C
      KTROC0ZU = KEND1                   
      KEND1   = KTROC0ZU + NTRAOC(ISINTR2ZU)!KTROC0ZU    <mu3|[H^ZU,T2^0]|HF>
C                                                              ====  
      KTROC3 = KEND1
      KEND1  = KTROC3   + NTRAOC(ISINTZU)!KTROC3- int. in <mu3|[H^{Z,U},T2]|HF>
C                                                               =======
      KT3OG1 = KEND1
      KEND1  = KT3OG1 + NTRAOC(ISINT1) ! <mu3|[H^0,T2^ZU]|HF>
C                                              ===
      KINTOC   = KEND1
      MAXOCC   = MAX(NTOTOC(ISINTR2Z),NTOTOC(ISINTZU))
      MMAXOCC  = MAX(NTOTOC(ISINTR2U),MAXOCC)
      MMMAXOCC = MAX(NTOTOC(ISINTR2ZU),MMAXOCC)
      KEND2    = KINTOC  + MAX(NTOTOC(ISINT1),MMMAXOCC)!KINTOC temporary storage
      LWRK2    = LWORK   - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in CC3_BMATT3ZU (2)')
      END IF
C
      CALL DZERO(WORK(KRBJIA),NT2SQ(ISYRES))
C
C-------------------------------------------------------------
C     occupied integrals for <mu3|[H^0,T2^ZU]|HF>
C-------------------------------------------------------------
C
      IOFF = 1
      IF (NTOTOC(ISINT1) .GT. 0) THEN
         CALL GETWA2(LUCKJD,FNCKJD,WORK(KINTOC),IOFF,NTOTOC(ISINT1))
      ENDIF

      !Transform (ai| delta j) integrals to (ai|kj) and sort as I(k,i,j,a)
      CALL CC3_TROCC(WORK(KINTOC),WORK(KT3OG1),WORK(KLAMDP),WORK(KEND2),
     *               LWRK2,ISINT1)
C
C-------------------------------------------------------------
C     Z-transformed occupied integrals for <mu3|[H^Z,T2^U]|HF>
C-------------------------------------------------------------
C
         IOFF = 1
         IF (NTOTOC(ISINTR2Z) .GT. 0) THEN
            CALL GETWA2(LUCKJDR,FNCKJDR,WORK(KINTOC),IOFF,
     *                  NTOTOC(ISINTR2Z))
         ENDIF
C
         IF (IPRINT .GT. 55) THEN
            XNORMVAL = DDOT(NTOTOC(ISINTR2Z),WORK(KINTOC),1,
     *                   WORK(KINTOC),1)
            WRITE(LUPRI,*) 'Norm of CCSDT_OC-INT (B transformed) ',
     *                      XNORMVAL
         ENDIF
C
         CALL CC3_TROCC(WORK(KINTOC),WORK(KTROC0Z),WORK(KLAMDP),
     *                  WORK(KEND2),LWRK2,ISINTR2Z)
C
C-------------------------------------------------------------
C     U-transformed occupied integrals for <mu3|[H^U,T2^Z]|HF>
C-------------------------------------------------------------
C
         IOFF = 1
         IF (NTOTOC(ISINTR2U) .GT. 0) THEN
            CALL GETWA2(LUCKJDR2,FNCKJDR2,WORK(KINTOC),IOFF,
     *                  NTOTOC(ISINTR2U))
         ENDIF
C
         IF (IPRINT .GT. 55) THEN
            XNORMVAL = DDOT(NTOTOC(ISINTR2U),WORK(KINTOC),1,
     *                   WORK(KINTOC),1)
            WRITE(LUPRI,*) 'Norm of CCSDT_OC-INT (C transformed) ',
     *                      XNORMVAL
         ENDIF
C
         CALL CC3_TROCC(WORK(KINTOC),WORK(KTROC0U),WORK(KLAMDP),
     *                  WORK(KEND2),LWRK2,ISINTR2U)
C
C--------------------------------------------------------------------------
C    Z,U-transformed occupied integrals for <mu3|[H^{Z,U},T2]|HF>
C    (double one-index transformed with first-order amplitudes)
C--------------------------------------------------------------------------
C
      IOFF = 1 
      IF (NTOTOC(ISINTZU) .GT. 0) THEN
         CALL GETWA2(LUCKJDR1,FNCKJDR1,WORK(KINTOC),IOFF,
     *               NTOTOC(ISINTZU))
      ENDIF
C
      CALL CC3_TROCC(WORK(KINTOC),WORK(KTROC3),WORK(KLAMDP),
     *                    WORK(KEND2),LWRK2,ISINTZU)
C
C---------------------------------------------------------------
C     ZU-transformed occupied integrals for <mu3|[H^ZU,T2^0]|HF>
C     (one-index transformed with second-order amplitudes)
C---------------------------------------------------------------
C
         IOFF = 1
         IF (NTOTOC(ISINTR2ZU) .GT. 0) THEN
            CALL GETWA2(LUCKJDR3,FNCKJDR3,WORK(KINTOC),IOFF,
     *                  NTOTOC(ISINTR2ZU))
         ENDIF
C
         IF (IPRINT .GT. 55) THEN
            XNORMVAL = DDOT(NTOTOC(ISINTR2ZU),WORK(KINTOC),1,
     *                   WORK(KINTOC),1)
            WRITE(LUPRI,*) 'Norm of CCSDT_OC-INT (BC transformed) ',
     *                      XNORMVAL
         ENDIF
C
         CALL CC3_TROCC(WORK(KINTOC),WORK(KTROC0ZU),WORK(KLAMDP),
     *                  WORK(KEND2),LWRK2,ISINTR2ZU)
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
         XNORMVAL  = DDOT(NTOTOC(ISINT1),WORK(KINTOC),1,WORK(KINTOC),1)
         WRITE(LUPRI,*) 'Norm of CCSDT_OC-INT ',XNORMVAL
      ENDIF
C
      !Transform (ia|j delta) integrals to (ia|j k)-tildeX and sort as (i,j,k,a)
      CALL CCLR_TROCC(WORK(KINTOC),WORK(KTROC),WORK(KLAMDHX),ISYMRX,
     *                  WORK(KEND2),LWRK2)
C
      !sort (i,j,k,a) as (a,i,j,k)
      CALL CCSDT_SRTOC2(WORK(KTROC),WORK(KTROC1),ISINT1X,
     *                  WORK(KEND2),LWRK2)
C
C---------------------
C     Start ISYMD-loop
C---------------------
C
      DO ISYMD = 1,NSYM
C
         ISCKB1  = MULD2H(ISINT1,ISYMD)      
         ISCKB1X  = MULD2H(ISINT1X,ISYMD)      
         ISCKB2Z = MULD2H(ISINTR2Z,ISYMD)
         ISCKB2U = MULD2H(ISINTR2U,ISYMD)
         ISCKB2ZU = MULD2H(ISINTR2ZU,ISYMD)
         ISCKBZU = MULD2H(ISINTZU,ISYMD)
C
C----------------------------------------
C        Read virtual integrals (fixed D)
C----------------------------------------
C
         !Use KEND1, because KINTOC is not needed any more
         KTRVI  = KEND1 
         KTRVI1 = KTRVI  + NCKATR(ISCKB1X)   !KTRVI - int. in contraction
         KEND2 = KTRVI1 + NCKATR(ISCKB1X)    !KTRVI1 - int. in contraction
         LWRK2  = LWORK  - KEND2
C
         KTRVI0Z  = KEND2
         KEND3 = KTRVI0Z + NCKATR(ISCKB2Z)!KTRVI0Z - int. in <mu3|[H^Z,T2^U]|HF>
         LWRK3    = LWORK    - KEND3      !                        ===
C
         KTRVI0U  = KEND3
         KEND3 = KTRVI0U + NCKATR(ISCKB2U)!KTRVI0U - int. in <mu3|[H^U,T2^Z]|HF>
         LWRK3    = LWORK    - KEND3      !                        ===
C
         KTRVI0ZU  = KEND3
         KEND3     = KTRVI0ZU + NCKATR(ISCKB2ZU)!KTRVI0ZU <mu3|[H^ZU,T2^0]|HF>
         LWRK3     = LWORK    - KEND3           !               ====   
C
         KTRVI7 = KEND3
         KEND3 = KTRVI7 + NCKATR(ISCKBZU)!KTRVI7 int. in <mu3|[H^{Z,U},T2^0]|HF>
         LWRK3  = LWORK  - KEND3         !                     ======
C
         KT3VDG1  = KEND3
         KEND3   = KT3VDG1  + NCKATR(ISCKB1)! <mu3|[H^0,T2^ZU]|HF>
C                                                   ===
         KINTVI = KEND3
         KEND4  = KINTVI + NCKA(ISCKB1)!KINTVI - temporary storage
         LWRK4  = LWORK  - KEND4
C
         IF (LWRK4 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND4
            CALL QUIT('Insufficient space in CC3_BMATT3ZU (3)')
         END IF
C
C--------------------
C        Start D-loop
C--------------------
C
         DO D = 1,NVIR(ISYMD)
C
C----------------------------------------------------------------------------
C           virtual integrals for <mu3|[H^0,T2^ZU]|HF> (fixed D)
C----------------------------------------------------------------------------
C
            IOFF = ICKBD(ISCKB1,ISYMD) + NCKATR(ISCKB1)*(D - 1) + 1
            IF (NCKATR(ISCKB1) .GT. 0) THEN
               CALL GETWA2(LUDKBC,FNDKBC,WORK(KT3VDG1),IOFF,
     *                     NCKATR(ISCKB1))
            ENDIF
C
C----------------------------------------------------------------------------
C           Z-transformed virtual integrals for <mu3|[H^Z,T2^U]|HF> (fixed D)
C----------------------------------------------------------------------------
C
            IOFF = ICKBD(ISCKB2Z,ISYMD) + NCKATR(ISCKB2Z)*(D - 1) + 1
            IF (NCKATR(ISCKB2Z) .GT. 0) THEN
               CALL GETWA2(LUDKBCR,FNDKBCR,WORK(KTRVI0Z),IOFF,
     &                     NCKATR(ISCKB2Z))
            ENDIF
C
C----------------------------------------------------------------------------
C           U-transformed virtual integrals for <mu3|[H^U,T2^Z]|HF> (fixed D)
C----------------------------------------------------------------------------
C
            IOFF = ICKBD(ISCKB2U,ISYMD) + NCKATR(ISCKB2U)*(D - 1) + 1
            IF (NCKATR(ISCKB2U) .GT. 0) THEN
               CALL GETWA2(LUDKBCR2,FNDKBCR2,WORK(KTRVI0U),IOFF,
     &                     NCKATR(ISCKB2U))
            ENDIF
C
C-----------------------------------------------------------------------------
C           Z,U-transformed virtual integrals for <mu3|[H^{Z,U},T2^0]|HF> 
C           (fixed D)
C           (doubly one-index transformed with first-order amplitudes)
C-----------------------------------------------------------------------------
C
            IOFF = ICKBD(ISCKBZU,ISYMD) + NCKATR(ISCKBZU)*(D - 1) + 1
            IF (NCKATR(ISCKBZU) .GT. 0) THEN
               CALL GETWA2(LUDKBCR1,FNDKBCR1,WORK(KTRVI7),IOFF,
     &                     NCKATR(ISCKBZU))
            ENDIF
C
C----------------------------------------------------------------------------
C           ZU-transformed virtual integrals for <mu3|[H^ZU,T2^0]|HF> (fixed D)
C           (one-index transformed with second-order amplitudes)
C----------------------------------------------------------------------------
C
            IOFF = ICKBD(ISCKB2ZU,ISYMD) + NCKATR(ISCKB2ZU)*(D - 1) + 1
            IF (NCKATR(ISCKB2ZU) .GT. 0) THEN
               CALL GETWA2(LUDKBCR3,FNDKBCR3,WORK(KTRVI0ZU),IOFF,
     &                     NCKATR(ISCKB2ZU))
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
            CALL CCLR_TRVIR(WORK(KINTVI),WORK(KTRVI),WORK(KLAMDPX),
     *                      ISYMRX,
     *                      ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
C
            IOFF = ICKAD(ISCKB1,ISYMD) + NCKA(ISCKB1)*(D - 1) + 1
            IF (NCKA(ISCKB1) .GT. 0) THEN
               CALL GETWA2(LU3VI,FN3VI,WORK(KINTVI),IOFF,
     *                     NCKA(ISCKB1))
            ENDIF
            CALL CCLR_TRVIR(WORK(KINTVI),WORK(KTRVI1),WORK(KLAMDPX),
     *                      ISYMRX,
     *                      ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
C
C---------------------------
C           Start ISYMB-loop
C---------------------------
C
            DO ISYMB = 1,NSYM
C
               ISYMZU  = MULD2H(ISYMRZ,ISYMRU)
               ISYMBD  = MULD2H(ISYMB,ISYMD)
               ISWMAT  = MULD2H(ISYMZU,ISYMBD)
               ISCKD2Z = MULD2H(ISINTR2Z,ISYMB)
               ISCKD2U = MULD2H(ISINTR2U,ISYMB)
               ISCKD2ZU = MULD2H(ISINTR2ZU,ISYMB)
               ISCKDZU = MULD2H(ISINTZU,ISYMB)
               ISCKD1  = MULD2H(ISINT1,ISYMB)
C
               !Use KEND3, because KINTVI is not needed any more
               KDIAGW  = KEND3
               KWMAT   = KDIAGW  + NCKIJ(ISWMAT)
               KINDSQW = KWMAT   + NCKIJ(ISWMAT)
               KTMAT   = KINDSQW + (6*NCKIJ(ISWMAT) - 1)/IRAT + 1
               KEND4   = KTMAT   + NCKIJ(ISWMAT)
C
               KTRVI8Z = KEND4
               KEND4   = KTRVI8Z + NCKATR(ISCKD2Z)!KTRVI8Z: <mu3|[H^Z,T2^U]|HF>
               LWRK4   = LWORK   - KEND4          !               ===
C
               KTRVI8U = KEND4
               KEND4   = KTRVI8U + NCKATR(ISCKD2U)!KTRVI8U: <mu3|[H^U,T2^Z]|HF>
               LWRK4   = LWORK   - KEND4          !               ===
C
               KTRVI8ZU = KEND4
               KEND4    = KTRVI8ZU + NCKATR(ISCKD2ZU)! <mu3|[H^ZU,T2^0]|HF>
               LWRK4    = LWORK    - KEND4           !       ====    
C
               KT3VBG1  = KEND4
               KEND4   = KT3VBG1  + NCKATR(ISCKD1)! <mu3|[H^0,T2^ZU]|HF>
C
               KTRVI9 = KEND4
               KEND4  = KTRVI9 + NCKATR(ISCKDZU)!KTRVI9: <mu3|[H^{Z,U},T2^0]|HF>
               LWRK4   = LWORK   - KEND4        !              =======
C
               IF (LWRK4 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Memory available : ',LWORK
                  WRITE(LUPRI,*) 'Memory needed    : ',KEND4
                  CALL QUIT('Insufficient space in CC3_BMATT3ZU (4)')
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
C                 virtual integrals for <mu3|[H^0,T2^ZU]|HF> (fixed B)
C----------------------------------------------------------------------------
C
                  IOFF = ICKBD(ISCKD1,ISYMB) + NCKATR(ISCKD1)*(B - 1)+ 1
                  IF (NCKATR(ISCKD1) .GT. 0) THEN
                     CALL GETWA2(LUDKBC,FNDKBC,WORK(KT3VBG1),IOFF,
     *                           NCKATR(ISCKD1))
                  ENDIF
C
C----------------------------------------------------------------------------
C                 Z-transformed virtual integrals for <mu3|[H^Z,T2^U]|HF> 
C                 (fixed B)
C----------------------------------------------------------------------------
C
                  IOFF = ICKBD(ISCKD2Z,ISYMB) + NCKATR(ISCKD2Z)*(B-1) +1
                  IF (NCKATR(ISCKD2Z) .GT. 0) THEN
                     CALL GETWA2(LUDKBCR,FNDKBCR,WORK(KTRVI8Z),IOFF,
     *                           NCKATR(ISCKD2Z))
                  ENDIF
C
C----------------------------------------------------------------------------
C                 U-transformed virtual integrals for <mu3|[H^U,T2^Z]|HF> 
C                 (fixed B)
C----------------------------------------------------------------------------
C
                  IOFF = ICKBD(ISCKD2U,ISYMB) + NCKATR(ISCKD2U)*(B-1) +1
                  IF (NCKATR(ISCKD2U) .GT. 0) THEN
                     CALL GETWA2(LUDKBCR2,FNDKBCR2,WORK(KTRVI8U),IOFF,
     *                           NCKATR(ISCKD2U))
                  ENDIF
C
C-----------------------------------------------------------------------------
C                 Z,U-transformed virtual integrals for <mu3|[H^{Z,U},T2^0]|HF> 
C                 (fixed B)
C                 (doubly one-index transformed with first-order amplitudes)
C-----------------------------------------------------------------------------
C
                  IOFF = ICKBD(ISCKDZU,ISYMB) + NCKATR(ISCKDZU)*(B-1) +1
                  IF (NCKATR(ISCKDZU) .GT. 0) THEN
                     CALL GETWA2(LUDKBCR1,FNDKBCR1,WORK(KTRVI9),IOFF,
     *                           NCKATR(ISCKDZU))
                  ENDIF
C
C----------------------------------------------------------------------------
C                 ZU-transformed virtual integrals for <mu3|[H^ZU,T2^0]|HF> 
C                 (fixed B)
C                 (one-index transformed with second-order amplitudes)
C----------------------------------------------------------------------------
C
                  IOFF = ICKBD(ISCKD2ZU,ISYMB)+ NCKATR(ISCKD2ZU)*(B-1)+1
                  IF (NCKATR(ISCKD2ZU) .GT. 0) THEN
                     CALL GETWA2(LUDKBCR3,FNDKBCR3,WORK(KTRVI8ZU),IOFF,
     *                           NCKATR(ISCKD2ZU))
                  ENDIF
C
C----------------------------------------------
C                 Calculate <mu3|[H^U,T2^Z]|HF> !b-matrix contribution
C----------------------------------------------
C
                  CALL WBD_GROUND(WORK(KT2Z),ISYMRZ,WORK(KTMAT),
     *                            WORK(KTRVI0U),WORK(KTRVI8U),
     *                            WORK(KTROC0U),ISINTR2U,WORK(KWMAT),
     *                            WORK(KEND4),LWRK4,
     *                            WORK(KINDSQW),LENSQW,
     *                            ISYMB,B,ISYMD,D)
C
C----------------------------------------------
C                 Calculate <mu3|[H^Z,T2^U]|HF> !b-matrix contribution
C----------------------------------------------
C
                  CALL WBD_GROUND(WORK(KT2U),ISYMRU,WORK(KTMAT),
     *                            WORK(KTRVI0Z),WORK(KTRVI8Z),
     *                            WORK(KTROC0Z),ISINTR2Z,WORK(KWMAT),
     *                            WORK(KEND4),LWRK4,
     *                            WORK(KINDSQW),LENSQW,
     *                            ISYMB,B,ISYMD,D)
C
C----------------------------------------------------------
C                 Calculate <mu3|[[[H,T1^Z],T1^U],T2^0]|HF> !b-matrix contrib.
C----------------------------------------------------------
C
                  CALL WBD_GROUND(WORK(KT2TP),ISYMT2,WORK(KTMAT),
     *                            WORK(KTRVI7),WORK(KTRVI9),
     *                            WORK(KTROC3),ISINTZU,WORK(KWMAT),
     *                            WORK(KEND4),LWRK4,
     *                            WORK(KINDSQW),LENSQW,
     *                            ISYMB,B,ISYMD,D)
C
C----------------------------------------------------
C                 Calculate <mu3|[[H,T1^ZU],T2^0]|HF> !jacobian contribution
C----------------------------------------------------
C
                  CALL WBD_GROUND(WORK(KT2TP),ISYMT2,WORK(KTMAT),
     *                            WORK(KTRVI0ZU),WORK(KTRVI8ZU),
     *                            WORK(KTROC0ZU),ISINTR2ZU,WORK(KWMAT),
     *                            WORK(KEND4),LWRK4,
     *                            WORK(KINDSQW),LENSQW,
     *                            ISYMB,B,ISYMD,D)
C
C-----------------------------------------------
C                 Calculate <mu3|[H^0,T2^ZU]|HF> !jacobian contribution
C-----------------------------------------------
C
                  CALL WBD_GROUND(WORK(KT2ZU),ISYMRZU,WORK(KTMAT),
     *                            WORK(KT3VDG1),WORK(KT3VBG1),
     *                            WORK(KT3OG1),ISINT1,WORK(KWMAT),
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
                  CALL T3_FORBIDDEN(WORK(KWMAT),ISYMZU,
     *                              ISYMB,B,ISYMD,D)

c                  call sum_pt3(work(kwmat),isymb,b,isymd,d,
c    *                          iswmat,work(kx3am),4)
c                  write(lupri,*) 'Total norm of WBD : ',
c    *            ddot(nckij(iswmat),work(kwmat),1,work(kwmat),1)

                  CALL WBD_DIA(B,ISYMB,D,ISYMD,FREQZU,
     *                         ISWMAT,WORK(KWMAT),
     *                         WORK(KDIAGW),WORK(KFOCKD))

C
C-----------------------------------------------
C                 Carry out the contractions...
C-----------------------------------------------
C
C
C
C------------------------------------------------------
C                 Calculate the  term <mu2|[H,W^BD(3)]|HF> 
C                 ( Fock matrix cont ) 
C                 added in OMEGA2EFF 
C------------------------------------------------------
C
                  CALL CC3_CY2F(OMEGA2EFF,ISYRES,WORK(KWMAT),ISWMAT,
     *                          WORK(KTMAT),WORK(KFOCKX),ISYFCKX,
     *                          WORK(KINDSQW),
     *                          LENSQW,WORK(KEND4),LWRK4,
     *                          ISYMB,B,ISYMD,D)
C                 
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
c      write(lupri,*) 'T3XY in CC3_BMATT3ZU  : '
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
c     write(lupri,*)'OMEGA2EFF after CC3_BMATT3ZU'
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
C-------------------------------------------
C     Close files (integrals in contraction)
C-------------------------------------------
C
      CALL WCLOSE2(LUTOC,FNTOC,'KEEP')
      CALL WCLOSE2(LU3VI,FN3VI,'KEEP')
      CALL WCLOSE2(LU3VI2,FN3VI2,'KEEP')
C
C------------------------------------------
C     Close files (H0 integrals for WBD)
C------------------------------------------
C
      CALL WCLOSE2(LUCKJD,FNCKJD,'KEEP')
      CALL WCLOSE2(LUDKBC,FNDKBC,'KEEP')
C
C------------------------------------------
C     Close temporary files (H^Z integrals)
C------------------------------------------
C
      CALL WCLOSE2(LU3SRTR,FN3SRTR,'DELETE')
      CALL WCLOSE2(LUCKJDR,FNCKJDR,'DELETE')
      CALL WCLOSE2(LUDELDR,FNDELDR,'DELETE')
      CALL WCLOSE2(LUDKBCR,FNDKBCR,'DELETE')
C
C------------------------------------------
C     Close temporary files (H^{Z,U} integrals)
C------------------------------------------
C
      CALL WCLOSE2(LUCKJDR1,FNCKJDR1,'DELETE')
      CALL WCLOSE2(LUDKBCR1,FNDKBCR1,'DELETE')
C
C------------------------------------------
C     Close temporary files (H^U integrals)
C------------------------------------------
C
      CALL WCLOSE2(LUCKJDR2,FNCKJDR2,'DELETE')
      CALL WCLOSE2(LUDKBCR2,FNDKBCR2,'DELETE')
C
C------------------------------------------
C     Close temporary files (H^ZU integrals)
C------------------------------------------
C
      CALL WCLOSE2(LUCKJDR3,FNCKJDR3,'DELETE')
      CALL WCLOSE2(LUDKBCR3,FNDKBCR3,'DELETE')
C
C-------------
C     End
C-------------
C
      CALL QEXIT('CC3_BMATT3ZU ')
C
      RETURN
      END
C
