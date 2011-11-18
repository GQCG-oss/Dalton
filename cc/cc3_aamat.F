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
      SUBROUTINE CC3_AAMATSD(LISTRX,IDLSTRX,
     *                                LISTRY,IDLSTRY,
     &                                OMEGA1,OMEGA2,
     &                                OMEGA1EFF,OMEGA2EFF,
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
      CHARACTER*3 LISTRX, LISTRY
      INTEGER IDLSTRX, IDLSTRY
C
      INTEGER ISYRES, LWORK
      integer isym
C
#if defined (SYS_CRAY)
      REAL OMEGA1(*),OMEGA2(*),OMEGA1EFF(*),OMEGA2EFF(*)
      REAL WORK(LWORK)
#else
      DOUBLE PRECISION OMEGA1(*),OMEGA2(*),OMEGA1EFF(*),OMEGA2EFF(*)
      DOUBLE PRECISION WORK(LWORK)
#endif
C
      CALL QENTER('AASD')
C
C----------------------
C     Lists check
C----------------------
C

      IF (LISTRX(1:3).EQ.'R1 ') THEN
        IF ((LISTRY(1:3).EQ.'R1 ') .OR. (LISTRY(1:3).EQ.'RE ')) THEN
          CONTINUE
        ELSE 
          WRITE(LUPRI,*)'LISTRX = ', LISTRX
          WRITE(LUPRI,*)'LISTRY = ', LISTRY
          WRITE(LUPRI,*)'Only implemented when '
     *    //'LISTRX = R1 and LISTRY = R1 or LISTRY = RE '
          CALL QUIT('Case not implemented in  CC3_AAMATSD')
        END IF
      ELSE
         WRITE(LUPRI,*)'LISTRX = ', LISTRX
         WRITE(LUPRI,*)'Only implemented for LISTRX = R1 '
         CALL QUIT('Case not implemented in  CC3_AAMATSD')
      END IF
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
      CALL CC3_AAMATSDVIR(LISTRX,IDLSTRX,
     *                                LISTRY,IDLSTRY,
     &                                OMEGA1,OMEGA2,
     &                                OMEGA1EFF,OMEGA2EFF,
     *                                ISYRES,WORK,LWORK,
     *                                LUCKJD,FNCKJD,LUTOC,FNTOC,
     *                                LUDELD,FNDELD,LUDKBC,FNDKBC,
     *                                LU3VI2,FN3VI2,LU3VI,FN3VI)
C
c
      IF (LISTRY(1:3).EQ.'R1 ') THEN
         CALL CC3_AAMATSDOCC(LISTRX,IDLSTRX,LISTRY,IDLSTRY,
     &                                   OMEGA1,OMEGA2,
     &                                   OMEGA1EFF,OMEGA2EFF,
     *                                   ISYRES,WORK,LWORK,
     *                                   LUCKJD,FNCKJD,LUTOC,FNTOC,
     *                                   LUDELD,FNDELD,LU3FOP,FN3FOP,
     *                                   LU3VI,FN3VI)
      END IF
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
      CALL QEXIT('AASD')
C
      RETURN
      END
C    /* Deck cc3_aamatsdvir */ 
      SUBROUTINE CC3_AAMATSDVIR(LISTRX,IDLSTRX,
     *                                LISTRY,IDLSTRY,
     &                                OMEGA1,OMEGA2,
     &                                OMEGA1EFF,OMEGA2EFF,
     *                                ISYRES,WORK,LWORK,
     *                                LUCKJD,FNCKJD,LUTOC,FNTOC,
     *                                LUDELD,FNDELD,LUDKBC,FNDKBC,
     *                                LU3VI2,FN3VI2,LU3VI,FN3VI)
C
*******************************************************************
*
* This routine calculates the virtual part of the right hand side vector T^XY
* for the paritioned second-order amplitude equations.
*
* Calculate the A^X matrix contributions (using W intermmediate): 
*
* W^BD =  <mu3|[[X,T1^Y],T3^0]|HF> + <mu3|[[X,T2^Y],T2^0]|HF>
*      +  <mu3|[X,T3^Y]|HF>
*
* (!!! To get the total contribution add the term from CC3_AAMATSDOCC !!!)

* projected up into the Singles-and-Doubles space:
*
* omega1eff = omega1eff - <mu1|[H,W^BD]|HF>
* omega2eff = omega2eff - <mu2|[H,W^BD]|HF>
*
*
*
C     Written by Poul Jorgensen and Filip Pawlowski, Aarhus, Spring 2003
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
      CHARACTER*12 FN3SRTR, FNCKJDRX, FNDELDRX, FNDKBCRX
      CHARACTER*13 FNCKJDRY,FNDELDRY,FNDKBCRY
      INTEGER LU3SRTR, LUCKJDRX, LUDELDRX, LUDKBCRX
      INTEGER LUCKJDRY,LUDELDRY,LUDKBCRY
C
      CHARACTER CDUMMY*1
C
      PARAMETER(FN3SRTR   = 'CCSDT_FBMAT1', FNCKJDRX  = 'CCSDT_FBMAT2', 
     *          FNDELDRX  = 'CCSDT_FBMAT3', FNDKBCRX  = 'CCSDT_FBMAT4')
C
      PARAMETER(FNCKJDRY  = 'CCSDT_FBMAT22',FNDELDRY  = 'CCSDT_FBMAT33',
     *          FNDKBCRY  = 'CCSDT_FBMAT44')
C
C
      CHARACTER*3 LISTRX, LISTRY
      CHARACTER*8 LABELRX,LABELRY
      INTEGER     IDLSTRX,IDLSTRY
      INTEGER     ISYMRX, ISYMRY
C
      LOGICAL T2XNET2Z
C
      INTEGER ISYRES, LWORK
      INTEGER ISYM0,ISINT1,ISINT2,ISYMT2,ISYMT3
      INTEGER KT2TP,KFOCKD,KLAMDP,KLAMDH,KEND00,LWRK00,KXIAJB,KFOCK0
      INTEGER IOPTTCME,IOPT,IOFF
      INTEGER ISYMWX,ISYMWY,ISYMWXY,ISYFCKX,ISYFCKY,ISYMXY
      INTEGER KT2X,KT2Y,KFCKX,KFCKY,KEND0,LWRK0
      INTEGER KFCKXYV,KFCKXYO,KFCKYXV,KFCKYXO,KT1X,KT1Y,KEND1,LWRK1
      INTEGER KLAMDPY,KLAMDHY,KLAMDPX,KLAMDHX
      INTEGER ISINTRX1,ISINTRX2,ISINTRY1,ISINTRY2
      INTEGER KRBJIAXY,KT3OG1,KT3OG2,KW3XOGX1,KW3YOGY1,KTROC,KTROC1
      INTEGER KINTOC
      INTEGER ISYMD,ISCKB1,ISCKB2,ISCKB2X,ISCKB2Y
      INTEGER KT3VDG3,KT3VDG2,KEND2,LWRK2
      INTEGER KT3VDG1,KEND3,LWRK3,KW3XVDGX1,KW3YVDGY1,KT3VDGX3,KT3VDGY3
      INTEGER KTRVI,KTRVI1,MAXXY
      INTEGER ISYMB,ISYALJ,ISYALJ2,ISYMBD,ISCKIJ,ISCKD2,ISCKD2X,ISCKD2Y
      INTEGER ISWMATX,ISWMATY,ISWMATXY
      INTEGER KSMAT,KSMAT3,KUMAT,KUMAT3,KDIAG,KT3MAT
      INTEGER KDIAGWX,KWMATX,KINDSQWX,KDIAGWY,KWMATY,KINDSQWY
      INTEGER KINDSQ,KINDEX,KINDEX2
      INTEGER MAXXXYY,KTMATW,KT3VBG1,KT3VBG2,KT3VBG3,KEND4,LWRK4
      INTEGER KW3XVDGX2,KW3YVDGY2,KT3VBGX3,KT3VBGY3,KINTVI
      INTEGER KINDSQWXY,KDIAGWXY,KWMATXY
      INTEGER LENSQ,LENSQWX,LENSQWY,LENSQWXY
      INTEGER AIBJCK_PERM
      INTEGER ILSTSYM
c
      integer kx3am
c
#if defined (SYS_CRAY)
      REAL OMEGA1(*),OMEGA2(*),OMEGA1EFF(*),OMEGA2EFF(*)
      REAL WORK(LWORK)
      REAL FREQRX,FREQRY,FREQXY
      REAL XNORMVAL,DDOT
#else
      DOUBLE PRECISION OMEGA1(*),OMEGA2(*),OMEGA1EFF(*),OMEGA2EFF(*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION FREQRX,FREQRY,FREQXY
      DOUBLE PRECISION XNORMVAL,DDOT
#endif
C
      CALL QENTER('AMATSDV')
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
COMMENT COMMENT COMMENT 
c
c    summing up contributions in omega1eff
c    must be changed when introduced in real implementation
c
c       CALL DZERO(OMEGA1EFF,NT1AM(1))
c
COMMENT COMMENT COMMENT 
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
      LUCKJDRX  = -1
      LUCKJDRY  = -1
      LUDELDRX  = -1
      LUDELDRY  = -1
      LUDKBCRX  = -1
      LUDKBCRY  = -1
C
      CALL WOPEN2(LU3SRTR,FN3SRTR,64,0)
      CALL WOPEN2(LUCKJDRX,FNCKJDRX,64,0)
      CALL WOPEN2(LUCKJDRY,FNCKJDRY,64,0)
      CALL WOPEN2(LUDELDRX,FNDELDRX,64,0)
      CALL WOPEN2(LUDELDRY,FNDELDRY,64,0)
      CALL WOPEN2(LUDKBCRX,FNDKBCRX,64,0)
      CALL WOPEN2(LUDKBCRY,FNDKBCRY,64,0)
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
      KFOCK0 = KXIAJB + NT2AM(ISINT1)
      KEND00 = KFOCK0 + N2BST(ISYM0)
      LWRK00 = LWORK  - KEND00
C
      IF (LWRK00 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_AAMATSDVIR (zeroth allo.')
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

*     write(lupri,*) 'XIAJB(x1) in cc3_aamat listry:i ',listry
*     call outpak(WORK(KXIAJB),NT1AM(1),1,lupri)
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
C------------------------------------------------------------------
C     Read in Fock matrix in AO basis (from file) and transform to    
C     Lambda_0 basis.             
C------------------------------------------------------------------
C
      CALL GET_FOCK0(WORK(KFOCK0),WORK(KLAMDP),WORK(KLAMDH),
     *               WORK(KEND00),LWRK00)
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
            call quit('Insufficient space (T3) in CC3_AAMATSDVIR')
         END IF
      endif
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
C
C determining R1 labels
C

      IF (LISTRX(1:3).EQ.'R1 ') THEN
            ISYMRX = ISYLRT(IDLSTRX)
            FREQRX  = FRQLRT(IDLSTRX)
            LABELRX = LRTLBL(IDLSTRX)
C
         IF (LISTRY(1:3).EQ.'R1 ') THEN
            ISYMRY = ISYLRT(IDLSTRY)
            FREQRY  = FRQLRT(IDLSTRY)
            LABELRY = LRTLBL(IDLSTRY)
         ELSE IF (LISTRY(1:3).EQ.'RE ') THEN
            ISYMRY = ILSTSYM(LISTRY,IDLSTRY)
            FREQRY  = EIGVAL(IDLSTRY)
            LABELRY = '- none -'
         ELSE
           WRITE(LUPRI,*) 'LISTRY = ', LISTRY
           CALL QUIT('Unnkown list for LISTRY in CC3_AAMATSDVIR.')
         END IF
C
      ELSE
        WRITE(LUPRI,*) 'LISTRX = ', LISTRX
        CALL QUIT('Unknown list for LISTRX in CC3_AAMATSDVIR.')
      END IF
C
      FREQXY = FREQRX + FREQRY
C
      ISYMWX   = MULD2H(ISYMT3,ISYMRX)
      ISYMWY   = MULD2H(ISYMT3,ISYMRY)
C
      ISYMWXY   = MULD2H(ISYMWX,ISYMRY)
C
      IF (ISYRES.NE.ISYMWXY) THEN
         WRITE(LUPRI,*) 'INCONSISTENT SYMMETRY SPECIFICATION'
         WRITE(LUPRI,*) 
     *      'ISYRES =',ISYRES,'ISYMWXY =',ISYMWXY
         STOP' CC3_AAMATSDVIR inconsistent symmetry specification' 
      END IF
C
      ISYFCKX = ISYMRX
      ISYFCKY = ISYMRY
      ISYMXY = MULD2H(ISYMRX,ISYMRY)
C
C-------------------------------------------------
C     Read T1^B and T2^B
C     Calculate (ck|de)-tilde(B) and (ck|lm)-tilde(B)
C     Used to construct T3^B
C-------------------------------------------------
C
      KT2X   = KEND00
      KT2Y   = KT2X   + NT2SQ(ISYMRX)
      KFCKX  = KT2Y   + NT2SQ(ISYMRY)
      KFCKY  = KFCKX  + N2BST(ISYFCKX)
      KEND0  = KFCKY  + N2BST(ISYFCKY)
      LWRK0  = LWORK  - KEND0
C
      KFCKXYV = KEND0
      KFCKXYO = KFCKXYV + N2BST(ISYMXY)
      KFCKYXV = KFCKXYO + N2BST(ISYMXY)
      KFCKYXO = KFCKYXV + N2BST(ISYMXY)
      KEND0   = KFCKYXO + N2BST(ISYMXY)
C
      KT1X   = KEND0
      KT1Y   = KT1X   + NT1AM(ISYMRX)
      KEND1  = KT1Y   + NT1AM(ISYMRY)
C
      KLAMDPY = KEND1
      KLAMDHY = KLAMDPY + NLAMDT
      KLAMDPX = KLAMDHY + NLAMDT
      KLAMDHX = KLAMDPX + NLAMDT
      KEND1   = KLAMDHX + NLAMDT
      LWRK1   = LWORK   - KEND1
C
      IF (LWRK1 .LT. NT2SQ(ISYMRX)) THEN
         CALL QUIT('Out of memory in CC3_AAMATSDVIR (TOT_T3Y) ')
      ENDIF
C
C     T1^X and T2^X amplitudes
C     
      IF (LISTRY(1:3).EQ.'R1 ') THEN
        IOPT = 3
        CALL GET_T1_T2(IOPT,.TRUE.,WORK(KT1X),WORK(KT2X),LISTRX,IDLSTRX,
     *                 ISYMRX,WORK(KEND1),LWRK1)
C
        IF (IPRINT .GT. 55) THEN
           XNORMVAL = DDOT(NT2SQ(ISYMRX),WORK(KT2X),1,WORK(KT2X),1)
           WRITE(LUPRI,*) 'Norm of T2X  ',XNORMVAL
           XNORMVAL = DDOT(NT1AM(ISYMRX),WORK(KT1X),1,WORK(KT1X),1)
           WRITE(LUPRI,*) 'Norm of T1X  ',XNORMVAL
        ENDIF
      END IF
C
C     T1^Y and T2^Y amplitudes
C     
      IOPT = 3
      CALL GET_T1_T2(IOPT,.TRUE.,WORK(KT1Y),WORK(KT2Y),LISTRY,IDLSTRY,
     *               ISYMRY,WORK(KEND1),LWRK1)
C
      IF (IPRINT .GT. 55) THEN
         XNORMVAL = DDOT(NT2SQ(ISYMRY),WORK(KT2Y),1,WORK(KT2Y),1)
         WRITE(LUPRI,*) 'Norm of T2Y  ',XNORMVAL
         XNORMVAL = DDOT(NT1AM(ISYMRY),WORK(KT1Y),1,WORK(KT1Y),1)
         WRITE(LUPRI,*) 'Norm of T1Y  ',XNORMVAL
      ENDIF
C
C     Integrals H^X in T3^X
C

      ISINTRX1 = MULD2H(ISINT1,ISYMRX)
      ISINTRX2 = MULD2H(ISINT2,ISYMRX)
C
      IF (LISTRY(1:3).EQ.'R1 ') THEN
        CALL CC3_BARINT(WORK(KT1X),ISYMRX,WORK(KLAMDP),
     *                  WORK(KLAMDH),WORK(KEND1),LWRK1,
     *                  LU3SRTR,FN3SRTR,LUCKJDRX,FNCKJDRX)
C
        CALL CC3_SORT1(WORK(KEND1),LWRK1,2,ISINTRX1,LU3SRTR,FN3SRTR,
     *                 LUDELDRX,FNDELDRX,IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *                 IDUMMY,CDUMMY)
C
        CALL CC3_SINT(WORK(KLAMDH),WORK(KEND1),LWRK1,ISINTRX1,
     *                LUDELDRX,FNDELDRX,LUDKBCRX,FNDKBCRX)
      END IF
C
C     Integrals H^Y in T3^Y
C

      ISINTRY1 = MULD2H(ISINT1,ISYMRY)
      ISINTRY2 = MULD2H(ISINT2,ISYMRY)
C
      CALL CC3_BARINT(WORK(KT1Y),ISYMRY,WORK(KLAMDP),
     *                WORK(KLAMDH),WORK(KEND1),LWRK1,
     *                LU3SRTR,FN3SRTR,LUCKJDRY,FNCKJDRY)
C
      CALL CC3_SORT1(WORK(KEND1),LWRK1,2,ISINTRY1,LU3SRTR,FN3SRTR,
     *               LUDELDRY,FNDELDRY,IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *               IDUMMY,CDUMMY)
C
      CALL CC3_SINT(WORK(KLAMDH),WORK(KEND1),LWRK1,ISINTRY1,
     *              LUDELDRY,FNDELDRY,LUDKBCRY,FNDKBCRY)
C
C------------------------------------------
C     Calculate the F^X matrix
C------------------------------------------
C
      CALL GET_FOCKX(WORK(KFCKX),LABELRX,IDLSTRX,ISYMRX,WORK(KLAMDP),
     *                  ISYM0,WORK(KLAMDH),ISYM0,WORK(KEND1),LWRK1)

C
C------------------------------------------
C     Calculate the F^Y matrix
C------------------------------------------
C
      IF (LISTRY(1:3).EQ.'R1 ') THEN
        CALL GET_FOCKX(WORK(KFCKY),LABELRY,IDLSTRY,ISYMRY,WORK(KLAMDP),
     *                    ISYM0,WORK(KLAMDH),ISYM0,WORK(KEND1),LWRK1)
      END IF
C
C------------------------------------------
C     Calculate the [X,T1^Y] matrix
C     Recall that we only need the occ-occ and vir-vir block.
C------------------------------------------
C
      CALL GET_LAMBDAX(WORK(KLAMDPY),WORK(KLAMDHY),LISTRY,IDLSTRY,
     *                 ISYMRY,WORK(KLAMDP),WORK(KLAMDH),WORK(KEND1),
     *                 LWRK1)
      ! get vir-vir block X_(c-,d)
      CALL GET_FOCKX(WORK(KFCKXYV),LABELRX,IDLSTRX,ISYMRX,WORK(KLAMDPY),
     *                  ISYMRY,WORK(KLAMDH),ISYM0,WORK(KEND1),LWRK1)
      ! get occ-occ block X_(l,k-)
      CALL GET_FOCKX(WORK(KFCKXYO),LABELRX,IDLSTRX,ISYMRX,WORK(KLAMDP),
     *                  ISYM0,WORK(KLAMDHY),ISYMRY,WORK(KEND1),LWRK1)
C
C------------------------------------------
C     Calculate the [Y,T1^X] matrix
C     Recall that we only need the occ-occ and vir-vir block.
C------------------------------------------
C
      IF (LISTRY(1:3).EQ.'R1 ') THEN
        CALL GET_LAMBDAX(WORK(KLAMDPX),WORK(KLAMDHX),LISTRX,IDLSTRX,
     *                   ISYMRX,WORK(KLAMDP),WORK(KLAMDH),WORK(KEND1),
     *                   LWRK1)
        ! get vir-vir block Y_(c-,d)
        CALL GET_FOCKX(WORK(KFCKYXV),LABELRY,IDLSTRY,ISYMRY,
     *                 WORK(KLAMDPX),
     *                    ISYMRX,WORK(KLAMDH),ISYM0,WORK(KEND1),LWRK1)
        ! get occ-occ block Y_(l,k-)
        CALL GET_FOCKX(WORK(KFCKYXO),LABELRY,IDLSTRY,ISYMRY,
     *                 WORK(KLAMDP),
     *                    ISYM0,WORK(KLAMDHX),ISYMRX,WORK(KEND1),LWRK1)
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
      KRBJIAXY = KEND0
      KT3OG1   = KRBJIAXY + NT2SQ(ISYRES)
      KT3OG2   = KT3OG1 + NTRAOC(ISINT2)
      KEND1    = KT3OG2+ NTRAOC(ISINT2)
      LWRK1    = LWORK  - KEND1
C
      KW3XOGX1 = KEND1
      KEND1   = KW3XOGX1 + NTRAOC(ISINTRX2)
C
      KW3YOGY1 = KEND1
      KEND1   = KW3YOGY1 + NTRAOC(ISINTRY2)
      LWRK1  = LWORK  - KEND1
C
      KTROC  = KEND1
      KTROC1 = KTROC  + NTRAOC(ISINT1)   !KTROC - int. in projection
      KEND1  = KTROC1 + NTRAOC(ISINT1)   !KTROC1 - int. in projection
      LWRK1  = LWORK  - KEND1
C
      KINTOC = KEND1
      KEND1  = KINTOC  + NTOTOC(ISINT1) !KINTOC - temporary storage
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in CC3_AAMATSDVIR')
      END IF
C
      CALL DZERO(WORK(KRBJIAXY),NT2SQ(ISYRES))
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
C     X transformed Occupied integrals.
C-----------------------------------------
C
      IF (LISTRY(1:3).EQ.'R1 ') THEN
        CALL INTOCC_T3X(LUCKJDRX,FNCKJDRX,WORK(KLAMDP),ISINTRX2,
     *                  WORK(KW3XOGX1),WORK(KEND1),LWRK1)
      END IF
C
C------------------------------------------
C     Y transformed Occupied integrals.
C-----------------------------------------
C
      CALL INTOCC_T3X(LUCKJDRY,FNCKJDRY,WORK(KLAMDP),ISINTRY2,
     *                WORK(KW3YOGY1),WORK(KEND1),LWRK1)
C
C-----------------------------------
C   Occupied integrals in projection
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
     *                  WORK(KEND1),LWRK1)
C
      !sort (i,j,k,a) as (a,i,j,k)
      CALL CCSDT_SRTOC2(WORK(KTROC),WORK(KTROC1),ISINT1,
     *                  WORK(KEND1),LWRK1)
c
C
C----------------------------
C     General loop structure.
C----------------------------
C
      DO ISYMD = 1,NSYM

         ISCKB1  = MULD2H(ISINT1,ISYMD)
         ISCKB2  = MULD2H(ISINT2,ISYMD)
         ISCKB2X = MULD2H(ISINTRX2,ISYMD)
         ISCKB2Y = MULD2H(ISINTRY2,ISYMD)
C
         IF (IPRINT .GT. 55) THEN
C
            WRITE(LUPRI,*) 'In CC3_AAMATSDVIR: ISCKB2 :',ISCKB2
            WRITE(LUPRI,*) 'In CC3_AAMATSDVIR: ISCKB2X:',ISCKB2X
            WRITE(LUPRI,*) 'In CC3_AAMATSDVIR: ISCKB2Y:',ISCKB2Y

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
         KW3XVDGX1  = KEND3
         KEND3    = KW3XVDGX1  + NCKATR(ISCKB2X)
         LWRK3    = LWORK    - KEND3
C
         KW3YVDGY1  = KEND3
         KEND3    = KW3YVDGY1  + NCKATR(ISCKB2Y)
         LWRK3    = LWORK    - KEND3
C
         KT3VDGX3 = KEND3
         KEND3    = KT3VDGX3 + NCKATR(ISCKB2X)
C
         KT3VDGY3 = KEND3
         KEND3    = KT3VDGY3 + NCKATR(ISCKB2Y)
C
         KTRVI  = KEND3
         KTRVI1 = KTRVI  + NCKATR(ISCKB1)   !KTRVI - int. in projection
         KEND3 = KTRVI1 + NCKATR(ISCKB1)    !KTRVI1 - int. in projection
         LWRK3  = LWORK  - KEND3
C
         KINTVI  = KEND3
         MAXXY   = MAX(NCKA(ISCKB2X),NCKA(ISCKB2Y))
         KEND3   = KINTVI  + MAX(MAXXY,NCKA(ISCKB1))
         LWRK3   = LWORK   - KEND3
C
         IF (LWRK3 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND3
            CALL QUIT('Insufficient space in CC3_AAMATSDVIR')
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
C        Read X transformed virtual integrals (H^X) used for W^X 
C-----------------------------------------------------------------------
C
            IF (LISTRY(1:3).EQ.'R1 ') THEN 
              IOFF = ICKBD(ISCKB2X,ISYMD) + NCKATR(ISCKB2X)*(D - 1) + 1
              IF (NCKATR(ISCKB2X) .GT. 0) THEN
                 CALL GETWA2(LUDKBCRX,FNDKBCRX,WORK(KW3XVDGX1),IOFF,
     &                       NCKATR(ISCKB2X))
              ENDIF
            END IF
C
C-----------------------------------------------------------------------
C        Read Y transformed virtual integrals (H^Y) used for W^Y 
C-----------------------------------------------------------------------
C
            IOFF = ICKBD(ISCKB2Y,ISYMD) + NCKATR(ISCKB2Y)*(D - 1) + 1
            IF (NCKATR(ISCKB2Y) .GT. 0) THEN
               CALL GETWA2(LUDKBCRY,FNDKBCRY,WORK(KW3YVDGY1),IOFF,
     &                     NCKATR(ISCKB2Y))
            ENDIF
C
C--------------------------------------------------------------------
C           Read virtual integrals [H,T1X] where X is LISTRX (used in W^X)
C--------------------------------------------------------------------
C
            IF (LISTRY(1:3).EQ.'R1 ') THEN
              IF (NCKA(ISCKB2X) .GT. 0) THEN
                 IOFF = ICKAD(ISCKB2X,ISYMD) +
     &                  NCKA(ISCKB2X)*(D - 1) + 1
                 CALL GETWA2(LUDELDRX,FNDELDRX,WORK(KINTVI),IOFF,
     *                NCKA(ISCKB2X))
              ENDIF
C
              CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KT3VDGX3),
     *                         WORK(KLAMDH),ISYMD,D,ISINTRX2,
     *                         WORK(KEND3),LWRK3)
            END IF
C
C--------------------------------------------------------------------
C           Read virtual integrals [H,T1Y] where Y is LISTRY (used in W^Y)
C--------------------------------------------------------------------
C
            IF (NCKA(ISCKB2Y) .GT. 0) THEN
               IOFF = ICKAD(ISCKB2Y,ISYMD) +
     &                NCKA(ISCKB2Y)*(D - 1) + 1
               CALL GETWA2(LUDELDRY,FNDELDRY,WORK(KINTVI),IOFF,
     *              NCKA(ISCKB2Y))
            ENDIF
C
            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KT3VDGY3),
     *                       WORK(KLAMDH),ISYMD,D,ISINTRY2,
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
            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI),WORK(KLAMDP),
     *                       ISYMD,D,ISYMOP,WORK(KEND3),LWRK3)
C
            IOFF = ICKAD(ISCKB1,ISYMD) + NCKA(ISCKB1)*(D - 1) + 1
            IF (NCKA(ISCKB1) .GT. 0) THEN
               CALL GETWA2(LU3VI,FN3VI,WORK(KINTVI),IOFF,
     *                     NCKA(ISCKB1))
            ENDIF
            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI1),WORK(KLAMDP),
     *                          ISYMD,D,ISYMOP,WORK(KEND3),LWRK3)
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
               ISCKD2X = MULD2H(ISINTRX2,ISYMB)
               ISCKD2Y = MULD2H(ISINTRY2,ISYMB)
               ISWMATX  = MULD2H(ISYMRX,ISCKIJ)
               ISWMATY  = MULD2H(ISYMRY,ISCKIJ)
               ISWMATXY  = MULD2H(ISWMATX,ISYMRY)

               IF (IPRINT .GT. 55) THEN

                  WRITE(LUPRI,*) 'In CC3_AAMATSDVIR: ISYMD  :',ISYMD
                  WRITE(LUPRI,*) 'In CC3_AAMATSDVIR: ISYMB  :',ISYMB
                  WRITE(LUPRI,*) 'In CC3_AAMATSDVIR: ISYALJ :',ISYALJ
                  WRITE(LUPRI,*) 'In CC3_AAMATSDVIR: ISYALJ2:',ISYALJ2
                  WRITE(LUPRI,*) 'In CC3_AAMATSDVIR: ISYMBD :',ISYMBD
                  WRITE(LUPRI,*) 'In CC3_AAMATSDVIR: ISCKIJ :',ISCKIJ
                  WRITE(LUPRI,*) 'In CC3_AAMATSDVIR: ISWMATX :',ISWMATX
                  WRITE(LUPRI,*) 'In CC3_AAMATSDVIR: ISWMATY :',ISWMATY

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
               KDIAGWX  = KEND3
               KWMATX   = KDIAGWX  + NCKIJ(ISWMATX)
               KINDSQWX = KWMATX   + NCKIJ(ISWMATX)
               KEND3  = KINDSQWX + (6*NCKIJ(ISWMATX) - 1)/IRAT + 1
C
               KDIAGWY  = KEND3
               KWMATY   = KDIAGWY  + NCKIJ(ISWMATY)
               KINDSQWY = KWMATY   + NCKIJ(ISWMATY)
               KEND3  = KINDSQWY + (6*NCKIJ(ISWMATY) - 1)/IRAT + 1
C
               KINDSQ  = KEND3
               KINDEX  = KINDSQ  + (6*NCKIJ(ISCKIJ) - 1)/IRAT + 1
               KINDEX2 = KINDEX  + (NCKI(ISYALJ) - 1)/IRAT + 1
               KEND3   = KINDEX2 + (NCKI(ISYALJ2) - 1)/IRAT + 1
C
               MAXXY    = MAX(NCKIJ(ISWMATX),NCKIJ(ISWMATY))
               MAXXXYY    = MAX(MAXXY,NCKIJ(ISWMATXY))
C
               KTMATW   = KEND3
               KT3VBG1  = KTMATW   + MAX(MAXXXYY,NCKIJ(ISCKIJ))
               KT3VBG2  = KT3VBG1  + NCKATR(ISCKD2)
               KT3VBG3 = KT3VBG2  + NCKATR(ISCKD2)
               KEND4   = KT3VBG3 + NCKATR(ISCKD2)
               LWRK4   = LWORK   - KEND4
C
               KW3XVDGX2 = KEND4
               KEND4   = KW3XVDGX2 + NCKATR(ISCKD2X)
C
               KW3YVDGY2 = KEND4
               KEND4   = KW3YVDGY2 + NCKATR(ISCKD2Y)
C
               KT3VBGX3 = KEND4
               KEND4    = KT3VBGX3 + NCKATR(ISCKD2X)
C
               KT3VBGY3 = KEND4
               KEND4    = KT3VBGY3 + NCKATR(ISCKD2Y)
C
               KINTVI  = KEND4
               KEND4   = KINTVI  + MAX(NCKA(ISCKD2X),NCKA(ISCKD2Y))
               LWRK4   = LWORK   - KEND4
C
               KINDSQWXY = KEND4
               KDIAGWXY  = KINDSQWXY + (6*NCKIJ(ISWMATXY) - 1)/IRAT + 1
               KWMATXY   = KDIAGWXY  + NCKIJ(ISWMATXY)
               KEND4     = KWMATXY   + NCKIJ(ISWMATXY)
               LWRK4     = LWORK     - KEND4
C
               IF (LWRK4 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Memory available : ',LWORK
                  WRITE(LUPRI,*) 'Memory needed    : ',KEND4
                  CALL QUIT('Insufficient space in CC3_AAMATSDVIR')
               END IF
C
C---------------------------------------------
C           Construct part of the diagonal.
C---------------------------------------------
C
               CALL CC3_DIAG(WORK(KDIAG),WORK(KFOCKD),ISCKIJ)
               CALL CC3_DIAG(WORK(KDIAGWX),WORK(KFOCKD),ISWMATX)
               CALL CC3_DIAG(WORK(KDIAGWY),WORK(KFOCKD),ISWMATY)
               CALL CC3_DIAG(WORK(KDIAGWXY),WORK(KFOCKD),ISWMATXY)
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
               LENSQWX = NCKIJ(ISWMATX)
               CALL CC3_INDSQ(WORK(KINDSQWX),LENSQWX,ISWMATX) 
               LENSQWY = NCKIJ(ISWMATY)
               CALL CC3_INDSQ(WORK(KINDSQWY),LENSQWY,ISWMATY) 
               LENSQWXY = NCKIJ(ISWMATXY)
               CALL CC3_INDSQ(WORK(KINDSQWXY),LENSQWXY,ISWMATXY) 


               DO B = 1,NVIR(ISYMB)
C
C------------------------------------
C           Initialise
C---------------------------------------
C
               CALL DZERO(WORK(KWMATX),NCKIJ(ISWMATX))
               CALL DZERO(WORK(KWMATY),NCKIJ(ISWMATY))
               CALL DZERO(WORK(KWMATXY),NCKIJ(ISWMATXY))
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
C           Read X transformed virtual integrals (H^X) used in W^X
C--------------------------------------------------------------------
C
                  IF (LISTRY(1:3).EQ.'R1 ') THEN
                    IOFF = ICKBD(ISCKD2X,ISYMB) + 
     &                     NCKATR(ISCKD2X)*(B - 1) + 1
                    IF (NCKATR(ISCKD2X) .GT. 0) THEN
                     CALL GETWA2(LUDKBCRX,FNDKBCRX,WORK(KW3XVDGX2),IOFF,
     &                           NCKATR(ISCKD2X))
                    ENDIF
                  END IF
C
C--------------------------------------------------------------------
C           Read Y transformed virtual integrals (H^Y) used in W^Y
C--------------------------------------------------------------------
C
                  IOFF = ICKBD(ISCKD2Y,ISYMB) + 
     &                   NCKATR(ISCKD2Y)*(B - 1) + 1
                  IF (NCKATR(ISCKD2Y) .GT. 0) THEN
                     CALL GETWA2(LUDKBCRY,FNDKBCRY,WORK(KW3YVDGY2),IOFF,
     &                           NCKATR(ISCKD2Y))
                  ENDIF
C
C--------------------------------------------------------------------
C           Read virtual integrals [H,T1X] where X is LISTRX (used in W^X)
C--------------------------------------------------------------------
C
                  IF (LISTRY(1:3).EQ.'R1 ') THEN
                    IF (NCKA(ISCKD2X) .GT. 0) THEN
                       IOFF = ICKAD(ISCKD2X,ISYMB) +
     &                        NCKA(ISCKD2X)*(B - 1) + 1
                       CALL GETWA2(LUDELDRX,FNDELDRX,WORK(KINTVI),IOFF,
     *                      NCKA(ISCKD2X))
                    ENDIF
C
                    CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KT3VBGX3),
     *                               WORK(KLAMDH),ISYMB,B,ISINTRX2,
     *                               WORK(KEND4),LWRK4)
                  END IF
C
C--------------------------------------------------------------------
C           Read virtual integrals [H,T1Y] where Y is LISTRY (used in W^Y)
C--------------------------------------------------------------------
C
                  IF (NCKA(ISCKD2Y) .GT. 0) THEN
                     IOFF = ICKAD(ISCKD2Y,ISYMB) +
     &                      NCKA(ISCKD2Y)*(B - 1) + 1
                     CALL GETWA2(LUDELDRY,FNDELDRY,WORK(KINTVI),IOFF,
     *                    NCKA(ISCKD2Y))
                  ENDIF
C
                  CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KT3VBGY3),
     *                             WORK(KLAMDH),ISYMB,B,ISINTRY2,
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
C Calculate WBD^X
C------------------------------------------------------
C

                  ! Copy T3^0 from KT3MAT to KTMATW. KTMATW used later
                  ! as help array and KTMATW overwritten.
c                 CALL DCOPY(NCKIJ(ISCKIJ),WORK(KT3MAT),1,
c    *                       WORK(KTMATW),1)
C

                  IF (LISTRY(1:3).EQ.'R1 ') THEN
C
C------------------------------------------------------
C     Calculate the  term <mu3|[X,T30]|HF> occupied contribution 
C     added in W^BD(a,i-,k-,j-), where i- means one-index transformation of i
C     (the three contributions are added: fixed BD, fixed aB, fixed Da 
C      added in W^BD).
C------------------------------------------------------
C
                     AIBJCK_PERM = 4 
C     write(lupri,*)'entering WXBD_O... '
                     CALL WXBD_O(AIBJCK_PERM,WORK(KT3MAT),ISCKIJ,
     *                           WORK(KFCKX),ISYMRX,
     *                           WORK(KWMATX),ISWMATX,WORK(KEND4),LWRK4)
C
C------------------------------------------------------
C     Calculate the  term <mu3|[[X,T2],T2]|HF> 
C     added in W^BD(a,i-,k-,j-).
C------------------------------------------------------
C
C     write(lupri,*)'entering WXBD_T2... '
                     CALL WXBD_T2(AIBJCK_PERM,B,ISYMB,D,ISYMD,
     *                            WORK(KT2TP),
     *                            ISYMT2,WORK(KFCKX),
     *                            ISYMRX,WORK(KINDSQWX),LENSQWX,
     *                            WORK(KWMATX),ISWMATX,WORK(KEND4),
     *                            LWRK4)

C
C----------------------------------------------------------------
C     Calculate the  terms <mu3|[H^X,T2]|HF> + <mu3|[H,T2^X]|HF>
C     added in W^BD(a,i-,k-,j-).
C----------------------------------------------------------------
C
C     write(lupri,*)'entering WXBD_GROUND(1)... '
                     CALL WXBD_GROUND(AIBJCK_PERM,
     *                        WORK(KT2X),ISYMRX,WORK(KTMATW),
     *                        WORK(KT3VDG1),WORK(KT3VBG1),WORK(KT3VBG3),
     *                        WORK(KT3VDG3),
     *                        WORK(KT3OG1),ISINT2,
     *                        WORK(KWMATX),WORK(KEND4),LWRK4,
     *                        WORK(KINDSQWX),LENSQWX,
     *                        ISYMB,B,ISYMD,D)
C
C     write(lupri,*)'entering WXBD_GROUND(2)... '
                     CALL WXBD_GROUND(AIBJCK_PERM,
     *                          WORK(KT2TP),ISYMT2,WORK(KTMATW),
     *                          WORK(KW3XVDGX1),WORK(KW3XVDGX2),
     *                          WORK(KT3VBGX3),WORK(KT3VDGX3),
     *                          WORK(KW3XOGX1),ISINTRX2,
     *                          WORK(KWMATX),WORK(KEND4),LWRK4,
     *                          WORK(KINDSQWX),LENSQWX,
     *                          ISYMB,B,ISYMD,D)
C
C------------------------------------------------
C     Divide by the energy difference and
C     remove the forbidden elements
C------------------------------------------------
C
                     CALL WBD_DIA(B,ISYMB,D,ISYMD,FREQRX,ISWMATX,
     *                          WORK(KWMATX),WORK(KDIAGWX),WORK(KFOCKD))
                     CALL T3_FORBIDDEN(WORK(KWMATX),ISYMRX,ISYMB,B,
     *                                 ISYMD,D)

c
c       call sum_pt3(work(KWMATX),isymb,b,isymd,d,
c    *             ISWMATX,work(kx3am),5)


C
c                 CALL GET_T3X_BD(ISYMRX,WORK(KWMATX),ISWMATX,
c    *                            WORK(KTMATW),ISCKIJ,WORK(KFCKX),
c    *                            ISYMRX,WORK(KT2TP),ISYMT2,
c    *                            WORK(KT2X),ISYMRX,WORK(KT3VDG1),
c    *                            WORK(KT3VBG1),WORK(KT3OG1),
c    *                            ISINT2,WORK(KW3XVDGX1),
c    *                            WORK(KW3XVDGX2),WORK(KW3XOGX1),
c    *                            ISINTRX2,FREQRX,WORK(KDIAGWX),
c    *                            WORK(KFOCKD),WORK(KINDSQWX),LENSQWX,
c    *                            B,ISYMB,D,ISYMD,WORK(KEND4),LWRK4)

C                 Calculate WYX intermmediate for <mu3|[Y,wX^{aBC}_{ijk}]|HF>
C                 where wX^{aBC}_{ijk} is part of T3^X with the virtual
C                 contribution (WXBD_V routine) taken out.

c
                     CALL WBD_O(WORK(KWMATX),ISWMATX,WORK(KFCKY),ISYMRY,
     *                    WORK(KWMATXY),ISWMATXY,WORK(KEND4),LWRK4)
C
                     CALL WBD_V(WORK(KWMATX),ISWMATX,WORK(KFCKY),ISYMRY,
     *                    WORK(KWMATXY),ISWMATXY,WORK(KEND4),LWRK4)
C
C                    Calculate <mu3|[[Y,T1X],T30]|HF>
C
                     CALL WBD_O(WORK(KT3MAT),ISCKIJ,WORK(KFCKYXO),
     *                    ISYMXY,
     *                    WORK(KWMATXY),ISWMATXY,WORK(KEND4),LWRK4)
C
                     CALL WBD_V(WORK(KT3MAT),ISCKIJ,WORK(KFCKYXV),
     *                    ISYMXY,
     *                    WORK(KWMATXY),ISWMATXY,WORK(KEND4),LWRK4)
C
C                    Calculate <mu3|[[Y,T2X],T20]|HF>
C
                     T2XNET2Z = .TRUE.
                     CALL WBD_T2(T2XNET2Z,B,ISYMB,D,ISYMD,
     *                           WORK(KT2X),ISYMRX,WORK(KT2TP),ISYMT2,
     *                           WORK(KFCKY),ISYMRY,
     *                           WORK(KINDSQWXY),LENSQWXY,WORK(KWMATXY),
     *                           ISWMATXY,WORK(KEND4),LWRK4)
C
C------------------------------------------------
C     Divide by the energy difference and
C     remove the forbidden elements (here only for debugging)
C------------------------------------------------
C
c                 CALL WBD_DIA(B,ISYMB,D,ISYMD,FREQXY,ISWMATXY,
c    *                        WORK(KWMATXY),WORK(KDIAGWXY),WORK(KFOCKD))
c                 CALL T3_FORBIDDEN(WORK(KWMATXY),ISYMXY,ISYMB,B,
c    *                              ISYMD,D)

 
c       call sum_pt3(work(KWMATXY),isymb,b,isymd,d,
c    *             ISWMATXY,work(kx3am),5)


                   END IF !LISTRY .eq. 'R1 '
C
C------------------------------------------------------
C Calculate WBD^Y
C------------------------------------------------------
C

                  ! Copy T3^0 from KT3MAT to KTMATW. KTMATW used later
                  ! as help array and KTMATW overwritten.
c                 CALL DCOPY(NCKIJ(ISCKIJ),WORK(KT3MAT),1,
c    *                       WORK(KTMATW),1)
C
C
C------------------------------------------------------
C     Calculate the  term <mu3|[Y,T30]|HF> occupied contribution 
C     added in W^BD(a,i-,k-,j-), where i- means one-index transformation of i
C     (the three contributions are added: fixed BD, fixed aB, fixed Da 
C      added in W^BD).
C------------------------------------------------------
C
                  AIBJCK_PERM = 4 ! transform all occupied indecies
c     write(lupri,*)'entering WXBD_O... Y'
                  IF (LISTRY(1:3).EQ.'R1 ') THEN
                     CALL WXBD_O(AIBJCK_PERM,WORK(KT3MAT),ISCKIJ,
     *                           WORK(KFCKY),ISYMRY,
     *                           WORK(KWMATY),ISWMATY,WORK(KEND4),LWRK4)
C
C------------------------------------------------------
C     Calculate the  term <mu3|[[Y,T2],T2]|HF> 
C     added in W^BD(a,i-,k-,j-).
C------------------------------------------------------
C
c     write(lupri,*)'entering WXBD_T2... Y'
                     CALL WXBD_T2(AIBJCK_PERM,B,ISYMB,D,ISYMD,
     *                           WORK(KT2TP),
     *                           ISYMT2,WORK(KFCKY),
     *                           ISYMRY,WORK(KINDSQWY),LENSQWY,
     *                           WORK(KWMATY),ISWMATY,WORK(KEND4),LWRK4)

C
                  END IF
C
C----------------------------------------------------------------
C     Calculate the  terms <mu3|[H^Y,T2]|HF> + <mu3|[H,T2^Y]|HF>
C     added in W^BD(a,i-,k-,j-).
C----------------------------------------------------------------
C
c     write(lupri,*)'entering WXBD_GROUND(1)... Y'
                  CALL WXBD_GROUND(AIBJCK_PERM,
     *                       WORK(KT2Y),ISYMRY,WORK(KTMATW),
     *                       WORK(KT3VDG1),WORK(KT3VBG1),WORK(KT3VBG3),
     *                       WORK(KT3VDG3),
     *                       WORK(KT3OG1),ISINT2,
     *                       WORK(KWMATY),WORK(KEND4),LWRK4,
     *                       WORK(KINDSQWY),LENSQWY,
     *                       ISYMB,B,ISYMD,D)
C
c     write(lupri,*)'entering WXBD_GROUND(2)... Y'
                  CALL WXBD_GROUND(AIBJCK_PERM,
     *                       WORK(KT2TP),ISYMT2,WORK(KTMATW),
     *                       WORK(KW3YVDGY1),WORK(KW3YVDGY2),
     *                       WORK(KT3VBGY3),WORK(KT3VDGY3),
     *                       WORK(KW3YOGY1),ISINTRY2,
     *                       WORK(KWMATY),WORK(KEND4),LWRK4,
     *                       WORK(KINDSQWY),LENSQWY,
     *                       ISYMB,B,ISYMD,D)
C
C------------------------------------------------
C     Divide by the energy difference and
C     remove the forbidden elements
C------------------------------------------------
C
                  CALL WBD_DIA(B,ISYMB,D,ISYMD,FREQRY,ISWMATY,
     *                         WORK(KWMATY),WORK(KDIAGWY),WORK(KFOCKD))
                  CALL T3_FORBIDDEN(WORK(KWMATY),ISYMRY,ISYMB,B,
     *                              ISYMD,D)

c
c       call sum_pt3(work(KWMATY),isymb,b,isymd,d,
c    *             ISWMATY,work(kx3am),5)


c                 CALL GET_T3X_BD(ISYMRY,WORK(KWMATY),ISWMATY,
c    *                            WORK(KTMATW),ISCKIJ,WORK(KFCKY),
c    *                            ISYMRY,WORK(KT2TP),ISYMT2,
c    *                            WORK(KT2Y),ISYMRY,WORK(KT3VDG1),
c    *                            WORK(KT3VBG1),WORK(KT3OG1),
c    *                            ISINT2,WORK(KW3YVDGY1),
c    *                            WORK(KW3YVDGY2),WORK(KW3YOGY1),
c    *                            ISINTRY2,FREQRY,WORK(KDIAGWY),
c    *                            WORK(KFOCKD),WORK(KINDSQWY),LENSQWY,
c    *                            B,ISYMB,D,ISYMD,WORK(KEND4),LWRK4)

C                 Calculate WXY intermmediate for <mu3|[X,wY^{aBC}_{ijk}]|HF>
C                 where wY^{aBC}_{ijk} is part of T3^Y with the virtual
C                 contribution (WXBD_V routine) taken out.

                  CALL WBD_O(WORK(KWMATY),ISWMATY,WORK(KFCKX),ISYMRX,
     *                 WORK(KWMATXY),ISWMATXY,WORK(KEND4),LWRK4)
 
                  CALL WBD_V(WORK(KWMATY),ISWMATY,WORK(KFCKX),ISYMRX,
     *                 WORK(KWMATXY),ISWMATXY,WORK(KEND4),LWRK4)
C
C                 Calculate <mu3|[[X,T1Y],T30]|HF>
C
                  CALL WBD_O(WORK(KT3MAT),ISCKIJ,WORK(KFCKXYO),ISYMXY,
     *                 WORK(KWMATXY),ISWMATXY,WORK(KEND4),LWRK4)
 
                  CALL WBD_V(WORK(KT3MAT),ISCKIJ,WORK(KFCKXYV),ISYMXY,
     *                 WORK(KWMATXY),ISWMATXY,WORK(KEND4),LWRK4)
C
C                 Calculate <mu3|[[X,T2Y],T20]|HF>
C

*     write(lupri,*) 'fock in cc3_aamat listry ',listry
*     call output(WORK(KFCKX),1,NORBT,1,NORBT,NORBT,NORBT,1,lupri)
*     write(lupri,*) 'E2 in cc3_aamat listry ',listry
*     call output(WORK(KT2Y),1,NVIRT*NRHFT*NRHFT,1,NVIRT,
*    * NVIRT*NRHFT*NRHFT,NVIRT,1,lupri)


                  T2XNET2Z = .TRUE.
                  CALL WBD_T2(T2XNET2Z,B,ISYMB,D,ISYMD,
     *                        WORK(KT2Y),ISYMRY,WORK(KT2TP),ISYMT2,
     *                        WORK(KFCKX),ISYMRX,
     *                        WORK(KINDSQWXY),LENSQWXY,WORK(KWMATXY),
     *                        ISWMATXY,WORK(KEND4),LWRK4)
C
C------------------------------------------------
C     Divide by the energy difference and
C     remove the forbidden elements 
C------------------------------------------------
C
                  CALL T3_FORBIDDEN(WORK(KWMATXY),ISYMXY,ISYMB,B,
     *                              ISYMD,D)

c                 call sum_pt3(work(KWMATXY),isymb,b,isymd,d,
c    *                       ISWMATXY,work(kx3am),4)

*     write(lupri,*) ' cc3_aamat: omega_f - omega_X = ', FREQXY

*                 call sum_pt3(work(KWMATXY),isymb,b,isymd,d,
*    *                       ISWMATXY,work(kx3am),4)

                  CALL WBD_DIA(B,ISYMB,D,ISYMD,FREQXY,ISWMATXY,
     *                        WORK(KWMATXY),WORK(KDIAGWXY),WORK(KFOCKD))

c                 call sum_pt3(work(KWMATXY),isymb,b,isymd,d,
c    *                       ISWMATXY,work(kx3am),4)

c

C
C-----------------------------------------------
C                 Carry out the contractions...
C-----------------------------------------------
C
C
C--------------------------------------------------------
C                 Calculate the  term <mu1|[H,WBD^XY(3)]|HF> 
C                 added in OMEGA1EFF 
C--------------------------------------------------------
C
c     write(lupri,*) 'XIAJB in cc3_aamat listry:i ',listry
c     call outpak(WORK(KXIAJB),NT1AM(1),1,lupri)


                  CALL CC3_CY1(OMEGA1EFF,ISYRES,WORK(KWMATXY),ISWMATXY,
     *                         WORK(KTMATW),WORK(KXIAJB),
     *                         ISINT1,WORK(KINDSQWXY),LENSQWXY,
     *                         WORK(KEND4),LWRK4,
     *                         ISYMB,B,ISYMD,D)

C
C------------------------------------------------------
C                 Calculate the  term <mu2|[H,WBD^XY(3)]|HF> 
C                 ( Fock matrix cont ) 
C                 added in OMEGA2EFF 
C------------------------------------------------------
C
                  CALL CC3_CY2F(OMEGA2EFF,ISYRES,WORK(KWMATXY),ISWMATXY,
     *                          WORK(KTMATW),WORK(KFOCK0),ISYM0,
     *                          WORK(KINDSQWXY),
     *                          LENSQWXY,WORK(KEND4),LWRK4,
     *                          ISYMB,B,ISYMD,D)
c                 
C------------------------------------------------------
C                 Calculate the  term <mu2|[H,WBD^XY(3)]|HF> 
C                 ( Occupied  cont ) 
C                 added in OMEGA2EFF 
C------------------------------------------------------
C
                 CALL CC3_CY2O(OMEGA2EFF,ISYRES,WORK(KWMATXY),ISWMATXY,
     *                          WORK(KTMATW),WORK(KTROC),
     *                          WORK(KTROC1),ISINT1,WORK(KEND4),LWRK4,
     *                          WORK(KINDSQWXY),LENSQWXY,
     *                          ISYMB,B,ISYMD,D)
C
C------------------------------------------------------
C                 Calculate the  term <mu2|[H,WBD^XY(3)]|HF> 
C                 ( Virtual  cont ) 
C                 added in OMEGA2EFF 
C------------------------------------------------------
C
                  CALL CC3_CY2V(OMEGA2EFF,ISYRES,WORK(KRBJIAXY),
     *                          WORK(KWMATXY),
     *                          ISWMATXY,WORK(KTMATW),WORK(KTRVI),
     *                          WORK(KTRVI1),ISINT1,WORK(KEND4),LWRK4,
     *                          WORK(KINDSQWXY),LENSQWXY,
     *                          ISYMB,B,ISYMD,D)
C
               END DO !B
            END DO    !ISYMB
C
         END DO       !D
      END DO          !ISYMD
c
c      

c

C
C
C------------------------------------------------------
C     Accumulate RBJIA from <mu2|[H,W^BD(3)]|HF> ( Virtual  cont ) 
C     in OMEGA2EFF 
C------------------------------------------------------
C
      CALL CC3_RBJIA(OMEGA2EFF,ISYRES,WORK(KRBJIAXY))
C
      IF (IPRINT .GT. 55) THEN
        XNORMVAL = DDOT(NT2AM(ISYRES),OMEGA2EFF,1,OMEGA2EFF,1)
        WRITE(LUPRI,*)'CC3_AAMATSDVIR:LISTRX,IDLSTRX,LISTRY,IDLSTRY ',
     *       LISTRX,IDLSTRX,LISTRY,IDLSTRY
        WRITE(LUPRI,*)'CC3_AAMATSDVIR:Norm of OMEGA2EFF final ',XNORMVAL
      ENDIF
C
      IF (IPRINT .GT. 55) THEN
         XNORMVAL = DDOT(NT1AM(ISYRES),OMEGA1EFF,1,OMEGA1EFF,1)
         WRITE(LUPRI,*) 'CC3_AAMATSDVIR:Norm of OMEGA1EFF final  ',
     *                   XNORMVAL
      ENDIF
C
c     write(lupri,*)' OMEGA1EFF after cc3_aamatvir: listry ', listry
c     call output(OMEGA1EFF,1,nvirt,1,nrhft,nvirt,nrhft,1,lupri)


C
C
COMMENT COMMENT COMMENT COMMENT 
COMMENT COMMENT COMMENT COMMENT 
C     write(lupri,*)'omega1 '
C     do i = 1,nt1am(isyres)
C        write(lupri,*) i , OMEGA1(i)
C     end do
COMMENT COMMENT COMMENT COMMENT 
COMMENT COMMENT COMMENT COMMENT 
COMMENT COMMENT
COMMENT COMMENT
c      write(lupri,*) 'cont to t3x in CC3_AAMATSDVIR'
c      call print_pt3(work(kx3am),1,4)
C      write(lupri,*) 'The summed S terms : '
C      call print_pt3(work(kx3am),1,1)
COMMENT COMMENT
COMMENT COMMENT
C
C--------------------------------
C     Close files for "response"
C--------------------------------
C
      CALL WCLOSE2(LU3SRTR,FN3SRTR,'DELETE')
      CALL WCLOSE2(LUCKJDRX,FNCKJDRX,'DELETE')
      CALL WCLOSE2(LUCKJDRY,FNCKJDRY,'DELETE')
      CALL WCLOSE2(LUDELDRX,FNDELDRX,'DELETE')
      CALL WCLOSE2(LUDELDRY,FNDELDRY,'DELETE')
      CALL WCLOSE2(LUDKBCRX,FNDKBCRX,'DELETE')
      CALL WCLOSE2(LUDKBCRY,FNDKBCRY,'DELETE')
C
C-------------
C     End
C-------------
C
      CALL QEXIT('AMATSDV')
C
      RETURN
      END
C  /* Deck tetax_jk_i */
      SUBROUTINE TETAX_JK_I(WJK,ISWJK,XOP,ISYMXOP,TETAXJK,ISTETAXJK,
     *                      WORK,LWORK)
C
C TETAXJK(bcai) = TETAXJK(bcai) 
C
C             - xop(li) wjk(bcal)
C                             

      IMPLICIT NONE
C
      INTEGER ISWJK, ISYMXOP, ISTETAXJK, LWORK
      INTEGER KFCLI,KEND0,LWRK0 
      INTEGER ISYMI,ISYML,KOFF1,KOFF2,ISYMBCA,KOFF3,NTOTBCA,NTOTL
C
#if defined (SYS_CRAY)
      REAL WJK(*), TETAXJK(*), XOP(*), WORK(LWORK)
      REAL ONE
      real xnormval,ddot
#else
      DOUBLE PRECISION WJK(*), TETAXJK(*), XOP(*), WORK(LWORK) 
      DOUBLE PRECISION ONE
      double precision xnormval,ddot
#endif
C
      PARAMETER (ONE = 1.0D0)
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccsdinp.h>
C
      CALL QENTER('TETJKI')

C
C RESORT OCC-OCC  FOCKY ELEMENTS (L,I)
C
C
      KFCLI  = 1
      KEND0  = KFCLI + NMATIJ(ISYMXOP)
      LWRK0  = LWORK  - KEND0
C
      IF (LWRK0 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWRK0
         WRITE(LUPRI,*) 'Memory needed    : ',KEND0
         CALL QUIT('Insufficient space in TETAX_JK_A')
      END IF
C   
      DO ISYMI = 1,NSYM
         ISYML = MULD2H(ISYMI,ISYMXOP)
         DO I = 1,NRHF(ISYMI)
             KOFF1 = IFCRHF(ISYML,ISYMI) + NORB(ISYML)*(I - 1) + 1
             KOFF2 = KFCLI + IMATIJ(ISYML,ISYMI) + NRHF(ISYML)*(I - 1)
             CALL DCOPY(NRHF(ISYML),XOP(KOFF1),1,WORK(KOFF2),1)
         END DO
      END DO
C
      DO ISYMI = 1,NSYM
         ISYML   = MULD2H(ISYMXOP,ISYMI)
         ISYMBCA = MULD2H(ISWJK,ISYML)
C
         KOFF1   = 1
     *           + IMAABCI(ISYMBCA,ISYML)
         KOFF2   = KFCLI 
     *           + IMATIJ(ISYML,ISYMI)
         KOFF3   = 1
     *           + IMAABCI(ISYMBCA,ISYMI)
C
         NTOTBCA = MAX(1,NMAABC(ISYMBCA))
         NTOTL   = MAX(1,NRHF(ISYML))
C
C TETAXJK(bcai) = TETAXJK(bcai) -  wjk(bcal) * xop(li)
C
         CALL DGEMM('N','N',NMAABC(ISYMBCA),NRHF(ISYMI),
     *              NRHF(ISYML),ONE,WJK(KOFF1),NTOTBCA,
     *              WORK(KOFF2),NTOTL,
     *              ONE,TETAXJK(KOFF3),NTOTBCA)
      END DO
C                             
      CALL QEXIT('TETJKI')
      RETURN
      END
C  /* Deck cc3_aamatsdocc */
      SUBROUTINE CC3_AAMATSDOCC(LISTRX,IDLSTRX,LISTRY,IDLSTRY,
     &                          OMEGA1,OMEGA2,
     &                          OMEGA1EFF,OMEGA2EFF,
     *                          ISYRES,WORK,LWORK,
     *                          LUCKJD,FNCKJD,LUTOC,FNTOC,
     *                          LUDELD,FNDELD,LU3FOP,FN3FOP,
     *                          LU3VI,FN3VI)
*
*******************************************************************
*
* This routine calculates the occupied part of the right hand side vector T^XY
* for the paritioned second-order amplitude equations.
*
* Calculate the A^X matrix contribution (using theta intermmediate): 
*
* THETAXY^LK =  <mu3|[X,T3^Y]|HF> (T3^Y also calculated using theta^LK)
*
* projected up into the Singles-and-Doubles space:
*
* omega1eff = omega1eff - <mu1|[H,THETAXY^LK]|HF>
* omega2eff = omega2eff - <mu2|[H,THETAXY^LK]|HF>
*
*
*
*     Written by Poul Jorgensen and Filip Pawlowski, Aarhus, Spring 2003
*
*******************************************************************
*
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
      CHARACTER*3 LISTRX, LISTRY
      CHARACTER*8 LABELRX,LABELRY
      INTEGER     IDLSTRX,IDLSTRY
      INTEGER     ISYMRX, ISYMRY
C
      INTEGER ISYRES, LWORK
      INTEGER ISYM0,ISYMT2,ISYMT3
      INTEGER KT2TP,KFOCKD,KLAMDP,KLAMDH,KEND00,LWRK00,KXIAJB,KFOCK0
      INTEGER IOPTTCME,IOPT
      INTEGER ISYMWX,ISYMWY,ISYMWXY
      INTEGER ISYFCKX,ISYFCKY,ISYMXY,KFCKX,KFCKY,KEND1,LWRK1
      INTEGER ISINT1,ISINT2,KTMP,KT3BOG2
      INTEGER KT3OG2,KT3VIJG1,KXGADCK
      INTEGER ISYMTETAX,ISYMTETAY,ISYMTETAXY
      INTEGER ISYMK,ISYML,ISYMKL,ISYT30KL,ISTETAXKL,ISTETAYKL,ISTETAXYKL
      INTEGER KT30KL,KTETAXKL,KTETAYKL,KTETAXYKL,KTMAT
      INTEGER KINTOC,KEND2,LWRK2,IOFF
c
      integer kx3am
c
#if defined (SYS_CRAY)
      REAL OMEGA1(*),OMEGA2(*),OMEGA1EFF(*),OMEGA2EFF(*)
      REAL WORK(LWORK)
      REAL FREQRX,FREQRY,FREQXY
      REAL XNORMVAL,DDOT
#else
      DOUBLE PRECISION OMEGA1(*),OMEGA2(*),OMEGA1EFF(*),OMEGA2EFF(*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION FREQRX,FREQRY,FREQXY
      DOUBLE PRECISION XNORMVAL,DDOT
#endif
C
      CALL QENTER('AASDOCC')
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
      KFOCK0 = KXIAJB + NT2AM(ISINT1)
      KEND00 = KFOCK0 + N2BST(ISYM0)
      LWRK00 = LWORK  - KEND00
C
      IF (LWRK00 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_AAMATSDOCC (zeroth allo.')
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
C------------------------------------------------------------------
C     Read in Fock matrix in AO basis (from file) and transform to    
C     Lambda_0 basis.             
C------------------------------------------------------------------
C
      CALL GET_FOCK0(WORK(KFOCK0),WORK(KLAMDP),WORK(KLAMDH),
     *               WORK(KEND00),LWRK00)
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
            call quit('Insufficient space (T3) in CC3_AAMATSDOCC')
         END IF
      endif
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
C
C determining R1 labels
C
      IF (LISTRX(1:3).EQ.'R1 ') THEN
         ISYMRX = ISYLRT(IDLSTRX)
         FREQRX  = FRQLRT(IDLSTRX)
         LABELRX = LRTLBL(IDLSTRX)
C
      ELSE
         CALL QUIT('Unknown list for LISTRX in CC3_AAMATSDOCC.') 
      END IF
C      
      IF (LISTRY(1:3).EQ.'R1 ') THEN
         ISYMRY = ISYLRT(IDLSTRY)
         FREQRY  = FRQLRT(IDLSTRY)
         LABELRY = LRTLBL(IDLSTRY)
C
      ELSE
         CALL QUIT('Unknown list for LISTRY in CC3_AAMATSDVIR.') 
      END IF
C
      FREQXY = FREQRX + FREQRY
C
      ISYMWX   = MULD2H(ISYMT3,ISYMRX)
      ISYMWY   = MULD2H(ISYMT3,ISYMRY)
C
      ISYMWXY   = MULD2H(ISYMWX,ISYMRY)
C
      IF (ISYRES.NE.ISYMWXY) THEN
         WRITE(LUPRI,*) 'INCONSISTENT SYMMETRY SPECIFICATION'
         WRITE(LUPRI,*) 
     *      'ISYRES =',ISYRES,'ISYMWXY =',ISYMWXY
         STOP' CC3_AAMATSDVIR inconsistent symmetry specification' 
      END IF
C
      ISYFCKX = ISYMRX
      ISYFCKY = ISYMRY
      ISYMXY = MULD2H(ISYMRX,ISYMRY)
C
C-------------------------------------------------
C     Read property matrices and Lambda matrices
C-------------------------------------------------
C
      KFCKX  = KEND00
      KFCKY  = KFCKX  + N2BST(ISYFCKX)
      KEND1  = KFCKY  + N2BST(ISYFCKY)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. NT2SQ(ISYMRX)) THEN
         CALL QUIT('Out of memory in CC3_AAMATSDOCC (TOT_T3Y) ')
      ENDIF
C
C------------------------------------------
C     Calculate the F^X matrix
C------------------------------------------
C
      CALL GET_FOCKX(WORK(KFCKX),LABELRX,IDLSTRX,ISYMRX,WORK(KLAMDP),
     *                  ISYM0,WORK(KLAMDH),ISYM0,WORK(KEND1),LWRK1)
C
C------------------------------------------
C     Calculate the F^Y matrix
C------------------------------------------
C
      CALL GET_FOCKX(WORK(KFCKY),LABELRY,IDLSTRY,ISYMRY,WORK(KLAMDP),
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
C
      KT3BOG2   = KEND1
      KT3OG2   = KT3BOG2   + NTRAOC(ISYM0)
      KEND1     = KT3OG2    + NTRAOC(ISINT2)
C
      KT3VIJG1  = KEND1
      KEND1     = KT3VIJG1  + NMAABCI(ISYM0)
      LWRK1     = LWORK     - KEND1
C
      KXGADCK   = KEND1 
      KEND1   = KXGADCK + NMAABCI(ISYM0)
      LWRK1     = LWORK     - KEND1
C
      KINTOC  = KEND1
      KTMP = KINTOC  + NTOTOC(ISYM0) !temporary storage
      KEND2   = KTMP + NTRAOC(ISYM0) !temporary storage
      LWRK2   = LWORK   - KEND2
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in CC3_AAMATSDOCC (3)')
      END IF
C
C-----------------------------------------------------------------
C     Construct occupied integrals which are required for occupied
C     part of projection (KT3BOG2)
C-----------------------------------------------------------------
C
c     This call has been replaced by explicit construction of KT3BOG2
c     integrals (the only ones needed) to safe the work space...
c
*     CALL INTOCC_T3BAR0(LUTOC,FNTOC,WORK(KLAMDH),ISYM0,WORK(KT3BOG1),
*    *                   WORK(KT3BOL1),WORK(KT3BOG2),WORK(KT3BOL2),
*    *                   WORK(KEND1),LWRK1)
C
      IOFF = 1
C    ! Read in integrals (ia|j delta) sitting as KINTOC(ah ip | jp delta)
      IF (NTOTOC(ISYM0) .GT. 0) THEN
         CALL GETWA2(LUTOC,FNTOC,WORK(KINTOC),IOFF,NTOTOC(ISYM0))
      ENDIF
C
C    ! Transform the integrals with hole index to (ia|jk) and sort as 
C    ! T3BOG1(kija)
      CALL CC3_TROCC(WORK(KINTOC),WORK(KTMP),WORK(KLAMDH),
     *               WORK(KEND2),LWRK2,ISYM0)
C
C    ! (ia|jk) sorted as T3BOG1(kija) and now sorted as T3BOG2(ajik)
      CALL CCFOP_SORT(WORK(KTMP),WORK(KT3BOG2),ISYM0,1)
C
C-----------------------------------------------------------------
C     Construct occupied integrals which are required to calculate    
C     t3_0 amplitudes
C-----------------------------------------------------------------
C
      CALL INTOCC_T30(LUCKJD,FNCKJD,WORK(KLAMDP),ISINT2,WORK(KTMP),
     *                WORK(KT3OG2),WORK(KEND2),LWRK2)
!
!     FROM NOW ON KEND1 CAN BE USED (don't need KTMP nor KINTOC anymore)
!

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
C     Get virtual integrals for virtual part of projection (KXGADCK)
C     KXGADCK g(kcad) = (kc ! ad) sorted as I(adck)
C-------------------------------------------------------------------
C
      IOPT = 1
      CALL INTVIR_T3B0_JK(IOPT,WORK(KXGADCK),DUMMY,ISYM0,
     *                    WORK(KLAMDP),ISYM0,LU3VI,FN3VI,LU3FOP,FN3FOP,
     *                    WORK(KEND1),LWRK1)
C
C----------------------------
C     Loop over K
C----------------------------
C
      ISYMTETAX = MULD2H(ISYM0,ISYMRX)
      ISYMTETAY = MULD2H(ISYM0,ISYMRY)
      ISYMTETAXY = MULD2H(ISYM0,ISYMXY)
      DO ISYMK = 1,NSYM

         DO K = 1,NRHF(ISYMK)
C
            DO ISYML = 1,NSYM
C
               ISYMKL = MULD2H(ISYMK,ISYML)
               ISYT30KL = MULD2H(ISYMKL,ISYMT3)
               ISTETAXKL  = MULD2H(ISYMKL,ISYMTETAX)
               ISTETAYKL  = MULD2H(ISYMKL,ISYMTETAY)
               ISTETAXYKL  = MULD2H(ISYMKL,ISYMTETAXY)
C
               KT30KL = KEND1
               KTETAXKL  = KT30KL + NMAABCI(ISYT30KL)
               KTETAYKL   = KTETAXKL + NMAABCI(ISTETAXKL)
               KTETAXYKL   = KTETAYKL + NMAABCI(ISTETAYKL)
               KEND2   = KTETAXYKL + NMAABCI(ISTETAXYKL)
               LWRK2  = LWORK  - KEND2
C
               KTMAT = KEND2
               KEND2 = KTMAT + NMAABCI(ISTETAXYKL)
               LWRK2  = LWORK  - KEND2
C
               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Memory available : ',LWORK
                  WRITE(LUPRI,*) 'Memory needed    : ',KEND2
                  CALL QUIT('Insufficient space in CC3_AAMATSDOCC (4)')
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
     *                              WORK(KEND2),LWRK2)
C
                  CALL GET_T30_IJ_V(WORK(KT30KL),ISYT30KL,WORK(KT2TP),
     *                              ISYM0,WORK(KT3VIJG1),
     *                              ISYM0,ISYML,L,ISYMK,K,
     *                              WORK(KEND2),LWRK2)
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
                    WRITE(LUPRI,*)'NORM OF KT30KL IN CC3_AAMATSDOCC ', 
     *                             XNORMVAL
                  END IF
C
C---------------------------------
C                 Calculate KTETAXKL
C---------------------------------
C
                  CALL DZERO(WORK(KTETAXKL),NMAABCI(ISTETAXKL))
C
                  CALL TETAX_JK_BC(WORK(KT30KL),ISYT30KL,WORK(KFCKX),
     *                             ISYFCKX,WORK(KTETAXKL),ISTETAXKL,
     *                             WORK(KEND2),LWRK2)
C
c                 call w3jk_dia(work(ktetaxkl),istetaxkl,isyml,l,isymk,
c    *                          k,work(kfockd),freqr1)
c                 call t3_forbidden_jk(work(ktetaxkl),isymtetax,isyml,l,
c    *                                  isymk,k)
c
c                 call sum_pt3_jk(work(KTETAXKL),isyml,l,isymk,k,
c    *                            ISTETAXKL,
c    *                            work(kx3am),1)
c
c                 xnormval = ddot(nmaabci(istetaxkl),work(ktetaxkl),1,
c    *                              work(ktetaxkl),1)
c                 write(lupri,*)'norm of ktetaxkl in cc3_AAMATSDOCC ',
c    *                             xnormval
C
                  CALL TETAX_JK_A(WORK(KT30KL),ISYT30KL,WORK(KFCKX),
     *                             ISYFCKX,WORK(KTETAXKL),ISTETAXKL,
     *                             WORK(KEND2),LWRK2)
C
                  CALL W3JK_DIA(WORK(KTETAXKL),ISTETAXKL,ISYML,L,ISYMK,
     *                           K,WORK(KFOCKD),FREQRX)
                  CALL T3_FORBIDDEN_JK(WORK(KTETAXKL),ISYMTETAX,ISYML,L,
     *                                  ISYMK,K)

c                 call sum_pt3_jk(work(KTETAXKL),isyml,l,isymk,k,
c    *                            ISTETAXKL,
c    *                            work(kx3am),4)

C
                  IF (IPRINT .GT. 55) THEN
                    XNORMVAL = DDOT(NMAABCI(ISTETAXKL),WORK(KTETAXKL),1,
     *                              WORK(KTETAXKL),1)
                    WRITE(LUPRI,*)'(FINAL) NORM OF KTETAXKL ', XNORMVAL
                  END IF
C
C---------------------------------
C                 Calculate KTETAYKL
C---------------------------------
C
                  CALL DZERO(WORK(KTETAYKL),NMAABCI(ISTETAYKL))
C
                  CALL TETAX_JK_BC(WORK(KT30KL),ISYT30KL,WORK(KFCKY),
     *                             ISYFCKY,WORK(KTETAYKL),ISTETAYKL,
     *                             WORK(KEND2),LWRK2)
C
C
                  CALL TETAX_JK_A(WORK(KT30KL),ISYT30KL,WORK(KFCKY),
     *                             ISYFCKY,WORK(KTETAYKL),ISTETAYKL,
     *                             WORK(KEND2),LWRK2)
C
                  CALL W3JK_DIA(WORK(KTETAYKL),ISTETAYKL,ISYML,L,ISYMK,
     *                           K,WORK(KFOCKD),FREQRY)
                  CALL T3_FORBIDDEN_JK(WORK(KTETAYKL),ISYMTETAY,ISYML,L,
     *                                  ISYMK,K)

c                 call sum_pt3_jk(work(KTETAYKL),isyml,l,isymk,k,
c    *                            ISTETAYKL,
c    *                            work(kx3am),1)

C
                  IF (IPRINT .GT. 55) THEN
                    XNORMVAL = DDOT(NMAABCI(ISTETAYKL),WORK(KTETAYKL),1,
     *                              WORK(KTETAYKL),1)
                    WRITE(LUPRI,*)'(FINAL) NORM OF KTETAYKL ', XNORMVAL
C
                  END IF
C
C-------------------------------------------------------------
C                 Calculate KTETAXYKL = <mu3|[Y,KTETAXKL]|HF>
C-------------------------------------------------------------
C
                  CALL DZERO(WORK(KTETAXYKL),NMAABCI(ISTETAXYKL))
C
                  ! virtual contribution:
                  ! KTETAXYKL = sum_d Y_ad KTETAXKL^{d-,b-,c-}_{iKL}
                  CALL TETAX_JK_A(WORK(KTETAXKL),ISTETAXKL,WORK(KFCKY),
     *                             ISYFCKY,WORK(KTETAXYKL),ISTETAXYKL,
     *                             WORK(KEND2),LWRK2)
C
                  ! occupied contribution:
                  ! KTETAXYKL = sum_m Y_mi KTETAXKL^{a-,b-,c-}_{mKL}
                  CALL TETAX_JK_I(WORK(KTETAXKL),ISTETAXKL,WORK(KFCKY),
     *                             ISYFCKY,WORK(KTETAXYKL),ISTETAXYKL,
     *                             WORK(KEND2),LWRK2)
C
C-------------------------------------------------------------
C                 Calculate KTETAXYKL = <mu3|[X,KTETAYKL]|HF>
C-------------------------------------------------------------
C
                  ! virtual contribution:
                  ! KTETAXYKL = sum_d X_ad KTETAYKL^{d-,b-,c-}_{iKL}
                  CALL TETAX_JK_A(WORK(KTETAYKL),ISTETAYKL,WORK(KFCKX),
     *                             ISYFCKX,WORK(KTETAXYKL),ISTETAXYKL,
     *                             WORK(KEND2),LWRK2)
C
                  ! occupied contribution:
                  ! KTETAXYKL = sum_m X_mi KTETAYKL^{a-,b-,c-}_{mKL}
                  CALL TETAX_JK_I(WORK(KTETAYKL),ISTETAYKL,WORK(KFCKX),
     *                             ISYFCKX,WORK(KTETAXYKL),ISTETAXYKL,
     *                             WORK(KEND2),LWRK2)
C
C---------------------------------------------------------------
C                 Divide by orbital energy difference and remove
C                 forbidden elements
C---------------------------------------------------------------
C
                  CALL T3_FORBIDDEN_JK(WORK(KTETAXYKL),ISYMTETAXY,ISYML,
     *                                 L,ISYMK,K)

c                 call sum_pt3_jk(work(KTETAXYKL),isyml,l,isymk,k,
c    *                            ISTETAXYKL,
c    *                            work(kx3am),4)

                  CALL W3JK_DIA(WORK(KTETAXYKL),ISTETAXYKL,ISYML,L,
     *                          ISYMK,K,WORK(KFOCKD),FREQXY)

c                 call sum_pt3_jk(work(KTETAXYKL),isyml,l,isymk,k,
c    *                            ISTETAXYKL,
c    *                            work(kx3am),4)


C
                  IF (IPRINT .GT. 55) THEN
                    XNORMVAL = DDOT(NMAABCI(ISTETAXYKL),WORK(KTETAXYKL),
     *                              1,WORK(KTETAXYKL),1)
                    WRITE(LUPRI,*)'(FINAL) NORM OF KTETAXYKL ', XNORMVAL
C
                  END IF
C
C-------------------------------------------------------
C                 Project up to singles space:
C                 omega1eff = - <mu1|[H,THETAXY^LK]|HF>
C-------------------------------------------------------
C
                  CALL CC3_CY1_JK(OMEGA1EFF,ISYRES,WORK(KTETAXYKL),
     *                            ISTETAXYKL,WORK(KTMAT),
     *                            WORK(KXIAJB),ISINT1,WORK(KEND2),LWRK2,
     *                            ISYML,L,ISYMK,K)
C
C-------------------------------------------------------
C                 Project up to doubles space:
C                 omega2eff = - <mu2|[H,THETAXY^LK]|HF>
C-------------------------------------------------------
C
                  ! Fock matrix contribution
                  CALL CC3_CY2F_JK(OMEGA2EFF,ISYRES,WORK(KTETAXYKL),
     *                             ISTETAXYKL,WORK(KTMAT),
     *                             WORK(KFOCK0),ISYM0,WORK(KEND2),LWRK2,
     *                             ISYML,L,ISYMK,K)
C
                  ! occupied integrals contribution
                  CALL CC3_CY2O_JK(OMEGA2EFF,ISYRES,WORK(KTETAXYKL),
     *                             ISTETAXYKL,WORK(KTMAT),WORK(KT3BOG2),
     *                             ISYM0,WORK(KEND2),LWRK2,
     *                             ISYML,L,ISYMK,K)
                  ! virtual integrals contribution
                  CALL CC3_CY2V_JK(OMEGA2EFF,ISYRES,WORK(KTETAXYKL),
     *                             ISTETAXYKL,WORK(KTMAT),WORK(KXGADCK),
     *                             ISYM0,WORK(KEND2),LWRK2,
     *                             ISYML,L,ISYMK,K)
C
               ENDDO   ! L
            ENDDO      ! ISYML
         ENDDO       ! K
      ENDDO          ! ISYMK 
C
c      write(lupri,*) 'T30KL in CC3_AAMATSDOCC'
c      call print_pt3(work(kx3am),isym0,4)
C
C
C------------------------------------------------------
C add: omega contributions to omegaeff     
C------------------------------------------------------
C
c     write(lupri,*)'OMEGA1EFF in CC3_AASDOCC'
c     call output(OMEGA1EFF,1,nvir(1),1,nrhf(1),nvir(1),nrhf(1),1,lupri)
 
c     write(lupri,*)'OMEGA2EFF in CC3_AASDOCC'
c     call outpak(OMEGA2EFF,NT1AM(isyres),1,lupri)
 
      IF (IPRINT .GT. 55) THEN
         XNORMVAL = DDOT(NT2AM(ISYRES),OMEGA2EFF,1,OMEGA2EFF,1)
         WRITE(LUPRI,*)'Norm of OMEGA2EFF final ',XNORMVAL
      ENDIF
C 
      IF (IPRINT .GT. 55) THEN
         XNORMVAL = DDOT(NT1AM(ISYRES),OMEGA1EFF,1,OMEGA1EFF,1)
         WRITE(LUPRI,*) 'Norm of OMEGA1EFF final, OMEGA1EFF   ',
     *                   XNORMVAL
      ENDIF
c
C-------------
C     End
C-------------
C
      CALL QEXIT('AASDOCC')
C
      RETURN
      END
C /* Deck cc3_cy1_jk */
      SUBROUTINE CC3_CY1_JK(OMEGA1,ISOMG1,WMAT,ISWMAT,TMAT,
     *                   XIAJB,ISYINT,WORK,LWORK,
     *                   ISYMIJ,IJ,ISYMIK,IK)
C
C     P. Jorgensen and F. Pawlowski.         Spring 2003
C
      IMPLICIT NONE
C
C     General symmetry: ISWMAT is symmetry of WMAT and TMAT intermdiates.  
C                       (exclusive isymj,isymk) 
C                       ISYINT = XIAJB
C                       ISOMG1 = ISWMAT*ISYINT
C 
C        Omega1(ai)   
C
C        = sum_{bjck} ( t^{abc}_{ijk}  -  t^{cba}_{ijk} ) * L_{jbkc}
C  (1):
C        =  ( W^{JK}_{bcai} - W^{JK}_{baci} ) * L_{bJcK} 
C  (2):
C        +  ( W^{IJ}_{abck} - W^{IJ}_{cbak} ) * L_{ckbJ} 
C  (3):
C        +  ( W^{KI}_{cabj} - W^{KI}_{acbj} ) * L_{bjcK} 
C
C
C
     
      INTEGER ISOMG1, ISWMAT, ISYINT, LWORK, ISYMIJ, IJ, ISYMIK, IK
      INTEGER ISYMJ,ISYMK,ISYMJK,ISYMAI,KBC,KWBACI,KEND1,LWRK1 
      INTEGER ISYMI,ISYMBCA,ISYMA,ISYMBC,KOFF1,KOFF2,KOFF3,NTOTBC 
      INTEGER ISYMBCK,KBCK,KWCBAK,ISYMC,ISYMCK,ISYMB,ISYMABC,ISYMAB
      INTEGER NTOTA,ISYMCBJ,KCBJ,ISYMBJ,ISYMACB,ISYMAC,ISYMCB
#if defined (SYS_CRAY)
      REAL OMEGA1(*), WMAT(*), TMAT(*), XIAJB(*), WORK(*)
      REAL ONE
#else
      DOUBLE PRECISION OMEGA1(*), WMAT(*), TMAT(*), XIAJB(*), WORK(*)
      DOUBLE PRECISION ONE
#endif
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccsdinp.h>
C
      PARAMETER ( ONE = 1.0D0)
C
      CALL QENTER('CC3_CY1_JK')
C
C----------------------------------
C     Calculate (1):
C        =  ( W^{JK}_{bcai} - W^{JK}_{baci} ) * L_{bJcK} 
C                T(bcai)                        I(bc)
C----------------------------------
C
      J = IJ
      K = IK
      ISYMJ = ISYMIJ
      ISYMK = ISYMIK
      ISYMJK  = MULD2H(ISYMJ,ISYMK)
      ISYMAI = ISOMG1
      ISYMBC  = MULD2H(ISOMG1,ISWMAT)
C
      KBC    = 1
      KWBACI = KBC     + NMATAB(ISYMBC)
      KEND1  = KWBACI  + NMAABCI(ISWMAT)
      LWRK1  = LWORK   - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space  in CC3_CY1_JK (1)')
      END IF
C
      CALL DZERO(WORK(KBC),NMATAB(ISYMBC))
      CALL DZERO(WORK(KWBACI),NMAABCI(ISWMAT))
C
C    I(bc) =  L_{bJcK}
C
      CALL SORT_L2_IJ_AB(J,ISYMJ,K,ISYMK,WORK(KBC),XIAJB,ISYINT)
C
C      T(bcai) = ( W^{JK}_{bcai} - W^{JK}_{baci} )  
C
      CALL DAXPY(NMAABCI(ISWMAT),-ONE,WMAT,1,WORK(KWBACI),1)
      CALL DCOPY(NMAABCI(ISWMAT),WMAT,1,TMAT,1)
C 

      CALL FACBI(TMAT,WORK(KWBACI),ISWMAT)
C
C        =  ( W^{JK}_{bcai} - W^{JK}_{baci} ) * L_{bJcK} 
C                T(bcai)                        I(bc)
C
      DO ISYMI = 1,NSYM
            DO I = 1,NRHF(ISYMI)
               ISYMBCA = MULD2H(ISWMAT,ISYMI)
               ISYMA = MULD2H(ISYMAI,ISYMI)
               ISYMBC = MULD2H(ISYMBCA,ISYMA)
C 
               KOFF1 = IMAABCI(ISYMBCA,ISYMI)
     *               + NMAABC(ISYMBCA)*(I-1)
     *               + IMAABC(ISYMBC,ISYMA)
     *               + 1
               KOFF2 = KBC
C
               KOFF3 = IT1AM(ISYMA,ISYMI)
     *               + NVIR(ISYMA)*(I-1)
     *               + 1

               NTOTBC = MAX(1,NMATAB(ISYMBC))
C 
               CALL DGEMV('T',NMATAB(ISYMBC),NVIR(ISYMA),-ONE,
     *                    TMAT(KOFF1),NTOTBC,WORK(KOFF2),1,
     *                    ONE,OMEGA1(KOFF3),1)
C
         END DO
      END DO
C
C  (2):
C        +  ( W^{IJ}_{abck} - W^{IJ}_{cbak} ) * L_{ckbJ} 
C                T(abck)                         I(bck)
C
      I = IJ
      J = IK
      ISYMI    = ISYMIJ
      ISYMJ    = ISYMIK
      ISYMAI   = ISOMG1
      ISYMA    = MULD2H(ISYMAI,ISYMI)
      ISYMBCK  = MULD2H(ISYINT,ISYMJ)
C
      KBCK    = 1
      KWCBAK  = KBCK     + NMAABI(ISYMBCK)
      KEND1   = KWCBAK   + NMAABCI(ISWMAT)
      LWRK1   = LWORK    - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in CC3_CY1_JK (2)')
      END IF
C
      CALL DZERO(WORK(KBCK),NMAABI(ISYMBCK))
      CALL DZERO(WORK(KWCBAK),NMAABCI(ISWMAT))
C
C    I(bck) = L_{ckbJ} 
C
      CALL SORT_L2_J_BAI(J,ISYMJ,WORK(KBCK),XIAJB,ISYINT)
C
C     T(abck) =  W^{IJ}_{abck} - W^{IJ}_{cbak}  
C
      CALL DAXPY(NMAABCI(ISWMAT),-ONE,WMAT,1,WORK(KWCBAK),1)
      CALL DCOPY(NMAABCI(ISWMAT),WMAT,1,TMAT,1)
C 

      CALL FCBAI(TMAT,WORK(KWCBAK),ISWMAT)
C 
C        +  ( W^{IJ}_{abck} - W^{IJ}_{cbak} ) * L_{ckbJ} 
C     omega1(aI) = omega1(aI) +   T(abck)     *  I(bck)
C
      DO ISYMK = 1,NSYM
         DO ISYMC = 1,NSYM
            ISYMCK = MULD2H(ISYMC,ISYMK)
            ISYMB  = MULD2H(ISYMBCK,ISYMCK) 
            ISYMBC = MULD2H(ISYMB,ISYMC)
            DO K = 1,NRHF(ISYMK)
               DO C = 1,NVIR(ISYMC)
                  ISYMABC = MULD2H(ISWMAT,ISYMK)
                  ISYMAB = MULD2H(ISYMABC,ISYMC)
C 
                  KOFF1 = IMAABCI(ISYMABC,ISYMK)
     *                  + NMAABC(ISYMABC)*(K-1)
     *                  + IMAABC(ISYMAB,ISYMC)
     *                  + NMATAB(ISYMAB)*(C-1)
     *                  + IMATAB(ISYMA,ISYMB)
     *                  + 1
                  KOFF2 = KBCK
     *                  + IMAABI(ISYMBC,ISYMK)
     *                  + NMATAB(ISYMBC)*(K-1)
     *                  + IMATAB(ISYMB,ISYMC)
     *                  + NVIR(ISYMB)*(C-1)
C
                  KOFF3 = IT1AM(ISYMA,ISYMI)
     *                  + NVIR(ISYMA)*(I-1)
     *                  + 1

                  NTOTA = MAX(1,NVIR(ISYMA))
C 
                  CALL DGEMV('N',NVIR(ISYMA),NVIR(ISYMB),-ONE,
     *                       TMAT(KOFF1),NTOTA,WORK(KOFF2),1,
     *                       ONE,OMEGA1(KOFF3),1)
C
               END DO
            END DO
         END DO
      END DO

C
C  (3):
C
C  Omega1(ai) = Omega1(ai)  
C        +  ( W^{KI}_{cabj} - W^{KI}_{acbj} ) * L_{bjcK} 
C
C        +           T(acbj)                  * I(cbj)
C
C
      K = IJ
      I = IK
      ISYMK    = ISYMIJ
      ISYMI    = ISYMIK
      ISYMAI   = ISOMG1
      ISYMA    = MULD2H(ISYMAI,ISYMI)
      ISYMCBJ  = MULD2H(ISYINT,ISYMK)
C
      KCBJ    = 1
      KEND1   = KCBJ     + NMAABI(ISYMCBJ)
      LWRK1   = LWORK    - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in CC3_CY1_JK (3)')
      END IF
C
      CALL DZERO(WORK(KCBJ),NMAABI(ISYMCBJ))
C
C     I(cbj) = L_{bjcK}
C
      CALL SORT_L2_J_BAI(K,ISYMK,WORK(KCBJ),XIAJB,ISYINT)
C
C     T(acbj) =  W^{KI}_{cabj} - W^{KI}_{acbj}
C
      CALL DZERO(TMAT,NMAABCI(ISWMAT))
      CALL FBACI(TMAT,WMAT,ISWMAT)
      CALL DAXPY(NMAABCI(ISWMAT),-ONE,WMAT,1,TMAT,1)
C
C        +  ( W^{KI}_{cabj} - W^{KI}_{acbj} ) * L_{bjcK} 
C
C  Omega1(ai) = Omega1(ai) +    T(acbj)       * I(cbj)
C
      DO ISYMJ = 1,NSYM
         DO ISYMB = 1,NSYM
            ISYMBJ = MULD2H(ISYMB,ISYMJ)
            ISYMC  = MULD2H(ISYMCBJ,ISYMBJ)
            DO J = 1,NRHF(ISYMJ)
               DO B = 1,NVIR(ISYMB)
                  ISYMACB = MULD2H(ISWMAT,ISYMJ)
                  ISYMAC  = MULD2H(ISYMACB,ISYMB) 
                  ISYMCB  = MULD2H(ISYMCBJ,ISYMJ)
C
                  KOFF1 = IMAABCI(ISYMACB,ISYMJ)
     *                  + NMAABC(ISYMACB)*(J-1)
     *                  + IMAABC(ISYMAC,ISYMB)
     *                  + NMATAB(ISYMAC)*(B-1)
     *                  + IMATAB(ISYMA,ISYMC)
     *                  + 1
                  KOFF2 = KCBJ
     *                  + IMAABI(ISYMCB,ISYMJ)
     *                  + NMATAB(ISYMCB)*(J-1)
     *                  + IMATAB(ISYMC,ISYMB)
     *                  + NVIR(ISYMC)*(B-1)
C
                  KOFF3 = IT1AM(ISYMA,ISYMI)
     *                  + NVIR(ISYMA)*(I-1)
     *                  + 1

                  NTOTA = MAX(1,NVIR(ISYMA))
C
                  CALL DGEMV('N',NVIR(ISYMA),NVIR(ISYMC),-ONE,
     *                       TMAT(KOFF1),NTOTA,WORK(KOFF2),1,
     *                       ONE,OMEGA1(KOFF3),1)
C
               END DO
            END DO
         END DO
      END DO
C
      CALL QEXIT('CC3_CY1_JK') 
C
      RETURN
      END
C  /* Deck CC3_CY2F_JK */
      SUBROUTINE CC3_CY2F_JK(XI2,ISYRES,WMAT,ISWMAT,TMAT,
     *                   FOCK0,ISYINT,WORK,LWORK,
     *                   ISYMIJ,IJ,ISYMIK,IK)
C
C  P. Jorgensen and F. Pawlowski.         Spring 2002
C
      IMPLICIT NONE
C
C  Calculate Fock part to XI2 
C
C
C  General symmetry: ISWMAT is symmetry of WMAT and TMAT intermdiates.  
C                       (exclusive isymj,isymk) 
C                       ISYINT is symmetry of FOCK0 
C                       ISYRES is symmetry of XI2 and XI2EFF 
C 
C        Omega2(aibj) = Omega2(aibj) +  
C
C
C           sum_{ck} ( t^{abc}_{ijk}  -  t^{cba}_{ijk} ) * F_{kc}
C  (1):
C        =  ( W^{JK}_{bcai} - W^{JK}_{baci} ) * F_{Kc}
C                 w(cbai)                     *   f(c)  = O^J(bai)
C  (2):
C        +  ( W^{IJ}_{abck} - W^{IJ}_{cbak} ) * F_{kc}
C                  w(abck)                    *  f(ck)  = O^IJ(ab)
C  (3):
C        +  ( W^{KI}_{cabj} - W^{KI}_{acbj} ) * F_{Kc}
C                  w(cabj)                    *   f(c)  = O^I(abj)


      INTEGER ISYRES, ISWMAT, ISYINT , LWORK
      INTEGER ISYMIJ, IJ, ISYMIK, IK, INDEX
      INTEGER KFMAT,KEND0,LEND0,ISYMC,ISYMK,KOFF1,KOFF2 
      INTEGER ISYMJ,ISYMBAI,KC,KBAI,KWBACI,KEND1,LWRK1
      INTEGER KOFF,ISYMI,ISYMBA,ISYMCBA,ISYMA,ISYMCB,ISYMB
      INTEGER KOFF3,NTOTC,ISYMAB,KAB,KWCBAK,ISYMABC,NTOTAB
      INTEGER ISYMABJ,KABJ,KWACBJ,ISYMCAB,ISYMCA
#if defined (SYS_CRAY)
      REAL XI2(*), WMAT(*), TMAT(*), FOCK0(*), WORK(*)
      REAL ZERO, ONE, TWO
#else
      DOUBLE PRECISION XI2(*), WMAT(*), TMAT(*), FOCK0(*), WORK(*)
      DOUBLE PRECISION ZERO, ONE, TWO
#endif

C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccsdinp.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C 
      CALL QENTER('CC3_CY2F_JK')

C       
C     work space allocations
C
C
      KFMAT    = 1
      KEND0    = KFMAT  + NT1AM(ISYINT)
      LEND0    = LWORK - KEND0
C
      IF (LEND0 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND0
         CALL QUIT('Insufficient space in CC3_CY2F_JK (0)')
      END IF
C
C
C   sort fock matrix 
C
C   FMAT(C,K) = FOCK0(K,C)
C
      DO ISYMC = 1,NSYM
         ISYMK = MULD2H(ISYMC,ISYINT)
         DO K = 1,NRHF(ISYMK)
            DO C = 1,NVIR(ISYMC)
               KOFF1 = IFCVIR(ISYMK,ISYMC) + NORB(ISYMK)*(C - 1) + K
               KOFF2 = KFMAT +  IT1AM(ISYMC,ISYMK)
     *               + NVIR(ISYMC)*(K - 1) + C -1
C
                  WORK(KOFF2) = FOCK0(KOFF1)
C
            END DO
         END DO
      END DO
C
C----------------------------------
C  (1):
C     Omega2(aibj) = Omega2(aibj)  
C        =  ( W^{JK}_{bcai} - W^{JK}_{baci} ) * F_{Kc}
C                 w(cbai)                     *   f(c)  = O^J(bai)
C----------------------------------
C
      J = IJ
      K = IK
      ISYMJ = ISYMIJ
      ISYMK = ISYMIK
      ISYMC = MULD2H(ISYINT,ISYMK)
      ISYMBAI = MULD2H(ISWMAT,ISYMC)
C
      KC      = KEND0
      KBAI    = KC + NVIR(ISYMC)
      KWBACI  = KBAI + NMAABI(ISYMBAI)
      KEND1   = KWBACI + NMAABCI(ISWMAT)
      LWRK1   = LWORK   - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in CC3_CY2F_JK (1)')
      END IF
C
      CALL DZERO(WORK(KBAI),NMAABI(ISYMBAI))
      CALL DZERO(WORK(KWBACI),NMAABCI(ISWMAT))


C
C  w(cbai) = W^{JK}_{bcai} - W^{JK}_{baci}
C
      CALL DZERO(TMAT,NMAABCI(ISWMAT))
      CALL FBACI(TMAT,WMAT,ISWMAT)
      CALL DAXPY(NMAABCI(ISWMAT),-ONE,WMAT,1,WORK(KWBACI),1)
C
      CALL FBCAI(TMAT,WORK(KWBACI),ISWMAT)
C
C f(c) =  F_{Kc}
C
      KOFF = KFMAT + IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K - 1) 
      CALL DCOPY(NVIR(ISYMC),WORK(KOFF),1,WORK(KC),1)
C
C      O^J(bai) = w(cbai) * f(c) 
C
      DO ISYMI = 1,NSYM
         ISYMBA = MULD2H(ISYMBAI,ISYMI)
         ISYMCBA = MULD2H(ISWMAT,ISYMI)
         DO ISYMA = 1,NSYM
            ISYMCB = MULD2H(ISYMCBA,ISYMA)
            ISYMB = MULD2H(ISYMBA,ISYMA)
            DO I = 1,NRHF(ISYMI)
               DO A = 1,NVIR(ISYMA)
                  KOFF1  =  IMAABCI(ISYMCBA,ISYMI)
     *               + NMAABC(ISYMCBA)*(I-1)
     *               + IMAABC(ISYMCB,ISYMA)
     *               + NMATAB(ISYMCB)*(A-1) 
     *               + IMATAB(ISYMC,ISYMB) + 1 
C
                  KOFF3  = KBAI + IMAABI(ISYMBA,ISYMI)
     *                   +  NMATAB(ISYMBA)*(I-1)
     *                   +  IMATAB(ISYMB,ISYMA)
     *                   +  NVIR(ISYMB)*(A-1)
C 
                  NTOTC = MAX(NVIR(ISYMC),1)
C
                  CALL DGEMV('T',NVIR(ISYMC),NVIR(ISYMB),-ONE,
     *                       TMAT(KOFF1), NTOTC,WORK(KC),1,ONE,
     *                       WORK(KOFF3),1)
C
               END DO
            END DO
         END DO
      END DO
C
C Distribute O^J(bai) in result vector
C
      CALL CC3_RI_ABJ(XI2,WORK(KBAI),ISYMJ,J,ISYRES)
C
C----------------------------------
C  (2):
C     Omega2(aibj) = Omega2(aibj)  
C        +  ( W^{IJ}_{abck} - W^{IJ}_{cbak} ) * F_{kc}
C                  w(abck)                    *  f(ck)  = O^IJ(ab)
C
      I = IJ
      J = IK 
      ISYMI = ISYMIJ
      ISYMJ = ISYMIK
C
C
      ISYMAB = MULD2H(ISWMAT,ISYINT)
C
      KAB    = KEND0 
      KWCBAK  = KAB + NMATAB(ISYMAB)
      KEND1   = KWCBAK + NMAABCI(ISWMAT) 
      LWRK1   = LWORK   - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in CC3_CY2F_JK (2)')
      END IF
C
      CALL DZERO(WORK(KAB),NMATAB(ISYMAB))
      CALL DZERO(WORK(KWCBAK),NMAABCI(ISWMAT))
C
C  w(abck) = W^{IJ}_{abck} - W^{IJ}_{cbak}
C
      CALL DCOPY(NMAABCI(ISWMAT),WMAT,1,TMAT,1)
      CALL DAXPY(NMAABCI(ISWMAT),-ONE,WMAT,1,WORK(KWCBAK),1)
C
      CALL FCBAI(TMAT,WORK(KWCBAK),ISWMAT)
C
C  O^IJ(ab) =  w(abck)  *  f(ck) 
C
      DO ISYMK = 1,NSYM
         ISYMABC = MULD2H(ISWMAT,ISYMK)
         ISYMC   = MULD2H(ISYINT,ISYMK)
         ISYMAB  = MULD2H(ISYMABC,ISYMC)
         DO K = 1,NRHF(ISYMK)
            KOFF1  =   IMAABCI(ISYMABC,ISYMK)
     *               + NMAABC(ISYMABC)*(K-1)
     *               + IMAABC(ISYMAB,ISYMC) + 1
C
            KOFF2  = KFMAT + IT1AM(ISYMC,ISYMK)
     *                     + NVIR(ISYMC)*(K-1)
C
            KOFF3  = KAB 
            NTOTAB = MAX(NMATAB(ISYMAB),1)
C
            CALL DGEMV('N',NMATAB(ISYMAB),NVIR(ISYMC),-ONE,
     *                 TMAT(KOFF1), NTOTAB,WORK(KOFF2),1,ONE,
     *                 WORK(KOFF3),1)
         END DO
      END DO
C
C Distribute O^J(ab) in result vector
C
      CALL CC3_RIJ_AB(XI2,WORK(KAB),ISYMI,I,ISYMJ,J,ISYRES)
C
C
C----------------------------------
C  (3):
C     Omega2(aibj) = Omega2(aibj)  
C        +  ( W^{KI}_{cabj} - W^{KI}_{acbj} ) * F_{Kc}
C                  w(cabj)                    *   f(c)  = O^I(abj)
C----------------------------------
C
      K = IJ
      I = IK
      ISYMK = ISYMIJ
      ISYMI = ISYMIK
      ISYMC = MULD2H(ISYINT,ISYMK)
      ISYMABJ = MULD2H(ISWMAT,ISYMC)
C
      KC      = KEND0
      KABJ    = KC + NVIR(ISYMC)
      KWACBJ  = KABJ + NMAABI(ISYMABJ)
      KEND1   = KWACBJ + NMAABCI(ISWMAT)
      LWRK1   = LWORK   - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in CC3_CY2F_JK (3)')
      END IF
C
      CALL DZERO(WORK(KABJ),NMAABI(ISYMABJ))
      CALL DZERO(WORK(KWACBJ),NMAABCI(ISWMAT))


C
C w(cabj)  = W^{KI}_{cabj} - W^{KI}_{acbj}
C
      CALL DCOPY(NMAABCI(ISWMAT),WMAT,1,TMAT,1)
      CALL DAXPY(NMAABCI(ISWMAT),-ONE,WMAT,1,WORK(KWACBJ),1)
C
      CALL FBACI(TMAT,WORK(KWACBJ),ISWMAT)
C
C f(c) =  F_{Kc}
C
      KOFF = KFMAT + IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K - 1)
      CALL DCOPY(NVIR(ISYMC),WORK(KOFF),1,WORK(KC),1)
C
C      O^I(abj) = w(cabj) * f(c)
C
      DO ISYMJ = 1,NSYM
         ISYMAB = MULD2H(ISYMABJ,ISYMJ)
         ISYMCAB = MULD2H(ISWMAT,ISYMJ)
         DO ISYMB = 1,NSYM
            ISYMCA = MULD2H(ISYMCAB,ISYMB)
            ISYMA = MULD2H(ISYMCA,ISYMC)
            DO J = 1,NRHF(ISYMJ)
               DO B = 1,NVIR(ISYMB)
                  KOFF1  =  IMAABCI(ISYMCAB,ISYMJ)
     *               + NMAABC(ISYMCAB)*(J-1)
     *               + IMAABC(ISYMCA,ISYMB)
     *               + NMATAB(ISYMCA)*(B-1)
     *               + IMATAB(ISYMC,ISYMA) + 1
C
                  KOFF3  = KABJ + IMAABI(ISYMAB,ISYMJ)
     *                   +  NMATAB(ISYMAB)*(J-1)
     *                   +  IMATAB(ISYMA,ISYMB)
     *                   +  NVIR(ISYMA)*(B-1)
C
                  NTOTC = MAX(NVIR(ISYMC),1)
C
                  CALL DGEMV('T',NVIR(ISYMC),NVIR(ISYMA),-ONE,
     *                       TMAT(KOFF1), NTOTC,WORK(KC),1,ONE,
     *                       WORK(KOFF3),1)
C
               END DO
            END DO
         END DO
      END DO
C
C Distribute O^J(bai) in result vector
C
      CALL CC3_RI_ABJ(XI2,WORK(KABJ),ISYMI,I,ISYRES)
C
      CALL QEXIT('CC3_CY2F_JK') 
C
C
      RETURN
      END
C  /* Deck cc3_cy2o_jk */
      SUBROUTINE CC3_CY2O_JK(XI2,ISYRES,WMAT,ISWMAT,TMAT,
     *                      T3BOG2,ISYINT,WORK,LWORK,
     *                      ISYMIJ,IJ,ISYMIK,IK)
C
C
C     General symmetry: ISWMAT is symmetry of WMAT and TMAT intermediates.
C                       (exclusive isymib,isymid)
C                       ISYINT is symmetry of integrals in TROCC and TROCC1.
C                       ISYRES = ISWMAT*ISYINT*ISYMIJ*ISYMIK
C
C   XI2(aibj) =  XI2(aibj)
C
C     + sum_{ckl} (2W^{abc}_{ikl}-W^{abc}_{lki}-W^{abc}_{ilk}) * (lc|kj)
C (1):
C     = (2W^{KL}_{bcai}-W^{KL}_{baci}-W^{KL}_{cbai}) * (Lc|Kj)
C                  w(cbai)                           *  g(cj)  = O(bjia)
C (2):
C     + (2W^{IK}_{abcl}-W^{IK}_{cbal}-W^{IK}_{acbl}) * (lc|Kj)
C                  w(abcl)                           *  g(jcl)  = O^I(abj)
C (3):
C    +  (2W^{LI}_{cabk}-W^{LI}_{acbk}-W^{LI}_{back}) * (Lc|kj)
C                  w(abck)                           *  g(jck)  = O^I(abj)
C
C                   (lc|kj) = g(cklj)
C
      IMPLICIT NONE
C
      INTEGER ISYRES, ISWMAT, ISYINT, LWORK
      INTEGER ISYMIJ, IJ, ISYMIK, IK 
      INTEGER ISYMK,ISYML,ISYMKL,ISYMCJ,ISYMI,ISYMCBA,ISYMA,ISYMIA 
      INTEGER ISYMBJ,ISYMBJI,ISYMCB,ISYMJ,ISYMC,ISYMB,ISYMCL,ISYMAB
      INTEGER ISYMABJ,ISYMJCK,ISYMCK
      INTEGER KCJ,KWBCAI,KRBJIA,KEND1,LWRK1,KOFF1,KOFF2,KOFF3
      INTEGER KJCL,KWCBAL,KRABJ,KJCK,KW
      INTEGER NTOTB,NTOTC,NTOTAB,NTOTJ,IOPT
      INTEGER ISYMJCL
C
#if defined (SYS_CRAY) 
      REAL XI2(*), WMAT(*), TMAT(*), T3BOG2(*), WORK(*) 
      REAL TWO, ONE
#else
      DOUBLE PRECISION XI2(*), WMAT(*), TMAT(*), T3BOG2(*), WORK(*) 
      DOUBLE PRECISION TWO, ONE 
#endif
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccsdinp.h>
C
      PARAMETER ( ONE = 1.0D0, TWO = 2.0D0)
C
      CALL QENTER('CC3_CY2O_JK')
C
C-------------------------
C     First occupied term.
C-------------------------
C (1):
C             XI2(aibj) =  XI2(aibj) 
C
C     = (2W^{KL}_{bcai}-W^{KL}_{baci}-W^{KL}_{cbai}) * (Lc|Kj)
C                  w(cbai)                           *  g(cj)  = O(bjia)
C
      K = IJ
      L = IK
      ISYMK = ISYMIJ
      ISYML = ISYMIK
C
      ISYMKL = MULD2H(ISYMK,ISYML)
      ISYMCJ = MULD2H(ISYINT,ISYMKL)
C
      KCJ      = 1
      KWBCAI   = KCJ    + NT1AM(ISYMCJ)
      KRBJIA   = KWBCAI + NMAABCI(ISWMAT)
      KEND1    = KRBJIA + NT2SQ(ISYRES) 
      LWRK1   = LWORK   - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in CC3_CY2O_JK (1)')
      END IF
C
      CALL DZERO(WORK(KCJ),NT1AM(ISYMCJ))
      CALL DZERO(WORK(KRBJIA),NT2SQ(ISYRES))


C
C  w(cbai) = 2W^{KL}_{bcai}-W^{KL}_{baci}-W^{KL}_{cbai}
C
      CALL DZERO(TMAT,NMAABCI(ISWMAT))
      CALL DAXPY(NMAABCI(ISWMAT),-ONE,WMAT,1,TMAT,1)
      CALL DZERO(WORK(KWBCAI),NMAABCI(ISWMAT))
      CALL DAXPY(NMAABCI(ISWMAT),-ONE,WMAT,1,WORK(KWBCAI),1)
      CALL FBCAI(TMAT,WORK(KWBCAI),ISWMAT)
      CALL DZERO(WORK(KWBCAI),NMAABCI(ISWMAT))
      CALL DAXPY(NMAABCI(ISWMAT),TWO,WMAT,1,WORK(KWBCAI),1)
      CALL FBACI(TMAT,WORK(KWBCAI),ISWMAT)
C
C Sort (Lc|Kj) = T3BOG2(C,k,l,j) as  g^LK(cj)
C
      CALL SORT_INT_AJ_IK(WORK(KCJ),T3BOG2,ISYINT,ISYML,L,ISYMK,K)
C
C   O(bjia) =  w(cbai) * g(cj) 
C
      DO ISYMI = 1,NSYM
         ISYMCBA = MULD2H(ISWMAT,ISYMI)
         DO ISYMA = 1,NSYM
            ISYMIA = MULD2H(ISYMI,ISYMA)
            ISYMBJ = MULD2H(ISYRES,ISYMIA)
            ISYMBJI = MULD2H(ISYMBJ,ISYMI)
            ISYMCB = MULD2H(ISYMCBA,ISYMA)
            DO ISYMJ = 1,NSYM
               ISYMC = MULD2H(ISYMCJ,ISYMJ)
               ISYMB = MULD2H(ISYMBJ,ISYMJ)
               DO I = 1,NRHF(ISYMI)
                  DO A = 1,NVIR(ISYMA)
C
                     KOFF1  =  IMAABCI(ISYMCBA,ISYMI)
     *                       + NMAABC(ISYMCBA)*(I-1)
     *                       + IMAABC(ISYMCB,ISYMA)
     *                       + NMATAB(ISYMCB)*(A-1)
     *                       + IMATAB(ISYMC,ISYMB) + 1
C
                     KOFF2  = KCJ + IT1AM(ISYMC,ISYMJ)
C
                     KOFF3  = KRBJIA + IT2SP(ISYMBJI,ISYMA) 
     *                       + NCKI(ISYMBJI)*(A-1)
     *                       + ISAIK(ISYMBJ,ISYMI) 
     *                       + NT1AM(ISYMBJ)*(I-1)
     *                       + IT1AM(ISYMB,ISYMJ)
C
                     NTOTB = MAX(NVIR(ISYMB),1)
                     NTOTC = MAX(NVIR(ISYMC),1)
C
C
                     CALL DGEMM('T','N',NVIR(ISYMB),NRHF(ISYMJ),
     *                          NVIR(ISYMC),ONE,TMAT(KOFF1),NTOTC,
     *                          WORK(KOFF2),NTOTC,
     *                          ONE,WORK(KOFF3),NTOTB)
C
                  END DO
               END DO
            END DO
         END DO
      END DO
C
C Distribute O(bjia) in result vector
C
      CALL CC3_RBJIA(XI2,ISYRES,WORK(KRBJIA))
C
C--------------------------
C     Second occupied term.
C--------------------------
C
C (2):
C
C   XI2(aibj) =  XI2(aibj)
C
C     + (2W^{IK}_{abcl}-W^{IK}_{cbal}-W^{IK}_{acbl}) * (lc|Kj)
C                  w(abcl)                           *  g(jcl)  = O^I(abj)
C 
      I = IJ
      K = IK
      ISYMI = ISYMIJ
      ISYMK = ISYMIK
C
      ISYMABJ = MULD2H(ISYRES,ISYMI)
      ISYMJCL = MULD2H(ISYINT,ISYMK)
C
      KJCL     = 1
      KWCBAL   = KJCL   + NCKI(ISYMJCL)
      KRABJ    = KWCBAL + NMAABCI(ISWMAT)
      KEND1    = KRABJ  + NMAABI(ISYMABJ)
      LWRK1   = LWORK   - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in CC3_CY2O_JK (2)')
      END IF
C
      CALL DZERO(WORK(KJCL),NCKI(ISYMJCL))
      CALL DZERO(WORK(KRABJ),NMAABI(ISYMABJ))
C
C  w(abcl  = 2W^{IK}_{abcl}-W^{IK}_{cbal}-W^{IK}_{acbl}
C
      CALL DZERO(TMAT,NMAABCI(ISWMAT))
      CALL DAXPY(NMAABCI(ISWMAT),TWO,WMAT,1,TMAT,1)
      CALL DZERO(WORK(KWCBAL),NMAABCI(ISWMAT))
      CALL DAXPY(NMAABCI(ISWMAT),-ONE,WMAT,1,WORK(KWCBAL),1)
      CALL FCBAI(TMAT,WORK(KWCBAL),ISWMAT)
      CALL DZERO(WORK(KWCBAL),NMAABCI(ISWMAT))
      CALL DAXPY(NMAABCI(ISWMAT),-ONE,WMAT,1,WORK(KWCBAL),1)
      CALL FACBI(TMAT,WORK(KWCBAL),ISWMAT)
C
      IF (NSYM.GT.1) THEN
         IOPT = 1
         CALL DZERO(WORK(KWCBAL),NMAABCI(ISWMAT))
         CALL FAB_CI(WORK(KWCBAL),TMAT,ISWMAT,IOPT)
         CALL DCOPY(NMAABCI(ISWMAT),WORK(KWCBAL),1,TMAT,1)
      END IF
C 
C Sort (lc|Kj) = T3BOG2(c,k,l,j) as  g^K(jcl) 
C
      CALL SORT_INT_JAI_K(WORK(KJCL),T3BOG2,ISYINT,ISYMK,K)
C
C O^I(abj) = w(abcl) *  g(jcl) 
C
      DO ISYMCL = 1,NSYM
         ISYMAB = MULD2H(ISWMAT,ISYMCL)
         ISYMJ  = MULD2H(ISYMJCL,ISYMCL)
C
         KOFF1 = 1
     *         + IMAAB_CI(ISYMAB,ISYMCL)
         KOFF2 = KJCL
     *         + IMAIAJ(ISYMJ,ISYMCL)
         KOFF3 = KRABJ
     *         + IMAABI(ISYMAB,ISYMJ) 
C
         NTOTAB = MAX(NMATAB(ISYMAB),1)
         NTOTJ  = MAX(NRHF(ISYMJ),1)
C
         CALL DGEMM('N','T',NMATAB(ISYMAB),NRHF(ISYMJ),NT1AM(ISYMCL),
     *              ONE,TMAT(KOFF1),NTOTAB,WORK(KOFF2),NTOTJ,
     *              ONE,WORK(KOFF3),NTOTAB)
C
      END DO ! ISYMN
C
C Distribute O^I(abj) in result vector
C
      CALL CC3_RI_ABJ(XI2,WORK(KRABJ),ISYMI,I,ISYRES)
C
C-------------------------
C     Third occupied term.
C-------------------------
C (3):
C             XI2(aibj) =  XI2(aibj) 
C
C    +  (2W^{LI}_{cabk}-W^{LI}_{acbk}-W^{LI}_{back}) * (Lc|kj)
C                  w(abck)                           *  g(jck)  = O^I(abj)
C
C----------------------------------
C
      L = IJ
      I = IK
      ISYML = ISYMIJ
      ISYMI = ISYMIK
C
      ISYMABJ = MULD2H(ISYRES,ISYMI)
      ISYMJCK = MULD2H(ISYINT,ISYML)
C
      KJCK     = 1
      KW       = KJCK   + NCKI(ISYMJCK)
      KRABJ    = KW     + NMAABCI(ISWMAT)
      KEND1    = KRABJ  + NMAABI(ISYMABJ)
      LWRK1   = LWORK   - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in CC3_CY2O_JK (3)')
      END IF
C
      CALL DZERO(WORK(KJCK),NCKI(ISYMJCK))
      CALL DZERO(WORK(KRABJ),NMAABI(ISYMABJ))
C
C  w(abck) = 2W^{LI}_{cabk}-W^{LI}_{acbk}-W^{LI}_{back}
C
      CALL DZERO(TMAT,NMAABCI(ISWMAT))
      CALL DZERO(WORK(KW),NMAABCI(ISWMAT))
      CALL DAXPY(NMAABCI(ISWMAT),TWO,WMAT,1,WORK(KW),1)
      CALL FCABI(TMAT,WORK(KW),ISWMAT)
      CALL DZERO(WORK(KW),NMAABCI(ISWMAT))
      CALL DAXPY(NMAABCI(ISWMAT),-ONE,WMAT,1,WORK(KW),1)
      CALL FACBI(TMAT,WORK(KW),ISWMAT)
      CALL DZERO(WORK(KW),NMAABCI(ISWMAT))
      CALL DAXPY(NMAABCI(ISWMAT),-ONE,WMAT,1,WORK(KW),1)
      CALL FBACI(TMAT,WORK(KW),ISWMAT)
C
      IF (NSYM.GT.1) THEN
         IOPT = 1
         CALL DZERO(WORK(KW),NMAABCI(ISWMAT))
         CALL FAB_CI(WORK(KW),TMAT,ISWMAT,IOPT)
         CALL DCOPY(NMAABCI(ISWMAT),WORK(KW),1,TMAT,1)
      END IF
C
C Sort (lc|Kj) = T3BOG2(c,k,l,j) as  g^L(jck)
C
      CALL SORT_INT_JAK_I(WORK(KJCK),T3BOG2,ISYINT,ISYML,L)
C
C O^I(abj) =  w(abck) * g(jck) 
C
      DO ISYMCK = 1,NSYM
         ISYMAB = MULD2H(ISWMAT,ISYMCK)
         ISYMJ  = MULD2H(ISYMJCK,ISYMCK)
C
         KOFF1 = 1
     *         + IMAAB_CI(ISYMAB,ISYMCK)
         KOFF2 = KJCK
     *         + IMAIAJ(ISYMJ,ISYMCK)
         KOFF3 = KRABJ
     *         + IMAABI(ISYMAB,ISYMJ)
C
         NTOTAB = MAX(NMATAB(ISYMAB),1)
         NTOTJ  = MAX(NRHF(ISYMJ),1)
C
         CALL DGEMM('N','T',NMATAB(ISYMAB),NRHF(ISYMJ),NT1AM(ISYMCK),
     *              ONE,TMAT(KOFF1),NTOTAB,WORK(KOFF2),NTOTJ,
     *              ONE,WORK(KOFF3),NTOTAB)
C
      END DO
C
C Distribute O^I(abj) in result vector
C
      CALL CC3_RI_ABJ(XI2,WORK(KRABJ),ISYMI,I,ISYRES)
C
      CALL QEXIT('CC3_CY2O_JK')
C
      RETURN
      END
C  /* Deck cc3_cy2v_jk */
      SUBROUTINE CC3_CY2V_JK(XI2,ISYRES,WMAT,ISWMAT,TMAT,XGACDK,
     *                    ISYINT,WORK,LWORK,
     *                    ISYMIJ,IJ,ISYMIK,IK)
C
C    
C General symmetry: ISWMAT is symmetry of WMAT and TMAT intermediates.
C                   (exclusive isymib,isymid)
C                   ISYINT is symmetry of integrals in TRVIR and TROVIR1.
C                   ISYRES = ISWMAT*ISYINT*ISYMIJ*ISYMIK
C                   RBJIA intermediate result vector
C
C
C   XI2(aibj) =  XI2(aibj)
C
C     + sum_{cdk} (2t^{cbd}_{ijk}-t^{cdb}_{ijk}-t^{dbc}_{ijk}) * (ac|kd)
C (1):
C     =  (2W^{JK}_{bdci}-W^{JK}_{dbci}-W^{JK}_{bcdi}) * (ac|Kd)
C                       w(dcbi)                       * g(dca)  = O^J(abi)
C (2):
C     +  (2W^{IJ}_{cbdk}-W^{IJ}_{cdbk}-W^{IJ}_{dbck}) * (ac|kd)
C                       w(bcdk)                       * g(acdk) = O^IJ(ab)
C (3):
C    +   (2W^{KI}_{dcbj}-W^{KI}_{bcdj}-W^{KI}_{cdbj}) * (ac|Kd)
C                       w(dcbj)                       * g(dca)  = O^I(abj)
C
C                   (ac|kd) = g(acdk)

      IMPLICIT NONE
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccsdinp.h>
C
      INTEGER ISYRES, ISWMAT, ISYINT , LWORK
      INTEGER ISYMIJ, IJ, ISYMIK, IK
      INTEGER ISYMJ,ISYMK,ISYMABI,ISYMDCA,ISYMI,ISYMDCB,ISYMA,ISYMDC 
      INTEGER ISYMB,ISYMIIJJ,ISYMAB,ISYMCDK
      INTEGER KDCA,KW,KRABI,KRABIS,KEND1,LWRK1,KOFF1,KOFF2,KOFF3
      INTEGER KRAB,KINT,KRABJ
      INTEGER NTOTA,NTOTDC 
      INTEGER NTOTB,ISYMABJ
#if defined (SYS_CRAY)
      REAL XI2(*), WMAT(*), TMAT(*), XGACDK(*), WORK(*)
      REAL  ONE, TWO
#else
      DOUBLE PRECISION XI2(*), WMAT(*), TMAT(*), XGACDK(*), WORK(*) 
      DOUBLE PRECISION  ONE, TWO
#endif

C
      PARAMETER ( ONE = 1.0D0, TWO = 2.0D0)
C
C
      CALL QENTER('CC3_CY2V_JK')
C
C------------------------
C (1):
C
C   XI2(aibj) =  XI2(aibj)
C
C     =  (2W^{JK}_{bdci}-W^{JK}_{dbci}-W^{JK}_{bcdi}) * (ac|Kd)
C                       w(dcbi)                       * g(dca)  = O^J(abi)
C
C------------------------
C
      J = IJ
      K = IK
      ISYMJ = ISYMIJ
      ISYMK = ISYMIK
C
      ISYMABI = MULD2H(ISYRES,ISYMJ)
      ISYMDCA = MULD2H(ISYINT,ISYMK)
C
      KDCA     = 1
      KW       = KDCA   + NMAABC(ISYMDCA)
      KRABI    = KW     + NMAABCI(ISWMAT)
      KRABIS   = KRABI  + NMAABI(ISYMABI)
      KEND1    = KRABIS + NMAABI(ISYMABI)
      LWRK1   = LWORK   - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in CC3_CY2V_JK (1)')
      END IF
C
      CALL DZERO(WORK(KDCA),NMAABC(ISYMDCA))
      CALL DZERO(WORK(KRABI),NMAABI(ISYMABI))
C
C  w(dcbi) = 2W^{JK}_{bdci}-W^{JK}_{dbci}-W^{JK}_{bcdi}
C
      CALL DZERO(TMAT,NMAABCI(ISWMAT))
      CALL DZERO(WORK(KW),NMAABCI(ISWMAT))
      CALL DAXPY(NMAABCI(ISWMAT),TWO,WMAT,1,WORK(KW),1)
      CALL FCABI(TMAT,WORK(KW),ISWMAT)
      CALL DZERO(WORK(KW),NMAABCI(ISWMAT))
      CALL DAXPY(NMAABCI(ISWMAT),-ONE,WMAT,1,WORK(KW),1)
      CALL FACBI(TMAT,WORK(KW),ISWMAT)
      CALL DZERO(WORK(KW),NMAABCI(ISWMAT))
      CALL DAXPY(NMAABCI(ISWMAT),-ONE,WMAT,1,WORK(KW),1)
      CALL FCBAI(TMAT,WORK(KW),ISWMAT)
C
C Sort (ac|Kd) =I(acdk) as g^K(dca)
C
      CALL SORT_INT_CBA(WORK(KDCA),ISYMK,K,XGACDK,ISYINT)
C
C O^J(abi) =  w(dcbi)  * g(dca)  
C
C
      DO ISYMI = 1,NSYM
         ISYMDCB = MULD2H(ISWMAT,ISYMI)
         DO ISYMA = 1,NSYM
            ISYMDC = MULD2H(ISYMDCA,ISYMA)
            ISYMB  = MULD2H(ISYMDCB,ISYMDC)
            ISYMAB = MULD2H(ISYMA,ISYMB)
            DO I = 1,NRHF(ISYMI)
C
               KOFF1  = KDCA 
     *                  + IMAABC(ISYMDC,ISYMA) 
C
               KOFF2  =   IMAABCI(ISYMDCB,ISYMI)
     *                  + NMAABC(ISYMDCB)*(I-1)
     *                  + IMAABC(ISYMDC,ISYMB) + 1
C
               KOFF3  = KRABI 
     *                  + IMAABI(ISYMAB,ISYMI)
     *                  + NMATAB(ISYMAB)*(I-1)
     *                  + IMATAB(ISYMA,ISYMB)
C
               NTOTDC = MAX(NMATAB(ISYMDC),1)
               NTOTA  = MAX(NVIR(ISYMA),1)
C
               CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),
     *                    NMATAB(ISYMDC),
     *                    -ONE,WORK(KOFF1),NTOTDC,TMAT(KOFF2),NTOTDC,
     *                    ONE,WORK(KOFF3),NTOTA)
C

C
            END DO
         END DO
      END DO
C
C Distribute O^J(abi) in result vector
C
C  sort O^J(abi) as O^J(bai)
C
      CALL CC3_MABI(WORK(KRABI),WORK(KRABIS),ISYMJ,J,ISYRES) 
C
C  put O^J(bai) in result vector
C
      CALL CC3_RI_ABJ(XI2,WORK(KRABIS),ISYMJ,J,ISYRES)
C
C (2):
C
C   XI2(aibj) =  XI2(aibj)
C
C     +  (2W^{IJ}_{cbdk}-W^{IJ}_{cdbk}-W^{IJ}_{dbck}) * (ac|kd)
C                       w(bcdk)                       * g(acdk) = O^IJ(ab)
C
C-----------------------------------
C
      I = IJ
      J = IK
      ISYMI = ISYMIJ
      ISYMJ = ISYMIK
C
      ISYMIIJJ  = MULD2H(ISYMI,ISYMJ)
      ISYMAB  = MULD2H(ISYRES,ISYMIIJJ)
C
      KW       = 1 
      KRAB     = KW     + NMAABCI(ISWMAT)
      KINT     = KRAB   + NMATAB(ISYMAB)
      KEND1    = KINT   + NMAABCI(ISYINT)
      LWRK1    = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in CC3_CY2V_JK (2)')
      END IF
C
      CALL DZERO(WORK(KRAB),NMATAB(ISYMAB))
C
C  w(bcdk) = 2W^{IJ}_{cbdk}-W^{IJ}_{cdbk}-W^{IJ}_{dbck}
C
      CALL DZERO(TMAT,NMAABCI(ISWMAT))
      CALL DZERO(WORK(KW),NMAABCI(ISWMAT))
      CALL DAXPY(NMAABCI(ISWMAT),TWO,WMAT,1,WORK(KW),1)
      CALL FBACI(TMAT,WORK(KW),ISWMAT)
      CALL DZERO(WORK(KW),NMAABCI(ISWMAT))
      CALL DAXPY(NMAABCI(ISWMAT),-ONE,WMAT,1,WORK(KW),1)
      CALL FBCAI(TMAT,WORK(KW),ISWMAT)
      CALL DZERO(WORK(KW),NMAABCI(ISWMAT))
      CALL DAXPY(NMAABCI(ISWMAT),-ONE,WMAT,1,WORK(KW),1)
      CALL FCABI(TMAT,WORK(KW),ISWMAT)

      IF (NSYM.GT.1) THEN
C
C  Sort integrals from XGADCK(d,c,b,i) to KDCBI(d,cbi)
C
         CALL DZERO(WORK(KINT),NMAABCI(ISYINT))
         CALL FA_BCI(WORK(KINT),XGACDK,ISYINT,2)
C
C Sort wmat
C
         CALL DZERO(WORK(KW),NMAABCI(ISWMAT))
         CALL FA_BCI(WORK(KW),TMAT,ISWMAT,2)
         CALL DCOPY(NMAABCI(ISWMAT),WORK(KW),1,TMAT,1)
      ELSE
         CALL DCOPY(NMAABCI(ISYINT),XGACDK,1,WORK(KINT),1)
      END IF
C  
C   O^IJ(ab) =  w(bcdk) * g(acdk) 
C
      DO ISYMCDK = 1,NSYM
         ISYMA   = MULD2H(ISYMCDK,ISYINT)
         ISYMB   = MULD2H(ISYMCDK,ISWMAT)
         KOFF1   = KINT  + IMAAOBCI(ISYMA,ISYMCDK)
         KOFF2   = 1     + IMAAOBCI(ISYMB,ISYMCDK)
         KOFF3   = KRAB   + IMATAB(ISYMA,ISYMB)
C
         NTOTB  = MAX(1,NVIR(ISYMB))
         NTOTA  = MAX(1,NVIR(ISYMA))
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),
     *                    NMAABI(ISYMCDK),-ONE,WORK(KOFF1),NTOTA,
     *                    TMAT(KOFF2),NTOTB,
     *                    ONE,WORK(KOFF3),NTOTA)
      END DO
C
C Distribute O^IJ(ab) in result vector
C
      CALL CC3_RIJ_AB(XI2,WORK(KRAB),ISYMI,I,ISYMJ,J,ISYRES)
C
C (3):
C
C   XI2(aibj) =  XI2(aibj)
C
C    +   (2W^{KI}_{dcbj}-W^{KI}_{bcdj}-W^{KI}_{cdbj}) * (ac|Kd)
C                       w(dcbj)                       * g(dca)  = O^I(abj)
C
C-----------------------------------
C
      K = IJ
      I = IK
      ISYMK = ISYMIJ
      ISYMI = ISYMIK
C
      ISYMABJ = MULD2H(ISYRES,ISYMI)
      ISYMDCA = MULD2H(ISYINT,ISYMK)
C
      KDCA     = 1
      KW       = KDCA   + NMAABC(ISYMDCA)
      KRABJ    = KW     + NMAABCI(ISWMAT)
      KEND1    = KRABJ  + NMAABI(ISYMABJ)
      LWRK1   = LWORK   - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in CC3_CY2V_JK (3)')
      END IF
C
      CALL DZERO(WORK(KDCA),NMAABC(ISYMDCA))
      CALL DZERO(WORK(KRABJ),NMAABI(ISYMABJ))
C
C  w(dcbj) = 2W^{KI}_{dcbj}-W^{KI}_{bcdj}-W^{KI}_{cdbj}
C
      CALL DZERO(TMAT,NMAABCI(ISWMAT))
      CALL DAXPY(NMAABCI(ISWMAT),TWO,WMAT,1,TMAT,1)
      CALL DZERO(WORK(KW),NMAABCI(ISWMAT))
      CALL DAXPY(NMAABCI(ISWMAT),-ONE,WMAT,1,WORK(KW),1)
      CALL FCBAI(TMAT,WORK(KW),ISWMAT)
      CALL DZERO(WORK(KW),NMAABCI(ISWMAT))
      CALL DAXPY(NMAABCI(ISWMAT),-ONE,WMAT,1,WORK(KW),1)
      CALL FBACI(TMAT,WORK(KW),ISWMAT)
C
C
C Sort (ac|Kd) =I(acdk) as g^K(dca)
C
      CALL SORT_INT_CBA(WORK(KDCA),ISYMK,K,XGACDK,ISYINT)
C
C O^I(abj) =  w(dcbj) * g(dca) 
C
C
      DO ISYMJ = 1,NSYM
         ISYMDCB = MULD2H(ISWMAT,ISYMJ)
         DO ISYMA = 1,NSYM
            ISYMDC = MULD2H(ISYMDCA,ISYMA)
            ISYMB  = MULD2H(ISYMDCB,ISYMDC)
            ISYMAB = MULD2H(ISYMA,ISYMB)
            DO J = 1,NRHF(ISYMJ)
C
               KOFF1  = KDCA
     *                  + IMAABC(ISYMDC,ISYMA)
C
               KOFF2  =   IMAABCI(ISYMDCB,ISYMJ)
     *                  + NMAABC(ISYMDCB)*(J-1)
     *                  + IMAABC(ISYMDC,ISYMB) + 1
C
               KOFF3  = KRABJ
     *                  + IMAABI(ISYMAB,ISYMJ)
     *                  + NMATAB(ISYMAB)*(J-1)
     *                  + IMATAB(ISYMA,ISYMB)
C
               NTOTDC = MAX(NMATAB(ISYMDC),1)
               NTOTA  = MAX(NVIR(ISYMA),1)
C
               CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),
     *                    NMATAB(ISYMDC),
     *                    -ONE,WORK(KOFF1),NTOTDC,TMAT(KOFF2),NTOTDC,
     *                    ONE,WORK(KOFF3),NTOTA)
C
            END DO
         END DO
      END DO
C
C
C Distribute O^I(abj) in result vector
C
      CALL CC3_RI_ABJ(XI2,WORK(KRABJ),ISYMI,I,ISYRES)
C
      CALL QEXIT('CC3_CY2V_JK')
      RETURN
      END
C  /* Deck sort_l2_ij_ab */
      SUBROUTINE SORT_L2_IJ_AB(I,ISYMI,J,ISYMJ,XAB,XAIBJ,ISYINT)
C
C-------------------------------
C     Sort from L2 integral stored as X(aI,bJ) to X^IJ(a,b)
C-------------------------------
C
      IMPLICIT NONE
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYMI,ISYMJ,ISYINT 
      INTEGER ISYMB,ISYMBIJ,ISYMA,ISYMAI,ISYMBJ,NBJ,NAB,NAI,NAIBJ 
      INTEGER ISYMIJ
      INTEGER INDEX
C
#if defined (SYS_CRAY)
      REAL XAIBJ(*),XAB(*)
#else
      DOUBLE PRECISION XAIBJ(*),XAB(*) 
#endif
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('L2_IJ_AB')
C
C
      ISYMIJ = MULD2H(ISYMI,ISYMJ)
C
      DO ISYMB = 1,NSYM
C
         ISYMBIJ = MULD2H(ISYMB,ISYMIJ)
         ISYMA   = MULD2H(ISYMBIJ,ISYINT)
         ISYMAI  = MULD2H(ISYMA,ISYMI)
         ISYMBJ  = MULD2H(ISYMB,ISYMJ)
C
         DO B = 1,NVIR(ISYMB)
C
            NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
C
            DO A = 1,NVIR(ISYMA)
C
               NAB = IMATAB(ISYMA,ISYMB)+ NVIR(ISYMA)*(B - 1) + A
               NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A
C
               NAIBJ = IT2AM(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
C
               XAB(NAB) = XAIBJ(NAIBJ)
C
            END DO
         END DO
      END DO
C
      CALL QEXIT('L2_IJ_AB')
C
      RETURN
      END
C  /* Deck sort_l2_j_bai */
      SUBROUTINE SORT_L2_J_BAI(J,ISYMJ,XBAI,XAIBJ,ISYINT)
C
C-------------------------------
C     Sort from L2 integral stored as X(ai,bJ) to X^J(bai)
C-------------------------------
C
      IMPLICIT NONE
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYMJ,ISYINT 
      INTEGER ISYMBAI,ISYMB,ISYMAI,ISYMA,ISYMI,ISYMBJ,NBJ,NBAI,NAI,NAIBJ
      INTEGER ISYMAB
      INTEGER INDEX
C
#if defined (SYS_CRAY)
      REAL XAIBJ(*),XBAI(*)
#else
      DOUBLE PRECISION XAIBJ(*),XBAI(*) 
#endif
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('L2_J_BAI')
C
C
C     Sort from L2 integral stored as X(ai,bJ) to X^J(bai)
      ISYMBAI = MULD2H(ISYINT,ISYMJ)
      DO ISYMB = 1,NSYM
C
         ISYMAI = MULD2H(ISYMB,ISYMBAI)
         DO ISYMA = 1,NSYM
C
            ISYMAB  = MULD2H(ISYMA,ISYMB)
            ISYMI   = MULD2H(ISYMAI,ISYMA)
            ISYMBJ  = MULD2H(ISYMB,ISYMJ)
C
            DO B = 1,NVIR(ISYMB)
C
               NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
C
               DO A = 1,NVIR(ISYMA)
C
                  DO I = 1,NRHF(ISYMI)

                     NBAI = IMAABI(ISYMAB,ISYMI)+ NMATAB(ISYMAB)*(I - 1)
     *                       +IMATAB(ISYMB,ISYMA) + NVIR(ISYMB)*(A-1)+ B
                     NAI  = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A
C
                     NAIBJ = IT2AM(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
C
                     XBAI(NBAI) = XAIBJ(NAIBJ)
C
                  END DO
               END DO
            END DO
         END DO
      END DO

C
      CALL QEXIT('L2_J_BAI')
C
      RETURN
      END
C  /* Deck CC3_RIJ_AB */
      SUBROUTINE CC3_RIJ_AB(OMEGA2,XAB,ISYMI,I,ISYMJ,J,ISYRES)
C
C distribute matrix elements in XAB(A,B) to result vector OMEGA2(AIBJ)
C
      IMPLICIT NONE 
C
      INTEGER ISYMA, ISYMB, ISYRES, INDEX
      INTEGER ISYAB, ISYIJ, ISYMJ, ISYMI, ISYMBJ, ISYMAI
      INTEGER NBJ, NAB, NAI, NAIBJ 
#if defined (SYS_CRAY) 
      REAL OMEGA2(*), XAB(*)
      REAL TWO
#else
      DOUBLE PRECISION OMEGA2(*), XAB(*)
      DOUBLE PRECISION TWO
#endif
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccsdinp.h>
C
      PARAMETER (TWO = 2.0D0)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC3_RIJ_AB')
C
      ISYIJ = MULD2H(ISYMI,ISYMJ)
      ISYAB = MULD2H(ISYIJ,ISYRES)
C
      DO 300 ISYMB = 1,NSYM
         ISYMA  = MULD2H(ISYAB,ISYMB)
         ISYMBJ = MULD2H(ISYMB,ISYMJ)
         ISYMAI = MULD2H(ISYMA,ISYMI)
C
         DO 310 B = 1,NVIR(ISYMB)
            NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
            IF (ISYMAI .EQ. ISYMBJ) THEN
C
               DO 320 A = 1,NVIR(ISYMA)
                  NAB = IMATAB(ISYMA,ISYMB) + NVIR(ISYMA)*(B - 1) + A
                  NAI = IT1AM(ISYMA,ISYMI)  + NVIR(ISYMA)*(I - 1) + A
                  IF (NAI .EQ. NBJ) XAB(NAB) = TWO*XAB(NAB)
                  NAIBJ = IT2AM(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
                  OMEGA2(NAIBJ) = OMEGA2(NAIBJ) + XAB(NAB)
  320          CONTINUE
C
            ELSE IF (ISYMAI .LT. ISYMBJ) THEN
C
               DO 330 A = 1,NVIR(ISYMA)
                  NAB = IMATAB(ISYMA,ISYMB) + NVIR(ISYMA)*(B - 1) + A
                  NAI = IT1AM(ISYMA,ISYMI)  + NVIR(ISYMA)*(I - 1) + A
                  NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                  + NT1AM(ISYMAI)*(NBJ-1) + NAI
                  OMEGA2(NAIBJ) = OMEGA2(NAIBJ) + XAB(NAB)
  330          CONTINUE
C
            ELSE IF (ISYMBJ .LT. ISYMAI) THEN
C
               DO 340 A = 1,NVIR(ISYMA)
                  NAB = IMATAB(ISYMA,ISYMB) + NVIR(ISYMA)*(B - 1) + A
                  NAI = IT1AM(ISYMA,ISYMI)  + NVIR(ISYMA)*(I - 1) + A
                  NAIBJ = IT2AM(ISYMBJ,ISYMAI)
     *                  + NT1AM(ISYMBJ)*(NAI-1) + NBJ
                  OMEGA2(NAIBJ) = OMEGA2(NAIBJ) + XAB(NAB)
  340          CONTINUE
C
            ENDIF
  310    CONTINUE
C
  300 CONTINUE
C
      CALL QEXIT('CC3_RIJ_AB') 
C
      RETURN
      END
C  /* Deck cc3_ri_abj */
      SUBROUTINE CC3_RI_ABJ(OMEGA2,XABJ,ISYMI,I,ISYRES)
C
C distribute matrix elements in XABJ(ABJ) to result vector OMEGA2(AIBJ)
C
      IMPLICIT NONE
C
      INTEGER ISYMI, ISYRES, INDEX
      INTEGER ISYABJ,ISYMA,ISYMBJ,ISYMAI,ISYMAA,NAI,NAAI,ISYMJ,ISYMB 
      INTEGER ISYMAB,NBJ,KOFF1,KOFF2,KOFF5,KOFF6
#if defined (SYS_CRAY)
      REAL OMEGA2(*), XABJ(*)
      REAL TWO
#else
      DOUBLE PRECISION OMEGA2(*), XABJ(*)
      DOUBLE PRECISION TWO
#endif
C
      PARAMETER (TWO = 2.0D0)
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC3_RI_ABJ')
C
      ISYABJ = MULD2H(ISYMI,ISYRES)
C
      DO 100 ISYMA = 1,NSYM
C
         ISYMBJ = MULD2H(ISYMA,ISYABJ)
         ISYMAI = MULD2H(ISYMA,ISYMI)
         ISYMAA = MULD2H(ISYMA,ISYMA)
C
         IF (NT1AM(ISYMAI) .LE. 0 ) GO TO 100
C
         DO 110 A = 1,NVIR(ISYMA)
C
            NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A
C
            IF (ISYMAI .EQ. ISYMBJ) THEN
C
               NAAI = IMAABI(ISYMAA,ISYMI)
     *              + NMATAB(ISYMAA)*(I-1)
     *              + IMATAB(ISYMA,ISYMA)
     *              + NVIR(ISYMA)*(A-1) +A
C
               XABJ(NAAI) = TWO * XABJ(NAAI)
C
               DO 230 ISYMJ = 1,NSYM
                  ISYMB  = MULD2H(ISYMBJ,ISYMJ)
                  ISYMAB = MULD2H(ISYMA,ISYMB)
                  DO 231 J = 1,NRHF(ISYMJ)
                     DO 232 B = 1,NVIR(ISYMB) 
                        NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J-1) + B
C
                        KOFF1 = IT2AM(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
                        KOFF2 = IMAABI(ISYMAB,ISYMJ)
     *                        + NMATAB(ISYMAB)*(J-1)
     *                        + IMATAB(ISYMA,ISYMB)
     *                        + NVIR(ISYMA)*(B-1) +A
C
C
                        OMEGA2(KOFF1) = OMEGA2(KOFF1) + XABJ(KOFF2)
C
  232                CONTINUE
  231             CONTINUE
  230          CONTINUE
C
            ENDIF
C
            IF (ISYMAI .LT. ISYMBJ) THEN
C
               DO 240 ISYMJ = 1,NSYM
                  ISYMB  = MULD2H(ISYMBJ,ISYMJ)
                  ISYMAB = MULD2H(ISYMA,ISYMB)
                  DO 241 J = 1,NRHF(ISYMJ)
                     DO 242 B = 1,NVIR(ISYMB)
                        NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J-1) + B
C
                        KOFF5 = IT2AM(ISYMAI,ISYMBJ)
     *                        + NT1AM(ISYMAI)*(NBJ-1) + NAI
                        KOFF6 = IMAABI(ISYMAB,ISYMJ)
     *                        + NMATAB(ISYMAB)*(J-1)
     *                        + IMATAB(ISYMA,ISYMB)
     *                        + NVIR(ISYMA)*(B-1) +A
C
                        OMEGA2(KOFF5) = OMEGA2(KOFF5) + XABJ(KOFF6)
C
  242                CONTINUE
  241             CONTINUE
  240          CONTINUE
C
            ENDIF
C
            IF (ISYMAI .GT. ISYMBJ) THEN
C
C
               DO 250 ISYMJ = 1,NSYM
                  ISYMB  = MULD2H(ISYMBJ,ISYMJ)
                  ISYMAB = MULD2H(ISYMA,ISYMB)
                  DO 251 J = 1,NRHF(ISYMJ)
                     DO 252 B = 1,NVIR(ISYMB)
                        NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J-1) + B
C
                        KOFF5 = IT2AM(ISYMAI,ISYMBJ)
     *                        + NT1AM(ISYMBJ)*(NAI-1) + NBJ
                        KOFF6 = IMAABI(ISYMAB,ISYMJ)
     *                        + NMATAB(ISYMAB)*(J-1)
     *                        + IMATAB(ISYMA,ISYMB)
     *                        + NVIR(ISYMA)*(B-1) +A

                        OMEGA2(KOFF5) = OMEGA2(KOFF5) + XABJ(KOFF6)
C
  252                CONTINUE
  251             CONTINUE
  250          CONTINUE
C
            ENDIF
C
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CC3_RI_ABJ')
C
      RETURN
      END
C  /* Deck CC3_MABI */
      SUBROUTINE CC3_MABI(RMAT1,RMAT2,ISYMJ,J,ISYRES) 
C
C   order matrix  M^J_(abi) as  M^J_(bai)
C
      IMPLICIT NONE
C
      INTEGER ISYMJ,ISYRES
      INTEGER ISYMABI,ISYMBAI,ISYMI,ISYMAB,ISYMBA,ISYMB,ISYMA,ISYAI 
      INTEGER KOFF1, KOFF2 
#if defined (SYS_CRAY)
      REAL RMAT1(*), RMAT2(*)
#else
      DOUBLE PRECISION RMAT1(*), RMAT2(*)
#endif
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccsdinp.h>

C
      CALL QENTER('CC3_MABI')
C
      ISYMABI   = MULD2H(ISYRES,ISYMJ)
      ISYMBAI  = ISYMABI
      DO ISYMI = 1,NSYM
         ISYMAB = MULD2H(ISYMI,ISYMABI)
         ISYMBA = ISYMAB
         DO ISYMB = 1,NSYM
            ISYMA  = MULD2H(ISYMAB,ISYMB)
            ISYAI  = MULD2H(ISYMA,ISYMI)
            DO I = 1,NRHF(ISYMI)
               DO B = 1,NVIR(ISYMB)
                  DO A = 1,NVIR(ISYMA)
                     KOFF1 = IMAABI(ISYMAB,ISYMI)
     *                     + NMATAB(ISYMAB)*(I-1)
     *                     + IMATAB(ISYMA,ISYMB)
     *                     + NVIR(ISYMA)*(B-1) + A 
                     KOFF2 = IMAABI(ISYMBA,ISYMI)
     *                     + NMATAB(ISYMBA)*(I-1)
     *                     + IMATAB(ISYMB,ISYMA)
     *                     + NVIR(ISYMB)*(A-1) + B 
                     RMAT2(KOFF2) = RMAT1(KOFF1)
                  END DO
               END DO
            END DO 
         END DO 
      END DO 
C
      CALL QEXIT('CC3_MABI')
C
      RETURN
      END
