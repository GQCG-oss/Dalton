C
C...   Copyright (c) 2001 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 1.2
C...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
C...   O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
C...   T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
C...   H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
C...   H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
C...   T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
C...   B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
C...   and O. Vahtras"
C...
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html
C
#include <single.h>
C
C  /* Deck cc3_aden */
      SUBROUTINE CC3_ADEN(LISTL,IDLSTL,LISTR,IDLSTR,
     *                            XLAMDP0,XLAMDH0,FOCK0,
     *                            DIJ,DAB,DIA,
     *                            WORK,LWORK,
     *                            LUDELD,FNDELD,LUCKJD,FNCKJD,LUDKBC,
     *                            FNDKBC,LUTOC,FNTOC,LU3VI,FN3VI,
     *                            LUDKBC3,FNDKBC3,LU3FOPX,FN3FOPX,
     *                            LU3FOP2X,FN3FOP2X)
C
      IMPLICIT NONE
#include <priunit.h>
#include <dummy.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      CHARACTER LISTL*3, LISTR*3
      CHARACTER*(*) FNTOC, FN3VI, FNDKBC3, FN3FOPX, FN3FOP2X
      CHARACTER*(*) FNDKBC,FNDELD,FNCKJD
      CHARACTER*5 FN3FOP
      CHARACTER*10 MODEL
C
      PARAMETER (FN3FOP  = 'PTFOP')
C
      INTEGER IDLSTL,IDLSTR,LWORK
      INTEGER LUTOC, LU3VI, LUDKBC3, LU3FOPX, LU3FOP2X
      INTEGER LUDKBC,LUDELD,LUCKJD
      INTEGER LU3FOP
      INTEGER ISYM0,KT1AMP,KLAMP0,KLAMH0,KEND1,LWRK1,IOPT
C
#if defined (SYS_CRAY)
      REAL XLAMDP0(*),XLAMDH0(*),FOCK0(*)
      REAL DAB(*),DIJ(*),DIA(*)
      REAL WORK(LWORK)
#else
      DOUBLE PRECISION XLAMDP0(*),XLAMDH0(*),FOCK0(*)
      DOUBLE PRECISION DAB(*),DIJ(*),DIA(*)
      DOUBLE PRECISION WORK(LWORK)
#endif
C
      CALL QENTER('CC3_ADEN')

      ISYM0 = 1
C
      KT1AMP = 1
      KLAMP0 = KT1AMP + NT1AM(ISYM0)
      KLAMH0 = KLAMP0 + NLAMDT
      KEND1 = KLAMH0 + NLAMDT
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CC3_ADEN (1)')
      ENDIF
C
*---------------------------------------------------------------------*
*     initialize 0.th-order Lambda:
*---------------------------------------------------------------------*
      IOPT = 1
      CALL CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KT1AMP),DUMMY)

      CALL LAMMAT(WORK(KLAMP0),WORK(KLAMH0),WORK(KT1AMP),
     &            WORK(KEND1),LWRK1)
C

C
C     Open the file
C
      LU3FOP  = 0
      CALL WOPEN2(LU3FOP,FN3FOP,64,0)
C
careful
      call dzero(DIJ,NRHFT*NRHFT)
      call dzero(DAB,NVIRT*NVIRT)
      call dzero(DIA,NVIRT*NRHFT)
careful
      CALL CC3_ADENOCC(LISTL,IDLSTL,LISTR,IDLSTR,
     *                            WORK(KLAMP0),WORK(KLAMH0),FOCK0,
     *                            DIJ,DAB,DIA,
     *                            WORK(KEND1),LWRK1,
     *                            LUDELD,FNDELD,LUCKJD,FNCKJD,LUDKBC,
     *                            FNDKBC,LUTOC,FNTOC,LU3VI,FN3VI,
     *                            LUDKBC3,FNDKBC3,LU3FOP,FN3FOP,
     *                            LU3FOPX,FN3FOPX,
     *                            LU3FOP2X,FN3FOP2X)
C
      write(lupri,*)'DIJ after CC3_ADENOCC '
      call output(DIJ,1,NRHFT,1,NRHFT,NRHFT,NRHFT,1,lupri)
      write(lupri,*)'DAB after CC3_ADENOCC '
      call output(DAB,1,NVIRT,1,NVIRT,NVIRT,NVIRT,1,lupri)
      write(lupri,*)'DIA after CC3_ADENOCC '
      call output(DIA,1,NVIRT,1,NRHFT,NVIRT,NRHFT,1,lupri)
c
      CALL WCLOSE2(LU3FOP,FN3FOP,'KEEP')
careful
       stop'after adenocc'
careful

C
C----------
C     End.
C----------
C
      CALL QEXIT('CC3_ADEN')
C
      RETURN
      END
C  /* Deck cc3_adenocc */
      SUBROUTINE CC3_ADENOCC(LISTL,IDLSTL,LISTR,IDLSTR,
     *                            XLAMDP0,XLAMDH0,FOCK0,
     *                            DIJ,DAB,DIA,
     *                            WORK,LWORK,
     *                            LUDELD,FNDELD,LUCKJD,FNCKJD,LUDKBC,
     *                            FNDKBC,LUTOC,FNTOC,LU3VI,FN3VI,
     *                            LUDKBC3,FNDKBC3,LU3FOP,FN3FOP,
     *                            LU3FOPX,FN3FOPX,
     *                            LU3FOP2X,FN3FOP2X)
*---------------------------------------------------------------------*
*
*    Purpose: 
*            
*    (1) Construct T30^LK
*
*    (2) Construct T3bar0^LK
*
*    (3) Construct W3bar0^LK
*             
*
*    Written by Poul Jorgensen and Filip Pawlowski, Fall 2002, Aarhus
*            
*=====================================================================*
C
      IMPLICIT NONE
#include<ccl1rsp.h>
#include<ccr1rsp.h>
#include<ccorb.h>
#include<ccsdsym.h>
#include<dummy.h>
#include<priunit.h>
#include<iratdef.h>
#include<ccinftap.h>
C
      INTEGER ISYM0
      PARAMETER(ISYM0 = 1)
C
      INTEGER IDLSTL,IDLSTR,IDLSTL0,LWORK
      INTEGER LUTOC, LU3VI, LUDKBC3, LU3FOPX, LU3FOP2X, LU3FOP
      INTEGER LUDKBC,LUDELD,LUCKJD
      INTEGER ISYML0,ISYML1,ISYMR1,ISINT1,ISINT2,ISINTR1,ISYFCKR1
      INTEGER ISYMK,ISYML,ISYMT3,ISYMKL,ISYT30KL
      INTEGER IOPT,LENGTH
      INTEGER KFOCKD,KFCKBA,KT2TP,KL1AM,KL2TP,KEND0,LWRK0
      INTEGER KL1,KL2,KFOCKL1,KR1,KR2,KFOCKR1,KEND1,LWRK1
      INTEGER KXIAJB,KT3BOG1,KT3BOL1,KT3BOG2,KT3BOL2,KT3OG1,KT3OG2
      INTEGER KLAMPR1,KLAMHR1,KT30KL
      INTEGER KFCKR1CK,KW3BXOG1,KW3BXOL1,KW3BXOGX1,KW3BXOLX1,KT3VIJG1
      INTEGER ISYMT3B,ISYT3B0KL,ISYW3BXKL
      INTEGER KXGADCK,KXLADCK
      INTEGER KT3B0KL,KW3BXKL,ISYMW3BX
      INTEGER KT3BOG2X,KT3BOL2X,KXGADCKX,KXLADCKX
      INTEGER ISYMTETAX,ISTETAXKL
      INTEGER KTETAXKL
C
      CHARACTER LISTL*3, LISTR*3, LISTL0*3
      CHARACTER*(*) FNTOC, FN3VI, FNDKBC3, FN3FOPX, FN3FOP2X, FN3FOP
      CHARACTER*(*) FNDKBC,FNDELD,FNCKJD
      CHARACTER LABELL1*8,LABELR1*8
C
      LOGICAL   LOCDBG,LORXL1
      PARAMETER (LOCDBG = .FALSE.)
C
C
      integer kx3am
C
C
#if defined (SYS_CRAY)
      REAL XLAMDP0(*),XLAMDH0(*),FOCK0(*) 
      REAL DAB(*),DIJ(*),DIA(*)
      REAL WORK(LWORK)
      REAL FREQL1,FREQR1
      REAL DDOT
#else
      DOUBLE PRECISION XLAMDP0(*),XLAMDH0(*),FOCK0(*) 
      DOUBLE PRECISION DAB(*),DIJ(*),DIA(*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION FREQL1,FREQR1
      DOUBLE PRECISION DDOT
#endif
C
      CALL QENTER('CC3_ADENOCC')
C
C------------------------------------------------------------
C     some initializations:
C------------------------------------------------------------
C
      LISTL0 = 'L0 '
      IDLSTL0 = 0
      ISYML0 = ISYM0
C
      ISYMT3 = ISYM0
      ISYMT3B = ISYM0

      IF (LISTL(1:3).EQ.'L1 ') THEN
         ! get symmetry, frequency and integral label for left list 
         ! from common blocks defined in ccl1rsp.h
         ISYML1  = ISYLRZ(IDLSTL)
         FREQL1  = FRQLRZ(IDLSTL)
         LABELL1 = LRZLBL(IDLSTL)
         LORXL1  = LORXLRZ(IDLSTL)

         IF (LORXL1) CALL QUIT('NO ORBITAL RELAX. IN CC3_ADENOCC')

      ELSE 
         CALL QUIT('Unkown left list in CC3_ADENOCC')
      END IF

      IF (LISTR(1:3).EQ.'R1 ') THEN
         ! get symmetry, frequency and integral label for right list 
         ! from common blocks defined in ccr1rsp.h
        ISYMR1  = ISYLRT(IDLSTR)
        FREQR1  = FRQLRT(IDLSTR)
        LABELR1 = LRTLBL(IDLSTR)
      ELSE
         CALL QUIT('Unkown right list in CC3_ADENOCC')
      END IF
C
C---------------------------------------------------------------------
C     initial allocations, orbital energy, fock matrix and T2 and L2 :
C---------------------------------------------------------------------
C
      KFOCKD  = 1
      KFCKBA  = KFOCKD  + NORBTS
      KT2TP   = KFCKBA  + NT1AMX 
      KL1AM   = KT2TP   + NT2SQ(ISYM0)
      KL2TP   = KL1AM   + NT1AM(ISYML0)
      KEND0   = KL2TP   + NT2SQ(ISYML0)
      LWRK0   = LWORK   - KEND0
C
      KL1     = KEND0
      KL2     = KL1     + NT1AM(ISYML1)
      KFOCKL1 = KL2     + NT2SQ(ISYML1)
      KR1     = KFOCKL1 + N2BST(ISYML1)
      KR2     = KR1     + NT1AM(ISYMR1)
      KFOCKR1 = KR2     + NT2SQ(ISYMR1)
      KEND1   = KFOCKR1 + N2BST(ISYMR1)
      LWRK1   = LWORK   - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CC3_ADENOCC (1)')
      ENDIF
C
C-------------------------------------
C     Read T2 amplitudes 
C-------------------------------------
C
      IOPT = 2
      CALL GET_T1_T2(IOPT,.FALSE.,DUMMY,WORK(KT2TP),'R0',0,ISYM0,
     *                WORK(KEND1),LWRK1)
C
      IF (LOCDBG) WRITE(LUPRI,*) 'Norm of T2TP ',
     *    DDOT(NT2SQ(ISYM0),WORK(KT2TP),1,WORK(KT2TP),1)
C
C-------------------------------------
C     Read L1 and L2 amplitudes 
C-------------------------------------
C
      IOPT = 3
      CALL GET_T1_T2(IOPT,.FALSE.,WORK(KL1AM),WORK(KL2TP),LISTL0,
     *               IDLSTL0,ISYML0,WORK(KEND1),LWRK1)
C
      IF (LOCDBG) WRITE(LUPRI,*) 'Norm of L2TP ',
     *    DDOT(NT2SQ(ISYML0),WORK(KL2TP),1,WORK(KL2TP),1)

C
C---------------------------------------------------------------
C     Read canonical orbital energies and delete frozen orbitals 
C     in Fock diagonal, if required
C---------------------------------------------------------------
C
      CALL GET_ORBEN(WORK(KFOCKD),WORK(KEND1),LWRK1)
C
C--------------------------------------------
C     Sort the Fock matrix to get F(ck) block
C--------------------------------------------
C
      CALL SORT_FOCKCK(WORK(KFCKBA),FOCK0,ISYM0)
C
C---------------------------------------------------------------------
C     Read information for L1 list
C---------------------------------------------------------------------
C
      IF (LISTL(1:3).EQ.'L1 ') THEN
C
C---------------------------------------------------------------------
C     Read the matrix the property integrals and trasform it to lambda 
C     basis (unsorted - need in WBX_JK_ETA)
C---------------------------------------------------------------------
C
         CALL GET_FOCKX(WORK(KFOCKL1),LABELL1,IDLSTL,ISYML1,XLAMDP0,
     *                  XLAMDH0,WORK(KEND1),LWRK1)
C
C-------------------------------------
C     Read L1 and L2 multipliers 
C-------------------------------------
C
         IOPT  = 3
         CALL GET_T1_T2(IOPT,.FALSE.,WORK(KL1),WORK(KL2),LISTL,
     *                  IDLSTL,ISYML1,WORK(KEND1),LWRK1)
      END IF
C
C---------------------------------------------------------------------
C     Read information for R1 list
C---------------------------------------------------------------------
C
      IF (LISTR(1:3).EQ.'R1 ') THEN

C---------------------------------------------------------------------
C     Read the matrix the property integrals and trasform it to lambda 
C     basis
C---------------------------------------------------------------------
C        
         CALL GET_FOCKX(WORK(KFOCKR1),LABELR1,IDLSTR,ISYMR1,XLAMDP0,
     *                  XLAMDH0,WORK(KEND1),LWRK1)
C
C-------------------------------------
C     Read R1 and R2 amplitudes 
C-------------------------------------
C     
         IOPT  = 3
         CALL GET_T1_T2(IOPT,.FALSE.,WORK(KR1),WORK(KR2),LISTR,
     *                  IDLSTR,ISYMR1,WORK(KEND1),LWRK1)
      END IF
C
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT If we want to sum the T3 amplitudes
COMMENT COMMENT COMMENT
C
      if (.false.) then
         kx3am  = kend1
         kend1 = kx3am + nrhft*nrhft*nrhft*nvirt*nvirt*nvirt
         call dzero(work(kx3am),nrhft*nrhft*nrhft*nvirt*nvirt*nvirt)
         lwrk0 = lwork - kend1
         if (lwrk0 .lt. 0) then
            write(lupri,*) 'Memory available : ',lwork
            write(lupri,*) 'Memory needed    : ',kend1
            call quit('Insufficient space (kx3am) in CC3_ADENOCC (2)')
         END IF
      endif
COMMENT COMMENT COMMENT
C      write(lupri,*) 'WMAT after dzero'
C      call print_pt3(work(kx3am),ISYMBY,4)
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
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
      ISINTR1   = MULD2H(ISINT1,ISYMR1)
      ISYFCKR1  = MULD2H(ISYMOP,ISYMR1)

      KXIAJB    = KEND1
      KEND1     = KXIAJB    + NT2AM(ISYM0)

      KT3BOG1   = KEND1
      KT3BOL1   = KT3BOG1   + NTRAOC(ISYM0)
      KT3BOG2   = KT3BOL1   + NTRAOC(ISYM0)
      KT3BOL2   = KT3BOG2   + NTRAOC(ISYM0)
      KT3OG1    = KT3BOL2   + NTRAOC(ISYM0)
      KT3OG2    = KT3OG1    + NTRAOC(ISINT2)
      KLAMPR1   = KT3OG2    + NTRAOC(ISINT2)
      KLAMHR1   = KLAMPR1   + NLAMDT
      KEND1     = KLAMHR1   + NLAMDT
C
      KFCKR1CK  = KEND1
      KW3BXOG1  = KFCKR1CK  + NT1AM(ISYFCKR1)
      KW3BXOL1  = KW3BXOG1  + NTRAOC(ISYM0)
      KW3BXOGX1 = KW3BXOL1  + NTRAOC(ISYM0)
      KW3BXOLX1 = KW3BXOGX1 + NTRAOC(ISINTR1)
      KT3VIJG1  = KW3BXOLX1 + NTRAOC(ISINTR1)
      KEND1     = KT3VIJG1  + NMAABCI(ISYM0)
      LWRK1     = LWORK     - KEND1
C
      KT3BOG2X   = KEND1
      KT3BOL2X   = KT3BOG2X + NTRAOC(ISINTR1)
      KEND1      = KT3BOL2X + NTRAOC(ISINTR1)
C
      KXGADCK   = KEND1 
      KXLADCK   = KXGADCK + NMAABCI(ISYM0)
      KEND1     = KXLADCK + NMAABCI(ISYM0)
      LWRK1     = LWORK     - KEND1
C
      KXGADCKX   = KEND1 
      KXLADCKX   = KXGADCKX + NMAABCI(ISINTR1)
      KEND1     = KXLADCKX + NMAABCI(ISINTR1)
      LWRK1     = LWORK     - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in CC3_ADENOCC (3)')
      END IF
C
C------------------------
C     Construct L(ia,jb).
C------------------------
C
      LENGTH = IRAT*NT2AM(ISYM0)

      REWIND(LUIAJB)
      CALL READI(LUIAJB,LENGTH,WORK(KXIAJB))

      CALL CCSD_TCMEPK(WORK(KXIAJB),1.0D0,ISYM0,1)

C
C--------------------------------------------------------------
C     Prepare to construct the integrals (occupied and virtual)
C--------------------------------------------------------------
C
C
C--------------------------
C     Get Lambda_Y matrices
C--------------------------
C
         CALL GET_LAMBDAX(WORK(KLAMPR1),WORK(KLAMHR1),LISTR,IDLSTR,
     *                    ISYMR1,XLAMDP0,XLAMDH0,WORK(KEND1),LWRK1)
C
C------------------------------------------------------------------
C        Calculate the F^Y matrix (kc elements evaluated and stored 
C        as ck)
C------------------------------------------------------------------
C
         CALL CC3LR_MFOCK(WORK(KFCKR1CK),WORK(KR1),WORK(KXIAJB),
     *                    ISYFCKR1)
C
C-----------------------------------------------------------------
C     Construct occupied integrals which are required to calculate    
C     t3bar_0 multipliers                                             
C-----------------------------------------------------------------
C
      CALL INTOCC_T3BAR0(LUTOC,FNTOC,XLAMDH0,ISYM0,WORK(KT3BOG1),
     *                   WORK(KT3BOL1),WORK(KT3BOG2),WORK(KT3BOL2),
     *                   WORK(KEND1),LWRK1)
C
C-----------------------------------------------------------------
C     Construct occupied integrals which are required to calculate    
C     t3bar_X multipliers                                             
C-----------------------------------------------------------------
C
      CALL INTOCC_T3BAR0(LUTOC,FNTOC,WORK(KLAMHR1),ISINTR1,
     *                   WORK(KT3BOG1),
     *                   WORK(KT3BOL1),WORK(KT3BOG2X),WORK(KT3BOL2X),
     *                   WORK(KEND1),LWRK1)
C
C-----------------------------------------------------------------
C     Construct occupied integrals which are required to calculate    
C     t3_0 amplitudes
C-----------------------------------------------------------------
C
      CALL INTOCC_T30(LUCKJD,FNCKJD,XLAMDP0,ISINT2,WORK(KT3OG1),
     *                WORK(KT3OG2),WORK(KEND1),LWRK1)
C
C-----------------------------------------------------------------
C     Construct occupied integrals which are required to calculate    
C     t3bar_Y multipliers                                             
C-----------------------------------------------------------------
C
      CALL INTOCC_T3BARX(LUTOC,FNTOC,ISYMOP,XLAMDH0,ISYM0,ISINT1,
     *                   WORK(KLAMHR1),ISYMR1,ISINTR1,WORK(KW3BXOG1),
     *                   WORK(KW3BXOL1),WORK(KW3BXOGX1),WORK(KW3BXOLX1),
     *                   WORK(KEND1),LWRK1)
C
C----------------------------------------------
C     Get virtual integrals for t30 amplitudes
C     KT3VIJG1 : (ck|da) sorted as I(ad|ck)
C----------------------------------------------
C
      CALL INTVIR_T30_IJ(WORK(KT3VIJG1),ISYM0,XLAMDH0,LUDELD,FNDELD,
     *                   WORK(KEND1),LWRK1)
C
C----------------------------------------------
C     Get virtual integrals for t3b0 multipliers
C     KXGADCK g(kcad) = (kc ! ad) sorted as I(adck)
C     KXLADCK L(kcad) sorted as I(adck)
C----------------------------------------------
C
      CALL INTVIR_T3B0_JK(WORK(KXGADCK),WORK(KXLADCK),ISYM0,XLAMDP0,
     *                         LU3VI,FN3VI,LU3FOP,FN3FOP,
     *                         WORK(KEND1),LWRK1)
C
C----------------------------------------------
C     Get virtual integrals for t3b0 multipliers
C----------------------------------------------
C
      CALL INTVIR_T3B0_JK(WORK(KXGADCKX),WORK(KXLADCKX),ISYMR1,
     *                    WORK(KLAMPR1),
     *                    LU3VI,FN3VI,LU3FOP,FN3FOP,
     *                    WORK(KEND1),LWRK1)
C
C----------------------------
C     Loop over K
C----------------------------
C
      ISYMW3BX = MULD2H(ISYM0,ISYML1)
      ISYMTETAX = MULD2H(ISYM0,ISYMR1)
      DO ISYMK = 1,NSYM

         DO K = 1,NRHF(ISYMK)
C
            DO ISYML = 1,NSYM
C
               ISYMKL = MULD2H(ISYMK,ISYML)
               ISYT30KL = MULD2H(ISYMKL,ISYMT3)
               ISYT3B0KL = MULD2H(ISYMKL,ISYMT3B)
               ISYW3BXKL  = MULD2H(ISYMKL,ISYMW3BX)
               ISTETAXKL  = MULD2H(ISYMKL,ISYMTETAX)
C
               KT30KL = KEND1
               KT3B0KL  = KT30KL + NMAABCI(ISYT30KL)
               KW3BXKL  = KT3B0KL + MAX( NMAABCI(ISYT3B0KL),
     *                                   NMAABCI(ISTETAXKL))
               KTETAXKL    = KW3BXKL + NMAABCI(ISYW3BXKL)
               KEND1   = KTETAXKL + NMAABCI(ISTETAXKL)
               LWRK1  = LWORK  - KEND1

               IF (LWRK1 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Memory available : ',LWORK
                  WRITE(LUPRI,*) 'Memory needed    : ',KEND1
                  CALL QUIT('Insufficient space in CC3_ADENOCC (4)')
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
C                Divide by orbital energy difference and remove 
C                forbidden elements
C
                 CALL T3JK_DIA(WORK(KT30KL),ISYT30KL,ISYML,L,ISYMK,K,
     *                         WORK(KFOCKD))
                 CALL T3_FORBIDDEN_JK(WORK(KT30KL),ISYMT3,ISYML,L,
     *                                ISYMK,K)
C
c                call sum_pt3_jk(work(kt30kl),isyml,l,isymk,k,isyt30kl,
c    *                           work(kx3am),1)
C
C-------------------------------------------
C                 Get T3BAR0^KL multipliers
C-------------------------------------------
C
                  CALL DZERO(WORK(KT3B0KL),NMAABCI(ISYT3B0KL))
C
                  CALL GET_T3B0_JK_O(WORK(KT3B0KL),ISYT3B0KL,
     *                           WORK(KL2TP),ISYML0,
     *                           WORK(KT3BOL2),WORK(KT3BOG2),ISYM0,
     *                           ISYML,L,ISYMK,K,
     *                           WORK(KEND1),LWRK1)
C
                  CALL GET_T3B0_JK_V(WORK(KT3B0KL),ISYT3B0KL,
     *                               WORK(KL2TP),ISYML0,
     *                               WORK(KXGADCK),WORK(KXLADCK),
     *                               ISYM0,ISYML,L,ISYMK,K,
     *                               WORK(KEND1),LWRK1)
C
                  CALL GET_T3B0_JK_L1F(WORK(KT3B0KL),ISYT3B0KL,
     *                            WORK(KL1AM),ISYML0,
     *                            WORK(KXIAJB),ISYM0,
     *                            WORK(KL2TP),ISYML0,
     *                            WORK(KFCKBA),ISYM0,
     *                            ISYML,L,ISYMK,K)
C
C                Divide by orbital energy difference and remove 
C                forbidden elements
C
                 CALL T3JK_DIA(WORK(KT3B0KL),ISYT3B0KL,ISYML,L,ISYMK,K,
     *                         WORK(KFOCKD))
                 CALL T3_FORBIDDEN_JK(WORK(KT3B0KL),ISYMT3B,ISYML,L,
     *                                ISYMK,K)
C
c               call sum_pt3_jk(work(kt3b0kl),isyml,l,isymk,k,isyt3b0kl,
c    *                          work(kx3am),7)
C
C-------------------------------------------
C                 Get W3BARX^KL multipliers
C-------------------------------------------
C
                  CALL DZERO(WORK(KW3BXKL),NMAABCI(ISYW3BXKL))
C           
C                 <L2|[Y,tau3]|HF> + <L3|[Y^,tau3]|HF>
C
                  CALL WBX_JK_ETA(WORK(KT3B0KL),ISYT3B0KL,WORK(KFOCKL1),
     *                            ISYML1,WORK(KW3BXKL),ISYW3BXKL,
     *                            WORK(KL2TP),ISYML0,ISYML,L,ISYMK,K,
     *                            WORK(KEND1),LWRK1)
C
C                 <L2Y|[H^,tau3]|HF>
C
                  CALL WBX_JK_FMAT(WORK(KW3BXKL),ISYW3BXKL,
     *                            WORK(KL2),ISYML1,
     *                            WORK(KFCKBA),ISYM0,
     *                            WORK(KT3BOL2),WORK(KT3BOG2),
     *                            WORK(KXGADCK),WORK(KXLADCK),ISYM0,
     *                            ISYML,L,ISYMK,K,
     *                            WORK(KEND1),LWRK1)
C
C                 <L2|[H^Y,tau3]|HF>
C
                  CALL WBX_JK_FMAT(WORK(KW3BXKL),ISYW3BXKL,
     *                            WORK(KL2TP),ISYML0,
     *                            WORK(KFCKR1CK),ISYFCKR1,
     *                            WORK(KT3BOL2X),WORK(KT3BOG2X),
     *                            WORK(KXGADCKX),WORK(KXLADCKX),ISINTR1,
     *                            ISYML,L,ISYMK,K,
     *                            WORK(KEND1),LWRK1)
C
C                 <L1Y|[H^,tau3]|HF>
C
                  CALL WBX_JK_L1(WORK(KW3BXKL),ISYW3BXKL,
     *                           WORK(KL1),ISYML1, 
     *                           WORK(KXIAJB),ISYM0, 
     *                           ISYML,L,ISYMK,K) 
C
C                Divide by orbital energy difference and remove 
C                forbidden elements
C
                 CALL W3JK_DIA(WORK(KW3BXKL),ISYW3BXKL,ISYML,L,ISYMK,K,
     *                         WORK(KFOCKD),FREQL1)
                 CALL T3_FORBIDDEN_JK(WORK(KW3BXKL),ISYMW3BX,ISYML,L,
     *                                ISYMK,K)
C
c               call sum_pt3_jk(work(kw3bxkl),isyml,l,isymk,k,isyw3bxkl,
c    *                          work(kx3am),4)
C
                CALL DZERO(WORK(KTETAXKL),NMAABCI(ISTETAXKL))
c
                CALL TETAX_JK_BC(WORK(KT30KL),ISYT30KL,WORK(KFOCKR1),
     *                           ISYMR1,WORK(KTETAXKL),ISTETAXKL,
     *                           WORK(KEND1),LWRK1)
C
                CALL DZERO(WORK(KT3B0KL),NMAABCI(ISTETAXKL))
C
                CALL DCOPY(NMAABCI(ISTETAXKL),WORK(KTETAXKL),1,
     *                     WORK(KT3B0KL),1)
C
                CALL W3JK_DIA(WORK(KTETAXKL),ISTETAXKL,ISYML,L,ISYMK,K,
     *                         WORK(KFOCKD),FREQR1)
                CALL T3_FORBIDDEN_JK(WORK(KTETAXKL),ISTETAXKL,ISYML,L,
     *                                ISYMK,K)

c
c               call sum_pt3_jk(work(KTETAXKL),isyml,l,isymk,k,
c    *                          ISTETAXKL,
c    *                          work(kx3am),1)
C
C
c
                CALL ADEN_DIJ_JK(DIJ,WORK(KTETAXKL),ISTETAXKL,
     *                           WORK(KW3BXKL),ISYW3BXKL)
C
                CALL TETAX_JK_A(WORK(KT30KL),ISYT30KL,WORK(KFOCKR1),
     *                           ISYMR1,WORK(KT3B0KL),ISTETAXKL,
     *                           WORK(KEND1),LWRK1)
C
                CALL W3JK_DIA(WORK(KT3B0KL),ISTETAXKL,ISYML,L,ISYMK,K,
     *                         WORK(KFOCKD),FREQR1)
                CALL T3_FORBIDDEN_JK(WORK(KT3B0KL),ISTETAXKL,ISYML,L,
     *                                ISYMK,K)

c
c               call sum_pt3_jk(work(KTETAXKL),isyml,l,isymk,k,
c    *                          ISTETAXKL,
c    *                          work(kx3am),4)
c
                 CALL ADEN_DAB_LM(DAB,WORK(KT3B0KL),ISTETAXKL,
     *                            WORK(KW3BXKL),ISYW3BXKL,
     *                            WORK(KEND1),LWRK1)
C
                 CALL ADEN_DAI_LM(DIA,
     *                       WORK(KL2),ISYML1,
     *                       WORK(KT3B0KL),ISTETAXKL,
     *                       ISYML,L,ISYMK,K,
     *                       WORK(KEND1),LWRK1)
C
               ENDDO   ! L
            ENDDO      ! ISYML
         ENDDO       ! K
      ENDDO          ! ISYMK 
C
c      write(lupri,*) 'T30KL in CC3_ADENOCC'
c      call print_pt3(work(kx3am),isym0,4)
C
C
C-------------
C     End
C-------------
C

      CALL QEXIT('CC3_ADENOCC')
C
      RETURN
      END
C  /* Deckntvir_t30_ij */
      SUBROUTINE INTVIR_T30_IJ(XINTADCK,ISYINT,XLAMH,
     *                         LUDELD,FNDELD,
     *                         WORK,LWORK) 
**********************************************************
*
*     Construvt the integrals used for t30^IJ calculation
*
*     Read virtual integrals I(ck | delta^h D^p) and transform and sort.
*
*     OUTPUT (XINTADCK) :  (ck|Da) sorted as I(aD|ck) 
*
*     P. Jorgensen, F. Pawlowski, 31-01-2003, Aarhus.
**********************************************************
C
      IMPLICIT NONE
#include <ccsdsym.h>
#include <ccorb.h>
#include <priunit.h>
C
      INTEGER ISYINT, LUDELD, LWORK
      INTEGER ISYMD, ISYCKA, KINTVI, KTRVI1, KEND1, LWRK1, IOFF
C
      CHARACTER*(*) FNDELD
C
#if defined (SYS_CRAY)
      REAL XINTADCK(*), XLAMH(*), WORK(LWORK)
#else
      DOUBLE PRECISION XINTADCK(*), XLAMH(*), WORK(LWORK)
#endif
C
      CALL QENTER('INTVIR_T30_IJ')
C

C
C--------------------------
C     Read MO coefficients.
C--------------------------
C

      DO ISYMD = 1, NSYM
         ISYCKA = MULD2H(ISYINT,ISYMD)
C
         KINTVI = 1
         KTRVI1 = KINTVI + NCKA(ISYCKA)
         KEND1  = KTRVI1 + NCKATR(ISYCKA)
         LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in INTVIR_T30_IJ ')
      ENDIF
C
         DO D = 1, NVIR(ISYMD)
C
C-------------------------------------------------------------------------
C           Read virtual integrals I(ck | delta^h D^p) 
C           and transform to I(ck | a^h D^p)
C-------------------------------------------------------------------------
C
            IOFF = ICKAD(ISYCKA,ISYMD) + NCKA(ISYCKA)*(D - 1) + 1
            IF (NCKA(ISYCKA) .GT. 0) THEN
               CALL GETWA2(LUDELD,FNDELD,WORK(KINTVI),IOFF,
     &                     NCKA(ISYCKA))
            ENDIF
C
            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI1),XLAMH,
     *                       ISYMD,D,ISYMOP,WORK(KEND1),LWRK1)
C
C---------------------------------------------------------
C           Integrals sitting as I(ck | a^h D^p) = (ck|Da)
C           Sort it now as I(aD|ck)
C---------------------------------------------------------
C
C
            CALL SORT_INTVIR_T30(XINTADCK,WORK(KTRVI1),
     *                           D,ISYMD,ISYCKA,WORK(KEND1),LWRK1)
         END DO !  D
      END DO !  ISYMD
C
      CALL QEXIT('INTVIR_T30_IJ')
C
      RETURN
      END
C  /* Deck sort_intvir_t30 */
      SUBROUTINE SORT_INTVIR_T30(XINTADCK,XINTCKAD,
     *                           D,ISYMD,ISYCKA,WORK,LWORK)
**************************************************************************
* Sort I(ck | a^h D^p) = (ck|Da) integrals (INTCKAD) as I(aD|ck) (INTADCK).
**************************************************************************
C
      IMPLICIT NONE
#include <ccsdsym.h>
#include <ccorb.h>
#include <priunit.h>
C
      INTEGER ISYMD,ISYCKA,LWORK
      INTEGER ISYMA,ISYMCK,ISYAD,ISYMC,ISYMK,ISYADC
      INTEGER KOFF1,KOFF2
C
#if defined (SYS_CRAY)
      REAL XINTADCK(*),XINTCKAD(*),WORK(LWORK)
#else
      DOUBLE PRECISION XINTADCK(*),XINTCKAD(*),WORK(LWORK)
#endif
C
      CALL QENTER('SORT_INTVIR_T30')


* Sort I^D(cka) = (ck|Da) integrals (INTCKAD) as I(aD|ck) (INTADCK).
      DO ISYMA = 1,NSYM
         ISYMCK = MULD2H(ISYCKA,ISYMA) 
         ISYAD  =  MULD2H(ISYMA,ISYMD)
         DO ISYMC = 1,NSYM
            ISYMK = MULD2H(ISYMCK,ISYMC)
            ISYADC = MULD2H(ISYAD,ISYMC)
            DO A = 1,NVIR(ISYMA)
               DO K = 1,NRHF(ISYMK)
                  DO C = 1,NVIR(ISYMC)
                     KOFF1 = ICKATR(ISYMCK,ISYMA)
     *                     + NT1AM(ISYMCK)*(A-1)
     *                     + IT1AM(ISYMC,ISYMK)
     *                     + NVIR(ISYMC)*(K-1)
     *                     + C
                     KOFF2 = IMAABCI(ISYADC,ISYMK)
     *                     + NMAABC(ISYADC)*(K-1)
     *                     + IMAABC(ISYAD,ISYMC)
     *                     + NMATAB(ISYAD)*(C-1)
     *                     + IMATAB(ISYMA,ISYMD)
     *                     + NVIR(ISYMA)*(D-1)
     *                     + A
                     XINTADCK(KOFF2) = XINTCKAD(KOFF1)
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      CALL QEXIT('SORT_INTVIR_T30')
      RETURN
      END
C  /* Deck get_t30_ij_o */
      SUBROUTINE GET_T30_IJ_O(T30JK,ISYT30JK,T2TP,
     *                        ISYMT2,
     *                        T3OG2,ISYINT,ISYMJ,J,ISYMK,K,
     *                        WORK,LWORK)
***********************************************************
*
*     T3OG2 : (ai | kj) sorted as I(a_1^p,j_2^h,k_2^p,i_1^h)
*                                 I(a,j,k,i)
*     T30JK sitting as (bcai)
***********************************************************
C
C     T30^(abc)_(iJK) = 
C     P(ai,bj,ck) (sum_d t^(ad)_(ij) (ck|bd) ) + 
C    - P(ai,bj,ck) (sum_l t^(ab)_(il) (ck|lj) )
C
C    In this routine we calculate the second (i.e. occupied) contribution:
C
C    T^JK(bcai) = - P(ai,bj,ck) (sum_l t^(ab)_(il) (ck|lj) )
C (1)
C               =    - sum_l t^(ab)_(il) (ck|lj)
C (2)                
C                    - sum_l t^(ba)_(jl) (ck|li)
C (3)
C                    - sum_l t^(bc)_(jl) (ai|lk)
C (4)
C                    - sum_l t^(ac)_(il) (bj|lk)
C (5)
C                    - sum_l t^(ca)_(kl) (bj|li)
C (6)
C                    - sum_l t^(cb)_(kl) (ai|lj)
C
C
C     Filip Pawlowski, Aarhus, Winter 2003
C
      IMPLICIT NONE
#include<ccsdsym.h>
#include<ccorb.h>
#include<priunit.h>
C
      INTEGER ISYT30JK,ISYMT2,ISYINT,ISYMJ,ISYMK,LWORK
      INTEGER ISYKJ,ISYCL,ISYML,ISYMC,ISYBAI
      INTEGER ISYT2BAL,ISYINTLCI,ISYBA,ISYCI
      INTEGER ISYT2BCL,ISYBC,ISYAI
      INTEGER ISYJK,ISYBL,ISYMB,ISYCAI
      INTEGER ISYT2CAL,ISYINTLBI,ISYCA,ISYBI
      INTEGER ISYT2CBL,ISYINTLAI,ISYCB
      INTEGER KT2LBAI,KINTCL,KCBAI,KEND1,LWRK1
      INTEGER KT2BAL,KINTLCI,KBACI
      INTEGER KT2BCL
      INTEGER KT2LCAI,KINTBL
      INTEGER KT2CAL,KINTLBI,KCABI
      INTEGER KT2CBL,KINTLAI
      INTEGER KOFF1,KOFF2,KOFF3
      INTEGER NTOTC,NTOTL,NTOTBA,NTOTBC,NTOTB,NTOTCA,NTOTCB
C
#if defined (SYS_CRAY)
      REAL T30JK(*),T2TP(*),T3OG2(*),WORK(LWORK)
      REAL ONE
#else
      DOUBLE PRECISION T30JK(*),T2TP(*),T3OG2(*),WORK(LWORK)
      DOUBLE PRECISION ONE
#endif
C
      PARAMETER (ONE = 1.0D0)
C
      CALL QENTER('GET_T30_IJ_O')
C
C=================================================
C     Calculate (1)   - sum_l t^(ab)_(il) (ck|lj)
C
C                             T(lbai) I^KJ(cl)
C=================================================
C
C-------------------------------
C     Sort T2TP(blia) as T(lbai)
C-------------------------------
C
      ISYKJ = MULD2H(ISYMK,ISYMJ)
      ISYCL = MULD2H(ISYINT,ISYKJ)
C
      KT2LBAI = 1
      KINTCL  = KT2LBAI + NT2SQ(ISYMT2)
      KCBAI   = KINTCL  + NT1AM(ISYCL)
      KEND1   = KCBAI   + NMAAOBCI(ISYT30JK)
      LWRK1   = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T30_IJ_O (1)')
      END IF
C
      CALL DZERO(WORK(KCBAI),NMAAOBCI(ISYT30JK))
C
      CALL SORT_T2_IABJ(WORK(KT2LBAI),T2TP,ISYMT2)
C
C-----------------------------
C     Sort (ck|lj) = T3OG2(c,j,k,l) as I^KJ(cl)
C-----------------------------
C
      CALL SORT_INT_AJ_IK(WORK(KINTCL),T3OG2,ISYINT,ISYMK,K,ISYMJ,J)
C
C------------------------------------------
C    Multiply I^KJ(cl) T(lbai) = T^JK(cbai)
C------------------------------------------
C
      DO ISYML = 1, NSYM
            ISYMC = MULD2H(ISYCL,ISYML)
            ISYBAI = MULD2H(ISYMT2,ISYML)
C      
            KOFF1 = KINTCL 
     *            + IT1AM(ISYMC,ISYML)
            KOFF2 = KT2LBAI
     *            + IMAJBAI(ISYML,ISYBAI)
            KOFF3 = KCBAI
     *            + IMAAOBCI(ISYMC,ISYBAI)
C
            NTOTC = MAX(NVIR(ISYMC),1)
            NTOTL = MAX(NRHF(ISYML),1)
C
            CALL DGEMM('N','N',NVIR(ISYMC),NMAABI(ISYBAI),NRHF(ISYML),
     *                 ONE,WORK(KOFF1),NTOTC,WORK(KOFF2),NTOTL,
     *                 ONE,WORK(KOFF3),NTOTC)
C
      END DO ! ISYML
C
C     T30JK(bcai) = T30JK(bcai) + T^JK(cbai)
C
C  add_occ(1)
C
      CALL FBACI(T30JK,WORK(KCBAI),ISYT30JK)
C
C=================================================
C     Calculate (2)  - sum_l t^(ba)_(jl) (ck|li) 
C
C                             T^J(bal) I^K(lci)
C=================================================
C
C-------------------------------
C     Sort T2TP(bjla) as T^J(bal)
C-------------------------------
C
      ISYT2BAL = MULD2H(ISYMT2,ISYMJ)
      ISYINTLCI = MULD2H(ISYINT,ISYMK)
C
      KT2BAL = KEND1
      KINTLCI  = KT2BAL + NMAABI(ISYT2BAL)
      KBACI   = KINTLCI + NCKI(ISYINTLCI)
      KEND1 = KBACI + NMAAB_CI(ISYT30JK)
      LWRK1  = LWORK - KEND1 
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T30_IJ_O (2)')
      END IF
C
      CALL DZERO(WORK(KBACI),NMAAB_CI(ISYT30JK))
C
      CALL SORT_T2_ABJ(WORK(KT2BAL),ISYMJ,J,T2TP,ISYMT2)
C
C-----------------------------
C     Sort (ck|li) = T3OG2(c,i,k,l) as I^K(lci)
C-----------------------------
C
      CALL SORT_INT_JAK_I(WORK(KINTLCI),T3OG2,ISYINT,ISYMK,K)
C
C------------------------------------------
C    Multiply T^J(bal) * I^K(lci) = T^JK(baci)
C------------------------------------------
C
      DO ISYML = 1, NSYM
         ISYBA = MULD2H(ISYT2BAL,ISYML)
         ISYCI = MULD2H(ISYINTLCI,ISYML)
C
         KOFF1 = KT2BAL
     *         + IMAABI(ISYBA,ISYML)
         KOFF2 = KINTLCI
     *         + IMAIAJ(ISYML,ISYCI)
         KOFF3 = KBACI
     *         + IMAAB_CI(ISYBA,ISYCI)
C
         NTOTBA = MAX(NMATAB(ISYBA),1)
         NTOTL  = MAX(NRHF(ISYML),1)
C
         CALL DGEMM('N','N',NMATAB(ISYBA),NT1AM(ISYCI),NRHF(ISYML),
     *              ONE,WORK(KOFF1),NTOTBA,WORK(KOFF2),NTOTL,
     *              ONE,WORK(KOFF3),NTOTBA)
C
      END DO ! ISYML
C
C-------------------------------------------
C     T30JK(bcai) = T30JK(bcai) + T^JK(baci)
C-------------------------------------------
C
C  add_occ(2)
C
      CALL FACBI(T30JK,WORK(KBACI),ISYT30JK)
C
C=================================================
C     Calculate (3)  - sum_l t^(bc)_(jl) (ai|lk) 
C
C                             T^J(bcl) I^K(lai)
C=================================================
C
C-------------------------------
C     Sort T2TP(bjlc) as T^J(bcl)
C-------------------------------
C
      ISYT2BCL = MULD2H(ISYMT2,ISYMJ)
      ISYINTLAI = MULD2H(ISYINT,ISYMK)
C
      KT2BCL = KEND1
      KINTLAI  = KT2BCL + NMAABI(ISYT2BCL)
      KEND1   = KINTLAI + NCKI(ISYINTLAI)
      LWRK1  = LWORK - KEND1 
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T30_IJ_O (3)')
      END IF
C
      CALL SORT_T2_ABJ(WORK(KT2BCL),ISYMJ,J,T2TP,ISYMT2)
C
C-----------------------------
C     Sort (ai|lk) = T3OG2(a,k,i,l) as I^K(lai)
C-----------------------------
C
      CALL SORT_INT_JAI_K(WORK(KINTLAI),T3OG2,ISYINT,ISYMK,K)
C
C------------------------------------------
C    Multiply T^J(bcl) * I^K(lai) = T^JK(bcai)
C------------------------------------------
C
      DO ISYML = 1, NSYM
         ISYBC = MULD2H(ISYT2BCL,ISYML)
         ISYAI = MULD2H(ISYINTLAI,ISYML)
C
         KOFF1 = KT2BCL
     *         + IMAABI(ISYBC,ISYML)
         KOFF2 = KINTLAI
     *         + IMAIAJ(ISYML,ISYAI)
         KOFF3 = IMAAB_CI(ISYBC,ISYAI) + 1
C
         NTOTBC = MAX(NMATAB(ISYBC),1)
         NTOTL  = MAX(NRHF(ISYML),1)
C
C  add_occ(3)
C
         CALL DGEMM('N','N',NMATAB(ISYBC),NT1AM(ISYAI),NRHF(ISYML),
     *              ONE,WORK(KOFF1),NTOTBC,WORK(KOFF2),NTOTL,
     *              ONE,T30JK(KOFF3),NTOTBC)
C
      END DO ! ISYML
C
C=================================================
C     Calculate (4)   - sum_l t^(ac)_(il) (bj|lk)
C
C                             T(lcai) I^JK(bl)
C=================================================
C
C-------------------------------
C     Sort T2TP(clia) as T(lcai)
C-------------------------------
C
      ISYJK = MULD2H(ISYMJ,ISYMK)
      ISYBL = MULD2H(ISYINT,ISYJK)
C
      KT2LCAI = KEND1
      KINTBL  = KT2LCAI + NT2SQ(ISYMT2)
      KEND1   = KINTBL  + NT1AM(ISYBL)
      LWRK1   = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T30_IJ_O (4)')
      END IF
C
      CALL SORT_T2_IABJ(WORK(KT2LCAI),T2TP,ISYMT2)
C
C-----------------------------
C     Sort (bj|lk) = T3OG2(b,k,j,l) as I^JK(bl)
C-----------------------------
C
      CALL SORT_INT_AJ_IK(WORK(KINTBL),T3OG2,ISYINT,ISYMJ,J,ISYMK,K)
C
C------------------------------------------
C    Multiply I^JK(bl) T(lcai) = T^JK(bcai)
C------------------------------------------
C
      DO ISYML = 1, NSYM
            ISYMB = MULD2H(ISYBL,ISYML)
            ISYCAI = MULD2H(ISYMT2,ISYML)
C      
            KOFF1 = KINTBL 
     *            + IT1AM(ISYMB,ISYML)
            KOFF2 = KT2LCAI
     *            + IMAJBAI(ISYML,ISYCAI)
            KOFF3 = IMAAOBCI(ISYMB,ISYCAI) + 1
C
            NTOTB = MAX(NVIR(ISYMB),1)
            NTOTL = MAX(NRHF(ISYML),1)
C
C  add_occ(4)
C
            CALL DGEMM('N','N',NVIR(ISYMB),NMAABI(ISYCAI),NRHF(ISYML),
     *                 ONE,WORK(KOFF1),NTOTB,WORK(KOFF2),NTOTL,
     *                 ONE,T30JK(KOFF3),NTOTB)
C
      END DO ! ISYML
C
C=================================================
C     Calculate (5)  - sum_l t^(ca)_(kl) (bj|li) 
C
C                             T^K(cal) I^J(lbi)
C=================================================
C
C-------------------------------
C     Sort T2TP(ckla) as T^K(cal)
C-------------------------------
C
      ISYT2CAL = MULD2H(ISYMT2,ISYMK)
      ISYINTLBI = MULD2H(ISYINT,ISYMJ)
C
      KT2CAL = KEND1
      KINTLBI  = KT2CAL + NMAABI(ISYT2CAL)
      KCABI   = KINTLBI + NCKI(ISYINTLBI)
      KEND1 = KCABI + NMAAB_CI(ISYT30JK)
      LWRK1  = LWORK - KEND1 
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T30_IJ_O (5)')
      END IF
C
      CALL DZERO(WORK(KCABI),NMAAB_CI(ISYT30JK))
C
      CALL SORT_T2_ABJ(WORK(KT2CAL),ISYMK,K,T2TP,ISYMT2)
C
C-----------------------------
C     Sort (bj|li) = T3OG2(b,i,j,l) as I^J(lbi)
C-----------------------------
C
      CALL SORT_INT_JAK_I(WORK(KINTLBI),T3OG2,ISYINT,ISYMJ,J)
C
C------------------------------------------
C    Multiply T^K(cal) * I^J(lbi) = T^KJ(cabi)
C------------------------------------------
C
      DO ISYML = 1, NSYM
         ISYCA = MULD2H(ISYT2CAL,ISYML)
         ISYBI = MULD2H(ISYINTLBI,ISYML)
C
         KOFF1 = KT2CAL
     *         + IMAABI(ISYCA,ISYML)
         KOFF2 = KINTLBI
     *         + IMAIAJ(ISYML,ISYBI)
         KOFF3 = KCABI
     *         + IMAAB_CI(ISYCA,ISYBI)
C
         NTOTCA = MAX(NMATAB(ISYCA),1)
         NTOTL  = MAX(NRHF(ISYML),1)
C
         CALL DGEMM('N','N',NMATAB(ISYCA),NT1AM(ISYBI),NRHF(ISYML),
     *              ONE,WORK(KOFF1),NTOTCA,WORK(KOFF2),NTOTL,
     *              ONE,WORK(KOFF3),NTOTBA)
C
      END DO ! ISYML
C
C-------------------------------------------
C     T30JK(bcai) = T30JK(bcai) + T^KJ(cabi)
C-------------------------------------------
C
C  add_occ(5)
C
      CALL FBCAI(T30JK,WORK(KCABI),ISYT30JK)
C
C=================================================
C     Calculate (6)  - sum_l t^(cb)_(kl) (ai|lj) 
C
C                             T^K(cbl) I^J(lai)
C=================================================
C
C-------------------------------
C     Sort T2TP(cklb) as T^K(cbl)
C-------------------------------
C
      ISYT2CBL = MULD2H(ISYMT2,ISYMK)
      ISYINTLAI = MULD2H(ISYINT,ISYMJ)
C
      KT2CBL = KEND1
      KINTLAI  = KT2BCL + NMAABI(ISYT2CBL)
      KCBAI   = KINTLAI + NCKI(ISYINTLAI)
      KEND1   = KCBAI   + NMAAB_CI(ISYT30JK)
      LWRK1  = LWORK - KEND1 
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T30_IJ_O (6)')
      END IF
C
      CALL DZERO(WORK(KCBAI),NMAAB_CI(ISYT30JK))
C
      CALL SORT_T2_ABJ(WORK(KT2CBL),ISYMK,K,T2TP,ISYMT2)
C
C-----------------------------
C     Sort (ai|lj) = T3OG2(a,j,i,l) as I^J(lai)
C-----------------------------
C
      CALL SORT_INT_JAI_K(WORK(KINTLAI),T3OG2,ISYINT,ISYMJ,J)
C
C------------------------------------------
C    Multiply T^K(cbl) * I^J(lai) = T^KJ(cbai)
C------------------------------------------
C
      DO ISYML = 1, NSYM
         ISYCB = MULD2H(ISYT2CBL,ISYML)
         ISYAI = MULD2H(ISYINTLAI,ISYML)
C
         KOFF1 = KT2CBL
     *         + IMAABI(ISYCB,ISYML)
         KOFF2 = KINTLAI
     *         + IMAIAJ(ISYML,ISYAI)
         KOFF3 = KCBAI  
     *         + IMAAB_CI(ISYCB,ISYAI) 
C
         NTOTCB = MAX(NMATAB(ISYCB),1)
         NTOTL  = MAX(NRHF(ISYML),1)
C
         CALL DGEMM('N','N',NMATAB(ISYCB),NT1AM(ISYAI),NRHF(ISYML),
     *              ONE,WORK(KOFF1),NTOTCB,WORK(KOFF2),NTOTL,
     *              ONE,WORK(KOFF3),NTOTCB)
C
      END DO ! ISYML
C
C     T30JK(bcai) = T30JK(bcai) + T^JK(cbai)
C
C  add_occ(6)
C
      CALL FBACI(T30JK,WORK(KCBAI),ISYT30JK)
C
C     Scale T30JK with minus sign
C
      CALL DSCAL(NMAABCI(ISYT30JK),-ONE,T30JK,1)
C
      CALL QEXIT('GET_T30_IJ_O')
C
      RETURN
      END 
C
C  /* Deck get_t30_ij_v */
      SUBROUTINE GET_T30_IJ_V(T30JK,ISYT30JK,T2TP,
     *                           ISYMT2,T3VIJG1,
     *                           ISYINT,ISYMJ,J,ISYMK,K,
     *                           WORK,LWORK)

***********************************************************
*    T3VIJG1 : (ck|da) sorted as I(ad|ck) 
*    
*     T30KL sitting as (bcai)
***********************************************************
C
C     T30^(abc)_(iJK) = 
C     P(ai,bj,ck) (sum_d t^(ad)_(ij) (ck|bd) ) + 
C    - P(ai,bj,ck) (sum_l t^(ab)_(il) (ck|lj) )
C
C    In this routine we calculate the first contribution:
C
C    T^JK(bcai) =  P(ai,bj,ck) (sum_d t^(ad)_(ij) (ck|bd) ) 
C
C 1)            +  t^ad_ij (ck|bd)
C
C 2)            +  t^bd_ji (ck!ad)
C
C 3)            +  t^bd_jk (ai!cd)
C
C 4)            +  t^ad_ik (bj!cd)
C
C 5)            +  t^cd_ki (bj!ad)
C
C 6)            +  t^cd_kj (ai!bd)
C
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYT30JK, ISYMT2, ISYINT, ISYMJ, ISYMK, LWORK
      INTEGER ISYMDAI, ISYMBCD, ISYMDA, ISYMBC, ISYMBCA, ISYMDBI  
      INTEGER ISYMACD, ISYMACBI, ISYMDB, ISYMAC, ISYMACB 
      INTEGER ISYMJK, ISYMBD, ISYMCAI, ISYMDCI, ISYMBAD, ISYMBACI 
      INTEGER ISYMDC, ISYMBA, ISYMBAC, ISYMKJ, ISYMCD, ISYMCBAI 
      INTEGER ISYMI, ISYMD, ISYMA, ISYMB, ISYMC
      INTEGER NTOTBC, NTOTD, NTOTAC, NTOTB, NTOTBA, NTOTC 
      INTEGER KDAI, KBCD, KEND1, LWRK1, KDBI, KACD, KACBI, KBD
      INTEGER KDCI, KBAD, KBACI, KCD, KCBAI 
      INTEGER KOFF1, KOFF2, KOFF3
      INTEGER ISYMBAI
C
#if defined (SYS_CRAY)
      REAL T30JK(*), T2TP(*), T3VIJG1(*), WORK(LWORK)
      REAL ONE
#else
      DOUBLE PRECISION T30JK(*), T2TP(*), T3VIJG1(*), WORK(LWORK) 
      DOUBLE PRECISION ONE
#endif
C
      PARAMETER (ONE = 1.0D0)
C
      CALL QENTER('GET_T30_IJ_V')
C
C***************************************************
C 1)               t^ad_ij    *   (ck|bd) 
C***************************************************
C
C t2tp(djia) =   I^J(dai) 
C
C (ck!bd) = I(dbck) =          I^K(bcd) 
C
C T^JK(bcai) = T^JK(bcai) + I^K(bcd)*I^J(dai) 
C
C symmetry and work allocation
C

      ISYMDAI = MULD2H(ISYMT2,ISYMJ)
      ISYMBCD = MULD2H(ISYINT,ISYMK)
      KDAI  = 1
      KBCD  = KDAI  + NMAABI(ISYMDAI)
      KEND1 = KBCD  + NMAABC(ISYMBCD)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T30_IJ_V (1)')
      END IF
C
C  sort t^ad_ij = t2tp(djia) as I^J(dai)
C
      CALL SORT_T2_ABJ(WORK(KDAI),ISYMJ,J,T2TP,ISYMT2)
C
C sort (ck|bd) = T3VIJG1(dbck) as I^K(bcd)
C 


      CALL SORT_INT_BCA(WORK(KBCD),ISYMK,K,T3VIJG1,ISYINT)
C
C T^JK(bcai) = T^JK(bcai) + I^K(bcd)*I^J(dai) 
C
      DO ISYMI = 1,NSYM
         ISYMDA = MULD2H(ISYMDAI,ISYMI)
         DO I =  1,NRHF(ISYMI)
            DO ISYMD = 1,NSYM
               ISYMBC  = MULD2H(ISYMBCD,ISYMD)
               ISYMA   = MULD2H(ISYMDA,ISYMD)
               ISYMBCA = MULD2H(ISYMBC,ISYMA)
               KOFF1   = KBCD + IMAABC(ISYMBC,ISYMD) 
               KOFF2   = KDAI 
     *                + IMAABI(ISYMDA,ISYMI)  
     *                + NMATAB(ISYMDA)*(I-1)
     *                + IMATAB(ISYMD,ISYMA)
               KOFF3   = 1 + IMAABCI(ISYMBCA,ISYMI)
     *                + NMAABC(ISYMBCA)*(I-1)
     *                + IMAABC(ISYMBC,ISYMA)
C  
               NTOTBC = MAX(1,NMATAB(ISYMBC))
               NTOTD  = MAX(1,NVIR(ISYMD))
C
C  add_vir(1)
C
               CALL DGEMM('N','N',NMATAB(ISYMBC),NVIR(ISYMA),
     *                    NVIR(ISYMD),ONE,WORK(KOFF1),NTOTBC,
     *                    WORK(KOFF2),NTOTD, 
     *                    ONE,T30JK(KOFF3),NTOTBC) 
            END DO
         END DO
      END DO

C
C***************************************************
C 2)            +  t^bd_ji (ck!ad)
C***************************************************
C
C t2tp(bjid)  =   I^J(dbi)
C
C (ck!ad) = I(dack)   =  I^K(acd)
C
C T^JK(bcai) = T^JK(bcai) +  I^K(acd) * I^J(dbi)
C
C symmetry and work allocation
C
      ISYMDBI = MULD2H(ISYMT2,ISYMJ)
      ISYMACD = MULD2H(ISYINT,ISYMK)
      ISYMACBI = MULD2H(ISYMACD,ISYMDBI)
      KDBI   = 1
      KACD   = KDBI  +  NMAABI(ISYMDBI)
      KACBI  = KACD  +  NMAABC(ISYMACD)
      KEND1  = KACBI +  NMAABCI(ISYMACBI)
      LWRK1  = LWORK -  KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T30_IJ_V (2)')
      END IF
C
      CALL DZERO(WORK(KACBI),NMAABCI(ISYMACBI))

C
C t2tp(bjid)  =   I^J(dbi)
C
      CALL SORT_T2_BAJ(WORK(KDBI),ISYMJ,J,T2TP,ISYMT2)
C
C (ck!ad) = I(dack)   =  I^K(acd)
C
      CALL SORT_INT_BCA(WORK(KACD),ISYMK,K,T3VIJG1,ISYINT)
C
C
      DO ISYMI = 1,NSYM
         ISYMDB = MULD2H(ISYMDBI,ISYMI)
         DO I =  1,NRHF(ISYMI)
            DO ISYMD = 1,NSYM
               ISYMB   = MULD2H(ISYMDB,ISYMD)
               ISYMAC  = MULD2H(ISYMACD,ISYMD)
               ISYMACB = MULD2H(ISYMAC,ISYMB)
               KOFF1   = KACD
     *                + IMAABC(ISYMAC,ISYMD)
               KOFF2   = KDBI
     *                + IMAABI(ISYMDB,ISYMI)
     *                + NMATAB(ISYMDB)*(I-1)
     *                + IMATAB(ISYMD,ISYMB)
               KOFF3   = KACBI
     *                + IMAABCI(ISYMACB,ISYMI)
     *                + NMAABC(ISYMACB)*(I-1)
     *                + IMAABC(ISYMAC,ISYMB)
C 
               NTOTAC = MAX(1,NMATAB(ISYMAC))
               NTOTD  = MAX(1,NVIR(ISYMD))
C
C  work(acbi) =  I^K(acd) * I^J(dbi)
C

               CALL DGEMM('N','N',NMATAB(ISYMAC),NVIR(ISYMB),
     *                    NVIR(ISYMD),ONE,WORK(KOFF1),NTOTAC,
     *                    WORK(KOFF2),NTOTD,
     *                    ONE,WORK(KOFF3),NTOTAC)

       
            END DO
         END DO
      END DO
C
C T^JK(bcai) = T^JK(bcai) +  work(acbi)
C
C
C  add_vir(2)
C
      CALL FCBAI(T30JK,WORK(KACBI),ISYT30JK)
C
C**************************************************    
C 3)            +  t^bd_jk (ai!cd)
C**************************************************    
C
C t2tp(bjkd)  =   I^JK(bd)
C
C (ai!cd) = I(dcai)   =  I(dcai)
C
C T^JK(bcai) = T^JK(bcai) +  I^JK(bd) * I(dcai) 
C
C symmetry and work allocation
C
      ISYMJK = MULD2H(ISYMJ,ISYMK)
      ISYMBD = MULD2H(ISYMT2,ISYMJK)
      KBD    = 1
      KEND1  = KBD + NMATAB(ISYMBD)
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T30_IJ_V (3)')
      END IF

C
C t2tp(bjkd)  =   I^JK(bd)
C
      CALL SORT_T2_AB(WORK(KBD),ISYMJ,J,ISYMK,K,T2TP,ISYMT2)
C
C T^JK(bcai) = T^JK(bcai) +  I^JK(bd) * I(dcai) 
C
      DO ISYMD = 1,NSYM
         ISYMB   = MULD2H(ISYMBD,ISYMD)
         ISYMCAI = MULD2H(ISYMD,ISYINT)
         KOFF1   = KBD + IMATAB(ISYMB,ISYMD)
         KOFF2   = 1  + IMAAOBCI(ISYMD,ISYMCAI)
         KOFF3   = 1  + IMAAOBCI(ISYMB,ISYMCAI)
C 
         NTOTD  = MAX(1,NVIR(ISYMD))
         NTOTB  = MAX(1,NVIR(ISYMB))
C
C
C  add_vir(3)
C
         CALL DGEMM('N','N',NVIR(ISYMB),NMAABI(ISYMCAI),
     *                    NVIR(ISYMD),ONE,WORK(KOFF1),NTOTB,
     *                    T3VIJG1(KOFF2),NTOTD,
     *                    ONE,T30JK(KOFF3),NTOTD)
      END DO
C
C****************************************************
C 4)            +  t^ad_ik (bj!cd)
C****************************************************
C
C t2tp(dkia) =   I^K(dai)
C
C (bj!cd) = I(dcbj) =      I^J(bcd)
C
C T^JK(bcai) = T^JK(bcai) + I^J(bcd)*I^K(dai)
C
C symmetry and work allocation
C
      ISYMDAI = MULD2H(ISYMT2,ISYMK)
      ISYMBCD = MULD2H(ISYINT,ISYMJ)
      KDAI  = 1
      KBCD  = KDAI  + NMAABI(ISYMDAI)
      KEND1 = KBCD  + NMAABC(ISYMBCD)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T30_IJ_V (4)')
      END IF
C
C t2tp(dkia) =   I^K(dai)
C
      CALL SORT_T2_ABJ(WORK(KDAI),ISYMK,K,T2TP,ISYMT2)
C
C sort (bj|cd) = T3VIJG1(dcbj) as I^J(bcd)
C
      CALL SORT_INT_CBA(WORK(KBCD),ISYMJ,J,T3VIJG1,ISYINT)
C
C T^JK(bcai) = T^JK(bcai) + I^K(bcd)*I^J(dai)
C
      DO ISYMI = 1,NSYM
         ISYMDA = MULD2H(ISYMDAI,ISYMI)
         DO I =  1,NRHF(ISYMI)
            DO ISYMD = 1,NSYM
               ISYMBC  = MULD2H(ISYMBCD,ISYMD)
               ISYMA   = MULD2H(ISYMDA,ISYMD)
               ISYMBCA = MULD2H(ISYMBC,ISYMA)
               KOFF1   = KBCD + IMAABC(ISYMBC,ISYMD)
               KOFF2   = KDAI
     *                + IMAABI(ISYMDA,ISYMI)
     *                + NMATAB(ISYMDA)*(I-1)
     *                + IMATAB(ISYMD,ISYMA)
               KOFF3   = 1 + IMAABCI(ISYMBCA,ISYMI)
     *                + NMAABC(ISYMBCA)*(I-1)
     *                + IMAABC(ISYMBC,ISYMA)
C 
               NTOTBC = MAX(1,NMATAB(ISYMBC))
               NTOTD  = MAX(1,NVIR(ISYMD))
C
C  add_vir(4)
C
               CALL DGEMM('N','N',NMATAB(ISYMBC),NVIR(ISYMA),
     *                    NVIR(ISYMD),ONE,WORK(KOFF1),NTOTBC,
     *                    WORK(KOFF2),NTOTD,
     *                    ONE,T30JK(KOFF3),NTOTBC)
            END DO
         END DO
      END DO
C**************************************************
C 5)            +  t^cd_ki (bj!ad)
C**************************************************
C
C t2tp(ckid)  =   I^K(dci)
C
C (bj!ad) = I(dabj)   =  I^J(bad)
C
C T^JK(bcai) = T^JK(bcai) + I^J(bad) * I^K(dci) 
C
C symmetry and work allocation
C
      ISYMDCI = MULD2H(ISYMT2,ISYMK)
      ISYMBAD = MULD2H(ISYINT,ISYMJ)
      ISYMBACI = MULD2H(ISYMBAD,ISYMDCI)
      KDCI   = 1
      KBAD   = KDCI  +  NMAABI(ISYMDCI)
      KBACI  = KBAD  +  NMAABC(ISYMBAD)
      KEND1  = KBACI +  NMAABCI(ISYMBACI)
      LWRK1  = LWORK -  KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T30_IJ_V (5)')
      END IF
C
      CALL DZERO(WORK(KBACI),NMAABCI(ISYMBACI))
C
C t2tp(ckid)  =   I^K(dci)
C
      CALL SORT_T2_BAJ(WORK(KDCI),ISYMK,K,T2TP,ISYMT2)
C
C (bj!ad) = I(dabj)   =  I^J(bad)
C
      CALL SORT_INT_CBA(WORK(KBAD),ISYMJ,J,T3VIJG1,ISYINT)
C
      DO ISYMI = 1,NSYM
         ISYMDC = MULD2H(ISYMDCI,ISYMI)
         DO I =  1,NRHF(ISYMI)
            DO ISYMD = 1,NSYM
               ISYMC   = MULD2H(ISYMDC,ISYMD)
               ISYMBA  = MULD2H(ISYMBAD,ISYMD)
               ISYMBAC = MULD2H(ISYMBA,ISYMC)
               KOFF1   = KBAD
     *                + IMAABC(ISYMBA,ISYMD)
               KOFF2   = KDCI
     *                + IMAABI(ISYMDC,ISYMI)
     *                + NMATAB(ISYMDC)*(I-1)
     *                + IMATAB(ISYMD,ISYMC)
               KOFF3   = KBACI
     *                + IMAABCI(ISYMBAC,ISYMI)
     *                + NMAABC(ISYMBAC)*(I-1)
     *                + IMAABC(ISYMBA,ISYMC)
C
               NTOTBA = MAX(1,NMATAB(ISYMBA))
               NTOTD  = MAX(1,NVIR(ISYMD))
C
C   WORK(baci) =  I^J(bad) * I^K(dci) 
               CALL DGEMM('N','N',NMATAB(ISYMBA),NVIR(ISYMC),
     *                    NVIR(ISYMD),ONE,WORK(KOFF1),NTOTBA,
     *                    WORK(KOFF2),NTOTD,
     *                    ONE,WORK(KOFF3),NTOTBA)
             END DO
         END DO
      END DO
C
C T^JK(bcai) = T^JK(bcai) + WORK(baci)  
C
C
C  add_vir(5)
C
      CALL FACBI(T30JK,WORK(KBACI),ISYT30JK)
C
C**************************************************
C 6)            +  t^cd_kj (ai!bd)
C****************************************
C
C t2tp(ckjd)  =   I^KJ(cd)
C
C (ai!bd) = I(dbai)   =  I(dbai)
C
C T^JK(bcai) = T^JK(bcai) + I^KJ(cd) * I(dbai) 
C
C symmetry and work allocation
C
      ISYMKJ  = MULD2H(ISYMK,ISYMJ)
      ISYMCD  = MULD2H(ISYMT2,ISYMKJ)
      ISYMCBAI = MULD2H(ISYINT,ISYMCD)
      KCD    = 1
      KCBAI  = KCD   +  NMATAB(ISYMCD)
      KEND1  = KCBAI +  NMAABCI(ISYMCBAI)
      LWRK1  = LWORK -  KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T30_IJ_V (6)')
      END IF
C
      CALL DZERO(WORK(KCBAI),NMAABCI(ISYMCBAI))
C
C t2tp(ckjd)  =   I^KJ(cd)
C
      CALL SORT_T2_AB(WORK(KCD),ISYMK,K,ISYMJ,J,T2TP,ISYMT2)
C
      DO ISYMD = 1,NSYM
         ISYMC   = MULD2H(ISYMCD,ISYMD)
         ISYMBAI = MULD2H(ISYMD,ISYINT)
         KOFF1   = KCD   + IMATAB(ISYMC,ISYMD)
         KOFF2   = 1     + IMAAOBCI(ISYMD,ISYMBAI)
         KOFF3   = KCBAI + IMAAOBCI(ISYMC,ISYMBAI)
C
         NTOTD  = MAX(1,NVIR(ISYMD))
         NTOTC  = MAX(1,NVIR(ISYMC))
C
C  work(cbai) =  I^KJ(cd) * I(dbai) 
         CALL DGEMM('N','N',NVIR(ISYMC),NMAABI(ISYMBAI),
     *                    NVIR(ISYMD),ONE,WORK(KOFF1),NTOTC,
     *                    T3VIJG1(KOFF2),NTOTD,
     *                    ONE,WORK(KOFF3),NTOTC)
      END DO
C
C T^JK(bcai) = T^JK(bcai) + work(cbai) 
C
C
C  add_vir(6)
C
      CALL FBACI(T30JK,WORK(KCBAI),ISYT30JK)
C
      CALL QEXIT('GET_T30_IJ_V')
C
      RETURN   
      END
C  /* Deck fbaci */
      SUBROUTINE FBACI(TABCI,TBACI,ISYMT)
C
C     TABCI = TABCI + TBACI
C
C     F. Pawlowski, Aarhus, Winter 2003
C
      IMPLICIT NONE
#include<ccsdsym.h>
#include<ccorb.h>
#include<priunit.h>
C
      INTEGER ISYMT,ISYMI,ISYABC,ISYBAC,ISYMC,ISYAB,ISYBA,ISYMB,ISYMA
      INTEGER KOFF1,KOFF2
C
#if defined (SYS_CRAY)
      REAL TABCI(*),TBACI(*)
#else
      DOUBLE PRECISION TABCI(*),TBACI(*)
#endif
C
      CALL QENTER('FBACI')
C
      DO  ISYMI = 1,NSYM
        ISYABC = MULD2H(ISYMT,ISYMI)
        ISYBAC = ISYABC
         DO ISYMC =  1,NSYM
            ISYAB = MULD2H(ISYABC,ISYMC)
            ISYBA = MULD2H(ISYBAC,ISYMC)
            DO ISYMB = 1,NSYM
               ISYMA =  MULD2H(ISYAB,ISYMB)
               DO I = 1,NRHF(ISYMI)
                  DO C = 1,NVIR(ISYMC)
                     DO B = 1,NVIR(ISYMB)
                        DO A = 1,NVIR(ISYMA)
                           KOFF1 = IMAABCI(ISYABC,ISYMI)
     *                           + NMAABC(ISYABC)*(I-1)
     *                           + IMAABC(ISYAB,ISYMC)
     *                           + NMATAB(ISYAB)*(C-1)
     *                           + IMATAB(ISYMA,ISYMB)
     *                           + NVIR(ISYMA)*(B-1)
     *                           + A
                           KOFF2 = IMAABCI(ISYBAC,ISYMI)
     *                           + NMAABC(ISYBAC)*(I-1)
     *                           + IMAABC(ISYBA,ISYMC)
     *                           + NMATAB(ISYBA)*(C-1)
     *                           + IMATAB(ISYMB,ISYMA)
     *                           + NVIR(ISYMB)*(A-1)
     *                           + B
C
                           TABCI(KOFF1) = TABCI(KOFF1) + TBACI(KOFF2)
C
                        END DO   
                     END DO   
                  END DO   
               END DO   
            END DO   
         END DO   
      END DO   
C
      CALL QEXIT('FBACI')
C
      RETURN
      END
C  /* Deck fbcai */
      SUBROUTINE FBCAI(TABCI,TBCAI,ISYMT)
C
C     TABCI = TABCI + TBCAI
C
C     F. Pawlowski, Aarhus, Winter 2003
C
      IMPLICIT NONE
#include<ccsdsym.h>
#include<ccorb.h>
#include<priunit.h>
C
      INTEGER ISYMT,ISYMI,ISYABC,ISYBCA,ISYMC,ISYAB,ISYMB,ISYMA,ISYBC
      INTEGER KOFF1,KOFF2
C
#if defined (SYS_CRAY)
      REAL TABCI(*),TBCAI(*)
#else
      DOUBLE PRECISION TABCI(*),TBCAI(*)
#endif
C
      CALL QENTER('FBCAI')
C
      DO  ISYMI = 1,NSYM
        ISYABC = MULD2H(ISYMT,ISYMI)
        ISYBCA = ISYABC
         DO ISYMC =  1,NSYM
            ISYAB = MULD2H(ISYABC,ISYMC)
            DO ISYMB = 1,NSYM
               ISYMA =  MULD2H(ISYAB,ISYMB)
               ISYBC = MULD2H(ISYBCA,ISYMA)
               DO I = 1,NRHF(ISYMI)
                  DO C = 1,NVIR(ISYMC)
                     DO B = 1,NVIR(ISYMB)
                        DO A = 1,NVIR(ISYMA)
                           KOFF1 = IMAABCI(ISYABC,ISYMI)
     *                           + NMAABC(ISYABC)*(I-1)
     *                           + IMAABC(ISYAB,ISYMC)
     *                           + NMATAB(ISYAB)*(C-1)
     *                           + IMATAB(ISYMA,ISYMB)
     *                           + NVIR(ISYMA)*(B-1)
     *                           + A
                           KOFF2 = IMAABCI(ISYBCA,ISYMI)
     *                           + NMAABC(ISYBCA)*(I-1)
     *                           + IMAABC(ISYBC,ISYMA)
     *                           + NMATAB(ISYBC)*(A-1)
     *                           + IMATAB(ISYMB,ISYMC)
     *                           + NVIR(ISYMB)*(C-1)
     *                           + B
C
                           TABCI(KOFF1) = TABCI(KOFF1) + TBCAI(KOFF2)
C
                        END DO   
                     END DO   
                  END DO   
               END DO   
            END DO   
         END DO   
      END DO   
C
      CALL QEXIT('FBCAI')
C
      RETURN
      END
C  /* Deck facbi */
      SUBROUTINE FACBI(TABCI,TACBI,ISYMT)
C
C     TABCI = TABCI + TACBI
C
C     F. Pawlowski, Aarhus, Winter 2003
C
      IMPLICIT NONE
#include<ccsdsym.h>  
#include<ccorb.h>
#include<priunit.h>
C
      INTEGER ISYMT,ISYMI,ISYABC,ISYACB,ISYMC,ISYAB,ISYMB,ISYAC,ISYMA
      INTEGER KOFF1,KOFF2 
C
#if defined (SYS_CRAY)
      REAL TABCI(*),TACBI(*)
#else
      DOUBLE PRECISION TABCI(*),TACBI(*)
#endif
C
      CALL QENTER('FACBI')

      DO  ISYMI = 1,NSYM
        ISYABC = MULD2H(ISYMT,ISYMI)
        ISYACB = ISYABC
         DO ISYMC =  1,NSYM
            ISYAB = MULD2H(ISYABC,ISYMC)
            DO ISYMB = 1,NSYM
               ISYAC = MULD2H(ISYACB,ISYMB)
               ISYMA =  MULD2H(ISYAB,ISYMB)
               DO I = 1,NRHF(ISYMI)
                  DO C = 1,NVIR(ISYMC)
                     DO B = 1,NVIR(ISYMB)
                        DO A = 1,NVIR(ISYMA)
                           KOFF1 = IMAABCI(ISYABC,ISYMI)
     *                           + NMAABC(ISYABC)*(I-1)
     *                           + IMAABC(ISYAB,ISYMC)
     *                           + NMATAB(ISYAB)*(C-1)
     *                           + IMATAB(ISYMA,ISYMB)
     *                           + NVIR(ISYMA)*(B-1)
     *                           + A
                           KOFF2 = IMAABCI(ISYACB,ISYMI)
     *                           + NMAABC(ISYACB)*(I-1)
     *                           + IMAABC(ISYAC,ISYMB)
     *                           + NMATAB(ISYAC)*(B-1)
     *                           + IMATAB(ISYMA,ISYMC)
     *                           + NVIR(ISYMA)*(C-1)
     *                           + A
C
                           TABCI(KOFF1) = TABCI(KOFF1) + TACBI(KOFF2)
C
                        END DO   
                     END DO   
                  END DO   
               END DO   
            END DO   
         END DO   
      END DO   
C
      CALL QEXIT('FACBI')
C
      RETURN
      END
C  /* Deck fcabi */
      SUBROUTINE FCABI(TABCI,TCABI,ISYMT)
C
C     TABCI = TABCI + TCABI
C
C     F. Pawlowski, Aarhus, Winter 2003
C
      IMPLICIT NONE
#include<ccsdsym.h>
#include<ccorb.h>
#include<priunit.h>
C
      INTEGER ISYMT,ISYMI,ISYABC,ISYCAB,ISYMC,ISYAB,ISYMB,ISYCA,ISYMA
      INTEGER KOFF1,KOFF2
C
#if defined (SYS_CRAY)
      REAL TABCI(*),TCABI(*)
#else
      DOUBLE PRECISION TABCI(*),TCABI(*)
#endif
C
      CALL QENTER('FCABI')
C
      DO  ISYMI = 1,NSYM
        ISYABC = MULD2H(ISYMT,ISYMI)
        ISYCAB = ISYABC
         DO ISYMC =  1,NSYM
            ISYAB = MULD2H(ISYABC,ISYMC)
            DO ISYMB = 1,NSYM
               ISYCA = MULD2H(ISYCAB,ISYMB)
               ISYMA =  MULD2H(ISYAB,ISYMB)
               DO I = 1,NRHF(ISYMI)
                  DO C = 1,NVIR(ISYMC)
                     DO B = 1,NVIR(ISYMB)
                        DO A = 1,NVIR(ISYMA)
                           KOFF1 = IMAABCI(ISYABC,ISYMI)
     *                           + NMAABC(ISYABC)*(I-1)
     *                           + IMAABC(ISYAB,ISYMC)
     *                           + NMATAB(ISYAB)*(C-1)
     *                           + IMATAB(ISYMA,ISYMB)
     *                           + NVIR(ISYMA)*(B-1)
     *                           + A
                           KOFF2 = IMAABCI(ISYCAB,ISYMI)
     *                           + NMAABC(ISYCAB)*(I-1)
     *                           + IMAABC(ISYCA,ISYMB)
     *                           + NMATAB(ISYCA)*(B-1)
     *                           + IMATAB(ISYMC,ISYMA)
     *                           + NVIR(ISYMC)*(A-1)
     *                           + C
C
                           TABCI(KOFF1) = TABCI(KOFF1) + TCABI(KOFF2)
C
                        END DO   
                     END DO   
                  END DO   
               END DO   
            END DO   
         END DO   
      END DO   
C
      CALL QEXIT('FCABI')
C
      RETURN
      END
C  /* Deck fcbai */
      SUBROUTINE FCBAI(TABCI,TCBAI,ISYMT)
C
C     TABCI = TABCI + TCBAI
C
C     F. Pawlowski, Aarhus, Winter 2003
C
      IMPLICIT NONE
#include<ccsdsym.h>  
#include<ccorb.h>
#include<priunit.h>
C
      INTEGER ISYMT,ISYMI,ISYABC,ISYCBA,ISYMC,ISYAB,ISYMB,ISYMA,ISYCB
      INTEGER KOFF1,KOFF2
C
#if defined (SYS_CRAY)
      REAL TABCI(*),TCBAI(*)
#else
      DOUBLE PRECISION TABCI(*),TCBAI(*)
#endif
C
      CALL QENTER('FCBAI')
C
      DO  ISYMI = 1,NSYM
        ISYABC = MULD2H(ISYMT,ISYMI)
        ISYCBA = ISYABC
         DO ISYMC =  1,NSYM
            ISYAB = MULD2H(ISYABC,ISYMC)
            DO ISYMB = 1,NSYM
               ISYMA =  MULD2H(ISYAB,ISYMB)
               ISYCB = MULD2H(ISYCBA,ISYMA)
               DO I = 1,NRHF(ISYMI)
                  DO C = 1,NVIR(ISYMC)
                     DO B = 1,NVIR(ISYMB)
                        DO A = 1,NVIR(ISYMA)
                           KOFF1 = IMAABCI(ISYABC,ISYMI)
     *                           + NMAABC(ISYABC)*(I-1)
     *                           + IMAABC(ISYAB,ISYMC)
     *                           + NMATAB(ISYAB)*(C-1)
     *                           + IMATAB(ISYMA,ISYMB)
     *                           + NVIR(ISYMA)*(B-1)
     *                           + A
                           KOFF2 = IMAABCI(ISYCBA,ISYMI)
     *                           + NMAABC(ISYCBA)*(I-1)
     *                           + IMAABC(ISYCB,ISYMA)
     *                           + NMATAB(ISYCB)*(A-1)
     *                           + IMATAB(ISYMC,ISYMB)
     *                           + NVIR(ISYMC)*(B-1)
     *                           + C
C
                           TABCI(KOFF1) = TABCI(KOFF1) + TCBAI(KOFF2)
C
                        END DO   
                     END DO   
                  END DO   
               END DO   
            END DO   
         END DO   
      END DO   
C
      CALL QEXIT('FCBAI')
C
      RETURN
      END
C  /* Deck fab_ci */
      SUBROUTINE FAB_CI(TAB_CI,TABCI,ISYMT)
C
C     TAB_CI = TAB_CI + TABCI
C  
C     which means: T(a,b,c,i) --> T(ab,ci)
C
C     F. Pawlowski, Aarhus, Winter 2003
C
      IMPLICIT NONE
#include<ccsdsym.h>
#include<ccorb.h>
#include<priunit.h>
C
      INTEGER ISYMT,ISYMI,ISYABC,ISYMC,ISYAB,ISYCI,ISYMB,ISYMA
      INTEGER NAB,NCI
      INTEGER KOFF1,KOFF2
C
#if defined (SYS_CRAY)
      REAL TABCI(*),TAB_CI(*)
#else
      DOUBLE PRECISION TABCI(*),TAB_CI(*)
#endif
C
      CALL QENTER('FAB_CI')
C
      DO  ISYMI = 1,NSYM
        ISYABC = MULD2H(ISYMT,ISYMI)
         DO ISYMC =  1,NSYM
            ISYAB = MULD2H(ISYABC,ISYMC)
            ISYCI = MULD2H(ISYMC,ISYMI)
            DO ISYMB = 1,NSYM
               ISYMA =  MULD2H(ISYAB,ISYMB)
               DO I = 1,NRHF(ISYMI)
                  DO C = 1,NVIR(ISYMC)
                     NCI = IT1AM(ISYMC,ISYMI)
     *                   + NVIR(ISYMC)*(I-1)
     *                   + C
                     DO B = 1,NVIR(ISYMB)
                        DO A = 1,NVIR(ISYMA)
                           NAB = IMATAB(ISYMA,ISYMB)
     *                         + NVIR(ISYMA)*(B-1)
     *                         + A
                           KOFF1 = IMAABCI(ISYABC,ISYMI)
     *                           + NMAABC(ISYABC)*(I-1)
     *                           + IMAABC(ISYAB,ISYMC)
     *                           + NMATAB(ISYAB)*(C-1)
     *                           + IMATAB(ISYMA,ISYMB)
     *                           + NVIR(ISYMA)*(B-1)
     *                           + A
                           KOFF2 = IMAAB_CI(ISYAB,ISYCI)
     *                           + NMATAB(ISYAB)*(NCI-1)
     *                           + NAB
C
                           TAB_CI(KOFF2) = TAB_CI(KOFF2) + TABCI(KOFF1)
C
                        END DO   
                     END DO   
                  END DO   
               END DO   
            END DO   
         END DO   
      END DO   
C
      CALL QEXIT('FAB_CI')
C
      RETURN
      END
C  /* Deck sort_int_aj_ik */
      SUBROUTINE SORT_INT_AJ_IK(XINTAJ,XINTAIJK,ISYINT,ISYMI,I,ISYMK,K)
C
C-----------------------------
C     Sort (ai|jk) = XINTAIJK(a,k,i,j) as I^IK(aj)
C-----------------------------
C
C     F. Pawlowski, Aarhus, Winter 2003
C
      IMPLICIT NONE
#include<ccsdsym.h>
#include<ccorb.h>
#include<priunit.h>
C
      INTEGER ISYINT,ISYMI,ISYMK,ISYAKI,ISYMJ,ISYAK,ISYMA
      INTEGER KOFF1,KOFF2
C
#if defined (SYS_CRAY)
      REAL XINTAJ(*),XINTAIJK(*)
#else
      DOUBLE PRECISION XINTAJ(*),XINTAIJK(*)
#endif
C
      CALL QENTER('SORT_INT_AJ_IK')
C
      DO ISYMJ = 1,NSYM
         ISYAKI = MULD2H(ISYINT,ISYMJ)
         ISYAK = MULD2H(ISYAKI,ISYMI)
         ISYMA = MULD2H(ISYAK,ISYMK)
         DO J = 1,NRHF(ISYMJ)
            DO A = 1,NVIR(ISYMA)

               KOFF1 = ISAIKJ(ISYAKI,ISYMJ)
     *               + NCKI(ISYAKI)*(J-1)
     *               + ICKI(ISYAK,ISYMI)
     *               + NT1AM(ISYAK)*(I-1)
     *               + IT1AM(ISYMA,ISYMK)
     *               + NVIR(ISYMA)*(K-1)
     *               + A
               KOFF2 = IT1AM(ISYMA,ISYMJ)
     *               + NVIR(ISYMA)*(J-1)
     *               + A
               XINTAJ(KOFF2) = XINTAIJK(KOFF1)
            END DO
         END DO
      END DO
      CALL QEXIT('SORT_INT_AJ_IK')
C
      RETURN
      END
C  /*  Deck sort_int_jak_i */
      SUBROUTINE SORT_INT_JAK_I(XINTJAK,XINTAIJK,ISYINT,ISYMI,I)
C
C-----------------------------
C     Sort (ai|jk) = T3OG2(a,k,i,j) as I^I(jak)
C-----------------------------
C
C     F. Pawlowski, Aarhus, Winter 2003
C
      IMPLICIT NONE
#include<ccsdsym.h>
#include<ccorb.h>
#include<priunit.h>
C
      INTEGER ISYINT,ISYMI,ISYAKI,ISYMJ,ISYAK,ISYMK,ISYMA,ISYJA
      INTEGER KOFF1,KOFF2
C
#if defined (SYS_CRAY)
      REAL XINTJAK(*),XINTAIJK(*)
#else
      DOUBLE PRECISION XINTJAK(*),XINTAIJK(*)
#endif
C
      CALL QENTER('SORT_INT_JAK_I')
C
      DO ISYMJ = 1,NSYM
         ISYAKI = MULD2H(ISYINT,ISYMJ)
         ISYAK = MULD2H(ISYAKI,ISYMI)
         DO ISYMK = 1,NSYM
            ISYMA = MULD2H(ISYAK,ISYMK)
            ISYJA = MULD2H(ISYMJ,ISYMA)
            DO J = 1,NRHF(ISYMJ)
               DO K = 1,NRHF(ISYMK)
                  DO A = 1,NVIR(ISYMA)
                     KOFF1 = ISAIKJ(ISYAKI,ISYMJ)
     *                     + NCKI(ISYAKI)*(J-1)
     *                     + ICKI(ISYAK,ISYMI)
     *                     + NT1AM(ISYAK)*(I-1)
     *                     + IT1AM(ISYMA,ISYMK)
     *                     + NVIR(ISYMA)*(K-1)
     *                     + A
                     KOFF2 = ICKI(ISYJA,ISYMK)
     *                     + NT1AM(ISYJA)*(K-1)
     *                     + IT1AMT(ISYMJ,ISYMA)
     *                     + NRHF(ISYMJ)*(A-1)
     *                     + J
C
                     XINTJAK(KOFF2) = XINTAIJK(KOFF1)
                  END DO
               END DO
            END DO
         END DO
      END DO
      CALL QEXIT('SORT_INT_JAK_I')
C
      RETURN
      END
C  /* Deck sort_int_jai_k */
      SUBROUTINE SORT_INT_JAI_K(XINTJAI,XINTAIJK,ISYINT,ISYMK,K)
C
C-----------------------------
C     Sort (ai|jk) = T3OG2(a,k,i,j) as I^K(jai)
C-----------------------------
C
C     F. Pawlowski, Aarhus, Winter 2003
C
      IMPLICIT NONE
#include<ccsdsym.h>
#include<ccorb.h>
#include<priunit.h>
C
      INTEGER ISYINT,ISYMK,ISYMI,ISYAKI,ISYMJ,ISYAK,ISYMA,ISYJA
      INTEGER KOFF1,KOFF2
C
#if defined (SYS_CRAY)
      REAL XINTJAI(*),XINTAIJK(*)
#else
      DOUBLE PRECISION XINTJAI(*),XINTAIJK(*)
#endif
C
      CALL QENTER('SORT_INT_JAI_K')
C
      DO ISYMJ = 1,NSYM
         ISYAKI = MULD2H(ISYINT,ISYMJ)
         DO ISYMI = 1,NSYM
            ISYAK = MULD2H(ISYAKI,ISYMI)
            ISYMA = MULD2H(ISYAK,ISYMK)
            ISYJA = MULD2H(ISYMJ,ISYMA)
            DO J = 1,NRHF(ISYMJ)
               DO I = 1,NRHF(ISYMI)
                  DO A = 1,NVIR(ISYMA)
                     KOFF1 = ISAIKJ(ISYAKI,ISYMJ)
     *                     + NCKI(ISYAKI)*(J-1)
     *                     + ICKI(ISYAK,ISYMI)
     *                     + NT1AM(ISYAK)*(I-1)
     *                     + IT1AM(ISYMA,ISYMK)
     *                     + NVIR(ISYMA)*(K-1)
     *                     + A
                     KOFF2 = ICKI(ISYJA,ISYMI)
     *                     + NT1AM(ISYJA)*(I-1)
     *                     + IT1AMT(ISYMJ,ISYMA)
     *                     + NRHF(ISYMJ)*(A-1)
     *                     + J
C
                     XINTJAI(KOFF2) = XINTAIJK(KOFF1)
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      CALL QEXIT('SORT_INT_JAI_K')
C
      RETURN
      END
C  /* Deck sort_int_bac */
      SUBROUTINE SORT_INT_BAC(XBAC,ISYMK,K,T3VIJG1,ISYINT)
C
C I^K(BAC) = T3VIJG1(abck) 
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYMK, ISYINT
      INTEGER ISYMBAC, ISYMABC, ISYMA, ISYMBA, ISYMC, ISYMB, ISYMAB
      INTEGER KOFF1, KOFF2
C
#if defined (SYS_CRAY)
      REAL XBAC(*), T3VIJG1(*)
#else
      DOUBLE PRECISION XBAC(*), T3VIJG1(*)
#endif
C
      CALL QENTER('SORT_INT_BAC')
C
      ISYMBAC = MULD2H(ISYINT,ISYMK)
      ISYMABC = ISYMBAC
      DO ISYMC = 1,NSYM
         ISYMBA = MULD2H(ISYMBAC,ISYMC)
         ISYMAB = ISYMBA
         DO ISYMA = 1,NSYM
            ISYMB = MULD2H(ISYMBA,ISYMA)
            ISYMAB = MULD2H(ISYMABC,ISYMC)
            DO A = 1,NVIR(ISYMA)
               DO B = 1,NVIR(ISYMB)
                  DO C = 1,NVIR(ISYMC)
                     KOFF2 = IMAABCI(ISYMABC,ISYMK)
     *                     + NMAABC(ISYMABC)*(K-1)
     *                     + IMAABC(ISYMAB,ISYMC)
     *                     + NMATAB(ISYMAB)*(C-1)
     *                     + IMATAB(ISYMA,ISYMB)
     *                     + NVIR(ISYMA)*(B-1)
     *                     + A
                     KOFF1 = IMAABC(ISYMBA,ISYMC)
     *                     + NMATAB(ISYMBA)*(C-1)
     *                     + IMATAB(ISYMB,ISYMA)
     *                     + NVIR(ISYMB)*(A-1)
     *                     + B
                     XBAC(KOFF1) = T3VIJG1(KOFF2)
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      CALL QEXIT('SORT_INT_BAC')
C
      RETURN
      END
C  /* Deck sort_int_cab */
      SUBROUTINE SORT_INT_CAB(XCAB,ISYMK,K,T3VIJG1,ISYINT)
C
C I^K(CAB) = T3VIJG1(abck) 
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYMK, ISYINT
      INTEGER ISYMCAB, ISYMABC, ISYMA, ISYMCA, ISYMC, ISYMB, ISYMAB
      INTEGER KOFF1, KOFF2
C
#if defined (SYS_CRAY)
      REAL XCAB(*), T3VIJG1(*)
#else
      DOUBLE PRECISION XCAB(*), T3VIJG1(*)
#endif
C
      CALL QENTER('SORT_INT_CAB')
C
      ISYMCAB = MULD2H(ISYINT,ISYMK)
      ISYMABC = ISYMCAB
      DO ISYMB = 1,NSYM
         ISYMCA = MULD2H(ISYMCAB,ISYMB)
         DO ISYMA = 1,NSYM
            ISYMC = MULD2H(ISYMCA,ISYMA)
            ISYMAB = MULD2H(ISYMABC,ISYMC)
            DO A = 1,NVIR(ISYMA)
               DO B = 1,NVIR(ISYMB)
                  DO C = 1,NVIR(ISYMC)
                     KOFF2 = IMAABCI(ISYMABC,ISYMK)
     *                     + NMAABC(ISYMABC)*(K-1)
     *                     + IMAABC(ISYMAB,ISYMC)
     *                     + NMATAB(ISYMAB)*(C-1)
     *                     + IMATAB(ISYMA,ISYMB)
     *                     + NVIR(ISYMA)*(B-1)
     *                     + A
                     KOFF1 = IMAABC(ISYMCA,ISYMB)
     *                     + NMATAB(ISYMCA)*(B-1)
     *                     + IMATAB(ISYMC,ISYMA)
     *                     + NVIR(ISYMC)*(A-1)
     *                     + C
                     XCAB(KOFF1) = T3VIJG1(KOFF2)
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      CALL QEXIT('SORT_INT_CAB')
C
      RETURN
      END
C  /* Deck sort_int_cba */
      SUBROUTINE SORT_INT_CBA(XCBA,ISYMK,K,T3VIJG1,ISYINT)
C
C I^K(cba) = T3VIJG1(abck) 
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYMK, ISYINT
      INTEGER ISYMCBA, ISYMABC, ISYMA, ISYMCB, ISYMC, ISYMB, ISYMAB
      INTEGER KOFF1, KOFF2
C
#if defined (SYS_CRAY)
      REAL XCBA(*), T3VIJG1(*)
#else
      DOUBLE PRECISION XCBA(*), T3VIJG1(*)
#endif
C
      CALL QENTER('SORT_INT_CBA')
C
      ISYMCBA = MULD2H(ISYINT,ISYMK)
      ISYMABC = ISYMCBA
      DO ISYMA = 1,NSYM
         ISYMCB = MULD2H(ISYMCBA,ISYMA)
         DO ISYMB = 1,NSYM
            ISYMC = MULD2H(ISYMCB,ISYMB)
            ISYMAB = MULD2H(ISYMABC,ISYMC)
            DO A = 1,NVIR(ISYMA)
               DO B = 1,NVIR(ISYMB)
                  DO C = 1,NVIR(ISYMC)
                     KOFF2 = IMAABCI(ISYMABC,ISYMK)
     *                     + NMAABC(ISYMABC)*(K-1)
     *                     + IMAABC(ISYMAB,ISYMC)
     *                     + NMATAB(ISYMAB)*(C-1)
     *                     + IMATAB(ISYMA,ISYMB)
     *                     + NVIR(ISYMA)*(B-1)
     *                     + A
                     KOFF1 = IMAABC(ISYMCB,ISYMA)
     *                     + NMATAB(ISYMCB)*(A-1)
     *                     + IMATAB(ISYMC,ISYMB)
     *                     + NVIR(ISYMC)*(B-1)
     *                     + C
                     XCBA(KOFF1) = T3VIJG1(KOFF2)
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      CALL QEXIT('SORT_INT_CBA')
C
      RETURN
      END
C  /* Deck sort_int_abc */
      SUBROUTINE SORT_INT_ABC(XABC,ISYMK,K,T3VIJG1,ISYINT)
C
C I^K(ABC) = T3VIJG1(abck) 
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYMK, ISYINT
      INTEGER ISYMABC, ISYMA, ISYMC, ISYMB, ISYMAB
      INTEGER KOFF1, KOFF2
C
#if defined (SYS_CRAY)
      REAL XABC(*), T3VIJG1(*)
#else
      DOUBLE PRECISION XABC(*), T3VIJG1(*)
#endif
C
      CALL QENTER('SORT_INT_ABC')
C
      ISYMABC = MULD2H(ISYINT,ISYMK)
      DO ISYMC = 1,NSYM
         ISYMAB = MULD2H(ISYMABC,ISYMC)
         DO ISYMB = 1,NSYM
            ISYMA = MULD2H(ISYMAB,ISYMB)
            DO A = 1,NVIR(ISYMA)
               DO B = 1,NVIR(ISYMB)
                  DO C = 1,NVIR(ISYMC)
                     KOFF2 = IMAABCI(ISYMABC,ISYMK)
     *                     + NMAABC(ISYMABC)*(K-1)
     *                     + IMAABC(ISYMAB,ISYMC)
     *                     + NMATAB(ISYMAB)*(C-1)
     *                     + IMATAB(ISYMA,ISYMB)
     *                     + NVIR(ISYMA)*(B-1)
     *                     + A
                     KOFF1 = IMAABC(ISYMAB,ISYMC)
     *                     + NMATAB(ISYMAB)*(C-1)
     *                     + IMATAB(ISYMA,ISYMB)
     *                     + NVIR(ISYMA)*(B-1)
     *                     + A
                     XABC(KOFF1) = T3VIJG1(KOFF2)
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      CALL QEXIT('SORT_INT_ABC')
C
      RETURN
      END
C  /* Deck sort_int_bca */
      SUBROUTINE SORT_INT_BCA(XBCA,ISYMK,K,T3VIJG1,ISYINT)
C
C I^K(BCA) = T3VIJG1(abck) 
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYMK, ISYINT
      INTEGER ISYMABC, ISYMBCA, ISYMBC, ISYMA, ISYMC, ISYMB, ISYMAB
      INTEGER KOFF1, KOFF2
C
#if defined (SYS_CRAY)
      REAL XBCA(*), T3VIJG1(*)
#else
      DOUBLE PRECISION XBCA(*), T3VIJG1(*)
#endif
C
      CALL QENTER('SORT_INT_BCA')
C
      ISYMBCA = MULD2H(ISYINT,ISYMK)
      ISYMABC = ISYMBCA
      DO ISYMA = 1,NSYM
         ISYMBC = MULD2H(ISYMBCA,ISYMA)
         DO ISYMC = 1,NSYM
            ISYMB = MULD2H(ISYMBC,ISYMC)
            ISYMAB = MULD2H(ISYMABC,ISYMC)
            DO A = 1,NVIR(ISYMA)
               DO B = 1,NVIR(ISYMB)
                  DO C = 1,NVIR(ISYMC)
                     KOFF2 = IMAABCI(ISYMABC,ISYMK)
     *                     + NMAABC(ISYMABC)*(K-1)
     *                     + IMAABC(ISYMAB,ISYMC)
     *                     + NMATAB(ISYMAB)*(C-1)
     *                     + IMATAB(ISYMA,ISYMB)
     *                     + NVIR(ISYMA)*(B-1)
     *                     + A
                     KOFF1 = IMAABC(ISYMBC,ISYMA)
     *                     + NMATAB(ISYMBC)*(A-1)
     *                     + IMATAB(ISYMB,ISYMC)
     *                     + NVIR(ISYMB)*(C-1)
     *                     + B
                     XBCA(KOFF1) = T3VIJG1(KOFF2)
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      CALL QEXIT('SORT_INT_BCA')
C
      RETURN
      END
C  /* Deck sort_int_acb */
      SUBROUTINE SORT_INT_ACB(XACB,ISYMK,K,T3VIJG1,ISYINT)
C
C I^K(acb) = T3VIJG1(abck) 
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYMK, ISYINT
      INTEGER ISYMACB, ISYMABC, ISYMA, ISYMAC, ISYMC, ISYMB, ISYMAB
      INTEGER KOFF1, KOFF2
C
#if defined (SYS_CRAY)
      REAL XACB(*), T3VIJG1(*)
#else
      DOUBLE PRECISION XACB(*), T3VIJG1(*)
#endif
C
      CALL QENTER('SORT_INT_ACB')
C
      ISYMACB = MULD2H(ISYINT,ISYMK)
      ISYMABC = ISYMACB
      DO ISYMB = 1,NSYM
         ISYMAC = MULD2H(ISYMACB,ISYMB)
         DO ISYMC = 1,NSYM
            ISYMA = MULD2H(ISYMAC,ISYMC)
            ISYMAB = MULD2H(ISYMABC,ISYMC)
            DO A = 1,NVIR(ISYMA)
               DO B = 1,NVIR(ISYMB)
                  DO C = 1,NVIR(ISYMC)
                     KOFF2 = IMAABCI(ISYMABC,ISYMK)
     *                     + NMAABC(ISYMABC)*(K-1)
     *                     + IMAABC(ISYMAB,ISYMC)
     *                     + NMATAB(ISYMAB)*(C-1)
     *                     + IMATAB(ISYMA,ISYMB)
     *                     + NVIR(ISYMA)*(B-1)
     *                     + A
                     KOFF1 = IMAABC(ISYMAC,ISYMB)
     *                     + NMATAB(ISYMAC)*(B-1)
     *                     + IMATAB(ISYMA,ISYMC)
     *                     + NVIR(ISYMA)*(C-1)
     *                     + A
                     XACB(KOFF1) = T3VIJG1(KOFF2)
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      CALL QEXIT('SORT_INT_ACB')
C
      RETURN
      END
C  /* Deck sort_t2_iabj */
      SUBROUTINE SORT_T2_IABJ(T2IABJ,T2TP,ISYMT2)
C
C-------------------------------
C     Sort T2TP(aijb) as T2(iabj)
C-------------------------------
C
C     F. Pawlowski, Aarhus, Winter 2003
C
      IMPLICIT NONE
#include<ccsdsym.h>
#include<ccorb.h>
#include<priunit.h>
C
      INTEGER ISYMT2,ISYMB,ISYAIJ,ISYMJ,ISYAI,ISYIA,ISYIAB,ISYMI,ISYMA
      INTEGER KOFF1,KOFF2
C
#if defined (SYS_CRAY)
      REAL T2IABJ(*),T2TP(*)
#else
      DOUBLE PRECISION T2IABJ(*),T2TP(*)
#endif
C
      CALL QENTER('SORT_T2_IABJ')
C
C     Sort T2TP(aijb) as T2(iabj)
      DO ISYMB = 1,NSYM
         ISYAIJ = MULD2H(ISYMT2,ISYMB)
         DO ISYMJ = 1,NSYM
            ISYAI = MULD2H(ISYAIJ,ISYMJ)
            ISYIA = ISYAI
            ISYIAB = MULD2H(ISYIA,ISYMB)
            DO ISYMI = 1,NSYM
               ISYMA = MULD2H(ISYAI,ISYMI)
               DO B = 1,NVIR(ISYMB)
                  DO J = 1,NRHF(ISYMJ)
                     DO I = 1,NRHF(ISYMI)
                        DO A = 1,NVIR(ISYMA)
                           KOFF1 = IT2SP(ISYAIJ,ISYMB)
     *                           + NCKI(ISYAIJ)*(B-1)
     *                           + ICKI(ISYAI,ISYMJ)
     *                           + NT1AM(ISYAI)*(J-1)
     *                           + IT1AM(ISYMA,ISYMI)
     *                           + NVIR(ISYMA)*(I-1)
     *                           + A
                           KOFF2 = IMAJBAIT(ISYIAB,ISYMJ)
     *                           + NMAIAB(ISYIAB)*(J-1)
     *                           + IMAIAB(ISYIA,ISYMB)
     *                           + NT1AM(ISYIA)*(B-1)
     *                           + IT1AMT(ISYMI,ISYMA)
     *                           + NRHF(ISYMI)*(A-1)
     *                           + I
C
                           T2IABJ(KOFF2) = T2TP(KOFF1)
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      CALL QEXIT('SORT_T2_IABJ')
C
      RETURN
      END
C  /* Deck sort_t2_abj */
      SUBROUTINE SORT_T2_ABJ(XABJ,ISYMI,I,T2TP,ISYMT2) 
C
C     t2tp(aijb) as I^I(abj)
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYMI, ISYMT2
      INTEGER ISYMABJ, ISYMJ, ISYMAB, ISYMB, ISYMA, ISYMAIJ, ISYMAI
      INTEGER KOFF1, KOFF2
#if defined (SYS_CRAY)
      REAL XABJ(*), T2TP(*) 
#else
      DOUBLE PRECISION XABJ(*), T2TP(*) 
#endif
C
      CALL QENTER('SORT_T2_ABJ')
C
      ISYMABJ = MULD2H(ISYMT2,ISYMI)
      DO ISYMJ = 1,NSYM
         ISYMAB = MULD2H(ISYMABJ,ISYMJ)
         DO ISYMB = 1,NSYM 
            ISYMA = MULD2H(ISYMAB,ISYMB)
            ISYMAIJ = MULD2H(ISYMT2,ISYMB)
            ISYMAI  = MULD2H(ISYMAIJ,ISYMJ)
            DO J = 1,NRHF(ISYMJ)
               DO B = 1,NVIR(ISYMB)
                  DO A = 1,NVIR(ISYMA)
                     KOFF1 = IMAABI(ISYMAB,ISYMJ)
     *                     + NMATAB(ISYMAB)*(J-1)
     *                     + IMATAB(ISYMA,ISYMB)
     *                     + NVIR(ISYMA)*(B-1)
     *                     + A
C
                     KOFF2 = IT2SP(ISYMAIJ,ISYMB)
     *                     + NCKI(ISYMAIJ)*(B-1)
     *                     + ICKI(ISYMAI,ISYMJ)
     *                     + NT1AM(ISYMAI)*(J-1)
     *                     + IT1AM(ISYMA,ISYMI)
     *                     + NVIR(ISYMA)*(I-1)
     *                     + A
C
                     XABJ(KOFF1) = T2TP(KOFF2)
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      CALL QEXIT('SORT_T2_ABJ')
C
      RETURN
      END
C  /* Deck sort_t2_abi */
      SUBROUTINE SORT_T2_ABI(XABI,ISYMJ,J,T2TP,ISYMT2) 
C
C     t2tp(aijb) as I^J(abi)
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYMI, ISYMT2
      INTEGER ISYMABI, ISYMJ, ISYMAB, ISYMB, ISYMA, ISYMAIJ, ISYMAI
      INTEGER KOFF1, KOFF2
#if defined (SYS_CRAY)
      REAL XABI(*), T2TP(*) 
#else
      DOUBLE PRECISION XABI(*), T2TP(*) 
#endif
C
      CALL QENTER('SORT_T2_ABI')
C
      ISYMABI = MULD2H(ISYMT2,ISYMJ)
      DO ISYMI = 1,NSYM
         ISYMAB = MULD2H(ISYMABI,ISYMI)
         DO ISYMB = 1,NSYM 
            ISYMA = MULD2H(ISYMAB,ISYMB)
            ISYMAI  = MULD2H(ISYMA,ISYMI)
            ISYMAIJ = MULD2H(ISYMAI,ISYMJ)
            DO I = 1,NRHF(ISYMI)
               DO B = 1,NVIR(ISYMB)
                  DO A = 1,NVIR(ISYMA)
                     KOFF1 = IMAABI(ISYMAB,ISYMI)
     *                     + NMATAB(ISYMAB)*(I-1)
     *                     + IMATAB(ISYMA,ISYMB)
     *                     + NVIR(ISYMA)*(B-1)
     *                     + A
C
                     KOFF2 = IT2SP(ISYMAIJ,ISYMB)
     *                     + NCKI(ISYMAIJ)*(B-1)
     *                     + ICKI(ISYMAI,ISYMJ)
     *                     + NT1AM(ISYMAI)*(J-1)
     *                     + IT1AM(ISYMA,ISYMI)
     *                     + NVIR(ISYMA)*(I-1)
     *                     + A
C
                     XABI(KOFF1) = T2TP(KOFF2)
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      CALL QEXIT('SORT_T2_ABI')
C
      RETURN
      END
C  /* Deck sort_t2_baj */
      SUBROUTINE  SORT_T2_BAJ(XBAJ,ISYMI,I,T2TP,ISYMT2) 
C
C     t2tp(aijb) as I^I(baj)
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYMI, ISYMT2
      INTEGER ISYMBAJ, ISYMJ, ISYMBA, ISYMB, ISYMA, ISYMAIJ, ISYMAI
      INTEGER KOFF1, KOFF2
#if defined (SYS_CRAY)
      REAL XBAJ(*), T2TP(*) 
#else
      DOUBLE PRECISION XBAJ(*), T2TP(*) 
#endif
C
      CALL QENTER('SORT_T2_BAJ')
C
      ISYMBAJ = MULD2H(ISYMT2,ISYMI)
      DO ISYMJ = 1,NSYM
         ISYMBA = MULD2H(ISYMBAJ,ISYMJ)
         DO ISYMA = 1,NSYM 
            ISYMB = MULD2H(ISYMBA,ISYMA)
            ISYMAIJ = MULD2H(ISYMT2,ISYMB)
            ISYMAI  = MULD2H(ISYMAIJ,ISYMJ)
            DO J = 1,NRHF(ISYMJ)
               DO A = 1,NVIR(ISYMA)
                  DO B = 1,NVIR(ISYMB)
                     KOFF1 = IMAABI(ISYMBA,ISYMJ)
     *                     + NMATAB(ISYMBA)*(J-1)
     *                     + IMATAB(ISYMB,ISYMA)
     *                     + NVIR(ISYMB)*(A-1)
     *                     + B
C
                     KOFF2 = IT2SP(ISYMAIJ,ISYMB)
     *                     + NCKI(ISYMAIJ)*(B-1)
     *                     + ICKI(ISYMAI,ISYMJ)
     *                     + NT1AM(ISYMAI)*(J-1)
     *                     + IT1AM(ISYMA,ISYMI)
     *                     + NVIR(ISYMA)*(I-1)
     *                     + A
C
                     XBAJ(KOFF1) = T2TP(KOFF2)
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      CALL QEXIT('SORT_T2_BAJ')
C
      RETURN
      END
C
C  /* Deck sort_t2_ab */
      SUBROUTINE SORT_T2_AB(XAB,ISYMI,I,ISYMJ,J,T2TP,ISYMT2)
C
C   XAB^IJ(ab) = t2tp(aijb)
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYMI, ISYMJ, ISYMT2
      INTEGER ISYMIJ, ISYMAB, ISYMA, ISYMB, ISYMAIJ, ISYMAI
      INTEGER KOFF1, KOFF2
#if defined (SYS_CRAY)
      REAL XAB(*), T2TP(*)
#else
      DOUBLE PRECISION XAB(*), T2TP(*)
#endif
C
      CALL QENTER('SORT_T2_AB')
C
      ISYMIJ = MULD2H(ISYMI,ISYMJ)
      ISYMAB = MULD2H(ISYMT2,ISYMIJ)
      DO ISYMA = 1,NSYM
         ISYMB = MULD2H(ISYMAB,ISYMA)
         ISYMAIJ = MULD2H(ISYMA,ISYMIJ)
         ISYMAI  = MULD2H(ISYMA,ISYMI)
         DO A = 1,NVIR(ISYMA)
            DO B = 1,NVIR(ISYMB)
               KOFF1 =  IMATAB(ISYMA,ISYMB)
     *                + NVIR(ISYMA)*(B-1)
     *                + A
C
               KOFF2 = IT2SP(ISYMAIJ,ISYMB)
     *                + NCKI(ISYMAIJ)*(B-1)
     *                + ICKI(ISYMAI,ISYMJ)
     *                + NT1AM(ISYMAI)*(J-1)
     *                + IT1AM(ISYMA,ISYMI)
     *                + NVIR(ISYMA)*(I-1)
     *                + A
C
                     XAB(KOFF1) = T2TP(KOFF2)
            END DO
         END DO
      END DO
C
      CALL QEXIT('SORT_T2_AB')
C
      RETURN
      END
C  /* Deck t3jk_dia */
      SUBROUTINE T3JK_DIA(TBCAI,ISYMBCAI,ISYMJ,J,ISYMK,K,FOCKD)
C
C  T3^abc_iJK = T3^JK(bcai) = TBCAI/(F(A)+F(B)+F(C)-F(I)-F(J)-F(K))
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
C
      INTEGER ISYMBCAI,ISYMJ,ISYMK
      INTEGER NJ, NK, ISYMI, ISYMBCA, ISYMC, ISYMBC, ISYMB,ISYMA
      INTEGER NI, NC, NB, NA, KOFF1 
C
#if defined (SYS_CRAY)
      REAL TBCAI(*), FOCKD(*) 
      REAL DEN 
#else
      DOUBLE PRECISION TBCAI(*), FOCKD(*) 
      DOUBLE PRECISION DEN 
#endif
C
      CALL QENTER('T3JK_DIA')
C

      NJ = IORB(ISYMJ) + J
      NK = IORB(ISYMK) + K
      DO  ISYMI = 1,NSYM
        ISYMBCA = MULD2H(ISYMBCAI,ISYMI)
         DO ISYMA =  1,NSYM
            ISYMBC = MULD2H(ISYMBCA,ISYMA)
            DO ISYMC = 1,NSYM
               ISYMB =  MULD2H(ISYMBC,ISYMC)
               DO I = 1,NRHF(ISYMI)
                  NI = IORB(ISYMI) + I
                  DO C = 1,NVIR(ISYMC)
                     NC = IORB(ISYMC) + NRHF(ISYMC) + C
                     DO B = 1,NVIR(ISYMB)
                        NB = IORB(ISYMB) + NRHF(ISYMB) + B
                        DO A = 1,NVIR(ISYMA)
                           NA = IORB(ISYMA) + NRHF(ISYMA) + A
                           KOFF1 = IMAABCI(ISYMBCA,ISYMI)
     *                           + NMAABC(ISYMBCA)*(I-1)
     *                           + IMAABC(ISYMBC,ISYMA)
     *                           + NMATAB(ISYMBC)*(A-1)
     *                           + IMATAB(ISYMB,ISYMC)
     *                           + NVIR(ISYMB)*(C-1)
     *                           + B
                            DEN =  FOCKD(NA) + FOCKD(NB) + FOCKD(NC)
     *                          -  FOCKD(NI) - FOCKD(NJ) - FOCKD(NK)
C
                           TBCAI(KOFF1) = TBCAI(KOFF1)/DEN
C
                        END DO   
                     END DO   
                  END DO   
               END DO   
            END DO   
         END DO   
      END DO   
C
      CALL QEXIT('T3JK_DIA')
C
      RETURN
      END
C  /* Deck sum_pt3_jk */
      SUBROUTINE SUM_PT3_JK(TMAT,ISYMJ,J,ISYMK,K,ISYMBCAI,T3SUM,IOPT)
C
C     Sum up the T3 amplitudes from the TMAT.
C
C     IOPT = 1: Start from (ai,bj,ck) 
C     IOPT = 2: Start from (bj,ai,ck)
C     IOPT = 3: Start from (bj,ck,ai)
C     IOPT = 4: w (=divide by 2)
C     IOPT = 5: Start from (ck,ai,bj)
C     IOPT = 6: t (=divide by 6)
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYMJ,ISYMK,ISYMBCAI,IOPT
      INTEGER ISYMBCA,ISYMBC
      INTEGER ISYMI,ISYMA,ISYMC,ISYMB
      INTEGER KOFF1,KOFF2,KOFF3,KOFF4,KOFF5,KOFF6,KOFF7,KH
C
#if defined (SYS_CRAY)
      REAL TMAT(*)
      REAL ONE,SIXTH,HALF,FACT
      REAL t3sum(nvirt,nvirt,nvirt,nrhft,nrhft,nrhft)
#else
      DOUBLE PRECISION TMAT(*)
      DOUBLE PRECISION ONE,SIXTH,HALF,FACT
      DOUBLE PRECISION t3sum(nvirt,nvirt,nvirt,nrhft,nrhft,nrhft)
#endif
C
      LOGICAL LDEBUG
C
      PARAMETER (ONE = 1.0D0, SIXTH = 1.0D0/6.0D0, HALF = 1.0D0/2.0D0)
C
      CALL QENTER('SUM_PT3_JK')
C
      LDEBUG = .true.
C
C
      IF ( (IOPT .GE. 1) .AND. (IOPT .LE. 3) ) THEN
         FACT = ONE
      ELSE IF ((IOPT .EQ. 4) ) THEN
         FACT = HALF
      ELSE IF ((IOPT .EQ. 6) ) THEN
         FACT = SIXTH
      ELSE
         CALL QUIT('Wrong IOPT in SUM_PT3_JK')
      ENDIF
C
      KOFF5 = 0
      DO KH = 1, ISYMJ-1
        KOFF5 = KOFF5 + NRHF(KH)
      ENDDO
      KOFF6 = 0
      DO KH = 1, ISYMK-1
        KOFF6 = KOFF6 + NRHF(KH)
      ENDDO
C
      DO ISYMI = 1, NSYM
         KOFF4 = 0
         DO KH = 1, ISYMI-1
           KOFF4 = KOFF4 + NRHF(KH)
         ENDDO
         ISYMBCA = MULD2H(ISYMI,ISYMBCAI)
         DO ISYMA = 1, NSYM
            KOFF1 = 0
            DO KH = 1, ISYMA-1
              KOFF1 = KOFF1 + NVIR(KH)
            ENDDO
            ISYMBC = MULD2H(ISYMBCA,ISYMA)
            DO ISYMC = 1, NSYM
C
               KOFF3 = 0
               DO KH = 1, ISYMC-1
                  KOFF3 = KOFF3 + NVIR(KH)
               ENDDO
               ISYMB = MULD2H(ISYMBC,ISYMC)
C
               KOFF2 = 0
               DO KH = 1, ISYMB-1
                  KOFF2 = KOFF2 + NVIR(KH)
               ENDDO
C
               DO I = 1, NRHF(ISYMI)
               DO A = 1, NVIR(ISYMA)
               DO C = 1, NVIR(ISYMC)
               DO B = 1, NVIR(ISYMB)
C
               IF ( (IOPT .GE. 1) .AND. (IOPT .LE. 6) ) THEN
                 KOFF7 = IMAABCI(ISYMBCA,ISYMI)
     *                 + NMAABC(ISYMBCA)*(I - 1)
     *                 + IMAABC(ISYMBC,ISYMA) + NMATAB(ISYMBC)*(A-1)
     *                 + IMATAB(ISYMB,ISYMC) + NVIR(ISYMB)*(C-1) + B
               ENDIF
C
         IF (LDEBUG) THEN
            IF (ABS(TMAT(KOFF7)) .GT. 1.0d-12) THEN
              IF ( (IOPT .GE. 1) .AND. (IOPT .LE. 6) ) THEN
                 WRITE(LUPRI,*) 'CONTRIBUTION FROM I = ',
     *                              KOFF7,' WITH TMAT = ',
     *                              TMAT(KOFF7)
              ENDIF
C
               WRITE(LUPRI,'(A,6I3)') ' A, B, C, I, J, K : ',A,B,C,I,J,K
               WRITE(LUPRI,*) 'KOFF1 = ',KOFF1,' KOFF2 = ',
     *                         KOFF2,'KOFF3 = ',KOFF3
               WRITE(LUPRI,*) 'KOFF4 = ',KOFF4,' KOFF5 = ',
     *                         KOFF5,'KOFF6 = ',KOFF6
            ENDIF
         ENDIF
C
         T3SUM(A+KOFF1,B+KOFF2,C+KOFF3,I+KOFF4,J+KOFF5,K+KOFF6) = 
     *   T3SUM(A+KOFF1,B+KOFF2,C+KOFF3,I+KOFF4,J+KOFF5,K+KOFF6) +
     *   FACT*TMAT(KOFF7)
C
         T3SUM(A+KOFF1,C+KOFF3,B+KOFF2,I+KOFF4,K+KOFF6,J+KOFF5) = 
     *   T3SUM(A+KOFF1,C+KOFF3,B+KOFF2,I+KOFF4,K+KOFF6,J+KOFF5) +
     *   FACT*TMAT(KOFF7)
C
         T3SUM(B+KOFF2,A+KOFF1,C+KOFF3,J+KOFF5,I+KOFF4,K+KOFF6) = 
     *   T3SUM(B+KOFF2,A+KOFF1,C+KOFF3,J+KOFF5,I+KOFF4,K+KOFF6) +
     *   FACT*TMAT(KOFF7)
C
         T3SUM(B+KOFF2,C+KOFF3,A+KOFF1,J+KOFF5,K+KOFF6,I+KOFF4) = 
     *   T3SUM(B+KOFF2,C+KOFF3,A+KOFF1,J+KOFF5,K+KOFF6,I+KOFF4) +
     *   FACT*TMAT(KOFF7)
C
         T3SUM(C+KOFF3,A+KOFF1,B+KOFF2,K+KOFF6,I+KOFF4,J+KOFF5) = 
     *   T3SUM(C+KOFF3,A+KOFF1,B+KOFF2,K+KOFF6,I+KOFF4,J+KOFF5) +
     *   FACT*TMAT(KOFF7)
C
         T3SUM(C+KOFF3,B+KOFF2,A+KOFF1,K+KOFF6,J+KOFF5,I+KOFF4) = 
     *   T3SUM(C+KOFF3,B+KOFF2,A+KOFF1,K+KOFF6,J+KOFF5,I+KOFF4) +
     *   FACT*TMAT(KOFF7)
C
                 ENDDO
                 ENDDO
                 ENDDO
                 ENDDO
            ENDDO
C
         ENDDO
C
      ENDDO
C
      CALL QEXIT('SUM_PT3_JK')
C
    1 FORMAT(1X,A8,I3,A1,I3,A1,I3,A1,I3,A1,I3,A1,I3,A4,E20.10)
      RETURN
      END
C  /* Deck w3jk_dia */
      SUBROUTINE W3JK_DIA(TBCAI,ISYMBCAI,ISYMJ,J,ISYMK,K,FOCKD,FREQ)
C
C  W3^abc_iJK = W3^JK(bcai) = WBCAI/(F(A)+F(B)+F(C)-F(I)-F(J)-F(K) - FREQ)
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
C
      INTEGER ISYMBCAI,ISYMJ,ISYMK
      INTEGER NJ, NK, ISYMI, ISYMBCA, ISYMC, ISYMBC, ISYMB,ISYMA
      INTEGER NI, NC, NB, NA, KOFF1 
C
#if defined (SYS_CRAY)
      REAL TBCAI(*), FOCKD(*) 
      REAL DEN,FREQ 
#else
      DOUBLE PRECISION TBCAI(*), FOCKD(*) 
      DOUBLE PRECISION DEN,FREQ 
#endif
C
      CALL QENTER('W3JK_DIA')
C

      NJ = IORB(ISYMJ) + J
      NK = IORB(ISYMK) + K
      DO  ISYMI = 1,NSYM
        ISYMBCA = MULD2H(ISYMBCAI,ISYMI)
         DO ISYMA =  1,NSYM
            ISYMBC = MULD2H(ISYMBCA,ISYMA)
            DO ISYMC = 1,NSYM
               ISYMB =  MULD2H(ISYMBC,ISYMC)
               DO I = 1,NRHF(ISYMI)
                  NI = IORB(ISYMI) + I
                  DO C = 1,NVIR(ISYMC)
                     NC = IORB(ISYMC) + NRHF(ISYMC) + C
                     DO B = 1,NVIR(ISYMB)
                        NB = IORB(ISYMB) + NRHF(ISYMB) + B
                        DO A = 1,NVIR(ISYMA)
                           NA = IORB(ISYMA) + NRHF(ISYMA) + A
                           KOFF1 = IMAABCI(ISYMBCA,ISYMI)
     *                           + NMAABC(ISYMBCA)*(I-1)
     *                           + IMAABC(ISYMBC,ISYMA)
     *                           + NMATAB(ISYMBC)*(A-1)
     *                           + IMATAB(ISYMB,ISYMC)
     *                           + NVIR(ISYMB)*(C-1)
     *                           + B
                            DEN =  FOCKD(NA) + FOCKD(NB) + FOCKD(NC)
     *                          -  FOCKD(NI) - FOCKD(NJ) - FOCKD(NK)
     *                          -  FREQ
C
                           TBCAI(KOFF1) = TBCAI(KOFF1)/DEN
C
                        END DO   
                     END DO   
                  END DO   
               END DO   
            END DO   
         END DO   
      END DO   
C
      CALL QEXIT('W3JK_DIA')
C
      RETURN
      END
C  /* Deck t3_forbidden_jk */
      SUBROUTINE T3_FORBIDDEN_JK(TMAT,ISYMIM,ISYMJ,J,ISYMK,K)
C
C     Purpose : Remove the forbidden t3/t3-bar amplitudes.
C
      IMPLICIT NONE
C
#include <ccsdsym.h>
#include <ccorb.h>
#include <priunit.h>
C
      INTEGER ISYMIM,ISYMJ,ISYMK,ISYABC,ISYMC,ISYMAB,ISYMB
      INTEGER ISYMA,KOFF1,ISYMAA,ISYMJK,ISYMAI,ISYMI
C
#if defined (SYS_CRAY)
      REAL TMAT(*),ZERO
#else
      DOUBLE PRECISION TMAT(*),ZERO
#endif
C
      PARAMETER (ZERO = 0.0D0)
C
      CALL QENTER('T3_FORBIDDEN_JK')
C
C---------------------------------------------------------
C     If J and K are the same remove all amplitudes
C     having an I which is the same as J and K.
C---------------------------------------------------------
C
      IF ((ISYMJ .EQ. ISYMK) .AND. (J .EQ. K) ) THEN
         ISYABC = MULD2H(ISYMJ,ISYMIM)
         DO ISYMC = 1, NSYM
            ISYMAB = MULD2H(ISYMC,ISYABC)
            DO ISYMB = 1, NSYM
               ISYMA  = MULD2H(ISYMB,ISYMAB)
C
               DO C = 1, NVIR(ISYMC)
               DO B = 1, NVIR(ISYMB)
               DO A = 1, NVIR(ISYMA)
                  KOFF1 = IMAABCI(ISYABC,ISYMJ)
     *                  + NMAABC(ISYABC)*(J-1)
     *                  + IMAABC(ISYMAB,ISYMC)
     *                  + NMATAB(ISYMAB)*(C-1)
     *                  + IMATAB(ISYMA,ISYMB)
     *                  + NVIR(ISYMA)*(B-1)
     *                  + A
C
                  TMAT(KOFF1)  = ZERO
C
               ENDDO
               ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF
C
C----------------------------------------------------------
C     Remove all amplitudes that has three indentical
C     virtual indices.
C----------------------------------------------------------
C
      ISYMAA = 1
      ISYMJK = MULD2H(ISYMJ,ISYMK)
      ISYMAI = MULD2H(ISYMIM,ISYMJK)
      DO ISYMA = 1, NSYM
         ISYMI = MULD2H(ISYMA,ISYMAI)
         DO A = 1, NVIR(ISYMA)
            DO I = 1, NRHF(ISYMI)
               KOFF1 = IMAABCI(ISYMA,ISYMI)
     *                 + NMAABC(ISYMA)*(I-1)
     *                 + IMAABC(ISYMAA,ISYMA)
     *                 + NMATAB(ISYMAA)*(A-1)
     *                 + IMATAB(ISYMA,ISYMA)
     *                 + NVIR(ISYMA)*(A-1)
     *                 + A

C
               TMAT(KOFF1)  = ZERO
C
            ENDDO
         ENDDO
      ENDDO
C
C-----------------------
C     End.
C-----------------------
C
      CALL QEXIT('T3_FORBIDDEN_JK')
C
      RETURN
      END
c  /* deck print_t3jk */
c
      subroutine print_t3jk(tmat,isymim,j,k,iopt)
c
c     Print t3^{JK}_{bcai} amplitudes which have been 
c     summed up (outside) to tmat array. (IOPT = 1)
c
c     Print W^{JK}_{bcai} "amplitudes" which have been 
c     summed up (outside) to wmat array. (IOPT = 2)
c
c     isymim is symmetry of (bcai)
c
c     P. Jorgensen and F. Pawlowski, Spring 2003.
c
      IMPLICIT NONE
c
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYMIM, KOFF1, KOFF4, KOFF5, KOFF6, KOFF7
      INTEGER ISYMA, ISYMI, ISYMB, ISYMC
      INTEGER KH, ISYCAI, ISYMBC, ISYMAI, ISYBCA
      INTEGER IOPT
C
#if defined (SYS_CRAY)
      REAL tmat(*)
#else
      DOUBLE PRECISION tmat(*)
#endif
C
      CALL QENTER('print_t3jk')
C
C----------------------------------------
C     Print the triples amplitudes.
C----------------------------------------
C
      DO ISYMB = 1, NSYM
C
         KOFF1 = 0
         DO KH = 1, ISYMB-1
           KOFF1 = KOFF1 + NVIR(KH)
         ENDDO
C
         ISYCAI = MULD2H(ISYMIM,ISYMB)
C
         DO ISYMC = 1, NSYM
C
            KOFF4 = 0
            DO KH = 1, ISYMC-1
              KOFF4 = KOFF4 + NVIR(KH)
            ENDDO
C
            ISYMAI = MULD2H(ISYCAI,ISYMC)
            ISYMBC = MULD2H(ISYMB,ISYMC)
C
            DO ISYMA = 1, NSYM
C
               KOFF5 = 0
               DO KH = 1, ISYMA-1
                 KOFF5 = KOFF5 + NVIR(KH)
               ENDDO
C
               ISYMI   = MULD2H(ISYMAI,ISYMA)
               ISYBCA = MULD2H(ISYMBC,ISYMA)
C
               KOFF6 = 0
               DO KH = 1, ISYMI-1
                 KOFF6 = KOFF6 + NRHF(KH)
               ENDDO
C
                DO B = 1, NVIR(ISYMB)
                DO C = 1, NVIR(ISYMC)
                DO A = 1, NVIR(ISYMA)
                DO I = 1, NRHF(ISYMI)
C
                KOFF7 = IMAABCI(ISYBCA,ISYMI)
     *                + NMAABC(ISYBCA)*(I - 1)
     *                + IMAABC(ISYMBC,ISYMA) + NMATAB(ISYMBC)*(A-1)
     *                + IMATAB(ISYMB,ISYMC) + NVIR(ISYMB)*(C-1) + B

C
        IF (ABS(TMAT(KOFF7)) 
     *                                    .GT. 1.0D-12) THEN
           IF (IOPT .EQ. 1) THEN
              write(lupri,1) 'T3JK(',a+koff5,',',b+koff1,',',
     *                               c+koff4,',',i+koff6,',',
     *                               j,',',k,') = ',
     *        TMAT(KOFF7) 
           ELSE IF (IOPT .EQ. 2) THEN
              write(lupri,1) 'WJK(',a+koff5,',',b+koff1,',',
     *                              c+koff4,',',i+koff6,',',
     *                              j,',',k,') = ',
     *        TMAT(KOFF7) 
           ELSE 
              call quit('Illegal value for IOPT in print_t3jk ')
           ENDIF
C
        ENDIF
C
                ENDDO 
                ENDDO
                ENDDO 
                ENDDO 
                ENDDO 
                ENDDO 
C
      ENDDO          
C
      CALL QEXIT('print_t3jk')
C
    1 FORMAT(1X,A6,I3,A1,I3,A1,I3,A1,I3,A1,I3,A1,I3,A4,E20.10)
      RETURN
      END
C  /* Deck get_t3b0_jk_o */
      SUBROUTINE GET_T3B0_JK_O(T3B0JK,ISYT3B0JK,
     *                           L2TP,ISYML2,
     *                           T3BOL2,T3BOG2,ISYINT,
     *                           ISYMJ,J,ISYMK,K,
     *                           WORK,LWORK)
********************************************************************
*
* T3BOL2 : L(ia|jk) sorted as T3BOL2(ajik)
* T3BOG2 : g(ia|jk) sorted as T3BOG2(ajik)
*
********************************************************************
*
* OBS !
*       t in the following comments of this routine denotes 
*         Lagrange multipliers 
*
********************************************************************
*
* In this routine we calculate these contributions to t3bar_0
* multipliers, which contain the integrals of the type (ooov).
*
* We thus calculate the following intermmediate (for two fixed 
* occupied index):
*
* T^JK(bcai) 
*       = P(ai,bj,ck) ( - t^(ab)_(in) L(jn|kc) + t^(ab)_(nk) g(in|jc) )
*
* (1)
*       = - t^(ab)_(in) L(jn|kc) + t^(ab)_(nk) g(in|jc)
* (2)
*         - t^(ba)_(jn) L(in|kc) + t^(ba)_(nk) g(jn|ic)
* (3)
*         - t^(bc)_(jn) L(kn|ia) + t^(bc)_(ni) g(jn|ka)
* (4)
*         - t^(ac)_(in) L(kn|jb) + t^(ac)_(nj) g(in|kb)
* (5)
*         - t^(ca)_(kn) L(in|jb) + t^(ca)_(nj) g(kn|ib)
* (6)
*         - t^(cb)_(kn) L(jn|ia) + t^(cb)_(ni) g(kn|ja)
*
* Filip Pawlowski, Aarhus, Winter 2003
*
********************************************************************
      IMPLICIT NONE
#include<ccsdsym.h>
#include<ccorb.h>
#include<priunit.h>
C
      INTEGER ISYT3B0JK,ISYML2,ISYINT,ISYMJ,ISYMK,LWORK
      INTEGER ISYKJ,ISYCN,ISYMN,ISYMC,ISYBAI
      INTEGER NTOTC,NTOTN
      INTEGER KL2NBAI,KINTCN,KCBAI,KEND1,LWRK1
      INTEGER ISYL2BAN,ISYINTNCI,ISYBA,ISYCI
      INTEGER NTOTBA
      INTEGER KL2BAN,KINTNCI,KBACI
      INTEGER ISYJK,ISYBN,ISYMB,ISYCAI
      INTEGER NTOTB
      INTEGER ISYL2CAN,ISYINTNBI,ISYCA,ISYBI
      INTEGER NTOTCA
      INTEGER KL2CAN,KINTNBI,KCABI
      INTEGER KL2NCAI,KINTBN
      INTEGER ISYL2BCN,ISYINTNAI,ISYBC,ISYAI
      INTEGER NTOTBC
      INTEGER KL2BCN,KINTNAI
      INTEGER KOFF1,KOFF2,KOFF3
      INTEGER ISYAN,ISYMA,ISYCBI,NTOTA
      INTEGER KL2NCBI,KINTAN,KACBI
      INTEGER ISYL2CBN,ISYCB,ISYBCI
      INTEGER NTOTCB
      INTEGER KL2CBN,KL2NBCI,KABCI
C
#if defined (SYS_CRAY)
      REAL T3B0JK(*),L2TP(*),T3BOL2(*),T3BOG2(*),WORK(LWORK)
      REAL ONE
#else
      DOUBLE PRECISION T3B0JK(*),L2TP(*),T3BOL2(*),T3BOG2(*),WORK(LWORK)
      DOUBLE PRECISION ONE
#endif
C
      PARAMETER (ONE = 1.0D0)
C
      CALL QENTER('GET_T3B0_JK_O')
C
C=================================================
C     Calculate (1a)   - t^(ab)_(in) L(jn|kc)
C                    = - t^(ba)_(ni) L(kc|jn)
C                             
C                      - T(nbai) I^KJ(cn)
C=================================================
C
C-------------------------------
C     Sort L2TP(bnia) as T(nbai)
C-------------------------------
C
      ISYKJ = MULD2H(ISYMK,ISYMJ)
      ISYCN = MULD2H(ISYINT,ISYKJ)
C
      KL2NBAI = 1
      KINTCN  = KL2NBAI + NT2SQ(ISYML2)
      KCBAI   = KINTCN  + NT1AM(ISYCN)
      KEND1   = KCBAI   + NMAAOBCI(ISYT3B0JK)
      LWRK1   = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3B0_JK_O (1a)')
      END IF
C
      CALL DZERO(WORK(KCBAI),NMAAOBCI(ISYT3B0JK))
C
      CALL SORT_T2_IABJ(WORK(KL2NBAI),L2TP,ISYML2)
C
C------------------------------------------------
C     Sort L(kc|jn) = T3BOL2(c,j,k,n) as I^KJ(cn)
C------------------------------------------------
C
      CALL SORT_INT_AJ_IK(WORK(KINTCN),T3BOL2,ISYINT,ISYMK,K,ISYMJ,J)
C
C------------------------------------------
C    Multiply I^KJ(cn) T(nbai) = T^JK(cbai)
C------------------------------------------
C
      DO ISYMN = 1, NSYM
            ISYMC = MULD2H(ISYCN,ISYMN)
            ISYBAI = MULD2H(ISYML2,ISYMN)
C      
            KOFF1 = KINTCN
     *            + IT1AM(ISYMC,ISYMN)
            KOFF2 = KL2NBAI
     *            + IMAJBAI(ISYMN,ISYBAI)
            KOFF3 = KCBAI
     *            + IMAAOBCI(ISYMC,ISYBAI)
C
            NTOTC = MAX(NVIR(ISYMC),1)
            NTOTN = MAX(NRHF(ISYMN),1)
C
            CALL DGEMM('N','N',NVIR(ISYMC),NMAABI(ISYBAI),NRHF(ISYMN),
     *                 -ONE,WORK(KOFF1),NTOTC,WORK(KOFF2),NTOTN,
     *                 ONE,WORK(KOFF3),NTOTC)
C
      END DO ! ISYMN
C
C     T3B0JK(bcai) = T3B0JK(bcai) + T^JK(cbai)
C
C  bar_occ(1a)
C
      CALL FBACI(T3B0JK,WORK(KCBAI),ISYT3B0JK)
C
C=================================================
C     Calculate (1b)   t^(ab)_(nk) g(in|jc)
C                    = t^(ba)_(kn) g(jc|in)
C                    
C                      T^K(ban) I^J(nci)
C=================================================
C
C-------------------------------
C     Sort L2TP(bkna) as T^K(ban)
C-------------------------------
C
      ISYL2BAN = MULD2H(ISYML2,ISYMK)
      ISYINTNCI = MULD2H(ISYINT,ISYMJ)
C
      KL2BAN = KEND1
      KINTNCI  = KL2BAN + NMAABI(ISYL2BAN)
      KBACI   = KINTNCI + NCKI(ISYINTNCI)
      KEND1 = KBACI + NMAAB_CI(ISYT3B0JK)
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3BAR0_IJ_O (1b)')
      END IF
C
      CALL DZERO(WORK(KBACI),NMAAB_CI(ISYT3B0JK))
C
      CALL SORT_T2_ABJ(WORK(KL2BAN),ISYMK,K,L2TP,ISYML2)
C
C----------------------------------------------
C     Sort g(jc|in) = T3OG2(c,i,j,n) as I^J(nci)
C----------------------------------------------
C
      CALL SORT_INT_JAK_I(WORK(KINTNCI),T3BOG2,ISYINT,ISYMJ,J)
C
C------------------------------------------
C    Multiply T^K(ban) * I^J(nci) = T^JK(baci)
C------------------------------------------
C
      DO ISYMN = 1, NSYM
         ISYBA = MULD2H(ISYL2BAN,ISYMN)
         ISYCI = MULD2H(ISYINTNCI,ISYMN)
C
         KOFF1 = KL2BAN
     *         + IMAABI(ISYBA,ISYMN)
         KOFF2 = KINTNCI
     *         + IMAIAJ(ISYMN,ISYCI)
         KOFF3 = KBACI
     *         + IMAAB_CI(ISYBA,ISYCI)
C
         NTOTBA = MAX(NMATAB(ISYBA),1)
         NTOTN  = MAX(NRHF(ISYMN),1)
C
         CALL DGEMM('N','N',NMATAB(ISYBA),NT1AM(ISYCI),NRHF(ISYMN),
     *              ONE,WORK(KOFF1),NTOTBA,WORK(KOFF2),NTOTN,
     *              ONE,WORK(KOFF3),NTOTBA)
C
      END DO ! ISYMN
C
C-------------------------------------------
C     T3B0JK(bcai) = T3B0JK(bcai) + T^JK(baci)
C-------------------------------------------
C
C  bar_occ(1b)
C
      CALL FACBI(T3B0JK,WORK(KBACI),ISYT3B0JK)
C
C=================================================
C     Calculate (2a)   - t^(ba)_(jn) L(in|kc)
C                    = - t^(ba)_(jn) L(kc|in)
C                    
C                      - T^J(ban) I^K(nci)
C=================================================
C
C-------------------------------
C     Sort L2TP(bjna) as T^J(ban)
C-------------------------------
C
      ISYL2BAN = MULD2H(ISYML2,ISYMJ)
      ISYINTNCI = MULD2H(ISYINT,ISYMK)
C
      KL2BAN = KEND1
      KINTNCI  = KL2BAN + NMAABI(ISYL2BAN)
      KBACI   = KINTNCI + NCKI(ISYINTNCI)
      KEND1 = KBACI + NMAAB_CI(ISYT3B0JK)
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3BAR0_IJ_O (2a)')
      END IF
C
      CALL DZERO(WORK(KBACI),NMAAB_CI(ISYT3B0JK))
C
      CALL SORT_T2_ABJ(WORK(KL2BAN),ISYMJ,J,L2TP,ISYML2)
C
C----------------------------------------------
C     Sort L(kc|in) = T3OL2(c,i,k,n) as I^K(nci)
C----------------------------------------------
C
      CALL SORT_INT_JAK_I(WORK(KINTNCI),T3BOL2,ISYINT,ISYMK,K)
C
C------------------------------------------
C    Multiply T^J(ban) * I^K(nci) = T^JK(baci)
C------------------------------------------
C
      DO ISYMN = 1, NSYM
         ISYBA = MULD2H(ISYL2BAN,ISYMN)
         ISYCI = MULD2H(ISYINTNCI,ISYMN)
C
         KOFF1 = KL2BAN
     *         + IMAABI(ISYBA,ISYMN)
         KOFF2 = KINTNCI
     *         + IMAIAJ(ISYMN,ISYCI)
         KOFF3 = KBACI
     *         + IMAAB_CI(ISYBA,ISYCI)
C
         NTOTBA = MAX(NMATAB(ISYBA),1)
         NTOTN  = MAX(NRHF(ISYMN),1)
C
         CALL DGEMM('N','N',NMATAB(ISYBA),NT1AM(ISYCI),NRHF(ISYMN),
     *              -ONE,WORK(KOFF1),NTOTBA,WORK(KOFF2),NTOTN,
     *              ONE,WORK(KOFF3),NTOTBA)
C
      END DO ! ISYMN
C
C-------------------------------------------
C     T3B0JK(bcai) = T3B0JK(bcai) + T^JK(baci)
C-------------------------------------------
C
C  bar_occ(2a)
C
      CALL FACBI(T3B0JK,WORK(KBACI),ISYT3B0JK)
C
C=================================================
C     Calculate (2b)   t^(ab)_(kn) g(jn|ic)
C                    = t^(ba)_(nk) g(ic|jn)
C                    
C                      T^K(ban) I^J(nci)
C=================================================
C
C-------------------------------
C     Sort L2TP(bnka) as T^K(ban)
C-------------------------------
C
      ISYL2BAN = MULD2H(ISYML2,ISYMK)
      ISYINTNCI = MULD2H(ISYINT,ISYMJ)
C
      KL2BAN = KEND1
      KINTNCI  = KL2BAN + NMAABI(ISYL2BAN)
      KBACI   = KINTNCI + NCKI(ISYINTNCI)
      KEND1 = KBACI + NMAAB_CI(ISYT3B0JK)
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3BAR0_IJ_O (2b)')
      END IF
C
      CALL DZERO(WORK(KBACI),NMAAB_CI(ISYT3B0JK))
C
      CALL SORT_T2_ABI(WORK(KL2BAN),ISYMK,K,L2TP,ISYML2)
C
C----------------------------------------------
C     Sort g(ic|jn) = T3OG2(c,j,i,n) as I^J(nci)
C----------------------------------------------
C
      CALL SORT_INT_JAI_K(WORK(KINTNCI),T3BOG2,ISYINT,ISYMJ,J)
C
C------------------------------------------
C    Multiply T^K(ban) * I^J(nci) = T^JK(baci)
C------------------------------------------
C
      DO ISYMN = 1, NSYM
         ISYBA = MULD2H(ISYL2BAN,ISYMN)
         ISYCI = MULD2H(ISYINTNCI,ISYMN)
C
         KOFF1 = KL2BAN
     *         + IMAABI(ISYBA,ISYMN)
         KOFF2 = KINTNCI
     *         + IMAIAJ(ISYMN,ISYCI)
         KOFF3 = KBACI
     *         + IMAAB_CI(ISYBA,ISYCI)
C
         NTOTBA = MAX(NMATAB(ISYBA),1)
         NTOTN  = MAX(NRHF(ISYMN),1)
C
         CALL DGEMM('N','N',NMATAB(ISYBA),NT1AM(ISYCI),NRHF(ISYMN),
     *              ONE,WORK(KOFF1),NTOTBA,WORK(KOFF2),NTOTN,
     *              ONE,WORK(KOFF3),NTOTBA)
C
      END DO ! ISYMN
C
C-------------------------------------------
C     T3B0JK(bcai) = T3B0JK(bcai) + T^JK(baci)
C-------------------------------------------
C
C  bar_occ(2b)
C
      CALL FACBI(T3B0JK,WORK(KBACI),ISYT3B0JK)
C
C=================================================
C     Calculate (3a)   - t^(bc)_(jn) L(kn|ia)
C                    = - t^(cb)_(nj) L(ia|kn)
C                    
C                      T^J(cbn) I^K(nai)
C=================================================
C
C-------------------------------
C     Sort L2TP(cnjb) as T^J(cbn)
C-------------------------------
C
      ISYL2CBN = MULD2H(ISYML2,ISYMJ)
      ISYINTNAI = MULD2H(ISYINT,ISYMK)
C
      KL2CBN = KEND1
      KINTNAI  = KL2CBN + NMAABI(ISYL2CBN)
      KCBAI   = KINTNAI + NCKI(ISYINTNAI)
      KEND1   = KCBAI   + NMAAB_CI(ISYT3B0JK)
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3BAR0_IJ_O (3a)')
      END IF
C
      CALL DZERO(WORK(KCBAI),NMAAB_CI(ISYT3B0JK))
C
      CALL SORT_T2_ABI(WORK(KL2CBN),ISYMJ,J,L2TP,ISYML2)
C
C----------------------------------------------
C     Sort L(ia|kn) = T3OL2(a,k,i,n) as I^K(nai)
C----------------------------------------------
C
      CALL SORT_INT_JAI_K(WORK(KINTNAI),T3BOL2,ISYINT,ISYMK,K)
C
C------------------------------------------
C    Multiply T^J(cbn) * I^K(nai) = T^JK(cbai)
C------------------------------------------
C
      DO ISYMN = 1, NSYM
         ISYCB = MULD2H(ISYL2CBN,ISYMN)
         ISYAI = MULD2H(ISYINTNAI,ISYMN)
C
         KOFF1 = KL2CBN
     *         + IMAABI(ISYCB,ISYMN)
         KOFF2 = KINTNAI
     *         + IMAIAJ(ISYMN,ISYAI)
         KOFF3 = KCBAI
     *         + IMAAB_CI(ISYCB,ISYAI) 
C
         NTOTCB = MAX(NMATAB(ISYCB),1)
         NTOTN  = MAX(NRHF(ISYMN),1)
C
         CALL DGEMM('N','N',NMATAB(ISYCB),NT1AM(ISYAI),NRHF(ISYMN),
     *              -ONE,WORK(KOFF1),NTOTCB,WORK(KOFF2),NTOTN,
     *              ONE,WORK(KOFF3),NTOTCB)
C
C
C     T30JK(bcai) = T30JK(bcai) + T^JK(cbai)
C
C  bar_occ(3a)
C
      CALL FBACI(T3B0JK,WORK(KCBAI),ISYT3B0JK)

      END DO ! ISYMN
C
C=================================================
C     Calculate (3b)   t^(cb)_(in) g(jn|ka)
C                    = t^(bc)_(ni) g(ka|jn)
C                             
C                      T(nbci) I^KJ(an)
C=================================================
C
C-------------------------------
C     Sort L2TP(bnic) as T(nbci)
C-------------------------------
C
      ISYKJ = MULD2H(ISYMK,ISYMJ)
      ISYAN = MULD2H(ISYINT,ISYKJ)
C
      KL2NBCI = KEND1
      KINTAN  = KL2NBCI + NT2SQ(ISYML2)
      KABCI   = KINTAN  + NT1AM(ISYAN)
      KEND1   = KABCI   + NMAAOBCI(ISYT3B0JK)
      LWRK1   = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3B0_JK_O (3b)')
      END IF
C
      CALL DZERO(WORK(KABCI),NMAAOBCI(ISYT3B0JK))
C
      CALL SORT_T2_IABJ(WORK(KL2NBCI),L2TP,ISYML2)
C
C------------------------------------------------
C     Sort g(ka|jn) = T3BOG2(a,j,k,n) as I^KJ(an)
C------------------------------------------------
C
      CALL SORT_INT_AJ_IK(WORK(KINTAN),T3BOG2,ISYINT,ISYMK,K,ISYMJ,J)
C
C------------------------------------------
C    Multiply I^KJ(an) T(nbci) = T^JK(abci)
C------------------------------------------
C
      DO ISYMN = 1, NSYM
            ISYMA = MULD2H(ISYAN,ISYMN)
            ISYBCI = MULD2H(ISYML2,ISYMN)
C      
            KOFF1 = KINTAN
     *            + IT1AM(ISYMA,ISYMN)
            KOFF2 = KL2NBCI
     *            + IMAJBAI(ISYMN,ISYBCI)
            KOFF3 = KABCI
     *            + IMAAOBCI(ISYMA,ISYBCI) 
C
            NTOTA = MAX(NVIR(ISYMA),1)
            NTOTN = MAX(NRHF(ISYMN),1)
C
            CALL DGEMM('N','N',NVIR(ISYMA),NMAABI(ISYBCI),NRHF(ISYMN),
     *                 ONE,WORK(KOFF1),NTOTA,WORK(KOFF2),NTOTN,
     *                 ONE,WORK(KOFF3),NTOTA)
C
      END DO ! ISYMN
C
C     T30JK(bcai) = T30JK(bcai) + T^JK(abci)
C
C  bar_occ(3b)
C
       CALL FCABI(T3B0JK,WORK(KABCI),ISYT3B0JK)
C
C=================================================
C     Calculate (4a)   - t^(ac)_(in) L(kn|jb)
C                    = - t^(ca)_(ni) L(jb|kn)
C                             
C                      - T(ncai) I^JK(bn)
C=================================================
C
C-------------------------------
C     Sort L2TP(cnia) as T(ncai)
C-------------------------------
C
      ISYJK = MULD2H(ISYMJ,ISYMK)
      ISYBN = MULD2H(ISYINT,ISYJK)
C
      KL2NCAI = KEND1
      KINTBN  = KL2NCAI + NT2SQ(ISYML2)
      KEND1   = KINTBN  + NT1AM(ISYBN)
      LWRK1   = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3B0_JK_O (4a)')
      END IF
C
      CALL SORT_T2_IABJ(WORK(KL2NCAI),L2TP,ISYML2)
C
C------------------------------------------------
C     Sort L(jb|kn) = T3BOL2(b,k,j,n) as I^JK(bn)
C------------------------------------------------
C
      CALL SORT_INT_AJ_IK(WORK(KINTBN),T3BOL2,ISYINT,ISYMJ,J,ISYMK,K)
C
C------------------------------------------
C    Multiply I^JK(bn) T(ncai) = T^JK(bcai)
C------------------------------------------
C
      DO ISYMN = 1, NSYM
            ISYMB = MULD2H(ISYBN,ISYMN)
            ISYCAI = MULD2H(ISYML2,ISYMN)
C      
            KOFF1 = KINTBN
     *            + IT1AM(ISYMB,ISYMN)
            KOFF2 = KL2NCAI
     *            + IMAJBAI(ISYMN,ISYCAI)
            KOFF3 = IMAAOBCI(ISYMB,ISYCAI) + 1
C
            NTOTB = MAX(NVIR(ISYMB),1)
            NTOTN = MAX(NRHF(ISYMN),1)
C
C  bar_occ(4a)
C
            CALL DGEMM('N','N',NVIR(ISYMB),NMAABI(ISYCAI),NRHF(ISYMN),
     *                 -ONE,WORK(KOFF1),NTOTB,WORK(KOFF2),NTOTN,
     *                 ONE,T3B0JK(KOFF3),NTOTB)
C
      END DO ! ISYMN
C
C=================================================
C     Calculate (4b)   t^(ac)_(nj) g(in|kb)
C                    = t^(ca)_(jn) g(kb|in)
C                    
C                      T^J(can) I^K(nbi)
C=================================================
C
C-------------------------------
C     Sort L2TP(cjna) as T^J(can)
C-------------------------------
C
      ISYL2CAN = MULD2H(ISYML2,ISYMJ)
      ISYINTNBI = MULD2H(ISYINT,ISYMK)
C
      KL2CAN = KEND1
      KINTNBI  = KL2CAN + NMAABI(ISYL2CAN)
      KCABI   = KINTNBI + NCKI(ISYINTNBI)
      KEND1 = KCABI + NMAAB_CI(ISYT3B0JK)
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3BAR0_IJ_O (4b)')
      END IF
C
      CALL DZERO(WORK(KCABI),NMAAB_CI(ISYT3B0JK))
C
      CALL SORT_T2_ABJ(WORK(KL2CAN),ISYMJ,J,L2TP,ISYML2)
C
C----------------------------------------------
C     Sort g(kb|in) = T3OG2(b,i,k,n) as I^K(nbi)
C----------------------------------------------
C
      CALL SORT_INT_JAK_I(WORK(KINTNBI),T3BOG2,ISYINT,ISYMK,K)
C
C---------------------------------------------
C    Multiply T^J(can) * I^K(nbi) = T^JK(cabi)
C---------------------------------------------
C
      DO ISYMN = 1, NSYM
         ISYCA = MULD2H(ISYL2CAN,ISYMN)
         ISYBI = MULD2H(ISYINTNBI,ISYMN)
C
         KOFF1 = KL2CAN
     *         + IMAABI(ISYCA,ISYMN)
         KOFF2 = KINTNBI
     *         + IMAIAJ(ISYMN,ISYBI)
         KOFF3 = KCABI
     *         + IMAAB_CI(ISYCA,ISYBI)
C
         NTOTCA = MAX(NMATAB(ISYCA),1)
         NTOTN  = MAX(NRHF(ISYMN),1)
C
         CALL DGEMM('N','N',NMATAB(ISYCA),NT1AM(ISYBI),NRHF(ISYMN),
     *              ONE,WORK(KOFF1),NTOTCA,WORK(KOFF2),NTOTN,
     *              ONE,WORK(KOFF3),NTOTCA)
C
      END DO ! ISYMN
C
C-------------------------------------------
C     T3B0JK(bcai) = T3B0JK(bcai) + T^JK(cabi)
C-------------------------------------------
C
C  bar_occ(4b)
C
      CALL FBCAI(T3B0JK,WORK(KCABI),ISYT3B0JK)
C
C=================================================
C     Calculate (5a)   - t^(ca)_(kn) L(in|jb)
C                    = - t^(ca)_(kn) L(jb|in)
C                    
C                      - T^K(can) I^J(nbi)
C=================================================
C
C-------------------------------
C     Sort L2TP(ckna) as T^K(can)
C-------------------------------
C
      ISYL2CAN = MULD2H(ISYML2,ISYMK)
      ISYINTNBI = MULD2H(ISYINT,ISYMJ)
C
      KL2CAN = KEND1
      KINTNBI  = KL2CAN + NMAABI(ISYL2CAN)
      KCABI   = KINTNBI + NCKI(ISYINTNBI)
      KEND1 = KCABI + NMAAB_CI(ISYT3B0JK)
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3BAR0_IJ_O (5a)')
      END IF
C
      CALL DZERO(WORK(KCABI),NMAAB_CI(ISYT3B0JK))
C
      CALL SORT_T2_ABJ(WORK(KL2CAN),ISYMK,K,L2TP,ISYML2)
C
C----------------------------------------------
C     Sort L(jb|in) = T3OL2(b,i,j,n) as I^J(nbi)
C----------------------------------------------
C
      CALL SORT_INT_JAK_I(WORK(KINTNBI),T3BOL2,ISYINT,ISYMJ,J)
C
C------------------------------------------
C    Multiply T^K(can) * I^J(nbi) = T^JK(cabi)
C------------------------------------------
C
      DO ISYMN = 1, NSYM
         ISYCA = MULD2H(ISYL2CAN,ISYMN)
         ISYBI = MULD2H(ISYINTNBI,ISYMN)
C
         KOFF1 = KL2CAN
     *         + IMAABI(ISYCA,ISYMN)
         KOFF2 = KINTNBI
     *         + IMAIAJ(ISYMN,ISYBI)
         KOFF3 = KCABI
     *         + IMAAB_CI(ISYCA,ISYBI)
C
         NTOTCA = MAX(NMATAB(ISYCA),1)
         NTOTN  = MAX(NRHF(ISYMN),1)
C
         CALL DGEMM('N','N',NMATAB(ISYCA),NT1AM(ISYBI),NRHF(ISYMN),
     *              -ONE,WORK(KOFF1),NTOTCA,WORK(KOFF2),NTOTN,
     *              ONE,WORK(KOFF3),NTOTCA)
C
      END DO ! ISYMN
C
C-------------------------------------------
C     T3B0JK(bcai) = T3B0JK(bcai) + T^JK(cabi)
C-------------------------------------------
C
C  bar_occ(5a)
C
      CALL FBCAI(T3B0JK,WORK(KCABI),ISYT3B0JK)
C
C=================================================
C     Calculate (5b)   t^(ac)_(jn) g(kn|ib)
C                    = t^(ca)_(nj) g(ib|kn)
C                    
C                      T^J(can) I^K(nbi)
C=================================================
C
C-------------------------------
C     Sort L2TP(cnja) as T^J(can)
C-------------------------------
C
      ISYL2CAN = MULD2H(ISYML2,ISYMJ)
      ISYINTNBI = MULD2H(ISYINT,ISYMK)
C
      KL2CAN = KEND1
      KINTNBI  = KL2CAN + NMAABI(ISYL2CAN)
      KCABI   = KINTNBI + NCKI(ISYINTNBI)
      KEND1 = KCABI + NMAAB_CI(ISYT3B0JK)
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3BAR0_IJ_O (5b)')
      END IF
C
      CALL DZERO(WORK(KCABI),NMAAB_CI(ISYT3B0JK))
C
      CALL SORT_T2_ABI(WORK(KL2CAN),ISYMJ,J,L2TP,ISYML2)
C
C----------------------------------------------
C     Sort g(ib|kn) = T3OG2(b,k,i,n) as I^K(nbi)
C----------------------------------------------
C
      CALL SORT_INT_JAI_K(WORK(KINTNBI),T3BOG2,ISYINT,ISYMK,K)
C
C------------------------------------------
C    Multiply T^J(can) * I^K(nbi) = T^JK(cabi)
C------------------------------------------
C
      DO ISYMN = 1, NSYM
         ISYCA = MULD2H(ISYL2CAN,ISYMN)
         ISYBI = MULD2H(ISYINTNBI,ISYMN)
C
         KOFF1 = KL2CAN
     *         + IMAABI(ISYCA,ISYMN)
         KOFF2 = KINTNBI
     *         + IMAIAJ(ISYMN,ISYBI)
         KOFF3 = KCABI
     *         + IMAAB_CI(ISYCA,ISYBI)
C
         NTOTCA = MAX(NMATAB(ISYCA),1)
         NTOTN  = MAX(NRHF(ISYMN),1)
C
         CALL DGEMM('N','N',NMATAB(ISYCA),NT1AM(ISYBI),NRHF(ISYMN),
     *              ONE,WORK(KOFF1),NTOTCA,WORK(KOFF2),NTOTN,
     *              ONE,WORK(KOFF3),NTOTCA)
C
      END DO ! ISYMN
C
C-------------------------------------------
C     T3B0JK(bcai) = T3B0JK(bcai) + T^JK(cabi)
C-------------------------------------------
C
C  bar_occ(5b)
C
      CALL FBCAI(T3B0JK,WORK(KCABI),ISYT3B0JK)
C
C=================================================
C     Calculate (6a)   - t^(cb)_(kn) L(jn|ia)
C                    = - t^(bc)_(nk) L(ia|jn)
C                    
C                      T^K(bcn) I^J(nai)
C=================================================
C
C-------------------------------
C     Sort L2TP(bnkc) as T^K(bcn)
C-------------------------------
C
      ISYL2BCN = MULD2H(ISYML2,ISYMK)
      ISYINTNAI = MULD2H(ISYINT,ISYMJ)
C
      KL2BCN = KEND1
      KINTNAI  = KL2BCN + NMAABI(ISYL2BCN)
      KEND1   = KINTNAI + NCKI(ISYINTNAI)
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3BAR0_IJ_O (6a)')
      END IF
C
      CALL SORT_T2_ABI(WORK(KL2BCN),ISYMK,K,L2TP,ISYML2)
C
C----------------------------------------------
C     Sort L(ia|jn) = T3OL2(a,j,i,n) as I^J(nai)
C----------------------------------------------
C
      CALL SORT_INT_JAI_K(WORK(KINTNAI),T3BOL2,ISYINT,ISYMJ,J)
C
C------------------------------------------
C    Multiply T^K(bcn) * I^J(nai) = T^JK(bcai)
C------------------------------------------
C
      DO ISYMN = 1, NSYM
         ISYBC = MULD2H(ISYL2BCN,ISYMN)
         ISYAI = MULD2H(ISYINTNAI,ISYMN)
C
         KOFF1 = KL2BCN
     *         + IMAABI(ISYBC,ISYMN)
         KOFF2 = KINTNAI
     *         + IMAIAJ(ISYMN,ISYAI)
         KOFF3 = IMAAB_CI(ISYBC,ISYAI) + 1
C
         NTOTBC = MAX(NMATAB(ISYBC),1)
         NTOTN  = MAX(NRHF(ISYMN),1)
C
C  bar_occ(6a)
C
         CALL DGEMM('N','N',NMATAB(ISYBC),NT1AM(ISYAI),NRHF(ISYMN),
     *              -ONE,WORK(KOFF1),NTOTBC,WORK(KOFF2),NTOTN,
     *              ONE,T3B0JK(KOFF3),NTOTBC)
C
      END DO ! ISYMN
C
C=================================================
C     Calculate (6b)   t^(bc)_(in) g(kn|ja)
C                    = t^(cb)_(ni) g(ja|kn)
C                             
C                      T(ncbi) I^JK(an)
C=================================================
C
C-------------------------------
C     Sort L2TP(cnib) as T(ncbi)
C-------------------------------
C
      ISYJK = MULD2H(ISYMJ,ISYMK)
      ISYAN = MULD2H(ISYINT,ISYJK)
C
      KL2NCBI = KEND1
      KINTAN  = KL2NCBI + NT2SQ(ISYML2)
      KACBI   = KINTAN  + NT1AM(ISYAN)
      KEND1   = KACBI   + NMAAOBCI(ISYT3B0JK)
      LWRK1   = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3B0_JK_O (6b)')
      END IF
C
      CALL DZERO(WORK(KACBI),NMAAOBCI(ISYT3B0JK))
C
      CALL SORT_T2_IABJ(WORK(KL2NCBI),L2TP,ISYML2)
C
C------------------------------------------------
C     Sort g(ja|kn) = T3BOG2(a,k,j,n) as I^JK(an)
C------------------------------------------------
C
      CALL SORT_INT_AJ_IK(WORK(KINTAN),T3BOG2,ISYINT,ISYMJ,J,ISYMK,K)
C
C------------------------------------------
C    Multiply I^JK(an) T(ncbi) = T^JK(acbi)
C------------------------------------------
C
      DO ISYMN = 1, NSYM
            ISYMA = MULD2H(ISYAN,ISYMN)
            ISYCBI = MULD2H(ISYML2,ISYMN)
C      
            KOFF1 = KINTAN
     *            + IT1AM(ISYMA,ISYMN)
            KOFF2 = KL2NCBI
     *            + IMAJBAI(ISYMN,ISYCBI)
            KOFF3 = KACBI
     *            + IMAAOBCI(ISYMA,ISYCBI) 
C
            NTOTA = MAX(NVIR(ISYMA),1)
            NTOTN = MAX(NRHF(ISYMN),1)
C
            CALL DGEMM('N','N',NVIR(ISYMA),NMAABI(ISYCBI),NRHF(ISYMN),
     *                 ONE,WORK(KOFF1),NTOTA,WORK(KOFF2),NTOTN,
     *                 ONE,WORK(KOFF3),NTOTA)
C
      END DO ! ISYMN
C
C     T30JK(bcai) = T30JK(bcai) + T^JK(acbi)
C
C  bar_occ(6b)
C
       CALL FCBAI(T3B0JK,WORK(KACBI),ISYT3B0JK)
C
C-------------
C     End.
C-------------
C
      CALL QEXIT('GET_T3B0_JK_O')
C
      RETURN
      END
C  /* Deckntvir_t3b0_jk */
      SUBROUTINE INTVIR_T3B0_JK(XGADCK,XLADCK,ISYINT,XLAMP,
     *                         LU3VI,FN3VI,LU3FOP,FN3FOP,
     *                         WORK,LWORK) 
**********************************************************
*
*     Construvt the integrals used for t3B0^JK calculation
*
*     Read virtual integrals (kc | delta D) stored as I^D(ckdelta)
*     Transform to           (kc ! ad)      stored as I^D(cka)
*     Final sort (kc ! ad)  as I(adck) 
*
*     OUTPUT (XGADCK) : g(kcad) = (kc ! ad) sorted as I(adck) 
*     OUTPUT (XLADCK) : L(kcad) sorted as I(adck) 
*
*     P. Jorgensen, F. Pawlowski, 31-01-2003, Aarhus.
**********************************************************
C Integrals (kc ! ad) stored as I^D(cka) sorted as I(adck)
C
      IMPLICIT NONE
#include <ccsdsym.h>
#include <ccorb.h>
#include <priunit.h>
C
      INTEGER ISYINT, LU3VI, LU3FOP, LWORK
      INTEGER ISYMD, ISYCKA, KINTVI, KTRVI1, KEND1, LWRK1, IOFF
C
      CHARACTER*(*) FN3VI, FN3FOP 
C
#if defined (SYS_CRAY)
      REAL XGADCK(*), XLADCK(*), XLAMP(*), WORK(LWORK)
#else
      DOUBLE PRECISION XGADCK(*), XLADCK(*), XLAMP(*), WORK(LWORK) 
#endif
C
      CALL QENTER('INTVIR_T3B0_JK')
C
C***********************************************************'
C     Get  (XGBDCK) : g(kcad) = (kc ! ad) sorted as I(adck) 
C***********************************************************'
C
      DO ISYMD = 1, NSYM
         ISYCKA = MULD2H(ISYINT,ISYMD)
C
         KINTVI = 1
         KTRVI1 = KINTVI + NCKA(ISYCKA)
         KEND1  = KTRVI1 + NCKATR(ISYCKA)
         LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in INTVIR_T3B0_JK ')
      ENDIF
C
         DO D = 1, NVIR(ISYMD)
C
C     Read virtual integrals (kc | delta D) stored as I^D(ckdelta)
C
            IOFF = ICKAD(ISYCKA,ISYMD) + NCKA(ISYCKA)*(D - 1) + 1
            IF (NCKA(ISYCKA) .GT. 0) THEN
               CALL GETWA2(LU3VI,FN3VI,WORK(KINTVI),IOFF,
     &                        NCKA(ISYCKA))
            ENDIF
C
C     Transform to           (kc ! ad)      stored as I^D(cka)
C 
            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI1),XLAMP,
     *                       ISYMD,D,ISYMOP,WORK(KEND1),LWRK1)
C
C     Final sort (kc ! ad)  as I(adck) 
C
            CALL SORT_INTVIR_T3B0(XGADCK,WORK(KTRVI1),
     *                           D,ISYMD,ISYCKA,WORK(KEND1),LWRK1)
C
C     OUTPUT (XGADCK) : g(kcad) = (kc ! ad) sorted as I(adck) 
C
         END DO !  D
      END DO !  ISYMD
C
C
C***********************************************************'
C     Get (XLBDCK) : L(kcad) sorted as I(adck) 
C***********************************************************'
C
      DO ISYMD = 1, NSYM
         ISYCKA = MULD2H(ISYINT,ISYMD)
C
         KINTVI = 1
         KTRVI1 = KINTVI + NCKA(ISYCKA)
         KEND1  = KTRVI1 + NCKATR(ISYCKA)
         LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in INTVIR_T3B0_JK ')
      ENDIF
C
         DO D = 1, NVIR(ISYMD)
C
C     Read virtual integrals L(kc | delta D) stored as L^D(ckdelta)
C
            IOFF = ICKAD(ISYCKA,ISYMD) + NCKA(ISYCKA)*(D - 1) + 1
            IF (NCKA(ISYCKA) .GT. 0) THEN
               CALL GETWA2(LU3FOP,FN3FOP,WORK(KINTVI),IOFF,
     &                        NCKA(ISYCKA))
            ENDIF
C
C     Transform to           L(kc ! ad)      stored as L^D(cka)
C 
            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI1),XLAMP,
     *                       ISYMD,D,ISYMOP,WORK(KEND1),LWRK1)
C
C     Final sort L(kc ! ad)  as I(adck) 
C
            CALL SORT_INTVIR_T3B0(XLADCK,WORK(KTRVI1),
     *                           D,ISYMD,ISYCKA,WORK(KEND1),LWRK1)
C
C     OUTPUT (XLADCK) : L(kcad) sorted as I(adck) 
C
         END DO !  D
      END DO !  ISYMD
C
      CALL QEXIT('INTVIR_T3B0_JK')
C
      RETURN
      END
C  /* Deck sort_intvir_t3B0 */
      SUBROUTINE SORT_INTVIR_T3B0(XADCK,XCKA,
     *                           D,ISYMD,ISYCKA,WORK,LWORK)
C
**************************************************************************
C Integrals (kc ! ad) stored as I^D(cka) (XCKA) sorted as I(adck) (XADCK)
**************************************************************************
C
      IMPLICIT NONE
#include <ccsdsym.h>
#include <ccorb.h>
#include <priunit.h>
C
      INTEGER ISYMD,ISYCKA,LWORK
      INTEGER ISYMA,ISYMCK,ISYAD,ISYMC,ISYMK,ISYADC
      INTEGER KOFF1,KOFF2
C
#if defined (SYS_CRAY)
      REAL XADCK(*),XCKA(*),WORK(LWORK)
#else
      DOUBLE PRECISION XADCK(*),XCKA(*),WORK(LWORK)
#endif
C
      CALL QENTER('SORT_INTVIR_T30')
C
      DO ISYMA = 1,NSYM
         ISYMCK = MULD2H(ISYCKA,ISYMA) 
         ISYAD  =  MULD2H(ISYMA,ISYMD)
         DO ISYMC = 1,NSYM
            ISYMK = MULD2H(ISYMCK,ISYMC)
            ISYADC = MULD2H(ISYAD,ISYMC)
            DO A = 1,NVIR(ISYMA)
               DO K = 1,NRHF(ISYMK)
                  DO C = 1,NVIR(ISYMC)
                     KOFF1 = ICKATR(ISYMCK,ISYMA)
     *                     + NT1AM(ISYMCK)*(A-1)
     *                     + IT1AM(ISYMC,ISYMK)
     *                     + NVIR(ISYMC)*(K-1)
     *                     + C
                     KOFF2 = IMAABCI(ISYADC,ISYMK)
     *                     + NMAABC(ISYADC)*(K-1)
     *                     + IMAABC(ISYAD,ISYMC)
     *                     + NMATAB(ISYAD)*(C-1)
     *                     + IMATAB(ISYMA,ISYMD)
     *                     + NVIR(ISYMA)*(D-1)
     *                     + A
                    XADCK(KOFF2) = XCKA(KOFF1)
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      CALL QEXIT('SORT_INTVIR_T30')
      RETURN
      END
C  /* Deck get_T3B0_JK_V */
      SUBROUTINE GET_T3B0_JK_V(T30JK,ISYT30JK,L2TP,
     *                           ISYML2,XGADCK,XLADCK,
     *                           ISYINT,ISYMJ,J,ISYMK,K,
     *                           WORK,LWORK)

***********************************************************
*     XGADCK : g(kcad) = (kc ! ad) sorted as I(adck) 
*     XLADCK : L(kcad) sorted as I(adck) 
*    
*     T30KL sitting as (bcai)
***********************************************************
C
C The T3B0^(abc)_(iJK) contribution with integrals having 3 virtual indices
C
C     T3B0^(abc)_(iJK) = T3B0^(abc)_(iJK) +
C     P(ai,bj,ck) (sum_d l^ad_ij L(dbkc) - l^ad_jk g(ibdc) 
C
C 1)  l^ad_ij L(dbkc) - l^ad_jk g(ibdc) 
C
C 2)  l^bd_ji L(dakc) - l^bd_ik g(jadc) 
C
C 3)  l^bd_jk L(dcia) - l^bd_ki g(jcda) 
C
C 4)  l^ad_ik L(dcjb) - l^ad_kj g(icdb) 
C
C 5)  l^cd_ki L(dajb) - l^cd_ij g(kadb) 
C
C 6)  l^cd_kj L(dbia) - l^cd_ji g(kbda) 
C
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYT30JK, ISYML2, ISYINT, ISYMJ, ISYMK, LWORK
      INTEGER ISYMDAI, ISYMBCD, ISYMBCAI, ISYMDA, ISYMCBI, ISYMBC 
      INTEGER ISYMBCA, ISYMAD, ISYMDBI , ISYMABCI 
      INTEGER ISYMACD, ISYMACBI, ISYMDB, ISYMAC, ISYMACB 
      INTEGER ISYMJK, ISYMBD, ISYMCAI, ISYMDCI, ISYMBAD, ISYMBACI 
      INTEGER ISYMDC, ISYMBA, ISYMBAC, ISYMKJ, ISYMCD, ISYMCBAI 
      INTEGER ISYMI, ISYMD, ISYMA, ISYMB, ISYMC, ISYMBCI
      INTEGER NTOTBC, NTOTD, NTOTAC, NTOTB, NTOTBA, NTOTC, NTOTA 
      INTEGER KDAI, KBCD, KBCAI, KEND1, LWRK1, KAD, KDBI, KACD, KACBI
      INTEGER KBD, KDCI, KBAD, KBACI, KCD, KCBAI, KABCI 
      INTEGER KOFF1, KOFF2, KOFF3
      INTEGER ISYMBAI
C
#if defined (SYS_CRAY)
      REAL T30JK(*), L2TP(*), XGADCK(*), XLADCK(*), WORK(LWORK)
      REAL ONE
#else
      DOUBLE PRECISION T30JK(*), L2TP(*), XGADCK(*), XLADCK(*)
      DOUBLE PRECISION WORK(LWORK) 
      DOUBLE PRECISION ONE
#endif
C
      PARAMETER (ONE = 1.0D0)
C
      CALL QENTER('GET_T3B0_JK_V')
C
C***************************************************
C 1)  l^ad_ij L(dbkc) - l^ad_jk g(ibdc) 
C***************************************************
C
C l2tp(djia) =   I^J(dai) 
C
C L(dbkc) = I(dbck)  sorted as I^K(bcd)
C
C T^JK(bcai) = T^JK(bcai) + I^K(bcd)*I^J(dai)
C
C symmetry and work allocation
C

      ISYMDAI = MULD2H(ISYML2,ISYMJ)
      ISYMBCD = MULD2H(ISYINT,ISYMK)
      ISYMBCAI = MULD2H(ISYMBCD,ISYMDAI)
      KDAI  = 1
      KBCD  = KDAI  + NMAABI(ISYMDAI)
      KBCAI = KBCD  + NMAABC(ISYMBCD)
      KEND1 = KBCAI + NMAABCI(ISYMBCAI)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3B0_JK_V (1L)')
      END IF
C
      CALL DZERO(WORK(KBCAI),NMAABCI(ISYMBCAI))
C
C
C  sort l^ad_ij = l2tp(djia) as I^J(dai) 
C
      CALL SORT_T2_ABJ(WORK(KDAI),ISYMJ,J,L2TP,ISYML2)
C
C L(dbkc) = I(dbck)  sorted as I^K(bcd)
C 
      CALL SORT_INT_BCA(WORK(KBCD),ISYMK,K,XLADCK,ISYINT)
C
C T^JK(bcai) = T^JK(bcai) + I^K(bcd)*I^J(dai) 
C
      DO ISYMI = 1,NSYM
         ISYMDA = MULD2H(ISYMDAI,ISYMI)
         DO I =  1,NRHF(ISYMI)
            DO ISYMD = 1,NSYM
               ISYMBC  = MULD2H(ISYMBCD,ISYMD)
               ISYMA   = MULD2H(ISYMDA,ISYMD)
               ISYMBCA = MULD2H(ISYMBC,ISYMA)
               KOFF1   = KBCD + IMAABC(ISYMBC,ISYMD) 
               KOFF2   = KDAI 
     *                + IMAABI(ISYMDA,ISYMI)  
     *                + NMATAB(ISYMDA)*(I-1)
     *                + IMATAB(ISYMD,ISYMA)
               KOFF3   = KBCAI + IMAABCI(ISYMBCA,ISYMI)
     *                + NMAABC(ISYMBCA)*(I-1)
     *                + IMAABC(ISYMBC,ISYMA)
C  
               NTOTBC = MAX(1,NMATAB(ISYMBC))
               NTOTD  = MAX(1,NVIR(ISYMD))
C
               CALL DGEMM('N','N',NMATAB(ISYMBC),NVIR(ISYMA),
     *                    NVIR(ISYMD),ONE,WORK(KOFF1),NTOTBC,
     *                    WORK(KOFF2),NTOTD, 
     *                    ONE,WORK(KOFF3),NTOTBC) 
            END DO
         END DO
      END DO
C
C  bar_vir(1L)
C
      CALL FAB_CI(T30JK,WORK(KBCAI),ISYT30JK)
C

C
C  - l^ad_jk g(ibdc) 
C
C  L2TP(aJKd) = I^JK(ad)
C 
C g(ibdc) = I(dcbi) 
C
C T^JK(bcai) = T^JK(bcai) - I^JK(ad) * I(dcbi)
C
C
C symmetry and work allocation
C
      ISYMJK = MULD2H(ISYMJ,ISYMK)
      ISYMAD = MULD2H(ISYML2,ISYMJK)
      ISYMACBI = MULD2H(ISYINT,ISYMAD)
      KAD    = 1
      KACBI  = KAD + NMATAB(ISYMAD)
      KEND1  = KACBI + NMAABCI(ISYMACBI) 
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3B0_JK_V (1g)')
      END IF
C
      CALL DZERO(WORK(KACBI),NMAABCI(ISYMACBI))
C
C  L2TP(aJKd) = I^JK(ad)
C
      CALL SORT_T2_AB(WORK(KAD),ISYMJ,J,ISYMK,K,L2TP,ISYML2)
C
C T^JK(bcai) = T^JK(bcai) +  I^JK(ad) * I(dcbi)
C
      DO ISYMD = 1,NSYM
         ISYMA   = MULD2H(ISYMAD,ISYMD)
         ISYMCBI = MULD2H(ISYMD,ISYINT)
         KOFF1   = KAD + IMATAB(ISYMA,ISYMD)
         KOFF2   = 1  + IMAAOBCI(ISYMD,ISYMCBI)
         KOFF3   = KACBI  + IMAAOBCI(ISYMA,ISYMCBI)
C
         NTOTD  = MAX(1,NVIR(ISYMD))
         NTOTA  = MAX(1,NVIR(ISYMA))
C
         CALL DGEMM('N','N',NVIR(ISYMA),NMAABI(ISYMCBI),
     *                    NVIR(ISYMD),-ONE,WORK(KOFF1),NTOTA,
     *                    XGADCK(KOFF2),NTOTD,
     *                    ONE,WORK(KOFF3),NTOTD)
      END DO
C
C T^JK(bcai) = T^JK(bcai) +  I^JK(ad) * I(dcbi)
C
C  bar_vir(1g)
C
      CALL FCBAI(T30JK,WORK(KACBI),ISYT30JK)
C
C***************************************************
C 2)  l^bd_ji L(dakc) - l^bd_ik g(jadc)
C***************************************************
C
C L2TP(bjid)  =   I^J(dbi)
C
C  L(dakc) = I(dack)  sorted as I^K(acd)
C
C T^JK(bcai) = T^JK(bcai) +  I^K(acd) * I^J(dbi)
C
C symmetry and work allocation
C
      ISYMDBI = MULD2H(ISYML2,ISYMJ)
      ISYMACD = MULD2H(ISYINT,ISYMK)
      ISYMACBI = MULD2H(ISYMACD,ISYMDBI)
      KDBI   = 1
      KACD   = KDBI  +  NMAABI(ISYMDBI)
      KACBI  = KACD  +  NMAABC(ISYMACD)
      KEND1  = KACBI +  NMAABCI(ISYMACBI)
      LWRK1  = LWORK -  KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3B0_JK_V (2L)')
      END IF
C
      CALL DZERO(WORK(KACBI),NMAABCI(ISYMACBI))
C
C L2TP(bjid)  =   I^J(dbi)
C
      CALL SORT_T2_BAJ(WORK(KDBI),ISYMJ,J,L2TP,ISYML2)
C
C  L(dakc) = I(dack)  sorted as I^K(acd)
C
      CALL SORT_INT_BCA(WORK(KACD),ISYMK,K,XLADCK,ISYINT)
C
C
      DO ISYMI = 1,NSYM
         ISYMDB = MULD2H(ISYMDBI,ISYMI)
         DO I =  1,NRHF(ISYMI)
            DO ISYMD = 1,NSYM
               ISYMB   = MULD2H(ISYMDB,ISYMD)
               ISYMAC  = MULD2H(ISYMACD,ISYMD)
               ISYMACB = MULD2H(ISYMAC,ISYMB)
               KOFF1   = KACD
     *                + IMAABC(ISYMAC,ISYMD)
               KOFF2   = KDBI
     *                + IMAABI(ISYMDB,ISYMI)
     *                + NMATAB(ISYMDB)*(I-1)
     *                + IMATAB(ISYMD,ISYMB)
               KOFF3   = KACBI
     *                + IMAABCI(ISYMACB,ISYMI)
     *                + NMAABC(ISYMACB)*(I-1)
     *                + IMAABC(ISYMAC,ISYMB)
C 
               NTOTAC = MAX(1,NMATAB(ISYMAC))
               NTOTD  = MAX(1,NVIR(ISYMD))
C
C  work(acbi) =  I^K(acd) * I^J(dbi)
C
               CALL DGEMM('N','N',NMATAB(ISYMAC),NVIR(ISYMB),
     *                    NVIR(ISYMD),ONE,WORK(KOFF1),NTOTAC,
     *                    WORK(KOFF2),NTOTD,
     *                    ONE,WORK(KOFF3),NTOTAC)
C
            END DO
         END DO
      END DO
C
C T^JK(bcai) = T^JK(bcai) +  work(acbi)
C
C
C  bar_vir(2L)
C
      CALL FCBAI(T30JK,WORK(KACBI),ISYT30JK)
C
C  - l^bd_ik g(jadc)
C
C L2TP(dKib) = I^K(dbi)
C
C g(jadc) = I(dcaj) sorted as I^J(acd)
C
C T^JK(bcai) = T^JK(bcai) +  I^J(acd) * I^K(dbi) 
C

      ISYMDBI = MULD2H(ISYML2,ISYMK)
      ISYMACD = MULD2H(ISYINT,ISYMJ)
      ISYMACBI = MULD2H(ISYMACD,ISYMDBI)
      KDBI  = 1
      KACD  = KDBI  + NMAABI(ISYMDBI)
      KACBI = KACD  + NMAABC(ISYMACD)
      KEND1  = KACBI + NMAABCI(ISYMACBI)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3B0_JK_V (2g)')
      END IF
C
      CALL DZERO(WORK(KACBI),NMAABCI(ISYMACBI))
C
C L2TP(dKib) = I^K(dbi)
C 
      CALL SORT_T2_ABJ(WORK(KDBI),ISYMK,K,L2TP,ISYML2)
C
C g(jadc) = I(dcaj) sorted as I^J(acd)
C
      CALL SORT_INT_CBA(WORK(KACD),ISYMJ,J,XGADCK,ISYINT)
C
      DO ISYMI = 1,NSYM
         ISYMDB = MULD2H(ISYMDBI,ISYMI)
         DO I =  1,NRHF(ISYMI)
            DO ISYMD = 1,NSYM
               ISYMAC  = MULD2H(ISYMACD,ISYMD)
               ISYMB   = MULD2H(ISYMDB,ISYMD)
               ISYMACB = MULD2H(ISYMAC,ISYMB)
               KOFF1   = KACD + IMAABC(ISYMAC,ISYMD)
               KOFF2   = KDBI
     *                + IMAABI(ISYMDB,ISYMI)
     *                + NMATAB(ISYMDB)*(I-1)
     *                + IMATAB(ISYMD,ISYMB)
               KOFF3   = KACBI + IMAABCI(ISYMACB,ISYMI)
     *                + NMAABC(ISYMACB)*(I-1)
     *                + IMAABC(ISYMAC,ISYMB)
C
               NTOTAC = MAX(1,NMATAB(ISYMAC))
               NTOTD  = MAX(1,NVIR(ISYMD))
C
C work(acbi) = -  I^J(acd) * I^K(dbi) 
C
               CALL DGEMM('N','N',NMATAB(ISYMAC),NVIR(ISYMB),
     *                    NVIR(ISYMD),-ONE,WORK(KOFF1),NTOTAC,
     *                    WORK(KOFF2),NTOTD,
     *                    ONE,WORK(KOFF3),NTOTAC)
            END DO
         END DO
      END DO
C
C T^JK(bcai) = T^JK(bcai) + work(acbi) 
C
C  bar_vir(2g)
C
      CALL FCBAI(T30JK,WORK(KACBI),ISYT30JK)
C
C**************************************************    
C 3)  l^bd_jk L(dcia) - l^bd_ki g(jcda)
C**************************************************    
C
C L2TP(bjkd)  =   I^JK(bd)
C
C  L(dcia) = I(dcai)
C
C T^JK(bcai) = T^JK(bcai) +  I^JK(bd) * I(dcai) 
C
C symmetry and work allocation
C
      ISYMJK = MULD2H(ISYMJ,ISYMK)
      ISYMBD = MULD2H(ISYML2,ISYMJK)
      KBD    = 1
      KEND1  = KBD + NMATAB(ISYMBD)
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3B0_JK_V (3L)')
      END IF

C
C L2TP(bjkd)  =   I^JK(bd)
C
      CALL SORT_T2_AB(WORK(KBD),ISYMJ,J,ISYMK,K,L2TP,ISYML2)
C
C T^JK(bcai) = T^JK(bcai) +  I^JK(bd) * I(dcai) 
C
      DO ISYMD = 1,NSYM
         ISYMB   = MULD2H(ISYMBD,ISYMD)
         ISYMCAI = MULD2H(ISYMD,ISYINT)
         KOFF1   = KBD + IMATAB(ISYMB,ISYMD)
         KOFF2   = 1  + IMAAOBCI(ISYMD,ISYMCAI)
         KOFF3   = 1  + IMAAOBCI(ISYMB,ISYMCAI)
C 
         NTOTD  = MAX(1,NVIR(ISYMD))
         NTOTB  = MAX(1,NVIR(ISYMB))
C
C
C  bar_vir(3L)
C
         CALL DGEMM('N','N',NVIR(ISYMB),NMAABI(ISYMCAI),
     *                    NVIR(ISYMD),ONE,WORK(KOFF1),NTOTB,
     *                    XLADCK(KOFF2),NTOTD,
     *                    ONE,T30JK(KOFF3),NTOTD)
      END DO
C
C - l^bd_ki g(jcda)
C
C L2TP(bkid) = I^K(dbi)
C
C g(jcda) = I(dacj)  sorted as I^J(acd)
C
C T^JK(bcai) = T^JK(bcai) - I^J(acd) * I^K(dbi)
C
C symmetry and work allocation
C
      ISYMDBI = MULD2H(ISYML2,ISYMK)
      ISYMACD = MULD2H(ISYINT,ISYMJ)
      ISYMACBI = MULD2H(ISYMACD,ISYMDBI)
      KDBI   = 1
      KACD   = KDBI  +  NMAABI(ISYMDBI)
      KACBI  = KACD  +  NMAABC(ISYMACD)
      KEND1  = KACBI +  NMAABCI(ISYMACBI)
      LWRK1  = LWORK -  KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3B0_JK_V (3g)')
      END IF
C
      CALL DZERO(WORK(KACBI),NMAABCI(ISYMACBI))
C
C L2TP(bKid)  =   I^K(dbi)
C
      CALL SORT_T2_BAJ(WORK(KDBI),ISYMK,K,L2TP,ISYML2)
C
C g(jcda) = I(dacj)  sorted as I^J(acd)
C
      CALL SORT_INT_BCA(WORK(KACD),ISYMJ,J,XGADCK,ISYINT)
C
C
C
      DO ISYMI = 1,NSYM
         ISYMDB = MULD2H(ISYMDBI,ISYMI)
         DO I =  1,NRHF(ISYMI)
            DO ISYMD = 1,NSYM
               ISYMB   = MULD2H(ISYMDB,ISYMD)
               ISYMAC  = MULD2H(ISYMACD,ISYMD)
               ISYMACB = MULD2H(ISYMAC,ISYMB)
               KOFF1   = KACD
     *                + IMAABC(ISYMAC,ISYMD)
               KOFF2   = KDBI
     *                + IMAABI(ISYMDB,ISYMI)
     *                + NMATAB(ISYMDB)*(I-1)
     *                + IMATAB(ISYMD,ISYMB)
               KOFF3   = KACBI
     *                + IMAABCI(ISYMACB,ISYMI)
     *                + NMAABC(ISYMACB)*(I-1)
     *                + IMAABC(ISYMAC,ISYMB)
C
               NTOTAC = MAX(1,NMATAB(ISYMAC))
               NTOTD  = MAX(1,NVIR(ISYMD))
C
C  work(acbi) =  I^J(acd) * I^K(dbi)
C
               CALL DGEMM('N','N',NMATAB(ISYMAC),NVIR(ISYMB),
     *                    NVIR(ISYMD),-ONE,WORK(KOFF1),NTOTAC,
     *                    WORK(KOFF2),NTOTD,
     *                    ONE,WORK(KOFF3),NTOTAC)
C
            END DO
         END DO
      END DO
C
C T^JK(bcai) = T^JK(bcai) +  work(acbi)
C
C  bar_vir(3g)
C
      CALL FCBAI(T30JK,WORK(KACBI),ISYT30JK)
C

C****************************************************
C 4)  l^ad_ik L(dcjb) - l^ad_kj g(icdb)
C****************************************************
C
C L2TP(dkia) =   I^K(dai)
C
C L(dcjb) = I(dcbj) stored as I^J(bcd)
C
C T^JK(bcai) = T^JK(bcai) + I^J(bcd)*I^K(dai)
C
C symmetry and work allocation
C
      ISYMDAI = MULD2H(ISYML2,ISYMK)
      ISYMBCD = MULD2H(ISYINT,ISYMJ)
      ISYMBCAI = MULD2H(ISYMBCD,ISYMDAI)
      KDAI  = 1
      KBCD  = KDAI  + NMAABI(ISYMDAI)
      KBCAI = KBCD  + NMAABC(ISYMBCD)
      KEND1 = KBCAI + NMAABCI(ISYMBCAI)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3B0_JK_V (4L)')
      END IF
C
      CALL DZERO(WORK(KBCAI),NMAABCI(ISYMBCAI))
C
C L2TP(dkia) =   I^K(dai)
C
      CALL SORT_T2_ABJ(WORK(KDAI),ISYMK,K,L2TP,ISYML2)
C
C L(dcjb) = I(dcbj) sorted as I^J(bcd)
C
      CALL SORT_INT_CBA(WORK(KBCD),ISYMJ,J,XLADCK,ISYINT)
C
C T^JK(bcai) = T^JK(bcai) + I^J(bcd)*I^K(dai)
C
      DO ISYMI = 1,NSYM
         ISYMDA = MULD2H(ISYMDAI,ISYMI)
         DO I =  1,NRHF(ISYMI)
            DO ISYMD = 1,NSYM
               ISYMBC  = MULD2H(ISYMBCD,ISYMD)
               ISYMA   = MULD2H(ISYMDA,ISYMD)
               ISYMBCA = MULD2H(ISYMBC,ISYMA)
               KOFF1   = KBCD + IMAABC(ISYMBC,ISYMD)
               KOFF2   = KDAI
     *                + IMAABI(ISYMDA,ISYMI)
     *                + NMATAB(ISYMDA)*(I-1)
     *                + IMATAB(ISYMD,ISYMA)
               KOFF3   = KBCAI + IMAABCI(ISYMBCA,ISYMI)
     *                + NMAABC(ISYMBCA)*(I-1)
     *                + IMAABC(ISYMBC,ISYMA)
C 
               NTOTBC = MAX(1,NMATAB(ISYMBC))
               NTOTD  = MAX(1,NVIR(ISYMD))
C
               CALL DGEMM('N','N',NMATAB(ISYMBC),NVIR(ISYMA),
     *                    NVIR(ISYMD),ONE,WORK(KOFF1),NTOTBC,
     *                    WORK(KOFF2),NTOTD,
     *                    ONE,WORK(KOFF3),NTOTBC)
            END DO
         END DO
      END DO
C
C  bar_vir(4L)
C
      CALL FAB_CI(T30JK,WORK(KBCAI),ISYT30JK)
C
C  - l^ad_kj g(icdb)
C
C  L2TP(akjd) = I^KJ(ad)
C
C g(icdb) = I(dbci)
C
C T^JK(bcai) = T^JK(bcai) + I^KJ(ad) * I(dbci)
C
C symmetry and work allocation
C
      ISYMJK = MULD2H(ISYMJ,ISYMK)
      ISYMAD = MULD2H(ISYML2,ISYMJK)
      ISYMABCI = MULD2H(ISYINT,ISYMAD)
      KAD    = 1
      KABCI  = KAD + NMATAB(ISYMAD)
      KEND1  = KABCI + NMAABCI(ISYMABCI)
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3B0_JK_V (4g)')
      END IF
C
      CALL DZERO(WORK(KABCI),NMAABCI(ISYMABCI))
C
C  L2TP(akjd) = I^KJ(ad)
C
      CALL SORT_T2_AB(WORK(KAD),ISYMK,K,ISYMJ,J,L2TP,ISYML2)
C
C T^JK(bcai) = T^JK(bcai) + I^KJ(ad) * I(dbci)
C
      DO ISYMD = 1,NSYM
         ISYMA   = MULD2H(ISYMAD,ISYMD)
         ISYMBCI = MULD2H(ISYMD,ISYINT)
         KOFF1   = KAD + IMATAB(ISYMA,ISYMD)
         KOFF2   = 1  + IMAAOBCI(ISYMD,ISYMBCI)
         KOFF3   = KABCI  + IMAAOBCI(ISYMA,ISYMBCI)
C
         NTOTD  = MAX(1,NVIR(ISYMD))
         NTOTA  = MAX(1,NVIR(ISYMA))
C
C  work(abci) = I^KJ(ad) * I(dbci)
C
         CALL DGEMM('N','N',NVIR(ISYMA),NMAABI(ISYMBCI),
     *                    NVIR(ISYMD),-ONE,WORK(KOFF1),NTOTA,
     *                    XGADCK(KOFF2),NTOTD,
     *                    ONE,WORK(KOFF3),NTOTA)
      END DO
C
C T^JK(bcai) = T^JK(bcai) + work(abci) 
C
C  bar_vir(4g)
C
      CALL FCABI(T30JK,WORK(KABCI),ISYT30JK)
C

C**************************************************
C 5)  l^cd_ki L(dajb) - l^cd_ij g(kadb)
C**************************************************
C
C L2TP(ckid)  =   I^K(dci)
C
C  L(dajb) = I(dabj) sorted as I^J(bad)
C
C T^JK(bcai) = T^JK(bcai) + I^J(bad) * I^K(dci) 
C
C symmetry and work allocation
C
      ISYMDCI = MULD2H(ISYML2,ISYMK)
      ISYMBAD = MULD2H(ISYINT,ISYMJ)
      ISYMBACI = MULD2H(ISYMBAD,ISYMDCI)
      KDCI   = 1
      KBAD   = KDCI  +  NMAABI(ISYMDCI)
      KBACI  = KBAD  +  NMAABC(ISYMBAD)
      KEND1  = KBACI +  NMAABCI(ISYMBACI)
      LWRK1  = LWORK -  KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3B0_JK_V (5L)')
      END IF
C
      CALL DZERO(WORK(KBACI),NMAABCI(ISYMBACI))
C
C L2TP(ckid)  =   I^K(dci)
C
      CALL SORT_T2_BAJ(WORK(KDCI),ISYMK,K,L2TP,ISYML2)
C
C  L(dajb) = I(dabj) sorted as I^J(bad)
C
      CALL SORT_INT_CBA(WORK(KBAD),ISYMJ,J,XLADCK,ISYINT)
C
      DO ISYMI = 1,NSYM
         ISYMDC = MULD2H(ISYMDCI,ISYMI)
         DO I =  1,NRHF(ISYMI)
            DO ISYMD = 1,NSYM
               ISYMC   = MULD2H(ISYMDC,ISYMD)
               ISYMBA  = MULD2H(ISYMBAD,ISYMD)
               ISYMBAC = MULD2H(ISYMBA,ISYMC)
               KOFF1   = KBAD
     *                + IMAABC(ISYMBA,ISYMD)
               KOFF2   = KDCI
     *                + IMAABI(ISYMDC,ISYMI)
     *                + NMATAB(ISYMDC)*(I-1)
     *                + IMATAB(ISYMD,ISYMC)
               KOFF3   = KBACI
     *                + IMAABCI(ISYMBAC,ISYMI)
     *                + NMAABC(ISYMBAC)*(I-1)
     *                + IMAABC(ISYMBA,ISYMC)
C
               NTOTBA = MAX(1,NMATAB(ISYMBA))
               NTOTD  = MAX(1,NVIR(ISYMD))
C
C   WORK(baci) =  I^J(bad) * I^K(dci) 
               CALL DGEMM('N','N',NMATAB(ISYMBA),NVIR(ISYMC),
     *                    NVIR(ISYMD),ONE,WORK(KOFF1),NTOTBA,
     *                    WORK(KOFF2),NTOTD,
     *                    ONE,WORK(KOFF3),NTOTBA)
             END DO
         END DO
      END DO
C
C T^JK(bcai) = T^JK(bcai) + WORK(baci)  
C
C
C  bar_vir(5L)
C
      CALL FACBI(T30JK,WORK(KBACI),ISYT30JK)
C
C  - l^cd_ij g(kadb)
C
C L2TP(dJic) = I^J(dci)
C
C g(kadb) = I(dbak) sorted as I^K(bad)
C
C  T^JK(bcai) = T^JK(bcai) - I^K(bad)*I^J(dci)
C 
C
C symmetry and work allocation
C
      ISYMDCI = MULD2H(ISYML2,ISYMJ)
      ISYMBAD = MULD2H(ISYINT,ISYMK)
      ISYMBACI = MULD2H(ISYMBAD,ISYMDCI)
      KDCI   = 1
      KBAD   = KDCI  +  NMAABI(ISYMDCI)
      KBACI  = KBAD  +  NMAABC(ISYMBAD)
      KEND1  = KBACI +  NMAABCI(ISYMBACI)
      LWRK1  = LWORK -  KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3B0_JK_V (5g)')
      END IF
C
      CALL DZERO(WORK(KBACI),NMAABCI(ISYMBACI))
C
C L2TP(dJic) = I^J(dci)
C
      CALL SORT_T2_ABJ(WORK(KDCI),ISYMJ,J,L2TP,ISYML2)
C
C g(kadb) = I(dbak) sorted as I^K(bad)
C
      CALL SORT_INT_BCA(WORK(KBAD),ISYMK,K,XGADCK,ISYINT)
C
      DO ISYMI = 1,NSYM
         ISYMDC = MULD2H(ISYMDCI,ISYMI)
         DO I =  1,NRHF(ISYMI)
            DO ISYMD = 1,NSYM
               ISYMC   = MULD2H(ISYMDC,ISYMD)
               ISYMBA  = MULD2H(ISYMBAD,ISYMD)
               ISYMBAC = MULD2H(ISYMBA,ISYMC)
               KOFF1   = KBAD
     *                + IMAABC(ISYMBA,ISYMD)
               KOFF2   = KDCI
     *                + IMAABI(ISYMDC,ISYMI)
     *                + NMATAB(ISYMDC)*(I-1)
     *                + IMATAB(ISYMD,ISYMC)
               KOFF3   = KBACI
     *                + IMAABCI(ISYMBAC,ISYMI)
     *                + NMAABC(ISYMBAC)*(I-1)
     *                + IMAABC(ISYMBA,ISYMC)
C
               NTOTBA = MAX(1,NMATAB(ISYMBA))
               NTOTD  = MAX(1,NVIR(ISYMD))
C
C   WORK(baci) =  I^K(bad) * I^J(dci)
               CALL DGEMM('N','N',NMATAB(ISYMBA),NVIR(ISYMC),
     *                    NVIR(ISYMD),-ONE,WORK(KOFF1),NTOTBA,
     *                    WORK(KOFF2),NTOTD,
     *                    ONE,WORK(KOFF3),NTOTBA)
             END DO
         END DO
      END DO
C
C T^JK(bcai) = T^JK(bcai) + WORK(baci)
C
C
C  bar_vir(5g)
C
      CALL FACBI(T30JK,WORK(KBACI),ISYT30JK)
C
C**************************************************
C 6)  l^cd_kj L(dbia) - l^cd_ji g(kbda)
C****************************************
C
C L2TP(ckjd)  =   I^KJ(cd)
C
C L(dbia) = I(dbai) 
C
C T^JK(bcai) = T^JK(bcai) + I^KJ(cd) * I(dbai) 
C
C symmetry and work allocation
C
      ISYMKJ  = MULD2H(ISYMK,ISYMJ)
      ISYMCD  = MULD2H(ISYML2,ISYMKJ)
      ISYMCBAI = MULD2H(ISYINT,ISYMCD)
      KCD    = 1
      KCBAI  = KCD   +  NMATAB(ISYMCD)
      KEND1  = KCBAI +  NMAABCI(ISYMCBAI)
      LWRK1  = LWORK -  KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3B0_JK_V (6L)')
      END IF
C
      CALL DZERO(WORK(KCBAI),NMAABCI(ISYMCBAI))
C
C L2TP(ckjd)  =   I^KJ(cd)
C
      CALL SORT_T2_AB(WORK(KCD),ISYMK,K,ISYMJ,J,L2TP,ISYML2)
C
      DO ISYMD = 1,NSYM
         ISYMC   = MULD2H(ISYMCD,ISYMD)
         ISYMBAI = MULD2H(ISYMD,ISYINT)
         KOFF1   = KCD   + IMATAB(ISYMC,ISYMD)
         KOFF2   = 1     + IMAAOBCI(ISYMD,ISYMBAI)
         KOFF3   = KCBAI + IMAAOBCI(ISYMC,ISYMBAI)
C
         NTOTD  = MAX(1,NVIR(ISYMD))
         NTOTC  = MAX(1,NVIR(ISYMC))
C
C  work(cbai) =  I^KJ(cd) * I(dbai) 
         CALL DGEMM('N','N',NVIR(ISYMC),NMAABI(ISYMBAI),
     *                    NVIR(ISYMD),ONE,WORK(KOFF1),NTOTC,
     *                    XLADCK(KOFF2),NTOTD,
     *                    ONE,WORK(KOFF3),NTOTC)
      END DO
C
C T^JK(bcai) = T^JK(bcai) + work(cbai) 
C
C  bar_vir(6L)
C
      CALL FBACI(T30JK,WORK(KCBAI),ISYT30JK)
C
C
C   - l^cd_ji g(kbda)
C
C L2TP(cJid) = I^J(dci)
C
C g(kbda) = I(dabk) sorted as I^K(bad)
C
C T^JK(bcai) = T^JK(bcai) - I^K(bad) * I^J(dci)
C
C
C symmetry and work allocation
C
      ISYMDCI = MULD2H(ISYML2,ISYMJ)
      ISYMBAD = MULD2H(ISYINT,ISYMK)
      ISYMBACI = MULD2H(ISYMBAD,ISYMDCI)
      KDCI   = 1
      KBAD   = KDCI  +  NMAABI(ISYMDCI)
      KBACI  = KBAD  +  NMAABC(ISYMBAD)
      KEND1  = KBACI +  NMAABCI(ISYMBACI)
      LWRK1  = LWORK -  KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3B0_JK_V (6g)')
      END IF
C
      CALL DZERO(WORK(KBACI),NMAABCI(ISYMBACI))
C
C L2TP(cJid)  =   I^J(dci)
C
      CALL SORT_T2_BAJ(WORK(KDCI),ISYMJ,J,L2TP,ISYML2)
C
C g(kbda) = I(dabk) sorted as I^K(bad)
C
      CALL SORT_INT_CBA(WORK(KBAD),ISYMK,K,XGADCK,ISYINT)
C
C
      DO ISYMI = 1,NSYM
         ISYMDC = MULD2H(ISYMDCI,ISYMI)
         DO I =  1,NRHF(ISYMI)
            DO ISYMD = 1,NSYM
               ISYMC   = MULD2H(ISYMDC,ISYMD)
               ISYMBA  = MULD2H(ISYMBAD,ISYMD)
               ISYMBAC = MULD2H(ISYMBA,ISYMC)
               KOFF1   = KBAD
     *                + IMAABC(ISYMBA,ISYMD)
               KOFF2   = KDCI
     *                + IMAABI(ISYMDC,ISYMI)
     *                + NMATAB(ISYMDC)*(I-1)
     *                + IMATAB(ISYMD,ISYMC)
               KOFF3   = KBACI
     *                + IMAABCI(ISYMBAC,ISYMI)
     *                + NMAABC(ISYMBAC)*(I-1)
     *                + IMAABC(ISYMBA,ISYMC)
C
               NTOTBA = MAX(1,NMATAB(ISYMBA))
               NTOTD  = MAX(1,NVIR(ISYMD))
C
C   WORK(baci) =  I^K(bad) * I^J(dci)
               CALL DGEMM('N','N',NMATAB(ISYMBA),NVIR(ISYMC),
     *                    NVIR(ISYMD),-ONE,WORK(KOFF1),NTOTBA,
     *                    WORK(KOFF2),NTOTD,
     *                    ONE,WORK(KOFF3),NTOTBA)
             END DO
         END DO
      END DO
C
C T^JK(bcai) = T^JK(bcai) + WORK(baci)
C
C
C  bar_vir(6g)
C
      CALL FACBI(T30JK,WORK(KBACI),ISYT30JK)
C
      CALL QEXIT('GET_T3B0_JK_V')
C
      RETURN   
      END
C  /* Deck wbx_jk_eta */
      SUBROUTINE WBX_JK_ETA(TB0JK,ISTB0JK,XOP,ISYMXOP,WBXJK,ISYMWBX,
     *                      L2TP,ISYML2,ISYMJ,J,ISYMK,K,WORK,LWORK)
C
C WBX_JK(bcai) = WBX_JK(bcai) 
C
C 1)            - xop(da) tb0_jk(bcdi)
C                             
C 2)            + xop(il) tb0_jk(bcal) 
C                             
C 3)            + xop(j,b)*t2(ai,ck) - xop(k,b)*t2(ai,cj) 
C               + xop(k,c)*t2(ai,bj) - xop(j,c)*t2(ai,bk) 
C
C 1)            - xop(da) tb0_jk(bcdi)
C
C SORT VIR-VIR  XOP ELEMENTS (D,A)
C

      IMPLICIT NONE
C
      INTEGER ISTB0JK, ISYMXOP, ISYMWBX, ISYML2, ISYMJ, ISYMK
      INTEGER KDA, KIL, KEND1 , LWRK1, KOFF1, KOFF2, KOFF3, LWORK
      INTEGER ISYMA, ISYMD, ISYMI, ISYMBCD, ISYMBCA, ISYML 
      INTEGER ISYMAI, ISYMBC, ISYMC, ISYMCJ, ISYMCK, ISYMB
      INTEGER ISYMBK, ISYMAICK, ISYMAICJ, ISYMAIBJ, ISYMAIBK, ISYMAIJ
      INTEGER ISYMAIK, ISYMBJ
      INTEGER NTOTBC, NTOTD, NTOTBCA, NTOTI 
C
#if defined (SYS_CRAY)
      REAL TB0JK(*), L2TP(*), WBXJK(*), XOP(*), WORK(LWORK)
      REAL ONE
#else
      DOUBLE PRECISION TB0JK(*), L2TP(*), WBXJK(*), XOP(*), WORK(LWORK) 
      DOUBLE PRECISION ONE
#endif
C
      PARAMETER (ONE = 1.0D0)
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccsdinp.h>
C
      CALL QENTER('WBX_JK_ETA')

      KDA  = 1
      KEND1  = KDA + NMATAB(ISYMXOP)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWRK1
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in WBX_JK_ETA')
      END IF
C
      DO ISYMA = 1,NSYM
         ISYMD = MULD2H(ISYMA,ISYMXOP)
         DO A = 1,NVIR(ISYMA)
            KOFF1 = IFCVIR(ISYMD,ISYMA) + NORB(ISYMD)*(A - 1)
     *                                  + NRHF(ISYMD) + 1
            KOFF2 = KDA + IMATAB(ISYMD,ISYMA) + NVIR(ISYMD)*(A - 1)
            CALL DCOPY(NVIR(ISYMD),XOP(KOFF1),1,WORK(KOFF2),1)
         END DO
      END DO
C
      DO ISYMI = 1,NSYM
         ISYMBCD = MULD2H(ISTB0JK,ISYMI)
         DO I = 1,NRHF(ISYMI)
            DO ISYMD = 1,NSYM
               ISYMA = MULD2H(ISYMD,ISYMXOP)
               ISYMBCA = MULD2H(ISYMXOP,ISYMBCD)
               ISYMBC  = MULD2H(ISYMBCD,ISYMD)
               KOFF1   = 1
     *                + IMAABCI(ISYMBCD,ISYMI)
     *                + NMAABC(ISYMBCD)*(I-1)
     *                + IMAABC(ISYMBC,ISYMD)
               KOFF2   = KDA
     *                + IMATAB(ISYMD,ISYMA)
               KOFF3   = 1
     *                + IMAABCI(ISYMBCA,ISYMI)
     *                + NMAABC(ISYMBCA)*(I-1)
     *                + IMAABC(ISYMBC,ISYMA)
C
               NTOTBC = MAX(1,NMATAB(ISYMBC))
               NTOTD  = MAX(1,NVIR(ISYMD))
C
C WBX_JK(bcai) = WBX_JK(bcai)  - xop(da) tb0_jk(bcdi)  
C
c add_wbx(1v)
               CALL DGEMM('N','N',NMATAB(ISYMBC),NVIR(ISYMA),
     *                    NVIR(ISYMD),-ONE,TB0JK(KOFF1),NTOTBC,
     *                    WORK(KOFF2),NTOTD,
     *                    ONE,WBXJK(KOFF3),NTOTBC)
            END DO
         END DO
      END DO
C                             
C 2)             xop(il) tb0_jk(bcal) 
C                             
C
C SORT OCC-OCC  XOP ELEMENTS (i,l)
C
C
      KIL  = 1
      KEND1  = KIL + NMATIJ(ISYMXOP)
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWRK1
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in WBX_JK_ETA')
      END IF
C   
      DO ISYML = 1,NSYM
         ISYMI = MULD2H(ISYML,ISYMXOP)
         DO L = 1,NRHF(ISYML)
             KOFF1 = IFCRHF(ISYMI,ISYML) + NORB(ISYMI)*(L - 1) + 1
             KOFF2 = KIL + IMATIJ(ISYMI,ISYML) + NRHF(ISYMI)*(L - 1)
             CALL DCOPY(NRHF(ISYMI),XOP(KOFF1),1,WORK(KOFF2),1)
         END DO
      END DO
C
      DO ISYML = 1,NSYM
         ISYMBCA = MULD2H(ISTB0JK,ISYML)
         ISYMI   = MULD2H(ISYML,ISYMXOP)
C
         KOFF1   = 1 + IMAABCI(ISYMBCA,ISYML)
         KOFF2   = KIL + IMATIJ(ISYMI,ISYML)
         KOFF3   = 1 + IMAABCI(ISYMBCA,ISYMI)
C
         NTOTBCA = MAX(1,NMAABC(ISYMBCA))
         NTOTI  = MAX(1,NRHF(ISYMI))
C
C WBX_JK(bcai) = WBX_JK(bcai) - tb0_jk(bcal) xop(il) 
C
c add_wbx(2o)
         CALL DGEMM('N','T',NMAABC(ISYMBCA),NRHF(ISYMI),
     *               NRHF(ISYML),ONE,TB0JK(KOFF1),NTOTBCA,
     *               WORK(KOFF2),NTOTI,
     *               ONE,WBXJK(KOFF3),NTOTBCA)
      END DO
C
C 3)
C
C WBX_JK(bcai) = WBX_JK(bcai) 
C               + xop(j,b)*t2(ai,ck) - xop(k,b)*t2(ai,cj)
C               + xop(k,c)*t2(ai,bj) - xop(j,c)*t2(ai,bk)
C
      DO ISYMI = 1,NSYM
         ISYMBCA = MULD2H(ISYMWBX,ISYMI)
         DO ISYMA = 1,NSYM
            ISYMAI = MULD2H(ISYMA,ISYMI)
            ISYMBC = MULD2H(ISYMBCA,ISYMA)
            DO ISYMC = 1,NSYM
               ISYMCJ = MULD2H(ISYMC,ISYMJ)
               ISYMCK = MULD2H(ISYMC,ISYMK)
               ISYMB = MULD2H(ISYMBC,ISYMC)
               ISYMBJ = MULD2H(ISYMB,ISYMJ)
               ISYMBK = MULD2H(ISYMB,ISYMK)
               ISYMAICK = MULD2H(ISYMAI,ISYMCK)
               ISYMAICJ = MULD2H(ISYMAI,ISYMCJ)
               ISYMAIBJ = MULD2H(ISYMAI,ISYMBJ)
               ISYMAIBK = MULD2H(ISYMAI,ISYMBK)
               ISYMAIJ  = MULD2H(ISYMAI,ISYMJ)
               ISYMAIK  = MULD2H(ISYMAI,ISYMK)
               DO I = 1,NRHF(ISYMI)
                  DO A = 1,NVIR(ISYMA)
                     DO C = 1,NVIR(ISYMC)
                        DO B = 1,NVIR(ISYMB)
                           KOFF1 = IMAABCI(ISYMBCA,ISYMI)
     *                           + NMAABC(ISYMBCA)*(I-1)
     *                           + IMAABC(ISYMBC,ISYMA)
     *                           + NMATAB(ISYMBC)*(A-1)
     *                           + IMATAB(ISYMB,ISYMC)
     *                           + NVIR(ISYMB)*(C-1)
     *                           + B
C
C                       xop(j,b)*t2(ai,ck)
C
                        IF (ISYMBJ .EQ. ISYMXOP .AND. 
     *                          ISYMAICK .EQ. ISYML2) THEN
C
                            KOFF2 = IFCVIR(ISYMJ,ISYMB) + 
     *                               NORB(ISYMJ)*(B - 1) + J
                            KOFF3 = IT2SP(ISYMAIK,ISYMC)
     *                            + NCKI(ISYMAIK)*(C-1)
     *                            + ISAIK(ISYMAI,ISYMK)
     *                            + NT1AM(ISYMAI)*(K-1)
     *                            + IT1AM(ISYMA,ISYMI)
     *                            + NVIR(ISYMA)*(I-1)
     *                            + A
C
c add_wbx(3a)
                            WBXJK(KOFF1) = WBXJK(KOFF1) + 
     *                                      XOP(KOFF2)* L2TP(KOFF3) 
                        END IF 
C
C                      - xop(k,b)*t2(ai,cj) 
C
                        IF (ISYMBK .EQ. ISYMXOP .AND. 
     *                          ISYMAICJ .EQ. ISYML2) THEN
C
                            KOFF2 = IFCVIR(ISYMK,ISYMB) +  
     *                               NORB(ISYMK)*(B - 1) + K
                            KOFF3 = IT2SP(ISYMAIJ,ISYMC)
     *                            + NCKI(ISYMAIJ)*(C-1)
     *                            + ISAIK(ISYMAI,ISYMJ)
     *                            + NT1AM(ISYMAI)*(J-1)
     *                            + IT1AM(ISYMA,ISYMI)
     *                            + NVIR(ISYMA)*(I-1)
     *                            + A
C
c add_wbx(3b)
                            WBXJK(KOFF1) = WBXJK(KOFF1) -
     *                                      XOP(KOFF2)* L2TP(KOFF3)
                        END IF 
C
C                       xop(k,c)*t2(ai,bj) 
C
                        IF (ISYMCK .EQ. ISYMXOP .AND. 
     *                          ISYMAIBJ .EQ. ISYML2) THEN
C
                            KOFF2 = IFCVIR(ISYMK,ISYMC) +
     *                               NORB(ISYMK)*(C - 1) + K
                            KOFF3 = IT2SP(ISYMAIJ,ISYMB)
     *                            + NCKI(ISYMAIJ)*(B-1)
     *                            + ISAIK(ISYMAI,ISYMJ)
     *                            + NT1AM(ISYMAI)*(J-1)
     *                            + IT1AM(ISYMA,ISYMI)
     *                            + NVIR(ISYMA)*(I-1)
     *                            + A
C
c add_wbx(3c)
                            WBXJK(KOFF1) = WBXJK(KOFF1) +
     *                                      XOP(KOFF2)* L2TP(KOFF3)
                        END IF 
C
C                        - xop(j,c)*t2(ai,bk) 
C
                        IF (ISYMCJ .EQ. ISYMXOP .AND. 
     *                          ISYMAIBK .EQ. ISYML2) THEN
C
                            KOFF2 = IFCVIR(ISYMJ,ISYMC) +
     *                               NORB(ISYMJ)*(C - 1) + J
                            KOFF3 = IT2SP(ISYMAIK,ISYMB)
     *                            + NCKI(ISYMAIK)*(B-1)
     *                            + ISAIK(ISYMAI,ISYMK)
     *                            + NT1AM(ISYMAI)*(K-1)
     *                            + IT1AM(ISYMA,ISYMI)
     *                            + NVIR(ISYMA)*(I-1)
     *                            + A
C
c add_wbx(3d)
                            WBXJK(KOFF1) = WBXJK(KOFF1) -
     *                                      XOP(KOFF2)* L2TP(KOFF3)
                        END IF 
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
      CALL QEXIT('WBX_JK_ETA')
      RETURN
      END
C  /* Deck GET_T3B0_JK_L1F */
      SUBROUTINE GET_T3B0_JK_L1F(T3B0JK,ISYT3B0JK,
     *                            T1AM,ISYMT1,
     *                            XIAJB,ISINT1,
     *                            T2AM,ISYMT2,
     *                            FOCKCK,ISYMFCK,
     *                            ISYMJ,J,ISYMK,K)
********************************************************************
*
* In this routine we calculate these contributions to t3bar_0
* multipliers, which contain the t1 multipliers
*
* We thus calculate the following intermmediate (for two fixed 
* occupied index):
*
* T^JK(bcai) 
*       = P(ai,bj,ck) ( t(ai)*L(jb|kc) - t(ak)*L(jb|ic) )
*
* (1) 
*       = t(ai)*L(jb|kc) - t(ak)*L(jb|ic)
* (2) 
*       + t(bj)*L(ia|kc) - t(bk)*L(ia|jc)
* (3) 
*       + t(bj)*L(kc|ia) - t(bi)*L(kc|ja)
* (4) 
*       + t(ai)*L(kc|jb) - t(aj)*L(kc|ib)
* (5) 
*       + t(ck)*L(ia|jb) - t(cj)*L(ia|kb)
* (6) 
*       + t(ck)*L(jb|ia) - t(ci)*L(jb|ka)
********************************************************************
*
* OBS !
*       t in the following comments of this routine denotes 
*         Lagrange multipliers 
*
********************************************************************
*
* Filip Pawlowski, Aarhus, Winter 2003
*
********************************************************************
C
      IMPLICIT NONE
C
      INTEGER ISYT3B0JK,ISYMT1,ISINT1,ISYMT2,ISYMFCK,ISYMJ,ISYMK
      INTEGER ISYMI,ISYMBCA,ISYMA,ISYMBC,ISYMC,ISYMB,ISYMCK,ISYMBJ
      INTEGER ISYMAJ,ISYMAK,ISYMAI,ISYMBK,ISYMCI,ISYMBI,ISYMCJ
      INTEGER NBJ,NBK,NCK,NCJ,NAK,NAJ,NAI,NCI,NBI
      INTEGER NCKBJ,NCKAJ,NAKBJ,NCIBJ,NCKAI,NBJAI,NCJAI,NCKBI,NBKAI
      INTEGER ISYMAIJ,ISYMAIK,ISYMAJK,ISYMAKJ,ISYMCKI,ISYMCKJ,ISYMCIJ
      INTEGER KKA,KJA,KIA,KKC,KJC,KIC,KKB,KJB,KIB
      INTEGER KAJKC,KAIKC,KAIJC,KAIJB,KAIKB,KAKJB,KCKIB,KCKJB,KCIJB
      INTEGER KOFF1
      INTEGER INDEX
C
#if defined (SYS_CRAY)
      REAL T3B0JK(*),T1AM(*),XIAJB(*),T2AM(*),FOCKCK(*)
      REAL TWO
#else
      DOUBLE PRECISION T3B0JK(*),T1AM(*),XIAJB(*),T2AM(*),FOCKCK(*)
      DOUBLE PRECISION TWO
#endif
C
      PARAMETER (TWO = 2.0D0)
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('GET_T3B0_JK_L1F')
C
C---------------------------------------
C     Contract the integrals with T1.
C---------------------------------------
C
      DO ISYMI = 1, NSYM
         ISYMBCA = MULD2H(ISYT3B0JK,ISYMI)
         DO ISYMA = 1,NSYM
            ISYMBC = MULD2H(ISYMBCA,ISYMA) 
            ISYMAJ = MULD2H(ISYMA,ISYMJ)
            ISYMAK = MULD2H(ISYMA,ISYMK)
            ISYMAI = MULD2H(ISYMA,ISYMI)
            DO ISYMC = 1,NSYM
               ISYMB = MULD2H(ISYMBC,ISYMC)
               ISYMCK = MULD2H(ISYMC,ISYMK)
               ISYMBK = MULD2H(ISYMB,ISYMK)
               ISYMCI = MULD2H(ISYMC,ISYMI)
               ISYMBI = MULD2H(ISYMB,ISYMI)
               ISYMCJ = MULD2H(ISYMC,ISYMJ)
               ISYMBJ = MULD2H(ISYMB,ISYMJ)
C
               DO B = 1, NVIR(ISYMB)
                  NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
                  NBK = IT1AM(ISYMB,ISYMK) + NVIR(ISYMB)*(K - 1) + B
C
                  DO C = 1, NVIR(ISYMC)
C
                     NCK = IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K - 1) + C
                     NCJ = IT1AM(ISYMC,ISYMJ) + NVIR(ISYMC)*(J - 1) + C
C
                     NCKBJ = IT2AM(ISYMCK,ISYMBJ) + INDEX(NCK,NBJ)
C
                     DO A = 1, NVIR(ISYMA)
                        NAK = IT1AM(ISYMA,ISYMK) + NVIR(ISYMA)*(K-1) + A
                        NAJ = IT1AM(ISYMA,ISYMJ) + NVIR(ISYMA)*(J-1) + A
                        NCKAJ = IT2AM(ISYMCK,ISYMAJ) + INDEX(NCK,NAJ)
                        NAKBJ = IT2AM(ISYMAK,ISYMBJ) + INDEX(NAK,NBJ)
                        DO I = 1, NRHF(ISYMI)
C
                          NAI=IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
                          NCI=IT1AM(ISYMC,ISYMI) + NVIR(ISYMC)*(I-1) + C
                          NBI=IT1AM(ISYMB,ISYMI) + NVIR(ISYMB)*(I-1) + B
                          NCIBJ = IT2AM(ISYMCI,ISYMBJ) + INDEX(NCI,NBJ)
                          NCKAI = IT2AM(ISYMCK,ISYMAI) + INDEX(NCK,NAI)
                          NBJAI = IT2AM(ISYMBJ,ISYMAI) + INDEX(NBJ,NAI)
                          NCJAI = IT2AM(ISYMCJ,ISYMAI) + INDEX(NCJ,NAI)
                          NCKBI = IT2AM(ISYMCK,ISYMBI) + INDEX(NCK,NBI)
                          NBKAI = IT2AM(ISYMBK,ISYMAI) + INDEX(NBK,NAI)
C
                          KOFF1 = IMAABCI(ISYMBCA,ISYMI)
     *                          + NMAABC(ISYMBCA)*(I-1)
     *                          + IMAABC(ISYMBC,ISYMA)
     *                          + NMATAB(ISYMBC)*(A-1) 
     *                          + IMATAB(ISYMB,ISYMC)
     *                          + NVIR(ISYMB)*(C-1)
     *                          + B
C
                           IF (ISYMA .EQ. ISYMI) THEN 
                             T3B0JK(KOFF1) = T3B0JK(KOFF1) 
     *                                     + TWO*T1AM(NAI)*XIAJB(NCKBJ)
                           END IF
                           IF (ISYMB .EQ. ISYMJ) THEN 
                             T3B0JK(KOFF1) = T3B0JK(KOFF1) 
     *                                     + TWO*T1AM(NBJ)*XIAJB(NCKAI)
                           END IF
                           IF (ISYMC .EQ. ISYMK) THEN 
                             T3B0JK(KOFF1) = T3B0JK(KOFF1) 
     *                                     + TWO*T1AM(NCK)*XIAJB(NBJAI)
                           END IF
                           IF (ISYMA .EQ. ISYMK) THEN 
                             T3B0JK(KOFF1) = T3B0JK(KOFF1) 
     *                                     - T1AM(NAK)*XIAJB(NCIBJ)
                           END IF
                           IF (ISYMB .EQ. ISYMK) THEN 
                             T3B0JK(KOFF1) = T3B0JK(KOFF1) 
     *                                     - T1AM(NBK)*XIAJB(NCJAI)
                           END IF
                           IF (ISYMB .EQ. ISYMI) THEN 
                             T3B0JK(KOFF1) = T3B0JK(KOFF1) 
     *                                     - T1AM(NBI)*XIAJB(NCKAJ)
                           END IF
                           IF (ISYMA .EQ. ISYMJ) THEN 
                             T3B0JK(KOFF1) = T3B0JK(KOFF1) 
     *                                     - T1AM(NAJ)*XIAJB(NCKBI)
                           END IF
                           IF (ISYMC .EQ. ISYMJ) THEN 
                             T3B0JK(KOFF1) = T3B0JK(KOFF1) 
     *                                     - T1AM(NCJ)*XIAJB(NBKAI)
                           END IF
                           IF (ISYMC .EQ. ISYMI) THEN 
                             T3B0JK(KOFF1) = T3B0JK(KOFF1) 
     *                                     - T1AM(NCI)*XIAJB(NAKBJ)
                           END IF

C
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
C
C---------------------------------------
C     Contract the Fock matrix with T2.
C---------------------------------------
C
      DO ISYMI = 1, NSYM
         ISYMBCA = MULD2H(ISYT3B0JK,ISYMI)
         DO ISYMA = 1,NSYM
            ISYMBC = MULD2H(ISYMBCA,ISYMA) 
            ISYMAJ = MULD2H(ISYMA,ISYMJ)
            ISYMAK = MULD2H(ISYMA,ISYMK)
            ISYMAI = MULD2H(ISYMA,ISYMI)
            ISYMAIJ = MULD2H(ISYMAI,ISYMJ)
            ISYMAIK = MULD2H(ISYMAI,ISYMK)
            ISYMAJK = MULD2H(ISYMAJ,ISYMK)
            ISYMAKJ = MULD2H(ISYMAK,ISYMJ)
            DO ISYMC = 1,NSYM
               ISYMB = MULD2H(ISYMBC,ISYMC)
               ISYMCK = MULD2H(ISYMC,ISYMK)
               ISYMBK = MULD2H(ISYMB,ISYMK)
               ISYMCI = MULD2H(ISYMC,ISYMI)
               ISYMBI = MULD2H(ISYMB,ISYMI)
               ISYMCJ = MULD2H(ISYMC,ISYMJ)
               ISYMBJ = MULD2H(ISYMB,ISYMJ)
               ISYMCKI = MULD2H(ISYMCK,ISYMI)
               ISYMCKJ = MULD2H(ISYMCK,ISYMJ)
               ISYMCIJ = MULD2H(ISYMCI,ISYMJ)
C
               DO I = 1, NRHF(ISYMI)
C
                  DO A = 1, NVIR(ISYMA)
                     KKA = IT1AM(ISYMA,ISYMK)  + NVIR(ISYMA)*(K-1)+ A
                     KJA = IT1AM(ISYMA,ISYMJ) +  NVIR(ISYMA)*(J-1)+ A
                     KIA = IT1AM(ISYMA,ISYMI) +  NVIR(ISYMA)*(I-1)+ A
C
                     DO C = 1, NVIR(ISYMC)
                        KKC = IT1AM(ISYMC,ISYMK)  + NVIR(ISYMC)*(K-1)+ C
                        KJC = IT1AM(ISYMC,ISYMJ) +  NVIR(ISYMC)*(J-1)+ C
                        KIC = IT1AM(ISYMC,ISYMI) +  NVIR(ISYMC)*(I-1)+ C
C
                        KAJKC = IT2SP(ISYMAJK,ISYMC)
     *                        + NCKI(ISYMAJK)*(C-1)
     *                        + ISAIK(ISYMAJ,ISYMK)
     *                        + NT1AM(ISYMAJ)*(K-1)
     *                        + IT1AM(ISYMA,ISYMJ)
     *                        + NVIR(ISYMA)*(J-1)
     *                        + A
C
                        KAIKC = IT2SP(ISYMAIK,ISYMC)
     *                        + NCKI(ISYMAIK)*(C-1)
     *                        + ISAIK(ISYMAI,ISYMK)
     *                        + NT1AM(ISYMAI)*(K-1)
     *                        + IT1AM(ISYMA,ISYMI)
     *                        + NVIR(ISYMA)*(I-1)
     *                        + A
C
                        KAIJC = IT2SP(ISYMAIJ,ISYMC)
     *                        + NCKI(ISYMAIJ)*(C-1)
     *                        + ISAIK(ISYMAI,ISYMJ)
     *                        + NT1AM(ISYMAI)*(J-1)
     *                        + IT1AM(ISYMA,ISYMI)
     *                        + NVIR(ISYMA)*(I-1)
     *                        + A
C

                        DO B = 1, NVIR(ISYMB)
                           KKB=IT1AM(ISYMB,ISYMK)  + NVIR(ISYMB)*(K-1)+B
                           KJB=IT1AM(ISYMB,ISYMJ) +  NVIR(ISYMB)*(J-1)+B
                           KIB=IT1AM(ISYMB,ISYMI) +  NVIR(ISYMB)*(I-1)+B
C
                           KAIJB = IT2SP(ISYMAIJ,ISYMB)
     *                           + NCKI(ISYMAIJ)*(B-1)
     *                           + ISAIK(ISYMAI,ISYMJ)
     *                           + NT1AM(ISYMAI)*(J-1)
     *                           + IT1AM(ISYMA,ISYMI)
     *                           + NVIR(ISYMA)*(I-1)
     *                           + A
C
                           KAIKB = IT2SP(ISYMAIK,ISYMB)
     *                           + NCKI(ISYMAIK)*(B-1)
     *                           + ISAIK(ISYMAI,ISYMK)
     *                           + NT1AM(ISYMAI)*(K-1)
     *                           + IT1AM(ISYMA,ISYMI)
     *                           + NVIR(ISYMA)*(I-1)
     *                           + A
C
                           KAKJB = IT2SP(ISYMAKJ,ISYMB)
     *                           + NCKI(ISYMAKJ)*(B-1)
     *                           + ISAIK(ISYMAK,ISYMJ)
     *                           + NT1AM(ISYMAK)*(J-1)
     *                           + IT1AM(ISYMA,ISYMK)
     *                           + NVIR(ISYMA)*(K-1)
     *                           + A
C
                           KCKIB = IT2SP(ISYMCKI,ISYMB)
     *                           + NCKI(ISYMCKI)*(B-1)
     *                           + ISAIK(ISYMCK,ISYMI)
     *                           + NT1AM(ISYMCK)*(I-1)
     *                           + IT1AM(ISYMC,ISYMK)
     *                           + NVIR(ISYMC)*(K-1)
     *                           + C
C
                           KCKJB = IT2SP(ISYMCKJ,ISYMB)
     *                           + NCKI(ISYMCKJ)*(B-1)
     *                           + ISAIK(ISYMCK,ISYMJ)
     *                           + NT1AM(ISYMCK)*(J-1)
     *                           + IT1AM(ISYMC,ISYMK)
     *                           + NVIR(ISYMC)*(K-1)
     *                           + C
C

                           KCIJB = IT2SP(ISYMCIJ,ISYMB)
     *                           + NCKI(ISYMCIJ)*(B-1)
     *                           + ISAIK(ISYMCI,ISYMJ)
     *                           + NT1AM(ISYMCI)*(J-1)
     *                           + IT1AM(ISYMC,ISYMI)
     *                           + NVIR(ISYMC)*(I-1)
     *                           + C
C
                          KOFF1 = IMAABCI(ISYMBCA,ISYMI)
     *                          + NMAABC(ISYMBCA)*(I-1)
     *                          + IMAABC(ISYMBC,ISYMA)
     *                          + NMATAB(ISYMBC)*(A-1) 
     *                          + IMATAB(ISYMB,ISYMC)
     *                          + NVIR(ISYMB)*(C-1)
     *                          + B
C
                           IF (ISYMC .EQ. ISYMK) THEN
                             T3B0JK(KOFF1) = T3B0JK(KOFF1)
     *                                     + TWO*T2AM(KAIJB)*FOCKCK(KKC)
                           END IF
                           IF (ISYMA .EQ. ISYMI) THEN
                             T3B0JK(KOFF1) = T3B0JK(KOFF1)
     *                                     + TWO*T2AM(KCKJB)*FOCKCK(KIA)
                           END IF
                           IF (ISYMB .EQ. ISYMJ) THEN
                             T3B0JK(KOFF1) = T3B0JK(KOFF1)
     *                                     + TWO*T2AM(KAIKC)*FOCKCK(KJB)
                           END IF
                           IF (ISYMC .EQ. ISYMJ) THEN
                             T3B0JK(KOFF1) = T3B0JK(KOFF1)
     *                                     - T2AM(KAIKB)*FOCKCK(KJC)
                           END IF
                           IF (ISYMC .EQ. ISYMI) THEN
                             T3B0JK(KOFF1) = T3B0JK(KOFF1)
     *                                     - T2AM(KAKJB)*FOCKCK(KIC)
                           END IF
                           IF (ISYMA .EQ. ISYMK) THEN
                             T3B0JK(KOFF1) = T3B0JK(KOFF1)
     *                                     - T2AM(KCIJB)*FOCKCK(KKA)
                           END IF
                           IF (ISYMB .EQ. ISYMK) THEN
                             T3B0JK(KOFF1) = T3B0JK(KOFF1)
     *                                     - T2AM(KAIJC)*FOCKCK(KKB)
                           END IF
                           IF (ISYMB .EQ. ISYMI) THEN
                             T3B0JK(KOFF1) = T3B0JK(KOFF1)
     *                                     - T2AM(KAJKC)*FOCKCK(KIB)
                           END IF
                           IF (ISYMA .EQ. ISYMJ) THEN
                             T3B0JK(KOFF1) = T3B0JK(KOFF1)
     *                                     - T2AM(KCKIB)*FOCKCK(KJA)
                           END IF
C
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
C
C-----------
C     End.
C-----------
C
      CALL QEXIT('GET_T3B0_JK_L1F')
C
      RETURN
      END
C  /* Deck wbx_jk_fmat */
      SUBROUTINE WBX_JK_FMAT(T3B0JK,ISYT3B0JK,
     *                            T2AM,ISYMT2,
     *                            FOCKCK,ISYMFCK,
     *                            T3BOL2,T3BOG2,
     *                            XGADCK,XLADCK,ISYINT,
     *                            ISYMJ,J,ISYMK,K,
     *                            WORK,LWORK)
********************************************************************
*
* In this routine we calculate the following contributions to t3bar_X
* multipliers:
*
*                     <L2|[H,tau3]|HF>
*
* We use W^JK(bcai) intermmediate.
*
* We thus calculate :
*
*     W^JK(bcai) = W^JK(bcai) + T2TP(aijb)*F(kc) 
*                             - T2TP(aikb)*F(jc)
*                             + T2TP(aikc)*F(jb)
*                             - T2TP(aijc)*F(kb)
*
*                   + T2TP(aijd)*L(d^bkc)
*                   - T2TP(ajkd)*g(ibd^c)
*                   + T2TP(aikd)*L(d^cjb)
*                   - T2TP(akjd)*g(icd^b)
*
*                   + T2TP(ailb)*L(jl^kc)
*                   - T2TP(alkb)*g(il^jc)
*                   + T2TP(ailc)*L(kl^jb)
*                   - T2TP(aljc)*g(il^kb)
*
********************************************************************
*
* OBS !
*       T2TP in the following comments of this routine denotes 
*         Lagrange multipliers 
*
********************************************************************
*
* Filip Pawlowski, Aarhus, Winter 2003
*
********************************************************************
C
      IMPLICIT NONE
C
      INTEGER ISYT3B0JK,ISYMT2,ISYMFCK,ISYINT,ISYMJ,ISYMK,LWORK
      INTEGER ISYMI,ISYMBCA,ISYMA,ISYMBC,ISYMC,ISYMB,ISYMCK,ISYMBJ
      INTEGER ISYMAJ,ISYMAK,ISYMAI,ISYMBK,ISYMCI,ISYMBI,ISYMCJ
      INTEGER ISYMAIJ,ISYMAIK,ISYMAJK,ISYMAKJ,ISYMCKI,ISYMCKJ,ISYMCIJ
      INTEGER KKA,KJA,KIA,KKC,KJC,KIC,KKB,KJB,KIB
      INTEGER KAJKC,KAIKC,KAIJC,KAIJB,KAIKB,KAKJB,KCKIB,KCKJB,KCIJB
      INTEGER ISYKJ,ISYCN
      INTEGER KL2NBAI,KINTCN,KCBAI,KEND1,LWRK1
      INTEGER ISYMN,ISYBAI
      INTEGER NTOTC,NTOTN
      INTEGER ISYL2BAN,ISYINTNCI
      INTEGER KL2BAN,KINTNCI,KBACI
      INTEGER ISYBA,ISYCI
      INTEGER NTOTBA
      INTEGER ISYJK,ISYBN
      INTEGER KL2NCAI,KINTBN
      INTEGER ISYCAI
      INTEGER NTOTB
      INTEGER ISYL2CAN,ISYINTNBI
      INTEGER KL2CAN,KINTNBI,KCABI
      INTEGER ISYCA,ISYBI
      INTEGER NTOTCA
      INTEGER ISYMDAI,ISYMBCD,ISYMBCAI
      INTEGER KDAI,KBCD,KBCAI
      INTEGER ISYMDA,ISYMD
      INTEGER NTOTBC
      INTEGER ISYMJK,ISYMAD,ISYMACBI
      INTEGER KACBI
      INTEGER ISYMCBI
      INTEGER NTOTD,NTOTA
      INTEGER KAD,KABCI
      INTEGER ISYMBCI
      INTEGER KOFF1,KOFF2,KOFF3
      INTEGER ISYMABCI
C
#if defined (SYS_CRAY)
      REAL T3B0JK(*),T2AM(*),FOCKCK(*),T3BOL2(*),T3BOG2(*)
      REAL XGADCK(*),XLADCK(*)
      REAL WORK(LWORK)
      REAL ONE
#else
      DOUBLE PRECISION T3B0JK(*),T2AM(*),FOCKCK(*),T3BOL2(*),T3BOG2(*)
      DOUBLE PRECISION XGADCK(*),XLADCK(*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION ONE
#endif
C
      PARAMETER (ONE = 1.0D0)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('WBX_JK_FMAT')
C
C---------------------------------------------------
C     Calculate all (four) Fock matrix contributions
C---------------------------------------------------
C
      DO ISYMI = 1, NSYM
         ISYMBCA = MULD2H(ISYT3B0JK,ISYMI)
         DO ISYMA = 1,NSYM
            ISYMBC = MULD2H(ISYMBCA,ISYMA) 
            ISYMAJ = MULD2H(ISYMA,ISYMJ)
            ISYMAK = MULD2H(ISYMA,ISYMK)
            ISYMAI = MULD2H(ISYMA,ISYMI)
            ISYMAIJ = MULD2H(ISYMAI,ISYMJ)
            ISYMAIK = MULD2H(ISYMAI,ISYMK)
            ISYMAJK = MULD2H(ISYMAJ,ISYMK)
            ISYMAKJ = MULD2H(ISYMAK,ISYMJ)
            DO ISYMC = 1,NSYM
               ISYMB = MULD2H(ISYMBC,ISYMC)
               ISYMCK = MULD2H(ISYMC,ISYMK)
               ISYMBK = MULD2H(ISYMB,ISYMK)
               ISYMCI = MULD2H(ISYMC,ISYMI)
               ISYMBI = MULD2H(ISYMB,ISYMI)
               ISYMCJ = MULD2H(ISYMC,ISYMJ)
               ISYMBJ = MULD2H(ISYMB,ISYMJ)
               ISYMCKI = MULD2H(ISYMCK,ISYMI)
               ISYMCKJ = MULD2H(ISYMCK,ISYMJ)
               ISYMCIJ = MULD2H(ISYMCI,ISYMJ)
C
               DO I = 1, NRHF(ISYMI)
C
                  DO A = 1, NVIR(ISYMA)
C
                     DO C = 1, NVIR(ISYMC)
                        KKC = IT1AM(ISYMC,ISYMK)  + NVIR(ISYMC)*(K-1)+ C
                        KJC = IT1AM(ISYMC,ISYMJ) +  NVIR(ISYMC)*(J-1)+ C
C
                        KAIKC = IT2SP(ISYMAIK,ISYMC)
     *                        + NCKI(ISYMAIK)*(C-1)
     *                        + ISAIK(ISYMAI,ISYMK)
     *                        + NT1AM(ISYMAI)*(K-1)
     *                        + IT1AM(ISYMA,ISYMI)
     *                        + NVIR(ISYMA)*(I-1)
     *                        + A
C
                        KAIJC = IT2SP(ISYMAIJ,ISYMC)
     *                        + NCKI(ISYMAIJ)*(C-1)
     *                        + ISAIK(ISYMAI,ISYMJ)
     *                        + NT1AM(ISYMAI)*(J-1)
     *                        + IT1AM(ISYMA,ISYMI)
     *                        + NVIR(ISYMA)*(I-1)
     *                        + A
C

                        DO B = 1, NVIR(ISYMB)
                           KKB=IT1AM(ISYMB,ISYMK)  + NVIR(ISYMB)*(K-1)+B
                           KJB=IT1AM(ISYMB,ISYMJ) +  NVIR(ISYMB)*(J-1)+B
C
                           KAIJB = IT2SP(ISYMAIJ,ISYMB)
     *                           + NCKI(ISYMAIJ)*(B-1)
     *                           + ISAIK(ISYMAI,ISYMJ)
     *                           + NT1AM(ISYMAI)*(J-1)
     *                           + IT1AM(ISYMA,ISYMI)
     *                           + NVIR(ISYMA)*(I-1)
     *                           + A
C
                           KAIKB = IT2SP(ISYMAIK,ISYMB)
     *                           + NCKI(ISYMAIK)*(B-1)
     *                           + ISAIK(ISYMAI,ISYMK)
     *                           + NT1AM(ISYMAI)*(K-1)
     *                           + IT1AM(ISYMA,ISYMI)
     *                           + NVIR(ISYMA)*(I-1)
     *                           + A
C
                          KOFF1 = IMAABCI(ISYMBCA,ISYMI)
     *                          + NMAABC(ISYMBCA)*(I-1)
     *                          + IMAABC(ISYMBC,ISYMA)
     *                          + NMATAB(ISYMBC)*(A-1) 
     *                          + IMATAB(ISYMB,ISYMC)
     *                          + NVIR(ISYMB)*(C-1)
     *                          + B
C
                           IF (ISYMC .EQ. ISYMK) THEN
                             T3B0JK(KOFF1) = T3B0JK(KOFF1)
     *                                     + T2AM(KAIJB)*FOCKCK(KKC)
                           END IF
                           IF (ISYMB .EQ. ISYMJ) THEN
                             T3B0JK(KOFF1) = T3B0JK(KOFF1)
     *                                     + T2AM(KAIKC)*FOCKCK(KJB)
                           END IF
                           IF (ISYMC .EQ. ISYMJ) THEN
                             T3B0JK(KOFF1) = T3B0JK(KOFF1)
     *                                     - T2AM(KAIKB)*FOCKCK(KJC)
                           END IF
                           IF (ISYMB .EQ. ISYMK) THEN
                             T3B0JK(KOFF1) = T3B0JK(KOFF1)
     *                                     - T2AM(KAIJC)*FOCKCK(KKB)
                           END IF
C
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
C
C---------------------------------------------------
C     Calculate all (four) occupied contributions
C---------------------------------------------------
C
C
C=================================================
C     Calculate (1L)   - t^(ab)_(in) L(jn|kc)
C                    = - t^(ba)_(ni) L(kc|jn)
C                             
C                      - T(nbai) I^KJ(cn)
C=================================================
C
C-------------------------------
C     Sort T2AM(bnia) as T(nbai)
C-------------------------------
C
      ISYKJ = MULD2H(ISYMK,ISYMJ)
      ISYCN = MULD2H(ISYINT,ISYKJ)
C
      KL2NBAI = 1
      KINTCN  = KL2NBAI + NT2SQ(ISYMT2)
      KCBAI   = KINTCN  + NT1AM(ISYCN)
      KEND1   = KCBAI   + NMAAOBCI(ISYT3B0JK)
      LWRK1   = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in WBX_JK_FMAT (1L)')
      END IF
C
      CALL DZERO(WORK(KCBAI),NMAAOBCI(ISYT3B0JK))
C
      CALL SORT_T2_IABJ(WORK(KL2NBAI),T2AM,ISYMT2)
C
C------------------------------------------------
C     Sort L(kc|jn) = T3BOL2(c,j,k,n) as I^KJ(cn)
C------------------------------------------------
C
      CALL SORT_INT_AJ_IK(WORK(KINTCN),T3BOL2,ISYINT,ISYMK,K,ISYMJ,J)
C
C------------------------------------------
C    Multiply I^KJ(cn) T(nbai) = T^JK(cbai)
C------------------------------------------
C
      DO ISYMN = 1, NSYM
            ISYMC = MULD2H(ISYCN,ISYMN)
            ISYBAI = MULD2H(ISYMT2,ISYMN)
C      
            KOFF1 = KINTCN
     *            + IT1AM(ISYMC,ISYMN)
            KOFF2 = KL2NBAI
     *            + IMAJBAI(ISYMN,ISYBAI)
            KOFF3 = KCBAI
     *            + IMAAOBCI(ISYMC,ISYBAI)
C
            NTOTC = MAX(NVIR(ISYMC),1)
            NTOTN = MAX(NRHF(ISYMN),1)
C
            CALL DGEMM('N','N',NVIR(ISYMC),NMAABI(ISYBAI),NRHF(ISYMN),
     *                 -ONE,WORK(KOFF1),NTOTC,WORK(KOFF2),NTOTN,
     *                 ONE,WORK(KOFF3),NTOTC)
C
      END DO ! ISYMN
C
C     T3B0JK(bcai) = T3B0JK(bcai) + T^JK(cbai)
C
C  w3x_occ(1L)
C
      CALL FBACI(T3B0JK,WORK(KCBAI),ISYT3B0JK)
C
C=================================================
C     Calculate (1g)   t^(ab)_(nk) g(in|jc)
C                    = t^(ba)_(kn) g(jc|in)
C                    
C                      T^K(ban) I^J(nci)
C=================================================
C
C-------------------------------
C     Sort T2AM(bkna) as T^K(ban)
C-------------------------------
C
      ISYL2BAN = MULD2H(ISYMT2,ISYMK)
      ISYINTNCI = MULD2H(ISYINT,ISYMJ)
C
      KL2BAN = KEND1
      KINTNCI  = KL2BAN + NMAABI(ISYL2BAN)
      KBACI   = KINTNCI + NCKI(ISYINTNCI)
      KEND1 = KBACI + NMAAB_CI(ISYT3B0JK)
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in WBX_JK_FMAT (1g)')
      END IF
C
      CALL DZERO(WORK(KBACI),NMAAB_CI(ISYT3B0JK))
C
      CALL SORT_T2_ABJ(WORK(KL2BAN),ISYMK,K,T2AM,ISYMT2)
C
C----------------------------------------------
C     Sort g(jc|in) = T3OG2(c,i,j,n) as I^J(nci)
C----------------------------------------------
C
      CALL SORT_INT_JAK_I(WORK(KINTNCI),T3BOG2,ISYINT,ISYMJ,J)
C
C------------------------------------------
C    Multiply T^K(ban) * I^J(nci) = T^JK(baci)
C------------------------------------------
C
      DO ISYMN = 1, NSYM
         ISYBA = MULD2H(ISYL2BAN,ISYMN)
         ISYCI = MULD2H(ISYINTNCI,ISYMN)
C
         KOFF1 = KL2BAN
     *         + IMAABI(ISYBA,ISYMN)
         KOFF2 = KINTNCI
     *         + IMAIAJ(ISYMN,ISYCI)
         KOFF3 = KBACI
     *         + IMAAB_CI(ISYBA,ISYCI)
C
         NTOTBA = MAX(NMATAB(ISYBA),1)
         NTOTN  = MAX(NRHF(ISYMN),1)
C
         CALL DGEMM('N','N',NMATAB(ISYBA),NT1AM(ISYCI),NRHF(ISYMN),
     *              ONE,WORK(KOFF1),NTOTBA,WORK(KOFF2),NTOTN,
     *              ONE,WORK(KOFF3),NTOTBA)
C
      END DO ! ISYMN
C
C-------------------------------------------
C     T3B0JK(bcai) = T3B0JK(bcai) + T^JK(baci)
C-------------------------------------------
C
C  w3x_occ(1g)
C
      CALL FACBI(T3B0JK,WORK(KBACI),ISYT3B0JK)
C
C=================================================
C     Calculate (2L)   - t^(ac)_(in) L(kn|jb)
C                    = - t^(ca)_(ni) L(jb|kn)
C                             
C                      - T(ncai) I^JK(bn)
C=================================================
C
C-------------------------------
C     Sort T2AM(cnia) as T(ncai)
C-------------------------------
C
      ISYJK = MULD2H(ISYMJ,ISYMK)
      ISYBN = MULD2H(ISYINT,ISYJK)
C
      KL2NCAI = KEND1
      KINTBN  = KL2NCAI + NT2SQ(ISYMT2)
      KEND1   = KINTBN  + NT1AM(ISYBN)
      LWRK1   = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in WBX_JK_FMAT (2L)')
      END IF
C
      CALL SORT_T2_IABJ(WORK(KL2NCAI),T2AM,ISYMT2)
C
C------------------------------------------------
C     Sort L(jb|kn) = T3BOL2(b,k,j,n) as I^JK(bn)
C------------------------------------------------
C
      CALL SORT_INT_AJ_IK(WORK(KINTBN),T3BOL2,ISYINT,ISYMJ,J,ISYMK,K)
C
C------------------------------------------
C    Multiply I^JK(bn) T(ncai) = T^JK(bcai)
C------------------------------------------
C
      DO ISYMN = 1, NSYM
            ISYMB = MULD2H(ISYBN,ISYMN)
            ISYCAI = MULD2H(ISYMT2,ISYMN)
C      
            KOFF1 = KINTBN
     *            + IT1AM(ISYMB,ISYMN)
            KOFF2 = KL2NCAI
     *            + IMAJBAI(ISYMN,ISYCAI)
            KOFF3 = IMAAOBCI(ISYMB,ISYCAI) + 1
C
            NTOTB = MAX(NVIR(ISYMB),1)
            NTOTN = MAX(NRHF(ISYMN),1)
C
C  w3x_occ(2L)
C
            CALL DGEMM('N','N',NVIR(ISYMB),NMAABI(ISYCAI),NRHF(ISYMN),
     *                 -ONE,WORK(KOFF1),NTOTB,WORK(KOFF2),NTOTN,
     *                 ONE,T3B0JK(KOFF3),NTOTB)
C
      END DO ! ISYMN
C
C=================================================
C     Calculate (2g)   t^(ac)_(nj) g(in|kb)
C                    = t^(ca)_(jn) g(kb|in)
C                    
C                      T^J(can) I^K(nbi)
C=================================================
C
C-------------------------------
C     Sort T2AM(cjna) as T^J(can)
C-------------------------------
C
      ISYL2CAN = MULD2H(ISYMT2,ISYMJ)
      ISYINTNBI = MULD2H(ISYINT,ISYMK)
C
      KL2CAN = KEND1
      KINTNBI  = KL2CAN + NMAABI(ISYL2CAN)
      KCABI   = KINTNBI + NCKI(ISYINTNBI)
      KEND1 = KCABI + NMAAB_CI(ISYT3B0JK)
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in WBX_JK_FMAT (2g)')
      END IF
C
      CALL DZERO(WORK(KCABI),NMAAB_CI(ISYT3B0JK))
C
      CALL SORT_T2_ABJ(WORK(KL2CAN),ISYMJ,J,T2AM,ISYMT2)
C
C----------------------------------------------
C     Sort g(kb|in) = T3OG2(b,i,k,n) as I^K(nbi)
C----------------------------------------------
C
      CALL SORT_INT_JAK_I(WORK(KINTNBI),T3BOG2,ISYINT,ISYMK,K)
C
C---------------------------------------------
C    Multiply T^J(can) * I^K(nbi) = T^JK(cabi)
C---------------------------------------------
C
      DO ISYMN = 1, NSYM
         ISYCA = MULD2H(ISYL2CAN,ISYMN)
         ISYBI = MULD2H(ISYINTNBI,ISYMN)
C
         KOFF1 = KL2CAN
     *         + IMAABI(ISYCA,ISYMN)
         KOFF2 = KINTNBI
     *         + IMAIAJ(ISYMN,ISYBI)
         KOFF3 = KCABI
     *         + IMAAB_CI(ISYCA,ISYBI)
C
         NTOTCA = MAX(NMATAB(ISYCA),1)
         NTOTN  = MAX(NRHF(ISYMN),1)
C
         CALL DGEMM('N','N',NMATAB(ISYCA),NT1AM(ISYBI),NRHF(ISYMN),
     *              ONE,WORK(KOFF1),NTOTCA,WORK(KOFF2),NTOTN,
     *              ONE,WORK(KOFF3),NTOTCA)
C
      END DO ! ISYMN
C
C-------------------------------------------
C     T3B0JK(bcai) = T3B0JK(bcai) + T^JK(cabi)
C-------------------------------------------
C
C  w3x_occ(2g)
C
      CALL FBCAI(T3B0JK,WORK(KCABI),ISYT3B0JK)
C
C---------------------------------------------------
C     Calculate all (four) virtual contributions
C---------------------------------------------------
C
C
C***************************************************
C 1)  l^ad_ij L(dbkc) - l^ad_jk g(ibdc) 
C***************************************************
C
C T2AM(djia) =   I^J(dai) 
C
C L(dbkc) = I(dbck)  sorted as I^K(bcd)
C
C T^JK(bcai) = T^JK(bcai) + I^K(bcd)*I^J(dai)
C
C symmetry and work allocation
C

      ISYMDAI = MULD2H(ISYMT2,ISYMJ)
      ISYMBCD = MULD2H(ISYINT,ISYMK)
      ISYMBCAI = MULD2H(ISYMBCD,ISYMDAI)
      KDAI  = KEND1
      KBCD  = KDAI  + NMAABI(ISYMDAI)
      KBCAI = KBCD  + NMAABC(ISYMBCD)
      KEND1 = KBCAI + NMAABCI(ISYMBCAI)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in WBX_JK_FMAT (1Lv)')
      END IF
C
      CALL DZERO(WORK(KBCAI),NMAABCI(ISYMBCAI))
C
C
C  sort l^ad_ij = T2AM(djia) as I^J(dai) 
C
      CALL SORT_T2_ABJ(WORK(KDAI),ISYMJ,J,T2AM,ISYMT2)
C
C L(dbkc) = I(dbck)  sorted as I^K(bcd)
C 
      CALL SORT_INT_BCA(WORK(KBCD),ISYMK,K,XLADCK,ISYINT)
C
C T^JK(bcai) = T^JK(bcai) + I^K(bcd)*I^J(dai) 
C
      DO ISYMI = 1,NSYM
         ISYMDA = MULD2H(ISYMDAI,ISYMI)
         DO I =  1,NRHF(ISYMI)
            DO ISYMD = 1,NSYM
               ISYMBC  = MULD2H(ISYMBCD,ISYMD)
               ISYMA   = MULD2H(ISYMDA,ISYMD)
               ISYMBCA = MULD2H(ISYMBC,ISYMA)
               KOFF1   = KBCD + IMAABC(ISYMBC,ISYMD) 
               KOFF2   = KDAI 
     *                + IMAABI(ISYMDA,ISYMI)  
     *                + NMATAB(ISYMDA)*(I-1)
     *                + IMATAB(ISYMD,ISYMA)
               KOFF3   = KBCAI + IMAABCI(ISYMBCA,ISYMI)
     *                + NMAABC(ISYMBCA)*(I-1)
     *                + IMAABC(ISYMBC,ISYMA)
C  
               NTOTBC = MAX(1,NMATAB(ISYMBC))
               NTOTD  = MAX(1,NVIR(ISYMD))
C
               CALL DGEMM('N','N',NMATAB(ISYMBC),NVIR(ISYMA),
     *                    NVIR(ISYMD),ONE,WORK(KOFF1),NTOTBC,
     *                    WORK(KOFF2),NTOTD, 
     *                    ONE,WORK(KOFF3),NTOTBC) 
            END DO
         END DO
      END DO
C
C  w3x_vir(1L)
C
      CALL FAB_CI(T3B0JK,WORK(KBCAI),ISYT3B0JK)
C

C
C  - l^ad_jk g(ibdc) 
C
C  T2AM(aJKd) = I^JK(ad)
C 
C g(ibdc) = I(dcbi) 
C
C T^JK(bcai) = T^JK(bcai) - I^JK(ad) * I(dcbi)
C
C
C symmetry and work allocation
C
      ISYMJK = MULD2H(ISYMJ,ISYMK)
      ISYMAD = MULD2H(ISYMT2,ISYMJK)
      ISYMACBI = MULD2H(ISYINT,ISYMAD)
      KAD    = KEND1
      KACBI  = KAD + NMATAB(ISYMAD)
      KEND1  = KACBI + NMAABCI(ISYMACBI) 
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in WBX_JK_FMAT (1gv)')
      END IF
C
      CALL DZERO(WORK(KACBI),NMAABCI(ISYMACBI))
C
C  T2AM(aJKd) = I^JK(ad)
C
      CALL SORT_T2_AB(WORK(KAD),ISYMJ,J,ISYMK,K,T2AM,ISYMT2)
C
C T^JK(bcai) = T^JK(bcai) +  I^JK(ad) * I(dcbi)
C
      DO ISYMD = 1,NSYM
         ISYMA   = MULD2H(ISYMAD,ISYMD)
         ISYMCBI = MULD2H(ISYMD,ISYINT)
         KOFF1   = KAD + IMATAB(ISYMA,ISYMD)
         KOFF2   = 1  + IMAAOBCI(ISYMD,ISYMCBI)
         KOFF3   = KACBI  + IMAAOBCI(ISYMA,ISYMCBI)
C
         NTOTD  = MAX(1,NVIR(ISYMD))
         NTOTA  = MAX(1,NVIR(ISYMA))
C
         CALL DGEMM('N','N',NVIR(ISYMA),NMAABI(ISYMCBI),
     *                    NVIR(ISYMD),-ONE,WORK(KOFF1),NTOTA,
     *                    XGADCK(KOFF2),NTOTD,
     *                    ONE,WORK(KOFF3),NTOTD)
      END DO
C
C T^JK(bcai) = T^JK(bcai) +  I^JK(ad) * I(dcbi)
C
C  w3x_vir(1g)
C
      CALL FCBAI(T3B0JK,WORK(KACBI),ISYT3B0JK)
C
C****************************************************
C 2)  l^ad_ik L(dcjb) - l^ad_kj g(icdb)
C****************************************************
C
C T2AM(dkia) =   I^K(dai)
C
C L(dcjb) = I(dcbj) stored as I^J(bcd)
C
C T^JK(bcai) = T^JK(bcai) + I^J(bcd)*I^K(dai)
C
C symmetry and work allocation
C
      ISYMDAI = MULD2H(ISYMT2,ISYMK)
      ISYMBCD = MULD2H(ISYINT,ISYMJ)
      ISYMBCAI = MULD2H(ISYMBCD,ISYMDAI)
      KDAI  = KEND1
      KBCD  = KDAI  + NMAABI(ISYMDAI)
      KBCAI = KBCD  + NMAABC(ISYMBCD)
      KEND1 = KBCAI + NMAABCI(ISYMBCAI)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in WBX_JK_FMAT (2Lv)')
      END IF
C
      CALL DZERO(WORK(KBCAI),NMAABCI(ISYMBCAI))
C
C T2AM(dkia) =   I^K(dai)
C
      CALL SORT_T2_ABJ(WORK(KDAI),ISYMK,K,T2AM,ISYMT2)
C
C L(dcjb) = I(dcbj) sorted as I^J(bcd)
C
      CALL SORT_INT_CBA(WORK(KBCD),ISYMJ,J,XLADCK,ISYINT)
C
C T^JK(bcai) = T^JK(bcai) + I^J(bcd)*I^K(dai)
C
      DO ISYMI = 1,NSYM
         ISYMDA = MULD2H(ISYMDAI,ISYMI)
         DO I =  1,NRHF(ISYMI)
            DO ISYMD = 1,NSYM
               ISYMBC  = MULD2H(ISYMBCD,ISYMD)
               ISYMA   = MULD2H(ISYMDA,ISYMD)
               ISYMBCA = MULD2H(ISYMBC,ISYMA)
               KOFF1   = KBCD + IMAABC(ISYMBC,ISYMD)
               KOFF2   = KDAI
     *                + IMAABI(ISYMDA,ISYMI)
     *                + NMATAB(ISYMDA)*(I-1)
     *                + IMATAB(ISYMD,ISYMA)
               KOFF3   = KBCAI + IMAABCI(ISYMBCA,ISYMI)
     *                + NMAABC(ISYMBCA)*(I-1)
     *                + IMAABC(ISYMBC,ISYMA)
C 
               NTOTBC = MAX(1,NMATAB(ISYMBC))
               NTOTD  = MAX(1,NVIR(ISYMD))
C
               CALL DGEMM('N','N',NMATAB(ISYMBC),NVIR(ISYMA),
     *                    NVIR(ISYMD),ONE,WORK(KOFF1),NTOTBC,
     *                    WORK(KOFF2),NTOTD,
     *                    ONE,WORK(KOFF3),NTOTBC)
            END DO
         END DO
      END DO
C
C  w3x_vir(2L)
C
      CALL FAB_CI(T3B0JK,WORK(KBCAI),ISYT3B0JK)
C
C  - l^ad_kj g(icdb)
C
C  T2AM(akjd) = I^KJ(ad)
C
C g(icdb) = I(dbci)
C
C T^JK(bcai) = T^JK(bcai) + I^KJ(ad) * I(dbci)
C
C symmetry and work allocation
C
      ISYMJK = MULD2H(ISYMJ,ISYMK)
      ISYMAD = MULD2H(ISYMT2,ISYMJK)
      ISYMABCI = MULD2H(ISYINT,ISYMAD)
      KAD    = KEND1
      KABCI  = KAD + NMATAB(ISYMAD)
      KEND1  = KABCI + NMAABCI(ISYMABCI)
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in WBX_JK_FMAT (2gv)')
      END IF
C
      CALL DZERO(WORK(KABCI),NMAABCI(ISYMABCI))
C
C  T2AM(akjd) = I^KJ(ad)
C
      CALL SORT_T2_AB(WORK(KAD),ISYMK,K,ISYMJ,J,T2AM,ISYMT2)
C
C T^JK(bcai) = T^JK(bcai) + I^KJ(ad) * I(dbci)
C
      DO ISYMD = 1,NSYM
         ISYMA   = MULD2H(ISYMAD,ISYMD)
         ISYMBCI = MULD2H(ISYMD,ISYINT)
         KOFF1   = KAD + IMATAB(ISYMA,ISYMD)
         KOFF2   = 1  + IMAAOBCI(ISYMD,ISYMBCI)
         KOFF3   = KABCI  + IMAAOBCI(ISYMA,ISYMBCI)
C
         NTOTD  = MAX(1,NVIR(ISYMD))
         NTOTA  = MAX(1,NVIR(ISYMA))
C
C  work(abci) = I^KJ(ad) * I(dbci)
C
         CALL DGEMM('N','N',NVIR(ISYMA),NMAABI(ISYMBCI),
     *                    NVIR(ISYMD),-ONE,WORK(KOFF1),NTOTA,
     *                    XGADCK(KOFF2),NTOTD,
     *                    ONE,WORK(KOFF3),NTOTA)
      END DO
C
C T^JK(bcai) = T^JK(bcai) + work(abci) 
C
C  w3x_vir(2g)
C
      CALL FCABI(T3B0JK,WORK(KABCI),ISYT3B0JK)
C
C
C-----------
C     End.
C-----------
C
      CALL QEXIT('WBX_JK_FMAT')
C
      RETURN
      END
C  /* Deck WBX_JK_L1 */
      SUBROUTINE WBX_JK_L1(T3B0JK,ISYT3B0JK,
     *                            T1AM,ISYMT1,
     *                            XIAJB,ISINT1,
     *                            ISYMJ,J,ISYMK,K)
********************************************************************
*
* In this routine we calculate the following contributions to t3bar_X
* multipliers:
*
*                     <L1Y|[H^,tau3]|HF>
*
* We use W^JK(bcai) intermmediate.
*
* We thus calculate :
*
*    WMAT^BC(aikj) = WMAT^BC(aikj) + T1(ai)*L(jBkC)             
*                                  - T1(ak)*L(jBiC)
*                                  + T1(ai)*L(kCjB)
*                                  - T1(aj)*L(kCiB)
*
********************************************************************
*
* OBS !
*       T1AM in the following comments of this routine denotes 
*         Lagrange multipliers 
*
********************************************************************
*
* Filip Pawlowski, Aarhus, Winter 2003
*
********************************************************************
C
      IMPLICIT NONE
C
      INTEGER ISYT3B0JK,ISYMT1,ISINT1,ISYMJ,ISYMK
      INTEGER ISYMI,ISYMBCA,ISYMA,ISYMBC,ISYMC,ISYMB,ISYMCK,ISYMBJ
      INTEGER ISYMAJ,ISYMAK,ISYMAI,ISYMBK,ISYMCI,ISYMBI,ISYMCJ
      INTEGER NBJ,NBK,NCK,NCJ,NAK,NAJ,NAI,NCI,NBI
      INTEGER NCKBJ,NCKAJ,NAKBJ,NCIBJ,NCKAI,NBJAI,NCJAI,NCKBI,NBKAI
      INTEGER KOFF1
      INTEGER INDEX
C
#if defined (SYS_CRAY)
      REAL T3B0JK(*),T1AM(*),XIAJB(*)
      REAL TWO
#else
      DOUBLE PRECISION T3B0JK(*),T1AM(*),XIAJB(*)
      DOUBLE PRECISION TWO
#endif
C
      PARAMETER (TWO = 2.0D0)
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('WBX_JK_L1')
C
C---------------------------------------
C     Contract the integrals with T1.
C---------------------------------------
C
      DO ISYMI = 1, NSYM
         ISYMBCA = MULD2H(ISYT3B0JK,ISYMI)
         DO ISYMA = 1,NSYM
            ISYMBC = MULD2H(ISYMBCA,ISYMA) 
            ISYMAJ = MULD2H(ISYMA,ISYMJ)
            ISYMAK = MULD2H(ISYMA,ISYMK)
            ISYMAI = MULD2H(ISYMA,ISYMI)
            DO ISYMC = 1,NSYM
               ISYMB = MULD2H(ISYMBC,ISYMC)
               ISYMCK = MULD2H(ISYMC,ISYMK)
               ISYMBK = MULD2H(ISYMB,ISYMK)
               ISYMCI = MULD2H(ISYMC,ISYMI)
               ISYMBI = MULD2H(ISYMB,ISYMI)
               ISYMCJ = MULD2H(ISYMC,ISYMJ)
               ISYMBJ = MULD2H(ISYMB,ISYMJ)
C
               DO B = 1, NVIR(ISYMB)
                  NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
C
                  DO C = 1, NVIR(ISYMC)
C
                     NCK = IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K - 1) + C
                     NCKBJ = IT2AM(ISYMCK,ISYMBJ) + INDEX(NCK,NBJ)
C
                     DO A = 1, NVIR(ISYMA)
                        NAJ = IT1AM(ISYMA,ISYMJ) + NVIR(ISYMA)*(J-1) + A
                        NAK = IT1AM(ISYMA,ISYMK) + NVIR(ISYMA)*(K-1) + A
C
                        DO I = 1, NRHF(ISYMI)
C
                          NAI=IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
                          NCI=IT1AM(ISYMC,ISYMI) + NVIR(ISYMC)*(I-1) + C
                          NBI=IT1AM(ISYMB,ISYMI) + NVIR(ISYMB)*(I-1) + B
                          NCIBJ = IT2AM(ISYMCI,ISYMBJ) + INDEX(NCI,NBJ)
                          NCKBI = IT2AM(ISYMCK,ISYMBI) + INDEX(NCK,NBI)
C
                          KOFF1 = IMAABCI(ISYMBCA,ISYMI)
     *                          + NMAABC(ISYMBCA)*(I-1)
     *                          + IMAABC(ISYMBC,ISYMA)
     *                          + NMATAB(ISYMBC)*(A-1) 
     *                          + IMATAB(ISYMB,ISYMC)
     *                          + NVIR(ISYMB)*(C-1)
     *                          + B
C
                           IF (ISYMA .EQ. ISYMI) THEN 
                             T3B0JK(KOFF1) = T3B0JK(KOFF1) 
     *                                     + TWO*T1AM(NAI)*XIAJB(NCKBJ)
                           END IF
                           IF (ISYMA .EQ. ISYMK) THEN 
                             T3B0JK(KOFF1) = T3B0JK(KOFF1) 
     *                                     - T1AM(NAK)*XIAJB(NCIBJ)
                           END IF
                           IF (ISYMA .EQ. ISYMJ) THEN 
                             T3B0JK(KOFF1) = T3B0JK(KOFF1) 
     *                                     - T1AM(NAJ)*XIAJB(NCKBI)
                           END IF
C
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
C
C-----------
C     End.
C-----------
C
      CALL QEXIT('WBX_JK_L1')
C
      RETURN
      END
C  /* Deck tetax_jk_bc */
      SUBROUTINE TETAX_JK_BC(T0JK,IST0JK,XOP,ISYMXOP,TETAXJK,ISTETAXJK,
     *                      WORK,LWORK)
C
C TETAXJK(bcai) = TETAXJK(bcai) 
C
C             - xop(cd) t0_jk(bdai)
C                             
C             - xop(bd) t0_jk(dcai)
C                             

      IMPLICIT NONE
C
      INTEGER IST0JK, ISYMXOP, ISTETAXJK, LWORK
      INTEGER KAD, KEND1, LWRK1, KOFF1, KOFF2, KOFF3
      INTEGER ISYMI, ISYMBD, ISYMBDA, ISYMD, ISYMA 
      INTEGER ISYMC, ISYMBCA, ISYMBC
      INTEGER ISYMB, ISYMDCA, ISYMDC 

      INTEGER NTOTC, NTOTB, NTOTD
C
#if defined (SYS_CRAY)
      REAL TB0JK(*), TETAXJK(*), XOP(*), WORK(LWORK)
      REAL ONE
#else
      DOUBLE PRECISION T0JK(*), TETAXJK(*), XOP(*), WORK(LWORK) 
      DOUBLE PRECISION ONE
#endif
C
      PARAMETER (ONE = 1.0D0)
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccsdinp.h>
C
      CALL QENTER('TETAX_JK_BC')

      KAD  = 1
      KEND1  = KAD + NMATAB(ISYMXOP)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWRK1
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in TETAX_JK_BC')
      END IF
C
C SORT VIR-VIR  XOP ELEMENTS (A,D)
C
C
      DO ISYMD = 1,NSYM
         ISYMA = MULD2H(ISYMD,ISYMXOP)
         DO D = 1,NVIR(ISYMD)
            KOFF1 = IFCVIR(ISYMA,ISYMD) + NORB(ISYMA)*(D - 1)
     *                                  + NRHF(ISYMA) + 1
            KOFF2 = KAD + IMATAB(ISYMA,ISYMD) + NVIR(ISYMA)*(D - 1)
            CALL DCOPY(NVIR(ISYMA),XOP(KOFF1),1,WORK(KOFF2),1)
         END DO
      END DO
C
C TETAXJK(bcai) = TETAXJK(bcai) - xop(cd) t0_jk(bdai) 
C
      DO ISYMI = 1,NSYM
         ISYMBDA = MULD2H(IST0JK,ISYMI)
         DO ISYMA = 1,NSYM
            ISYMBD = MULD2H(ISYMBDA,ISYMA)
            DO ISYMD = 1,NSYM
            ISYMC = MULD2H(ISYMD,ISYMXOP)
            ISYMB = MULD2H(ISYMBD,ISYMD)
            ISYMBC  = MULD2H(ISYMB,ISYMC)
            ISYMBCA = MULD2H(ISYMBC,ISYMA)
               DO I = 1,NRHF(ISYMI)
                  DO A = 1,NVIR(ISYMA)
C
                     KOFF1   = 1
     *                        + IMAABCI(ISYMBDA,ISYMI)
     *                        + NMAABC(ISYMBDA)*(I-1)
     *                        + IMAABC(ISYMBD,ISYMA)
     *                        + NMATAB(ISYMBD)*(A-1)
     *                        + IMATAB(ISYMB,ISYMD)
                     KOFF2   = KAD
     *                        + IMATAB(ISYMC,ISYMD)
                     KOFF3   = 1
     *                        + IMAABCI(ISYMBCA,ISYMI)
     *                        + NMAABC(ISYMBCA)*(I-1)
     *                        + IMAABC(ISYMBC,ISYMA)
     *                        + NMATAB(ISYMBC)*(A-1)
     *                        + IMATAB(ISYMB,ISYMC)
C
                     NTOTB = MAX(1,NVIR(ISYMB))
                     NTOTC  = MAX(1,NVIR(ISYMC))
C
C TETAXJK(bcai) = TETAXJK(bcai) - xop(cd) t0_jk(bdai) 
C
                     CALL DGEMM('N','T',NVIR(ISYMB),NVIR(ISYMC),
     *                          NVIR(ISYMD),-ONE,T0JK(KOFF1),NTOTB,
     *                          WORK(KOFF2),NTOTC,
     *                          ONE,TETAXJK(KOFF3),NTOTB)
                  END DO
               END DO
            END DO
         END DO
      END DO
C
C TETAXJK(bcai) = TETAXJK(bcai) - xop(bd) t0_jk(dcai) 
C     
      DO ISYMI = 1,NSYM
         ISYMDCA = MULD2H(IST0JK,ISYMI)
         DO ISYMA = 1,NSYM
            ISYMDC = MULD2H(ISYMDCA,ISYMA)
            DO ISYMC = 1,NSYM
               ISYMD = MULD2H(ISYMDC,ISYMC)
               ISYMB = MULD2H(ISYMD,ISYMXOP)
               ISYMBC  = MULD2H(ISYMB,ISYMC)
               ISYMBCA = MULD2H(ISYMBC,ISYMA)
                  DO I = 1,NRHF(ISYMI)
                     DO A = 1,NVIR(ISYMA)
C
                     KOFF1   = KAD
     *                        + IMATAB(ISYMB,ISYMD)
                     KOFF2   = 1
     *                        + IMAABCI(ISYMDCA,ISYMI)
     *                        + NMAABC(ISYMDCA)*(I-1)
     *                        + IMAABC(ISYMDC,ISYMA)
     *                        + NMATAB(ISYMDC)*(A-1)
     *                        + IMATAB(ISYMD,ISYMC)
                     KOFF3   = 1
     *                        + IMAABCI(ISYMBCA,ISYMI)
     *                        + NMAABC(ISYMBCA)*(I-1)
     *                        + IMAABC(ISYMBC,ISYMA)
     *                        + NMATAB(ISYMBC)*(A-1)
     *                        + IMATAB(ISYMB,ISYMC)
C     
                     NTOTB = MAX(1,NVIR(ISYMB))
                     NTOTD  = MAX(1,NVIR(ISYMD))
C
C TETAXJK(bcai) = TETAXJK(bcai) - xop(bd) t0_jk(dcai) 
C     
                     CALL DGEMM('N','N',NVIR(ISYMB),NVIR(ISYMC),
     *                          NVIR(ISYMD),-ONE,WORK(KOFF1),NTOTB,
     *                          T0JK(KOFF2),NTOTD,
     *                          ONE,TETAXJK(KOFF3),NTOTB)
                  END DO
               END DO
            END DO
         END DO
      END DO

      CALL QEXIT('TETAX_JK_BC')
      RETURN
      END
C  /* Deck tetax_jk_a */
      SUBROUTINE TETAX_JK_A(T0JK,IST0JK,XOP,ISYMXOP,TETAXJK,ISTETAXJK,
     *                      WORK,LWORK)
C
C TETAXJK(bcai) = TETAXJK(bcai) 
C
C             - xop(ad) t0_jk(bcdi)
C                             

      IMPLICIT NONE
C
      INTEGER IST0JK, ISYMXOP, ISTETAXJK, LWORK
      INTEGER KAD, KEND1, LWRK1, KOFF1, KOFF2, KOFF3
      INTEGER ISYMI, ISYMBCD, ISYMD, ISYMA, ISYMBCA, ISYMBC
      INTEGER NTOTBC, NTOTA
C
#if defined (SYS_CRAY)
      REAL TB0JK(*), TETAXJK(*), XOP(*), WORK(LWORK)
      REAL ONE
#else
      DOUBLE PRECISION T0JK(*), TETAXJK(*), XOP(*), WORK(LWORK) 
      DOUBLE PRECISION ONE
#endif
C
      PARAMETER (ONE = 1.0D0)
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccsdinp.h>
C
      CALL QENTER('TETAX_JK_A')

      KAD  = 1
      KEND1  = KAD + NMATAB(ISYMXOP)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWRK1
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in TETAX_JK_A')
      END IF
C
C SORT VIR-VIR  XOP ELEMENTS (A,D)
C
C
      DO ISYMD = 1,NSYM
         ISYMA = MULD2H(ISYMD,ISYMXOP)
         DO D = 1,NVIR(ISYMD)
            KOFF1 = IFCVIR(ISYMA,ISYMD) + NORB(ISYMA)*(D - 1)
     *                                  + NRHF(ISYMA) + 1
            KOFF2 = KAD + IMATAB(ISYMA,ISYMD) + NVIR(ISYMA)*(D - 1)
            CALL DCOPY(NVIR(ISYMA),XOP(KOFF1),1,WORK(KOFF2),1)
         END DO
      END DO
C
      DO ISYMI = 1,NSYM
         ISYMBCD = MULD2H(IST0JK,ISYMI)
         DO I = 1,NRHF(ISYMI)
            DO ISYMD = 1,NSYM
               ISYMA = MULD2H(ISYMD,ISYMXOP)
               ISYMBCA = MULD2H(ISYMXOP,ISYMBCD)
               ISYMBC  = MULD2H(ISYMBCD,ISYMD)
               KOFF1   = 1
     *                + IMAABCI(ISYMBCD,ISYMI)
     *                + NMAABC(ISYMBCD)*(I-1)
     *                + IMAABC(ISYMBC,ISYMD)
               KOFF2   = KAD
     *                + IMATAB(ISYMA,ISYMD)
               KOFF3   = 1
     *                + IMAABCI(ISYMBCA,ISYMI)
     *                + NMAABC(ISYMBCA)*(I-1)
     *                + IMAABC(ISYMBC,ISYMA)
C
               NTOTBC = MAX(1,NMATAB(ISYMBC))
               NTOTA  = MAX(1,NVIR(ISYMA))
C
C TETAXJK(bcai) = TETAXJK(bcai)  - xop(ad) tb0_jk(bcdi)  
C
               CALL DGEMM('N','T',NMATAB(ISYMBC),NVIR(ISYMA),
     *                    NVIR(ISYMD),-ONE,T0JK(KOFF1),NTOTBC,
     *                    WORK(KOFF2),NTOTA,
     *                    ONE,TETAXJK(KOFF3),NTOTBC)
            END DO
         END DO
      END DO
C                             
      CALL QEXIT('TETAX_JK_A')
      RETURN
      END
C  /* Deck aden_dij_jk */
      SUBROUTINE ADEN_DIJ_JK(DIJ,THLM,ISYMTHLM,WLM,ISYMWLM)
C
      IMPLICIT NONE
#include <priunit.h>
#include <dummy.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYMTHLM,ISYMWLM
      INTEGER ISYMJ,ISYMDEF,ISYMI,KOFF1,KOFF2,KOFF3,NTOTDEF,NTOTI
C
#if defined (SYS_CRAY)
      REAL DIJ(*),THLM(*),WLM(*)
      REAL ONE,HALF
#else
      DOUBLE PRECISION DIJ(*),THLM(*),WLM(*)
      DOUBLE PRECISION ONE,HALF
#endif
C
      PARAMETER (HALF = 0.5D0, ONE = 1.0D0)
C
      CALL QENTER('ADEN_DIJ_JK')

C     D(ij) = THETA^LM(defi)*W^LM(defj) 
C
      DO ISYMJ = 1,NSYM
         ISYMDEF = MULD2H(ISYMWLM,ISYMJ)
         ISYMI = MULD2H(ISYMTHLM,ISYMDEF)
C
         KOFF1 = IMAABCI(ISYMDEF,ISYMI) + 1
         KOFF2 = IMAABCI(ISYMDEF,ISYMJ) + 1
         KOFF3 = IMATIJ(ISYMI,ISYMJ)    + 1
C
         NTOTDEF = MAX(NMAABC(ISYMDEF),1)
         NTOTI   = MAX(NRHF(ISYMI),1)
C
c        write(lupri,*)'THLM(vvv,o) in ADEN_DIJ_JK'
c        call output(THLM(KOFF1),1,NMAABC(ISYMDEF),1,NRHF(ISYMI),
c    *               NMAABC(ISYMDEF),NRHF(ISYMI),1,lupri)
c        write(lupri,*)'WLM(vvv,o) in ADEN_DIJ_JK'
c        call output(WLM(KOFF2),1,NMAABC(ISYMDEF),1,NRHF(ISYMI),
c    *               NMAABC(ISYMDEF),NRHF(ISYMI),1,lupri)
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NMAABC(ISYMDEF),
     *              HALF,THLM(KOFF1),NTOTDEF,WLM(KOFF2),NTOTDEF,
     *              ONE,DIJ(KOFF3),NTOTI)
C
c        write(lupri,*)'DIJ(o,o) in ADEN_DIJ_JK'
c        call output(DIJ(KOFF3),1,NRHF(ISYMI),1,NRHF(ISYMJ),
c    *               NRHF(ISYMI),NRHF(ISYMI),1,lupri)
C
      END DO
C
      CALL QEXIT('ADEN_DIJ_JK')
C
      RETURN
      END
C  /* Deck aden_dab_lm */
      SUBROUTINE ADEN_DAB_LM(DAB,THLM,ISYMTHLM,WLM,ISYMWLM,WORK,LWORK)
C
      IMPLICIT NONE
#include <priunit.h>
#include <dummy.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYMTHLM,ISYMWLM,LWORK
      INTEGER ISYMN,ISYMDEB,ISYMDEA,ISYMB,ISYMDE,ISYMA
      INTEGER KOFF1,KOFF2,KOFF3
      INTEGER NTOTDE,NTOTA
      INTEGER KWDAEN,KTHDBEN,KEND1,LWRK1
      INTEGER ISYMDBE,ISYMDAE,ISYME,ISYMDB,ISYMDA,ISYMEN,ISYMD
      INTEGER NTOTD
C
#if defined (SYS_CRAY)
      REAL DAB(*),THLM(*),WLM(*)
      REAL ONE,HALF
      REAL WORK(LWORK)
#else
      DOUBLE PRECISION DAB(*),THLM(*),WLM(*)
      DOUBLE PRECISION ONE,HALF
      DOUBLE PRECISION WORK(LWORK)
#endif
C
      PARAMETER (HALF = 0.5D0, ONE = 1.0D0)
C
      CALL QENTER('ADEN_DAB_LM')

C
C     D(ab) = W^LM(dean) * THETA^LM(debn)
C
      DO ISYMN = 1,NSYM
         ISYMDEB = MULD2H(ISYMTHLM,ISYMN)
         ISYMDEA = MULD2H(ISYMWLM,ISYMN)
         DO ISYMB = 1,NSYM
            ISYMDE = MULD2H(ISYMDEB,ISYMB)
            ISYMA  = MULD2H(ISYMDEA,ISYMDE)
            DO N = 1,NRHF(ISYMN)
C
               KOFF1 = IMAABCI(ISYMDEA,ISYMN) 
     *               + NMAABC(ISYMDEA)*(N-1)
     *               + IMAABC(ISYMDE,ISYMA)
     *               + 1
               KOFF2 = IMAABCI(ISYMDEB,ISYMN) 
     *               + NMAABC(ISYMDEB)*(N-1)
     *               + IMAABC(ISYMDE,ISYMB)
     *               + 1
               KOFF3 = IMATAB(ISYMA,ISYMB) + 1
C
               NTOTDE = MAX(NMATAB(ISYMDE),1)
               NTOTA   = MAX(NVIR(ISYMA),1)
C
c     write(lupri,*)'THLM(vv,v) in ADEN_DAB_LM'
c     call output(THLM(KOFF2),1,NMATAB(ISYMDE),1,NVIR(ISYMB),
c    *            NMATAB(ISYMDE),NVIR(ISYMB),1,lupri)
c     write(lupri,*)'WLM(vv,v) in ADEN_DAB_LM'
c     call output(WLM(KOFF1),1,NMATAB(ISYMDE),1,NVIR(ISYMA),
c    *            NMATAB(ISYMDE),NVIR(ISYMA),1,lupri)
C
               CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),
     *                    NMATAB(ISYMDE),-HALF,WLM(KOFF1),NTOTDE,
     *                    THLM(KOFF2),NTOTDE,ONE,DAB(KOFF3),NTOTA)
C
c     write(lupri,*)'DAB(v,v) in ADEN_DAB_LM'
c     call output(DAB(KOFF3),1,NVIR(ISYMA),1,NVIR(ISYMB),
c    *            NVIR(ISYMA),NVIR(ISYMB),1,lupri)
C
            END DO   ! N
         END DO      ! ISYMB
      END DO         ! ISYMN
C
C     Calculate second contribution to D(ab)
C
      KWDAEN = 1
      KTHDBEN  = KWDAEN + NMAABCI(ISYMWLM)
      KEND1   = KTHDBEN + NMAABCI(ISYMTHLM)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWRK1
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in ADEN_DAB_LM')
      END IF
C
      CALL DZERO(WORK(KWDAEN),NMAABCI(ISYMWLM))
      CALL DZERO(WORK(KTHDBEN),NMAABCI(ISYMTHLM))
C
C     Sort W^LM(dean) to W^LM(daen)
C
      CALL FACBI(WORK(KWDAEN),WLM,ISYMWLM)
C
C     Sort THETA^LM(debn) to THETA^LM(dben)
C
      CALL FACBI(WORK(KTHDBEN),THLM,ISYMTHLM)
C
C     D(ab) = W^LM(daen) * THETA^LM(dben)
C
      DO ISYMN = 1,NSYM
         ISYMDEB = MULD2H(ISYMTHLM,ISYMN)
         ISYMDEA = MULD2H(ISYMWLM,ISYMN)
         DO ISYMB = 1,NSYM
            ISYMDE = MULD2H(ISYMDEB,ISYMB)
            ISYMA  = MULD2H(ISYMDEA,ISYMDE)
            DO N = 1,NRHF(ISYMN)
C
               KOFF1 = IMAABCI(ISYMDEA,ISYMN) 
     *               + NMAABC(ISYMDEA)*(N-1)
     *               + IMAABC(ISYMDE,ISYMA)
     *               + KWDAEN
               KOFF2 = IMAABCI(ISYMDEB,ISYMN) 
     *               + NMAABC(ISYMDEB)*(N-1)
     *               + IMAABC(ISYMDE,ISYMB)
     *               + KTHDBEN
               KOFF3 = IMATAB(ISYMA,ISYMB) + 1
C
               NTOTDE = MAX(NMATAB(ISYMDE),1)
               NTOTA   = MAX(NVIR(ISYMA),1)
C
c     write(lupri,*)'THLM(vv,v) in ADEN_DAB_LM'
c     call output(WORK(KOFF2),1,NMATAB(ISYMDE),1,NVIR(ISYMB),
c    *            NMATAB(ISYMDE),NVIR(ISYMB),1,lupri)
c     write(lupri,*)'WLM(vv,v) in ADEN_DAB_LM'
c     call output(WORK(KOFF1),1,NMATAB(ISYMDE),1,NVIR(ISYMA),
c    *            NMATAB(ISYMDE),NVIR(ISYMA),1,lupri)
C
               CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),
     *                    NMATAB(ISYMDE),-ONE,WORK(KOFF1),NTOTDE,
     *                    WORK(KOFF2),NTOTDE,ONE,DAB(KOFF3),NTOTA)
C
c     write(lupri,*)'DAB(v,v) in ADEN_DAB_LM'
c     call output(DAB(KOFF3),1,NVIR(ISYMA),1,NVIR(ISYMB),
c    *            NVIR(ISYMA),NVIR(ISYMB),1,lupri)
C
            END DO   ! N
         END DO      ! ISYMB
      END DO         ! ISYMN
C
      CALL QEXIT('ADEN_DAB_LM')
C
      RETURN
      END
C  /* Deck aden_dai_lm */
      SUBROUTINE ADEN_DAI_LM(DAI,
     *                       T2TP,ISYMT2,
     *                       TETAXKL,ISTETAXKL,
     *                       ISYML,L,ISYMM,M,
     *                       WORK,LWORK)
C
C d(ia) = d(ia) + t2tp(dLMe) * ( teta^LM(deai) - teta^LM(daei) )
C                    I^LM(de)         I^LM(deai)
C
      IMPLICIT NONE
#include <priunit.h>
#include <dummy.h>
#include <ccsdsym.h> 
#include <ccorb.h>
C
      INTEGER ISYMT2,ISTETAXKL,ISYML,ISYMM,LWORK
      INTEGER KT2DE,KTHDAEI,KEND1,LWRK1
      INTEGER ISYMI,ISYMDEA,ISYMDE,ISYMA,ISYMD,ISYME
      INTEGER KOFF1,KOFF2,KOFF3
      INTEGER NTOTDE
C
#if defined (SYS_CRAY)
      REAL DAI(*),T2TP(*),TETAXKL(*),WORK(LWORK)
      REAL ONE
#else
      DOUBLE PRECISION DAI(*),T2TP(*),TETAXKL(*),WORK(LWORK)
      DOUBLE PRECISION ONE
#endif
C
      PARAMETER (ONE = 1.0D0)
C
      CALL QENTER('ADEN_DAI_LM')
C
      KT2DE = 1
      KTHDAEI = KT2DE + NT2SQ(ISYMT2)
      KEND1   = KTHDAEI + NMAABCI(ISTETAXKL)
      LWRK1   = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in ADEN_DAI_LM')
      END IF
C
      CALL DZERO(WORK(KT2DE),NT2SQ(ISYMT2))
      CALL DZERO(WORK(KTHDAEI),NMAABCI(ISTETAXKL))      
C
      CALL SORT_T2_AB(WORK(KT2DE),ISYML,L,ISYMM,M,T2TP,ISYMT2)
C
      CALL DAXPY(NMAABCI(ISTETAXKL),-ONE,TETAXKL,1,WORK(KTHDAEI),1)
C 
      CALL FACBI(TETAXKL,WORK(KTHDAEI),ISTETAXKL)
C
      DO ISYMI = 1,NSYM
         ISYMDEA = MULD2H(ISTETAXKL,ISYMI)
         DO ISYMA = 1,NSYM
            ISYMDE = MULD2H(ISYMDEA,ISYMA)
            DO ISYMD = 1,NSYM
               ISYME = MULD2H(ISYMDE,ISYMD)
               DO I = 1,NRHF(ISYMI)
 
                  KOFF1 = IMAABCI(ISYMDEA,ISYMI)
     *                  + NMAABC(ISYMDEA)*(I-1)
     *                  + IMAABC(ISYMDE,ISYMA)
     *                  + 1
                  KOFF2 = KT2DE 
     *                  + IMATAB(ISYMD,ISYME)
 
                  KOFF3 = IT1AM(ISYMA,ISYMI)
     *                  + 1
 
                  NTOTDE = MAX(1,NMATAB(ISYMDE))
 
                  CALL DGEMV('T',NMATAB(ISYMDE),NVIR(ISYMA),-ONE,
     *                       TETAXKL(KOFF1),NTOTDE,WORK(KOFF2),1,
     *                       ONE,DAI(KOFF3),1)
               END DO
            END DO
         END DO
      END DO
C
      CALL QEXIT('ADEN_DAI_LM')
C
      RETURN
      END

