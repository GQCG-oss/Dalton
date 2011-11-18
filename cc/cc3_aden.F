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
C  /* Deck cc3_aden */
      SUBROUTINE CC3_ADEN(LISTL,IDLSTL,LISTR,IDLSTR,
     *                            XLAMDP0,XLAMDH0,FOCK0,
     *                            DIJ,DAB,DO_DIA,DIA,ISYDEN,
     *                            DO_YMMAT,YMMAT,
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
#include <ccsdinp.h>
C
      LOGICAL DO_YMMAT, DO_DIA
C
      CHARACTER LISTL*3, LISTR*3
      CHARACTER*(*) FNTOC, FN3VI, FNDKBC3, FN3FOPX, FN3FOP2X
      CHARACTER*(*) FNDKBC,FNDELD,FNCKJD
      CHARACTER*5 FN3FOP
      CHARACTER*8 FN3VI2
      CHARACTER*6 FN3FOP2
      CHARACTER*10 MODEL
C
      PARAMETER (FN3FOP  = 'PTFOP')
      PARAMETER (FN3VI2  = 'CC3_VI12')
      PARAMETER (FN3FOP2 = 'PTFOP2')
C
      INTEGER ISYDEN,IDLSTL,IDLSTR,LWORK
      INTEGER LUTOC, LU3VI, LUDKBC3, LU3FOPX, LU3FOP2X
      INTEGER LUDKBC,LUDELD,LUCKJD
      INTEGER LU3FOP
      INTEGER LU3VI2, LU3FOP2
      INTEGER ISYM0,KT1AMP,KLAMP0,KLAMH0,KEND1,LWRK1,IOPT
C
#if defined (SYS_CRAY)
      REAL XLAMDP0(*),XLAMDH0(*),FOCK0(*)
      REAL DAB(*),DIJ(*),DIA(*)
      REAL YMMAT(*)
      REAL WORK(LWORK)
#else
      DOUBLE PRECISION XLAMDP0(*),XLAMDH0(*),FOCK0(*)
      DOUBLE PRECISION DAB(*),DIJ(*),DIA(*)
      DOUBLE PRECISION YMMAT(*)
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
      CALL DZERO(DAB,NMATAB(ISYDEN))
      CALL DZERO(DIJ,NMATIJ(ISYDEN))
      IF (DO_DIA) THEN
         CALL DZERO(DIA,NT1AM(ISYDEN))
      END IF
C
C     Open the file
C
      LU3FOP  = -1
      LU3VI2  = -1
      LU3FOP2 = -1
      CALL WOPEN2(LU3FOP,FN3FOP,64,0)
      CALL WOPEN2(LU3VI2,FN3VI2,64,0)
      CALL WOPEN2(LU3FOP2,FN3FOP2,64,0)
C
      CALL CC3_ADENVIR(DIJ,DAB,DO_DIA,DIA,ISYDEN,
     *                 DO_YMMAT,YMMAT,
     *                   LISTL,IDLSTL,LISTR,IDLSTR,
     *                   LUTOC,FNTOC,LU3VI,FN3VI,LU3VI2,FN3VI2,
     *                   LUDKBC3,FNDKBC3,
     *                   LU3FOPX,FN3FOPX,LU3FOP2X,FN3FOP2X,
     *                   LU3FOP,FN3FOP,LU3FOP2,FN3FOP2,
     *                   LUCKJD,FNCKJD,LUDKBC,FNDKBC,LUDELD,FNDELD,
     *                   WORK(KEND1),LWRK1)
C
      IF (IPRINT .GT. 55) THEN
         WRITE(LUPRI,*)'DAB density after CC3_ADENVIR '
         CALL PRINT_MATAB(DAB,ISYDEN)
         WRITE(LUPRI,*)'DIJ density after CC3_ADENVIR '
         CALL PRINT_MATIJ(DIJ,ISYDEN)
         IF (DO_DIA) THEN
            WRITE(LUPRI,*)'DIA density after CC3_ADENVIR '
            CALL PRINT_MATAI(DIA,ISYDEN)
         END IF
      END IF
C
      IF (LISTR(1:3).EQ.'R1 ') THEN
         CALL CC3_ADENOCC(LISTL,IDLSTL,LISTR,IDLSTR,
     *                               WORK(KLAMP0),WORK(KLAMH0),FOCK0,
     *                               DIJ,DAB,DO_DIA,DIA,ISYDEN,
     *                               WORK(KEND1),LWRK1,
     *                               LUDELD,FNDELD,LUCKJD,FNCKJD,LUDKBC,
     *                               FNDKBC,LUTOC,FNTOC,LU3VI,FN3VI,
     *                               LUDKBC3,FNDKBC3,LU3FOP,FN3FOP,
     *                               LU3FOPX,FN3FOPX,
     *                               LU3FOP2X,FN3FOP2X)
 
         IF (IPRINT .GT. 55) THEN
            WRITE(LUPRI,*)'DAB density after CC3_ADENOCC '
            CALL PRINT_MATAB(DAB,ISYDEN)
            WRITE(LUPRI,*)'DIJ density after CC3_ADENOCC '
            CALL PRINT_MATIJ(DIJ,ISYDEN)
            IF (DO_DIA) THEN 
               WRITE(LUPRI,*)'DIA density after CC3_ADENOCC '
               CALL PRINT_MATAI(DIA,ISYDEN)
            END IF
         END IF
 
      END IF
C
      CALL WCLOSE2(LU3FOP,FN3FOP,'KEEP')
      CALL WCLOSE2(LU3VI2,FN3VI2,'KEEP')
      CALL WCLOSE2(LU3FOP2,FN3FOP2,'KEEP')
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
     *                            DIJ,DAB,DO_DIA,DIA,ISYDEN,
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
#include <ccl1rsp.h>
#include <ccr1rsp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <dummy.h>
#include <priunit.h>
#include <iratdef.h>
#include <ccinftap.h>
#include <ccsdinp.h>
#include <cclrmrsp.h>
#include <ccexci.h>
#include <ccn2rsp.h>
C
      LOGICAL DO_DIA
      LOGICAL LSKIPL1R
C
      INTEGER ISYM0
      PARAMETER(ISYM0 = 1)
C
      INTEGER ISYDEN,IDLSTL,IDLSTR,IDLSTL0,LWORK
      INTEGER LUTOC, LU3VI, LUDKBC3, LU3FOPX, LU3FOP2X, LU3FOP
      INTEGER LUDKBC,LUDELD,LUCKJD
      INTEGER ISYML0,ISYML1,ISYMR1,ISINT1,ISINT2,ISINTR1,ISYFCKL1R
      INTEGER ISYMK,ISYML,ISYMT3,ISYMKL,ISYT30KL
      INTEGER IOPT,LENGTH
      INTEGER KFOCKD,KFCKBA,KT2TP,KL1AM,KL2TP,KEND0,LWRK0
      INTEGER KL1,KL2,KFOCKL1,KR1,KR2,KFOCKR1,KEND1,LWRK1
      INTEGER KXIAJB,KT3BOG1,KT3BOL1,KT3BOG2,KT3BOL2,KT3OG1,KT3OG2
      INTEGER KLAMPL1R,KLAMHL1R,KT30KL
      INTEGER KFOCKL1RCK,KT3VIJG1
      INTEGER ISYMT3B,ISYT3B0KL,ISYW3BXKL
      INTEGER KXGADCK,KXLADCK
      INTEGER KT3B0KL,KW3BXKL,ISYMW3BX
      INTEGER KT3BOG2X,KT3BOL2X,KXGADCKX,KXLADCKX
      INTEGER ISYMTETAX,ISTETAXKL
      INTEGER KTETAXKL
      INTEGER IDLSTL1R,ISYML1R
      INTEGER ISINT2L1R,KT1L1R
      INTEGER KEND2,LWRK2
C
      INTEGER IR1TAMP
      INTEGER ILSTSYM
C
      CHARACTER LISTL*3, LISTR*3, LISTL0*3, LISTL1R*3
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
      integer isymi,isyabc,koff1
C
#if defined (SYS_CRAY)
      REAL XLAMDP0(*),XLAMDH0(*),FOCK0(*) 
      REAL DAB(*),DIJ(*),DIA(*)
      REAL WORK(LWORK)
      REAL FREQL1,FREQR1,FREQL1R
      REAL DDOT,XNORMVAL,ONE
#else
      DOUBLE PRECISION XLAMDP0(*),XLAMDH0(*),FOCK0(*) 
      DOUBLE PRECISION DAB(*),DIJ(*),DIA(*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION FREQL1,FREQR1,FREQL1R
      DOUBLE PRECISION DDOT,XNORMVAL,ONE
#endif
C
      PARAMETER (ONE = 1.0D0)
C
      CALL QENTER('CC3_ADENOCC')
C
C------------------------------------------------------------
C     some initializations:
C------------------------------------------------------------
C
*     LISTL0 = 'L0 '
*     IDLSTL0 = 0
*     ISYML0 = ISYM0
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

        LISTL1R  = 'R1 '
        IDLSTL1R = IR1TAMP(LABELL1,LORXL1,FREQL1,ISYML1)
        ! get symmetry and frequency from common blocks
        ! defined in ccl1rsp.h
        ISYML1R  = ISYLRT(IDLSTL1R)
        FREQL1R  = FRQLRT(IDLSTL1R)
C
        !LITSL0 corresponding to LISTL
        LISTL0 = 'L0 '
        IDLSTL0 = 0
        ISYML0 = ISYM0
C
        IF (FREQL1R .NE. FREQL1) THEN
           WRITE(LUPRI,*)'FREQL1R: ', FREQL1R
           WRITE(LUPRI,*)'FREQL1: ', FREQL1
           CALL QUIT('Frequency mismatch in CC3_ADENOCC(L1)')
        END IF

      ELSE IF (LISTL(1:3).EQ.'M1 ') THEN
        ISYML1 = ILSTSYM(LISTL,IDLSTL)
        FREQL1 = FRQLRM(IDLSTL)
        LABELL1 = '- none -'
C
        ! find corresponding right eigenvector
        LISTL1R = 'RE '
        IDLSTL1R = ILRM(IDLSTL)
        ISYML1R = ISYML1
        FREQL1R = EIGVAL(IDLSTL1R)
C
        !LITSL0 corresponding to LISTL
        LISTL0 = 'L0 '
        IDLSTL0 = 0
        ISYML0 = ISYM0
C
        IF (FREQL1R .NE. FREQL1) THEN
           WRITE(LUPRI,*)'FREQL1R: ', FREQL1R
           WRITE(LUPRI,*)'FREQL1: ', FREQL1
           CALL QUIT('Frequency mismatch in CC3_ADENOCC(M1)')
        END IF
C
      ELSE IF (LISTL(1:3).EQ.'N2 ') THEN
        ISYML1 = ILSTSYM(LISTL,IDLSTL)
        FREQL1 = FRQIN2(IDLSTL) + FRQFN2(IDLSTL)
        LABELL1 = '- none -'
C
        ! find corresponding right eigenvector
        LISTL1R = 'RE '
        IDLSTL1R = IFN2(IDLSTL)
        ISYML1R = ILSTSYM(LISTL1R,IDLSTL1R)
        FREQL1R = FRQFN2(IDLSTL)
C
        !LITSL0 corresponding to LISTL
        LISTL0 = 'LE '
        IDLSTL0 = IIN2(IDLSTL)
        ISYML0 = ILSTSYM(LISTL0,IDLSTL0)
C
      ELSE IF (LISTL(1:3).EQ.'LE ') THEN
        ISYML1 = ILSTSYM(LISTL,IDLSTL)
        FREQL1 = -EIGVAL(IDLSTL)
        LABELL1 = '- none -'
C
        !we don't have any "right" vector entering a right hand side
        LISTL1R = '---'
        IDLSTL1R = -99
C
        !LITSL0 corresponding to LISTL (not used for LE)
        LISTL0 = 'L0 '
        IDLSTL0 = 0
        ISYML0 = ISYM0
C
      ELSE 
         CALL QUIT('Unknown left list in CC3_ADENOCC')
      END IF

      IF (LISTR(1:3).EQ.'R1 ') THEN
         ! get symmetry, frequency and integral label for right list 
         ! from common blocks defined in ccr1rsp.h
        ISYMR1  = ISYLRT(IDLSTR)
        FREQR1  = FRQLRT(IDLSTR)
        LABELR1 = LRTLBL(IDLSTR)
      ELSE
         CALL QUIT('Unknown right list in CC3_ADENOCC')
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
     *                  ISYM0,XLAMDH0,ISYM0,WORK(KEND1),LWRK1)

      END IF
C
C-------------------------------------
C     Read L1 and L2 multipliers 
C-------------------------------------
C
      IOPT  = 3
      CALL GET_T1_T2(IOPT,.FALSE.,WORK(KL1),WORK(KL2),LISTL,
     *               IDLSTL,ISYML1,WORK(KEND1),LWRK1)
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
     *                  ISYM0,XLAMDH0,ISYM0,WORK(KEND1),LWRK1)
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
C---------------------------------------------------
C If we want to sum the T3 amplitudes (for debugging)
C---------------------------------------------------
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
      IF ( (LISTL(1:3).EQ.'L1 ') .OR. (LISTL(1:3).EQ.'M1 ') 
     *                           .OR. (LISTL(1:3).EQ.'N2 ') ) THEN
         ISINT2L1R = MULD2H(ISYML1R,ISINT2)
         ISYFCKL1R  = MULD2H(ISYMOP,ISYML1R)
      END IF

      KXIAJB    = KEND1
      KEND1     = KXIAJB    + NT2AM(ISYM0)

      KT3BOG1   = KEND1
      KT3BOL1   = KT3BOG1   + NTRAOC(ISYM0)
      KT3BOG2   = KT3BOL1   + NTRAOC(ISYM0)
      KT3BOL2   = KT3BOG2   + NTRAOC(ISYM0)
      KT3OG1    = KT3BOL2   + NTRAOC(ISYM0)
      KT3OG2    = KT3OG1    + NTRAOC(ISINT2)
      KLAMPL1R   = KT3OG2    + NTRAOC(ISINT2)
      KLAMHL1R   = KLAMPL1R   + NLAMDT
      KEND1     = KLAMHL1R   + NLAMDT
C
      KT3VIJG1 = KEND1
      KEND1    = KT3VIJG1  + NMAABCI(ISYM0)
      LWRK1     = LWORK     - KEND1
C
      IF ( (LISTL(1:3).EQ.'L1 ') .OR. (LISTL(1:3).EQ.'M1 ') 
     *                           .OR. (LISTL(1:3).EQ.'N2 ') ) THEN
         KFOCKL1RCK  = KEND1
         KEND1  = KFOCKL1RCK  + NT1AM(ISYFCKL1R)
         LWRK1     = LWORK     - KEND1
C
         KT3BOG2X   = KEND1
         KT3BOL2X   = KT3BOG2X + NTRAOC(ISINT2L1R)
         KEND1      = KT3BOL2X + NTRAOC(ISINT2L1R)
C
         KXGADCKX   = KEND1 
         KXLADCKX   = KXGADCKX + NMAABCI(ISINT2L1R)
         KEND1     = KXLADCKX + NMAABCI(ISINT2L1R)
         LWRK1     = LWORK     - KEND1
C
         KT1L1R  = KEND1
         KEND1  = KT1L1R + NT1AM(ISYML1R)
         LWRK1   = LWORK  - KEND1
C
      END IF
C
      KXGADCK   = KEND1 
      KXLADCK   = KXGADCK + NMAABCI(ISYM0)
      KEND1     = KXLADCK + NMAABCI(ISYM0)
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
C----------------------------------------------------------
C     Get Lambda for right list depended on left LISTL list
C----------------------------------------------------------
C
      IF ( (LISTL(1:3).EQ.'L1 ') .OR. (LISTL(1:3).EQ.'M1 ') 
     *                           .OR. (LISTL(1:3).EQ.'N2 ') ) THEN
         CALL GET_LAMBDAX(WORK(KLAMPL1R),WORK(KLAMHL1R),LISTL1R,
     *                    IDLSTL1R,
     *                    ISYML1R,XLAMDP0,XLAMDH0,WORK(KEND1),LWRK1)
C
C------------------------------------------------------------------
C        Calculate the F^L1R matrix (kc elements evaluated and stored 
C        as ck)
C------------------------------------------------------------------
C
         IOPT = 1
         CALL GET_T1_T2(IOPT,.FALSE.,WORK(KT1L1R),DUMMY,LISTL1R,
     *                  IDLSTL1R,
     *                  ISYML1R,WORK(KEND1),LWRK1)
         CALL CC3LR_MFOCK(WORK(KFOCKL1RCK),WORK(KT1L1R),WORK(KXIAJB),
     *                    ISYFCKL1R)
C
      END IF
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

C     KW3BXOG1 and KW3BXOL1 have to be taken out
C     and ccfop_sort included
C
!     CALL INTOCC_T3BARX(.FALSE.,
!    *                   LUTOC,FNTOC,ISYMOP,XLAMDH0,ISYM0,ISINT2,
!    *                   WORK(KLAMHL1R),ISYML1R,ISINT2L1R,
!    *                   WORK(KW3BXOG1),WORK(KW3BXOL1),
!    *                   WORK(KT3BOG2X1),WORK(KT3BOL2X1),
!    *                   WORK(KEND1),LWRK1)
C
c     CALL CCFOP_SORT(WORK(KT3BOG2X1),WORK(KT3BOG2X),ISINT2L1R,1)
c     CALL CCFOP_SORT(WORK(KT3BOL2X1),WORK(KT3BOL2X),ISINT2L1R,1)
c
      IF ( (LISTL(1:3).EQ.'L1 ') .OR. (LISTL(1:3).EQ.'M1 ') 
     *                           .OR. (LISTL(1:3).EQ.'N2 ') ) THEN
         CALL INTOCC_T3BARX_JK(LUTOC,FNTOC,ISYMOP,
     *                      WORK(KLAMHL1R),ISYML1R,ISINT2L1R,
     *                      DUMMY,DUMMY,.TRUE.,
     *                      WORK(KT3BOG2X),WORK(KT3BOL2X),
     *                      WORK(KEND1),LWRK1)
      END IF
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
      CALL INTVIR_T3B0_JK(2,WORK(KXGADCK),WORK(KXLADCK),ISYM0,XLAMDP0,
     *                    ISYM0,
     *                         LU3VI,FN3VI,LU3FOP,FN3FOP,
     *                         WORK(KEND1),LWRK1)
C
C----------------------------------------------
C     Get virtual integrals for t3b0 multipliers
C----------------------------------------------
C
      IF ( (LISTL(1:3).EQ.'L1 ') .OR. (LISTL(1:3).EQ.'M1 ') 
     *                           .OR. (LISTL(1:3).EQ.'N2 ') ) THEN
         CALL INTVIR_T3BX_JK(WORK(KXGADCKX),WORK(KXLADCKX),ISINT2L1R,
     *                       WORK(KLAMPL1R),ISYML1R,
     *                       LU3VI,FN3VI,LU3FOP,FN3FOP,
     *                       WORK(KEND1),LWRK1)
      END IF
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
               KEND2   = KTETAXKL + NMAABCI(ISTETAXKL)
               LWRK2  = LWORK  - KEND2

               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Memory available : ',LWORK
                  WRITE(LUPRI,*) 'Memory needed    : ',KEND2
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
c    *                           work(kx3am),1)
C
                  IF (IPRINT .GT. 55) THEN
                    WRITE(LUPRI,*)'ISYML,L,ISYMK,K ', ISYML,L,ISYMK,K
                    XNORMVAL = DDOT(NMAABCI(ISYT30KL),WORK(KT30KL),1,
     *                              WORK(KT30KL),1)
                    WRITE(LUPRI,*)'NORM OF KT30KL IN CC3_ADENOCC ', 
     *                             XNORMVAL
                  END IF
C
C-------------------------------------------
C                 Get T3BAR0^KL multipliers
C-------------------------------------------
C
                  CALL DZERO(WORK(KW3BXKL),NMAABCI(ISYW3BXKL))

                  IF  (LISTL(1:3).EQ.'L1 ') THEN
                    CALL DZERO(WORK(KT3B0KL),NMAABCI(ISYT3B0KL))
C
                    CALL GET_T3B0_JK_O(WORK(KT3B0KL),ISYT3B0KL,
     *                             WORK(KL2TP),ISYML0,
     *                             WORK(KT3BOL2),WORK(KT3BOG2),ISYM0,
     *                             ISYML,L,ISYMK,K,
     *                             WORK(KEND2),LWRK2)
                    CALL GET_T3B0_JK_V(WORK(KT3B0KL),ISYT3B0KL,
     *                                 WORK(KL2TP),ISYML0,
     *                                 WORK(KXGADCK),WORK(KXLADCK),
     *                                 ISYM0,ISYML,L,ISYMK,K,
     *                                 WORK(KEND2),LWRK2)
C
                    CALL GET_T3B0_JK_L1F(WORK(KT3B0KL),ISYT3B0KL,
     *                              WORK(KL1AM),ISYML0,
     *                              WORK(KXIAJB),ISYM0,
     *                              WORK(KL2TP),ISYML0,
     *                              WORK(KFCKBA),ISYM0,
     *                              ISYML,L,ISYMK,K)
C
C----------------------------------------------------------------
C                   Divide by orbital energy difference and remove 
C                   forbidden elements
C----------------------------------------------------------------
C
                    CALL T3JK_DIA(WORK(KT3B0KL),ISYT3B0KL,ISYML,L,ISYMK,
     *                            K,WORK(KFOCKD))
                    CALL T3_FORBIDDEN_JK(WORK(KT3B0KL),ISYMT3B,ISYML,L,
     *                                   ISYMK,K)
c
c                c  all sum_pt3_jk(work(kt3b0kl),isyml,l,isymk,k,isyt3b0kl,
c    *                             work(kx3am),7)
C
                    IF (IPRINT .GT. 55) THEN
                     XNORMVAL = DDOT(NMAABCI(ISYT3B0KL),WORK(KT3B0KL),1,
     *                                 WORK(KT3B0KL),1)
                       WRITE(LUPRI,*)'NORM OF KT3B0KL CC3_ADENOCC ',
     *                                XNORMVAL
                    END IF
C
C---------------------------------------------
C                   Get W3BARX^KL multipliers
C---------------------------------------------
C
C           
C                   <L2|[Y,tau3]|HF> + <L3|[Y^,tau3]|HF>
C
                    CALL WBX_JK_ETA(WORK(KT3B0KL),ISYT3B0KL,
     *                             WORK(KFOCKL1),
     *                             ISYML1,WORK(KW3BXKL),ISYW3BXKL,
     *                             WORK(KL2TP),ISYML0,ISYML,L,ISYMK,K,
     *                             WORK(KEND2),LWRK2)
                   END IF
C
C                 <L2|[H^Y,tau3]|HF>
C
                   IF ( (LISTL(1:3).EQ.'L1 ') .OR. (LISTL(1:3).EQ.'M1 ')
     *                 .OR. (LISTL(1:3).EQ.'N2 ') ) THEN
                     CALL WBX_JK_FMAT(WORK(KW3BXKL),ISYW3BXKL,
     *                               WORK(KL2TP),ISYML0,
     *                               WORK(KFOCKL1RCK),ISYFCKL1R,
     *                               WORK(KT3BOL2X),WORK(KT3BOG2X),
     *                               WORK(KXGADCKX),WORK(KXLADCKX),
     *                               ISINT2L1R,
     *                               ISYML,L,ISYMK,K,
     *                               WORK(KEND2),LWRK2)
                   END IF
C
C                 <L2Y|[H^,tau3]|HF>
C
                  CALL WBX_JK_FMAT(WORK(KW3BXKL),ISYW3BXKL,
     *                            WORK(KL2),ISYML1,
     *                            WORK(KFCKBA),ISYM0,
     *                            WORK(KT3BOL2),WORK(KT3BOG2),
     *                            WORK(KXGADCK),WORK(KXLADCK),ISYM0,
     *                            ISYML,L,ISYMK,K,
     *                            WORK(KEND2),LWRK2)
C
C                 <L1Y|[H^,tau3]|HF>
C
                  CALL WBX_JK_L1(WORK(KW3BXKL),ISYW3BXKL,
     *                           WORK(KL1),ISYML1, 
     *                           WORK(KXIAJB),ISYM0, 
     *                           ISYML,L,ISYMK,K) 
C
C--------------------------------------------------------------
C                 Divide by orbital energy difference and remove 
C                 forbidden elements
C--------------------------------------------------------------
C
                  CALL W3JK_DIA(WORK(KW3BXKL),ISYW3BXKL,ISYML,L,ISYMK,K,
     *                          WORK(KFOCKD),-FREQL1)
                  CALL T3_FORBIDDEN_JK(WORK(KW3BXKL),ISYMW3BX,ISYML,L,
     *                                ISYMK,K)
C
                  !To conform with real sign of t3b multipliers
                  !(noddy code definition)
                  CALL DSCAL(NMAABCI(ISYW3BXKL),-ONE,WORK(KW3BXKL),1)

c                 call sum_pt3_jk(work(kw3bxkl),isyml,l,isymk,k,isyw3bxkl,
c    *                           work(kx3am),4)
C
                  IF (IPRINT .GT. 55) THEN
                     XNORMVAL = DDOT(NMAABCI(ISYW3BXKL),WORK(KW3BXKL),1,
     *                               WORK(KW3BXKL),1)
                     WRITE(LUPRI,*)'NORM OF KW3BXKL IN CC3_ADENOCC ', 
     *                              XNORMVAL
                  END IF
C
C---------------------------------
C                 Calculate KTETAXKL
C---------------------------------
C
                  CALL DZERO(WORK(KTETAXKL),NMAABCI(ISTETAXKL))
C
                  CALL TETAX_JK_BC(WORK(KT30KL),ISYT30KL,WORK(KFOCKR1),
     *                             ISYMR1,WORK(KTETAXKL),ISTETAXKL,
     *                             WORK(KEND2),LWRK2)
C
                  CALL DZERO(WORK(KT3B0KL),NMAABCI(ISTETAXKL))
C
                  CALL DCOPY(NMAABCI(ISTETAXKL),WORK(KTETAXKL),1,
     *                       WORK(KT3B0KL),1)
C
                  CALL W3JK_DIA(WORK(KTETAXKL),ISTETAXKL,ISYML,L,ISYMK,
     *                          K,WORK(KFOCKD),FREQR1)
                  CALL T3_FORBIDDEN_JK(WORK(KTETAXKL),ISYMTETAX,ISYML,L,
     *                                  ISYMK,K)
c
c                 call sum_pt3_jk(work(KTETAXKL),isyml,l,isymk,k,
c    *                            ISTETAXKL,
c    *                            work(kx3am),1)
C
                  IF (IPRINT .GT. 55) THEN
                    XNORMVAL = DDOT(NMAABCI(ISTETAXKL),WORK(KTETAXKL),1,
     *                              WORK(KTETAXKL),1)
                    WRITE(LUPRI,*)'NORM OF KTETAXKL IN CC3_ADENOCC ',
     *                             XNORMVAL
                  END IF 
C
                  CALL ADEN_DIJ_JK(DIJ,WORK(KTETAXKL),ISTETAXKL,
     *                             WORK(KW3BXKL),ISYW3BXKL)
c
c                 write(lupri,*)'DIJ after ADEN_DIJ_JK '
c                 call PRINT_MATIJ(DIJ,ISYDEN)
C
                  IF (IPRINT .GT. 55) THEN
                     XNORMVAL = DDOT(NMATIJ(ISYDEN),DIJ,1,DIJ,1)
                     WRITE(LUPRI,*)'NORM OF DIJ AFTER ADEN_DIJ_JK ',
     *                              XNORMVAL
                  END IF

C
                  CALL TETAX_JK_A(WORK(KT30KL),ISYT30KL,WORK(KFOCKR1),
     *                             ISYMR1,WORK(KT3B0KL),ISTETAXKL,
     *                             WORK(KEND2),LWRK2)
C
                  CALL W3JK_DIA(WORK(KT3B0KL),ISTETAXKL,ISYML,L,ISYMK,K,
     *                           WORK(KFOCKD),FREQR1)
C
                  CALL T3_FORBIDDEN_JK(WORK(KT3B0KL),ISYMTETAX,ISYML,L,
     *                                  ISYMK,K)
C
                  IF (IPRINT .GT. 55) THEN
                    XNORMVAL = DDOT(NMAABCI(ISTETAXKL),WORK(KT3B0KL),1,
     *                              WORK(KT3B0KL),1)
                    WRITE(LUPRI,*)'(FINAL) NORM OF KT3B0KL ', XNORMVAL
C
                    XNORMVAL = DDOT(NMAABCI(ISYW3BXKL),WORK(KW3BXKL),1,
     *                              WORK(KW3BXKL),1)
                    WRITE(LUPRI,*)'(FINAL) NORM OF KW3BXKL ', XNORMVAL
                  END IF
C
                  CALL ADEN_DAB_LM(DAB,WORK(KT3B0KL),ISTETAXKL,
     *                             WORK(KW3BXKL),ISYW3BXKL,
     *                             WORK(KEND2),LWRK2)
C
                  IF (IPRINT .GT. 55) THEN
                    XNORMVAL = DDOT(NMATAB(ISYDEN),DAB,1,DAB,1)
                    WRITE(LUPRI,*)'NORM OF DAB AFTER ADEN_DAB_LM ',
     *                             XNORMVAL
                  END IF
C
                  IF (DO_DIA) THEN
                     CALL ADEN_DAI_LM(DIA,WORK(KL2),ISYML1,
     *                          WORK(KT3B0KL),ISTETAXKL,
     *                          ISYML,L,ISYMK,K,
     *                          WORK(KEND2),LWRK2)
C
                     IF (IPRINT .GT. 55) THEN
                       XNORMVAL = DDOT(NT1AM(ISYDEN),DIA,1,DIA,1)
                       WRITE(LUPRI,*)'NORM OF DIA AFTER ADEN_DAI_LM ',
     *                                XNORMVAL
                     END IF
                 END IF
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
C  /* Deck intvir_t30_ij */
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
C  /* Deck intvir_t3bx_jk */
      SUBROUTINE INTVIR_T3BX_JK(XGADCK,XLADCK,ISYINT,XLAMP,ISYLAMP,
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
      INTEGER ISYINT, ISYLAMP, LU3VI, LU3FOP, LWORK
      INTEGER ISYMD, ISYCKA, KINTVI, KTRVI1, KEND1, LWRK1, IOFF
      INTEGER ISYM0,ISYCKA0
C
      CHARACTER*(*) FN3VI, FN3FOP 
C
#if defined (SYS_CRAY)
      REAL XGADCK(*), XLADCK(*), XLAMP(*), WORK(LWORK)
      real xnormval,ddot
#else
      DOUBLE PRECISION XGADCK(*), XLADCK(*), XLAMP(*), WORK(LWORK) 
      double precision xnormval,ddot
#endif
C
      CALL QENTER('INTVIR_T3BX_JK')
C
C***********************************************************'
C     Get  (XGBDCK) : g(kcad) = (kc ! ad) sorted as I(adck) 
C***********************************************************'
C
      ISYM0 = 1
C
      DO ISYMD = 1, NSYM
         ISYCKA0 = MULD2H(ISYM0,ISYMD)
         ISYCKA = MULD2H(ISYINT,ISYMD)
C
         KINTVI = 1
         KTRVI1 = KINTVI + NCKA(ISYCKA0)
         KEND1  = KTRVI1 + NCKATR(ISYCKA)
         LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in INTVIR_T3BX_JK ')
      ENDIF
C
         DO D = 1, NVIR(ISYMD)
C
C     Read virtual integrals (kc | delta D) stored as I^D(ckdelta)
C
            IOFF = ICKAD(ISYCKA0,ISYMD) + NCKA(ISYCKA0)*(D - 1) + 1
            IF (NCKA(ISYCKA0) .GT. 0) THEN
               CALL GETWA2(LU3VI,FN3VI,WORK(KINTVI),IOFF,
     &                        NCKA(ISYCKA0))
            ENDIF
C
C     Transform to           (kc ! ad)      stored as I^D(cka)
C 
c           CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI1),XLAMP,
c    *                       ISYMD,D,ISYMOP,WORK(KEND1),LWRK1)
c
            CALL CCLR_TRVIR(WORK(KINTVI),WORK(KTRVI1),
     *                      XLAMP,ISYLAMP,
     *                      ISYMD,D,ISYMOP,WORK(KEND1),LWRK1)

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
      ISYM0 = 1
C
      DO ISYMD = 1, NSYM
         ISYCKA0 = MULD2H(ISYM0,ISYMD)
         ISYCKA = MULD2H(ISYINT,ISYMD)
C
         KINTVI = 1
         KTRVI1 = KINTVI + NCKA(ISYCKA0)
         KEND1  = KTRVI1 + NCKATR(ISYCKA)
         LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in INTVIR_T3BX_JK ')
      ENDIF
C
         DO D = 1, NVIR(ISYMD)
C
C     Read virtual integrals L(kc | delta D) stored as L^D(ckdelta)
C
            IOFF = ICKAD(ISYCKA0,ISYMD) + NCKA(ISYCKA0)*(D - 1) + 1
            IF (NCKA(ISYCKA0) .GT. 0) THEN
               CALL GETWA2(LU3FOP,FN3FOP,WORK(KINTVI),IOFF,
     &                        NCKA(ISYCKA0))
            ENDIF
C
C     Transform to           L(kc ! ad)      stored as L^D(cka)
C 
c           CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI1),XLAMP,
c    *                       ISYMD,D,ISYMOP,WORK(KEND1),LWRK1)
c
            CALL CCLR_TRVIR(WORK(KINTVI),WORK(KTRVI1),
     *                      XLAMP,ISYLAMP,
     *                      ISYMD,D,ISYMOP,WORK(KEND1),LWRK1)
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
      CALL QEXIT('INTVIR_T3BX_JK')
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
*
*     The memory bug fixed on 27-Oct-2003, Aarhus, FP.
C
      IMPLICIT NONE
#include <ccsdsym.h>
#include <ccorb.h>
#include <priunit.h>
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
      INTEGER KBCAI
      INTEGER KTEMP,KEND2,LWRK2
c
      integer isymci,isyab,isymi,isyabc
C
#if defined (SYS_CRAY)
      REAL T30JK(*),T2TP(*),T3OG2(*),WORK(LWORK)
      REAL ONE
      real xnormval,ddot
#else
      DOUBLE PRECISION T30JK(*),T2TP(*),T3OG2(*),WORK(LWORK)
      DOUBLE PRECISION ONE
      double precision xnormval,ddot
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
*     KT2LBAI = 1
*     KINTCL  = KT2LBAI + NT2SQ(ISYMT2)
*     KCBAI   = KINTCL  + NT1AM(ISYCL)
*     KEND1   = KCBAI   + NMAAOBCI(ISYT30JK)
*     LWRK1   = LWORK - KEND1
C
      KCBAI   = 1
      KEND1   = KCBAI   + NMAAOBCI(ISYT30JK)
      LWRK1   = LWORK - KEND1
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T30_IJ_O (1)')
      END IF
C
      KT2LBAI = KEND1
      KINTCL  = KT2LBAI + NT2SQ(ISYMT2)
      KEND2   = KINTCL  + NT1AM(ISYCL)
      LWRK2   = LWORK - KEND2
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in GET_T30_IJ_O (1x)')
      END IF
C
      CALL DZERO(WORK(KCBAI),NMAAOBCI(ISYT30JK))
C
      CALL SORT_T2_I_ABJ(WORK(KT2LBAI),T2TP,ISYMT2)
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
      IF (NSYM .GT. 1) THEN
C
         KTEMP = KEND1 
         KEND2    = KTEMP + NMAABCI(ISYT30JK)
         LWRK2    = LWORK - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND2
            CALL QUIT('Insufficient space in GET_T30_IJ_O (1a)')
         END IF
C
         CALL DZERO(WORK(KTEMP),NMAABCI(ISYT30JK))
C
         ! Sort from KCBAI(c,bai) to KTEMP(c,b,a,i)
         CALL FA_BCI(WORK(KTEMP),WORK(KCBAI),ISYT30JK,1)
         CALL DCOPY(NMAAOBCI(ISYT30JK),WORK(KTEMP),1,WORK(KCBAI),1)
      END IF

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
*     KT2BAL = 1
*     KINTLCI  = KT2BAL + NMAABI(ISYT2BAL)
*     KBACI   = KINTLCI + NCKI(ISYINTLCI)
*     KEND1 = KBACI + NMAAB_CI(ISYT30JK)
*     LWRK1  = LWORK - KEND1 
C
      KBACI   = 1
      KEND1 = KBACI + NMAAB_CI(ISYT30JK)
      LWRK1  = LWORK - KEND1 
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T30_IJ_O (2)')
      END IF
C
      KT2BAL = KEND1
      KINTLCI  = KT2BAL + NMAABI(ISYT2BAL)
      KEND2   = KINTLCI + NCKI(ISYINTLCI)
      LWRK2  = LWORK - KEND2 
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in GET_T30_IJ_O (2x)')
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
C-----------------------------------------------
C     Sort from KBACI(ba,ci) to KTEMP(b,a,c,i)
C     and copy back to KBACI
C-----------------------------------------------
C
      IF (NSYM .GT. 1) THEN
C
         KTEMP = KEND1
         KEND2    = KTEMP + NMAAB_CI(ISYT30JK)
         LWRK2  = LWORK - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND2
            CALL QUIT('Insufficient space in GET_T30_IJ_O (2a)')
         END IF
C
         CALL DZERO(WORK(KTEMP),NMAAB_CI(ISYT30JK))
C
         CALL FAB_CI(WORK(KTEMP),WORK(KBACI),ISYT30JK,2)
         CALL DCOPY(NMAABCI(ISYT30JK),WORK(KTEMP),1,WORK(KBACI),1)
      END IF
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
*     KT2BCL = 1
*     KINTLAI  = KT2BCL + NMAABI(ISYT2BCL)
*     KBCAI   = KINTLAI + NCKI(ISYINTLAI)
*     KEND1  = KBCAI + NMAABCI(ISYT30JK)
*     LWRK1  = LWORK - KEND1 
C
      KBCAI   = 1
      KEND1  = KBCAI + NMAABCI(ISYT30JK)
      LWRK1  = LWORK - KEND1 
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T30_IJ_O (3)')
      END IF
C
      KT2BCL = KEND1
      KINTLAI  = KT2BCL + NMAABI(ISYT2BCL)
      KEND2   = KINTLAI + NCKI(ISYINTLAI)
      LWRK2  = LWORK - KEND2 
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in GET_T30_IJ_O (3x)')
      END IF
C
      CALL DZERO(WORK(KBCAI),NMAABCI(ISYT30JK))
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
         KOFF3 = KBCAI
     *         + IMAAB_CI(ISYBC,ISYAI) 
C
         NTOTBC = MAX(NMATAB(ISYBC),1)
         NTOTL  = MAX(NRHF(ISYML),1)
C
C  add_occ(3)
C
         CALL DGEMM('N','N',NMATAB(ISYBC),NT1AM(ISYAI),NRHF(ISYML),
     *              ONE,WORK(KOFF1),NTOTBC,WORK(KOFF2),NTOTL,
     *              ONE,WORK(KOFF3),NTOTBC)
C
      END DO ! ISYML
c
C
         CALL FAB_CI(T30JK,WORK(KBCAI),ISYT30JK,2)
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
*     KT2LCAI = 1
*     KINTBL  = KT2LCAI + NT2SQ(ISYMT2)
*     KBCAI   = KINTBL  + NT1AM(ISYBL)
*     KEND1 = KBCAI + NMAABCI(ISYT30JK)
*     LWRK1   = LWORK - KEND1
C
      KBCAI   = 1
      KEND1 = KBCAI + NMAABCI(ISYT30JK)
      LWRK1   = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T30_IJ_O (4)')
      END IF
C
      KT2LCAI = KEND1
      KINTBL  = KT2LCAI + NT2SQ(ISYMT2)
      KEND2   = KINTBL  + NT1AM(ISYBL)
      LWRK2   = LWORK - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in GET_T30_IJ_O (4x)')
      END IF
C
      CALL DZERO(WORK(KBCAI),NMAABCI(ISYT30JK))
C
      CALL SORT_T2_I_ABJ(WORK(KT2LCAI),T2TP,ISYMT2)
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
            KOFF3 = KBCAI
     *            + IMAAOBCI(ISYMB,ISYCAI) 
C
            NTOTB = MAX(NVIR(ISYMB),1)
            NTOTL = MAX(NRHF(ISYML),1)
C
C  add_occ(4)
C
            CALL DGEMM('N','N',NVIR(ISYMB),NMAABI(ISYCAI),NRHF(ISYML),
     *                 ONE,WORK(KOFF1),NTOTB,WORK(KOFF2),NTOTL,
     *                 ONE,WORK(KOFF3),NTOTB)
C
      END DO ! ISYML
      IF (NSYM .GT. 1) THEN
C
         KTEMP = KEND1
         KEND2    = KTEMP + NMAABCI(ISYT30JK)
         LWRK2    = LWORK - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND2
            CALL QUIT('Insufficient space in GET_T30_IJ_O (4a)')
         END IF
C
         ! Sort from KCBAI(b,cai) to KTEMP(b,c,a,i)
         CALL FA_BCI(WORK(KTEMP),WORK(KBCAI),ISYT30JK,1)
         CALL DCOPY(NMAAOBCI(ISYT30JK),WORK(KTEMP),1,WORK(KBCAI),1)
      END IF
C
      DO I = 1,NMAAOBCI(ISYT30JK)
         T30JK(I) = T30JK(I) + WORK(KBCAI+I-1)
      END DO
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
*     KT2CAL = 1
*     KINTLBI  = KT2CAL + NMAABI(ISYT2CAL)
*     KCABI   = KINTLBI + NCKI(ISYINTLBI)
*     KEND1 = KCABI + NMAAB_CI(ISYT30JK)
*     LWRK1  = LWORK - KEND1 
C
      KCABI   = 1
      KEND1 = KCABI + NMAAB_CI(ISYT30JK)
      LWRK1  = LWORK - KEND1 
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T30_IJ_O (5)')
      END IF
C
      KT2CAL = KEND1
      KINTLBI  = KT2CAL + NMAABI(ISYT2CAL)
      KEND2   = KINTLBI + NCKI(ISYINTLBI)
      LWRK2  = LWORK - KEND2 
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in GET_T30_IJ_O (5x)')
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
     *              ONE,WORK(KOFF3),NTOTCA)
C
      END DO ! ISYML
C
      IF (NSYM .GT. 1) THEN
C
         KTEMP = KEND1
         KEND2    = KTEMP + NMAAB_CI(ISYT30JK)
         LWRK2  = LWORK - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND2
            CALL QUIT('Insufficient space in GET_T30_IJ_O (5a)')
         END IF
C
         CALL DZERO(WORK(KTEMP),NMAAB_CI(ISYT30JK))
C
         CALL FAB_CI(WORK(KTEMP),WORK(KCABI),ISYT30JK,2)
         CALL DCOPY(NMAABCI(ISYT30JK),WORK(KTEMP),1,WORK(KCABI),1)
      END IF

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
*     KT2CBL = 1
*     KINTLAI  = KT2CBL + NMAABI(ISYT2CBL)
*     KCBAI   = KINTLAI + NCKI(ISYINTLAI)
*     KEND1   = KCBAI   + NMAAB_CI(ISYT30JK)
*     LWRK1  = LWORK - KEND1 
C
      KCBAI   = 1
      KEND1   = KCBAI   + NMAAB_CI(ISYT30JK)
      LWRK1  = LWORK - KEND1 
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T30_IJ_O (6)')
      END IF
C
      KT2CBL = KEND1
      KINTLAI  = KT2CBL + NMAABI(ISYT2CBL)
      KEND2   = KINTLAI + NCKI(ISYINTLAI)
      LWRK2  = LWORK - KEND2 
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in GET_T30_IJ_O (6x)')
      END IF
C
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
      IF (NSYM .GT. 1) THEN
C
         KTEMP = KEND1
         KEND2    = KTEMP + NMAAB_CI(ISYT30JK)
         LWRK2  = LWORK - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND2
            CALL QUIT('Insufficient space in GET_T30_IJ_O (6a)')
         END IF
C
         CALL DZERO(WORK(KTEMP),NMAAB_CI(ISYT30JK))
C
         CALL FAB_CI(WORK(KTEMP),WORK(KCBAI),ISYT30JK,2)
         CALL DCOPY(NMAABCI(ISYT30JK),WORK(KTEMP),1,WORK(KCBAI),1)
      END IF

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
*
* Fixed for memory problems, 29-Oct-2003, Aarhus, FP.
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
      INTEGER KBCAI,KTEMP,KEND2,LWRK2
      INTEGER KDCAI,KDBAI
      integer isyabc
C
#if defined (SYS_CRAY)
      REAL T30JK(*), T2TP(*), T3VIJG1(*), WORK(LWORK)
      REAL ONE
      real xnormval,ddot
#else
      DOUBLE PRECISION T30JK(*), T2TP(*), T3VIJG1(*), WORK(LWORK) 
      DOUBLE PRECISION ONE
      double precision xnormval,ddot
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
C
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
C
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
C
*     KBD    = 1
*     KBCAI  = KBD + NMATAB(ISYMBD)
*     KDCAI  = KBCAI + NMAABCI(ISYT30JK)
*     KEND1  = KDCAI + NMAABCI(ISYINT)
*     LWRK1  = LWORK - KEND1
C
      KBCAI  = 1
      KEND1  = KBCAI + NMAABCI(ISYT30JK)
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T30_IJ_V (3)')
      END IF
C
      KBD    = KEND1
      KDCAI  = KBD + NMATAB(ISYMBD)
      KEND2  = KDCAI + NMAABCI(ISYINT)
      LWRK2  = LWORK - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in GET_T30_IJ_V (3x)')
      END IF
C
      CALL DZERO(WORK(KBCAI),NMAABCI(ISYT30JK))

      CALL FA_BCI(WORK(KDCAI),T3VIJG1,ISYINT,2)
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
         KOFF2   = KDCAI  + IMAAOBCI(ISYMD,ISYMCAI)
         KOFF3   = KBCAI  + IMAAOBCI(ISYMB,ISYMCAI)
C 
         NTOTD  = MAX(1,NVIR(ISYMD))
         NTOTB  = MAX(1,NVIR(ISYMB))
C
C
C
         CALL DGEMM('N','N',NVIR(ISYMB),NMAABI(ISYMCAI),
     *                    NVIR(ISYMD),ONE,WORK(KOFF1),NTOTB,
     *                    WORK(KOFF2),NTOTD,
     *                    ONE,WORK(KOFF3),NTOTB)
      END DO
C
      IF (NSYM .GT. 1) THEN
C    
         KTEMP = KEND1
         KEND2    = KTEMP + NMAABCI(ISYT30JK)
         LWRK2    = LWORK - KEND2
C    
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND2
            CALL QUIT('Insufficient space in GET_T30_IJ_V (3a)')
         END IF
C
         ! Sort from KCBAI(b,cai) to KTEMP(b,c,a,i)
         CALL FA_BCI(WORK(KTEMP),WORK(KBCAI),ISYT30JK,1)
         CALL DCOPY(NMAABCI(ISYT30JK),WORK(KTEMP),1,WORK(KBCAI),1)
      END IF
C     
C  add_vir(3)
C
      DO I = 1,NMAABCI(ISYT30JK)
         T30JK(I) = T30JK(I) + WORK(KBCAI+I-1)
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
C
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
C
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
C
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
C
*     KCD    = 1
*     KCBAI  = KCD   +  NMATAB(ISYMCD)
*     KDBAI  = KCBAI +  NMAABCI(ISYMCBAI)
*     KEND1  = KDBAI +  NMAABCI(ISYINT)
*     LWRK1  = LWORK -  KEND1
C
      KCBAI  = 1
      KEND1  = KCBAI + NMAABCI(ISYMCBAI)
      LWRK1  = LWORK -  KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T30_IJ_V (6)')
      END IF
C
      KCD    = KEND1
      KDBAI  = KCD   +  NMATAB(ISYMCD)
      KEND2  = KDBAI +  NMAABCI(ISYINT)
      LWRK2  = LWORK -  KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in GET_T30_IJ_V (6x)')
      END IF
C
      CALL DZERO(WORK(KCBAI),NMAABCI(ISYMCBAI))
C
      CALL FA_BCI(WORK(KDBAI),T3VIJG1,ISYINT,2)
C
C t2tp(ckjd)  =   I^KJ(cd)
C
      CALL SORT_T2_AB(WORK(KCD),ISYMK,K,ISYMJ,J,T2TP,ISYMT2)
C
      DO ISYMD = 1,NSYM
         ISYMC   = MULD2H(ISYMCD,ISYMD)
         ISYMBAI = MULD2H(ISYMD,ISYINT)
         KOFF1   = KCD   + IMATAB(ISYMC,ISYMD)
         KOFF2   = KDBAI     + IMAAOBCI(ISYMD,ISYMBAI)
         KOFF3   = KCBAI + IMAAOBCI(ISYMC,ISYMBAI)
C
         NTOTD  = MAX(1,NVIR(ISYMD))
         NTOTC  = MAX(1,NVIR(ISYMC))
C
C  work(cbai) =  I^KJ(cd) * I(dbai) 
         CALL DGEMM('N','N',NVIR(ISYMC),NMAABI(ISYMBAI),
     *                    NVIR(ISYMD),ONE,WORK(KOFF1),NTOTC,
     *                    WORK(KOFF2),NTOTD,
     *                    ONE,WORK(KOFF3),NTOTC)
      END DO
      IF (NSYM .GT. 1) THEN
C    
         KTEMP = KEND1
         KEND2    = KTEMP + NMAABCI(ISYT30JK)
         LWRK2    = LWORK - KEND2
C    
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND2
            CALL QUIT('Insufficient space in GET_T30_IJ_V (6a)')
         END IF
C
         ! Sort from KCBAI(b,cai) to KTEMP(b,c,a,i)
         CALL FA_BCI(WORK(KTEMP),WORK(KCBAI),ISYT30JK,1)
         CALL DCOPY(NMAABCI(ISYT30JK),WORK(KTEMP),1,WORK(KCBAI),1)
      END IF

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
#include <ccsdsym.h>
#include <ccorb.h>
#include <priunit.h>
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
#include <ccsdsym.h>
#include <ccorb.h>
#include <priunit.h>
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
#include <ccsdsym.h>  
#include <ccorb.h>
#include <priunit.h>
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
#include <ccsdsym.h>
#include <ccorb.h>
#include <priunit.h>
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
#include <ccsdsym.h>  
#include <ccorb.h>
#include <priunit.h>
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
      SUBROUTINE FAB_CI(TAB_CI,TABCI,ISYMT,IOPT)
C
C     TAB_CI = TAB_CI + TABCI
C  
C     if (iopt .eq. 1) then 
C        T(a,b,c,i) --> T(ab,ci)
C     else if (iopt .eq. 2) then
C        T(ab,ci)   --> T(a,b,c,i)
C     end if
C
C     F. Pawlowski, Aarhus, Winter 2003
C
      IMPLICIT NONE
#include <ccsdsym.h>
#include <ccorb.h>
#include <priunit.h>
C
      INTEGER ISYMT,ISYMI,ISYABC,ISYMC,ISYAB,ISYCI,ISYMB,ISYMA
      INTEGER NAB,NCI
      INTEGER KOFF1,KOFF2
      INTEGER IOPT
C
#if defined (SYS_CRAY)
      REAL TABCI(*),TAB_CI(*)
#else
      DOUBLE PRECISION TABCI(*),TAB_CI(*)
#endif
C
      CALL QENTER('FAB_CI')
C
      DO  ISYCI = 1,NSYM
        ISYAB = MULD2H(ISYMT,ISYCI)
         DO ISYMI =  1,NSYM
            ISYMC = MULD2H(ISYCI,ISYMI)
            ISYABC = MULD2H(ISYAB,ISYMC)
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

                           IF (IOPT .EQ. 1) THEN
                              ! T(a,b,c,i) --> T(ab,ci)
                              TAB_CI(KOFF2)=TAB_CI(KOFF2) + TABCI(KOFF1)
                           ELSE IF (IOPT .EQ. 2) THEN
                              ! T(ab,ci)   --> T(a,b,c,i)
                              TAB_CI(KOFF1)=TAB_CI(KOFF1) + TABCI(KOFF2)
                           ELSE
                              WRITE(LUPRI,*)'IOPT = ', IOPT
                              CALL QUIT('Illegal option in FAB_CI')
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
      CALL QEXIT('FAB_CI')
C
      RETURN
      END
C  /* Deck fab_ci */
      SUBROUTINE FA_BCI(TABCI,TAPBCI,ISYMT,IOPT)
C
C     if (iopt .eq. 1) then
C          TAPBCI(a,bci) -->  TABCI(a,b,c,i)
C     else if (iopt .eq. 2) then
C          TAPBCI(a,b,c,i)) --> TABCI(a,bci)
C     end if
C
C     F. Pawlowski, Aarhus, Winter 2003
C
      IMPLICIT NONE
#include <ccsdsym.h>
#include <ccorb.h>
#include <priunit.h>
C
      INTEGER IOPT
      INTEGER ISYMT,ISYMI,ISYABC,ISYMC,ISYAB,ISYMB,ISYMA
      INTEGER ISYBC,ISYBCI
      INTEGER NBCI
      INTEGER KOFF1,KOFF2
C
#if defined (SYS_CRAY)
      REAL TABCI(*),TAPBCI(*)
#else
      DOUBLE PRECISION TABCI(*),TAPBCI(*)
#endif
C
      CALL QENTER('FA_BCI')
C
      DO  ISYMA = 1,NSYM
        ISYBCI = MULD2H(ISYMT,ISYMA)
         DO ISYMI =  1,NSYM
            ISYBC = MULD2H(ISYBCI,ISYMI)
            DO ISYMC = 1,NSYM
               ISYMB =  MULD2H(ISYBC,ISYMC)
               ISYAB = MULD2H(ISYMA,ISYMB)
               ISYABC = MULD2H(ISYAB,ISYMC)
               DO I = 1,NRHF(ISYMI)
                  DO C = 1,NVIR(ISYMC)
                     DO B = 1,NVIR(ISYMB)
C
                        NBCI  = ICKASR(ISYBC,ISYMI)
     *                        + NMATAB(ISYBC)*(I-1)
     *                        + IMATAB(ISYMB,ISYMC)
     *                        + NVIR(ISYMB)*(C-1)
     *                        + B
C
                           KOFF1 = IMAABCI(ISYABC,ISYMI)
     *                           + NMAABC(ISYABC)*(I-1)
     *                           + IMAABC(ISYAB,ISYMC)
     *                           + NMATAB(ISYAB)*(C-1)
     *                           + IMATAB(ISYMA,ISYMB)  
     *                           + NVIR(ISYMA)*(B-1) +  1
                           KOFF2 = IMAAOBCI(ISYMA,ISYBCI) 
     *                           + NVIR(ISYMA)*(NBCI-1) + 1

C
                           IF (IOPT .EQ. 1) THEN
                              ! TAPBCI(a,bci) -->  TABCI(a,b,c,i)
                              CALL DCOPY(NVIR(ISYMA),TAPBCI(KOFF2),1,
     *                                   TABCI(KOFF1),1)
                           ELSE IF (IOPT .EQ. 2) THEN
                              ! TAPBCI(a,b,c,i)) --> TABCI(a,bci)
                              CALL DCOPY(NVIR(ISYMA),TAPBCI(KOFF1),1,
     *                                   TABCI(KOFF2),1)
                           ELSE
                              WRITE(LUPRI,*)'IOPT = ', IOPT
                              CALL QUIT('illegal option in FA_BCI ')
                           END IF
C
                     END DO   
                  END DO   
               END DO   
            END DO   
         END DO   
      END DO   
C
      CALL QEXIT('FA_BCI')
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
#include <ccsdsym.h>
#include <ccorb.h>
#include <priunit.h>
C
      INTEGER ISYINT,ISYMI,ISYMK,ISYAKI,ISYMJ,ISYAK,ISYMA
      INTEGER KOFF1,KOFF2
      integer isyaj
C
#if defined (SYS_CRAY)
      REAL XINTAJ(*),XINTAIJK(*)
      real xnormval,ddot
#else
      DOUBLE PRECISION XINTAJ(*),XINTAIJK(*)
      double precision xnormval,ddot
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
C
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
#include <ccsdsym.h>
#include <ccorb.h>
#include <priunit.h>
C
      INTEGER ISYINT,ISYMI,ISYAKI,ISYMJ,ISYAK,ISYMK,ISYMA,ISYJA
      INTEGER KOFF1,KOFF2,NAK
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

                     NAK   = IT1AM(ISYMA,ISYMK)
     *                     + NVIR(ISYMA)*(K-1)
     *                     + A
                     KOFF2 = IMAIAJ(ISYMJ,ISYAK)
     *                     + NRHF(ISYMJ)*(NAK-1)
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
#include <ccsdsym.h>
#include <ccorb.h>
#include <priunit.h>
C
      INTEGER ISYINT,ISYMK,ISYMI,ISYAKI,ISYMJ,ISYAK,ISYMA,ISYAI
      INTEGER KOFF1,KOFF2,NAI
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
            ISYAI = MULD2H(ISYMA,ISYMI)
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
                     NAI   = IT1AM(ISYMA,ISYMI)
     *                     + NVIR(ISYMA)*(I-1)
     *                     + A
                     KOFF2 = IMAIAJ(ISYMJ,ISYAI)
     *                     + NRHF(ISYMJ)*(NAI-1)
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
#include <ccsdsym.h>
#include <ccorb.h>
#include <priunit.h>
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
C  /* Deck sort_t2_i_abj */
      SUBROUTINE SORT_T2_I_ABJ(T2IABJ,T2TP,ISYMT2)
C
C-------------------------------
C     Sort from T2TP(a,i,j,b) to T2(i,abj)
C-------------------------------
C
C     F. Pawlowski, Aarhus, Winter 2003
C
      IMPLICIT NONE
#include <ccsdsym.h>
#include <ccorb.h>
#include <priunit.h>
C
      INTEGER ISYMT2,ISYMB,ISYAIJ,ISYMJ,ISYAI,ISYMI,ISYMA,ISYAB,ISYABJ
      INTEGER KOFF1,KOFF2,NABJ
C
#if defined (SYS_CRAY)
      REAL T2IABJ(*),T2TP(*)
#else
      DOUBLE PRECISION T2IABJ(*),T2TP(*)
#endif
C
      CALL QENTER('SORT_T2_I_ABJ')
C
C     Sort from T2TP(a,i,j,b) to T2(i,abj)
      DO ISYMB = 1,NSYM
         ISYAIJ = MULD2H(ISYMT2,ISYMB)
         DO ISYMJ = 1,NSYM
            ISYAI = MULD2H(ISYAIJ,ISYMJ)
            DO ISYMI = 1,NSYM
               ISYMA = MULD2H(ISYAI,ISYMI)
               ISYAB = MULD2H(ISYMA,ISYMB)
               ISYABJ = MULD2H(ISYAB,ISYMJ)
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
                           NABJ  = ICKASR(ISYAB,ISYMJ)
     *                           + NMATAB(ISYAB)*(J-1)
     *                           + IMATAB(ISYMA,ISYMB)
     *                           + NVIR(ISYMA)*(B-1)
     *                           + A
                           KOFF2 = IMAJBAI(ISYMI,ISYABJ)
     *                           + NRHF(ISYMI)*(NABJ-1)
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
      CALL QEXIT('SORT_T2_I_ABJ')
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
C     IOPT = 7: Sum up the T3 amplitudes only with aibjck_perm
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
      LDEBUG = .false.
C
C
      IF ( (IOPT .GE. 1) .AND. (IOPT .LE. 3) ) THEN
         FACT = ONE
      ELSE IF ((IOPT .EQ. 4) ) THEN
         FACT = HALF
      ELSE IF ((IOPT .EQ. 6) ) THEN
         FACT = SIXTH
      ELSE IF ((IOPT .EQ. 7) ) THEN
         FACT = ONE
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
               IF ( (IOPT .GE. 1) .AND. (IOPT .LE. 7) ) THEN
                 KOFF7 = IMAABCI(ISYMBCA,ISYMI)
     *                 + NMAABC(ISYMBCA)*(I - 1)
     *                 + IMAABC(ISYMBC,ISYMA) + NMATAB(ISYMBC)*(A-1)
     *                 + IMATAB(ISYMB,ISYMC) + NVIR(ISYMB)*(C-1) + B
               ENDIF
C
         IF (LDEBUG) THEN
            IF (ABS(TMAT(KOFF7)) .GT. 1.0d-12) THEN
              IF ( (IOPT .GE. 1) .AND. (IOPT .LE. 7) ) THEN
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
         IF (IOPT .NE. 7) THEN
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
         END IF
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
      INTEGER ISYMABCI,ISYMABC
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
      ISYMJK = MULD2H(ISYMJ,ISYMK)
      ISYMABCI = MULD2H(ISYMIM,ISYMJK)
      DO ISYMI = 1, NSYM
         ISYMABC = MULD2H(ISYMI,ISYMABCI)
         DO ISYMC = 1,NSYM
            ISYMAB = MULD2H(ISYMABC,ISYMC)
            DO ISYMB = 1,NSYM
               ISYMA = MULD2H(ISYMAB,ISYMB)
               IF ( (ISYMA .EQ. ISYMB) .AND. (ISYMB .EQ. ISYMC) ) THEN
                  ISYMAA = MULD2H(ISYMB,ISYMC)
                  DO I = 1, NRHF(ISYMI)
                     DO A = 1, NVIR(ISYMA)
                        KOFF1 = IMAABCI(ISYMA,ISYMI)
     *                        + NMAABC(ISYMA)*(I-1)
     *                        + IMAABC(ISYMAA,ISYMA)
     *                        + NMATAB(ISYMAA)*(A-1)
     *                        + IMATAB(ISYMA,ISYMA)
     *                        + NVIR(ISYMA)*(A-1)
     *                        + A
C
                              TMAT(KOFF1)  = ZERO
C
                     END DO   
                  END DO   
               END IF
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
C
C Fixed for memory problems, 29-Oct-2003, Aarhus, FP.
*
********************************************************************
      IMPLICIT NONE
#include <ccsdsym.h>
#include <ccorb.h>
#include <priunit.h>
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
      INTEGER KBCAI
      INTEGER KTEMP
c
      integer isyabc,isymi,ioff
C
      INTEGER KEND2,LWRK2
c
C
#if defined (SYS_CRAY)
      REAL T3B0JK(*),L2TP(*),T3BOL2(*),T3BOG2(*),WORK(LWORK)
      REAL ONE
      REAL DDOT,XNORMVAL
#else
      DOUBLE PRECISION T3B0JK(*),L2TP(*),T3BOL2(*),T3BOG2(*),WORK(LWORK)
      DOUBLE PRECISION ONE
      DOUBLE PRECISION DDOT,XNORMVAL
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
*     KL2NBAI = 1
*     KINTCN  = KL2NBAI + NT2SQ(ISYML2)
*     KCBAI   = KINTCN  + NT1AM(ISYCN)
*     KTEMP   = KCBAI   + NMAAOBCI(ISYT3B0JK)
*     KEND1   = KTEMP + NMAAOBCI(ISYT3B0JK)
*     LWRK1   = LWORK - KEND1
C
      KCBAI   = 1
      KEND1   = KCBAI   + NMAAOBCI(ISYT3B0JK)
      LWRK1   = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3B0_JK_O (1a)')
      END IF
C
      KL2NBAI = KEND1
      KINTCN  = KL2NBAI + NT2SQ(ISYML2)
      KEND2   = KINTCN  + NT1AM(ISYCN)
      LWRK2   = LWORK   - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in GET_T3B0_JK_O (1x)')
      END IF
C
      CALL DZERO(WORK(KCBAI),NMAAOBCI(ISYT3B0JK))
C
      CALL SORT_T2_I_ABJ(WORK(KL2NBAI),L2TP,ISYML2)
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
      IF (NSYM .GT. 1) THEN
C
         KTEMP   = KEND1
         KEND1   = KTEMP + NMAAOBCI(ISYT3B0JK)
         LWRK1   = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND1
            CALL QUIT('Insufficient space in GET_T3B0_JK_O (1y)')
         END IF
C
         ! Sort from KCBAI(c,bai) to TABCI(c,b,a,i)
         CALL FA_BCI(WORK(KTEMP),WORK(KCBAI),ISYT3B0JK,1)
         CALL DCOPY(NMAAOBCI(ISYT3B0JK),WORK(KTEMP),1,WORK(KCBAI),1)
      END IF
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
*     KL2BAN = 1
*     KINTNCI  = KL2BAN + NMAABI(ISYL2BAN)
*     KBACI   = KINTNCI + NCKI(ISYINTNCI)
*     KTEMP = KBACI + NMAAB_CI(ISYT3B0JK)
*     KEND1 = KTEMP + NMAAB_CI(ISYT3B0JK)
*     LWRK1  = LWORK - KEND1
C
      KBACI   = 1
      KEND1   = KBACI + NMAAB_CI(ISYT3B0JK)
      LWRK1   = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3BAR0_IJ_O (1b)')
      END IF
C
      KL2BAN   = KEND1
      KINTNCI  = KL2BAN  + NMAABI(ISYL2BAN)
      KEND2    = KINTNCI + NCKI(ISYINTNCI)
      LWRK2    = LWORK   - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in GET_T3BAR0_IJ_O (1bx)')
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
      IF (NSYM .GT. 1) THEN
C
         KTEMP = KEND1
         KEND1 = KTEMP + NMAAB_CI(ISYT3B0JK)
         LWRK1  = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND1
            CALL QUIT('Insufficient space in GET_T3BAR0_IJ_O (1by)')
         END IF
C
         CALL DZERO(WORK(KTEMP),NMAAB_CI(ISYT3B0JK))
C
         CALL FAB_CI(WORK(KTEMP),WORK(KBACI),ISYT3B0JK,2)
         CALL DCOPY(NMAABCI(ISYT3B0JK),WORK(KTEMP),1,WORK(KBACI),1)
      END IF
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
*     KL2BAN = 1
*     KINTNCI  = KL2BAN + NMAABI(ISYL2BAN)
*     KBACI   = KINTNCI + NCKI(ISYINTNCI)
*     KTEMP = KBACI + NMAAB_CI(ISYT3B0JK)
*     KEND1 = KTEMP + NMAABCI(ISYT3B0JK)
*     LWRK1  = LWORK - KEND1
C
      KBACI = 1
      KEND1 = KBACI + NMAAB_CI(ISYT3B0JK)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3BAR0_IJ_O (2a)')
      END IF
C
      KL2BAN   = KEND1
      KINTNCI  = KL2BAN  + NMAABI(ISYL2BAN)
      KEND2    = KINTNCI + NCKI(ISYINTNCI)
      LWRK2    = LWORK   - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in GET_T3BAR0_IJ_O (2ax)')
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
      IF (NSYM .GT. 1) THEN
C
         KTEMP = KEND1
         KEND1 = KTEMP + NMAABCI(ISYT3B0JK)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND1
            CALL QUIT('Insufficient space in GET_T3BAR0_IJ_O (2ay)')
         END IF
C
         CALL DZERO(WORK(KTEMP),NMAABCI(ISYT3B0JK))
C
         CALL FAB_CI(WORK(KTEMP),WORK(KBACI),ISYT3B0JK,2)
         CALL DCOPY(NMAABCI(ISYT3B0JK),WORK(KTEMP),1,WORK(KBACI),1)
      END IF
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
*     KL2BAN = 1
*     KINTNCI  = KL2BAN + NMAABI(ISYL2BAN)
*     KBACI   = KINTNCI + NCKI(ISYINTNCI)
*     KTEMP = KBACI + NMAAB_CI(ISYT3B0JK)
*     KEND1 = KTEMP + NMAABCI(ISYT3B0JK)
*     LWRK1  = LWORK - KEND1
C
      KBACI = 1
      KEND1 = KBACI + NMAAB_CI(ISYT3B0JK)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3BAR0_IJ_O (2b)')
      END IF
C
      KL2BAN   = KEND1
      KINTNCI  = KL2BAN  + NMAABI(ISYL2BAN)
      KEND2    = KINTNCI + NCKI(ISYINTNCI)
      LWRK2    = LWORK   - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in GET_T3BAR0_IJ_O (2bx)')
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
      IF (NSYM .GT. 1) THEN
C
         KTEMP = KEND1
         KEND1 = KTEMP + NMAABCI(ISYT3B0JK)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND1
            CALL QUIT('Insufficient space in GET_T3BAR0_IJ_O (2by)')
         END IF
C
         CALL DZERO(WORK(KTEMP),NMAABCI(ISYT3B0JK))
C
         CALL FAB_CI(WORK(KTEMP),WORK(KBACI),ISYT3B0JK,2)
         CALL DCOPY(NMAABCI(ISYT3B0JK),WORK(KTEMP),1,WORK(KBACI),1)
      END IF
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
*     KL2CBN = 1
*     KINTNAI  = KL2CBN + NMAABI(ISYL2CBN)
*     KCBAI   = KINTNAI + NCKI(ISYINTNAI)
*     KTEMP   = KCBAI   + NMAAB_CI(ISYT3B0JK)
*     KEND1 = KTEMP + NMAABCI(ISYT3B0JK)
*     LWRK1  = LWORK - KEND1
C
      KCBAI = 1
      KEND1 = KCBAI + NMAAB_CI(ISYT3B0JK)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3BAR0_IJ_O (3a)')
      END IF
C
      KL2CBN   = KEND1
      KINTNAI  = KL2CBN  + NMAABI(ISYL2CBN)
      KEND2    = KINTNAI + NCKI(ISYINTNAI)
      LWRK2    = LWORK   - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in GET_T3BAR0_IJ_O (3ax)')
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
      END DO ! ISYMN
C
      IF (NSYM .GT. 1) THEN
C
         KTEMP = KEND1
         KEND1 = KTEMP + NMAABCI(ISYT3B0JK)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND1
            CALL QUIT('Insufficient space in GET_T3BAR0_IJ_O (3ay)')
         END IF
C
         CALL DZERO(WORK(KTEMP),NMAABCI(ISYT3B0JK))
C
         CALL FAB_CI(WORK(KTEMP),WORK(KCBAI),ISYT3B0JK,2)
         CALL DCOPY(NMAABCI(ISYT3B0JK),WORK(KTEMP),1,WORK(KCBAI),1)
      END IF
C
C     T30JK(bcai) = T30JK(bcai) + T^JK(cbai)
C
C  bar_occ(3a)
C
      CALL FBACI(T3B0JK,WORK(KCBAI),ISYT3B0JK)
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
*     KL2NBCI = 1
*     KINTAN  = KL2NBCI + NT2SQ(ISYML2)
*     KABCI   = KINTAN  + NT1AM(ISYAN)
*     KTEMP   = KABCI   + NMAAOBCI(ISYT3B0JK)
*     KEND1   = KTEMP   + NMAAOBCI(ISYT3B0JK)
*     LWRK1   = LWORK - KEND1
C
      KABCI = 1
      KEND1 = KABCI + NMAAOBCI(ISYT3B0JK)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3B0_JK_O (3b)')
      END IF
C
      KL2NBCI = KEND1
      KINTAN  = KL2NBCI + NT2SQ(ISYML2)
      KEND2   = KINTAN  + NT1AM(ISYAN)
      LWRK2   = LWORK   - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in GET_T3B0_JK_O (3bx)')
      END IF
C
      CALL DZERO(WORK(KABCI),NMAAOBCI(ISYT3B0JK))
C
      CALL SORT_T2_I_ABJ(WORK(KL2NBCI),L2TP,ISYML2)
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

      IF (NSYM .GT. 1) THEN
C
         KTEMP   = KEND1
         KEND1   = KTEMP   + NMAAOBCI(ISYT3B0JK)
         LWRK1   = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND1
            CALL QUIT('Insufficient space in GET_T3B0_JK_O (3by)')
         END IF
C
         ! Sort from KABCI(a,bci) to KTEMP(a,b,c,i)
         CALL FA_BCI(WORK(KTEMP),WORK(KABCI),ISYT3B0JK,1)
         CALL DCOPY(NMAAOBCI(ISYT3B0JK),WORK(KTEMP),1,WORK(KABCI),1)
      END IF

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
*     KL2NCAI = 1
*     KINTBN  = KL2NCAI + NT2SQ(ISYML2)
*     KBCAI   = KINTBN  + NT1AM(ISYBN)
*     KTEMP   = KBCAI   + NMAABCI(ISYT3B0JK)
*     KEND1   = KTEMP   + NMAABCI(ISYT3B0JK)
*     LWRK1   = LWORK - KEND1
C
      KBCAI = 1
      KEND1 = KBCAI + NMAABCI(ISYT3B0JK)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3B0_JK_O (4a)')
      END IF
C
      KL2NCAI = KEND1
      KINTBN  = KL2NCAI + NT2SQ(ISYML2)
      KEND2   = KINTBN  + NT1AM(ISYBN)
      LWRK2   = LWORK   - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in GET_T3B0_JK_O (4ax)')
      END IF
C
      CALL DZERO(WORK(KBCAI),NMAABCI(ISYT3B0JK))
C
      CALL SORT_T2_I_ABJ(WORK(KL2NCAI),L2TP,ISYML2)
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
            KOFF3 = KBCAI 
     *            + IMAAOBCI(ISYMB,ISYCAI) 
C
            NTOTB = MAX(NVIR(ISYMB),1)
            NTOTN = MAX(NRHF(ISYMN),1)
C
C  bar_occ(4a)
C
            CALL DGEMM('N','N',NVIR(ISYMB),NMAABI(ISYCAI),NRHF(ISYMN),
     *                 -ONE,WORK(KOFF1),NTOTB,WORK(KOFF2),NTOTN,
     *                 ONE,WORK(KOFF3),NTOTB)
C
      END DO ! ISYMN
C
      IF (NSYM .GT. 1) THEN
C
         KTEMP   = KEND1
         KEND1   = KTEMP   + NMAABCI(ISYT3B0JK)
         LWRK1   = LWORK   - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND1
            CALL QUIT('Insufficient space in GET_T3B0_JK_O (4ay)')
         END IF
C
         ! Sort from KBCAI(b,cai) to KTEMP(b,c,a,i)
         CALL FA_BCI(WORK(KTEMP),WORK(KBCAI),ISYT3B0JK,1)
         CALL DCOPY(NMAAOBCI(ISYT3B0JK),WORK(KTEMP),1,WORK(KBCAI),1)
      END IF
C
      DO I = 1,NMAABCI(ISYT3B0JK)
         T3B0JK(I) = T3B0JK(I) + WORK(KBCAI+I-1)
      END DO
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
*     KL2CAN = 1
*     KINTNBI  = KL2CAN + NMAABI(ISYL2CAN)
*     KCABI   = KINTNBI + NCKI(ISYINTNBI)
*     KTEMP = KCABI + NMAAB_CI(ISYT3B0JK)
*     KEND1 = KTEMP + NMAABCI(ISYT3B0JK)
*     LWRK1  = LWORK - KEND1
C
      KCABI = 1
      KEND1 = KCABI + NMAAB_CI(ISYT3B0JK)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3BAR0_IJ_O (4b)')
      END IF
C
      KL2CAN   = KEND1
      KINTNBI  = KL2CAN  + NMAABI(ISYL2CAN)
      KEND2    = KINTNBI + NCKI(ISYINTNBI)
      LWRK2    = LWORK   - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in GET_T3BAR0_IJ_O (4bx)')
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
      IF (NSYM .GT. 1) THEN
C
         KTEMP = KEND1
         KEND1 = KTEMP + NMAABCI(ISYT3B0JK)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND1
            CALL QUIT('Insufficient space in GET_T3BAR0_IJ_O (4by)')
         END IF
C
         CALL DZERO(WORK(KTEMP),NMAABCI(ISYT3B0JK))
C
         CALL FAB_CI(WORK(KTEMP),WORK(KCABI),ISYT3B0JK,2)
         CALL DCOPY(NMAABCI(ISYT3B0JK),WORK(KTEMP),1,WORK(KCABI),1)
      END IF
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
*     KL2CAN = 1
*     KINTNBI  = KL2CAN + NMAABI(ISYL2CAN)
*     KCABI   = KINTNBI + NCKI(ISYINTNBI)
*     KTEMP = KCABI + NMAAB_CI(ISYT3B0JK)
*     KEND1 = KTEMP + NMAABCI(ISYT3B0JK)
*     LWRK1  = LWORK - KEND1
C
      KCABI = 1
      KEND1 = KCABI + NMAAB_CI(ISYT3B0JK)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3BAR0_IJ_O (5a)')
      END IF
C
      KL2CAN  = KEND1
      KINTNBI = KL2CAN  + NMAABI(ISYL2CAN)
      KEND2   = KINTNBI + NCKI(ISYINTNBI)
      LWRK2   = LWORK   - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in GET_T3BAR0_IJ_O (5ax)')
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
      IF (NSYM .GT. 1) THEN
C
         KTEMP = KEND1
         KEND1 = KTEMP + NMAABCI(ISYT3B0JK)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND1
            CALL QUIT('Insufficient space in GET_T3BAR0_IJ_O (5ay)')
         END IF
C
         CALL DZERO(WORK(KTEMP),NMAABCI(ISYT3B0JK))
C
         CALL FAB_CI(WORK(KTEMP),WORK(KCABI),ISYT3B0JK,2)
         CALL DCOPY(NMAABCI(ISYT3B0JK),WORK(KTEMP),1,WORK(KCABI),1)
      END IF
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
*     KL2CAN = 1
*     KINTNBI  = KL2CAN + NMAABI(ISYL2CAN)
*     KCABI   = KINTNBI + NCKI(ISYINTNBI)
*     KTEMP = KCABI + NMAAB_CI(ISYT3B0JK)
*     KEND1 = KTEMP + NMAABCI(ISYT3B0JK)
*     LWRK1  = LWORK - KEND1
C
      KCABI = 1
      KEND1 = KCABI + NMAAB_CI(ISYT3B0JK)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3BAR0_IJ_O (5b)')
      END IF
C
      KL2CAN  = KEND1
      KINTNBI = KL2CAN  + NMAABI(ISYL2CAN)
      KEND2   = KINTNBI + NCKI(ISYINTNBI)
      LWRK2   = LWORK   - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in GET_T3BAR0_IJ_O (5bx)')
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
      IF (NSYM .GT. 1) THEN
C
         KTEMP = KEND1
         KEND1 = KTEMP + NMAABCI(ISYT3B0JK)
         LWRK1  = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND1
            CALL QUIT('Insufficient space in GET_T3BAR0_IJ_O (5by)')
         END IF
C
         CALL DZERO(WORK(KTEMP),NMAABCI(ISYT3B0JK))
C
         CALL FAB_CI(WORK(KTEMP),WORK(KCABI),ISYT3B0JK,2)
         CALL DCOPY(NMAABCI(ISYT3B0JK),WORK(KTEMP),1,WORK(KCABI),1)
      END IF
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
*     KL2BCN = 1
*     KINTNAI  = KL2BCN + NMAABI(ISYL2BCN)
*     KBCAI   = KINTNAI + NCKI(ISYINTNAI)
*     KTEMP = KBCAI + NMAABCI(ISYT3B0JK)
*     KEND1 = KTEMP + NMAABCI(ISYT3B0JK)
*     LWRK1  = LWORK - KEND1
C
      KBCAI = 1
      KEND1 = KBCAI + NMAABCI(ISYT3B0JK)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3BAR0_IJ_O (6a)')
      END IF
C
      KL2BCN  = KEND1
      KINTNAI = KL2BCN  + NMAABI(ISYL2BCN)
      KEND2   = KINTNAI + NCKI(ISYINTNAI)
      LWRK2   = LWORK   - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in GET_T3BAR0_IJ_O (6ax)')
      END IF
C
      CALL DZERO(WORK(KBCAI),NMAABCI(ISYT3B0JK))
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
         KOFF3 = KBCAI
     *         + IMAAB_CI(ISYBC,ISYAI) 
C
         NTOTBC = MAX(NMATAB(ISYBC),1)
         NTOTN  = MAX(NRHF(ISYMN),1)
C
C  bar_occ(6a)
C
         CALL DGEMM('N','N',NMATAB(ISYBC),NT1AM(ISYAI),NRHF(ISYMN),
     *              -ONE,WORK(KOFF1),NTOTBC,WORK(KOFF2),NTOTN,
     *              ONE,WORK(KOFF3),NTOTBC)
C
      END DO ! ISYMN
C
      IF (NSYM .GT. 1) THEN
C
         KTEMP = KEND1
         KEND1 = KTEMP + NMAABCI(ISYT3B0JK)
         LWRK1  = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND1
            CALL QUIT('Insufficient space in GET_T3BAR0_IJ_O (6ay)')
         END IF
C
         CALL DZERO(WORK(KTEMP),NMAABCI(ISYT3B0JK))
C
         CALL FAB_CI(WORK(KTEMP),WORK(KBCAI),ISYT3B0JK,2)
         CALL DCOPY(NMAABCI(ISYT3B0JK),WORK(KTEMP),1,WORK(KBCAI),1)
      END IF
C
      DO I = 1,NMAABCI(ISYT3B0JK)
         T3B0JK(I) = T3B0JK(I) + WORK(KBCAI+I-1)
      END DO
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
*     KL2NCBI = 1
*     KINTAN  = KL2NCBI + NT2SQ(ISYML2)
*     KACBI   = KINTAN  + NT1AM(ISYAN)
*     KTEMP   = KACBI   + NMAAOBCI(ISYT3B0JK)
*     KEND1   = KTEMP   + NMAAOBCI(ISYT3B0JK)
*     LWRK1   = LWORK - KEND1
C
      KACBI   = 1
      KEND1   = KACBI + NMAAOBCI(ISYT3B0JK)
      LWRK1   = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3B0_JK_O (6b)')
      END IF
C
      KL2NCBI = KEND1
      KINTAN  = KL2NCBI + NT2SQ(ISYML2)
      KEND2   = KINTAN  + NT1AM(ISYAN)
      LWRK2   = LWORK - KEND2
C
       IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in GET_T3B0_JK_O (6bx)')
      END IF
C
      CALL DZERO(WORK(KACBI),NMAAOBCI(ISYT3B0JK))
C
      CALL SORT_T2_I_ABJ(WORK(KL2NCBI),L2TP,ISYML2)
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
      IF (NSYM .GT. 1) THEN
C
         KTEMP   = KEND1
         KEND1   = KTEMP   + NMAAOBCI(ISYT3B0JK)
         LWRK1   = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND1
            CALL QUIT('Insufficient space in GET_T3B0_JK_O (6by)')
         END IF
C

         CALL FA_BCI(WORK(KTEMP),WORK(KACBI),ISYT3B0JK,1)
         CALL DCOPY(NMAAOBCI(ISYT3B0JK),WORK(KTEMP),1,WORK(KACBI),1)
      END IF
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
C  /* Deck intvir_t3b0_jk */
      SUBROUTINE INTVIR_T3B0_JK(IOPT,XGADCK,XLADCK,ISYINT,XLAMP,ISYMLP,
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
*     IF IOPT = 1 do only g integrals
*     IF IOPT = 2 do both 
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
      INTEGER ISYINT, LU3VI, LU3FOP, LWORK, ISYMLP
      INTEGER ISYMD, ISYCKA, KINTVI, KTRVI1, KEND1, LWRK1, IOFF
      INTEGER ISYCKALP
      INTEGER IOPT
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
      !check the option
      IF ( (IOPT .LT. 1) .OR. (IOPT .GT. 2) ) THEN
         WRITE(LUPRI,*)'IOPT = ', IOPT
         WRITE(LUPRI,*)'Only 1 or 2 are allowed values for IOPT'
         CALL QUIT('Wrong option IOPT in INTVIR_T3B0_JK')
      END IF
C
C***********************************************************'
C     Get  (XGBDCK) : g(kcad) = (kc ! ad) sorted as I(adck) 
C***********************************************************'
C
      DO ISYMD = 1, NSYM
         ISYCKA = MULD2H(ISYINT,ISYMD)
         ISYCKALP = MULD2H(ISYCKA,ISYMLP)
C
         KINTVI = 1
         KTRVI1 = KINTVI + NCKA(ISYCKA)
         KEND1  = KTRVI1 + NCKATR(ISYCKALP)
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
c           CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI1),XLAMP,
c    *                       ISYMD,D,ISYMOP,WORK(KEND1),LWRK1)
C
            CALL CCLR_TRVIR(WORK(KINTVI),WORK(KTRVI1),XLAMP,ISYMLP,
     *                       ISYMD,D,ISYINT,WORK(KEND1),LWRK1)
C
C     Final sort (kc ! ad)  as I(adck) 
C
            CALL SORT_INTVIR_T3B0(XGADCK,WORK(KTRVI1),
     *                           D,ISYMD,ISYCKALP,WORK(KEND1),LWRK1)
C
C     OUTPUT (XGADCK) : g(kcad) = (kc ! ad) sorted as I(adck) 
C
         END DO !  D
      END DO !  ISYMD
C
      IF (IOPT .EQ. 2) THEN
C
C***********************************************************'
C     Get (XLBDCK) : L(kcad) sorted as I(adck) 
C***********************************************************'
C
      DO ISYMD = 1, NSYM
         ISYCKA = MULD2H(ISYINT,ISYMD)
         ISYCKALP = MULD2H(ISYCKA,ISYMLP)
C
         KINTVI = 1
         KTRVI1 = KINTVI + NCKA(ISYCKA)
         KEND1  = KTRVI1 + NCKATR(ISYCKALP)
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
c           CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI1),XLAMP,
c    *                       ISYMD,D,ISYMOP,WORK(KEND1),LWRK1)
            CALL CCLR_TRVIR(WORK(KINTVI),WORK(KTRVI1),XLAMP,ISYMLP,
     *                       ISYMD,D,ISYINT,WORK(KEND1),LWRK1)
C
C     Final sort L(kc ! ad)  as I(adck) 
C
            CALL SORT_INTVIR_T3B0(XLADCK,WORK(KTRVI1),
     *                           D,ISYMD,ISYCKALP,WORK(KEND1),LWRK1)
C
C     OUTPUT (XLADCK) : L(kcad) sorted as I(adck) 
C
         END DO !  D
      END DO !  ISYMD
C
      END IF
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
C  /* Deck get_t3b0_jk_v */
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
*
* Fixed for memory problems, 29-Oct-2003, Aarhus, FP.
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
      INTEGER KDCBI,KDCAI,KDBCI,KDBAI
      INTEGER KTEMP
C
      INTEGER KEND2,LWRK2
C
#if defined (SYS_CRAY)
      REAL T30JK(*), L2TP(*), XGADCK(*), XLADCK(*), WORK(LWORK)
      REAL ONE
      real ddot,xnormval
#else
      DOUBLE PRECISION T30JK(*), L2TP(*), XGADCK(*), XLADCK(*)
      DOUBLE PRECISION WORK(LWORK) 
      DOUBLE PRECISION ONE
      double precision ddot,xnormval
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
C
      KDAI  = 1
      KBCD  = KDAI  + NMAABI(ISYMDAI)
      KEND1 = KBCD  + NMAABC(ISYMBCD)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3B0_JK_V (1L)')
      END IF
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
               KOFF3   =  IMAABCI(ISYMBCA,ISYMI)
     *                + NMAABC(ISYMBCA)*(I-1)
     *                + IMAABC(ISYMBC,ISYMA)
     *                + 1
C  
               NTOTBC = MAX(1,NMATAB(ISYMBC))
               NTOTD  = MAX(1,NVIR(ISYMD))
C
C  bar_vir(1L)
C
               CALL DGEMM('N','N',NMATAB(ISYMBC),NVIR(ISYMA),
     *                    NVIR(ISYMD),ONE,WORK(KOFF1),NTOTBC,
     *                    WORK(KOFF2),NTOTD, 
     *                    ONE,T30JK(KOFF3),NTOTBC) 
            END DO
         END DO
      END DO
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
C
*     KAD    = 1
*     KACBI  = KAD + NMATAB(ISYMAD)
*     KTEMP  = KACBI + NMAABCI(ISYT30JK) 
*     KDCBI   = KTEMP   + NMAAOBCI(ISYT30JK)
*     KEND1  = KDCBI + NMAABCI(ISYINT)
*     LWRK1  = LWORK - KEND1
C
      KACBI  = 1
      KEND1  = KACBI + NMAABCI(ISYT30JK) 
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3B0_JK_V (1g)')
      END IF
C
      KAD    = KEND1
      KDCBI  = KAD + NMATAB(ISYMAD)
      KEND2  = KDCBI + NMAABCI(ISYINT)
      LWRK2  = LWORK - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in GET_T3B0_JK_V (1gx)')
      END IF
C
      CALL DZERO(WORK(KACBI),NMAABCI(ISYT30JK))
C
C  L2TP(aJKd) = I^JK(ad)
C
      CALL SORT_T2_AB(WORK(KAD),ISYMJ,J,ISYMK,K,L2TP,ISYML2)
C
C  Sort integrals from XGADCK(d,c,b,i) to KDCBI(d,cbi)
C
      CALL FA_BCI(WORK(KDCBI),XGADCK,ISYINT,2)
C
C T^JK(bcai) = T^JK(bcai) +  I^JK(ad) * I(dcbi)
C
      DO ISYMD = 1,NSYM
         ISYMA   = MULD2H(ISYMAD,ISYMD)
         ISYMCBI = MULD2H(ISYMD,ISYINT)
         KOFF1   = KAD + IMATAB(ISYMA,ISYMD)
         KOFF2   = KDCBI  + IMAAOBCI(ISYMD,ISYMCBI)
         KOFF3   = KACBI  + IMAAOBCI(ISYMA,ISYMCBI)
C
         NTOTD  = MAX(1,NVIR(ISYMD))
         NTOTA  = MAX(1,NVIR(ISYMA))
C
         CALL DGEMM('N','N',NVIR(ISYMA),NMAABI(ISYMCBI),
     *                    NVIR(ISYMD),-ONE,WORK(KOFF1),NTOTA,
     *                    WORK(KOFF2),NTOTD,
     *                    ONE,WORK(KOFF3),NTOTA)
      END DO
C
      IF (NSYM .GT. 1) THEN
C
      KTEMP  = KEND1
      KEND1  = KTEMP + NMAAOBCI(ISYT30JK)
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3B0_JK_V (1gy)')
      END IF
C
         CALL FA_BCI(WORK(KTEMP),WORK(KACBI),ISYT30JK,1)
         CALL DCOPY(NMAAOBCI(ISYT30JK),WORK(KTEMP),1,WORK(KACBI),1)
      END IF
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
C
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
C
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
*
*     KBD    = 1
*     KDCAI  = KBD + NMATAB(ISYMBD)
*     KBCAI  = KDCAI + NMAABCI(ISYINT)
*     KTEMP  = KBCAI + NMAABCI(ISYT30JK)
*     KEND1   = KTEMP   + NMAAOBCI(ISYT30JK)
*     LWRK1  = LWORK - KEND1
C
      KBCAI  = 1
      KEND1  = KBCAI + NMAABCI(ISYT30JK)
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3B0_JK_V (3L)')
      END IF
C
      KBD    = KEND1
      KDCAI  = KBD + NMATAB(ISYMBD)
      KEND2  = KDCAI + NMAABCI(ISYINT)
      LWRK2  = LWORK - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in GET_T3B0_JK_V (3Lx)')
      END IF
C
      CALL DZERO(WORK(KBCAI),NMAABCI(ISYT30JK))
C
C  Sort integrals from XGADCK(d,c,b,i) to KDCAI(d,cai)
C
      CALL FA_BCI(WORK(KDCAI),XLADCK,ISYINT,2)
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
         KOFF2   = KDCAI  + IMAAOBCI(ISYMD,ISYMCAI)
         KOFF3   = KBCAI  + IMAAOBCI(ISYMB,ISYMCAI)
C 
         NTOTD  = MAX(1,NVIR(ISYMD))
         NTOTB  = MAX(1,NVIR(ISYMB))
C
C
C  bar_vir(3L)
C
         CALL DGEMM('N','N',NVIR(ISYMB),NMAABI(ISYMCAI),
     *                    NVIR(ISYMD),ONE,WORK(KOFF1),NTOTB,
     *                    WORK(KOFF2),NTOTD,
     *                    ONE,WORK(KOFF3),NTOTB)
      END DO
C
      IF (NSYM .GT. 1) THEN
C
         KTEMP  = KEND1
         KEND1  = KTEMP   + NMAAOBCI(ISYT30JK)
         LWRK1  = LWORK   - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND1
            CALL QUIT('Insufficient space in GET_T3B0_JK_V (3Ly)')
         END IF
C
*        CALL DZERO(WORK(KTEMP),NMAABCI(ISYT30JK))
C
         CALL FA_BCI(WORK(KTEMP),WORK(KBCAI),ISYT30JK,1)
         CALL DCOPY(NMAAOBCI(ISYT30JK),WORK(KTEMP),1,WORK(KBCAI),1)
      END IF
C
      DO I = 1,NMAABCI(ISYT30JK)
         T30JK(I) = T30JK(I) + WORK(KBCAI+I-1)
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
C
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
C
      KDAI  = 1
      KBCD  = KDAI  + NMAABI(ISYMDAI)
      KEND1 = KBCD  + NMAABC(ISYMBCD)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3B0_JK_V (4L)')
      END IF
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
               KOFF3   = 1 + IMAABCI(ISYMBCA,ISYMI)
     *                + NMAABC(ISYMBCA)*(I-1)
     *                + IMAABC(ISYMBC,ISYMA)
C 
               NTOTBC = MAX(1,NMATAB(ISYMBC))
               NTOTD  = MAX(1,NVIR(ISYMD))
C
C
C  bar_vir(4L)
               CALL DGEMM('N','N',NMATAB(ISYMBC),NVIR(ISYMA),
     *                    NVIR(ISYMD),ONE,WORK(KOFF1),NTOTBC,
     *                    WORK(KOFF2),NTOTD,
     *                    ONE,T30JK(KOFF3),NTOTBC)
            END DO
         END DO
      END DO
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
C
*     KAD    = 1
*     KABCI  = KAD + NMATAB(ISYMAD)
*     KDBCI  = KABCI + NMAABCI(ISYMABCI)
*     KTEMP  = KDBCI + NMAABCI(ISYINT)
*     KEND1 = KTEMP + NMAABCI(ISYT30JK)
*     LWRK1  = LWORK - KEND1
C
      KABCI = 1
      KEND1 = KABCI + NMAABCI(ISYMABCI)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3B0_JK_V (4g)')
      END IF
C
      KAD    = KEND1
      KDBCI  = KAD + NMATAB(ISYMAD)
      KEND2  = KDBCI + NMAABCI(ISYINT)
      LWRK2  = LWORK - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in GET_T3B0_JK_V (4gx)')
      END IF
C
      CALL DZERO(WORK(KABCI),NMAABCI(ISYMABCI))
C
C  Sort integrals from XGADCK(d,b,c,i) to KDBCI(d,bci)
C
      CALL FA_BCI(WORK(KDBCI),XGADCK,ISYINT,2)
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
         KOFF2   = KDBCI  + IMAAOBCI(ISYMD,ISYMBCI)
         KOFF3   = KABCI  + IMAAOBCI(ISYMA,ISYMBCI)
C
         NTOTD  = MAX(1,NVIR(ISYMD))
         NTOTA  = MAX(1,NVIR(ISYMA))
C
C  work(abci) = I^KJ(ad) * I(dbci)
C
         CALL DGEMM('N','N',NVIR(ISYMA),NMAABI(ISYMBCI),
     *                    NVIR(ISYMD),-ONE,WORK(KOFF1),NTOTA,
     *                    WORK(KOFF2),NTOTD,
     *                    ONE,WORK(KOFF3),NTOTA)
      END DO
C
      IF (NSYM .GT. 1) THEN
C
         KTEMP  = KEND1
         KEND1 = KTEMP + NMAABCI(ISYT30JK)
         LWRK1  = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND1
            CALL QUIT('Insufficient space in GET_T3B0_JK_V (4gy)')
         END IF
C
*        CALL DZERO(WORK(KTEMP),NMAABCI(ISYT30JK))
C
         CALL FA_BCI(WORK(KTEMP),WORK(KABCI),ISYT30JK,1)
         CALL DCOPY(NMAAOBCI(ISYT30JK),WORK(KTEMP),1,WORK(KABCI),1)
      END IF
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
C
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
C
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
C
*     KCD    = 1
*     KCBAI  = KCD   +  NMATAB(ISYMCD)
*     KDBAI  = KCBAI +  NMAABCI(ISYMCBAI)
*     KTEMP  = KDBAI + NMAABCI(ISYINT)
*     KEND1  = KTEMP + NMAAOBCI(ISYT30JK)
*     LWRK1  = LWORK -  KEND1
C
      KCBAI = 1
      KEND1 = KCBAI +  NMAABCI(ISYMCBAI)
      LWRK1  = LWORK -  KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_T3B0_JK_V (6L)')
      END IF
C
      KCD    = KEND1
      KDBAI  = KCD   +  NMATAB(ISYMCD)
      KEND2  = KDBAI + NMAABCI(ISYINT)
      LWRK2  = LWORK -  KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in GET_T3B0_JK_V (6Lx)')
      END IF
C
      CALL DZERO(WORK(KCBAI),NMAABCI(ISYMCBAI))
C
C  Sort integrals from XLADCK(d,b,a,i) to KDBAI(d,bai)
C
      CALL FA_BCI(WORK(KDBAI),XLADCK,ISYINT,2)
C
C L2TP(ckjd)  =   I^KJ(cd)
C
      CALL SORT_T2_AB(WORK(KCD),ISYMK,K,ISYMJ,J,L2TP,ISYML2)
C
      DO ISYMD = 1,NSYM
         ISYMC   = MULD2H(ISYMCD,ISYMD)
         ISYMBAI = MULD2H(ISYMD,ISYINT)
         KOFF1   = KCD   + IMATAB(ISYMC,ISYMD)
         KOFF2   = KDBAI     + IMAAOBCI(ISYMD,ISYMBAI)
         KOFF3   = KCBAI + IMAAOBCI(ISYMC,ISYMBAI)
C
         NTOTD  = MAX(1,NVIR(ISYMD))
         NTOTC  = MAX(1,NVIR(ISYMC))
C
C  work(cbai) =  I^KJ(cd) * I(dbai) 
         CALL DGEMM('N','N',NVIR(ISYMC),NMAABI(ISYMBAI),
     *                    NVIR(ISYMD),ONE,WORK(KOFF1),NTOTC,
     *                    WORK(KOFF2),NTOTD,
     *                    ONE,WORK(KOFF3),NTOTC)
      END DO
C
      IF (NSYM .GT. 1) THEN
C
         KTEMP  = KEND1
         KEND1  = KTEMP + NMAAOBCI(ISYT30JK)
         LWRK1  = LWORK -  KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND1
            CALL QUIT('Insufficient space in GET_T3B0_JK_V (6Ly)')
         END IF
C
*        CALL DZERO(WORK(KTEMP),NMAABCI(ISYT30JK))
C
         CALL FA_BCI(WORK(KTEMP),WORK(KCBAI),ISYT30JK,1)
         CALL DCOPY(NMAAOBCI(ISYT30JK),WORK(KTEMP),1,WORK(KCBAI),1)
      END IF
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
C
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
      real xnormval,ddot
#else
      DOUBLE PRECISION TB0JK(*), L2TP(*), WBXJK(*), XOP(*), WORK(LWORK) 
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
C
      CALL QEXIT('WBX_JK_ETA')
      RETURN
      END
C  /* Deck get_t3b0_jk_l1f */
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
C Fixed memory problem, 29-Oct-2003, Aarhus, FP.
C
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
      INTEGER KDCBI,KDBCI
      INTEGER KTEMP
C
      INTEGER KEND2,LWRK2
C
#if defined (SYS_CRAY)
      REAL T3B0JK(*),T2AM(*),FOCKCK(*),T3BOL2(*),T3BOG2(*)
      REAL XGADCK(*),XLADCK(*)
      REAL WORK(LWORK)
      REAL ONE
      real ddot,xnormval
#else
      DOUBLE PRECISION T3B0JK(*),T2AM(*),FOCKCK(*),T3BOL2(*),T3BOG2(*)
      DOUBLE PRECISION XGADCK(*),XLADCK(*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION ONE
      double precision ddot,xnormval
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
                           IF (MULD2H(ISYMC,ISYMK) .EQ. ISYMFCK) THEN
                             T3B0JK(KOFF1) = T3B0JK(KOFF1)
     *                                     + T2AM(KAIJB)*FOCKCK(KKC)
                           END IF
                           IF (MULD2H(ISYMB,ISYMJ) .EQ. ISYMFCK) THEN
                             T3B0JK(KOFF1) = T3B0JK(KOFF1)
     *                                     + T2AM(KAIKC)*FOCKCK(KJB)
                           END IF
                           IF (MULD2H(ISYMC,ISYMJ) .EQ. ISYMFCK) THEN
                             T3B0JK(KOFF1) = T3B0JK(KOFF1)
     *                                     - T2AM(KAIKB)*FOCKCK(KJC)
                           END IF
                           IF (MULD2H(ISYMB,ISYMK) .EQ. ISYMFCK) THEN
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
*     KL2NBAI = 1
*     KINTCN  = KL2NBAI + NT2SQ(ISYMT2)
*     KCBAI   = KINTCN  + NT1AM(ISYCN)
*     KTEMP   = KCBAI   + NMAAOBCI(ISYT3B0JK)
*     KEND1   = KTEMP   + NMAAOBCI(ISYT3B0JK)
*     LWRK1   = LWORK - KEND1
C
      KCBAI   = 1
      KEND1   = KCBAI   + NMAAOBCI(ISYT3B0JK)
      LWRK1   = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in WBX_JK_FMAT (1L)')
      END IF
C
      KL2NBAI = KEND1
      KINTCN  = KL2NBAI + NT2SQ(ISYMT2)
      KEND2   = KINTCN  + NT1AM(ISYCN)
      LWRK2   = LWORK - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in WBX_JK_FMAT (1Lx)')
      END IF
C
      CALL DZERO(WORK(KCBAI),NMAAOBCI(ISYT3B0JK))
C
      CALL SORT_T2_I_ABJ(WORK(KL2NBAI),T2AM,ISYMT2)
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
      IF (NSYM .GT. 1) THEN
C
         KTEMP   = KEND1
         KEND1   = KTEMP   + NMAAOBCI(ISYT3B0JK)
         LWRK1   = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND1
            CALL QUIT('Insufficient space in WBX_JK_FMAT (1Ly)')
         END IF
C
         CALL FA_BCI(WORK(KTEMP),WORK(KCBAI),ISYT3B0JK,1)
         CALL DCOPY(NMAAOBCI(ISYT3B0JK),WORK(KTEMP),1,WORK(KCBAI),1)
      END IF
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
*     KL2BAN = 1
*     KINTNCI  = KL2BAN + NMAABI(ISYL2BAN)
*     KBACI   = KINTNCI + NCKI(ISYINTNCI)
*     KTEMP = KBACI + NMAAB_CI(ISYT3B0JK)
*     KEND1 = KTEMP + NMAABCI(ISYT3B0JK)
*     LWRK1  = LWORK - KEND1
C
      KBACI   = 1
      KEND1   = KBACI + NMAAB_CI(ISYT3B0JK)
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in WBX_JK_FMAT (1g)')
      END IF
C
      KL2BAN = KEND1
      KINTNCI  = KL2BAN + NMAABI(ISYL2BAN)
      KEND2   = KINTNCI + NCKI(ISYINTNCI)
      LWRK2  = LWORK - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in WBX_JK_FMAT (1gx)')
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
      IF (NSYM .GT. 1) THEN
C
         KTEMP = KEND1
         KEND1 = KTEMP + NMAABCI(ISYT3B0JK)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND1
            CALL QUIT('Insufficient space in WBX_JK_FMAT (1gy)')
         END IF
C
         CALL DZERO(WORK(KTEMP),NMAABCI(ISYT3B0JK))
C
         CALL FAB_CI(WORK(KTEMP),WORK(KBACI),ISYT3B0JK,2)
         CALL DCOPY(NMAABCI(ISYT3B0JK),WORK(KTEMP),1,WORK(KBACI),1)
      END IF
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
*     KL2NCAI = 1
*     KINTBN  = KL2NCAI + NT2SQ(ISYMT2)
*     KBCAI   = KINTBN  + NT1AM(ISYBN)
*     KTEMP   = KBCAI   + NMAAOBCI(ISYT3B0JK)
*     KEND1   = KTEMP   + NMAAOBCI(ISYT3B0JK)
*     LWRK1   = LWORK - KEND1
C
      KBCAI = 1
      KEND1 = KBCAI   + NMAAOBCI(ISYT3B0JK)
      LWRK1   = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in WBX_JK_FMAT (2L)')
      END IF
C
      KL2NCAI = KEND1
      KINTBN  = KL2NCAI + NT2SQ(ISYMT2)
      KEND2   = KINTBN  + NT1AM(ISYBN)
      LWRK2   = LWORK - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in WBX_JK_FMAT (2Lx)')
      END IF
C
      CALL DZERO(WORK(KBCAI),NMAAOBCI(ISYT3B0JK))
C
      CALL SORT_T2_I_ABJ(WORK(KL2NCAI),T2AM,ISYMT2)
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
            KOFF3 = KBCAI 
     *            + IMAAOBCI(ISYMB,ISYCAI) 
C
            NTOTB = MAX(NVIR(ISYMB),1)
            NTOTN = MAX(NRHF(ISYMN),1)
C
C  w3x_occ(2L)
C
            CALL DGEMM('N','N',NVIR(ISYMB),NMAABI(ISYCAI),NRHF(ISYMN),
     *                 -ONE,WORK(KOFF1),NTOTB,WORK(KOFF2),NTOTN,
     *                 ONE,WORK(KOFF3),NTOTB)
C
      END DO ! ISYMN
C
      IF (NSYM .GT. 1) THEN
C
         KTEMP   = KEND1
         KEND1   = KTEMP   + NMAAOBCI(ISYT3B0JK)
         LWRK1   = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND1
            CALL QUIT('Insufficient space in WBX_JK_FMAT (2Ly)')
         END IF
C
         CALL FA_BCI(WORK(KTEMP),WORK(KBCAI),ISYT3B0JK,1)
         CALL DCOPY(NMAAOBCI(ISYT3B0JK),WORK(KTEMP),1,WORK(KBCAI),1)
      END IF
C
      DO I = 1,NMAABCI(ISYT3B0JK)
         T3B0JK(I) = T3B0JK(I) + WORK(KBCAI+I-1)
      END DO
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
*     KL2CAN = 1
*     KINTNBI  = KL2CAN + NMAABI(ISYL2CAN)
*     KCABI   = KINTNBI + NCKI(ISYINTNBI)
*     KTEMP = KCABI + NMAAB_CI(ISYT3B0JK)
*     KEND1 = KTEMP + NMAABCI(ISYT3B0JK)
*     LWRK1  = LWORK - KEND1
C
      KCABI = 1
      KEND1 = KCABI + NMAAB_CI(ISYT3B0JK)
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in WBX_JK_FMAT (2g)')
      END IF
C
      KL2CAN = KEND1
      KINTNBI  = KL2CAN + NMAABI(ISYL2CAN)
      KEND2   = KINTNBI + NCKI(ISYINTNBI)
      LWRK2  = LWORK - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in WBX_JK_FMAT (2gx)')
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
      IF (NSYM .GT. 1) THEN
C
         KTEMP = KEND1
         KEND1 = KTEMP + NMAABCI(ISYT3B0JK)
         LWRK1  = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND1
            CALL QUIT('Insufficient space in WBX_JK_FMAT (2gy)')
         END IF
C
         CALL DZERO(WORK(KTEMP),NMAABCI(ISYT3B0JK))
C
         CALL FAB_CI(WORK(KTEMP),WORK(KCABI),ISYT3B0JK,2)
         CALL DCOPY(NMAABCI(ISYT3B0JK),WORK(KTEMP),1,WORK(KCABI),1)
      END IF
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
C
      KDAI  = 1
      KBCD  = KDAI  + NMAABI(ISYMDAI)
      KEND1 = KBCD  + NMAABC(ISYMBCD)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in WBX_JK_FMAT (1Lv)')
      END IF
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
               KOFF3   = 1 + IMAABCI(ISYMBCA,ISYMI)
     *                + NMAABC(ISYMBCA)*(I-1)
     *                + IMAABC(ISYMBC,ISYMA)
C  
               NTOTBC = MAX(1,NMATAB(ISYMBC))
               NTOTD  = MAX(1,NVIR(ISYMD))
C
C  w3x_vir(1L)
C
               CALL DGEMM('N','N',NMATAB(ISYMBC),NVIR(ISYMA),
     *                    NVIR(ISYMD),ONE,WORK(KOFF1),NTOTBC,
     *                    WORK(KOFF2),NTOTD, 
     *                    ONE,T3B0JK(KOFF3),NTOTBC) 
            END DO
         END DO
      END DO
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
C
*     KAD    = 1
*     KACBI  = KAD + NMATAB(ISYMAD)
*     KDCBI  = KACBI + NMAABCI(ISYMACBI) 
*     KTEMP  = KDCBI + NMAABCI(ISYINT)
*     KEND1  = KTEMP + NMAAOBCI(ISYT3B0JK)
*     LWRK1  = LWORK - KEND1
C
      KACBI = 1
      KEND1 = KACBI + NMAABCI(ISYMACBI)
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in WBX_JK_FMAT (1gv)')
      END IF
C
      KAD    = KEND1
      KDCBI  = KAD + NMATAB(ISYMAD)
      KEND2  = KDCBI + NMAABCI(ISYINT)
      LWRK2  = LWORK - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in WBX_JK_FMAT (1gvx)')
      END IF
C
      CALL DZERO(WORK(KACBI),NMAABCI(ISYMACBI))
C
C  Sort integrals from XGADCK(d,c,b,i) to KDCBI(d,cbi)
C
      CALL FA_BCI(WORK(KDCBI),XGADCK,ISYINT,2)
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
         KOFF2   = KDCBI  + IMAAOBCI(ISYMD,ISYMCBI)
         KOFF3   = KACBI  + IMAAOBCI(ISYMA,ISYMCBI)
C
         NTOTD  = MAX(1,NVIR(ISYMD))
         NTOTA  = MAX(1,NVIR(ISYMA))
C
         CALL DGEMM('N','N',NVIR(ISYMA),NMAABI(ISYMCBI),
     *                    NVIR(ISYMD),-ONE,WORK(KOFF1),NTOTA,
     *                    WORK(KOFF2),NTOTD,
     *                    ONE,WORK(KOFF3),NTOTA)
      END DO
C
      IF (NSYM .GT. 1) THEN
C
         KTEMP  = KEND1
         KEND1  = KTEMP + NMAAOBCI(ISYT3B0JK)
         LWRK1  = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND1
            CALL QUIT('Insufficient space in WBX_JK_FMAT (1gvy)')
         END IF
C
         CALL FA_BCI(WORK(KTEMP),WORK(KACBI),ISYT3B0JK,1)
         CALL DCOPY(NMAAOBCI(ISYT3B0JK),WORK(KTEMP),1,WORK(KACBI),1)
      END IF
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
C
      KDAI  = 1
      KBCD  = KDAI  + NMAABI(ISYMDAI)
      KEND1 = KBCD  + NMAABC(ISYMBCD)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in WBX_JK_FMAT (2Lv)')
      END IF
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
               KOFF3   = 1 + IMAABCI(ISYMBCA,ISYMI)
     *                + NMAABC(ISYMBCA)*(I-1)
     *                + IMAABC(ISYMBC,ISYMA)
C 
               NTOTBC = MAX(1,NMATAB(ISYMBC))
               NTOTD  = MAX(1,NVIR(ISYMD))
C
C  w3x_vir(2L)
C
               CALL DGEMM('N','N',NMATAB(ISYMBC),NVIR(ISYMA),
     *                    NVIR(ISYMD),ONE,WORK(KOFF1),NTOTBC,
     *                    WORK(KOFF2),NTOTD,
     *                    ONE,T3B0JK(KOFF3),NTOTBC)
            END DO
         END DO
      END DO
C
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
C
*     KAD    = 1
*     KABCI  = KAD + NMATAB(ISYMAD)
*     KDBCI  = KABCI + NMAABCI(ISYMABCI)
*     KTEMP  = KDBCI + NMAABCI(ISYINT)
*     KEND1  = KTEMP + NMAAOBCI(ISYT3B0JK)
*     LWRK1  = LWORK - KEND1
C
      KABCI = 1
      KEND1 = KABCI + NMAABCI(ISYMABCI)
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in WBX_JK_FMAT (2gv)')
      END IF
C
      KAD = KEND1
      KDBCI = KAD + NMATAB(ISYMAD)
      KEND2 = KDBCI + NMAABCI(ISYINT)
      LWRK2 = LWORK - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in WBX_JK_FMAT (2gvx)')
      END IF
C
      CALL DZERO(WORK(KABCI),NMAABCI(ISYMABCI))
C
C  Sort integrals from XGADCK(d,b,c,i) to KDBCI(d,bci)
C
      CALL FA_BCI(WORK(KDBCI),XGADCK,ISYINT,2)
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
         KOFF2   = KDBCI  + IMAAOBCI(ISYMD,ISYMBCI)
         KOFF3   = KABCI  + IMAAOBCI(ISYMA,ISYMBCI)
C
         NTOTD  = MAX(1,NVIR(ISYMD))
         NTOTA  = MAX(1,NVIR(ISYMA))
C
C  work(abci) = I^KJ(ad) * I(dbci)
C
         CALL DGEMM('N','N',NVIR(ISYMA),NMAABI(ISYMBCI),
     *                    NVIR(ISYMD),-ONE,WORK(KOFF1),NTOTA,
     *                    WORK(KOFF2),NTOTD,
     *                    ONE,WORK(KOFF3),NTOTA)
      END DO
C
      IF (NSYM .GT. 1) THEN
C
         KTEMP  = KEND1
         KEND1  = KTEMP + NMAAOBCI(ISYT3B0JK)
         LWRK1  = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND1
            CALL QUIT('Insufficient space in WBX_JK_FMAT (2gvy)')
         END IF
C
         CALL FA_BCI(WORK(KTEMP),WORK(KABCI),ISYT3B0JK,1)
         CALL DCOPY(NMAAOBCI(ISYT3B0JK),WORK(KTEMP),1,WORK(KABCI),1)
      END IF
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
C  /* Deck wbx_jk_l1 */
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
*    WMAT^JK(bcai) = WMAT^JK(bcai) + T1(ai)*L(JbKc)             
*                                  - T1(aK)*L(Jbic)
*                                  + T1(ai)*L(KcJb)
*                                  - T1(aJ)*L(Kcib)
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
                           IF (MULD2H(ISYMA,ISYMI) .EQ. ISYMT1) THEN 
                             T3B0JK(KOFF1) = T3B0JK(KOFF1) 
     *                                     + TWO*T1AM(NAI)*XIAJB(NCKBJ)
                           END IF
                           IF (MULD2H(ISYMA,ISYMK) .EQ. ISYMT1) THEN 
                             T3B0JK(KOFF1) = T3B0JK(KOFF1) 
     *                                     - T1AM(NAK)*XIAJB(NCIBJ)
                           END IF
                           IF (MULD2H(ISYMA,ISYMJ) .EQ. ISYMT1) THEN 
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
      real ddot,xnormval
#else
      DOUBLE PRECISION T0JK(*), TETAXJK(*), XOP(*), WORK(LWORK) 
      DOUBLE PRECISION ONE
      double precision ddot,xnormval
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
C
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
      real xnormval,ddot
#else
      DOUBLE PRECISION T0JK(*), TETAXJK(*), XOP(*), WORK(LWORK) 
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
C TETAXJK(bcai) = TETAXJK(bcai) 
C
C             - xop(ad) t0_jk(bcdi)
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
      real xnormval,ddot
#else
      DOUBLE PRECISION DIJ(*),THLM(*),WLM(*)
      DOUBLE PRECISION ONE,HALF
      double precision xnormval,ddot
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
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NMAABC(ISYMDEF),
     *              HALF,THLM(KOFF1),NTOTDEF,WLM(KOFF2),NTOTDEF,
     *              ONE,DIJ(KOFF3),NTOTI)
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
               CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),
     *                    NMATAB(ISYMDE),-HALF,WLM(KOFF1),NTOTDE,
     *                    THLM(KOFF2),NTOTDE,ONE,DAB(KOFF3),NTOTA)
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
               CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),
     *                    NMATAB(ISYMDE),-ONE,WORK(KOFF1),NTOTDE,
     *                    WORK(KOFF2),NTOTDE,ONE,DAB(KOFF3),NTOTA)
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
      INTEGER ISYDE,ISYAI,ISYLM
C
#if defined (SYS_CRAY)
      REAL DAI(*),T2TP(*),TETAXKL(*),WORK(LWORK)
      REAL ONE
      real xnormval,ddot
#else
      DOUBLE PRECISION DAI(*),T2TP(*),TETAXKL(*),WORK(LWORK)
      DOUBLE PRECISION ONE
      double precision xnormval,ddot
#endif
C
      PARAMETER (ONE = 1.0D0)
C
      CALL QENTER('ADEN_DAI_LM')
C
      ISYLM = MULD2H(ISYML,ISYMM)
      ISYDE = MULD2H(ISYMT2,ISYLM)
      ISYAI = MULD2H(ISYDE,ISTETAXKL)
      KT2DE = 1
      KTHDAEI = KT2DE + NMATAB(ISYDE)
      KEND1   = KTHDAEI + NMAABCI(ISTETAXKL)
      LWRK1   = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in ADEN_DAI_LM')
      END IF
C
      CALL DZERO(WORK(KT2DE),NMATAB(ISYDE))
      CALL DZERO(WORK(KTHDAEI),NMAABCI(ISTETAXKL))      
C
      CALL SORT_T2_AB(WORK(KT2DE),ISYML,L,ISYMM,M,T2TP,ISYMT2)
C
      CALL DAXPY(NMAABCI(ISTETAXKL),-ONE,TETAXKL,1,WORK(KTHDAEI),1)
C 
      CALL FACBI(TETAXKL,WORK(KTHDAEI),ISTETAXKL)
C
C d(ia) = d(ia) + t2tp(dLMe) * ( teta^LM(deai) - teta^LM(daei) )
C                    I^LM(de)         I^LM(deai)
      DO ISYMI = 1,NSYM
            DO I = 1,NRHF(ISYMI)
               ISYMDEA = MULD2H(ISTETAXKL,ISYMI)
               ISYMA = MULD2H(ISYAI,ISYMI)
               ISYMDE = MULD2H(ISYMDEA,ISYMA)
C 
               KOFF1 = IMAABCI(ISYMDEA,ISYMI)
     *               + NMAABC(ISYMDEA)*(I-1)
     *               + IMAABC(ISYMDE,ISYMA)
     *               + 1
               KOFF2 = KT2DE 
 
               KOFF3 = IT1AM(ISYMA,ISYMI)
     *               + NVIR(ISYMA)*(I-1)
     *               + 1
 
               NTOTDE = MAX(1,NMATAB(ISYMDE))
C 
               CALL DGEMV('T',NMATAB(ISYMDE),NVIR(ISYMA),-ONE,
     *                    TETAXKL(KOFF1),NTOTDE,WORK(KOFF2),1,
     *                    ONE,DAI(KOFF3),1)
C
         END DO
      END DO
C
      CALL QEXIT('ADEN_DAI_LM')
C
      RETURN
      END
C  /* Deck cc3_adenvir */
      SUBROUTINE CC3_ADENVIR(DIJ,DAB,DO_DIA,DIA,ISYDEN,
     *                   DO_YMMAT,YMMAT,
     *                   LISTL,IDLSTL,LISTR,IDLSTR,
     *                   LUTOC,FNTOC,LU3VI,FN3VI,LU3VI2,FN3VI2,
     *                   LUDKBC3,FNDKBC3,
     *                   LU3FOPX,FN3FOPX,LU3FOP2X,FN3FOP2X,
     *                   LU3FOP,FN3FOP,LU3FOP2,FN3FOP2,
     *                   LUCKJD,FNCKJD,LUDKBC,FNDKBC,LUDELD,FNDELD,
     *                   WORK,LWORK)

*---------------------------------------------------------------------*
*
*    Purpose: compute triples component of eta vector 
*             projected into the singles and doubles space
*             
*    W3BAR^Y = ( <L2|[Y,tau3]|HF> + <L3|[Y^,tau3]|HF>
*
*            +   <L2|[H^Y,tau3]|HF> 
*
*            +    <L1^Y + L2^Y|[H^,tau3]|HF> )/ (w+epsiln(tau3))
*
*
*    Written by Poul Jorgensen and Filip Pawlowski, Fall 2002, Aarhus
*            
*=====================================================================*
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <dummy.h>
#include <iratdef.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccinftap.h>
#include <inftap.h>
#include <cc3t3d.h>
#include <ccl1rsp.h>
#include <ccr1rsp.h>
#include <cclrmrsp.h>
#include <ccexci.h>
#include <ccn2rsp.h>
C
      CHARACTER*10 MODEL
C
      LOGICAL DO_YMMAT,DO_DIA
      LOGICAL LSKIPL1R
C
      INTEGER ISYM0
      PARAMETER(ISYM0 = 1)
      CHARACTER LISTL0*3, LISTL*3,LISTR*3,LISTL1R*3,LABELL1*8,LABELR1*8
      CHARACTER*(*) FNTOC, FN3VI, FNDKBC3, FN3FOPX, FN3FOP2X 
      CHARACTER*(*) FNDKBC,FNDELD,FN3VI2,FN3FOP,FN3FOP2,FNCKJD
C
      CHARACTER*10 FNT3, FNWBMAT
      CHARACTER*14 FNTHETA
      PARAMETER(FNT3 = 'CC3_T3_TMP', FNWBMAT = 'CC3_W3_TMP', 
     *          FNTHETA = 'CC3_THETA3_TMP')
C
      CHARACTER*12 FN3SRTR, FNCKJDR, FNDELDR, FNDKBCR
      PARAMETER(FN3SRTR  = 'CCSDT_FBMAT1', FNCKJDR  = 'CCSDT_FBMAT2',
     *          FNDELDR  = 'CCSDT_FBMAT3', FNDKBCR  = 'CCSDT_FBMAT4')
      INTEGER LU3SRTR, LUCKJDR, LUDELDR, LUDKBCR
C
      INTEGER LUTOC, LU3VI, LUDKBC3, LU3FOPX, LU3FOP2X
      INTEGER LUDKBC,LUDELD,LU3VI2,LU3FOP,LU3FOP2,LUCKJD
      INTEGER LUT3,LUWBMAT,LUTHETA
C
      CHARACTER*12 FNT2Y
      PARAMETER ( FNT2Y = 'CC3_R2TP_TMP' )
      INTEGER LUT2Y
C
      LOGICAL   LOCDBG,LORXL1
      PARAMETER (LOCDBG = .FALSE.)
C
      INTEGER  AIBJCK_PERM
      LOGICAL QUADR
      LOGICAL T2XNET2Y
C
      CHARACTER CDUMMY*1 
      PARAMETER (CDUMMY = ' ')

      INTEGER   ISYDEN,IDLSTL,IDLSTR,LWORK
C
      INTEGER IDLSTL0,IDLSTL1R
      INTEGER ISYML1,ISYML1R,ISYMR1
      INTEGER ISINT1,ISINT2
      INTEGER KLAMP0,KLAMH0,KFOCKD,KFOCK0CK,KT2TP,KL1AM,KL2TP
      INTEGER KEND0,LWRK0
      INTEGER KL1L1,KL2L1,KT1R1,KT2R1,KFOCK0,KFOCKL1,KFOCKR1
      INTEGER KEND1,LWRK1
      INTEGER IOPT
      INTEGER ISINT1R1,ISINT2R1,ISINT2L1R,ISYFCKL1R
      INTEGER KXIAJB,KT3BOG1,KT3BOL1,KT3BOG2,KT3BOL2,KT3OG1,KT3OG2
      INTEGER KLAMPL1R,KLAMHL1R,KW3XOGX1,KFOCKL1RCK,KW3BXOG1
      INTEGER KW3BXOL1,KW3BXOGX1,KW3BXOLX1,KT1L1R
      INTEGER KEND2,LWRK2
      INTEGER LENGTH
      INTEGER ISINT1L1R
      INTEGER ISYMD,ISYCKBD0,ISYCKBDL1R,ISYCKBDR1
      INTEGER KT3VDG1,KT3VDG2,KT3VDG3,KT3BVDL1,KT3BVDL2,KT3BVDL3
      INTEGER KEND3,LWRK3
      INTEGER KT3BVDG1,KT3BVDG2,KT3BVDG3,KW3BXVDG1,KW3BXVDG2
      INTEGER KW3BXVDL1,KW3BXVDL2,KW3BXVDGX1,KW3BXVDGX2,KW3BXVDLX1
      INTEGER KW3BXVDLX2,KW3XVDGX1,KINTVI,KTRVI6
      INTEGER KEND4,LWRK4
      INTEGER IOFF
      INTEGER ISYMB,ISYALJB0,ISYALJD0,ISYALJBL1,ISYALJDL1,ISYMBD
      INTEGER ISCKIJ,ISWBMAT,ISWMAT,ISYCKD,ISYCKDBR1
      INTEGER KSMAT2,KUMAT2,KDIAG,KDIAGWB,KDIAGW,KINDSQ,KINDSQWB
      INTEGER KINDSQW,KINDEXB,KINDEXD,KINDEXBL1,KINDEXDL1,KTMAT
      INTEGER KT3MAT,KW3BMAT,KW3MAT,KWTEMP,KS3MAT,KU3MAT,KS3MAT3
      INTEGER KU3MAT3,KT3VBG1,KT3VBG2,KT3VBG3,KT3BVBG1,KT3BVBG2
      INTEGER KT3BVBG3,KSMAT4,KUMAT4,KT3BVBL1,KT3BVBL2,KT3BVBL3
      INTEGER KW3XVDGX2
      INTEGER KEND5,LWRK5
      INTEGER LENSQ,LENSQWB,LENSQW
      INTEGER ISYML,ISYMDL,ISAIBJ,ISYMJ,ISYMBJ,ISYMAI,ISYAIL
      INTEGER KOFF1,NBJ,IADR
      INTEGER KW3MATD
      INTEGER KDAB0,KDIJ0
      INTEGER KT3VBGX3
      INTEGER KR2TP
      INTEGER ISYML0
      INTEGER ISYALJBL0,ISYALJDL0,KINDEXBL0,KINDEXDL0
C
      INTEGER IR1TAMP
      INTEGER ILSTSYM
C
      integer kx3am
C

#if defined (SYS_CRAY)
      REAL      FREQL1,FREQL1R,FREQR1
      REAL      WORK(LWORK)
      REAL      XNORMVAL
      REAL      DAB(*),DIJ(*),DIA(*)
      REAL      YMMAT(*)
      REAL      DDOT,HALF,ONE,TWO
#else
      DOUBLE PRECISION      FREQL1,FREQL1R,FREQR1
      DOUBLE PRECISION      WORK(LWORK)
      DOUBLE PRECISION      XNORMVAL
      DOUBLE PRECISION      DAB(*),DIJ(*),DIA(*)
      DOUBLE PRECISION      YMMAT(*)
      DOUBLE PRECISION      DDOT,HALF,ONE,TWO
#endif
C
      PARAMETER(HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      CALL QENTER('CC3_ADENVIR')
C
C--------------------------------
C     Open temporary files
C--------------------------------
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
C------------------------------------------------------------
C     some initializations:
C------------------------------------------------------------
C
*     LISTL0 = 'L0 '
*     IDLSTL0 = 0 

      IF (LISTL(1:3).EQ.'L1 ') THEN

         ! get symmetry, frequency and integral label from common blocks
         ! defined in ccl1rsp.h
         ISYML1  = ISYLRZ(IDLSTL)
         FREQL1  = FRQLRZ(IDLSTL)
         LABELL1 = LRZLBL(IDLSTL)
         LORXL1  = LORXLRZ(IDLSTL)

         IF (LORXL1) CALL QUIT('NO ORBITAL RELAX. IN CC3_ADENVIR')

        LISTL1R  = 'R1 '
        IDLSTL1R = IR1TAMP(LABELL1,LORXL1,FREQL1,ISYML1)
        ! get symmetry and frequency from common blocks
        ! defined in ccl1rsp.h
        ISYML1R  = ISYLRT(IDLSTL1R)
        FREQL1R  = FRQLRT(IDLSTL1R)
C
        !LITSL0 corresponding to LISTL
        LISTL0 = 'L0 '
        IDLSTL0 = 0
        ISYML0 = 1
C
        IF (ISYML1 .NE. ISYML1R) THEN
           WRITE(LUPRI,*)'ISYML1: ', ISYML1
           WRITE(LUPRI,*)'ISYML1R: ', ISYML1R
           CALL QUIT('Symmetry mismatch in CC3_ADENVIR')
        END IF
C
        IF (FREQL1R .NE. FREQL1) THEN
           WRITE(LUPRI,*)'FREQL1R: ', FREQL1R
           WRITE(LUPRI,*)'FREQL1: ', FREQL1
           CALL QUIT('Frequency mismatch in CC3_ADENVIR(L1)')
        END IF
C
      ELSE IF (LISTL(1:3).EQ.'M1 ') THEN
        ISYML1 = ILSTSYM(LISTL,IDLSTL)
        FREQL1 = FRQLRM(IDLSTL)
        LABELL1 = '- none -'
C
        ! find corresponding right eigenvector
        LISTL1R = 'RE '
        IDLSTL1R = ILRM(IDLSTL)
        ISYML1R = ISYML1
        FREQL1R = EIGVAL(IDLSTL1R)
C
        !LITSL0 corresponding to LISTL
        LISTL0 = 'L0 '
        IDLSTL0 = 0
        ISYML0 = 1
C
        IF (FREQL1R .NE. FREQL1) THEN
           WRITE(LUPRI,*)'FREQL1R: ', FREQL1R
           WRITE(LUPRI,*)'FREQL1: ', FREQL1
           CALL QUIT('Frequency mismatch in CC3_ADENVIR(M1)')
        END IF
C
      ELSE IF (LISTL(1:3).EQ.'N2 ') THEN
        ISYML1 = ILSTSYM(LISTL,IDLSTL)
        FREQL1 = FRQIN2(IDLSTL) + FRQFN2(IDLSTL)
        LABELL1 = '- none -'
C
        ! find corresponding right eigenvector
        LISTL1R = 'RE '
        IDLSTL1R = IFN2(IDLSTL)
        ISYML1R = ILSTSYM(LISTL1R,IDLSTL1R)
        FREQL1R = FRQFN2(IDLSTL)
C
        !LITSL0 corresponding to LISTL
        LISTL0 = 'LE '
        IDLSTL0 = IIN2(IDLSTL)
        ISYML0 = ILSTSYM(LISTL0,IDLSTL0)
C
      ELSE IF (LISTL(1:3).EQ.'LE ') THEN
        ISYML1 = ILSTSYM(LISTL,IDLSTL)
        FREQL1 = -EIGVAL(IDLSTL)
        LABELL1 = '- none -'
C
        !we don't have any "right" vector entering a right hand side
        LISTL1R = '---'
        IDLSTL1R = -99
C
        !LITSL0 corresponding to LISTL (not used for LE)
        LISTL0 = 'L0 '
        IDLSTL0 = 0
        ISYML0 = 1
C
      ELSE
         CALL QUIT('Unknown left list in CC3_ADENVIR')
      END IF

      IF (LISTR(1:3).EQ.'R1 ') THEN
         ! get symmetry, frequency and integral label for right list 
         ! from common blocks defined in ccr1rsp.h
        ISYMR1  = ISYLRT(IDLSTR)
        FREQR1  = FRQLRT(IDLSTR)
        LABELR1 = LRTLBL(IDLSTR)
c
      ELSE IF (LISTR(1:3).EQ.'RE ') THEN
         ISYMR1 = ILSTSYM(LISTR,IDLSTR)
         FREQR1 = EIGVAL(IDLSTR)
         LABELR1 = '- none -'
      ELSE
         WRITE(LUPRI,*)'LISTR ',LISTR
         CALL QUIT('Unknown right list in CC3_ADENVIR')
      END IF
C
C-------------------------------------------------------
C     initial allocations, orbital energy, fock matrix and T2 and L2 :
C-------------------------------------------------------
C
C     Symmetry of integrals in contraction:
C
      ISINT1 = ISYM0
      ISINT2 = ISYM0
C
      KLAMP0 = 1
      KLAMH0  = KLAMP0  + NLAMDT
      KFOCKD  = KLAMH0  + NLAMDT
      KFOCK0CK  = KFOCKD  + NORBTS
      KT2TP   = KFOCK0CK  + NT1AMX 
      KL1AM   = KT2TP   + NT2SQ(ISYM0)
      KL2TP   = KL1AM   + NT1AM(ISYML0)
      KEND0   = KL2TP   + NT2SQ(ISYML0)
      LWRK0   = LWORK   - KEND0
C
      KL1L1   = KEND0
      KL2L1   = KL1L1   + NT1AM(ISYML1)
      KT1R1   = KL2L1   + NT2SQ(ISYML1)
      KT2R1   = KT1R1   + NT1AM(ISYMR1)
      KFOCK0  = KT2R1   + NT2SQ(ISYMR1)
      KFOCKL1  = KFOCK0    + N2BST(ISYM0)
      KFOCKR1   = KFOCKL1    + N2BST(ISYML1)
      KEND1    = KFOCKR1 + N2BST(ISYMR1)
      LWRK1    = LWORK - KEND1
C
      IF (DO_DIA) THEN
         KDAB0 = KEND1
         KDIJ0 = KDAB0 + NMATAB(ISYML1)
         KEND1 = KDIJ0 + NMATIJ(ISYML1)
         LWRK1 = LWORK - KEND1
C
         CALL DZERO(WORK(KDAB0),NMATAB(ISYML1))
         CALL DZERO(WORK(KDIJ0),NMATIJ(ISYML1))
      END IF
C
C-------------------------------------
C     Read in lamdap and lamdh
C-------------------------------------
C
      CALL GET_LAMBDA0(WORK(KLAMP0),WORK(KLAMH0),WORK(KEND1),LWRK1)
C
C---------------------------------------------------------------------
C     Read zeroth-order AO Fock matrix from file and trasform it to 
C     lambda basis
C---------------------------------------------------------------------
C
      CALL GET_FOCK0(WORK(KFOCK0),WORK(KLAMP0),WORK(KLAMH0),WORK(KEND1),
     *               LWRK1)
C
C---------------------------------------------------------------------
C     Read the matrix the property integrals and trasform it to lambda 
C     basis for L1 list and R1 list
C---------------------------------------------------------------------
C
      IF (LISTL(1:3).EQ.'L1 ') THEN
         CALL GET_FOCKX(WORK(KFOCKL1),LABELL1,IDLSTL,ISYML1,
     *                  WORK(KLAMP0),ISYM0,
     *                  WORK(KLAMH0),ISYM0,WORK(KEND1),LWRK1)
      END IF
C
      IF (LISTR(1:3).EQ.'R1 ') THEN
         CALL GET_FOCKX(WORK(KFOCKR1),LABELR1,IDLSTR,ISYMR1,
     *                  WORK(KLAMP0),ISYM0,
     *                  WORK(KLAMH0),ISYM0,WORK(KEND1),LWRK1)
      END IF
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
C     Read L1 and L2 zeroth-order multipliers 
C-------------------------------------
C
      IOPT = 3
      CALL GET_T1_T2(IOPT,.FALSE.,WORK(KL1AM),WORK(KL2TP),LISTL0,
     *                IDLSTL0,
     *               ISYML0,WORK(KEND1),LWRK1)
C
      IF (LOCDBG) WRITE(LUPRI,*) 'Norm of L2TP ',
     *    DDOT(NT2SQ(ISYML0),WORK(KL2TP),1,WORK(KL2TP),1)

C
C-------------------------------------
C     Read L1L1 and L2L1 multipliers 
C-------------------------------------
C
      IOPT  = 3
      CALL GET_T1_T2(IOPT,.FALSE.,WORK(KL1L1),WORK(KL2L1),LISTL,
     *               IDLSTL,ISYML1,WORK(KEND1),LWRK1)
C
C-------------------------------------
C     Read T1R1 and L2R1 amplitudes 
C-------------------------------------
C
      IOPT  = 3
      CALL GET_T1_T2(IOPT,.TRUE.,WORK(KT1R1),WORK(KT2R1),LISTR,
     *               IDLSTR,ISYMR1,WORK(KEND1),LWRK1)
C
C------------------------------------------------------
C     If calculating M(imfN) 
C     get T2 amplitudes for R1 list in a special manner
C------------------------------------------------------
C     
      IF (DO_YMMAT) THEN
          LUT2Y = -1
          CALL WOPEN2(LUT2Y,FNT2Y,64,0)

          KR2TP = KEND1
          KEND1 = KR2TP + NT2SQ(ISYMR1)
          LWRK1 = LWORK - KEND1
          IF (LWRK1.LT.NT2AM(ISYMR1))
     &      CALL QUIT('Out of memory in CC3_ADENVIR (special T2)')

          IOPT = 2
          CALL CC_RDRSP(LISTR,IDLSTR,ISYMR1,IOPT,MODEL,
     &                  DUMMY,WORK(KEND1))
          CALL CCLR_DIASCL(WORK(KEND1),TWO,ISYMR1)
          CALL CC_T2SQ(WORK(KEND1),WORK(KR2TP),ISYMR1)

          CALL PUTWA2(LUT2Y,FNT2Y,WORK(KR2TP),1,NT2SQ(ISYMR1))
      END IF
C
C----------------------------------------
C     Integrals [H,T1Y] where Y is LISTR
C----------------------------------------
C
      ISINT1R1 = MULD2H(ISINT1,ISYMR1)
      ISINT2R1 = MULD2H(ISINT2,ISYMR1)
C
      CALL CC3_BARINT(WORK(KT1R1),ISYMR1,WORK(KLAMP0),
     *                WORK(KLAMH0),WORK(KEND1),LWRK1,
     *                LU3SRTR,FN3SRTR,LUCKJDR,FNCKJDR)
C
      CALL CC3_SORT1(WORK(KEND1),LWRK1,2,ISINT1R1,LU3SRTR,FN3SRTR,
     *               LUDELDR,FNDELDR,IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *               IDUMMY,CDUMMY)
C
      CALL CC3_SINT(WORK(KLAMH0),WORK(KEND1),LWRK1,ISINT1R1,
     *              LUDELDR,FNDELDR,LUDKBCR,FNDKBCR)
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
      CALL SORT_FOCKCK(WORK(KFOCK0CK),WORK(KFOCK0),ISYM0)
C
C----------------------------------------
C     If we want to sum the T3 amplitudes
C----------------------------------------
C
      if (.false.) then
         kx3am  = kend1
         kend1 = kx3am + nrhft*nrhft*nrhft*nvirt*nvirt*nvirt
         call dzero(work(kx3am),nrhft*nrhft*nrhft*nvirt*nvirt*nvirt)
         lwrk0 = lwork - kend1
         if (lwrk0 .lt. 0) then
            write(lupri,*) 'Memory available : ',lwork
            write(lupri,*) 'Memory needed    : ',kend1
            call quit('Insufficient space (T3) in CC3_ADENVIR')
         END IF
      endif
C
C      write(lupri,*) 'WBMAT after dzero'
C      call print_pt3(work(kx3am),ISYML1,4)
C
C-----------------------------
C     Memory allocation.
C-----------------------------
C
      IF ( (LISTL(1:3).EQ.'L1 ') .OR. (LISTL(1:3).EQ.'M1 ') 
     *                           .OR. (LISTL(1:3).EQ.'N2 ') ) THEN
         ISINT2L1R = MULD2H(ISYML1R,ISINT2)
         ISYFCKL1R = MULD2H(ISYMOP,ISYML1R)
      END IF

      KXIAJB   = KEND1
      KEND1   = KXIAJB  + NT2AM(ISYM0)

      KT3BOG1 = KEND1
      KT3BOL1 = KT3BOG1 + NTRAOC(ISYM0)
      KT3BOG2 = KT3BOL1 + NTRAOC(ISYM0)
      KT3BOL2 = KT3BOG2 + NTRAOC(ISYM0)
      KT3OG1  = KT3BOL2 + NTRAOC(ISYM0)
      KT3OG2 = KT3OG1  + NTRAOC(ISINT2)
      KEND1  = KT3OG2 + NTRAOC(ISINT2)
      LWRK1   = LWORK   - KEND1
C
      KW3XOGX1 = KEND1
      KEND1   = KW3XOGX1 + NTRAOC(ISINT2R1)
C
      IF ( (LISTL(1:3).EQ.'L1 ') .OR. (LISTL(1:3).EQ.'M1 ') 
     *                           .OR. (LISTL(1:3).EQ.'N2 ') ) THEN
         KFOCKL1RCK    = KEND1
         KW3BXOGX1   = KFOCKL1RCK    + NT1AM(ISYFCKL1R)
         KW3BXOLX1   = KW3BXOGX1   + NTRAOC(ISINT2L1R)
         KEND1      = KW3BXOLX1   + NTRAOC(ISINT2L1R)
         LWRK1      = LWORK      - KEND1
      END IF
C
      KW3BXOG1   = KEND1
      KW3BXOL1   = KW3BXOG1   + NTRAOC(ISYM0)
      KEND1      = KW3BXOL1   +  NTRAOC(ISYM0)
      LWRK1      = LWORK      - KEND1
C
      IF ( (LISTL(1:3).EQ.'L1 ') .OR. (LISTL(1:3).EQ.'M1 ') 
     *                           .OR. (LISTL(1:3).EQ.'N2 ') ) THEN
         KLAMPL1R  = KEND1
         KLAMHL1R  = KLAMPL1R  + NLAMDT
         KEND1   = KLAMHL1R  + NLAMDT
         LWRK1   = LWORK   - KEND1
C
         KT1L1R  = KEND1
         KEND2  = KT1L1R + NT1AM(ISYML1R)
         LWRK2   = LWORK  - KEND2
      ELSE
         KEND2 = KEND1
         LWRK2 = LWRK1
      END IF
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in CC3_ADENVIR')
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
C---------------------------------------------------
C     Prepare to construct the occupied integrals...
C---------------------------------------------------
C
C        isint1  - symmetry of integrals in standard H, transformed
C                  with LambdaH_0
C        ISINT1L1R - symmetry of integrals in standard H, transformed
C                  with LambdaH_L1R

      ISINT1  = 1

      IF ( (LISTL(1:3).EQ.'L1 ') .OR. (LISTL(1:3).EQ.'M1 ') 
     *                           .OR. (LISTL(1:3).EQ.'N2 ') ) THEN
         ISINT1L1R = MULD2H(ISINT1,ISYML1R)
C
C--------------------------
C     Get Lambda for right list depended on left LISTL list
C--------------------------
C
         CALL GET_LAMBDAX(WORK(KLAMPL1R),WORK(KLAMHL1R),LISTL1R,
     *                    IDLSTL1R,
     *                    ISYML1R,
     *                    WORK(KLAMP0),WORK(KLAMH0),WORK(KEND2),LWRK2)
C
C------------------------------------------------------------------
C        Calculate the F^L1R matrix (kc elements evaluated and stored 
C        as ck)
C------------------------------------------------------------------
C
         IOPT = 1
         CALL GET_T1_T2(IOPT,.FALSE.,WORK(KT1L1R),DUMMY,LISTL1R,
     *                  IDLSTL1R,
     *                  ISYML1R,WORK(KEND2),LWRK2)
         CALL CC3LR_MFOCK(WORK(KFOCKL1RCK),WORK(KT1L1R),WORK(KXIAJB),
     *                    ISYFCKL1R)
C
      END IF
C
         ! From now on WORK(KEND1) is used again, since we do not need
         ! KT1L1R amplitudes any more...
C
C-----------------------------------------------------------------
C     Construct occupied integrals which are required to calculate    
C     t3bar_0 multipliers                                             
C-----------------------------------------------------------------
C
      CALL INTOCC_T3BAR0(LUTOC,FNTOC,WORK(KLAMH0),ISYM0,WORK(KT3BOG1),
     *                   WORK(KT3BOL1),WORK(KT3BOG2),WORK(KT3BOL2),
     *                   WORK(KEND1),LWRK1)
C
C-----------------------------------------------------------------
C     Construct occupied integrals which are required to calculate    
C     t3_0 amplitudes
C-----------------------------------------------------------------
C
      CALL INTOCC_T30(LUCKJD,FNCKJD,WORK(KLAMP0),ISINT2,WORK(KT3OG1),
     *                WORK(KT3OG2),WORK(KEND1),LWRK1)
C
C-----------------------------------------------------------------
C     Construct occupied integrals which are required to calculate    
C     t3bar_Y multipliers                                             
C-----------------------------------------------------------------
C
       
      IF ( (LISTL(1:3).EQ.'L1 ') .OR. (LISTL(1:3).EQ.'M1 ') 
     *                           .OR. (LISTL(1:3).EQ.'N2 ') ) THEN
       LSKIPL1R = .FALSE.
       CALL INTOCC_T3BARX(LSKIPL1R,
     *                   LUTOC,FNTOC,ISYMOP,WORK(KLAMH0),ISYM0,ISINT1,
     *                   WORK(KLAMHL1R),ISYML1R,ISINT1L1R,
     *                   WORK(KW3BXOG1),
     *                   WORK(KW3BXOL1),WORK(KW3BXOGX1),WORK(KW3BXOLX1),
     *                   WORK(KEND1),LWRK1)
      ELSE IF (LISTL(1:3).EQ.'LE ') THEN
       LSKIPL1R = .TRUE.
       CALL INTOCC_T3BARX(LSKIPL1R,
     *                   LUTOC,FNTOC,ISYMOP,WORK(KLAMH0),ISYM0,ISINT1,
     *                   DUMMY,IDUMMY,IDUMMY,
     *                   WORK(KW3BXOG1),
     *                   WORK(KW3BXOL1),DUMMY,DUMMY,
     *                   WORK(KEND1),LWRK1)
      END IF

C
C------------------------------------------------------------------
C     Read occupied integrals [H,T1Y] where Y is LISTR1 (used in W)
C-----------------------------------------------------------------
C
      CALL INTOCC_T3X(LUCKJDR,FNCKJDR,WORK(KLAMP0),ISINT2R1,
     *                WORK(KW3XOGX1),WORK(KEND1),LWRK1)
C
C---------------------------------------------
C     Open files for Tbar and W intermediates:
C---------------------------------------------
C
      LUT3 = -1
      LUWBMAT = -1
      LUTHETA = -1

      CALL WOPEN2(LUT3,FNT3,64,0)
      CALL WOPEN2(LUWBMAT,FNWBMAT,64,0)
      CALL WOPEN2(LUTHETA,FNTHETA,64,0)
C
C----------------------------
C     Loop over D
C----------------------------
C
      DO ISYMD = 1,NSYM

         ISYCKBD0  = MULD2H(ISYMD,ISYM0)
         ISYCKBDR1  = MULD2H(ISYMD,ISINT2R1)
C
         IF ( (LISTL(1:3).EQ.'L1 ') .OR. (LISTL(1:3).EQ.'M1 ') 
     *                              .OR. (LISTL(1:3).EQ.'N2 ') ) THEN
            ISYCKBDL1R  = MULD2H(ISYMD,ISYML1R)
         END IF
C
         DO D = 1,NVIR(ISYMD)
C
C           ------------------
C           Memory allocation.
C           ------------------
            KT3VDG1  = KEND1
            KT3VDG2  = KT3VDG1  + NCKATR(ISYCKBD0)
            KT3VDG3   = KT3VDG2  + NCKATR(ISYCKBD0)
            KEND1   = KT3VDG3 + NCKATR(ISYCKBD0)
C
            KT3BVDL1  = KEND1
            KT3BVDL2  = KT3BVDL1 + NCKATR(ISYCKBD0)
            KT3BVDL3  = KT3BVDL2 + NCKATR(ISYCKBD0)
            KEND3   = KT3BVDL3 + NCKATR(ISYCKBD0)
            LWRK3   = LWORK  - KEND3
           
            KT3BVDG1 = KEND3
            KT3BVDG2 = KT3BVDG1 + NCKATR(ISYCKBD0)
            KT3BVDG3 = KT3BVDG2 + NCKATR(ISYCKBD0)
            KEND3   = KT3BVDG3 + NCKATR(ISYCKBD0)
            LWRK3   = LWORK  - KEND3

            KW3BXVDG1  = KEND3
            KW3BXVDG2  = KW3BXVDG1  + NCKATR(ISYCKBD0)
            KW3BXVDL1  = KW3BXVDG2  + NCKATR(ISYCKBD0)
            KW3BXVDL2  = KW3BXVDL1  + NCKATR(ISYCKBD0)
            KEND3     = KW3BXVDL2  + NCKATR(ISYCKBD0)
            LWRK3     = LWORK     - KEND3

            IF ( (LISTL(1:3).EQ.'L1 ') .OR. (LISTL(1:3).EQ.'M1 ') 
     *                                 .OR. (LISTL(1:3).EQ.'N2 ') ) THEN
              KW3BXVDGX1  = KEND3
              KW3BXVDGX2  = KW3BXVDGX1  + NCKATR(ISYCKBDL1R)
              KW3BXVDLX1  = KW3BXVDGX2  + NCKATR(ISYCKBDL1R)
              KW3BXVDLX2  = KW3BXVDLX1  + NCKATR(ISYCKBDL1R)
              KEND3     = KW3BXVDLX2  + NCKATR(ISYCKBDL1R)
              LWRK3     = LWORK     - KEND3
            END IF
C
            KW3XVDGX1  = KEND3
            KEND3    = KW3XVDGX1  + NCKATR(ISYCKBDR1)
            LWRK3    = LWORK    - KEND3
C
            KINTVI = KEND3
            KTRVI6 = KINTVI + MAX(NCKA(ISYCKBD0),NCKA(ISYCKBDL1R))
            KEND4  = KTRVI6 + NCKATR(ISYCKBD0) 
            LWRK4  = LWORK  - KEND4
           
            IF (LWRK4 .LT. 0) THEN
               WRITE(LUPRI,*) 'Memory available : ',LWORK
               WRITE(LUPRI,*) 'Memory needed    : ',KEND4
               CALL QUIT('Insufficient space in CC3_ADENVIR')
            END IF
C
C-----------------------------------------------------------------------
C           Construct virtual integrals (for fixed D) which are required 
C           to calculate t3_0 amplitudes
C-----------------------------------------------------------------------
C
            CALL INTVIR_T30_D(LUDKBC,FNDKBC,LUDELD,FNDELD,ISINT2,
     *                        WORK(KT3VDG1),WORK(KT3VDG2),WORK(KT3VDG3),
     *                        WORK(KLAMH0),ISYMD,D,WORK(KEND4),LWRK4)
C
C-----------------------------------------------------------------------
C           Construct virtual integrals (for fixed D) which are required 
C           to calculate t3bar_0 multipliers
C-----------------------------------------------------------------------
C
            CALL INTVIR_T3BAR0_D(LU3FOPX,FN3FOPX,LU3FOP2X,FN3FOP2X,
     *                           LUDKBC3,FNDKBC3,LU3VI,FN3VI,ISYM0,
     *                           WORK(KT3BVDL1),WORK(KT3BVDG1),
     *                           WORK(KT3BVDG2),WORK(KT3BVDL2),
     *                           WORK(KT3BVDG3),WORK(KT3BVDL3),
     *                           WORK(KLAMP0),ISYMD,D,WORK(KEND4),LWRK4)
C
C-----------------------------------------------------------------------
C           Construct virtual integrals (for fixed D) which are required 
C           to calculate t3bar_X multipliers
C-----------------------------------------------------------------------
C
            IF ( (LISTL(1:3).EQ.'L1 ') .OR. (LISTL(1:3).EQ.'M1 ') 
     *                                 .OR. (LISTL(1:3).EQ.'N2 ') ) THEN
             LSKIPL1R = .FALSE.
             CALL INTVIR_T3BARX_D(LSKIPL1R,
     *                            ISYMOP,LU3VI,FN3VI,LU3VI2,FN3VI2,
     *                            LU3FOP,FN3FOP,LU3FOP2,FN3FOP2,
     *                            WORK(KW3BXVDGX1),WORK(KW3BXVDG1),
     *                            WORK(KW3BXVDGX2),WORK(KW3BXVDG2),
     *                            WORK(KW3BXVDLX1),WORK(KW3BXVDL1),
     *                            WORK(KW3BXVDLX2),WORK(KW3BXVDL2),
     *                            WORK(KLAMPL1R),ISYML1R,WORK(KLAMP0),
     *                            ISYM0,ISYMD,D,WORK(KEND4),LWRK4)
            ELSE IF (LISTL(1:3).EQ.'LE ') THEN
             LSKIPL1R = .TRUE.
             CALL INTVIR_T3BARX_D(LSKIPL1R,
     *                            ISYMOP,LU3VI,FN3VI,LU3VI2,FN3VI2,
     *                            LU3FOP,FN3FOP,LU3FOP2,FN3FOP2,
     *                            DUMMY,WORK(KW3BXVDG1),
     *                            DUMMY,WORK(KW3BXVDG2),
     *                            DUMMY,WORK(KW3BXVDL1),
     *                            DUMMY,WORK(KW3BXVDL2),
     *                            DUMMY,IDUMMY,WORK(KLAMP0),
     *                            ISYM0,ISYMD,D,WORK(KEND4),LWRK4)
            END IF
C
C-----------------------------------------------------------------------
C        Read virtual integrals [H,T1Y] where Y is LISTR1 (used in W)
C-----------------------------------------------------------------------
C
            IOFF = ICKBD(ISYCKBDR1,ISYMD) + NCKATR(ISYCKBDR1)*(D - 1) 
     *           + 1
            IF (NCKATR(ISYCKBDR1) .GT. 0) THEN
               CALL GETWA2(LUDKBCR,FNDKBCR,WORK(KW3XVDGX1),IOFF,
     &                     NCKATR(ISYCKBDR1))
            ENDIF
C
            DO ISYMB = 1,NSYM

               ISYALJB0  = MULD2H(ISYMB,ISYM0)
               ISYALJD0 = MULD2H(ISYMD,ISYM0)
               ISYALJBL0  = MULD2H(ISYMB,ISYML0)
               ISYALJDL0 = MULD2H(ISYMD,ISYML0)
               ISYALJBL1  = MULD2H(ISYMB,ISYML1)
               ISYALJDL1 = MULD2H(ISYMD,ISYML1)
               ISYMBD  = MULD2H(ISYMD,ISYMB)
               ISCKIJ  = MULD2H(ISYMBD,ISYM0)
               ISWBMAT  = MULD2H(ISCKIJ,ISYML1)
               ISWMAT  = MULD2H(ISCKIJ,ISYMR1)
               ISYCKD  = MULD2H(ISYM0,ISYMB)
C
               ISYCKDBR1  = MULD2H(ISYMB,ISINT2R1)

C              Can use kend3 since we do not need the integrals anymore.
               KSMAT2     = KEND3
               KUMAT2     = KSMAT2    + NCKIJ(ISCKIJ)
               KDIAG      = KUMAT2    + NCKIJ(ISCKIJ)
               KDIAGWB     = KDIAG     + NCKIJ(ISCKIJ)
               KDIAGW     = KDIAGWB    + NCKIJ(ISWBMAT)
               KINDSQ     = KDIAGW    + NCKIJ(ISWMAT)
               KINDSQWB    = KINDSQ    + (6*NCKIJ(ISCKIJ) - 1)/IRAT + 1
               KINDSQW     = KINDSQWB  + (6*NCKIJ(ISWBMAT) - 1)/IRAT + 1
               KINDEXB     = KINDSQW   + (6*NCKIJ(ISWMAT) - 1)/IRAT + 1
               KINDEXD    = KINDEXB    + (NCKI(ISYALJB0)  - 1)/IRAT + 1
               KINDEXBL1   = KINDEXD   + (NCKI(ISYALJD0) - 1)/IRAT + 1
               KINDEXDL1   = KINDEXBL1 + (NCKI(ISYALJBL1)  - 1)/IRAT + 1
               KTMAT      = KINDEXDL1  + (NCKI(ISYALJDL1) - 1)/IRAT + 1
               KT3MAT     = KTMAT    + MAX(NCKIJ(ISCKIJ),NCKIJ(ISWBMAT))
               KW3BMAT      = KT3MAT    + NCKIJ(ISCKIJ)
               KW3MAT      = KW3BMAT     + NCKIJ(ISWBMAT)
               KWTEMP     = KW3MAT      + NCKIJ(ISWMAT)
               KW3MATD      = KWTEMP     + NCKIJ(ISWMAT)
               KEND4        = KW3MATD    + NCKIJ(ISWMAT)
               LWRK4      = LWORK     - KEND4
C
               KINDEXBL0     = KEND4
               KINDEXDL0    = KINDEXBL0 + (NCKI(ISYALJBL0) - 1)/IRAT + 1
               KEND4   = KINDEXDL0   + (NCKI(ISYALJDL0) - 1)/IRAT + 1
               LWRK4      = LWORK     - KEND4
C
               KS3MAT   = KEND4
               KU3MAT   = KS3MAT  + NCKIJ(ISCKIJ)
               KS3MAT3  = KU3MAT  + NCKIJ(ISCKIJ)
               KU3MAT3  = KS3MAT3 + NCKIJ(ISCKIJ)
               KEND4    = KU3MAT3 + NCKIJ(ISCKIJ)

               KT3VBG1  = KEND4
               KT3VBG2  = KT3VBG1  + NCKATR(ISYCKD)
               KT3VBG3   = KT3VBG2  + NCKATR(ISYCKD)
               KEND4   = KT3VBG3 + NCKATR(ISYCKD)

               KT3BVBG1 = KEND4
               KT3BVBG2 = KT3BVBG1 + NCKATR(ISYCKD)
               KT3BVBG3 = KT3BVBG2 + NCKATR(ISYCKD)
               KEND4   = KT3BVBG3 + NCKATR(ISYCKD)
               LWRK4   = LWORK   - KEND4

               KSMAT4  = KEND4
               KUMAT4  = KSMAT4  + NCKIJ(ISCKIJ)
               KT3BVBL1 = KUMAT4  + NCKIJ(ISCKIJ)
               KT3BVBL2 = KT3BVBL1 + NCKATR(ISYCKD)
               KT3BVBL3 = KT3BVBL2 + NCKATR(ISYCKD)
               KEND4   = KT3BVBL3 + NCKATR(ISYCKD)
               LWRK4   = LWORK   - KEND4
C
               KW3XVDGX2 = KEND4
               KEND4   = KW3XVDGX2 + NCKATR(ISYCKDBR1)
C
               KT3VBGX3 = KEND4
               KEND4    = KT3VBGX3 + NCKATR(ISYCKDBR1)
C
               KINTVI  = KEND4
               KEND5   = KINTVI  + NCKA(ISYCKDBR1)
               LWRK5   = LWORK   - KEND5

               IF (LWRK5 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Memory available : ',LWORK
                  WRITE(LUPRI,*) 'Memory needed    : ',KEND5
                  CALL QUIT('Insufficient space in CC3_ADENVIR')
               END IF
C
C
C              -------------------------------
C              Construct part of the diagonal.
C              -------------------------------
C
               CALL CC3_DIAG(WORK(KDIAG), WORK(KFOCKD),ISCKIJ)
               CALL CC3_DIAG(WORK(KDIAGWB),WORK(KFOCKD),ISWBMAT)
               CALL CC3_DIAG(WORK(KDIAGW),WORK(KFOCKD),ISWMAT)

C
C              -----------------------
C              Construct index arrays.
C              -----------------------
C
               LENSQ  = NCKIJ(ISCKIJ)
               CALL CC3_INDSQ(WORK(KINDSQ),LENSQ,ISCKIJ)
               LENSQWB  = NCKIJ(ISWBMAT)
               CALL CC3_INDSQ(WORK(KINDSQWB),LENSQWB,ISWBMAT)
               LENSQW = NCKIJ(ISWMAT)
               CALL CC3_INDSQ(WORK(KINDSQW),LENSQW,ISWMAT) 
C
               CALL CC3_INDEX(WORK(KINDEXB),ISYALJB0)
               CALL CC3_INDEX(WORK(KINDEXD),ISYALJD0)
               CALL CC3_INDEX(WORK(KINDEXBL0),ISYALJBL0)
               CALL CC3_INDEX(WORK(KINDEXDL0),ISYALJDL0)
               CALL CC3_INDEX(WORK(KINDEXBL1),ISYALJBL1)
               CALL CC3_INDEX(WORK(KINDEXDL1),ISYALJDL1)

               DO B = 1,NVIR(ISYMB)
C
                  CALL DZERO(WORK(KW3BMAT),NCKIJ(ISWBMAT))
                  CALL DZERO(WORK(KW3MAT),NCKIJ(ISWMAT))
                  CALL DZERO(WORK(KW3MATD),NCKIJ(ISWMAT))
C
C-----------------------------------------------------------------------
C           Construct virtual integrals (for fixed B) which are required 
C           to calculate t3_0 amplitudes
C           (the same routine as in d-loop is used)
C-----------------------------------------------------------------------
C
                 CALL INTVIR_T30_D(LUDKBC,FNDKBC,LUDELD,FNDELD,ISINT2,
     *                             WORK(KT3VBG1),WORK(KT3VBG2),
     *                             WORK(KT3VBG3),WORK(KLAMH0),ISYMB,B,
     *                             WORK(KEND5),LWRK5)

C
C-----------------------------------------------------------------------
C           Construct virtual integrals (for fixed B) which are required 
C           to calculate t3bar_0 multipliers
C           (the same routine as in d-loop is used)
C-----------------------------------------------------------------------
C
                  CALL INTVIR_T3BAR0_D(LU3FOPX,FN3FOPX,LU3FOP2X,
     *                                 FN3FOP2X,LUDKBC3,FNDKBC3,
     *                                 LU3VI,FN3VI,ISYM0,WORK(KT3BVBL1),
     *                                 WORK(KT3BVBG1),WORK(KT3BVBG2),
     *                                 WORK(KT3BVBL2),WORK(KT3BVBG3),
     *                                 WORK(KT3BVBL3),WORK(KLAMP0),
     *                                 ISYMB,B,WORK(KEND5),LWRK5)
C
C--------------------------------------------------------------------
C           Read virtual integrals [H,T1Y] where Y is LISTR1 (used in W)
C--------------------------------------------------------------------
C
                  IOFF = ICKBD(ISYCKDBR1,ISYMB) +
     &                   NCKATR(ISYCKDBR1)*(B - 1) + 1
                  IF (NCKATR(ISYCKDBR1) .GT. 0) THEN
                     CALL GETWA2(LUDKBCR,FNDKBCR,WORK(KW3XVDGX2),IOFF,
     &                           NCKATR(ISYCKDBR1))
                  ENDIF
C
                  IOFF = ICKAD(ISYCKDBR1,ISYMB) +
     &                   NCKA(ISYCKDBR1)*(B - 1) + 1
                  IF (NCKA(ISYCKDBR1) .GT. 0) THEN
                     CALL GETWA2(LUDELDR,FNDELDR,WORK(KINTVI),IOFF,
     *                    NCKA(ISYCKDBR1))
                  ENDIF
C
                  CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KT3VBGX3),
     *                             WORK(KLAMH0),ISYMB,B,ISINT2R1,
     *                             WORK(KEND5),LWRK5)



C
C-----------------------------------------------------
C                 Get T3_BD amplitudes (using S and U)
C-----------------------------------------------------
C
                  CALL GET_T30_BD(ISYM0,ISINT2,WORK(KT2TP),ISYM0,
     *                            WORK(KT3MAT),WORK(KFOCKD),WORK(KDIAG),
     *                            WORK(KINDSQ),LENSQ,WORK(KS3MAT),
     *                            WORK(KT3VDG1),WORK(KT3VDG2),
     *                            WORK(KT3OG1),WORK(KINDEXB),
     *                            WORK(KS3MAT3),WORK(KT3VBG1),
     *                            WORK(KT3VBG2),WORK(KINDEXD),
     *                            WORK(KU3MAT),WORK(KT3VDG3),
     *                            WORK(KT3OG2),WORK(KU3MAT3),
     *                            WORK(KT3VBG3),ISYMB,B,ISYMD,D,ISCKIJ,
     *                            WORK(KEND5),LWRK5)
c
c       call sum_pt3(work(KT3MAT),isymb,b,isymd,d,
c    *             ISYM0,work(kx3am),4)
C
C---------------------------------------------------------
C                 Get T3bar_BD multipliers (using S and U)
C---------------------------------------------------------
C
                  IF (LISTL(1:3).EQ.'L1 ') THEN
                    CALL GET_T3BAR0_BD(ISYM0,WORK(KL1AM),ISYML0,
     *                                 WORK(KL2TP),ISYML0,WORK(KTMAT),
     *                                 WORK(KFOCK0CK),WORK(KFOCKD),
     *                                 WORK(KDIAG),WORK(KXIAJB),ISYM0,
     *                                 ISYM0,WORK(KINDSQ),LENSQ,
     *                                 WORK(KSMAT2),WORK(KT3BVDG1),
     *                                 WORK(KT3BVDG2),WORK(KT3BVDL1),
     *                                 WORK(KT3BVDL2),WORK(KT3BOG1),
     *                                 WORK(KT3BOL1),WORK(KINDEXBL0),
     *                                 WORK(KSMAT4),WORK(KT3BVBG1),
     *                                 WORK(KT3BVBG2),WORK(KT3BVBL1),
     *                                 WORK(KT3BVBL2),WORK(KINDEXDL0),
     *                                 WORK(KUMAT2),WORK(KT3BVDG3),
     *                                 WORK(KT3BVDL3),WORK(KT3BOG2),
     *                                 WORK(KT3BOL2),WORK(KUMAT4),
     *                                 WORK(KT3BVBG3),WORK(KT3BVBL3),
     *                                 ISYMB,B,ISYMD,D,ISCKIJ,
     *                                 WORK(KEND5),LWRK5)
C Solely for debugging !!!
*                 IF (DO_YMMAT) THEN
*                    CALL WRITE_T3_DL(LUWBMATX,FNWBMATX,WORK(KTMAT),
*    *                                ISYM0,ISYMD,ISYMB,B)
*                 END IF
C
C
c
c       call sum_pt3(work(KTMAT),isymb,b,isymd,d,
c    *             ISYM0,work(kx3am),4)
C
C----------------------------------------------------
C                 Get T3barX_BD multipliers (using W)
C----------------------------------------------------
C
                   CALL GET_T3BARX_BD(.FALSE.,
     *                              WORK(KTMAT),ISCKIJ,WORK(KFOCKL1),
     *                              ISYML1,WORK(KW3BMAT),ISWBMAT,
     *                              WORK(KL2TP),ISYML0,WORK(KFOCKL1RCK),
     *                              ISYFCKL1R,WORK(KW3BXVDLX2),
     *                              WORK(KW3BXVDLX1),WORK(KW3BXVDGX2),
     *                              WORK(KW3BXVDGX1),WORK(KW3BXOLX1),
     *                              WORK(KW3BXOGX1),ISINT2L1R,
     *                              WORK(KINDEXBL0),WORK(KINDEXDL0),
     *                              WORK(KINDSQWB),LENSQWB,WORK(KL2L1),
     *                              ISYML1,WORK(KFOCK0CK),ISYM0,
     *                              WORK(KW3BXVDL2),WORK(KW3BXVDL1),
     *                              WORK(KW3BXVDG2),WORK(KW3BXVDG1),
     *                              WORK(KW3BXOL1),WORK(KW3BXOG1),
     *                              ISINT2,WORK(KINDEXBL1),
     *                              WORK(KINDEXDL1),WORK(KL1L1),ISYML1,
     *                              WORK(KXIAJB),ISINT1,-FREQL1,
     *                              WORK(KDIAGWB),WORK(KFOCKD),B,ISYMB,
     *                              D,ISYMD,ISYML1,WORK(KEND5),LWRK5)
c       call sum_pt3(work(KW3BMAT),isymb,b,isymd,d,
c    *             ISWBMAT,work(kx3am),4)
                  ELSE IF ( (LISTL(1:3).EQ.'M1 ') 
     *                     .OR. (LISTL(1:3).EQ.'N2 ') ) THEN
                     CALL GET_M3BAR_BD(WORK(KTMAT),WORK(KW3BMAT),
     *                      ISWBMAT,WORK(KL2TP),ISYML0,WORK(KFOCKL1RCK),
     *                      ISYFCKL1R,WORK(KW3BXVDLX2),
     *                      WORK(KW3BXVDLX1),WORK(KW3BXVDGX2),
     *                      WORK(KW3BXVDGX1),WORK(KW3BXOLX1),
     *                      WORK(KW3BXOGX1),ISINT2L1R,
     *                      WORK(KINDEXBL0),WORK(KINDEXDL0),
     *                      WORK(KINDSQWB),LENSQWB,WORK(KL2L1),
     *                      ISYML1,WORK(KFOCK0CK),ISYM0,
     *                      WORK(KW3BXVDL2),WORK(KW3BXVDL1),
     *                      WORK(KW3BXVDG2),WORK(KW3BXVDG1),
     *                      WORK(KW3BXOL1),WORK(KW3BXOG1),
     *                      ISINT2,WORK(KINDEXBL1),
     *                      WORK(KINDEXDL1),WORK(KL1L1),ISYML1,
     *                      WORK(KXIAJB),ISINT1,-FREQL1,
     *                      WORK(KDIAGWB),WORK(KFOCKD),B,ISYMB,
     *                      D,ISYMD,ISYML1,WORK(KEND5),LWRK5)

c
c       call sum_pt3(work(KW3BMAT),isymb,b,isymd,d,
c    *             ISWBMAT,work(kx3am),4)
                  ELSE IF (LISTL(1:3).EQ.'LE ') THEN
                     CALL GET_L3BAR_BD(WORK(KTMAT),WORK(KW3BMAT),
     *                       ISWBMAT,
     *                       WORK(KINDSQWB),LENSQWB,WORK(KL2L1),
     *                       ISYML1,WORK(KFOCK0CK),ISYM0,
     *                       WORK(KW3BXVDL2),WORK(KW3BXVDL1),
     *                       WORK(KW3BXVDG2),WORK(KW3BXVDG1),
     *                       WORK(KW3BXOL1),WORK(KW3BXOG1),
     *                       ISINT2,WORK(KINDEXBL1),
     *                       WORK(KINDEXDL1),WORK(KL1L1),ISYML1,
     *                       WORK(KXIAJB),ISINT1,-FREQL1,
     *                       WORK(KDIAGWB),WORK(KFOCKD),B,ISYMB,
     *                       D,ISYMD,ISYML1,WORK(KEND5),LWRK5)
c
c       call sum_pt3(work(KW3BMAT),isymb,b,isymd,d,
c    *             ISWBMAT,work(kx3am),4)
                  END IF
C
                  !To conform with real sign of t3b multipliers
                  !(noddy code definition)
                  CALL DSCAL(NCKIJ(ISWBMAT),-ONE,WORK(KW3BMAT),1)
c
c       call sum_pt3(work(KW3BMAT),isymb,b,isymd,d,
c    *             ISWBMAT,work(kx3am),4)
C
C------------------------------------------------------
C     Calculate the  term <mu3|[X,T3]|HF> occupied contribution 
C     added in W^BD(KW3MAT)
C------------------------------------------------------
C
                  AIBJCK_PERM = 1
C
                  IF (LISTR(1:3).EQ.'R1 ') THEN
                   CALL WXBD_O(AIBJCK_PERM,WORK(KT3MAT),ISCKIJ,
     *                        WORK(KFOCKR1),ISYMR1,
     *                        WORK(KW3MAT),ISWMAT,WORK(KEND5),LWRK5)
C
C------------------------------------------------------
C     Calculate the  term <mu3|[[X,T2],T2]|HF> 
C     added in W^BD(KW3MAT)
C------------------------------------------------------
C
                   CALL WXBD_T2(AIBJCK_PERM,B,ISYMB,D,ISYMD,WORK(KT2TP),
     *                         ISYM0,
     *                         WORK(KFOCKR1),ISYMR1,
     *                         WORK(KINDSQW),LENSQW,WORK(KW3MAT),ISWMAT,
     *                         WORK(KEND5),LWRK5)
                  END IF
C
C------------------------------------------------------
C     To get the entire T3^X add the two terms
C------------------------------------------------------
C
                  CALL WXBD_GROUND(AIBJCK_PERM,WORK(KT2R1),ISYMR1,
     *                            WORK(KWTEMP),
     *                            WORK(KT3VDG1),WORK(KT3VBG1),DUMMY,
     *                            DUMMY,
     *                            WORK(KT3OG1),ISINT2,
     *                            WORK(KW3MAT),WORK(KEND5),LWRK5,
     *                            WORK(KINDSQW),LENSQW,ISYMB,B,ISYMD,D)
 
                  CALL WXBD_GROUND(AIBJCK_PERM,WORK(KT2TP),ISYM0,
     *                            WORK(KWTEMP),
     *                            WORK(KW3XVDGX1),WORK(KW3XVDGX2),
     *                            DUMMY,DUMMY,
     *                            WORK(KW3XOGX1),ISINT2R1,
     *                            WORK(KW3MAT),WORK(KEND5),LWRK5,
     *                            WORK(KINDSQW),LENSQW,ISYMB,B,ISYMD,D)
C
C------------------------------------------------
C     Divide by the energy difference and
C     remove the forbidden elements
C------------------------------------------------
C
                  CALL WBD_DIA(B,ISYMB,D,ISYMD,FREQR1,ISWMAT,
     *                         WORK(KW3MAT),WORK(KDIAGW),WORK(KFOCKD))
                  CALL T3_FORBIDDEN(WORK(KW3MAT),ISYMR1,ISYMB,B,
     *                              ISYMD,D)
c
c       call sum_pt3(work(KW3MAT),isymb,b,isymd,d,
c    *             ISWMAT,work(kx3am),5)
C
C------------------------------------------------------
C     Calculate the  term <mu3|[X,T3]|HF> occupied contribution 
C     added in W^BD(KW3MATD)
C------------------------------------------------------
C
                  AIBJCK_PERM = 3
C
                  IF (LISTR(1:3).EQ.'R1 ') THEN
                   CALL WXBD_O(AIBJCK_PERM,WORK(KT3MAT),ISCKIJ,
     *                        WORK(KFOCKR1),ISYMR1,
     *                        WORK(KW3MATD),ISWMAT,WORK(KEND5),LWRK5)
C
C------------------------------------------------------
C     Calculate the  term <mu3|[[X,T2],T2]|HF> 
C     added in W^BD(KW3MATD)
C------------------------------------------------------
C
                   CALL WXBD_T2(AIBJCK_PERM,B,ISYMB,D,ISYMD,WORK(KT2TP),
     *                         ISYM0,
     *                         WORK(KFOCKR1),ISYMR1,
     *                        WORK(KINDSQW),LENSQW,WORK(KW3MATD),ISWMAT,
     *                         WORK(KEND5),LWRK5)
                  END IF
C
C------------------------------------------------------
C     To get the entire T3^X add the two terms
C------------------------------------------------------
C
                  CALL WXBD_GROUND(AIBJCK_PERM,WORK(KT2R1),ISYMR1,
     *                            WORK(KWTEMP),
     *                            WORK(KT3VDG1),WORK(KT3VBG1),
     *                            WORK(KT3VBG3),DUMMY,
     *                            WORK(KT3OG1),ISINT2,
     *                            WORK(KW3MATD),WORK(KEND5),LWRK5,
     *                            WORK(KINDSQW),LENSQW,ISYMB,B,ISYMD,D)
 
                  CALL WXBD_GROUND(AIBJCK_PERM,WORK(KT2TP),ISYM0,
     *                            WORK(KWTEMP),
     *                            WORK(KW3XVDGX1),WORK(KW3XVDGX2),
     *                            WORK(KT3VBGX3),DUMMY,
     *                            WORK(KW3XOGX1),ISINT2R1,
     *                            WORK(KW3MATD),WORK(KEND5),LWRK5,
     *                            WORK(KINDSQW),LENSQW,ISYMB,B,ISYMD,D)
C
C------------------------------------------------
C     Divide by the energy difference and
C     remove the forbidden elements
C------------------------------------------------
C
                  CALL WBD_DIA(B,ISYMB,D,ISYMD,FREQR1,ISWMAT,
     *                         WORK(KW3MATD),WORK(KDIAGW),WORK(KFOCKD))
 
                  CALL T3_FORBIDDEN(WORK(KW3MATD),ISYMR1,ISYMB,B,
     *                              ISYMD,D)
c
c       call sum_pt3(work(KW3MATD),isymb,b,isymd,d,
c    *             ISWMAT,work(kx3am),5)
C
C--------------------------------------------
C    Add WORK(KW3MAT) + 0.5*WORK(KW3MATD) ...
C--------------------------------------------
C
                  CALL DAXPY(NCKIJ(ISWMAT),HALF,WORK(KW3MATD),1,
     *                       WORK(KW3MAT),1)
c
c       call sum_pt3(work(KW3MAT),isymb,b,isymd,d,
c    *             ISWMAT,work(kx3am),6)
 
C--------------------------------------------------------
C                 Write WBMAT as WBMAT^D(ai,bj,l) to disc
C--------------------------------------------------------
                  CALL WRITE_T3_DL(LUWBMAT,FNWBMAT,WORK(KW3BMAT),ISYML1,
     *                             ISYMD,ISYMB,B)
C
                  !To conform with real sign of t3 amplitudes
                  CALL DSCAL(NCKIJ(ISCKIJ),-ONE,WORK(KT3MAT),1)
C-------------------------------------------------------------
C                 Write T3 amplitudes as T3^D(ai,bj,l) to disc
C-------------------------------------------------------------
                  CALL WRITE_T3_DL(LUT3,FNT3,WORK(KT3MAT),ISYM0,
     *                             ISYMD,ISYMB,B)
C
C                 
                  !To conform with noddy code
                  CALL DSCAL(NCKIJ(ISWMAT),-ONE,WORK(KW3MAT),1)
C-----------------------------------------------------------------------
C    Write WORK(KW3MAT) + 0.5*WORK(KW3MATD) as KW3MAT^D(ai,bj,l) to disc
C-----------------------------------------------------------------------
                  CALL WRITE_T3_DL(LUTHETA,FNTHETA,WORK(KW3MAT),ISYMR1,
     *                             ISYMD,ISYMB,B)
C
C--------------------------------------------------------
C                 Get T2X T2Y contribution to DIA density 
C--------------------------------------------------------
C
                  IF (DO_DIA) THEN
                   T2XNET2Y = .TRUE.
                   CALL CC_XI_DEN_IA(T2XNET2Y,DIA,WORK(KW3BMAT),ISWBMAT,
     *                               WORK(KT2R1),ISYMR1,
     *                                WORK(KT2TP),ISYM0,WORK(KINDSQWB),
     *                                LENSQWB,
     *                                B,ISYMB,D,ISYMD,WORK(KEND5),LWRK5)
                  END IF
C
               ENDDO   ! B
            ENDDO      ! ISYMB
C
C--------------------------------------------------------------------
C          y^M(fn,im) <- M(fn,im)+ sum_em y^T2(dl,ei) tbar(dl,em,fn)
C--------------------------------------------------------------------
C
           IF (DO_YMMAT) THEN
              
              !if tbar30 (all W intermediates are contained in t3bar)
              !for t3barx it is only a part of W intermediates
              CALL MMATX(YMMAT,LISTL,
     *                   LUWBMAT,FNWBMAT,ISYML1,
     *                   LUT2Y,FNT2Y,ISYMR1,
     *                   ISYMD,D,
     *                   WORK(KEND5),LWRK5)
C
              IF ( (LISTL(1:3) .EQ. 'L1 ') .OR. (LISTL(1:3) .EQ. 'M1 ')
     *            .OR. (LISTL(1:3) .EQ. 'N2 ')
     *            .OR. (LISTL(1:3) .EQ. 'LE ') ) THEN
                 ! for t3barx calculate the remaining one W intermediate
                 CALL MMATXL1(YMMAT,LISTL,LISTR,IDLSTR,
     *                        LUWBMAT,FNWBMAT,ISYML1,
     *                        ISYMD,D,
     *                        WORK(KEND5),LWRK5)
              END IF
           END IF

C
C-------------------------------------------------------
C          Get DAB0 and DIJ0 densities
C-------------------------------------------------------
C
            IF (DO_DIA) THEN
               QUADR = .FALSE.
               CALL CC_XI_DEN_ABIJ(QUADR,LISTR,WORK(KDAB0),WORK(KDIJ0),
     *                             .FALSE.,DUMMY,
     *                             DUMMY,IDUMMY,
     *                             IDUMMY,DUMMY,
     *                             ISYM0,ISYML1,IDUMMY,
     *                             LUT3,FNT3,LUWBMAT,FNWBMAT,
     *                             IDUMMY,CDUMMY,
     *                             DUMMY,DUMMY,
     *                             WORK(KEND5),LWRK5,ISYMD,D)
            END IF
C
C--------------------------------------------------------------------------
C           Get DAB,DIJ,DIA densities contributions for quadratice response
C--------------------------------------------------------------------------
C
            QUADR = .TRUE.
            CALL CC_XI_DEN_ABIJ(QUADR,LISTR,DAB,DIJ,DO_DIA,DIA,
     *                          WORK(KL2L1),ISYML1,
     *                          ISYMR1,WORK(KFOCKR1),
     *                          ISYM0,ISYML1,ISYMR1,
     *                          LUT3,FNT3,LUWBMAT,FNWBMAT,
     *                          LUTHETA,FNTHETA,
     *                          WORK(KFOCKD),FREQR1,
     *                          WORK(KEND5),LWRK5,ISYMD,D)
C
         ENDDO       ! D
      ENDDO          ! ISYMD 
C
      IF (DO_DIA) THEN
        CALL CC3_XI_DEN_AI_T1(DIA,ISYDEN,WORK(KDAB0),WORK(KDIJ0),ISYML1,
     *                              WORK(KT1R1),ISYMR1)
      END IF
C
COMMENT COMMENT
c      write(lupri,*) 'w3x (usual) in CC3_ADENVIR'
c      write(lupri,*) 'w3xD  in CC3_ADENVIR'
c      write(lupri,*) 'w3bx  in CC3_ADENVIR'
c      write(lupri,*) 'DO_YMMAT = ', DO_YMMAT
c      write(lupri,*) 'w3x + 0.5w3xD in CC3_ADENVIR'
c      call print_pt3(work(kx3am),ISYM0,4)
cOMMENT COMMENT
C
C---------------------------------
C     Close the files
C---------------------------------
C
      CALL WCLOSE2(LUT3,FNT3,'DELETE')
      CALL WCLOSE2(LUWBMAT,FNWBMAT,'DELETE')
      CALL WCLOSE2(LUTHETA,FNTHETA,'DELETE')
C
      IF (DO_YMMAT) THEN 
         CALL WCLOSE2(LUT2Y,FNT2Y,'DELETE')
      END IF
C
C--------------------------------
C     Close files for "response"
C--------------------------------
C
      CALL WCLOSE2(LU3SRTR,FN3SRTR,'DELETE')
      CALL WCLOSE2(LUCKJDR,FNCKJDR,'DELETE')
      CALL WCLOSE2(LUDELDR,FNDELDR,'DELETE')
      CALL WCLOSE2(LUDKBCR,FNDKBCR,'DELETE')
C
C-------------
C     End
C-------------
C
C
      CALL QEXIT('CC3_ADENVIR')
C
      RETURN
      END
C  /* Deck wxbd_o */
      SUBROUTINE WXBD_O(AIBJCK_PERM,TMAT,ISTMAT,FOCKY,ISYFKY,
     *                 WMAT,ISWMAT,WRK,LWRK)
C 
C If (AIBJCK_PERM.eq.1) then (bjdk) permutation symmetry 
C
C     WBD(aikj) = WBD(aikj) - t(aBD,ljk) * fock(li)
C
C                            tmatBD(alkj)
C
C else if (AIBJCK_PERM.eq.2) then (aidk) permutation symmetry
C
C     WBD(aikj) = WBD(aikj) - t(aBD,ilk) * fock(lj)
C
C                            tmatBD(aikl)
C
C else if (AIBJCK_PERM.eq.3) then (aibj) permutation symmetry
C
C     WBD(aikj) = WBD(aikj) - t(aBD,ijl) * fock(lk)
C
C                            tmatBD(ailj)
C
C else if (AIBJCK_PERM.eq.4) then calculate all terms
C 
C
C Written by P. Jorgensen and F. Pawlowski, Spring 2002.
C (modyfied for AIBJCK_PERM flag - spring 2003.)
C

      IMPLICIT NONE
C
      INTEGER AIBJCK_PERM
C
      INTEGER LWRK, KFCLI, KEND0, LWRK0, KOFF1, KOFF2 
      INTEGER NL, KOFFY, KOFFT, KOFFW 
      INTEGER ISTMAT, ISYFKY, ISWMAT, ISALKJ
      INTEGER ISYMA, ISYAI, ISYAIK, ISYALK, ISYAL, NA
      INTEGER ISYMJ, ISYMK, ISYMI, ISYML, ISYFI
      INTEGER ISYAIL,NAI,NAIK
C
#if defined (SYS_CRAY)
      REAL TMAT(*), FOCKY(*), WMAT(*), WRK(*)
      REAL HALF, ONE 
      REAL XNORMVAL,DDOT
#else
      DOUBLE PRECISION TMAT(*), FOCKY(*), WMAT(*), WRK(*)
      DOUBLE PRECISION HALF, ONE
      DOUBLE PRECISION XNORMVAL,DDOT
#endif
C
#include <priunit.h>
#include <dummy.h>
#include <iratdef.h>
#include <ccsdsym.h>
#include <inftap.h>
#include <ccinftap.h>
#include <ccorb.h>
#include <ccsdinp.h>
C
      PARAMETER (HALF = 0.5D0, ONE = 1.0D0)
C
      CALL QENTER('WXBDO')
C
C---------------------------------------
C     Initial test of AIBJCK_PERM option
C---------------------------------------
C
      IF ( (AIBJCK_PERM .LT. 1) .OR. (AIBJCK_PERM .GT. 4) ) THEN
         WRITE(LUPRI,*)'AIBJCK_PERM = ',AIBJCK_PERM
         WRITE(LUPRI,*)'should be between 1 and 4 '
         CALL QUIT('Illegal value of AIBJCK_PERM option in WXBD_O')
      END IF
C
C RESORT OCC-OCC  FOCKY ELEMENTS (L,I)
C
C
      KFCLI  = 1
      KEND0  = KFCLI + NMATIJ(ISYFKY)
      LWRK0  = LWRK  - KEND0
C
      IF (LWRK0 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWRK0
         WRITE(LUPRI,*) 'Memory needed    : ',KEND0
         CALL QUIT('Insufficient space in WXBD_O')
      END IF
C   
      DO ISYMI = 1,NSYM
         ISYML = MULD2H(ISYMI,ISYFKY)
         DO I = 1,NRHF(ISYMI)
             KOFF1 = IFCRHF(ISYML,ISYMI) + NORB(ISYML)*(I - 1) + 1
             KOFF2 = KFCLI + IMATIJ(ISYML,ISYMI) + NRHF(ISYML)*(I - 1)
             CALL DCOPY(NRHF(ISYML),FOCKY(KOFF1),1,WRK(KOFF2),1) 
         END DO
      END DO
C
      IF ((AIBJCK_PERM.EQ.1) .OR. (AIBJCK_PERM.EQ.4)) THEN
C
C CARRY OUT MATRIX MULTIPLICATION
C WBD(a,i,k,j) = WBD(a,i-,k,j) - sum (l) tmatBD(a,l,k,j)*focky(l,i)
C 
         ISALKJ = ISTMAT
         DO ISYMJ = 1,NSYM
            ISYALK =MULD2H(ISYMJ,ISALKJ)
            DO J = 1,NRHF(ISYMJ)
               DO ISYMK = 1,NSYM
                  ISYAL = MULD2H(ISYMK,ISYALK)
                  DO K  = 1,NRHF(ISYMK)
                     DO ISYML = 1,NSYM
                        ISYMA = MULD2H(ISYAL,ISYML)   
                        ISYMI = MULD2H(ISYFKY,ISYML)
                        ISYAIK = MULD2H(ISWMAT,ISYMJ)
                        ISYAI = MULD2H(ISYAIK,ISYMK)
                        NA    = MAX(1,NVIR(ISYMA))
                        NL    = MAX(1,NRHF(ISYML))
                        KOFFY = KFCLI + IMATIJ(ISYML,ISYMI) 
                        KOFFT = ISAIKJ(ISYALK,ISYMJ)
     *                        + NCKI(ISYALK)*(J-1)
     *                        + ISAIK(ISYAL,ISYMK) 
     *                        + NT1AM(ISYAL)*(K-1)
     *                        + IT1AM(ISYMA,ISYML) + 1
                        KOFFW = ISAIKJ(ISYAIK,ISYMJ)
     *                        + NCKI(ISYAIK)*(J-1)
     *                        + ISAIK(ISYAI,ISYMK) 
     *                        + NT1AM(ISYAI)*(K-1)
     *                        + IT1AM(ISYMA,ISYMI) + 1
C
C SYMMETRY BETWEEN BJ AND CK INTRODUCE A FACTOR TWO
C DENOTE t3 IS CALCULATED WITH NEGATIVE SIGN
C
                        CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),
     *                             NRHF(ISYML),ONE,TMAT(KOFFT),NA,
     *                             WRK(KOFFY),NL,ONE,WMAT(KOFFW),NA)
                     END DO
                  END DO
               END DO
            END DO
         END DO
C
      END IF
      IF ((AIBJCK_PERM.EQ.2) .OR. (AIBJCK_PERM.EQ.4)) THEN
C
C     WBD(aikj) = WBD(aikj) - t(aBD,ilk) * fock(lj)
C
C                            tmatBD(aikl)
C
C WBD(a,i,k,j) = WBD(a,i,k-,j) - sum (l) tmatBD(a,i,k,l)*focky(l,j)
C
         DO ISYMJ = 1,NSYM
            ISYML = MULD2H(ISYFKY,ISYMJ)
            ISYAIK =MULD2H(ISTMAT,ISYML)
            NAIK    = MAX(1,NCKI(ISYAIK))
            NL    = MAX(1,NRHF(ISYML))
            KOFFY = KFCLI + IMATIJ(ISYML,ISYMJ)
            KOFFT = ISAIKJ(ISYAIK,ISYML)
     *            + 1
            KOFFW = ISAIKJ(ISYAIK,ISYMJ)
     *            + 1
C
C DENOTE t3 IS CALCULATED WITH NEGATIVE SIGN
C
            CALL DGEMM('N','N',NCKI(ISYAIK),NRHF(ISYMJ),
     *                 NRHF(ISYML),ONE,TMAT(KOFFT),NAIK,
     *                 WRK(KOFFY),NL,ONE,WMAT(KOFFW),NAIK)
         END DO

      END IF
      IF ((AIBJCK_PERM.EQ.3) .OR. (AIBJCK_PERM.EQ.4)) THEN
C
C CARRY OUT MATRIX MULTIPLICATION
C WBD(a,i,k,j) = WBD(a,i,k-,j) - sum (l) tmatBD(a,i,l,j)*focky(l,k)
C 
         DO ISYMJ = 1,NSYM
            ISYAIL =MULD2H(ISTMAT,ISYMJ)
            ISYAIK = MULD2H(ISWMAT,ISYMJ)
            DO J = 1,NRHF(ISYMJ)
               DO ISYMK = 1,NSYM
                  ISYAI = MULD2H(ISYAIK,ISYMK)
                  ISYML = MULD2H(ISYFKY,ISYMK)
                     NAI    = MAX(1,NT1AM(ISYAI))
                     NL    = MAX(1,NRHF(ISYML))
                     KOFFY = KFCLI + IMATIJ(ISYML,ISYMK) 
                     KOFFT = ISAIKJ(ISYAIL,ISYMJ)
     *                     + NCKI(ISYAIL)*(J-1)
     *                     + ISAIK(ISYAI,ISYML) 
     *                     + 1
                     KOFFW = ISAIKJ(ISYAIK,ISYMJ)
     *                     + NCKI(ISYAIK)*(J-1)
     *                     + ISAIK(ISYAI,ISYMK) 
     *                     + 1
C
C DENOTE t3 IS CALCULATED WITH NEGATIVE SIGN
C
                     CALL DGEMM('N','N',NT1AM(ISYAI),NRHF(ISYMK),
     *                          NRHF(ISYML),ONE,TMAT(KOFFT),NAI,
     *                          WRK(KOFFY),NL,ONE,WMAT(KOFFW),NAI)
               END DO   
            END DO   
         END DO   
C
      END IF
C
      CALL QEXIT('WXBDO')
C
      RETURN
      END
C
C  /* Deck wxbd_t2 */
      SUBROUTINE WXBD_T2(AIBJCK_PERM,B,ISYMB,D,ISYMD,T2TP,ISYMT2,FOCKY,
     *                   ISYFKY,INDSQ,LENSQ,WMAT,ISWMAT,WRK,LWRK)
C 
C 
C If (AIBJCK_PERM.eq.1) then (aibjdk) + (aidkbj)  permutation 
C
C     WBD(aikj) = WBD(aikj)  -  focky(l,f)* t2(ai,dl)*t2(fk,bj) 
C
C                            - focky(l,f)*t2(ai,bl)*t2(fj,dk)
C
C else (AIBJCK_PERM.eq.2)  then  (bjdkai) + (bjaidk)  permutation
C
C     WBD(aikj) = WBD(aikj)  -  focky(l,f)* t2(bj,al)*t2(fi,dk)
C 
C                            -  focky(l,f)* t2(bj,dl)*t2(fk,ai)
C 

C else (AIBJCK_PERM.eq.3)  then  (dkbjai) + (dkaibj)  permutation
C
C     WBD(aikj) = WBD(aikj)  -  focky(l,f)* t2(dk,al)*t2(fi,bj)
C 
C                            -  focky(l,f)* t2(dk,bl)*t2(fj,ai)
C
C else if (AIBJCK_PERM.eq.4) then calculate all terms
C 
C     Written by P. Jorgensen and F. Pawlowski, Spring 2002.
C

      IMPLICIT NONE
C
      INTEGER AIBJCK_PERM
C
      INTEGER LENSQ
      INTEGER INDSQ(LENSQ,6)
      INTEGER LWRK,KFCLF, KEND0, LWRK0, KTB, KEND1, LWRK1 
      INTEGER NL, NF, KOFFY, KOFFT2, KOFFT, KOFFW, KTD, KW
      INTEGER ISYMB, ISYMD, ISYMT2, ISYFKY, ISWMAT 
      INTEGER ISYAIL, ISYAI, ISYAIK, NA, NAI, LENGTH
      INTEGER ISYFIJ,ISYLIJ,ISYAKL,ISYMJ,ISYFI,ISYMI,ISYMF,ISYML
      INTEGER ISYFKJ,ISYTB,ISYMK,ISYFJK,ISYTD,ISYFJ,ISYLJ
      INTEGER ISYAK,ISYAKI,NAK
      INTEGER ISYBD,ISYLK,ISYFK,ISYAIJ,NAIJ,ISYLI
      INTEGER KLIJ,KAKL,KLK,KFK
      INTEGER ISYFIK,ISYLIK,ISYAJL,KLIK,KAJL
      INTEGER ISYAJ,ISYAJI,NAJ
      INTEGER KLJ,KFJ,NAIK
C
      INTEGER KOFF1,KOFF2,KOFF3
C
#if defined (SYS_CRAY)
      REAL T2TP(*), FOCKY(*), WMAT(*), WRK(*)
      REAL HALF, ONE, ZERO
#else
      DOUBLE PRECISION T2TP(*), FOCKY(*), WMAT(*), WRK(*)
      DOUBLE PRECISION HALF, ONE, ZERO
#endif
C
#include <priunit.h>
#include <dummy.h>
#include <iratdef.h>
#include <ccsdsym.h>
#include <inftap.h>
#include <ccinftap.h>
#include <ccorb.h>
#include <ccsdinp.h>
C
      PARAMETER (HALF = 0.5D0, ONE = 1.0D0, ZERO = 0.0D0)
C
      CALL QENTER('WXBDT2')
C
C---------------------------------------
C     Initial test of AIBJCK_PERM option
C---------------------------------------
C
      IF ( (AIBJCK_PERM .LT. 1) .OR. (AIBJCK_PERM .GT. 4) ) THEN
         WRITE(LUPRI,*)'AIBJCK_PERM = ',AIBJCK_PERM
         WRITE(LUPRI,*)'should be between 1 and 4 '
         CALL QUIT('Illegal value of AIBJCK_PERM option in WXBD_T2')
      END IF
C
C
C RESORT VIR-OCC  FOCKY ELEMENTS (l,f)
C
C
      KW = 1
      KFCLF = KW + NCKIJ(ISWMAT)
      KEND0  = KFCLF + NT1AM(ISYFKY)
      LWRK0  = LWRK  - KEND0
      CALL DZERO(WRK(KW),NCKIJ(ISWMAT))
C
      IF (LWRK0 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWRK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND0
         CALL QUIT('Insufficient space in WXBD_T2 (1)')
      END IF
C
      DO ISYMF = 1,NSYM
         ISYML = MULD2H(ISYMF,ISYFKY)
         DO L = 1,NRHF(ISYML)
            DO F = 1,NVIR(ISYMF)
               KOFF1 = IFCVIR(ISYML,ISYMF) + NORB(ISYML)*(F - 1) + L
               KOFF2 = KFCLF +  IT1AMT(ISYML,ISYMF) 
     *               + NRHF(ISYML)*(F - 1) + L -1 
C
                  WRK(KOFF2) = FOCKY(KOFF1)
C
            END DO
         END DO
      END DO
C
      IF ((AIBJCK_PERM.EQ.1) .OR. (AIBJCK_PERM.EQ.4)) THEN
C
C    calculate first t2 contribution to W matrix
C
C construct tB(l,k,j) = sum (f) focky(l,f)*t2tp(f,k,j,B)
C
         ISYFKJ   = MULD2H(ISYMT2,ISYMB) 
         ISYTB    = MULD2H(ISYFKY,ISYFKJ) 
         KTB      = KEND0
         KEND1    = KTB  + NMAIJK(ISYTB)
         LWRK1    = LWRK  - KEND1
C
         CALL DZERO(WRK(KTB),NMAIJK(ISYTB))
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWRK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND1
            CALL QUIT('Insufficient space in WXBD_T2 (2)')
         END IF
C
         DO ISYMJ = 1,NSYM
            ISYFK  = MULD2H(ISYFKJ,ISYMJ) 
            DO J  = 1,NRHF(ISYMJ)
               DO ISYMK = 1,NSYM
                  ISYMF = MULD2H(ISYFK,ISYMK)
                  ISYML = MULD2H(ISYFKY,ISYMF)
                  ISYLK  = MULD2H(ISYML,ISYMK)
                  NL = MAX(1,NRHF(ISYML))
                  NF = MAX(1,NVIR(ISYMF))
                  KOFFY  = KFCLF + IT1AMT(ISYML,ISYMF)  
                  KOFFT2 = IT2SP(ISYFKJ,ISYMB) + NCKI(ISYFKJ)*(B-1)
     *                    + ISAIK(ISYFK,ISYMJ) + NT1AM(ISYFK)*(J-1)
     *                    + IT1AM(ISYMF,ISYMK) + 1
                  KOFFT =  KTB + IMAIJK(ISYLK,ISYMJ)  
     *                         + NMATIJ(ISYLK)*(J-1) 
     *                         + IMATIJ(ISYML,ISYMK)
C
                  CALL DGEMM('N','N',NRHF(ISYML),NRHF(ISYMK),
     *                    NVIR(ISYMF),ONE,WRK(KOFFY),NL,
     *                    T2TP(KOFFT2),NF,ONE,WRK(KOFFT),NL)
C
               END DO
            END DO
         END DO
C
C         WBD(a,i,k,j) = WBD(a,i,k,j) -
C                        sum (f,l) focky(l,f)* t2(ai,Dl)*t2(fk,Bj) 
C                      = WBD(a,i,k,j) - 
C                        sum(l) t2tp(a,i,l,D) * tB(l,k,j)
C
         ISYAIL = MULD2H(ISYMT2,ISYMD)
         DO ISYMJ = 1,NSYM
            ISYLK  = MULD2H(ISYTB,ISYMJ)
            DO J  = 1,NRHF(ISYMJ)
               DO ISYMK = 1,NSYM
                  ISYML = MULD2H(ISYLK,ISYMK)
                  ISYAI = MULD2H(ISYAIL,ISYML) 
                  ISYAIK = MULD2H(ISYAI,ISYMK)
                  NAI = MAX(1,NT1AM(ISYAI))
                  NL = MAX(1,NRHF(ISYML)) 
                  KOFFT2 = IT2SP(ISYAIL,ISYMD) + NCKI(ISYAIL)*(D-1)
     *                    + ISAIK(ISYAI,ISYML) + 1
                  KOFFT =  KTB + IMAIJK(ISYLK,ISYMJ) 
     *                         + NMATIJ(ISYLK)*(J-1) 
     *                         + IMATIJ(ISYML,ISYMK)
                  KOFFW  = ISAIKJ(ISYAIK,ISYMJ) + NCKI(ISYAIK)*(J-1)
     *                    + ISAIK(ISYAI,ISYMK) + 1
                  CALL DGEMM('N','N',NT1AM(ISYAI),NRHF(ISYMK),
     *                       NRHF(ISYML),-ONE,T2TP(KOFFT2),NAI,
     *                       WRK(KOFFT),NL,ONE,WMAT(KOFFW),NAI)

C
               END DO
            END DO
         END DO
C
C    calculate second t2 contribution to W matrix
C
C
C construct tD(l,j,k) = sum (f) focky(l,f)*t2tp(f,j,k,D)
C
         ISYFJK   = MULD2H(ISYMT2,ISYMD) 
         ISYTD    = MULD2H(ISYFKY,ISYFJK) 
         KTD      = KEND0
         KEND1    = KTD  + NMAIJK(ISYTD)
         LWRK1    = LWRK  - KEND1
C
         CALL DZERO(WRK(KTD),NMAIJK(ISYTD))
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWRK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND1
            CALL QUIT('Insufficient space in WXBD_T2 (3)')
         END IF
C

         DO ISYMK = 1,NSYM
            ISYFJ  = MULD2H(ISYFJK,ISYMK) 
            DO K  = 1,NRHF(ISYMK)
               DO ISYMJ = 1,NSYM
                  ISYMF = MULD2H(ISYFJ,ISYMJ)
                  ISYML = MULD2H(ISYFKY,ISYMF)
                  ISYLJ  = MULD2H(ISYML,ISYMJ)
                  NL = MAX(1,NRHF(ISYML))
                  NF = MAX(1,NVIR(ISYMF))
                  KOFFY  = KFCLF + IT1AMT(ISYML,ISYMF)  
                  KOFFT2 = IT2SP(ISYFJK,ISYMD) + NCKI(ISYFJK)*(D-1)
     *                    + ISAIK(ISYFJ,ISYMK) + NT1AM(ISYFJ)*(K-1)
     *                    + IT1AM(ISYMF,ISYMJ) + 1
                  KOFFT =  KTD + IMAIJK(ISYLJ,ISYMK) 
     *                    + NMATIJ(ISYLJ)*(K-1) 
     *                    + IMATIJ(ISYML,ISYMJ)
                  CALL DGEMM('N','N',NRHF(ISYML),NRHF(ISYMJ),
     *                       NVIR(ISYMF),ONE,WRK(KOFFY),NL,
     *                       T2TP(KOFFT2),NF,ONE,WRK(KOFFT),NL)
C
               END DO
            END DO
         END DO
C
C      WBD(a,i,k,j) = WBD(a,i,k,j) -
C                        sum (f,l) focky(l,f)*t2(ai,Bl)*t2(fj,Dk) )
C                   = WBD(a,i,k,j) -
C                        sum(l) t2tp(a,i,l,B) * tD(l,j,k)
C
         ISYAIL = MULD2H(ISYMT2,ISYMB)
         DO ISYMK = 1,NSYM
            ISYLJ  = MULD2H(ISYTD,ISYMK)
            DO K  = 1,NRHF(ISYMK)
               DO ISYMJ = 1,NSYM
                  ISYML = MULD2H(ISYLJ,ISYMJ)
                  ISYAI = MULD2H(ISYAIL,ISYML) 
                  ISYAIJ = MULD2H(ISYAI,ISYMJ) 
                  NAI = MAX(1,NT1AM(ISYAI))
                  NL = MAX(1,NRHF(ISYML)) 
                  KOFFT2 = IT2SP(ISYAIL,ISYMB) + NCKI(ISYAIL)*(B-1)
     *                    + ISAIK(ISYAI,ISYML) + 1
                  KOFFT =  KTD + IMAIJK(ISYLJ,ISYMK) 
     *                         + NMATIJ(ISYLJ)*(K-1) 
     *                         + IMATIJ(ISYML,ISYMJ)
                  KOFFW  = KW  + ISAIKJ(ISYAIJ,ISYMK) 
     *                         + NCKI(ISYAIJ)*(K-1)
     *                         + ISAIK(ISYAI,ISYMJ) 
                  CALL DGEMM('N','N',NT1AM(ISYAI),NRHF(ISYMJ),
     *                       NRHF(ISYML),-ONE,T2TP(KOFFT2),NAI,
     *                       WRK(KOFFT),NL,ONE,WRK(KOFFW),NAI)

C
               END DO
            END DO
         END DO
C
C     change order aijk to aikj
C
         DO I = 1,NCKIJ(ISWMAT)
            WMAT(I) = WMAT(I) + WRK(INDSQ(I,3))
         END DO
C
C
      END IF
      IF ((AIBJCK_PERM.EQ.2) .OR. (AIBJCK_PERM.EQ.4)) THEN
C
C
C     WBD(aikj) = WBD(aikj)  -  focky(l,f)* t2(bj,al)*t2(fi,dk)
C
C                                          T^B(ajl)    t(fikD)
C                  
C                work(lik) = focky(lf) * t(fikD) 
C
C                work(ajik) = T^B(ajl) * work(lik)

C
C                work(lij) = focky(lf) * T^B(fij)
C
         ISYFIK     = MULD2H(ISYMT2,ISYMD)
         ISYLIK     = MULD2H(ISYFKY,ISYFIK)
         ISYAJL     = MULD2H(ISYMT2,ISYMB)
         KLIK       = KEND0
         KAJL       = KLIK  + NMAIJK(ISYLIK)
         KEND1      = KAJL  + NCKI(ISYAJL)
         LWRK1      = LWRK   - KEND1
C
         CALL DZERO(WRK(KLIK),NMAIJK(ISYLIK))
C
         CALL DZERO(WRK(KW),NCKIJ(ISWMAT))
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWRK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND1
            CALL QUIT('Insufficient space in WXBD_T2 (4)')
         END IF
C
C

         DO ISYMK = 1,NSYM
            ISYFI  = MULD2H(ISYFIK,ISYMK)
            DO K  = 1,NRHF(ISYMK)
               DO ISYMI = 1,NSYM
                  ISYMF = MULD2H(ISYFI,ISYMI)
                  ISYML = MULD2H(ISYFKY,ISYMF)
                  ISYLI  = MULD2H(ISYMI,ISYML)
                  NL = MAX(1,NRHF(ISYML))
                  NF = MAX(1,NVIR(ISYMF))
                  KOFF1  = KFCLF + IT1AMT(ISYML,ISYMF)
                  KOFF2 = IT2SP(ISYFIK,ISYMD) + NCKI(ISYFIK)*(D-1)
     *                    + ISAIK(ISYFI,ISYMK) + NT1AM(ISYFI)*(K-1)
     *                    + IT1AM(ISYMF,ISYMI) + 1
                  KOFF3 =  KLIK + IMAIJK(ISYLI,ISYMK)
     *                    + NMATIJ(ISYLI)*(K-1)
     *                    + IMATIJ(ISYML,ISYMI)
C
C                work(lik) = focky(lf) * t(fikD) 
C
                  CALL DGEMM('N','N',NRHF(ISYML),NRHF(ISYMI),
     *                       NVIR(ISYMF),ONE,WRK(KOFF1),NL,
     *                       T2TP(KOFF2),NF,ONE,WRK(KOFF3),NL)
C
               END DO
            END DO
         END DO
C
C        T^B(ajl) =     t2(bj,al)
C
         CALL SORT_T2_AJI(WRK(KAJL),ISYMB,B,T2TP,ISYMT2)
C
C
         DO ISYMK = 1,NSYM
            ISYLI  = MULD2H(ISYLIK,ISYMK)
            DO K  = 1,NRHF(ISYMK)
               DO ISYMI = 1,NSYM
                  ISYML = MULD2H(ISYLI,ISYMI)
                  ISYAJ = MULD2H(ISYAJL,ISYML)
                  ISYAJI = MULD2H(ISYAJ,ISYMI)
                  NAJ = MAX(1,NT1AM(ISYAJ))
                  NL = MAX(1,NRHF(ISYML))
                  KOFF1 = KAJL + ISAIK(ISYAJ,ISYML)
                  KOFF2 = KLIK + IMAIJK(ISYLI,ISYMK)
     *                         + NMATIJ(ISYLI)*(K-1)
     *                         + IMATIJ(ISYML,ISYMI)
                  KOFF3  = KW  + ISAIKJ(ISYAJI,ISYMK)
     *                         + NCKI(ISYAJI)*(K-1)
     *                         + ISAIK(ISYAJ,ISYMI)
C
C                work(ajik) = T^B(ajl) * work(lik)
C
                  CALL DGEMM('N','N',NT1AM(ISYAJ),NRHF(ISYMI),
     *                       NRHF(ISYML),-ONE,WRK(KOFF1),NAJ,
     *                       WRK(KOFF2),NL,ONE,WRK(KOFF3),NAJ)

C
               END DO
            END DO
         END DO
C
C     change order ajik to aikj
C
         DO I = 1,NCKIJ(ISWMAT)
            WMAT(I) = WMAT(I) + WRK(INDSQ(I,4))
         END DO
C
C     WBD(aikj) = WBD(aikj) - focky(l,f)* t2(bj,dl)*t2(fk,ai)
C                                                     
C                                           T^DB(lj)  t(aikf)
C             
C                work(fj) = focky(l,f) * T^DB(lj)
C                
C              WMAT(aikj) = WMAT(aikj) - t(aikf) * work(fj) 
C
         ISYBD    = MULD2H(ISYMB,ISYMD)
         ISYLJ    = MULD2H(ISYBD,ISYMT2)
         ISYFJ    = MULD2H(ISYFKY,ISYLJ)
C
         KLJ      = KEND0
         KFJ      = KLJ  + NMATIJ(ISYLJ)
         KEND1    = KFJ  + NT1AM(ISYFJ)
         LWRK1    = LWRK  - KEND1
C
         CALL DZERO(WRK(KFJ),NT1AM(ISYFJ))
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWRK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND1
            CALL QUIT('Insufficient space in WXBD_T2 (5)')
         END IF
C
C
         CALL  SORT_T2_IJ(WRK(KLJ),ISYMD,D,ISYMB,B,T2TP,ISYMT2)
C
C                work(fj) = focky(l,f) * T^DB(lj)
C
         DO ISYMJ = 1,NSYM
            ISYML = MULD2H(ISYMJ,ISYLJ)
            ISYMF = MULD2H(ISYFKY,ISYML)
C
            KOFF1  = KFCLF + IT1AMT(ISYML,ISYMF)
            KOFF2  = KLJ   + IMATIJ(ISYML,ISYMJ)
            KOFF3 =  KFJ   + IT1AM(ISYMF,ISYMJ)
C
            NF = MAX(1,NVIR(ISYMF))
            NL = MAX(1,NRHF(ISYML))
C
            CALL DGEMM('T','N',NVIR(ISYMF),NRHF(ISYMJ),
     *                 NRHF(ISYML),ONE,WRK(KOFF1),NL,
     *                 WRK(KOFF2),NL,ONE,WRK(KOFF3),NF)
C
         END DO
C                
C              WMAT(aikj) = WMAT(aikj) - t(aikf) * work(fj) 
C
         DO ISYMJ = 1,NSYM
            ISYMF = MULD2H(ISYMJ,ISYFJ)
            ISYAIK = MULD2H(ISYMT2,ISYMF)
C
            KOFF1 =  IT2SP(ISYAIK,ISYMF) + 1
            KOFF2 =  KFJ +  IT1AM(ISYMF,ISYMJ)
            KOFF3 =  ISAIKJ(ISYAIK,ISYMJ) + 1
C
            NAIK = MAX(1,NCKI(ISYAIK))
            NF = MAX(1,NVIR(ISYMF))

            CALL DGEMM('N','N',NCKI(ISYAIK),NRHF(ISYMJ),
     *                 NVIR(ISYMF),-ONE,T2TP(KOFF1),NAIK,
     *                 WRK(KOFF2),NF,ONE,WMAT(KOFF3),NAIK)
C
         END DO
C
      END IF
      IF ((AIBJCK_PERM.EQ.3) .OR. (AIBJCK_PERM.EQ.4)) THEN
C
C     WBD(aikj) = WBD(aikj)  -  focky(l,f)* t2(dk,al)*t2(fi,bj)
C 
C                                            I^D(alk)  T^B(fij)
C
C                work(lij) = focky(lf) * T^B(fij)
C                 
C                work(akij) = T^D(akl) * work(lij)
C
C     WBD(aikj) = WBD(aikj)  -  work(akij)
C
C
C                work(lij) = focky(lf) * T^B(fij)
C
         ISYFIJ     = MULD2H(ISYMT2,ISYMB)
         ISYLIJ     = MULD2H(ISYFKY,ISYFIJ)
         ISYAKL     = MULD2H(ISYMT2,ISYMD)
         KLIJ       = KEND0
         KAKL       = KLIJ  + NMAIJK(ISYLIJ)
         KEND1      = KAKL  + NCKI(ISYAKL)
         LWRK1      = LWRK   - KEND1
C
         CALL DZERO(WRK(KLIJ),NMAIJK(ISYLIJ))
C
         CALL DZERO(WRK(KW),NCKIJ(ISWMAT))
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWRK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND1
            CALL QUIT('Insufficient space in WXBD_T2 (6)')
         END IF
C

         DO ISYMJ = 1,NSYM
            ISYFI  = MULD2H(ISYFIJ,ISYMJ)
            DO J  = 1,NRHF(ISYMJ) 
               DO ISYMI = 1,NSYM
                  ISYMF = MULD2H(ISYFI,ISYMI)
                  ISYML = MULD2H(ISYFKY,ISYMF)
                  ISYLI  = MULD2H(ISYMI,ISYML)
                  NL = MAX(1,NRHF(ISYML))
                  NF = MAX(1,NVIR(ISYMF))
                  KOFF1  = KFCLF + IT1AMT(ISYML,ISYMF)
                  KOFF2 = IT2SP(ISYFIJ,ISYMB) + NCKI(ISYFIJ)*(B-1)
     *                    + ISAIK(ISYFI,ISYMJ) + NT1AM(ISYFI)*(J-1)
     *                    + IT1AM(ISYMF,ISYMI) + 1
                  KOFF3 =  KLIJ + IMAIJK(ISYLI,ISYMJ)
     *                    + NMATIJ(ISYLI)*(J-1)
     *                    + IMATIJ(ISYML,ISYMI)
                  CALL DGEMM('N','N',NRHF(ISYML),NRHF(ISYMI),
     *                       NVIR(ISYMF),ONE,WRK(KOFF1),NL,
     *                       T2TP(KOFF2),NF,ONE,WRK(KOFF3),NL)
C
               END DO
            END DO
         END DO
C
C                work(akij) = T^D(akl) * work(lij)
C
         CALL SORT_T2_AJI(WRK(KAKL),ISYMD,D,T2TP,ISYMT2)
C
C
         DO ISYMJ = 1,NSYM
            ISYLI  = MULD2H(ISYLIJ,ISYMJ)
            DO J  = 1,NRHF(ISYMJ)
               DO ISYMI = 1,NSYM
                  ISYML = MULD2H(ISYLI,ISYMI)
                  ISYAK = MULD2H(ISYAKL,ISYML)
                  ISYAKI = MULD2H(ISYAK,ISYMI)
                  NAK = MAX(1,NT1AM(ISYAK))
                  NL = MAX(1,NRHF(ISYML))
                  KOFF1 = KAKL + ISAIK(ISYAK,ISYML) 
                  KOFF2 = KLIJ + IMAIJK(ISYLI,ISYMJ)
     *                         + NMATIJ(ISYLI)*(J-1)
     *                         + IMATIJ(ISYML,ISYMI)
                  KOFF3  = KW  + ISAIKJ(ISYAKI,ISYMJ) 
     *                         + NCKI(ISYAKI)*(J-1)
     *                         + ISAIK(ISYAK,ISYMI) 
                  CALL DGEMM('N','N',NT1AM(ISYAK),NRHF(ISYMI),
     *                       NRHF(ISYML),-ONE,WRK(KOFF1),NAK,
     *                       WRK(KOFF2),NL,ONE,WRK(KOFF3),NAK)

C
               END DO
            END DO
         END DO
C
C     change order akij to aikj
C
         DO I = 1,NCKIJ(ISWMAT)
            WMAT(I) = WMAT(I) + WRK(INDSQ(I,1))
         END DO
C
C     WBD(aikj) = WBD(aikj)  -  focky(l,f)* t2(dk,bl)*t2(fj,ai)
C                                            
C                                           I^BD(lk)  I(aijf)
C
C                  work(fk) = focky(lf) * I^BD(lk)
C                  
C                  work(aijk) = I(aijf) * work(fk)
C
C     WBD(aikj) = WBD(aikj)  - work(aijk)
C
         ISYBD    = MULD2H(ISYMB,ISYMD)
         ISYLK    = MULD2H(ISYBD,ISYMT2)
         ISYFK    = MULD2H(ISYFKY,ISYLK)
C
         KLK      = KEND0
         KFK      = KLK  + NMATIJ(ISYLK) 
         KEND1    = KFK  + NT1AM(ISYFK)
         LWRK1    = LWRK  - KEND1
C
         CALL DZERO(WRK(KFK),NT1AM(ISYFK))
         CALL DZERO(WRK(KW),NCKIJ(ISWMAT))
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWRK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND1
            CALL QUIT('Insufficient space in WXBD_T2 (7)')
         END IF
C
         CALL  SORT_T2_IJ(WRK(KLK),ISYMB,B,ISYMD,D,T2TP,ISYMT2)
C
C                  work(fk) = focky(lf) * I^BD(lk)

         DO ISYMK = 1,NSYM
            ISYML = MULD2H(ISYMK,ISYLK)
            ISYMF = MULD2H(ISYFKY,ISYML)
C
            KOFF1  = KFCLF + IT1AMT(ISYML,ISYMF)  
            KOFF2  = KLK   + IMATIJ(ISYML,ISYMK)
            KOFF3 =  KFK   + IT1AM(ISYMF,ISYMK)
C
            NF = MAX(1,NVIR(ISYMF))
            NL = MAX(1,NRHF(ISYML))        
C
            CALL DGEMM('T','N',NVIR(ISYMF),NRHF(ISYMK),
     *                 NRHF(ISYML),ONE,WRK(KOFF1),NL,
     *                 WRK(KOFF2),NL,ONE,WRK(KOFF3),NF)
C
         END DO
C
C                  work(aijk) = I(aijf) * work(fk)
C
         DO ISYMK = 1,NSYM
            ISYMF = MULD2H(ISYMK,ISYFK)
            ISYAIJ = MULD2H(ISYMT2,ISYMF)
C
            KOFF1 =  IT2SP(ISYAIJ,ISYMF) + 1 
            KOFF2 =  KFK +  IT1AM(ISYMF,ISYMK)
            KOFF3 =  KW  +  ISAIKJ(ISYAIJ,ISYMK) 
C
            NAIJ = MAX(1,NCKI(ISYAIJ))
            NF = MAX(1,NVIR(ISYMF))

            CALL DGEMM('N','N',NCKI(ISYAIJ),NRHF(ISYMK),
     *                 NVIR(ISYMF),-ONE,T2TP(KOFF1),NAIJ,
     *                 WRK(KOFF2),NF,ONE,WRK(KOFF3),NAIJ)
C
         END DO
C
C     change order aijk to aikj
C
         DO I = 1,NCKIJ(ISWMAT)
            WMAT(I) = WMAT(I) + WRK(INDSQ(I,3))
         END DO
C
      END IF
C

      CALL QEXIT('WXBDT2')
C
      RETURN
      END
C  /* Deck wxbd_ground */
      SUBROUTINE WXBD_GROUND(AIBJCK_PERM,
     *                       T2TP,ISYMT2,TMAT,
     *                       TRVIR,TRVIR2,TRVIR3,TRVIR4,
     *                       TROCC,ISYINT,
     *                       WMAT,WORK,LWORK,
     *                       INDSQ,LENSQ,
     *                       ISYMB,B,ISYMC,C)
C
C     Kasper Hald, Summer 2002
C
C     Modified: Filip Pawlowski, Spring 2003. (AIBJCK_PERM)
C
C
C     ISYINT is the symmetry of the integrals.
C     ISYMT2 is the symmetry of T2TP.
C
C
C
C If (AIBJCK_PERM.eq.1) then (aibjck) + (aickbj)  permutation 
C
C     WBD(ai,k,j) = WBD(ai,k,j)
C                 + t(ai,dj) (Ck|Bd) 
C                 + t(ai,dj) (Bk|Cd)   ! right array is obtained using INDSQ(3)
C                                      ! after DGEMM
C                 - t(ai,Bl) (Ck|lj) 
C                 - t(ai,Cl) (Bk|lj)   ! right array is obtained using INDSQ(3)
C                                      ! after DGEMM
C
C else if (AIBJCK_PERM.eq.2) then (bjckai) + (bjaick)  permutation
C
C     WBD(ai,k,j) = WBD(ai,k,j)
C                 + t(Bj,dk) (ai|Cd) 
C                 + t(Bj,di) (Ck|ad)   
C                 - t(Bj,Cl) (ai|lk) 
C                 - t(Bj,al) (Ck|li)   
C
C else if (AIBJCK_PERM.eq.3) then (ckbjai) + (ckaibj)  permutation
C
C     WBD(ai,k,j) = WBD(ai,k,j)
C                 + t(Ck,dj) (ai|Bd) 
C                 + t(Ck,di) (Bj|ad)   
C                 - t(Ck,Bl) (ai|lj) 
C                 - t(Ck,al) (Bj|li)   
C
C else if (AIBJCK_PERM.eq.4) then calculate all terms:
C
C      WBD(ai,k,j) = (aibjck) + (aickbj) + (bjckai) + (bjaick)
C                  + (ckbjai) + (ckaibj)
C
C
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C#include <cclr.h>
C
      INTEGER AIBJCK_PERM
C
      INTEGER ISYMT2, ISYINT, LWORK, LENSQ, ISYMB, ISYMC
      INTEGER INDSQ(LENSQ,6)
      INTEGER ISYMBC, ISYRES, JSAIKJ, ISYMDK, LENGTH, ISYMK
      INTEGER ISYMD, ISYAIJ, KOFF1, KOFF2, KOFF3, NTOAIJ, NVIRD
      INTEGER ISYAIL, ISYLKJ, ISYMJ, ISYMLK, ISYML, ISYMAI
      INTEGER ISYAIK, NTOTAI, NRHFL
C
      INTEGER ISYAID,ISYDJK,ISYAIJK,ISYMDJ
      INTEGER NTOTD
      INTEGER ISYAJD,ISYDIK,ISYAJIK,ISYMDI,ISYMI,ISYMAJ,ISYAJI
      INTEGER KAJD,KEND1,LWRK1
      INTEGER NTOTAJ
      INTEGER ISYAJIL,ISYBLK,ISYMBL
      INTEGER KLK,KAJIL
      INTEGER NTOTAJI,NTOTL
      INTEGER ISYALK,ISYLJI,ISYMAL,ISYMLJ,ISYMA
      INTEGER NTOTA
      INTEGER ISYDKJ
      INTEGER ISYAKD,ISYDIJ,ISYAKIJ,KAKD
      INTEGER ISYMAK,ISYAKI,NTOTAK,ISYAKIL,ISYMCB,KLJ,KAKIL,NTOTAKI
      INTEGER ISYALJ,ISYLKI
C
#if defined (SYS_CRAY)
      REAL T2TP(*), TRVIR(*), TRVIR2(*), TRVIR3(*), TROCC(*)
      REAL TRVIR4(*)
      REAL TMAT(*), WMAT(*), WORK(LWORK)
      REAL ZERO, ONE, DDOT, XNORM
#else
      DOUBLE PRECISION T2TP(*), TRVIR(*), TRVIR2(*), TRVIR3(*), TROCC(*)
      DOUBLE PRECISION TRVIR4(*)
      DOUBLE PRECISION TMAT(*), WMAT(*), WORK(LWORK)
      DOUBLE PRECISION ZERO, ONE, DDOT, XNORM
#endif
C
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
C
      CALL QENTER('WXBDGR')
C
C---------------------------------------
C     Initial test of AIBJCK_PERM option
C---------------------------------------
C
      IF ( (AIBJCK_PERM .LT. 1) .OR. (AIBJCK_PERM .GT. 4) ) THEN
         WRITE(LUPRI,*)'AIBJCK_PERM = ',AIBJCK_PERM
         WRITE(LUPRI,*)'should be between 1 and 4 '
         CALL QUIT('Illegal value of AIBJCK_PERM option in WXBD_GROUND')
      END IF
C
C************************
C************************
      IF ((AIBJCK_PERM.EQ.1) .OR. (AIBJCK_PERM.EQ.4)) THEN 
C************************
C************************

C
C--------------------------------
C        First virtual contribution.
C--------------------------------
C
         ISYMBC = MULD2H(ISYMB,ISYMC)
         ISYRES = MULD2H(ISYINT,ISYMT2)
         JSAIKJ = MULD2H(ISYMBC,ISYRES)
         ISYMDK = MULD2H(ISYMBC,ISYINT)
C
         LENGTH = NCKIJ(JSAIKJ)
C
         IF (LWORK .LT. LENGTH) THEN
            CALL QUIT('Insufficient core in WXBD_GROUND (1)')
         ENDIF
C
         DO ISYMK = 1,NSYM
C
            ISYMD  = MULD2H(ISYMK,ISYMDK)
            ISYAIJ = MULD2H(ISYMK,JSAIKJ)
C
            KOFF1 = IT2SP(ISYAIJ,ISYMD)  + 1
            KOFF2 = ICKATR(ISYMDK,ISYMB) + NT1AM(ISYMDK)*(B - 1)
     *            + IT1AM(ISYMD,ISYMK)   + 1
            KOFF3 = ISAIKJ(ISYAIJ,ISYMK) + 1
C
            NTOAIJ = MAX(1,NCKI(ISYAIJ))
            NVIRD  = MAX(NVIR(ISYMD),1)
C
            CALL DGEMM('N','N',NCKI(ISYAIJ),NRHF(ISYMK),
     *                 NVIR(ISYMD),ONE,T2TP(KOFF1),NTOAIJ,
     *                 TRVIR(KOFF2),NVIRD,ZERO,
     *                 WORK(KOFF3),NTOAIJ)
C
         ENDDO
C
c add_wbdx_1v
         DO I = 1,LENGTH
            WMAT(I) = WMAT(I) + WORK(INDSQ(I,3))
         ENDDO
C
         IF (IPRINT .GT. 55) THEN
            XNORM = DDOT(NCKIJ(JSAIKJ),WMAT,1,WMAT,1)
            WRITE(LUPRI,*) 'In WXBD_GROUND: 1. Norm of WMAT ',XNORM
         ENDIF
C
C---------------------------------
C        Second virtual contribution.
C---------------------------------
C
         ISYMBC = MULD2H(ISYMB,ISYMC)
         ISYRES = MULD2H(ISYINT,ISYMT2)
         JSAIKJ = MULD2H(ISYMBC,ISYRES)
         ISYMDK = MULD2H(ISYMBC,ISYINT)
C
         LENGTH = NCKIJ(JSAIKJ)
C
         IF (LWORK .LT. LENGTH) THEN
            CALL QUIT('Insufficient core in WXBD_GROUND (2)')
         ENDIF
C
         DO ISYMK = 1,NSYM
C
            ISYMD  = MULD2H(ISYMK,ISYMDK)
            ISYAIJ = MULD2H(ISYMK,JSAIKJ)
C
            KOFF1 = IT2SP(ISYAIJ,ISYMD)  + 1
            KOFF2 = ICKATR(ISYMDK,ISYMC) + NT1AM(ISYMDK)*(C - 1)
     *            + IT1AM(ISYMD,ISYMK)   + 1
            KOFF3 = ISAIKJ(ISYAIJ,ISYMK) + 1
C
            NTOAIJ = MAX(1,NCKI(ISYAIJ))
            NVIRD  = MAX(NVIR(ISYMD),1)
C
            CALL DGEMM('N','N',NCKI(ISYAIJ),NRHF(ISYMK),
     *                 NVIR(ISYMD),ONE,T2TP(KOFF1),NTOAIJ,
     *                 TRVIR2(KOFF2),NVIRD,ZERO,
     *                 WORK(KOFF3),NTOAIJ)
C
         ENDDO
C
c add_wbdx_2v
         DO I = 1,LENGTH
            WMAT(I) = WMAT(I) + WORK(I)
         ENDDO
C
         IF (IPRINT .GT. 55) THEN
            XNORM = DDOT(NCKIJ(JSAIKJ),WMAT,1,WMAT,1)
            WRITE(LUPRI,*) 'In WXBD_GROUND: 2. Norm of WMAT ',XNORM
         ENDIF
C
C---------------------------------
C        First occupied contribution.
C---------------------------------
C
         ISYAIL = MULD2H(ISYMB,ISYMT2)
         ISYLKJ = MULD2H(ISYMC,ISYINT)
C
         DO ISYMJ = 1,NSYM
C
            ISYMLK = MULD2H(ISYMJ,ISYLKJ)
C
            DO J = 1,NRHF(ISYMJ)
C
               DO ISYMK = 1,NSYM
C
                  ISYML  = MULD2H(ISYMK,ISYMLK)
                  ISYMAI = MULD2H(ISYAIL,ISYML)
                  ISYAIK = MULD2H(ISYMAI,ISYMK)
C
                  KOFF1 = IT2SP(ISYAIL,ISYMB)
     *                  + NCKI(ISYAIL)*(B - 1)
     *                  + ICKI(ISYMAI,ISYML) + 1
                  KOFF2 = ISJIKA(ISYLKJ,ISYMC)
     *                  + NMAJIK(ISYLKJ)*(C - 1)
     *                  + ISJIK(ISYMLK,ISYMJ)
     *                  + NMATIJ(ISYMLK)*(J - 1)
     *                  + IMATIJ(ISYML,ISYMK) + 1
                  KOFF3 = ISAIKJ(ISYAIK,ISYMJ)
     *                  + NCKI(ISYAIK)*(J - 1)
     *                  + ICKI(ISYMAI,ISYMK) + 1
C
                  NTOTAI = MAX(1,NT1AM(ISYMAI))
                  NRHFL  = MAX(1,NRHF(ISYML))
C
c add_wbdx_1o
                  CALL DGEMM('N','N',NT1AM(ISYMAI),NRHF(ISYMK),
     *                       NRHF(ISYML),-ONE,T2TP(KOFF1),NTOTAI,
     *                       TROCC(KOFF2),NRHFL,ONE,WMAT(KOFF3),
     *                       NTOTAI)
C
               ENDDO
            ENDDO
         ENDDO
C
         IF (IPRINT .GT. 55) THEN
            XNORM = DDOT(NCKIJ(JSAIKJ),WMAT,1,WMAT,1)
            WRITE(LUPRI,*) 'In WXBD_GROUND: 3. Norm of WMAT ',XNORM
         ENDIF
C
C----------------------------------
C        Second occupied contribution.
C----------------------------------
C
         ISYAIL = MULD2H(ISYMC,ISYMT2)
         ISYLKJ = MULD2H(ISYMB,ISYINT)
C
         DO ISYMJ = 1,NSYM
C
            ISYMLK = MULD2H(ISYMJ,ISYLKJ)
C
            DO J = 1,NRHF(ISYMJ)
C
               DO ISYMK = 1,NSYM
C
                  ISYML  = MULD2H(ISYMK,ISYMLK)
                  ISYMAI = MULD2H(ISYAIL,ISYML)
                  ISYAIK = MULD2H(ISYMAI,ISYMK)
C
                  KOFF1 = IT2SP(ISYAIL,ISYMC)
     *                  + NCKI(ISYAIL)*(C - 1)
     *                  + ICKI(ISYMAI,ISYML) + 1
                  KOFF2 = ISJIKA(ISYLKJ,ISYMB)
     *                  + NMAJIK(ISYLKJ)*(B - 1)
     *                  + ISJIK(ISYMLK,ISYMJ)
     *                  + NMATIJ(ISYMLK)*(J - 1)
     *                  + IMATIJ(ISYML,ISYMK) + 1
                  KOFF3 = ISAIKJ(ISYAIK,ISYMJ)
     *                  + NCKI(ISYAIK)*(J - 1)
     *                  + ICKI(ISYMAI,ISYMK) + 1
C
                  NTOTAI = MAX(1,NT1AM(ISYMAI))
                  NRHFL  = MAX(1,NRHF(ISYML))
C
                  CALL DGEMM('N','N',NT1AM(ISYMAI),NRHF(ISYMK),
     *                       NRHF(ISYML),-ONE,T2TP(KOFF1),NTOTAI,
     *                       TROCC(KOFF2),NRHFL,ZERO,TMAT(KOFF3),
     *                       NTOTAI)
C
               ENDDO
            ENDDO
         ENDDO
C
c add_wbdx_2o
         DO I = 1,LENGTH
            WMAT(I) = WMAT(I) + TMAT(INDSQ(I,3))
         ENDDO
C
         IF (IPRINT .GT. 55) THEN
            XNORM = DDOT(NCKIJ(JSAIKJ),WMAT,1,WMAT,1)
            WRITE(LUPRI,*) 'In WXBD_GROUND: 4. Norm of WMAT ',XNORM
         ENDIF
C
C**********
C**********
      END IF
      IF ( (AIBJCK_PERM .EQ. 2) .OR. (AIBJCK_PERM .EQ. 4)) THEN
C**********
C**********
C    (bjckai) + (bjaick)  permutation
C
C     WBD(ai,k,j) = WBD(ai,k,j)
C                 + t(Bj,dk) (ai|Cd) 
C                 + t(Bj,di) (Ck|ad)   
C                 - t(Bj,Cl) (ai|lk) 
C                 - t(Bj,al) (Ck|li)   
C-----------------------------------
C        First virtual contribution
C-----------------------------------
C
C WBD(ai,k,j) = WBD(ai,k,j) + t(Cj,dk) (ai|Cd)
C
C                         tmat(aikj) = TRVIR4^C(aid)*T2TP^B(dkj)
C
C WBD(ai,k,j) = WBD(ai,k,j) + tmat(aikj)
C
         ISYAID  = MULD2H(ISYINT,ISYMC)
         ISYDKJ  = MULD2H(ISYMT2,ISYMB)
C
C--------------------------------------------------------
C        tmat(aikj) = TRVIR4^C(aid)*T2TP^B(dkj)
C--------------------------------------------------------
C
         DO ISYMJ = 1,NSYM
            ISYMDK = MULD2H(ISYDKJ,ISYMJ)
            DO ISYMD = 1,NSYM
               ISYMK  = MULD2H(ISYMDK,ISYMD)
               ISYMAI = MULD2H(ISYAID,ISYMD)
               ISYAIK = MULD2H(ISYMAI,ISYMK)
               DO J = 1,NRHF(ISYMJ)
C
                  KOFF1 = ICKATR(ISYMAI,ISYMD)
     *                  + 1
                  KOFF2 = IT2SP(ISYDKJ,ISYMB)
     *                  + NCKI(ISYDKJ)*(B-1)
     *                  + ICKI(ISYMDK,ISYMJ)
     *                  + NT1AM(ISYMDK)*(J-1)
     *                  + IT1AM(ISYMD,ISYMK)
     *                  + 1
                  KOFF3 = ISAIKJ(ISYAIK,ISYMJ)
     *                  + NCKI(ISYAIK)*(J-1)
     *                  + ICKI(ISYMAI,ISYMK)
     *                  + 1
C
                  NTOTAI = MAX(NT1AM(ISYMAI),1)
                  NTOTD  = MAX(NVIR(ISYMD),1)
C
                  CALL DGEMM('N','N',NT1AM(ISYMAI),NRHF(ISYMK),
     *                       NVIR(ISYMD),ONE,TRVIR4(KOFF1),NTOTAI,
     *                       T2TP(KOFF2),NTOTD,ONE,WMAT(KOFF3),
     *                       NTOTAI)
C
               END DO
            END DO
         END DO
C
C-----------------------------------
C        Second virtual contribution
C-----------------------------------
C
C WBD(ai,k,j) = WBD(ai,k,j) + t(Bj,di) (Ck|ad)
C
C                         tmat(akij) = work(akd)*T2TP^B(dij)
C
C WBD(ai,k,j) = WBD(ai,k,j) + tmat(akij)
C
         ISYAKD  = MULD2H(ISYINT,ISYMC)
         ISYDIJ  = MULD2H(ISYMT2,ISYMB)
         ISYAKIJ = MULD2H(ISYAKD,ISYDIJ)
C
         KAKD  = 1
         KEND1 = KAKD  + NCKATR(ISYAKD)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*)'Memory available : ', LWORK
            WRITE(LUPRI,*)'Memory needed    : ', KEND1
            CALL QUIT('Insufficient core in WXBD_GROUND (3)')
         ENDIF
C
         CALL DZERO(TMAT,NCKIJ(ISYAKIJ))
C
C-------------------------------------------------------
C        Sort (Ck|ad) integrals sitting as TRVIR^C(dka)
C        to WORK(akd)
C-------------------------------------------------------
C
         CALL CCSDT_SRTVIR(TRVIR,WORK(KAKD),WORK(KEND1),LWRK1,
     *                     ISYMC,ISYINT)
C
C--------------------------------------------------------
C        tmat(akij) = work(akd)*T2TP^B(dij)
C--------------------------------------------------------
C
         DO ISYMJ = 1,NSYM
            ISYMDI = MULD2H(ISYDIJ,ISYMJ)
            DO ISYMD = 1,NSYM
               ISYMI  = MULD2H(ISYMDI,ISYMD)
               ISYMAK = MULD2H(ISYAKD,ISYMD)
               ISYAKI = MULD2H(ISYMAK,ISYMI)
               DO J = 1,NRHF(ISYMJ)
C
                  KOFF1 = KAKD
     *                  + ICKATR(ISYMAK,ISYMD)
                  KOFF2 = IT2SP(ISYDIJ,ISYMB)
     *                  + NCKI(ISYDIJ)*(B-1)
     *                  + ICKI(ISYMDI,ISYMJ)
     *                  + NT1AM(ISYMDI)*(J-1)
     *                  + IT1AM(ISYMD,ISYMI)
     *                  + 1
                  KOFF3 = ISAIKJ(ISYAKI,ISYMJ)
     *                  + NCKI(ISYAKI)*(J-1)
     *                  + ICKI(ISYMAK,ISYMI)
     *                  + 1
C
                  NTOTAK = MAX(NT1AM(ISYMAK),1)
                  NTOTD  = MAX(NVIR(ISYMD),1)
C
                  CALL DGEMM('N','N',NT1AM(ISYMAK),NRHF(ISYMI),
     *                       NVIR(ISYMD),ONE,WORK(KOFF1),NTOTAK,
     *                       T2TP(KOFF2),NTOTD,ONE,TMAT(KOFF3),
     *                       NTOTAK)
C
               END DO
            END DO
         END DO
C
C----------------------------------------------
C        WBD(ai,k,j) = WBD(ai,k,j) + tmat(akij)
C----------------------------------------------
C
c add_wbdx_2va
         DO I = 1,NCKIJ(ISYAKIJ)
            WMAT(I) = WMAT(I) + TMAT(INDSQ(I,1))
         ENDDO
C
C-----------------------------------
C        First occupied contribution
C-----------------------------------
C
C WBD(ai,k,j) = WBD(ai,k,j) - t(Bj,Cl) (ai|lk)
C
C               tmat(akij) = -         work(akil) * work(lj) 
C
C WBD(ai,k,j) = WBD(ai,k,j) + tmat(akij)
C
         ISYAKIL = ISYINT
         ISYMCB  = MULD2H(ISYMC,ISYMB)
         ISYMLJ  = MULD2H(ISYMT2,ISYMCB)
         ISYAKIJ = MULD2H(ISYAKIL,ISYMLJ)
C
         KLJ   = 1
         KAKIL = KLJ   + NMATIJ(ISYMLJ)
         KEND1 = KAKIL + NTRAOC(ISYAKIL)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*)'Memory available : ', LWORK
            WRITE(LUPRI,*)'Memory needed    : ', KEND1
            CALL QUIT('Insufficient core in WXBD_GROUND (4)')
         ENDIF
C
         CALL DZERO(TMAT,NCKIJ(ISYAKIJ))
C
C------------------------------------------------------
C        Sort integrals (ai|lk) sitting as TROCC(lik,a)
C        to WORK(aki,l)
C------------------------------------------------------
C
         CALL CCFOP_SORT(TROCC,WORK(KAKIL),ISYINT,1)
C
C-------------------------------------------------------
C        Sort T2 amplitudes from T2TP(Cl,jB) to WORK(lj)
C-------------------------------------------------------
C
         CALL SORT_T2_IJ(WORK(KLJ),ISYMC,C,ISYMB,B,T2TP,ISYMT2)
C
C--------------------------------------------
C        tmat(akij) = - work(akil) * work(lj) 
C--------------------------------------------
C
         DO ISYML = 1,NSYM
            ISYAKI = MULD2H(ISYAKIL,ISYML)
            ISYMJ  = MULD2H(ISYMLJ,ISYML)
C
            KOFF1 = KAKIL
     *            + ISAIKJ(ISYAKI,ISYML)
            KOFF2 = KLJ
     *            + IMATIJ(ISYML,ISYMJ)
            KOFF3 = ISAIKJ(ISYAKI,ISYMJ)
     *            + 1
C
            NTOTAKI = MAX(NCKI(ISYAKI),1)
            NTOTL   = MAX(NRHF(ISYML),1)
C
            CALL DGEMM('N','N',NCKI(ISYAKI),NRHF(ISYMJ),NRHF(ISYML),
     *                 -ONE,WORK(KOFF1),NTOTAKI,WORK(KOFF2),NTOTL,
     *                 ONE,TMAT(KOFF3),NTOTAKI)
C
         END DO
C
C----------------------------------------------
C        WBD(ai,k,j) = WBD(ai,k,j) + tmat(akij)
C----------------------------------------------
C
c add_wbdx_1oa
         DO I = 1,NCKIJ(ISYAKIJ)
            WMAT(I) = WMAT(I) + TMAT(INDSQ(I,1))
         ENDDO 
C
C-----------------------------------
C        Second occupied contribution
C-----------------------------------
C
C WBD(ai,k,j) = WBD(ai,k,j) - t(Bj,al) (Ck|li)
C
C               tmat(akij) = - T2TP^B(alj) * TROCC(lki,C)
C
C WBD(ai,k,j) = WBD(ai,k,j) + tmat(akij)
C
         ISYALJ  = MULD2H(ISYMT2,ISYMB)
         ISYLKI  = MULD2H(ISYINT,ISYMC)
         ISYAKIJ = MULD2H(ISYALJ,ISYLKI)
C
         CALL DZERO(TMAT,NCKIJ(ISYAKIJ))
C
         DO ISYMJ = 1,NSYM
            ISYMAL = MULD2H(ISYALJ,ISYMJ)
            ISYAKI = MULD2H(ISYAKIJ,ISYMJ)
            DO ISYMI = 1,NSYM
               ISYMLK = MULD2H(ISYLKI,ISYMI)
               ISYMAK = MULD2H(ISYAKI,ISYMI)
               DO ISYML = 1,NSYM
                  ISYMA = MULD2H(ISYMAL,ISYML)
                  ISYMK = MULD2H(ISYMLK,ISYML)
                  DO J = 1,NRHF(ISYMJ)
                     DO I = 1,NRHF(ISYMI)
C
                        KOFF1 = IT2SP(ISYALJ,ISYMB)
     *                        + NCKI(ISYALJ)*(B-1)
     *                        + ICKI(ISYMAL,ISYMJ)
     *                        + NT1AM(ISYMAL)*(J-1)
     *                        + IT1AM(ISYMA,ISYML)
     *                        + 1
                        KOFF2 = ISJIKA(ISYLKI,ISYMC)
     *                        + NMAJIK(ISYLKI)*(C-1)
     *                        + ISJIK(ISYMLK,ISYMI)
     *                        + NMATIJ(ISYMLK)*(I-1)
     *                        + IMATIJ(ISYML,ISYMK)
     *                        + 1
                        KOFF3 = ISAIKJ(ISYAKI,ISYMJ)
     *                        + NCKI(ISYAKI)*(J-1)
     *                        + ICKI(ISYMAK,ISYMI)
     *                        + NT1AM(ISYMAK)*(I-1)
     *                        + IT1AM(ISYMA,ISYMK)
     *                        + 1
C
                        NTOTA = MAX(NVIR(ISYMA),1)
                        NTOTL = MAX(NRHF(ISYML),1)
C
                        CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMK),
     *                             NRHF(ISYML),-ONE,T2TP(KOFF1),
     *                             NTOTA,TROCC(KOFF2),NTOTL,ONE,
     *                             TMAT(KOFF3),NTOTA)
C
                     END DO
                  END DO
               END DO
            END DO
         END DO
C
C----------------------------------------------
C        WBD(ai,k,j) = WBD(ai,k,j) + tmat(akij)
C----------------------------------------------
C
c add_wbdx_2oa
         DO I = 1,NCKIJ(ISYAKIJ)
            WMAT(I) = WMAT(I) + TMAT(INDSQ(I,1))
         ENDDO
C
C**********
C**********
      END IF
      IF ( (AIBJCK_PERM .EQ. 3) .OR. (AIBJCK_PERM .EQ. 4)) THEN
C**********
C**********
C    (ckbjai) + (ckaibj)  permutation
C
C     WBD(ai,k,j) = WBD(ai,k,j)
C                 + t(Ck,dj) (ai|Bd) 
C                 + t(Ck,di) (Bj|ad)   
C                 - t(Ck,Bl) (ai|lj) 
C                 - t(Ck,al) (Bj|li)   
C-----------------------------------
C        First virtual contribution
C-----------------------------------
C
C WBD(ai,k,j) = WBD(ai,k,j) + t(Ck,dj) (ai|Bd)
C
C                         tmat(aijk) = TRVIR3^B(aid)*T2TP^C(djk)
C
C WBD(ai,k,j) = WBD(ai,k,j) + tmat(aijk)
C
         ISYAID  = MULD2H(ISYINT,ISYMB)
         ISYDJK  = MULD2H(ISYMT2,ISYMC)
         ISYAIJK = MULD2H(ISYAID,ISYDJK)
C
         CALL DZERO(TMAT,NCKIJ(ISYAIJK))
C
C--------------------------------------------------------
C        tmat(aijk) = TRVIR3^B(aid)*T2TP^C(djk)
C--------------------------------------------------------
C
         DO ISYMK = 1,NSYM
            ISYMDJ = MULD2H(ISYDJK,ISYMK)
            DO ISYMD = 1,NSYM
               ISYMJ  = MULD2H(ISYMDJ,ISYMD)
               ISYMAI = MULD2H(ISYAID,ISYMD)
               ISYAIJ = MULD2H(ISYMAI,ISYMJ)
               DO K = 1,NRHF(ISYMK)
C
                  KOFF1 = ICKATR(ISYMAI,ISYMD)
     *                  + 1
                  KOFF2 = IT2SP(ISYDJK,ISYMC)
     *                  + NCKI(ISYDJK)*(C-1)
     *                  + ICKI(ISYMDJ,ISYMK)
     *                  + NT1AM(ISYMDJ)*(K-1)
     *                  + IT1AM(ISYMD,ISYMJ)
     *                  + 1
                  KOFF3 = ISAIKJ(ISYAIJ,ISYMK)
     *                  + NCKI(ISYAIJ)*(K-1)
     *                  + ICKI(ISYMAI,ISYMJ)
     *                  + 1
C
                  NTOTAI = MAX(NT1AM(ISYMAI),1)
                  NTOTD  = MAX(NVIR(ISYMD),1)
C
                  CALL DGEMM('N','N',NT1AM(ISYMAI),NRHF(ISYMJ),
     *                       NVIR(ISYMD),ONE,TRVIR3(KOFF1),NTOTAI,
     *                       T2TP(KOFF2),NTOTD,ONE,TMAT(KOFF3),
     *                       NTOTAI)
C
               END DO
            END DO
         END DO
C
C----------------------------------------------
C        WBD(ai,k,j) = WBD(ai,k,j) + tmat(aijk)
C----------------------------------------------
C
c add_wbdx_1va
         DO I = 1,NCKIJ(ISYAIJK)
            WMAT(I) = WMAT(I) + TMAT(INDSQ(I,3))
         ENDDO
C
C-----------------------------------
C        Second virtual contribution
C-----------------------------------
C
C WBD(ai,k,j) = WBD(ai,k,j) + t(Ck,di) (Bj|ad)
C
C                         tmat(ajik) = work(ajd)*T2TP^C(dik)
C
C WBD(ai,k,j) = WBD(ai,k,j) + tmat(ajik)
C
         ISYAJD  = MULD2H(ISYINT,ISYMB)
         ISYDIK  = MULD2H(ISYMT2,ISYMC)
         ISYAJIK = MULD2H(ISYAJD,ISYDIK)
C
         KAJD  = 1
         KEND1 = KAJD  + NCKATR(ISYAJD)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*)'Memory available : ', LWORK
            WRITE(LUPRI,*)'Memory needed    : ', KEND1
            CALL QUIT('Insufficient core in WXBD_GROUND (5)')
         ENDIF
C
         CALL DZERO(TMAT,NCKIJ(ISYAJIK))
C
C-------------------------------------------------------
C        Sort (Bj|ad) integrals sitting as TRVIR2^B(dja)
C        to WORK(ajd)
C-------------------------------------------------------
C
         CALL CCSDT_SRTVIR(TRVIR2,WORK(KAJD),WORK(KEND1),LWRK1,
     *                     ISYMB,ISYINT)
C
C--------------------------------------------------------
C        tmat(ajik) = work(ajd)*T2TP^C(dik)
C--------------------------------------------------------
C
         DO ISYMK = 1,NSYM
            ISYMDI = MULD2H(ISYDIK,ISYMK)
            DO ISYMD = 1,NSYM
               ISYMI  = MULD2H(ISYMDI,ISYMD)
               ISYMAJ = MULD2H(ISYAJD,ISYMD)
               ISYAJI = MULD2H(ISYMAJ,ISYMI)
               DO K = 1,NRHF(ISYMK)
C
                  KOFF1 = KAJD
     *                  + ICKATR(ISYMAJ,ISYMD)
                  KOFF2 = IT2SP(ISYDIK,ISYMC)
     *                  + NCKI(ISYDIK)*(C-1)
     *                  + ICKI(ISYMDI,ISYMK)
     *                  + NT1AM(ISYMDI)*(K-1)
     *                  + IT1AM(ISYMD,ISYMI)
     *                  + 1
                  KOFF3 = ISAIKJ(ISYAJI,ISYMK)
     *                  + NCKI(ISYAJI)*(K-1)
     *                  + ICKI(ISYMAJ,ISYMI)
     *                  + 1
C
                  NTOTAJ = MAX(NT1AM(ISYMAJ),1)
                  NTOTD  = MAX(NVIR(ISYMD),1)
C
                  CALL DGEMM('N','N',NT1AM(ISYMAJ),NRHF(ISYMI),
     *                       NVIR(ISYMD),ONE,WORK(KOFF1),NTOTAJ,
     *                       T2TP(KOFF2),NTOTD,ONE,TMAT(KOFF3),
     *                       NTOTAJ)
C
               END DO
            END DO
         END DO
C
C----------------------------------------------
C        WBD(ai,k,j) = WBD(ai,k,j) + tmat(ajik)
C----------------------------------------------
C
c add_wbdx_2va
         DO I = 1,NCKIJ(ISYAJIK)
            WMAT(I) = WMAT(I) + TMAT(INDSQ(I,4))
         ENDDO
C
C-----------------------------------
C        First occupied contribution
C-----------------------------------
C
C WBD(ai,k,j) = WBD(ai,k,j) - t(Ck,Bl) (ai|lj)
C
C               tmat(ajik) = -         work(ajil) * work(lk) 
C
C WBD(ai,k,j) = WBD(ai,k,j) + tmat(ajik)
C
         ISYAJIL = ISYINT
         ISYMBC  = MULD2H(ISYMB,ISYMC)
         ISYMLK  = MULD2H(ISYMT2,ISYMBC)
         ISYAJIK = MULD2H(ISYAJIL,ISYMLK)
C
         KLK   = 1
         KAJIL = KLK   + NMATIJ(ISYMLK)
         KEND1 = KAJIL + NTRAOC(ISYAJIL)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*)'Memory available : ', LWORK
            WRITE(LUPRI,*)'Memory needed    : ', KEND1
            CALL QUIT('Insufficient core in WXBD_GROUND (6)')
         ENDIF
C
         CALL DZERO(TMAT,NCKIJ(ISYAJIK))
C
C------------------------------------------------------
C        Sort integrals (ai|lj) sitting as TROCC(lij,a)
C        to WORK(aji,l)
C------------------------------------------------------
C
         CALL CCFOP_SORT(TROCC,WORK(KAJIL),ISYINT,1)
C
C-------------------------------------------------------
C        Sort T2 amplitudes from T2TP(Bl,kC) to WORK(lk)
C-------------------------------------------------------
C
         CALL SORT_T2_IJ(WORK(KLK),ISYMB,B,ISYMC,C,T2TP,ISYMT2)
C
C--------------------------------------------
C        tmat(ajik) = - work(ajil) * work(lk) 
C--------------------------------------------
C
         DO ISYML = 1,NSYM
            ISYAJI = MULD2H(ISYAJIL,ISYML)
            ISYMK  = MULD2H(ISYMLK,ISYML)
C
            KOFF1 = KAJIL
     *            + ISAIKJ(ISYAJI,ISYML)
            KOFF2 = KLK
     *            + IMATIJ(ISYML,ISYMK)
            KOFF3 = ISAIKJ(ISYAJI,ISYMK)
     *            + 1
C
            NTOTAJI = MAX(NCKI(ISYAJI),1)
            NTOTL   = MAX(NRHF(ISYML),1)
C
            CALL DGEMM('N','N',NCKI(ISYAJI),NRHF(ISYMK),NRHF(ISYML),
     *                 -ONE,WORK(KOFF1),NTOTAJI,WORK(KOFF2),NTOTL,
     *                 ONE,TMAT(KOFF3),NTOTAJI)
C
         END DO
C
C----------------------------------------------
C        WBD(ai,k,j) = WBD(ai,k,j) + tmat(ajik)
C----------------------------------------------
C
c add_wbdx_1oa
         DO I = 1,NCKIJ(ISYAJIK)
            WMAT(I) = WMAT(I) + TMAT(INDSQ(I,4))
         ENDDO 
C
C-----------------------------------
C        Second occupied contribution
C-----------------------------------
C
C WBD(ai,k,j) = WBD(ai,k,j) - t(Ck,al) (Bj|li)
C
C               tmat(ajik) = - T2TP^C(alk) * TROCC(lji,B)
C
C WBD(ai,k,j) = WBD(ai,k,j) + tmat(ajik)
C
         ISYALK  = MULD2H(ISYMT2,ISYMC)
         ISYLJI  = MULD2H(ISYINT,ISYMB)
         ISYAJIK = MULD2H(ISYALK,ISYLJI)
C
         CALL DZERO(TMAT,NCKIJ(ISYAJIK))
C
         DO ISYMK = 1,NSYM
            ISYMAL = MULD2H(ISYALK,ISYMK)
            ISYAJI = MULD2H(ISYAJIK,ISYMK)
            DO ISYMI = 1,NSYM
               ISYMLJ = MULD2H(ISYLJI,ISYMI)
               ISYMAJ = MULD2H(ISYAJI,ISYMI)
               DO ISYML = 1,NSYM
                  ISYMA = MULD2H(ISYMAL,ISYML)
                  ISYMJ = MULD2H(ISYMLJ,ISYML)
                  DO K = 1,NRHF(ISYMK)
                     DO I = 1,NRHF(ISYMI)
C
                        KOFF1 = IT2SP(ISYALK,ISYMC)
     *                        + NCKI(ISYALK)*(C-1)
     *                        + ICKI(ISYMAL,ISYMK)
     *                        + NT1AM(ISYMAL)*(K-1)
     *                        + IT1AM(ISYMA,ISYML)
     *                        + 1
                        KOFF2 = ISJIKA(ISYLJI,ISYMB)
     *                        + NMAJIK(ISYLJI)*(B-1)
     *                        + ISJIK(ISYMLJ,ISYMI)
     *                        + NMATIJ(ISYMLJ)*(I-1)
     *                        + IMATIJ(ISYML,ISYMJ)
     *                        + 1
                        KOFF3 = ISAIKJ(ISYAJI,ISYMK)
     *                        + NCKI(ISYAJI)*(K-1)
     *                        + ICKI(ISYMAJ,ISYMI)
     *                        + NT1AM(ISYMAJ)*(I-1)
     *                        + IT1AM(ISYMA,ISYMJ)
     *                        + 1
C
                        NTOTA = MAX(NVIR(ISYMA),1)
                        NTOTL = MAX(NRHF(ISYML),1)
C
                        CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMJ),
     *                             NRHF(ISYML),-ONE,T2TP(KOFF1),
     *                             NTOTA,TROCC(KOFF2),NTOTL,ONE,
     *                             TMAT(KOFF3),NTOTA)
C
                     END DO
                  END DO
               END DO
            END DO
         END DO
C
C----------------------------------------------
C        WBD(ai,k,j) = WBD(ai,k,j) + tmat(ajik)
C----------------------------------------------
C
c add_wbdx_2oa
         DO I = 1,NCKIJ(ISYAJIK)
            WMAT(I) = WMAT(I) + TMAT(INDSQ(I,4))
         ENDDO
C
C**********
C**********
      END IF
C**********
C**********
C
C-----------
C     End.
C-----------
C
      CALL QEXIT('WXBDGR')
C
      RETURN
      END
C  /* Deck aden_dai_t2_d */
      SUBROUTINE ADEN_DAI_T2_D(DAI,ISYMDAI,T2TP,ISYMT2,
     *                         TETA,ISYMTETA,ISYMD,D,
     *                         ISYML,L,WORK,LWORK)
C
C     Calculate the contributions to Omega1 in CCSD(T) unrelaxed f.o.p.
C
C     DAI(ai) = DAI(ai) + (T{Dea}_{Lmi} - T{Dea}_{Lim}) * T2TP{emLD})
C
C      TETA contain on input TETA = T(Dea,Lmi) = T^DL(em,ai)
C      work(kemai) = T^DL(em,ai) - T^DL(ei,am) 
C
      IMPLICIT NONE
C
      INTEGER ISYMDAI,ISYMT2,ISYMTETA,ISYMD,ISYML,LWORK 
      INTEGER ISYMEM,ISYMEML,NTOTEM,KOFF1,KOFF2
      INTEGER KIJMJI,KEND1,LWRK1
C
#if defined (SYS_CRAY)
      REAL DAI(*), T2TP(*), TETA(*), WORK(LWORK), ONE, DM1
#else
      DOUBLE PRECISION DAI(*), T2TP(*), TETA(*), WORK(LWORK), ONE, DM1 
#endif
C
      PARAMETER ( ONE = 1.0D0, DM1 =-1.0D0)
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
C
      CALL QENTER('ADEN_DAI_T2_D')
C
C      TETA contain on input TETA = T(Dea,Lmi) = T^DL(em,ai)
C      work(kemai) = T^DL(em,ai) - T^DL(ei,am) 
C
      KIJMJI = 1
      KEND1  = KIJMJI + NT2SQ(ISYMTETA)
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in ADEN_DAI_T2_D')
      END IF
C
      CALL T2_AIBJ_EXC(WORK(KIJMJI),TETA,ISYMTETA,ONE,DM1)
C
C     DAI(ai) = DAI(ai) +  work(kemai) * T2TP{emLD}
C
      ISYMEM = MULD2H(ISYMDAI,ISYMTETA) 
      ISYMEML = MULD2H(ISYMEM,ISYML)
C
      NTOTEM = MAX(NT1AM(ISYMEM),1)
C
      KOFF1 = KIJMJI + IT2SQ(ISYMEM,ISYMDAI) 
C
      KOFF2 = 1 + IT2SP(ISYMEML,ISYMD)
     *       + NCKI(ISYMEML) *(D-1) 
     *       + ISAIK(ISYMEM,ISYML) 
     *       + NT1AM(ISYMEM)*(L-1) 
C 
c
      CALL DGEMV('T',NT1AM(ISYMEM),NT1AM(ISYMDAI),ONE,WORK(KOFF1),
     *           NTOTEM,T2TP(KOFF2),1,ONE,DAI,1)
c
C
      CALL QEXIT('ADEN_DAI_T2_D')
      RETURN
      END
C  /* Deck cc3_xi_den_ia */
      SUBROUTINE CC_XI_DEN_IA(T2XNET2Y,DIA,WMAT,ISWMAT,T2TPX,ISYMT2X,
     *                         T2TPY,ISYMT2Y,INDSQ,LENSQ,
     *                         IB,ISYMIB,ID,ISYMID,WORK,LWORK)
C
C            D(ia) = <L3Z|[[E_ia,T2X],T2Y]|HF>
C
C            D(ia) = D(ia)
C 
C                    - tbarZ^(dbc)_(jlk) t^x(db)_(ji) t^y(ca)_(kl)
C
C                    IF (T2X .NE. T2Y) THEN ADD
C
C                    - tbarZ^(dbc)_(jlk) t^y(db)_(ji) t^x(ca)_(kl)
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      LOGICAL T2XNET2Y
      INTEGER ISWMAT,ISYMT2X,ISYMT2Y,IB,ISYMIB,ID,ISYMID,LWORK
      INTEGER LENSQ,INDSQ(LENSQ,6)
C
#if defined (SYS_CRAY)
      REAL DIA(*),WMAT(*),T2TPX(*),T2TPY(*),WORK(LWORK)
#else
      DOUBLE PRECISION DIA(*),WMAT(*),T2TPX(*),T2TPY(*),WORK(LWORK)
#endif
C
      CALL QENTER('CC_XI_DEN_IA')
C
      CALL XI_DEN_IA(DIA,WMAT,ISWMAT,T2TPX,ISYMT2X,
     *               T2TPY,ISYMT2Y,INDSQ,LENSQ,
     *               IB,ISYMIB,ID,ISYMID,WORK,LWORK)

      IF (T2XNET2Y) THEN
         CALL XI_DEN_IA(DIA,WMAT,ISWMAT,T2TPY,ISYMT2Y,
     *                  T2TPX,ISYMT2X,INDSQ,LENSQ,
     *                  IB,ISYMIB,ID,ISYMID,WORK,LWORK)
      END IF
C
      CALL QEXIT('CC_XI_DEN_IA')
C
      RETURN
      END
C  /* Deck xi_den_ia */
      SUBROUTINE XI_DEN_IA(DIA,WMAT,ISWMAT,T2TPX,ISYMT2X,
     *                         T2TPY,ISYMT2Y,INDSQ,LENSQ,
     *                         IB,ISYMIB,ID,ISYMID,WORK,LWORK)
C
C=========================================================================
C    Calculate the density 
C
C
C            D(ia) = D(ia)
C 
C                  = - tbarZ^(dbc)_(jlk) t^x(db)_(ji) t^y(ca)_(kl) 
C
C                  = ( - W^BC(djkl) T2TPX(djiB) T2TPY(alkC)
C                      - W^DB(cklj) T2TPX(DjiB) T2TPY(ckla)
C                      - W^CD(bljk) T2TPX(bijD) T2TPY(alkC) )  
C
C            Density stored as AI:  DS(ai) = D(ia)
C
C            Outside the density DS(ai) is to be contracted with X_(ia)
C
C
C    Written by Filip Pawlowski, Fall 2002, Aarhus
C=========================================================================
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISWMAT,ISYMT2X,ISYMT2Y,IB,ISYMIB,ID,ISYMID,LWORK
      INTEGER ISYMB,ISYMC,ISYDJI,ISYALK,ISYKLI,ISYMDJ,ISYMKL,ISYMI
      INTEGER ISYMK,ISYMAL,ISYMA,ISYML,ISYMD,ISYMBD,ISYMJI,ISYCKLI
      INTEGER ISYMJ,ISYCKL,ISYBIJ,ISYLKI,ISYMBI,ISYMBJ,ISYMLK,ISYMAK
      INTEGER LENGTH,KTMAT,KEND1,LWRK1,KGKLI,KT2KLA
      INTEGER KT2JI,KGCKLI,KGLKI,KT2BJI,KEND2,LWRK2,KT2ALK
      INTEGER KOFF1,KOFF2,KOFF3
      INTEGER NTOTDJ,NTOTKL,NTOTA,NTOTJ,NTOTCKL,NTOTBJ,NTOTLK
      INTEGER LENSQ,INDSQ(LENSQ,6)
      INTEGER ISYMIBID,ISYDEN
C
#if defined (SYS_CRAY)
      REAL DIA(*),WMAT(*),T2TPX(*),T2TPY(*),WORK(LWORK),ONE,TWO
#else
      DOUBLE PRECISION DIA(*),WMAT(*),T2TPX(*),T2TPY(*),WORK(*),ONE,TWO
#endif
C
      PARAMETER(ONE = 1.0D0, TWO = 2.0D0)
C
      CALL QENTER('XI_DEN_IA')
c
      LENGTH = NCKIJ(ISWMAT)
C
      IF (LENSQ .NE. LENGTH) THEN
         WRITE(LUPRI,*) 'LENSQ = ', LENSQ
         WRITE(LUPRI,*) 'LENGTH = ', LENGTH
         CALL QUIT('Inconsistency in length in XI_DEN_IA')
      END IF
C
      KTMAT  = 1
      KEND1  = KTMAT + NCKIJ(ISWMAT)
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory in XI_DEN_IA (0)')
      ENDIF
C
      ISYMIBID = MULD2H(ISYMIB,ISYMID)
      ISYDEN   = MULD2H(ISWMAT,ISYMIBID)
C
C
C=========================================================================
C    calculate  first contribution: 
C                 DS(ai) = DS(ai) - W^BC(djkl) T2TPX(djiB) T2TPY(alkC)
C                                                         
C
C                                           Gx(kl,i)        Ty^C(kl,a)
C=========================================================================
C
      B = IB
      C = ID
      ISYMB = ISYMIB
      ISYMC = ISYMID
C
C------------------------
C Set symmetry flags
C------------------------
C
      ISYDJI = MULD2H(ISYMT2X,ISYMB)
      ISYALK = MULD2H(ISYMT2Y,ISYMC)
      ISYKLI = MULD2H(ISWMAT,ISYDJI)
C
C-----------------------
C Memory allocation
C-----------------------
C
      KGKLI  = KEND1
      KT2KLA = KGKLI  + NMAIJK(ISYKLI)
      KEND2  = KT2KLA + NCKI(ISYALK)
      LWRK2  = LWORK  - KEND2
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory in XI_DEN_IA (1)')
      ENDIF
C
C---------------------
C Initialize
C---------------------
C
      CALL DZERO(WORK(KTMAT),LENGTH)
      CALL DZERO(WORK(KT2KLA),NCKI(ISYALK))
      CALL DZERO(WORK(KGKLI),NMAIJK(ISYKLI))

C
C-------------------------------------------------------------------------
C     If symmetry, sort W^BC(djkl) as W^BC(dj,kl)
C-------------------------------------------------------------------------
C
      DO I = 1, LENGTH
         WORK(KTMAT-1+I) = WMAT(I)
      ENDDO
C
      IF (NSYM .GT. 1) THEN
         IF (LWRK2 .LT. LENGTH) THEN
            CALL QUIT('Out of memory in XI_DEN_IA'
     *                 //' (CC_GATHER-1)')
         ENDIF
         CALL CC_GATHER(LENGTH,WORK(KEND2),WORK(KTMAT),INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK(KEND2),1,WORK(KTMAT),1)
      ENDIF
C
C-------------------------------------------------------------------------
C     Calculate G(kl,i) = - W^BC(dj,kl) T2TPX^B(dj,i)
C-------------------------------------------------------------------------
C
      DO ISYMDJ = 1,NSYM
C
         ISYMKL = MULD2H(ISWMAT,ISYMDJ)
         ISYMI  = MULD2H(ISYDJI,ISYMDJ)
C
         KOFF1  = KTMAT + ISAIKL(ISYMDJ,ISYMKL)  
         KOFF2  = IT2SP(ISYDJI,ISYMB)
     *          + NCKI(ISYDJI)*(B-1)
     *          + ISAIK(ISYMDJ,ISYMI)
     *          + 1
         KOFF3 = KGKLI  + IMAIJK(ISYMKL,ISYMI)
C
         NTOTDJ = MAX(NT1AM(ISYMDJ),1)
         NTOTKL = MAX(NMATIJ(ISYMKL),1)
C
         CALL DGEMM('T','N',NMATIJ(ISYMKL),NRHF(ISYMI),
     *              NT1AM(ISYMDJ),-ONE,WORK(KOFF1),NTOTDJ,
     *              T2TPX(KOFF2),NTOTDJ,ONE,WORK(KOFF3),NTOTKL)
C
      END DO ! ISYMDJ
C
C-------------------------------------
C     Sort T2   T^C(kl,a) = T2TPY(alkC)
C-------------------------------------
C
      DO ISYMK = 1, NSYM
         ISYMAL = MULD2H(ISYALK,ISYMK)
         DO ISYMA = 1, NSYM
            ISYML = MULD2H(ISYMAL,ISYMA)
            ISYMKL = MULD2H(ISYMK,ISYML)
C
            DO K = 1, NRHF(ISYMK)
               DO L = 1, NRHF(ISYML)
                  KOFF1 = IT2SP(ISYALK,ISYMC)
     *                  + NCKI(ISYALK)*(C-1)
     *                  + ICKI(ISYMAL,ISYMK)
     *                  + NT1AM(ISYMAL)*(K-1)
     *                  + IT1AM(ISYMA,ISYML)
     *                  + NVIR(ISYMA)*(L-1)
     *                  + 1
                  KOFF2 = KT2KLA - 1           
     *                  + IMAIJA(ISYMKL,ISYMA)
     *                  + IMATIJ(ISYMK,ISYML)
     *                  + NRHF(ISYMK)*(L-1)
     *                  + K
C
                  CALL DCOPY(NVIR(ISYMA),T2TPY(KOFF1),1,
     *                       WORK(KOFF2),NMATIJ(ISYMKL))
               ENDDO ! L
            ENDDO    ! K
         ENDDO       ! ISYMA
      ENDDO          ! ISYMK
C
C-------------------------------------------------------------------------
C     Calculate D(ia) = T^C(kl,a)  G(kl,i) 
C-------------------------------------------------------------------------
C
      DO ISYMA = 1, NSYM
         ISYMKL = MULD2H(ISYALK,ISYMA)
         ISYMI  = MULD2H(ISYKLI,ISYMKL)
C
         KOFF1 = KT2KLA + IMAIJA(ISYMKL,ISYMA)
         KOFF2 = KGKLI  + IMAIJK(ISYMKL,ISYMI)
         KOFF3 = IT1AM(ISYMA,ISYMI) + 1
C
         NTOTA  = MAX(NVIR(ISYMA),1)
         NTOTKL = MAX(NMATIJ(ISYMKL),1)
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NMATIJ(ISYMKL),
     *              ONE,WORK(KOFF1),NTOTKL,
     *              WORK(KOFF2),NTOTKL,ONE,DIA(KOFF3),NTOTA)
C
      END DO ! ISYMA

C
C=========================================================================
C    calculate  second contribution: 
C                 DS(ai) = DS(ai) - W^DB(cklj) T2TPX(DjiB) T2TPY(ckla)
C                                                         
C                                          G(ckl,i)        T(ckl,a)
C=========================================================================
C
      D = IB
      B = ID
      ISYMD = ISYMIB
      ISYMB = ISYMID
C
C------------------------
C Set symmetry flags
C------------------------
C
      ISYMBD  = MULD2H(ISYMB,ISYMD)
      ISYMJI  = MULD2H(ISYMBD,ISYMT2X)
      ISYCKLI = MULD2H(ISWMAT,ISYMJI)
C
C-----------------------
C Memory allocation
C-----------------------
C
      KT2JI  = KEND1
      KGCKLI   = KT2JI + NMATIJ(ISYMJI)
      KEND2   = KGCKLI  + NCKIJ(ISYCKLI)
      LWRK2   = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         CALL QUIT('Out of memory in XI_DEN_IA (2)')
      ENDIF
C
C---------------------
C Initialize
C---------------------
C
      CALL DZERO(WORK(KT2JI),NMATIJ(ISYMJI))
      CALL DZERO(WORK(KGCKLI),NCKIJ(ISYCKLI))
C
C ---------------------------------
C     Sort T2TPX(DjiB) as T^{DB}_{ji} 
C ---------------------------------
C
      ISYDJI = MULD2H(ISYMJI,ISYMD)
      DO ISYMI = 1,NSYM
         ISYMDJ = MULD2H(ISYDJI,ISYMI)
         ISYMJ = MULD2H(ISYMJI,ISYMI)
         DO I = 1,NRHF(ISYMI)
            DO J = 1,NRHF(ISYMJ)
               KOFF1 = IT2SP(ISYDJI,ISYMB)
     *                  + NCKI(ISYDJI)*(B-1)
     *                  + ICKI(ISYMDJ,ISYMI)
     *                  + NT1AM(ISYMDJ)*(I-1)
     *                  + IT1AM(ISYMD,ISYMJ)
     *                  + NVIR(ISYMD)*(J-1)
     *                  + D
               KOFF2 = IMATIJ(ISYMJ,ISYMI)
     *                  + NRHF(ISYMJ)*(I-1)
     *                  + J
               WORK(KT2JI-1+KOFF2) = T2TPX(KOFF1)
            ENDDO
         ENDDO
      ENDDO
C
C-------------------------------------------------------------------------
C     Calculate G(ckl,i) = - W^DB(cklj) T^{DB}_{ji}
C-------------------------------------------------------------------------
C
      DO ISYMI = 1,NSYM
C
         ISYMJ   = MULD2H(ISYMJI,ISYMI)
         ISYCKL  = MULD2H(ISWMAT,ISYMJ)
C
         KOFF1  = ISAIKJ(ISYCKL,ISYMJ) + 1
         KOFF2  = IMATIJ(ISYMJ,ISYMI)  + KT2JI
         KOFF3  = ISAIKJ(ISYCKL,ISYMI) + KGCKLI
C
         NTOTJ      =  MAX(1,NRHF(ISYMJ))
         NTOTCKL    =  MAX(1,NCKI(ISYCKL))
C
         CALL DGEMM('N','N',NCKI(ISYCKL),NRHF(ISYMI),
     *               NRHF(ISYMJ),-ONE,WMAT(KOFF1),NTOTCKL,
     *               WORK(KOFF2),NTOTJ,ONE,WORK(KOFF3),NTOTCKL)

      END DO ! ISYMI
C
C-------------------------------------------------------------------------
C       DS(ai) = DS(ai) + Ty(ckl,a)  G(ckl,i) 
C-------------------------------------------------------------------------
C
      DO ISYMA = 1,NSYM
C
         ISYCKL  = MULD2H(ISYMT2Y,ISYMA)
         ISYMI  = MULD2H(ISYCKLI,ISYCKL)
C
         KOFF1  = IT2SP(ISYCKL,ISYMA) + 1
         KOFF2  = ISAIKJ(ISYCKL,ISYMI) + KGCKLI
         KOFF3  = IT1AM(ISYMA,ISYMI)  + 1
C
         NTOTCKL   = MAX(1,NCKI(ISYCKL))
         NTOTA  = MAX(1,NVIR(ISYMA))
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),
     *               NCKI(ISYCKL),ONE,T2TPY(KOFF1),NTOTCKL,
     *               WORK(KOFF2),NTOTCKL,ONE,DIA(KOFF3),NTOTA)
      END DO ! ISYMA

C
C========================================================================= C    calculate  third contribution: 
C                 DS(ai) = DS(ai) - W^CD(bljk) T2TPX(bijD) T2TPY(alkC)
C
C                                          G(lk,i)        T^C(lk,a)
C=========================================================================
C
      C = IB
      D = ID
      ISYMC = ISYMIB
      ISYMD = ISYMID
C
C------------------------
C Set symmetry flags
C------------------------
C
      ISYBIJ = MULD2H(ISYMT2X,ISYMD)
      ISYALK = MULD2H(ISYMT2Y,ISYMC)
      ISYLKI = MULD2H(ISWMAT,ISYBIJ)
C
C-----------------------
C Memory allocation
C-----------------------
C
      KGLKI  = KEND1
      KEND2  = KGLKI  + NMAIJK(ISYLKI)
C
      KT2BJI = KEND2
      KEND2  = KT2BJI + NCKI(ISYBIJ)
      LWRK2  = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         CALL QUIT('Out of memory in XI_DEN_IA (3)')
      ENDIF
C
C---------------------
C Initialize
C---------------------
C
      CALL DZERO(WORK(KTMAT),LENGTH)
      CALL DZERO(WORK(KT2BJI),NCKI(ISYBIJ))
      CALL DZERO(WORK(KGLKI),NMAIJK(ISYLKI))
C
C-------------------------------------------------------------------------
C     Sort W^CD(bljk) as W^CD(bjlk)
C-------------------------------------------------------------------------
C
      DO I = 1, LENGTH
         WORK(KTMAT-1+I) = WMAT(INDSQ(I,1))
      ENDDO
C
C-------------------------------------------------------------------------
C     If symmetry, sort W^CD(bjlk) as W^CD(bj,lk)
C-------------------------------------------------------------------------
C
      IF (NSYM .GT. 1) THEN
         IF (LWRK2 .LT. LENGTH) THEN
            CALL QUIT('Out of memory in XI_DEN_IA'
     *                 //' (CC_GATHER-2)')
         ENDIF
         CALL CC_GATHER(LENGTH,WORK(KEND2),WORK(KTMAT),INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK(KEND2),1,WORK(KTMAT),1)
      ENDIF
C
C-------------------------------------
C     Sort T2   T^XD(bj,i) = T2TPX(bijD)
C-------------------------------------
C
      DO ISYMJ = 1, NSYM
         ISYMBI = MULD2H(ISYBIJ,ISYMJ)
         DO ISYMB = 1, NSYM
            ISYMI = MULD2H(ISYMBI,ISYMB)
            ISYMBJ = MULD2H(ISYMB,ISYMJ)
C
            DO J = 1, NRHF(ISYMJ)
               DO I = 1, NRHF(ISYMI)
                  KOFF1 = IT2SP(ISYBIJ,ISYMD)
     *                  + NCKI(ISYBIJ)*(D-1)
     *                  + ICKI(ISYMBI,ISYMJ)
     *                  + NT1AM(ISYMBI)*(J-1)
     *                  + IT1AM(ISYMB,ISYMI)
     *                  + NVIR(ISYMB)*(I-1)
     *                  + 1
                  KOFF2 = KT2BJI        
     *                  + ICKI(ISYMBJ,ISYMI)
     *                  + NT1AM(ISYMBJ)*(I-1)
     *                  + IT1AM(ISYMB,ISYMJ)
     *                  + NVIR(ISYMB)*(J-1)
C
                  CALL DCOPY(NVIR(ISYMB),T2TPX(KOFF1),1,
     *                       WORK(KOFF2),1) 
               ENDDO ! I
            ENDDO    ! J
         ENDDO       ! ISYMB
      ENDDO          ! ISYMJ
C    
C-------------------------------------------------------------------------
C     Calculate G(lk,i) = - W^CD(bj,lk) T^XD(bj,i)     
C-------------------------------------------------------------------------
C
      DO ISYMBJ = 1,NSYM
C
         ISYMLK = MULD2H(ISWMAT,ISYMBJ)
         ISYMI  = MULD2H(ISYBIJ,ISYMBJ)
C
         KOFF1  = KTMAT  + ISAIKL(ISYMBJ,ISYMLK) 
         KOFF2  = KT2BJI + ICKI(ISYMBJ,ISYMI)
         KOFF3  = KGLKI  + IMAIJK(ISYMLK,ISYMI)
C
         NTOTBJ = MAX(NT1AM(ISYMBJ),1)
         NTOTLK = MAX(NMATIJ(ISYMLK),1)
C
         CALL DGEMM('T','N',NMATIJ(ISYMLK),NRHF(ISYMI),
     *              NT1AM(ISYMBJ),-ONE,WORK(KOFF1),NTOTBJ,
     *              WORK(KOFF2),NTOTBJ,ONE,WORK(KOFF3),NTOTLK)
C

      END DO ! ISYMBJ
C
C-------------------------------------
C     Sort T2   T^yC(lk,a) = T2TPY(alkC)
C-------------------------------------
C
      KT2ALK = KEND2
      KEND2  = KT2ALK + NCKI(ISYALK)
      LWRK2  = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         CALL QUIT('Out of memory in XI_DEN_IA (4)')
      ENDIF
C
      CALL DZERO(WORK(KT2ALK),NCKI(ISYALK))
C
      DO ISYML = 1, NSYM
         ISYMAK = MULD2H(ISYALK,ISYML)
         DO ISYMA = 1, NSYM
            ISYMK = MULD2H(ISYMAK,ISYMA)
            ISYMLK = MULD2H(ISYML,ISYMK)
            ISYMAL = MULD2H(ISYML,ISYMA)
C
            DO L = 1, NRHF(ISYML)
               DO K = 1, NRHF(ISYMK)
                  KOFF1 = IT2SP(ISYALK,ISYMC)
     *                  + NCKI(ISYALK)*(C-1)
     *                  + ICKI(ISYMAL,ISYMK)
     *                  + NT1AM(ISYMAL)*(K-1)
     *                  + IT1AM(ISYMA,ISYML)
     *                  + NVIR(ISYMA)*(L-1)
     *                  + 1
                  KOFF2 = KT2ALK - 1
     *                  + IMAIJA(ISYMLK,ISYMA)
     *                  + IMATIJ(ISYML,ISYMK)
     *                  + NRHF(ISYML)*(K-1)
     *                  + L
C
                  CALL DCOPY(NVIR(ISYMA),T2TPY(KOFF1),1,
     *                       WORK(KOFF2),NMATIJ(ISYMLK))
               ENDDO ! K
            ENDDO    ! L
         ENDDO       ! ISYMA
      ENDDO          ! ISYML
C
C-------------------------------------------------------------------------
C     Calculate DS(ai) = T^YC(lk,a)  G(lk,i) 
C-------------------------------------------------------------------------
C
      DO ISYMA = 1, NSYM
         ISYMLK = MULD2H(ISYALK,ISYMA)
         ISYMI  = MULD2H(ISYLKI,ISYMLK)
C
         KOFF1 = KT2ALK + IMAIJA(ISYMLK,ISYMA)
         KOFF2 = KGLKI  + IMAIJK(ISYMLK,ISYMI)
         KOFF3 = IT1AM(ISYMA,ISYMI) + 1
C
         NTOTA  = MAX(NVIR(ISYMA),1)
         NTOTLK = MAX(NMATIJ(ISYMLK),1)
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NMATIJ(ISYMLK),
     *              ONE,WORK(KOFF1),NTOTLK,
     *              WORK(KOFF2),NTOTLK,ONE,DIA(KOFF3),NTOTA)
C
      END DO ! ISYMA

      CALL QEXIT('XI_DEN_IA')
C
      RETURN
      END
C  /* Deck cc3_xi_den_abij */
      SUBROUTINE CC_XI_DEN_ABIJ(QUADR,LISTR,DAB,DIJ,DO_DIA,DAI,
     *                           L2L1,ISYML1,
     *                           ISYFCK,FOCK,
     *                           ISYMT3,ISWMAT,ISTHETA,
     *                           LUT3,FNT3,LUWBMAT,FNWBMAT,
     *                           LUTHETA,FNTHETA,
     *                           FOCKD,FREQ,
     *                           WORK,LWORK,ISYMD,D)
C
C=========================================================================
C
C    QUADR has to be .FALSE. for linear response calculations
C
C    QUADR has to be .TRUE.  for quadratic response calculations
C=========================================================================
C
C    if  quadr   calculate 
C
C            DAI(ai) = DAI(ai) +  L2L1{emLD}*(THETA{Dea}_{Lmi} - THETA{Dea}_{Lim}) 
C
C    Calculate the densities 
C
C     (1)    D(ab) = <L3Y|[E_ab,T3]|HF> 
C
C            D(ab) = 1/2 tbarY^(dea)_(lmn) t^(deb)_(lmn) 
C
C                  = ( 1/2 W^de(anml) + W^da(emnl) ) t^(deb)_(lmn)
C                  =       Wtilde_bar^DL(em,aN)      T3^DL(em,bN)
C 
C   if quadr              t^(deb)_(lmn) = theta^(deb)_(l-m-n-)
C                                       = THETA^DL-(em-,bN-)
C
C                         THETA^DL-(em-,bN-) is symmetrized inside
C                         the routine
C
C------------------------------------------------------------------------
C
C
C     (2)    D(ij) = <L3Y|[E_ij,T3]|HF> 
C
C            D(ij) = - 1/2 tbarY^(def)_(lmj) t^(def)_(lmi) 
C
C                  = ( 1/2 W^de(fjml) + W^df(emjl) ) t^(def)_(lmi)
C                  =       Wtilde_bar^DL(em,aJ)      T3^DL(em,fi)
C
C   if quadr    t^(def)_(lmi) = [theta^(def)_(l-m-i-) + theta^(def-)_(lmi)]
C                             = THETA^DL(em,fi)
C
C------------------------------------------------------------------------
C
C   if quadr   calculate additionally :
C
C    D(ij) = D(ij) + W^df(emlj) [ theta^(def-)_(iml) + theta^(de-f)_(iml) ]
C
C
C
C
C     Common notation for (1) and (2):
C
C                   a -> f     j -> n
C                   ======     ======
C
C     Outside the densities D(ab) snd D(ij) are to be contracted 
C     with X_(ab) and X_(ij), respectively
C
C
C    Written by Filip Pawlowski, Fall 2002, Aarhus
C=========================================================================
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <cc3t3d.h>
C
      LOGICAL DO_DIA

      CHARACTER LISTR*3

      INTEGER ISYMT3,ISWMAT,LUT3,LUWBMAT,LWORK,ISYMD
      INTEGER ISYML,ISYMDL,ISWMATDL,ISYMT3DL,ISYMN,ISYEMF,ISYMBN,ISYMEM
      INTEGER ISYMB,ISYMF,ISYMFI,ISYMI,ISYEMB,ISYMFN
      INTEGER KT3,KWMAT,KEND1,LWRK1
      INTEGER KOFF1,KOFF2,KOFF3,KBN,KFN
      INTEGER NTOTEM,NTOTF,NNEMF
      INTEGER IADR
C
      INTEGER ISTHETA,ISYFCK,LUTHETA
      INTEGER ISTHETADL,ISTHETAF
      INTEGER KTHETA,KTHETAF
      INTEGER KFI
      INTEGER IOPT
C
      INTEGER ISYMJ,ISYMFJ,KFJ
      INTEGER ISYMM,ISYME
C
      INTEGER ISYMDAI,ISYML1
C
      LOGICAL QUADR
C
      CHARACTER*(*) FNT3,FNWBMAT,FNTHETA
C
#if defined (SYS_CRAY)
      REAL DAB(*),DIJ(*),DAI(*),WORK(LWORK),ONE,HALF
      REAL FOCK(*),L2L1(*),FOCKD(*),FREQ
      REAL XNORMVAL,DDOT
#else
      DOUBLE PRECISION DAB(*),DIJ(*),DAI(*),WORK(LWORK),ONE,HALF
      DOUBLE PRECISION FOCK(*),L2L1(*),FOCKD(*),FREQ
      DOUBLE PRECISION XNORMVAL,DDOT
#endif
C
      PARAMETER(ONE = 1.0D0, HALF = 0.5D0)
C
      CALL QENTER('XIDENABIJ')
C
      DO ISYML = 1,NSYM
C
         ISYMDL = MULD2H(ISYMD,ISYML)
         ISWMATDL = MULD2H(ISWMAT,ISYMDL)
         ISYMT3DL = MULD2H(ISYMT3,ISYMDL)
         IF (QUADR) THEN
            ISTHETADL = MULD2H(ISTHETA,ISYMDL)
            ISTHETAF  = MULD2H(ISYMT3DL,ISYFCK)
            ! symmetry check
            IF (ISTHETADL .NE. ISTHETAF) THEN
               WRITE(LUPRI,*)'ISTHETADL : ', ISTHETADL
               WRITE(LUPRI,*)'ISTHETAF  : ', ISTHETAF
               CALL QUIT('Symmetry mismatch in CC3_XI_DEN_ABIJ')
            END IF
         END IF
C
         KT3  = 1
         KWMAT  = KT3 + NT2SQ(ISYMT3DL)
         KEND1 = KWMAT + NT2SQ(ISWMATDL)
         LWRK1  = LWORK - KEND1
         IF (QUADR) THEN
            KTHETA  = KEND1 
            KTHETAF = KTHETA  + NT2SQ(ISTHETADL) 
            KEND1   = KTHETAF + NT2SQ(ISTHETAF)
            LWRK1   = LWORK   - KEND1
         END IF
C
         IF ( LWRK1 .LT. 0 ) THEN
           CALL QUIT('Out of memory in CC3_XI_DEN_ABIJ (x)')
         ENDIF
C
         DO L = 1, NRHF(ISYML)
C
C           --------------------------------------------
C           Read T3 amplitudes from file:
C           --------------------------------------------
C
            IADR = ISWTL(ISYMT3DL,ISYML) + NT2SQ(ISYMT3DL)*(L-1) + 1
            CALL GETWA2(LUT3,FNT3,WORK(KT3),IADR,NT2SQ(ISYMT3DL))
c
C
            IF (QUADR) THEN
C           ----------------------------------------------------
C           Read THETA amplitudes from file and symmetrize them:
C           ----------------------------------------------------
C
               IADR = ISWTL(ISTHETADL,ISYML) + NT2SQ(ISTHETADL)*(L-1) 
     *              + 1
               CALL GETWA2(LUTHETA,FNTHETA,WORK(KTHETA),IADR,
     *                     NT2SQ(ISTHETADL))
C
               CALL CC_T2MOD(WORK(KTHETA),ISTHETADL,ONE)
C
C-----------------------------------------------------------------------
C  DAI(ai) = DAI(ai) +  L2L1{emLD}*(THETA{Dea}_{Lmi} - THETA{Dea}_{Lim})
C-----------------------------------------------------------------------
C
               ISYMDAI = MULD2H(ISTHETA,ISYML1)
C
               IF (DO_DIA) THEN
                  CALL ADEN_DAI_T2_D(DAI,ISYMDAI,L2L1,ISYML1,
     *                               WORK(KTHETA),ISTHETADL,ISYMD,D,
     *                               ISYML,L,WORK(KEND1),LWRK1)
               END IF
            END IF
C
C           ------------------------------------------------
C           Read WMAT_bar from file and generate WMAT-tilde:
C           ------------------------------------------------
C
            IADR = ISWTL(ISWMATDL,ISYML) + NT2SQ(ISWMATDL)*(L-1) + 1
            CALL GETWA2(LUWBMAT,FNWBMAT,WORK(KWMAT),IADR,
     *                  NT2SQ(ISWMATDL))
C
            CALL CC_T2MOD(WORK(KWMAT),ISWMATDL,HALF)

C
C           -----------------------------------
C           Loop over outermost occupied index:
C           -----------------------------------
            DO ISYMN = 1, NSYM
               ISYEMF = MULD2H(ISWMATDL,ISYMN)
               IF (QUADR) THEN
                  ISYEMB = MULD2H(ISTHETADL,ISYMN)
               ELSE
                  ISYEMB = MULD2H(ISYMT3DL,ISYMN)
               END IF
C
               DO N = 1, NRHF(ISYMN)
C
C                 -------------------------------------------------------
C                 D(fb) <- D(fb)+ sum_em Wtilde_bar^DL(em,fN) T3^DL(em,bN):
                  ! FOR QUADRU = .TRUE. T3^DL(em,bN) becomes THETA_Z^DL(em,bN)
C                 -------------------------------------------------------
                  DO ISYMEM = 1, NSYM
                     ISYMB  = MULD2H(ISYEMB,ISYMEM)
                     ISYMF  = MULD2H(ISYEMF,ISYMEM)
                     ISYMFN = MULD2H(ISYMF,ISYMN)
                     ISYMBN = MULD2H(ISYMB,ISYMN)

                     KFN    = IT1AM(ISYMF,ISYMN)+NVIR(ISYMF)*(N-1)+1
                     KOFF1  = KWMAT + IT2SQ(ISYMEM,ISYMFN)
     *                              + NT1AM(ISYMEM)*(KFN-1)
                     KBN    = IT1AM(ISYMB,ISYMN)+NVIR(ISYMB)*(N-1)+1
C                     
                     IF (QUADR) THEN
                        KOFF2  = KTHETA + IT2SQ(ISYMEM,ISYMBN)
     *                                  + NT1AM(ISYMEM)*(KBN-1)
                     ELSE
                        KOFF2  = KT3    + IT2SQ(ISYMEM,ISYMBN)
     *                                  + NT1AM(ISYMEM)*(KBN-1)
                     END IF
C
                     KOFF3  = IMATAB(ISYMF,ISYMB) + 1


                     NTOTEM = MAX(NT1AM(ISYMEM),1)
                     NTOTF  = MAX(NVIR(ISYMF),1)

                     CALL DGEMM('T','N',NVIR(ISYMF),NVIR(ISYMB),
     *                          NT1AM(ISYMEM),ONE,WORK(KOFF1),NTOTEM,
     *                          WORK(KOFF2),NTOTEM,ONE,DAB(KOFF3),NTOTF)

                  END DO ! ISYMEM
C
               END DO ! N
            END DO    ! ISYMN
C
            IF (QUADR .AND. (LISTR(1:3).EQ.'R1 ') ) THEN
               ! Construct THETA^DL(em,f- i)
               CALL DZERO(WORK(KTHETAF),NT2SQ(ISTHETAF))


               IOPT = 1
               CALL WXDL_V(IOPT,WORK(KT3),ISYMT3DL,FOCK,ISYFCK,
     *                 WORK(KTHETAF),ISTHETAF,WORK(KEND1),LWRK1)
C
               CALL W3DL_DIA(WORK(KTHETAF),ISTHETAF,ISYML,L,ISYMD,D,
     *                       FOCKD,FREQ)
               CALL T3_FORBIDDEN_DL(WORK(KTHETAF),ISTHETA,ISYMD,D,
     *                              ISYML,L)

C
               ! add THETA^DL(em,f- i) to THETA^DL-(em-,fi-) sitting in 
               ! KTHETA
               DO ISYMI = 1,NSYM
                  DO ISYMF = 1,NSYM
                     ISYMFI = MULD2H(ISYMF,ISYMI)
                     ISYMEM = MULD2H(ISTHETADL,ISYMFI)
                     DO ISYMM = 1,NSYM
                        ISYME = MULD2H(ISYMEM,ISYMM)
                        DO I = 1,NRHF(ISYMI)
                           DO F = 1,NVIR(ISYMF)
                              DO M = 1,NRHF(ISYMM)
                                 DO E = 1,NVIR(ISYME)
                                    KFI   = IT1AM(ISYMF,ISYMI)
     *                                    + NVIR(ISYMF)*(I-1)
     *                                    + F
                                    KOFF1 = IT2SQ(ISYMEM,ISYMFI)
     *                                    + NT1AM(ISYMEM)*(KFI-1)
     *                                    + IT1AM(ISYME,ISYMM)
     *                                    + NVIR(ISYME)*(M-1)
     *                                    + E -1
C
                                  WORK(KTHETA+KOFF1)=WORK(KTHETA+KOFF1)
     *                                              +WORK(KTHETAF+KOFF1)
                                 END DO
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
C



            END IF
C
            DO ISYMN = 1, NSYM
               ISYEMF = MULD2H(ISWMATDL,ISYMN)
C
               DO N = 1, NRHF(ISYMN)

C                   -------------------------------------------------------
C                   D(iN) <- D(iN)- sum_emf Wtilde_bar^DL(em,fN) t^DL(em,fi):
                  ! FOR QUADRU = .TRUE. t^DL(em,fi) becomes THETA_Z^DL(em,fi)
C                   -------------------------------------------------------
                  DO ISYMEM = 1, NSYM
                     IF (QUADR) THEN
                        ISYMFI = MULD2H(ISTHETADL,ISYMEM)
                     ELSE
                        ISYMFI = MULD2H(ISYMT3DL,ISYMEM)
                     END IF
                     ISYMF  = MULD2H(ISYEMF,ISYMEM)
                     ISYMI  = MULD2H(ISYMFI,ISYMF)
                     ISYMFN = MULD2H(ISYMF,ISYMN)

C
                     IF (QUADR) THEN
                        KOFF1 = KTHETA+ IT2SQ(ISYMEM,ISYMFI)
     *                                + NT1AM(ISYMEM)*IT1AM(ISYMF,ISYMI)
                     ELSE
                        KOFF1  = KT3  + IT2SQ(ISYMEM,ISYMFI) 
     *                                + NT1AM(ISYMEM)*IT1AM(ISYMF,ISYMI)
                     END IF
C
                     KFN    = IT1AM(ISYMF,ISYMN)+NVIR(ISYMF)*(N-1)+1
                     KOFF2  = KWMAT + IT2SQ(ISYMEM,ISYMFN)
     *                              + NT1AM(ISYMEM)*(KFN-1)

                     KOFF3  = IMATIJ(ISYMI,ISYMN) 
     *                              + NRHF(ISYMI)*(N-1) + 1

                     NNEMF  = MAX(NT1AM(ISYMEM)*NVIR(ISYMF),1)
C



                     CALL DGEMV('T',NT1AM(ISYMEM)*NVIR(ISYMF),
     *                         NRHF(ISYMI),-ONE,WORK(KOFF1),NNEMF,
     *                         WORK(KOFF2),1,ONE,DIJ(KOFF3),1)
C
                  END DO ! ISYMFI

               END DO ! N 
            END DO    ! ISYMN 
C
            IF (QUADR .AND. (LISTR(1:3).EQ.'R1 ') ) THEN
       ! Calculate the extra contribution to D(ij) density:
       ! D(ij) = D(ij) + W^Df(emlj) * [theta^(Def-)_(iml) + theta^(De-f)_(iml)]
C
C  ----------------------------
C  Read T3^DL(em,fi) amplitudes
C  ----------------------------
C
               ! KT3 is recycled here
               CALL READ_T3_AIBL(LUT3,FNT3,ISYMT3,WORK(KT3),
     *                           ISYMT3DL,L,ISYML,ISYMD)
C
C ----------------------------------------------
C Contract T3^DL(em,fi) with X operator 
C to get THDL(em,fi) = [ THETA^DL(em,f-i) + THETA^DL(e-m,fi) ]
C ----------------------------------------------
C
               ! KTHETAF is recycled here
               CALL DZERO(WORK(KTHETAF),NT2SQ(ISTHETAF))
               IOPT = 2
               CALL WXDL_V(IOPT,WORK(KT3),ISYMT3DL,FOCK,ISYFCK,
     *                     WORK(KTHETAF),ISTHETAF,WORK(KEND1),LWRK1)
               ! Divide it by orbital energy difference and remove the 
               ! forbidden elements
               CALL W3DL_DIA(WORK(KTHETAF),ISTHETAF,ISYML,L,ISYMD,D,
     *                       FOCKD,FREQ)
               CALL T3_FORBIDDEN_DL(WORK(KTHETAF),ISTHETA,ISYMD,D,
     *                              ISYML,L)
C
C ------------------------------------------------
C Read WBMAT^DL(em,fj) from the file
C ------------------------------------------------
C
               ! KWMAT is recycled here
               CALL READ_T3_AIBL(LUWBMAT,FNWBMAT,ISWMAT,WORK(KWMAT),
     *                           ISWMATDL,L,ISYML,ISYMD)

C
C------------------------------------------------
C Contract D(ij) <- WBMAT^DL(em,fj) * THDL(em,fi) 
C------------------------------------------------
C
               DO ISYMJ = 1,NSYM
                  ISYEMF = MULD2H(ISWMATDL,ISYMJ)
                  DO J = 1,NRHF(ISYMJ)
                     DO ISYMEM = 1, NSYM
                        ISYMFI = MULD2H(ISTHETAF,ISYMEM)
                        ISYMF  = MULD2H(ISYEMF,ISYMEM)
                        ISYMI  = MULD2H(ISYMFI,ISYMF)
                        ISYMFJ = MULD2H(ISYMF,ISYMJ)

C
                        KOFF1 = KTHETAF+ IT2SQ(ISYMEM,ISYMFI)
     *                                + NT1AM(ISYMEM)*IT1AM(ISYMF,ISYMI)

                        KFJ    = IT1AM(ISYMF,ISYMJ)+NVIR(ISYMF)*(J-1)+1
                        KOFF2  = KWMAT + IT2SQ(ISYMEM,ISYMFJ)
     *                              + NT1AM(ISYMEM)*(KFJ-1)

                        KOFF3  = IMATIJ(ISYMI,ISYMJ)
     *                              + NRHF(ISYMI)*(J-1) + 1

                        NNEMF  = MAX(NT1AM(ISYMEM)*NVIR(ISYMF),1)
C
                        CALL DGEMV('T',NT1AM(ISYMEM)*NVIR(ISYMF),
     *                            NRHF(ISYMI),-ONE,WORK(KOFF1),NNEMF,
     *                            WORK(KOFF2),1,ONE,DIJ(KOFF3),1)
C

                     END DO ! ISYMEM

               END DO ! J 
            END DO    ! ISYMJ 
C
            END IF 
C
         END DO       ! L 
      END DO          ! ISYML 
C
      CALL QEXIT('XIDENABIJ')
C
      RETURN
      END
C  /* Deck wxdl_v */
      SUBROUTINE WXDL_V(IOPT,TMAT,ISTMAT,FOCKY,ISYFKY, 
     *                 WMAT,ISWMAT,WRK,LWRK)
*
*********************************************************
*
* if (iopt .eq. 1) then calculate ONLY : 
*    WDL(em,f-i) = - TDL(em,di) * X(f,d)
* else if (iopt .eq. 2) then calculate ADDITIONALLY :
*    WDL(e-m,fi) = - TDL(dm,fi) * X(e,d)
* else if (iopt .eq. 3) then calculate ONLY :
*    WDL(e-m,fi) = - TDL(dm,fi) * X(e,d)
* end if
*
*     F. Pawlowski, Winter 2003, Aarhus.
*     (modified IOPT = 3, 14-Sep-2003, Aarhus.)
********************************************************
*
      IMPLICIT NONE
C
#include <priunit.h>
#include <dummy.h>
#include <iratdef.h>
#include <ccsdsym.h>
#include <inftap.h>
#include <ccinftap.h>
#include <ccorb.h>
#include <ccsdinp.h>
C
      INTEGER IOPT,ISTMAT,ISYFKY,ISWMAT,LWRK
      INTEGER KFCFD,KEND0,LWRK0
      INTEGER ISYMD,ISYMF
      INTEGER KOFF1,KOFF2
      INTEGER ISYMI,ISYMDI,ISYMEM,ISYMFI
      INTEGER NDI,NFI,KFI
      INTEGER KOFF3
      INTEGER NTOTEM,NTOTF
C
      INTEGER ISYMDMF,ISYME,ISYMDM,ISYMM,NTOTE,NTOTD
C
#if defined (SYS_CRAY)
      REAL TMAT(*),FOCKY(*),WMAT(*),WRK(LWRK)
      REAL ONE
      real xnormval,ddot
#else
      DOUBLE PRECISION TMAT(*),FOCKY(*),WMAT(*),WRK(LWRK)
      DOUBLE PRECISION ONE
      double precision xnormval,ddot
#endif
C
      PARAMETER (ONE = 1.0D0)
C
      CALL QENTER('WXDL_V')
C
C Initial test of iopt
C
      IF ( (IOPT .NE. 1) .AND. (IOPT .NE. 2) .AND. (IOPT .NE. 3)) THEN
         WRITE(LUPRI,*)'IOPT : ', IOPT
         CALL QUIT('Wrong option in WXDL_V')
      END IF
C
C
C RESORT VIR-VIR FOCKY ELEMENTS (F,D)
C
C
      KFCFD  = 1
      KEND0  = KFCFD + NMATAB(ISYFKY)
      LWRK0  = LWRK  - KEND0
C
      IF (LWRK0 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWRK0
         WRITE(LUPRI,*) 'Memory needed    : ',KEND0
         CALL QUIT('Insufficient space in WXDL_V')
      END IF
C
      DO ISYMD = 1,NSYM
         ISYMF = MULD2H(ISYMD,ISYFKY)
         DO D = 1,NVIR(ISYMD)
            KOFF1 = IFCVIR(ISYMF,ISYMD) + NORB(ISYMF)*(D - 1)
     *                                  + NRHF(ISYMF) + 1
            KOFF2 = KFCFD + IMATAB(ISYMF,ISYMD) + NVIR(ISYMF)*(D - 1) 
            CALL DCOPY(NVIR(ISYMF),FOCKY(KOFF1),1,WRK(KOFF2),1)
         END DO
      END DO


      IF (IOPT .LE. 2) THEN
C
C
C CARRY OUT MATRIX MULTIPLICATION
C WDL(em,fi) = TDL(em,di) * X(f,d)
C              I(em,d)      X(f,d)
C 
         DO ISYMI = 1,NSYM
            DO ISYMD = 1,NSYM
               ISYMDI = MULD2H(ISYMD,ISYMI)
               ISYMEM = MULD2H(ISYMDI,ISTMAT)
               ISYMF  = MULD2H(ISYFKY,ISYMD)
               ISYMFI = MULD2H(ISYMF,ISYMI)
               DO I = 1,NRHF(ISYMI)
C
                  NDI   = IT1AM(ISYMD,ISYMI)+NVIR(ISYMD)*(I-1)+1
                  NFI   = IT1AM(ISYMF,ISYMI)+NVIR(ISYMF)*(I-1)+1
C
                  KOFF1 = IT2SQ(ISYMEM,ISYMDI)
     *                  + NT1AM(ISYMEM)*(NDI-1)
     *                  + 1
                  KOFF2 = KFCFD 
     *                  + IMATAB(ISYMF,ISYMD) 
                  KOFF3 = IT2SQ(ISYMEM,ISYMFI)
     *                  + NT1AM(ISYMEM)*(NFI-1)
     *                  + 1
C
                  NTOTEM = MAX(NT1AM(ISYMEM),1)
                  NTOTF  = MAX(NVIR(ISYMF),1)
C                       
                  CALL DGEMM('N','T',NT1AM(ISYMEM),NVIR(ISYMF), 
     *                       NVIR(ISYMD),-ONE,TMAT(KOFF1),NTOTEM,
     *                       WRK(KOFF2),NTOTF,ONE,WMAT(KOFF3),NTOTEM)
c
               END DO
            END DO
         END DO
C
      END IF
C
      IF (IOPT .GE. 2) THEN
C
C CARRY OUT MATRIX MULTIPLICATION
C WDL(em,fi) = TDL(dm,fi) * X(e,d)
C
C              X(e,d)   I(d,m)      
C 
         DO ISYMI = 1,NSYM
            ISYMDMF = MULD2H(ISTMAT,ISYMI)
            DO ISYMD = 1,NSYM
               ISYME  = MULD2H(ISYFKY,ISYMD)
               DO ISYMF  = 1,NSYM
                  ISYMDM = MULD2H(ISYMDMF,ISYMF)
                  ISYMFI = MULD2H(ISYMI,ISYMF)
                  ISYMM  = MULD2H(ISYMDM,ISYMD)
                  ISYMEM = MULD2H(ISYME,ISYMM)
                  DO I = 1,NRHF(ISYMI)
                     DO F = 1,NVIR(ISYMF)
C
                        KFI   = IT1AM(ISYMF,ISYMI)+NVIR(ISYMF)*(I-1)+F
C
                        KOFF1 = KFCFD 
     *                        + IMATAB(ISYME,ISYMD) 
                        KOFF2 = IT2SQ(ISYMDM,ISYMFI)
     *                        + NT1AM(ISYMDM)*(KFI-1)
     *                        + IT1AM(ISYMD,ISYMM)
     *                        + 1
                        KOFF3 = IT2SQ(ISYMEM,ISYMFI)
     *                        + NT1AM(ISYMEM)*(KFI-1)
     *                        + IT1AM(ISYME,ISYMM)
     *                        + 1
C
                        NTOTE = MAX(NVIR(ISYME),1)
                        NTOTD  = MAX(NVIR(ISYMD),1)
C                    
                        CALL DGEMM('N','N',NVIR(ISYME),NRHF(ISYMM), 
     *                             NVIR(ISYMD),-ONE,WRK(KOFF1),NTOTE,
     *                             TMAT(KOFF2),NTOTD,ONE,WMAT(KOFF3),
     *                             NTOTE)
C
                     END DO
                  END DO
               END DO
            END DO
         END DO
C
      END IF
C
      CALL QEXIT('WXDL_V')
C
      RETURN
      END
C  /* Deck w3dl_dia */
      SUBROUTINE W3DL_DIA(TEMFI,ISYMEMFI,ISYML,L,ISYMD,D,FOCKD,FREQ)
C
C  W3^DL(em,fi) = W3^DL(em,fi) / (F(D)+F(E)+F(F)-F(L)-F(M)-F(I) - FREQ)
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
C
      INTEGER ISYMEMFI,ISYML,ISYMD
      INTEGER ISYMI,ISYMEMF,ISYMF,ISYMEM,ISYMFI,ISYMM,ISYME
      INTEGER NL,ND,NI,NF,NM,NE
      INTEGER KFI,KOFF1
C
#if defined (SYS_CRAY)
      REAL TEMFI(*),FOCKD(*)
      REAL DEN,FREQ
#else
      DOUBLE PRECISION TEMFI(*),FOCKD(*)
      DOUBLE PRECISION DEN,FREQ  
#endif
C
      CALL QENTER('W3DL_DIA')
C
      NL = IORB(ISYML) + L
      ND = IORB(ISYMD) + NRHF(ISYMD) + D
      DO ISYMI = 1,NSYM
         ISYMEMF = MULD2H(ISYMEMFI,ISYMI)
         DO ISYMF = 1,NSYM
            ISYMEM = MULD2H(ISYMEMF,ISYMF)
            ISYMFI = MULD2H(ISYMF,ISYMI)
            DO ISYMM = 1,NSYM
               ISYME = MULD2H(ISYMEM,ISYMM)
               DO I = 1,NRHF(ISYMI)
                  NI = IORB(ISYMI) + I
                  DO F = 1,NVIR(ISYMF)
                     NF = IORB(ISYMF) + NRHF(ISYMF) + F
C
                     KFI   = IT1AM(ISYMF,ISYMI)+NVIR(ISYMF)*(I-1)+F
C
                     DO M = 1,NRHF(ISYMM)
                        NM = IORB(ISYMM) + M
                        DO E = 1,NVIR(ISYME)
                           NE = IORB(ISYME) + NRHF(ISYME) + E
C
                           KOFF1 = IT2SQ(ISYMEM,ISYMFI)
     *                           + NT1AM(ISYMEM)*(KFI-1)
     *                           + IT1AM(ISYME,ISYMM)
     *                           + NVIR(ISYME)*(M-1)
     *                           + E
C
                           DEN   = FOCKD(ND) + FOCKD(NE) + FOCKD(NF)
     *                           - FOCKD(NL) - FOCKD(NM) - FOCKD(NI)
     *                           - FREQ
C
                           TEMFI(KOFF1) = TEMFI(KOFF1)/DEN

C
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
C
C     End.
C
      CALL QEXIT('W3DL_DIA')
C
      RETURN
      END
C
C  /* Deck t3_forbidden_dl */
      SUBROUTINE T3_FORBIDDEN_DL(TMAT,ISYMIM,ISYMD,D,ISYML,L)
C
C     Purpose : Remove the forbidden t3/t3-bar amplitudes
C     sitting as T3^DL(em,fi)
C
      IMPLICIT NONE
C
#include <ccsdsym.h>
#include <ccorb.h>
#include <priunit.h>
C
      INTEGER ISYMIM,ISYMD,ISYML
      INTEGER ISYMM,ISYDEF,ISYMF,ISYDE,ISYME,ISYMEM,ISYMFI
      INTEGER NFI,KOFF1
      INTEGER ISYLMI,ISYMI,ISYLM
C
#if defined (SYS_CRAY)
      REAL TMAT(*),ZERO
#else
      DOUBLE PRECISION TMAT(*),ZERO
#endif
C
      PARAMETER (ZERO = 0.0D0)
C
      CALL QENTER('T3_FORBIDDEN_DL')
C
C---------------------------------------------------------
C     If L and M are the same remove all amplitudes
C     having an I which is the same as L and M.
C---------------------------------------------------------
C
      DO ISYMM = 1,NSYM
         DO M = 1,NRHF(ISYMM)
            IF ((ISYML .EQ. ISYMM) .AND. (L .EQ. M) ) THEN
               ISYDEF = MULD2H(ISYML,ISYMIM)
               DO ISYMF = 1,NSYM
                  ISYDE = MULD2H(ISYDEF,ISYMF)
                  DO F = 1,NVIR(ISYMF)
C
                     NFI   = IT1AM(ISYMF,ISYMM)+NVIR(ISYMF)*(M-1)+F
C
                     ISYME = MULD2H(ISYDE,ISYMD)
                     ISYMEM = MULD2H(ISYME,ISYMM)
                     ISYMFI = MULD2H(ISYMF,ISYMM)
C
                     DO E =1,NVIR(ISYME)
C
                        KOFF1 = IT2SQ(ISYMEM,ISYMFI)
     *                        + NT1AM(ISYMEM)*(NFI-1)
     *                        + IT1AM(ISYME,ISYMM)
     *                        + NVIR(ISYME)*(M-1)
     *                        + E
C
                        TMAT(KOFF1)  = ZERO
C
                    END DO
                 END DO
              END DO
           END IF
        END DO
      END DO
C
C---------------------------------------------------------
C     If D and E are the same remove all amplitudes
C     having an F which is the same as D and E.
C---------------------------------------------------------
C
      DO ISYME = 1,NSYM
         DO E = 1,NVIR(ISYME)
            IF ((ISYMD .EQ. ISYME) .AND. (D .EQ. E) ) THEN
               ISYLMI = MULD2H(ISYMD,ISYMIM)
               DO ISYMI = 1,NSYM
                  ISYLM = MULD2H(ISYLMI,ISYMI)
                  DO I = 1,NRHF(ISYMI)
C
                  NFI   = IT1AM(ISYME,ISYMI)+NVIR(ISYME)*(I-1)+E
C
                  ISYMM = MULD2H(ISYLM,ISYML)
                  ISYMEM = MULD2H(ISYME,ISYMM)
                  ISYMFI = MULD2H(ISYME,ISYMI)
C
                  DO M =1,NRHF(ISYMM)
C
                     KOFF1 = IT2SQ(ISYMEM,ISYMFI)
     *                     + NT1AM(ISYMEM)*(NFI-1)
     *                     + IT1AM(ISYME,ISYMM)
     *                     + NVIR(ISYME)*(M-1)
     *                     + E
C
                     TMAT(KOFF1)  = ZERO
C
                    END DO
                 END DO
              END DO
           END IF
        END DO
      END DO
C
C-----------------------
C     End.
C-----------------------
C
      CALL QEXIT('T3_FORBIDDEN_DL')
C
      RETURN
      END
C  /* Deck read_t3_aibl */
      SUBROUTINE READ_T3_AIBL(LUFILE,FNFILE,ISYMT3,T2SQ,ISYMT2,
     *                        J,ISYMJ,ISYMD)

      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <cc3t3d.h>
#include <ccorb.h>
C
      CHARACTER*(*) FNFILE
C
      INTEGER LUFILE,ISYMT3,ISYMT2,ISYMJ,ISYMD
C
      INTEGER ISYMAIBJL,ISYMAIBL,ISYML,ISYMAIB,ISYMAIBJ,ISYMB
      INTEGER ISYMBJ,ISYMBL,ISYMAI,NBJ,NBL,KOFFT2,IADR
C
#if defined (SYS_CRAY)
      REAL T2SQ(*) 
#else
      DOUBLE PRECISION T2SQ(*)
#endif
C
      CALL QENTER('READ_T3_AIBL')
C
      ISYMAIBJL = MULD2H(ISYMT3,ISYMD)
      ISYMAIBL  = MULD2H(ISYMAIBJL,ISYMJ)
      DO ISYML = 1,NSYM
         ISYMAIB = MULD2H(ISYMAIBL,ISYML)
         ISYMAIBJ = MULD2H(ISYMAIB,ISYMJ)
         DO ISYMB = 1,NSYM
            ISYMBL = MULD2H(ISYMB,ISYML)
            ISYMBJ = MULD2H(ISYMB,ISYMJ)
            ISYMAI = MULD2H(ISYMAIB,ISYMB)
            DO L = 1,NRHF(ISYML)
               DO B = 1,NVIR(ISYMB)
C
                  NBL  = IT1AM(ISYMB,ISYML)+NVIR(ISYMB)*(L-1)+B
C
                  KOFFT2 = IT2SQ(ISYMAI,ISYMBL)+NT1AM(ISYMAI)*(NBL-1)+1 
C
                  NBJ  = IT1AM(ISYMB,ISYMJ)+NVIR(ISYMB)*(J-1)+B
C
                  IADR = ISWTL(ISYMAIBJ,ISYML)+NT2SQ(ISYMAIBJ)*(L-1)+
     *                    IT2SQ(ISYMAI,ISYMBJ)+NT1AM(ISYMAI)*(NBJ-1)+1

                  CALL GETWA2(LUFILE,FNFILE,T2SQ(KOFFT2),
     *                        IADR,NT1AM(ISYMAI))
C
               END DO
            END DO
         END DO
      END DO
C
      CALL QEXIT('READ_T3_AIBL')
C
      RETURN
      END
C  /* Deck cc3_xi_den_ai_t1 */
      SUBROUTINE CC3_XI_DEN_AI_T1(DIA,ISYDIA,DAB0,DIJ0,ISYDEN0,
     *                            T1R1,ISYMR1)
C
C     D_ia(TL1,TR1) = D_ia(TL1,TR1) - DAB0(ba) * TR1(bi)
C                                   + DIJ0(ij) * TR1(aj)
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYDIA,ISYDEN0,ISYMR1
      INTEGER ISYMA,ISYMI,ISYMB,ISYMJ
      INTEGER KOFFDBA,KOFFDIJ,KOFFDIA,KOFFTBI,KOFFTAJ
      INTEGER NVIRA,NVIRB,NRHFI
C
#if defined (SYS_CRAY)
      REAL DIA(*),DAB0(*),DIJ0(*),T1R1(*)
      REAL ONE
#else
      DOUBLE PRECISION DIA(*),DAB0(*),DIJ0(*),T1R1(*)
      DOUBLE PRECISION ONE
#endif
C
      PARAMETER (ONE = 1.0D0)
C
      CALL QENTER('CC3XIT1')
C
C     D_ia(TL1,TR1) = D_ia(TL1,TR1) - DAB0(ba) * TR1(bi)
C                                   + DIJ0(ij) * TR1(aj)
      DO ISYMA = 1, NSYM
        ISYMI = MULD2H(ISYDIA,ISYMA)
        ISYMB = MULD2H(ISYDEN0,ISYMA)
        ISYMJ = MULD2H(ISYDEN0,ISYMI)
C
C       Symmetry check
C
        IF ( (MULD2H(ISYMR1,ISYMA) .NE. ISYMJ)
     *       .OR. (MULD2H(ISYMR1,ISYMI) .NE. ISYMB) ) THEN
           WRITE(LUPRI,*)' Symmetry mismatch in CC3_XI_DEN_AI_T1 '
           WRITE(LUPRI,*)' ISYMR1,ISYMI,MULD2H(ISYMR1,ISYMI) ',
     *                     ISYMR1,ISYMI,MULD2H(ISYMR1,ISYMI)
           WRITE(LUPRI,*)' NOT equal ISYMB ', ISYMB
           WRITE(LUPRI,*)' OR ISYMR1,ISYMA,MULD2H(ISYMR1,ISYMA) ',
     *                        ISYMR1,ISYMA,MULD2H(ISYMR1,ISYMA)
           WRITE(LUPRI,*)' NOT equal ISYMJ ', ISYMJ
           CALL QUIT('Symmetry mismatch in CC3_XI_DEN_AI_T1')
        END IF
        

        KOFFDBA = IMATAB(ISYMB,ISYMA) + 1
        KOFFDIJ = IMATIJ(ISYMI,ISYMJ) + 1
        KOFFDIA = IT1AM(ISYMA,ISYMI)  + 1
        KOFFTBI = IT1AM(ISYMB,ISYMI)  + 1
        KOFFTAJ = IT1AM(ISYMA,ISYMJ)  + 1

        NVIRA   = MAX(NVIR(ISYMA),1)
        NVIRB   = MAX(NVIR(ISYMB),1)
        NRHFI   = MAX(NRHF(ISYMI),1)

        CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMB),
     &             -ONE,DAB0(KOFFDBA),NVIRB,T1R1(KOFFTBI),NVIRB,
     &              ONE,DIA(KOFFDIA),NVIRA)

        CALL DGEMM('N','T',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     &              ONE,T1R1(KOFFTAJ),NVIRA,DIJ0(KOFFDIJ),NRHFI,
     &              ONE,DIA(KOFFDIA),NVIRA)

      END DO
C
      CALL QEXIT('CC3XIT1')
C
      RETURN
      END
C /* t2_aibj_exc */
      SUBROUTINE T2_AIBJ_EXC(XIJMJI,XAIBJ,ISYAIBJ,FAIBJ,FAJBI)
C
C     Calculate: XIJMJI = FAIBJ*XAIBJ(aibj) + FAJBI*XAIBJ(ajbi)
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYAIBJ
      INTEGER ISYMBJ,ISYMAI,ISYMJ,ISYMB,ISYMI,ISYMA,ISYMAJ,ISYMBI
      INTEGER NBJ,NBI
      INTEGER KAIBJ,KAJBI
C
#if defined (SYS_CRAY)
      REAL XIJMJI(*),XAIBJ(*)
      REAL FAIBJ,FAJBI
#else
      DOUBLE PRECISION XIJMJI(*),XAIBJ(*)
      DOUBLE PRECISION FAIBJ,FAJBI
#endif
C
      CALL QENTER('T2_AIBJ_EXC')
C
      DO ISYMBJ = 1,NSYM
         ISYMAI = MULD2H(ISYAIBJ,ISYMBJ)
         DO ISYMJ = 1,NSYM
            ISYMB = MULD2H(ISYMBJ,ISYMJ)
            DO ISYMI = 1,NSYM
               ISYMA = MULD2H(ISYMAI,ISYMI)
               ISYMAJ = MULD2H(ISYMA,ISYMJ)
               ISYMBI = MULD2H(ISYMB,ISYMI)
               DO J = 1,NRHF(ISYMJ)
                  DO B = 1,NVIR(ISYMB)
C
                     NBJ=IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J-1) + B
C
                     DO I = 1,NRHF(ISYMI)
C
                        NBI=IT1AM(ISYMB,ISYMI) + NVIR(ISYMB)*(I-1) + B
C
                        DO A = 1,NVIR(ISYMA)
C
                           KAIBJ = IT2SQ(ISYMAI,ISYMBJ)
     *                           + NT1AM(ISYMAI)*(NBJ-1)
     *                           + IT1AM(ISYMA,ISYMI)
     *                           + NVIR(ISYMA)*(I-1)
     *                           + A 
                           KAJBI = IT2SQ(ISYMAJ,ISYMBI)
     *                           + NT1AM(ISYMAJ)*(NBI-1)
     *                           + IT1AM(ISYMA,ISYMJ)
     *                           + NVIR(ISYMA)*(J-1)
     *                           + A 
C
                           XIJMJI(KAIBJ) = FAIBJ*XAIBJ(KAIBJ)
     *                                   + FAJBI*XAIBJ(KAJBI)
C
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
C
C     End.
C
      CALL QEXIT('T2_AIBJ_EXC')
C
      RETURN
      END
C  /* Deck intocc_t3barx_jk */
      SUBROUTINE INTOCC_T3BARX_JK(LUTOC,FNTOC,ISYHAM,
     *                         LAMHY,ISYMY,ISYINTY,
     *                         W3BXOGX1,W3BXOLX1,INTERNAL,
     *                         W3BXOGX2,W3BXOLX2,
     *                         WORK,LWORK)
*
***********************************************************************
*                                                                     *
*     Construct occupied integrals which are required to calculate    *
*     t3bar_X multipliers in JK-LOOP (for constructing A density)     *
*     (otherwise use intocc_t3barx routine)                           *
*                                                                     *
*     g(ia|j k-)    and     L(ia|j k-)                                *
*                                                                     *
*     ISYHAM  - symmetry of Hamiltonian                               *
*     ISYMY   - symmetry of operator in LambdaH_Y transformation      *
*     ISYINTY - symmetry of LambdaH_Y-transformed integrals           *
*                                                                     *
*     IF INTERNAL = .TRUE. W3BXOGX1,W3BXOLX1 integrals are used only  *
*     as internal intermmediates and do not need to be declarted      *
*     outside.                                                        *
*                                                                     *
*     Filip Pawlowski, 15-Apr-2003, Aarhus                            *
*                                                                     *
***********************************************************************
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      LOGICAL INTERNAL
C
      INTEGER ISYHAM,ISYMY,ISYINTY,LWORK,LUTOC
      INTEGER KINTOC,KEND1,LWRK1,IOFF
      INTEGER KW3BXOGX1,KW3BXOLX1
C
      CHARACTER*(*) FNTOC
C
#if defined (SYS_CRAY) 
      REAL LAMHY(*),W3BXOGX1(*),W3BXOGX2(*)
      REAL W3BXOLX1(*),W3BXOLX2(*),WORK(LWORK)
#else
      DOUBLE PRECISION LAMHY(*),W3BXOGX2(*)
      DOUBLE PRECISION W3BXOGX1(*),W3BXOLX1(*),W3BXOLX2(*)
      DOUBLE PRECISION WORK(LWORK)
#endif
C
      CALL QENTER('INTOCC_T3BARX_JK')
C
C--------------------------
C     Initial symmetry check
C--------------------------
C
      IF (ISYINTY .NE. MULD2H(ISYMY,ISYHAM) ) THEN
         WRITE(LUPRI,*)'ISYINTY = ', ISYINTY
         WRITE(LUPRI,*)'should be equal to '
         WRITE(LUPRI,*)'MULD2H(ISYMY,ISYHAM) = ',MULD2H(ISYMY,ISYHAM)
         WRITE(LUPRI,*)
         WRITE(LUPRI,*)'ISYMY,ISYHAM ',ISYMY,ISYHAM
         CALL QUIT('Symmetry mismatch in INTOCC_T3BARX_JK')
      END IF
C
C---------------
C     Allocation
C---------------
C
      KINTOC = 1
      KEND1  = KINTOC + NTOTOC(ISYHAM)
      LWRK1  = LWORK  - KEND1
      IF (INTERNAL) THEN
         KW3BXOGX1 = KEND1
         KW3BXOLX1 = KW3BXOGX1 + NTRAOC(ISYINTY)
         KEND1     = KW3BXOLX1 + NTRAOC(ISYINTY)
         LWRK1     = LWORK  - KEND1
      END IF
C      
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in INTOCC_T3BARX_JK')
      END IF
C
      CALL DZERO(WORK(KINTOC),NTOTOC(ISYHAM))
      IF (INTERNAL) THEN
         CALL DZERO(WORK(KW3BXOGX1),NTRAOC(ISYINTY))
         CALL DZERO(WORK(KW3BXOLX1),NTRAOC(ISYINTY))
      END IF
C
C-------------------------------
C     Read in occupied integrals 
C-------------------------------
C
      IOFF = 1
      IF (NTOTOC(ISYHAM) .GT. 0) THEN
         CALL GETWA2(LUTOC,FNTOC,WORK(KINTOC),IOFF,NTOTOC(ISYHAM))
      ENDIF
C
C-----------------------------
C     LambdaH_Y transformation
C-----------------------------
C
      IF (INTERNAL) THEN
         CALL INTOCC_T3BARX2(WORK(KINTOC),LAMHY,ISYMY,ISYINTY,
     *                       WORK(KW3BXOGX1),WORK(KW3BXOLX1),
     *                       WORK(KEND1),LWRK1)
      ELSE 
         CALL INTOCC_T3BARX2(WORK(KINTOC),LAMHY,ISYMY,ISYINTY,W3BXOGX1,
     *                       W3BXOLX1,WORK(KEND1),LWRK1)
      END IF
C
C-------------------
C     Sort 
C-------------------
C
      IF (INTERNAL) THEN
         CALL CCFOP_SORT(WORK(KW3BXOGX1),W3BXOGX2,ISYINTY,1)
         CALL CCFOP_SORT(WORK(KW3BXOLX1),W3BXOLX2,ISYINTY,1)
      ELSE 
         CALL CCFOP_SORT(W3BXOGX1,W3BXOGX2,ISYINTY,1)
         CALL CCFOP_SORT(W3BXOLX1,W3BXOLX2,ISYINTY,1)
      END IF
C
C-------------
C     End
C-------------
C
      CALL QEXIT('INTOCC_T3BARX_JK')
C
      RETURN
      END
C  /* Deck mmatx */
      SUBROUTINE MMATX(YMMAT,LISTL,
     *                LUWBAR,FNWBAR,ISWBAR,
     *                LUT2Y,FNT2Y,IST2Y,
     *                ISYMD,D,
     *                WORK,LWORK)
*
***************************************************************************
* Prupose: calculate M intermediate:
*
* y^M(fn,im) <- M(fn,im)+ sum_em y^t^DL(ei) tbar^DL(em,fn)
*
* F. Pawlowski, 20-Oct-2003, Aarhus.
***************************************************************************
*
      IMPLICIT NONE
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <cc3t3d.h>
C
      CHARACTER LISTL*3
      CHARACTER FNWBAR*(*),FNT2Y*(*)
C
      INTEGER LUWBAR,LUT2Y
      INTEGER ISWBAR,IST2Y
      INTEGER ISYMD,LWORK
C
      INTEGER ISYML,ISYMDL,IST2YDL,ISWBARDL,KT2Y,KWBAR,KEND1,LWRK1
      INTEGER NDL,IADR,ISYMF,ISYMN,ISYMFN,ISYMEM,ISYMM,ISYME,ISYMI
      INTEGER ISYMIM
      INTEGER KOFF1,KFN,KOFF2,KOFF3,NVIRE,NRHFI
#if defined (SYS_CRAY)
      REAL YMMAT(*)
      REAL ONE,FAC
      REAL WORK(LWORK)
#else
      DOUBLE PRECISION YMMAT(*)
      DOUBLE PRECISION ONE,FAC
      DOUBLE PRECISION WORK(LWORK)
#endif
C
      PARAMETER (ONE = 1.0D0)

C
      CALL QENTER('MMATX')
C
      !check the list and set the factor
      IF (LISTL(1:3) .EQ. 'L0 ') THEN
         FAC = ONE
      ELSE IF ( (LISTL(1:3) .EQ. 'L1 ') .OR. (LISTL(1:3) .EQ. 'M1 ') 
     *         .OR. (LISTL(1:3) .EQ. 'N2 ')
     *         .OR. (LISTL(1:3) .EQ. 'LE ') ) THEN
         FAC = -ONE
      ELSE
         WRITE(LUPRI,*)'LISTL = ', LISTL(1:3)
         WRITE(LUPRI,*)'Case not implemented in MMATX'
         CALL QUIT('Wrong LISTL in MMATX')
      END IF
C
      DO ISYML = 1,NSYM
C
         ISYMDL   = MULD2H(ISYMD,ISYML)
         IST2YDL  = MULD2H(IST2Y,ISYMDL)
         ISWBARDL = MULD2H(ISWBAR,ISYMDL)
C
         KT2Y = 1
         KWBAR = KT2Y  + NT1AM(IST2YDL)
         KEND1 = KWBAR + NT2SQ(ISWBARDL)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*)'Memory available: ', LWORK
            WRITE(LUPRI,*)'Memory needed   : ', KEND1
            CALL QUIT('Insufficient memory in MMATX (1)')
         END IF
C
         DO L = 1,NRHF(ISYML)

C           --------------------------------------------
C           Read response doubles amplitudes T2^y,DL:
C           --------------------------------------------
            NDL  = IT1AM(ISYMD,ISYML) + NVIR(ISYMD)*(L-1) + D
            IADR = IT2SQ(IST2YDL,ISYMDL) + NT1AM(IST2YDL)*(NDL-1) + 1

            CALL GETWA2(LUT2Y,FNT2Y,WORK(KT2Y),IADR,NT1AM(IST2YDL))

C           --------------------------------------------
C           Read WMAT from file and generate WMAT-tilde:
C           --------------------------------------------
            IADR = ISWTL(ISWBARDL,ISYML) + NT2SQ(ISWBARDL)*(L-1) + 1
            CALL GETWA2(LUWBAR,FNWBAR,WORK(KWBAR),IADR,NT2SQ(ISWBARDL))

            IF ( (LISTL(1:3) .EQ. 'L1 ') .OR. (LISTL(1:3) .EQ. 'M1 ') 
     *         .OR. (LISTL(1:3) .EQ. 'N2 ')
     *         .OR. (LISTL(1:3) .EQ. 'LE ') ) THEN
               CALL CC_T2MOD(WORK(KWBAR),ISWBARDL,ONE)
            END IF

C           -------------------------------------------------------
C           y^M(imfn) <- M(imfn)+ sum_em y^t^DL(ei) tbar^DL(em,fn)
C           -------------------------------------------------------

            DO ISYMF = 1,NSYM
               DO ISYMN = 1,NSYM
                  ISYMFN = MULD2H(ISYMF,ISYMN)
                  ISYMEM = MULD2H(ISWBARDL,ISYMFN)
C
                  DO ISYMM  = 1,NSYM
                     ISYME  = MULD2H(ISYMEM,ISYMM)
                     ISYMI  = MULD2H(IST2YDL,ISYME)
                     ISYMIM = MULD2H(ISYMI,ISYMM)
C
                     DO N = 1,NRHF(ISYMN)
                        DO F = 1,NVIR(ISYMF)
                           KOFF1 = KT2Y 
     *                           + IT1AM(ISYME,ISYMI)
C
                           KFN    = IT1AM(ISYMF,ISYMN)+NVIR(ISYMF)*(N-1)
     *                            + F 
C
                           KOFF2 = KWBAR
     *                           + IT2SQ(ISYMEM,ISYMFN)
     *                           + NT1AM(ISYMEM)*(KFN-1)
     *                           + IT1AM(ISYME,ISYMM)
C
                           KOFF3 = 1
     *                           + ISAIKL(ISYMFN,ISYMIM)
     *                           + NMATIJ(ISYMIM)*(KFN-1)
     *                           + IMATIJ(ISYMI,ISYMM)
C
                           NVIRE  = MAX(NVIR(ISYME),1)
                           NRHFI  = MAX(NRHF(ISYMI),1)
C
                           CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMM),
     *                                NVIR(ISYME),FAC,WORK(KOFF1),NVIRE,
     *                                WORK(KOFF2),NVIRE,
     *                                ONE,YMMAT(KOFF3),NRHFI)
C
                        END DO
                     END DO
                  END DO
               END DO
            END DO
C
         END DO !L
      END DO !ISYML
C
      CALL QEXIT('MMATX')
C
      RETURN
      END
C  /* Deck mmatxl1 */
      SUBROUTINE MMATXL1(YMMAT,LISTL,LISTR,IDLSTR,
     *                LUWBAR,FNWBAR,ISWBAR,
     *                ISYMD,D,
     *                WORK,LWORK)
*
***************************************************************************
* Prupose: calculate a special contribution to M intermediate
* (i.e. the one which occurs when LISTL = 'L1 '):
*
* y^M(fn,im) <- M(fn,im)+ sum_em y^T2(dl,ei) tbar(dl,em,fn)
*
* for 'L1 ' (first-order multipliers) we have:
* tbar(dl,em,fn) = W^ef(dlnm) + W^df(emnl) + W^de(fnml)
*
* The contractions of last two W's (W^df and W^de) with T2(dl,ei) are
* corried out in MMATX routine.
*
* In THIS ROUTINE we calculate the contraction of the first W (W^ef)
* with T2(dl,ei):
*
* W^ef(dlnm)*T2(dl,ei) 
* = W^df(elnm)*T2(el,di)     -> M(fn,im).
*
* This is calculated as:
*
*  M(fn,iM) <- M(fn,iM) + I^DM(el,fn) * I^D(el,i)
*
* F. Pawlowski, 20-Oct-2003, Aarhus.
***************************************************************************
*
      IMPLICIT NONE
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <cc3t3d.h>
#include <ccr1rsp.h>
#include <dummy.h>
C
      CHARACTER LISTL*3, LISTR*3
      CHARACTER FNWBAR*(*)
C
      INTEGER LUWBAR
      INTEGER ISWBAR
      INTEGER ISYMD,LWORK
      INTEGER IDLSTR
C
      INTEGER ISYMR1,KT2R1,KEND1,LWRK1
      INTEGER IOPT,ISYMM,ISYMDM,ISWBARDM,KWBAR,KEND2,LWRK2
      INTEGER IADR,ISYMFN,ISYMEL,ISYMELI,ISYMI,ISYMIM,KOFF1,KOFF2
      INTEGER NMI,KOFF3,NEL,NFN
      integer kim,kfn,Kel,NI
      INTEGER ISYMF,ISYMN
#if defined (SYS_CRAY)
      REAL YMMAT(*)
      REAL ONE
      REAL WORK(LWORK)
      real xtmp
#else
      DOUBLE PRECISION YMMAT(*)
      DOUBLE PRECISION ONE
      DOUBLE PRECISION WORK(LWORK)
      double precision xtmp
#endif
C
      PARAMETER (ONE = 1.0D0)
C
      CALL QENTER('MMATXL1')
C
      !check the list 
      IF ( (LISTL(1:3) .NE. 'L1 ') .AND. (LISTL(1:3) .NE. 'M1 ') 
     *     .AND. (LISTL(1:3) .NE. 'N2 ')
     *     .AND. (LISTL(1:3) .NE. 'LE ') ) THEN
         WRITE(LUPRI,*)'LISTL = ', LISTL(1:3)
         WRITE(LUPRI,*)'MMATXL1 routine is designed for LISTL = '
     *   //'L1 or M1 or N2 or LE '
         CALL QUIT('Wrong LISTL in MMATXL1')
      END IF
C
      !get symmetry of right list (doubles amplitudes)
      ISYMR1  = ISYLRT(IDLSTR)
C
      KT2R1 = 1
      KEND1 = KT2R1 + NT2SQ(ISYMR1)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*)'Memory available: ', LWORK
         WRITE(LUPRI,*)'Memory needed   : ', KEND1
         CALL QUIT('Insufficient memory in MMATXL1 (1)')
      END IF
C     --------------------
C     Construct T2TP(elid)
C     --------------------
      IOPT  = 2
      CALL GET_T1_T2(IOPT,.TRUE.,DUMMY,WORK(KT2R1),LISTR,
     *               IDLSTR,ISYMR1,WORK(KEND1),LWRK1)
C
      DO ISYMM = 1,NSYM
C
         ISYMDM   = MULD2H(ISYMD,ISYMM)
         ISWBARDM = MULD2H(ISWBAR,ISYMDM)
C
         KWBAR = KEND1
         KEND2 = KWBAR + NT2SQ(ISWBARDM)
         LWRK2 = LWORK - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*)'Memory available: ', LWORK
            WRITE(LUPRI,*)'Memory needed   : ', KEND2
            CALL QUIT('Insufficient memory in MMATXL1 (2)')
         END IF
C
         DO M = 1,NRHF(ISYMM)

C           --------------------
C           Read WMAT from file 
C           --------------------
            IADR = ISWTL(ISWBARDM,ISYMM) + NT2SQ(ISWBARDM)*(M-1) + 1
            CALL GETWA2(LUWBAR,FNWBAR,WORK(KWBAR),IADR,NT2SQ(ISWBARDM))
C
C           We have now sitting the following arrays:
C           KT2R1(elid) of symmetry ISYMR1
C           KWBAR^DM(el,fn) of symmetry ISWBARDM
C

C           ----------------------------------------------
C           Calculate:
C           M(iM,fn) <- M(iM,fn) +  I^D(el,i) * I^DM(el,fn)
C           ----------------------------------------------

            DO ISYMN = 1,NSYM
              DO ISYMF = 1,NSYM
               ISYMFN = MULD2H(ISYMF,ISYMN)
               ISYMEL = MULD2H(ISWBARDM,ISYMFN)
               ISYMELI = MULD2H(ISYMR1,ISYMD)
               ISYMI   = MULD2H(ISYMELI,ISYMEL)
               ISYMIM  = MULD2H(ISYMI,ISYMM)
C
               DO N = 1,NRHF(ISYMN)
                 DO F = 1,NVIR(ISYMF)
C
                  KOFF1 = KT2R1
     *                  + IT2SP(ISYMELI,ISYMD)
     *                  + NCKI(ISYMELI)*(D-1)
     *                  + ISAIK(ISYMEL,ISYMI)
C
                  KFN = IT1AM(ISYMF,ISYMN) + NVIR(ISYMF)*(N-1) + F
C
                  KOFF2 = KWBAR
     *                  + IT2SQ(ISYMEL,ISYMFN)
     *                  + NT1AM(ISYMEL)*(KFN-1)
C
                  KOFF3  = 1 + ISAIKL(ISYMFN,ISYMIM) 
     *                   + NMATIJ(ISYMIM)*(KFN-1)
     *                   + IMATIJ(ISYMI,ISYMM) + NRHF(ISYMI)*(M-1)

C
                  NEL = MAX(NT1AM(ISYMEL),1)
                  NI = MAX(NRHF(ISYMI),1)
C
                  CALL DGEMM('T','N',NRHF(ISYMI),1,
     *                       NT1AM(ISYMEL),-ONE,WORK(KOFF1),NEL,
     *                       WORK(KOFF2),NEL,ONE,YMMAT(KOFF3),NI)
C
               END DO
              END DO
             END DO
            END DO

         END DO !M
      END DO !ISYMM

C
      CALL QEXIT('MMATXL1')
C
      RETURN
      END

