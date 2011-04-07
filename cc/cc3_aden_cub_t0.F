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
C  /* Deck cc3_aden_cub_t0 */
      SUBROUTINE CC3_ADEN_CUB_T0(LISTL,IDLSTL,LISTR,IDLSTR,
     *                            XLAMDP0,XLAMDH0,FOCK0,
     *                            DIJ,DAB,ISYDEN,
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
      LOGICAL DO_YMMAT
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
#if defined (SYS_CRAY)
      REAL XLAMDP0(*),XLAMDH0(*),FOCK0(*)
      REAL DAB(*),DIJ(*)
      REAL YMMAT(*)
      REAL WORK(LWORK)
#else
      DOUBLE PRECISION XLAMDP0(*),XLAMDH0(*),FOCK0(*)
      DOUBLE PRECISION DAB(*),DIJ(*)
      DOUBLE PRECISION YMMAT(*)
      DOUBLE PRECISION WORK(LWORK)
#endif
C
      CALL QENTER('CC3DENCB')

      ISYM0 = 1
C
      KT1AMP = 1
      KLAMP0 = KT1AMP + NT1AM(ISYM0)
      KLAMH0 = KLAMP0 + NLAMDT
      KEND1 = KLAMH0 + NLAMDT
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CC3_ADEN_CUB (1)')
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
c     IF (IPRINT .GT. 55) THEN
c        WRITE(LUPRI,*)'DAB density before CC3_ADENVIR_CUB '
c        CALL PRINT_MATAB(DAB,ISYDEN)
c        WRITE(LUPRI,*)'DIJ density after CC3_ADENVIR_CUB '
c        CALL PRINT_MATIJ(DIJ,ISYDEN)
c     END IF
      CALL CC3_ADENVIR_CUB_T0(DIJ,DAB,ISYDEN,
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
         WRITE(LUPRI,*)'DAB density after CC3_ADENVIR_CUB '
         CALL PRINT_MATAB(DAB,ISYDEN)
         WRITE(LUPRI,*)'DIJ density after CC3_ADENVIR_CUB '
         CALL PRINT_MATIJ(DIJ,ISYDEN)
      END IF
C
      CALL CC3_ADENOCC_CUB_T0(LISTL,IDLSTL,LISTR,IDLSTR,
     *                            WORK(KLAMP0),WORK(KLAMH0),FOCK0,
     *                            DIJ,DAB,ISYDEN,
     *                            WORK(KEND1),LWRK1,
     *                            LUDELD,FNDELD,LUCKJD,FNCKJD,LUDKBC,
     *                            FNDKBC,LUTOC,FNTOC,LU3VI,FN3VI,
     *                            LUDKBC3,FNDKBC3,LU3FOP,FN3FOP,
     *                            LU3FOPX,FN3FOPX,
     *                            LU3FOP2X,FN3FOP2X)
C
      IF (IPRINT .GT. 55) THEN
         WRITE(LUPRI,*)'DAB density after CC3_ADENOCC_CUB '
         CALL PRINT_MATAB(DAB,ISYDEN)
         WRITE(LUPRI,*)'DIJ density after CC3_ADENOCC_CUB '
         CALL PRINT_MATIJ(DIJ,ISYDEN)
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
      CALL QEXIT('CC3DENCB')
C
      RETURN
      END
C  /* Deck cc3_adenocc_cub_t0 */
      SUBROUTINE CC3_ADENOCC_CUB_t0(LISTL,IDLSTL,LISTR,IDLSTR,
     *                            XLAMDP0,XLAMDH0,FOCK0,
     *                            DIJ,DAB,ISYDEN,
     *                            WORK,LWORK,
     *                            LUDELD,FNDELD,LUCKJD,FNCKJD,LUDKBC,
     *                            FNDKBC,LUTOC,FNTOC,LU3VI,FN3VI,
     *                            LUDKBC3,FNDKBC3,LU3FOP,FN3FOP,
     *                            LU3FOPX,FN3FOPX,
     *                            LU3FOP2X,FN3FOP2X)
*---------------------------------------------------------------------*
*
*    The same as CC3_ADENOCC_CUB, but simplified since tbarX is
*    replaced by tbar0; thus some terms are droipped.
*
*    Written by Filip Pawlowski, Fall 2003, Aarhus
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
#include <ccr2rsp.h>
C
      INTEGER ISYM0
      PARAMETER(ISYM0 = 1)
C
      CHARACTER CDUMMY*1 
      PARAMETER (CDUMMY = ' ')
C
      CHARACTER*14 FN3SRTR, FNCKJDRZ, FNDELDRZ, FNDKBCRZ
      PARAMETER(FN3SRTR = 'CCSDT_FBMAT1_Z',FNCKJDRZ = 'CCSDT_FBMAT2_Z',
     *          FNDELDRZ = 'CCSDT_FBMAT3_Z',FNDKBCRZ = 'CCSDT_FBMAT4_Z')
      INTEGER LU3SRTR, LUCKJDRZ, LUDELDRZ, LUDKBCRZ
C
      CHARACTER*14 FNCKJDRU, FNDELDRU, FNDKBCRU
      PARAMETER(FNCKJDRU = 'CCSDT_FBMAT2_U',
     *          FNDELDRU = 'CCSDT_FBMAT3_U',FNDKBCRU = 'CCSDT_FBMAT4_U')
      INTEGER LUCKJDRU, LUDELDRU, LUDKBCRU
C
      INTEGER ISYDEN,IDLSTL,IDLSTR,IDLSTL0,LWORK
      INTEGER LUTOC, LU3VI, LUDKBC3, LU3FOPX, LU3FOP2X, LU3FOP
      INTEGER LUDKBC,LUDELD,LUCKJD
      INTEGER ISYML0,ISYML1,ISYMRZ,ISINT1,ISINT2,ISINT1RZ,ISYFCKL1R
      INTEGER ISYMK,ISYML,ISYMT3,ISYMKL,ISYT30KL
      INTEGER IOPT,LENGTH
      INTEGER KFOCKD,KFCKBA,KT2TP,KL1AM,KL2TP,KEND0,LWRK0
      INTEGER KT1RZ,KT2RZ,KFOCKRZ,KEND1,LWRK1
      INTEGER KXIAJB,KT3BOG1,KT3BOL1,KT3BOG2,KT3BOL2,KT3OG1,KT3OG2
      INTEGER KT30KL
      INTEGER KT3VIJG1
      INTEGER ISYMT3B,ISYT3B0KL,ISYW3BXKL
      INTEGER KXGADCK,KXLADCK
      INTEGER KT3B0KL,KW3BXKL,ISYMW3BX
      INTEGER ISYMTETAZ,ISTETAZKL
      INTEGER KTETAXKL
      INTEGER IDLSTL1R,ISYML1R
      INTEGER ISINT2L1R
C
      INTEGER IDLSTZU,IDLSTRZ,IDLSTRU,ISYMRU
      INTEGER KFOCKRU,ISYMZU,ISYMTETAU,ISYMTETAZU,ISTETAUKL,ISTETAZUKL
      INTEGER MAXX1
      INTEGER KABCI,KT3B0KL3
      INTEGER KFCKZUV,KFCKUZV,KLAMDPZ,KLAMDHTMP,KLAMDPU
C
      INTEGER KGBCDK
      INTEGER KT1RU,KT2RU
C
      INTEGER ISINT2RZ,ISINT1RU,ISINT2RU,KT3OG2Z
      INTEGER KT3OG2U,KGBCDKZ,KGBCDKU
      INTEGER KEND2,LWRK2
C
      INTEGER IR1TAMP
C
      CHARACTER LISTL*3, LISTR*3, LISTL0*3, LISTL1R*3
      CHARACTER*(*) FNTOC, FN3VI, FNDKBC3, FN3FOPX, FN3FOP2X, FN3FOP
      CHARACTER*(*) FNDKBC,FNDELD,FNCKJD
      CHARACTER LABELL1*8,LABELRZ*8,LABELRU*8
C
      CHARACTER LISTRZ*3,LISTRU*3
C
      LOGICAL   LOCDBG,LORXL1
      PARAMETER (LOCDBG = .FALSE.)
      LOGICAL   LORXRZ,LORXRU
C
C
      integer kx3am
C
#if defined (SYS_CRAY)
      REAL XLAMDP0(*),XLAMDH0(*),FOCK0(*) 
      REAL DAB(*),DIJ(*)
      REAL WORK(LWORK)
      REAL FREQL1,FREQRZ,FREQL1R,FREQRU,FREQZU
      REAL DDOT,XNORMVAL,ONE
#else
      DOUBLE PRECISION XLAMDP0(*),XLAMDH0(*),FOCK0(*) 
      DOUBLE PRECISION DAB(*),DIJ(*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION FREQL1,FREQRZ,FREQL1R,FREQRU,FREQZU
      DOUBLE PRECISION DDOT,XNORMVAL,ONE
#endif
C
      PARAMETER (ONE = 1.0D0)
C
      CALL QENTER('CC3AOCB')
C
C--------------------------------
C     Open temporary files
C--------------------------------
C
      LU3SRTR  = -1
      LUCKJDRZ = -1
      LUDELDRZ = -1
      LUDKBCRZ = -1
C
      CALL WOPEN2(LU3SRTR,FN3SRTR,64,0)
      CALL WOPEN2(LUCKJDRZ,FNCKJDRZ,64,0)
      CALL WOPEN2(LUDELDRZ,FNDELDRZ,64,0)
      CALL WOPEN2(LUDKBCRZ,FNDKBCRZ,64,0)
C
      LUCKJDRU = -1
      LUDELDRU = -1
      LUDKBCRU = -1
C
      CALL WOPEN2(LUCKJDRU,FNCKJDRU,64,0)
      CALL WOPEN2(LUDELDRU,FNDELDRU,64,0)
      CALL WOPEN2(LUDKBCRU,FNDKBCRU,64,0)
C
C------------------------------------------------------------
C     some initializations:
C------------------------------------------------------------
C
      ISINT1 = 1
      ISINT2 = 1
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

         IF (LORXL1) CALL QUIT('NO ORBITAL RELAX. IN CC3_ADENOCC_CUB')

        LISTL1R  = 'R1 '
        IDLSTL1R = IR1TAMP(LABELL1,LORXL1,FREQL1,ISYML1)
        ! get symmetry and frequency from common blocks
        ! defined in ccl1rsp.h
        ISYML1R  = ISYLRT(IDLSTL1R)
        FREQL1R  = FRQLRT(IDLSTL1R)
C
        IF (FREQL1R .NE. FREQL1) THEN
           WRITE(LUPRI,*)'FREQL1R: ', FREQL1R
           WRITE(LUPRI,*)'FREQL1: ', FREQL1
           CALL QUIT('Frequency mismatch in CC3_ADENOCC_CUB')
        END IF

      ELSE 
         CALL QUIT('Unknown left list in CC3_ADENOCC_CUB')
      END IF

      IF (LISTR(1:3).EQ.'R2 ') THEN
         IDLSTZU = IDLSTR 
         ! get symmetry, frequency and integral label for right list 
         ! from common blocks defined in ccr1rsp.h
         LISTRZ  = 'R1 '
         LABELRZ = LBLR2T(IDLSTZU,1)
         ISYMRZ  = ISYR2T(IDLSTZU,1)
         FREQRZ  = FRQR2T(IDLSTZU,1)
         LORXRZ  = LORXR2T(IDLSTZU,1)
         IDLSTRZ = IR1TAMP(LABELRZ,LORXRZ,FREQRZ,ISYMRZ)

         LISTRU  = 'R1 '   
         LABELRU = LBLR2T(IDLSTZU,2)
         ISYMRU  = ISYR2T(IDLSTZU,2)
         FREQRU  = FRQR2T(IDLSTZU,2)
         LORXRU  = LORXR2T(IDLSTZU,2)
         IDLSTRU = IR1TAMP(LABELRU,LORXRU,FREQRU,ISYMRU)

C
         IF (LORXRZ.OR.LORXRU) THEN
          CALL QUIT('Orbital relaxation not allowed in CC3_ADENVIR_CUB')
         END IF
C
      ELSE
         WRITE(LUPRI,*)'LISTR = ',LISTR(1:3)
         WRITE(LUPRI,*)'CC3_ADENVIR_CUB is designed for LISTR = R2'
         CALL QUIT('Unknown right list in CC3_ADENOCC_CUB')
      END IF

      FREQZU = FREQRZ + FREQRU
      ISYMZU = MULD2H(ISYMRZ,ISYMRU)

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
      KT1RZ     = KEND0
      KT2RZ     = KT1RZ     + NT1AM(ISYMRZ)
      KFOCKRZ = KT2RZ     + NT2SQ(ISYMRZ)
      KEND1   = KFOCKRZ + N2BST(ISYMRZ)
      LWRK1   = LWORK   - KEND1
C
      KT1RU     = KEND1
      KT2RU     = KT1RU     + NT1AM(ISYMRU)
      KEND1 = KT2RU     + NT2SQ(ISYMRU)
      LWRK1   = LWORK   - KEND1
C
      KFOCKRU = KEND1
      KEND1   = KFOCKRU + N2BST(ISYMRU)
      LWRK1   = LWORK   - KEND1
C
      KFCKZUV  = KEND1 + N2BST(ISYMZU)
      KFCKUZV  = KFCKZUV + N2BST(ISYMZU)
      KEND1   = KFCKUZV + N2BST(ISYMZU)
      LWRK1   = LWORK   - KEND1
C
      KLAMDPZ = KEND1
      KLAMDPU = KLAMDPZ + NLAMDT
      KLAMDHTMP = KLAMDPU + NLAMDT
      KEND1   = KLAMDHTMP + NLAMDT
      LWRK1   = LWORK   - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CC3_ADENOCC_CUB (1)')
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
C      WRITE(LUPRI,*) 'Norm of L2TP (after reading)',
C    *    DDOT(NT2SQ(ISYML0),WORK(KL2TP),1,WORK(KL2TP),1)

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
C     Read the matrix the property integrals and trasform it to lambda 
C     basis (Z operator)
C---------------------------------------------------------------------
C        
         CALL GET_FOCKX(WORK(KFOCKRZ),LABELRZ,IDLSTRZ,ISYMRZ,XLAMDP0,
     *                  ISYM0,XLAMDH0,ISYM0,WORK(KEND1),LWRK1)
C
C---------------------------------------------------------------------
C     Read the matrix the property integrals and trasform it to lambda 
C     basis (U operator)
C---------------------------------------------------------------------
C        
         CALL GET_FOCKX(WORK(KFOCKRU),LABELRU,IDLSTRU,ISYMRU,XLAMDP0,
     *                  ISYM0,XLAMDH0,ISYM0,WORK(KEND1),LWRK1)
C
C------------------------------------------
C     Calculate the [U,T1^Z] matrix
C     Recall that we only need vir-vir block.
C------------------------------------------
C
      CALL GET_LAMBDAX(WORK(KLAMDPZ),WORK(KLAMDHTMP),LISTRZ,IDLSTRZ,
     *                 ISYMRZ,XLAMDP0,XLAMDH0,WORK(KEND1),
     *                 LWRK1)
      ! get vir-vir block U_(c-,d)
      CALL GET_FOCKX(WORK(KFCKUZV),LABELRU,IDLSTRU,ISYMRU,WORK(KLAMDPZ),
     *                  ISYMRZ,XLAMDH0,ISYM0,WORK(KEND1),LWRK1)
C
C------------------------------------------
C     Calculate the [Z,T1^U] matrix
C     Recall that we only need the vir-vir block.
C------------------------------------------
C
      CALL GET_LAMBDAX(WORK(KLAMDPU),WORK(KLAMDHTMP),LISTRU,IDLSTRU,
     *                 ISYMRU,XLAMDP0,XLAMDH0,WORK(KEND1),
     *                 LWRK1)
      ! get vir-vir block Z_(c-,d)
      CALL GET_FOCKX(WORK(KFCKZUV),LABELRZ,IDLSTRZ,ISYMRZ,WORK(KLAMDPU),
     *                  ISYMRU,XLAMDH0,ISYM0,WORK(KEND1),LWRK1)
C
C-------------------------------------
C     Read R1 and R2 amplitudes 
C-------------------------------------
C     
         IOPT  = 3
         CALL GET_T1_T2(IOPT,.TRUE.,WORK(KT1RZ),WORK(KT2RZ),LISTRZ,
     *                  IDLSTRZ,ISYMRZ,WORK(KEND1),LWRK1)
C
C-------------------------------------
C     Read R1 and R2 amplitudes 
C-------------------------------------
C     
         IOPT  = 3
         CALL GET_T1_T2(IOPT,.TRUE.,WORK(KT1RU),WORK(KT2RU),LISTRU,
     *                  IDLSTRU,ISYMRU,WORK(KEND1),LWRK1)
C
C
C----------------------------------------
C     Integrals [H,T1Z] where Z is LISTRZ
C----------------------------------------
C
      ISINT1RZ = MULD2H(ISINT1,ISYMRZ)
      ISINT2RZ = MULD2H(ISINT2,ISYMRZ)
C
      CALL CC3_BARINT(WORK(KT1RZ),ISYMRZ,XLAMDP0,
     *                XLAMDH0,WORK(KEND1),LWRK1,
     *                LU3SRTR,FN3SRTR,LUCKJDRZ,FNCKJDRZ)
C
      CALL CC3_SORT1(WORK(KEND1),LWRK1,2,ISINT1RZ,LU3SRTR,FN3SRTR,
     *               LUDELDRZ,FNDELDRZ,IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *               IDUMMY,CDUMMY)
C
      CALL CC3_SINT(XLAMDH0,WORK(KEND1),LWRK1,ISINT1RZ,
     *              LUDELDRZ,FNDELDRZ,LUDKBCRZ,FNDKBCRZ)

C
C----------------------------------------
C     Integrals [H,T1U] where U is LISTRU
C----------------------------------------
C
      ISINT1RU = MULD2H(ISINT1,ISYMRU)
      ISINT2RU = MULD2H(ISINT2,ISYMRU)
C
      CALL CC3_BARINT(WORK(KT1RU),ISYMRU,XLAMDP0,
     *                XLAMDH0,WORK(KEND1),LWRK1,
     *                LU3SRTR,FN3SRTR,LUCKJDRU,FNCKJDRU)
C
      CALL CC3_SORT1(WORK(KEND1),LWRK1,2,ISINT1RU,LU3SRTR,FN3SRTR,
     *               LUDELDRU,FNDELDRU,IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *               IDUMMY,CDUMMY)
C
      CALL CC3_SINT(XLAMDH0,WORK(KEND1),LWRK1,ISINT1RU,
     *              LUDELDRU,FNDELDRU,LUDKBCRU,FNDKBCRU)

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
           call quit('Insufficient space(kx3am) in CC3_ADENOCC_CUB (2)')
         END IF
      endif
C
C-----------------------------
C     Memory allocation.
C-----------------------------
C
C        isint1, isint2  - symmetry of integrals in standard H, transformed
C                  with LambdaH_0
C        ISINT1RZ - symmetry of integrals in standard H, transformed
C                  with LambdaH_R1

      ISINT1    = 1
      ISINT2    = 1
      ISINT1RZ   = MULD2H(ISINT1,ISYMRZ)
      ISINT2L1R = MULD2H(ISYML1R,ISINT2)
      ISYFCKL1R  = MULD2H(ISYMOP,ISYML1R)

      KXIAJB    = KEND1
      KEND1     = KXIAJB    + NT2AM(ISYM0)
C
      MAXX1 = MAX(NTRAOC(ISINT2RZ),NTRAOC(ISINT2RU))
C

      KT3BOG1   = KEND1
      KT3BOL1   = KT3BOG1   + MAX(NTRAOC(ISINT2L1R),NTRAOC(ISYM0))
      KT3BOG2   = KT3BOL1   + MAX(NTRAOC(ISINT2L1R),NTRAOC(ISYM0))
      KT3BOL2   = KT3BOG2   + NTRAOC(ISYM0)
      KT3OG1    = KT3BOL2   + NTRAOC(ISYM0)
      KT3OG2    = KT3OG1    + MAX(NTRAOC(ISINT2),MAXX1)
      KEND1   = KT3OG2    + NTRAOC(ISINT2)
C
      KT3OG2Z    = KEND1
      KEND1      = KT3OG2Z       + NTRAOC(ISINT2RZ)
      LWRK1     = LWORK     - KEND1
C
      KT3OG2U    = KEND1
      KEND1      = KT3OG2U       + NTRAOC(ISINT2RU)
      LWRK1     = LWORK     - KEND1
C
      KT3VIJG1  = KEND1
      KEND1     = KT3VIJG1  + NMAABCI(ISYM0)
      LWRK1     = LWORK     - KEND1
C
      KXGADCK   = KEND1 
      KXLADCK   = KXGADCK + NMAABCI(ISYM0)
      KEND1     = KXLADCK + NMAABCI(ISYM0)
      LWRK1     = LWORK     - KEND1
C
      KGBCDK    = KEND1
      KEND1     = KGBCDK + NMAABCI(ISYM0)
      LWRK1     = LWORK     - KEND1
C
      KGBCDKZ    = KEND1
      KEND1     = KGBCDKZ + NMAABCI(ISYMRZ)
      LWRK1     = LWORK     - KEND1
C
      KGBCDKU    = KEND1
      KEND1     = KGBCDKU + NMAABCI(ISYMRU)
      LWRK1     = LWORK     - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in CC3_ADENOCC_CUB (3)')
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
C     t3_x amplitudes
C-----------------------------------------------------------------
C
      CALL INTVIR_T3X_JK(WORK(KGBCDK),ISYM0,LUDKBC,FNDKBC,
     *                   WORK(KEND1),LWRK1)
C
      CALL INTVIR_T3X_JK(WORK(KGBCDKZ),ISYMRZ,LUDKBCRZ,FNDKBCRZ,
     *                   WORK(KEND1),LWRK1)
C
      CALL INTVIR_T3X_JK(WORK(KGBCDKU),ISYMRU,LUDKBCRU,FNDKBCRU,
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
C     t3_x amplitudes
C-----------------------------------------------------------------
C
      CALL INTOCC_T30(LUCKJDRZ,FNCKJDRZ,XLAMDP0,ISINT2RZ,WORK(KT3OG1),
     *                WORK(KT3OG2Z),WORK(KEND1),LWRK1)
C
      CALL INTOCC_T30(LUCKJDRU,FNCKJDRU,XLAMDP0,ISINT2RU,WORK(KT3OG1),
     *                WORK(KT3OG2U),WORK(KEND1),LWRK1)
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
C----------------------------
C     Loop over K
C----------------------------
C
      ISYMW3BX = MULD2H(ISYM0,ISYML1)
      ISYMTETAZ = MULD2H(ISYM0,ISYMRZ)
      ISYMTETAU = MULD2H(ISYM0,ISYMRU)
      ISYMTETAZU = MULD2H(ISYM0,ISYMZU)
      DO ISYMK = 1,NSYM

         DO K = 1,NRHF(ISYMK)
C
            DO ISYML = 1,NSYM
C
               ISYMKL = MULD2H(ISYMK,ISYML)
               ISYT30KL = MULD2H(ISYMKL,ISYMT3)
               ISYT3B0KL = MULD2H(ISYMKL,ISYMT3B)
               ISYW3BXKL  = MULD2H(ISYMKL,ISYMW3BX)
               ISTETAZKL  = MULD2H(ISYMKL,ISYMTETAZ)
               ISTETAUKL  = MULD2H(ISYMKL,ISYMTETAU)
               ISTETAZUKL  = MULD2H(ISYMKL,ISYMTETAZU)
C
               MAXX1 = MAX(NMAABCI(ISTETAZKL),NMAABCI(ISTETAUKL))
C
               KT30KL = KEND1
               KT3B0KL  = KT30KL + NMAABCI(ISYT30KL)
               KW3BXKL  = KT3B0KL + MAX( NMAABCI(ISYT3B0KL),MAXX1)
               KTETAXKL = KW3BXKL 
     *                  + MAX(NMAABCI(ISYW3BXKL),NMAABCI(ISTETAZKL))
               KEND2   = KTETAXKL + MAX(MAXX1,NMAABCI(ISTETAZUKL))
               LWRK2  = LWORK  - KEND2
C
               KABCI = KEND2
               KT3B0KL3 = KABCI + NMAABCI(ISTETAUKL)
               KEND2 = KT3B0KL3 + NMAABCI(ISYT3B0KL)
               LWRK2  = LWORK  - KEND2

               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Memory available : ',LWORK
                  WRITE(LUPRI,*) 'Memory needed    : ',KEND2
                  CALL QUIT('Insufficient space in CC3_ADENOCC_CUB (4)')
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
                    WRITE(LUPRI,*)'NORM OF KT30KL IN CC3_ADENOCC_CUB ', 
     *                             XNORMVAL
                  END IF
C
C-------------------------------------------
C                 Get T3BAR0^KL multipliers
C-------------------------------------------
C
                  CALL DZERO(WORK(KT3B0KL3),NMAABCI(ISYT3B0KL))
C
                  CALL GET_T3B0_JK_O(WORK(KT3B0KL3),ISYT3B0KL,
     *                           WORK(KL2TP),ISYML0,
     *                           WORK(KT3BOL2),WORK(KT3BOG2),ISYM0,
     *                           ISYML,L,ISYMK,K,
     *                           WORK(KEND2),LWRK2)
                  CALL GET_T3B0_JK_V(WORK(KT3B0KL3),ISYT3B0KL,
     *                               WORK(KL2TP),ISYML0,
     *                               WORK(KXGADCK),WORK(KXLADCK),
     *                               ISYM0,ISYML,L,ISYMK,K,
     *                               WORK(KEND2),LWRK2)
C
                  CALL GET_T3B0_JK_L1F(WORK(KT3B0KL3),ISYT3B0KL,
     *                            WORK(KL1AM),ISYML0,
     *                            WORK(KXIAJB),ISYM0,
     *                            WORK(KL2TP),ISYML0,
     *                            WORK(KFCKBA),ISYM0,
     *                            ISYML,L,ISYMK,K)

C
C--------------------------------------------------------------
C                 Divide by orbital energy difference and remove 
C                 forbidden elements
C--------------------------------------------------------------
C
                 CALL T3JK_DIA(WORK(KT3B0KL3),ISYT3B0KL,ISYML,L,ISYMK,K,
     *                          WORK(KFOCKD))
                 CALL T3_FORBIDDEN_JK(WORK(KT3B0KL3),ISYMT3B,ISYML,L,
     *                                 ISYMK,K)
                 CALL DSCAL(NMAABCI(ISYT3B0KL),-1.0D0/3.0D0,
     *                      WORK(KT3B0KL3),1)

c               call sum_pt3_jk(work(KT3B0KL3),isyml,l,isymk,k,isyt3b0kl,
c    *                          work(kx3am),6)
C
                  IF (IPRINT .GT. 55) THEN
                    XNORMVAL = DDOT(NMAABCI(ISYT3B0KL),WORK(KT3B0KL3),1,
     *                               WORK(KT3B0KL3),1)
                    WRITE(LUPRI,*)'NORM OF KT3B0KL3 CC3_ADENOCC_CUB ',
     *                              XNORMVAL
                  END IF
C
                  !KT3B0KL is reused here
                  CALL DZERO(WORK(KT3B0KL),NMAABCI(ISTETAZKL))

                  IOPT = 2
                  CALL TETAX_JK_BC_CUB(.TRUE.,.FALSE.,
     *                                 IOPT,WORK(KT30KL),ISYT30KL,
     *                                 WORK(KFOCKRZ),ISYMRZ,
     *                                 WORK(KT3B0KL),ISTETAZKL,
     *                                 WORK(KEND2),LWRK2)

C
C                 KT3B0KL(d- e- b-)_(LMn) = KT3B0KL(d- e- b)_(LMn)
C                                        + KT30KL(deb)_(LMn) * FOCKZ
C
                  CALL TETAX_JK_A(WORK(KT30KL),ISYT30KL,WORK(KFOCKRZ),
     *                                 ISYMRZ,WORK(KT3B0KL),ISTETAZKL,
     *                                 WORK(KEND2),LWRK2)

C
                  CALL W3JK_DIA(WORK(KT3B0KL),ISTETAZKL,ISYML,L,ISYMK,K,
     *                           WORK(KFOCKD),FREQRZ)
C
                  CALL T3_FORBIDDEN_JK(WORK(KT3B0KL),ISYMTETAZ,ISYML,L,
     *                                  ISYMK,K)

C
C                 KTETAXKL = KT3B0KL(d- e- b-)_(LMn) * FOCKU
C

                  CALL DZERO(WORK(KTETAXKL),NMAABCI(ISTETAZUKL))
C
                  IOPT = 2
                  CALL TETAX_JK_BC_CUB(.TRUE.,.FALSE.,
     *                                 IOPT,WORK(KT3B0KL),ISTETAZKL,
     *                                 WORK(KFOCKRU),ISYMRU,
     *                                 WORK(KTETAXKL),ISTETAZUKL,
     *                                 WORK(KEND2),LWRK2)
C
                  CALL TETAX_JK_A(WORK(KT3B0KL),ISTETAZKL,WORK(KFOCKRU),
     *                                 ISYMRU,WORK(KTETAXKL),ISTETAZUKL,
     *                                 WORK(KEND2),LWRK2)
C
C
C                 INCLUDE P(ZU) permutation
C
C

C
C                  KT3B0KL(d- e- b)_(LMn) = KT30KL(deb)_(LMn) * FOCKU
C
                  !KT3B0KL is reused here
                  CALL DZERO(WORK(KT3B0KL),NMAABCI(ISTETAUKL))

                  IOPT = 2
                  CALL TETAX_JK_BC_CUB(.TRUE.,.FALSE.,
     *                                 IOPT,WORK(KT30KL),ISYT30KL,
     *                                 WORK(KFOCKRU),ISYMRU,
     *                                 WORK(KT3B0KL),ISTETAUKL,
     *                                 WORK(KEND2),LWRK2)

C
C                 KT3B0KL(d- e- b-)_(LMn) = KT3B0KL(d- e- b)_(LMn)
C                                        + KT30KL(deb)_(LMn) * FOCKU
C
                  CALL TETAX_JK_A(WORK(KT30KL),ISYT30KL,WORK(KFOCKRU),
     *                                 ISYMRU,WORK(KT3B0KL),ISTETAUKL,
     *                                 WORK(KEND2),LWRK2)

C
                  CALL W3JK_DIA(WORK(KT3B0KL),ISTETAUKL,ISYML,L,ISYMK,K,
     *                           WORK(KFOCKD),FREQRU)
C
                  CALL T3_FORBIDDEN_JK(WORK(KT3B0KL),ISYMTETAU,ISYML,L,
     *                                  ISYMK,K)

C
C                 KTETAXKL = KT3B0KL(d- e- b-)_(LMn) * FOCKZ
C
C
                  IOPT = 2
                  CALL TETAX_JK_BC_CUB(.TRUE.,.FALSE.,
     *                                 IOPT,WORK(KT3B0KL),ISTETAUKL,
     *                                 WORK(KFOCKRZ),ISYMRZ,
     *                                 WORK(KTETAXKL),ISTETAZUKL,
     *                                 WORK(KEND2),LWRK2)
C
                  CALL TETAX_JK_A(WORK(KT3B0KL),ISTETAUKL,WORK(KFOCKRZ),
     *                                 ISYMRZ,WORK(KTETAXKL),ISTETAZUKL,
     *                                 WORK(KEND2),LWRK2)
C
                  CALL W3JK_DIA(WORK(KTETAXKL),ISTETAZUKL,ISYML,L,
     *                          ISYMK,K,WORK(KFOCKD),FREQZU)
C
                  CALL T3_FORBIDDEN_JK(WORK(KTETAXKL),ISYMTETAZU,
     *                                 ISYML,L,ISYMK,K)
 
c                 call sum_pt3_jk(work(KTETAXKL),isyml,l,isymk,k,
c    *                            ISYMTETAZU,
c    *                            work(kx3am),4)



                  !CONTRACTION: 3rd line of Eq. 61
                

                  IOPT = 2
                  CALL ADEN_DAB_LM_CUB(IOPT,DAB,
     *                             WORK(KTETAXKL),ISTETAZUKL,
     *                             WORK(KT3B0KL3),ISYT3B0KL,
     *                             WORK(KEND2),LWRK2)
                  !CONTRACTION: last line of Eq. 62 (1st term)
 
                  CALL ADEN_DIJ_JK(DIJ,WORK(KTETAXKL),ISTETAZUKL,
     *                             WORK(KT3B0KL3),ISYT3B0KL)

                  call dzero(work(KW3BXKL),NMAABCI(ISTETAZKL))
                  call dzero(work(KABCI),NMAABCI(ISTETAUKL))
                  !the real construction of wJK(abci-)

                  CALL TETAX_JK_I(WORK(KT30KL),ISYT30KL,
     *                            WORK(KFOCKRZ),ISYMRZ,
     *                            WORK(KW3BXKL),ISTETAZKL,
     *                            WORK(KEND2),LWRK2)
C
                  CALL WJK_T2(ONE,L,ISYML,K,ISYMK,WORK(KT2TP),ISYM0,
     *                        WORK(KT2TP),
     *                        ISYM0,
     *                        WORK(KFOCKRZ),ISYMRZ,
     *                        WORK(KW3BXKL),ISTETAZKL,
     *                        WORK(KEND2),LWRK2)

c                call sum_pt3_jk(work(KW3BXKL),isyml,l,isymk,k,ISTETAZKL,
c    *                                       work(kx3am),7)
C
                  CALL WJK_GROUND(WORK(KW3BXKL),ISTETAZKL,
     *                            WORK(KT2RZ),ISYMRZ,
     *                            WORK(KGBCDK),ISYM0,
     *                            ISYML,L,ISYMK,K,
     *                            WORK(KEND2),LWRK2)
                  CALL WJK_GROUND_OCC(WORK(KW3BXKL),ISTETAZKL,
     *                                WORK(KT2RZ),ISYMRZ,
     *                                WORK(KT3OG2),ISYM0,
     *                                ISYML,L,ISYMK,K,
     *                                WORK(KEND2),LWRK2)

                  CALL WJK_GROUND(WORK(KW3BXKL),ISTETAZKL,
     *                            WORK(KT2TP),ISYM0,
     *                            WORK(KGBCDKZ),ISYMRZ,
     *                            ISYML,L,ISYMK,K,
     *                            WORK(KEND2),LWRK2)
                  CALL WJK_GROUND_OCC(WORK(KW3BXKL),ISTETAZKL,
     *                                WORK(KT2TP),ISYM0,
     *                                WORK(KT3OG2Z),ISYMRZ,
     *                                ISYML,L,ISYMK,K,
     *                                WORK(KEND2),LWRK2)

c                 call sum_pt3_jk(work(KW3BXKL),isyml,l,isymk,k,ISTETAZKL,
c    *                           work(kx3am),7)

C
                  CALL W3JK_DIA(WORK(KW3BXKL),ISTETAZKL,ISYML,L,ISYMK,K,
     *                           WORK(KFOCKD),FREQRZ)
C
                  CALL T3_FORBIDDEN_JK(WORK(KW3BXKL),ISYMTETAZ,ISYML,L,
     *                                  ISYMK,K)

c                 call sum_pt3_jk(work(KW3BXKL),isyml,l,isymk,k,ISTETAZKL,
c    *                           work(kx3am),7)

      

                 !to include P(ZU) permutation in KABCI

                 CALL TETAX_JK_I(WORK(KT30KL),ISYT30KL,
     *                           WORK(KFOCKRU),ISYMRU,
     *                           WORK(KABCI),ISTETAUKL,
     *                           WORK(KEND2),LWRK2)

                 CALL WJK_T2(ONE,L,ISYML,K,ISYMK,WORK(KT2TP),ISYM0,
     *                       WORK(KT2TP),
     *                       ISYM0,
     *                       WORK(KFOCKRU),ISYMRU,
     *                       WORK(KABCI),ISTETAUKL,
     *                       WORK(KEND2),LWRK2)
C
                 CALL WJK_GROUND(WORK(KABCI),ISTETAUKL,
     *                           WORK(KT2RU),ISYMRU,
     *                           WORK(KGBCDK),ISYM0,
     *                           ISYML,L,ISYMK,K,
     *                           WORK(KEND2),LWRK2)
                 CALL WJK_GROUND_OCC(WORK(KABCI),ISTETAUKL,
     *                               WORK(KT2RU),ISYMRU,
     *                               WORK(KT3OG2),ISYM0,
     *                               ISYML,L,ISYMK,K,
     *                               WORK(KEND2),LWRK2)
                 CALL WJK_GROUND(WORK(KABCI),ISTETAUKL,
     *                           WORK(KT2TP),ISYM0,
     *                           WORK(KGBCDKU),ISYMRU,
     *                           ISYML,L,ISYMK,K,
     *                           WORK(KEND2),LWRK2)
                 CALL WJK_GROUND_OCC(WORK(KABCI),ISTETAUKL,
     *                               WORK(KT2TP),ISYM0,
     *                               WORK(KT3OG2U),ISYMRU,
     *                               ISYML,L,ISYMK,K,
     *                               WORK(KEND2),LWRK2)

C
                  CALL W3JK_DIA(WORK(KABCI),ISTETAUKL,ISYML,L,ISYMK,K,
     *                           WORK(KFOCKD),FREQRU)
C
                  CALL T3_FORBIDDEN_JK(WORK(KABCI),ISYMTETAU,ISYML,L,
     *                                  ISYMK,K)

c                 call sum_pt3_jk(work(KW3BXKL),isyml,l,isymk,k,ISTETAZKL,
c    *                           work(kx3am),7)

                  CALL DZERO(WORK(KTETAXKL),NMAABCI(ISTETAZUKL))
C
                 !1st cont to 57
                  CALL TETAX_JK_A(WORK(KW3BXKL),ISTETAZKL,
     *                                 WORK(KFOCKRU),ISYMRU,
     *                                 WORK(KTETAXKL),ISTETAZUKL,
     *                                 WORK(KEND2),LWRK2)
C
                 !2nd cont to 57
                  CALL TETAX_JK_A(WORK(KT30KL),ISYT30KL,
     *                                 WORK(KFCKUZV),ISYMZU,
     *                                 WORK(KTETAXKL),ISTETAZUKL,
     *                                 WORK(KEND2),LWRK2)
                  !3rd cont to 57
                  CALL DZERO(WORK(KT3B0KL),NMAABCI(ISTETAZKL))
C
                  ! thetaZ(deb-)_(LMn)
                  IOPT = 3
                  CALL TETAX_JK_A(WORK(KT30KL),ISYT30KL,
     *                                 WORK(KFOCKRZ),ISYMRZ,
     *                                 WORK(KT3B0KL),ISTETAZKL,
     *                                 WORK(KEND2),LWRK2)
C
                  CALL W3JK_DIA(WORK(KT3B0KL),ISTETAZKL,ISYML,L,ISYMK,K,
     *                           WORK(KFOCKD),FREQRZ)
C
                  CALL T3_FORBIDDEN_JK(WORK(KT3B0KL),ISYMTETAZ,ISYML,L,
     *                                  ISYMK,K)
C
                  ! thetaZU(deb-)_(LMn-) = thetaZ(deb- )_(LMk) *FOCKU(k,n)
C
                  CALL TETAX_JK_I(WORK(KT3B0KL),ISTETAZKL,
     *                            WORK(KFOCKRU),ISYMRU,
     *                            WORK(KTETAXKL),ISTETAZUKL,
     *                            WORK(KEND2),LWRK2)
C
C      INCLUDE P(ZU) permutation
C
                 !1st cont to 57
                  CALL TETAX_JK_A(WORK(KABCI),ISTETAUKL,
     *                                 WORK(KFOCKRZ),ISYMRZ,
     *                                 WORK(KTETAXKL),ISTETAZUKL,
     *                                 WORK(KEND2),LWRK2)

                 !2nd cont to 57
                  CALL TETAX_JK_A(WORK(KT30KL),ISYT30KL,
     *                                 WORK(KFCKZUV),ISYMZU,
     *                                 WORK(KTETAXKL),ISTETAZUKL,
     *                                 WORK(KEND2),LWRK2)
                  !3rd cont to 57
                  CALL DZERO(WORK(KT3B0KL),NMAABCI(ISTETAUKL))
C
                  ! thetaU(deb-)_(LMn)
                  IOPT = 3
                  CALL TETAX_JK_A(WORK(KT30KL),ISYT30KL,
     *                                 WORK(KFOCKRU),ISYMRU,
     *                                 WORK(KT3B0KL),ISTETAUKL,
     *                                 WORK(KEND2),LWRK2)
C
                  CALL W3JK_DIA(WORK(KT3B0KL),ISTETAUKL,ISYML,L,ISYMK,K,
     *                           WORK(KFOCKD),FREQRU)
C
                  CALL T3_FORBIDDEN_JK(WORK(KT3B0KL),ISYMTETAU,ISYML,L,
     *                                  ISYMK,K)
C
                  ! thetaZU(deb-)_(LMn-) = thetaU(deb- )_(LMk) *FOCKZ(k,n)
C
                  CALL TETAX_JK_I(WORK(KT3B0KL),ISTETAUKL,
     *                            WORK(KFOCKRZ),ISYMRZ,
     *                            WORK(KTETAXKL),ISTETAZUKL,
     *                            WORK(KEND2),LWRK2)
C
                  CALL W3JK_DIA(WORK(KTETAXKL),ISTETAZUKL,ISYML,L,
     *                          ISYMK,K,WORK(KFOCKD),FREQZU)
 
                  CALL T3_FORBIDDEN_JK(WORK(KTETAXKL),ISYMTETAZU,
     *                                 ISYML,L,ISYMK,K)

                  !last line in 61 (term 1)
                  iopt = 3
                  CALL ADEN_DAB_LM_CUB(IOPT,DAB,
     *                             WORK(KTETAXKL),ISTETAZUKL,
     *                             WORK(KT3B0KL3),ISYT3B0KL,
     *                             WORK(KEND2),LWRK2)

c--------
                  !construct theta for last line in 61 (term 2 & 3)
                  !KT3B0KL is reused
                  CALL DZERO(WORK(KT3B0KL),NMAABCI(ISTETAZKL))
C
                  ! thetaZ(de- b)_(LMn)
                  IOPT = 3
                  CALL TETAX_JK_BC_CUB(.TRUE.,.FALSE.,
     *                                 IOPT,WORK(KT30KL),ISYT30KL,
     *                                 WORK(KFOCKRZ),ISYMRZ,
     *                                 WORK(KT3B0KL),ISTETAZKL,
     *                                 WORK(KEND2),LWRK2)
C
                  CALL W3JK_DIA(WORK(KT3B0KL),ISTETAZKL,ISYML,L,ISYMK,K,
     *                           WORK(KFOCKD),FREQRZ)
C
                  CALL T3_FORBIDDEN_JK(WORK(KT3B0KL),ISYMTETAZ,ISYML,L,
     *                                  ISYMK,K)
C
                  ! thetaZU(de- b)_(LMn-) = thetaZ(de- b)_(LMk) *FOCKU(k,n)
                  CALL DZERO(WORK(KTETAXKL),NMAABCI(ISTETAZUKL))
C
                  CALL TETAX_JK_I(WORK(KT3B0KL),ISTETAZKL,
     *                            WORK(KFOCKRU),ISYMRU,
     *                            WORK(KTETAXKL),ISTETAZUKL,
     *                            WORK(KEND2),LWRK2)
C
C                INCLUDE P(ZU) permutation now
C

                  CALL DZERO(WORK(KT3B0KL),NMAABCI(ISTETAUKL))
C
                  ! thetaU(de- b)_(LMn)
                  IOPT = 3
                  CALL TETAX_JK_BC_CUB(.TRUE.,.FALSE.,
     *                                 IOPT,WORK(KT30KL),ISYT30KL,
     *                                 WORK(KFOCKRU),ISYMRU,
     *                                 WORK(KT3B0KL),ISTETAUKL,
     *                                 WORK(KEND2),LWRK2)
C
                  CALL W3JK_DIA(WORK(KT3B0KL),ISTETAUKL,ISYML,L,ISYMK,K,
     *                           WORK(KFOCKD),FREQRU)
C
                  CALL T3_FORBIDDEN_JK(WORK(KT3B0KL),ISYMTETAU,ISYML,L,
     *                                  ISYMK,K)
C
                  ! thetaZU(de- b)_(LMn-) = thetaU(de- b)_(LMk) *FOCKZ(k,n)
                  CALL TETAX_JK_I(WORK(KT3B0KL),ISTETAUKL,
     *                            WORK(KFOCKRZ),ISYMRZ,
     *                            WORK(KTETAXKL),ISTETAZUKL,
     *                            WORK(KEND2),LWRK2)
C

                  IOPT = 3
                  CALL TETAX_JK_BC_CUB(.TRUE.,.FALSE.,
     *                                 IOPT,WORK(KW3BXKL),ISTETAZKL,
     *                                 WORK(KFOCKRU),ISYMRU,
     *                                 WORK(KTETAXKL),ISTETAZUKL,
     *                                 WORK(KEND2),LWRK2)
C
C      INCLUDE P(ZU) permutation
C
                  IOPT = 3
                  CALL TETAX_JK_BC_CUB(.TRUE.,.FALSE.,
     *                                 IOPT,WORK(KABCI),ISTETAUKL,
     *                                 WORK(KFOCKRZ),ISYMRZ,
     *                                 WORK(KTETAXKL),ISTETAZUKL,
     *                                 WORK(KEND2),LWRK2)

C
                  CALL W3JK_DIA(WORK(KTETAXKL),ISTETAZUKL,ISYML,L,
     *                          ISYMK,K,WORK(KFOCKD),FREQZU)
 
                  CALL T3_FORBIDDEN_JK(WORK(KTETAXKL),ISYMTETAZU,
     *                                 ISYML,L,ISYMK,K)


                  !last line in 61 (term 3)
                  iopt = 0
                  CALL ADEN_DAB_LM_CUB(IOPT,DAB,
     *                             WORK(KTETAXKL),ISTETAZUKL,
     *                             WORK(KT3B0KL3),ISYT3B0KL,
     *                             WORK(KEND2),LWRK2)
                  !last line in 61 (term 2) 
                  iopt = 1
                  CALL ADEN_DAB_LM_CUB_t0(IOPT,DAB,
     *                             WORK(KTETAXKL),ISTETAZUKL,
     *                             WORK(KT3B0KL3),ISYT3B0KL,
     *                             WORK(KEND2),LWRK2)

                  !intermmediates for last line of Eq. 62 (2nd term)
                  CALL DZERO(WORK(KT3B0KL),NMAABCI(ISTETAZKL))
C
                  ! thetaZ(de- b)_(LMn)
                  IOPT = 3
                  CALL TETAX_JK_BC_CUB(.TRUE.,.FALSE.,
     *                                 IOPT,WORK(KT30KL),ISYT30KL,
     *                                 WORK(KFOCKRZ),ISYMRZ,
     *                                 WORK(KT3B0KL),ISTETAZKL,
     *                                 WORK(KEND2),LWRK2)
C
                  CALL W3JK_DIA(WORK(KT3B0KL),ISTETAZKL,ISYML,L,ISYMK,K,
     *                           WORK(KFOCKD),FREQRZ)
C
                  CALL T3_FORBIDDEN_JK(WORK(KT3B0KL),ISYMTETAZ,ISYML,L,
     *                                  ISYMK,K)
C
                  CALL DZERO(WORK(KTETAXKL),NMAABCI(ISTETAZUKL))

                  ! thetaZU(de- b-)_(LMn)
                  CALL TETAX_JK_A(WORK(KT3B0KL),ISTETAZKL,
     *                                 WORK(KFOCKRU),ISYMRU,
     *                                 WORK(KTETAXKL),ISTETAZUKL,
     *                                 WORK(KEND2),LWRK2)
C
                  CALL DZERO(WORK(KT3B0KL),NMAABCI(ISTETAZKL))
C
                  ! thetaZ(deb- )_(LMn)
                  CALL TETAX_JK_A(WORK(KT30KL),ISYT30KL,
     *                                 WORK(KFOCKRZ),ISYMRZ,
     *                                 WORK(KT3B0KL),ISTETAZKL,
     *                                 WORK(KEND2),LWRK2)
C
                  CALL W3JK_DIA(WORK(KT3B0KL),ISTETAZKL,ISYML,L,ISYMK,K,
     *                           WORK(KFOCKD),FREQRZ)
C
                  CALL T3_FORBIDDEN_JK(WORK(KT3B0KL),ISYMTETAZ,ISYML,L,
     *                                  ISYMK,K)
C
                  ! thetaZU(de- b-)_(LMn)
                  CALL TETAX_JK_BC_CUB(.TRUE.,.FALSE.,
     *                                 IOPT,WORK(KT3B0KL),ISTETAZKL,
     *                                 WORK(KFOCKRU),ISYMRU,
     *                                 WORK(KTETAXKL),ISTETAZUKL,
     *                                 WORK(KEND2),LWRK2)


C
C                 INCLUDE P(ZU) permutation
C

                  CALL DZERO(WORK(KT3B0KL),NMAABCI(ISTETAUKL))
C
                  ! thetaU(de- b)_(LMn)
                  IOPT = 3
                  CALL TETAX_JK_BC_CUB(.TRUE.,.FALSE.,
     *                                 IOPT,WORK(KT30KL),ISYT30KL,
     *                                 WORK(KFOCKRU),ISYMRU,
     *                                 WORK(KT3B0KL),ISTETAUKL,
     *                                 WORK(KEND2),LWRK2)
C
                  CALL W3JK_DIA(WORK(KT3B0KL),ISTETAUKL,ISYML,L,ISYMK,K,
     *                           WORK(KFOCKD),FREQRU)
C
                  CALL T3_FORBIDDEN_JK(WORK(KT3B0KL),ISYMTETAU,ISYML,L,
     *                                  ISYMK,K)
C
                  ! thetaZU(de- b-)_(LMn)
                  CALL TETAX_JK_A(WORK(KT3B0KL),ISTETAUKL,
     *                                 WORK(KFOCKRZ),ISYMRZ,
     *                                 WORK(KTETAXKL),ISTETAZUKL,
     *                                 WORK(KEND2),LWRK2)
C
                  CALL DZERO(WORK(KT3B0KL),NMAABCI(ISTETAUKL))
C
                  ! thetaZ(deb- )_(LMn)
                  CALL TETAX_JK_A(WORK(KT30KL),ISYT30KL,
     *                                 WORK(KFOCKRU),ISYMRU,
     *                                 WORK(KT3B0KL),ISTETAUKL,
     *                                 WORK(KEND2),LWRK2)
C
                  CALL W3JK_DIA(WORK(KT3B0KL),ISTETAUKL,ISYML,L,ISYMK,K,
     *                           WORK(KFOCKD),FREQRU)
C
                  CALL T3_FORBIDDEN_JK(WORK(KT3B0KL),ISYMTETAU,ISYML,L,
     *                                  ISYMK,K)
C
                  ! thetaZU(de- b-)_(LMn)
                  CALL TETAX_JK_BC_CUB(.TRUE.,.FALSE.,
     *                                 IOPT,WORK(KT3B0KL),ISTETAUKL,
     *                                 WORK(KFOCKRZ),ISYMRZ,
     *                                 WORK(KTETAXKL),ISTETAZUKL,
     *                                 WORK(KEND2),LWRK2)
C
                  CALL W3JK_DIA(WORK(KTETAXKL),ISTETAZUKL,ISYML,L,
     *                          ISYMK,K,WORK(KFOCKD),FREQZU)

                  CALL T3_FORBIDDEN_JK(WORK(KTETAXKL),ISYMTETAZU,
     *                                 ISYML,L,ISYMK,K)



                  !CONTRACTION: last line of Eq. 62 (2nd term)
C
                  CALL ADEN_DIJ_JK_t0(DIJ,WORK(KTETAXKL),ISTETAZUKL,
     *                             WORK(KT3B0KL3),ISYT3B0KL)

C
                  IF (IPRINT .GT. 55) THEN
                    XNORMVAL = DDOT(NMATAB(ISYDEN),DAB,1,DAB,1)
                    WRITE(LUPRI,*)'NORM OF DAB AFTER ADEN_DAB_LM ',
     *                             XNORMVAL
                  END IF
C
               ENDDO   ! L
            ENDDO      ! ISYML
         ENDDO       ! K
      ENDDO          ! ISYMK 
C
c      write(lupri,*) 'T30KL in CC3_ADENOCC_CUB_t0'
c      call print_pt3(work(kx3am),isym0,4)
C
C--------------------------------
C     Close files for "response"
C--------------------------------
C
      CALL WCLOSE2(LU3SRTR,FN3SRTR,'DELETE')
      CALL WCLOSE2(LUCKJDRZ,FNCKJDRZ,'DELETE')
      CALL WCLOSE2(LUDELDRZ,FNDELDRZ,'DELETE')
      CALL WCLOSE2(LUDKBCRZ,FNDKBCRZ,'DELETE')
C
      CALL WCLOSE2(LUCKJDRU,FNCKJDRU,'DELETE')
      CALL WCLOSE2(LUDELDRU,FNDELDRU,'DELETE')
      CALL WCLOSE2(LUDKBCRU,FNDKBCRU,'DELETE')
C
C
C-------------
C     End
C-------------
C

      CALL QEXIT('CC3AOCB')
C
      RETURN
      END
C  /* Deck cc3_adenvir_cub_t0 */
      SUBROUTINE CC3_ADENVIR_CUB_T0(DIJ,DAB,ISYDEN,
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
*    Simplification of CC3_ADENVIR_CUB.
*    tbarX are replaced byt tbar0, so some terms are dropped.
*
*    Written by Filip Pawlowski, Fall 2003, Aarhus
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
#include <ccr2rsp.h>

C
      CHARACTER*10 MODEL
C
      LOGICAL DO_YMMAT
C
      INTEGER ISYM0
      PARAMETER(ISYM0 = 1)
      CHARACTER LISTL0*3, LISTL*3,LISTR*3,LISTL1R*3,LABELL1*8,LABELRZ*8
      CHARACTER LABELRU*8
      CHARACTER LISTRZ*3,LISTRU*3
      CHARACTER*(*) FNTOC, FN3VI, FNDKBC3, FN3FOPX, FN3FOP2X 
      CHARACTER*(*) FNDKBC,FNDELD,FN3VI2,FN3FOP,FN3FOP2,FNCKJD
C
      CHARACTER*11 FNWBMATX
      PARAMETER(FNWBMATX = 'CC3_W3_TMPX')
      INTEGER LUWBMATX
C
      CHARACTER*12 FNT2Y
      PARAMETER ( FNT2Y = 'CC3_R2TP_TMP' )
      INTEGER LUT2Y
C
      CHARACTER*10 FNT3, FNWBMAT
      CHARACTER*14 FNTHETA,FNWZU
      PARAMETER(FNT3 = 'CC3_T3_TMP', FNWBMAT = 'CC3_W3_TMP',  
     *          FNTHETA = 'CC3_THETA3_TMP',FNWZU = 'CC3_WZU____TMP')
C
      CHARACTER*14 FN3SRTR, FNCKJDRZ, FNDELDRZ, FNDKBCRZ
      PARAMETER(FN3SRTR = 'CCSDT_FBMAT1_Z',FNCKJDRZ = 'CCSDT_FBMAT2_Z',
     *          FNDELDRZ = 'CCSDT_FBMAT3_Z',FNDKBCRZ = 'CCSDT_FBMAT4_Z')
      INTEGER LU3SRTR, LUCKJDRZ, LUDELDRZ, LUDKBCRZ
C
      CHARACTER*14 FNCKJDRU, FNDELDRU, FNDKBCRU
      PARAMETER(FNCKJDRU = 'CCSDT_FBMAT2_U',
     *          FNDELDRU = 'CCSDT_FBMAT3_U',FNDKBCRU = 'CCSDT_FBMAT4_U')
      INTEGER LUCKJDRU, LUDELDRU, LUDKBCRU
C
      ![[H,T1Z],T1U]
      CHARACTER*14 FNCKJDRZU, FNDELDRZU, FNDKBCRZU
      PARAMETER(FNCKJDRZU ='CCSDT_FBMAT2ZU',
     *          FNDELDRZU ='CCSDT_FBMAT3ZU',FNDKBCRZU ='CCSDT_FBMAT4ZU')
      INTEGER LUCKJDRZU, LUDELDRZU, LUDKBCRZU
C
      ![H,T1ZU]
      CHARACTER*14 FNCKJDR2, FNDELDR2, FNDKBCR2
      PARAMETER(FNCKJDR2 = 'CCSDT_FBMAT2R2',
     *          FNDELDR2 = 'CCSDT_FBMAT3R2',FNDKBCR2 = 'CCSDT_FBMAT4R2')
      INTEGER LUCKJDR2, LUDELDR2, LUDKBCR2
C
C
      INTEGER LUTOC, LU3VI, LUDKBC3, LU3FOPX, LU3FOP2X
      INTEGER LUDKBC,LUDELD,LU3VI2,LU3FOP,LU3FOP2,LUCKJD
      INTEGER LUT3,LUWBMAT,LUTHETA,LUWZU
C
      LOGICAL   LOCDBG,LORXL1
      LOGICAL   LORXRZ,LORXRU
      PARAMETER (LOCDBG = .FALSE.)
C
      INTEGER  AIBJCK_PERM
      LOGICAL QUADR
      LOGICAL T2XNET2Y
      LOGICAL T2XNET2Z,NOVIRT
C
      CHARACTER CDUMMY*1 
      PARAMETER (CDUMMY = ' ')

      INTEGER   ISYDEN,IDLSTL,IDLSTR,LWORK
C
      INTEGER IDLSTL0,IDLSTL1R
      INTEGER ISYML1,ISYML1R,ISYMRZ,ISYMRU
      INTEGER ISINT1,ISINT2
      INTEGER KLAMP0,KLAMH0,KFOCKD,KFOCK0CK,KT2TP,KL1AM,KL2TP
      INTEGER KEND0,LWRK0
      INTEGER KL1L1,KL2L1,KT1RZ,KT2RZ,KFOCK0,KFOCKL1,KFOCKRZ
      INTEGER KEND1,LWRK1
      INTEGER IOPT
      INTEGER ISINT1RZ,ISINT2RZ,ISINT2L1R,ISYFCKL1R
      INTEGER KXIAJB,KT3BOG1,KT3BOL1,KT3BOG2,KT3BOL2,KT3OG1,KT3OG2
      INTEGER KLAMPL1R,KLAMHL1R,KW3ZOGZ1,KFOCKL1RCK,KW3BXOG1
      INTEGER KW3BXOL1,KW3BXOGX1,KW3BXOLX1,KT1L1R,KT2L1R
      INTEGER KEND2,LWRK2
      INTEGER LENGTH
      INTEGER ISINT1L1R
      INTEGER ISYMD,ISYCKBD0,ISYCKBDL1R,ISYCKBDR1Z
      INTEGER KT3VDG1,KT3VDG2,KT3VDG3,KT3BVDL1,KT3BVDL2,KT3BVDL3
      INTEGER KEND3,LWRK3
      INTEGER KT3BVDG1,KT3BVDG2,KT3BVDG3,KW3BXVDG1,KW3BXVDG2
      INTEGER KW3BXVDL1,KW3BXVDL2,KW3BXVDGX1,KW3BXVDGX2,KW3BXVDLX1
      INTEGER KW3BXVDLX2,KW3ZVDGZ1,KINTVI,KTRVI6
      INTEGER KEND4,LWRK4
      INTEGER IOFF
      INTEGER ISYMB,ISYALJB0,ISYALJD0,ISYALJBL1,ISYALJDL1,ISYMBD
      INTEGER ISCKIJ,ISWBMAT,ISWMATZ,ISYCKD,ISYCKDBR1Z
      INTEGER KSMAT2,KUMAT2,KDIAG,KDIAGWB,KDIAGWZ,KINDSQ,KINDSQWB
      INTEGER KINDSQWZ,KINDEX,KINDEX2,KINDEXBL1,KINDEXDL1,KTMAT
      INTEGER KT3MAT,KW3MATZ,KWTEMP,KS3MAT,KU3MAT,KS3MAT3
      INTEGER KU3MAT3,KT3VBG1,KT3VBG2,KT3VBG3,KT3BVBG1,KT3BVBG2
      INTEGER KT3BVBG3,KSMAT4,KUMAT4,KT3BVBL1,KT3BVBL2,KT3BVBL3
      INTEGER KW3ZVDGZ2
      INTEGER KEND5,LWRK5
      INTEGER LENSQ,LENSQWB,LENSQWZ
      INTEGER ISYML,ISYMDL,ISAIBJ,ISYMJ,ISYMBJ,ISYMAI,ISYAIL
      INTEGER KOFF1,NBJ,IADR
      INTEGER KT3VBGZ3
      INTEGER IDLSTZU,IDLSTRZ,IDLSTRU
      INTEGER KT3VDGZ3,KFOCKRU,KWMATZU,KFCKUZO,KLAMDPZ,KLAMDHZ,KINDSQWZU
      INTEGER LENSQWZU,KDIAGWZU,ISYMZU,MAXX1,MAXX2,ISWMATZU
      INTEGER KW3MATU,ISWMATU,KINDSQWU,LENSQWU,KDIAGWU,KT2RU,KT1RU
      INTEGER KW3UVDGU1,KW3UVDGU2,KT3VBGU3,KT3VDGU3,KW3UOGU1
      INTEGER ISINT1RU,ISINT2RU,ISYCKBDR1U,MAXX3,ISYCKDBR1U
      INTEGER KFCKZUO,KLAMDPU,KLAMDHU
      INTEGER KWMATZUD
      INTEGER ISINT1RZU,ISINT2RZU,KW3ZUOGZU1,ISYCKBDR1ZU,KW3ZUVDGZU1
      INTEGER MAXX4,ISYCKDBR1ZU,KW3ZUVDGZU2,KT3VBGZU3,KT3VDGZU3
      INTEGER KT1ZU,KT2ZU
      INTEGER KWZUVDGR21,KWZUVDGR22,KWZUVBGR23,KWZUVDGR23,KWZUOGR21
      INTEGER ISINT1R2,ISINT2R2
      INTEGER KFCKZUV,KFCKUZV
      INTEGER KR2TP
C
      INTEGER IR1TAMP
C
      integer kx3am
C
      INTEGER FKW3BXVDG1,FKW3BXVDG2,FKW3BXVDL1,FKW3BXVDL2
      INTEGER FKW3BXVDGX1,FKW3BXVDGX2,FKW3BXVDLX1,FKW3BXVDLX2

#if defined (SYS_CRAY)
      REAL      FREQL1,FREQL1R,FREQRZ,FREQRU,FREQZU
      REAL      WORK(LWORK)
      REAL      XNORMVAL
      REAL      DAB(*),DIJ(*)
      REAL      YMMAT(*)
      REAL      DDOT,HALF,ONE,TWO
#else
      DOUBLE PRECISION      FREQL1,FREQL1R,FREQRZ,FREQRU,FREQZU
      DOUBLE PRECISION      WORK(LWORK)
      DOUBLE PRECISION      XNORMVAL
      DOUBLE PRECISION      DAB(*),DIJ(*)
      DOUBLE PRECISION      YMMAT(*)
      DOUBLE PRECISION      DDOT,HALF,ONE,TWO
#endif
C
      PARAMETER(HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      CALL QENTER('CC3DENVCB')
C
C--------------------------------
C     Open temporary files
C--------------------------------
C
      IF (DO_YMMAT) THEN
         LUWBMATX = -1
         CALL WOPEN2(LUWBMATX,FNWBMATX,64,0)
      END IF
C
      LU3SRTR   = -1
      LUCKJDRZ  = -1
      LUDELDRZ  = -1
      LUDKBCRZ  = -1
C
      CALL WOPEN2(LU3SRTR,FN3SRTR,64,0)
      CALL WOPEN2(LUCKJDRZ,FNCKJDRZ,64,0)
      CALL WOPEN2(LUDELDRZ,FNDELDRZ,64,0)
      CALL WOPEN2(LUDKBCRZ,FNDKBCRZ,64,0)
C
      LUCKJDRU  = -1
      LUDELDRU  = -1
      LUDKBCRU  = -1
C
      CALL WOPEN2(LUCKJDRU,FNCKJDRU,64,0)
      CALL WOPEN2(LUDELDRU,FNDELDRU,64,0)
      CALL WOPEN2(LUDKBCRU,FNDKBCRU,64,0)
C
      ![[H,T1Z],T1U]
      LUCKJDRZU = -1
      LUDELDRZU = -1
      LUDKBCRZU = -1
C
      CALL WOPEN2(LUCKJDRZU,FNCKJDRZU,64,0)
      CALL WOPEN2(LUDELDRZU,FNDELDRZU,64,0)
      CALL WOPEN2(LUDKBCRZU,FNDKBCRZU,64,0)
C
      ![H,T1ZU]
      LUCKJDR2  = -1
      LUDELDR2  = -1
      LUDKBCR2  = -1
C
      CALL WOPEN2(LUCKJDR2,FNCKJDR2,64,0)
      CALL WOPEN2(LUDELDR2,FNDELDR2,64,0)
      CALL WOPEN2(LUDKBCR2,FNDKBCR2,64,0)
C
C------------------------------------------------------------
C     some initializations:
C------------------------------------------------------------
C
      LISTL0 = 'L0 '
      IDLSTL0 = 0 

      IF (LISTL(1:3).EQ.'L1 ') THEN

         ! get symmetry, frequency and integral label from common blocks
         ! defined in ccl1rsp.h
         ISYML1  = ISYLRZ(IDLSTL)
         FREQL1  = FRQLRZ(IDLSTL)
         LABELL1 = LRZLBL(IDLSTL)
         LORXL1  = LORXLRZ(IDLSTL)
c

         IF (LORXL1) CALL QUIT('NO ORBITAL RELAX. IN CC3_ADENVIR_CUB')

        LISTL1R  = 'R1 '
        IDLSTL1R = IR1TAMP(LABELL1,LORXL1,FREQL1,ISYML1)
        ! get symmetry and frequency from common blocks
        ! defined in ccl1rsp.h
        ISYML1R  = ISYLRT(IDLSTL1R)
        FREQL1R  = FRQLRT(IDLSTL1R)
C
        IF (ISYML1 .NE. ISYML1R) THEN
           WRITE(LUPRI,*)'ISYML1: ', ISYML1
           WRITE(LUPRI,*)'ISYML1R: ', ISYML1R
           CALL QUIT('Symmetry mismatch in CC3_ADENVIR_CUB')
        END IF
C
        IF (FREQL1R .NE. FREQL1) THEN
           WRITE(LUPRI,*)'FREQL1R: ', FREQL1R
           WRITE(LUPRI,*)'FREQL1: ', FREQL1
           CALL QUIT('Frequency mismatch in CC3_ADENVIR_CUB')
        END IF
C
      ELSE
         WRITE(LUPRI,*)'LISTL ',LISTL
         CALL QUIT('Unknown left list in CC3_ADENVIR_CUB_T0')
      END IF

      IF (LISTR(1:3).EQ.'R2 ') THEN
         IDLSTZU = IDLSTR
         ! get symmetry, frequency and integral label for right list 
         ! from common blocks defined in ccr1rsp.h
         LISTRZ  = 'R1 '
         LABELRZ = LBLR2T(IDLSTZU,1)
         ISYMRZ  = ISYR2T(IDLSTZU,1)
         FREQRZ  = FRQR2T(IDLSTZU,1)
         LORXRZ  = LORXR2T(IDLSTZU,1)
         IDLSTRZ = IR1TAMP(LABELRZ,LORXRZ,FREQRZ,ISYMRZ)

         LISTRU  = 'R1 '
         LABELRU = LBLR2T(IDLSTZU,2)
         ISYMRU  = ISYR2T(IDLSTZU,2)
         FREQRU  = FRQR2T(IDLSTZU,2)
         LORXRU  = LORXR2T(IDLSTZU,2)
         IDLSTRU = IR1TAMP(LABELRU,LORXRU,FREQRU,ISYMRU)
 
C
         IF (LORXRZ.OR.LORXRU) THEN
          CALL QUIT('Orbital relaxation not allowed in CC3_ADENVIR_CUB')
         END IF
C
      ELSE
         WRITE(LUPRI,*)'LISTR = ',LISTR(1:3)
         WRITE(LUPRI,*)'CC3_ADENVIR_CUB is designed for LISTR = R2'
         CALL QUIT('Unknown right list in CC3_ADENVIR_CUB')
      END IF

      FREQZU = FREQRZ + FREQRU
C
C-------------------------------------------------------
C     initial allocations, orbital energy, fock matrix and T2 and L2 :
C-------------------------------------------------------
C
C     Symmetry of integrals in contraction:
C
      ISINT1 = ISYM0
      ISINT2 = ISYM0
      ISYMZU = MULD2H(ISYMRZ,ISYMRU)
C
      KLAMP0 = 1
      KLAMH0  = KLAMP0  + NLAMDT
      KFOCKD  = KLAMH0  + NLAMDT
      KFOCK0CK  = KFOCKD  + NORBTS
      KT2TP   = KFOCK0CK  + NT1AMX 
      KL1AM   = KT2TP   + NT2SQ(ISYM0)
      KL2TP   = KL1AM   + NT1AM(ISYM0)
      KEND0   = KL2TP   + NT2SQ(ISYM0)
      LWRK0   = LWORK   - KEND0
C
      KL1L1   = KEND0
      KL2L1   = KL1L1   + NT1AM(ISYML1)
      KT1RZ   = KL2L1   + NT2SQ(ISYML1)
      KT2RZ   = KT1RZ   + NT1AM(ISYMRZ)
      KFOCK0  = KT2RZ   + NT2SQ(ISYMRZ)
      KFOCKL1  = KFOCK0    + N2BST(ISYM0)
      KFOCKRZ   = KFOCKL1    + N2BST(ISYML1)
      KEND1    = KFOCKRZ + N2BST(ISYMRZ)
      LWRK1    = LWORK - KEND1
C
      KT2RU   = KEND1
      KT1RU   = KT2RU + NT2SQ(ISYMRU)
      KEND1   = KT1RU + NT1AM(ISYMRU)
      LWRK1    = LWORK - KEND1
C
      KT2ZU = KEND1
      KEND1 = KT2ZU + NT2SQ(ISYMZU)
      LWRK1    = LWORK - KEND1
C
      KT1ZU = KEND1
      KEND1 = KT1ZU + NT1AM(ISYMZU)
      LWRK1    = LWORK - KEND1
C
      KFOCKRU = KEND1
      KEND1   = KFOCKRU + N2BST(ISYMRU)
      LWRK1    = LWORK - KEND1
C
      KFCKUZO  = KEND1
      KFCKZUO  = KFCKUZO + N2BST(ISYMZU)
      KFCKZUV  = KFCKZUO + N2BST(ISYMZU)
      KFCKUZV  = KFCKZUV + N2BST(ISYMZU)
      KEND1   = KFCKUZV + N2BST(ISYMZU)
C
      KLAMDPZ = KEND1
      KLAMDHZ = KLAMDPZ + NLAMDT
      KLAMDPU = KLAMDHZ + NLAMDT
      KLAMDHU = KLAMDPU + NLAMDT
      KEND1   = KLAMDHU + NLAMDT 
      LWRK1   = LWORK   - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_ADENVIR_CUB (00) ')
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
      ! FZ
      IF (LISTRZ(1:3).EQ.'R1 ') THEN
         CALL GET_FOCKX(WORK(KFOCKRZ),LABELRZ,IDLSTRZ,ISYMRZ,
     *                  WORK(KLAMP0),ISYM0,
     *                  WORK(KLAMH0),ISYM0,WORK(KEND1),LWRK1)
      END IF

      ! FU
      IF (LISTRU(1:3).EQ.'R1 ') THEN
         CALL GET_FOCKX(WORK(KFOCKRU),LABELRU,IDLSTRU,ISYMRU,
     *                  WORK(KLAMP0),ISYM0,
     *                  WORK(KLAMH0),ISYM0,WORK(KEND1),LWRK1)
      END IF
C
C------------------------------------------
C     Calculate the [U,T1^Z] matrix
C     Recall that we only need the occ-occ and vir-vir block.
C------------------------------------------
C
      CALL GET_LAMBDAX(WORK(KLAMDPZ),WORK(KLAMDHZ),LISTRZ,IDLSTRZ,
     *                 ISYMRZ,WORK(KLAMP0),WORK(KLAMH0),WORK(KEND1),
     *                 LWRK1)
      ! get vir-vir block U_(c-,d)
      CALL GET_FOCKX(WORK(KFCKUZV),LABELRU,IDLSTRU,ISYMRU,WORK(KLAMDPZ),
     *                  ISYMRZ,WORK(KLAMH0),ISYM0,WORK(KEND1),LWRK1)
      ! get occ-occ block U_(l,k-)
      CALL GET_FOCKX(WORK(KFCKUZO),LABELRU,IDLSTRU,ISYMRU,WORK(KLAMP0),
     *                  ISYM0,WORK(KLAMDHZ),ISYMRZ,WORK(KEND1),LWRK1)
C
C------------------------------------------
C     Calculate the [Z,T1^U] matrix
C     Recall that we only need the occ-occ and vir-vir block.
C------------------------------------------
C
      CALL GET_LAMBDAX(WORK(KLAMDPU),WORK(KLAMDHU),LISTRU,IDLSTRU,
     *                 ISYMRU,WORK(KLAMP0),WORK(KLAMH0),WORK(KEND1),
     *                 LWRK1)
      ! get vir-vir block Z_(c-,d)
      CALL GET_FOCKX(WORK(KFCKZUV),LABELRZ,IDLSTRZ,ISYMRZ,WORK(KLAMDPU),
     *                  ISYMRU,WORK(KLAMH0),ISYM0,WORK(KEND1),LWRK1)
      ! get occ-occ block Z_(l,k-)
      CALL GET_FOCKX(WORK(KFCKZUO),LABELRZ,IDLSTRZ,ISYMRZ,WORK(KLAMP0),
     *                  ISYM0,WORK(KLAMDHU),ISYMRU,WORK(KEND1),LWRK1)

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
     *               ISYM0,WORK(KEND1),LWRK1)
C
      IF (LOCDBG) WRITE(LUPRI,*) 'Norm of L2TP ',
     *    DDOT(NT2SQ(ISYM0),WORK(KL2TP),1,WORK(KL2TP),1)

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
C     Read T1Z and T2Z amplitudes 
C-------------------------------------
C
      IOPT  = 3
      CALL GET_T1_T2(IOPT,.TRUE.,WORK(KT1RZ),WORK(KT2RZ),LISTRZ,
     *               IDLSTRZ,ISYMRZ,WORK(KEND1),LWRK1)
C
C-------------------------------------
C     Read T1U and T2U amplitudes 
C-------------------------------------
C
      IOPT  = 3
      CALL GET_T1_T2(IOPT,.TRUE.,WORK(KT1RU),WORK(KT2RU),LISTRU,
     *               IDLSTRU,ISYMRU,WORK(KEND1),LWRK1)
C
C-------------------------------------------------------
C     Read in T1^ZU and T2^ZU   !second-order amplitudes
C-------------------------------------------------------
C
      IOPT = 3
      CALL GET_T1_T2(IOPT,.TRUE.,WORK(KT1ZU),WORK(KT2ZU),LISTR,IDLSTR,
     *               ISYMZU,WORK(KEND1),LWRK1)
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
          KEND1 = KR2TP + NT2SQ(ISYMZU)
          LWRK1 = LWORK - KEND1
          IF (LWRK1.LT.NT2AM(ISYMZU))
     &     CALL QUIT('Out of memory in CC3_ADENVIR_CUB_T0 (special T2)')

          IOPT = 2
          CALL CC_RDRSP(LISTR,IDLSTR,ISYMZU,IOPT,MODEL,
     &                  DUMMY,WORK(KEND1))
          CALL CCLR_DIASCL(WORK(KEND1),TWO,ISYMZU)
          CALL CC_T2SQ(WORK(KEND1),WORK(KR2TP),ISYMZU)

          CALL PUTWA2(LUT2Y,FNT2Y,WORK(KR2TP),1,NT2SQ(ISYMZU))
      END IF

C
C----------------------------------------
C     Integrals [H,T1Z] where Z is LISTRZ
C----------------------------------------
C
      ISINT1RZ = MULD2H(ISINT1,ISYMRZ)
      ISINT2RZ = MULD2H(ISINT2,ISYMRZ)
C
      CALL CC3_BARINT(WORK(KT1RZ),ISYMRZ,WORK(KLAMP0),
     *                WORK(KLAMH0),WORK(KEND1),LWRK1,
     *                LU3SRTR,FN3SRTR,LUCKJDRZ,FNCKJDRZ)
C
      CALL CC3_SORT1(WORK(KEND1),LWRK1,2,ISINT1RZ,LU3SRTR,FN3SRTR,
     *               LUDELDRZ,FNDELDRZ,IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *               IDUMMY,CDUMMY)
C
      CALL CC3_SINT(WORK(KLAMH0),WORK(KEND1),LWRK1,ISINT1RZ,
     *              LUDELDRZ,FNDELDRZ,LUDKBCRZ,FNDKBCRZ)

C
C----------------------------------------
C     Integrals [H,T1U] where U is LISTRU
C----------------------------------------
C
      ISINT1RU = MULD2H(ISINT1,ISYMRU)
      ISINT2RU = MULD2H(ISINT2,ISYMRU)
C
      CALL CC3_BARINT(WORK(KT1RU),ISYMRU,WORK(KLAMP0),
     *                WORK(KLAMH0),WORK(KEND1),LWRK1,
     *                LU3SRTR,FN3SRTR,LUCKJDRU,FNCKJDRU)
C
      CALL CC3_SORT1(WORK(KEND1),LWRK1,2,ISINT1RU,LU3SRTR,FN3SRTR,
     *               LUDELDRU,FNDELDRU,IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *               IDUMMY,CDUMMY)
C
      CALL CC3_SINT(WORK(KLAMH0),WORK(KEND1),LWRK1,ISINT1RU,
     *              LUDELDRU,FNDELDRU,LUDKBCRU,FNDKBCRU)

C
C------------------------------------------------------
C     Calculate the (ck|de)-{Z,U}tilde and (ck|lm)-{Z,U}tilde
C     (double one-index transformed with first-order amplitudes)
C------------------------------------------------------
C
      ISINT1RZU = MULD2H(ISINT1,ISYMZU)
      ISINT2RZU = MULD2H(ISINT2,ISYMZU)

      CALL CC3_3BARINT(ISYMRZ,LISTRZ,IDLSTRZ,ISYMRU,LISTRU,IDLSTRU,
     *                 IDUMMY,CDUMMY,IDUMMY,.FALSE.,
     *                 WORK(KLAMP0),WORK(KLAMH0),WORK(KEND1),LWRK1,
     *                 LU3SRTR,FN3SRTR,LUCKJDRZU,FNCKJDRZU)
C
      CALL CC3_SORT1(WORK(KEND1),LWRK1,2,ISINT1RZU,LU3SRTR,FN3SRTR,
     *               LUDELDRZU,FNDELDRZU,IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *               IDUMMY,CDUMMY)
C
      CALL CC3_SINT(WORK(KLAMH0),WORK(KEND1),LWRK1,ISINT1RZU,
     *              LUDELDRZU,FNDELDRZU,LUDKBCRZU,FNDKBCRZU)

C
C----------------------------------------
C     Integrals [H,T1ZU] where ZU is LISTR
C     (one-index transformed with second-order amplitudes)
C----------------------------------------
C
      ISINT1R2 = MULD2H(ISINT1,ISYMZU)
      ISINT2R2 = MULD2H(ISINT2,ISYMZU)
C
      CALL CC3_BARINT(WORK(KT1ZU),ISYMZU,WORK(KLAMP0),
     *                WORK(KLAMH0),WORK(KEND1),LWRK1,
     *                LU3SRTR,FN3SRTR,LUCKJDR2,FNCKJDR2)
C
      CALL CC3_SORT1(WORK(KEND1),LWRK1,2,ISINT1R2,LU3SRTR,FN3SRTR,
     *               LUDELDR2,FNDELDR2,IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *               IDUMMY,CDUMMY)
C
      CALL CC3_SINT(WORK(KLAMH0),WORK(KEND1),LWRK1,ISINT1R2,
     *              LUDELDR2,FNDELDR2,LUDKBCR2,FNDKBCR2)
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
            call quit('Insufficient space (T3) in CC3_ADENVIR_CUB')
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
      ISINT2L1R = MULD2H(ISYML1R,ISINT2)
      ISYFCKL1R = MULD2H(ISYMOP,ISYML1R)

      KXIAJB   = KEND1
      KEND1   = KXIAJB  + NT2AM(ISYM0)

      KT3BOG1 = KEND1
      KT3BOL1 = KT3BOG1 + NTRAOC(ISYM0)
      KT3BOG2 = KT3BOL1 + NTRAOC(ISYM0)
      KT3BOL2 = KT3BOG2 + NTRAOC(ISYM0)
      KT3OG1  = KT3BOL2 + NTRAOC(ISYM0)
      KT3OG2 = KT3OG1  + NTRAOC(ISINT2)
      KLAMPL1R  = KT3OG2 + NTRAOC(ISINT2)
      KLAMHL1R  = KLAMPL1R  + NLAMDT
      KEND1   = KLAMHL1R  + NLAMDT
      LWRK1   = LWORK   - KEND1
C
      KW3ZOGZ1 = KEND1
      KEND1   = KW3ZOGZ1 + NTRAOC(ISINT2RZ)
C
      KWZUOGR21 = KEND1
      KEND1     = KWZUOGR21 + NTRAOC(ISINT2RZU)
C
      KW3UOGU1 = KEND1
      KEND1    = KW3UOGU1 + NTRAOC(ISINT2RU)
C
      KW3ZUOGZU1 = KEND1
      KEND1      = KW3ZUOGZU1 + NTRAOC(ISINT2RZU)
C
      KFOCKL1RCK    = KEND1
      KW3BXOG1   = KFOCKL1RCK    + NT1AM(ISYFCKL1R)
      KW3BXOL1   = KW3BXOG1   + NTRAOC(ISYM0)
      KW3BXOGX1   = KW3BXOL1   + NTRAOC(ISYM0)
      KW3BXOLX1   = KW3BXOGX1   + NTRAOC(ISINT2L1R)
      KEND1      = KW3BXOLX1   + NTRAOC(ISINT2L1R)
      LWRK1      = LWORK      - KEND1
C
      KT2L1R = KEND1 
      KEND1  = KT2L1R + NT2SQ(ISYML1R)
      LWRK1      = LWORK      - KEND1
C
      KT1L1R  = KEND1
      KEND2  = KT1L1R + NT1AM(ISYML1R)
      LWRK2   = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in CC3_ADENVIR_CUB')
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
         IOPT = 3
         CALL GET_T1_T2(IOPT,.TRUE.,WORK(KT1L1R),WORK(KT2L1R),LISTL1R,
     *                  IDLSTL1R,
     *                  ISYML1R,WORK(KEND2),LWRK2)
         CALL CC3LR_MFOCK(WORK(KFOCKL1RCK),WORK(KT1L1R),WORK(KXIAJB),
     *                    ISYFCKL1R)
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
      CALL INTOCC_T3BARX(.FALSE.,
     *                   LUTOC,FNTOC,ISYMOP,WORK(KLAMH0),ISYM0,ISINT1,
     *                   WORK(KLAMHL1R),ISYML1R,ISINT1L1R,
     *                   WORK(KW3BXOG1),
     *                   WORK(KW3BXOL1),WORK(KW3BXOGX1),WORK(KW3BXOLX1),
     *                   WORK(KEND1),LWRK1)
C
C------------------------------------------------------------------
C     Read occupied integrals [H,T1Z] where Z is LISTRZ (used in WZ)
C-----------------------------------------------------------------
C
      CALL INTOCC_T3X(LUCKJDRZ,FNCKJDRZ,WORK(KLAMP0),ISINT2RZ,
     *                WORK(KW3ZOGZ1),WORK(KEND1),LWRK1)

C
C------------------------------------------------------------------
C     Read occupied integrals [H,T1ZU] (used in WZU)
C-----------------------------------------------------------------
C
      CALL INTOCC_T3X(LUCKJDR2,FNCKJDR2,WORK(KLAMP0),ISINT2RZU,
     *                WORK(KWZUOGR21),WORK(KEND1),LWRK1)
C
C------------------------------------------------------------------
C     Read occupied integrals [H,T1U] where U is LISTRU (used in WU)
C-----------------------------------------------------------------
C
      CALL INTOCC_T3X(LUCKJDRU,FNCKJDRU,WORK(KLAMP0),ISINT2RU,
     *                WORK(KW3UOGU1),WORK(KEND1),LWRK1)

C
C------------------------------------------------------------------
C     Read occupied integrals [[H,T1Z],T1U] (used in WZU)
C-----------------------------------------------------------------
C
      CALL INTOCC_T3X(LUCKJDRZU,FNCKJDRZU,WORK(KLAMP0),ISINT2RZU,
     *                WORK(KW3ZUOGZU1),WORK(KEND1),LWRK1)

C
C---------------------------------------------
C     Open files for Tbar and W intermediates:
C---------------------------------------------
C
      LUT3 = -1
      LUWBMAT = -1
      LUTHETA = -1
      LUWZU = -1

      CALL WOPEN2(LUT3,FNT3,64,0)
      CALL WOPEN2(LUWBMAT,FNWBMAT,64,0)
      CALL WOPEN2(LUTHETA,FNTHETA,64,0)
      CALL WOPEN2(LUWZU,FNWZU,64,0)
C
C----------------------------
C     Loop over D
C----------------------------
C
      DO ISYMD = 1,NSYM

         ISYCKBD0  = MULD2H(ISYMD,ISYM0)
         ISYCKBDL1R  = MULD2H(ISYMD,ISYML1R)
         ISYCKBDR1Z  = MULD2H(ISYMD,ISINT2RZ)
         ISYCKBDR1U  = MULD2H(ISYMD,ISINT2RU)
         ISYCKBDR1ZU  = MULD2H(ISYMD,ISINT2RZU)
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

            KW3BXVDGX1  = KEND3
            KW3BXVDGX2  = KW3BXVDGX1  + NCKATR(ISYCKBDL1R)
            KW3BXVDLX1  = KW3BXVDGX2  + NCKATR(ISYCKBDL1R)
            KW3BXVDLX2  = KW3BXVDLX1  + NCKATR(ISYCKBDL1R)
            KEND3     = KW3BXVDLX2  + NCKATR(ISYCKBDL1R)
            LWRK3     = LWORK     - KEND3
C
            KW3ZVDGZ1  = KEND3
            KEND3    = KW3ZVDGZ1  + NCKATR(ISYCKBDR1Z)
            LWRK3    = LWORK    - KEND3
C
            KWZUVDGR21 = KEND3
            KEND3      = KWZUVDGR21 + NCKATR(ISYCKBDR1ZU)
            LWRK3    = LWORK    - KEND3
C
            KT3VDGZ3 = KEND3
            KEND3    = KT3VDGZ3 + NCKATR(ISYCKBDR1Z)
            LWRK3    = LWORK    - KEND3
C
            KWZUVDGR23 = KEND3
            KEND3      = KWZUVDGR23 + NCKATR(ISYCKBDR1ZU)
            LWRK3    = LWORK    - KEND3
C
            KW3UVDGU1 = KEND3
            KEND3     = KW3UVDGU1 + NCKATR(ISYCKBDR1U)
            LWRK3    = LWORK    - KEND3
C
            KW3ZUVDGZU1 = KEND3
            KEND3       = KW3ZUVDGZU1 + NCKATR(ISYCKBDR1ZU)
            LWRK3    = LWORK    - KEND3
C
            KT3VDGU3 = KEND3
            KEND3    = KT3VDGU3 + NCKATR(ISYCKBDR1U)
            LWRK3    = LWORK    - KEND3
C
            KT3VDGZU3 = KEND3
            KEND3     = KT3VDGZU3 + NCKATR(ISYCKBDR1ZU)
            LWRK3    = LWORK    - KEND3
C
            MAXX1 = MAX(NCKA(ISYCKBD0),NCKA(ISYCKBDL1R))
            MAXX2 = MAX(MAXX1,NCKA(ISYCKBDR1Z))
            MAXX3 = MAX(MAXX2,NCKA(ISYCKBDR1U))
            MAXX4 = MAX(MAXX3,NCKA(ISYCKBDR1ZU))
C
            KINTVI = KEND3
            KTRVI6 = KINTVI + MAXX4
            KEND4  = KTRVI6 + NCKATR(ISYCKBD0) 
            LWRK4  = LWORK  - KEND4
           
            IF (LWRK4 .LT. 0) THEN
               WRITE(LUPRI,*) 'Memory available : ',LWORK
               WRITE(LUPRI,*) 'Memory needed    : ',KEND4
               CALL QUIT('Insufficient space in CC3_ADENVIR_CUB')
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
            CALL INTVIR_T3BARX_D(.FALSE.,
     *                           ISYMOP,LU3VI,FN3VI,LU3VI2,FN3VI2,
     *                           LU3FOP,FN3FOP,LU3FOP2,FN3FOP2,
     *                           WORK(KW3BXVDGX1),WORK(KW3BXVDG1),
     *                           WORK(KW3BXVDGX2),WORK(KW3BXVDG2),
     *                           WORK(KW3BXVDLX1),WORK(KW3BXVDL1),
     *                           WORK(KW3BXVDLX2),WORK(KW3BXVDL2),
     *                           WORK(KLAMPL1R),ISYML1R,WORK(KLAMP0),
     *                           ISYM0,ISYMD,D,WORK(KEND4),LWRK4)
C
C-----------------------------------------------------------------------
C        Read virtual integrals [H,T1Z] where Z is LISTRZ (used in WZ)
C-----------------------------------------------------------------------
C
            IOFF = ICKBD(ISYCKBDR1Z,ISYMD) + NCKATR(ISYCKBDR1Z)*(D - 1) 
     *           + 1
            IF (NCKATR(ISYCKBDR1Z) .GT. 0) THEN
               CALL GETWA2(LUDKBCRZ,FNDKBCRZ,WORK(KW3ZVDGZ1),IOFF,
     &                     NCKATR(ISYCKBDR1Z))
            ENDIF

C
C-----------------------------------------------------------------------
C        Read virtual integrals [H,T1ZU] (used in WZU)
C-----------------------------------------------------------------------
C
            IOFF = ICKBD(ISYCKBDR1ZU,ISYMD) + NCKATR(ISYCKBDR1ZU)*(D-1) 
     *           + 1
            IF (NCKATR(ISYCKBDR1ZU) .GT. 0) THEN
               CALL GETWA2(LUDKBCR2,FNDKBCR2,WORK(KWZUVDGR21),IOFF,
     &                     NCKATR(ISYCKBDR1ZU))
            ENDIF
C
C-----------------------------------------------------------------------
C        Read virtual integrals [H,T1U] where U is LISTRU (used in WU)
C-----------------------------------------------------------------------
C
            IOFF = ICKBD(ISYCKBDR1U,ISYMD) + NCKATR(ISYCKBDR1U)*(D - 1)
     *           + 1
            IF (NCKATR(ISYCKBDR1U) .GT. 0) THEN
               CALL GETWA2(LUDKBCRU,FNDKBCRU,WORK(KW3UVDGU1),IOFF,
     &                     NCKATR(ISYCKBDR1U))
            ENDIF

C
C-----------------------------------------------------------------------
C        Read virtual integrals [[H,T1Z],T1U] (used in WZU)
C-----------------------------------------------------------------------
C
            IOFF = ICKBD(ISYCKBDR1ZU,ISYMD) + NCKATR(ISYCKBDR1ZU)*(D-1)
     *           + 1
            IF (NCKATR(ISYCKBDR1ZU) .GT. 0) THEN
               CALL GETWA2(LUDKBCRZU,FNDKBCRZU,WORK(KW3ZUVDGZU1),IOFF,
     &                     NCKATR(ISYCKBDR1ZU))
            ENDIF

C
C--------------------------------------------------------------------
C           Read virtual integrals [H,T1Z] where Z is LISTRZ (used in W^Z)
C--------------------------------------------------------------------
C
            IF (NCKA(ISYCKBDR1Z) .GT. 0) THEN
               IOFF = ICKAD(ISYCKBDR1Z,ISYMD) +
     &                NCKA(ISYCKBDR1Z)*(D - 1) + 1
               CALL GETWA2(LUDELDRZ,FNDELDRZ,WORK(KINTVI),IOFF,
     *              NCKA(ISYCKBDR1Z))
            ENDIF
C
            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KT3VDGZ3),
     *                       WORK(KLAMH0),ISYMD,D,ISINT2RZ,
     *                       WORK(KEND4),LWRK4)

C
C--------------------------------------------------------------------
C           Read virtual integrals [H,T1ZU] (used in W^ZU)
C--------------------------------------------------------------------
C
            IF (NCKA(ISYCKBDR1ZU) .GT. 0) THEN
               IOFF = ICKAD(ISYCKBDR1ZU,ISYMD) +
     &                NCKA(ISYCKBDR1ZU)*(D - 1) + 1
               CALL GETWA2(LUDELDR2,FNDELDR2,WORK(KINTVI),IOFF,
     *              NCKA(ISYCKBDR1ZU))
            ENDIF
C
            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KWZUVDGR23),
     *                       WORK(KLAMH0),ISYMD,D,ISINT2RZU,
     *                       WORK(KEND4),LWRK4)

C
C--------------------------------------------------------------------
C           Read virtual integrals [H,T1U] where U is LISTRU (used in W^U)
C--------------------------------------------------------------------
C
            IF (NCKA(ISYCKBDR1U) .GT. 0) THEN
               IOFF = ICKAD(ISYCKBDR1U,ISYMD) +
     &                NCKA(ISYCKBDR1U)*(D - 1) + 1
               CALL GETWA2(LUDELDRU,FNDELDRU,WORK(KINTVI),IOFF,
     *              NCKA(ISYCKBDR1U))
            ENDIF    
C
            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KT3VDGU3),
     *                       WORK(KLAMH0),ISYMD,D,ISINT2RU,
     *                       WORK(KEND4),LWRK4)
C
C--------------------------------------------------------------------
C           Read virtual integrals [[H,T1Z],T1U] (used in W^ZU)
C--------------------------------------------------------------------
C
            IF (NCKA(ISYCKBDR1ZU) .GT. 0) THEN
               IOFF = ICKAD(ISYCKBDR1ZU,ISYMD) +
     &                NCKA(ISYCKBDR1ZU)*(D - 1) + 1
               CALL GETWA2(LUDELDRZU,FNDELDRZU,WORK(KINTVI),IOFF,
     *              NCKA(ISYCKBDR1ZU))
            ENDIF    
C
            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KT3VDGZU3),
     *                       WORK(KLAMH0),ISYMD,D,ISINT2RZU,
     *                       WORK(KEND4),LWRK4)



C
            DO ISYMB = 1,NSYM

               ISYALJB0  = MULD2H(ISYMB,ISYM0)
               ISYALJD0 = MULD2H(ISYMD,ISYM0)
               ISYALJBL1  = MULD2H(ISYMB,ISYML1)
               ISYALJDL1 = MULD2H(ISYMD,ISYML1)
               ISYMBD  = MULD2H(ISYMD,ISYMB)
               ISCKIJ  = MULD2H(ISYMBD,ISYM0)
               ISWBMAT  = MULD2H(ISCKIJ,ISYML1)
               ISWMATZ  = MULD2H(ISCKIJ,ISYMRZ)
               ISWMATU  = MULD2H(ISCKIJ,ISYMRU)
               ISWMATZU  = MULD2H(ISWMATZ,ISYMRU)
               ISYCKD  = MULD2H(ISYM0,ISYMB)
C
               ISYCKDBR1Z  = MULD2H(ISYMB,ISINT2RZ)
               ISYCKDBR1U  = MULD2H(ISYMB,ISINT2RU)
               ISYCKDBR1ZU  = MULD2H(ISYMB,ISINT2RZU)

C              Can use kend3 since we do not need the integrals anymore.
               KSMAT2     = KEND3
               KUMAT2     = KSMAT2    + NCKIJ(ISCKIJ)
               KDIAG      = KUMAT2    + NCKIJ(ISCKIJ)
               KDIAGWB     = KDIAG     + NCKIJ(ISCKIJ)
               KDIAGWZ     = KDIAGWB    + NCKIJ(ISWBMAT)
               KINDSQ     = KDIAGWZ    + NCKIJ(ISWMATZ)
               KINDSQWB    = KINDSQ    + (6*NCKIJ(ISCKIJ) - 1)/IRAT + 1
               KINDSQWZ    = KINDSQWB  + (6*NCKIJ(ISWBMAT) - 1)/IRAT + 1
               KINDEX     = KINDSQWZ   + (6*NCKIJ(ISWMATZ) - 1)/IRAT + 1
               KINDEX2    = KINDEX    + (NCKI(ISYALJB0)  - 1)/IRAT + 1
               KINDEXBL1   = KINDEX2   + (NCKI(ISYALJD0) - 1)/IRAT + 1
               KINDEXDL1   = KINDEXBL1 + (NCKI(ISYALJBL1)  - 1)/IRAT + 1
               KTMAT      = KINDEXDL1  + (NCKI(ISYALJDL1) - 1)/IRAT + 1
               KT3MAT     = KTMAT    + MAX(NCKIJ(ISCKIJ),NCKIJ(ISWBMAT))
               KW3MATZ      = KT3MAT    + NCKIJ(ISCKIJ)
C
               KWTEMP     = KW3MATZ      + NCKIJ(ISWMATZ)
               KEND4      = KWTEMP       + NCKIJMAX
               LWRK4      = LWORK        - KEND4
C
               KW3MATU  = KEND4
               KINDSQWU = KW3MATU  + NCKIJ(ISWMATU)
               KDIAGWU  = KINDSQWU + (6*NCKIJ(ISWMATU) - 1)/IRAT + 1
               KEND4    = KDIAGWU  + NCKIJ(ISWMATU)
               LWRK4    = LWORK    - KEND4

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
               KWMATZU = KEND4
               KEND4   = KWMATZU + NCKIJ(ISWMATZU)
               LWRK4   = LWORK   - KEND4
C
               KWMATZUD = KEND4
               KEND4    = KWMATZUD + NCKIJ(ISWMATZU)
               LWRK4   = LWORK   - KEND4
C
               KINDSQWZU = KEND4
               KDIAGWZU  = KINDSQWZU + (6*NCKIJ(ISWMATZU) - 1)/IRAT + 1
               KEND4     = KDIAGWZU  + NCKIJ(ISWMATZU)
               LWRK4   = LWORK   - KEND4
C
               KW3ZVDGZ2 = KEND4
               KEND4   = KW3ZVDGZ2 + NCKATR(ISYCKDBR1Z)
C
               KWZUVDGR22 = KEND4
               KEND4      = KWZUVDGR22 + NCKATR(ISYCKDBR1ZU)
C
               KW3UVDGU2 = KEND4
               KEND4     = KW3UVDGU2 + NCKATR(ISYCKDBR1U)
C
               KW3ZUVDGZU2 = KEND4
               KEND4       = KW3ZUVDGZU2 + NCKATR(ISYCKDBR1ZU)
C
               KT3VBGZ3 = KEND4
               KEND4    = KT3VBGZ3 + NCKATR(ISYCKDBR1Z)
C
               KWZUVBGR23 = KEND4
               KEND4      = KWZUVBGR23 + NCKATR(ISYCKDBR1ZU)
C
               KT3VBGU3 = KEND4
               KEND4    = KT3VBGU3 + NCKATR(ISYCKDBR1U)
C
               KT3VBGZU3 = KEND4
               KEND4     = KT3VBGZU3 + NCKATR(ISYCKDBR1ZU)
C
               MAXX1   = MAX(NCKA(ISYCKDBR1Z),NCKA(ISYCKDBR1U))
               MAXX2   = MAX(MAXX1,NCKA(ISYCKDBR1ZU))
C
               FKW3BXVDG1  = KEND4
               FKW3BXVDG2  = FKW3BXVDG1  + NCKATR(ISYALJB0)
               FKW3BXVDL1  = FKW3BXVDG2  + NCKATR(ISYALJB0)
               FKW3BXVDL2  = FKW3BXVDL1  + NCKATR(ISYALJB0)
               KEND4     = FKW3BXVDL2  + NCKATR(ISYALJB0)
               LWRK4     = LWORK     - KEND4
C
               ISYCKBDL1R  = MULD2H(ISYMB,ISYML1R)
C
               FKW3BXVDGX1  = KEND4
               FKW3BXVDGX2  = FKW3BXVDGX1  + NCKATR(ISYCKBDL1R)
               FKW3BXVDLX1  = FKW3BXVDGX2  + NCKATR(ISYCKBDL1R)
               FKW3BXVDLX2  = FKW3BXVDLX1  + NCKATR(ISYCKBDL1R)
               KEND4     = FKW3BXVDLX2  + NCKATR(ISYCKBDL1R)
               LWRK4     = LWORK     - KEND4
C
               KINTVI  = KEND4
               KEND5   = KINTVI + MAXX2
               LWRK5   = LWORK   - KEND5

               IF (LWRK5 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Memory available : ',LWORK
                  WRITE(LUPRI,*) 'Memory needed    : ',KEND5
                  CALL QUIT('Insufficient space in CC3_ADENVIR_CUB')
               END IF
C
C
C              -------------------------------
C              Construct part of the diagonal.
C              -------------------------------
C
               CALL CC3_DIAG(WORK(KDIAG), WORK(KFOCKD),ISCKIJ)
               CALL CC3_DIAG(WORK(KDIAGWB),WORK(KFOCKD),ISWBMAT)
               CALL CC3_DIAG(WORK(KDIAGWZ),WORK(KFOCKD),ISWMATZ)
               CALL CC3_DIAG(WORK(KDIAGWU),WORK(KFOCKD),ISWMATU)
               CALL CC3_DIAG(WORK(KDIAGWZU),WORK(KFOCKD),ISWMATZU)

C
C              -----------------------
C              Construct index arrays.
C              -----------------------
C
               LENSQ  = NCKIJ(ISCKIJ)
               CALL CC3_INDSQ(WORK(KINDSQ),LENSQ,ISCKIJ)
               LENSQWB  = NCKIJ(ISWBMAT)
               CALL CC3_INDSQ(WORK(KINDSQWB),LENSQWB,ISWBMAT)
               LENSQWZ = NCKIJ(ISWMATZ)
               CALL CC3_INDSQ(WORK(KINDSQWZ),LENSQWZ,ISWMATZ) 
               LENSQWU = NCKIJ(ISWMATU)
               CALL CC3_INDSQ(WORK(KINDSQWU),LENSQWU,ISWMATU) 
               LENSQWZU = NCKIJ(ISWMATZU)
               CALL CC3_INDSQ(WORK(KINDSQWZU),LENSQWZU,ISWMATZU)

C
               CALL CC3_INDEX(WORK(KINDEX),ISYALJB0)
               CALL CC3_INDEX(WORK(KINDEX2),ISYALJD0)
               CALL CC3_INDEX(WORK(KINDEXBL1),ISYALJBL1)
               CALL CC3_INDEX(WORK(KINDEXDL1),ISYALJDL1)

               DO B = 1,NVIR(ISYMB)
C
                  CALL DZERO(WORK(KW3MATZ),NCKIJ(ISWMATZ))
                  CALL DZERO(WORK(KW3MATU),NCKIJ(ISWMATU))
                  CALL DZERO(WORK(KWMATZU),NCKIJ(ISWMATZU))
                  CALL DZERO(WORK(KWMATZUD),NCKIJ(ISWMATZU))
C
                  CALL INTVIR_T3BARX_D(.FALSE.,
     *                           ISYMOP,LU3VI,FN3VI,LU3VI2,FN3VI2,
     *                           LU3FOP,FN3FOP,LU3FOP2,FN3FOP2,
     *                           WORK(FKW3BXVDGX1),WORK(FKW3BXVDG1),
     *                           WORK(FKW3BXVDGX2),WORK(FKW3BXVDG2),
     *                           WORK(FKW3BXVDLX1),WORK(FKW3BXVDL1),
     *                           WORK(FKW3BXVDLX2),WORK(FKW3BXVDL2),
     *                           WORK(KLAMPL1R),ISYML1R,WORK(KLAMP0),
     *                           ISYM0,ISYMB,B,WORK(KEND5),LWRK5)
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
C           Read virtual integrals [H,T1Z] where Z is LISTRZ (used in WZ)
C--------------------------------------------------------------------
C
                  IOFF = ICKBD(ISYCKDBR1Z,ISYMB) +
     &                   NCKATR(ISYCKDBR1Z)*(B - 1) + 1
                  IF (NCKATR(ISYCKDBR1Z) .GT. 0) THEN
                     CALL GETWA2(LUDKBCRZ,FNDKBCRZ,WORK(KW3ZVDGZ2),IOFF,
     &                           NCKATR(ISYCKDBR1Z))
                  ENDIF
C
                  IOFF = ICKAD(ISYCKDBR1Z,ISYMB) +
     &                   NCKA(ISYCKDBR1Z)*(B - 1) + 1
                  IF (NCKA(ISYCKDBR1Z) .GT. 0) THEN
                     CALL GETWA2(LUDELDRZ,FNDELDRZ,WORK(KINTVI),IOFF,
     *                    NCKA(ISYCKDBR1Z))
                  ENDIF
C
                  CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KT3VBGZ3),
     *                             WORK(KLAMH0),ISYMB,B,ISINT2RZ,
     *                             WORK(KEND5),LWRK5)

C
C--------------------------------------------------------------------
C           Read virtual integrals [H,T1ZU] (used in WZU)
C--------------------------------------------------------------------
C
                  IOFF = ICKBD(ISYCKDBR1ZU,ISYMB) +
     &                   NCKATR(ISYCKDBR1ZU)*(B - 1) + 1
                  IF (NCKATR(ISYCKDBR1ZU) .GT. 0) THEN
                     CALL GETWA2(LUDKBCR2,FNDKBCR2,WORK(KWZUVDGR22),
     *                           IOFF,NCKATR(ISYCKDBR1ZU))
                  ENDIF
C
                  IOFF = ICKAD(ISYCKDBR1ZU,ISYMB) +
     &                   NCKA(ISYCKDBR1ZU)*(B - 1) + 1
                  IF (NCKA(ISYCKDBR1ZU) .GT. 0) THEN
                     CALL GETWA2(LUDELDR2,FNDELDR2,WORK(KINTVI),IOFF,
     *                    NCKA(ISYCKDBR1ZU))
                  ENDIF
C
                  CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KWZUVBGR23),
     *                             WORK(KLAMH0),ISYMB,B,ISINT2RZU,
     *                             WORK(KEND5),LWRK5)
C
C--------------------------------------------------------------------
C           Read virtual integrals [H,T1U] where U is LISTRU (used in WU)
C--------------------------------------------------------------------
C
                  IOFF = ICKBD(ISYCKDBR1U,ISYMB) +
     &                   NCKATR(ISYCKDBR1U)*(B - 1) + 1
                  IF (NCKATR(ISYCKDBR1U) .GT. 0) THEN
                     CALL GETWA2(LUDKBCRU,FNDKBCRU,WORK(KW3UVDGU2),IOFF,
     &                           NCKATR(ISYCKDBR1U))
                  ENDIF
C
                  IOFF = ICKAD(ISYCKDBR1U,ISYMB) +
     &                   NCKA(ISYCKDBR1U)*(B - 1) + 1
                  IF (NCKA(ISYCKDBR1U) .GT. 0) THEN
                     CALL GETWA2(LUDELDRU,FNDELDRU,WORK(KINTVI),IOFF,
     *                    NCKA(ISYCKDBR1U))
                  ENDIF
C
                  CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KT3VBGU3),
     *                             WORK(KLAMH0),ISYMB,B,ISINT2RU,
     *                             WORK(KEND5),LWRK5)

C
C--------------------------------------------------------------------
C           Read virtual integrals [[H,T1Z],T1U] (used in WZU)
C--------------------------------------------------------------------
C
                  IOFF = ICKBD(ISYCKDBR1ZU,ISYMB) +
     &                   NCKATR(ISYCKDBR1ZU)*(B - 1) + 1
                  IF (NCKATR(ISYCKDBR1ZU) .GT. 0) THEN
                     CALL GETWA2(LUDKBCRZU,FNDKBCRZU,WORK(KW3ZUVDGZU2),
     *                           IOFF,NCKATR(ISYCKDBR1ZU))
                  ENDIF
C
                  IOFF = ICKAD(ISYCKDBR1ZU,ISYMB) +
     &                   NCKA(ISYCKDBR1ZU)*(B - 1) + 1
                  IF (NCKA(ISYCKDBR1ZU) .GT. 0) THEN
                     CALL GETWA2(LUDELDRZU,FNDELDRZU,WORK(KINTVI),IOFF,
     *                    NCKA(ISYCKDBR1ZU))
                  ENDIF
C
                  CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KT3VBGZU3),
     *                             WORK(KLAMH0),ISYMB,B,ISINT2RZU,
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
     *                            WORK(KT3OG1),WORK(KINDEX),
     *                            WORK(KS3MAT3),WORK(KT3VBG1),
     *                            WORK(KT3VBG2),WORK(KINDEX2),
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
                  CALL GET_T3BAR0_BD(ISYM0,WORK(KL1AM),ISYM0,
     *                               WORK(KL2TP),ISYM0,WORK(KTMAT),
     *                               WORK(KFOCK0CK),WORK(KFOCKD),
     *                               WORK(KDIAG),WORK(KXIAJB),ISYM0,
     *                               ISYM0,WORK(KINDSQ),LENSQ,
     *                               WORK(KSMAT2),WORK(KT3BVDG1),
     *                               WORK(KT3BVDG2),WORK(KT3BVDL1),
     *                               WORK(KT3BVDL2),WORK(KT3BOG1),
     *                               WORK(KT3BOL1),WORK(KINDEX),
     *                               WORK(KSMAT4),WORK(KT3BVBG1),
     *                               WORK(KT3BVBG2),WORK(KT3BVBL1),
     *                               WORK(KT3BVBL2),WORK(KINDEX2),
     *                               WORK(KUMAT2),WORK(KT3BVDG3),
     *                               WORK(KT3BVDL3),WORK(KT3BOG2),
     *                               WORK(KT3BOL2),WORK(KUMAT4),
     *                               WORK(KT3BVBG3),WORK(KT3BVBL3),
     *                               ISYMB,B,ISYMD,D,ISCKIJ,
     *                               WORK(KEND5),LWRK5)
                  IF (DO_YMMAT) THEN
                     CALL WRITE_T3_DL(LUWBMATX,FNWBMATX,WORK(KTMAT),
     *                                ISYM0,ISYMD,ISYMB,B)
                  END IF

c
c       call sum_pt3(work(KTMAT),isymb,b,isymd,d,
c    *             ISYM0,work(kx3am),4)
C
                  CALL DSCAL(NCKIJ(ISCKIJ),1.0D0/3.0D0,WORK(KTMAT),1)
c
c       call sum_pt3(work(KTMAT),isymb,b,isymd,d,
c    *             ISYM0,work(kx3am),5)
                  CALL DSCAL(NCKIJ(ISCKIJ),-ONE,WORK(KTMAT),1)
C
c
c       call sum_pt3(work(KTMAT),isymb,b,isymd,d,
c    *             ISYM0,work(kx3am),4)
C
C--------------------------------------------------------
C                 Write WBMAT as KTMAT^D(ai,bj,l) to disc
C--------------------------------------------------------
                  CALL WRITE_T3_DL(LUWBMAT,FNWBMAT,WORK(KTMAT),
     *                             ISYM0,ISYMD,ISYMB,B)

                  CALL DZERO(WORK(KW3MATZ),NCKIJ(ISWMATZ))
                  CALL DZERO(WORK(KW3MATU),NCKIJ(ISWMATU))
C
*****************************************************************
*****************************************************************
*
* Now we prepare 
* theta^{abc}_{i-- j-- k--} = C^{abc}_{ijk} wZU^{abc}_{i-- j- k-}
*
*****************************************************************
*****************************************************************
C

C
C=====================================================================
C Start with wZU^{abc}_{i-- j- k-} =
C
C     = - [ P(ZU) {   U_{li} wZ^{abc}_{l- j- k-}                    (1)
C
C                   + U(Z)_{li} t{abc}_{ljk}                        (2)
C                
C                   + b^{abc}_{ijk}(U, t2Z, t20)                    (3)
C                  
C                   + A^{abc}_{ijk} (t2UZ)                          (4)
C
C                   + B^{abc}_{ijk} (t2U, t2Z) ]                    (5)
C
C        * 1 / (epsilon^{abc}_{ijk} - omega_Z - omega_U)
C
C  Permutation P(ZU) is explicit in the following !
C
C=====================================================================
C

C                 ---------
C                  TERM (1)
C                 ---------

C
C------------------------------------------------------
C Calculate wZ^{abc}_{l- j- k-}
C------------------------------------------------------
C
                  AIBJCK_PERM = 4 ! means that we transform ALL occupied
                                  ! indeces

                  ! <mu3|[Z,T30]|HF> occupied contribution 
                  CALL WXBD_O(AIBJCK_PERM,WORK(KT3MAT),ISCKIJ,
     *                        WORK(KFOCKRZ),ISYMRZ,
     *                        WORK(KW3MATZ),ISWMATZ,WORK(KEND5),LWRK5)
C
                  ! <mu3|[[Z,T2],T2]|HF> 
                  CALL WXBD_T2(AIBJCK_PERM,B,ISYMB,D,ISYMD,WORK(KT2TP),
     *                         ISYM0,WORK(KFOCKRZ),
     *                         ISYMRZ,WORK(KINDSQWZ),LENSQWZ,
     *                         WORK(KW3MATZ),ISWMATZ,WORK(KEND5),LWRK5)

                  !<mu3|[H^Z,T2]|HF> + <mu3|[H,T2^Z]|HF>
                  CALL WXBD_GROUND(AIBJCK_PERM,
     *                       WORK(KT2RZ),ISYMRZ,WORK(KWTEMP),
     *                       WORK(KT3VDG1),WORK(KT3VBG1),WORK(KT3VBG3),
     *                       WORK(KT3VDG3),
     *                       WORK(KT3OG1),ISINT2,
     *                       WORK(KW3MATZ),WORK(KEND5),LWRK5,
     *                       WORK(KINDSQWZ),LENSQWZ,
     *                       ISYMB,B,ISYMD,D)
C
                  CALL WXBD_GROUND(AIBJCK_PERM,
     *                       WORK(KT2TP),ISYM0,WORK(KWTEMP),
     *                       WORK(KW3ZVDGZ1),WORK(KW3ZVDGZ2),
     *                       WORK(KT3VBGZ3),WORK(KT3VDGZ3),
     *                       WORK(KW3ZOGZ1),ISINT2RZ,
     *                       WORK(KW3MATZ),WORK(KEND5),LWRK5,
     *                       WORK(KINDSQWZ),LENSQWZ,
     *                       ISYMB,B,ISYMD,D)

                  !Divide by the energy difference and
                  !remove the forbidden elements
                  CALL WBD_DIA(B,ISYMB,D,ISYMD,FREQRZ,ISWMATZ,
     *                         WORK(KW3MATZ),WORK(KDIAGWZ),WORK(KFOCKD))
                  CALL T3_FORBIDDEN(WORK(KW3MATZ),ISYMRZ,ISYMB,B,
     *                              ISYMD,D)

 
c       call sum_pt3(work(KW3MATZ),isymb,b,isymd,d,
c    *             ISWMATZ,work(kx3am),5)

C
C-------------------------------------------------------------------
C                 Contract wZ with U operator:
C                 wZU^{abc}_{i-- j- k-} = wZU^{abc}_{i-- j- k-} 
C                                       + U_{li} wZ^{abc}_{l- j- k-}
C-------------------------------------------------------------------
C
                  CALL WBD_O(WORK(KW3MATZ),ISWMATZ,WORK(KFOCKRU),ISYMRU,
     *                 WORK(KWMATZU),ISWMATZU,WORK(KEND5),LWRK5)

c do the same to get wZU^{abc}_{i- j- k--} (put in KWMATZUD)

                  AIBJCK_PERM = 3
                  CALL WXBD_O(AIBJCK_PERM,WORK(KW3MATZ),ISWMATZ,
     *                 WORK(KFOCKRU),ISYMRU,
     *                 WORK(KWMATZUD),ISWMATZU,WORK(KEND5),LWRK5)



C                 ---------
C                  TERM (2)
C                 ---------

C
C----------------------------------------------------------------
C                 wZU^{abc}_{i-- j- k-} = wZU^{abc}_{i-- j- k-}
C                                       +  U(Z)_{li} t{abc}_{ljk}
C----------------------------------------------------------------
C
                  !Calculate <mu3|[[U,T1Z],T30]|HF>
                  CALL WBD_O(WORK(KT3MAT),ISCKIJ,WORK(KFCKUZO),ISYMZU,
     *                 WORK(KWMATZU),ISWMATZU,WORK(KEND5),LWRK5)

c do the same to get wZU^{abc}_{i- j- k--} (put in KWMATZUD)

                  AIBJCK_PERM = 3
                  CALL WXBD_O(AIBJCK_PERM,WORK(KT3MAT),ISCKIJ,
     *                 WORK(KFCKUZO),ISYMZU,
     *                 WORK(KWMATZUD),ISWMATZU,WORK(KEND5),LWRK5)

C                 ---------
C                  TERM (3)
C                 ---------

C
C--------------------------------------------------------------------
C                 wZU^{abc}_{i-- j- k-} = wZU^{abc}_{i-- j- k-}
C                                       +  b^{abc}_{ijk}(U, t2Z, t20)
C--------------------------------------------------------------------
C
                  !Calculate <mu3|[[U,T2Z],T20]|HF>
                  T2XNET2Z = .TRUE.
                  CALL WBD_T2(T2XNET2Z,B,ISYMB,D,ISYMD,
     *                        WORK(KT2RZ),ISYMRZ,WORK(KT2TP),ISYM0,
     *                        WORK(KFOCKRU),ISYMRU,
     *                        WORK(KINDSQWZU),LENSQWZU,WORK(KWMATZU),
     *                        ISWMATZU,WORK(KEND5),LWRK5)

c do the same to get wZU^{abc}_{i- j- k--} (put in KWMATZUD)

                  T2XNET2Z = .TRUE.
                  AIBJCK_PERM = 3 
                  CALL WXBD_T2_CUB(T2XNET2Z,AIBJCK_PERM,B,ISYMB,D,ISYMD,
     *                        WORK(KT2RZ),ISYMRZ,WORK(KT2TP),ISYM0,
     *                        WORK(KFOCKRU),ISYMRU,
     *                        WORK(KINDSQWZU),LENSQWZU,WORK(KWMATZUD),
     *                        ISWMATZU,WORK(KEND5),LWRK5)

C                 ---------
C                  TERM (4) 
C                 ---------
                  !P(ZU) permutation does not apply here: see the formula
C
C------------------------------------------------------
C Calculate A^{abc}_{ijk} (t2UZ)
C------------------------------------------------------
C
                  !<mu3|[[H,T1^ZU],T2^0]|HF>
                  AIBJCK_PERM = 1 ! means wZU^{abc}_{i-- j- k-}
                  CALL WXBD_GROUND(AIBJCK_PERM,
     *                       WORK(KT2TP),ISYM0,WORK(KWTEMP),
     *                       WORK(KWZUVDGR21),WORK(KWZUVDGR22),
     *                       WORK(KWZUVBGR23),WORK(KWZUVDGR23),
     *                       WORK(KWZUOGR21),ISINT2RZU,
     *                       WORK(KWMATZU),WORK(KEND5),LWRK5,
     *                       WORK(KINDSQWZU),LENSQWZU,
     *                       ISYMB,B,ISYMD,D)

                  !<mu3|[H^0,T2^ZU]|HF>
                  AIBJCK_PERM = 1 ! means wZU^{abc}_{i-- j- k-}
                  CALL WXBD_GROUND(AIBJCK_PERM,
     *                       WORK(KT2ZU),ISYMZU,WORK(KWTEMP),
     *                       WORK(KT3VDG1),WORK(KT3VBG1),WORK(KT3VBG3),
     *                       WORK(KT3VDG3),
     *                       WORK(KT3OG1),ISINT2,
     *                       WORK(KWMATZU),WORK(KEND5),LWRK5,
     *                       WORK(KINDSQWZU),LENSQWZU,
     *                       ISYMB,B,ISYMD,D)

c do the same to get wZU^{abc}_{i- j- k--} (put in KWMATZUD)

                  !<mu3|[[H,T1^ZU],T2^0]|HF>
                  AIBJCK_PERM = 3 ! means wZU^{abc}_{i- j- k--}
                  CALL WXBD_GROUND(AIBJCK_PERM,
     *                       WORK(KT2TP),ISYM0,WORK(KWTEMP),
     *                       WORK(KWZUVDGR21),WORK(KWZUVDGR22),
     *                       WORK(KWZUVBGR23),WORK(KWZUVDGR23),
     *                       WORK(KWZUOGR21),ISINT2RZU,
     *                       WORK(KWMATZUD),WORK(KEND5),LWRK5,
     *                       WORK(KINDSQWZU),LENSQWZU,
     *                       ISYMB,B,ISYMD,D)

                  !<mu3|[H^0,T2^ZU]|HF>
                  AIBJCK_PERM = 3 ! means wZU^{abc}_{i- j- k--}
                  CALL WXBD_GROUND(AIBJCK_PERM,
     *                       WORK(KT2ZU),ISYMZU,WORK(KWTEMP),
     *                       WORK(KT3VDG1),WORK(KT3VBG1),WORK(KT3VBG3),
     *                       WORK(KT3VDG3),
     *                       WORK(KT3OG1),ISINT2,
     *                       WORK(KWMATZUD),WORK(KEND5),LWRK5,
     *                       WORK(KINDSQWZU),LENSQWZU,
     *                       ISYMB,B,ISYMD,D)

C                 ---------
C                  TERM (5)
C                 ---------

C
C------------------------------------------------------
C Calculate B^{abc}_{ijk} (t2U, t2Z)
C------------------------------------------------------
C

                  !<mu3|[H^U,T2^Z]|HF>
                  AIBJCK_PERM = 1 ! means wZU^{abc}_{i-- j- k-}
                  CALL WXBD_GROUND(AIBJCK_PERM,
     *                       WORK(KT2RZ),ISYMRZ,WORK(KWTEMP),
     *                       WORK(KW3UVDGU1),WORK(KW3UVDGU2),
     *                       WORK(KT3VBGU3),
     *                       WORK(KT3VDGU3),
     *                       WORK(KW3UOGU1),ISINT2RU,
     *                       WORK(KWMATZU),WORK(KEND5),LWRK5,
     *                       WORK(KINDSQWZU),LENSQWZU,
     *                       ISYMB,B,ISYMD,D)

                  !<mu3|[[[H,T1^Z],T1^U],T2^0]|HF> 

                  !P(ZU) permutation taken into account here simply by
                  ! skipping the factor 1/2 from the formula.
                  ! Thus there is no need to have this term again in the
                  ! "permutation" part of this routine.

                  AIBJCK_PERM = 1 ! means wZU^{abc}_{i-- j- k-}
                  CALL WXBD_GROUND(AIBJCK_PERM,
     *                       WORK(KT2TP),ISYM0,WORK(KWTEMP),
     *                       WORK(KW3ZUVDGZU1),WORK(KW3ZUVDGZU2),
     *                       WORK(KT3VBGZU3),WORK(KT3VDGZU3),
     *                       WORK(KW3ZUOGZU1),ISINT2RZU,
     *                       WORK(KWMATZU),WORK(KEND5),LWRK5,
     *                       WORK(KINDSQWZU),LENSQWZU,
     *                       ISYMB,B,ISYMD,D)


c do the same to get wZU^{abc}_{i- j- k--} (put in KWMATZUD)

                  !<mu3|[H^U,T2^Z]|HF>
                  AIBJCK_PERM = 3 !means wZU^{abc}_{i- j- k--} (put in KWMATZUD)
                  CALL WXBD_GROUND(AIBJCK_PERM,
     *                       WORK(KT2RZ),ISYMRZ,WORK(KWTEMP),
     *                       WORK(KW3UVDGU1),WORK(KW3UVDGU2),
     *                       WORK(KT3VBGU3),
     *                       WORK(KT3VDGU3),
     *                       WORK(KW3UOGU1),ISINT2RU,
     *                       WORK(KWMATZUD),WORK(KEND5),LWRK5,
     *                       WORK(KINDSQWZU),LENSQWZU,
     *                       ISYMB,B,ISYMD,D)
                  !<mu3|[[[H,T1^Z],T1^U],T2^0]|HF> 
                  AIBJCK_PERM = 3 !means wZU^{abc}_{i- j- k--} (put in KWMATZUD)
                  CALL WXBD_GROUND(AIBJCK_PERM,
     *                       WORK(KT2TP),ISYM0,WORK(KWTEMP),
     *                       WORK(KW3ZUVDGZU1),WORK(KW3ZUVDGZU2),
     *                       WORK(KT3VBGZU3),WORK(KT3VDGZU3),
     *                       WORK(KW3ZUOGZU1),ISINT2RZU,
     *                       WORK(KWMATZUD),WORK(KEND5),LWRK5,
     *                       WORK(KINDSQWZU),LENSQWZU,
     *                       ISYMB,B,ISYMD,D)







                  !Divide by the energy difference and
                  !remove the forbidden elements (here only for debugging)
 
c                 call wbd_dia(b,isymb,d,isymd,freqzu,iswmatzu,
c    *                        work(kwmatzu),work(kdiagwzu),work(kfockd))
c                 call t3_forbidden(work(kwmatzu),isymzu,isymb,b,
c    *                              isymd,d)

 
c       call sum_pt3(work(KWMATZU),isymb,b,isymd,d,
c    *             ISWMATZU,work(kx3am),5)


C                -------------------------------------
C                 Repeat the TERMS (1)--(3) to include 
C                 P(ZU) PERMUTATION explicitly
C                -------------------------------------


C                 ---------
C                  TERM (1) (permuted)
C                 ---------

C
C------------------------------------------------------
C Calculate wU^{abc}_{l- j- k-}
C------------------------------------------------------
C
                  AIBJCK_PERM = 4 ! means that we transform ALL occupied
                                  ! indeces

                  ! <mu3|[U,T30]|HF> occupied contribution 
                  CALL WXBD_O(AIBJCK_PERM,WORK(KT3MAT),ISCKIJ,
     *                        WORK(KFOCKRU),ISYMRU,
     *                        WORK(KW3MATU),ISWMATU,WORK(KEND5),LWRK5)
C
                  ! <mu3|[[U,T2],T2]|HF> 
                  CALL WXBD_T2(AIBJCK_PERM,B,ISYMB,D,ISYMD,WORK(KT2TP),
     *                         ISYM0,WORK(KFOCKRU),
     *                         ISYMRU,WORK(KINDSQWU),LENSQWU,
     *                         WORK(KW3MATU),ISWMATU,WORK(KEND5),LWRK5)

                  !<mu3|[H^U,T2]|HF> + <mu3|[H,T2^U]|HF>
                  CALL WXBD_GROUND(AIBJCK_PERM,
     *                       WORK(KT2RU),ISYMRU,WORK(KWTEMP),
     *                       WORK(KT3VDG1),WORK(KT3VBG1),WORK(KT3VBG3),
     *                       WORK(KT3VDG3),
     *                       WORK(KT3OG1),ISINT2,
     *                       WORK(KW3MATU),WORK(KEND5),LWRK5,
     *                       WORK(KINDSQWU),LENSQWU,
     *                       ISYMB,B,ISYMD,D)
C
                  CALL WXBD_GROUND(AIBJCK_PERM,
     *                       WORK(KT2TP),ISYM0,WORK(KWTEMP),
     *                       WORK(KW3UVDGU1),WORK(KW3UVDGU2),
     *                       WORK(KT3VBGU3),WORK(KT3VDGU3),
     *                       WORK(KW3UOGU1),ISINT2RU,
     *                       WORK(KW3MATU),WORK(KEND5),LWRK5,
     *                       WORK(KINDSQWU),LENSQWU,
     *                       ISYMB,B,ISYMD,D)

                  !Divide by the energy difference and
                  !remove the forbidden elements
                  CALL WBD_DIA(B,ISYMB,D,ISYMD,FREQRU,ISWMATU,
     *                         WORK(KW3MATU),WORK(KDIAGWU),WORK(KFOCKD))
                  CALL T3_FORBIDDEN(WORK(KW3MATU),ISYMRU,ISYMB,B,
     *                              ISYMD,D)

 
c       call sum_pt3(work(KW3MATU),isymb,b,isymd,d,
c    *             ISWMATU,work(kx3am),5)

C
C-------------------------------------------------------------------
C                 Contract wU with Z operator:
C                 wZU^{abc}_{i-- j- k-} = wZU^{abc}_{i-- j- k-} 
C                                       + Z_{li} wU^{abc}_{l- j- k-}
C-------------------------------------------------------------------
C
                  CALL WBD_O(WORK(KW3MATU),ISWMATU,WORK(KFOCKRZ),ISYMRZ,
     *                 WORK(KWMATZU),ISWMATZU,WORK(KEND5),LWRK5)


c do the same to get wZU^{abc}_{i- j- k--} (put in KWMATZUD)

                  AIBJCK_PERM = 3
                  CALL WXBD_O(AIBJCK_PERM,WORK(KW3MATU),ISWMATU,
     *                 WORK(KFOCKRZ),ISYMRZ,
     *                 WORK(KWMATZUD),ISWMATZU,WORK(KEND5),LWRK5)

C                 ---------
C                  TERM (2) (permuted)
C                 ---------

C
C----------------------------------------------------------------
C                 wZU^{abc}_{i-- j- k-} = wZU^{abc}_{i-- j- k-}
C                                       +  Z(U)_{li} t{abc}_{ljk}
C----------------------------------------------------------------
C
                  !Calculate <mu3|[[Z,T1U],T30]|HF>
                  CALL WBD_O(WORK(KT3MAT),ISCKIJ,WORK(KFCKZUO),ISYMZU,
     *                 WORK(KWMATZU),ISWMATZU,WORK(KEND5),LWRK5)


c do the same to get wZU^{abc}_{i- j- k--} (put in KWMATZUD)

                  AIBJCK_PERM = 3
                  CALL WXBD_O(AIBJCK_PERM,WORK(KT3MAT),ISCKIJ,
     *                 WORK(KFCKZUO),ISYMZU,
     *                 WORK(KWMATZUD),ISWMATZU,WORK(KEND5),LWRK5)

C                 ---------
C                  TERM (3) (permuted)
C                 ---------

C
C--------------------------------------------------------------------
C                 wZU^{abc}_{i-- j- k-} = wZU^{abc}_{i-- j- k-}
C                                       +  b^{abc}_{ijk}(Z, t2U, t20)
C--------------------------------------------------------------------
C
                  !Calculate <mu3|[[Z,T2U],T20]|HF>
                  T2XNET2Z = .TRUE.
                  CALL WBD_T2(T2XNET2Z,B,ISYMB,D,ISYMD,
     *                        WORK(KT2RU),ISYMRU,WORK(KT2TP),ISYM0,
     *                        WORK(KFOCKRZ),ISYMRZ,
     *                        WORK(KINDSQWZU),LENSQWZU,WORK(KWMATZU),
     *                        ISWMATZU,WORK(KEND5),LWRK5)


c do the same to get wZU^{abc}_{i- j- k--} (put in KWMATZUD)
                  T2XNET2Z = .TRUE.
                  AIBJCK_PERM = 3 
                  CALL WXBD_T2_CUB(T2XNET2Z,AIBJCK_PERM,B,ISYMB,D,ISYMD,
     *                     WORK(KT2RU),ISYMRU,
     *                     WORK(KT2TP),ISYM0,WORK(KFOCKRZ),ISYMRZ,
     *                     WORK(KINDSQWZU),LENSQWZU,WORK(KWMATZUD),
     *                     ISWMATZU,WORK(KEND5),LWRK5)

C                 ---------
C                  TERM (5) (permuted)
C                 ---------

C
C------------------------------------------------------
C Calculate B^{abc}_{ijk} (t2Z, t2U)
C------------------------------------------------------
C

                  !<mu3|[H^Z,T2^U]|HF>
                  AIBJCK_PERM = 1 ! means wZU^{abc}_{i-- j- k-}
                  CALL WXBD_GROUND(AIBJCK_PERM,
     *                       WORK(KT2RU),ISYMRU,WORK(KWTEMP),
     *                       WORK(KW3ZVDGZ1),WORK(KW3ZVDGZ2),
     *                       WORK(KT3VBGZ3),
     *                       WORK(KT3VDGZ3),
     *                       WORK(KW3ZOGZ1),ISINT2RZ,
     *                       WORK(KWMATZU),WORK(KEND5),LWRK5,
     *                       WORK(KINDSQWZU),LENSQWZU,
     *                       ISYMB,B,ISYMD,D)

c do the same to get wZU^{abc}_{i- j- k--} (put in KWMATZUD)

                  !<mu3|[H^Z,T2^U]|HF>
                  AIBJCK_PERM = 3 !means wZU^{abc}_{i- j- k--} (put in KWMATZUD)
                  CALL WXBD_GROUND(AIBJCK_PERM,
     *                       WORK(KT2RU),ISYMRU,WORK(KWTEMP),
     *                       WORK(KW3ZVDGZ1),WORK(KW3ZVDGZ2),
     *                       WORK(KT3VBGZ3),
     *                       WORK(KT3VDGZ3),
     *                       WORK(KW3ZOGZ1),ISINT2RZ,
     *                       WORK(KWMATZUD),WORK(KEND5),LWRK5,
     *                       WORK(KINDSQWZU),LENSQWZU,
     *                       ISYMB,B,ISYMD,D)




                  !Divide by the energy difference and
                  !remove the forbidden elements 
                  CALL WBD_DIA(B,ISYMB,D,ISYMD,FREQZU,ISWMATZU,
     *                        WORK(KWMATZU),WORK(KDIAGWZU),WORK(KFOCKD))
                  CALL T3_FORBIDDEN(WORK(KWMATZU),ISYMZU,ISYMB,B,
     *                              ISYMD,D)


c       call sum_pt3(work(KWMATZU),isymb,b,isymd,d,
c    *             ISWMATZU,work(kx3am),5)


c do the same for wZU^{abc}_{i- j- k--} (put in KWMATZUD)

                  !Divide by the energy difference and
                  !remove the forbidden elements 
                  CALL WBD_DIA(B,ISYMB,D,ISYMD,FREQZU,ISWMATZU,
     *                       WORK(KWMATZUD),WORK(KDIAGWZU),WORK(KFOCKD))
                  CALL T3_FORBIDDEN(WORK(KWMATZUD),ISYMZU,ISYMB,B,
     *                              ISYMD,D)

c       call sum_pt3(work(KWMATZUD),isymb,b,isymd,d,
c    *             ISWMATZU,work(kx3am),5)


c get now wtildeU^{abc}_{ijk} = (1 + 0.5 P(ck,ai) ) wZU^{abc}_{i-- j- k-}

                  CALL DAXPY(NCKIJ(ISWMATZU),HALF,WORK(KWMATZUD),1,
     *                       WORK(KWMATZU),1)

C-----------------------------------------------------------------------
C    Write WORK(KWMATZU) + 0.5*WORK(KWMATZUD) as KW3MATZU^D(ai,bj,l) to disc
C-----------------------------------------------------------------------
                  !To conform with noddy code
                  CALL DSCAL(NCKIJ(ISWMATZU),-ONE,WORK(KWMATZU),1)
C
                  CALL WRITE_T3_DL(LUTHETA,FNTHETA,WORK(KWMATZU),ISYMZU,
     *                             ISYMD,ISYMB,B)

c       call sum_pt3(work(KWMATZU),isymb,b,isymd,d,
c    *             ISWMATZU,work(kx3am),4)

C                 ...now KWMATZU and KWMATZUD can be reused...


C
*****************************************************************
*****************************************************************
*
* Now we prepare 
* wZU^{a- bc}_{i- j- k-} = w^{a- bc}_{i- j- k-} = theta^{a- bc}_{i- j- k-}
*
*****************************************************************
*****************************************************************
C

C
C=====================================================================
C wZU^{a- bc}_{i- j- k-} =
C
C     = - [ P(ZU) {   U_{ad} wZ^{dbc}_{i- j- k-}                    (1)
C
C                   + U(Z)_{ad} t{dbc}_{ijk}                        (2)
C                
C                   + U_{li} thetaZ^{a- bc}_{ljk}                   (3)
C                  
C                   + U_{lj} thetaZ^{a- bc}_{ilk}                   (4)
C
C                   + U_{lk} thetaZ^{a- bc}_{ijl} ]                 (5)
C
C        * 1 / (epsilon^{abc}_{ijk} - omega_Z - omega_U)
C
C  Permutation P(ZU) is explicit in the following !
C
C=====================================================================
C


C                 We will reuse here KWMATZU
                  CALL DZERO(WORK(KWMATZU),NCKIJ(ISWMATZU))

C                 ---------
C                  TERM (1)
C                 ---------

C
C--------------------------------------------------------------------
C                 wZU^{a- bc}_{i- j- k-} = wZU^{a- bc}_{i- j- k-}
C                                        + U_{ad} wZ^{dbc}_{i- j- k-}
C--------------------------------------------------------------------
C

C                 wZ^{abc}_{l- j- k-} is already there sitting in
C                 KW3MATZ array.

                  CALL WBD_V(WORK(KW3MATZ),ISWMATZ,WORK(KFOCKRU),ISYMRU,
     *                 WORK(KWMATZU),ISWMATZU,WORK(KEND5),LWRK5)

C                 ---------
C                  TERM (1) (permuted)
C                 ---------

C
C--------------------------------------------------------------------
C                 wZU^{a- bc}_{i- j- k-} = wZU^{a- bc}_{i- j- k-}
C                                        + Z_{ad} wU^{dbc}_{i- j- k-}
C--------------------------------------------------------------------
C

C                 wU^{abc}_{l- j- k-} is already there sitting in
C                 KW3MATU array.

                  CALL WBD_V(WORK(KW3MATU),ISWMATU,WORK(KFOCKRZ),ISYMRZ,
     *                 WORK(KWMATZU),ISWMATZU,WORK(KEND5),LWRK5)


C                 ---------
C                  TERM (2)
C                 ---------

C
C--------------------------------------------------------------------
C                 wZU^{a- bc}_{i- j- k-} = wZU^{a- bc}_{i- j- k-}
C                                        + U(Z)_{ad} t{dbc}_{ijk}
C--------------------------------------------------------------------
C

C                 t{dbc}_{ijk} is already there sitting in
C                 KT3MAT array.

                  !Calculate <mu3|[[U,T1Z],T30]|HF>
                  CALL WBD_V(WORK(KT3MAT),ISCKIJ,WORK(KFCKZUV),ISYMZU,
     *                 WORK(KWMATZU),ISWMATZU,WORK(KEND5),LWRK5)

C                 ---------
C                  TERM (2) (permuted)
C                 ---------

C
C--------------------------------------------------------------------
C                 wZU^{a- bc}_{i- j- k-} = wZU^{a- bc}_{i- j- k-}
C                                        + Z(U)_{ad} t{dbc}_{ijk}
C--------------------------------------------------------------------
C

C                 t{dbc}_{ijk} is already there sitting in
C                 KT3MAT array.

                  !Calculate <mu3|[[Z,T1U],T30]|HF>
                  CALL WBD_V(WORK(KT3MAT),ISCKIJ,WORK(KFCKUZV),ISYMZU,
     *                 WORK(KWMATZU),ISWMATZU,WORK(KEND5),LWRK5)

C                 ---------------------
C                  TERM (3) + (4) + (5)
C                 ---------------------

C
C-------------------------------------------------------------------------
C                 wZU^{a- bc}_{i- j- k-} = wZU^{a- bc}_{i- j- k-}
C                                        + U_{li} thetaZ^{a- bc}_{ljk} (3)
C                                        + U_{lj} thetaZ^{a- bc}_{ilk} (4)
C                                        + U_{lk} thetaZ^{a- bc}_{ijl} (5)
C-------------------------------------------------------------------------
C

C-------------------------------------------------------
C                 First we need thetaZ^{a- bc}_{ijk}...
C-------------------------------------------------------

C                 Let's reuse KW3MATZ array
                  CALL DZERO(WORK(KW3MATZ),NCKIJ(ISWMATZ))

                  ! <mu3|[Z,T30]|HF> virtual contribution 
                  CALL WBD_V(WORK(KT3MAT),ISCKIJ,
     *                        WORK(KFOCKRZ),ISYMRZ,
     *                        WORK(KW3MATZ),ISWMATZ,WORK(KEND5),LWRK5)

                  !Divide by the energy difference and
                  !remove the forbidden elements
                  CALL WBD_DIA(B,ISYMB,D,ISYMD,FREQRZ,ISWMATZ,
     *                         WORK(KW3MATZ),WORK(KDIAGWZ),WORK(KFOCKD))
                  CALL T3_FORBIDDEN(WORK(KW3MATZ),ISYMRZ,ISYMB,B,
     *                                ISYMD,D)


c       call sum_pt3(work(KW3MATZ),isymb,b,isymd,d,
c    *             ISWMATZ,work(kx3am),5)

C
C-------------------------------------------------------------------------
C                 Now contract thetaZ with U operator:
C                 wZU^{a- bc}_{i- j- k-} = wZU^{a- bc}_{i- j- k-}
C                                        + U_{li} thetaZ^{a- bc}_{ljk} 
C                                        + U_{lj} thetaZ^{a- bc}_{ilk} 
C                                        + U_{lk} thetaZ^{a- bc}_{ijl} 
C-------------------------------------------------------------------------
C

                  AIBJCK_PERM = 4 ! transform all occ indeces simultanously
                  CALL WXBD_O(AIBJCK_PERM,WORK(KW3MATZ),ISWMATZ,
     *                 WORK(KFOCKRU),ISYMRU,
     *                 WORK(KWMATZU),ISWMATZU,WORK(KEND5),LWRK5)

C                 ---------------------
C                  TERM (3) + (4) + (5) (permuted)
C                 ---------------------

C
C-------------------------------------------------------------------------
C                 wZU^{a- bc}_{i- j- k-} = wZU^{a- bc}_{i- j- k-}
C                                        + Z_{li} thetaU^{a- bc}_{ljk} (3)
C                                        + Z_{lj} thetaU^{a- bc}_{ilk} (4)
C                                        + Z_{lk} thetaU^{a- bc}_{ijl} (5)
C-------------------------------------------------------------------------
C

C-------------------------------------------------------
C                 First we need thetaU^{a- bc}_{ijk}...
C-------------------------------------------------------

C                 Let's reuse KW3MATU array
                  CALL DZERO(WORK(KW3MATU),NCKIJ(ISWMATU))

                  ! <mu3|[U,T30]|HF> virtual contribution 
                  CALL WBD_V(WORK(KT3MAT),ISCKIJ,
     *                        WORK(KFOCKRU),ISYMRU,
     *                        WORK(KW3MATU),ISWMATU,WORK(KEND5),LWRK5)

                  !Divide by the energy difference and
                  !remove the forbidden elements
                  CALL WBD_DIA(B,ISYMB,D,ISYMD,FREQRU,ISWMATU,
     *                         WORK(KW3MATU),WORK(KDIAGWU),WORK(KFOCKD))
                  CALL T3_FORBIDDEN(WORK(KW3MATU),ISYMRU,ISYMB,B,
     *                                ISYMD,D)


c       call sum_pt3(work(KW3MATU),isymb,b,isymd,d,
c    *             ISWMATU,work(kx3am),5)

C
C-------------------------------------------------------------------------
C                 Now contract thetaU with Z operator:
C                 wZU^{a- bc}_{i- j- k-} = wZU^{a- bc}_{i- j- k-}
C                                        + Z_{li} thetaU^{a- bc}_{ljk} 
C                                        + Z_{lj} thetaU^{a- bc}_{ilk} 
C                                        + Z_{lk} thetaU^{a- bc}_{ijl} 
C-------------------------------------------------------------------------
C

                  AIBJCK_PERM = 4 ! transform all occ indeces simultanously
                  CALL WXBD_O(AIBJCK_PERM,WORK(KW3MATU),ISWMATU,
     *                 WORK(KFOCKRZ),ISYMRZ,
     *                 WORK(KWMATZU),ISWMATZU,WORK(KEND5),LWRK5)

                  !Divide by the energy difference and
                  !remove the forbidden elements 
                  CALL WBD_DIA(B,ISYMB,D,ISYMD,FREQZU,ISWMATZU,
     *                        WORK(KWMATZU),WORK(KDIAGWZU),WORK(KFOCKD))
                  CALL T3_FORBIDDEN(WORK(KWMATZU),ISYMZU,ISYMB,B,
     *                              ISYMD,D)


c       call sum_pt3(work(KWMATZU),isymb,b,isymd,d,
c    *             ISWMATZU,work(kx3am),5)

C-----------------------------------------------------------------------
C    Write wZU^{a- bc}_{i- j- k-} to file
C-----------------------------------------------------------------------
                  !To conform with noddy code
                  CALL DSCAL(NCKIJ(ISWMATZU),-ONE,WORK(KWMATZU),1)
C
                  CALL WRITE_T3_DL(LUWZU,FNWZU,WORK(KWMATZU),ISYMZU,
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

               ENDDO   ! B
            ENDDO      ! ISYMB
C--------------------------------------------------------------------
C          Get y^M(imfN) <- y^M(imfN)+ sum_em y^t^DL(ei) tbar^DLN(em,f)
C--------------------------------------------------------------------
C
           IF (DO_YMMAT) THEN
              CALL MMATX(YMMAT,LISTL0,
     *                   LUWBMATX,FNWBMATX,ISYM0,
     *                   LUT2Y,FNT2Y,ISYMZU,
     *                   ISYMD,D,
     *                   WORK(KEND5),LWRK5)
           END IF

      !LUWBMAT contains now T30bar 
c            so sym of WBAR (11th keyword) is set to ISYM0
            QUADR = .TRUE.
            CALL CC_XI_DEN_ABIJ_CUB_T0(QUADR,DAB,DIJ,ISYDEN,
     *                          WORK(KL2L1),ISYML1,
     *                          ISYMRZ,WORK(KFOCKRZ),
     *                          ISYMRU,WORK(KFOCKRU),
     *                          ISYM0,ISYM0,ISYMZU,
     *                          LUT3,FNT3,LUWBMAT,FNWBMAT,
     *                          LUTHETA,FNTHETA,
     *                          LUWZU,FNWZU,
     *                          WORK(KFOCKD),FREQRZ,FREQRU,
     *                          WORK(KEND5),LWRK5,ISYMD,D)
C
         ENDDO       ! D
      ENDDO          ! ISYMD 

c      write(lupri,*) 'w3x (usual) in CC3_ADENVIR_CUB'
c      write(lupri,*) 'w3xD  in CC3_ADENVIR_CUB_t0'
c      write(lupri,*) 'w3x  in CC3_ADENVIR_CUB_t0'
c      write(lupri,*) 'w3x + 0.5w3xD in CC3_ADENVIR_CUB'
c      write(lupri,*) 'w3zu  in CC3_ADENVIR_CUB'
c      call print_pt3(work(kx3am),ISWMATZU,4)

C
C---------------------------------
C     Close the files
C---------------------------------
C
      IF (DO_YMMAT) THEN
         CALL WCLOSE2(LUT2Y,FNT2Y,'DELETE')
         CALL WCLOSE2(LUWBMATX,FNWBMATX,'DELETE')
      END IF
C
      CALL WCLOSE2(LUT3,FNT3,'DELETE')
      CALL WCLOSE2(LUWBMAT,FNWBMAT,'DELETE')
      CALL WCLOSE2(LUTHETA,FNTHETA,'DELETE')
      CALL WCLOSE2(LUWZU,FNWZU,'DELETE')
C
C--------------------------------
C     Close files for "response"
C--------------------------------
C
      CALL WCLOSE2(LU3SRTR,FN3SRTR,'DELETE')
      CALL WCLOSE2(LUCKJDRZ,FNCKJDRZ,'DELETE')
      CALL WCLOSE2(LUDELDRZ,FNDELDRZ,'DELETE')
      CALL WCLOSE2(LUDKBCRZ,FNDKBCRZ,'DELETE')
C
      CALL WCLOSE2(LUCKJDRU,FNCKJDRU,'DELETE')
      CALL WCLOSE2(LUDELDRU,FNDELDRU,'DELETE')
      CALL WCLOSE2(LUDKBCRU,FNDKBCRU,'DELETE')
C
      CALL WCLOSE2(LUCKJDRZU,FNCKJDRZU,'DELETE')
      CALL WCLOSE2(LUDELDRZU,FNDELDRZU,'DELETE')
      CALL WCLOSE2(LUDKBCRZU,FNDKBCRZU,'DELETE')
C
      CALL WCLOSE2(LUCKJDR2,FNCKJDR2,'DELETE')
      CALL WCLOSE2(LUDELDR2,FNDELDR2,'DELETE')
      CALL WCLOSE2(LUDKBCR2,FNDKBCR2,'DELETE')
C
C-------------
C     End
C-------------
C
C
      CALL QEXIT('CC3DENVCB')
C
      RETURN
      END
C  /* Deck cc3_xi_den_abij_cub_t0 */
      SUBROUTINE CC_XI_DEN_ABIJ_CUB_T0(QUADR,DAB,DIJ,ISYDEN,
     *                           L2L1,ISYML1,
     *                           ISYFCKX,FOCKX,
     *                           ISYFCKY,FOCKY,
     *                           ISYMT3,ISWMAT,ISTHETA,
     *                           LUT3,FNT3,LUWBMAT,FNWBMAT,
     *                           LUTHETA,FNTHETA,
     *                           LUWZU,FNWZU,
     *                           FOCKD,FREQX,FREQY,
     *                           WORK,LWORK,ISYMD,D)
C
C=========================================================================
C
C    QUADR has to be .FALSE. for linear response calculations
C
C    QUADR has to be .TRUE.  for cubic response calculations
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
C    Written by Filip Pawlowski, Fall 2003, Aarhus
C=========================================================================
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <cc3t3d.h>
C
      INTEGER ISYMT3,ISWMAT,LUT3,LUWBMAT,LWORK,ISYMD,ISYDEN
      INTEGER ISYML,ISYMDL,ISWMATDL,ISYMT3DL,ISYMN,ISYEMF,ISYMBN,ISYMEM
      INTEGER ISYMB,ISYMF,ISYMFI,ISYMI,ISYEMB,ISYMFN
      INTEGER KT3,KWMAT,KEND1,LWRK1
      INTEGER KOFF1,KOFF2,KOFF3,KBN,KFN
      INTEGER NTOTEM,NTOTF,NNEMF
      INTEGER IADR
C
      INTEGER ISTHETA,ISYFCKX,ISYFCKY,LUTHETA,LUWZU
      INTEGER ISTHETADL,ISTHETAFX,ISTHETAFY
      INTEGER KTHETA,KTHETAF,KWZU
      INTEGER KFI
      INTEGER IOPT
      INTEGER MAXX1
C
      INTEGER ISYMJ,ISYMFJ,KFJ
      INTEGER ISYMM,ISYME
C
      INTEGER ISYML1
C
      LOGICAL QUADR
      LOGICAL TRANSPOSEW
C
      CHARACTER*(*) FNT3,FNWBMAT,FNTHETA,FNWZU
C
#if defined (SYS_CRAY)
      REAL DAB(*),DIJ(*),WORK(LWORK),ONE,HALF
      REAL FOCKX(*),FOCKY(*),L2L1(*),FOCKD(*),FREQX,FREQY
      REAL XNORMVAL,DDOT,FREQXY
#else
      DOUBLE PRECISION DAB(*),DIJ(*),WORK(LWORK),ONE,HALF
      DOUBLE PRECISION FOCKX(*),FOCKY(*),L2L1(*),FOCKD(*),FREQX,FREQY
      DOUBLE PRECISION XNORMVAL,DDOT,FREQXY
#endif
C
      PARAMETER(ONE = 1.0D0, HALF = 0.5D0)
C
      CALL QENTER('DENABIJC')
C
      DO ISYML = 1,NSYM
C
         ISYMDL = MULD2H(ISYMD,ISYML)
         ISWMATDL = MULD2H(ISWMAT,ISYMDL)
         ISYMT3DL = MULD2H(ISYMT3,ISYMDL)
         IF (QUADR) THEN
            ISTHETADL = MULD2H(ISTHETA,ISYMDL)
            ISTHETAFX  = MULD2H(ISYMT3DL,ISYFCKX)
            ISTHETAFY  = MULD2H(ISYMT3DL,ISYFCKY)
         END IF
C
         KT3  = 1
         KWMAT  = KT3 + NT2SQ(ISYMT3DL)
         KEND1 = KWMAT + NT2SQ(ISWMATDL)
         LWRK1  = LWORK - KEND1
C
         IF (QUADR) THEN
C
            MAXX1 = MAX(NT2SQ(ISTHETAFX),NT2SQ(ISTHETAFY))
C
            KTHETA  = KEND1 
            KTHETAF = KTHETA  + NT2SQ(ISTHETADL)
            KEND1   = KTHETAF + MAX(MAXX1,NT2SQ(ISTHETADL))
            LWRK1   = LWORK   - KEND1
C
            KWZU    = KEND1
            KEND1   = KWZU + NT2SQ(ISTHETADL)
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
C            ---------------------------------------------
C            4ht line of Eq. 62 (second cont)
C            ---------------------------------------------

C
C               KTHETAF(De- f)_(lmi) = KT3 * FOCKX
C
               ! KTHETAF is recycled here
               CALL DZERO(WORK(KTHETAF),NT2SQ(ISTHETAFX))
               IOPT = 3
               CALL WXDL_V(IOPT,WORK(KT3),ISYMT3DL,FOCKX,ISYFCKX,
     *                     WORK(KTHETAF),ISTHETAFX,WORK(KEND1),LWRK1)
               ! Divide it by orbital energy difference and remove the 
               ! forbidden elements
               CALL W3DL_DIA(WORK(KTHETAF),ISTHETAFX,ISYML,L,ISYMD,D,
     *                       FOCKD,FREQX)
               CALL T3_FORBIDDEN_DL(WORK(KTHETAF),ISYFCKX,ISYMD,D,
     *                              ISYML,L)

C
C               KTHETA(De- f-)_(lmi) = KTHETAF(De- f)_(lmi) * FOCKY
C
               ! KTHETA is recycled here
               CALL DZERO(WORK(KTHETA),NT2SQ(ISTHETADL))
               IOPT = 1
               CALL WXDL_V(IOPT,WORK(KTHETAF),ISTHETAFX,FOCKY,ISYFCKY,
     *                     WORK(KTHETA),ISTHETADL,WORK(KEND1),LWRK1)

C
C               KTHETAF(Def- )_(lmi) = KT3 * FOCKX
C
               ! KTHETAF is recycled here
               CALL DZERO(WORK(KTHETAF),NT2SQ(ISTHETAFX))
               IOPT = 1
               CALL WXDL_V(IOPT,WORK(KT3),ISYMT3DL,FOCKX,ISYFCKX,
     *                     WORK(KTHETAF),ISTHETAFX,WORK(KEND1),LWRK1)
               ! Divide it by orbital energy difference and remove the 
               ! forbidden elements
               CALL W3DL_DIA(WORK(KTHETAF),ISTHETAFX,ISYML,L,ISYMD,D,
     *                       FOCKD,FREQX)
               CALL T3_FORBIDDEN_DL(WORK(KTHETAF),ISYFCKX,ISYMD,D,
     *                              ISYML,L)

C
C               KTHETA(De- f-)_(lmi) = KTHETAF(Def- )_(lmi) * FOCKY
C
               IOPT = 3
               CALL WXDL_V(IOPT,WORK(KTHETAF),ISTHETAFX,FOCKY,ISYFCKY,
     *                     WORK(KTHETA),ISTHETADL,WORK(KEND1),LWRK1)


C
C              Include P(XY) permutation
C
C
C               KTHETAF(De- f)_(lmi) = KT3 * FOCKY
C
               ! KTHETAF is recycled here
               CALL DZERO(WORK(KTHETAF),NT2SQ(ISTHETAFY))
               IOPT = 3
               CALL WXDL_V(IOPT,WORK(KT3),ISYMT3DL,FOCKY,ISYFCKY,
     *                     WORK(KTHETAF),ISTHETAFY,WORK(KEND1),LWRK1)
               ! Divide it by orbital energy difference and remove the 
               ! forbidden elements
               CALL W3DL_DIA(WORK(KTHETAF),ISTHETAFY,ISYML,L,ISYMD,D,
     *                       FOCKD,FREQY)
               CALL T3_FORBIDDEN_DL(WORK(KTHETAF),ISYFCKY,ISYMD,D,
     *                              ISYML,L)

C
C               KTHETA(De- f-)_(lmi) = KTHETAF(De- f)_(lmi) * FOCKX
C
               IOPT = 1
               CALL WXDL_V(IOPT,WORK(KTHETAF),ISTHETAFY,FOCKX,ISYFCKX,
     *                     WORK(KTHETA),ISTHETADL,WORK(KEND1),LWRK1)

C
C               KTHETAF(Def- )_(lmi) = KT3 * FOCKY
C
               ! KTHETAF is recycled here
               CALL DZERO(WORK(KTHETAF),NT2SQ(ISTHETAFY))
               IOPT = 1
               CALL WXDL_V(IOPT,WORK(KT3),ISYMT3DL,FOCKY,ISYFCKY,
     *                     WORK(KTHETAF),ISTHETAFY,WORK(KEND1),LWRK1)
               ! Divide it by orbital energy difference and remove the 
               ! forbidden elements
               CALL W3DL_DIA(WORK(KTHETAF),ISTHETAFY,ISYML,L,ISYMD,D,
     *                       FOCKD,FREQY)
               CALL T3_FORBIDDEN_DL(WORK(KTHETAF),ISYFCKY,ISYMD,D,
     *                              ISYML,L)

C
C               KTHETA(De- f-)_(lmi) = KTHETAF(Def- )_(lmi) * FOCKX
C
               IOPT = 3
               CALL WXDL_V(IOPT,WORK(KTHETAF),ISTHETAFY,FOCKX,ISYFCKX,
     *                     WORK(KTHETA),ISTHETADL,WORK(KEND1),LWRK1)

               ! Divide it by orbital energy difference and remove the 
               ! forbidden elements
C
               FREQXY = FREQX + FREQY
C
               CALL W3DL_DIA(WORK(KTHETA),ISTHETADL,ISYML,L,ISYMD,D,
     *                       FOCKD,FREQXY)
               CALL T3_FORBIDDEN_DL(WORK(KTHETA),ISTHETA,ISYMD,D,
     *                              ISYML,L)

            !4th line in Eq. (62) (Dij) (second cont)
             TRANSPOSEW = .TRUE.
C          ---------------------------------------------
C          3rd line of Eq. (62)
C          ---------------------------------------------
 
C
C               KTHETAF(Def- )_(lmi) = KT3 * FOCKX
C
               ! KTHETAF is recycled here
               CALL DZERO(WORK(KTHETAF),NT2SQ(ISTHETAFX))
               IOPT = 1
               CALL WXDL_V(IOPT,WORK(KT3),ISYMT3DL,FOCKX,ISYFCKX,
     *                     WORK(KTHETAF),ISTHETAFX,WORK(KEND1),LWRK1)
               ! Divide it by orbital energy difference and remove the 
               ! forbidden elements
               CALL W3DL_DIA(WORK(KTHETAF),ISTHETAFX,ISYML,L,ISYMD,D,
     *                       FOCKD,FREQX)
               CALL T3_FORBIDDEN_DL(WORK(KTHETAF),ISYFCKX,ISYMD,D,
     *                              ISYML,L)
 
C
C               KWZU(Def-- )_(lmi) = KTHETAF(Def- )_(lmi) * FOCKY
C
               ! KTHETA is recycled here
               CALL DZERO(WORK(KWZU),NT2SQ(ISTHETADL))
               IOPT = 1
               CALL WXDL_V(IOPT,WORK(KTHETAF),ISTHETAFX,FOCKY,ISYFCKY,
     *                     WORK(KWZU),ISTHETADL,WORK(KEND1),LWRK1)

C
C              Include P(XY) permutation
C

C
C               KTHETAF(Def- )_(lmi) = KT3 * FOCKY
C
               ! KTHETAF is recycled here
               CALL DZERO(WORK(KTHETAF),NT2SQ(ISTHETAFY))
               IOPT = 1
               CALL WXDL_V(IOPT,WORK(KT3),ISYMT3DL,FOCKY,ISYFCKY,
     *                     WORK(KTHETAF),ISTHETAFY,WORK(KEND1),LWRK1)
               ! Divide it by orbital energy difference and remove the 
               ! forbidden elements
               CALL W3DL_DIA(WORK(KTHETAF),ISTHETAFY,ISYML,L,ISYMD,D,
     *                       FOCKD,FREQY)
               CALL T3_FORBIDDEN_DL(WORK(KTHETAF),ISYFCKY,ISYMD,D,
     *                              ISYML,L)

C
C               KWZU(Def-- )_(lmi) = KTHETAF(Def- )_(lmi) * FOCKX
C
               IOPT = 1
               CALL WXDL_V(IOPT,WORK(KTHETAF),ISTHETAFY,FOCKX,ISYFCKX,
     *                     WORK(KWZU),ISTHETADL,WORK(KEND1),LWRK1)


               ! Divide it by orbital energy difference and remove the 
               ! forbidden elements
C
               FREQXY = FREQX + FREQY
C
               CALL W3DL_DIA(WORK(KWZU),ISTHETADL,ISYML,L,ISYMD,D,
     *                       FOCKD,FREQXY)
               CALL T3_FORBIDDEN_DL(WORK(KWZU),ISTHETA,ISYMD,D,
     *                              ISYML,L)

               ! add KWZU(Def-- )_(lmi) + KTHETA(De- f-)_(lmi)
               CALL DAXPY(NT2SQ(ISTHETADL),ONE,WORK(KWZU),1,
     *                    WORK(KTHETA),1)
C           ------------------------------------------------
C           Read WMAT_bar from file 
C           ------------------------------------------------
C
            IADR = ISWTL(ISWMATDL,ISYML) + NT2SQ(ISWMATDL)*(L-1) + 1
            CALL GETWA2(LUWBMAT,FNWBMAT,WORK(KWMAT),IADR,
     *                  NT2SQ(ISWMATDL))

C
C               contract... (3rd line of Eq. (62))
C

            TRANSPOSEW = .FALSE.
            CALL DIJ_CONT_CUB(TRANSPOSEW,DIJ,ONE,WORK(KWMAT),ISWMATDL,
     *                        WORK(KTHETA),ISTHETADL,WORK(KEND1),LWRK1)
                

C
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

C           ----------------------------------------------------
C           Read wZU^{Deb-}_{l- m- n-} from file...
C           ----------------------------------------------------
C
               IADR = ISWTL(ISTHETADL,ISYML) + NT2SQ(ISTHETADL)*(L-1) 
     *              + 1
               CALL GETWA2(LUWZU,FNWZU,WORK(KWZU),IADR,
     *                     NT2SQ(ISTHETADL))
             !second contribution to Dab (second line in Eq. (61))
             TRANSPOSEW = .TRUE.

             CALL DAB_CONT_CUB_t0(TRANSPOSEW,DAB,WORK(KWMAT),ISWMATDL,
     *                         WORK(KWZU),ISTHETADL,WORK(KEND1),LWRK1)

            !4th line in Eq. (62) (Dij)
             TRANSPOSEW = .TRUE.
             CALL DIJ_CONT_CUB(TRANSPOSEW,DIJ,HALF,WORK(KWMAT),ISWMATDL,
     *                         WORK(KWZU),ISTHETADL,WORK(KEND1),LWRK1)


C
C           -----------------------------------------------------------
C           ...and create wZU^{Deb-}_{l- m- n-} + wZU^{Dbe-}_{l- n- m-}
C           -----------------------------------------------------------

               CALL CC_T2MOD(WORK(KWZU),ISTHETADL,ONE)

C            ---------------------------------------------------------
C            Get THETA + wZU^{Deb-}_{l- m- n-} + wZU^{Dbe-}_{l- n- m-}
C            ---------------------------------------------------------

               CALL DAXPY(NT2SQ(ISTHETADL),ONE,WORK(KWZU),1,
     *                       WORK(KTHETA),1)

            END IF
C


            !generate WMAT-tilde:
            CALL CC_T2MOD(WORK(KWMAT),ISWMATDL,HALF)
C
C            -------------------------------------------------------
C            D(fb) <- D(fb)+ sum_em Wtilde_bar^DL(em,fN) T3^DL(em,bN):
             ! FOR QUADRU = .TRUE. T3^DL(em,bN) becomes THETA_Z^DL(em,bN)
C            -------------------------------------------------------
             TRANSPOSEW = .FALSE.


             CALL DAB_CONT_CUB(TRANSPOSEW,DAB,WORK(KWMAT),ISWMATDL,
     *                         WORK(KTHETA),ISTHETADL,WORK(KEND1),LWRK1)
C
            IF (QUADR) THEN
               ! Construct THETA^DL(em,f- i)
               CALL DZERO(WORK(KTHETAF),NT2SQ(ISTHETAFX))
               IOPT = 1
               CALL WXDL_V(IOPT,WORK(KT3),ISYMT3DL,FOCKX,ISYFCKX,
     *                 WORK(KTHETAF),ISTHETAFX,WORK(KEND1),LWRK1)
C
               CALL W3DL_DIA(WORK(KTHETAF),ISTHETAFX,ISYML,L,ISYMD,D,
     *                       FOCKD,FREQX)
               CALL T3_FORBIDDEN_DL(WORK(KTHETAF),ISTHETA,ISYMD,D,
     *                              ISYML,L)
C
            END IF

            TRANSPOSEW = .FALSE.
            CALL DIJ_CONT_CUB(TRANSPOSEW,DIJ,ONE,WORK(KWMAT),ISWMATDL,
     *                        WORK(KTHETA),ISTHETADL,WORK(KEND1),LWRK1)
C
            IF (QUADR) THEN
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
       !Now it  becomes the second line of Eq (62) (part of it).
c



C
C               (1) KTHETAF(X) = KT3 * FOCKX
C
               ! KTHETAF is recycled here
               CALL DZERO(WORK(KTHETAF),NT2SQ(ISTHETAFX))
               IOPT = 2
               CALL WXDL_V(IOPT,WORK(KT3),ISYMT3DL,FOCKX,ISYFCKX,
     *                     WORK(KTHETAF),ISTHETAFX,WORK(KEND1),LWRK1)
               ! Divide it by orbital energy difference and remove the 
               ! forbidden elements
               CALL W3DL_DIA(WORK(KTHETAF),ISTHETAFX,ISYML,L,ISYMD,D,
     *                       FOCKD,FREQX)
               CALL T3_FORBIDDEN_DL(WORK(KTHETAF),ISYFCKX,ISYMD,D,
     *                              ISYML,L)
C
C               (2) KTHETA = KTHETAF(X) * FOCKY
C
               !KTHETA is reused here
               CALL DZERO(WORK(KTHETA),NT2SQ(ISTHETADL))
C
               IOPT = 2
               CALL WXDL_V(IOPT,WORK(KTHETAF),ISTHETAFX,FOCKY,ISYFCKY,
     *                     WORK(KTHETA),ISTHETADL,WORK(KEND1),LWRK1)

C
C               Apply P(XY) permutation
C

C
C               (3) KTHETAF(Y) = KT3 * FOCKY
C
               ! KTHETAF is recycled here
               CALL DZERO(WORK(KTHETAF),NT2SQ(ISTHETAFY))
               IOPT = 2
               CALL WXDL_V(IOPT,WORK(KT3),ISYMT3DL,FOCKY,ISYFCKY,
     *                     WORK(KTHETAF),ISTHETAFY,WORK(KEND1),LWRK1)
               ! Divide it by orbital energy difference and remove the 
               ! forbidden elements
               CALL W3DL_DIA(WORK(KTHETAF),ISTHETAFY,ISYML,L,ISYMD,D,
     *                       FOCKD,FREQY)
               CALL T3_FORBIDDEN_DL(WORK(KTHETAF),ISYFCKY,ISYMD,D,
     *                              ISYML,L)
C
C               (4) KTHETA = KTHETAF(Y) * FOCKX
C
               IOPT = 2
               CALL WXDL_V(IOPT,WORK(KTHETAF),ISTHETAFY,FOCKX,ISYFCKX,
     *                     WORK(KTHETA),ISTHETADL,WORK(KEND1),LWRK1)


               ! Divide it by orbital energy difference and remove the 
               ! forbidden elements
C
               FREQXY = FREQX + FREQY
C
               CALL W3DL_DIA(WORK(KTHETA),ISTHETADL,ISYML,L,ISYMD,D,
     *                       FOCKD,FREQXY)
               CALL T3_FORBIDDEN_DL(WORK(KTHETA),ISTHETA,ISYMD,D,
     *                              ISYML,L)
C
c      Now we construct wXY(Def-)_(i- m- l-)

               CALL DZERO(WORK(KWZU),NT2SQ(ISTHETADL))
               CALL DZERO(WORK(KTHETAF),NT2SQ(ISTHETADL))
C
               CALL READ_T3_ALBJ(LUWZU,FNWZU,ISTHETA,WORK(KWZU),
     *                           ISTHETADL,L,ISYML,ISYMD)
c
               !transpose and accumalte and accumulate
               CALL TRANS_AIBJ_BJAI(WORK(KWZU),WORK(KTHETAF),ISTHETADL)
               CALL DAXPY(NT2SQ(ISTHETADL),ONE,WORK(KTHETAF),1,
     *                    WORK(KTHETA),1)
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
                        ISYMFI = MULD2H(ISTHETADL,ISYMEM)
                        ISYMF  = MULD2H(ISYEMF,ISYMEM)
                        ISYMI  = MULD2H(ISYMFI,ISYMF)
                        ISYMFJ = MULD2H(ISYMF,ISYMJ)

C
                        KOFF1 = KTHETA+ IT2SQ(ISYMEM,ISYMFI)
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
      CALL QEXIT('DENABIJC')
C
      RETURN
      END
C  /* Deck dab_cont_cub_t0 */
      SUBROUTINE DAB_CONT_CUB_t0(TRANSPOSEW,DAB,WBARDL,ISWMATDL,THETADL,
     *                        ISTHETADL,WORK,LWORK)
*
**********************************************************************
*
* Calculate the contribution to the DAB density (cubic response) of 
* following type:
*
* Wbar^{Da}(emnl) * theta^{Deb}_{lmn}.
*
* The multiplication is carried out for fixed DL:
*
*
* IF (.NOT. TRANSPOSEW) THEN
*
*    D(ab) = D(ab) + WBARDL(em,aN) * THETADL(em,bN)
*
* ELSE
*
*    D(ab) = D(ab) + WBARDL(aN,em) * THETADL(em,bN)
*
* END IF
*
* Filip Pawlowski, 05-Sep-2003, Aarhus.
**********************************************************************
*
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      LOGICAL TRANSPOSEW
C
      INTEGER ISWMATDL,ISTHETADL
      INTEGER LWORK
      INTEGER ISYMN,ISYEMA,ISYEMB,ISYMEM,ISYMB,ISYMA,ISYMAN,ISYMBN
      INTEGER KAN,KOFF1,KBN,KOFF2,KOFF3,NTOTEM,NTOTA
      INTEGER ISYMM,ISYANE,ISYBNE,ISYME,KEM,NTOTB
      INTEGER NTOTAN
      INTEGER KWBARTR,KEND1,LWRK1
C
#if defined (SYS_CRAY)
      REAL DAB(*),WBARDL(*),THETADL(*)
      REAL WORK(LWORK)
      REAL ONE,HALF
#else
      DOUBLE PRECISION DAB(*),WBARDL(*),THETADL(*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION ONE,HALF
#endif
C
      PARAMETER (ONE = 1.0D0, HALF = 0.5D0)
C
      CALL QENTER('DABCUB')
C

      IF (TRANSPOSEW) THEN 
      !transpose Wbar^DL(em,an) to Wbar^DL(an,em)
         KWBARTR = 1
         KEND1   = KWBARTR + NT2SQ(ISWMATDL)
         LWRK1   = LWORK   - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*)'Memory available: ', LWORK
            WRITE(LUPRI,*)'Memory needed   : ', KEND1
            CALL QUIT('Insufficient memory in DAB_CONT_CUB ')
         END IF
C
         CALL TRANS_AIBJ_BJAI(WBARDL,WORK(KWBARTR),ISWMATDL)
C
      END IF

C--------------------------------------------------------------
C     Calculate D(ab) = D(ab) + WBARDL(em,aN) * THETADL(em,bN)
C--------------------------------------------------------------
C
C     -----------------------------------
C     Loop over outermost occupied index:
C     -----------------------------------
C
      DO ISYMN = 1, NSYM
         ISYEMA = MULD2H(ISWMATDL,ISYMN)
         ISYEMB = MULD2H(ISTHETADL,ISYMN)
C
         DO N = 1, NRHF(ISYMN)
C
C           -------------------------------------------------------
C           D(ab) <- D(ab)+ sum_em Wbar^DL(em,aN) THETA^DL(em,bN):
C           -------------------------------------------------------
            DO ISYMEM = 1, NSYM
               ISYMB  = MULD2H(ISYEMB,ISYMEM)
               ISYMA  = MULD2H(ISYEMA,ISYMEM)
               ISYMAN = MULD2H(ISYMA,ISYMN)
               ISYMBN = MULD2H(ISYMB,ISYMN)

               KAN    = IT1AM(ISYMA,ISYMN)+NVIR(ISYMA)*(N-1)+1
C
               IF (.NOT.TRANSPOSEW) THEN
                  KOFF1  = 1       + IT2SQ(ISYMEM,ISYMAN)
     *                             + NT1AM(ISYMEM)*(KAN-1)
               ELSE
                  KOFF1  = KWBARTR + IT2SQ(ISYMEM,ISYMAN)
     *                             + NT1AM(ISYMEM)*(KAN-1)
               END IF
C               
               KBN    = IT1AM(ISYMB,ISYMN)+NVIR(ISYMB)*(N-1)+1
               KOFF2  = 1     + IT2SQ(ISYMEM,ISYMBN)
     *                        + NT1AM(ISYMEM)*(KBN-1)
C
               KOFF3  = IMATAB(ISYMA,ISYMB) + 1
C
               NTOTEM = MAX(NT1AM(ISYMEM),1)
               NTOTA  = MAX(NVIR(ISYMA),1)

               IF (.NOT.TRANSPOSEW) THEN
                  CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),
     *                       NT1AM(ISYMEM),half,WBARDL(KOFF1),NTOTEM,
     *                       THETADL(KOFF2),NTOTEM,ONE,DAB(KOFF3),
     *                       NTOTA)
               ELSE
                  CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),
     *                       NT1AM(ISYMEM),half,WORK(KOFF1),NTOTEM,
     *                       THETADL(KOFF2),NTOTEM,ONE,DAB(KOFF3),
     *                       NTOTA)
               END IF


            END DO ! ISYMEM
C
         END DO ! N
      END DO    ! ISYMN
C
      CALL QEXIT('DABCUB')
C
      RETURN
      END
C  /* Deck aden_dij_jk_t0 */
      SUBROUTINE ADEN_DIJ_JK_t0(DIJ,THLM,ISYMTHLM,WLM,ISYMWLM)
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
     *              ONE,THLM(KOFF1),NTOTDEF,WLM(KOFF2),NTOTDEF,
     *              ONE,DIJ(KOFF3),NTOTI)
C
      END DO
C
      CALL QEXIT('ADEN_DIJ_JK')
C
      RETURN
      END
C  /* Deck aden_dab_lm_cub_t0 */
      SUBROUTINE ADEN_DAB_LM_CUB_t0(IOPT,DAB,THLM,ISYMTHLM,WLM,
     *                           ISYMWLM,
     *                           WORK,LWORK)
C
      IMPLICIT NONE
#include <priunit.h>
#include <dummy.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER IOPT
      INTEGER ISYMTHLM,ISYMWLM,LWORK
      INTEGER ISYMN,ISYMDEB,ISYMDEA,ISYMB,ISYMDE,ISYMA
      INTEGER KOFF1,KOFF2,KOFF3
      INTEGER NTOTDE,NTOTA
      INTEGER KWDAEN,KTHDBEN,KEND1,LWRK1
      INTEGER ISYMDBE,ISYMDAE,ISYME,ISYMDB,ISYMDA,ISYMEN,ISYMD
      INTEGER NTOTD
      INTEGER ISYMEDN,NTOTB
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
      CALL QENTER('DABLMCB')
C
      IF (IOPT .GT. 3) 
     *   CALL QUIT('Wrong IOPT value in ADEN_DAB_LM_CUB')
C
      IF ((IOPT .EQ. 1) .OR. (IOPT .EQ. 2)) THEN

C
C        D(ab) = W^LM(dean) * THETA^LM(debn)
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
     *                  + NMAABC(ISYMDEA)*(N-1)
     *                  + IMAABC(ISYMDE,ISYMA)
     *                  + 1
                  KOFF2 = IMAABCI(ISYMDEB,ISYMN) 
     *                  + NMAABC(ISYMDEB)*(N-1)
     *                  + IMAABC(ISYMDE,ISYMB)
     *                  + 1
                  KOFF3 = IMATAB(ISYMA,ISYMB) + 1
C
                  NTOTDE = MAX(NMATAB(ISYMDE),1)
                  NTOTA   = MAX(NVIR(ISYMA),1)
C
                  CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),
     *                       NMATAB(ISYMDE),-one,WLM(KOFF1),NTOTDE,
     *                       THLM(KOFF2),NTOTDE,ONE,DAB(KOFF3),NTOTA)
C
               END DO   ! N
            END DO      ! ISYMB
         END DO         ! ISYMN
C
      END IF
C
      IF ((IOPT .EQ. 2).OR.(IOPT .EQ. 3)) THEN   
C
C        Calculate second contribution to D(ab)
C
         KWDAEN = 1
         KTHDBEN  = KWDAEN + NMAABCI(ISYMWLM)
         KEND1   = KTHDBEN + NMAABCI(ISYMTHLM)
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWRK1
            WRITE(LUPRI,*) 'Memory needed    : ',KEND1
            CALL QUIT('Insufficient space in ADEN_DAB_LM_CUB')
         END IF
C
         CALL DZERO(WORK(KWDAEN),NMAABCI(ISYMWLM))
         CALL DZERO(WORK(KTHDBEN),NMAABCI(ISYMTHLM))
 
         !Sort W^LM(dean) to W^LM(daen)
         CALL FACBI(WORK(KWDAEN),WLM,ISYMWLM)
 
         !Sort THETA^LM(debn) to THETA^LM(dben)
         CALL FACBI(WORK(KTHDBEN),THLM,ISYMTHLM)
C
C        D(ab) = W^LM(daen) * THETA^LM(dben)
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
     *                  + NMAABC(ISYMDEA)*(N-1)
     *                  + IMAABC(ISYMDE,ISYMA)
     *                  + KWDAEN
                  KOFF2 = IMAABCI(ISYMDEB,ISYMN) 
     *                  + NMAABC(ISYMDEB)*(N-1)
     *                  + IMAABC(ISYMDE,ISYMB)
     *                  + KTHDBEN
                  KOFF3 = IMATAB(ISYMA,ISYMB) + 1
C
                  NTOTDE = MAX(NMATAB(ISYMDE),1)
                  NTOTA   = MAX(NVIR(ISYMA),1)
C
                  CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),
     *                       NMATAB(ISYMDE),-ONE,WORK(KOFF1),NTOTDE,
     *                       WORK(KOFF2),NTOTDE,ONE,DAB(KOFF3),NTOTA)
C
               END DO   ! N
            END DO      ! ISYMB
         END DO         ! ISYMN
C
      END IF
C
      IF (IOPT .EQ. 0) THEN
C
         KWDAEN = 1
         KTHDBEN  = KWDAEN + NMAABCI(ISYMWLM)
         KEND1   = KTHDBEN + NMAABCI(ISYMTHLM)
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWRK1
            WRITE(LUPRI,*) 'Memory needed    : ',KEND1
            CALL QUIT('Insufficient space in ADEN_DAB_LM_CUB (2)')
         END IF
C
         CALL DZERO(WORK(KWDAEN),NMAABCI(ISYMWLM))
         CALL DZERO(WORK(KTHDBEN),NMAABCI(ISYMTHLM))
 
         !Sort W^LM(aedn) to W^LM(a,edn)
         CALL FA_BCI(WORK(KWDAEN),WLM,ISYMWLM,2)
 
         !Sort THETA^LM(bedn) to THETA^LM(b,edn)
         CALL FA_BCI(WORK(KTHDBEN),THLM,ISYMTHLM,2)

C
C        D(ab) = W^LM(aedn) * THETA^LM(bedn)
C
         DO ISYMEDN = 1,NSYM
            ISYMB = MULD2H(ISYMTHLM,ISYMEDN)
            ISYMA = MULD2H(ISYMWLM,ISYMEDN)
C
            KOFF1 = IMAAOBCI(ISYMA,ISYMEDN)
     *            + KWDAEN
            KOFF2 = IMAAOBCI(ISYMB,ISYMEDN)
     *            + KTHDBEN
            KOFF3 = IMATAB(ISYMA,ISYMB) + 1
C
            NTOTB = MAX(NVIR(ISYMB),1)
            NTOTA   = MAX(NVIR(ISYMA),1)
C
            CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),
     *               NMAABI(ISYMEDN),-ONE,WORK(KOFF1),NTOTA,
     *               WORK(KOFF2),NTOTB,ONE,DAB(KOFF3),NTOTA)
C
         END DO         ! ISYMEDN
C
      END IF
C
      CALL QEXIT('DABLMCB')
C
      RETURN
      END
