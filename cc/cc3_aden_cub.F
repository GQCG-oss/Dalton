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
C  /* Deck cc3_aden_cub */
      SUBROUTINE CC3_ADEN_CUB(LISTL,IDLSTL,LISTR,IDLSTR,
     *                            XLAMDP0,XLAMDH0,FOCK0,
     *                            DIJ,DAB,DIA,ISYDEN,
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
      REAL WORK(LWORK)
#else
      DOUBLE PRECISION XLAMDP0(*),XLAMDH0(*),FOCK0(*)
      DOUBLE PRECISION DAB(*),DIJ(*),DIA(*)
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
      CALL DZERO(DIA,NT1AM(ISYDEN))
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
      CALL CC3_ADENVIR_CUB(DIJ,DAB,DIA,ISYDEN,LISTL,IDLSTL,LISTR,IDLSTR,
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
         WRITE(LUPRI,*)'DIA density after CC3_ADENVIR_CUB '
         CALL PRINT_MATAI(DIA,ISYDEN)
      END IF
C
      IF (LISTR(1:3).EQ.'R2 ') THEN
         CALL CC3_ADENOCC_CUB(LISTL,IDLSTL,LISTR,IDLSTR,
     *                               WORK(KLAMP0),WORK(KLAMH0),FOCK0,
     *                               DIJ,DAB,DIA,ISYDEN,
     *                               WORK(KEND1),LWRK1,
     *                               LUDELD,FNDELD,LUCKJD,FNCKJD,LUDKBC,
     *                               FNDKBC,LUTOC,FNTOC,LU3VI,FN3VI,
     *                               LUDKBC3,FNDKBC3,LU3FOP,FN3FOP,
     *                               LU3FOPX,FN3FOPX,
     *                               LU3FOP2X,FN3FOP2X)
C
         IF (IPRINT .GT. 55) THEN
            WRITE(LUPRI,*)'DAB density after CC3_ADENOCC_CUB '
            CALL PRINT_MATAB(DAB,ISYDEN)
            WRITE(LUPRI,*)'DIJ density after CC3_ADENOCC_CUB '
            CALL PRINT_MATIJ(DIJ,ISYDEN)
            WRITE(LUPRI,*)'DIA density after CC3_ADENOCC_CUB '
            CALL PRINT_MATAI(DIA,ISYDEN)
         END IF
C
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
C  /* Deck cc3_adenocc_cub */
      SUBROUTINE CC3_ADENOCC_CUB(LISTL,IDLSTL,LISTR,IDLSTR,
     *                            XLAMDP0,XLAMDH0,FOCK0,
     *                            DIJ,DAB,DIA,ISYDEN,
     *                            WORK,LWORK,
     *                            LUDELD,FNDELD,LUCKJD,FNCKJD,LUDKBC,
     *                            FNDKBC,LUTOC,FNTOC,LU3VI,FN3VI,
     *                            LUDKBC3,FNDKBC3,LU3FOP,FN3FOP,
     *                            LU3FOPX,FN3FOPX,
     *                            LU3FOP2X,FN3FOP2X)
*
*************************************************************************
*             
*  Calculate all those contractions for A density which should be 
*  evaluated for 2 fixed occupied indeces:
*
*  1) Contributions to Dab density:
*
*     Dab <-- 1/2 Wbar^LM(naed) theta^{d--e--b--}_{LMn} 
*               + Wbar^LM(nead) theta^{d--b--e--}_{LMn}
*               + thetabar^{d-ae}_{LMn} theta^{dbe-}_{LMn-}
*               + thetabar^{d-ea}_{LMn} theta^{de-b}_{LMn-}
*               + thetabar^{aed-}_{LMn} theta^{be-d}_{LMn-}
*
*  2) Contributions to Dij density:
*
*     Dij <-- 1/2 Wbar^LM(fjed) theta^{d--e--f--}_{LMi}
*               + thetabar^{d-ef}_{LMj} theta^{de-f-}_{LMi}
*
*  3) Contributions to Dai density:
*
*     Dai <-- T2bar^{de}_{LM} ( theta^{d--e--a--}_{LMi} 
*                               - theta^{d--a--e--}_{LMi} )
*
*
*  where single bar "-" denotes single transformation of an index
*  and  double bar "--" denotes double transformation of an index. 
*
*************************************************************************
*
*  Before the contractions are carried out the following intermediates
*  must be evaluated in this routine ("eps" denotes orbital energy
*  difference and "ome_X" denotes a frequency associated with perturnation X):
*
*
*
*  1) Intermediates for the first-order triples multipliers tbarZ:
*
*    1a) Wbar^LM(naed) intermadiate known from quadratic response densities:
*
*    1b) thetabar^{d-ae}_{LMn} = 
*                        - (Z_dc tbar0^{cae}_{LMn}) / (eps^{dae}_{LMn}+omega_Z)
*
*    1c) thetabar^{d-ea}_{LMn} defined in the same way as 1b)
*    
*    1d) thetabar^{aed-}_{LMn} defined in the same way as 1b), but requiring
*        the transformation of the last index in tbar0
*
*-------------------------------------------------------------------------
*
*
*  2) Intermediates for the second-order triples amplitudes tXY:
*
*    2a) theta^{d--e--b--}_{LMn} = theta^{d--eb}_{LMn} + theta^{d-e-b}_{LMn}
*
*    where
*
*      2aa) theta^{d--eb}_{LMn} = 
*             PXY (X_dc theta^{c-eb}_{LMn}) / (eps^{deb}_{LMn} - ome_X - ome_Y)
*
*      2ab) theta^{d-e-b}_{LMn} =
*             PXY (X_dc theta^{ce-b}_{LMn} +  X_ec theta^{d-cb}_{LMn})
*                  / (eps^{deb}_{LMn} - ome_X - ome_Y)
*
*      where
*
*        2aaa) theta^{c-eb}_{LMn} = 
*                  (Y_ca t0^{aeb}_{LMn}) / (eps^{ceb}_{LMn} - ome_Y)
*
*        and
*
*        PXY is the permutation operator (permutes X and Y perturbations)
*
*
*    2b) theta^{dbe-}_{LMn-} =
*          PXY ( X_ec w^{dbc}_{LMn-} + X^{Y}_ec t0^{dbc}_{LMn}
*                - X_jn theta^{dbe-}_{LMj} ) / (eps^{dbe}_{LMn} - ome_X - ome_Y)
*
*    where
*
*        2ba) w^{dbc}_{lmn-} = Wdb(cnlm) - theta^{dbe-}_{lmn}
*
*        2bb) theta^{dbe-}_{LMj} is calculated like in 2aaa)
*
*        2bc) X^{Y} = [X,T1Y]
*
*
*    2c) theta^{be-d}_{LMn-} = 
*          PXY ( X_ec w^{bcd}_{LMn-} - X_jn theta^{be-d}_{LMj} )
*              / (eps^{deb}_{LMn} - ome_X - ome_Y)
*
*    where w^{bcd}_{LMn-} and theta^{be-d}_{LMj} have been defined 
*    in 2ba) and 2bb) respectively.
*                     
*
*************************************************************************
*    Written by Filip Pawlowski, Fall 2003, Aarhus
*************************************************************************
*            
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
      PARAMETER(FN3SRTR  = 'CCSDT_FBMAT1_Z',FNCKJDRZ = 'CCSDT_FBMAT2_Z',
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
      INTEGER KL1,KL2,KFOCKL1,KT1RZ,KT2RZ,KFOCKRZ,KEND1,LWRK1
      INTEGER KXIAJB,KT3BOG1,KT3BOL1,KT3BOG2,KT3BOL2,KT3OG1,KT3OG2
      INTEGER KLAMPL1R,KLAMHL1R,KT30KL
      INTEGER KFOCKL1RCK,KT3VIJG1
      INTEGER ISYMT3B,ISYT3B0KL,ISYW3BXKL
      INTEGER KXGADCK,KXLADCK
      INTEGER KT3B0KL,KW3BXKL,ISYMW3BX
      INTEGER KT3BOG2X,KT3BOL2X,KXGADCKX,KXLADCKX
      INTEGER ISYMTETAZ,ISTETAZKL
      INTEGER KTETAXKL
      INTEGER IDLSTL1R,ISYML1R
      INTEGER ISINT2L1R,KT1L1R
C
      INTEGER IDLSTZU,IDLSTRZ,IDLSTRU,ISYMRU
      INTEGER KFOCKRU,ISYMZU,ISYMTETAU,ISYMTETAZU,ISTETAUKL,ISTETAZUKL
      INTEGER MAXX1
      INTEGER K1,K1X,KABCI
      INTEGER KFCKZUV,KFCKUZV,KLAMDPZ,KLAMDHTMP,KLAMDPU
C
      INTEGER KGBCDK
      INTEGER KT1RU,KT2RU
C
      INTEGER ISINT2RZ,ISINT1RU,ISINT2RU,KT3OG2Z
      INTEGER KT3OG2U,KGBCDKZ,KGBCDKU
      INTEGER KEND2,LWRK2
      INTEGER KEND3,LWRK3
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
      integer kx3am
C
#if defined (SYS_CRAY)
      REAL XLAMDP0(*),XLAMDH0(*),FOCK0(*) 
      REAL DAB(*),DIJ(*),DIA(*)
      REAL WORK(LWORK)
      REAL FREQL1,FREQRZ,FREQL1R,FREQRU,FREQZU
      REAL DDOT,XNORMVAL,ONE
#else
      DOUBLE PRECISION XLAMDP0(*),XLAMDH0(*),FOCK0(*) 
      DOUBLE PRECISION DAB(*),DIJ(*),DIA(*)
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
      KL1     = KEND0
      KL2     = KL1     + NT1AM(ISYML1)
      KFOCKL1 = KL2     + NT2SQ(ISYML1)
      KT1RZ     = KFOCKL1 + N2BST(ISYML1)
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
C      WRITE(LUPRI,*) 'Norm of L2TP (after readeing)',
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

      KT3BOG1   = KEND1
      KT3BOL1   = KT3BOG1   + MAX(NTRAOC(ISINT2L1R),NTRAOC(ISYM0))
      KT3BOG2   = KT3BOL1   + MAX(NTRAOC(ISINT2L1R),NTRAOC(ISYM0))
      KT3BOL2   = KT3BOG2   + NTRAOC(ISYM0)
      KT3OG1    = KT3BOL2   + NTRAOC(ISYM0)
      KT3OG2    = KT3OG1    + MAX(NTRAOC(ISINT2),MAXX1)
      KLAMPL1R   = KT3OG2    + NTRAOC(ISINT2)
      KLAMHL1R   = KLAMPL1R   + NLAMDT
      KEND1     = KLAMHL1R   + NLAMDT
C
      KT3OG2Z    = KEND1
      KEND1      = KT3OG2Z       + NTRAOC(ISINT2RZ)
      LWRK1     = LWORK     - KEND1
C
      KT3OG2U    = KEND1
      KEND1      = KT3OG2U       + NTRAOC(ISINT2RU)
      LWRK1     = LWORK     - KEND1
C
      KFOCKL1RCK  = KEND1
      KT3VIJG1  = KFOCKL1RCK  + NT1AM(ISYFCKL1R)
      KEND1     = KT3VIJG1  + NMAABCI(ISYM0)
      LWRK1     = LWORK     - KEND1
C
      KT3BOG2X   = KEND1
      KT3BOL2X   = KT3BOG2X + NTRAOC(ISINT2L1R)
      KEND1      = KT3BOL2X + NTRAOC(ISINT2L1R)
C
      KXGADCK   = KEND1 
      KXLADCK   = KXGADCK + NMAABCI(ISYM0)
      KEND1     = KXLADCK + NMAABCI(ISYM0)
      LWRK1     = LWORK     - KEND1
C
      KXGADCKX   = KEND1 
      KXLADCKX   = KXGADCKX + NMAABCI(ISINT2L1R)
      KEND1     = KXLADCKX + NMAABCI(ISINT2L1R)
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
      KT1L1R  = KEND1
      KEND1  = KT1L1R + NT1AM(ISYML1R)
      LWRK1   = LWORK  - KEND1
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
C----------------------------------------------------------
C     Get Lambda for right list depended on left LISTL list
C----------------------------------------------------------
C
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
C-----------------------------------------------------------------
C     Construct occupied integrals which are required to calculate    
C     t3bar_Y multipliers                                             
C-----------------------------------------------------------------
C
      CALL INTOCC_T3BARX_JK(LUTOC,FNTOC,ISYMOP,
     *                   WORK(KLAMHL1R),ISYML1R,ISINT2L1R,
     *                   DUMMY,DUMMY,.TRUE.,
     *                   WORK(KT3BOG2X),WORK(KT3BOL2X),
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
      CALL INTVIR_T3B0_JK(2,WORK(KXGADCK),WORK(KXLADCK),ISYM0,XLAMDP0,
     *                    ISYM0,
     *                         LU3VI,FN3VI,LU3FOP,FN3FOP,
     *                         WORK(KEND1),LWRK1)
C
C----------------------------------------------
C     Get virtual integrals for t3b0 multipliers
C----------------------------------------------
C
      CALL INTVIR_T3BX_JK(WORK(KXGADCKX),WORK(KXLADCKX),ISINT2L1R,
     *                    WORK(KLAMPL1R),ISYML1R,
     *                    LU3VI,FN3VI,LU3FOP,FN3FOP,
     *                    WORK(KEND1),LWRK1)
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
 
                  !Divide by orbital energy difference and remove 
                  !forbidden elements
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
C---------------------------------------------------------------------------
C                 Calculate KT3B0KL(d- e- b-)_(LMn) = KT3B0KL(d- e- b)_(LMn)
C                           + KT30KL(deb)_(LMn) * FOCKZ
C---------------------------------------------------------------------------
C
                  !KT3B0KL is used here first time
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
C                 ------------------------------------------
C                 KTETAXKL = KT3B0KL(d- e- b-)_(LMn) * FOCKU
C                 ------------------------------------------
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
C                 ------------------------------------------
C                  KT3B0KL(d- e- b)_(LMn) = KT30KL(deb)_(LMn) * FOCKU
C                 ------------------------------------------
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
C                 ------------------------------------------
C                 KT3B0KL(d- e- b-)_(LMn) = KT3B0KL(d- e- b)_(LMn)
C                                        + KT30KL(deb)_(LMn) * FOCKU
C                 ------------------------------------------
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
C                 ------------------------------------------
C                 KTETAXKL = KT3B0KL(d- e- b-)_(LMn) * FOCKZ
C                 ------------------------------------------
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


C
C-------------------------------------------
C                 Get T3BAR0^KL multipliers
C-------------------------------------------
C
                  !KT3B0KL is reused here
                  CALL DZERO(WORK(KT3B0KL),NMAABCI(ISYT3B0KL))
C
                  CALL GET_T3B0_JK_O(WORK(KT3B0KL),ISYT3B0KL,
     *                           WORK(KL2TP),ISYML0,
     *                           WORK(KT3BOL2),WORK(KT3BOG2),ISYM0,
     *                           ISYML,L,ISYMK,K,
     *                           WORK(KEND2),LWRK2)
                  CALL GET_T3B0_JK_V(WORK(KT3B0KL),ISYT3B0KL,
     *                               WORK(KL2TP),ISYML0,
     *                               WORK(KXGADCK),WORK(KXLADCK),
     *                               ISYM0,ISYML,L,ISYMK,K,
     *                               WORK(KEND2),LWRK2)
C
                  CALL GET_T3B0_JK_L1F(WORK(KT3B0KL),ISYT3B0KL,
     *                            WORK(KL1AM),ISYML0,
     *                            WORK(KXIAJB),ISYM0,
     *                            WORK(KL2TP),ISYML0,
     *                            WORK(KFCKBA),ISYM0,
     *                            ISYML,L,ISYMK,K)

 
                  !Divide by orbital energy difference and remove 
                  !forbidden elements
                  CALL T3JK_DIA(WORK(KT3B0KL),ISYT3B0KL,ISYML,L,ISYMK,K,
     *                          WORK(KFOCKD))
                  CALL T3_FORBIDDEN_JK(WORK(KT3B0KL),ISYMT3B,ISYML,L,
     *                                 ISYMK,K)
 
c               call sum_pt3_jk(work(kt3b0kl),isyml,l,isymk,k,isyt3b0kl,
c    *                          work(kx3am),6)
C
                  IF (IPRINT .GT. 55) THEN
                     XNORMVAL = DDOT(NMAABCI(ISYT3B0KL),WORK(KT3B0KL),1,
     *                               WORK(KT3B0KL),1)
                     WRITE(LUPRI,*)'NORM OF KT3B0KL CC3_ADENOCC_CUB ',
     *                              XNORMVAL
                  END IF

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
     *                            WORK(KEND2),LWRK2)
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
C                 <L2|[H^Y,tau3]|HF>
C
                  CALL WBX_JK_FMAT(WORK(KW3BXKL),ISYW3BXKL,
     *                            WORK(KL2TP),ISYML0,
     *                            WORK(KFOCKL1RCK),ISYFCKL1R,
     *                            WORK(KT3BOL2X),WORK(KT3BOG2X),
     *                            WORK(KXGADCKX),WORK(KXLADCKX),
     *                            ISINT2L1R,
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

c               call sum_pt3_jk(work(kw3bxkl),isyml,l,isymk,k,isyw3bxkl,
c    *                          work(kx3am),4)
C
                  IF (IPRINT .GT. 55) THEN
                    XNORMVAL = DDOT(NMAABCI(ISYW3BXKL),WORK(KW3BXKL),1,
     *                              WORK(KW3BXKL),1)
                    WRITE(LUPRI,*)'NORM OF KW3BXKL IN CC3_ADENOCC_CUB ',
     *                             XNORMVAL
                  END IF


                  !CONTRACTION: 3rd line of Eq. 61:

                  !1/2 Wbar^LM(naed) theta^{d--e--b--}_{LMn} 
                  ! + Wbar^LM(nead) theta^{d--b--e--}_{LMn}
                  !(-- denotes double transformation of an index)
                  IOPT = 2
                  CALL ADEN_DAB_LM_CUB(IOPT,DAB,
     *                             WORK(KTETAXKL),ISTETAZUKL,
     *                             WORK(KW3BXKL),ISYW3BXKL,
     *                             WORK(KEND2),LWRK2)
C
                  !CONTRACTION: last line of Eq. 62 (1st term):

                  !1/2 Wbar^LM(fjed) theta^{d--e--f--}_{LMi}
                  CALL ADEN_DIJ_JK(DIJ,WORK(KTETAXKL),ISTETAZUKL,
     *                             WORK(KW3BXKL),ISYW3BXKL)

C
                  !CONTRACTION: last term of Eq. 63 

                  !T2bar^{de}_{LM} ( theta^{d--e--a--}_{LMi}
                  !                  - theta^{d--a--e--}_{LMi} )
                  CALL ADEN_DAI_LM(DIA,WORK(KL2),ISYML1,
     *                       WORK(KTETAXKL),ISTETAZUKL,
     *                       ISYML,L,ISYMK,K,
     *                       WORK(KEND2),LWRK2)



                 CALL DZERO(WORK(KW3BXKL),NMAABCI(ISTETAZKL))
C
                 CALL WJK_GROUND_OCC(WORK(KW3BXKL),ISTETAZKL,
     *                               WORK(KT2RZ),ISYMRZ,
     *                               WORK(KT3OG2),ISYM0,
     *                               ISYML,L,ISYMK,K,
     *                               WORK(KEND2),LWRK2)
C
                 CALL WJK_GROUND_OCC(WORK(KW3BXKL),ISTETAZKL,
     *                               WORK(KT2TP),ISYM0,
     *                               WORK(KT3OG2Z),ISYMRZ,
     *                               ISYML,L,ISYMK,K,
     *                               WORK(KEND2),LWRK2)
C
                 !allocation !!!
                 KABCI = KEND2
                 KEND3  = KABCI + NMAABCI(ISTETAUKL)
                 LWRK3  = LWORK  - KEND3
C
                 IF (LWRK3 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Memory available : ',LWORK
                  WRITE(LUPRI,*) 'Memory needed    : ',KEND3
                  CALL QUIT('Insufficient space in CC3_ADENOCC_CUB (5)')
                 END IF

                 CALL DZERO(WORK(KABCI),NMAABCI(ISTETAUKL))
C
                 CALL WJK_GROUND_OCC(WORK(KABCI),ISTETAUKL,
     *                               WORK(KT2RU),ISYMRU,
     *                               WORK(KT3OG2),ISYM0,
     *                               ISYML,L,ISYMK,K,
     *                               WORK(KEND3),LWRK3)
C
                 CALL WJK_GROUND_OCC(WORK(KABCI),ISTETAUKL,
     *                               WORK(KT2TP),ISYM0,
     *                               WORK(KT3OG2U),ISYMRU,
     *                               ISYML,L,ISYMK,K,
     *                               WORK(KEND3),LWRK3)
C
                 !the real construction of wJK(abci-)
                 CALL TETAX_JK_I(WORK(KT30KL),ISYT30KL,
     *                           WORK(KFOCKRZ),ISYMRZ,
     *                           WORK(KW3BXKL),ISTETAZKL,
     *                           WORK(KEND3),LWRK3)
C
                 CALL WJK_T2(ONE,L,ISYML,K,ISYMK,WORK(KT2TP),ISYM0,
     *                       WORK(KT2TP),
     *                       ISYM0,
     *                       WORK(KFOCKRZ),ISYMRZ,
     *                       WORK(KW3BXKL),ISTETAZKL,
     *                       WORK(KEND3),LWRK3)

c                call sum_pt3_jk(work(KW3BXKL),isyml,l,isymk,k,ISTETAZKL,
c    *                          work(kx3am),7)
C
                 CALL WJK_GROUND(WORK(KW3BXKL),ISTETAZKL,
     *                           WORK(KT2RZ),ISYMRZ,
     *                           WORK(KGBCDK),ISYM0,
     *                           ISYML,L,ISYMK,K,
     *                           WORK(KEND3),LWRK3)
C

                 CALL WJK_GROUND(WORK(KW3BXKL),ISTETAZKL,
     *                           WORK(KT2TP),ISYM0,
     *                           WORK(KGBCDKZ),ISYMRZ,
     *                           ISYML,L,ISYMK,K,
     *                           WORK(KEND3),LWRK3)
C


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
     *                            WORK(KFOCKRU),ISYMRU,
     *                            WORK(KABCI),ISTETAUKL,
     *                            WORK(KEND3),LWRK3)

                  CALL WJK_T2(ONE,L,ISYML,K,ISYMK,WORK(KT2TP),ISYM0,
     *                        WORK(KT2TP),
     *                        ISYM0,
     *                        WORK(KFOCKRU),ISYMRU,
     *                        WORK(KABCI),ISTETAUKL,
     *                        WORK(KEND3),LWRK3)
C
                  CALL WJK_GROUND(WORK(KABCI),ISTETAUKL,
     *                            WORK(KT2RU),ISYMRU,
     *                            WORK(KGBCDK),ISYM0,
     *                            ISYML,L,ISYMK,K,
     *                            WORK(KEND3),LWRK3)
C
                  CALL WJK_GROUND(WORK(KABCI),ISTETAUKL,
     *                            WORK(KT2TP),ISYM0,
     *                            WORK(KGBCDKU),ISYMRU,
     *                            ISYML,L,ISYMK,K,
     *                            WORK(KEND3),LWRK3)
C

C
                  CALL W3JK_DIA(WORK(KABCI),ISTETAUKL,ISYML,L,ISYMK,K,
     *                           WORK(KFOCKD),FREQRU)
C
                  CALL T3_FORBIDDEN_JK(WORK(KABCI),ISYMTETAU,ISYML,L,
     *                                  ISYMK,K)


                  !allocation !!!
                  K1 = KEND3
                  K1X = K1 + NMAABCI(ISYW3BXKL)
                  KEND3 =  K1X + NMAABCI(ISYW3BXKL)
                  LWRK3 = LWORK  - KEND3
C
                  IF (LWRK3 .LT. 0) THEN
                   WRITE(LUPRI,*) 'Memory available : ',LWORK
                   WRITE(LUPRI,*) 'Memory needed    : ',KEND3
                   CALL QUIT('Insufficient space in CC3_ADENOCC_CUB(6)')
                  END IF

c get extra thetaBAR(d-ea)_(LMn) intermediate (special for cubic)
                  CALL DZERO(WORK(K1),NMAABCI(ISYW3BXKL))

                  IOPT = 1
                  CALL TETAX_JK_BC_CUB(.FALSE.,.TRUE.,
     *                                 IOPT,WORK(KT3B0KL),ISYT3B0KL,
     *                                 WORK(KFOCKL1),ISYML1,
     *                                 WORK(K1),ISYW3BXKL,
     *                                 WORK(KEND3),LWRK3)
C
                  CALL W3JK_DIA(WORK(K1),ISYW3BXKL,ISYML,L,ISYMK,K,
     *                          WORK(KFOCKD),-FREQL1)
                  CALL T3_FORBIDDEN_JK(WORK(K1),ISYMW3BX,ISYML,L,
     *                                ISYMK,K)
C
                  !To conform with real sign of t3b multipliers
                  !(noddy code definition)
                  CALL DSCAL(NMAABCI(ISYW3BXKL),-ONE,WORK(K1),1)

c get extra thetaBAR(aed-)_(LMn) intermediate (special for cubic) 
                  CALL DZERO(WORK(K1X),NMAABCI(ISYW3BXKL))

                  CALL TETAX_JK_A_CUB(.FALSE.,.TRUE.,
     *                                 WORK(KT3B0KL),ISYT3B0KL,
     *                                 WORK(KFOCKL1),ISYML1,
     *                                 WORK(K1X),ISYW3BXKL,
     *                                 WORK(KEND3),LWRK3)
C
                  CALL W3JK_DIA(WORK(K1X),ISYW3BXKL,ISYML,L,ISYMK,K,
     *                          WORK(KFOCKD),-FREQL1)
                  CALL T3_FORBIDDEN_JK(WORK(K1X),ISYMW3BX,ISYML,L,
     *                                ISYMK,K)
C
                  !To conform with real sign of t3b multipliers
                  !(noddy code definition)
                  CALL DSCAL(NMAABCI(ISYW3BXKL),-ONE,WORK(K1X),1)

c                 call sum_pt3_jk(work(KW3BXKL),isyml,l,isymk,k,ISTETAZKL,
c    *                           work(kx3am),7)

                  CALL DZERO(WORK(KTETAXKL),NMAABCI(ISTETAZUKL))
C
                 !1st cont to 57
                 !theta^{dbe-}_{LMn-} <-- U_ec w^{dbc}_{LMn-}
                  CALL TETAX_JK_A(WORK(KW3BXKL),ISTETAZKL,
     *                                 WORK(KFOCKRU),ISYMRU,
     *                                 WORK(KTETAXKL),ISTETAZUKL,
     *                                 WORK(KEND3),LWRK3)
C
                 !2nd cont to 57
                 !theta^{dbe-}_{LMn-} <-- U^{Z}_ec t0^{dbc}_{LMn}
                  CALL TETAX_JK_A(WORK(KT30KL),ISYT30KL,
     *                                 WORK(KFCKUZV),ISYMZU,
     *                                 WORK(KTETAXKL),ISTETAZUKL,
     *                                 WORK(KEND3),LWRK3)
                  !3rd cont to 57
                 !theta^{dbe-}_{LMn-} <-- - U_jn theta^{dbe-}_{LMj}
                  CALL DZERO(WORK(KT3B0KL),NMAABCI(ISTETAZKL))
C
                  ! thetaZ(deb-)_(LMn)
                  IOPT = 3
                  CALL TETAX_JK_A(WORK(KT30KL),ISYT30KL,
     *                                 WORK(KFOCKRZ),ISYMRZ,
     *                                 WORK(KT3B0KL),ISTETAZKL,
     *                                 WORK(KEND3),LWRK3)
C
                  CALL W3JK_DIA(WORK(KT3B0KL),ISTETAZKL,ISYML,L,ISYMK,K,
     *                           WORK(KFOCKD),FREQRZ)
C
                  CALL T3_FORBIDDEN_JK(WORK(KT3B0KL),ISYMTETAZ,ISYML,L,
     *                                  ISYMK,K)
C
                  ! thetaZU(deb-)_(LMn-) = thetaZ(deb- )_(LMk) *FOCKU(k,n)
                  CALL TETAX_JK_I(WORK(KT3B0KL),ISTETAZKL,
     *                            WORK(KFOCKRU),ISYMRU,
     *                            WORK(KTETAXKL),ISTETAZUKL,
     *                            WORK(KEND3),LWRK3)
C
C      INCLUDE P(ZU) permutation
C
                 !1st cont to 57
                 !theta^{dbe-}_{LMn-} <-- Z_ec w^{dbc}_{LMn-}
                  CALL TETAX_JK_A(WORK(KABCI),ISTETAUKL,
     *                                 WORK(KFOCKRZ),ISYMRZ,
     *                                 WORK(KTETAXKL),ISTETAZUKL,
     *                                 WORK(KEND3),LWRK3)

                 !2nd cont to 57
                 !theta^{dbe-}_{LMn-} <-- Z^{U}_ec t0^{dbc}_{LMn}
                  CALL TETAX_JK_A(WORK(KT30KL),ISYT30KL,
     *                                 WORK(KFCKZUV),ISYMZU,
     *                                 WORK(KTETAXKL),ISTETAZUKL,
     *                                 WORK(KEND3),LWRK3)
                  !3rd cont to 57
                 !theta^{dbe-}_{LMn-} <-- - Z_jn theta^{dbe-}_{LMj}
                  CALL DZERO(WORK(KT3B0KL),NMAABCI(ISTETAUKL))
C
                  ! thetaU(deb-)_(LMn)
                  IOPT = 3
                  CALL TETAX_JK_A(WORK(KT30KL),ISYT30KL,
     *                                 WORK(KFOCKRU),ISYMRU,
     *                                 WORK(KT3B0KL),ISTETAUKL,
     *                                 WORK(KEND3),LWRK3)
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
     *                            WORK(KEND3),LWRK3)
                  CALL W3JK_DIA(WORK(KTETAXKL),ISTETAZUKL,ISYML,L,
     *                          ISYMK,K,WORK(KFOCKD),FREQZU)
 
                  CALL T3_FORBIDDEN_JK(WORK(KTETAXKL),ISYMTETAZU,
     *                                 ISYML,L,ISYMK,K)

                  !last line in Eq. 61 (term 1)
                  !thetabar^{d-ae}_{LMn} theta^{dbe-}_{LMn-}
                  IOPT = 3
                  CALL ADEN_DAB_LM_CUB(IOPT,DAB,
     *                             WORK(KTETAXKL),ISTETAZUKL,
     *                             WORK(k1),ISYW3BXKL,
     *                             WORK(KEND3),LWRK3)


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
     *                                 WORK(KEND3),LWRK3)
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
     *                            WORK(KEND3),LWRK3)
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
     *                                 WORK(KEND3),LWRK3)
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
     *                            WORK(KEND3),LWRK3)
C

                  IOPT = 3
                  CALL TETAX_JK_BC_CUB(.TRUE.,.FALSE.,
     *                                 IOPT,WORK(KW3BXKL),ISTETAZKL,
     *                                 WORK(KFOCKRU),ISYMRU,
     *                                 WORK(KTETAXKL),ISTETAZUKL,
     *                                 WORK(KEND3),LWRK3)
C
C      INCLUDE P(ZU) permutation
C
                  IOPT = 3
                  CALL TETAX_JK_BC_CUB(.TRUE.,.FALSE.,
     *                                 IOPT,WORK(KABCI),ISTETAUKL,
     *                                 WORK(KFOCKRZ),ISYMRZ,
     *                                 WORK(KTETAXKL),ISTETAZUKL,
     *                                 WORK(KEND3),LWRK3)

                  CALL W3JK_DIA(WORK(KTETAXKL),ISTETAZUKL,ISYML,L,
     *                          ISYMK,K,WORK(KFOCKD),FREQZU)
 
                  CALL T3_FORBIDDEN_JK(WORK(KTETAXKL),ISYMTETAZU,
     *                                 ISYML,L,ISYMK,K)



                  !last line in Eq. 61 (term 2)
                  !thetabar^{d-ea}_{LMn} theta^{de-b}_{LMn-}
                  IOPT = 1
                  CALL DSCAL(NMAABCI(ISYW3BXKL),2.0D0,WORK(K1),1)
                  CALL ADEN_DAB_LM_CUB(IOPT,DAB,
     *                             WORK(KTETAXKL),ISTETAZUKL,
     *                             WORK(k1),ISYW3BXKL,
     *                             WORK(KEND3),LWRK3)
                  !last line in Eq. 61 (term 3)
                  !thetabar^{aed-}_{LMn} theta^{be-d}_{LMn-}
                  IOPT = 0
                  CALL ADEN_DAB_LM_CUB(IOPT,DAB,
     *                             WORK(KTETAXKL),ISTETAZUKL,
     *                             WORK(k1x),ISYW3BXKL,
     *                             WORK(KEND3),LWRK3)

                  !intermmediates for last line of Eq. 62 (2nd term)
                  CALL DZERO(WORK(KT3B0KL),NMAABCI(ISTETAZKL))
C
                  ! thetaZ(de- b)_(LMn)
                  IOPT = 3
                  CALL TETAX_JK_BC_CUB(.TRUE.,.FALSE.,
     *                                 IOPT,WORK(KT30KL),ISYT30KL,
     *                                 WORK(KFOCKRZ),ISYMRZ,
     *                                 WORK(KT3B0KL),ISTETAZKL,
     *                                 WORK(KEND3),LWRK3)
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
     *                                 WORK(KEND3),LWRK3)
C
                  CALL DZERO(WORK(KT3B0KL),NMAABCI(ISTETAZKL))
C
                  ! thetaZ(deb- )_(LMn)
                  CALL TETAX_JK_A(WORK(KT30KL),ISYT30KL,
     *                                 WORK(KFOCKRZ),ISYMRZ,
     *                                 WORK(KT3B0KL),ISTETAZKL,
     *                                 WORK(KEND3),LWRK3)
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
     *                                 WORK(KEND3),LWRK3)


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
     *                                 WORK(KEND3),LWRK3)
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
     *                                 WORK(KEND3),LWRK3)
C
                  CALL DZERO(WORK(KT3B0KL),NMAABCI(ISTETAUKL))
C
                  ! thetaZ(deb- )_(LMn)
                  CALL TETAX_JK_A(WORK(KT30KL),ISYT30KL,
     *                                 WORK(KFOCKRU),ISYMRU,
     *                                 WORK(KT3B0KL),ISTETAUKL,
     *                                 WORK(KEND3),LWRK3)
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
     *                                 WORK(KEND3),LWRK3)
C
                  CALL W3JK_DIA(WORK(KTETAXKL),ISTETAZUKL,ISYML,L,
     *                          ISYMK,K,WORK(KFOCKD),FREQZU)

                  CALL T3_FORBIDDEN_JK(WORK(KTETAXKL),ISYMTETAZU,
     *                                 ISYML,L,ISYMK,K)



                  !CONTRACTION: last line of Eq. 62 (2nd term)
                  !thetabar^{d-ef}_{LMj} theta^{de-f-}_{LMi}
C
                  CALL ADEN_DIJ_JK(DIJ,WORK(KTETAXKL),ISTETAZUKL,
     *                             WORK(k1),ISYW3BXKL)

C
                  IF (IPRINT .GT. 55) THEN
                    XNORMVAL = DDOT(NMATAB(ISYDEN),DAB,1,DAB,1)
                    WRITE(LUPRI,*)'NORM OF DAB AFTER ADEN_DAB_LM ',
     *                             XNORMVAL
                  END IF
C
                  IF (IPRINT .GT. 55) THEN
                    XNORMVAL = DDOT(NT1AM(ISYDEN),DIA,1,DIA,1)
                    WRITE(LUPRI,*)'NORM OF DIA AFTER ADEN_DAI_LM ',
     *                             XNORMVAL
                  END IF
C
               ENDDO   ! L
            ENDDO      ! ISYML
         ENDDO       ! K
      ENDDO          ! ISYMK 
C
c      write(lupri,*) 'W3BAR in CC3_ADENOCC_CUB'
c      write(lupri,*) 'T30KL in CC3_ADENOCC_CUB'
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
C  /* Deck tetax_jk_a_cub */
      SUBROUTINE TETAX_JK_A_CUB(LAMP,LMUL,T0JK,IST0JK,XOP,ISYMXOP,
     *                          TETAXJK,ISTETAXJK,WORK,LWORK)
C
C TETAXJK(bcai) = TETAXJK(bcai) 
C
C             - xop(ad) t0_jk(bcdi)
C                             
C LAMP = .TRUE. : carry out amplitudes-like transformations
C LMUL = .TRUE. : carry out multipliers-like transformations


      IMPLICIT NONE
C
      LOGICAL LAMP,LMUL
      INTEGER IST0JK, ISYMXOP, ISTETAXJK, LWORK
      INTEGER KAD, KEND1, LWRK1, KOFF1, KOFF2, KOFF3
      INTEGER ISYMI, ISYMBCD, ISYMD, ISYMA, ISYMBCA, ISYMBC
      INTEGER NTOTBC, NTOTA
      INTEGER NTOTD
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
      CALL QENTER('TETACB')
C
      !initial test of logic
      IF (LAMP .EQV. LMUL) THEN
         WRITE(LUPRI,*)'LAMP = ', LAMP
         WRITE(LUPRI,*)'LMUL = ', LMUL
         WRITE(LUPRI,*)'LAMP and LMUL must have opposite values '
         CALL QUIT('Logic fault in TETAX_JK_A_CUB')
      END IF

      KAD  = 1
      KEND1  = KAD + NMATAB(ISYMXOP)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWRK1
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in TETAX_JK_A_CUB')
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
C
               IF (LAMP) THEN
                  KOFF2   = KAD
     *                   + IMATAB(ISYMA,ISYMD)
               ELSE
                  KOFF2   = KAD
     *                   + IMATAB(ISYMD,ISYMA)
               END IF
C
               KOFF3   = 1
     *                + IMAABCI(ISYMBCA,ISYMI)
     *                + NMAABC(ISYMBCA)*(I-1)
     *                + IMAABC(ISYMBC,ISYMA)
C
               NTOTBC = MAX(1,NMATAB(ISYMBC))
               IF (LAMP) THEN
                  NTOTA  = MAX(1,NVIR(ISYMA))
               ELSE 
                  NTOTD  = MAX(1,NVIR(ISYMD))
               END IF
C
C TETAXJK(bcai) = TETAXJK(bcai)  - xop(ad) tb0_jk(bcdi)  
C
               IF (LAMP) THEN
                  CALL DGEMM('N','T',NMATAB(ISYMBC),NVIR(ISYMA),
     *                       NVIR(ISYMD),-ONE,T0JK(KOFF1),NTOTBC,
     *                       WORK(KOFF2),NTOTA,
     *                       ONE,TETAXJK(KOFF3),NTOTBC)
               ELSE
                  CALL DGEMM('N','N',NMATAB(ISYMBC),NVIR(ISYMA),
     *                       NVIR(ISYMD),-ONE,T0JK(KOFF1),NTOTBC,
     *                       WORK(KOFF2),NTOTD,
     *                       ONE,TETAXJK(KOFF3),NTOTBC)

               END IF
C
            END DO
         END DO
      END DO
C                             
      CALL QEXIT('TETACB')
      RETURN
      END
C  /* Deck wjk_ground */
      SUBROUTINE WJK_GROUND(T30JK,ISYT30JK,T2TP,
     *                           ISYMT2,T3VIJG1,
     *                           ISYINT,ISYMJ,J,ISYMK,K,
     *                           WORK,LWORK)

***********************************************************
*    T3VIJG1 : g(ck|bd) sitting as I(bcd,k)
*    
*     T30KL sitting as (bcai)
***********************************************************
C
C     T3X^(abc)_(iJK) = 
C     P(ai,bj,ck) (sum_d t^(ad)_(ij) (ck|bd) ) + 
C    - P(ai,bj,ck) (sum_l t^(ab)_(il) (ck|lj) )
C
C    In this routine we calculate the first contribution in terms of 
C    W intermediate:
C
C    W^JK(bcai) =  W^JK(bcai)
C
C 1)            +  t^ad_ij (ck|bd)
C
C 4)            +  t^ad_ik (bj|cd)
C
C F. Pawlowski, 02-10-2003, Aarhus.
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
      CALL QENTER('WJKGR')
C
C***************************************************
C 1)               t^ad_ij    *   (ck|bd) 
C***************************************************
C
C t2tp(djia) =   I^J(dai) 
C
C (ck|bd) = I(bcd,k) =  I^K(bcd) 
C
C W^JK(bcai) = W^JK(bcai) + I^K(bcd)*I^J(dai) 
C
C symmetry and work allocation
C

      ISYMDAI = MULD2H(ISYMT2,ISYMJ)
      ISYMBCD = MULD2H(ISYINT,ISYMK)
C
      KDAI  = 1
      KEND1  = KDAI  + NMAABI(ISYMDAI)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in WJK_GROUND (1)')
      END IF
C
C  sort t^ad_ij = t2tp(djia) as I^J(dai)
C
      CALL SORT_T2_ABJ(WORK(KDAI),ISYMJ,J,T2TP,ISYMT2)
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
               KOFF1   = 1 + IMAABCI(ISYMBCD,ISYMK)
     *                 + NMAABC(ISYMBCD)*(K-1)
     *                 + IMAABC(ISYMBC,ISYMD) 
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
     *                    NVIR(ISYMD),ONE,T3VIJG1(KOFF1),NTOTBC,
     *                    WORK(KOFF2),NTOTD, 
     *                    ONE,T30JK(KOFF3),NTOTBC) 
            END DO
         END DO
      END DO
C
C****************************************************
C 4)            +  t^ad_ik (bj|cd)
C****************************************************
C
C t2tp(dkia) =   I^K(dai)
C
C (bj|cd) = I(cbdj) =      I^J(cbd)
C
C T^JK(bcai) = T^JK(bcai) + I^J(bcd)*I^K(dai)
C
C symmetry and work allocation
C
      ISYMDAI = MULD2H(ISYMT2,ISYMK)
      ISYMBCD = MULD2H(ISYINT,ISYMJ)
C
      KDAI  = 1
      KEND1  = KDAI  + NMAABI(ISYMDAI)
      LWRK1 = LWORK - KEND1
C
      KCBAI = KEND1
      KEND1 = KCBAI + NMAABCI(ISYT30JK)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in WJK_GROUND (4)')
      END IF
C
      CALL DZERO(WORK(KCBAI),NMAABCI(ISYT30JK))
C
C t2tp(dkia) =   I^K(dai)
C
      CALL SORT_T2_ABJ(WORK(KDAI),ISYMK,K,T2TP,ISYMT2)
C
C T^JK(bcai) = T^JK(bcai) + I^K(cbd)*I^J(dai)
C
      DO ISYMI = 1,NSYM
         ISYMDA = MULD2H(ISYMDAI,ISYMI)
         DO I =  1,NRHF(ISYMI)
            DO ISYMD = 1,NSYM
               ISYMBC  = MULD2H(ISYMBCD,ISYMD)
               ISYMA   = MULD2H(ISYMDA,ISYMD)
               ISYMBCA = MULD2H(ISYMBC,ISYMA)
               KOFF1   = 1 + IMAABCI(ISYMBCD,ISYMJ)
     *                 + NMAABC(ISYMBCD)*(J-1)
     *                 + IMAABC(ISYMBC,ISYMD)
               KOFF2   = KDAI
     *                + IMAABI(ISYMDA,ISYMI)
     *                + NMATAB(ISYMDA)*(I-1)
     *                + IMATAB(ISYMD,ISYMA)
               KOFF3   = KCBAI + IMAABCI(ISYMBCA,ISYMI)
     *                + NMAABC(ISYMBCA)*(I-1)
     *                + IMAABC(ISYMBC,ISYMA)
C 
               NTOTBC = MAX(1,NMATAB(ISYMBC))
               NTOTD  = MAX(1,NVIR(ISYMD))
C
C  add_vir(4)
C
               CALL DGEMM('N','N',NMATAB(ISYMBC),NVIR(ISYMA),
     *                    NVIR(ISYMD),ONE,T3VIJG1(KOFF1),NTOTBC,
     *                    WORK(KOFF2),NTOTD,
     *                    ONE,WORK(KOFF3),NTOTBC)
            END DO
         END DO
      END DO
C
      !put W(cbai) to W(bcai)
      CALL  FBACI(T30JK,WORK(KCBAI),ISYT30JK)
C
      CALL QEXIT('WJKGR')
C
      RETURN   
      END
C  /* Deck tetax_jk_bc_cub */
      SUBROUTINE TETAX_JK_BC_CUB(LAMP,LMUL,IOPT,T0JK,IST0JK,XOP,ISYMXOP,
     *                           TETAXJK,ISTETAXJK,WORK,LWORK)
C
C TETAXJK(bcai) = TETAXJK(bcai) 
C
C             - xop(bd) t0_jk(dcai) (1)
C                             
C             - xop(cd) t0_jk(bdai) (2)
C                             
C IOPT = 1 : calculate only term (1)
C IOPT = 2 : calculate both terms 
C IOPT = 3 : calculate only term (2)
C
C LAMP = .TRUE. : carry out amplitudes-like transformations
C LMUL = .TRUE. : carry out multipliers-like transformations

      IMPLICIT NONE
C
      LOGICAL LAMP,LMUL
      INTEGER IOPT
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
      CALL QENTER('TETJKCB')
C
      !initial test of logic
      IF (LAMP .EQV. LMUL) THEN
         WRITE(LUPRI,*)'LAMP = ', LAMP
         WRITE(LUPRI,*)'LMUL = ', LMUL
         WRITE(LUPRI,*)'LAMP and LMUL must have opposite values '
         CALL QUIT('Logic fault in TETAX_JK_BC_CUB')
      END IF
C
      KAD  = 1
      KEND1  = KAD + NMATAB(ISYMXOP)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWRK1
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in TETAX_JK_BC_CUB')
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

      IF (IOPT .GE. 2) THEN
C
C        TETAXJK(bcai) = TETAXJK(bcai) - xop(cd) t0_jk(bdai)  (term 2)
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
     *                           + IMAABCI(ISYMBDA,ISYMI)
     *                           + NMAABC(ISYMBDA)*(I-1)
     *                           + IMAABC(ISYMBD,ISYMA)
     *                           + NMATAB(ISYMBD)*(A-1)
     *                           + IMATAB(ISYMB,ISYMD)
C
                        IF (LAMP) THEN
                           KOFF2   = KAD
     *                             + IMATAB(ISYMC,ISYMD)
                        ELSE 
                           KOFF2   = KAD
     *                             + IMATAB(ISYMD,ISYMC)
                        END IF 
C
                        KOFF3   = 1
     *                           + IMAABCI(ISYMBCA,ISYMI)
     *                           + NMAABC(ISYMBCA)*(I-1)
     *                           + IMAABC(ISYMBC,ISYMA)
     *                           + NMATAB(ISYMBC)*(A-1)
     *                           + IMATAB(ISYMB,ISYMC)
C
                        NTOTB = MAX(1,NVIR(ISYMB))
C
                        IF (LAMP) THEN
                           NTOTC  = MAX(1,NVIR(ISYMC))
                        ELSE
                           NTOTD  = MAX(1,NVIR(ISYMD))
                        END IF
C
C                       TETAXJK(bcai) = TETAXJK(bcai) - xop(cd) t0_jk(bdai) 
C
                        IF (LAMP) THEN
                           CALL DGEMM('N','T',NVIR(ISYMB),NVIR(ISYMC),
     *                               NVIR(ISYMD),-ONE,T0JK(KOFF1),NTOTB,
     *                               WORK(KOFF2),NTOTC,
     *                               ONE,TETAXJK(KOFF3),NTOTB)
                        ELSE
                           CALL DGEMM('N','N',NVIR(ISYMB),NVIR(ISYMC),
     *                               NVIR(ISYMD),-ONE,T0JK(KOFF1),NTOTB,
     *                               WORK(KOFF2),NTOTD,
     *                               ONE,TETAXJK(KOFF3),NTOTB)
                        END IF
C
                     END DO
                  END DO
               END DO
            END DO
         END DO
C
      END IF
C
      IF (IOPT .LE. 2) THEN

C
C        TETAXJK(bcai) = TETAXJK(bcai) - xop(bd) t0_jk(dcai) (term 1)
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
                        IF (LAMP) THEN
                           KOFF1   = KAD
     *                              + IMATAB(ISYMB,ISYMD)
                         ELSE
                           KOFF1   = KAD
     *                              + IMATAB(ISYMD,ISYMB)
                         END IF
C
                        KOFF2   = 1
     *                           + IMAABCI(ISYMDCA,ISYMI)
     *                           + NMAABC(ISYMDCA)*(I-1)
     *                           + IMAABC(ISYMDC,ISYMA)
     *                           + NMATAB(ISYMDC)*(A-1)
     *                           + IMATAB(ISYMD,ISYMC)
                        KOFF3   = 1
     *                           + IMAABCI(ISYMBCA,ISYMI)
     *                           + NMAABC(ISYMBCA)*(I-1)
     *                           + IMAABC(ISYMBC,ISYMA)
     *                           + NMATAB(ISYMBC)*(A-1)
     *                           + IMATAB(ISYMB,ISYMC)
C        
                        NTOTB = MAX(1,NVIR(ISYMB))
                        NTOTD  = MAX(1,NVIR(ISYMD))
C
C                       TETAXJK(bcai) = TETAXJK(bcai) - xop(bd) t0_jk(dcai) 
C        
                        IF (LAMP) THEN
                           CALL DGEMM('N','N',NVIR(ISYMB),NVIR(ISYMC),
     *                               NVIR(ISYMD),-ONE,WORK(KOFF1),NTOTB,
     *                               T0JK(KOFF2),NTOTD,
     *                               ONE,TETAXJK(KOFF3),NTOTB)
                        ELSE
                           CALL DGEMM('T','N',NVIR(ISYMB),NVIR(ISYMC),
     *                               NVIR(ISYMD),-ONE,WORK(KOFF1),NTOTD,
     *                               T0JK(KOFF2),NTOTD,
     *                               ONE,TETAXJK(KOFF3),NTOTB)
                        END IF
C
                     END DO
                  END DO
               END DO
            END DO
         END DO
C
      END IF
C
      CALL QEXIT('TETJKCB')
      RETURN
      END
C  /* Deck aden_dab_lm_cub */
      SUBROUTINE ADEN_DAB_LM_CUB(IOPT,DAB,THLM,ISYMTHLM,WLM,
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
     *                       NMATAB(ISYMDE),-HALF,WLM(KOFF1),NTOTDE,
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
            WRITE(LUPRI,*) 'Memory available : ',LWORK
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
            WRITE(LUPRI,*) 'Memory available : ',LWORK
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
C  /* Deck cc3_adenvir_cub */
      SUBROUTINE CC3_ADENVIR_CUB(DIJ,DAB,DIA,ISYDEN,LISTL,IDLSTL,LISTR,
     *                   IDLSTR,
     *                   LUTOC,FNTOC,LU3VI,FN3VI,LU3VI2,FN3VI2,
     *                   LUDKBC3,FNDKBC3,
     *                   LU3FOPX,FN3FOPX,LU3FOP2X,FN3FOP2X,
     *                   LU3FOP,FN3FOP,LU3FOP2,FN3FOP2,
     *                   LUCKJD,FNCKJD,LUDKBC,FNDKBC,LUDELD,FNDELD,
     *                   WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Calculate these terms to A density for cubic reponse that
*    should be calculated for 2 fixed virtual indeces.
*
*    (see comments in this routines for the formulas)
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
#include <ccer1rsp.h>

C
      INTEGER ISYM0
      PARAMETER(ISYM0 = 1)
      CHARACTER LISTL0*3, LISTL*3,LISTR*3,LISTL1R*3,LABELL1*8,LABELRZ*8
      CHARACTER LABELRU*8
      CHARACTER LISTRZ*3,LISTRU*3
      CHARACTER*(*) FNTOC, FN3VI, FNDKBC3, FN3FOPX, FN3FOP2X 
      CHARACTER*(*) FNDKBC,FNDELD,FN3VI2,FN3FOP,FN3FOP2,FNCKJD
C
      CHARACTER*10 FNT3, FNWBMAT,FNWBZU
      CHARACTER*14 FNTHETA,FNWZU
      PARAMETER(FNT3 = 'CC3_T3_TMP', FNWBMAT = 'CC3_W3_TMP',  
     *          FNWBZU = 'CC3_WZUTMP',
     *          FNTHETA = 'CC3_THETA3_TMP',FNWZU = 'CC3_WZU____TMP')
C
      CHARACTER*14 FN3SRTR, FNCKJDRZ, FNDELDRZ, FNDKBCRZ
      PARAMETER(FN3SRTR  = 'CCSDT_FBMAT1_Z',FNCKJDRZ = 'CCSDT_FBMAT2_Z',
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
      INTEGER LUT3,LUWBMAT,LUTHETA,LUWZU,LUWBZU
C
      LOGICAL   LOCDBG,LORXL1
      LOGICAL   LORXRZ,LORXRU
      PARAMETER (LOCDBG = .FALSE.)
C
      INTEGER  AIBJCK_PERM
      LOGICAL QUADR
      LOGICAL CUBIC
      LOGICAL T2XNET2Y
      LOGICAL T2XNET2Z,NOVIRT
      LOGICAL LSKIPL1R
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
      INTEGER ISYMD,ISYCKBD0,ISYCKBL1R,ISYCKBDR1Z
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
      INTEGER KT3MAT,KW3BMAT,KW3MATZ,KWTEMP,KS3MAT,KU3MAT,KS3MAT3
      INTEGER KU3MAT3,KT3VBG1,KT3VBG2,KT3VBG3,KT3BVBG1,KT3BVBG2
      INTEGER KT3BVBG3,KSMAT4,KUMAT4,KT3BVBL1,KT3BVBL2,KT3BVBL3
      INTEGER KW3ZVDGZ2
      INTEGER KEND5,LWRK5
      INTEGER LENSQ,LENSQWB,LENSQWZ
      INTEGER ISYML,ISYMDL,ISAIBJ,ISYMJ,ISYMBJ,ISYMAI,ISYAIL
      INTEGER KOFF1,NBJ,IADR
      INTEGER KDAB0,KDIJ0
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
C
      INTEGER IR1TAMP
      INTEGER ILSTSYM
C
      integer kx3am
C
      INTEGER FKW3BXVDG1,FKW3BXVDG2,FKW3BXVDL1,FKW3BXVDL2
      INTEGER FKW3BXVDGX1,FKW3BXVDGX2,FKW3BXVDLX1,FKW3BXVDLX2
      INTEGER ISYCKDL1R

#if defined (SYS_CRAY)
      REAL      FREQL1,FREQL1R,FREQRZ,FREQRU,FREQZU
      REAL      WORK(LWORK)
      REAL      XNORMVAL
      REAL      DAB(*),DIJ(*),DIA(*)
      REAL      DDOT,HALF,ONE
#else
      DOUBLE PRECISION      FREQL1,FREQL1R,FREQRZ,FREQRU,FREQZU
      DOUBLE PRECISION      WORK(LWORK)
      DOUBLE PRECISION      XNORMVAL
      DOUBLE PRECISION      DAB(*),DIJ(*),DIA(*)
      DOUBLE PRECISION      DDOT,HALF,ONE
#endif
C
      PARAMETER(HALF = 0.5D0, ONE = 1.0D0)
C
      CALL QENTER('CC3DENVCB')
C--------------------------------
C     Open temporary files
C--------------------------------
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
      LUCKJDRZU  = -1
      LUDELDRZU  = -1
      LUDKBCRZU  = -1
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
      ELSE IF (LISTL(1:3).EQ.'LE ') THEN
        ISYML1 = ILSTSYM(LISTL,IDLSTL)
        FREQL1 = -EIGVAL(IDLSTL)
        LABELL1 = '- none -'
C
        !we don't have any "right" vector entering a right hand side
        LISTL1R = '---'
        IDLSTL1R = -99
        ISYML1R = IDUMMY
        FREQL1R = DUMMY
C
      ELSE
         CALL QUIT('Unknown left list in CC3_ADENVIR_CUB')
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
      ELSE IF (LISTR(1:3).EQ.'ER1') THEN
        IDLSTZU = IDLSTR
C
         LISTRZ  = 'R1 '
         LABELRZ = lbler1(IDLSTZU)
         ISYMRZ  = isyoer1(IDLSTZU)
         FREQRZ  = frqer1(IDLSTZU)
         LORXRZ  = lorxer1(IDLSTZU)
         IDLSTRZ = IR1TAMP(LABELRZ,LORXRZ,FREQRZ,ISYMRZ)
C
         LISTRU  = 'RE '
         LABELRU = '-- XX --'
         ISYMRU  = isyser1(IDLSTZU)
         FREQRU  = eiger1(IDLSTZU)
         LORXRU  = .FALSE.
         IDLSTRU = ister1(IDLSTZU)
C
      ELSE
       WRITE(LUPRI,*)'LISTR = ',LISTR(1:3)
       WRITE(LUPRI,*)'CC3_ADENVIR_CUB is designed for LISTR = R2 or ER1'
       CALL QUIT('Unknown right list in CC3_ADENVIR_CUB')
      END IF
C
      IF (LORXRZ.OR.LORXRU) THEN
       CALL QUIT('Orbital relaxation not allowed in CC3_ADENVIR_CUB')
      END IF
C
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
      KDAB0 = KEND1
      KDIJ0 = KDAB0 + NMATAB(ISYML1)
      KEND1 = KDIJ0 + NMATIJ(ISYML1)
      LWRK1 = LWORK - KEND1
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
      CALL DZERO(WORK(KDAB0),NMATAB(ISYML1))
      CALL DZERO(WORK(KDIJ0),NMATIJ(ISYML1))
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
      IF (LISTRU(1:3).EQ.'R1 ') THEN
        CALL GET_LAMBDAX(WORK(KLAMDPZ),WORK(KLAMDHZ),LISTRZ,IDLSTRZ,
     *                   ISYMRZ,WORK(KLAMP0),WORK(KLAMH0),WORK(KEND1),
     *                   LWRK1)
        ! get vir-vir block U_(c-,d)
        CALL GET_FOCKX(WORK(KFCKUZV),LABELRU,IDLSTRU,ISYMRU,
     *                    WORK(KLAMDPZ),
     *                    ISYMRZ,WORK(KLAMH0),ISYM0,WORK(KEND1),LWRK1)
        ! get occ-occ block U_(l,k-)
        CALL GET_FOCKX(WORK(KFCKUZO),LABELRU,IDLSTRU,ISYMRU,
     *                    WORK(KLAMP0),
     *                    ISYM0,WORK(KLAMDHZ),ISYMRZ,WORK(KEND1),LWRK1)
      END IF
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
      IF (LISTL(1:3).EQ.'L1 ') THEN
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
      KW3BXOG1   = KEND1
      KW3BXOL1   = KW3BXOG1   + NTRAOC(ISYM0)
      KEND1   = KW3BXOL1   + NTRAOC(ISYM0)
      LWRK1      = LWORK      - KEND1
C
      IF (LISTL(1:3).EQ.'L1 ') THEN
         KFOCKL1RCK    = KEND1
         KW3BXOGX1   = KFOCKL1RCK    + NT1AM(ISYFCKL1R)
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
      ELSE IF (LISTL(1:3).EQ.'LE ') THEN
         KEND2 = KEND1
         LWRK2 = LWRK1
      END IF
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
C
      IF (LISTL(1:3).EQ.'L1 ') THEN
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
      END IF
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
      IF (LISTL(1:3).EQ.'L1 ') THEN
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
      LUT3    = -1
      LUWBMAT = -1
      LUWBZU  = -1
      LUTHETA = -1
      LUWZU   = -1

      CALL WOPEN2(LUT3,FNT3,64,0)
      CALL WOPEN2(LUWBMAT,FNWBMAT,64,0)
      CALL WOPEN2(LUWBZU,FNWBZU,64,0)
      CALL WOPEN2(LUTHETA,FNTHETA,64,0)
      CALL WOPEN2(LUWZU,FNWZU,64,0)
C
C----------------------------
C     Loop over D
C----------------------------
C
      DO ISYMD = 1,NSYM

         ISYCKBD0  = MULD2H(ISYMD,ISYM0)
         ISYCKBDR1Z  = MULD2H(ISYMD,ISINT2RZ)
         ISYCKBDR1U  = MULD2H(ISYMD,ISINT2RU)
         ISYCKBDR1ZU  = MULD2H(ISYMD,ISINT2RZU)
         IF (LISTL(1:3).EQ.'L1 ') THEN
            ISYCKBL1R  = MULD2H(ISYMD,ISYML1R)
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

            IF (LISTL(1:3).EQ.'L1 ') THEN
               KW3BXVDGX1  = KEND3
               KW3BXVDGX2  = KW3BXVDGX1  + NCKATR(ISYCKBL1R)
               KW3BXVDLX1  = KW3BXVDGX2  + NCKATR(ISYCKBL1R)
               KW3BXVDLX2  = KW3BXVDLX1  + NCKATR(ISYCKBL1R)
               KEND3     = KW3BXVDLX2  + NCKATR(ISYCKBL1R)
               LWRK3     = LWORK     - KEND3
            END IF
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
            IF (LISTL(1:3).EQ.'L1 ') THEN
               MAXX1 = MAX(NCKA(ISYCKBD0),NCKA(ISYCKBL1R))
            ELSE IF (LISTL(1:3).EQ.'LE ') THEN
               MAXX1 = NCKA(ISYCKBD0)
            END IF
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
            IF (LISTL(1:3).EQ.'L1 ') THEN
              LSKIPL1R = .FALSE.
              CALL INTVIR_T3BARX_D(LSKIPL1R,
     *                             ISYMOP,LU3VI,FN3VI,LU3VI2,FN3VI2,
     *                             LU3FOP,FN3FOP,LU3FOP2,FN3FOP2,
     *                             WORK(KW3BXVDGX1),WORK(KW3BXVDG1),
     *                             WORK(KW3BXVDGX2),WORK(KW3BXVDG2),
     *                             WORK(KW3BXVDLX1),WORK(KW3BXVDL1),
     *                             WORK(KW3BXVDLX2),WORK(KW3BXVDL2),
     *                             WORK(KLAMPL1R),ISYML1R,WORK(KLAMP0),
     *                             ISYM0,ISYMD,D,WORK(KEND4),LWRK4)
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
               KW3BMAT      = KT3MAT    + NCKIJ(ISCKIJ)
               KW3MATZ      = KW3BMAT     + NCKIJ(ISWBMAT)
c
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
c
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

               IF (LISTL(1:3).EQ.'L1 ') THEN
                  ISYCKDL1R  = MULD2H(ISYMB,ISYML1R)
                  FKW3BXVDGX1  = KEND4
                  FKW3BXVDGX2  = FKW3BXVDGX1  + NCKATR(ISYCKDL1R)
                  FKW3BXVDLX1  = FKW3BXVDGX2  + NCKATR(ISYCKDL1R)
                  FKW3BXVDLX2  = FKW3BXVDLX1  + NCKATR(ISYCKDL1R)
                  KEND4     = FKW3BXVDLX2  + NCKATR(ISYCKDL1R)
                  LWRK4     = LWORK     - KEND4
               END IF
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

               CALL CC3_INDEX(WORK(KINDEX),ISYALJB0)
               CALL CC3_INDEX(WORK(KINDEX2),ISYALJD0)
               CALL CC3_INDEX(WORK(KINDEXBL1),ISYALJBL1)
               CALL CC3_INDEX(WORK(KINDEXDL1),ISYALJDL1)

               DO B = 1,NVIR(ISYMB)
                  CALL DZERO(WORK(KW3BMAT),NCKIJ(ISWBMAT))
                  CALL DZERO(WORK(KW3MATZ),NCKIJ(ISWMATZ))
                  CALL DZERO(WORK(KW3MATU),NCKIJ(ISWMATU))
                  CALL DZERO(WORK(KWMATZU),NCKIJ(ISWMATZU))
                  CALL DZERO(WORK(KWMATZUD),NCKIJ(ISWMATZU))
C
                  IF (LISTL(1:3).EQ.'L1 ') THEN
                    LSKIPL1R = .FALSE.
                    CALL INTVIR_T3BARX_D(LSKIPL1R,
     *                             ISYMOP,LU3VI,FN3VI,LU3VI2,FN3VI2,
     *                             LU3FOP,FN3FOP,LU3FOP2,FN3FOP2,
     *                             WORK(FKW3BXVDGX1),WORK(FKW3BXVDG1),
     *                             WORK(FKW3BXVDGX2),WORK(FKW3BXVDG2),
     *                             WORK(FKW3BXVDLX1),WORK(FKW3BXVDL1),
     *                             WORK(FKW3BXVDLX2),WORK(FKW3BXVDL2),
     *                             WORK(KLAMPL1R),ISYML1R,WORK(KLAMP0),
     *                             ISYM0,ISYMB,B,WORK(KEND5),LWRK5)
                  ELSE IF (LISTL(1:3).EQ.'LE ') THEN
                    LSKIPL1R = .TRUE.
                    CALL INTVIR_T3BARX_d(LSKIPL1R,
     *                             ISYMOP,LU3VI,FN3VI,LU3VI2,FN3VI2,
     *                             LU3FOP,FN3FOP,LU3FOP2,FN3FOP2,
     *                             DUMMY,WORK(FKW3BXVDG1),
     *                             DUMMY,WORK(FKW3BXVDG2),
     *                             DUMMY,WORK(FKW3BXVDL1),
     *                             DUMMY,WORK(FKW3BXVDL2),
     *                             DUMMY,IDUMMY,WORK(KLAMP0),
     *                             ISYM0,ISYMB,B,WORK(KEND5),LWRK5)
                  END IF
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
c
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
C
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
c
c       call sum_pt3(work(KTMAT),isymb,b,isymd,d,
c    *             ISYM0,work(kx3am),4)
C
                  IF (LISTL(1:3).EQ.'L1 ') THEN
                    !<L3|[Y^,tau3]|HF> (virt. part)
                    CALL WBARBD_V(WORK(KTMAT),ISCKIJ,
     *                          WORK(KFOCKL1),
     *                                 ISYML1,WORK(KW3BMAT),ISWBMAT,
     *                                 WORK(KEND5),LWRK5)
C
                    !<L3|[Y^,tau3]|HF> (occ. part)
                    CALL WX_BD_O(1,.FALSE.,.TRUE.,WORK(KTMAT),ISCKIJ,
     *                          WORK(KFOCKL1),
     *                                 ISYML1,WORK(KW3BMAT),ISWBMAT,
     *                                 WORK(KEND5),LWRK5)
 
                    ! <L2|[Y,tau3]|HF>
                    CALL WBARXBD_T2(1,B,ISYMB,D,ISYMD,WORK(KL2TP),ISYM0,
     *                             WORK(KFOCKL1),
     *                 ISYML1,WORK(KW3BMAT),ISWBMAT)
 
C
                    !<L2|[H^Y,tau3]|HF>
                    CALL WBARXBD_TMAT(1,
     *                    WORK(KL2TP),ISYM0,WORK(KW3BMAT),WORK(KWTEMP),
     *                    ISWBMAT,WORK(KFOCKL1RCK),ISYFCKL1R,
     *                    WORK(KW3BXVDLX2),WORK(KW3BXVDLX1),
     *                    WORK(KW3BXVDGX2),
     *                    WORK(KW3BXVDGX1),WORK(KW3BXOLX1),
     *                    WORK(KW3BXOGX1),ISINT2L1R,
     *                    WORK(KEND5),LWRK5,
     *                    WORK(KINDEX),WORK(KINDEX2),
     *                    WORK(KINDSQWB),LENSQWB,
     *                    ISYMB,B,ISYMD,D)
                  END IF
C
                  !<L2Y|[H^,tau3]|HF>
                  CALL WBARXBD_TMAT(1,
     *                  WORK(KL2L1),ISYML1,WORK(KW3BMAT),WORK(KWTEMP),
     *                  ISWBMAT,WORK(KFOCK0CK),ISYM0,
     *                  WORK(KW3BXVDL2),WORK(KW3BXVDL1),
     *                  WORK(KW3BXVDG2),WORK(KW3BXVDG1),
     *                  WORK(KW3BXOL1),WORK(KW3BXOG1),
     *                  ISINT2,
     *                  WORK(KEND5),LWRK5,WORK(KINDEXBL1),
     *                  WORK(KINDEXDL1),WORK(KINDSQWB),LENSQWB,
     *                  ISYMB,B,ISYMD,D)
C
                  !<L1Y|[H^,tau3]|HF>
                  CALL WBARXBD_L1(1,WORK(KL1L1),ISYML1,WORK(KWTEMP),
     *                           WORK(KXIAJB),
     *                ISINT1,WORK(KW3BMAT),WORK(KEND5),LWRK5,
     *                WORK(KINDSQWB),LENSQWB,ISYMB,B,ISYMD,D)
C
                  CALL WBD_DIA(B,ISYMB,D,ISYMD,-FREQL1,ISWBMAT,
     *                         WORK(KW3BMAT),WORK(KDIAGWB),WORK(KFOCKD))
                  CALL T3_FORBIDDEN(WORK(KW3BMAT),ISYML1,ISYMB,B,
     *                              ISYMD,D)
C
                  !To conform with real sign of t3b multipliers
                  !(noddy code definition)
                  CALL DSCAL(NCKIJ(ISWBMAT),-ONE,WORK(KW3BMAT),1)
c
c       call sum_pt3(work(KW3BMAT),isymb,b,isymd,d,
c    *             ISWBMAT,work(kx3am),4)
c
C
C--------------------------------------------------------
C                 Write WBMAT as WBMAT^D(ai,bj,l) to disc
C--------------------------------------------------------
                  CALL WRITE_T3_DL(LUWBMAT,FNWBMAT,WORK(KW3BMAT),ISYML1,
     *                             ISYMD,ISYMB,B)

C
C--------------------------------------------------------
C                 Get T2ZU T20 contribution to DIA density 
C                 (comes from tX * A{Y} * tZU )
C--------------------------------------------------------
C
                  T2XNET2Y = .TRUE.
                  CALL CC_XI_DEN_IA(T2XNET2Y,DIA,WORK(KW3BMAT),ISWBMAT,
     *                              WORK(KT2ZU),ISYMZU,
     *                               WORK(KT2TP),ISYM0,WORK(KINDSQWB),
     *                               LENSQWB,
     *                               B,ISYMB,D,ISYMD,WORK(KEND5),LWRK5)
C
C----------------------------------------------------------
C                 Get again T3barX_BD multipliers (using W)
C                 but now without virtual contribution:
C----------------------------------------------------------
C

                  !reuse KW3BMAT array
                  CALL DZERO(WORK(KW3BMAT),NCKIJ(ISWBMAT))
C
                  IF (LISTL(1:3).EQ.'L1 ') THEN

                    !<L3|[Y^,tau3]|HF> (occ. part)
                    CALL WX_BD_O(3,.FALSE.,.TRUE.,WORK(KTMAT),ISCKIJ,
     *                           WORK(KFOCKL1),
     *                                  ISYML1,WORK(KW3BMAT),ISWBMAT,
     *                                  WORK(KEND5),LWRK5)
 
                    ! <L2|[Y,tau3]|HF>
                    CALL WBARXBD_T2(3,B,ISYMB,D,ISYMD,WORK(KL2TP),ISYM0,
     *                              WORK(KFOCKL1),
     *                  ISYML1,WORK(KW3BMAT),ISWBMAT)
C

                    !<L2|[H^Y,tau3]|HF>
                    CALL WBARXBD_TMAT(3,
     *                    WORK(KL2TP),ISYM0,WORK(KW3BMAT),WORK(KTMAT),
     *                    ISWBMAT,WORK(KFOCKL1RCK),ISYFCKL1R,
     *                    WORK(fKW3BXVDLX2),WORK(fKW3BXVDLX1),
     *                    WORK(fKW3BXVDGX2),
     *                    WORK(fKW3BXVDGX1),WORK(KW3BXOLX1),
     *                    WORK(KW3BXOGX1),ISINT2L1R,
     *                    WORK(KEND5),LWRK5,
     *                    WORK(KINDEX),WORK(KINDEX2),
     *                    WORK(KINDSQWB),LENSQWB,
     *                    ISYMB,B,ISYMD,D)
                  END IF
C
                  !<L2Y|[H^,tau3]|HF>
                  CALL WBARXBD_TMAT(3,
     *                  WORK(KL2L1),ISYML1,WORK(KW3BMAT),WORK(KTMAT),
     *                  ISWBMAT,WORK(KFOCK0CK),ISYM0,
     *                  WORK(fKW3BXVDL2),WORK(fKW3BXVDL1),
     *                  WORK(fKW3BXVDG2),WORK(fKW3BXVDG1),
     *                  WORK(KW3BXOL1),WORK(KW3BXOG1),
     *                  ISINT2,
     *                  WORK(KEND5),LWRK5,WORK(KINDEXBL1),
     *                  WORK(KINDEXDL1),WORK(KINDSQWB),LENSQWB,
     *                  ISYMB,B,ISYMD,D)
C
                  !<L1Y|[H^,tau3]|HF>
                  CALL WBARXBD_L1(3,WORK(KL1L1),ISYML1,WORK(KTMAT),
     *                           WORK(KXIAJB),
     *                ISINT1,WORK(KW3BMAT),WORK(KEND5),LWRK5,
     *                WORK(KINDSQWB),LENSQWB,ISYMB,B,ISYMD,D)
C
                  CALL WBD_DIA(B,ISYMB,D,ISYMD,-FREQL1,ISWBMAT,
     *                         WORK(KW3BMAT),WORK(KDIAGWB),WORK(KFOCKD))
                  CALL T3_FORBIDDEN(WORK(KW3BMAT),ISYML1,ISYMB,B,
     *                              ISYMD,D)
                  !To conform with real sign of t3b multipliers
                  !(noddy code definition)
                  CALL DSCAL(NCKIJ(ISWBMAT),-ONE,WORK(KW3BMAT),1)
c
c       call sum_pt3(work(KW3BMAT),isymb,b,isymd,d,
c    *             ISWBMAT,work(kx3am),4)
c
C
C--------------------------------------------------------
C                 Write WBMAT as WBMAT^D(ai,bj,l) to disc
C--------------------------------------------------------

                  CALL WRITE_T3_DL(LUWBZU,FNWBZU,WORK(KW3BMAT),ISYML1,
     *                             ISYMD,ISYMB,B)

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
                  IF (LISTRU(1:3).EQ.'R1 ') THEN
                    AIBJCK_PERM = 4 ! means that we transform ALL occupied
                                    ! indeces

                    ! <mu3|[Z,T30]|HF> occupied contribution 

                    CALL WXBD_O(AIBJCK_PERM,WORK(KT3MAT),ISCKIJ,
     *                          WORK(KFOCKRZ),ISYMRZ,
     *                          WORK(KW3MATZ),ISWMATZ,WORK(KEND5),LWRK5)
C
                    ! <mu3|[[Z,T2],T2]|HF> 
                    CALL WXBD_T2(AIBJCK_PERM,B,ISYMB,D,ISYMD,
     *                          WORK(KT2TP),
     *                          ISYM0,WORK(KFOCKRZ),
     *                          ISYMRZ,WORK(KINDSQWZ),LENSQWZ,
     *                          WORK(KW3MATZ),ISWMATZ,WORK(KEND5),LWRK5)

                    !<mu3|[H^Z,T2]|HF> + <mu3|[H,T2^Z]|HF>
                    CALL WXBD_GROUND(AIBJCK_PERM,
     *                        WORK(KT2RZ),ISYMRZ,WORK(KWTEMP),
     *                        WORK(KT3VDG1),WORK(KT3VBG1),WORK(KT3VBG3),
     *                        WORK(KT3VDG3),
     *                        WORK(KT3OG1),ISINT2,
     *                        WORK(KW3MATZ),WORK(KEND5),LWRK5,
     *                        WORK(KINDSQWZ),LENSQWZ,
     *                        ISYMB,B,ISYMD,D)
C
                    CALL WXBD_GROUND(AIBJCK_PERM,
     *                         WORK(KT2TP),ISYM0,WORK(KWTEMP),
     *                         WORK(KW3ZVDGZ1),WORK(KW3ZVDGZ2),
     *                         WORK(KT3VBGZ3),WORK(KT3VDGZ3),
     *                         WORK(KW3ZOGZ1),ISINT2RZ,
     *                         WORK(KW3MATZ),WORK(KEND5),LWRK5,
     *                         WORK(KINDSQWZ),LENSQWZ,
     *                         ISYMB,B,ISYMD,D)

                    !Divide by the energy difference and
                    !remove the forbidden elements
                    CALL WBD_DIA(B,ISYMB,D,ISYMD,FREQRZ,ISWMATZ,
     *                         WORK(KW3MATZ),WORK(KDIAGWZ),WORK(KFOCKD))
                    CALL T3_FORBIDDEN(WORK(KW3MATZ),ISYMRZ,ISYMB,B,
     *                                ISYMD,D)

 
c                  call sum_pt3(work(KW3MATZ),isymb,b,isymd,d,
c    *                        ISWMATZ,work(kx3am),5)

C
C-------------------------------------------------------------------
C                 Contract wZ with U operator:
C                 wZU^{abc}_{i-- j- k-} = wZU^{abc}_{i-- j- k-} 
C                                       + U_{li} wZ^{abc}_{l- j- k-}
C-------------------------------------------------------------------
C

                    CALL WBD_O(WORK(KW3MATZ),ISWMATZ,WORK(KFOCKRU),
     *                   ISYMRU,
     *                   WORK(KWMATZU),ISWMATZU,WORK(KEND5),LWRK5)

c do the same to get wZU^{abc}_{i- j- k--} (put in KWMATZUD)

                    AIBJCK_PERM = 3
                    CALL WXBD_O(AIBJCK_PERM,WORK(KW3MATZ),ISWMATZ,
     *                   WORK(KFOCKRU),ISYMRU,
     *                   WORK(KWMATZUD),ISWMATZU,WORK(KEND5),LWRK5)

                  END IF

C                 ---------
C                  TERM (2)
C                 ---------

C
C----------------------------------------------------------------
C                 wZU^{abc}_{i-- j- k-} = wZU^{abc}_{i-- j- k-}
C                                       +  U(Z)_{li} t{abc}_{ljk}
C----------------------------------------------------------------
C
                  IF (LISTRU(1:3).EQ.'R1 ') THEN
                    !Calculate <mu3|[[U,T1Z],T30]|HF>
                    CALL WBD_O(WORK(KT3MAT),ISCKIJ,WORK(KFCKUZO),ISYMZU,
     *                   WORK(KWMATZU),ISWMATZU,WORK(KEND5),LWRK5)

c do the same to get wZU^{abc}_{i- j- k--} (put in KWMATZUD)

                    AIBJCK_PERM = 3
                    CALL WXBD_O(AIBJCK_PERM,WORK(KT3MAT),ISCKIJ,
     *                   WORK(KFCKUZO),ISYMZU,
     *                   WORK(KWMATZUD),ISWMATZU,WORK(KEND5),LWRK5)
                  END IF

C                 ---------
C                  TERM (3)
C                 ---------

C
C--------------------------------------------------------------------
C                 wZU^{abc}_{i-- j- k-} = wZU^{abc}_{i-- j- k-}
C                                       +  b^{abc}_{ijk}(U, t2Z, t20)
C--------------------------------------------------------------------
C
                  IF (LISTRU(1:3).EQ.'R1 ') THEN
                    !Calculate <mu3|[[U,T2Z],T20]|HF>
                    T2XNET2Z = .TRUE.
                    CALL WBD_T2(T2XNET2Z,B,ISYMB,D,ISYMD,
     *                          WORK(KT2RZ),ISYMRZ,WORK(KT2TP),ISYM0,
     *                          WORK(KFOCKRU),ISYMRU,
     *                          WORK(KINDSQWZU),LENSQWZU,WORK(KWMATZU),
     *                          ISWMATZU,WORK(KEND5),LWRK5)

c do the same to get wZU^{abc}_{i- j- k--} (put in KWMATZUD)

                    T2XNET2Z = .TRUE.
                    AIBJCK_PERM = 3 
                    CALL WXBD_T2_CUB(T2XNET2Z,AIBJCK_PERM,B,ISYMB,D,
     *                          ISYMD,
     *                          WORK(KT2RZ),ISYMRZ,WORK(KT2TP),ISYM0,
     *                          WORK(KFOCKRU),ISYMRU,
     *                          WORK(KINDSQWZU),LENSQWZU,WORK(KWMATZUD),
     *                          ISWMATZU,WORK(KEND5),LWRK5)
                  END IF

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

                  IF (LISTRU(1:3).EQ.'R1 ') THEN
                    ! <mu3|[U,T30]|HF> occupied contribution 
                    CALL WXBD_O(AIBJCK_PERM,WORK(KT3MAT),ISCKIJ,
     *                          WORK(KFOCKRU),ISYMRU,
     *                          WORK(KW3MATU),ISWMATU,WORK(KEND5),LWRK5)
C
                    ! <mu3|[[U,T2],T2]|HF> 
                    CALL WXBD_T2(AIBJCK_PERM,B,ISYMB,D,ISYMD,
     *                          WORK(KT2TP),
     *                          ISYM0,WORK(KFOCKRU),
     *                          ISYMRU,WORK(KINDSQWU),LENSQWU,
     *                          WORK(KW3MATU),ISWMATU,WORK(KEND5),LWRK5)
                  END IF

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

                  IF (LISTRU(1:3).EQ.'R1 ') THEN
                    CALL WBD_V(WORK(KW3MATZ),ISWMATZ,WORK(KFOCKRU),
     *                   ISYMRU,
     *                   WORK(KWMATZU),ISWMATZU,WORK(KEND5),LWRK5)
                  END IF

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
C                                        + Z(U)_{ad} t{dbc}_{ijk}
C--------------------------------------------------------------------
C

C                 t{dbc}_{ijk} is already there sitting in
C                 KT3MAT array.

                  !Calculate <mu3|[[Z,T1U],T30]|HF>
                  CALL WBD_V(WORK(KT3MAT),ISCKIJ,WORK(KFCKZUV),ISYMZU,
     *                 WORK(KWMATZU),ISWMATZU,WORK(KEND5),LWRK5)

C                 ---------
C                  TERM (2) (permuted)
C                 ---------

C
C--------------------------------------------------------------------
C                 wZU^{a- bc}_{i- j- k-} = wZU^{a- bc}_{i- j- k-}
C                                        + U(Z)_{ad} t{dbc}_{ijk}
C--------------------------------------------------------------------
C

C                 t{dbc}_{ijk} is already there sitting in
C                 KT3MAT array.

                  IF (LISTRU(1:3).EQ.'R1 ') THEN
                    !Calculate <mu3|[[U,T1Z],T30]|HF>
                    CALL WBD_V(WORK(KT3MAT),ISCKIJ,WORK(KFCKUZV),ISYMZU,
     *                   WORK(KWMATZU),ISWMATZU,WORK(KEND5),LWRK5)
                  END IF

                  IF (LISTRU(1:3).EQ.'R1 ') THEN
C                    ---------------------
C                     TERM (3) + (4) + (5)
C                    ---------------------

C
C-------------------------------------------------------------------------
C                    wZU^{a- bc}_{i- j- k-} = wZU^{a- bc}_{i- j- k-}
C                                           + U_{li} thetaZ^{a- bc}_{ljk} (3)
C                                           + U_{lj} thetaZ^{a- bc}_{ilk} (4)
C                                           + U_{lk} thetaZ^{a- bc}_{ijl} (5)
C-------------------------------------------------------------------------
C

C-------------------------------------------------------
C                    First we need thetaZ^{a- bc}_{ijk}...
C-------------------------------------------------------

C                    Let's reuse KW3MATZ array
                     CALL DZERO(WORK(KW3MATZ),NCKIJ(ISWMATZ))

                     ! <mu3|[Z,T30]|HF> virtual contribution 
                     CALL WBD_V(WORK(KT3MAT),ISCKIJ,
     *                          WORK(KFOCKRZ),ISYMRZ,
     *                          WORK(KW3MATZ),ISWMATZ,WORK(KEND5),LWRK5)

                     !Divide by the energy difference and
                     !remove the forbidden elements
                     CALL WBD_DIA(B,ISYMB,D,ISYMD,FREQRZ,ISWMATZ,
     *                         WORK(KW3MATZ),WORK(KDIAGWZ),WORK(KFOCKD))
                     CALL T3_FORBIDDEN(WORK(KW3MATZ),ISYMRZ,ISYMB,B,
     *                                   ISYMD,D)


c                    call sum_pt3(work(KW3MATZ),isymb,b,isymd,d,
c    *                          ISWMATZ,work(kx3am),5)

C
C-------------------------------------------------------------------------
C                    Now contract thetaZ with U operator:
C                    wZU^{a- bc}_{i- j- k-} = wZU^{a- bc}_{i- j- k-}
C                                           + U_{li} thetaZ^{a- bc}_{ljk} 
C                                           + U_{lj} thetaZ^{a- bc}_{ilk} 
C                                           + U_{lk} thetaZ^{a- bc}_{ijl} 
C-------------------------------------------------------------------------
C

                     AIBJCK_PERM = 4 ! transform all occ indeces simultanously
                     CALL WXBD_O(AIBJCK_PERM,WORK(KW3MATZ),ISWMATZ,
     *                    WORK(KFOCKRU),ISYMRU,
     *                    WORK(KWMATZU),ISWMATZU,WORK(KEND5),LWRK5)

C                    ---------------------
C                     TERM (3) + (4) + (5) (permuted)
C                    ---------------------

C
C-------------------------------------------------------------------------
C                    wZU^{a- bc}_{i- j- k-} = wZU^{a- bc}_{i- j- k-}
C                                           + Z_{li} thetaU^{a- bc}_{ljk} (3)
C                                           + Z_{lj} thetaU^{a- bc}_{ilk} (4)
C                                           + Z_{lk} thetaU^{a- bc}_{ijl} (5)
C-------------------------------------------------------------------------
C

C-------------------------------------------------------
C                    First we need thetaU^{a- bc}_{ijk}...
C-------------------------------------------------------

C                    Let's reuse KW3MATU array
                     CALL DZERO(WORK(KW3MATU),NCKIJ(ISWMATU))

                     ! <mu3|[U,T30]|HF> virtual contribution 
                     CALL WBD_V(WORK(KT3MAT),ISCKIJ,
     *                          WORK(KFOCKRU),ISYMRU,
     *                          WORK(KW3MATU),ISWMATU,WORK(KEND5),LWRK5)

                     !Divide by the energy difference and
                     !remove the forbidden elements
                     CALL WBD_DIA(B,ISYMB,D,ISYMD,FREQRU,ISWMATU,
     *                         WORK(KW3MATU),WORK(KDIAGWU),WORK(KFOCKD))
                     CALL T3_FORBIDDEN(WORK(KW3MATU),ISYMRU,ISYMB,B,
     *                                   ISYMD,D)


c                    call sum_pt3(work(KW3MATU),isymb,b,isymd,d,
c    *                          ISWMATU,work(kx3am),5)

C
C-------------------------------------------------------------------------
C                    Now contract thetaU with Z operator:
C                    wZU^{a- bc}_{i- j- k-} = wZU^{a- bc}_{i- j- k-}
C                                           + Z_{li} thetaU^{a- bc}_{ljk} 
C                                           + Z_{lj} thetaU^{a- bc}_{ilk} 
C                                           + Z_{lk} thetaU^{a- bc}_{ijl} 
C-------------------------------------------------------------------------
C

                     AIBJCK_PERM = 4 ! transform all occ indeces simultanously
                     CALL WXBD_O(AIBJCK_PERM,WORK(KW3MATU),ISWMATU,
     *                    WORK(KFOCKRZ),ISYMRZ,
     *                    WORK(KWMATZU),ISWMATZU,WORK(KEND5),LWRK5)

                  END IF ! LISTRU .EQ. 'R1 '
C
                  !Divide by the energy difference and
                  !remove the forbidden elements 
                  CALL WBD_DIA(B,ISYMB,D,ISYMD,FREQZU,ISWMATZU,
     *                        WORK(KWMATZU),WORK(KDIAGWZU),WORK(KFOCKD))
                  CALL T3_FORBIDDEN(WORK(KWMATZU),ISYMZU,ISYMB,B,
     *                              ISYMD,D)


c                 call sum_pt3(work(KWMATZU),isymb,b,isymd,d,
c    *                       ISWMATZU,work(kx3am),5)

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
C
C-------------------------------------------------------
C          Get DAB0 and DIJ0 densities
C-------------------------------------------------------
C
            QUADR = .FALSE.
            CALL CC_XI_DEN_ABIJ(QUADR,LISTR,WORK(KDAB0),WORK(KDIJ0),
     *                          .FALSE.,DUMMY,
     *                          DUMMY,IDUMMY,
     *                          IDUMMY,DUMMY,
     *                          ISYM0,ISYML1,IDUMMY,
     *                          LUT3,FNT3,LUWBMAT,FNWBMAT,
     *                          IDUMMY,CDUMMY,
     *                          DUMMY,DUMMY,
     *                          WORK(KEND5),LWRK5,ISYMD,D)

C
C-------------------------------------------------------
C          Get DAB and DIJ densities
C-------------------------------------------------------
C
            CUBIC = .TRUE.
            CALL CC_XI_DEN_ABIJ_CUB(CUBIC,LISTL,LISTRZ,LISTRU,
     *                          DAB,DIJ,DIA,ISYDEN,
     *                          WORK(KL2L1),ISYML1,
     *                          ISYMRZ,WORK(KFOCKRZ),
     *                          ISYMRU,WORK(KFOCKRU),
     *                          ISYM0,ISYML1,ISYMZU,
     *                          LUT3,FNT3,LUWBMAT,FNWBMAT,
     *                          LUTHETA,FNTHETA,
     *                          LUWZU,FNWZU,
     *                          LUWBZU,FNWBZU,
     *                          WORK(KFOCKD),FREQRZ,FREQRU,
     *                          WORK(KEND5),LWRK5,ISYMD,D)
C
         ENDDO       ! D
      ENDDO          ! ISYMD 
C
      CALL CC3_XI_DEN_AI_T1(DIA,ISYDEN,WORK(KDAB0),WORK(KDIJ0),ISYML1,
     *                            WORK(KT1ZU),ISYMZU)
C
c      write(lupri,*) 'w3x (usual) in CC3_ADENVIR_CUB'
c      write(lupri,*) 'w3xD  in CC3_ADENVIR_CUB'
c      write(lupri,*) 'w3bx  in CC3_ADENVIR_CUB'
c      write(lupri,*) 'w3x + 0.5w3xD in CC3_ADENVIR_CUB'
c      write(lupri,*) 'w3zu  in CC3_ADENVIR_CUB'
c      call print_pt3(work(kx3am),ISYM0,4)
C
C---------------------------------
C     Close the files
C---------------------------------
C
      CALL WCLOSE2(LUT3,FNT3,'DELETE')
      CALL WCLOSE2(LUWBMAT,FNWBMAT,'DELETE')
      CALL WCLOSE2(LUWBZU,FNWBZU,'DELETE')
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
C  /* Deck wbarxbd_t2 */
      SUBROUTINE WBARXBD_T2(AIBJCK_PERM,B,ISYMB,D,ISYMD,T2TP,ISYMT2,
     *                      FOCKY,ISYFKY,WMAT,ISWMAT)
C 
C IF (AIBJCK_PERM = 1) THEN  (aibjdk + aidkbj permutation)
C
C    WBD(a,i,k,j) = WBD(a,i,k,j) +
C       focky(j,B)*t2(ai,Dk) - focky(k,B)*t2(ai,Dj) 
C       focky(k,D)*t2(ai,Bj) - focky(j,D)*t2(ai,Bk) 
C
C ELSE IF (AIBJCK_PERM = 3) THEN (dkbjai + dkaibj permutation)
C       focky(j,B)*t2(Dk,ai) - focky(i,B)*t2(Dk,aj) 
C       focky(i,a)*t2(Dk,Bj) - focky(j,a)*t2(Dk,Bi) 
C
C ELSE quit with the error message.
C        
C
C     Written by P. Jorgensen and F. Pawlowski, Spring 2002.
C     (modified for AIBJCK_PERM option, Autumn 2003.)
C

      IMPLICIT NONE
C
      INTEGER AIBJCK_PERM
      INTEGER ISYMB, ISYMD, ISYMT2, ISYFKY, ISWMAT 
      INTEGER ISYMJ, KJB, KJD, ISYMK, KKB, KKD, ISYMI, ISYIJ, ISYIK
      INTEGER ISYMA, ISYAI, ISYAIK, ISYAIJ, KAIKD, KAIJD, KAIJB
      INTEGER KAIKB, KAIKJ
C
      INTEGER ISYDK,ISYDKI,KDKIA,ISYAJK,ISYAK,ISYDKJ,KIB,KDKJA
      INTEGER ISYBD,ISYKJ,KIA,KDKJB,ISYKI,KJA,KDKIB
C
#if defined (SYS_CRAY)
      REAL T2TP(*), FOCKY(*), WMAT(*)
#else
      DOUBLE PRECISION T2TP(*), FOCKY(*), WMAT(*)
#endif
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccsdinp.h>
C
      CALL QENTER('WBXT2')

      IF (AIBJCK_PERM .EQ. 1) THEN
C
C        focky(j,B)*t2(ai,Dk) - focky(k,B)*t2(ai,Dj) 
C        focky(k,D)*t2(ai,Bj) - focky(j,D)*t2(ai,Bk) 
C
         
C
C        (1)   wmat(aikj) = wmat(aikj) +  focky(j,B)*t2(ai,Dk) 
C
         ISYMJ = MULD2H(ISYFKY,ISYMB)
         ISYAIK = MULD2H(ISYMT2,ISYMD)
         DO ISYMK = 1,NSYM
            ISYAI = MULD2H(ISYAIK,ISYMK)
            DO ISYMI = 1,NSYM
               ISYMA = MULD2H(ISYAI,ISYMI)
               DO J = 1,NRHF(ISYMJ)
                  KJB = IFCVIR(ISYMJ,ISYMB) + NORB(ISYMJ)*(B - 1) + J
                  DO K = 1,NRHF(ISYMK)
                     DO I = 1,NRHF(ISYMI)
                        DO A = 1,NVIR(ISYMA)
                           KAIKD = IT2SP(ISYAIK,ISYMD)
     *                           + NCKI(ISYAIK)*(D-1)
     *                           + ISAIK(ISYAI,ISYMK)
     *                           + NT1AM(ISYAI)*(K-1)
     *                           + IT1AM(ISYMA,ISYMI)
     *                           + NVIR(ISYMA)*(I-1)
     *                           + A

                           KAIKJ = ISAIKJ(ISYAIK,ISYMJ)
     *                           + NCKI(ISYAIK)*(J-1)
     *                           + ISAIK(ISYAI,ISYMK)
     *                           + NT1AM(ISYAI)*(K-1)
     *                           + IT1AM(ISYMA,ISYMI)
     *                           + NVIR(ISYMA)*(I-1)
     *                           + A

                           WMAT(KAIKJ) = WMAT(KAIKJ) 
     *                                 + FOCKY(KJB)*T2TP(KAIKD)    
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
             
C
C        (2)  wmat(aikj) = wmat(aikj) - focky(k,B)*t2(ai,Dj) 
C
         ISYMK = MULD2H(ISYFKY,ISYMB)
         ISYAIJ = MULD2H(ISYMT2,ISYMD)
         DO ISYMJ = 1,NSYM
            ISYAI = MULD2H(ISYAIJ,ISYMJ)
            ISYAIK = MULD2H(ISYAI,ISYMK)
            DO ISYMI = 1,NSYM
               ISYMA = MULD2H(ISYAI,ISYMI)
               DO J = 1,NRHF(ISYMJ)
                  DO K = 1,NRHF(ISYMK)
                  KKB = IFCVIR(ISYMK,ISYMB) + NORB(ISYMK)*(B - 1) + K
                     DO I = 1,NRHF(ISYMI)
                        DO A = 1,NVIR(ISYMA)

                           KAIJD = IT2SP(ISYAIJ,ISYMD)
     *                           + NCKI(ISYAIJ)*(D-1)
     *                           + ISAIK(ISYAI,ISYMJ)
     *                           + NT1AM(ISYAI)*(J-1)
     *                           + IT1AM(ISYMA,ISYMI)
     *                           + NVIR(ISYMA)*(I-1)
     *                           + A

                           KAIKJ = ISAIKJ(ISYAIK,ISYMJ)
     *                           + NCKI(ISYAIK)*(J-1)
     *                           + ISAIK(ISYAI,ISYMK)
     *                           + NT1AM(ISYAI)*(K-1)
     *                           + IT1AM(ISYMA,ISYMI)
     *                           + NVIR(ISYMA)*(I-1)
     *                           + A

                           WMAT(KAIKJ) = WMAT(KAIKJ) 
     *                                 - FOCKY(KKB)*T2TP(KAIJD)
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO

         
C
C        (3)  wmat(aikj) = wmat(aikj) + focky(k,D)*t2(ai,Bj) 
C
         ISYMK = MULD2H(ISYFKY,ISYMD)
         ISYAIJ = MULD2H(ISYMT2,ISYMB)
         DO ISYMJ = 1,NSYM
            ISYAI = MULD2H(ISYAIJ,ISYMJ)
            ISYAIK = MULD2H(ISYAI,ISYMK)
            DO ISYMI = 1,NSYM
               ISYMA = MULD2H(ISYAI,ISYMI)
               DO J = 1,NRHF(ISYMJ)
                  DO K = 1,NRHF(ISYMK)
                     KKD = IFCVIR(ISYMK,ISYMD) + NORB(ISYMK)*(D - 1) + K
                     DO I = 1,NRHF(ISYMI)
                        DO A = 1,NVIR(ISYMA)

                           KAIJB = IT2SP(ISYAIJ,ISYMB)
     *                           + NCKI(ISYAIJ)*(B-1)
     *                           + ISAIK(ISYAI,ISYMJ)
     *                           + NT1AM(ISYAI)*(J-1)
     *                           + IT1AM(ISYMA,ISYMI)
     *                           + NVIR(ISYMA)*(I-1)
     *                           + A

                           KAIKJ = ISAIKJ(ISYAIK,ISYMJ)
     *                           + NCKI(ISYAIK)*(J-1)
     *                           + ISAIK(ISYAI,ISYMK)
     *                           + NT1AM(ISYAI)*(K-1)
     *                           + IT1AM(ISYMA,ISYMI)
     *                           + NVIR(ISYMA)*(I-1)
     *                           + A

                           WMAT(KAIKJ) = WMAT(KAIKJ) 
     *                                 + FOCKY(KKD)*T2TP(KAIJB)
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO

C
C        (4)  wmat(aikj) = wmat(aikj) -  focky(j,D)*t2(ai,Bk) 
C
         ISYMJ = MULD2H(ISYFKY,ISYMD)
         ISYAIK = MULD2H(ISYMT2,ISYMB)
         DO ISYMK = 1,NSYM
            ISYAI = MULD2H(ISYAIK,ISYMK)
            DO ISYMI = 1,NSYM
               ISYMA = MULD2H(ISYAI,ISYMI)
               DO J = 1,NRHF(ISYMJ)
                  KJD = IFCVIR(ISYMJ,ISYMD) + NORB(ISYMJ)*(D - 1) + J
                  DO K = 1,NRHF(ISYMK)
                     DO I = 1,NRHF(ISYMI)
                        DO A = 1,NVIR(ISYMA)

                           KAIKB = IT2SP(ISYAIK,ISYMB)
     *                           + NCKI(ISYAIK)*(B-1)
     *                           + ISAIK(ISYAI,ISYMK)
     *                           + NT1AM(ISYAI)*(K-1)
     *                           + IT1AM(ISYMA,ISYMI)
     *                           + NVIR(ISYMA)*(I-1)
     *                           + A

                           KAIKJ = ISAIKJ(ISYAIK,ISYMJ)
     *                           + NCKI(ISYAIK)*(J-1)
     *                           + ISAIK(ISYAI,ISYMK)
     *                           + NT1AM(ISYAI)*(K-1)
     *                           + IT1AM(ISYMA,ISYMI)
     *                           + NVIR(ISYMA)*(I-1)
     *                           + A

                           WMAT(KAIKJ) = WMAT(KAIKJ)
     *                                 - FOCKY(KJD)*T2TP(KAIKB)
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
C
      ELSE IF (AIBJCK_PERM .EQ. 3) THEN
C
C       focky(j,B)*t2(Dk,ai) - focky(i,B)*t2(Dk,aj) 
C       focky(i,a)*t2(Dk,Bj) - focky(j,a)*t2(Dk,Bi) 
C

C
C        (1)   wmat(aikj) = wmat(aikj) +  focky(j,B)*t2(Dk,ai) 
C
         ISYMJ = MULD2H(ISYFKY,ISYMB)
         ISYAIK = MULD2H(ISYMT2,ISYMD)
         DO ISYMK = 1,NSYM
            ISYDK = MULD2H(ISYMD,ISYMK)
            ISYAI = MULD2H(ISYAIK,ISYMK)
            DO ISYMI = 1,NSYM
               ISYDKI = MULD2H(ISYDK,ISYMI)
               ISYMA = MULD2H(ISYAI,ISYMI)
               DO J = 1,NRHF(ISYMJ)
                  KJB = IFCVIR(ISYMJ,ISYMB) + NORB(ISYMJ)*(B - 1) + J
                  DO K = 1,NRHF(ISYMK)
                     DO I = 1,NRHF(ISYMI)
                        DO A = 1,NVIR(ISYMA)
                           KDKIA = IT2SP(ISYDKI,ISYMA)
     *                           + NCKI(ISYDKI)*(A-1)
     *                           + ISAIK(ISYDK,ISYMI)
     *                           + NT1AM(ISYDK)*(I-1)
     *                           + IT1AM(ISYMD,ISYMK)
     *                           + NVIR(ISYMD)*(K-1)
     *                           + D

                           KAIKJ = ISAIKJ(ISYAIK,ISYMJ)
     *                           + NCKI(ISYAIK)*(J-1)
     *                           + ISAIK(ISYAI,ISYMK)
     *                           + NT1AM(ISYAI)*(K-1)
     *                           + IT1AM(ISYMA,ISYMI)
     *                           + NVIR(ISYMA)*(I-1)
     *                           + A

                           WMAT(KAIKJ) = WMAT(KAIKJ)
     *                                 + FOCKY(KJB)*T2TP(KDKIA)
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO

C        (2)  wmat(aikj) = wmat(aikj) - focky(i,B)*t2(Dk,aj)
C
         ISYMI = MULD2H(ISYFKY,ISYMB)
         ISYAJK = MULD2H(ISYMT2,ISYMD)
         DO ISYMJ = 1,NSYM
            ISYAK = MULD2H(ISYAJK,ISYMJ)
            ISYAIK = MULD2H(ISYAK,ISYMI)
            DO ISYMK = 1,NSYM
               ISYDK = MULD2H(ISYMK,ISYMD)
               ISYDKJ = MULD2H(ISYDK,ISYMJ)
               ISYMA = MULD2H(ISYAK,ISYMK)
               ISYAI = MULD2H(ISYAIK,ISYMK)
               DO J = 1,NRHF(ISYMJ)
                  DO K = 1,NRHF(ISYMK)
                     DO I = 1,NRHF(ISYMI)
                        KIB = IFCVIR(ISYMI,ISYMB)+NORB(ISYMI)*(B-1) + I
                        DO A = 1,NVIR(ISYMA)

                           KDKJA = IT2SP(ISYDKJ,ISYMA)
     *                           + NCKI(ISYDKJ)*(A-1)
     *                           + ISAIK(ISYDK,ISYMJ)
     *                           + NT1AM(ISYDK)*(J-1)
     *                           + IT1AM(ISYMD,ISYMK)
     *                           + NVIR(ISYMD)*(K-1)
     *                           + D

                           KAIKJ = ISAIKJ(ISYAIK,ISYMJ)
     *                           + NCKI(ISYAIK)*(J-1)
     *                           + ISAIK(ISYAI,ISYMK)
     *                           + NT1AM(ISYAI)*(K-1)
     *                           + IT1AM(ISYMA,ISYMI)
     *                           + NVIR(ISYMA)*(I-1)
     *                           + A

                           WMAT(KAIKJ) = WMAT(KAIKJ)
     *                                 - FOCKY(KIB)*T2TP(KDKJA)
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
C
C        (3)  wmat(aikj) = wmat(aikj) + focky(i,a)*t2(Dk,Bj)
C
         ISYBD = MULD2H(ISYMD,ISYMB)
         ISYKJ = MULD2H(ISYMT2,ISYBD)
         ISYDKJ = MULD2H(ISYMD,ISYKJ)
         DO ISYMJ = 1,NSYM
            ISYMK = MULD2H(ISYKJ,ISYMJ)
            ISYDK = MULD2H(ISYMD,ISYMK)
            DO ISYMI = 1,NSYM
               ISYMA = MULD2H(ISYFKY,ISYMI)
               ISYAI = MULD2H(ISYMA,ISYMI)
               ISYAIK = MULD2H(ISYAI,ISYMK)
               DO J = 1,NRHF(ISYMJ)
                  DO K = 1,NRHF(ISYMK)
                     DO I = 1,NRHF(ISYMI)
                        DO A = 1,NVIR(ISYMA)
                           KIA = IFCVIR(ISYMI,ISYMA)+NORB(ISYMI)*(A-1)+I

                           KDKJB = IT2SP(ISYDKJ,ISYMB)
     *                           + NCKI(ISYDKJ)*(B-1)
     *                           + ISAIK(ISYDK,ISYMJ)
     *                           + NT1AM(ISYDK)*(J-1)
     *                           + IT1AM(ISYMD,ISYMK)
     *                           + NVIR(ISYMD)*(K-1)
     *                           + D

                           KAIKJ = ISAIKJ(ISYAIK,ISYMJ)
     *                           + NCKI(ISYAIK)*(J-1)
     *                           + ISAIK(ISYAI,ISYMK)
     *                           + NT1AM(ISYAI)*(K-1)
     *                           + IT1AM(ISYMA,ISYMI)
     *                           + NVIR(ISYMA)*(I-1)
     *                           + A

                           WMAT(KAIKJ) = WMAT(KAIKJ)
     *                                 + FOCKY(KIA)*T2TP(KDKJB)
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO

C
C        (4)  wmat(aikj) = wmat(aikj) -  focky(j,a)*t2(Dk,Bi)
C
         ISYDKI = MULD2H(ISYMT2,ISYMB)
         ISYKI  = MULD2H(ISYDKI,ISYMD)
         DO ISYMJ = 1,NSYM
            ISYMA = MULD2H(ISYFKY,ISYMJ)
            DO ISYMK = 1,NSYM
               ISYMI = MULD2H(ISYKI,ISYMK)
               ISYAI = MULD2H(ISYMA,ISYMI)
               ISYDK = MULD2H(ISYDKI,ISYMI)
               ISYAIK = MULD2H(ISYAI,ISYMK)
               DO J = 1,NRHF(ISYMJ)
                  DO K = 1,NRHF(ISYMK)
                     DO I = 1,NRHF(ISYMI)
                        DO A = 1,NVIR(ISYMA)
                           KJA = IFCVIR(ISYMJ,ISYMA)+NORB(ISYMJ)*(A-1)+J

                           KDKIB = IT2SP(ISYDKI,ISYMB)
     *                           + NCKI(ISYDKI)*(B-1)
     *                           + ISAIK(ISYDK,ISYMI)
     *                           + NT1AM(ISYDK)*(I-1)
     *                           + IT1AM(ISYMD,ISYMK)
     *                           + NVIR(ISYMD)*(K-1)
     *                           + D

                           KAIKJ = ISAIKJ(ISYAIK,ISYMJ)
     *                           + NCKI(ISYAIK)*(J-1)
     *                           + ISAIK(ISYAI,ISYMK)
     *                           + NT1AM(ISYAI)*(K-1)
     *                           + IT1AM(ISYMA,ISYMI)
     *                           + NVIR(ISYMA)*(I-1)
     *                           + A

                           WMAT(KAIKJ) = WMAT(KAIKJ)
     *                                 - FOCKY(KJA)*T2TP(KDKIB)
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
C
      ELSE
         WRITE(LUPRI,*) 'AIBJCK_PERM = ', AIBJCK_PERM
         WRITE(LUPRI,*) 'WBARXBD_T2 works for AIBJCK_PERM = 1 or 3'
         CALL QUIT('Wrong AIBJCK_PERM option in WBARXBD_T2')
      END IF

      CALL QEXIT('WBXT2')
C
      RETURN
      END
C  /* Deck wbarxbd_l1 */
      SUBROUTINE WBARXBD_L1(AIBJCK_PERM,T1AM,ISYMT1,TMAT,XIAJB,
     *                      ISINT1,
     *                      WMAT,WORK,LWORK,
     *                      INDSQ,LENSQ,ISYMB,B,ISYMC,C)
*---------------------------------------------------------------------*
*
*    Purpose: compute Tbar1^Y contribution to triples component of 
*    first-order multipliers vector:
*
*    <Tbar1^Y|[H_0^,tau3]|HF> = P^(abc)_(ijk) (  t1bar^Y(ai)*L(jbkc) 
*                                              - t1bar^Y(ak)*L(jbic)  )
*
*             Use W intermmediates: 
*
* IF (AIBJCK_PERM .EQ. 1) THEN (aibjck + aickbj permutation)
*
*    WMAT^BC(aikj) = WMAT^BC(aikj) + T1(ai)*L(jBkC)             
*                                  - T1(ak)*L(jBiC)
*                                  + T1(ai)*L(kCjB)
*                                  - T1(aj)*L(kCiB)
*
* ELSE IF (AIBJCK_PERM = 3) THEN (ckbjai + ckaibj permutation)
*
*    WMAT^BC(aikj) = WMAT^BC(aikj) + T1(Ck)*L(jBia)
*                                  - T1(Ci)*L(jBka)
*                                  + T1(Ck)*L(iajB)
*                                  - T1(Cj)*L(iakB)
*
* ELSE quit with the error message.
*
*    Written by Filip Pawlowski, Fall 2002, Aarhus
*    (modified for AIBJCK_PERM option, Fall 2003.)
*            
*=====================================================================*
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INTEGER AIBJCK_PERM
      INTEGER ISYMT1, ISINT1, LENSQ, ISYMB, ISYMC
      INTEGER ISYMBC, ISYRES, JSAIKJ, LENGTH, ISYMK, ISYMJ 
      INTEGER ISYMAI, ISYAIK, ISYMJK, ISYMCK, NBJ, NCK, ISYMBJ
      INTEGER NCKBJ, NBJCK, NAI, NAIKJ 
      INTEGER INDEX, INDSQ(LENSQ,6)
      INTEGER LWORK
C
      INTEGER ISYBIA,NAIBJ,ISYAIB,NBJAI
C
#if defined (SYS_CRAY)
      REAL T1AM(*), TMAT(*), XIAJB(*) 
      REAL WMAT(*),WORK(LWORK)
      real xnormval,ddot
#else
      DOUBLE PRECISION T1AM(*), TMAT(*), XIAJB(*) 
      DOUBLE PRECISION WMAT(*),WORK(LWORK)  
      double precision xnormval,ddot
#endif
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('WBX_L1')
C
      ISYRES  = MULD2H(ISYMT1,ISINT1)
C
      ISYMBC = MULD2H(ISYMB,ISYMC)
      JSAIKJ = MULD2H(ISYRES,ISYMBC)
      LENGTH = NCKIJ(JSAIKJ)
C
      IF (AIBJCK_PERM .EQ. 1) THEN
C
C-----------------------------------------------
C        First contribution from both T1 terms
C
C        WMAT^BC(aikj) = WMAT^BC(aikj) + T1(ai)*L(jBkC)             
C                                     - T1(ak)*L(jBiC)
C
C-----------------------------------------------
C
         ISYMJK = MULD2H(ISYMBC,ISINT1)
C
C------------------------------------------
C        Contract the integrals with T1.
C------------------------------------------
C
         CALL DZERO(TMAT,LENGTH)
C
         ISYMAI = ISYMT1
         DO ISYMJ = 1, NSYM
            ISYMK  = MULD2H(ISYMJK,ISYMJ)
            ISYAIK = MULD2H(ISYMK,ISYMAI)
            ISYMCK = MULD2H(ISYMC,ISYMK)
            ISYMBJ = MULD2H(ISYMB,ISYMJ)
C
            DO J = 1, NRHF(ISYMJ)
               NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
C
               DO K = 1, NRHF(ISYMK)
C
                  NCK = IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K - 1) + C
C
                  NCKBJ = IT2AM(ISYMCK,ISYMBJ) + INDEX(NCK,NBJ)
C
                  DO NAI = 1, NT1AM(ISYMAI)
C
                     NAIKJ = ISAIKJ(ISYAIK,ISYMJ)
     *                     + NCKI(ISYAIK)*(J - 1)
     *                     + ICKI(ISYMAI,ISYMK)
     *                     + NT1AM(ISYMAI)*(K - 1) + NAI
C
                     TMAT(NAIKJ) = T1AM(NAI)*XIAJB(NCKBJ)
C
                  ENDDO
               ENDDO
            ENDDO
C
         ENDDO
C
C-------------------------------------------
C        Sum the result into WMAT.
C-------------------------------------------
C
         JSAIKJ = MULD2H(ISYMAI,ISYMJK)
         DO I = 1, NCKIJ(JSAIKJ)
C            First :  WMAT^BC(aikj) = WMAT^BC(aikj) + T1(ai)*L(jBkC)             
             WMAT(I) = WMAT(I) + TMAT(I)
C            Second : WMAT^BC(aikj) = WMAT^BC(aikj) - T1(ak)*L(jBiC)
             WMAT(I) = WMAT(I) - TMAT(INDSQ(I,1))
         ENDDO
C
C-----------------------------------------------
C        Second contribution from both T1 terms
C
C        WMAT^BC(aikj) = WMAT^BC(aikj) + T1(ai)*L(kCjB)
C                                     - T1(aj)*L(kCiB)
C
C
C-----------------------------------------------
C
         ISYMJK = MULD2H(ISYMBC,ISINT1)
C
C------------------------------------------
C        Contract the integrals with T1.
C------------------------------------------
C
         CALL DZERO(TMAT,LENGTH)
C
         ISYMAI = ISYMT1
         DO ISYMK = 1, NSYM
            ISYAIK = MULD2H(ISYMK,ISYMAI)
            ISYMJ  = MULD2H(ISYMJK,ISYMK)
            ISYMCK = MULD2H(ISYMC,ISYMK)
            ISYMBJ = MULD2H(ISYMB,ISYMJ)
C
            DO K = 1, NRHF(ISYMK)
C
               NCK = IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K - 1) + C
C
               DO J = 1, NRHF(ISYMJ)
                  NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
C
                  NBJCK = IT2AM(ISYMBJ,ISYMCK) + INDEX(NBJ,NCK)
C
                  DO NAI = 1, NT1AM(ISYMAI)
C
                     NAIKJ = ISAIKJ(ISYAIK,ISYMJ)
     *                     + NCKI(ISYAIK)*(J - 1)
     *                     + ICKI(ISYMAI,ISYMK)
     *                     + NT1AM(ISYMAI)*(K-1) + NAI
C
                     TMAT(NAIKJ) = T1AM(NAI)*XIAJB(NBJCK)
C
                  ENDDO
               ENDDO
            ENDDO
C
         ENDDO
c
C
C-------------------------------------------
C        Sum the result into WMAT.
C-------------------------------------------
C
         JSAIKJ = MULD2H(ISYMAI,ISYMJK)
         DO I = 1, NCKIJ(JSAIKJ)
C            First :  WMAT^BC(aikj) = WMAT^BC(aikj) + T1(ai)*L(kCjB)
             WMAT(I) = WMAT(I) + TMAT(I)
C            Second : WMAT^BC(aikj) = WMAT^BC(aikj) - T1(aj)*L(kCiB)
             WMAT(I) = WMAT(I) - TMAT(INDSQ(I,5))
         ENDDO
C
      ELSE IF (AIBJCK_PERM .EQ. 3) THEN
C
C-----------------------------------------------
C        First contribution from both T1 terms
C
C        WMAT^BC(aikj) = WMAT^BC(aikj) + T1(Ck)*L(jBia)
C                                      - T1(Ci)*L(jBka)
C
C-----------------------------------------------
C
C
C------------------------------------------
C        Contract the integrals with T1.
C------------------------------------------
C
         CALL DZERO(TMAT,LENGTH)
C
         ISYMCK = ISYMT1
         ISYMK  = MULD2H(ISYMCK,ISYMC)
         DO ISYMJ = 1, NSYM
            ISYBIA = MULD2H(ISINT1,ISYMJ)
            ISYMAI = MULD2H(ISYBIA,ISYMB)
            ISYMBJ = MULD2H(ISYMB,ISYMJ)
            ISYAIK = MULD2H(ISYMAI,ISYMK)
C
            DO J = 1, NRHF(ISYMJ)
               NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
C
               DO K = 1, NRHF(ISYMK)
C
                  NCK = IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K - 1) + C
C
                  DO NAI = 1, NT1AM(ISYMAI)
C
                     NAIBJ = IT2AM(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
C
                     NAIKJ = ISAIKJ(ISYAIK,ISYMJ)
     *                     + NCKI(ISYAIK)*(J - 1)
     *                     + ICKI(ISYMAI,ISYMK)
     *                     + NT1AM(ISYMAI)*(K - 1) + NAI
C
                     TMAT(NAIKJ) = T1AM(NCK)*XIAJB(NAIBJ)
C
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
C
C-------------------------------------------
C        Sum the result into WMAT.
C-------------------------------------------
C
c        JSAIKJ = MULD2H(ISYMAI,ISYMJK)
c        DO I = 1, NCKIJ(JSAIKJ)
         DO I = 1, LENGTH
C            First :  WMAT^BC(aikj) = WMAT^BC(aikj) + T1(Ck)*L(jBia)
 
             WMAT(I) = WMAT(I) + TMAT(I)
C            Second : WMAT^BC(aikj) = WMAT^BC(aikj) - T1(Ci)*L(jBka)
             WMAT(I) = WMAT(I) - TMAT(INDSQ(I,1))
         ENDDO
C
C-----------------------------------------------
C        Second contribution from both T1 terms
C
C        WMAT^BC(aikj) = WMAT^BC(aikj) + T1(Ck)*L(iajB)
C                                      - T1(Cj)*L(iakB)
C
C-----------------------------------------------
C

C
C------------------------------------------
C        Contract the integrals with T1.
C------------------------------------------
C
         CALL DZERO(TMAT,LENGTH)
C
         ISYMCK = ISYMT1
         ISYMK = MULD2H(ISYMCK,ISYMC)
         DO ISYMJ = 1, NSYM
            ISYAIB = MULD2H(ISINT1,ISYMJ)
            ISYMAI = MULD2H(ISYAIB,ISYMB)
            ISYAIK = MULD2H(ISYMAI,ISYMK)
            ISYMBJ = MULD2H(ISYMB,ISYMJ)
C
            DO K = 1, NRHF(ISYMK)
C
               NCK = IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K - 1) + C
C
               DO J = 1, NRHF(ISYMJ)
                  NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
C
                  DO NAI = 1, NT1AM(ISYMAI)
C
                     NBJAI = IT2AM(ISYMBJ,ISYMAI) + INDEX(NBJ,NAI)
C
                     NAIKJ = ISAIKJ(ISYAIK,ISYMJ)
     *                     + NCKI(ISYAIK)*(J - 1)
     *                     + ICKI(ISYMAI,ISYMK)
     *                     + NT1AM(ISYMAI)*(K-1) + NAI
C
                     TMAT(NAIKJ) = T1AM(NCK)*XIAJB(NBJAI)
C
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
C
C-------------------------------------------
C        Sum the result into WMAT.
C-------------------------------------------
C
c        JSAIKJ = MULD2H(ISYMAI,ISYMJK)
c        DO I = 1, NCKIJ(JSAIKJ)
         DO I = 1, LENGTH
C            First :  WMAT^BC(aikj) = WMAT^BC(aikj) + T1(Ck)*L(iajB)
             WMAT(I) = WMAT(I) + TMAT(I)
C            Second : WMAT^BC(aikj) = WMAT^BC(aikj) - T1(Cj)*L(iakB)
             WMAT(I) = WMAT(I) - TMAT(INDSQ(I,3))
         ENDDO
C
      ELSE 
         WRITE(LUPRI,*) 'AIBJCK_PERM = ', AIBJCK_PERM
         WRITE(LUPRI,*) 'WBARXBD_T2 works for AIBJCK_PERM = 1 or 3'
         CALL QUIT('Wrong AIBJCK_PERM option in WBARXBD_L1')
      END IF

C
      CALL QEXIT('WBX_L1')
C
      RETURN
      END
C  /* Deck wbarxbd_tmat */
      SUBROUTINE WBARXBD_TMAT(AIBJCK_PERM,
     *                      T2TP,ISYMT2,WMAT,TMAT,ISWMAT,FOCK,
     *                      ISYFOCK,VLDKBC,VLDKCB,VGDKBC,VGDKCB,TROCCL,
     *                      TROCCG,ISINT2,WORK,LWORK,INDAJLB,
     *                      INDAJLC,INDSQ,LENSQ,ISYMB,B,ISYMC,C)
C
C     Written by Kasper Hald, Fall 2001.
C     (generalized for AIBJCK_PERM, F. Pawlowski, Fall 2003.)
C
C     General symmetry: ISINT2 is symmetry of integrals 
C                       ISYMT2 is symmetry of T2TP
C
C     Virtual integrals stored as:
C          L(kcd^b) -> IC(d^kB):  VLDKBC
C          L(kcd^b) -> IB(d^kC):  VLDKCB
C          g(kcd^b) -> IC(d^kB):  VGDKBC
C          g(kcd^b) -> IB(d^kC):  VGDKCB

C     Occupied integrals stored as:
C          L(ia|j k-) -> I(k-,i,j,a): TROCCL 
C          g(ia|j k-) -> I(k-,i,j,a): TROCCG 
C
C
C IF (AIBJCK_PERM .EQ. 1) THEN (aibjck + aickbj permutation)
C
C     WMAT^BC(aikj) = WMAT^BC(aikj) + T2TP(aijB)*F(kC) 
C                                   - T2TP(aikB)*F(jC)
C                                   + T2TP(aikC)*F(jB)
C                                   - T2TP(aijC)*F(kB)
C
C                   + T2TP(aijd)*L(d^BkC)
C                   - T2TP(ajkd)*g(iBd^C)
C                   + T2TP(aikd)*L(d^CjB)
C                   - T2TP(akjd)*g(iCd^B)
C
C                   + T2TP(ailB)*L(jl^kC)
C                   - T2TP(alkB)*g(il^jC)
C                   + T2TP(ailC)*L(kl^jB)
C                   - T2TP(aljC)*g(il^kB)
C
C ELSE IF (AIBJCK_PERM = 3) THEN (ckbjai + ckaibj permutation)
C
C     WMAT^BC(aikj) = WMAT^BC(aikj) + T2TP(CkjB)*F(ia) 
C                                   - T2TP(CkiB)*F(ja)
C                                   + T2TP(Ckia)*F(jB)
C                                   - T2TP(Ckja)*F(iB)
C
C                   + T2TP(Ckjd)*L(d^Bia)
C                   - T2TP(Cjid)*g(kBd^a)
C                   + T2TP(Ckid)*L(d^ajB)
C                   - T2TP(Cijd)*g(kad^B)
C
C                   + T2TP(CklB)*L(jl^ia)
C                   - T2TP(CliB)*g(kl^ja)
C                   + T2TP(Ckla)*L(il^jB)
C                   - T2TP(Clja)*g(kl^iB)
C
C ELSE quit with the error message.

C

      IMPLICIT NONE 
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INTEGER AIBJCK_PERM
      INTEGER ISYMT2,ISWMAT,ISINT2,ISYMB,ISYMC,ISYRES,ISYMBC
      INTEGER JSAIKJ,ISYMK,ISYAIJ,ISYMJ,ISYMAI,ISYAIK,ISYMI
      INTEGER ISYMA,ISYMDK,ISYMD,ISYMDI,ISYAJK,ISYMDJ,ISYAIL
      INTEGER ISYLKJ,ISYMLK,ISYML,ISYALK,ISYLJI,ISYAKJ,ISYMLJ
      INTEGER ISYMAK,ISYAJL,ISYLKI,ISYMAJ,ISYFOCK
      INTEGER NAI,NAIJB,NCK,NAIKJ,NCJ,NAIKB,NBJ,NAIKC,NAIJC,NBK
      INTEGER NTOAIJ,NVIRD,NTOAJK,NTOAIK,NTOAKJ,NTOTAI,NRHFL
      INTEGER NTOTAK,NTOTAJ
      INTEGER INDAJLB,INDAJLC,LENSQ,INDSQ(LENSQ,6),INDEX
      INTEGER KOFF1,KOFF2,KOFF3,KALK,KEND1,KALJ,KOFF
      INTEGER LWORK,LENGTH,LWRK1
C
      INTEGER ISYCKJ,ISYKJ,ISYMCK,NCKJB
      INTEGER ISYAID,ISYDJK,KAID,ISYDJ,ISYAI
      INTEGER ILOOP
      INTEGER ISYLK,KLK,KAJIL,NTOTL,ISYAJI,NTOTAJI

#if defined (SYS_CRAY)
      REAL T2TP(*),WMAT(*),TMAT(*),FOCK(*)
      REAL VLDKBC(*),VLDKCB(*),VGDKBC(*),VGDKCB(*),TROCCL(*),TROCCG(*)
      REAL WORK(*)
      REAL XWMAT,ONE,DDOT
#else
      DOUBLE PRECISION T2TP(*),WMAT(*),TMAT(*),FOCK(*)
      DOUBLE PRECISION VLDKBC(*),VLDKCB(*),VGDKBC(*),VGDKCB(*),TROCCL(*)
      DOUBLE PRECISION TROCCG(*),WORK(*)
      DOUBLE PRECISION XWMAT,ONE,DDOT
#endif
C
      PARAMETER(ONE = 1.0D0)
C
      CALL QENTER('WBXTMT')
C
      ISYRES = MULD2H(ISYMT2,ISINT2)
C
      ISYMBC = MULD2H(ISYMB,ISYMC)
      JSAIKJ = MULD2H(ISYRES,ISYMBC)
C
      IF (JSAIKJ .NE. ISWMAT) THEN
         WRITE(LUPRI,*)'JSAIKJ ', JSAIKJ
         WRITE(LUPRI,*)'ISWMAT ', ISWMAT
         WRITE(LUPRI,*)'ISWMAT and JSAIKJ should be equal '
         CALL QUIT('Symmetry inconsistency in WBARXBD_TMAT')
      END IF
C
      LENGTH = NCKIJ(JSAIKJ)
C
      IF (AIBJCK_PERM .EQ. 1) THEN
C
C--------------------------------------------------------------------------
C    
C        WMAT^BC(aikj) = WMAT^BC(aikj) + T2TP(aijB)*F(kC) 
C                                      - T2TP(aikB)*F(jC)
C                                      + T2TP(aikC)*F(jB)
C                                      - T2TP(aijC)*F(kB)
C
C--------------------------------------------------------------------------
C        Contribution from both Fock terms:
C--------------------------------------------------------------------------
C
C
C        WMAT^BC(aikj) = WMAT^BC(aikj) + T2TP(aijB)*F(kC) 
C
         ISYMK  = MULD2H(ISYFOCK,ISYMC)
         ISYAIJ = MULD2H(ISYMT2,ISYMB)
C
         DO ISYMJ = 1, NSYM
            ISYMAI = MULD2H(ISYAIJ,ISYMJ)
            ISYAIK = MULD2H(ISYMK,ISYMAI)
            DO ISYMI = 1, NSYM
               ISYMA = MULD2H(ISYMAI,ISYMI)
C
               DO J = 1, NRHF(ISYMJ)
C
                  DO I = 1, NRHF(ISYMI)
                  DO A = 1, NVIR(ISYMA)
C
                     NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
C
C                    Index for sorted T2 amplitudes.
C
                     NAIJB = IT2SP(ISYAIJ,ISYMB)
     *                     + NCKI(ISYAIJ)*(B - 1)
     *                     + ICKI(ISYMAI,ISYMJ)
     *                     + NT1AM(ISYMAI)*(J - 1) + NAI
C
                     DO K = 1, NRHF(ISYMK)
C
                        NCK = IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K-1) +C 
                        NAIKJ = ISAIKJ(ISYAIK,ISYMJ)
     *                        + NCKI(ISYAIK)*(J - 1)
     *                        + ICKI(ISYMAI,ISYMK)
     *                        + NT1AM(ISYMAI)*(K-1)
     *                        + NAI
                               
C
C                       Fock 1.0 contribution addWMAT
C

                        WMAT(NAIKJ) = WMAT(NAIKJ)+T2TP(NAIJB)*FOCK(NCK)
C
                     ENDDO
                  ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
C
C        WMAT^BC(aikj) = WMAT^BC(aikj) - T2TP(aikB)*F(jC)
C
         ISYMJ  = MULD2H(ISYFOCK,ISYMC)
         ISYAIK = MULD2H(ISYMT2,ISYMB)
C
         DO ISYMK = 1, NSYM
            ISYMAI = MULD2H(ISYAIK,ISYMK)
            DO ISYMI = 1, NSYM
               ISYMA = MULD2H(ISYMAI,ISYMI)
C
               DO J = 1, NRHF(ISYMJ)
                  NCJ = IT1AM(ISYMC,ISYMJ) + NVIR(ISYMC)*(J-1) + C
C
                  DO I = 1, NRHF(ISYMI)
                  DO A = 1, NVIR(ISYMA)
C
                     NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
C
C                    Index for sorted T2 amplitudes.
C
                     DO K = 1, NRHF(ISYMK)
C
                        NAIKB = IT2SP(ISYAIK,ISYMB)
     *                        + NCKI(ISYAIK)*(B - 1)
     *                        + ICKI(ISYMAI,ISYMK) 
     *                        + NT1AM(ISYMAI)*(K - 1) + NAI
C
                        NAIKJ = ISAIKJ(ISYAIK,ISYMJ)
     *                        + NCKI(ISYAIK)*(J - 1)
     *                        + ICKI(ISYMAI,ISYMK)
     *                        + NT1AM(ISYMAI)*(K-1)
     *                        + NAI

C
C                       Fock 2.0 contribution addWMAT
C
                        WMAT(NAIKJ) = WMAT(NAIKJ)-T2TP(NAIKB)*FOCK(NCJ)
C
                     ENDDO
                  ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

C
C        WMAT^BC(aikj) = WMAT^BC(aikj) + T2TP(aikC)*F(jB)
C                                     
         ISYMJ  = MULD2H(ISYFOCK,ISYMB)
         ISYAIK = MULD2H(ISYMT2,ISYMC)
C
         DO ISYMK = 1, NSYM
            ISYMAI = MULD2H(ISYAIK,ISYMK)
            DO ISYMI = 1, NSYM
               ISYMA = MULD2H(ISYMAI,ISYMI)
C
               DO J = 1, NRHF(ISYMJ)
                  NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J-1) + B
C
                  DO I = 1, NRHF(ISYMI)
                  DO A = 1, NVIR(ISYMA)
C 
                     NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
C
C                    Index for sorted T2 amplitudes.
C 
                     DO K = 1, NRHF(ISYMK)
C                       
                        NAIKC = IT2SP(ISYAIK,ISYMC)
     *                        + NCKI(ISYAIK)*(C - 1)
     *                        + ICKI(ISYMAI,ISYMK) 
     *                        + NT1AM(ISYMAI)*(K - 1) + NAI
C
                        NAIKJ = ISAIKJ(ISYAIK,ISYMJ)
     *                        + NCKI(ISYAIK)*(J - 1)
     *                        + ICKI(ISYMAI,ISYMK)
     *                        + NT1AM(ISYMAI)*(K-1)
     *                        + NAI

C                       
C                       Fock 3.0 contribution addWMAT
C
                        WMAT(NAIKJ) = WMAT(NAIKJ)+T2TP(NAIKC)*FOCK(NBJ)
C                    
                     ENDDO
                  ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

C
C        WMAT^BC(aikj) = WMAT^BC(aikj) - T2TP(aijC)*F(kB)
C
         ISYMK  = MULD2H(ISYFOCK,ISYMB)
         ISYAIJ = MULD2H(ISYMT2,ISYMC)
C
         DO ISYMJ = 1, NSYM
            ISYMAI = MULD2H(ISYAIJ,ISYMJ)
            ISYAIK = MULD2H(ISYMK,ISYMAI)
            DO ISYMI = 1, NSYM
               ISYMA = MULD2H(ISYMAI,ISYMI)
C
               DO J = 1, NRHF(ISYMJ)
C
                  DO I = 1, NRHF(ISYMI)
                  DO A = 1, NVIR(ISYMA)
C
                     NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
C
C                    Index for sorted T2 amplitudes.
C
                     NAIJC = IT2SP(ISYAIJ,ISYMC)
     *                     + NCKI(ISYAIJ)*(C - 1)
     *                     + ICKI(ISYMAI,ISYMJ)
     *                     + NT1AM(ISYMAI)*(J - 1) + NAI
C
                     DO K = 1, NRHF(ISYMK)
C
                        NBK = IT1AM(ISYMB,ISYMK) + NVIR(ISYMB)*(K-1) + B
                        NAIKJ = ISAIKJ(ISYAIK,ISYMJ)
     *                        + NCKI(ISYAIK)*(J - 1)
     *                        + ICKI(ISYMAI,ISYMK)
     *                        + NT1AM(ISYMAI)*(K-1)
     *                        + NAI

C
C                       Fock 4.0 contribution addWMAT
C
                        WMAT(NAIKJ) = WMAT(NAIKJ)-T2TP(NAIJC)*FOCK(NBK)
C
                     ENDDO
                  ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
C
C------------------------------------------------------------
C        First virtual contribution of L term.
C
C         WMAT^BC(aikj) = WMAT^BC(aikj) + T2TP(aijd)*L(d^BkC)
C------------------------------------------------------------
C
         ISYMBC = MULD2H(ISYMB,ISYMC)
         ISYRES = MULD2H(ISINT2,ISYMT2)
         JSAIKJ = MULD2H(ISYMBC,ISYRES)
         ISYMDK = MULD2H(ISYMBC,ISINT2)
C
         LENGTH = NCKIJ(JSAIKJ)
C
         CALL DZERO(TMAT,LENGTH)
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
C           Virtual-L 1.0 contribution addWMAT
C
*     write(lupri,*)'T2TP(voo,v), isymb,b,isymc,c,isymk ',
*    * isymb,b,isymc,c,isymk,
*    * ddot(NCKI(ISYAIJ)*NVIR(ISYMD),T2TP(KOFF1),1,T2TP(KOFF1),1)
c     call output(T2TP(KOFF1),1,NCKI(ISYAIJ),1,NVIR(ISYMD),
c    * NCKI(ISYAIJ),NVIR(ISYMD),1,lupri)
*     write(lupri,*)'VLDKBC(v,o) ',
*    * ddot(NVIR(ISYMD)*NRHF(ISYMK),VLDKBC(KOFF2),1,VLDKBC(KOFF2),1)
c     call output(VLDKBC(KOFF2),1,NVIR(ISYMD),1,NRHF(ISYMK),
c    * NVIR(ISYMD),NRHF(ISYMK),1,lupri)
            CALL DGEMM('N','N',NCKI(ISYAIJ),NRHF(ISYMK),
     *                 NVIR(ISYMD),ONE,T2TP(KOFF1),NTOAIJ,
     *                 VLDKBC(KOFF2),NVIRD,ONE,
     *                 TMAT(KOFF3),NTOAIJ)
*     write(lupri,*)'TMAT(voo,o) ',
*    *ddot(NCKI(ISYAIJ)*NRHF(ISYMK),TMAT(KOFF3),1,TMAT(KOFF3),1)
c     call output(TMAT(KOFF3),1,NCKI(ISYAIJ),1,NRHF(ISYMK),
c    * NCKI(ISYAIJ),NRHF(ISYMK),1,lupri)
C
         ENDDO
C
C         CALL CC_GATHER(LENGTH,SMAT,WORK,INDSQ(1,3))
         DO I = 1,LENGTH
            WMAT(I) = WMAT(I) + TMAT(INDSQ(I,3))
         ENDDO
C
         IF (IPRINT .GT. 55) THEN
            XWMAT = DDOT(NCKIJ(JSAIKJ),WMAT,1,WMAT,1)
            WRITE(LUPRI,*) 'In WBARXBD_TMAT: vir-L1 Norm of WMAT ',XWMAT
         ENDIF


C------------------------------------------------------------
C        First virtual contribution of g term.
C
C         WMAT^BC(aikj) = WMAT^BC(aikj) - T2TP(ajkd)*g(iBd^C)
C------------------------------------------------------------
C
         ISYMBC = MULD2H(ISYMB,ISYMC)
         ISYRES = MULD2H(ISINT2,ISYMT2)
         JSAIKJ = MULD2H(ISYMBC,ISYRES)
         ISYMDI = MULD2H(ISYMBC,ISINT2)
C
         LENGTH = NCKIJ(JSAIKJ)
C
         CALL DZERO(TMAT,LENGTH)
C
         DO ISYMI = 1,NSYM
C
            ISYMD  = MULD2H(ISYMI,ISYMDI)
            ISYAJK = MULD2H(ISYMI,JSAIKJ)
C
            KOFF1 = IT2SP(ISYAJK,ISYMD)  + 1
            KOFF2 = ICKATR(ISYMDI,ISYMB) + NT1AM(ISYMDI)*(B - 1)
     *            + IT1AM(ISYMD,ISYMI)   + 1
            KOFF3 = ISAIKJ(ISYAJK,ISYMI) + 1
C
            NTOAJK = MAX(1,NCKI(ISYAJK))
            NVIRD  = MAX(NVIR(ISYMD),1)
C
C           Virtual-g 1.0 contribution addWMAT
C
            CALL DGEMM('N','N',NCKI(ISYAJK),NRHF(ISYMI),
     *                 NVIR(ISYMD),ONE,T2TP(KOFF1),NTOAJK,
     *                 VGDKCB(KOFF2),NVIRD,ONE,
     *                 TMAT(KOFF3),NTOAJK)
C
         ENDDO
C        
         DO I = 1,LENGTH
            WMAT(I) = WMAT(I) - TMAT(INDSQ(I,5))
         ENDDO
C
         IF (IPRINT .GT. 55) THEN
            XWMAT = DDOT(NCKIJ(JSAIKJ),WMAT,1,WMAT,1)
            WRITE(LUPRI,*) 'In WBARXBD_TMAT: vir-g1 Norm of WMAT ',XWMAT
         ENDIF
C
C------------------------------------------------------------
C        Second virtual contribution of L term.
C
C         WMAT^BC(aikj) = WMAT^BC(aikj) + T2TP(aikd)*L(d^CjB)
C------------------------------------------------------------
C
         ISYMBC = MULD2H(ISYMB,ISYMC)
         ISYRES = MULD2H(ISINT2,ISYMT2)
         JSAIKJ = MULD2H(ISYMBC,ISYRES)
         ISYMDJ = MULD2H(ISYMBC,ISINT2)
C
         LENGTH = NCKIJ(JSAIKJ)
C
         CALL DZERO(TMAT,LENGTH)
C
         DO ISYMJ = 1,NSYM
C
            ISYMD  = MULD2H(ISYMJ,ISYMDJ)
            ISYAIK = MULD2H(ISYMJ,JSAIKJ)
C
            KOFF1 = IT2SP(ISYAIK,ISYMD)  + 1
            KOFF2 = ICKATR(ISYMDJ,ISYMB) + NT1AM(ISYMDJ)*(B - 1)
     *            + IT1AM(ISYMD,ISYMJ)   + 1
            KOFF3 = ISAIKJ(ISYAIK,ISYMJ) + 1
C
            NTOAIK = MAX(1,NCKI(ISYAIK))
            NVIRD  = MAX(NVIR(ISYMD),1)
C
C           Virtual-L 2.0 contribution addWMAT
C
            CALL DGEMM('N','N',NCKI(ISYAIK),NRHF(ISYMJ),
     *                 NVIR(ISYMD),ONE,T2TP(KOFF1),NTOAIK,
     *                 VLDKCB(KOFF2),NVIRD,ONE,
     *                 TMAT(KOFF3),NTOAIK)
C
         ENDDO
C
         DO I = 1,LENGTH
            WMAT(I) = WMAT(I) + TMAT(I)
         ENDDO
C        
         IF (IPRINT .GT. 55) THEN
            XWMAT = DDOT(NCKIJ(JSAIKJ),WMAT,1,WMAT,1)
            WRITE(LUPRI,*) 'In WBARXBD_TMAT: vir-L2 Norm of WMAT ',XWMAT
         ENDIF
C
C------------------------------------------------------------
C        Second virtual contribution of g term.
C
C         WMAT^BC(aikj) = WMAT^BC(aikj) - T2TP(akjd)*g(iCd^B) 
C------------------------------------------------------------
C
         ISYMBC = MULD2H(ISYMB,ISYMC)
         ISYRES = MULD2H(ISINT2,ISYMT2)
         JSAIKJ = MULD2H(ISYMBC,ISYRES)
         ISYMDI = MULD2H(ISYMBC,ISINT2)
C
         LENGTH = NCKIJ(JSAIKJ)
C
         CALL DZERO(TMAT,LENGTH)
C
         DO ISYMI = 1,NSYM
C
            ISYMD  = MULD2H(ISYMI,ISYMDI)
            ISYAKJ = MULD2H(ISYMI,JSAIKJ)
C
            KOFF1 = IT2SP(ISYAKJ,ISYMD)  + 1
            KOFF2 = ICKATR(ISYMDI,ISYMB) + NT1AM(ISYMDI)*(B - 1)
     *            + IT1AM(ISYMD,ISYMI)   + 1
            KOFF3 = ISAIKJ(ISYAKJ,ISYMI) + 1
C
            NTOAKJ = MAX(1,NCKI(ISYAKJ))
            NVIRD  = MAX(NVIR(ISYMD),1)
C
C           Virtual-g 2.0 contribution addWMAT
C
            CALL DGEMM('N','N',NCKI(ISYAKJ),NRHF(ISYMI),
     *                 NVIR(ISYMD),ONE,T2TP(KOFF1),NTOAKJ,
     *                 VGDKBC(KOFF2),NVIRD,ONE,
     *                 TMAT(KOFF3),NTOAKJ)
C
         ENDDO
C
C        
         DO I = 1,LENGTH
            WMAT(I) = WMAT(I) - TMAT(INDSQ(I,2))
         ENDDO
C
         IF (IPRINT .GT. 55) THEN
            XWMAT = DDOT(NCKIJ(JSAIKJ),WMAT,1,WMAT,1)
            WRITE(LUPRI,*) 'In WBARXBD_TMAT: vir-g2 Norm of WMAT ',XWMAT
         ENDIF
C
C-------------------------------------------
C        First occupied L contribution.
C
C         WMAT^BC(aikj) = WMAT^BC(aikj)  
C                       + T2TP(ailB)*L(jl^kC) 
C                                         
C                         TB(ail)*LC(l^kj) = R(aikj)
C-------------------------------------------
C
         ISYMBC = MULD2H(ISYMB,ISYMC)
         ISYRES = MULD2H(ISINT2,ISYMT2)
         JSAIKJ = MULD2H(ISYMBC,ISYRES)
C
         ISYAIL = MULD2H(ISYMB,ISYMT2)
         ISYLKJ = MULD2H(ISYMC,ISINT2)
C
         CALL DZERO(TMAT,LENGTH)
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
C                 Occupied-L 1.0 contribution addWMAT
C
                  CALL DGEMM('N','N',NT1AM(ISYMAI),NRHF(ISYMK),
     *                       NRHF(ISYML),ONE,T2TP(KOFF1),NTOTAI,
     *                       TROCCL(KOFF2),NRHFL,ONE,TMAT(KOFF3),
     *                       NTOTAI)
C
               ENDDO
            ENDDO
         ENDDO
C
         DO I = 1,NCKIJ(JSAIKJ)
            WMAT(I) = WMAT(I) - TMAT(I)
         ENDDO
         IF (IPRINT .GT. 55) THEN
            XWMAT = DDOT(NCKIJ(JSAIKJ),WMAT,1,WMAT,1)
            WRITE(LUPRI,*) 'In WBARXBD_TMAT: occ-L1 Norm of WMAT ',XWMAT
         ENDIF
C
C-------------------------------------------
C        Second occupied L contribution.
C
C         WMAT^BC(aikj) = WMAT^BC(aikj)  
C                       + T2TP(ailC)*L(kl^jB)
C                                         
C                         TC(ail)*LB(l^jk) = R(aijk)
C
C-------------------------------------------
C

         ISYMBC = MULD2H(ISYMB,ISYMC)
         ISYRES = MULD2H(ISINT2,ISYMT2)
         JSAIKJ = MULD2H(ISYMBC,ISYRES)
C
         ISYAIL = MULD2H(ISYMC,ISYMT2)
         ISYLKJ = MULD2H(ISYMB,ISINT2)
C
         CALL DZERO(TMAT,LENGTH)
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
C                 Occupied-L 2.0 contribution addWMAT
C
                  CALL DGEMM('N','N',NT1AM(ISYMAI),NRHF(ISYMK),
     *                       NRHF(ISYML),ONE,T2TP(KOFF1),NTOTAI,
     *                       TROCCL(KOFF2),NRHFL,ONE,TMAT(KOFF3),
     *                       NTOTAI)
C
               ENDDO
            ENDDO
         ENDDO
C
         DO I = 1,NCKIJ(JSAIKJ)
            WMAT(I) = WMAT(I) - TMAT(INDSQ(I,3))
         ENDDO
         IF (IPRINT .GT. 55) THEN
            XWMAT = DDOT(NCKIJ(JSAIKJ),WMAT,1,WMAT,1)
            WRITE(LUPRI,*) 'In WBARXBD_TMAT: occ-L2 Norm of WMAT ',XWMAT
         ENDIF

C
C
C-------------------------------------------
C        First occupied g contribution.
C
C         WMAT^BC(aikj) = WMAT^BC(aikj)  
C                       - T2TP(alkB)*g(il^jC)
C                                         
C                         TB(akl)*gC(l^ji) = R(akji)
C
C-------------------------------------------
C
C
         ISYMBC = MULD2H(ISYMB,ISYMC)
         ISYRES = MULD2H(ISINT2,ISYMT2)
         JSAIKJ = MULD2H(ISYMBC,ISYRES)
C
         ISYALK = MULD2H(ISYMB,ISYMT2)
         ISYLJI = MULD2H(ISYMC,ISINT2)
C
         KALK = 1
         KEND1  = KALK   + NCKI(ISYALK)
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Not enough space in WBARXBD_TMAT (1)')
         END IF
C
         CALL DZERO(TMAT,NCKIJ(JSAIKJ))
C
C
C        T2TP(alkB) put in WORK(akl)
C
         KOFF = IT2SP(ISYALK,ISYMB) + NCKI(ISYALK)*(B - 1) + 1
         CALL CC_GATHER(NCKI(ISYALK),WORK(KALK),T2TP(KOFF),INDAJLB)
C
         DO ISYMI = 1,NSYM
C
            ISYAKJ = MULD2H(JSAIKJ,ISYMI)
            ISYMLJ = MULD2H(ISYLJI,ISYMI)
            DO I = 1,NRHF(ISYMI)
C
               DO ISYML = 1,NSYM
C
                  ISYMAK = MULD2H(ISYALK,ISYML)
                  ISYMJ  = MULD2H(ISYMLJ,ISYML)
C
                  KOFF1 = KALK
     *                  + ICKI(ISYMAK,ISYML) 
                  KOFF2 = ISJIKA(ISYLJI,ISYMC)
     *                  + NMAJIK(ISYLJI)*(C - 1)
     *                  + ISJIK(ISYMLJ,ISYMI)
     *                  + NMATIJ(ISYMLJ)*(I - 1)
     *                  + IMATIJ(ISYML,ISYMJ) + 1
                  KOFF3 = ISAIKJ(ISYAKJ,ISYMI)
     *                  + NCKI(ISYAKJ)*(I - 1)
     *                  + ICKI(ISYMAK,ISYMJ) + 1
C
                  NTOTAK = MAX(1,NT1AM(ISYMAK))
                  NRHFL  = MAX(1,NRHF(ISYML))
C
C                 Occupied-g 1.0 contribution addWMAT
C
                  CALL DGEMM('N','N',NT1AM(ISYMAK),NRHF(ISYMJ),
     *                       NRHF(ISYML),ONE,WORK(KOFF1),NTOTAK,
     *                       TROCCG(KOFF2),NRHFL,ONE,TMAT(KOFF3),
     *                       NTOTAK)
C
               ENDDO
            ENDDO
         ENDDO
C
         DO I = 1,NCKIJ(JSAIKJ)
            WMAT(I) = WMAT(I) + TMAT(INDSQ(I,2))
         ENDDO
C
         IF (IPRINT .GT. 55) THEN
            XWMAT = DDOT(NCKIJ(JSAIKJ),WMAT,1,WMAT,1)
            WRITE(LUPRI,*) 'In WBARXBD_TMAT: occ-g1 Norm of WMAT ',XWMAT
         ENDIF
C

C-------------------------------------------
C        Second occupied g contribution.
C
C         WMAT^BC(aikj) = WMAT^BC(aikj)  
C                       - T2TP(aljC)*g(il^kB)
C                                         
C                         TC(ajl)*gB(l^ki) = R(ajki)
C-------------------------------------------
C
         ISYMBC = MULD2H(ISYMB,ISYMC)
         ISYRES = MULD2H(ISINT2,ISYMT2)
         JSAIKJ = MULD2H(ISYMBC,ISYRES)
C
         ISYAJL = MULD2H(ISYMC,ISYMT2)
         ISYLKI = MULD2H(ISYMB,ISINT2)
C
         KALJ = 1
         KEND1  = KALJ   + NCKI(ISYAJL)
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Not enough space in WBARXBD_TMAT (2)')
         END IF
C
         CALL DZERO(TMAT,NCKIJ(JSAIKJ))
C
C
C        T2TP(aljC) put in WORK(ajl)
C
         KOFF = IT2SP(ISYAJL,ISYMC) + NCKI(ISYAJL)*(C - 1) + 1
         CALL CC_GATHER(NCKI(ISYAJL),WORK(KALJ),T2TP(KOFF),INDAJLC)
C
         DO ISYMI = 1,NSYM
C
            ISYAJK = MULD2H(JSAIKJ,ISYMI)
            ISYMLK = MULD2H(ISYLKI,ISYMI)
            DO I = 1,NRHF(ISYMI)
C
               DO ISYML = 1,NSYM
C
                  ISYMAJ = MULD2H(ISYAJL,ISYML)
                  ISYMK  = MULD2H(ISYMLK,ISYML)
C
                  KOFF1 = KALJ
     *                  + ICKI(ISYMAJ,ISYML) 
                  KOFF2 = ISJIKA(ISYLKI,ISYMB)
     *                  + NMAJIK(ISYLKI)*(B - 1)
     *                  + ISJIK(ISYMLK,ISYMI)
     *                  + NMATIJ(ISYMLK)*(I - 1)
     *                  + IMATIJ(ISYML,ISYMK) + 1
                  KOFF3 = ISAIKJ(ISYAJK,ISYMI)
     *                  + NCKI(ISYAJK)*(I - 1)
     *                  + ICKI(ISYMAJ,ISYMK) + 1
C
                  NTOTAJ = MAX(1,NT1AM(ISYMAJ))
                  NRHFL  = MAX(1,NRHF(ISYML))
C
C                 Occupied-g 2.0 contribution addWMAT
C
                  CALL DGEMM('N','N',NT1AM(ISYMAJ),NRHF(ISYMK),
     *                       NRHF(ISYML),ONE,WORK(KOFF1),NTOTAJ,
     *                       TROCCG(KOFF2),NRHFL,ONE,TMAT(KOFF3),
     *                       NTOTAJ)
C
               ENDDO
            ENDDO
         ENDDO
C
         DO I = 1,NCKIJ(JSAIKJ)
            WMAT(I) = WMAT(I) + TMAT(INDSQ(I,5))
         ENDDO
C
         IF (IPRINT .GT. 55) THEN
            XWMAT = DDOT(NCKIJ(JSAIKJ),WMAT,1,WMAT,1)
            WRITE(LUPRI,*) 'In WBARXBD_TMAT: occ-g2 Norm of WMAT ',XWMAT
         ENDIF
C
      ELSE IF (AIBJCK_PERM .EQ. 3) THEN
C
C--------------------------------------------------------------------------
C    
C        WMAT^BC(aikj) = WMAT^BC(aikj) + T2TP(CkjB)*F(ia) 
C                                      - T2TP(CkiB)*F(ja)
C                                      + T2TP(Ckia)*F(jB)
C                                      - T2TP(Ckja)*F(iB)
C
C--------------------------------------------------------------------------
C        Contribution from both Fock terms:
C--------------------------------------------------------------------------
C
         CALL DZERO(TMAT,LENGTH)
C
C        WMAT^BC(aikj) = WMAT^BC(aikj) + T2TP(CkjB)*F(ia) 
C        WMAT^BC(aikj) = WMAT^BC(aikj) - T2TP(CkiB)*F(ja)
C
         ISYMAI = ISYFOCK
         ISYCKJ = MULD2H(ISYMT2,ISYMB)
         ISYKJ  = MULD2H(ISYCKJ,ISYMC)
C
         DO ISYMJ = 1, NSYM
            ISYMK = MULD2H(ISYKJ,ISYMJ)
            ISYMCK = MULD2H(ISYMC,ISYMK)
            ISYAIK = MULD2H(ISYMAI,ISYMK)
            DO ISYMI = 1, NSYM
               ISYMA = MULD2H(ISYMAI,ISYMI)
C
               DO J = 1, NRHF(ISYMJ)
C
                  DO I = 1, NRHF(ISYMI)
                  DO A = 1, NVIR(ISYMA)
C
                     NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
C
C                    Index for sorted T2 amplitudes.
C
                     DO K = 1, NRHF(ISYMK)
C
                        NCK = IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K-1) +C
                        NCKJB = IT2SP(ISYCKJ,ISYMB)
     *                        + NCKI(ISYCKJ)*(B - 1)
     *                        + ICKI(ISYMCK,ISYMJ)
     *                        + NT1AM(ISYMCK)*(J - 1) + NCK
C
                        NAIKJ = ISAIKJ(ISYAIK,ISYMJ)
     *                        + NCKI(ISYAIK)*(J - 1)
     *                        + ICKI(ISYMAI,ISYMK)
     *                        + NT1AM(ISYMAI)*(K-1)
     *                        + NAI

C

                        TMAT(NAIKJ) = T2TP(NCKJB)*FOCK(NAI)
C
                     ENDDO
                  ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
C
C-------------------------------------------
C        Sum the result into WMAT.
C-------------------------------------------
C
         DO I = 1, LENGTH
C            Fock 1.0 contribution addWMAT
C            First :  WMAT^BC(aikj) = WMAT^BC(aikj) + T2TP(CkjB)*F(ia)
             WMAT(I) = WMAT(I) + TMAT(I)
C            Fock 2.0 contribution addWMAT
C            Second : WMAT^BC(aikj) = WMAT^BC(aikj) - T2TP(CkiB)*F(ja)
             WMAT(I) = WMAT(I) - TMAT(INDSQ(I,5))
         ENDDO
C
C------------------------------------------
C        Third and fourth Fock term 
C------------------------------------------
C
         CALL DZERO(TMAT,LENGTH)
C
C        WMAT^BC(aikj) = WMAT^BC(aikj) + T2TP(Ckia)*F(jB)
                                         !note that T2TP(Ckia) = T2TP(aikC)
C        WMAT^BC(aikj) = WMAT^BC(aikj) - T2TP(Ckja)*F(iB)
C                                     
         ISYMJ  = MULD2H(ISYFOCK,ISYMB)
         ISYAIK = MULD2H(ISYMT2,ISYMC)
C
         DO ISYMK = 1, NSYM
            ISYMAI = MULD2H(ISYAIK,ISYMK)
            DO ISYMI = 1, NSYM
               ISYMA = MULD2H(ISYMAI,ISYMI)
C
               DO J = 1, NRHF(ISYMJ)
                  NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J-1) + B
C
                  DO I = 1, NRHF(ISYMI)
                  DO A = 1, NVIR(ISYMA)
C 
                     NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
C
C                    Index for sorted T2 amplitudes.
C 
                     DO K = 1, NRHF(ISYMK)
C                       
                        !note that T2TP(Ckia) = T2TP(aikC)
                        NAIKC = IT2SP(ISYAIK,ISYMC)
     *                        + NCKI(ISYAIK)*(C - 1)
     *                        + ICKI(ISYMAI,ISYMK)
     *                        + NT1AM(ISYMAI)*(K - 1) + NAI
C
                        NAIKJ = ISAIKJ(ISYAIK,ISYMJ)
     *                        + NCKI(ISYAIK)*(J - 1)
     *                        + ICKI(ISYMAI,ISYMK)
     *                        + NT1AM(ISYMAI)*(K-1)
     *                        + NAI

C
                        TMAT(NAIKJ) = T2TP(NAIKC)*FOCK(NBJ)
C                    
                     ENDDO
                  ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
C
C-------------------------------------------
C        Sum the result into WMAT.
C-------------------------------------------
C
         DO I = 1, LENGTH
C            Fock 3.0 contribution addWMAT
C            First :  WMAT^BC(aikj) = WMAT^BC(aikj) + T2TP(Ckia)*F(jB)
             WMAT(I) = WMAT(I) + TMAT(I)
C            Fock 4.0 contribution addWMAT
C            Second : WMAT^BC(aikj) = WMAT^BC(aikj) - T2TP(Ckja)*F(iB)
             WMAT(I) = WMAT(I) - TMAT(INDSQ(I,5))
         ENDDO
C
C--------------------------------------------------------------------------
C     Calculate ALL virtual contributions here (in ILOOP = 1,4 loop)
C     Can be done, because: g(kad^B) = g(d^Bka), etc.
C     At the end use appropriate INDSQ.
C
C        First virtual contribution of L term.
C         WMAT^BC(aikj) = WMAT^BC(aikj) + T2TP(Ckjd)*L(d^Bia) (ILOOP = 1)
C
C        First virtual contribution of g term.
C         WMAT^BC(aikj) = WMAT^BC(aikj) - T2TP(Cjid)*g(kBd^a) (ILOOP = 2)
C
C        Second virtual contribution of L term.
C         WMAT^BC(aikj) = WMAT^BC(aikj) + T2TP(Ckid)*L(d^ajB) (ILOOP = 3)
C
C        Second virtual contribution of g term.
C         WMAT^BC(aikj) = WMAT^BC(aikj) - T2TP(Cijd)*g(kad^B) (ILOOP = 4)
C
C--------------------------------------------------------------------------
C

         ISYMBC = MULD2H(ISYMB,ISYMC)
         ISYRES = MULD2H(ISINT2,ISYMT2)
         JSAIKJ = MULD2H(ISYMBC,ISYRES)
C
         ISYAID = MULD2H(ISINT2,ISYMB)
         ISYDJK = MULD2H(ISYMT2,ISYMC)
C
         DO ILOOP = 1,4
            !sort integrals from VLDKCB(di,a) to KAID(ai,d)
            KAID = 1
            KEND1 = KAID + NCKATR(ISYAID)
            LWRK1  = LWORK  - KEND1
C
            IF (LWRK1 .LT. NCKATR(ISYAID)) THEN
               WRITE(LUPRI,*)'Memory available : ', LWRK1
               WRITE(LUPRI,*)'Memory needed    : ', NCKATR(ISYAID)
               IF (ILOOP .EQ. 1) THEN
                  CALL QUIT('Not enough space in WBARXBD_TMAT (3a)')
               ELSE IF (ILOOP .EQ. 2) THEN
                  CALL QUIT('Not enough space in WBARXBD_TMAT (3b)')
               ELSE IF (ILOOP .EQ. 3) THEN
                  CALL QUIT('Not enough space in WBARXBD_TMAT (3c)')
               ELSE IF (ILOOP .EQ. 4) THEN
                  CALL QUIT('Not enough space in WBARXBD_TMAT (3d)')
               END IF
            END IF
C
            IF (ILOOP .EQ. 1) THEN
               CALL DCOPY(NCKATR(ISYAID),VLDKCB,1,WORK(KAID),1)
            ELSE IF (ILOOP .EQ. 2) THEN
               CALL DCOPY(NCKATR(ISYAID),VGDKBC,1,WORK(KAID),1)
            ELSE IF (ILOOP .EQ. 3) THEN
               CALL DCOPY(NCKATR(ISYAID),VLDKBC,1,WORK(KAID),1)
            ELSE IF (ILOOP .EQ. 4) THEN
               CALL DCOPY(NCKATR(ISYAID),VGDKCB,1,WORK(KAID),1)
            END IF
            CALL CCSDT_SRVIR3(WORK(KAID),WORK(KEND1),ISYMB,B,ISINT2)
C
            LENGTH = NCKIJ(JSAIKJ)
C
            CALL DZERO(TMAT,LENGTH)
C
            DO ISYMK = 1,NSYM
               ISYDJ = MULD2H(ISYDJK,ISYMK)
               DO ISYMJ = 1,NSYM
                  ISYMD = MULD2H(ISYDJ,ISYMJ)
                  ISYAI = MULD2H(ISYAID,ISYMD)
                  ISYAIJ = MULD2H(ISYAI,ISYMJ)
                  DO K = 1,NRHF(ISYMK)
C
                     KOFF1 = KAID
     *                     + ICKATR(ISYAI,ISYMD)  
                     KOFF2 = IT2SP(ISYDJK,ISYMC) 
     *                     + NCKI(ISYDJK)*(C-1)
     *                     + ICKI(ISYDJ,ISYMK)
     *                     + NT1AM(ISYDJ)*(K-1)
     *                     + IT1AM(ISYMD,ISYMJ)
     *                     + 1
                     KOFF3 = ISAIKJ(ISYAIJ,ISYMK) 
     *                     + NCKI(ISYAIJ)*(K-1)
     *                     + ICKI(ISYAI,ISYMJ)
     *                     + 1
C
                     NTOTAI = MAX(1,NT1AM(ISYAI))
                     NVIRD  = MAX(NVIR(ISYMD),1)
C                                
C                    Virtual-L 1.0 contribution addWMAT
C
                     CALL DGEMM('N','N',NT1AM(ISYAI),NRHF(ISYMJ),
     *                          NVIR(ISYMD),ONE,WORK(KOFF1),NTOTAI,
     *                          T2TP(KOFF2),NVIRD,ONE,
     *                          TMAT(KOFF3),NTOTAI)
C
                  ENDDO
               ENDDO
            ENDDO
C
C            CALL CC_GATHER(LENGTH,SMAT,WORK,INDSQ(1,3))
            DO I = 1,LENGTH
               IF (ILOOP .EQ. 1) THEN   
                  WMAT(I) = WMAT(I) + TMAT(INDSQ(I,3))
               ELSE IF (ILOOP .EQ. 2) THEN
                  WMAT(I) = WMAT(I) - TMAT(INDSQ(I,1))
               ELSE IF (ILOOP .EQ. 3) THEN
                  WMAT(I) = WMAT(I) + TMAT(INDSQ(I,4))
               ELSE IF (ILOOP .EQ. 4) THEN
                  WMAT(I) = WMAT(I) - TMAT(INDSQ(I,2))
               END IF
            ENDDO
C
            IF (IPRINT .GT. 55) THEN
               XWMAT = DDOT(NCKIJ(JSAIKJ),WMAT,1,WMAT,1)
               IF (ILOOP .EQ. 1) THEN
                 WRITE(LUPRI,*) 'WBARXBD_TMAT virL1 Norm of WMAT ',XWMAT
               ELSE IF (ILOOP .EQ. 2) THEN
                 WRITE(LUPRI,*) 'WBARXBD_TMAT virg2 Norm of WMAT ',XWMAT
               ELSE IF (ILOOP .EQ. 3) THEN
                 WRITE(LUPRI,*) 'WBARXBD_TMAT virL2 Norm of WMAT ',XWMAT
               ELSE IF (ILOOP .EQ. 4) THEN
                 WRITE(LUPRI,*) 'WBARXBD_TMAT virg2 Norm of WMAT ',XWMAT
               END IF
            ENDIF
C
         END DO ! ILOOP
C
C---------------------------------------------------------------------
C        First occupied L contribution.
C         WMAT^BC(aikj) = WMAT^BC(aikj) - T2TP(CklB)*L(jl^ia) (ILOOP = 1)
C
C        First occupied g contribution.
C         WMAT^BC(aikj) = WMAT^BC(aikj) + T2TP(BilC)*g(kl^ja) (ILOOP = 2)
C
C---------------------------------------------------------------------
C
         ISYMBC = MULD2H(ISYMB,ISYMC)
         ISYRES = MULD2H(ISINT2,ISYMT2)
         JSAIKJ = MULD2H(ISYMBC,ISYRES)
C
         !T2TP(BlkC) = T2TP(CklB)
         ISYLK = MULD2H(ISYMT2,ISYMBC)
C
         DO ILOOP = 1,2
C
            KLK = 1
            KAJIL = KLK + NMATIJ(ISYLK)
            KEND1 = KAJIL + NTRAOC(ISINT2)
            LWRK1  = LWORK  - KEND1
C
            IF (LWRK1 .LT. 0) THEN
               WRITE(LUPRI,*)'Memory available : ', LWORK
               WRITE(LUPRI,*)'Memory needed    : ', KEND1
               IF (ILOOP .EQ. 1) THEN
                  CALL QUIT('Not enough space in WBARXBD_TMAT (4a)')
               ELSE IF (ILOOP .EQ. 2) THEN
                  CALL QUIT('Not enough space in WBARXBD_TMAT (4b)')
               END IF
            END IF
C
            CALL DZERO(TMAT,LENGTH)
C
            !sort from T2TP(BlkC) to KLK(lk)
            IF (ILOOP .EQ. 1) THEN
               CALL SORT_T2_IJ(WORK(KLK),ISYMB,B,ISYMC,C,T2TP,ISYMT2)
            ELSE IF (ILOOP .EQ. 2) THEN
               CALL SORT_T2_IJ(WORK(KLK),ISYMC,C,ISYMB,B,T2TP,ISYMT2)
            END IF
            !sort from TROCCL(lij,a) to KAJIL(ajil)
            IF (ILOOP .EQ. 1) THEN
               CALL CCFOP_SORT(TROCCL,WORK(KAJIL),ISINT2,1)
            ELSE IF (ILOOP .EQ. 2) THEN
               CALL CCFOP_SORT(TROCCG,WORK(KAJIL),ISINT2,1)
            END IF
C
            !multiply KAJIL(aji,l)*KLK(l,k) --> TMAT(aji,k) (ILOOP = 1)
            !multiply KAJIL(akj,l)*KLK(l,i) --> TMAT(akj,i) (ILOOP = 2)
C
            DO ISYML = 1,NSYM
               ISYMK  = MULD2H(ISYLK,ISYML)
               ISYAJI = MULD2H(JSAIKJ,ISYMK)
C
               KOFF1 = KAJIL
     *               + ISAIKJ(ISYAJI,ISYML)
C 
               KOFF2 = KLK
     *               + IMATIJ(ISYML,ISYMK)
C
               KOFF3 = ISAIKJ(ISYAJI,ISYMK)
     *               + 1
C
               NTOTAJI = MAX(1,NCKI(ISYAJI))
               NTOTL  = MAX(NRHF(ISYML),1)
C
               CALL DGEMM('N','N',NCKI(ISYAJI),NRHF(ISYMK),
     *                    NRHF(ISYML),ONE,WORK(KOFF1),NTOTAJI,
     *                    WORK(KOFF2),NTOTL,ONE,
     *                    TMAT(KOFF3),NTOTAJI)
C
            END DO
C
            DO I = 1,LENGTH
               IF (ILOOP .EQ. 1) THEN
                  !TMAT(ajik) --> WMAT(aikj)
                  WMAT(I) = WMAT(I) - TMAT(INDSQ(I,4))
               ELSE IF (ILOOP .EQ. 2) THEN
                  !TMAT(akji) --> WMAT(aikj)
                  WMAT(I) = WMAT(I) + TMAT(INDSQ(I,2))
               END IF
            END DO
C
         END DO !ILOOP
C
C------------------------------------------------------------
C        Second occupied L contribution.
C         WMAT^BC(aikj) = WMAT^BC(aikj) - T2TP(Ckla)*L(il^jB) 
C
C------------------------------------------------------------
C

C
C------------------------------------------------------------
C        Second occupied g contribution.
C         WMAT^BC(aikj) = WMAT^BC(aikj) + T2TP(ajlC)*g(kl^iB) 
C
C        This is part of code is actually strongly based on 
C         the second occupied L contribution in "AIBJCK_PERM = 1" 
C         part; this means that the nomenclature might be a bit
C         confusing, because:
C
C         i <--> j
C
C------------------------------------------------------------
C
         ISYMBC = MULD2H(ISYMB,ISYMC)
         ISYRES = MULD2H(ISINT2,ISYMT2)
         JSAIKJ = MULD2H(ISYMBC,ISYRES)
C
         ISYAIL = MULD2H(ISYMC,ISYMT2)
         ISYLKJ = MULD2H(ISYMB,ISINT2)
C
         CALL DZERO(TMAT,LENGTH)
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
     *                       NRHF(ISYML),ONE,T2TP(KOFF1),NTOTAI,
     *                       TROCCG(KOFF2),NRHFL,ONE,TMAT(KOFF3),
     *                       NTOTAI)
C
               ENDDO
            ENDDO
         ENDDO
C
         !TMAT(ajik) --> WMAT(aikj)
         DO I = 1,NCKIJ(JSAIKJ)
            WMAT(I) = WMAT(I) + TMAT(INDSQ(I,4))
         ENDDO
C
C------------------------------------------------------------
C        Second occupied L contribution.
C         WMAT^BC(aikj) = WMAT^BC(aikj) - T2TP(alkC)*L(il^jB)
C
C        This is part of code is actually strongly based on 
C         the second occupied g contribution in "AIBJCK_PERM = 1" 
C         part; this means that the nomenclature might be a bit
C         confusing, because:
C
C         j <--> k
C
C------------------------------------------------------------
C
         ISYMBC = MULD2H(ISYMB,ISYMC)
         ISYRES = MULD2H(ISINT2,ISYMT2)
         JSAIKJ = MULD2H(ISYMBC,ISYRES)
C
         ISYAJL = MULD2H(ISYMC,ISYMT2)
         ISYLKI = MULD2H(ISYMB,ISINT2)
C
         KALJ = 1
         KEND1  = KALJ   + NCKI(ISYAJL)
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Not enough space in WBARXBD_TMAT (5)')
         END IF
C
         CALL DZERO(TMAT,NCKIJ(JSAIKJ))
C
C
C        T2TP(aljC) put in WORK(ajl)
C
         KOFF = IT2SP(ISYAJL,ISYMC) + NCKI(ISYAJL)*(C - 1) + 1
         CALL CC_GATHER(NCKI(ISYAJL),WORK(KALJ),T2TP(KOFF),INDAJLC)
C
         DO ISYMI = 1,NSYM
C
            ISYAJK = MULD2H(JSAIKJ,ISYMI)
            ISYMLK = MULD2H(ISYLKI,ISYMI)
            DO I = 1,NRHF(ISYMI)
C
               DO ISYML = 1,NSYM
C
                  ISYMAJ = MULD2H(ISYAJL,ISYML)
                  ISYMK  = MULD2H(ISYMLK,ISYML)
C
                  KOFF1 = KALJ
     *                  + ICKI(ISYMAJ,ISYML) 
                  KOFF2 = ISJIKA(ISYLKI,ISYMB)
     *                  + NMAJIK(ISYLKI)*(B - 1)
     *                  + ISJIK(ISYMLK,ISYMI)
     *                  + NMATIJ(ISYMLK)*(I - 1)
     *                  + IMATIJ(ISYML,ISYMK) + 1
                  KOFF3 = ISAIKJ(ISYAJK,ISYMI)
     *                  + NCKI(ISYAJK)*(I - 1)
     *                  + ICKI(ISYMAJ,ISYMK) + 1
C
                  NTOTAJ = MAX(1,NT1AM(ISYMAJ))
                  NRHFL  = MAX(1,NRHF(ISYML))
C
                  CALL DGEMM('N','N',NT1AM(ISYMAJ),NRHF(ISYMK),
     *                       NRHF(ISYML),ONE,WORK(KOFF1),NTOTAJ,
     *                       TROCCL(KOFF2),NRHFL,ONE,TMAT(KOFF3),
     *                       NTOTAJ)
C
               ENDDO
            ENDDO
         ENDDO
C
         DO I = 1,NCKIJ(JSAIKJ)
            WMAT(I) = WMAT(I) - TMAT(INDSQ(I,2))
         ENDDO
C
      ELSE
         WRITE(LUPRI,*) 'AIBJCK_PERM = ', AIBJCK_PERM
         WRITE(LUPRI,*) 'WBARXBD_T2 works for AIBJCK_PERM = 1 or 3'
         CALL QUIT('Wrong AIBJCK_PERM option in WBARXBD_TMAT')
      END IF

C
      CALL QEXIT('WBXTMT')
C
      RETURN
      END
C  /* Deck wx_bd_o */
      SUBROUTINE WX_BD_O(AIBJCK_PERM,LW,LWBAR,TMAT,ISTMAT,FOCKY,ISYFKY,
     *                 WMAT,ISWMAT,WRK,LWRK)
* 
* If (AIBJCK_PERM.eq.1) then (bjdk) permutation symmetry 
*
*     WBD(aikj) = WBD(aikj) - t(aBD,ljk) * fock(li)
*
*                            tmatBD(alkj)
*
* else if (AIBJCK_PERM.eq.2) then (aidk) permutation symmetry
*
*     WBD(aikj) = WBD(aikj) - t(aBD,ilk) * fock(lj)
*
*                            tmatBD(aikl)
*
* else if (AIBJCK_PERM.eq.3) then (aibj) permutation symmetry
*
*     WBD(aikj) = WBD(aikj) - t(aBD,ijl) * fock(lk)
*
*                            tmatBD(ailj)
*
* else if (AIBJCK_PERM.eq.4) then calculate all terms
* 
*
* Written by P. Jorgensen and F. Pawlowski, Spring 2002.
* (modyfied for AIBJCK_PERM flag - spring 2003.)
*
* Autumn 2003, F. Pawlowski:
*
* Generalized to treat either the triples amplitudes 
* (LW = .TRUE., LWBAR = .FALSE.) or 
* the triplees multipliers (LW = .FALSE., LWBAR = .TRUE.).
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
C
      LOGICAL LW,LWBAR
C
      INTEGER AIBJCK_PERM
C
      INTEGER LWRK, KFCLI, KEND0, LWRK0, KOFF1, KOFF2 
      INTEGER NL, KOFFY, KOFFT, KOFFW 
      INTEGER ISTMAT, ISYFKY, ISWMAT, ISALKJ
      INTEGER ISYMA, ISYAI, ISYAIK, ISYALK, ISYAL, NA
      INTEGER ISYMJ, ISYMK, ISYMI, ISYML, ISYFI
      INTEGER ISYAIL,NAI,NAIK
      INTEGER NI,NJ,NK
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
      PARAMETER (HALF = 0.5D0, ONE = 1.0D0)
C
      CALL QENTER('WX_BDO')
C
C---------------------------------------
C     Initial test of AIBJCK_PERM option
C---------------------------------------
C
      IF ( (AIBJCK_PERM .LT. 1) .OR. (AIBJCK_PERM .GT. 4) ) THEN
         WRITE(LUPRI,*)'AIBJCK_PERM = ',AIBJCK_PERM
         WRITE(LUPRI,*)'should be between 1 and 4 '
         CALL QUIT('Illegal value of AIBJCK_PERM option in WX_BD_O')
      END IF
C
C Initial test of logical flags
C
      IF (LW .AND. .NOT.LWBAR) THEN
         CONTINUE
      ELSE IF (.NOT.LW .AND. LWBAR) THEN
         CONTINUE
      ELSE
         WRITE(LUPRI,*) 'LW = ', LW
         WRITE(LUPRI,*) 'LWBAR = ', LWBAR
         WRITE(LUPRI,*) 'LW and LWBAR flags must have opposite values '
         CALL QUIT('Logic inconsistency in WX_BD_O')
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
         CALL QUIT('Insufficient space in WX_BD_O')
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
C IF (LW) THEN
C    WBD(a,i,k,j) = WBD(a,i-,k,j) - sum (l) tmatBD(a,l,k,j)*focky(l,i)
C ELSE
C    WBD(a,i,k,j) = WBD(a,i-,k,j) - sum (l) tmatBD(a,l,k,j)*focky(i,l)
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
C
                        IF (LW) THEN
                           NL    = MAX(1,NRHF(ISYML))
                        ELSE
                           NI    = MAX(1,NRHF(ISYMI))
                        END IF
C
                        IF (LW) THEN
                           KOFFY = KFCLI + IMATIJ(ISYML,ISYMI) 
                        ELSE
                           KOFFY = KFCLI + IMATIJ(ISYMI,ISYML)
                        END IF
C
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
                        IF (LW) THEN
                           CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),
     *                                NRHF(ISYML),ONE,TMAT(KOFFT),NA,
     *                                WRK(KOFFY),NL,ONE,WMAT(KOFFW),NA)
                        ELSE
                           CALL DGEMM('N','T',NVIR(ISYMA),NRHF(ISYMI),
     *                                NRHF(ISYML),ONE,TMAT(KOFFT),NA,
     *                                WRK(KOFFY),NI,ONE,WMAT(KOFFW),NA)
                        END IF
C
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
C IF (LW) THEN
C    WBD(a,i,k,j) = WBD(a,i,k,j-) - sum (l) tmatBD(a,i,k,l)*focky(l,j)
C ELSE
C    WBD(a,i,k,j) = WBD(a,i,k,j-) - sum (l) tmatBD(a,i,k,l)*focky(j,l)
C
         DO ISYMJ = 1,NSYM
            ISYML = MULD2H(ISYFKY,ISYMJ)
            ISYAIK =MULD2H(ISTMAT,ISYML)
            NAIK    = MAX(1,NCKI(ISYAIK))
C
            IF (LW) THEN
               NL    = MAX(1,NRHF(ISYML))
            ELSE 
               NJ    = MAX(1,NRHF(ISYMJ))
            END IF
C
            IF (LW) THEN
               KOFFY = KFCLI + IMATIJ(ISYML,ISYMJ)
            ELSE
               KOFFY = KFCLI + IMATIJ(ISYMJ,ISYML)
            END IF
C
            KOFFT = ISAIKJ(ISYAIK,ISYML)
     *            + 1
            KOFFW = ISAIKJ(ISYAIK,ISYMJ)
     *            + 1
C
C DENOTE t3 IS CALCULATED WITH NEGATIVE SIGN
C
            IF (LW) THEN
               CALL DGEMM('N','N',NCKI(ISYAIK),NRHF(ISYMJ),
     *                    NRHF(ISYML),ONE,TMAT(KOFFT),NAIK,
     *                    WRK(KOFFY),NL,ONE,WMAT(KOFFW),NAIK)
            ELSE
               CALL DGEMM('N','T',NCKI(ISYAIK),NRHF(ISYMJ),
     *                    NRHF(ISYML),ONE,TMAT(KOFFT),NAIK,
     *                    WRK(KOFFY),NJ,ONE,WMAT(KOFFW),NAIK)
            END IF               
C
         END DO

      END IF
      IF ((AIBJCK_PERM.EQ.3) .OR. (AIBJCK_PERM.EQ.4)) THEN
C
C CARRY OUT MATRIX MULTIPLICATION
C IF (LW) THEN
C    WBD(a,i,k,j) = WBD(a,i,k-,j) - sum (l) tmatBD(a,i,l,j)*focky(l,k)
C ELSE
C    WBD(a,i,k,j) = WBD(a,i,k-,j) - sum (l) tmatBD(a,i,l,j)*focky(k,l)
C 
         DO ISYMJ = 1,NSYM
            ISYAIL =MULD2H(ISTMAT,ISYMJ)
            ISYAIK = MULD2H(ISWMAT,ISYMJ)
            DO J = 1,NRHF(ISYMJ)
               DO ISYMK = 1,NSYM
                  ISYAI = MULD2H(ISYAIK,ISYMK)
                  ISYML = MULD2H(ISYFKY,ISYMK)
                     NAI    = MAX(1,NT1AM(ISYAI))
C
                     IF (LW) THEN
                        NL    = MAX(1,NRHF(ISYML))
                     ELSE
                        NK    = MAX(1,NRHF(ISYMK))
                     END IF
C
                     IF (LW) THEN
                        KOFFY = KFCLI + IMATIJ(ISYML,ISYMK) 
                     ELSE
                        KOFFY = KFCLI + IMATIJ(ISYMK,ISYML)
                     END IF
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
                     IF (LW) THEN
                        CALL DGEMM('N','N',NT1AM(ISYAI),NRHF(ISYMK),
     *                             NRHF(ISYML),ONE,TMAT(KOFFT),NAI,
     *                             WRK(KOFFY),NL,ONE,WMAT(KOFFW),NAI)
                     ELSE
                        CALL DGEMM('N','T',NT1AM(ISYAI),NRHF(ISYMK),
     *                             NRHF(ISYML),ONE,TMAT(KOFFT),NAI,
     *                             WRK(KOFFY),NK,ONE,WMAT(KOFFW),NAI)
                     END IF
C
               END DO   
            END DO   
         END DO   
C
      END IF
C
      CALL QEXIT('WX_BDO')
C
      RETURN
      END
C  /* Deck cc3_xi_den_abij_cub */
      SUBROUTINE CC_XI_DEN_ABIJ_CUB(CUBIC,LISTL,LISTRZ,LISTRU,
     *                           DAB,DIJ,DAI,ISYDEN,
     *                           L2L1,ISYML1,
     *                           ISYFCKX,FOCKX,
     *                           ISYFCKY,FOCKY,
     *                           ISYMT3,ISWMAT,ISTHETA,
     *                           LUT3,FNT3,LUWBMAT,FNWBMAT,
     *                           LUTHETA,FNTHETA,
     *                           LUWZU,FNWZU,
     *                           LUWBZU,FNWBZU,
     *                           FOCKD,FREQX,FREQY,
     *                           WORK,LWORK,ISYMD,D)
C
C=========================================================================
C
C    CUBIC has to be .TRUE.  for CUBIC response calculations
C=========================================================================
C
C     Dab, Dij and Dai densities for cc3 cubic response ( A{Y} matrix).
C
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <cc3t3d.h>
C
      CHARACTER LISTRZ*3,LISTRU*3,LISTL*3
C
      INTEGER ISYMT3,ISWMAT,LUT3,LUWBMAT,LWORK,ISYMD,ISYDEN
      INTEGER ISYML,ISYMDL,ISWMATDL,ISYMT3DL,ISYMN,ISYEMF,ISYMBN,ISYMEM
      INTEGER ISYMB,ISYMF,ISYMFI,ISYMI,ISYEMB,ISYMFN
      INTEGER KT3,KWMAT,KEND1,LWRK1
      INTEGER KOFF1,KOFF2,KOFF3,KBN,KFN
      INTEGER NTOTEM,NTOTF,NNEMF
      INTEGER IADR
C
      INTEGER ISTHETA,ISYFCKX,ISYFCKY,LUTHETA,LUWZU,LUWBZU
      INTEGER ISTHETADL,ISTHETAFX,ISTHETAFY
      INTEGER KTHETA,KTHETAF,KWZU
      INTEGER KFI
      INTEGER IOPT
      INTEGER MAXX1
C
      INTEGER KWBZU
C
      INTEGER ISYMJ,ISYMFJ,KFJ
      INTEGER ISYMM,ISYME
C
      INTEGER ISYMDAI,ISYML1
C
      LOGICAL CUBIC
      LOGICAL TRANSPOSEW
C
      CHARACTER*(*) FNT3,FNWBMAT,FNTHETA,FNWZU,FNWBZU
C
#if defined (SYS_CRAY)
      REAL DAB(*),DIJ(*),DAI(*),WORK(LWORK),ONE,HALF
      REAL FOCKX(*),FOCKY(*),L2L1(*),FOCKD(*),FREQX,FREQY
      REAL XNORMVAL,DDOT,FREQXY
#else
      DOUBLE PRECISION DAB(*),DIJ(*),DAI(*),WORK(LWORK),ONE,HALF
      DOUBLE PRECISION FOCKX(*),FOCKY(*),L2L1(*),FOCKD(*),FREQX,FREQY
      DOUBLE PRECISION XNORMVAL,DDOT,FREQXY
#endif
C
      PARAMETER(ONE = 1.0D0, HALF = 0.5D0)
C
      CALL QENTER('DENABIJC')
C
      ISYMDAI = MULD2H(ISTHETA,ISYML1)
      !symmetry check
      IF (ISYMDAI .NE. ISYDEN) THEN
         WRITE(LUPRI,*)'ISYMDAI ', ISYMDAI
         WRITE(LUPRI,*)'ISYDEN ', ISYDEN
         WRITE(LUPRI,*)'These symmetries should be the same '
         CALL QUIT('Symmetry inconsistency in CC_XI_DEN_ABIJ_CUB')
      END IF
C
      DO ISYML = 1,NSYM
C
         ISYMDL = MULD2H(ISYMD,ISYML)
         ISWMATDL = MULD2H(ISWMAT,ISYMDL)
         ISYMT3DL = MULD2H(ISYMT3,ISYMDL)
         ISTHETADL = MULD2H(ISTHETA,ISYMDL)
         IF (LISTRU(1:3).EQ.'R1 ') THEN
            ISTHETAFX  = MULD2H(ISYMT3DL,ISYFCKX)
            ISTHETAFY  = MULD2H(ISYMT3DL,ISYFCKY)
         END IF
C
         KT3  = 1
         KWMAT  = KT3 + NT2SQ(ISYMT3DL)
         KWBZU = KWMAT + NT2SQ(ISWMATDL)
         KEND1 = KWBZU + NT2SQ(ISWMATDL)
         LWRK1  = LWORK - KEND1
C
         IF (CUBIC) THEN
C
            MAXX1 = 0
            IF (LISTRU(1:3).EQ.'R1 ') THEN
               MAXX1 = MAX(NT2SQ(ISTHETAFX),NT2SQ(ISTHETAFY))
            END IF
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
C
C           ------------------------------------------------
C           Read wMAT_bar from file 
C           ------------------------------------------------
C
            IADR = ISWTL(ISWMATDL,ISYML) + NT2SQ(ISWMATDL)*(L-1) + 1
            CALL GETWA2(LUWBZU,FNWBZU,WORK(KWBZU),IADR,
     *                  NT2SQ(ISWMATDL))

            IF (LISTRU(1:3).EQ.'R1 ') THEN
C              ---------------------------------------------
C              4ht line of Eq. 62 (second cont)
C              ---------------------------------------------

C
C              KTHETAF(De- f)_(lmi) = KT3 * FOCKX
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
               CALL DIJ_CONT_CUB(TRANSPOSEW,DIJ,ONE,WORK(KWBZU),
     *                         ISWMATDL,
     *                         WORK(KTHETA),ISTHETADL,WORK(KEND1),LWRK1)

            END IF ! LISTRU .EQ. R1
C
C
C           ---------------------------------------------
C           3rd line of Eq. (62)
C           ---------------------------------------------
 
C           ------------------------------------------------
C           Read WMAT_bar from file 
C           ------------------------------------------------
C
            IADR = ISWTL(ISWMATDL,ISYML) + NT2SQ(ISWMATDL)*(L-1) + 1
            CALL GETWA2(LUWBMAT,FNWBMAT,WORK(KWMAT),IADR,
     *                  NT2SQ(ISWMATDL))
C
            IF (LISTRU(1:3).EQ.'R1 ') THEN
C
C              KTHETAF(Def- )_(lmi) = KT3 * FOCKX
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
C
C              contract... (3rd line of Eq. (62))
C

               TRANSPOSEW = .FALSE.
               CALL DIJ_CONT_CUB(TRANSPOSEW,DIJ,ONE,WORK(KWMAT),
     *                         ISWMATDL,
     *                         WORK(KTHETA),ISTHETADL,WORK(KEND1),LWRK1)
C
            END IF ! LISTRU .EQ. R1
                
C
C           ----------------------------------------------------
C           Read THETA(Deb)_(l-- m-- n--)  amplitudes from file 
C           and symmetrize them:
C           ----------------------------------------------------
C
            IADR = ISWTL(ISTHETADL,ISYML) + NT2SQ(ISTHETADL)*(L-1) 
     *           + 1
            CALL GETWA2(LUTHETA,FNTHETA,WORK(KTHETA),IADR,
     *                  NT2SQ(ISTHETADL))
C
            CALL CC_T2MOD(WORK(KTHETA),ISTHETADL,ONE)

C           ----------------------------------------------------
C           Read wZU^{Deb-}_{l- m- n-} from file...
C           ----------------------------------------------------
C
            IADR = ISWTL(ISTHETADL,ISYML) + NT2SQ(ISTHETADL)*(L-1) 
     *           + 1
            CALL GETWA2(LUWZU,FNWZU,WORK(KWZU),IADR,
     *                  NT2SQ(ISTHETADL))
            !second contribution to Dab (second line in Eq. (61))

            CALL DAXPY(NT2SQ(ISWMATDL),HALF,WORK(KWMAT),1,
     *                 WORK(KWBZU),1)
            TRANSPOSEW = .TRUE.
            CALL DAB_CONT_CUB(TRANSPOSEW,DAB,WORK(KWBZU),ISWMATDL,
     *                        WORK(KWZU),ISTHETADL,WORK(KEND1),LWRK1)

            !4th line in Eq. (62) (Dij) (first cont)
            TRANSPOSEW = .TRUE.
            CALL DIJ_CONT_CUB(TRANSPOSEW,DIJ,HALF,WORK(KWMAT),ISWMATDL,
     *                        WORK(KWZU),ISTHETADL,WORK(KEND1),LWRK1)


C
C          ---------------------------------------------------------
C          DAI(ai) = DAI(ai) +  L2L1{emld}*(w{Aed-}_{i-m-l-} - (w{Aed-}_{m-i-l-}
C          ---------------------------------------------------------
C
            CALL ADEN_DAI_T2_D_CUB(DAI,ISYMDAI,L2L1,ISYML1,
     *                         WORK(KWZU),ISTHETADL,ISYMD,D,
     *                         ISYML,L,WORK(KEND1),LWRK1)
C
C           -----------------------------------------------------------
C           ...and create wZU^{Deb-}_{l- m- n-} + wZU^{Dbe-}_{l- n- m-}
C           -----------------------------------------------------------

            CALL CC_T2MOD(WORK(KWZU),ISTHETADL,ONE)

C           ---------------------------------------------------------
C           Get THETA + wZU^{Deb-}_{l- m- n-} + wZU^{Dbe-}_{l- n- m-}
C           ---------------------------------------------------------

            CALL DAXPY(NT2SQ(ISTHETADL),ONE,WORK(KWZU),1,
     *                    WORK(KTHETA),1)

C
C-----------------------------------------------------------------------
C  DAI(ai) = DAI(ai) +  L2L1{emLD}*(THETA{Dea}_{Lmi} - THETA{Dea}_{Lim})
C-----------------------------------------------------------------------
C
C
            CALL ADEN_DAI_T2_D(DAI,ISYMDAI,L2L1,ISYML1,
     *                         WORK(KTHETA),ISTHETADL,ISYMD,D,
     *                         ISYML,L,WORK(KEND1),LWRK1)
C

            !generate WMAT-tilde:
            CALL CC_T2MOD(WORK(KWMAT),ISWMATDL,HALF)
C
C           -------------------------------------------------------
C           D(fb) <- D(fb)+ sum_em Wtilde_bar^DL(em,fN) T3^DL(em,bN):
            ! FOR CUBIC = .TRUE. T3^DL(em,bN) becomes THETA_Z^DL(em,bN)
C           -------------------------------------------------------

            TRANSPOSEW = .FALSE.
            CALL DAB_CONT_CUB(TRANSPOSEW,DAB,WORK(KWMAT),ISWMATDL,
     *                        WORK(KTHETA),ISTHETADL,WORK(KEND1),LWRK1)
C
C
c             -------------------------------------------------------
c             D(iN) <- D(iN)- sum_emf Wtilde_bar^DL(em,fN) t^DL(em,fi):
c           ! FOR CUBIC = .TRUE. t^DL(em,fi) becomes THETA_Z^DL(em,fi)
c             -------------------------------------------------------

            TRANSPOSEW = .FALSE.
            CALL DIJ_CONT_CUB(TRANSPOSEW,DIJ,ONE,WORK(KWMAT),ISWMATDL,
     *                        WORK(KTHETA),ISTHETADL,WORK(KEND1),LWRK1)
C
       ! Calculate the extra contribution to D(ij) density:
       ! D(ij) = D(ij) + W^Df(emlj) * [theta^(Def-)_(iml) + theta^(De-f)_(iml)]
C
C  ----------------------------
C  Read T3^DL(em,fi) amplitudes
C  ----------------------------
C
            ! KT3 is recycled here
            CALL READ_T3_AIBL(LUT3,FNT3,ISYMT3,WORK(KT3),
     *                        ISYMT3DL,L,ISYML,ISYMD)
C
C ----------------------------------------------
C Contract T3^DL(em,fi) with X operator 
C to get THDL(em,fi) = [ THETA^DL(em,f-i) + THETA^DL(e-m,fi) ]
C ----------------------------------------------
C
cNow it  becomes the second line of Eq (62) (part of it).
c

            IF (LISTRU(1:3).EQ.'R1 ') THEN
C
C              (1) KTHETAF(X) = KT3 * FOCKX
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
            ELSE IF (LISTRU(1:3).EQ.'RE ') THEN
               CALL DZERO(WORK(KTHETA),NT2SQ(ISTHETADL))
            END IF
C
c           Now we construct wXY(Def-)_(i- m- l-)

            CALL DZERO(WORK(KWZU),NT2SQ(ISTHETADL))
            CALL DZERO(WORK(KTHETAF),NT2SQ(ISTHETADL))
C
            CALL READ_T3_ALBJ(LUWZU,FNWZU,ISTHETA,WORK(KWZU),
     *                        ISTHETADL,L,ISYML,ISYMD)
c
            !transpose and accumalte 
            CALL TRANS_AIBJ_BJAI(WORK(KWZU),WORK(KTHETAF),ISTHETADL)
            CALL DAXPY(NT2SQ(ISTHETADL),ONE,WORK(KTHETAF),1,
     *                 WORK(KTHETA),1)
C
C ------------------------------------------------
C Read WBMAT^DL(em,fj) from the file
C ------------------------------------------------
C
            ! KWBMAT is recycled here
            CALL READ_T3_AIBL(LUWBMAT,FNWBMAT,ISWMAT,WORK(KWMAT),
     *                        ISWMATDL,L,ISYML,ISYMD)

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
     *                             + NT1AM(ISYMEM)*IT1AM(ISYMF,ISYMI)

                     KFJ    = IT1AM(ISYMF,ISYMJ)+NVIR(ISYMF)*(J-1)+1
                     KOFF2  = KWMAT + IT2SQ(ISYMEM,ISYMFJ)
     *                           + NT1AM(ISYMEM)*(KFJ-1)

                     KOFF3  = IMATIJ(ISYMI,ISYMJ)
     *                           + NRHF(ISYMI)*(J-1) + 1

                     NNEMF  = MAX(NT1AM(ISYMEM)*NVIR(ISYMF),1)
C
                     CALL DGEMV('T',NT1AM(ISYMEM)*NVIR(ISYMF),
     *                         NRHF(ISYMI),-ONE,WORK(KOFF1),NNEMF,
     *                         WORK(KOFF2),1,ONE,DIJ(KOFF3),1)
C

                  END DO ! ISYMEM

               END DO ! J 
            END DO    ! ISYMJ 
C
C
C
         END DO       ! L 
      END DO          ! ISYML 
C
      CALL QEXIT('DENABIJC')
C
      RETURN
      END
C  /* Deck dij_cont_cub */
      SUBROUTINE DIJ_CONT_CUB(TRANSPOSEW,DIJ,FACTOR,WBARDL,ISWMATDL,
     *                        THETADL,ISTHETADL,WORK,LWORK)
*
**********************************************************************
*
* Calculate the contribution to the DIJ density (cubic response) of 
* following type:
*
* FACTOR*Wbar^{Df}(emjl) * theta^{Def}_{lmi}.
*
* The multiplication is carried out for fixed DL:
*
*
* IF (.NOT. TRANSPOSEW) THEN
*
*    D(ij) = D(ij) + FACTOR*WBARDL(em,fj) * THETADL(em,fi)
*
* ELSE
*
*    D(ij) = D(ij) + FACTOR*WBARDL(fj,em) * THETADL(em,fi)
*
* END IF
*
* Filip Pawlowski, 11-Sep-2003, Aarhus.
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
      INTEGER ISWMATDL,ISTHETADL,LWORK
      INTEGER KWBARTR,KEND1,LWRK1,ISYMJ,ISYEMF,ISYMEM,ISYMFI,ISYMF
      INTEGER ISYMI,ISYMFJ,KOFF1,KFJ,KOFF2,KOFF3,NNEMF
C
#if defined (SYS_CRAY)
      REAL DIJ(*),WBARDL(*),THETADL(*)
      REAL WORK(LWORK)
      REAL ONE, FACTOR
#else
      DOUBLE PRECISION DIJ(*),WBARDL(*),THETADL(*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION ONE, FACTOR
#endif
C
      PARAMETER (ONE = 1.0D0)
C
      CALL QENTER('DIJCUB')
C

      IF (TRANSPOSEW) THEN
      !transpose Wbar^DL(em,fj) to Wbar^DL(fj,em))
         KWBARTR = 1
         KEND1   = KWBARTR + NT2SQ(ISWMATDL)
         LWRK1   = LWORK   - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*)'Memory available: ', LWORK
            WRITE(LUPRI,*)'Memory needed   : ', KEND1
            CALL QUIT('Insufficient memory in DIJ_CONT_CUB ')
         END IF
C
         CALL TRANS_AIBJ_BJAI(WBARDL,WORK(KWBARTR),ISWMATDL)
C
      END IF

C---------------------------------------------------------------
C     Calculate D(ij) = D(ij) - Wtilde_bar^DL(em,fj) t^DL(em,fi)
C---------------------------------------------------------------
C
C     -----------------------------------
C     Loop over outermost occupied index:
C     -----------------------------------
C
      DO ISYMJ = 1, NSYM
         ISYEMF = MULD2H(ISWMATDL,ISYMJ)
C
         DO J = 1, NRHF(ISYMJ)

C             -----------------------------------------------------
C             D(ij) <- D(ij)- sum_emf Wbar^DL(em,fj) THETA^DL(em,fi)
C             -----------------------------------------------------
            DO ISYMEM = 1, NSYM
               ISYMFI = MULD2H(ISTHETADL,ISYMEM)
               ISYMF  = MULD2H(ISYEMF,ISYMEM)
               ISYMI  = MULD2H(ISYMFI,ISYMF)
               ISYMFJ = MULD2H(ISYMF,ISYMJ)

C
               KOFF1 = 1     + IT2SQ(ISYMEM,ISYMFI)
     *                       + NT1AM(ISYMEM)*IT1AM(ISYMF,ISYMI)
C
               KFJ    = IT1AM(ISYMF,ISYMJ)+NVIR(ISYMF)*(J-1)+1
C
               IF (.NOT.TRANSPOSEW) THEN
                  KOFF2  = 1       + IT2SQ(ISYMEM,ISYMFJ)
     *                             + NT1AM(ISYMEM)*(KFJ-1)
               ELSE
                  KOFF2  = KWBARTR + IT2SQ(ISYMEM,ISYMFJ)
     *                             + NT1AM(ISYMEM)*(KFJ-1)
               END IF
C
               KOFF3  = IMATIJ(ISYMI,ISYMJ) 
     *                        + NRHF(ISYMI)*(J-1) + 1

               NNEMF  = MAX(NT1AM(ISYMEM)*NVIR(ISYMF),1)
C
               IF (.NOT.TRANSPOSEW) THEN
                  CALL DGEMV('T',NT1AM(ISYMEM)*NVIR(ISYMF),
     *                      NRHF(ISYMI),-FACTOR,THETADL(KOFF1),NNEMF,
     *                      WBARDL(KOFF2),1,ONE,DIJ(KOFF3),1)
               ELSE
                  CALL DGEMV('T',NT1AM(ISYMEM)*NVIR(ISYMF),
     *                      NRHF(ISYMI),-FACTOR,THETADL(KOFF1),NNEMF,
     *                      WORK(KOFF2),1,ONE,DIJ(KOFF3),1)
               END IF
C
            END DO ! ISYMFI

         END DO ! J 
      END DO    ! ISYMJ 
C
      CALL QEXIT('DIJCUB')
C
      RETURN
      END
C  /* Deck dab_cont_cub */
      SUBROUTINE DAB_CONT_CUB(TRANSPOSEW,DAB,WBARDL,ISWMATDL,THETADL,
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
      REAL ONE
#else
      DOUBLE PRECISION DAB(*),WBARDL(*),THETADL(*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION ONE
#endif
C
      PARAMETER (ONE = 1.0D0)
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
     *                       NT1AM(ISYMEM),ONE,WBARDL(KOFF1),NTOTEM,
     *                       THETADL(KOFF2),NTOTEM,ONE,DAB(KOFF3),
     *                       NTOTA)
               ELSE
                  CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),
     *                       NT1AM(ISYMEM),ONE,WORK(KOFF1),NTOTEM,
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
C  /* Deck read_t3_albj */
      SUBROUTINE READ_T3_ALBJ(LUFILE,FNFILE,ISYMT3,T2SQ,ISYMT2,
     *                        I,ISYMI,ISYMD)

      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <cc3t3d.h>
#include <ccorb.h>
C
      CHARACTER*(*) FNFILE
C
      INTEGER LUFILE,ISYMT3,ISYMT2,ISYMI,ISYMD
C
      INTEGER ISYMAIBJL,ISYMBJ,ISYMAIL,ISYMAL,ISYML,ISYMAI,ISYMAIBJ
      INTEGER ISYMJ,ISYMB,NBJ
      INTEGER KOFFT2,IADR
      INTEGER ISYMA
C
#if defined (SYS_CRAY)
      REAL T2SQ(*) 
#else
      DOUBLE PRECISION T2SQ(*)
#endif
C
      CALL QENTER('RDALBJ')
C
      ISYMAIBJL = MULD2H(ISYMT3,ISYMD)
      DO ISYMBJ = 1,NSYM
         ISYMAIL  = MULD2H(ISYMAIBJL,ISYMBJ)
         ISYMAL   = MULD2H(ISYMAIL,ISYMI)
         DO ISYML = 1,NSYM
            ISYMAI = MULD2H(ISYMAIL,ISYML)
            ISYMA  = MULD2H(ISYMAI,ISYMI)
            ISYMAIBJ = MULD2H(ISYMAI,ISYMBJ)
            DO ISYMJ = 1,NSYM
               ISYMB = MULD2H(ISYMBJ,ISYMJ)
               DO L = 1,NRHF(ISYML)
                  DO J = 1,NRHF(ISYMJ)
                     DO B = 1,NVIR(ISYMB)
            
C
                        NBJ  = IT1AM(ISYMB,ISYMJ)+NVIR(ISYMB)*(J-1)+B
C
                        KOFFT2 = IT2SQ(ISYMAL,ISYMBJ)
     *                         + NT1AM(ISYMAL)*(NBJ-1)
     *                         + IT1AM(ISYMA,ISYML)
     *                         + NVIR(ISYMA)*(L-1)
     *                         + 1
C
                        IADR = ISWTL(ISYMAIBJ,ISYML)
     *                       + NT2SQ(ISYMAIBJ)*(L-1)
     *                       + IT2SQ(ISYMAI,ISYMBJ)
     *                       + NT1AM(ISYMAI)*(NBJ-1)
     *                       + IT1AM(ISYMA,ISYMI)
     *                       + NVIR(ISYMA)*(I-1)
     *                       + 1

                        CALL GETWA2(LUFILE,FNFILE,T2SQ(KOFFT2),
     *                              IADR,NVIR(ISYMA))
C
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
C
C
      CALL QEXIT('RDALBJ')
C
      RETURN
      END
C  /* aden_dai_t2_d_cub */
      SUBROUTINE ADEN_DAI_T2_D_CUB(DAI,ISYMDAI,T2TP,ISYMT2,
     *                         TETA,ISYMTETA,ISYMD,D,
     *                         ISYML,L,WORK,LWORK)
*
************************************************************************
*
* Calculate contribution to the virtual part of Dai density for cubic
* response:
*
* DAI(ai) = DAI(ai) + T2^{de}_{lm} * (w^{Aed-}_{iml} - w^{Aed-}_{mil}).
*
* w^{Aed-}_{iml} is actually sitting as TETA^AI(em,dl) and, therfore,
* the two terms are calculated separetly.
*
*
*
* 1) DAI(ai) = DAI(ai) + T2^{de}_{lm} * w^{Aed-}_{iml} is calculated as:
*
*    DAI(AI) = DAI(AI) + KT2AM(em,dl) * TETA^AI(em,dl), which requires
*    sorting of T2 amplitudes first.
*
* 2) DAI(ai) = DAI(ai) - T2^{de}_{lm} * w^{Aed-}_{mil} is calculated as:
*
*    DAI(Ai) = DAI(Ai) - TETA^AM(dl,ei) * KT2AM(dl,eM)
*
************************************************************************
*     Written by F. Pawlowski, Fall 2003, Aarhus.
************************************************************************
*
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INTEGER ISYMDAI,ISYMT2,ISYMTETA,ISYMD,ISYML,LWORK
      INTEGER KT2AM,KEND1,LWRK1
      INTEGER ISYMA,ISYMI,ISYMDL,ISYMEM,KOFF1,KOFF2
      INTEGER ISYMM,ISYMEI,ISYME,KEM,KOFF3,NDLE
      INTEGER KT2AMTR,KTETATTR
      INTEGER ISYMAI
C
#if defined (SYS_CRAY)
      REAL DAI(*),T2TP(*),TETA(*),WORK(LWORK)
      REAL DDOT,ONE
#else
      DOUBLE PRECISION DAI(*),T2TP(*),TETA(*),WORK(LWORK)
      DOUBLE PRECISION DDOT,ONE
#endif
C
      PARAMETER (ONE = 1.0D0)
C
      CALL QENTER('DAIT2C')
C
C----------------------------------------------------------
C     1) DAI(ai) = DAI(ai) + T2^{de}_{lm} * w^{Aed-}_{iml}
C----------------------------------------------------------
C
      KT2AM = 1
      KEND1 = KT2AM + NT2SQ(ISYMT2)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in ADEN_DAI_T2_D_CUB (1)')
      END IF

      !first sort the T2 amplitudes:
      !T2^{de}_{lm} -> KT2AM(em,dl)
      CALL SORT_T2_AI_BJ(WORK(KT2AM),T2TP,ISYMT2)
C
* 1) DAI(ai) = DAI(ai) + T2^{de}_{lm} * w^{Aed-}_{iml} is calculated as:
*
*    DAI(AI) = DAI(AI) + KT2AM(em,dl) * TETA^AI(em,dl), which requires
*    sorting of T2 amplitudes first.

      !multiply KT2AM(em,dl) * TETA^AI(em,dl)
      ISYMA = ISYMD
      ISYMI = ISYML
      A     = D
      I     = L
C
      ISYMAI = MULD2H(ISYMA,ISYMI)
      IF (ISYMAI .EQ. ISYMDAI) THEN
C
         IF (ISYMT2 .EQ. ISYMTETA) THEN
C
            KOFF1  = IT1AM(ISYMA,ISYMI)
     *             + NVIR(ISYMA)*(I-1)
     *             + A
C
            DAI(KOFF1) = DAI(KOFF1) + DDOT(NT2SQ(ISYMTETA),TETA,1,
     *                                     WORK(KT2AM),1)
         END IF
      END IF
C
C
C----------------------------------------------------------
C     2) DAI(ai) = DAI(ai) - T2^{de}_{lm} * w^{Aed-}_{mil}
C----------------------------------------------------------
C
      ISYMA = ISYMD 
      ISYMM = ISYML
      A     = D
      M     = L
      DO ISYMDL = 1,NSYM
         ISYMEM = MULD2H(ISYMT2,ISYMDL)
         ISYMEI = MULD2H(ISYMTETA,ISYMDL)
         ISYME  = MULD2H(ISYMEM,ISYMM)
         ISYMI  = MULD2H(ISYMEI,ISYME)
         KOFF1 = 1 + IT2SQ(ISYMDL,ISYMEI)
     *             + NT1AM(ISYMDL)*IT1AM(ISYME,ISYMI)
C
         KEM   = IT1AM(ISYME,ISYMM) + NVIR(ISYME)*(M-1) + 1
C
         KOFF2 = KT2AM + IT2SQ(ISYMDL,ISYMEM)
     *                 + NT1AM(ISYMDL)*(KEM-1)
         KOFF3 = IT1AM(ISYMA,ISYMI)
     *         + A
C
         NDLE  = MAX(NT1AM(ISYMDL)*NVIR(ISYME),1)
C
         CALL DGEMV('T',NT1AM(ISYMDL)*NVIR(ISYME),
     *              NRHF(ISYMI),-ONE,TETA(KOFF1),NDLE,
     *              WORK(KOFF2),1,ONE,DAI(KOFF3),NVIR(ISYMA))
C
      END DO

      CALL QEXIT('DAIT2C')
C
      RETURN
      END
C  /* Deck sort_t2_ai_bj */
      SUBROUTINE SORT_T2_AI_BJ(T2AM,T2TP,ISYMT2)
C
C     Reorder t2 amplitudes as:
C
C     t2am(ai,bj) = t2tp(aijb)
C 
C     F. Pawlowski, Fall 2003, Aarhus.
C
#include <implicit.h>
C
      DIMENSION T2AM(*),T2TP(*)
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('T2AI_BJ')
C
      DO 100 ISYMB = 1,NSYM
C
         ISYAIJ = MULD2H(ISYMB,ISYMT2)
C
         DO 110 ISYMJ = 1,NSYM
C
            ISYMBJ = MULD2H(ISYMB,ISYMJ)
            ISYMAI = MULD2H(ISYMBJ,ISYMT2)
C
            DO 120 J = 1,NRHF(ISYMJ)
C
               DO 130 B = 1,NVIR(ISYMB)
C
                  NBJ   = IT1AM(ISYMB,ISYMJ)
     *                  + NVIR(ISYMB)*(J - 1) + B
C
                  KOFF1 = IT2SQ(ISYMAI,ISYMBJ)
     *                  + NT1AM(ISYMAI)*(NBJ - 1) + 1
C
                  KOFF2 = IT2SP(ISYAIJ,ISYMB)
     *                  + NCKI(ISYAIJ)*(B - 1)
     *                  + ISAIK(ISYMAI,ISYMJ)
     *                  + NT1AM(ISYMAI)*(J - 1) + 1
C
                  CALL DCOPY(NT1AM(ISYMAI),T2TP(KOFF2),1,T2AM(KOFF1),1)
C
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('T2AI_BJ')
C
      RETURN
      END
C  /* Deck wjk_ground_occ */
      SUBROUTINE WJK_GROUND_OCC(T30JK,ISYT30JK,T2TP,
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
C (4)
C                    - sum_l t^(ac)_(il) (bj|lk)
C
C     Filip Pawlowski, Aarhus, Winter 2003
*
*     Fixed for memory problems, 29-Oct-2003, Aarhus, FP.
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
      INTEGER ILOOP
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
      CALL QENTER('WJKGRO')
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
      KCBAI = 1
      KEND1 = KCBAI   + NMAAOBCI(ISYT30JK)
      LWRK1   = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in WJK_GROUND_OCC (1)')
      END IF
C
      DO ILOOP = 1,2
C
        KT2LBAI = KEND1
        KINTCL  = KT2LBAI + NT2SQ(ISYMT2)
        KEND2   = KINTCL  + NT1AM(ISYCL)
        LWRK2   = LWORK - KEND2
C
        IF (LWRK2 .LT. 0) THEN
           WRITE(LUPRI,*) 'Memory available : ',LWORK
           WRITE(LUPRI,*) 'Memory needed    : ',KEND2
           IF (ILOOP .EQ. 1) THEN
              CALL QUIT('Insufficient space in WJK_GROUND_OCC (1x)')
           ELSE
              CALL QUIT('Insufficient space in WJK_GROUND_OCC (1xx)')
           END IF
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
        IF (ILOOP .EQ. 1) THEN
          CALL SORT_INT_AJ_IK(WORK(KINTCL),T3OG2,ISYINT,ISYMK,K,ISYMJ,J)
        ELSE
          CALL SORT_INT_AJ_IK(WORK(KINTCL),T3OG2,ISYINT,ISYMJ,J,ISYMK,K)
        ENDIF
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
     *              + IT1AM(ISYMC,ISYML)
              KOFF2 = KT2LBAI
     *              + IMAJBAI(ISYML,ISYBAI)
              KOFF3 = KCBAI
     *              + IMAAOBCI(ISYMC,ISYBAI)
C
              NTOTC = MAX(NVIR(ISYMC),1)
              NTOTL = MAX(NRHF(ISYML),1)
C
              CALL DGEMM('N','N',NVIR(ISYMC),NMAABI(ISYBAI),NRHF(ISYML),
     *                   -ONE,WORK(KOFF1),NTOTC,WORK(KOFF2),NTOTL,
     *                   ONE,WORK(KOFF3),NTOTC)
C
        END DO ! ISYML
C
C       T30JK(bcai) = T30JK(bcai) + T^JK(cbai)
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
              CALL QUIT('Insufficient space in WJK_GROUND_OCC (1a)')
           END IF
C
           CALL DZERO(WORK(KTEMP),NMAABCI(ISYT30JK))
C
           ! Sort from KCBAI(c,bai) to KTEMP(c,b,a,i)
           CALL FA_BCI(WORK(KTEMP),WORK(KCBAI),ISYT30JK,1)
           CALL DCOPY(NMAAOBCI(ISYT30JK),WORK(KTEMP),1,WORK(KCBAI),1)
        END IF

        IF (ILOOP .EQ. 1) THEN
          !sort W(cbai) as W(bcai)
          CALL FBACI(T30JK,WORK(KCBAI),ISYT30JK)
        ELSE
          !put W(bcai) in the final array
          CALL FABCI_COLLECT(T30JK,WORK(KCBAI),ISYT30JK)
        ENDIF

      ENDDO
C
      CALL QEXIT('WJKGRO')
C
      RETURN
      END 
C  /* Deck wjk_t2 */
      SUBROUTINE WJK_T2(FAC,J,ISYMJ,K,ISYMK,T2TPX,ISYMT2X,T2TPZ,
     *                 ISYMT2Z,
     *                 FOCKY,ISYFKY,
     *                 WMAT,ISWMAT,WRK,LWRK)
C 
C WJK(bda,i) = WBD(bda,i) -
C      sum (f,l) focky(l,f)*( t2X(ai,dl)*t2Z(fk,bj) + t2X(ai,bl)*t2Z(fj,dk) )
C
C
C
C     Written by F. Pawlowski, Fall 2003.
C
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
      INTEGER LWRK, KFCLF, KEND0, LWRK0, KOFF1, KOFF2, KTB, KEND1, LWRK1
      INTEGER NL, NF, KOFFY, KOFFT2, KOFFT, KOFFW, KTD, KW
      INTEGER ISYMB, ISYMD, ISYMT2X, ISYFKY, ISWMAT 
      INTEGER ISYAIL, ISYAI, ISYAIK, NA, NAI, LENGTH
      INTEGER ISYMF, ISYML, ISYFKJ, ISYTB, ISYMJ, ISYFK, ISYMK, ISYLK
      INTEGER ISYFJK, ISYTD, ISYLJ, ISYFJ, ISYAIJ 
      INTEGER ISYMT2Z
C
      INTEGER ISYMKJ,ISYFB,ISYTJK,KTJK,ISYLB,ISYMBDAI,KWTEMP,ISYDIA
      INTEGER ISYMI,ISYMLI,ISYMBDA,ISYMA,ISYMLIA,ISYMBD,ISYMDL
      INTEGER ISYMDLI,ND,NB
      INTEGER KT2FB
      INTEGER ILOOP
C
#if defined (SYS_CRAY)
      REAL T2TPX(*), FOCKY(*), WMAT(*), WRK(*)
      REAL HALF, ONE, ZERO
      REAL T2TPZ(*)
      REAL FAC 
#else
      DOUBLE PRECISION T2TPX(*), FOCKY(*), WMAT(*), WRK(*)
      DOUBLE PRECISION HALF, ONE, ZERO
      DOUBLE PRECISION T2TPZ(*)
      DOUBLE PRECISION FAC 
#endif
C
      PARAMETER (HALF = 0.5D0, ONE = 1.0D0, ZERO = 0.0D0)
C
      CALL QENTER('WJKT2')
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
         WRITE(LUPRI,*) 'Memory available : ',LWRK0
         WRITE(LUPRI,*) 'Memory needed    : ',KEND0
         CALL QUIT('Insufficient space in WJK_T2 (1)')
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
C    calculate first t2 contribution to W matrix
C
C construct tZJK(l,b) = sum (f) focky(l,f)*t2tpZ(f,K,J,b)
C
C calculated as:
C            tZJK(l,b) = sum (f) focky(l,f)*t2tpZJK(f,b)
C
      DO ILOOP = 1,2
C
        ISYMKJ = MULD2H(ISYMK,ISYMJ)
        ISYFB  = MULD2H(ISYMT2Z,ISYMKJ) 
        ISYTJK    = MULD2H(ISYFKY,ISYFB) 
C
        KTJK      = KEND0
        KEND1    = KTJK  + NT1AM(ISYTJK)
        LWRK1    = LWRK  - KEND1
C
        KT2FB    = KEND1
        KEND1    = KT2FB + NMATAB(ISYFB)
        LWRK1    = LWRK  - KEND1
C
        IF (LWRK1 .LT. 0) THEN
           WRITE(LUPRI,*) 'Memory available : ',LWRK
           WRITE(LUPRI,*) 'Memory needed    : ',KEND1
           CALL QUIT('Insufficient space in WJK_T2 (2)')
        END IF
C
        CALL DZERO(WRK(KTJK),NT1AM(ISYTJK))
C
        IF (ILOOP .EQ. 1) THEn
        !sort t2tpZ(f,K,J,b) to KT2FB(f,b)
        CALL SORT_T2_AB(WRK(KT2FB),ISYMK,K,ISYMJ,J,T2TPZ,ISYMT2Z)
        ELSE
        !sort t2tpZ(f,J,K,b) to KT2FB(f,b)
        CALL SORT_T2_AB(WRK(KT2FB),ISYMJ,J,ISYMK,K,T2TPZ,ISYMT2Z)
        ENDIF
C
        !tZJK(l,b) = sum (f) focky(l,f)*t2tpZJK(f,b)
        DO ISYMF = 1,NSYM
           ISYML = MULD2H(ISYFKY,ISYMF)
           ISYMB = MULD2H(ISYFB,ISYMF)
C
           KOFFY  = KFCLF + IT1AMT(ISYML,ISYMF)
           KOFFT2 = KT2FB +IMATAB(ISYMF,ISYMB)
           KOFFT  = KTJK  + IT1AMT(ISYML,ISYMB)
C
           NL = MAX(NRHF(ISYML),1)
           NF = MAX(NVIR(ISYMF),1)
C
           CALL DGEMM('N','N',NRHF(ISYML),NVIR(ISYMB),
     *                NVIR(ISYMF),ONE,WRK(KOFFY),NL,
     *                WRK(KOFFT2),NF,ONE,WRK(KOFFT),NL)

C
        END DO !ISYMF

C 
C WJK(bda,i) = WBD(bda,i) -
C      sum (f,l) focky(l,f)*t2X(ai,dl)*t2Z(fk,bj) 
C            = WBD(bda,i) -
C      sum (l) t2tpX(a,i,l,d) * tZJK(l,b)
C

C Multiply as tZJK(l,b) * t2tpX(d,l,i,a) --> WJK(bda,i)

C

        ISYMBDAI = MULD2H(ISYMT2X,ISYTJK)
C
        !symmmetry check
        IF (ISYMBDAI .NE. ISWMAT) THEN
           WRITE(LUPRI,*) 'ISYMBDAI = ', ISYMBDAI
           WRITE(LUPRI,*) 'ISWMAT = ', ISWMAT
           WRITE(LUPRI,*) 'These symmetries should be EQUAL!'
           CALL QUIT('Symmetry inconsistency in WJK_T2')
        END IF
C
        KWTEMP   = KEND1
        KEND1    = KWTEMP  + NMAABCI(ISYMBDAI)
        LWRK1    = LWRK  - KEND1
C
        IF (LWRK1 .LT. 0) THEN
           WRITE(LUPRI,*) 'Memory available : ',LWRK
           WRITE(LUPRI,*) 'Memory needed    : ',KEND1
           CALL QUIT('Insufficient space in WJK_T2 (3)')
        END IF
C
        CALL DZERO(WRK(KWTEMP),NMAABCI(ISYMBDAI))
C
        !Multiply as tZJK(l,b) * t2tpX(d,l,i,a) --> WJK(bda,i)
        DO ISYML = 1,NSYM
           ISYDIA = MULD2H(ISYMT2X,ISYML)
           ISYMB = MULD2H(ISYTJK,ISYML)
           DO ISYMI = 1,NSYM
              ISYMLI = MULD2H(ISYML,ISYMI)
              ISYMBDA = MULD2H(ISYMBDAI,ISYMI)
              DO ISYMA = 1,NSYM
                 ISYMLIA = MULD2H(ISYMLI,ISYMA)
                 ISYMBD  = MULD2H(ISYMBDA,ISYMA)
                 ISYMD = MULD2H(ISYMT2X,ISYMLIA)
                 ISYMDL = MULD2H(ISYMD,ISYML)
                 ISYMDLI = MULD2H(ISYMD,ISYMLI)
                 DO I = 1,NRHF(ISYMI)
                    DO A = 1,NVIR(ISYMA)
C
                     KOFFT =  KTJK + IT1AMT(ISYML,ISYMB)
                     KOFFT2 = IT2SP(ISYMDLI,ISYMA) + NCKI(ISYMDLI)*(A-1)
     *                       + ISAIK(ISYMDL,ISYMI) + NT1AM(ISYMDL)*(I-1)
     *                       + IT1AM(ISYMD,ISYML)  + 1
                     KOFFW  = KWTEMP + IMAABCI(ISYMBDA,ISYMI) 
     *                     + NMAABC(ISYMBDA)*(I-1)
     *                     + IMAABC(ISYMBD,ISYMA) + NMATAB(ISYMBD)*(A-1)
     *                     + IMATAB(ISYMB,ISYMD)  
C
                     NL = MAX(NRHF(ISYML),1)
                     ND = MAX(NVIR(ISYMD),1)
                     NB = MAX(NVIR(ISYMB),1)
C
                     CALL DGEMM('T','T',NVIR(ISYMB),NVIR(ISYMD),
     *                          NRHF(ISYML),-FAC,WRK(KOFFT),NL,
     *                          T2TPX(KOFFT2),ND,ONE,WRK(KOFFW),NB)

C
                    END DO
                 END DO
              END DO
           END DO
        END DO
C
        IF (ILOOP .EQ. 1) THEn
           !First contribution
           CALL FABCI_COLLECT(WMAT,WRK(KWTEMP),ISYMBDAI)
        ELSE
           !Second contribution( (bj) <-> (dk) permutation)
           CALL FBACI(WMAT,WRK(KWTEMP),ISYMBDAI)
        ENDIF

      ENDDO
C
      CALL QEXIT('WJKT2')
C
      RETURN
      END
C  /* Deck fabci */
      SUBROUTINE FABCI_COLLECT(TABCI,TABCITMP,ISYMT)
C
C     TABCI = TABCI + TABCITMP
C
C     F. Pawlowski, Aarhus, Fall 2003
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
      REAL TABCI(*),TABCITMP(*)
#else
      DOUBLE PRECISION TABCI(*),TABCITMP(*)
#endif
C
      CALL QENTER('FABCICLL')
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
C
                           TABCI(KOFF1) = TABCI(KOFF1) + TABCITMP(KOFF1)
C
                        END DO   
                     END DO   
                  END DO   
               END DO   
            END DO   
         END DO   
      END DO   
C
      CALL QEXIT('FABCICLL')
C
      RETURN
      END
C  /* Deck wxbd_t2_cub */
      SUBROUTINE WXBD_T2_CUB(T2XNET2Z,AIBJCK_PERM,B,ISYMB,D,ISYMD,
     *                     T2TPX,ISYMT2X,
     *                     T2TPZ,ISYMT2Z,FOCKY,
     *                     ISYFKY,INDSQ,LENSQ,WMAT,ISWMAT,WRK,LWRK)

      IMPLICIT NONE

#include <priunit.h>
#include <ccsdsym.h>

      LOGICAL T2XNET2Z

      INTEGER AIBJCK_PERM

      INTEGER LENSQ, INDSQ(LENSQ,6), LWRK
      INTEGER ISYMB, ISYMD, ISYMT2X, ISYFKY, ISWMAT, ISYMT2Z

#if defined (SYS_CRAY)
      REAL T2TPX(*), FOCKY(*), WMAT(*), WRK(*), T2TPZ(*)
#else
      DOUBLE PRECISION T2TPX(*), FOCKY(*), WMAT(*), WRK(*), T2TPZ(*)
#endif

      CALL QENTER('WXBDT2C')

      CALL WXBD_T2_1(AIBJCK_PERM,B,ISYMB,D,ISYMD,T2TPX,ISYMT2X,
     *                     T2TPZ,ISYMT2Z,FOCKY,
     *                     ISYFKY,INDSQ,LENSQ,WMAT,ISWMAT,WRK,LWRK)

      IF (T2XNET2Z) THEN
C
         CALL WXBD_T2_1(AIBJCK_PERM,B,ISYMB,D,ISYMD,T2TPZ,ISYMT2Z,
     *                     T2TPX,ISYMT2X,FOCKY,
     *                     ISYFKY,INDSQ,LENSQ,WMAT,ISWMAT,WRK,LWRK)

      END IF

      CALL QEXIT('WXBDT2C')

      RETURN
      END


C  /* Deck wxbd_t2_1 */
      SUBROUTINE WXBD_T2_1(AIBJCK_PERM,B,ISYMB,D,ISYMD,T2TPX,ISYMT2X,
     *                     T2TPZ,ISYMT2Z,FOCKY,
     *                     ISYFKY,INDSQ,LENSQ,WMAT,ISWMAT,WRK,LWRK)


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
C     Written by F. Pawlowski, Spring 2003.
C

      IMPLICIT NONE
C
      INTEGER AIBJCK_PERM
C
      INTEGER LENSQ
      INTEGER INDSQ(LENSQ,6)
      INTEGER LWRK,KFCLF, KEND0, LWRK0, KTB, KEND1, LWRK1 
      INTEGER NL, NF, KOFFY, KOFFT2, KOFFT, KOFFW, KTD, KW
      INTEGER ISYMB, ISYMD, ISYMT2X, ISYMT2Z, ISYFKY, ISWMAT 
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
      REAL T2TPX(*), T2TPZ(*), FOCKY(*), WMAT(*), WRK(*)
      REAL HALF, ONE, ZERO
#else
      DOUBLE PRECISION T2TPX(*), T2TPZ(*), FOCKY(*), WMAT(*), WRK(*)
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
      CALL QENTER('WXBDT21')
C
C---------------------------------------
C     Initial test of AIBJCK_PERM option
C---------------------------------------
C
      IF ( (AIBJCK_PERM .LT. 1) .OR. (AIBJCK_PERM .GT. 4) ) THEN
         WRITE(LUPRI,*)'AIBJCK_PERM = ',AIBJCK_PERM
         WRITE(LUPRI,*)'should be between 1 and 4 '
         CALL QUIT('Illegal value of AIBJCK_PERM option in WXBD_T2_1')
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
         CALL QUIT('Insufficient space in WXBD_T2_1 (1)')
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
C construct tZB(l,k,j) = sum (f) focky(l,f)*t2tpZ(f,k,j,B)
C
         ISYFKJ   = MULD2H(ISYMT2Z,ISYMB) 
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
            CALL QUIT('Insufficient space in WXBD_T2_1 (2)')
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
     *                    T2TPZ(KOFFT2),NF,ONE,WRK(KOFFT),NL)
C
               END DO
            END DO
         END DO
C
C         WBD(a,i,k,j) = WBD(a,i,k,j) -
C                        sum (f,l) focky(l,f)* t2X(ai,Dl)*t2Z(fk,Bj) 
C                      = WBD(a,i,k,j) - 
C                        sum(l) t2tpX(a,i,l,D) * tZB(l,k,j)
C
         ISYAIL = MULD2H(ISYMT2X,ISYMD)
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
     *                       NRHF(ISYML),-ONE,T2TPX(KOFFT2),NAI,
     *                       WRK(KOFFT),NL,ONE,WMAT(KOFFW),NAI)

C
               END DO
            END DO
         END DO
C
C    calculate second t2 contribution to W matrix
C
C
C construct tD(l,j,k) = sum (f) focky(l,f)*t2tpZ(f,j,k,D)
C
         ISYFJK   = MULD2H(ISYMT2Z,ISYMD) 
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
            CALL QUIT('Insufficient space in WXBD_T2_1 (3)')
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
     *                       T2TPZ(KOFFT2),NF,ONE,WRK(KOFFT),NL)
C
               END DO
            END DO
         END DO
C
C      WBD(a,i,k,j) = WBD(a,i,k,j) -
C                        sum (f,l) focky(l,f)*t2X(ai,Bl)*t2Z(fj,Dk) )
C                   = WBD(a,i,k,j) -
C                        sum(l) t2tpX(a,i,l,B) * tZD(l,j,k)
C
         ISYAIL = MULD2H(ISYMT2X,ISYMB)
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
     *                       NRHF(ISYML),-ONE,T2TPX(KOFFT2),NAI,
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
C     WBD(aikj) = WBD(aikj)  -  focky(l,f)* t2X(bj,al)*t2Z(fi,dk)
C
C                                          TX^B(ajl)    tZ(fikD)
C                  
C                work(lik) = focky(lf) * tZ(fikD) 
C
C                work(ajik) = TX^B(ajl) * work(lik)

C
C                work(lik) = focky(lf) * tZ(fikD) 
C
         ISYFIK     = MULD2H(ISYMT2Z,ISYMD)
         ISYLIK     = MULD2H(ISYFKY,ISYFIK)
         ISYAJL     = MULD2H(ISYMT2X,ISYMB)
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
            CALL QUIT('Insufficient space in WXBD_T2_1 (4)')
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
     *                       T2TPZ(KOFF2),NF,ONE,WRK(KOFF3),NL)
C
               END DO
            END DO
         END DO
C
C        TX^B(ajl) =     t2X(bj,al)
C
         CALL SORT_T2_AJI(WRK(KAJL),ISYMB,B,T2TPX,ISYMT2X)
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
C                work(ajik) = TX^B(ajl) * work(lik)
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
C     WBD(aikj) = WBD(aikj) - focky(l,f)* t2X(bj,dl)*t2Z(fk,ai)
C                                                     
C                                           TX^DB(lj)  tZ(aikf)
C             
C                work(fj) = focky(l,f) * TX^DB(lj)
C                
C              WMAT(aikj) = WMAT(aikj) - tZ(aikf) * work(fj) 
C
         ISYBD    = MULD2H(ISYMB,ISYMD)
         ISYLJ    = MULD2H(ISYBD,ISYMT2X)
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
            CALL QUIT('Insufficient space in WXBD_T2_1 (5)')
         END IF
C
C
         CALL  SORT_T2_IJ(WRK(KLJ),ISYMD,D,ISYMB,B,T2TPX,ISYMT2X)
C
C                work(fj) = focky(l,f) * TX^DB(lj)
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
C              WMAT(aikj) = WMAT(aikj) - tZ(aikf) * work(fj) 
C
         DO ISYMJ = 1,NSYM
            ISYMF = MULD2H(ISYMJ,ISYFJ)
            ISYAIK = MULD2H(ISYMT2Z,ISYMF)
C
            KOFF1 =  IT2SP(ISYAIK,ISYMF) + 1
            KOFF2 =  KFJ +  IT1AM(ISYMF,ISYMJ)
            KOFF3 =  ISAIKJ(ISYAIK,ISYMJ) + 1
C
            NAIK = MAX(1,NCKI(ISYAIK))
            NF = MAX(1,NVIR(ISYMF))

            CALL DGEMM('N','N',NCKI(ISYAIK),NRHF(ISYMJ),
     *                 NVIR(ISYMF),-ONE,T2TPZ(KOFF1),NAIK,
     *                 WRK(KOFF2),NF,ONE,WMAT(KOFF3),NAIK)
C
         END DO
C
      END IF
      IF ((AIBJCK_PERM.EQ.3) .OR. (AIBJCK_PERM.EQ.4)) THEN
C
C     WBD(aikj) = WBD(aikj)  -  focky(l,f)* t2X(dk,al)*t2Z(fi,bj)
C 
C                                            IX^D(alk)  TZ^B(fij)
C
C                work(lij) = focky(lf) * TZ^B(fij)
C                 
C                work(akij) = TX^D(akl) * work(lij)
C
C     WBD(aikj) = WBD(aikj)  -  work(akij)
C
C
C                work(lij) = focky(lf) * TZ^B(fij)
C
         ISYFIJ     = MULD2H(ISYMT2Z,ISYMB)
         ISYLIJ     = MULD2H(ISYFKY,ISYFIJ)
         ISYAKL     = MULD2H(ISYMT2X,ISYMD)
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
            CALL QUIT('Insufficient space in WXBD_T2_1 (6)')
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
     *                       T2TPZ(KOFF2),NF,ONE,WRK(KOFF3),NL)
C
               END DO
            END DO
         END DO
C
C                work(akij) = TX^D(akl) * work(lij)
C
         CALL SORT_T2_AJI(WRK(KAKL),ISYMD,D,T2TPX,ISYMT2X)
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
C     WBD(aikj) = WBD(aikj)  -  focky(l,f)* t2X(dk,bl)*t2Z(fj,ai)
C                                            
C                                           IX^BD(lk)  IZ(aijf)
C
C                  work(fk) = focky(lf) * IX^BD(lk)
C                  
C                  work(aijk) = IZ(aijf) * work(fk)
C
C     WBD(aikj) = WBD(aikj)  - work(aijk)
C
         ISYBD    = MULD2H(ISYMB,ISYMD)
         ISYLK    = MULD2H(ISYBD,ISYMT2X)
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
            CALL QUIT('Insufficient space in WXBD_T2_1 (7)')
         END IF
C
         CALL  SORT_T2_IJ(WRK(KLK),ISYMB,B,ISYMD,D,T2TPX,ISYMT2X)
C
C                  work(fk) = focky(lf) * IX^BD(lk)

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
C                  work(aijk) = IZ(aijf) * work(fk)
C
         DO ISYMK = 1,NSYM
            ISYMF = MULD2H(ISYMK,ISYFK)
            ISYAIJ = MULD2H(ISYMT2Z,ISYMF)
C
            KOFF1 =  IT2SP(ISYAIJ,ISYMF) + 1 
            KOFF2 =  KFK +  IT1AM(ISYMF,ISYMK)
            KOFF3 =  KW  +  ISAIKJ(ISYAIJ,ISYMK) 
C
            NAIJ = MAX(1,NCKI(ISYAIJ))
            NF = MAX(1,NVIR(ISYMF))

            CALL DGEMM('N','N',NCKI(ISYAIJ),NRHF(ISYMK),
     *                 NVIR(ISYMF),-ONE,T2TPZ(KOFF1),NAIJ,
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

      CALL QEXIT('WXBDT21')
C
      RETURN
      END
C  /* Deck intvir_t3x_jk */
      SUBROUTINE INTVIR_T3X_JK(XGBCDK,ISYINT,
     *                         LUFIL,FNFIL,
     *                         WORK,LWORK) 
**********************************************************
*
*     Construct the integrals used for t3x^JK calculation
*
*     Read virtual integrals (Ck|bd) stored as I^C(dk,b)
*     Final sort (Ck|bd)  as I(bcd,k) 
*
*     OUTPUT (XGBCDK) : g(ckbd) = (ck|bd) sorted as I(bcd,k) 
*
*     F. Pawlowski, 02-10-2003, Aarhus.
**********************************************************
C
      IMPLICIT NONE
#include <ccsdsym.h>
#include <ccorb.h>
#include <priunit.h>
C
      INTEGER ISYINT, LUFIL, LWORK
      INTEGER ISYMD, ISYCKA, KINTVI, KEND1, LWRK1, IOFF
C
      CHARACTER*(*) FNFIL
C
#if defined (SYS_CRAY)
      REAL XGBCDK(*), WORK(LWORK)
#else
      DOUBLE PRECISION XGBCDK(*), WORK(LWORK) 
#endif
C
      CALL QENTER('INTV3XJK')
C
C***********************************************************'
C     Get  (XGBDCK) : g(ckbd) = (ck|bd) sorted as I(bcd,k) 
C***********************************************************'
C
      DO ISYMD = 1, NSYM
         ISYCKA = MULD2H(ISYINT,ISYMD)
C
         KINTVI = 1
         KEND1 = KINTVI + NCKATR(ISYCKA)
         LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in INTVIR_T3X_JK ')
      ENDIF
C
         DO D = 1, NVIR(ISYMD)
C
C     Read virtual integrals (ck|bd) sorted as I^C(dk,b)
C
            IOFF = ICKBD(ISYCKA,ISYMD) + NCKATR(ISYCKA)*(D - 1) + 1
            IF (NCKATR(ISYCKA) .GT. 0) THEN
               CALL GETWA2(LUFIL,FNFIL,WORK(KINTVI),IOFF,
     &                        NCKATR(ISYCKA))
            ENDIF
C
C     Final sort (ck|bd) as I(bcd,k)
C
            CALL SORT_INTVIR_T3B0(XGBCDK,WORK(KINTVI),
     *                           D,ISYMD,ISYCKA,WORK(KEND1),LWRK1)
         END DO !  D
      END DO !  ISYMD
C
      CALL QEXIT('INTV3XJK')
C
      RETURN
      END
