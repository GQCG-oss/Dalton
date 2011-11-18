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
      SUBROUTINE CC3_GMAT(LISTL,IDLSTL,LISTB,IDLSTB,
     &                    LISTC,IDLSTC,OMEGA1,OMEGA2,ISYRES,WORK,LWORK,
     *                    LUCKJD,FNCKJD,LUTOC,FNTOC,
     *                    LU3VI,FN3VI,LUDKBC3,FNDKBC3,
     *                    LU3FOPX,FN3FOPX,LU3FOP2X,FN3FOP2X)
C
*    Purpose: compute triples contribution to G matrix transformation
*
*  (G T^B T^C)^eff_1,2 = (G T^B T^C)_1,2(CCSD) 
*                            + (G T^B T^C)_1,2(L_3)
*        
C    (G T^B T^C)_1,2(L_3) =
C 
C                      +  <L3|[H^BC,T^{0}_{2}],tau_{1}]|HF> 
C
C                      +  <L3|[H^B,T^{C}_{2}],tau_{1}]|HF> 
C
C                      +  <L3|[H^C,T^{B}_{2}],tau_{1}]|HF> 
C
C                      +  <L3|[H^BC,tau_{2}]|HF>
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <dummy.h>
#include <ccsdsym.h>
#include <inftap.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <iratdef.h>
#include <ccinftap.h>
C
      INTEGER IDLSTL,IDLSTB,IDLSTC,LWORK,LUCKJD,LUTOC,LU3VI,LUDKBC3
      INTEGER LU3FOPX,LU3FOP2X
      INTEGER LU4VBC,LU4VC,LU4VB,LU3SRTR,LUCKJDR,LUDELDR,LUDKBCR
      INTEGER LUDKBCR4
      INTEGER ISINT1,ISINT2,ISYMT1,ISYMT2,ISYML,ISYML1
      INTEGER ISYML2,ISYMB,ISYMB1,ISYMB2,ISYMC1,ISYMC2
      INTEGER ISINTB,ISINTC,ISYMBC,ISINTBC,ISYMIM,ISYRES
      INTEGER KFOCKD,KFCKBA,KT1AM,KL1AM,KL2TP,KB1AM,KB2TP
      INTEGER KC1AM,KC2TP,KT2TP,KLAMDP,KLAMDH,KEND0,LWRK0
      INTEGER IOPT
      INTEGER KBIOOOO,KBIOVVO,KBIOOVV,KCIOOOO,KCIOVVO,KCIOOVV
      INTEGER KBCIOOOO,KBCIOVVO,KBCIOOVV
      INTEGER ISYMC,ISYMK,KOFF1,KOFF2
      INTEGER KTROC,KTROC0,KTROC1,KTROC2,KTROC3,KTROC4,KXIAJB
      INTEGER KEND1,LWRK1,KINTOC,MAXOCC,KEND2,LWRK2,LENGTH
      INTEGER ISYOPE,IOPTTCME,IOFF
      INTEGER ISYMD,ISAIJ1,ISYCKB,ISCKB1,ISCKB2,ISCKBB,ISCKBC,ISCKBBC
      INTEGER KTRVI7,KTRVI8,KTRVI2,KRMAT1,KTRVI0,KTRVI3,KTRVI4,KTRVI5
      INTEGER KTRVI6,KVVVVB,KVVVVC,KVVVVBC,KEND3,LWRK3,KINTVI
      INTEGER KEND4,LWRK4
      INTEGER ISYALJ,ISAIJ2,ISYMBD,ISCKIJ,KSMAT,KQMAT,KDIAG
      INTEGER KINDSQ,KINDEX,KTMAT,KRMAT2,LENSQ,ISYRES1,KOME1,KOME2
      INTEGER LUFCK
C
COMMENT COMMENT
COMMENT COMMENT
      integer kt3am
COMMENT COMMENT
COMMENT COMMENT
C
C     Functions : 
      INTEGER ILSTSYM
C
#if defined (SYS_CRAY)
      REAL OMEGA1(*), OMEGA2(*)
      REAL WORK(LWORK)
      REAL XNORMVAL,DDOT
      REAL HALF , ONE , TWO 
#else
      DOUBLE PRECISION OMEGA1(*), OMEGA2(*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION XNORMVAL,DDOT
      DOUBLE PRECISION HALF , ONE , TWO 
#endif
C
C
      CHARACTER*(*) FNCKJD, FNTOC, FN3VI, FNDKBC3, FN3FOPX
      CHARACTER*(*) FN3FOP2X
      CHARACTER*1 CDUMMY
      CHARACTER*3 LISTL, LISTB, LISTC
      CHARACTER*10 MODEL
      CHARACTER*13 FN4VC, FN4VB, FN4VBC, FN3SRTR, FNCKJDR
      CHARACTER*13 FNDELDR, FNDKBCR, FNDKBCR4
C
      PARAMETER(HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      CALL QENTER('CC3_GMAT')
C
C----------------------------------------------------
C     Initialise character strings and open files
C----------------------------------------------------
C
      CDUMMY = ' '
C
      LU4VBC   = -1
      LU4VC    = -1
      LU4VB    = -1
      LU3SRTR  = -1
      LUCKJDR  = -1
      LUDELDR  = -1
      LUDKBCR  = -1
      LUDKBCR4 = -1
C
      FN4VBC   = 'CC3_GMAT_TMP1'
      FN4VC    = 'CC3_GMAT_TMP2'
      FN4VB    = 'CC3_GMAT_TMP3'
      FN3SRTR  = 'CC3_GMAT_TMP4'
      FNCKJDR  = 'CC3_GMAT_TMP5'
      FNDELDR  = 'CC3_GMAT_TMP6'
      FNDKBCR  = 'CC3_GMAT_TMP7'
      FNDKBCR4 = 'CC3_GMAT_TMP8'
C
      CALL WOPEN2(LU4VBC,FN4VBC,64,0)
      CALL WOPEN2(LU4VC,FN4VC,64,0)
      CALL WOPEN2(LU4VB,FN4VB,64,0)
      CALL WOPEN2(LU3SRTR,FN3SRTR,64,0)
      CALL WOPEN2(LUCKJDR,FNCKJDR,64,0)
      CALL WOPEN2(LUDELDR,FNDELDR,64,0)
      CALL WOPEN2(LUDKBCR,FNDKBCR,64,0)
      CALL WOPEN2(LUDKBCR4,FNDKBCR4,64,0)
C
C-------------------------------------------------------------
C     Set symmetry flags.
C
C     omega = int1*T2* ISYMB * ISYMC
C     isymres is symmetry of result(omega)
C     isint1 is symmetry of integrals in contraction.(int1)
C     the int1 integrals are one index transformed with  B or C or BC
C     isint2 is symmetry of integrals in the triples equation.(int2)
C     isymim is symmetry of S and Q intermediates.(t2*int2)
C      (sym is for all index of S and Q (cbd,klj)
C       thus cklj=b*d*isymim)
C-------------------------------------------------------------
C
      IPRCC   = IPRINT
      ISINT1  = ISYMOP
      ISINT2  = ISYMOP
      ISYMT1  = 1
      ISYMT2  = 1
      ISYML   = ILSTSYM(LISTL,IDLSTL)
      ISYML1  = ISYML
      ISYML2  = ISYML
      ISYMB   = ILSTSYM(LISTB,IDLSTB)
      ISYMB1  = ISYMB
      ISYMB2  = ISYMB
      ISYMC   = ILSTSYM(LISTC,IDLSTC)
      ISYMC1  = ISYMC
      ISYMC2  = ISYMC

      ISINTB   = MULD2H(ISINT1,ISYMB)
      ISINTC   = MULD2H(ISINT1,ISYMC)
      ISYMBC   = MULD2H(ISYMB,ISYMC)
      ISINTBC  = MULD2H(ISINT1,ISYMBC)
      ISYMIM   = MULD2H(ISYML,ISYMOP)
      ISYRES1  = MULD2H(ISYMIM,ISYMBC)
      IF (ISYRES .NE. ISYRES1) THEN
         WRITE(LUPRI,*) ' SYMMETRY MISMATCH IN CC3_GMAT'
         WRITE(LUPRI,*) ' ISYRES =', ISYRES
         WRITE(LUPRI,*) ' ISYRES1 =', ISYRES1
         CALL QUIT('Symmetry mismatch IN CC3_GMAT')
      END IF
C
C-----------------------------------
C     Work space allocation
C-----------------------------------
C
      KOME1   = 1
      KOME2   = KOME1  + NT1AM(ISYRES)
      KFOCKD  = KOME2  + NT2SQ(ISYRES) 
      KFCKBA  = KFOCKD + NORBTS
      KT1AM   = KFCKBA + N2BST(ISYMOP)
      KL1AM   = KT1AM  + NT1AM(ISYMT1)
      KL2TP   = KL1AM  + NT1AM(ISYML1)
      KB1AM   = KL2TP  + NT2SQ(ISYML2)
      KB2TP   = KB1AM  + NT1AM(ISYMB1)
      KC1AM   = KB2TP  + NT2SQ(ISYMB2) 
      KC2TP   = KC1AM  + NT1AM(ISYMC1)
      KT2TP   = KC2TP  + NT2SQ(ISYMC2)
      KLAMDP  = KT2TP  + NT2SQ(ISYMT2)
      KLAMDH  = KLAMDP + NLAMDT
      KEND0   = KLAMDH + NLAMDT
      LWRK0   = LWORK  - KEND0
C
      CALL DZERO(WORK(KOME1),NT1AM(ISYRES))
      CALL DZERO(WORK(KOME2),NT2AM(ISYRES))
COMMENT COMMENT
COMMENT COMMENT
      if (.false.) then
         kt3am  = kend0
         kend0  = kt3am + nt1amx*nt1amx*nt1amx
         lwrk0 = lwork - kend0
         call dzero(work(kt3am),nt1amx*nt1amx*nt1amx)
      endif
COMMENT COMMENT
COMMENT COMMENT
C
      IF (LWRK0 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND0
         CALL QUIT('Insufficient space in CC3_GMAT')
      END IF
C
C---------------------------------------------------------
C     Read the zero'th order amplitudes and calculate
C     the zero'th order Lambda matrices
C---------------------------------------------------------
C
      IOPT   = 1
      CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KT1AM),DUMMY)
C
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),
     &            WORK(KEND0),LWRK0)
C
C--------------------------------------
C     Reorder the t2-amplitudes i T2TP.
C--------------------------------------
C
      IF (LWORK .LT. NT2SQ(1)) THEN
        CALL QUIT('Not enough memory to construct T2TP in CC3_GMAT')
      ENDIF
C
      IOPT = 2
      CALL CC_RDRSP('R0',0,1,IOPT,MODEL,
     *              DUMMY,WORK(KT2TP))
      CALL CC_T2SQ(WORK(KT2TP),WORK(KEND0),ISYMT2)
      CALL CC3_T2TP(WORK(KT2TP),WORK(KEND0),ISYMT2)
C
      IF (IPRINT .GT. 55) THEN
         XNORMVAL = DDOT(NT2SQ(ISYMT2),WORK(KT2TP),1,WORK(KT2TP),1)
         WRITE(LUPRI,*) 'Norm of T2TP ',XNORMVAL
      ENDIF
C
C--------------------------------------
C     Reorder the l2-amplitudes i L2TP.
C--------------------------------------
C
      IF (LWORK .LT. NT2SQ(ISYML2)) THEN
        CALL QUIT('Not enough memory to construct L2TP in CC3_GMAT')
      ENDIF
C
      IOPT = 3
      CALL CC_RDRSP(LISTL,IDLSTL,ISYML2,IOPT,MODEL,
     *              WORK(KL1AM),WORK(KL2TP))
      CALL CC_T2SQ(WORK(KL2TP),WORK(KEND0),ISYML2)
      CALL CC3_T2TP(WORK(KL2TP),WORK(KEND0),ISYML2)
C
      IF (IPRINT .GT. 55) THEN
         XNORMVAL = DDOT(NT2SQ(ISYML2),WORK(KL2TP),1,WORK(KL2TP),1)
         WRITE(LUPRI,*) 'Norm of L2TP ',XNORMVAL
      ENDIF
C
C--------------------------------------
C     Reorder the B2-amplitudes i B2TP.
C--------------------------------------
C
      IF (LWORK .LT. NT2SQ(ISYMB2)) THEN
        CALL QUIT('Not enough memory to construct B2TP in CC3_GMAT')
      ENDIF
C
      IOPT = 3
      CALL CC_RDRSP(LISTB,IDLSTB,ISYMB2,IOPT,MODEL,
     &              WORK(KB1AM),WORK(KB2TP))
      CALL CCLR_DIASCL(WORK(KB2TP),TWO,ISYMB2)
      CALL CC_T2SQ(WORK(KB2TP),WORK(KEND0),ISYMB2)
      CALL CC3_T2TP(WORK(KB2TP),WORK(KEND0),ISYMB2)
C
      IF (IPRINT .GT. 55) THEN
         XNORMVAL = DDOT(NT2SQ(ISYMB2),WORK(KB2TP),1,WORK(KB2TP),1)
         WRITE(LUPRI,*) 'Norm of B2TP ',XNORMVAL
      ENDIF
C
C--------------------------------------
C     Reorder the C2-amplitudes i C2TP.
C--------------------------------------
C
      IF (LWORK .LT. NT2SQ(ISYMC2)) THEN
        CALL QUIT('Not enough memory to construct B2TP in CC3_GMAT')
      ENDIF
C
      IOPT = 3
      CALL CC_RDRSP(LISTC,IDLSTC,ISYMC2,IOPT,MODEL,
     &              WORK(KC1AM),WORK(KC2TP))
      CALL CCLR_DIASCL(WORK(KC2TP),TWO,ISYMC2)
      CALL CC_T2SQ(WORK(KC2TP),WORK(KEND0),ISYMC2)
      CALL CC3_T2TP(WORK(KC2TP),WORK(KEND0),ISYMC2)
C
      IF (IPRINT .GT. 55) THEN
         XNORMVAL = DDOT(NT2SQ(ISYMC2),WORK(KC2TP),1,WORK(KC2TP),1)
         WRITE(LUPRI,*) 'Norm of C2TP ',XNORMVAL
      ENDIF
C

C
C------------------------------------------------------
C     Calculate the (ck|de)-BCtilde and (ck|lm)-BCtilde
C------------------------------------------------------
C
      CALL CC3_3BARINT(ISYMB1,LISTB,IDLSTB,ISYMC1,LISTC,IDLSTC,
     *                IDUMMY,CDUMMY,IDUMMY,.FALSE.,
     *                WORK(KLAMDP),WORK(KLAMDH),WORK(KEND0),LWRK0,
     *                LU3SRTR,FN3SRTR,LUCKJDR,FNCKJDR)
C
      CALL CC3_SORT1(WORK(KEND0),LWRK0,2,ISYMBC,LU3SRTR,FN3SRTR,
     *               LUDELDR,FNDELDR,IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *               IDUMMY,CDUMMY)
C
      CALL CC3_SINT(WORK(KLAMDH),WORK(KEND0),LWRK0,ISYMBC,
     *              LUDELDR,FNDELDR,LUDKBCR,FNDKBCR)
C
      CALL CC3_TCME(WORK(KLAMDP),ISYMBC,WORK(KEND0),LWRK0,
     *              IDUMMY,CDUMMY,LUDKBCR,FNDKBCR,
     *              IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *              IDUMMY,CDUMMY,LUDKBCR4,FNDKBCR4,2)
C
C---------------------------------------------------------------
C     Calculate the g^B integrals for
C     OOOO, OVVO and OOVV integrals. VVVV is stored on file
C---------------------------------------------------------------
C
      KBIOOOO = KEND0
      KBIOVVO = KBIOOOO + N3ORHF(ISINTB) 
      KBIOOVV = KBIOVVO + NT2SQ(ISINTB)
      KEND0   = KBIOOVV + NT2SQ(ISINTB)
      LWRK0   = LWORK   - KEND0
C
      IF (LWRK0 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_GMAT (g^B[2o2v] kind)')
      ENDIF
C
      CALL DZERO(WORK(KBIOOOO),N3ORHF(ISINTB))
      CALL DZERO(WORK(KBIOVVO),NT2SQ(ISINTB))
      CALL DZERO(WORK(KBIOOVV),NT2SQ(ISINTB))
C
      CALL CC3_INT(WORK(KBIOOOO),WORK(KBIOOVV),WORK(KBIOVVO),
     *             ISINTB,LU4VB,FN4VB,
     *             .TRUE.,LISTB,IDLSTB,ISYMB1,
     *             .FALSE.,'DUMMY',IDUMMY,IDUMMY,
     *             WORK(KEND0),LWRK0)
C
C---------------------------------------------------------------
C     Calculate the g^C integrals for
C     OOOO, OVVO and OOVV integrals. VVVV is stored on file
C---------------------------------------------------------------
C
      KCIOOOO = KEND0
      KCIOVVO = KCIOOOO + N3ORHF(ISINTC) 
      KCIOOVV = KCIOVVO + NT2SQ(ISINTC)
      KEND0   = KCIOOVV + NT2SQ(ISINTC)
      LWRK0   = LWORK   - KEND0
C
      IF (LWRK0 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_GMAT (g^C[2o2v] kind)')
      ENDIF
C
      CALL DZERO(WORK(KCIOOOO),N3ORHF(ISINTC))
      CALL DZERO(WORK(KCIOVVO),NT2SQ(ISINTC))
      CALL DZERO(WORK(KCIOOVV),NT2SQ(ISINTC))
C
      CALL CC3_INT(WORK(KCIOOOO),WORK(KCIOOVV),WORK(KCIOVVO),
     *             ISINTC,LU4VC,FN4VC,
     *             .TRUE.,LISTC,IDLSTC,ISYMC1,
     *             .FALSE.,'DUMMY',IDUMMY,IDUMMY,
     *             WORK(KEND0),LWRK0)
C
C---------------------------------------------------------------
C     Calculate the g^BC integrals for
C     OOOO, OVVO and OOVV integrals. VVVV is stored on file
C---------------------------------------------------------------
C
      KBCIOOOO = KEND0
      KBCIOVVO = KBCIOOOO + N3ORHF(ISINTBC)
      KBCIOOVV = KBCIOVVO + NT2SQ(ISINTBC)
      KEND0   = KBCIOOVV + NT2SQ(ISINTBC)
      LWRK0   = LWORK   - KEND0
C
      IF (LWRK0 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_GMAT (g^B[2o2v] kind)')
      ENDIF
C
      CALL DZERO(WORK(KBCIOOOO),N3ORHF(ISINTBC))
      CALL DZERO(WORK(KBCIOVVO),NT2SQ(ISINTBC))
      CALL DZERO(WORK(KBCIOOVV),NT2SQ(ISINTBC))
C
      CALL CC3_INT(WORK(KBCIOOOO),WORK(KBCIOOVV),WORK(KBCIOVVO),
     *             ISINTBC,LU4VBC,FN4VBC,
     *             .TRUE.,LISTB,IDLSTB,ISYMB1,
     *             .TRUE.,LISTC,IDLSTC,ISYMC1,
     *             WORK(KEND0),LWRK0)
C
C---------------------------------------------------------
C     Read canonical orbital energies and MO coefficients.
C---------------------------------------------------------
C
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
C
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KFOCKD+I-1), I=1,NORBTS)
C
      CALL GPCLOSE(LUSIFC,'KEEP')
C
C---------------------------------------------
C     Delete frozen orbitals in Fock diagonal.
C---------------------------------------------
C
      IF (FROIMP .OR. FROEXP)
     *   CALL CCSD_DELFRO(WORK(KFOCKD),WORK(KEND0),LWRK0)
C
C-----------------------------------------------------
C     Construct the transformed Fock matrix
C-----------------------------------------------------
C
      LUFCK = -1
C     This AO Fock matrix is constructed from the T1 transformed density
      CALL GPOPEN(LUFCK,'CC_FCKH','UNKNOWN',' ','UNFORMATTED',
     *              IDUMMY,.FALSE.)
C     This AO Fock matrix is constructed from the CMO transformed density
C      CALL GPOPEN(LUFCK,'CC_FCKREF','UNKNOWN',' ','UNFORMATTED',
C     *            IDUMMY,.FALSE.)
      REWIND(LUFCK)
      READ(LUFCK)(WORK(KFCKBA + I-1),I = 1,N2BST(ISYMOP))
      CALL GPCLOSE(LUFCK,'KEEP' )
C
      IF (IPRINT .GT. 140) THEN
         CALL AROUND( 'Usual Fock AO matrix' )
         CALL CC_PRFCKAO(WORK(KFCKBA),ISYMOP)
      ENDIF
C
      CALL CC_FCKMO(WORK(KFCKBA),WORK(KLAMDP),WORK(KLAMDH),
     *              WORK(KEND0),LWRK0,1,1,1)
C
      IF (IPRINT .GT. 50) THEN
         CALL AROUND( 'In CC3_GMAT: Triples Fock MO matrix' )
         CALL CC_PRFCKMO(WORK(KFCKBA),ISYMOP)
      ENDIF
C
C     Sort the fock matrix
C
C
      CALL DCOPY(N2BST(ISINT1),WORK(KFCKBA),1,WORK(KEND0),1)
C
      DO ISYMC = 1,NSYM
C
         ISYMK = MULD2H(ISYMC,ISINT1)
C
         DO K = 1,NRHF(ISYMK)
C
            DO C = 1,NVIR(ISYMC)
C
               KOFF1 = KEND0 + IFCVIR(ISYMK,ISYMC) + 
     *                 NORB(ISYMK)*(C - 1) + K - 1
               KOFF2 = KFCKBA + IT1AM(ISYMC,ISYMK)
     *               + NVIR(ISYMC)*(K - 1) + C - 1
C
               WORK(KOFF2) = WORK(KOFF1)
C
            ENDDO
         ENDDO
      ENDDO
C
      IF (IPRINT .GT. 50) THEN
         CALL AROUND('In CC3_GMAT: Triples Fock MO matrix (sort)')
         CALL CC_PRFCKMO(WORK(KFCKBA),ISYMOP)
      ENDIF
C
C-----------------------------
C     Read occupied integrals.
C-----------------------------
C
C     Memory allocation.
C
      KTROC  = KEND0
      KTROC0 = KTROC  + NTRAOC(ISINT2)
      KTROC1 = KTROC0 + NTRAOC(ISINT2)
      KTROC2 = KTROC1 + NTRAOC(ISINT2)
      KTROC3 = KTROC2 + NTRAOC(ISINT2)
      KTROC4 = KTROC3 + NTRAOC(ISINTBC)
      KXIAJB = KTROC4 + NTRAOC(ISINTBC)
      KEND1  = KXIAJB + NT2AM(ISYMOP)
      LWRK1  = LWORK  - KEND1
C
      KINTOC = KEND1
      MAXOCC = MAX(NTOTOC(ISINT2),NTOTOC(ISINTBC))
      KEND2  = KINTOC + MAX(NTOTOC(ISYMOP),MAXOCC)
      LWRK2  = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in CC3_LHTR_L3')
      END IF
C
C------------------------
C     Construct L(ia,jb).
C------------------------
C
      LENGTH = IRAT*NT2AM(ISYMOP)
C
      REWIND(LUIAJB)
      CALL READI(LUIAJB,LENGTH,WORK(KXIAJB))
C
      ISYOPE = ISYMOP
      IOPTTCME = 1
      CALL CCSD_TCMEPK(WORK(KXIAJB),1.0D0,ISYOPE,IOPTTCME)
C
      IF ( IPRINT .GT. 55) THEN
         XNORMVAL = DDOT(NT2AM(ISYMOP),WORK(KXIAJB),1,
     *                WORK(KXIAJB),1)
         WRITE(LUPRI,*) 'Norm of IAJB ',XNORMVAL
      ENDIF
C
C--------------------------------------------------------------------------
C     Read in BC-transformed integrals used in contractions and transform.
C--------------------------------------------------------------------------
C
      IOFF = 1
      IF (NTOTOC(ISINTBC) .GT. 0) THEN
         CALL GETWA2(LUCKJDR,FNCKJDR,WORK(KINTOC),IOFF,NTOTOC(ISINTBC))
      ENDIF
C
      CALL CC3_TROCC(WORK(KINTOC),WORK(KTROC3),WORK(KLAMDP),
     *                    WORK(KEND2),LWRK2,ISINTBC)
C
      CALL CCFOP_SORT(WORK(KTROC3),WORK(KTROC4),ISINTBC,1)
C
      CALL CC3_LSORT1(WORK(KTROC3),ISINTBC,WORK(KEND2),LWRK2,5)
C
C------------------------------------------------------------
C     Read in integrals used in contractions and transform.
C------------------------------------------------------------
C
      IOFF = 1
      IF (NTOTOC(ISINT2) .GT. 0) THEN
         CALL GETWA2(LUCKJD,FNCKJD,WORK(KINTOC),IOFF,NTOTOC(ISINT2))
      ENDIF
C
      CALL CC3_TROCC(WORK(KINTOC),WORK(KTROC),WORK(KLAMDP),
     *                    WORK(KEND2),LWRK2,ISINT2)
C
      CALL CCFOP_SORT(WORK(KTROC),WORK(KTROC1),ISINT2,1)
C
      CALL CC3_LSORT1(WORK(KTROC),ISINT2,WORK(KEND2),LWRK2,5)
C
C------------------------------------------------------------
C     Read in integrals used in L3 and transform.
C------------------------------------------------------------
C
      IOFF = 1
      IF (NTOTOC(ISINT2) .GT. 0) THEN
         CALL GETWA2(LUTOC,FNTOC,WORK(KINTOC),IOFF,NTOTOC(ISINT2))
      ENDIF
C
      CALL CC3_TROCC(WORK(KINTOC),WORK(KTROC0),WORK(KLAMDH),
     *                 WORK(KEND2),LWRK2,ISINT2)
C
C-----------------------------------------------------------
C     Construct 2*C-E of the integrals.
C-----------------------------------------------------------
C
      CALL CCSDT_TCMEOCC(WORK(KTROC0),WORK(KTROC2),ISINT2)
C
C-------------------------------
C     Write out norms of arrays.
C-------------------------------
C
      IF (IPRINT .GT. 55) THEN
         XNORMVAL = DDOT(NTOTOC(ISINT2),WORK(KINTOC),1,
     *                WORK(KINTOC),1)
         WRITE(LUPRI,*) 'Norm of CKJDEL-INT  ',XNORMVAL
      ENDIF
C
      IF (IPRINT .GT. 55) THEN
         XNORMVAL = DDOT(NTRAOC(ISINT2),WORK(KTROC0),1,
     *                WORK(KTROC0),1)
         WRITE(LUPRI,*) 'Norm of TROC0 CC3_GMAT : ',XNORMVAL
      ENDIF
C
      IF (IPRINT .GT. 55) THEN
         XNORMVAL = DDOT(NTRAOC(ISINT2),WORK(KTROC2),1,
     *                WORK(KTROC2),1)
         WRITE(LUPRI,*) 'Norm of TROC2 CC3_GMAT : ',XNORMVAL
      ENDIF
C
C----------------------------
C     General loop structure.
C----------------------------
C
      DO ISYMD = 1,NSYM
C
         ISAIJ1 = MULD2H(ISYMD,ISYRES)
         ISYCKB = MULD2H(ISYMD,ISYMOP)
         ISCKB1 = MULD2H(ISINT1,ISYMD)
         ISCKB2 = MULD2H(ISINT2,ISYMD)

         ISCKBB = MULD2H(ISINTB,ISYMD)
         ISCKBC = MULD2H(ISINTC,ISYMD)
         ISCKBBC = MULD2H(ISINTBC,ISYMD)
C
         IF (IPRINT .GT. 55) THEN
C
            WRITE(LUPRI,*) 'In CC3_GMAT: ISAIJ1 :',ISAIJ1
            WRITE(LUPRI,*) 'In CC3_GMAT: ISYCKB :',ISYCKB
            WRITE(LUPRI,*) 'In CC3_GMAT: ISCKB2 :',ISCKB2

            WRITE(LUPRI,*) 'In CC3_GMAT: ISCKBB :',ISCKBB
            WRITE(LUPRI,*) 'In CC3_GMAT: ISCKBC :',ISCKBC
            WRITE(LUPRI,*) 'In CC3_GMAT: ISCKBBC :',ISCKBBC
C
         ENDIF
C
C--------------------------
C        Memory allocation.
C--------------------------
C 
         KTRVI7 = KEND1 
         KTRVI8 = KTRVI7 + NCKATR(ISCKBBC)
         KTRVI2 = KTRVI8 + NCKATR(ISCKBBC)
         KRMAT1 = KTRVI2 + NCKATR(ISCKB2)
         KEND2  = KRMAT1 + NCKI(ISAIJ1)
         LWRK2  = LWORK  - KEND2
C
         KTRVI0  = KEND2
         KTRVI3  = KTRVI0  + NCKATR(ISCKB2)
         KTRVI4  = KTRVI3  + NCKATR(ISCKB2)
         KTRVI5  = KTRVI4  + NCKATR(ISCKB2)
         KTRVI6  = KTRVI5  + NCKATR(ISCKB2)
         KVVVVB  = KTRVI6  + NCKATR(ISCKB2)
         KVVVVC  = KVVVVB  + NMAABC(ISCKBB)
         KVVVVBC = KVVVVC  + NMAABC(ISCKBC)
         KEND3   = KVVVVBC  + NMAABC(ISCKBBC)
         LWRK3   = LWORK  - KEND3
C
         KINTVI = KEND3
         KEND4  = KINTVI + MAX(NCKA(ISYMD),NCKA(ISCKB2))
         LWRK4  = LWORK  - KEND4
C
         IF (LWRK4 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND4
            CALL QUIT('Insufficient space in CC3_LHTR_L3')
         END IF
C
C---------------------
C        Sum over D
C---------------------
C
         DO D = 1,NVIR(ISYMD)
C
C------------------------------------
C           Initialize the R1 matrix.
C------------------------------------
C
            CALL DZERO(WORK(KRMAT1),NCKI(ISAIJ1))
C
C--------------------------------------------
C           Read in g^B_{vvvv} for a given D
C--------------------------------------------
C
            IF (NMAABC(ISCKBB) .GT. 0) THEN
               IOFF = I3VVIR(ISCKBB,ISYMD)
     *              + NMAABC(ISCKBB)*(D-1)
     *              + 1
               CALL GETWA2(LU4VB,FN4VB,WORK(KVVVVB),IOFF,NMAABC(ISCKBB))
            ENDIF
C
C--------------------------------------------
C           Read in g^B_{vvvv} for a given D
C--------------------------------------------
C
            IF (NMAABC(ISCKBC) .GT. 0) THEN
               IOFF = I3VVIR(ISCKBC,ISYMD)
     *              + NMAABC(ISCKBC)*(D-1)
     *              + 1
               CALL GETWA2(LU4VC,FN4VC,WORK(KVVVVC),IOFF,NMAABC(ISCKBC))
            ENDIF
C
C
C--------------------------------------------
C           Read in g^B_{vvvv} for a given D
C--------------------------------------------
C
            IF (NMAABC(ISCKBBC) .GT. 0) THEN
               IOFF = I3VVIR(ISCKBBC,ISYMD)
     *              + NMAABC(ISCKBBC)*(D-1)
     *              + 1
               CALL GETWA2(LU4VBC,FN4VBC,WORK(KVVVVBC),IOFF,
     *                     NMAABC(ISCKBBC))
            ENDIF
C
C--------------------------------------------------------------
C           Read B-transformed integrals used in contraction
C--------------------------------------------------------------
C
            IF (NCKATR(ISCKBBC) .GT. 0) THEN
               IOFF = ICKBD(ISCKBBC,ISYMD) + NCKATR(ISCKBBC)*(D - 1) + 1
               CALL GETWA2(LUDKBCR,FNDKBCR,WORK(KTRVI7),IOFF,
     &                     NCKATR(ISCKBBC))
            ENDIF
C
            IF (LWRK4 .LT. NCKATR(ISCKBBC)) THEN
               CALL QUIT('Insufficient space for allocation in '//
     &                   'CC3_GMAT (TRVI7)')
            END IF
C
            CALL CCSDT_SRVIR3(WORK(KTRVI7),WORK(KEND4),ISYMD,D,ISINTBC)
C
            IF (NCKATR(ISCKBBC) .GT. 0) THEN
               IOFF = ICKBD(ISCKBBC,ISYMD) + NCKATR(ISCKBBC)*(D - 1) + 1
               CALL GETWA2(LUDKBCR4,FNDKBCR4,WORK(KTRVI8),IOFF,
     &                     NCKATR(ISCKBBC))
            ENDIF
C
            IF (LWRK4 .LT. NCKATR(ISCKBBC)) THEN
               CALL QUIT('Insufficient space for allocation in '//
     &                   'CC3_GMAT (TRVI8)')
            END IF
C
            CALL CCSDT_SRVIR3(WORK(KTRVI8),WORK(KEND4),ISYMD,D,ISINTBC)
C
C-----------------------------------------------
C           Integrals used in L3am.
C-----------------------------------------------
C
            IOFF = ICKBD(ISCKB2,ISYMD) + NCKATR(ISCKB2)*(D - 1) + 1
            IF (NCKATR(ISCKB2) .GT. 0) THEN
               CALL GETWA2(LUDKBC3,FNDKBC3,WORK(KTRVI0),IOFF,
     &                     NCKATR(ISCKB2))
            ENDIF
C
            CALL CCSDT_SRVIR3(WORK(KTRVI0),WORK(KEND4),ISYMD,D,ISINT2)
C
C------------------------------------------------------
C           Read 2*C-E of integral used for L3
C------------------------------------------------------
C
            IOFF = ICKBD(ISCKB2,ISYMD) + NCKATR(ISCKB2)*(D - 1) + 1
            IF (NCKATR(ISCKB2) .GT. 0) THEN
               CALL GETWA2(LU3FOP2X,FN3FOP2X,WORK(KTRVI4),IOFF,
     &                     NCKATR(ISCKB2))
            ENDIF
C
C-----------------------------------------------------------
C           Sort the integrals for s3am and for t3-bar
C-----------------------------------------------------------
C
            CALL CCSDT_SRTVIR(WORK(KTRVI0),WORK(KTRVI2),WORK(KEND4),
     *                        LWRK4,ISYMD,ISINT2)
C
            CALL CCSDT_SRTVIR(WORK(KTRVI4),WORK(KTRVI5),WORK(KEND4),
     *                        LWRK4,ISYMD,ISINT2)
C
C
            IF (IPRINT .GT. 55) THEN
               XNORMVAL = DDOT(NCKATR(ISCKB2),WORK(KTRVI0),1,
     *                      WORK(KTRVI0),1)
               WRITE(LUPRI,*) 'Norm of TRVI0 ',XNORMVAL
            ENDIF
C
            IF (IPRINT .GT. 55) THEN
               XNORMVAL = DDOT(NCKATR(ISCKB2),WORK(KTRVI2),1,
     *                      WORK(KTRVI2),1)
               WRITE(LUPRI,*) 'Norm of TRVI2 ',XNORMVAL
            ENDIF
C
            IF (IPRINT .GT. 55) THEN
               XNORMVAL = DDOT(NCKATR(ISCKB2),WORK(KTRVI4),1,
     *                      WORK(KTRVI4),1)
               WRITE(LUPRI,*) 'Norm of TRVI4 ',XNORMVAL
            ENDIF
C
            IF (IPRINT .GT. 55) THEN
               XNORMVAL = DDOT(NCKATR(ISCKB2),WORK(KTRVI5),1,
     *                      WORK(KTRVI5),1)
               WRITE(LUPRI,*) 'Norm of TRVI5 ',XNORMVAL
            ENDIF
C
C------------------------------------------------------
C           Calculate virtual integrals used in q3am.
C------------------------------------------------------
C
            IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D - 1) + 1
            IF (NCKA(ISYCKB) .GT. 0) THEN
               CALL GETWA2(LU3VI,FN3VI,WORK(KINTVI),IOFF,
     *                     NCKA(ISYCKB))
            ENDIF
C
            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI3),WORK(KLAMDP),
     *                       ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
C
            IF (LWRK3 .LT. NCKATR(ISYCKB)) THEN
               CALL QUIT('Insufficient space for allocation in '//
     *                   'CC3_GMAT  (TRVI3)')
            END IF
C
C           Can use kend3 since dont need the integrals anymore
            CALL CCSDT_SRVIR3(WORK(KTRVI3),WORK(KEND3),ISYMD,D,ISINT2)
C
C---------------------------------------------------------------
C           Read virtual integrals used in q3am/u3am for t3-bar.
C---------------------------------------------------------------
C
            IOFF = ICKBD(ISYCKB,ISYMD) + NCKATR(ISYCKB)*(D - 1) + 1
            IF (NCKATR(ISYCKB) .GT. 0) THEN
               CALL GETWA2(LU3FOPX,FN3FOPX,WORK(KTRVI6),IOFF,
     *                     NCKATR(ISYCKB))
            ENDIF
C
            IF (LWRK3 .LT. NCKATR(ISYCKB)) THEN
               CALL QUIT('Insufficient space for allocation in '//
     &                   'CC3_GMAT (TRVI6)')
            END IF
C
C           Can use kend3 since dont need the integrals anymore
            CALL CCSDT_SRVIR3(WORK(KTRVI6),WORK(KEND4),ISYMD,D,ISINT2)
C
            IF (IPRINT .GT. 55) THEN
               XNORMVAL = DDOT(NCKATR(ISCKB2),WORK(KTRVI6),1,
     *                      WORK(KTRVI6),1)
               WRITE(LUPRI,*) 'Norm of TRVI6 ',XNORMVAL
            ENDIF
C
C---------------------
C           Calculate.
C---------------------
C
            DO ISYMB = 1,NSYM
C
               ISYALJ  = MULD2H(ISYMB,ISYML2) 
               ISAIJ2  = MULD2H(ISYMB,ISYRES) 
               ISYMBD  = MULD2H(ISYMB,ISYMD)
               ISCKIJ  = MULD2H(ISYMBD,ISYMIM)
C
               IF (IPRINT .GT. 55) THEN
C
                  WRITE(LUPRI,*) 'In CC3_LHTR_L3: ISYMD :',ISYMD
                  WRITE(LUPRI,*) 'In CC3_LHTR_L3: ISYMB :',ISYMB
                  WRITE(LUPRI,*) 'In CC3_LHTR_L3: ISYALJ:',ISYALJ
                  WRITE(LUPRI,*) 'In CC3_LHTR_L3: ISAIJ2:',ISAIJ2
                  WRITE(LUPRI,*) 'In CC3_LHTR_L3: ISYMBD:',ISYMBD
                  WRITE(LUPRI,*) 'In CC3_LHTR_L3: ISCKIJ:',ISCKIJ
C
               ENDIF
C
C              Can use kend3 since we do not need the integrals anymore.
               KSMAT   = KEND3
               KQMAT   = KSMAT   + NCKIJ(ISCKIJ)
               KDIAG   = KQMAT   + NCKIJ(ISCKIJ)
               KINDSQ  = KDIAG   + NCKIJ(ISCKIJ)
               KINDEX  = KINDSQ  + (6*NCKIJ(ISCKIJ) - 1)/IRAT + 1
               KTMAT   = KINDEX  + (NCKI(ISYALJ) - 1)/IRAT + 1
               KRMAT2  = KTMAT   + NCKIJ(ISCKIJ)
               KEND4   = KRMAT2  + NCKI(ISAIJ2)
               LWRK4   = LWORK   - KEND4
C
               IF (LWRK4 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Memory available : ',LWORK
                  WRITE(LUPRI,*) 'Memory needed    : ',KEND4
                  CALL QUIT('Insufficient space in CC3_LHTR_L3 (inner)')
               END IF
C
C---------------------------------------------
C              Construct part of the diagonal.
C---------------------------------------------
C
               CALL CC3_DIAG(WORK(KDIAG),WORK(KFOCKD),ISCKIJ)
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
C
               DO B = 1,NVIR(ISYMB)
C
C-----------------------------------------
C                 Initialize the R2 matrix.
C-----------------------------------------
C
                  CALL DZERO(WORK(KRMAT2),NCKI(ISAIJ2))
C
C---------------------------------------------------------------------------
C                 Calculate the S(ci,bk,dj) matrix for for B,D for L3
C                 Scale the smat with a half since have a factor
C                 of two in the routine.
C---------------------------------------------------------------------------
C
C
                  CALL DZERO(WORK(KSMAT),NCKIJ(ISCKIJ))
C
                  CALL CCFOP_SMAT(0.0D0,WORK(KL1AM),ISYML1,WORK(KL2TP),
     *                            ISYML2,WORK(KTMAT),
     *                            WORK(KFCKBA),WORK(KXIAJB),ISINT1,
     *                            WORK(KTRVI0),WORK(KTRVI2),
     *                            WORK(KTRVI4),WORK(KTRVI5),
     *                            WORK(KTROC0),WORK(KTROC2),
     *                            ISINT2,WORK(KFOCKD),WORK(KDIAG),
     *                            WORK(KSMAT),WORK(KEND4),LWRK4,
     *                            WORK(KINDEX),WORK(KINDSQ),LENSQ,
     *                            ISYMB,B,ISYMD,D)
C
                  CALL DSCAL(NCKIJ(ISCKIJ),HALF,WORK(KSMAT),1)
C
                  CALL T3_FORBIDDEN(WORK(KSMAT),ISYMIM,ISYMB,B,ISYMD,D)
COMMENT COMMENT
COMMENT COMMENT
C      call sum_pt3(work(ksmat),isymb,b,isymd,d,
C     *             isckij,work(kt3am),1)
COMMENT COMMENT
COMMENT COMMENT
C
C
                  IF (IPRINT .GT. 55) THEN
                     XNORMVAL = DDOT(NCKIJ(ISCKIJ),WORK(KSMAT),1,
     *                       WORK(KSMAT),1)
                     WRITE(LUPRI,*) 'Norm of SMAT-L3    ',XNORMVAL
                  ENDIF
C
C-------------------------------------------------------------------
C                 Calculate Q(ci,jk) for fixed b,d for t3-bar.
C-------------------------------------------------------------------
C
C
                  CALL DZERO(WORK(KQMAT),NCKIJ(ISCKIJ))
C
                  ! g integrals
c                 call dzero(WORK(KTRVI3),NCKATR(ISCKB2))
                  ! L integrals
c                 call dzero(WORK(KTRVI6),NCKATR(ISCKB2))
C
                  ! g integrals
c                 call dzero(WORK(KTROC0),NTRAOC(ISINT2))
c                 ! L integrals
c                 call dzero(WORK(KTROC2),NTRAOC(ISINT2))
c
c                 call dzero(WORK(KXIAJB),NT2AM(ISYMOP))
c                 call dzero(WORK(KFCKBA),NT1AM(1))
                  CALL CCFOP_QMAT(0.0D0,WORK(KL1AM),ISYML1,
     *                            WORK(KL2TP),ISYML2,
     *                            WORK(KTMAT),WORK(KFCKBA),
     *                            WORK(KXIAJB),ISINT1,WORK(KTRVI3),
     *                            WORK(KTRVI6),WORK(KTROC0),
     *                            WORK(KTROC2),ISINT2,WORK(KFOCKD),
     *                            WORK(KDIAG),WORK(KQMAT),
     *                            WORK(KEND4),LWRK4,WORK(KINDSQ),
     *                            LENSQ,ISYMB,B,ISYMD,D)
C
                  CALL DSCAL(NCKIJ(ISCKIJ),HALF,WORK(KQMAT),1)
C
                  CALL T3_FORBIDDEN(WORK(KQMAT),ISYMIM,ISYMB,B,ISYMD,D)
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
c      call sum_pt3(work(kqmat),isymb,b,isymd,d,
c    *             isckij,work(kt3am),2)
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
C
C
                  IF (IPRINT .GT. 55) THEN
                     XNORMVAL = DDOT(NCKIJ(ISCKIJ),WORK(KQMAT),1,
     *                       WORK(KQMAT),1)
                     WRITE(LUPRI,*) 'Norm of QMAT-L3  ',XNORMVAL
                  ENDIF
C
C-----------------------------------------------------------------------
C                 Calculate the contributions to omega2
C                         <L3|[H^BC,tau_{2}]|HF>
C-----------------------------------------------------------------------
c     WRITE(lupri,*)'ome2 before CFOP_CONOCC'
c     call output(WORK(KOME2),1,nt1am(1),1,nt1am(1),
c    *            nt1am(1),nt1am(1),1,lupri)
C
                  CALL CCFOP_CONVIR(WORK(KRMAT2),WORK(KSMAT),
     *                              WORK(KQMAT),WORK(KTMAT),ISYMIM,
     *                              WORK(KTRVI7),WORK(KTRVI8),ISINTBC,
     *                              WORK(KEND4),LWRK4,WORK(KINDSQ),
     *                              LENSQ,ISYMB,B,ISYMD,D)
c      write(lupri,*)'ISYMB,B,ISYMD,D,ISINTBC',ISYMB,B,ISYMD,D,ISINTBC
c      xnormval = ddot(NCKIJ(ISCKIJ),WORK(KSMAT),1,WORK(KSMAT),1)
c      write(lupri,*)'KSMAT',xnormval
c      xnormval = ddot(NCKIJ(ISCKIJ),WORK(KQMAT),1,WORK(KQMAT),1)
c      write(lupri,*)'KQMAT',xnormval
c      xnormval = ddot(NCKATR(ISCKBBC),WORK(KTRVI7),1,WORK(KTRVI7),1)
c      write(lupri,*)'KTRVI7,ISINTBC',xnormval,ISINTBC
c      xnormval = ddot(NCKATR(ISCKBBC),WORK(KTRVI8),1,WORK(KTRVI8),1)
c      write(lupri,*)'KTRVI8',xnormval
c      xnormval =  ddot(NTRAOC(ISINTBC),WORK(KTROC3),1,WORK(KTROC3),1)
c      write(lupri,*)'KTROC3',xnormval
c      xnormval =  ddot(NTRAOC(ISINTBC),WORK(KTROC4),1,WORK(KTROC4),1)
c      write(lupri,*)'KTROC4',xnormval
C
                  CALL CCFOP_CONOCC(WORK(KOME2),WORK(KRMAT1),
     *                              WORK(KRMAT2),WORK(KSMAT),
     *                              WORK(KTMAT),ISYMIM,
     *                              WORK(KTROC3),WORK(KTROC4),ISINTBC,
     *                              WORK(KEND4),LWRK4,WORK(KINDSQ),
     *                              LENSQ,ISYMB,B,ISYMD,D)
c              XNORMVAL = DDOT(NT2AM(ISYRES),WORK(KOME2),1,
c    *                         WORK(KOME2),1)
c              WRITE(LUPRI,*) 'Norm  Rho22-after CCFOP_CONOCC',XNORMVAL
c     WRITE(lupri,*)'ome2 after CFOP_CONOCC'
c     call output(WORK(KOME2),1,nt1am(1),1,nt1am(1),
c    *            nt1am(1),nt1am(1),1,lupri)

C
C------------------------------------------------------------------------
C                 Calculate the L3 contribution to omega1
C------------------------------------------------------------------------
C
                  CALL DSCAL(NCKIJ(ISCKIJ),-ONE,WORK(KSMAT),1)
                  CALL DSCAL(NCKIJ(ISCKIJ),-ONE,WORK(KQMAT),1)
C
C                        <L3|[H^BC,T^{0}_{2}],tau_{1}]|HF>
c             write(lupri,*)'ISYMB,B,ISYMD,D,ISINTBC'
c             write(lupri,*)ISYMB,B,ISYMD,D,ISINTBC
C   
c             XNORMVAL = DDOT(N3ORHF(ISINTBC),WORK(KBCIOOOO),1,
c    *                                        WORK(KBCIOOOO),1)
c             write(lupri,*)'KBCIOOOO',xnormval
c             XNORMVAL = DDOT(NT2SQ(ISINTBC),WORK(KBCIOVVO),1,
c    *                                        WORK(KBCIOVVO),1)
c             write(lupri,*)'KBCIOVVO',xnormval
c             XNORMVAL = DDOT(NT2SQ(ISINTBC),WORK(KBCIOOVV),1,
c    *                                        WORK(KBCIOOVV),1)
c             write(lupri,*)'KBCIOOVV',xnormval
c             XNORMVAL = DDOT(NMAABC(ISCKBBC),WORK(KVVVVBC),1,
c    *                                        WORK(KVVVVBC),1)
c             write(lupri,*)'KVVVVBC',xnormval
                  !<L3|[H^BC,T^{0}_{2}],tau_{1}]|HF>
                  CALL CC3_L3_OMEGA1(WORK(KOME1),ISYRES,WORK(KSMAT),
     *                               WORK(KQMAT),WORK(KTMAT),ISYMIM,
     *                               WORK(KBCIOOOO),WORK(KBCIOVVO),
     *                               WORK(KBCIOOVV),WORK(KVVVVBC),
     *                               ISINTBC,
     *                               WORK(KT2TP),ISYMT2,WORK(KEND4),
     *                               LWRK4,LENSQ,WORK(KINDSQ),
     *                               ISYMB,B,ISYMD,D)
C
c     WRITE(lupri,*)'ome1 after CC3_L3_OMEGA2_NODDY t20 cont'
c     call output(WORK(KOME1),1,nvir(1),1,nrhf(1),
c    *            nvir(1),nrhf(1),1,lupri)

                  IF (IPRINT .GT. 55) THEN
                     XNORMVAL = DDOT(NT1AM(ISYRES),OMEGA1,1,OMEGA1,1)
                     WRITE(LUPRI,*) 'Norm CC3_GMAT after CC3_L3_OMEGA1',
     *                               XNORMVAL
                  ENDIF
C
C                        <L3|[H^B,T^{C}_{2}],tau_{1}]|HF>
C   
C
                  CALL CC3_L3_OMEGA1(WORK(KOME1),ISYRES,WORK(KSMAT),
     *                               WORK(KQMAT),WORK(KTMAT),ISYMIM,
     *                               WORK(KBIOOOO),WORK(KBIOVVO),
     *                               WORK(KBIOOVV),WORK(KVVVVB),ISINTB,
     *                               WORK(KC2TP),ISYMC2,WORK(KEND4),
     *                               LWRK4,LENSQ,WORK(KINDSQ),
     *                               ISYMB,B,ISYMD,D)
c     WRITE(lupri,*)'ome1 after CC3_L3_OMEGA2_NODDY t2c cont'
c     call output(WORK(KOME1),1,nvir(1),1,nrhf(1),
c    *            nvir(1),nrhf(1),1,lupri)

C
                  IF (IPRINT .GT. 55) THEN
                     XNORMVAL = DDOT(NT1AM(ISYRES),OMEGA1,1,OMEGA1,1)
                     WRITE(LUPRI,*) 'Norm CC3_GMAT after CC3_L3_OMEGA1',
     *                               XNORMVAL
                  END IF

C

C                        <L3|[H^C,T^{B}_{2}],tau_{1}]|HF>
C 
                  CALL CC3_L3_OMEGA1(WORK(KOME1),ISYRES,WORK(KSMAT),
     *                               WORK(KQMAT),WORK(KTMAT),ISYMIM,
     *                               WORK(KCIOOOO),WORK(KCIOVVO),
     *                               WORK(KCIOOVV),WORK(KVVVVC),ISINTC,
     *                               WORK(KB2TP),ISYMB2,WORK(KEND4),
     *                               LWRK4,LENSQ,WORK(KINDSQ),
     *                               ISYMB,B,ISYMD,D)
c     WRITE(lupri,*)'ome1 after CC3_L3_OMEGA2_NODDY t2b cont'
c     call output(WORK(KOME1),1,nvir(1),1,nrhf(1),
c    *            nvir(1),nrhf(1),1,lupri)

C
                  IF (IPRINT .GT. 55) THEN
                     XNORMVAL = DDOT(NT1AM(ISYRES),OMEGA1,1,OMEGA1,1)
                     WRITE(LUPRI,*) 'Norm CC3_GMAT after CC3_L3_OMEGA1',
     *                               XNORMVAL
                  ENDIF
C
C-------------------------------------------
C                 Accumulate R2 in Omega2
C-------------------------------------------
C
                  CALL CC3_RACC(WORK(KOME2),WORK(KRMAT2),ISYMB,B,
     *                          ISYRES)
c     WRITE(lupri,*)'ome2 after cc3_racc  rmat2'
c     call output(WORK(KOME2),1,nt1am(1),1,nt1am(1),
c    *            nt1am(1),nt1am(1),1,lupri)
c     WRITE(lupri,*)'ome2 triangular form '
c     call outpak(WORK(KOME2),nt1amx,1,lupri)
C
                  IF (IPRINT .GT. 55) THEN
                     XNORMVAL = DDOT(NT2AM(ISYRES),WORK(KOME2),1,
     *                               WORK(KOME2),1)
                     WRITE(LUPRI,*) 'Norm of Rho22-after CC3_RACC',
     *               XNORMVAL
                  ENDIF
C
                  IF (IPRINT .GT. 220) THEN
                     CALL AROUND('After CC3_RACC: ')
                     CALL CC_PRP(DUMMY,OMEGA2,ISYRES,0,1)
                  ENDIF
C
               ENDDO   ! B
            ENDDO      ! ISYMB
C
C---------------------------------------
C           Accumulate R1 in Omega2.
C---------------------------------------
C
            CALL CC3_RACC(WORK(KOME2),WORK(KRMAT1),ISYMD,D,ISYRES)
c     WRITE(lupri,*)'ome2 after cc3_racc  rmat1'
c     call output(WORK(KOME2),1,nt1am(1),1,nt1am(1),
c    *            nt1am(1),nt1am(1),1,lupri)
C
            IF (IPRINT .GT. 55) THEN
               XNORMVAL = DDOT(NT2AM(ISYRES),WORK(KOME2),1,
     *                         WORK(KOME2),1)
               WRITE(LUPRI,*) 'Norm of Rho22-after CC3_RACC-2',XNORMVAL
            ENDIF
C
            IF (IPRINT .GT. 220) THEN
               CALL AROUND('After CC3_RACC-2: ')
               CALL CC_PRP(DUMMY,WORK(KOME2),ISYRES,0,1)
            ENDIF
C
         ENDDO       ! D
      ENDDO          ! ISYMD
C
COMMENT COMMENT COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT COMMENT COMMENT
c      write(lupri,*) 'L3 Amplitudes in cc3_Gmat : '
C      call dscal(nt1amx*nt1amx*nt1amx,-1.0D0,work(kt3am),1)
c      call print_pt3(work(kt3am),isckij,3)
COMMENT COMMENT COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT COMMENT COMMENT
C-------------------------------
C     Close and delete files
C-------------------------------
C
      CALL WCLOSE2(LU4VBC,FN4VBC,'DELETE')
      CALL WCLOSE2(LU4VC,FN4VC,'DELETE')
      CALL WCLOSE2(LU4VB,FN4VB,'DELETE')
      CALL WCLOSE2(LU3SRTR,FN3SRTR,'DELETE')
      CALL WCLOSE2(LUCKJDR,FNCKJDR,'DELETE')
      CALL WCLOSE2(LUDELDR,FNDELDR,'DELETE')
      CALL WCLOSE2(LUDKBCR,FNDKBCR,'DELETE')
      CALL WCLOSE2(LUDKBCR4,FNDKBCR4,'DELETE')
C
COMMENT COMMENT COMMENT COMMENT 
C       write(lupri,*) 'OMEGA1EFF in cc3_Gmat before added'
C       do i = 1,nt1am(isyres)
C           write(lupri,*) i,OMEGA1EFF(i)
C       end do
COMMENT COMMENT COMMENT COMMENT 
C     write(lupri,*) 'OMEGA2EFF in cc3_Gmat before added'
C     call outpak(OMEGA2EFF,NT2AM(isyres),1,lupri)
COMMENT COMMENT COMMENT COMMENT 

C
C     WRITE(LUPRI,*)' OME1 FINAL'
C     CALL PRINT_MATAI(WORK(KOME1),ISYRES)
C     XNORMVAL = DDOT(NT1AM(ISYRES),WORK(KOME1),1,
C    *                WORK(KOME1),1)
C     WRITE(LUPRI,*) 'Norm of OME1',XNORMVAL
C
      DO I = 1,NT1AM(ISYRES)
         OMEGA1(I) = OMEGA1(I) + WORK(KOME1-1+I)
      END DO
C
c     WRITE(LUPRI,*)' OMEGA1 FINAL'
c     CALL PRINT_MATAI(OMEGA1,ISYRES)
c     XNORMVAL = DDOT(NT1AM(ISYRES),OMEGA1,1,
c    *                     OMEGA1,1)
c     WRITE(LUPRI,*) 'Norm of OMEGA1',XNORMVAL
 
c     WRITE(lupri,*)'ome2 triangular form final '
c     call outpak(WORK(KOME2),nt1amx,1,lupri)
C     XNORMVAL = DDOT(NT2AM(ISYRES),WORK(KOME2),1,
C    *                WORK(KOME2),1)
C     WRITE(LUPRI,*) 'Norm of ome2 ',XNORMVAL
C
      DO I = 1,NT2AM(ISYRES)
         OMEGA2(I) = OMEGA2(I) + WORK(KOME2-1+I)
      END DO
C
c     WRITE(lupri,*)'omega2 triangular form final '
c     call outpak(OMEGA2,nt1amx,1,lupri)
C     XNORMVAL = DDOT(NT2AM(ISYRES),OMEGA2,1,
C    *                     OMEGA2,1)
C     WRITE(LUPRI,*) 'Norm of OMEGA2',XNORMVAL
C
COMMENT COMMENT COMMENT COMMENT 
C       write(lupri,*) 'OMEGA1EFF in cc3_Gmat after added'
C       do i = 1,nt1am(isyres)
C           write(lupri,*) i,OMEGA1EFF(i)
C       end do
COMMENT COMMENT COMMENT COMMENT 
C     write(lupri,*) 'OMEGA2EFF in cc3_Gmat after added'
C     call outpak(OMEGA2EFF,NT2AM(isyres),1,lupri)
COMMENT COMMENT COMMENT COMMENT 

C
C-------------
C     End
C-------------
C
      CALL QEXIT('CC3_GMAT')
C
      RETURN
C
    1 FORMAT(7X,'Time used in',2X,A12,F12.2,' seconds')
C
      END
C  /* Deck cc3_3barint */
      SUBROUTINE CC3_3BARINT(ISYMX,LISTX,IDLSTX,
     *                       ISYMY,LISTY,IDLSTY,
     *                       ISYMZ1,LISTZ,IDLSTZ,Z1EXIST,
     *                       XLAMDP,XLAMDH,WORK,LWORK,
     *                       LU3SRTR,FN3SRTR,LUCKJDR,FNCKJDR)
C
C
C     Interface to calculate 
C     the (ck|de)-X,Y,Z bar and (ck|lm)-X,Y,Z bar integrals
C
C     e and l index is normal T1-transformed (and can be transformed outside)
C
C     (alfa beta /gamma  del) is transformed to  (ck/d del) and (ck/del m)
C     
C     if Z1EXIST eq true
C       one indes transform with z1am
C     else
C       use XLAMDP,XLAMDH as transformation matrix
C     endif 
C
C      The atomic integrals (alfa beta /gamma  del) have symmetry isymop
C
      IMPLICIT NONE
C
#include <ccsdsym.h>
#include <priunit.h>
#include <dummy.h>
#include <maxorb.h>
#include <maxash.h>
#include <mxcent.h>
#include <aovec.h>
#include <iratdef.h>
#include <ccorb.h>
#include <ccisao.h>
#include <r12int.h>
#include <blocks.h>
#include <ccsdinp.h>
#include <distcl.h>
#include <cbieri.h>
#include <eribuf.h>
C#include <cclr.h>
C
      INTEGER ISYMX,ISYMY,ISYMZ,ISYMZ1,LWORK,LU3SRTR,LUCKJDR 
      INTEGER KLAMPX,KLAMHX,KLAMPY,KLAMHY,KLAMPZ,KEND1,LWRK1
      INTEGER KCCFB1, KINDXB, KODCL1, KODCL2, KODBC1, KODBC2
      INTEGER KRDBC1, KRDBC2, KODPP1, KODPP2, KRDPP1, KRDPP2
      INTEGER KFREE,LFREE
      INTEGER NTOSYM,NTOT,ILLL,NUMDIS,KRECNR,IDEL2,IDEL,ISYMD
      INTEGER ISYDIS,ISYAIK,ISYCKJ,KXINT,KCKJ,KEND2,LWRK2
      INTEGER KLAMHZ,ISYMD1
      INTEGER IDLSTX,IDLSTY,IDLSTZ
      INTEGER IDUM
C
      INTEGER INDEXA(MXCORB)
#if defined (SYS_CRAY)
      REAL XLAMDP(*),XLAMDH(*)
      REAL WORK(LWORK)
      real ddot
#else
      DOUBLE PRECISION XLAMDP(*),XLAMDH(*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION ddot 
#endif
C
      CHARACTER*3 LISTX, LISTY, LISTZ
      CHARACTER*(*) FN3SRTR, FNCKJDR
      CHARACTER*1 CDUMMY
C
      LOGICAL Z1EXIST
C
      CALL QENTER('CC3_3BARINT')
C
      CDUMMY = ' '
C
      KLAMPX = 1
      KLAMHX = KLAMPX + NLAMDT
      KLAMPY = KLAMHX + NLAMDT
      KLAMHY = KLAMPY + NLAMDT
      KLAMPZ = KLAMHY + NLAMDT
      KLAMHZ = KLAMPZ + NLAMDT
      KEND1  = KLAMHZ + NLAMDT
      LWRK1  = LWORK  - KEND1
C
C
C--------------------------------------
C     Calculate lambda-X, lambda-Y, lambda-Z
C--------------------------------------
C
c     CALL CCLR_LAMTRA(XLAMDP,WORK(KLAMPX),XLAMDH,
c    *                 WORK(KLAMHX),X1AM,ISYMX)
c     CALL CCLR_LAMTRA(XLAMDP,WORK(KLAMPY),XLAMDH,
c    *                 WORK(KLAMHY),Y1AM,ISYMY)
      CALL GET_LAMBDAX(WORK(KLAMPX),WORK(KLAMHX),LISTX,IDLSTX,
     *                    ISYMX,XLAMDP,XLAMDH,WORK(KEND1),LWRK1)
      CALL GET_LAMBDAX(WORK(KLAMPY),WORK(KLAMHY),LISTY,IDLSTY,
     *                    ISYMY,XLAMDP,XLAMDH,WORK(KEND1),LWRK1)

      IF (Z1EXIST) THEN
         ISYMZ = ISYMZ1
      CALL GET_LAMBDAX(WORK(KLAMPZ),WORK(KLAMHZ),LISTZ,IDLSTZ,
     *                    ISYMZ,XLAMDP,XLAMDH,WORK(KEND1),LWRK1)
C        CALL CCLR_LAMTRA(XLAMDP,WORK(KLAMPZ),XLAMDH,
C    *                    WORK(KLAMHZ),Z1AM,ISYMZ)
      ELSE
         CALL DZERO(WORK(KLAMPZ),NLAMDT)
         CALL DZERO(WORK(KLAMHZ),NLAMDT)
         ISYMZ = 1
         CALL DCOPY(NLAMDT,XLAMDP,1,WORK(KLAMPZ),1)
         CALL DCOPY(NLAMDT,XLAMDH,1,WORK(KLAMHZ),1)
      END IF
C
C--------------------------------
C     Calculate integrals
C--------------------------------
C
      IF (DIRECT) THEN
         IF (HERDIR) THEN
            CALL HERDI1(WORK(KEND1),LWRK1,IPRERI)
         ELSE
            KCCFB1 = KEND1
            KINDXB = KCCFB1 + MXPRIM*MXCONT
            KEND1  = KINDXB + (8*MXSHEL*MXCONT + 1)/IRAT
            LWRK1  = LWORK  - KEND1
            CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     &                  KODPP1,KODPP2,KRDPP1,KRDPP2,
     &                  KFREE,LFREE,KEND1,WORK(KCCFB1),WORK(KINDXB),
     &                  WORK(KEND1),LWRK1,IPRERI)
            KEND1 = KFREE
            LWRK1 = LFREE
         ENDIF
         NTOSYM = 1
      ELSE
         NTOSYM = NSYM
      ENDIF
C
      DO ISYMD1 = 1,NTOSYM
C
         IF (DIRECT) THEN
            IF (HERDIR) THEN
               NTOT = MAXSHL
            ELSE
               NTOT = MXCALL
            ENDIF
         ELSE
            NTOT = NBAS(ISYMD1)
         ENDIF
C
         DO ILLL = 1,NTOT
C
C------------------------------------------
C        If direct calculate the integrals.
C------------------------------------------
C
            IF (DIRECT) THEN
C
               IF (HERDIR) THEN
                  CALL HERDI2(WORK(KEND1),LWRK1,INDEXA,ILLL,NUMDIS,
     &                        IPRERI)
               ELSE
                  CALL ERIDI2(ILLL,INDEXA,NUMDIS,0,0,
     &                        WORK(KODCL1),WORK(KODCL2),
     &                        WORK(KODBC1),WORK(KODBC2),
     &                        WORK(KRDBC1),WORK(KRDBC2),
     &                        WORK(KODPP1),WORK(KODPP2),
     &                        WORK(KRDPP1),WORK(KRDPP2),
     &                        WORK(KCCFB1),WORK(KINDXB),
     &                        WORK(KEND1), LWRK1,IPRERI)
               ENDIF
C
               KRECNR = KEND1
               KEND1  = KRECNR + (NBUFX(0) - 1)/IRAT + 1
               LWRK1  = LWORK  - KEND1
               IF (LWRK1 .LT. 0) THEN
                  CALL QUIT('Insufficient core in CC3_3BARINT')
               END IF
C
            ELSE
               NUMDIS = 1
            ENDIF
C
C--------------------------------------------------
C           Loop over number of distributions in disk.
C--------------------------------------------------
C
            DO IDEL2 = 1,NUMDIS
C
               IF (DIRECT) THEN
                  IDEL  = INDEXA(IDEL2)
                  IF (NOAUXB) THEN
                     IDUM = 1
                     CALL IJKAUX(IDEL,IDUM,IDUM,IDUM)
                  END IF
                  ISYMD = ISAO(IDEL)
               ELSE
                  IDEL  = IBAS(ISYMD1) + ILLL
                  ISYMD = ISYMD1
               ENDIF
C
               ISYDIS = MULD2H(ISYMD,ISYMOP)
               ISYAIK = MULD2H(ISYDIS,ISYMOP)
               ISYCKJ = ISYDIS
C
C------------------------------------------
C              Work space allocation no. 2.
C------------------------------------------
C
               KXINT  = KEND1
               KCKJ   = KXINT + NDISAO(ISYDIS)
               KEND2  = KCKJ  + NCKI(ISYCKJ)
               LWRK2  = LWORK - KEND2
C
               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
                  CALL QUIT('Insufficient space in CC3_GMAT ')
               ENDIF
C
               CALL DZERO(WORK(KCKJ),NCKI(ISYCKJ))
C
C-----------------------------------------
C              Read in batch of integrals.
C-----------------------------------------
C
               CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND2),LWRK2,
     *                     WORK(KRECNR),DIRECT)
C
C-----------------------------------------------------------
C              Calculate transformed X,Y,Z bar integrals 
C-----------------------------------------------------------
C
               CALL CC3_BARXYZINT(WORK(KXINT),XLAMDP,XLAMDH
     *                        ,WORK(KLAMPX),WORK(KLAMHX),ISYMX
     *                        ,WORK(KLAMPY),WORK(KLAMHY),ISYMY
     *                        ,WORK(KLAMPZ),WORK(KLAMHZ),ISYMZ,
     *                        WORK(KEND2),LWRK2,
     *                        IDEL,ISYMD,LU3SRTR,FN3SRTR,
     *                        LUCKJDR,FNCKJDR)
c     write(lupri,*)'IDEL,ISYMD',IDEL,ISYMD
c     write(lupri,*)'KXINT',ddot(NDISAO(ISYDIS),WORK(KXINT),
c    *                                      1,WORK(KXINT),1)
C
            ENDDO     ! IDEL2
         ENDDO        ! ILLL
      ENDDO           ! ISYMD1
C
      CALL QEXIT('CC3_3BARINT')
      RETURN
      END
C  /* Deck cc3_barxyzint */
      SUBROUTINE CC3_BARXYZINT(XINT,XLAMDP,XLAMDH,XLAMPX,XLAMHX,ISYMX,
     *                         XLAMPY,XLAMHY,ISYMY,XLAMPZ,XLAMHZ,ISYMZ,
     *                         WORK,LWORK,IDEL,ISYDEL,
     *                         LUOUT1,FNOUT1,LUOUT2,FNOUT2)
C
C     Purpose: Calculate X,Y,Z BAR integrals used in CC3.
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccinftap.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INTEGER ISYMX,ISYMY,ISYMZ,LWORK,IDEL,ISYDEL,LUOUT1,LUOUT2
      INTEGER KXCKD,KXCKJ,ID,LENGTH,IOFF
      INTEGER KEND1,LWRK1
      INTEGER ISYMXY,ISYMXYZ,ISYCKD,ISYCKJ
#if defined (SYS_CRAY)
      REAL XINT(*), XLAMDP(*),XLAMDH(*),XLAMPX(*),XLAMHX(*) 
      REAL XLAMPY(*),XLAMHY(*),XLAMPZ(*),XLAMHZ(*)
      REAL WORK(LWORK)
      real xnormval,ddot
#else
      DOUBLE PRECISION XINT(*), XLAMDP(*),XLAMDH(*),XLAMPX(*),XLAMHX(*) 
      DOUBLE PRECISION XLAMPY(*),XLAMHY(*),XLAMPZ(*),XLAMHZ(*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION xnormval,ddot
#endif

C
      CHARACTER*(*) FNOUT1, FNOUT2
C
      CALL QENTER('CC3_BARXYZINT')
C
c     xnormval = ddot(nlamdt,XLAMDP,1,XLAMDP,1)
c     write(lupri,*)'XLAMDP',xnormval
c     xnormval = ddot(nlamdt,XLAMDH,1,XLAMDH,1)
c     write(lupri,*)'XLAMDH',xnormval
c     xnormval = ddot(nglmdt(isymx),XLAMPX,1,XLAMPX,1)
c     write(lupri,*)'XLAMPX',xnormval,isymx
c     xnormval = ddot(nglmdt(isymx),XLAMHX,1,XLAMHX,1)
c     write(lupri,*)'XLAMHX',xnormval,isymx
c     xnormval = ddot(nglmdt(isymY),XLAMPY,1,XLAMPY,1)
c     write(lupri,*)'XLAMPY',xnormval,isymY
c     xnormval = ddot(nglmdt(isymY),XLAMHY,1,XLAMHY,1)
c     write(lupri,*)'XLAMHY',xnormval,isymY
c     xnormval = ddot(nglmdt(isymZ),XLAMPZ,1,XLAMPZ,1)
c     write(lupri,*)'XLAMPZ',xnormval,isymz
c     xnormval = ddot(nglmdt(isymz),XLAMHZ,1,XLAMHZ,1)
c     write(lupri,*)'XLAMHZ',xnormval,isymz
C
      ISYMXY  = MULD2H(ISYMX,ISYMY)
      ISYMXYZ = MULD2H(ISYMXY,ISYMZ)
      ISYCKD  = MULD2H(ISYMXYZ,MULD2H(ISYDEL,ISYMOP))
      ISYCKJ = ISYCKD
C
C---------------------------------
C        Allocation of work space.
C---------------------------------
C
      KXCKD  = 1
      KXCKJ  = KXCKD  + NCKATR(ISYCKD)
      KEND1  = KXCKJ  + NCKI(ISYCKJ)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient core in CC3_BARXYZINT')
      ENDIF
C
C----------------------------------------
C        Calculate transformed integrals.
C----------------------------------------
C
      CALL DZERO(WORK(KXCKD),NCKATR(ISYCKD))
      CALL DZERO(WORK(KXCKJ),NCKI(ISYCKJ))
C
      CALL CC3_CKD1(XINT,WORK(KXCKD),WORK(KXCKJ),
     *             XLAMPX(NT1AO(ISYMX)+1),XLAMHY,XLAMPZ(NT1AO(ISYMZ)+1),
     *             XLAMHZ,ISYMX,ISYMY,ISYMZ,ISYMZ,
     *             WORK(KEND1),LWRK1,IDEL,ISYDEL)
C
      CALL CC3_CKD1(XINT,WORK(KXCKD),WORK(KXCKJ),
     *             XLAMPX(NT1AO(ISYMX)+1),XLAMHZ,XLAMPY(NT1AO(ISYMY)+1),
     *             XLAMHY,ISYMX,ISYMZ,ISYMY,ISYMY,
     *             WORK(KEND1),LWRK1,IDEL,ISYDEL)
C
      CALL CC3_CKD1(XINT,WORK(KXCKD),WORK(KXCKJ),
     *             XLAMPY(NT1AO(ISYMY)+1),XLAMHX,XLAMPZ(NT1AO(ISYMZ)+1),
     *             XLAMHZ,ISYMY,ISYMX,ISYMZ,ISYMZ,
     *             WORK(KEND1),LWRK1,IDEL,ISYDEL)
C
      CALL CC3_CKD1(XINT,WORK(KXCKD),WORK(KXCKJ),
     *             XLAMPZ(NT1AO(ISYMZ)+1),XLAMHX,XLAMPY(NT1AO(ISYMY)+1),
     *             XLAMHY,ISYMZ,ISYMX,ISYMY,ISYMY,
     *             WORK(KEND1),LWRK1,IDEL,ISYDEL)
C
      CALL CC3_CKD1(XINT,WORK(KXCKD),WORK(KXCKJ),
     *             XLAMPY(NT1AO(ISYMY)+1),XLAMHZ,XLAMPX(NT1AO(ISYMX)+1),
     *             XLAMHX,ISYMY,ISYMZ,ISYMX,ISYMX,
     *             WORK(KEND1),LWRK1,IDEL,ISYDEL)
C
      CALL CC3_CKD1(XINT,WORK(KXCKD),WORK(KXCKJ),
     *             XLAMPZ(NT1AO(ISYMZ)+1),XLAMHY,XLAMPX(NT1AO(ISYMX)+1),
     *             XLAMHX,ISYMZ,ISYMY,ISYMX,ISYMX,
     *             WORK(KEND1),LWRK1,IDEL,ISYDEL)
C
C--------------------------------
C     Write to disk (ck|d alpha).
C--------------------------------
C
      ID     = IDEL - IBAS(ISYDEL)
C
      LENGTH = NCKATR(ISYCKD)
C
      IOFF = ICKDAO(ISYCKD,ISYDEL) + NCKATR(ISYCKD)*(ID - 1) + 1
C
      IF (LENGTH .GT. 0) THEN
         CALL PUTWA2(LUOUT1,FNOUT1,WORK(KXCKD),IOFF,LENGTH)
      ENDIF
c     write(lupri,*)'id,isydel',id,isydel
C     xnormval = ddot(NCKATR(ISYCKD),WORK(KXCKD),1,WORK(KXCKD),1)
C     write(lupri,*)'KXCKD',xnormval
C
      LENGTH = NCKI(ISYCKJ)
C
      IOFF  = ICKID(ISYCKJ,ISYDEL) + NCKI(ISYCKJ)*(ID - 1) + 1
C
      IF (LENGTH .GT. 0) THEN
         CALL PUTWA2(LUOUT2,FNOUT2,WORK(KXCKJ),IOFF,LENGTH)
      ENDIF
C      xnormval = ddot(NCKI(ISYCKJ),WORK(KXCKJ),1,WORK(KXCKJ),1)
C      write(lupri,*)'KXCKJ',xnormval
C
      CALL QEXIT('CC3_BARXYZINT')
      RETURN
      END
