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
C  /* Deck cc3_lhtr_l3 */
      SUBROUTINE CC3_FMAT(LISTL,IDLSTL,LISTB,IDLSTB,IOPTRES,
     *                    OMEGA1,OMEGA2,OMEGA1EFF,OMEGA2EFF,
     *                    IDOTS,DOTPROD,LISTDP,ITRAN,
     *                    NFTRAN,MXVEC,WORK,LWORK,
     *                    LUCKJD,FNCKJD,LUDKBC,FNDKBC,LUTOC,FNTOC,
     *                    LU3VI,FN3VI,LUDKBC3,FNDKBC3,
     *                    LU3FOPX,FN3FOPX,LU3FOP2X,FN3FOP2X)
C
C     Written by K. Hald, Spring 2002.
C
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
#include <second.h>
C
      INTEGER IDLSTL, IDLSTB, IOPTRES, ITRAN, NFTRAN, MXVEC, LWORK
      INTEGER LUCKJD, LUDKBC, LUTOC, LU3VI, LUDKBC3, LU3FOPX, LU3FOP2X
      INTEGER IDOTS(MXVEC,NFTRAN)
      INTEGER LUDKBC4, LU4V, ISYMTR, ISINT1, ISINT2, ISYMIM
      INTEGER ISYMT1, ISYMT2, ISYML, ISYML1, ISYML2, KT1AM
      INTEGER KL1AM, KL2TP, KLAMDP, KLAMDH, IOPT
      INTEGER KFOCKD, KFCKBA, KEND0, LWRK0, KEND1, LWRK1
      INTEGER KEND2, LWRK2, KEND3, LWRK3, KEND4, LWRK4
      INTEGER LUFCK, ISYMC, ISYMK, KOFF1, KOFF2, KTROC, KTROC0, KTROC1
      INTEGER KTROC2, KXIAJB, KINTOC, LENGTH, ISYOPE, IOPTTCME
      INTEGER IOFF, ISYMD, ISAIJ1, ISYCKB, ISCKB1, ISCKB2, KTRVI
      INTEGER KTRVI1, KTRVI2, KTRVI0, KTRVI3, KTRVI4
      INTEGER KTRVI5, KTRVI6, KINTVI, ISYMB, ISYALJ, ISAIJ2, ISYMBD
      INTEGER ISCKIJ, KSMAT, KQMAT, KDIAG, KINDSQ, KINDEX, KTMAT
      INTEGER KRMAT1, KRMAT2, LENSQ, ISYRES
      INTEGER ISYMB1, ISYMB2, KB1AM, KB2TP, KT2TP
      INTEGER KVVVVO, KOIOOOO, KOIOVVO, KOIOOVV
      INTEGER KVVVVB, KBIOOOO, KBIOVVO, KBIOOVV, LU4VB
      INTEGER LU3SRTR, LUCKJDR, LUDELDR, LUDKBCR, LUDKBCR4
      INTEGER ISINTR1
      INTEGER KTROC3, KTROC4, ISCKBR, KTRVI7, KTRVI8
c
      integer kt3am
c
C     Functions : 
      INTEGER ILSTSYM
C
#if defined (SYS_CRAY)
      REAL OMEGA1(*), OMEGA2(*), OMEGA1EFF(*), OMEGA2EFF(*)
      REAL DOTPROD(MXVEC,NFTRAN), WORK(LWORK)
      REAL TITRAN, TISORT, TISMAT, TIQMAT, TICONT, TIOME1
      REAL HALF, ONE, TWO
      REAL DTIME, DDOT, XNORM
#else
      DOUBLE PRECISION OMEGA1(*), OMEGA2(*), OMEGA1EFF(*), OMEGA2EFF(*)
      DOUBLE PRECISION DOTPROD(MXVEC,NFTRAN), WORK(LWORK)
      DOUBLE PRECISION TITRAN, TISORT, TISMAT, TIQMAT, TICONT, TIOME1
      DOUBLE PRECISION HALF, ONE, TWO
      DOUBLE PRECISION DTIME, DDOT, XNORM
#endif
C
      CHARACTER*(*) FNCKJD, FNDKBC, FNTOC, FN3VI, FNDKBC3, FN3FOPX
      CHARACTER*(*) FN3FOP2X
      CHARACTER*1 CDUMMY
      CHARACTER*3 LISTL, LISTB, LISTDP
      CHARACTER*10 MODEL
      CHARACTER*13 FNDKBC4, FN4V, FN4VB, FN3SRTR, FNCKJDR
      CHARACTER*13 FNDELDR, FNDKBCR, FNDKBCR4
C
      PARAMETER(HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      CALL QENTER('CC3_FMAT')
C
C----------------------------------------------------
C     Initialise character strings and open files
C----------------------------------------------------
C
      CDUMMY = ' '
C
      LUDKBC4  = -1
      LU4V     = -1
      LU4VB    = -1
      LU3SRTR  = -1
      LUCKJDR  = -1
      LUDELDR  = -1
      LUDKBCR  = -1
      LUDKBCR4 = -1
C
      FNDKBC4  = 'CC3_FMAT_TMP1'
      FN4V     = 'CC3_FMAT_TMP2'
      FN4VB    = 'CC3_FMAT_TMP3'
      FN3SRTR  = 'CC3_FMAT_TMP4'
      FNCKJDR  = 'CC3_FMAT_TMP5'
      FNDELDR  = 'CC3_FMAT_TMP6'
      FNDKBCR  = 'CC3_FMAT_TMP7'
      FNDKBCR4 = 'CC3_FMAT_TMP8'
C
      CALL WOPEN2(LUDKBC4,FNDKBC4,64,0)
      CALL WOPEN2(LU4V,FN4V,64,0)
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
C     omega = int1*T2*int2
C     isymres is symmetry of result(omega)
C     isint1 is symmetry of integrals in contraction.(int1)
C     isint2 is symmetry of integrals in the triples equation.(int2)
C     isymim is symmetry of S and Q intermediates.(t2*int2)
C      (sym is for all index of S and Q (cbd,klj)
C       thus cklj=b*d*isymim)
C-------------------------------------------------------------
C
      IPRCC   = IPRINT
      ISYMTR  = ISYMOP
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
      ISINTR1 = MULD2H(ISINT1,ISYMB1)
      ISYMIM  = MULD2H(ISYML,ISYMOP)
      ISYRES  = MULD2H(ISYMIM,ISYMB)
C
C--------------------
C     Time variables.
C--------------------
C
      TITRAN = 0.0D0
      TISORT = 0.0D0
      TISMAT = 0.0D0
      TIQMAT = 0.0D0
      TICONT = 0.0D0
      TIOME1 = 0.0D0
C
C-----------------------------------
C     Work space allocation
C-----------------------------------
C
      KFOCKD  = 1
      KFCKBA  = KFOCKD + NORBTS
      KT1AM   = KFCKBA + N2BST(ISYMOP)
      KL1AM   = KT1AM  + NT1AM(ISYMT1)
      KL2TP   = KL1AM  + NT1AM(ISYML1)
      KB1AM   = KL2TP  + NT2SQ(ISYML2)
      KB2TP   = KB1AM  + NT1AM(ISYMB1)
      KT2TP   = KB2TP  + NT2SQ(ISYMB2)
      KLAMDP  = KT2TP  + NT2SQ(ISYMT2)
      KLAMDH  = KLAMDP + NLAMDT
      KEND0   = KLAMDH + NLAMDT
      LWRK0   = LWORK  - KEND0
C
      if (.false.) then
         kt3am  = kend0
         kend0  = kt3am + nt1amx*nt1amx*nt1amx
         lwrk0 = lwork - kend0
         call dzero(work(kt3am),nt1amx*nt1amx*nt1amx)
      endif
C
      IF (LWRK0 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND0
         CALL QUIT('Insufficient space in CC3_FMAT')
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
C-----------------------------------------------------------
C     Resort DKBC -> DKBC4
C-----------------------------------------------------------
C
      CALL CC3_TCME(WORK(KLAMDP),ISINT1,WORK(KEND0),LWRK0,IDUMMY,CDUMMY,
     *              LUDKBC,FNDKBC,IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *              IDUMMY,CDUMMY,LUDKBC4,FNDKBC4,2)
C
C--------------------------------------
C     Reorder the t2-amplitudes i T2TP.
C--------------------------------------
C
      IF (LWORK .LT. NT2SQ(1)) THEN
        CALL QUIT('Not enough memory to construct T2TP in CC3_FMAT')
      ENDIF
C
      IOPT = 2
      CALL CC_RDRSP('R0',0,1,IOPT,MODEL,
     *              DUMMY,WORK(KT2TP))
      CALL CC_T2SQ(WORK(KT2TP),WORK(KEND0),ISYMT2)
      CALL CC3_T2TP(WORK(KT2TP),WORK(KEND0),ISYMT2)
C
      IF (IPRINT .GT. 55) THEN
         XNORM = DDOT(NT2SQ(ISYMT2),WORK(KT2TP),1,WORK(KT2TP),1)
         WRITE(LUPRI,*) 'Norm of T2TP ',XNORM
      ENDIF
C
C--------------------------------------
C     Reorder the l2-amplitudes i L2TP.
C--------------------------------------
C
      IF (LWORK .LT. NT2SQ(ISYML2)) THEN
        CALL QUIT('Not enough memory to construct L2TP in CC3_FMAT')
      ENDIF
C
      IOPT = 3
      CALL CC_RDRSP(LISTL,IDLSTL,ISYML2,IOPT,MODEL,
     *              WORK(KL1AM),WORK(KL2TP))
      CALL CC_T2SQ(WORK(KL2TP),WORK(KEND0),ISYML2)
      CALL CC3_T2TP(WORK(KL2TP),WORK(KEND0),ISYML2)
C
      IF (IPRINT .GT. 55) THEN
         XNORM = DDOT(NT2SQ(ISYML2),WORK(KL2TP),1,WORK(KL2TP),1)
         WRITE(LUPRI,*) 'Norm of L2TP ',XNORM
      ENDIF
C
C--------------------------------------
C     Reorder the B2-amplitudes i B2TP.
C--------------------------------------
C
      IF (LWORK .LT. NT2SQ(ISYMB2)) THEN
        CALL QUIT('Not enough memory to construct B2TP in CC3_FMAT')
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
         XNORM = DDOT(NT2SQ(ISYMB2),WORK(KB2TP),1,WORK(KB2TP),1)
         WRITE(LUPRI,*) 'Norm of B2TP ',XNORM
      ENDIF
C
C------------------------------------------------------
C     Calculate the (ck|de)-tilde and (ck|lm)-tilde
C------------------------------------------------------
C
      CALL CC3_BARINT(WORK(KB1AM),ISYMB1,WORK(KLAMDP),
     *                WORK(KLAMDH),WORK(KEND0),LWRK0,
     *                LU3SRTR,FN3SRTR,LUCKJDR,FNCKJDR)
C
      CALL CC3_SORT1(WORK(KEND0),LWRK0,2,ISINTR1,LU3SRTR,FN3SRTR,
     *               LUDELDR,FNDELDR,IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *               IDUMMY,CDUMMY)
C
      CALL CC3_SINT(WORK(KLAMDH),WORK(KEND0),LWRK0,ISINTR1,
     *              LUDELDR,FNDELDR,LUDKBCR,FNDKBCR)
C
      CALL CC3_TCME(WORK(KLAMDP),ISINTR1,WORK(KEND0),LWRK0,
     *              IDUMMY,CDUMMY,LUDKBCR,FNDKBCR,
     *              IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *              IDUMMY,CDUMMY,LUDKBCR4,FNDKBCR4,2)
C
C--------------------------------------------------------------
C     Calculate the normal g^0 integrals for
C     OOOO, OVVO and OOVV integrals. VVVV is stored on file
C--------------------------------------------------------------
C
      KOIOOOO = KEND0
      KOIOVVO = KOIOOOO + N3ORHF(ISYMOP)
      KOIOOVV = KOIOVVO + NT2SQ(ISYMOP)
      KEND0   = KOIOOVV + NT2SQ(ISYMOP)
      LWRK0   = LWORK   - KEND0
C
      IF (LWRK0 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_FMAT (g^0[2o2v] kind)')
      ENDIF
C
      CALL DZERO(WORK(KOIOVVO),NT2SQ(ISYMOP))
      CALL DZERO(WORK(KOIOOVV),NT2SQ(ISYMOP))
C
      CALL CC3_INT(WORK(KOIOOOO),WORK(KOIOOVV),WORK(KOIOVVO),
     *             ISYMOP,LU4V,FN4V,
     *             .FALSE.,'DUMMY',IDUMMY,IDUMMY,
     *             .FALSE.,'DUMMY',IDUMMY,IDUMMY,
     *             WORK(KEND0),LWRK0)
C
C---------------------------------------------------------------
C     Calculate the g^B integrals for
C     OOOO, OVVO and OOVV integrals. VVVV is stored on file
C---------------------------------------------------------------
C
      KBIOOOO = KEND0
      KBIOVVO = KBIOOOO + N3ORHF(ISINTR1)
      KBIOOVV = KBIOVVO + NT2SQ(ISINTR1)
      KEND0   = KBIOOVV + NT2SQ(ISINTR1)
      LWRK0   = LWORK   - KEND0
C
      IF (LWRK0 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_FMAT (g^B[2o2v] kind)')
      ENDIF
C
      CALL DZERO(WORK(KBIOOOO),N3ORHF(ISINTR1))
      CALL DZERO(WORK(KBIOVVO),NT2SQ(ISINTR1))
      CALL DZERO(WORK(KBIOOVV),NT2SQ(ISINTR1))
C
      CALL CC3_INT(WORK(KBIOOOO),WORK(KBIOOVV),WORK(KBIOVVO),
     *             ISINTR1,LU4VB,FN4VB,
     *             .TRUE.,LISTB,IDLSTB,ISYMB1,
     *             .FALSE.,'DUMMY',IDUMMY,IDUMMY,
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
         CALL AROUND( 'In CC3_FMAT: Triples Fock MO matrix' )
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
         CALL AROUND('In CC3_FMAT: Triples Fock MO matrix (sort)')
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
      KTROC4 = KTROC3 + NTRAOC(ISINTR1)
      KXIAJB = KTROC4 + NTRAOC(ISINTR1)
      KEND1  = KXIAJB + NT2AM(ISYMOP)
      LWRK1  = LWORK  - KEND1
C
      KINTOC = KEND1
      KEND2  = KINTOC + MAX(NTOTOC(ISYMOP),NTOTOC(ISINT2))
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
         XNORM = DDOT(NT2AM(ISYMOP),WORK(KXIAJB),1,
     *                WORK(KXIAJB),1)
         WRITE(LUPRI,*) 'Norm of IAJB ',XNORM
      ENDIF
C
C--------------------------------------------------------------------------
C     Read in B-transformed integrals used in contractions and transform.
C--------------------------------------------------------------------------
C
      IOFF = 1
      IF (NTOTOC(ISINTR1) .GT. 0) THEN
         CALL GETWA2(LUCKJDR,FNCKJDR,WORK(KINTOC),IOFF,NTOTOC(ISINTR1))
      ENDIF
C
      CALL CC3_TROCC(WORK(KINTOC),WORK(KTROC3),WORK(KLAMDP),
     *                    WORK(KEND2),LWRK2,ISINTR1)
C
      CALL CCFOP_SORT(WORK(KTROC3),WORK(KTROC4),ISINTR1,1)
C
      CALL CC3_LSORT1(WORK(KTROC3),ISINTR1,WORK(KEND2),LWRK2,5)
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
         XNORM = DDOT(NTOTOC(ISINT2),WORK(KINTOC),1,
     *                WORK(KINTOC),1)
         WRITE(LUPRI,*) 'Norm of CKJDEL-INT  ',XNORM
      ENDIF
C
      IF (IPRINT .GT. 55) THEN
         XNORM = DDOT(NTRAOC(ISINT2),WORK(KTROC0),1,
     *                WORK(KTROC0),1)
         WRITE(LUPRI,*) 'Norm of TROC0 CC3_FMAT : ',XNORM
      ENDIF
C
      IF (IPRINT .GT. 55) THEN
         XNORM = DDOT(NTRAOC(ISINT2),WORK(KTROC2),1,
     *                WORK(KTROC2),1)
         WRITE(LUPRI,*) 'Norm of TROC2 CC3_FMAT : ',XNORM
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
         ISCKBR = MULD2H(ISINTR1,ISYMD)
C
         IF (IPRINT .GT. 55) THEN
C
            WRITE(LUPRI,*) 'In CC3_FMAT: ISAIJ1 :',ISAIJ1
            WRITE(LUPRI,*) 'In CC3_FMAT: ISYCKB :',ISYCKB
            WRITE(LUPRI,*) 'In CC3_FMAT: ISCKB2 :',ISCKB2
            WRITE(LUPRI,*) 'In CC3_FMAT: ISCKBR :',ISCKBR
C
         ENDIF
C
C--------------------------
C        Memory allocation.
C--------------------------
C
         KTRVI  = KEND1
         KTRVI1 = KTRVI  + NCKATR(ISCKB1)
         KTRVI7 = KTRVI1 + NCKATR(ISCKB1)
         KTRVI8 = KTRVI7 + NCKATR(ISCKBR)
         KTRVI2 = KTRVI8 + NCKATR(ISCKBR)
         KRMAT1 = KTRVI2 + NCKATR(ISCKB2)
         KEND2  = KRMAT1 + NCKI(ISAIJ1)
         LWRK2  = LWORK  - KEND2
C
         KTRVI0  = KEND2
         KTRVI3  = KTRVI0  + NCKATR(ISCKB2)
         KTRVI4  = KTRVI3  + NCKATR(ISCKB2)
         KTRVI5  = KTRVI4  + NCKATR(ISCKB2)
         KTRVI6  = KTRVI5  + NCKATR(ISCKB2)
         KVVVVO  = KTRVI6  + NCKATR(ISCKB2)
         KVVVVB  = KVVVVO  + NMAABC(ISCKB2)
         KEND3   = KVVVVB  + NMAABC(ISCKBR)
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
C           Read in g^0_{vvvv} for a given D
C--------------------------------------------
C
            IF (NMAABC(ISCKB2) .GT. 0) THEN
               IOFF = I3VVIR(ISCKB2,ISYMD)
     *              + NMAABC(ISCKB2)*(D-1)
     *              + 1
               CALL GETWA2(LU4V,FN4V,WORK(KVVVVO),IOFF,NMAABC(ISCKB2))
            ENDIF
C
C--------------------------------------------
C           Read in g^B_{vvvv} for a given D
C--------------------------------------------
C
            IF (NMAABC(ISCKBR) .GT. 0) THEN
               IOFF = I3VVIR(ISCKBR,ISYMD)
     *              + NMAABC(ISCKBR)*(D-1)
     *              + 1
               CALL GETWA2(LU4VB,FN4VB,WORK(KVVVVB),IOFF,NMAABC(ISCKBR))
            ENDIF
C
C--------------------------------------------------------------
C           Read B-transformed integrals used in contraction
C--------------------------------------------------------------
C
            IF (NCKATR(ISCKBR) .GT. 0) THEN
               IOFF = ICKBD(ISCKBR,ISYMD) + NCKATR(ISCKBR)*(D - 1) + 1
               CALL GETWA2(LUDKBCR,FNDKBCR,WORK(KTRVI7),IOFF,
     &                     NCKATR(ISCKBR))
            ENDIF
C
            IF (LWRK4 .LT. NCKATR(ISCKBR)) THEN
               CALL QUIT('Insufficient space for allocation in '//
     &                   'CC3_FMAT (TRVI7)')
            END IF
C
            DTIME = SECOND()
            CALL CCSDT_SRVIR3(WORK(KTRVI7),WORK(KEND4),ISYMD,D,ISINTR1)
C
            IF (NCKATR(ISCKBR) .GT. 0) THEN
               IOFF = ICKBD(ISCKBR,ISYMD) + NCKATR(ISCKBR)*(D - 1) + 1
               CALL GETWA2(LUDKBCR4,FNDKBCR4,WORK(KTRVI8),IOFF,
     &                     NCKATR(ISCKBR))
            ENDIF
C
            IF (LWRK4 .LT. NCKATR(ISCKBR)) THEN
               CALL QUIT('Insufficient space for allocation in '//
     &                   'CC3_FMAT (TRVI8)')
            END IF
C
            DTIME = SECOND()
            CALL CCSDT_SRVIR3(WORK(KTRVI8),WORK(KEND4),ISYMD,D,ISINTR1)
C
C------------------------------------------------------------
C           Read and transform integrals used in contraction.
C------------------------------------------------------------
C
            IF (NCKATR(ISCKB1) .GT. 0) THEN
               IOFF = ICKBD(ISCKB1,ISYMD) + NCKATR(ISCKB1)*(D - 1) + 1
               CALL GETWA2(LUDKBC,FNDKBC,WORK(KTRVI),IOFF,
     &                     NCKATR(ISCKB1))
            ENDIF
C
            IF (LWRK4 .LT. NCKATR(ISCKB1)) THEN
               CALL QUIT('Insufficient space for allocation in '//
     &                   'CC3_L3 (TRVI)')
            END IF
C
            DTIME = SECOND()
            CALL CCSDT_SRVIR3(WORK(KTRVI),WORK(KEND4),ISYMD,D,ISINT1)
C
            DTIME  = SECOND() - DTIME
            TISORT = TISORT   + DTIME
C
            IF (NCKATR(ISCKB1) .GT. 0) THEN
               IOFF = ICKBD(ISCKB1,ISYMD) + NCKATR(ISCKB1)*(D - 1) + 1
               CALL GETWA2(LUDKBC4,FNDKBC4,WORK(KTRVI1),IOFF,
     &                     NCKATR(ISCKB1))
            ENDIF
C
            IF (LWRK4 .LT. NCKATR(ISCKB1)) THEN
               CALL QUIT('Insufficient space for allocation in '//
     &                   'CC3_L3 (TRVI1)')
            END IF
C
            DTIME = SECOND()
            CALL CCSDT_SRVIR3(WORK(KTRVI1),WORK(KEND4),ISYMD,D,ISINT1)
C
            DTIME  = SECOND() - DTIME
            TISORT = TISORT   + DTIME
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
            DTIME = SECOND()
            CALL CCSDT_SRTVIR(WORK(KTRVI0),WORK(KTRVI2),WORK(KEND4),
     *                        LWRK4,ISYMD,ISINT2)
C
            CALL CCSDT_SRTVIR(WORK(KTRVI4),WORK(KTRVI5),WORK(KEND4),
     *                        LWRK4,ISYMD,ISINT2)
C
            DTIME  = SECOND() - DTIME
            TISORT = TISORT   + DTIME
C
            IF (IPRINT .GT. 55) THEN
               XNORM = DDOT(NCKATR(ISCKB2),WORK(KTRVI0),1,
     *                      WORK(KTRVI0),1)
               WRITE(LUPRI,*) 'Norm of TRVI0 ',XNORM
            ENDIF
C
            IF (IPRINT .GT. 55) THEN
               XNORM = DDOT(NCKATR(ISCKB2),WORK(KTRVI2),1,
     *                      WORK(KTRVI2),1)
               WRITE(LUPRI,*) 'Norm of TRVI2 ',XNORM
            ENDIF
C
            IF (IPRINT .GT. 55) THEN
               XNORM = DDOT(NCKATR(ISCKB2),WORK(KTRVI4),1,
     *                      WORK(KTRVI4),1)
               WRITE(LUPRI,*) 'Norm of TRVI4 ',XNORM
            ENDIF
C
            IF (IPRINT .GT. 55) THEN
               XNORM = DDOT(NCKATR(ISCKB2),WORK(KTRVI5),1,
     *                      WORK(KTRVI5),1)
               WRITE(LUPRI,*) 'Norm of TRVI5 ',XNORM
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
     *                   'CC3_FMAT  (TRVI3)')
            END IF
C
C           Can use kend3 since dont need the integrals anymore
            DTIME = SECOND()
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
     &                   'CC3_FMAT (TRVI6)')
            END IF
C
C           Can use kend3 since dont need the integrals anymore
            DTIME = SECOND()
            CALL CCSDT_SRVIR3(WORK(KTRVI6),WORK(KEND4),ISYMD,D,ISINT2)
C
            IF (IPRINT .GT. 55) THEN
               XNORM = DDOT(NCKATR(ISCKB2),WORK(KTRVI6),1,
     *                      WORK(KTRVI6),1)
               WRITE(LUPRI,*) 'Norm of TRVI6 ',XNORM
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
                  XNORM = DDOT(NCKIJ(ISCKIJ),WORK(KDIAG),1,
     *                    WORK(KDIAG),1)
                  WRITE(LUPRI,*) 'Norm of DIA  ',XNORM
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
                  DTIME = SECOND()
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
                  DTIME  = SECOND() - DTIME
                  TISMAT = TISMAT   + DTIME
C
                  IF (IPRINT .GT. 55) THEN
                     XNORM = DDOT(NCKIJ(ISCKIJ),WORK(KSMAT),1,
     *                       WORK(KSMAT),1)
                     WRITE(LUPRI,*) 'Norm of SMAT-L3    ',XNORM
                  ENDIF
C
C-------------------------------------------------------------------
C                 Calculate Q(ci,jk) for fixed b,d for t3-bar.
C-------------------------------------------------------------------
C
                  DTIME = SECOND()
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
                  DTIME  = SECOND() - DTIME
                  TIQMAT = TIQMAT   + DTIME
C
                  IF (IPRINT .GT. 55) THEN
                     XNORM = DDOT(NCKIJ(ISCKIJ),WORK(KQMAT),1,
     *                       WORK(KQMAT),1)
                     WRITE(LUPRI,*) 'Norm of QMAT-L3  ',XNORM
                  ENDIF
C
C-----------------------------------------------------------------------
C                 Calculate the contributions to omega2
C-----------------------------------------------------------------------
C
                  CALL CCFOP_CONVIR(WORK(KRMAT2),WORK(KSMAT),
     *                              WORK(KQMAT),WORK(KTMAT),ISYMIM,
     *                              WORK(KTRVI7),WORK(KTRVI8),ISINTR1,
     *                              WORK(KEND4),LWRK4,WORK(KINDSQ),
     *                              LENSQ,ISYMB,B,ISYMD,D)
C
                  CALL CCFOP_CONOCC(OMEGA2,WORK(KRMAT1),
     *                              WORK(KRMAT2),WORK(KSMAT),
     *                              WORK(KTMAT),ISYMIM,
     *                              WORK(KTROC3),WORK(KTROC4),ISINTR1,
     *                              WORK(KEND4),LWRK4,WORK(KINDSQ),
     *                              LENSQ,ISYMB,B,ISYMD,D)
C
C------------------------------------------------------------------------
C                 Calculate the L3 contribution to omega1
C------------------------------------------------------------------------
C
                  CALL DSCAL(NCKIJ(ISCKIJ),-ONE,WORK(KSMAT),1)
                  CALL DSCAL(NCKIJ(ISCKIJ),-ONE,WORK(KQMAT),1)
C
C     t^{B}_{1} <L3|[H^0,T^{C}_{2}],tau_{1}]|HF>
                  CALL CC3_L3_OMEGA1(OMEGA1,ISYRES,WORK(KSMAT),
     *                               WORK(KQMAT),WORK(KTMAT),ISYMIM,
     *                               WORK(KOIOOOO),WORK(KOIOVVO),
     *                               WORK(KOIOOVV),WORK(KVVVVO),ISINT1,
     *                               WORK(KB2TP),ISYMB2,WORK(KEND4),
     *                               LWRK4,LENSQ,WORK(KINDSQ),
     *                               ISYMB,B,ISYMD,D)
C
                  IF (IPRINT .GT. 55) THEN
                     XNORM = DDOT(NT1AM(ISYRES),OMEGA1,1,OMEGA1,1)
                     WRITE(LUPRI,*) 'Norm CC3_FMAT after CC3_L3_OMEGA1',
     *                               XNORM
                  ENDIF
C
C     t^{B}_{1} <L3|[H^B,T^{0}_{2}],tau_{1}]|HF>
                  CALL CC3_L3_OMEGA1(OMEGA1,ISYRES,WORK(KSMAT),
     *                               WORK(KQMAT),WORK(KTMAT),ISYMIM,
     *                               WORK(KBIOOOO),WORK(KBIOVVO),
     *                               WORK(KBIOOVV),WORK(KVVVVB),ISINTR1,
     *                               WORK(KT2TP),ISYMT2,WORK(KEND4),
     *                               LWRK4,LENSQ,WORK(KINDSQ),
     *                               ISYMB,B,ISYMD,D)
C
                  IF (IPRINT .GT. 55) THEN
                     XNORM = DDOT(NT1AM(ISYRES),OMEGA1,1,OMEGA1,1)
                     WRITE(LUPRI,*) 'Norm CC3_FMAT after CC3_L3_OMEGA1',
     *                               XNORM
                  ENDIF
C
C-------------------------------------------
C                 Accumulate R2 in Omega2
C-------------------------------------------
C
                  CALL CC3_RACC(OMEGA2,WORK(KRMAT2),ISYMB,B,
     *                          ISYRES)
C
                  IF (IPRINT .GT. 55) THEN
                     XNORM = DDOT(NT2AM(ISYRES),OMEGA2,1,OMEGA2,1)
                     WRITE(LUPRI,*) 'Norm of Rho22-after CC3_RACC',XNORM
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
            CALL CC3_RACC(OMEGA2,WORK(KRMAT1),ISYMD,D,ISYRES)
C
            IF (IPRINT .GT. 55) THEN
               XNORM = DDOT(NT2AM(ISYRES),OMEGA2,1,OMEGA2,1)
               WRITE(LUPRI,*) 'Norm of Rho22-after CC3_RACC-2',XNORM
            ENDIF
C
            IF (IPRINT .GT. 220) THEN
               CALL AROUND('After CC3_RACC-2: ')
               CALL CC_PRP(DUMMY,OMEGA2,ISYRES,0,1)
            ENDIF
C
         ENDDO       ! D
      ENDDO          ! ISYMD
C
COMMENT COMMENT COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT COMMENT COMMENT
c      write(lupri,*) 'L3 Amplitudes in cc3_fmat : '
C      call dscal(nt1amx*nt1amx*nt1amx,-1.0D0,work(kt3am),1)
c      call print_pt3(work(kt3am),isckij,3)
COMMENT COMMENT COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT COMMENT COMMENT
C-------------------------------
C     Close and delete files
C-------------------------------
C
      CALL WCLOSE2(LUDKBC4,FNDKBC4,'DELETE')
      CALL WCLOSE2(LU4V,FN4V,'DELETE')
      CALL WCLOSE2(LU4VB,FN4VB,'DELETE')
      CALL WCLOSE2(LU3SRTR,FN3SRTR,'DELETE')
      CALL WCLOSE2(LUCKJDR,FNCKJDR,'DELETE')
      CALL WCLOSE2(LUDELDR,FNDELDR,'DELETE')
      CALL WCLOSE2(LUDKBCR,FNDKBCR,'DELETE')
      CALL WCLOSE2(LUDKBCR4,FNDKBCR4,'DELETE')
C
COMMENT COMMENT COMMENT COMMENT 
C       write(lupri,*) 'OMEGA1EFF in cc3_fmat before added'
C       do i = 1,nt1am(isyres)
C           write(lupri,*) i,OMEGA1EFF(i)
C       end do
COMMENT COMMENT COMMENT COMMENT 
C     write(lupri,*) 'OMEGA2EFF in cc3_fmat before added'
C     call outpak(OMEGA2EFF,NT2AM(isyres),1,lupri)
COMMENT COMMENT COMMENT COMMENT 

C     DO I = 1,NT1AM(ISYRES)
C        OMEGA1EFF(I) = OMEGA1EFF(I) + OMEGA1(I)
C     END DO
C
C     DO I = 1,NT2AM(ISYRES)
C        OMEGA2EFF(I) = OMEGA2EFF(I) + OMEGA2(I)
C     END DO
COMMENT COMMENT COMMENT COMMENT 
C       write(lupri,*) 'OMEGA1EFF in cc3_fmat after added'
C       do i = 1,nt1am(isyres)
C           write(lupri,*) i,OMEGA1EFF(i)
C       end do
COMMENT COMMENT COMMENT COMMENT 
C     write(lupri,*) 'OMEGA2EFF in cc3_fmat after added'
C     call outpak(OMEGA2EFF,NT2AM(isyres),1,lupri)
COMMENT COMMENT COMMENT COMMENT 

C
C
C-------------------
C     Print timings.
C-------------------
C
      IF (IPRINT .GT. 9) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*)
         WRITE(LUPRI,1) 'CC3_TRAN  : ',TITRAN
         WRITE(LUPRI,1) 'CC3_SORT  : ',TISORT
         WRITE(LUPRI,1) 'CC3_SMAT  : ',TISMAT
         WRITE(LUPRI,1) 'CC3_QMAT  : ',TIQMAT
         WRITE(LUPRI,1) 'CC3_OME1  : ',TIOME1
         WRITE(LUPRI,*)
      END IF
C
C-------------
C     End
C-------------
C
      CALL QEXIT('CC3_FMAT')
C
      RETURN
C
    1 FORMAT(7X,'Time used in',2X,A12,F12.2,' seconds')
C
      END
C  /* Deck cc3_int */
      SUBROUTINE CC3_INT(XOOOO,XOOVV,XOVVO,ISYINT,LU4V,FN4V,
     *                   DO_LAMB,LISTB,IDLSTB,ISYMB,
     *                   DO_LAMC,LISTC,IDLSTC,ISYMC,
     *                   WORK,LWORK)
C
C     Written by Kasper Hald, Summer 2002
C
C     Purpose:
C
C     Calculate g_{oooo}, g_{ovvo} and g_{oovv}
C
C     IF DO_LAMB = TRUE we calculate the T^{B}_{1} transformed integrals.
C     IF DO_LAMB = TRUE and DO_LAMC = TRUE calculated T^{B}_{1} T^{C}_{1}
C     (doubly) transformed integrals.
C     Otherwise we calculate the ordinary (lambda^{(0)}) integrals.
C
C     DO_LAMC can only be TRUE if DO_LAMB is TRUE !!!
C
C     IF LVVVV = .TRUE. calculate VVVV integrals and store them on LU4V file.
C     (LVVVV is sitting in a common block in the ccsdinp.h file and defined 
C     through the input).)
C
      IMPLICIT NONE
C
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
#include <ccsdsym.h>
#include <distcl.h>
#include <cbieri.h>
#include <eribuf.h>
#include <second.h>
C
      INTEGER ISYINT, LU4V, IDLSTB, ISYMB, LWORK
      INTEGER KLAMDP, KLAMDH, KEND1, LWRK1, KEND2, LWRK2, KENDS2, LWRKS2
      INTEGER KENDSV, LWRKSV, KEND3, LWRK3
      INTEGER LUIN3O, LUIN3O2, LU3V, LU3VB
      INTEGER KCCFB1, KINDXB, KODCL1, KODCL2, KODBC1, KODBC2
      INTEGER KRDBC1, KRDBC2, KODPP1, KODPP2, KRDPP1, KRDPP2
      INTEGER KFREE, LFREE, NTOSYM, ISYMD1, ILLL, KRECNR
      INTEGER IDEL2, IDEL, ISYMD, ISYDIS, KXINT, ISYMTI, ISYMTI2
      INTEGER KDSRHF, K3OINT, NTOT, NUMDIS, IOPT
      INTEGER KLAMPB, KLAMHB, KT1AM, ISYMLP, K3OINT2
      INTEGER KOOOO, ISYMB2
      INTEGER INDEXA(MXCORB)
      INTEGER ISYMC2,ISYMC,IDLSTC,ISYMBC,MAXBC,KLAMPC,KLAMHC,IDUM
C
#if defined (SYS_CRAY)
      REAL XOOOO(*), XOOVV(*), XOVVO(*), WORK(LWORK)
      REAL TIMTOT, TIRDAO, TIME1O, TIME3O
      REAL TIMEINTDEL, TIME2O2V, TIMHE1, TIMHE2
      REAL DTIME, ZERO, ONE
      real ddot ,xnormval
#else
      DOUBLE PRECISION XOOOO(*), XOOVV(*), XOVVO(*),WORK(LWORK)
      DOUBLE PRECISION TIMTOT, TIRDAO, TIME1O, TIME3O
      DOUBLE PRECISION TIMEINTDEL, TIME2O2V, TIMHE1, TIMHE2
      DOUBLE PRECISION DTIME, ZERO, ONE
      DOUBLE PRECISION ddot ,xnormval
#endif
C
      LOGICAL DO_LAMB,DO_LAMC
C
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
C
      CHARACTER*(*) FN4V, LISTB, LISTC
      CHARACTER*12 FNIN3O, FNIN3O2, FN3V, FN3VB
      CHARACTER*10 MODEL
C
      CALL QENTER('CC3_INT')
C
C     DO_LAMC can only be TRUE if DO_LAMB is TRUE !!!
C
      IF ( (DO_LAMC) .AND. (.NOT.DO_LAMB) ) THEN
         WRITE(LUPRI,*)'DO_LAMC can only be TRUE if DO_LAMB is TRUE !!!'
         WRITE(LUPRI,*)'DO_LAMC = ', DO_LAMC
         WRITE(LUPRI,*)'DO_LAMB = ', DO_LAMB
         CALL QUIT('Incorrect logic in CC3_INT')
      END IF
    
      TIMTOT     = SECOND()
      TIRDAO     = ZERO
      TIME1O     = ZERO
      TIME3O     = ZERO
      TIMEINTDEL = ZERO
      TIME2O2V   = ZERO
      TIME2O2V   = ZERO
      TIMHE1     = ZERO
      TIMHE2     = ZERO
C
      IF (.NOT. DO_LAMB) THEN
         ISYMB2 = ISYMOP
      ELSE
         ISYMB2 = ISYMB
      ENDIF
C
      IF (.NOT. DO_LAMC) THEN
         ISYMC2 = ISYMOP
      ELSE
         ISYMC2 = ISYMC
      ENDIF
C
      ISYMBC = MULD2H(ISYMB2,ISYMC2)
      IF (ISYMBC .NE. ISYINT) THEN
         WRITE(LUPRI,*)'ISYMBC and ISYINT should be the same'
         WRITE(LUPRI,*)'ISYMBC = ',ISYMBC
         WRITE(LUPRI,*)'ISYINT = ',ISYINT
         CALL QUIT('Symmetry mismatch in CC3_INT')
      END IF
C
C-----------------------------
C     Work-space allocation 1.
C-----------------------------
C
      MAXBC  = MAX(NT1AM(ISYMB2),NT1AM(ISYMC2))
      KLAMDP = 1
      KLAMDH = KLAMDP + NLAMDT
      KLAMPB = KLAMDH + NLAMDT
      KLAMHB = KLAMPB + NGLMDT(ISYMB2)
      KLAMPC = KLAMHB + NGLMDT(ISYMB2)
      KLAMHC = KLAMPC + NGLMDT(ISYMC2)
      KT1AM  = KLAMHC + NGLMDT(ISYMC2)
      KEND1  = KT1AM  + MAX(NT1AM(ISYMOP),MAXBC)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CC3_OVVO_INT')
      ENDIF
C
C-------------------------------
C     Open scratch files
C-------------------------------
C
      LUIN3O  = -1
      LU3V    = -1
C
      FNIN3O = 'CC3_INT_TMP1'
      FN3V   = 'CC3_INT_TMP2'
C
      CALL WOPEN2(LUIN3O,FNIN3O,64,0)
      IF (LVVVV) THEN
         CALL WOPEN2(LU3V,FN3V,64,0)
      END IF
C
      IF (DO_LAMB) THEN
         LUIN3O2 = -1
         FNIN3O2 = 'CC3_INT_TMP3'
         CALL WOPEN2(LUIN3O2,FNIN3O2,64,0)
      ENDIF
C
C--------------------------------------------
C     Calculate the zero'th lamda matrices.
C--------------------------------------------
C
      IOPT = 1
      CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KT1AM),DUMMY)
C
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),WORK(KEND1),
     &            LWRK1)
C
C--------------------------------------------
C     Calculate the "B" lamda matrices.
C--------------------------------------------
C
      IF (DO_LAMB) THEN
         IOPT = 1
         CALL CC_RDRSP(LISTB,IDLSTB,ISYMB2,IOPT,MODEL,
     *                 WORK(KT1AM),DUMMY)
C

         CALL CCLR_LAMTRA(WORK(KLAMDP),WORK(KLAMPB),WORK(KLAMDH),
     *                    WORK(KLAMHB),WORK(KT1AM),ISYMB2)
      ELSE
         CALL DCOPY(NLAMDT,WORK(KLAMDP),1,WORK(KLAMPB),1)
         CALL DCOPY(NLAMDT,WORK(KLAMDH),1,WORK(KLAMHB),1)
      ENDIF
C
C--------------------------------------------
C     Calculate the "C" lamda matrices.
C--------------------------------------------
C
      IF (DO_LAMC) THEN
         IOPT = 1
         CALL CC_RDRSP(LISTC,IDLSTC,ISYMC2,IOPT,MODEL,
     *                 WORK(KT1AM),DUMMY)
C

         CALL CCLR_LAMTRA(WORK(KLAMDP),WORK(KLAMPC),WORK(KLAMDH),
     *                    WORK(KLAMHC),WORK(KT1AM),ISYMC2)
      ELSE
         CALL DCOPY(NLAMDT,WORK(KLAMDP),1,WORK(KLAMPC),1)
         CALL DCOPY(NLAMDT,WORK(KLAMDH),1,WORK(KLAMHC),1)
      ENDIF
C
C------------------------------------------
C     Regain work space from T1-amplitudes.
C------------------------------------------
C
      KEND1 = KT1AM
      LWRK1 = LWORK - KEND1
C
C-----------------------------------
C     Start the loop over integrals.
C-----------------------------------
C
      KENDS2 = KEND1
      LWRKS2 = LWRK1
C
      IF (DIRECT) THEN
         DTIME  = SECOND()
         IF (HERDIR) THEN
            CALL HERDI1(WORK(KEND1),LWRK1,IPRERI)
         ELSE
            KCCFB1 = KEND1
            KINDXB = KCCFB1 + MXPRIM*MXCONT
            KEND1  = KINDXB + (8*MXSHEL*MXCONT + 1)/IRAT
            LWRK1  = LWORK  - KEND1
            CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     *                  KODPP1,KODPP2,KRDPP1,KRDPP2,
     *                  KFREE,LFREE,KEND1,WORK(KCCFB1),WORK(KINDXB),
     *                  WORK(KEND1),LWRK1,IPRERI)
            KEND1 = KFREE
            LWRK1 = LFREE
         ENDIF
         DTIME  = SECOND() - DTIME
         TIMHE1 = TIMHE1 + DTIME
         NTOSYM = 1
      ELSE
         NTOSYM = NSYM
      ENDIF
C
      KENDSV = KEND1
      LWRKSV = LWRK1
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
C-----------------------------------------------------------------
C           If direct calculate the integrals and transposed CTR2.
C-----------------------------------------------------------------
C
            IF (DIRECT) THEN
C
               KEND1 = KENDSV
               LWRK1 = LWRKSV
C
               DTIME  = SECOND()
               IF (HERDIR) THEN
                  CALL HERDI2(WORK(KEND1),LWRK1,INDEXA,ILLL,NUMDIS,
     &                        IPRINT)
               ELSE
                  CALL ERIDI2(ILLL,INDEXA,NUMDIS,0,0,
     *                        WORK(KODCL1),WORK(KODCL2),
     *                        WORK(KODBC1),WORK(KODBC2),
     *                        WORK(KRDBC1),WORK(KRDBC2),
     *                        WORK(KODPP1),WORK(KODPP2),
     *                        WORK(KRDPP1),WORK(KRDPP2),
     *                        WORK(KCCFB1),WORK(KINDXB),
     *                        WORK(KEND1),LWRK1,IPRERI)
               ENDIF
               DTIME  = SECOND() - DTIME
               TIMHE2 = TIMHE2 + DTIME
C
               KRECNR = KEND1
               KEND1  = KRECNR + (NBUFX(0) - 1)/IRAT + 1
               LWRK1  = LWORK  - KEND1
               IF (LWRK1 .LT. 0) THEN
                  CALL QUIT('Insufficient core in CC3_INT')
               END IF
C
            ELSE
               NUMDIS = 1
            ENDIF
C
C-----------------------------------------------------
C           Loop over number of distributions in disk.
C-----------------------------------------------------
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
C
C------------------------------------------
C              Work space allocation no. 2.
C------------------------------------------
C
               KXINT  = KEND1
               KEND2  = KXINT + NDISAO(ISYDIS)
               LWRK2  = LWORK - KEND2
C
               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
                  CALL QUIT('Insufficient space in CC3_INT')
               ENDIF
C
C--------------------------------------------
C              Read AO integral distribution.
C--------------------------------------------
C
               DTIME = SECOND()
               CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND2),LWRK2,
     *                     WORK(KRECNR),DIRECT)
               DTIME = SECOND() - DTIME
               TIRDAO = TIRDAO + DTIME
C
C------------------------------------------
C              Work space allocation no. 3.
C------------------------------------------
C
               ISYMTI  = MULD2H(ISYMD,ISYMC2)
               ISYMTI2 = MULD2H(ISYMD,ISYMB2)
C
               KDSRHF = KEND2
               K3OINT = KDSRHF + NDSRHF(ISYMD)
               KEND3  = K3OINT + NMAIJK(ISYMTI)
               IF (DO_LAMB) THEN
                  K3OINT2 = KEND3
                  KEND3   = K3OINT2 + NMAIJK(ISYMTI2)
               ENDIF
               LWRK3  = LWORK  - KEND3
C
               IF (LWRK3 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND3,'Available : ',LWORK
                  CALL QUIT('Insufficient space in CC3_INT')
               ENDIF
C
C-----------------------------------------------------------------------
C              Transform one index in the integral batch.
C              (alpha beta | j delta) meaning that we only need lam^0
C-----------------------------------------------------------------------
C
               DTIME = SECOND()
               ISYMLP = 1
               CALL CCTRBT(WORK(KXINT),WORK(KDSRHF),WORK(KLAMDP),
     *                    ISYMLP,WORK(KEND3),LWRK3,ISYDIS)
               DTIME = SECOND() - DTIME
               TIME1O = TIME1O + DTIME
C
C-------------------------------------------------------------------
C              Calculate integral batch with three occupied indices.
C-------------------------------------------------------------------
C
               DTIME = SECOND()
               CALL CC_INT3O(WORK(K3OINT),WORK(KDSRHF),WORK(KLAMHC),
     *                       ISYMC2,WORK(KLAMDP),WORK(KEND3),LWRK3,
     *                       IDEL,ISYMD,LUIN3O,FNIN3O)
               DTIME = SECOND() - DTIME
               TIME3O = TIME3O + DTIME
C
               IF (DO_LAMB) THEN
                  DTIME = SECOND()
                  CALL CC_INT3O(WORK(K3OINT2),WORK(KDSRHF),WORK(KLAMHB),
     *                          ISYMB2,WORK(KLAMDP),WORK(KEND3),LWRK3,
     *                          IDEL,ISYMD,LUIN3O2,FNIN3O2)
c          xnormval = ddot(NMAIJK(ISYMTI2),WORK(K3OINT2),1,
c    *                                     WORK(K3OINT2),1)
c          write(lupri,*)'IDEL,ISYMD,luin3o2',IDEL,ISYMD
c          write(lupri,*)'WORK(K3OINT2)',xnormval

                  DTIME = SECOND() - DTIME
                  TIME3O = TIME3O + DTIME
               ENDIF
C
C--------------------------------------------------------------------
C              Calculate integrals for CC3.
C              Integrals with 3 virtual indices are stored on disc
C--------------------------------------------------------------------
C
               IF (LVVVV) THEN
                 DTIME = SECOND()
                 IF (.NOT. DO_LAMB) THEN
                   CALL CC3_INTDEL(WORK(KXINT),1,LU3V,FN3V,WORK(KLAMDP),
     *                             1,WORK(KLAMDH),1,1,
     *                             WORK(KEND3),LWRK3,IDEL,ISYMD)
                 ELSE
                   CALL CC3_INTDEL2(WORK(KXINT),ISYMOP,LU3V,FN3V,
     *                              WORK(KLAMPC),ISYMC2,
     *                              WORK(KLAMDH),ISYMOP,
     *                              WORK(KLAMPB),ISYMB2,
     *                              WORK(KLAMHB),ISYMB2,
     *                              ISYINT,WORK(KEND3),LWRK3,IDEL,ISYMD)
                 ENDIF
                 DTIME = SECOND() - DTIME
                 TIMEINTDEL = TIMEINTDEL + DTIME
               END IF
C
               DTIME = SECOND()
               CALL CC3_2O2V(WORK(KXINT),ISYMOP,WORK(KDSRHF),ISYMOP,
     *                       XOVVO,XOOVV,WORK(KLAMDP),ISYMOP,
     *                       WORK(KLAMDH),ISYMOP,
     *                       WORK(KLAMPB),ISYMB2,
     *                       WORK(KLAMHC),ISYMC2,ISYINT,
     *                       WORK(KEND3),LWRK3,IDEL,ISYMD)
C
               IF (DO_LAMB) THEN
                  CALL CC3_2O2V(WORK(KXINT),ISYMOP,WORK(KDSRHF),ISYMOP,
     *                          XOVVO,XOOVV,WORK(KLAMDP),ISYMOP,
     *                          WORK(KLAMDH),ISYMOP,
     *                          WORK(KLAMPC),ISYMC2,
     *                          WORK(KLAMHB),ISYMB2,ISYINT,
     *                          WORK(KEND3),LWRK3,IDEL,ISYMD)
               ENDIF
C
               DTIME = SECOND() - DTIME
               TIME2O2V = TIME2O2V + DTIME
C              
C---------------------------------------
C     END ALL LOOPS
C---------------------------------------
C
            ENDDO    ! IDEL2
         ENDDO       ! ILLL
      ENDDO          ! ISYMD1
C
C------------------------
C     Recover work space.
C------------------------
C
      KEND1 = KENDS2
      LWRK1 = LWRKS2
C
      CALL DZERO(XOOOO,N3ORHF(ISYINT))
C
C     for DO_LAMB lamH^B is in KLAMHB otherwise lamH^0
C
      IF (LVVVV) THEN
         IOPT = 3
      ELSE
         IOPT = 1
      END IF
      CALL CC3_INTSTORE(LUIN3O,FNIN3O,XOOOO,ISYINT,
     *                  WORK(KLAMHB),ISYMB2,WORK(KLAMDH),ISYMOP,
     *                  LU3V,FN3V,LU4V,FN4V,ISYINT,
     *                  WORK(KEND1),LWRK1,IOPT)
C
      IF (DO_LAMB) THEN
         KOOOO = KEND1
         KEND1 = KOOOO + N3ORHF(ISYINT)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Out of memory (koooo) in CC3_INT')
         ENDIF
C
         CALL DZERO(WORK(KOOOO),N3ORHF(ISYINT))
C
         IOPT = 1
         CALL CC3_INTSTORE(LUIN3O2,FNIN3O2,WORK(KOOOO),ISYINT,
     *                     WORK(KLAMHC),ISYMC2,WORK(KLAMDH),ISYMOP,
     *                     LU3V,FN3V,LU4V,FN4V,ISYINT,
     *                     WORK(KEND1),LWRK1,IOPT)
C
         CALL DAXPY(N3ORHF(ISYINT),ONE,WORK(KOOOO),1,XOOOO,1)
      ENDIF
C
      CALL CC3_SORT4O(XOOOO,ISYINT,WORK(KEND1),LWRK1)
C
C-------------------------------
C     Delete intermediate files.
C-------------------------------
C
      CALL WCLOSE2(LUIN3O,FNIN3O,'DELETE')
      IF (LVVVV) THEN
         CALL WCLOSE2(LU3V,FN3V,'DELETE')
      END IF
      IF (DO_LAMB) THEN
         CALL WCLOSE2(LUIN3O2,FNIN3O2,'DELETE')
      ENDIF
C
      TIMTOT = SECOND() - TIMTOT
C
C-------------------------------
C     Write out program timings.
C-------------------------------
C
      IF (IPRINT .GT. 3) THEN
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) 'Time used in CC3_INT :',TIMTOT
         WRITE(LUPRI,*) ' '
      ENDIF
      IF (IPRINT .GT. 9) THEN
         WRITE(LUPRI,*) 'Time used as follows:'
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) 'Time used in CCINT3O:',TIME3O
      ENDIF
C
   1  FORMAT(1x,A35,1X,E20.10)
C
C
      CALL QEXIT('CC3_INT')
      RETURN
      END
C  /* Deck cc3_intdel2 */
      SUBROUTINE CC3_INTDEL2(AOINT,ISYMAO,LUINT,FNINT,
     *                       XLAMP0,ISYMLP0,XLAMH0,ISYMLH0,
     *                       XLAMPB,ISYMLPB,XLAMHB,ISYMLHB,
     *                       ISYINT,WORK,LWORK,IDEL,ISYMD)
C
C     Written by K. Hald, Summer 2002.
C
C     Calculate integrals that are needed for the CC3 left hand side,
C     and store on file.
C
C     VVV,delta (V=vir.) are needed ... 
C     equals (v-bar v|v delta) + (vv|v-bar delta)
C
C
      IMPLICIT NONE
C
      INTEGER ISYMAO, LUINT, ISYMLP0, ISYMLH0, ISYMLPB, ISYMLHB
      INTEGER ISYINT, LWORK, IDEL, ISYMD
      INTEGER ISYABG, ISYTMP, ISYABC, KVVVV, KEND1, KEND2, LWRK1, LWRK2
      INTEGER ISYMG, ISYMC, ISALBE, ISYMAB, KINT, KSCR1, KSCR2
      INTEGER KOFF1, KOFF2, KOFF3, ISYMB, ISYMBE, ISYMAL, ISYMA
      INTEGER NBASAL, NBASBE, NVIRA, NAB, NBASG, IOFF
C
#if defined (SYS_CRAY)
      REAL AOINT(*), XLAMP0(*), XLAMH0(*)
      REAL XLAMPB(*), XLAMHB(*)
      REAL WORK(LWORK), ZERO, ONE
#else
      DOUBLE PRECISION AOINT(*), XLAMP0(*), XLAMH0(*)
      DOUBLE PRECISION XLAMPB(*), XLAMHB(*)
      DOUBLE PRECISION WORK(LWORK), ZERO, ONE
#endif
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      CHARACTER FNINT*(*)
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C
      CALL QENTER('CC3_INTDEL')
C
C-------------------------------------------
C     Work space allocation.
C-------------------------------------------
C
      ISYABG = MULD2H(ISYMAO,ISYMD)
C
      ISYTMP = MULD2H(ISYINT,ISYMD)
      ISYABC = MULD2H(ISYTMP,ISYMLH0)
C
      KVVVV  = 1
      KEND1  = KVVVV  + NMAABC(ISYABC)
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory in CC3_INTDEL')
      ENDIF
C
      CALL DZERO(WORK(KVVVV),NMAABC(ISYABC))
C
C-----------------------------------------------------
C     Transform AO-integrals to g_{v-barvv,delta}
C-----------------------------------------------------
C
      DO ISYMG = 1, NSYM
         ISYMC  = MULD2H(ISYMG,ISYMLP0)
         ISALBE = MULD2H(ISYABG,ISYMG)
         ISYMAB = MULD2H(ISYABC,ISYMC)
         ISYTMP = MULD2H(ISYMAB,ISYMLH0)
C
         KINT   = KEND1
         KSCR1  = KINT  + NMATAB(ISYMAB)*NBAS(ISYMG)
         KSCR2  = KSCR1 + N2BST(ISALBE)
         KEND2  = KSCR2 + NEMAT1(ISYTMP)
         LWRK2  = LWORK - KEND2
COMMENT
COMMENT  allocate to much space for kscr2 at the moment
COMMENT
C
         IF (LWRK2 .LT. 0) THEN
            CALL QUIT('Out of memory in CC3_INTDEL (2)')
         ENDIF
C
         DO G = 1, NBAS(ISYMG)
C
            KOFF1 = IDSAOG(ISYMG,ISYMD) + NNBST(ISALBE)*(G-1) + 1
            CALL CCSD_SYMSQ(AOINT(KOFF1),ISALBE,WORK(KSCR1))
C
            DO ISYMB = 1,NSYM
C
               ISYMBE = MULD2H(ISYMB,ISYMLH0)
               ISYMAL = MULD2H(ISYMBE,ISALBE)
               ISYMA  = MULD2H(ISYMAL,ISYMLPB)
C
               KOFF1 = KSCR1 
     *               + IAODIS(ISYMAL,ISYMBE)
               KOFF2 = IGLMVI(ISYMBE,ISYMB) + 1
               KOFF3 = KSCR2
C
               NBASAL = MAX(NBAS(ISYMAL),1)
               NBASBE = MAX(NBAS(ISYMBE),1)
C
               CALL DGEMM('N','N',NBAS(ISYMAL),NVIR(ISYMB),NBAS(ISYMBE),
     *                    ONE,WORK(KOFF1),NBASAL,XLAMH0(KOFF2),NBASBE,
     *                    ZERO,WORK(KOFF3),NBASAL)
C
               KOFF1 = IGLMVI(ISYMAL,ISYMA) + 1
               KOFF2 = KSCR2
               KOFF3 = KINT 
     *               + NMATAB(ISYMAB)*(G - 1)
     *               + IMATAB(ISYMA,ISYMB)
C
               NBASAL = MAX(NBAS(ISYMAL),1)
               NVIRA  = MAX(NVIR(ISYMA),1)
C
               CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),NBAS(ISYMAL),
     *                    ONE,XLAMPB(KOFF1),NBASAL,WORK(KOFF2),NBASAL,
     *                    ZERO,WORK(KOFF3),NVIRA)
C
            ENDDO
C
         ENDDO
C
         KOFF2 = IGLMVI(ISYMG,ISYMC)  + 1
         KOFF3 = KVVVV 
     *         + IMAABC(ISYMAB,ISYMC)
C
         NAB    = MAX(NMATAB(ISYMAB),1)
         NBASG  = MAX(NBAS(ISYMG),1)
C
         CALL DGEMM('N','N',NMATAB(ISYMAB),NVIR(ISYMC),NBAS(ISYMG),
     *              ONE,WORK(KINT),NAB,XLAMP0(KOFF2),NBASG,
     *              ONE,WORK(KOFF3),NAB)
C
      ENDDO
C
C-----------------------------------------------------
C     Transform AO-integrals to g_{vvv-bar,delta}
C-----------------------------------------------------
C
      DO ISYMG = 1, NSYM
         ISYMC  = MULD2H(ISYMG,ISYMLPB)
         ISALBE = MULD2H(ISYABG,ISYMG)
         ISYMAB = MULD2H(ISYABC,ISYMC)
         ISYTMP = MULD2H(ISYMAB,ISYMLH0)
C
         KINT   = KEND1
         KSCR1  = KINT  + NMATAB(ISYMAB)*NBAS(ISYMG)
         KSCR2  = KSCR1 + N2BST(ISALBE)
         KEND2  = KSCR2 + NEMAT1(ISYTMP)
         LWRK2  = LWORK - KEND2
COMMENT
COMMENT  allocate to much space for kscr2 at the moment
COMMENT
C
         IF (LWRK2 .LT. 0) THEN
            CALL QUIT('Out of memory in CC3_INTDEL (2)')
         ENDIF
C
         DO G = 1, NBAS(ISYMG)
C
            KOFF1 = IDSAOG(ISYMG,ISYMD) + NNBST(ISALBE)*(G-1) + 1
            CALL CCSD_SYMSQ(AOINT(KOFF1),ISALBE,WORK(KSCR1))
C
            DO ISYMB = 1,NSYM
C
               ISYMBE = MULD2H(ISYMB,ISYMLH0)
               ISYMAL = MULD2H(ISYMBE,ISALBE)
               ISYMA  = MULD2H(ISYMAL,ISYMLP0)
C
               KOFF1 = KSCR1 
     *               + IAODIS(ISYMAL,ISYMBE)
               KOFF2 = IGLMVI(ISYMBE,ISYMB) + 1
               KOFF3 = KSCR2
C
               NBASAL = MAX(NBAS(ISYMAL),1)
               NBASBE = MAX(NBAS(ISYMBE),1)
C
               CALL DGEMM('N','N',NBAS(ISYMAL),NVIR(ISYMB),NBAS(ISYMBE),
     *                    ONE,WORK(KOFF1),NBASAL,XLAMH0(KOFF2),NBASBE,
     *                    ZERO,WORK(KOFF3),NBASAL)
C
               KOFF1 = IGLMVI(ISYMAL,ISYMA) + 1
               KOFF2 = KSCR2
               KOFF3 = KINT 
     *               + NMATAB(ISYMAB)*(G - 1)
     *               + IMATAB(ISYMA,ISYMB)
C
               NBASAL = MAX(NBAS(ISYMAL),1)
               NVIRA  = MAX(NVIR(ISYMA),1)
C
               CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),NBAS(ISYMAL),
     *                    ONE,XLAMP0(KOFF1),NBASAL,WORK(KOFF2),NBASAL,
     *                    ZERO,WORK(KOFF3),NVIRA)
C
            ENDDO
C
         ENDDO
C
         KOFF2 = IGLMVI(ISYMG,ISYMC)  + 1
         KOFF3 = KVVVV 
     *         + IMAABC(ISYMAB,ISYMC)
C
         NAB    = MAX(NMATAB(ISYMAB),1)
         NBASG  = MAX(NBAS(ISYMG),1)
C
         CALL DGEMM('N','N',NMATAB(ISYMAB),NVIR(ISYMC),NBAS(ISYMG),
     *              ONE,WORK(KINT),NAB,XLAMPB(KOFF2),NBASG,
     *              ONE,WORK(KOFF3),NAB)
C
      ENDDO
C
C----------------------------------------
C     Save the g_{vvv,delta} to disc.
C----------------------------------------
C
      IF (NMAABC(ISYABC) .GT. 0) THEN
         KOFF1 = IDEL - IBAS(ISYMD)
         IOFF = I3VDEL(ISYABC,ISYMD) + NMAABC(ISYABC)*(KOFF1-1) + 1
         CALL PUTWA2(LUINT,FNINT,WORK(KVVVV),IOFF,NMAABC(ISYABC))
      ENDIF
C
C--------------------------
C     End.
C--------------------------
C
      CALL QEXIT('CC3_INTDEL')
C
      RETURN
      END
C  /* Deck cc3_intdel2 */
      SUBROUTINE CC3_BARINT(B1AM,ISYMB,XLAMDP,XLAMDH,WORK,LWORK,
     *                      LU3SRTR,FN3SRTR,LUCKJDR,FNCKJDR)
C
C     Written by Kasper Hald, Summer 2002.
C
C     Interface to calculate the (ck|de)-bar and (ck|lm)-bar integrals
C
C     NSIMTR = 1 !!!!!!!
C
      IMPLICIT NONE
C
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
#include <ccsdsym.h>
#include <distcl.h>
#include <cbieri.h>
#include <eribuf.h>
C#include <cclr.h>
C
      INTEGER ISYMB, LWORK, LU3SRTR, LUCKJDR
      INTEGER KEND1, LWRK1, KEND2, LWRK2, KENDSV, LWRKSV
      INTEGER KCCFB1, KINDXB, KODCL1, KODCL2, KODBC1, KODBC2
      INTEGER KRDBC1, KRDBC2, KODPP1, KODPP2, KRDPP1, KRDPP2
      INTEGER KFREE, LFREE, NTOSYM, ISYMD1, NTOT, ILLL
      INTEGER KRECNR, NUMDIS, IDEL2, IDEL, ISYMD, ISYDIS
      INTEGER ISYAIK, KXINT, KLAMPB, KLAMHB, ISYCKJ, KCKJ
      INTEGER INDEXA(MXCORB),IDUM
C
#if defined (SYS_CRAY)
#else
      DOUBLE PRECISION B1AM(*), XLAMDP(*), XLAMDH(*), WORK(LWORK)
#endif
C
      CHARACTER*(*) FN3SRTR, FNCKJDR
      CHARACTER*1 CDUMMY
C
      CALL QENTER('CC3_BARINT')
C
      CDUMMY = ' '
C
      KLAMPB = 1
      KLAMHB = KLAMPB + NLAMDT
      KEND1  = KLAMHB + NLAMDT
      LWRK1 = LWORK
C
      KENDSV = KEND1
      LWRKSV = LWRK1
C
C--------------------------------------
C     Calculate lambda-B
C--------------------------------------
C
      CALL CCLR_LAMTRA(XLAMDP,WORK(KLAMPB),XLAMDH,
     *                 WORK(KLAMHB),B1AM,ISYMB)
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
               KEND1 = KENDSV
               LWRK1 = LWRKSV
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
                  CALL QUIT('Insufficient core in CC3_BARINT')
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
                  CALL QUIT('Insufficient space in CC3_FMAT ')
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
C              Calculate transformed integrals used in t3am.
C-----------------------------------------------------------
C
               CALL CC3_T3INT(WORK(KXINT),XLAMDP,XLAMDH,
     *                        B1AM,ISYMB,WORK(KEND2),LWRK2,
     *                        IDEL,ISYMD,2,LU3SRTR,FN3SRTR,
     *                        LUCKJDR,FNCKJDR)
C
            ENDDO     ! IDEL2
         ENDDO        ! ILLL
      ENDDO           ! ISYMD1
C
      CALL QEXIT('CC3_BARINT')
      RETURN
      END
C  /* Deck wbd_ground */
      SUBROUTINE WBD_GROUND(T2TP,ISYMT2,TMAT,TRVIR,TRVIR2,TROCC,ISYINT,
     *                      WMAT,WORK,LWORK,INDSQ,LENSQ,
     *                      ISYMB,B,ISYMC,C)
C
C     Kasper Hald, Summer 2002
C
C     WBD(ai,k,j) = WBD(ai,k,j)
C                 + t(ai,ej) (Dk|Be) + t(ai,ej) (Bk|De)
C                 - t(ai,Bl) (Dk|lj) - t(ai,Dl) (Bk|lj)
C
C     ISYINT is the symmetry of the integrals.
C     ISYMT2 is the symmetry of T2TP.
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C#include <cclr.h>
C
      INTEGER ISYMT2, ISYINT, LWORK, LENSQ, ISYMB, ISYMC
      INTEGER INDSQ(LENSQ,6)
      INTEGER ISYMBC, ISYRES, JSAIKJ, ISYMDK, LENGTH, ISYMK
      INTEGER ISYMD, ISYAIJ, KOFF1, KOFF2, KOFF3, NTOAIJ, NVIRD
      INTEGER ISYAIL, ISYLKJ, ISYMJ, ISYMLK, ISYML, ISYMAI
      INTEGER ISYAIK, NTOTAI, NRHFL
C
      integer isymdkb
C
#if defined (SYS_CRAY)
      REAL T2TP(*), TRVIR(*), TRVIR2(*), TROCC(*)
      REAL TMAT(*), WMAT(*), WORK(LWORK)
      REAL ZERO, ONE, DDOT, XNORM
#else
      DOUBLE PRECISION T2TP(*), TRVIR(*), TRVIR2(*), TROCC(*)
      DOUBLE PRECISION TMAT(*), WMAT(*), WORK(LWORK)
      DOUBLE PRECISION ZERO, ONE, DDOT, XNORM
#endif
C
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
C
      CALL QENTER('WBD_GROUND')

C
C--------------------------------
C     First virtual contribution.
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
         CALL QUIT('Insufficient core in WBD_GROUND (1)')
      ENDIF
C
      DO ISYMK = 1,NSYM
C
         ISYMD  = MULD2H(ISYMK,ISYMDK)
         ISYAIJ = MULD2H(ISYMK,JSAIKJ)
C
         KOFF1 = IT2SP(ISYAIJ,ISYMD)  + 1
         KOFF2 = ICKATR(ISYMDK,ISYMB) + NT1AM(ISYMDK)*(B - 1)
     *         + IT1AM(ISYMD,ISYMK)   + 1
         KOFF3 = ISAIKJ(ISYAIJ,ISYMK) + 1
C
         NTOAIJ = MAX(1,NCKI(ISYAIJ))
         NVIRD  = MAX(NVIR(ISYMD),1)
C
         CALL DGEMM('N','N',NCKI(ISYAIJ),NRHF(ISYMK),
     *              NVIR(ISYMD),ONE,T2TP(KOFF1),NTOAIJ,
     *              TRVIR(KOFF2),NVIRD,ZERO,
     *              WORK(KOFF3),NTOAIJ)
C
      ENDDO
C

      DO I = 1,LENGTH
         WMAT(I) = WMAT(I) + WORK(INDSQ(I,3))
      ENDDO

C
      IF (IPRINT .GT. 55) THEN
         XNORM = DDOT(NCKIJ(JSAIKJ),WMAT,1,WMAT,1)
         WRITE(LUPRI,*) 'In WBD_GROUND: 1. Norm of WMAT ',XNORM
      ENDIF
C
C---------------------------------
C     Second virtual contribution.
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
         CALL QUIT('Insufficient core in WBD_GROUND (2)')
      ENDIF
C
      DO ISYMK = 1,NSYM
C
         ISYMD  = MULD2H(ISYMK,ISYMDK)
         ISYAIJ = MULD2H(ISYMK,JSAIKJ)
C
         KOFF1 = IT2SP(ISYAIJ,ISYMD)  + 1
         KOFF2 = ICKATR(ISYMDK,ISYMC) + NT1AM(ISYMDK)*(C - 1)
     *         + IT1AM(ISYMD,ISYMK)   + 1
         KOFF3 = ISAIKJ(ISYAIJ,ISYMK) + 1
C
         NTOAIJ = MAX(1,NCKI(ISYAIJ))
         NVIRD  = MAX(NVIR(ISYMD),1)
C
         CALL DGEMM('N','N',NCKI(ISYAIJ),NRHF(ISYMK),
     *              NVIR(ISYMD),ONE,T2TP(KOFF1),NTOAIJ,
     *              TRVIR2(KOFF2),NVIRD,ZERO,
     *              WORK(KOFF3),NTOAIJ)
C
      ENDDO
C
      DO I = 1,LENGTH
         WMAT(I) = WMAT(I) + WORK(I)
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XNORM = DDOT(NCKIJ(JSAIKJ),WMAT,1,WMAT,1)
         WRITE(LUPRI,*) 'In WBD_GROUND: 2. Norm of WMAT ',XNORM
      ENDIF
C
C---------------------------------
C     First occupied contribution.
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
     *               + NCKI(ISYAIL)*(B - 1)
     *               + ICKI(ISYMAI,ISYML) + 1
               KOFF2 = ISJIKA(ISYLKJ,ISYMC)
     *               + NMAJIK(ISYLKJ)*(C - 1)
     *               + ISJIK(ISYMLK,ISYMJ)
     *               + NMATIJ(ISYMLK)*(J - 1)
     *               + IMATIJ(ISYML,ISYMK) + 1
               KOFF3 = ISAIKJ(ISYAIK,ISYMJ)
     *               + NCKI(ISYAIK)*(J - 1)
     *               + ICKI(ISYMAI,ISYMK) + 1
C
               NTOTAI = MAX(1,NT1AM(ISYMAI))
               NRHFL  = MAX(1,NRHF(ISYML))
C
               CALL DGEMM('N','N',NT1AM(ISYMAI),NRHF(ISYMK),
     *                    NRHF(ISYML),-ONE,T2TP(KOFF1),NTOTAI,
     *                    TROCC(KOFF2),NRHFL,ONE,WMAT(KOFF3),
     *                    NTOTAI)
C
            ENDDO
         ENDDO
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XNORM = DDOT(NCKIJ(JSAIKJ),WMAT,1,WMAT,1)
         WRITE(LUPRI,*) 'In WBD_GROUND: 3. Norm of WMAT ',XNORM
      ENDIF
C
C----------------------------------
C     Second occupied contribution.
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
     *               + NCKI(ISYAIL)*(C - 1)
     *               + ICKI(ISYMAI,ISYML) + 1
               KOFF2 = ISJIKA(ISYLKJ,ISYMB)
     *               + NMAJIK(ISYLKJ)*(B - 1)
     *               + ISJIK(ISYMLK,ISYMJ)
     *               + NMATIJ(ISYMLK)*(J - 1)
     *               + IMATIJ(ISYML,ISYMK) + 1
               KOFF3 = ISAIKJ(ISYAIK,ISYMJ)
     *               + NCKI(ISYAIK)*(J - 1)
     *               + ICKI(ISYMAI,ISYMK) + 1
C
               NTOTAI = MAX(1,NT1AM(ISYMAI))
               NRHFL  = MAX(1,NRHF(ISYML))
C
               CALL DGEMM('N','N',NT1AM(ISYMAI),NRHF(ISYMK),
     *                    NRHF(ISYML),-ONE,T2TP(KOFF1),NTOTAI,
     *                    TROCC(KOFF2),NRHFL,ZERO,TMAT(KOFF3),
     *                    NTOTAI)
C
            ENDDO
         ENDDO
      ENDDO
C
      DO I = 1,LENGTH
         WMAT(I) = WMAT(I) + TMAT(INDSQ(I,3))
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XNORM = DDOT(NCKIJ(JSAIKJ),WMAT,1,WMAT,1)
         WRITE(LUPRI,*) 'In WBD_GROUND: 4. Norm of WMAT ',XNORM
      ENDIF
C
C---------------------------------------
C     End.
C---------------------------------------
C
      CALL QEXIT('WBD_GROUND')
C
      RETURN
      END
C  /* Deck cclr_trvir */
      SUBROUTINE CCLR_TRVIR(XINT,TRINT,XLAMDP,ISYMLP,ISYMD,D,ISYINT,
     *                      WORK,LWORK)
C
C     Henrik Koch and Alfredo Sanchez.         Dec 1994
C
C     Transform (ck|alpha d) integrals to (ck|a d)
C
C     Ove Christiansen 11-1-1996: ISYINT is symmetry of (ck|alpha d)
C
C     Kasper Hald, 01062002 : General symmetry of Lambda
C
#include <implicit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION XINT(*),TRINT(*),XLAMDP(*),WORK(LWORK)
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('CCLR_TRVIR')
C
      ISYDIS = MULD2H(ISYMD,ISYINT)
C
      DO 100 ISYMA = 1,NSYM
C
COMMENT COMMENT
COMMENT COMMENT
         ISYMAL = MULD2H(ISYMA,ISYMLP)
C         ISYMCK = MULD2H(ISYMA,ISYDIS)
         ISYMCK = MULD2H(ISYMAL,ISYDIS)
COMMENT COMMENT
COMMENT COMMENT
C
         NTOTCK = MAX(NT1AM(ISYMCK),1)
COMMENT COMMENT
C         NBASA  = MAX(NBAS(ISYMA),1)
         NBASA  = MAX(NBAS(ISYMAL),1)
COMMENT COMMENT
C
COMMENT COMMENT COMMENT
C         KOFF1  = ICKALP(ISYMCK,ISYMA) + 1
         KOFF1  = ICKALP(ISYMCK,ISYMAL) + 1
C         KOFF2  = ILMVIR(ISYMA) + 1
         KOFF2  = IGLMVI(ISYMAL,ISYMA)
     *          + 1
COMMENT COMMENT COMMENT
         KOFF3  = ICKATR(ISYMCK,ISYMA) + 1
C
COMMENT COMMENT
C         CALL DGEMM('N','N',NT1AM(ISYMCK),NVIR(ISYMA),NBAS(ISYMA),
         CALL DGEMM('N','N',NT1AM(ISYMCK),NVIR(ISYMA),NBAS(ISYMAL),
COMMENT COMMENT
     *              ONE,XINT(KOFF1),NTOTCK,XLAMDP(KOFF2),NBASA,
     *              ZERO,TRINT(KOFF3),NTOTCK)
C
  100 CONTINUE
C
      CALL QEXIT('CCLR_TRVIR')
C
      RETURN
      END
C  /* Deck cclr_trocc */
      SUBROUTINE CCLR_TROCC(XINT,TRINT,XLAMDH,ISYMLH,WORK,LWORK)
C
C     Henrik Koch and Alfredo Sanchez.         Dec 1994
C
C     Transform (ia|j delta) integrals to (ia|j k) and sort as (ij,k,a)
C
C     Kasper Hald, Summer 2002 -> generalisation to nonsymmetric lambda
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INTEGER ISYMLH, LWORK
      INTEGER ISYMK, ISYINT, ISYAIJ, LENGTH, ISYMD
      INTEGER NTOIAJ, NBASD, KOFF1, KOFF2, ISYMJ, ISYMAI
      INTEGER ISYMI, ISYMA, ISYMJI, ISYJIK, NTOJIK
C
#if defined (SYS_CRAY)
#else
      DOUBLE PRECISION XINT(*), TRINT(*), XLAMDH(*), WORK(LWORK)
      DOUBLE PRECISION ZERO, ONE, DDOT
#endif
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C
C
      CALL QENTER('CCLR_TROCC')
C
      DO ISYMD = 1,NSYM
C
COMMENT COMMENT
C         ISYMK  = ISYMD
C         ISYAIJ = MULD2H(ISYMD,ISYMOP)
         ISYMK  = MULD2H(ISYMD,ISYMLH)
C         ISYINT = MULD2H(ISYMLH,ISYMOP)
C         ISYAIJ = MULD2H(ISYMK,ISYINT)
         ISYAIJ = MULD2H(ISYMD,ISYMOP)
COMMENT COMMENT
C
         LENGTH = NCKI(ISYAIJ) * NRHF(ISYMK)
         IF (LENGTH .GT. LWORK) THEN
            CALL QUIT('Not enough space in CCLR_TROCC')
         END IF
C
         NTOIAJ = MAX(NCKI(ISYAIJ),1)
         NBASD  = MAX(NBAS(ISYMD),1)
C
         KOFF1 = ICKID(ISYAIJ,ISYMD)  + 1
COMMENT COMMENT
COMMENT COMMENT
C         KOFF2 = ILMRHF(ISYMD) + 1
         KOFF2 = IGLMRH(ISYMD,ISYMK) + 1
COMMENT COMMENT
COMMENT COMMENT
C
         CALL DGEMM('N','N',NCKI(ISYAIJ),NRHF(ISYMK),NBAS(ISYMD),
     *              ONE,XINT(KOFF1),NTOIAJ,XLAMDH(KOFF2),NBASD,
     *              ZERO,WORK,NTOIAJ)
C
C        Sort
C
         DO ISYMJ = 1,NSYM
C
            ISYMAI = MULD2H(ISYAIJ,ISYMJ)
C
            DO ISYMI = 1,NSYM
C
               ISYMA  = MULD2H(ISYMAI,ISYMI)
               ISYMJI = MULD2H(ISYMJ,ISYMI)
               ISYJIK = MULD2H(ISYMJI,ISYMK)
C
               DO K = 1,NRHF(ISYMK)
C
                  DO J = 1,NRHF(ISYMJ)
C
                     DO I = 1,NRHF(ISYMI)
C
                        NTOJIK = NMAJIK(ISYJIK)
C
                        KOFF1 = NCKI(ISYAIJ)*(K - 1)
     *                        + ISAIK(ISYMAI,ISYMJ)
     *                        + NT1AM(ISYMAI)*(J - 1)
     *                        + IT1AM(ISYMA,ISYMI)
     *                        + NVIR(ISYMA)*(I - 1) + 1
                        KOFF2 = ISJIKA(ISYJIK,ISYMA)
     *                        + ISJIK(ISYMJI,ISYMK)
     *                        + NMATIJ(ISYMJI)*(K - 1)
     *                        + IMATIJ(ISYMJ,ISYMI)
     *                        + NRHF(ISYMJ)*(I - 1) + J
C
                        CALL DCOPY(NVIR(ISYMA),WORK(KOFF1),1,
     *                             TRINT(KOFF2),NTOJIK)
C
                     ENDDO    ! I
                  ENDDO       ! J
               ENDDO          ! K
            ENDDO             ! ISYMI
         ENDDO                ! ISYMJ
C
      ENDDO                   ! ISYMD
C
      CALL QEXIT('CCLR_TROCC')
C
      RETURN
      END
