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
C  /* Deck ccsdpt_dens2_sc */
      SUBROUTINE CCSDPT_DENS2_SC(T1AM,ISYMT1,T2TP,ISYMT2,MODEL,
     *                        L1AM,ISYML1,L2TP,ISYML2,WORK,LWORK,
     *                        LUDELD,FNDELD,LUCKJD,FNCKJD,LUDKBC,FNDKBC,
     *                        LUTOC,FNTOC,LU3VI,FN3VI,LUDKBC3,FNDKBC3,
     *                        LU3FOP,FN3FOP,LU3FOP2,FN3FOP2,
     *                        LU3FOPX,FN3FOPX,LU3FOP2X,FN3FOP2X)
C
C     Written by K. Hald, Fall 2001.
C     Modified (temporary) version for CCSD(T) gradients with symmetry. 
C     S. Coriani, Fall 2002 
C
C     Calculate the triples contribution to the electronic densities
C     in the MO basis and store them on file.
C     Calculate also the diagonal kappabar multipliers if (RELORB).
C
C     ISYMT2 is symmetry of T2TP
C     ISYMT1 is symmetry of T1AM
C     Isyres = isymt1*isymt2*isymop
C
C     For CCSD(T) LUDKBC3, FNDKBC3 is actually LU3VI2, FN3VI2
C     For CCSD(T) we do not use LUDELD,FNDELD,LUCKJD,FNCKJD,LUDKBC,FNDKBC
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
#include <ccfop.h>
#include <second.h>
C
      INTEGER ISYMTR, ISYMT1, ISYMT2, ISYML1, ISYML2, LWORK
      INTEGER ISYRES, ISINT1, ISINT2, ISYMIM, KFOCKD
      INTEGER KOMG22, KCMO, KEND0, LWRK0, KTROC2
      INTEGER KTROC0, KXIAJB, KEND1, LWRK1, KINTOC, KEND2, LWRK2
      INTEGER LENGTH, ISYOPE, IOPTTCME, IOFF, ISYMD, ISAIJ1, ISYCKB
      INTEGER ISCKB1, ISCKB2, KTRVI1, KTRVI2, KRMAT1, KTRVI0
      INTEGER KTRVI3, KEND3, LWRK3, KINTVI, KEND4, LWRK4, ISYMB
      INTEGER ISYALJ, ISAIJ2, ISYMBD, ISCKIJ, KSMAT2, KSMAT, KQMAT
      INTEGER KDIAG, ISYMC, ISYMK, KOFF1, KOFF2, KOFF3
      INTEGER KINDSQ, KINDEX, KTMAT, KRMAT2, LENSQ
      INTEGER LUFCK, KFCKBA, KT2TCME, IOPTT2, KTRVI4, KTRVI5
      INTEGER KTRVI6, KQMAT2, KVIR1, KVIR2, KVIR3, KVIR4, LUPTIA
      INTEGER LUPTIAJB, LUABI1, LUABI2, LUABI3, LUABI4, ISYAIB
      INTEGER ISYMAI, ISYAID, KOCC1, KOCC2, KOMG1, KUMAT, KUMAT2
      INTEGER LUAIJK, LUIAJK, LUPTAB, LUPTIJ, LUPTIA2
      INTEGER KTROC02, KTROC22, KTRVI7, KTRVI8, KTRVI9, KTRVI10
      INTEGER KTRVI11, KTRVI12, KTRVI13, KDENSAB, KDENSIJ
      INTEGER ISYCKD, ISCKD2, KSMAT3, KUMAT3, KEND5, LWRK5
      INTEGER KKAPAA, KKAPII, LUKAPAB, LUKAPIJ, ISYALJ2, KINDEX2
      INTEGER KSMAT4, KUMAT4, KOMG12, KLAMDP, KLAMDH
      INTEGER KTRVI14, KTRVI15, KTRVI16, KTRVI17, KTRVI18, KTRVI19
      INTEGER KTRVI20, ISYTMP, KTROC01, KTROC21, KTROC03, KTROC23
      INTEGER LUDELD, LUCKJD, LUDKBC, LUTOC, LU3VI, LUDKBC3, LU3FOP
      INTEGER LU3FOP2, LU3FOPX, LU3FOP2X
      integer KENDSV,NTOT
C
#if defined (SYS_CRAY)
      REAL T1AM(*), T2TP(*)
      REAL L1AM(*), L2TP(*)
      REAL WORK(LWORK), ONE
      REAL TITRAN, TISORT, TISMAT, TIQMAT, TIOME1
      REAL RHO1N, RHO2N
      REAL XT2TP, DDOT, XIAJB, XINT, XTROC, XTROC1, XTROC0
      REAL XTRVI0, XTRVI2, XTRVI3, XTRVI, XTRVI1, XDIA
      REAL XSMAT, XTMAT, XQMAT, XRMAT, ZERO, TWO, HALF
      REAL DTIME
#else
      DOUBLE PRECISION T1AM(*), T2TP(*)
      DOUBLE PRECISION L1AM(*), L2TP(*)
      DOUBLE PRECISION WORK(LWORK), ONE
      DOUBLE PRECISION TITRAN, TISORT, TISMAT, TIQMAT, TIOME1
      DOUBLE PRECISION RHO1N, RHO2N
      DOUBLE PRECISION XT2TP, DDOT, XIAJB, XINT, XTROC, XTROC1, XTROC0
      DOUBLE PRECISION XTRVI0, XTRVI2, XTRVI3, XTRVI, XTRVI1, XDIA
      DOUBLE PRECISION XSMAT, XTMAT, XQMAT, XRMAT, ZERO, TWO, HALF
      DOUBLE PRECISION DTIME
#endif
C
      LOGICAL   C3LRSV, CC1ASV, CC1BSV, LDEBUG
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      CHARACTER*(*) FNDELD, FNCKJD, FNDKBC, FNTOC, FN3VI, FNDKBC3
      CHARACTER*(*) FN3FOP, FN3FOP2
      CHARACTER*(*) FN3FOPX, FN3FOP2X
      CHARACTER*5 FNDPTIA, FNDPTAB, FNDPTIJ, FNKAPAB, FNKAPIJ
      CHARACTER*6 FNDPTIA2
      CHARACTER*7 FNDIAJB, FNDAIJK, FNDIAJK
      CHARACTER*8 FNDABI1, FNDABI2, FNDABI3, FNDABI4
      CHARACTER*10 MODEL
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      CALL QENTER('CCSDPT_DENS2_SC')
C
      LDEBUG = .FALSE.
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
      IPRCC = IPRINT
      ISYMTR = MULD2H(ISYMT1,ISYMT2)
      ISYRES = MULD2H(ISYMTR,ISYMOP)
      ISINT1 = ISYMOP
      ISINT2 = MULD2H(ISYMT1,ISYMOP)
      ISYMIM = MULD2H(ISYMTR,ISYMOP)
C
C--------------------
C     Time variables.
C--------------------
C
      TITRAN = 0.0D0
      TISORT = 0.0D0
      TISMAT = 0.0D0
      TIQMAT = 0.0D0
      TIOME1 = 0.0D0
C
C--------------------------------------
C     Reorder the t2-amplitudes i T2TP.
C--------------------------------------
C
      IF (LWORK .LT. NT2SQ(ISYMT2)) THEN
         CALL QUIT('Not enough memory to construct T2TP (CCSDPT_DENS2)')
      ENDIF
C
      CALL DCOPY(NT2SQ(ISYMT2),T2TP,1,WORK,1)
      CALL CC3_T2TP(T2TP,WORK,ISYMT2)
C
      IF (IPRINT .GT. 55) THEN
         XT2TP = DDOT(NT2SQ(ISYMT2),T2TP,1,T2TP,1)
         WRITE(LUPRI,*) 'Norm of T2TP ',XT2TP
      ENDIF
C
C-----------------------------------------------
C     Reorder the l2-amplitudes i L2TP if CC3.
C-----------------------------------------------
C
      IF (CC3) THEN
C
         IF (LWORK .LT. NT2SQ(ISYML2)) THEN
            CALL QUIT('Not enough memory to construct L2TP')
         ENDIF
C
         CALL DCOPY(NT2SQ(ISYML2),L2TP,1,WORK,1)
         CALL CC3_T2TP(L2TP,WORK,ISYML2)
C
         IF (IPRINT .GT. 55) THEN
            XT2TP = DDOT(NT2SQ(ISYML2),L2TP,1,L2TP,1)
            WRITE(LUPRI,*) 'Norm of L2TP ',XT2TP
         ENDIF
C
      ENDIF
C
C---------------------------------------------------------
C     Read canonical orbital energies and MO coefficients.
C---------------------------------------------------------
C
      KFOCKD = 1
      KOMG1  = KFOCKD + NORBTS
      KOMG22 = KOMG1  + NT1AM(ISYMOP)
      KFCKBA = KOMG22 + NT2AM(ISYMOP)
      KEND0  = KFCKBA + N2BST(ISYMOP)
C
      IF (CC3) THEN
         KLAMDP = KEND0
         KLAMDH = KLAMDP + NLAMDT
         KEND0  = KLAMDH + NLAMDT
      ELSE
         KCMO = KEND0
         KEND0 = KCMO + NLAMDS
      ENDIF
C
      LWRK0  = LWORK  - KEND0
C
      IF (LWRK0 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND0
         CALL QUIT('Insufficient space in CCSDPT_DENS2')
      END IF
C
      CALL DZERO(WORK(KOMG1),NT1AM(ISYMOP))
      CALL DZERO(WORK(KOMG22),NT2AM(ISYMOP))
C
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
C
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KFOCKD+I-1), I=1,NORBTS)
C
      IF (.NOT. CC3) THEN
         READ (LUSIFC) (WORK(KCMO+I-1), I=1,NLAMDS)
         CALL CMO_REORDER(WORK(KCMO),WORK(KEND0),LWRK0)
      ENDIF
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
C----------------------------------------------
C     Calculate the lambda matrices for CC3
C----------------------------------------------
C
      IF (CC3) THEN
         CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),T1AM,
     *               WORK(KEND0),LWRK0)
C
         IF (IPRINT .GT.100) THEN
            CALL AROUND('Usual Lambda matrices ')
            CALL CC_PRLAM(WORK(KLAMDP),WORK(KLAMDH),1)
         ENDIF
      ENDIF
C
C-----------------------------------------------------
C     Construct the transformed Fock matrix
C-----------------------------------------------------
C
      LUFCK = -1
C
      IF (CC3) THEN
C     This AO Fock matrix is constructed from the T1 transformed density
         CALL GPOPEN(LUFCK,'CC_FCKH','UNKNOWN',' ','UNFORMATTED',
     *                 IDUMMY,.FALSE.)
      ELSE
C     This AO Fock matrix is constructed from the CMO transformed density
         CALL GPOPEN(LUFCK,'CC_FCKREF','UNKNOWN',' ','UNFORMATTED',
     *               IDUMMY,.FALSE.)
      ENDIF
C
      REWIND(LUFCK)
      READ(LUFCK)(WORK(KFCKBA + I-1),I = 1,N2BST(ISYMOP))
      CALL GPCLOSE(LUFCK,'KEEP' )
C
      IF (IPRINT .GT. 140) THEN
         CALL AROUND( 'Usual Fock AO matrix' )
         CALL CC_PRFCKAO(WORK(KFCKBA),ISYMOP)
      ENDIF
C
      ! SCF Fock matrix in transformed using CMO vector
      IF (CC3) THEN
         CALL CC_FCKMO(WORK(KFCKBA),WORK(KLAMDP),WORK(KLAMDH),
     *                 WORK(KEND0),LWRK0,1,1,1)
      ELSE
         CALL CC_FCKMO(WORK(KFCKBA),WORK(KCMO),WORK(KCMO),
     *                 WORK(KEND0),LWRK0,1,1,1)
      ENDIF
C
      IF (IPRINT .GT. 50) THEN
         CALL AROUND( 'In CC_ETA: Triples Fock MO matrix' )
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
         CALL AROUND('In CCSDPT_DENS2: Triples Fock MO matrix (sort)')
         CALL CC_PRFCKMO(WORK(KFCKBA),ISYMOP)
      ENDIF
C
C------------------------------------------------------------------
C     Read in another T2 amplitude, and transform it to 2*C-E
C     Square up to full matrix and reorder the index
C------------------------------------------------------------------
C
      IF (.NOT. CC3) THEN
         KT2TCME = KEND0
         KEND0   = KT2TCME + NT2SQ(1)
         LWRK0   = LWORK - KEND0
C
         IF (LWRK0 .LT. NT2SQ(1)) 
     *        CALL QUIT('Too litlle workspace CCSDPT_DENS2 T2TCME')
C
         IOPTT2 = 2
         CALL CC_RDRSP('R0',0,1,IOPTT2,MODEL,DUMMY,WORK(KEND0))
C
         ISYOPE = ISYMOP
         IOPTT2 = 1
         CALL CCSD_TCMEPK(WORK(KEND0),1.0D0,ISYOPE,IOPTT2)
C
         CALL CC_T2SQ(WORK(KEND0),WORK(KT2TCME),1)
C
         CALL DCOPY(NT2SQ(1),WORK(KT2TCME),1,WORK(KEND0),1)
         CALL CC3_T2TP(WORK(KT2TCME),WORK(KEND0),1)
C
         IF (IPRINT .GT. 55) THEN
            XT2TP = DDOT(NT2SQ(1),WORK(KT2TCME),1,WORK(KT2TCME),1)
            WRITE(LUPRI,*) 'Norm of 2*C-E T2 amplitudes after resort ',
     *                       XT2TP
         ENDIF
      ENDIF
C
C-----------------------------
C     Read occupied integrals.
C-----------------------------
C
C     Memory allocation.
C
      KTROC0 = KEND0
      KTROC02= KTROC0 + NTRAOC(ISINT2)
      KTROC2 = KTROC02+ NTRAOC(ISINT2)
      KTROC22= KTROC2 + NTRAOC(ISINT2)
      KXIAJB = KTROC22+ NTRAOC(ISINT2)
      KOCC1  = KXIAJB + NT2AM(ISYMOP)
      KOCC2  = KOCC1  + NCKIJ(ISYRES)
      KKAPAA = KOCC2  + NCKIJ(ISYRES)
      KKAPII = KKAPAA + NVIRT
      KEND1  = KKAPII + NRHFT
      LWRK1  = LWORK  - KEND1
C
      IF (CC3) THEN
         KTROC01 = KEND1
         KTROC21 = KTROC01 + NTRAOC(ISINT2)
         KTROC03 = KTROC21 + NTRAOC(ISINT2)
         KTROC23 = KTROC03 + NTRAOC(ISINT2)
         KEND1   = KTROC23 + NTRAOC(ISINT2)
         LWRK1   = LWORK  - KEND1
      ENDIF
C
      IF (.NOT. RELORB) THEN
         KOMG12  = KEND1
         KDENSAB = KOMG12  + NT1AM(ISYRES)
         KDENSIJ = KDENSAB + NMATAB(ISYRES)
         KEND1   = KDENSIJ + NMATIJ(ISYRES)
         LWRK1   = LWORK - KEND1
C
         CALL DZERO(WORK(KOMG12),NT1AM(ISYRES))
         CALL DZERO(WORK(KDENSAB),NMATAB(ISYRES))
         CALL DZERO(WORK(KDENSIJ),NMATIJ(ISYRES))
      ENDIF
C
      KINTOC = KEND1
      KEND2  = KINTOC + MAX(NTOTOC(ISYMOP),NTOTOC(ISINT2))
      LWRK2  = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in CCSDPT_DENS2')
      END IF
C
C
C----------------------------------
C     Initialize result vectors
C----------------------------------
C
      CALL DZERO(WORK(KOCC1),NCKIJ(ISYRES))
      CALL DZERO(WORK(KOCC2),NCKIJ(ISYRES))
      CALL DZERO(WORK(KKAPAA),NVIRT)
      CALL DZERO(WORK(KKAPII),NRHFT)
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
         XIAJB = DDOT(NT2AM(ISYMOP),WORK(KXIAJB),1,
     *                WORK(KXIAJB),1)
         WRITE(LUPRI,*) 'Norm of IAJB ',XIAJB
      ENDIF
C
C------------------------
C     Occupied integrals.
C------------------------
C
      IF (CC3) THEN
         IOFF = 1
         IF (NTOTOC(ISYMOP) .GT. 0) THEN
            CALL GETWA2(LUCKJD,FNCKJD,WORK(KINTOC),IOFF,NTOTOC(ISYMOP))
         ENDIF
      ELSE
         IOFF = 1
         IF (NTOTOC(ISYMOP) .GT. 0) THEN
            CALL GETWA2(LUTOC,FNTOC,WORK(KINTOC),IOFF,NTOTOC(ISYMOP))
         ENDIF
      ENDIF
C
C----------------------------------
C     Write out norms of Integrals.
C----------------------------------
C
      IF (IPRINT .GT. 55) THEN
         XINT  = DDOT(NTOTOC(ISYMOP),WORK(KINTOC),1,
     *                WORK(KINTOC),1)
         WRITE(LUPRI,*) 'Norm of OCC-INT ',XINT
      ENDIF
C
C----------------------------------------------------------------------
C     Transform (ia|j delta) integrals to (ia|j k) and sort as (ij,k,a)
C----------------------------------------------------------------------
C
      DTIME = SECOND()
      IF (CC3) THEN
         CALL CC3_TROCC(WORK(KINTOC),WORK(KTROC0),WORK(KLAMDP),
     *                    WORK(KEND2),LWRK2,ISINT2)
      ELSE
         CALL CCSDT_TROCC(WORK(KINTOC),WORK(KTROC0),WORK(KCMO),
     *                    WORK(KEND2),LWRK2)
      ENDIF
C
      DTIME  = SECOND() - DTIME
      TITRAN = TITRAN   + DTIME

C
      DTIME = SECOND()
C
      DTIME  = SECOND() - DTIME
      TISORT = TISORT   + DTIME
C
C-----------------------------------------------------------
C     Construct 2*C-E of the integrals.
C     Have integral for both (ij,k,a) and (a,k,j,i)
C-----------------------------------------------------------
C
      CALL CCSDT_TCMEOCC(WORK(KTROC0),WORK(KTROC2),ISINT2)
C
      IF (CC3) THEN
         IOFF = 1
         IF (NTOTOC(ISINT2) .GT. 0) THEN
            CALL GETWA2(LUTOC,FNTOC,WORK(KINTOC),IOFF,NTOTOC(ISINT2))
         ENDIF
C
         CALL CC3_TROCC(WORK(KINTOC),WORK(KTROC01),WORK(KLAMDH),
     *                    WORK(KEND2),LWRK2,ISINT2)
C
         CALL CCSDT_TCMEOCC(WORK(KTROC01),WORK(KTROC21),ISINT2)
C
         CALL CCFOP_SORT(WORK(KTROC01),WORK(KTROC03),ISINT2,1)
C
         CALL CCFOP_SORT(WORK(KTROC21),WORK(KTROC23),ISINT2,1)
      ENDIF
C
      CALL CCFOP_SORT(WORK(KTROC0),WORK(KTROC02),ISINT2,1)
C
      CALL CCFOP_SORT(WORK(KTROC2),WORK(KTROC22),ISINT2,1)
C
C-------------------------------
C     Write out norms of arrays.
C-------------------------------
C
      IF (IPRINT .GT. 55) THEN
         XINT  = DDOT(NTOTOC(ISINT2),WORK(KINTOC),1,
     *                WORK(KINTOC),1)
         WRITE(LUPRI,*) 'Norm of CKJDEL-INT  ',XINT
      ENDIF
C
      IF (IPRINT .GT. 55) THEN
         XTROC0 = DDOT(NTRAOC(ISINT2),WORK(KTROC0),1,
     *                WORK(KTROC0),1)
         WRITE(LUPRI,*) 'Norm of TROC0 ',XTROC0
      ENDIF
C
      IF (IPRINT .GT. 55) THEN
         XTROC0 = DDOT(NTRAOC(ISINT2),WORK(KTROC2),1,
     *                WORK(KTROC2),1)
         WRITE(LUPRI,*) 'Norm of TROC2 ',XTROC0
      ENDIF
C
C--------------------------------------------------------
C     Open files to the one and two electron densities.
C--------------------------------------------------------
C
      LUPTIA   = -1
      FNDPTIA  = 'DPTIA'
C     d_{ia}
      CALL WOPEN2(LUPTIA,FNDPTIA,64,0)
C
      IF ((.NOT. CC3) .AND. (RELORB)) THEN
         LUPTIAJB = -1
         LUABI1   = -1
         LUABI2   = -1
         LUABI3   = -1
         LUABI4   = -1
         LUAIJK   = -1
         LUIAJK   = -1
         FNDIAJB  = 'DPTIAJB'
         FNDABI1  = 'DPTABIC1'
         FNDABI2  = 'DPTABIC2'
         FNDABI3  = 'DPTABCI1'
         FNDABI4  = 'DPTABCI2'
         FNDAIJK  = 'DPTAIJK'
         FNDIAJK  = 'DPTIAJK'
C
C        d_{iajb}
         CALL WOPEN2(LUPTIAJB,FNDIAJB,64,0)
C        d_{abic_1}
         CALL WOPEN2(LUABI1,FNDABI1,64,0)
C        d_{abic_2}
         CALL WOPEN2(LUABI2,FNDABI2,64,0)
C        d_{abci_1}
         CALL WOPEN2(LUABI3,FNDABI3,64,0)
C        d_{abci_2}
         CALL WOPEN2(LUABI4,FNDABI4,64,0)
C        d_{aijk}
         CALL WOPEN2(LUAIJK,FNDAIJK,64,0)
C        d_{iajk}
         CALL WOPEN2(LUIAJK,FNDIAJK,64,0)
      ELSE
         LUPTIA2  = -1
         LUPTAB   = -1
         LUPTIJ   = -1
         FNDPTIA2 = 'DPTIA2'
         FNDPTAB  = 'DPTAB'
         FNDPTIJ  = 'DPTIJ'
C        d_{ia}
         CALL WOPEN2(LUPTIA2,FNDPTIA2,64,0)
C        d_{ab}
         CALL WOPEN2(LUPTAB,FNDPTAB,64,0)
C        d_{ij}
         CALL WOPEN2(LUPTIJ,FNDPTIJ,64,0)
      ENDIF
C
C-----------------------------------------
!     Sonia: old code!!!!
C     Save workspace adress to regain workspace
C     in the end.
C-----------------------------------------
C
      KENDSV = KEND1
C
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
C
         IF (IPRINT .GT. 55) THEN
C
            WRITE(LUPRI,*) 'In CCSDPT_DENS2: ISAIJ1 :',ISAIJ1
            WRITE(LUPRI,*) 'In CCSDPT_DENS2: ISYCKB :',ISYCKB
            WRITE(LUPRI,*) 'In CCSDPT_DENS2: ISCKB1 :',ISCKB1
            WRITE(LUPRI,*) 'In CCSDPT_DENS2: ISCKB2 :',ISCKB2
C
         ENDIF
C
C--------------------------
C        Memory allocation.
C--------------------------
C
         KTRVI1 = KEND1
         KTRVI2 = KTRVI1 + NCKATR(ISCKB2)
         KRMAT1 = KTRVI2 + NCKATR(ISCKB2)
         KEND2  = KRMAT1 + NCKI(ISAIJ1)
         LWRK2  = LWORK  - KEND2
C
         KTRVI0  = KEND2
         KTRVI3  = KTRVI0  + NCKATR(ISCKB2)
         KTRVI4  = KTRVI3  + NCKATR(ISCKB2)
         KTRVI5  = KTRVI4  + NCKATR(ISCKB2)
         KTRVI6  = KTRVI5  + NCKATR(ISCKB2)
         KTRVI7  = KTRVI6  + NCKATR(ISCKB2)
         KVIR1   = KTRVI7  + NCKATR(ISCKB2)
         KVIR2   = KVIR1   + NCKATR(ISAIJ1)
         KVIR3   = KVIR2   + NCKATR(ISAIJ1)
         KVIR4   = KVIR3   + NCKATR(ISAIJ1)
         KEND3   = KVIR4   + NCKATR(ISAIJ1)
         LWRK3   = LWORK  - KEND3
C
         IF (CC3) THEN
            KTRVI14 = KEND3
            KTRVI15 = KTRVI14 + NCKATR(ISCKB2)
            KTRVI18 = KTRVI15 + NCKATR(ISCKB2)
            KTRVI19 = KTRVI18 + NCKATR(ISCKB2)
            KEND3   = KTRVI19 + NCKATR(ISCKB2)
            LWRK3   = LWORK  - KEND3
         ENDIF
C
         KINTVI = KEND3
         KEND4  = KINTVI + MAX(NCKA(ISYMD),NCKA(ISCKB2))
         LWRK4  = LWORK  - KEND4
C
         IF (LWRK4 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND4
            CALL QUIT('Insufficient space in CCSDPT_DENS2')
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
            CALL DZERO(WORK(KVIR1),NCKATR(ISAIJ1))
            CALL DZERO(WORK(KVIR2),NCKATR(ISAIJ1))
            CALL DZERO(WORK(KVIR3),NCKATR(ISAIJ1))
            CALL DZERO(WORK(KVIR4),NCKATR(ISAIJ1))
C
C-----------------------------------------------
C           Integrals used in s3am.
C-----------------------------------------------
C
            IF (CC3) THEN
               IOFF = ICKBD(ISCKB2,ISYMD) + NCKATR(ISCKB2)*(D - 1) + 1
               IF (NCKATR(ISCKB2) .GT. 0) THEN
                  CALL GETWA2(LUDKBC,FNDKBC,WORK(KTRVI0),IOFF,
     &                        NCKATR(ISCKB2))
               ENDIF
            ELSE
               IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D - 1) + 1
               IF (NCKA(ISYCKB) .GT. 0) THEN
                  CALL GETWA2(LUDKBC3,FNDKBC3,WORK(KINTVI),IOFF,
     &                        NCKA(ISYCKB))
               ENDIF
C
               CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI0),WORK(KCMO),
     *                          ISYMD,D,ISINT2,WORK(KEND4),LWRK4)
            ENDIF
C
C------------------------------------------------------
C           Read 2*C-E of integral used for t3-bar
C------------------------------------------------------
C
            IF (CC3) THEN
               IOFF = ICKBD(ISCKB2,ISYMD) + NCKATR(ISCKB2)*(D - 1) + 1
               IF (NCKATR(ISCKB2) .GT. 0) THEN
                  CALL GETWA2(LU3FOP2X,FN3FOP2X,WORK(KTRVI4),IOFF,
     &                        NCKATR(ISCKB2))
               ENDIF
            ELSE
               IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D - 1) + 1
               IF (NCKA(ISYCKB) .GT. 0) THEN
                  CALL GETWA2(LU3FOP2,FN3FOP2,WORK(KINTVI),IOFF,
     *                        NCKA(ISYCKB))
               ENDIF
C
               CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI4),WORK(KCMO),
     *                          ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
            ENDIF
C
C------------------------------------------------------------
C           Integrals used for t3-bar for cc3
C------------------------------------------------------------
C
            IF (CC3) THEN
               IOFF = ICKBD(ISCKB2,ISYMD) + NCKATR(ISCKB2)*(D - 1) + 1
               IF (NCKATR(ISCKB2) .GT. 0) THEN
                  CALL GETWA2(LUDKBC3,FNDKBC3,WORK(KTRVI14),IOFF,
     &                        NCKATR(ISCKB2))
               ENDIF
               CALL CCSDT_SRVIR3(WORK(KTRVI14),WORK(KEND4),
     *                           ISYMD,D,ISINT2)
               CALL CCSDT_SRTVIR(WORK(KTRVI14),WORK(KTRVI15),WORK(KEND4)
     *                           ,LWRK4,ISYMD,ISINT2)
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
               XTRVI0= DDOT(NCKATR(ISCKB2),WORK(KTRVI0),1,
     *                      WORK(KTRVI0),1)
               WRITE(LUPRI,*) 'Norm of TRVI0 ',XTRVI0
            ENDIF
C
            IF (IPRINT .GT. 55) THEN
               XTRVI2= DDOT(NCKATR(ISCKB2),WORK(KTRVI2),1,
     *                      WORK(KTRVI2),1)
               WRITE(LUPRI,*) 'Norm of TRVI2 ',XTRVI2
            ENDIF
C
            IF (IPRINT .GT. 55) THEN
               XTRVI0= DDOT(NCKATR(ISCKB2),WORK(KTRVI4),1,
     *                      WORK(KTRVI4),1)
               WRITE(LUPRI,*) 'Norm of TRVI4 ',XTRVI0
            ENDIF
C
            IF (IPRINT .GT. 55) THEN
               XTRVI2= DDOT(NCKATR(ISCKB2),WORK(KTRVI5),1,
     *                      WORK(KTRVI5),1)
               WRITE(LUPRI,*) 'Norm of TRVI5 ',XTRVI2
            ENDIF
C
C------------------------------------------------------
C           Read virtual integrals used in contraction.
C------------------------------------------------------
C
            IF (CC3) THEN
               IOFF = ICKAD(ISCKB2,ISYMD) + NCKA(ISCKB2)*(D - 1) + 1
               IF (NCKA(ISCKB2) .GT. 0) THEN
                  CALL GETWA2(LUDELD,FNDELD,WORK(KINTVI),IOFF,
     *                        NCKA(ISCKB2))
               ENDIF
C
               CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI1),WORK(KLAMDH),
     *                          ISYMD,D,ISINT2,WORK(KEND4),LWRK4)
C
            ELSE
               IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D - 1) + 1
               IF (NCKA(ISYCKB) .GT. 0) THEN
                  CALL GETWA2(LU3VI,FN3VI,WORK(KINTVI),IOFF,
     &                        NCKA(ISYCKB))
               ENDIF
C
               CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI1),WORK(KCMO),
     *                          ISYMD,D,ISINT2,WORK(KEND4),LWRK4)
            ENDIF
C
C--------------------------------------------------------
C           Calculate virtual integrals used in q3am.
C--------------------------------------------------------
C
            CALL DCOPY(NCKATR(ISCKB2),WORK(KTRVI1),1,WORK(KTRVI3),1)
C
            IF (LWRK3 .LT. NCKATR(ISYCKB)) THEN
               CALL QUIT('Insufficient space for allocation in '//
     &                   'CCSDPT_DENS2 (1)')
            END IF
C
            DTIME = SECOND()
            CALL CCSDT_SRVIR3(WORK(KTRVI3),WORK(KEND4),ISYMD,D,ISINT2)
C
            DTIME  = SECOND() - DTIME
            TISORT = TISORT   + DTIME
C
            IF (IPRINT .GT. 55) THEN
               XTRVI3= DDOT(NCKATR(ISCKB2),WORK(KTRVI3),1,
     *                      WORK(KTRVI3),1)
               WRITE(LUPRI,*) 'Norm of TRVI3 ',XTRVI3
            ENDIF
C
C---------------------------------------------------------------
C           Read virtual integrals used in q3am/u3am for t3-bar.
C---------------------------------------------------------------
C
            IF (CC3) THEN
               IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D - 1) + 1
               IF (NCKA(ISYCKB) .GT. 0) THEN
                  CALL GETWA2(LU3VI,FN3VI,WORK(KINTVI),IOFF,
     *                        NCKA(ISYCKB))
               ENDIF
C
               CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI19),WORK(KLAMDP),
     *                          ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
C
               IF (LWRK4 .LT. NCKATR(ISYCKB)) THEN
                  CALL QUIT('Insufficient space for allocation in '//
     *                      'CCSDPT_DENS2  (CC3 TRVI)')
               END IF
C
               CALL CCSDT_SRTVIR(WORK(KTRVI19),WORK(KTRVI18),WORK(KEND4)
     *                           ,LWRK4,ISYMD,ISINT2)
            ENDIF
C
            IF (CC3) THEN
               IOFF = ICKBD(ISYCKB,ISYMD) + NCKATR(ISYCKB)*(D - 1) + 1
               IF (NCKATR(ISYCKB) .GT. 0) THEN
                  CALL GETWA2(LU3FOPX,FN3FOPX,WORK(KTRVI6),IOFF,
     *                        NCKATR(ISYCKB))
               ENDIF
            ELSE
               IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D - 1) + 1
               IF (NCKA(ISYCKB) .GT. 0) THEN
                  CALL GETWA2(LU3FOP,FN3FOP,WORK(KINTVI),IOFF,
     *                        NCKA(ISYCKB))
               ENDIF
C
               CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI6),WORK(KCMO),
     *                          ISYMD,D,ISINT2,WORK(KEND4),LWRK4)
            ENDIF
C
            IF (LWRK3 .LT. NCKATR(ISYCKB)) THEN
               CALL QUIT('Insufficient space for allocation in '//
     &                   'CCSDPT_DENS2 (2)')
            END IF
C
            CALL DCOPY(NCKATR(ISCKB2),WORK(KTRVI6),1,WORK(KTRVI7),1)
C
            DTIME = SECOND()
            CALL CCSDT_SRVIR3(WORK(KTRVI6),WORK(KEND4),ISYMD,D,ISINT2)
C
            DTIME  = SECOND() - DTIME
            TISORT = TISORT   + DTIME
C
            IF (IPRINT .GT. 55) THEN
               XTRVI3= DDOT(NCKATR(ISCKB2),WORK(KTRVI6),1,
     *                      WORK(KTRVI6),1)
               WRITE(LUPRI,*) 'Norm of TRVI6 ',XTRVI3
               XTRVI3= DDOT(NCKATR(ISCKB2),WORK(KTRVI7),1,
     *                      WORK(KTRVI7),1)
               WRITE(LUPRI,*) 'Norm of TRVI7 ',XTRVI3
            ENDIF
C
            IF (IPRINT .GT. 55) THEN
               XTRVI1= DDOT(NCKATR(ISCKB2),WORK(KTRVI1),1,
     *                      WORK(KTRVI1),1)
               WRITE(LUPRI,*) 'Norm of TRVI1 ',XTRVI1
            ENDIF
C
C---------------------
C           Calculate.
C---------------------
C
            DO ISYMB = 1,NSYM
C
               ISYALJ  = MULD2H(ISYMB,ISYMT2)
               ISYALJ2 = MULD2H(ISYMD,ISYMT2)
               ISAIJ2  = MULD2H(ISYMB,ISYRES)
               ISYMBD  = MULD2H(ISYMB,ISYMD)
               ISCKIJ  = MULD2H(ISYMBD,ISYMIM)
               ISYCKD  = MULD2H(ISYMOP,ISYMB)
               ISCKD2  = MULD2H(ISINT2,ISYMB)
C
               IF ((IPRINT .GT. 55)) THEN
C
                  WRITE(LUPRI,*) 'In CCSDPT_DENS2: ISYMD :',ISYMD
                  WRITE(LUPRI,*) 'In CCSDPT_DENS2: ISYMB :',ISYMB
                  WRITE(LUPRI,*) 'In CCSDPT_DENS2: ISYALJ:',ISYALJ
                  WRITE(LUPRI,*) 'In CCSDPT_DENS2: ISAIJ2:',ISAIJ2
                  WRITE(LUPRI,*) 'In CCSDPT_DENS2: ISYMBD:',ISYMBD
                  WRITE(LUPRI,*) 'In CCSDPT_DENS2: ISCKIJ:',ISCKIJ
C
               ENDIF
C
C              Can use kend3 since we do not need the integrals anymore.
               KSMAT   = KEND3
               KQMAT   = KSMAT   + NCKIJ(ISCKIJ)
               KSMAT2  = KQMAT   + NCKIJ(ISCKIJ)
               KSMAT3  = KSMAT2  + NCKIJ(ISCKIJ)
               KQMAT2  = KSMAT3  + NCKIJ(ISCKIJ)
               KUMAT   = KQMAT2  + NCKIJ(ISCKIJ)
               KUMAT2  = KUMAT   + NCKIJ(ISCKIJ)
               KUMAT3  = KUMAT2  + NCKIJ(ISCKIJ)
               KDIAG   = KUMAT3  + NCKIJ(ISCKIJ)
               KINDSQ  = KDIAG   + NCKIJ(ISCKIJ)
               KINDEX  = KINDSQ  + (6*NCKIJ(ISCKIJ) - 1)/IRAT + 1
               KINDEX2 = KINDEX  + (NCKI(ISYALJ) - 1)/IRAT + 1
               KTMAT   = KINDEX2 + (NCKI(ISYALJ2) - 1)/IRAT + 1
               KRMAT2  = KTMAT   + NCKIJ(ISCKIJ)
               KTRVI8  = KRMAT2  + NCKI(ISAIJ2)
               KTRVI9  = KTRVI8  + NCKATR(ISCKD2)
               KTRVI10 = KTRVI9  + NCKATR(ISCKD2)
               KEND4   = KTRVI10 + NCKATR(ISCKD2)
               LWRK4   = LWORK   - KEND4
C
               IF (CC3) THEN
                  KTRVI16 = KEND4
                  KTRVI17 = KTRVI16 + NCKATR(ISCKD2)
                  KTRVI20 = KTRVI17 + NCKATR(ISCKD2)
                  KEND4   = KTRVI20 + NCKATR(ISCKD2)
                  LWRK4   = LWORK  - KEND4
               ENDIF
C
               IF (.NOT. RELORB) THEN
                  KSMAT4  = KEND4
                  KUMAT4  = KSMAT4 + NCKIJ(ISCKIJ)
                  KTRVI11 = KUMAT4 + NCKIJ(ISCKIJ)
                  KTRVI12 = KTRVI11 + NCKATR(ISCKD2)
                  KTRVI13 = KTRVI12 + NCKATR(ISCKD2)
                  KEND4   = KTRVI13 + NCKATR(ISCKD2)
                  LWRK4   = LWORK-KEND4
               ENDIF
C
               KINTVI  = KEND4
COMMENT COMMENT
C               KEND5   = KINTVI  + NCKA(ISCKD2)
               KEND5   = KINTVI  + MAX(NCKA(ISYMB),NCKA(ISCKD2))
COMMENT COMMENT
               LWRK5   = LWORK   - KEND5
C
               IF (LWRK5 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Memory available : ',LWORK
                  WRITE(LUPRI,*) 'Memory needed    : ',KEND5
                  CALL QUIT('Insufficient space in CCSDPT_DENS2')
               END IF
C
C---------------------------------------------
C              Construct part of the diagonal.
C---------------------------------------------
C
               CALL CC3_DIAG(WORK(KDIAG),WORK(KFOCKD),ISCKIJ)
C
               IF ((IPRINT .GT. 55)) THEN
                  XDIA  = DDOT(NCKIJ(ISCKIJ),WORK(KDIAG),1,
     *                    WORK(KDIAG),1)
                  WRITE(LUPRI,*) 'Norm of DIA  ',XDIA
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
C
               DO B = 1,NVIR(ISYMB)
C
C-----------------------------------------
C                 Initialize the R2 matrix.
C-----------------------------------------
C
                  CALL DZERO(WORK(KRMAT2),NCKI(ISAIJ2))
C
C-------------------------------------------------------------
C           Read and transform integrals used in second S
C-------------------------------------------------------------
C
                  IF (CC3) THEN
                     IOFF = ICKBD(ISYCKD,ISYMB) 
     *                    + NCKATR(ISYCKD)*(B - 1) + 1
                     IF (NCKATR(ISYCKD) .GT. 0) THEN
                        CALL GETWA2(LUDKBC,FNDKBC,WORK(KTRVI8),IOFF,
     *                              NCKATR(ISYCKD))
                     ENDIF
                  ELSE
C
                     IOFF = ICKAD(ISYCKD,ISYMB) 
     *                    + NCKA(ISYCKD)*(B - 1) + 1
                     IF (NCKA(ISYCKD) .GT. 0) THEN
                        CALL GETWA2(LUDKBC3,FNDKBC3,WORK(KINTVI),IOFF,
     *                             NCKA(ISYCKD))
                     ENDIF
C
                     CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI8),
     *                                WORK(KCMO),ISYMB,B,ISINT2,
     *                                WORK(KEND5),LWRK5)
                  ENDIF
C
                  CALL CCSDT_SRTVIR(WORK(KTRVI8),WORK(KTRVI9),
     *                              WORK(KEND4),LWRK4,ISYMB,ISINT2)
C
                  IF (CC3) THEN
                     IOFF = ICKBD(ISYCKD,ISYMB) 
     *                    + NCKATR(ISYCKD)*(B - 1) + 1
                     IF (NCKATR(ISYCKD) .GT. 0) THEN
                        CALL GETWA2(LUDKBC3,FNDKBC3,WORK(KTRVI16),IOFF,
     *                              NCKATR(ISYCKD))
                     ENDIF
                     CALL CCSDT_SRVIR3(WORK(KTRVI16),WORK(KEND5),
     *                                 ISYMB,B,ISINT2)
                     CALL CCSDT_SRTVIR(WORK(KTRVI16),WORK(KTRVI17),
     *                                 WORK(KEND4),LWRK4,ISYMB,ISINT2)
                  ENDIF
C
C----------------------------------------------------------
C           Read virtual integrals used in second U
C----------------------------------------------------------
C
C
                  IF (CC3) THEN
                     IOFF = ICKAD(ISCKD2,ISYMB) 
     *                    + NCKA(ISCKD2)*(B - 1) + 1
                     IF (NCKA(ISYCKD) .GT. 0) THEN
                        CALL GETWA2(LUDELD,FNDELD,WORK(KINTVI),IOFF,
     *                              NCKA(ISCKD2))
                     ENDIF
C
                     CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI10),
     *                                WORK(KLAMDH),ISYMB,B,ISINT2,
     *                                WORK(KEND5),LWRK5)
C
                  ELSE
C
                     IOFF = ICKAD(ISYCKD,ISYMB) 
     *                    + NCKA(ISYCKD)*(B - 1) + 1
                     IF (NCKA(ISYCKD) .GT. 0) THEN
                        CALL GETWA2(LU3VI,FN3VI,WORK(KINTVI),IOFF,
     *                              NCKA(ISYCKD))
                     ENDIF
C
                     CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI10),
     *                                WORK(KCMO),ISYMB,B,ISYMOP,
     *                                WORK(KEND5),LWRK5)
                  ENDIF
C
C------------------------------------------------------------------------
C           Read and transform integrals used in second S-bar and U-bar
C           NOT used for CC3
C------------------------------------------------------------------------
C
                  IF (.NOT. RELORB) THEN
C
                     IF (CC3) THEN
                        IOFF = ICKBD(ISYCKD,ISYMB) 
     *                       + NCKATR(ISYCKD)*(B-1) + 1
                        IF (NCKATR(ISYCKD) .GT. 0) THEN
                           CALL GETWA2(LU3FOP2X,FN3FOP2X,WORK(KTRVI11),
     *                                 IOFF,NCKATR(ISYCKD))
                        ENDIF
                     ELSE
                        IOFF = ICKAD(ISYCKD,ISYMB) 
     *                       + NCKA(ISYCKD)*(B-1) + 1
                        IF (NCKA(ISYCKD) .GT. 0) THEN
                           CALL GETWA2(LU3FOP2,FN3FOP2,WORK(KINTVI),
     *                                 IOFF,NCKA(ISYCKD))
                        ENDIF
C
                        CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI11),
     *                                   WORK(KCMO),ISYMB,B,ISYMOP,
     *                                   WORK(KEND5),LWRK5)
C
                     ENDIF
C
                     CALL CCSDT_SRTVIR(WORK(KTRVI11),WORK(KTRVI12),
     *                                 WORK(KEND5),LWRK5,ISYMB,
     *                                 ISINT2)
C
                     IF (CC3) THEN
                        IOFF = ICKBD(ISYCKD,ISYMB) 
     *                       + NCKATR(ISYCKD)*(B - 1) + 1
                        IF (NCKATR(ISYCKD) .GT. 0) THEN
                           CALL GETWA2(LU3FOPX,FN3FOPX,WORK(KTRVI13),
     *                                 IOFF,NCKATR(ISYCKD))
                        ENDIF
C
                        IOFF = ICKAD(ISYCKD,ISYMB) 
     *                       + NCKA(ISYCKD)*(B - 1) + 1
                        IF (NCKA(ISYCKD) .GT. 0) THEN
                           CALL GETWA2(LU3VI,FN3VI,WORK(KINTVI),IOFF,
     *                                 NCKA(ISYCKD))
                        ENDIF
C
                        CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI20),
     *                                   WORK(KLAMDP),ISYMB,B,ISYMOP,
     *                                   WORK(KEND4),LWRK4)
                     ELSE
                        IOFF = ICKAD(ISYCKD,ISYMB) 
     *                       + NCKA(ISYCKD)*(B-1) + 1
                        IF (NCKA(ISYCKD) .GT. 0) THEN
                           CALL GETWA2(LU3FOP,FN3FOP,WORK(KINTVI),IOFF,
     *                                 NCKA(ISYCKD))
                        ENDIF
C
                        CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI13),
     *                                   WORK(KCMO),ISYMB,B,ISINT2,
     *                                   WORK(KEND5),LWRK5)
                     ENDIF
                  ENDIF
C
C-------------------------------------------------------------------
C                 Calculate the S(ci,bk,dj) matrix for T3 for B,D.
C-------------------------------------------------------------------
C
                  DTIME = SECOND()
                  CALL CC3_SMAT(0.0D0,T2TP,ISYMT2,WORK(KTMAT),
     *                          WORK(KTRVI0),
     *                          WORK(KTRVI2),WORK(KTROC0),ISINT2,
     *                          WORK(KFOCKD),WORK(KDIAG),
     *                          WORK(KSMAT),WORK(KEND4),LWRK4,
     *                          WORK(KINDEX),WORK(KINDSQ),LENSQ,
     *                          ISYMB,B,ISYMD,D)
C
                  CALL T3_FORBIDDEN(WORK(KSMAT),ISYMIM,ISYMB,B,ISYMD,D)
C
                  DTIME  = SECOND() - DTIME
                  TISMAT = TISMAT   + DTIME
C
                  IF (IPRINT .GT. 55) THEN
                     XSMAT = DDOT(NCKIJ(ISCKIJ),WORK(KSMAT),1,
     *                       WORK(KSMAT),1)
                     WRITE(LUPRI,*) 'Norm of SMAT     ',XSMAT
                  ENDIF
C
C-------------------------------------------------------------------
C                 Calculate the S(ci,bk,dj) matrix for T3 for D,B.
C-------------------------------------------------------------------
C
                  DTIME = SECOND()
                  CALL CC3_SMAT(0.0D0,T2TP,ISYMT2,WORK(KTMAT),
     *                          WORK(KTRVI8),
     *                          WORK(KTRVI9),WORK(KTROC0),ISINT2,
     *                          WORK(KFOCKD),WORK(KDIAG),
     *                          WORK(KSMAT3),WORK(KEND4),LWRK4,
     *                          WORK(KINDEX2),WORK(KINDSQ),LENSQ,
     *                          ISYMD,D,ISYMB,B)
C
                  CALL T3_FORBIDDEN(WORK(KSMAT3),ISYMIM,ISYMD,D,ISYMB,B)
C
                  DTIME  = SECOND() - DTIME
                  TISMAT = TISMAT   + DTIME
C
                  IF (IPRINT .GT. 55) THEN
                     XSMAT = DDOT(NCKIJ(ISCKIJ),WORK(KSMAT3),1,
     *                       WORK(KSMAT3),1)
                     WRITE(LUPRI,*) 'Norm of SMAT3    ',XSMAT
                  ENDIF
C
C---------------------------------------------------------------------------
C                 Calculate the S(ci,bk,dj) matrix for for B,D for T3-BAR.
C---------------------------------------------------------------------------
C
                  DTIME = SECOND()
C
                  CALL DZERO(WORK(KSMAT2),NCKIJ(ISCKIJ))
C
                  IF (CC3) THEN
                     CALL CCFOP_SMAT(0.0D0,L1AM,ISYML1,L2TP,ISYML2,
     *                               WORK(KTMAT),
     *                               WORK(KFCKBA),WORK(KXIAJB),ISINT1,
     *                               WORK(KTRVI14),WORK(KTRVI15),
     *                               WORK(KTRVI4),WORK(KTRVI5),
     *                               WORK(KTROC01),WORK(KTROC21),
     *                               ISINT2,WORK(KFOCKD),WORK(KDIAG),
     *                               WORK(KSMAT2),WORK(KEND4),LWRK4,
     *                               WORK(KINDEX),WORK(KINDSQ),LENSQ,
     *                               ISYMB,B,ISYMD,D)
C
                     CALL DSCAL(NCKIJ(ISCKIJ),HALF,WORK(KSMAT2),1)
C
                  ELSE
C
                     CALL CCFOP_SMAT(0.0D0,T1AM,ISYMT1,WORK(KT2TCME),
     *                               ISYMT2,WORK(KTMAT),
     *                               WORK(KFCKBA),WORK(KXIAJB),ISINT1,
     *                               WORK(KTRVI0),WORK(KTRVI2),
     *                               WORK(KTRVI4),WORK(KTRVI5),
     *                               WORK(KTROC0),WORK(KTROC2),
     *                               ISINT2,WORK(KFOCKD),WORK(KDIAG),
     *                               WORK(KSMAT2),WORK(KEND4),LWRK4,
     *                               WORK(KINDEX),WORK(KINDSQ),LENSQ,
     *                               ISYMB,B,ISYMD,D)
                  ENDIF
C
                  CALL T3_FORBIDDEN(WORK(KSMAT2),ISYMIM,ISYMB,B,ISYMD,D)
C
                  DTIME  = SECOND() - DTIME
                  TISMAT = TISMAT   + DTIME
C
                  IF (IPRINT .GT. 55) THEN
                     XSMAT = DDOT(NCKIJ(ISCKIJ),WORK(KSMAT2),1,
     *                       WORK(KSMAT2),1)
                     WRITE(LUPRI,*) 'Norm of SMAT-BAR-1 ',XSMAT
                  ENDIF
C
C-----------------------------------------------------------------------------
C                 Calculate the S(ci,bk,dj) matrix for for D,B for T3-BAR.
C-----------------------------------------------------------------------------
C
                  IF (.NOT. RELORB) THEN
C
                     DTIME = SECOND()
C
                     CALL DZERO(WORK(KSMAT4),NCKIJ(ISCKIJ))
C
                     IF (CC3) THEN
                        CALL CCFOP_SMAT(0.0D0,L1AM,ISYML1,L2TP,
     *                                  ISYML2,WORK(KTMAT),WORK(KFCKBA),
     *                                  WORK(KXIAJB),ISINT1,
     *                                  WORK(KTRVI16),WORK(KTRVI17),
     *                                  WORK(KTRVI11),WORK(KTRVI12),
     *                                  WORK(KTROC01),WORK(KTROC21),
     *                                  ISINT2,WORK(KFOCKD),WORK(KDIAG),
     *                                  WORK(KSMAT4),WORK(KEND4),LWRK4,
     *                                  WORK(KINDEX2),WORK(KINDSQ),
     *                                  LENSQ,ISYMD,D,ISYMB,B)
C
                        CALL DSCAL(NCKIJ(ISCKIJ),HALF,WORK(KSMAT4),1)
C
                     ELSE
C
                        CALL CCFOP_SMAT(0.0D0,T1AM,ISYMT1,WORK(KT2TCME),
     *                                  ISYMT2,WORK(KTMAT),WORK(KFCKBA),
     *                                  WORK(KXIAJB),ISINT1,
     *                                  WORK(KTRVI8),WORK(KTRVI9),
     *                                  WORK(KTRVI11),WORK(KTRVI12),
     *                                  WORK(KTROC0),WORK(KTROC2),
     *                                  ISINT2,WORK(KFOCKD),WORK(KDIAG),
     *                                  WORK(KSMAT4),WORK(KEND4),LWRK4,
     *                                  WORK(KINDEX2),WORK(KINDSQ),
     *                                  LENSQ,ISYMD,D,ISYMB,B)
                     ENDIF
C
                     CALL T3_FORBIDDEN(WORK(KSMAT4),ISYMIM,
     *                                 ISYMD,D,ISYMB,B)
C
                     DTIME  = SECOND() - DTIME
                     TISMAT = TISMAT   + DTIME
C
                     IF (IPRINT .GT. 55) THEN
                        XSMAT = DDOT(NCKIJ(ISCKIJ),WORK(KSMAT4),1,
     *                          WORK(KSMAT4),1)
                        WRITE(LUPRI,*) 'Norm of SMAT-BAR-2 ',XSMAT
                     ENDIF
C
                  ENDIF
C
C--------------------------------------------------
C                 Calculate Q(ci,jk) for fixed b,d.
C--------------------------------------------------
C
                  DTIME = SECOND()
                  CALL CC3_QMAT(0.0D0,T2TP,ISYMT2,WORK(KTRVI3),
     *                          WORK(KTROC0),ISINT2,WORK(KFOCKD),
     *                          WORK(KDIAG),WORK(KQMAT),
     *                          WORK(KEND4),LWRK4,WORK(KINDSQ),LENSQ,
     *                          ISYMB,B,ISYMD,D)
C
                  CALL T3_FORBIDDEN(WORK(KQMAT),ISYMIM,ISYMB,B,ISYMD,D)
C
                  DTIME  = SECOND() - DTIME
                  TIQMAT = TIQMAT   + DTIME
C
                  IF (IPRINT .GT. 55) THEN
                     XQMAT = DDOT(NCKIJ(ISCKIJ),WORK(KQMAT),1,
     *                       WORK(KQMAT),1)
                     WRITE(LUPRI,*) 'Norm of QMAT     ',XQMAT
                  ENDIF
C
C-------------------------------------------------------------------
C                 Calculate Q(ci,jk) for fixed b,d for t3-bar.
C-------------------------------------------------------------------
C
                  DTIME = SECOND()
C
                  CALL DZERO(WORK(KQMAT2),NCKIJ(ISCKIJ))
C
                  IF (CC3) THEN
                     CALL CCFOP_QMAT(0.0D0,L1AM,ISYML1,L2TP,ISYML2,
     *                               WORK(KTMAT),WORK(KFCKBA),
     *                               WORK(KXIAJB),ISINT1,WORK(KTRVI18),
     *                               WORK(KTRVI6),WORK(KTROC01),
     *                               WORK(KTROC21),ISINT2,WORK(KFOCKD),
     *                               WORK(KDIAG),WORK(KQMAT2),
     *                               WORK(KEND4),LWRK4,WORK(KINDSQ),
     *                               LENSQ,ISYMB,B,ISYMD,D)
C
                     CALL DSCAL(NCKIJ(ISCKIJ),HALF,WORK(KQMAT2),1)
C
                  ELSE
                     CALL CCFOP_QMAT(0.0D0,T1AM,ISYMT1,
     *                               WORK(KT2TCME),ISYMT2,
     *                               WORK(KTMAT),WORK(KFCKBA),
     *                               WORK(KXIAJB),ISINT1,WORK(KTRVI3),
     *                               WORK(KTRVI6),WORK(KTROC0),
     *                               WORK(KTROC2),ISINT2,WORK(KFOCKD),
     *                               WORK(KDIAG),WORK(KQMAT2),
     *                               WORK(KEND4),LWRK4,WORK(KINDSQ),
     *                               LENSQ,ISYMB,B,ISYMD,D)
                  ENDIF
C
                  CALL T3_FORBIDDEN(WORK(KQMAT2),ISYMIM,ISYMB,B,ISYMD,D)
C
                  DTIME  = SECOND() - DTIME
                  TIQMAT = TIQMAT   + DTIME
C
                  IF (IPRINT .GT. 55) THEN
                     XQMAT = DDOT(NCKIJ(ISCKIJ),WORK(KQMAT2),1,
     *                       WORK(KQMAT2),1)
                     WRITE(LUPRI,*) 'Norm of QMAT-BAR ',XQMAT
                  ENDIF
C
C--------------------------------------------------
C                 Calculate U(ci,jk) for fixed b,d.
C--------------------------------------------------
C
                  DTIME = SECOND()
                  CALL CC3_UMAT(0.0D0,T2TP,ISYMT2,WORK(KTRVI1),
     *                          WORK(KTROC02),ISINT2,WORK(KFOCKD),
     *                          WORK(KDIAG),WORK(KUMAT),WORK(KTMAT),
     *                          WORK(KEND4),LWRK4,WORK(KINDSQ),LENSQ,
     *                          ISYMB,B,ISYMD,D)
C
                  CALL T3_FORBIDDEN(WORK(KUMAT),ISYMIM,ISYMB,B,ISYMD,D)
C
                  DTIME  = SECOND() - DTIME
                  TIQMAT = TIQMAT   + DTIME
C
                  IF (IPRINT .GT. 55) THEN
                     XQMAT = DDOT(NCKIJ(ISCKIJ),WORK(KUMAT),1,
     *                       WORK(KUMAT),1)
                     WRITE(LUPRI,*) 'Norm of UMAT     ',XQMAT
                  ENDIF
C
C--------------------------------------------------
C                 Calculate U(ci,jk) for fixed d,b.
C--------------------------------------------------
C
                  DTIME = SECOND()
                  CALL CC3_UMAT(0.0D0,T2TP,ISYMT2,WORK(KTRVI10),
     *                          WORK(KTROC02),ISINT2,WORK(KFOCKD),
     *                          WORK(KDIAG),WORK(KUMAT3),WORK(KTMAT),
     *                          WORK(KEND4),LWRK4,WORK(KINDSQ),LENSQ,
     *                          ISYMD,D,ISYMB,B)
C
                  CALL T3_FORBIDDEN(WORK(KUMAT3),ISYMIM,ISYMD,D,ISYMB,B)
C
                  DTIME  = SECOND() - DTIME
                  TIQMAT = TIQMAT   + DTIME
C
                  IF (IPRINT .GT. 55) THEN
                     XQMAT = DDOT(NCKIJ(ISCKIJ),WORK(KUMAT),1,
     *                       WORK(KUMAT),1)
                     WRITE(LUPRI,*) 'Norm of UMAT3    ',XQMAT
                  ENDIF
C
C-----------------------------------------------------------------
C                 Calculate U(ci,jk) for fixed b,d for t3-bar.
C-----------------------------------------------------------------
C
                  DTIME = SECOND()
C
                  CALL DZERO(WORK(KUMAT2),NCKIJ(ISCKIJ))
C
                  IF (CC3) THEN
                     CALL CCFOP_UMAT(0.0D0,L1AM,ISYML1,L2TP,ISYML2,
     *                               WORK(KXIAJB),ISINT1,WORK(KFCKBA),
     *                               WORK(KTRVI19),WORK(KTRVI7),
     *                               WORK(KTROC03),WORK(KTROC23),ISINT2,
     *                               WORK(KFOCKD),WORK(KDIAG),
     *                               WORK(KUMAT2),
     *                               WORK(KTMAT),WORK(KEND4),LWRK4,
     *                               WORK(KINDSQ),LENSQ,ISYMB,B,ISYMD,D)
C
                     CALL DSCAL(NCKIJ(ISCKIJ),HALF,WORK(KUMAT2),1)
C
                  ELSE
                     CALL CCFOP_UMAT(0.0D0,T1AM,ISYMT1,WORK(KT2TCME),
     *                               ISYMT2,
     *                               WORK(KXIAJB),ISINT1,WORK(KFCKBA),
     *                               WORK(KTRVI1),WORK(KTRVI7),
     *                               WORK(KTROC02),WORK(KTROC22),ISINT2,
     *                               WORK(KFOCKD),WORK(KDIAG),
     *                               WORK(KUMAT2),
     *                               WORK(KTMAT),WORK(KEND4),LWRK4,
     *                               WORK(KINDSQ),LENSQ,ISYMB,B,ISYMD,D)
                  ENDIF
C
                  CALL T3_FORBIDDEN(WORK(KUMAT2),ISYMIM,ISYMB,B,ISYMD,D)
C
                  DTIME  = SECOND() - DTIME
                  TIQMAT = TIQMAT   + DTIME
C
                  IF (IPRINT .GT. 55) THEN
                     XQMAT = DDOT(NCKIJ(ISCKIJ),WORK(KUMAT2),1,
     *                       WORK(KUMAT2),1)
                     WRITE(LUPRI,*) 'Norm of UMAT-BAR-1 ',XQMAT
                  ENDIF
C
C-----------------------------------------------------------------
C                 Calculate U(ci,jk) for fixed d,b for t3-bar.
C-----------------------------------------------------------------
C
                  IF (.NOT. RELORB) THEN
C
                     DTIME = SECOND()
C
                     CALL DZERO(WORK(KUMAT4),NCKIJ(ISCKIJ))
C
                     IF (CC3) THEN
                        CALL CCFOP_UMAT(0.0D0,L1AM,ISYML1,L2TP,ISYML2,
     *                                  WORK(KXIAJB),ISINT1,
     *                                  WORK(KFCKBA),WORK(KTRVI20),
     *                                  WORK(KTRVI13),WORK(KTROC03),
     *                                  WORK(KTROC23),ISINT2,
     *                                  WORK(KFOCKD),WORK(KDIAG),
     *                                  WORK(KUMAT4),WORK(KTMAT),
     *                                  WORK(KEND4),LWRK4,WORK(KINDSQ),
     *                                  LENSQ,ISYMD,D,ISYMB,B)
C
                     CALL DSCAL(NCKIJ(ISCKIJ),HALF,WORK(KUMAT4),1)
C
                     ELSE
                        CALL CCFOP_UMAT(0.0D0,T1AM,ISYMT1,WORK(KT2TCME),
     *                                  ISYMT2,WORK(KXIAJB),ISINT1,
     *                                  WORK(KFCKBA),WORK(KTRVI10),
     *                                  WORK(KTRVI13),WORK(KTROC02),
     *                                  WORK(KTROC22),ISINT2,
     *                                  WORK(KFOCKD),WORK(KDIAG),
     *                                  WORK(KUMAT4),WORK(KTMAT),
     *                                  WORK(KEND4),LWRK4,WORK(KINDSQ),
     *                                  LENSQ,ISYMD,D,ISYMB,B)
                     ENDIF
C
                     CALL T3_FORBIDDEN(WORK(KUMAT4),ISYMIM,
     *                                 ISYMD,D,ISYMB,B)
C
                     DTIME  = SECOND() - DTIME
                     TIQMAT = TIQMAT   + DTIME
C
                     IF (IPRINT .GT. 55) THEN
                        XQMAT = DDOT(NCKIJ(ISCKIJ),WORK(KUMAT4),1,
     *                          WORK(KUMAT4),1)
                        WRITE(LUPRI,*) 'Norm of UMAT-BAR-2 ',XQMAT
                     ENDIF
C
                  ENDIF
C
C-----------------------------------------------------------
C                 Construct Kappabar_{aa} and Kappabar_{ii}
C-----------------------------------------------------------
C
                  IF ((.NOT. CC3) .AND. (RELORB)) THEN
C
                     CALL CCSDT_KAPPADIAG(WORK(KKAPAA),WORK(KKAPII),
     *                                    WORK(KSMAT2),WORK(KSMAT),
     *                                    WORK(KSMAT3),WORK(KUMAT2),
     *                                    WORK(KUMAT),WORK(KUMAT3),
     *                                    WORK(KTMAT),WORK(KINDSQ),
     *                                    LENSQ,ISCKIJ,
     *                                    WORK(KEND4),LWRK4)
C
                  ENDIF
C
C----------------------------------------------------------------
C                 Calculate the three extra contributions to the 
C                 one-electron density if nonrelaxed
C----------------------------------------------------------------
C
                  IF (.NOT. RELORB) THEN
                     CALL CCFOP_NONREL(WORK(KOMG12),WORK(KDENSAB),
     *                                 WORK(KDENSIJ),ISCKIJ,
     *                                 WORK(KSMAT),WORK(KSMAT3),
     *                                 WORK(KSMAT2),WORK(KSMAT4),
     *                                 WORK(KUMAT),WORK(KUMAT3),
     *                                 WORK(KUMAT2),WORK(KUMAT4),
     *                                 WORK(KTMAT),T2TP,ISYMT2,
     *                                 WORK(KINDSQ),LENSQ,
     *                                 ISYMB,B,ISYMD,D,
     *                                 WORK(KEND4),LWRK4)
                  ENDIF
C
C---------------------------------------------
C                 Contract with integrals.
C---------------------------------------------
C
                  DTIME = SECOND()
C
                  IF ((.NOT. CC3) .AND. (RELORB)) THEN
C
                     CALL CCFOP_DENVIR_SC(WORK(KVIR1),WORK(KVIR2),
     *                                 WORK(KSMAT),WORK(KQMAT),
     *                                 WORK(KTMAT),ISYMIM,
     *                                 WORK(KT2TCME),ISYMT2,
     &                                 WORK(KEND4),LWRK4,
     &                                 WORK(KINDSQ),LENSQ,
     *                                 ISYMB,B,ISYMD,D,1)
C
                     IF ((IPRINT .GT. 55)) THEN
                        XRMAT = DDOT(NCKATR(ISAIJ1),WORK(KVIR1),1,
     *                          WORK(KVIR1),1)
                        WRITE(LUPRI,*) 'Norm DENS1 - CCFOP_DENVIR',XRMAT
                        XRMAT = DDOT(NCKATR(ISAIJ1),WORK(KVIR2),1,
     *                          WORK(KVIR2),1)
                        WRITE(LUPRI,*) 'Norm DENS2 - CCFOP_CONVIR',XRMAT
                     ENDIF
C
                     CALL CCFOP_DENVIR_SC(WORK(KVIR3),WORK(KVIR4),
     *                                 WORK(KSMAT2),WORK(KQMAT2),
     *                                 WORK(KTMAT),ISYMIM,
     *                                 T2TP,ISYMT2,WORK(KEND4),
     *                                 LWRK4,WORK(KINDSQ),LENSQ,
     *                                 ISYMB,B,ISYMD,D,2)
C
                     IF ((IPRINT .GT. 55)) THEN
                        XRMAT = DDOT(NCKATR(ISAIJ1),WORK(KVIR1),1,
     *                          WORK(KVIR1),1)
                        WRITE(LUPRI,*) 'Norm DENS1 - CCFOP_DENVIR',XRMAT
                        XRMAT = DDOT(NCKATR(ISAIJ1),WORK(KVIR2),1,
     *                          WORK(KVIR2),1)
                        WRITE(LUPRI,*) 'Norm DENS2 - CCFOP_CONVIR',XRMAT
                     ENDIF
C
C
                     CALL CCFOP_DENOCC_SC(WORK(KOCC1),WORK(KSMAT),
     *                                 WORK(KQMAT),WORK(KTMAT),ISYMIM,
     *                                 WORK(KT2TCME),ISYMT2,WORK(KEND4),
     *                                 LWRK4,WORK(KINDSQ),LENSQ,
     *                                 ISYMB,B,ISYMD,D,1)

C
                     CALL CCFOP_DENOCC_SC(WORK(KOCC2),WORK(KSMAT2),
     *                                 WORK(KQMAT2),WORK(KTMAT),ISYMIM,
     *                                 T2TP,ISYMT2,WORK(KEND4),
     *                                 LWRK4,WORK(KINDSQ),LENSQ,
     *                                 ISYMB,B,ISYMD,D,2)
C

                     IF ((IPRINT .GT. 55)) THEN
                        XRMAT = DDOT(NCKIJ(ISYRES),WORK(KOCC1),1,
     *                          WORK(KOCC1),1)
                        WRITE(LUPRI,*) 'Norm DENS1 - CCFOP_DENOCC',XRMAT
                        XRMAT = DDOT(NCKIJ(ISYRES),WORK(KOCC2),1,
     *                          WORK(KOCC2),1)
                        WRITE(LUPRI,*) 'Norm DENS2 - CCFOP_DENOCC',XRMAT
                     END IF
C---------------------------------------
C                 Calculate Omega22.
C---------------------------------------
C
                     DTIME = SECOND()
C
                     CALL CCFOP_ONEL(WORK(KOMG22),WORK(KRMAT1),
     *                               WORK(KRMAT2),T1AM,WORK(KSMAT),
     *                               WORK(KTMAT),ISYMIM,ISINT1,
     *                               WORK(KINDSQ),LENSQ,
     *                               WORK(KEND4),LWRK4,ISYMB,B,ISYMD,D)
C
                     IF ((IPRINT .GT. 55)) THEN
                        RHO2N = DDOT(NT2AM(ISYRES),WORK(KOMG22),1,
     *                                             WORK(KOMG22),1)
                        WRITE(LUPRI,*) 'Norm of Rho22 after CC3_ONEL',
     *                                  RHO2N
                     ENDIF
C
                     IF (IPRINT .GT. 220) THEN
                        CALL AROUND('After CC3_ONEL: ')
                        CALL CC_PRP(DUMMY,WORK(KOMG22),ISYRES,0,1)
                     ENDIF
C
                     IF (IPRINT .GT. 55) THEN
                        XRMAT = DDOT(NCKI(ISAIJ1),WORK(KRMAT1),1,
     *                          WORK(KRMAT1),1)
                        WRITE(LUPRI,*) 'Norm of RMAT1 -after ONEL',XRMAT
                        XRMAT = DDOT(NCKI(ISAIJ2),WORK(KRMAT2),1,
     *                          WORK(KRMAT2),1)
                        WRITE(LUPRI,*) 'Norm of RMAT2 -after ONEL',XRMAT
                        XTMAT = DDOT(NCKIJ(ISCKIJ),WORK(KSMAT),1,
     *                          WORK(KSMAT),1)
                        WRITE(LUPRI,*) 'Norm of SMAT -after ONEL',XSMAT
                        XTMAT = DDOT(NCKIJ(ISCKIJ),WORK(KTMAT),1,
     *                          WORK(KTMAT),1)
                        WRITE(LUPRI,*) 'Norm of TMAT -after ONEL',XTMAT
                     ENDIF
C
                  ENDIF   ! RELORB
C
C---------------------------------------------------
C                 Calculate Omega1.
C---------------------------------------------------
C
                  DTIME  = SECOND() - DTIME
                  TIOME1 = TIOME1   + DTIME
C
                  IF (CC3) THEN
                     CALL CCFOP_ONED(WORK(KOMG1),L2TP,ISYML2,
     *                               WORK(KSMAT),WORK(KTMAT),ISYMIM,
     *                               WORK(KINDSQ),LENSQ,
     *                               WORK(KEND4),LWRK4,ISYMB,B,ISYMD,D)
                  ELSE
                     CALL CCFOP_ONED(WORK(KOMG1),WORK(KT2TCME),ISYMT2,
     *                               WORK(KSMAT),WORK(KTMAT),ISYMIM,
     *                               WORK(KINDSQ),LENSQ,
     *                               WORK(KEND4),LWRK4,ISYMB,B,ISYMD,D)
                  ENDIF
C
                  IF ((IPRINT .GT. 55)) THEN
                     XT2TP = DDOT(NT1AM(ISYMOP),WORK(KOMG1),1,
     *                            WORK(KOMG1),1)
                     WRITE(LUPRI,*) 'Norm of 1 e- density : ',XT2TP
                  ENDIF
C
                  IF (IPRINT .GT. 220) THEN
                     CALL AROUND('After CCFOP_ONED: ')
                     CALL CC_PRP(WORK(KOMG1),DUMMY,ISYRES,1,0)
                  ENDIF
C
                  DTIME  = SECOND() - DTIME
                  TIOME1 = TIOME1   + DTIME
C
C---------------------------------------------------------
C                 Accumulate the R2 matrix in Omega22
C---------------------------------------------------------
C
                  IF ((.NOT. CC3) .AND. (RELORB)) THEN
C
                     CALL CC3_RACC(WORK(KOMG22),WORK(KRMAT2),ISYMB,B,
     *                             ISYRES)
C
                     IF ((IPRINT .GT. 55)) THEN
                        RHO2N = DDOT(NT2AM(ISYRES),WORK(KOMG22),1,
     *                                             WORK(KOMG22),1)
                        WRITE(LUPRI,*) 'Norm of Rho22-after CC3_RACC',
     *                                  RHO2N
                     ENDIF
C
                     IF (IPRINT .GT. 220) THEN
                        CALL AROUND('After CC3_RACC: ')
                        CALL CC_PRP(DUMMY,WORK(KOMG22),ISYRES,0,1)
                     ENDIF
C
                  ENDIF
C
               ENDDO   ! B
            ENDDO      ! ISYMB
C
C---------------------------------------------------
C           Accumulate the R1 matrix in Omega22.
C---------------------------------------------------
C
            IF ((.NOT. CC3) .AND. (RELORB)) THEN
C
               CALL CC3_RACC(WORK(KOMG22),WORK(KRMAT1),ISYMD,D,ISYRES)
C
               IF (IPRINT .GT. 55) THEN
                  RHO2N = DDOT(NT2AM(ISYRES),WORK(KOMG22),1,
     *                                       WORK(KOMG22),1)
                  WRITE(LUPRI,*) 'Norm of Rho22-after CC3_RACC-2',RHO2N
               ENDIF
C
               IF (IPRINT .GT. 220) THEN
                  CALL AROUND('After CC3_RACC-2: ')
                  CALL CC_PRP(DUMMY,WORK(KOMG22),ISYRES,0,1)
               ENDIF
C
C--------------------------------------------------------------
C        Sort the two electron densities from T3 for a constant
C        d and write them to file.
C--------------------------------------------------------------
C
               IF (LWRK4 .LT. NCKATR(ISAIJ1)) THEN
                  CALL QUIT('Exceeded memory in CCSDPT_DENS2 (sort)')
               ENDIF
!
!Sonia: following not working with symmetry
C
!               CALL DEN_AIBSORT(WORK(KVIR1),WORK(KEND4),ISAIJ1)
C
!               CALL DEN_AIBSORT(WORK(KVIR2),WORK(KEND4),ISAIJ1)
C
                CALL DEN_BIASORT(WORK(KVIR1),WORK(KEND4),ISAIJ1)
C
                CALL DEN_BIASORT(WORK(KVIR2),WORK(KEND4),ISAIJ1)
C
               IOFF = ICKBD(ISAIJ1,ISYMD)
     *              + NCKATR(ISAIJ1)*(D-1)
     *              + 1
               CALL PUTWA2(LUABI2,FNDABI2,WORK(KVIR1),IOFF,
     *                     NCKATR(ISAIJ1))
C
               CALL PUTWA2(LUABI1,FNDABI1,WORK(KVIR2),IOFF,
     *                     NCKATR(ISAIJ1))
C
               IF (IPRINT .GT. 55) THEN
                  RHO1N = DDOT(NCKATR(ISAIJ1),WORK(KVIR1),1,
     *                         WORK(KVIR1),1)
                  WRITE(LUPRI,*) 'Norm of VIR1 : ',RHO1N
                  RHO1N = DDOT(NCKATR(ISAIJ1),WORK(KVIR2),1,
     *                         WORK(KVIR2),1)
                  WRITE(LUPRI,*) 'Norm of VIR2 : ',RHO1N
               ENDIF
C
C----------------------------------------------------------------------
C        Sort the two electron densities from T3-bar for a constant
C        d and write them to file.
C----------------------------------------------------------------------
C
               IF (LWRK4 .LT. NCKATR(ISAIJ1)) THEN
                  CALL QUIT('Exceeded memory in CCSDPT_DENS2 (sort)')
               ENDIF
C
!
!Sonia: following not working with symmetry
C
!               CALL DEN_AIBSORT(WORK(KVIR3),WORK(KEND4),ISAIJ1)
C
!               CALL DEN_AIBSORT(WORK(KVIR4),WORK(KEND4),ISAIJ1)
C
                CALL DEN_BIASORT(WORK(KVIR3),WORK(KEND4),ISAIJ1)
C
                CALL DEN_BIASORT(WORK(KVIR4),WORK(KEND4),ISAIJ1)
C
               IOFF = ICKBD(ISAIJ1,ISYMD)
     *              + NCKATR(ISAIJ1)*(D-1)
     *              + 1
               CALL PUTWA2(LUABI4,FNDABI4,WORK(KVIR3),IOFF,
     *                     NCKATR(ISAIJ1))
C
               CALL PUTWA2(LUABI3,FNDABI3,WORK(KVIR4),IOFF,
     *                     NCKATR(ISAIJ1))
C
               IF (IPRINT .GT. 55) THEN
                  RHO1N = DDOT(NCKATR(ISAIJ1),WORK(KVIR3),1,
     *                         WORK(KVIR3),1)
                  WRITE(LUPRI,*) 'Norm of VIR3 : ',RHO1N
                  RHO1N = DDOT(NCKATR(ISAIJ1),WORK(KVIR4),1,
     *                         WORK(KVIR4),1)
                  WRITE(LUPRI,*) 'Norm of VIR4 : ',RHO1N
               ENDIF
C
            ENDIF    ! RELORB
C
         ENDDO       ! D
      ENDDO          ! ISYMD
C
C -------------------------------------------------------------
C The following section is taken from old code, as the new does
C not seem to work with symmetry
C Sonia, Dec 2002
C -------------------------------------------------------------
C
      IF ((.NOT. CC3)) THEN

         KEND4 = KENDSV
         LWRK4 = LWORK - KEND4

         IF (RELORB) THEN
C
            NTOT = 0
            DO ISYMD = 1, NSYM
               IF (NT1AM(ISYMD) .GT. NTOT) THEN
                  NTOT = NT1AM(ISYMD)
               ENDIF
            ENDDO
C
            KVIR1 = KEND4
            KVIR2 = KVIR1 + NTOT
            KEND4 = KVIR2 + NTOT
            LWRK4 = LWORK - KEND4

C--------------------------------------------
C     Resort the densities occ1 and occ2
C--------------------------------------------
      IF ((IPRINT.GT.55).or.(locdbg)) THEN
        RHO2N = ddot(nckij(isyres),work(kocc1),1,work(kocc1),1)
        WRITE(lupri,*) 'Norm of Occ1 (before cc3_lsort1) : ',rho2n
        RHO2N = ddot(nckij(isyres),work(kocc2),1,work(kocc2),1)
        WRITE(lupri,*) 'Norm of Occ2 (before cc3_lsort1) : ',rho2n
      END IF
C
            CALL CC3_LSORT1(WORK(KOCC1),ISYRES,WORK(KEND4),LWRK4,3)
C
            CALL CC3_LSORT1(WORK(KOCC2),ISYRES,WORK(KEND4),LWRK4,3)
C
      IF ((IPRINT.GT.55).or.(locdbg)) THEN
        rho2n = ddot(nckij(isyres),work(kocc1),1,work(kocc1),1)
        write(lupri,*) 'Norm of Occ1 (after  cc3_lsort1) : ',rho2n
        rho2n = ddot(nckij(isyres),work(kocc2),1,work(kocc2),1)
        write(lupri,*) 'Norm of Occ2 (after  cc3_lsort1) : ',rho2n
      END IF

         ENDIF
      END IF
C
C-------------------------------------- End of old part....
C---------------------------------------------------------
C     Construct 2*C-E of work(komg22) and write to file.
C---------------------------------------------------------
C
      IF ((.NOT. CC3) .AND. (RELORB)) THEN
         IOPTTCME = 1
         ISYOPE   = ISYRES
         CALL CCSD_TCMEPK(WORK(KOMG22),1.0D0,ISYOPE,IOPTTCME)
C
         IF ((IPRINT .GT. 55)) THEN
            RHO2N = DDOT(NT2AM(ISYRES),WORK(KOMG22),1,WORK(KOMG22),1)
            WRITE(LUPRI,*) 'Norm of Rho22 at the end     ',RHO2N
         ENDIF
C
         IF (IPRINT .GT. 100) THEN
            CALL AROUND('TWO ELECTRON DENSITY : D_{IAJB}')
            CALL CC_PRP(T1AM,WORK(KOMG22),ISYRES,0,1)
         ENDIF
C
         IF (NT2AM(ISYRES) .GT. 0) THEN
            IOFF = 1
            CALL PUTWA2(LUPTIAJB,FNDIAJB,WORK(KOMG22),IOFF,
     *                  NT2AM(ISYRES))
         ENDIF
C
         IF (LDEBUG .AND. (.NOT. CC3)) THEN
            LENGTH = IRAT*NT2AM(ISYRES)
C
            REWIND(LUIAJB)
            CALL READI(LUIAJB,LENGTH,WORK(KXIAJB))
            CALL CCLR_DIASCL(WORK(KXIAJB),0.5D0,ISYMop)
            CALL DSCAL(NT2AM(ISYMOP),2.0D0,WORK(KOMG22),1)
C
            XQMAT = DDOT(NT2AM(ISYRES),WORK(KOMG22),1,WORK(KXIAJB),1)
            WRITE(LUPRI,*) 'DEBUGGING CCSD(T) : E5 = ',XQMAT
         ENDIF
C
      ENDIF
C
C---------------------------------------
C     Scale and store the 1 e- density :
C---------------------------------------
C
      IF (CC3) THEN
         CALL DSCAL(NT1AM(ISYRES),-ONE,WORK(KOMG1),1)
      ELSE
         CALL DSCAL(NT1AM(ISYRES),-TWO,WORK(KOMG1),1)
      ENDIF
C
      IF (IPRINT .GT. 55) THEN
         RHO1N = DDOT(NT1AM(ISYRES),WORK(KOMG1),1,WORK(KOMG1),1)
         WRITE(LUPRI,*) 'Norm of OMEG1 at the end : ',RHO1N
      ENDIF
C
      IF (IPRINT .GT. 100) THEN
         CALL AROUND('1 e- in CCSDPT_DENS2 : ')
         CALL CC_PRP(WORK(KOMG1),DUMMY,ISYRES,1,0)
      ENDIF
C
      IF (NT1AM(ISYRES) .GT. 0) THEN
         IOFF = 1
         CALL PUTWA2(LUPTIA,FNDPTIA,WORK(KOMG1),IOFF,NT1AM(ISYRES))
      ENDIF
C
C----------------------------------------------------------
C      Add the 1e- to the d_{iajk} for j=k and for i=k
C----------------------------------------------------------
C
      IF ((IPRINT .GT. 55).or.locdbg) THEN
         RHO1N = DDOT(NCKIJ(ISYRES),WORK(KOCC1),1,WORK(KOCC1),1)
         WRITE(LUPRI,*) 'Norm of OCC1 (iajk) (before dens1to2) =',RHO1N
      ENDIF
C
      IF ((.NOT. CC3) .AND. (RELORB)) THEN
         CALL DENS1TO2(WORK(KOMG1),WORK(KOCC1),ISYRES)
      ENDIF
C
      IF ((IPRINT .GT. 55).or.locdbg) THEN
         RHO1N = DDOT(NCKIJ(ISYRES),WORK(KOCC1),1,WORK(KOCC1),1)
         WRITE(LUPRI,*) 'Norm of OCC1 (iajk) (after dens1to2)  =',RHO1N
      ENDIF
C
C-----------------------------------------------------------
C     If nonrel store the three extra terms on disc
C-----------------------------------------------------------
C
      IF (.NOT. RELORB) THEN
C
         IF (NMATAB(ISYRES) .GT. 0) THEN
           IOFF = 1
           CALL PUTWA2(LUPTAB,FNDPTAB,WORK(KDENSAB),IOFF,NMATAB(ISYRES))
         ENDIF
C
         IF (NMATIJ(ISYRES) .GT. 0) THEN
           IOFF = 1
           CALL PUTWA2(LUPTIJ,FNDPTIJ,WORK(KDENSIJ),IOFF,NMATIJ(ISYRES))
         ENDIF
C
         CALL DSCAL(NT1AM(ISYRES),-TWO,WORK(KOMG12),1)
         IF (NT1AM(ISYRES) .GT. 0) THEN
           IOFF = 1
           CALL PUTWA2(LUPTIA2,FNDPTIA2,WORK(KOMG12),IOFF,NT1AM(ISYRES))
         ENDIF
C
      ENDIF
C
C---------------------------------------------------------------
C     Construct the total d(ab,ic) density stored as (ai,b,c)
C     from the T3 amplitudes.
C---------------------------------------------------------------
C
      IF ((.NOT. CC3) .AND. (RELORB)) THEN
C
         CALL DENSTORE(WORK(KVIR2),LUABI1,FNDABI1,
     *                 WORK(KVIR1),LUABI2,FNDABI2,ISYRES)
C
C
C---------------------------------------------------------------
C     Construct the total d(ab,ci) density stored as (bi,a,c)
C     from the T3-bar amplitudes.
C---------------------------------------------------------------
C
         CALL DENSTORE(WORK(KVIR4),LUABI3,FNDABI3,
     *                 WORK(KVIR3),LUABI4,FNDABI4,ISYRES)
C
C-----------------------------------------
C     Store the d_{iajk} as kjia
C-----------------------------------------
C
         IF (NCKIJ(ISYRES) .GT. 0) THEN
            IOFF = 1
            CALL PUTWA2(LUIAJK,FNDIAJK,WORK(KOCC1),IOFF,NCKIJ(ISYRES))
         ENDIF
C
C-----------------------------------------
C     Store the d_{aijk} as jkia
C-----------------------------------------
C
      IF ((IPRINT .GT. 55).or.locdbg) THEN
         RHO1N = DDOT(NCKIJ(ISYRES),WORK(KOCC2),1,WORK(KOCC2),1)
         WRITE(LUPRI,*) 'Norm of OCC2 (aijk) = ',RHO1N
      ENDIF
C
         IF (NCKIJ(ISYRES) .GT. 0) THEN
            IOFF = 1
            CALL PUTWA2(LUAIJK,FNDAIJK,WORK(KOCC2),IOFF,NCKIJ(ISYRES))
         ENDIF
C
C------------------------------------------
C     Store kappabar_{aa} and kappabar_{ii}
C------------------------------------------
C
         IF ((IPRINT .GT. 55).or.locdbg) THEN
            RHO1N = DDOT(NRHFT,WORK(KKAPII),1,WORK(KKAPII),1)
            WRITE(LUPRI,*) 'Norm of KAPII : ',RHO1N
            RHO1N = DDOT(NVIRT,WORK(KKAPAA),1,WORK(KKAPAA),1)
            WRITE(LUPRI,*) 'Norm of KAPAA : ',RHO1N
         ENDIF
C
         LUKAPAB  = -1
         LUKAPIJ  = -1
         FNKAPAB  = 'KAPAB'
         FNKAPIJ  = 'KAPIJ'
         CALL WOPEN2(LUKAPAB,FNKAPAB,64,0)
         CALL WOPEN2(LUKAPIJ,FNKAPIJ,64,0)
C
         IF (NVIRT .GT. 0) THEN
            IOFF = 1
            CALL PUTWA2(LUKAPAB,FNKAPAB,WORK(KKAPAA),IOFF,NVIRT)
         ENDIF
C
         IF (NRHFT .GT. 0) THEN
            IOFF = 1
            CALL PUTWA2(LUKAPIJ,FNKAPIJ,WORK(KKAPII),IOFF,NRHFT)
         ENDIF
C
         CALL WCLOSE2(LUKAPAB,FNKAPAB,'KEEP')
         CALL WCLOSE2(LUKAPIJ,FNKAPIJ,'KEEP')
C
      ENDIF   ! RELORB
C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C     SONIA: symmetrize/reorder some two-electron densities
C     before closing, and generate backtransformed ones!!!
C     Symmetrize the vir.vir.occ.vir and vir.vir.vir.occ
C     Symmetrize the occ.occ.occ.vir and occ.occ.vir.occ
C     and backtransform last index to delta
C----------------------------------------------------------
C
      IF ((.NOT. CC3) .AND. RELORB) THEN
         if (.false.) then
            CALL SYMMBACK(MODEL,LUIAJK,FNDIAJK,LUAIJK,FNDAIJK,
     *                 LUABI1,FNDABI1,LUABI3,FNDABI3,
     *                 LUPTIAJB,FNDIAJB,
     *                 ISYRES,WORK(KEND4),LWRK4)
         else
            CALL SYMMBACK_1(MODEL,LUIAJK,FNDIAJK,LUAIJK,FNDAIJK,
     *                 LUABI1,FNDABI1,LUABI3,FNDABI3,
     *                 LUPTIAJB,FNDIAJB,
     *                 ISYRES,WORK(KEND4),LWRK4)
         end if
      ENDIF
C
C---------------------------------------
C     Close files.
C---------------------------------------
C
      CALL WCLOSE2(LUPTIA,FNDPTIA,'KEEP')
C
      IF ((.NOT. CC3) .AND. (RELORB)) THEN
         CALL WCLOSE2(LUPTIAJB,FNDIAJB,'KEEP')
         CALL WCLOSE2(LUABI1,FNDABI1,'KEEP')
         CALL WCLOSE2(LUABI2,FNDABI2,'DELETE')
         CALL WCLOSE2(LUABI3,FNDABI3,'KEEP')
         CALL WCLOSE2(LUABI4,FNDABI4,'DELETE')
         CALL WCLOSE2(LUAIJK,FNDAIJK,'KEEP')
         CALL WCLOSE2(LUIAJK,FNDIAJK,'KEEP')
      ELSE
         CALL WCLOSE2(LUPTIA2,FNDPTIA2,'KEEP')
         CALL WCLOSE2(LUPTAB,FNDPTAB,'KEEP')
         CALL WCLOSE2(LUPTIJ,FNDPTIJ,'KEEP')
      ENDIF
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
      CALL QEXIT('CCSDPT_DENS2_SC')
C
      RETURN
C
    1 FORMAT(7X,'Time used in',2X,A12,F12.2,' seconds')
C
      END
C  /* Deck ccfop_denvir_sc */
      SUBROUTINE CCFOP_DENVIR_SC(RINTE1,RINTE2,SMAT,QMAT,TMAT,ISYMIM,
     *                        T2TCME,ISYMT2,WORK,LWORK,INDSQ,LENSQ,
     *                        ISYMB,B,ISYMD,D,IOPT)
C
C     Kasper Hald, Fall 2001.
C     OLd version working for CCSD(T) gradient, Sonia Coriani
C
C     Calculate the two electron density (abic) for a constant index D,
C     and add to the density RINTE1 and RINTE2.
C
C     ISYMIM is the symmetry of the SMAT and TMAT intermdiates.
C     ISYMT2 is the symmetry of the T2 amplitudes.
C
C     IOPT = 1. Calculate the terms from T3AM.
C     IOPT = 2. Calculate the terms from T3BAR.
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INTEGER ISYMIM, ISYMT2, LWORK, LENSQ, ISYMB, ISYMD, IOPT
      INTEGER INDSQ(LENSQ,6)
      INTEGER INDEX, ISYRES, ISYMBD, ISCKIJ, LENGTH, ISYAIJ, ISYMAI
      INTEGER ISYMA, ISYMIJ, ISYMI, ISYMJ, NAI, KOFF1, KOFF2, KOFF3
      INTEGER ISYMK, ISYCIJ, ISYMC, ISYMAB, ISYMCK, NTOTCK, NTOTIJ
      INTEGER ISYBIJ, ISYMBJ, NTOTA, KEND1, LWRK1
C
#if defined (SYS_CRAY)
      REAL RINTE1(*), RINTE2(*), SMAT(*), QMAT(*)
      REAL TMAT(*), T2TCME(*), WORK(LWORK)
      REAL ZERO, ONE, TWO, HALF
#else
      DOUBLE PRECISION RINTE1(*), RINTE2(*), SMAT(*), QMAT(*)
      DOUBLE PRECISION TMAT(*), T2TCME(*), WORK(LWORK)
      DOUBLE PRECISION ZERO, ONE, TWO, HALF
#endif
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, HALF = 0.5D0)
C
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CCFOP_DENVIR_SC')
C
C-----------------------------------------------
C     Sanity check and symmetry calculation.
C-----------------------------------------------
C
      IF (IOPT .NE. 1 .AND. IOPT .NE. 2) THEN
         CALL QUIT('Wrong IOPT in CCFOP_DENVIR_SC')
      ENDIF
C
      ISYRES = MULD2H(ISYMIM,ISYMT2)
C
      ISYMBD = MULD2H(ISYMB,ISYMD)
      ISCKIJ = MULD2H(ISYMBD,ISYMIM)
C
      LENGTH = NCKIJ(ISCKIJ)
C
C
C----------------------------------------
C     Sort the T2 for a constant B
C----------------------------------------
C
      ISYAIJ = MULD2H(ISYMT2,ISYMB)
C
      IF (LWORK .LT. NCKATR(ISYAIJ)) THEN
         CALL QUIT('Exceeded work memory in CCFOP_DENVIR')
      ENDIF
C
      DO ISYMA = 1, NSYM
         ISYMIJ = MULD2H(ISYAIJ,ISYMA)
         ISYBIJ = MULD2H(ISYMIJ,ISYMB)
         DO ISYMI = 1, NSYM
C
            ISYMJ  = MULD2H(ISYMIJ,ISYMI)
            ISYMBJ = MULD2H(ISYMJ,ISYMB)
            ISYMAI = MULD2H(ISYMA,ISYMI)
C
            DO A = 1, NVIR(ISYMA)
               DO I = 1, NRHF(ISYMI)
C
!Sonia: NAI is NOT used anywhere....
!
!                  NAI = IT1AM(ISYMA,ISYMI)
!     *                + NVIR(ISYMA)*(I-1) + A
C
                  KOFF1 =  IT2SP(ISYBIJ,ISYMA)
     *                  +  NCKI(ISYBIJ)*(A - 1)
     *                  +  ICKI(ISYMBJ,ISYMI)
     *                  +  NT1AM(ISYMBJ)*(I-1)
     *                  +  IT1AM(ISYMB,ISYMJ)
     *                  +  B
C
!                  KOFF2 =  ISAIK(ISYMAI,ISYMJ)
!     *                  +  IT1AM(ISYMA,ISYMI)
!     *                  +  NVIR(ISYMA)*(I-1)
!     *                  +  A
C
!                  CALL DCOPY(NRHF(ISYMJ),T2TCME(KOFF1),NVIR(ISYMB),
!     *                       WORK(KOFF2),NT1AM(ISYMAI))
!Sonia: replace what above with the following:
!
                  KOFF2 = IMAIJA(ISYMIJ,ISYMA)
     *                  + NMATIJ(ISYMIJ)*(A-1)
     *                  + IMATIJ(ISYMI,ISYMJ)
C     *                  + NRHF(ISYMI)*(J-1)
     *                  + I
C
                  CALL DCOPY(NRHF(ISYMJ),T2TCME(KOFF1),NVIR(ISYMB),
     *                       WORK(KOFF2),NRHF(ISYMI))
C
               ENDDO    ! I
            ENDDO       ! A
         ENDDO          ! ISYMI
      ENDDO             ! ISYMA
C
C------------------------
C     First term.
C------------------------
C
      DO I = 1,LENGTH
C
         IF (IOPT .EQ. 1) THEN
C
            TMAT(I) =  TWO*SMAT(I)
     *              -      SMAT(INDSQ(I,1))
     *              -      SMAT(INDSQ(I,5))
     *              +  TWO*QMAT(INDSQ(I,3))
     *              -      QMAT(INDSQ(I,2))
     *              -      QMAT(INDSQ(I,4))
C
         ELSE
            TMAT(I) =-HALF*SMAT(I)
     *               -HALF*QMAT(INDSQ(I,3))
         ENDIF
C
      ENDDO
!
!Sonia: this IF-block is also from old routine....
!
      IF (NSYM .GT. 1) THEN
         KEND1 = 1 + NCKI(ISYAIJ)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. LENGTH) THEN
            CALL QUIT('Out of memory in CCFOP_DENVIR (GATHER-1)')
         ENDIF
C
         CALL CC_GATHER(LENGTH,WORK(KEND1),TMAT,INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK(KEND1),1,TMAT,1)
      ENDIF
!
C------------------------------
C     Contract with T2
C------------------------------
C
!
!Sonia: replace following with old version.....
!
!      DO ISYMK = 1,NSYM
C
!         ISYCIJ = MULD2H(ISCKIJ,ISYMK)
C
!         DO ISYMC = 1, NSYM
C
!            ISYMIJ = MULD2H(ISYCIJ,ISYMC)
!            ISYMAB = MULD2H(ISYMT2,ISYMIJ)
!            ISYMA  = MULD2H(ISYMAB,ISYMB)
!            ISYMCK = MULD2H(ISYMC,ISYMK)
C
!            NTOTCK = MAX(NT1AM(ISYMCK),1)
!            NTOTIJ = MAX(NMATIJ(ISYMIJ),1)
!            NTOTA  = MAX(NVIR(ISYMA),1)
C
!            KOFF1 = ISAIK(ISYMAI,ISYMJ) + 1
!            KOFF2 = ISAIKL(ISYMCK,ISYMIJ) + 1
!            KOFF3 = ICKATR(ISYMCK,ISYMA)  + 1
C
!            CALL DGEMM('N','T',NVIR(ISYMA),NT1AM(ISYMCK),
!     *                 NMATIJ(ISYMIJ),TWO,WORK(KOFF1),NTOTA,
!     *                 TMAT(KOFF2),NTOTCK,ONE,RINTE1(KOFF3),
!     *                 NTOTA)
C
!         ENDDO         ! ISYMC
!      ENDDO            ! ISYMK

      DO ISYMCK = 1,NSYM
C
         ISYMIJ = MULD2H(ISCKIJ,ISYMCK)
         ISYMAB = MULD2H(ISYMT2,ISYMIJ)
         ISYMA  = MULD2H(ISYMAB,ISYMB)
C
         NTOTCK = MAX(NT1AM(ISYMCK),1)
         NTOTIJ = MAX(NMATIJ(ISYMIJ),1)
         NTOTA  = MAX(NVIR(ISYMA),1)
C
         KOFF1 = ISAIKL(ISYMCK,ISYMIJ) + 1
         KOFF2 = IMAIJA(ISYMIJ,ISYMA) + 1
         KOFF3 = ICKATR(ISYMCK,ISYMA)  + 1
C
         CALL DGEMM('N','N',NT1AM(ISYMCK),NVIR(ISYMA),
     *              NMATIJ(ISYMIJ),TWO,TMAT(KOFF1),NTOTCK,
     *              WORK(KOFF2),NTOTIJ,ONE,RINTE1(KOFF3),
     *              NTOTCK)
C
      ENDDO
C
C-------------------------
C     Second term.
C-------------------------
C
      DO I = 1,LENGTH
C
         IF (IOPT .EQ. 1) THEN
C
            TMAT(I) =  TWO*SMAT(INDSQ(I,1))
     *              -      SMAT(I) 
     *              -      SMAT(INDSQ(I,2))
     *              +  TWO*QMAT(INDSQ(I,2))
     *              -      QMAT(INDSQ(I,3))
     *              -      QMAT(INDSQ(I,1))
C
         ELSE
            TMAT(I) =-HALF*SMAT(INDSQ(I,1))
     *               -HALF*QMAT(INDSQ(I,2))
         ENDIF
      ENDDO
!
!Sonia: This IF-block is also from old version....
!
      IF (NSYM .GT. 1) THEN
         KEND1 = 1 + NCKI(ISYAIJ)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. LENGTH) THEN
            CALL QUIT('Out of memory in CCFOP_DENVIR (GATHER-1)')
         ENDIF
C
         CALL CC_GATHER(LENGTH,WORK(KEND1),TMAT,INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK(KEND1),1,TMAT,1)
      ENDIF
C
C------------------------------
C     Contract with T2
C------------------------------
C
!
!Sonia: As before, replace following with old version.....
!
!      DO ISYMK = 1,NSYM
C
!         ISYCIJ = MULD2H(ISCKIJ,ISYMK)
C
!         DO ISYMC = 1, NSYM
C
!            ISYMIJ = MULD2H(ISYCIJ,ISYMC)
!            ISYMAB = MULD2H(ISYMT2,ISYMIJ)
!            ISYMA  = MULD2H(ISYMAB,ISYMB)
!            ISYMCK = MULD2H(ISYMC,ISYMK)
C
!            NTOTCK = MAX(NT1AM(ISYMCK),1)
!            NTOTIJ = MAX(NMATIJ(ISYMIJ),1)
!            NTOTA  = MAX(NVIR(ISYMA),1)
C
!            KOFF1 = ISAIK(ISYMAI,ISYMJ) + 1
!            KOFF2 = ISAIKL(ISYMCK,ISYMIJ) + 1
!            KOFF3 = ICKATR(ISYMCK,ISYMA)  + 1
C
!            CALL DGEMM('N','T',NVIR(ISYMA),NT1AM(ISYMCK),
!     *                 NMATIJ(ISYMIJ),TWO,WORK(KOFF1),NTOTA,
!     *                 TMAT(KOFF2),NTOTCK,ONE,RINTE2(KOFF3),
!     *                 NTOTA)
C
!         ENDDO         ! ISYMC
!      ENDDO            ! ISYMK
C
      DO ISYMCK = 1,NSYM
C
          ISYMIJ = MULD2H(ISCKIJ,ISYMCK)
          ISYMAB = MULD2H(ISYMT2,ISYMIJ)
          ISYMA  = MULD2H(ISYMAB,ISYMB)
C
          NTOTCK = MAX(NT1AM(ISYMCK),1)
          NTOTIJ = MAX(NMATIJ(ISYMIJ),1)
          NTOTA  = MAX(NVIR(ISYMA),1)
C
          KOFF1 = ISAIKL(ISYMCK,ISYMIJ) + 1
          KOFF2 = IMAIJA(ISYMIJ,ISYMA) + 1
          KOFF3 = ICKATR(ISYMCK,ISYMA)  + 1
C
          CALL DGEMM('N','N',NT1AM(ISYMCK),NVIR(ISYMA),
     *               NMATIJ(ISYMIJ),TWO,TMAT(KOFF1),NTOTCK,
     *               WORK(KOFF2),NTOTIJ,ONE,RINTE2(KOFF3),
     *               NTOTCK)
C
      ENDDO
C
      CALL QEXIT('CCFOP_DENVIR_SC')
C
      RETURN
      END
C  /* Deck ccfop_denocc_sc */
      SUBROUTINE CCFOP_DENOCC_SC(OCC,SMAT,QMAT,TMAT,ISYMIM,T2AM,ISYMT2,
     *                        WORK,LWORK,INDSQ,LENSQ,ISYMB,B,
     *                        ISYMD,D,IOPT)
C
C     Written by Kasper Hald, Fall 2001.
C     Original version working with symmetry, Sonia Coriani 2002
C
C     Purpose : Calculate the contributions to the t3 and t3-bar
C               densities d_{iajk} and d_{aijk} respectively.
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYMIM, ISYMT2, LWORK,LENSQ, ISYMB, ISYMD, IOPT
      INTEGER INDSQ(LENSQ,6)
      INTEGER ISYMBD, ISELJI, ISYELK, ISYIJK, ISYMK, ISYMEL, ISYMIJ
      INTEGER NTOTEL, NTOTK, KOFF1, KOFF2, KOFF3, ISYML, ISYME
      INTEGER ISYMEK, NTOTIJ, isyres
C
#if defined (SYS_CRAY)
      REAL OCC(*), SMAT(*), QMAT(*), TMAT(*), T2AM(*)
      REAL WORK(LWORK), TWO, ONE, HALF
#else
      DOUBLE PRECISION OCC(*), SMAT(*), QMAT(*), TMAT(*), T2AM(*)
      DOUBLE PRECISION WORK(LWORK), TWO, ONE, HALF
      double precision xnorm, ddot
#endif
C
      PARAMETER (ONE = 1.0D0, TWO = 2.0D0, HALF = 0.5D0)
C
      CALL QENTER('CCFOP_DENOCC_SC')
C
C--------------------------------------
C     Symmetries and sanity check.
C--------------------------------------
C
      IF (IOPT .NE. 1 .AND. IOPT .NE. 2) THEN
         CALL QUIT('Wrong IOPT in CCFOP_DENOCC_SC')
      ENDIF
C
!      ISYMBD = MULD2H(ISYMB,ISYMD)
!      ISYELK = MULD2H(ISYMT2,ISYMD)
!      ISELJI = MULD2H(ISYMIM,ISYMBD)
!      ISYIJK = MULD2H(ISYELK,ISELJI)
!

      ISYRES = MULD2H(ISYMT2,ISYMIM)
      ISYMBD = MULD2H(ISYMB,ISYMD)
      ISYELK = MULD2H(ISYMT2,ISYMD)
      ISELJI = MULD2H(ISYMIM,ISYMBD)
      ISYIJK = MULD2H(ISYRES,ISYMB)
C
      IF (LWORK .LT. NCKI(ISYELK)) THEN
         CALL QUIT('Not enough memory in CCFOP_DENOCC_SC')
      ENDIF
C
C------------------------------------
C     Sort T2 to first term.
C------------------------------------
C
C      DO ISYMK = 1, NSYM
C        ISYMEL = MULD2H(ISYELK,ISYMK)
C        DO ISYML = 1, NSYM
C           ISYME = MULD2H(ISYMEL,ISYML)
CC
C          DO K = 1, NRHF(ISYMK)
C            DO L = 1, NRHF(ISYML)
CC
C               KOFF1 = IT2SP(ISYELK,ISYMD)
C     *                      + NCKI(ISYELK)*(D - 1)
C     *                      + ISAIK(ISYMEL,ISYMK)
C     *                      + NT1AM(ISYMEL)*(K - 1)
C     *                      + IT1AM(ISYME,ISYML)
C     *                      + NVIR(ISYME)*(L-1)
C     *                      + E
CC
C               KOFF2 = ISAIK(ISYMEL,ISYMK)
C     *                      + NT1AM(ISYMEL)*(K - 1)
C     *                      + IT1AM(ISYME,ISYML)
C     *                      + NVIR(ISYME)*(L-1)
C     *                      + E
CC
C               CALL DCOPY(NVIR(ISYME),T2AM(KOFF1),1,WORK(KOFF2),1)
CC
C            ENDDO
C          ENDDO
C        ENDDO
C      ENDDO
C
C--------------------------------------------
C     Contract with S and Q intermediates.
C--------------------------------------------
C
      DO I = 1, NCKIJ(ISELJI)
C
         IF (IOPT .EQ. 1) THEN
            TMAT(I) =       SMAT(INDSQ(I,5))
     *                - TWO*SMAT(I)
     *                +     SMAT(INDSQ(I,3))
     *                +     QMAT(INDSQ(I,4))
     *                - TWO*QMAT(INDSQ(I,3))
     *                +     QMAT(I)
         ELSE
            TMAT(I) = -HALF*SMAT(I)
     *                -HALF*QMAT(INDSQ(I,3))
         ENDIF
C
      ENDDO
!
!Sonia: following bit comes from old code
!
      IF (NSYM .GT. 1) THEN
         IF (LWORK .LT. NCKIJ(ISELJI)) THEN
            CALL QUIT('Out of memory in CCFOP_DENOCC_SC (Gather-1)')
         ENDIF
C
         CALL CC_GATHER(NCKIJ(ISELJI),WORK,TMAT,INDSQ(1,6))
         CALL DCOPY(NCKIJ(ISELJI),WORK,1,TMAT,1)
      ENDIF

!      xnorm = ddot(NCKIJ(ISELJI),TMAT,1,TMAT,1)
!      write(lupri,*) ' TMAT-1 in FOP_DENOC ', xnorm
!
!Sonia: end of it
!
C
      DO ISYMK = 1, NSYM
         ISYMEL = MULD2H(ISYELK,ISYMK)
         ISYMIJ = MULD2H(ISYIJK,ISYMK)
C
         NTOTEL = MAX(NT1AM(ISYMEL),1)
         NTOTK  = MAX(NRHF(ISYMK),1)
C
!
!Sonia: removed following and change to old version. 
!       It does not work with symmetry!
!
!         KOFF1  = IT2SP(ISYELK,ISYMD)
!     *          + NCKI(ISYELK)*(D-1)
!     *          + ISAIK(ISYMEL,ISYMK) + 1
!         KOFF2  = ISAIKL(ISYMEL,ISYMIJ) + 1
!         KOFF3  = I3OVIR(ISYIJK,ISYMB)
!     *          + NMAIJK(ISYIJK)*(B-1)
!     *          + IMAIJK(ISYMIJ,ISYMK)
!     *          + 1
C
!         CALL DGEMM('T','N',NRHF(ISYMK),NMATIJ(ISYMIJ),
!     *              NT1AM(ISYMEL),TWO,T2AM(KOFF1),NTOTEL,
!     *              TMAT(KOFF2),NTOTEL,ONE,OCC(KOFF3),
!     *              NTOTK)

         NTOTIJ = MAX(NMATIJ(ISYMIJ),1)

         KOFF1  = ISAIKL(ISYMEL,ISYMIJ) + 1
         KOFF2  = IT2SP(ISYELK,ISYMD)
     *          + NCKI(ISYELK)*(D-1)
     *          + ISAIK(ISYMEL,ISYMK) + 1
         KOFF3  = I3OVIR(ISYIJK,ISYMB)
     *          + NMAIJK(ISYIJK)*(B-1)
     *          + IMAIJK(ISYMIJ,ISYMK)
     *          + 1
C
         CALL DGEMM('T','N',NMATIJ(ISYMIJ),NRHF(ISYMK),
     *              NT1AM(ISYMEL),TWO,TMAT(KOFF1),NTOTEL,
     *              T2AM(KOFF2),NTOTEL,ONE,OCC(KOFF3),
     *              NTOTIJ)
C
      ENDDO
!      xnorm = ddot(NCKIJ(isyres),OCC,1,OCC,1)
!      write(lupri,*) ' OCC in FOP_DENOC ', xnorm
!      xnorm = ddot(NT2SQ(isyres),T2AM,1,T2AM,1)
!      write(lupri,*) ' T2AM in FOP_DENOC ', xnorm
C
C--------------------------------------------
C     Contract with S and Q intermediates.
C     1. Build second TMAT first
C     2. Contract TMAT with sorted T2
C--------------------------------------------
C
      DO I = 1, NCKIJ(ISELJI)
C
         IF (IOPT .EQ. 1) THEN
            TMAT(I) =       SMAT(INDSQ(I,2))
     *                - TWO*SMAT(INDSQ(I,1))
     *                +     SMAT(INDSQ(I,4))
     *                +     QMAT(INDSQ(I,1))
     *                - TWO*QMAT(INDSQ(I,2))
     *                +     QMAT(INDSQ(I,5))
         ELSE
            TMAT(I) = -HALF*SMAT(INDSQ(I,1))
     *                -HALF*QMAT(INDSQ(I,2))
         ENDIF
C
      ENDDO
!
!Sonia: following bit from old code
!
      IF (NSYM .GT. 1) THEN
         IF (LWORK .LT. NCKIJ(ISELJI)) THEN
            CALL QUIT('Out of memory in CCFOP_DENOCC (Gather-1)')
         ENDIF
C
         CALL CC_GATHER(NCKIJ(ISELJI),WORK,TMAT,INDSQ(1,6))
         CALL DCOPY(NCKIJ(ISELJI),WORK,1,TMAT,1)
      ENDIF
!
!Sonia: end of it
!
C
C------------------------------------
C     Sort T2 to second term.
C------------------------------------
C
      DO ISYMK = 1, NSYM
        ISYMEL = MULD2H(ISYELK,ISYMK)
        DO ISYML = 1, NSYM
           ISYME  = MULD2H(ISYMEL,ISYML)
           ISYMEK = MULD2H(ISYME,ISYMK)
C
          DO K = 1, NRHF(ISYMK)
            DO L = 1, NRHF(ISYML)
C
               KOFF1 = IT2SP(ISYELK,ISYMD)
     *                      + NCKI(ISYELK)*(D - 1)
     *                      + ISAIK(ISYMEL,ISYMK)
     *                      + NT1AM(ISYMEL)*(K - 1)
     *                      + IT1AM(ISYME,ISYML)
     *                      + NVIR(ISYME)*(L-1)
     *                      + 1
C
               KOFF2 = ISAIK(ISYMEK,ISYML)
     *                      + NT1AM(ISYMEK)*(L - 1)
     *                      + IT1AM(ISYME,ISYMK)
     *                      + NVIR(ISYME)*(K-1)
     *                      + 1
C
               CALL DCOPY(NVIR(ISYME),T2AM(KOFF1),1,WORK(KOFF2),1)
C
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C
      DO ISYMK = 1, NSYM
         ISYMEL = MULD2H(ISYELK,ISYMK)
         ISYMIJ = MULD2H(ISYIJK,ISYMK)
C
         NTOTEL = MAX(NT1AM(ISYMEL),1)
         NTOTK  = MAX(NRHF(ISYMK),1)
C
!
!Sonia: following bit does not work with symmetry!
!       Replaced with older version
!
!         KOFF1  = ISAIK(ISYMEL,ISYMK) + 1
!         KOFF2  = ISAIKL(ISYMEL,ISYMIJ) + 1
!         KOFF3  = I3OVIR(ISYIJK,ISYMB)
!     *          + NMAIJK(ISYIJK)*(B-1)
!     *          + IMAIJK(ISYMIJ,ISYMK)
!     *          + 1
C
!         CALL DGEMM('T','N',NRHF(ISYMK),NMATIJ(ISYMIJ),
!     *              NT1AM(ISYMEL),TWO,WORK(KOFF1),NTOTEL,
!     *              TMAT(KOFF2),NTOTEL,ONE,OCC(KOFF3),
!     *              NTOTK)
!
         NTOTIJ = MAX(NMATIJ(ISYMIJ),1)
         KOFF1  = ISAIKL(ISYMEL,ISYMIJ) + 1
         KOFF2  = ISAIK(ISYMEL,ISYMK) + 1
         KOFF3  = I3OVIR(ISYIJK,ISYMB)
     *          + NMAIJK(ISYIJK)*(B-1)
     *          + IMAIJK(ISYMIJ,ISYMK)
     *          + 1
C
         CALL DGEMM('T','N',NMATIJ(ISYMIJ),NRHF(ISYMK),
     *              NT1AM(ISYMEL),TWO,TMAT(KOFF1),NTOTEL,
     *              WORK(KOFF2),NTOTEL,ONE,OCC(KOFF3),
     *              NTOTIJ)
      ENDDO
!      xnorm = ddot(NCKIJ(isyres),OCC,1,OCC,1)
!      write(lupri,*) ' OCC in FOP_DENOC ', xnorm
C
C-----------------------
C     End.
C-----------------------
C
      CALL QEXIT('CCFOP_DENOCC_SC')
C
      RETURN
      END
C  /* Deck den_biasort */
      SUBROUTINE DEN_BIASORT(VIRREAL,VIRTMP,ISYVIR)
C
C     Written by Kasper Hald, 2001.
C
C     Purpose : Sort the two electron densities d(bia) -> d(aib)
C               where the densities have a constant C.
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYVIR, ISYMB, ISYMA, ISYMI, ISYMAI
      INTEGER KOFF1, KOFF2, ISYMBI
COMMENT COMMENT
      integer khcount
COMMENT COMMENT
C
#if defined (SYS_CRAY)
      REAL VIRREAL(*), VIRTMP(*)
#else
      DOUBLE PRECISION VIRREAL(*), VIRTMP(*), tmp
#endif
C
      CALL QENTER('DEN_BIASORT')
C
C----------------------------------------
C     Sort matrix.
C----------------------------------------
C
      DO ISYMB = 1, NSYM
         ISYMAI = MULD2H(ISYMB,ISYVIR)
         DO ISYMA = 1, NSYM
            ISYMI  = MULD2H(ISYMAI,ISYMA)
            ISYMBI = MULD2H(ISYMI,ISYMB)
            DO B = 1, NVIR(ISYMB)
               DO I = 1, NRHF(ISYMI)
C
                  KOFF1 = ICKATR(ISYMAI,ISYMB)
     *                  + NT1AM(ISYMAI)*(B-1)
     *                  + IT1AM(ISYMA,ISYMI)
     *                  + NVIR(ISYMA)*(I-1)
     *                  + 1
C
                  KOFF2 = ICKATR(ISYMBI,ISYMA)
     *                  + IT1AM(ISYMB,ISYMI)
     *                  + NVIR(ISYMB)*(I-1)
     *                  + B
C
                  CALL DCOPY(NVIR(ISYMA),VIRREAL(KOFF1),1,
     *                       VIRTMP(KOFF2),NT1AM(ISYMBI))
C
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C------------------------------------
C     Copy back to original matrix
C------------------------------------
C
      CALL DCOPY(NCKATR(ISYVIR),VIRTMP(1),1,VIRREAL(1),1)
C
C-----------------------
C     End.
C-----------------------
C
      CALL QEXIT('DEN_BIASORT')
C
      RETURN
      END

