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
C  /* Deck ccsdpt_eta */
      SUBROUTINE CCSDPT_ETA(OMEGA1,OMEGA2,T1AM,ISYMT1,T2TP,
     *                      ISYMT2,MODEL,WORK,LWORK,
     *                      LUTOC,FNTOC,LU3VI,FN3VI,LU3VI2,FN3VI2,
     *                      LU3FOP,FN3FOP,LU3FOP2,FN3FOP2)
C
C     Written by K. Hald, Fall 2001.
C
C     Add the triples contribution to the eta vector for solving
C     for the zero order single and double amplitude multipliers.
C
C         ISYMT2 is symmetry of T2TP
C         ISYMT1 is symmetry of T1AM
C
C         Isyres = isymt1*isymt2*isymop
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
#include <second.h>
C
      INTEGER ISYMT1, ISYMT2, LWORK
      INTEGER ISYMTR, ISYRES, ISINT1, ISINT2, ISYMIM, KFOCKD, KOMG1
      INTEGER KOMG22, KCMO, KEND0, LWRK0, KTROC, KTROC1, KTROC2
      INTEGER KTROC0, KXIAJB, KEND1, LWRK1, KINTOC, KEND2, LWRK2
      INTEGER LENGTH, ISYOPE, IOPTTCME, IOFF, ISYMD, ISAIJ1, ISYCKB
      INTEGER ISCKB1, ISCKB2, KTRVI, KTRVI1, KTRVI2, KRMAT1, KTRVI0
      INTEGER KTRVI3, KEND3, LWRK3, KINTVI, KEND4, LWRK4, ISYMB
      INTEGER ISYALJ, ISAIJ2, ISYMBD, ISCKIJ, KSMAT2, KSMAT, KQMAT
      INTEGER KDIAG, ISYMC, ISYMK, KOFF1, KOFF2
      INTEGER KINDSQ, KINDEX, KTMAT, KRMAT2, KRMAT4, LENSQ
      INTEGER LUFCK, KFCKBA, KT2TCME, IOPTT2, KTRVI4, KTRVI5
      INTEGER KTRVI6, KQMAT2, KRMAT3
      INTEGER LUTOC, LU3VI, LU3VI2, LU3FOP, LU3FOP2
C
#if defined (SYS_CRAY)
      REAL OMEGA1(*), OMEGA2(*), T1AM(*), T2TP(*)
      REAL WORK(LWORK), ONE
      REAL TITRAN, TISORT, TISMAT, TIQMAT, TIOME1
      REAL TICONV, TICONO, RHO1N, RHO2N
      REAL XT2TP, DDOT, XIAJB, XINT, XTROC, XTROC1, XTROC0
      REAL XTRVI0, XTRVI2, XTRVI3, XTRVI, XTRVI1, XDIA
      REAL XSMAT, XTMAT, XQMAT, XRMAT, TWO, ZERO, HALF
      REAL DTIME
#else
      DOUBLE PRECISION OMEGA1(*), OMEGA2(*), T1AM(*), T2TP(*)
      DOUBLE PRECISION WORK(LWORK), ONE
      DOUBLE PRECISION TITRAN, TISORT, TISMAT, TIQMAT, TIOME1
      DOUBLE PRECISION TICONV, TICONO, RHO1N, RHO2N
      DOUBLE PRECISION XT2TP, DDOT, XIAJB, XINT, XTROC, XTROC1, XTROC0
      DOUBLE PRECISION XTRVI0, XTRVI2, XTRVI3, XTRVI, XTRVI1, XDIA
      DOUBLE PRECISION XSMAT, XTMAT, XQMAT, XRMAT, TWO, ZERO, HALF
      DOUBLE PRECISION DTIME
#endif
C
      LOGICAL   C3LRSV, CC1ASV, CC1BSV
      CHARACTER*10 MODEL
      CHARACTER*(*) FNTOC, FN3VI, FN3VI2, FN3FOP, FN3FOP2
C
      PARAMETER (HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0, ZERO = 0.0D0)
C
      CALL QENTER('CCSDPT_ETA')
C
C--------------------------
C     Save and set flags.
C--------------------------
C
      C3LRSV = CC3LR
      CC1ASV = CC1A
      CC1BSV = CC1B
      CC3LR  = .FALSE.
      CC1A   = .TRUE.
      CC1B   = .TRUE.
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
      TICONO = 0.0D0
      TICONV = 0.0D0
      TIOME1 = 0.0D0
C
C--------------------------------------
C     Reorder the t2-amplitudes i T2TP.
C--------------------------------------
C
      IF (LWORK .LT. NT2SQ(ISYMT2)) THEN
         CALL QUIT('Not enough memory to construct T2TP in CC3_OMEG')
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
C---------------------------------------------------------
C     Read canonical orbital energies and MO coefficients.
C---------------------------------------------------------
C
      KFOCKD = 1
      KOMG1  = KFOCKD + NORBTS
      KOMG22 = KOMG1  + NT1AM(ISYRES)
      KCMO   = KOMG22 + NT2AM(ISYRES)
      KFCKBA = KCMO   + NLAMDS
      KEND0  = KFCKBA + N2BST(ISYMOP)
      LWRK0  = LWORK  - KEND0
C
      CALL DZERO(WORK(KOMG1),NT1AM(ISYRES))
      CALL DZERO(WORK(KOMG22),NT2AM(ISYRES))
C
      IF (LWRK0 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    > ',KEND0
         CALL QUIT('Insufficient space in CCSDT_OMEG')
      END IF
C
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
C
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KFOCKD+I-1), I=1,NORBTS)
      READ (LUSIFC) (WORK(KCMO+I-1), I=1,NLAMDS)
C
      CALL GPCLOSE(LUSIFC,'KEEP')
C
      CALL CMO_REORDER(WORK(KCMO),WORK(KEND0),LWRK0)
C
C---------------------------------------------
C     Delete frozen orbitals in Fock diagonal.
C---------------------------------------------
C
      IF (FROIMP .OR. FROEXP)
     *   CALL CCSD_DELFRO(WORK(KFOCKD),WORK(KEND0),LWRK0)
C
C-----------------------------------------------------
C     Construct the CMO transformed Fock matrix
C-----------------------------------------------------
C
        LUFCK = -1
C       This AO Fock matrix is constructed from the T1 transformed density
C        CALL GPOPEN(LUFCK,'CC_FCKH','UNKNOWN',' ','UNFORMATTED',
C     *              IDUMMY,.FALSE.)
C       This AO Fock matrix is constructed from the CMO transformed density
        CALL GPOPEN(LUFCK,'CC_FCKREF','UNKNOWN',' ','UNFORMATTED',
     *              IDUMMY,.FALSE.)
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
        CALL CC_FCKMO(WORK(KFCKBA),WORK(KCMO),WORK(KCMO),
     *                WORK(KEND0),LWRK0,1,1,1)
C
        IF (IPRINT .GT. 50) THEN
           CALL AROUND( 'In CC_ETA: Triples Fock MO matrix' )
           CALL CC_PRFCKMO(WORK(KFCKBA),ISYMOP)
        ENDIF
C
C        Sort the fock matrix
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
     *                    NORB(ISYMK)*(C - 1) + K - 1
                  KOFF2 = KFCKBA + IT1AM(ISYMC,ISYMK)
     *                  + NVIR(ISYMC)*(K - 1) + C - 1
C
                  WORK(KOFF2) = WORK(KOFF1)
C
               ENDDO
            ENDDO
         ENDDO
C
        IF (IPRINT .GT. 50) THEN
           CALL AROUND( 'In CC_ETA: Triples Fock MO matrix (sort)' )
           CALL CC_PRFCKMO(WORK(KFCKBA),ISYMOP)
        ENDIF
C
C------------------------------------------------------------------
C     Read in another T2 amplitude, and transform it to 2*C-E
C     Square up to full matrix and reorder the indexing
C------------------------------------------------------------------
C
      KT2TCME = KEND0
      KEND0   = KT2TCME + NT2SQ(ISYMT2)
      LWRK0   = LWORK - KEND0
C
      IF (LWRK0 .LT. NT2SQ(1)) 
     *     CALL QUIT('Too litlle workspace CCSDPT_ETA T2')
C
      IOPTT2 = 2
      CALL CC_RDRSP('R0',0,1,IOPTT2,MODEL,DUMMY,WORK(KEND0))
C
      ISYOPE = ISYMT2
      IOPTT2 = 1
      CALL CCSD_TCMEPK(WORK(KEND0),1.0D0,ISYOPE,IOPTT2)
C
      CALL CC_T2SQ(WORK(KEND0),WORK(KT2TCME),ISYMT2)
C
      CALL DCOPY(NT2SQ(ISYMT2),WORK(KT2TCME),1,WORK(KEND0),1)
      CALL CC3_T2TP(WORK(KT2TCME),WORK(KEND0),1)
C
      IF (IPRINT .GT. 55) THEN
         XT2TP = DDOT(NT2SQ(ISYMT2),WORK(KT2TCME),1,WORK(KT2TCME),1)
         WRITE(LUPRI,*) 'Norm of 2*C-E T2 amplitudes after resort ',
     *                    XT2TP
      ENDIF
C
C-----------------------------
C     Read occupied integrals.
C-----------------------------
C
C     Memory allocation.
C
      KTROC  = KEND0
      KTROC1 = KTROC  + NTRAOC(ISINT1)
      KTROC0 = KTROC1 + NTRAOC(ISINT1)
      KTROC2 = KTROC0 + NTRAOC(ISINT2)
      KXIAJB = KTROC2 + NTRAOC(ISINT2)
      KEND1  = KXIAJB + NT2AM(ISYMOP)
      LWRK1  = LWORK  - KEND1
C
      KINTOC = KEND1
      KEND2  = KINTOC + MAX(NTOTOC(ISYMOP),NTOTOC(ISINT2))
      LWRK2  = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    > ',KEND2
         CALL QUIT('Insufficient space in CCSDT_OMEG')
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
         XIAJB = DDOT(NT2AM(ISYMOP),WORK(KXIAJB),1,
     *                WORK(KXIAJB),1)
         WRITE(LUPRI,*) 'Norm of IAJB ',XIAJB
      ENDIF
C
C------------------------
C     Occupied integrals.
C------------------------
C
      DTIME = SECOND()
      IOFF = 1
      IF (NTOTOC(ISYMOP) .GT. 0) THEN
         CALL GETWA2(LUTOC,FNTOC,WORK(KINTOC),IOFF,NTOTOC(ISYMOP))
      ENDIF
C
C----------------------------------
C     Write out norms of Integrals.
C----------------------------------
C
      IF (IPRINT .GT. 55) THEN
         XINT  = DDOT(NTOTOC(ISYMOP),WORK(KINTOC),1,
     *                WORK(KINTOC),1)
         WRITE(LUPRI,*) 'Norm of CCSDT_OC-INT ',XINT
      ENDIF
C
C----------------------------------------------------------------------
C     Transform (ia|j delta) integrals to (ia|j k) and sort as (ij,k,a)
C----------------------------------------------------------------------
C
      DTIME  = SECOND() - DTIME
C
      CALL CCSDT_TROCC(WORK(KINTOC),WORK(KTROC),WORK(KCMO),
     *                 WORK(KEND2),LWRK2)
C
C----------------------------------------------------------------------
C     Transform (ia|j delta) integrals to (ia|j k) and sort as (ij,k,a)
C----------------------------------------------------------------------
C
      DTIME = SECOND()
      CALL CCSDT_TROCC(WORK(KINTOC),WORK(KTROC0),WORK(KCMO),
     *                 WORK(KEND2),LWRK2)
C
      DTIME  = SECOND() - DTIME
      TITRAN = TITRAN   + DTIME

C
      DTIME = SECOND()
C
      CALL CCSDT_SRTOC2(WORK(KTROC),WORK(KTROC1),ISINT1,
     *                  WORK(KEND2),LWRK2)
C
      DTIME  = SECOND() - DTIME
      TISORT = TISORT   + DTIME
C
C-----------------------------------------------------------
C     Construct 2*C-E for the occupied integrals.
C-----------------------------------------------------------
C
      CALL CCSDT_TCMEOCC(WORK(KTROC0),WORK(KTROC2),ISINT2)
C
C-------------------------------
C     Write out norms of arrays.
C-------------------------------
C
      IF (IPRINT .GT. 55) THEN
         XTROC = DDOT(NTRAOC(ISINT1),WORK(KTROC),1,
     *                WORK(KTROC),1)
         WRITE(LUPRI,*) 'Norm of TROC ',XTROC
      ENDIF
C
      IF (IPRINT .GT. 55) THEN
         XINT  = DDOT(NTOTOC(ISINT2),WORK(KINTOC),1,
     *                WORK(KINTOC),1)
         WRITE(LUPRI,*) 'Norm of CKJDEL-INT  ',XINT
      ENDIF
C
      IF (IPRINT .GT. 55) THEN
         XTROC1 = DDOT(NTRAOC(ISINT1),WORK(KTROC1),1,
     *                WORK(KTROC1),1)
         WRITE(LUPRI,*) 'Norm of TROC1 ',XTROC1
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
            WRITE(LUPRI,*) 'In CCSDPT_ETA: ISAIJ1:',ISAIJ1
            WRITE(LUPRI,*) 'In CCSDPT_ETA: ISYCKB:',ISYCKB
            WRITE(LUPRI,*) 'In CCSDPT_ETA: ISCKB1:',ISCKB1
            WRITE(LUPRI,*) 'In CCSDPT_ETA: ISCKB2:',ISCKB2
C
         ENDIF
C
C--------------------------
C        Memory allocation.
C--------------------------
C
         KTRVI  = KEND1
         KTRVI1 = KTRVI  + NCKATR(ISCKB1)
         KTRVI2 = KTRVI1 + NCKATR(ISCKB1)
         KRMAT1 = KTRVI2 + NCKATR(ISCKB2)
         KRMAT3 = KRMAT1 + NCKI(ISAIJ1)
         KEND2  = KRMAT3 + NCKI(ISAIJ1)
         LWRK2  = LWORK  - KEND2
C
         KTRVI0 = KEND2
         KTRVI3 = KTRVI0 + NCKATR(ISCKB2)
         KTRVI4 = KTRVI3 + NCKATR(ISCKB2)
         KTRVI5 = KTRVI4 + NCKATR(ISCKB2)
         KTRVI6 = KTRVI5 + NCKATR(ISCKB2)
         KEND3  = KTRVI6 + NCKATR(ISCKB2)
         LWRK3  = LWORK  - KEND3
C
         KINTVI = KEND3
         KEND4  = KINTVI + MAX(NCKA(ISYCKB),NCKA(ISYMD),NCKA(ISCKB2))
         LWRK4  = LWORK  - KEND4
C
         IF (LWRK4 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    > ',KEND4
            CALL QUIT('Insufficient space in CCSDT_OMEG')
         END IF
C
         DO D = 1,NVIR(ISYMD)
C
C------------------------------------
C           Initialize the R1 matrix.
C------------------------------------
C
            CALL DZERO(WORK(KRMAT1),NCKI(ISAIJ1))
            CALL DZERO(WORK(KRMAT3),NCKI(ISAIJ1))
C
C------------------------------------------------------------
C           Read and transform integrals used in contraction.
C------------------------------------------------------------
C
            IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D - 1) + 1
            IF (NCKA(ISYCKB) .GT. 0) THEN
               CALL GETWA2(LU3VI2,FN3VI2,WORK(KINTVI),IOFF,
     &                     NCKA(ISYCKB))
            ENDIF
C
            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI),WORK(KCMO),
     *                       ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
C
C-------------------------------------------------------
C           Calculate virtual integrals used in s3am.
C-------------------------------------------------------
C
            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI0),WORK(KCMO),
     *                       ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
C
C------------------------------------------------------
C           Read 2*C-E of integral used for t3-bar
C------------------------------------------------------
C
            IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D - 1) + 1
            IF (NCKA(ISYCKB) .GT. 0) THEN
               CALL GETWA2(LU3FOP2,FN3FOP2,WORK(KINTVI),IOFF,
     &                     NCKA(ISYCKB))
            ENDIF
C
            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI4),WORK(KCMO),
     *                       ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
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
            IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D - 1) + 1
            IF (NCKA(ISYCKB) .GT. 0) THEN
               CALL GETWA2(LU3VI,FN3VI,WORK(KINTVI),IOFF,
     &                     NCKA(ISYCKB))
            ENDIF
C
            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI1),WORK(KCMO),
     *                       ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
C
C--------------------------------------------------------
C           Calculate virtual integrals used in q3am.
C--------------------------------------------------------
C
C
            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI3),WORK(KCMO),
     *                       ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
C
            IF (LWRK3 .LT. NCKATR(ISYCKB)) THEN
               CALL QUIT('Insufficient space for allocation in '//
     &                   'CCSDPT_ETA (1)')
            END IF
C
C           Can use kend3 since dont need the integrals anymore
            DTIME = SECOND()
            CALL CCSDT_SRVIR3(WORK(KTRVI3),WORK(KEND3),ISYMD,D,ISINT2)
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
C           Read virtual integrals used in q3am for t3-bar.
C---------------------------------------------------------------
C
            IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D - 1) + 1
            IF (NCKA(ISYCKB) .GT. 0) THEN
               CALL GETWA2(LU3FOP,FN3FOP,WORK(KINTVI),IOFF,
     &                     NCKA(ISYCKB))
            ENDIF
C
            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI6),WORK(KCMO),
     *                       ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
C
            IF (LWRK3 .LT. NCKATR(ISYCKB)) THEN
               CALL QUIT('Insufficient space for allocation in '//
     &                   'CCSDPT_ETA (2)')
            END IF
C
C           Can use kend3 since dont need the integrals anymore
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
            ENDIF
C
C---------------------------------------------
C           Construct integrals used in CC3LR.
C---------------------------------------------
C
            IF (IPRINT .GT. 55) THEN
               XTRVI= DDOT(NCKATR(ISCKB1),WORK(KTRVI),1,
     *                      WORK(KTRVI),1)
               WRITE(LUPRI,*) 'Norm of TRVI ',XTRVI
            ENDIF
C
            IF (IPRINT .GT. 55) THEN
               XTRVI1= DDOT(NCKATR(ISCKB1),WORK(KTRVI1),1,
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
               ISYALJ = MULD2H(ISYMB,ISYMT2)
               ISAIJ2 = MULD2H(ISYMB,ISYRES)
               ISYMBD = MULD2H(ISYMB,ISYMD)
               ISCKIJ = MULD2H(ISYMBD,ISYMIM)
C
               IF (IPRINT .GT. 55) THEN
C
                  WRITE(LUPRI,*) 'In CC3_OMEG: ISYMD :',ISYMD
                  WRITE(LUPRI,*) 'In CC3_OMEG: ISYMB :',ISYMB
                  WRITE(LUPRI,*) 'In CC3_OMEG: ISYALJ:',ISYALJ
                  WRITE(LUPRI,*) 'In CC3_OMEG: ISAIJ2:',ISAIJ2
                  WRITE(LUPRI,*) 'In CC3_OMEG: ISYMBD:',ISYMBD
                  WRITE(LUPRI,*) 'In CC3_OMEG: ISCKIJ:',ISCKIJ
C
               ENDIF
C
               KSMAT  = KEND3
               KQMAT  = KSMAT  + NCKIJ(ISCKIJ)
               KSMAT2 = KQMAT  + NCKIJ(ISCKIJ)
               KQMAT2 = KSMAT2 + NCKIJ(ISCKIJ)
               KDIAG  = KQMAT2 + NCKIJ(ISCKIJ)
               KINDSQ = KDIAG  + NCKIJ(ISCKIJ)
               KINDEX = KINDSQ + (6*NCKIJ(ISCKIJ) - 1)/IRAT + 1
               KTMAT  = KINDEX + (NCKI(ISYALJ) - 1)/IRAT + 1
               KRMAT2 = KTMAT  + NCKIJ(ISCKIJ)
               KRMAT4 = KRMAT2 + NCKI(ISAIJ2)
               KEND4  = KRMAT4 + NCKI(ISAIJ2)
               LWRK4  = LWORK  - KEND4
C
               IF (LWRK4 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Memory available : ',LWORK
                  WRITE(LUPRI,*) 'Memory needed    > ',KEND4
                  CALL QUIT('Insufficient space in CCSDT_OMEG')
               END IF
C
C---------------------------------------------
C              Construct part of the diagonal.
C---------------------------------------------
C
               CALL CC3_DIAG(WORK(KDIAG),WORK(KFOCKD),ISCKIJ)
C
               IF (IPRINT .GT. 55) THEN
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
C
               DO B = 1,NVIR(ISYMB)
C
C-----------------------------------------
C                 Initialize the R2 matrix.
C-----------------------------------------
C
                  CALL DZERO(WORK(KRMAT2),NCKI(ISAIJ2))
                  CALL DZERO(WORK(KRMAT4),NCKI(ISAIJ2))
C
C-------------------------------------------------------------
C                 Calculate the S(ci,bk,dj) matrix for T3.
C-------------------------------------------------------------
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
                  DTIME  = SECOND() - DTIME
                  TISMAT = TISMAT   + DTIME
C
                  IF (IPRINT .GT. 55) THEN
                     XSMAT = DDOT(NCKIJ(ISCKIJ),WORK(KSMAT),1,
     *                       WORK(KSMAT),1)
                     WRITE(LUPRI,*) 'Norm of SMAT  ',XSMAT
                  ENDIF
C
                  IF (IPRINT .GT. 55) THEN
                     XTMAT = DDOT(NCKIJ(ISCKIJ),WORK(KTMAT),1,
     *                       WORK(KTMAT),1)
                     WRITE(LUPRI,*) 'Norm of TMAT  ',XTMAT
                  ENDIF
C
C---------------------------------------------------------------
C                 Calculate the S(ci,bk,dj) matrix for T3-BAR.
C---------------------------------------------------------------
C
                  DTIME = SECOND()
C
                  CALL DZERO(WORK(KSMAT2),NCKIJ(ISCKIJ))
C
                  CALL CCFOP_SMAT(0.0D0,T1AM,ISYMT1,WORK(KT2TCME),
     *                            ISYMT2,WORK(KTMAT),
     *                            WORK(KFCKBA),WORK(KXIAJB),ISINT1,
     *                            WORK(KTRVI0),WORK(KTRVI2),
     *                            WORK(KTRVI4),WORK(KTRVI5),
     *                            WORK(KTROC0),WORK(KTROC2),
     *                            ISINT2,WORK(KFOCKD),
     *                            WORK(KDIAG),WORK(KSMAT2),WORK(KEND4),
     *                            LWRK4,WORK(KINDEX),WORK(KINDSQ),LENSQ,
     *                            ISYMB,B,ISYMD,D)
C
                  DTIME  = SECOND() - DTIME
                  TISMAT = TISMAT   + DTIME
C
                  IF (IPRINT .GT. 55) THEN
                     XSMAT = DDOT(NCKIJ(ISCKIJ),WORK(KSMAT2),1,
     *                       WORK(KSMAT2),1)
                     WRITE(LUPRI,*) 'Norm of SMAT2 ',XSMAT
                  ENDIF
C
                  IF (IPRINT .GT. 55) THEN
                     XTMAT = DDOT(NCKIJ(ISCKIJ),WORK(KTMAT),1,
     *                       WORK(KTMAT),1)
                     WRITE(LUPRI,*) 'Norm of TMAT  ',XTMAT
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
                  DTIME  = SECOND() - DTIME
                  TIQMAT = TIQMAT   + DTIME
C
                  IF (IPRINT .GT. 55) THEN
                     XQMAT = DDOT(NCKIJ(ISCKIJ),WORK(KQMAT),1,
     *                       WORK(KQMAT),1)
                     WRITE(LUPRI,*) 'Norm of QMAT  ',XQMAT
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
                  CALL CCFOP_QMAT(0.0D0,T1AM,ISYMT1,WORK(KT2TCME),
     *                            ISYMT2,WORK(KTMAT),WORK(KFCKBA),
     *                            WORK(KXIAJB),ISINT1,WORK(KTRVI3),
     *                            WORK(KTRVI6),WORK(KTROC0),
     *                            WORK(KTROC2),ISINT2,WORK(KFOCKD),
     *                            WORK(KDIAG),WORK(KQMAT2),
     *                            WORK(KEND4),LWRK4,WORK(KINDSQ),
     *                            LENSQ,ISYMB,B,ISYMD,D)
C
                  DTIME  = SECOND() - DTIME
                  TIQMAT = TIQMAT   + DTIME
C
                  IF (IPRINT .GT. 55) THEN
                     XQMAT = DDOT(NCKIJ(ISCKIJ),WORK(KQMAT),1,
     *                       WORK(KQMAT),1)
                     WRITE(LUPRI,*) 'Norm of QMAT  ',XQMAT
                  ENDIF
C
C-----------------------------------------
C                 Contract with integrals.
C-----------------------------------------
C
                  DTIME = SECOND()
                  CALL CC3_CONVIR(WORK(KRMAT2),WORK(KSMAT),
     *                            WORK(KQMAT),WORK(KTMAT),ISYMIM,
     *                            WORK(KTRVI),WORK(KTRVI1),ISINT1,
     *                            WORK(KEND4),LWRK4,
     *                            WORK(KINDSQ),LENSQ,ISYMB,B,ISYMD,D)
                  IF (IPRINT .GT. 55) THEN
                     XRMAT = DDOT(NCKI(ISAIJ2),WORK(KRMAT2),1,
     *                       WORK(KRMAT2),1)
                     WRITE(LUPRI,*) 'Norm RMAT2 - CC3_CONVIR',XRMAT
                  ENDIF
C
                  CALL CCFOP_CONVIR(WORK(KRMAT4),WORK(KSMAT2),
     *                              WORK(KQMAT2),WORK(KTMAT),ISYMIM,
     *                              WORK(KTRVI),WORK(KTRVI1),ISINT1,
     *                              WORK(KEND4),LWRK4,WORK(KINDSQ),
     *                              LENSQ,ISYMB,B,ISYMD,D)
C
                  IF (IPRINT .GT. 55) THEN
                     XRMAT = DDOT(NCKI(ISAIJ2),WORK(KRMAT4),1,
     *                       WORK(KRMAT4),1)
                     WRITE(LUPRI,*) 'Norm RMAT4 - CCFOP_CONVIR ',XRMAT
                  ENDIF
C
                  DTIME  = SECOND() - DTIME
                  TICONV = TICONV   + DTIME
C
                  DTIME = SECOND()
                  CALL CC3_CONOCC(WORK(KOMG22),WORK(KRMAT1),WORK(KRMAT2)
     *                           ,WORK(KSMAT),WORK(KTMAT),ISYMIM,
     *                            WORK(KTROC),WORK(KTROC1),ISINT1,
     *                            WORK(KEND4),LWRK4,
     *                            WORK(KINDSQ),LENSQ,ISYMB,B,ISYMD,D)
C
                  IF (IPRINT .GT. 55) THEN
                     XRMAT = DDOT(NCKI(ISAIJ1),WORK(KRMAT1),1,
     *                       WORK(KRMAT1),1)
                     WRITE(LUPRI,*) 'Norm of RMAT1 - CC3_CONOCC ',XRMAT
                     XRMAT = DDOT(NCKI(ISAIJ2),WORK(KRMAT2),1,
     *                       WORK(KRMAT2),1)
                     WRITE(LUPRI,*) 'Norm of RMAT2 - CC3_CONOCC',XRMAT
                     RHO2N = DDOT(NT2AM(ISYRES),WORK(KOMG22),1,
     *                                          WORK(KOMG22),1)
                     WRITE(LUPRI,*) 'Norm of Rho2 -after CC3_CONOCC',
     *                                                 RHO2N
                  ENDIF
C
                  IF (IPRINT .GT. 220) THEN
                     CALL AROUND('After CC3_CONOCC: ')
                     CALL CC_PRP(OMEGA1,WORK(KOMG22),ISYRES,0,1)
                  ENDIF
C
                  DTIME = SECOND()
                  CALL CCFOP_CONOCC(OMEGA2,WORK(KRMAT3),
     *                              WORK(KRMAT4),WORK(KSMAT2),
     *                              WORK(KTMAT),ISYMIM,
     *                              WORK(KTROC),WORK(KTROC1),ISINT1,
     *                              WORK(KEND4),LWRK4,WORK(KINDSQ),
     *                              LENSQ,ISYMB,B,ISYMD,D)
C
                  IF (IPRINT .GT. 55) THEN
                     XRMAT = DDOT(NCKI(ISAIJ1),WORK(KRMAT3),1,
     *                       WORK(KRMAT3),1)
                     WRITE(LUPRI,*) 'Norm RMAT3 - CCFOP_CONOCC ',XRMAT
                     XRMAT = DDOT(NCKI(ISAIJ2),WORK(KRMAT4),1,
     *                       WORK(KRMAT4),1)
                     WRITE(LUPRI,*) 'Norm RMAT4 - CCFOP_CONOCC',XRMAT
                     RHO2N = DDOT(NT2AM(ISYRES),WORK(KOMG22),1,
     *                                          WORK(KOMG22),1)
                     WRITE(LUPRI,*) 'Norm of Rho2 -after CCFOP_CONOCC',
     *                                                 RHO2N
                  ENDIF
C
                  IF (IPRINT .GT. 220) THEN
                     CALL AROUND('After CCFOP_CONOCC: ')
                     CALL CC_PRP(OMEGA1,OMEGA2,ISYRES,0,1)
                  ENDIF
C
                  DTIME  = SECOND() - DTIME
                  TICONO = TICONO   + DTIME
C
C----------------------------------
C                 Calculate Omega1.
C----------------------------------
C
                  DTIME = SECOND()
C
                  CALL CC3_ONEL(WORK(KOMG1),WORK(KOMG22),WORK(KRMAT1),
     *                          WORK(KRMAT2),WORK(KFCKBA),WORK(KSMAT),
     *                          WORK(KTMAT),ISYMIM,WORK(KXIAJB),ISINT1,
     *                          WORK(KINDSQ),LENSQ,WORK(KEND4),LWRK4,
     *                          ISYMB,B,ISYMD,D)
C
                  IF (IPRINT .GT. 55) THEN
                     RHO1N = DDOT(NT1AM(ISYRES),OMEGA1,1,OMEGA1,1)
                     RHO2N = DDOT(NT2AM(ISYRES),WORK(KOMG22),1,
     *                                          WORK(KOMG22),1)
                     WRITE(LUPRI,*) 'Norm of Rho1 -after CC3_ONEL',RHO1N
                     WRITE(LUPRI,*) 'Norm of Rho2 -after CC3_ONEL',RHO2N
                  ENDIF
C
                  IF (IPRINT .GT. 220) THEN
                     CALL AROUND('After CC3_ONEL: ')
                     CALL CC_PRP(OMEGA1,WORK(KOMG22),ISYRES,1,1)
                  ENDIF
C
                  IF (IPRINT .GT. 55) THEN
                     XRMAT = DDOT(NCKI(ISAIJ1),WORK(KRMAT1),1,
     *                       WORK(KRMAT1),1)
                     WRITE(LUPRI,*) 'Norm of RMAT1 -after ONEL',XRMAT
                     XRMAT = DDOT(NCKI(ISAIJ2),WORK(KRMAT2),1,
     *                       WORK(KRMAT2),1)
                     WRITE(LUPRI,*) 'Norm of RMAT2 -after ONEL',XRMAT
                  ENDIF
C
C
                  DTIME  = SECOND() - DTIME
                  TIOME1 = TIOME1   + DTIME
C
C---------------------------------------------------------
C                 Accumulate the R2 matrix in Omg22 and
C                                R4        in OMEGA2.
C---------------------------------------------------------
C
                  CALL CC3_RACC(WORK(KOMG22),WORK(KRMAT2),ISYMB,B,
     *                          ISYRES)
                  CALL CC3_RACC(OMEGA2,WORK(KRMAT4),ISYMB,B,
     *                          ISYRES)
C
                  IF (IPRINT .GT. 55) THEN
                     RHO1N = DDOT(NT1AM(ISYRES),OMEGA1,1,OMEGA1,1)
                     RHO2N = DDOT(NT2AM(ISYRES),OMEGA2,1,OMEGA2,1)
                     WRITE(LUPRI,*) 'Norm of Rho1 -after CC3_RACC',RHO1N
                     WRITE(LUPRI,*) 'Norm of Rho2 -after CC3_RACC',RHO2N
                     RHO2N = DDOT(NT2AM(ISYRES),WORK(KOMG22),1,
     *                                          WORK(KOMG22),1)
                     WRITE(LUPRI,*) 'Norm of Rho22-after CC3_RACC',RHO2N
                  ENDIF
C
                  IF (IPRINT .GT. 220) THEN
                     CALL AROUND('After CC3_RACC: ')
                     CALL CC_PRP(OMEGA1,WORK(KOMG22),ISYRES,1,1)
                  ENDIF
C
               ENDDO        ! B
            ENDDO           ! ISYMB
C
C---------------------------------------------------
C           Accumulate the R1 matrix in Omega22 and
C                          R3        in Omega2.
C---------------------------------------------------
C
            CALL CC3_RACC(WORK(KOMG22),WORK(KRMAT1),ISYMD,D,ISYRES)
            CALL CC3_RACC(OMEGA2,WORK(KRMAT3),ISYMD,D,ISYRES)
C
            IF (IPRINT .GT. 55) THEN
               RHO1N = DDOT(NT1AM(ISYRES),OMEGA1,1,OMEGA1,1)
               RHO2N = DDOT(NT2AM(ISYRES),OMEGA2,1,OMEGA2,1)
               WRITE(LUPRI,*) 'Norm of Rho1 -after CC3_RACC-2',RHO1N
               WRITE(LUPRI,*) 'Norm of Rho2 -after CC3_RACC-2',RHO2N
               RHO2N = DDOT(NT2AM(ISYRES),WORK(KOMG22),1,WORK(KOMG22),1)
               WRITE(LUPRI,*) 'Norm of Rho22-after CC3_RACC-2',RHO2N
            ENDIF
C
            IF (IPRINT .GT. 220) THEN
               CALL AROUND('After CC3_RACC-2: ')
               CALL CC_PRP(OMEGA1,WORK(KOMG22),ISYRES,1,1)
            ENDIF
C
         ENDDO     ! D
      ENDDO        ! ISYMD
C
C-----------------------------------------------
C     Scale the one-electron result vector with
C     two since t1* = 2 t1 and add to real
C     result vector.
C-----------------------------------------------
C
      CALL DAXPY(NT1AM(ISYRES),TWO,WORK(KOMG1),1,OMEGA1,1)
C
C----------------------------------------------------------------
C     Take two times (two coulomb minus exchange) in
C     the double result vector from the normal T3 amplitudes 
C     and sum up in Omega2
C----------------------------------------------------------------
C
      CALL DSCAL(NT2AM(ISYRES),TWO,WORK(KOMG22),1)
C
      IOPTTCME = 1
      ISYOPE   = ISYRES
      CALL CCSD_TCMEPK(WORK(KOMG22),1.0D0,ISYOPE,IOPTTCME)
C
      CALL DAXPY(NT2AM(ISYRES),ONE,WORK(KOMG22),1,OMEGA2,1)
C
      IF (IPRINT .GT. 110) THEN
         CALL AROUND('Omega1 and Omega2 at the end of CCSDPT_ETA')
         CALL CC_PRP(OMEGA1,OMEGA2,ISYRES,1,1)
      ENDIF
C
C-----------------------
C     Restore flags.
C-----------------------
C
      CC3LR = C3LRSV
      CC1A  = CC1ASV
      CC1B  = CC1BSV
C
C-------------------
C     Print timings.
C-------------------
C
      IF (IPRINT .GT. 9) THEN
COMMENT COMMENT
COMMENT COMMENT  Have a look at the timings in this routine.
COMMENT COMMENT
         WRITE(LUPRI,*)
         WRITE(LUPRI,*)
         WRITE(LUPRI,1) 'CC3_TRAN  : ',TITRAN
         WRITE(LUPRI,1) 'CC3_SORT  : ',TISORT
         WRITE(LUPRI,1) 'CC3_SMAT  : ',TISMAT
         WRITE(LUPRI,1) 'CC3_QMAT  : ',TIQMAT
         WRITE(LUPRI,1) 'CC3_CONV  : ',TICONV
         WRITE(LUPRI,1) 'CC3_CONO  : ',TICONO
         WRITE(LUPRI,1) 'CC3_OME1  : ',TIOME1
         WRITE(LUPRI,*)
      END IF
C
C-------------
C     End
C-------------
C
      CALL QEXIT('CCSDPT_ETA')
C
      RETURN
C
    1 FORMAT(7X,'Time used in',2X,A12,F12.2,' seconds')
C
      END
C  /* Deck ccfop_smat */
      SUBROUTINE CCFOP_SMAT(ECURR,T1AM,ISYMT1,T2TCME,ISYMT2,TMAT,FOCK,
     *                      XIAJB,
     *                      ISINT1,TRVIR,TRVIR2,TRVIR4,TRVIR5,TROCC,
     *                      TROCC2,ISINT2,FOCKD,DIAG,SMAT,WORK,LWORK,
     *                      INDAJL,INDSQ,LENSQ,ISYMB,B,ISYMC,C)
C
C     Written by Kasper Hald, Fall 2001.
C
C     Calculate the S matrix for the t3-bar coefficients
C     in the CCSD(T) model.
C
C     S is stored as S(ai,k,j) for fixed b and c
C     (kc|bd) is stored as I(dk,b,c)
C
C     T2TCME is two times coulomb minus exchange of the T2 ampl.
C     IN TRVIR,  TRVIR2 and TROCC  are the normal integrals
C     IN TRVIR4, TRVIR5 and TROCC2 are the integrals transformed to 2*C-E
C
C     General symmetry: ISINT1 / ISINT2 is symmetry of integrals 
C                       ISYMT1 / ISYMT2 is symmetry of T1AM / (T2TP,T2TCME).
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INTEGER ISYMT1, ISYMT2, ISINT1, LWORK, LENSQ, ISYMB, ISYMC
      INTEGER ISINT2, ISYMBC, ISYRES, JSAIKJ, ISYMDK, LENGTH, ISYMK
      INTEGER ISYMD, ISYAIJ, KOFF1, KOFF2, KOFF3, NTOAIJ, NVIRD, ISYAKD
      INTEGER ISYDIJ, ISYMJ, ISYMDI, ISYMI, ISYMAK, ISYAKI, NTOTAK
      INTEGER ISYAIL, ISYLKJ, ISYMLK, ISYML, ISYMAI, ISYAIK, NTOTAI
      INTEGER NRHFL, ISYAJL, ISYLKI, KOFF, ISYMAJ, ISYAJK, NTOTAJ
      INTEGER NB, NC, ISYMA, ISYRES2
      INTEGER ISYMJK, ISYMCK, NBJ, NKJ, NCK, NCKBJ, ISYMBJ
      INTEGER NAI, NAIKJ, NAIK, NAICK, NAIBJ, NCKTEMP, NAIKJTEMP
      INTEGER ISYMBI, NBI, NAJ, NAJBI
      INTEGER INDEX, INDAJL, INDSQ(LENSQ,6)
C
#if defined (SYS_CRAY)
      REAL T1AM(*), TMAT(*), XIAJB(*), FOCK(*)
      REAL TRVIR(*), TRVIR2(*), TROCC(*), FOCKD(*), DIAG(*)
      REAL SMAT(*), WORK(LWORK), XSMAT, EPSIBC, T2TCME(*)
      REAL TROCC2(*), TRVIR4(*), TRVIR5(*)
      REAL ZERO, ONE, TWO, FOUR, ECURR
      REAL DDOT
#else
      DOUBLE PRECISION T1AM(*), TMAT(*), XIAJB(*), FOCK(*)
      DOUBLE PRECISION TRVIR(*), TRVIR2(*), TROCC(*), FOCKD(*), DIAG(*)
      DOUBLE PRECISION SMAT(*), WORK(LWORK), XSMAT, EPSIBC, T2TCME(*)
      DOUBLE PRECISION TROCC2(*), TRVIR4(*), TRVIR5(*)
      DOUBLE PRECISION ZERO, ONE, TWO, FOUR, ECURR
      DOUBLE PRECISION DDOT
#endif
C
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, FOUR = 4.0D0)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CCFOP_SMAT')
C
      ISYRES  = MULD2H(ISYMT1,ISINT1)
      ISYRES2 = MULD2H(ISYMT2,ISINT2)
C
      IF (ISYRES .NE. ISYRES2)
     *   CALL QUIT('Symmetry mismatch in CCFOP_SMAT')
C
      ISYMBC = MULD2H(ISYMB,ISYMC)
      JSAIKJ = MULD2H(ISYRES,ISYMBC)
      LENGTH = NCKIJ(JSAIKJ)
C
C--------------------------------------------
C     First contribution from both T1 terms
C--------------------------------------------
C
      ISYMJK = MULD2H(ISYMBC,ISINT1)
C
      if (.true.) then
C----------------------------------------------
C     Sort integrals for constant B and C
C----------------------------------------------
C
      IF (LWORK .LT. NMATIJ(ISYMJK)) THEN
         CALL QUIT('Too little workspace in CCFOP_SMAT')
      ENDIF
C
      DO ISYMJ = 1, NSYM
C
         ISYMK  = MULD2H(ISYMJK,ISYMJ)
         ISYMCK = MULD2H(ISYMC,ISYMK)
         ISYMBJ = MULD2H(ISYMB,ISYMJ)
C
         DO J = 1, NRHF(ISYMJ)
C
            NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J-1) + B
C
            DO K = 1, NRHF(ISYMK)
C
               NKJ = IMATIJ(ISYMK,ISYMJ)+ NRHF(ISYMK)*(J - 1) + K
               NCK = IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K - 1) + C
C
               NCKBJ = IT2AM(ISYMCK,ISYMBJ) + INDEX(NCK,NBJ)
C
               WORK(NKJ) = XIAJB(NCKBJ)
C
            ENDDO ! K
         ENDDO    ! J
      ENDDO       ! ISYMJ
C
C---------------------------------------
C     Contract the integrals with T1.
C---------------------------------------
C
      CALL DZERO(TMAT,LENGTH)
C
      ISYMAI = ISYMT1
      DO ISYMJ = 1, NSYM
         ISYMK  = MULD2H(ISYMJK,ISYMJ)
         ISYAIK = MULD2H(ISYMK,ISYMAI)
C
         DO J = 1, NRHF(ISYMJ)
            DO K = 1, NRHF(ISYMK)
C
               NKJ = IMATIJ(ISYMK,ISYMJ)+ NRHF(ISYMK)*(J - 1) + K
C
               DO NAI = 1, NT1AM(ISYMAI)
C
                  NAIKJ = ISAIKJ(ISYAIK,ISYMJ)
     *                  + NCKI(ISYAIK)*(J - 1)
     *                  + ICKI(ISYMAI,ISYMK)
     *                  + NT1AM(ISYMAI)*(K-1) + NAI
C
                  TMAT(NAIKJ) = TWO*T1AM(NAI)*WORK(NKJ)
C
               ENDDO
            ENDDO
         ENDDO
C
      ENDDO
C
C----------------------------------------
C     Sum the result into SMAT.
C----------------------------------------
C
      JSAIKJ = MULD2H(ISYMAI,ISYMJK)
      DO I = 1, NCKIJ(JSAIKJ)
C         First :
          SMAT(I) = SMAT(I) + TMAT(I)
C         Second :
          SMAT(I) = SMAT(I) - TMAT(INDSQ(I,1))
      ENDDO
C
C-----------------------------------------------
C     Second contribution from both T1 terms
C-----------------------------------------------
C
      ISYAIK = MULD2H(ISINT1,ISYMC)
C
C------------------------------------
C     Sort integrals for constant C
C------------------------------------
C
      IF (LWORK .LT. NCKI(ISYAIK)) THEN
         CALL QUIT('Too little workspace in CCFOP_SMAT (2)')
      ENDIF
C
      DO ISYMK = 1, NSYM
         ISYMAI = MULD2H(ISYAIK,ISYMK)
         ISYMCK = MULD2H(ISYMC,ISYMK)
         DO K = 1, NRHF(ISYMK)
            NCK = IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K - 1) + C
            DO NAI = 1, NT1AM(ISYMAI)
C
               NAIK  = ICKI(ISYMAI,ISYMK)+NT1AM(ISYMAI)*(K - 1)+NAI
               NAICK = IT2AM(ISYMAI,ISYMCK) + INDEX(NAI,NCK)
C
               WORK(NAIK) = XIAJB(NAICK)
C
            ENDDO
         ENDDO
      ENDDO
C
C----------------------------------
C     Contract integrals with T1.
C----------------------------------
C
      CALL DZERO(TMAT,LENGTH)
C
      ISYMJ = MULD2H(ISYMT1,ISYMB)
C
      DO ISYMK = 1, NSYM
         ISYMAI = MULD2H(ISYAIK,ISYMK)
         ISYMCK = MULD2H(ISYMC,ISYMK)
C
         DO K = 1, NRHF(ISYMK)
            NCK = IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K - 1) + C
            DO J = 1, NRHF(ISYMJ)
C
               NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
C
               DO NAI = 1, NT1AM(ISYMAI)
C
                  NAIKJ = ISAIKJ(ISYAIK,ISYMJ)
     *                  + NCKI(ISYAIK)*(J - 1)
     *                  + ICKI(ISYMAI,ISYMK) 
     *                  + NT1AM(ISYMAI)*(K-1) + NAI
C
                  NAIK  = ICKI(ISYMAI,ISYMK)+ NT1AM(ISYMAI)*(K - 1)+ NAI
C
                  TMAT(NAIKJ) = TWO*T1AM(NBJ)*WORK(NAIK)
C
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C--------------------------------------
C     Sum the result into SMAT.
C--------------------------------------
C
      JSAIKJ = MULD2H(ISYAIK,ISYMJ)
      DO I = 1, NCKIJ(JSAIKJ)
C         First :
          SMAT(I) = SMAT(I) + TMAT(I)
C         Second :
          SMAT(I) = SMAT(I) - TMAT(INDSQ(I,3))
      ENDDO
C
      endif ! The end of the if (if .false.) then
C
C-----------------------------------------------------------------------
C     Contribution from both Fock terms
C-----------------------------------------------------------------------
C
      if (.true.) then
C
      CALL DZERO(TMAT,LENGTH)
C
      ISYMK  = MULD2H(ISINT2,ISYMC)
      ISYAIJ = MULD2H(ISYMT2,ISYMB)
      NCKTEMP = IT1AM(ISYMC,ISYMK) + C
C
      DO ISYMJ = 1, NSYM
         ISYMAI = MULD2H(ISYAIJ,ISYMJ)
         ISYMBJ = MULD2H(ISYMB,ISYMJ)
         ISYAIK = MULD2H(ISYMK,ISYMAI)
         DO ISYMI = 1, NSYM
            ISYMA = MULD2H(ISYMAI,ISYMI)
            ISYMAJ = MULD2H(ISYMA,ISYMJ)
            ISYMBI = MULD2H(ISYMB,ISYMI)
C
            DO J = 1, NRHF(ISYMJ)
               NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J-1) + B
C
               DO I = 1, NRHF(ISYMI)
               DO A = 1, NVIR(ISYMA)
C
                  NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
C
C                 Index for sorted T2 amplitudes.
C
                  NAIBJ = IT2SP(ISYAIJ,ISYMB)
     *                  + NCKI(ISYAIJ)*(B - 1)
     *                  + ICKI(ISYMAI,ISYMJ)
     *                  + NT1AM(ISYMAI)*(J - 1) + NAI
C
                  NAIKJTEMP = ISAIKJ(ISYAIK,ISYMJ)
     *                      + NCKI(ISYAIK)*(J - 1)
     *                      + ICKI(ISYMAI,ISYMK)
     *                      + NAI
C
                  DO K = 1, NRHF(ISYMK)
C
                     NCK = NCKTEMP + NVIR(ISYMC)*(K-1)
                     NAIKJ = NAIKJTEMP
     *                     + NT1AM(ISYMAI)*(K-1)
C
                     TMAT(NAIKJ) = TWO*T2TCME(NAIBJ)*FOCK(NCK)
C
                  ENDDO
               ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C------------------------------------
C     Sum the result into SMAT.
C------------------------------------
C
      ISYMBC = MULD2H(ISYMB,ISYMC)
      ISYRES = MULD2H(ISINT2,ISYMT2)
      JSAIKJ = MULD2H(ISYMBC,ISYRES)
C
      DO I = 1, NCKIJ(JSAIKJ)
         ! First term
         SMAT(I) = SMAT(I) + TWO*TMAT(I)
         ! Second term
         SMAT(I) = SMAT(I) - TMAT(INDSQ(I,3))
         SMAT(I) = SMAT(I) - TMAT(INDSQ(I,1))
      ENDDO
C
      endif
C
C----------------------------------------------
C     First virtual contribution of L term.
C----------------------------------------------
C
      ISYMBC = MULD2H(ISYMB,ISYMC)
      ISYRES = MULD2H(ISINT2,ISYMT2)
      JSAIKJ = MULD2H(ISYMBC,ISYRES)
      ISYMDK = MULD2H(ISYMBC,ISINT2)
C
      LENGTH = NCKIJ(JSAIKJ)
C
      IF (LWORK .LT. LENGTH) THEN
         CALL QUIT('Insufficient core in CCFOP_SMAT')
      ENDIF
C
      if (.true.) then
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
     *              NVIR(ISYMD),ONE,T2TCME(KOFF1),NTOAIJ,
     *              TRVIR5(KOFF2),NVIRD,ZERO,
     *              WORK(KOFF3),NTOAIJ)
C
      ENDDO
C
C      CALL CC_GATHER(LENGTH,SMAT,WORK,INDSQ(1,3))
      DO I = 1,LENGTH
         SMAT(I) = SMAT(I) + TWO*WORK(INDSQ(I,3))
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XSMAT = DDOT(NCKIJ(JSAIKJ),SMAT,1,SMAT,1)
         WRITE(LUPRI,*) 'In CCFOP_SMAT: 1. Norm of SMAT ',XSMAT
      ENDIF
C
C-------------------------------------------------
C     Second virtual contribution of L term.
C-------------------------------------------------
C
      ISYAKD = MULD2H(ISYMC,ISINT2)
      ISYDIJ = MULD2H(ISYMB,ISYMT2)
C
      DO ISYMJ = 1,NSYM
C
         ISYMDI = MULD2H(ISYMJ,ISYDIJ)
C
         DO J = 1,NRHF(ISYMJ)
C
            DO ISYMI = 1,NSYM
C
               ISYMD  = MULD2H(ISYMDI,ISYMI)
               ISYMAK = MULD2H(ISYMD,ISYAKD)
               ISYAKI = MULD2H(ISYMAK,ISYMI)
C
               KOFF1 = ICKATR(ISYMAK,ISYMD) + 1
C
               KOFF2 = IT2SP(ISYDIJ,ISYMB)
     *               + NCKI(ISYDIJ)*(B - 1)
     *               + ISAIK(ISYMDI,ISYMJ)
     *               + NT1AM(ISYMDI)*(J - 1)
     *               + IT1AM(ISYMD,ISYMI) + 1
C
               KOFF3 = ISAIKJ(ISYAKI,ISYMJ)
     *               + NCKI(ISYAKI)*(J - 1)
     *               + ISAIK(ISYMAK,ISYMI) + 1
C
               NVIRD  = MAX(NVIR(ISYMD),1)
               NTOTAK = MAX(NT1AM(ISYMAK),1)
C
               CALL DGEMM('N','N',NT1AM(ISYMAK),NRHF(ISYMI),
     *                    NVIR(ISYMD),ONE,TRVIR4(KOFF1),NTOTAK,
     *                    T2TCME(KOFF2),NVIRD,ZERO,TMAT(KOFF3),
     *                    NTOTAK)
C
            ENDDO
         ENDDO
      ENDDO
C
c     CALL CC_GATHER(LENGTH,WORK,TMAT,INDSQ(1,1))
c     CALL DAXPY(LENGTH,ONE,WORK,1,SMAT,1)
C
      DO I = 1,LENGTH
         SMAT(I) = SMAT(I) + TWO*TMAT(INDSQ(I,1))
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XSMAT = DDOT(NCKIJ(JSAIKJ),SMAT,1,SMAT,1)
         WRITE(LUPRI,*) 'In CCFOP_SMAT: 2. Norm of SMAT ',XSMAT
      ENDIF
C
      endif
C
C------------------------------------------------
C     First virtual contribution of g term.
C------------------------------------------------
C
      ISYMBC = MULD2H(ISYMB,ISYMC)
      ISYRES = MULD2H(ISINT2,ISYMT2)
      JSAIKJ = MULD2H(ISYMBC,ISYRES)
      ISYMDK = MULD2H(ISYMBC,ISINT2)
C
      LENGTH = NCKIJ(JSAIKJ)
C
      IF (LWORK .LT. LENGTH) THEN
         CALL QUIT('Insufficient core in CCSDT_SMAT')
      ENDIF
C
      if (.true.) then
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
     *              NVIR(ISYMD),ONE,T2TCME(KOFF1),NTOAIJ,
     *              TRVIR(KOFF2),NVIRD,ZERO,
     *              WORK(KOFF3),NTOAIJ)
C
      ENDDO
C
      DO I = 1,LENGTH
         SMAT(I) = SMAT(I) - TWO*WORK(INDSQ(I,2))
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XSMAT = DDOT(NCKIJ(JSAIKJ),SMAT,1,SMAT,1)
         WRITE(LUPRI,*) 'In CCFOP_SMAT: 3. Norm of SMAT ',XSMAT
      ENDIF
C
C-------------------------------------------------
C     Second virtual contribution of g term.
C-------------------------------------------------
C
      ISYAKD = MULD2H(ISYMC,ISINT2)
      ISYDIJ = MULD2H(ISYMB,ISYMT2)
C
      DO ISYMJ = 1,NSYM
C
         ISYMDI = MULD2H(ISYMJ,ISYDIJ)
C
         DO J = 1,NRHF(ISYMJ)
C
            DO ISYMI = 1,NSYM
C
               ISYMD  = MULD2H(ISYMDI,ISYMI)
               ISYMAK = MULD2H(ISYMD,ISYAKD)
               ISYAKI = MULD2H(ISYMAK,ISYMI)
C
               KOFF1 = ICKATR(ISYMAK,ISYMD) + 1
C
               KOFF2 = IT2SP(ISYDIJ,ISYMB)
     *               + NCKI(ISYDIJ)*(B - 1)
     *               + ISAIK(ISYMDI,ISYMJ)
     *               + NT1AM(ISYMDI)*(J - 1)
     *               + IT1AM(ISYMD,ISYMI) + 1
C
               KOFF3 = ISAIKJ(ISYAKI,ISYMJ)
     *               + NCKI(ISYAKI)*(J - 1)
     *               + ISAIK(ISYMAK,ISYMI) + 1
C
               NVIRD  = MAX(NVIR(ISYMD),1)
               NTOTAK = MAX(NT1AM(ISYMAK),1)
C
               CALL DGEMM('N','N',NT1AM(ISYMAK),NRHF(ISYMI),
     *                    NVIR(ISYMD),ONE,TRVIR2(KOFF1),NTOTAK,
     *                    T2TCME(KOFF2),NVIRD,ZERO,TMAT(KOFF3),
     *                    NTOTAK)
C
            ENDDO
         ENDDO
      ENDDO
C
      DO I = 1,LENGTH
         SMAT(I) = SMAT(I) - TWO*TMAT(INDSQ(I,4))
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XSMAT = DDOT(NCKIJ(JSAIKJ),SMAT,1,SMAT,1)
         WRITE(LUPRI,*) 'In CCFOP_SMAT: 4. Norm of SMAT ',XSMAT
      ENDIF
C
      endif
C
C----------------------------------------
C     First occupied L contribution.
C----------------------------------------
C
      ISYAIL = MULD2H(ISYMB,ISYMT2)
      ISYLKJ = MULD2H(ISYMC,ISINT2)
C
      if (.true.) then
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
     *                    NRHF(ISYML),-TWO,T2TCME(KOFF1),NTOTAI,
     *                    TROCC2(KOFF2),NRHFL,ONE,SMAT(KOFF3),
     *                    NTOTAI)
C
            ENDDO
         ENDDO
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XSMAT = DDOT(NCKIJ(JSAIKJ),SMAT,1,SMAT,1)
         WRITE(LUPRI,*) 'In CC3_SMAT: 3. Norm of SMAT ',XSMAT
      ENDIF
C
C----------------------------------------
C     Second occupied L contribution.
C----------------------------------------
C
      ISYAJL = MULD2H(ISYMB,ISYMT2)
      ISYLKI = MULD2H(ISYMC,ISINT2)
C
      IF (LWORK .LT. NCKI(ISYAJL)) THEN
         CALL QUIT('Not enough space in CCSDT_SMAT')
      END IF
C
      KOFF = IT2SP(ISYAJL,ISYMB) + NCKI(ISYAJL)*(B - 1) + 1
      CALL CC_GATHER(NCKI(ISYAJL),WORK,T2TCME(KOFF),INDAJL)
C
      DO ISYMI = 1,NSYM
C
         ISYMLK = MULD2H(ISYMI,ISYLKI)
C
         DO I = 1,NRHF(ISYMI)
C
            DO ISYMK = 1,NSYM
C
               ISYML  = MULD2H(ISYMK,ISYMLK)
               ISYMAJ = MULD2H(ISYAJL,ISYML)
               ISYAJK = MULD2H(ISYMAJ,ISYMK)
C
               KOFF1 = ICKI(ISYMAJ,ISYML) + 1
C
               KOFF2 = ISJIKA(ISYLKI,ISYMC)
     *               + NMAJIK(ISYLKI)*(C - 1)
     *               + ISJIK(ISYMLK,ISYMI)
     *               + NMATIJ(ISYMLK)*(I - 1)
     *               + IMATIJ(ISYML,ISYMK) + 1
C
               KOFF3 = ISAIKJ(ISYAJK,ISYMI)
     *               + NCKI(ISYAJK)*(I - 1)
     *               + ICKI(ISYMAJ,ISYMK) + 1
C
               NTOTAJ = MAX(1,NT1AM(ISYMAJ))
               NRHFL  = MAX(1,NRHF(ISYML))
C
               CALL DGEMM('N','N',NT1AM(ISYMAJ),NRHF(ISYMK),
     *                    NRHF(ISYML),TWO,WORK(KOFF1),NTOTAJ,
     *                    TROCC2(KOFF2),NRHFL,ZERO,TMAT(KOFF3),
     *                    NTOTAJ)
C
            ENDDO
         ENDDO
      ENDDO
C
c     CALL CC_GATHER(NCKIJ(JSAIKJ),WORK,TMAT,INDSQ(1,5))
c     CALL DAXPY(NCKIJ(JSAIKJ),-ONE,WORK,1,SMAT,1)
C
      DO I = 1,NCKIJ(JSAIKJ)
         SMAT(I) = SMAT(I) - TMAT(INDSQ(I,5))
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XSMAT = DDOT(NCKIJ(JSAIKJ),SMAT,1,SMAT,1)
         WRITE(LUPRI,*) 'In CC3_SMAT: 4. Norm of SMAT ',XSMAT
      ENDIF
C
      endif
C
C---------------------------------------
C     First occupied g contribution.
C---------------------------------------
C
      ISYAIL = MULD2H(ISYMB,ISYMT2)
      ISYLKJ = MULD2H(ISYMC,ISINT2)
C
      if (.true.) then
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
     *                    NRHF(ISYML),-TWO,T2TCME(KOFF1),NTOTAI,
     *                    TROCC(KOFF2),NRHFL,ZERO,TMAT(KOFF3),
     *                    NTOTAI)
C
            ENDDO
         ENDDO
      ENDDO
C
      DO I = 1,NCKIJ(JSAIKJ)
         SMAT(I) = SMAT(I) - TMAT(INDSQ(I,1))
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XSMAT = DDOT(NCKIJ(JSAIKJ),SMAT,1,SMAT,1)
         WRITE(LUPRI,*) 'In CC3_SMAT: 3. Norm of SMAT ',XSMAT
      ENDIF
C
C----------------------------------------
C     Second occupied g contribution.
C----------------------------------------
C
      ISYAJL = MULD2H(ISYMB,ISYMT2)
      ISYLKI = MULD2H(ISYMC,ISINT2)
C
      IF (LWORK .LT. NCKI(ISYAJL)) THEN
         CALL QUIT('Not enough space in CCSDT_SMAT')
      END IF
C
      KOFF = IT2SP(ISYAJL,ISYMB) + NCKI(ISYAJL)*(B - 1) + 1
      CALL CC_GATHER(NCKI(ISYAJL),WORK,T2TCME(KOFF),INDAJL)
C
      DO ISYMI = 1,NSYM
C
         ISYMLK = MULD2H(ISYMI,ISYLKI)
C
         DO I = 1,NRHF(ISYMI)
C
            DO ISYMK = 1,NSYM
C
               ISYML  = MULD2H(ISYMK,ISYMLK)
               ISYMAJ = MULD2H(ISYAJL,ISYML)
               ISYAJK = MULD2H(ISYMAJ,ISYMK)
C
               KOFF1 = ICKI(ISYMAJ,ISYML) + 1
C
               KOFF2 = ISJIKA(ISYLKI,ISYMC)
     *               + NMAJIK(ISYLKI)*(C - 1)
     *               + ISJIK(ISYMLK,ISYMI)
     *               + NMATIJ(ISYMLK)*(I - 1)
     *               + IMATIJ(ISYML,ISYMK) + 1
C
               KOFF3 = ISAIKJ(ISYAJK,ISYMI)
     *               + NCKI(ISYAJK)*(I - 1)
     *               + ICKI(ISYMAJ,ISYMK) + 1
C
               NTOTAJ = MAX(1,NT1AM(ISYMAJ))
               NRHFL  = MAX(1,NRHF(ISYML))
C
               CALL DGEMM('N','N',NT1AM(ISYMAJ),NRHF(ISYMK),
     *                    NRHF(ISYML),TWO,WORK(KOFF1),NTOTAJ,
     *                    TROCC(KOFF2),NRHFL,ZERO,TMAT(KOFF3),
     *                    NTOTAJ)
C
            ENDDO
         ENDDO
      ENDDO
C
c     CALL CC_GATHER(NCKIJ(JSAIKJ),WORK,TMAT,INDSQ(1,5))
c     CALL DAXPY(NCKIJ(JSAIKJ),-ONE,WORK,1,SMAT,1)
C
      DO I = 1,NCKIJ(JSAIKJ)
         SMAT(I) = SMAT(I) + TMAT(INDSQ(I,2))
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XSMAT = DDOT(NCKIJ(JSAIKJ),SMAT,1,SMAT,1)
         WRITE(LUPRI,*) 'In CC3_SMAT: 4. Norm of SMAT ',XSMAT
      ENDIF
C
      endif
C
C-----------------------------------------
C     Divide by the Fock matrix diagonals.
C-----------------------------------------
C
      NB = IORB(ISYMB) + NRHF(ISYMB) + B
      NC = IORB(ISYMC) + NRHF(ISYMC) + C
C
      EPSIBC = FOCKD(NB) + FOCKD(NC) - ECURR
C
      DO L = 1,NCKIJ(JSAIKJ)
         SMAT(L) = SMAT(L)/(DIAG(L) + EPSIBC)
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XSMAT = DDOT(NCKIJ(JSAIKJ),SMAT,1,SMAT,1)
         WRITE(LUPRI,*) 'In CC3_SMAT: 5. Norm of SMAT ',XSMAT
      ENDIF
C
C----------------------
C     Print if desired.
C----------------------
C
      IF (IPRCC .GT. 75) THEN
C
         CALL AROUND('The S matrix')
         WRITE(LUPRI,*)
         WRITE(LUPRI,'(2X,A,I4)')  'JSAIKJ ', JSAIKJ
         WRITE(LUPRI,'(2X,A,4I4)') 'isymb,b,isymc,c',ISYMB,B,ISYMC,C
         WRITE(LUPRI,*)
C
         DO ISYMJ = 1,NSYM
C
            ISYAIK = MULD2H(JSAIKJ,ISYMJ)
C
            DO J = 1,NRHF(ISYMJ)
C
               WRITE(LUPRI,'(5X,A,2I4)') 'isymj,j',ISYMJ,J
               WRITE(LUPRI,*)
C
               DO ISYMK = 1,NSYM
C
                  ISYMAI = MULD2H(ISYAIK,ISYMK)
C
                  DO K = 1,NRHF(ISYMK)
C
                     WRITE(LUPRI,'(8X,A,2I4)') 'isymk,k',ISYMK,K
                     WRITE(LUPRI,*)
C
                     DO ISYMI = 1,NSYM
C
                        ISYMA = MULD2H(ISYMAI,ISYMI)
C
                        KOFF = ISAIKJ(ISYAIK,ISYMJ)
     *                       + NCKI(ISYAIK)*(J - 1)
     *                       + ICKI(ISYMAI,ISYMK)
     *                       + NT1AM(ISYMAI)*(K - 1)
     *                       + IT1AM(ISYMA,ISYMI) + 1
C
                        CALL OUTPUT(SMAT(KOFF),1,NVIR(ISYMA),1,
     *                              NRHF(ISYMI),NVIR(ISYMA),
     *                              NRHF(ISYMI),1,LUPRI)
C
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
C
      END IF
C
      CALL QEXIT('CCFOP_SMAT')
C
      RETURN
      END
C  /* Deck ccsdt_tcmeocc */
      SUBROUTINE CCSDT_TCMEOCC(TRINP,TROUT,ISYINT)
C
C     Kasper Hald, Fall 2001
C
C     Calculate 2*C-E from the occupied integrals
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INTEGER ISYINT
      INTEGER ISYMA, ISYJIK, ISYMK, ISYMJI, ISYMJ, ISYMI, KOFF1, KOFF2
      integer isymjk
C
#if defined (SYS_CRAY)
      REAL TRINP(*), TROUT(*), TWO
#else
      DOUBLE PRECISION TRINP(*), TROUT(*), TWO
#endif
C
      PARAMETER(TWO = 2.0D0)
C
      CALL QENTER('CCSDT_TCMEOCC')
C
      DO ISYMA = 1, NSYM
C
         ISYJIK = MULD2H(ISYINT,ISYMA)
C
         DO ISYMK = 1, NSYM
            ISYMJI = MULD2H(ISYJIK,ISYMK)
            DO ISYMJ = 1, NSYM
               ISYMI = MULD2H(ISYMJI,ISYMJ)
               ISYMJK = MULD2H(ISYMJ,ISYMK)
               DO K = 1, NRHF(ISYMK)
               DO A = 1, NVIR(ISYMA)
                  DO J = 1, NRHF(ISYMJ)
                  DO I = 1, NRHF(ISYMI)
C
                     KOFF1 = ISJIKA(ISYJIK,ISYMA)
     *                     + NMAJIK(ISYJIK)*(A-1)
     *                     + ISJIK(ISYMJI,ISYMK)
     *                     + NMATIJ(ISYMJI)*(K - 1)
     *                     + IMATIJ(ISYMJ,ISYMI)
     *                     + NRHF(ISYMJ)*(I - 1) + J
C
C     First :
C                     KOFF2 = ISJIKA(ISYJIK,ISYMA)
C     *                     + NMAJIK(ISYJIK)*(A-1)
C     *                     + ISJIK(ISYMJI,ISYMK)
C     *                     + NMATIJ(ISYMJI)*(K - 1)
C     *                     + IMATIJ(ISYMI,ISYMJ)
C     *                     + NRHF(ISYMI)*(J - 1) + I
C
                     KOFF2 = ISJIKA(ISYJIK,ISYMA)
     *                     + NMAJIK(ISYJIK)*(A-1)
     *                     + ISJIK(ISYMJK,ISYMI)
     *                     + NMATIJ(ISYMJK)*(I - 1)
     *                     + IMATIJ(ISYMJ,ISYMK)
     *                     + NRHF(ISYMJ)*(K - 1) + J
C
C
                     TROUT(KOFF1) = TWO*TRINP(KOFF1) - TRINP(KOFF2)
C
                  ENDDO   ! I
                  ENDDO   ! J
               ENDDO      ! A
               ENDDO      ! K
            ENDDO         ! ISYMA
         ENDDO            ! ISYMJ
      ENDDO               ! ISYMK
C
      CALL QEXIT('CCSDT_TCMEOCC')
C
      RETURN
      END
C  /* Deck ccfop_qmat */
      SUBROUTINE CCFOP_QMAT(ECURR,T1AM,ISYMT1,T2TCME,ISYMT2,TMAT,FOCK,
     *                      XIAJB,ISINT1,TRVIR3,TRVIR6,TROCC,TROCC2,
     *                      ISINT2,FOCKD,DIAG,QMAT,WORK,LWORK,
     *                      INDSQ,LENSQ,ISYMB,B,ISYMD,D)
C
C     Written by K. Hald, Fall 2001.
C
C     Calculate QMAT for t3-bar.
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INTEGER ISYMT1, ISYMT2, ISINT1, ISINT2, LWORK, LENSQ
      INTEGER ISYMB, ISYMD
      INTEGER INDSQ(LENSQ,6), INDEX
      INTEGER ISYMBD, ISYMJK, ISYMJ, ISYMK, ISYMBK, ISYMDJ, NDJ, NKJ
      INTEGER NBK, NBKDJ, ISYMAI, ISYAIK, ISYMA, ISYMI, NAI, NAIKJ
      INTEGER JSAIKJ, ISYAIJ, NBKTEMP, ISYMAJ, ISYMDI
      INTEGER NAIDJ, NAIKJTEMP, ISYRES, ISYMFK, ISYMF, NVIRF
      INTEGER LENGTH, KOFF1, KOFF2, KOFF3, NTOAIJ, ISYAIL, ISYLKJ
      INTEGER ISYMLK, ISYML, NTOTAI, NRHFL, NB, ND, NAIK, NAIBK
C
#if defined (SYS_CRAY)
      REAL T1AM(*), T2TCME(*), FOCK(*), XIAJB(*), TRVIR3(*)
      REAL TRVIR6(*), TROCC(*), TROCC2(*), FOCKD(*), DIAG(*)
      REAL TMAT(*), QMAT(*), WORK(LWORK), EPSIBD, ECURR
      REAL ZERO, ONE, TWO, FOUR, XQMAT, DDOT
#else
      DOUBLE PRECISION T1AM(*), T2TCME(*), FOCK(*), XIAJB(*), TRVIR3(*)
      DOUBLE PRECISION TRVIR6(*), TROCC(*), TROCC2(*), FOCKD(*),DIAG(*)
      DOUBLE PRECISION TMAT(*), QMAT(*), WORK(LWORK), EPSIBD, ECURR
      DOUBLE PRECISION ZERO, ONE, TWO, FOUR, XQMAT, DDOT
#endif
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, FOUR = 4.0D0)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CCFOP_QMAT')
C
C------------------------------------------
C     Contribution from the two T1 terms.
C------------------------------------------
      if (.true.) then
C
      ISYAIK = MULD2H(ISINT1,ISYMB)
      ISYMJ = MULD2H(ISYMT1,ISYMD)
      JSAIKJ = MULD2H(ISYAIK,ISYMJ)
C
      CALL DZERO(TMAT,NCKIJ(JSAIKJ))
C
C------------------------------------
C     Sort integrals for constant B
C------------------------------------
C
      IF (LWORK .LT. NCKI(ISYAIK)) THEN
         CALL QUIT('Too little workspace in CCFOP_QMAT (1)')
      ENDIF
C
      DO ISYMK = 1, NSYM
         ISYMAI = MULD2H(ISYAIK,ISYMK)
         ISYMBK = MULD2H(ISYMB,ISYMK)
         DO K = 1, NRHF(ISYMK)
            NBK = IT1AM(ISYMB,ISYMK) + NVIR(ISYMB)*(K - 1) + B
            DO NAI = 1, NT1AM(ISYMAI)
C
               NAIK  = ICKI(ISYMAI,ISYMK)+NT1AM(ISYMAI)*(K - 1)+NAI
               NAIBK = IT2AM(ISYMAI,ISYMBK) + INDEX(NAI,NBK)
C
               WORK(NAIK) = XIAJB(NAIBK)
C
            ENDDO
         ENDDO
      ENDDO
C
C----------------------------------
C     Contract integrals with T1.
C----------------------------------
C
      DO ISYMK = 1, NSYM
         ISYMAI = MULD2H(ISYAIK,ISYMK)
C
         DO K = 1, NRHF(ISYMK)
            DO J = 1, NRHF(ISYMJ)
C
               NDJ = IT1AM(ISYMD,ISYMJ) + NVIR(ISYMD)*(J - 1) + D
C
               DO NAI = 1, NT1AM(ISYMAI)
C
                  NAIKJ = ISAIKJ(ISYAIK,ISYMJ)
     *                  + NCKI(ISYAIK)*(J - 1)
     *                  + ICKI(ISYMAI,ISYMK) 
     *                  + NT1AM(ISYMAI)*(K-1) + NAI
C
                  NAIK  = ICKI(ISYMAI,ISYMK)+ NT1AM(ISYMAI)*(K - 1)+ NAI
C
                  TMAT(NAIKJ) = TWO*T1AM(NDJ)*WORK(NAIK)
C
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C--------------------------------------
C     Sum the result into SMAT.
C--------------------------------------
C
      DO I = 1, NCKIJ(JSAIKJ)
C         First :
          QMAT(I) = QMAT(I) + TMAT(I)
C         Second :
          QMAT(I) = QMAT(I) - TMAT(INDSQ(I,3))
      ENDDO
C
      endif
C
C-----------------------------------------------------------------------
C     Contribution from both Fock terms
C-----------------------------------------------------------------------
C
      if (.true.) then
C
      ISYMK  = MULD2H(ISINT2,ISYMB)
      ISYAIJ = MULD2H(ISYMT2,ISYMD)
      NBKTEMP = IT1AM(ISYMB,ISYMK) + B
      JSAIKJ  = MULD2H(ISYMK,ISYAIJ)
C
      CALL DZERO(TMAT,NCKIJ(JSAIKJ))
C
      DO ISYMJ = 1, NSYM
         ISYMAI = MULD2H(ISYAIJ,ISYMJ)
         ISYMDJ = MULD2H(ISYMD,ISYMJ)
         ISYAIK = MULD2H(ISYMK,ISYMAI)
         DO ISYMI = 1, NSYM
            ISYMA = MULD2H(ISYMAI,ISYMI)
            ISYMAJ = MULD2H(ISYMA,ISYMJ)
            ISYMDI = MULD2H(ISYMD,ISYMI)
C
            DO J = 1, NRHF(ISYMJ)
               NDJ = IT1AM(ISYMD,ISYMJ) + NVIR(ISYMD)*(J-1) + D
C
               DO I = 1, NRHF(ISYMI)
               DO A = 1, NVIR(ISYMA)
C
                  NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
C
C                 Index for sorted T2 amplitudes.
C
                  NAIDJ = IT2SP(ISYAIJ,ISYMD)
     *                  + NCKI(ISYAIJ)*(D - 1)
     *                  + ICKI(ISYMAI,ISYMJ)
     *                  + NT1AM(ISYMAI)*(J - 1) + NAI
C
                  NAIKJTEMP = ISAIKJ(ISYAIK,ISYMJ)
     *                      + NCKI(ISYAIK)*(J - 1)
     *                      + ICKI(ISYMAI,ISYMK)
     *                      + NAI
C
                  DO K = 1, NRHF(ISYMK)
C
                     NBK = NBKTEMP + NVIR(ISYMB)*(K-1)
                     NAIKJ = NAIKJTEMP
     *                     + NT1AM(ISYMAI)*(K-1)
C
                     TMAT(NAIKJ) = TWO*T2TCME(NAIDJ)*FOCK(NBK)
C
                  ENDDO
               ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C------------------------------------
C     Sum the result into SMAT.
C------------------------------------
C
      DO I = 1, NCKIJ(JSAIKJ)
         ! First term
         QMAT(I) = QMAT(I) + TMAT(I)
         ! Second term
         QMAT(I) = QMAT(I) - TMAT(INDSQ(I,3))
      ENDDO
C
      endif
C----------------------------------------------
C     Virtual contribution of L term.
C----------------------------------------------
C
      ISYMBD = MULD2H(ISYMB,ISYMD)
      ISYRES = MULD2H(ISINT2,ISYMT2)
      JSAIKJ = MULD2H(ISYMBD,ISYRES)
      ISYMFK = MULD2H(ISYMBD,ISINT2)
C
      LENGTH = NCKIJ(JSAIKJ)
C
      IF (LWORK .LT. LENGTH) THEN
         CALL QUIT('Insufficient core in CCSDT_SMAT')
      ENDIF
C
      if (.true.) then
C
      DO ISYMK = 1,NSYM
C
         ISYMF  = MULD2H(ISYMK,ISYMFK)
         ISYAIJ = MULD2H(ISYMK,JSAIKJ)
C
         KOFF1 = IT2SP(ISYAIJ,ISYMF)  + 1
         KOFF2 = ICKATR(ISYMFK,ISYMB) + NT1AM(ISYMFK)*(B - 1)
     *         + IT1AM(ISYMF,ISYMK)   + 1
         KOFF3 = ISAIKJ(ISYAIJ,ISYMK) + 1
C
         NTOAIJ = MAX(1,NCKI(ISYAIJ))
         NVIRF  = MAX(NVIR(ISYMF),1)
C
         CALL DGEMM('N','N',NCKI(ISYAIJ),NRHF(ISYMK),
     *              NVIR(ISYMF),ONE,T2TCME(KOFF1),NTOAIJ,
     *              TRVIR6(KOFF2),NVIRF,ZERO,
     *              WORK(KOFF3),NTOAIJ)
C
      ENDDO
C
C      CALL CC_GATHER(LENGTH,SMAT,WORK,INDSQ(1,3))
      DO I = 1,LENGTH
         QMAT(I) = QMAT(I) + TWO*WORK(INDSQ(I,3))
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XQMAT = DDOT(NCKIJ(JSAIKJ),QMAT,1,QMAT,1)
         WRITE(LUPRI,*) 'In CCFOP_QMAT: 1. Norm of QMAT ',XQMAT
      ENDIF
C
      endif
C
C----------------------------------------------
C     Virtual contribution of g term.
C----------------------------------------------
C
      ISYMBD = MULD2H(ISYMB,ISYMD)
      ISYRES = MULD2H(ISINT2,ISYMT2)
      JSAIKJ = MULD2H(ISYMBD,ISYRES)
      ISYMFK = MULD2H(ISYMBD,ISINT2)
C
      LENGTH = NCKIJ(JSAIKJ)
C
      IF (LWORK .LT. LENGTH) THEN
         CALL QUIT('Insufficient core in CCSDT_SMAT')
      ENDIF
C
      if (.true.) then
C
      DO ISYMK = 1,NSYM
C
         ISYMF  = MULD2H(ISYMK,ISYMFK)
         ISYAIJ = MULD2H(ISYMK,JSAIKJ)
C
         KOFF1 = IT2SP(ISYAIJ,ISYMF)  + 1
         KOFF2 = ICKATR(ISYMFK,ISYMB) + NT1AM(ISYMFK)*(B - 1)
     *         + IT1AM(ISYMF,ISYMK)   + 1
         KOFF3 = ISAIKJ(ISYAIJ,ISYMK) + 1
C
         NTOAIJ = MAX(1,NCKI(ISYAIJ))
         NVIRF  = MAX(NVIR(ISYMF),1)
C
         CALL DGEMM('N','N',NCKI(ISYAIJ),NRHF(ISYMK),
     *              NVIR(ISYMF),ONE,T2TCME(KOFF1),NTOAIJ,
     *              TRVIR3(KOFF2),NVIRF,ZERO,
     *              WORK(KOFF3),NTOAIJ)
C
      ENDDO
C
      DO I = 1,LENGTH
         QMAT(I) = QMAT(I) - TWO*WORK(INDSQ(I,2))
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XQMAT = DDOT(NCKIJ(JSAIKJ),QMAT,1,QMAT,1)
         WRITE(LUPRI,*) 'In CCFOP_QMAT: 2. Norm of QMAT ',XQMAT
      ENDIF
C
C
      endif
C
C----------------------------------------
C     Occupied L contribution.
C----------------------------------------
C
      ISYAIL = MULD2H(ISYMD,ISYMT2)
      ISYLKJ = MULD2H(ISYMB,ISINT2)
C
      if (.true.) then
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
               KOFF1 = IT2SP(ISYAIL,ISYMD)
     *               + NCKI(ISYAIL)*(D - 1)
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
     *                    NRHF(ISYML),-TWO,T2TCME(KOFF1),NTOTAI,
     *                    TROCC2(KOFF2),NRHFL,ONE,QMAT(KOFF3),
     *                    NTOTAI)
C
            ENDDO
         ENDDO
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XQMAT = DDOT(NCKIJ(JSAIKJ),QMAT,1,QMAT,1)
         WRITE(LUPRI,*) 'In CCFOP_QMAT: 3. Norm of QMAT ',XQMAT
      ENDIF
C
      endif
C---------------------------------------
C     Occupied g contribution.
C---------------------------------------
C
      ISYAIL = MULD2H(ISYMD,ISYMT2)
      ISYLKJ = MULD2H(ISYMB,ISINT2)
C
      if (.true.) then
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
               KOFF1 = IT2SP(ISYAIL,ISYMD)
     *               + NCKI(ISYAIL)*(D - 1)
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
     *                    NRHF(ISYML),-TWO,T2TCME(KOFF1),NTOTAI,
     *                    TROCC(KOFF2),NRHFL,ZERO,TMAT(KOFF3),
     *                    NTOTAI)
C
            ENDDO
         ENDDO
      ENDDO
C
      DO I = 1,NCKIJ(JSAIKJ)
         QMAT(I) = QMAT(I) - TMAT(INDSQ(I,1))
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XQMAT = DDOT(NCKIJ(JSAIKJ),QMAT,1,QMAT,1)
         WRITE(LUPRI,*) 'In CCFOP_QMAT: 4. Norm of QMAT ',XQMAT
      ENDIF
C
      endif
C
C-----------------------------------------
C     Divide by the Fock matrix diagonals.
C-----------------------------------------
C
      NB = IORB(ISYMB) + NRHF(ISYMB) + B
      ND = IORB(ISYMD) + NRHF(ISYMD) + D
C
      EPSIBD = FOCKD(NB) + FOCKD(ND) - ECURR
C
      DO L = 1,NCKIJ(JSAIKJ)
         QMAT(L) = QMAT(L)/(DIAG(L) + EPSIBD)
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XQMAT = DDOT(NCKIJ(JSAIKJ),QMAT,1,QMAT,1)
         WRITE(LUPRI,*) 'In CCFOP_SMAT: 5. Norm of QMAT ',XQMAT
      ENDIF
C
C----------------------
C     Print if desired.
C----------------------
CC
C      IF (IPRCC .GT. 75) THEN
CC
C         CALL AROUND('The S matrix')
C         WRITE(LUPRI,*)
C         WRITE(LUPRI,'(2X,A,I4)')  'JSAIKJ ', JSAIKJ
C         WRITE(LUPRI,'(2X,A,4I4)') 'isymb,b,isymd,d',ISYMB,B,ISYMD,D
C         WRITE(LUPRI,*)
CC
C         DO ISYMJ = 1,NSYM
CC
C            ISYAIK = MULD2H(JSAIKJ,ISYMJ)
CC
C            DO J = 1,NRHF(ISYMJ)
CC
C               WRITE(LUPRI,'(5X,A,2I4)') 'isymj,j',ISYMJ,J
C               WRITE(LUPRI,*)
CC
C               DO ISYMK = 1,NSYM
CC
C                  ISYMAI = MULD2H(ISYAIK,ISYMK)
CC
C                  DO K = 1,NRHF(ISYMK)
CC
C                     WRITE(LUPRI,'(8X,A,2I4)') 'isymk,k',ISYMK,K
C                     WRITE(LUPRI,*)
CC
C                     DO ISYMI = 1,NSYM
CC
C                        ISYMA = MULD2H(ISYMAI,ISYMI)
CC
C                        KOFF1 = ISAIKJ(ISYAIK,ISYMJ)
C     *                        + NCKI(ISYAIK)*(J - 1)
C     *                        + ICKI(ISYMAI,ISYMK)
C     *                        + NT1AM(ISYMAI)*(K - 1)
C     *                        + IT1AM(ISYMA,ISYMI) + 1
CC
C                        CALL OUTPUT(QMAT(KOFF1),1,NVIR(ISYMA),1,
C     *                              NRHF(ISYMI),NVIR(ISYMA),
C     *                              NRHF(ISYMI),1,LUPRI)
CC
C                     ENDDO
C                  ENDDO
C               ENDDO
C            ENDDO
C         ENDDO
CC
C      END IF
CC
C---------------------------------
C     Finish
C---------------------------------
C
      CALL QEXIT('CCFOP_QMAT')
C
      RETURN
      END
C  /* Deck ccfop_convir */
      SUBROUTINE CCFOP_CONVIR(RMAT,SMAT,QMAT,TMAT,ISYMIM,TRVIR,
     *                        TRVIR1,ISYINT,WORK,LWORK,INDSQ,LENSQ,
     *                        ISYMB,B,ISYMD,D)
C
C     K. Hald, Fall 2001, Nearly identical with cc3_convir.
C     Need only one of the amplidudes.
C
C     General symmetry: ISYMIM is the symmetry of the SMAT, QMAT 
C                       and TMAT intermdiates.
C                       ISYINT is symmetry of FOCKAK and XIAJB
C                       ISYRES = ISYMIM*ISYINT
C
#include <implicit.h>
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION RMAT(*),SMAT(*),QMAT(*)
      DIMENSION TMAT(*),TRVIR(*),TRVIR1(*)
      DIMENSION WORK(LWORK),INDSQ(LENSQ,6)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CCFOP_CONVIR')
C
      ISYRES = MULD2H(ISYMIM,ISYINT)
C
      ISYMBD = MULD2H(ISYMB,ISYMD)
      ISCKIJ = MULD2H(ISYMBD,ISYMIM)
C
      LENGTH = NCKIJ(ISCKIJ)
C
C------------------------
C     First virtual term.
C------------------------
C
      IF (LWORK .LT. NCKIJ(ISCKIJ)) THEN
         CALL QUIT('Insufficient core in CCSDT_CONVIR')
      ENDIF
C
C
      DO I = 1,LENGTH
C
         TMAT(I) = - SMAT(INDSQ(I,1))
     *           -   QMAT(INDSQ(I,2))
C
      ENDDO
C
C---------------------------
C     Contract with (ac|kd).
C---------------------------
C
      DO ISYMJ = 1,NSYM
C
         ISYMBJ = MULD2H(ISYMB,ISYMJ)
         ISYMAI = MULD2H(ISYMBJ,ISYRES)
         ISYCKI = MULD2H(ISCKIJ,ISYMJ)
C
         KSCR1  = 1
         KEND1  = KSCR1 + NT1AM(ISYMAI)
         LWRK1  = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CCSDT_CONVIR')
         ENDIF
C
         DO J = 1,NRHF(ISYMJ)
C
            DO ISYMI = 1,NSYM
C
               ISYMCK = MULD2H(ISYCKI,ISYMI)
               ISYMA  = MULD2H(ISYMAI,ISYMI)
C
               NTOTCK = MAX(NT1AM(ISYMCK),1)
               NVIRA  = MAX(NVIR(ISYMA),1)
C
               KOFF1 = ICKATR(ISYMCK,ISYMA) + 1
               KOFF2 = ISAIKJ(ISYCKI,ISYMJ)
     *               + NCKI(ISYCKI)*(J - 1)
     *               + ISAIK(ISYMCK,ISYMI)  + 1
               KOFF3 = ISAIK(ISYMAI,ISYMJ)
     *               + NT1AM(ISYMAI)*(J - 1)
     *               + IT1AM(ISYMA,ISYMI) + 1
C
               CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NT1AM(ISYMCK),
     *                    ONE,TRVIR(KOFF1),NTOTCK,TMAT(KOFF2),NTOTCK,
     *                    ONE,RMAT(KOFF3),NVIRA)
C
            ENDDO
         ENDDO
      ENDDO
C
C-------------------------
C     Second virtual term.
C-------------------------
C
C
      DO I = 1,LENGTH
C
         TMAT(I) = - SMAT(I) 
     *             - QMAT(INDSQ(I,3))
C
      ENDDO
C
C---------------------------
C     Contract with (ad|kc).
C---------------------------
C
      DO ISYMJ = 1,NSYM
C
         ISYMBJ = MULD2H(ISYMB,ISYMJ)
         ISYMAI = MULD2H(ISYMBJ,ISYRES)
         ISYCKI = MULD2H(ISCKIJ,ISYMJ)
C
         KSCR1  = 1
         KEND1  = KSCR1 + NT1AM(ISYMAI)
         LWRK1  = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CCSDT_CONVIR')
         ENDIF
C
         DO J = 1,NRHF(ISYMJ)
C
            DO ISYMI = 1,NSYM
C
               ISYMCK = MULD2H(ISYCKI,ISYMI)
               ISYMA  = MULD2H(ISYMAI,ISYMI)
C
               NTOTCK = MAX(NT1AM(ISYMCK),1)
               NVIRA  = MAX(NVIR(ISYMA),1)
C
               KOFF1 = ICKATR(ISYMCK,ISYMA) + 1
               KOFF2 = ISAIKJ(ISYCKI,ISYMJ)
     *               + NCKI(ISYCKI)*(J - 1)
     *               + ISAIK(ISYMCK,ISYMI)  + 1
               KOFF3 = ISAIK(ISYMAI,ISYMJ)
     *               + NT1AM(ISYMAI)*(J - 1)
     *               + IT1AM(ISYMA,ISYMI) + 1
C
               CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NT1AM(ISYMCK),
     *                    ONE,TRVIR1(KOFF1),NTOTCK,TMAT(KOFF2),NTOTCK,
     *                    ONE,RMAT(KOFF3),NVIRA)
C
            ENDDO
         ENDDO
      ENDDO
C
      CALL QEXIT('CCFOP_CONVIR')
C
      RETURN
      END
C  /* Deck ccfop_conocc */
      SUBROUTINE CCFOP_CONOCC(OMEGA2,RMAT1,RMAT2,SMAT,TMAT,ISYMIM,
     *                        TROCC,TROCC1,ISYINT,WORK,LWORK,INDSQ,
     *                        LENSQ,ISYMIB,IB,ISYMID,ID)
C
C     Kasper Hald, Fall 2001. 
C     Nearly identical to cc3_conocc by
C     Henrik Koch and Alfredo Sanchez.         Dec 1994
C     Ove Christiansen 9-1-1996
C
C     Set up combinations of S's and contract with integrals.
C
C     General symmetry: ISYMIM is symmetry of SMAT and TMAT intermediates.
C                       (including isymib*isymid)
C                       ISYINT is symmetry of integrals in TROCC and TROCC1.
C                       ISYRES = ISYMIM*ISYINT
C
#include <implicit.h>
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION OMEGA2(*),RMAT1(*),RMAT2(*),SMAT(*),TMAT(*)
      DIMENSION TROCC(*),TROCC1(*),WORK(LWORK),INDSQ(LENSQ,6)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CCFOP_CONOCC')
C
      IF (LWORK .LT. LENSQ) THEN
         CALL QUIT('Insufficient core in CONOCC')
      ENDIF
C
      ISYRES = MULD2H(ISYMIM,ISYINT)
C
C-------------------------
C     First occupied term.
C-------------------------
C
      C = ID
      B = IB
C
      ISYMC = ISYMID
      ISYMB = ISYMIB
C
      ISYMBC = MULD2H(ISYMB,ISYMC)
      JSAIKL = MULD2H(ISYMBC,ISYMIM)
C
      LENGTH = NCKIJ(JSAIKL)
C
C----------------------------------
C     Setup combinations of smat's.
C----------------------------------
C
      DO I = 1,LENGTH
C
C         TMAT(I) =       SMAT(I) 
C     *             - TWO*SMAT(INDSQ(I,3))
C     *             +     SMAT(INDSQ(I,4))
C
         TMAT(I) =  - SMAT(INDSQ(I,3))
C
      ENDDO
C
C----------------------------------
C     Symmetry sorting if symmetry.
C----------------------------------
C
      IF (NSYM .GT. 1) THEN
         CALL CC_GATHER(LENGTH,WORK,TMAT,INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK,1,TMAT,1)
      ENDIF
C
      IF (IPRINT .GT. 55) THEN
         XTMAT = DDOT(NCKIJ(JSAIKL),TMAT,1,TMAT,1)
         WRITE(LUPRI,*) 'In CC3_CONOCC: 1. Norm of TMAT = ',XTMAT
      ENDIF
C
C-----------------------
C     First contraction.
C-----------------------
C
      DO ISYMJ = 1,NSYM
C
         ISYMBJ = MULD2H(ISYMB,ISYMJ)
         ISYMAI = MULD2H(ISYMBJ,ISYRES)
         ISYMKL = MULD2H(JSAIKL,ISYMAI)
         ISYKLJ = MULD2H(ISYMKL,ISYMJ)
C
         NTOTAI = MAX(NT1AM(ISYMAI),1)
         NTOTKL = MAX(NMATIJ(ISYMKL),1)
C
         KOFF1  = ISAIKL(ISYMAI,ISYMKL) + 1
         KOFF2  = ISJIKA(ISYKLJ,ISYMC)
     *          + NMAJIK(ISYKLJ)*(C - 1)
     *          + ISJIK(ISYMKL,ISYMJ) + 1
         KOFF3  = ISAIK(ISYMAI,ISYMJ) + 1
C
         CALL DGEMM('N','N',NT1AM(ISYMAI),NRHF(ISYMJ),NMATIJ(ISYMKL),
     *              -ONE,TMAT(KOFF1),NTOTAI,TROCC(KOFF2),NTOTKL,
     *              ONE,RMAT2(KOFF3),NTOTAI)
C
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XRMAT = DDOT(NCKI(ISYRES),RMAT2,1,RMAT2,1)
         WRITE(LUPRI,*) 'In CC3_CONOCC: Norm of RMAT2 =  ',XRMAT
      ENDIF
C
C--------------------------
C     Second occupied term.
C--------------------------
C
      B = ID
      C = IB
C
      ISYMB = ISYMID
      ISYMC = ISYMIB
C
      ISYMBC = MULD2H(ISYMB,ISYMC)
      JSAIKL = MULD2H(ISYMBC,ISYMIM)
C
      LENGTH = NCKIJ(JSAIKL)
C
C----------------------------------
C     Setup combinations of smat's.
C----------------------------------
C
      DO I = 1,LENGTH
C
C         TMAT(I) = - TWO*SMAT(I) 
C     *             +     SMAT(INDSQ(I,3))
C     *             +     SMAT(INDSQ(I,5))
C
         TMAT(I) = - SMAT(I)
C
      ENDDO
C
C----------------------------------
C     Symmetry sorting if symmetry.
C----------------------------------
C
      IF (NSYM .GT. 1) THEN
         CALL CC_GATHER(LENGTH,WORK,TMAT,INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK,1,TMAT,1)
      ENDIF
C
      IF (IPRINT .GT. 55) THEN
         XTMAT = DDOT(NCKIJ(JSAIKL),TMAT,1,TMAT,1)
         WRITE(LUPRI,*) 'In CC3_CONOCC: 2. Norm of TMAT = ',XTMAT
      ENDIF
C
C------------------------
C     Second contraction.
C------------------------
C
      DO ISYMJ = 1,NSYM
C
         ISYMBJ = MULD2H(ISYMB,ISYMJ)
         ISYMAI = MULD2H(ISYMBJ,ISYRES)
         ISYMKL = MULD2H(JSAIKL,ISYMAI)
         ISYKLJ = MULD2H(ISYMKL,ISYMJ)
C
         NTOTAI = MAX(NT1AM(ISYMAI),1)
         NTOTKL = MAX(NMATIJ(ISYMKL),1)
C
         KOFF1  = ISAIKL(ISYMAI,ISYMKL) + 1
         KOFF2  = ISJIKA(ISYKLJ,ISYMC)
     *          + NMAJIK(ISYKLJ)*(C - 1)
     *          + ISJIK(ISYMKL,ISYMJ) + 1
         KOFF3  = ISAIK(ISYMAI,ISYMJ) + 1
C
         CALL DGEMM('N','N',NT1AM(ISYMAI),NRHF(ISYMJ),NMATIJ(ISYMKL),
     *              -ONE,TMAT(KOFF1),NTOTAI,TROCC(KOFF2),NTOTKL,
     *              ONE,RMAT1(KOFF3),NTOTAI)
C
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XRMAT = DDOT(NCKI(ISYRES),RMAT1,1,RMAT1,1)
         WRITE(LUPRI,*) 'In CC3_CONOCC: Norm of RMAT1 =  ',XRMAT
      ENDIF
C
C-------------------------
C     Third occupied term.
C-------------------------
C
      A = ID
      B = IB
C
      ISYMA = ISYMID
      ISYMB = ISYMIB
C
      ISYMAB = MULD2H(ISYMA,ISYMB)
      JSCKLI = MULD2H(ISYMAB,ISYMIM)
C
      LENGTH = NCKIJ(JSCKLI)
C
C----------------------------------
C     Setup combinations of smat's.
C----------------------------------
C
      DO I = 1,LENGTH
C
C         TMAT(I) =       SMAT(INDSQ(I,5)) 
C     *             - TWO*SMAT(INDSQ(I,2))
C     *             +     SMAT(INDSQ(I,3))
C
         TMAT(I) = - SMAT(INDSQ(I,2))
C
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XTMAT = DDOT(NCKIJ(JSAIKL),TMAT,1,TMAT,1)
         WRITE(LUPRI,*) 'In CC3_CONOCC: 3. Norm of TMAT = ',XTMAT
      ENDIF
C
C-----------------------
C     Third contraction.
C-----------------------
C
      DO ISYMJ = 1,NSYM
C
         ISYMBJ = MULD2H(ISYMB,ISYMJ)
         ISYMAI = MULD2H(ISYMBJ,ISYRES)
         ISYMI  = MULD2H(ISYMAI,ISYMA)
         ISYCKL = MULD2H(ISYMI,JSCKLI)
C
         IF (LWORK .LT. NRHF(ISYMI)*NRHF(ISYMJ)) THEN
            CALL QUIT('Insufficient memory in CCSDT_CONOCC')
         END IF
C
         NTOCKL = MAX(NCKI(ISYCKL),1)
         NRHFI  = MAX(NRHF(ISYMI),1)
C
         KOFF1  = ISAIKJ(ISYCKL,ISYMI) + 1
         KOFF2  = ISAIKJ(ISYCKL,ISYMJ) + 1
         KOFF3  = 1
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NCKI(ISYCKL),
     *              ONE,TMAT(KOFF1),NTOCKL,TROCC1(KOFF2),NTOCKL,
     *              ZERO,WORK(KOFF3),NRHFI)
C
         DO J = 1,NRHF(ISYMJ)
C
            NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
C
            IF (ISYMAI.EQ.ISYMBJ) THEN
C
               DO I = 1,NRHF(ISYMI)
C
                  NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A
C
                  KOFF5 = IT2AM(ISYMAI,ISYMBJ)
     *                 + INDEX(NAI,NBJ)
C
                  KOFF6 = NRHF(ISYMI)*(J - 1) + I
C
                  IF (NAI .EQ. NBJ) WORK(KOFF6) = TWO*WORK(KOFF6)
C
                  OMEGA2(KOFF5) = OMEGA2(KOFF5) - WORK(KOFF6)
C
               ENDDO
C
            ELSE IF (ISYMAI .LT. ISYMBJ) THEN
C
               DO I = 1,NRHF(ISYMI)
C
                  NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A
C
                  KOFF5 = IT2AM(ISYMAI,ISYMBJ)
     *                  + NT1AM(ISYMAI)*(NBJ-1) + NAI
C
                  KOFF6 = NRHF(ISYMI)*(J - 1) + I
                  OMEGA2(KOFF5) = OMEGA2(KOFF5) - WORK(KOFF6)
C
               ENDDO
C
            ELSE IF (ISYMBJ .LT. ISYMAI) THEN
C
               DO I = 1,NRHF(ISYMI)
C
                  NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A
C
                  KOFF5 = IT2AM(ISYMAI,ISYMBJ)
     *                  + NT1AM(ISYMBJ)*(NAI-1) + NBJ
C
                  KOFF6 = NRHF(ISYMI)*(J - 1) + I
                  OMEGA2(KOFF5) = OMEGA2(KOFF5) - WORK(KOFF6)
C
               ENDDO
C
            ENDIF
C
         ENDDO
C
      ENDDO
C
      CALL QEXIT('CCFOP_CONOCC')
C
      RETURN
      END
C------------------------------------------------------------------------- Sonia
C  /* Deck ccsdpt_dens2 */
      SUBROUTINE CCSDPT_DENS2(T1AM,ISYMT1,T2TP,ISYMT2,MODEL,
     *                        L1AM,ISYML1,L2TP,ISYML2,WORK,LWORK,
     *                        LUDELD,FNDELD,LUCKJD,FNCKJD,LUDKBC,FNDKBC,
     *                        LUTOC,FNTOC,LU3VI,FN3VI,LUDKBC3,FNDKBC3,
     *                        LU3FOP,FN3FOP,LU3FOP2,FN3FOP2,
     *                        LU3FOPX,FN3FOPX,LU3FOP2X,FN3FOP2X)
C
C     Written by K. Hald, Fall 2001.
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
C
#if defined (SYS_CRAY)
      REAL T1AM(*), T2TP(*)
      REAL L1AM(*), L2TP(*)
      REAL WORK(LWORK), ONE
      REAL TITRAN, TISORT, TISMAT, TIQMAT, TIOME1
      REAL TICONV, TICONO, RHO1N, RHO2N
      REAL XT2TP, DDOT, XIAJB, XINT, XTROC, XTROC1, XTROC0
      REAL XTRVI0, XTRVI2, XTRVI3, XTRVI, XTRVI1, XDIA
      REAL XSMAT, XTMAT, XQMAT, XRMAT, ZERO, TWO, HALF
      REAL DTIME
#else
      DOUBLE PRECISION T1AM(*), T2TP(*)
      DOUBLE PRECISION L1AM(*), L2TP(*)
      DOUBLE PRECISION WORK(LWORK), ONE
      DOUBLE PRECISION TITRAN, TISORT, TISMAT, TIQMAT, TIOME1
      DOUBLE PRECISION TICONV, TICONO, RHO1N, RHO2N
      DOUBLE PRECISION XT2TP, DDOT, XIAJB, XINT, XTROC, XTROC1, XTROC0
      DOUBLE PRECISION XTRVI0, XTRVI2, XTRVI3, XTRVI, XTRVI1, XDIA
      DOUBLE PRECISION XSMAT, XTMAT, XQMAT, XRMAT, ZERO, TWO, HALF
      DOUBLE PRECISION DTIME
#endif
C
      LOGICAL   C3LRSV, CC1ASV, CC1BSV, LDEBUG
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
      PARAMETER (LDEBUG = .FALSE.)
C
      CALL QENTER('CCSDPT_DENS2')
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
      TICONO = 0.0D0
      TICONV = 0.0D0
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
C     Calculate the lamda matrices for cc3
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
                     CALL CCFOP_QMAT(0.0D0,T1AM,ISYMT1,WORK(KT2TCME),
     *                               ISYMT2,WORK(KTMAT),WORK(KFCKBA),
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
                     CALL CCFOP_DENVIR(WORK(KVIR1),WORK(KVIR2),
     *                                 WORK(KSMAT),WORK(KQMAT),
     *                                 WORK(KTMAT),ISYMIM,
     *                                 WORK(KT2TCME),ISYMT2,WORK(KEND4),
     *                                 LWRK4,WORK(KINDSQ),LENSQ,
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
                     CALL CCFOP_DENVIR(WORK(KVIR3),WORK(KVIR4),
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
                     CALL CCFOP_DENOCC(WORK(KOCC1),WORK(KSMAT),
     *                                 WORK(KQMAT),WORK(KTMAT),ISYMIM,
     *                                 WORK(KT2TCME),ISYMT2,WORK(KEND4),
     *                                 LWRK4,WORK(KINDSQ),LENSQ,
     *                                 ISYMB,B,ISYMD,D,1)
C
                     CALL CCFOP_DENOCC(WORK(KOCC2),WORK(KSMAT2),
     *                                 WORK(KQMAT2),WORK(KTMAT),ISYMIM,
     *                                 T2TP,ISYMT2,WORK(KEND4),
     *                                 LWRK4,WORK(KINDSQ),LENSQ,
     *                                 ISYMB,B,ISYMD,D,2)
C
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
C
               CALL DEN_AIBSORT(WORK(KVIR1),WORK(KEND4),ISAIJ1)
C
               CALL DEN_AIBSORT(WORK(KVIR2),WORK(KEND4),ISAIJ1)
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
               CALL DEN_AIBSORT(WORK(KVIR3),WORK(KEND4),ISAIJ1)
C
               CALL DEN_AIBSORT(WORK(KVIR4),WORK(KEND4),ISAIJ1)
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
      IF ((IPRINT .GT. 55)) THEN
         RHO1N = DDOT(NCKIJ(ISYRES),WORK(KOCC1),1,WORK(KOCC1),1)
         WRITE(LUPRI,*) 'Norm of OCC1 (iajk) (before dens1to2) =',RHO1N
      ENDIF
C
      IF ((.NOT. CC3) .AND. (RELORB)) THEN
         CALL DENS1TO2(WORK(KOMG1),WORK(KOCC1),ISYRES)
      ENDIF
C
      IF ((IPRINT .GT. 55)) THEN
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
      IF ((IPRINT .GT. 55)) THEN
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
         IF ((IPRINT .GT. 55)) THEN
            RHO1N = DDOT(NRHFT,WORK(KKAPII),1,WORK(KKAPII),1)
            WRITE(LUPRI,*) 'Norm of KAPII : ',RHO1N
            RHO1N = DDOT(NVIRT,WORK(KKAPAA),1,WORK(KKAPAA),1)
            WRITE(LUPRI,*) 'Norm of KAPAA : ',RHO1N
         ENDIF
C
         LUKAPAB = -1
         LUKAPIJ = -1
         FNKAPAB = 'KAPAB'
         FNKAPIJ = 'KAPIJ'
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
         CALL SYMMBACK(MODEL,LUIAJK,FNDIAJK,LUAIJK,FNDAIJK,
     *                 LUABI1,FNDABI1,LUABI3,FNDABI3,
     *                 LUPTIAJB,FNDIAJB,
     *                 ISYRES,WORK(KEND4),LWRK4)
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
         WRITE(LUPRI,1) 'CC3_CONV  : ',TICONV
         WRITE(LUPRI,1) 'CC3_CONO  : ',TICONO
         WRITE(LUPRI,1) 'CC3_OME1  : ',TIOME1
         WRITE(LUPRI,*)
      END IF
C
C-------------
C     End
C-------------
C
      CALL QEXIT('CCSDPT_DENS2')
C
      RETURN
C
    1 FORMAT(7X,'Time used in',2X,A12,F12.2,' seconds')
C
      END
C  /* Deck ccfop_onel */
      SUBROUTINE CCFOP_ONEL(OMEGA2,RMAT1,RMAT2,T1AM,SMAT,TMAT,
     *                      ISYMIM,ISYINT,INDSQ,LENSQ,WORK,LWORK,
     *                      ISYMIB,IB,ISYMID,ID)
C
C     Kasper Hald, Fall 2001.
C     Based on cc3_onel by
C     Henrik Koch and Alfredo Sanchez.         Dec 1994
C     Ove Christiansen 9-1-1996:
C
C     Calculate 2 electon density iajb in CCSD(T)
C
C
C     General symmetry: ISYMIM is symmetry of SMAT and TMAT 
C                       intermdiates.(incl isymd,isymb)
C                       ISYINT is symmetry of T1AM
C                       ISYRES = ISYMIM*ISYINT
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INTEGER ISYMIM, ISYINT, LENSQ, LWORK, ISYMIB, IB, ISYMID, ID
      INTEGER ISYRES, ISYMB, ISYMC, ISYMK, ISYMBC, JSAIKJ
      INTEGER ISYAIJ, ISYMCK, LENGTH, NCK, KOFF1, NTOAIJ
      INTEGER NTOTC, ISYMBK, NBK, NTOTB, JSAKIJ, ISYMIJ, ISYMAK
      INTEGER NTOTAK, NTOTIJ, ISYMJ, ISYMI, ISYMCI, NBJ, NIJ, NCI, NCIBJ
      INTEGER ISYMBJ
      INTEGER INDEX, INDSQ(LENSQ,6)
C
#if defined (SYS_CRAY)
      REAL OMEGA2(*), RMAT1(*), RMAT2(*), T1AM(*), SMAT(*)
      REAL TMAT(*), WORK(LWORK)
      REAL ZERO, ONE, TWO
#else
      DOUBLE PRECISION OMEGA2(*), RMAT1(*), RMAT2(*), T1AM(*), SMAT(*)
      DOUBLE PRECISION TMAT(*), WORK(LWORK)
      DOUBLE PRECISION ZERO, ONE, TWO
#endif
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CCFOP_ONEL')
C
      ISYRES = MULD2H(ISYMIM,ISYINT)
C
      B = IB
      C = ID
C
      ISYMB = ISYMIB
      ISYMC = ISYMID
C
C----------------------------------
C     First contribution to Omega2.
C----------------------------------
C
      ISYMK  = MULD2H(ISYMC,ISYINT)
      ISYMBC = MULD2H(ISYMB,ISYMC)
      JSAIKJ = MULD2H(ISYMBC,ISYMIM)
      ISYAIJ = MULD2H(ISYMK,JSAIKJ)
      ISYMCK = MULD2H(ISYMC,ISYMK)
C
      LENGTH = NCKIJ(JSAIKJ)
C
      IF (LWORK .LT. LENGTH) THEN
         CALL QUIT('Not enough core in CCSDT_ONEL')
      END IF
C
      DO I = 1,LENGTH
         TMAT(I) =   SMAT(INDSQ(I,4)) 
     *             - SMAT(INDSQ(I,3))
      ENDDO
C
      NCK = IT1AM(ISYMC,ISYMK) + C
C
      KOFF1 = ISAIKJ(ISYAIJ,ISYMK) + 1
C
      NTOAIJ = MAX(NCKI(ISYAIJ),1)
      NTOTC  = MAX(NVIR(ISYMC),1)
C
      CALL DGEMV('N',NCKI(ISYAIJ),NRHF(ISYMK),ONE,TMAT(KOFF1),NTOAIJ,
     *           T1AM(NCK),NTOTC,ONE,RMAT2,1)
C
C-----------------------------------
C     Second contribution to Omega2.
C-----------------------------------
C
      ISYMK  = MULD2H(ISYMB,ISYINT)
      ISYMBC = MULD2H(ISYMB,ISYMC)
      JSAIKJ = MULD2H(ISYMBC,ISYMIM)
      ISYAIJ = MULD2H(ISYMK,JSAIKJ)
      ISYMBK = MULD2H(ISYMB,ISYMK)
C
      LENGTH = NCKIJ(JSAIKJ)
C
      IF (LWORK .LT. LENGTH) THEN
         CALL QUIT('Not enough core in CCFOP_ONEL')
      END IF
C
      DO I = 1,LENGTH
         TMAT(I) =   SMAT(INDSQ(I,5)) 
     *             - SMAT(I)
      ENDDO
C
      NBK = IT1AM(ISYMB,ISYMK) + B
C
      KOFF1 = ISAIKJ(ISYAIJ,ISYMK) + 1
C
      NTOAIJ = MAX(NCKI(ISYAIJ),1)
      NTOTB  = MAX(NVIR(ISYMB),1)
C
      CALL DGEMV('N',NCKI(ISYAIJ),NRHF(ISYMK),ONE,TMAT(KOFF1),NTOAIJ,
     *           T1AM(NBK),NTOTB,ONE,RMAT1,1)
C
C----------------------------------
C     Third contribution to Omega2.
C----------------------------------
C
      ISYMBC = MULD2H(ISYMB,ISYMC)
      JSAKIJ = MULD2H(ISYMBC,ISYMIM)
      ISYMIJ = MULD2H(ISYMBC,ISYRES)
      ISYMAK = MULD2H(JSAKIJ,ISYMIJ)
C
      LENGTH = NCKIJ(JSAKIJ)
C
      IF (LWORK .LT. LENGTH) THEN
         CALL QUIT('Not enough core in CCSDT_ONEL')
      END IF
C
      DO I = 1,LENGTH
         TMAT(I) =   SMAT(INDSQ(I,1)) 
     *             - SMAT(I)
      ENDDO
C
C     Symmetry sorting if symmetry
C     ----------------------------
C
      IF (NSYM .GT. 1) THEN
         CALL CC_GATHER(LENGTH,WORK,TMAT,INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK,1,TMAT,1)
      ENDIF
C
      NTOTAK = MAX(NT1AM(ISYMAK),1)
      NTOTIJ = MAX(NMATIJ(ISYMIJ),1)
C
      KOFF1 = ISAIKL(ISYMAK,ISYMIJ) + 1
C
      CALL DGEMV('T',NT1AM(ISYMAK),NMATIJ(ISYMIJ),ONE,TMAT(KOFF1),
     *           NTOTAK,T1AM,1,ZERO,WORK,1)
C
      DO ISYMJ = 1,NSYM
C
         ISYMI  = MULD2H(ISYMIJ,ISYMJ)
C
         ISYMBJ = MULD2H(ISYMB,ISYMJ)
         ISYMCI = MULD2H(ISYMC,ISYMI)
C
         DO J = 1,NRHF(ISYMJ)
C
            NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
C
            IF (ISYMCI .EQ. ISYMBJ) THEN
C
               DO I = 1,NRHF(ISYMI)
C
                  NIJ = IMATIJ(ISYMI,ISYMJ) + NRHF(ISYMI)*(J - 1) + I
                  NCI = IT1AM(ISYMC,ISYMI)  + NVIR(ISYMC)*(I - 1) + C
C
                  IF (NCI .EQ. NBJ) WORK(NIJ) = TWO*WORK(NIJ)
C
                  NCIBJ = IT2AM(ISYMCI,ISYMBJ) + INDEX(NCI,NBJ)
C
                  OMEGA2(NCIBJ) = OMEGA2(NCIBJ) + WORK(NIJ)
C
               ENDDO
C
            ELSE IF (ISYMCI .LT. ISYMBJ) THEN
C
               DO I = 1,NRHF(ISYMI)
C
                  NIJ = IMATIJ(ISYMI,ISYMJ) + NRHF(ISYMI)*(J - 1) + I
                  NCI = IT1AM(ISYMC,ISYMI)  + NVIR(ISYMC)*(I - 1) + C
C
                  NCIBJ = IT2AM(ISYMCI,ISYMBJ)
     *                  + NT1AM(ISYMCI)*(NBJ-1) + NCI
C
                  OMEGA2(NCIBJ) = OMEGA2(NCIBJ) + WORK(NIJ)
C
               ENDDO
C
            ELSE IF (ISYMBJ .LT. ISYMCI) THEN
C
               DO I = 1,NRHF(ISYMI)
C
                  NIJ = IMATIJ(ISYMI,ISYMJ) + NRHF(ISYMI)*(J - 1) + I
                  NCI = IT1AM(ISYMC,ISYMI)  + NVIR(ISYMC)*(I - 1) + C
C
                  NCIBJ = IT2AM(ISYMBJ,ISYMCI)
     *                  + NT1AM(ISYMBJ)*(NCI-1) + NBJ
C
                  OMEGA2(NCIBJ) = OMEGA2(NCIBJ) + WORK(NIJ)
C
               ENDDO
C
            ENDIF
C
         ENDDO
C
      ENDDO
C
      CALL QEXIT('CCFOP_ONEL')
C
      RETURN
      END
C  /* Deck ccfop_denvir */
      SUBROUTINE CCFOP_DENVIR(RINTE1,RINTE2,SMAT,QMAT,TMAT,ISYMIM,
     *                        T2TCME,ISYMT2,WORK,LWORK,INDSQ,LENSQ,
     *                        ISYMB,B,ISYMD,D,IOPT)
C
C     Kasper Hald, Fall 2001.
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
      INTEGER ISYBIJ, ISYMBJ, NTOTA
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
      CALL QENTER('CCFOP_DENVIR')
C
C-----------------------------------------------
C     Sanity check and symmetry calculation.
C-----------------------------------------------
C
      IF (IOPT .NE. 1 .AND. IOPT .NE. 2) THEN
         CALL QUIT('Wrong IOPT in CCFOP_DENVIR')
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
                  NAI = IT1AM(ISYMA,ISYMI)
     *                + NVIR(ISYMA)*(I-1) + A
C
                  KOFF1 =  IT2SP(ISYBIJ,ISYMA)
     *                  +  NCKI(ISYBIJ)*(A - 1)
C     *                  +  ISAIK(ISYMBJ,ISYMI)
     *                  +  ICKI(ISYMBJ,ISYMI)
     *                  +  NT1AM(ISYMBJ)*(I-1)
     *                  +  IT1AM(ISYMB,ISYMJ)
     *                  +  B
C
                  KOFF2 =  ISAIK(ISYMAI,ISYMJ)
     *                  +  IT1AM(ISYMA,ISYMI)
     *                  +  NVIR(ISYMA)*(I-1)
     *                  +  A
C
                  CALL DCOPY(NRHF(ISYMJ),T2TCME(KOFF1),NVIR(ISYMB),
     *                       WORK(KOFF2),NT1AM(ISYMAI))
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
C------------------------------
C     Contract with T2
C------------------------------
C
      DO ISYMK = 1,NSYM
C
         ISYCIJ = MULD2H(ISCKIJ,ISYMK)
C
         DO ISYMC = 1, NSYM
C
            ISYMIJ = MULD2H(ISYCIJ,ISYMC)
            ISYMAB = MULD2H(ISYMT2,ISYMIJ)
            ISYMA  = MULD2H(ISYMAB,ISYMB)
            ISYMCK = MULD2H(ISYMC,ISYMK)
C
            NTOTCK = MAX(NT1AM(ISYMCK),1)
            NTOTIJ = MAX(NMATIJ(ISYMIJ),1)
            NTOTA  = MAX(NVIR(ISYMA),1)
C
            KOFF1 = ISAIK(ISYMAI,ISYMJ) + 1
            KOFF2 = ISAIKL(ISYMCK,ISYMIJ) + 1
            KOFF3 = ICKATR(ISYMCK,ISYMA)  + 1
C
            CALL DGEMM('N','T',NVIR(ISYMA),NT1AM(ISYMCK),
     *                 NMATIJ(ISYMIJ),TWO,WORK(KOFF1),NTOTA,
     *                 TMAT(KOFF2),NTOTCK,ONE,RINTE1(KOFF3),
     *                 NTOTA)
C
         ENDDO         ! ISYMC
      ENDDO            ! ISYMK
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
C
C------------------------------
C     Contract with T2
C------------------------------
C
      DO ISYMK = 1,NSYM
C
         ISYCIJ = MULD2H(ISCKIJ,ISYMK)
C
         DO ISYMC = 1, NSYM
C
            ISYMIJ = MULD2H(ISYCIJ,ISYMC)
            ISYMAB = MULD2H(ISYMT2,ISYMIJ)
            ISYMA  = MULD2H(ISYMAB,ISYMB)
            ISYMCK = MULD2H(ISYMC,ISYMK)
C
            NTOTCK = MAX(NT1AM(ISYMCK),1)
            NTOTIJ = MAX(NMATIJ(ISYMIJ),1)
            NTOTA  = MAX(NVIR(ISYMA),1)
C
            KOFF1 = ISAIK(ISYMAI,ISYMJ) + 1
            KOFF2 = ISAIKL(ISYMCK,ISYMIJ) + 1
            KOFF3 = ICKATR(ISYMCK,ISYMA)  + 1
C
            CALL DGEMM('N','T',NVIR(ISYMA),NT1AM(ISYMCK),
     *                 NMATIJ(ISYMIJ),TWO,WORK(KOFF1),NTOTA,
     *                 TMAT(KOFF2),NTOTCK,ONE,RINTE2(KOFF3),
     *                 NTOTA)
C
         ENDDO         ! ISYMC
      ENDDO            ! ISYMK
C
C
      CALL QEXIT('CCFOP_DENVIR')
C
      RETURN
      END
C  /* Deck t3_forbidden */
      SUBROUTINE T3_FORBIDDEN(SMAT,ISYMIM,ISYMB,B,ISYMD,D)
C
C     Written by Kasper Hald, Fall 2001.
C
C     Purpose : Remove the forbidden t3/t3-bar amplitudes.
C
      IMPLICIT NONE
C
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYMIM, ISYMB, ISYMD, ISYIJK, ISYMK, ISYMIJ
      INTEGER ISYBIJ, ISYMJ, ISYMI, ISYMBI, KOFF1, ISYMBD
      INTEGER ISYMAI, ISYMA
C
#if defined (SYS_CRAY)
      REAL SMAT(*), ZERO
#else
      DOUBLE PRECISION SMAT(*), ZERO
#endif
C
      PARAMETER (ZERO = 0.0D0)
C
      CALL QENTER('T3_FORBIDDEN')
C
C---------------------------------------------------------
C     If B and D are the same remove all amplitudes
C     having an A which is the same as B and D.
C---------------------------------------------------------
C
      IF ((ISYMB .EQ. ISYMD) .AND. (B .EQ. D) ) THEN
         ISYIJK = MULD2H(ISYMB,ISYMIM)
         DO ISYMK = 1, NSYM
            ISYMIJ = MULD2H(ISYMK,ISYIJK)
            ISYBIJ = MULD2H(ISYMIJ,ISYMB)
            DO ISYMJ = 1, NSYM
               ISYMI  = MULD2H(ISYMJ,ISYMIJ)
               ISYMBI = MULD2H(ISYMB,ISYMI)
C
               DO K = 1, NRHF(ISYMK)
               DO J = 1, NRHF(ISYMJ)
               DO I = 1, NRHF(ISYMI)
                  KOFF1 = ISAIKJ(ISYBIJ,ISYMK)
     *                  + NCKI(ISYBIJ)*(K - 1)
     *                  + ISAIK(ISYMBI,ISYMJ)
     *                  + NT1AM(ISYMBI)*(J-1)
     *                  + IT1AM(ISYMB,ISYMI)
     *                  + NVIR(ISYMB)*(I-1)
     *                  + B
C
                  SMAT(KOFF1)  = ZERO
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
C     occupied indices.
C----------------------------------------------------------
C
      ISYMBD = MULD2H(ISYMB,ISYMD)
      ISYMAI = MULD2H(ISYMIM,ISYMBD)
      DO ISYMA = 1, NSYM
         ISYMI = MULD2H(ISYMA,ISYMAI)
         DO A = 1, NVIR(ISYMA)
            DO I = 1, NRHF(ISYMI)
               KOFF1 = ISAIKJ(ISYMA,ISYMI)
     *               + NCKI(ISYMA)*(I - 1)
     *               + ISAIK(ISYMAI,ISYMI)
     *               + NT1AM(ISYMAI)*(I-1)
     *               + IT1AM(ISYMA,ISYMI)
     *               + NVIR(ISYMA)*(I-1)
     *               + A
C
               SMAT(KOFF1)  = ZERO
C
            ENDDO
         ENDDO
      ENDDO
C
C-----------------------
C     End.
C-----------------------
C
      CALL QEXIT('T3_FORBIDDEN')
C
      RETURN
      END
C  /* Deck den_aibsort */
      SUBROUTINE DEN_AIBSORT(VIRREAL,VIRTMP,ISYVIR)
C
C     Written by Kasper Hald, 2001.
C
C     Purpose : Sort the two electron densities d(abi) -> d(aib)
C               where the densities have a constant C.
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYVIR, ISYMB, ISYMA, ISYMI, ISYMAI, ISYMBI
      INTEGER KOFF1, KOFF2, ISYMAB
C
#if defined (SYS_CRAY)
      REAL VIRREAL(*), VIRTMP(*)
#else
      DOUBLE PRECISION VIRREAL(*), VIRTMP(*), tmp
#endif
C
      CALL QENTER('DEN_AIBSORT')
C
C----------------------------------------
C     Sort matrix.
C----------------------------------------
C
      DO ISYMB = 1, NSYM
         ISYMAI = MULD2H(ISYMB,ISYVIR)
         DO ISYMA = 1, NSYM
            ISYMAB = MULD2H(ISYMA,ISYMB)
            ISYMI  = MULD2H(ISYMAI,ISYMA)
            DO B = 1, NVIR(ISYMB)
               DO A = 1, NVIR(ISYMA)
                  DO I = 1, NRHF(ISYMI)
                     KOFF1 = ICKASR(ISYMAB,ISYMI)
     *                     + NMATAB(ISYMAB)*(I-1)
     *                     + IMATAB(ISYMA,ISYMB)
     *                     + NVIR(ISYMA)*(B-1)
     *                     + 1
C
                     KOFF2 = ICKATR(ISYMAI,ISYMB)
     *                     + NT1AM(ISYMAI)*(B-1)
     *                     + IT1AM(ISYMA,ISYMI)
     *                     + NVIR(ISYMA)*(I-1)
     *                     + 1
C
                     CALL DCOPY(NVIR(ISYMA),VIRREAL(KOFF1),1,
     *                          VIRTMP(KOFF2),1)
C
                  ENDDO
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
      CALL QEXIT('DEN_AIBSORT')
C
      RETURN
      END
C  /* Deck denpt */
      SUBROUTINE DENSTORE(VIR1,LUVIR1,FNVIR1,VIR2,LUVIR2,FNVIR2,ISYRES)
C
C     Written by K. Hald, Fall 2001.
C
C     Purpose : Get the two different densities from file
C               and sum them up with correct index.
C               Store the total density on file.
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER LUVIR1, LUVIR2, ISYRES
      INTEGER ISYMD, ISYAIB, ISYMB, ISYMAI, ISYAID, KOFF1, KOFF2
C
#if defined (SYS_CRAY)
      DOUBLE PRECISION VIR1(*), VIR2(*), ONE
#else
      DOUBLE PRECISION VIR1(*), VIR2(*), ONE
#endif
C
      CHARACTER*(*) FNVIR1, FNVIR2
C
      PARAMETER (ONE = 1.0D0)
C
      CALL QENTER('DENSTORE')
C
      DO ISYMD = 1, NSYM
         ISYAIB = MULD2H(ISYMD,ISYRES)
         DO ISYMB = 1, NSYM
            ISYMAI  = MULD2H(ISYAIB,ISYMB)
            ISYAID  = MULD2H(ISYMAI,ISYMD)
            DO D = 1, NVIR(ISYMD)
               DO B = 1, NVIR(ISYMB)
C
                  KOFF1 = ICKBD(ISYAIB,ISYMD)
     *                  + NCKATR(ISYAIB)*(D-1)
     *                  + ICKATR(ISYMAI,ISYMB)
     *                  + NT1AM(ISYMAI)*(B-1)
     *                  + 1
C
                  CALL GETWA2(LUVIR2,FNVIR2,VIR2,KOFF1,
     *                        NT1AM(ISYMAI))
C
                  KOFF2 = ICKBD(ISYAID,ISYMB)
     *                  + NCKATR(ISYAID)*(B-1)
     *                  + ICKATR(ISYMAI,ISYMD)
     *                  + NT1AM(ISYMAI)*(D-1)
     *                  + 1
C
                  CALL GETWA2(LUVIR1,FNVIR1,VIR1,
     *                        KOFF2,NT1AM(ISYMAI))
C
                  CALL DAXPY(NT1AM(ISYMAI),ONE,VIR2,1,
     *                       VIR1,1)
C
                  CALL PUTWA2(LUVIR1,FNVIR1,VIR1,
     *                        KOFF2,NT1AM(ISYMAI))
C
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C
C-----------------------
C     End.
C-----------------------
C
      CALL QEXIT('DENSTORE')
C
      RETURN
      END
C  /* Deck ccfop_denocc */
      SUBROUTINE CCFOP_DENOCC(OCC,SMAT,QMAT,TMAT,ISYMIM,T2AM,ISYMT2,
     *                        WORK,LWORK,INDSQ,LENSQ,ISYMB,B,
     *                        ISYMD,D,IOPT)
C
C     Written by Kasper Hald, Fall 2001.
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
      INTEGER ISYMEK
C
#if defined (SYS_CRAY)
      REAL OCC(*), SMAT(*), QMAT(*), TMAT(*), T2AM(*)
      REAL WORK(LWORK), TWO, ONE, HALF
#else
      DOUBLE PRECISION OCC(*), SMAT(*), QMAT(*), TMAT(*), T2AM(*)
      DOUBLE PRECISION WORK(LWORK), TWO, ONE, HALF
#endif
C
      PARAMETER (ONE = 1.0D0, TWO = 2.0D0, HALF = 0.5D0)
C
      CALL QENTER('CCFOP_DENOCC')
C
C--------------------------------------
C     Symmetries and sanity check.
C--------------------------------------
C
      IF (IOPT .NE. 1 .AND. IOPT .NE. 2) THEN
         CALL QUIT('Wrong IOPT in CCFOP_DENOCC')
      ENDIF
C
      ISYMBD = MULD2H(ISYMB,ISYMD)
      ISELJI = MULD2H(ISYMIM,ISYMBD)
      ISYELK = MULD2H(ISYMT2,ISYMD)
      ISYIJK = MULD2H(ISYELK,ISELJI)
C
      IF (LWORK .LT. NCKI(ISYELK)) THEN
         CALL QUIT('Not enough memory in CCFOP_DENOCC')
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
C
      DO ISYMK = 1, NSYM
         ISYMEL = MULD2H(ISYELK,ISYMK)
         ISYMIJ = MULD2H(ISYIJK,ISYMK)
C
         NTOTEL = MAX(NT1AM(ISYMEL),1)
         NTOTK  = MAX(NRHF(ISYMK),1)
C
         KOFF1  = IT2SP(ISYELK,ISYMD)
     *          + NCKI(ISYELK)*(D-1)
     *          + ISAIK(ISYMEL,ISYMK) + 1
         KOFF2  = ISAIKL(ISYMEL,ISYMIJ) + 1
         KOFF3  = I3OVIR(ISYIJK,ISYMB)
     *          + NMAIJK(ISYIJK)*(B-1)
     *          + IMAIJK(ISYMIJ,ISYMK)
     *          + 1
!         KOFF3  = I3OVIR(ISYIJK,ISYMB)
!     *          + NMAIJK(ISYIJK)*(B-1)
!     *          + Itestd(ISYMk,ISYMIJ)
!     *          + 1
C
         CALL DGEMM('T','N',NRHF(ISYMK),NMATIJ(ISYMIJ),
     *              NT1AM(ISYMEL),TWO,T2AM(KOFF1),NTOTEL,
     *              TMAT(KOFF2),NTOTEL,ONE,OCC(KOFF3),
     *              NTOTK)
      ENDDO
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
C--------------------------------------------
C     Contract with S and Q intermediates.
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
C
      DO ISYMK = 1, NSYM
         ISYMEL = MULD2H(ISYELK,ISYMK)
         ISYMIJ = MULD2H(ISYIJK,ISYMK)
C
         NTOTEL = MAX(NT1AM(ISYMEL),1)
         NTOTK  = MAX(NRHF(ISYMK),1)
C
         KOFF1  = ISAIK(ISYMEL,ISYMK) + 1
         KOFF2  = ISAIKL(ISYMEL,ISYMIJ) + 1
         KOFF3  = I3OVIR(ISYIJK,ISYMB)
     *          + NMAIJK(ISYIJK)*(B-1)
     *          + IMAIJK(ISYMIJ,ISYMK)
     *          + 1
!         KOFF3  = I3OVIR(ISYIJK,ISYMB)
!     *          + NMAIJK(ISYIJK)*(B-1)
!     *          + Itestd(ISYMK,ISYMIJ)
!     *          + 1
C
         CALL DGEMM('T','N',NRHF(ISYMK),NMATIJ(ISYMIJ),
     *              NT1AM(ISYMEL),TWO,WORK(KOFF1),NTOTEL,
     *              TMAT(KOFF2),NTOTEL,ONE,OCC(KOFF3),
     *              NTOTK)
      ENDDO
C
C-----------------------
C     End.
C-----------------------
C
      CALL QEXIT('CCFOP_DENOCC')
C
      RETURN
      END
C  /* Deck ccfop_oned */
      SUBROUTINE CCFOP_ONED(OMEGA1,T2AM,ISYMT2,SMAT,TMAT,ISYMIM,INDSQ,
     *                      LENSQ,WORK,LWORK,ISYMIB,IB,ISYMID,ID)
C
C     Written by K. Hald, Fall 2001.
C
C     Based on cc3_onel by
C     Henrik Koch and Alfredo Sanchez.         Dec 1994
C     Ove Christiansen 9-1-1996:
C
C     Calculate the contributions to Omega1 in CCSD(T) unrelaxed f.o.p.
C
C     omega1(ai) = (t^{dea}_{lmi} - t^{dea}_{lim}) * (t^{* (0)}_{dl,em})
C
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INTEGER ISYMT2, ISYMIM, LENSQ, LWORK, ISYMIB, IB, ISYMID, ID
      INTEGER ISYRES, ISYMB, ISYMC, ISYMI, ISYAKJ, ISYMJ, ISYMBJ, ISYMAK
      INTEGER NBJ, NAK, NAKBJ, NAKJ, NTOTC, NTOAKJ, KOFF1, KOFF2
      INTEGER ISYMBC, JSAIKJ, LENGTH, ISYCKJ, ISYMAI, ISYMKJ, ISYMCK
      INTEGER NKJ, NCK, NCKBJ, NTOTAI, ISYMK, JSAKIJ
      INTEGER INDEX, INDSQ(LENSQ,6)
C
#if defined (SYS_CRAY)
      REAL OMEGA1(*), T2AM(*), SMAT(*), TMAT(*), WORK(LWORK)
      REAL ZERO, ONE, TWO
#else
      DOUBLE PRECISION OMEGA1(*), T2AM(*), SMAT(*), TMAT(*), WORK(LWORK)
      DOUBLE PRECISION ZERO, ONE, TWO
#endif
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CCFOP_ONED')
C
      ISYRES = MULD2H(ISYMIM,ISYMT2)
C
      B = IB
      C = ID
C
      ISYMB = ISYMIB
      ISYMC = ISYMID
C
C-----------------------------------
C     First SMAT / TMAT magic.
C-----------------------------------
C
      ISYMBC = MULD2H(ISYMB,ISYMC)
      JSAIKJ = MULD2H(ISYMBC,ISYMIM)
      LENGTH = NCKIJ(JSAIKJ)
C
      DO I = 1,LENGTH
C
         TMAT(I) =   SMAT(INDSQ(I,3))
     *             - SMAT(INDSQ(I,4))
C
      ENDDO
C
C----------------------------------
C     First contribution to Omega1.
C----------------------------------
C
      ISYMI  = MULD2H(ISYMC,ISYRES)
      ISYAKJ = MULD2H(ISYMB,ISYMT2)
C
      IF (NRHF(ISYMI) .NE. 0) THEN
C
         IF (LWORK .LT. NCKI(ISYAKJ)) THEN
            CALL QUIT('Not enough core in CCFOP_ONED')
         END IF
C
C        Construct M(ak,j) = T(ak,bj)
C        ---------------------------
C
         DO ISYMJ = 1,NSYM
C
            ISYMAK = MULD2H(ISYMJ,ISYAKJ)
C
            DO J = 1,NRHF(ISYMJ)
C
               NAKBJ = IT2SP(ISYAKJ,ISYMB)
     *               + NCKI(ISYAKJ)*(B-1)
     *               + ICKI(ISYMAK,ISYMJ)
     *               + NT1AM(ISYMAK)*(J - 1)
     *               + 1
C
               NAKJ  = ICKI(ISYMAK,ISYMJ)
     *               + NT1AM(ISYMAK)*(J - 1) 
     *               + 1
C
               CALL DCOPY(NT1AM(ISYMAK),T2AM(NAKBJ),1,
     *                    WORK(NAKJ),1)
C
            ENDDO
         ENDDO
C
         NTOTC  = MAX(NVIR(ISYMC),1)
         NTOAKJ = MAX(NCKI(ISYAKJ),1)
C
         KOFF1 = ISAIKJ(ISYAKJ,ISYMI) + 1
         KOFF2 = IT1AM(ISYMC,ISYMI) + C
C
         CALL DGEMV('T',NCKI(ISYAKJ),NRHF(ISYMI),ONE,TMAT(KOFF1),
     *              NTOAKJ,WORK,1,ONE,OMEGA1(KOFF2),NTOTC)
C
      ENDIF
C
C---------------------------------------
C     Second contribution to Omega1
C---------------------------------------
C
      ISYMBC = MULD2H(ISYMB,ISYMC)
      JSAIKJ = MULD2H(ISYMBC,ISYMIM)
C
      LENGTH = NCKIJ(JSAIKJ)
C
      IF (LWORK .LT. LENGTH) THEN
         CALL QUIT('Not enough core in CCFOP_ONED')
      END IF
C
      DO I = 1,LENGTH
         TMAT(I) =   SMAT(I)
     *             - SMAT(INDSQ(I,5))
      ENDDO
C
      IF (NSYM .GT. 1) THEN
         CALL CC_GATHER(LENGTH,WORK,TMAT,INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK,1,TMAT,1)
      ENDIF
C
      ISYMKJ = MULD2H(ISYMBC,ISYMT2)
      ISYMAI = ISYRES
C
C        Construct M(k,j) = T(ck,bj)
C        ---------------------------
C
      DO ISYMJ = 1,NSYM
C
         ISYMK  = MULD2H(ISYMJ,ISYMKJ)
         ISYMCK = MULD2H(ISYMC,ISYMK)
         ISYCKJ = MULD2H(ISYMCK,ISYMJ)
C
         DO J = 1,NRHF(ISYMJ)
C
            NKJ   = IMATIJ(ISYMK,ISYMJ)
     *            + NRHF(ISYMK)*(J - 1) 
     *            + 1
C
            NCKBJ = IT2SP(ISYCKJ,ISYMB)
     *            + NCKI(ISYCKJ)*(B-1)
     *            + ICKI(ISYMCK,ISYMJ)
     *            + NT1AM(ISYMCK)*(J-1)
     *            + IT1AM(ISYMC,ISYMK)
     *            + C
C
            CALL DCOPY(NRHF(ISYMK),T2AM(NCKBJ),NVIR(ISYMC),WORK(NKJ),1)
C
         ENDDO
      ENDDO
C
      NTOTAI = MAX(NT1AM(ISYMAI),1)
C
      KOFF1 = ISAIKL(ISYMAI,ISYMKJ) + 1
C
      CALL DGEMV('N',NT1AM(ISYMAI),NMATIJ(ISYMKJ),ONE,TMAT(KOFF1),
     *           NTOTAI,WORK,1,ONE,OMEGA1,1)
C
C--------------------------------------------
C     Third contribution to omega1
C--------------------------------------------
C
      ISYMBC = MULD2H(ISYMB,ISYMC)
      JSAKIJ = MULD2H(ISYMBC,ISYMIM)
C
      LENGTH = NCKIJ(JSAKIJ)
C
      IF (LWORK .LT. LENGTH) THEN
         CALL QUIT('Not enough core in CCFOP_ONED')
      END IF
C
      DO I = 1,LENGTH
         TMAT(I) =   SMAT(I)
     *             - SMAT(INDSQ(I,1))
      ENDDO
C
C     Symmetry sorting if symmetry
C     ----------------------------
C
      IF (NSYM .GT. 1) THEN
         CALL CC_GATHER(LENGTH,WORK,TMAT,INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK,1,TMAT,1)
      ENDIF
C
      ISYMKJ = MULD2H(ISYMBC,ISYMT2)
      ISYMAI = ISYRES
C
C        Construct M(k,j) = T(ck,bj)
C        ---------------------------
C
      DO ISYMJ = 1,NSYM
C
         ISYMK  = MULD2H(ISYMJ,ISYMKJ)
         ISYMCK = MULD2H(ISYMC,ISYMK)
         ISYCKJ = MULD2H(ISYMCK,ISYMJ)
C
         DO J = 1,NRHF(ISYMJ)
C
            NKJ = IMATIJ(ISYMK,ISYMJ)
     *          + NRHF(ISYMK)*(J - 1) 
     *          + 1
C
            NCKBJ = IT2SP(ISYCKJ,ISYMB)
     *            + NCKI(ISYCKJ)*(B-1)
     *            + ICKI(ISYMCK,ISYMJ)
     *            + NT1AM(ISYMCK)*(J-1)
     *            + IT1AM(ISYMC,ISYMK)
     *            + C
C
            CALL DCOPY(NRHF(ISYMK),T2AM(NCKBJ),NVIR(ISYMC),WORK(NKJ),1)
C
         ENDDO
      ENDDO
C
      NTOTAI = MAX(NT1AM(ISYMAI),1)
C
      KOFF1 = ISAIKL(ISYMAI,ISYMKJ) + 1
C
      CALL DGEMV('N',NT1AM(ISYMAI),NMATIJ(ISYMKJ),ONE,TMAT(KOFF1),
     *           NTOTAI,WORK,1,ONE,OMEGA1,1)
C
C----------------------------------------
C     End.
C----------------------------------------
C
      CALL QEXIT('CCFOP_ONED')
C
C
      RETURN
      END
      SUBROUTINE CCPT_TEST(AOINT,ISYAO,CMO,MOINT,LOCAMO,ABINT,IJINT,
     *                     WORK,LWORK)
C
C     Written by K. Hald, Fall 2001
C
C     Purpose : Calculate ia block of MO integrals from corresponding AO,
C               and transpose the matrix to (ai) from (ia)
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER LWORK, ISYAO
      INTEGER AL, BE
C
#if defined (SYS_CRAY)
      REAL AOINT(NBAST,NBAST), CMO(NBAST,NORBT)
      REAL MOINT(NVIRT,NRHFT)
      REAL LOCAMO(NRHFT,NVIRT)
      REAL ABINT(NVIRT,NVIRT)
      REAL IJINT(NRHFT,NRHFT)
      REAL WORK(LWORK), ZERO, ONE, RHO1N, DDOT
#else
      DOUBLE PRECISION AOINT(NBAST,NBAST), CMO(NBAST,NORBT)
      DOUBLE PRECISION MOINT(NVIRT,NRHFT)
      DOUBLE PRECISION LOCAMO(NRHFT,NVIRT)
      DOUBLE PRECISION ABINT(NVIRT,NVIRT)
      DOUBLE PRECISION IJINT(NRHFT,NRHFT)
      DOUBLE PRECISION WORK(LWORK), ZERO, ONE, RHO1N, DDOT
#endif
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C
      CALL QENTER('CCPT_TEST')
C
      CALL DZERO(MOINT,NVIRT*NRHFT)
      CALL DZERO(LOCAMO,NVIRT*NRHFT)
      CALL DZERO(ABINT,NVIRT*NVIRT)
      CALL DZERO(IJINT,NRHFT*NRHFT)
C
C-------------------------------------------
C     Print the integrals in the AO basis
C-------------------------------------------
C
      DO AL = 1, NBAST
         DO BE = 1, NBAST
           if (abs(aoint(al,be)) .gt. 1.0D-9) then
               write(lupri,*) 'X^[ao}_{ia}(',al,',',be,') = ',
     *             aoint(al,be)
           endif
         enddo
      enddo
C
C---------------------------------------------------------
C     Transform the AO integrals to the MO basis (ia)
C     The transformed vector is put in (ai) though
C     it is the (ia) block we calculate.
C---------------------------------------------------------
C
      DO AL = 1, NBAST
         DO BE = 1, NBAST
            DO A = 1, NVIRT
               DO I = 1, NRHFT
                 locamo(i,a) = locamo(i,a) 
     *                      + aoint(al,be)*CMO(al,i)*CMO(be,nrhft+a)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C-------------------------------
C     Print integrals.
C-------------------------------
C
      DO AL = 1, NRHFT
         DO BE = 1, NVIRT
           if (abs(locamo(al,be)) .gt. 1.0D-9) then
               write(lupri,*) 'X^{mo}_{ia}(',al,',',be,') = ',
     *             locamo(al,be)
           endif
         enddo
      enddo
C
C----------------------------
C     Transform integrals
C----------------------------
C
      DO AL = 1, NRHFT
         DO BE = 1, NVIRT
              moint(be,al) = locamo(al,be)
         enddo
      enddo
C
C----------------------------------------------------
C     Transform the AO integrals to MO (ij block)
C----------------------------------------------------
C
      DO AL = 1, NBAST
         DO BE = 1, NBAST
            DO I = 1, NRHFT
               DO J = 1, NRHFT
                 ijint(i,j) = ijint(i,j)
     *                      + aoint(al,be)*CMO(al,i)*CMO(be,j)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C-------------------------------
C     Print integrals.
C-------------------------------
C
      DO AL = 1, NRHFT
         DO BE = 1, NRHFT
           if (abs(ijint(al,be)) .gt. 1.0D-9) then
               write(lupri,*) 'X^{mo}_{ij}(',al,',',be,') = ',
     *             ijint(al,be)
           endif
         enddo
      enddo
C
C----------------------------------------------------
C     Transform the AO integrals to MO (ab part)
C----------------------------------------------------
C
      DO AL = 1, NBAST
         DO BE = 1, NBAST
            DO A = 1, NVIRT
               DO B = 1, NVIRT
                  abint(a,b) = abint(a,b)
     *                   + aoint(al,be)*CMO(al,nrhft+a)*CMO(be,nrhft+b)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C-------------------------------
C     Print integrals.
C-------------------------------
C
      DO AL = 1, NVIRT
         DO BE = 1, NVIRT
           if (abs(abint(al,be)) .gt. 1.0D-9) then
               write(lupri,*) 'X^{mo}_{ab}(',al,',',be,') = ',
     *             abint(al,be)
           endif
         enddo
      enddo
C
C---------------
C     END
C---------------
C
      CALL QEXIT('CCPT_TEST')
C
      RETURN
      END
C  /* DECK SUM_PT3 */
      SUBROUTINE SUM_PT3(SMAT,ISYMB,B,ISYMD,D,IAIKJ,T3SUM,IOPT)
C
C     Sum up the T3 amplitudes from the S-MAT alone (IOPT = 1).
C     Sum up the T3 amplitudes from the Q-MAT alone (IOPT = 2).
C     Sum up the T3 amplitudes from the U-MAT alone (IOPT = 3).
C     Sum up the T3 amplitudes from the W-MAT alone (IOPT = 4).
C     Sum up the T3 amplitudes only with aibjck_perm (IOPT = 5).
C     Sum up the T3 amplitudes only with aibjck_perm + bjaick_perm (IOPT = 6).
C
C     HOWEVER : PLEASE DO NOT REMOVE THIS ROUTINE
C
C     K. Hald, Fall 2001.
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYMB, ISYMD, IAIKJ, IOPT
      INTEGER KOFF1, KOFF2, ISYMJ, ISYMK, ISYMA, ISYMI, ISYAIK, ISYMBJ
      INTEGER ISYMAI, ISYBJI, KOFF3, KOFF4, KOFF5, KOFF6, KOFF7, ISYAIJ
      INTEGER ISYMDJ, ISYMDK, ISYMBK, KH
C
#if defined (SYS_CRAY)
      REAL SMAT(*), FACT
      REAL t3sum(nvirt,nvirt,nvirt,nrhft,nrhft,nrhft), HALF
#else
      DOUBLE PRECISION SMAT(*), FACT
      DOUBLE PRECISION t3sum(nvirt,nvirt,nvirt,nrhft,nrhft,nrhft), HALF
#endif
C
      LOGICAL LDEBUG
C
      PARAMETER (HALF = 0.5D0)
      PARAMETER (LDEBUG = .FALSE.)
C
      CALL QENTER('SUM_PT3')
C
C
C
      IF (IOPT .EQ. 1) THEN
         FACT = HALF
      ELSE IF ((IOPT .EQ. 2) .OR. (IOPT .EQ. 3)) THEN
         FACT = 1.0D0
      ELSE IF (IOPT .EQ. 4) THEN
         FACT = HALF
      ELSE IF (IOPT .EQ. 5) THEN
         FACT = 1.0D0
      ELSE IF (IOPT .EQ. 6) THEN
         FACT = 1.0D0
      ELSE
         CALL QUIT('Wrong IOPT in sum_pt3')
      ENDIF
C
      KOFF2 = 0
      DO KH = 1, ISYMB-1
        KOFF2 = KOFF2 + NVIR(KH)
      ENDDO
      KOFF3 = 0
      DO KH = 1, ISYMD-1
        KOFF3 = KOFF3 + NVIR(KH)
      ENDDO
C
      DO ISYMJ = 1, NSYM
         KOFF5 = 0
         DO KH = 1, ISYMJ-1
           KOFF5 = KOFF5 + NRHF(KH)
         ENDDO
         ISYAIK = MULD2H(ISYMJ,IAIKJ)
         ISYMBJ = MULD2H(ISYMB,ISYMJ)
         ISYMDJ = MULD2H(ISYMD,ISYMJ)
         DO ISYMK = 1, NSYM
            KOFF6 = 0
            DO KH = 1, ISYMK-1
              KOFF6 = KOFF6 + NRHF(KH)
            ENDDO
            ISYMAI = MULD2H(ISYAIK,ISYMK)
            ISYAIJ = MULD2H(ISYMAI,ISYMJ)
            ISYMBK = MULD2H(ISYMB,ISYMK)
            ISYMDK = MULD2H(ISYMD,ISYMK)
            DO ISYMA = 1, NSYM
C
               KOFF1 = 0
               DO KH = 1, ISYMA-1
                  KOFF1 = KOFF1 + NVIR(KH)
               ENDDO
C
               ISYMI  = MULD2H(ISYMAI,ISYMA)
C
               KOFF4 = 0
               DO KH = 1, ISYMI-1
                  KOFF4 = KOFF4 + NRHF(KH)
               ENDDO
               ISYBJI = MULD2H(ISYMBJ,ISYMI)
C
               DO I = 1, NRHF(ISYMI)
               DO J = 1, NRHF(ISYMJ)
               DO K = 1, NRHF(ISYMK)
               DO A = 1, NVIR(ISYMA)
C
               IF (IOPT .EQ. 1) THEN
                 KOFF7 = ISAIKJ(ISYAIK,ISYMJ)
     *                 + NCKI(ISYAIK)*(J - 1)
     *                 + ISAIK(ISYMAI,ISYMK) + NT1AM(ISYMAI)*(K-1)
     *                 + IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
               ELSE IF (IOPT .EQ. 2) THEN
                 KOFF7 = ISAIKJ(ISYAIJ,ISYMK)
     *                 + NCKI(ISYAIJ)*(K - 1)
     *                 + ISAIK(ISYMAI,ISYMJ) + NT1AM(ISYMAI)*(J-1)
     *                 + IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
               ELSE IF (IOPT .EQ. 3) THEN
                 KOFF7 = ISAIKJ(ISYAIK,ISYMJ)
     *                 + NCKI(ISYAIK)*(J - 1)
     *                 + ISAIK(ISYMAI,ISYMK) + NT1AM(ISYMAI)*(K-1)
     *                 + IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
               ELSE IF (IOPT .EQ. 4) THEN
                 KOFF7 = ISAIKJ(ISYAIK,ISYMJ)
     *                 + NCKI(ISYAIK)*(J - 1)
     *                 + ISAIK(ISYMAI,ISYMK) + NT1AM(ISYMAI)*(K-1)
     *                 + IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
               ELSE IF (IOPT .EQ. 5) THEN
                 KOFF7 = ISAIKJ(ISYAIK,ISYMJ)
     *                 + NCKI(ISYAIK)*(J - 1)
     *                 + ISAIK(ISYMAI,ISYMK) + NT1AM(ISYMAI)*(K-1)
     *                 + IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
               ELSE IF (IOPT .EQ. 6) THEN
                 KOFF7 = ISAIKJ(ISYAIK,ISYMJ)
     *                 + NCKI(ISYAIK)*(J - 1)
     *                 + ISAIK(ISYMAI,ISYMK) + NT1AM(ISYMAI)*(K-1)
     *                 + IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
               ENDIF
C
         IF (LDEBUG) THEN
            IF (ABS(SMAT(KOFF7)) .GT. 1.0d-12) THEN
              IF (IOPT .EQ. 1) THEN
                 WRITE(LUPRI,*) 'CONTRIBUTION FROM I = ',
     *                              KOFF7,' WITH SMAT = ',
     *                              SMAT(KOFF7)
              ELSE IF (IOPT .EQ. 2) THEN
                 WRITE(LUPRI,*) 'CONTRIBUTION FROM I = ',
     *                              KOFF7,' WITH QMAT = ',
     *                              SMAT(KOFF7)
              ELSE IF (IOPT .EQ. 3) THEN
                 WRITE(LUPRI,*) 'CONTRIBUTION FROM I = ',
     *                              KOFF7,' WITH UMAT = ',
     *                              SMAT(KOFF7)
              ELSE IF (IOPT .EQ. 4) THEN
                 WRITE(LUPRI,*) 'CONTRIBUTION FROM I = ',
     *                              KOFF7,' WITH WMAT = ',
     *                              SMAT(KOFF7)
              ELSE IF (IOPT .EQ. 5) THEN
                 WRITE(LUPRI,*) 'CONTRIBUTION FROM I = ',
     *                              KOFF7,' WITH WMAT = ',
     *                              SMAT(KOFF7)
              ELSE IF (IOPT .EQ. 6) THEN
                 WRITE(LUPRI,*) 'CONTRIBUTION FROM I = ',
     *                              KOFF7,' WITH WMAT = ',
     *                              SMAT(KOFF7)
              ENDIF
C
               WRITE(LUPRI,'(A,6I3)') ' A, B, D, I, J, K : ',A,B,D,I,J,K
               WRITE(LUPRI,*) 'KOFF1 = ',KOFF1,' KOFF2 = ',
     *                         KOFF2,'KOFF3 = ',KOFF3
               WRITE(LUPRI,*) 'KOFF4 = ',KOFF4,' KOFF5 = ',
     *                         KOFF5,'KOFF6 = ',KOFF6
            endif
         ENDIF
C
         T3SUM(A+KOFF1,B+KOFF2,D+KOFF3,I+KOFF4,J+KOFF5,K+KOFF6) = 
     *   T3SUM(A+KOFF1,B+KOFF2,D+KOFF3,I+KOFF4,J+KOFF5,K+KOFF6) +
     *   FACT*SMAT(KOFF7)
C
         IF (IOPT .NE. 5) THEN
C
         T3SUM(B+KOFF2,A+KOFF1,D+KOFF3,J+KOFF5,I+KOFF4,K+KOFF6) = 
     *   T3SUM(B+KOFF2,A+KOFF1,D+KOFF3,J+KOFF5,I+KOFF4,K+KOFF6) +
     *   FACT*SMAT(KOFF7)
C
         IF (IOPT .NE. 6) THEN
C
         T3SUM(A+KOFF1,D+KOFF3,B+KOFF2,I+KOFF4,K+KOFF6,J+KOFF5) = 
     *   T3SUM(A+KOFF1,D+KOFF3,B+KOFF2,I+KOFF4,K+KOFF6,J+KOFF5) +
     *   FACT*SMAT(KOFF7)
C
         T3SUM(B+KOFF2,D+KOFF3,A+KOFF1,J+KOFF5,K+KOFF6,I+KOFF4) = 
     *   T3SUM(B+KOFF2,D+KOFF3,A+KOFF1,J+KOFF5,K+KOFF6,I+KOFF4) +
     *   FACT*SMAT(KOFF7)
C
         T3SUM(D+KOFF3,A+KOFF1,B+KOFF2,K+KOFF6,I+KOFF4,J+KOFF5) = 
     *   T3SUM(D+KOFF3,A+KOFF1,B+KOFF2,K+KOFF6,I+KOFF4,J+KOFF5) +
     *   FACT*SMAT(KOFF7)
C
         T3SUM(D+KOFF3,B+KOFF2,A+KOFF1,K+KOFF6,J+KOFF5,I+KOFF4) = 
     *   T3SUM(D+KOFF3,B+KOFF2,A+KOFF1,K+KOFF6,J+KOFF5,I+KOFF4) +
     *   FACT*SMAT(KOFF7)
C
         END IF
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
      CALL QEXIT('SUM_PT3')
C
    1 FORMAT(1X,A8,I3,A1,I3,A1,I3,A1,I3,A1,I3,A1,I3,A4,E20.10)
      RETURN
      END
C  /* DECK PRINT_PT3 */
      SUBROUTINE PRINT_PT3(T3SUM,ISYMIM,IOPT)
C
C     Remove the forbidden T3 amplitudes and print.
C
C     HOWEVER : PLEASE DO NOT REMOVE THIS ROUTINE
C
C     K. Hald, Fall 2001.
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYMIM, IOPT, KOFF1, KOFF2, KOFF3, KOFF4, KOFF5, KOFF6
      INTEGER ISYMA, ISYMB, ISYMAB, IAIKJ, ISYMC, ISYMI, ISYMJ, ISYMK
      INTEGER KH, ISYABC, ISYIJK, ISYMJK, ISYMAI, ISYMBC
C
#if defined (SYS_CRAY)
      REAL t3sum(nvirt,nvirt,nvirt,nrhft,nrhft,nrhft), zero
#else
      DOUBLE PRECISION t3sum(nvirt,nvirt,nvirt,nrhft,nrhft,nrhft), zero
#endif
C
      LOGICAL LDEBUG
      PARAMETER ( LDEBUG = .FALSE. )
C
      PARAMETER (ZERO = 0.0D0)
C
      CALL QENTER('PRINT_PT3')
C
C
C--------------------------------------------------------------
C     Remove forbidden elements of the triples amplitudes.
C--------------------------------------------------------------
C
      DO ISYMA = 1, NSYM
         KOFF1 = 0
         DO KH = 1, ISYMA-1
            KOFF1 = KOFF1 + NVIR(KH)
         ENDDO
      DO ISYMI = 1, NSYM
         KOFF4 = 0
         DO KH = 1, ISYMI-1
            KOFF4 = KOFF4 + NRHF(KH)
         ENDDO
C
         ISYMAI = MULD2H(ISYMA,ISYMI)
         ISYMBC = MULD2H(ISYMAI,ISYMIM)
C
         DO ISYMB = 1, NSYM
            KOFF2 = 0
            DO KH = 1, ISYMB-1
               KOFF2 = KOFF2 + NVIR(KH)
            ENDDO
C
            ISYMC = MULD2H(ISYMBC,ISYMB)
               KOFF3 = 0
               DO KH = 1, ISYMC-1
                  KOFF3 = KOFF3 + NVIR(KH)
               ENDDO
C
            DO A = 1, NVIR(ISYMA)
            DO I = 1, NRHF(ISYMI)
C
               DO B = 1, NVIR(ISYMB)
               DO C = 1, NVIR(ISYMC)
C
         T3SUM(A+KOFF1,B+KOFF2,C+KOFF3,I+KOFF4,I+KOFF4,I+KOFF4) = ZERO
         T3SUM(A+KOFF1,C+KOFF3,B+KOFF2,I+KOFF4,I+KOFF4,I+KOFF4) = ZERO
         T3SUM(B+KOFF2,A+KOFF1,C+KOFF3,I+KOFF4,I+KOFF4,I+KOFF4) = ZERO
         T3SUM(B+KOFF2,C+KOFF3,A+KOFF1,I+KOFF4,I+KOFF4,I+KOFF4) = ZERO
         T3SUM(C+KOFF3,A+KOFF1,B+KOFF2,I+KOFF4,I+KOFF4,I+KOFF4) = ZERO
         T3SUM(C+KOFF3,B+KOFF2,A+KOFF1,I+KOFF4,I+KOFF4,I+KOFF4) = ZERO
C
               ENDDO ! C
               ENDDO ! B
            ENDDO    ! I
            ENDDO    ! A
         ENDDO       ! ISYMB
C
         DO ISYMJ = 1, NSYM
            KOFF5 = 0
            DO KH = 1, ISYMJ-1
               KOFF5 = KOFF5 + NRHF(KH)
            ENDDO
C
            ISYMK = MULD2H(ISYMBC,ISYMJ)
               KOFF6 = 0
               DO KH = 1, ISYMK-1
                  KOFF6 = KOFF6 + NRHF(KH)
               ENDDO
C
            DO A = 1,NVIR(ISYMA)
            DO I = 1,NRHF(ISYMI)
               DO J = 1,NRHF(ISYMJ)
               DO K = 1,NRHF(ISYMK)
C
         T3SUM(A+KOFF1,A+KOFF1,A+KOFF1,I+KOFF4,J+KOFF5,K+KOFF6) = ZERO
         T3SUM(A+KOFF1,A+KOFF1,A+KOFF1,I+KOFF4,K+KOFF6,J+KOFF5) = ZERO
         T3SUM(A+KOFF1,A+KOFF1,A+KOFF1,J+KOFF5,I+KOFF4,K+KOFF6) = ZERO
         T3SUM(A+KOFF1,A+KOFF1,A+KOFF1,J+KOFF5,K+KOFF6,I+KOFF4) = ZERO
         T3SUM(A+KOFF1,A+KOFF1,A+KOFF1,K+KOFF6,I+KOFF4,J+KOFF5) = ZERO
         T3SUM(A+KOFF1,A+KOFF1,A+KOFF1,K+KOFF6,J+KOFF5,I+KOFF4) = ZERO
C
               ENDDO ! K
               ENDDO ! J
            ENDDO    ! I
            ENDDO    ! A
         ENDDO
      ENDDO          ! ISYMI
      ENDDO          ! ISYMA
C
C----------------------------------------
C     Print the triples amplitudes.
C----------------------------------------
C
      DO ISYMA = 1, NSYM
C
         KOFF1 = 0
         DO KH = 1, ISYMA-1
           KOFF1 = KOFF1 + NVIR(KH)
         ENDDO
C
         DO ISYMB = 1, NSYM
C
             ISYMAB = MULD2H(ISYMA,ISYMB)
C
             KOFF2 = 0
             DO KH = 1, ISYMB-1
               KOFF2 = KOFF2 + NVIR(KH)
             ENDDO
C
            DO ISYMC = 1, NSYM
C
                ISYABC = MULD2H(ISYMAB,ISYMC)
                ISYIJK = MULD2H(ISYMIM,ISYABC)
C
                KOFF3 = 0
                DO KH = 1, ISYMC-1
                  KOFF3 = KOFF3 + NVIR(KH)
                ENDDO
C
                DO ISYMI = 1, NSYM
C
                KOFF4 = 0
                DO KH = 1, ISYMI-1
                  KOFF4 = KOFF4 + NRHF(KH)
                ENDDO
C
                ISYMJK = MULD2H(ISYIJK,ISYMI)
C
                DO ISYMJ = 1, NSYM
C
                KOFF5 = 0
                DO KH = 1, ISYMJ-1
                  KOFF5 = KOFF5 + NRHF(KH)
                ENDDO
C
                ISYMK = MULD2H(ISYMJK,ISYMJ)
C
                KOFF6 = 0
                DO KH = 1, ISYMK-1
                  KOFF6 = KOFF6 + NRHF(KH)
                ENDDO
C
                DO A = 1, NVIR(ISYMA)
                DO B = 1, NVIR(ISYMB)
                DO C = 1, NVIR(ISYMC)
                DO I = 1, NRHF(ISYMI)
                DO J = 1, NRHF(ISYMJ)
                DO K = 1, NRHF(ISYMK)
C
        IF (ABS(T3SUM(A+KOFF1,B+KOFF2,C+KOFF3,I+KOFF4,J+KOFF5,K+KOFF6)) 
     *                                    .GT. 1.0D-12) THEN
           IF (IOPT .EQ. 1) THEN
              write(lupri,1) 'T3AM(',a+koff1,',',b+koff2,',',
     *                               c+koff3,',',i+koff4,',',
     *                               j+koff5,',',k+koff6,') = ',
     *        T3SUM(A+KOFF1,B+KOFF2,C+KOFF3,I+KOFF4,J+KOFF5,K+KOFF6) 
C
           ELSE IF (IOPT .EQ. 2) THEN
              write(lupri,1) 'T3-BAR(',a+koff1,',',b+koff2,',',
     *                                 c+koff3,',',i+koff4,',',
     *                                 j+koff5,',',k+koff6,') = ',
     *        T3SUM(A+KOFF1,B+KOFF2,C+KOFF3,I+KOFF4,J+KOFF5,K+KOFF6)
C
           ELSE IF (IOPT .EQ. 3) THEN
              write(lupri,1) 'L3AM(',a+koff1,',',b+koff2,',',
     *                               c+koff3,',',i+koff4,',',
     *                               j+koff5,',',k+koff6,') = ',
     *        T3SUM(A+KOFF1,B+KOFF2,C+KOFF3,I+KOFF4,J+KOFF5,K+KOFF6) 
C
           ELSE IF (IOPT .EQ. 4) THEN
              write(lupri,1) 'WMAT(',a+koff1,',',b+koff2,',',
     *                               c+koff3,',',i+koff4,',',
     *                               j+koff5,',',k+koff6,') = ',
     *        T3SUM(A+KOFF1,B+KOFF2,C+KOFF3,I+KOFF4,J+KOFF5,K+KOFF6) 
           ELSE IF (IOPT .EQ. 5) THEN
              write(lupri,1) 'WBDD(',a+koff1,',',b+koff2,',',
     *                               c+koff3,',',i+koff4,',',
     *                               j+koff5,',',k+koff6,') = ',
     *        T3SUM(A+KOFF1,B+KOFF2,C+KOFF3,I+KOFF4,J+KOFF5,K+KOFF6) 
C
           ELSE
             CALL QUIT('Wrong IOPT in PRINT_PT3')
           ENDIF
        ENDIF
C
                ENDDO ! K
                ENDDO ! J
                ENDDO ! I
                ENDDO ! C
                ENDDO ! B
                ENDDO ! A
                ENDDO ! ISYMJ
                ENDDO ! ISYMI
                ENDDO ! ISYMC
         ENDDO        ! ISYMB
C
      ENDDO           ! ISYMA
C
      CALL QEXIT('PRINT_PT3')
C
    1 FORMAT(1X,A6,I3,A1,I3,A1,I3,A1,I3,A1,I3,A1,I3,A4,E20.10)
      RETURN
      END
C  /* deck dens1to2 */
      SUBROUTINE DENS1TO2(DENS1,DENS2,ISYRES)
C
C     Written by K. Hald, Fall 2001
C
C     Purpose : Calculate the contributions to the 2 electron
C               density which is identical to the 1 electron
C               density with a delta function.
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYRES, ISYMA, ISYMI, ISYMJ, ISYMIJ, ISYMJJ, KOFF1
      INTEGER KOFF2
C
#if defined (SYS_CRAY)
      REAL DENS1(*), DENS2(*), TWO
#else
      DOUBLE PRECISION DENS1(*), DENS2(*), TWO
#endif
C
      PARAMETER (TWO = 2.0D0)
C
      CALL QENTER('DENS1TO2')
C
C--------------------------
C     First contribution
C--------------------------
C
      DO ISYMA = 1, NSYM
         ISYMI = MULD2H(ISYMA,ISYRES)
         DO ISYMJ = 1, NSYM
            ISYMIJ = MULD2H(ISYMI,ISYMJ)
            ISYMJJ = MULD2H(ISYMJ,ISYMJ)   ! Now that is a tricky one
            DO A = 1, NVIR(ISYMA)
               DO I = 1, NRHF(ISYMI)
                  DO J = 1, NRHF(ISYMJ)
                     KOFF1 = I3OVIR(ISYMI,ISYMA)
     *                     + NMAIJK(ISYMI)*(A-1)
     *                     + IMAIJK(ISYMJJ,ISYMI)
     *                     + NMATIJ(ISYMJJ)*(I-1)
     *                     + IMATIJ(ISYMJ,ISYMJ)
     *                     + NRHF(ISYMJ)*(J-1)
     *                     + J
C
                     KOFF2 = IT1AM(ISYMA,ISYMI)
     *                     + NVIR(ISYMA)*(I-1)
     *                     + A
C
                     DENS2(KOFF1) = DENS2(KOFF1) - TWO*DENS1(KOFF2)
C
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C--------------------------
C     Second contribution
C--------------------------
C
      DO ISYMA = 1, NSYM
         ISYMI = MULD2H(ISYMA,ISYRES)
         DO ISYMJ = 1, NSYM
            ISYMIJ = MULD2H(ISYMI,ISYMJ)
            ISYMJJ = MULD2H(ISYMJ,ISYMJ)   ! Now that is a tricky one
            DO A = 1, NVIR(ISYMA)
               DO I = 1, NRHF(ISYMI)
                  DO J = 1, NRHF(ISYMJ)
                     KOFF1 = I3OVIR(ISYMI,ISYMA)
     *                     + NMAIJK(ISYMI)*(A-1)
     *                     + IMAIJK(ISYMIJ,ISYMJ)
     *                     + NMATIJ(ISYMIJ)*(J-1)
     *                     + IMATIJ(ISYMJ,ISYMI)
     *                     + NRHF(ISYMJ)*(I-1)
     *                     + J
C
                     KOFF2 = IT1AM(ISYMA,ISYMI)
     *                     + NVIR(ISYMA)*(I-1)
     *                     + A
C
                     DENS2(KOFF1) = DENS2(KOFF1) + DENS1(KOFF2)
C
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C--------------------------
C     End.
C--------------------------
C
      CALL QEXIT('DENS1TO2')
C
      RETURN
      END
C  /* Deck cc3_umat */
      SUBROUTINE CC3_UMAT(ECURR,T2TP,ISYMT2,TRVIR,TROCC,ISYINT,FOCKD,
     *                    DIAG,UMAT,TMAT,WORK,LWORK,INDSQ,LENSQ,
     *                    ISYMB,B,ISYMD,D)
C
C     Written by K. Hald, Fall 2001.
C
C     Purpose : Calculate the U-intermediate which is used in the
C               multiplication of T3 with T3-BAR.
C
C               U^{bd}(ck,i,j) = t2(bi,cl) (dj|lk) - t2(bi,dk) (dj|cd)
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccsdinp.h>
C
      INTEGER ISYMT2, ISYINT, LWORK, LENSQ, ISYMB, ISYMD
      INTEGER INDSQ(LENSQ,6)
      INTEGER ISYRES, ISYMBD, ISCIKJ, ISYCJF, ISYFIK, LENGTH
      INTEGER ISYMF, ISYMCJ, ISYMIK, ISYMC, ISYMJ, ISYMFJ, ISYCIK
      INTEGER KOFF1, KOFF2, KOFF3, NVIRF, NTOTCJ, ISYCIL, ISYLKJ
      INTEGER ISYMIL, ISYMI, ISYML, ISYMLK, ISYMK, ISYMFK, ISYCJK
      INTEGER NTOTCI, NRHFL, NB, ND, ISYMCI, ISYMCL, ISYCIJ, NTOTIK
      INTEGER ISYMJL, ISYBJL, ISCIKL, ISYMBJ, NTOCIK
C
#if defined (SYS_CRAY)
      REAL T2TP(*), TRVIR(*), TROCC(*), FOCKD(*), DIAG(*)
      REAL UMAT(*), TMAT(*), WORK(LWORK), EPSIBD, XUMAT
      DOUBLE PRECISION DDOT, ZERO, ONE,ECURR
#else
      DOUBLE PRECISION T2TP(*), TRVIR(*), TROCC(*), FOCKD(*), DIAG(*)
      DOUBLE PRECISION UMAT(*), TMAT(*), WORK(LWORK), EPSIBD, XUMAT
      DOUBLE PRECISION DDOT, ZERO, ONE,ECURR
#endif
      LOGICAL LDEBUG
C
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      PARAMETER(LDEBUG = .FALSE.)
C
      CALL QENTER('CC3_UMAT')
C
C
C------------------------
C     Symmetries :
C------------------------
C
      ISYRES = MULD2H(ISYMT2,ISYINT)
      ISYMBD = MULD2H(ISYMB,ISYMD)
      ISCIKJ = MULD2H(ISYMBD,ISYRES)
      ISYCJF = MULD2H(ISYMD,ISYINT)
      ISYFIK = MULD2H(ISYMB,ISYMT2)
C
C--------------------------
C     Virtual contribution.
C--------------------------
C
      LENGTH = NCKIJ(ISCIKJ)
C
      IF (LWORK .LT. LENGTH) THEN
         CALL QUIT('Insufficient memory in CC3_UMAT')
      ENDIF
C
C--------------------
C     Sort the T2.
C--------------------
C
      DO ISYMF = 1, NSYM
        ISYMIK = MULD2H(ISYMF,ISYFIK)
        DO ISYMI = 1, NSYM
           ISYMK = MULD2H(ISYMIK,ISYMI)
           ISYMFK = MULD2H(ISYMK,ISYMF)
           DO F = 1, NVIR(ISYMF)
              DO I = 1, NRHF(ISYMI)
                 DO K = 1, NRHF(ISYMK)
C
                    KOFF1 = IT2SP(ISYFIK,ISYMB)
     *                    + NCKI(ISYFIK)*(B-1)
     *                    + ISAIK(ISYMFK,ISYMI)
     *                    + NT1AM(ISYMFK)*(I-1)
     *                    + IT1AM(ISYMF,ISYMK)
     *                    + NVIR(ISYMF)*(K-1)
     *                    + F
                    KOFF2 = IMAIJA(ISYMIK,ISYMF)
     *                    + NMATIJ(ISYMIK)*(F-1)
     *                    + IMATIJ(ISYMI,ISYMK)
     *                    + NRHF(ISYMI)*(K-1)
     *                    + I
C
                    WORK(KOFF2) = T2TP(KOFF1)
C
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
      ENDDO
C
      DO ISYMF = 1, NSYM
         ISYMCJ = MULD2H(ISYMF,ISYCJF)
         ISYMIK = MULD2H(ISYFIK,ISYMF)
C
            KOFF1  = ICKATR(ISYMCJ,ISYMF)
     *             + 1
            KOFF2  = IMAIJA(ISYMIK,ISYMF)
     *             + 1
            KOFF3  = ISAIKL(ISYMCJ,ISYMIK)
     *             + 1
C
            NVIRF  = MAX(NVIR(ISYMF),1)
            NTOTCJ = MAX(NT1AM(ISYMCJ),1)
            NTOTIK = MAX(NMATIJ(ISYMIK),1)
C
            CALL DGEMM('N','T',NT1AM(ISYMCJ),NMATIJ(ISYMIK),
     *                 NVIR(ISYMF),ONE,TRVIR(KOFF1),NTOTCJ,
     *                 WORK(KOFF2),NTOTIK,ZERO,
     *                 TMAT(KOFF3),NTOTCJ)
C
      ENDDO
C
C--------------------------------------------
C     Sort if symmetry and add to umat.
C--------------------------------------------
C
      IF (NSYM .GT. 1) THEN
         CALL CCSDPT_SYMSORT(TMAT,ISCIKJ,WORK(1),LWORK)
      ENDIF
C
      DO I = 1, LENGTH
         UMAT(I) = TMAT(INDSQ(I,3))
      ENDDO
C
      IF ((IPRINT .GT. 55) .OR. LDEBUG) THEN
         XUMAT = DDOT(LENGTH,UMAT,1,UMAT,1)
         WRITE(LUPRI,*) 'In CC3_UMAT: 1. Norm of UMAT ',XUMAT
      ENDIF
C
C---------------------------
C     Occupied contribution.
C---------------------------
C
      ISYMJL = MULD2H(ISYMBD,ISYMT2)
      ISYBJL = MULD2H(ISYMJL,ISYMB)
      ISCIKL = ISYINT
C
C-------------------------
C     Sort T2
C-------------------------
C
      IF (LWORK .LT. NCKI(ISYBJL))
     *   CALL QUIT('Not enough memory in CC3_UMAT')
C
      DO ISYMJ = 1, NSYM
         ISYML  = MULD2H(ISYMJL,ISYMJ)
         ISYMBJ = MULD2H(ISYMJ,ISYMB)
C
         DO J = 1, NRHF(ISYMJ)
C
            KOFF1 = IT2SP(ISYBJL,ISYMD)
     *            + NCKI(ISYBJL)*(D - 1)
     *            + ISAIK(ISYMBJ,ISYML)  
     *            + IT1AM(ISYMB,ISYMJ)
     *            + NVIR(ISYMB)*(J-1)
     *            + B
            KOFF2 = IMATIJ(ISYML,ISYMJ)
     *            + NRHF(ISYML)*(J-1)
     *            + 1
C
            CALL DCOPY(NRHF(ISYML),T2TP(KOFF1),NT1AM(ISYMBJ),
     *                 WORK(KOFF2),1)
         ENDDO
      ENDDO
C
      DO ISYML = 1,NSYM
C
         ISYMJ  = MULD2H(ISYMJL,ISYML)
         ISYCIK = MULD2H(ISCIKL,ISYML)
C
         NTOCIK = MAX(NCKI(ISYCIK),1)
         NRHFL  = MAX(NRHF(ISYML),1)
C
         KOFF1 = ISAIKJ(ISYCIK,ISYML)
     *         + 1
         KOFF2 = IMATIJ(ISYML,ISYMJ)  
     *         + 1
         KOFF3 = ISAIKJ(ISYCIK,ISYMJ)
     *         + 1
C
         CALL DGEMM('N','N',NCKI(ISYCIK),NRHF(ISYMJ),NRHF(ISYML),
     *              -ONE,TROCC(KOFF1),NTOCIK,WORK(KOFF2),NRHFL,
     *              ZERO,TMAT(KOFF3),NTOCIK)
C
      ENDDO
C
      DO I = 1, LENGTH
         UMAT(I) = UMAT(I) + TMAT(INDSQ(I,1))
      ENDDO
C
      IF ((IPRINT .GT. 55) .OR. LDEBUG) THEN
         XUMAT = DDOT(LENGTH,UMAT,1,UMAT,1)
         WRITE(LUPRI,*) 'In CC3_UMAT: 2. Norm of UMAT ',XUMAT
      ENDIF
C
C-----------------------------------------
C     Divide by the Fock matrix diagonals.
C-----------------------------------------
C
      NB = IORB(ISYMB) + NRHF(ISYMB) + B
      ND = IORB(ISYMD) + NRHF(ISYMD) + D
C
      EPSIBD = FOCKD(NB) + FOCKD(ND) - ECURR
C
      DO L = 1,LENGTH
C
         UMAT(L) = UMAT(L)/(DIAG(L) + EPSIBD)
C
      ENDDO
C
      IF ((IPRINT .GT. 55) .OR. LDEBUG) THEN
         XUMAT = DDOT(LENGTH,UMAT,1,UMAT,1)
         WRITE(LUPRI,*) 'In CC3_UMAT: 3. Norm of UMAT ',XUMAT
      ENDIF
C
C--------------------------
C     End.
C--------------------------
C
      CALL QEXIT('CC3_UMAT')
C
      RETURN
      END
C  /* Deck ccfop_umat */
      SUBROUTINE CCFOP_UMAT(ECURR,T1AM,ISYMT1,T2TCME,ISYMT2,XIAJB,
     *                      ISINT1,FOCK,TRVIR,TRVIR7,TROCC,TROCC2,
     *                      ISYINT,FOCKD,DIAG,
     *                      UMAT,TMAT,WORK,LWORK,INDSQ,LENSQ,
     *                      ISYMB,B,ISYMD,D)
C
C     Written by K. Hald, Fall 2001.
C
C     Purpose : Calculate the U-intermediate which is used in the
C               multiplication of T3-BAR with T3.
C
C               U^{bd}(ck,i,j) = t2(bi,cl) (dj|lk) - t2(bi,dk) (dj|cd)
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccsdinp.h>
C
      INTEGER ISYMT1, ISYMT2, ISINT1, ISYINT, LWORK, LENSQ, ISYMB, ISYMD
      INTEGER INDSQ(LENSQ,6), INDEX
      INTEGER ISYRES, ISYMBD, ISCIKJ, ISYCJF, ISYFIK, LENGTH
      INTEGER ISYMF, ISYMCJ, ISYMIK, ISYMC, ISYMJ, ISYMFJ, ISYCIK
      INTEGER KOFF1, KOFF2, KOFF3, NVIRF, NTOTCJ, ISYCIL, ISYLKJ
      INTEGER ISYMIL, ISYMI, ISYML, ISYMLK, ISYMK, ISYMFK, ISYCJK
      INTEGER NTOTCI, NRHFL, NB, ND, ISYMCI, ISYMCL, ISYCIJ, NTOTIK
      INTEGER ISYAIK, ISYMAI, ISYMBK, NBK, NAI, NAIK, NAIBK, NDJ, NAIKJ
      INTEGER ISYBKJ, ISYMKJ, NBKDJ, NAIKJTEMP, NAISUM, ISYMA
      INTEGER ISYBIL, ISYMBI, ISYAJK, ISYMJK, NTOAJK, ISYAFK, ISYFIJ
      INTEGER ISYMAK, ISYMIJ, NTOTAK, ISYRES2, ISYMFI, NTOTIJ
C
#if defined (SYS_CRAY)
      REAL T1AM(*), T2TCME(*), XIAJB(*), TRVIR(*), TRVIR7(*)
      REAL TROCC(*), TROCC2(*), FOCK(*), FOCKD(*), DIAG(*)
      REAL UMAT(*), TMAT(*), WORK(LWORK), EPSIBD, XUMAT, ECURR
      REAL DDOT, ZERO, ONE, TWO
#else
      DOUBLE PRECISION T1AM(*), T2TCME(*), XIAJB(*), TRVIR(*), TRVIR7(*)
      DOUBLE PRECISION TROCC(*), TROCC2(*), FOCK(*), FOCKD(*), DIAG(*)
      DOUBLE PRECISION UMAT(*), TMAT(*), WORK(LWORK), EPSIBD, XUMAT
      DOUBLE PRECISION DDOT, ZERO, ONE, TWO, ECURR
#endif
C
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CCFOP_UMAT')
C
C------------------------
C     Symmetries :
C------------------------
C
      ISYRES  = MULD2H(ISYMT2,ISYINT)
      ISYRES2 = MULD2H(ISYMT1,ISINT1)
C
      IF (ISYRES .NE. ISYRES2)
     *       CALL QUIT('Symmetry mismatch in CCFOP_UMAT')
C
      ISYMBD  = MULD2H(ISYMB,ISYMD)
      ISCIKJ  = MULD2H(ISYMBD,ISYRES)
C
      LENGTH = NCKIJ(ISCIKJ)
C
C------------------------------------------
C     Contribution from the two T1 terms.
C------------------------------------------
C
      if (.true.) then
C
      ISYAIK = MULD2H(ISINT1,ISYMB)
C
C------------------------------------
C     Sort integrals for constant B
C------------------------------------
C
      IF (LWORK .LT. NCKI(ISYAIK)) THEN
         CALL QUIT('Too little workspace in CCFOP_UMAT (1)')
      ENDIF
C
      DO ISYMK = 1, NSYM
         ISYMAI = MULD2H(ISYAIK,ISYMK)
         ISYMBK = MULD2H(ISYMB,ISYMK)
         DO K = 1, NRHF(ISYMK)
            NBK = IT1AM(ISYMB,ISYMK) + NVIR(ISYMB)*(K - 1) + B
            DO NAI = 1, NT1AM(ISYMAI)
C
               NAIK  = ICKI(ISYMAI,ISYMK)+NT1AM(ISYMAI)*(K - 1)+NAI
               NAIBK = IT2AM(ISYMAI,ISYMBK) + INDEX(NAI,NBK)
C
               WORK(NAIK) = XIAJB(NAIBK)
C
            ENDDO
         ENDDO
      ENDDO
C
C----------------------------------
C     Contract integrals with T1.
C----------------------------------
C
      CALL DZERO(TMAT,LENGTH)
C
      ISYMJ = MULD2H(ISYMT1,ISYMD)
C
      DO ISYMK = 1, NSYM
         ISYMAI = MULD2H(ISYAIK,ISYMK)
         ISYMBK = MULD2H(ISYMB,ISYMK)
C
         DO K = 1, NRHF(ISYMK)
            NBK = IT1AM(ISYMB,ISYMK) + NVIR(ISYMB)*(K - 1) + B
            DO J = 1, NRHF(ISYMJ)
C
               NDJ = IT1AM(ISYMD,ISYMJ) + NVIR(ISYMD)*(J - 1) + D
C
               DO NAI = 1, NT1AM(ISYMAI)
C
                  NAIKJ = ISAIKJ(ISYAIK,ISYMJ)
     *                  + NCKI(ISYAIK)*(J - 1)
     *                  + ICKI(ISYMAI,ISYMK)
     *                  + NT1AM(ISYMAI)*(K-1) + NAI
C
                  NAIK  = ICKI(ISYMAI,ISYMK)+ NT1AM(ISYMAI)*(K - 1)+ NAI
C
                  TMAT(NAIKJ) = TWO*T1AM(NDJ)*WORK(NAIK)
C
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C--------------------------------------
C     Sum the result into UMAT.
C--------------------------------------
C
      DO I = 1, LENGTH
C         First :
          UMAT(I) = UMAT(I) + TMAT(INDSQ(I,3))
CCCCCC          UMAT(I) = UMAT(I) + TMAT(INDSQ(I,4))
C         Second :
          UMAT(I) = UMAT(I) - TMAT(INDSQ(I,2))
CCCCCC          UMAT(I) = UMAT(I) - TMAT(I)
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XUMAT = DDOT(LENGTH,UMAT,1,UMAT,1)
         WRITE(LUPRI,*) 'In CCFOP_UMAT: 1. Norm of UMAT ',XUMAT
      ENDIF
C
      endif
C
C-----------------------------------------------------------------------
C     Contribution from both Fock terms
C-----------------------------------------------------------------------
C
      if (.true.) then
C
      CALL DZERO(TMAT,LENGTH)
C
      ISYBKJ = MULD2H(ISYMT2,ISYMD)
      ISYMKJ = MULD2H(ISYBKJ,ISYMB)
      ISYMAI = ISYINT
C
      DO ISYMJ = 1, NSYM
         ISYMK  = MULD2H(ISYMKJ,ISYMJ)
         ISYMBK = MULD2H(ISYMB,ISYMK)
         ISYAIK = MULD2H(ISYMAI,ISYMK)
C
         DO ISYMI = 1, NSYM
C
         ISYMA = MULD2H(ISYMAI,ISYMI)
C
            DO J = 1, NRHF(ISYMJ)
               NDJ = IT1AM(ISYMD,ISYMJ) + NVIR(ISYMD)*(J-1) + D
C
               DO K = 1, NRHF(ISYMK)
C
                  NBKDJ = IT2SP(ISYBKJ,ISYMD)
     *                  + NCKI(ISYBKJ)*(D - 1)
     *                  + ICKI(ISYMBK,ISYMJ)
     *                  + NT1AM(ISYMBK)*(J - 1)
     *                  + IT1AM(ISYMB,ISYMK)
     *                  + NVIR(ISYMB)*(K-1)
     *                  + B
C
                  NAIKJTEMP = ISAIKJ(ISYAIK,ISYMJ)
     *                      + NCKI(ISYAIK)*(J - 1)
     *                      + ICKI(ISYMAI,ISYMK)
     *                      + NT1AM(ISYMAI)*(K-1)
     *                      + IT1AM(ISYMA,ISYMI)
C
                  DO A = 1, NVIR(ISYMA)
                  DO I = 1, NRHF(ISYMI)
C
                     NAIKJ = NAIKJTEMP
     *                     + NVIR(ISYMA)*(I-1)
     *                     + A
                     NAI   = IT1AM(ISYMA,ISYMI)
     *                     + NVIR(ISYMA)*(I-1)
     *                     + A
C
                     TMAT(NAIKJ) = TWO*T2TCME(NBKDJ)*FOCK(NAI)
C
                  ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C------------------------------------
C     Sum the result into SMAT.
C------------------------------------
C
      DO I = 1, LENGTH
         ! First term
         UMAT(I) = UMAT(I) + TMAT(INDSQ(I,3))
         ! Second term
         UMAT(I) = UMAT(I) - TMAT(INDSQ(I,4))
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XUMAT = DDOT(LENGTH,UMAT,1,UMAT,1)
         WRITE(LUPRI,*) 'In CCFOP_UMAT: 2. Norm of UMAT ',XUMAT
      ENDIF
C
      endif
C
C----------------------------------------------
C     Virtual contribution of L term.
C----------------------------------------------
C
      ISYAFK = MULD2H(ISYMD,ISYINT)
      ISYFIJ = MULD2H(ISYMB,ISYMT2)
C
      DO ISYMF = 1, NSYM
         ISYMIJ = MULD2H(ISYMF,ISYFIJ)
         DO ISYMI = 1, NSYM
            ISYMFI= MULD2H(ISYMF,ISYMI)
            ISYMJ = MULD2H(ISYMIJ,ISYMI)
C
            DO F = 1, NVIR(ISYMF)
               DO J = 1, NRHF(ISYMJ)
C
                      KOFF1 = IT2SP(ISYFIJ,ISYMB)
     *                      + NCKI(ISYFIJ)*(B-1)
     *                      + ICKI(ISYMFI,ISYMJ)
     *                      + NT1AM(ISYMFI)*(J-1)
     *                      + IT1AM(ISYMF,ISYMI)
     *                      + F
                      KOFF2 = IMAIJA(ISYMIJ,ISYMF)
     *                      + NMATIJ(ISYMIJ)*(F-1)
     *                      + IMATIJ(ISYMI,ISYMJ)
     *                      + NRHF(ISYMI)*(J-1)
     *                      + 1
C
                      CALL DCOPY(NRHF(ISYMI),T2TCME(KOFF1),NVIR(ISYMF),
     *                           WORK(KOFF2),1)
C
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C
      if (.true.) then
C
      DO ISYMF = 1,NSYM
C
         ISYMAK = MULD2H(ISYMF,ISYAFK)
         ISYMIJ = MULD2H(ISYMF,ISYFIJ)
C
         KOFF1 = ICKATR(ISYMAK,ISYMF)
     *         + 1
C
         KOFF2 = IMAIJA(ISYMIJ,ISYMF)
     *         + 1
         KOFF3 = ISAIKL(ISYMAK,ISYMIJ) + 1
C
         NTOTAK = MAX(1,NT1AM(ISYMAK))
         NTOTIJ = MAX(1,NMATIJ(ISYMIJ))
C
         CALL DGEMM('N','T',NT1AM(ISYMAK),NMATIJ(ISYMIJ),
     *              NVIR(ISYMF),TWO,TRVIR7(KOFF1),NTOTAK,
     *              WORK(KOFF2),NTOTIJ,ZERO,
     *              TMAT(KOFF3),NTOTAK)
C
      ENDDO
C
      IF (NSYM .GT. 1) THEN
C        Do not destroy the sorted T2!
         KOFF1 = NCKI(ISYFIJ) + 1
         CALL CCSDPT_SYMSORT(TMAT,ISCIKJ,WORK(KOFF1),LWORK)
      ENDIF
C
      DO I = 1,LENGTH
         UMAT(I) = UMAT(I) + TMAT(I)
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XUMAT = DDOT(LENGTH,UMAT,1,UMAT,1)
         WRITE(LUPRI,*) 'In CCFOP_UMAT: 3. Norm of UMAT ',XUMAT
      ENDIF
C
      endif
C
C----------------------------------------------
C     Virtual contribution of g term.
C----------------------------------------------
C
      ISYAFK = MULD2H(ISYMD,ISYINT)
      ISYFIJ = MULD2H(ISYMB,ISYMT2)
C
      if (.true.) then
C
      DO ISYMF = 1,NSYM
C
         ISYMAK = MULD2H(ISYMF,ISYAFK)
         ISYMIJ = MULD2H(ISYMF,ISYFIJ)
C
         KOFF1 = ICKATR(ISYMAK,ISYMF)
     *         + 1
         KOFF2 = IMAIJA(ISYMIJ,ISYMF)
     *         + 1
         KOFF3 = ISAIKL(ISYMAK,ISYMIJ) + 1
C
         NTOTAK = MAX(1,NT1AM(ISYMAK))
         NTOTIJ = MAX(1,NMATIJ(ISYMIJ))
C
         CALL DGEMM('N','T',NT1AM(ISYMAK),NMATIJ(ISYMIJ),
     *              NVIR(ISYMF),TWO,TRVIR(KOFF1),NTOTAK,
     *              WORK(KOFF2),NTOTIJ,ZERO,
     *              TMAT(KOFF3),NTOTAK)
C
      ENDDO
C
      IF (NSYM .GT. 1) THEN
         CALL CCSDPT_SYMSORT(TMAT,ISCIKJ,WORK(1),LWORK)
      ENDIF
C
      DO I = 1,LENGTH
         UMAT(I) = UMAT(I) - TMAT(INDSQ(I,5))
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XUMAT = DDOT(LENGTH,UMAT,1,UMAT,1)
         WRITE(LUPRI,*) 'In CCFOP_UMAT: 4. Norm of UMAT ',XUMAT
      ENDIF
C
      endif
C
C----------------------------------------
C     Occupied L contribution.
C----------------------------------------
C
      ISYBIL = MULD2H(ISYMD,ISYMT2)
      ISYMIL = MULD2H(ISYMB,ISYBIL)
C
C     Sort the T2 amplitudes for given B and D
C
      DO ISYML = 1, NSYM
         ISYMI  = MULD2H(ISYMIL,ISYML)
         ISYMBI = MULD2H(ISYMI,ISYMB)
C
         DO I = 1, NRHF(ISYMI)
            KOFF1 = IT2SP(ISYBIL,ISYMD)
     *            + NCKI(ISYBIL)*(D - 1)
     *            + ICKI(ISYMBI,ISYML)
     *            + IT1AM(ISYMB,ISYMI)
     *            + NVIR(ISYMB)*(I-1)
     *            + B
            KOFF2 = IMATIJ(ISYML,ISYMI)
     *            + NRHF(ISYML)*(I-1)
     *            + 1
C
            CALL DCOPY(NRHF(ISYML),T2TCME(KOFF1),NT1AM(ISYMBI),
     *                 WORK(KOFF2),1)
         ENDDO
      ENDDO
C
      if (.true.) then
C
      DO ISYML = 1,NSYM
C
         ISYMI = MULD2H(ISYML,ISYMIL)
         ISYAJK = MULD2H(ISYML,ISYINT)
C
         KOFF1 = ISAIKJ(ISYAJK,ISYML)
     *         + 1
         KOFF2 = IMATIJ(ISYML,ISYMI)
     *         + 1
         KOFF3 = ISAIKJ(ISYAJK,ISYMI)
     *         + 1
C
         NTOAJK = MAX(1,NCKI(ISYAJK))
         NRHFL  = MAX(1,NRHF(ISYML))
C
         CALL DGEMM('N','N',NCKI(ISYAJK),NRHF(ISYMI),
     *              NRHF(ISYML),-TWO,TROCC2(KOFF1),NTOAJK,
     *              WORK(KOFF2),NRHFL,ZERO,TMAT(KOFF3),
     *              NTOAJK)
C
      ENDDO
C
      DO I = 1, NCKIJ(ISCIKJ)
         UMAT(I) = UMAT(I) + TMAT(INDSQ(I,1))
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XUMAT = DDOT(LENGTH,UMAT,1,UMAT,1)
         WRITE(LUPRI,*) 'In CCFOP_UMAT: 5. Norm of UMAT ',XUMAT
      ENDIF
C
      endif
C------------------------------------------
C     Occupied g contribution ... use the
C     sorted T2 from occ. L term.
C------------------------------------------
C
      ISYBIL = MULD2H(ISYMD,ISYMT2)
      ISYMIL = MULD2H(ISYMB,ISYBIL)
C
      if (.true.) then
C
      DO ISYML = 1,NSYM
C
         ISYMI = MULD2H(ISYML,ISYMIL)
         ISYAJK = MULD2H(ISYML,ISYINT)
C
         KOFF1 = ISAIKJ(ISYAJK,ISYML)
     *         + 1
         KOFF2 = IMATIJ(ISYML,ISYMI)
     *         + 1
         KOFF3 = ISAIKJ(ISYAJK,ISYMI)
     *         + 1
C
         NTOAJK = MAX(1,NCKI(ISYAJK))
         NRHFL  = MAX(1,NRHF(ISYML))
C
         CALL DGEMM('N','N',NCKI(ISYAJK),NRHF(ISYMI),
     *              NRHF(ISYML),-TWO,TROCC(KOFF1),NTOAJK,
     *              WORK(KOFF2),NRHFL,ZERO,TMAT(KOFF3),
     *              NTOAJK)
C
      ENDDO
C
      DO I = 1, LENGTH
         UMAT(I) = UMAT(I) - TMAT(INDSQ(I,2))
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XUMAT = DDOT(LENGTH,UMAT,1,UMAT,1)
         WRITE(LUPRI,*) 'In CCFOP_UMAT: 6. Norm of UMAT ',XUMAT
      ENDIF
C
      endif
C-----------------------------------------
C     Divide by the Fock matrix diagonals.
C-----------------------------------------
C
      NB = IORB(ISYMB) + NRHF(ISYMB) + B
      ND = IORB(ISYMD) + NRHF(ISYMD) + D
C
      EPSIBD = FOCKD(NB) + FOCKD(ND) - ECURR
C
      DO L = 1,LENGTH
C
         UMAT(L) = UMAT(L)/(DIAG(L) + EPSIBD)
C
      ENDDO
C
      IF (IPRINT .GT. 55) THEN
         XUMAT = DDOT(LENGTH,UMAT,1,UMAT,1)
         WRITE(LUPRI,*) 'In CCFOP_UMAT: 7. Norm of UMAT ',XUMAT
      ENDIF
C
C--------------------------
C     End.
C--------------------------
C
      CALL QEXIT('CCFOP_UMAT')
C
      RETURN
      END
C  /* Deck ccfop_sort */
      SUBROUTINE CCFOP_SORT(INT1,INT2,ISYINT,IOPT)
C
C     Written by K. Hald, Fall 2001.
C
C     Sort INT1 and place it in INT2.
C     IOPT = 1 : Sort FROM ljka TO akjl
C     IOPT = 2 : Sort FROM akjl TO ljka
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYINT, IOPT, KOFF1, KOFF2
      INTEGER ISYMA, ISYLJK, ISYMK, ISYMLJ, ISYMAK, ISYMJ
      INTEGER ISYML, ISYAKJ
      INTEGER ISYMAL, ISYALJ
C
#if defined (SYS_CRAY)
      REAL INT1(*), INT2(*)
#else
      DOUBLE PRECISION INT1(*), INT2(*)
#endif
C
      CALL QENTER('CCFOP_SORT')
C
C--------------------------
C     Sort.
C--------------------------
C
      DO ISYMA = 1, NSYM
         ISYLJK = MULD2H(ISYINT,ISYMA)
         DO ISYMK = 1, NSYM
            ISYMLJ = MULD2H(ISYLJK,ISYMK)
            ISYMAK = MULD2H(ISYMA,ISYMK)
            DO ISYMJ = 1, NSYM
               ISYML  = MULD2H(ISYMLJ,ISYMJ)
               ISYMAL = MULD2H(ISYMA,ISYML)
               ISYAKJ = MULD2H(ISYMAK,ISYMJ)
               ISYALJ = MULD2H(ISYMAL,ISYMJ)
C
               DO A = 1, NVIR(ISYMA)
               DO K = 1, NRHF(ISYMK)
               DO J = 1, NRHF(ISYMJ)
C
                  KOFF1 = ISJIKA(ISYLJK,ISYMA)
     *                  + NMAJIK(ISYLJK)*(A-1)
     *                  + ISJIK(ISYMLJ,ISYMK)
     *                  + NMATIJ(ISYMLJ)*(K - 1)
     *                  + IMATIJ(ISYML,ISYMJ)
     *                  + NRHF(ISYML)*(J - 1) 
     *                  + 1
C
                  KOFF2 = ISAIKJ(ISYAKJ,ISYML)
     *                  + ICKI(ISYMAK,ISYMJ)
     *                  + NT1AM(ISYMAK)*(J-1)
     *                  + IT1AM(ISYMA,ISYMK)
     *                  + NVIR(ISYMA)*(K-1)
     *                  + A
C
                  IF (IOPT .EQ. 1) THEN
                     CALL DCOPY(NRHF(ISYML),INT1(KOFF1),1,
     *                          INT2(KOFF2),NCKI(ISYAKJ))
                  ELSE IF (IOPT .EQ. 2) THEN
                     CALL DCOPY(NRHF(ISYML),INT1(KOFF2),NCKI(ISYAKJ),
     *                          INT2(KOFF1),1)
                  ELSE
                    CALL QUIT('Wrong IOPT in CCFOP_SORT')
                  ENDIF
C
               ENDDO
               ENDDO
               ENDDO
C
            ENDDO
         ENDDO
      ENDDO
C
C--------------------------
C     End.
C--------------------------
C
      CALL QEXIT('CCFOP_SORT')
C
      RETURN
      END
C  /* Deck ccsdt_kappadiag */
      SUBROUTINE CCSDT_KAPPADIAG(KAPAA,KAPII,SMATBAR,SMAT,SMAT3,
     *                           UMATBAR,UMAT,UMAT3,TMAT,INDSQ,
     *                           LENSQ,ISSMAT,WORK,LWORK)
C
C     Written by K. Hald, Fall 2001.
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER LENSQ, ISSMAT, LWORK
      INTEGER INDSQ(LENSQ,6)
      INTEGER KTMAT, KEND1, LWRK1
      INTEGER LENGTH, ISYMK, ISYAJL, ISYML, ISYMAJ, ISYMA, ISYMJ, KOFF1
      INTEGER KOFF2, KOFF3, KHCOUNT
C
#if defined (SYS_CRAY)
      REAL KAPAA(*), KAPII(*), SMATBAR(*), SMAT(*)
      REAL SMAT3(*), UMATBAR(*), UMAT(*), UMAT3(*)
      REAL TMAT(*), DDOT, WORK(LWORK)
#else
      DOUBLE PRECISION KAPAA(*), KAPII(*), SMATBAR(*), SMAT(*)
      DOUBLE PRECISION SMAT3(*), UMATBAR(*), UMAT(*), UMAT3(*)
      DOUBLE PRECISION TMAT(*), DDOT, WORK(LWORK)
#endif
C
      LOGICAL LDEBUG
      PARAMETER(LDEBUG = .FALSE.)
C
      CALL QENTER('CCSDT_KAPPADIAG')
C
C
      LENGTH = NCKIJ(ISSMAT)
C
      KTMAT = 1
      KEND1 = KTMAT + LENGTH
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0 ) THEN
         CALL QUIT('Out of memory in CCSDT_KAPPADIAG')
      ENDIF
C
      DO I = 1, LENGTH
         TMAT(I) = 
     *             SMAT(I)
     *           + UMAT(I)
     *           + SMAT3(INDSQ(I,3))
     *           + UMAT3(INDSQ(I,3))
      ENDDO
C
C--------------------------------------------------------
C     Calculate kappa_{aa} and kappa_{ii} from smatbar
C--------------------------------------------------------
C
      DO I = 1, LENGTH
         WORK(KTMAT-1+I) = 
     *                   +  SMATBAR(I)
     *                   +  UMATBAR(I)
      ENDDO
C
      DO ISYMK = 1, NSYM
         ISYAJL = MULD2H(ISSMAT,ISYMK)
         DO ISYML = 1, NSYM
            ISYMAJ = MULD2H(ISYML,ISYAJL)
            DO ISYMJ = 1, NSYM
               ISYMA = MULD2H(ISYMJ,ISYMAJ)
C              Do not use ivir(isyma) since this appears to
C              be broken. IRHF seems to be fine, but ....
               KOFF2 = 0
               DO KHCOUNT = 1, ISYMA-1
                  KOFF2 = KOFF2 + NVIR(KHCOUNT)
               ENDDO
               KOFF3 = 0
               DO KHCOUNT = 1, ISYMJ-1
                  KOFF3 = KOFF3 + NRHF(KHCOUNT)
               ENDDO
               DO A = 1, NVIR(ISYMA)
                  DO K = 1, NRHF(ISYMK)
                     DO L = 1, NRHF(ISYML)
                        DO J = 1, NRHF(ISYMJ)
C
                        KOFF1 = ISAIKJ(ISYAJL,ISYMK)
     *                        + NCKI(ISYAJL)*(K-1)
     *                        + ICKI(ISYMAJ,ISYML)
     *                        + NT1AM(ISYMAJ)*(L-1)
     *                        + IT1AM(ISYMA,ISYMJ)
     *                        + NVIR(ISYMA)*(J-1)
     *                        + A
C
                        KAPAA(KOFF2+A) = KAPAA(KOFF2+A) 
     *                                 + WORK(KTMAT+KOFF1-1)*TMAT(KOFF1)
                        KAPII(KOFF3+J) = KAPII(KOFF3+J) 
     *                                 - WORK(KTMAT+KOFF1-1)*TMAT(KOFF1)
C
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      IF (LDEBUG) THEN
         DO A = 1, NVIRT
            IF (ABS(KAPAA(A)) .GT. 1.0D-20) THEN
               WRITE(LUPRI,*) 'KAPAA(',A,') = ',KAPAA(A)
            ENDIF
         ENDDO
         DO J = 1, NRHFT
            IF (ABS(KAPII(J)) .GT. 1.0D-20) THEN
              WRITE(LUPRI,*) 'KAPII(',j,') = ',KAPII(J)
            ENDIF
         ENDDO
      ENDIF
C
C--------------------------
C     End.
C--------------------------
C
      CALL QEXIT('CCSDT_KAPPADIAG')
C
      RETURN
      END
C  /* Deck ccsdpt_symsort */
      SUBROUTINE CCSDPT_SYMSORT(UMAT,ISSMAT,WORK,LWORK)
C
C     Written by K. Hald, Fall 2001.
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISSMAT, LWORK
      INTEGER ISYMK, ISYAIJ, ISYMJ, ISYMAI, ISYMJK, NJK
      INTEGER KOFF1, KOFF2
C
#if defined (SYS_CRAY)
      REAL UMAT(*), WORK(LWRK)
#else
      DOUBLE PRECISION UMAT(*), WORK(LWORK)
#endif
C
      CALL QENTER('CCSDPT_SYMSORT')
C
      IF (LWORK .LT. NCKIJ(ISSMAT)) THEN
         CALL QUIT('Exceeded workspace in CCSDPT_SYMSORT')
      ENDIF
C
      CALL DCOPY(NCKIJ(ISSMAT),UMAT,1,WORK,1)
C
C---------------------
C     Sort.
C---------------------
C
      DO ISYMK = 1, NSYM
         ISYAIJ = MULD2H(ISSMAT,ISYMK)
         DO ISYMJ = 1, NSYM
            ISYMAI = MULD2H(ISYAIJ,ISYMJ)
            ISYMJK = MULD2H(ISYMJ,ISYMK)
            DO K = 1, NRHF(ISYMK)
               DO J = 1, NRHF(ISYMJ)
C
                  NJK   = IMATIJ(ISYMJ,ISYMK)
     *                    + NRHF(ISYMJ)*(K - 1) + J
                  KOFF1 = ISAIKL(ISYMAI,ISYMJK)
     *                  + NT1AM(ISYMAI)*(NJK - 1)
     *                  + 1
                  KOFF2 = ISAIKJ(ISYAIJ,ISYMK)
     *                  + NCKI(ISYAIJ)*(K - 1)
     *                  + ISAIK(ISYMAI,ISYMJ)
     *                  + NT1AM(ISYMAI)*(J - 1)
     *                  + 1
C
                  CALL DCOPY(NT1AM(ISYMAI),WORK(KOFF1),1,
     *                       UMAT(KOFF2),1)
C
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C---------------------
C     End.
C---------------------
C
      CALL QEXIT('CCSDPT_SYMSORT')
C
      RETURN
      END
C  /* Deck ccfop_nonrel */
      SUBROUTINE CCFOP_NONREL(OMEGA1,DENSAB,DENSIJ,ISSMAT,SMAT,SMAT2,
     *                        SMATBAR,SMATBAR2,UMAT,UMAT2,
     *                        UMATBAR,UMATBAR2,TMAT,T2TP,ISYMT2,
     *                        INDSQ,LENSQ,ISYMB,B,ISYMD,D,WORK,LWORK)
C
C     Written by K. Hald, Winter 2001/2002.
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INTEGER ISSMAT, ISYMT2, ISYMB, ISYMD, LWORK
      INTEGER ISYFIK, ISYANM, ISYRES, KT2SOR, KT2SOR2, KEND1, LWRK1
      INTEGER ISYMF, ISYMIK, ISYMI, ISYMK, ISYMFK, ISYMFI
      INTEGER KOFF1, KOFF2, KOFF3, ISYMA, ISYMNM, ISYMN, ISYMM
      INTEGER ISYMAN, ISYMAM, LENGTH, NTOTFK, NTOTNM, NTOTA
      INTEGER KTMAT, ISYMJ, ISYCKL, ISYIJK, ISYMC, ISYMBD
      INTEGER LENSQ
      INTEGER INDSQ(LENSQ,6)
C
#if defined (SYS_CRAY)
      REAL OMEGA1(*), DENSAB(*), DENSIJ(*)
      REAL SMAT(*), SMAT2(*), SMATBAR(*), SMATBAR2(*)
      REAL UMAT(*), UMAT2(*), UMATBAR(*), UMATBAR2(*)
      REAL TMAT(*), T2TP(*), WORK(LWORK), ONE, ZERO, HALF
      REAL DDOT, TEMP
#else
      DOUBLE PRECISION OMEGA1(*), DENSAB(*), DENSIJ(*)
      DOUBLE PRECISION SMAT(*), SMAT2(*), SMATBAR(*), SMATBAR2(*)
      DOUBLE PRECISION UMAT(*), UMAT2(*), UMATBAR(*), UMATBAR2(*)
      DOUBLE PRECISION TMAT(*), T2TP(*), WORK(LWORK), ONE, ZERO, HALF
      DOUBLE PRECISION DDOT, TEMP
#endif
C
      PARAMETER(ONE = 1.0D0, ZERO = 0.0D0, HALF = 0.5D0)
C
      CALL QENTER('CCFOP_NONREL')
C
      LENGTH = NCKIJ(ISSMAT)
C
C------------------------------------------------------------
C     Calculate d_{ia} from <\mu_{3}| [[V,T2],T2]|HF >
C------------------------------------------------------------
C
      if (.true.) then
C
      ISYFIK = MULD2H(ISYMT2,ISYMB)
      ISYANM = MULD2H(ISYMT2,ISYMD)
      ISYMBD = MULD2H(ISYMB,ISYMD)
      ISYRES = MULD2H(ISSMAT,ISYMBD)      ! *isymt2*isymt2
C
      KT2SOR  = 1
      KT2SOR2 = KT2SOR  + NCKI(ISYFIK)
      KEND1   = KT2SOR2 + NCKI(ISYANM)
      LWRK1   = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('No more workspace in CCFOP_NONREL (T2,T2)')
      ENDIF
C
C-------------------------------
C     Sort T2 ... to two diff.
C-------------------------------
C
      DO ISYMF = 1, NSYM
         ISYMIK = MULD2H(ISYFIK,ISYMF)
         DO ISYMI = 1, NSYM
            ISYMK  = MULD2H(ISYMIK,ISYMI)
            ISYMFK = MULD2H(ISYMF,ISYMK)
            ISYMFI = MULD2H(ISYMF,ISYMI)
C
            DO K = 1, NRHF(ISYMK)
               DO I = 1, NRHF(ISYMI)
C
                  KOFF1 = IT2SP(ISYFIK,ISYMB)
     *                  + NCKI(ISYFIK)*(B-1)
     *                  + ISAIK(ISYMFI,ISYMK)
     *                  + NT1AM(ISYMFI)*(K-1)
     *                  + IT1AM(ISYMF,ISYMI)
     *                  + NVIR(ISYMF)*(I-1)
     *                  + 1
                  KOFF2 = KT2SOR
     *                  + ISAIK(ISYMFK,ISYMI)
     *                  + NT1AM(ISYMFK)*(I-1)
     *                  + IT1AM(ISYMF,ISYMK)
     *                  + NVIR(ISYMF)*(K-1)
C
                  CALL DCOPY(NVIR(ISYMF),T2TP(KOFF1),1,WORK(KOFF2),1)
C
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      DO ISYMA = 1, NSYM
         ISYMNM = MULD2H(ISYANM,ISYMA)
         DO ISYMN = 1, NSYM
            ISYMM = MULD2H(ISYMNM,ISYMN)
            ISYMAN = MULD2H(ISYMA,ISYMN)
            ISYMAM = MULD2H(ISYMA,ISYMM)
C
            DO N = 1, NRHF(ISYMN)
               DO M = 1, NRHF(ISYMM)
C
                  KOFF1 = IT2SP(ISYANM,ISYMD)
     *                  + NCKI(ISYANM)*(D-1)
     *                  + ISAIK(ISYMAN,ISYMM)
     *                  + NT1AM(ISYMAN)*(M-1)
     *                  + IT1AM(ISYMA,ISYMN)
     *                  + NVIR(ISYMA)*(N-1)
     *                  + 1
                  KOFF2 = KT2SOR2-1
     *                  + IMAIJA(ISYMNM,ISYMA)
     *                  + IMATIJ(ISYMN,ISYMM)
     *                  + NRHF(ISYMN)*(M-1)
     *                  + N
C
                  CALL DCOPY(NVIR(ISYMA),T2TP(KOFF1),1,WORK(KOFF2),
     *                       NMATIJ(ISYMNM))
C
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      DO I = 1, LENGTH
         TMAT(I) =  SMATBAR(INDSQ(I,2))
     *           +  UMATBAR(INDSQ(I,2))
     *           + SMATBAR2(INDSQ(I,1))
     *           + UMATBAR2(INDSQ(I,1))
      ENDDO
C
      IF (NSYM .GT. 1) THEN
         CALL CC_GATHER(LENGTH,WORK(KEND1),TMAT,INDSQ(1,6))
         CALL DCOPY(LENGTH,WORK(KEND1),1,TMAT,1)
      ENDIF
C
      DO ISYMA = 1, NSYM
         ISYMI = MULD2H(ISYRES,ISYMA)
         ISYMFK = MULD2H(ISYFIK,ISYMI)
         ISYMNM = MULD2H(ISSMAT,ISYMFK)
C
         KOFF1 = ISAIKL(ISYMFK,ISYMNM)
     *         + 1
         KOFF2 = KT2SOR
     *         + ISAIK(ISYMFK,ISYMI)
         KOFF3 = KEND1
     *         + IMAIJK(ISYMNM,ISYMI)
C
         NTOTFK = MAX(1,NT1AM(ISYMFK))
         NTOTNM = MAX(1,NMATIJ(ISYMNM))
C
         CALL DGEMM('T','N',NMATIJ(ISYMNM),NRHF(ISYMI),NT1AM(ISYMFK),
     *              ONE,TMAT(KOFF1),NTOTFK,WORK(KOFF2),NTOTFK,ZERO,
     *              WORK(KOFF3),NTOTNM)
C
         KOFF1 = KT2SOR2
     *         + IMAIJA(ISYMNM,ISYMA)
         KOFF2 = KEND1
     *         + IMAIJK(ISYMNM,ISYMI)
         KOFF3 = IT1AM(ISYMA,ISYMI)
     *         + 1
C
         NTOTNM = MAX(1,NMATIJ(ISYMNM))
         NTOTA  = MAX(1,NVIR(ISYMA))
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NMATIJ(ISYMNM),
     *              -HALF,WORK(KOFF1),NTOTNM,WORK(KOFF2),NTOTNM,
     *              ONE,OMEGA1(KOFF3),NTOTA)
C
      ENDDO
C
      endif ! of if (.true.) then
C-----------------------------------------------------
C     Calculate d_{ij} from <\mu_{3}|[V,T3]|HF >
C-----------------------------------------------------
C
      KTMAT = 1
      KEND1 = KTMAT + LENGTH
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory in CCFOP_NONREL (ab-part)')
      ENDIF
C
      if (.true.) then
C
      DO I = 1, LENGTH
         TMAT(I) = 
     *             SMAT(I)
     *           + UMAT(I)
     *           + SMAT2(INDSQ(I,3))
     *           + UMAT2(INDSQ(I,3))
      ENDDO
C
      DO I = 1, LENGTH
         WORK(KTMAT-1+I) = 
     *                     SMATBAR(I)
     *                   + UMATBAR(I)
     *                   + SMATBAR2(INDSQ(I,3))
     *                   + UMATBAR2(INDSQ(I,3))
      ENDDO
C
      DO ISYMJ = 1, NSYM
         ISYMI = MULD2H(ISYRES,ISYMJ)
         ISYCKL = MULD2H(ISSMAT,ISYMI)
         DO J = 1, NRHF(ISYMJ)
            KOFF1 = KTMAT
     *            + ISAIKJ(ISYCKL,ISYMJ)
     *            + NCKI(ISYCKL)*(J-1)
            DO I = 1, NRHF(ISYMI)
               KOFF2 = ISAIKJ(ISYCKL,ISYMI)
     *               + NCKI(ISYCKL)*(I-1)
     *               + 1
               KOFF3 = IMATIJ(ISYMI,ISYMJ)
     *               + NRHF(ISYMI)*(J-1)
     *               + I
C
               TEMP  = DDOT(NCKI(ISYCKL),WORK(KOFF1),1,TMAT(KOFF2),1)
C
               DENSIJ(KOFF3) = DENSIJ(KOFF3) - HALF*TEMP
            ENDDO
         ENDDO
      ENDDO
C
      endif
C-----------------------------------------------------
C     Calculate d_{ab} from <\mu_{3}|[V,T3]|HF >
C-----------------------------------------------------
C
      if (.true.) then
C
      DO I = 1, LENGTH
         TMAT(I) = 
     *             SMAT(I)
     *           + UMAT(I)
     *           + SMAT2(INDSQ(I,3))
     *           + UMAT2(INDSQ(I,3))
      ENDDO
C
      DO I = 1, LENGTH
         WORK(KTMAT-1+I) =  
     *                     SMATBAR(I)
     *                   + UMATBAR(I)
     *                   + SMATBAR2(INDSQ(I,3))
     *                   + UMATBAR2(INDSQ(I,3))
      ENDDO
C
      CALL DCOPY(LENGTH,TMAT,1,WORK(KEND1),1)
      CALL CCFOP_SORT(WORK(KEND1),TMAT,ISSMAT,2)
C
      CALL DCOPY(LENGTH,WORK(KTMAT),1,WORK(KEND1),1)
      CALL CCFOP_SORT(WORK(KEND1),WORK(KTMAT),ISSMAT,2)
C
      DO ISYMA = 1, NSYM
         ISYMC = MULD2H(ISYRES,ISYMA)
         ISYIJK = MULD2H(ISSMAT,ISYMA)
         DO A = 1, NVIR(ISYMA)
            KOFF1 = KTMAT
     *            + ISJIKA(ISYIJK,ISYMA)
     *            + NMAJIK(ISYIJK)*(A-1)
            DO C = 1, NVIR(ISYMC)
               KOFF2 = ISJIKA(ISYIJK,ISYMC)
     *               + NMAJIK(ISYIJK)*(C-1)
     *               + 1
               KOFF3 = IMATAB(ISYMA,ISYMC)
     *               + NVIR(ISYMA)*(C-1)
     *               + A
C
               TEMP = DDOT(NMAJIK(ISYIJK),WORK(KOFF1),1,TMAT(KOFF2),1)
C
               DENSAB(KOFF3) = DENSAB(KOFF3) + HALF*TEMP
C
            ENDDO
         ENDDO
      ENDDO
C
      endif
C
C---------------------
C     End.
C---------------------
C
      CALL QEXIT('CCFOP_NONREL')
C
      RETURN
      END
