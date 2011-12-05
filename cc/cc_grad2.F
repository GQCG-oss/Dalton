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
C  /* Deck cc_grad2 */
      SUBROUTINE CC_GRAD2(I1DXORD,WORK,LWORK)
C
C     Written by Asger Halkier & Christof Haettig May 1998.
C
C     Version: 1.01
C
C     Purpose: calculate expectation values of one- and two-electron
C              operators and effective Fock matrices
C              using the Coupled Cluster density matrices!
C
C              the information about which expectation values and 
C              which expectation values are to be calculated is read
C              from a common block (ccexpfck.h)
C              expectation values are returned on this common block,
C              the eff. Fock matrices are written to a direct access
C              file and only the start address are returned on the 
C              on the common block.
C
C              I1DXORD = 0   --  calculate usual effective Fock matrices
C                                and expectation values 
C           
C              I1DXORD = 1   --  calculate effective Fock matrices
C                                from one-index transformed densities
C
C     Current models: CCD, CCSD
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxash.h>
#include <maxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <aovec.h>
#include <iratdef.h>
#include <nuclei.h>
#include <symmet.h>
#include <chrnos.h>
#include <ccorb.h>
CCN #include <infind.h> not compatible with R12!
#include <ccisao.h>
#include <r12int.h>
#include <blocks.h>
#include <ccfield.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <distcl.h>
#include <cbieri.h>
#include <eribuf.h>
#include <inftap.h>
#include <cclr.h>
#include <ccroper.h>
#include <ccropr2.h>
#include <ccfop.h>
#include <ccfro.h>
#include <ccexpfck.h>
#include <cc1dxfck.h>
#include <ccr1rsp.h>
#include <chrxyz.h>
#include <ccinftap.h>
#include <cch2d.h>
C
      LOGICAL LOCDBG 
      PARAMETER( LOCDBG = .FALSE. )
      INTEGER INDEXA(MXCORB)
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (FOUR = 4.0D0)
      PARAMETER (ISYM0 = 1)
      INTEGER LWORK
      DIMENSION WORK(LWORK)
      LOGICAL LDERINT(8,3)
      LOGICAL DIRINT, SYM1ONLY, DIRSAV
      CHARACTER*8 LABEL1, LABEL1H
C
      CHARACTER LSTRLX*(3), FILFCK*(8), MODEL*(10)
      LOGICAL LDOTHISATOM, LANYTHING, LEXPEC, LFOCK, LTWO, SQRINT
C
      CALL QENTER('CC_GRAD2')
 
      ! save the present value of the DIRECT flag
      DIRSAV = DIRECT
C
      IF (LOCDBG) WRITE(LUPRI,*) 'WORK:',LWORK,WORK(1),WORK(LWORK)
C
C------------------------------------------
C     Open file for effective Fock matrices 
C------------------------------------------
C
      IF      (I1DXORD.EQ.0) THEN
        FILFCK = FILFCKEFF
      ELSE IF (I1DXORD.EQ.1) THEN
        FILFCK = FIL1DXFCK
      ELSE
        CALL QUIT('Illegal value of I1DXORD in CC_GRAD2.')
      END IF
C
      LUFCK = -1
      CALL WOPEN2(LUFCK,FILFCK,64,0)
C
C-----------------------------------------
C     Initialization of timing parameters.
C-----------------------------------------
C
      TIMTOT = ZERO
      TIMTOT = SECOND()
      TIMDEN = ZERO
      TIMRES = ZERO
      TIRDAO = ZERO
      TIMHE2 = ZERO
      TIMONE = ZERO
      TIMONE = SECOND()
C
C----------------------------------------------------
C     check if anything to do:
C----------------------------------------------------
C
      LANYTHING = .FALSE.
    
      IF (I1DXORD.EQ.0) THEN
         DO IDLST = 1, NEXPFCKLBL
            IF ( LEXPFCK(1,IDLST) .OR. LEXPFCK(2,IDLST) ) THEN
               LANYTHING = .TRUE.
            END IF
         END DO
      ELSE IF (I1DXORD.EQ.1) THEN
         IF (N1DXFLBL.GT.0) LANYTHING = .TRUE.
      END IF

      IF (.NOT. LANYTHING) GOTO 9000
C
      IF (LOCDBG) THEN
         WRITE (LUPRI,*) 'Entering new CC_GRAD2 with I1DXORD = ',I1DXORD
         WRITE (LUPRI,*) 'DIRECT = ',DIRECT
      END IF
C
C---------------------------------------------------------------
C     for I1DXORD = 1, find list indeces for relaxation vectors:
C---------------------------------------------------------------
C
      DO IDX = 1, N1DXFLBL
         IF (LST1DXFCK(IDX)(1:2).EQ.'R1') THEN
            IRELAX1DX(IDX) = IR1KAPPA(LABRLX(IDX),FRQRLX(IDX),
     &                               ISYMRLX(IDX))
         ELSE
            CALL QUIT('Unknown List in CC_GRAD2.')
         END IF
      END DO
C
C----------------------------------------------------
C     Both zeta- and t-vectors are totally symmetric.
C----------------------------------------------------
C
      ISYMTR = 1
      ISYMOP = 1
C
C-----------------------------------
C     Initial work space allocation.
C-----------------------------------
C
      N2BSTM  = 0
      NALLAIM = 0
      NGLMDTM = 0
      DO ISYM = 1, NSYM
        N2BSTM  = MAX(N2BSTM,N2BST(ISYM))
        NALLAIM = MAX(NALLAIM,NALLAI(ISYM))
        NGLMDTM = MAX(NGLMDTM,NGLMDT(ISYM))
      END DO

      KAODEN   = 1
      KZ2AM    = KAODEN   + N2BSTM
      KT2AM    = KZ2AM    + NT2SQ(1)
      KT2AMT   = KT2AM    + NT2AMX
      KLAMDP   = KT2AMT   + NT2AMX
      KLAMDH   = KLAMDP   + NLAMDT
      KT1AM    = KLAMDH   + NLAMDT
      KEND1    = KT1AM    + NT1AMX
C
      KDNS1D   = KEND1                ! derivative inactive Fock matrix
      KRMAT    = KDNS1D   + N2BSTM    ! connection matrix
      KKAPPA   = KRMAT    + N2BSTM    ! orbital relaxation vector
      KCMOPQ   = KKAPPA   + 2*NALLAIM ! derivative MO vector
      KCMOHQ   = KCMOPQ   + NGLMDTM   ! alpha-idx transf. HF den
      KB1DHFAO = KCMOHQ   + NGLMDTM   ! alpha-idx transf. HF den
      KB2DHFAO = KB1DHFAO + N2BSTM    ! beta-idx transf. HF den
      KLAMDPQ  = KB2DHFAO + N2BSTM    ! derivative Lambda^p
      KLAMDHQ  = KLAMDPQ  + NGLMDTM   ! derivative Lambda^h
      KEND1    = KLAMDHQ  + NGLMDTM  
C
      IF (I1DXORD.GT.0) THEN
         KOVERLP  = KEND1               ! overlap matrix
         KQAOS    = KOVERLP  + N2BST(1) ! Q^ao x S^AO
         KB1KABAO = KQAOS    + N2BSTM   ! symmetrized 1-idx. trnsf. den.
         KB2KABAO = KB1KABAO + N2BSTM   ! beta-idx transf. relax. con.
         KLAMDPQ2 = KB2KABAO + N2BSTM   
         KLAMDHQ2 = KLAMDPQ2 + NGLMDTM  ! derivative lambda matrices
         KEND1    = KLAMDHQ2 + NGLMDTM  ! with left & right Q exchanged
      END IF
C
      KZ1AM  = KEND1
      KEND1  = KZ1AM  + NT1AMX
C
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for initial allocation '//
     &             'in CC_GRAD2')
      ENDIF
C
      IF (LOCDBG) THEN
         CALL AROUND("OUTPUT FROM CC_GRAD2")
         WRITE (LUPRI,*) 'CCS :',CCS
         WRITE (LUPRI,*) 'CC2 :',CC2
         WRITE (LUPRI,*) 'CCSD:',CCSD
      END IF
C
C------------------------------------------------------------------
C     if we are going to calculate one-index transformed effective
C     Fock matrices, we read here the overlap matrix:
C------------------------------------------------------------------
C
      IF (I1DXORD.EQ.1) THEN
         IF (LWRK1.LT.NBAST) 
     &        CALL QUIT('Insufficient work space in CC_GRAD2.')
         CALL RDONEL('OVERLAP ',.TRUE.,WORK(KEND1),NBAST)
         CALL CCSD_SYMSQ(WORK(KEND1),ISYM0,WORK(KOVERLP))
      END IF 
C
C----------------------------------------
C     Read zero'th order zeta amplitudes.
C----------------------------------------
C
      IF (.NOT.CCS) THEN
        IOPT   = 3
        CALL CC_RDRSP('L0',0,1,IOPT,MODEL,WORK(KZ1AM),WORK(KZ2AM))
      ELSE
        CALL DZERO(WORK(KZ1AM),NT1AMX)
        CALL DZERO(WORK(KZ2AM),NT2AMX)
      END IF
C
C--------------------------------
C     Square up zeta2 amplitudes.
C--------------------------------
C
      CALL DCOPY(NT2AMX,WORK(KZ2AM),1,WORK(KT2AM),1)
      CALL CC_T2SQ(WORK(KT2AM),WORK(KZ2AM),1)
C
C-------------------------------------------
C     Read zero'th order cluster amplitudes.
C-------------------------------------------
C
      IF (.NOT.CCS) THEN
         IOPT = 3
         CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KT1AM),WORK(KT2AM))
      ELSE
        CALL DZERO(WORK(KT1AM),NT1AMX)
        CALL DZERO(WORK(KT2AM),NT2AMX)
      END IF
C
C------------------------------------------------
C     Zero out single vectors in CCD-calculation.
C------------------------------------------------
C
      IF (CCD) THEN
         CALL DZERO(WORK(KT1AM),NT1AMX)
         CALL DZERO(WORK(KZ1AM),NT1AMX)
      ENDIF
C
C----------------------------------
C     Calculate the lamda matrices.
C----------------------------------
C
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),WORK(KEND1),
     *            LWRK1)
C
C---------------------------------------
C     Set up 2C-E of cluster amplitudes.
C---------------------------------------
C
      ISYOPE = 1
C
      CALL DCOPY(NT2AMX,WORK(KT2AM),1,WORK(KT2AMT),1)
      IOPTTCME = 1
      CALL CCSD_TCMEPK(WORK(KT2AMT),1.0D0,ISYOPE,IOPTTCME)
C
C-------------------------------
C     Work space allocation one.
C     Note that D(ai) = ZETA(ai)
C     and both D(ia) and h(ia) 
C     are stored transposed!
C-------------------------------
C
      LENBAR = 2*NT1AMX + NMATIJ(1) + NMATAB(1)
C
      KONEAI = KZ1AM
      KONEAB = KONEAI + NT1AMX
      KONEIJ = KONEAB + NMATAB(1)
      KONEIA = KONEIJ + NMATIJ(1)
      KXMAT  = KONEIA + NT1AMX
      KYMAT  = KXMAT  + NMATIJ(1)
      KMINT  = KYMAT  + NMATAB(1)
      KMIRES = KMINT  + N3ORHF(1)
      KD1ABT = KMIRES + N3ORHF(1)
      KD1IJT = KD1ABT + NMATAB(1)
      KKABAR = KD1IJT + NMATIJ(1)
      KDHFAO = KKABAR + LENBAR
      KKABAO = KDHFAO + N2BST(1)
      KCMO   = KKABAO + N2BST(1)
      KEND1  = KCMO   + NLAMDS
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for allocation 1 CC_GRAD2')
      ENDIF
C
C------------------------------------------------------
C     Initialize remaining one electron density arrays.
C------------------------------------------------------
C
      CALL DZERO(WORK(KONEAB),NMATAB(1))
      CALL DZERO(WORK(KONEIJ),NMATIJ(1))
      CALL DZERO(WORK(KONEIA),NT1AMX)
C
C--------------------------------------------------------
C     Calculate X-intermediate of zeta- and t-amplitudes.
C--------------------------------------------------------
C
      CALL CC_XI(WORK(KXMAT),WORK(KZ2AM),ISYMTR,WORK(KT2AM),ISYMOP,
     *             WORK(KEND1),LWRK1)
C
C--------------------------------------------------------
C     Calculate Y-intermediate of zeta- and t-amplitudes.
C--------------------------------------------------------
C
      CALL CC_YI(WORK(KYMAT),WORK(KZ2AM),ISYMTR,WORK(KT2AM),ISYMOP,
     *           WORK(KEND1),LWRK1)
C
C--------------------------------------------------------------
C     Construct three remaining blocks of one electron density.
C--------------------------------------------------------------
C
      CALL DCOPY(NMATAB(1),WORK(KYMAT),1,WORK(KONEAB),1)
      CALL CC_EITR(WORK(KONEAB),WORK(KONEIJ),WORK(KEND1),LWRK1,1)
      CALL DIJGEN(WORK(KONEIJ),WORK(KXMAT))
      CALL DIAGEN(WORK(KONEIA),WORK(KT2AMT),WORK(KONEAI))
C
C---------------------------------
C     Set up transposed densities.
C---------------------------------
C
      CALL DCOPY(NMATAB(1),WORK(KONEAB),1,WORK(KD1ABT),1)
      CALL DCOPY(NMATIJ(1),WORK(KONEIJ),1,WORK(KD1IJT),1)
      CALL CC_EITR(WORK(KD1ABT),WORK(KD1IJT),WORK(KEND1),LWRK1,1)
C
C----------------------------------------------
C     Read orbital relaxation vector from disc.
C----------------------------------------------
C
      CALL DZERO(WORK(KKABAR),LENBAR)
C
      IF (.NOT.RELORB) THEN
        WRITE (LUPRI,*) 'CC_GRAD2 needs RELORB=.TRUE.'
        CALL QUIT('CC_GRAD2 needs RELORB=.TRUE.')
      END IF
C
      IF (.NOT. (CCS.OR.CC2)) THEN
         LUBAR0 = -765
         CALL GPOPEN(LUBAR0,'CCKABAR0','OLD',' ','UNFORMATTED',IDUMMY,
     *               .FALSE.)
         REWIND(LUBAR0)
         READ(LUBAR0) (WORK(KKABAR+I-1), I = 1,LENBAR)
         CALL GPCLOSE(LUBAR0,'KEEP')
      END IF
C
C----------------------------------------------------------
C     Read MO-coefficients from interface file and reorder.
C----------------------------------------------------------
C
      CALL CC_GET_CMO(WORK(KCMO))
C
      CALL CMO_REORDER(WORK(KCMO),WORK(KEND1),LWRK1)
C
C--------------------------------------------------------------
C     Calculate ao-transformed zeta-kappa-bar-0 and HF density.
C--------------------------------------------------------------
C
      KOFDIJ = KKABAR
      KOFDAB = KOFDIJ + NMATIJ(1)
      KOFDAI = KOFDAB + NMATAB(1)
      KOFDIA = KOFDAI + NT1AMX
C
      ISDEN = 1
      CALL DZERO(WORK(KKABAO),N2BST(1))
      CALL CC_DENAO(WORK(KKABAO),ISDEN,WORK(KOFDAI),WORK(KOFDAB),
     *              WORK(KOFDIJ),WORK(KOFDIA),ISDEN,WORK(KCMO),1,
     *              WORK(KCMO),1,WORK(KEND1),LWRK1)
C
      CALL CCS_D1AO(WORK(KDHFAO),WORK(KEND1),LWRK1)
      CALL DSCAL(N2BST(1),HALF,WORK(KDHFAO),1) 
C
C---------------------------------------------------------
C     Add orbital relaxation for effective density matrix.
C---------------------------------------------------------
C
      CALL DCOPY(N2BST(1),WORK(KKABAO),1,WORK(KAODEN),1)
C
C------------------------------------------------------------
C     Backtransform the one electron density to AO-basis.
C     We thus have the entire effective one-electron density.
C------------------------------------------------------------
C
      ISDEN = 1
      CALL CC_DENAO(WORK(KAODEN),ISDEN,WORK(KONEAI),WORK(KONEAB),
     *              WORK(KONEIJ),WORK(KONEIA),ISDEN,WORK(KLAMDP),1,
     *              WORK(KLAMDH),1,WORK(KEND1),LWRK1)
C
C--------------------------------------------------------
C     Calculate M-intermediate of zeta- and t-amplitudes.
C--------------------------------------------------------
C
      CALL CC_MI(WORK(KMINT),WORK(KZ2AM),ISYMTR,WORK(KT2AM),ISYMOP,
     *           WORK(KEND1),LWRK1)
C
C-----------------------------------------------------------------------
C     loop over list of Fock matrices / expectation values and calculate
C     the 1el contributions to the exp. values and eff. Fock matrices:
C     Close file since opened in CCEFFFOCK1 with 'another' LUFCK
C-----------------------------------------------------------------------
C
      CALL DSCAL(N2BST(1),TWO,WORK(KDHFAO),1)
C
      CALL WCLOSE2(LUFCK,FILFCK,'KEEP')
C
      CALL CCEFFFOCK1(I1DXORD,WORK(KAODEN),WORK(KDHFAO),
     *                WORK(KEND1),LWRK1)
C
      CALL DSCAL(N2BST(1),HALF,WORK(KDHFAO),1)
C
C--------------------------------------------------------
C     Calculate resorted M-intermediate M(imjk)->M(mkij). 
C--------------------------------------------------------
C
      CALL CC_MIRS(WORK(KMIRES),WORK(KMINT))
C
      TIMONE = SECOND() - TIMONE
C
C----------------------------------------------------------------
C     Start the loop over atoms (integrals evaluated atom wise):
C     for each atom compute the contributions of the operators
C     which are associated with this atom:
C     (IATOM = 0 for HAM0)
C----------------------------------------------------------------
C
      IF (LUFCK .LE. 0) THEN
        CALL WOPEN2(LUFCK,FILFCK,64,0)
      ENDIF
C
      KEND1SV = KEND1
      LWRK1SV = LWRK1
C
      DO 90 IATOM = 0, NUCIND
    
         KEND1 = KEND1SV
         LWRK1 = LWRK1SV

         CALL CCEFFFOCKHLP1(IATOM,I1DXORD,LDOTHISATOM,LABEL1H)
            
         IF (.NOT.LDOTHISATOM) GOTO 90
         
         IF ( LABEL1H(1:5).EQ.'1DHAM' .OR. LABEL1H(1:5).EQ.'dh/dB') THEN
            ! for derivative integrals use the DIRGRD flag to
            ! switch between direct/non-direct mode
            DIRECT = DIRGRD 

            ! set dorps.h common blocks... 
            ! this has to be adapted for symmetry...
            SYM1ONLY = .FALSE.
            CALL CC_SETDORPS(LABEL1H,SYM1ONLY,0)

            IF (LOCDBG) THEN
               WRITE (LUPRI,'(A,I3,2X,A,2X,L3)') 
     &           'CC_GRAD2> IATOM,LABEL1H,DIRECT:',IATOM,LABEL1H,DIRECT
            END IF

            ! precalculate and sort derivative integrals for this atom
            IF (.NOT.DIRECT) THEN
              CALL CCDER1(IATOM,LABEL1H,LDERINT,WORK(KEND1),LWRK1)
            END IF
         ELSE
            DIRECT = DIRSAV
         END IF

         IF (DIRECT) THEN
            NTOSYM = 1
C           IF (HERDIR) THEN
C             CALL HERDI1(WORK(KEND),LWRK,IPRERI)
C           ELSE
              KCCFB1 = KEND1
              KINDXB = KCCFB1 + MXPRIM*MXCONT
              KEND   = KINDXB + (8*MXSHEL*MXCONT + 1)/IRAT
              LWRK   = LWORK  - KEND
              IF (LOCDBG) THEN
                WRITE(LUPRI,*) 'call ERIDI1> IATOM = ',IATOM
              END IF
              CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     *                    KODPP1,KODPP2,KRDPP1,KRDPP2,
     *                    KFREE,LFREE,KEND,WORK(KCCFB1),WORK(KINDXB),
     *                    WORK(KEND),LWRK,IPRERI)

C           END IF
            KEND1 = KFREE
            LWRK1 = LFREE
         ELSE
            NTOSYM = NSYM
         END IF
C
C--------------------------------------------------------------------
C        for I1DXORD=1 loop here over operator and find the 
C        ones which belong to this atoms:
C--------------------------------------------------------------------
C
         IF (I1DXORD.EQ.0) THEN
           NFOCK1LBL = 1
         ELSE IF (I1DXORD.EQ.1) THEN
           NFOCK1LBL = N1DXFLBL
         END IF
               
         DO IDLST1 = 1, NFOCK1LBL
      
            IF (I1DXORD.EQ.1) THEN
               IDLST   = IDLST1

               CALL CCEFFFOCKHLP2(IDLST,JATOM,IATOM,ISYMOPR,ISYRLX,
     &                            ISYFCK,ISCOOR,ICOOR,ICORSY,MXCOMP,
     &                            LEXPEC,LFOCK,LTWO,LABEL1,IOPER,
     &                            WORK(KKAPPA),WORK(KRMAT),WORK(KQAOS),
     &                            WORK(KT1AM),WORK(KOVERLP),
     &                            WORK(KLAMDPQ),WORK(KLAMDHQ),
     &                            WORK(KCMOPQ),WORK(KCMOHQ),
     &                            WORK(KB1DHFAO),WORK(KB2DHFAO),
     &                            WORK(KB1KABAO),WORK(KB2KABAO), 
     &                            WORK(KLAMDPQ2),WORK(KLAMDHQ2), 
     &                            WORK(KCMO),WORK(KKABAR),
     &                            WORK(KEND1),LWRK1)
               CALL DSCAL(N2BST(ISYRLX),HALF,WORK(KB1DHFAO),1) 
               CALL DSCAL(N2BST(ISYRLX),HALF,WORK(KB2DHFAO),1) 
           
               IF (LOCDBG) THEN
                  WRITE (LUPRI,*) 'IDXORD.EQ.1 case:'
                  WRITE (LUPRI,*) 'two-electron part for:'
                  WRITE (LUPRI,*) 'IDLST:',IDLST
                  WRITE (LUPRI,*) 'LABEL1:',LABEL1
                  WRITE (LUPRI,*) 'IOPER:',IOPER
                  WRITE (LUPRI,*) 'IATOM:',JATOM
                  WRITE (LUPRI,*) 'LTWO :',LTWO
C
                  XNORM=DNRM2(N2BST(ISYM0),WORK(KDHFAO),1)
                  WRITE (LUPRI,*) 'norm^2(dhfao):',xnorm
c                 WRITE (LUPRI,*) 'DHFAO:'
c                 CALL CC_PRONELAO(WORK(KDHFAO),ISYM0)

                  XNORM=DNRM2(N2BST(ISYRLX),WORK(KB1DHFAO),1)
                  WRITE (LUPRI,*) 'norm^2(b1dhfao):',xnorm
c                 WRITE (LUPRI,*) 'B1DHFAO:'
c                 CALL CC_PRONELAO(WORK(KB1DHFAO),ISYRLX)

                  XNORM=DNRM2(N2BST(ISYRLX),WORK(KB2DHFAO),1)
                  WRITE (LUPRI,*) 'norm^2(b2dhfao):',xnorm
c                 WRITE (LUPRI,*) 'B2DHFAO:'
c                 CALL CC_PRONELAO(WORK(KB2DHFAO),ISYRLX)

                  XNORM=DNRM2(N2BST(ISYM0),WORK(KKABAO),1)
                  WRITE (LUPRI,*) 'norm^2(kabao):',xnorm
c                 WRITE (LUPRI,*) 'KABAO:'
c                 CALL CC_PRONELAO(WORK(KKABAO),ISYM0)

                  XNORM=DNRM2(N2BST(ISYRLX),WORK(KB1KABAO),1)
                  WRITE (LUPRI,*) 'norm^2(b1kabao):',xnorm
c                 WRITE (LUPRI,*) 'B1KABAO:'
c                 CALL CC_PRONELAO(WORK(KB1KABAO),ISYRLX)

                  XNORM=DNRM2(N2BST(ISYRLX),WORK(KB2KABAO),1)
                  WRITE (LUPRI,*) 'norm^2(b2kabao):',xnorm
c                 WRITE (LUPRI,*) 'B2KABAO:'
c                 CALL CC_PRONELAO(WORK(KB2KABAO),ISYRLX)
               END IF

            END IF

            IF ( I1DXORD.EQ.0 .OR. (JATOM.EQ.IATOM .AND. LTWO) ) THEN
C
C--------------------------------------
C        Start the loop over integrals.
C--------------------------------------
C
         DO 100 ISYMD1 = 1,NTOSYM
C
            IF (DIRECT) THEN
C             IF (HERDIR) THEN
C               NTOT = MAXSHL
C             ELSE
                NTOT = MXCALL
C             ENDIF
            ELSE
              NTOT  = NBAS(ISYMD1)
            END IF
C
         DO 110 ILLL = 1,NTOT
C
            IF (DIRECT) THEN

              IF      (LABEL1H(1:5).EQ.'HAM0 ') THEN
                NGDER  = 0
                NBDER  = 0
                NFILES = 1
              ELSE IF (LABEL1H(1:5).EQ.'1DHAM') THEN
                NGDER  = 1
                NBDER  = 0
                NFILES = 1 + 3*NUCDEP
              ELSE IF (LABEL1H(1:5).EQ.'dh/dB') THEN
                NGDER  = 0
                NBDER  = 1
                NFILES = 1 + 3*2
              ELSE
                CALL QUIT('Unknown 2e- integral type in CC_GRAD2.')
              END IF

C             IF (HERDIR) THEN
C                CALL HERDI2(WORK(KEND1),LWRK1,INDEXA,ILLL,NUMDIS,
C    &                       IPRINT)
C             ELSE
                 IF (LOCDBG) THEN
                   WRITE(LUPRI,*) 'call ERIDI2> ILLL = ',ILLL,NTOT
                   WRITE(LUPRI,*)'NGDER,NBDER,NFIL:',NGDER,NBDER,NFILES
C                  IF ((NGDER+NBDER).EQ.1) THEN
C                    INTPRI = .true.
C                    WRITE(LUPRI,*) 'INTPRI switched on:',INTPRI
C                  END IF
                 END IF
                 CALL ERIDI2(ILLL,INDEXA,NUMDIS,NGDER,NBDER,
     *                       WORK(KODCL1),WORK(KODCL2),
     *                       WORK(KODBC1),WORK(KODBC2),
     *                       WORK(KRDBC1),WORK(KRDBC2),
     *                       WORK(KODPP1),WORK(KODPP2),
     *                       WORK(KRDPP1),WORK(KRDPP2),
     *                       WORK(KCCFB1),WORK(KINDXB),
     *                       WORK(KEND1), LWRK1,IPRERI)
C             END IF

              ! allocate memory for orbital/record LABEL1s
              NBUFMX = NBUFX(0)
              DO I = 1, NFILES-1
                NBUFMX = MAX(NBUFMX,NBUFX(I))
              END DO
              KRECNR  = KEND1
              KNRECS  = KRECNR + (NBUFMX*NFILES - 1)/IRAT + 1
              KEND1A  = KNRECS + (NFILES - 1)/IRAT + 1
              LWRK1A  = LWORK - KEND1A

              IF (LWRK1A .LT. 0) THEN
               CALL QUIT('Insufficient work space in CC_GRAD2 (ERIDI2)')
              END IF
 
              CALL RDERILBS(WORK(KRECNR),WORK(KNRECS),NBUFMX,NFILES)

              IF (LOCDBG) THEN
                WRITE(LUPRI,*) 'NUMDIS,NBUFMX:',NUMDIS,NBUFMX
                WRITE(LUPRI,*) 'INDEXA:',(INDEXA(I),I=1,NUMDIS)
                CALL FLSHFO(LUPRI)
              END IF
  
            ELSE
              NUMDIS = 1
              KEND1A = KEND1
              LWRK1A = LWRK1
            END IF
C
            DO IDEL2 = 1, NUMDIS
C
C              set orbital index and symmetry class for next
C              AO integral distribution:
C
               IF (DIRECT) THEN
                 IDEL  = INDEXA(IDEL2)
                  IF (NOAUXB) THEN
                     IDUM = 1
                     CALL IJKAUX(IDEL,IDUM,IDUM,IDUM)
                  END IF
                 ISYMD = ISAO(IDEL)
                 WRITE(LUPRI,*)'IDEL2,INDEXA:',IDEL2,INDEXA(IDEL2)
                 WRITE(LUPRI,*)'IDEL,ISAO:',IDEL,ISAO(IDEL)
               ELSE
                 IDEL  = IBAS(ISYMD1) + ILLL
                 ISYMD = ISYMD1
               END IF

               IF (LOCDBG) THEN
                 WRITE(LUPRI,*) 'IATOM:',IATOM
                 WRITE(LUPRI,*) 'ILLL,NUMDIS:',ILLL,NUMDIS
                 WRITE(LUPRI,*) 'INDEXA:',(INDEXA(I),I=1,NUMDIS)
                 WRITE(LUPRI,*) 'IDEL2:',IDEL2
                 WRITE(LUPRI,*) 'IDEL:',IDEL
                 WRITE(LUPRI,*) 'ISYMD:',ISYMD
               END IF
C
C----------------------------------------
C              Work space allocation two.
C----------------------------------------
C
               ISYDEN = ISYMD
C
               KD2IJG = KEND1A
               KD2AIG = KD2IJG + ND2IJG(ISYDEN)
               KD2IAG = KD2AIG + ND2AIG(ISYDEN)
               KD2ABG = KD2IAG + ND2AIG(ISYDEN)
               KEND2  = KD2ABG + ND2ABG(ISYDEN)
C
               IF (I1DXORD.EQ.1) THEN
                  ISYDNB = MULD2H(ISYDEN,ISYRLX)

                  ! for lambda^q transformed
                  KDB2IJG = KEND2
                  KDB2AIG = KDB2IJG + ND2IJG(ISYDNB)
                  KDB2IAG = KDB2AIG + ND2AIG(ISYDNB)
                  KDB2ABG = KDB2IAG + ND2AIG(ISYDNB)
                  KEND2   = KDB2ABG + ND2ABG(ISYDNB)

                  ! for lambda^q2 transformed
                  KDB22IJG = KEND2
                  KDB22AIG = KDB22IJG + ND2IJG(ISYDNB)
                  KDB22IAG = KDB22AIG + ND2AIG(ISYDNB)
                  KDB22ABG = KDB22IAG + ND2AIG(ISYDNB)
                  KEND2    = KDB22ABG + ND2ABG(ISYDNB)
               END IF
C
               LWRK2  = LWORK  - KEND2
               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:',KEND2
                  CALL QUIT('Insufficient space in CC_GRAD2 (2).')
               ENDIF
C
C-------------------------------------------------------
C              Initialize 4 two electron density arrays.
C-------------------------------------------------------
C
               CALL DZERO(WORK(KD2IJG),ND2IJG(ISYDEN))
               CALL DZERO(WORK(KD2AIG),ND2AIG(ISYDEN))
               CALL DZERO(WORK(KD2IAG),ND2AIG(ISYDEN))
               CALL DZERO(WORK(KD2ABG),ND2ABG(ISYDEN))
C
               IF (I1DXORD.EQ.1) THEN
                  CALL DZERO(WORK(KDB2IJG),ND2IJG(ISYDNB))
                  CALL DZERO(WORK(KDB2AIG),ND2AIG(ISYDNB))
                  CALL DZERO(WORK(KDB2IAG),ND2AIG(ISYDNB))
                  CALL DZERO(WORK(KDB2ABG),ND2ABG(ISYDNB))

                  CALL DZERO(WORK(KDB22IJG),ND2IJG(ISYDNB))
                  CALL DZERO(WORK(KDB22AIG),ND2AIG(ISYDNB))
                  CALL DZERO(WORK(KDB22IAG),ND2AIG(ISYDNB))
                  CALL DZERO(WORK(KDB22ABG),ND2ABG(ISYDNB))
               END IF
C
C-------------------------------------------------------------------
C              Calculate the two electron density d(pq,gamma;delta).
C-------------------------------------------------------------------
C
               AUTIME = SECOND()
C
               CALL CC_DEN2(WORK(KD2IJG),WORK(KD2AIG),WORK(KD2IAG),
     *                      WORK(KD2ABG),WORK(KZ2AM),WORK(KT2AM),
     *                      WORK(KT2AMT),WORK(KMINT),WORK(KXMAT),
     *                      WORK(KYMAT),WORK(KONEAB),WORK(KONEAI),
     *                      WORK(KONEIA),WORK(KMIRES),WORK(KLAMDH),1,
     *                      WORK(KLAMDP),1,WORK(KEND2),LWRK2,IDEL,ISYMD)
C
               IF (I1DXORD.EQ.1) THEN
                 CALL CC_DEN2(WORK(KDB2IJG),WORK(KDB2AIG),WORK(KDB2IAG),
     *                        WORK(KDB2ABG),WORK(KZ2AM),WORK(KT2AM),
     *                        WORK(KT2AMT),WORK(KMINT),WORK(KXMAT),
     *                        WORK(KYMAT),WORK(KONEAB),WORK(KONEAI),
     *                        WORK(KONEIA),WORK(KMIRES),
     *                        WORK(KLAMDHQ),ISYRLX,WORK(KLAMDP),1,
     *                        WORK(KEND2),LWRK2,IDEL,ISYMD)
C
                 CALL CC_DEN2(WORK(KDB2IJG),WORK(KDB2AIG),WORK(KDB2IAG),
     *                        WORK(KDB2ABG),WORK(KZ2AM),WORK(KT2AM),
     *                        WORK(KT2AMT),WORK(KMINT),WORK(KXMAT),
     *                        WORK(KYMAT),WORK(KONEAB),WORK(KONEAI),
     *                        WORK(KONEIA),WORK(KMIRES),
     *                        WORK(KLAMDH),1,WORK(KLAMDPQ),ISYRLX,
     *                        WORK(KEND2),LWRK2,IDEL,ISYMD)
C
                 CALL CC_DEN2(WORK(KDB22IJG),WORK(KDB22AIG),
     *                        WORK(KDB22IAG),WORK(KDB22ABG),
     *                        WORK(KZ2AM),WORK(KT2AM),
     *                        WORK(KT2AMT),WORK(KMINT),WORK(KXMAT),
     *                        WORK(KYMAT),WORK(KONEAB),WORK(KONEAI),
     *                        WORK(KONEIA),WORK(KMIRES),
     *                        WORK(KLAMDHQ2),ISYRLX,WORK(KLAMDP),1,
     *                        WORK(KEND2),LWRK2,IDEL,ISYMD)
C
                 CALL CC_DEN2(WORK(KDB22IJG),WORK(KDB22AIG),
     *                        WORK(KDB22IAG),WORK(KDB22ABG),
     *                        WORK(KZ2AM),WORK(KT2AM),
     *                        WORK(KT2AMT),WORK(KMINT),WORK(KXMAT),
     *                        WORK(KYMAT),WORK(KONEAB),WORK(KONEAI),
     *                        WORK(KONEIA),WORK(KMIRES),
     *                        WORK(KLAMDH),1,WORK(KLAMDPQ2),ISYRLX,
     *                        WORK(KEND2),LWRK2,IDEL,ISYMD)
               END IF
C
               AUTIME = SECOND() - AUTIME
               TIMDEN = TIMDEN + AUTIME
C
C----------------------------------------------------------------------
C              for I1DXORD=0 loop here over the operators and find the 
C              ones which belong to this atoms:
C----------------------------------------------------------------------
C
               IF (I1DXORD.EQ.0) THEN
                 NFOCK0LBL = NEXPFCKLBL
               ELSE 
                 NFOCK0LBL = 1
               END IF
               
               DO IDLST0 = 1, NFOCK0LBL
      
                  IF (I1DXORD.EQ.0) THEN
                     IDLST   = IDLST0

                   CALL CCEFFFOCKHLP3(IDLST,IATOM,JATOM,ISYMOPR,ISYFCK,
     &                                LABEL1,IOPER,IORDER,
     &                                MXCOMP,ISCOOR,ICORSY,ICOOR,
     &                                LEXPEC,LFOCK,LTWO,WORK(KDNS1D),
     &                                WORK(KEND2),LWRK2)
                  END IF

                  IF (JATOM.EQ.IATOM .AND. LTWO .AND.
     &                (LABEL1(1:5).EQ.'HAM0 ' .OR. DIRECT
     &                 .OR.LDERINT(ICORSY,ICOOR)) ) THEN 

C
C--------------------------------------------------------------------
C                 read the integrals and calculate the contributions:
C--------------------------------------------------------------------
C
                     ISYDIS = MULD2H(ISYMD,ICORSY)
C
C------------------------------------------------------------
C                    Start loop over second AO-index (gamma).
C                    and read AO integral distribution.
C------------------------------------------------------------
C
                     DO 130 ISYMG = 1, NSYM
                     DO 140 G = 1, NBAS(ISYMG)

                        ISYMAB = MULD2H(ISYMG,ISYDIS)
                        IGAM   = G + IBAS(ISYMG)
                        NUMD   = 1
                        NUMG   = 1
C
C-------------------------------------------------------------
C                       Set addresses for 2-electron densities
C-------------------------------------------------------------
C
                        AUTIME = SECOND()
C
                        ISYMPQ = MULD2H(ISYMG,ISYDEN)
C
                        KD2GIJ = KD2IJG + ID2IJG(ISYMPQ,ISYMG)
     *                         + NMATIJ(ISYMPQ)*(G - 1) 
                        KD2GAI = KD2AIG + ID2AIG(ISYMPQ,ISYMG)
     *                         + NT1AM(ISYMPQ)*(G - 1)
                        KD2GAB = KD2ABG + ID2ABG(ISYMPQ,ISYMG)
     *                         + NMATAB(ISYMPQ)*(G - 1)
                        KD2GIA = KD2IAG + ID2AIG(ISYMPQ,ISYMG)
     *                         + NT1AM(ISYMPQ)*(G - 1)
C
                        IF (I1DXORD.EQ.1) THEN
                           ISYPQ1  = MULD2H(ISYMPQ,ISYRLX)

                           KDB2GIJ = KDB2IJG + ID2IJG(ISYPQ1,ISYMG)
     *                             + NMATIJ(ISYPQ1)*(G - 1) 
                           KDB2GAI = KDB2AIG + ID2AIG(ISYPQ1,ISYMG)
     *                             + NT1AM(ISYPQ1)*(G - 1)
                           KDB2GAB = KDB2ABG + ID2ABG(ISYPQ1,ISYMG)
     *                             + NMATAB(ISYPQ1)*(G - 1)
                           KDB2GIA = KDB2IAG + ID2AIG(ISYPQ1,ISYMG)
     *                             + NT1AM(ISYPQ1)*(G - 1)

                           KDB22GIJ = KDB22IJG + ID2IJG(ISYPQ1,ISYMG)
     *                              + NMATIJ(ISYPQ1)*(G - 1) 
                           KDB22GAI = KDB22AIG + ID2AIG(ISYPQ1,ISYMG)
     *                              + NT1AM(ISYPQ1)*(G - 1)
                           KDB22GAB = KDB22ABG + ID2ABG(ISYPQ1,ISYMG)
     *                              + NMATAB(ISYPQ1)*(G - 1)
                           KDB22GIA = KDB22IAG + ID2AIG(ISYPQ1,ISYMG)
     *                              + NT1AM(ISYPQ1)*(G - 1)
                        END IF
C
C--------------------------------------------------------------------
C                       compute the contributions to the expectation
C                       values and effective Fock matrices:
C--------------------------------------------------------------------
C
                        D = IDEL - IBAS(ISYMD)
C
                        CALL CCEFFFOCK2(I1DXORD,LABEL1,
     &                                  ISCOOR,ICOOR,ICORSY,
     &                                  MXCOMP,ISYRLX,
     &                                  G,ISYMG,D,ISYMD,
     &                                  LUFCK,FILFCK,ISYFCK,IDLST,
     &                                  LFOCK,LEXPEC,ISYDEN,ISYMOPR,
     &                                  WORK(KLAMDP),WORK(KLAMDH),
     &                                  WORK(KQAOS),WORK(KDHFAO),
     &                                  WORK(KDNS1D),WORK(KKABAO),
     &                                  WORK(KB1KABAO),WORK(KB1DHFAO),
     &                                  WORK(KB2KABAO),WORK(KB2DHFAO),
     &                                  WORK(KD2GAI),  WORK(KD2GAB),  
     &                                  WORK(KD2GIJ),  WORK(KD2GIA),
     &                                  WORK(KDB2GAI),WORK(KDB2GAB),
     &                                  WORK(KDB2GIJ),WORK(KDB2GIA),
     &                                  WORK(KDB22GAI),WORK(KDB22GAB),
     &                                  WORK(KDB22GIJ),WORK(KDB22GIA),
     &                                  DIRECT,LDERINT,WORK(KRECNR),
     &                                  WORK(KNRECS),NBUFMX,NFILES,
     &                                  WORK(KEND2),LWRK2)
C
  140                CONTINUE
  130                CONTINUE
C
                  END IF ! ( (IATOM.EQ.0) .OR. LDERINT(ICORSY,ICOOR) ) 
C
               END DO ! IDLST
C
               AUTIME = SECOND() - AUTIME
               TIMRES = TIMRES + AUTIME
C
            END DO ! IDEL2
C
  110    END DO ! ILLL
  100    END DO ! ISYMD1
         END IF
         END DO
   90 CONTINUE
C
C----------------------------------------------------------------
C    some debug print out:
C----------------------------------------------------------------
C
      IF (LOCDBG .OR. DEBUG .OR. IPRINT.GE.2) THEN
         WRITE (LUPRI,'(//,2X,A)') 'Results of CC_GRAD2:'
         WRITE (LUPRI,'(2X,A,/)') '===================='

         WRITE (LUPRI,'(2X,A,I5)') 'one-index transformation order:',
     &            I1DXORD 

         WRITE (LUPRI,'(/2X,3A)') 'LABEL1     Sym. ',
     &    ' <1el> CC    <2el> CC    <1el> HF    <2el> HF',
     &    '    Fock idx '
         WRITE (LUPRI,'(2X,72("-"))')

         IF (I1DXORD.EQ.0) THEN
            NFOCKLBL = NEXPFCKLBL
         ELSE IF (I1DXORD.EQ.1) THEN
            NFOCKLBL = N1DXFLBL
         ELSE
            CALL QUIT('Illegal value of I1DXORD in CC_GRAD2.')
         END IF

         DO IDLST = 1, NFOCKLBL
            IF (I1DXORD.EQ.0) THEN
              LABEL1   = LBLEXPFCK(IDLST)
              ISYMOPR = ISYEXPFCK(IDLST)
              LEXPEC  = LEXPFCK(1,IDLST)
              LFOCK   = LEXPFCK(2,IDLST)
              ISYFCK  = ISYMOPR
            ELSE
              LABEL1  = LBL1DXFCK(IDLST)
              IOPER  = IROPER(LABEL1,ISYMOPR)
              LEXPEC = .FALSE.
              LFOCK  = .TRUE.
              LSTRLX = LST1DXFCK(IDLST)
              IDXRLX = IRELAX1DX(IDLST)
              ISYRLX = ILSTSYMRLX(LSTRLX,IDXRLX)
              ISYFCK = MULD2H(ISYMOPR,ISYRLX)
            END IF
            IF (LEXPEC.AND.LFOCK) THEN
             WRITE (LUPRI,'(2X,A8,2I3,4G16.8,2I5)') LABEL1, ISYMOPR, 
     &           ISYFCK,
     &           (EXPVALUE(I,IDLST),I=1,4),
     &           (IADRFCK(I,IDLST),I=1,2)
            ELSE IF (LEXPEC) THEN
             WRITE (LUPRI,'(2X,A8,2I3,4G16.8,A10)') LABEL1,ISYMOPR,
     &           ISYFCK,
     &           (EXPVALUE(I,IDLST),I=1,4), 
     &           '  ---  ---'
            ELSE IF (LFOCK)  THEN
             WRITE (LUPRI,'(2X,A8,2I3,A48,2I5)') LABEL1,ISYMOPR,ISYFCK,
     &           '  --.-----    --.-----    --.-----    --.-----  ',
     &             (IADRFCK(I,IDLST),I=1,2)
            ELSE
             WRITE (LUPRI,'(2X,A8,2I3,A48,A10)')LABEL1,ISYMOPR,ISYFCK,
     &           '  --.-----    --.-----    --.-----    --.-----  ',
     &           '  ---  ---'
            END IF
         END DO
C
         WRITE (LUPRI,'(2X,72("-"),/)')
C
         TIMTOT = SECOND() - TIMTOT
C
         WRITE(LUPRI,'((2X,A,F12.2))')
     &        'Total time used in CC_GRAD2:           ', TIMTOT,
     &        'Time used for setting up 2 e- density :', TIMDEN,
     &        'Time used for contraction with integr.:', TIMRES,
     &        'Time used for reading 2 e- AO-integr. :', TIRDAO,
     &        'Time used for calc. 2 e- AO-integrals :', TIMHE2,
     &        'Time used for 1 e- density & intermed.:', TIMONE

      ENDIF
C
C----------------------------------------------------------------
C     continue here if no expectation values or Fock matrices are
C     to be calculated (But still you have to close the file!):
C----------------------------------------------------------------
C
9000  CONTINUE
C
C----------------------------
C     close fock matrix file:
C----------------------------
C
      CALL WCLOSE2(LUFCK,FILFCK,'KEEP')
C
CCH  switch of print of integrals
C     INTPRI = .false.
CCH
C
      ! restore the DIRECT flag 
      DIRECT = DIRSAV

      CALL QEXIT('CC_GRAD2')
      RETURN
      END
