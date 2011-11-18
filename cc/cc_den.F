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
C  /* Deck cc_den */
      SUBROUTINE CC_DEN(POTNUC,AODEN,ZKDIA,WORK,LWORK,IOPT)
C
C     Written by Asger Halkier February/March 1996.
C
C     Version: 1.0
C
C     Purpose: To calculate various different matrices,
C              depending on the value of IOPT, all constructed
C              from the Coupled Cluster density matrices!
C
C     Current models: CCD, CCSD
C
C     For IOPT = 1 the right hand side of the equation for the
C     zero'th order orbital rotation lagrange multiplier response
C     (Zeta-Kappa-0) is returned in the array AODEN!
C
C     For IOPT = 2 both the one and two electron Coupled Cluster lamda
C     density is set up and used to calculate the ccsd-energy ECCSD!
C     (Used to debug the construction of the density-matrices)
C
C     For IOPT = 3 both the one and two electron effective Coupled Cluster
C     density is set up and used to calculate the ccsd-energy ECCSD!
C     (Used to debug the construction of the density-matrices)
C
C     Modified Summer 1998 to incorporate canonical surface, so the
C     diagonal blocks of Kappa-bar-0 is returned in ZKDIA for 
C     frozen core calculations.
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxash.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
#include <iratdef.h>
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (FOUR = 4.0D0)
      DIMENSION INDEXA(MXCORB)
      DIMENSION AODEN(*), ZKDIA(*), WORK(LWORK)
#include <ccorb.h>
#include <ccisao.h>
#include <r12int.h>
#include <inftap.h>
#include <blocks.h>
#include <ccfield.h>
#include <ccsdinp.h>
#include <ccinftap.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <distcl.h>
#include <cbieri.h>
#include <eribuf.h>
#include <ccfro.h>
C
      CHARACTER MODEL*10
      CHARACTER NAME1*8
      CHARACTER NAME2*8
C
      CALL QENTER('CC_DEN')
C
      IF ((FROIMP) .AND. (IOPT .EQ. 1)) THEN
C
         NAME1 = 'CCFRO1IN'
         NAME2 = 'CCFRO2IN'
C
         LUN1  = -1
         LUN2  = -1
C
         CALL WOPEN2(LUN1,NAME1,64,0)
         CALL WOPEN2(LUN2,NAME2,64,0)
C
      ENDIF
C
      IF (IOPT .LT. 2) THEN
         CALL HEADER('Constructing right-hand-side for CCSD-kappa-0',-1)
      ENDIF
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
      KZ2AM  = 1
      KT2AM  = KZ2AM  + NT2SQ(1)
      KT2AMT = KT2AM  + NT2AMX
      KLAMDP = KT2AMT + NT2AMX
      KLAMDH = KLAMDP + NLAMDT
      KT1AM  = KLAMDH + NLAMDT
      KZ1AM  = KT1AM  + NT1AMX
      KEND1  = KZ1AM  + NT1AMX
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'CC_DEN: Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for initial allocation '//
     &             'in CC_DEN')
      ENDIF
C
C----------------------------------------
C     Read zero'th order zeta amplitudes.
C----------------------------------------
C
      IOPTRW   = 3
      CALL CC_RDRSP('L0',0,1,IOPTRW,MODEL,WORK(KZ1AM),WORK(KZ2AM))
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
      IOPTRW = 3
      CALL CC_RDRSP('R0',0,1,IOPTRW,MODEL,WORK(KT1AM),WORK(KT2AM))
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
      KONEAI = KZ1AM
      KONEAB = KONEAI + NT1AMX
      KONEIJ = KONEAB + NMATAB(1)
      KONEIA = KONEIJ + NMATIJ(1)
      KXMAT  = KONEIA + NT1AMX
      KYMAT  = KXMAT  + NMATIJ(1)
      KMINT  = KYMAT  + NMATAB(1)
      KONINT = KMINT  + N3ORHF(1)
      KMIRES = KONINT + N2BST(ISYMOP)
      IF (IOPT .GT. 1) THEN
         IF (IOPT .EQ. 2) THEN
            KEND1  = KMIRES + N3ORHF(1)
            LWRK1  = LWORK  - KEND1
         ELSE IF (IOPT .EQ. 3) THEN
            LENBAR = 2*NT1AMX    + NMATIJ(1)   + NMATAB(1)
     *             + 2*NCOFRO(1) + 2*NT1FRO(1)
            KKABAR = KMIRES + N3ORHF(1)
            KDHFAO = KKABAR + LENBAR
            KKABAO = KDHFAO + N2BST(1)
            KCMO   = KKABAO + N2BST(1)
            KEND1  = KCMO   + NLAMDS
            LWRK1  = LWORK  - KEND1
         ENDIF
      ELSE IF (IOPT .LT. 2) THEN
         KINTAI = KMIRES + N3ORHF(1)
         KINTAB = KINTAI + NT1AMX
         KINTIJ = KINTAB + NMATAB(1)
         KINTIA = KINTIJ + NMATIJ(1)
         KINABT = KINTIA + NT1AMX
         KINIJT = KINABT + NMATAB(1)
         KD1ABT = KINIJT + NMATIJ(1)
         KD1IJT = KD1ABT + NMATAB(1)
         KEND1  = KD1IJT + NMATIJ(1)
         LWRK1  = LWORK  - KEND1
         IF (FROIMP) THEN
            KFROII = KEND1
            KFROIJ = KFROII + NFROFR(1)
            KFROJI = KFROIJ + NCOFRO(1)
            KFROAI = KFROJI + NCOFRO(1)
            KFROIA = KFROAI + NT1FRO(1)
            KFD1II = KFROIA + NT1FRO(1)
            KEND1  = KFD1II + NFROFR(1)
            LWRK1  = LWORK  - KEND1
         ENDIF
      ENDIF
C
      IF (FROIMP) THEN
         KCMOF = KEND1
         KEND1 = KCMOF + NLAMDS
         LWRK1 = LWORK - KEND1
      ENDIF
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'CC_DEN: Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for allocation 1 CC_DEN')
      ENDIF
C
      IF (FROIMP) THEN
C
C-------------------------------------------
C        Get the FULL MO coefficient matrix.
C-------------------------------------------
C
         CALL CMO_ALL(WORK(KCMOF),WORK(KEND1),LWRK1)
C
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
C     Read one-electron integrals.
C---------------------------------
C
      CALL CCRHS_ONEAO(WORK(KONINT),WORK(KEND1),LWRK1)
C
      IF (IOPT .LT. 2) THEN
C
C---------------------------------------------------------
C        Ove 24-20-1996
C        Read one-electron integrals into Fock-matrix for
C        finite field.
C---------------------------------------------------------
C
         DO 13 IF = 1, NFIELD
c           DTIME  = SECOND()
            FF = EFIELD(IF)
            CALL CC_ONEP(WORK(KONINT),WORK(KEND1),LWRK1,FF,1,LFIELD(IF))
c           DTIME  = DTIME - SECOND()
c           TIMFCK = TIMFCK + DTIME
 13      CONTINUE
C
C--------------------------------------------------
C     Transform one electron integrals to MO-basis.
C--------------------------------------------------
C
         ISYM = 1
         CALL CCDINTMO(WORK(KINTIJ),WORK(KINTIA),WORK(KINTAB),
     *                 WORK(KINTAI),WORK(KONINT),WORK(KLAMDP),
     *                 WORK(KLAMDH),WORK(KEND1),LWRK1,ISYM)
C
         IF (IPRINT .GT. 50) THEN
C
            CALL AROUND('One electron integrals in MO-basis')
C
            DO 10 ISYM = 1,NSYM
C
               WRITE(LUPRI,333) 'Symmetry block number:', ISYM
  333          FORMAT(/3X,A22,2X,I1)
C
               KOFFIJ = KINTIJ + IMATIJ(ISYM,ISYM)
               KOFFAI = KINTAI + IT1AM(ISYM,ISYM)
               KOFFIA = KINTIA + IT1AM(ISYM,ISYM)
               KOFFAB = KINTAB + IMATAB(ISYM,ISYM)
C
               CALL AROUND('occ.occ part')
               CALL OUTPUT(WORK(KOFFIJ),1,NRHF(ISYM),1,NRHF(ISYM),
     *                     NRHF(ISYM),NRHF(ISYM),1,LUPRI)
C
               CALL AROUND('vir.occ part')
               CALL OUTPUT(WORK(KOFFAI),1,NVIR(ISYM),1,NRHF(ISYM),
     *                     NVIR(ISYM),NRHF(ISYM),1,LUPRI)
C
               CALL AROUND('occ.vir part')
               CALL OUTPUT(WORK(KOFFIA),1,NVIR(ISYM),1,NRHF(ISYM),
     *                     NVIR(ISYM),NRHF(ISYM),1,LUPRI)
C
               CALL AROUND('vir.vir part')
               CALL OUTPUT(WORK(KOFFAB),1,NVIR(ISYM),1,NVIR(ISYM),
     *                     NVIR(ISYM),NVIR(ISYM),1,LUPRI)
C
  10        CONTINUE
C
            HIJNO = DDOT(NMATIJ(1),WORK(KINTIJ),1,WORK(KINTIJ),1)
            HAINO = DDOT(NT1AMX,WORK(KINTAI),1,WORK(KINTAI),1)
            HIANO = DDOT(NT1AMX,WORK(KINTIA),1,WORK(KINTIA),1)
            HABNO = DDOT(NMATAB(1),WORK(KINTAB),1,WORK(KINTAB),1)
C
            WRITE(LUPRI,*) ' '
            WRITE(LUPRI,*) 'Norm of one electron integrals in MO-basis:'
            WRITE(LUPRI,*) HIJNO, HAINO, HIANO, HABNO
C
         ENDIF
C
         IF (FROIMP) THEN
C
            ISYM = 1
            CALL CCIFROMO(WORK(KFROIJ),WORK(KFROJI),WORK(KFROAI),
     *                    WORK(KFROIA),WORK(KFROII),WORK(KONINT),
     *                    WORK(KLAMDP),WORK(KLAMDH),WORK(KCMOF),
     *                    WORK(KEND1),LWRK1,ISYM)
C
            CALL CCFD1II(WORK(KFD1II))
C
         ENDIF
C
C--------------------------------------------------
C        Set up transposed integrals and densities.
C--------------------------------------------------
C
         CALL DCOPY(NMATAB(1),WORK(KINTAB),1,WORK(KINABT),1)
         CALL DCOPY(NMATIJ(1),WORK(KINTIJ),1,WORK(KINIJT),1)
         CALL DCOPY(NMATAB(1),WORK(KONEAB),1,WORK(KD1ABT),1)
         CALL DCOPY(NMATIJ(1),WORK(KONEIJ),1,WORK(KD1IJT),1)
C
         CALL CC_EITR(WORK(KINABT),WORK(KINIJT),WORK(KEND1),LWRK1,1)
         CALL CC_EITR(WORK(KD1ABT),WORK(KD1IJT),WORK(KEND1),LWRK1,1)
C
C------------------------------------------------------------
C        Calculate one electron contribution to Zeta-kappa-0.
C------------------------------------------------------------
C
         CALL DZERO(AODEN,NT1AM(1))
C
         ISYM = 1
         CALL CCDENZK0(AODEN,WORK(KINTIJ),WORK(KINTAI),WORK(KINTIA),
     *                 WORK(KINTAB),WORK(KONEIJ),WORK(KONEAI),
     *                 WORK(KONEIA),WORK(KONEAB),WORK(KINIJT),
     *                 WORK(KINABT),WORK(KD1IJT),WORK(KD1ABT),
     *                 WORK(KT1AM),WORK(KEND1),LWRK1,ISYM)
C
         IF (FROIMP) THEN
C
C--------------------------------------------------------
C           Calculate one-electron contribution to right-
C           hand-side of correlated-frozen multipliers.
C--------------------------------------------------------
C
            ISYM = 1
            ICON = 1
            KOFRES = 2*NT1AMX + NMATIJ(1) + NMATAB(1) + 2*NT1FRO(1) + 1
            CALL CC_ETIJF(ZKDIA(KOFRES),WORK(KONEIJ),WORK(KONEAB),
     *                    WORK(KONEAI),WORK(KONEIA),WORK(KD1IJT),
     *                    WORK(KFD1II),SK1,SK2,SK3,SK4,WORK(KINTIJ),
     *                    WORK(KINTAI),WORK(KINTIA),WORK(KINIJT),
     *                    WORK(KINABT),WORK(KFROIJ),WORK(KFROJI),
     *                    WORK(KFROAI),WORK(KFROIA),WORK(KFROII),
     *                    WORK(KT1AM),WORK(KEND1),LWRK1,ISYM,ICON)
C
C--------------------------------------------------------
C           Calculate one-electron contribution to right-
C           hand-side of virtual-frozen multipliers.
C--------------------------------------------------------
C
            ISYM = 1
            ICON = 1
            KOFRES = 2*NT1AMX + NMATIJ(1) + NMATAB(1) + 1
            CALL CC_ETAIF(ZKDIA(KOFRES),WORK(KONEAB),WORK(KONEAI),
     *                    WORK(KONEIA),WORK(KD1IJT),WORK(KD1ABT),
     *                    WORK(KFD1II),SK1,SK2,SK3,SK4,WORK(KINTIJ),
     *                    WORK(KINTAB),WORK(KINTAI),WORK(KINTIA),
     *                    WORK(KINABT),WORK(KFROIJ),WORK(KFROJI),
     *                    WORK(KFROAI),WORK(KFROIA),WORK(KFROII),
     *                    WORK(KT1AM),WORK(KEND1),LWRK1,ISYM,ICON)
C
C--------------------------------------------------------
C           Calculate one-electron contribution to right-
C           hand-side of frozen-frozen multipliers.
C--------------------------------------------------------
C
c           ISYM = 1
c           ICON = 1
c           KOFRES = 2*NT1AMX + NMATIJ(1) + NMATAB(1) + 2*NT1FRO(1)
c    *             + 2*NCOFRO(1) + 1
c           CALL CC_ETIFJF(ZKDIA(KOFRES),WORK(KFD1II),SK1,SK2,SK3,SK4,
c    *                     WORK(KFROIJ),WORK(KFROJI),WORK(KFROAI),
c    *                     WORK(KFROIA),WORK(KFROII),ISYM,ICON)
C
C----------------------------------------------------
C           Calculate one-electron contribution to
C           right-hand-side for diagonal multipliers.
C----------------------------------------------------
C
c           ISYM = 1
c           CALL CCDIAZK0(ZKDIA,WORK(KINTIJ),WORK(KINTAI),WORK(KINTIA),
c    *                    WORK(KINTAB),WORK(KONEIJ),WORK(KONEAI),
c    *                    WORK(KONEIA),WORK(KONEAB),WORK(KINIJT),
c    *                    WORK(KINABT),WORK(KD1IJT),WORK(KD1ABT),
c    *                    WORK(KT1AM),WORK(KEND1),LWRK1,ISYM)
C
         ENDIF
      ENDIF
C
      IF (IOPT .GT. 1) THEN
C
C-------------------------------------------------
C        Read orbital relaxation vector from disc.
C-------------------------------------------------
C
         IF (IOPT .EQ. 3) THEN
C
            CALL DZERO(WORK(KKABAR),LENBAR)
C
            LUBAR0 = -978
            CALL GPOPEN(LUBAR0,'CCKABAR0','UNKNOWN',' ','UNFORMATTED',
     *                  IDUMMY,.FALSE.)
            REWIND(LUBAR0)
            READ(LUBAR0) (WORK(KKABAR+I-1), I = 1,LENBAR)
            CALL GPCLOSE(LUBAR0,'KEEP')
C
C----------------------------------------------------------------
C           Read MO-coefficients from interface file and reorder.
C----------------------------------------------------------------
C
            CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &                  .FALSE.)
            REWIND LUSIFC
            CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
            READ (LUSIFC)
            READ (LUSIFC)
            READ (LUSIFC) (WORK(KCMO+I-1), I=1,NLAMDS)
            CALL GPCLOSE(LUSIFC,'KEEP')
C
            CALL CMO_REORDER(WORK(KCMO),WORK(KEND1),LWRK1)
C
C--------------------------------------------------------------------
C           Calculate ao-transformed zeta-kappa-bar-0 and HF density.
C--------------------------------------------------------------------
C
            KOFDIJ = KKABAR
            KOFDAB = KOFDIJ + NMATIJ(1)
            KOFDAI = KOFDAB + NMATAB(1)
            KOFDIA = KOFDAI + NT1AMX
C
            ISDEN = 1
            CALL DZERO(WORK(KKABAO),N2BST(1))
            CALL CC_DENAO(WORK(KKABAO),ISDEN,WORK(KOFDAI),WORK(KOFDAB),
     *                    WORK(KOFDIJ),WORK(KOFDIA),ISDEN,
     *                    WORK(KCMO),1,WORK(KCMO),1,WORK(KEND1),LWRK1)
C
            CALL CCS_D1AO(WORK(KDHFAO),WORK(KEND1),LWRK1)
            IF (FROIMP .OR. FROEXP) THEN
               MODEL = 'DUMMY'
               CALL CC_FCD1AO(WORK(KDHFAO),WORK(KEND1),LWRK1,MODEL)
            ENDIF
C
C---------------------------------------------------------------
C           Add orbital relaxation for effective density matrix.
C---------------------------------------------------------------
C
            CALL DCOPY(N2BST(1),WORK(KKABAO),1,AODEN,1)
C
         ENDIF
C
C-----------------------------------------------------------
C        Backtransform the one electron density to AO-basis.
C-----------------------------------------------------------
C
         IF (IOPT .EQ. 2) CALL DZERO(AODEN,N2BST(1))
C
         ISDEN = 1
         CALL CC_DENAO(AODEN,ISDEN,WORK(KONEAI),WORK(KONEAB),
     *                 WORK(KONEIJ),WORK(KONEIA),ISDEN,WORK(KLAMDP),1,
     *                 WORK(KLAMDH),1,WORK(KEND1),LWRK1)
C
C------------------------------------------------------
C        Add frozen core contributions to AO densities.
C------------------------------------------------------
C
         IF (FROIMP) THEN
            IF (IOPT .EQ. 3) THEN
C
               KOFFAI = KKABAR + NMATIJ(1) + NMATAB(1) + 2*NT1AMX
               KOFFIA = KOFFAI + NT1FRO(1)
               KOFFIJ = KOFFIA + NT1FRO(1)
               KOFFJI = KOFFIJ + NCOFRO(1)
C
               ISDEN = 1
               ICON  = 1
               CALL CC_D1FCB(AODEN,WORK(KOFFIJ),WORK(KOFFJI),
     *                       WORK(KOFFAI),WORK(KOFFIA),WORK(KEND1),
     *                       LWRK1,ISDEN,ICON)
C
               ISDEN = 1
               ICON  = 2
               CALL CC_D1FCB(WORK(KKABAO),WORK(KOFFIJ),WORK(KOFFJI),
     *                       WORK(KOFFAI),WORK(KOFFIA),WORK(KEND1),
     *                       LWRK1,ISDEN,ICON)
C
            ELSE
C
               MODEL = 'DUMMY'
               CALL CC_FCD1AO(AODEN,WORK(KEND1),LWRK1,MODEL)
C
            ENDIF
         ENDIF
C
C-----------------------------------------------------------
C        Calculate one electron contribution to ccsd-energy.
C-----------------------------------------------------------
C
         WRITE (LUPRI,*) 'AODEN in CC_DEN for IOPT=',IOPT
         CALL CC_PRONELAO(AODEN,1)
C
         ECCSD1 = DDOT(N2BST(ISYMOP),AODEN,1,WORK(KONINT),1)
         ECCSD2 = ZERO
C
      ENDIF
C
C--------------------------------------------------------
C     Calculate M-intermediate of zeta- and t-amplitudes.
C--------------------------------------------------------
C
      CALL CC_MI(WORK(KMINT),WORK(KZ2AM),ISYMTR,WORK(KT2AM),ISYMOP,
     *           WORK(KEND1),LWRK1)
C
C--------------------------------------------------------
C     Calculate resorted M-intermediate M(imjk)->M(mkij). 
C--------------------------------------------------------
C
      CALL CC_MIRS(WORK(KMIRES),WORK(KMINT))
C
      TIMONE = SECOND() - TIMONE
C
C-----------------------------------
C     Start the loop over integrals.
C-----------------------------------
C
      KENDS2 = KEND1
      LWRKS2 = LWRK1
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
     *                  KODPP1,KODPP2,KRDPP1,KRDPP2,
     *                  KFREE,LFREE,KEND1,WORK(KCCFB1),WORK(KINDXB),
     *                  WORK(KEND1),LWRK1,IPRERI)
            KEND1 = KFREE
            LWRK1 = LFREE
         ENDIF
         NTOSYM = 1
      ELSE
         NTOSYM = NSYM
      ENDIF
C
      KENDSV = KEND1
      LWRKSV = LWRK1
C
      ICDEL1 = 0
      DO 100 ISYMD1 = 1,NTOSYM
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
         DO 110 ILLL = 1,NTOT
C
C---------------------------------------------
C           If direct calculate the integrals.
C---------------------------------------------
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
     *                        WORK(KEND1), LWRK1,IPRERI)
               ENDIF
               DTIME  = SECOND() - DTIME
               TIMHE2 = TIMHE2   + DTIME
C
               KRECNR = KEND1
               KEND1  = KRECNR + (NBUFX(0) - 1)/IRAT + 1
               LWRK1  = LWORK  - KEND1
               IF (LWRK1 .LT. 0) THEN
                  CALL QUIT('Insufficient core in CC_DEN')
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
            DO 120 IDEL2 = 1,NUMDIS
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
C----------------------------------------
C              Work space allocation two.
C----------------------------------------
C
               ISYDEN = ISYMD
C
               KD2IJG = KEND1
               KD2AIG = KD2IJG + ND2IJG(ISYDEN)
               KD2IAG = KD2AIG + ND2AIG(ISYDEN)
               KD2ABG = KD2IAG + ND2AIG(ISYDEN)
               KEND2  = KD2ABG + ND2ABG(ISYDEN)
               LWRK2  = LWORK  - KEND2
C
               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:',KEND2
                  CALL QUIT('Insufficient space for allocation '//
     &                      '2 in CC_DEN')
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
               AUTIME = SECOND() - AUTIME
               TIMDEN = TIMDEN + AUTIME
C
C------------------------------------------
C              Work space allocation three.
C------------------------------------------
C
               ISYDIS = MULD2H(ISYMD,ISYMOP)
C
               KXINT  = KEND2
               KEND3  = KXINT  + NDISAO(ISYDIS)
               LWRK3  = LWORK  - KEND3
C
               IF (LWRK3 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:',KEND3
                  CALL QUIT('Insufficient space for allocation '//
     &                      '3 in CC_DEN')
               ENDIF
C
C--------------------------------------------
C              Read AO integral distribution.
C--------------------------------------------
C
               AUTIME = SECOND()
               CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND3),LWRK3,
     *                     WORK(KRECNR),DIRECT)
               AUTIME = SECOND() - AUTIME
               TIRDAO = TIRDAO + AUTIME
C
C----------------------------------------------------------------------
C              Calculate integral intermediates needed for frozen core.
C----------------------------------------------------------------------
C
               IF ((FROIMP) .AND. (IOPT .EQ. 1)) THEN
C
                  KDSRHF = KEND3
                  KOFOIN = KDSRHF + NDSRHF(ISYMD)
                  KDSFRO = KOFOIN + NOFROO(ISYDIS)
                  KSCRAI = KDSFRO + NDSFRO(ISYDIS)
                  KSCAIF = KSCRAI + NOFROO(ISYDIS)
                  KEND3  = KSCAIF + NCOFRF(ISYDIS)
                  LWRK3  = LWORK  - KEND3
C
                  IF (LWRK3 .LT. 0) THEN
                     WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:',KEND3
                     CALL QUIT('Insufficient space for allocation '//
     &                         'in CC_DEN')
                  ENDIF
C
                  CALL DZERO(WORK(KSCRAI),NOFROO(ISYDIS))
                  CALL DZERO(WORK(KSCAIF),NCOFRF(ISYDIS))
C
C-------------------------------------------------------------------------
C                 Transform one index in the integral batch to correlated.
C-------------------------------------------------------------------------
C
                  CALL CCTRBT(WORK(KXINT),WORK(KDSRHF),WORK(KLAMDP),
     *                        ISYMOP,WORK(KEND3),LWRK3,ISYDIS)
C
C---------------------------------------------------------------------
C                 Transform one index in the integral batch to frozen.
C---------------------------------------------------------------------
C
                  CALL CC_GTOFRO(WORK(KXINT),WORK(KDSFRO),WORK(KCMOF),
     *                           WORK(KEND3),LWRK3,ISYDIS)
C
C--------------------------------------------------------------
C                 Calculate integral batch (cor fro | cor del).
C--------------------------------------------------------------
C
                  CALL CC_OFROIN(WORK(KDSRHF),WORK(KOFOIN),WORK(KCMOF),
     *                           WORK(KEND3),LWRK3,ISYDIS)
C
C------------------------------------------------------------------
C                 Calculate terms to ai-part from KOFOIN integrals.
C------------------------------------------------------------------
C
                  CALL CC_FRCOC1(WORK(KSCRAI),WORK(KOFOIN),ISYDIS)
C
C----------------------------------------------------------------
C                 Calculate exchange parts with KDSFRO integrals.
C----------------------------------------------------------------
C
                  CALL CC_FRCOMI(WORK(KSCRAI),WORK(KSCAIF),
     *                           WORK(KDSFRO),WORK(KCMOF),
     *                           WORK(KEND3),LWRK3,ISYMD)
C
C----------------------------------------------------
C                 Calculate coulomb part of aI block.
C----------------------------------------------------
C
                  CALL CC_FRCOF1(WORK(KSCAIF),WORK(KDSFRO),WORK(KCMOF),
     *                           WORK(KEND3),LWRK3,ISYMD)
C
C-----------------------------------------------------
C                 Calculate exchange part of aI block.
C-----------------------------------------------------
C
                  CALL CC_FRCOF2(WORK(KSCAIF),WORK(KDSRHF),WORK(KCMOF),
     *                           WORK(KEND3),LWRK3,ISYMD)
C
C----------------------------------------------------------
C                 Dump intermediates to disc for later use.
C----------------------------------------------------------
C
                  CALL CC_FRINDU(WORK(KSCRAI),WORK(KSCAIF),IDEL,ISYMD,
     &                           LUN1,LUN2)
C
               ENDIF
C
C------------------------------------------------------
C              Start loop over second AO-index (gamma).
C------------------------------------------------------
C
               AUTIME = SECOND()
C
               DO 130 ISYMG = 1,NSYM
C
                  ISYMPQ = MULD2H(ISYMG,ISYDEN)
C
                  DO 140 G = 1,NBAS(ISYMG)
C
                     KD2GIJ = KD2IJG + ID2IJG(ISYMPQ,ISYMG)
     *                      + NMATIJ(ISYMPQ)*(G - 1) 
                     KD2GAI = KD2AIG + ID2AIG(ISYMPQ,ISYMG)
     *                      + NT1AM(ISYMPQ)*(G - 1)
                     KD2GAB = KD2ABG + ID2ABG(ISYMPQ,ISYMG)
     *                      + NMATAB(ISYMPQ)*(G - 1)
                     KD2GIA = KD2IAG + ID2AIG(ISYMPQ,ISYMG)
     *                      + NT1AM(ISYMPQ)*(G - 1)
C
C-----------------------------------------------
C                    Work space allocation four.
C-----------------------------------------------
C
                     KINTAO = KEND3
                     KD2AOB = KINTAO + N2BST(ISYMPQ)
                     KEND4  = KD2AOB + N2BST(ISYMPQ)
                     LWRK4  = LWORK  - KEND4
C
                     IF (LWRK4 .LT. 0) THEN
                        WRITE(LUPRI,*) 'Available:', LWORK
                        WRITE(LUPRI,*) 'Needed:', KEND4
                        CALL QUIT('Insufficient work space in CC_DEN')
                     ENDIF
C
                     CALL DZERO(WORK(KD2AOB),N2BST(ISYMPQ))
C
C-------------------------------------------------------------
C                    Calculate frozen core contributions to d.
C-------------------------------------------------------------
C
                     IF (FROIMP) THEN
C
                        KFD2IJ = KEND4
                        KFD2JI = KFD2IJ + NCOFRO(ISYMPQ)
                        KFD2AI = KFD2JI + NCOFRO(ISYMPQ)
                        KFD2IA = KFD2AI + NT1FRO(ISYMPQ)
                        KFD2II = KFD2IA + NT1FRO(ISYMPQ)
                        KEND4  = KFD2II + NFROFR(ISYMPQ)
                        LWRK4  = LWORK  - KEND4
C
                        IF (LWRK4 .LT. 0) THEN
                           WRITE (LUPRI,*) 'Available:', LWORK
                           WRITE (LUPRI,*) 'Needed:', KEND4
                           CALL QUIT(
     *                          'Insufficient work space in CC_DEN')
                        ENDIF
C
                        CALL DZERO(WORK(KFD2IJ),NCOFRO(ISYMPQ))
                        CALL DZERO(WORK(KFD2JI),NCOFRO(ISYMPQ))
                        CALL DZERO(WORK(KFD2AI),NT1FRO(ISYMPQ))
                        CALL DZERO(WORK(KFD2IA),NT1FRO(ISYMPQ))
                        CALL DZERO(WORK(KFD2II),NFROFR(ISYMPQ))
C
                        CALL CC_FD2BL(WORK(KFD2II),WORK(KFD2IJ),
     *                                WORK(KFD2JI),WORK(KFD2AI),
     *                                WORK(KFD2IA),WORK(KONEIJ),
     *                                WORK(KONEAB),WORK(KONEAI),
     *                                WORK(KONEIA),WORK(KCMOF),
     *                                WORK(KLAMDH),WORK(KLAMDP),
     *                                WORK(KEND4),LWRK4,IDEL,
     *                                ISYMD,G,ISYMG)
C
                        CALL CC_FD2AO(WORK(KD2AOB),WORK(KFD2II),
     *                                WORK(KFD2IJ),WORK(KFD2JI),
     *                                WORK(KFD2AI),WORK(KFD2IA),
     *                                WORK(KCMOF),WORK(KLAMDH),
     *                                WORK(KLAMDP),WORK(KEND4),LWRK4,
     *                                ISYMPQ)
C
                        CALL CC_D2GAF(WORK(KD2GIJ),WORK(KD2GAB),
     *                                WORK(KD2GAI),WORK(KD2GIA),
     *                                WORK(KONEIJ),WORK(KONEAB),
     *                                WORK(KONEAI),WORK(KONEIA),
     *                                WORK(KCMOF),IDEL,ISYMD,G,ISYMG)
C
                     ENDIF
C
C-------------------------------------------------------
C                    Square up AO-integral distribution.
C-------------------------------------------------------
C
                     KOFFIN = KXINT + IDSAOG(ISYMG,ISYDIS) 
     *                      + NNBST(ISYMPQ)*(G - 1) 
C
                     CALL CCSD_SYMSQ(WORK(KOFFIN),ISYMPQ,
     *                               WORK(KINTAO))
C
C------------------------------------------------------------
C                    Backtransform density fully to AO basis.
C------------------------------------------------------------
C
                     IF (IOPT .GE. 2) THEN
C
                        CALL CC_DENAO(WORK(KD2AOB),ISYMPQ,
     *                                WORK(KD2GAI),WORK(KD2GAB),
     *                                WORK(KD2GIJ),WORK(KD2GIA),ISYMPQ,
     *                                WORK(KLAMDP),1,WORK(KLAMDH),1,
     *                                WORK(KEND4),LWRK4)
C
C---------------------------------------------------------------------
C                    Add relaxation terms to set up effective density.
C---------------------------------------------------------------------
C
                        IF (IOPT .EQ. 3) THEN
C
                           ICON = 1
                           CALL CC_D2EFF(WORK(KD2AOB),G,ISYMG,IDEL,
     *                                   ISYMD,WORK(KKABAO),
     *                                   WORK(KDHFAO),ICON)
C
                        ENDIF
C
C----------------------------------------------------------------------
C                    Calcalate the 2 e- density contribution to E-ccsd.
C----------------------------------------------------------------------
C
                        ECCSD2 = ECCSD2 + HALF*DDOT(N2BST(ISYMPQ),
     *                                    WORK(KD2AOB),1,WORK(KINTAO),1)
C
                     ELSE
C
C-----------------------------------------------
C                    Work space allocation five.
C-----------------------------------------------
C
                        KIJINT = KEND4
                        KAIINT = KIJINT + NMATIJ(ISYMPQ)
                        KIAINT = KAIINT + NT1AM(ISYMPQ)
                        KABINT = KIAINT + NT1AM(ISYMPQ)
                        KABTIN = KABINT + NMATAB(ISYMPQ)
                        KIJTIN = KABTIN + NMATAB(ISYMPQ)
                        KD2TAB = KIJTIN + NMATIJ(ISYMPQ)
                        KD2TIJ = KD2TAB + NMATAB(ISYMPQ)
                        KEND5  = KD2TIJ + NMATIJ(ISYMPQ)
                        LWRK5  = LWORK  - KEND5
                        IF (FROIMP) THEN
                           KIIFRO = KEND5
                           KIJFRO = KIIFRO + NFROFR(ISYMPQ)
                           KJIFRO = KIJFRO + NCOFRO(ISYMPQ)
                           KAIFRO = KJIFRO + NCOFRO(ISYMPQ)
                           KIAFRO = KAIFRO + NT1FRO(ISYMPQ)
                           KEND5  = KIAFRO + NT1FRO(ISYMPQ)
                           LWRK5  = LWORK  - KEND5
                        ENDIF
C
                        IF (LWRK5 .LT. 0) THEN
                           WRITE(LUPRI,*) 'Available:', LWORK
                           WRITE(LUPRI,*) 'Needed:', KEND5
                           CALL QUIT('Insufficient work space '//
     &                               'in CC_DEN')
                        ENDIF
C
C----------------------------------------------------------------
C                       Transform 2 integral indices to MO-basis.
C----------------------------------------------------------------
C
                        ISYM = ISYMPQ
                        CALL CCDINTMO(WORK(KIJINT),WORK(KIAINT),
     *                                WORK(KABINT),WORK(KAIINT),
     *                                WORK(KINTAO),WORK(KLAMDP),
     *                                WORK(KLAMDH),WORK(KEND5),
     *                                LWRK5,ISYM)
C
                        IF (FROIMP) THEN
C
                           ISYM = ISYMPQ
                           CALL CCIFROMO(WORK(KIJFRO),WORK(KJIFRO),
     *                                   WORK(KAIFRO),WORK(KIAFRO),
     *                                   WORK(KIIFRO),WORK(KINTAO),
     *                                   WORK(KLAMDP),WORK(KLAMDH),
     *                                   WORK(KCMOF),WORK(KEND5),
     *                                   LWRK5,ISYM)
C
                        ENDIF
C
C-----------------------------------------------------------------
C                       Set up transposed integrals and densities.
C-----------------------------------------------------------------
C
                        CALL DCOPY(NMATAB(ISYMPQ),WORK(KABINT),1,
     *                             WORK(KABTIN),1)
                        CALL DCOPY(NMATIJ(ISYMPQ),WORK(KIJINT),1,
     *                             WORK(KIJTIN),1)
                        CALL DCOPY(NMATAB(ISYMPQ),WORK(KD2GAB),1,
     *                             WORK(KD2TAB),1)
                        CALL DCOPY(NMATIJ(ISYMPQ),WORK(KD2GIJ),1,
     *                             WORK(KD2TIJ),1)
C
                        CALL CC_EITR(WORK(KABTIN),WORK(KIJTIN),
     *                               WORK(KEND5),LWRK5,ISYMPQ)
                        CALL CC_EITR(WORK(KD2TAB),WORK(KD2TIJ),
     *                               WORK(KEND5),LWRK5,ISYMPQ)
C
C-------------------------------------------------------------------
C                       Calculate 2 e- contribution to Zeta-Kappa-0.
C-------------------------------------------------------------------
C
                        ISYM = ISYMPQ
                        CALL CCDENZK0(AODEN,WORK(KIJINT),WORK(KAIINT),
     *                                WORK(KIAINT),WORK(KABINT),
     *                                WORK(KD2GIJ),WORK(KD2GAI),
     *                                WORK(KD2GIA),WORK(KD2GAB),
     *                                WORK(KIJTIN),WORK(KABTIN),
     *                                WORK(KD2TIJ),WORK(KD2TAB),
     *                                WORK(KT1AM),WORK(KEND5),LWRK5,
     *                                ISYM)
C
                        IF (FROIMP) THEN
C
                           ISYM = ISYMPQ
C
                           CALL CCFRETAI(AODEN,WORK(KIJFRO),
     *                                   WORK(KJIFRO),WORK(KAIFRO),
     *                                   WORK(KIAFRO),WORK(KFD2IJ),
     *                                   WORK(KFD2JI),WORK(KFD2AI),
     *                                   WORK(KFD2IA),WORK(KT1AM),
     *                                   WORK(KEND5),LWRK5,ISYM)
C
C-----------------------------------------------------------------------
C                          Calculate two-electron contribution to right-
C                          hand-side of correlated-frozen multipliers.
C-----------------------------------------------------------------------
C
                           ICON = 2
                           KOFRES = 2*NT1AMX + NMATIJ(1) + NMATAB(1)
     *                            + 2*NT1FRO(1) + 1
C
                           CALL CC_ETIJF(ZKDIA(KOFRES),WORK(KD2GIJ),
     *                                   WORK(KD2GAB),WORK(KD2GAI),
     *                                   WORK(KD2GIA),WORK(KD2TIJ),
     *                                   WORK(KFD2II),WORK(KFD2IJ),
     *                                   WORK(KFD2JI),WORK(KFD2AI),
     *                                   WORK(KFD2IA),WORK(KIJINT),
     *                                   WORK(KAIINT),WORK(KIAINT),
     *                                   WORK(KIJTIN),WORK(KABTIN),
     *                                   WORK(KIJFRO),WORK(KJIFRO),
     *                                   WORK(KAIFRO),WORK(KIAFRO),
     *                                   WORK(KIIFRO),WORK(KT1AM),
     *                                   WORK(KEND5),LWRK5,ISYM,ICON)
C
C-----------------------------------------------------------------------
C                          Calculate two-electron contribution to right-
C                          hand-side of virtual-frozen multipliers.
C-----------------------------------------------------------------------
C
                           ICON = 2
                           KOFRE = 2*NT1AMX + NMATIJ(1) + NMATAB(1) + 1
C
                           CALL CC_ETAIF(ZKDIA(KOFRE),WORK(KD2GAB),
     *                                   WORK(KD2GAI),WORK(KD2GIA),
     *                                   WORK(KD2TIJ),WORK(KD2TAB),
     *                                   WORK(KFD2II),WORK(KFD2IJ),
     *                                   WORK(KFD2JI),WORK(KFD2AI),
     *                                   WORK(KFD2IA),WORK(KIJINT),
     *                                   WORK(KABINT),WORK(KAIINT),
     *                                   WORK(KIAINT),WORK(KABTIN),
     *                                   WORK(KIJFRO),WORK(KJIFRO),
     *                                   WORK(KAIFRO),WORK(KIAFRO),
     *                                   WORK(KIIFRO),WORK(KT1AM),
     *                                   WORK(KEND5),LWRK5,ISYM,ICON)
C
C-----------------------------------------------------------------------
C                          Calculate two-electron contribution to right-
C                          hand-side of frozen-frozen multipliers.
C-----------------------------------------------------------------------
C
                           ICON = 2
                           KOFRES = 2*NT1AMX + NMATIJ(1) + NMATAB(1)
     *                            + 2*NT1FRO(1) + 2*NCOFRO(1) + 1
C
c                          CALL CC_ETIFJF(ZKDIA(KOFRES),WORK(KFD2II),
c    *                                   WORK(KFD2IJ),WORK(KFD2JI),
c    *                                   WORK(KFD2AI),WORK(KFD2IA),
c    *                                   WORK(KIJFRO),WORK(KJIFRO),
c    *                                   WORK(KAIFRO),WORK(KIAFRO),
c    *                                   WORK(KIIFRO),ISYM,ICON)
C
c                          CALL CCDIAZK0(ZKDIA,WORK(KIJINT),
c    *                                   WORK(KAIINT),
c    *                                   WORK(KIAINT),WORK(KABINT),
c    *                                   WORK(KD2GIJ),WORK(KD2GAI),
c    *                                   WORK(KD2GIA),WORK(KD2GAB),
c    *                                   WORK(KIJTIN),WORK(KABTIN),
c    *                                   WORK(KD2TIJ),WORK(KD2TAB),
c    *                                   WORK(KT1AM),WORK(KEND5),
c    *                                   LWRK5,ISYM)
C
c                          CALL CCFRETAB(ZKDIA,WORK(KIJFRO),
c    *                                   WORK(KJIFRO),WORK(KAIFRO),
c    *                                   WORK(KIAFRO),WORK(KFD2IJ),
c    *                                   WORK(KFD2JI),WORK(KFD2AI),
c    *                                   WORK(KFD2IA),WORK(KT1AM),
c    *                                   WORK(KEND5),LWRK5,ISYM)
C
c                          CALL CCFRETIJ(ZKDIA,WORK(KIJFRO),
c    *                                   WORK(KJIFRO),WORK(KAIFRO),
c    *                                   WORK(KIAFRO),WORK(KFD2IJ),
c    *                                   WORK(KFD2JI),WORK(KFD2AI),
c    *                                   WORK(KFD2IA),WORK(KT1AM),
c    *                                   WORK(KEND5),LWRK5,ISYM)
C
                        ENDIF
                     ENDIF
C
  140                CONTINUE
  130             CONTINUE
C
                  AUTIME = SECOND() - AUTIME
                  TIMRES = TIMRES + AUTIME
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
C-----------------------
C     Regain work space.
C-----------------------
C
      KEND1 = KENDS2
      LWRK1 = LWRKS2
C
C---------------------------------------------------------------
C     Add nuclear nuclear repulsion energy and write out E-ccsd.
C---------------------------------------------------------------
C
      IF (IOPT .GT. 1) THEN
C
         ECCSD = ECCSD1 + ECCSD2 + POTNUC
C
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) 'Coupled Cluster energy constructed'
         WRITE(LUPRI,*) 'from density matrices:'
         IF (CCD)  WRITE(LUPRI,*) 'CCD-energy:',  ECCSD
         IF (CCSD) WRITE(LUPRI,*) 'CCSD-energy:', ECCSD
         WRITE(LUPRI,*) 'H1 energy, ECCSD1 = ',ECCSD1
         WRITE(LUPRI,*) 'H2 energy, ECCSD2 = ',ECCSD2
         WRITE(LUPRI,*) 'Nuc. Pot. energy  = ',POTNUC
C
      ENDIF
C
C------------------------------------------------
C     Write out frozen block of eta-kappa-bar-0.
C------------------------------------------------
C
      IF ((IPRINT .GT. 9) .AND. (FROIMP) .AND. (IOPT .EQ. 1)) THEN
         CALL AROUND('Eta-kappa-bar-0 correlated-frozen block')
         KOFRES = NMATIJ(1) + NMATAB(1) + 2*NT1AMX + 2*NT1FRO(1) + 1
         ZKAPF1 = DDOT(NCOFRO(1),ZKDIA(KOFRES),1,
     *                 ZKDIA(KOFRES),1)
         ZKAPFR = ZKAPF1**0.5
         WRITE(LUPRI,*) 'Norm of correlated-frozen block:', ZKAPFR
      ENDIF
C
C--------------------------------------------------
C     Calculate the diagonal blocks of kappa-bar-0.
C--------------------------------------------------
C
c     IF ((FROIMP) .AND. (IOPT .EQ. 1)) THEN
c        CALL CCSD_ZKBLO(ZKDIA,WORK(KEND1),LWRK1)
c     ENDIF
C
C-----------------------------------------------------------
C     Calculate the correlated-frozen blocks of kappa-bar-0.
C-----------------------------------------------------------
C
      IF ((FROIMP) .AND. (IOPT .EQ. 1)) THEN
         KOFFIJ = 2*NT1AMX + NMATIJ(1) + NMATAB(1) + 2*NT1FRO(1) + 1
         CALL CC_ZKFCB(ZKDIA(KOFFIJ),WORK(KEND1),LWRK1)
C
C----------------------------------------------------------------
C        Calculate correction terms from correlated-frozen block.
C----------------------------------------------------------------
C
         KOFFAI = 2*NT1AMX + NMATIJ(1) + NMATAB(1) + 1
         CALL CC_ETACOR(AODEN,ZKDIA(KOFFAI),ZKDIA(KOFFIJ),
     *                  WORK(KCMOF),LUN1,LUN2,WORK(KEND1),LWRK1)
C
      ENDIF
C
C---------------------
C     Reorder results.
C---------------------
C
      KOFFAI = 2*NT1AMX + NMATIJ(1) + NMATAB(1) + 1
      CALL CC_ETARE(AODEN,ZKDIA(KOFFAI),WORK(KEND1),LWRK1)
C
C---------------------------------
C     Write out eta-ai and eta-aI.
C---------------------------------
C
      IF (IPRINT .GT. 20) THEN
C
         CALL AROUND('Eta-kappa-bar-0-ai vector exiting CC_DEN')
C
         DO 20 ISYM = 1,NSYM
C
            WRITE(LUPRI,*) ' '
            WRITE(LUPRI,444) 'Sub-symmetry block number:', ISYM
            WRITE(LUPRI,555) '--------------------------'
  444       FORMAT(3X,A26,2X,I1)
  555       FORMAT(3X,A25)
C
            KOFF = IALLAI(ISYM,ISYM) + 1
            CALL OUTPUT(AODEN(KOFF),1,NVIR(ISYM),1,NRHFS(ISYM),
     *                  NVIR(ISYM),NRHFS(ISYM),1,LUPRI)
C
            IF ((NVIR(ISYM) .EQ. 0) .OR. (NRHFS(ISYM) .EQ. 0)) THEN
               WRITE(LUPRI,*) 'This sub-symmetry is empty'
            ENDIF
C
  20     CONTINUE
      ENDIF
C
      IF (IPRINT .GT. 9) THEN
         ETAKAN = DDOT(NALLAI(1),AODEN,1,AODEN,1)
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) 'Norm of occupied-virtual block:', ETAKAN
      ENDIF
C
C-------------------------------------------
C     Write out frozen block of kappa-bar-0.
C-------------------------------------------
C
      IF ((IPRINT .GT. 9) .AND. (FROIMP) .AND. (IOPT .EQ. 1)) THEN
         CALL AROUND('Zeta-kappa-0 correlated-frozen block')
         KOFRES = 2*NT1AMX + NMATIJ(1) + NMATAB(1) + 2*NT1FRO(1) + 1
         ZKAPF1 = DDOT(NCOFRO(1),ZKDIA(KOFRES),1,
     *                 ZKDIA(KOFRES),1)
         ZKAPFR = ZKAPF1**0.5
         WRITE(LUPRI,*) 'Norm of correlated-frozen block:', ZKAPFR
      ENDIF
C
C------------------------------
C     Close intermediate files.
C------------------------------
C
      IF ((FROIMP) .AND. (IOPT .EQ. 1)) THEN
C
         CALL WCLOSE2(LUN1,NAME1,'KEEP')
         CALL WCLOSE2(LUN2,NAME2,'KEEP')
      ENDIF
C
C-----------------------
C     Write out timings.
C-----------------------
C
  99  TIMTOT = SECOND() - TIMTOT
C
      IF (IPRINT .GT. 3) THEN
         WRITE (LUPRI,*) ' '
         WRITE (LUPRI,*) 'Requested density dependent '//
     &        'quantities calculated'
         WRITE (LUPRI,*) 'Total time used in CC_DEN:', TIMTOT
      ENDIF
      IF (IPRINT .GT. 9) THEN
         WRITE (LUPRI,*) 'Time used for setting up 2 e- density:',TIMDEN
         WRITE (LUPRI,*) 'Time used for contraction with integrals:',
     &        TIMRES
         WRITE (LUPRI,*) 'Time used for reading 2 e- AO-integrals:',
     &        TIRDAO
         WRITE (LUPRI,*) 'Time used for calculating 2 e- AO-integrals:',
     *              TIMHE2
         WRITE (LUPRI,*) 'Time used for 1 e- density & intermediates:',
     *              TIMONE
      ENDIF
C
      CALL QEXIT('CC_DEN')
      RETURN
      END
C  /* Deck cc_den2 */
      SUBROUTINE CC_DEN2(D2IJG,D2AIG,D2IAG,D2ABG,Z2AM,T2AM,T2AMT,XMINT,
     *                   XMAT,YMAT,DAB,DAI,DIA,XMIRES,
     *                   XLAMDH,ISYMLH,XLAMDP,ISYMLP,
     *                   WORK,LWORK,IDEL,ISYMD)
C
C     Written by Asger Halkier 19/2 - 1996
C
C     Version: 1.0
C
C     Purpose: Directs the calculation of the 2 electron CC density
C              d(pq,gam;del) for a given del (IDEL). The 4 blocks pq
C              of the result is returned in D2IJG through D2ABG!
C
C     Modifications for CC2 by A. Halkier and S. Coriani 13/01-2000.
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      INTEGER ISYMLP, ISYMLH
      DIMENSION D2IJG(*), D2AIG(*), D2IAG(*), D2ABG(*), Z2AM(*)
      DIMENSION T2AM(*), T2AMT(*), XMINT(*), XMAT(*), YMAT(*)
      DIMENSION DAB(*), DAI(*), DIA(*), XMIRES(*), XLAMDH(*)
      DIMENSION XLAMDP(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
C
      CALL QENTER('CC_DEN2')
      IF (DEBUG) WRITE(LUPRI,*) 'Entered CC_DEN2'
C
C-------------------------------
C     set some symmetries:
C-------------------------------
C
      ISYD1  = 1                     ! one-particle density
      ISYMTR = 1                     ! Z2AM
      ISYD2H = MULD2H(ISYMD,ISYMLH)  ! lambdah backtransformed density
      ISYDEN = MULD2H(ISYD2H,ISYMLP) ! lambdah + lambdap transformed 
      ISYMTI = MULD2H(ISYMLH,ISYMD)
C
C-------------------------------
C     Work space allocation one.
C-------------------------------
C
      KT2DEL = 1
      KEND1  = KT2DEL + NT2BCD(ISYMTI)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient space for first allocation in CC_DEN2')
      ENDIF
C
C------------------------------------------------------
C     Calculate the delta backtransformed t2-amplitude.
C------------------------------------------------------
C
      CALL CC_TI(WORK(KT2DEL),ISYMTI,T2AM,ISYMOP,XLAMDH,ISYMLH,
     *           WORK(KEND1),LWRK1,IDEL,ISYMD)
C
C---------------------------------------------------
C     Calculate the (occ.occ,vir;del) density block.
C---------------------------------------------------
C
      KD2IJA = KEND1
      KEND2  = KD2IJA + NMAIJA(ISYD2H)
      LWRK2  = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE (LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
         CALL QUIT(
     *        'Insufficient space for second allocation in CC_DEN2')
      ENDIF
C
      CALL DZERO(WORK(KD2IJA),NMAIJA(ISYD2H))
C
C------------------------------------------------
C     Contributions from projection against <u1|.
C     Remember that D_ai = zeta_ai.
C------------------------------------------------
C
      CALL DIJA11(WORK(KD2IJA),ISYD2H,DAI,ISYD1,XLAMDH,ISYMLH,
     &            WORK(KEND2),LWRK2,IDEL,ISYMD)
C
C------------------------------------------------
C     Contributions from projection against <u2|.
C------------------------------------------------
C
      IF (CCSD) THEN
         CALL DIJA21(WORK(KD2IJA),ISYD2H,DAB,ISYD1,XLAMDH,ISYMLH,
     *               WORK(KEND2),LWRK2,IDEL,ISYMD)
         CALL DIJA22(WORK(KD2IJA),ISYD2H,Z2AM,WORK(KT2DEL),ISYMTI,
     *               WORK(KEND2),LWRK2)
         CALL DIJA23(WORK(KD2IJA),ISYD2H,Z2AM,WORK(KT2DEL),ISYMTI,
     *               WORK(KEND2),LWRK2)
      ENDIF
C
C-----------------------------------------------------------------
C     Backtransform third virtual index to AO and store in result.
C-----------------------------------------------------------------
C
      ICON = 4
      CALL CCD2_PQAO(D2IJG,ISYDEN,WORK(KD2IJA),ISYD2H,
     &               XLAMDP,ISYMLP,ICON)
C
      IF (CCSD) THEN
C
C-------------------------------------------
C     Calculate the (vir.vir,vir;del) blocks
C     contribution to D2ABG directly.
C-------------------------------------------
C
         DO 100 ISYMG = 1,NSYM
C
            ISYZ2I = MULD2H(MULD2H(ISYMG,ISYMTR),ISYMLP)
C
            KZ2GAM = KEND1
            KEND2  = KZ2GAM + NT2BCD(ISYZ2I)
            LWRK2  = LWORK  - KEND2
C
            IF (LWRK2 .LT. 0) THEN
               WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
               CALL QUIT
     &          ('Insufficient space for 3. allocation in CC_DEN2')
            ENDIF
C
            DO 110 G = 1,NBAS(ISYMG)
C
C---------------------------------------------------------------
C           Calculate the gamma backtransformed zeta2-amplitude.
C---------------------------------------------------------------
C
               ICON = 2
               CALL CC_T2AO(Z2AM,XLAMDP,ISYMLP,WORK(KZ2GAM),
     *                      WORK(KEND2),LWRK2,
     *                      G,ISYMG,ISYMTR,ICON)
C
C--------------------------------------
C           Calculate the contribution.
C--------------------------------------
C
               CALL DABGAO(D2ABG,ISYDEN,WORK(KT2DEL),ISYMTI,
     *                     WORK(KZ2GAM),ISYZ2I,WORK(KEND2),LWRK2,
     *                     G,ISYMG)
C
  110       CONTINUE
  100    CONTINUE
C
      ENDIF
C
C------------------------------------------------------------
C     Calculate terms of (occ.vir,occ;del) block with t2-del.
C     Note that the Z-intermediate is overwritten by the 
C     last calculation!
C------------------------------------------------------------
C
      ISYMZI = MULD2H(ISYMTI,ISYMTR)
C
      KD2IAJ = KEND1
      KZINT  = KD2IAJ + NT2BCD(ISYD2H)
      KEND2  = KZINT  + NT2BCD(ISYMZI)
      LWRK2  = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
         CALL QUIT(
     *        'Insufficient space for fourth allocation in CC_DEN2')
      ENDIF
C
      CALL DZERO(WORK(KD2IAJ),NT2BCD(ISYD2H))
C
      IF (CCSD) THEN
C
         ICON = 1
         CALL CC_ZWVI(WORK(KZINT),Z2AM,ISYMTR,WORK(KT2DEL),ISYMTI,
     *                  WORK(KEND2),LWRK2,ICON)
C
         CALL DIAJ25(WORK(KD2IAJ),ISYD2H,WORK(KT2DEL),ISYMTI,XMIRES,
     *               WORK(KEND2),LWRK2)
         CALL DIAJ27(WORK(KD2IAJ),ISYD2H,T2AM,WORK(KZINT),ISYMZI,
     *               WORK(KEND2),LWRK2)
C
         ICON = 2
         CALL CCSD_T2TP(Z2AM,WORK(KEND2),LWRK2,ISYMTR)
         CALL CC_ZWVI(WORK(KZINT),Z2AM,ISYMTR,WORK(KT2DEL),ISYMTI,
     *                  WORK(KEND2),LWRK2,ICON)
         CALL DIAJ27(WORK(KD2IAJ),ISYD2H,T2AM,WORK(KZINT),ISYMZI,
     *               WORK(KEND2),LWRK2)
         CALL CCSD_T2TP(Z2AM,WORK(KEND2),LWRK2,ISYMTR)
C
      ENDIF
C
C--------------------------------------------------------------
C     Calculate delta backtransformed T2AMT & concomitant ZINT.
C     Note that these overwrite the old intermediates, since 
C     all terms in need of the old ones have been calculated!
C--------------------------------------------------------------
C
      CALL CC_TI(WORK(KT2DEL),ISYMTI,T2AMT,ISYMOP,XLAMDH,ISYMLH,
     *           WORK(KEND2),LWRK2,IDEL,ISYMD)
C
      IF (CCSD) THEN
C
         ICON = 1
         CALL CC_ZWVI(WORK(KZINT),Z2AM,ISYMTR,WORK(KT2DEL),ISYMTI,
     *                  WORK(KEND2),LWRK2,ICON)
C
C--------------------------------------------------------------------
C     Calculate terms of (occ.vir,occ;del) block with Zbar (in ZINT).
C--------------------------------------------------------------------
C
         CALL DIAJ26(WORK(KD2IAJ),ISYD2H,T2AM,WORK(KZINT),ISYMZI,
     *               WORK(KEND2),LWRK2)
         CALL DIAJ28(WORK(KD2IAJ),ISYD2H,T2AMT,WORK(KZINT),ISYMZI,
     *               WORK(KEND2),LWRK2)
C
      ENDIF
C
C-----------------------------------------------------------------
C     Calculate the remaining terms of (occ.vir,occ;del) block
C     and calculate the remainder of the (vir.occ,occ;del) block.
C     Note that Zbar is the first contribution to this last block.
C-----------------------------------------------------------------
C
      KD2AIJ = KZINT
      KEND2  = KD2AIJ + NT2BCD(ISYD2H)
      LWRK2  = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
         CALL QUIT(
     *        'Insufficient space for fourth allocation in CC_DEN2')
      ENDIF
C
      IF (CC2) CALL DZERO(WORK(KD2AIJ),NT2BCD(ISYD2H))
C
C------------------------------------------------
C     Contributions from projection against <HF|.
C------------------------------------------------
C
      IF (ISYMTI.NE.ISYD2H) CALL QUIT('Symmetry mismatch in CC_DEN2')
      CALL DAXPY(NT2BCD(ISYMTI),TWO,WORK(KT2DEL),1,WORK(KD2IAJ),1)
C
C------------------------------------------------
C     Contributions from projection against <u1|.
C------------------------------------------------
C
      CALL DAIJ11(WORK(KD2AIJ),ISYD2H,DAI,ISYD1,XLAMDH,ISYMLH,
     *            IDEL,ISYMD)
      CALL DAIJ12(WORK(KD2AIJ),ISYD2H,DAI,ISYD1,XLAMDH,ISYMLH,
     *            WORK(KEND2),LWRK2,IDEL,ISYMD)
      CALL DIAJ11(WORK(KD2IAJ),ISYD2H,DIA,ISYD1,XLAMDH,ISYMLH,
     *            IDEL,ISYMD)
      CALL DIAJ13(WORK(KD2IAJ),ISYD2H,DAI,ISYD1,T2AMT,XLAMDH,ISYMLH,
     *            WORK(KEND2),LWRK2,IDEL,ISYMD)
C
C------------------------------------------------
C     Contributions from projection against <u2|.
C------------------------------------------------
C
      IF (CCSD) THEN
C
         CALL DAIJ22(WORK(KD2AIJ),ISYD2H,DAB,ISYD1,XLAMDH,ISYMLH,
     *               WORK(KEND2),LWRK2,IDEL,ISYMD)
         CALL DIAJ21(WORK(KD2IAJ),ISYD2H,YMAT,WORK(KT2DEL),ISYMTI)
         CALL DIAJ22(WORK(KD2IAJ),ISYD2H,XMAT,WORK(KT2DEL),ISYMTI)
         CALL DIAJ23(WORK(KD2IAJ),ISYD2H,XMAT,WORK(KT2DEL),ISYMTI)
         CALL DIAJ24(WORK(KD2IAJ),ISYD2H,T2AMT,DAB,ISYD1,XLAMDH,ISYMLH,
     *               WORK(KEND2),LWRK2,IDEL,ISYMD)
C
      ENDIF
C
C-------------------------------------------------------------
C     Backtransform occupied index to AO and store in results.
C-------------------------------------------------------------
C
      ICON = 2
      CALL CCD2_PQAO(D2IAG,ISYDEN,WORK(KD2IAJ),ISYD2H,
     *               XLAMDP,ISYMLP,ICON)
      CALL CCD2_PQAO(D2AIG,ISYDEN,WORK(KD2AIJ),ISYD2H,
     *               XLAMDP,ISYMLP,ICON)
C
C--------------------------------------------------
C     Calculate the three blocks: (occ.vir,vir;del)
C     (vir.occ,vir;del) & (vir.vir,occ;del).
C--------------------------------------------------
C
      KD2IAB = KEND1
      KD2ABI = KD2IAB + NCKATR(ISYD2H)
      KD2AIB = KD2ABI + NCKASR(ISYD2H)
      KEND2  = KD2AIB + NCKATR(ISYD2H)
      LWRK2  = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
         CALL QUIT('Insufficient space for fifth allocation in CC_DEN2')
      ENDIF
C
      CALL DZERO(WORK(KD2IAB),NCKATR(ISYD2H))
      CALL DZERO(WORK(KD2ABI),NCKASR(ISYD2H))
      CALL DZERO(WORK(KD2AIB),NCKATR(ISYD2H))
C
C--------------------------------------------------------------
C     Calculate backtransformed zeta2-amplitude zeta(ai,b;del).
C     Note that this equals the d(ai,b;del) density block.
C--------------------------------------------------------------
C
      CALL DAIB21(WORK(KD2AIB),ISYD2H,Z2AM,XLAMDH,ISYMLH,
     *            WORK(KEND2),LWRK2,IDEL,ISYMD)
C
C-------------------------------------------------------
C     Backtransform third index of the (vir.occ,vir;del)
C     block to AO-basis and store in result.
C-------------------------------------------------------
C
      ICON = 5
      CALL CCD2_PQAO(D2AIG,ISYDEN,WORK(KD2AIB),ISYD2H,
     *               XLAMDP,ISYMLP,ICON)
C
C-----------------------------------------------
C     Calculate the zeta(ai,b;del)-containing 
C     contributions to the remaining two blocks.
C-----------------------------------------------
C
      KZ2AO   = KD2AIB
      ISYZ2AO = ISYD2H
C
      IF (CCSD) THEN
C
         CALL DIAC22(WORK(KD2IAB),ISYD2H,T2AMT,WORK(KZ2AO),ISYZ2AO,
     *               WORK(KEND2),LWRK2)
         CALL DABI22(WORK(KD2ABI),ISYD2H,T2AM,WORK(KZ2AO),ISYZ2AO,
     *               WORK(KEND2),LWRK2)
         CALL DABI23(WORK(KD2ABI),ISYD2H,T2AM,WORK(KZ2AO),ISYZ2AO,
     *               WORK(KEND2),LWRK2)
C
      ENDIF
C
C--------------------------------------------------------------------
C     Calculate remaining contributions from projection against <u1|.
C--------------------------------------------------------------------
C
      CALL DABI11(WORK(KD2ABI),ISYD2H,DAI,ISYD1,WORK(KT2DEL),ISYMTI)
      CALL DIAC11(WORK(KD2IAB),ISYD2H,DAI,ISYD1,WORK(KT2DEL),ISYMTI)
C
C--------------------------------------------------------------------
C     Calculate remaining contributions from projection against <u2|.
C--------------------------------------------------------------------
C
      IF (CCSD) THEN
C
         CALL DABI21(WORK(KD2ABI),ISYD2H,DAB,ISYD1,XLAMDH,ISYMLH,
     &               IDEL,ISYMD)
         CALL DIAC21(WORK(KD2IAB),ISYD2H,YMAT,ISYD1,XLAMDH,ISYMLH,
     &               IDEL,ISYMD)
C
      ENDIF
C
C---------------------------------------------------
C     Backtransform third index of the two remaining
C     blocks to AO and store in results.
C---------------------------------------------------
C
      ICON = 5
      CALL CCD2_PQAO(D2IAG,ISYDEN,WORK(KD2IAB),ISYD2H,
     *               XLAMDP,ISYMLP,ICON)
      ICON = 3
      CALL CCD2_PQAO(D2ABG,ISYDEN,WORK(KD2ABI),ISYD2H,
     *               XLAMDP,ISYMLP,ICON)
C
C---------------------------------------------------
C     Calculate the (occ.occ,occ;del) density block.
C---------------------------------------------------
C
      KD2IJK = KEND1
      KEND2  = KD2IJK + NMAIJK(ISYD2H)
      LWRK2  = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
         CALL QUIT('Insufficient space for seventh allocation '//
     &             'in CC_DEN2')
      ENDIF
C
      CALL DZERO(WORK(KD2IJK),NMAIJK(ISYD2H))
C
C------------------------------------------------
C     Contributions from projection against <HF|.
C------------------------------------------------
C
      CALL DIJK01(WORK(KD2IJK),ISYD2H,XLAMDH,ISYMLH,IDEL,ISYMD)
      CALL DIJK02(WORK(KD2IJK),ISYD2H,XLAMDH,ISYMLH,IDEL,ISYMD)
C
C------------------------------------------------
C     Contributions from projection against <u1|.
C------------------------------------------------
C
      CALL DIJK11(WORK(KD2IJK),ISYD2H,XLAMDH,ISYMLH,DIA,ISYD1,
     &            WORK(KEND2),LWRK2,IDEL,ISYMD)
      CALL DIJK13(WORK(KD2IJK),ISYD2H,DAI,ISYD1,WORK(KT2DEL),ISYMTI)
C
C------------------------------------------------
C     Contributions from projection against <u2|.
C------------------------------------------------
C
      IF (CCSD) THEN
C
         CALL DIJK21(WORK(KD2IJK),ISYD2H,XMAT,XLAMDH,ISYMLH,
     *               WORK(KEND2),LWRK2,IDEL,ISYMD)
         CALL DIJK23(WORK(KD2IJK),ISYD2H,XMAT,XLAMDH,ISYMLH,IDEL,ISYMD)
         CALL DIJK24(WORK(KD2IJK),ISYD2H,XMAT,XLAMDH,ISYMLH,IDEL,ISYMD)
         CALL DIJK25(WORK(KD2IJK),ISYD2H,XMINT,XLAMDH,ISYMLH,IDEL,ISYMD)
C
      ENDIF
C
C------------------------------------------------------------------
C     Backtransform third occupied index to AO and store in result.
C------------------------------------------------------------------
C
      ICON = 1
      CALL CCD2_PQAO(D2IJG,ISYDEN,WORK(KD2IJK),ISYD2H,
     *               XLAMDP,ISYMLP,ICON)
C
      IF (DEBUG) WRITE(LUPRI,*) 'Leave CC_DEN2'
      CALL QEXIT('CC_DEN2')
C
      RETURN
      END
C  /* Deck ccd2_pqao */
      SUBROUTINE CCD2_PQAO(D2PQG,ISYPQG,D2PQMO,ISYDMO,
     *                     XLAMDP,ISYMLP,ICON)
C
C     Written by Asger Halkier 19/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the backtransformation of the CC density
C              d(pq,r;del) to d(pq,gam;del), where p, q, & r are 
C              passed through the variable ICON!
C
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>

#if defined (SYS_CRAY)
      REAL ZERO, ONE
#else
      DOUBLE PRECISION ZERO, ONE
#endif
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)

      INTEGER ISYPQG, ISYDMO, ISYMLP, ICON

#if defined (SYS_CRAY)
      REAL D2PQG(*), D2PQMO(*), XLAMDP(*)
#else
      DOUBLE PRECISION D2PQG(*), D2PQMO(*), XLAMDP(*)
#endif
    
      INTEGER ISYMK, ISYMG, ISYMA, ISYMB, KOFF1, KOFF2, KOFF3,
     *        NTOTIJ, NTOTAI, NTOTAB, NTOTG, ISYMIJ, ISYMAB, ISYMAI
C
      CALL QENTER('CCD2_PQAO')
C
      IF (MULD2H(ISYDMO,ISYMLP).NE.ISYPQG) THEN
       CALL QUIT('Symmetry mismatch in CCD2_PQAO')
      END IF
C
C-----------------------------------------------------------------
C     Calculate the transformation of the (occ.occ,occ;del) block.
C-----------------------------------------------------------------
C
      IF (ICON .EQ. 1) THEN
C
         DO 100 ISYMIJ = 1,NSYM
C
            ISYMK  = MULD2H(ISYMIJ,ISYDMO)
            ISYMG  = MULD2H(ISYMK,ISYMLP)
C
            KOFF1  = IMAIJK(ISYMIJ,ISYMK) + 1
            KOFF2  = IGLMRH(ISYMG,ISYMK) + 1
            KOFF3  = ID2IJG(ISYMIJ,ISYMG) + 1
C
            NTOTIJ = MAX(NMATIJ(ISYMIJ),1)
            NTOTG  = MAX(NBAS(ISYMG),1)
C
            CALL DGEMM('N','T',NMATIJ(ISYMIJ),NBAS(ISYMG),NRHF(ISYMK),
     *                 ONE,D2PQMO(KOFF1),NTOTIJ,XLAMDP(KOFF2),NTOTG,
     *                 ONE,D2PQG(KOFF3),NTOTIJ)
C
  100    CONTINUE
C
C---------------------------------------------------------------
C     Calculate the transformation of both the (occ.vir,occ;del) 
C     block and the (vir.occ,occ;del) block (stored equally).
C---------------------------------------------------------------
C
      ELSE IF (ICON .EQ. 2) THEN
C
         DO 110 ISYMAI = 1,NSYM
C
            ISYMK  = MULD2H(ISYMAI,ISYDMO)
            ISYMG  = MULD2H(ISYMK,ISYMLP)
C
            KOFF1  = IT2BCD(ISYMAI,ISYMK) + 1
            KOFF2  = IGLMRH(ISYMG,ISYMK) + 1
            KOFF3  = ID2AIG(ISYMAI,ISYMG) + 1
C
            NTOTAI = MAX(NT1AM(ISYMAI),1)
            NTOTG  = MAX(NBAS(ISYMG),1)
C
            CALL DGEMM('N','T',NT1AM(ISYMAI),NBAS(ISYMG),NRHF(ISYMK),
     *                 ONE,D2PQMO(KOFF1),NTOTAI,XLAMDP(KOFF2),NTOTG,
     *                 ONE,D2PQG(KOFF3),NTOTAI)
C
  110    CONTINUE
C
C-----------------------------------------------------------------
C     Calculate the transformation of the (vir.vir,occ;del) block.
C-----------------------------------------------------------------
C
      ELSE IF (ICON .EQ. 3) THEN
C
         DO 120 ISYMAB = 1,NSYM
C
            ISYMK  = MULD2H(ISYMAB,ISYDMO)
            ISYMG  = MULD2H(ISYMK,ISYMLP)
C
            KOFF1  = ICKASR(ISYMAB,ISYMK) + 1
            KOFF2  = IGLMRH(ISYMG,ISYMK) + 1
            KOFF3  = ID2ABG(ISYMAB,ISYMG) + 1
C
            NTOTAB = MAX(NMATAB(ISYMAB),1)
            NTOTG  = MAX(NBAS(ISYMG),1)
C
            CALL DGEMM('N','T',NMATAB(ISYMAB),NBAS(ISYMG),NRHF(ISYMK),
     *                 ONE,D2PQMO(KOFF1),NTOTAB,XLAMDP(KOFF2),NTOTG,
     *                 ONE,D2PQG(KOFF3),NTOTAB)
C
  120    CONTINUE
C
C-----------------------------------------------------------------
C     Calculate the transformation of the (occ.occ,vir;del) block.
C-----------------------------------------------------------------
C
      ELSE IF (ICON .EQ. 4) THEN
C
         DO 130 ISYMIJ = 1,NSYM
C
            ISYMA  = MULD2H(ISYMIJ,ISYDMO)
            ISYMG  = MULD2H(ISYMA,ISYMLP)
C
            KOFF1  = IMAIJA(ISYMIJ,ISYMA) + 1
            KOFF2  = IGLMVI(ISYMG,ISYMA) + 1
            KOFF3  = ID2IJG(ISYMIJ,ISYMG) + 1
C
            NTOTIJ = MAX(NMATIJ(ISYMIJ),1)
            NTOTG  = MAX(NBAS(ISYMG),1)
C
            CALL DGEMM('N','T',NMATIJ(ISYMIJ),NBAS(ISYMG),NVIR(ISYMA),
     *                 ONE,D2PQMO(KOFF1),NTOTIJ,XLAMDP(KOFF2),NTOTG,
     *                 ONE,D2PQG(KOFF3),NTOTIJ)
C
  130    CONTINUE
C
C-----------------------------------------------------------------
C     Calculate the transformation of the (occ.vir,vir;del) block.
C-----------------------------------------------------------------
C
      ELSE IF (ICON .EQ. 5) THEN
C
         DO 140 ISYMAI = 1,NSYM
C
            ISYMB  = MULD2H(ISYMAI,ISYDMO)
            ISYMG  = MULD2H(ISYMB,ISYMLP)
C
            KOFF1  = ICKATR(ISYMAI,ISYMB) + 1
            KOFF2  = IGLMVI(ISYMG,ISYMB) + 1
            KOFF3  = ID2AIG(ISYMAI,ISYMG) + 1
C
            NTOTAI = MAX(NT1AM(ISYMAI),1)
            NTOTG  = MAX(NBAS(ISYMG),1)
C
            CALL DGEMM('N','T',NT1AM(ISYMAI),NBAS(ISYMG),NVIR(ISYMB),
     *                 ONE,D2PQMO(KOFF1),NTOTAI,XLAMDP(KOFF2),NTOTG,
     *                 ONE,D2PQG(KOFF3),NTOTAI)
C
  140    CONTINUE
C
      ENDIF
C
      CALL QEXIT('CCD2_PQAO')
C
      RETURN
      END
C  /* Deck diac11 */
      SUBROUTINE DIAC11(D2IAB,ISYIAB,DAI,ISYD1,T2TIN,ISYMTI)
C
C     Written by Asger Halkier 20/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the one and only contribution to D2IAB
C              from projection against singles!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION D2IAB(*), DAI(*), T2TIN(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('DIAC11')
C
      IF (MULD2H(ISYMTI,ISYD1).NE.ISYIAB)
     &  CALL QUIT('Symmetry mismatch in DIAC11')
C
C----------------------------
C     Calculate contribution.
C----------------------------
C
      DO 100 ISYMAI = 1,NSYM
C
         ISYMM  = MULD2H(ISYMAI,ISYMTI)
         ISYMC  = MULD2H(ISYMM,ISYD1)
C
         KOFF1  = IT2BCD(ISYMAI,ISYMM) + 1
         KOFF2  = IT1AM(ISYMC,ISYMM) + 1
         KOFF3  = ICKATR(ISYMAI,ISYMC) + 1
C
         NTOTAI = MAX(NT1AM(ISYMAI),1)
         NTOTC  = MAX(NVIR(ISYMC),1)
C
         CALL DGEMM('N','T',NT1AM(ISYMAI),NVIR(ISYMC),NRHF(ISYMM),ONE,
     *              T2TIN(KOFF1),NTOTAI,DAI(KOFF2),NTOTC,ONE,
     *              D2IAB(KOFF3),NTOTAI)
C
  100 CONTINUE
C
      CALL QEXIT('DIAC11')
C
      RETURN
      END
C  /* Deck dija11 */
      SUBROUTINE DIJA11(D2IJA,ISYDEN,DAI,ISYDAI,XLAMDH,ISYMLH,
     &                  WORK,LWORK,IDEL,ISYMD)
C
C     Written by Asger Halkier 20/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the first and second contribution to D2IJA
C              from projection against singles!
C
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>

#if defined (SYS_CRAY)
      REAL ZERO, ONE, TWO
#else
      DOUBLE PRECISION ZERO, ONE, TWO
#endif
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)

      INTEGER ISYDEN, ISYDAI, ISYMLH, LWORK, IDEL, ISYMD
      INTEGER ISYMK, ISYMA, ISYMIJ, ISYMI, ISYMJ, ISYMAJ,
     &        KOFF1, KOFF2, NTOTA, NIJA, NDEI, NAJ

#if defined (SYS_CRAY)
      REAL D2IJA(*), DAI(*), XLAMDH(*), WORK(LWORK)
#else
      DOUBLE PRECISION D2IJA(*), DAI(*), XLAMDH(*), WORK(LWORK)
#endif

      CALL QENTER('DIJA11')
C
C
      IF (ISYDAI.NE.1) CALL QUIT('Illegal ISYDAI in DIJA11')
C
      ISYMK  = MULD2H(ISYMLH,ISYMD)
      ISYMA  = MULD2H(ISYDAI,ISYMK)
C
      IF (ISYMA.NE.ISYDEN) THEN
        CALL QUIT('Symmetry mismatch in DIJA11. (2)')
      END IF
C
C--------------------------------------------------------
C     Calculate contraction of D(ai) & lamda hole matrix.
C     (back transform occupied index of D(ai) with lambda:
c         d(a d) = sum_k D(ak) XL(dk)
C--------------------------------------------------------
C
      IF (LWORK .LT. NVIR(ISYMA)) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', NVIR(ISYMA)
         CALL QUIT('Insufficient space for allocation in DIJA11')
      ENDIF
C
      CALL DZERO(WORK,NVIR(ISYMA))
C
      KOFF1 = IT1AM(ISYMA,ISYMK) + 1
      KOFF2 = IGLMRH(ISYMD,ISYMK) + IDEL - IBAS(ISYMD)
C
      NTOTA = MAX(NVIR(ISYMA),1)
C
      CALL DGEMV('N',NVIR(ISYMA),NRHF(ISYMK),ONE,DAI(KOFF1),NTOTA,
     *           XLAMDH(KOFF2),NBAS(ISYMD),ZERO,WORK,1)
C
C--------------------------------------
C     Calculate the first contribution.
C--------------------------------------
C
      ISYMIJ = 1
C
      DO 100 A = 1,NVIR(ISYMA)
C
         DO 110 ISYMJ = 1,NSYM
C
            ISYMI = ISYMJ
C
            DO 120 J = 1,NRHF(ISYMJ)
C
               I    = J
               NIJA = IMAIJA(ISYMIJ,ISYMA) + NMATIJ(ISYMIJ)*(A - 1)
     *              + IMATIJ(ISYMI,ISYMJ)  + NRHF(ISYMI)*(J - 1) + I
C
               D2IJA(NIJA) = D2IJA(NIJA) + TWO*WORK(A)
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
C---------------------------------------
C     Calculate the second contribution.
C---------------------------------------
C
      ISYMAJ = ISYDAI
      ISYMI  = MULD2H(ISYMLH,ISYMD)
C
      IF (ISYMI.NE.ISYDEN) THEN
        CALL QUIT('Symmetry mismatch in DIJA11. (1)')
      END IF
C
      DO 130 I = 1,NRHF(ISYMI)
C
         NDEI = IGLMRH(ISYMD,ISYMI) + NBAS(ISYMD)*(I - 1)
     *        + IDEL - IBAS(ISYMD)
C
         DO 140 ISYMA = 1,NSYM
C
            ISYMJ  = MULD2H(ISYMA,ISYMAJ)
            ISYMIJ = MULD2H(ISYMI,ISYMJ)
C
            DO 150 J = 1,NRHF(ISYMJ)
C
               DO 160 A = 1,NVIR(ISYMA)
C
                  NAJ  = IT1AM(ISYMA,ISYMJ) + NVIR(ISYMA)*(J - 1) + A
                  NIJA = IMAIJA(ISYMIJ,ISYMA) + NMATIJ(ISYMIJ)*(A - 1)
     *              + IMATIJ(ISYMI,ISYMJ)  + NRHF(ISYMI)*(J - 1) + I
C
                  D2IJA(NIJA) = D2IJA(NIJA) - DAI(NAJ)*XLAMDH(NDEI)
C
  160          CONTINUE
  150       CONTINUE
  140    CONTINUE
  130 CONTINUE
C
      CALL QEXIT('DIJA11')
C
      RETURN
      END
C  /* Deck diac21 */
      SUBROUTINE DIAC21(D2IAC,ISYIAC,YMAT,ISYMY,XLAMDH,ISYMLH,
     *                  IDEL,ISYMD)
C
C     Written by Asger Halkier 21/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the first contribution to D2IAC
C              from projection against doubles!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION D2IAC(*), YMAT(*), XLAMDH(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('DIAC21')
C
      ISYMI  = MULD2H(ISYMD,ISYMLH)
      ISYMAC = ISYMY
C
      IF (MULD2H(ISYMI,ISYMAC).NE.ISYIAC)
     *  CALL QUIT('Symmetry mismatch in DIAC21')
C
C--------------------------------
C     Calculate the contribution.
C--------------------------------
C
      DO 100 I = 1,NRHF(ISYMI)
C
         NDEI = IGLMRH(ISYMD,ISYMI) + NBAS(ISYMD)*(I - 1) 
     *        + IDEL - IBAS(ISYMD)
C
         FACT = -XLAMDH(NDEI)
C
         DO 110 ISYMA = 1,NSYM
C
            ISYMC  = MULD2H(ISYMA,ISYMAC)
            ISYMAI = MULD2H(ISYMA,ISYMI)
C
            DO 120 C = 1,NVIR(ISYMC)
C
               NAC  = IMATAB(ISYMA,ISYMC)  + NVIR(ISYMA)*(C - 1) + 1
               NAIC = ICKATR(ISYMAI,ISYMC) + NT1AM(ISYMAI)*(C - 1)
     *              + IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + 1
C
               CALL DAXPY(NVIR(ISYMA),FACT,YMAT(NAC),1,
     *                    D2IAC(NAIC),1)
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('DIAC21')
C
      RETURN
      END
C  /* Deck dija21 */
      SUBROUTINE DIJA21(D2IJA,ISYDEN,DAB,ISYDAB,XLAMDH,ISYMLH,
     &                  WORK,LWORK,IDEL,ISYMD)
C
C     Written by Asger Halkier 21/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the first contribution to D2IJA
C              from projection against doubles!
C
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
#if defined (SYS_CRAY)
      REAL ZERO, ONE, TWO
#else
      DOUBLE PRECISION ZERO, ONE, TWO
#endif
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)

      INTEGER ISYDEN, ISYDAB, ISYMLH, IDEL, ISYMD, LWORK
      INTEGER ISYMB, ISYMA, ISYMIJ, ISYMI, ISYMJ,
     &        KOFF1, KOFF2, NTOTA, NIJA

#if defined (SYS_CRAY)
      REAL D2IJA(*), DAB(*), XLAMDH(*), WORK(LWORK)
#else
      DOUBLE PRECISION D2IJA(*), DAB(*), XLAMDH(*), WORK(LWORK)
#endif
C
      CALL QENTER('DIJA21')
C
      IF (ISYDAB.NE.1) CALL QUIT('Illegal ISYDAI in DIJA21')
C
      ISYMB  = MULD2H(ISYMLH,ISYMD)
      ISYMA  = MULD2H(ISYDAB,ISYMB)
      ISYMIJ = 1
C
      IF (ISYMA.NE.ISYDEN) THEN
        CALL QUIT('Symmetry mismatch in DIJA21. (2)')
      END IF
C
C-------------------------------------------------------------------
C     Calculate contraction of transposed Y-matrix and lamda matrix.
C-------------------------------------------------------------------
C
      IF (LWORK .LT. NVIR(ISYMA)) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', NVIR(ISYMA)
         CALL QUIT('Insufficient space for allocation in DIJA21')
      ENDIF
C
      CALL DZERO(WORK,NVIR(ISYMA))
C
      KOFF1 = IMATAB(ISYMA,ISYMB) + 1
      KOFF2 = IGLMVI(ISYMD,ISYMB) + IDEL - IBAS(ISYMD)
C
      NTOTA = MAX(NVIR(ISYMA),1)
C
      CALL DGEMV('N',NVIR(ISYMA),NVIR(ISYMB),ONE,DAB(KOFF1),NTOTA,
     *           XLAMDH(KOFF2),NBAS(ISYMD),ZERO,WORK,1)
C
      DO 100 A = 1,NVIR(ISYMA)
C
         DO 110 ISYMI = 1,NSYM
C
            ISYMJ = ISYMI
C
            DO 120 I = 1,NRHF(ISYMI)
C
            J    = I
            NIJA = IMAIJA(ISYMIJ,ISYMA) + NMATIJ(ISYMIJ)*(A - 1)
     *           + IMATIJ(ISYMI,ISYMJ)  + NRHF(ISYMI)*(J - 1) + I
C
            D2IJA(NIJA) = D2IJA(NIJA)   + TWO*WORK(A)
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('DIJA21')
C
      RETURN
      END
C  /* Deck dija22 */
      SUBROUTINE DIJA22(D2IJA,ISYDEN,Z2AM,T2INT,ISYMTI,WORK,LWORK)
C
C     Written by Asger Halkier 21/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the second contribution to D2IJA
C              from projection against doubles!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION D2IJA(*), Z2AM(*), T2INT(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('DIJA22')
C
      IF (ISYDEN.NE.ISYMTI) CALL QUIT('Symmetry mismatch in DIJA22')
C
      ISYCKI = ISYDEN
      ISYIJA = ISYDEN
C
      DO 100 ISYMI = 1,NSYM
C
         ISYMCK = MULD2H(ISYMI,ISYCKI)
         ISYMAJ = ISYMCK
C
C------------------------------
C        Work space allocation.
C------------------------------
C
         KSCR = 1
         KEND = KSCR  + NT1AM(ISYMAJ)*NRHF(ISYMI)
         LWRK = LWORK - KEND
C
         IF (LWRK .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND
            CALL QUIT('Insufficient work space for '//
     &                'allocation in DIJA22')
         ENDIF
C
         CALL DZERO(WORK(KSCR),NT1AM(ISYMAJ)*NRHF(ISYMI))
C
C--------------------------------------------
C        Calculate contraction of zeta and t.
C--------------------------------------------
C
         KOFF1  = IT2SQ(ISYMAJ,ISYMCK) + 1
         KOFF2  = IT2BCD(ISYMCK,ISYMI) + 1
C
         NTOTAJ = MAX(NT1AM(ISYMAJ),1)
         NTOTCK = MAX(NT1AM(ISYMCK),1)
C
         CALL DGEMM('N','N',NT1AM(ISYMAJ),NRHF(ISYMI),NT1AM(ISYMCK),
     *              ONE,Z2AM(KOFF1),NTOTAJ,T2INT(KOFF2),NTOTCK,ZERO,
     *              WORK(KSCR),NTOTAJ)
C
C--------------------------------------------------------- 
C        Store properly reordered scratch-array in result.
C--------------------------------------------------------- 
C
         DO 110 ISYMA = 1,NSYM
C
            ISYMJ  = MULD2H(ISYMA,ISYMAJ)
            ISYMIJ = MULD2H(ISYMJ,ISYMI)
C
            DO 120 I = 1,NRHF(ISYMI)
C
               DO 130 J = 1,NRHF(ISYMJ)
C
                  DO 140 A = 1,NVIR(ISYMA)
C
                     NAJI = NT1AM(ISYMAJ)*(I - 1) + IT1AM(ISYMA,ISYMJ)
     *                    + NVIR(ISYMA)*(J - 1)   + A
                     NIJA = IMAIJA(ISYMIJ,ISYMA) 
     *                    + NMATIJ(ISYMIJ)*(A - 1) 
     *                    + IMATIJ(ISYMI,ISYMJ) 
     *                    + NRHF(ISYMI)*(J - 1) + I
C
                     D2IJA(NIJA) = D2IJA(NIJA) - WORK(KSCR + NAJI - 1)
C
  140             CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('DIJA22')
C
      RETURN
      END
C  /* Deck dija23 */
      SUBROUTINE DIJA23(D2IJA,ISYDEN,Z2AM,T2INT,ISYMTI,WORK,LWORK)
C
C     Written by Asger Halkier 21/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the third contribution to D2IJA
C              from projection against doubles!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION D2IJA(*), Z2AM(*), T2INT(*), WORK(LWORK)
#include <priunit.h>      
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('DIJA23')
C
      IF (ISYDEN.NE.ISYMTI) CALL QUIT('Symmetry mismatch in DIJA23')
C
      ISYCIK = ISYDEN
      ISYIJA = ISYDEN
C
      DO 100 ISYMI = 1,NSYM
C
         ISYMCK = MULD2H(ISYMI,ISYCIK)
         ISYMAJ = ISYMCK
C
C----------------------------------
C        Work space allocation one.
C----------------------------------
C
         KSCR1 = 1
         KEND1 = KSCR1 + NRHF(ISYMI)*NT1AM(ISYMCK)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT(
     *           'Insufficient work space for allocation in DIJA23')
         ENDIF
C
         CALL DZERO(WORK(KSCR1),NRHF(ISYMI)*NT1AM(ISYMCK))
C
C-------------------------------------------------------------
C        Reorder backtransformed t-amplitude in scratch-array.
C-------------------------------------------------------------
C
         DO 110 ISYMC = 1,NSYM
C
            ISYMK  = MULD2H(ISYMC,ISYMCK)
            ISYMCI = MULD2H(ISYMC,ISYMI)
C
            DO 120 K = 1,NRHF(ISYMK)
C
               DO 130 I = 1,NRHF(ISYMI)
C
                  DO 140 C = 1,NVIR(ISYMC)
C
                     NCIK = IT2BCD(ISYMCI,ISYMK) 
     *                    + NT1AM(ISYMCI)*(K - 1) + IT1AM(ISYMC,ISYMI)
     *                    + NVIR(ISYMC)*(I - 1) + C
                     NCK  = IT1AM(ISYMC,ISYMK) 
     *                    + NVIR(ISYMC)*(K - 1) + C
                     NICK = KSCR1 + NRHF(ISYMI)*(NCK - 1) + I - 1
C
                     WORK(NICK) = T2INT(NCIK)
C
  140             CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
C
         DO 150 ISYMA = 1,NSYM
C
            ISYMJ  = MULD2H(ISYMA,ISYMAJ)
            ISYMIJ = MULD2H(ISYMA,ISYDEN) 
C
C-------------------------------------
C           Work space allocation two.
C-------------------------------------
C
            KSCR2 = KEND1
            KEND2 = KSCR2 + NT1AM(ISYMCK)*NRHF(ISYMJ)
            LWRK2 = LWORK - KEND2
C
            IF (LWRK2 .LT. 0) THEN
               WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
               CALL QUIT('Insufficient work space for '//
     &                   'allocation in DIJA23')
            ENDIF
C
            DO 160 A = 1,NVIR(ISYMA)
C
               CALL DZERO(WORK(KSCR2),NT1AM(ISYMCK)*NRHF(ISYMJ))
C
C---------------------------------------------------
C              Copy submatrix out of zeta-amplitude.
C---------------------------------------------------
C
               DO 170 ISYMC = 1,NSYM
C
                  ISYMCJ = MULD2H(ISYMC,ISYMJ)
                  ISYMAK = ISYMCJ
                  ISYMK  = MULD2H(ISYMC,ISYMCK)
C
                  DO 180 K = 1,NRHF(ISYMK)
C
                     NAK = IT1AM(ISYMA,ISYMK) + NVIR(ISYMA)*(K - 1) + A
C
                     DO 190 J = 1,NRHF(ISYMJ)
C
                        NCJAK = IT2SQ(ISYMCJ,ISYMAK)
     *                        + NT1AM(ISYMCJ)*(NAK - 1)
     *                        + IT1AM(ISYMC,ISYMJ)
     *                        + NVIR(ISYMC)*(J - 1) + 1
                        NCKJ  = KSCR2 + NT1AM(ISYMCK)*(J - 1)
     *                        + IT1AM(ISYMC,ISYMK) 
     *                        + NVIR(ISYMC)*(K - 1) 
C
                        CALL DCOPY(NVIR(ISYMC),Z2AM(NCJAK),1,
     *                             WORK(NCKJ),1)
C
  190                CONTINUE
  180             CONTINUE
  170          CONTINUE
C
C------------------------------------------------------
C              Final contraction and storage in result.
C------------------------------------------------------
C
               KOFF1  = IMAIJA(ISYMIJ,ISYMA) + NMATIJ(ISYMIJ)*(A - 1)
     *                + IMATIJ(ISYMI,ISYMJ) + 1
C
               NTOTI  = MAX(NRHF(ISYMI),1)
               NTOTCK = MAX(NT1AM(ISYMCK),1)
C
               CALL DGEMM('N','N',NRHF(ISYMI),NRHF(ISYMJ),
     *                    NT1AM(ISYMCK),-ONE,WORK(KSCR1),NTOTI,
     *                    WORK(KSCR2),NTOTCK,ONE,D2IJA(KOFF1),NTOTI)
C
  160       CONTINUE
  150    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('DIJA23')
C
      RETURN
      END
C  /* Deck daib21 */
      SUBROUTINE DAIB21(D2AIB,ISYAIB,Z2AM,XLAMDH,ISYMLH,WORK,LWORK,
     *                  IDEL,ISYMD)
C
C     Written by Asger Halkier 22/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the one and only contribution to D2AIB
C              at all which arises from projection against doubles!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION D2AIB(*), Z2AM(*), XLAMDH(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('DAIB21')
C
      ISYMJ  = MULD2H(ISYMD,ISYMLH)
C
      IF (MULD2H(ISYMLH,ISYMD).NE.ISYAIB) 
     &  CALL QUIT('Symmetry mismatch in DAIB21')
C
C-------------------------------
C     Work space allocation one.
C-------------------------------
C
      KJVEC = 1
      KEND1 = KJVEC + NRHF(ISYMJ)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space for allocation in DAIB21')
      ENDIF
C
      CALL DZERO(WORK(KJVEC),NRHF(ISYMJ))
C
C------------------------------------------
C     Copy vector out of lamda hole matrix.
C------------------------------------------
C
      KOFF1 = IGLMRH(ISYMD,ISYMJ) + IDEL - IBAS(ISYMD)
C
      CALL DCOPY(NRHF(ISYMJ),XLAMDH(KOFF1),NBAS(ISYMD),WORK(KJVEC),1)
C
      DO 100 ISYMB = 1,NSYM
C
         ISYMAI = MULD2H(ISYMB,ISYAIB)
         ISYMBJ = MULD2H(ISYMB,ISYMJ)
C
C----------------------------------
C        Work space allocation two.
C----------------------------------
C
         KZSUB = KEND1
         KEND2 = KZSUB + NT1AM(ISYMAI)*NRHF(ISYMJ)
         LWRK2 = LWORK - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
            CALL QUIT(
     *           'Insufficient work space for allocation in DAIB21')
         ENDIF
C
         CALL DZERO(WORK(KZSUB),NT1AM(ISYMAI)*NRHF(ISYMJ))
C
         DO 110 B = 1,NVIR(ISYMB)
C
C------------------------------------------------
C           Copy submatrix out of zeta-amplitude.
C------------------------------------------------
C
            DO 120 J = 1,NRHF(ISYMJ)
C
               NBJ   = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
               NAIBJ = IT2SQ(ISYMAI,ISYMBJ) 
     *               + NT1AM(ISYMAI)*(NBJ - 1) + 1
               NAIJ  = KZSUB + NT1AM(ISYMAI)*(J - 1)
C
               CALL DCOPY(NT1AM(ISYMAI),Z2AM(NAIBJ),1,WORK(NAIJ),1)
C
  120       CONTINUE
C
C---------------------------------------------------
C           Final contraction and storage in result.
C---------------------------------------------------
C
            KOFF2  = ICKATR(ISYMAI,ISYMB) + NT1AM(ISYMAI)*(B - 1) + 1
C
            NTOTAI = MAX(NT1AM(ISYMAI),1)
C
            CALL DGEMV('N',NT1AM(ISYMAI),NRHF(ISYMJ),ONE,WORK(KZSUB),
     *                 NTOTAI,WORK(KJVEC),1,ONE,D2AIB(KOFF2),1)
C
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('DAIB21')
C
      RETURN
      END
C  /* Deck diac22 */
      SUBROUTINE DIAC22(D2IAB,ISYIAC,T2AMT,Z2AO,ISYZ2AO,WORK,LWORK)
C
C     Written by Asger Halkier 24/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the second contribution to D2IAC
C              from projection against doubles!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION D2IAB(*), T2AMT(*), Z2AO(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('DIAC22')
C
      DO 100 ISYMI = 1,NSYM
C
         ISYMAC = MULD2H(ISYMI,ISYIAC)
C
         DO 110 I = 1,NRHF(ISYMI)
C
            DO 120 ISYMA = 1,NSYM
C
               ISYMC  = MULD2H(ISYMA,ISYMAC)
               ISYMDK = MULD2H(ISYMC,ISYZ2AO)
               ISYMAI = MULD2H(ISYMA,ISYMI)
C
C------------------------------------
C              Work space allocation.
C------------------------------------
C
               KT2SUB = 1
               KACRES = KT2SUB + NVIR(ISYMA)*NT1AM(ISYMDK)
               KEND1  = KACRES + NVIR(ISYMA)*NVIR(ISYMC)
               LWRK1  = LWORK  - KEND1
C
               IF (LWRK1 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
                  CALL QUIT('Insufficient work space in routine DIAC22')
               ENDIF
C
               CALL DZERO(WORK(KT2SUB),NVIR(ISYMA)*NT1AM(ISYMDK))
               CALL DZERO(WORK(KACRES),NVIR(ISYMA)*NVIR(ISYMC))
C
C------------------------------------------
C              Copy submatrix out of t2amt.
C------------------------------------------
C
               DO 130 NDK = 1,NT1AM(ISYMDK)
C
                  DO 140 A = 1,NVIR(ISYMA)
C
                     NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A
                     NAIDK = IT2AM(ISYMAI,ISYMDK) + INDEX(NAI,NDK)
                     NADK  = KT2SUB + NVIR(ISYMA)*(NDK - 1) + A - 1
C
                     WORK(NADK) = T2AMT(NAIDK)
C
  140             CONTINUE
  130          CONTINUE
C
C-----------------------------------------------------------
C              Calculate contraction of t2amt(sub) and z2ao.
C-----------------------------------------------------------
C
               KOFF1  = ICKATR(ISYMDK,ISYMC) + 1
C
               NTOTA  = MAX(NVIR(ISYMA),1)
               NTOTDK = MAX(NT1AM(ISYMDK),1)
C
               CALL DGEMM('N','N',NVIR(ISYMA),NVIR(ISYMC),
     *                    NT1AM(ISYMDK),ONE,WORK(KT2SUB),NTOTA,
     *                    Z2AO(KOFF1),NTOTDK,ZERO,WORK(KACRES),NTOTA)
C
C--------------------------------------
C              Final storage in result.
C--------------------------------------
C
               DO 150 C = 1,NVIR(ISYMC)
C
                  NAC  = KACRES + NVIR(ISYMA)*(C - 1)
                  NAIC = ICKATR(ISYMAI,ISYMC) + NT1AM(ISYMAI)*(C - 1)
     *                 + IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + 1
C
                  CALL DAXPY(NVIR(ISYMA),ONE,WORK(NAC),1,D2IAB(NAIC),1)
C
  150          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('DIAC22')
C
      RETURN
      END
C  /* Deck dabi22 */
      SUBROUTINE DABI22(D2ABI,ISYABI,T2AM,Z2AO,ISYZ2AO,WORK,LWORK)
C
C     Written by Asger Halkier 24/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the second contribution to D2AIB
C              from projection against doubles!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION D2ABI(*), T2AM(*), Z2AO(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('DABI22')
C
      DO 100 ISYMI = 1,NSYM
C
         ISYMAB = MULD2H(ISYMI,ISYABI)
C
         DO 110 I = 1,NRHF(ISYMI)
C
            DO 120 ISYMA = 1,NSYM
C
               ISYMB  = MULD2H(ISYMA,ISYMAB)
               ISYMCK = MULD2H(ISYMA,ISYZ2AO)
               ISYMBI = MULD2H(ISYMA,ISYABI)
C
C------------------------------------
C              Work space allocation.
C------------------------------------
C
               KT2SUB = 1
               KBARES = KT2SUB + NVIR(ISYMB)*NT1AM(ISYMCK)
               KEND1  = KBARES + NVIR(ISYMB)*NVIR(ISYMA)
               LWRK1  = LWORK  - KEND1
C
               IF (LWRK1 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
                  CALL QUIT('Insufficient work space in routine DABI22')
               ENDIF
C
               CALL DZERO(WORK(KT2SUB),NVIR(ISYMB)*NT1AM(ISYMCK))
               CALL DZERO(WORK(KBARES),NVIR(ISYMB)*NVIR(ISYMA))
C
C------------------------------------------
C              Copy submatrix out of t2amt.
C------------------------------------------
C
               DO 130 NCK = 1,NT1AM(ISYMCK)
C
                  DO 140 B = 1,NVIR(ISYMB)
C
                     NBI = IT1AM(ISYMB,ISYMI) + NVIR(ISYMB)*(I - 1) + B
                     NBICK = IT2AM(ISYMBI,ISYMCK) + INDEX(NBI,NCK)
                     NBCK  = KT2SUB + NVIR(ISYMB)*(NCK - 1) + B - 1
C
                     WORK(NBCK) = T2AM(NBICK)
C
  140             CONTINUE
  130          CONTINUE
C
C----------------------------------------------------------
C              Calculate contraction of t2am(sub) and z2ao.
C----------------------------------------------------------
C
               KOFF1  = ICKATR(ISYMCK,ISYMA) + 1
C
               NTOTB  = MAX(NVIR(ISYMB),1)
               NTOTCK = MAX(NT1AM(ISYMCK),1)
C
               CALL DGEMM('N','N',NVIR(ISYMB),NVIR(ISYMA),
     *                    NT1AM(ISYMCK),ONE,WORK(KT2SUB),NTOTB,
     *                    Z2AO(KOFF1),NTOTCK,ZERO,WORK(KBARES),NTOTB)
C
C--------------------------------------
C              Final storage in result.
C--------------------------------------
C
               DO 150 A = 1,NVIR(ISYMA)
C
                  NBA = KBARES + NVIR(ISYMB)*(A - 1)
                  NAB = ICKASR(ISYMAB,ISYMI) + NMATAB(ISYMAB)*(I - 1)
     *                + IMATAB(ISYMA,ISYMB) + A
C
                  CALL DAXPY(NVIR(ISYMB),-ONE,WORK(NBA),1,D2ABI(NAB),
     *                       NVIR(ISYMA))
C
  150          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('DABI22')
C
      RETURN
      END
C  /* Deck dabi23 */
      SUBROUTINE DABI23(D2ABI,ISYABI,T2AM,Z2AO,ISYZ2AO,WORK,LWORK)
C
C     Written by Asger Halkier 24/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the third contribution to D2AIB
C              from projection against doubles!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION D2ABI(*), T2AM(*), Z2AO(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('DABI23')
C
      DO 100 ISYMA = 1,NSYM
C
         ISYMCK = MULD2H(ISYMA,ISYZ2AO)
         ISYMBI = ISYMCK
C
C----------------------------------
C        Work space allocation one.
C----------------------------------
C
         KZREOR = 1
         KEND1  = KZREOR + NVIR(ISYMA)*NT1AM(ISYMCK)
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT('Insufficient work space in routine DABI23')
         ENDIF
C
         CALL DZERO(WORK(KZREOR),NVIR(ISYMA)*NT1AM(ISYMCK))
C
C-------------------------------------------------
C        Reorder z2ao(ak,c;del) to z2ao(a,ck;del).
C-------------------------------------------------
C
         DO 110 ISYMC = 1,NSYM
C
            ISYMK  = MULD2H(ISYMC,ISYMCK)
            ISYMAK = MULD2H(ISYMC,ISYZ2AO)
C
            DO 120 C = 1,NVIR(ISYMC)
C
               DO 130 K = 1,NRHF(ISYMK)
C
                  NAKC = ICKATR(ISYMAK,ISYMC) + NT1AM(ISYMAK)*(C - 1)
     *                 + IT1AM(ISYMA,ISYMK) + NVIR(ISYMA)*(K - 1) + 1
                  NCK  = IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K - 1) + C
                  NACK = KZREOR + NVIR(ISYMA)*(NCK - 1)
C
                  CALL DCOPY(NVIR(ISYMA),Z2AO(NAKC),1,WORK(NACK),1)
C
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
C
         DO 140 ISYMI = 1,NSYM
C
            ISYMB  = MULD2H(ISYMI,ISYMBI)
            ISYMAB = MULD2H(ISYMI,ISYABI)
C
C-------------------------------------
C           Work space allocation two.
C-------------------------------------
C
            KT2SUB = KEND1
            KEND2  = KT2SUB + NT1AM(ISYMCK)*NVIR(ISYMB)
            LWRK2  = LWORK  - KEND2
C
            IF (LWRK2 .LT. 0) THEN
               WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
               CALL QUIT('Insufficient work space in routine DABI23')
            ENDIF
C
            DO 150 I = 1,NRHF(ISYMI)
C
               CALL DZERO(WORK(KT2SUB),NT1AM(ISYMCK)*NVIR(ISYMB))
C
C-----------------------------------------
C              Copy submatrix out of t2am.
C-----------------------------------------
C
               DO 160 ISYMK = 1,NSYM
C
                  ISYMC  = MULD2H(ISYMK,ISYMCK)
                  ISYMCI = MULD2H(ISYMC,ISYMI)
                  ISYMBK = ISYMCI
C
                  DO 170 K = 1,NRHF(ISYMK)
C
                     DO 180 B = 1,NVIR(ISYMB)
C
                        NBK = IT1AM(ISYMB,ISYMK) + NVIR(ISYMB)*(K - 1)
     *                      + B
C
                        DO 190 C = 1,NVIR(ISYMC)
C
                           NCI   = IT1AM(ISYMC,ISYMI) 
     *                           + NVIR(ISYMC)*(I - 1) + C
                           NCIBK = IT2AM(ISYMCI,ISYMBK) 
     *                           + INDEX(NCI,NBK)
                           NCKB  = KT2SUB + NT1AM(ISYMCK)*(B - 1)
     *                           + IT1AM(ISYMC,ISYMK) 
     *                           + NVIR(ISYMC)*(K - 1) + C - 1
C
                           WORK(NCKB) = T2AM(NCIBK) 
C
  190                   CONTINUE
  180                CONTINUE
  170             CONTINUE
  160          CONTINUE
C
C------------------------------------------------------
C              Final contraction and storage in result.
C------------------------------------------------------
C
               KOFF1  = ICKASR(ISYMAB,ISYMI) + NMATAB(ISYMAB)*(I - 1) 
     *                + IMATAB(ISYMA,ISYMB)  + 1
C
               NTOTA  = MAX(NVIR(ISYMA),1)
               NTOTCK = MAX(NT1AM(ISYMCK),1)
C
               CALL DGEMM('N','N',NVIR(ISYMA),NVIR(ISYMB),
     *                    NT1AM(ISYMCK),-ONE,WORK(KZREOR),NTOTA,
     *                    WORK(KT2SUB),NTOTCK,ONE,D2ABI(KOFF1),NTOTA)
C
  150       CONTINUE
  140    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('DABI23')
C
      RETURN
      END
C  /* Deck diaj26 */
      SUBROUTINE DIAJ26(D2IAJ,ISYIAJ,T2AM,ZINT,ISYMZI,WORK,LWORK)
C
C     Written by Asger Halkier 25/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the sixth contribution to D2IAJ
C              from projection against doubles!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION D2IAJ(*), T2AM(*), ZINT(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('DIAJ26')
C
      DO 100 ISYMJ = 1,NSYM
C
         ISYMAI = MULD2H(ISYMJ,ISYIAJ)
C
         DO 110 J = 1,NRHF(ISYMJ)
C
            DO 120 ISYMA = 1,NSYM
C
               ISYMI  = MULD2H(ISYMA,ISYMAI)
               ISYMEM = MULD2H(ISYMI,ISYMZI)
               ISYMAJ = ISYMEM
C
C------------------------------------
C              Work space allocation.
C------------------------------------
C
               KT2SUB = 1
               KEND1  = KT2SUB + NVIR(ISYMA)*NT1AM(ISYMEM)
               LWRK1  = LWORK  - KEND1
C
               IF (LWRK1 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
                  CALL QUIT('Insufficient work space in routine DIAJ26')
               ENDIF
C
               CALL DZERO(WORK(KT2SUB),NVIR(ISYMA)*NT1AM(ISYMEM))
C
C-------------------------------------------------
C              Copy submatrix out of t2-amplitude.
C-------------------------------------------------
C
               DO 130 NEM = 1,NT1AM(ISYMEM)
C
                  DO 140 A = 1,NVIR(ISYMA)
C
                     NAJ = IT1AM(ISYMA,ISYMJ) + NVIR(ISYMA)*(J - 1) + A
                     NAJEM = IT2AM(ISYMAJ,ISYMEM) + INDEX(NAJ,NEM)
                     NAEM  = KT2SUB + NVIR(ISYMA)*(NEM - 1) + A - 1
C
                     WORK(NAEM) = T2AM(NAJEM)
C
  140             CONTINUE
  130          CONTINUE
C
C------------------------------------------------------
C              Final contraction and storage in result.
C------------------------------------------------------
C
               KOFF1  = IT2BCD(ISYMEM,ISYMI) + 1
               KOFF2  = IT2BCD(ISYMAI,ISYMJ) + NT1AM(ISYMAI)*(J - 1) 
     *                + IT1AM(ISYMA,ISYMI)   + 1
C
               NTOTA  = MAX(NVIR(ISYMA),1)
               NTOTEM = MAX(NT1AM(ISYMEM),1)
C
               CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),
     *                    NT1AM(ISYMEM),-ONE,WORK(KT2SUB),NTOTA,
     *                    ZINT(KOFF1),NTOTEM,ONE,D2IAJ(KOFF2),NTOTA)
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('DIAJ26')
C
      RETURN
      END
C  /* Deck diaj27 */
      SUBROUTINE DIAJ27(D2IAJ,ISYIAJ,T2AM,ZINT,ISYMZI,WORK,LWORK)
C
C     Written by Asger Halkier 25/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the seventh contribution to D2IAJ
C              from projection against doubles!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION D2IAJ(*), T2AM(*), ZINT(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('DIAJ27')
C
      IF (ISYIAJ.NE.ISYMZI) CALL QUIT('Symmetry mismatch in DIAJ27')
C
      DO 100 ISYMJ = 1,NSYM
C
         ISYMAI = MULD2H(ISYMJ,ISYIAJ)
C
         DO 110 J = 1,NRHF(ISYMJ)
C
            DO 120 ISYMA = 1,NSYM
C
               ISYMI  = MULD2H(ISYMA,ISYMAI)
               ISYMEM = MULD2H(ISYMI,ISYMZI)
               ISYMAJ = ISYMEM
C
C------------------------------------
C              Work space allocation.
C------------------------------------
C
               KT2SUB = 1
               KEND1  = KT2SUB + NVIR(ISYMA)*NT1AM(ISYMEM)
               LWRK1  = LWORK  - KEND1
C
               IF (LWRK1 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
                  CALL QUIT('Insufficient work space in routine DIAJ27')
               ENDIF
C
               CALL DZERO(WORK(KT2SUB),NVIR(ISYMA)*NT1AM(ISYMEM))
C
C-------------------------------------------------
C              Copy submatrix out of t2-amplitude.
C-------------------------------------------------
C
               DO 130 ISYME = 1,NSYM
C
                  ISYMM  = MULD2H(ISYME,ISYMEM)
                  ISYMAM = MULD2H(ISYMM,ISYMA)
                  ISYMEJ = ISYMAM
C
                  DO 140 E = 1,NVIR(ISYME)
C
                     NEJ = IT1AM(ISYME,ISYMJ) + NVIR(ISYME)*(J - 1) + E
C
                     DO 150 M = 1,NRHF(ISYMM)
C
                        NEM = IT1AM(ISYME,ISYMM)
     *                      + NVIR(ISYME)*(M - 1) + E
C
                        DO 160 A = 1,NVIR(ISYMA)
C
                           NAM   = IT1AM(ISYMA,ISYMM)
     *                           + NVIR(ISYMA)*(M - 1) + A
                           NAMEJ = IT2AM(ISYMAM,ISYMEJ) 
     *                           + INDEX(NAM,NEJ)
                           NAEM  = KT2SUB + NVIR(ISYMA)*(NEM - 1) 
     *                           + A - 1
C
                           WORK(NAEM) = T2AM(NAMEJ)
C
  160                   CONTINUE
  150                CONTINUE
  140             CONTINUE
  130          CONTINUE
C
C------------------------------------------------------
C              Final contraction and storage in result.
C------------------------------------------------------
C
               KOFF1  = IT2BCD(ISYMEM,ISYMI) + 1
               KOFF2  = IT2BCD(ISYMAI,ISYMJ) + NT1AM(ISYMAI)*(J - 1)
     *                + IT1AM(ISYMA,ISYMI)   + 1
C
               NTOTA  = MAX(NVIR(ISYMA),1)
               NTOTEM = MAX(NT1AM(ISYMEM),1)
C
               CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),
     *                    NT1AM(ISYMEM),ONE,WORK(KT2SUB),NTOTA,
     *                    ZINT(KOFF1),NTOTEM,ONE,D2IAJ(KOFF2),NTOTA)
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('DIAJ27')
C
      RETURN
      END
C  /* Deck diaj28 */
      SUBROUTINE DIAJ28(D2IAJ,ISYIAJ,T2AMT,ZINT,ISYMZI,WORK,LWORK)
C
C     Written by Asger Halkier 26/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the eighth contribution to D2IAJ
C              from projection against doubles!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION D2IAJ(*), T2AMT(*), ZINT(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('DIAJ28')
C
      DO 100 ISYMI = 1,NSYM
C
         ISYMAJ = MULD2H(ISYMI,ISYIAJ)
C
         DO 110 I = 1,NRHF(ISYMI)
C
            DO 120 ISYMA = 1,NSYM
C
               ISYMJ  = MULD2H(ISYMA,ISYMAJ)
               ISYMEM = MULD2H(ISYMJ,ISYMZI)
               ISYMAI = ISYMEM
C
C------------------------------------
C              Work space allocation.
C------------------------------------
C
               KT2SUB = 1
               KEND1  = KT2SUB + NVIR(ISYMA)*NT1AM(ISYMEM)
               LWRK1  = LWORK  - KEND1
C
               IF (LWRK1 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
                  CALL QUIT('Insufficient work space in routine DIAJ28')
               ENDIF
C
               CALL DZERO(WORK(KT2SUB),NVIR(ISYMA)*NT1AM(ISYMEM))
C
C------------------------------------------
C              Copy submatrix out of t2amt.
C------------------------------------------
C
               DO 130 NEM = 1,NT1AM(ISYMEM)
C
                  DO 140 A = 1,NVIR(ISYMA)
C
                     NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A
                     NAIEM = IT2AM(ISYMAI,ISYMEM) + INDEX(NAI,NEM)
                     NAEM  = KT2SUB + NVIR(ISYMA)*(NEM - 1) + A - 1
C
                     WORK(NAEM) = T2AMT(NAIEM)
C
  140             CONTINUE
  130          CONTINUE
C
C---------------------------------------------------------
C              Calculate final contraction in loop over j.
C---------------------------------------------------------
C
               DO 150 J = 1,NRHF(ISYMJ)
C
                  KOFF1 = IT2BCD(ISYMEM,ISYMJ) 
     *                  + NT1AM(ISYMEM)*(J - 1) + 1
                  KOFF2 = IT2BCD(ISYMAI,ISYMJ)
     *                  + NT1AM(ISYMAI)*(J - 1) 
     *                  + IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + 1
C
                  NTOTA = MAX(NVIR(ISYMA),1)
C
                  CALL DGEMV('N',NVIR(ISYMA),NT1AM(ISYMEM),ONE,
     *                       WORK(KT2SUB),NTOTA,ZINT(KOFF1),1,ONE,
     *                       D2IAJ(KOFF2),1)
C
  150          CONTINUE
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('DIAJ28')
C
      RETURN
      END
C  /* Deck cc_mirs */
      SUBROUTINE CC_MIRS(XMIRES,XMINT)
C
C     Written by Asger Halkier 26/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To resort the M-intermediate XMINT (im,j;k)
C              to (mk,i;j) and store the result in XMIRES! 
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION XMIRES(*), XMINT(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('CC_MIRS')
C
      CALL DZERO(XMIRES,N3ORHF(1))
C
C----------------------------------------------
C     Do the reordering through loads of loops.
C----------------------------------------------
C
      DO 100 ISYMK = 1,NSYM
C
         ISYIMJ = ISYMK
C
         DO 110 K = 1,NRHF(ISYMK)
C
            DO 120 ISYMJ = 1,NSYM
C
               ISYMIM = MULD2H(ISYMJ,ISYIMJ)
               ISYMKI = ISYMJ
C
               DO 130 J = 1,NRHF(ISYMJ)
C
                  DO 140 ISYMM = 1,NSYM
C
                     ISYMI  = MULD2H(ISYMM,ISYMIM)
                     ISYMMK = MULD2H(ISYMI,ISYMKI)
C
                     DO 150 M = 1,NRHF(ISYMM)
C
                        NIMJK = I3ORHF(ISYIMJ,ISYMK)
     *                        + NMAIJK(ISYIMJ)*(K - 1)
     *                        + IMAIJK(ISYMIM,ISYMJ)
     *                        + NMATIJ(ISYMIM)*(J - 1)
     *                        + IMATIJ(ISYMI,ISYMM) 
     *                        + NRHF(ISYMI)*(M - 1) + 1
                        NMK   = IMATIJ(ISYMM,ISYMK) 
     *                        + NRHF(ISYMM)*(K - 1) + M
                        NMKIJ = I3ORHF(ISYMKI,ISYMJ)
     *                        + NMAIJK(ISYMKI)*(J - 1)
     *                        + IMAIJK(ISYMMK,ISYMI) + NMK
C
                        CALL DCOPY(NRHF(ISYMI),XMINT(NIMJK),1,
     *                             XMIRES(NMKIJ),NMATIJ(ISYMMK))
C
  150                CONTINUE
  140             CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CC_MIRS')
C
      RETURN
      END
C  /* Deck diaj25 */
      SUBROUTINE DIAJ25(D2IAJ,ISYIAJ,T2INT,ISYMTI,XMIRES,WORK,LWORK)
C
C     Written by Asger Halkier 26/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the fifth contribution to D2IAJ
C              from projection against doubles!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION D2IAJ(*), T2INT(*), XMIRES(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('DIAJ25')
C
      IF (ISYMTI.NE.ISYIAJ) CALL QUIT('Symmetry mismatch in DIAJ25')
C
      DO 100 ISYMA = 1,NSYM
C
         ISYMMK = MULD2H(ISYMA,ISYMTI)
         ISYMIJ = MULD2H(ISYMA,ISYIAJ)
C
C------------------------------
C        Work space allocation.
C------------------------------
C
         KTREOR = 1
         KEND1  = KTREOR + NVIR(ISYMA)*NMATIJ(ISYMMK)
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT('Insufficient work space in routine DIAJ25')
         ENDIF
C
         CALL DZERO(WORK(KTREOR),NVIR(ISYMA)*NMATIJ(ISYMMK))
C
C--------------------------------------------
C        Reorder t2int in scratch-work-array.
C--------------------------------------------
C
         DO 110 ISYMM = 1,NSYM
C
            ISYMK  = MULD2H(ISYMM,ISYMMK)
            ISYMAM = MULD2H(ISYMK,ISYMTI)
C
            DO 120 K = 1,NRHF(ISYMK)
C
               DO 130 M = 1,NRHF(ISYMM)
C
                  NMK   = IMATIJ(ISYMM,ISYMK) + NRHF(ISYMM)*(K - 1) + M
                  NAMKN = IT2BCD(ISYMAM,ISYMK) + NT1AM(ISYMAM)*(K - 1)
     *                  + IT1AM(ISYMA,ISYMM) + NVIR(ISYMA)*(M - 1) + 1
                  NAMKR = KTREOR + NVIR(ISYMA)*(NMK - 1)
C
                  CALL DCOPY(NVIR(ISYMA),T2INT(NAMKN),1,WORK(NAMKR),1)
C
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
C
C---------------------------------------------------------
C        Calculate contraction with xmires in loop over j.
C---------------------------------------------------------
C
         DO 140 ISYMJ = 1,NSYM
C
            ISYMI  = MULD2H(ISYMJ,ISYMIJ)
            ISYMAI = MULD2H(ISYMJ,ISYIAJ)
            ISYMKI = ISYMJ
C
            DO 150 J = 1,NRHF(ISYMJ)
C
               KOFF1  = I3ORHF(ISYMKI,ISYMJ) 
     *                + NMAIJK(ISYMKI)*(J - 1) 
     *                + IMAIJK(ISYMMK,ISYMI) + 1
               KOFF2  = IT2BCD(ISYMAI,ISYMJ) + NT1AM(ISYMAI)*(J - 1)
     *                + IT1AM(ISYMA,ISYMI)   + 1
C
               NTOTA  = MAX(NVIR(ISYMA),1)
               NTOTMK = MAX(NMATIJ(ISYMMK),1)
C
               CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),
     *                    NMATIJ(ISYMMK),ONE,WORK(KTREOR),NTOTA,
     *                    XMIRES(KOFF1),NTOTMK,ONE,D2IAJ(KOFF2),NTOTA)
C
  150       CONTINUE
  140    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('DIAJ25')
C
      RETURN
      END
C  /* Deck dabgao */
      SUBROUTINE DABGAO(D2ABG,ISYDEN,T2INT,ISYMTI,Z2INT,ISYMZI,
     *                  WORK,LWORK,G,ISYMG)
C
C     Written by Asger Halkier 27/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the only contribution from D2ABC
C              (originating from projection against doubles) 
C              to D2ABG directly!
C
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
#if defined (SYS_CRAY)
      REAL ZERO, ONE
#else
      DOUBLE PRECISION ZERO, ONE
#endif
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)

      INTEGER ISYDEN, ISYMTI, ISYMZI, ISYMG, LWORK

#if defined (SYS_CRAY)
      REAL D2ABG(*), T2INT(*), Z2INT(*), WORK(LWORK)
#else
      DOUBLE PRECISION D2ABG(*), T2INT(*), Z2INT(*), WORK(LWORK)
#endif

      INTEGER ISYMA,ISYMAB,ISYMB,ISYMM, ISYMAM, ISYMBM, ISYMMN, ISYMN, 
     &        KZREOR, KTREOR, KEND1, LWRK1, kOFF1, NTOTA, NTOTMN, 
     &        NMN, NAMNN, NAMNR, NBMN, NMNB
C
      CALL QENTER('DABGAO')
C
      IF (MULD2H(ISYMTI,ISYMZI).NE.MULD2H(ISYDEN,ISYMG)) THEN
        WRITE(LUPRI,*) 'ISYMTI,ISYMZI,ISYMG,ISYDEN:',
     &                  ISYMTI,ISYMZI,ISYMG,ISYDEN
        CALL QUIT('Symmetry mismatch in DABGAO. (1)')
      END IF
C
      ISYMAB = MULD2H(ISYDEN,ISYMG)
C
      DO 100 ISYMA = 1,NSYM
C
         ISYMB  = MULD2H(ISYMA,ISYMAB)
         ISYMMN = MULD2H(ISYMA,ISYMZI)
         IF (MULD2H(ISYMB,ISYMTI).NE.ISYMMN) THEN
           CALL QUIT('Symmetry mismatch in DABGAO. (2)')
         END IF
C
C------------------------------
C        Work space allocation.
C------------------------------
C
         KZREOR = 1
         KTREOR = KZREOR + NVIR(ISYMA)*NMATIJ(ISYMMN)
         KEND1  = KTREOR + NMATIJ(ISYMMN)*NVIR(ISYMB)
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT('Insufficient work space in routine DABGAO')
         ENDIF
C
         CALL DZERO(WORK(KZREOR),NVIR(ISYMA)*NMATIJ(ISYMMN))
         CALL DZERO(WORK(KTREOR),NMATIJ(ISYMMN)*NVIR(ISYMB))
C
         DO 110 ISYMM = 1,NSYM
C
            ISYMN  = MULD2H(ISYMM,ISYMMN)
C
C-----------------------------------------------
C           Reorder z2int in scratch-work-array.
C-----------------------------------------------
C
            ISYMAM = MULD2H(ISYMN,ISYMZI)

            DO 120 N = 1,NRHF(ISYMN)
C
               DO 130 M = 1,NRHF(ISYMM)
C
                  NMN   = IMATIJ(ISYMM,ISYMN) + NRHF(ISYMM)*(N - 1) + M
                  NAMNN = IT2BCD(ISYMAM,ISYMN) + NT1AM(ISYMAM)*(N - 1)
     *                  + IT1AM(ISYMA,ISYMM) + NVIR(ISYMA)*(M - 1) + 1
                  NAMNR = KZREOR + NVIR(ISYMA)*(NMN - 1)
C
                  CALL DCOPY(NVIR(ISYMA),Z2INT(NAMNN),1,WORK(NAMNR),1)
C
  130          CONTINUE
  120       CONTINUE
C
C-----------------------------------------------
C           Reorder t2int in scratch-work-array.
C-----------------------------------------------
C
            ISYMBM = MULD2H(ISYMN,ISYMTI)

            DO 140 N = 1,NRHF(ISYMN)
C
               DO 150 M = 1,NRHF(ISYMM)
C
                  NMN  = IMATIJ(ISYMM,ISYMN)  + NRHF(ISYMM)*(N - 1) + M
                  NBMN = IT2BCD(ISYMBM,ISYMN) + NT1AM(ISYMBM)*(N - 1)
     *                 + IT1AM(ISYMB,ISYMM)   + NVIR(ISYMB)*(M - 1) + 1
                  NMNB = KTREOR + NMN - 1
C
                  CALL DCOPY(NVIR(ISYMB),T2INT(NBMN),1,WORK(NMNB),
     *                       NMATIJ(ISYMMN))
C
  150          CONTINUE
  140       CONTINUE
  110    CONTINUE
C
C------------------------------------------------------------------
C        Contraction of two scratch-matrices and storage in result.
C------------------------------------------------------------------
C
         KOFF1  = ID2ABG(ISYMAB,ISYMG) + NMATAB(ISYMAB)*(G - 1)
     *          + IMATAB(ISYMA,ISYMB)  + 1
C
         NTOTA  = MAX(NVIR(ISYMA),1)
         NTOTMN = MAX(NMATIJ(ISYMMN),1)
C
         CALL DGEMM('N','N',NVIR(ISYMA),NVIR(ISYMB),NMATIJ(ISYMMN),
     *              ONE,WORK(KZREOR),NTOTA,WORK(KTREOR),NTOTMN,ONE,
     *              D2ABG(KOFF1),NTOTA)
C
  100 CONTINUE
C
      CALL QEXIT('DABGAO')
C
      RETURN
      END
C  /* Deck cc1intmo */
      SUBROUTINE CCDINTMO(XINTIJ,XINTIA,XINTAB,XINTAI,XAOINT,
     *                    XLAMDP,XLAMDH,WORK,LWORK,ISYM)
C
C     Written by Asger Halkier 19/3 - 1996
C
C     Version: 1.0
C
C     Purpose: To transform the AO integrals in (XAOINT)
C              to MO-basis (stored in XINTIJ through XINTAI)!
C              Note that the AO integrals either can be the one
C              electron integrals or a submatrix alfa,beta of the two
C              electron integrals (alfa beta|gamma delta)!
C              ISYM is the integral symmetry!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION XINTIJ(*), XINTIA(*), XINTAB(*), XINTAI(*) 
      DIMENSION XAOINT(*), XLAMDP(*), XLAMDH(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('CCDINTMO')
C
C------------------------------
C     Initialize result arrays.
C------------------------------
C
      CALL DZERO(XINTIJ,NMATIJ(ISYM))
      CALL DZERO(XINTIA,NT1AM(ISYM))
      CALL DZERO(XINTAB,NMATAB(ISYM))
      CALL DZERO(XINTAI,NT1AM(ISYM))
C
      DO 100 ISYMAL = 1,NSYM
C
         ISYMBE = MULD2H(ISYMAL,ISYM)
         ISYMP  = ISYMAL
         ISYMQ  = ISYMBE
C
C------------------------------
C        Work space allocation.
C------------------------------
C
         LESCR = MAX(NBAS(ISYMAL)*NRHF(ISYMQ),NBAS(ISYMAL)*NVIR(ISYMQ))
C
         KSCR  = 1
         KEND1 = KSCR  + LESCR
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT('Insufficient memory for allocation in CCDINTMO')
         ENDIF
C
         CALL DZERO(WORK(KSCR),NBAS(ISYMAL)*NRHF(ISYMQ))
C
C----------------------------------------------------------
C        Transform second integral index to occupied space.
C----------------------------------------------------------
C
         ISYMI  = ISYMQ
C
         KOFF1  = IAODIS(ISYMAL,ISYMBE) + 1
         KOFF2  = IGLMRH(ISYMBE,ISYMI) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTBE = MAX(NBAS(ISYMBE),1)
C
         CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMI),NBAS(ISYMBE),ONE,
     *              XAOINT(KOFF1),NTOTAL,XLAMDH(KOFF2),NTOTBE,ZERO,
     *              WORK(KSCR),NTOTAL)
C
C------------------------------------------
C        Calculate the (occ.occ) integrals.
C------------------------------------------
C
         ISYMJ  = ISYMP
C
         KOFF3  = IGLMRH(ISYMAL,ISYMJ) + 1
         KOFF4  = IMATIJ(ISYMJ,ISYMI) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTJ  = MAX(NRHF(ISYMJ),1)
C
         CALL DGEMM('T','N',NRHF(ISYMJ),NRHF(ISYMI),NBAS(ISYMAL),ONE,
     *              XLAMDP(KOFF3),NTOTAL,WORK(KSCR),NTOTAL,ZERO,
     *              XINTIJ(KOFF4),NTOTJ)
C
C------------------------------------------
C        Calculate the (vir.occ) integrals.
C------------------------------------------
C
         ISYMA  = ISYMP
C
         KOFF5  = IGLMVI(ISYMAL,ISYMA) + 1
         KOFF6  = IT1AM(ISYMA,ISYMI) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTA  = MAX(NVIR(ISYMA),1)
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMAL),ONE,
     *              XLAMDP(KOFF5),NTOTAL,WORK(KSCR),NTOTAL,ZERO,
     *              XINTAI(KOFF6),NTOTA)
C
C---------------------------------------------------------
C        Transform second integral index to virtual space.
C---------------------------------------------------------
C
         CALL DZERO(WORK(KSCR),NBAS(ISYMAL)*NVIR(ISYMQ))
C
         ISYMA  = ISYMQ
C
         KOFF7  = IAODIS(ISYMAL,ISYMBE) + 1
         KOFF8  = IGLMVI(ISYMBE,ISYMA) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTBE = MAX(NBAS(ISYMBE),1)
C
         CALL DGEMM('N','N',NBAS(ISYMAL),NVIR(ISYMA),NBAS(ISYMBE),ONE,
     *              XAOINT(KOFF7),NTOTAL,XLAMDH(KOFF8),NTOTBE,ZERO,
     *              WORK(KSCR),NTOTAL)
C
C------------------------------------------
C        Calculate the (occ.vir) integrals.
C        Note that these are stored trans-
C        posed, i.e. (vir.occ) like a t1!
C------------------------------------------
C
         ISYMI  = ISYMP
C
         KOFF9  = IGLMRH(ISYMAL,ISYMI) + 1
         KOFF10 = IT1AM(ISYMA,ISYMI) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTA  = MAX(NVIR(ISYMA),1)
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMAL),ONE,
     *              WORK(KSCR),NTOTAL,XLAMDP(KOFF9),NTOTAL,ZERO,
     *              XINTIA(KOFF10),NTOTA)
C
C------------------------------------------
C        Calculate the (vir.vir) integrals.
C------------------------------------------
C
         ISYMB  = ISYMP
C
         KOFF11 = IGLMVI(ISYMAL,ISYMB) + 1
         KOFF12 = IMATAB(ISYMB,ISYMA) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTB  = MAX(NVIR(ISYMB),1)
C
         CALL DGEMM('T','N',NVIR(ISYMB),NVIR(ISYMA),NBAS(ISYMAL),ONE,
     *              XLAMDP(KOFF11),NTOTAL,WORK(KSCR),NTOTAL,ZERO,
     *              XINTAB(KOFF12),NTOTB)
C
  100 CONTINUE
C
      CALL QEXIT('CCDINTMO')
C
      RETURN
      END
C  /* Deck ccdenzk0 */
      SUBROUTINE CCDENZK0(ETAKA,XINTIJ,XINTAI,XINTIA,XINTAB,DIJ,DAI,
     *                    DIA,DAB,XIJT,XABT,DIJT,DABT,T1AM,
     *                    WORK,LWORK,ISYM)
C
C     Written by Asger Halkier 20/3 - 1996
C
C     Version: 1.0
C
C     Purpose: To set up the right hand side of the equation for
C              zeta-kappa-0 (ETAKA) from MO-integrals (XI*), Coupled 
C              Cluster densities (D*) and t1-amplitudes (T1AM)!
C              Note that due to the symmetry in the formulas, this
C              routine is able to handle both the one- and the two-
C              electron contributions!
C              ISYM is the symmetry of both the density and the 
C              integrals!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION ETAKA(*), XINTIJ(*), XINTAI(*), XINTIA(*), XINTAB(*)
      DIMENSION DIJ(*), DAI(*), DIA(*), DAB(*), XIJT(*), XABT(*)
      DIMENSION DIJT(*), DABT(*), T1AM(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('CCDENZK0')
C
      DO 100 ISYMA = 1,NSYM
C
C-------------------------------------------------------
C        Calculate terms originating from [H(t1),E(ai)].
C-------------------------------------------------------
C
         ISYMI  = ISYMA
         ISYMB  = MULD2H(ISYMA,ISYM)
         ISYMJ  = MULD2H(ISYMA,ISYM)
C
         KOFFRE = IT1AM(ISYMA,ISYMI)  + 1
C
         NTOTRE = MAX(NVIR(ISYMA),1)
         NTOTI  = MAX(NRHF(ISYMI),1)
         NTOTB  = MAX(NVIR(ISYMB),1)
         NTOTJ  = MAX(NRHF(ISYMJ),1)
C
         KOFF1  = IMATAB(ISYMA,ISYMB) + 1
         KOFF2  = IT1AM(ISYMB,ISYMI)  + 1
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMB),
     *              ONE,XABT(KOFF1),NTOTRE,DAI(KOFF2),NTOTB,ONE,
     *              ETAKA(KOFFRE),NTOTRE)
C
         KOFF3  = IMATAB(ISYMA,ISYMB) + 1
         KOFF4  = IT1AM(ISYMB,ISYMI)  + 1
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMB),
     *              -ONE,DAB(KOFF3),NTOTRE,XINTIA(KOFF4),NTOTB,ONE,
     *              ETAKA(KOFFRE),NTOTRE)
C
         KOFF5  = IT1AM(ISYMA,ISYMJ)  + 1
         KOFF6  = IMATIJ(ISYMJ,ISYMI) + 1
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     *              ONE,XINTIA(KOFF5),NTOTRE,DIJ(KOFF6),NTOTJ,ONE,
     *              ETAKA(KOFFRE),NTOTRE)
C
         KOFF7  = IT1AM(ISYMA,ISYMJ)  + 1
         KOFF8  = IMATIJ(ISYMJ,ISYMI) + 1
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     *              -ONE,DAI(KOFF7),NTOTRE,XIJT(KOFF8),NTOTJ,ONE,
     *              ETAKA(KOFFRE),NTOTRE)
C
C-------------------------------------------------------
C        Calculate terms originating from [H(t1),E(ia)].
C-------------------------------------------------------
C
         KOFF9  = IMATAB(ISYMA,ISYMB) + 1
         KOFF10 = IT1AM(ISYMB,ISYMI)  + 1
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMB),
     *              -ONE,DABT(KOFF9),NTOTRE,XINTAI(KOFF10),NTOTB,
     *              ONE,ETAKA(KOFFRE),NTOTRE)
C
         KOFF11 = IMATAB(ISYMA,ISYMB) + 1
         KOFF12 = IT1AM(ISYMB,ISYMI)  + 1
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMB),
     *              ONE,XINTAB(KOFF11),NTOTRE,DIA(KOFF12),NTOTB,
     *              ONE,ETAKA(KOFFRE),NTOTRE)
C
         KOFF13 = IT1AM(ISYMA,ISYMJ)  + 1
         KOFF14 = IMATIJ(ISYMJ,ISYMI) + 1
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     *              -ONE,DIA(KOFF13),NTOTRE,XINTIJ(KOFF14),NTOTJ,
     *              ONE,ETAKA(KOFFRE),NTOTRE)
C
         KOFF15 = IT1AM(ISYMA,ISYMJ)  + 1
         KOFF16 = IMATIJ(ISYMJ,ISYMI) + 1
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMJ),
     *              ONE,XINTAI(KOFF15),NTOTRE,DIJT(KOFF16),NTOTJ,
     *              ONE,ETAKA(KOFFRE),NTOTRE)
C
  100 CONTINUE
C
C----------------------------------------------------
C     Calculate terms originating from [H(t1),E(ik)].
C----------------------------------------------------
C
      DO 110 ISYMA = 1,NSYM
C
         ISYMI = ISYMA
         ISYMK = ISYMA
         ISYMB = MULD2H(ISYMI,ISYM)
         ISYMJ = MULD2H(ISYMI,ISYM)
C
         KOFFR = IT1AM(ISYMA,ISYMI) + 1
         KOFFT = IT1AM(ISYMA,ISYMK) + 1
C
         NTOTR = MAX(NVIR(ISYMA),1)
         NTOTB = MAX(NVIR(ISYMB),1)
         NTOTK = MAX(NRHF(ISYMK),1)
         NTOTI = MAX(NRHF(ISYMI),1)
         NTOTJ = MAX(NRHF(ISYMJ),1)
C
C----------------------------------
C        Work space allocation one.
C----------------------------------
C
         KSCR1 = 1
         KEND1 = KSCR1 + NRHF(ISYMK)*NRHF(ISYMI)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT('Insufficient memory for work '//
     &                'allocation in CCDENZK0')
         ENDIF
C
         CALL DZERO(WORK(KSCR1),NRHF(ISYMK)*NRHF(ISYMI))
C
         KOFF17 = IT1AM(ISYMB,ISYMK) + 1
         KOFF18 = IT1AM(ISYMB,ISYMI) + 1
C
         CALL DGEMM('T','N',NRHF(ISYMK),NRHF(ISYMI),NVIR(ISYMB),ONE,
     *              XINTIA(KOFF17),NTOTB,DIA(KOFF18),NTOTB,ONE,
     *              WORK(KSCR1),NTOTK)
C
         KOFF19 = IT1AM(ISYMB,ISYMK) + 1
         KOFF20 = IT1AM(ISYMB,ISYMI) + 1
C
         CALL DGEMM('T','N',NRHF(ISYMK),NRHF(ISYMI),NVIR(ISYMB),-ONE,
     *              DAI(KOFF19),NTOTB,XINTAI(KOFF20),NTOTB,ONE,
     *              WORK(KSCR1),NTOTK)
C
         KOFF21 = IMATIJ(ISYMK,ISYMJ) + 1
         KOFF22 = IMATIJ(ISYMJ,ISYMI) + 1
C
         CALL DGEMM('N','N',NRHF(ISYMK),NRHF(ISYMI),NRHF(ISYMJ),ONE,
     *              XINTIJ(KOFF21),NTOTK,DIJT(KOFF22),NTOTJ,ONE,
     *              WORK(KSCR1),NTOTK)
C
         KOFF23 = IMATIJ(ISYMK,ISYMJ) + 1
         KOFF24 = IMATIJ(ISYMJ,ISYMI) + 1
C
         CALL DGEMM('N','N',NRHF(ISYMK),NRHF(ISYMI),NRHF(ISYMJ),-ONE,
     *              DIJT(KOFF23),NTOTK,XINTIJ(KOFF24),NTOTJ,ONE,
     *              WORK(KSCR1),NTOTK)
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMK),ONE,
     *              T1AM(KOFFT),NTOTR,WORK(KSCR1),NTOTK,ONE,
     *              ETAKA(KOFFR),NTOTR)
C
  110 CONTINUE
C
C----------------------------------------------------
C     Calculate terms originating from [H(t1),E(ca)].
C----------------------------------------------------
C
      DO 120 ISYMA = 1,NSYM
C
         ISYMI = ISYMA
         ISYMC = ISYMI
         ISYMB = MULD2H(ISYMA,ISYM)
         ISYMJ = MULD2H(ISYMA,ISYM)
C
         KOFFR = IT1AM(ISYMA,ISYMI) + 1
         KOFFT = IT1AM(ISYMC,ISYMI) + 1
C
         NTOTR = MAX(NVIR(ISYMA),1)
         NTOTC = MAX(NVIR(ISYMC),1)
         NTOTB = MAX(NVIR(ISYMB),1)
C
C----------------------------------
C        Work space allocation two.
C----------------------------------
C
         KSCR2 = 1
         KEND2 = KSCR2 + NVIR(ISYMA)*NVIR(ISYMC)
         LWRK2 = LWORK - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
            CALL QUIT('Insufficient memory for work '//
     &                'allocation in CCDENZK0')
         ENDIF
C
         CALL DZERO(WORK(KSCR2),NVIR(ISYMA)*NVIR(ISYMC))
C
         KOFF25 = IT1AM(ISYMA,ISYMJ) + 1
         KOFF26 = IT1AM(ISYMC,ISYMJ) + 1
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMC),NRHF(ISYMJ),ONE,
     *              DIA(KOFF25),NTOTR,XINTIA(KOFF26),NTOTC,ONE,
     *              WORK(KSCR2),NTOTR)
C
         KOFF27 = IT1AM(ISYMA,ISYMJ) + 1
         KOFF28 = IT1AM(ISYMC,ISYMJ) + 1
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMC),NRHF(ISYMJ),-ONE,
     *              XINTAI(KOFF27),NTOTR,DAI(KOFF28),NTOTC,ONE,
     *              WORK(KSCR2),NTOTR)
C
         KOFF29 = IMATAB(ISYMA,ISYMB) + 1
         KOFF30 = IMATAB(ISYMB,ISYMC) + 1
C
         CALL DGEMM('N','N',NVIR(ISYMA),NVIR(ISYMC),NVIR(ISYMB),ONE,
     *              DABT(KOFF29),NTOTR,XINTAB(KOFF30),NTOTB,ONE,
     *              WORK(KSCR2),NTOTR)
C
         KOFF31 = IMATAB(ISYMA,ISYMB) + 1
         KOFF32 = IMATAB(ISYMB,ISYMC) + 1
C
         CALL DGEMM('N','N',NVIR(ISYMA),NVIR(ISYMC),NVIR(ISYMB),-ONE,
     *              XINTAB(KOFF31),NTOTR,DABT(KOFF32),NTOTB,ONE,
     *              WORK(KSCR2),NTOTR)
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMC),ONE,
     *              WORK(KSCR2),NTOTR,T1AM(KOFFT),NTOTC,ONE,
     *              ETAKA(KOFFR),NTOTR)
C
  120 CONTINUE
C
C----------------------------------------------------
C     Calculate terms originating from [H(t1),E(ck)].
C----------------------------------------------------
C
      DO 130 ISYMA = 1,NSYM
C
         ISYMI = ISYMA
         ISYMK = ISYMA
         ISYMC = ISYMI
         ISYMB = MULD2H(ISYMC,ISYM)
         ISYMJ = MULD2H(ISYMC,ISYM)
C
         KOFFR = IT1AM(ISYMA,ISYMI) + 1
         KOFTO = IT1AM(ISYMA,ISYMK) + 1
         KOFTI = IT1AM(ISYMC,ISYMI) + 1
C
         NTOTR = MAX(NVIR(ISYMA),1)
         NTOTC = MAX(NVIR(ISYMC),1)
         NTOTB = MAX(NVIR(ISYMB),1)
         NTOTK = MAX(NRHF(ISYMK),1)
         NTOTJ = MAX(NRHF(ISYMJ),1)
C
C------------------------------------
C        Work space allocation three.
C------------------------------------
C
         KSCR3 = 1
         KSCR4 = KSCR3 + NVIR(ISYMC)*NRHF(ISYMK)
         KEND3 = KSCR4 + NRHF(ISYMK)*NRHF(ISYMI)
         LWRK3 = LWORK - KEND3
C
         IF (LWRK3 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND3
            CALL QUIT('Insufficient memory for work allocation '//
     &                'in CCDENZK0')
         ENDIF
C
         CALL DZERO(WORK(KSCR3),NVIR(ISYMC)*NRHF(ISYMK))
         CALL DZERO(WORK(KSCR4),NRHF(ISYMK)*NRHF(ISYMI))
C
         KOFF33 = IMATAB(ISYMC,ISYMB) + 1
         KOFF34 = IT1AM(ISYMB,ISYMK)  + 1
C
         CALL DGEMM('N','N',NVIR(ISYMC),NRHF(ISYMK),NVIR(ISYMB),ONE,
     *              XABT(KOFF33),NTOTC,DAI(KOFF34),NTOTB,ONE,
     *              WORK(KSCR3),NTOTC)
C
         KOFF35 = IMATAB(ISYMC,ISYMB) + 1
         KOFF36 = IT1AM(ISYMB,ISYMK)  + 1
C
         CALL DGEMM('N','N',NVIR(ISYMC),NRHF(ISYMK),NVIR(ISYMB),-ONE,
     *              DAB(KOFF35),NTOTC,XINTIA(KOFF36),NTOTB,ONE,
     *              WORK(KSCR3),NTOTC)
C
         KOFF37 = IT1AM(ISYMC,ISYMJ)  + 1
         KOFF38 = IMATIJ(ISYMJ,ISYMK) + 1
C
         CALL DGEMM('N','N',NVIR(ISYMC),NRHF(ISYMK),NRHF(ISYMJ),ONE,
     *              XINTIA(KOFF37),NTOTC,DIJ(KOFF38),NTOTJ,ONE,
     *              WORK(KSCR3),NTOTC)
C
         KOFF39 = IT1AM(ISYMC,ISYMJ)  + 1
         KOFF40 = IMATIJ(ISYMJ,ISYMK) + 1
C
         CALL DGEMM('N','N',NVIR(ISYMC),NRHF(ISYMK),NRHF(ISYMJ),-ONE,
     *              DAI(KOFF39),NTOTC,XIJT(KOFF40),NTOTJ,ONE,
     *              WORK(KSCR3),NTOTC)
C
         CALL DGEMM('T','N',NRHF(ISYMK),NRHF(ISYMI),NVIR(ISYMC),ONE,
     *              WORK(KSCR3),NTOTC,T1AM(KOFTI),NTOTC,ZERO,
     *              WORK(KSCR4),NTOTK)
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMK),ONE,
     *              T1AM(KOFTO),NTOTR,WORK(KSCR4),NTOTK,ONE,
     *              ETAKA(KOFFR),NTOTR)
C
  130 CONTINUE
C
      CALL QEXIT('CCDENZK0')
C
      RETURN
      END
C  /* Deck cc_d2eff */
      SUBROUTINE CC_D2EFF(AODEN,G,ISYMG,IDEL,ISYMD,ZKABAO,DHFAO,ICON)
C
C     Written by Asger Halkier 12/5 - 1998
C
C     Version: 1.0
C
C     Purpose: To add the extra terms that add to the 2-elctron
C              coupled cluster lamda density to form the effective
C              2-electron density.
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION AODEN(*), ZKABAO(*), DHFAO(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('CC_D2EFF')
C
      FACI = ONE
      IF (ICON .EQ. 2) FACI = ONE/TWO
C
C-----------------------
C     Add coulomb terms.
C-----------------------
C
      ISALBE = MULD2H(ISYMG,ISYMD)
      D      = IDEL - IBAS(ISYMD)
C
      IF (ISYMG .EQ. ISYMD) THEN
         KOFFGD = IAODIS(ISYMG,ISYMD) + NBAS(ISYMG)*(D - 1) + G
         FAC1 = TWO*DHFAO(KOFFGD)*FACI
         CALL DAXPY(N2BST(ISALBE),FAC1,ZKABAO,1,AODEN,1)
      ENDIF
C
C------------------------
C     Add exchange terms.
C------------------------
C
      ISYMB = ISYMG
      ISYMA = ISYMD
C
      DO 100 B  = 1,NBAS(ISYMB)
C
         KOFFGB = IAODIS(ISYMG,ISYMB) + NBAS(ISYMG)*(B - 1) + G
         KOFFAD = IAODIS(ISYMA,ISYMD) + NBAS(ISYMA)*(D - 1) + 1
         KOFFAB = IAODIS(ISYMA,ISYMB) + NBAS(ISYMA)*(B - 1) + 1
C
         FAC2 = -DHFAO(KOFFGB)*FACI
         CALL DAXPY(NBAS(ISYMA),FAC2,ZKABAO(KOFFAD),1,
     *              AODEN(KOFFAB),1)
C
  100 CONTINUE
C
      CALL QEXIT('CC_D2EFF')
C
      RETURN
      END
C  /* Deck ccdiazk0 */
      SUBROUTINE CCDIAZK0(ZKDIA,XINTIJ,XINTAI,XINTIA,XINTAB,DIJ,DAI,
     *                    DIA,DAB,XIJT,XABT,DIJT,DABT,T1AM,
     *                    WORK,LWORK,ISYM)
C
C     Written by Asger Halkier 17/6 - 1998
C
C     Purpose: To set up the right hand side for the diagonal blocks
C              of kappa-bar-0 (ZKDIA) from MO-integrals (XI*), Coupled
C              Cluster densities (D*) and t1-amplitudes (T1AM).
C
C              Note that due to the symmetry in the formulas, this
C              routine is able to handle both the one- and the two-
C              electron contributions.
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION ZKDIA(*), XINTIJ(*), XINTAI(*), XINTIA(*), XINTAB(*)
      DIMENSION DIJ(*), DAI(*), DIA(*), DAB(*), XIJT(*), XABT(*)
      DIMENSION DIJT(*), DABT(*), T1AM(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('CCDIAZK0')
C
C===================================
C     First we do the virtual block.
C===================================
C
      DO 100 ISYMA = 1,NSYM
C
         ISYMB = ISYMA
         ISYMC = MULD2H(ISYMA,ISYM)
         ISYMI = MULD2H(ISYMA,ISYM)
C
         KOFFR = NMATIJ(1) + IMATAB(ISYMA,ISYMB) + 1
C
         NTOTA = MAX(NVIR(ISYMA),1)
         NTOTB = MAX(NVIR(ISYMB),1)
         NTOTC = MAX(NVIR(ISYMC),1)
C
C-----------------------------------------------------------------
C        Direct contributions with intermediate loop over virtual.
C-----------------------------------------------------------------
C
         KOFF1 = IMATAB(ISYMA,ISYMC) + 1
         KOFF2 = IMATAB(ISYMC,ISYMB) + 1
C
         CALL DGEMM('N','N',NVIR(ISYMA),NVIR(ISYMB),NVIR(ISYMC),ONE,
     *              XABT(KOFF1),NTOTA,DAB(KOFF2),NTOTC,ONE,
     *              ZKDIA(KOFFR),NTOTA)
C
         KOFF3 = IMATAB(ISYMA,ISYMC) + 1
         KOFF4 = IMATAB(ISYMC,ISYMB) + 1
C
         CALL DGEMM('N','N',NVIR(ISYMA),NVIR(ISYMB),NVIR(ISYMC),-ONE,
     *              DAB(KOFF3),NTOTA,XABT(KOFF4),NTOTC,ONE,
     *              ZKDIA(KOFFR),NTOTA)
C
         KOFF5 = IMATAB(ISYMA,ISYMC) + 1
         KOFF6 = IMATAB(ISYMC,ISYMB) + 1
C
         CALL DGEMM('N','N',NVIR(ISYMA),NVIR(ISYMB),NVIR(ISYMC),-ONE,
     *              DABT(KOFF5),NTOTA,XINTAB(KOFF6),NTOTC,ONE,
     *              ZKDIA(KOFFR),NTOTA)
C
         KOFF7 = IMATAB(ISYMA,ISYMC) + 1
         KOFF8 = IMATAB(ISYMC,ISYMB) + 1
C
         CALL DGEMM('N','N',NVIR(ISYMA),NVIR(ISYMB),NVIR(ISYMC),ONE,
     *              XINTAB(KOFF7),NTOTA,DABT(KOFF8),NTOTC,ONE,
     *              ZKDIA(KOFFR),NTOTA)
C
C------------------------------------------------------------------
C        Direct contributions with intermediate loop over occupied.
C------------------------------------------------------------------
C
         KOFF9  = IT1AM(ISYMA,ISYMI) + 1
         KOFF10 = IT1AM(ISYMB,ISYMI) + 1
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMI),ONE,
     *              XINTIA(KOFF9),NTOTA,DIA(KOFF10),NTOTB,ONE,
     *              ZKDIA(KOFFR),NTOTA)
C
         KOFF11 = IT1AM(ISYMA,ISYMI) + 1
         KOFF12 = IT1AM(ISYMB,ISYMI) + 1
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMI),-ONE,
     *              DAI(KOFF11),NTOTA,XINTAI(KOFF12),NTOTB,ONE,
     *              ZKDIA(KOFFR),NTOTA)
C
         KOFF13 = IT1AM(ISYMA,ISYMI) + 1
         KOFF14 = IT1AM(ISYMB,ISYMI) + 1
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMI),-ONE,
     *              DIA(KOFF13),NTOTA,XINTIA(KOFF14),NTOTB,ONE,
     *              ZKDIA(KOFFR),NTOTA)
C
         KOFF15 = IT1AM(ISYMA,ISYMI) + 1
         KOFF16 = IT1AM(ISYMB,ISYMI) + 1
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMI),ONE,
     *              XINTAI(KOFF15),NTOTA,DAI(KOFF16),NTOTB,ONE,
     *              ZKDIA(KOFFR),NTOTA)
C
  100 CONTINUE
C
      DO 110 ISYMA = 1,NSYM
C
         ISYMB = ISYMA
         ISYMK = ISYMA
         ISYMC = MULD2H(ISYMA,ISYM)
         ISYMI = MULD2H(ISYMA,ISYM)
C
         KOFFR = NMATIJ(1) + IMATAB(ISYMA,ISYMB) + 1
         KOFFT = IT1AM(ISYMA,ISYMK) + 1
C
         NTOTA = MAX(NVIR(ISYMA),1)
         NTOTB = MAX(NVIR(ISYMB),1)
         NTOTC = MAX(NVIR(ISYMC),1)
         NTOTI = MAX(NRHF(ISYMI),1)
C
C----------------------------------
C        Work space allocation one.
C----------------------------------
C
         KSCR1 = 1
         KEND1 = KSCR1 + NVIR(ISYMA)*NRHF(ISYMK)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT('Insufficient memory for work allocation '//
     &                'in CCDIAZK0')
         ENDIF
C
         LEN = NVIR(ISYMA)*NRHF(ISYMK)
C
         CALL DZERO(WORK(KSCR1),LEN)
C
C-------------------------------------------------------------------
C        Indirect contributions with intermediate loop over virtual.
C-------------------------------------------------------------------
C
         KOFF1 = IMATAB(ISYMA,ISYMC) + 1
         KOFF2 = IT1AM(ISYMC,ISYMK)  + 1
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMK),NVIR(ISYMC),ONE,
     *              XABT(KOFF1),NTOTA,DAI(KOFF2),NTOTC,ZERO,
     *              WORK(KSCR1),NTOTA)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMK),ONE,
     *              WORK(KSCR1),NTOTA,T1AM(KOFFT),NTOTB,ONE,
     *              ZKDIA(KOFFR),NTOTA)
C
         CALL DZERO(WORK(KSCR1),LEN)
C
         KOFF3 = IMATAB(ISYMB,ISYMC) + 1
         KOFF4 = IT1AM(ISYMC,ISYMK)  + 1
C
         CALL DGEMM('N','N',NVIR(ISYMB),NRHF(ISYMK),NVIR(ISYMC),ONE,
     *              XABT(KOFF3),NTOTB,DAI(KOFF4),NTOTC,ZERO,
     *              WORK(KSCR1),NTOTB)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMK),-ONE,
     *              T1AM(KOFFT),NTOTA,WORK(KSCR1),NTOTB,ONE,
     *              ZKDIA(KOFFR),NTOTA)
C
         CALL DZERO(WORK(KSCR1),LEN)
C
         KOFF5 = IMATAB(ISYMA,ISYMC) + 1
         KOFF6 = IT1AM(ISYMC,ISYMK)  + 1
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMK),NVIR(ISYMC),ONE,
     *              DAB(KOFF5),NTOTA,XINTIA(KOFF6),NTOTC,ZERO,
     *              WORK(KSCR1),NTOTA)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMK),-ONE,
     *              WORK(KSCR1),NTOTA,T1AM(KOFFT),NTOTB,ONE,
     *              ZKDIA(KOFFR),NTOTA)
C
         CALL DZERO(WORK(KSCR1),LEN)
C
         KOFF7 = IMATAB(ISYMB,ISYMC) + 1
         KOFF8 = IT1AM(ISYMC,ISYMK)  + 1
C
         CALL DGEMM('N','N',NVIR(ISYMB),NRHF(ISYMK),NVIR(ISYMC),ONE,
     *              DAB(KOFF7),NTOTB,XINTIA(KOFF8),NTOTC,ZERO,
     *              WORK(KSCR1),NTOTB)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMK),ONE,
     *              T1AM(KOFFT),NTOTA,WORK(KSCR1),NTOTB,ONE,
     *              ZKDIA(KOFFR),NTOTA)
C
         CALL DZERO(WORK(KSCR1),LEN)
C
C--------------------------------------------------------------------
C        Indirect contributions with intermediate loop over occupied.
C--------------------------------------------------------------------
C
         KOFF9  = IT1AM(ISYMA,ISYMI)  + 1
         KOFF10 = IMATIJ(ISYMI,ISYMK) + 1
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMK),NRHF(ISYMI),ONE,
     *              XINTIA(KOFF9),NTOTA,DIJ(KOFF10),NTOTI,ZERO,
     *              WORK(KSCR1),NTOTA)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMK),ONE,
     *              WORK(KSCR1),NTOTA,T1AM(KOFFT),NTOTB,ONE,
     *              ZKDIA(KOFFR),NTOTA)
C
         CALL DZERO(WORK(KSCR1),LEN)
C
         KOFF11 = IT1AM(ISYMB,ISYMI)  + 1
         KOFF12 = IMATIJ(ISYMI,ISYMK) + 1
C
         CALL DGEMM('N','N',NVIR(ISYMB),NRHF(ISYMK),NRHF(ISYMI),ONE,
     *              XINTIA(KOFF11),NTOTB,DIJ(KOFF12),NTOTI,ZERO,
     *              WORK(KSCR1),NTOTB)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMK),-ONE,
     *              T1AM(KOFFT),NTOTA,WORK(KSCR1),NTOTB,ONE,
     *              ZKDIA(KOFFR),NTOTA)
C
         CALL DZERO(WORK(KSCR1),LEN)
C
         KOFF13 = IT1AM(ISYMA,ISYMI)  + 1
         KOFF14 = IMATIJ(ISYMI,ISYMK) + 1
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMK),NRHF(ISYMI),ONE,
     *              DAI(KOFF13),NTOTA,XIJT(KOFF14),NTOTI,ZERO,
     *              WORK(KSCR1),NTOTA)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMK),-ONE,
     *              WORK(KSCR1),NTOTA,T1AM(KOFFT),NTOTB,ONE,
     *              ZKDIA(KOFFR),NTOTA)
C
         CALL DZERO(WORK(KSCR1),LEN)
C
         KOFF15 = IT1AM(ISYMB,ISYMI)  + 1
         KOFF16 = IMATIJ(ISYMI,ISYMK) + 1
C
         CALL DGEMM('N','N',NVIR(ISYMB),NRHF(ISYMK),NRHF(ISYMI),ONE,
     *              DAI(KOFF15),NTOTB,XIJT(KOFF16),NTOTI,ZERO,
     *              WORK(KSCR1),NTOTB)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMK),ONE,
     *              T1AM(KOFFT),NTOTA,WORK(KSCR1),NTOTB,ONE,
     *              ZKDIA(KOFFR),NTOTA)
C
  110 CONTINUE
C
C===================================
C     Then we do the occupied block.
C===================================
C
      DO 120 ISYMI = 1,NSYM
C
         ISYMJ = ISYMI
         ISYMC = MULD2H(ISYMI,ISYM)
         ISYMK = MULD2H(ISYMI,ISYM)
C
         KOFFR = IMATIJ(ISYMI,ISYMJ) + 1
C
         NTOTI = MAX(NRHF(ISYMI),1)
         NTOTC = MAX(NVIR(ISYMC),1)
         NTOTK = MAX(NRHF(ISYMK),1)
C
C-----------------------------------------------------------------
C        Direct contributions with intermediate loop over virtual.
C-----------------------------------------------------------------
C
         KOFF17 = IT1AM(ISYMC,ISYMI) + 1
         KOFF18 = IT1AM(ISYMC,ISYMJ) + 1
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMC),ONE,
     *              XINTAI(KOFF17),NTOTC,DAI(KOFF18),NTOTC,ONE,
     *              ZKDIA(KOFFR),NTOTI) 
C
         KOFF19 = IT1AM(ISYMC,ISYMI) + 1
         KOFF20 = IT1AM(ISYMC,ISYMJ) + 1
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMC),-ONE,
     *              DAI(KOFF19),NTOTC,XINTAI(KOFF20),NTOTC,ONE,
     *              ZKDIA(KOFFR),NTOTI) 
C
         KOFF21 = IT1AM(ISYMC,ISYMI) + 1
         KOFF22 = IT1AM(ISYMC,ISYMJ) + 1
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMC),-ONE,
     *              DIA(KOFF21),NTOTC,XINTIA(KOFF22),NTOTC,ONE,
     *              ZKDIA(KOFFR),NTOTI) 
C
         KOFF23 = IT1AM(ISYMC,ISYMI) + 1
         KOFF24 = IT1AM(ISYMC,ISYMJ) + 1
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMC),ONE,
     *              XINTIA(KOFF23),NTOTC,DIA(KOFF24),NTOTC,ONE,
     *              ZKDIA(KOFFR),NTOTI) 
C
C------------------------------------------------------------------
C        Direct contributions with intermediate loop over occupied.
C------------------------------------------------------------------
C
         KOFF25 = IMATIJ(ISYMI,ISYMK) + 1
         KOFF26 = IMATIJ(ISYMK,ISYMJ) + 1
C
         CALL DGEMM('N','N',NRHF(ISYMI),NRHF(ISYMJ),NRHF(ISYMK),ONE,
     *              XIJT(KOFF25),NTOTI,DIJ(KOFF26),NTOTK,ONE,
     *              ZKDIA(KOFFR),NTOTI)
C
         KOFF27 = IMATIJ(ISYMI,ISYMK) + 1
         KOFF28 = IMATIJ(ISYMK,ISYMJ) + 1
C
         CALL DGEMM('N','N',NRHF(ISYMI),NRHF(ISYMJ),NRHF(ISYMK),-ONE,
     *              DIJT(KOFF27),NTOTI,XINTIJ(KOFF28),NTOTK,ONE,
     *              ZKDIA(KOFFR),NTOTI)
C
         KOFF29 = IMATIJ(ISYMI,ISYMK) + 1
         KOFF30 = IMATIJ(ISYMK,ISYMJ) + 1
C
         CALL DGEMM('N','N',NRHF(ISYMI),NRHF(ISYMJ),NRHF(ISYMK),-ONE,
     *              DIJ(KOFF29),NTOTI,XIJT(KOFF30),NTOTK,ONE,
     *              ZKDIA(KOFFR),NTOTI)
C
         KOFF31 = IMATIJ(ISYMI,ISYMK) + 1
         KOFF32 = IMATIJ(ISYMK,ISYMJ) + 1
C
         CALL DGEMM('N','N',NRHF(ISYMI),NRHF(ISYMJ),NRHF(ISYMK),ONE,
     *              XINTIJ(KOFF31),NTOTI,DIJT(KOFF32),NTOTK,ONE,
     *              ZKDIA(KOFFR),NTOTI)
C
  120 CONTINUE
C
      DO 130 ISYMI = 1,NSYM
C
         ISYMJ = ISYMI
         ISYMC = ISYMI
         ISYMD = MULD2H(ISYMI,ISYM)
         ISYMK = MULD2H(ISYMI,ISYM)
C
         KOFFR = IMATIJ(ISYMI,ISYMJ) + 1
         KOFFT = IT1AM(ISYMC,ISYMI)  + 1
C
         NTOTI = MAX(NRHF(ISYMI),1)
         NTOTK = MAX(NRHF(ISYMK),1)
         NTOTC = MAX(NVIR(ISYMC),1)
         NTOTD = MAX(NVIR(ISYMD),1)
C
C----------------------------------
C        Work space allocation two.
C----------------------------------
C
         KSCR2 = 1
         KEND2 = KSCR2 + NVIR(ISYMC)*NRHF(ISYMI)
         LWRK2 = LWORK - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
            CALL QUIT('Insufficient memory for work allocation '//
     &                'in CCDIAZK0')
         ENDIF
C
         LEN = NVIR(ISYMC)*NRHF(ISYMI)
C
         CALL DZERO(WORK(KSCR2),LEN)
C
C-------------------------------------------------------------------
C        Indirect contributions with intermediate loop over virtual.
C-------------------------------------------------------------------
C
         KOFF33 = IMATAB(ISYMC,ISYMD) + 1
         KOFF34 = IT1AM(ISYMD,ISYMI)  + 1
C
         CALL DGEMM('N','N',NVIR(ISYMC),NRHF(ISYMI),NVIR(ISYMD),ONE,
     *              XABT(KOFF33),NTOTC,DAI(KOFF34),NTOTD,ZERO,
     *              WORK(KSCR2),NTOTC)
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMC),ONE,
     *              WORK(KSCR2),NTOTC,T1AM(KOFFT),NTOTC,ONE,
     *              ZKDIA(KOFFR),NTOTI)
C
         CALL DZERO(WORK(KSCR2),LEN)
C
         KOFF35 = IMATAB(ISYMC,ISYMD) + 1
         KOFF36 = IT1AM(ISYMD,ISYMJ)  + 1
C
         CALL DGEMM('N','N',NVIR(ISYMC),NRHF(ISYMJ),NVIR(ISYMD),ONE,
     *              XABT(KOFF35),NTOTC,DAI(KOFF36),NTOTD,ZERO,
     *              WORK(KSCR2),NTOTC)
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMC),-ONE,
     *              T1AM(KOFFT),NTOTC,WORK(KSCR2),NTOTC,ONE,
     *              ZKDIA(KOFFR),NTOTI)
C
         CALL DZERO(WORK(KSCR2),LEN)
C
         KOFF37 = IMATAB(ISYMC,ISYMD) + 1
         KOFF38 = IT1AM(ISYMD,ISYMI)  + 1
C
         CALL DGEMM('N','N',NVIR(ISYMC),NRHF(ISYMI),NVIR(ISYMD),ONE,
     *              DAB(KOFF37),NTOTC,XINTIA(KOFF38),NTOTD,ZERO,
     *              WORK(KSCR2),NTOTC)
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMC),-ONE,
     *              WORK(KSCR2),NTOTC,T1AM(KOFFT),NTOTC,ONE,
     *              ZKDIA(KOFFR),NTOTI)  
C
         CALL DZERO(WORK(KSCR2),LEN)
C
         KOFF39 = IMATAB(ISYMC,ISYMD) + 1
         KOFF40 = IT1AM(ISYMD,ISYMJ)  + 1
C
         CALL DGEMM('N','N',NVIR(ISYMC),NRHF(ISYMJ),NVIR(ISYMD),ONE,
     *              DAB(KOFF39),NTOTC,XINTIA(KOFF40),NTOTD,ZERO,
     *              WORK(KSCR2),NTOTC)
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMC),ONE,
     *              T1AM(KOFFT),NTOTC,WORK(KSCR2),NTOTC,ONE,
     *              ZKDIA(KOFFR),NTOTI)  
C
C-------------------------------------------------------------------
C        Indirect contributions with intermediate loop over virtual.
C-------------------------------------------------------------------
C
         CALL DZERO(WORK(KSCR2),LEN)
C
         KOFF41 = IT1AM(ISYMC,ISYMK)  + 1
         KOFF42 = IMATIJ(ISYMK,ISYMI) + 1
C
         CALL DGEMM('N','N',NVIR(ISYMC),NRHF(ISYMI),NRHF(ISYMK),ONE,
     *              XINTIA(KOFF41),NTOTC,DIJ(KOFF42),NTOTK,ZERO,
     *              WORK(KSCR2),NTOTC)
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMC),ONE,
     *              WORK(KSCR2),NTOTC,T1AM(KOFFT),NTOTC,ONE,
     *              ZKDIA(KOFFR),NTOTI)  
C
         CALL DZERO(WORK(KSCR2),LEN)
C
         KOFF43 = IT1AM(ISYMC,ISYMK)  + 1
         KOFF44 = IMATIJ(ISYMK,ISYMJ) + 1
C
         CALL DGEMM('N','N',NVIR(ISYMC),NRHF(ISYMJ),NRHF(ISYMK),ONE,
     *              XINTIA(KOFF43),NTOTC,DIJ(KOFF44),NTOTK,ZERO,
     *              WORK(KSCR2),NTOTC)
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMC),-ONE,
     *              T1AM(KOFFT),NTOTC,WORK(KSCR2),NTOTC,ONE,
     *              ZKDIA(KOFFR),NTOTI)
C
         CALL DZERO(WORK(KSCR2),LEN)
C
         KOFF45 = IT1AM(ISYMC,ISYMK)  + 1
         KOFF46 = IMATIJ(ISYMK,ISYMI) + 1
C
         CALL DGEMM('N','N',NVIR(ISYMC),NRHF(ISYMI),NRHF(ISYMK),ONE,
     *              DAI(KOFF45),NTOTC,XIJT(KOFF46),NTOTK,ZERO,
     *              WORK(KSCR2),NTOTC)
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMC),-ONE,
     *              WORK(KSCR2),NTOTC,T1AM(KOFFT),NTOTC,ONE,
     *              ZKDIA(KOFFR),NTOTI)  
C
         CALL DZERO(WORK(KSCR2),LEN)
C
         KOFF47 = IT1AM(ISYMC,ISYMK)  + 1
         KOFF48 = IMATIJ(ISYMK,ISYMJ) + 1
C
         CALL DGEMM('N','N',NVIR(ISYMC),NRHF(ISYMJ),NRHF(ISYMK),ONE,
     *              DAI(KOFF47),NTOTC,XIJT(KOFF48),NTOTK,ZERO,
     *              WORK(KSCR2),NTOTC)
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMC),ONE,
     *              T1AM(KOFFT),NTOTC,WORK(KSCR2),NTOTC,ONE,
     *              ZKDIA(KOFFR),NTOTI)  
C
  130 CONTINUE
C
      CALL QEXIT('CCDIAZK0')
C
      RETURN
      END
C  /* Deck ccsd_zkblo */
      SUBROUTINE CCSD_ZKBLO(ZKDIA,WORK,LWORK)
C
C     Written by Asger Halkier 4/8 - 1998
C
C     Version: 1.0
C
C     Purpose: To calculate the ab & ij parts of the CCSD kappa-bar-0,
C              from right-hand-sides (ZKDIA on input) and canonical
C              orbital energies.
C     Control numerical instabilities. Sonia, 2002
C
#include <implicit.h>
#include <dummy.h>
      PARAMETER(ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (THRDNM = 1.0D-12)
      DIMENSION ZKDIA(*), WORK(LWORK)
#include <priunit.h>
#include <maxorb.h>
#include <ccorb.h>
#include <iratdef.h>
#include <inftap.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccsdinp.h>
C
      CALL QENTER('CCSD_ZKBLO')
C
C---------------------------
C     Work space allocation.
C---------------------------
C
      KFOCKD = 1
      KETAIJ = KFOCKD + NORBTS
      KETAAB = KETAIJ + NMATIJ(1)
      KEND1  = KETAAB + NMATAB(1)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Need:', KEND1, 'Available:', LWORK
         CALL QUIT('Insufficient memory for allocation in CCSD_ZKBLO')
      ENDIF
C
      CALL DZERO(WORK,KEND1)
      CALL DCOPY(NMATIJ(1),ZKDIA(1),1,WORK(KETAIJ),1)
      CALL DCOPY(NMATAB(1),ZKDIA(NMATIJ(1)+1),1,WORK(KETAAB),1)
      CALL DZERO(ZKDIA,NMATIJ(1)+NMATAB(1))
C
C-------------------------------------
C     Read canonical orbital energies.
C-------------------------------------
C
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND (LUSIFC)
C
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KFOCKD + I - 1), I = 1,NORBTS)
C
      CALL GPCLOSE(LUSIFC,'KEEP')
C
C----------------------------------------------------------------
C     Change symmetry ordering of the canonical orbital energies.
C----------------------------------------------------------------
C
      IF (FROIMP)
     *    CALL CCSD_DELFRO(WORK(KFOCKD),WORK(KEND1),LWRK1)
C
      CALL FOCK_REORDER(WORK(KFOCKD),WORK(KEND1),LWRK1)
C
C---------------------------
C     Calculate the results:
C     Occupied block:
C---------------------------
C
      DO 100 ISYMI = 1,NSYM
         ISYMJ = ISYMI
         DO 110 J = 1,NRHF(ISYMJ)
            KOFFJ = KFOCKD + IRHF(ISYMJ) + J - 1
            DO 120 I = J+1,NRHF(ISYMI)
               KOFFI = KFOCKD + IRHF(ISYMI) + I - 1
               DENOM = WORK(KOFFJ) - WORK(KOFFI)
               IF (ABS(DENOM) .GT. THRDNM) THEN
                  NIJ = IMATIJ(ISYMI,ISYMJ) + NRHF(ISYMI)*(J - 1) + I
                  NJI = IMATIJ(ISYMJ,ISYMI) + NRHF(ISYMJ)*(I - 1) + J
C
                  ZKDIA(NIJ) = HALF*WORK(KETAIJ+NIJ-1)/DENOM
                  ZKDIA(NJI) = ZKDIA(NIJ)
               END IF
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
C-------------------
C     Virtual block:
C-------------------
C
      DO 130 ISYMA = 1,NSYM
         ISYMB = ISYMA
         DO 140 B = 1,NVIR(ISYMB)
            KOFFB = KFOCKD + IVIR(ISYMB) + B - 1
            DO 150 A = B+1,NVIR(ISYMA)
               KOFFA = KFOCKD + IVIR(ISYMA) + A - 1
               DENOM = WORK(KOFFB) - WORK(KOFFA)
               IF (ABS(DENOM) .GT. THRDNM) THEN
                  NAB = IMATAB(ISYMA,ISYMB) + NVIR(ISYMA)*(B - 1) + A
                  NBA = IMATAB(ISYMB,ISYMA) + NVIR(ISYMB)*(A - 1) + B
C
                  ZKDIA(NMATIJ(1) + NAB) = HALF*WORK(KETAAB+NAB-1)/
     *                                     DENOM
                  ZKDIA(NMATIJ(1) + NBA) = ZKDIA(NMATIJ(1) + NAB)
               END IF
C
  150       CONTINUE
  140    CONTINUE
  130 CONTINUE
C
      CALL QEXIT('CCSD_ZKBLO')
C
      RETURN
      END
C  /* Deck cc_d2gaf */
      SUBROUTINE CC_D2GAF(D2GIJ,D2GAB,D2GAI,D2GIA,DIJ,DAB,DAI,DIA,
     *                    CMO,IDEL,ISYMD,G,ISYMG)
C
C     Written by Asger Halkier 12/8 - 1998
C
C     Purpose: To calculate the contributions to d(pq,gam;del) where
C              gamma has been backtransformed from a frozen index.
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION D2GIJ(*), D2GAB(*), D2GAI(*), D2GIA(*)
      DIMENSION DIJ(*), DAB(*), DAI(*), DIA(*), CMO(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('CC_D2GAF')
C
      IF (ISYMG .EQ. ISYMD) THEN
C
         ISYML = ISYMD
C
         ND = ILRHSI(ISYMD) + IDEL - IBAS(ISYMD)
         NG = ILRHSI(ISYMG) + G
C
         FACT = TWO*DDOT(NRHFFR(ISYML),CMO(ND),NBAS(ISYMD),CMO(NG),
     *               NBAS(ISYMG))
C
         CALL DAXPY(NMATIJ(1),FACT,DIJ,1,D2GIJ,1)
         CALL DAXPY(NMATAB(1),FACT,DAB,1,D2GAB,1)
         CALL DAXPY(NT1AMX,FACT,DAI,1,D2GAI,1)
         CALL DAXPY(NT1AMX,FACT,DIA,1,D2GIA,1)
C 
      ENDIF
C
      CALL QEXIT('CC_D2GAF')
C
      RETURN
      END
C  /* Deck cc_fd2ao */
      SUBROUTINE CC_FD2AO(D2AO,D2II,D2IJ,D2JI,D2AI,D2IA,CMO,XLAMDH,
     *                    XLAMDP,WORK,LWORK,ISYMPQ)
C
C     Written by Asger Halkier 12/8 - 1998
C
C     Purpose: To calculate the contributions to D2AO from d(pq,gam;del) 
C              where at least one of the two indices p & q is frozen.
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION D2AO(*), D2II(*), D2IJ(*), D2JI(*), D2AI(*)
      DIMENSION D2IA(*), CMO(*), XLAMDH(*), XLAMDP(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
      CALL QENTER('CC_FD2AO')
C
C--------------------------------
C     Do the frozen-frozen block.
C--------------------------------
C
      DO 100 ISYMAL = 1,NSYM
C
         ISYMBE = MULD2H(ISYMPQ,ISYMAL)
         ISYMI  = ISYMAL
         ISYMJ  = ISYMBE
C
         IF ((NRHFFR(ISYMI) .EQ. 0).OR.(NRHFFR(ISYMJ) .EQ. 0)) GOTO 100
C
         KSCR1  = 1
         KEND1  = KSCR1 + NBAS(ISYMAL)*NRHFFR(ISYMJ)
         LWRK1  = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT('Insufficient work space in CC_FD2AO')
         ENDIF
C
         CALL DZERO(WORK(KSCR1),KEND1)
C
         KOFF1  = ILRHSI(ISYMAL) + 1
         KOFF2  = IFROFR(ISYMI,ISYMJ) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTI  = MAX(NRHFFR(ISYMI),1)
C
         CALL DGEMM('N','N',NBAS(ISYMAL),NRHFFR(ISYMJ),NRHFFR(ISYMI),
     *              ONE,CMO(KOFF1),NTOTAL,D2II(KOFF2),NTOTI,ZERO,
     *              WORK(KSCR1),NTOTAL)
C
         KOFF3  = ILRHSI(ISYMBE) + 1
         KOFF4  = IAODIS(ISYMAL,ISYMBE) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTBE = MAX(NBAS(ISYMBE),1)
C
         CALL DGEMM('N','T',NBAS(ISYMAL),NBAS(ISYMBE),NRHFFR(ISYMJ),
     *              ONE,WORK(KSCR1),NTOTAL,CMO(KOFF3),NTOTBE,ONE,
     *              D2AO(KOFF4),NTOTAL)
C
  100 CONTINUE
C
C------------------------------------
C     Do the correlated-frozen block.
C------------------------------------
C
      DO 110 ISYMAL = 1,NSYM
C
         ISYMBE = MULD2H(ISYMPQ,ISYMAL)
         ISYMI  = ISYMAL
         ISYMJ  = ISYMBE
C
         IF (NRHFFR(ISYMJ) .EQ. 0) GOTO 110
C
         KSCR1  = 1
         KEND1  = KSCR1 + NBAS(ISYMAL)*NRHFFR(ISYMJ)
         LWRK1  = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT('Insufficient work space in CC_FD2AO')
         ENDIF
C
         CALL DZERO(WORK(KSCR1),KEND1)
C
         KOFF5  = IGLMRH(ISYMAL,ISYMI) + 1
         KOFF6  = ICOFRO(ISYMI,ISYMJ) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTI  = MAX(NRHF(ISYMI),1)
C
         CALL DGEMM('N','N',NBAS(ISYMAL),NRHFFR(ISYMJ),NRHF(ISYMI),
     *              ONE,XLAMDP(KOFF5),NTOTAL,D2IJ(KOFF6),NTOTI,ZERO,
     *              WORK(KSCR1),NTOTAL)
C
         KOFF7  = ILRHSI(ISYMBE) + 1
         KOFF8  = IAODIS(ISYMAL,ISYMBE) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTBE = MAX(NBAS(ISYMBE),1)
C
         CALL DGEMM('N','T',NBAS(ISYMAL),NBAS(ISYMBE),NRHFFR(ISYMJ),
     *              ONE,WORK(KSCR1),NTOTAL,CMO(KOFF7),NTOTBE,ONE,
     *              D2AO(KOFF8),NTOTAL)
C
  110 CONTINUE
C
C------------------------------------
C     Do the frozen-correlated block.
C------------------------------------
C
      DO 120 ISYMAL = 1,NSYM
C
         ISYMBE = MULD2H(ISYMPQ,ISYMAL)
         ISYMI  = ISYMAL
         ISYMJ  = ISYMBE
C
         IF (NRHFFR(ISYMI) .EQ. 0) GOTO 120
C
         KSCR1  = 1
         KEND1  = KSCR1 + NBAS(ISYMBE)*NRHFFR(ISYMI)
         LWRK1  = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT('Insufficient work space in CC_FD2AO')
         ENDIF
C
         CALL DZERO(WORK(KSCR1),KEND1)
C
         KOFF9  = IGLMRH(ISYMBE,ISYMJ) + 1
         KOFF10 = ICOFRO(ISYMJ,ISYMI) + 1
C
         NTOTBE = MAX(NBAS(ISYMBE),1)
         NTOTJ  = MAX(NRHF(ISYMJ),1)
C
         CALL DGEMM('N','N',NBAS(ISYMBE),NRHFFR(ISYMI),NRHF(ISYMJ),
     *              ONE,XLAMDH(KOFF9),NTOTBE,D2JI(KOFF10),NTOTJ,
     *              ZERO,WORK(KSCR1),NTOTBE)
C
         KOFF11 = ILRHSI(ISYMAL) + 1
         KOFF12 = IAODIS(ISYMAL,ISYMBE) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTBE = MAX(NBAS(ISYMBE),1)
C
         CALL DGEMM('N','T',NBAS(ISYMAL),NBAS(ISYMBE),NRHFFR(ISYMI),
     *              ONE,CMO(KOFF11),NTOTAL,WORK(KSCR1),NTOTBE,ONE,
     *              D2AO(KOFF12),NTOTAL)
C
  120 CONTINUE
C
C---------------------------------
C     Do the virtual-frozen block.
C---------------------------------
C
      DO 130 ISYMAL = 1,NSYM
C
         ISYMBE = MULD2H(ISYMPQ,ISYMAL)
         ISYMA  = ISYMAL
         ISYMI  = ISYMBE
C
         IF (NRHFFR(ISYMI) .EQ. 0) GOTO 130
C
         KSCR1  = 1
         KEND1  = KSCR1 + NBAS(ISYMAL)*NRHFFR(ISYMI)
         LWRK1  = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT('Insufficient work space in CC_FD2AO')
         ENDIF
C
         CALL DZERO(WORK(KSCR1),KEND1)
C
         KOFF13 = IGLMVI(ISYMAL,ISYMA) + 1
         KOFF14 = IT1FRO(ISYMA,ISYMI) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTA  = MAX(NVIR(ISYMA),1)
C
         CALL DGEMM('N','N',NBAS(ISYMAL),NRHFFR(ISYMI),NVIR(ISYMA),
     *              ONE,XLAMDP(KOFF13),NTOTAL,D2AI(KOFF14),NTOTA,
     *              ZERO,WORK(KSCR1),NTOTAL)
C
         KOFF15 = ILRHSI(ISYMBE) + 1
         KOFF16 = IAODIS(ISYMAL,ISYMBE) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTBE = MAX(NBAS(ISYMBE),1)
C
         CALL DGEMM('N','T',NBAS(ISYMAL),NBAS(ISYMBE),NRHFFR(ISYMI),
     *              ONE,WORK(KSCR1),NTOTAL,CMO(KOFF15),NTOTBE,ONE,
     *              D2AO(KOFF16),NTOTAL)
C
  130 CONTINUE
C
C---------------------------------
C     Do the frozen-virtual block.
C---------------------------------
C
      DO 140 ISYMAL = 1,NSYM
C
         ISYMBE = MULD2H(ISYMPQ,ISYMAL)
         ISYMI  = ISYMAL
         ISYMA  = ISYMBE
C
         IF (NRHFFR(ISYMI) .EQ. 0) GOTO 140
C
         KSCR1  = 1
         KEND1  = KSCR1 + NBAS(ISYMBE)*NRHFFR(ISYMI)
         LWRK1  = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT('Insufficient work space in CC_FD2AO')
         ENDIF
C
         CALL DZERO(WORK(KSCR1),KEND1)
C
         KOFF17 = IGLMVI(ISYMBE,ISYMA) + 1
         KOFF18 = IT1FRO(ISYMA,ISYMI) + 1
C
         NTOTBE = MAX(NBAS(ISYMBE),1)
         NTOTA  = MAX(NVIR(ISYMA),1)
C
         CALL DGEMM('N','N',NBAS(ISYMBE),NRHFFR(ISYMI),NVIR(ISYMA),
     *              ONE,XLAMDH(KOFF17),NTOTBE,D2IA(KOFF18),NTOTA,
     *              ZERO,WORK(KSCR1),NTOTBE)
C
         KOFF19 = ILRHSI(ISYMAL) + 1
         KOFF20 = IAODIS(ISYMAL,ISYMBE) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTBE = MAX(NBAS(ISYMBE),1)
C
         CALL DGEMM('N','T',NBAS(ISYMAL),NBAS(ISYMBE),NRHFFR(ISYMI),
     *              ONE,CMO(KOFF19),NTOTAL,WORK(KSCR1),NTOTBE,ONE,
     *              D2AO(KOFF20),NTOTAL)
C
  140 CONTINUE
C
      CALL QEXIT('CC_FD2AO')
C
      RETURN
      END
C  /* Deck cc_fd2bl */
      SUBROUTINE CC_FD2BL(D2II,D2IJ,D2JI,D2AI,D2IA,DIJ,DAB,DAI,DIA,
     *                    CMO,XLAMDH,XLAMDP,WORK,LWORK,IDEL,ISYMD,
     *                    G,ISYMG)
C
C     Written by Asger Halkier 12/8 - 1998
C
C     Purpose: To calculate the contributions to d(pq,gam;del) 
C              where at least one of the two indices p & q is frozen.
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION D2II(*), D2IJ(*), D2JI(*), D2AI(*), D2IA(*), DIJ(*)
      DIMENSION DAB(*), DAI(*), DIA(*), CMO(*), XLAMDH(*)
      DIMENSION XLAMDP(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
      CALL QENTER('CC_FD2BL')
C
C-------------------------------
C     Work space allocation one.
C-------------------------------
C
      ISYMB  = ISYMD
      ISYMK  = ISYMD
      ISYMA  = ISYMB
      ISYMI  = ISYMG
C
      KVECA = 1
      KVECB = KVECA + NVIR(ISYMA)
      KVECK = KVECB + NVIR(ISYMB)
      KEND1 = KVECK + NRHF(ISYMK)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space in CC_FD2BL')
      ENDIF
C
      CALL DZERO(WORK(KVECA),KEND1)
C
      KOFF1 = IGLMVI(ISYMD,ISYMB) + IDEL - IBAS(ISYMD)
      KOFF2 = IGLMRH(ISYMD,ISYMK) + IDEL - IBAS(ISYMD)
      CALL DCOPY(NVIR(ISYMB),XLAMDH(KOFF1),NBAS(ISYMD),WORK(KVECB),1)
      CALL DCOPY(NRHF(ISYMK),XLAMDH(KOFF2),NBAS(ISYMD),WORK(KVECK),1)
C
C--------------------------------------
C     Calculate intermediate vector Va.
C--------------------------------------
C
      KOFF3 = IMATAB(ISYMA,ISYMB) + 1
      KOFF4 = IT1AM(ISYMA,ISYMK)  + 1
C
      NTOTA = MAX(NVIR(ISYMA),1)
C
      CALL DGEMV('N',NVIR(ISYMA),NVIR(ISYMB),ONE,DAB(KOFF3),NTOTA,
     *           WORK(KVECB),1,ZERO,WORK(KVECA),1)
C
      CALL DGEMV('N',NVIR(ISYMA),NRHF(ISYMK),ONE,DAI(KOFF4),NTOTA,
     *           WORK(KVECK),1,ONE,WORK(KVECA),1)
C
C------------------------------------
C     Calculate virtual-frozen block.
C------------------------------------
C
      DO 100 I = 1,NRHFFR(ISYMI)
C
         KOFF5 = ILRHSI(ISYMG) + NBAS(ISYMG)*(I - 1) + G
         KOFF6 = IT1FRO(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + 1
         CALL DAXPY(NVIR(ISYMA),-CMO(KOFF5),WORK(KVECA),1,D2AI(KOFF6),1)
C
  100 CONTINUE
C
C-----------------------------------------------------
C     Add contribution to frozen-frozen block from Va.
C-----------------------------------------------------
C
      IF (ISYMG .EQ. ISYMD) THEN
C
         KOFF7 = IGLMVI(ISYMG,ISYMA) + G
         FAC = DDOT(NVIR(ISYMA),WORK(KVECA),1,XLAMDP(KOFF7),NBAS(ISYMG))
C
         DO 110 ISYMJ = 1,NSYM
            DO 120 J = 1,NRHFFR(ISYMJ)
               KOFF8 = IFROFR(ISYMJ,ISYMJ) + NRHFFR(ISYMJ)*(J - 1) + J
               D2II(KOFF8) = D2II(KOFF8) + TWO*FAC
  120       CONTINUE
  110    CONTINUE
      ENDIF
C
C-------------------------------
C     Work space allocation two.
C-------------------------------
C
      ISYMA = ISYMD
      ISYML = ISYMD
      ISYMI = ISYMA
      ISYMJ = ISYMG
C
      KVECI = 1
      KVECA = KVECI + NRHF(ISYMI)
      KVECL = KVECA + NVIR(ISYMA)
      KEND1 = KVECL + NRHF(ISYML)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space in CC_FD2BL')
      ENDIF
C
      CALL DZERO(WORK(KVECI),KEND1)
C
      KOFF9  = IGLMVI(ISYMD,ISYMA) + IDEL - IBAS(ISYMD)
      KOFF10 = IGLMRH(ISYMD,ISYML) + IDEL - IBAS(ISYMD)
      CALL DCOPY(NVIR(ISYMA),XLAMDH(KOFF9),NBAS(ISYMD),WORK(KVECA),1)
      CALL DCOPY(NRHF(ISYML),XLAMDH(KOFF10),NBAS(ISYMD),WORK(KVECL),1)
C
C--------------------------------------
C     Calculate intermediate vector Vi.
C--------------------------------------
C
      KOFF11 = IT1AM(ISYMA,ISYMI)  + 1
      KOFF12 = IMATIJ(ISYMI,ISYML) + 1
C
      NTOTA  = MAX(NVIR(ISYMA),1)
      NTOTI  = MAX(NRHF(ISYMI),1)
C
      CALL DGEMV('T',NVIR(ISYMA),NRHF(ISYMI),ONE,DIA(KOFF11),NTOTA,
     *           WORK(KVECA),1,ZERO,WORK(KVECI),1)
C
      CALL DGEMV('N',NRHF(ISYMI),NRHF(ISYML),ONE,DIJ(KOFF12),NTOTI,
     *           WORK(KVECL),1,ONE,WORK(KVECI),1)
C
C---------------------------------------
C     Calculate correlated-frozen block.
C---------------------------------------
C
      DO 130 J = 1,NRHFFR(ISYMJ)
C
         KOFF13 = ILRHSI(ISYMG) + NBAS(ISYMG)*(J - 1) + G
         KOFF14 = ICOFRO(ISYMI,ISYMJ) + NRHF(ISYMI)*(J - 1) + 1
         CALL DAXPY(NRHF(ISYMI),-CMO(KOFF13),WORK(KVECI),1,
     *              D2IJ(KOFF14),1)
C
  130 CONTINUE
C
C-----------------------------------------------------
C     Add contribution to frozen-frozen block from Va.
C-----------------------------------------------------
C
      IF (ISYMG .EQ. ISYMD) THEN
C
         ISYMK = ISYMD
         KVECK = KVECI
C
         KOFF15 = IGLMRH(ISYMG,ISYMK) + G
         FAC = DDOT(NRHF(ISYMK),WORK(KVECK),1,XLAMDP(KOFF15),
     *              NBAS(ISYMG))
C
         DO 140 ISYMM = 1,NSYM
            DO 150 M = 1,NRHFFR(ISYMM)
               KOFF16 = IFROFR(ISYMM,ISYMM) + NRHFFR(ISYMM)*(M - 1) + M
               D2II(KOFF16) = D2II(KOFF16) + TWO*FAC
  150       CONTINUE
  140    CONTINUE
      ENDIF
C
C----------------------------------------------------------------------
C     Calculate Hartree-Fock like contributions to frozen-frozen block.
C----------------------------------------------------------------------
C
      IF (ISYMG .EQ. ISYMD) THEN
C
         KOFF17 = ILRHSI(ISYMG) + G
         KOFF18 = ILRHSI(ISYMD) + IDEL - IBAS(ISYMD)
         FAC = TWO*DDOT(NRHFFR(ISYMG),CMO(KOFF17),NBAS(ISYMG),
     *              CMO(KOFF18),NBAS(ISYMD))
C
         DO 160 ISYMI = 1,NSYM
            DO 170 I = 1,NRHFFR(ISYMI)
               KOFF19 = IFROFR(ISYMI,ISYMI) + NRHFFR(ISYMI)*(I - 1) + I
               D2II(KOFF19) = D2II(KOFF19) + TWO*FAC
  170       CONTINUE
  160    CONTINUE
      ENDIF
C
      ISYMI = ISYMD
      ISYMJ = ISYMG
      D     = IDEL - IBAS(ISYMD)
C
      DO 180 J = 1,NRHFFR(ISYMJ)
         DO 190 I = 1,NRHFFR(ISYMI)
            KOFF20 = IFROFR(ISYMI,ISYMJ) + NRHFFR(ISYMI)*(J - 1) + I
            KOFF21 = ILRHSI(ISYMD) + NBAS(ISYMD)*(I - 1) + D
            KOFF22 = ILRHSI(ISYMG) + NBAS(ISYMG)*(J - 1) + G
            D2II(KOFF20) = D2II(KOFF20) - TWO*CMO(KOFF21)*CMO(KOFF22)
  190    CONTINUE
  180 CONTINUE
C
C---------------------------------
C     Work space allocation three.
C---------------------------------
C
      ISYMB = ISYMG
      ISYMJ = ISYMG
      ISYMA = ISYMB
      ISYMI = ISYMD
C
      KVECA = 1
      KVECB = KVECA + NVIR(ISYMA)
      KVECJ = KVECB + NVIR(ISYMB)
      KEND1 = KVECJ + NRHF(ISYMJ)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space in CC_FD2BL')
      ENDIF
C
      CALL DZERO(WORK(KVECI),KEND1)
C
      KOFF23 = IGLMVI(ISYMG,ISYMB) + G
      KOFF24 = IGLMRH(ISYMG,ISYMJ) + G
      CALL DCOPY(NVIR(ISYMB),XLAMDP(KOFF23),NBAS(ISYMG),WORK(KVECB),1)
      CALL DCOPY(NRHF(ISYMJ),XLAMDP(KOFF24),NBAS(ISYMG),WORK(KVECJ),1)
C
C--------------------------------------
C     Calculate intermediate vector Wa.
C--------------------------------------
C
      KOFF25 = IMATAB(ISYMB,ISYMA) + 1
      KOFF26 = IT1AM(ISYMA,ISYMJ)  + 1
C
      NTOTB  = MAX(NVIR(ISYMB),1)
      NTOTA  = MAX(NVIR(ISYMA),1)
C
      CALL DGEMV('T',NVIR(ISYMB),NVIR(ISYMA),ONE,DAB(KOFF25),NTOTB,
     *           WORK(KVECB),1,ZERO,WORK(KVECA),1)
C
      CALL DGEMV('N',NVIR(ISYMA),NRHF(ISYMJ),ONE,DIA(KOFF26),NTOTA,
     *           WORK(KVECJ),1,ONE,WORK(KVECA),1)
C
C------------------------------------
C     Calculate frozen-virtual block.
C------------------------------------
C
      DO 200 I = 1,NRHFFR(ISYMI)
C
         KOFF27 = ILRHSI(ISYMD) + NBAS(ISYMD)*(I - 1)
     *          + IDEL - IBAS(ISYMD)
         KOFF28 = IT1FRO(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + 1
         CALL DAXPY(NVIR(ISYMA),-CMO(KOFF27),WORK(KVECA),1,
     *              D2IA(KOFF28),1)
C
  200 CONTINUE
C
C--------------------------------
C     Work space allocation four.
C--------------------------------
C
      ISYMA = ISYMG
      ISYMK = ISYMG
      ISYMJ = ISYMA
      ISYMI = ISYMD
C
      KVECJ = 1
      KVECA = KVECJ + NRHF(ISYMJ)
      KVECK = KVECA + NVIR(ISYMA)
      KEND1 = KVECK + NRHF(ISYMK)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space in CC_FD2BL')
      ENDIF
C
      CALL DZERO(WORK(KVECI),KEND1)
C
      KOFF29 = IGLMVI(ISYMG,ISYMA) + G
      KOFF30 = IGLMRH(ISYMG,ISYMk) + G
      CALL DCOPY(NVIR(ISYMA),XLAMDP(KOFF29),NBAS(ISYMG),WORK(KVECA),1)
      CALL DCOPY(NRHF(ISYMK),XLAMDP(KOFF30),NBAS(ISYMG),WORK(KVECK),1)
C
C--------------------------------------
C     Calculate intermediate vector Wj.
C--------------------------------------
C
      KOFF31 = IT1AM(ISYMA,ISYMJ)  + 1
      KOFF32 = IMATIJ(ISYMK,ISYMJ) + 1
C
      NTOTA  = MAX(NVIR(ISYMA),1)
      NTOTK  = MAX(NRHF(ISYMK),1)
C
      CALL DGEMV('T',NVIR(ISYMA),NRHF(ISYMJ),ONE,DAI(KOFF31),NTOTA,
     *           WORK(KVECA),1,ZERO,WORK(KVECJ),1)
C
      CALL DGEMV('T',NRHF(ISYMK),NRHF(ISYMJ),ONE,DIJ(KOFF32),NTOTK,
     *           WORK(KVECK),1,ONE,WORK(KVECJ),1)
C
C---------------------------------------
C     Calculate frozen-correlated block.
C---------------------------------------
C
      DO 210 I = 1,NRHFFR(ISYMI)
C
         KOFF33 = ILRHSI(ISYMD) + NBAS(ISYMD)*(I - 1)
     *          + IDEL - IBAS(ISYMD)
         KOFF34 = ICOFRO(ISYMJ,ISYMI) + NRHF(ISYMJ)*(I - 1) + 1
         CALL DAXPY(NRHF(ISYMJ),-CMO(KOFF33),WORK(KVECJ),1,
     *              D2JI(KOFF34),1)
C
  210 CONTINUE
C
      CALL QEXIT('CC_FD2BL')
C
      RETURN
      END
C  /* Deck ccs_den2 */
      SUBROUTINE CCS_DEN2(D2IJG,CMO,WORK,LWORK,IDEL,ISYMD)
C
C     Written by Asger Halkier 4/6 - 1998.
C
C     Purpose: Calculate the 2 electron CCS (=SCF) density
C              d(pq,gam;del) for a given delta (IDEL).
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION D2IJG(*), CMO(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
      CALL QENTER('CCS_DEN2')
C
      ISYDEN = ISYMD
C
C---------------------------------------------------
C     Calculate the (occ.occ,occ;del) density block.
C---------------------------------------------------
C
      KD2IJK = 1
      KEND1  = KD2IJK + NFRIJK(ISYDEN)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient space for allocation in CCS_DEN2')
      ENDIF
C
      CALL DZERO(WORK(KD2IJK),NFRIJK(ISYDEN))
C
C---------------------------------
C     Calculate the contributions.
C---------------------------------
C
      CALL FIJK01(WORK(KD2IJK),CMO,IDEL,ISYMD)
      CALL FIJK02(WORK(KD2IJK),CMO,IDEL,ISYMD)
C
      CALL F2_PQAO(D2IJG,WORK(KD2IJK),CMO,ISYDEN)
C
      CALL QEXIT('CCS_DEN2')
C
      RETURN
      END
C  /* Deck mp2_den2 */
      SUBROUTINE MP2_DEN2(D2IAG,T2AMM,CMO,WORK,LWORK,IDEL,ISYMD)
C
C     Written by Asger Halkier 5/6 - 1998.
C
C     Purpose: Calculate the non-SCF part of the 2 electron MP2
C              density d(pq,gam;del) for a given delta (IDEL).
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION D2IAG(*), T2AMM(*), CMO(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('MP2_DEN2')
C
      ISYDEN = ISYMD
C
C----------------------------------------------------
C     Calculate terms of the (occ.vir,occ;del) block.
C     Note that we exploit the full permutational
C     symmetry of the 2-electron integrals later on.
C----------------------------------------------------
C
      KD2IAJ = 1
      KEND1  = KD2IAJ + NT2BCD(ISYDEN)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT(
     *        'Insufficient space for fourth allocation in MP2_DEN2')
      ENDIF
C
      CALL DZERO(WORK(KD2IAJ),NT2BCD(ISYDEN))
C
C---------------------------------
C     Calculate the contributions.
C---------------------------------
C
      ISYMT2 = 1
      CALL CC_TI(WORK(KD2IAJ),ISYMD,T2AMM,ISYMT2,CMO,1,
     *           WORK(KEND1),LWRK1,IDEL,ISYMD)
C
      ICON = 2
      CALL CCD2_PQAO(D2IAG,ISYDEN,WORK(KD2IAJ),ISYDEN,CMO,1,ICON)
C
      CALL QEXIT('MP2_DEN2')
C
      RETURN
      END
C  /* Deck ccmp_dao */
      SUBROUTINE CCMP_DAO(AODEN,DENIJ,DENIA,CMO,XCMO,WORK,LWORK,ISDEN)
C
C     Written by Asger Halkier 8/6 - 1998
C
C     Version: 1.0
C
C     Purpose: To transform the blocks of the CCS (1 block) and
C              MP2 (2 blocks) one electron density to AO-basis
C              and add the results!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION AODEN(*), DENIJ(*), DENIA(*)
      DIMENSION CMO(*), XCMO(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <ccfro.h>
C
      CALL QENTER('CCMP_DAO')
C
      DO 100 ISYM1 = 1,NSYM
C
         ISYM2 = MULD2H(ISDEN,ISYM1)
C
         LNEED = NBAS(ISYM1)*NRHFS(ISYM2)
C
         IF (LWORK .LT. LNEED) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', LNEED
            CALL QUIT('Insufficient work space in CCMP_DAO')
         ENDIF
C
         CALL DZERO(WORK,LNEED)
C
C------------------------------------------
C        Transform occupied-occupied block.
C------------------------------------------
C
         KOFF1  = ILRHSI(ISYM1) + 1
         KOFF2  = IFROIJ(ISYM1,ISYM2) + 1
C
         NTOBA1 = MAX(NBAS(ISYM1),1)
         NTORH1 = MAX(NRHFS(ISYM1),1)
C
         CALL DGEMM('N','N',NBAS(ISYM1),NRHFS(ISYM2),NRHFS(ISYM1),
     *              ONE,CMO(KOFF1),NTOBA1,DENIJ(KOFF2),NTORH1,
     *              ZERO,WORK,NTOBA1)
C
         KOFF3  = ILRHSI(ISYM2) + 1
         KOFF4  = IAODIS(ISYM1,ISYM2) + 1
C
         NTOBA1 = MAX(NBAS(ISYM1),1)
         NTOBA2 = MAX(NBAS(ISYM2),1)
C
         CALL DGEMM('N','T',NBAS(ISYM1),NBAS(ISYM2),NRHFS(ISYM2),
     *              ONE,WORK,NTOBA1,CMO(KOFF3),NTOBA2,ONE,
     *              AODEN(KOFF4),NTOBA1)
C
  100 CONTINUE
C
      IF (MP2) THEN
         DO 110 ISYM1 = 1,NSYM
C
            ISYM2 = MULD2H(ISDEN,ISYM1)
C
            LNEED = NBAS(ISYM2)*NRHF(ISYM1)
C
            IF (LWORK .LT. LNEED) THEN
               WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', LNEED
               CALL QUIT('Insufficient work space in CCMP_DAO')
            ENDIF
C
            CALL DZERO(WORK,LNEED)
C
C-------------------------------------------------------------
C        Transform occupied-virtual block (stored transposed).
C-------------------------------------------------------------
C
         KOFF1  = IGLMVI(ISYM2,ISYM2) + 1
         KOFF2  = IT1AM(ISYM2,ISYM1) + 1
C
         NTOBA2 = MAX(NBAS(ISYM2),1)
         NTOVI2 = MAX(NVIR(ISYM2),1)
C
         CALL DGEMM('N','N',NBAS(ISYM2),NRHF(ISYM1),NVIR(ISYM2),
     *              ONE,XCMO(KOFF1),NTOBA2,DENIA(KOFF2),NTOVI2,
     *              ZERO,WORK,NTOBA2)
C
         KOFF3  = IGLMRH(ISYM1,ISYM1) + 1
         KOFF4  = IAODIS(ISYM1,ISYM2) + 1
C
         NTOBA1 = MAX(NBAS(ISYM1),1)
         NTOBA2 = MAX(NBAS(ISYM2),1)
C
         CALL DGEMM('N','T',NBAS(ISYM1),NBAS(ISYM2),NRHF(ISYM1),
     *              ONE,XCMO(KOFF3),NTOBA1,WORK,NTOBA2,ONE,
     *              AODEN(KOFF4),NTOBA1)
C
  110    CONTINUE
      ENDIF
C
      CALL QEXIT('CCMP_DAO')
C
      RETURN
      END
C  /* Deck fijk01 */
      SUBROUTINE FIJK01(D2IJK,CMO,IDEL,ISYMD)
C
C     Written by Asger Halkier 11/6 - 1998
C
C     Purpose: To calculate the first contribution to the
C              CCS 2-electron density including frozen-core
C              contributions - based on DIJK01.
C
#include <implicit.h>
      PARAMETER(FOUR = 4.0D0)
      DIMENSION D2IJK(*), CMO(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
      CALL QENTER('FIJK01')
C
      ISYDEN = ISYMD
      ISYMK  = ISYMD
      ISYMIJ = MULD2H(ISYDEN,ISYMK)
C
C------------------------------------------
C     Calculate adresses and add to result.
C------------------------------------------
C
      DO 100 K = 1,NRHFS(ISYMK)
C
         DO 110 ISYMJ = 1,NSYM
C
            DO 120 J = 1,NRHFS(ISYMJ)
C
               NIJK = IFRIJK(ISYMIJ,ISYMK) + NFROIJ(ISYMIJ)*(K - 1)
     *              + IFROIJ(ISYMJ,ISYMJ) + NRHFS(ISYMJ)*(J - 1) + J
               NDEL = ILRHSI(ISYMK) + NBAS(ISYMD)*(K - 1)
     *              + IDEL - IBAS(ISYMD)
C
               D2IJK(NIJK) = D2IJK(NIJK) + FOUR*CMO(NDEL)
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('FIJK01')
C
      RETURN
      END
C  /* Deck fijk02 */
      SUBROUTINE FIJK02(D2IJK,CMO,IDEL,ISYMD)
C
C     Written by Asger Halkier 11/6 - 1998
C
C     Purpose: To calculate the second contribution to the
C              CCS 2-electron density including frozen-core
C              contributions - based on DIJK02.
C
#include <implicit.h>
      PARAMETER(TWO = 2.0D0)
      DIMENSION D2IJK(*), CMO(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
      CALL QENTER('FIJK02')
C
      ISYDEN = ISYMD
      ISYMI  = ISYMD
C
C------------------------------------------
C     Calculate adresses and add to result.
C------------------------------------------
C
      DO 100 ISYMK = 1,NSYM
C
         ISYMIJ = MULD2H(ISYMI,ISYMK)
C
         DO 110 K = 1,NRHFS(ISYMK)
C
            DO 120 I = 1,NRHFS(ISYMI)
C
               NIJK = IFRIJK(ISYMIJ,ISYMK) + NFROIJ(ISYMIJ)*(K - 1)
     *              + IFROIJ(ISYMI,ISYMK) + NRHFS(ISYMI)*(K - 1) + I
               NDEL = ILRHSI(ISYMI) + NBAS(ISYMD)*(I - 1)
     *              + IDEL - IBAS(ISYMD)
C
               D2IJK(NIJK) = D2IJK(NIJK) - TWO*CMO(NDEL)
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('FIJK02')
C
      RETURN
      END
C  /* Deck f2_pqao */
      SUBROUTINE F2_PQAO(D2IJG,D2IJK,CMO,ISYDEN)
C
C     Written by Asger Halkier 11/6 - 1998
C
C     Purpose: To calculate the backtransformation of the CCS
C              two-electron density d(ij,k;del) to d(ij,gam;del)
C              including the contributions from frozen core orbitals.
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION D2IJG(*), D2IJK(*), CMO(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
      CALL QENTER('F2_PQAO')
C
C----------------------------------
C     Calculate the transformation.
C----------------------------------
C
         DO 100 ISYMIJ = 1,NSYM
C
            ISYMK  = MULD2H(ISYMIJ,ISYDEN)
            ISYMG  = ISYMK
C
            KOFF1  = IFRIJK(ISYMIJ,ISYMK) + 1
            KOFF2  = ILRHSI(ISYMK) + 1
            KOFF3  = IF2IJG(ISYMIJ,ISYMG) + 1
C
            NTOTIJ = MAX(NFROIJ(ISYMIJ),1)
            NTOTG  = MAX(NBAS(ISYMG),1)
C
            CALL DGEMM('N','T',NFROIJ(ISYMIJ),NBAS(ISYMG),NRHFS(ISYMK),
     *                 ONE,D2IJK(KOFF1),NTOTIJ,CMO(KOFF2),NTOTG,
     *                 ONE,D2IJG(KOFF3),NTOTIJ)
C
  100    CONTINUE
C
      CALL QEXIT('F2_PQAO')
C
      RETURN
      END
C  /* Deck cc_d2HFao */
      SUBROUTINE CC_D2HFAO(AODEN,DENIJ,CMO,WORK,LWORK,ISDEN)
C
C     Written by Sonia Coriani 20/11/2008
C     Used to obtain the correlation gradient in Incremental Scheme
C     Version: 1.0, based on CCMP_DAO
C
C     Purpose: To transform the blocks of the CCS/HF (1 block) 
C              two electron density (IJ_GAMMA_DELTA) to AO-basis
C              and subtract the results!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION AODEN(*), DENIJ(*)
      DIMENSION CMO(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <ccfro.h>
      CALL QENTER('CC_D2HFAO')
C
      DO 100 ISYM1 = 1,NSYM
C
         ISYM2 = MULD2H(ISDEN,ISYM1)
C
         LNEED = NBAS(ISYM1)*NRHFS(ISYM2)
C
         IF (LWORK .LT. LNEED) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', LNEED
            CALL QUIT('Insufficient work space in CCMP_DAO')
         ENDIF
C
         CALL DZERO(WORK,LNEED)
C
C------------------------------------------
C        Transform occupied-occupied block.
C------------------------------------------
C
         KOFF1  = ILRHSI(ISYM1) + 1
         KOFF2  = IFROIJ(ISYM1,ISYM2) + 1
C
         NTOBA1 = MAX(NBAS(ISYM1),1)
         NTORH1 = MAX(NRHFS(ISYM1),1)
C
         CALL DGEMM('N','N',NBAS(ISYM1),NRHFS(ISYM2),NRHFS(ISYM1),
     *              ONE,CMO(KOFF1),NTOBA1,DENIJ(KOFF2),NTORH1,
     *              ZERO,WORK,NTOBA1)
C
         KOFF3  = ILRHSI(ISYM2) + 1
         KOFF4  = IAODIS(ISYM1,ISYM2) + 1
C
         NTOBA1 = MAX(NBAS(ISYM1),1)
         NTOBA2 = MAX(NBAS(ISYM2),1)
C
         CALL DGEMM('N','T',NBAS(ISYM1),NBAS(ISYM2),NRHFS(ISYM2),
     *              -ONE,WORK,NTOBA1,CMO(KOFF3),NTOBA2,ONE,
     *              AODEN(KOFF4),NTOBA1)
C
  100 CONTINUE
C
      CALL QEXIT('CC_D2HFAO')
C
      RETURN
      END

