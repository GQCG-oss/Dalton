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
C   /* Deck cc_den_ptfc */
      SUBROUTINE CC_DEN_PTFC(POTNUC,ETAAI,ZKDIA,WORK,LWORK,
     &                     IOPT,IMODEL,LTSTE)
C
C     Written by S. Coriani, based on CC_DEN
C     Fall 2001
C
C     Version: 2.0
C
C     Purpose: 
C
C     For IOPT = 1 calculate the (non corrected) right hand side of the 
C     equation for the zeroth order orbital rotation lagrange multiplier 
C     response (Zeta-Kappa-0)_ai and return it in the array ETAAI
C     For test purposes ...
C     The (P) densities are read in from file
C
C     For IOPT = 2 calculate ALSO the "diagonal" ZKAPPA_pp (returned in
C     ZKDIA), as well as the "diagonal" RHS eta_ij, eta_ab (local arrays)
C     and the corresponding ZKAPPA_ij, ZKAPPA_ab (also returned in ZKDIA)
C     and use them to correct ETAAI (I moved inside the call to ETACORR!!!)
C
C     ZKAPPA_ij, ZKAPPA_ab are simply equal to eta_ij eta_ab scaled by
C     orbital energies.....
C
C     Current models: CCSD (IMODEL=0, test purposes), CCSD(T) (IMODEL=1)
C
C     LTSTE = .true. test densities via Energy calculation
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxorb.h>
#include <maxash.h>
#include <mxcent.h>
#include <aovec.h>
#include <iratdef.h>
      PARAMETER (ZERO = 0.0D0,HALF=0.5D0,ONE = 1.0D0,TWO = 2.0D0)
      PARAMETER (TRE = 3.0D0, FOUR = 4.0D0)
      DIMENSION INDEXA(MXCORB)
      DIMENSION ETAAI(*), ZKDIA(*), WORK(LWORK)
      LOGICAL LTSTE, LETAFI, LETIFJ
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
CTEST
C#include <ccfop.h>

      CHARACTER MODEL*10
      CHARACTER NAME1*8
      CHARACTER NAME2*8
      CHARACTER*5 FNKAPAB, FNKAPIJ, FNDPTIA

      LOGICAL LOCDBG
      PARAMETER (LOCDBG=.FALSE.)
C
      CALL QENTER('CC_DEN_PTFC')

C
      IF (FROIMP) THEN
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
      IF (IOPT .LE. 2) THEN
         CALL HEADER('Constructing RHS for CCSD(T)-kapbar-0',-1)
         CALL HEADER('Inside ccpt_den()',-1)
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
      LUNGO = 2*NT1AMX    + NMATIJ(1)   + NMATAB(1)
     *          + 2*NCOFRO(1) + 2*NT1FRO(1)
C
C-----------------------------------
C     Initial work space allocation.
C-----------------------------------
C
      IF (LTSTE) THEN
        KD1AOB = 1
        KSTART = KD1AOB + N2BST(1)
      ELSE 
        KSTART = 1
      END IF

      KZ2AM  = KSTART
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
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for initial allocation '//
     &             'in CC_DEN_PTFC')
      ENDIF
C
C----------------------------------------
C     Read zero-th order zeta amplitudes.
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
C     Read zero-th order cluster amplitudes.
C-------------------------------------------
C
      IOPTRW = 3
      CALL CC_RDRSP('R0',0,1,IOPTRW,MODEL,WORK(KT1AM),WORK(KT2AM))
C
C----------------------------------
C     Calculate the lambda matrices.
C     Used in CCSD-like part
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
C----------------------------------------------
C     Work space allocation one. CCSD-like part
C     Note that D(ai) = ZETA(ai), and both 
C     D(ia) and h(ia) are stored transposed!
C----------------------------------------------
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
C
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

      KCMO  = KEND1
      KEND1 = KCMO + NLAMDS
      LWRK1 = LWORK - KEND1

      IF (FROIMP) THEN
         KCMOF = KEND1
         KEND1 = KCMOF + NLAMDS
         LWRK1 = LWORK - KEND1
      ENDIF
C
      IF (LWRK1 .LT. 0) THEN
        WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
        CALL QUIT('Insuff. memory for allocation 1 CC_DEN_PTFC')
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
C-------------------------------------------------
C     Get the non-frozen MO coefficient matrix
C     for (T) part and reorder.
C-------------------------------------------------
C
      CALL CC_GET_CMO(WORK(KCMO))
      CALL CMO_REORDER(WORK(KCMO),WORK(KEND1),LWRK1)
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
C------------------------------------------------------------------------
C     Construct three remaining blocks of CCSD-like one electron density.
C------------------------------------------------------------------------
C
      CALL DCOPY(NMATAB(1),WORK(KYMAT),1,WORK(KONEAB),1)
      CALL CC_EITR(WORK(KONEAB),WORK(KONEIJ),WORK(KEND1),LWRK1,1)
      CALL DIJGEN(WORK(KONEIJ),WORK(KXMAT))
      CALL DIAGEN(WORK(KONEIA),WORK(KT2AMT),WORK(KONEAI))

      IF (LOCDBG) THEN
         DIJNO = DDOT(NMATIJ(1),WORK(KONEIJ),1,WORK(KONEIJ),1)
         DAINO = DDOT(NT1AMX,WORK(KONEAI),1,WORK(KONEAI),1)
         DIANO = DDOT(NT1AMX,WORK(KONEIA),1,WORK(KONEIA),1)
         DABNO = DDOT(NMATAB(1),WORK(KONEAB),1,WORK(KONEAB),1)
         WRITE(LUPRI,*) 'CC_DEN_PTFC: IOPT =  ', IOPT
         WRITE(LUPRI,*) 
     &   'Norm of CCSD-like one electron densities in MO-basis:'
         WRITE(LUPRI,*) DIJNO, DAINO, DIANO, DABNO
      ENDIF
C
C---------------------------------
C     Read one-electron integrals.
C---------------------------------
C
      CALL CCRHS_ONEAO(WORK(KONINT),WORK(KEND1),LWRK1)

      IF (LTSTE) THEN

         CALL DZERO(WORK(KD1AOB),N2BST(1))
         ISDEN = 1
         CALL CC_DENAO(WORK(KD1AOB),ISDEN,WORK(KONEAI),WORK(KONEAB),
     *                 WORK(KONEIJ),WORK(KONEIA),ISDEN,WORK(KLAMDP),1,
     *                 WORK(KLAMDH),1,WORK(KEND1),LWRK1)
C
         IF (FROIMP) THEN
            MODEL = 'DUMMY'
            CALL CC_FCD1AO(WORK(KD1AOB),WORK(KEND1),LWRK1,MODEL)
         END IF

         ECCSD1 = DDOT(N2BST(ISYMOP),WORK(KD1AOB),1,WORK(KONINT),1)
         ECCSD2 = ZERO

      END IF
C
C---------------------------------------------------------
C        Ove 24-20-1996
C        Read one-electron integrals into Fock-matrix for
C        finite field.
C---------------------------------------------------------
C
        DO 13 IF = 1, NFIELD
c          DTIME  = SECOND()
           FF = EFIELD(IF)
           CALL CC_ONEP(WORK(KONINT),WORK(KEND1),LWRK1,FF,1,LFIELD(IF))
c          DTIME  = DTIME - SECOND()
c          TIMFCK = TIMFCK + DTIME
 13     CONTINUE
C
C--------------------------------------------------
C       Transform one electron integrals to MO-basis.
C--------------------------------------------------
C
        ISYM = 1
        CALL CCDINTMO(WORK(KINTIJ),WORK(KINTIA),WORK(KINTAB),
     *                 WORK(KINTAI),WORK(KONINT),WORK(KLAMDP),
     *                 WORK(KLAMDH),WORK(KEND1),LWRK1,ISYM)
C
        IF (LOCDBG) THEN
 
           HIJNO = DDOT(NMATIJ(1),WORK(KINTIJ),1,WORK(KINTIJ),1)
           HAINO = DDOT(NT1AMX,WORK(KINTAI),1,WORK(KINTAI),1)
           HIANO = DDOT(NT1AMX,WORK(KINTIA),1,WORK(KINTIA),1)
           HABNO = DDOT(NMATAB(1),WORK(KINTAB),1,WORK(KINTAB),1)
           WRITE(LUPRI,*) 'CC_DEN_PTFC: IOPT = ', IOPT
           WRITE(LUPRI,*) 
     &     'UFFA: Norm of Lambda one electron integrals in MO-basis:'
           WRITE(LUPRI,*) HIJNO, HAINO, HIANO, HABNO
 
           TAMNOR = DDOT(NT1AM(1),WORK(KT1AM),1,WORK(KT1AM),1)
           WRITE(LUPRI,*) 'Norm of T1AM amplitudes', TAMNOR
 
        ENDIF
C
        IF (FROIMP) THEN
C
            ISYM = 1
            !obtain integrals with frozen indices
            ! h_Ij h_jI h_aJ h_Ia h_IJ
            !
            CALL CCIFROMO(WORK(KFROIJ),WORK(KFROJI),WORK(KFROAI),
     *                    WORK(KFROIA),WORK(KFROII),WORK(KONINT),
     *                    WORK(KLAMDP),WORK(KLAMDH),WORK(KCMOF),
     *                    WORK(KEND1),LWRK1,ISYM)
C
            !calculare D_II = 2 delta_IJ
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
        ISYM = 1

        IF (IOPT.EQ.2) THEN

          KOFFIJ = 1
          CALL CC2_ETIJ(ZKDIA(KOFFIJ),WORK(KINTIJ),WORK(KINTAI),
     &                WORK(KINTIA),WORK(KINTAB),WORK(KONEIJ),
     &                WORK(KONEAI),WORK(KONEIA),WORK(KONEAB),
     &                WORK(KT1AM),WORK(KEND1),LWRK1,ISYM)


          KOFFAB = NMATIJ(1) + 1
          CALL CC2_ETAB(ZKDIA(KOFFAB),WORK(KINTIJ),WORK(KINTAI),
     &                WORK(KINTIA),WORK(KINTAB),WORK(KONEIJ),
     &                WORK(KONEAI),WORK(KONEIA),WORK(KONEAB),
     &                WORK(KT1AM),WORK(KEND1),LWRK1,ISYM)


        END IF

C------------------------------------------------------------

        CALL DZERO(ETAAI,NALLAI(1))
        CALL CCDENZK0(ETAAI,WORK(KINTIJ),WORK(KINTAI),WORK(KINTIA),
     *                WORK(KINTAB),WORK(KONEIJ),WORK(KONEAI),
     *                WORK(KONEIA),WORK(KONEAB),WORK(KINIJT),
     *                WORK(KINABT),WORK(KD1IJT),WORK(KD1ABT),
     *                WORK(KT1AM),WORK(KEND1),LWRK1,ISYM)
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
            !
            !eta_iJ (one el)
            !
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
            !
            !eta_aI
            !
            CALL CC_ETAIF(ZKDIA(KOFRES),WORK(KONEAB),WORK(KONEAI),
     *                    WORK(KONEIA),WORK(KD1IJT),WORK(KD1ABT),
     *                    WORK(KFD1II),SK1,SK2,SK3,SK4,WORK(KINTIJ),
     *                    WORK(KINTAB),WORK(KINTAI),WORK(KINTIA),
     *                    WORK(KINABT),WORK(KFROIJ),WORK(KFROJI),
     *                    WORK(KFROAI),WORK(KFROIA),WORK(KFROII),
     *                    WORK(KT1AM),WORK(KEND1),LWRK1,ISYM,ICON)

        ENDIF

C
C-----------------------------------------------------------------------
C     Transform now one-electron integrals with regular CMO coefficients
C-----------------------------------------------------------------------
C
      IF (IMODEL.EQ.1) THEN
        ISYM = 1
        CALL CCDINTMO(WORK(KINTIJ),WORK(KINTIA),WORK(KINTAB),
     *                 WORK(KINTAI),WORK(KONINT),WORK(KCMO),
     *                 WORK(KCMO),WORK(KEND1),LWRK1,ISYM)
C
        IF (LOCDBG) THEN
 
           HIJNO = DDOT(NMATIJ(1),WORK(KINTIJ),1,WORK(KINTIJ),1)
           HAINO = DDOT(NT1AMX,WORK(KINTAI),1,WORK(KINTAI),1)
           HIANO = DDOT(NT1AMX,WORK(KINTIA),1,WORK(KINTIA),1)
           HABNO = DDOT(NMATAB(1),WORK(KINTAB),1,WORK(KINTAB),1)
           WRITE(LUPRI,*) 'CC_DEN_PTFC: IOPT = ', IOPT
           WRITE(LUPRI,*) 
     &     'Norm of CMO one electron integrals in MO-basis:'
           WRITE(LUPRI,*) HIJNO, HAINO, HIANO, HABNO
 
        ENDIF
C
C--------------------------------------------------
C     Read in D^(T)_ia from file, contract with h_pq
C     and add to result 
C--------------------------------------------------

        KONEIA_PT = KEND1
        KEND1     = KONEIA_PT + NT1AMX
        LWRK1     = LWORK - KEND1
C
        LUPTIA = -1
        FNDPTIA = 'DPTIA'
        CALL WOPEN2(LUPTIA,FNDPTIA,64,0)
C
        IOFF = 1
        CALL GETWA2(LUPTIA,FNDPTIA,WORK(KONEIA_PT),IOFF,NT1AMX)
        CALL WCLOSE2(LUPTIA,FNDPTIA,'KEEP')

        if (locdbg) then
           write(lupri,*)'Norm of (T) D_ia:', 
     &     ddot(NT1AMX,WORK(KONEIA_PT),1,WORK(KONEIA_PT),1)
        end if

        IF (LTSTE) THEN
           !
           !D(T)_ia contribution to the energy: sum_ia D_ia h_ia
           !
           EN1PT = DDOT(NT1AMX,WORK(KONEIA_PT),1,WORK(KINTIA),1)

        END IF

C
C------------------------------------------------------------
C      Add one electron contribution to Zeta-kappa-0.
C      from the (T) one electron density...
C------------------------------------------------------------
C
        ISYM = 1

        IF (IOPT.EQ.2) THEN

           KOFFIJ = 1
           KOFFAB = 1 + NMATIJ(1)
           CALL CCPT_ETARS_1E(ZKDIA(KOFFIJ),ZKDIA(KOFFAB),
     &                        WORK(KINTIJ),WORK(KINTAI),
     &                WORK(KINTIA),WORK(KINTAB),WORK(KONEIA_PT),
     &                WORK(KEND1),LWRK1,ISYM)


        END IF
C
C------------------------------------------------------------
C Calculate (T) contribution to ETA_ai and transfer outside
C------------------------------------------------------------
C
        CALL CCPT_ETAAI_1E(ETAAI,WORK(KINTIJ),WORK(KINTAI),
     &                   WORK(KINTIA),WORK(KINTAB),
     &                   WORK(KONEIA_PT),WORK(KEND1),
     &                   LWRK1,ISYM)

        LUNGO = NMATIJ(1)   + NMATAB(1) + 2*NT1AMX  
     *                      + 2*NCOFRO(1) + 2*NT1FRO(1)
        XZKDIA = DDOT(LUNGO,ZKDIA,1,ZKDIA,1)
        WRITE(LUPRI,*) 
     &     'CC_DEN_PTFC: ZKDIA before TWO electron part', XZKDIA

        IF (FROIMP) THEN
C
            ISYM = 1
            !obtain integrals with frozen indices on regular MO basis
            ! h_Ij h_jI h_aJ h_Ia h_IJ
            !
            CALL CCIFROMO(WORK(KFROIJ),WORK(KFROJI),WORK(KFROAI),
     *                    WORK(KFROIA),WORK(KFROII),WORK(KONINT),
     *                    WORK(KCMO),WORK(KCMO),WORK(KCMOF),
     *                    WORK(KEND1),LWRK1,ISYM)
           IF (LOCDBG) THEN
 
             HIJFN = DDOT(NCOFRO(1),WORK(KFROIJ),1,WORK(KFROIJ),1)
             HJIFN = DDOT(NCOFRO(1),WORK(KFROJI),1,WORK(KFROJI),1)
             HAIFN = DDOT(NT1FRO(1),WORK(KFROAI),1,WORK(KFROAI),1)
             HIAFN = DDOT(NT1FRO(1),WORK(KFROIA),1,WORK(KFROIA),1)
             WRITE(LUPRI,*) 'CC_DEN_PTFC: IOPT = ', IOPT
             WRITE(LUPRI,*) 
     &       'Norm of CMO one ele fro integrals in MO-basis:'
             WRITE(LUPRI,*) HIJFN, HJIFN,HAIFN, HIAFN
 
           ENDIF
 
            KOFRE1 = 2*NT1AMX + NMATIJ(1) + NMATAB(1) + 1
            KOFRE2 = 2*NT1AMX + NMATIJ(1) + NMATAB(1) + 2*NT1FRO(1) + 1
 
            CALL CC_ETFRO1E(ZKDIA(KOFRE1),ZKDIA(KOFRE2),
     &                        WORK(KONEIA_PT),
     &                        WORK(KFROIJ),WORK(KFROJI),
     &                        WORK(KFROAI),WORK(KFROIA),
     &                        WORK(KEND1),LWRK1,ISYM)
C
         END IF
      END IF

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
C--------------------------------------------
C     Start the loop over integrals.
C     Salva tutto quanto definito fino ad ora
C--------------------------------------------
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
                  CALL QUIT('Insufficient core in CC_DEN_PTFC')
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
     &                      '2 in CC_DEN_PTFC')
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
C-----------------------------------------------------------------------------
C              Calculate the CCSD-like two electron density d(pq,gamma;delta).
C-----------------------------------------------------------------------------
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
     &                      '3 in CC_DEN_PTFC')
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
               IF (FROIMP) THEN

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
     &                         'in CC_DEN_PTFC')
                  ENDIF
C 
                  CALL DZERO(WORK(KSCRAI),NOFROO(ISYDIS))
                  CALL DZERO(WORK(KSCAIF),NCOFRF(ISYDIS))
C
C-------------------------------------------------------------------------
C                 Transform one index in the integral batch to correlated.
C-------------------------------------------------------------------------
C
                  !(alp-bet,i,delta)
                  CALL CCTRBT(WORK(KXINT),WORK(KDSRHF),WORK(KLAMDP),
     *                        ISYMOP,WORK(KEND3),LWRK3,ISYDIS)
C
C---------------------------------------------------------------------
C                 Transform one index in the integral batch to frozen.
C---------------------------------------------------------------------
C
                  !(alp-bet,i,delta)
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
C-------------------------------------------------------------------------
C                    Work space allocation four.
C                    Note: d2aob is only used for the ltest case!!!!!!!!!!
C-------------------------------------------------------------------------
C
                     KINTAO = KEND3
                     KD2AOB = KINTAO + N2BST(ISYMPQ)
                     KEND4  = KD2AOB + N2BST(ISYMPQ)
                     LWRK4  = LWORK  - KEND4
C
                     IF (LWRK4 .LT. 0) THEN
                        WRITE(LUPRI,*) 'Available:', LWORK
                        WRITE(LUPRI,*) 'Needed:', KEND4
                        CALL QUIT('Insufficient  space in CC_DEN_PTFC')
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
     *                       'Insufficient work space in CC_DEN_PTFC')
                        ENDIF
C
                        CALL DZERO(WORK(KFD2IJ),NCOFRO(ISYMPQ))
                        CALL DZERO(WORK(KFD2JI),NCOFRO(ISYMPQ))
                        CALL DZERO(WORK(KFD2AI),NT1FRO(ISYMPQ))
                        CALL DZERO(WORK(KFD2IA),NT1FRO(ISYMPQ))
                        CALL DZERO(WORK(KFD2II),NFROFR(ISYMPQ))
C
C              To calculate the contributions to d(pq,gam;del)
C              where at least one of the two indices p & q is frozen.
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
C              ! calculate the contributions to D2AO from d(pq,gam;del)
C              ! where at least one of the two indices p & q is frozen

                        CALL CC_FD2AO(WORK(KD2AOB),WORK(KFD2II),
     *                                WORK(KFD2IJ),WORK(KFD2JI),
     *                                WORK(KFD2AI),WORK(KFD2IA),
     *                                WORK(KCMOF),WORK(KLAMDH),
     *                          WORK(KLAMDP),WORK(KEND4),LWRK4,
     *                                ISYMPQ)
C
C     Purpose: To calculate the contributions to d(pq,gam;del) where
C              gamma has been backtransformed from a frozen index.
C
                        CALL CC_D2GAF(WORK(KD2GIJ),WORK(KD2GAB),
     *                                WORK(KD2GAI),WORK(KD2GIA),
     *                                WORK(KONEIJ),WORK(KONEAB),
     *                                WORK(KONEAI),WORK(KONEIA),
     *                           WORK(KCMOF),IDEL,ISYMD,G,ISYMG)
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
C---------------------------------------------------------------------------
C                    If energy test backtransform density fully to AO basis.
C---------------------------------------------------------------------------
C
                     IF (LTSTE) THEN
C
                        CALL CC_DENAO(WORK(KD2AOB),ISYMPQ,
     *                                WORK(KD2GAI),WORK(KD2GAB),
     *                                WORK(KD2GIJ),WORK(KD2GIA),ISYMPQ,
     *                                WORK(KLAMDP),1,WORK(KLAMDH),1,
     *                                WORK(KEND4),LWRK4)
C
C---------------------------------------------------------------------
C                       Add relaxation terms to set up effective density.
C---------------------------------------------------------------------
C
!                        IF (IOPT .EQ. 3) THEN
C
!                           ICON = 1
!                           CALL CC_D2EFF(WORK(KD2AOB),G,ISYMG,IDEL,
!     *                                   ISYMD,WORK(KKABAO),
!     *                                   WORK(KDHFAO),ICON)
C
!                        ENDIF
C
C----------------------------------------------------------------------
C                    Calculate the 2 e- density contribution to E-ccsd.
C----------------------------------------------------------------------
C
                        ECCSD2 = ECCSD2 + HALF*DDOT(N2BST(ISYMPQ),
     *                                    WORK(KD2AOB),1,WORK(KINTAO),1)
C
                      end if
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
     &                               'in CC_DEN_PTFC')
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
C Prepare integrals g_pq (gam,del) where one index is frozen
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
                        IF (IOPT.EQ.2) THEN

                          KOFFIJ = 1
                          CALL CC2_ETIJ(ZKDIA(KOFFIJ),
     &                                  WORK(KIJINT),WORK(KAIINT),
     &                                  WORK(KIAINT),WORK(KABINT),
     &                                  WORK(KD2GIJ),WORK(KD2GAI),
     &                                  WORK(KD2GIA),WORK(KD2GAB),
     &                                  WORK(KT1AM),WORK(KEND5),LWRK5,
     &                                  ISYM)

                          KOFFAB = NMATIJ(1) + 1
                          CALL CC2_ETAB(ZKDIA(KOFFAB),
     &                                  WORK(KIJINT),WORK(KAIINT),
     *                                  WORK(KIAINT),WORK(KABINT),
     *                                  WORK(KD2GIJ),WORK(KD2GAI),
     *                                  WORK(KD2GIA),WORK(KD2GAB),
     *                                  WORK(KT1AM),WORK(KEND5),LWRK5,
     *                                  ISYM)
                        END IF

                        CALL CCDENZK0(ETAAI,WORK(KIJINT),WORK(KAIINT),
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
                           !
                           ! contributions to eta_ai from loop over frozen
                           !
                           CALL CCFRETAI(ETAAI,WORK(KIJFRO),
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
                           !
                           ! eta_iJ
                           !
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
                           !
                           ! eta_aI
                           !
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
c                          !
c                          ! eta_ab contribs from loop over frozen I
c                          !
                           CALL CCFRETAB(ZKDIA,WORK(KIJFRO),
     *                                   WORK(KJIFRO),WORK(KAIFRO),
     *                                   WORK(KIAFRO),WORK(KFD2IJ),
     *                                   WORK(KFD2JI),WORK(KFD2AI),
     *                                   WORK(KFD2IA),WORK(KT1AM),
     *                                   WORK(KEND5),LWRK5,ISYM)
                           !
                           ! eta_ij contribs from loop over frozen L
                           !
                           CALL CCFRETIJ(ZKDIA,WORK(KIJFRO),
     *                                   WORK(KJIFRO),WORK(KAIFRO),
     *                                   WORK(KIAFRO),WORK(KFD2IJ),
     *                                   WORK(KFD2JI),WORK(KFD2AI),
     *                                   WORK(KFD2IA),WORK(KT1AM),
     *                                   WORK(KEND5),LWRK5,ISYM)
c
C
                        ENDIF  !froimp
!                     END IF    !ltste
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
C--------------------------------------------
C     Add contribution from 2-e (T) densities
C--------------------------------------------
C
      if (locdbg) then
          KOFFIJ = 1
          KOFFAB = 1 + NMATIJ(1)
          KOFAFI = 2*NT1AMX + NMATIJ(1) + NMATAB(1) + 1
          KOFIFJ = 2*NT1AMX + NMATIJ(1) + NMATAB(1) + 2*NT1FRO(1) + 1
          write(lupri,*) '                         '
          write(lupri,*) 'Before call to CC_DEN_PT2'
          xtest=ddot(nmatij(1),zkdia(koffij),1,zkdia(koffij),1)
          write(lupri,*) 'Norm of ETAIJ: ', xtest
          xtest=ddot(2*NCOFRO(1),zkdia(kofifj),1,zkdia(kofifj),1)
          write(lupri,*) 'Norm of ETIFJ: ', xtest
          xtest=ddot(nmatab(1),zkdia(koffab),1,zkdia(koffab),1)
          write(lupri,*) 'Norm of ETAAB: ', xtest
          xtest=ddot(2*nt1amx,etaai(1),1,etaai(1),1)
          write(lupri,*) 'Norm of ETAAI: ', xtest
          xtest=ddot(2*nt1fro(1),zkdia(kofafi),1,zkdia(kofafi),1)
          write(lupri,*) 'Norm of ETAFI: ', xtest
      end if

      IF (IMODEL.EQ.1) THEN

         IOPT1 = 2
         en2pt = zero

         if (.false.) then
         KOFFIJ = 1
         KOFFAB = 1 + NMATIJ(1)
         CALL CC_DEN_PT2(ETAAI,ZKDIA(KOFFIJ),
     &                   ZKDIA(KOFFAB),WORK(KEND1),LWRK1,IOPT1,
     &                   ltste,en2pt)
         else

         CALL CCPT_DEN2FC(ETAAI,ZKDIA,
     &                     WORK(KONEIA_PT),
     &                     WORK(KEND1),LWRK1,
     &                     IOPT1,LTSTE,EN2PT)
         end if
        
      END IF
      if (locdbg) then
          KOFFIJ = 1
          KOFFAB = 1 + NMATIJ(1)
          KOFIFJ = 2*NT1AMX + NMATIJ(1) + NMATAB(1) + 2*NT1FRO(1) + 1
          KOFAFI = 2*NT1AMX + NMATIJ(1) + NMATAB(1) + 1
          write(lupri,*) '                         '
          write(lupri,*) 'After call to CC_DEN_PT2'
          xtest=ddot(nmatij(1),zkdia(koffij),1,zkdia(koffij),1)
          write(lupri,*) 'Norm of ETAIJ: ', xtest
          xtest=ddot(NCOFRO(1),zkdia(kofifj),1,zkdia(kofifj),1)
          write(lupri,*) 'Norm of ETIFJ: ', xtest
          xtest=ddot(nmatab(1),zkdia(koffab),1,zkdia(koffab),1)
          write(lupri,*) 'Norm of ETAAB: ', xtest
          xtest=ddot(nt1amx,etaai(1),1,etaai(1),1)
          write(lupri,*) 'Norm of ETAAI: ', xtest
          xtest=ddot(nt1fro(1),zkdia(kofafi),1,zkdia(kofafi),1)
          write(lupri,*) 'Norm of ETAFI: ', xtest
      end if
C
C------------------------------------------------------------------------
C Calculate the ZK0(ij) and ZK0(ab) blocks (to be used to correct eta_ai)
C from eta_ij/e_i-e_j and eta_ab/e_a-e_b
C------------------------------------------------------------------------
C
      IF (IOPT.EQ.2) THEN

         CALL CCSD_ZKBLO(ZKDIA,WORK(KEND1),LWRK1)

C
C------------------------------------------------------------------------
C Add the diagonal ZK0(ii) and ZK0(aa) elements in the proper places
C------------------------------------------------------------------------
C
         if (.true.) then
         IF (IMODEL.EQ.1) THEN

           LUKAPAB  = -1
           LUKAPIJ  = -1
           FNKAPAB  = 'KAPAB'
           FNKAPIJ  = 'KAPIJ'


           CALL WOPEN2(LUKAPAB,FNKAPAB,64,0)
           CALL WOPEN2(LUKAPIJ,FNKAPIJ,64,0)
C
           KKAPII = KEND1
           KEND1  = KKAPII + NRHFT
           LWRK1  = LWORK  - KEND1

           IF (NRHFT .GT. 0) THEN
              IOFF = 1
              CALL GETWA2(LUKAPIJ,FNKAPIJ,WORK(KKAPII),IOFF,NRHFT)
           ENDIF

           DO  ISYMI = 1, NSYM
             DO I = 1, NRHF(ISYMI)
               NI  = IRHF(ISYMI) + I
               NII = IMATIJ(ISYMI,ISYMI) + NRHF(ISYMI)*(I-1) + I
               ZKDIA(NII) = WORK(KKAPII+NI-1)
             END DO
           END DO

           KKAPAA = KEND1
           KEND1  = KKAPAA + NVIRT
           LWRK1  = LWORK  - KEND1
C
           IF (NVIRT .GT. 0) THEN
              IOFF = 1
              CALL GETWA2(LUKAPAB,FNKAPAB,WORK(KKAPAA),IOFF,NVIRT)
           ENDIF

! Bemaerk, IVIR(ISYM) initialized to NRHFT

           DO  ISYMA = 1, NSYM
             DO A = 1, NVIR(ISYMA)
               NA  = IVIR(ISYMA) + A - NRHFT
               NAA = IMATAB(ISYMA,ISYMA) + NVIR(ISYMA)*(A-1) + A
               ZKDIA(NMATIJ(1) + NAA) = WORK(KKAPAA+NA-1)
             END DO
           END DO

           if (locdbg) then
              XZKDIA = DDOT(LUNGO,ZKDIA,1,ZKDIA,1)
              WRITE(LUPRI,*) 'CC_DEN_PTFC: ZKDIA after ii, aa ', XZKDIA

              XZKDIA = DDOT(LUNGO,ZKDIA,1,ZKDIA,1)
              WRITE(LUPRI,*) 'CC_DEN_PTFC: ZKDIA after ii, aa ', XZKDIA
           end if

           CALL WCLOSE2(LUKAPAB,FNKAPAB,'KEEP')
           CALL WCLOSE2(LUKAPIJ,FNKAPIJ,'KEEP')


           IF (LTSTE) THEN

             !multiply by epsilon_p and sum over p to get the energy

             KFOCKDIA = KEND1
             KEND1    = KFOCKDIA + NORBTS
             LWRK1    = LWORK-KEND1
C
C-------------------------------------
C     Read canonical orbital energies.
C-------------------------------------
C
             CALL DZERO(WORK(KFOCKDIA),NORBTS)
             CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
             REWIND (LUSIFC)
C
             CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
             READ (LUSIFC)
             READ (LUSIFC) (WORK(KFOCKDIA + I - 1), I = 1,NORBTS)
C
             CALL GPCLOSE(LUSIFC,'KEEP')
C
C----------------------------------------------------------------
C     Change symmetry ordering of the canonical orbital energies.
C----------------------------------------------------------------
C
             IF (FROIMP)
     *           CALL CCSD_DELFRO(WORK(KFOCKDIA),WORK(KEND1),LWRK1)
C
                 CALL FOCK_REORDER(WORK(KFOCKDIA),WORK(KEND1),LWRK1)
C
C--------------------------------------------
C     Calculate sum_p kappabar_pp * epsilon_p
C     Occupied block:
C--------------------------------------------
C
            SKAPEP = DDOT(NORBT,WORK(KKAPII),1,WORK(KFOCKDIA),1)
             END IF
           END IF
c
      END IF
      end if
C
C------------------------------------------------
C Correct Eta_ai with ZK0(ij) and ZK0(ab) blocks 
C------------------------------------------------
C
      IF ((IOPT.EQ.2).OR.(FROIMP)) THEN
         IF (FROIMP) THEN

             KOFIFJ = 2*NT1AMX + NMATIJ(1) + NMATAB(1) + 2*NT1FRO(1) + 1
             !
             !calculate kappabar_iJ = eta_iJ/e_i-e_J
             !
             CALL CC_ZKFCB(ZKDIA(KOFIFJ),WORK(KEND1),LWRK1)
             !xnorm = ddot(ncofro(1),zkdia(kofifj),1,zkdia(kofifj),1)
             !write(lupri,*) ' Norm of ZKbar_iJ : ', xnorm
         END IF

         !if (.true.) then 
           write(lupri,*)'CC_DEN_PTFC using CCETACOR'
           CALL CCETACOR(ETAAI,ZKDIA,WORK(KEND1),LWRK1)
         !else
         !  write(lupri,*)'CC_DEN_PTFC using CCETACORNEW'
         !  CALL CCETACORNEW(ETAAI,ZKDIA,WORK(KEND1),LWRK1)
         !end if
 
          KOFFIJ = 1
          KOFFAB = NMATIJ(1) + 1
          KOFAFI = 2*NT1AMX + NMATIJ(1) + NMATAB(1) + 1
          KOFIFJ = 2*NT1AMX + NMATIJ(1) + NMATAB(1) + 2*NT1FRO(1) + 1

          !write(lupri,*) 'After call to CCETACORNEW'
          !xtest=ddot(nmatij(1),zkdia(koffij),1,zkdia(koffij),1)
          !write(lupri,*) 'Norm of ETAIJ: ', xtest
          !xtest=ddot(NCOFRO(1),zkdia(kofifj),1,zkdia(kofifj),1)
          !write(lupri,*) 'Norm of ETIFJ: ', xtest
          !xtest=ddot(nmatab(1),zkdia(koffab),1,zkdia(koffab),1)
          !write(lupri,*) 'Norm of ETAAB: ', xtest
          !xtest=ddot(nt1amx,etaai(1),1,etaai(1),1)
          !write(lupri,*) 'Norm of ETAAI: ', xtest
          !xtest=ddot(nt1fro(1),zkdia(kofafi),1,zkdia(kofafi),1)
          !write(lupri,*) 'Norm of ETAFI: ', xtest

      END IF

      XZKDIA = DDOT(LUNGO,ZKDIA,1,ZKDIA,1)
      WRITE(LUPRI,*) 'CC_DEN_PTFC: ZKDIA before CCETACOR', XZKDIA
C
C------------------------------------------------
C     Write out frozen block of eta-kappa-bar-0.
C     eta_iJ
C------------------------------------------------
C
      IF (((IPRINT .GT. 9).OR.(LOCDBG)) .AND. (FROIMP)) THEN
         CALL AROUND('Eta-kappa-bar-0 correlated-frozen block')
         KOFRES = NMATIJ(1) + NMATAB(1) + 2*NT1AMX + 2*NT1FRO(1) + 1
         ZKAPF1 = DDOT(NCOFRO(1),ZKDIA(KOFRES),1,
     *                 ZKDIA(KOFRES),1)
         ZKAPFR = ZKAPF1**0.5
         WRITE(LUPRI,*) 'Norm of correlated-frozen block:', ZKAPFR
      ENDIF
C
C-----------------------------------------------------------
C     Calculate the correlated-frozen blocks of kappa-bar-0.
C-----------------------------------------------------------
C
      IF (FROIMP) THEN

         !KOFIFJ = 2*NT1AMX + NMATIJ(1) + NMATAB(1) + 2*NT1FRO(1) + 1
         !
         !calculate kappabar_iJ = eta_iJ/e_i-e_J
         !
         !CALL CC_ZKFCB(ZKDIA(KOFIFJ),WORK(KEND1),LWRK1)
         !xnorm = ddot(ncofro(1),zkdia(kofifj),1,zkdia(kofifj),1)
         !write(lupri,*) ' Norm of ZKbar_iJ : ', xnorm
C
C----------------------------------------------------------------
C        Calculate correction terms from correlated-frozen block 
C        eta_iJ to both eta_ai and eta_aI
C----------------------------------------------------------------
C
         KOFAFI = 2*NT1AMX + NMATIJ(1) + NMATAB(1) + 1
         CALL CC_ETACOR(ETAAI,ZKDIA(KOFAFI),ZKDIA(KOFIFJ),
     *                  WORK(KCMOF),LUN1,LUN2,WORK(KEND1),LWRK1)
C
      ENDIF
C
C---------------------
C     Reorder results.
C---------------------
C
      KOFAFI = 2*NT1AMX + NMATIJ(1) + NMATAB(1) + 1
      CALL CC_ETARE(ETAAI,ZKDIA(KOFAFI),WORK(KEND1),LWRK1)
C
C---------------------------------
C     Write out eta-ai and eta-aI.
C---------------------------------
C
      IF ((IPRINT .GT. 20).OR.(LOCDBG)) THEN
C
         CALL AROUND('Eta-kappa-bar-0-ai vector exiting CC_DEN_PT')
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
            CALL OUTPUT(ETAAI(KOFF),1,NVIR(ISYM),1,NRHFS(ISYM),
     *                  NVIR(ISYM),NRHFS(ISYM),1,LUPRI)
C
            IF ((NVIR(ISYM) .EQ. 0) .OR. (NRHFS(ISYM) .EQ. 0)) THEN
               WRITE(LUPRI,*) 'This sub-symmetry is empty'
            ENDIF
C
  20     CONTINUE
      ENDIF
C
      IF ((IPRINT .GT. 9).OR.(LOCDBG)) THEN
         ETAKAN = DDOT(NALLAI(1),ETAAI,1,ETAAI,1)
         WRITE(LUPRI,*) 'CC_DEN_PTFC '
         WRITE(LUPRI,*) 'Norm of occupied-virtual block:', ETAKAN
      ENDIF
C
C-------------------------------------------
C     Write out frozen block of kappa-bar-0.
C-------------------------------------------
C
      IF (((IPRINT .GT. 9).OR.(LOCDBG)) .AND. (FROIMP)) THEN
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
      IF (FROIMP) THEN
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
         WRITE (LUPRI,*) 'Total time used in CC_DEN_PTFC:', TIMTOT
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
C---------------------------------------------------------------
C If energy test add nuclear nuclear repulsion energy and write out E-ccsd.
C---------------------------------------------------------------
C
      IF (ltste) THEN
C
         ECCSD = ECCSD1 + ECCSD2 + POTNUC
C
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) 'Coupled Cluster energy constructed'
         WRITE(LUPRI,*) 'from density matrices:'
         !IF (CCSD) WRITE(LUPRI,*) 'CCSD-(type) energy:', ECCSD
         WRITE(LUPRI,'(A,f15.10)') 'H1 energy, ECCSD1(type)  = ',ECCSD1
         WRITE(LUPRI,'(A,f15.10)') 'H2 energy, ECCSD2(type)  = ',ECCSD2
         WRITE(lupri,'(A,f15.10)') 'sum_p e_p kbar_pp        = ',SKAPEP
         WRITE(lupri,'(A,f15.10)') 'D(T)_ia contrib inc      = ',EN1PT
         WRITE(lupri,'(A,f15.10)') 'd(T)_pqrs contribution   = ',EN2PT
         WRITE(LUPRI,'(A,f15.10)') 'Nuc. Pot. energy         = ',POTNUC
         WRITE(lupri,'(A,f15.10)') 'CCSD(T) energy ?         = ',
     &        ECCSD1+EN1PT+ECCSD2+EN2PT+SKAPEP + POTNUC
         !WRITE(LUPRI,*) 'OBS POTNUC is missing!!! '
         
C
      ENDIF
C----------------------------------------------------------------------
      CALL QEXIT('CC_DEN_PTFC')
      RETURN
      END
C----------------------------------------------------------------------
C  /* Deck ccpt_den2fc */
      SUBROUTINE CCPT_DEN2FC(ETAAI,ZKDIA,
     &                        D1PTIA,
     &                        WORK,LWORK,
     &                        IOPT,LTSTEN,en2pt)
C
C     Written by S. Coriani, based on CC_DEN_PT
C     January 2002
C
C     Version: 1.0
C
C     Purpose: 
C     drive the calculation of the "pure d_pqrs(T)" contributions to the 
C     ^kappabar-eta_pq RHS of the orbital multipliers
C     LTSTEN = true, test densities via energy calculation
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxorb.h>
#include <maxash.h>
#include <mxcent.h>
#include <aovec.h>
#include <iratdef.h>
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (FOUR = 4.0D0)
      DIMENSION INDEXA(MXCORB)
      DIMENSION ETAAI(*), ZKDIA(*), WORK(LWORK)
      DIMENSION D1PTIA(*)
      LOGICAL LTSTEN
#include <ccorb.h>
CCN #include <infind.h> not compatible with R12!
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

      LOGICAL LETAB, LETIJ, LETAI
      LOGICAL LOCDBG
      PARAMETER (LOCDBG=.FALSE.)
C
      CALL QENTER('CCPT_DEN2FC')
C
      CALL HEADER('Construct part of rhs for CCSD(T)-kappa-0',-1)
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


      IF (LTSTEN) EN2PT=ZERO
C
C----------------------------------------------------
C     Both zeta- and t-vectors are totally symmetric.
C----------------------------------------------------
C
      ISYMTR = 1
      ISYMOP = 1
C
C----------------------------------------
C     Get CMO coefficients
C----------------------------------------
C
      KT1AM   = 1
      KD1PTAI = KT1AM + NT1AMX
      KD1PTAB = KD1PTAI + NT1AM(1)
      KD1PTIJ = KD1PTAB + NMATAB(1)
      KEND0   = KD1PTIJ + NMATIJ(1)

      call dzero(work(kd1ptai),nt1am(1))
      call dzero(work(kd1ptij),nmatij(1))
      call dzero(work(kd1ptab),nmatab(1))
      call dzero(work(kt1am),nt1amx)

      KCMO  = KEND0
      KEND1 = KCMO  + NLAMDS
      LWRK1 = LWORK - KEND1
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for allocation 1 CCPT_DEN2')
      ENDIF
C
      IF (FROIMP) THEN
         KCMOF = KEND1
         KEND1 = KCMOF + NLAMDS
         LWRK1 = LWORK - KEND1
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT('Insufficient memory for allocation 2 CCPT_DEN2')
         ENDIF
C
C-------------------------------------------
C        Get the FULL MO coefficient matrix.
C-------------------------------------------
C
         CALL CMO_ALL(WORK(KCMOF),WORK(KEND1),LWRK1)
C
      END IF
C
C----------------------------------------------------------
C     Read MO-coefficients from interface file and reorder.
C----------------------------------------------------------
C
      LUSIFC = -1
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',
     *            IDUMMY,.FALSE.)
      REWIND LUSIFC
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KCMO+I-1), I=1,NLAMDS)
      CALL GPCLOSE (LUSIFC,'KEEP')
C
      CALL CMO_REORDER(WORK(KCMO),WORK(KEND1),LWRK1)

C
C-------------------------------------
C  Two-electron part starts here.....
C-------------------------------------
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
                CALL QUIT('Insufficient core in CC_DEN_PT2')
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
C---------------------------------------------------------
C              Sonia
C              Work space allocation for the (T) densities
C              with third index backtransformed to gamma
C              All gammas together
C---------------------------------------------------------
C
               ISYDEN = ISYMD
C
               KD2IJG_PT = KEND1
               KD2AIG_PT = KD2IJG_PT + ND2IJG(ISYDEN)
               KD2IAG_PT = KD2AIG_PT + ND2AIG(ISYDEN)
               KD2ABG_PT = KD2IAG_PT + ND2AIG(ISYDEN)
               KEND2     = KD2ABG_PT + ND2ABG(ISYDEN)
               LWRK2     = LWORK  - KEND2
C
               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:',KEND2
                  CALL QUIT('Insufficient space for allocation '//
     &                      '2and1/2 in CCPT_DEN2')
               ENDIF
C
C-------------------------------------------------------
C              Initialize 4 two electron density arrays.
C-------------------------------------------------------
C
               CALL DZERO(WORK(KD2IJG_PT),ND2IJG(ISYDEN))
               CALL DZERO(WORK(KD2AIG_PT),ND2AIG(ISYDEN))
               CALL DZERO(WORK(KD2IAG_PT),ND2AIG(ISYDEN))
               CALL DZERO(WORK(KD2ABG_PT),ND2ABG(ISYDEN))
C
C-------------------------------------------------------------------
C              Calculate the two electron density d(pq,gamma;delta).
C-------------------------------------------------------------------
C
               AUTIME = SECOND()
C
               CALL CC_DEN2_PT(WORK(KD2IJG_PT),WORK(KD2AIG_PT),
     &                         WORK(KD2IAG_PT),WORK(KD2ABG_PT),
     &                         WORK(KCMO),1,
     &                         WORK(KEND2),LWRK2,
     &                         IDEL,ISYMD)
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
     &                      '3 in CC_DEN_PT2')
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
                     KD2GIJ = KD2IJG_PT + ID2IJG(ISYMPQ,ISYMG)
     *                      + NMATIJ(ISYMPQ)*(G - 1) 
                     KD2GAI = KD2AIG_PT + ID2AIG(ISYMPQ,ISYMG)
     *                      + NT1AM(ISYMPQ)*(G - 1)
                     KD2GIA = KD2IAG_PT + ID2AIG(ISYMPQ,ISYMG)
     *                      + NT1AM(ISYMPQ)*(G - 1)
                     KD2GAB = KD2ABG_PT + ID2ABG(ISYMPQ,ISYMG)
     *                      + NMATAB(ISYMPQ)*(G - 1)
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
                        CALL QUIT('Insufficient space in CC_DEN_PT2')
                     ENDIF
C
                     CALL DZERO(WORK(KD2AOB),N2BST(ISYMPQ))
c
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
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
                        CALL CCPT_FD2BL(WORK(KFD2II),WORK(KFD2IJ),
     *                                  WORK(KFD2JI),WORK(KFD2AI),
     *                                  WORK(KFD2IA),D1PTIA,
     &                                  WORK(KCMO),WORK(KCMOF),
     *                                  WORK(KEND4),LWRK4,
     &                                  IDEL,ISYMD,G,ISYMG)

                        CALL CC_FD2AO(WORK(KD2AOB),WORK(KFD2II),
     *                                WORK(KFD2IJ),WORK(KFD2JI),
     *                                WORK(KFD2AI),WORK(KFD2IA),
     *                                WORK(KCMOF),WORK(KCMO),
     *                                WORK(KCMO),WORK(KEND4),LWRK4,
     *                                ISYMPQ)

                        CALL CCPT_D2GAF(WORK(KD2GIA),
     *                                  D1PTIA,WORK(KCMOF),
     &                                  IDEL,ISYMD,G,ISYMG)
                     END IF

cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
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
C---------------------------------------------------------
C
C---------------------------------------------------------
C
                     if (LTSTEN) then

                        CALL CC_DENAO(WORK(KD2AOB),ISYMPQ,
     *                                WORK(KD2GAI),WORK(KD2GAB),
     *                                WORK(KD2GIJ),WORK(KD2GIA),ISYMPQ,
     *                                WORK(KCMO),1,WORK(KCMO),1,
     *                                WORK(KEND4),LWRK4)
C
C----------------------------------------------------------------------
C                    Calculate the 2 e- density contribution to E-ccsd.
C----------------------------------------------------------------------
C
                        EN2PT = EN2PT + HALF*DDOT(N2BST(ISYMPQ),
     *                                    WORK(KD2AOB),1,WORK(KINTAO),1)
C
                     end if
C
C-----------------------------------------------
C                    Work space allocation five.
C-----------------------------------------------
C
                        KIJINT = KEND4
                        KAIINT = KIJINT + NMATIJ(ISYMPQ)
                        KIAINT = KAIINT + NT1AM(ISYMPQ)
                        KABINT = KIAINT + NT1AM(ISYMPQ)
                        KEND5  = KABINT + NMATAB(ISYMPQ)
                        LWRK5  = LWORK  - KEND5
!Sonia: skod allocations
                        KABTIN = KEND5
                        KIJTIN = KABTIN + NMATAB(ISYMPQ)
                        KD2TAB = KIJTIN + NMATIJ(ISYMPQ)
                        KD2TIJ = KD2TAB + NMATAB(ISYMPQ)
                        KEND5  = KD2TIJ + NMATIJ(ISYMPQ)
                        LWRK5  = LWORK  - KEND5
c
                        IF (FROIMP) THEN
                           KIIFRO = KEND5
                           KIJFRO = KIIFRO + NFROFR(ISYMPQ)
                           KJIFRO = KIJFRO + NCOFRO(ISYMPQ)
                           KAIFRO = KJIFRO + NCOFRO(ISYMPQ)
                           KIAFRO = KAIFRO + NT1FRO(ISYMPQ)
                           KEND5  = KIAFRO + NT1FRO(ISYMPQ)
                           LWRK5  = LWORK  - KEND5
                        END IF
c
                        IF (LWRK5 .LT. 0) THEN
                           WRITE(LUPRI,*) 'Available:', LWORK
                           WRITE(LUPRI,*) 'Needed:', KEND5
                           CALL QUIT('Insufficient work space '//
     &                               'in CC_DEN_PT2')
                        ENDIF
C
C----------------------------------------------------------------
C                       Transform 2 integral indices to MO-basis.
C----------------------------------------------------------------
C
                        ISYM = ISYMPQ
                        CALL CCDINTMO(WORK(KIJINT),WORK(KIAINT),
     *                                WORK(KABINT),WORK(KAIINT),
     *                                WORK(KINTAO),WORK(KCMO),
     *                                WORK(KCMO),WORK(KEND5),
     *                                LWRK5,ISYM)

                        IF (FROIMP) THEN
C
C Prepare integrals g_pq (gam,del) where one index is frozen
C
                           ISYM = ISYMPQ
                           CALL CCIFROMO(WORK(KIJFRO),WORK(KJIFRO),
     *                                   WORK(KAIFRO),WORK(KIAFRO),
     *                                   WORK(KIIFRO),WORK(KINTAO),
     *                                   WORK(KCMO),WORK(KCMO),
     *                                   WORK(KCMOF),WORK(KEND5),
     *                                   LWRK5,ISYM)
C
                        ENDIF
C
C-----------------------------------------------------------------
C  A lot of unnecessary allocs!!!! 
C                Set up transposed integrals and densities.
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

                        IF (IOPT.EQ.2) THEN

                           KETAIJ = 1
                           KETAAB = 1 + NMATIJ(1)

                           CALL CCPT_ETARS_2E(
     &                                  ZKDIA(KETAIJ),ZKDIA(KETAAB),
     &                                  WORK(KIJINT),WORK(KAIINT),
     &                                  WORK(KIAINT),WORK(KABINT),
     &                                  WORK(KD2GIJ),WORK(KD2GAI),
     &                                  WORK(KD2GIA),WORK(KD2GAB),
     &                                  WORK(KEND5),LWRK5,ISYM)
                        END IF

                        CALL CCPT_ETAAI_2E(ETAAI,
     &                                WORK(KIJINT),WORK(KAIINT),
     &                                WORK(KIAINT),WORK(KABINT),
     &                                WORK(KD2GIJ),WORK(KD2GAI),
     &                                WORK(KD2GIA),WORK(KD2GAB),
     &                                WORK(KEND5),LWRK5,
     &                                ISYM)

                        IF (FROIMP) THEN

                           ISYM = ISYMPQ
!
!Contributi a eta_ai dai loop sui frozen
!
                           CALL CCFRETAI(ETAAI,WORK(KIJFRO),
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
C                          hand-side of virtual-frozen multipliers. et_aI
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

!                          CALL CCDIAZK0(ZKDIA,WORK(KIJINT),
!     *                                   WORK(KAIINT),
!     *                                   WORK(KIAINT),WORK(KABINT),
!     *                                   WORK(KD2GIJ),WORK(KD2GAI),
!     *                                   WORK(KD2GIA),WORK(KD2GAB),
!     *                                   WORK(KIJTIN),WORK(KABTIN),
!     *                                   WORK(KD2TIJ),WORK(KD2TAB),
!     *                                   WORK(KT1AM),WORK(KEND5),
!     *                                   LWRK5,ISYM)

                          CALL CCFRETAB(ZKDIA,WORK(KIJFRO),
     *                                   WORK(KJIFRO),WORK(KAIFRO),
     *                                   WORK(KIAFRO),WORK(KFD2IJ),
     *                                   WORK(KFD2JI),WORK(KFD2AI),
     *                                   WORK(KFD2IA),WORK(KT1AM),
     *                                   WORK(KEND5),LWRK5,ISYM)

                          CALL CCFRETIJ(ZKDIA,WORK(KIJFRO),
     *                                   WORK(KJIFRO),WORK(KAIFRO),
     *                                   WORK(KIAFRO),WORK(KFD2IJ),
     *                                   WORK(KFD2JI),WORK(KFD2AI),
     *                                   WORK(KFD2IA),WORK(KT1AM),
     *                                   WORK(KEND5),LWRK5,ISYM)
C
                        END IF
  140                CONTINUE
  130             CONTINUE
C
                  AUTIME = SECOND() - AUTIME
                  TIMRES = TIMRES + AUTIME

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
C------------------------
C
C------------------------
C
      if (ltsten) then
         write(lupri,*)'CC_DEN_PT2--> EN2PT: ', EN2PT
      end if
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
         WRITE (LUPRI,*) 'Total time used in CC_DEN_PT2:', TIMTOT
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
      CALL QEXIT('CCPT_DEN2FC')
      RETURN
      END
