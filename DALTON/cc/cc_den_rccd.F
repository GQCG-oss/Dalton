!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
C
C   /* Deck cc_den_rccd */
      SUBROUTINE CC_DEN_RCCD(POTNUC,ETAAI,ZKDIA,WORK,LWORK,
     &                     IOPT,IMODEL,LTSTE)
C
C     Written by S. Coriani, based on CC_DEN_PTFC
C     Debugged version using particle-symmetrized densities
C
C     Version: 3.0
C
C     Current models: RCCD & DRCCD & SOSEX
C
C     LTSTE = .true. test densities via Energy calculation
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxorb.h"
#include "maxash.h"
#include "mxcent.h"
#include "aovec.h"
#include "iratdef.h"
      PARAMETER (ZERO = 0.0D0,HALF=0.5D0,ONE = 1.0D0,TWO = 2.0D0)
      PARAMETER (TRE = 3.0D0, FOUR = 4.0D0)
      DIMENSION INDEXA(MXCORB_CC)
      DIMENSION ETAAI(*), ZKDIA(*), WORK(LWORK)
      LOGICAL LTSTE, LETAFI, LETIFJ
#include "ccorb.h"
#include "ccisao.h"
#include "r12int.h"
#include "inftap.h"
#include "blocks.h"
#include "ccfield.h"
#include "ccsdinp.h"
#include "ccinftap.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "distcl.h"
#include "cbieri.h"
#include "eritap.h"
#include "ccfro.h"
CAMT
C#include "dftcom.h"
C#include "oepopt.h"
C#include "ccandy.h"

CTEST
C#include "ccfop.h"

      CHARACTER MODEL*10
      CHARACTER NAME1*8
      CHARACTER NAME2*8

      LOGICAL LOCDBG
      PARAMETER (LOCDBG=.FALSE.)
C
      CALL QENTER('CC_DEN_RCCD')

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
        !IF (LPRNCC) THEN
         CALL HEADER('CC_DEN_RCCD: constructing RHS for RCCD-kapbar-0',
     &               -1)
         call flshfo(lupri)
        !ENDIF
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
      KZ2TIL = KLAMDH + NLAMDT   !2C-E of multipliers
      KZ2PCK = KZ2TIL + NT2SQ(1)
      KT2SQ  = KZ2PCK + NT2AMX
      KEND1  = KT2SQ  + NT2SQ(1)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for initial allocation '//
     &             'in CC_DEN_RCCD')
      ENDIF
C
C----------------------------------------
C     Read zero-th order zeta amplitudes.
C----------------------------------------
C
      IOPTRW   = 2
      CALL CC_RDRSP('L0',0,1,IOPTRW,MODEL,WORK(KEND1),WORK(KZ2AM))
      call flshfo(lupri)
C
C--------------------------------------------------------
C     Calculate tbar_tilde = 2C-E of Tbar for RCCD and
C     for dRCCD just 2*tbar in squared form
C     and save a packed copy of Tbar in KZ2PCK
C--------------------------------------------------------
C
      CALL DCOPY(NT2AMX,WORK(KZ2AM),1,WORK(KZ2PCK),1)
      CALL DCOPY(NT2AMX,WORK(KZ2AM),1,WORK(KT2AM),1)
      if (RCCD) then
        ISYOPE = 1
        IOPTTCME = 1
        CALL CCSD_TCMEPK(WORK(KT2AM),1.0D0,ISYOPE,IOPTTCME)
      else
        CALL DSCAL(NT2AMX,TWO,WORK(KT2AM),1)
      end if
      CALL CC_T2SQ(WORK(KT2AM),WORK(KZ2TIL),1)
C
C-------------------------------------------------------------
C     Square up zeta2 amplitudes.
C-------------------------------------------------------------
C
      CALL DCOPY(NT2AMX,WORK(KZ2AM),1,WORK(KT2AM),1)
      CALL CC_T2SQ(WORK(KT2AM),WORK(KZ2AM),1)
C
C-------------------------------------------
C     Read zero'th order cluster amplitudes.
C-------------------------------------------
C
      IOPTRW = 2
      CALL CC_RDRSP('R0',0,1,IOPTRW,MODEL,WORK(KEND1),WORK(KT2AM))
      CALL CC_T2SQ(WORK(KT2AM),WORK(KT2SQ),1)
C
C-------------------------------------------------
C     Set up 2C-E of cluster amplitudes (T2 tilde).
C     for RCCD and SOSEX, otherwise just 2*ampl
C-------------------------------------------------
C
      ISYOPE = 1
      CALL DCOPY(NT2AMX,WORK(KT2AM),1,WORK(KT2AMT),1)
      if (DRCCD) then
         if (SOSEX) then
            IOPTTCME = 1
            CALL CCSD_TCMEPK(WORK(KT2AMT),1.0D0,ISYOPE,IOPTTCME)
         else
           CALL DSCAL(NT2AMX,TWO,WORK(KT2AMT),1)
         end if
      else !if (RCCD) then
         IOPTTCME = 1
         CALL CCSD_TCMEPK(WORK(KT2AMT),1.0D0,ISYOPE,IOPTTCME)
      end if
C
C----------------------------------------------------------------
C     Calculate the lambda matrices.
C     Redundant, it's just CMO, but I let it go to avoid problems
C----------------------------------------------------------------
C
      KT1AM = KEND1
      KEND1 = KT1AM + NT1AMX
      LWRK1 = LWORK-KEND1
      CALL DZERO(WORK(KT1AM),NT1AMX)
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),WORK(KEND1),
     *               LWRK1)
C
C----------------------------------------------
C     Work space allocation one. CCSD-like part
C     Note that D(ai) = D(ia) = 0, and both 
C     D(ia) and h(ia) are stored transposed!
C----------------------------------------------
C
      KONEAI = KEND1
      KONEAB = KONEAI + NT1AMX
      KONEIJ = KONEAB + NMATAB(1)
      KONEIA = KONEIJ + NMATIJ(1)
      KXMAT  = KONEIA + NT1AMX
      KYMAT  = KXMAT  + NMATIJ(1)
      KONINT = KYMAT  + NMATAB(1)
C
      KINTAI = KONINT + N2BST(ISYMOP)
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
        CALL QUIT('Insuff. memory for allocation 1 CC_DEN_RCCD')
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
C     Get the non-frozen MO coefficient matrix reorder.
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
      CALL DZERO(WORK(KONEAI),NT1AMX)
C
C--------------------------------------------------------
C     Calculate X-intermediate of zeta- and t-amplitudes.
C--------------------------------------------------------
C
      CALL CC_XI(WORK(KXMAT),WORK(KZ2AM),ISYMTR,WORK(KT2AM),ISYMOP,
     *             WORK(KEND1),LWRK1)
      CALL DSCAL(NMATIJ(1),TWO,WORK(KXMAT),1)
C
C--------------------------------------------------------
C     Calculate Y-intermediate of zeta- and t-amplitudes.
C--------------------------------------------------------
C
      CALL CC_YI(WORK(KYMAT),WORK(KZ2AM),ISYMTR,WORK(KT2AM),ISYMOP,
     *           WORK(KEND1),LWRK1)
      CALL DSCAL(NMATAB(1),TWO,WORK(KYMAT),1)
C
C------------------------------------------------------------------------
C     Construct non-zero blocks of one electron density.
C     Note that X and Y are actually 2*X and 2*Y
C------------------------------------------------------------------------
C
      CALL DCOPY(NMATAB(1),WORK(KYMAT),1,WORK(KONEAB),1)
      CALL CC_EITR(WORK(KONEAB),WORK(KONEIJ),WORK(KEND1),LWRK1,1)
      CALL DIJGEN(WORK(KONEIJ),WORK(KXMAT))

C rescale X and Y back to "true" X and Y value
      CALL DSCAL(NMATIJ(1),HALF,WORK(KXMAT),1)
      CALL DSCAL(NMATAB(1),HALF,WORK(KYMAT),1)

      IF (LOCDBG) THEN
         DIJNO = DDOT(NMATIJ(1),WORK(KONEIJ),1,WORK(KONEIJ),1)
         DAINO = DDOT(NT1AMX,WORK(KONEAI),1,WORK(KONEAI),1)
         DIANO = DDOT(NT1AMX,WORK(KONEIA),1,WORK(KONEIA),1)
         DABNO = DDOT(NMATAB(1),WORK(KONEAB),1,WORK(KONEAB),1)
         WRITE(LUPRI,*) 'CC_DEN_RCCD: IOPT =  ', IOPT
         WRITE(LUPRI,*) 
     &   'Norm of ',MODEL(1:5),' one electron densities in MO-basis:'
         WRITE(LUPRI,*) DIJNO, DAINO, DIANO, DABNO
         call flshfo(lupri)

         write(lupri,*)'The IJ density of ', MODEL(1:5)
         CALL OUTPUT(WORK(KONEIJ),1,NRHF(1),1,NRHF(1),
     *                  NRHF(1),NRHF(1),1,LUPRI)
         write(lupri,*)'The AI density of ', MODEL(1:5)
         CALL OUTPUT(WORK(KONEAI),1,NVIR(1),1,NRHF(1),
     *                  NVIR(1),NRHF(1),1,LUPRI)
         write(lupri,*)'The IA density of ', MODEL(1:5)
         CALL OUTPUT(WORK(KONEIA),1,NVIR(1),1,NRHF(1),
     *                  NVIR(1),NRHF(1),1,LUPRI)
         write(lupri,*)'The AB density of ', MODEL(1:5)
         CALL OUTPUT(WORK(KONEAB),1,NVIR(1),1,NVIR(1),
     *                  NVIR(1),NVIR(1),1,LUPRI)

      ENDIF !locdbg
C
C---------------------------------
C     Read one-electron integrals.
C---------------------------------
C
      CALL CCRHS_ONEAO(WORK(KONINT),WORK(KEND1),LWRK1)

      IF (LTSTE) THEN
         !IF (LPRNCC) write(lupri,*)'LTSTE=', LTSTE
         !call flshfo(lupri)

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

      END IF !LTSTE
C
C---------------------------------------------------------
C        Ove 24-20-1996
C        Read one-electron integrals into Fock-matrix for
C        finite field.
C---------------------------------------------------------
C
      DO 13 IF = 1, NFIELD
         FF = EFIELD(IF)
         CALL CC_ONEP(WORK(KONINT),WORK(KEND1),LWRK1,FF,1,LFIELD(IF))
 13   CONTINUE
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
C
      IF (FROIMP) THEN
C
         ISYM = 1
         !obtain integrals with frozen indices
         ! h_Ij h_jI h_aJ h_Ia h_IJ
         !
         CALL CCIFROMO(WORK(KFROIJ),WORK(KFROJI),WORK(KFROAI),
     *                 WORK(KFROIA),WORK(KFROII),WORK(KONINT),
     *                    WORK(KLAMDP),WORK(KLAMDH),WORK(KCMOF),
     *                    WORK(KEND1),LWRK1,ISYM)
C
         !calculate D_II = 2 delta_IJ
         CALL CCFD1II(WORK(KFD1II))
C
      ENDIF !froimp
C
C--------------------------------------------------
C     Set up transposed integrals and densities.
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
C     Calculate one electron contribution to Zeta-kappa-0.
C------------------------------------------------------------
C
      ISYM = 1

      IF (IOPT.EQ.2) THEN

         !I let it go thru this unaltered as T1AM is set to zero
         !compute eta_ij
         KOFFIJ = 1
         !CALL DZERO(NMATIJ(1),ZKDIA(KOFFIJ))
         CALL RCCD_ETIJ(ZKDIA(KOFFIJ),WORK(KINTIJ),WORK(KINTAI),
     &                WORK(KINTIA),WORK(KINTAB),WORK(KONEIJ),
     &                WORK(KONEAI),WORK(KONEIA),WORK(KONEAB),
     &                WORK(KEND1),LWRK1,ISYM)
c        write(lupri,*)'Norm of eta_ij=',
c    & SQRT(ABS(DDOT(NMATIJ(1),ZKDIA(KOFFIJ),1,ZKDIA(KOFFIJ),1)))

         !I let it go thru this unaltered as T1AM is set to zero
         !compute eta_ab
         KOFFAB = NMATIJ(1) + 1
         !CALL DZERO(NMATAB(1),ZKDIA(KOFFAB))
         CALL RCCD_ETAB(ZKDIA(KOFFAB),WORK(KINTIJ),WORK(KINTAI),
     &                WORK(KINTIA),WORK(KINTAB),WORK(KONEIJ),
     &                WORK(KONEAI),WORK(KONEIA),WORK(KONEAB),
     &                WORK(KEND1),LWRK1,ISYM)
c        write(lupri,*)'Norm of eta_ab=',
c    & SQRT(ABS(DDOT(NMATAB(1),ZKDIA(KOFFAB),1,ZKDIA(KOFFAB),1)))

      END IF !iopt2

C------------------------------------------------------------

      CALL DZERO(ETAAI,NALLAI(1))
      CALL DZERO(WORK(KT1AM),NT1AMX)
      !I let it go thru this unaltered as T1AM is set to zero
      !eta_ai
      CALL CCDENZK0(ETAAI,WORK(KINTIJ),WORK(KINTAI),WORK(KINTIA),
     *              WORK(KINTAB),WORK(KONEIJ),WORK(KONEAI),
     *              WORK(KONEIA),WORK(KONEAB),WORK(KINIJT),
     *              WORK(KINABT),WORK(KD1IJT),WORK(KD1ABT),
     *              WORK(KT1AM),WORK(KEND1),LWRK1,ISYM)
C
!      write(lupri,*)'The eta_IJ one-electron of ', MODEL(1:5)
!      CALL OUTPUT(ZKDIA(KOFFIJ),1,NRHF(1),1,NRHF(1),
!     *                 NRHF(1),NRHF(1),1,LUPRI)
!      write(lupri,*)'The eta_AI one-electron of ', MODEL(1:5)
!      CALL OUTPUT(ETAAI,1,NVIR(1),1,NRHF(1),
!     *                  NVIR(1),NRHF(1),1,LUPRI)
!         write(lupri,*)'The eta_IA one-electron of RPA'
!         CALL OUTPUT(WORK(KINTIA),1,NVIR(1),1,NRHF(1),
!     *                  NVIR(1),NRHF(1),1,LUPRI)
!       write(lupri,*)'The eta_AB one-electron of ', MODEL(1:5)
!       CALL OUTPUT(ZKDIA(KOFFAB),1,NVIR(1),1,NVIR(1),
!     *                  NVIR(1),NVIR(1),1,LUPRI)


      IF (FROIMP) THEN
C
C--------------------------------------------------------
C       Calculate one-electron contribution to right-
C       hand-side of correlated-frozen multipliers.
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
C       Calculate one-electron contribution to right-
C       hand-side of virtual-frozen multipliers.
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

      ENDIF !froimp
C
      TIMONE = SECOND() - TIMONE
C
C--------------------------------------------
C     Start the loop over 2-e integrals.
C     Salva tutto quanto definito fino ad ora
C--------------------------------------------
C
      ! IF (LPRNCC) 
      !&  
      WRITE(LUPRI,*)'DONE WITH 1-E, START 2-E, 1e Energy=', ECCSD1
      call flshfo(lupri)

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
      ENDIF  !direct
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
                  CALL QUIT('Insufficient core in CC_DEN_RCCD')
               END IF
C
            ELSE
               NUMDIS = 1
            ENDIF !direct
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
     &                      '2 in CC_DEN_RCCD')
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
C              Calculate the RCCD  two electron density d(pq,gamma;delta).
C-----------------------------------------------------------------------------
C
               AUTIME = SECOND()
C
               CALL CC_DEN2_RCCD(WORK(KD2IJG),WORK(KD2AIG),
     &                      WORK(KD2IAG),WORK(KD2ABG),
     &                      WORK(KZ2PCK),WORK(KZ2AM),
     &                      WORK(KT2AM),WORK(KT2AMT),WORK(KT2SQ),
     &                      WORK(KZ2TIL),WORK(KXMAT),
     &                      WORK(KYMAT),WORK(KONEAB),WORK(KONEAI),
     &                      WORK(KONEIA),WORK(KLAMDH),1,
     &                      WORK(KLAMDP),1,WORK(KEND2),LWRK2,IDEL,ISYMD)
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

               ENDIF !froimp
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
                     ENDIF !froimp
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
C                    Calculate the 2 e- density contribution to E_rccd
C----------------------------------------------------------------------
C
                        ECCSD2 = ECCSD2 + HALF*DDOT(N2BST(ISYMPQ),
     *                                    WORK(KD2AOB),1,WORK(KINTAO),1)
C
                     END IF !ltste
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
     &                               'in CC_DEN_RCCD')
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
                        ENDIF !froimp
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
                        END IF !iopt2

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
                        ENDIF   !froimp
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
      if (locdbg) then
          KOFFIJ = 1
          KOFFAB = 1 + NMATIJ(1)
          KOFAFI = 2*NT1AMX + NMATIJ(1) + NMATAB(1) + 1
          KOFIFJ = 2*NT1AMX + NMATIJ(1) + NMATAB(1) + 2*NT1FRO(1) + 1
          write(lupri,*) '                         '
          write(lupri,*) 'Before call to CCSD_ZKBLO'
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
          call flshfo(lupri)
      end if !locdbg
C
C------------------------------------------------------------------------
C Calculate the ZK0(ij) and ZK0(ab) blocks (to be used to correct eta_ai)
C from eta_ij/e_i-e_j and eta_ab/e_a-e_b
C------------------------------------------------------------------------
C
      IF (IOPT.EQ.2) THEN

         CALL CCSD_ZKBLO(ZKDIA,WORK(KEND1),LWRK1)
         write(lupri,*) 'I came out of ccsd_zkblo'

         if (locdbg) then
             KOFFIJ = 1
             KOFFAB = 1 + NMATIJ(1)
             KOFAFI = 2*NT1AMX + NMATIJ(1) + NMATAB(1) + 1
             KOFIFJ = 2*NT1AMX + NMATIJ(1) + NMATAB(1) + 2*NT1FRO(1) + 1
             write(lupri,*) 'after call to CCSD_ZKBLO'
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
             call flshfo(lupri)
          end if
      END IF !locdbg
C
C------------------------------------------------------------------------
C Add the diagonal ZK0(ii) and ZK0(aa) elements in the proper places
C------------------------------------------------------------------------
C
!         if (.false.) then
!
!           IF (LTSTE) THEN
!
!             !multiply by epsilon_p and sum over p to get the energy
!
!             KFOCKDIA = KEND1
!             KEND1    = KFOCKDIA + NORBTS
!             LWRK1    = LWORK-KEND1
C
C-------------------------------------
C     Read canonical orbital energies.
C-------------------------------------
C
!             CALL DZERO(WORK(KFOCKDIA),NORBTS)
!             CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
!     &            .FALSE.)
!             REWIND (LUSIFC)
C
!             CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
!             READ (LUSIFC)
!             READ (LUSIFC) (WORK(KFOCKDIA + I - 1), I = 1,NORBTS)
C
!             CALL GPCLOSE(LUSIFC,'KEEP')
C
C----------------------------------------------------------------
C     Change symmetry ordering of the canonical orbital energies.
C----------------------------------------------------------------
C
!             IF (FROIMP)
!     *           CALL CCSD_DELFRO(WORK(KFOCKDIA),WORK(KEND1),LWRK1)
C
!                 CALL FOCK_REORDER(WORK(KFOCKDIA),WORK(KEND1),LWRK1)
C
C--------------------------------------------
C     Calculate sum_p kappabar_pp * epsilon_p
C     Occupied block:
C--------------------------------------------
C
!            SKAPEP = DDOT(NORBT,WORK(KKAPII),1,WORK(KFOCKDIA),1)
!             END IF
!           END IF
c
!      end if
!      END IF
C
C------------------------------------------------
C Correct Eta_ai with ZK0(ij) and ZK0(ab) blocks 
C------------------------------------------------
C-tbp:
c     CALL AROUND('Eta-kappa-bar-0-ai vector before modification')
c     do ISYM = 1,NSYM
c           WRITE(LUPRI,*) ' '
c           WRITE(LUPRI,444) 'Sub-symmetry block number:', ISYM
c           WRITE(LUPRI,555) '--------------------------'
c           KOFF = IALLAI(ISYM,ISYM) + 1
c           write(lupri,'(A,1P,D25.16)') 'Norm=',
c    *      sqrt(ddot(NVIR(ISYM)*NRHFS(ISYM),ETAAI(KOFF),1,
c    *                                       ETAAI(KOFF),1))
c           CALL OUTPUT(ETAAI(KOFF),1,NVIR(ISYM),1,NRHFS(ISYM),
c    *                  NVIR(ISYM),NRHFS(ISYM),1,LUPRI)
c           IF ((NVIR(ISYM) .EQ. 0) .OR. (NRHFS(ISYM) .EQ. 0)) THEN
c              WRITE(LUPRI,*) 'This sub-symmetry is empty'
c           ENDIF
c     end do

C
      IF (IOPT.EQ.2) THEN
          !SONIA: NB!! I removed the froimp stuff 
           !IF (LPRNCC) write(lupri,*)'CC_DEN_RCCD using CCETACOR'
          call flshfo(lupri)
          CALL CCETACOR(ETAAI,ZKDIA,WORK(KEND1),LWRK1)
           !IF (LPRNCC) write(lupri,*)'CC_DEN_RCCD out of etacor'
          call flshfo(lupri)
 
          KOFFIJ = 1
          KOFFAB = NMATIJ(1) + 1
          KOFAFI = 2*NT1AMX + NMATIJ(1) + NMATAB(1) + 1
          KOFIFJ = 2*NT1AMX + NMATIJ(1) + NMATAB(1) + 2*NT1FRO(1) + 1

      END IF

      XZKDIA = DDOT(LUNGO,ZKDIA,1,ZKDIA,1)
!      IF (LPRNCC) 
!     &   WRITE(LUPRI,*) 'CC_DEN_RCCD: ZKDIA before CCETACOR', XZKDIA
!           call flshfo(lupri)
C
C---------------------
C     Reorder results.
C     it does nothing if it's NOT froimp
C---------------------
C
!      KOFAFI = 2*NT1AMX + NMATIJ(1) + NMATAB(1) + 1
!      CALL CC_ETARE(ETAAI,ZKDIA(KOFAFI),WORK(KEND1),LWRK1)
C
C---------------------------------
C     Write out eta-ai and eta-aI.
C---------------------------------
C
      IF ((IPRINT .GT. 20).OR.(LOCDBG)) THEN
C
         CALL AROUND('Eta-kappa-bar-0-ai vector exiting CC_DEN_RCCD')
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
            write(lupri,'(A,1P,D25.16)') 'Norm=',
     *      sqrt(ddot(NVIR(ISYM)*NRHFS(ISYM),ETAAI(KOFF),1,
     *                                       ETAAI(KOFF),1))
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
         WRITE(LUPRI,*) 'CC_DEN_RCCD '
         WRITE(LUPRI,*) 'Norm of occupied-virtual block:', ETAKAN
         call flshfo(lupri)
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
         WRITE (LUPRI,*) 'Total time used in CC_DEN_RCCD:', TIMTOT
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
        !IF (LPRNCC) THEN
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) 'RPA energy constructed'
         WRITE(LUPRI,*) 'from density matrices:'
         !IF (CCSD) WRITE(LUPRI,*) 'CCSD-(type) energy:', ECCSD
         WRITE(LUPRI,'(A,f15.10)') 'H1 energy, ECCSD1(type)  = ',ECCSD1
         WRITE(LUPRI,'(A,f15.10)') 'H2 energy, ECCSD2(type)  = ',ECCSD2
         !WRITE(lupri,'(A,f15.10)') 'sum_p e_p kbar_pp        = ',SKAPEP
         WRITE(LUPRI,'(A,f15.10)') 'Nuc. Pot. energy         = ',POTNUC
         WRITE(lupri,'(A,f15.10)') 'CCSD energy ?         = ',
     &        ECCSD1+ECCSD2+ POTNUC
           call flshfo(lupri)
         WRITE(lupri,'(A,f15.10)') 'CCSD energy (2) ?     = ',ECCSD
         !WRITE(LUPRI,*) 'OBS POTNUC is missing!!! '
        !ENDIF
C
      ENDIF
C----------------------------------------------------------------------
      CALL QEXIT('CC_DEN_RCCD')
      RETURN
      END
C----------------------------------------------------------------------


C  /* Deck rccd_etij */
      SUBROUTINE RCCD_ETIJ(ETAIJ,XINTIJ,XINTAI,XINTIA,XINTAB,DIJ,DAI,
     *                    DIA,DAB,WORK,LWORK,ISYM)
C
C     Written by S. Coriani 2010
C
C     Version: 1.0
C
C     Purpose: To set up the right hand side of the equation for
C              zeta-kappa-0_ij (ETAIJ) from MO-integrals (XIN*) and RCCD
C              densities (D*) 
C              Note that due to the symmetry in the formulas, this
C              routine is able to handle both the one- and the two-
C              electron contributions!
C              ISYM is the symmetry of both the density and the
C              integrals!
C
#include "implicit.h"
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION ETAIJ(*), XINTIJ(*), XINTAI(*), XINTIA(*), XINTAB(*)
      DIMENSION DIJ(*), DAI(*), DIA(*), DAB(*), WORK(LWORK)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "cclr.h"
C
      CALL QENTER('RCCD_ETIJ')
C
      DO 100 ISYMI = 1,NSYM
C
C----------------------------------------------------------------
C        Calculate direct terms to eta_ij.
C----------------------------------------------------------------
C
         ISYMJ  = ISYMI
         ISYMK  = MULD2H(ISYMI,ISYM)
         ISYMC  = MULD2H(ISYMI,ISYM)
C
         KOFFRE = IMATIJ(ISYMI,ISYMJ) + 1
C
         NTOTRE = MAX(NRHF(ISYMI),1)
         NTOTI  = MAX(NRHF(ISYMI),1)
         NTOTJ  = MAX(NRHF(ISYMJ),1)
         NTOTK  = MAX(NRHF(ISYMK),1)
         NTOTC  = MAX(NVIR(ISYMC),1)
C
         KOFF1  = IMATIJ(ISYMK,ISYMI) + 1
         KOFF2  = IMATIJ(ISYMK,ISYMJ) + 1
         KOFF3  = IMATIJ(ISYMI,ISYMK) + 1
         KOFF4  = IMATIJ(ISYMJ,ISYMK) + 1
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NRHF(ISYMK),ONE,
     *              XINTIJ(KOFF1),NTOTK,DIJ(KOFF2),NTOTK,ONE,
     *              ETAIJ(KOFFRE),NTOTRE)
C
         CALL DGEMM('N','T',NRHF(ISYMI),NRHF(ISYMJ),NRHF(ISYMK),ONE,
     *              XINTIJ(KOFF3),NTOTI,DIJ(KOFF4),NTOTJ,ONE,
     *              ETAIJ(KOFFRE),NTOTRE)
C
         CALL DGEMM('N','T',NRHF(ISYMI),NRHF(ISYMJ),NRHF(ISYMK),-ONE,
     *              DIJ(KOFF3),NTOTI,XINTIJ(KOFF4),NTOTJ,ONE,
     *              ETAIJ(KOFFRE),NTOTRE)
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NRHF(ISYMK),-ONE,
     *              DIJ(KOFF1),NTOTK,XINTIJ(KOFF2),NTOTK,ONE,
     *              ETAIJ(KOFFRE),NTOTRE)
C
         KOFF5  = IT1AM(ISYMC,ISYMI) + 1
         KOFF6  = IT1AM(ISYMC,ISYMJ) + 1
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMC),ONE,
     *              XINTAI(KOFF5),NTOTC,DAI(KOFF6),NTOTC,ONE,
     *              ETAIJ(KOFFRE),NTOTRE)
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMC),ONE,
     *              XINTIA(KOFF5),NTOTC,DIA(KOFF6),NTOTC,ONE,
     *              ETAIJ(KOFFRE),NTOTRE)
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMC),-ONE,
     *              DIA(KOFF5),NTOTC,XINTIA(KOFF6),NTOTC,ONE,
     *              ETAIJ(KOFFRE),NTOTRE)
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMC),-ONE,
     *              DAI(KOFF5),NTOTC,XINTAI(KOFF6),NTOTC,ONE,
     *              ETAIJ(KOFFRE),NTOTRE)
C
  100 CONTINUE
C
      CALL QEXIT('RCCD_ETIJ')
C
      RETURN
      END
C  /* Deck rccd_etab */
      SUBROUTINE RCCD_ETAB(ETAAB,XINTIJ,XINTAI,XINTIA,XINTAB,DIJ,DAI,
     *                    DIA,DAB,WORK,LWORK,ISYM)
C
C     Written by S. Coriani 2010
C
C     Version: 1.0
C
C     Purpose: To set up the right hand side of the equation for
C              zeta-kappa-0_ab (ETAAB) from MO-integrals (XIN*) & RCCD
C              densities (D*) 
C              Note that due to the symmetry in the formulas, this
C              routine is able to handle both the one- and the two-
C              electron contributions!
C              ISYM is the symmetry of both the density and the
C              integrals!
C
#include "implicit.h"
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION ETAAB(*), XINTIJ(*), XINTAI(*), XINTIA(*), XINTAB(*)
      DIMENSION DIJ(*), DAI(*), DIA(*), DAB(*), WORK(LWORK)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "cclr.h"
C
      CALL QENTER('RCCD_ETAB')
C
      DO 100 ISYMA = 1,NSYM
C
C----------------------------------------------------------------
C        Calculate direct terms to eta_ab.
C----------------------------------------------------------------
C
         ISYMB  = ISYMA
         ISYMK  = MULD2H(ISYMA,ISYM)
         ISYMC  = MULD2H(ISYMA,ISYM)
C
         KOFFRE = IMATAB(ISYMA,ISYMB) + 1
C
         NTOTRE = MAX(NVIR(ISYMA),1)
         NTOTA  = MAX(NVIR(ISYMA),1)
         NTOTB  = MAX(NVIR(ISYMB),1)
         NTOTC  = MAX(NVIR(ISYMC),1)
         NTOTK  = MAX(NRHF(ISYMK),1)
C
         KOFF1  = IT1AM(ISYMA,ISYMK) + 1
         KOFF2  = IT1AM(ISYMB,ISYMK) + 1
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMK),ONE,
     *              XINTIA(KOFF1),NTOTA,DIA(KOFF2),NTOTB,ONE,
     *              ETAAB(KOFFRE),NTOTRE)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMK),-ONE,
     *              DAI(KOFF1),NTOTA,XINTAI(KOFF2),NTOTB,ONE,
     *              ETAAB(KOFFRE),NTOTRE)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMK),-ONE,
     *              DIA(KOFF1),NTOTA,XINTIA(KOFF2),NTOTB,ONE,
     *              ETAAB(KOFFRE),NTOTRE)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMK),ONE,
     *              XINTAI(KOFF1),NTOTA,DAI(KOFF2),NTOTB,ONE,
     *              ETAAB(KOFFRE),NTOTRE)
C
         KOFF3  = IMATAB(ISYMC,ISYMA) + 1
         KOFF4  = IMATAB(ISYMC,ISYMB) + 1
         KOFF5  = IMATAB(ISYMA,ISYMC) + 1
         KOFF6  = IMATAB(ISYMB,ISYMC) + 1
C
         CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),NVIR(ISYMC),ONE,
     *              XINTAB(KOFF3),NTOTC,DAB(KOFF4),NTOTC,ONE,
     *              ETAAB(KOFFRE),NTOTRE)
C 
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NVIR(ISYMC),-ONE,
     *              DAB(KOFF5),NTOTA,XINTAB(KOFF6),NTOTB,ONE,
     *              ETAAB(KOFFRE),NTOTRE)
C
         CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),NVIR(ISYMC),-ONE,
     *              DAB(KOFF3),NTOTC,XINTAB(KOFF4),NTOTC,ONE,
     *              ETAAB(KOFFRE),NTOTRE)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NVIR(ISYMC),ONE,
     *              XINTAB(KOFF5),NTOTA,DAB(KOFF6),NTOTB,ONE,
     *              ETAAB(KOFFRE),NTOTRE)
C
  100 CONTINUE
C
C
      CALL QEXIT('RCCD_ETAB')
C
      RETURN
      END
C  /* Deck Rccdenzk0 */
      SUBROUTINE RCCDENZK0(ETAKA,XINTIJ,XINTAI,XINTIA,XINTAB,DIJ,DAI,
     *                    DIA,DAB,XIJT,XABT,DIJT,DABT,
     *                    WORK,LWORK,ISYM)
C
C     Written by Sonia Halkier 28/1 - 2011 - based on CCDENZK0
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
#include "implicit.h"
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION ETAKA(*), XINTIJ(*), XINTAI(*), XINTIA(*), XINTAB(*)
      DIMENSION DIJ(*), DAI(*), DIA(*), DAB(*), XIJT(*), XABT(*)
      DIMENSION DIJT(*), DABT(*), WORK(LWORK)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
      CALL QENTER('RCCDENZK0')
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
      CALL QEXIT('RCCDENZK0')
C
      RETURN
      END

C  /* Deck rccd_zkblo */
      SUBROUTINE RCCD_ZKBLO(ZKDIA,WORK,LWORK)
C
C     Written by Sonia Halkier  2011
C
C     Version: 1.0
C
C     Purpose: To calculate the ab & ij parts of the CCSD kappa-bar-0,
C              from right-hand-sides (ZKDIA on input) and canonical
C              orbital energies.
C     Control numerical instabilities. Sonia, 2002
C
#include "implicit.h"
#include "dummy.h"
      PARAMETER(ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER(FOUR = 4.0D0, EIGHT=8.0D0)
      PARAMETER (THRDNM = 1.0D-12)
      DIMENSION ZKDIA(*), WORK(LWORK)
#include "priunit.h"
#include "maxorb.h"
#include "ccorb.h"
#include "iratdef.h"
#include "inftap.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
C
      CALL QENTER('RCCD_ZKBLO')
C-tbp:
      write(lupri,*)
      write(lupri,*) 'RCCD_ZKBLO WARNING: I''m not numerically stable!'
      write(lupri,*) 'RCCD_ZKBLO WARNING: use CCSD_ZKBLO instead...'
      write(lupri,*)
      call flshfo(lupri)
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
!      do i=1,NORBTS
!         write(lupri,*) 'Epsilon_',i,'=',WORK(KFOCKD+i-1)
!      end do
C
C---------------------------
C     Calculate the results:
C     Occupied block:
C---------------------------
      DO 100 ISYMI = 1,NSYM
         ISYMJ = ISYMI
         DO 110 J = 1,NRHF(ISYMJ)
            KOFFJ = KFOCKD + IRHF(ISYMJ) + J - 1
            DO 120 I = J+1,NRHF(ISYMI)
               KOFFI = KFOCKD + IRHF(ISYMI) + I - 1
               DENOM = WORK(KOFFJ) - WORK(KOFFI)
               !write(lupri,*) 'Denom IJ in RCCD_ZKBLO is:', DENOM
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
               !write(lupri,*) 'Denom AB in RCCD_ZKBLO is:', DENOM
               IF (ABS(DENOM) .GT. THRDNM) THEN
                  NAB = IMATAB(ISYMA,ISYMB) + NVIR(ISYMA)*(B - 1) + A
                  NBA = IMATAB(ISYMB,ISYMA) + NVIR(ISYMB)*(A - 1) + B
C
                  ZKDIA(NMATIJ(1) + NAB)= HALF*WORK(KETAAB+NAB-1)/
     *                                     DENOM
                  ZKDIA(NMATIJ(1) + NBA)= ZKDIA(NMATIJ(1) + NAB)
               END IF
C
  150       CONTINUE
  140    CONTINUE
  130 CONTINUE
C
      CALL QEXIT('RCCD_ZKBLO')
C
      RETURN
      END

!============
C  /* Deck cc_den2_rccd */
      SUBROUTINE CC_DEN2_RCCD(D2IJG,D2AIG,D2IAG,D2ABG,
     &                   Z2PK,Z2AM,T2AM,T2TILDE,T2SQ,
     &                   Z2TILDE,
     *                   XMAT,YMAT,DAB,DAI,DIA,
     *                   XLAMDH,ISYMLH,XLAMDP,ISYMLP,
     *                   WORK,LWORK,IDEL,ISYMD)
C
C     Written by  Sonia & Francesca, 18/3/2010, based on CC_DEN2
C     Debugged by Sonia & Thomas BP, 18/6/2012
C     Purpose: Directs the calculation of the 2 electron RCCD density
C              d(pq,gam;del) for a given del (IDEL). The 4 blocks pq
C              of the result is returned in D2IJG through D2ABG!
C     Densities must be particle-symmetrized!!!!!!!!!!!!!!!!!!!!!!!!!
C
#include "implicit.h"
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      INTEGER   ISYMLP, ISYMLH
      DIMENSION D2IJG(*), D2AIG(*), D2IAG(*), D2ABG(*), Z2AM(*)
      DIMENSION T2AM(*), T2TILDE(*), Z2TILDE(*), XMAT(*), YMAT(*)
      DIMENSION DAB(*), DAI(*), DIA(*), Z2PK(*), XLAMDH(*)
      DIMENSION T2SQ(*)
      DIMENSION XLAMDP(*), WORK(LWORK)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdinp.h"
C
      CALL QENTER('CC_DEN2_RCCD')
      IF (DEBUG) WRITE(LUPRI,*) 'Entered CC_DEN2_RCCD'
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
      KZ2DEL = KT2DEL + NT2BCD(ISYMTI)
      KEND1  = KZ2DEL + NT2BCD(ISYMTI)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insuff. space for 1st allocation, CC_DEN2_RCCD')
      ENDIF

C
C------------------------------------------------------
C     Calculate the delta backtransformed amplitude
C     t_ai,j;delta and multiplier tbar_ai,j;delta
C------------------------------------------------------
C
      CALL CC_TI(WORK(KT2DEL),ISYMTI,T2AM,ISYMOP,XLAMDH,ISYMLH,
     *           WORK(KEND1),LWRK1,IDEL,ISYMD)

      CALL CC_TI(WORK(KZ2DEL),ISYMTI,Z2PK,ISYMOP,XLAMDH,ISYMLH,
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
C     Contributions from projection against <u2|.
C------------------------------------------------
C
      !Y part, note that D_ab=2*Y_ba 
      FACTOR=2.0d0
      CALL RPADIJA21(WORK(KD2IJA),ISYD2H,DAB,ISYD1,XLAMDH,ISYMLH,
     *               WORK(KEND2),LWRK2,IDEL,ISYMD,FACTOR)
      IF (RCCD) THEN
         FACTOR=2.0d0
         CALL RPADIJA22(WORK(KD2IJA),ISYD2H,T2SQ,WORK(KZ2DEL),ISYMTI,
     *               WORK(KEND2),LWRK2,FACTOR)
      END IF
C
C-----------------------------------------------------------------
C     Backtransform third virtual index to AO and store in result.
C-----------------------------------------------------------------
C
      ICON = 4
      CALL CCD2_PQAO(D2IJG,ISYDEN,WORK(KD2IJA),ISYD2H,
     &               XLAMDP,ISYMLP,ICON)
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
      KZBIN  = KZINT  + NT2BCD(ISYMZI)
      KEND2  = KZBIN  + NT2BCD(ISYMZI)
      LWRK2  = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
         CALL QUIT('Insuff space for 4th allocation in CC_DEN2_RCCD')
      ENDIF
C
      CALL DZERO(WORK(KD2IAJ),NT2BCD(ISYD2H))
C
C--------------------------------------------------------------
C     Calculate ZINT of t^delta amplitudes.
C--------------------------------------------------------------
C
      ICON = 1
      CALL CC_ZWVI(WORK(KZINT),Z2AM,ISYMTR,WORK(KT2DEL),ISYMTI,
     *                  WORK(KEND2),LWRK2,ICON)
C
C--------------------------------------------------------------------
C     Calculate terms of (occ.vir,occ;del) block with Zbar (in ZINT).
C--------------------------------------------------------------------
C
      CALL DSCAL(NT2BCD(ISYMZI),4.0d0,WORK(KZINT),1)  !to get -4*t*Z
      IF (RCCD) THEN
         CALL DIAJ26(WORK(KD2IAJ),ISYD2H,T2AM,WORK(KZINT),ISYMZI,
     *               WORK(KEND2),LWRK2)
      END IF

      CALL DSCAL(NT2BCD(ISYMZI),2.0d0,WORK(KZINT),1)  !to get 8*t*Z
      CALL DIAJ28(WORK(KD2IAJ),ISYD2H,T2AM,WORK(KZINT),ISYMZI,
     *               WORK(KEND2),LWRK2)
C
C-----------------------------------------------------------------
C     Add the backtransformed 2C-E of t ampl to (occ.vir,occ;del) 
C     (= Contribution from projection against <HF|)
C-----------------------------------------------------------------
C
      !Compute backtranformed 2C-E of amplitudes and add to density
      CALL CC_TI(WORK(KT2DEL),ISYMTI,T2TILDE,ISYMOP,XLAMDH,ISYMLH,
     *           WORK(KEND2),LWRK2,IDEL,ISYMD)
      IF (ISYMTI.NE.ISYD2H) CALL QUIT('Symmetry mismatch in CC_DEN2')
      !scale by 2 
      CALL DAXPY(NT2BCD(ISYMTI),TWO,WORK(KT2DEL),1,WORK(KD2IAJ),1)
C
C--------------------------------------------------------------
C     Calculate ZINT of tbar^delta multipliers (ZBIN)
C     and place it in of the (vir.occ,occ;del) block.
C     (first contribution from proj against <u2|?)
C--------------------------------------------------------------
C
      ICON = 1
!      recycle later on ZINT
!      CALL CC_ZWVI(WORK(KZINT),T2SQ,ISYMTR,WORK(KZ2DEL),ISYMTI,
!     *                  WORK(KEND2),LWRK2,ICON)
      call dzero(WORK(KZBIN),NT2BCD(ISYMZI))
      CALL CC_ZWVI(WORK(KZBIN),T2SQ,ISYMTR,WORK(KZ2DEL),ISYMTI,
     *                  WORK(KEND2),LWRK2,ICON)
C
      !KD2AIJ = KZINT
      KD2AIJ = KZBIN
      KEND2  = KD2AIJ + NT2BCD(ISYD2H)
      LWRK2  = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
         CALL QUIT(
     *        'Insufficient space for fourth allocation in CC_DEN2')
      ENDIF
      !4*Zbar term
      CALL DSCAL(NT2BCD(ISYMZI),4.0d0,WORK(KD2AIJ),1)
C
C------------------------------------------------
C     Second contribution from projection against <u2|.
C------------------------------------------------
C
      FACTOR=1.0d0 
      CALL RPADAIJ22(WORK(KD2AIJ),ISYD2H,DAB,ISYD1,XLAMDH,ISYMLH,
     *               WORK(KEND2),LWRK2,IDEL,ISYMD,FACTOR)
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
         CALL QUIT('Insuff space for 5th alloc in CC_DEN2_RCCD')
      ENDIF
C
      CALL DZERO(WORK(KD2IAB),NCKATR(ISYD2H))
      CALL DZERO(WORK(KD2ABI),NCKASR(ISYD2H))
      CALL DZERO(WORK(KD2AIB),NCKATR(ISYD2H))
C
C--------------------------------------------------------------
C     Backtransf 2C-E of tbar: tbartilde(ai,b;del).
C     Note that this equals the d(ai,b;del) density block.
C--------------------------------------------------------------
C
      CALL DAIB21(WORK(KD2AIB),ISYD2H,Z2TILDE,XLAMDH,ISYMLH,
     &            WORK(KEND2),LWRK2,IDEL,ISYMD)
      CALL DSCAL(NCKATR(ISYD2H),1.0d0,WORK(KD2AIB),1)
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
      KZ2AO   = KD2AIB  !ricicla lo spazio
      ISYZ2AO = ISYD2H
      !j-backtr of tam placed on Z2AO
      CALL DZERO(WORK(KZ2AO),NCKATR(ISYZ2AO))
      CALL DAIB21(WORK(KZ2AO),ISYD2H,T2SQ,XLAMDH,ISYMLH,
     *            WORK(KEND2),LWRK2,IDEL,ISYMD)
!
! THIS IS THE CORRECT ONE FOR IABJ
!
      CALL DIAC22(WORK(KD2IAB),ISYD2H,Z2PK,WORK(KZ2AO),ISYZ2AO,
     &            WORK(KEND2),LWRK2)
      CALL DSCAL(NCKATR(ISYD2H),4.0d0,WORK(KD2IAB),1)
!
      IF (RCCD) THEN
        CALL DABI22(WORK(KD2ABI),ISYD2H,Z2PK,WORK(KZ2AO),ISYZ2AO,
     *               WORK(KEND2),LWRK2)
        CALL DSCAL(NCKASR(ISYD2H),2.0d0,WORK(KD2ABI),1)
      END IF
C
C--------------------------------------------------------------------
C     Calculate remaining contributions from projection against <u2|.
C--------------------------------------------------------------------
C
      FACTOR=2.0d0
      CALL RPADABI21(WORK(KD2ABI),ISYD2H,DAB,ISYD1,XLAMDH,ISYMLH,
     &               IDEL,ISYMD,FACTOR)
      
      FACTOR=2.0d0
      CALL RPADIAC21(WORK(KD2IAB),ISYD2H,YMAT,ISYD1,XLAMDH,ISYMLH,
     &               IDEL,ISYMD,FACTOR)
C
C---------------------------------------------------
C     Backtransform third index of the two remaining
C     blocks to AO and store in results.
C---------------------------------------------------

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
C     Contributions from projection against <u2|.
C     Note that X is ''true'' X, not scaled by two!!!!
C     Density must be defined symmetrized, so IJKL=KLIJ
C------------------------------------------------
C
!    Use CCSD code directly

      CALL DIJK21(WORK(KD2IJK),ISYD2H,XMAT,XLAMDH,ISYMLH,
     &            WORK(KEND2),LWRK2,IDEL,ISYMD)
      CALL DIJK23(WORK(KD2IJK),ISYD2H,XMAT,XLAMDH,ISYMLH,IDEL,ISYMD)
      CALL DIJK24(WORK(KD2IJK),ISYD2H,XMAT,XLAMDH,ISYMLH,IDEL,ISYMD)
      !scale by 2 to match RPA definition
      CALL DSCAL(NMAIJK(ISYD2H),2.0d0,WORK(KD2IJK),1)  
C
C------------------------------------------------
C     Contributions from projection against <HF|.
C------------------------------------------------
C    
      CALL DIJK01(WORK(KD2IJK),ISYD2H,XLAMDH,ISYMLH,IDEL,ISYMD)
      CALL DIJK02(WORK(KD2IJK),ISYD2H,XLAMDH,ISYMLH,IDEL,ISYMD)
C
C------------------------------------------------------------------
C     Backtransform third occupied index to AO and store in result.
C------------------------------------------------------------------
C
      ICON = 1
      CALL CCD2_PQAO(D2IJG,ISYDEN,WORK(KD2IJK),ISYD2H,
     *               XLAMDP,ISYMLP,ICON)
C
      CALL QEXIT('CC_DEN2_RCCD')
C
      RETURN
      END


C  /* Deck rpadija21 */
      SUBROUTINE RPADIJA21(D2IJA,ISYDEN,DAB,ISYDAB,XLAMDH,ISYMLH,
     &                  WORK,LWORK,IDEL,ISYMD,FACTOR)
C
C     MEMO: D_ab is 2Y already!!!!
C     Modified version of DIJA21
C     Purpose: To calculate the first contribution to D2IJA
C              from projection against doubles in RPA/RCCD/DRCCD!
C     Sonia
C
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
#if defined (SYS_CRAY)
      REAL ZERO, ONE, TWO, FACTOR
#else
      DOUBLE PRECISION ZERO, ONE, TWO, FACTOR
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
      CALL QENTER('RPADIJA21')
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
C     Calculate contraction of transposed Y-matrix and lambda matrix.
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
      CALL DGEMV('N',NVIR(ISYMA),NVIR(ISYMB),FACTOR,DAB(KOFF1),NTOTA,
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
            D2IJA(NIJA) = D2IJA(NIJA)   + ONE*WORK(A)
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('RPADIJA21')
C
      RETURN
      END

C  /* Deck rpadija22 */
      SUBROUTINE rpaDIJA22(D2IJA,ISYDEN,T2SQ,Z2INT,ISYMTI,WORK,
     &                     LWORK,FACTOR)
C
C     Purpose: To calculate the second contribution to D2IJA
C              from projection against doubles!
C              in RPA, -2 sum_ck t_aj,ck tbar_ck,i;delta 
C
#include "implicit.h"
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, FOUR = 4.0d0)
      DIMENSION D2IJA(*), T2SQ(*), Z2INT(*), WORK(LWORK)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
      CALL QENTER('RPADIJA22')
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
     *              ONE,T2SQ(KOFF1),NTOTAJ,Z2INT(KOFF2),NTOTCK,ZERO,
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
                     D2IJA(NIJA) = D2IJA(NIJA) 
     &                             -FACTOR*WORK(KSCR + NAJI - 1)
C
  140             CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('RPADIJA22')
C
      RETURN
      END

C  /* Deck rpadijk21 */
      SUBROUTINE rpaDIJK21(D2IJK,ISYIJK,XMAT,XLAMDH,ISYMLH,
     *                  WORK,LWORK,IDEL,ISYMD,FACTOR)
C
C     Purpose: To calculate the first contribution to D2IJK
C              in RCCD/DRCCD, (projection against doubles)!
C
#include "implicit.h"
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION D2IJK(*), XMAT(*), XLAMDH(*), WORK(LWORK)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
      CALL QENTER('rpaDIJK21')
C
      ISYML  = MULD2H(ISYMD,ISYMLH)
      ISYMX1 = ISYML
C
      IF (ISYML.NE.ISYIJK) CALL QUIT('Symmetry mismatch in DIJK21')
C
C-----------------------------------------------------------------
C     Calculate contraction of X-intermediate & lambda hole matrix.
C-----------------------------------------------------------------
C
      IF (LWORK .LT. NRHF(ISYMX1)) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', NRHF(ISYMX1)
         CALL QUIT('Insufficient work space in DIJK21')
      ENDIF
C
      CALL DZERO(WORK,NRHF(ISYMX1))
C
      KOFF1 = IMATIJ(ISYMX1,ISYML) + 1
      KOFF2 = IGLMRH(ISYMD,ISYML) + IDEL - IBAS(ISYMD)
C
      NTOTI = MAX(NRHF(ISYMX1),1)
C
      CALL DGEMV('N',NRHF(ISYMX1),NRHF(ISYML),FACTOR,XMAT(KOFF1),NTOTI,
     *           XLAMDH(KOFF2),NBAS(ISYMD),ZERO,WORK,1)
C
C----------------------------------
C     Calculate first contribution.
C----------------------------------
C
      ISYMI = ISYMX1
C
      DO 100 ISYMK = 1,NSYM
C
         ISYMJ  = ISYMK
         ISYMIJ = MULD2H(ISYMJ,ISYMI)
C
         DO 110 K = 1,NRHF(ISYMK)
C
            J    = K
            NIJK = IMAIJK(ISYMIJ,ISYMK) + NMATIJ(ISYMIJ)*(K - 1)
     *           + IMATIJ(ISYMI,ISYMJ)  + NRHF(ISYMI)*(J - 1) + 1
C
            CALL DAXPY(NRHF(ISYMI),ONE,WORK,1,D2IJK(NIJK),1)
C
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('rpaDIJK21')
C
      RETURN
      END

C  /* Deck rpadaij22 */
      SUBROUTINE RPADAIJ22(D2AIJ,ISYAIJ,DAB,ISYD1,XLAMDH,ISYMLH,
     *                     WORK,LWORK,IDEL,ISYMD,FACTOR)
C
C     Written by Asger Halkier 10/2 - 1996
C     FACTOR added for RPA, Sonia & Francesca, 2010
C     Version: 1.0
C
C     Purpose: To calculate the second contribution to D2AIJ
C              from projection against doubles!
C
#include "implicit.h"
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION D2AIJ(*), DAB(*), XLAMDH(*), WORK(LWORK)
#include "priunit.h"      
#include "ccorb.h"
#include "ccsdsym.h"
C
      CALL QENTER('RPADAIJ22')
C
      ISYMB  = MULD2H(ISYMD,ISYMLH)
      ISYMA  = MULD2H(ISYMB,ISYD1)
      ISYMIJ = 1
C
      IF (ISYMA.NE.ISYAIJ) CALL QUIT('Symmetry mismatch in DAIJ22')
C
C-----------------------------------------------------------------
C     Calculate contraction of Y-intermediate & lamda hole matrix.
C-----------------------------------------------------------------
C
      IF (LWORK .LT. NVIR(ISYMA)) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', NVIR(ISYMA)
         CALL QUIT('Insufficient work space in DAIJ22')
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
C-------------------------------------
C     Calculate contribution to D2AIJ.
C-------------------------------------
C
      DO 100 ISYMJ = 1,NSYM
C
         ISYMI  = ISYMJ
         ISYMAI = MULD2H(ISYMI,ISYMA)
C
         DO 110 J = 1,NRHF(ISYMJ)
C
            I    = J
            NAIJ = IT2BCD(ISYMAI,ISYMJ) + NT1AM(ISYMAI)*(J - 1)
     *           + IT1AM(ISYMA,ISYMI)   + NVIR(ISYMA)*(I - 1) + 1
C
            CALL DAXPY(NVIR(ISYMA),-FACTOR,WORK,1,D2AIJ(NAIJ),1)
C
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('RPADAIJ22')
C
      RETURN
      END
!--
C  /* Deck rpadabi21 */
      SUBROUTINE rpaDABI21(D2ABI,ISYABI,DAB,ISYD1,XLAMDH,ISYMLH,
     &                     IDEL,ISYMD,FACTOR)
C
C     Written by Asger Halkier 11/2 - 1996
C     FACTOR added for RPA, Sonia e Francesca, 2010
C     Version: 1.0
C
C     Purpose: To calculate the first contribution to D2ABI
C              from projection against doubles!
C
#include "implicit.h"
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION D2ABI(*), DAB(*), XLAMDH(*)
#include "ccorb.h"
#include "ccsdsym.h"
C
      CALL QENTER('RPADABI21')
C
      ISYMI  = MULD2H(ISYMD,ISYMLH)
      ISYMAB = ISYD1
C
      IF (MULD2H(ISYMI,ISYMAB).NE.ISYABI)
     *  CALL QUIT('Symmetry mismatch in DABI21')
C
C--------------------------------
C     Calculate the contribution.
C--------------------------------
C
      DO 100 I = 1,NRHF(ISYMI)
C
         NDEI = IGLMRH(ISYMD,ISYMI) + NBAS(ISYMD)*(I - 1) 
     *        + IDEL - IBAS(ISYMD)
         NABI = ICKASR(ISYMAB,ISYMI) + NMATAB(ISYMAB)*(I - 1) + 1
C
         FACT = FACTOR*XLAMDH(NDEI)
C
         CALL DAXPY(NMATAB(ISYMAB),FACT,DAB,1,D2ABI(NABI),1)
C
  100 CONTINUE
C
      CALL QEXIT('RPADABI21')
C
      RETURN
      END
!---
C  /* Deck rpadiac21 */
      SUBROUTINE rpaDIAC21(D2IAC,ISYIAC,YMAT,ISYMY,XLAMDH,ISYMLH,
     *                  IDEL,ISYMD,FACTOR)
C
C     Written by Asger Halkier 21/2 - 1996
C     FACTOR introduced for RPA, Sonia & Francesca, 2010
C     Version: 1.0
C
C     Purpose: To calculate the first contribution to D2IAC
C              from projection against doubles!
C
#include "implicit.h"
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION D2IAC(*), YMAT(*), XLAMDH(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
      CALL QENTER('RPADIAC21')
C
      ISYMI  = MULD2H(ISYMD,ISYMLH)
      ISYMAC = ISYMY
C
      IF (MULD2H(ISYMI,ISYMAC).NE.ISYIAC)
     *  CALL QUIT('Symmetry mismatch in RPADIAC21')
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
         FACT = -XLAMDH(NDEI)*FACTOR
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
      CALL QEXIT('RPADIAC21')
C
      RETURN
      END
