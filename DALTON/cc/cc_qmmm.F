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
!
C
C
C***********************************************************************
C  /* Deck ccmm_addg */
      SUBROUTINE CCMM_ADDG(FOCK,WORK,LWORK)
C**************************************************************
C
C     Purpose: Construct the EFFective G solvent operator
C     and add to Fock
C     In the new QMMM implementation it replaces CCMM_RHSTG
C
C     The general strategy of the routine is the following:
C      IF (LGFILE==FALSE)
C      1) Get induced dipoles from file (preceeding HF/MM or CC/MM 
C         calc depending on first (HF/MM) or later(CC/MM) outer iteration
C         number.
C      2) Construct effective operator G^pol=-sum_a\mu^ind_aF^el_a (polarization
C         contribution)
C      3) Get G^elecstat from file (multipole integrals calculated in 
C        preceeding HF/MM calc)
C      ELSE IF (LFGILE==TRUE)
C     Get G operator from file dumped in the first iteration.
C     ENDIF
C     Add effective operator G to FOCK
C
C     The LGFILE is reset in TINE CC_QMMME
C     KS oct 09
C**************************************************************
C
#include "implicit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "qmmm.h"
#include "dummy.h"
#include "iratdef.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
#include "ccfield.h"
#include "exeinf.h"
#include "ccfdgeo.h"
#include "qm3.h"
#include "ccinftap.h"
#include "nuclei.h"
#include "orgcom.h"
C
      PARAMETER (ZERO = 0.0D0, MONE=-1.0D0 ,ONE = 1.0D0,
     &             IZERO = 0 , TWO = 2.0D0)
      DIMENSION WORK(LWORK),FOCK(*)
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE,TRIMAT,EXP1VL,LOCDEB
      PARAMETER (LOCDEB=.FALSE.)
C
      CHARACTER*8 LABEL
C
C
      CALL QENTER('CCMM_ADDG')
C-------------------------
C       Init parameters
C-------------------------
      NDIM=3*NNZAL
C--------------------------
C       Dynamical allocation
C--------------------------
C
      KTAO  =  1 
      KINT  = KTAO + NNBASX
      KWRK1 = KINT + N2BST(ISYMOP)
      LWRK1 = LWORK   - KWRK1

      IF (LWRK1.LT.0) CALL QUIT( 'Too little work in CCMM_ADDG' )
        
      CALL DZERO(WORK(KTAO),NNBASX)
      CALL DZERO(WORK(KINT),N2BST(ISYMOP))
C
C
C---------------------------------------
C     First the polarization contribution
C---------------------------------------

      IF (.NOT. LGFILE) THEN
        IF (IPOLTP .GT. 0 .AND. NNZAL.NE.0) THEN 
C
          CALL CCMM_POL(WORK(KTAO),WORK(KWRK1),LWRK1)
C
        END IF
C
C--------------------------------------
C     Then the electrostatic contribution 
C--------------------------------------
      CALL CCMM_ELSTAT(WORK(KTAO),WORK(KWRK1),LWRK1)
C
C
C
        IF (IPQMMM .GT. 14) THEN
          WRITE(LUPRI,*) 'NORM of G after electrostatic cont: ',
     *    DDOT(NNBASX,WORK(KTAO),1,WORK(KTAO),1)
          WRITE (LUPRI,*) ' G matrix: '
          CALL OUTPAK(WORK(KTAO),NBAST,1,LUPRI)
        ENDIF
C
C
C--------------------------------------------------
C     Add total contribution to effective AO fock matrix:
C--------------------------------------------------

        LABEL= 'GIVE INT'
        CALL CCMMINT(LABEL,-1,WORK(KINT),WORK(KTAO),
     &             IRREP,ISYM,IERR,WORK(KWRK1),LWRK1)
C
C Print G matrix to file
        IF (FFIRST ) THEN  !needed for PECC2 to ensure that we dump also the HF G to fil
          CALL PUT_TO_FILE('G_MATRIXHF',N2BST(ISYMOP),WORK(KINT))
          FFIRST=.FALSE.
        ELSE
          CALL PUT_TO_FILE('G_MATRIX',N2BST(ISYMOP),WORK(KINT))

C Put LGFILE to TRUE to indicate that we read it from file in future
C inner iterations. This keyword is reset in TINE QMMME
          LGFILE=.TRUE.
        ENDIF

      ELSE IF (LGFILE) THEN

        CALL GET_FROM_FILE('G_MATRIX',N2BST(ISYMOP),WORK(KINT))

      END IF

      IF (IPQMMM .GT.14) THEN
         CALL AROUND( 'CCMM cont. to AO matrix' )
         CALL OUTPUT(WORK(KINT),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      ENDIF
C
      IF (IPQMMM .GT.14) THEN
         CALL AROUND( 'Usual Fock AO matrix' )
         CALL CC_PRFCKAO(FOCK,ISYMOP)
      ENDIF
C
      CALL DAXPY(N2BST(ISYMOP),ONE,WORK(KINT),1,FOCK,1)
C
      IF (IPQMMM .GT.14) THEN
         CALL AROUND( 'Modified Fock AO matrix when leaving ADDG')
         CALL CC_PRFCKAO(FOCK,ISYMOP)
      ENDIF
C
      CALL QEXIT('CCMM_ADDG')
      END
C
C***********************************************************************
C  /* Deck ccmm_addg */
      SUBROUTINE CCMM_ADDGHF(FOCK,WORK,LWORK)
C**************************************************************
C
C     Purpose: Construct the Effective G solvent operator FROM
C     THE HF INDUCED DIPOLES and add to Fock
C
C     The general strategy of the routine is the following:
C      IF (LGFILE==FALSE)
C      1) Get induced dipoles from file (preceeding HF/MM) 
C         calc 
C      2) Construct effective operator G^pol=-sum_a\mu^ind_aF^el_a (polarization
C         contribution)
C      3) Get G^elecstat from file (multipole integrals calculated in 
C        preceeding HF/MM calc)
C      ELSE IF (LFGILE==TRUE)
C     Get G operator from file dumped in the first iteration.
C     ENDIF
C     Add effective operator G to FOCK
C
C     The LGFILE is reset in TINE CC_QMMME
C     KS oct 09
C**************************************************************
C
#include "implicit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "qmmm.h"
#include "dummy.h"
#include "iratdef.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
#include "ccfield.h"
#include "exeinf.h"
#include "ccfdgeo.h"
#include "qm3.h"
#include "ccinftap.h"
#include "nuclei.h"
#include "orgcom.h"
C
      PARAMETER (ZERO = 0.0D0, MONE=-1.0D0 ,ONE = 1.0D0,
     &             IZERO = 0 , TWO = 2.0D0)
      PARAMETER (D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0)
      DIMENSION WORK(LWORK),FOCK(*)
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE,TRIMAT,EXP1VL,LOCDEB
      PARAMETER (LOCDEB=.FALSE.)
C
      CHARACTER*8 LABEL
C
C
      CALL QENTER('CCMM_ADDGHF')
C-------------------------
C       Init parameters
C-------------------------
      NDIM=3*NNZAL
C--------------------------
C       Dynamical allocation
C--------------------------
C
      KTAO  =  1 
      KINT = KTAO + NNBASX
      KWRK1 = KINT + N2BST(ISYMOP)
      LWRK1   = LWORK   - KWRK1

      IF (LWRK1.LT.0) CALL QUIT( 'Too little work in CCMM_ADDGHF' )
        
      CALL DZERO(WORK(KTAO),NNBASX)
      CALL DZERO(WORK(KINT),N2BST(ISYMOP))
C
C---------------------------------------
C     First the polarization contribution
C---------------------------------------

      IF (.NOT. LHFFIL ) THEN
        IF (IPOLTP .GT. 0 .AND. NNZAL.NE.0) THEN 
C
          CALL CCMM_POL(WORK(KTAO),WORK(KWRK1),LWRK1)
C
        END IF
C
C--------------------------------------
C     Then the electrostatic contribution 
C--------------------------------------
      CALL CCMM_ELSTAT(WORK(KTAO),WORK(KWRK1),LWRK1)
C
C
C
        IF (IPQMMM .GT. 14 ) THEN
          WRITE(LUPRI,*) 'NORM of G after electrostatic cont: ',
     *    DDOT(NNBASX,WORK(KTAO),1,WORK(KTAO),1)
          WRITE (LUPRI,*) ' G matrix: '
          CALL OUTPAK(WORK(KTAO),NBAST,1,LUPRI)
        ENDIF
C--------------------------------------------------
C     Add total contribution to effective AO fock matrix:
C--------------------------------------------------

        LABEL= 'GIVE INT'
        CALL CCMMINT(LABEL,-1,WORK(KINT),WORK(KTAO),
     &             IRREP,ISYM,IERR,WORK(KWRK1),LWRK1)
C
C Print G matrix to file
        CALL PUT_TO_FILE('G_MATRIXHF',N2BST(ISYMOP),WORK(KINT))
C Put LHFFIL to TRUE to indicate that we read it from file in all future
C inner iterations. Unlike LGFILE this keyword is never reset.
        LHFFIL=.TRUE.

      ELSE IF (LHFFIL) THEN

        CALL GET_FROM_FILE('G_MATRIXHF',N2BST(ISYMOP),WORK(KINT))

      END IF

      IF (IPQMMM .GT.14 ) THEN
         CALL AROUND( 'CCMM cont. to AO matrix' )
         CALL OUTPUT(WORK(KINT),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      ENDIF
C
      IF (IPQMMM .GT.14) THEN
         CALL AROUND( 'Usual Fock AO matrix' )
         CALL CC_PRFCKAO(FOCK,ISYMOP)
      ENDIF
C
      CALL DAXPY(N2BST(ISYMOP),ONE,WORK(KINT),1,FOCK,1)
C
      IF (IPQMMM .GT.14) THEN
         CALL AROUND( 'Modified Fock AO matrix when leaving ADDGHF')
         CALL CC_PRFCKAO(FOCK,ISYMOP)
      ENDIF
C

      CALL QEXIT('CCMM_ADDGHF')
      END
C
C***********************************************************************
C
C  /* Deck ccmm_transformer */
      SUBROUTINE CCMM_TRANSFORMER(RHO1,RHO2,CTR1,CTR2,MODEL,
     *                          ISYMTR,LR,WORK,LWORK)
C
C-----------------------------------------------------------------------------
C
C   Construct effective operators from contracted densities and do 
C   response transformation. The type of density and transformation
C   is controlled by the LR flag. 
C
C   IF (LR NE '0')
C    1) Construct auxiliary density 
C    4) Construct effective operator,e.g G^C=-sum_a\∆mu^{DC}_a eps_a
C   ELSE IF (LR EQ '0')
C    Simply read GCC matrix from file constructed in the preceeding ADDG routine.
C   END IF
C   4) Do xi/eta-transformation depending on LR
C   5) Add to rho vector
C
C   Sneskov, oct 09 
C-----------------------------------------------------------------------------
C
#include "implicit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "qmmm.h"
#include "dummy.h"
#include "iratdef.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
#include "ccfield.h"
#include "exeinf.h"
#include "ccfdgeo.h"
#include "ccslvinf.h"
#include "qm3.h"
#include "ccinftap.h"
#include "nuclei.h"
#include "orgcom.h"

      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
      PARAMETER (HALF = 0.5D0)
      INTEGER NDIM
C
      DIMENSION WORK(LWORK),RHO1(*),RHO2(*),CTR1(*),CTR2(*)
      CHARACTER LISDUM*(2)
      INTEGER IDLDUM, ISYDUM
C
      CHARACTER MODEL*10 
      LOGICAL LOCDEB
      PARAMETER (LOCDEB=.FALSE.)
C
      CHARACTER*8 LABEL,LR*(1), LIST*(2)
C
C
      CALL QENTER('CCMM_TRANSFORMER')

      IF (IPQMMM.GT.10) THEN
        WRITE(LUPRI,*)'CCMM_TRANS : ISYMTR:', ISYMTR,
     &               ' and LR: ', LR
      ENDIF
C
C---------------------
C     Init parameters.
C---------------------
C
      NTAMP1  = NT1AM(ISYMTR)
      NTAMP2  = NT2AM(ISYMTR)
      NTAMP   = NTAMP1 + NTAMP2
      NDIM = 3*NNZAL
C
      IF (CCS)  CALL QUIT('CCMM_TRANSFORMER: Not implemented for CCS')
C
      IF (DISCEX) CALL QUIT('CCMM_TRANSFORM:Not implemented for DISCEX')
C fixed field approximation - no derivative response terms for fixed
C ccfld
      IF (CCFIXF .AND. (LR.NE.'0') ) THEN
        CALL QEXIT('CCMM_TRANSFORMER')
        RETURN 
      END IF
C
C------------------------------------
C     Dynamical allocation for CCMM :
C------------------------------------
C
      KDENS = 1
      KGMAT = KDENS + N2BST(ISYMTR)
      KETA = KGMAT  + N2BST(ISYMTR)
      KWRK = KETA   + NTAMP
      LWRK = LWORK  - KWRK

      IF (LWRK .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KWRK
         CALL QUIT( 'Too little work in CCMM_TRANSFORMER')
      END IF
C
C Initialize the working matrices
      CALL DZERO(WORK(KDENS),N2BST(ISYMTR))
      CALL DZERO(WORK(KETA),NTAMP) 
      CALL DZERO(WORK(KGMAT),N2BST(ISYMTR))
C
C Print section
      IF ( IPQMMM .GT. 10 .OR. LOCDEB) THEN
        RHO1N = DDOT(NT1AM(ISYMTR),RHO1,1,RHO1,1)
        RHO2N = DDOT(NT2AM(ISYMTR),RHO2,1,RHO2,1)
        WRITE(LUPRI,*) ' Norm of RHO1 in CCMM_TRANS on input:', RHO1N
        WRITE(LUPRI,*) ' Norm of RHO2 in CCMM_TRANS on input:', RHO2N
        RHO1N = DDOT(NT1AM(ISYMTR),CTR1,1,CTR1,1)
        RHO2N = DDOT(NT2AM(ISYMTR),CTR2,1,CTR2,1)
        WRITE(LUPRI,*) ' Norm af C1AM in CCMM_TRANS on input:', RHO1N
        WRITE(LUPRI,*) ' Norm af C2AM in CCMM_TRANS on input:', RHO2N
      ENDIF

      IF ( (LR.NE.'0' .AND. IPOLTP .GT. 0) .AND. (.NOT. HFFLD) ) THEN

C----------------------------
C     Construct auxiliary densities
C----------------------------
        CALL CCMM_D1AO(WORK(KDENS),CTR1,CTR2,MODEL,LR,
     *                 LISDUM,IDLDUM,ISYDUM,
     *                 WORK(KWRK),LWRK)
C----------------------------
C     Construct effective G operator
C----------------------------
        CALL CCMM_GEFF(WORK(KDENS),WORK(KGMAT),WORK(KWRK),LWRK)

C

      ELSE IF (LR .EQ. '0' ) THEN
        IF (HFFLD ) THEN
          CALL GET_FROM_FILE('G_MATRIXHF',N2BST(ISYMOP),WORK(KGMAT))
        ELSE
          CALL GET_FROM_FILE('G_MATRIX',N2BST(ISYMOP),WORK(KGMAT))
        END IF
      END IF
     
C Transform the effective operator according to LR
C---------------------------------------------------------------------
C             (LR.EQ.'L').OR.(LR.EQ.'F')                Calculate eta^Geff 
C             (LR.EQ.'R').OR.(LR.EQ.'P').OR.(LR.EQ.'0') Calculate xi^Geff
C---------------------------------------------------------------------
C
      KETA1 = KETA
      KETA2 = KETA1 + NTAMP1
C      

      IF ( (LR.EQ.'L') .OR. (LR.EQ.'F') ) THEN
         LABEL = 'GIVE INT'
         LIST  = 'L0'
         IDLINO = 1
C
         CALL CC_ETAC(ISYMTR,LABEL,WORK(KETA),
     *             LIST,IDLINO,0,WORK(KGMAT),WORK(KWRK),LWRK)
C
      ELSE IF ( (LR.EQ.'R') .OR. (LR.EQ.'P') .OR. (LR.EQ.'0')) THEN 
         LABEL = 'GIVE INT'
C
         CALL CC_XKSI(WORK(KETA),LABEL,ISYMTR,0,WORK(KGMAT),
     *                 WORK(KWRK),LWRK) 
C
         IF (LR .EQ. 'R') THEN
            CALL CCLR_DIASCL(WORK(KETA2),TWO,ISYMTR)
         END IF
      END IF
C

      IF(LOCDEB .OR. IPQMMM .GT. 14) THEN
         TAL1= DDOT(NT1AM(ISYMTR),WORK(KETA1),1,WORK(KETA1),1)
         TAL2= DDOT(NT2AM(ISYMTR),WORK(KETA2),1,WORK(KETA2),1)
         WRITE(LUPRI,*) 'Printing TRANS contribution. 
     &                     Norm2 of singles: ', TAL1,
     &                    'Norm2 of doubles: ', TAL2
      END IF

         CALL DAXPY(NT1AM(ISYMTR),ONE,WORK(KETA1),1,RHO1,1)
         CALL DAXPY(NT2AM(ISYMTR),ONE,WORK(KETA2),1,RHO2,1)
C
      IF (LOCDEB .OR. IPQMMM .GT. 14) THEN
         TAL1= DDOT(NT1AM(ISYMTR),RHO1,1,RHO1,1)
         TAL2= DDOT(NT2AM(ISYMTR),RHO2,1,RHO2,1)
         WRITE(LUPRI,*) 'Printing RHO:
     &                     Norm2 of singles: ', TAL1,
     &                    'Norm2 of doubles: ', TAL2
      END IF
C
      CALL QEXIT('CCMM_TRANSFORMER')
         
      END
C
C***********************************************************************
C  /* Deck ccmm_d1ao */
      SUBROUTINE CCMM_D1AO(AODEN,CTR1,CTR2,MODEL,LR,
     *                    LISTB,IDLSTB,ISYMTB,
     *                    WORK,LWORK)
C
C     Kristian Sneskov July 2009 inspired by CCX_D1AO
C
C     Purpose: To calculate contributions to the one electron 
C              auxiliary density matrix appearing in TINE CCMM_TRANSFORMER 
C              and return it backtransformed to AO-basis in AODEN.
C
C     IF(LR=R.OR.LR=F) -> Calc D^C
C     IF(LR=L.OR.LR=P) -> Calc D^B
C     Current models: CCS, CC2 and CCSD    
C     KS sep 10: Modified to calculate QR density when LR=Q
!     Note that the list arguments are only dummy when LR=R,L
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxash.h"
#include "maxorb.h"
#include "mxcent.h"
#include "aovec.h"
#include "iratdef.h"
#include "ccorb.h"
#include "infind.h"
#include "blocks.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccexci.h"
#include "ccsdio.h"
#include "distcl.h"
#include "cbieri.h"
#include "eribuf.h"
#include "cclr.h"
#include "qm3.h"
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION AODEN(*), WORK(LWORK), CTR1(*), CTR2(*)

      CHARACTER*(*) LISTB
      INTEGER IDLSTB, ISYMTB
      LOGICAL LOCDEB
      PARAMETER (LOCDEB=.FALSE.)
      CHARACTER MODEL*10,MODDUM*10
      CHARACTER LR*1
C
      CALL QENTER('CCMM_D1AO')
C
C---------------------------------------------
C     Find symmetry of left and right vectors. 
C---------------------------------------------
C
        ISYML = ILSTSYM('L0',1) ! Symmetry of multipliers and b-trial
        ISYMR = ILSTSYM('R0',1) ! symmetry of (t2)amplitudes
        IF (LR .EQ. 'Q' ) ISYML=ISYMTB

        IF (LR .NE. 'Q') THEN
          IF (ISYMR .NE. ISYML) 
     &      CALL QUIT('CCMM_D1AO: Density not total sym.')
        END IF
C
       
C-----------------------------------
C     Initial work space allocation.
C-----------------------------------
C
      KONEAI = 1
      KONEAB = KONEAI + NT1AMX
      KONEIJ = KONEAB + NMATAB(1)
      KONEIA = KONEIJ + NMATIJ(1)
      KL1AM  = KONEIA + NT1AMX
      KL2AM  = KL1AM  + NT1AM(ISYML) ! stores t2 when D^B and tbar2 when D^C
      KT1AM  = KL2AM  + NT2SQ(ISYML)
      KR2AM  = KT1AM  + NT1AM(ISYMR)
      KEND0  = KR2AM  + NT2AM(ISYMTR)
C
      KR2AMT = KEND0
      KEND1  = KR2AMT  + NT2AM(ISYMTR)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('CCMM_D1AO: Not enough memory')
      ENDIF
C
C--------------------------------------------
C     Initialize one electron density arrays.
C--------------------------------------------
C
      CALL DZERO(WORK(KONEAB),NMATAB(1))
      CALL DZERO(WORK(KONEIJ),NMATIJ(1))
      CALL DZERO(WORK(KONEAI),NT1AMX)
      CALL DZERO(WORK(KONEIA),NT1AMX)
C--------------------------------------------
C     Initialize working arrays.
C     NOTE if D^B -> KL1AM is just a dummy array
C--------------------------------------------
      CALL DZERO(WORK(KL1AM),NT1AM(ISYML))
      CALL DZERO(WORK(KL2AM),NT2SQ(ISYML))
      CALL DZERO(WORK(KT1AM),NT1AM(ISYMR))
      CALL DZERO(WORK(KR2AM),NT2AM(ISYMTR))
      CALL DZERO(WORK(KR2AMT),NT2AM(ISYMTR))
C
C----------------------
C     Read left vector.
C----------------------
C
C
      IF((LR.EQ.'R').OR.(LR.EQ.'F').OR.(LR.EQ.'Q')) THEN
         IOPT = 3
         CALL CC_RDRSP('L0',1,ISYML,IOPT,MODDUM,
     *             WORK(KL1AM),WORK(KR2AMT))
C when LR=Q overwrite with B vector
         IF (LR.EQ.'Q') THEN
           CALL CC_RDRSP(LISTB,IDLSTB,ISYMTB,IOPT,MODDUM,
     *             WORK(KL1AM),WORK(KR2AMT))
         END IF
C Read in t1 amplitudes
         IOPT = 1
         CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KT1AM),DUMMY)
C     Copy CTR2 to R2AM since we will scale diagonal
         CALL DCOPY(NT2AM(ISYMR),CTR2,1,WORK(KR2AM),1) 

         IF (ISYMR.EQ.1) CALL CCLR_DIASCL(WORK(KR2AM),TWO,ISYMR)
         IF (LOCDEB) THEN
            TAL1 = DDOT(NT1AMX,WORK(KL1AM),1,WORK(KL1AM),1)
            TAL2 = DDOT(NT2AMX,WORK(KR2AMT),1,WORK(KR2AMT),1)
            WRITE(LUPRI,*) 'CCMM_TGDEN: Norm of L1: ', TAL1, 
     &                     'Norm of L2: ', TAL2
         END IF
      ELSE IF ( (LR .EQ. 'L') .OR. (LR .EQ. 'P') ) THEN
         CALL DCOPY(NT2AM(ISYMTR),CTR2,1,WORK(KR2AMT),1)
         IOPT = 3
         CALL CC_RDRSP('R0',1,ISYMR,IOPT,MODDUM,
     *              WORK(KT1AM),WORK(KR2AM)) 
         IF (LOCDEB) THEN
            TAL1 = DDOT(NT1AMX,WORK(KT1AM),1,WORK(KT1AM),1)
            TAL2 = DDOT(NT2AMX,WORK(KR2AM),1,WORK(KR2AM),1)
            WRITE(LUPRI,*) 'CCMM_TGDEN: Norm of T1: ', TAL1, 
     &                     'Norm of T2: ', TAL2
         END IF
      END IF
C
C----------------------
C     Square up tbar_2 multipliers or b_2 elements
C----------------------
C
      CALL CC_T2SQ(WORK(KR2AMT),WORK(KL2AM),ISYML)
C
C---------------------------------------
C     Set up 2C-E of cluster amplitudes.
C---------------------------------------
C
C Note that diascl amplitudes are expected 
      IF (.NOT. MP2) THEN
         CALL DCOPY(NT2AM(ISYMR),WORK(KR2AM),1,WORK(KR2AMT),1) 
         IOPTTCME = 1
         CALL CCSD_TCMEPK(WORK(KR2AMT),1.0D0,ISYMR,IOPTTCME)
      ENDIF
C
C-------------------------------
C     Work space allocation one.
C     Note that D(ia) is stored transposed
C-------------------------------
C
      KXMAT  = KEND1
      KYMAT  = KXMAT  + NMATIJ(1)
      KEND2  = KYMAT  + NMATAB(1)
      LWRK2  = LWORK  - KEND1
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for 2. alloc. in CCMM_D1AO')
      ENDIF
      CALL DZERO(WORK(KXMAT),NMATIJ(1))
      CALL DZERO(WORK(KYMAT),NMATAB(1))
C
C-----------------------------------------------------------
C     Calculate X-intermediate of L2AM- and R2AM-amplitudes.
C-----------------------------------------------------------
C Note that KL2AM must be squared when input in these routines
C
      CALL CC_XI(WORK(KXMAT),WORK(KL2AM),ISYML,WORK(KR2AM),ISYMR,
     *           WORK(KEND2),LWRK2) 
C
C-----------------------------------------------------------
C     Calculate Y-intermediate of L2AM- and R2AM-amplitudes.
C-----------------------------------------------------------
C
      CALL CC_YI(WORK(KYMAT),WORK(KL2AM),ISYML,WORK(KR2AM),ISYMR,
     *           WORK(KEND2),LWRK2) 
C
C--------------------------------------------------------------
C     Construct <L2|[Emn,R2]|HF> contribution to DXaa and DXii.
C--------------------------------------------------------------
C
      CALL DCOPY(NMATAB(1),WORK(KYMAT),1,WORK(KONEAB),1)
      CALL CC_EITR(WORK(KONEAB),WORK(KONEIJ),WORK(KEND2),LWRK2,1)
      CALL DAXPY(NMATIJ(1),-ONE,WORK(KXMAT),1,WORK(KONEIJ),1)
C
      IF ( (LR .EQ. 'R') .OR. (LR .EQ. 'F') .OR. (LR .EQ. 'Q')) THEN
C
C--------------------------------------------------------------
C     Construct <HF|[Emn,R1]|HF> contribution.
C--------------------------------------------------------------
C 
         IF ( (LR .EQ. 'R') .OR. (LR .EQ. 'F') ) THEN
           CALL DAXPY(NT1AMX,TWO,CTR1,1,WORK(KONEIA),1)
         END IF
C
C
C--------------------------------------------------------------
C     Construct <L1|[Emn,R1]|HF> contribution to DXaa and DXii.
C--------------------------------------------------------------
         CALL CC_DXIJ(WORK(KONEIJ),WORK(KL1AM),ISYML,CTR1,ISYMR)
         CALL CC_DXAB(WORK(KONEAB),WORK(KL1AM),ISYML,CTR1,ISYMR)
C
C--------------------------------------------------------------
C     Construct <L1|[Eia,R2]|HF> contribution to DXia
C--------------------------------------------------------------
C
         CALL CC_DXIA12(WORK(KONEIA),WORK(KL1AM),ISYML,
     *         WORK(KR2AMT),ISYMR)
C
C----------------------------------------------------------
C     Construct <L2|[[Eia,R1],T2]|HF> contribution to DXia.
C----------------------------------------------------------
         KT2AM  = KEND0
         KEND3  = KT2AM  + NT2AM(ISYMOP) 
         LWRK3  = LWORK  - KEND3 
C read t2 amplitudes here rather than in CC_DXIA21
         IOPT = 2
         CALL CC_RDRSP('R0',0,1,IOPT,MODEL,DUMMY,WORK(KT2AM))
C

         IF (LWRK3 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT('Insufficient memory for 3. alloc. in CCMM_D1AO')
         ENDIF
         CALL CC_DXIA21(WORK(KONEIA),WORK(KL2AM),ISYML,
     *               CTR1,ISYMR,WORK(KT2AM),ISYMR,
     *               WORK(KEND3),LWRK3)

      ELSE IF ((LR .EQ. 'L') .OR. (LR .EQ. 'P')) THEN
C
C----------------------------------------------------------
C     Construct <B1|[Eai|HF> contribution to DXai.
C----------------------------------------------------------
         CALL DAXPY(NT1AMX,ONE,CTR1,1,WORK(KONEAI),1)

C--------------------------------------------------------------
C     Construct <B1|[Eia,T2]|HF> contribution to DXia
C--------------------------------------------------------------
         CALL CC_DXIA12(WORK(KONEIA),CTR1,ISYMTR,WORK(KR2AMT),ISYMR)
      END IF 
C
C--------------------------
C     Write out test norms.
C--------------------------
C
      IF (LOCDEB) THEN
         XIJ = DDOT(NMATIJ(1),WORK(KONEIJ),1,WORK(KONEIJ),1)
         XAB = DDOT(NMATAB(1),WORK(KONEAB),1,WORK(KONEAB),1)
         XAI = DDOT(NT1AMX,WORK(KONEAI),1,WORK(KONEAI),1)
         XIA = DDOT(NT1AMX,WORK(KONEIA),1,WORK(KONEIA),1)
         WRITE(LUPRI,*) 'Norms: DXIJ',XIJ
         WRITE(LUPRI,*) 'Norms: DXAB',XAB
         WRITE(LUPRI,*) 'Norms: DXAI',XAI
         WRITE(LUPRI,*) 'Norms: DXIA',XIA
      ENDIF
C
C----------------------------------
C     Calculate the lambda matrices.
C----------------------------------
C
      KLAMDP = KEND0
      KLAMDH = KLAMDP + NLAMDT
      KEND4  = KLAMDH + NLAMDT
      LWRK4  = LWORK  - KEND4 
      IF (LWRK4 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for 4. alloc. in CCMM_D1AO')
      ENDIF
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),WORK(KEND4),
     *            LWRK4)
C
C--------------------------------------------------------
C     Add the one electron density in the AO-basis.
C--------------------------------------------------------
CC
      ISDEN = 1
      CALL CC_DENAO(AODEN,ISDEN,WORK(KONEAI),WORK(KONEAB),
     *              WORK(KONEIJ),WORK(KONEIA),ISDEN,WORK(KLAMDP),1,
     *              WORK(KLAMDH),1,WORK(KEND4),LWRK4)
C
      CALL QEXIT('CCMM_D1AO')
C
   1  FORMAT(1x,A35,1X,E20.10)
      RETURN
      END
C
C***********************************************************************
C  /* Deck ccmm_f2dip */
      SUBROUTINE CCMM_F2DIP(ELF,INDMOM,WORK,LWORK)
C
C     April 2010, KS 
C
C     Purpose: Transform electric field to induced dipole
C     by either direct inversion (MMMAT) or by iterative means (MMITER).
C
#include "implicit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "dummy.h"
#include "iratdef.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
#include "ccfield.h"
#include "exeinf.h"
#include "ccfdgeo.h"
#include "ccslvinf.h"
#include "ccinftap.h"
#include "nuclei.h"
#include "orgcom.h"
#include "qm3.h"
#include "qmmm.h"
C
      DIMENSION WORK(LWORK)
      DOUBLE PRECISION ELF(*),INDMOM(*)
      LOGICAL LOCDEB,FNDLAB
      PARAMETER (LOCDEB=.FALSE.)
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
      INTEGER NDIM

      CALL QENTER('CCMM_F2DIP')

      NDIM=3*NNZAL
      NDIM2P = NDIM * (NDIM+1)/2!packed response matrix

      KRELAY=1
      KWRK1=KRELAY+NDIM2P 
      LWRK1 =LWORK-KWRK1
      IF (LWRK1.LT.0) CALL QUIT('Too little work in CCMM_F2DIP')

      CALL DZERO(WORK(KRELAY),NDIM2P)
C
C
      IF (MMMAT) THEN
         LUQMMM = -1
         IF(LUQMMM .LT. 0) THEN
            CALL GPOPEN(LUQMMM,'QMMMIM','UNKNOWN','SEQUENTIAL',
     &                    'UNFORMATTED',IDUMMY,.FALSE.)
         ENDIF
         REWIND(LUQMMM)
C      
         IF (FNDLAB('QQMMMMAT',LUQMMM)) THEN
            CALL READT(LUQMMM,NDIM2P,WORK(KRELAY))
         ELSE 
            CALL QUIT('Problem reading the relay 
     &                 matrix from QMMMIM file')
         ENDIF

         CALL GPCLOSE(LUQMMM,'KEEP')
C
         IF (IPQMMM .GE. 5) THEN
            WRITE(LUPRI,*) 'The relay matrix is read from file'
            CALL OUTPAK(WORK(KRELAY),NDIM,1,LUPRI)
         ENDIF
C
         IF (IPRINT .GT. 1) THEN
            WRITE(LUPRI,*)
            WRITE(LUPRI,1051)
            WRITE(LUPRI,1050)
            WRITE(LUPRI,1051)
            WRITE(LUPRI,*)
         ENDIF
C
         IF(LOCDEB) THEN
            WRITE(LUPRI,*) 'F vector to be transformed to dipole: '
            DO I=1,NDIM
               WRITE(LUPRI,*) ELF(I)
            END DO
         END IF

C
C calculate the corresponding induced moments
         CALL DSPMV('L', NDIM, ONE, WORK(KRELAY), ELF, 1,ZERO, INDMOM,1)

C 
      ELSE IF (MMITER) THEN
C convert the F vector into induced moments in an iterative fashion               
         IOPT=1
         CALL F2QMMM(ELF,NNZAL,INDMOM,WORK(KRELAY),
     *                 LWORK,IOPT,IPQMMM)
C
         IF (IPQMMM .GT. 1) THEN
            WRITE(LUPRI,*)
            WRITE(LUPRI,1030)
            WRITE(LUPRI,*)
            WRITE(LUPRI,1000)
            WRITE(LUPRI,1010)
            WRITE(LUPRI,1000)
         ENDIF
C
         IINIM = 0

         DO I=1,MMCENT
            IF (ZEROAL(I) .EQ. -1) THEN
               DIPX = 0.0D0
               DIPY = 0.0D0
               DIPZ = 0.0D0
            ELSE
               DIPX = INDMOM(IINIM) 
               DIPY = INDMOM(IINIM +1)
               DIPZ = INDMOM(IINIM +2)
               IINIM = IINIM +3
            ENDIF
            IF (IPQMMM .GT. 1) THEN
               WRITE(LUPRI,1020) I, DIPX, DIPY, DIPZ
            END IF
         END DO

         IF(IPQMMM .GT. 1) THEN
            WRITE(LUPRI,1000)
            WRITE(LUPRI,*)
         ENDIF
      END IF
C
C
 1050 FORMAT('  Induced dipole moments   ')
 1051 FORMAT(2X,'=',22('-'),'=',2X)
 1000 FORMAT(1X,51('='))
 1010 FORMAT(' | Site |     X     |     Y     |     Z     |')
 1020 FORMAT(1X,I6,3(4X,F10.6))
 1030 FORMAT(' Total induced dipole moments: ')
 
      CALL QEXIT('CCMM_F2DIP')
      END
CC
C***********************************************************************
C  /* Deck ccmm_elstat */
      SUBROUTINE CCMM_ELSTAT(NSAO,WORK,LWORK)
C
C     April 2010, KS 
C
C     Purpose: Calculate electrostatic contribution to G operator. In
C     old qm3 formulation we refered to this contribution 
C     as Ns and it is included in the
C     effective operator (G) by setting OLDTG=TRUE
C     We get the multipole integrals from the HF/MM
C     calculation available from file.
C     We do all the multipole corrections simultaneosly. We will analyze
C     the HF/MM results if we want to consider each multipole contribution
C     separately.
C
C     NSAO is a vector containing the sum (of all orders) of 
C     the multipole integrals. This is the one we need at output. Note that
C     it is returned packed i.e. LENGTH=NNBASX 
C     Note that all signs are carried in the integrals themselves (cf.
C     sirqmmm.F) and daxpy is always with a +
C
#include "implicit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "dummy.h"
#include "iratdef.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
#include "ccfield.h"
#include "exeinf.h"
#include "ccfdgeo.h"
#include "ccslvinf.h"
#include "ccinftap.h"
#include "nuclei.h"
#include "orgcom.h"
#include "qm3.h"
#include "qmmm.h"
C
      DOUBLE PRECISION NSAO(*)
      DIMENSION WORK(LWORK)
      PARAMETER (ONE=1.0D0)
      PARAMETER (MONE=-1.0D0)
      CALL QENTER('CCMM_ELSTAT')

      KMULTI = 1
      KWRK = KMULTI + NNBASX
      LWRK = LWORK - KWRK
      IF (LWRK.LT.0) CALL QUIT('Too little work in CCMM_ELSTAT')
      CALL DZERO(WORK(KMULTI),NNBASX)

      IF (NMULT .GE. 0) THEN
C
         LUQMMM = -1
         CALL GPOPEN(LUQMMM,'MU0INT','OLD',' ',
     &            'UNFORMATTED',IDUMMY,.FALSE.)
         REWIND (LUQMMM)
         CALL READT(LUQMMM,NNBASX,WORK(KMULTI))
         CALL GPCLOSE(LUQMMM,'KEEP')
         IF (IPQMMM .GT. 14) THEN
            WRITE(LUPRI,*) 'CCMM_ELSTAT: Print the stored 
     &         charge integrals:'
            CALL OUTPAK(WORK(KMULTI),NBAST,1,LUPRI)
         END IF

         CALL DAXPY(NNBASX,ONE,WORK(KMULTI),1,NSAO,1)

      END IF
C 2) The dipole contribution
      IF (NMULT .GE. 1) THEN
         CALL DZERO(WORK(KMULTI),NNBASX)
         CALL GPOPEN(LUQMMM,'MU1INT','OLD',' ',
     &            'UNFORMATTED',IDUMMY,.FALSE.)
         REWIND (LUQMMM)
         CALL READT(LUQMMM,NNBASX,WORK(KMULTI))
         CALL GPCLOSE(LUQMMM,'KEEP')
         IF (IPQMMM .GT. 14) THEN
            WRITE(LUPRI,*) 'CCMM_ELSTAT: Print the stored 
     &         dipole integrals:'
            CALL OUTPAK(WORK(KMULTI),NBAST,1,LUPRI)
         END IF

         CALL DAXPY(NNBASX,ONE,WORK(KMULTI),1,NSAO,1)

      END IF
C 3) The quadrupole contribution
      IF (NMULT .GE. 2) THEN
         CALL DZERO(WORK(KMULTI),NNBASX)
         CALL GPOPEN(LUQMMM,'MU2INT','OLD',' ',
     &            'UNFORMATTED',IDUMMY,.FALSE.)
         REWIND (LUQMMM)
         CALL READT(LUQMMM,NNBASX,WORK(KMULTI))
         CALL GPCLOSE(LUQMMM,'KEEP')
         IF (IPQMMM .GT. 14) THEN
            WRITE(LUPRI,*) 'CCMM_ELSTAT: Print the stored 
     &        quadrupole integrals:'
            CALL OUTPAK(WORK(KMULTI),NBAST,1,LUPRI)
         END IF

         CALL DAXPY(NNBASX,ONE,WORK(KMULTI),1,NSAO,1)

      END IF
C

      CALL QEXIT('CCMM_ELSTAT')
      END
C***********************************************************************
C  /* Deck ccmm_dens2f */
      SUBROUTINE CCMM_DENS2F(DENS,ELF,WORK,LWORK)
C
C     April 2010, KS 
C
C     Purpose: Transform aux density to electric field.
C     As of now (april-10) it only handles the aux density
C     and not the construction of the normal induced dipoles.
C     These are instead handled in FSTATIC
C     KEPSAO contains the packed electric field integrals.
C     KEPSx/y/z  contains the unpacked field integrals.
C
#include "implicit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "dummy.h"
#include "iratdef.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
#include "ccfield.h"
#include "exeinf.h"
#include "ccfdgeo.h"
#include "ccslvinf.h"
#include "ccinftap.h"
#include "nuclei.h"
#include "orgcom.h"
#include "qm3.h"
#include "qmmm.h"
#include "cclr.h"
C
      DIMENSION WORK(LWORK), DENS(*), ELF(*)
      INTEGER NDIM
      LOGICAL LOCDEB, LSKIP 
      CHARACTER*8 LABEL
C
      CALL QENTER('CCMM_DENS2F')
C
      NDIM=3*NNZAL

C Memory allocation
      KEPSAO = 1
      KEPSx = KEPSAO + 3*NNBASX
      KEPSy = KEPSx + N2BST(ISYMTR)
      KEPSz = KEPSy + N2BST(ISYMTR)
      KWRK  = KEPSz + N2BST(ISYMTR)
      LWRK  = LWORK   - KWRK
      IF (LWRK.LT.0) CALL QUIT('Too little work in CCMM_DENS2F')
C

      LRI=1
C
      DO 200 I = 1,MMCENT

        LSKIP = .FALSE.

        CALL CCMM_EPSAO(WORK(KEPSAO),I,LSKIP,WORK(KWRK),LWRK)

        IF (LSKIP) THEN
          IF (ZEROAL(I) .EQ. -1) THEN
            GOTO 200
          ELSE  
            LRI = LRI + 3
            GOTO 200
          END IF
        END IF

C Unpack the integrals
        LABEL = 'GIVE INT'

        CALL DZERO(WORK(KEPSx),3*N2BST(ISYMTR))

        CALL CCMMINT(LABEL,LUQMMM,WORK(KEPSx),WORK(KEPSAO),
     &          IRREP,ISYM,IERR,WORK(KWRK),LWRK)
        CALL CCMMINT(LABEL,LUQMMM,WORK(KEPSy),WORK(KEPSAO+NNBASX),
     &             IRREP,ISYM,IERR,WORK(KWRK),LWRK)
        CALL CCMMINT(LABEL,LUQMMM,WORK(KEPSz),WORK(KEPSAO+2*NNBASX),
     &          IRREP,ISYM,IERR,WORK(KWRK),LWRK)
C 
C Do the trace operation of the product matrices
C
         PRODx = DDOT(N2BST(ISYMTR),WORK(KEPSx),1,DENS,1)
         PRODy = DDOT(N2BST(ISYMTR),WORK(KEPSy),1,DENS,1)
         PRODz = DDOT(N2BST(ISYMTR),WORK(KEPSz),1,DENS,1)

         ELF(LRI + 0) = PRODx
         ELF(LRI + 1) = PRODy
         ELF(LRI + 2) = PRODz 
         
         LRI = LRI + 3

 200  CONTINUE


      IF (IPQMMM .GE. 5) THEN
         WRITE(LUPRI,*) 'Transformed electric fields and NDIM', NDIM
         CALL OUTPUT(ELF,0,NDIM,1,1,NDIM,1,1,LUPRI)
         WRITE(LUPRI,*)
      END IF
      CALL QEXIT('CCMM_DENS2F')

      END
C***********************************************************************
C  /* Deck ccmm_pol */
      SUBROUTINE CCMM_POL(GPOL,WORK,LWORK)
C
C     April 2010, KS 
C
C     Purpose: Polarization contribution to effective operator G
C     The first time we enter this routine we read the converged 
C     induced dipoles from the preceeding HF/MM calculation.
C     In later iterations we read the induced dipoles written in 
C     CC_QMMME
C     G=-sum_a mu^ind_a eps_a
C
C
#include "implicit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "dummy.h"
#include "iratdef.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
#include "ccfield.h"
#include "exeinf.h"
#include "ccfdgeo.h"
#include "ccslvinf.h"
#include "ccinftap.h"
#include "nuclei.h"
#include "orgcom.h"
#include "qm3.h"
#include "qmmm.h"
#include "cclr.h"
C
      DIMENSION WORK(LWORK)
      DOUBLE PRECISION GPOL(*)
      INTEGER NDIM
      LOGICAL LSKIP

      CALL QENTER('CCMM_POL')
      
      NDIM=3*NNZAL
      
      KEPSAO  = 1
      KINDIP = KEPSAO + 3*NNBASX
      KWRK  = KINDIP + NDIM
      LWRK = LWORK - KWRK

      IF (LWRK.LT.0) CALL QUIT('Too little work in CCMM_POL')
      CALL DZERO(WORK(KINDIP),NDIM)
      
C Read the induced moments from file 
      CALL GET_FROM_FILE_1('INDUCED_DIPOLES',NNZAL,WORK(KINDIP))
C
C-------------------
C     Print section:
C-------------------
C
      
      IF (IPQMMM .GT. 10) THEN
         WRITE(LUPRI,'(//,A)')
     *        ' +============================================'
     *       //'==============+'
         WRITE(LUPRI,'(A)')
     *        ' | MMcent | <mu^ind>_x  | <mu^ind>_y      |'
     *       //' <mu^ind>_z      |'
         WRITE(LUPRI,'(1X,A)')
     *        '+============================================'
     *       //'==============+'

         LRI = 0

         DO 100 I = 1, MMCENT

            IF (ZEROAL(I) .EQ. -1) GOTO 100


            WRITE(LUPRI,'(1X,A,I3,A,F16.10,A,
     &                       F16.10,A,F16.10,A)')
     *         '| ', I,'|', WORK(KINDIP+LRI+0),' |',
     *              WORK(KINDIP+LRI+1),' |', 
     *              WORK(KINDIP+LRI+2),' |'
                    WRITE(LUPRI,'(1X,A)')
     *           '+---------------------------------------------'
     *           //'-------------+'
            LRI = LRI + 3
 100     CONTINUE 
         WRITE(LUPRI,'(//,A)')
      END IF

C
C---------------------------------------
C Calculate electric field due to electrons at s in ao basis 
C---------------------------------------
      
      LRI = 0
C
      DO 200 I = 1, MMCENT

         LSKIP =.FALSE.

         CALL CCMM_EPSAO(WORK(KEPSAO),I,LSKIP,WORK(KWRK),LWRK)

         IF (LSKIP) THEN
            IF (ZEROAL(I) .EQ. -1) THEN
               GOTO 200
            ELSE  
               LRI = LRI + 3
               GOTO 200
            END IF
         END IF

         FACx =  - WORK(KINDIP + LRI +0 )  
         FACy =  - WORK(KINDIP + LRI +1 )  
         FACz =  - WORK(KINDIP + LRI +2 )  
C
         CALL DAXPY(NNBASX,FACx,WORK(KEPSAO),
     *                   1,GPOL,1)
C
         CALL DAXPY(NNBASX,FACy,WORK(KEPSAO+NNBASX),
     *                     1,GPOL,1)
C 
         CALL DAXPY(NNBASX,FACz,WORK(KEPSAO+2*NNBASX),
     *                  1,GPOL,1)

         LRI = LRI + 3
 200  CONTINUE
C
      IF (IPQMMM.GT.14) THEN
         WRITE (LUPRI,*) ' NORM of G matrix after pol. contr.: ',
     *   DDOT(NNBASX,GPOL,1,GPOL,1)
         WRITE (LUPRI,'(/A)') ' G matrix: '
         CALL OUTPAK(GPOL,NBAST,1,LUPRI)
      ENDIF

      CALL QEXIT('CCMM_POL')
      END

C***********************************************************************
C  /* Deck ccmm_fstatic */
      SUBROUTINE CCMM_FSTATIC(ELF,AODEN,WORK,LWORK)
C
C     April 2010, KS 
C
C     Purpose: Calculate the total static field as
C     F^sta=F^mul+F^el+F^nuc
C
C
#include "implicit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "dummy.h"
#include "iratdef.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
#include "ccfield.h"
#include "exeinf.h"
#include "ccfdgeo.h"
#include "ccslvinf.h"
#include "ccinftap.h"
#include "nuclei.h"
#include "orgcom.h"
#include "qm3.h"
#include "qmmm.h"
#include "cclr.h"
C
      DIMENSION WORK(LWORK)
      DIMENSION ELF(*),AODEN(*)
      INTEGER NDIM
      LOGICAL LOCDEB, LSKIP
      PARAMETER (LOCDEB=.FALSE.)
      CHARACTER LABEL*8
   
      CALL QENTER('CCMM_FSTATIC')

      NDIM=3*NNZAL
      ISYMOP=1

      KEPSAO=1
      KEPSx=KEPSAO+3*NNBASX
      KEPSy=KEPSx+N2BST(ISYMOP)
      KEPSz=KEPSy+N2BST(ISYMOP)
      KWRK=KEPSz+N2BST(ISYMOP)
      LWRK=LWORK-KWRK

      LRI = 1 ! Row-index in supervector

      DO 100 I = 1,MMCENT

C First we calculate the field due to the QM electrons by calculating
C the electric field operator and dotting with the AO CC density to 
C obtain the CC expectation value
        LSKIP = .FALSE.

        CALL CCMM_EPSAO(WORK(KEPSAO),I,LSKIP,WORK(KWRK),LWRK)

        IF (LSKIP) THEN
           IF (ZEROAL(I) .EQ. -1) THEN
              GOTO 100
           ELSE  
              LRI = LRI + 3
              GOTO 100
           END IF
        END IF
C

C Unpack the integrals
         CALL DZERO(WORK(KEPSx),3*N2BST(ISYMOP))
         LABEL = 'GIVE INT'
         CALL CCMMINT(LABEL,LUQMMM,WORK(KEPSx),WORK(KEPSAO),
     *                  IRREP,ISYM,IERR,WORK(KWRK),LWRK)
         CALL CCMMINT(LABEL,LUQMMM,WORK(KEPSy),
     *                   WORK(KEPSAO+NNBASX),IRREP,ISYM,IERR,
     *                   WORK(KWRK),LWRK)
         CALL CCMMINT(LABEL,LUQMMM,WORK(KEPSz),
     *                   WORK(KEPSAO+2*NNBASX),IRREP,ISYM,IERR,
     *                   WORK(KWRK),LWRK)
C Dot with the AO density

         EXELCO = DDOT(N2BST(ISYMOP),AODEN,1,WORK(KEPSx),1)
         EYELCO = DDOT(N2BST(ISYMOP),AODEN,1,WORK(KEPSy),1)
         EZELCO = DDOT(N2BST(ISYMOP),AODEN,1,WORK(KEPSz),1)

         IF (LOCDEB) THEN
            WRITE(LUPRI,*) 'At MM site ', I, ' the field due to the QM',
     &      ' electrons: ',EXELCO,EYELCO,EZELCO
         END IF
 
C Add to the electric field from the multipoles 
         ELF(LRI+0) = ELF(LRI+0) + EXELCO
         ELF(LRI+1) = ELF(LRI+1) + EYELCO
         ELF(LRI+2) = ELF(LRI+2) + EZELCO

C Second we add the field due to the distributed multipoles
         CALL CCMM_FMUL(ELF,LRI,I)

C Finally add the field from the QM nuclei
         CALL CCMM_FNUC(ELF,LRI,I)              

         LRI = LRI + 3

 100  CONTINUE !Sum over all sites - done with the static field
 
      CALL QEXIT('CCMM_FSTATIC')
      END
C**************************************************************************************
C  /* Deck put_to_file */
      SUBROUTINE PUT_TO_FILE(FLNAME,NULOOP,DDATA)
C**************************************************************
C
#include "implicit.h"
#include "dummy.h"
C
      CHARACTER*(*) FLNAME
      INTEGER   NMBU,NULOOP
      DIMENSION DDATA(*)
C
      NMBU = -1
      CALL GPOPEN(NMBU,FLNAME,'UNKNOWN',' ','FORMATTED',IDUMMY,.FALSE.)
C
      REWIND (NMBU)
      DO 820 L = 1,NULOOP
        WRITE(NMBU,'(I5,3E25.15)') L,DDATA(L)
  820 CONTINUE
C
      CALL GPCLOSE(NMBU,'KEEP')
C
      END
C
C**************************************************************
C  /* Deck get_from_file */
      SUBROUTINE GET_FROM_FILE(FLNAME,NULOOP,DDATA)
C**************************************************************
C
#include "implicit.h"
#include "dummy.h"
C
      CHARACTER*(*) FLNAME
      INTEGER   NMBU,NULOOP
      DIMENSION DDATA(*)
C
      NMBU = -1
      CALL GPOPEN(NMBU,FLNAME,'UNKNOWN',' ','FORMATTED',IDUMMY,.FALSE.)
C
      REWIND (NMBU)
      DO 820 L = 1,NULOOP
        READ(NMBU,'(I5,3E25.15)') LK,DDATA(L)
  820 CONTINUE
C
      IF (LK.NE.NULOOP) THEN
        CALL QUIT('Problem in dimension in GET_FROM_FILE')
      ENDIF

      CALL GPCLOSE(NMBU,'KEEP')
C
      END
C
C***********************************************************************
C  /* Deck ccmm_epsao */
      SUBROUTINE CCMM_EPSAO(EPSAO,SITE,LSKIP,WORK,LWORK)
C
C     July 2010, KS 
C
C     Purpose: Calculate the electric field at a site SITE due to
C     the electrons. LSKIP keeps track of whether a site is skipped
C     either due to having no polarizability or because of distance
C     cutoff. The vector is returned in the AO basis.
C     Note we are doing everything integral-direct to avoid
C     a possible storage-problem when working with large systems.

#include "implicit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "dummy.h"
#include "iratdef.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
#include "ccfield.h"
#include "exeinf.h"
#include "ccfdgeo.h"
#include "ccslvinf.h"
#include "ccinftap.h"
#include "nuclei.h"
#include "orgcom.h"
#include "qmmm.h"
#include "qm3.h"
#include "cclr.h"

      PARAMETER (D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0)
      DIMENSION WORK(LWORK)
      DIMENSION EPSAO(*)
      INTEGER SITE
      LOGICAL LSKIP
      LOGICAL TOFILE,TRIMAT,EXP1VL,LOCDEB 
      PARAMETER (LOCDEB=.FALSE.)
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      CHARACTER*8 LABINT(9*MXCENT)

      KWRK = 1
      LWRK = LWORK - KWRK

      IF (ZEROAL(SITE) .EQ. -1) THEN 
        LSKIP = .TRUE.
        RETURN
      END IF

C Backup dipole origin
      OBKPX = DIPORG(1)
      OBKPY = DIPORG(2)
      OBKPZ = DIPORG(3)

      DIPORG(1) = MMCORD(1,SITE)
      DIPORG(2) = MMCORD(2,SITE)
      DIPORG(3) = MMCORD(3,SITE)
 
      DIST2 = (MMCORD(1,SITE)-QMCOM(1))**2 + 
     &        (MMCORD(2,SITE)-QMCOM(2))**2 + 
     &        (MMCORD(3,SITE)-QMCOM(3))**2  
      DIST = SQRT(DIST2)

      IF (DIST .GT. RCUTMM) THEN
        LSKIP =.TRUE.
        IF (LOCDEB) THEN
          WRITE(LUPRI,*) 'Skipping site in mu^ind vector 
     &                    in CCMM_POL', SITE
        ENDIF
        RETURN
      ENDIF
      
      CALL DZERO(EPSAO,3*NNBASX)
      KPATOM = 0
      NOCOMP = 3
      TOFILE = .FALSE.
      TRIMAT = .TRUE.
      EXP1VL = .FALSE.   

      RUNQM3 = .TRUE.
      CALL GET1IN(EPSAO,'NEFIELD',NOCOMP,WORK(KWRK),
     *           LWRK,LABINT,INTREP,INTADR,SITE,TOFILE,KPATOM,
     *           TRIMAT,DUMMY,EXP1VL,DUMMY,IPQMMM)
      RUNQM3 = .FALSE.

      IF (QMDAMP) THEN
        IF ( (IDAMP .EQ. 3) .AND. (NQMNUC .NE. NUCIND) ) THEN
          CALL QUIT('ERROR in no. of assigned QM 
     &                  polarizabilities')
        ENDIF 
        IF ( (IDAMP .EQ. 1) .OR. (IDAMP .EQ. 3) ) THEN
          DIST = 9.99D+99
          MHIT = 0
          DO 205 M=1,NUCIND
            DISTC = (DIPORG(1)-CORD(1,M))**2 +
     &              (DIPORG(2)-CORD(2,M))**2 +
     &              (DIPORG(3)-CORD(3,M))**2
            IF (DISTC .LE. DIST) THEN
              DIST = DISTC
              MHIT = M
            ENDIF
  205     CONTINUE
        ELSE IF (IDAMP .EQ. 2) THEN
          DIST = (DIPORG(1)-QMCOM(1))**2 +
     &           (DIPORG(2)-QMCOM(2))**2 +
     &           (DIPORG(3)-QMCOM(3))**2
        ENDIF
        DIST = SQRT(DIST)
C
        IF (IDAMP .EQ. 3) THEN
          IF (IPOLTP .EQ. 2) THEN
            TEMPI =  D3I*(POLMM(1,SITE)+POLMM(4,SITE)+POLMM(6,SITE))
          ELSE IF (IPOLTP .EQ. 1) THEN
            TEMPI = POLIMM(SITE)
          ENDIF 
          TEMP = (TEMPI*QMPOL(MHIT))**(D6I)
          SIJ = 2.1304D0*DIST/TEMP
          DFACT = 1.0D0 - (1.0D0+SIJ+0.50D0*SIJ*SIJ)*exp(-SIJ)
        ELSE
          DFACT = (1-exp(-ADAMP*DIST))**3
        ENDIF
        CALL DSCAL(3*NNBASX,DFACT,EPSAO,1)
      ENDIF

C Print section
      IF ( (IPQMMM.GT.15) .OR. (LOCDEB) ) THEN
        WRITE(LUPRI,*) 'Packed AO electric field_x matrix due to QM
     &  electrons at MM site ',SITE
        CALL OUTPAK(EPSAO(1),NBAST,1,LUPRI)

        WRITE(LUPRI,*) 'Packed AO electric field_y matrix due to QM
     &  electrons at MM site ',SITE
        CALL OUTPAK(EPSAO(1+NNBASX),NBAST,1,LUPRI)

        WRITE(LUPRI,*) 'Packed AO electric field_z matrix due to QM
     &  electrons at MM site ',SITE
        CALL OUTPAK(EPSAO(1+2*NNBASX),NBAST,1,LUPRI)
      ENDIF
C Restore dipole origin
      DIPORG(1) = OBKPX
      DIPORG(2) = OBKPY
      DIPORG(3) = OBKPZ

      END
C**************************************************************
C  /* Deck cc_qmmme */
      SUBROUTINE CC_QMMME(ECCMM,EPOL,EEL,EELMUL,EREP,AODEN,WORK,LWORK)
C**************************************************************
C
C     KS, 09
C     Purpose: Updates the induced dipoles besides calculating E(QM/MM)
C
C     In new QMMM implementation (NYQMMM) the strategy is as follows
C     1) Calculate F^sta=F^el+F^mul+F^nuc
C     2) Calculate induced dipole by transforming with relay 
C        matrix (MMMAT) or by iterative means (MMITER)
C     3) Put induced moments to file
C     4) Calculate Epol=-1/2sum_a\mu^ind_aF^sta_a
C     5) Calculate Eelesta=-sum_s<N_s>
C
C     Also the LGFILE keyword is reset such that a new G matrix is
C     written to file when TGT is called next time.
C
C**************************************************************
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "qmmm.h"
#include "qm3.h"
#include "iratdef.h"
#include "maxorb.h"
#include "orgcom.h"
#include "nuclei.h"
#include "ccinftap.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
#include "ccfield.h"
#include "exeinf.h"
#include "ccfdgeo.h"
#include "thrzer.h"

C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
      PARAMETER (HALF = 1.0D0/2.0D0)
      PARAMETER (D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0)
C
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE, TRIMAT, EXP1VL, LOCDEB, FNDLAB, EXCENT
      PARAMETER (LOCDEB=.FALSE.)
      DIMENSION WORK(LWORK),AODEN(*)
      CHARACTER LABEL*8
      INTEGER NTOTI
C
      CALL QENTER('CC_QMMME')
C------------------------------------------------------------------
C       Init parameters
C------------------------------------------------------------------

      ISYMOP = 1
C
C------------------------------------------------------------------
C       Dynamical allocation for CCMM 
C------------------------------------------------------------------
C
      NDIM = 3*NNZAL
      LGFILE=.FALSE. ! reset G matrix construction
      KELF = 1
      KINDIP = KELF + NDIM
      KNSAO = KINDIP + NDIM
      KINT = KNSAO + NNBASX
      KWRK1 = KINT + N2BST(ISYMOP)
      LWRK1 = LWORK - KWRK1
C
      CALL DZERO(WORK(KELF),NDIM)
      CALL DZERO(WORK(KINDIP),NDIM)
C
C
      IF (LWRK1.LT.0) CALL QUIT( 'Too little work in CC_QMMME')
C 
C SPLDIP option not implemented for CC
      IF (SPLDIP) WRITE(LUPRI,*) 'WARNING: SPLDIP option
     &   not implemented for CC/MM. Keyword ignored!'
      IF (FIXDIP) WRITE(LUPRI,*) 'WARNING: FIXDIP option
     &   not implemented for CC/MM. Keyword ignored!'
      IF (LGSPOL) CALL QUIT( 'LGSPOL not implemented for CC/MM')

C
      IF ( ( IPOLTP .GT. 0 ) .AND. (.NOT. HFFLD) ) THEN
         CALL CCMM_FSTATIC(WORK(KELF),AODEN,WORK(KWRK1),LWRK1)
C
CC  Transform the field to induced dipoles
         CALL CCMM_F2DIP(WORK(KELF),WORK(KINDIP),WORK(KWRK1),LWRK1)
C
CC Write the updated induced moments to file 
         CALL PUT_TO_FILE_1('INDUCED_DIPOLES',NNZAL,WORK(KINDIP))
      ELSE IF ( HFFLD .AND. (IPOLTP .GT. 0 ) ) THEN
         CALL CCMM_FSTATIC(WORK(KELF),AODEN,WORK(KWRK1),LWRK1)
         CALL GET_FROM_FILE_1('INDUCED_DIPOLES',NNZAL,WORK(KINDIP))
CC
      END IF

      ECCMM = ZERO
      EEL = ZERO
      EELMUL = ZERO
      EPOL = ZERO
      ETEMP = ZERO
C-------------------------------------
C     Add up the energy contributions:
C     1) Epol:
C-------------------------------------

      ETEMP = DDOT(NDIM,WORK(KINDIP),1,WORK(KELF),1)
     
      EPOL = -HALF*ETEMP

C Now add the electrostatic energy contribution
C we do this by getting the multipole integrals from file
C and taking the expectation value. The integrals are at our disposal
C from the preceeding HF/MM calculation handled by sirius
C--------
C 2) The electrostatic contribution:
C--------
C
      CALL DZERO(WORK(KNSAO),NNBASX)
      CALL CCMM_ELSTAT(WORK(KNSAO),WORK(KWRK1),LWRK1)

C Unpack the multipole integrals and dot with the density to
C obtain the expectation value

      LABEL = 'GIVE INT'
      CALL CCMMINT(LABEL,LUQMMM,WORK(KINT),WORK(KNSAO),
     *               IRREP,ISYM,IERR,WORK(KWRK1),LWRK1)
         
      EELMUL = DDOT(N2BST(ISYMOP),WORK(KINT),1,AODEN,1)

C Finally we add the interaction energy between QM nuclei and MM multipoles. This 
C was found in the previous SCF/MM calculation and is readily available 
C in the COMMON block. 
      EEL = EELMUL + ENUMUL

      ECCMM = EPOL + EEL
      IF (LOCDEB .OR. IPQMMM .GT. 4) THEN
         WRITE(LUPRI,*) 'EELMUL: ', EELMUL
         WRITE(LUPRI,*) 'EELSTAT: ', EEL
         WRITE(LUPRI,*) 'EPOL: ', EPOL
         WRITE(LUPRI,*) 'ECCMM: ', ECCMM
         WRITE(LUPRI,*) 'The nuclear-multipole contribution: ', ENUMUL
      END IF
C
C
C-------------
C 3) Edisp
C-------------
C
C sne - New qmmm not implemented for dispersion yet
C sne - Anyhow it is just a static energy contribution (for now) so
C straightforward to do - copy how it is done in dft
C      ECCMM = ECCMM + ECLVDW
C
C-----------------------------------------------
C sne New qmmm not implemented for explicit repulsion yet
C 4) E_repulsion
C
C      CALL CCMM_REP2(EREP,AODEN,WORK(KWRK1),LWRK1)
C      ECCMM = ECCMM + EREP
C
C-----------------------------------------------
C   Writing out the final energy contributions one by one:
C-----------------------------------------------
C
 1050 FORMAT('  Induced dipole moments   ')
 1051 FORMAT(2X,'=',22('-'),'=',2X)
 1000 FORMAT(1X,51('='))
 1020 FORMAT(1X,I6,3(4X,F10.6))
 1010 FORMAT(' | Site |     X     |     Y     |     Z     |')
 1030 FORMAT(' Total induced dipole moments: ')

      CALL QEXIT('CC_QMMME')
C
      END
C
C***********************************************************************
C  /* Deck ccmm_fmul */
      SUBROUTINE CCMM_FMUL(ELF,LRI,SITEI)
C
C     Aug 2010, KS 
C
C     Purpose: Calculate the electric field due to the MM multipoles
C     at SITEI
C     

#include "implicit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "dummy.h"
#include "iratdef.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
#include "ccfield.h"
#include "exeinf.h"
#include "ccfdgeo.h"
#include "ccslvinf.h"
#include "ccinftap.h"
#include "nuclei.h"
#include "orgcom.h"
#include "qm3.h"
#include "qmmm.h"
#include "cclr.h"

      DIMENSION ELF(*)
      INTEGER SITEI, SITEJ, LRI
      LOGICAL LOCDEB, EXCENT
      PARAMETER (LOCDEB=.FALSE.)

      CALL QENTER('CCMM_FMUL')

      DO 200 SITEJ=1,MMCENT

        IF (SITEJ .NE. SITEI) THEN     
          DIST2 = (MMCORD(1,SITEI)-MMCORD(1,SITEJ))**2 +
     &            (MMCORD(2,SITEI)-MMCORD(2,SITEJ))**2 +
     &            (MMCORD(3,SITEI)-MMCORD(3,SITEJ))**2
          DIST = SQRT(DIST2)

          IF(DIST .GT. RCUTMM) THEN
            IF (LOCDEB) THEN
              WRITE(LUPRI,*)
     &           'Skipping element in T2 CC', SITEI,SITEJ
            ENDIF
            GOTO 200
          ENDIF
               
          EXCENT=.FALSE.
               
C Check if any centers were in the same LoProp calculation
          IF (NEWEXC) THEN
            DO L = 1, NEXLST
              IF (EXLIST(1,SITEI) .EQ. EXLIST(L,SITEJ)) EXCENT = .TRUE.
            ENDDO
          ELSE
            DO L = 1, NEXLST
              IF (EXLIST(L,SITEI) .EQ. EXLIST(1,SITEJ)) EXCENT = .TRUE.
            ENDDO
          ENDIF

C Form F vector due to the permament moments

          IF (.NOT. (EXCENT) ) THEN

C              A) The point-charge contribution
            IF((ABS(MUL0MM(SITEJ)) .GT. THRMM) 
     &                .AND. (NMULT.GE.0) ) THEN
        
              ELFLDX = 0.0D0 
              ELFLDY = 0.0D0
              ELFLDZ = 0.0D0 

              CALL GET_CHARGE_ELFLD(MUL0MM(SITEJ),
     *                 MMCORD(1,SITEJ),MMCORD(2,SITEJ),MMCORD(3,SITEJ),
     *                 MMCORD(1,SITEI),MMCORD(2,SITEI),MMCORD(3,SITEI),
     *                 ELFLDX,ELFLDY,ELFLDZ)
              ELF(LRI+0) = ELF(LRI+0)+ELFLDX
              ELF(LRI+1) = ELF(LRI+1)+ELFLDY
              ELF(LRI+2) = ELF(LRI+2)+ELFLDZ

            ENDIF
            
C              B) The dipole contribution
            IF (NMULT .GE. 1) THEN
              ELFLDX = 0.0D0
              ELFLDY = 0.0D0
              ELFLDZ = 0.0D0
              CALL GET_DIPOLE_ELFLD(MUL1MM(1,SITEJ),
     *                  MUL1MM(2,SITEJ),MUL1MM(3,SITEJ),
     *                  MMCORD(1,SITEJ),MMCORD(2,SITEJ),MMCORD(3,SITEJ),
     *                  MMCORD(1,SITEI),MMCORD(2,SITEI),MMCORD(3,SITEI),
     *                  ELFLDX,ELFLDY,ELFLDZ)
              ELF(LRI+0) = ELF(LRI+0)+ELFLDX
              ELF(LRI+1) = ELF(LRI+1)+ELFLDY
              ELF(LRI+2) = ELF(LRI+2)+ELFLDZ
            ENDIF

C              C) The quadrupole contribution
            IF (NMULT .GE. 2) THEN
              ELFLDX = 0.0D0
              ELFLDY = 0.0D0
              ELFLDZ = 0.0D0
              CALL GET_QUADRUPOLE_ELFLD(
     *                  MUL2MM(1,SITEJ),MUL2MM(2,SITEJ),MUL2MM(3,SITEJ),
     *                  MUL2MM(4,SITEJ),MUL2MM(5,SITEJ),MUL2MM(6,SITEJ),
     *                  MMCORD(1,SITEJ),MMCORD(2,SITEJ),MMCORD(3,SITEJ),
     *                  MMCORD(1,SITEI),MMCORD(2,SITEI),MMCORD(3,SITEI),
     *                  ELFLDX,ELFLDY,ELFLDZ)
              ELF(LRI+0) = ELF(LRI+0)+ELFLDX
              ELF(LRI+1) = ELF(LRI+1)+ELFLDY
              ELF(LRI+2) = ELF(LRI+2)+ELFLDZ
            ENDIF 
          ENDIF
        END IF
 200  CONTINUE

      CALL QEXIT('CCMM_FMUL')
      END
C***********************************************************************
C  /* Deck ccmm_fnuc */
      SUBROUTINE CCMM_FNUC(ELF,LRI,SITEI)
C
C     Aug 2010, KS 
C
C     Purpose: Calculate the electric field due to the QM nuclei
C     at SITEI
C     

#include "implicit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "dummy.h"
#include "iratdef.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
#include "ccfield.h"
#include "exeinf.h"
#include "ccfdgeo.h"
#include "ccslvinf.h"
#include "ccinftap.h"
#include "nuclei.h"
#include "orgcom.h"
#include "qm3.h"
#include "qmmm.h"
#include "cclr.h"

      DIMENSION ELF(*)
      INTEGER SITEI, LRI, SITEJ
      LOGICAL LOCDEB
      PARAMETER (LOCDEB=.FALSE.)
      PARAMETER (D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0)

      TJKX = 0.0D0
      TJKY = 0.0D0
      TJKZ = 0.0D0

      DO 145 SITEJ=1,NUCIND ! Dalton keyword for nr of QM nuclei
        ELFLDX = 0.0D0
        ELFLDY = 0.0D0
        ELFLDZ = 0.0D0
        CALL GET_CHARGE_ELFLD(CHARGE(SITEJ),
     *             CORD(1,SITEJ),CORD(2,SITEJ),CORD(3,SITEJ),
     *             MMCORD(1,SITEI),MMCORD(2,SITEI),MMCORD(3,SITEI),
     *             ELFLDX,ELFLDY,ELFLDZ)

        IF (QMDAMP) THEN
          IF ((IDAMP .EQ. 3).AND.(NQMNUC .NE. NUCIND)) THEN
            CALL QUIT('ERROR in no. of assigned QM 
     &                  polarizabilities')
          ENDIF
          IF ( (IDAMP .EQ. 1) .OR. (IDAMP .EQ. 3) ) THEN
            DIQM = 9.99D+99
            MHIT = 0
            DO 155 M=1,NUCIND
              DIQMC = (MMCORD(1,SITEI)-CORD(1,M))**2 +
     &                (MMCORD(2,SITEI)-CORD(2,M))**2 +
     &                (MMCORD(3,SITEI)-CORD(3,M))**2
              IF (DIQMC .LE. DIQM) THEN
                DIQM = DIQMC
                MHIT = M
              ENDIF
  155       CONTINUE
          ELSE IF (IDAMP .EQ. 2) THEN
            DIQM = (MMCORD(1,SITEI)-QMCOM(1))**2 +
     &             (MMCORD(2,SITEI)-QMCOM(2))**2 +
     &             (MMCORD(3,SITEI)-QMCOM(3))**2
          ENDIF
               DIQM = SQRT(DIQM)

          IF (IDAMP .EQ. 3) THEN
            IF (IPOLTP .EQ. 2) THEN
              TEMPI =D3I*(POLMM(1,SITEI)+POLMM(4,SITEI)+POLMM(6,SITEI))
            ELSE IF (IPOLTP .EQ. 1) THEN
              TEMPI = POLIMM(SITEI)
            ENDIF
            TEMP = (TEMPI*QMPOL(MHIT))**(D6I)
            SIJ = 2.1304D0*DIQM/TEMP
            DFACT=1.0D0-(1.0D0+SIJ+0.50D0*SIJ*SIJ)*exp(-SIJ)
          ELSE
            DFACT = (1-exp(-ADAMP*DIQM))**3
          ENDIF

          ELFLDX = ELFLDX*DFACT
          ELFLDY = ELFLDY*DFACT
          ELFLDZ = ELFLDZ*DFACT
        END IF

        TJKX = TJKX + ELFLDX
        TJKY = TJKY + ELFLDY
        TJKZ = TJKZ + ELFLDZ
 145  CONTINUE
      ELF(LRI+0)=ELF(LRI+0)+TJKX
      ELF(LRI+1)=ELF(LRI+1)+TJKY
      ELF(LRI+2)=ELF(LRI+2)+TJKZ


      IF (LOCDEB) THEN
        WRITE(LUPRI,*) 'Nuclear field: CC ', TJKX,TJKY,TJKZ
      ENDIF

      END
C***********************************************************************
C
C  /* Deck ccmm_qrtransformer */
      SUBROUTINE CCMM_QRTRANSFORMER(RHO1,RHO2,ISYRES,
     *                          LISTB,IDLSTB,ISYMTB,
     *                          LISTC,IDLSTC,ISYMTC,
     *                          MODEL,RSPTYP,WORK,LWORK)
C
C-----------------------------------------------------------------------------
C
C   Construct effective operators from contracted densities and do 
C   response transformation. This is the QR analog to the TINE
C   CCMM_LRTRANSFORMER. 
C   For B and G we have 3 contributions 
C   For K transformation 2 contributions

C   A: Lagrangian multipliers to be used for F transformation
C   B: B trial vector
C   C: C trial vector
C
C   The strategy is as follows:
C    1) Construct auxiliary density 
C    2) Transform density to electric field
C    3) Transform electric field to induced dipoles
C    4) Construct effective operator. G^C=-sum_a\∆mu^C_a eps_a
C    5) Do response-transformation the type of which determined by QR flag
C    6) Add to rho vector
C
C   KS, aug 10 
C   KS, Sep 10: Modified to be called from FMATNEW to generate TPA pseudo
C   B matrix contributions
C-----------------------------------------------------------------------------
C
#include "implicit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "qmmm.h"
#include "dummy.h"
#include "iratdef.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
#include "ccfield.h"
#include "exeinf.h"
#include "ccfdgeo.h"
#include "ccslvinf.h"
#include "qm3.h"
#include "ccinftap.h"
#include "nuclei.h"
#include "orgcom.h"

      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, IZERO = 0 , TWO = 2.0D0)
      PARAMETER (HALF = 0.5D0)
C
      DIMENSION WORK(LWORK),RHO1(*),RHO2(*)
C
      CHARACTER*(*) LISTB, LISTC
      INTEGER IDLSTB, IDLSTC
      INTEGER ISYCB, ISYRES, ISYMTB, ISYMTC
      CHARACTER MODEL*10 
      LOGICAL LOCDEB, LSAME
      SAVE LOCDEB
      DATA LOCDEB /.FALSE./
      CHARACTER*(2) LISTL
C
      CHARACTER*8 LABEL,DENTYP*(1),RSPTYP*(1) 
C
C
      CALL QENTER('CCMM_QRTRANSFORMER')
      LISTL='L0'

      IF (CCS)  CALL QUIT('CCMM_QRTRANS: Not implemented for CCS')
C
      IF (DISCEX) CALL QUIT('CCMM_QRTRANS:Not implemented for DISCEX')

      IF (IPQMMM.GT.10 .OR. LOCDEB) THEN
        WRITE(LUPRI,*) 'CCMM_QRTRANSFORMER: RSPTYP: ',RSPTYP
        WRITE(LUPRI,*)'CCMM_QRTRANSFORMER: ISYRES= ',ISYRES 
        WRITE(LUPRI,*)'CCMM_QRTRANSFORMER: ISYMTB= ',ISYMTB 
        WRITE(LUPRI,*)'CCMM_QRTRANSFORMER: ISYMTC= ',ISYMTC 
        WRITE(LUPRI,*)'CCMM_QRTRANSFORMER: LISTB= ',LISTB 
        WRITE(LUPRI,*)'CCMM_QRTRANSFORMER: LISTC= ',LISTC 
        WRITE(LUPRI,*)'CCMM_QRTRANSFORMER: IDLSTB= ',IDLSTB 
        WRITE(LUPRI,*)'CCMM_QRTRANSFORMER: IDLSTC= ',IDLSTC 
        CALL FLSHFO(LUPRI)
      ENDIF

      IF (IPOLTP .EQ. 0) THEN
        WRITE(LUPRI,*) 'No QM/MM polarization so no contribution
     &                 from quadratic response transformer.'
        CALL QEXIT('CCMM_QRTRANSFORMER')
        RETURN
      END IF
C
C     Check the symmetry
      ISYCB=MULD2H(ISYMTB,ISYMTC)
      IF (ISYCB .NE. ISYRES) THEN
        CALL QUIT('Symmetry problem in CCMM_QRTRANSFORMER')
      END IF
C
C Check what kind of response.
C Find out which density and hence which operator to construct
C For F start with left density and change flag when done with this cont.
      IF ( (RSPTYP.EQ.'G') .OR. (RSPTYP.EQ.'B') ) THEN
        DENTYP='R'
      ELSE IF (RSPTYP.EQ.'K' .OR. (RSPTYP .EQ. 'F')) THEN
        DENTYP='L'
      END IF
C
C---------------------
C     Init parameters.
C---------------------
C
      NB1  = NT1AM(ISYMTB)
      NB2  = NT2AM(ISYMTB)
      NDIM = 3*NNZAL
C
C------------------------------------
C     Dynamical allocation for CCMM :
C------------------------------------
      KB1=1
      KB2=KB1+NB1
      KDENS=KB2+NB2
      KGBMAT=KDENS+N2BST(ISYMTB)
      KWRK1=KGBMAT+N2BST(ISYMTB)
      LWRK1=LWORK-KWRK1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KWRK1
         CALL QUIT( 'Too little work in CCMM_QRTRANSFORMER, 1')
      END IF
C
C Initialize the working matrices
      CALL DZERO(WORK(KB1),NB1)
      CALL DZERO(WORK(KB2),NB2)
      CALL DZERO(WORK(KDENS),N2BST(ISYMTB))
      CALL DZERO(WORK(KGBMAT),N2BST(ISYMTB))
C
C Read C vector from file
      IOPT = 3
      CALL CC_RDRSP(LISTB,IDLSTB,ISYMTB,IOPT,MODEL,
     *            WORK(KB1),WORK(KB2))
C
C Print section
      IF ( IPQMMM .GT. 10 .OR. LOCDEB) THEN
        RHO1N = DDOT(NT1AM(ISYCB),RHO1,1,RHO1,1)
        RHO2N = DDOT(NT2AM(ISYCB),RHO2,1,RHO2,1)
        WRITE(LUPRI,*)'Norm of RHO1 in CCMM_QRTRANS on input:,1', RHO1N
        WRITE(LUPRI,*)'Norm of RHO2 in CCMM_QRTRANS on input:,1', RHO2N
        RHO1N = DDOT(NB1,WORK(KB1),1,WORK(KB1),1)
        RHO2N = DDOT(NB2,WORK(KB2),1,WORK(KB2),1)
        WRITE(LUPRI,*)'Norm af B1AM in CCMM_QRTRANS on input:,1', RHO1N
        WRITE(LUPRI,*)'Norm af B2AM in CCMM_QRTRANS on input:,1', RHO2N
      ENDIF
C----------------------------
C     Construct auxiliary densities
C----------------------------
      CALL CCMM_D1AO(WORK(KDENS),WORK(KB1),WORK(KB2),MODEL,DENTYP,
     *                 LISTC,IDLSTC,ISYMTC,
     *                 WORK(KWRK1),LWRK1)

C----------------------------
C     Construct effective operator
C----------------------------
      CALL CCMM_GEFF(WORK(KDENS),WORK(KGBMAT),WORK(KWRK1),LWRK1)
C
C Transform the effective operator according to RESPTYP flag
      IF(LOCDEB .OR. IPQMMM .GT. 14) THEN
         TAL1= DDOT(N2BST(ISYCB),WORK(KGBMAT),1,WORK(KGBMAT),1)
         WRITE(LUPRI,*) 'Print Norm2 GBMAT: ', TAL1
      END IF
C
C Dynamical memory allocation to store the response transformed operator
      KETA1=KWRK1
      KETA2=KETA1+NT1AM(ISYCB)
      KWRK2=KETA2+NT2AM(ISYCB)
      LWRK2=LWORK-KWRK2

      LABEL='GIVE INT'
      IF ( (RSPTYP.EQ.'G') .OR. (RSPTYP.EQ.'F') ) THEN
        CALL CCLR_FA(LABEL,ISYMTB,LISTC,IDLSTC,LISTL,0,
     *       WORK(KGBMAT),WORK(KETA1),LWRK2)
      ELSE IF (RSPTYP.EQ.'B') THEN
        CALL CCCR_AA(LABEL,ISYMTB,LISTC,IDLSTC,WORK(KGBMAT),
     *       WORK(KETA1),LWRK2)
        CALL CCLR_DIASCL(WORK(KETA2),TWO,ISYCB)
      ELSE IF (RSPTYP.EQ.'K') THEN
        CALL CC_ETAC(ISYMTB,LABEL,WORK(KETA1),LISTC,IDLSTC,0,
     *                WORK(KGBMAT),WORK(KWRK2),LWRK2)
      END IF
        

C Print section
      IF (LOCDEB .OR. IPQMMM .GT. 14) THEN
         TAL1= DDOT(NT1AM(ISYCB),WORK(KETA1),1,WORK(KETA1),1)
         TAL2= DDOT(NT2AM(ISYCB),WORK(KETA2),1,WORK(KETA2),1)
         WRITE(LUPRI,*) 'Printing transformed G^B*C contribution. 
     &                     Norm2 of singles: ', TAL1,
     &                    'Norm2 of doubles: ', TAL2
      END IF

C Note that if B=C then we just multiply the already found contribution by 2
C If we are dealing with the pseudo B matrix contribution no contributions
C will be the same and we prepare to make a right density.
      IF (RSPTYP .EQ. 'F') THEN 
        FACTKS=ONE
        LSAME=.FALSE.
        DENTYP='R'
      ELSE 
        LSAME= (LISTC.EQ.LISTB .AND. IDLSTC.EQ.IDLSTB)
        IF (LSAME) THEN
          FACTKS=TWO
        ELSE
          FACTKS=ONE
        END IF  
      END IF

      CALL DAXPY(NT1AM(ISYCB),FACTKS,WORK(KETA1),1,RHO1,1)
      CALL DAXPY(NT2AM(ISYCB),FACTKS,WORK(KETA2),1,RHO2,1)
C
      IF (LOCDEB .OR. IPQMMM .GT. 14) THEN
        TAL1= DDOT(NT1AM(ISYCB),RHO1,1,RHO1,1)
        TAL2= DDOT(NT2AM(ISYCB),RHO2,1,RHO2,1)
        WRITE(LUPRI,*) 'Printing RHO:
     &                     Norm2 of singles: ', TAL1,
     &                    'Norm2 of doubles: ', TAL2
      END IF

C Then the second contribution ^BG*C. 
      IF ( .NOT. LSAME) THEN

C Initial parameters
        NC1  = NT1AM(ISYMTC)
        NC2  = NT2AM(ISYMTC)

C Renaming (for readability) and resetting
        KC1=KB1
        KC2=KC1+NC1
        KDENS=KC2+NC2
        KGCMAT=KDENS+N2BST(ISYMTC)
        KWRK1=KGCMAT+N2BST(ISYMTC)
        LWRK1=LWORK-KWRK1

C Zero the workspace
        CALL DZERO(WORK(KC1),NC1)
        CALL DZERO(WORK(KC2),NC2)
        CALL DZERO(WORK(KDENS),N2BST(ISYMTC))
        CALL DZERO(WORK(KGCMAT),N2BST(ISYMTC))
C
C Read C vector from file
        IOPT = 3
        CALL CC_RDRSP(LISTC,IDLSTC,ISYMTC,IOPT,MODEL,
     *            WORK(KC1),WORK(KC2))
C
C Print section
        IF ( IPQMMM .GT. 10 .OR. LOCDEB) THEN
          RHO1N = DDOT(NT1AM(ISYCB),RHO1,1,RHO1,1)
          RHO2N = DDOT(NT2AM(ISYCB),RHO2,1,RHO2,1)
          WRITE(LUPRI,*)'Norm of RHO1 in CCMM_QRTRANS on input,2:',RHO1N
          WRITE(LUPRI,*)'Norm of RHO2 in CCMM_QRTRANS on input,2:',RHO2N
          RHO1N = DDOT(NC1,WORK(KC1),1,WORK(KC1),1)
          RHO2N = DDOT(NC2,WORK(KC2),1,WORK(KC2),1)
          WRITE(LUPRI,*)'Norm af C1AM in CCMM_QRTRANS on input,2:',RHO1N
          WRITE(LUPRI,*)'Norm af C2AM in CCMM_QRTRANS on input,2:',RHO2N
        ENDIF
C----------------------------
C     Construct auxiliary densities
C----------------------------
        CALL CCMM_D1AO(WORK(KDENS),WORK(KC1),WORK(KC2),MODEL,DENTYP,
     *                 LISTB,IDLSTB,ISYMTB,
     *                 WORK(KWRK1),LWRK1)
C
C----------------------------
C     Construct effective operator 
C----------------------------
        CALL CCMM_GEFF(WORK(KDENS),WORK(KGCMAT),WORK(KWRK1),LWRK1)

C Transform the effective operator according to RSPTYP flag
C
C Dynamical memory allocation to store the response transformed operator
        KETA1=KWRK1
        KETA2=KETA1+NT1AM(ISYCB)
        KWRK2=KETA2+NT2AM(ISYCB)
        LWRK2=LWORK-KWRK2

        IF (LWRK2 .LT. 0) THEN
          WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KWRK2
          CALL QUIT( 'Too little work in CCMM_QRTRANSFORMER, 2')
        END IF
        LABEL='GIVE INT'
        IF (RSPTYP.EQ.'G') THEN
          CALL CCLR_FA(LABEL,ISYMTC,LISTB,IDLSTB,LISTL,0,
     *         WORK(KGCMAT),WORK(KETA1),LWRK2)
        ELSE IF (RSPTYP.EQ.'B') THEN
          CALL CCCR_AA(LABEL,ISYMTC,LISTB,IDLSTB,WORK(KGCMAT),
     *       WORK(KETA1),LWRK2)
          CALL CCLR_DIASCL(WORK(KETA2),TWO,ISYCB)
        ELSE IF ( (RSPTYP.EQ.'K') .OR. (RSPTYP.EQ.'F')) THEN
          CALL CC_ETAC(ISYMTC,LABEL,WORK(KETA1),LISTB,IDLSTB,0,
     *                WORK(KGCMAT),WORK(KWRK2),LWRK2)
          DENTYP='Q' ! only effect for F term since K contains no more conts
        END IF

C Print section
        IF(LOCDEB .OR. IPQMMM .GT. 14) THEN
          TAL1= DDOT(NT1AM(ISYCB),WORK(KETA1),1,WORK(KETA1),1)
          TAL2= DDOT(NT2AM(ISYCB),WORK(KETA2),1,WORK(KETA2),1)
          WRITE(LUPRI,*) 'Printing transformed G^C*B contribution. 
     &                     Norm2 of singles: ', TAL1,
     &                    'Norm2 of doubles: ', TAL2
        END IF

        CALL DAXPY(NT1AM(ISYCB),ONE,WORK(KETA1),1,RHO1,1)
        CALL DAXPY(NT2AM(ISYCB),ONE,WORK(KETA2),1,RHO2,1)
C
        IF (LOCDEB .OR. IPQMMM .GT. 14) THEN
          TAL1= DDOT(NT1AM(ISYCB),RHO1,1,RHO1,1)
          TAL2= DDOT(NT2AM(ISYCB),RHO2,1,RHO2,1)
          WRITE(LUPRI,*) 'Printing RHO:
     &                     Norm2 of singles: ', TAL1,
     &                    'Norm2 of doubles: ', TAL2
        END IF
      END IF

C Finally the third contribution GBC
      IF ( (RSPTYP.EQ.'G') .OR. (RSPTYP.EQ.'B') 
     &                                   .OR. (RSPTYP.EQ.'F')) THEN
C Note we do not delete C vector from work
        KGBC=KDENS+N2BST(ISYCB)
        KWRK1=KGBC+N2BST(ISYCB)
        LWRK1=LWORK-KWRK1     

C Reset the primary workspace
        CALL DZERO(WORK(KDENS),N2BST(ISYCB))
        CALL DZERO(WORK(KGBC),N2BST(ISYCB))

C Construct the quadratic CCMM auxiliary density
        IF ((RSPTYP .EQ. 'G') .OR. (RSPTYP .EQ. 'B') ) THEN
          CALL CCMM_D2AO(WORK(KDENS),ISYRES,
     *             LISTB,IDLSTB,ISYMTB,
     *             LISTC,IDLSTC,ISYMTC,
     *             MODEL,WORK(KWRK1),LWRK1)
        ELSE IF (RSPTYP .EQ. 'F') THEN
          CALL CCMM_D1AO(WORK(KDENS),WORK(KC1),WORK(KC2),MODEL,DENTYP,
     *                 LISTB,IDLSTB,ISYMTB,
     *                 WORK(KWRK1),LWRK1)
        END IF

C Construct effective operator
        IF(LOCDEB .OR. IPQMMM .GT. 14) THEN
          TAL1= DDOT(N2BST(ISYCB),WORK(KGBC),1,WORK(KGBC),1)
          WRITE(LUPRI,*) 'Printing {GBC} before build. 
     &                     Norm2 of operator: ', TAL1
        END IF
        CALL CCMM_GEFF(WORK(KDENS),WORK(KGBC),WORK(KWRK1),LWRK1)
C
        IF(LOCDEB .OR. IPQMMM .GT. 14) THEN
          TAL1= DDOT(N2BST(ISYCB),WORK(KGBC),1,WORK(KGBC),1)
          WRITE(LUPRI,*) 'Printing {GBC} after build. 
     &                     Norm2 of operator: ', TAL1
        END IF

C Then do response transformation 
        KETA1=KWRK1
        KETA2=KETA1+NT1AM(ISYCB)
        KWRK2=KETA2+NT2AM(ISYCB)
        LWRK2=LWORK-KWRK2
        IF (LWRK2 .LT. 0) THEN
          WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KWRK2
          CALL QUIT( 'Too little work in CCMM_QRTRANSFORMER, 2')
        END IF
        CALL DZERO(WORK(KETA1),NT1AM(ISYCB)+NT2AM(ISYCB))
    
        LABEL='GIVE INT'
        IF ( (RSPTYP .EQ. 'G') .OR. (RSPTYP .EQ. 'F')) THEN
          CALL CC_ETAC(ISYCB,LABEL,WORK(KETA1),LISTL,0,0,
     *                WORK(KGBC),WORK(KWRK2),LWRK2)
        ELSE IF (RSPTYP .EQ. 'B') THEN
          CALL CC_XKSI(WORK(KETA1),LABEL,ISYCB,0,WORK(KGBC),
     *                 WORK(KWRK2),LWRK2)
          CALL CCLR_DIASCL(WORK(KETA2),TWO,ISYCB)
        END IF

      LOCDEB=.TRUE.
        IF(LOCDEB .OR. IPQMMM .GT. 14) THEN
          TAL1= DDOT(NT1AM(ISYCB),WORK(KETA1),1,WORK(KETA1),1)
          TAL2= DDOT(NT2AM(ISYCB),WORK(KETA2),1,WORK(KETA2),1)
          WRITE(LUPRI,*) 'Printing transformed G^{BC} contribution. 
     &                     Norm2 of singles: ', TAL1,
     &                    'Norm2 of doubles: ', TAL2
        END IF
        LOCDEB=.FALSE.
 
        CALL DAXPY(NT1AM(ISYCB),ONE,WORK(KETA1),1,RHO1,1)
        CALL DAXPY(NT2AM(ISYCB),ONE,WORK(KETA2),1,RHO2,1)
      END IF
C
      CALL QEXIT('CCMM_QRTRANSFORMER')
         
      END
C
C***********************************************************************
C
C  /* Deck ccmm_geff */
      SUBROUTINE CCMM_GEFF(DENS,GEFF,WORK,LWORK)
C
C-----------------------------------------------------------------------------
C
C   Construct effective operators from contracted densities.
C
C   KS, aug 10 
C-----------------------------------------------------------------------------
C
#include "implicit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "qmmm.h"
#include "dummy.h"
#include "iratdef.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
#include "ccfield.h"
#include "exeinf.h"
#include "ccfdgeo.h"
#include "ccslvinf.h"
#include "qm3.h"
#include "ccinftap.h"
#include "nuclei.h"
#include "orgcom.h"

      DIMENSION WORK(LWORK)
      DIMENSION DENS(*), GEFF(*)
      INTEGER NDIM
      LOGICAL LSKIP
      CHARACTER*8 LABEL
C
C Initial parameters
      NDIM=3*NNZAL

C Workspace initialization
      KEPSAO = 1
      KEFIEL = KEPSAO + 3*NNBASX
      KINDIP = KEFIEL + NDIM
      KGAO   = KINDIP + NDIM
      KWRK   = KGAO   + NNBASX
      LWRK   = LWORK - KWRK
      IF (LWRK .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KWRK
         CALL QUIT( 'Too little work in CCMM_GEFF')
      END IF
      
C Convert auxiliary density to electric field
      CALL CCMM_DENS2F(DENS,WORK(KEFIEL),WORK(KWRK),LWRK)

C Transform the field to induced dipoles
      CALL CCMM_F2DIP(WORK(KEFIEL),WORK(KINDIP),WORK(KWRK),LWRK)

      CALL DZERO(WORK(KGAO),NNBASX)

C Construct response G from the induced moments - eg. G^C=-sum_a ∆mu^C eps_a
      LRI = 0

      DO 100 I=1,MMCENT

        LSKIP = .FALSE. 

        CALL CCMM_EPSAO(WORK(KEPSAO),I,LSKIP,WORK(KWRK),LWRK)

        IF (LSKIP) THEN
          IF (ZEROAL(I) .EQ. -1) THEN
            GOTO 100
          ELSE  
            LRI = LRI + 3
            GOTO 100
          END IF
        END IF
C       
        FACx = -WORK(KINDIP + LRI + 0)
 
        FACy = -WORK(KINDIP + LRI + 1)

        FACz = -WORK(KINDIP + LRI + 2)
        
             
        CALL DAXPY(NNBASX,FACx,WORK(KEPSAO),
     *               1,WORK(KGAO),1)

        CALL DAXPY(NNBASX,FACy,WORK(KEPSAO+NNBASX),
     *               1,WORK(KGAO),1)

        CALL DAXPY(NNBASX,FACz,WORK(KEPSAO+2*NNBASX),
     *               1,WORK(KGAO),1)

        LRI = LRI + 3
 100  CONTINUE

C Unpack the integrals and store in GEFF
      LABEL = 'GIVE INT'
      CALL CCMMINT(LABEL,LUQMMM,GEFF,WORK(KGAO),
     &         IRREP,ISYM,IERR,WORK(KWRK),LWRK)

      END
C***********************************************************************
C  /* Deck ccmm_d2ao */
      SUBROUTINE CCMM_D2AO(AODEN,ISYRES,
     *                    LISTB,IDLSTB,ISYMB,
     *                    LISTC,IDLSTC,ISYMC,
     *                    MODEL,WORK,LWORK)
C
C     Purpose: To calculate contributions to the one electron 
C              quadratic auxiliary density matrix to be used
C              in CCMM quadratic response theory.
C              
C
C     Current models: CCSD    
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxash.h"
#include "maxorb.h"
#include "mxcent.h"
#include "aovec.h"
#include "iratdef.h"
#include "ccorb.h"
#include "infind.h"
#include "blocks.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccexci.h"
#include "ccsdio.h"
#include "distcl.h"
#include "cbieri.h"
#include "eribuf.h"
#include "cclr.h"
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION AODEN(*), WORK(LWORK)
      CHARACTER*(10) MODDUM, MODEL
      CHARACTER*(2) LISTB, LISTC
      CHARACTER DENTYP*(1)
      INTEGER ISYRES, IDLSTB, IDLSTC, ISYMB, ISYMC
      LOGICAL LOCDEB
      PARAMETER (LOCDEB=.FALSE.)

      CALL QENTER('CCMM_D2AO')

      ISYML = ILSTSYM('L0',1) ! Symmetry of multipliers 
      ISYMR = ILSTSYM('R0',1) ! Symmetry of amplitudes
C-----------------------------------
C     Initial work space allocation.
C-----------------------------------
C
      KONEIA=1
      KL1AM =KONEIA+NT1AMX
      KL2AM =KL1AM +NT1AMX
      KB1   =KL2AM +NT2SQ(ISYML)
      KB2   =KB1   +NT1AM(ISYMB)
      KC1   =KB2   +NT2AM(ISYMB)
      KC2   =KC1   +NT1AM(ISYMC)
      KWRK1 =KC2   +NT2AM(ISYMC)
      LWRK1 =LWORK -KWRK1
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KWRK1
         CALL QUIT('CCMM_D2AO: Not enough memory')
      ENDIF

C-----------------------------------
C     Reset the working matrices
C-----------------------------------
      CALL DZERO(WORK(KONEIA),NT1AMX)
      CALL DZERO(WORK(KL1AM),NT1AM(ISYML))
      CALL DZERO(WORK(KL2AM),NT2SQ(ISYML))
      CALL DZERO(WORK(KB1),NT1AM(ISYMB))
      CALL DZERO(WORK(KB2),NT2AM(ISYMB))
      CALL DZERO(WORK(KC1),NT1AM(ISYMC))
      CALL DZERO(WORK(KC2),NT2AM(ISYMC))
C-----------------------------------
C     Read in multipliers and singles part of trial vectors
C-----------------------------------
C      WRITE(LUPRI,*) 'NT1AMX: ', NT1AMX
C      WRITE(LUPRI,*) 'NT2AMX: ', NT2AMX
C      WRITE(LUPRI,*) 'NT2AM(ISYML): ', NT2AM(ISYML)
C      WRITE(LUPRI,*) 'NB1: ', NT1AM(ISYMB)
C      WRITE(LUPRI,*) 'NB2: ', NT2AM(ISYMB)
C      WRITE(LUPRI,*) 'NC1: ', NT1AM(ISYMC)
C      WRITE(LUPRI,*) 'NC2: ', NT2AM(ISYMC)
 
      IOPT = 3
      CALL CC_RDRSP('L0',1,ISYML,IOPT,MODDUM,
     *             WORK(KL1AM),WORK(KWRK1))
C      CALL DZERO(WORK(KL1AM),NT1AM(ISYMB))
C      CALL DZERO(WORK(KWRK1),NT2AM(ISYMB))

      CALL CC_RDRSP(LISTB,IDLSTB,ISYMB,IOPT,MODDUM,
     *             WORK(KB1),WORK(KB2))
      CALL CC_RDRSP(LISTC,IDLSTC,ISYMC,IOPT,MODDUM,
     *             WORK(KC1),WORK(KC2))

C scale doubles diagonal
      IF (ISYMB .EQ. 1) CALL CCLR_DIASCL(WORK(KB2),TWO,ISYMB)
      IF (ISYMC .EQ. 1) CALL CCLR_DIASCL(WORK(KC2),TWO,ISYMC)

C Square up the L0_2 block.
      CALL CC_T2SQ(WORK(KWRK1),WORK(KL2AM),ISYML)

      IF (LOCDEB) THEN
        TAL1 = DDOT(NT1AMX,WORK(KL1AM),1,WORK(KL1AM),1)
        TAL2 = DDOT(NT2SQ(ISYML),WORK(KL2AM),1,WORK(KL2AM),1)
        WRITE(LUPRI,*) 'CCMM_D2AO: Norm of L0_1: ', TAL1, 
     &                     'Norm of L0_2: ', TAL2
        TAL1 = DDOT(NT1AMX,WORK(KB1),1,WORK(KB1),1)
        TAL2 = DDOT(NT2AMX,WORK(KB2),1,WORK(KB2),1)
        WRITE(LUPRI,*) 'CCMM_D2AO: Norm of B_1: ', TAL1, 
     &                     'Norm of B_2: ', TAL2
        TAL1 = DDOT(NT1AMX,WORK(KC1),1,WORK(KC1),1)
        TAL2 = DDOT(NT2AMX,WORK(KC2),1,WORK(KC2),1)
        WRITE(LUPRI,*) 'CCMM_D2AO: Norm of C_1: ', TAL1, 
     &                     'Norm of C_2: ', TAL2
      END IF

C----------------------------------------------------------
C     Construct <L2|[[Eia,B1],C2]|HF> contribution to DXia.
C----------------------------------------------------------
      IF (LOCDEB) THEN
        TAL1 = DDOT(NT1AMX,WORK(KONEIA),1,WORK(KONEIA),1)
        WRITE(LUPRI,*) 'CCMM_D2AO: Norm of DIA: ', TAL1 
      END IF
C
      CALL CC_DXIA21(WORK(KONEIA),WORK(KL2AM),ISYML,
     *               WORK(KB1),ISYMB,WORK(KC2),ISYMC,
     *               WORK(KWRK1),LWRK1)
      IF (LOCDEB) THEN
        TAL1 = DDOT(NT1AMX,WORK(KONEIA),1,WORK(KONEIA),1)
        WRITE(LUPRI,*) 'CCMM_D2AO: Norm of DIA after 1: ', TAL1 
      END IF

C
C----------------------------------------------------------
C     Construct <L2|[[Eia,C1],B2]|HF> contribution to DXia.
C----------------------------------------------------------
      CALL CC_DXIA21(WORK(KONEIA),WORK(KL2AM),ISYML,
     *               WORK(KC1),ISYMC,WORK(KB2),ISYMB,
     *               WORK(KWRK1),LWRK1)
      IF (LOCDEB) THEN
        TAL1 = DDOT(NT1AMX,WORK(KONEIA),1,WORK(KONEIA),1)
        WRITE(LUPRI,*) 'CCMM_D2AO: Norm of DIA after 2: ', TAL1 
      END IF

      KONEAI=KWRK1
      KONEAB=KONEAI+NT1AMX
      KONEIJ=KONEAB+NMATAB(1)
      KWRK2=KONEIJ+NMATIJ(1) 
      LWRK2=LWORK-KWRK2
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KWRK2
         CALL QUIT('Insufficient memory for 3. alloc. in CCMM_D1AO')
      ENDIF
      CALL DZERO(WORK(KONEAI),NT1AMX)
      CALL DZERO(WORK(KONEAB),NMATAB(1))
      CALL DZERO(WORK(KONEIJ),NMATIJ(1))
C
C----------------------------------------------------------
C     Construct <L1|[[Eia,B1],C1]|HF> contribution to DXia.
C----------------------------------------------------------
C We do this by doing two contractions of the type 
C sum_ia t_ia b_iq c_pa = sum_i zeta_ip b_iq = D_pq (=D_ia)
C For efficiency we start by performing the sum over virtual orbitals
      CALL CC_DXIJ(WORK(KONEIJ),WORK(KL1AM),ISYML,WORK(KB1),ISYMB)
      CALL CCMM_DIA(WORK(KONEIA),WORK(KONEIJ),ISYRES,WORK(KC1),ISYMC)

      IF (LOCDEB) THEN
        TAL1 = DDOT(NT1AMX,WORK(KONEIA),1,WORK(KONEIA),1)
        WRITE(LUPRI,*) 'CCMM_D2AO: Norm of DIA after 3: ', TAL1 
      END IF
C and the related contribution
      CALL DZERO(WORK(KONEIJ),NMATIJ(1))

      CALL CC_DXIJ(WORK(KONEIJ),WORK(KL1AM),ISYML,WORK(KC1),ISYMC)
      CALL CCMM_DIA(WORK(KONEIA),WORK(KONEIJ),ISYRES,WORK(KB1),ISYMB)

      IF (LOCDEB) THEN
        TAL1 = DDOT(NT1AMX,WORK(KONEIA),1,WORK(KONEIA),1)
        WRITE(LUPRI,*) 'CCMM_D2AO: Norm of DIA after 4: ', TAL1 
      END IF

      KT1AM  = KWRK2
      KLAMDP = KT1AM  + NT1AMX
      KLAMDH = KLAMDP + NLAMDT
      KWRK3  = KLAMDH + NLAMDT
      LWRK3  = LWORK  - KWRK3 
      IF (LWRK3 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KWRK3
         CALL QUIT('Insufficient memory for 2. alloc. in CCMM_D1AO')
      ENDIF

      CALL DZERO(WORK(KT1AM),NT1AMX)
      IOPT = 1

C Read singles to be used for construction of Lambda matrices
      KDUM = KWRK3
      CALL CC_RDRSP('R0',1,ISYMR,IOPT,MODDUM,
     *             WORK(KT1AM),WORK(KDUM))

      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),WORK(KWRK3),
     *            LWRK3)
C
C--------------------------------------------------------
C     Add the one electron density in the AO-basis.
C--------------------------------------------------------
CC
      CALL DZERO(WORK(KONEIJ),NMATIJ(1))
      
      ISDEN = 1
      CALL CC_DENAO(AODEN,ISDEN,WORK(KONEAI),WORK(KONEAB),
     *              WORK(KONEIJ),WORK(KONEIA),ISDEN,WORK(KLAMDP),1,
     *              WORK(KLAMDH),1,WORK(KWRK3),LWRK3)
C
      CALL QEXIT('CCMM_D2AO')

      END

C***********************************************************************
C  /* Deck ccmm_addgdiff */
      SUBROUTINE CCMM_ADDGDIFF(FOCK,WORK,LWORK)
C**************************************************************
C
C     Purpose: For CC2 fullfld Construct the EFFective Diff G solvent operator FROM
C     the diff between cc induced dipoles and the HF INDUCED DIPOLES and add to Fock
C     (in AO basis)
C
C     
C     Get G operators from file dumped 
C     Subtract G^CC-G^HF
C     Add effective operator G to FOCK
C
C     KS 11
C**************************************************************
C
#include "implicit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "qmmm.h"
#include "dummy.h"
#include "iratdef.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
#include "ccfield.h"
#include "exeinf.h"
#include "ccfdgeo.h"
#include "qm3.h"
#include "ccinftap.h"
#include "nuclei.h"
#include "orgcom.h"
C
      PARAMETER (ZERO = 0.0D0, MONE=-1.0D0 ,ONE = 1.0D0,
     &             IZERO = 0 , TWO = 2.0D0)
      PARAMETER (D3I = 1.0D0/3.0D0, D6I = 1.0D0/6.0D0)
      DIMENSION WORK(LWORK),FOCK(*)
      DIMENSION INTREP(9*MXCENT), INTADR(9*MXCENT)
      CHARACTER*8 LABINT(9*MXCENT)
      LOGICAL TOFILE,TRIMAT,EXP1VL,LOCDEB
      PARAMETER (LOCDEB=.FALSE.)
C
      CHARACTER*8 LABEL
C
C
      CALL QENTER('CCMM_ADDGDIFF')
C-------------------------
C       Init parameters
C-------------------------
      NDIM=3*NNZAL
C--------------------------
C       Dynamical allocation
C--------------------------
C
      KCC =  1 
      KHF = KCC + N2BST(ISYMOP)
      KWRK = KHF + N2BST(ISYMOP)
      LWRK1   = LWORK   - KWRK

      IF (LWRK1.LT.0) CALL QUIT( 'Too little work in CCMM_ADDGDIFF' )
        
      CALL DZERO(WORK(KCC),N2BST(ISYMOP))
      CALL DZERO(WORK(KHF),N2BST(ISYMOP))
C

      WRITE(LUPRI,*) 'Reading G from file'
      CALL GET_FROM_FILE('G_MATRIX',N2BST(ISYMOP),WORK(KCC))

      WRITE(LUPRI,*) 'Reading GHF from file'
      CALL GET_FROM_FILE('G_MATRIXHF',N2BST(ISYMOP),WORK(KHF))


      IF (IPQMMM .GT.14 ) THEN
         CALL AROUND( 'CCMM cont. to AO matrix' )
         CALL OUTPUT(WORK(KCC),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      ENDIF
C
      IF (IPQMMM .GT.14) THEN
         CALL AROUND( 'Usual Fock AO matrix' )
         CALL CC_PRFCKAO(FOCK,ISYMOP)
      ENDIF

      CALL DAXPY(N2BST(ISYMOP),-ONE,WORK(KHF),1,WORK(KCC),1)

      TAL=DDOT(N2BST(ISYMOP),WORK(KCC),1,WORK(KCC),1)
      WRITE(LUPRI,*) 'Difference density cont norm ', TAL
C
      CALL DAXPY(N2BST(ISYMOP),ONE,WORK(KCC),1,FOCK,1)

      IF (IPQMMM .GT.14) THEN
         CALL AROUND( 'Modified Fock AO matrix when leaving ADDGDIFF')
         CALL CC_PRFCKAO(FOCK,ISYMOP)
      ENDIF
C
      CALL QEXIT('CCMM_ADDGDIFF')
      END
