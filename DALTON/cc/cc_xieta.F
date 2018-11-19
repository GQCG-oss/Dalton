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
*=====================================================================*
       SUBROUTINE CC_XIETA( IXETRAN, NXETRAN, IOPTRES, IORDER, LISTL, 
     &                      FILXI,   IXDOTS, XCONS, 
     &                      FILETA,  IEDOTS, ECONS, 
     &                      LCAUCHY,  
     &                      MXVEC,   WORK,    LWORK         )
*---------------------------------------------------------------------*
*
*    Purpose: batched loop over Xi and Eta vector calculations
*
*             for more detailed documentation see: CC_XIETA1
*        
*     Written by Christof Haettig, April 1999, based on CC_FMATRIX.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "maxorb.h"
#include "cclists.h"
#include "ccsdinp.h"
#include "ccnoddy.h"
Cholesky
#include "ccdeco.h"
Cholesky

      LOGICAL LOCDBG
      PARAMETER (LOCDBG=.FALSE.)
 
      INTEGER MAXXETRAN, MXVEC
      PARAMETER (MAXXETRAN = 100)

      CHARACTER*8 FN3VI2, FNTOC
      CHARACTER*4 FNDKBC
      CHARACTER*5 FNDKBC3
      CHARACTER*6 FNDELD, FNCKJD, FN3VI, FN3FOPX
      CHARACTER*7 FN3FOP2X

      PARAMETER (FNDELD  = 'CKDELD'  , FNCKJD  = 'CKJDEL'  )
      PARAMETER (FNDKBC  = 'DKBC'    , FNTOC   = 'CCSDT_OC')
      PARAMETER (FN3VI   = 'CC3_VI'  , FN3VI2  = 'CC3_VI12')
      PARAMETER (FN3FOPX = 'PTFOPX'  , FN3FOP2X= 'PTFOP2X' )
      PARAMETER (FNDKBC3 = 'DKBC3'   )


      CHARACTER*(*) LISTL, FILXI, FILETA
      LOGICAL LCAUCHY
      INTEGER IOPTRES, IORDER
      INTEGER NXETRAN, LWORK
      INTEGER IXETRAN(MXDIM_XEVEC,NXETRAN)
      INTEGER IXDOTS(MXVEC,NXETRAN), IEDOTS(MXVEC,NXETRAN)
      
#if defined (SYS_CRAY) 
      REAL WORK(LWORK) 
      REAL XCONS(MXVEC,NXETRAN), ECONS(MXVEC,NXETRAN)
#else 
      DOUBLE PRECISION WORK(LWORK) 
      DOUBLE PRECISION XCONS(MXVEC,NXETRAN), ECONS(MXVEC,NXETRAN)
#endif 

      INTEGER NTRAN, IFIRST, IBATCH, NBATCH, IDUM
      INTEGER ISTART, IEND, IADRX0, IADRX1, IADRE0, IADRE1, IADRPQ0
      INTEGER IADRPQ1, IADRPQMO, IT2DEL0, IT2DELB, IADRH0, IADRH1
      INTEGER KCMOPQ, KCMOHQ, KLAMDPQ, KLAMDHQ, KDNSHFB, KCHI, KCHIQ
      INTEGER KDPKHFB, KDENSB, KDPCKB, KONEHB, KFOCKB, KFCKHFB
      INTEGER KRBIM, KGBIM, KRHO1, KRHO2, KXINT, KYINT, KBFZI, KFBIM
      INTEGER KBFZ0, KF0IM, KRZ0I, KRZBI, KEND,  LEND
      INTEGER KZDPK0, KZDEN0, KZFCK0, KZDPKB, KZDENB, KZFCKB
      INTEGER KIDXL_XIDEN, KIDXL_EADEN, KIDXR_EADEN, MXXIDEN, MXEALEFT
      INTEGER MXEADEN, KIEASTEND, NXIDEN, NEADEN, NEALEFT

      CALL QENTER('CC_XIETA')

*---------------------------------------------------------------------*
* Main section:   CC_XIETA1, driven by a loop over Xi/Eta vectors
* -------------
*   singles and doubles models:
*               calculation of the single and double excitation
*               parts of the Xi and Eta vectors and respective
*               dot products (IOPTRES=5).
*
*   triple models: 
*               calculation of the effective Xi and Eta vectors,
*               for IOPTRES=5 only singles and doubles contributions
*               are computed in CC_XIETA1, the triples contributions
*               to the dot products are computed in CCSDT_XIETA_CONT
*               (see below).
*---------------------------------------------------------------------*


Cholesky
C
      IF (CHOINT .AND. CC2) THEN
         CALL CC_CHOXIETA(IXETRAN,NXETRAN,IOPTRES,IORDER,LISTL,
     &                    FILXI,IXDOTS,XCONS,FILETA,IEDOTS,ECONS,
     &                    MXVEC,WORK,LWORK)
         GOTO 1234
      ENDIF
C
Cholesky
C

      NBATCH = (NXETRAN+MAXXETRAN-1)/MAXXETRAN

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'Batching over Xi/Eta vector calculations:'
        WRITE (LUPRI,*) 'nb. of batches needed:', NBATCH
        CALL FLSHFO(LUPRI)
      END IF
  
      DO IBATCH = 1, NBATCH
        IFIRST = (IBATCH-1) * MAXXETRAN + 1
        NTRAN  = MIN(NXETRAN-(IFIRST-1),MAXXETRAN)

        ISTART      = 1
        IEND        = ISTART   + NTRAN
        IADRX0      = IEND     + NTRAN
        IADRX1      = IADRX0   + MXCORB_CC*1
        IADRE0      = IADRX1   + MXCORB_CC*NTRAN
        IADRE1      = IADRE0   + MXCORB_CC*NTRAN
        IADRPQ0     = IADRE1   + MXCORB_CC*NTRAN
        IADRPQ1     = IADRPQ0  + MXCORB_CC*NTRAN
        IADRPQMO    = IADRPQ1  + MXCORB_CC*NTRAN
        IT2DEL0     = IADRPQMO + MXCORB_CC*NTRAN
        IT2DELB     = IT2DEL0  + MXCORB_CC*1
        IADRH0      = IT2DELB  + MXCORB_CC*NTRAN
        IADRH1      = IADRH0   + MXCORB_CC*NTRAN
        KCMOPQ      = IADRH1   + MXCORB_CC*NTRAN
        KCMOHQ      = KCMOPQ   + NTRAN
        KLAMDPQ     = KCMOHQ   + NTRAN
        KLAMDHQ     = KLAMDPQ  + NTRAN
        KDNSHFB     = KLAMDHQ  + NTRAN
        KDPKHFB     = KDNSHFB  + NTRAN
        KDENSB      = KDPKHFB  + NTRAN
        KDPCKB      = KDENSB   + NTRAN
        KONEHB      = KDPCKB   + NTRAN
        KFOCKB      = KONEHB   + NTRAN
        KFCKHFB     = KFOCKB   + NTRAN
        KRBIM       = KFCKHFB  + NTRAN
        KGBIM       = KRBIM    + NTRAN
        KRHO1       = KGBIM    + NTRAN
        KRHO2       = KRHO1    + NTRAN
        KXINT       = KRHO2    + NTRAN
        KYINT       = KXINT    + NTRAN
        KBFZ0       = KYINT    + NTRAN
        KBFZI       = KBFZ0    + NTRAN
        KF0IM       = KBFZI    + NTRAN
        KFBIM       = KF0IM    + NTRAN
        KRZ0I       = KFBIM    + NTRAN
        KRZBI       = KRZ0I    + NTRAN
        KZDPK0      = KRZBI    + NTRAN
        KZDEN0      = KZDPK0   + NTRAN
        KZFCK0      = KZDEN0   + NTRAN
        KZDPKB      = KZFCK0   + NTRAN
        KZDENB      = KZDPKB   + NTRAN
        KZFCKB      = KZDENB   + NTRAN
        KCHI        = KZFCKB   + NTRAN
        KCHIQ       = KCHI     + NTRAN
        KEND        = KCHIQ    + NTRAN
        LEND        = LWORK    - KEND

        IF (LEND .LT. 0) THEN
           WRITE (LUPRI,*) 'Insufficient work space in CC_XIETA.'
           WRITE (LUPRI,*) 'Available    :',LWORK,' words,'
           WRITE (LUPRI,*) 'Need at least:',KEND, ' words.'
           CALL QUIT('Insufficient work space in CC_XIETA.')
        END IF

        IF (LOCDBG) THEN
          WRITE (LUPRI,*) 'Batch No.:',IBATCH
          WRITE (LUPRI,*) 'start at :',IFIRST
          WRITE (LUPRI,*) '# transf.:',NTRAN
          WRITE (LUPRI,*) 'kend     :',KEND
          IDUM = 0 
          WRITE (LUPRI,*) 'work(0)  :',WORK(IDUM)
        END IF

        CALL CC_XIETA1( IXETRAN(1,IFIRST),NTRAN,IOPTRES,
     &                  IORDER,LISTL,MXVEC,LCAUCHY,
     &                  FILXI, IXDOTS(1,IFIRST),XCONS(1,IFIRST), 
     &                  FILETA,IEDOTS(1,IFIRST),ECONS(1,IFIRST), 
     &                  WORK(ISTART),   WORK(IEND),    WORK(IADRX0),
     &                  WORK(IADRX1),   WORK(IADRE0),  WORK(IADRE1),
     &                  WORK(IADRPQ0),  WORK(IADRPQ1), WORK(IADRPQMO),
     &                  WORK(IT2DEL0),  WORK(IT2DELB),
     &                  WORK(IADRH0),   WORK(IADRH1),
     &                  WORK(KCMOPQ),   WORK(KCMOHQ),
     &                  WORK(KLAMDPQ),  WORK(KLAMDHQ),
     &                  WORK(KDNSHFB),  WORK(KDPKHFB), WORK(KDENSB),
     &                  WORK(KDPCKB),   WORK(KONEHB),  WORK(KFOCKB),
     &                  WORK(KFCKHFB),  WORK(KRBIM),   WORK(KGBIM),
     &                  WORK(KRHO1),WORK(KRHO2),WORK(KXINT),WORK(KYINT),
     &                  WORK(KBFZ0),WORK(KBFZI),WORK(KF0IM),WORK(KFBIM),
     &                  WORK(KRZ0I),WORK(KRZBI),WORK(KCHI),WORK(KCHIQ),
     &                  WORK(KZDPK0),WORK(KZDEN0),WORK(KZFCK0),
     &                  WORK(KZDPKB),WORK(KZDENB),WORK(KZFCKB),
     &                  WORK(KEND), LEND )

        IF (LOCDBG) THEN
          WRITE (LUPRI,*) 'returned from CC_XIETA1:'
          IDUM = 0
          WRITE (LUPRI,*) 'work(0)  :',WORK(IDUM)
        END IF

      END DO

*---------------------------------------------------------------------*
* special triples section:
* ------------------------
*     the calculation of the triples contributions to dot products
*     between a Xi/Eta vector and another response vector are computed
*     via effective one-electron densities, which requires for
*     efficiency a different outermost loop structure, which is setup
*     in CC_XIETA_DENPREP.
*---------------------------------------------------------------------*
      IF (IOPTRES.EQ.5 .AND. CCSDT) THEN

        MXXIDEN  = MXVEC
        MXEALEFT = NXETRAN
        MXEADEN  = MXVEC * MXEALEFT

        KIDXL_XIDEN = 1
        KIDXL_EADEN = KIDXL_XIDEN + MXXIDEN
        KIDXR_EADEN = KIDXL_EADEN + MXEADEN
        KIEASTEND   = KIDXR_EADEN + MXEADEN
        KEND        = KIEASTEND   + 2*MXEALEFT
        LEND        = LWORK       - KEND
        IF (LEND.LT.0) CALL QUIT('Out of memory in CC_XIETA.') 

        CALL CC_XIETA_DENPREP(IXETRAN,NXETRAN,MXVEC,
     &                        IXDOTS,FILXI,
     &                        IEDOTS,FILETA,LISTL,
     &                        WORK(KIDXL_XIDEN),NXIDEN,MXXIDEN,
     &                        WORK(KIDXL_EADEN),WORK(KIDXR_EADEN),
     &                        NEADEN,MXEADEN,
     &                        WORK(KIEASTEND),NEALEFT,MXEALEFT)

        IF (NXIDEN.GT.0) THEN
          IF (NODDY_XI .AND. NODDY_XI_ALTER) THEN

            CALL CCSDT_XI_CONT_NODDY(FILXI,XCONS,
     &                               WORK(KIDXL_XIDEN),NXIDEN,
     &                               IXDOTS,IXETRAN,MXVEC,NXETRAN,
     &                               .FALSE.,'X','X','X','X',
     &                               WORK(KEND),LEND)

          ELSEIF (.NOT.NODDY_XI) THEN
            CALL CC3_XI_CONT(FILXI,XCONS,
     &                           WORK(KIDXL_XIDEN),NXIDEN,
     &                           IXDOTS,IXETRAN,MXVEC,NXETRAN,
     &                           WORK(KEND),LEND)


          END IF
        END IF

        ! call CCSDT_ETA_CONT if NEADEN is nonzero and
        ! none or both of the noddy flags is set:
        !   none of the noddy flag set --> use true code
        !   both of the noddy flag set --> use noddy density
        IF (NEADEN.GT.0 .AND. 
     &      (.NOT.(NODDY_ETA .AND. .NOT.NODDY_ETA_ALTER))) THEN

          IF (LOCDBG) WRITE(LUPRI,*)
     &     'call CCSDT_ETA_CONT... NODD_ETA_ALTER = ',NODDY_ETA_ALTER

          CALL CCSDT_ETA_CONT(LISTL,FILETA,NODDY_ETA_ALTER,
     &                        IEDOTS,ECONS,
     &                        NXETRAN,IXETRAN,MXVEC,
     &                        NEALEFT,WORK(KIEASTEND),
     &                        NEADEN,
     &                        WORK(KIDXL_EADEN),WORK(KIDXR_EADEN),
     &                        WORK(KEND),LEND,
     &                        FNDELD,FNCKJD,FNDKBC,FNTOC,
     &                        FN3VI,FNDKBC3,FN3FOPX,FN3FOP2X)
        END IF
      END IF

 1234 CONTINUE      ! From Cholesky

      CALL QEXIT('CC_XIETA')

      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC_XIETA                             *
*---------------------------------------------------------------------*
c /* deck cc_xieta1 */
*=====================================================================*
       SUBROUTINE CC_XIETA1( IXETRAN, NXETRAN, IOPTRES, 
     &                       IORDER,   LISTL,   MXVEC,  LCAUCHY,
     &                       FILXI,   IXDOTS,  XCONS,  
     &                       FILETA,  IEDOTS,  ECONS, 
     &                       ISTART,  IEND,    IADRX0,
     &                       IADRX1,  IADRE0,  IADRE1,
     &                       IADRPQ0, IADRPQ1, IADRPQMO,
     &                       IT2DEL0, IT2DELB, IADRH0,  IADRH1,
     &                       KCMOPQ,  KCMOHQ,  KLAMDPQ, KLAMDHQ,
     &                       KDNSHFB, KDPKHFB, KDENSB,
     &                       KDPCKB,  KONEHB,  KFOCKB,
     &                       KFCKHFB, KRBIM,   KGBIM,
     &                       KRHO1,   KRHO2,   KXINT,   KYINT,   
     &                       KBFZ0,   KBFZI,   KF0IM,   KFBIM,
     &                       KRZ0I,   KRZBI,   KCHI,    KCHIQ,
     &                       KZDPK0,  KZDEN0,  KZFCK0,  KZDPKB,  
     &                       KZDENB,  KZFCKB,  WORK,    LWORK )
*---------------------------------------------------------------------*
*
*   Purpose: Calculation of the CC vector function or the Jacobian
*            left transformation with a generalized Hamiltonian
*            describing a (first-order) perturbation including the
*            contributions of derivative integrals, orbital relaxation
*            and reorthonormalization (orbital connection).
*
*            Used to calculate right-hand-side vectors for:
*              -- general one-electron perturbations T response "O1"
*              -- geometrical first-derivatives T response "O1"
*              -- magnetic first-derivatives T response "O1"
*              -- general one-electron perturbations T-bar response "X1"
*              -- geometrical first-derivatives T-bar response "X1"
*              -- magnetic first-derivatives T-bar response "X1"
*              -- Cauchy expansion of one-electron T response "CO1"
*
*  <mu|exp(-T^0) Hbar^(1) |CC> and/or <L|exp(-T^0)[Hbar^(1),mu]|CC>
*
*       Hbar^(1) = hbar^(1)_pq E_pq + 1/2 gbar^(1)_pqrs e_pqrs
*
*       hbar^(1)_pq =  h^[1]_pq + 
*                      h^[0]_alp,q LQ^p_alp,p + h^[0]_p,bet LQ^h_bet,q
*       gbar_pqrs defined analogously
*
*    IXETRAN(1,*)  --  operator indices in LBLOPR array 
*    IXETRAN(2,*)  --  indices for left vectors for ETA result vector
*    IXETRAN(3,*)  --  indices for output Xi vector on FILXI list
*    IXETRAN(4,*)  --  indices for output Eta vector on FILETA list
*    IXETRAN(5,*)  --  indices for the 1. orbital relaxation vectors
*                      (for unrelaxed use 0)
*    IXETRAN(6,*)  --  indices for the 2. orbital relaxation vectors
*                      (for unrelaxed use 0, only partially implemented)
*    IXETRAN(7,*)  --  indices for the 3. orbital relaxation vectors
*                      (not yet used....)
*    IXETRAN(8,*)  --  indices for the 4. orbital relaxation vectors
*                      (not yet used....)
*
*    NXETRAN  -- number of requested transformations/vectors
*
*    FILXI    -- list type of the Xi result vectors (IOPTRES=3,4) or
*                of the vectors which are dotted on the Xi vectors
*                (IOPTRES=5)
*    FILETA   -- list type of the Eta result vectors (IOPTRES=3,4) or
*                of the vectors which are dotted on the Eta vectors
*                (IOPTRES=5)
*
*    LCAUCHY  -- compute rhs vectors for the Cauchy expansion of T1
*                (used for CC3 where "effective" rhs vectors needed)
*
*    IXCONS   -- indices of vectors to be dotted on the Xi vectors
*    IECONS   -- indices of vectors to be dotted on the Eta vectors
*
*    XCONS    -- contains for IOPTRES=5 the Xi dot products on return
*    ECONS    -- contains for IOPTRES=5 the Eta dot products on return
*
*    IOPTRES = 0 : all result vectors are written to direct access
*                  files. FILXI and FILETA are used as file names,
*                  the start addresses of the vectors are returned
*                  in IXETRAN(3,*) and IXETRAN(4,*) 
*
*    IOPTRES = 3 : each result vector is written to its own file by
*                  a CALL to CC_WRRSP using FILXI/FILETA as list type
*                  and IXETRAN(3,*)/IXETRAN(4,*) as list index
*
*    IOPTRES = 4 : each result vector is added to a vector on file by
*                  a CALL to CC_WARSP using FILXI/FILETA as list type
*                  and IXETRAN(3,*)/IXETRAN(4,*) as list index
*
*    IOPTRES = 5 : the result vectors are dotted on an array of vectors,
*                  the type of the arrays is given by FILXI for the
*                  Xi result vectors and by FILETA for the Eta result
*                  vectors and the indices are taken from the matrices
*                  IXDOTS and IEDOTS, respectively. the resulting
*                  scalar results are returned in the matrices XCONS
*                  and ECONS.
*                  
*
*
*    IT2DEL0,IT2DELB -- integer scratch arrays of dimensions
*                       MXCORB_CC and MXCORB_CC x NXETRAN
*
*    operator labels:
*           HAM0     : unperturbed Hamiltonian (1-el & 2-el part)
*                      (for test purposes)
*           1DHAMxxx : geometrical first derivatives (1-el & 2-el part)
*
*           ELSE the label is intepreted as a one-electron operator
*                and the reorthonormalization terms are skipped
*                (see subroutine CC_GET_RMAT for details about the
*                 treatment of the connection matrix)
*
*    Written by Christof Haettig, May 1998 -- Jan 1999.
*
*    CC-R12 contributions by Christian Neiss  spring 2005
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "aovec.h"
#include "maxorb.h"
#include "blocks.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "nuclei.h"
#include "symmet.h"
#include "inftap.h"
#include "ccorb.h"
#include "ccfield.h"
#include "ccsdsym.h"
#include "ccsdinp.h"
#include "iratdef.h"
#include "distcl.h"
#include "eritap.h"
#include "ccisao.h"
#include "ccroper.h"
#include "ccfro.h"
#include "cco1rsp.h"
#include "ccx1rsp.h"
#include "ccr1rsp.h"
#include "cclists.h"
#include "chrxyz.h"
#include "dummy.h"
#include "ccsections.h"
#include "qm3.h"
!#include "qmmm.h"
#include "ccnoddy.h"
#include "cch2d.h"
#include "ccrc1rsp.h"
#include "r12int.h"
#include "ccr12int.h"

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      INTEGER ISYM0
      PARAMETER (ISYM0 = 1) ! symmetry of the reference state
      INTEGER LU0IAJB, LU0IABJ, LU0IJBA, LU0AIBJ, LUHINT0
      INTEGER LU1IAJB, LU1IABJ, LU1IJBA, LU1AIBJ, LUHINT1
      INTEGER LUBFDX0, LUBFDX1, LUBFDE0, LUBFDE1, LUPQMO, LUPQ0, LUPQ1
      INTEGER LUDIM, LUCIM, LUBFIM, LUBFZI, LUXI, LUETA,  LUFOCK
      INTEGER LUDKBC, LUCKJD, LUDELD, LUTOC, LU3VI, LU3VI2
      CHARACTER*8 FN0IAJB, FN0IABJ, FN0IJBA, FN0AIBJ, FNHINT0
      CHARACTER*8 FN1IAJB, FN1IABJ, FN1IJBA, FN1AIBJ, FNHINT1
      CHARACTER*8 FNBFDX0, FNBFDX1, FNBFDE0, FNBFDE1, FNBFIM, FNBFZI
      CHARACTER*8 FILPQMO, FILPQ0, FILPQ1, FNDIM, FNCIM, FNTOC
      CHARACTER*8 FN3VI2
      CHARACTER*4 FNDKBC
      CHARACTER*5 FNDKBC3
      CHARACTER*6 FNDELD, FNCKJD, FN3VI, FN3FOPX
      CHARACTER*7 FN3FOP2X
      PARAMETER (FN0IAJB = 'CCXIAJB0', FN0IABJ = 'CCXIABJ0')
      PARAMETER (FN1IAJB = 'CCXIAJB1', FN1IABJ = 'CCXIABJ1')
      PARAMETER (FN0IJBA = 'CCXIJBA0', FN0AIBJ = 'CCXAIBJ0')
      PARAMETER (FN1IJBA = 'CCXIJBA1', FN1AIBJ = 'CCXAIBJ1')
      PARAMETER (FNHINT0 = 'CCHIAJB0', FNHINT1 = 'CCHIAJB1')
      PARAMETER (FNBFDX0 = 'CCBFDEN1', FNBFDX1 = 'CCBFDEN5')
      PARAMETER (FNBFDE0 = 'CCBFDEN2', FNBFDE1 = 'CCBFDEN6')
      PARAMETER (FILPQMO = 'CCPQIMMO', FILPQ0  = 'CCPQINT0')
      PARAMETER (FNCIM   = 'CCXIETAC', FNDIM   = 'CCXIETAD')
      PARAMETER (FNBFIM  = 'CCXEBFIM', FNBFZI  = 'CCXEBFZI')
      PARAMETER (FILPQ1  = 'CCPQINT1')
      PARAMETER (FNDELD  = 'CKDELD'  , FNCKJD  = 'CKJDEL'  )
      PARAMETER (FNDKBC  = 'DKBC'    , FNTOC   = 'CCSDT_OC')
      PARAMETER (FN3VI   = 'CC3_VI'  , FN3VI2  = 'CC3_VI12')
      PARAMETER (FN3FOPX  = 'PTFOPX'   , FN3FOP2X = 'PTFOP2X'  )
      PARAMETER (FNDKBC3 = 'DKBC3'   )


#if defined (SYS_CRAY)
      REAL DDOT, ONE, TWO, THREE, HALF, ZERO
#else
      DOUBLE PRECISION DDOT, ONE, TWO, THREE, HALF, ZERO
#endif
      PARAMETER (ONE = 1.0d0, TWO = 2.0d0, THREE = 3.0d0, HALF = 0.5d0)
      PARAMETER (ZERO = 0.0d0)

      CHARACTER*10 MODEL
      

      CHARACTER*(*) LISTL, FILXI, FILETA
      CHARACTER*8 LABELH, LAB1, LABEL1, LABEL2
      LOGICAL LCAUCHY
      INTEGER LWORK, NXETRAN, IOPTRES, MAXSIM, MXVEC, IORDER
      INTEGER IXETRAN(MXDIM_XEVEC,NXETRAN)
      INTEGER IXDOTS(MXVEC,NXETRAN), IEDOTS(MXVEC,NXETRAN)
      INTEGER IT2DEL0(MXCORB_CC), IT2DELB(MXCORB_CC,NXETRAN)
      INTEGER ISTART(NXETRAN), IEND(NXETRAN)
      INTEGER KRBIM(NXETRAN), KGBIM(NXETRAN), KFBIM(NXETRAN)
      INTEGER KRHO1(NXETRAN), KRHO2(NXETRAN), KBFZI(NXETRAN)
      INTEGER KBFZ0(NXETRAN), KF0IM(NXETRAN), KRZ0I(NXETRAN)
      INTEGER KFOCKB(NXETRAN), KFCKHFB(NXETRAN), KRZBI(NXETRAN)
      INTEGER KONEHB(NXETRAN), KCHI(NXETRAN),   KCHIQ(NXETRAN)
      INTEGER KZDPK0(NXETRAN), KZDEN0(NXETRAN), KZFCK0(NXETRAN)
      INTEGER KZDPKB(NXETRAN), KZDENB(NXETRAN), KZFCKB(NXETRAN)
      INTEGER KLAMDPQ(NXETRAN), KLAMDHQ(NXETRAN)
      INTEGER KDENSB(NXETRAN),  KDPCKB(NXETRAN)
      INTEGER KDNSHFB(NXETRAN), KDPKHFB(NXETRAN)
      INTEGER KXINT(NXETRAN),   KYINT(NXETRAN)
      INTEGER KCMOPQ(NXETRAN),  KCMOHQ(NXETRAN)
      INTEGER IADRX0(MXCORB_CC),          IADRX1(MXCORB_CC,NXETRAN)
      INTEGER IADRH0(MXCORB_CC,NXETRAN),  IADRH1(MXCORB_CC,NXETRAN)

      INTEGER IADRE0(MXCORB_CC,NXETRAN),  IADRE1(MXCORB_CC,NXETRAN)
      INTEGER IADRPQ0(MXCORB_CC,NXETRAN), IADRPQ1(MXCORB_CC,NXETRAN)
      INTEGER IADRPQMO(MXCORB_CC,NXETRAN)

      INTEGER INDEXA(MXCORB_CC), IGAM(MXCORB_CC)

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL AVERAGE, XNORM, FREQ, FREQLST, KTEST
      REAL XCONS(MXVEC,NXETRAN), ECONS(MXVEC,NXETRAN)
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION AVERAGE, XNORM, FREQ, FREQLST, KTEST
      DOUBLE PRECISION XCONS(MXVEC,NXETRAN), ECONS(MXVEC,NXETRAN)
#endif

      CHARACTER*10 MODELW
      LOGICAL LDERINT(8,3), LNEWXI, SKIPXI, SYM1ONLY, SKIPETA
      LOGICAL LTWOEL, LRELAX, LLTWOEL, LLRELAX, NEWTYPE, LZERO
      LOGICAL LZEROIMDONE, LZEROLFT, SQRINT, LPDBS1, LPDBS2, DIRSAV
      INTEGER MT2BGD, MT2ORT, MT2BCD, MT2AM, MT2SQ, M2BST, MEMAT1,MT1AO
      INTEGER MDISAO, MGAMMA, MSCRATCH0, MSCRATCH1, MWORK, NBATCH,MT1AM
      INTEGER ISYM, IATOM, ITRAN, IOPER, IBATCH, IOPT, IFIELD, ILLL
      INTEGER ISYHOP, ISYRES, ISYMD1, IPRERI, IOPER_OLD, IRELAX_OLD
      INTEGER IADRI0, IADRIB, NTOSYM, NTOT, NUMDIS, IDEL2, IDEL, IOPTYP
      INTEGER ISY0DIS, ISY1DIS, ISYDEL, ISYMG, NUMG,ICOOR, ICORSY
      INTEGER NGDER, NBDER, NBUFMX, NFILES, KNRECS
      INTEGER KCCFB1,KINDXB, KEND,LWRK, KODCL1, KODCL2, KODBC1, KODBC2
      INTEGER KRDBC1,KRDBC2, KODPP1,KODPP2, KRDPP1,KRDPP2,KFREE,LFREE
      INTEGER KRECNR, ITYPE, NUMD, MXCOMP, JGAM, IADRF_XI, IADRF_ETA
      INTEGER IFILE, KEND1SV, LWRK1SV, ITST, KT2AM, LEN0, LEN1
      INTEGER IADRBFE0, IADRBFE1, IADRBFX0, IADRBFX1, IADR, IADRCI
      INTEGER IADRPQ, IADRPQI0, IADRPQI1, IADR0, IADRB, IADRDI
      INTEGER IADRZ0, IADRZ1
      INTEGER KONEH0, KFOCK0, KXI1, KXI2, KETA1, KETA2, IOPTW, ISYMJ
      INTEGER IRELAX, IRELAX1, IRELAX2, IRELAX2_OLD, LEN, IOPTWE
      INTEGER KGAMMA, KEMAT1, KEMAT2, ISY0IAJ, ISY1IAJ, KXIAJB, KDPK0
      INTEGER K1XINT, K0XINT, KEND3, LWRK3, IDLSTL, KXAIBJ
      INTEGER KD1PRHF, KD0PRHF, KD1HRHF, KD0HRHF, ICAU
      INTEGER KZETA1, KZETA2, KEND2, LWRK2, ISYCTR, ISYETA
      INTEGER KEND0, LWRK0, KEND1, LWRK1, KG0IM, KR0IM, KFOCKBAO
      INTEGER KMINT, KEND4, LWRK4, NGIM, IDLSTL_OLD, MAXDIF, IVEC
      INTEGER ICORE, KLAMH0, KLAMP0, KDENS0, KT1AMP0, KCMO0, KDNSHF0
      INTEGER KKAPPA, KRMAT, IDXR, IERR, KDPKHF0, JCOOR, JSCOOR, ISCOOR
      INTEGER INUM, ISIGNH, IOPER1, IOPER2, ISYOP1, ISYOP2, IREAL
      INTEGER KRBIMOLD, KFOCK0MO, KXI1EFF, KXI2EFF, KETA1EFF,KETA2EFF
      INTEGER LUDKBC3, LU3FOPX, LU3FOP2X, IOPTWR12, KXIR12, KETAR12
      INTEGER ktr12,ktr12sq,kxmatsq,kxmat,kvxintsq,
     &        kprpao,koffx,kxir12sq,ketar12sq,kctrr12,kctrr12sq
      INTEGER luxint,ian,isymv
      INTEGER IDUM
      INTEGER KSCR1,KSCR2,KSCR3,ISYM1,ISYM2
    

* external functions:
      INTEGER ILSTSYM
      INTEGER IROPER
      INTEGER IROPER2
      INTEGER ILRCAMP

      CALL QENTER('CC_XIETA1')

*---------------------------------------------------------------------*
* begin: check wavefunction model and flush print unit
*---------------------------------------------------------------------*
C     IF (CCSDT) THEN
C       WRITE (LUPRI,'(/1x,a)')
C    &        'CC_XIETA not yet implemented for triples ...'
C       CALL QUIT('Triples not implemented in CC_XIETA.')
C     END IF

      IF ( .NOT. (CCS .OR. CC2 .OR. CCSD .OR. CC3) ) THEN
        WRITE (LUPRI,'(/1x,a)') 'CC_XIETA Called for a Coupled Cluster '
     &          //'method not implemented in CC_XIETA...'
        CALL QUIT('Unk own CC method in CC_XIETA.')
      END IF

      CALL FLSHFO(LUPRI)

      ! save the present value of the DIRECT flag
      DIRSAV = DIRECT

      IF (LOCDBG) THEN
        ITST = 0
        WRITE (LUPRI,'(/1x,a,i15)') 'Work space in CC_XIETA:',LWORK
        WRITE (LUPRI,*) 'FILXI  = ',FILXI(1:3)
        WRITE (LUPRI,*) 'FILETA = ',FILETA(1:3)
        WRITE (LUPRI,*) 'NXETRAN, MXVEC:',NXETRAN,MXVEC
        WRITE (LUPRI,*) 'WORK(0) = ',WORK(ITST)
        CALL FLSHFO(LUPRI)
      END IF

*---------------------------------------------------------------------*
* check return option for the result vectors and initialize output:
*---------------------------------------------------------------------*
      LUXI  = -1
      LUETA = -1
      IF (IOPTRES.EQ.0 .OR. IOPTRES.EQ.1) THEN
         CALL WOPEN2(LUXI, FILXI, 64,0)
         CALL WOPEN2(LUETA,FILETA,64,0)
         IADRF_XI  = 1
         IADRF_ETA = 1
         IF (CCSDT) CALL QUIT('Problem in CC_XIETA.')
      ELSE IF (IOPTRES.EQ.2) THEN
         CALL QUIT('IOPTRES=2 option not implemented in CC_XIETA.')
      ELSE IF (IOPTRES.EQ.3 .OR. IOPTRES.EQ.4) THEN
         CONTINUE
      ELSE IF (IOPTRES.EQ.5) THEN
         IF (MXVEC*NXETRAN .LE. 0) THEN
            WRITE (LUPRI,*)
     &           'WARNING: CC_XIETA called, but nothing to do!'
            RETURN
         END IF
      ELSE
         CALL QUIT('Illegel value of IOPTRES in CC_XIETA.')
      END IF

      IOPTWR12 = 0
      IF (CCS) THEN
         MODELW = 'CCS       '
         IOPTW  = 1
      ELSE IF (CC2) THEN
         MODELW = 'CC2       '
         IF (CCR12) THEN
           MODELW = 'CC2-R12   '
           IOPTWR12 = 32
         END IF
         IOPTW  = 3
      ELSE IF (CCSD) THEN
         MODELW = 'CCSD      '
         IF (CCR12) THEN
           MODELW = 'CCSD-R12  '
           IOPTWR12 = 32
         END IF
         IOPTW  = 3
      ELSE IF (CC3) THEN
         MODELW = 'CC3       '
         IF (CCR12) THEN
           MODELW = 'CC3-R12   '
           IOPTWR12 = 32
CCN        CALL QUIT('No CC3-R12 yet in CC_XIETA!')
         END IF
         IOPTW  = 3
         IOPTWE = 24
      ELSE
         CALL QUIT('Unk own coupled cluster model in CC_XIETA.')
      END IF


*---------------------------------------------------------------------*
* open files for the different integrals:
*---------------------------------------------------------------------*
      LU0AIBJ = -1
      LU0IABJ = -1
      LU0IAJB = -1
      LU0IJBA = -1
      LU1AIBJ = -1
      LU1IABJ = -1
      LU1IAJB = -1
      LU1IJBA = -1
      LUHINT0 = -1
      LUHINT1 = -1
      CALL WOPEN2(LU0IAJB, FN0IAJB, 64, 0)
      CALL WOPEN2(LU1IAJB, FN1IAJB, 64, 0)
      IF (CC2) THEN
        CALL WOPEN2(LU0AIBJ, FN0AIBJ, 64, 0)
        CALL WOPEN2(LU1AIBJ, FN1AIBJ, 64, 0)
        CALL WOPEN2(LUHINT0, FNHINT0, 64, 0)
        CALL WOPEN2(LUHINT1, FNHINT1, 64, 0)
      ELSE
        CALL WOPEN2(LU0IABJ, FN0IABJ, 64, 0)
        CALL WOPEN2(LU1IABJ, FN1IABJ, 64, 0)
        CALL WOPEN2(LU0IJBA, FN0IJBA, 64, 0)
        CALL WOPEN2(LU1IJBA, FN1IJBA, 64, 0)
      END IF

*---------------------------------------------------------------------*
* open files for P and Q intermediates: 
*---------------------------------------------------------------------*
      LUPQMO = -1
      LUPQ0  = -1
      LUPQ1  = -1
      CALL WOPEN2(LUPQMO, FILPQMO, 64, 0)
      CALL WOPEN2(LUPQ0,  FILPQ0,  64, 0)
      CALL WOPEN2(LUPQ1,  FILPQ1,  64, 0)

*=====================================================================*
* prepare batching:
*=====================================================================*
*     estimate memory requirements:
      MT2BGD = 0
      MT2BCD = 0
      MT2ORT = 0
      MT2AM  = 0
      MT2SQ  = 0
      MDISAO = 0
      MEMAT1 = 0
      M2BST  = 0
      MGAMMA = 0
      MT1AO  = 0
      MT1AM  = 0
      DO ISYM = 1, NSYM
        MT2BCD = MAX(MT2BCD,NT2BCD(ISYM))
        MT2BGD = MAX(MT2BGD,NT2BGD(ISYM))
        MT2ORT = MAX(MT2ORT,NT2ORT(ISYM))
        MT2AM  = MAX(MT2AM ,NT2AM(ISYM))
        MT2SQ  = MAX(MT2SQ ,NT2SQ(ISYM))
        MDISAO = MAX(MDISAO,NDISAO(ISYM))
        MEMAT1 = MAX(MEMAT1,NEMAT1(ISYM))
        M2BST  = MAX(M2BST ,N2BST(ISYM) )
        MGAMMA = MAX(MGAMMA,NGAMMA(ISYM))
        MT1AO  = MAX(MT1AO, NT1AO(ISYM))
        MT1AM  = MAX(MT1AM, NT1AM(ISYM))
      END DO

* fixed requirement, independent of batch length:
*     maximum of MT2SQ and MGAMMA (for A, C and D terms)
*     a NT2BGD array (CC_BF) and 6 NT2BCD arrays (transformations)
*     + NT2SQ (T2AMP0, and some other intermediates) + NEMAT1(1) (R0IM)
*     + 2 (**|*del) integral arrays and some extra reserve

      IF (CCS) THEN
        MSCRATCH0 = 2*MDISAO + 10*N2BASX
      ELSE IF (CC2 .OR. CCSD .OR. CCSDT) THEN
        MSCRATCH0 = MAX(MT2SQ,MGAMMA,MT2BGD,6*MT2BCD) + 2*MDISAO + 
     &        MAX(MT2SQ,2*MT2ORT) + NEMAT1(ISYM0) + 10*N2BASX
      ELSE
        CALL QUIT('Unk own CC model in CC_XIETA.')
      END IF

* requirements per simultaneous transformation:
*     BF and R intermediate and Fock matrix and one-electron integrals

      IF (CCS) THEN
        MSCRATCH1 = 2*M2BST
      ELSE IF (CC2 .OR. CCSD .OR. CCSDT) THEN
        MSCRATCH1 = MAX(2*MT2ORT,MT2AM) + MEMAT1 + 2*M2BST
      ELSE
        CALL QUIT('Unk own CC model in CC_XIETA.')
      END IF

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'CC_XIETA> scratch space estimates:'
        WRITE (LUPRI,*) 'CC_XIETA> fixed (MSCRATCH0)   : ',MSCRATCH0
        WRITE (LUPRI,*) 'CC_XIETA> per simult. transf. : ',MSCRATCH1
        WRITE (LUPRI,*) 'WORK(0) = ',WORK(ITST)
      END IF

* prepare batching, make sure that for derivatives all transf. 
* in a batch belong to coordinates of the same atom,
* do not mix derivatives with other perturbations
* -> IOPTYP = 0 : 1el perturbations
*    IOPTYP = 1 : first geom. derivatives
*    IOPTYP = 2 : first magn. derivatives (London integrals)
*
      NBATCH  = 1
      MWORK   = 0
      IATOM   = 0
      IOPTYP  = 0
      ISTART(NBATCH) = 1
      DO ITRAN = 1, NXETRAN

        IOPER  = IXETRAN(1,ITRAN)
        LABELH = LBLOPR(IOPER)

        IF (LABELH(1:5).EQ.'1DHAM') THEN

           IF (DIRECT) THEN
             NEWTYPE = (IOPTYP.NE.1)
           ELSE
             NEWTYPE = (IATOPR(IOPER).NE.IATOM .OR. IOPTYP.NE.1)
           END IF
           IATOM   = IATOPR(IOPER)
           IOPTYP  = 1
           IF ( (NEWTYPE .AND. ITRAN.GT.1) .OR. 
     &          ((MWORK+MSCRATCH1+MSCRATCH0).GT.LWORK) ) THEN
              IEND(NBATCH)   = ITRAN-1
              NBATCH         = NBATCH + 1
              ISTART(NBATCH) = ITRAN
              MWORK = 0
           END IF

        ELSE IF (LABELH(1:5).EQ.'dh/dB') THEN

           NEWTYPE = (IOPTYP.NE.2)
           IATOM   = 0
           IOPTYP  = 2
           IF ( (NEWTYPE .AND. ITRAN.GT.1) .OR. 
     &          ((MWORK+MSCRATCH1+MSCRATCH0).GT.LWORK) ) THEN
              IEND(NBATCH)   = ITRAN-1
              NBATCH         = NBATCH + 1
              ISTART(NBATCH) = ITRAN
              MWORK = 0
           END IF

        ELSE 

           NEWTYPE = (IOPTYP.NE.0)
           IATOM   = 0
           IOPTYP  = 0
           IF ( (NEWTYPE .AND. ITRAN.GT.1) .OR.
     &          ((MWORK+MSCRATCH1+MSCRATCH0).GT.LWORK) ) THEN
              IEND(NBATCH)   = ITRAN-1
              NBATCH         = NBATCH + 1
              ISTART(NBATCH) = ITRAN
              MWORK = 0
           END IF

        END IF

        MWORK = MWORK + MSCRATCH1
        IF ((MWORK+MSCRATCH0) .GT. LWORK) THEN
           WRITE (LUPRI,*) 'Insufficient work space in CC_XIETA.'
           WRITE (LUPRI,*) 'Available:',LWORK,' words.'
           WRITE (LUPRI,*) 'Need at least:',MWORK+MSCRATCH0,' words.'
           CALL FLSHFO(LUPRI)
           CALL QUIT('Insufficient work space in CC_XIETA.')
        END IF

      END DO
      IEND(NBATCH) = NXETRAN

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'CC_XIETA> batching statistics:'
        WRITE (LUPRI,*) 'CC_XIETA> NBATCH : ',NBATCH
        WRITE (LUPRI,*) 'CC_XIETA> ISTART : ',
     &       (ISTART(IBATCH),IBATCH=1,NBATCH)
        WRITE (LUPRI,*) 'CC_XIETA> IEND   : ',
     &       (IEND(IBATCH),  IBATCH=1,NBATCH)
        WRITE (LUPRI,*) 'WORK(0) = ',WORK(ITST)
      END IF

C     -------------------------------------------------------
C     set start address for 0.th-order intermediates on file:
C     -------------------------------------------------------
      IADRI0   = 1
      IADRPQ   = 1
      IADRPQI0 = 1
      IADRBFX0 = 1
      IADRBFE0 = 1

C     --------------------------------------------------------------
C     allocate work space for 0.th-order Lambda matrices, MO vector,
C     densities, Fock matrices, etc., which are kept in memory:
C     --------------------------------------------------------------
      KDENS0  = 1
      KDPK0   = KDENS0  + N2BST(ISYM0)
      KT1AMP0 = KDPK0   + NNBST(ISYM0)
      KLAMP0  = KT1AMP0 + NT1AMX
      KLAMH0  = KLAMP0  + NLAMDT
      KONEH0  = KLAMH0  + NLAMDT
      KFOCK0  = KONEH0  + N2BST(ISYM0)
      KEND0   = KFOCK0  + N2BST(ISYM0)
       
      IF (CC2) THEN
        KDNSHF0 = KEND0
        KDPKHF0 = KDNSHF0 + N2BST(ISYM0)
        KCMO0   = KDPKHF0 + NNBST(ISYM0)
        KEND0   = KCMO0   + NLAMDS
      END IF

      IF (CCSDT) THEN
        KFOCK0MO = KEND0
        KEND0    = KFOCK0MO + N2BST(ISYM0)
      END IF

      LWRK0   = LWORK   - KEND0

      IF (LWRK0 .LT. 0) THEN
        CALL QUIT('Insufficient work space in CC_XIETA. (0)')
      END IF

C     ----------------------------------------------------------------
C     initialize 0.th-order Lambda matrices, densities, etc.:
C     ----------------------------------------------------------------
      IOPT = 1
      CALL CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KT1AMP0),DUMMY)

      CALL LAMMAT(WORK(KLAMP0),WORK(KLAMH0),WORK(KT1AMP0),
     &            WORK(KEND0),LWRK0)

      ICORE = 1
      CALL CC_AODENS(WORK(KLAMP0),WORK(KLAMH0),WORK(KDENS0),
     &               ISYM0,ICORE,WORK(KEND0),LWRK0)
      CALL CC_DNSPK(WORK(KDENS0),WORK(KDPK0),ISYM0)

C     ----------------------------------------------------------------
C     for CC2 we need also the SCF MO coefficients and density matrix:
C     ----------------------------------------------------------------
      IF (CC2) THEN
        CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',
     *              IDUMMY,.FALSE.)
        REWIND LUSIFC
        CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
        READ(LUSIFC)
        READ(LUSIFC)
        READ(LUSIFC) (WORK(KCMO0+I-1),I=1,NLAMDS)
        CALL GPCLOSE(LUSIFC,'KEEP')

        CALL CMO_REORDER(WORK(KCMO0),WORK(KEND0),LWRK0)

        ICORE = 1
        CALL CC_AODENS(WORK(KCMO0),WORK(KCMO0),WORK(KDNSHF0),
     &                 ISYM0,ICORE,WORK(KEND0),LWRK0)
        CALL CC_DNSPK(WORK(KDNSHF0),WORK(KDPKHF0),ISYM0)
      END  IF
C
C     ----------------------------------------------------------------
C     get the zeroth-order one-electron Hamiltonian and Fock matrices:
C     (for CCS we have to calculate the Fock matrix and do only
C      initialize it here with the one-electron part)
C     ----------------------------------------------------------------
*     read zeroth-order one-electron hamiltonian from file 
*     and add finite fields contributions (if any)
      CALL CCRHS_ONEAO(WORK(KONEH0),WORK(KEND0),LWRK0)
      DO IFIELD = 1, NFIELD
        CALL CC_ONEP(WORK(KONEH0),WORK(KEND0),LWRK0,
     &               EFIELD(IFIELD),ISYM0,LFIELD(IFIELD)) 
      END DO
C
C-------------------------------------
C     Solvent contribution.
C     Put into one-electron integrals.
C SLV98,OC
C CCMM02,JA+AO
C-------------------------------------
C
      IF (CCSLV .AND. (.NOT. CCMM )) THEN
         CALL CCSL_RHSTG(WORK(KONEH0),WORK(KEND0),LWRK0)
      ENDIF
      IF (CCMM) THEN
         IF (.NOT. NYQMMM) THEN
            CALL CCMM_RHSTG(WORK(KONEH0),WORK(KEND0),LWRK0)
         ELSE IF (NYQMMM) THEN
            IF (HFFLD) THEN 
              CALL CCMM_ADDGHF(WORK(KONEH0),WORK(KEND0),LWRK0)
            ELSE 
              CALL CCMM_ADDG(WORK(KONEH0),WORK(KEND0),LWRK0)
            END IF
         END IF
      ENDIF
C
C---------------------------------------------
C
      IF (CCS .OR. CC2) THEN
        CALL DCOPY(N2BST(ISYM0),WORK(KONEH0),1,WORK(KFOCK0),1)
      ELSE
*       read zeroth-order AO Fock matrix from file: 
        LUFOCK = -1
        CALL GPOPEN(LUFOCK,'CC_FCKH','OLD',' ','UNFORMATTED',IDUMMY,
     &              .FALSE.)
        REWIND(LUFOCK)
        READ(LUFOCK) (WORK(KFOCK0-1+I),I=1,N2BST(ISYM0))
        CALL GPCLOSE(LUFOCK,'KEEP')
      END IF

C     ------------------------------------------------
C     for CC3 precalculate zeroth-order CC Fock matrix 
C     in Lambda MO basis:
C     ------------------------------------------------
      IF (CCSDT) THEN
        CALL DCOPY(N2BST(1),WORK(KFOCK0),1,WORK(KFOCK0MO),1)
        CALL CC_FCKMO(WORK(KFOCK0MO),WORK(KLAMP0),WORK(KLAMH0),
     &                WORK(KEND0),LWRK0,1,1,1)
      END IF

C     -----------------------------------------------------------------
C     allocate and initialize zero order R & G intermediates with zero:
C     -----------------------------------------------------------------
      IF (.NOT. CCS) THEN
         KR0IM = KEND0
         KG0IM = KR0IM + NEMAT1(ISYM0)
         KEND0 = KG0IM + MT1AO
         LWRK0 = LWORK - KEND0

         IF (LWRK0 .LT. 0) THEN
           CALL QUIT('Insufficient work space in CC_XIETA. (0b)')
         END IF

         CALL DZERO(WORK(KR0IM),NEMAT1(ISYM0))
         CALL DZERO(WORK(KG0IM),NT1AO(ISYM0))
      END IF

*=====================================================================*
* loop over the batches of transformations:
*=====================================================================*
      DO IBATCH = 1, NBATCH

C       -----------------------------------------
C       set start address for intermediate files:
C       -----------------------------------------
        IADRIB   = 1
        IADRPQI1 = 1
        IADRBFX1 = 1
        IADRBFE1 = 1
        
        IADRZ0   = 1
        IADRZ1   = 1

*---------------------------------------------------------------------*
* loop over transformations, allocate work space and 
* initialze R, G, F, BF & BFZ intermediates: 
*---------------------------------------------------------------------*
        LLRELAX     = .FALSE.
        LLTWOEL     = .FALSE.
        LZEROIMDONE = .FALSE.

        KEND1 = KEND0

        IOPER_OLD   = -1
        IRELAX_OLD  = -1
        IDLSTL_OLD  = -1

        IRELAX2     = -1
        IRELAX2_OLD = -1

        DO ITRAN = ISTART(IBATCH), IEND(IBATCH)

          IOPER  = IXETRAN(1,ITRAN)   ! operator index
          IDLSTL = IXETRAN(2,ITRAN)   ! ZETA vector index
          IRELAX = IXETRAN(5,ITRAN)   ! flag for 1. orbital relaxation

          LABELH = LBLOPR(IOPER)      ! operator label
          ISYHOP = ISYOPR(IOPER)      ! symmetry of the operator
          IREAL  = ISYMAT(IOPER)      ! flag for real/imag. operators
          LTWOEL = LPDBSOP(IOPER)     ! two-electron contribution
          LRELAX = LTWOEL .OR. (IRELAX.GE.1)

          SKIPETA = ( IXETRAN(4,ITRAN) .EQ. -1 )

          IF (IORDER.GE.2) THEN
             IRELAX2 = IXETRAN(6,ITRAN) 

             ! get info about the different contributing perturbations
             LABEL1 = '?'
             IF (IRELAX  .GE. 1) LABEL1 = LRTHFLBL(IRELAX)
             LABEL2 = '?'
             IF (IRELAX2 .GE. 1) LABEL2 = LRTHFLBL(IRELAX2)
              
             IF (IORDER.EQ.2) THEN
                INUM   = IROPER2(LABEL1,LABEL2,LABELH,ISIGNH,ISYM)
                IF (LABEL1(1:1).NE.'?') THEN
                   IOPER1 = IROPER(LABEL1,ISYOP1)
                   LPDBS1 = LPDBSOP(IOPER1) 
                END IF
                IF (LABEL2(1:1).NE.'?') THEN
                   IOPER2 = IROPER(LABEL2,ISYOP2)
                   LPDBS2 = LPDBSOP(IOPER2) 
                END IF
                ! for second-order operators we set the LRELAX flag
                ! only if both perturbations are either relaxed or
                ! have field-dep. basis sets, i.e. when we get a
                ! contrib. from the undifferentiated 2e- integrals
                LRELAX = (LPDBS1 .OR. IRELAX .GE.1) .AND.
     &                   (LPDBS2 .OR. IRELAX2.GE.1) 
             ELSE
                CALL QUIT('IORDER out of range in CC_XIETA.')
             END IF

          END IF

          LLTWOEL = (LLTWOEL .OR. LTWOEL)  ! any two-el. perturbation
          LLRELAX = (LLRELAX .OR. LRELAX)  ! any relaxed perturbation

C         LZERO   = (ITRAN.EQ.ISTART(IBATCH)) ! 0.-order intermediates
          IF ( (.NOT.LZEROIMDONE) .AND. LRELAX) THEN
            LZERO = .TRUE.
          ELSE
            LZERO = .FALSE.
          END IF
C
          IF (.NOT.SKIPETA) THEN
            ISYCTR = ILSTSYM(LISTL,IDLSTL) ! sym. of ZETA vector
            ISYETA = MULD2H(ISYHOP,ISYCTR) ! sym. of ETA result vect.
          END IF

C         --------------------------------------------------------------
C         intermediates for Xi & Eta: they depend only on the integrals
C         and the Lambda/LambdaQ matrices and have only to be calculated
C         when IOPER or IRELAX has changed, not if only IDLSTL changed
C         --------------------------------------------------------------
          IF ( (IOPER.NE.IOPER_OLD) .OR. 
     &         (IRELAX.NE.IRELAX_OLD).OR.(IRELAX2.NE.IRELAX2_OLD)) THEN

            IF (LRELAX) THEN
              KCMOPQ(ITRAN)   = KEND1
              KCMOHQ(ITRAN)   = KCMOPQ(ITRAN)   + NGLMDT(ISYHOP)
              KLAMDPQ(ITRAN)  = KCMOHQ(ITRAN)   + NGLMDT(ISYHOP)
              KLAMDHQ(ITRAN)  = KLAMDPQ(ITRAN)  + NGLMDT(ISYHOP)
              KDNSHFB(ITRAN)  = KLAMDHQ(ITRAN)  + NGLMDT(ISYHOP)
              KDPKHFB(ITRAN)  = KDNSHFB(ITRAN)  + N2BST(ISYHOP)
              KDENSB(ITRAN)   = KDPKHFB(ITRAN)  + NNBST(ISYHOP)
              KDPCKB(ITRAN)   = KDENSB(ITRAN)   + N2BST(ISYHOP)
              KEND1           = KDPCKB(ITRAN)   + NNBST(ISYHOP)
            END IF

            IF ( (.NOT. CCS) .AND. (LTWOEL.OR.LRELAX) ) THEN
              NGIM  = MAX(NT1AO(ISYHOP),NT1AO(ISYETA))
              KRBIM(ITRAN) = KEND1
              KGBIM(ITRAN) = KRBIM(ITRAN) + NEMAT1(ISYHOP)
              KEND1        = KGBIM(ITRAN) + NGIM
              IF (.NOT.CC2) THEN
                KRHO2(ITRAN) = KEND1
                KEND1        = KRHO2(ITRAN) + 2 * NT2ORT(ISYHOP)
              END IF
            END IF

          ELSE
            KLAMDPQ(ITRAN)  = KLAMDPQ(ITRAN-1)
            KLAMDHQ(ITRAN)  = KLAMDHQ(ITRAN-1) 
            KDENSB(ITRAN)   = KDENSB(ITRAN-1) 
            KDPCKB(ITRAN)   = KDPCKB(ITRAN-1)
            KDNSHFB(ITRAN)  = KDNSHFB(ITRAN-1)
            KDPKHFB(ITRAN)  = KDPKHFB(ITRAN-1)
            KCMOPQ(ITRAN)   = KCMOPQ(ITRAN-1)
            KCMOHQ(ITRAN)   = KCMOHQ(ITRAN-1)
            KRBIM(ITRAN)    = KRBIM(ITRAN-1) 
            KGBIM(ITRAN)    = KGBIM(ITRAN-1)
            KRHO2(ITRAN)    = KRHO2(ITRAN-1)
          END IF

          KONEHB(ITRAN)  = KEND1
          KFOCKB(ITRAN)  = KONEHB(ITRAN)  + N2BST(ISYHOP)
          KFCKHFB(ITRAN) = KFOCKB(ITRAN)  + N2BST(ISYHOP)
          KEND1          = KFCKHFB(ITRAN) + N2BST(ISYHOP)
          LWRK1          = LWORK - KEND1
          IF (LWRK1 .LT. 0) THEN
             CALL QUIT('Insufficient work space in CC_XIETA. (A2)')
          END IF
C         -----------------------------------------------------------
C         get the orbital relaxation vector and the connection matrix 
C         R, and calculate the derivative Lambda matrices "LambdaQ" 
C         and the derivative density matrices:
C         Note, that we have a orbital relaxation vectors if one of
C         the IRELAX  variables is > 0
C         -----------------------------------------------------------
          IF (LRELAX) THEN

             KKAPPA = KEND1
             KRMAT  = KKAPPA + 2*NALLAI(ISYHOP)
             KEND2  = KRMAT  + N2BST(ISYHOP)
             LWRK2  = LWORK  - KEND2
       
             IF (LWRK2 .LT. 0) THEN
                CALL QUIT('Insufficient work space in CC_XIETA. (2)')
             END IF

             IF (LABELH.EQ.'HAM0    ' .AND. IRELAX.GE.1) THEN
               IRELAX = 0         
               WRITE (LUPRI,*) 'Test case "HAM0"... no '//
     &              'relaxation vector used.'
             END IF

             IF (IORDER.EQ.1 .AND. IRELAX.GE.1) THEN
               CALL CC_RDHFRSP('R1 ',IRELAX,ISYHOP,WORK(KKAPPA))
             ELSE IF (IORDER.GT.1) THEN
               CALL QUIT(
     &            'IORDER > 1 not fully implemented in CC_XIETA.')
             ELSE 
               CALL DZERO(WORK(KKAPPA),2*NALLAI(ISYHOP))
             END IF
 
             IF (LOCDBG) THEN
               WRITE (LUPRI,*) 'orbital relaxation vector:',
     &                         IRELAX,'  ',LABELH
               CALL OUTPUT(WORK(KKAPPA),1,2*NALLAI(ISYHOP),1,1,
     &                                    2*NALLAI(ISYHOP),1,1,LUPRI)
             END IF

             CALL CC_GET_RMAT(WORK(KRMAT),IOPER,IORDER,ISYHOP,
     &                        WORK(KEND2),LWRK2)

             IOPT = 1
             CALL CC_LAMBDAQ(WORK(KLAMDPQ(ITRAN)),WORK(KLAMDHQ(ITRAN)),
     &                       WORK(KCMOPQ(ITRAN)), WORK(KCMOHQ(ITRAN)),
     &                       ISYHOP,WORK(KT1AMP0),
     &                       WORK(KKAPPA),WORK(KRMAT),IREAL,IOPT,
     &                       WORK(KEND2),LWRK2)

C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C         THIS is still unfinished...
C         frozen core contributions have to be introduced!
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ICORE = 1
            IOPT  = 2
            CALL CC_AODENS2(WORK(KLAMP0),WORK(KLAMH0),ISYM0,
     &               WORK(KLAMDPQ(ITRAN)),WORK(KLAMDHQ(ITRAN)),ISYHOP,
     &               WORK(KDENSB(ITRAN)),ICORE,IOPT,WORK(KEND2),LWRK2)

            IF (CC2) THEN
             ICORE = 1
             IOPT  = 2
             CALL CC_AODENS2(WORK(KCMO0),WORK(KCMO0),ISYM0,
     &                WORK(KCMOPQ(ITRAN)),WORK(KCMOHQ(ITRAN)),ISYHOP,
     &                WORK(KDNSHFB(ITRAN)),ICORE,IOPT,WORK(KEND2),LWRK2)
            END IF

          END IF

C         ------------------------------------------------------------
C         initialize the remaining intermediates with zero:
C         ------------------------------------------------------------
          IF ( (.NOT. CCS) .AND. (LTWOEL.OR.LRELAX) ) THEN
             CALL DZERO(WORK(KGBIM(ITRAN)),NT1AO(ISYHOP))
          END IF
          IF ( CCSD .AND. (LTWOEL.OR.LRELAX) ) THEN
             CALL DZERO(WORK(KRHO2(ITRAN)),2*NT2ORT(ISYHOP))
          END IF


C         ----------------------------
C         extra intermediates for Eta:
C         ----------------------------
          IF (.NOT.CCS) THEN
            ! X & Y intermediates are only calculate if IDLSTL changed
            ! else we reuse the previous one
            IF ((IDLSTL .NE. IDLSTL_OLD).AND.(.NOT.SKIPETA)) THEN
              KXINT(ITRAN) = KEND1
              KYINT(ITRAN) = KXINT(ITRAN) + NMATIJ(ISYCTR)
              KEND1        = KYINT(ITRAN) + NMATAB(ISYCTR)
            ELSE
              KXINT(ITRAN) = KXINT(ITRAN-1)
              KYINT(ITRAN) = KYINT(ITRAN-1) 
            END IF

            IF (LTWOEL.OR.LRELAX) THEN
               KFBIM(ITRAN) = KEND1
               KEND1        = KFBIM(ITRAN) + NT1AO(ISYETA)
               IF (CC2) THEN
                 KZDPKB(ITRAN) = KEND1
                 KZDENB(ITRAN) = KZDPKB(ITRAN) + NNBST(ISYETA)
                 KZFCKB(ITRAN) = KZDENB(ITRAN) + N2BST(ISYETA)
                 KCHIQ(ITRAN)  = KZFCKB(ITRAN) + N2BST(ISYETA)
                 KEND1         = KCHIQ(ITRAN)  + NGLMDT(ISYETA)
               ELSE
                 KBFZI(ITRAN) = KEND1
                 KEND1        = KBFZI(ITRAN) + 2 * NT2ORT(ISYETA)
               END IF

C              IF (IDLSTL .NE. IDLSTL_OLD) THEN
                 KF0IM(ITRAN) = KEND1
                 KEND1        = KF0IM(ITRAN) + NT1AO(ISYCTR)
                 IF (CC2) THEN
                   KZDPK0(ITRAN) = KEND1
                   KZDEN0(ITRAN) = KZDPK0(ITRAN) + NNBST(ISYCTR)
                   KZFCK0(ITRAN) = KZDEN0(ITRAN) + N2BST(ISYCTR)
                   KCHI(ITRAN)   = KZFCK0(ITRAN) + N2BST(ISYCTR)
                   KEND1         = KCHI(ITRAN)   + NGLMDT(ISYCTR)
                 ELSE
                   KBFZ0(ITRAN) = KEND1
                   KEND1        = KBFZ0(ITRAN) + 2 * NT2ORT(ISYCTR)
                 END IF
C              ELSE
C                KF0IM(ITRAN) = KF0IM(ITRAN-1)
C                IF (.NOT.CC2) KBFZ0(ITRAN) = KBFZ0(ITRAN-1) 
C              END IF

            END IF
          END IF

          LWRK1   = LWORK - KEND1
          IF (LWRK1 .LT. 0) THEN
             CALL QUIT('Insufficient work space in CC_XIETA. (B)')
          END IF


          IF ((.NOT.CCS) .AND. (LTWOEL.OR.LRELAX) ) THEN
             CALL DZERO(WORK(KFBIM(ITRAN)),NT1AO(ISYETA))
             IF (CC2) THEN
                 CALL DZERO(WORK(KZFCKB(ITRAN)),N2BST(ISYETA))
             ELSE
                 CALL DZERO(WORK(KBFZI(ITRAN)),2*NT2ORT(ISYETA))
             END IF
C            IF (IDLSTL .NE. IDLSTL_OLD) THEN
               CALL DZERO(WORK(KF0IM(ITRAN)),NT1AO(ISYCTR))
               IF (CC2) THEN
                 CALL DZERO(WORK(KZFCK0(ITRAN)),N2BST(ISYCTR))
               ELSE
                 CALL DZERO(WORK(KBFZ0(ITRAN)),2*NT2ORT(ISYCTR))
               END IF
C            END IF
          END IF

          IOPER_OLD  = IOPER
          IRELAX_OLD = IRELAX
          IF (SKIPETA) THEN
            IDLSTL_OLD = -1
          ELSE
            IDLSTL_OLD = IDLSTL
          END IF
          

          IF (LZERO) LZEROIMDONE = .TRUE.

          IF (LOCDBG) THEN
            WRITE (LUPRI,*) 'CC_XIETA> memory alloc. for AO section:'
            WRITE (LUPRI,*) 'CC_XIETA> LABELH,ISYHOP = ',LABELH,ISYHOP
            WRITE (LUPRI,*) 'CC_XIETA> IOPER,IRELAX,IDLSTL=',
     *                       IOPER,IRELAX,IDLSTL
            WRITE (LUPRI,*) 'CC_XIETA> LTWOEL,LRELAX = ',LTWOEL,LRELAX
            WRITE (LUPRI,*) 'CC_XIETA> KEND1, LWRK1:',KEND1,LWRK1
          END IF

        END DO ! ITRAN

*---------------------------------------------------------------------*
* open files for several effective densities: 
*---------------------------------------------------------------------*
        LUBFDX0 = -1
        LUBFDX1 = -1
        LUBFDE0 = -1
        LUBFDE1 = -1
        CALL WOPEN2(LUBFDX0, FNBFDX0, 64, 0)
        CALL WOPEN2(LUBFDX1, FNBFDX1, 64, 0)
        CALL WOPEN2(LUBFDE0, FNBFDE0, 64, 0)
        CALL WOPEN2(LUBFDE1, FNBFDE1, 64, 0)
     
*=====================================================================*
* Precalculate X, Y, P, and Q intermediates and the effective density
* for the BFZ intermediate for the Eta vector calculations.
* For these intermediates we need the amplituded packed in core
* and the Lagrangian multipliers squared. 
* (The Zeta vector is read in CCETAINT1 into WORK(KZETA2) )
*=====================================================================*
      IF (CCS) THEN
         KEND2 = KEND1
         LWRK2 = LWRK1
         KT2AM = IDUMMY
      ELSE
         KT2AM = KEND1
         KEND2 = KT2AM  + NT2AM(ISYM0)
         LWRK2 = LWORK  - KEND2

         IF (LWRK2 .LT. 0) THEN
           CALL QUIT('Insufficient work space in CC_XIETA. (CCETAINT1)')
         END IF

         IOPT = 2
         CALL CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,DUMMY,WORK(KT2AM))
      END IF

      DO ITRAN = ISTART(IBATCH), IEND(IBATCH)

         IOPER  = IXETRAN(1,ITRAN) ! operator index
         IDLSTL = IXETRAN(2,ITRAN) ! ZETA vector index
         IRELAX = IXETRAN(5,ITRAN) ! flag for orbital relaxation
         LABELH = LBLOPR(IOPER)    ! operator label
         ISYHOP = ISYOPR(IOPER)    ! symmetry of hamiltonian
         LTWOEL = LPDBSOP(IOPER)   ! two-electron contribution
         LRELAX = LTWOEL .OR. (IRELAX.GE.1)

         SKIPETA = ( IXETRAN(4,ITRAN) .EQ. -1 )

         IF (.NOT. SKIPETA) THEN
           ISYCTR = ILSTSYM(LISTL,IDLSTL)  ! symmetry of ZETA vector
         END IF


         IF (IORDER.GE.2) THEN
            IRELAX2 = IXETRAN(6,ITRAN) 

            ! get info about the different contributing perturbations
            LABEL1 = '?'
            IF (IRELAX  .GE. 1) LABEL1 = LRTHFLBL(IRELAX)
            LABEL2 = '?'
            IF (IRELAX2 .GE. 1) LABEL2 = LRTHFLBL(IRELAX2)
             
            IF (IORDER.EQ.2) THEN
               INUM   = IROPER2(LABEL1,LABEL2,LABELH,ISIGNH,ISYM)
               IF (LABEL1(1:1).NE.'?') THEN
                  IOPER1 = IROPER(LABEL1,ISYOP1)
                  LPDBS1 = LPDBSOP(IOPER1) 
               END IF
               IF (LABEL2(1:1).NE.'?') THEN
                  IOPER2 = IROPER(LABEL2,ISYOP2)
                  LPDBS2 = LPDBSOP(IOPER2) 
               END IF
               ! for second-order operators we set the LRELAX flag
               ! only if both perturbations are either relaxed or
               ! have field-dep. basis sets, i.e. when we get a
               ! contrib. from the undifferentiated 2e- integrals
               LRELAX = (LPDBS1 .OR. IRELAX .GE.1) .AND.
     &                  (LPDBS2 .OR. IRELAX2.GE.1) 
            ELSE
               CALL QUIT('IORDER out of range in CC_XIETA.')
            END IF

         END IF


         IF (.NOT. SKIPETA) THEN
            KZETA1 = KEND2
            KZETA2 = KZETA1 + NT1AM(ISYCTR)
            KEND3  = KZETA2 + NT2SQ(ISYCTR)
            IF (.NOT.(CCS.OR.CC2)) THEN
              KMINT  = KEND3 
              KEND3  = KMINT  + N3ORHF(ISYCTR)
            END IF
            LWRK3  = LWORK  - KEND3

            CALL CCETAINT1(ITRAN,ISTART(IBATCH),LISTL,IDLSTL,
     &           WORK(KZDPK0(ITRAN)),WORK(KZDEN0(ITRAN)),
     &           WORK(KZDPKB(ITRAN)),WORK(KZDENB(ITRAN)),
     &           WORK(KXINT(ITRAN)),WORK(KYINT(ITRAN)),WORK(KMINT),
     &           WORK(KCHI(ITRAN)),WORK(KCHIQ(ITRAN)),
     &           WORK(KZETA1),WORK(KZETA2),ISYCTR, 
     &           WORK(KT2AM),WORK(KLAMP0), WORK(KLAMH0), ISYM0,
     &           WORK(KLAMDPQ(ITRAN)),WORK(KLAMDHQ(ITRAN)), ISYHOP,
     &           FNBFDE0, LUBFDE0, IADRE0,  IADRBFE0,
     &           FNBFDE1, LUBFDE1, IADRE1,  IADRBFE1,
     &           FILPQMO, LUPQMO,  IADRPQMO,IADRPQ,
     &           FILPQ0,  LUPQ0,   IADRPQ0, IADRPQI0,
     &           FILPQ1,  LUPQ1,   IADRPQ1, IADRPQI1,
     &           LRELAX,  LTWOEL,  WORK(KEND3),LWRK3)

         ELSE IF (ITRAN.GT.1) THEN
           DO IDEL = 1, NBAST
             IADRE0(IDEL,ITRAN)   = IADRE0(IDEL,ITRAN-1)
             IADRPQMO(IDEL,ITRAN) = IADRPQMO(IDEL,ITRAN-1)
             IADRPQ0(IDEL,ITRAN)  = IADRPQ0(IDEL,ITRAN-1)
             IADRPQ1(IDEL,ITRAN)  = -1
           END DO
         END IF

      END DO

*---------------------------------------------------------------------*
* for the following we need the zeroth-order cluster amplitudes
* squared in memory:
*---------------------------------------------------------------------*
      IF (CCS) THEN
         KEND2 = KEND1
         LWRK2 = LWRK1
         KT2AM = IDUMMY
      ELSE
         KT2AM = KEND1
         KEND2 = KT2AM + NT2SQ(ISYM0)
         LWRK2 = LWORK - KEND2

         IF (LWRK2 .LT. NT2AM(ISYM0)) THEN
           CALL QUIT('Insufficient work space in CC_XIETA. (CCXIINT1)')
         END IF

         IOPT = 2
         CALL CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,DUMMY,WORK(KEND2))
         CALL CC_T2SQ(WORK(KEND2),WORK(KT2AM),ISYM0)
      END IF

*---------------------------------------------------------------------*
* initialize some zero-order intermediates
*---------------------------------------------------------------------*
      IF (LLTWOEL .OR. LLRELAX) THEN
        IF (CCS) THEN
          CALL DCOPY(N2BST(ISYM0),WORK(KONEH0),1,WORK(KFOCK0),1)
        ELSE IF (CC2) THEN 
          CALL DCOPY(N2BST(ISYM0),WORK(KONEH0),1,WORK(KFOCK0),1)
          CALL DZERO(WORK(KG0IM),NT1AO(ISYM0))
        END IF
      END IF

*=====================================================================*
* Precalculate some intermediates for Xi & Eta vector calculation,
* including the effective densities for the BF intermediates:
* - some zeroth-order intermediates need only be computed once
* - the others have to be computed whenever IOPER or IRELAX changes
* - for these intermediates we need the amplitudes squared, but
*   no Lagrangian multiplier vectors.
*=====================================================================*

      IOPER_OLD   = -1
      IRELAX_OLD  = -1
      LZEROIMDONE = .FALSE.

      DO ITRAN = ISTART(IBATCH), IEND(IBATCH)

         IOPER  = IXETRAN(1,ITRAN) ! operator index
         IRELAX = IXETRAN(5,ITRAN) ! flag for orbital relaxation
         ISYHOP = ISYOPR(IOPER)    ! symmetry of hamiltonian
         LABELH = LBLOPR(IOPER)    ! label for 1-el integrals
         LTWOEL = LPDBSOP(IOPER)   ! two-electron contribution
         LRELAX = LTWOEL .OR. (IRELAX.GE.1)

C        LZERO  = (ITRAN.EQ.ISTART(IBATCH)) ! 0.-order intermediates
         IF ( (.NOT.LZEROIMDONE) .AND. LRELAX) THEN
           LZERO = .TRUE.
         ELSE
           LZERO = .FALSE.
         END IF

         LNEWXI=(IOPER.NE.IOPER_OLD).OR.(IRELAX.NE.IRELAX_OLD).OR.LZERO

         IF (IORDER.GE.2) THEN
            IRELAX2 = IXETRAN(6,ITRAN) 

            ! get info about the different contributing perturbations
            LABEL1 = '?'
            IF (IRELAX  .GE. 1) LABEL1 = LRTHFLBL(IRELAX)
            LABEL2 = '?'
            IF (IRELAX2 .GE. 1) LABEL2 = LRTHFLBL(IRELAX2)
             
            IF (IORDER.EQ.2) THEN
               INUM   = IROPER2(LABEL1,LABEL2,LABELH,ISIGNH,ISYM)
               IF (LABEL1(1:1).NE.'?') THEN
                  IOPER1 = IROPER(LABEL1,ISYOP1)
                  LPDBS1 = LPDBSOP(IOPER1) 
               END IF
               IF (LABEL2(1:1).NE.'?') THEN
                  IOPER2 = IROPER(LABEL2,ISYOP2)
                  LPDBS2 = LPDBSOP(IOPER2) 
               END IF
               ! for second-order operators we set the LRELAX flag
               ! only if both perturbations are either relaxed or
               ! have field-dep. basis sets, i.e. when we get a
               ! contrib. from the undifferentiated 2e- integrals
               LRELAX = (LPDBS1 .OR. IRELAX .GE.1) .AND.
     &                  (LPDBS2 .OR. IRELAX2.GE.1) 
            ELSE
               CALL QUIT('IORDER out of range in CC_XIETA.')
            END IF

         END IF

         CALL CCXIINT1(ITRAN, LABELH, IORDER, WORK(KT2AM),
     &                 WORK(KDPCKB(ITRAN)),WORK(KDPKHFB(ITRAN)),
     &                 WORK(KFOCKB(ITRAN)),WORK(KFCKHFB(ITRAN)),
     &                 WORK(KLAMP0), WORK(KLAMH0), ISYM0,
     &                 WORK(KLAMDPQ(ITRAN)),WORK(KLAMDHQ(ITRAN)), 
     &                 WORK(KDENSB(ITRAN)), WORK(KDNSHFB(ITRAN)),
     &                 WORK(KONEHB(ITRAN)), ISYHOP,
     &                 FNBFDX0, LUBFDX0, IADRX0,  IADRBFX0,
     &                 FNBFDX1, LUBFDX1, IADRX1,  IADRBFX1,
     &                 LRELAX,  LTWOEL,  LZERO,   LNEWXI,
     &                 WORK(KEND2),LWRK2)

         IOPER_OLD  = IOPER
         IRELAX_OLD = IRELAX

         IF (LZERO) LZEROIMDONE = .TRUE.

      END DO

*---------------------------------------------------------------------*
* initialize integral program:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
         WRITE (LUPRI,*) 'Initialize now integral program:'
         WRITE (LUPRI,*) 'LLTWOEL = ',LLTWOEL
         WRITE (LUPRI,*) 'LABELH  = ',LABELH
      END IF

      IF (LLTWOEL .AND. 
     &    (LABELH(1:5).EQ.'1DHAM' .OR. LABELH(1:5).EQ.'dh/dB') ) THEN
         
         ! for derivative integrals use the DIRGRD flag to
         ! switch between direct/non-direct mode
         DIRECT = DIRGRD 

         SYM1ONLY = .FALSE.
         CALL CC_SETDORPS(LABELH,SYM1ONLY,0)

         IOPER   = IXETRAN(1,ISTART(IBATCH)) 
         IATOM   = IATOPR(IOPER)

         IF (.NOT.DIRECT) THEN
           CALL CCDER1(IATOM,LABELH,LDERINT,WORK(KEND2),LWRK2)
         END IF

      ELSE
         DIRECT = DIRSAV
      END IF

      IF (LLTWOEL .OR. LLRELAX .OR. LABELH(1:5).EQ.'HAM0 ') THEN

         IF (DIRECT) THEN
            NTOSYM = 1
C           IF (HERDIR) THEN
C             CALL HERDI1(WORK(KEND),LWRK,IPRERI)
C           ELSE
              KCCFB1 = KEND2
              KINDXB = KCCFB1 + MXPRIM*MXCONT
              KEND   = KINDXB + (8*MXSHEL*MXCONT + 1)/IRAT
              LWRK   = LWORK  - KEND
              CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     *                    KODPP1,KODPP2,KRDPP1,KRDPP2,
     *                    KFREE,LFREE,KEND,WORK(KCCFB1),WORK(KINDXB),
     *                    WORK(KEND),LWRK,IPRERI)
 
C           END IF
            KEND2 = KFREE
            LWRK2 = LFREE
         ELSE
            NTOSYM = NSYM
         END IF

      END IF

*---------------------------------------------------------------------*
* loop over integral distributions:
*---------------------------------------------------------------------*
      IF (LLTWOEL .OR. LLRELAX) THEN

        DO ISYMD1 = 1, NTOSYM

          IF (DIRECT) THEN
C           IF (HERDIR) THEN
C              NTOT = MAXSHL
C           ELSE
               NTOT = MXCALL
C           ENDIF                 
          ELSE
            NTOT = NBAS(ISYMD1)
          END IF

          DO ILLL = 1, NTOT
     
            IF (DIRECT) THEN

              IF      (LABELH(1:5).EQ.'HAM0 ' .OR. (.NOT.LLTWOEL)) THEN
                NGDER = 0
                NBDER = 0
                NFILES = 1
              ELSE IF (LABELH(1:5).EQ.'1DHAM') THEN
                NGDER = 1
                NBDER = 0
                NFILES = 1 + 3*NUCDEP
              ELSE IF (LABELH(1:5).EQ.'dh/dB') THEN
                NGDER = 0
                NBDER = 1
                NFILES = 1 + 3*2
              ELSE
                CALL QUIT('Unknown 2e- integral type in CC_XIETA.')
              END IF

              KEND = KEND2
              LWRK = LWRK2

C             IF (HERDIR) THEN
C                CALL HERDI2(WORK(KEND),LWRK,INDEXA,ILLL,NUMDIS,
C    &                       IPRINT)
C             ELSE
                 CALL ERIDI2(ILLL,INDEXA,NUMDIS,NGDER,NBDER,
     *                       WORK(KODCL1),WORK(KODCL2),
     *                       WORK(KODBC1),WORK(KODBC2),
     *                       WORK(KRDBC1),WORK(KRDBC2),
     *                       WORK(KODPP1),WORK(KODPP2),
     *                       WORK(KRDPP1),WORK(KRDPP2),
     *                       WORK(KCCFB1),WORK(KINDXB),
     *                       WORK(KEND), LWRK,IPRERI)
C             END IF

              NBUFMX = NBUFX(0)
              DO I = 1, NFILES-1
                NBUFMX = MAX(NBUFMX,NBUFX(I))
              END DO

              ! allocate memory for orbital/record labels
              KRECNR = KEND
              KNRECS = KRECNR + (NBUFMX*NFILES - 1)/IRAT + 1
              KEND3  = KNRECS + (NFILES - 1)/IRAT + 1
              LWRK3  = LWORK - KEND3

              IF (LWRK3 .LT. 0) THEN
               CALL QUIT('Insufficient work space in CC_XIETA (ERIDI2)')
              END IF
 
              CALL RDERILBS(WORK(KRECNR),WORK(KNRECS),NBUFMX,NFILES)
  
            ELSE
              NUMDIS = 1
              KEND3  = KEND2
              LWRK3  = LWRK2
            END IF

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'CC_XIETA> after integral evaluation '//
     &        'for ILLL:',ILLL
        WRITE (LUPRI,*) 'CC_XIETA> KEND3, LWRK3:',KEND3,LWRK3
        WRITE (LUPRI,*) 'WORK(0) = ',WORK(ITST)
      END IF

*=====================================================================*
* Calculate intermediates depending on the two-electron AO integrals,
* in the loop over AO integral shells. We put the loop over the result
* vectors here, because the required derivatives integrals can change
* for each ITRAN. We assume that all Xi & Eta vectors for operators
* coming from the same IATOM are grouped together. In this way we
* avoid recalculation of the derivative integrals for each vector, but
* accept that the AO integrals are reread from file for each ITRAN.
*=====================================================================*
      IOPER_OLD   = -1
      IDLSTL_OLD  = -1
      IRELAX_OLD  = -1
      LZEROIMDONE = .FALSE.

      DO ITRAN = ISTART(IBATCH), IEND(IBATCH)

        IOPER  = IXETRAN(1,ITRAN)   ! operator index
        IDLSTL = IXETRAN(2,ITRAN)   ! ZETA vector index
        IRELAX = IXETRAN(5,ITRAN)   ! flag for orbital relaxation

        ISYHOP = ISYOPR(IOPER)      ! symmetry of hamiltonian
        IREAL  = ISYMAT(IOPER)      ! flag for real/imag. operators
        IATOM  = IATOPR(IOPER)      ! associated atom
        LABELH = LBLOPR(IOPER)      ! operator label
        LTWOEL = LPDBSOP(IOPER)     ! pert.-dep. basis set
        LRELAX = LTWOEL .OR. (IRELAX.GE.1)

        SKIPETA = ( IXETRAN(4,ITRAN) .EQ. -1 )

        ! number of components (directions) in total for this operator,
        ! coordinate index (ICOOR) and coordinate symmetry (ICORSY)
        MXCOMP = 0
        IF      ( LABELH(1:5).EQ.'HAM0 ' ) THEN
           MXCOMP = 1
           ICOOR  = 1
           ICORSY = 1
           ISCOOR = 0
        ELSE IF ( LABELH(1:5).EQ.'dh/dB' ) THEN
           MXCOMP = 3
           ICORSY = ISYHOP
           DO JCOOR = 1, MXCOMP
              IF (CHRXYZ(JCOOR).EQ.LABELH(6:6)) THEN
                 ICOOR = JCOOR
              END IF
           END DO
           ISCOOR = ICOOR
        ELSE IF ( LABELH(1:5).EQ.'1DHAM' ) THEN
           MXCOMP = 3
           ICORSY = ISYHOP
           READ(LABELH,'(A5,I3)') LAB1,ISCOOR

           DO JCOOR = 1, MXCOMP
              JSCOOR = IPTCNT(3*(IATOM-1)+JCOOR,ICORSY-1,1)
              IF (JSCOOR.EQ.ISCOOR) THEN
                 ICOOR = JCOOR
              END IF
           END DO
        ELSE IF (LTWOEL) THEN
           WRITE (LUPRI,*) 'Unk own 2el operator label in CC_XIETA.'
           CALL QUIT('Unk own 2el operator label in CC_XIETA.')
        END IF

        IF (LOCDBG) THEN
          WRITE (LUPRI,*) 'CC_XIETA> IOPER         = ',IOPER
          WRITE (LUPRI,*) 'CC_XIETA> LTWOEL,LRELAX = ',LTWOEL,LRELAX
          WRITE (LUPRI,*) 'CC_XIETA> LABELH,ISYHOP = ',LABELH,ISYHOP
          WRITE (LUPRI,*) 'CC_XIETA> IATOM,ISCOOR  = ',IATOM,ISCOOR
          WRITE (LUPRI,*) 'CC_XIETA> ICORSY,ICOOR  = ',ICORSY,ICOOR
        END IF

*---------------------------------------------------------------------*
*     loop over number of distributions in this shell:
*     (we enter this loop only for those ITRAN where we have
*      a two-electron contribution, i.e. if LRELAX is true. )
*---------------------------------------------------------------------*
      IF (LRELAX) THEN

        DO IDEL2 = 1, NUMDIS
C
C         set orbital index and symmetry class for next 
C         AO integral distribution:
C
          IF (DIRECT) THEN
            IDEL   = INDEXA(IDEL2)
            IF (NOAUXB) THEN
              IDUM = 1
              CALL IJKAUX(IDEL,IDUM,IDUM,IDUM)
            END IF
            ISYDEL = ISAO(IDEL)
          ELSE
            IDEL   = IBAS(ISYMD1) + ILLL
            ISYDEL = ISYMD1
          END IF

          ISY0DIS = MULD2H(ISYDEL,ISYM0)
          ISY1DIS = MULD2H(ISYDEL,ISYHOP)
C
C         read undifferentiated AO integral distribution,
C         and transform gamma index to occupied:
C
          K0XINT  = KEND3
          KD0PRHF = K0XINT  + NDISAO(ISY0DIS)
          KD0HRHF = KD0PRHF + NDSRHF(ISY0DIS)
          KEND4   = KD0HRHF + NDSRHF(ISY0DIS)
          LWRK4   = LWORK   - KEND4

          IF (LWRK3 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CC_XIETA. (CCRDAOD)')
          END IF

COLD
C         CALL CCRDAO(WORK(K0XINT),IDEL,IDEL2,
C    &                WORK(KEND4),LWRK4,WORK(KRECNR),DIRECT)
COLD

C
C         read all gamma for given delta (``3-index'' approach)
C
          NUMD = 1
          NUMG = 0
          DO ISYMG = 1, NSYM
            DO G = 1, NBAS(ISYMG)
              NUMG       = NUMG + 1
              IGAM(NUMG) = G + IBAS(ISYMG)
            END DO
          END DO

          CALL CC_RDAOD(WORK(K0XINT),0,1,1,IGAM,IDEL,NUMG,NUMD,
     &                  WORK(KEND4),LWRK4,WORK(KRECNR),WORK(KNRECS),
     &                  NBUFMX,NFILES,DIRECT,0,1,LDERINT)

          IOPT = 0
          CALL CCTRBT2(WORK(K0XINT),WORK(KD0PRHF),WORK(KLAMP0),ISYM0,
     &              WORK(KEND4),LWRK4,ISY0DIS,IOPT,.FALSE.,.FALSE.,ONE)

          IOPT = 0
          CALL CCTRBT2(WORK(K0XINT),WORK(KD0HRHF),WORK(KLAMH0),ISYM0,
     &              WORK(KEND4),LWRK4,ISY0DIS,IOPT,.FALSE.,.FALSE.,ONE)


          IF (LTWOEL) THEN
C
C            read first derivative AO integral distribution:
C            and transform gamma index to occupied:
C
             IF       ( LABELH(1:5).EQ.'HAM0 ' ) THEN
               ITYPE  = 0
               SQRINT = .FALSE.
             ELSE IF ( LABELH(1:5).EQ.'1DHAM' ) THEN
               ITYPE  = 1
               SQRINT = .FALSE.
             ELSE IF ( LABELH(1:5).EQ.'dh/dB' ) THEN
               ITYPE  = 5
               SQRINT = .TRUE.
             ELSE
               CALL QUIT('ITYPE unknown in CC_XIETA.')
             END IF

             K1XINT  = KEND4
             IF (SQRINT) THEN
               KD1PRHF = K1XINT  + NDISAOSQ(ISY1DIS)
               KD1HRHF = KD1PRHF + NDSRHFSQ(ISY1DIS)
               KEND4   = KD1HRHF + NDSRHFSQ(ISY1DIS)
             ELSE
               KD1PRHF = K1XINT  + NDISAO(ISY1DIS)
               KD1HRHF = KD1PRHF + NDSRHF(ISY1DIS)
               KEND4   = KD1HRHF + NDSRHF(ISY1DIS)
             END IF
             LWRK4   = LWORK   - KEND4

             IF (LWRK3 .LT. 0) THEN
               CALL QUIT('Insufficient work space in '//
     &               'CC_XIETA. (CCRDAOD)')
             END IF
C
C            read all gamma for given delta (``3-index'' approach)
C
             NUMD = 1
             NUMG = 0
             DO ISYMG = 1, NSYM
               DO G = 1, NBAS(ISYMG)
                 NUMG       = NUMG + 1
                 IGAM(NUMG) = G + IBAS(ISYMG)
               END DO
             END DO

             CALL CC_RDAOD(WORK(K1XINT),ISCOOR,ICOOR,ICORSY,
     &                     IGAM,IDEL,NUMG,NUMD,
     &                     WORK(KEND4),LWRK4,WORK(KRECNR),WORK(KNRECS),
     &                     NBUFMX,NFILES,DIRECT,ITYPE,MXCOMP,LDERINT)

             IOPT = 0
             CALL CCTRBT2(WORK(K1XINT),WORK(KD1PRHF),WORK(KLAMP0),ISYM0,
     &                WORK(KEND4),LWRK4,ISY0DIS,IOPT,SQRINT,.FALSE.,ONE)

             IOPT = 0
             CALL CCTRBT2(WORK(K1XINT),WORK(KD1HRHF),WORK(KLAMH0),ISYM0,
     &                WORK(KEND4),LWRK4,ISY0DIS,IOPT,SQRINT,.FALSE.,ONE)


             IF (LOCDBG) THEN
               WRITE (LUPRI,*) 'Norm^2 of undifferentiated integrals:',
     &          DDOT(NDISAO(ISY0DIS),WORK(K0XINT),1,WORK(K0XINT),1)
              IF (SQRINT) THEN
               WRITE (LUPRI,*) 
     &            'Norm^2 of derivative integrals (squared):',
     &            DDOT(NDISAOSQ(ISY1DIS),WORK(K1XINT),1,WORK(K1XINT),1)
              ELSE
               WRITE (LUPRI,*) 
     &            'Norm^2 of derivative integrals (packed):',
     &            DDOT(NDISAO(ISY1DIS),WORK(K1XINT),1,WORK(K1XINT),1)
              END IF
             END IF

          END IF

*---------------------------------------------------------------------*
*         if in this batch we have not yet calculated the zero-order
*         intermediates for relaxed/derivative vectors, we test if
*         can (have to) calculated them now.
*
*         the intermediates needed for the Eta vector which depend
*         only on the left vector and the (zero-order) amplitudes
*         have to recalculated each time the left vector changes
*---------------------------------------------------------------------*
          IF ( (.NOT.LZEROIMDONE) .AND. LRELAX) THEN
            LZERO = .TRUE.
          ELSE
            LZERO = .FALSE.
          END IF

C         LZEROLFT = (.NOT.SKIPETA) .AND. (IDLSTL .NE. IDLSTL_OLD)
          LZEROLFT = LRELAX .OR. LTWOEL

*---------------------------------------------------------------------*
*        calculate intermediates needed for XI: 
*         - some zeroth-order intermediates need only be computed once
*           i.e. when the LZERO flag is set
*         - the others have to be computed whenever IOPER or IRELAX 
*           change
*---------------------------------------------------------------------*

          LNEWXI = (IOPER.NE.IOPER_OLD) .OR. (IRELAX.NE.IRELAX_OLD)
     &             .OR. LZERO

          IF (LOCDBG) THEN
            WRITE (LUPRI,*) 'ITRAN      :',ITRAN
            WRITE (LUPRI,*) 'IOPER      :',IOPER
            WRITE (LUPRI,*) 'LABELH     :',LABELH
            WRITE (LUPRI,*) 'IRELAX     :',IRELAX
            WRITE (LUPRI,*) 'LZEROIMDONE:',LZEROIMDONE
            WRITE (LUPRI,*) 'LZERO      :',LZERO
            WRITE (LUPRI,*) 'LRELAX     :',LRELAX
            WRITE (LUPRI,*) 'LNEWXI     :',LNEWXI
            WRITE (LUPRI,*) 
     &            'KG0IM :',KG0IM,       KG0IM       +NT1AO(ISYM0) 
            WRITE (LUPRI,*) 
     &            'KG0IM :',KGBIM(ITRAN),KGBIM(ITRAN)+NT1AO(ISYHOP)
          END IF

          CALL CCXIINTAO(WORK(K0XINT),WORK(KD0PRHF),ISY0DIS,
     &                   WORK(K1XINT),WORK(KD1PRHF),ISY1DIS,SQRINT,
     &                   WORK(KLAMP0),WORK(KLAMH0),
     &                   WORK(KLAMDPQ(ITRAN)),WORK(KLAMDHQ(ITRAN)),
     &                   WORK(KDENS0),        WORK(KDPK0),
     &                   WORK(KDNSHF0),       WORK(KDPKHF0),
     &                   WORK(KDENSB(ITRAN)), WORK(KDPCKB(ITRAN)),
     &                   WORK(KDNSHFB(ITRAN)),WORK(KDPKHFB(ITRAN)),
     &                   WORK(KFOCK0),
     &                   WORK(KFOCKB(ITRAN)), WORK(KFCKHFB(ITRAN)),
     &                   WORK(KRHO2(ITRAN)),WORK(KG0IM),
     &                   WORK(KGBIM(ITRAN)),
     &                   LUBFDX0,FNBFDX0,IADRX0,
     &                   LUBFDX1,FNBFDX1,IADRX1(1,ITRAN),
     &                   LU0IAJB,LU0IABJ,LU0IJBA,LU0AIBJ,
     &                   FN0IAJB,FN0IABJ,FN0IJBA,FN0AIBJ,
     &                   LU1IAJB,LU1IABJ,LU1IJBA,LU1AIBJ,
     &                   FN1IAJB,FN1IABJ,FN1IJBA,FN1AIBJ,
     &                   IT2DEL0,IADRI0,IT2DELB,IADRIB,
     &                   ITRAN,IDEL,LZERO,LNEWXI,LRELAX,LTWOEL,
     &                   IREAL,ISYHOP,WORK(KEND4), LWRK4)

*---------------------------------------------------------------------*
*        calculate additional intermediates needed for ETA: 
*         - some zeroth-order intermediates need only be computed once
*         - the others have to be computed whenever IOPER, IRELAX or
*           IDLSTL change, which should be the case for every ITRAN
*           for which a ETA vector is requested
*---------------------------------------------------------------------*

          IF (.NOT. SKIPETA) THEN
            ISYCTR = ILSTSYM(LISTL,IDLSTL)  ! symmetry of ZETA vector

            CALL CCETAINT2(WORK(K0XINT),WORK(KD0HRHF), ISY0DIS,
     &                     WORK(K1XINT),WORK(KD1HRHF), ISY1DIS,SQRINT,
     &                     WORK(KBFZ0(ITRAN)), WORK(KBFZI(ITRAN)), 
     &                     WORK(KF0IM(ITRAN)), WORK(KFBIM(ITRAN)),
     &                     WORK(KZFCK0(ITRAN)),WORK(KZFCKB(ITRAN)),
     &                     WORK(KLAMP0), WORK(KLAMH0), 
     &                     WORK(KLAMDPQ(ITRAN)), WORK(KLAMDHQ(ITRAN)), 
     &                     WORK(KCHI(ITRAN)),  WORK(KCHIQ(ITRAN)),
     &                     WORK(KZDPK0(ITRAN)),WORK(KZDEN0(ITRAN)),
     &                     WORK(KZDPKB(ITRAN)),WORK(KZDENB(ITRAN)),
     &                     IREAL,   ISYHOP,  ISYCTR,
     &                     FNBFDE0, LUBFDE0, IADRE0(1,ITRAN),  
     &                     FNBFDE1, LUBFDE1, IADRE1(1,ITRAN),  
     &                     FILPQ0,  LUPQ0,   IADRPQ0(1,ITRAN),  
     &                     FILPQ1,  LUPQ1,   IADRPQ1(1,ITRAN),  
     &                     FNHINT0, LUHINT0, IADRH0(1,ITRAN), IADRZ0,
     &                     FNHINT1, LUHINT1, IADRH1(1,ITRAN), IADRZ1,
     &                     LRELAX,  LTWOEL,  LZEROLFT,
     &                     IDEL,    WORK(KEND4), LWRK4)

          END IF

*---------------------------------------------------------------------*
*     close loop over AO two-electron integrals and ITRAN
*---------------------------------------------------------------------*
        END DO ! IDEL2

          IF (LZERO) LZEROIMDONE = .TRUE.

        END IF ! LRELAX

        IOPER_OLD  = IOPER
        IF (SKIPETA) THEN
           IDLSTL_OLD = -1
        ELSE
           IDLSTL_OLD = IDLSTL
        END IF
        IRELAX_OLD = IRELAX

      END DO ! ITRAN 
      END DO ! ILLL
      END DO ! ISYMD1
      END IF ! (LLTWOEL .OR. LLRELAX) 

*---------------------------------------------------------------------*
* finished with intermediates depending two-electron AO integrals
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'CC_XIETA> after loop over AO'
        CALL FLSHFO(LUPRI)
      END IF

*---------------------------------------------------------------------*
* close & delete the files with eff. densities:
*---------------------------------------------------------------------*
      CALL WCLOSE2(LUBFDX0, FNBFDX0, 'KEEP'  )
      CALL WCLOSE2(LUBFDE0, FNBFDE0, 'KEEP'  )
      CALL WCLOSE2(LUBFDX1, FNBFDX1, 'DELETE')
      CALL WCLOSE2(LUBFDE1, FNBFDE1, 'DELETE')

*---------------------------------------------------------------------*
*     calculate R intermediate from (ia|j del) integrals
*     and R-bar intermediate from (ia|j del)-bar integrals
*---------------------------------------------------------------------*
      IF ( (LLTWOEL.OR.LLRELAX) .AND. (.NOT.CCS) )THEN

         KT2AM = KEND1
         KEND2 = KT2AM + NT2SQ(ISYM0)
         LWRK2 = LWORK - KEND2

         KXIAJB = KEND2

         IF (LWRK2 .LT. MAX(NT2AM(ISYM0),MT2BCD) ) THEN
           CALL QUIT('Insufficient work space in CC_XIETA. '//
     &           '(RIM section)')
         END IF

C        ---------------------------------------------------
C        read zeroth-order cluster amplitudes, square up and
C        calculate 2T2(ia,jb) - T2(ja,ib) in place:
C        ---------------------------------------------------
         IOPT = 2
         CALL CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,DUMMY,WORK(KEND2))
         CALL CC_T2SQ(WORK(KEND2),WORK(KT2AM),ISYM0)
         CALL CCRHS_T2TR(WORK(KT2AM),WORK(KEND2),LWRK2,ISYM0)

C        ---------------------------------------------------
C        initialize R intermediates:
C        ---------------------------------------------------
         CALL DZERO(WORK(KR0IM),NEMAT1(ISYM0))

         KRBIMOLD = -1
         DO ITRAN = ISTART(IBATCH), IEND(IBATCH)
            IOPER   = IXETRAN(1,ITRAN)   ! operator index
            ISYHOP = ISYOPR(IOPER)      ! symmetry of hamiltonian
            LTWOEL  = LPDBSOP(IOPER)     ! pert.-dep. basis set
            IRELAX  = IXETRAN(5,ITRAN)   ! flag for orbital relax.
            IF ( LTWOEL .OR. (IRELAX.GE.1)) THEN
              IF (KRBIM(ITRAN).NE.KRBIMOLD) THEN
                KRBIMOLD = KRBIM(ITRAN)
                CALL DZERO(WORK(KRBIM(ITRAN)),NEMAT1(ISYHOP))
              END IF
            END IF
         END DO

         DO IDEL = 1, NBAST

            ISYDEL  = ISAO(IDEL)
            ISY0IAJ = MULD2H(ISYDEL,ISYM0)
            LEN0    = NT2BCD(ISY0IAJ)
            IADR0   = IT2DEL0(IDEL)

C           ------------------------------------------------------
C           calculate R intermediates form zeroth-order integrals:
C           ------------------------------------------------------
C           IF (IBATCH.EQ.1) THEN
               CALL GETWA2(LU0IAJB,FN0IAJB,WORK(KXIAJB),IADR0,LEN0)
               CALL CC_RIM(WORK(KR0IM),WORK(KT2AM),ISYM0,
     &                     WORK(KXIAJB),ISY0IAJ,IDEL,ISYDEL)
C           END IF

            KRBIMOLD = -1
            DO ITRAN = ISTART(IBATCH), IEND(IBATCH)
              IF (KRBIM(ITRAN).NE.KRBIMOLD) THEN
               KRBIMOLD = KRBIM(ITRAN)

               IOPER   = IXETRAN(1,ITRAN)   ! operator index
               LTWOEL  = LPDBSOP(IOPER)     ! pert.-dep. basis set
               IRELAX  = IXETRAN(5,ITRAN)   ! flag for orbital relax.

               IF ( LTWOEL .OR. (IRELAX.GE.1)) THEN
                 ISYHOP  = ISYOPR(IOPER)      
                 ISY1IAJ = MULD2H(ISYDEL,ISYHOP)
                 LEN1    = NT2BCD(ISY1IAJ)
                 IADRB   = IT2DELB(IDEL,ITRAN) 

C                ----------------------------------------------------
C                calculate R intermediates from derivative integrals:
C                ----------------------------------------------------
                 CALL GETWA2(LU1IAJB,FN1IAJB,WORK(KXIAJB),IADRB,LEN1)
                 CALL CC_RIM(WORK(KRBIM(ITRAN)),WORK(KT2AM),ISYM0,
     &                       WORK(KXIAJB),ISY1IAJ,IDEL,ISYDEL)

               END IF 
              END IF 

            END DO

         END DO

C           DO ITRAN = ISTART(IBATCH), IEND(IBATCH)
C              IOPER   = IXETRAN(1,ITRAN)   ! operator index
C              LABELH = LBLOPR(IOPER)      ! operator label
C              LTWOEL  = LPDBSOP(IOPER)     ! pert.-dep. basis set
C              IRELAX  = IXETRAN(5,ITRAN)   ! flag for orbital relax.
C              WRITE (LUPRI,*) 'IOPER,LABELH:',IOPER,LABELH
C              WRITE (LUPRI,*) 'LTWOEL,IRELAX:',LTWOEL,IRELAX
C              IF ( LTWOEL .OR. (IRELAX.GE.1)) THEN
C                 WRITE (LUPRI,*) 'RBIM matrix:'
C                 CALL OUTPUT(WORK(KRBIM(ITRAN)),1,NVIRT,1,NBAST,
C    &                        NVIRT,NBAST,1,LUPRI)
C              END IF
C           END DO

      END IF ! (LLTWOEL .OR. LLRELAX) 

*---------------------------------------------------------------------*
* write the BF and BZeta intermediates to file: 
* (this overwrites the eff. 'densities', which are no longer needed)
*---------------------------------------------------------------------*
      IF ( (LLTWOEL.OR.LLRELAX) .AND. CCSD )THEN

         LUBFIM = -1
         CALL WOPEN2(LUBFIM, FNBFIM, 64, 0)

         IADR = 1
         DO ITRAN = ISTART(IBATCH), IEND(IBATCH)
           IOPER  = IXETRAN(1,ITRAN)   ! operator index
           LABELH = LBLOPR(IOPER)      ! operator label
           ISYHOP = ISYOPR(IOPER)      ! symmetry of hamiltonian

           LTWOEL = LPDBSOP(IOPER)     ! two-electron contribution
           IRELAX = IXETRAN(5,ITRAN)   ! relaxation flag
           LRELAX = LTWOEL .OR. (IRELAX.GE.1)

           IF (IORDER.GE.2) THEN
             IRELAX2 = IXETRAN(6,ITRAN)
             LRELAX  = LTWOEL .OR. ((IRELAX.GE.1) .AND. (IRELAX.GE.2))
           END IF

           IF (LRELAX) THEN
             LEN    = 2 * NT2ORT(ISYHOP) ! length of BF intermediate
             CALL PUTWA2(LUBFIM,FNBFIM,WORK(KRHO2(ITRAN)),IADR,LEN)
             IADRX1(1,ITRAN) = IADR
             IADR = IADR + LEN
           END IF
         END DO


         ! open the file for the BFZ intermediates
         LUBFZI = -1
         CALL WOPEN2(LUBFZI,FNBFZI,64,0)

         IADR = 1
C
         DO ITRAN = ISTART(IBATCH), IEND(IBATCH)
           IOPER  = IXETRAN(1,ITRAN)   ! operator index
           IDLSTL = IXETRAN(2,ITRAN)   ! ZETA vector index
           ISYHOP = ISYOPR(IOPER)      ! symmetry of hamiltonian
           LTWOEL = LPDBSOP(IOPER)     ! two-electron contribution

           IRELAX = IXETRAN(5,ITRAN)   ! relaxation flag
           LRELAX = LTWOEL .OR. (IRELAX.GE.1)
           IF (IORDER.GE.2) THEN
             IRELAX2 = IXETRAN(6,ITRAN)
             LRELAX  = LTWOEL .OR. ((IRELAX.GE.1) .AND. (IRELAX.GE.2))
           END IF

           SKIPETA = ( IXETRAN(4,ITRAN) .EQ. -1 )
    
C          IF (.NOT. SKIPETA) THEN
C          END IF
   
           IF (LRELAX .AND. .NOT.SKIPETA) THEN
             ISYCTR = ILSTSYM(LISTL,IDLSTL) ! symmetry of ZETA vector
             LEN    = 2 * NT2ORT(ISYCTR)    ! length of BFZ0 intermed.
             CALL PUTWA2(LUBFZI,FNBFZI,WORK(KBFZ0(ITRAN)),IADR,LEN)
             IADRE0(1,ITRAN) = IADR
             IADR = IADR + LEN

             ISYETA = MULD2H(ISYHOP,ISYCTR) ! symmetry of ETA vector
             LEN    = 2 * NT2ORT(ISYETA)    ! length of BFZ1 intermed.
             CALL PUTWA2(LUBFZI,FNBFZI,WORK(KBFZI(ITRAN)),IADR,LEN)
             IADRE1(1,ITRAN) = IADR
             IADR = IADR + LEN
           END IF
         END DO


      END IF ! (LLTWOEL .OR. LLRELAX) 

*---------------------------------------------------------------------*
*     for CC2 calculate RZeta intermediate from (ai|del j) integrals
*     and RZeta-bar intermediate from (ai|del j)-bar integrals
*---------------------------------------------------------------------*
      IF ( (LLTWOEL.OR.LLRELAX) .AND. CC2 )THEN

C        ---------------------------------------------------
C        allocate and initialize RZeta intermediates:
C        ---------------------------------------------------
         IDLSTL_OLD = -1
         DO ITRAN = ISTART(IBATCH), IEND(IBATCH)
            IOPER   = IXETRAN(1,ITRAN)   ! operator index
            IDLSTL  = IXETRAN(2,ITRAN)   ! ZETA vector index
            IRELAX  = IXETRAN(5,ITRAN)   ! flag for orbital relax.
            ISYHOP  = ISYOPR(IOPER)      ! symmetry of hamiltonian
            LTWOEL  = LPDBSOP(IOPER)     ! pert.-dep. basis set
            SKIPETA = ( IXETRAN(4,ITRAN) .EQ. -1 )
            LRELAX  = LTWOEL .OR. (IRELAX.GE.1)

            IF ( LRELAX .AND. (.NOT.SKIPETA) ) THEN
              ISYCTR = ILSTSYM(LISTL,IDLSTL) ! symmetry of ZETA vector
              ISYETA = MULD2H(ISYCTR,ISYHOP) ! symmetry of RZBI interm.
              IF (IDLSTL .NE. IDLSTL_OLD) THEN
                 IDLSTL_OLD   = IDLSTL
                 KRZ0I(ITRAN) = KEND1
                 KEND1        = KRZ0I(ITRAN) + NEMAT1(ISYCTR)
                 IF (KEND1.GT.LWORK) THEN
                   CALL QUIT('Insufficient memory in CC_XIETA. (RZ0)')
                 END IF
                 CALL DZERO(WORK(KRZ0I(ITRAN)),NEMAT1(ISYCTR))
              ELSE
                 KRZ0I(ITRAN) = KRZ0I(ITRAN-1) ! promote last pointer
              END IF
              KRZBI(ITRAN) = KEND1
              KEND1        = KRZBI(ITRAN) + NEMAT1(ISYETA)
              IF (KEND1.GT.LWORK) THEN
                CALL QUIT('Insufficient memory in CC_XIETA. (RZI)')
              END IF
              CALL DZERO(WORK(KRZBI(ITRAN)),NEMAT1(ISYETA))
            ELSE
              IDLSTL_OLD   = -1
              KRZ0I(ITRAN) = KRZ0I(ITRAN-1) ! promote last pointer
            END IF

         END DO

C        ---------------------------------------------------
C        start new loop to do the actual calculations:
C        ---------------------------------------------------
         IDLSTL_OLD = -1
         DO ITRAN = ISTART(IBATCH), IEND(IBATCH)
            IOPER   = IXETRAN(1,ITRAN)   ! operator index
            IDLSTL  = IXETRAN(2,ITRAN)   ! ZETA vector index
            IRELAX  = IXETRAN(5,ITRAN)   ! flag for orbital relax.
            ISYHOP  = ISYOPR(IOPER)      ! symmetry of hamiltonian
            LTWOEL  = LPDBSOP(IOPER)     ! pert.-dep. basis set
            SKIPETA = ( IXETRAN(4,ITRAN) .EQ. -1 )
            LRELAX  = LTWOEL .OR. (IRELAX.GE.1)

            IF ( LRELAX .AND. (.NOT.SKIPETA) ) THEN
              ISYCTR = ILSTSYM(LISTL,IDLSTL) ! symmetry of ZETA vector
              ISYETA = MULD2H(ISYCTR,ISYHOP) ! symmetry of RZBI interm.

              ! if index of zeta vector changed, get the new one
              IF (IDLSTL .NE. IDLSTL_OLD) THEN
                KZETA2 = KEND1
                KEND2  = KZETA2 + NT2SQ(ISYCTR)
                LWRK2  = LWORK  - KEND2

                KXAIBJ = KEND2

                IF ( LWRK2 .LT. MAX(NT2AM(ISYCTR),MT2BCD) ) THEN
                  CALL QUIT('Insufficient memory in CC_XIETA. (RZ)')
                END IF
   
                IOPT = 2
                CALL CC_RDRSP(LISTL,IDLSTL,ISYCTR,IOPT,MODEL,
     &                        DUMMY,WORK(KEND2))
                CALL CC_T2SQ(WORK(KEND2),WORK(KZETA2),ISYCTR)
              END IF

              DO IDEL = 1, NBAST

C                ------------------------------------------------------
C                calculate RZeta intermed. form zeroth-order integrals:
C                ------------------------------------------------------
                 IF (IDLSTL .NE. IDLSTL_OLD) THEN
                    ISYDEL  = ISAO(IDEL)
                    ISY0IAJ = MULD2H(ISYDEL,ISYM0)
                    LEN0    = NT2BCD(ISY0IAJ)
                    IADR0   = IT2DEL0(IDEL)
                    CALL GETWA2(LU0AIBJ,FN0AIBJ,WORK(KXAIBJ),IADR0,LEN0)
                    CALL CC_RIM(WORK(KRZ0I(ITRAN)),WORK(KZETA2),ISYCTR,
     &                          WORK(KXAIBJ),ISY0IAJ,IDEL,ISYDEL)
                 END IF

C                ----------------------------------------------------
C                calculate RZeta intermed. from derivative integrals:
C                ----------------------------------------------------
                 ISY1IAJ = MULD2H(ISYDEL,ISYHOP)
                 LEN1    = NT2BCD(ISY1IAJ)
                 IADRB   = IT2DELB(IDEL,ITRAN) 
                 CALL GETWA2(LU1AIBJ,FN1AIBJ,WORK(KXAIBJ),IADRB,LEN1)
                 CALL CC_RIM(WORK(KRZBI(ITRAN)),WORK(KZETA2),ISYETA,
     &                       WORK(KXAIBJ),ISY1IAJ,IDEL,ISYDEL)

              END DO

              IDLSTL_OLD = IDLSTL

            ELSE
              IDLSTL_OLD = -1
            END IF ! ( LRELAX .AND. (.NOT.SKIPETA) ) 

         END DO

         IF (LOCDBG .AND. NSYM.EQ.1) THEN
            DO ITRAN = ISTART(IBATCH), IEND(IBATCH)
               IOPER   = IXETRAN(1,ITRAN)   ! operator index
               IDLSTL  = IXETRAN(2,ITRAN)   ! ZETA vector index
               LABELH  = LBLOPR(IOPER)      ! operator label
               LTWOEL  = LPDBSOP(IOPER)     ! pert.-dep. basis set
               IRELAX  = IXETRAN(5,ITRAN)   ! flag for orbital relax.
               SKIPETA = ( IXETRAN(4,ITRAN) .EQ. -1 )
               LRELAX  = LTWOEL .OR. (IRELAX.GE.1)
               IF ( LRELAX .AND. (.NOT.SKIPETA) ) THEN
                  ISYCTR = ILSTSYM(LISTL,IDLSTL) 
                  ISYETA = MULD2H(ISYCTR,ISYHOP) 
                  WRITE (LUPRI,*) 'IOPER,LABELH  :',IOPER,LABELH
                  WRITE (LUPRI,*) 'LTWOEL,IRELAX :',LTWOEL,IRELAX
                  WRITE (LUPRI,*) 'IDLSTL,SKIPETA:',IDLSTL,SKIPETA
                  WRITE (LUPRI,*) 'RZeta0 matrix:',KRZ0I(ITRAN)
                  CALL OUTPUT(WORK(KRZ0I(ITRAN)),1,NVIRT,1,NBAST,
     &                        NVIRT,NBAST,1,LUPRI)
                  WRITE (LUPRI,*) 'RZetaB matrix:',KRZBI(ITRAN)
                  CALL OUTPUT(WORK(KRZBI(ITRAN)),1,NVIRT,1,NBAST,
     &                        NVIRT,NBAST,1,LUPRI)
               END IF
            END DO
         END IF

      END IF ! ( (LLTWOEL.OR.LLRELAX) .AND. CC2 )

*---------------------------------------------------------------------*
* open files for C and D intermediates:
*---------------------------------------------------------------------*
      LUCIM = -1
      LUDIM = -1
      CALL WOPEN2(LUCIM, FNCIM, 64, 0)
      CALL WOPEN2(LUDIM, FNDIM, 64, 0)

*---------------------------------------------------------------------*
* start a new loop over the transformations and calculate now the
* complete vectors from the intermediates:
*---------------------------------------------------------------------*
      KEND1SV = KEND1
      LWRK1SV = LWRK1

      DO ITRAN = ISTART(IBATCH), IEND(IBATCH)
          
        IOPER  = IXETRAN(1,ITRAN)  ! operator index
        IDLSTL = IXETRAN(2,ITRAN)  ! ZETA vector index
        IRELAX = IXETRAN(5,ITRAN)  ! flag for relax. contrib.

        ISYHOP = ISYOPR(IOPER)     ! symmetry of hamiltonian
        LABELH = LBLOPR(IOPER)     ! operator label
        LTWOEL = LPDBSOP(IOPER)    ! two-electron contribution

        SKIPETA = ( IXETRAN(4,ITRAN) .EQ. -1 )
        LRELAX  = LTWOEL .OR. (IRELAX.GE.1)

        IF (IORDER.GE.2) THEN
           IRELAX1 = IXETRAN(5,ITRAN)
           IRELAX2 = IXETRAN(6,ITRAN) 

           ! get info about the different contributing perturbations
           LABEL1 = '?'
           IF (IRELAX1 .GE. 1) LABEL1 = LRTHFLBL(IRELAX1)
           LABEL2 = '?'
           IF (IRELAX2 .GE. 1) LABEL2 = LRTHFLBL(IRELAX2)
            
           IF (IORDER.EQ.2) THEN
              INUM   = IROPER2(LABEL1,LABEL2,LABELH,ISIGNH,ISYM)
              IF (LABEL1(1:1).NE.'?') THEN
                 IOPER1 = IROPER(LABEL1,ISYOP1)
                 LPDBS1 = LPDBSOP(IOPER1) 
              END IF
              IF (LABEL2(1:1).NE.'?') THEN
                 IOPER2 = IROPER(LABEL2,ISYOP2)
                 LPDBS2 = LPDBSOP(IOPER2) 
              END IF
              ! for second-order operators we set the LRELAX flag
              ! only if both perturbations are either relaxed or
              ! have field-dep. basis sets, i.e. when we get a
              ! contrib. from the undifferentiated 2e- integrals
              LRELAX  = (LPDBS1 .OR. (IRELAX1.GE.1)) .AND.
     &                  (LPDBS2 .OR. (IRELAX2.GE.1)) 
           ELSE
              CALL QUIT('IORDER out of range in CC_XIETA.')
           END IF

        END IF


C       SKIPXI = ( IXETRAN(3,ITRAN) .EQ. -1 )
C       SKIPXI = .FALSE.

        ! for the Cauchy vectors we need not to compute the usual
        ! CCSD part of the Xi vectors but just some intermediates
        ! as e.g. F^B Fock matrix
        SKIPXI = LCAUCHY  
C

        IF (.NOT. SKIPETA) THEN
          ISYCTR = ILSTSYM(LISTL,IDLSTL) ! sym. of ZETA vector
          ISYETA = MULD2H(ISYHOP,ISYCTR) ! sym. of ETA result vector
        END IF

        IF (LOCDBG) THEN
          WRITE (LUPRI,*) 'CC_XIETA> IOPER, IRELAX = ',IOPER,IRELAX
          WRITE (LUPRI,*) 'CC_XIETA> LTWOEL,LRELAX = ',LTWOEL,LRELAX
          WRITE (LUPRI,*) 'CC_XIETA> LABELH,ISYHOP = ',LABELH,ISYHOP
          IF (IORDER.GE.2) THEN
            WRITE (LUPRI,*) 'CC_XIETA> IRELAX2 = ', IRELAX2
            WRITE (LUPRI,*) 'CC_XIETA> LABEL1,LABEL2 = ',LABEL1,LABEL2
            WRITE (LUPRI,*) 'CC_XIETA> IOPER1,IOPER2 = ',IOPER1,IOPER2
          END IF
          WRITE (LUPRI,*) 'CC_XIETA> SKIPXI        = ',SKIPXI
          WRITE (LUPRI,*) 'WORK(0) = ',WORK(ITST)
          CALL FLSHFO(LUPRI)
        END IF

        KEND1 = KEND1SV

        KEMAT1 = KEND1
        KEMAT2 = KEMAT1 + NMATAB(ISYHOP)
        KEND1  = KEMAT2 + NMATIJ(ISYHOP)
        LWRK1  = LWORK - KEND1
        
        ISYRES = ISYHOP ! symmetry of result vector

        IF (CCSD .AND. (LTWOEL.OR.LRELAX)) THEN
           KGAMMA = KEND1
           KEND1  = KGAMMA + NGAMMA(ISYHOP)
        END IF
        IF (CCSDT) THEN
           KFOCKBAO = KEND1
           KEND1    = KFOCKBAO + N2BST(ISYRES)
        END IF

        KXI1  = KEND1
        KEND2 = KXI1 + NT1AM(ISYRES)
        IF ( CC2 .OR. CCSD .OR. CCSDT) THEN
           KXI2   = KEND2
           KEND2  = KXI2  + NT2AM(ISYRES)
        END IF

        LWRK2 = LWORK - KEND2
        IF (LWRK2 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CC_XIETA. (1b)')
        END IF

        IF (CCSDT) THEN
           ! save F^B in AO basis needed for test in CCSDT Eta
           CALL DCOPY(N2BST(ISYRES),WORK(KFOCKB(ITRAN)),1,
     &                              WORK(KFOCKBAO),1)
        END IF

        CALL DZERO(WORK(KXI1),NT1AM(ISYRES))

        IADRDI = 1
        IADRCI = 1

        CALL CCXI2(WORK(KXI1),WORK(KXI2), IORDER, SKIPXI,
     &             WORK(KGAMMA),WORK(KEMAT1), WORK(KEMAT2),
     &             WORK(KG0IM), WORK(KGBIM(ITRAN)),
     &             WORK(KR0IM), WORK(KRBIM(ITRAN)), 
     &             WORK(KLAMP0), WORK(KLAMH0), WORK(KCMO0), 
     &             WORK(KONEH0), WORK(KFOCK0),
     &             WORK(KLAMDPQ(ITRAN)),WORK(KLAMDHQ(ITRAN)),
     &             WORK(KCMOPQ(ITRAN)), WORK(KCMOHQ(ITRAN)),
     &             WORK(KONEHB(ITRAN)), WORK(KFOCKB(ITRAN)), 
     &             WORK(KFCKHFB(ITRAN)),ISYHOP,
     &             LABEL1,IRELAX1,ISYOP1,
     &             LABEL2,IRELAX2,ISYOP2,
     &             FNBFIM,  LUBFIM,  IADRX1(1,ITRAN),
     &             FN0IAJB, FN0IJBA, FN0IABJ, FN0AIBJ,
     &             LU0IAJB, LU0IJBA, LU0IABJ, LU0AIBJ,
     &             FN1IAJB, FN1IJBA, FN1IABJ, FN1AIBJ,
     &             LU1IAJB, LU1IJBA, LU1IABJ, LU1AIBJ,
     &             LUDIM,FNDIM,IADRDI,
     &             LUCIM,FNCIM,IADRCI,
     &             IT2DEL0, IT2DELB(1,ITRAN), IADRI0,
     &             LRELAX,  LTWOEL,  WORK(KEND2), LWRK2)

        ! redefine SKIPXI: things are only skip if the result vectors
        ! for this section are not needed/wanted.
        SKIPXI = ( IXETRAN(3,ITRAN) .EQ. -1 )

        IF (CCR12 .AND. CCS) THEN
          CONTINUE
        ELSE IF (CCR12 .AND. (.NOT.SKIPXI)) THEN
          IF (LRELAX) THEN
            CALL QUIT('CC-R12 response can only handle unrelaxed '//
     &                'orbitals')
          END IF

C------------------------
C         allocate memory
C------------------------
          KXIR12 = KEND2
          KEND2  = KXIR12 + NTR12AM(ISYRES)

          kxir12sq= kend2
          ktr12   = kxir12sq + ntr12sq(isyres)
          ktr12sq = ktr12 + ntr12am(1)
          kxmat   = ktr12sq + ntr12sq(1)
          kxmatsq = kxmat + nr12r12p(1)
          kvxintsq= kxmatsq + nr12r12sq(1) 
          kprpao  = kvxintsq + nr12r12sq(ISYRES)
          kend3   = kprpao + N2BST(ISYRES)
          LWRK3   = LWORK - KEND3
          IF (LWRK3 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CC_XIETA. (1b-R12)')
          END IF

C------------------------------------------------------------------
C         Read R12 amplitudes from disk and reorder to full square:
C------------------------------------------------------------------
          iopt=32
          call cc_rdrsp('R0 ',0,1,iopt,model,dummy,work(ktr12))
         
C         if (locdbg) then
C           call around('R12 amplitudes in CC_XIETA1 before unpacking')
C           call cc_prpr12(work(ktr12),1,1,.false.)
C         end if
        
          iopt = 1 
          call ccr12unpck2(work(ktr12),1,work(ktr12sq),'N',iopt)

C-----------------------------------------------------------------------
C         read X-integrals (R12 overlap matrix) from disk and reorder
C         to full square:
C-----------------------------------------------------------------------
          luxint = -1
          call gpopen(luxint,fccr12x,'old',' ','unformatted',idummy,
     &                .false.)
          rewind(luxint)
 9999     read(luxint) ian
          read(luxint) (work(kxmat+i), i=0, nr12r12p(1)-1 )
          if (ian.ne.ianr12) goto 9999
          call gpclose(luxint,'KEEP')
          iopt = 2
          call ccr12unpck2(work(kxmat),1,work(kxmatsq),'N',iopt)

C         if (locdbg) then
C           call  around('R12 overlap integrals in CC_XIETA1')
C           call cc_prsqr12(work(kxmatsq),1,'N',1,.false.)
C         end if

C--------------------------------------------------
C         read in VXINT and reorder to full square:
C--------------------------------------------------
          call dzero(work(kvxintsq),nr12r12sq(isyres))
          call cc_r12rdvxint(work(kvxintsq),work(kend3),lwrk3,one,
     &                       isyres,labelh)

C         if (locdbg) then
C           call around('VXINT in CC_XIETA1')
C           call cc_prsqr12(work(kvxintsq),isyres,'N',1,.false.)
C         end if

C--------------------------------------------------------------
C         read in V (perturbation operator) matrix in AO-basis:
C--------------------------------------------------------------
          call ccprpao(labelh,.TRUE.,work(kprpao),isymv,isym,
     &                 ierr,work(kend3),lwrk3) 
          IF ((IERR.GT.0).OR.((IERR.EQ.0).AND.(isymv.NE.isyres))) THEN
           CALL QUIT('CC_XIETA1: error while reading operator '//LABELH)
          ELSE IF (IERR.LT.0) THEN
           CALL DZERO(work(kprpao),N2BST(isyres))
          END IF 

C----------------------------------
C         calculate R12 part of Xi:
C----------------------------------
          call cc_r12xi(work(kxir12sq),isyres,'N',work(ktr12sq),1,
     &                  work(kxmatsq),work(kvxintsq),isyres,
     &                  work(kprpao),WORK(KLAMP0),WORK(KLAMH0),'N',
     &                  work(kend3),lwrk3)
          ! pack Xi to triangular format
          iopt = 1
          call ccr12pck2(work(kxir12),isyres,.false.,work(kxir12sq),
     &                   'N',iopt)
          call cclr_diasclr12(work(kxir12),brascl,isyres)

        END IF
        


        IF (LCAUCHY) THEN
          ! for the right Cauchy vectors "CO1" the situation is 
          ! slightly different: the right hand side vectors are then
          ! identical with the solution vectors with the Cauchy order
          ! decreased by 1 --> we don't need extra files for the right
          ! hand side vectors but use those for the solution vectors

          ! decremented Cauchy order
          ICAU = ILRCAU(IXETRAN(3,ITRAN))-1 

          ! list number for the vector with Cauchy order decreased by 1
          IFILE = ILRCAMP(LABELH,ICAU,ISYRES)


          ! initialize the SD part of the vector with C(n-1)
          IOPT = 3
          CALL CC_RDRSP(FILXI,IFILE,ISYRES,IOPT,MODEL,
     &                  WORK(KXI1),WORK(KXI2))
          CALL CCLR_DIASCL(WORK(KXI2),TWO,ISYRES)

        END IF


        IF (CCSDT .AND. (.NOT.SKIPXI)) THEN

          ! allocate extra memory for the 'effective' CC3 rhs vector
          KXI1EFF = KEND2
          KXI2EFF = KXI1EFF + NT1AM(ISYRES)
          KEND2   = KXI2EFF + NT2AM(ISYRES)
          LWRK2 = LWORK - KEND2
          IF (LWRK2 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CC_XIETA. (1c)')
          END IF

          ! find the associated frequency
          IF (IOPTRES.EQ.5 .OR. LCAUCHY) THEN
            FREQ  = 0.0D0
          ELSE
            IFILE = IXETRAN(3,ITRAN)
            FREQ  = FREQLST(FILXI,IFILE)
          END IF

C
          IF (IOPTRES.NE.5 .AND. (.NOT.NODDY_XI)) THEN
C           -----------
C           Open files:
C           -----------
C
            LUDELD  = -1
            LUCKJD  = -1
            LUDKBC  = -1
            LUTOC   = -1
            LU3VI   = -1
            LU3VI2  = -1
C
            CALL WOPEN2(LUDELD,FNDELD,64,0)
            CALL WOPEN2(LUCKJD,FNCKJD,64,0)
            CALL WOPEN2(LUDKBC,FNDKBC,64,0)
            CALL WOPEN2(LUTOC,FNTOC,64,0)
            CALL WOPEN2(LU3VI,FN3VI,64,0)
            CALL WOPEN2(LU3VI2,FN3VI2,64,0)
C
            CALL DZERO(WORK(KXI2EFF),NT2AM(ISYRES))

            CALL CC3_XISD(WORK(KXI1),WORK(KXI2),
     *                  WORK(KXI1EFF),WORK(KXI2EFF),
     *                  ISYRES,
     *                  WORK(KFOCKB(ITRAN)),ISYHOP,
     *                  FREQ,ICAU,LCAUCHY,LABELH,
     *                  WORK(KLAMP0),WORK(KLAMH0), 
     *                  WORK(KFOCK0MO),
     *                  LUCKJD,FNCKJD,LUDKBC,FNDKBC,LUDELD,FNDELD,
     *                  LUTOC,FNTOC,LU3VI,FN3VI,LU3VI2,FN3VI2,
     *                  WORK(KEND2),LWRK2)
 
            CALL WCLOSE2(LUDELD,FNDELD,'KEEP')
            CALL WCLOSE2(LUCKJD,FNCKJD,'KEEP')
            CALL WCLOSE2(LUTOC,FNTOC,'KEEP')
            CALL WCLOSE2(LUDKBC,FNDKBC,'KEEP')
            CALL WCLOSE2(LU3VI,FN3VI,'KEEP')
            CALL WCLOSE2(LU3VI2,FN3VI2,'KEEP')
C
          ELSE IF (NODDY_XI .AND. 
     &               (.NOT.(IOPTRES.EQ.5 .AND. NODDY_XI_ALTER))) THEN
C
            CALL CCSDT_XI_NODDY(IOPTRES,WORK(KXI1),WORK(KXI2),
     &                          WORK(KXI1EFF),WORK(KXI2EFF),
     &                          LABELH,WORK(KFOCKB(ITRAN)),
     &                          FREQ,ICAU,LCAUCHY,
     &                          WORK(KLAMP0),WORK(KLAMH0),
     &                          IXDOTS,XCONS,FILXI,ITRAN,
     &                          NXETRAN,MXVEC,WORK(KEND2),LWRK2)
          END IF
c
        END IF

*---------------------------------------------------------------------*
* transformation finished: scale diagonal of double excitation part
* and write result vector to file:
*---------------------------------------------------------------------*
        SKIPXI = ( IXETRAN(3,ITRAN) .EQ. -1 )

        IF (.NOT. SKIPXI) THEN

          IF (.NOT.CCS) THEN
            CALL CCLR_DIASCL(WORK(KXI2),HALF,ISYRES)
            IF (CCSDT) THEN
              CALL CCLR_DIASCL(WORK(KXI2EFF),HALF,ISYRES)
            END IF
          END IF

          IF (IOPTRES.EQ.0 .OR. IOPTRES.EQ.1) THEN
            IXETRAN(3,ITRAN) = IADRF_XI
            CALL PUTWA2(LUXI,FILXI,WORK(KXI1),IADRF_XI,NT1AM(ISYRES))
            IADRF_XI = IADRF_XI + NT1AM(ISYRES)
            IF (.NOT.CCS) THEN
              CALL PUTWA2(LUXI,FILXI,WORK(KXI2),IADRF_XI,NT2AM(ISYRES))
              IADRF_XI = IADRF_XI + NT2AM(ISYRES)
            END IF
            IF (CCR12) THEN
              CALL PUTWA2(LUXI,FILXI,WORK(KXIR12),IADRF_XI,
     &                    NTR12AM(ISYRES))
              IADRF_XI = IADRF_XI + NTR12AM(ISYRES)
            END IF
            IF (CCSDT) CALL QUIT('Problem in CC_XIETA')
          ELSE IF (IOPTRES.EQ.3) THEN
            IFILE  = IXETRAN(3,ITRAN)
            IF (LCAUCHY) IFILE = ILRCAMP(LABELH,ICAU,ISYRES)
            IF (ILSTSYM(FILXI,IFILE).NE.ISYRES) THEN
              CALL QUIT('Symmetry mismatch for Xi vector in CC_XIETA.')
            END IF
            IF (.NOT.LCAUCHY) THEN
              CALL CC_WRRSP(FILXI,IFILE,ISYRES,IOPTW,MODELW,DUMMY,
     &                      WORK(KXI1),WORK(KXI2),WORK(KEND2),LWRK2)
            END IF
            IF (CCR12) THEN
              CALL CC_WRRSP(FILXI,IFILE,ISYRES,IOPTWR12,MODELW,DUMMY,
     &                   DUMMY,WORK(KXIR12),WORK(KEND2),LWRK2)
            END IF
            IF (CCSDT) THEN
              CALL CC_WRRSP(FILXI,IFILE,ISYRES,IOPTWE,MODELW,DUMMY,
     &                   WORK(KXI1EFF),WORK(KXI2EFF),WORK(KEND2),LWRK2)
            END IF
          ELSE IF (IOPTRES.EQ.4) THEN
            IFILE  = IXETRAN(3,ITRAN)
            IF (LCAUCHY) IFILE = ILRCAMP(LABELH,ICAU,ISYRES)
            IF (ILSTSYM(FILXI,IFILE).NE.ISYRES) THEN
              CALL QUIT('Symmetry mismatch for Xi vector in CC_XIETA.')
            END IF
            IF (.NOT.LCAUCHY) THEN
              CALL CC_WARSP(FILXI,IFILE,ISYRES,IOPTW,MODELW,DUMMY,
     &                      WORK(KXI1),WORK(KXI2),WORK(KEND2),LWRK2)
            END IF
            IF (CCR12) THEN
              CALL CC_WARSP(FILXI,IFILE,ISYRES,IOPTWR12,MODELW,DUMMY,
     &                   DUMMY,WORK(KXIR12),WORK(KEND2),LWRK2)
            END IF
            IF (CCSDT) THEN
              CALL CC_WARSP(FILXI,IFILE,ISYRES,IOPTWE,MODELW,DUMMY,
     &                   WORK(KXI1EFF),WORK(KXI2EFF),WORK(KEND2),LWRK2)
            END IF
          ELSE IF (IOPTRES.EQ.5) THEN

            IF (LOCDBG.AND.CCSDT) THEN
             WRITE(LUPRI,*) 'XCONS TRIPLES CONTRIBUTION:'
             IVEC = 1
             DO WHILE (IXDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
              WRITE(LUPRI,*)'XCONS:',IVEC,ITRAN,XCONS(IVEC,ITRAN),IOPTW
              IVEC = IVEC + 1
             END DO
            END IF

            IF (.NOT.CCS) CALL CCLR_DIASCL(WORK(KXI2),TWO,ISYRES)
            CALL CCDOTRSP(IXDOTS,XCONS,IOPTW,FILXI,ITRAN,NXETRAN,MXVEC,
     &                    WORK(KXI1),WORK(KXI2),ISYRES,
     &                    WORK(KEND2),LWRK2)

            IF (CCR12) THEN
              CALL CCDOTRSP(IXDOTS,XCONS,IOPTWR12,FILXI,ITRAN,NXETRAN,
     &                      MXVEC,DUMMY,WORK(KXIR12),ISYRES,
     &                      WORK(KEND2),LWRK2)
            END IF  

            IF (LOCDBG) THEN
             IVEC = 1
             DO WHILE (IXDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
              WRITE(LUPRI,*)'XCONS:',IVEC,ITRAN,XCONS(IVEC,ITRAN),IOPTW
              IVEC = IVEC + 1
             END DO
            END IF
          ELSE
            CALL QUIT('Illegal value for IOPTRES in CC_XIETA.')
          END IF
     
          IF (LOCDBG) THEN
            WRITE (LUPRI,*) 'Final result of CC_XIETA:'
            WRITE (LUPRI,*) 'operator label:      ',LBLOPR(IOPER)
            WRITE (LUPRI,*) 'output file type:    ',FILXI(1:3)
            WRITE (LUPRI,*) 'index of output file:',IFILE
            WRITE (LUPRI,*) 'two-electron contr.: ',LTWOEL
            WRITE (LUPRI,*) 'relax/reorth. contr.:',LRELAX
            I = 1
            IF (CCS) I = 0
            CALL CC_PRP(WORK(KXI1),WORK(KXI2),ISYRES,1,I)
            IF (CCR12) THEN
              CALL CC_PRPR12(WORK(KXIR12),ISYRES,1,.TRUE.)
              WRITE(LUPRI,*) 'Norm^2 of Xi: ',
     &          DDOT(NT1AM(ISYRES),WORK(KXI1),1,WORK(KXI1),1) +
     &          DDOT(NT2AM(ISYRES),WORK(KXI2),1,WORK(KXI2),1) +
     &          DDOT(NTR12AM(ISYRES),WORK(KXIR12),1,WORK(KXIR12),1)
            END IF
            CALL FLSHFO(LUPRI)
            WRITE (LUPRI,*) 'WORK(0) = ',WORK(ITST)
          END IF

        END IF

*=====================================================================*
* calculate now the ETA vector:
*=====================================================================*
        IF (LOCDBG) THEN
          WRITE (LUPRI,*) 'SKIPETA:',SKIPETA
          WRITE (LUPRI,*) 'IDLSTL :',IDLSTL
          WRITE (LUPRI,*) 'IXETRAN(4,ITRAN):',IXETRAN(2,ITRAN)  
        END IF
 
        IF (.NOT. SKIPETA) THEN

          ISYRES = ISYETA ! symmetry of result vector

          KZETA1 = KEND1
          KETA1  = KZETA1 + NT1AM(ISYCTR)
          KEND2  = KETA1  + NT1AM(ISYRES)
          IF (CC2 .OR. CCSD .OR. CCSDT) THEN
             KETA2 = KEND2
             KEND2 = KETA2  + NT2AM(ISYRES)
          END IF
          LWRK2 = LWORK - KEND2

          IF (LWRK2 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CC_XIETA. (CCETA2)')
          END IF

          IOPT = 1
          CALL CC_RDRSP(LISTL,IDLSTL,ISYCTR,IOPT,MODEL,
     &                  WORK(KZETA1),DUMMY)

          IADRDI = 1
          IADRCI = 1
cch
c          write(lupri,*) 'kxint,kyint:',kxint(itran),kyint(itran)
cch

          CALL CCETA2(WORK(KETA1),WORK(KETA2),AVERAGE,
     &                WORK(KG0IM),WORK(KGBIM(ITRAN)),
     &                WORK(KR0IM),WORK(KRBIM(ITRAN)),
     &                WORK(KF0IM(ITRAN)),WORK(KFBIM(ITRAN)),
     &                WORK(KRZ0I(ITRAN)),WORK(KRZBI(ITRAN)),
     &                WORK(KFOCKB(ITRAN)),WORK(KONEHB(ITRAN)),
     &                WORK(KFCKHFB(ITRAN)),
     &                WORK(KZFCK0(ITRAN)),WORK(KZFCKB(ITRAN)),
     &                WORK(KGAMMA),WORK(KXINT(ITRAN)),
     &                WORK(KYINT(ITRAN)),WORK(KZETA1),
     &                LUBFIM,FNBFIM,IADRX1(1,ITRAN),
     &                LUBFZI,FNBFZI,IADRE0(1,ITRAN),IADRE1(1,ITRAN),
     &                LU0IAJB,FN0IAJB, IADRI0,
     &                LUDIM,FNDIM,IADRDI,
     &                LUCIM,FNCIM,IADRCI,
     &                LU0IAJB,FN0IAJB,IT2DEL0,
     &                LU1IAJB,FN1IAJB,IT2DELB(1,ITRAN),
     &                LUHINT0,FNHINT0,IADRH0(1,ITRAN),
     &                LUHINT1,FNHINT1,IADRH1(1,ITRAN),
     &                LUPQMO, FILPQMO,IADRPQMO(1,ITRAN),
     &                WORK(KLAMP0),WORK(KLAMH0),
     &                WORK(KLAMDPQ(ITRAN)),WORK(KLAMDHQ(ITRAN)),
     &                LISTL,IDLSTL,LABELH,ISYHOP,
     &                LTWOEL,LRELAX,WORK(KEND2), LWRK2)

C--------------------------------------
C         calculate R12 part of Eta{O}:
C--------------------------------------
          IF (CCR12 .AND. CCS) THEN
            CONTINUE
          ELSE IF (CCR12) THEN
            KETAR12 = KEND2
            KEND2   = KETAR12 + NTR12AM(ISYRES)
            LWRK2   = LWORK - KEND2

            KETAR12SQ = KEND2
            KTR12     = KETAR12SQ + NTR12SQ(ISYRES)
            KTR12SQ   = KTR12 + NTR12AM(1)
            KCTRR12   = KTR12SQ + NTR12SQ(1)
            KCTRR12SQ = KCTRR12 + NTR12AM(ISYCTR)
            KXMAT     = KCTRR12SQ + NTR12SQ(ISYCTR)
            KXMATSQ   = KXMAT + NR12R12P(1)
            KPRPAO    = KXMATSQ + NR12R12SQ(1)
            KEND3     = KPRPAO + N2BAST
            LWRK3   = LWORK - KEND3
            IF (LWRK3 .LT. 0) THEN
              CALL QUIT('Insufficient work space in CC_XIETA.(ETA-R12)')
            END IF

C--------------------------------------------------------------------
C         Read R12 amplitudes from disk and reorder to square:
C--------------------------------------------------------------------
            iopt=32
            call cc_rdrsp('R0 ',0,1,iopt,model,dummy,work(ktr12))
            iopt = 1
            call ccr12unpck2(work(ktr12),1,work(ktr12sq),'N',iopt)

C--------------------------------------------------------------------
C         Read R12 Lagr. multipliers from disk and reorder to square:
C--------------------------------------------------------------------
            iopt=32
            call cc_rdrsp(listl,idlstl,isyctr,iopt,model,dummy,
     &                    work(kctrr12))
            call cclr_diasclr12(work(kctrr12),2.0D0*brascl,isyctr)
            iopt = 1
            call ccr12unpck2(work(kctrr12),isyctr,work(kctrr12sq),
     &                       'N',iopt)

C-----------------------------------------------------------------------
C         read X-integrals (R12 overlap matrix) from disk and reorder
C         to full square:
C-----------------------------------------------------------------------
            luxint = -1
            call gpopen(luxint,fccr12x,'old',' ','unformatted',idummy,
     &                  .false.)
            rewind(luxint)
 7777       read(luxint) ian
            read(luxint) (work(kxmat+i), i=0, nr12r12p(1)-1 )
            if (ian.ne.ianr12) goto 7777
            call gpclose(luxint,'KEEP')
            iopt = 2
            call ccr12unpck2(work(kxmat),1,work(kxmatsq),'N',iopt)

C--------------------------------------------------------------
C         read in V (perturbation operator) matrix in AO-basis:
C--------------------------------------------------------------
            call ccprpao(labelh,.TRUE.,work(kprpao),isymv,isym,
     &                   ierr,work(kend3),lwrk3)
            IF (N2BST(ISYMV).GT.N2BAST) THEN
              CALL QUIT('Memory allocation error in CC_XIETA1')
            ELSE IF ((IERR.GT.0)) THEN
              CALL QUIT('CC_XIETA1: error while reading operator '//
     &                  LABELH)
            ELSE IF (IERR.LT.0) THEN
              CALL DZERO(work(kprpao),N2BST(isymv))
            END IF
 
C--------------------------------------------------
C         read in VXINT and reorder to full square:
C--------------------------------------------------
            KVXINTSQ = KEND3
            KEND3    = KVXINTSQ + NR12R12SQ(ISYMV)
            LWRK3   = LWORK - KEND3
            IF (LWRK3 .LT. 0) THEN
              CALL QUIT('Insufficient work space in CC_XIETA.(ETA-R12)')
            END IF
            call dzero(work(kvxintsq),nr12r12sq(isymv))
            call cc_r12rdvxint(work(kvxintsq),work(kend3),lwrk3,one,
     &                         isymv,labelh)

C--------------------------------------------------
C         calculate R12 contribution to Eta{O}_(ai):
C--------------------------------------------------
            call cc_r12etaa(work(keta1),isyres,work(kctrr12sq),isyctr,
     &                      work(ktr12sq),1,work(kxmatsq),work(kprpao),
     &                      isymv,WORK(KLAMP0),WORK(KLAMH0),.false.,
     &                      work(kend3),lwrk3)

C-------------------------------------
C         calculate R12 part of Eta{O}:
C-------------------------------------
            call cc_r12xi(work(ketar12sq),isyres,'N',work(kctrr12sq),
     &                    isyctr,work(kxmatsq),work(kvxintsq),isymv,
     &                    work(kprpao),WORK(KLAMP0),WORK(KLAMH0),'T',
     &                    work(kend3),lwrk3)

            ! pack Eta{O} to triangular format
            iopt = 1
            call ccr12pck2(work(ketar12),isyres,.false.,work(ketar12sq),
     &                     'N',iopt)
            call cclr_diasclr12(work(ketar12),0.5D0*ketscl,isyres)

            if (locdbg) then
              if (isyres.eq.1) then
                write(lupri,*) 'propr12 in CC_XIETA: ',
     &            ddot(ntr12am(1),work(ketar12),1,work(ktr12),1)
              else
                write(lupri,*) 'propr12 in CC_XIETA: zero by symmetry'
              end if
            end if 

          END IF
  
          IF (CCSDT) THEN
            ! allocate extra memory for the 'effective' CC3 rhs vector
            KETA1EFF = KEND2
            KETA2EFF = KETA1EFF + NT1AM(ISYRES)
            KEND2    = KETA2EFF + NT2AM(ISYRES)
            LWRK2    = LWORK - KEND2
            IF (LWRK2 .LT. 0) THEN
              CALL QUIT('Insufficient work space in CC_XIETA.(1c-eta)')
            END IF
C
            CALL DZERO(WORK(KETA1EFF),NT1AM(ISYRES))
            CALL DZERO(WORK(KETA2EFF),NT2AM(ISYRES))
C
            ! find the associated frequency
            IF (IOPTRES.EQ.5) THEN
              FREQ  = 0.0D0
            ELSE
              IFILE = IXETRAN(4,ITRAN)
              FREQ  = FREQLST(FILETA,IFILE)
            END IF
 
            IF (IOPTRES.EQ.5 .AND. (.NOT.NODDY_ETA)) THEN

              CONTINUE

            ELSE IF (.NOT.NODDY_ETA) THEN

C             -----------
C             Open files:
C             -----------
 
              LUDKBC   = -1
              LUCKJD   = -1
              LUTOC    = -1
              LU3VI    = -1
              LUDKBC3  = -1
              LU3FOPX  = -1
              LU3FOP2X = -1
 
              CALL WOPEN2(LUDKBC,FNDKBC,64,0)
              CALL WOPEN2(LUCKJD,FNCKJD,64,0)
              CALL WOPEN2(LUTOC,FNTOC,64,0)
              CALL WOPEN2(LU3VI,FN3VI,64,0)
              CALL WOPEN2(LUDKBC3,FNDKBC3,64,0)
              CALL WOPEN2(LU3FOPX,FN3FOPX,64,0)
              CALL WOPEN2(LU3FOP2X,FN3FOP2X,64,0)
 
              CALL CC3_ETASD(LISTL,IDLSTL,WORK(KFOCKB(ITRAN)),
     *                       ISYHOP,FREQ,WORK(KFOCK0MO),
     *                       WORK(KETA1),WORK(KETA2),WORK(KETA1EFF),
     *                       WORK(KETA2EFF),ISYRES,WORK(KEND2),LWRK2,
     *                       LUDKBC,FNDKBC,LUCKJD,FNCKJD,LUTOC,FNTOC,
     *                       LU3VI,FN3VI,LUDKBC3,FNDKBC3,
     *                       LU3FOPX,FN3FOPX,LU3FOP2X,FN3FOP2X)
 

C             ------------
C             Close files:
C             ------------
 
              CALL WCLOSE2(LUDKBC,FNDKBC,'KEEP')
              CALL WCLOSE2(LUCKJD,FNCKJD,'KEEP')
              CALL WCLOSE2(LUTOC,FNTOC,'KEEP')
              CALL WCLOSE2(LU3VI,FN3VI,'KEEP')
              CALL WCLOSE2(LUDKBC3,FNDKBC3,'KEEP')
              CALL WCLOSE2(LU3FOPX,FN3FOPX,'KEEP')
              CALL WCLOSE2(LU3FOP2X,FN3FOP2X,'KEEP')
 
            ELSE IF (NODDY_ETA .AND.
     &                (.NOT.(IOPTRES.EQ.5 .AND. NODDY_ETA_ALTER))) THEN

              CALL CCSDT_ETA_NODDY(LISTL,IDLSTL,IOPTRES,
     &                             WORK(KFOCKB(ITRAN)),WORK(KFOCKBAO),
     &                             FREQ,WORK(KFOCK0MO),
     &                             WORK(KETA1),WORK(KETA2),
     &                             WORK(KETA1EFF),WORK(KETA2EFF),
     &                             IEDOTS,ECONS,FILETA,ITRAN,
     &                             NXETRAN,MXVEC,WORK(KEND2),LWRK2 )

            END IF

          END IF


          IF (IOPTRES.EQ.0 .OR. IOPTRES.EQ.1) THEN
            IXETRAN(4,ITRAN) = IADRF_ETA
            CALL PUTWA2(LUETA,FILETA,WORK(KETA1),IADRF_ETA,
     &                  NT1AM(ISYRES))
            IADRF_ETA = IADRF_ETA + NT1AM(ISYRES)
            IF (.NOT.CCS) THEN
              CALL PUTWA2(LUETA,FILETA,WORK(KETA2),IADRF_ETA,
     &                    NT2AM(ISYRES))
              IADRF_ETA = IADRF_ETA + NT2AM(ISYRES)
            END IF
            IF (CCR12) THEN
              CALL PUTWA2(LUETA,FILETA,WORK(KETAR12),IADRF_ETA,
     &                    NTR12AM(ISYRES))
              IADRF_ETA = IADRF_ETA + NTR12AM(ISYRES)
            END IF
            IF (CCSDT) CALL QUIT('Problem in CC_XIETA')
          ELSE IF (IOPTRES.EQ.3) THEN
           IFILE  = IXETRAN(4,ITRAN)
           IF (ILSTSYM(FILETA,IFILE).NE.ISYRES) THEN
             CALL QUIT('Symmetry mismatch for Eta vector in CC_XIETA.')
           END IF
           CALL CC_WRRSP(FILETA,IFILE,ISYRES,IOPTW,MODELW,DUMMY,
     &                   WORK(KETA1),WORK(KETA2),WORK(KEND2),LWRK2)
           IF (CCR12 .AND. IOPTWR12.EQ.32) THEN
             CALL CC_WRRSP(FILETA,IFILE,ISYRES,IOPTWR12,MODELW,DUMMY,
     &                  DUMMY,WORK(KETAR12),WORK(KEND2),LWRK2)
           END IF
           IF (CCSDT) THEN
             CALL CC_WRRSP(FILETA,IFILE,ISYRES,IOPTWE,MODELW,DUMMY,
     &                 WORK(KETA1EFF),WORK(KETA2EFF),WORK(KEND2),LWRK2)
           END IF
          ELSE IF (IOPTRES.EQ.4) THEN
           IFILE  = IXETRAN(4,ITRAN)
           IF (ILSTSYM(FILETA,IFILE).NE.ISYRES) THEN
             CALL QUIT('Symmetry mismatch for Eta vector in CC_XIETA.')
           END IF
           CALL CC_WARSP(FILETA,IFILE,ISYRES,IOPTW,MODELW,DUMMY,
     &                   WORK(KETA1),WORK(KETA2),WORK(KEND2),LWRK2)
           IF (CCR12 .AND. IOPTWR12.EQ.32) THEN
             CALL CC_WARSP(FILETA,IFILE,ISYRES,IOPTWR12,MODELW,DUMMY,
     &                  DUMMY,WORK(KETAR12),WORK(KEND2),LWRK2)
           END IF
           IF (CCSDT) THEN
            CALL CC_WARSP(FILETA,IFILE,ISYRES,IOPTWE,MODELW,DUMMY,
     &                 WORK(KETA1EFF),WORK(KETA2EFF),WORK(KEND2),LWRK2)
           END IF
          ELSE IF (IOPTRES.EQ.5) THEN
           IF (LOCDBG) THEN
             IVEC = 1
             WRITE(LUPRI,*) 'ECONS TRIPLES CONTRIBUTION:'
             DO WHILE (IEDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
                WRITE (LUPRI,*) 
     &                'ECONS:',IVEC,ITRAN,ECONS(IVEC,ITRAN),IOPTW
                IVEC = IVEC + 1
             END DO
           END IF                
C
           IF (.NOT.CCS) CALL CCLR_DIASCL(WORK(KETA2),TWO,ISYRES)
           CALL CCDOTRSP(IEDOTS,ECONS,IOPTW,FILETA,ITRAN,NXETRAN,MXVEC,
     &                   WORK(KETA1),WORK(KETA2),ISYRES,
     &                   WORK(KEND2),LWRK2)
           IF (LOCDBG) THEN
             IVEC = 1
             DO WHILE (IEDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
                WRITE (LUPRI,*) 
     &                'ECONS:',IVEC,ITRAN,ECONS(IVEC,ITRAN),IOPTW
                IVEC = IVEC + 1
             END DO
           END IF
C
           IF (CCR12) THEN
             CALL CCDOTRSP(IEDOTS,ECONS,IOPTWR12,FILETA,ITRAN,NXETRAN,
     &                     MXVEC,DUMMY,WORK(KETAR12),ISYRES,
     &                     WORK(KEND2),LWRK2)
           END IF
C
           IF (LOCDBG) THEN
             write(lupri,*) 'ECONS R12 doubles contrib.:'
             IVEC = 1
             DO WHILE (IEDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
                WRITE (LUPRI,*) 
     &                'ECONS:',IVEC,ITRAN,ECONS(IVEC,ITRAN),IOPTWR12
                IVEC = IVEC + 1
             END DO
           END IF                
               
          ELSE
           CALL QUIT('Illegal value for IOPTRES in CC_XIETA.')
          END IF

          ! X1 vectors store <HF|J^(1)|CC> average on common block
          IF (FILETA(1:3).EQ.'X1 ') THEN
            AVEX1(IFILE) = AVERAGE
          END IF
     
          IF (LOCDBG) THEN
            WRITE (LUPRI,*) 'Final result of CC_XIETA:'
            WRITE (LUPRI,*) 'operator label:      ',LBLOPR(IOPER)
            WRITE (LUPRI,*) 'output file type:    ',FILETA
            WRITE (LUPRI,*) 'index of output file:',IFILE
            WRITE (LUPRI,*) 'two-electron contr.: ',LTWOEL
            WRITE (LUPRI,*) 'relax/reorth. contr.:',LRELAX
            I = 1
            IF (CCS) I = 0
            CALL CC_PRP(WORK(KETA1),WORK(KETA2),ISYRES,1,I)
            IF (CCR12) THEN
              CALL CC_PRPR12(WORK(KETAR12),ISYRES,I,.TRUE.)
              WRITE(LUPRI,*) 'Norm^2 of Eta: ',
     &          DDOT(NT1AM(ISYRES),WORK(KETA1),1,WORK(KETA1),1) +
     &          DDOT(NT2AM(ISYRES),WORK(KETA2),1,WORK(KETA2),1) +
     &          DDOT(NTR12AM(ISYRES),WORK(KETAR12),1,WORK(KETAR12),1)
            END IF
            IF (CCSDT.AND.IOPTRES.LT.5) THEN
              XNORM = DDOT(NT1AM(ISYRES),WORK(KETA1EFF),1,
     &                                   WORK(KETA1EFF),1)
     &              + DDOT(NT2AM(ISYRES),WORK(KETA2EFF),1,
     &                                   WORK(KETA2EFF),1)
              WRITE (LUPRI,*) 'effective Eta vector:',XNORM
              CALL CC_PRP(WORK(KETA1EFF),WORK(KETA2EFF),ISYRES,1,I)
            END IF
            CALL FLSHFO(LUPRI)
          END IF

        END IF
*=====================================================================*
* end of loop over transformations within this batch:
*=====================================================================*

      END DO ! ITRAN

*---------------------------------------------------------------------*
* delete the files with the BF and BFZ intermediates,
* close files for C and D intermediates
*---------------------------------------------------------------------*
         IF ( (LLTWOEL.OR.LLRELAX) .AND. (CCSD) )THEN
           CALL WCLOSE2(LUBFIM, FNBFIM,   'DELETE') 
C          IF (IBATCH.EQ.NBATCH) THEN
              CALL WCLOSE2(LUBFZI, FNBFZI,'DELETE') 
C          END IF
         END IF

         CALL WCLOSE2(LUCIM,  FNCIM,    'DELETE')
         CALL WCLOSE2(LUDIM,  FNDIM,    'DELETE')

*---------------------------------------------------------------------*
* close the loop over batches
*---------------------------------------------------------------------*
      END DO ! IBATCH 

*---------------------------------------------------------------------*
* close and delete integral and P, Q, C & D intermediate files:
*---------------------------------------------------------------------*
      CALL WCLOSE2(LU0IAJB, FN0IAJB, 'DELETE')
      CALL WCLOSE2(LU1IAJB, FN1IAJB, 'DELETE')
      IF (CC2) THEN
        CALL WCLOSE2(LU0AIBJ, FN0AIBJ, 'DELETE')
        CALL WCLOSE2(LU1AIBJ, FN1AIBJ, 'DELETE')
        CALL WCLOSE2(LUHINT0, FNHINT0, 'DELETE')
        CALL WCLOSE2(LUHINT1, FNHINT1, 'DELETE')
      ELSE
        CALL WCLOSE2(LU0IABJ, FN0IABJ, 'DELETE')
        CALL WCLOSE2(LU1IABJ, FN1IABJ, 'DELETE')
        CALL WCLOSE2(LU0IJBA, FN0IJBA, 'DELETE')
        CALL WCLOSE2(LU1IJBA, FN1IJBA, 'DELETE')
      END IF

      CALL WCLOSE2(LUPQMO, FILPQMO,  'DELETE')
      CALL WCLOSE2(LUPQ0,  FILPQ0,   'DELETE')
      CALL WCLOSE2(LUPQ1,  FILPQ1,   'DELETE')
     
      IF (IOPTRES.EQ.0) THEN
         CALL WCLOSE2(LUXI, FILXI, 'KEEP')
         CALL WCLOSE2(LUETA,FILETA,'KEEP')
      ELSE IF (IOPTRES.EQ.1) THEN
         CALL WCLOSE2(LUXI, FILXI, 'DELETE')
         CALL WCLOSE2(LUETA,FILETA,'DELETE')
         CALL QUIT('IOPTRES=1 not yet implemented in CC_XIETA.')
      END IF

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'CC_XIETA ended successfully (?) ... '//
     &        'return now.'
        CALL FLSHFO(LUPRI)
      END IF

      ! restore the DIRECT flag 
      DIRECT = DIRSAV 

      CALL FLSHFO(LUPRI)
      CALL QEXIT('CC_XIETA1')
      RETURN
      END
*=====================================================================*
*                  END OF SUBROUTINE CC_XIETA                         *
*=====================================================================*
*=====================================================================*
      SUBROUTINE CCXIINT1(ITRAN,   LABELH,  IORDER,  T2AMP, 
     &                    DENSPKQ, DPKHFQ,  FOCKQ,   FCKHFQ, 
     &                    XLAMDP,  XLAMDH,  ISYAMP,
     &                    XLAMDPQ, XLAMDHQ, DENSQ,   DNSHFQ,
     &                    ONEHQ,   ISYHOP,
     &                    FNBFDX0, LUBFDX0, IADRX0,  IADRBFX0,
     &                    FNBFDX1, LUBFDX1, IADRX1,  IADRBFX1,
     &                    LRELAX,  LTWOEL,  LZERO,   LNEWXI,
     &                    WORK,    LWORK)
*---------------------------------------------------------------------*
*
* Purpose:
*
*     Precalculate some intermediates for XI vector calculation
*     including the effective densities for the BF intermediate.
*
*     zeroth-order BF density only computed for LZERO and if
*        LTWOEL or LRELAX are set
*
*     derivative BF density only computed for LNEWXI and LRELAX
*
*     Fock and OneHam intermediates computed always
*
*     Christof Haettig 5-11-1998
*
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include "priunit.h"
#include "ccsdsym.h"
#include "maxorb.h"
#include "ccorb.h"
#include "ccfield.h"
#include "ccsdinp.h"

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      LOGICAL LRELAX, LTWOEL, LZERO, LNEWXI
      CHARACTER*(*) FNBFDX0, FNBFDX1
      CHARACTER*(8) LABELH
      INTEGER ITRAN, ISYAMP, ISYHOP, IORDER, LWORK
      INTEGER LUBFDX0, IADRBFX0, IADRX0(MXCORB_CC)
      INTEGER LUBFDX1, IADRBFX1, IADRX1(MXCORB_CC,*)
      
#if defined (SYS_CRAY)
      REAL DENSPKQ(*), DPKHFQ(*), FOCKQ(*), FCKHFQ(*)
      REAL XLAMDP(*), XLAMDH(*), XLAMDPQ(*), XLAMDHQ(*)
      REAL DENSQ(*), DNSHFQ(*), ONEHQ(*), T2AMP(*), WORK(*)
      REAL ZERO, THREE, DUMMY
#else
      DOUBLE PRECISION DENSPKQ(*), DPKHFQ(*), FOCKQ(*), FCKHFQ(*)
      DOUBLE PRECISION XLAMDP(*), XLAMDH(*), XLAMDPQ(*), XLAMDHQ(*)
      DOUBLE PRECISION DENSQ(*), DNSHFQ(*), ONEHQ(*), T2AMP(*), WORK(*)
      DOUBLE PRECISION ZERO, THREE, DUMMY
#endif
      PARAMETER (ZERO = 0.0D0, THREE = 3.0D0)

      CHARACTER MODEL*(10), LABEL1*(8), LABEL2*(8)
      INTEGER IOPT,IDEL,IDUMMY,IFIELD,IRREP,ISYM,IERR,ISIGNH,INUM

* external function:
      INTEGER IROPER2

      CALL QENTER('CCXIINT1')

*---------------------------------------------------------------------*
* generate lower triangular packed density matrix:
*---------------------------------------------------------------------*
      IF (LRELAX) THEN
        CALL CC_DNSPK(DENSQ,DENSPKQ,ISYHOP)
        IF (CC2) CALL CC_DNSPK(DNSHFQ,DPKHFQ,ISYHOP)
      END IF

*---------------------------------------------------------------------*
* get AO one-electron integrals h^X
*---------------------------------------------------------------------*
      IF ( LABELH(1:8) .EQ. 'HAM0    ' ) THEN

        CALL CCRHS_ONEAO(ONEHQ,WORK,LWORK)
*       for zeroth-order Hamiltonian add finite fields:
        DO IFIELD = 1, NFIELD
          CALL CC_ONEP(ONEHQ,WORK,LWORK,
     &                 EFIELD(IFIELD),ISYHOP,LFIELD(IFIELD)  ) 
        END DO

        ! scale the one-electron integrals with three:
        IF (LRELAX) THEN
           CALL DSCAL(N2BST(ISYHOP),THREE,ONEHQ,1)
           WRITE (LUPRI,*) 'Warning: multiply ONEHQ with 3 ...'
        END IF

      ELSE IF ( LABELH(1:8) .EQ. 'DUMMYOP ' ) THEN
        CALL DZERO(ONEHQ,N2BST(ISYHOP))
      ELSE IF ( IORDER.EQ.2 .AND. LABELH(1:2) .EQ. '->' ) THEN
        ! `second-order operator' with zero AO integrals
        CALL DZERO(ONEHQ,N2BST(ISYHOP))
      ELSE 

        CALL CCPRPAO(LABELH,.TRUE.,ONEHQ,IRREP,ISYM,IERR,WORK,LWORK)
        IF ((IERR.GT.0) .OR. (IERR.EQ.0 .AND. IRREP.NE.ISYHOP)) THEN
          CALL QUIT('CCXIINT1: error while reading operator '//LABELH)
        ELSE IF (IERR.LT.0) THEN
          CALL DZERO(ONEHQ,N2BST(ISYHOP))
        END IF

        ! check if second-order operator and in case fix the sign
        IF (IORDER.EQ.2) THEN
          LABEL1 = '?'
          LABEL2 = '?'
          INUM = IROPER2(LABEL1,LABEL2,LABELH,ISIGNH,ISYM)
          IF (ISIGNH.NE.1) THEN
            WRITE (LUPRI,*) 'CCXIINT1>', LABELH, ISIGNH, DBLE(ISIGNH)
            CALL DSCAL(N2BST(ISYHOP),DBLE(ISIGNH),ONEHQ,1)
          END IF
        END IF

      END IF

*---------------------------------------------------------------------*
* initialize derivative AO Fock matrix with h^x integrals,
*---------------------------------------------------------------------*
      CALL DCOPY(N2BST(ISYHOP),ONEHQ,1,FOCKQ,1)
      IF (CC2) CALL DCOPY(N2BST(ISYHOP),ONEHQ,1,FCKHFQ,1)

*---------------------------------------------------------------------*
* calculate effective density matrices for the BF and G intermediates:
*---------------------------------------------------------------------*
      IF (.NOT. CCS) THEN

         IF (LZERO .AND. (LRELAX.OR.LTWOEL) ) THEN
           IOPT = 1
           CALL CC_BFDEN(T2AMP,  ISYAMP, DUMMY,  IDUMMY,
     *                   XLAMDH, ISYAMP, XLAMDH, ISYAMP,
     *                   XLAMDH, ISYAMP, DUMMY,  IDUMMY,
     *                   FNBFDX0,LUBFDX0,IADRX0, IADRBFX0,
     *                   1, IOPT, .FALSE., WORK, LWORK)
         END IF

         IF (LNEWXI .AND. LRELAX) THEN
            IOPT = 5
            CALL CC_BFDEN(T2AMP,  ISYAMP, DUMMY,   IDUMMY,
     *                    XLAMDHQ,ISYHOP, XLAMDH,  ISYAMP,
     *                    XLAMDH, ISYAMP, XLAMDHQ, ISYHOP,
     *                    FNBFDX1,LUBFDX1,IADRX1, IADRBFX1,
     *                    ITRAN, IOPT, .FALSE., WORK, LWORK)
         ELSE IF (LRELAX) THEN
            DO IDEL = 1, NBAST
              IADRX1(IDEL,ITRAN) = IADRX1(IDEL,ITRAN-1)
            END DO
         ELSE
            DO IDEL = 1, NBAST
              IADRX1(IDEL,ITRAN) = -999999
            END DO
         END IF

      END IF

*---------------------------------------------------------------------*
* that's it; return:
*---------------------------------------------------------------------*
      CALL QEXIT('CCXIINT1')

      RETURN

      END 
*=====================================================================*
*=====================================================================*
      SUBROUTINE CCETAINT1(ITRAN,   ISTART,  LISTL,   IDLSTL,
     &                     ZDPK0,   ZDEN0,   ZDPKB,   ZDENB,
     &                     XINT,    YINT,    MINT,    CHI, CHIQ,
     &                     ZETA1,   ZETA2,   ISYCTR,
     &                     T2AMP,   XLAMDP,  XLAMDH,  ISYAMP,
     &                     XLAMDPQ, XLAMDHQ, ISYHOP,
     &                     FNBFDE0, LUBFDE0, IADRE0,  IADRBFE0,
     &                     FNBFDE1, LUBFDE1, IADRE1,  IADRBFE1,
     &                     FILPQMO, LUPQMO,  IADRPQMO,IADRPQ,
     &                     FILPQ0,  LUPQ0,   IADRPQ0, IADRPQI0, 
     &                     FILPQ1,  LUPQ1,   IADRPQ1, IADRPQI1,
     &                     LRELAX,  LTWOEL,  WORK,    LWORK)
*---------------------------------------------------------------------*
*
* Purpose:
*
*     Precalculate some intermediates for Eta vector calculation:
*
*     CCSD: the X and Y intermediates, the effective density for the
*           BFZeta intermediate and the P and Q intermediates
*
*     CC2:  the X and Y intermediates and the two-index backtransf.
*           Zeta2 vector
*
*     Christof Haettig 5-11-1998, CC2 added in spring 2000
*
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include "priunit.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "maxorb.h"
#include "ccorb.h"

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)

      LOGICAL LRELAX, LTWOEL
      CHARACTER*(*) LISTL, FNBFDE0, FNBFDE1, FILPQMO, FILPQ0, FILPQ1
      INTEGER ITRAN, ISTART, IDLSTL, ISYAMP, ISYCTR, ISYHOP, LWORK
      INTEGER LUBFDE0, IADRBFE0, IADRE0(MXCORB_CC,*)
      INTEGER LUBFDE1, IADRBFE1, IADRE1(MXCORB_CC,*)
      INTEGER LUPQMO,  IADRPQ,   IADRPQMO(MXCORB_CC,*)
      INTEGER LUPQ0,   IADRPQI0, IADRPQ0(MXCORB_CC,*)
      INTEGER LUPQ1,   IADRPQI1, IADRPQ1(MXCORB_CC,*)
      
#if defined (SYS_CRAY)
      REAL XINT(*), YINT(*), MINT(*), ZETA1(*), ZETA2(*)
      REAL ZDPK0(*), ZDEN0(*), ZDPKB(*), ZDENB(*), CHI(*), CHIQ(*)
      REAL XLAMDP(*), XLAMDH(*), XLAMDPQ(*), XLAMDHQ(*)
      REAL T2AMP(*), WORK(*)
      REAL DUMMY
#else
      DOUBLE PRECISION XINT(*), YINT(*), MINT(*), ZETA1(*), ZETA2(*)
      DOUBLE PRECISION ZDPK0(*), ZDEN0(*), ZDPKB(*), ZDENB(*)
      DOUBLE PRECISION CHI(*), CHIQ(*)
      DOUBLE PRECISION XLAMDP(*), XLAMDH(*), XLAMDPQ(*), XLAMDHQ(*)
      DOUBLE PRECISION T2AMP(*), WORK(*)
      DOUBLE PRECISION DUMMY
#endif

      CHARACTER MODEL*(10)
      INTEGER IOPT, ISYINT, ISYINTQ, KCHI, KCHIQ, KEND1, LWRK1, IDEL
      INTEGER IDUMMY, N2VEC
      INTEGER IDLSTOLD
      LOGICAL DONE_BFDEN, DONE_PQIM
      SAVE IDLSTOLD, DONE_BFDEN, DONE_PQIM

      CALL QENTER('CCETAINT1')

*---------------------------------------------------------------------*
* do some tests and set symmetries:
*---------------------------------------------------------------------*

* nothing to do for unrelaxed CCS:
      IF (CCS .AND. (.NOT.LRELAX)) THEN
         CALL QEXIT('CCETAINT1')
         RETURN
      ENDIF

* if ITRAN=ISTART, make sure that everything will be initialized:
CCH
         IDLSTOLD   = IDLSTL - 1
CCH
      IF (ITRAN.EQ.ISTART) THEN
         IDLSTOLD   = IDLSTL - 1
         DONE_BFDEN = .FALSE.
         DONE_PQIM  = .FALSE.
      END IF

* set symmetries for intermediates:
      ISYINT  = MULD2H(ISYCTR,ISYAMP)
      ISYINTQ = MULD2H(ISYINT,ISYHOP)

* flag for double excitation vector
      N2VEC = 1
      IF (CCS) N2VEC = 0

*---------------------------------------------------------------------*
* if left vector (IDLSTL) changed, read new multipliers into memory:
*---------------------------------------------------------------------*
      IF (IDLSTL .NE. IDLSTOLD) THEN

         DONE_BFDEN = .FALSE.
         DONE_PQIM  = .FALSE.

         IF (CCS) THEN
           IOPT = 1
           CALL CC_RDRSP(LISTL,IDLSTL,ISYCTR,IOPT,MODEL,ZETA1,WORK)
         ELSE
           IF (LWORK .LT. NT2AM(ISYCTR)) THEN
              CALL QUIT('Insufficient memroy in CCETAINT1 (a)')
           END IF 

           IOPT = 3
           CALL CC_RDRSP(LISTL,IDLSTL,ISYCTR,IOPT,MODEL,ZETA1,WORK)

           CALL CC_T2SQ(WORK,ZETA2,ISYCTR)
         END IF

         IF (LOCDBG) THEN
            WRITE (LUPRI,*) 'CCETAINT1> the zeta2 vector:'
            CALL CC_PRP(WORK,WORK,ISYCTR,1,N2VEC)
            WRITE (LUPRI,*) 'CCETAINT1> the t2amp vector:'
            CALL CC_PRP(WORK,T2AMP,ISYAMP,1,N2VEC)
         END IF

      END IF

*---------------------------------------------------------------------*
* for relaxed CC2 calculate Zeta densities:
*---------------------------------------------------------------------*
      IF ( (CCS .OR. CC2) .AND. LRELAX ) THEN

         IOPT = 1
         CALL DZERO(ZDEN0,N2BST(ISYINT))
         CALL CC_ZETADEN(ZDEN0,ZETA1,ISYCTR,XLAMDP,XLAMDH,ISYM0,
     &                   XLAMDP,XLAMDH,ISYM0,IOPT,WORK,LWORK)
         CALL CC_DNSPK(ZDEN0,ZDPK0,ISYINT)

         IOPT = 2
         CALL DZERO(ZDENB,N2BST(ISYINTQ))
         CALL CC_ZETADEN(ZDENB,ZETA1,ISYCTR,XLAMDP,XLAMDH,ISYM0,
     &                   XLAMDPQ,XLAMDHQ,ISYHOP,IOPT,WORK,LWORK)
         CALL CC_DNSPK(ZDENB,ZDPKB,ISYINTQ)

      END IF

      ! nothing more to do for CCS:
      IF (CCS) THEN
        IDLSTOLD = IDLSTL
        RETURN
      END IF

*---------------------------------------------------------------------*
* if left vector changed, calculate new X and Y intermediates:
*---------------------------------------------------------------------*
      IF ( IDLSTL .NE. IDLSTOLD ) THEN

          CALL CC_XI(XINT,ZETA2,ISYCTR,T2AMP,ISYAMP,WORK,LWORK)  

          CALL CC_YI(YINT,ZETA2,ISYCTR,T2AMP,ISYAMP,WORK,LWORK)  

          IF (LOCDBG) THEN
             WRITE (LUPRI,'(//A)') 'CCETAINT1> X-intermediate:'
             WRITE (LUPRI,'(5G15.6)') (XINT(I),I=1,NMATIJ(ISYINT))
             WRITE (LUPRI,'(//A)') 'CCETAINT1> Y-intermediate:'
             WRITE (LUPRI,'(5G15.6)') (YINT(I),I=1,NMATAB(ISYINT))
          END IF

      END IF

*---------------------------------------------------------------------*
* calculate Chi matrices (backtransformed Zeta1 + X intermediate):
*     --  KCHI  : backtransf. with XLAMDP
*     --  KCHIQ : backtransf. with XLAMDPQ
*---------------------------------------------------------------------*
      IF (CCSD .AND. LRELAX) THEN
        KCHI  = 1
        KCHIQ = KCHI  + NGLMDT(ISYINT)
        KEND1 = KCHIQ + NGLMDT(ISYINTQ)
        LWRK1 = LWORK - KEND1
        
        IF (LWORK .LT. 0) THEN
          CALL QUIT('Insufficient memroy in CCETAINT1 (b)')
        END IF 
  
        CALL CCLT_CHI(ZETA1,XINT,ISYINT,XLAMDP, ISYAMP,WORK(KCHI), 1)
       
        CALL CCLT_CHI(ZETA1,XINT,ISYINT,XLAMDPQ,ISYHOP,WORK(KCHIQ),1)

      ELSE IF (CC2 .AND. LRELAX) THEN
        KEND1 = 1
        LWRK1 = LWORK - KEND1

        CALL CCLT_CHI(ZETA1,DUMMY,ISYINT,XLAMDP, ISYAMP,CHI, 0)
       
        CALL CCLT_CHI(ZETA1,DUMMY,ISYINT,XLAMDPQ,ISYHOP,CHIQ,0)

      END IF

*---------------------------------------------------------------------*
* calculate the effective BFZ density which does not depend on the
* external perturbation and the M intermediate only when the left 
* vector has changed and they have not been calculate before for in 
* this batch. For CC2 we calculate here just the two-index backtransf.
* Zeta2 vector, i.e. the contributions from Chi and M are excluded.
*---------------------------------------------------------------------*
      IF (LRELAX .AND. (.NOT.DONE_BFDEN) ) THEN

        IF (CC2) THEN
          IOPT = 0
        ELSE
          CALL CC_MI(MINT,ZETA2,ISYCTR,T2AMP,ISYAMP,WORK(KEND1),LWRK1)  
          IOPT = 2
        END IF

        CALL CC_BFDEN(ZETA2,     ISYCTR, MINT,   ISYINT,
     *                XLAMDP,    ISYAMP, XLAMDP, ISYAMP,
     *                WORK(KCHI),ISYINT, DUMMY,  IDUMMY,
     *                FNBFDE0, LUBFDE0, IADRE0, IADRBFE0,
     *                ITRAN, IOPT, .FALSE., WORK(KEND1), LWRK1)

        DONE_BFDEN = .TRUE.

      ELSE IF (ITRAN.GT.1) THEN
        DO IDEL = 1, NBAST
           IADRE0(IDEL,ITRAN) = IADRE0(IDEL,ITRAN-1)
        END DO
      END IF

*---------------------------------------------------------------------*
* calculate the effective BFZ density which does depend on the
* external perturbation for each transformation. For CC2 the contribs.
* from ChiQ and Mint are skipped.
*---------------------------------------------------------------------*
      IF (LRELAX) THEN
        IOPT = 6
        IF (CC2) IOPT = 9
        CALL CC_BFDEN(ZETA2,     ISYCTR, MINT,        ISYINT,
     *                XLAMDPQ,   ISYHOP, XLAMDP,      ISYAMP,
     *                WORK(KCHI),ISYINT, WORK(KCHIQ), ISYINTQ,
     *                FNBFDE1,   LUBFDE1, IADRE1, IADRBFE1,
     *                ITRAN, IOPT, .FALSE., WORK(KEND1), LWRK1)
      END IF

*---------------------------------------------------------------------*
* calculate one-index backtransformed P and Q intermediates used 
* in CC_21I for the G term contribution to the Eta vector:
*---------------------------------------------------------------------*
      IF (CCSD .AND. LRELAX) THEN

        IF ( .NOT. DONE_PQIM ) THEN

          IOPT = 1
          CALL CC_PQIMO(ZETA2,ISYCTR,T2AMP,ISYAMP,YINT,IOPT,
     *                  FILPQMO,LUPQMO,IADRPQMO,IADRPQ,
     *                  ITRAN,WORK(KEND1),LWRK1)

          CALL CC_PQIAO(FILPQMO,LUPQMO,IADRPQMO,ISYINT,
     *                  FILPQ0, LUPQ0, IADRPQ0, IADRPQI0,
     *                  ITRAN, XLAMDH, ISYAMP,WORK(KEND1), LWRK1)

          DONE_PQIM = .TRUE.

        ELSE IF (ITRAN.GT.1) THEN
          DO IDEL = 1, NBAST
            IADRPQMO(IDEL,ITRAN) = IADRPQMO(IDEL,ITRAN-1)
            IADRPQ0(IDEL,ITRAN)  = IADRPQ0(IDEL,ITRAN-1)
          END DO
        END IF


        CALL CC_PQIAO(FILPQMO,LUPQMO,IADRPQMO,ISYINT,
     *                FILPQ1, LUPQ1, IADRPQ1, IADRPQI1,
     *                ITRAN,XLAMDHQ,ISYHOP,WORK(KEND1), LWRK1)

      ELSE IF (ITRAN.GT.1) THEN
        DO IDEL = 1, NBAST
          IADRPQMO(IDEL,ITRAN) = IADRPQMO(IDEL,ITRAN-1)
          IADRPQ0(IDEL,ITRAN)  = IADRPQ0(IDEL,ITRAN-1)
          IADRPQ1(IDEL,ITRAN)  = -1
        END DO
      END IF
*---------------------------------------------------------------------*
* save present IDLSTL in IDLSTOLD and return:
*---------------------------------------------------------------------*
      IDLSTOLD = IDLSTL

      CALL QEXIT('CCETAINT1')

      RETURN

      END 
*=====================================================================*
*=====================================================================*
      SUBROUTINE CCETAINT2(X0INT, D0SRHF, ISY0DIS, 
     &                     X1INT, D1SRHF, ISY1DIS, SQRINT,
     &                     BFZ0, BFZI, F0IM, FBIM, ZFCK0, ZFCKB,
     &                     XLAMDP0,XLAMDH0,XLAMDPQ,XLAMDHQ,CHI,CHIQ,
     &                     ZDPK0,   ZDEN0,   ZDPKB,   ZDENB, 
     &                     IREAL,   ISYHOP,  ISYCTR,
     &                     FNBFDE0, LUBFDE0, IADRE0,  
     &                     FNBFDE1, LUBFDE1, IADRE1,  
     &                     FILPQ0,  LUPQ0,   IADRPQ0,  
     &                     FILPQ1,  LUPQ1,   IADRPQ1,  
     &                     FNHINT0, LUHINT0, IADRH0, IADRZ0,
     &                     FNHINT1, LUHINT1, IADRH1, IADRZ1,
     &                     LRELAX,  LTWOEL,  LZERO,
     &                     IDEL,    WORK,    LWORK)
*---------------------------------------------------------------------*
*
* Purpose:
*
*     Precalculate some intermediates for Eta vector calculation
*
*     deriv. contrib. depending on X1INT are only computed for LTWOEL
*     deriv. contrib. depending on X0INT are only computed for LRELAX
*     zeroth-order intermediates are only computed for LZERO
*
*     Christof Haettig 5-11-1998, CC2 added in spring 2000
*
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include "priunit.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "maxorb.h"
#include "ccorb.h"
#include "ccisao.h"

      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)

      LOGICAL LRELAX, LTWOEL, LZERO, SQRINT
      CHARACTER*(*) FNBFDE0, FNBFDE1, FILPQ0, FILPQ1, FNHINT0, FNHINT1
      INTEGER ISYCTR,ISYAMP,ISYHOP,ISY0DIS,ISY1DIS,IDEL,IREAL,LWORK
      INTEGER IADRZ0, IADRZ1
      INTEGER LUBFDE0, IADRE0(MXCORB_CC)
      INTEGER LUBFDE1, IADRE1(MXCORB_CC)
      INTEGER LUPQ0, IADRPQ0(MXCORB_CC)
      INTEGER LUPQ1, IADRPQ1(MXCORB_CC)
      INTEGER LUHINT0, IADRH0(MXCORB_CC)
      INTEGER LUHINT1, IADRH1(MXCORB_CC)
      
#if defined (SYS_CRAY)
      REAL XLAMDP0(*), XLAMDH0(*), XLAMDPQ(*), XLAMDHQ(*)
      REAL BFZ0(*), BFZI(*), F0IM(*), FBIM(*), CHI(*), CHIQ(*)
      REAL ZFCK0(*), ZFCKB(*), ZDPK0(*), ZDEN0(*), ZDPKB(*), ZDENB(*)
      REAL X0INT(*), X1INT(*), D0SRHF(*), D1SRHF(*)
      REAL WORK(LWORK), DUMMY, SIGN, ONE, DDOT, XNORM
#else
      DOUBLE PRECISION XLAMDP0(*), XLAMDH0(*), XLAMDPQ(*), XLAMDHQ(*)
      DOUBLE PRECISION BFZ0(*), BFZI(*), F0IM(*), FBIM(*)
      DOUBLE PRECISION ZFCK0(*), ZFCKB(*), CHI(*), CHIQ(*)
      DOUBLE PRECISION ZDPK0(*), ZDEN0(*), ZDPKB(*), ZDENB(*)
      DOUBLE PRECISION X0INT(*), X1INT(*), D0SRHF(*), D1SRHF(*)
      DOUBLE PRECISION WORK(LWORK), DUMMY, SIGN, ONE, DDOT, XNORM
#endif
      PARAMETER (ONE = 1.0D0)
     
      INTEGER ISYDEL, ISYZT1, ISYZT2, NMGD, IADR, KMGD, KEND1, LWRK1
      INTEGER ISYM, LEN1, LEN0, KPINT0, KQINT0, KPINT1, KQINT1, IOPT
      INTEGER ISYETA, IDUMMY, NSRHF, KDSRHF, KX0IAJB, KX1IAJB
      INTEGER KEND2, LWRK2, ISY0IAJ, ISY1IAJ, ISY0ALBE, ISY1ALBE
      INTEGER JGAM, ISYGAM, KOFF0, KOFF1

      CALL QENTER('CCETAINT2')

*---------------------------------------------------------------------*
*     set some symmetries:
*---------------------------------------------------------------------*
      ISYDEL = ISAO(IDEL)
      D      = IDEL - IBAS(ISYDEL)

      ISYETA = MULD2H(ISYCTR,ISYHOP)
      ISYZT1 = MULD2H(ISYDEL,ISYCTR)
      ISYZT2 = MULD2H(ISYDEL,ISYETA)

*---------------------------------------------------------------------*
*     for relaxed CCS/CC2 calculate the Zeta-Fock matrices:
*---------------------------------------------------------------------*
      IF ( (CC2 .OR. CCS) .AND. LRELAX) THEN
   
        CALL CC_AOFOCK2(X0INT,ZDEN0,ZDPK0,ZFCK0,WORK,LWORK,
     &                  IDEL,ISY0DIS,ISYDEL,ISYCTR,.FALSE.)

        CALL CC_AOFOCK2(X0INT,ZDENB,ZDPKB,ZFCKB,WORK,LWORK,
     &                  IDEL,ISY0DIS,ISYDEL,ISYETA,.FALSE.)

        IF (LTWOEL) THEN
          CALL CC_AOFOCK2(X1INT,ZDEN0,ZDPK0,ZFCKB,WORK,LWORK,
     &                    IDEL,ISY1DIS,ISYDEL,ISYCTR,SQRINT)
        END IF

      END IF

      ! this was everything for CCS
      IF (CCS) RETURN

*---------------------------------------------------------------------*
*     for CCSD calculate the BZeta intermediates:
*     for CC2 calculate the GZeta intermediates:
*---------------------------------------------------------------------*
      IF (.NOT.CCS) THEN

*        ----------------------------------------------
*        allocate an array for the effective densities:
*        ----------------------------------------------
         NMGD = 0
         DO ISYM = 1, NSYM
           NMGD = MAX(NMGD,NT2BGD(ISYM))
         END DO

         NSRHF = NDSRHF(ISY1DIS)
         IF (SQRINT) NSRHF = NDSRHFSQ(ISY1DIS)
       
         KMGD  = 1
         KEND1 = KMGD  + NMGD
         IF (CC2) THEN
           KDSRHF = KEND1
           KEND1  = KDSRHF + NSRHF
         END IF
         LWRK1 = LWORK - KEND1

         IF (LWRK1 .LT. 0) THEN
           CALL QUIT('Insufficient work space in CC_XIETA. (4)')
         END IF


         IF (LTWOEL .OR. LRELAX .OR. LZERO) THEN
*           -------------------------------------------------------
*           read zeroth-order BFZ effective density 
*           -------------------------------------------------------
            IADR = IADRE0(IDEL)
            NMGD = NT2BGD(ISYZT1)
            CALL GETWA2(LUBFDE0,FNBFDE0,WORK(KMGD),IADR,NMGD)

            IF (CC2) THEN
*             -------------------------------------------------------
*             for CC2 the first contributions to the GZeta intermeds.
*             which are stored in the F0IM/FBIM arrays
*             this uses complex conjug. integrals --> sign for LAO's
*             -------------------------------------------------------
              CALL CCTRBT2(X0INT,WORK(KDSRHF),XLAMDHQ,ISYHOP,
     &               WORK(KEND1),LWRK1,ISY0DIS,0,.FALSE.,SQRINT,ONE)

              IF (LTWOEL) THEN 
                SIGN = DBLE(IREAL)
                CALL DAXPY(NSRHF,SIGN,D1SRHF,1,WORK(KDSRHF),1)
              END IF

              IOPT = 0
              CALL CC_GIM1(WORK(KDSRHF),ISY1DIS,WORK(KMGD),ISYZT1,
     &                     FBIM,IOPT,SQRINT,WORK(KEND1),LWRK1)

              IF (LZERO) THEN
                IOPT = 0
                CALL CC_GIM1(D0SRHF,ISY0DIS,WORK(KMGD),ISYZT1,
     &                       F0IM,IOPT,.FALSE.,WORK(KEND1),LWRK1)
              END IF 

            ELSE
*             -------------------------------------------------------
*             for CCSD the first contribution to the BFZ intermeds.
*             Note: the calculation of the BFZ intermediate uses
*                   complex conjugated integrals --> sign for LAO's
*             -------------------------------------------------------
              IF (LTWOEL) THEN
                SIGN = DBLE(IREAL)
                CALL CC_BFI(BFZI,X1INT,ISY1DIS,WORK(KMGD),ISYZT1,
     &                      D,ISYDEL,SQRINT,SIGN,WORK(KEND1),LWRK1)
              END IF

              IF (LZERO) THEN
                 CALL CC_BFI(BFZ0,X0INT,ISY0DIS,WORK(KMGD),ISYZT1,
     &                       D,ISYDEL,.FALSE.,ONE,WORK(KEND1),LWRK1)
              END IF

            END IF

         END IF


         IF (LRELAX) THEN
*           -------------------------------------------------------
*           read response contribution to BFZ effective density
*           -------------------------------------------------------
            IADR = IADRE1(IDEL)
            NMGD = NT2BGD(ISYZT2)
            CALL GETWA2(LUBFDE1,FNBFDE1,WORK(KMGD),IADR,NMGD)

            IF (CC2) THEN
*               ----------------------------------------------------
*               calculate 2. contribution to the GZeta intermediate
*               which is stored in the FBIM array
*               ----------------------------------------------------
                IOPT = 0
                CALL CC_GIM1(D0SRHF,ISY0DIS,WORK(KMGD),ISYZT2,
     &                       FBIM,IOPT,.FALSE.,WORK(KEND1),LWRK1)
            ELSE
*               ----------------------------------------------------
*               calculate 2. contribution to the BFZ intermediate
*               ----------------------------------------------------
                CALL CC_BFI(BFZI,X0INT,ISY0DIS,WORK(KMGD),ISYZT2,
     &                      D,ISYDEL,.FALSE.,ONE,WORK(KEND1),LWRK1)
            END IF

         END IF


      END IF

*---------------------------------------------------------------------*
*     calculate the 21F term contributions for Eta
*---------------------------------------------------------------------*
      IF (.NOT.CC2) THEN 
         
         LEN0   = NT2BCD(ISYZT1)
         LEN1   = NT2BCD(ISYZT2)

         KPINT0 = 1
         KQINT0 = KPINT0 + LEN0
         KPINT1 = KQINT0 + LEN0
         KQINT1 = KPINT1 + LEN1
         KEND1  = KQINT1 + LEN1
         LWRK1  = LWORK  - KEND1

         IF (LWRK1 .LT. 0) THEN
           CALL QUIT('Insufficient work space in CC_XIETA. (4c)')
         END IF

         IADR = IADRPQ0(IDEL)
         CALL GETWA2(LUPQ0,FILPQ0,WORK(KPINT0),IADR,LEN0)
         IADR = IADRPQ0(IDEL) + LEN0
         CALL GETWA2(LUPQ0,FILPQ0,WORK(KQINT0),IADR,LEN0)

         IADR = IADRPQ1(IDEL)
         CALL GETWA2(LUPQ1,FILPQ1,WORK(KPINT1),IADR,LEN1)
         IADR = IADRPQ1(IDEL) + LEN1
         CALL GETWA2(LUPQ1,FILPQ1,WORK(KQINT1),IADR,LEN1)
 
         IF (LZERO) THEN
            IOPT = 1
            CALL CC_21I2(F0IM,X0INT,ISY0DIS,DUMMY,IDUMMY,
     &                   WORK(KPINT0),WORK(KQINT0),ISYZT1,
     &                   DUMMY,       DUMMY,       IDUMMY,
     &                   XLAMDP0,XLAMDH0,ISYM0,XLAMDP0,ISYM0,
     &                   WORK(KEND1),LWRK1,IOPT,
     &                   .TRUE.,.FALSE.,.FALSE.)
         END IF

         IOPT = 3
         CALL CC_21I2(FBIM,X0INT,ISY0DIS,X1INT,ISY1DIS,
     &                WORK(KPINT0),WORK(KQINT0),ISYZT1,
     &                WORK(KPINT1),WORK(KQINT1),ISYZT2,
     &                XLAMDP0,XLAMDH0,ISYM0,XLAMDPQ,ISYHOP,
     &                WORK(KEND1),LWRK1,IOPT,
     &                .TRUE.,LTWOEL,SQRINT)
             
      END IF

*---------------------------------------------------------------------*
*     for CC2 compute modified integrals needed for neta^H contrib.
*---------------------------------------------------------------------*
      IF (CC2 .AND. LRELAX) THEN

         ISY0IAJ = MULD2H(ISY0DIS,ISYCTR)
         ISY1IAJ = MULD2H(ISY0DIS,ISYETA)

C        ---------------------------------------
C        allocate memory for modified integrals:
C        ---------------------------------------
         KX0IAJB = KEND1 
         KX1IAJB = KX0IAJB + NT2BCD(ISY0IAJ)
         KEND2   = KX1IAJB + NT2BCD(ISY1IAJ)
         LWRK2   = LWORK   - KEND2

         IF (LWRK2 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CCETAINT2. (iajb)')
         END IF

         CALL DZERO(WORK(KX0IAJB),NT2BCD(ISY0IAJ))
         CALL DZERO(WORK(KX1IAJB),NT2BCD(ISY1IAJ))

C        ---------------------------------------------------
C        do the 3-index transformation in a loop over gamma:
C        ---------------------------------------------------
         DO ISYGAM = 1, NSYM
           
           ISY0ALBE = MULD2H(ISY0DIS,ISYGAM)
           ISY1ALBE = MULD2H(ISY1DIS,ISYGAM)

           DO G = 1, NBAS(ISYGAM)
             JGAM = G + IBAS(ISYGAM)

             KOFF0  =  IDSAOG(ISYGAM,ISY0DIS)  +NNBST(ISY0ALBE)*(G-1)+1
             IF (SQRINT) THEN
               KOFF1 = IDSAOGSQ(ISYGAM,ISY1DIS)+N2BST(ISY1ALBE)*(G-1)+1
             ELSE
               KOFF1 = IDSAOG(ISYGAM,ISY1DIS)  +NNBST(ISY1ALBE)*(G-1)+1
             END IF

             IOPT = 0
             CALL CC_IAJB(X0INT(KOFF0),ISY0ALBE,X1INT(KOFF1),ISY1ALBE,
     &                    IDEL,JGAM,.FALSE.,IDUMMY,
     &                    WORK(KX0IAJB),DUMMY,DUMMY,
     &                    WORK(KX1IAJB),DUMMY,DUMMY,
     &                    XLAMDP0,XLAMDH0,ISYM0,XLAMDPQ,XLAMDHQ,ISYHOP,
     &                    CHI, DUMMY, ISYCTR, CHIQ, DUMMY, ISYETA,
     &                    WORK(KEND2),LWRK2,IOPT,LTWOEL,LRELAX,LZERO,
     &                    .FALSE.,SQRINT,IREAL)
           END DO 
         END DO

C        --------------------------------------------
C        write 3-index transformed integrals to disk:
C        --------------------------------------------
         LEN0 = NT2BCD(ISY0IAJ)
         CALL PUTWA2(LUHINT0,FNHINT0,WORK(KX0IAJB),IADRZ0,LEN0)
         IADRH0(IDEL) = IADRZ0
         IADRZ0       = IADRZ0 + LEN0

         LEN1 = NT2BCD(ISY1IAJ)
         CALL PUTWA2(LUHINT1,FNHINT1,WORK(KX1IAJB),IADRZ1,LEN1)
         IADRH1(IDEL) = IADRZ1
         IADRZ1       = IADRZ1 + LEN1

      END IF

*----------------------------------------------------------------------*
*     That's it; return:
*----------------------------------------------------------------------*

      CALL QEXIT('CCETAINT2')

      RETURN
      END 
*======================================================================*
*                   END OF SUBROUTINE CCETAINT2                        *
*======================================================================*
*======================================================================*
c /* deck ccxiintao */
*======================================================================*
      SUBROUTINE CCXIINTAO(X0INT, D0SRHF, ISY0DIS,
     &                     X1INT, D1SRHF, ISY1DIS, SQRINT,
     &                     XLAMDP0,XLAMDH0,XLAMDQP,XLAMDQH,
     &                     DENS0,DPCK0,DNSHF0,DPKHF0,
     &                     DENSQ,DPCKQ,DNSHFQ,DPKHFQ,
     &                     FOCK0,FOCKB,FCKHFB, 
     &                     RHO2,G0IM,GBIM,
     &                     LUBFDX0,FNBFDX0,IADRX0,
     &                     LUBFDX1,FNBFDX1,IADRX1,
     &                     LU0IAJB,LU0IABJ,LU0IJBA,LU0AIBJ,
     &                     FN0IAJB,FN0IABJ,FN0IJBA,FN0AIBJ,
     &                     LU1IAJB,LU1IABJ,LU1IJBA,LU1AIBJ,
     &                     FN1IAJB,FN1IABJ,FN1IJBA,FN1AIBJ,
     &                     IT2DEL0,IADR0,IT2DELB,IADRB,
     &                     ITRAN,IDEL,LZERO,LNEWXI,LRELAX,LTWOEL,
     &                     IREAL,ISYHOP,WORK, LWORK)
*----------------------------------------------------------------------*
*
*   Purpose: calculate intermediates for Xi right hand side vector
*            which depend on the AO integrals
*
*     contrib. depending on X1INT/D1SRHF are only computed for LTWOEL
*     contrib. depending on X0INT/D0SRHF are only computed for LRELAX
*     (only exception: the (ia|jb), (ij|ba), (ia|bj) integrals)
*     zeroth-order intermediates are only computed for LZERO
*
*   Written by Christof Haettig, November 1998
*   LAO's, CC2 and symmetry in winter 1999/2000
*
*----------------------------------------------------------------------*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccorb.h"
#include "maxorb.h"
#include "ccsdsym.h"
#include "ccsdinp.h"
#include "ccisao.h"

      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

* variables:
      LOGICAL LZERO, LNEWXI, LRELAX, LTWOEL, SQRINT
      CHARACTER*(*) FNBFDX0, FNBFDX1
      CHARACTER*(*) FN0IAJB,FN0IABJ,FN0IJBA,FN0AIBJ
      CHARACTER*(*) FN1IAJB,FN1IABJ,FN1IJBA,FN1AIBJ
      INTEGER ITRAN, IDEL, IREAL, ISYHOP, LWORK
      INTEGER ISY0DIS, ISY1DIS, LUBFDX0, LUBFDX1
      INTEGER LU0IAJB, LU0IABJ, LU0IJBA, LU0AIBJ
      INTEGER LU1IAJB, LU1IABJ, LU1IJBA, LU1AIBJ
      INTEGER IADRX0(*), IADRX1(*)
      INTEGER IT2DEL0(MXCORB_CC)
      INTEGER IT2DELB(MXCORB_CC,*)

#if defined (SYS_CRAY)
      REAL X0INT(*), X1INT(*), D0SRHF(*), D1SRHF(*)
      REAL XLAMDP0(*), XLAMDH0(*), XLAMDQP(*), XLAMDQH(*)
      REAL DENS0(*), DPCK0(*), DNSHF0(*), DPKHF0(*)
      REAL DENSQ(*), DPCKQ(*), DNSHFQ(*), DPKHFQ(*)
      REAL RHO2(*), G0IM(*), GBIM(*)
      REAL FOCK0(*), FOCKB(*), FCKHFB(*)
      REAL WORK(LWORK)
      REAL ONE, ZERO, TWO, XNORM, DDOT, DUMMY
#else
      DOUBLE PRECISION X0INT(*), X1INT(*), D0SRHF(*), D1SRHF(*)
      DOUBLE PRECISION XLAMDP0(*), XLAMDH0(*), XLAMDQP(*), XLAMDQH(*)
      DOUBLE PRECISION DENS0(*), DPCK0(*), DNSHF0(*), DPKHF0(*)
      DOUBLE PRECISION DENSQ(*), DPCKQ(*), DNSHFQ(*), DPKHFQ(*)
      DOUBLE PRECISION RHO2(*), G0IM(*), GBIM(*)
      DOUBLE PRECISION FOCK0(*), FOCKB(*), FCKHFB(*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION ONE, ZERO, TWO, XNORM, DDOT, DUMMY
#endif
      PARAMETER (ONE = 1.0d0, ZERO = 0.0d0, TWO = 2.0d0)

      LOGICAL LLAO
      INTEGER ISYDEL, ISYMM1, ISYMM2, NMGD, KEND4, LWRK4, IADR, KMGD
      INTEGER KSCRCM2, KSCRCM1, ISY0IAJ, ISY1IAJ, KDSRHF, NSRHF
      INTEGER KX1IAJB, KX1IABJ, KX1IJBA, KX0AIBJ
      INTEGER KX0IAJB, KX0IABJ, KX0IJBA, KX1AIBJ
      INTEGER LEN0, LEN1, IADR0, IADRB, ISYGAM, ISY0ALBE, ISY1ALBE
      INTEGER JGAM, KOFF0, KOFF1, ISYSRH1, KEND5, LWRK5, IOPT, ISYM
      INTEGER IDUMMY

      CALL QENTER('CCXIINTAO')

*---------------------------------------------------------------------*
*     begin:
*---------------------------------------------------------------------*
      IF (LOCDBG) WRITE (LUPRI,*) 'entered from CCXIINTAO...'

* nothing to do for unrelaxed CCS
      IF (CCS .AND. (.NOT.LRELAX)) RETURN

      ISYDEL = ISAO(IDEL)
      D      = IDEL - IBAS(ISYDEL)

*---------------------------------------------------------------------*
*     add contribution to Fock matrices Fbar and F^0:
*---------------------------------------------------------------------*

C     IF (LNEWXI .AND. LRELAX) THEN
      IF (LRELAX) THEN
         CALL CC_AOFOCK2(X0INT,DENSQ,DPCKQ,FOCKB,WORK,LWORK,
     &                   IDEL,ISY0DIS,ISYDEL,ISYHOP,.FALSE.)

         IF ((CCS .OR. CC2) .AND. LZERO) THEN
            CALL CC_AOFOCK2(X0INT,DENS0,DPCK0,FOCK0,WORK,LWORK,
     &                      IDEL,ISY0DIS,ISYDEL,ISYM0,.FALSE.)
         END IF

      END IF

C     IF (LNEWXI .AND. LTWOEL) THEN
      IF (LTWOEL) THEN
         CALL CC_AOFOCK2(X1INT,DENS0,DPCK0,FOCKB,WORK,LWORK,
     &                   IDEL,ISY1DIS,ISYDEL,ISYM0,SQRINT)
      END IF


      IF (CC2 .AND. LRELAX) THEN
         CALL CC_AOFOCK2(X0INT,DNSHFQ,DPKHFQ,FCKHFB,WORK,LWORK,
     &                   IDEL,ISY0DIS,ISYDEL,ISYHOP,.FALSE.)

         IF (LTWOEL) THEN
            CALL CC_AOFOCK2(X1INT,DNSHF0,DPKHF0,FCKHFB,WORK,LWORK,
     &                      IDEL,ISY1DIS,ISYDEL,ISYM0,SQRINT)
         END IF
      END IF

*---------------------------------------------------------------------*
*     for CCS we are done ...
*---------------------------------------------------------------------*
      IF (CCS) RETURN

*---------------------------------------------------------------------*
*     for CCSD calculate the first-order BF intermediate,
*     for CC2 the zero- and first-order G intermediates
*---------------------------------------------------------------------*
      ISYMM1 = MULD2H(ISYDEL,ISYM0)
      ISYMM2 = MULD2H(ISYDEL,ISYHOP)

*     --------------------------------------------------------
*     allocate an array for the different effective densities:
*     --------------------------------------------------------
      NMGD = 0
      DO ISYM = 1, NSYM
         NMGD = MAX(NMGD,NT2BGD(ISYM))
      END DO
       
      NSRHF = NDSRHF(ISY1DIS)  
      IF (SQRINT) NSRHF = NDSRHFSQ(ISY1DIS) 

      KMGD  = 1
      KEND4 = KMGD  + NMGD
      IF (CC2) THEN
        KDSRHF = KEND4
        KEND4  = KDSRHF + NSRHF
      END IF
      LWRK4 = LWORK - KEND4

      IF (LWRK4 .LT. 0) THEN
        CALL QUIT('Insufficient work space in CCXIINTAO. (4)')
      END IF

      IF (LNEWXI .AND. (LTWOEL .OR. LRELAX)) THEN
*        --------------------------------------
*        read zeroth-order BF effective density
*        --------------------------------------
         IADR = IADRX0(IDEL)
         NMGD = NT2BGD(ISYMM1)
         CALL GETWA2(LUBFDX0,FNBFDX0,WORK(KMGD),IADR,NMGD)

         IF (CC2) THEN
*           -----------------------------------------------------
*           CC2 contribution: the G intermediates...
*           -----------------------------------------------------
            CALL CCTRBT2(X0INT,WORK(KDSRHF),XLAMDQP,ISYHOP,
     &              WORK(KEND4),LWRK4,ISY0DIS,0,.FALSE.,SQRINT,ONE)

            IF (LTWOEL) CALL DAXPY(NSRHF,ONE,D1SRHF,1,WORK(KDSRHF),1)

            IOPT = 1
            CALL CC_GIM1(WORK(KDSRHF),ISY1DIS,WORK(KMGD),ISYMM1, 
     &                   GBIM,1,SQRINT,WORK(KEND4),LWRK4)

            IF (LZERO) THEN ! calculate zero-order G intermediate... 
               ! note that BF-density has already been transformed
               ! to 2Cou-Exc form --> IOPT = 0
               IOPT = 0
               CALL CC_GIM1(D0SRHF,ISY0DIS,WORK(KMGD),ISYMM1,G0IM,IOPT,
     &                      .FALSE.,WORK(KEND4),LWRK4)
            END IF

         ELSE 
*           -----------------------------------------------------
*           CCSD contribution: the BF intermediate...
*           -----------------------------------------------------
            IF (LTWOEL) THEN
              CALL CC_BFI(RHO2,X1INT,ISY1DIS,WORK(KMGD),ISYMM1,
     &                    D,ISYDEL,SQRINT,ONE,WORK(KEND4),LWRK4)
            END IF

         END IF

      END IF

      IF (LNEWXI .AND. LRELAX) THEN
*        --------------------------------------------------------
*        read response contribution to BF effective density
*        --------------------------------------------------------
         IADR = IADRX1(IDEL)
         NMGD = NT2BGD(ISYMM2)
         CALL GETWA2(LUBFDX1,FNBFDX1,WORK(KMGD),IADR,NMGD)

         IF (CC2) THEN
*           -----------------------------------------------------
*           CC2 contribution: the G intermediate...
*           -----------------------------------------------------
            IOPT = 1
            CALL CC_GIM1(D0SRHF,ISY0DIS,WORK(KMGD),ISYMM2,GBIM,IOPT,
     &                   .FALSE.,WORK(KEND4),LWRK4)

         ELSE
*           -----------------------------------------------------
*           CCSD contribution: the BF intermediate...
*           -----------------------------------------------------
            CALL CC_BFI(RHO2,X0INT,ISY0DIS,WORK(KMGD),ISYMM2,
     &                  D,ISYDEL,.FALSE.,ONE,WORK(KEND4),LWRK4)
         END IF

      END IF

*---------------------------------------------------------------------*
*    calculate 3-index transformed integrals: 
*          (ia|j del), (ia|j del)-bar  (for CCSD & CC2)
*          (ia|del j), (ia|del j)-bar  (for CCSD only)
*          (ij|del a), (ij|del a)-bar  (for CCSD only)
*          (ai|del j), (ai|del j)-bar  (for CC2  only)
*---------------------------------------------------------------------*

      ISY0IAJ = MULD2H(ISY0DIS,ISYM0)
      ISY1IAJ = MULD2H(ISY0DIS,ISYHOP)

C     -------------------------------------
C     allocate memory for integral batches:
C     -------------------------------------
      KEND4 = 1

      IF (LZERO) THEN
         KX0IAJB = KEND4
         KEND4   = KX0IAJB + NT2BCD(ISY0IAJ)
         IF (CC2) THEN
           KX0AIBJ = KEND4
           KEND4   = KX0AIBJ + NT2BCD(ISY0IAJ)
         ELSE
           KX0IABJ = KEND4
           KX0IJBA = KX0IABJ + NT2BCD(ISY0IAJ)
           KEND4   = KX0IJBA + NT2BCD(ISY0IAJ)
         END IF
      END IF

      IF (LNEWXI) THEN
         KX1IAJB = KEND4
         KEND4   = KX1IAJB + NT2BCD(ISY1IAJ)
         IF (CC2) THEN
           KX1AIBJ = KEND4
           KEND4   = KX1AIBJ + NT2BCD(ISY0IAJ)
         ELSE
           KX1IABJ = KEND4
           KX1IJBA = KX1IABJ + NT2BCD(ISY1IAJ)
           KEND4   = KX1IJBA + NT2BCD(ISY1IAJ)
         END IF
      END IF

      LWRK4   = LWORK - KEND4

      IF (LWRK4 .LT. 0) THEN
         CALL QUIT('Insufficient work space in CCXIINTAO. (4b)')
      END IF

      IF (LZERO) THEN
         CALL DZERO(WORK(KX0IAJB),NT2BCD(ISY0IAJ))
         IF (.NOT.CC2) CALL DZERO(WORK(KX0IABJ),NT2BCD(ISY0IAJ))
         IF (.NOT.CC2) CALL DZERO(WORK(KX0IJBA),NT2BCD(ISY0IAJ))
         IF (CC2)      CALL DZERO(WORK(KX0AIBJ),NT2BCD(ISY0IAJ))
      END IF

      IF (LNEWXI) THEN
         CALL DZERO(WORK(KX1IAJB),NT2BCD(ISY1IAJ))
         IF (.NOT.CC2) CALL DZERO(WORK(KX1IABJ),NT2BCD(ISY1IAJ))
         IF (.NOT.CC2) CALL DZERO(WORK(KX1IJBA),NT2BCD(ISY1IAJ))
         IF (CC2)      CALL DZERO(WORK(KX1AIBJ),NT2BCD(ISY1IAJ))
      END IF

C     ---------------------------------------------------
C     do the 3-index transformation in a loop over gamma:
C     ---------------------------------------------------
      IF (LNEWXI) THEN
        DO ISYGAM = 1, NSYM

          ISY0ALBE = MULD2H(ISY0DIS,ISYGAM)
          ISY1ALBE = MULD2H(ISY1DIS,ISYGAM)

          DO G = 1, NBAS(ISYGAM)
            JGAM = G + IBAS(ISYGAM)
                   
            KOFF0  =  IDSAOG(ISYGAM,ISY0DIS)  +NNBST(ISY0ALBE)*(G-1)+1

            IF (SQRINT) THEN
              KOFF1 = IDSAOGSQ(ISYGAM,ISY1DIS)+N2BST(ISY1ALBE)*(G-1)+1
            ELSE
              KOFF1 = IDSAOG(ISYGAM,ISY1DIS)  +NNBST(ISY1ALBE)*(G-1)+1
            END IF


C           ------------------------------------------------------
C           CCSD: compute (i a|j del), (i a|del j), and (ij|del a) 
C           CC2 : compute only (i a|j del) integrals
C           ------------------------------------------------------
            IOPT = 2
            IF (CC2) IOPT = 0

            CALL CC_IAJB(X0INT(KOFF0),ISY0ALBE,X1INT(KOFF1),ISY1ALBE,
     &                   IDEL,JGAM,.FALSE.,IDUMMY,
     &                   WORK(KX0IAJB),WORK(KX0IABJ),WORK(KX0IJBA),
     &                   WORK(KX1IAJB),WORK(KX1IABJ),WORK(KX1IJBA),
     &                   XLAMDP0,XLAMDH0,ISYM0,XLAMDQP,XLAMDQH,ISYHOP,
     &                   XLAMDP0,XLAMDH0,ISYM0,XLAMDQP,XLAMDQH,ISYHOP,
     &                   WORK(KEND4),LWRK4,IOPT,LTWOEL,LRELAX,LZERO,
     &                   .FALSE.,SQRINT,IREAL)

            
C           ---------------------------------------------------------
C           for CC2 compute (a i|del j) integrals as (i a|j del) with
C           particle and hole Lambda matrices interchanged --> for
C           imaginary perturbations we must invert the sign of the 
C           AO derivative integrals before the transformation
C           ---------------------------------------------------------
            IF (CC2) THEN
              IF (IREAL .EQ. -1) THEN
                 IF (.NOT.SQRINT) 
     *              CALL QUIT('Illegal IREAL/SQRINT combin.')
                 CALL DSCAL(N2BST(ISY1ALBE),-ONE,X1INT(KOFF1),1)
              END IF

              IOPT = 0
              CALL CC_IAJB(X0INT(KOFF0),ISY0ALBE,X1INT(KOFF1),ISY1ALBE,
     &                     IDEL,JGAM,.FALSE.,IDUMMY,
     &                     WORK(KX0AIBJ),DUMMY,DUMMY,
     &                     WORK(KX1AIBJ),DUMMY,DUMMY,
     &                     XLAMDH0,XLAMDP0,ISYM0,XLAMDQH,XLAMDQP,ISYHOP,
     &                     XLAMDH0,XLAMDP0,ISYM0,XLAMDQH,XLAMDQP,ISYHOP,
     &                     WORK(KEND4),LWRK4,IOPT,LTWOEL,LRELAX,LZERO,
     &                     .FALSE.,SQRINT,IREAL)

              IF (IREAL .EQ. -1) THEN
                 CALL DSCAL(N2BST(ISY1ALBE),-ONE,X1INT(KOFF1),1)
              END IF
            END IF

          END DO

        END DO
      END IF

C     ------------------------------------
C     transform (ia|del j) to L(ia|del j):
C     ------------------------------------
      IF (LZERO  .AND. (.NOT.CC2)) THEN
        CALL DSCAL(NT2BCD(ISY0IAJ), TWO,WORK(KX0IABJ),1)
        CALL DAXPY(NT2BCD(ISY0IAJ),-ONE,WORK(KX0IJBA),1,
     *             WORK(KX0IABJ),1)
      END IF

      IF (LNEWXI .AND. (.NOT.CC2)) THEN
        CALL DSCAL(NT2BCD(ISY1IAJ), TWO,WORK(KX1IABJ),1)
        CALL DAXPY(NT2BCD(ISY1IAJ),-ONE,WORK(KX1IJBA),1,
     *             WORK(KX1IABJ),1)
      END IF


C     --------------------------------------------
C     write 3-index transformed integrals to disk:
C     --------------------------------------------
      IF (LZERO) THEN
         LEN0 = NT2BCD(ISY0IAJ)

         CALL PUTWA2(LU0IAJB, FN0IAJB, WORK(KX0IAJB), IADR0, LEN0)
         IF (.NOT.CC2) THEN
           CALL PUTWA2(LU0IABJ, FN0IABJ, WORK(KX0IABJ), IADR0, LEN0)
           CALL PUTWA2(LU0IJBA, FN0IJBA, WORK(KX0IJBA), IADR0, LEN0)
         ELSE
           CALL PUTWA2(LU0AIBJ, FN0AIBJ, WORK(KX0AIBJ), IADR0, LEN0)
         END  IF

         IT2DEL0(IDEL) = IADR0
         IADR0 = IADR0 + LEN0
      END IF


      IF (LNEWXI) THEN
         LEN1 = NT2BCD(ISY1IAJ)

         CALL PUTWA2(LU1IAJB, FN1IAJB, WORK(KX1IAJB), IADRB, LEN1)
         IF (.NOT.CC2) THEN
           CALL PUTWA2(LU1IABJ, FN1IABJ, WORK(KX1IABJ), IADRB, LEN1)
           CALL PUTWA2(LU1IJBA, FN1IJBA, WORK(KX1IJBA), IADRB, LEN1)
         ELSE
           CALL PUTWA2(LU1AIBJ, FN1AIBJ, WORK(KX1AIBJ), IADRB, LEN1)
         END  IF

         IT2DELB(IDEL,ITRAN) = IADRB
         IADRB = IADRB + LEN1
      ELSE
         IT2DELB(IDEL,ITRAN) = IT2DELB(IDEL,ITRAN-1)
      ENDIF

      IF (IT2DELB(IDEL,ITRAN).LE.0) THEN
        WRITE (LUPRI,*) 
     *        'IDEL,ITRAN,IT2DELB:',IDEL,ITRAN,IT2DELB(IDEL,ITRAN)
        CALL QUIT('Fatal error in CCXIINTAO!')
      END IF

*---------------------------------------------------------------------*
* that's it;  return:
*---------------------------------------------------------------------*
      IF (LOCDBG) WRITE (LUPRI,*) 'return from CCXIINTAO...'

      CALL QEXIT('CCXIINTAO')

      RETURN
      END 
*=====================================================================*
*                  END OF SUBROUTINE CCXIINTAO                        *
*=====================================================================*
*=====================================================================*
      SUBROUTINE CCXI2(RHO1,    RHO2,    IORDER,  SKIPXI,
     &                 GAMMA,   EMAT1,   EMAT2,
     &                 G0IM,    GBIM,    R0IM,    RBIM, 
     &                 XLAMDP0, XLAMDH0, CMO0,  
     &                 ONEH0,   FOCK0,   
     &                 XLAMDQP, XLAMDQH, CMOPQ,   CMOHQ, 
     &                 ONEHB,   FOCKB,   FCKHFB,  ISYHOP,
     &                 LABEL1,  IRELAX1, ISYOP1,
     &                 LABEL2,  IRELAX2, ISYOP2,
     &                 FNBFDX1, LUBFDX1, IADRX1,
     &                 FN0IAJB, FN0IJBA, FN0IABJ, FN0AIBJ,
     &                 LU0IAJB, LU0IJBA, LU0IABJ, LU0AIBJ,
     &                 FN1IAJB, FN1IJBA, FN1IABJ, FN1AIBJ,
     &                 LU1IAJB, LU1IJBA, LU1IABJ, LU1AIBJ,
     &                 LUDIM,FNDIM,IADRDI,
     &                 LUCIM,FNCIM,IADRCI,
     &                 IT2DEL0, IT2DELB, IADRINT,
     &                 LRELAX,  LTWOEL,  
     &                 WORK,    LWORK)
*---------------------------------------------------------------------*
*
* Purpose: calculate Xi vector from intermediates
*          (if SKIPXI=TRUE, calculate only intermediates)
*
*     Christof Haettig 26-11-1998
*
*---------------------------------------------------------------------*
       IMPLICIT NONE
#include "priunit.h"
#include "dummy.h"
#include "ccsdsym.h"
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccfield.h"


      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER LUBF0, LUFCK

      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)

      CHARACTER*(*) FNBFDX1, FNCIM,   FNDIM
      CHARACTER*(*) FN0IAJB, FN0IJBA, FN0IABJ, FN0AIBJ
      CHARACTER*(*) FN1IAJB, FN1IJBA, FN1IABJ, FN1AIBJ
      CHARACTER*(8) LABEL1, LABEL2
      LOGICAL LRELAX, LTWOEL, SKIPXI
      INTEGER LWORK, ISYHOP, LUBFDX1, LUCIM, LUDIM, IORDER
      INTEGER LU0IAJB, LU0IJBA, LU0IABJ, LU0AIBJ
      INTEGER LU1IAJB, LU1IJBA, LU1IABJ, LU1AIBJ
      INTEGER IT2DEL0(*), IT2DELB(*), IADRX1, IADRINT, IADRCI, IADRDI
      INTEGER ISYOP1, ISYOP2, IRELAX1, IRELAX2

#if defined (SYS_CRAY)
      REAL XLAMDP0(*), XLAMDH0(*), XLAMDQP(*), XLAMDQH(*)
      REAL CMO0(*), CMOPQ(*), CMOHQ(*), ONEHB(*), ONEH0(*)
      REAL FOCK0(*), FOCKB(*), FCKHFB(*)
      REAL GAMMA(*), G0IM(*), GBIM(*)
      REAL RHO1(*), RHO2(*), EMAT1(*), EMAT2(*)
      REAL R0IM(*), RBIM(*), WORK(LWORK)
      REAL HALF, ONE, FAC, XNORM, DDOT, DNRM2
#else
      DOUBLE PRECISION XLAMDP0(*), XLAMDH0(*), XLAMDQP(*), XLAMDQH(*)
      DOUBLE PRECISION CMO0(*), CMOPQ(*), CMOHQ(*), ONEHB(*), ONEH0(*)
      DOUBLE PRECISION FOCK0(*), FOCKB(*), FCKHFB(*)
      DOUBLE PRECISION GAMMA(*), G0IM(*), GBIM(*)
      DOUBLE PRECISION RHO1(*), RHO2(*), EMAT1(*), EMAT2(*)
      DOUBLE PRECISION R0IM(*), RBIM(*), WORK(LWORK)
      DOUBLE PRECISION HALF, ONE, FAC, XNORM, DDOT, DNRM2
#endif
      PARAMETER (HALF = 0.5D0, ONE = 1.0D0)

      CHARACTER CDUM, MODEL*(10)
      LOGICAL LGAMMA, LO3BF, LBFZETA, LRCON, LGCON, FCKCON, LOCC
      LOGICAL LRELAX1, LRELAX2, LTRSPF
      INTEGER KEND1, LWRK1, KEND2, LWRK2, LENBF0, LENBF1, KBF, IOPTG
      INTEGER ICON, ISYMBF, IOPT, KT2AMP0, ISYIAJB, KXIAJB, KT2AM
      INTEGER KCBAR, KDBAR, IDUM, KEND3, LWRK3, ISYRES, IF, IOPTR12
      INTEGER KONEH1, KONEH2, KCMOP1, KCMOP2, KLAMDP1, KLAMDP2
      INTEGER KCMOH1, KCMOH2, KLAMDH1, KLAMDH2, KT1AMP0, KFCKHF0
      INTEGER IRREP, IERR, ISYM, IOPER1, IOPER2, KXAIBJ, KFIELD

* external functions:
      INTEGER IROPER

      CALL QENTER('CCXI2')

*---------------------------------------------------------------------*
*     some initializations:
*---------------------------------------------------------------------*
      KEND1  = 1
      LWRK1  = LWORK

      ISYRES = ISYHOP

      IF (LOCDBG) THEN
         WRITE (LUPRI,*) 'CCXI2> IORDER, SKIPXI:',IORDER,SKIPXI
         WRITE (LUPRI,*) 'CCXI2> LRELAX, LTWOEL:',LRELAX,LTWOEL
         WRITE (LUPRI,*) 'CCXI2> IRELAX1, IRELAX2:',IRELAX1,IRELAX2
      END IF

*---------------------------------------------------------------------*
*     transform one-electron hamiltonian and ``barred'' Fock matrix 
*     to MO basis and add extra relaxation contributions:
*---------------------------------------------------------------------*
      KEND2 = KEND1
   
      IF (CC2 .AND. LRELAX) THEN
        KFCKHF0 = KEND2
        KEND2   = KFCKHF0 + N2BST(ISYM0)
        IF (NFIELD.GT.0) THEN
          KFIELD  = KEND2
          KEND2   = KFIELD  + N2BST(ISYM0)
        END IF
      END IF
      LWRK2 = LWORK - KEND2
      IF (LWRK2 .LT. 0) THEN
         CALL QUIT('Insufficient work space in CCXI2. (00)')
      END IF

      IF (CC2 .AND. LRELAX) THEN
        ! read SCF AO-Fock matrix (includes all relaxed fields)
        LUFCK = -1
        CALL GPOPEN(LUFCK,'CC_FCKREF','UNKNOWN',' ','UNFORMATTED',
     *              IDUMMY,.FALSE.)
        REWIND(LUFCK)
        READ(LUFCK)(WORK(KFCKHF0-1+I),I=1,N2BST(ISYM0))
        CALL GPCLOSE(LUFCK,'KEEP')

        IF (NFIELD.GT.0) THEN
          ! collect unrelaxed external fields in WORK(KFIELD)
          CALL DZERO(WORK(KFIELD),N2BST(ISYM0))
          DO IF = 1, NFIELD
            IF ( NHFFIELD(IF) ) THEN
              CALL CC_ONEP(WORK(KFIELD),WORK(KEND2),LWRK2,EFIELD(IF),
     &                     ISYM0,LFIELD(IF))
            END IF
          END DO
        END IF
      END IF

      IF (LOCDBG) THEN
         IF ( CC2 .AND. LRELAX ) THEN 
           WRITE (LUPRI,*) 'CCXI2> FCKHFB matrix in AO:'
           CALL CC_PRFCKMO(FCKHFB,ISYHOP)
           WRITE (LUPRI,*) 'CCXI2> FCKHF0 matrix in AO:'
           CALL CC_PRFCKMO(WORK(KFCKHF0),ISYM0)
         END IF
      END IF

      LRELAX1 = .FALSE.
      LRELAX2 = .FALSE.

      IF (IORDER.GE.2) THEN
          KONEH1   = KEND2  
          KONEH2   = KONEH1  + N2BST(ISYOP1)
          KCMOP1   = KONEH2  + N2BST(ISYOP2)
          KCMOP2   = KCMOP1  + NGLMDT(ISYOP1)
          KCMOH1   = KCMOP2  + NGLMDT(ISYOP2)
          KCMOH2   = KCMOH1  + NGLMDT(ISYOP1)
          KLAMDP1  = KCMOH2  + NGLMDT(ISYOP2)
          KLAMDP2  = KLAMDP1 + NGLMDT(ISYOP1)
          KLAMDH1  = KLAMDP2 + NGLMDT(ISYOP2)
          KLAMDH2  = KLAMDH1 + NGLMDT(ISYOP1)
          KT1AMP0  = KLAMDH2 + NGLMDT(ISYOP2)
          KEND2    = KT1AMP0 + NT1AM(ISYM0)
          LWRK2  = LWORK - KEND2

          IF (LWRK2 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CCXI2. (0a)')
          END IF

          IOPT = 1
          CALL CC_RDRSP('R0 ',0,ISYM0,IOPT,MODEL,
     &                  WORK(KT1AMP0),DUMMY)

          LRELAX1 = (IRELAX1 .GE. 1)
          LRELAX2 = (IRELAX2 .GE. 1)

          IF (LRELAX1) THEN
            CALL CCPRPAO(LABEL2,.TRUE.,WORK(KONEH2),IRREP,
     &                   ISYM,IERR,WORK(KEND2),LWRK2)
            IF (IERR.LT.0) THEN
              CALL DZERO(WORK(KONEH2),N2BST(ISYOP2))
            END IF
          END IF

          IF (LRELAX2) THEN
            CALL CCPRPAO(LABEL1,.TRUE.,WORK(KONEH1),IRREP,
     &                   ISYM,IERR,WORK(KEND2),LWRK2)
            IF (IERR.LT.0) THEN
              CALL DZERO(WORK(KONEH1),N2BST(ISYOP1))
            END IF
          END IF

          IOPT   = 1
          IOPER1 = IROPER(LABEL1,ISYOP1)
          CALL CC_GET_LAMBDA1(IRELAX1,IOPER1,IOPT,
     &                        WORK(KLAMDP1),WORK(KLAMDH1),
     &                        WORK(KCMOP1), WORK(KCMOH1),
     &                        ISYOP1,WORK(KT1AMP0),WORK(KEND2),LWRK2)

          IOPT   = 1
          IOPER2 = IROPER(LABEL2,ISYOP2)
          CALL CC_GET_LAMBDA1(IRELAX2,IOPER2,IOPT,
     &                        WORK(KLAMDP2),WORK(KLAMDH2),
     &                        WORK(KCMOP2), WORK(KCMOH2),
     &                        ISYOP2,WORK(KT1AMP0),WORK(KEND2),LWRK2)

      END IF


      CALL CC_FCKRLX(ONEH0,ONEHB,WORK(KONEH1),WORK(KONEH2),
     &               XLAMDP0,XLAMDH0,ISYM0,
     &               XLAMDQP,XLAMDQH,ISYHOP,LRELAX,
     &               WORK(KLAMDP1),WORK(KLAMDH1),ISYOP1,LRELAX1,
     &               WORK(KLAMDP2),WORK(KLAMDH2),ISYOP2,LRELAX2,
     &               IORDER,WORK(KEND2),LWRK2)


      CALL CC_FCKRLX(FOCK0,FOCKB,WORK(KONEH1),WORK(KONEH2),
     &               XLAMDP0,XLAMDH0,ISYM0,
     &               XLAMDQP,XLAMDQH,ISYHOP,LRELAX,
     &               WORK(KLAMDP1),WORK(KLAMDH1),ISYOP1,LRELAX1,
     &               WORK(KLAMDP2),WORK(KLAMDH2),ISYOP2,LRELAX2,
     &               IORDER,WORK(KEND2),LWRK2)

      IF ( CC2 .AND. LRELAX ) THEN 
         ! for RELAXED CC2 use CMO coefficients to avoid contributions
         ! from t1 amplitudes included in the Lambda matrices in E-term.
         CALL CC_FCKRLX(WORK(KFCKHF0),FCKHFB,WORK(KONEH1),WORK(KONEH2),
     &                  CMO0,CMO0,ISYM0,
     &                  CMOPQ,CMOHQ,ISYHOP,LRELAX,
     &                  WORK(KCMOP1),WORK(KCMOH1),ISYOP1,LRELAX1,
     &                  WORK(KCMOP2),WORK(KCMOH2),ISYOP2,LRELAX2,
     &                  IORDER,WORK(KEND2),LWRK2)

         ! add contribution from unrelaxed external fields, which
         ! has to be transformed using the Lambda matrices
         IF (NFIELD.GT.0) THEN
            CALL CC_FCKMO2(FCKHFB,WORK(KFIELD),ISYM0,
     &                     XLAMDP0,XLAMDH0,ISYM0,
     &                     XLAMDQP,XLAMDQH,ISYHOP,WORK(KEND2),LWRK2)

            ! for 2. order operators there might be something missing?
            IF (LRELAX1.OR.LRELAX2) THEN
              WRITE (LUPRI,*) 
     &            'CCXI2> finite diff. for relaxed CC2 unfinished.'
              CALL QUIT(
     &            'CCXI2> finite diff. for relaxed CC2 unfinished.')
            END IF
         END IF
      END IF

      IF (LOCDBG) THEN
         WRITE (LUPRI,*) 'CCXI2> ONEHB matrix in MO:'
         CALL CC_PRFCKMO(ONEHB,ISYHOP)
         WRITE (LUPRI,*) 'CCXI2> FOCKB matrix in MO:'
         CALL CC_PRFCKMO(FOCKB,ISYHOP)
      END IF

*---------------------------------------------------------------------*
*     initialize single & double excitation result vectors:
*---------------------------------------------------------------------*
      CALL DZERO(RHO1,NT1AM(ISYRES))
      IF (.NOT.CCS) CALL DZERO(RHO2,NT2AM(ISYRES))

*---------------------------------------------------------------------*
*     for CCSD calculate BF-terms and Gamma and G intermediates:
*---------------------------------------------------------------------*
      IF ( (.NOT. (CCS .OR. CC2)) .AND. (LTWOEL.OR.LRELAX)) THEN

C         ---------------------------------------------------
C         allocate work space for BF intermediates:
C         ---------------------------------------------------
          LENBF0 = 2 * NT2ORT(ISYM0)
          LENBF1 = 2 * NT2ORT(ISYHOP)

          KBF    = KEND1
          KEND2  = KBF    + MAX(LENBF0,LENBF1)
          LWRK2  = LWORK - KEND2

          IF (LWRK2 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CCXI2. (1)')
          END IF

C         ---------------------------------------------------------
C         read first-order BF intermediate from file and initialize
C         GAMMA and G intermediates with zero:
C         ---------------------------------------------------------
          CALL GETWA2(LUBFDX1,FNBFDX1,WORK(KBF),IADRX1,LENBF1)

          CALL DZERO(GAMMA,NGAMMA(ISYHOP))
          CALL DZERO(GBIM,NT1AO(ISYHOP))
          CALL DZERO(G0IM,NT1AO(ISYM0))

          IF (LOCDBG) THEN
             XNORM = DDOT(LENBF1,WORK(KBF),1,WORK(KBF),1)
             WRITE (LUPRI,*) 'address/norm^2 of BF1 intermediate:',
     &                        IADRX1,XNORM
          END IF

C         ---------------------------------------------------------
C         transform BF to MO representation using XLAMDP0 matrices:
C         ---------------------------------------------------------
          ICON    = 1
          IOPTG   = 2
          ISYMBF  = ISYRES
          LGAMMA  = .TRUE.
          LO3BF   = .FALSE.
          LBFZETA = .FALSE.

          CALL CC_T2MO3(DUMMY,DUMMY,1,WORK(KBF),RHO2,
     *                  GAMMA,GBIM,DUMMY,
     *                  XLAMDP0,ISYM0,XLAMDP0,ISYM0,WORK(KEND2),LWRK2,
     *                  ISYMBF,ICON,LGAMMA,IOPTG,LO3BF,LBFZETA)

   
          IF ( LRELAX ) THEN

C           -----------------------------------------------------------
C           read zero-order BF intermediate (overwriting the other BF):
C           -----------------------------------------------------------
            LUBF0 = -1
            CALL GPOPEN(LUBF0,'CC_BFIM','OLD',' ','UNFORMATTED',
     &                  IDUMMY,.FALSE.)
            READ(LUBF0) (WORK(KBF-1+I),I=1,LENBF0)
            CALL GPCLOSE(LUBF0,'KEEP')

            IF (LOCDBG) THEN
               XNORM = DDOT(LENBF0,WORK(KBF),1,WORK(KBF),1)
               WRITE (LUPRI,*) 'norm^2 of BF0 intermediate:',XNORM
            END IF

C           ---------------------------------------------------
C           transform to MO using XLAMDP0 and XLAMDQP matrices:
C           ---------------------------------------------------
            ICON    = 4
            IOPTG   = 2
            ISYMBF  = ISYM0
            LGAMMA  = .TRUE.
            LO3BF   = .FALSE.
            LBFZETA = .FALSE.
            CALL CC_T2MO3(DUMMY,DUMMY,1,WORK(KBF),RHO2,
     *                    GAMMA,G0IM,GBIM,
     *                    XLAMDP0,ISYM0,XLAMDQP,ISYHOP,
     *                    WORK(KEND2),LWRK2,ISYMBF,
     *                    ICON,LGAMMA,IOPTG,LO3BF,LBFZETA)

          END IF

          IF (LOCDBG) THEN
             WRITE (LUPRI,*) 'CCXI2> RHO after BF-term section:'
             CALL CC_PRP(RHO1,RHO2,ISYRES,1,1)
          END IF

      END IF
   
*---------------------------------------------------------------------*
*     for CC2 calculate F term (i.e. the (ai|bj)-bar integrals):
*---------------------------------------------------------------------*
      IF ( CC2 .AND. (LTWOEL.OR.LRELAX)) THEN

         KXAIBJ = KEND1
         KEND2  = KXAIBJ + NT2SQ(ISYHOP)
         LWRK2  = LWORK  - KEND2
         IF (LWRK2 .LT. 0) THEN
           CALL QUIT('Insufficient work space in CCXI2. (F-term)')
         END IF

         CALL DZERO(WORK(KXAIBJ),NT2SQ(ISYHOP))

         IOPT = 2
         CALL CC_IAJB2(WORK(KXAIBJ),ISYHOP,IOPT,LRELAX,.FALSE.,.FALSE.,
     &                 LU0AIBJ,FN0AIBJ,IT2DEL0,XLAMDP0,ISYM0,
     &                 LU1AIBJ,FN1AIBJ,IT2DELB,XLAMDQP,ISYHOP,
     &                 WORK(KEND2),LWRK2)

         IOPT = 0
         CALL CC_T2PK(RHO2,WORK(KXAIBJ),ISYHOP,IOPT)

         IF (LOCDBG) THEN
           WRITE (LUPRI,*) 'CCXI2> RHO after F-term section:'
           CALL CC_PRP(RHO1,RHO2,ISYRES,0,1)
           XNORM = DDOT(NT2AM(ISYRES),RHO2,1,RHO2,1)
           WRITE (LUPRI,*) 'Norm^2(RHO2):',XNORM
         END IF

      END IF

*---------------------------------------------------------------------*
*     calculate E1 and E2 intermediates from the R, G & Fock matrices:
*---------------------------------------------------------------------*
      IF (.NOT. CCS) THEN

         IF ( (.NOT. CC2) .AND. (LTWOEL.OR.LRELAX)) THEN
             LRCON  = .TRUE.  ! include contrib. from R 
             LGCON  = .TRUE.  ! include contrib. from G
             FCKCON = .TRUE.  ! include contrib. from Fock 
             IOPT   =   2
             CALL CC_EIM(EMAT1,EMAT2,R0IM,RBIM,G0IM,GBIM,FOCKB,ONEHB,
     *                   XLAMDH0,XLAMDP0,ISYM0,XLAMDQH,XLAMDQP,ISYHOP,
     *                   FCKCON,LRCON,LGCON,LRELAX,IOPT,ISYRES)

         ELSE IF ( CC2 .AND. (LTWOEL.OR.LRELAX)) THEN

c            XNORM = DNRM2(N2BST(ISYRES),FCKHFB,1)
c            WRITE (LUPRI,*) 'norm^2(fckhfb):',xnorm
c            XNORM = DNRM2(NEMAT1(ISYM0),R0IM,1)
c            WRITE (LUPRI,*) 'norm^2(r0im):',xnorm
c            XNORM = DNRM2(NEMAT1(ISYRES),RBIM,1)
c            WRITE (LUPRI,*) 'norm^2(rbim):',xnorm
c            XNORM = DNRM2(NT1AO(ISYM0),G0IM,1)
c            WRITE (LUPRI,*) 'norm^2(g0im):',xnorm
c            XNORM = DNRM2(NT1AO(ISYRES),GBIM,1)
c            WRITE (LUPRI,*) 'norm^2(gbim):',xnorm

             LRCON  = .FALSE. ! skip contrib. from R 
             LGCON  = .FALSE. ! skip contrib. from G
             FCKCON = .TRUE.  ! include contrib. from SCF Fock matrix
             IOPT   =   2
             CALL CC_EIM(EMAT1,EMAT2,R0IM,RBIM,G0IM,GBIM,FCKHFB,FCKHFB,
     *                   XLAMDH0,XLAMDP0,ISYM0,XLAMDQH,XLAMDQP,ISYHOP,
     *                   FCKCON,LRCON,LGCON,LRELAX,IOPT,ISYRES)

         ELSE ! unrelaxed CC2 or CCSD 
             LRCON  = .FALSE. ! skip contrib. from R 
             LGCON  = .FALSE. ! skip contrib. from G
             FCKCON = .TRUE.  ! include contrib. from ONEHB
             IOPT   =   2
             CALL CC_EIM(EMAT1,EMAT2,R0IM,RBIM,G0IM,GBIM,ONEHB,ONEHB,
     *                   XLAMDH0,XLAMDP0,ISYM0,XLAMDQH,XLAMDQP,ISYHOP,
     *                   FCKCON,LRCON,LGCON,LRELAX,IOPT,ISYRES)
         END IF

         IF (LOCDBG) THEN
            CALL AROUND('E intermediates in CCXI2:')
            XNORM = DDOT(NMATAB(ISYHOP),EMAT1,1,EMAT1,1)
            WRITE (LUPRI,*) 'Norm of EMAT1:',XNORM
            XNORM = DDOT(NMATIJ(ISYHOP),EMAT2,1,EMAT2,1)
            WRITE (LUPRI,*) 'Norm of EMAT2:',XNORM
            WRITE (LUPRI,*) 'LRCON,LGCON,FCKCON:',LRCON,LGCON,FCKCON
            CALL CC_PREI(EMAT1,EMAT2,ISYRES,1)
         END IF

      END IF

*---------------------------------------------------------------------*
*     read zeroth-order T2 amplitudes and square them up:
*---------------------------------------------------------------------*
      IF (.NOT. CCS) THEN

          KT2AMP0 = KEND1
          KEND2   = KT2AMP0 + NT2SQ(ISYM0)
          LWRK2   = LWORK   - KEND2

          IF (LWRK2 .LT. NT2AM(ISYM0)) THEN
            CALL QUIT('Insufficient work space in CCXI2. (A-term)')
          END IF

          IOPT = 2
          CALL CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,DUMMY,
     &                  WORK(KEND2))
          CALL CC_T2SQ(WORK(KEND2),WORK(KT2AMP0),ISYM0)

      END IF

*---------------------------------------------------------------------*
*     If CC2 or CCSD and no two-electron integral contributions 
*     (i.e. the E1 & E2 terms not included in A & B terms) then
*     calculate the E terms here:
*---------------------------------------------------------------------*
      IF ( (.NOT.SKIPXI) .AND. (.NOT.CCS) ) THEN
        IF ( CC2 .OR. (.NOT.(LTWOEL.OR.LRELAX)) ) THEN

           CALL CCRHS_E(RHO2,WORK(KT2AMP0),EMAT1,EMAT2,
     &                  WORK(KEND2),LWRK2,ISYM0,ISYHOP)

           IF (LOCDBG) THEN
             WRITE (LUPRI,*) 'CCXI2> RHO2 after E-term section:'
             CALL CC_PRP(RHO1,RHO2,ISYRES,0,1)
             XNORM = DDOT(NT2AM(ISYRES),RHO2,1,RHO2,1)
             WRITE (LUPRI,*) 'Norm^2(RHO2) after E-term section:',XNORM
           END IF

        END IF
      END IF

*---------------------------------------------------------------------*
*     calculate A-term contribution:
*       requires:  -- Gamma intermediate in core
*                  -- T2 amplitudes squared in core
*---------------------------------------------------------------------*
      IF ( (.NOT. (CCS .OR. CC2)) .AND. (LTWOEL.OR.LRELAX)) THEN

          ! add E2 intermediate to the diagonals of Gamma:
          CALL CC_GAMMA2(GAMMA,EMAT2,ISYRES)

          ! calculate the A-term:
          IF (.NOT. SKIPXI) THEN
            IOPT = 1
            CALL CCRHS_A(RHO2,WORK(KT2AMP0),GAMMA,
     *                   WORK(KEND2),LWRK2,ISYRES,ISYM0,IOPT)

            IF (LOCDBG) THEN
               WRITE (LUPRI,*) 'CCXI2> Rho after A-term section:'
               CALL CC_PRP(RHO1,RHO2,ISYRES,1,1)
            END IF
          END IF

      END IF

*---------------------------------------------------------------------*
*     calculate I-term contributions:
*       requires:  -- MO Fock matrix in core
*                  -- T2 amplitudes squared in core
*---------------------------------------------------------------------*
      IF ( (.NOT. SKIPXI) .AND. (.NOT. CCS) ) THEN

         CALL CCRHS_T2TR(WORK(KT2AMP0),WORK(KEND2),LWRK2,1)

         CALL CCRHS_I(RHO1,WORK(KT2AMP0),FOCKB, 
     &                WORK(KEND2),LWRK2,ISYM0,ISYHOP)

         CALL CCRHS_T2BT(WORK(KT2AMP0),WORK(KEND2),LWRK2,1)

         IF (LOCDBG) THEN
             WRITE (LUPRI,*) 'RHO1 after CCRHS_I:'
             CALL CC_PRP(RHO1,RHO2,ISYRES,1,0)
             XNORM = DDOT(NT1AM(ISYRES),RHO1,1,RHO1,1)
             WRITE (LUPRI,*) 'Norm of RHO1 after CCRHS_I:',XNORM
         END IF
      END IF

*---------------------------------------------------------------------*
*     calculate J, G and H-term contributions:
*       requires:  -- Fock matrix 
*                  -- G intermediates
*                  -- R intermediates
*---------------------------------------------------------------------*
      IF (.NOT. SKIPXI) THEN

        IF ( (.NOT. CCS) .AND. (LTWOEL.OR.LRELAX)) THEN
           LRCON  = .TRUE.   ! include contrib. from R 
           LGCON  = .TRUE.   ! include contrib. from G
           FCKCON = .TRUE.
           LTRSPF = .FALSE.
           IOPT   =   2
           CALL CC_GHJ(RHO1,R0IM,RBIM,G0IM,GBIM,ONEHB,
     *                 XLAMDH0,XLAMDP0,ISYM0,
     *                 XLAMDQH,XLAMDQP,ISYHOP,
     *                 FCKCON,LRCON,LGCON,LRELAX,LTRSPF,IOPT,ISYRES)

        ELSE IF ( CCS .AND. (LTWOEL.OR.LRELAX)) THEN
           LRCON  = .FALSE.  ! skip contrib. from R 
           LGCON  = .FALSE.  ! skip contrib. from G
           FCKCON = .TRUE.
           LTRSPF = .FALSE.
           IOPT   =   2
           CALL CC_GHJ(RHO1,R0IM,RBIM,G0IM,GBIM,FOCKB,
     *                 XLAMDH0,XLAMDP0,ISYM0,
     *                 XLAMDQH,XLAMDQP,ISYHOP,
     *                 FCKCON,LRCON,LGCON,LRELAX,LTRSPF,IOPT,ISYRES)

        ELSE ! unrelaxed CCS, CC2 or CCSD
           LRCON  = .FALSE.  ! skip contrib. from R 
           LGCON  = .FALSE.  ! skip contrib. from G
           FCKCON = .TRUE.
           LTRSPF = .FALSE.
           IOPT   =   2
           CALL CC_GHJ(RHO1,R0IM,RBIM,G0IM,GBIM,ONEHB,
     *                 XLAMDH0,XLAMDP0,ISYM0,
     *                 XLAMDQH,XLAMDQP,ISYHOP,
     *                 FCKCON,LRCON,LGCON,LRELAX,LTRSPF,IOPT,ISYRES)
        END IF
     
        IF (LOCDBG) THEN
          WRITE (LUPRI,*) 'Finished GHJ-term section... Rho1:'
          CALL CC_PRP(RHO1,RHO2,ISYRES,1,0)
          XNORM = DDOT(NT1AM(ISYRES),RHO1,1,RHO1,1)
          WRITE (LUPRI,*) 'Norm of RHO1 after CC_GHJ:',XNORM
          XNORM = DDOT(NT2AM(ISYRES),RHO2,1,RHO2,1)
          WRITE (LUPRI,*) 'Norm of RHO2 after CC_GHJ:',XNORM
        END IF

      END IF

*----------------------------------------------------------------------*
*     read (ia|j del) and (ia|j del)-bar integrals and transform then
*     to (ia|jb)-bar and read packed T2 amplitudes:
*----------------------------------------------------------------------*
      IF ( (.NOT.CCS) .AND. (LTWOEL.OR.LRELAX)) THEN

         ISYIAJB = ISYHOP

         KXIAJB = KEND1
         KEND1  = KXIAJB + NT2SQ(ISYIAJB)
         IF (.NOT.CC2) THEN
           KT2AM  = KEND1
           KEND2  = KT2AM  + NT2AM(ISYM0)
         END IF
         LWRK2  = LWORK  - KEND2

         IF (LWRK2 .LT. 0) THEN
           CALL QUIT('Insufficient work space in CCXI2. (IAJB)')
         END IF

         CALL DZERO(WORK(KXIAJB),NT2SQ(ISYIAJB))

         IOPT = 2
         LOCC = .FALSE.
         CALL CC_IAJB2(WORK(KXIAJB),ISYIAJB,IOPT,LRELAX,LOCC,.FALSE.,
     &                 LU0IAJB,FN0IAJB,IT2DEL0,XLAMDH0,ISYM0,
     &                 LU1IAJB,FN1IAJB,IT2DELB,XLAMDQH,ISYHOP,
     &                 WORK(KEND2),LWRK2)

* write (ia|jb)-bar integrals back to file for later use:
* (append after (ia|jdelta) integrals)
         CALL PUTWA2(LU0IAJB,FN0IAJB,WORK(KXIAJB),
     &               IADRINT,NT2SQ(ISYIAJB))

         IF (LOCDBG) THEN
            WRITE (LUPRI,*) '(ia|jb) integrals:'
            CALL CC_PRSQ(DUMMY,WORK(KXIAJB),ISYRES,0,1)
         END IF

         IF (.NOT.CC2) THEN
C          --------------------------
C          read amplitudes from file:
C          --------------------------
           IOPT = 2
           CALL CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,DUMMY,
     *                   WORK(KT2AM))
         END IF

      END IF
*----------------------------------------------------------------------*
*       calculate C-term contribution:
*         requires:  -- (ia|jb)-bar integrals squared in core
*                    -- T2 amplitudes packed in core
*                    -- C intermediate in core (squared matrix)
*                    -- result vector in core (packed)
*----------------------------------------------------------------------*
      IF ( (.NOT. (CCS .OR. CC2)) .AND. (LTWOEL.OR.LRELAX)) THEN
       
          KCBAR   = KEND2
          KEND3   = KCBAR + NT2SQ(ISYIAJB)  
          LWRK3   = LWORK - KEND3
    
          IF (LWRK3 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CCX2. (C-term)')
          END IF

* contract (ia|jb)-bar integrals with T2 amplitudes:
          IOPT = 1
          CALL CCB_CDBAR('C',WORK(KXIAJB),ISYIAJB, WORK(KT2AM),ISYM0,
     &                       WORK(KCBAR),ISYIAJB, WORK(KEND3),LWRK3,
     &                       CDUM, IDUM, IDUM, IOPT)


* for XKSI vector scale with factor 1/2:
          CALL DSCAL(NT2SQ(ISYIAJB),HALF,WORK(KCBAR),1)

* save this intermediate result for calculation of the C intermediate
* needed for the ETA vector:
          CALL DCOPY(NT2SQ(ISYIAJB),WORK(KCBAR),1,WORK(KXIAJB),1)

C         -----------------------------------------------------------
C         read (ki|del c) integrals into memory, transform to (ki|ac)
C         using XLAMDP0/QP matrices and add to the CBAR intermediate:
C         -----------------------------------------------------------
          IOPT = 2
          LOCC = .FALSE.
          CALL CC_IAJB2(WORK(KCBAR),ISYIAJB,IOPT,LRELAX,LOCC,.FALSE.,
     &                  LU0IJBA,FN0IJBA,IT2DEL0,XLAMDP0,ISYM0,
     &                  LU1IJBA,FN1IJBA,IT2DELB,XLAMDQP,ISYHOP,
     &                  WORK(KEND3),LWRK3)


C         ------------------------------------------------------------
C         add E1 intermediate to the diagonal of the C intermediate:
C         ------------------------------------------------------------
          FAC = -0.5D0
          CALL CC_CDBAR2(WORK(KCBAR),EMAT1,FAC,.FALSE.,ISYRES)


C         ------------------------------------------------------------
C         save the C intermediate needed for the ETA vector on file:
C         ------------------------------------------------------------
          CALL DAXPY(NT2SQ(ISYIAJB),ONE,WORK(KCBAR),1,WORK(KXIAJB),1)
          CALL PUTWA2(LUCIM,FNCIM,WORK(KXIAJB),IADRCI,NT2SQ(ISYIAJB))

C         ------------------------------------------
C         calculate C term and add to result vector:
C         ------------------------------------------
          IF (.NOT. SKIPXI) THEN
             ioptr12 = 0
             CALL CC_CD('C',+1,ioptr12,RHO2,ISYRES,WORK(KT2AM),
     &                  ISYM0,WORK(KCBAR),ISYIAJB,WORK(KEND3),LWRK3)

             IF (LOCDBG) THEN
                WRITE (LUPRI,*) 'ISYRES,ISYIAJB:',ISYRES,ISYIAJB
                WRITE (LUPRI,*) 'Finished C-term section... RHO:'
                CALL CC_PRP(RHO1,RHO2,ISYRES,1,1)
                XNORM = DDOT(NT2AM(ISYRES),RHO2,1,RHO2,1)
                WRITE (LUPRI,*) 'Norm of RHO2 after C-term:',XNORM
             END IF
          END IF

      END IF
*---------------------------------------------------------------------*
*       calculate D term contribution:
*         requires:  -- (ia|jb)-bar integrals squared in core
*                    -- T2 amplitudes packed in core
*                    -- D intermediate in core (squared matrix)
*                    -- result vector in core (packed)
*---------------------------------------------------------------------*
      IF ( (.NOT. (CCS .OR. CC2)) .AND. (LTWOEL.OR.LRELAX)) THEN

         KDBAR   = KEND2
         KEND3   = KDBAR + NT2SQ(ISYIAJB)  
         LWRK3   = LWORK - KEND3
    
         IF (LWRK3 .LT. 0) THEN
           CALL QUIT('Insufficient work space in CCXI2. (D-term)')
         END IF

* read (ia|jb)-bar integrals from file:
         CALL GETWA2(LU0IAJB,FN0IAJB,WORK(KXIAJB),
     &               IADRINT,NT2SQ(ISYIAJB))

* contract (ia|jb)-bar integrals with T2 amplitudes:
* (this overwrites T2 with 2T2(ia;jb) - T2(ja;ib)...)
         IOPT = 1
         CALL CCB_CDBAR('D',WORK(KXIAJB),ISYIAJB, WORK(KT2AM),ISYM0,
     &                      WORK(KDBAR),ISYIAJB, WORK(KEND3),LWRK3,
     &                      CDUM, IDUM, IDUM, IOPT)

* for XKSI vector scale with factor 1/2:
         CALL DSCAL(NT2SQ(ISYIAJB),HALF,WORK(KDBAR),1)

* save this intermediate result for calculation of the D intermediate
* needed for the ETA vector:
         CALL DCOPY(NT2SQ(ISYIAJB),WORK(KDBAR),1,WORK(KXIAJB),1)

C        -------------------------------------------------------------
C        read L(kc|del i) integrals into memory, transform to L(kc|ai)
C        using XLAMDP0/QP matrices and add to the CBAR intermediate:
C        -------------------------------------------------------------
         IOPT = 2
         LOCC = .FALSE.
         CALL CC_IAJB2(WORK(KDBAR),ISYIAJB,IOPT,LRELAX,LOCC,.FALSE.,
     &                 LU0IABJ,FN0IABJ,IT2DEL0,XLAMDP0,ISYM0,
     &                 LU1IABJ,FN1IABJ,IT2DELB,XLAMDQP,ISYHOP,
     &                 WORK(KEND3),LWRK3)

C        ------------------------------------------------------------
C        add E1 intermediate to the diagonal of the D intermediate:
C        ------------------------------------------------------------
         FAC = 0.5D0
         CALL CC_CDBAR2(WORK(KDBAR),EMAT1,FAC,.FALSE.,ISYRES)

C        ------------------------------------------------------------
C        save the D intermediate on file:
C        (this overwrites the (ia|del j)-bar integrals...)
C        ------------------------------------------------------------
         CALL DAXPY(NT2SQ(ISYIAJB),ONE,WORK(KDBAR),1,WORK(KXIAJB),1)
         CALL PUTWA2(LUDIM,FNDIM,WORK(KXIAJB),IADRDI,NT2SQ(ISYIAJB))

C        ------------------------------------------------------------
C        calculate D term and add to result vector:
C        (note, that this assumes 2T(ia;jb)-T(ja;ib) in T amplitudes)
C        ------------------------------------------------------------
         IF (.NOT. SKIPXI) THEN
            ioptr12 = 0
            CALL CC_CD('D',+1,ioptr12,RHO2,ISYRES,WORK(KT2AM),
     &              ISYM0,WORK(KDBAR),ISYIAJB,WORK(KEND3),LWRK3)

            IF (LOCDBG) THEN
               WRITE (LUPRI,*) 'Finished D-term section...RHO1:'
               CALL CC_PRP(RHO1,RHO2,ISYRES,1,1)
               XNORM = DDOT(NT2AM(ISYRES),RHO2,1,RHO2,1)
               WRITE (LUPRI,*) 'Norm of RHO2 after D-term:',XNORM
            END IF
         END IF

      END IF

*----------------------------------------------------------------------*
*     Xi vector completed; print debug output and return:
*----------------------------------------------------------------------*
      IF (LOCDBG .AND. .NOT. SKIPXI) THEN
          I = 1
          IF (CCS) I = 0
          WRITE (LUPRI,*) 'Final result of CCXI2:'
          CALL CC_PRP(RHO1,RHO2,ISYRES,1,I)
          CALL FLSHFO(LUPRI)
      END IF

      CALL QEXIT('CCXI2')

      RETURN
      END  
*======================================================================*
*                        END OF SUBROUTINE CCXI2                       *
*======================================================================*
*======================================================================*
c /* deck cceta2 */
*======================================================================*
      SUBROUTINE CCETA2(RHO1,RHO2,AVERAGE,
     &                  G0IM,GBIM,R0IM,RBIM,F0IM,FBIM,RZ0I,RZBI,
     &                  FOCK,ONEH,FOCKHF,ZFCK0,ZFCKB,
     &                  XGAMMA,XINT,YINT,ZETA1,
     &                  LUBFX1,FNBFX1,IADRX1,
     &                  LUBFZI,FNBFZI,IADRE0,IADRE1,
     &                  LUMO,  FNMO,  IADRMO,
     &                  LUDIM, FNDIM, IDRDIM,
     &                  LUCIM, FNCIM, IDRCIM,
     &                  LU0IAJB,FN0IAJB,IT2DEL0,
     &                  LU1IAJB,FN1IAJB,IT2DELB,
     &                  LUHINT0,FNHINT0,IADRH0,
     &                  LUHINT1,FNHINT1,IADRH1,
     &                  LUPQMO, FILPQMO,IADRPQMO,
     &                  XLAMDP0,XLAMDH0,XLAMDQP,XLAMDQH,
     &                  LISTL,IDLSTL,LABELH,ISYHOP,
     &                  LTWOEL,LRELAX,WORK, LWORK)
*----------------------------------------------------------------------*
*
*   Purpose: calculate contributions for ETA vector from intermediates
*
*   Written by Christof Haettig, November 1998
*
*----------------------------------------------------------------------*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdinp.h"
#include "maxorb.h"
#include "ccisao.h"
#include "ccfield.h"
#include "dummy.h"

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER ISYM0, LUBF0
      PARAMETER (ISYM0 = 1)

* variables:
      CHARACTER*(*) LISTL, LABELH
      CHARACTER*(*) FNBFX1, FNBFZI, FNMO, FNDIM ,FNCIM
      CHARACTER*(*) FN1IAJB,FN0IAJB, FILPQMO, FNHINT0, FNHINT1
      LOGICAL LRELAX, LTWOEL
      INTEGER LWORK, IDLSTL, ISYHOP, LU1IAJB, LU0IAJB, LUPQMO
      INTEGER LUBFX1,IADRX1,LUBFZI,IADRE0,IADRE1,LUMO,IADRMO
      INTEGER LUCIM,IDRCIM,LUDIM,IDRDIM, LUHINT0, LUHINT1
      INTEGER IT2DEL0(*), IT2DELB(*), IADRPQMO(*), IADRH0(*), IADRH1(*)

#if defined (SYS_CRAY)
      REAL XLAMDP0(*), XLAMDH0(*), XLAMDQP(*), XLAMDQH(*)
      REAL RHO1(*), RHO2(*), G0IM(*), GBIM(*)
      REAL R0IM(*), RBIM(*), F0IM(*), FBIM(*)
      REAL RZ0I(*), RZBI(*), ZFCK0(*), ZFCKB(*), FOCKHF(*)
      REAL FOCK(*), ONEH(*), XGAMMA(*), XINT(*), YINT(*)
      REAL WORK(LWORK), ZETA1(*)
      REAL HALF, ONE, ZERO, DDOT, XNORM
      REAL AVE1, AVE2, AVERAGE
#else
      DOUBLE PRECISION XLAMDP0(*), XLAMDH0(*), XLAMDQP(*), XLAMDQH(*)
      DOUBLE PRECISION RHO1(*), RHO2(*), G0IM(*), GBIM(*)
      DOUBLE PRECISION R0IM(*), RBIM(*), F0IM(*), FBIM(*)
      DOUBLE PRECISION RZ0I(*), RZBI(*), ZFCK0(*), ZFCKB(*), FOCKHF(*)
      DOUBLE PRECISION FOCK(*), ONEH(*), XGAMMA(*), XINT(*), YINT(*)
      DOUBLE PRECISION WORK(LWORK), ZETA1(*)
      DOUBLE PRECISION HALF, ONE, ZERO, TWO, DDOT, XNORM
      DOUBLE PRECISION AVE1, AVE2, AVERAGE
#endif
      PARAMETER (HALF = 0.5d0, ONE = 1.0d0, TWO = 2.0d0, ZERO = 0.0d0)

      CHARACTER*(10) MODEL
      LOGICAL LGAMMA, LO3BF, LBFZETA, LRCON, LGCON, FCKCON, LOCC, LTRSPF
      INTEGER ISYCTR, ISYRES, ISYETA, KZETA2, KZETPK, KBF, KEND3A
      INTEGER ISYMI, ISYMA, KOFF1, KOFF2, LENBF0, LENBF1, KEND3, LWRK3
      INTEGER IOPT, ICON, KA2IM, KA2RS, KXPCK, KXIAJB, KT2AMP, IOPTG
      INTEGER KB2CON, KEMAT1, KEMAT2, KCTR2, KEND2, LWRK2, MT2BCD
      INTEGER ISYIFJL, KXIFJL, KYJLIF, KPINT0, KQINT0, KWINT0
      INTEGER ISYMF, ISYALJ, ISYJLI, IVIRF, KYPS0, KYPS1, LEN
      INTEGER NTAIKJ(9), ITAIKJ(8,8), ICOUNT, ISYM, ISYMAIK, ISYMJ
      INTEGER ISYMBF, KCDBAR, IADR, LENBFZ0, LENBFZ1, KBFZ, IORD
      INTEGER ILSTSYM, IOPTTCME, IOPTR12

      CALL QENTER('CCETA2')

*----------------------------------------------------------------------*
* begin:
*----------------------------------------------------------------------*
      ISYCTR = ILSTSYM(LISTL,IDLSTL)
      ISYETA = MULD2H(ISYHOP,ISYCTR)
      ISYRES = ISYETA
     
* intialize single excitation part of the ETA vector:
      CALL DZERO(RHO1,NT1AM(ISYETA))

* set NTAIKJ and ITAIKJ arrays:
      DO ISYM = 1, NSYM
        ICOUNT = 0
        DO ISYMAIK = 1, NSYM
           ISYMJ  = MULD2H(ISYMAIK,ISYM)
           ITAIKJ(ISYMAIK,ISYMJ) = ICOUNT
           ICOUNT = ICOUNT + NT2BCD(ISYMAIK)*NRHF(ISYMJ)
        END DO
        NTAIKJ(ISYM) = ICOUNT
      END DO

* initialize AVERAGE:
      AVERAGE = ZERO

*----------------------------------------------------------------------*
* for relaxed CC2 transform Zeta Fock matrix to MO basis:
*----------------------------------------------------------------------*
      IF ( (CCS .OR. CC2) .AND. (LRELAX .OR. LTWOEL) ) THEN

        CALL CC_FCKRLX1(ZFCKB,ZFCK0,ISYETA,ISYCTR,
     &                  XLAMDP0,XLAMDH0,ISYM0,ISYM0,
     &                  XLAMDQP,XLAMDQH,ISYHOP,ISYHOP,
     &                  LRELAX,WORK,LWORK)

        IF (LOCDBG) THEN
          WRITE (LUPRI,*) 'ZFCKB matrix in MO:'
          CALL CC_PRONELAO(ZFCKB,ISYETA)
        END IF
      END IF

*----------------------------------------------------------------------*
*     for 'L0' calculate projection of singles onto <HF| state:
*     !!! check that FOCK contains the correct Fock matrix for CC2!!!
*     calculate first contribution to <HF|J^(1)|CC> average
*----------------------------------------------------------------------*
      IF ( LISTL(1:2).EQ.'L0' ) THEN

         DO ISYMI = 1, NSYM
            ISYMA = MULD2H(ISYMI,ISYHOP)
            DO I = 1, NRHF(ISYMI)
               DO A = 1, NVIR(ISYMA)
                  KOFF1 = IT1AM(ISYMA,ISYMI)  + NVIR(ISYMA)*(I-1)+A
                  KOFF2 = IFCVIR(ISYMI,ISYMA) + NORB(ISYMI)*(A-1)+I
                  RHO1(KOFF1) = TWO * FOCK(KOFF2)
               END DO
            END DO
         END DO

         IF (ISYHOP.EQ.1) THEN
            AVE1 = ZERO
            AVE2 = ZERO
            DO ISYMI = 1, NSYM
               DO I = 1, NRHF(ISYMI)
                  KOFF2   = IFCRHF(ISYMI,ISYMI) + NORB(ISYMI)*(I-1)+I
                  AVE1    = AVE1 + ONEH(KOFF2) + FOCK(KOFF2)
               END DO
            END DO
            AVERAGE = AVE1 
         END IF

         IF (LOCDBG) THEN
            WRITE (LUPRI,*) 'Rho after singles <HF| term:'
            CALL CC_PRP(RHO1,RHO2,ISYETA,1,0)
            WRITE (LUPRI,*) '<HF|Jhat^(1)|HF> term:',AVE1
            WRITE (LUPRI,*) 'average after <HF|Jhat^(1)|HF> term:',
     &           AVERAGE
         END IF

      ELSE
          IF (LOCDBG) THEN
             WRITE (LUPRI,*) 'no singles <HF| term included.'
          END IF
      END IF

*---------------------------------------------------------------------*
*     calculate single excitation E2 terms form BF intermediates:
*---------------------------------------------------------------------*
      IF ( (.NOT. (CCS .OR. CC2)) .AND. (LTWOEL.OR.LRELAX) ) THEN

        LENBF0  = 2 * NT2ORT(ISYM0)
        LENBF1  = 2 * NT2ORT(ISYHOP)

        KZETA2  = 1
        KBF     = KZETA2 + NT2SQ(ISYCTR)
        KZETPK  = KBF    + MAX(LENBF0,LENBF1)
        KEND3   = KZETPK + NT2AM(ISYCTR)
        LWRK3   = LWORK  - KEND3

        IF (LWRK3 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CCETA2. (CTR2)')
        END IF

C       --------------------------------------------------
C       read lagrange multipliers from file and square up:
C       --------------------------------------------------
        IOPT = 2
        CALL CC_RDRSP(LISTL,IDLSTL,ISYCTR,IOPT,MODEL,DUMMY,
     &                WORK(KZETPK))

        CALL CC_T2SQ(WORK(KZETPK),WORK(KZETA2),ISYCTR)

C       ---------------------------------
C       read first-order BF intermediate:
C       ---------------------------------
        CALL GETWA2(LUBFX1,FNBFX1,WORK(KBF),IADRX1,LENBF1)

C       ------------------------------------------------------
C       transform to MO representation using XLAMDP0 matrices:
C       ------------------------------------------------------
        ICON    = 3
        IOPTG   = 0
        LGAMMA  = .FALSE.
        LO3BF   = .FALSE.
        LBFZETA = .FALSE.
        CALL CC_T2MO3(RHO1,WORK(KZETA2),ISYCTR,
     &                WORK(KBF),DUMMY,DUMMY,DUMMY,DUMMY,
     &                XLAMDP0,ISYM0,XLAMDP0,ISYM0,
     &                WORK(KEND3),LWRK3,ISYHOP,
     &                ICON,LGAMMA,IOPTG,LO3BF,LBFZETA)

        IF (LOCDBG) THEN
          CALL AROUND("RHO1 after E2 term contributions (part1):")
          CALL CC_PRP(RHO1,WORK,ISYETA,1,0)
        END IF

        IF ( LRELAX ) THEN

C            ----------------------------------
C            read zeroth-order BF intermediate:
C            ----------------------------------
             LUBF0 = -1
             CALL GPOPEN(LUBF0,'CC_BFIM','OLD',' ','UNFORMATTED',
     &                   IDUMMY,.FALSE.)
             READ(LUBF0) (WORK(KBF-1+I),I=1,2*NT2ORT(ISYM0))
             CALL GPCLOSE(LUBF0,'KEEP')

C            ------------------------------------
C            transform to MO representation using 
C            XLAMDP0 and XLAMDQP matrices:
C            ------------------------------------
             ICON    = 6
             IOPTG   = 0
             LGAMMA  = .FALSE.
             LO3BF   = .FALSE.
             LBFZETA = .FALSE.
             CALL CC_T2MO3(RHO1,WORK(KZETA2),ISYCTR,
     &                     WORK(KBF),DUMMY,DUMMY,DUMMY,DUMMY,
     &                     XLAMDP0,ISYM0,XLAMDQP,ISYHOP,
     &                     WORK(KEND3),LWRK3,ISYM0,
     &                     ICON,LGAMMA,IOPTG,LO3BF,LBFZETA)

        END IF

        IF (LOCDBG) THEN
          CALL AROUND("RHO1 after E2 term contributions:")
          CALL CC_PRP(RHO1,WORK,ISYETA,1,0)
        END IF

      END IF

*---------------------------------------------------------------------*
*     calculate A2 contributions:
*       requires:  -- T2 amplitudes or (ia|jb) packed in core
*---------------------------------------------------------------------*
      IF ( (.NOT.CCS) .AND. (LTWOEL.OR.LRELAX) ) THEN

          KA2IM  = 1
          KA2RS  = KA2IM  + NT1AM(ISYCTR)
          KXPCK  = KA2RS  + NT1AM(ISYRES)
          KXIAJB = KXPCK  + NT2AM(ISYHOP)
          KEND3  = KXIAJB + NT2SQ(ISYHOP)

          KT2AMP = KXPCK
          KEND3  = MAX(KEND3,KT2AMP + NT2AM(ISYM0))

          LWRK3  = LWORK  - KEND3

          IF (LWRK3 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CCETA2. (CTR2)')
          END IF

C         ---------------------------------------------------
C         read double excitation amplitudes; construct 2T - T
C         ---------------------------------------------------
          IOPT = 2
          CALL CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,
     &                  DUMMY,WORK(KT2AMP))
          IOPTTCME = 1
          CALL CCSD_TCMEPK(WORK(KT2AMP),ONE,ISYM0,IOPTTCME)

C         -----------------------------------------
C         contract: A2IM = sum_bj (2T-T)_aibj L1_bj
C         -----------------------------------------
          CALL CCG_LXD(WORK(KA2IM),ISYCTR,ZETA1,ISYCTR, 
     &                 WORK(KT2AMP),ISYM0,0)

C         ----------------------------------------------------------
C         read (ia|jb) derivative integrals; pack and calc. L(ia,jb)
C         !!! should be done with squared integrals... CCG_LXD is
C         !!! now able to do this!
C         ----------------------------------------------------------
          IOPT   = 0
          CALL GETWA2(LUMO,FNMO,WORK(KXIAJB),IADRMO,NT2SQ(ISYHOP))
          CALL CC_T2PK(WORK(KXPCK),WORK(KXIAJB),ISYHOP,IOPT)
          IOPTTCME = 1
          CALL CCSD_TCMEPK(WORK(KXPCK),ONE,ISYHOP,IOPTTCME)

C         -------------------------------------------
C         contract: eta_ai = sum_bj L(ia,jb) A2IM(jb)
C         -------------------------------------------
          CALL CCG_LXD(WORK(KA2RS),ISYETA,WORK(KA2IM),ISYCTR, 
     &                 WORK(KXPCK),ISYHOP,0)
          CALL DAXPY(NT1AM(ISYRES),ONE,WORK(KA2RS),1,RHO1,1)

          IF (LOCDBG) THEN
              WRITE (LUPRI,*) 'RHO1 after A2 contribution:'
              CALL CC_PRP(RHO1,RHO2,ISYRES,1,0)
              XNORM = DDOT(NT1AM(ISYRES),RHO1,1,RHO1,1)
              WRITE (LUPRI,*) 
     &             'Norm of RHO1 after A2 contribution:',XNORM
          END IF

      END IF

*----------------------------------------------------------------------*
*     for 'L0' and relaxed or two-electron perturbation calculate 
*     projection of doubles onto <HF| state; in any case initialize
*     the double excitation part of the ETA vector:
*----------------------------------------------------------------------*
      IF ( LISTL(1:2).EQ.'L0' .AND. (.NOT.CCS) .AND. 
     &     (LTWOEL.OR.LRELAX)                        ) THEN

          KXIAJB = 1
          KEND3  = KXIAJB + NT2SQ(ISYHOP)

          IF (ISYHOP.EQ.1) THEN
            KT2AMP = 1
            KEND3A = KT2AMP + NT2AM(ISYM0)
            KEND3  = MAX(KEND3,KEND3A)
          END IF

          LWRK3  = LWORK  - KEND3
          IF (LWRK3.LT.0) THEN
            CALL QUIT('Insufficient work space in CCETA2. (B2)')
          END IF

C         ----------------------------------------------------------
C         read (ia|jb) derivative integrals; pack and calc. L(ia,jb)
C         ----------------------------------------------------------
          IOPT   = 0
          CALL GETWA2(LUMO,FNMO,WORK(KXIAJB),IADRMO,NT2SQ(ISYHOP))
          CALL CC_T2PK(RHO2,WORK(KXIAJB),ISYHOP,IOPT)
          IOPTTCME = 1
          CALL CCSD_TCMEPK(RHO2,ONE,ISYHOP,IOPTTCME)
          CALL DSCAL(NT2AM(ISYETA),TWO,RHO2,1)

          IF (ISYHOP.EQ.1) THEN
             IOPT   = 2
             CALL CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,
     &                     DUMMY,WORK(KT2AMP))
             CALL CCLR_DIASCL(WORK(KT2AMP),HALF,ISYM0)
             AVE2 = DDOT(NT2AMX,WORK(KT2AMP),1,RHO2,1)
             AVERAGE = AVERAGE + AVE2
          END IF

          IF (LOCDBG) THEN
             WRITE (LUPRI,*) 'Rho after doubles <HF| term:'
             CALL CC_PRP(RHO1,RHO2,ISYETA,1,1)
             WRITE (LUPRI,*) '<HF|J^(1) T_2|HF> term:',AVE2
             WRITE (LUPRI,*) 
     &             'average after <HF|J^(1) T_2|HF> term:',AVERAGE
          END IF

      ELSE IF ( .NOT. CCS ) THEN
          CALL DZERO(RHO2,NT2AM(ISYETA))
      END IF


*----------------------------------------------------------------------*
*     calculate B2 contributions:
*----------------------------------------------------------------------*
      KB2CON = 1
      KEMAT1 = KB2CON + NT1AM(ISYETA)
      KEMAT2 = KEMAT1 + NMATAB(ISYHOP)
      KEND3  = KEMAT2 + NMATIJ(ISYHOP)
      LWRK3  = LWORK  - KEND3

      IF (LWRK3 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CC_XIETA. (B2)')
      END IF

      IF ( (.NOT.(CCS.OR.CC2)) .AND. (LTWOEL.OR.LRELAX) ) THEN
           LRCON  = .TRUE.  ! include contrib. from R 
           LGCON  = .TRUE.  ! include contrib. from G
           FCKCON = .TRUE.  ! include contrib. from Fock
           IOPT   =   2
           CALL CC_EIM(WORK(KEMAT1),WORK(KEMAT2),
     *                 R0IM,RBIM,G0IM,GBIM,ONEH,ONEH,
     *                 XLAMDH0,XLAMDP0,ISYM0,
     *                 XLAMDQH,XLAMDQP,ISYHOP,
     *                 FCKCON,LRCON,LGCON,LRELAX,IOPT,ISYHOP)

      ELSE IF ( CC2 .AND. (LTWOEL.OR.LRELAX) ) THEN
           LRCON  = .TRUE.  ! include contrib. from R 
           LGCON  = .TRUE.  ! include contrib. from G
           FCKCON = .TRUE.  ! include contrib. from Fock
           IOPT   =   2
           CALL CC_EIM(WORK(KEMAT1),WORK(KEMAT2),
     *                 R0IM,RBIM,G0IM,GBIM,FOCK,ONEH,
     *                 XLAMDH0,XLAMDP0,ISYM0,
     *                 XLAMDQH,XLAMDQP,ISYHOP,
     *                 FCKCON,LRCON,LGCON,LRELAX,IOPT,ISYHOP)

      ELSE 
           LRCON  = .FALSE. ! skip contrib. from R 
           LGCON  = .FALSE. ! skip contrib. from G
           FCKCON = .TRUE.  ! include contrib. from Fock
           IOPT   =   2
           CALL CC_EIM(WORK(KEMAT1),WORK(KEMAT2),
     *                 R0IM,RBIM,G0IM,GBIM,FOCK,FOCK,
     *                 XLAMDH0,XLAMDP0,ISYM0,
     *                 XLAMDQH,XLAMDQP,ISYHOP,
     *                 FCKCON,LRCON,LGCON,LRELAX,IOPT,ISYHOP)

      END IF

      IOPT = 2
      CALL CCG_1ITRVO(WORK(KB2CON),ISYETA,WORK(KEMAT2),WORK(KEMAT1),
     *                ISYHOP,ZETA1,ISYCTR,IOPT)
      CALL DAXPY(NT1AM(ISYETA),ONE,WORK(KB2CON),1,RHO1,1)
 
      IF (LOCDBG) THEN
          CALL AROUND('E intermediates:')
          CALL CC_PREI(WORK(KEMAT1),WORK(KEMAT2),ISYHOP,1)
          WRITE (LUPRI,*) 'RHO1 after B2 contribution:'
          CALL CC_PRP(RHO1,RHO2,ISYRES,1,0)
          XNORM = DDOT(NT1AM(ISYRES),RHO1,1,RHO1,1)
          WRITE (LUPRI,*) 'Norm of RHO1 after B2 contribution:',XNORM
          CALL FLSHFO(LUPRI)
      END IF

*----------------------------------------------------------------------*
*     For relaxed CC2 replace the E-intermediates by SCF Fock matrix,
*     so that in the following step the correct eta^G_aibj is obtained:
*----------------------------------------------------------------------*
      IF (CC2 .AND. (LTWOEL.OR.LRELAX) ) THEN

         IOPT = 1
         CALL CC_EIM(WORK(KEMAT1),WORK(KEMAT2),DUMMY,DUMMY,DUMMY,DUMMY,
     &             FOCKHF,FOCKHF,DUMMY,DUMMY,IDUMMY,DUMMY,DUMMY,IDUMMY,
     &            .TRUE.,.FALSE.,.FALSE.,.FALSE.,IOPT,ISYHOP)
 
         IF (NFIELD .NE. 0) THEN
           CALL QUIT('relaxed CC2 with finite field not yet finished.')
         END IF

      END IF

*----------------------------------------------------------------------*
*     for non-orbital-relaxed one-electron perturbations calculate
*     here the contribution of the E intermediates to the RHO2
*     (for relaxed CC2 this becomes the eta^G_aibj contribution)
*----------------------------------------------------------------------*
      IF ( (.NOT.(LTWOEL.OR.LRELAX)) .AND. (.NOT.CCS) .OR. CC2) THEN

         KCTR2  = KEND3
         KZETA2 = KCTR2  + NT2AM(ISYCTR)
         KEND3  = KZETA2 + NT2SQ(ISYCTR)
         LWRK3  = LWORK  - KEND3
   
         IF (LWRK3 .LT. 0) THEN
             CALL QUIT('Insufficient work space in CC_XIETA. (E)')
         END IF

         IOPT = 2
         CALL CC_RDRSP(LISTL,IDLSTL,ISYCTR,IOPT,MODEL,DUMMY,
     &                 WORK(KCTR2))
         CALL CC_T2SQ(WORK(KCTR2),WORK(KZETA2),ISYCTR)

         CALL CC_EITR(WORK(KEMAT1),WORK(KEMAT2),
     &                WORK(KEND3),LWRK3,ISYHOP)

         CALL CCRHS_E(RHO2,WORK(KZETA2),WORK(KEMAT1),WORK(KEMAT2),
     &                WORK(KEND3),LWRK3,ISYCTR,ISYHOP)

         IF (LOCDBG) THEN
C           WRITE (LUPRI,*) 'E term contribution from CCRHS_E:'
C           CALL CC_PRP(RHO1,RHO2,ISYETA,0,1)
            XNORM = DDOT(NT2AM(ISYRES),RHO2,1,RHO2,1)
            WRITE (LUPRI,*) 'Norm of Rho2 after CCRHS_E:',XNORM
         END IF

      END IF

*----------------------------------------------------------------------*
*     calculate D2 contributions:
*----------------------------------------------------------------------*
      IF (.NOT. CCS) THEN

         IF (LOCDBG) THEN
           WRITE (LUPRI,'(//A)') 'CCETA2> X-intermediate:'
           WRITE (LUPRI,'(5G15.6)') (XINT(I),I=1,NMATIJ(ISYCTR))
           WRITE (LUPRI,'(//A)') 'CCETA2> Y-intermediate:'
           WRITE (LUPRI,'(5G15.6)') (YINT(I),I=1,NMATAB(ISYCTR))
         END IF

         IF (CC2 .AND. LRELAX .AND. (.NOT.NONHF)) THEN
           CONTINUE
         ELSE IF (CC2 .AND. LRELAX .AND. NONHF) THEN
           CALL QUIT('relaxed CC2 with finite field not yet completed.')
         ELSE ! ccsd or unrelaxed cc2
           CALL CC_21EFM(RHO1,ONEH,ISYHOP,XINT,YINT,ISYCTR,WORK,LWORK)
         END IF

         IF (LOCDBG) THEN
           WRITE (LUPRI,*) 'D2 term contribution from CC_21EFM:'
           CALL CC_PRP(RHO1,WORK,ISYETA,1,0)
           XNORM = DDOT(NT1AM(ISYRES),RHO1,1,RHO1,1)
           WRITE (LUPRI,*) 'Norm of RHO1 after D2 contribution:',XNORM
         END IF

      END IF

*----------------------------------------------------------------------*
*     calculate BFZeta-term and GZeta intermediate:
*----------------------------------------------------------------------*
      IF ( (.NOT. (CCS .OR. CC2)) .AND. (LTWOEL.OR.LRELAX)) THEN

          CALL DZERO(GBIM,NT1AO(ISYETA))
          CALL DZERO(G0IM,NT1AO(ISYCTR))

C         ---------------------------------------------------
C         allocate work space for BFZ intermediates:
C         ---------------------------------------------------
          LENBFZ0 = 2 * NT2ORT(ISYCTR)
          LENBFZ1 = 2 * NT2ORT(ISYETA)

          KBFZ   = 1
          KEND2  = KBFZ  + MAX(LENBFZ0,LENBFZ1)
          LWRK2  = LWORK - KEND2

          IF (LWRK2 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CCETA2. (1)')
          END IF


* read & transform BFZ intermediate to MO representation using 
* XLAMDH0 matrices:
          LGAMMA  = .FALSE.
          LO3BF   = .FALSE.
          LBFZETA = .FALSE.

          CALL GETWA2(LUBFZI,FNBFZI,WORK(KBFZ),IADRE1,LENBFZ1)
          IF (LOCDBG) THEN
             XNORM = DDOT(LENBFZ1,WORK(KBFZ),1,WORK(KBFZ),1)
             WRITE (LUPRI,*) 'norm^2 of BFZ1 intermediate:',XNORM
          END IF


          ICON   = 1
          IOPTG  = 1
          ISYMBF = ISYETA
          CALL CC_T2MO3(DUMMY,DUMMY,1,WORK(KBFZ),RHO2,
     *                  DUMMY,GBIM,DUMMY,
     *                  XLAMDH0,ISYM0,XLAMDH0,ISYM0,WORK(KEND2),LWRK2,
     *                  ISYMBF,ICON,LGAMMA,IOPTG,LO3BF,LBFZETA)
   
          IF ( LRELAX ) THEN

             CALL GETWA2(LUBFZI,FNBFZI,WORK(KBFZ),IADRE0,LENBFZ0)

             IF (LOCDBG) THEN
                XNORM = DDOT(LENBFZ0,WORK(KBFZ),1,WORK(KBFZ),1)
                WRITE (LUPRI,*) 'norm^2 of BFZ0 intermediate:',XNORM
             END IF

             ICON   = 4
             ISYMBF = ISYCTR
             CALL CC_T2MO3(DUMMY,DUMMY,1,WORK(KBFZ),RHO2,
     *                     DUMMY,G0IM,GBIM,
     *                     XLAMDH0,ISYM0,XLAMDQH,ISYHOP,
     *                     WORK(KEND2),LWRK2,ISYMBF,
     *                     ICON,LGAMMA,IOPTG,LO3BF,LBFZETA)

          END IF
   
           
          IF (LOCDBG) THEN
            WRITE (LUPRI,*) 'CCETA2> RHO after B-term section:'
            CALL CC_PRP(RHO1,RHO2,ISYETA,1,1)
            CALL FLSHFO(LUPRI)
          END IF

      END IF

*----------------------------------------------------------------------*
*     calculate E1 term from GZeta intermediate:
*     (note that for relaxed CC2 the GZeta is stored in F0IM/FBIM 
*      and we calculate then also the E2 term from RZ0I/RZBI and the
*      A1 term from ZFCKB)
*----------------------------------------------------------------------*
      IF ( (LTWOEL.OR.LRELAX) .AND. .NOT.(CCS.OR.CC2) ) THEN

           LGCON  = .TRUE.
           LRCON  = .FALSE.
           FCKCON = .FALSE.
           LTRSPF = .FALSE.
           IOPT   =   2
           CALL CC_GHJ(RHO1,DUMMY,DUMMY,G0IM,GBIM,DUMMY,
     *                 XLAMDP0,XLAMDH0,ISYM0,
     *                 XLAMDQP,XLAMDQH,ISYHOP,
     *                 FCKCON,LRCON,LGCON,LRELAX,LTRSPF,IOPT,ISYETA)
C          note that use of hole & particle matrices is exchanged
C          compared to the call in the XI vector part...

C       ... the above should be changed to:
C       CALL CC_T1AM(RHO1,ISYETA,G0IM,ISYCTR,XLAMDQH,ISYHOP,ONE)
C       CALL CC_T1AM(RHO1,ISYETA,GBIM,ISYETA,XLAMDH0,ISYM0,ONE)

      ELSE IF ( (LTWOEL.OR.LRELAX) .AND. CC2 ) THEN

c       XNORM = DNRM2(NEMAT1(ISYCTR),RZ0I,1)
c       WRITE (LUPRI,*) 'norm(rz0i):',xnorm
c       XNORM = DNRM2(NEMAT1(ISYETA),RZBI,1)
c       WRITE (LUPRI,*) 'norm(rzbi):',xnorm
c       XNORM = DNRM2(NT1AO(ISYCTR),F0IM,1)
c       WRITE (LUPRI,*) 'norm(f0im):',xnorm
c       XNORM = DNRM2(NT1AO(ISYETA),FBIM,1)
c       WRITE (LUPRI,*) 'norm(fbim):',xnorm
c       XNORM = DNRM2(N2BST(ISYETA),ZFCKB,1)
c       WRITE (LUPRI,*) 'norm(zfckb):',xnorm

        LGCON  = .TRUE.
        LRCON  = .TRUE.
        FCKCON = .TRUE.
        LTRSPF = .TRUE.
        IOPT   =   2
        CALL CC_GHJ(RHO1,RZ0I,RZBI,F0IM,FBIM,ZFCKB,
     *              XLAMDP0,XLAMDH0,ISYM0,
     *              XLAMDQP,XLAMDQH,ISYHOP,
     *              FCKCON,LRCON,LGCON,LRELAX,LTRSPF,IOPT,ISYETA)
C        note that use of hole & particle matrices is exchanged
C        compared to the call in the XI vector part...

      ELSE IF ( (LTWOEL.OR.LRELAX) .AND. CCS ) THEN

        LGCON  = .FALSE.
        LRCON  = .FALSE.
        FCKCON = .TRUE.
        LTRSPF = .TRUE.
        IOPT   =   2
        CALL CC_GHJ(RHO1,RZ0I,RZBI,F0IM,FBIM,ZFCKB,
     *              XLAMDP0,XLAMDH0,ISYM0,
     *              XLAMDQP,XLAMDQH,ISYHOP,
     *              FCKCON,LRCON,LGCON,LRELAX,LTRSPF,IOPT,ISYETA)
C        note that use of hole & particle matrices is exchanged
C        compared to the call in the XI vector part...

      END IF

      IF (LOCDBG .AND. (LTWOEL.OR.LRELAX)) THEN
        WRITE (LUPRI,*) 'Finished E1-term section... RHO1:'
        CALL CC_PRP(RHO1,WORK,ISYRES,1,0)
        XNORM = DDOT(NT1AM(ISYRES),RHO1,1,RHO1,1)
        WRITE (LUPRI,*) 'Norm of RHO1 after E1 term:',XNORM
        CALL FLSHFO(LUPRI)
      ENDIF

*---------------------------------------------------------------------*
* allocate work space for packed ZETA amplitudes and read them in core:
*---------------------------------------------------------------------*
      IF ( (.NOT. (CCS .OR. CC2)) .AND. (LTWOEL.OR.LRELAX)) THEN

          KCTR2 = 1
          KEND2 = KCTR2 + NT2AM(ISYCTR)
          LWRK2 = LWORK - KEND2

          IF (LWRK2 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CC_XIETA. (CTR2)')
          END IF

* read lagrange multipliers from file:
          IOPT = 2
          CALL CC_RDRSP(LISTL,IDLSTL,ISYCTR,IOPT,MODEL,DUMMY,
     &                  WORK(KCTR2))

          IF (LOCDBG) THEN
             WRITE (LUPRI,*) 'KCTR2:',KCTR2
             WRITE (LUPRI,*) 'ISYCTR:',ISYCTR
             WRITE (LUPRI,*) 'CTR2:'
             CALL CC_PRP(WORK,WORK(KCTR2),ISYCTR,0,1)
             CALL FLSHFO(LUPRI)
          END IF

      END IF

*---------------------------------------------------------------------*
*     calculate A-term contribution:
*       requires:  -- Gamma intermediate in core
*                  -- ZETA2 amplitudes squared in core
*---------------------------------------------------------------------*
      IF ( (.NOT. (CCS .OR. CC2)) .AND. (LTWOEL.OR.LRELAX)) THEN

* square ZETA2 amplitudes up (overwritting T2AMP0 array):
          KZETA2 = KEND2
          KEND3  = KZETA2 + NT2SQ(ISYCTR)
          LWRK3 = LWORK - KEND3

          IF (LWRK3 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CCETA2. (A)')
          END IF

          CALL CC_T2SQ(WORK(KCTR2),WORK(KZETA2),ISYCTR)

* calculate the A-term:
          IOPT = 2
          CALL CCRHS_A(RHO2,WORK(KZETA2),XGAMMA,
     *                 WORK(KEND3),LWRK3,ISYHOP,ISYCTR,IOPT)

          IF (LOCDBG) THEN
             WRITE (LUPRI,*) 'CC_XIETA> after A-term section:'
             WRITE (LUPRI,*) 'RHO:'
             CALL CC_PRP(RHO1,RHO2,ISYRES,1,1)
             WRITE (LUPRI,*) 'ZETA2:'
             CALL CC_PRSQ(WORK,WORK(KZETA2),ISYRES,0,1)
             CALL FLSHFO(LUPRI)
             WRITE (LUPRI,*) 'GAMMA:',(XGAMMA(I),I=1,NGAMMA(ISYHOP))
          END IF

      END IF

*---------------------------------------------------------------------*
*       calculate C- and D-term contributions:
*         requires:  -- ZETA amplitudes packed in core
*                    -- C/D intermediate in core (squared matrix)
*                    -- result vector in core (packed)
*       (note that ZETA2 vector is overwritten in this section)
*---------------------------------------------------------------------*
      IF ( (.NOT. (CCS .OR. CC2)) .AND. (LTWOEL.OR.LRELAX)) THEN
       
          KCDBAR  = KEND2
          KEND3   = KCDBAR + NT2SQ(ISYHOP)  
          LWRK3   = LWORK  - KEND3
    
          IF (LWRK3 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CCETA2. (C-term)')
          END IF

C         ------------------------------------------------------------
C         read D intermediate from file:
C         ------------------------------------------------------------
          CALL GETWA2(LUDIM,FNDIM,WORK(KCDBAR),IDRDIM,NT2SQ(ISYHOP))

C         ------------------------------------------------------------
C         calculate D term and add to result vector:
C         ------------------------------------------------------------
          ioptr12 = 0
          CALL CC_CD('D',-1,ioptr12,RHO2,ISYETA,WORK(KCTR2),
     &               ISYCTR,WORK(KCDBAR),ISYHOP,WORK(KEND3),LWRK3)

          IF (LOCDBG) THEN
             WRITE (LUPRI,*) 'Finished D-term section...'
             WRITE (LUPRI,*) 'RHO2:'
             CALL CC_PRP(WORK,RHO2,ISYRES,0,1)
             CALL FLSHFO(LUPRI)
          END IF

C         ------------------------------------------------------------
C         read C intermediate from file:
C         ------------------------------------------------------------
          CALL GETWA2(LUCIM,FNCIM,WORK(KCDBAR),IDRCIM,NT2SQ(ISYHOP))

C         ------------------------------------------
C         calculate C term and add to result vector:
C         ------------------------------------------
          ioptr12 = 0
          CALL CC_CD('C',-1,ioptr12,RHO2,ISYETA,WORK(KCTR2),
     &               ISYCTR,WORK(KCDBAR),ISYHOP,WORK(KEND3),LWRK3)

          IF (LOCDBG) THEN
             WRITE (LUPRI,*) 'Finished C-term section... RHO:'
             CALL CC_PRP(RHO1,RHO2,ISYETA,1,1)
             CALL FLSHFO(LUPRI)
          END IF

      END IF

*---------------------------------------------------------------------*
*       calculate 21F and 21G contributions:
*---------------------------------------------------------------------*
      IF ( (LTWOEL.OR.LRELAX) .AND. .NOT.(CCS .OR. CC2) ) THEN

        CALL CC_T1AM(RHO1,ISYETA,F0IM,ISYCTR,XLAMDQH,ISYHOP,ONE)

        CALL CC_T1AM(RHO1,ISYETA,FBIM,ISYETA,XLAMDH0,ISYM0,ONE)

        IF (LOCDBG) THEN
           WRITE (LUPRI,*) 'Finished 21F term section... RHO1:'
           CALL CC_PRP(RHO1,WORK,ISYETA,1,0)
           CALL FLSHFO(LUPRI)
        END IF

        MT2BCD = 0
        DO ISYM = 1, NSYM
          MT2BCD = MAX(MT2BCD,NT2BCD(ISYM))
        END DO

        ISYIFJL = ISYHOP
        KXIFJL  = KEND2
        KYJLIF  = KXIFJL + NTAIKJ(ISYIFJL)
        KPINT0  = KYJLIF + NTAIKJ(ISYIFJL)
        KQINT0  = KPINT0 + MT2BCD
        KWINT0  = KQINT0 + MT2BCD
        KEND3   = KWINT0 + MT2BCD
        LWRK3   = LWORK  - KEND3

        IF (LWRK3 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CC_XIETA. (21G-term)')
        END IF


C       -------------------------------------------------
C       calculate (if|jl) integrals, stored as X(fi,l,j),
C       and resort to Y(jli,f):
C       -------------------------------------------------
        CALL DZERO(WORK(KXIFJL),NTAIKJ(ISYIFJL))

        IOPT = 2
        LOCC = .TRUE.
        CALL CC_IAJB2(WORK(KXIFJL),ISYIFJL,IOPT,LRELAX,LOCC,.FALSE.,
     &                LU0IAJB,FN0IAJB,IT2DEL0,XLAMDH0,ISYM0,
     &                LU1IAJB,FN1IAJB,IT2DELB,XLAMDQH,ISYHOP,
     &                WORK(KEND3),LWRK3)

        CALL CCG_SORT1(WORK(KXIFJL),WORK(KYJLIF),ISYIFJL,5)

C       --------------------------------------------------------
C       contract with P and Q intermediates to 21G contribution:
C       --------------------------------------------------------
        CALL DZERO(WORK(KWINT0),MT2BCD)

        DO ISYMF = 1, NSYM

           ISYALJ = MULD2H(ISYCTR,ISYMF)
           ISYJLI = MULD2H(ISYIFJL,ISYMF)

           DO F = 1, NVIR(ISYMF)
              IVIRF  = IVIR(ISYMF) + F
              IADR   = IADRPQMO(IVIRF)
              LEN    = NT2BCD(ISYALJ)
              CALL GETWA2(LUPQMO,FILPQMO,WORK(KPINT0),IADR    ,LEN)
              CALL GETWA2(LUPQMO,FILPQMO,WORK(KQINT0),IADR+LEN,LEN)

              KOFF1 = KYJLIF + ISJIKA(ISYJLI,ISYMF) 
     &                       + NMAIJK(ISYJLI)*(F-1)

              CALL CC_21H(RHO1,ISYETA,WORK(KQINT0),
     &                    WORK(KWINT0),WORK(KPINT0),ISYALJ,
     &                    WORK(KOFF1),ISYHOP,WORK(KEND3),LWRK3,ISYMF)
           END DO

        END DO 

        IF (LOCDBG) THEN
           WRITE (LUPRI,*) 'Finished 21G term section... RHO1:'
           CALL CC_PRP(RHO1,WORK,ISYETA,1,1)
           CALL FLSHFO(LUPRI)
        END IF

      END IF

*---------------------------------------------------------------------*
*       calculate the H contribution:
*---------------------------------------------------------------------*
      IF (.NOT. CCS) THEN

         KYPS0  = 1
         KYPS1  = KYPS0  + NGLMDT(ISYCTR)
         KXPCK  = KYPS1  + NGLMDT(ISYETA)
         KXIAJB = KXPCK  + NT2AM(ISYETA)
         KEND3  = KXIAJB + NT2SQ(ISYETA)
         LWRK3  = LWORK  - KEND3
    
         IF (LWRK3 .LT. 0) THEN
           CALL QUIT('Insufficient work space in CC_XIETA. (C-term)')
         END IF
C
C        the CTR1 x Fock contribution:
C
         CALL CC_L1FCK(RHO2,ZETA1,FOCK,ISYCTR,ISYHOP,
     &                 WORK(KEND3),LWRK3)

         IF (LOCDBG) THEN
            WRITE (LUPRI,*) 'after CC_L1FCK contribution to RHO2:'
            CALL CC_PRP(WORK,RHO2,ISYETA,0,1)
            XNORM = DDOT(NT1AM(ISYCTR),ZETA1,1,ZETA1,1)
            WRITE (LUPRI,*) 'NORM^2(ZETA1):',XNORM
            XNORM = DDOT(N2BST(ISYHOP),FOCK,1,FOCK,1)
            WRITE (LUPRI,*) 'NORM^2(FOCK):',XNORM
         END IF

         IF (LTWOEL.OR.LRELAX) THEN
C
C           calculate the Ypsilon matrices:
C           (for CC2 skip contribution from YINT)
C
            IOPT = 1
            IF (CC2) IOPT = 0
            CALL CCLT_YPS(ZETA1,YINT,ISYCTR,XLAMDH0,ISYM0, 
     &                    WORK(KYPS0),IOPT)
            CALL CCLT_YPS(ZETA1,YINT,ISYCTR,XLAMDQH,ISYHOP,
     &                    WORK(KYPS1),IOPT)
C
C           transform the (ia|j del) integrals with Ypsilon matrices:
C
            CALL DZERO(WORK(KXIAJB),NT2SQ(ISYETA))

            IOPT = 2
            LOCC = .FALSE.
            CALL CC_IAJB2(WORK(KXIAJB),ISYETA,IOPT,LRELAX,LOCC,.FALSE.,
     &                    LU0IAJB,FN0IAJB,IT2DEL0,WORK(KYPS0),ISYCTR,
     &                    LU1IAJB,FN1IAJB,IT2DELB,WORK(KYPS1),ISYETA,
     &                    WORK(KEND3),LWRK3)

C
C           for cc2 add sum_c Zeta_cj (ia|cb) contribution:
C
            IF (CC2) THEN
              CALL DSCAL(NT2SQ(ISYETA),-ONE,WORK(KXIAJB),1)

              IOPT = 2
              LOCC = .FALSE.
              CALL CC_IAJB2(WORK(KXIAJB),ISYETA,IOPT,
     &                      LRELAX,LOCC,.FALSE.,
     &                      LUHINT0,FNHINT0,IADRH0,XLAMDH0,ISYM0,
     &                      LUHINT1,FNHINT1,IADRH1,XLAMDQH,ISYHOP,
     &                      WORK(KEND3),LWRK3)
            END IF
C
C           apply P^ab_ij, pack the integrals and add to result vector:
C
            IOPT = 1
            CALL CC_T2PK(WORK(KXPCK),WORK(KXIAJB),ISYETA,1) 
            IOPTTCME = 1
            CALL CCSD_TCMEPK(WORK(KXPCK),ONE,ISYETA,IOPTTCME)

            IF (CC2) THEN
              CALL DAXPY(NT2AM(ISYETA),+ONE,WORK(KXPCK),1,RHO2,1)
            ELSE
              CALL DAXPY(NT2AM(ISYETA),-ONE,WORK(KXPCK),1,RHO2,1)
            END IF

            IF (LOCDBG) THEN
               WRITE (LUPRI,*) 'L(jb|ia^Y) integrals in H-term:'
               CALL CC_PRP(WORK,WORK(KXPCK),ISYETA,0,1)
               WRITE (LUPRI,*) 'RHO after H-term:' 
               CALL CC_PRP(RHO1,RHO2,ISYETA,1,1)
            END IF

         END IF

      END IF

*---------------------------------------------------------------------*
*     That's it; return
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
         I = 1
         IF (CC2) I = 0
         CALL AROUND('final result vector in CCETA2')
         CALL CC_PRP(RHO1,RHO2,ISYETA,1,1) 
      END IF
     
      CALL QEXIT('CCETA2')

      RETURN
      END 
*=====================================================================*
*                      END OF SUBROUTINE CCETA2                       *
*=====================================================================*
*======================================================================
      SUBROUTINE CC_RIM(RIM,T2,ISYMTR,XLDM,ISYLDM,IDEL,ISYDEL)
*----------------------------------------------------------------------
*
*     Purpose: Calculate the R intermediate which is a precursor
*              for one of the E intermediates and for the H contrib. 
*              to the vector function and relaxed right transformations
*
*              RIM_{b,del} += sum_{dlm} T2_{dl,bm} (ld|m del)
*
*              T2_{ai,bj} = 2 T_{ai,bj} - T_{aj,bi}
*              XLDM_{ld,m}  = (ld|m del)
*
*              
*              Symmetries:     ISYMTR  --  T2 amplitudes
*                              ISYLDM  --  XLDM integrals
*                              ISYDEL  --  IDEL delta index
*
*
*     Christof Haettig, May 1998
*
*======================================================================
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"

      INTEGER ISYMTR, ISYDEL, ISYLDM, IDEL

#if defined (SYS_CRAY)
      REAL ONE 
      REAL RIM(*), T2(*), XLDM(*)
#else
      DOUBLE PRECISION ONE
      DOUBLE PRECISION RIM(*), T2(*), XLDM(*)
#endif
      PARAMETER(ONE=1.0D0)

      INTEGER ISYMM,ISYMDL,ISYMBM,ISYMB,NT1DL,KOFF1,KOFF2,KOFF3,KBM
 
      CALL QENTER('CC_RIM')

 
*----------------------------------------------------------------------
*     Contract the integrals in XLDM with T2 amplitudes.
*----------------------------------------------------------------------
 
      DO ISYMM = 1,NSYM
 
         ISYMDL = MULD2H(ISYMM,ISYLDM)
         ISYMBM = MULD2H(ISYMDL,ISYMTR)
         ISYMB  = MULD2H(ISYMBM,ISYMM)
 
         DO M = 1,NRHF(ISYMM)
 
            D = IDEL - IBAS(ISYDEL) 
 
            NT1DL = MAX(NT1AM(ISYMDL),1)

            KBM   = IT1AM(ISYMB,ISYMM) + NVIR(ISYMB)*(M - 1) + 1

            KOFF1 = IT2SQ(ISYMDL,ISYMBM) + NT1AM(ISYMDL)*(KBM - 1) + 1
            KOFF2 = IT2BCD(ISYMDL,ISYMM) + NT1AM(ISYMDL)*(M - 1) + 1
            KOFF3 = IEMAT1(ISYMB,ISYDEL) + (D - 1)*NVIR(ISYMB) + 1
 
            CALL DGEMV('T',NT1AM(ISYMDL),NVIR(ISYMB),ONE,T2(KOFF1),
     &                 NT1DL,XLDM(KOFF2),1,ONE,RIM(KOFF3),1)
 
         END DO
 
      END DO
 
      CALL QEXIT('CC_RIM')

      RETURN
      END
*======================================================================
*=====================================================================*
      SUBROUTINE CC_GHJ(OMEGA1,RAIM,RBIM,GAIM,GBIM,ONEHAM,
     *                  XLAMDHA,XLAMDPA,ISYLMA,XLAMDHB,XLAMDPB,ISYLMB,
     *                  ONHCON,LRCON,LGCON,LRELAX,TRSPOH,IOPT,ISYRES)
*---------------------------------------------------------------------*
*
*     Calculate G-, H- and J-terms from R and G intermediates
*     and the MO one-electron hamiltonian integrals
*
*     IOPT = 1:     OMEGA1 = ONEHAM + XLAMDHA * RAIM + XLAMDPA * GAIM
*                   
*                   RBIM,GBIM,XLAMDHB,XLAMDPB,ISYMB are dummy input
*
*
*     IOPT = 2:     OMEGA1 = ONEHAM + XLAMDHA * RBIM + XLAMDHB * RAIM
*                                   + XLAMDPA * GBIM + XLAMDPB * GAIM
*
*     ONHCON = .FALSE.  :  skip contribution from ONEHAM
*     LRCON  = .FALSE.  :  skip contribution from R intermediates
*     LGCON  = .FALSE.  :  skip contribution from G intermediates
*
*     LRELAX = .FALSE.  :  skip relaxation contributions from 
*                          XLAMDPB, XLAMDHB
*     TRSPOH = .TRUE.   :  add transposed ONEHAM matrix
*
*     Symmetries:    ISYMRES  --  OMEGA1,  ONEHAM 
*                    ISYLMA   --  XLAMDHA, XLAMDPA   
*                    ISYLMB   --  XLAMDHB, XLAMDPB   
*
*     Christof Haettig 20-6-1998
*
*---------------------------------------------------------------------*
C
#include "implicit.h"
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
      DIMENSION OMEGA1(*),ONEHAM(*)
      DIMENSION RAIM(*),RBIM(*),GAIM(*),GBIM(*)
      DIMENSION XLAMDHA(*),XLAMDHB(*),XLAMDPA(*),XLAMDPB(*)
      LOGICAL ONHCON, LRCON, LGCON, LRELAX, TRSPOH
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
      INTEGER IOPT
C
      CALL QENTER('CC_GHJ')
C
C---------------------------------------------------------
C     Transform the AO index of R intermediate(s) to i.
C---------------------------------------------------------
C
      IF (LRCON) THEN
C
         DO ISYMD = 1,NSYM
C
            ISYMI = MULD2H(ISYMD,ISYLMA)
            ISYMA = MULD2H(ISYMI,ISYRES)
C
            NVIRA = MAX(NVIR(ISYMA),1)
            NBASD = MAX(NBAS(ISYMD),1)
C
            KOFF1 = IEMAT1(ISYMA,ISYMD) + 1
            KOFF2 = IGLMRH(ISYMD,ISYMI) + 1
            KOFF3 = IT1AM(ISYMA,ISYMI)  + 1
C
            IF ( IOPT .EQ. 1) THEN

              CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMD),
     *                   -ONE,RAIM(KOFF1),NVIRA,XLAMDHA(KOFF2),NBASD,
     *                   ONE,OMEGA1(KOFF3),NVIRA)

            ELSE

              CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMD),
     *                   -ONE,RBIM(KOFF1),NVIRA,XLAMDHA(KOFF2),NBASD,
     *                   ONE,OMEGA1(KOFF3),NVIRA)

              IF (LRELAX) THEN

                ISYMI = MULD2H(ISYMD,ISYLMB)
                ISYMA = MULD2H(ISYMI,ISYRES)

                NVIRA = MAX(NVIR(ISYMA),1)
                NBASD = MAX(NBAS(ISYMD),1)
C
                KOFF1 = IEMAT1(ISYMA,ISYMD) + 1
                KOFF2 = IGLMRH(ISYMD,ISYMI) + 1
                KOFF3 = IT1AM(ISYMA,ISYMI)  + 1
C
                CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMD),
     *                     -ONE,RAIM(KOFF1),NVIRA,XLAMDHB(KOFF2),NBASD,
     *                     ONE,OMEGA1(KOFF3),NVIRA)
              END IF

            END IF
C
         END DO
C
      END IF 
C
C
C---------------------------------------------------------
C     Transform the AO index of G intermediate(s) to a.
C---------------------------------------------------------
C
      IF (LGCON) THEN
C
         DO ISYMD = 1,NSYM
C
            ISYMA = MULD2H(ISYMD,ISYLMA)
            ISYMI = MULD2H(ISYMA,ISYRES)
C
            NVIRA = MAX(NVIR(ISYMA),1)
            NBASD = MAX(NBAS(ISYMD),1)
C
            KOFF1 = IT1AO(ISYMD,ISYMI)  + 1
            KOFF2 = IGLMVI(ISYMD,ISYMA) + 1
            KOFF3 = IT1AM(ISYMA,ISYMI) + 1
C
            IF ( IOPT .EQ. 1) THEN

              CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMD),
     *                   ONE,XLAMDPA(KOFF2),NBASD,GAIM(KOFF1),NBASD,
     *                   ONE,OMEGA1(KOFF3),NVIRA)

            ELSE

              CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMD),
     *                   ONE,XLAMDPA(KOFF2),NBASD,GBIM(KOFF1),NBASD,
     *                   ONE,OMEGA1(KOFF3),NVIRA)
C
              IF (LRELAX) THEN
C
                ISYMA = MULD2H(ISYMD,ISYLMB)
                ISYMI = MULD2H(ISYMA,ISYRES)
C
                NVIRA = MAX(NVIR(ISYMA),1)
                NBASD = MAX(NBAS(ISYMD),1)
C
                KOFF1 = IT1AO(ISYMD,ISYMI)  + 1
                KOFF2 = IGLMVI(ISYMD,ISYMA) + 1
                KOFF3 = IT1AM(ISYMA,ISYMI)  + 1
C
                CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMD),
     *                     ONE,XLAMDPB(KOFF2),NBASD,GAIM(KOFF1),NBASD,
     *                     ONE,OMEGA1(KOFF3),NVIRA)
              END IF
C
            END IF
C
         END DO
C
      END IF 
C
C--------------------------------------
C     Add the Fock matrix contribution:
C--------------------------------------
C
      IF (ONHCON) THEN
C
       IF (TRSPOH) THEN
C
         ! omega1_ai = omega1_ai + h_ia
         DO ISYMI = 1,NSYM
           ISYMA = MULD2H(ISYMI,ISYRES)
           DO I = 1,NRHF(ISYMI)
             DO A = 1,NVIR(ISYMA)
               KOFF1 = IT1AM(ISYMA,ISYMI) +NVIR(ISYMA)*(I-1)+A
               KOFF2 = IFCVIR(ISYMI,ISYMA)+NORB(ISYMI)*(A-1)+I
               OMEGA1(KOFF1) = OMEGA1(KOFF1) + ONEHAM(KOFF2)
             END DO
           END DO
         END DO
C
       ELSE
C
         ! omega1_ai = omega1_ai + h_ai
         DO ISYMI = 1,NSYM
           ISYMA = MULD2H(ISYMI,ISYRES)
           DO I = 1,NRHF(ISYMI)
            KOFF1 = IT1AM(ISYMA,ISYMI) +NVIR(ISYMA)*(I-1)+1
            KOFF2 = IFCRHF(ISYMA,ISYMI)+NORB(ISYMA)*(I-1)+NRHF(ISYMA)+1
            CALL DAXPY(NVIR(ISYMA),ONE,ONEHAM(KOFF2),1,OMEGA1(KOFF1),1)
           END DO
         END DO
C
       END IF
C
      END IF
C
      CALL QEXIT('CC_GHJ')
C
      RETURN
      END
*======================================================================*
C  /* Deck cclt_yps */
      SUBROUTINE CCLT_YPS(CTR1,YI,ISYCTR,XLAMDH,ISYLAM,YPS,IOPTY)
C
C     Purpose: To calculate the Ypsilon intermediate:
C
C     Yps(alpha a)  =   sum_k XLAMDH(alpha k) CTR1(a k) 
C                     + sum_f XLAMDH(alpha f)   YI(f a) 
C
C     if IOPTY <> 1 the contribution from YI is skipped
C
C     ISYCTR : symmetry of CTR1, XI
C     ISYLAM : symmetry of XLAMDH
C
C     Christof Haettig, October 1998
C
#include "implicit.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION CTR1(*),XLAMDH(*),YPS(*),YI(*)
C
      CALL QENTER('CCLT_YPS')
C
C---------------------------------------------------
C     Half-transformation to AO-basis of CTR1 and YI
C---------------------------------------------------
C
      ISYMAO = MULD2H(ISYCTR,ISYLAM)
C
      CALL DZERO(YPS,NGLMDT(ISYMAO))
C
      DO ISYMAL = 1,NSYM
C
         ISYMA = MULD2H(ISYMAL,ISYMAO)
         ISYMK = MULD2H(ISYMA,ISYCTR)
         ISYMF = MULD2H(ISYMA,ISYCTR)
C
         KOFF1 = IGLMRH(ISYMAL,ISYMK) + 1
         KOFF2 = IT1AM(ISYMA,ISYMK)   + 1
         KOFF3 = IGLMVI(ISYMAL,ISYMA) + 1
C
         NTOTBA = MAX(NBAS(ISYMAL),1)
         NTOTVI = MAX(NVIR(ISYMA),1)
C
         CALL DGEMM('N','T',NBAS(ISYMAL),NVIR(ISYMA),NRHF(ISYMK),
     *               ONE,XLAMDH(KOFF1),NTOTBA,CTR1(KOFF2),NTOTVI,
     *               ONE,YPS(KOFF3),NTOTBA)
C
         IF (IOPTY.EQ.1) THEN
            KOFF4 = IMATAB(ISYMF,ISYMA)  + 1
            KOFF5 = IGLMVI(ISYMAL,ISYMF) + 1
C
            NTOTVI = MAX(NVIR(ISYMF),1)
C
            CALL DGEMM('N','N',NBAS(ISYMAL),NVIR(ISYMA),NVIR(ISYMF),
     *                  ONE,XLAMDH(KOFF5),NTOTBA,YI(KOFF4),NTOTVI,
     *                  ONE,YPS(KOFF3),NTOTBA)
         END IF
C
      END DO
C
      CALL QEXIT('CCLT_YPS')
C
      RETURN
      END
*======================================================================
*======================================================================
      SUBROUTINE CC_FCKRLX(FOCK0,FOCKD,FOCK1,FOCK2,
     &                     XLAMDP0,XLAMDH0,ISYM0,
     &                     XLAMDPD,XLAMDHD,ISYMD,LRELAX,
     &                     XLAMDP1,XLAMDH1,ISYM1,LRELAX1,
     &                     XLAMDP2,XLAMDH2,ISYM2,LRELAX2,
     &                     IORDER,WORK,LWORK)
*----------------------------------------------------------------------
*
*     Purpose: transform derivative AO fock matrix to MO basis and
*              add relaxation contributions coming from the 
*              derivatives of the transformation matrices
* 
*      if a LRELAX flag is false, skip corresp. relaxation contrib.
*
*      FOCKD : derivative fock matrix, replaced on output
*      FOCK0 : zeroth-order fock matrix, unchanged on output
*      FOCK1 : first-order fock for operator 1, unchanged on output
*      FOCK2 : first-order fock for operator 2, unchanged on output
*
*     Christof Haettig, July 1998
*     second-order stuff for one operator, Christof Haettig, June 1999
*
*======================================================================
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      LOGICAL LRELAX, LRELAX1, LRELAX2
      INTEGER ISYM0, ISYMD, ISYM1, ISYM2
      INTEGER LWORK, ISYRES, KEND1, KSCR, LWRK1, IORDER

#if defined (SYS_CRAY)
      REAL FOCKD(*), FOCK0(*), FOCK1(*), FOCK2(*), WORK(*)
      REAL XLAMDP0(*), XLAMDH0(*), XLAMDPD(*), XLAMDHD(*)
      REAL XLAMDP1(*), XLAMDH1(*), XLAMDP2(*), XLAMDH2(*)
      REAL ONE, XNORM, DDOT
#else
      DOUBLE PRECISION FOCKD(*), FOCK0(*), FOCK1(*), FOCK2(*), WORK(*)
      DOUBLE PRECISION XLAMDP0(*), XLAMDH0(*), XLAMDPD(*), XLAMDHD(*)
      DOUBLE PRECISION XLAMDP1(*), XLAMDH1(*), XLAMDP2(*), XLAMDH2(*)
      DOUBLE PRECISION ONE, XNORM, DDOT
#endif
      PARAMETER(ONE=1.0D0)

      CALL QENTER('CC_FCKRLX')

 
*---------------------------------------------------------------------*
*       if debug flag set, print input matrices in AO:
*---------------------------------------------------------------------*
        IF (LOCDBG) THEN
          WRITE (LUPRI,*) 'CC_FCKRLX> IORDER,LRELAX:',IORDER,LRELAX
          WRITE (LUPRI,*) 'CC_FCKRLX> LRELAX1,LRELAX2:',LRELAX1,LRELAX2
          WRITE (LUPRI,*) 'CC_FCKRLX> FOCKD in AO:'
          CALL CC_PRFCKAO(FOCKD,ISYMD)
          XNORM = DDOT(N2BST(ISYMD),FOCKD,1,FOCKD,1)
          WRITE (LUPRI,*) 'Norm of AO FOCKD matrix:', XNORM
          WRITE (LUPRI,*) 'CC_FCKRLX> FOCK0 in AO:'
          CALL CC_PRFCKAO(FOCK0,ISYM0)
          XNORM = DDOT(N2BST(ISYM0),FOCK0,1,FOCK0,1)
          WRITE (LUPRI,*) 'Norm of AO FOCK0 matrix:', XNORM
          CALL FLSHFO(LUPRI)
        END IF

*---------------------------------------------------------------------*
*       transform derivative AO Fock matrix to MO using XLAMDP0/XLAMDH0
*---------------------------------------------------------------------*
        ISYRES = MULD2H(ISYMD,ISYM0)

        CALL CC_FCKMO(FOCKD,XLAMDP0,XLAMDH0,
     *                WORK,LWORK,ISYMD,ISYM0,ISYM0)

*---------------------------------------------------------------------*
*       transform zeroth-order Fock matrix to MO using XLAMDPD/XLAMDHD
*---------------------------------------------------------------------*
        IF (LRELAX) THEN
 
          CALL CC_FCKMO2(FOCKD,FOCK0,ISYM0,XLAMDP0,XLAMDH0,ISYM0,
     &                   XLAMDPD,XLAMDHD,ISYMD,WORK,LWORK)

        END IF

*---------------------------------------------------------------------*
*       add additional contributions for IORDER = 2:
*---------------------------------------------------------------------*
        IF (IORDER.EQ.2) THEN

          IF (LRELAX2) THEN
            CALL CC_FCKMO2(FOCKD,FOCK1,ISYM1,XLAMDP0,XLAMDH0,ISYM0,
     &                     XLAMDP2,XLAMDH2,ISYM2,WORK,LWORK)
          END IF

          IF (LRELAX1) THEN
            CALL CC_FCKMO2(FOCKD,FOCK2,ISYM2,XLAMDP0,XLAMDH0,ISYM0,
     &                     XLAMDP1,XLAMDH1,ISYM1,WORK,LWORK)
          END IF

          IF (LRELAX1 .AND. LRELAX2) THEN
            CALL CC_FCKMO2(FOCKD,FOCK0,ISYM0,XLAMDP1,XLAMDH1,ISYM1,
     &                     XLAMDP2,XLAMDH2,ISYM2,WORK,LWORK)
          END IF

        END IF


*---------------------------------------------------------------------*
*       print debug output and return:
*---------------------------------------------------------------------*

        IF (LOCDBG) THEN
          WRITE (LUPRI,*) 'CC_FCKRLX> IORDER:',IORDER
          WRITE (LUPRI,*) 'CC_FCKRLX> LRELAX, LRELAX1, LRELAX2:',
     &                        LRELAX, LRELAX1, LRELAX2       
          WRITE (LUPRI,*) 'CC_FCKRLX> FOCKD in MO:'
          CALL CC_PRFCKMO(FOCKD,ISYRES)
          XNORM = DDOT(N2BST(ISYRES),FOCKD,1,FOCKD,1)
          WRITE (LUPRI,*) 'Norm of MO FOCKD matrix:', XNORM
          CALL FLSHFO(LUPRI)
        END IF

      CALL QEXIT('CC_FCKRLX')

      RETURN

      END

*======================================================================
*======================================================================
      SUBROUTINE CC_FCKMO2(FOCKMO,FOCKAO,ISYFAO,
     &                     XLAMDP1,XLAMDH1,ISYM1,
     &                     XLAMDP2,XLAMDH2,ISYM2,
     &                     WORK,LWORK)
*----------------------------------------------------------------------
*
*     Purpose: transform a AO fock matrix to MO basis using two
*              two different sets of XLAMDA matrices as:
*
*     FOCKMO = FOCKMO + XLP1 x FOCKAO x XLH2 + XLP2 x FOCKAO x XLH1
*
*     Note, that FOCKMO is not initialized here!
* 
*     Christof Haettig, June 1999
*
*======================================================================
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER ISYFAO, ISYM2, ISYM1, LWORK, ISYRES, KEND1, KSCR, LWRK1

#if defined (SYS_CRAY)
      REAL FOCKMO(*), FOCKAO(*), WORK(*)
      REAL XLAMDP2(*), XLAMDH2(*), XLAMDP1(*), XLAMDH1(*)
      REAL ONE, XNORM, DDOT
#else
      DOUBLE PRECISION FOCKMO(*), FOCKAO(*), WORK(*)
      DOUBLE PRECISION XLAMDP2(*), XLAMDH2(*), XLAMDP1(*), XLAMDH1(*)
      DOUBLE PRECISION ONE, XNORM, DDOT
#endif
      PARAMETER(ONE=1.0D0)

      CALL QENTER('CC_FCKMO2')

*---------------------------------------------------------------------*
*       transform derivative AO Fock matrix to MO:
*---------------------------------------------------------------------*
      ISYRES = MULD2H(ISYFAO,MULD2H(ISYM1,ISYM2))

      KSCR  = 1
      KEND1 = KSCR + MAX(N2BST(ISYFAO),N2BST(ISYRES))
      LWRK1 = LWORK - KEND1

      IF ( LWRK1 .LT. 0 ) THEN
        CALL QUIT('Insufficient work space in CC_FCKMO2.')
      END IF

*     duplicate AO Fock matrix in WORK:
      CALL DCOPY(N2BST(ISYFAO),FOCKAO,1,WORK(KSCR),1)

*     transform AO FOCK with XLAMDP1 and XLAMDH2,
*     and add to transformed derivative Fock matrix:
      CALL CC_FCKMO(WORK(KSCR),XLAMDP1,XLAMDH2,
     &                WORK(KEND1),LWRK1,ISYFAO,ISYM1,ISYM2)

      CALL DAXPY(N2BST(ISYRES),ONE,WORK(KSCR),1,FOCKMO,1)

*     create a new duplicate of the AO Fock matrix in WORK:
      CALL DCOPY(N2BST(ISYFAO),FOCKAO,1,WORK(KSCR),1)

*     transform zeroth-order AO FOCK with XLAMDP2 and XLAMDH1,
*     and add to transformed derivative Fock matrix:
      CALL CC_FCKMO(WORK(KSCR),XLAMDP2,XLAMDH1,
     &              WORK(KEND1),LWRK1,ISYFAO,ISYM2,ISYM1)

      CALL DAXPY(N2BST(ISYRES),ONE,WORK(KSCR),1,FOCKMO,1)

      CALL QEXIT('CC_FCKMO2')

      RETURN

      END

*======================================================================
