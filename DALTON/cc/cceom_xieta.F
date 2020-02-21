!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2015 (2015), see http://daltonprogram.org"
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
C
*=====================================================================*
       SUBROUTINE CCEOM_XIETA( IXETRAN, NXETRAN, IOPTRES, IORDER, LISTL,
     &                      FILXI,   IXDOTS, XCONS, 
     &                      FILETA,  IEDOTS, ECONS, 
     &                      MXVEC,   WORK,    LWORK         )
*---------------------------------------------------------------------*
*
*    Purpose: batched loop over Xi and Eta vector calculations
*    for EOM/CCCI
*
*             for more detailed documentation see: CC_XIETA1
*        
*     Written by Sonia and Filip. 2015,
*     Based on CC_XIETA by Christof Haettig
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
!sonia
#include "second.h"

      LOGICAL LOCDBG
      PARAMETER (LOCDBG=.false.)
 
      INTEGER MAXXETRAN, MXVEC
      PARAMETER (MAXXETRAN = 100)

      CHARACTER*(*) LISTL, FILXI, FILETA
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
      INTEGER ISTART, IEND, KEND, LEND

      CALL QENTER('CCEOM_XIETA')

*---------------------------------------------------------------------*
* Main section:   CCEOM_XIETA1, driven by a loop over Xi/Eta vectors
* -------------
*   singles and doubles models:
*               calculation of the single and double excitation
*               parts of the Xi and Eta vectors and respective
*               dot products (IOPTRES=5).
*---------------------------------------------------------------------*

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
        KEND        = IEND     + NTRAN
        LEND        = LWORK    - KEND
        if (locdbg) then
        WRITE (LUPRI,*) 'istart,iend,kend,lend,lwork,ntran',
     &  istart,iend,kend,lend,lwork,ntran
        end if

        IF (LEND .LT. 0) THEN
           WRITE (LUPRI,*) 'Insufficient work space in CCEOM_XIETA.'
           WRITE (LUPRI,*) 'Available    :',LWORK,' words,'
           WRITE (LUPRI,*) 'Need at least:',KEND, ' words.'
           CALL QUIT('Insufficient work space in CCEOM_XIETA.')
        END IF

        IF (LOCDBG) THEN
          WRITE (LUPRI,*) 'Batch No.:',IBATCH
          WRITE (LUPRI,*) 'start at :',IFIRST
          WRITE (LUPRI,*) '# transf.:',NTRAN
          WRITE (LUPRI,*) 'kend     :',KEND
          IDUM = 0
          WRITE (LUPRI,*) 'work(0)  :',WORK(IDUM)
        END IF

        CALL CCEOM_XIETA1(IXETRAN(1,IFIRST),NTRAN,IOPTRES,
     &                  IORDER,LISTL,MXVEC,
     &                  FILXI, IXDOTS(1,IFIRST),XCONS(1,IFIRST),
     &                  FILETA,IEDOTS(1,IFIRST),ECONS(1,IFIRST),
     &                  WORK(ISTART),   WORK(IEND),
     &                  WORK(KEND), LEND )

        IF (LOCDBG) THEN
          WRITE (LUPRI,*) 'returned from CCEOM_XIETA1:'
          IDUM = 0
          WRITE (LUPRI,*) 'work(0)  :',WORK(IDUM)
        END IF

      END DO

      CALL QEXIT('CCEOM_XIETA')

      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCEOM_XIETA                             *
*---------------------------------------------------------------------*
c /* deck ccEOM_xieta1 */
*=====================================================================*
       SUBROUTINE CCEOM_XIETA1( IXETRAN, NXETRAN, IOPTRES,
     &                       IORDER,   LISTL,   MXVEC,
     &                       FILXI,   IXDOTS,  XCONS,
     &                       FILETA,  IEDOTS,  ECONS,
     &                       ISTART,  IEND,
     &                       WORK,    LWORK )
*---------------------------------------------------------------------*
*   NODDY VERSION OF EOM XIETA
*   Purpose: Calculation of the Jacobian
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
#include "eribuf.h"
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
#include "ccexci.h"
#include "ccnoddy.h"
#include "cch2d.h"
#include "ccrc1rsp.h"
#include "r12int.h"
#include "ccr12int.h"

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .false.)
      INTEGER ISYM0
      PARAMETER (ISYM0 = 1) ! symmetry of the reference state
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
      INTEGER LWORK, NXETRAN, IOPTRES, MAXSIM, MXVEC, IORDER
      INTEGER IXETRAN(MXDIM_XEVEC,NXETRAN)
      INTEGER IXDOTS(MXVEC,NXETRAN), IEDOTS(MXVEC,NXETRAN)
      INTEGER ISTART(NXETRAN), IEND(NXETRAN)
!Sonia
      INTEGER IADRF_ETA,IADRF_XI,IATOM, IBATCH, IDLSTL,IDLSTL_OLD
      INTEGER IFILE,IOPER,IOPER_OLD,IOPT,IOPTW,IOPTWE,IOPTWR12
      integer IOPTYP,IRELAX,IRELAX2_OLD,IRELAX2,IRELAX_OLD,ISYCTR
      integer ISYETA,isyhop,isyres, itran,itst
      integer IVEC,KEND0,KEND1,KEND1SV,KEND2,KETA1,KETA2
      logical LLrelax,lltwoel,lrelax,ltwoel,LZEROIMDONE,newtype
      logical skipeta,skipxi
      integer lueta,luxi,LWRK1,LWRK1SV ,LWRK2, kzeta1,kzeta2
      integer mscratch0,mscratch1,mwork, nbatch
      CHARACTER*10 MODELW
      INTEGER :: ISPNETA, ISPNOP, ISPNTR
      INTEGER :: KEND3, LEND3

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL AVERAGE, FREQ, FREQLST, KTEST
      REAL XCONS(MXVEC,NXETRAN), ECONS(MXVEC,NXETRAN)
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION AVERAGE, FREQ, FREQLST, KTEST
      DOUBLE PRECISION XCONS(MXVEC,NXETRAN), ECONS(MXVEC,NXETRAN)
#endif

* external functions:
      INTEGER ILSTSYM
      INTEGER IROPER

      CALL QENTER('CCEOM_XIETA1')

*---------------------------------------------------------------------*
* begin: check wavefunction model and flush print unit
*---------------------------------------------------------------------*

      IF ( .NOT. (CCS .OR. CC2 .OR. CCSD .OR. CC3) ) THEN
        WRITE (LUPRI,'(/1x,a)') 'CCEOM_XIETA Called for a CC '
     &          //'method not implemented in CCEOM_XIETA...'
        CALL QUIT('Unknown CC method in CCEOM_XIETA.')
      END IF

      CALL FLSHFO(LUPRI)

      ! save the present value of the DIRECT flag
      !DIRSAV = DIRECT

      IF (LOCDBG) THEN
        ITST = 0
        WRITE (LUPRI,'(/1x,a,i15)') 'Work space in CCEOM_XIETA:',LWORK
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
         IF (CCSDT) CALL QUIT('Problem in CCEOM_XIETA.')
      ELSE IF (IOPTRES.EQ.2) THEN
         CALL QUIT('IOPTRES=2 option not implemented in CC_XIETA.')
      ELSE IF (IOPTRES.EQ.3 .OR. IOPTRES.EQ.4) THEN
         CONTINUE
      ELSE IF (IOPTRES.EQ.5) THEN
         IF (MXVEC*NXETRAN .LE. 0) THEN
            WRITE (LUPRI,*)
     &      'WARNING: CCEOM_XIETA called, but nothing to do!'
            RETURN
         END IF
      ELSE
         CALL QUIT('Illegal value IOPTRES in CCEOM_XIETA.')
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
         CALL QUIT('Unknown CC model in CCEOM_XIETA.')
      END IF
*---------------------------------------------------------------------*
* prepare batching, make sure that for derivatives all transf. 
* in a batch belong to coordinates of the same atom,
* do not mix derivatives with other perturbations
* -> IOPTYP = 0 : 1el perturbations

      NBATCH  = 1
      MWORK   = 0
      IATOM   = 0
      IOPTYP  = 0
      ISTART(NBATCH) = 1
      DO ITRAN = 1, NXETRAN

        IOPER  = IXETRAN(1,ITRAN)
        LABELH = LBLOPR(IOPER)

        NEWTYPE = (IOPTYP.NE.0)
        IATOM   = 0
        IOPTYP  = 0
!        IF ( (NEWTYPE .AND. ITRAN.GT.1) .OR.
!     &          ((MWORK+MSCRATCH1+MSCRATCH0).GT.LWORK) ) THEN
!              IEND(NBATCH)   = ITRAN-1
!              NBATCH         = NBATCH + 1
!              ISTART(NBATCH) = ITRAN
!              MWORK = 0
!        END IF

!        MWORK = MWORK + MSCRATCH1
!        IF ((MWORK+MSCRATCH0) .GT. LWORK) THEN
!           WRITE (LUPRI,*) 'Insufficient work space in CC_XIETA.'
!           WRITE (LUPRI,*) 'Available:',LWORK,' words.'
!           WRITE (LUPRI,*) 'Need at least:',MWORK+MSCRATCH0,' words.'
!           CALL FLSHFO(LUPRI)
!           CALL QUIT('Insufficient work space in CC_XIETA.')
!        END IF

      END DO
      IEND(NBATCH) = NXETRAN

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'CCEOM_XIETA> batching statistics:'
        WRITE (LUPRI,*) 'CCEOM_XIETA> NBATCH : ',NBATCH
        WRITE (LUPRI,*) 'CCEOM_XIETA> ISTART : ',
     &       (ISTART(IBATCH),IBATCH=1,NBATCH)
        WRITE (LUPRI,*) 'CCEOM_XIETA> IEND   : ',
     &       (IEND(IBATCH),  IBATCH=1,NBATCH)
        WRITE (LUPRI,*) 'WORK(0) = ',WORK(ITST)
      END IF

*=====================================================================*
* loop over the batches of transformations:
*=====================================================================*


      DO IBATCH = 1, NBATCH

*---------------------------------------------------------------------*
* loop over transformations, allocate work space and 
* initialze R, G, F, BF & BFZ intermediates: 
*---------------------------------------------------------------------*
        LLRELAX     = .FALSE.
        LLTWOEL     = .FALSE.
        LZEROIMDONE = .FALSE.

        KEND0 = 1 !Sonia ??? here???
        KEND1 = KEND0
        LWRK1 = LWORK - KEND1

        IOPER_OLD   = -1
        IRELAX_OLD  = -1
        IDLSTL_OLD  = -1

        IRELAX2     = -1
        IRELAX2_OLD = -1

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
        ISPNOP = 1                 ! Only singlet operators for now

        SKIPETA = ( IXETRAN(4,ITRAN) .EQ. -1 )
        LRELAX  = LTWOEL .OR. (IRELAX.GE.1)

        SKIPXI = .true.
C
        IF (.NOT. SKIPETA) THEN
          ISYCTR = ILSTSYM(LISTL,IDLSTL) ! sym. of ZETA vector
          ISYETA = MULD2H(ISYHOP,ISYCTR) ! sym. of ETA result vector
          IF (LISTL(1:2).EQ.'LE') THEN
            ISPNTR = IMULTE(IDLSTL)
          ELSE
            ISPNTR = 1
          ENDIF
          ISPNETA = ISPNTR ! Only singlet operators work for now
        END IF


          KZETA1 = KEND1
          KETA1  = KZETA1 + NT1AM(ISYCTR)
          KEND2  = KETA1  + NT1AM(ISYETA)
          IF (CC2 .OR. CCSD .OR. CCSDT) THEN
             !write(lupri,*)' EOM_XIETA'
             !call flshfo(lupri)
             KETA2 = KEND2
             KEND2 = KETA2  + NT2AM(ISYETA)
             IF (ISPNETA.EQ.3) THEN ! Eta is triplet: twice as long!
                KEND2 = KEND2 + NT2AM(ISYETA)
             END IF
          END IF
          LWRK2 = LWORK - KEND2

          ! Perform self-test: Calculate Xi and compare!!
C          IF (ISPNTR.EQ.3) THEN
          IF (.FALSE.) THEN
             CALL CC_XIC13(ISYMOP,LABELH,WORK(KETA1),FILETA,IDLSTL,
     &                     .TRUE.,DUMMY,WORK(KEND2),LWRK2)
C
             KEND3 = KEND2 + MXVEC
             LEND3 = LWRK2 - MXVEC
             IOPT = 5
C
             CALL CCDOTRSP(IEDOTS,WORK(KEND2),IOPT,LISTL,ITRAN,MXVEC,
     &                     NXETRAN,
     &                     WORK(KETA1),WORK(KETA2),ISYCTR,
     &                     WORK(KEND3),LEND3)

             IVEC = 1
             DO WHILE (IEDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
                WRITE (LUPRI,*)
     &              ' 2 ECONS:',ITRAN,IVEC,WORK(KEND2-1+IVEC),IOPT
                IVEC = IVEC + 1
             END DO

          END IF
          !write(lupri,*)'SONIA: KETA1:', KETA1, ' KETA2:', KETA2
          !write(lupri,*)'SONIA: KEND2:', KEND2
          !call dzero(WORK(KETA1),NT1AM(ISYETA)+NT2AM(ISYETA))

          IF (LWRK2 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CC_XIETA. (CCETA2)')
          END IF

          IF (ISPNTR.EQ.3.AND.ISPNOP.EQ.1) THEN

             CALL CC_ETAC13(ISYHOP,LABELH,WORK(KETA1),LISTL,IDLSTL,
     &                     .TRUE., DUMMY,WORK(KEND2),LWRK2)
          ELSE
             call CCCI_ETAC(ISYHOP,LABELH,WORK(KETA1),LISTL,IDLSTL,
     &                      0, DUMMY,WORK(KEND2),LWRK2)
          END IF

        IF (LOCDBG) THEN
          WRITE (LUPRI,*) 'CCEOM_XIETA> IOPER, IRELAX = ',IOPER,IRELAX
          WRITE (LUPRI,*) 'CCEOM_XIETA> LTWOEL,LRELAX = ',LTWOEL,LRELAX
          WRITE (LUPRI,*) 'CCEOM_XIETA> LABELH,ISYHOP = ',LABELH,ISYHOP
          WRITE (LUPRI,*) 'CCEOM_XIETA> SKIPXI        = ',SKIPXI
          WRITE (LUPRI,*) 'WORK(0) = ',WORK(ITST)
          CALL FLSHFO(LUPRI)
        END IF

        KEND1 = KEND1SV

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

          IF (LWRK2 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CC_XIETA. (CCETA2)')
          END IF

          IOPT = 1
          !useless
          CALL CC_RDRSP(LISTL,IDLSTL,ISYCTR,IOPT,MODEL,
     &                  WORK(KZETA1),DUMMY)

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
            IF (CCSDT) CALL QUIT('Problem in CC_XIETA')
          ELSE IF (IOPTRES.EQ.3) THEN
           IFILE  = IXETRAN(4,ITRAN)
           IF (ILSTSYM(FILETA,IFILE).NE.ISYRES) THEN
             CALL QUIT('Symmetry mismatch for Eta vector in CC_XIETA.')
           END IF
           CALL CC_WRRSP(FILETA,IFILE,ISYRES,IOPTW,MODELW,DUMMY,
     &                   WORK(KETA1),WORK(KETA2),WORK(KEND2),LWRK2)
          ELSE IF (IOPTRES.EQ.4) THEN
           IFILE  = IXETRAN(4,ITRAN)
           IF (ILSTSYM(FILETA,IFILE).NE.ISYRES) THEN
             CALL QUIT('Symmetry mismatch for Eta vector in CC_XIETA.')
           END IF
           CALL CC_WARSP(FILETA,IFILE,ISYRES,IOPTW,MODELW,DUMMY,
     &                   WORK(KETA1),WORK(KETA2),WORK(KEND2),LWRK2)
          ELSE IF (IOPTRES.EQ.5) THEN
           IF (LOCDBG) THEN
             IVEC = 1
             WRITE(LUPRI,*) 'ECONS TRIPLES CONTRIBUTION:'
             DO WHILE (IEDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
                WRITE (LUPRI,*)
     &              '1 ECONS:',IVEC,ITRAN,ECONS(IVEC,ITRAN),IOPTW
                IVEC = IVEC + 1
             END DO
           END IF
C
           !shall we do it?   EOM, SONIA?
           IF ((.NOT.CCS).AND.(ISPNETA.EQ.1))
     &                    CALL CCLR_DIASCL(WORK(KETA2),TWO,ISYRES)
           iopt = ioptw
           if (ispneta.eq.3) iopt = 5
           !
           CALL CCDOTRSP(IEDOTS,ECONS,IOPT,FILETA,ITRAN,NXETRAN,MXVEC,
     &                   WORK(KETA1),WORK(KETA2),ISYRES,
     &                   WORK(KEND2),LWRK2)
           IF (LOCDBG) THEN
             IVEC = 1
             DO WHILE (IEDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
                WRITE (LUPRI,*)
     &              ' 2 ECONS:',IVEC,ITRAN,ECONS(IVEC,ITRAN),IOPT
                IVEC = IVEC + 1
             END DO
           END IF
C
          ELSE
           CALL QUIT('Illegal value for IOPTRES in CCEOM_XIETA.')
          END IF
C
          IF (LOCDBG) THEN
            WRITE (LUPRI,*) 'Final result of CCEOM_XIETA:'
            WRITE (LUPRI,*) 'operator label:      ',LBLOPR(IOPER)
            WRITE (LUPRI,*) 'output file type:    ',FILETA
            WRITE (LUPRI,*) 'index of output file:',IFILE
            WRITE (LUPRI,*) 'two-electron contr.: ',LTWOEL
            WRITE (LUPRI,*) 'relax/reorth. contr.:',LRELAX
            I = 1
            IF (CCS) I = 0
            CALL CC_PRP(WORK(KETA1),WORK(KETA2),ISYRES,1,I)
            CALL FLSHFO(LUPRI)
          END IF

        END IF
*=====================================================================*
* end of loop over transformations within this batch:
*=====================================================================*

      END DO ! ITRAN

*---------------------------------------------------------------------*
* close the loop over batches
*---------------------------------------------------------------------*
      END DO ! IBATCH 

      IF (IOPTRES.EQ.0) THEN
         CALL WCLOSE2(LUXI, FILXI, 'KEEP')
         CALL WCLOSE2(LUETA,FILETA,'KEEP')
      ELSE IF (IOPTRES.EQ.1) THEN
         CALL WCLOSE2(LUXI, FILXI, 'DELETE')
         CALL WCLOSE2(LUETA,FILETA,'DELETE')
         CALL QUIT('IOPTRES=1 not yet implemented in CCEOM_XIETA.')
      END IF

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'CCEOM_XIETA ended successfully (?) ... '//
     &        'return now.'
        CALL FLSHFO(LUPRI)
      END IF

      ! restore the DIRECT flag 
      !DIRECT = DIRSAV 

      CALL FLSHFO(LUPRI)
      CALL QEXIT('CCEOM_XIETA1')
      RETURN
      END
*=====================================================================*
*                  END OF SUBROUTINE CCEOM_XIETA1                     *
*=====================================================================*

c*DECK CCCI_ETAC
      SUBROUTINE CCCI_ETAC(ISYMC,LBLC,        !Operator stuff
     &                   ETAC,                !Result vector
     &                   LIST,ILSTNR,         !Left vector
     &                   IOPTCC2,
     &                   XINT,WORK,LWORK)
C
C-----------------------------------------------------------------------
C     Purpose: Calculate the CCCI/EOMCC etaC(L0) or
C     etaC(L) vector (Modified version of the CCSD etaC(l0/l1) code
C
C     Important note: Requires work space of dimension of
C             NT2AM + NT2SQ in addition to ETAC, so please take care.
C
C     eta(tau_nu)= (<HF| + Sum(mu)L(0 or 1)<mu|)
C                         exp(-t)[C,tau_nu]exp(T)|HF>
C
C     LIST= 'L0' for zeroth order left amplitudes.
C                ISYML should be ISYMOP in this case.
C
C           'L1' for first order left amplitudes, read in from file
C                In this case the vector is found according to its list
C                number ILSTNR.
C
C                For L1 HF contribution is skipped.
C
C     IOPTCC2 = 1 -- transform for CC2 the Fock matrix entering the
C                    E term contribution with CMO vector instead with
C                    Lambda matrices
C
C     IOPTL0 = 1 -- The L0 vectors are passed in memory
C
C     C property integrals read according to LBLC
C
C     SLV98,OC: Allow for input of integrals if
C               LBLC.eq.'GIVE INT'
C
C     Sonia & Filip, Maj 2015
C
C-----------------------------------------------------------------------
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxorb.h"
#include "ccorb.h"
#include "iratdef.h"
#include "cclr.h"
#include "ccexci.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
C
!Sonia & Filip: find a better place
#include "ccsections.h"
!
      PARAMETER( TWO = 2.0D00, XHALF = 0.5D00 )
      DIMENSION ETAC(*),WORK(LWORK),XINT(*)
      CHARACTER LBLC*(*),LIST*(*),MODEL*10
      INTEGER IOPTCC2
      LOGICAL  FCKCON,ETRAN, LOCDBG, LEOM
      PARAMETER( LOCDBG = .false., LEOM = .TRUE. )
C
!      CALL AROUND( 'IN CCCI_ETAC: Constructing ETAC(L1) vector ')
      IF ( (IPRINT .GT. 10).or.locdbg ) THEN
         CALL AROUND( 'IN CCCI_ETAC: Constructing ETAC(L1) vector ')
      ENDIF
C
C--------------------------------
C     find symmetry of D operator.
C--------------------------------
C
      ISYML = ILSTSYM(LIST,ILSTNR)
C
      ISYRES = MULD2H(ISYML,ISYMC)
      IF (( LIST .EQ. 'L0').AND.(ISYML.NE.1)) THEN
         CALL QUIT('Misuse of CCCI_ETAC')
      ENDIF
C
      TIMEC = SECOND()
C
      MODEL = 'CCSD      '
      IF (CCS) MODEL = 'CCS       '
      IF (CC2) MODEL = 'CC2       '
C
C--------------------
C     Allocate space.
C--------------------
C
      KCTMO  = 1
      KT1AM  = KCTMO  + N2BST(ISYMC)
      KLAMDP = KT1AM  + NT1AM(ISYMOP)
      KLAMDH = KLAMDP + NLAMDT
      KEND1  = KLAMDH + NLAMDT
C
      IF (CC2 .AND. IOPTCC2.EQ.1) THEN
        KCMO   = KEND1
        KFCKHF = KCMO   + NLAMDT
        KEND1  = KFCKHF + N2BST(ISYMC)
      END IF
C
      LEND1  = LWORK  - KEND1
C
      IF ( .NOT. CCS) THEN

         !if (EOMCCSD.or.CCCISD) then
           !Integrals in MO basis
           !highly redundant, but I need a second check
           !write(lupri,*)'CCCI_ETAC: EOM alloc for X_mo integrals'
         KINTIA = KEND1
         KEND1  = KINTIA + NT1AM(ISYMC)
         LEND1  = LWORK - KEND1
         CALL DZERO(WORK(KintIA),NT1AM(isymc))
         !end if
C
         KL1AM = KEND1
         KL2AM = KL1AM + NT1AM(ISYML)
         KEND2 = KL2AM + NT2SQ(ISYML)
         LEND2 = LWORK - KEND2
         KT2AM = KEND2
         KEND21= KT2AM + MAX(NT2AM(ISYML),NT2AM(1))
         LEND21= LWORK - KEND21
C
      ELSE
C
         KL1AM = KEND1
         KEND2 = KL1AM + NT1AM(ISYML)
         LEND2 = LEND1
         KEND21= KEND1
         LEND21= LEND1
C
      ENDIF
C
      IF (LEND21.LT. 0 ) CALL QUIT('TOO LITTLE WORKSPACE IN CCCI_ETAC')
C
C-----------------------
C     get T1 amplitudes.
C-----------------------
C
      CALL DZERO(WORK(KT1AM),NT1AM(1))
      IF ( .NOT. CCS) THEN
         IOPT = 1
         CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KT1AM),DUMMY)
         !skod = DDOT(NT1AM(1),WORK(KT1AM),1,WORK(KT1AM),1)
         !WRITE(LUPRI,1) 'Norm of T0_1:       ',skod
      ENDIF
C
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),
     *            WORK(KEND21),LEND21)
C
      IF (CC2 .AND. IOPTCC2.EQ.1) THEN
        LUSIFC = -1
        CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',
     *              IDUMMY,.FALSE.)
        REWIND(LUSIFC)
        CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
        READ(LUSIFC)
        READ(LUSIFC)
        READ(LUSIFC) (WORK(KCMO+I-1),I=1,NLAMDS)
        CALL GPCLOSE(LUSIFC,'KEEP')
        CALL CMO_REORDER(WORK(KCMO),WORK(KEND21),LEND21)
      END IF
C
C-------------------------------
C     get AO property integrals.
C-------------------------------
C
      CALL DZERO(WORK(KCTMO),N2BST(ISYMC))
      FF = 1.0D0
C SLV98,OC give integrals option
      IF (LBLC.EQ.'GIVE INT') THEN
        CALL DCOPY(N2BST(ISYMC),XINT(1),1,WORK(KCTMO),1)
      ELSE
        FF = 1.0D0
        CALL CC_ONEP(WORK(KCTMO),WORK(KEND21),LEND21,FF,ISYMC,LBLC)
      ENDIF
C
      IF (CC2 .AND. IOPTCC2.EQ.1) THEN
        CALL DCOPY(N2BST(ISYMC),WORK(KCTMO),1,WORK(KFCKHF),1)
        CALL CC_FCKMO(WORK(KFCKHF),WORK(KCMO),WORK(KCMO),
     *                WORK(KEND21),LEND21,ISYMC,1,1)
      END IF
C
C-----------------------------------------------
C     Make MO T1-transformed property integrals.
C-----------------------------------------------
C
      !X_pq
      CALL CC_FCKMO(WORK(KCTMO),WORK(KLAMDP),WORK(KLAMDH),
     *              WORK(KEND21),LEND21,ISYMC,1,1)
C
C----------------------------------------------------------
C     Calculate 2Cia (stored ia) Hartree-Fock contribution.
C----------------------------------------------------------
C
      CALL DZERO(ETAC,NT1AM(ISYRES))

      DO 100 ISYMI = 1,NSYM
C
         ISYMA = MULD2H(ISYMI,ISYMC)
C
         DO 110 A = 1,NVIR(ISYMA)
C
            DO 120 I = 1,NRHF(ISYMI)
C
               KOFF1 = KINTIA + IT1AM(ISYMA,ISYMI)
     &               + NVIR(ISYMA)*(I - 1) + A-1
               KOFF2 = KCTMO + IFCVIR(ISYMI,ISYMA)
     *               + NORB(ISYMI)*(A - 1) + I - 1
C
               WORK(KOFF1) = WORK(KOFF2)
C
  120       CONTINUE
  110    CONTINUE
C
  100 CONTINUE
C
      IF (LIST .EQ. 'L0') THEN
         CALL DAXPY(NT1AM(ISYRES),TWO,WORK(KINTIA),1,ETAC(1),1)
      ENDIF
C----------------------------------------
      IF ( DEBUG.or.locdbg ) THEN
         ETA1 = DDOT(NT1AM(ISYRES),ETAC(1),1,ETAC(1),1)
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,1) 'Norm of ETAC - 2Cia contribution:',ETA1
      ENDIF

      !return
C
C------------------------
C     IF CCS then return.
C------------------------
C
      IF ( CCS .AND. (LIST .EQ. 'L0')) RETURN
C
C----------------------------------------------
C     Read zero'th order multipliers.
C     Tbar2 are SQUARED
C----------------------------------------------
C
      IOPT = 3
      CALL CC_RDRSP(LIST,ILSTNR,ISYML,IOPT,MODEL,
     *              WORK(KL1AM),WORK(KT2AM))
C      skod = DDOT(NT1AM(ISYML),WORK(KL1AM),1,WORK(KL1AM),1)
C      skod = skod+DDOT(NT2AM(ISYML),WORK(KT2AM),1,WORK(KT2AM),1)
      !WRITE(LUPRI,1) 'Norm of Tbar (full):       ',skod
      IF (.NOT. CCS) CALL CC_T2SQ(WORK(KT2AM),WORK(KL2AM),ISYML)
C
C--------------------------------
C     Put T2 (packed) amplitudes in etac2.
C--------------------------------
C
      IF (.NOT. CCS) THEN
         !read in T2AM (packed)
         IOPT = 2
         CALL CC_RDRSP('R0',0,1,IOPT,MODEL,DUMMY,WORK(KT2AM))
C         skod = DDOT(NT2AM(1),WORK(KT2AM),1,WORK(KT2AM),1)
         ! WRITE(LUPRI,1) 'Norm of T2(0) (full):   ',skod
      ENDIF
C
C--------------------------------
C     Make X and Y intermediates.
C--------------------------------
C
      IF (.NOT. CCS) THEN
         KXMAT = KEND21
         KYMAT = KXMAT + NMATIJ(ISYML)
         KEND3 = KYMAT + NMATAB(ISYML)
         LEND3 = LWORK - KEND3
         IF (LEND3.LT. 0 )
     &        CALL QUIT(' TOO LITTLE WORKSPACE IN CC_ETAC-2')
C
         IF ( DEBUG.or.LOCDBG ) THEN
            XYI   = DDOT(NT1AM(ISYML),WORK(KL1AM),1,WORK(KL1AM),1)
            WRITE(LUPRI,1) 'CC_ETAC: L1AM vector:              ',XYI
            XYI   = DDOT(NT2SQ(ISYML),WORK(KL2AM),1,WORK(KL2AM),1)
            WRITE(LUPRI,1) 'CC_ETAC: L2AM vector:              ',XYI
            XXI   = DDOT(NT2AM(ISYMOP),WORK(KT2AM),1,WORK(KT2AM),1)
            WRITE(LUPRI,1) 'T2AM vector :                      ',XXI
         ENDIF
         CALL CC_XI(WORK(KXMAT),WORK(KL2AM),ISYML,WORK(KT2AM),1,
     *              WORK(KEND3),LEND3)
         CALL CC_YI(WORK(KYMAT),WORK(KL2AM),ISYML,WORK(KT2AM),1,
     *              WORK(KEND3),LEND3)
         IF ( DEBUG.or.LOCDBG ) THEN
            XYI   = DDOT(NMATAB(ISYML),WORK(KYMAT),1,WORK(KYMAT),1)
            WRITE(LUPRI,1) 'CC_ETAC: YI  intermediate is:      ',XYI
            XXI   = DDOT(NMATIJ(ISYML),WORK(KXMAT),1,WORK(KXMAT),1)
            WRITE(LUPRI,1) 'CC_ETAC: XI  intermediate is:      ',XXI
         ENDIF
      ELSE
         KEND3 = KEND2
         LEND3 = LEND2
      ENDIF
C
C----------------------------------------------
C     Calculate X and Y contributions to etac1.
C     etac1 = -sum(e)Cie*Yae - sum(l)Cla*Xli
C----------------------------------------------
C
      IF ( (.NOT.CCS) .AND. (.NOT.(CC2.AND.IOPTCC2.EQ.1)) ) THEN
C      DN = DDOT(N2BST(ISYMC),WORK(KCTMO),1,WORK(KCTMO),1)
      !WRITE(LUPRI,*) 'Integral norm before CC_21EFM out norm:',DN
         CALL CC_21EFM(ETAC,WORK(KCTMO),ISYMC,WORK(KXMAT),
     *                 WORK(KYMAT),ISYML,WORK(KEND3),LEND3)
C      DN = DDOT(N2BST(ISYMC),WORK(KCTMO),1,WORK(KCTMO),1)
      !WRITE(LUPRI,*) 'Integral norm after CC_21EFM out norm:',DN
C
         IF ( DEBUG.or.locdbg ) THEN
            ETA1 = DDOT(NT1AM(ISYRES),ETAC(1),1,ETAC(1),1)
            WRITE(LUPRI,1) 'Norm of eta1-after X&Y cont:       ',ETA1
         ENDIF
      ENDIF

      IF (LEOM) THEN
C
C---------------------------------------------------------------
C     EOM contribution to ETAC:
C                     ~                  ~
C     etac = 2sum(ei) L_{ckei}*(C_{ei} + t_{ei,fn}*C_{nf})
C---------------------------------------------------------------
         KCINT = KEND21
         KETAC = KCINT
         KINTAI = KCINT + MAX(NT1AM(ISYMC),NT1AM(ISYRES))
C        Extract C_{ai}
         DO ISYMI = 1, NSYM
            ISYMA = MULD2H(ISYMI,ISYMC)
            DO I = 1, NRHF(ISYMI)
               KOFF1 = KINTAI + IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1)
               KOFF2 = KCTMO + IFCRHF(ISYMA,ISYMI) +
     &                 NORB(ISYMA)*(I-1) + NRHF(ISYMA)
               CALL DCOPY(NVIR(ISYMA),WORK(KOFF2),1,WORK(KOFF1),1)
            END DO
         END DO
         !get ttilde (overwrites T2am)
         CALL CCSD_TCMEPK(WORK(KT2AM),1.0D0,1,1)
C            ~
C        Add t_{em,fn} *C_{nf} to C_{em}
         CALL CCG_LXD(WORK(KCINT),ISYMC,WORK(KINTIA),ISYMC,
     &                WORK(KT2AM),1,0)
         CALL DAXPY(NT1AM(ISYMC),1.D0,WORK(KCINT),1,WORK(KINTAI),1)
C        Calculate the term as L_{ai,em} * C'_{em}
         CALL CCG_LXD(WORK(KETAC),ISYRES,WORK(KINTAI),ISYMC,
     &                                   WORK(KL2AM),ISYML,1)
         IF ( DEBUG.or.locdbg ) THEN
            ETA1 = DDOT(NT1AM(ISYRES),WORK(KETAC),1,WORK(KETAC),1)
            WRITE(LUPRI,1) 'Norm of alone Tbar_ck,ei*Cei: ',ETA1
         ENDIF
         !removed factor 2 to get FCI limit!
         CALL DAXPY(NT1AM(ISYRES),1.0D0,WORK(KETAC),1,ETAC(1),1)
      END IF
C
C------------------------------------------------
C     Workspace for T2AM and X and Y is now free.
C     etac2 = P(ab,ij)(2l(ai)*Cjb - l(aj)*c(ib))
C------------------------------------------------
C
      IF (.NOT. CCS) THEN
         CALL DZERO(ETAC(1+NT1AM(ISYRES)),NT2AM(ISYRES))

C      DN = DDOT(N2BST(ISYMC),WORK(KCTMO),1,WORK(KCTMO),1)
      !WRITE(LUPRI,*) 'Integral norm before CC_L1FCK out norm:',DN

         CALL CC_L1FCK(ETAC(1+NT1AM(ISYRES)),WORK(KL1AM),WORK(KCTMO),
     *                 ISYML,ISYMC,WORK(KEND2),LEND2)

C      DN = DDOT(N2BST(ISYMC),WORK(KCTMO),1,WORK(KCTMO),1)
!      WRITE(LUPRI,*) 'Integral norm after CC_L1FCK out norm:',DN
!     *                 ISYML,ISYMC,WORK(KEND2),LEND2)
C
         IF ( DEBUG.or.locdbg ) THEN
            ETA1 = DDOT(NT1AM(ISYRES),ETAC(1),1,ETAC(1),1)
            ETA2 = DDOT(NT2AM(ISYRES),ETAC(1+NT1AM(ISYRES)),1,
     *                  ETAC(1+NT1AM(ISYRES)),1)
            WRITE(LUPRI,1) 'Norm of eta1-after L1c cont:       ',ETA1
            WRITE(LUPRI,1) 'Norm of eta2-after L1c cont:       ',ETA2
         ENDIF
      ENDIF
C
      KEI1   = KEND2
      KEI2   = KEI1   + NEMAT1(ISYMC)
      KEND3  = KEI2   + NMATIJ(ISYMC)
      LEND3  = LWORK  - KEND3
      IF (LEND3.LT. 0 ) CALL QUIT(' TOO LITTLE WORKSPACE IN CC_ETAC-3')
C
C--------------------------------
C     Put A into E matrix format.
C--------------------------------
C
      FCKCON = .TRUE.
      ETRAN  = .FALSE.
      CALL CCRHS_EFCK(WORK(KEI1),WORK(KEI2),WORK(KLAMDH),
     *                WORK(KCTMO),WORK(KEND3),LEND3,FCKCON,
     *                ETRAN,ISYMC)
C
C--------------------------------------------
C     etac1 =  sum(b)Lbi*Cba - sum(j)Laj*Cij.
C--------------------------------------------
C
      IF ( DEBUG.or.locdbg ) THEN
         XE1 = DDOT(NMATAB(ISYMC),WORK(KEI1),1,WORK(KEI1),1)
         XE2 = DDOT(NMATIJ(ISYMC),WORK(KEI2),1,WORK(KEI2),1)
         WRITE(LUPRI,1) 'Norm of EI1  -after EFCK:          ',XE1
         WRITE(LUPRI,1) 'Norm of EI2  -after EFCK:          ',XE2
         XE2 = DDOT(NMATIJ(ISYMC),WORK(KINTIJ),1,WORK(KINTIJ),1)
         WRITE(LUPRI,1) 'Norm of Cij integrals  -after EFCK: ',XE2
         ETA1 = DDOT(NT1AM(ISYML),WORK(KL1AM),1,WORK(KL1AM),1)
         WRITE(LUPRI,1) 'Norm of L1AM before  CCLR_E1C1:    ',ETA1
      ENDIF
C
      CALL CCLR_E1C1(ETAC,WORK(KL1AM),WORK(KEI1),WORK(KEI2),
     *               WORK(KEND3),LEND3,ISYML,ISYMC,'T')
C
      IF ( DEBUG.or.locdbg ) THEN
         ETA1 = DDOT(NT1AM(ISYRES),ETAC(1),1,ETAC(1),1)
         WRITE(LUPRI,1) 'Norm of eta1 - after CCLR_E1C1:    ',ETA1
      ENDIF
C
C---------------------------------------------------------------
C     etac2 = P(ab,ij)(sum(e)2L(aiej)*Ceb - sym(k)L(aibk)*C(jk))
C---------------------------------------------------------------
C
      IF (.NOT. CCS) THEN
C
         IF (CC2 .AND. IOPTCC2.EQ.1) THEN
           FCKCON = .TRUE.
           ETRAN  = .FALSE.
           CALL CCRHS_EFCK(WORK(KEI1),WORK(KEI2),WORK(KCMO),
     *             WORK(KFCKHF),WORK(KEND3),LEND3,FCKCON,ETRAN,ISYMC)
         END IF

         CALL CC_EITR(WORK(KEI1),WORK(KEI2),WORK(KEND3),LEND3,
     *                ISYMC)
C
         CALL CCRHS_E(ETAC(1+NT1AM(ISYRES)),WORK(KL2AM),
     *                WORK(KEI1),WORK(KEI2),WORK(KEND3),
     *                LEND3,ISYML,ISYMC)
C
         IF (IPRINT .GT. 40 ) THEN
            CALL AROUND( 'In CCCI_ETAC:  EtaC vector ' )
            CALL CC_PRP(ETAC(1),ETAC(1+NT1AM(ISYRES)),ISYMC,1,1)
         ENDIF
C
         IF (DEBUG .OR. ( IPRINT .GT. 20 )) THEN
            ETA1 = DDOT(NT1AM(ISYRES),ETAC(1),1,ETAC(1),1)
            ETA2 = DDOT(NT2AM(ISYRES),ETAC(1+NT1AM(ISYRES)),1,
     *                  ETAC(1+NT1AM(ISYRES)),1)
            WRITE(LUPRI,1) 'Norm of eta1 - end of CCCI_ETAC:   ',ETA1
            WRITE(LUPRI,1) 'Norm of eta2 - end of CCCI_ETAC:   ',ETA2
            CALL AROUND( 'END OF CC_ETAC ')
         ENDIF
      ENDIF

      !end if
      IF (IPRINT .GT. 5 ) THEN
         TIMEC = SECOND() - TIMEC
         WRITE(LUPRI,9999) 'CCCI_ETA^C      ', TIMEC
      ENDIF
C
   1  FORMAT(1x,A35,1X,E20.10)
9999  FORMAT(1x,'Time used in',2x,A18,2x,': ',f10.2,' seconds')
C
      END
c*DECK CCCI_ETAC2
      SUBROUTINE CCCI_ETAC2(ISYMC,LBLC,        !Operator stuff
     &                   ETAC,                !Result vector
     &                   LIST,ILSTNR,         !Left vector
     &                   IOPTCC2,
     &                   ioptl0,xl0,
     &                   XINT,WORK,LWORK)
C
C-----------------------------------------------------------------------
C     Purpose: Calculate the CCCI/EOMCC etaC(L0) or
C     etaC(L) vector (Modified version of the CCSD etaC(l0/l1) code
C
C     Important note: Requires work space of dimension of
C             NT2AM + NT2SQ in addition to ETAC, so please take care.
C
C     eta(tau_nu)= (<HF| + Sum(mu)L(0 or 1)<mu|)
C                         exp(-t)[C,tau_nu]exp(T)|HF>
C
C     LIST= 'L0' for zeroth order left amplitudes.
C                ISYML should be ISYMOP in this case.
C
C           'L1' for first order left amplitudes, read in from file
C                In this case the vector is found according to its list
C                number ILSTNR.
C
C                For L1 HF contribution is skipped.
C
C     IOPTCC2 = 1 -- transform for CC2 the Fock matrix entering the
C                    E term contribution with CMO vector instead with
C                    Lambda matrices
C
C     IOPTL0 = 1 -- The L0 vectors are passed in memory and Hartree-Fock part skipped.
C     IOPTL0 = 2 -- only Hartree-Fock part calculated.
C
C     C property integrals read according to LBLC
C
C     SLV98,OC: Allow for input of integrals if
C               LBLC.eq.'GIVE INT'
C
C     Sonia & Filip, Maj 2015
C
C-----------------------------------------------------------------------
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxorb.h"
#include "ccorb.h"
#include "iratdef.h"
#include "cclr.h"
#include "ccexci.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
C
!Sonia & Filip: find a better place
#include "ccsections.h"
#include "fdeom.h"
!
      PARAMETER( TWO = 2.0D00, XHALF = 0.5D00 )
      DIMENSION ETAC(*),WORK(LWORK),XINT(*), Xl0(*)
      CHARACTER LBLC*(*),LIST*(*),MODEL*10
      INTEGER IOPTCC2, ioptl0
      LOGICAL  FCKCON,ETRAN, LOCDBG
      PARAMETER( LOCDBG = .true. )
C
      IF ( (IPRINT .GT. 10).or.locdbg ) THEN
         CALL AROUND( 'IN CCCI_ETAC2: Constructing ETAC(L1) vector ')
      ENDIF
C
      if ( (ioptl0.ne.1) .and. (ioptl0.ne.2) ) then
         call quit('Unknown value of ioptl0 in CCCI_ETAC2')
      end if
C
      if ( (ioptl0.eq.2) .and. (LIST(1:2).ne.'L0') ) then
         call quit('Wrong setup in CCCI_ETAC2')
      end if
C
C--------------------------------
C     find symmetry of D operator.
C--------------------------------
C
      ISYML = ILSTSYM(LIST,ILSTNR)
C
      ISYRES = MULD2H(ISYML,ISYMC)
      IF (( LIST .EQ. 'L0').AND.(ISYML.NE.1)) THEN
         CALL QUIT('Misuse of CCCI_ETAC')
      ENDIF
C
      TIMEC = SECOND()
C
      MODEL = 'CCSD      '
      IF (CCS) MODEL = 'CCS       '
      IF (CC2) MODEL = 'CC2       '
C
C--------------------
C     Allocate space.
C--------------------
C
      KCTMO  = 1
      KT1AM  = KCTMO  + N2BST(ISYMC)
      KLAMDP = KT1AM  + NT1AM(ISYMOP)
      KLAMDH = KLAMDP + NLAMDT
      KEND1  = KLAMDH + NLAMDT
C
      IF (CC2 .AND. IOPTCC2.EQ.1) THEN
        KCMO   = KEND1
        KFCKHF = KCMO   + NLAMDT
        KEND1  = KFCKHF + N2BST(ISYMC)
      END IF
C
      LEND1  = LWORK  - KEND1
C
      IF ( .NOT. CCS) THEN

         !if (EOMCCSD.or.CCCISD.or.FDEOM) then
           !Integrals in MO basis
           !highly redundant, but I need a second check
           write(lupri,*)'CCCI_ETAC2: EOM alloc for X_mo integrals'
           KINTAI = KEND1
           KINTAB = KINTAI + NT1AM(ISYMC)
           KINTIJ = KINTAB + NMATAB(ISYMC)
           KINTIA = KINTIJ + NMATIJ(ISYMC)
           KEND1  = KINTIA + NT1AM(ISYMC)
           LEND1  = LWORK - KEND1
           CALL DZERO(WORK(KINTAB),NMATAB(isymc))
           CALL DZERO(WORK(KintIJ),NMATIJ(isymc))
           CALL DZERO(WORK(KintIA),NT1AM(isymc))
           CALL DZERO(WORK(KintAI),NT1AM(isymc))
         !end if
C
         KL1AM = KEND1
         KL2AM = KL1AM + NT1AM(ISYML)
         KEND2 = KL2AM + NT2SQ(ISYML)
         LEND2 = LWORK - KEND2
         KT2AM = KEND2
         KEND21= KT2AM + MAX(NT2AM(ISYML),NT2AM(1))
         LEND21= LWORK - KEND21
C
      ELSE
C
         KL1AM = KEND1
         KEND2 = KL1AM + NT1AM(ISYML)
         LEND2 = LEND1
         KEND21= KEND1
         LEND21= LEND1
C
      ENDIF
C
      IF (LEND21.LT. 0 ) CALL QUIT('TOO LITTLE WORKSPACE IN CCCI_ETAC2')
C
C-----------------------
C     get T1 amplitudes.
C-----------------------
C
      CALL DZERO(WORK(KT1AM),NT1AM(1))
      IF ( .NOT. CCS) THEN
         IOPT = 1
         CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KT1AM),DUMMY)
         skod = DDOT(NT1AM(1),WORK(KT1AM),1,WORK(KT1AM),1)
         WRITE(LUPRI,1) 'Norm of T0_1:       ',skod
      ENDIF
C
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),
     *            WORK(KEND21),LEND21)
C
      IF (CC2 .AND. IOPTCC2.EQ.1) THEN
        LUSIFC = -1
        CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',
     *              IDUMMY,.FALSE.)
        REWIND(LUSIFC)
        CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
        READ(LUSIFC)
        READ(LUSIFC)
        READ(LUSIFC) (WORK(KCMO+I-1),I=1,NLAMDS)
        CALL GPCLOSE(LUSIFC,'KEEP')
        CALL CMO_REORDER(WORK(KCMO),WORK(KEND21),LEND21)
      END IF
C
C-------------------------------
C     get AO property integrals.
C-------------------------------
C
      CALL DZERO(WORK(KCTMO),N2BST(ISYMC))
      FF = 1.0D0
C SLV98,OC give integrals option
      IF (LBLC.EQ.'GIVE INT') THEN
        CALL DCOPY(N2BST(ISYMC),XINT(1),1,WORK(KCTMO),1)
      ELSE
        FF = 1.0D0
        CALL CC_ONEP(WORK(KCTMO),WORK(KEND21),LEND21,FF,ISYMC,LBLC)
      ENDIF
C
      IF (CC2 .AND. IOPTCC2.EQ.1) THEN
        CALL DCOPY(N2BST(ISYMC),WORK(KCTMO),1,WORK(KFCKHF),1)
        CALL CC_FCKMO(WORK(KFCKHF),WORK(KCMO),WORK(KCMO),
     *                WORK(KEND21),LEND21,ISYMC,1,1)
      END IF
C
C-----------------------------------------------
C     Make MO T1-transformed property integrals.
C-----------------------------------------------
C
      !if (eomccsd.or.cccisd.or.fdeom) then

         DN = DDOT(N2BST(ISYMC),WORK(KCTMO),1,WORK(KCTMO),1)
         WRITE(LUPRI,*) 'IN ONEP: FOCK out norm:',DN

         write(lupri,*) 'X_mo integrals of EOM from CCDINTMO'
         ISYM = ISYMC
         CALL CCDINTMO(WORK(KINTIJ),WORK(KINTIA),WORK(KINTAB),
     *                 WORK(KINTAI),WORK(KCTMO),WORK(KLAMDP),
     *                 WORK(KLAMDH),WORK(KEND21),LEND21,ISYM)
         IF ( DEBUG.or.locdbg ) THEN
            ETA1 = DDOT(NT1AM(ISYMC),WORK(KINTIA),1,WORK(KINTIA),1)
            WRITE(LUPRI,*) ' '
            WRITE(LUPRI,1) '4*Norm of C_ie for EOM:',ETA1*4.0d0
         ENDIF
         IF ((IPRINT .GT. 50).or.locdbg) THEN
C
            CALL AROUND('One electron integrals in MO-basis')
            DO 10 ISYM = 1,NSYM
               WRITE(LUPRI,333) 'Symmetry block number:', ISYM
  333          FORMAT(/3X,A22,2X,I1)
               KOFFIJ = KINTIJ + IMATIJ(ISYM,ISYM)
               KOFFAI = KINTAI + IT1AM(ISYM,ISYM)
               KOFFIA = KINTIA + IT1AM(ISYM,ISYM)
               KOFFAB = KINTAB + IMATAB(ISYM,ISYM)
!               CALL AROUND('occ.occ part')
!               CALL OUTPUT(WORK(KOFFIJ),1,NRHF(ISYM),1,NRHF(ISYM),
!     *                     NRHF(ISYM),NRHF(ISYM),1,LUPRI)
!               CALL AROUND('vir.occ part')
!               CALL OUTPUT(WORK(KOFFAI),1,NVIR(ISYM),1,NRHF(ISYM),
!     *                     NVIR(ISYM),NRHF(ISYM),1,LUPRI)
!               CALL AROUND('occ.vir part')
!               CALL OUTPUT(WORK(KOFFIA),1,NVIR(ISYM),1,NRHF(ISYM),
!     *                     NVIR(ISYM),NRHF(ISYM),1,LUPRI)
!               CALL AROUND('vir.vir part')
!               CALL OUTPUT(WORK(KOFFAB),1,NVIR(ISYM),1,NVIR(ISYM),
!     *                     NVIR(ISYM),NVIR(ISYM),1,LUPRI)
C
  10        CONTINUE
            HIJNO = DDOT(NMATIJ(isymc),WORK(KINTIJ),1,WORK(KINTIJ),1)
            HAINO = DDOT(NT1AM(isymc),WORK(KINTAI),1,WORK(KINTAI),1)
            HIANO = DDOT(NT1AM(isymc),WORK(KINTIA),1,WORK(KINTIA),1)
            HABNO = DDOT(NMATAB(isymc),WORK(KINTAB),1,WORK(KINTAB),1)
            WRITE(LUPRI,*) ' '
            WRITE(LUPRI,*) 'Norm of one electron integrals in MO-basis:'
            WRITE(LUPRI,*) HIJNO, HAINO, HIANO, HABNO
            WRITE(LUPRI,*) 'Integrals sum: ', HIJNO+HAINO+HIANO+HABNO
         ENDIF
      !end if

       DN = DDOT(N2BST(ISYMC),WORK(KCTMO),1,WORK(KCTMO),1)
       WRITE(LUPRI,*) 'IN ONEP: FOCK out norm:',DN

      !X_pq
      CALL CC_FCKMO(WORK(KCTMO),WORK(KLAMDP),WORK(KLAMDH),
     *              WORK(KEND21),LEND21,ISYMC,1,1)
       DN = DDOT(N2BST(ISYMC),WORK(KCTMO),1,WORK(KCTMO),1)
       WRITE(LUPRI,*) 'Integral norm after FCKMO out norm:',DN
C
C----------------------------------------------------------
C     Calculate 2Cia (stored ia) Hartree-Fock contribution.
C----------------------------------------------------------
C
      CALL DZERO(ETAC,NT1AM(ISYRES))
C
      IF ((LIST .EQ. 'L0').and.(ioptl0.ne.1)) THEN
         DO 100 ISYMI = 1,NSYM
C
            ISYMA = MULD2H(ISYMI,ISYMC)
C
            DO 110 A = 1,NVIR(ISYMA)
C
               DO 120 I = 1,NRHF(ISYMI)
C
                  KOFF1 = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + A
                  KOFF2 = KCTMO + IFCVIR(ISYMI,ISYMA)
     *                  + NORB(ISYMI)*(A - 1) + I - 1
C
                  ETAC(KOFF1) = TWO*WORK(KOFF2)
C
  120          CONTINUE
  110       CONTINUE
  100    CONTINUE
C
      ENDIF
!
!     For ioptl0 = 2 we want only the Hartree-Fock contribution (which has just been
!     (calculated above), so skip the rest:
      if (ioptl0.eq.2) then
         go to 8888
      end if
C----------------------------------------
      IF ( DEBUG.or.locdbg ) THEN
         ETA1 = DDOT(NT1AM(ISYRES),ETAC(1),1,ETAC(1),1)
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,1) 'Norm of ETAC - 2Cia contribution:',ETA1
      ENDIF

      !return
C
C------------------------
C     IF CCS then return.
C------------------------
C
      IF ( CCS .AND. (LIST .EQ. 'L0')) RETURN
C
C----------------------------------------------
C     Read zero'th order multipliers.
C     Tbar2 are SQUARED
C----------------------------------------------
C
      if (ioptl0.eq.1) then
        call dcopy(NT1AM(ISYML),XL0(1),1,WORK(KL1AM),1)
        call dcopy(NT2AM(ISYML),XL0(1+NT1AM(ISYML)),1,WORK(KT2AM),1)
      else
        IOPT = 3
        CALL CC_RDRSP(LIST,ILSTNR,ISYML,IOPT,MODEL,
     *              WORK(KL1AM),WORK(KT2AM))
      end if
      skod = DDOT(NT1AM(ISYML),WORK(KL1AM),1,WORK(KL1AM),1)
      skod = skod+DDOT(NT2AM(ISYML),WORK(KT2AM),1,WORK(KT2AM),1)
      WRITE(LUPRI,1) 'Norm of Tbar (full):       ',skod
      IF (.NOT. CCS) CALL CC_T2SQ(WORK(KT2AM),WORK(KL2AM),ISYML)
C
C--------------------------------
C     Put T2 (packed) amplitudes in etac2.
C--------------------------------
C
      IF (.NOT. CCS) THEN
         !read in T2AM (packed)
         IOPT = 2
         CALL CC_RDRSP('R0',0,1,IOPT,MODEL,DUMMY,WORK(KT2AM))
         skod = DDOT(NT2AM(1),WORK(KT2AM),1,WORK(KT2AM),1)
          WRITE(LUPRI,1) 'Norm of T2(0) (full):   ',skod
      ENDIF
C
C--------------------------------
C     Make X and Y intermediates.
C--------------------------------
C
      IF (.NOT. CCS) THEN
         KXMAT = KEND21
         KYMAT = KXMAT + NMATIJ(ISYML)
         KEND3 = KYMAT + NMATAB(ISYML)
         LEND3 = LWORK - KEND3
         IF (LEND3.LT. 0 )
     &        CALL QUIT(' TOO LITTLE WORKSPACE IN CC_ETAC-2')
C
         !IF ( DEBUG ) THEN
            XYI   = DDOT(NT1AM(ISYML),WORK(KL1AM),1,WORK(KL1AM),1)
            WRITE(LUPRI,1) 'CC_ETAC: L1AM vector:              ',XYI
            XYI   = DDOT(NT2SQ(ISYML),WORK(KL2AM),1,WORK(KL2AM),1)
            WRITE(LUPRI,1) 'CC_ETAC: L2AM vector:              ',XYI
            XXI   = DDOT(NT2AM(ISYMOP),WORK(KT2AM),1,WORK(KT2AM),1)
            WRITE(LUPRI,1) 'T2AM vector :                      ',XXI
         !ENDIF
         CALL CC_XI(WORK(KXMAT),WORK(KL2AM),ISYML,WORK(KT2AM),1,
     *              WORK(KEND3),LEND3)
         CALL CC_YI(WORK(KYMAT),WORK(KL2AM),ISYML,WORK(KT2AM),1,
     *              WORK(KEND3),LEND3)
         !IF ( DEBUG ) THEN
            XYI   = DDOT(NMATAB(ISYML),WORK(KYMAT),1,WORK(KYMAT),1)
            WRITE(LUPRI,1) 'CC_ETAC: YI  intermediate is:      ',XYI
            XXI   = DDOT(NMATIJ(ISYML),WORK(KXMAT),1,WORK(KXMAT),1)
            WRITE(LUPRI,1) 'CC_ETAC: XI  intermediate is:      ',XXI
         !ENDIF
      ELSE
         KEND3 = KEND2
         LEND3 = LEND2
      ENDIF
C
C----------------------------------------------
C     Calculate X and Y contributions to etac1.
C     etac1 = -sum(e)Cie*Yae - sum(l)Cla*Xli
C----------------------------------------------
C
      IF ( (.NOT.CCS) .AND. (.NOT.(CC2.AND.IOPTCC2.EQ.1)) ) THEN
      DN = DDOT(N2BST(ISYMC),WORK(KCTMO),1,WORK(KCTMO),1)
      WRITE(LUPRI,*) 'Integral norm before CC_21EFM out norm:',DN
         CALL CC_21EFM(ETAC,WORK(KCTMO),ISYMC,WORK(KXMAT),
     *                 WORK(KYMAT),ISYML,WORK(KEND3),LEND3)
      DN = DDOT(N2BST(ISYMC),WORK(KCTMO),1,WORK(KCTMO),1)
      WRITE(LUPRI,*) 'Integral norm after CC_21EFM out norm:',DN
C
         IF ( DEBUG.or.locdbg ) THEN
            ETA1 = DDOT(NT1AM(ISYRES),ETAC(1),1,ETAC(1),1)
            WRITE(LUPRI,1) 'Norm of eta1-after X&Y cont:       ',ETA1
         ENDIF
      ENDIF

C
C------------------------------------------------
C     Workspace for T2AM and X and Y is now free.
C     etac2 = P(ab,ij)(2l(ai)*Cjb - l(aj)*c(ib))
C------------------------------------------------
C
      IF (.NOT. CCS) THEN
         CALL DZERO(ETAC(1+NT1AM(ISYRES)),NT2AM(ISYRES))

      DN = DDOT(N2BST(ISYMC),WORK(KCTMO),1,WORK(KCTMO),1)
      WRITE(LUPRI,*) 'Integral norm before CC_L1FCK out norm:',DN

         CALL CC_L1FCK(ETAC(1+NT1AM(ISYRES)),WORK(KL1AM),WORK(KCTMO),
     *                 ISYML,ISYMC,WORK(KEND2),LEND2)

      DN = DDOT(N2BST(ISYMC),WORK(KCTMO),1,WORK(KCTMO),1)
      WRITE(LUPRI,*) 'Integral norm after CC_L1FCK out norm:',DN
!     *                 ISYML,ISYMC,WORK(KEND2),LEND2)
C
         IF ( DEBUG.or.locdbg ) THEN
            ETA1 = DDOT(NT1AM(ISYRES),ETAC(1),1,ETAC(1),1)
            ETA2 = DDOT(NT2AM(ISYRES),ETAC(1+NT1AM(ISYRES)),1,
     *                  ETAC(1+NT1AM(ISYRES)),1)
            WRITE(LUPRI,1) 'Norm of eta1-after L1c cont:       ',ETA1
            WRITE(LUPRI,1) 'Norm of eta2-after L1c cont:       ',ETA2
         ENDIF
      ENDIF
C
      KEI1   = KEND2
      KEI2   = KEI1   + NEMAT1(ISYMC)
      KEND3  = KEI2   + NMATIJ(ISYMC)
      LEND3  = LWORK  - KEND3
      IF (LEND3.LT. 0 ) CALL QUIT(' TOO LITTLE WORKSPACE IN CC_ETAC-3')
C
C--------------------------------
C     Put A into E matrix format.
C--------------------------------
C
      FCKCON = .TRUE.
      ETRAN  = .FALSE.
      DN = DDOT(N2BST(ISYMC),WORK(KCTMO),1,WORK(KCTMO),1)
      WRITE(LUPRI,*) 'Integral norm before CCRHS_EFCK out norm:',DN
      CALL CCRHS_EFCK(WORK(KEI1),WORK(KEI2),WORK(KLAMDH),
     *                WORK(KCTMO),WORK(KEND3),LEND3,FCKCON,
     *                ETRAN,ISYMC)
C
C--------------------------------------------
C     etac1 =  sum(b)Lbi*Cba - sum(j)Laj*Cij.
C--------------------------------------------
C
      IF ( DEBUG.or.locdbg ) THEN
         XE1 = DDOT(NMATAB(ISYMC),WORK(KEI1),1,WORK(KEI1),1)
         XE2 = DDOT(NMATIJ(ISYMC),WORK(KEI2),1,WORK(KEI2),1)
         WRITE(LUPRI,1) 'Norm of EI1  -after EFCK:          ',XE1
         WRITE(LUPRI,1) 'Norm of EI2  -after EFCK:          ',XE2
         XE2 = DDOT(NMATIJ(ISYMC),WORK(KINTIJ),1,WORK(KINTIJ),1)
         WRITE(LUPRI,1) 'Norm of Cij integrals  -after EFCK: ',XE2
         ETA1 = DDOT(NT1AM(ISYML),WORK(KL1AM),1,WORK(KL1AM),1)
         WRITE(LUPRI,1) 'Norm of L1AM before  CCLR_E1C1:    ',ETA1
      ENDIF
C
c test
c     kei11= kend3
c     kei21= kei11+ NMATAB(ISYMC)
c     kend3 = kei21+ NMATIJ(ISYMC)
c     lend3 = lwork -kend3
c     call dzero(work(kei11),NMATAB(ISYMC))
c     call dzero(work(kei21),NMATIJ(ISYMC))
c     call dcopy(NMATAB(ISYMC),work(kei1),1,work(kei11),1)
c     call dcopy(NMATIJ(ISYMC),work(kei2),1,work(kei21),1)
c     CALL CCLR_E1C1(ETAC,WORK(KL1AM),WORK(KEI11),WORK(KEI21),
c    *               WORK(KEND3),LEND3,ISYML,ISYMC,'T')
c test
C
      CALL CCLR_E1C1(ETAC,WORK(KL1AM),WORK(KEI1),WORK(KEI2),
     *               WORK(KEND3),LEND3,ISYML,ISYMC,'T')
C
      IF ( DEBUG.or.locdbg ) THEN
         ETA1 = DDOT(NT1AM(ISYRES),ETAC(1),1,ETAC(1),1)
         WRITE(LUPRI,1) 'Norm of eta1 - after CCLR_E1C1:    ',ETA1
      ENDIF
C
C---------------------------------------------------------------
C     etac2 = P(ab,ij)(sum(e)2L(aiej)*Ceb - sym(k)L(aibk)*C(jk))
C---------------------------------------------------------------
C
      IF (.NOT. CCS) THEN
C
         IF (CC2 .AND. IOPTCC2.EQ.1) THEN
           FCKCON = .TRUE.
           ETRAN  = .FALSE.
           CALL CCRHS_EFCK(WORK(KEI1),WORK(KEI2),WORK(KCMO),
     *             WORK(KFCKHF),WORK(KEND3),LEND3,FCKCON,ETRAN,ISYMC)
         END IF

         CALL CC_EITR(WORK(KEI1),WORK(KEI2),WORK(KEND3),LEND3,
     *                ISYMC)
C
         CALL CCRHS_E(ETAC(1+NT1AM(ISYRES)),WORK(KL2AM),
     *                WORK(KEI1),WORK(KEI2),WORK(KEND3),
     *                LEND3,ISYML,ISYMC)
C
         IF (IPRINT .GT. 40 ) THEN
            CALL AROUND( 'In CC_ETAC:  EtaC vector ' )
            CALL CC_PRP(ETAC(1),ETAC(1+NT1AM(ISYRES)),ISYMC,1,1)
         ENDIF
C
         IF (DEBUG .OR. ( IPRINT .GT. 20 )) THEN
            ETA1 = DDOT(NT1AM(ISYRES),ETAC(1),1,ETAC(1),1)
            ETA2 = DDOT(NT2AM(ISYRES),ETAC(1+NT1AM(ISYRES)),1,
     *                  ETAC(1+NT1AM(ISYRES)),1)
            WRITE(LUPRI,1) 'Norm of eta1 - end of CC_ETAC:     ',ETA1
            WRITE(LUPRI,1) 'Norm of eta2 - end of CC_ETAC:     ',ETA2
            CALL AROUND( 'END OF CC_ETAC ')
         ENDIF
      ENDIF
C
      !if (cccisd.or.eomccsd.or.fdeom) then

         !restart from zero to make sure not to mess up
         KETAC = KEND21
         KT2AM = KETAC + NT1AM(ISYRES)
         KT2SQ = KT2AM + NT2AM(1)
         KL2AM = KT2SQ + NT2SQ(1)
         KL2SQ = KL2AM + NT2AM(ISYML)
         KZint = KL2SQ + NT2SQ(ISYML)
         KEND3 = KZint + NT2SQ(isyml)
         LEND3 = LWORK - KEND3
         IF (LEND3 .LT. 0 )
     &     CALL QUIT('INSUFFICIENT WORK SPACE IN CC_EOMETA')

C
C---------------------------------------------------------------
C     etac1 = 2sum(ei)L(ckei)*Cei
C     requires packed L2, so I reread it into T2AM
C     can also used it square, if last argument in called routine
C     is set = 1
C---------------------------------------------------------------
Css
      if (ioptl0.eq.1) then
        call dcopy(NT2AM(ISYML),XL0(1+NT1AM(ISYML)),1,WORK(KL2AM),1)
      else
         IOPT = 2
         CALL CC_RDRSP(LIST,ILSTNR,ISYML,IOPT,MODEL,
     *                 WORK(KEND3),WORK(KL2AM))
      end if
         !CALL DSCAL(NT1AM(isymc),TWO,WORK(KINTAI),1)
         !
!         CALL CCG_LXD(ETAC(1),ISYRES,WORK(KINTAI),ISYMC,
!     &                               WORK(KL2AM),ISYML,0)

         !call dzero(WORK(KETAC),NT1AM(ISYRES))
         !warning: result vector is zeroed inside!
         CALL CCG_LXD(WORK(KETAC),ISYRES,WORK(KINTAI),ISYMC,
     &                                   WORK(KL2AM),ISYML,0)
         IF ( DEBUG.or.locdbg ) THEN
            ETA1 = DDOT(NT1AM(ISYRES),WORK(KETAC),1,WORK(KETAC),1)
            WRITE(LUPRI,1) 'Norm of alone Tbar_ck,ei*Cei: ',ETA1
         ENDIF
         !removed factor 2 to get FCI limit!
         CALL DAXPY(NT1AM(ISYRES),1.0D0,WORK(KETAC),1,ETAC(1),1)
         ! add to total tensor
         IF ( DEBUG.or.locdbg ) THEN
            ETA1 = DDOT(NT1AM(ISYRES),ETAC(1),1,ETAC(1),1)
            WRITE(LUPRI,1)'Norm Eta^C after +Tbar_ck,ei*Cei:',ETA1
         ENDIF
C
C---------------------------------------------------------------
C     etac2 = 2(sum(nf)*Cnf(sum(bj)tbar(ck,bj)*ttilde(fn,bj))
C     the sum(bj)tbar(ck,bj)*ttilde(fn,bj) term is part of the
C     d_aijb density!!!
C---------------------------------------------------------------
C
         !I reread T2AM 
         IOPT = 2
         CALL CC_RDRSP('R0',0,1,IOPT,MODEL,DUMMY,WORK(KT2AM))
         ISYOPE = 1
         IOPTTCME = 1
         !get ttilde (overwrites T2am)
         CALL CCSD_TCMEPK(WORK(KT2AM),1.0D0,ISYOPE,IOPTTCME)
         write(lupri,*) ' norm of ttilde ', 
     &         ddot(NT2AM(1),WORK(KT2AM),1,WORK(KT2AM),1)
         !square it up
         CALL CC_T2SQ(WORK(KT2AM),WORK(KT2SQ),1)
         !Square up L2 as well
         IF (.NOT. CCS) CALL CC_T2SQ(WORK(KL2AM),WORK(KL2SQ),ISYML)

         !zero-ed inside the subroutine called
         call CC_ZIsq(work(KZint),work(kl2sq),ISYML,
     &                work(kT2sq),1,WORK(kend3),LEND3)
         write(lupri,*) ' Norm of Z intermediate ', 
     &   ddot(NT2SQ(isyml),WORK(KZINT),1,WORK(KZINT),1)

         !if (iprint.ge.45) then
         !  write(lupri,*) 'CC_EOMETA: Zint intermediate'
         !  call cc_prsq(work(kend3),WORK(KZint),1,0,1)
         !end if

         call dzero(WORK(KETAC),NT1AM(ISYRES))
         CALL CCG_LXD(WORK(KETAC),ISYRES,WORK(KINTIA),ISYMC,
     &                               WORK(KZINT),ISYML,1)

         IF ( DEBUG.or.locdbg ) THEN
            ETA1 = DDOT(NT1AM(ISYRES),WORK(KETAC),1,WORK(KETAC),1)
            !ETA1 = DDOT(NT1AM(ISYRES),ETAC(1),1,ETAC(1),1)
            WRITE(LUPRI,1) 'Norm of eta1 (alone) Zint*C: ',ETA1
         ENDIF
         !removed factor 2 to get FCI limit
         CALL DAXPY(NT1AM(ISYRES),1.0D0,WORK(KETAC),1,ETAC(1),1)  
         !
         IF ( DEBUG.or.locdbg ) THEN
            ETA1 = DDOT(NT1AM(ISYRES),ETAC(1),1,ETAC(1),1)
            WRITE(LUPRI,1) 'Norm of eta after 2*Zint*Cei: ',ETA1
         ENDIF

      !end if
      IF (IPRINT .GT. 5 ) THEN
         TIMEC = SECOND() - TIMEC
         WRITE(LUPRI,9999) 'CCCI_ETA^C      ', TIMEC
      ENDIF
C
   1  FORMAT(1x,A35,1X,E20.10)
9999  FORMAT(1x,'Time used in',2x,A18,2x,': ',f10.2,' seconds')
C
8888  CONTINUE
      END

C  /* Deck cc_zisq */
      SUBROUTINE CC_ZIsq(ZMAT,CLTR2,ISYMLV,T2AM,ISYMRV,
     &                   WORK,LWORK)
C
C     Written by Sonia Coriani Halkier 27/04 - 2015
C
C     Version: 1.0
C
C     Purpose: To calculate the Z-intermediate entering the EOM
C              etaC.
C
C     It is assumed that both CLTR2 and T2am are squared up
C
C     General symmetry of both CLTR2 and T2AM for F-mat and etaC.
C
#include "implicit.h"
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .false.)
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION ZMAT(*),CLTR2(*),T2AM(*),WORK(LWORK)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "cclr.h"
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      ISYRES = MULD2H(ISYMLV,ISYMRV)
C
C-----------------------------
C     Initialize result array.
C-----------------------------
C
      CALL DZERO(ZMAT,NT2SQ(ISYRES))
C
      DO 100 ISYMCK = 1,NSYM
C
         ISYMBJ = MULD2H(ISYMCK,ISYMLV)  !left vector (tbar)
         ISYMFN = MULD2H(ISYMBJ,ISYMRV)  !right vector (ttilde)
C
C-------------------------------------
C        Contract the two matrices.
C-------------------------------------------
C
         KOFF1  = IT2SQ(ISYMCK,ISYMBJ) + 1
         KOFF2  = IT2SQ(ISYMBJ,ISYMFN) + 1
         KOFF3  = IT2SQ(ISYMCK,ISYMFN) + 1
C
         NTOTCK = MAX(NT1AM(ISYMCK),1)
         NTOTBJ = MAX(NT1AM(ISYMBJ),1)
         NTOTFN = MAX(NT1AM(ISYMFN),1)
C
         CALL DGEMM('N','N',NTOTCK,NTOTFN,NTOTBJ,
     &               ONE,CLTR2(KOFF1),NTOTCK,T2AM(KOFF2),NTOTFN,
     &               ONE,ZMAT(KOFF3),NTOTCK)
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      IF (LOCDBG) THEN
        WRITE(LUPRI,*) '>>>>>CC_ZI: Norm of Z-intermediate:',
     &    DDOT(NT2SQ(ISYRES),ZMAT,1,ZMAT,1)
      END IF

!      DO 100 ISYMCK = 1,NSYM
C
!         ISYMBJ = MULD2H(ISYMCK,ISYMLV)  !left vector (tbar)
!         ISYMFN = MULD2H(ISYMBJ,ISYMRV)  !right vector (ttilde)
C
C-------------------------------------
C        Contract the two matrices.
C-------------------------------------------
C
!         KOFF1  = IMATAB(ISYME,ISYMF) + 1
C
!         NTOTCK = MAX(NT1AM(ISYMCK),1)
!         NTOTBJ = MAX(NT1AM(ISYMBJ),1)
!         NTOTFN = MAX(NT1AM(ISYMFN),1)
C
!         CALL DGEMM('N','N',NTOTCK,NTOTFN,NT1AM(ISYMBJ),
!     &               ONE,CLTR2,NTOTCK,KT2AM,NTOTFN,
!     &               ONE,ZMAT,NTOTBJ)
C
!  120       CONTINUE
!  110    CONTINUE
!  100 CONTINUE
C
      RETURN
      END
