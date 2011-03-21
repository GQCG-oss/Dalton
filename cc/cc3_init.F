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
C  /* Deck get_lambdax */
      SUBROUTINE GET_LAMBDAX(LAMPX,LAMHX,LISTX,IDLSTX,ISYMX,LAMP0,LAMH0,
     *                       WORK,LWORK)
*
***********************************************************************
*                                                                     *
*     Construct Lambda^X matrices.                                    *
*                                                                     *
*     X may actually be X,Y,Z,U,..., etc. operator.                   *
*                                                                     *
*     Filip Pawlowski, 05-Dec-2002, Aarhus                            *
*                                                                     *
***********************************************************************
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccr1rsp.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <dummy.h>
C
      INTEGER IDLSTX,ISYMX
      INTEGER KT1X
      INTEGER LWORK,KEND1,LWRK1
      INTEGER IOPT
C
      CHARACTER LISTX*3, MODEL*10
C
#if defined (SYS_CRAY) 
      REAL LAMPX(*),LAMHX(*),LAMP0(*),LAMH0(*),WORK(LWORK)
#else
      DOUBLE PRECISION LAMPX(*),LAMHX(*),LAMP0(*),LAMH0(*),WORK(LWORK)
#endif
C
      CALL QENTER('GET_LAMBDAX')
C
C----------------------------------------
C     Initial test of symmetry consitency
C----------------------------------------
C
C     IF (ISYMX .NE. ISYLRT(IDLSTX)) THEN
C        WRITE(LUPRI,*) 'ISYLRT(IDLSTX) = ',ISYLRT(IDLSTX)
C        WRITE(LUPRI,*) 'ISYMX = ',ISYMX
C        CALL QUIT('Symmetry mismatch for operator in GET_LAMBDAX')
C     END IF
C
C-------------------
C     Initialization
C-------------------
C
      IOPT = 1
C
C----------------
C     Allocation
C----------------
C
      KT1X  = 1
      KEND1 = KT1X  + NT1AM(ISYMX)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_LAMBDAX (1)')
      END IF
C
      CALL DZERO(WORK(KT1X),NT1AM(ISYMX))
C
C-----------------
C     Get Lambda_X
C-----------------
C
      CALL CC_RDRSP(LISTX,IDLSTX,ISYMX,IOPT,MODEL,WORK(KT1X),DUMMY)
C
      CALL CCLR_LAMTRA(LAMP0,LAMPX,LAMH0,LAMHX,WORK(KT1X),ISYMX)
C
C-------------
C     End
C-------------
C
      CALL QEXIT('GET_LAMBDAX')
C
      RETURN
      END
C  /* Deck get_fock0 */
      SUBROUTINE GET_FOCK0(FOCK0,LAMP0,LAMH0,WORK,LWORK)
*
***********************************************************************
*                                                                     *
*     Read in Fock matrix in AO basis (from file) and transform to    *
*     Lambda_0 basis.                                                 *
*                                                                     *
*     Could be generalized to get MO-transformed Fock matrix or       *
*     f.ex. Fa,i-bar matrix, but we do not need that in CC3.          *
*                                                                     *
*     Filip Pawlowski, 05-Dec-2002, Aarhus                            *
*                                                                     *
***********************************************************************
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <dummy.h>
C
      INTEGER ISYM0,LWORK,LUFOCK
C
#if defined (SYS_CRAY) 
      REAL FOCK0(*),LAMP0(*),LAMH0(*),WORK(LWORK)
#else
      DOUBLE PRECISION FOCK0(*),LAMP0(*),LAMH0(*),WORK(LWORK)
#endif
C
      CALL QENTER('GET_FOCK0')
C
C-------------------------
C     Some initializations
C-------------------------
C
      ISYM0  = 1
C
C-----------------------------------------------
C     Read zeroth-order AO Fock matrix from file
C-----------------------------------------------
C
      LUFOCK = -1
      CALL GPOPEN(LUFOCK,'CC_FCKH','OLD',' ','UNFORMATTED',IDUMMY,
     *            .FALSE.)
      REWIND(LUFOCK)
      READ(LUFOCK) (FOCK0(I),I=1,N2BST(ISYM0))
      CALL GPCLOSE(LUFOCK,'KEEP')
C
C-------------------------------------------
C     Trasform Fock matrix to Lambda_0 basis
C-------------------------------------------
C
      CALL CC_FCKMO(FOCK0,LAMP0,LAMH0,WORK,LWORK,ISYM0,ISYM0,ISYM0)
C
C-------------
C     End
C-------------
C
      CALL QEXIT('GET_FOCK0')
C
      RETURN
      END
C  /* Deck get_fockx */
      SUBROUTINE GET_FOCKX(FOCKX,LABELX,IDLSTX,ISYMX,LAMPY,ISYMPY,
     *                         LAMHZ,ISYMHZ,WORK,LWORK)
*
***********************************************************************
*                                                                     *
*     Read in matrix with property integrals (FOCKX) in AO basis      *
*     and transform using Lambda_pY and Lambda_hZ.
*                                                                     *
*     Poul Jorgensen, Filip Pawlowski, Spring-2003, Aarhus
*                                                                     *
***********************************************************************
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccl1rsp.h>
#include <ccorb.h>
C
      INTEGER ISYMX,IRREP,ISYM
      INTEGER LWORK,KEND1,LWRK1
      INTEGER IERR
      INTEGER IDLSTX
      INTEGER ISYMPY,ISYMHZ
C
      INTEGER ISYMYZ,ISYMXYZ,LENFCKX,KFOCKX
C
      CHARACTER LABELX*8
C
#if defined (SYS_CRAY) 
      REAL FOCKX(*),LAMPY(*),LAMHZ(*),WORK(LWORK)
#else
      DOUBLE PRECISION FOCKX(*),LAMPY(*),LAMHZ(*),WORK(LWORK)
#endif
C
      CALL QENTER('FCKX')
C
      ISYMYZ  = MULD2H(ISYMPY,ISYMHZ)
      ISYMXYZ = MULD2H(ISYMX,ISYMYZ)
C
      LENFCKX = MAX(N2BST(ISYMX),N2BST(ISYMXYZ))
C
      KFOCKX = 1
      KEND1  = KFOCKX + LENFCKX
      LWRK1  = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*)'Memory available: LWORK = ', LWORK
         WRITE(LUPRI,*)'Memory needed:    KEND1 = ', KEND1
         CALL QUIT('Insufficient memory in GET_FOCKX')
      END IF
C
C--------------------------------------
C     Read property integrals (FOCKX) in AO basis from file
C--------------------------------------
C
      CALL CCPRPAO(LABELX,.TRUE.,WORK(KFOCKX),IRREP,ISYM,IERR,
     *             WORK(KEND1),LWRK1)
C
      IF ((IERR.GT.0) .OR. (IERR.EQ.0 .AND. IRREP.NE.ISYMX)) THEN
          CALL QUIT('GET_FOCKX: error reading oper. '//LABELX)
      ELSE IF (IERR.LT.0) THEN
          CALL DZERO(WORK(KFOCKX),LENFCKX)
      END IF
C
C---------------------------------------------------
C     Transform property integrals (FOCKX) using Lambda_pY and Lambda_hZ
C---------------------------------------------------
C
      CALL CC_FCKMO(WORK(KFOCKX),LAMPY,LAMHZ,WORK(KEND1),LWRK1,
     *              ISYMX,ISYMPY,ISYMHZ)
C
C---------------------------------
C     Copy result to FOCKX array
C---------------------------------
C
      CALL DCOPY(N2BST(ISYMXYZ),WORK(KFOCKX),1,FOCKX,1)
C
C-------------
C     End
C-------------
C
      CALL QEXIT('FCKX')
C
      RETURN
      END
C  /* Deck get_lambda0 */
      SUBROUTINE GET_LAMBDA0(LAMP0,LAMH0,WORK,LWORK)
*
***********************************************************************
*                                                                     *
*     Construct Lambda^0 matrices.                                    *
*                                                                     *
*     Filip Pawlowski, 05-Dec-2002, Aarhus                            *
*                                                                     *
***********************************************************************
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccr1rsp.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <dummy.h>
C
      INTEGER IDLST0,ISYM0
      INTEGER KT1AMP0
      INTEGER LWORK,KEND1,LWRK1
      INTEGER IOPT
C
      CHARACTER LIST0*2, MODEL*10
C
#if defined (SYS_CRAY) 
      REAL LAMP0(*),LAMH0(*),WORK(LWORK)
#else
      DOUBLE PRECISION LAMP0(*),LAMH0(*),WORK(LWORK)
#endif
C
      CALL QENTER('GET_LAMBDA0')
C
C-------------------------
C     Some initializations
C-------------------------
C
      IOPT   = 1
      LIST0  = 'R0'
      IDLST0 = 0
      ISYM0  = 1
C
C----------------
C     Allocation
C----------------
C
      KT1AMP0 = 1
      KEND1   = KT1AMP0 + NT1AM(ISYM0)
      LWRK1   = LWORK   - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in GET_LAMBDA0')
      END IF
C
      CALL DZERO(WORK(KT1AMP0),NT1AM(ISYM0))
C
      CALL CC_RDRSP(LIST0,IDLST0,ISYM0,IOPT,MODEL,WORK(KT1AMP0),DUMMY)
C
      CALL LAMMAT(LAMP0,LAMH0,WORK(KT1AMP0),WORK(KEND1),LWRK1)
C
C-------------
C     End
C-------------
C
      CALL QEXIT('GET_LAMBDA0')
C
      RETURN
      END
C  /* Deck get_t1_t2 */
      SUBROUTINE GET_T1_T2(IOPT,DIASCL,T1,T2,LIST,IDLST,ISYM,WORK,LWORK)
*
***********************************************************************
*                                                                     *
*     Get T1 and/or T2 amplitudes/multipliers for a given list.       *
*                                                                     *
*     IOPT = 1 :  Get only T1                                         *
*     IOPT = 2 :  Get only T2                                         *
*     IOPT = 3 :  Get both T1 and T2                                  *
*                                                                     *
*     Filip Pawlowski, 06-Dec-2002, Aarhus                            *
*                                                                     *
***********************************************************************
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccr1rsp.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <dummy.h>
C
      INTEGER IDLST,ISYM,LWORK,IOPT
C
      CHARACTER LIST*(*), MODEL*10
C
      LOGICAL DIASCL
C
#if defined (SYS_CRAY) 
      REAL T1(*),T2(*),WORK(LWORK)
      REAL TWO
#else
      DOUBLE PRECISION T1(*),T2(*),WORK(LWORK)
      DOUBLE PRECISION TWO
#endif
C
      PARAMETER (TWO = 2.0D0)
C
      CALL QENTER('GET_T1_T2')
C
C-------------------------
C     Initial test of IOPT
C-------------------------
C
      IF ((IOPT.NE.1).AND.(IOPT.NE.2).AND.(IOPT.NE.3)) THEN
         WRITE(LUPRI,*)'IOPT = ', IOPT
         CALL QUIT('Illegal option in GET_T1_T2')
      END IF
C
C-------------------------------------------------
C     If iopt=1 get only T1 amplitudes/multipliers 
C-------------------------------------------------
C
      IF (IOPT.EQ.1) THEN
C
         CALL CC_RDRSP(LIST,IDLST,ISYM,IOPT,MODEL,T1,DUMMY)
C
      END IF
C
C-------------------------------------------------
C     If iopt=2 get only T2 amplitudes/multipliers 
C-------------------------------------------------
C
      IF (IOPT.EQ.2) THEN
C
         CALL CC_RDRSP(LIST,IDLST,ISYM,IOPT,MODEL,DUMMY,T2)
C
         IF (DIASCL) THEN
            CALL CCLR_DIASCL(T2,TWO,ISYM)
         END IF
C
         CALL CC_T2SQ(T2,WORK,ISYM)
C
         CALL CC3_T2TP(T2,WORK,ISYM)
C
      END IF
C
C---------------------------------------------------
C     If iopt=3 get T1 and T2 amplitudes/multipliers 
C---------------------------------------------------
C
      IF (IOPT.EQ.3) THEN
C
         CALL CC_RDRSP(LIST,IDLST,ISYM,IOPT,MODEL,T1,T2)
C
         IF (DIASCL) THEN
            CALL CCLR_DIASCL(T2,TWO,ISYM)
         END IF
C
         CALL CC_T2SQ(T2,WORK,ISYM)
C
         CALL CC3_T2TP(T2,WORK,ISYM)
C
      END IF
C
C-------------
C     End
C-------------
C
      CALL QEXIT('GET_T1_T2')
C
      RETURN
      END
C  /* Deck sort_fockck */
      SUBROUTINE SORT_FOCKCK(FOCKCK,FOCK0,ISYFCK)
*
***********************************************************************
*                                                                     *
*     Sort the Fock matrix to get F(ck) block                         *
*                                                                     *
*     Filip Pawlowski, 06-Dec-2002, Aarhus                            *
*                                                                     *
***********************************************************************
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYFCK,ISYMC,ISYMK,KOFF1,KOFF2
C
#if defined (SYS_CRAY) 
      REAL FOCK0(*),FOCKCK(*)
#else
      DOUBLE PRECISION FOCK0(*),FOCKCK(*)
#endif
C
      CALL QENTER('SORT_FOCKCK')
C
      DO ISYMC = 1,NSYM
C
         ISYMK = MULD2H(ISYMC,ISYFCK)
C
         DO K = 1,NRHF(ISYMK)
C
            DO C = 1,NVIR(ISYMC)
C
               KOFF1 = IFCVIR(ISYMK,ISYMC) +
     *                 NORB(ISYMK)*(C - 1) + K
               KOFF2 = IT1AM(ISYMC,ISYMK)
     *               + NVIR(ISYMC)*(K - 1) + C
C
               FOCKCK(KOFF2) = FOCK0(KOFF1)
C
            ENDDO ! C
         ENDDO    ! K
      ENDDO       ! ISYMC
C
      CALL QEXIT('SORT_FOCKCK')
C
      RETURN
      END
C  /* Deck get_orben */
      SUBROUTINE GET_ORBEN(FOCKD,WORK,LWORK)
*
***********************************************************************
*                                                                     *
*     Get orbital energies.                                           *
*                                                                     *
*     Filip Pawlowski, 06-Dec-2002, Aarhus                            *
*                                                                     *
***********************************************************************
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <dummy.h>
#include <inftap.h>
#include <ccorb.h>
#include <ccsdinp.h>
C
      INTEGER LWORK
C
#if defined (SYS_CRAY) 
      REAL FOCKD(*),WORK(LWORK)
#else
      DOUBLE PRECISION FOCKD(*),WORK(LWORK)
#endif
C
      CALL QENTER('GET_ORBEN')
C
C-------------------------------------
C     Read canonical orbital energies:
C-------------------------------------
C
      LUSIFC = -1
C
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     *            .FALSE.)
      REWIND LUSIFC
C 
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (FOCKD(I), I=1,NORBTS)
C
      CALL GPCLOSE(LUSIFC,'KEEP')
C
C---------------------------------------------
C     Delete frozen orbitals in Fock diagonal.
C---------------------------------------------
C
      IF (FROIMP .OR. FROEXP)
     *   CALL CCSD_DELFRO(FOCKD,WORK,LWORK)
C
C-------------
C     End
C-------------
C
      CALL QEXIT('GET_ORBEN')
C
      RETURN
      END
C  /* Deck intocc_t3bar0 */
      SUBROUTINE INTOCC_T3BAR0(LUTOC,FNTOC,LAMH0,ISYINT,T3BOG1,T3BOL1,
     *                         T3BOG2,T3BOL2,WORK,LWORK)
*
***********************************************************************
*                                                                     *
*     Construct occupied integrals which are required to calculate    *
*     t3bar_0 multipliers                                             *
*                                                                     *
*     Construct 2*C-E of the integrals.                               *
*     Have integral for both (ij,k,a) and (a,k,j,i)                   *
*                                                                     *
*                                                                     *
*     Filip Pawlowski, 09-Dec-2002, Aarhus                            *
*                                                                     *
***********************************************************************
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
C
      INTEGER ISYINT,LWORK,KINTOC,KEND1,LWRK1,IOFF
      INTEGER LUTOC
C
      CHARACTER*(*) FNTOC
C
#if defined (SYS_CRAY) 
      REAL LAMH0(*),T3BOG1(*),T3BOL1(*),T3BOG2(*),T3BOL2(*),WORK(LWORK)
#else
      DOUBLE PRECISION LAMH0(*),T3BOG1(*),T3BOL1(*),T3BOG2(*),T3BOL2(*)
      DOUBLE PRECISION WORK(LWORK)
#endif
C
      CALL QENTER('INTOCC_T3BAR0')
C
C---------------
C     Allocation
C---------------
C
      KINTOC = 1
      KEND1  = KINTOC + NTOTOC(ISYINT)
      LWRK1  = LWORK  - KEND1
C      
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in INTOCC_T3BAR0')
      END IF
C
      CALL DZERO(WORK(KINTOC),NTOTOC(ISYINT))
C
C-----------------------------------------------------------
C     Construct 2*C-E of the integrals.
C-----------------------------------------------------------
C
      IOFF = 1
C    ! Read in integrals (ia|j delta) sitting as KINTOC(ah ip | jp delta)
      IF (NTOTOC(ISYINT) .GT. 0) THEN
         CALL GETWA2(LUTOC,FNTOC,WORK(KINTOC),IOFF,NTOTOC(ISYINT))
      ENDIF
C
C    ! Transform the integrals with hole index to (ia|jk) and sort as 
C    ! T3BOG1(kija)
      CALL CC3_TROCC(WORK(KINTOC),T3BOG1,LAMH0,WORK(KEND1),LWRK1,ISYINT)

C     ! Calculate 2*(ia|jk) - (ka|ji) sorted as T3BOL1(kija) 
      CALL CCSDT_TCMEOCC(T3BOG1,T3BOL1,ISYINT)
C
C    ! (ia|jk) sorted as T3BOG1(kija) and now sorted as T3BOG2(ajik)
      CALL CCFOP_SORT(T3BOG1,T3BOG2,ISYINT,1)
C
      CALL CCFOP_SORT(T3BOL1,T3BOL2,ISYINT,1)
C
C-------------
C     End
C-------------
C
      CALL QEXIT('INTOCC_T3BAR0')
C
      RETURN
      END
C  /* Deck intocc_t30 */
      SUBROUTINE INTOCC_T30(LUCKJD,FNCKJD,LAMP0,ISYINT,T3OG1,T3OG2,
     *                      WORK,LWORK)
*
***********************************************************************
*                                                                     *
*     Construct occupied integrals which are required to calculate    *
*     t3_0 amplitudes (in terms of SMAT and UMAT).                    *
*                                                                     *
*                                                                     *
*     Filip Pawlowski, 09-Dec-2002, Aarhus                            *
*                                                                     *
***********************************************************************
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
C
      INTEGER ISYINT,LWORK,KINTOC,KEND1,LWRK1,IOFF
      INTEGER LUCKJD
C
      CHARACTER*(*) FNCKJD
C
#if defined (SYS_CRAY) 
      REAL LAMP0(*),T3OG1(*),T3OG2(*),WORK(LWORK)
#else
      DOUBLE PRECISION LAMP0(*),T3OG1(*),T3OG2(*),WORK(LWORK)
#endif
C
      CALL QENTER('INTOCC_T30')
C
C---------------
C     Allocation
C---------------
C
      KINTOC = 1
      KEND1  = KINTOC + NTOTOC(ISYINT)
      LWRK1  = LWORK  - KEND1
C      
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in INTOCC_T30')
      END IF
C
      CALL DZERO(WORK(KINTOC),NTOTOC(ISYINT))
C
C------------------------
C     Occupied integrals for S3MAT and U3MAT.
C-----------------------
C
      IOFF = 1
      IF (NTOTOC(ISYINT) .GT. 0) THEN
         CALL GETWA2(LUCKJD,FNCKJD,WORK(KINTOC),IOFF,NTOTOC(ISYINT))
      ENDIF
C
C----------------------------------------------------------------------
C     Transform (ai| delta j) integrals to (ai|kj) and sort as I(k,i,j,a)
C----------------------------------------------------------------------
C     here the xlamdp transformation need to be used 
C     (delta is particle index) 
C
      CALL CC3_TROCC(WORK(KINTOC),T3OG1,LAMP0,WORK(KEND1),LWRK1,ISYINT)

C
C----------------------------------------------------------------------
C     (ai|kj) sorted as T3OG1(kija) becomes T3OG2(ajik) needed for U3MAT
C----------------------------------------------------------------------
C
      CALL CCFOP_SORT(T3OG1,T3OG2,ISYINT,1)
C
C-------------
C     End
C-------------
C
      CALL QEXIT('INTOCC_T30')
C
      RETURN
      END
C  /* Deck intocc_t3barx */
      SUBROUTINE INTOCC_T3BARX(LSKIPL1R,
     *                         LUTOC,FNTOC,ISYHAM,LAMH0,ISYM0,ISYINT0,
     *                         LAMHY,ISYMY,ISYINTY,W3BXOG1,W3BXOL1,
     *                         W3BXOGX1,W3BXOLX1,WORK,LWORK)
*
***********************************************************************
*                                                                     *
*     Construct occupied integrals which are required to calculate    *
*     t3bar_X multipliers:                                            *
*                                                                     *
*     g(ia|j k-)    and     L(ia|j k-)                                *
*                                                                     *
*     Both LambdaH_X and LambdaH_0 transformed integrals are needed   *
*                                                                     *
*     ISYHAM  - symmetry of Hamiltonian                               *
*     ISYM0   - symmetry of operator in LambdaH_0 transformation      *
*     ISYINT0 - symmetry of LambdaH_0-transformed integrals           *
*     ISYMY   - symmetry of operator in LambdaH_Y transformation      *
*     ISYINTY - symmetry of LambdaH_Y-transformed integrals           *
*                                                                     *
*                                                                     *
*     Filip Pawlowski, 09-Dec-2002, Aarhus                            *
*                                                                     *
***********************************************************************
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
C
      LOGICAL LSKIPL1R
C
      INTEGER ISYHAM,ISYM0,ISYINT0,ISYMY,ISYINTY,LWORK,LUTOC
      INTEGER KINTOC,KEND1,LWRK1,IOFF
C
      CHARACTER*(*) FNTOC
C
#if defined (SYS_CRAY) 
      REAL LAMH0(*),LAMHY(*),W3BXOG1(*),W3BXOL1(*),W3BXOGX1(*)
      REAL W3BXOLX1(*),WORK(LWORK)
#else
      DOUBLE PRECISION LAMH0(*),LAMHY(*),W3BXOG1(*),W3BXOL1(*)
      DOUBLE PRECISION W3BXOGX1(*),W3BXOLX1(*)
      DOUBLE PRECISION WORK(LWORK)
#endif
C
      CALL QENTER('INTOCC_T3BARX')
C
C---------------
C     Allocation
C---------------
C
      KINTOC = 1
      KEND1  = KINTOC + NTOTOC(ISYHAM)
      LWRK1  = LWORK  - KEND1
C      
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in INTOCC_T3BARX2')
      END IF
C
      CALL DZERO(WORK(KINTOC),NTOTOC(ISYHAM))
C
C-------------------------------
C     Read in occupied integrals 
C-------------------------------
C
      IOFF = 1
      IF (NTOTOC(ISYHAM) .GT. 0) THEN
         CALL GETWA2(LUTOC,FNTOC,WORK(KINTOC),IOFF,NTOTOC(ISYHAM))
      ENDIF
C
C-----------------------------
C     LambdaH_0 transformation
C-----------------------------
C
      CALL INTOCC_T3BARX2(WORK(KINTOC),LAMH0,ISYM0,ISYINT0,W3BXOG1,
     *                    W3BXOL1,WORK(KEND1),LWRK1)
C
C-----------------------------
C     LambdaH_Y transformation
C-----------------------------
C
      IF (.NOT. LSKIPL1R) THEN
        CALL INTOCC_T3BARX2(WORK(KINTOC),LAMHY,ISYMY,ISYINTY,W3BXOGX1,
     *                      W3BXOLX1,WORK(KEND1),LWRK1)
      END IF
C
C-------------
C     End
C-------------
C
      CALL QEXIT('INTOCC_T3BARX')
C
      RETURN
      END
C  /* Deck intocc_t3barx2 */
      SUBROUTINE INTOCC_T3BARX2(OCCINT,LAMH,ISYMX,ISYINT,TROCCG,TROCCL,
     *                          WORK,LWORK)
***********************************************************************
*                                                                     *
*     Construct occupied integrals which are required to calculate    *
*     t3bar_X multipliers.                                            *
*                                                                     *
*     ISYMX  - symmetry of an operator used in the transformation     *
*     ISYINT - symmetry of transformed integrals                      *
*                                                                     *
*                                                                     *
*     Filip Pawlowski, 09-Dec-2002, Aarhus                            *
*                                                                     *
***********************************************************************
C
      IMPLICIT NONE
C
#include <priunit.h>
C
      INTEGER ISYMX,ISYINT,LWORK
C
#if defined (SYS_CRAY) 
      REAL OCCINT(*),LAMH(*),TROCCG(*),TROCCL(*),WORK(LWORK)
#else
      DOUBLE PRECISION OCCINT(*),LAMH(*),TROCCG(*),TROCCL(*),WORK(LWORK)
#endif
C
      CALL QENTER('INTOCC_T3BARX2')
C
C----------------------------------------------------------------------
C     Transform (ia|j delta) integrals to (ia|j k-) and sort  G(j,i,k-,a)
C----------------------------------------------------------------------
C
      CALL CCLR_TROCC(OCCINT,TROCCG,LAMH,ISYMX,WORK,LWORK)
C
C----------------------------------------------------------------------
C integrals g(ia|j k-) sorted as G(j,i,k-,a) resorted as G(k-,i,j,a)
C                       on input KTROCCG   on output     KTROCCG
C           L(ia|j k-)                                   KTROCCL
C----------------------------------------------------------------------
C
      CALL CC3_SRTOCC(TROCCG,TROCCL,ISYINT)
C
C
C-------------
C     End
C-------------
C
      CALL QEXIT('INTOCC_T3BARX2')
C
      RETURN
      END
C  /* Deck intocc_t3x */
      SUBROUTINE INTOCC_T3X(LUCKJDR,FNCKJDR,LAMP0,ISYINT,T3OGX,WORK,
     *                      LWORK)
*
***********************************************************************
*                                                                     *
*     Besides two occupied integrals constructed in INTOCC_T30        *
*     routine the additional integrals T3OGX are needed to calculate  *
*     t3_X amplitudes.                                                *
*     Those additional ocupied integrals (T3OGX) are constructed here.*
*                                                                     *
*                                                                     *
*     Filip Pawlowski, 06-Jan-2003, Aarhus                            *
*                                                                     *
***********************************************************************
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
C
      INTEGER ISYINT,LWORK,KINTOC,KEND1,LWRK1,IOFF
      INTEGER LUCKJDR
C
      CHARACTER*(*) FNCKJDR
C
#if defined (SYS_CRAY) 
      REAL LAMP0(*),T3OGX(*),WORK(LWORK)
#else
      DOUBLE PRECISION LAMP0(*),T3OGX(*),WORK(LWORK)
#endif
C
      CALL QENTER('INTOCC_T3X')
C
C---------------
C     Allocation
C---------------
C
      KINTOC = 1
      KEND1  = KINTOC + NTOTOC(ISYINT)
      LWRK1  = LWORK  - KEND1
C      
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in INTOCC_T3X')
      END IF
C
      CALL DZERO(WORK(KINTOC),NTOTOC(ISYINT))
C
C------------------------
C     Occupied integrals for S3MAT and U3MAT.
C-----------------------
C
      IOFF = 1
      IF (NTOTOC(ISYINT) .GT. 0) THEN
         CALL GETWA2(LUCKJDR,FNCKJDR,WORK(KINTOC),IOFF,NTOTOC(ISYINT))
      ENDIF
C
C----------------------------------------------------------------------
C     Transform (ai|j delta) integrals to (ai|j k) and sort as (ij,k,a)
C----------------------------------------------------------------------
C     here the xlamdp transformation need to be used 
C     (delta is particle index) 
C
      CALL CC3_TROCC(WORK(KINTOC),T3OGX,LAMP0,WORK(KEND1),LWRK1,ISYINT)
C
C-------------
C     End
C-------------
C
      CALL QEXIT('INTOCC_T3X')
C
      RETURN
      END
C  /* Deck intvir_t30_d */
      SUBROUTINE INTVIR_T30_D(LUDKBC,FNDKBC,LUDELD,FNDELD,ISYINT,T3VDG1,
     *                        T3VDG2,T3VDG3,LAMH,ISYMD,D,WORK,LWORK)
*
***********************************************************************
*                                                                     *
*     Construct virtual integrals which are required to calculate     *
*     t3_0 amplitudes (in terms of SMAT and UMAT).                    *
*                                                                     *
*     ISYINT - symmetry of integrals                                  *
*     ISYMD  - symmetry of external (fixed) D index)                  *
*                                                                     *
*     ****************************************************            *
*     !!! THIS ROUTINE IS TO BE USED INSIDE THE D-LOOP !!!            * 
*     ****************************************************            *
*                                                                     *
*                                                                     *
*     Filip Pawlowski, 09-Dec-2002, Aarhus                            *
*                                                                     *
***********************************************************************
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INTEGER ISYINT,ISYMD,LWORK,KINTVI,KEND1,LWRK1,IOFF,ISYCKB
      INTEGER LUDKBC,LUDELD
C
      CHARACTER*(*) FNDKBC,FNDELD
C
#if defined (SYS_CRAY) 
      REAL T3VDG1(*),T3VDG2(*),T3VDG3(*),LAMH(*),WORK(LWORK)
#else
      DOUBLE PRECISION T3VDG1(*),T3VDG2(*),T3VDG3(*),LAMH(*),WORK(LWORK)
#endif
C
      CALL QENTER('INTVIR_T30_D')
C
      ISYCKB = MULD2H(ISYINT,ISYMD)
C
C---------------
C     Allocation
C---------------
C
      KINTVI = 1
      KEND1  = KINTVI + NCKA(ISYCKB)
      LWRK1  = LWORK  - KEND1
C      
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in INTVIR_T30_D')
      END IF
C
      CALL DZERO(WORK(KINTVI),NCKA(ISYCKB))
C
C-----------------------------------------------------
C     Read virtual integrals used in first S3MAT
C-----------------------------------------------------
C
      IOFF = ICKBD(ISYCKB,ISYMD) + NCKATR(ISYCKB)*(D - 1) + 1
      IF (NCKATR(ISYCKB) .GT. 0) THEN
         CALL GETWA2(LUDKBC,FNDKBC,T3VDG1,IOFF,NCKATR(ISYCKB))
      ENDIF
C
C-----------------------------------------------------------
C     Sort the integrals for S3MAT
C-----------------------------------------------------------
C
      CALL CCSDT_SRTVIR(T3VDG1,T3VDG2,WORK(KEND1),LWRK1,ISYMD,ISYINT)
C
C------------------------------------------------------
C     Read virtual integrals used in first U3MAT.
C------------------------------------------------------
C
      IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D - 1) + 1
      IF (NCKA(ISYCKB) .GT. 0) THEN
         CALL GETWA2(LUDELD,FNDELD,WORK(KINTVI),IOFF,NCKA(ISYCKB))
      ENDIF
C
      CALL CCSDT_TRVIR(WORK(KINTVI),T3VDG3,LAMH,ISYMD,D,ISYINT,
     *                 WORK(KEND1),LWRK1)
C
C-------------
C     End
C-------------
C
      CALL QEXIT('INTVIR_T30_D')
C
      RETURN
      END
C  /* Deck intvir_t3bar0_d */
      SUBROUTINE INTVIR_T3BAR0_D(LU3FOPX,FN3FOPX,LU3FOP2X,FN3FOP2X,
     *                           LUDKBC3,FNDKBC3,LU3VI,FN3VI,ISYINT,
     *                           T3BVDL1,T3BVDG1,T3BVDG2,T3BVDL2,
     *                           T3BVDG3,T3BVDL3,LAMP,ISYMD,D,
     *                           WORK,LWORK)
***********************************************************************
*                                                                     *
*     Construct virtual integrals which are required to calculate     *
*     t3bar_0 amplitudes (in terms of SMAT and UMAT).                 *
*                                                                     *
*     ISYINT - symmetry of integrals                                  *
*     ISYMD  - symmetry of external (fixed) D index)                  *
*                                                                     *
*                                                                     *
*     **********************************************************      *
*     !!! THIS ROUTINE IS TO BE USED INSIDE THE D- or B-LOOP !!!      * 
*     **********************************************************      *
*                                                                     *
*                                                                     *
*     Filip Pawlowski, 09-Dec-2002, Aarhus                            *
*                                                                     *
***********************************************************************
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INTEGER ISYINT,ISYMD,LWORK
      INTEGER LU3FOPX,LU3FOP2X,LUDKBC3,LU3VI
      INTEGER ISYCKB,IOFF
      INTEGER KINTVI,KTRVI6,KEND1,LWRK1
C
      CHARACTER*(*) FN3FOPX,FN3FOP2X,FNDKBC3,FN3VI
C
#if defined (SYS_CRAY) 
      REAL T3BVDL1(*),T3BVDG1(*),T3BVDG2(*),T3BVDL2(*),T3BVDG3(*)
      REAL T3BVDL3(*),WORK(LWORK)
      REAL LAMP(*)
#else
      DOUBLE PRECISION T3BVDL1(*),T3BVDG1(*),T3BVDG2(*),T3BVDL2(*)
      DOUBLE PRECISION T3BVDG3(*),T3BVDL3(*),WORK(LWORK)
      DOUBLE PRECISION LAMP(*)
#endif
C
      CALL QENTER('INTVIR_T3BAR0_D')
C
      ISYCKB = MULD2H(ISYINT,ISYMD)
C
C     -------------------------------------------
C     Read 2*C-E of integral used for t3-bar
C     -------------------------------------------
C
      IOFF = ICKBD(ISYCKB,ISYMD) + NCKATR(ISYCKB)*(D - 1) + 1
      IF (NCKATR(ISYCKB) .GT. 0) THEN
         CALL GETWA2(LU3FOP2X,FN3FOP2X,T3BVDL1,IOFF,NCKATR(ISYCKB))
      ENDIF
C
C     -------------------------------------------------
C     Integrals used for t3-bar for cc3
C     -------------------------------------------------
C
      IOFF = ICKBD(ISYCKB,ISYMD) + NCKATR(ISYCKB)*(D - 1) + 1
      IF (NCKATR(ISYCKB) .GT. 0) THEN
         CALL GETWA2(LUDKBC3,FNDKBC3,T3BVDG1,IOFF,NCKATR(ISYCKB))
      ENDIF
C
      CALL CCSDT_SRVIR3(T3BVDG1,WORK,ISYMD,D,ISYINT) !Probably not neccesary
C
      CALL CCSDT_SRTVIR(T3BVDG1,T3BVDG2,WORK,LWORK,ISYMD,ISYINT)
C
C     ------------------------------------------------
C     Sort the integrals for t3-bar
C     ------------------------------------------------
C
      CALL CCSDT_SRTVIR(T3BVDL1,T3BVDL2,WORK,LWORK,ISYMD,ISYINT)
C
C---------------
C     Allocation
C---------------
C
      KINTVI = 1
      KEND1 = KINTVI + NCKA(ISYCKB)
      LWRK1 = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space for allocation in '//
     *             'INTVIR_T3BAR0_D (1)')
      END IF
C
      CALL DZERO(WORK(KINTVI),NCKATR(ISYCKB))
C
C     ----------------------------------------------------
C     Read virtual integrals used in u3am for t3-bar.
C     ----------------------------------------------------
C
      IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D - 1) + 1
      IF (NCKA(ISYCKB) .GT. 0) THEN
         CALL GETWA2(LU3VI,FN3VI,WORK(KINTVI),IOFF,
     *               NCKA(ISYCKB))
      ENDIF   
C
      CALL CCSDT_TRVIR(WORK(KINTVI),T3BVDG3,LAMP,ISYMD,D,ISYINT,
     *                 WORK(KEND1),LWRK1)
C
C--------------------------------------
C     Allocation (regain the workspace)
C--------------------------------------
C
      KTRVI6 = 1
      KEND1 = KTRVI6 + NCKATR(ISYCKB)
      LWRK1 = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space for allocation in '//
     *             'INTVIR_T3BAR0_D (2)')
      END IF
C
      CALL DZERO(WORK(KTRVI6),NCKATR(ISYCKB))
C
C----------------------
C     Read in integrals 
C----------------------
C
      IOFF = ICKBD(ISYCKB,ISYMD) + NCKATR(ISYCKB)*(D - 1) + 1
      IF (NCKATR(ISYCKB) .GT. 0) THEN
         CALL GETWA2(LU3FOPX,FN3FOPX,WORK(KTRVI6),IOFF,NCKATR(ISYCKB))
      ENDIF
C
      IF (LWRK1 .LT. NCKATR(ISYCKB)) THEN
         CALL QUIT('Insufficient space for allocation in '//
     &             'INTVIR_T3BAR0_D (2)')
      END IF

      CALL DCOPY(NCKATR(ISYCKB),WORK(KTRVI6),1,T3BVDL3,1)
C
C-------------
C     End
C-------------
C
      CALL QEXIT('INTVIR_T3BAR0_D')
C
      RETURN
      END
C  /* Deck intvir_t3barx_d */
      SUBROUTINE INTVIR_T3BARX_D(LSKIPL1R,
     *                           ISYINT,LU3VI,FN3VI,LU3VI2,FN3VI2,
     *                           LU3FOP,FN3FOP,LU3FOP2,FN3FOP2,
     *                           W3BXVDGX1,W3BXVDG1,W3BXVDGX2,W3BXVDG2,
     *                           W3BXVDLX1,W3BXVDL1,W3BXVDLX2,W3BXVDL2,
     *                           LAMPY,ISYMY,LAMP0,ISYM0,ISYMD,D,
     *                           WORK,LWORK)
*
***********************************************************************
*                                                                     *
*     Construct virtual integrals which are required to calculate    *
*     t3bar_X multipliers:                                            *
*                                                                     *
*                                                                     *
*     Both LambdaH_X and LambdaH_0 transformed integrals are needed   *
*                                                                     *
*     ISYM0   - symmetry of operator in LambdaH_0 transformation      *
*     ISYINT  - symmetry of integrals before transformation
*     ISYMY   - symmetry of operator in LambdaH_Y transformation      *
*                                                                     *
*                                                                     *
*     USE ONLY INSIDE D-LOOP
*                                                                     *
*     Filip Pawlowski, 09-Dec-2002, Aarhus                            *
*                                                                     *
***********************************************************************
C
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      LOGICAL LSKIPL1R
C
      INTEGER ISYINT,ISYM0,ISYMY,ISYMD,LWORK
      INTEGER LU3VI2,LU3FOP,LU3FOP2,LU3VI
      INTEGER KINTVI,KEND1,LWRK1,IOFF
      INTEGER ISYCKB,ISYCKBY
C
      CHARACTER*(*) FN3VI2,FN3FOP,FN3FOP2,FN3VI
C
#if defined (SYS_CRAY) 
      REAL LAMP0(*),LAMPY(*),W3BXVDGX1(*),W3BXVDG1(*),W3BXVDGX2(*)
      REAL W3BXVDG2(*)
      REAL W3BXVDLX1(*),W3BXVDL1(*),W3BXVDLX2(*),W3BXVDL2(*)
      REAL WORK(LWORK)
#else
      DOUBLE PRECISION LAMP0(*),LAMPY(*),W3BXVDGX1(*),W3BXVDG1(*)
      DOUBLE PRECISION W3BXVDLX1(*),W3BXVDL1(*),W3BXVDLX2(*),W3BXVDL2(*)
      DOUBLE PRECISION W3BXVDGX2(*),W3BXVDG2(*)
      DOUBLE PRECISION WORK(LWORK)
#endif
C
      CALL QENTER('INTVIR_T3BARX_D')
C
      ISYCKB  = MULD2H(ISYINT,ISYMD)
      IF (.NOT. LSKIPL1R) THEN
         ISYCKBY = MULD2H(ISYMY,ISYMD)
      END IF
C
      KINTVI = 1
      IF (.NOT. LSKIPL1R) THEN
         KEND1  = KINTVI + MAX(NCKA(ISYCKB),NCKA(ISYCKBY))
      ELSE
         KEND1  = KINTVI + NCKA(ISYCKB)
      ENDIF
      LWRK1  = LWORK  - KEND1 
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in INTVIR_T3BARX_D')
      END IF
C
C------------------------------------------------------------------------
C
            IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D - 1) + 1
            IF (NCKA(ISYCKB) .GT. 0) THEN
               CALL GETWA2(LU3VI,FN3VI,WORK(KINTVI),IOFF,
     *                     NCKA(ISYCKB))
            ENDIF

C
C       -----------------
C       LAMPY TRANSFORMED
C       -----------------
C
C     Transform g(c1^h k1^p | delta b2^h) integrals
C     to g(c1^h k1^p | d2^pY b2^h) with lampY
C
      IF (.NOT. LSKIPL1R) THEN
        CALL INTVIR_T3BARX_D2(WORK(KINTVI),LAMPY,ISYMY,ISYINT,W3BXVDGX1,
     *                        ISYMD,D,WORK(KEND1),LWRK1)
      END IF
C
C       -----------------
C       LAMP0 TRANSFORMED
C       -----------------
C
C       Transform g(c1^h k1^p | delta b2^h) integrals
C       to g(c1^h k1^p | d2^p b2^h) with lamp0
C
      CALL INTVIR_T3BARX_D2(WORK(KINTVI),LAMP0,ISYM0,ISYINT,W3BXVDG1,
     *                      ISYMD,D,WORK(KEND1),LWRK1)
C------------------------------------------------------------------------
C
      IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D - 1) + 1
      IF (NCKA(ISYCKB) .GT. 0) THEN
C
C     Read in g(b2^h k1^p | delta c1^h) integrals and 
C     transform and sort --- see above
C
         CALL GETWA2(LU3VI2,FN3VI2,WORK(KINTVI),IOFF,NCKA(ISYCKB))
      ENDIF
C
C       -----------------
C       LAMPY TRANSFORMED
C       -----------------
C
C
      IF (.NOT. LSKIPL1R) THEN
       CALL INTVIR_T3BARX_D2(WORK(KINTVI),LAMPY,ISYMY,ISYINT,W3BXVDGX2,
     *                       ISYMD,D,WORK(KEND1),LWRK1)
      END IF
C
C       -----------------
C       LAMP0 TRANSFORMED
C       -----------------
C
      CALL INTVIR_T3BARX_D2(WORK(KINTVI),LAMP0,ISYM0,ISYINT,W3BXVDG2,
     *                      ISYMD,D,WORK(KEND1),LWRK1)
C------------------------------------------------------------------------
C
C
C           Read in L(c1^h k1^p | delta b2^h) integrals
C
               CALL GETWA2(LU3FOP,FN3FOP,WORK(KINTVI),IOFF,
     &                        NCKA(ISYCKB))
C
C                     -----------------
C                     LAMPY TRANSFORMED
C                     -----------------
C
C          Transform L(c1^h k1^p | delta b2^h) integrals
C          to L(c1^h k1^p | d2^pY b2^h) 
C
      IF (.NOT. LSKIPL1R) THEN
       CALL INTVIR_T3BARX_D2(WORK(KINTVI),LAMPY,ISYMY,ISYINT,W3BXVDLX1,
     *                       ISYMD,D,WORK(KEND1),LWRK1)
      END IF
C
C                     -----------------
C                     LAMP0 TRANSFORMED
C                     -----------------
C
C          Transform L(c1^h k1^p | delta b2^h) integrals
C          to L(c1^h k1^p | d2^p b2^h) 
C
      CALL INTVIR_T3BARX_D2(WORK(KINTVI),LAMP0,ISYM0,ISYINT,W3BXVDL1,
     *                      ISYMD,D,WORK(KEND1),LWRK1)
C------------------------------------------------------------------------
C
C     Read in L(b2^h k1^p | delta c1^h) integrals and 
C     transform and sort --- see above
C
      CALL GETWA2(LU3FOP2,FN3FOP2,WORK(KINTVI),IOFF,NCKA(ISYCKB))
C
C
C                     -----------------
C                     LAMPY TRANSFORMED
C                     -----------------
C
      IF (.NOT. LSKIPL1R) THEN
       CALL INTVIR_T3BARX_D2(WORK(KINTVI),LAMPY,ISYMY,ISYINT,W3BXVDLX2,
     *                       ISYMD,D,WORK(KEND1),LWRK1)
      END IF
C
C
C                     -----------------
C                     LAMP0 TRANSFORMED
C                     -----------------
C
      CALL INTVIR_T3BARX_D2(WORK(KINTVI),LAMP0,ISYM0,ISYINT,W3BXVDL2,
     *                      ISYMD,D,WORK(KEND1),LWRK1)
C------------------------------------------------------------------------
C
C-------------
C     End
C-------------
C
      CALL QEXIT('INTVIR_T3BARX_D')
C
      RETURN
      END
C  /* Deck intvir_t3barx_d2 */
      SUBROUTINE INTVIR_T3BARX_D2(VIRINT,LAMP,ISYMX,ISYINT,TRVIR,ISYMD,
     *                            D,WORK,LWORK)
***********************************************************************
*                                                                     *
*     Construct virtual integrals which are required to calculate    *
*     t3bar_X multipliers.                                            *
*                                                                     *
*     ISYINT  - symmetry of integrals before transformation           *
*     ISYMX   - symmetry of operator                                  *
*     ISYINTX - symmetry of integrals after  transformation           *
*                                                                     *
*     Filip Pawlowski, 09-Dec-2002, Aarhus                            *
*                                                                     *
***********************************************************************
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYMX,ISYINT,ISYMD,LWORK
      INTEGER ISYINTX
C
#if defined (SYS_CRAY) 
      REAL VIRINT(*),LAMP(*),TRVIR(*),WORK(LWORK)
#else
      DOUBLE PRECISION VIRINT(*),LAMP(*),TRVIR(*),WORK(LWORK)
#endif
C
      CALL QENTER('INTVIR_T3BARX_D2')
C
      ISYINTX = MULD2H(ISYINT,ISYMX)
C
C----------------------------------------------------------------------
C     Transform virtual g-integrals with LambdaP
C----------------------------------------------------------------------
C
      CALL CCLR_TRVIR(VIRINT,TRVIR,LAMP,ISYMX,ISYMD,D,ISYINT,WORK,LWORK)
C
C-------------------
C     Sort integrals 
C-------------------
C
      CALL CCSDT_SRVIR3(TRVIR,WORK,ISYMD,D,ISYINTX)
C
C
C-------------
C     End
C-------------
C
      CALL QEXIT('INTVIR_T3BARX_D2')
C
      RETURN
      END
C  /* Deck get_t30_bd */
      SUBROUTINE GET_T30_BD(ISYMT3,ISYINT,T2TP,ISYMT2,T3MAT,FOCKD,DIAG,
     *                      INDSQ,LENSQ,S3MAT,T3VDG1,T3VDG2,T3OG1,
     *                      IINDEX,S3MAT3,T3VBG1,T3VBG2,INDEX2,U3MAT,
     *                      T3VDG3,T3OG2,U3MAT3,T3VBG3,ISYMB,B,ISYMD,D,
     *                      ISCKIJ,WORK,LWORK)
*
***********************************************************************
*                                                                     *
*     Construct t3_0 ampllitudes for fixed B an D usin S and U        *
*     intermediates 
*                                                                     *
*     USE ONLY INSIDE BD-LOOPS
*                                                                     *
*     Filip Pawlowski, 10-Dec-2002, Aarhus                            *
*                                                                     *
***********************************************************************
C
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYMT3,ISYINT,ISYMT2,LENSQ,INDSQ(LENSQ,6),IINDEX,INDEX2
      INTEGER ISYMB,ISYMD,ISCKIJ,LWORK
      INTEGER ISYMBD
C
#if defined (SYS_CRAY) 
      REAL T2TP(*),T3MAT(*),FOCKD(*),DIAG(*)
      REAL S3MAT(*),T3VDG1(*),T3VDG2(*),T3OG1(*)
      REAL S3MAT3(*),T3VBG1(*),T3VBG2(*)
      REAL U3MAT(*),T3VDG3(*),T3OG2(*)
      REAL U3MAT3(*),T3VBG3(*)
      REAL WORK(LWORK)
#else
      DOUBLE PRECISION T2TP(*),T3MAT(*),FOCKD(*),DIAG(*)
      DOUBLE PRECISION S3MAT(*),T3VDG1(*),T3VDG2(*),T3OG1(*)
      DOUBLE PRECISION S3MAT3(*),T3VBG1(*),T3VBG2(*)
      DOUBLE PRECISION U3MAT(*),T3VDG3(*),T3OG2(*)
      DOUBLE PRECISION U3MAT3(*),T3VBG3(*)
      DOUBLE PRECISION WORK(LWORK)
#endif
C
      CALL QENTER('GET_T30_BD')
C
C------------------------------------------
C     Initial check of symmetry consistency
C------------------------------------------
C
      ISYMBD = MULD2H(ISYMB,ISYMD)
C
      IF (ISCKIJ .NE. MULD2H(ISYMBD,ISYMT3)) THEN 
         WRITE(LUPRI,*)'ISCKIJ = ', ISCKIJ
         WRITE(LUPRI,*)'MULD2H(ISYMBD,ISYMT3) = ', MULD2H(ISYMBD,ISYMT3)
         CALL QUIT('Symmetry inconsistency in GET_T30_BD')
      END IF
C
C------------------------------------------------------------------------
C     Calculate the S(ci,bk,dj) matrix for T3 for B,D.  (S3MAT)
C-------------------------------------------------------------------
C
      CALL CC3_SMAT(0.0d0,T2TP,ISYMT2,T3MAT,T3VDG1,T3VDG2,T3OG1,ISYINT,
     *              FOCKD,DIAG,S3MAT,WORK,LWORK,IINDEX,INDSQ,LENSQ,
     *              ISYMB,B,ISYMD,D)
C
      CALL T3_FORBIDDEN(S3MAT,ISYMT3,ISYMB,B,ISYMD,D)
C
C-------------------------------------------------------------------
C     Calculate the S(ci,bk,dj) matrix for T3 for D,B.  (S3MAT3)
C-------------------------------------------------------------------
C
      CALL CC3_SMAT(0.0d0,T2TP,ISYMT2,T3MAT,T3VBG1,T3VBG2,T3OG1,
     *              ISYINT,FOCKD,DIAG,S3MAT3,WORK,LWORK,INDEX2,
     *              INDSQ,LENSQ,ISYMD,D,ISYMB,B)
C
      CALL T3_FORBIDDEN(S3MAT3,ISYMT3,ISYMD,D,ISYMB,B)
C
C---------------------------------------------------------------------------
C                 Calculate U(ci,jk) for fixed b,d.               (U3MAT)
C--------------------------------------------------
C
      CALL CC3_UMAT(0.0D0,T2TP,ISYMT2,T3VDG3,T3OG2,ISYINT,FOCKD,DIAG,
     *              U3MAT,T3MAT,WORK,LWORK,INDSQ,LENSQ,ISYMB,B,ISYMD,D)
C
      CALL T3_FORBIDDEN(U3MAT,ISYMT3,ISYMB,B,ISYMD,D)
C
C--------------------------------------------------
C     Calculate U(ci,jk) for fixed d,b.               (U3MAT3)
C--------------------------------------------------
C
      CALL CC3_UMAT(0.0D0,T2TP,ISYMT2,T3VBG3,T3OG2,ISYINT,FOCKD,DIAG,
     *              U3MAT3,T3MAT,WORK,LWORK,INDSQ,LENSQ,ISYMD,D,ISYMB,B)
C
      CALL T3_FORBIDDEN(U3MAT3,ISYMT3,ISYMD,D,ISYMB,B)
C
C--------------------------------------------------
C     Sum up S and U intermediates to get T3 BD amplitudes
C--------------------------------------------------
C
      CALL CC3_T3BD(ISCKIJ,S3MAT,S3MAT3,U3MAT,U3MAT3,T3MAT,INDSQ,LENSQ)
C
C-------------
C     End
C-------------
C
      CALL QEXIT('GET_T30_BD')
C
      RETURN
      END
C  /* Deck get_t3bar0_bd */
      SUBROUTINE GET_T3BAR0_BD(ISYMT3,T1AM,ISYMT1,T2TP,ISYMT2,TMAT,
     *                         FCKBA,FOCKD,DIAG,XIAJB,ISINT1,ISINT2,
     *                         INDSQ,LENSQ,SMAT2,T3BVDG1,T3BVDG2,
     *                         T3BVDL1,T3BVDL2,T3BOG1,T3BOL1,IINDEX,
     *                         SMAT4,T3BVBG1,T3BVBG2,T3BVBL1,T3BVBL2,
     *                         INDEX2,UMAT2,T3BVDG3,T3BVDL3,T3BOG2,
     *                         T3BOL2,UMAT4,T3BVBG3,T3BVBL3,ISYMB,B,
     *                         ISYMD,D,ISCKIJ,WORK,LWORK)
*
***********************************************************************
*                                                                     *
*     Construct t3bar_0 multipliers for fixed B an D usin S and U        *
*     intermediates 
*                                                                     *
*     USE ONLY INSIDE BD-LOOPS
*                                                                     *
*     Filip Pawlowski, 10-Dec-2002, Aarhus                            *
*                                                                     *
***********************************************************************
C
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYMT3,ISYMT1,ISYMT2,ISINT1,ISINT2
      INTEGER LENSQ,INDSQ(LENSQ,6),IINDEX,INDEX2
      INTEGER ISYMB,ISYMD,ISCKIJ,LWORK
      INTEGER ISYMBD
C
#if defined (SYS_CRAY) 
      REAL T1AM(*),T2TP(*),TMAT(*)
      REAL FCKBA(*),FOCKD(*),DIAG(*),XIAJB(*)
      REAL SMAT2(*),T3BVDG1(*),T3BVDG2(*),T3BVDL1(*),T3BVDL2(*)
      REAL T3BOG1(*),T3BOL1(*)
      REAL SMAT4(*),T3BVBG1(*),T3BVBG2(*),T3BVBL1(*),T3BVBL2(*)
      REAL UMAT2(*),T3BVDG3(*),T3BVDL3(*),T3BOG2(*),T3BOL2(*)
      REAL UMAT4(*),T3BVBG3(*),T3BVBL3(*)
      REAL WORK(LWORK)
      REAL HALF
#else
      DOUBLE PRECISION  T1AM(*),T2TP(*),TMAT(*)
      DOUBLE PRECISION  FCKBA(*),FOCKD(*),DIAG(*),XIAJB(*)
      DOUBLE PRECISION  SMAT2(*),T3BVDG1(*),T3BVDG2(*),T3BVDL1(*)
      DOUBLE PRECISION  T3BVDL2(*),T3BOG1(*),T3BOL1(*)
      DOUBLE PRECISION  SMAT4(*),T3BVBG1(*),T3BVBG2(*),T3BVBL1(*)
      DOUBLE PRECISION  T3BVBL2(*)
      DOUBLE PRECISION  UMAT2(*),T3BVDG3(*),T3BVDL3(*),T3BOG2(*)
      DOUBLE PRECISION  T3BOL2(*)
      DOUBLE PRECISION  UMAT4(*),T3BVBG3(*),T3BVBL3(*)
      DOUBLE PRECISION  WORK(LWORK)
      DOUBLE PRECISION  HALF
#endif
C
      PARAMETER(HALF = 0.5D0)
C
      CALL QENTER('GET_T3BAR0_BD')
C
C------------------------------------------
C     Initial check of symmetry consistency
C------------------------------------------
C
      ISYMBD = MULD2H(ISYMB,ISYMD)
C
      IF (ISCKIJ .NE. MULD2H(ISYMBD,ISYMT3)) THEN
         WRITE(LUPRI,*)'ISCKIJ = ', ISCKIJ
         WRITE(LUPRI,*)'MULD2H(ISYMBD,ISYMT3) = ', MULD2H(ISYMBD,ISYMT3)
         CALL QUIT('Symmetry inconsistency in GET_T3BAR0_BD')
      END IF
C     ----------------------------------------------------
C     Calculate the S(ci,bk,dj) matrix for B,D for T3-BAR:
C     ----------------------------------------------------
      CALL DZERO(SMAT2,NCKIJ(ISCKIJ))
C
      CALL CCFOP_SMAT(0.0D0,T1AM,ISYMT1,T2TP,ISYMT2,TMAT,FCKBA,XIAJB,
     *                ISINT1,T3BVDG1,T3BVDG2,T3BVDL1,T3BVDL2,T3BOG1,
     *                T3BOL1,ISINT2,FOCKD,DIAG,SMAT2,WORK,LWORK,IINDEX,
     *                INDSQ,LENSQ,ISYMB,B,ISYMD,D)
C
      CALL DSCAL(NCKIJ(ISCKIJ),HALF,SMAT2,1)
C
      CALL T3_FORBIDDEN(SMAT2,ISYMT3,ISYMB,B,ISYMD,D)
C
C     ----------------------------------------------------
C     Calculate the S(ci,bk,dj) matrix for D,B for T3-BAR:
C     ----------------------------------------------------
      CALL DZERO(SMAT4,NCKIJ(ISCKIJ))

      CALL CCFOP_SMAT(0.0D0,T1AM,ISYMT1,T2TP,ISYMT2,TMAT,FCKBA,XIAJB,
     *                ISINT1,T3BVBG1,T3BVBG2,T3BVBL1,T3BVBL2,T3BOG1,
     *                T3BOL1,ISINT2,FOCKD,DIAG,SMAT4,WORK,LWORK,INDEX2,
     *                INDSQ,LENSQ,ISYMD,D,ISYMB,B)
C
      CALL DSCAL(NCKIJ(ISCKIJ),HALF,SMAT4,1)
C
      CALL T3_FORBIDDEN(SMAT4,ISYMT3,ISYMD,D,ISYMB,B)
C
C     ------------------------------------------------
C     Calculate U(ci,jk) for fixed b,d for t3-bar.
C     ------------------------------------------------
      CALL DZERO(UMAT2,NCKIJ(ISCKIJ))

      CALL CCFOP_UMAT(0.0D0,T1AM,ISYMT1,T2TP,ISYMT2,XIAJB,ISINT1,
     *                FCKBA,T3BVDG3,T3BVDL3,T3BOG2,T3BOL2,ISINT2,FOCKD,
     *                DIAG,UMAT2,TMAT,WORK,LWORK,INDSQ,LENSQ,ISYMB,B,
     *                ISYMD,D)
C
      CALL DSCAL(NCKIJ(ISCKIJ),HALF,UMAT2,1)
C
      CALL T3_FORBIDDEN(UMAT2,ISYMT3,ISYMB,B,ISYMD,D)
C
C     ------------------------------------------------
C     Calculate U(ci,jk) for fixed d,b for t3-bar.
C     ------------------------------------------------
      CALL DZERO(UMAT4,NCKIJ(ISCKIJ))
C
      CALL CCFOP_UMAT(0.0D0,T1AM,ISYMT1,T2TP,ISYMT2,XIAJB,ISINT1,
     *                FCKBA,T3BVBG3,T3BVBL3,T3BOG2,T3BOL2,ISINT2,FOCKD,
     *                DIAG,UMAT4,TMAT,WORK,LWORK,INDSQ,LENSQ,ISYMD,D,
     *                ISYMB,B)
C
      CALL DSCAL(NCKIJ(ISCKIJ),HALF,UMAT4,1)

      CALL T3_FORBIDDEN(UMAT4,ISYMT3,ISYMD,D,ISYMB,B)
C
C     -------------------------------------------
C     Sum up the S-bar and U-bar to get a real T3-bar
C     -------------------------------------------
      CALL CC3_T3BD(ISCKIJ,SMAT2,SMAT4,UMAT2,UMAT4,TMAT,INDSQ,LENSQ)
C
C-------------
C     End
C-------------
C
      CALL QEXIT('GET_T3BAR0_BD')
C
      RETURN
      END
C  /* Deck get_t3barx_bd */
      SUBROUTINE GET_T3BARX_BD(NOVIRT,
     *                         TMAT,ISTMAT,FOCKY,ISYFKY,WMAT,ISWMAT,
     *                         L2TP,ISYML2,FCKYCK,ISYFCKY,W3BXVDLX2,
     *                         W3BXVDLX1,W3BXVDGX2,W3BXVDGX1,W3BXOLX1,
     *                         W3BXOGX1,ISINT2Y,INDAJLB,INDAJLC,INDSQ,
     *                         LENSQ,L2TPY,ISYML2Y,FCKBA,ISYFCKBA,
     *                         W3BXVDL2,W3BXVDL1,W3BXVDG2,W3BXVDG1,
     *                         W3BXOL1,W3BXOG1,ISINT2,INDAJLBY,INDAJLCY,
     *                         L1Y,ISYML1Y,XIAJB,ISINT1,FREQ,DIAGW,
     *                         FOCKD,B,ISYMB,D,ISYMD,ISYMT3,WORK,LWORK)
*
***********************************************************************
*                                                                     *
*     Construct t3bar_X multipliers for fixed B an D using W          *
*     intermediates 
*                                                                     *
*     USE ONLY INSIDE BD-LOOPS
*                                                                     *
*     Filip Pawlowski, 10-Dec-2002, Aarhus                            *
*                                                                     *
***********************************************************************
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      LOGICAL NOVIRT
      INTEGER ISTMAT,ISYFKY,ISWMAT,ISYML2,ISYFCKY,ISINT2Y,ISYML2Y
      INTEGER ISYFCKBA,ISINT2,ISYML1Y,ISINT1,ISYMB,ISYMD,ISYMT3
      INTEGER INDAJLB,INDAJLC,LENSQ,INDSQ(LENSQ,6),INDAJLBY,INDAJLCY
      INTEGER LWORK
      INTEGER ISYMBD
C
#if defined (SYS_CRAY) 
      REAL TMAT(*),FOCKY(*),WMAT(*),L2TP(*),FCKYCK(*)
      REAL W3BXVDLX2(*),W3BXVDLX1(*),W3BXVDGX2(*),W3BXVDGX1(*)
      REAL W3BXOLX1(*),W3BXOGX1(*),L2TPY(*),FCKBA(*)
      REAL W3BXVDL2(*),W3BXVDL1(*),W3BXVDG2(*),W3BXVDG1(*)
      REAL W3BXOL1(*),W3BXOG1(*),L1Y(*),XIAJB(*)
      REAL DIAGW(*),FOCKD(*),WORK(LWORK)
      REAL FREQ
#else
      DOUBLE PRECISION TMAT(*),FOCKY(*),WMAT(*),L2TP(*),FCKYCK(*)
      DOUBLE PRECISION W3BXVDLX2(*),W3BXVDLX1(*)
      DOUBLE PRECISION W3BXVDGX2(*),W3BXVDGX1(*)
      DOUBLE PRECISION W3BXOLX1(*),W3BXOGX1(*),L2TPY(*),FCKBA(*)
      DOUBLE PRECISION W3BXVDL2(*),W3BXVDL1(*),W3BXVDG2(*),W3BXVDG1(*)
      DOUBLE PRECISION W3BXOL1(*),W3BXOG1(*),L1Y(*),XIAJB(*)
      DOUBLE PRECISION DIAGW(*),FOCKD(*),WORK(LWORK)
      DOUBLE PRECISION FREQ
#endif
C
      CALL QENTER('GET_T3BARX_BD')
C
C------------------------------------------
C     Initial check of symmetry consistency
C------------------------------------------
C
      ISYMBD = MULD2H(ISYMB,ISYMD)
C
      IF (ISWMAT .NE. MULD2H(ISYMBD,ISYMT3)) THEN
         WRITE(LUPRI,*)'ISWMAT = ', ISWMAT
         WRITE(LUPRI,*)'MULD2H(ISYMBD,ISYMT3) = ', MULD2H(ISYMBD,ISYMT3)
         CALL QUIT('Symmetry inconsistency in GET_T3BARX_BD')
      END IF
C
C    calculate     <L3|[Y^,tau3]|HF>
C
      IF (.NOT.NOVIRT) THEN
         !virtual part
         CALL WBARBD_V(TMAT,ISTMAT,FOCKY,ISYFKY,WMAT,ISWMAT,WORK,LWORK)
      END IF
C
      !occupied part
      CALL WBARBD_O(TMAT,ISTMAT,FOCKY,ISYFKY,WMAT,ISWMAT,WORK,LWORK)
C 
C    calculate     <L2|[Y,tau3]|HF>
C
      CALL WBARBD_T2(B,ISYMB,D,ISYMD,L2TP,ISYML2,FOCKY,ISYFKY,WMAT,
     *               ISWMAT)
C
C    calculate     <L2|[H^Y,tau3]|HF>
C
      CALL WBARBD_TMAT(L2TP,ISYML2,WMAT,TMAT,ISWMAT,FCKYCK,ISYFCKY,
     *                 W3BXVDLX2,W3BXVDLX1,W3BXVDGX2,W3BXVDGX1,W3BXOLX1,
     *                 W3BXOGX1,ISINT2Y,WORK,LWORK,INDAJLB,INDAJLC,
     *                 INDSQ,LENSQ,ISYMB,B,ISYMD,D)
C
C    calculate     <L2Y|[H^,tau3]|HF>
C              
      CALL WBARBD_TMAT(L2TPY,ISYML2Y,WMAT,TMAT,ISWMAT,FCKBA,ISYFCKBA,
     *                 W3BXVDL2,W3BXVDL1,W3BXVDG2,W3BXVDG1,W3BXOL1,
     *                 W3BXOG1,ISINT2,WORK,LWORK,INDAJLBY,INDAJLCY,
     *                 INDSQ,LENSQ,ISYMB,B,ISYMD,D)
C
C    calculate     <L1Y|[H^,tau3]|HF>
C
      CALL WBARBD_L1(L1Y,ISYML1Y,TMAT,XIAJB,ISINT1,WMAT,WORK,LWORK,
     *               INDSQ,LENSQ,ISYMB,B,ISYMD,D)
C
C------------------------------------------------
C     Divide by the energy difference and
C     remove the forbidden elements
C------------------------------------------------
C
      CALL WBD_DIA(B,ISYMB,D,ISYMD,FREQ,ISWMAT,WMAT,DIAGW,FOCKD)
C
      CALL T3_FORBIDDEN(WMAT,ISYMT3,ISYMB,B,ISYMD,D)
C
C
C-------------
C     End
C-------------
C
      CALL QEXIT('GET_T3BARX_BD')
C
      RETURN
      END
C  /* Deck get_t3x_bd */
      SUBROUTINE GET_T3X_BD(ISYMT3,WMAT,ISWMAT,TMAT,ISTMAT,FOCKX,ISYFKX,
     *                      T2TP,ISYMT2,T2TPX,ISYMT2X,
     *                      W3XVDG1,W3XVDG2,W3XOG1,ISYINT,
     *                      W3XVDGX1,W3XVDGX2,W3XOGX1,ISYINTX,
     *                      FREQX,DIAG,FOCKD,
     *                      INDSQ,LENSQ,B,ISYMB,D,ISYMD,WORK,LWORK)
*
***********************************************************************
*                                                                     *
*     Construct t3_X amplitudes for fixed B an D using W              *
*     intermediates 
*                                                                     *
*     USE ONLY INSIDE BD-LOOPS
*                                                                     *
*     Filip Pawlowski, 13-Dec-2002, Aarhus                            *
*                                                                     *
***********************************************************************
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYMT3,ISWMAT,ISTMAT,ISYFKX,ISYMT2,ISYMT2X,ISYINT
      INTEGER ISYINTX,LENSQ,INDSQ(LENSQ,6),ISYMB,ISYMD,LWORK
      INTEGER ISYMBD
C
#if defined (SYS_CRAY) 
      REAL WMAT(*),TMAT(*),FOCKX(*),T2TP(*),T2TPX(*),DIAG(*),FOCKD(*)
      REAL W3XVDG1(*),W3XVDG2(*),W3XOG1(*)
      REAL W3XVDGX1(*),W3XVDGX2(*),W3XOGX1(*)
      REAL WORK(LWORK)
      REAL FREQX
#else
      DOUBLE PRECISION WMAT(*),TMAT(*),FOCKX(*),T2TP(*),T2TPX(*)
      DOUBLE PRECISION DIAG(*),FOCKD(*)
      DOUBLE PRECISION W3XVDG1(*),W3XVDG2(*),W3XOG1(*)
      DOUBLE PRECISION W3XVDGX1(*),W3XVDGX2(*),W3XOGX1(*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION FREQX
#endif
C
      CALL QENTER('GET_T3X_BD')
C
C------------------------------------------
C     Initital test of symmetry consistency
C------------------------------------------
C
      ISYMBD = MULD2H(ISYMB,ISYMD)
C
      IF (MULD2H(ISYMT3,ISYMBD) .NE. ISWMAT) THEN
         WRITE(LUPRI,*)'ISYMT3 = ', ISYMT3
         WRITE(LUPRI,*)'ISWMAT = ', ISWMAT
         CALL QUIT('Symmetry inconsistency in GET_T3X_BD')
      END IF
C
C------------------------------------------------------
C     Calculate the  term <mu3|[X,T3]|HF> virtual contribution 
C     added in W^BD(KWMAT)
C------------------------------------------------------
C
      CALL WBD_V(TMAT,ISTMAT,FOCKX,ISYFKX,WMAT,ISWMAT,WORK,LWORK)
C
C------------------------------------------------------
C     Calculate the  term <mu3|[X,T3]|HF> occupied contribution 
C     added in W^BD(KWMAT)
C------------------------------------------------------
C
      CALL WBD_O(TMAT,ISTMAT,FOCKX,ISYFKX,WMAT,ISWMAT,WORK,LWORK)
C
C------------------------------------------------------
C     Calculate the  term <mu3|[[X,T2],T2]|HF> 
C     added in W^BD(KWMAT)
C------------------------------------------------------
C
      CALL WBD_T2(.FALSE.,B,ISYMB,D,ISYMD,T2TP,ISYMT2,T2TP,ISYMT2,
     *                 FOCKX,ISYFKX,
     *                 INDSQ,LENSQ,WMAT,ISWMAT,WORK,LWORK)
C
C------------------------------------------------------
C     To get the entire T3^X add the two terms
C------------------------------------------------------
C
      CALL WBD_GROUND(T2TPX,ISYMT2X,TMAT,W3XVDG1,W3XVDG2,W3XOG1,ISYINT,
     *                WMAT,WORK,LWORK,INDSQ,LENSQ,ISYMB,B,ISYMD,D)
C
      CALL WBD_GROUND(T2TP,ISYMT2,TMAT,W3XVDGX1,W3XVDGX2,W3XOGX1,
     *                ISYINTX,WMAT,WORK,LWORK,INDSQ,LENSQ,ISYMB,B,
     *                ISYMD,D)
C
C------------------------------------------------
C     Divide by the energy difference and
C     remove the forbidden elements
C------------------------------------------------
C
      CALL WBD_DIA(B,ISYMB,D,ISYMD,FREQX,ISWMAT,WMAT,DIAG,FOCKD) 
C
      CALL T3_FORBIDDEN(WMAT,ISYMT3,ISYMB,B,ISYMD,D) 
C
C-------------
C     End
C-------------
C
      CALL QEXIT('GET_T3X_BD')
C
      RETURN
      END
C  /* Deck get_e3_bd */
      SUBROUTINE GET_E3_BD(ISYMT3,WMAT,ISWMAT,TMAT,
     *                      T2TP,ISYMT2,T2TPX,ISYMT2X,
     *                      W3XVDG1,W3XVDG2,W3XOG1,ISYINT,
     *                      W3XVDGX1,W3XVDGX2,W3XOGX1,ISYINTX,
     *                      FREQX,DIAG,FOCKD,
     *                      INDSQ,LENSQ,B,ISYMB,D,ISYMD,WORK,LWORK)
*
***********************************************************************
*                                                                     *
*     Construct e3 eigenvectors for fixed B an D using W              *
*     intermediates 
*                                                                     *
*     USE ONLY INSIDE BD-LOOPS
*                                                                     *
*     Filip Pawlowski, 12-Nov-2003, Aarhus                            *
*                                                                     *
***********************************************************************
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYMT3,ISWMAT,ISYMT2,ISYMT2X,ISYINT
      INTEGER ISYINTX,LENSQ,INDSQ(LENSQ,6),ISYMB,ISYMD,LWORK
      INTEGER ISYMBD
C
#if defined (SYS_CRAY) 
      REAL WMAT(*),TMAT(*),T2TP(*),T2TPX(*),DIAG(*),FOCKD(*)
      REAL W3XVDG1(*),W3XVDG2(*),W3XOG1(*)
      REAL W3XVDGX1(*),W3XVDGX2(*),W3XOGX1(*)
      REAL WORK(LWORK)
      REAL FREQX
#else
      DOUBLE PRECISION WMAT(*),TMAT(*),T2TP(*),T2TPX(*)
      DOUBLE PRECISION DIAG(*),FOCKD(*)
      DOUBLE PRECISION W3XVDG1(*),W3XVDG2(*),W3XOG1(*)
      DOUBLE PRECISION W3XVDGX1(*),W3XVDGX2(*),W3XOGX1(*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION FREQX
#endif
C
      CALL QENTER('E3BD')
C
C------------------------------------------
C     Initital test of symmetry consistency
C------------------------------------------
C
      ISYMBD = MULD2H(ISYMB,ISYMD)
C
      IF (MULD2H(ISYMT3,ISYMBD) .NE. ISWMAT) THEN
         WRITE(LUPRI,*)'ISYMT3 = ', ISYMT3
         WRITE(LUPRI,*)'ISWMAT = ', ISWMAT
         CALL QUIT('Symmetry inconsistency in GET_E3_BD')
      END IF
C
C------------------------------------------------------
C     <mu3|[H,E2]|HF> + <mu3|[[H,E1],T20]|HF> 
C------------------------------------------------------
C
      CALL WBD_GROUND(T2TPX,ISYMT2X,TMAT,W3XVDG1,W3XVDG2,W3XOG1,ISYINT,
     *                WMAT,WORK,LWORK,INDSQ,LENSQ,ISYMB,B,ISYMD,D)
C
      CALL WBD_GROUND(T2TP,ISYMT2,TMAT,W3XVDGX1,W3XVDGX2,W3XOGX1,
     *                ISYINTX,WMAT,WORK,LWORK,INDSQ,LENSQ,ISYMB,B,
     *                ISYMD,D)
C
C------------------------------------------------
C     Divide by the energy difference and
C     remove the forbidden elements
C------------------------------------------------
C
      CALL WBD_DIA(B,ISYMB,D,ISYMD,FREQX,ISWMAT,WMAT,DIAG,FOCKD) 
C
      CALL T3_FORBIDDEN(WMAT,ISYMT3,ISYMB,B,ISYMD,D) 
C
C-------------
C     End
C-------------
C
      CALL QEXIT('E3BD')
C
      RETURN
      END
C  /* Deck get_m3bar_bd */
      SUBROUTINE GET_M3BAR_BD(TMAT,WMAT,ISWMAT,L2TP,ISYML2,FCKYCK,
     *                        ISYFCKY,W3BXVDLX2,W3BXVDLX1,W3BXVDGX2,
     *                        W3BXVDGX1,W3BXOLX1,W3BXOGX1,ISINT2Y,
     *                        INDAJLB,INDAJLC,INDSQ,LENSQ,M2TP,ISYMM2,
     *                        FCKBA,ISYFCKBA,W3BXVDL2,W3BXVDL1,W3BXVDG2,
     *                        W3BXVDG1,W3BXOL1,W3BXOG1,ISINT2,INDAJLBY,
     *                        INDAJLCY,M1,ISYMM1,XIAJB,ISINT1,FREQ,
     *                        DIAGW,FOCKD,B,ISYMB,D,ISYMD,ISYMT3,WORK,
     *                        LWORK)
*
***********************************************************************
*                                                                     *
*     Construct the triples part of the Mbar auxilary vector (used in *
*     the calculation of the transition moments) fixed B an D using W *
*     intermediates                                                   *
*                                                                     *
*     USE ONLY INSIDE BD-LOOPS                                        *
*                                                                     *
*     Filip Pawlowski, 07-Jan-2002, Aarhus                            *
*                                                                     *
***********************************************************************
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISWMAT,ISYML2,ISYFCKY,ISINT2Y,ISYMM2
      INTEGER ISYFCKBA,ISINT2,ISYMM1,ISINT1,ISYMB,ISYMD,ISYMT3
      INTEGER INDAJLB,INDAJLC,LENSQ,INDSQ(LENSQ,6),INDAJLBY,INDAJLCY
      INTEGER LWORK
      INTEGER ISYMBD
C
#if defined (SYS_CRAY) 
      REAL TMAT(*),WMAT(*),L2TP(*),FCKYCK(*)
      REAL W3BXVDLX2(*),W3BXVDLX1(*),W3BXVDGX2(*),W3BXVDGX1(*)
      REAL W3BXOLX1(*),W3BXOGX1(*),M2TP(*),FCKBA(*)
      REAL W3BXVDL2(*),W3BXVDL1(*),W3BXVDG2(*),W3BXVDG1(*)
      REAL W3BXOL1(*),W3BXOG1(*),M1(*),XIAJB(*)
      REAL DIAGW(*),FOCKD(*),WORK(LWORK)
      REAL FREQ
#else
      DOUBLE PRECISION TMAT(*),WMAT(*),L2TP(*),FCKYCK(*)
      DOUBLE PRECISION W3BXVDLX2(*),W3BXVDLX1(*)
      DOUBLE PRECISION W3BXVDGX2(*),W3BXVDGX1(*)
      DOUBLE PRECISION W3BXOLX1(*),W3BXOGX1(*),M2TP(*),FCKBA(*)
      DOUBLE PRECISION W3BXVDL2(*),W3BXVDL1(*),W3BXVDG2(*),W3BXVDG1(*)
      DOUBLE PRECISION W3BXOL1(*),W3BXOG1(*),M1(*),XIAJB(*)
      DOUBLE PRECISION DIAGW(*),FOCKD(*),WORK(LWORK)
      DOUBLE PRECISION FREQ
#endif
C
      CALL QENTER('GET_M3BAR_BD')
C
C------------------------------------------
C     Initial check of symmetry consistency
C------------------------------------------
C
      ISYMBD = MULD2H(ISYMB,ISYMD)
C
      IF (ISWMAT .NE. MULD2H(ISYMBD,ISYMT3)) THEN
         WRITE(LUPRI,*)'ISWMAT = ', ISWMAT
         WRITE(LUPRI,*)'MULD2H(ISYMBD,ISYMT3) = ', MULD2H(ISYMBD,ISYMT3)
         CALL QUIT('Symmetry inconsistency in GET_M3BAR_BD')
      END IF
C
C    calculate     <L2|[[H,R1],tau3]|HF>
C
      CALL WBARBD_TMAT(L2TP,ISYML2,WMAT,TMAT,ISWMAT,FCKYCK,ISYFCKY,
     *                 W3BXVDLX2,W3BXVDLX1,W3BXVDGX2,W3BXVDGX1,W3BXOLX1,
     *                 W3BXOGX1,ISINT2Y,WORK,LWORK,INDAJLB,INDAJLC,
     *                 INDSQ,LENSQ,ISYMB,B,ISYMD,D)
C
C    calculate     <M2|[H^,tau3]|HF>
C              
      CALL WBARBD_TMAT(M2TP,ISYMM2,WMAT,TMAT,ISWMAT,FCKBA,ISYFCKBA,
     *                 W3BXVDL2,W3BXVDL1,W3BXVDG2,W3BXVDG1,W3BXOL1,
     *                 W3BXOG1,ISINT2,WORK,LWORK,INDAJLBY,INDAJLCY,
     *                 INDSQ,LENSQ,ISYMB,B,ISYMD,D)
C
C    calculate     <M1|[H^,tau3]|HF>
C
      CALL WBARBD_L1(M1,ISYMM1,TMAT,XIAJB,ISINT1,WMAT,WORK,LWORK,
     *               INDSQ,LENSQ,ISYMB,B,ISYMD,D)

C
C------------------------------------------------
C     Divide by the energy difference and
C     remove the forbidden elements
C------------------------------------------------
C
      CALL WBD_DIA(B,ISYMB,D,ISYMD,FREQ,ISWMAT,WMAT,DIAGW,FOCKD)
C
      CALL T3_FORBIDDEN(WMAT,ISYMT3,ISYMB,B,ISYMD,D)
C
C-------------
C     End
C-------------
C
      CALL QEXIT('GET_M3BAR_BD')
C
      RETURN
      END
C  /* Deck get_l3bar_bd */
      SUBROUTINE GET_L3BAR_BD(TMAT,WMAT,ISWMAT,INDSQ,LENSQ,L2TP,ISYML2,
     *                        FCKBA,ISYFCKBA,W3BXVDL2,W3BXVDL1,W3BXVDG2,
     *                        W3BXVDG1,W3BXOL1,W3BXOG1,ISINT2,INDAJLBY,
     *                        INDAJLCY,L1,ISYML1,XIAJB,ISINT1,FREQ,
     *                        DIAGW,FOCKD,B,ISYMB,D,ISYMD,ISYMT3,WORK,
     *                        LWORK)
*
***********************************************************************
*                                                                     *
*     Construct the triples part of the left eigenvector of the       *
*     Jacobian for fixed B an D using W intermediates                 *
*                                                                     *
*     USE ONLY INSIDE BD-LOOPS                                        *
*                                                                     *
*     Filip Pawlowski, 07-Jan-2002, Aarhus                            *
*                                                                     *
***********************************************************************
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISWMAT,ISYML2
      INTEGER ISYFCKBA,ISINT2,ISYML1,ISINT1,ISYMB,ISYMD,ISYMT3
      INTEGER LENSQ,INDSQ(LENSQ,6),INDAJLBY,INDAJLCY
      INTEGER LWORK
      INTEGER ISYMBD
C
#if defined (SYS_CRAY) 
      REAL TMAT(*),WMAT(*)
      REAL L2TP(*),FCKBA(*)
      REAL W3BXVDL2(*),W3BXVDL1(*),W3BXVDG2(*),W3BXVDG1(*)
      REAL W3BXOL1(*),W3BXOG1(*),L1(*),XIAJB(*)
      REAL DIAGW(*),FOCKD(*),WORK(LWORK)
      REAL FREQ
#else
      DOUBLE PRECISION TMAT(*),WMAT(*)
      DOUBLE PRECISION L2TP(*),FCKBA(*)
      DOUBLE PRECISION W3BXVDL2(*),W3BXVDL1(*),W3BXVDG2(*),W3BXVDG1(*)
      DOUBLE PRECISION W3BXOL1(*),W3BXOG1(*),L1(*),XIAJB(*)
      DOUBLE PRECISION DIAGW(*),FOCKD(*),WORK(LWORK)
      DOUBLE PRECISION FREQ
#endif
C
      CALL QENTER('GET_L3BAR_BD')
C
C------------------------------------------
C     Initial check of symmetry consistency
C------------------------------------------
C
      ISYMBD = MULD2H(ISYMB,ISYMD)
C
      IF (ISWMAT .NE. MULD2H(ISYMBD,ISYMT3)) THEN
         WRITE(LUPRI,*)'ISWMAT = ', ISWMAT
         WRITE(LUPRI,*)'MULD2H(ISYMBD,ISYMT3) = ', MULD2H(ISYMBD,ISYMT3)
         CALL QUIT('Symmetry inconsistency in GET_L3BAR_BD')
      END IF
C
C    calculate     <L2|[H^,tau3]|HF>
C              
      CALL WBARBD_TMAT(L2TP,ISYML2,WMAT,TMAT,ISWMAT,FCKBA,ISYFCKBA,
     *                 W3BXVDL2,W3BXVDL1,W3BXVDG2,W3BXVDG1,W3BXOL1,
     *                 W3BXOG1,ISINT2,WORK,LWORK,INDAJLBY,INDAJLCY,
     *                 INDSQ,LENSQ,ISYMB,B,ISYMD,D)
C
C    calculate     <L1|[H^,tau3]|HF>
C
      CALL WBARBD_L1(L1,ISYML1,TMAT,XIAJB,ISINT1,WMAT,WORK,LWORK,
     *               INDSQ,LENSQ,ISYMB,B,ISYMD,D)
C
C------------------------------------------------
C     Divide by the energy difference and
C     remove the forbidden elements
C------------------------------------------------
C
      CALL WBD_DIA(B,ISYMB,D,ISYMD,FREQ,ISWMAT,WMAT,DIAGW,FOCKD)
C
      CALL T3_FORBIDDEN(WMAT,ISYMT3,ISYMB,B,ISYMD,D)
C
C-------------
C     End
C-------------
C
      CALL QEXIT('GET_L3BAR_BD')
C
      RETURN
      END
C  /* Deck write_t3_dl */
      SUBROUTINE WRITE_T3_DL(LUFILE,FNFILE,T3,ISYMT3,ISYMD,ISYMB,B)
***************************************************
*  
*  Write T3 amplitudes as T3^D(ai,bj,l) to disc.
*  
*  F. Pawlowski, 25-02-2003, Aarhus.
***************************************************
C
      IMPLICIT NONE
C
#include <ccsdsym.h>
#include <ccorb.h>
#include <cc3t3d.h>
C
      INTEGER LUFILE
      INTEGER ISYMT3,ISYMD,ISYMB
      INTEGER ISYML,ISYMDL,ISAIBJ,ISYMJ,ISYMBJ,ISYMAI,ISYAIL
      INTEGER KOFF1,NBJ,IADR
C
      CHARACTER*(*) FNFILE
C
#if defined (SYS_CRAY)
      REAL T3(*)
#else
      DOUBLE PRECISION T3(*)
#endif
C
      CALL QENTER('WRITE_T3_DL')
C
      DO ISYML = 1, NSYM
         ISYMDL = MULD2H(ISYMD,ISYML)
         ISAIBJ = MULD2H(ISYMT3,ISYMDL)
         DO L = 1, NRHF(ISYML)
            DO ISYMJ = 1, NSYM
               ISYMBJ = MULD2H(ISYMJ,ISYMB)
               ISYMAI = MULD2H(ISAIBJ,ISYMBJ)
               ISYAIL = MULD2H(ISYMAI,ISYML)
               DO J = 1, NRHF(ISYMJ)

                 KOFF1 = ISAIKJ(ISYAIL,ISYMJ)+ NCKI(ISYAIL)*(J-1)
     *                 + ISAIK(ISYMAI, ISYML)+NT1AM(ISYMAI)*(L-1)
     *                 + 1

                 NBJ  = IT1AM(ISYMB,ISYMJ)+NVIR(ISYMB)*(J-1)+B
                 IADR = ISWTL(ISAIBJ,ISYML)+NT2SQ(ISAIBJ)*(L-1)
     *                + IT2SQ(ISYMAI,ISYMBJ)+NT1AM(ISYMAI)*(NBJ-1)
     *                + 1

                 CALL PUTWA2(LUFILE,FNFILE,T3(KOFF1),IADR,NT1AM(ISYMAI))
C
               END DO
            END DO
         END DO
      END DO
C
C--------------
C     End.
C--------------
C
      CALL QEXIT('WRITE_T3_DL')
C
      RETURN
      END
C
C  /* Deck sort_t2_ij */
      SUBROUTINE SORT_T2_IJ(XIJ,ISYMA,A,ISYMB,B,T2TP,ISYMT2)
C
C   XIJ^AB(ij) = t2tp(aijb)
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYMA, ISYMB, ISYMT2
      INTEGER ISYMAB, ISYMIJ, ISYMI, ISYMJ, ISYMAIJ, ISYMAI
      INTEGER KOFF1, KOFF2
#if defined (SYS_CRAY)
      REAL XIJ(*), T2TP(*)
#else
      DOUBLE PRECISION XIJ(*), T2TP(*)
#endif
C
      CALL QENTER('SORT_T2_IJ')
C
      ISYMAB = MULD2H(ISYMA,ISYMB)
      ISYMIJ = MULD2H(ISYMT2,ISYMAB)
      DO ISYMI = 1,NSYM
         ISYMJ = MULD2H(ISYMIJ,ISYMI)
         ISYMAIJ = MULD2H(ISYMA,ISYMIJ)
         ISYMAI  = MULD2H(ISYMA,ISYMI)
         DO I = 1,NRHF(ISYMI)
            DO J = 1,NRHF(ISYMJ)
               KOFF1 =  IMATIJ(ISYMI,ISYMJ)
     *                + NRHF(ISYMI)*(J-1)
     *                + I
C
               KOFF2 = IT2SP(ISYMAIJ,ISYMB)
     *                + NCKI(ISYMAIJ)*(B-1)
     *                + ICKI(ISYMAI,ISYMJ)
     *                + NT1AM(ISYMAI)*(J-1)
     *                + IT1AM(ISYMA,ISYMI)
     *                + NVIR(ISYMA)*(I-1)
     *                + A
C
                     XIJ(KOFF1) = T2TP(KOFF2)
            END DO
         END DO
      END DO
C
      CALL QEXIT('SORT_T2_IJ')
C
      RETURN
      END
C  /* Deck sort_t2_aji */
      SUBROUTINE  SORT_T2_AJI(XAJI,ISYMB,B,T2TP,ISYMT2) 
C
C     t2tp(aijb) as I^I(AJI)
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYMI, ISYMT2
      INTEGER ISYMAJI, ISYMJ, ISYMAJ, ISYMB, ISYMA, ISYMAIJ, ISYMAI
      INTEGER KOFF1, KOFF2
#if defined (SYS_CRAY)
      REAL XAJI(*), T2TP(*) 
#else
      DOUBLE PRECISION XAJI(*), T2TP(*) 
#endif
C
      CALL QENTER('SORT_T2_AJI')
C
      ISYMAIJ = MULD2H(ISYMT2,ISYMB)
      ISYMAJI = MULD2H(ISYMT2,ISYMB)
      DO ISYMJ = 1,NSYM
         ISYMAI = MULD2H(ISYMAJI,ISYMJ)
         DO ISYMA = 1,NSYM 
            ISYMAJ = MULD2H(ISYMJ,ISYMA)
            ISYMI = MULD2H(ISYMAI,ISYMA)
            DO J = 1,NRHF(ISYMJ)
               DO I = 1,NRHF(ISYMI)
                  DO A = 1,NVIR(ISYMA)
                     KOFF1 = ISAIK(ISYMAJ,ISYMI)
     *                      + NT1AM(ISYMAJ)*(I-1)
     *                      + IT1AM(ISYMA,ISYMJ)
     *                      + NVIR(ISYMA)*(J-1) + A
C
                     KOFF2 = IT2SP(ISYMAIJ,ISYMB)
     *                     + NCKI(ISYMAIJ)*(B-1)
     *                     + ICKI(ISYMAI,ISYMJ)
     *                     + NT1AM(ISYMAI)*(J-1)
     *                     + IT1AM(ISYMA,ISYMI)
     *                     + NVIR(ISYMA)*(I-1)
     *                     + A
C
                     XAJI(KOFF1) = T2TP(KOFF2)
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      CALL QEXIT('SORT_T2_AJI')
C
      RETURN
      END
