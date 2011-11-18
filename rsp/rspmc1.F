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
#ifdef REV_LOG
C===========================================================================
CRevision 1.2  2000/05/24 19:02:25  hjj
Cnew GETREF calls, fixing CSF problem for triplet
C===========================================================================
#endif
      SUBROUTINE A3INIT(KZYVR,KZYV1,KZYV2,IGRSYM,ISYMV1,ISYMV2,
     *                 OPLBL,IOPSYM,VEC1,VEC2,A3TRS,XINDX,UDV,
     *                 CMO,MJWOP,WRK,LWRK)
C
C     Layout the core for the calculation of A3 times two vectors
C
#include <implicit.h>
#include <priunit.h>
#include <infdim.h>
#include <inforb.h>
#include <maxorb.h>
#include <maxash.h>
#include <infvar.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <infpri.h>
C
      CHARACTER*8 OPLBL
C
      DIMENSION WRK(*)
      DIMENSION A3TRS(KZYVR), MJWOP(2,MAXWOP,8)
      DIMENSION VEC1(KZYV1)
      DIMENSION VEC2(KZYV2)
      DIMENSION XINDX(*)
      DIMENSION CMO(*)
      DIMENSION UDV(NASHDI,NASHDI)
C
C     Initialise the gradient to zero.
C
      CALL DZERO(A3TRS,KZYVR)
C
      KOPMAT  = 1
      KDEN1 = KOPMAT + NORBT * NORBT
      KFREE = KDEN1  + NASHT * NASHT
      LWRKF   = LWRK - KFREE + 1
      IF (LWRKF.LT.0) CALL ERRWRK('A3INIT',KFREE-1,LWRK)
C
      IF (IPRRSP .GT. 50 ) THEN
         WRITE(LUPRI,'(//A)')  ' Characteristics of A3 gradient'
         WRITE(LUPRI,'(A)')    ' =============================='
         WRITE(LUPRI,'(A,I8)') ' Gradient symmetry :',IGRSYM
         WRITE(LUPRI,'(A,I8)') ' Gradient length   :',KZYVR
         WRITE(LUPRI,'(A,I8)') ' Vector1 symmetry   :',ISYMV1
         WRITE(LUPRI,'(A,I8)') ' Vector1 length     :',KZYV1
         WRITE(LUPRI,'(A,I8)') ' Vector2 symmetry   :',ISYMV2
         WRITE(LUPRI,'(A,I8)') ' Vector2 length     :',KZYV2
         WRITE(LUPRI,'(A,I8)') ' Operator symmetry :',IOPSYM
         WRITE(LUPRI,'(A,A8)') ' Operator label    :',OPLBL
         WRITE(LUPRI,'(A//)')  ' =============================='
      END IF
C
      CALL A3DRV(KZYVR,KZYV1,KZYV2,IGRSYM,ISYMV1,ISYMV2,
     *                 OPLBL,IOPSYM,VEC1,VEC2,A3TRS,WRK(KOPMAT),
     *                 WRK(KDEN1),UDV,XINDX,CMO,MJWOP,WRK(KFREE),LWRKF)
C
      RETURN
      END
      SUBROUTINE A3DRV(KZYVR,KZYV1,KZYV2,IGRSYM,ISYMV1,ISYMV2,
     *            OPLBL,IOPSYM,VEC1,VEC2,A3TRS,OPMAT,
     *            DEN1,UDV,XINDX,CMO,MJWOP,WRK,LWRK)
C
C      Purpose:
C      Outer driver routine for X[3] times two vectors. This subroutine
C      calls the setup of the operator transformation, constructs a den-
C      sity matrix if necessary  and calls RSP1GR to compute the gradient.
C
#include <implicit.h>
C
#include <maxorb.h>
#include <infdim.h>
#include <inforb.h>
#include <infvar.h>
#include <qrinf.h>
#include <infrsp.h>
C
      CHARACTER*8 OPLBL
C
      DIMENSION A3TRS(KZYVR), MJWOP(2,MAXWOP,8)
      DIMENSION XINDX(*)
      DIMENSION OPMAT(NORBT,NORBT)
      DIMENSION DEN1(NASHDI,NASHDI)
      DIMENSION UDV(NASHDI,NASHDI)
      DIMENSION WRK(*)
      DIMENSION VEC1(KZYV1)
      DIMENSION VEC2(KZYV2)
      DIMENSION CMO(*)
C
C     Get the operator matrix
C
      KSYMP = -1
      CALL PRPGET(OPLBL,CMO,OPMAT,KSYMP,ANTSYM,WRK,LWRK,IPRRSP)
      IF (KSYMP.NE.IOPSYM) CALL QUIT(
     &   'A3DRV: unexpected symmetry of operator matrix')
      CALL DGETRN(OPMAT,NORBT,NORBT)
C
      CALL ACASE1(KZYVR,KZYV1,KZYV2,IGRSYM,ISYMV1,ISYMV2,
     *            IOPSYM,VEC1,VEC2,A3TRS,OPMAT,
     *            DEN1,UDV,XINDX,WRK,LWRK)
C
      CALL ACASE2(KZYVR,KZYV1,KZYV2,IGRSYM,ISYMV1,ISYMV2,
     *            IOPSYM,VEC1,VEC2,A3TRS,OPMAT,
     *            DEN1,UDV,XINDX,MJWOP,WRK,LWRK)
C
      CALL ACASE3(KZYVR,KZYV1,KZYV2,IGRSYM,ISYMV1,ISYMV2,
     *            IOPSYM,VEC1,VEC2,A3TRS,OPMAT,
     *            DEN1,UDV,XINDX,MJWOP,WRK,LWRK)
C
      CALL ACASE1(KZYVR,KZYV2,KZYV1,IGRSYM,ISYMV2,ISYMV1,
     *            IOPSYM,VEC2,VEC1,A3TRS,OPMAT,
     *            DEN1,UDV,XINDX,WRK,LWRK)
C
      CALL ACASE2(KZYVR,KZYV2,KZYV1,IGRSYM,ISYMV2,ISYMV1,
     *            IOPSYM,VEC2,VEC1,A3TRS,OPMAT,
     *            DEN1,UDV,XINDX,MJWOP,WRK,LWRK)
C
      CALL ACASE3(KZYVR,KZYV2,KZYV1,IGRSYM,ISYMV2,ISYMV1,
     *            IOPSYM,VEC2,VEC1,A3TRS,OPMAT,
     *            DEN1,UDV,XINDX,MJWOP,WRK,LWRK)
C
C     Swap the final result to conform with the notation for A[3]
C
      CALL DSWAP(MZVAR(IGRSYM),A3TRS,1,
     *           A3TRS(MZVAR(IGRSYM)+1),1)
C
      RETURN
      END
      SUBROUTINE ACASE1(KZYVR,KZYV1,KZYV2,IGRSYM,ISYMV1,ISYMV2,
     *                 IOPSYM,VEC1,VEC2,A3TRS,OPMAT,
     *                 DEN1,UDV,XINDX,WRK,LWRK)
C
#include <implicit.h>
C
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0, D3 = 3D0, D6 = 6.0D0 )
C
#include <maxorb.h>
#include <maxash.h>
#include <inforb.h>
#include <infind.h>
#include <infvar.h>
#include <infdim.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <rspprp.h>
#include <infhyp.h>
#include <qrinf.h>
#include <infpri.h>
#include <infspi.h>
C
      DIMENSION A3TRS(KZYVR)
      DIMENSION XINDX(*)
      DIMENSION OPMAT(NORBT,NORBT)
      DIMENSION DEN1(NASHDI,NASHDI)
      DIMENSION UDV(NASHDI,NASHDI)
      DIMENSION WRK(*)
      DIMENSION VEC1(KZYV1)
      DIMENSION VEC2(KZYV2)
C
      LOGICAL   LCON, LORB
      LOGICAL   TDM, NORHO2, LREF
C
      IF (MZCONF(ISYMV1) .EQ. 0 .OR. MZCONF(ISYMV2) .EQ. 0) RETURN
C
C     Initialise variables and layout some workspace
C
      TDM    = .TRUE.
      NORHO2 = .TRUE.
      NSIM = 1
      INTSYM = 1
      ISPIN = 0
      IPRONE = 75
C
      KCREF  = 1
C     KRES   = KCREF + NCREF
      KRES   = KCREF + MZCONF(1)
      KFREE  = KRES  + KZYVR
      LFREE  = LWRK  - KFREE
      IF (LFREE.LT.0) CALL ERRWRK('ACASE1',KFREE-1,LWRK)
C
      CALL GETREF(WRK(KCREF),MZCONF(1))
C
C     Case 1a corresponding to X3 does not exist
C     =======
C
C     Case 1b
C     =======
C     /               0                         \
C     | { -1/3<02L|A|0> - 1/6<0|A|02R> }*Sj(1)' |
C     |               0                         |
C     \ { -1/6<02L|A|0> - 1/3<0|A|02R> }*Sj(1)  /
C
C     F2R=<0|A|-02R>
C
      ILSYM  = IREFSY
      IRSYM  = MULD2H(IREFSY,ISYMV2)
      NCL    = MZCONF(1)
      NCR    = MZCONF(ISYMV2)
      IOFF   = 1
      CALL DZERO(DEN1,NASHT*NASHT)
      CALL RSPDM(ILSYM,IRSYM,NCL,NCR,WRK(KCREF),VEC2(IOFF),
     *           DEN1,DUMMY,ISPIN,ISPIN,TDM,NORHO2,XINDX,WRK,
     *           KFREE,LFREE)
      OVLAP = D0
      IF (ILSYM.EQ.IRSYM) 
     *   OVLAP = DDOT(NCL,WRK(KCREF),1,VEC2(IOFF),1)
C
      CALL MELONE(OPMAT,IOPSYM,DEN1,OVLAP,F2R,IPRONE,'F2R in ACASE1')
C
C     F2L=<02L|A|0>
C
      ILSYM  = MULD2H(IREFSY,ISYMV2)
      IRSYM  = IREFSY
      NCL    = MZCONF(ISYMV2)
      NCR    = MZCONF(1)
      IOFF   = MZVAR(ISYMV2) + 1
      CALL DZERO(DEN1,NASHT*NASHT)
      CALL RSPDM(ILSYM,IRSYM,NCL,NCR,VEC2(IOFF),WRK(KCREF),
     *           DEN1,DUMMY,ISPIN,ISPIN,TDM,NORHO2,XINDX,WRK,
     *           KFREE,LFREE)
      OVLAP = D0
      IF (ILSYM.EQ.IRSYM) 
     *   OVLAP = DDOT(NCL,WRK(KCREF),1,VEC2(IOFF),1)
C
      CALL MELONE(OPMAT,IOPSYM,DEN1,OVLAP,F2L,IPRONE,'F2L in ACASE1')
C
      NZCONF = MZCONF(IGRSYM)
      NZVAR  = MZVAR(IGRSYM)
      IF (IGRSYM.EQ.ISYMV1) THEN
         FACT   = - F2L/D3 + F2R/D6
         CALL DAXPY(NZCONF,FACT,VEC1,1,A3TRS,1)
         FACT   = - F2L/D6 + F2R/D3 
         CALL DAXPY(NZCONF,FACT,VEC1(NZVAR+1),1,A3TRS(NZVAR+1),1)
      END IF
C
C     Case 1c
C     =======
C     /      0       \
C     | -<0| A |j>   | * -1/6 * S(1)S(2)'
C     |      0       |
C     \  <j| A |0>   /
C
      IF (ISYMV1.NE.ISYMV2) RETURN
C
      NZCONF = MZCONF(ISYMV1)
      NZVAR = MZVAR(ISYMV1)
      FACT = DDOT(NZCONF,VEC1,1,VEC2(NZVAR+1),1)
      FACT = -FACT/D6
C
      ISYMST = MULD2H(IGRSYM,IREFSY)
      IF(ISYMST .EQ. IREFSY ) THEN
         LCON = ( MZCONF(IGRSYM) .GT. 1 )
      ELSE
         LCON = ( MZCONF(IGRSYM) .GT. 0 )
      END IF
      IF (LCON) THEN
         LREF = .TRUE.
         CALL DZERO(WRK(KRES),KZYVR)
         CALL CONSX(NSIM,KZYVR,IGRSYM,OPMAT,WRK(KCREF),MZCONF(1),
     *              MZCONF(1),IREFSY,MZCONF(IGRSYM),ISYMST,LREF,
     *              WRK(KRES),XINDX,WRK(KFREE),LFREE)
         CALL DAXPY(KZYVR,FACT,WRK(KRES),1,A3TRS,1)
      END IF
C
      RETURN
      END
      SUBROUTINE ACASE2(KZYVR,KZYV1,KZYV2,IGRSYM,ISYMV1,ISYMV2,
     *                 IOPSYM,VEC1,VEC2,A3TRS,OPMAT,
     *                 DEN1,UDV,XINDX,MJWOP,WRK,LWRK)
C
#include <implicit.h>
C
      PARAMETER ( D1 = 1.0D0, DMH = -0.5D0 )
C
#include <maxorb.h>
#include <maxash.h>
#include <inforb.h>
#include <infind.h>
#include <infvar.h>
#include <infdim.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <rspprp.h>
#include <infhyp.h>
#include <qrinf.h>
#include <infpri.h>
#include <infspi.h>
C
      DIMENSION A3TRS(KZYVR), MJWOP(2,MAXWOP,8)
      DIMENSION XINDX(*)
      DIMENSION OPMAT(NORBT,NORBT)
      DIMENSION DEN1(NASHDI,NASHDI)
      DIMENSION UDV(NASHDI,NASHDI)
      DIMENSION WRK(*)
      DIMENSION VEC1(KZYV1), VEC2(KZYV2)
C
      LOGICAL   LCON, LREF
C
C     Initialise variables and layout some workspace
C
      NSIM = 1
      IPRONE = 75
C
      KZYM   = 1
      KA3X   = KZYM  + N2ORBX
      KFREE  = KA3X  + N2ORBX
      LFREE  = LWRK  - KFREE
      IF (LFREE.LT.0) CALL ERRWRK('ACASE1',KFREE-1,LWRK)
C
      IF (MZCONF(ISYMV2) .EQ. 0) RETURN
      IF ( MZWOPT(ISYMV1) .EQ. 0 ) RETURN
C
C     Transform the operator with kappa and put the factor -1/2
C     present in these terms into the operator
C
      CALL GTZYMT(NSIM,VEC1,KZYV1,ISYMV1,WRK(KZYM),MJWOP)
      CALL DZERO(WRK(KA3X),N2ORBX)
      CALL OITH1(ISYMV1,WRK(KZYM),OPMAT,WRK(KA3X),IOPSYM)
      CALL DSCAL(NORBT*NORBT,DMH,WRK(KA3X),1)
C
C     Case 2a
C     =======
C     /           0          \
C     |  1/2<02L| A(k1) |j>  |
C     |           0          |
C     \ -1/2<j| A(k1) |02R>  /
C
C     Make the gradient
C
      ISYMST = MULD2H(IGRSYM,IREFSY)
      IF ( ISYMST .EQ. IREFSY ) THEN
         LCON = ( MZCONF(IGRSYM) .GT. 1 )
      ELSE
         LCON = ( MZCONF(IGRSYM) .GT. 0 )
      END IF
      LREF = .FALSE.
      NZYVEC = MZYVAR(ISYMV2)
      NZCVEC = MZCONF(ISYMV2)
      IF (LCON) THEN
         CALL CONSX(NSIM,KZYVR,IGRSYM,WRK(KA3X),VEC2,NZYVEC,NZCVEC,
     *              ISYMV2,MZCONF(IGRSYM),ISYMST,LREF,A3TRS,
     *              XINDX,WRK(KFREE),LFREE)
      END IF
C
C     Case 2b
C     =======
C     /   0    \
C     | Sj(2)' | * -1/2<0| A(k1) |0>
C     |   0    |
C     \ Sj(2)  /
C
      IF (IGRSYM.EQ.ISYMV2) THEN
         OVLAP = D1
         CALL MELONE(WRK(KA3X),1,UDV,OVLAP,FACT,IPRONE,'FACT in ACASE2')
         NZCONF = MZCONF(IGRSYM)
         NZVAR  = MZVAR(IGRSYM)
         CALL DAXPY(NZCONF,FACT,VEC2,1,A3TRS,1)
         CALL DAXPY(NZCONF,FACT,VEC2(NZVAR+1),1,A3TRS(NZVAR+1),1)
      END IF
C
      RETURN
      END
      SUBROUTINE ACASE3(KZYVR,KZYV1,KZYV2,IGRSYM,ISYMV1,ISYMV2,
     *                 IOPSYM,VEC1,VEC2,A3TRS,OPMAT,
     *                 DEN1,UDV,XINDX,MJWOP,WRK,LWRK)
C
#include <implicit.h>
C
      PARAMETER ( D0 = 0.0D0, DMH = -0.5D0 , D1 = 1.0D0, D3 = 3.0D0 )
C
#include <maxorb.h>
#include <maxash.h>
#include <inforb.h>
#include <infind.h>
#include <infvar.h>
#include <infdim.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <rspprp.h>
#include <infhyp.h>
#include <qrinf.h>
#include <infpri.h>
#include <infspi.h>
C
      DIMENSION A3TRS(KZYVR), MJWOP(2,MAXWOP,8)
      DIMENSION XINDX(*)
      DIMENSION OPMAT(NORBT,NORBT)
      DIMENSION DEN1(NASHDI,NASHDI)
      DIMENSION UDV(NASHDI,NASHDI)
      DIMENSION WRK(*)
      DIMENSION VEC1(KZYV1), VEC2(KZYV2)
C
      LOGICAL   LCON, LORB, LREF
C
C     Initialise variables and layout some workspace
C
      ISPIN  = 0
C
      KCREF  = 1
      KRES   = KCREF + MZCONF(1)
      KZYM   = KRES  + KZYVR
      KA3X   = KZYM  + N2ORBX
      KA3XX  = KA3X  + N2ORBX
      KFREE  = KA3XX + N2ORBX
      LFREE  = LWRK  - KFREE
      IF (LFREE.LT.0) CALL ERRWRK('ACASE3',KFREE-1,LWRK)
C
      IF ((MZWOPT(ISYMV2) .EQ. 0).OR.(MZWOPT(ISYMV1) .EQ. 0)) RETURN
C
C     / 1/3<0| [qj+,A(k1,k2)] |0> \
C     |   -<0| A(k1,k2) |j>       | * -1/2 
C     | 1/3<0| [qj ,A(k1,k2)] |0> |
C     \    <j| A(k1,k2) |0>       /
C
C     Transform the operator with 2k,1k
C
C     Put the factor of minus one half present in this term into the
C     ZY matrix used for transforming the integrals
C
      NSIM = 1
      CALL GTZYMT(NSIM,VEC1,KZYV1,ISYMV1,WRK(KZYM),MJWOP)
      CALL DSCAL(NORBT*NORBT,DMH,WRK(KZYM),1)
      CALL DZERO(WRK(KA3X),N2ORBX)
      CALL OITH1(ISYMV1,WRK(KZYM),OPMAT,WRK(KA3X),IOPSYM)
      CALL GTZYMT(NSIM,VEC2,KZYV2,ISYMV2,WRK(KZYM),MJWOP)
      CALL DZERO(WRK(KA3XX),N2ORBX)
      ISYMX = MULD2H(IOPSYM,ISYMV1)
      CALL OITH1(ISYMV2,WRK(KZYM),WRK(KA3X),WRK(KA3XX),ISYMX)
C
      CALL GETREF(WRK(KCREF),MZCONF(1))
C
C     We have the density matrix over the reference state already
C
      ISYMDN = 1
      OVLAP  = D1
C
C     Make the gradient
C
      ISYMST = MULD2H(IGRSYM,IREFSY)
      IF ( ISYMST .EQ. IREFSY ) THEN
         LCON = ( MZCONF(IGRSYM) .GT. 1 )
      ELSE
         LCON = ( MZCONF(IGRSYM) .GT. 0 )
      END IF
      LORB   = ( MZWOPT(IGRSYM) .GT. 0 )
      LREF = .TRUE.
      IF (LCON) THEN
         CALL DZERO(WRK(KRES),KZYVR)
         CALL CONSX(NSIM,KZYVR,IGRSYM,WRK(KA3XX),WRK(KCREF),MZCONF(1),
     *              MZCONF(1),IREFSY,MZCONF(IGRSYM),ISYMST,LREF,
     *              WRK(KRES),XINDX,WRK(KFREE),LFREE)
         CALL DAXPY(KZYVR,D1,WRK(KRES),1,A3TRS,1)
      END IF
      IF (LORB) THEN
         CALL DZERO(WRK(KRES),KZYVR)
         CALL ORBSX(NSIM,IGRSYM,KZYVR,WRK(KRES),WRK(KA3XX),OVLAP,
     *              ISYMDN,UDV,MJWOP,WRK(KFREE),LFREE)
         CALL DAXPY(KZYVR,1/D3,WRK(KRES),1,A3TRS,1)
      END IF
C
      RETURN
      END
      SUBROUTINE X3INIT(KZYVR,KZYV1,KZYV2,IGRSYM,ISYMV1,ISYMV2,
     *                 OPLBL,IOPSYM,VEC1,VEC2,X3TRS,XINDX,UDV,
     *                 CMO,MJWOP,WRK,LWRK)
C
C     Layout the core for the calculation of X3 times two vectors
C
#include <implicit.h>
#include <priunit.h>
#include <infdim.h>
#include <inforb.h>
#include <maxorb.h>
#include <maxash.h>
#include <infvar.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <rspprp.h>
#include <infhyp.h>
#include <qrinf.h>
#include <infpri.h>
#include <infspi.h>
#include <infcr.h>
C
      CHARACTER*8 OPLBL
C
      DIMENSION WRK(*)
      DIMENSION X3TRS(KZYVR), MJWOP(2,MAXWOP,8)
      DIMENSION VEC1(KZYV1), VEC2(KZYV2)
      DIMENSION XINDX(*)
      DIMENSION CMO(*)
      DIMENSION UDV(NASHDI,NASHDI)
C
C     Initialise the gradient to zero.
C
      CALL DZERO(X3TRS,KZYVR)
C
      KOPMAT  = 1
      KDEN1 = KOPMAT + NORBT * NORBT
      KFREE = KDEN1  + NASHT * NASHT
      LWRKF   = LWRK - KFREE + 1
      IF (LWRKF.LT.0) CALL ERRWRK('X3INIT',KFREE-1,LWRK)
C
      IF (IPRRSP .GT. 100 ) THEN
         WRITE(LUPRI,'(//A)')  ' Characteristics of X3 gradient'
         WRITE(LUPRI,'(A)')    ' =============================='
         WRITE(LUPRI,'(A,I8)') ' Gradient symmetry :',IGRSYM
         WRITE(LUPRI,'(A,I8)') ' Gradient length   :',KZYVR
         WRITE(LUPRI,'(A,I8)') ' Vector1 symmetry   :',ISYMV1
         WRITE(LUPRI,'(A,I8)') ' Vector1 length     :',KZYV1
         WRITE(LUPRI,'(A,I8)') ' Vector2 symmetry   :',ISYMV2
         WRITE(LUPRI,'(A,I8)') ' Vector2 length     :',KZYV2
         WRITE(LUPRI,'(A,I8)') ' Operator symmetry :',IOPSYM
         WRITE(LUPRI,'(A,A8)') ' Operator label    :',OPLBL
         WRITE(LUPRI,'(A//)')  ' =============================='
      END IF
C
      CALL X3DRV(KZYVR,KZYV1,KZYV2,IGRSYM,ISYMV1,ISYMV2,
     *                 OPLBL,IOPSYM,VEC1,VEC2,X3TRS,WRK(KOPMAT),
     *                 WRK(KDEN1),UDV,XINDX,CMO,MJWOP,WRK(KFREE),LWRKF)
C
      RETURN
      END
      SUBROUTINE X3DRV(KZYVR,KZYV1,KZYV2,IGRSYM,ISYMV1,ISYMV2,
     *            OPLBL,IOPSYM,VEC1,VEC2,X3TRS,OPMAT,
     *            DEN1,UDV,XINDX,CMO,MJWOP,WRK,LWRK)
C
C      Purpose:
C      Outer driver routine for X[3] times two vectors. This subroutine
C      calls the setup of the operator transformation, constructs a den-
C      sity matrix if necessary  and calls RSP1GR to compute the gradient.
C
#include <implicit.h>
C
#include <maxorb.h>
#include <infdim.h>
#include <inforb.h>
#include <infvar.h>
#include <infrsp.h>
C
      CHARACTER*8 OPLBL
C
      DIMENSION X3TRS(KZYVR), MJWOP(2,MAXWOP,8)
      DIMENSION XINDX(*)
      DIMENSION OPMAT(NORBT,NORBT)
      DIMENSION DEN1(NASHDI,NASHDI), UDV(NASHDI,NASHDI)
      DIMENSION WRK(*)
      DIMENSION VEC1(KZYV1), VEC2(KZYV2)
      DIMENSION CMO(*)
C
C
C     Get the operator matrix
C
      KSYMP = -1
      CALL PRPGET(OPLBL,CMO,OPMAT,KSYMP,ANTSYM,WRK,LWRK,IPRRSP)
      IF (KSYMP.NE.IOPSYM) CALL QUIT(
     &   'X3DRV: unexpected symmetry of operator matrix')
C
      CALL XCASE1(KZYVR,KZYV1,KZYV2,IGRSYM,ISYMV1,ISYMV2,
     *            IOPSYM,VEC1,VEC2,X3TRS,OPMAT,
     *            DEN1,UDV,XINDX,MJWOP,WRK,LWRK)
C
      CALL XCASE2(KZYVR,KZYV1,KZYV2,IGRSYM,ISYMV1,ISYMV2,
     *            IOPSYM,VEC1,VEC2,X3TRS,OPMAT,
     *            DEN1,UDV,XINDX,MJWOP,WRK,LWRK)
C
      CALL XCASE3(KZYVR,KZYV1,KZYV2,IGRSYM,ISYMV1,ISYMV2,
     *            IOPSYM,VEC1,VEC2,X3TRS,OPMAT,
     *            DEN1,UDV,XINDX,MJWOP,WRK,LWRK)
C
      CALL XCASE2(KZYVR,KZYV2,KZYV1,IGRSYM,ISYMV2,ISYMV1,
     *            IOPSYM,VEC2,VEC1,X3TRS,OPMAT,
     *            DEN1,UDV,XINDX,MJWOP,WRK,LWRK)
C
      CALL XCASE3(KZYVR,KZYV2,KZYV1,IGRSYM,ISYMV2,ISYMV1,
     *            IOPSYM,VEC2,VEC1,X3TRS,OPMAT,
     *            DEN1,UDV,XINDX,MJWOP,WRK,LWRK)
C
      RETURN
      END
      SUBROUTINE XCASE1(KZYVR,KZYV1,KZYV2,IGRSYM,ISYMV1,ISYMV2,
     *                 IOPSYM,VEC1,VEC2,X3TRS,OPMAT,
     *                 DEN1,UDV,XINDX,MJWOP,WRK,LWRK)
C
#include <implicit.h>
C
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0, DH = 0.5D0 )
C
#include <maxorb.h>
#include <maxash.h>
#include <inforb.h>
#include <infind.h>
#include <infvar.h>
#include <infdim.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <rspprp.h>
#include <infhyp.h>
#include <qrinf.h>
#include <infpri.h>
#include <infspi.h>
C
      DIMENSION X3TRS(KZYVR), MJWOP(2,MAXWOP,8)
      DIMENSION XINDX(*)
      DIMENSION OPMAT(NORBT,NORBT)
      DIMENSION DEN1(NASHDI,NASHDI)
      DIMENSION UDV(NASHDI,NASHDI)
      DIMENSION WRK(*)
      DIMENSION VEC1(KZYV1), VEC2(KZYV2)
C
      LOGICAL   LCON, LORB, TDM, NORHO2, LREF
C
C     Initialise variables and layout some workspace
C
      TDM    = .TRUE.
      NORHO2 = .TRUE.
      NSIM = 1
      ISPIN = 0
      IPRONE = 75
C
      KCREF  = 1
      KRES   = KCREF + MZCONF(1)
      KFREE  = KRES  + KZYVR
      LFREE  = LWRK  - KFREE
      IF (LFREE.LT.0) CALL ERRWRK('XCASE1',KFREE-1,LWRK)
C
      IF (MZCONF(ISYMV1) .EQ. 0 .OR. MZCONF(ISYMV2) .EQ. 0) RETURN
C
      CALL GETREF(WRK(KCREF),MZCONF(1))
C
C     Case 1a
C     =======
C     /   <01L| [qj,X] |02R>  + <02L| [qj,X] |01R>  \
C     |                       0                      |
C     |   <01L| [qj+,X] |02R> + <02L| [qj+,X] |01R>  |
C     \                       0                     /
C
C     Construct <01L|..|02R> + <02L|..|01R> density
C
      ILSYM  = MULD2H(IREFSY,ISYMV1)
      IRSYM  = MULD2H(IREFSY,ISYMV2)
      NCL    = MZCONF(ISYMV1)
      NCR    = MZCONF(ISYMV2)
      KZVARL = MZYVAR(ISYMV1)
      KZVARR = MZYVAR(ISYMV2)
      LREF   = .FALSE.
      ISYMDN = MULD2H(ILSYM,IRSYM)
      CALL DZERO(DEN1,N2ASHX)
      CALL RSPGDM(NSIM,ILSYM,IRSYM,NCL,NCR,KZVARL,KZVARR,
     *         VEC1,VEC2,OVLAP,DEN1,DUMMY,ISPIN,ISPIN,TDM,NORHO2,
     *         XINDX,WRK,KFREE,LFREE,LREF)
C
C     Make the gradient
C
C
      IF ( MZWOPT(IGRSYM) .GT. 0 ) THEN
         CALL ORBSX(NSIM,IGRSYM,KZYVR,X3TRS,OPMAT,OVLAP,ISYMDN,
     *              DEN1,MJWOP,WRK(KFREE),LFREE)
      END IF
C
C     Case 1b
C     =======
C     /               0                     \
C     | { 1/2<02L|X|0> + <0|X|02R> }*Sj(1)  |
C     |               0                     | + permutation
C     \ { <02L|X|0> + 1/2<0|X|02R> }*Sj(1)' /
C
C     F1R=<0|X|-01R>
C
      ILSYM  = IREFSY
      IRSYM  = MULD2H(IREFSY,ISYMV1)
      NCL    = MZCONF(1)
      NCR    = MZCONF(ISYMV1)
      IOFF   = 1
      CALL DZERO(DEN1,NASHT*NASHT)
      CALL RSPDM(ILSYM,IRSYM,NCL,NCR,WRK(KCREF),VEC1(IOFF),
     *           DEN1,DUMMY,ISPIN,ISPIN,TDM,NORHO2,XINDX,WRK,
     *           KFREE,LFREE)
      OVLAP = D0
      IF (ILSYM.EQ.IRSYM) 
     *   OVLAP = DDOT(NCL,WRK(KCREF),1,VEC1(IOFF),1)
C
      CALL MELONE(OPMAT,IOPSYM,DEN1,OVLAP,F1R,IPRONE,'F1R in XCASE1')
C
C     F2R=<0|X|-02R>
C
      IRSYM  = MULD2H(IREFSY,ISYMV2)
      NCR    = MZCONF(ISYMV2)
      CALL DZERO(DEN1,NASHT*NASHT)
      CALL RSPDM(ILSYM,IRSYM,NCL,NCR,WRK(KCREF),VEC2(IOFF),
     *           DEN1,DUMMY,ISPIN,ISPIN,TDM,NORHO2,XINDX,WRK,
     *           KFREE,LFREE)
      OVLAP = D0
      IF (ILSYM.EQ.IRSYM) 
     *   OVLAP = DDOT(NCL,WRK(KCREF),1,VEC2(IOFF),1)
C
      CALL MELONE(OPMAT,IOPSYM,DEN1,OVLAP,F2R,IPRONE,'F2R in XCASE1')
C
C     F1L=<01L|X|0>
C
      ILSYM  = MULD2H(IREFSY,ISYMV1)
      IRSYM  = IREFSY
      NCL    = MZCONF(ISYMV1)
      NCR    = MZCONF(1)
      IOFF   = MZVAR(ISYMV1) + 1
      CALL DZERO(DEN1,NASHT*NASHT)
      CALL RSPDM(ILSYM,IRSYM,NCL,NCR,VEC1(IOFF),WRK(KCREF),
     *           DEN1,DUMMY,ISPIN,ISPIN,TDM,NORHO2,XINDX,WRK,
     *           KFREE,LFREE)
      OVLAP = D0
      IF (ILSYM.EQ.IRSYM) 
     *   OVLAP = DDOT(NCL,WRK(KCREF),1,VEC1(IOFF),1)
C
      CALL MELONE(OPMAT,IOPSYM,DEN1,OVLAP,F1L,IPRONE,'F1L in XCASE1')
C
C     F2L=<02L|X|0>
C
      ILSYM  = MULD2H(IREFSY,ISYMV2)
      NCL    = MZCONF(ISYMV2)
      IOFF   = MZVAR(ISYMV2) + 1
      CALL DZERO(DEN1,NASHT*NASHT)
      CALL RSPDM(ILSYM,IRSYM,NCL,NCR,VEC2(IOFF),WRK(KCREF),
     *           DEN1,DUMMY,ISPIN,ISPIN,TDM,NORHO2,XINDX,WRK,
     *           KFREE,LFREE)
      OVLAP = D0
      IF (ILSYM.EQ.IRSYM) 
     *   OVLAP = DDOT(NCL,WRK(KCREF),1,VEC2(IOFF),1)
C
      CALL MELONE(OPMAT,IOPSYM,DEN1,OVLAP,F2L,IPRONE,'F2L in XCASE1')
C
      NZCONF = MZCONF(IGRSYM)
      NZVAR  = MZVAR(IGRSYM)
      IF (IGRSYM.EQ.ISYMV1) THEN 
         FACT   = DH*F2L - F2R
         CALL DAXPY(NZCONF,FACT,VEC1,1,X3TRS,1)
         FACT   = -DH*F2R + F2L
         CALL DAXPY(NZCONF,FACT,VEC1(NZVAR+1),1,X3TRS(NZVAR+1),1)
      END IF
      IF (IGRSYM.EQ.ISYMV2) THEN 
         FACT   = DH*F1L - F1R
         CALL DAXPY(NZCONF,FACT,VEC2,1,X3TRS,1)
         FACT   = -DH*F1R + F1L
         CALL DAXPY(NZCONF,FACT,VEC2(NZVAR+1),1,X3TRS(NZVAR+1),1)
      END IF
C
C     Case 1c
C     =======
C     /  <0| [qj,X] |0>   \
C     | 1/2<j| X |0>      | * ( S(1)S(2)' + S(1)'S(2) )
C     |  <0| [qj+ ,X] |0> |
C     \ -1/2<0| X |j>     /
C
      IF (ISYMV1.NE.ISYMV2) RETURN
C
      NZCONF = MZCONF(ISYMV1)
      NZVAR = MZVAR(ISYMV1)
      FACT = DDOT(NZCONF,VEC1,1,VEC2(NZVAR+1),1) +
     *       DDOT(NZCONF,VEC2,1,VEC1(NZVAR+1),1)
C
      ISYMDN = 1
      OVLAP  = D1
      ISYMST = MULD2H(IGRSYM,IREFSY)
      IF(ISYMST .EQ. IREFSY ) THEN
         LCON = ( MZCONF(IGRSYM) .GT. 1 )
      ELSE
         LCON = ( MZCONF(IGRSYM) .GT. 0 )
      END IF
      LORB   = ( MZWOPT(IGRSYM) .GT. 0 )
      LREF = .TRUE.
      IF (LCON) THEN
         CALL DZERO(WRK(KRES),KZYVR)
         CALL CONSX(NSIM,KZYVR,IGRSYM,OPMAT,WRK(KCREF),
     *              MZCONF(1),MZCONF(1),IREFSY,MZCONF(IGRSYM),ISYMST,
     *              LREF,WRK(KRES),XINDX,WRK(KFREE),LFREE)
         CALL DSCAL(KZYVR,DH,WRK(KRES),1)
         CALL DAXPY(KZYVR,FACT,WRK(KRES),1,X3TRS,1)
      END IF
      IF (LORB) THEN
         CALL DZERO(WRK(KRES),KZYVR)
         CALL ORBSX(NSIM,IGRSYM,KZYVR,WRK(KRES),OPMAT,OVLAP,ISYMDN,
     *              UDV,MJWOP,WRK(KFREE),LFREE)
         CALL DAXPY(KZYVR,FACT,WRK(KRES),1,X3TRS,1)
      END IF
C
      RETURN
      END
      SUBROUTINE XCASE2(KZYVR,KZYV1,KZYV2,IGRSYM,ISYMV1,ISYMV2,
     *                 IOPSYM,VEC1,VEC2,X3TRS,OPMAT,
     *                 DEN1,UDV,XINDX,MJWOP,WRK,LWRK)
C
#include <implicit.h>
C
      PARAMETER ( D1 = 1.0D0 )
C
#include <maxorb.h>
#include <maxash.h>
#include <inforb.h>
#include <infind.h>
#include <infvar.h>
#include <infdim.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <rspprp.h>
#include <infhyp.h>
#include <qrinf.h>
#include <infpri.h>
#include <infspi.h>
C
      DIMENSION X3TRS(KZYVR), MJWOP(2,MAXWOP,8)
      DIMENSION XINDX(*)
      DIMENSION OPMAT(NORBT,NORBT)
      DIMENSION DEN1(NASHDI,NASHDI)
      DIMENSION UDV(NASHDI,NASHDI)
      DIMENSION WRK(*)
      DIMENSION VEC1(KZYV1)
      DIMENSION VEC2(KZYV2)
C
      LOGICAL   LCON, LORB
      LOGICAL   TDM, NORHO2, LREF
C
C     Initialise variables and layout some workspace
C
      ISPIN  = 0
      TDM    = .TRUE.
      NORHO2 = .TRUE.
      IPRONE = 75
C
      KCREF  = 1
      KZYM   = KCREF + MZCONF(1)
      KX3X   = KZYM  + N2ORBX
      KFREE  = KX3X  + N2ORBX
      LFREE  = LWRK  - KFREE
      IF (LFREE.LT.0) CALL ERRWRK('XCASE1',KFREE-1,LWRK)
C
      IF ( MZCONF(ISYMV2) .EQ. 0 ) RETURN
      IF ( MZWOPT(ISYMV1) .EQ. 0 ) RETURN
C
C     Transform the operator with kappa
C
      NSIM = 1
      CALL GTZYMT(NSIM,VEC1,KZYV1,ISYMV1,WRK(KZYM),MJWOP)
      CALL DZERO(WRK(KX3X),N2ORBX)
      CALL OITH1(ISYMV1,WRK(KZYM),OPMAT,WRK(KX3X),IOPSYM)
C
      CALL GETREF(WRK(KCREF),MZCONF(1))
C
C     Case 2a
C     =======
C     /   <0| [qj,X(k1)] |02R>  + <02L| [qj,X(k1)] |0>  \
C     |   <j| X(k1) |02R>                                |
C     |   <0| [qj+,X(k1)] |02R> + <02L| [qj+,X(k1)] |0>  |
C     \  -<02L| X(k1) |j>                               /
C
C     Construct the density matrix <0L|..|0> + <0|..|0R>
C
      ILSYM  = IREFSY
      IRSYM  = MULD2H(IREFSY,ISYMV2)
      NCL    = MZCONF(1)
      NCR    = MZCONF(ISYMV2)
      KZVARL = MZCONF(1)
      KZVARR = MZYVAR(ISYMV2)
      LREF   = .TRUE.
      ISYMDN = MULD2H(ILSYM,IRSYM)
      CALL DZERO(DEN1,NASHT*NASHT)
      CALL RSPGDM(NSIM,ILSYM,IRSYM,NCL,NCR,KZVARL,KZVARR,
     *         WRK(KCREF),VEC2,OVLAP,DEN1,DUMMY,ISPIN,ISPIN,TDM,
     *         NORHO2,XINDX,WRK,KFREE,LFREE,LREF)
C
C     Make the gradient
C
      ISYMST = MULD2H(IGRSYM,IREFSY)
      IF ( ISYMST .EQ. IREFSY ) THEN
         LCON = ( MZCONF(IGRSYM) .GT. 1 )
      ELSE
         LCON = ( MZCONF(IGRSYM) .GT. 0 )
      END IF
      LORB   = ( MZWOPT(IGRSYM) .GT. 0 )
      LREF = .FALSE.
      NZYVEC = MZYVAR(ISYMV2)
      NZCVEC = MZCONF(ISYMV2)
      CALL RSP1GR(NSIM,KZYVR,IDUM,ISPIN,IGRSYM,ISPIN,ISYMV2,X3TRS,
     *            VEC2,NZYVEC,NZCVEC,OVLAP,ISYMDN,DEN1,WRK(KX3X),
     *            XINDX,MJWOP,WRK(KFREE),LFREE,LORB,LCON,LREF)
C
C     Case 2b
C     =======
C     /   0    \
C     | Sj(2)  | * <0| X(k1) |0>
C     |   0    |
C     \ Sj(2)' /
C
      IF (IGRSYM.EQ.ISYMV2) THEN
         OVLAP = D1
         CALL MELONE(WRK(KX3X),1,UDV,OVLAP,FACT,IPRONE,'FACT in XCASE2')
         NZCONF = MZCONF(IGRSYM)
         NZVAR  = MZVAR(IGRSYM)
         CALL DAXPY(NZCONF,FACT,VEC2,1,X3TRS,1)
         CALL DAXPY(NZCONF,FACT,VEC2(NZVAR+1),1,X3TRS(NZVAR+1),1)
      END IF
C
      RETURN
      END
      SUBROUTINE XCASE3(KZYVR,KZYV1,KZYV2,IGRSYM,ISYMV1,ISYMV2,
     *                 IOPSYM,VEC1,VEC2,X3TRS,OPMAT,
     *                 DEN1,UDV,XINDX,MJWOP,WRK,LWRK)
C
#include <implicit.h>
C
      PARAMETER ( D0 = 0.0D0, DH = 0.5D0 , D1 = 1.0D0 )
C
#include <maxorb.h>
#include <maxash.h>
#include <inforb.h>
#include <infind.h>
#include <infvar.h>
#include <infdim.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <rspprp.h>
#include <infhyp.h>
#include <qrinf.h>
#include <infpri.h>
#include <infspi.h>
C
      DIMENSION X3TRS(KZYVR), MJWOP(2,MAXWOP,8)
      DIMENSION XINDX(*)
      DIMENSION OPMAT(NORBT,NORBT)
      DIMENSION DEN1(NASHDI,NASHDI)
      DIMENSION UDV(NASHDI,NASHDI)
      DIMENSION WRK(*)
      DIMENSION VEC1(KZYV1), VEC2(KZYV2)
C
      LOGICAL   LCON, LORB, TDM, LREF
C
C     Initialise variables and layout some workspace
C
      ISPIN  = 0
      TDM    = .TRUE.
C
      KCREF  = 1
      KZYM   = KCREF + MZCONF(1)
      KX3X   = KZYM  + N2ORBX
      KX3XX  = KX3X  + N2ORBX
      KFREE  = KX3XX + N2ORBX
      LFREE  = LWRK  - KFREE
      IF (LFREE.LT.0) CALL ERRWRK('XCASE1',KFREE-1,LWRK)
C
      IF ((MZWOPT(ISYMV2) .EQ. 0).OR.(MZWOPT(ISYMV1) .EQ. 0)) RETURN
C
C     / <0| [qj ,X(k1,k2)] |0> \
C     | <j| X(k1,k2) |0>       | * 1/2
C     | <0| [qj+,X(k1,k2)] |0> |
C     \ -<0| X(k1,k2) |j>      /
C
C     Transform the operator with 2k,1k
C
C     Put the factor of one half present in this term into the
C     ZY matrix used for transforming the integrals
C
      NSIM = 1
      CALL GTZYMT(NSIM,VEC1,KZYV1,ISYMV1,WRK(KZYM),MJWOP)
      CALL DSCAL(NORBT*NORBT,DH,WRK(KZYM),1)
      CALL DZERO(WRK(KX3X),N2ORBX)
      CALL OITH1(ISYMV1,WRK(KZYM),OPMAT,WRK(KX3X),IOPSYM)
      CALL GTZYMT(NSIM,VEC2,KZYV2,ISYMV2,WRK(KZYM),MJWOP)
      CALL DZERO(WRK(KX3XX),N2ORBX)
      ISYMX = MULD2H(IOPSYM,ISYMV1)
      CALL OITH1(ISYMV2,WRK(KZYM),WRK(KX3X),WRK(KX3XX),ISYMX)
C
      CALL GETREF(WRK(KCREF),MZCONF(1))
C
C     We have the density matrices over the reference state already
C
      ISYMDN = 1
      OVLAP  = D1
C
C     Make the gradient
C
      ISYMV  = IREFSY
      ISYMST = MULD2H(IGRSYM,IREFSY)
      IF ( ISYMST .EQ. IREFSY ) THEN
         LCON = ( MZCONF(IGRSYM) .GT. 1 )
      ELSE
         LCON = ( MZCONF(IGRSYM) .GT. 0 )
      END IF
      LORB   = ( MZWOPT(IGRSYM) .GT. 0 )
      LREF = .TRUE.
      NZYVEC = MZCONF(1)
      NZCVEC = MZCONF(1)
      CALL RSP1GR(NSIM,KZYVR,IDUM,ISPIN,IGRSYM,ISPIN,ISYMV,X3TRS,
     *            WRK(KCREF),NZYVEC,NZCVEC,OVLAP,ISYMDN,UDV,WRK(KX3XX),
     *            XINDX,MJWOP,WRK(KFREE),LFREE,LORB,LCON,LREF)
C
      RETURN
      END
