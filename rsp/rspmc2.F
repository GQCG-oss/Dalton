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
===========================================================================
Revision 1.2  2000/05/24 19:03:29  hjj
new getref calls with appropriate NDREF instead of NCREF
(fixing error for triplet with CSF)
===========================================================================
#endif
      SUBROUTINE E4DRV(VECA,VEC1,VEC2,VEC3,IBCDEQ,
     *                 E4TRS,XINDX,H1,ZYM1,ZYM2,ZYM3,DEN1,DEN2,
     *                 UDV,PV,WRK,LWRK,KZYVR,KZYV1,KZYV2,KZYV3,
     *                 IGRSYM,ISYMV1,ISYMV2,ISYMV3,CMO,FC,MJWOP)
C
C      Purpose:
C      Outer driver routine for E[4] times three vectors.
C
#include <implicit.h>
C
#include <maxorb.h>
#include <infdim.h>
#include <inforb.h>
#include <infpri.h>
#include <infrsp.h>
#include <infvar.h>
#include <wrkrsp.h>
C
      PARAMETER ( DM1 = -1.0D0 , D2 = 2.0D0 , D6 = 6.0D0 )
C
      DIMENSION VECA(KZYVR),VEC1(KZYV1),VEC2(KZYV2),VEC3(KZYV3)
      DIMENSION E4TRS(KZYVR),MJWOP(2,MAXWOP,8),XINDX(*),H1(*)
      DIMENSION ZYM1(*),ZYM2(*),ZYM3(*)
      DIMENSION DEN1(*),DEN2(*),UDV(*),PV(*)
      DIMENSION WRK(*),CMO(*),FC(*)
C
      CALL QENTER('E4DRV')
      IDAX1 = 1
C
      IF (IBCDEQ.EQ.3)  GO TO 200
      IF (IBCDEQ.EQ.24) GO TO 400
C
      CALL ECASE1(VECA,VEC1,VEC2,VEC3,
     *           E4TRS,XINDX,H1,ZYM1,ZYM2,ZYM3,DEN1,DEN2,
     *           UDV,PV,WRK,LWRK,KZYVR,KZYV1,KZYV2,KZYV3,
     *           IGRSYM,ISYMV1,ISYMV2,ISYMV3,CMO,FC,MJWOP)
C
      IF (IBCDEQ.EQ.2) CALL DSCAL(KZYVR,D2,E4TRS,1)
      IF (IBCDEQ.EQ.4) CALL DSCAL(KZYVR,1/D2,E4TRS,1)
C
      CALL ECASE2(VECA,VEC1,VEC2,VEC3,
     *           E4TRS,XINDX,H1,ZYM1,ZYM2,ZYM3,DEN1,DEN2,
     *           UDV,PV,WRK,LWRK,KZYVR,KZYV1,KZYV2,KZYV3,
     *           IGRSYM,ISYMV1,ISYMV2,ISYMV3,CMO,FC,MJWOP)
C
      IF (IBCDEQ.EQ.2) CALL DSCAL(KZYVR,1/D2,E4TRS,1)
      IF (IBCDEQ.EQ.2) GO TO 100
C
      CALL ECASE3(VECA,VEC1,VEC2,VEC3,
     *           E4TRS,XINDX,H1,ZYM1,ZYM2,ZYM3,DEN1,DEN2,
     *           UDV,PV,WRK,LWRK,KZYVR,KZYV1,KZYV2,KZYV3,
     *           IGRSYM,ISYMV1,ISYMV2,ISYMV3,CMO,FC,MJWOP)
C
      CALL ECASE4(VECA,VEC1,VEC2,VEC3,
     *           E4TRS,XINDX,H1,ZYM1,ZYM2,ZYM3,DEN1,DEN2,
     *           UDV,PV,WRK,LWRK,KZYVR,KZYV1,KZYV2,KZYV3,
     *           IGRSYM,ISYMV1,ISYMV2,ISYMV3,CMO,FC,MJWOP)
C
 100  CALL ECASE3(VECA,VEC2, VEC1, VEC3,
     *           E4TRS,XINDX,H1,ZYM2,ZYM1,ZYM3,DEN1,DEN2,
     *           UDV,PV,WRK,LWRK,KZYVR,KZYV2,KZYV1,KZYV3,
     *           IGRSYM,ISYMV2,ISYMV1,ISYMV3,CMO,FC,MJWOP)
C
      CALL ECASE4(VECA,VEC2, VEC1, VEC3,
     *           E4TRS,XINDX,H1,ZYM2,ZYM1,ZYM3,DEN1,DEN2,
     *           UDV,PV,WRK,LWRK,KZYVR,KZYV2,KZYV1,KZYV3,
     *           IGRSYM,ISYMV2,ISYMV1,ISYMV3,CMO,FC,MJWOP)
C
      IF (.NOT.DIROIT) THEN
         CALL OITCLO(IDAX1,'DELETE')
      END IF
C
      IF (IBCDEQ.EQ.2) CALL DSCAL(KZYVR,D2,E4TRS,1)
      IF (IBCDEQ.EQ.2) GO TO 400
C
 200  CALL ECASE1(VECA,VEC2, VEC3, VEC1,
     *           E4TRS,XINDX,H1,ZYM2,ZYM3,ZYM1,DEN1,DEN2,
     *           UDV,PV,WRK,LWRK,KZYVR,KZYV2,KZYV3,KZYV1,
     *           IGRSYM,ISYMV2,ISYMV3,ISYMV1,CMO,FC,MJWOP)
C
      IF (IBCDEQ.EQ.3) CALL DSCAL(KZYVR,1/D2,E4TRS,1)
      IF (IBCDEQ.EQ.4) CALL DSCAL(KZYVR,D2,E4TRS,1)
C
      CALL ECASE2(VECA,VEC2, VEC3, VEC1,
     *           E4TRS,XINDX,H1,ZYM2,ZYM3,ZYM1,DEN1,DEN2,
     *           UDV,PV,WRK,LWRK,KZYVR,KZYV2,KZYV3,KZYV1,
     *           IGRSYM,ISYMV2,ISYMV3,ISYMV1,CMO,FC,MJWOP)
C
      IF (IBCDEQ.EQ.4) CALL DSCAL(KZYVR,1/D2,E4TRS,1)
      IF (IBCDEQ.EQ.4) GO TO 300
C
      CALL ECASE3(VECA,VEC2, VEC3, VEC1,
     *           E4TRS,XINDX,H1,ZYM2,ZYM3,ZYM1,DEN1,DEN2,
     *           UDV,PV,WRK,LWRK,KZYVR,KZYV2,KZYV3,KZYV1,
     *           IGRSYM,ISYMV2,ISYMV3,ISYMV1,CMO,FC,MJWOP)
C
      CALL ECASE4(VECA,VEC2, VEC3, VEC1,
     *           E4TRS,XINDX,H1,ZYM2,ZYM3,ZYM1,DEN1,DEN2,
     *           UDV,PV,WRK,LWRK,KZYVR,KZYV2,KZYV3,KZYV1,
     *           IGRSYM,ISYMV2,ISYMV3,ISYMV1,CMO,FC,MJWOP)
C
 300  CALL ECASE3(VECA,VEC3, VEC2, VEC1,
     *           E4TRS,XINDX,H1,ZYM3,ZYM2,ZYM1,DEN1,DEN2,
     *           UDV,PV,WRK,LWRK,KZYVR,KZYV3,KZYV2,KZYV1,
     *           IGRSYM,ISYMV3,ISYMV2,ISYMV1,CMO,FC,MJWOP)
C
      CALL ECASE4(VECA,VEC3, VEC2, VEC1,
     *           E4TRS,XINDX,H1,ZYM3,ZYM2,ZYM1,DEN1,DEN2,
     *           UDV,PV,WRK,LWRK,KZYVR,KZYV3,KZYV2,KZYV1,
     *           IGRSYM,ISYMV3,ISYMV2,ISYMV1,CMO,FC,MJWOP)
C
      IF (.NOT.DIROIT) THEN
         CALL OITCLO(IDAX1,'DELETE')
      END IF
C
      IF (IBCDEQ.EQ.4) GO TO 600
C
 400  CALL ECASE1(VECA,VEC3, VEC1, VEC2,
     *           E4TRS,XINDX,H1,ZYM3,ZYM1,ZYM2,DEN1,DEN2,
     *           UDV,PV,WRK,LWRK,KZYVR,KZYV3,KZYV1,KZYV2,
     *           IGRSYM,ISYMV3,ISYMV1,ISYMV2,CMO,FC,MJWOP)
C
      IF (IBCDEQ.EQ.2) CALL DSCAL(KZYVR,1/D2,E4TRS,1)
      IF (IBCDEQ.EQ.3) CALL DSCAL(KZYVR,D2,E4TRS,1)
C
      CALL ECASE2(VECA,VEC3, VEC1, VEC2,
     *           E4TRS,XINDX,H1,ZYM3,ZYM1,ZYM2,DEN1,DEN2,
     *           UDV,PV,WRK,LWRK,KZYVR,KZYV3,KZYV1,KZYV2,
     *           IGRSYM,ISYMV3,ISYMV1,ISYMV2,CMO,FC,MJWOP)
C
      IF (IBCDEQ.EQ.3) CALL DSCAL(KZYVR,1/D2,E4TRS,1)
      IF (IBCDEQ.EQ.3) GO TO 500
      IF (IBCDEQ.EQ.24) CALL DSCAL(KZYVR,1/D2,E4TRS,1)
      IF (IBCDEQ.EQ.24) GO TO 500
C
      CALL ECASE3(VECA,VEC3, VEC1, VEC2,
     *           E4TRS,XINDX,H1,ZYM3,ZYM1,ZYM2,DEN1,DEN2,
     *           UDV,PV,WRK,LWRK,KZYVR,KZYV3,KZYV1,KZYV2,
     *           IGRSYM,ISYMV3,ISYMV1,ISYMV2,CMO,FC,MJWOP)
C
      CALL ECASE4(VECA,VEC3, VEC1, VEC2,
     *           E4TRS,XINDX,H1,ZYM3,ZYM1,ZYM2,DEN1,DEN2,
     *           UDV,PV,WRK,LWRK,KZYVR,KZYV3,KZYV1,KZYV2,
     *           IGRSYM,ISYMV3,ISYMV1,ISYMV2,CMO,FC,MJWOP)
C
 500  CALL ECASE3(VECA,VEC1, VEC3, VEC2,
     *           E4TRS,XINDX,H1,ZYM1,ZYM3,ZYM2,DEN1,DEN2,
     *           UDV,PV,WRK,LWRK,KZYVR,KZYV1,KZYV3,KZYV2,
     *           IGRSYM,ISYMV1,ISYMV3,ISYMV2,CMO,FC,MJWOP)
C
      CALL ECASE4(VECA,VEC1, VEC3, VEC2,
     *           E4TRS,XINDX,H1,ZYM1,ZYM3,ZYM2,DEN1,DEN2,
     *           UDV,PV,WRK,LWRK,KZYVR,KZYV1,KZYV3,KZYV2,
     *           IGRSYM,ISYMV1,ISYMV3,ISYMV2,CMO,FC,MJWOP)
C
      IF (.NOT.DIROIT) THEN
         CALL OITCLO(IDAX1,'DELETE')
      END IF
C
      IF (IBCDEQ.EQ.2) CALL DSCAL(KZYVR,D2,E4TRS,1)
      IF (IBCDEQ.EQ.3) CALL DSCAL(KZYVR,D2,E4TRS,1)
 600  IF (IBCDEQ.EQ.4) CALL DSCAL(KZYVR,D2,E4TRS,1)
      IF (IBCDEQ.EQ.24) CALL DSCAL(KZYVR,D6,E4TRS,1)
C
      CALL QEXIT('E4DRV')
      RETURN
      END
      SUBROUTINE ECASE1(VECA,VEC1,VEC2,VEC3,
     *                 E4TRS,XINDX,H1,ZYM1,ZYM2,ZYM3,DEN1,DEN2,
     *                 UDV,PV,WRK,LWRK,KZYVR,KZYV1,KZYV2,KZYV3,
     *                 IGRSYM,ISYMV1,ISYMV2,ISYMV3,CMO,FC,MJWOP)
C
C
#include <implicit.h>
C
      PARAMETER ( DH= 0.5D0, D0 = 0.0D0, D1 = 1.0D0 )
      PARAMETER ( D2 = 2.0D0, D3 = 3.0D0, D6 = 6.0D0 )
C
#include <maxorb.h>
#include <infdim.h>
#include <inforb.h>
#include <wrkrsp.h>
#include <infrsp.h>
#include <infpri.h>
#include <infvar.h>
#include <qrinf.h>
#include <infinp.h>
#include <infopt.h>
#include <infden.h>
C
C
      DIMENSION VECA(KZYVR),VEC1(KZYV1),VEC2(KZYV2),VEC3(KZYV3)
      DIMENSION E4TRS(KZYVR),MJWOP(2,MAXWOP,8),XINDX(*),H1(*)
      DIMENSION ZYM1(*),ZYM2(*),ZYM3(*)
      DIMENSION DEN1(*),DEN2(*),UDV(*),PV(*)
      DIMENSION WRK(*),CMO(*),FC(*)
C
      LOGICAL   TDM, NORHO2, LREF
      LOGICAL   LCON, LORB
C
      IF (MZCONF(ISYMV1).EQ.0 .OR. 
     *    MZCONF(ISYMV2).EQ.0 .OR.
     *    MZCONF(ISYMV3).EQ.0) RETURN
C
      IF (ISYMV2.NE.ISYMV3) RETURN
      CALL QENTER('ECASE1')
C
C     Initialise variables
C
      INTSYM = 1
      NSIM   = 1
      ISPIN  = 0
      TDM    = .TRUE.
      NORHO2 = .FALSE.
      IDAX0 = 0
      IKLVL = 0
C
      KCREF = 1
      KRES  = KCREF + MZCONF(1)
      KFI   = KRES  + KZYVR
      KFREE = KFI   + NORBT*NORBT
      LFREE = LWRK - KFREE
      IF (LFREE.LT.0) CALL ERRWRK('ECASE1',KFREE-1,LWRK)
C
C     Get one electron matrix and reference state
C
      CALL RSPH1(H1,CMO,WRK,LWRK)
      CALL GETREF(WRK(KCREF),MZCONF(1))
C
C     Case 1a
C     =======
C     /  0     \
C     | Sj(1)  | * { <02L|H|03R> + 2/3<0|H|0>(S(2)S(3)'+S(3)S(2)') }
C     |  0     | 
C     \ Sj(1)' /
C
C     F1 = S(2)S(3)' + S(2)'S(3)
C
      NZCONF = MZCONF(ISYMV2)
      NZVAR  = MZVAR(ISYMV2)
      F1 = DDOT(NZCONF,VEC2,1,VEC3(NZVAR+1),1) +
     *     DDOT(NZCONF,VEC3,1,VEC2(NZVAR+1),1)
C
C     F2 = <02L|H|03R> + <03L|H|02R> = 2*<02L|H|03R>
C
C     Construct density matrix
C
      ILSYM  = MULD2H(IREFSY,ISYMV2)
      IRSYM  = MULD2H(IREFSY,ISYMV3)
      NCL    = MZCONF(ISYMV2)
      NCR    = MZCONF(ISYMV3)
      KZVARL = MZYVAR(ISYMV2)
      KZVARR = MZYVAR(ISYMV3)
      LREF   = .FALSE.
      ISYMDN = MULD2H(ILSYM,IRSYM)
      CALL RSPGDM(NSIM,ILSYM,IRSYM,NCL,NCR,KZVARL,KZVARR,VEC2,VEC3,
     *            OVLAP,DEN1,DEN2,ISPIN,ISPIN,TDM,NORHO2,XINDX,
     *            WRK,KFREE,LWRK,LREF)
C
      CALL MELTWO(H1,DEN1,DEN2,OVLAP,ISYMDN,IDAX0,INTSYM,
     *            ISPIN,ISPIN,ISPIN,IKLVL,
     *            DUMMY,IDUM,DUMMY,IDUM,DUMMY,IDUM,
     *            F2,CMO,WRK(KFREE),LFREE)
      FACT = DH*F2 + D2/D3*F1*(EMCSCF-POTNUC)
      NZCONF = MZCONF(ISYMV1)
      NZVAR  = MZVAR(ISYMV1)
      CALL DAXPY(NZCONF,FACT,VEC1,1,E4TRS,1)
      CALL DAXPY(NZCONF,FACT,VEC1(NZVAR+1),1,E4TRS(NZVAR+1),1)
C
      CALL PRIRES(E4TRS,VECA,IGRSYM,'ECASE1a')
C
C     Case 1b
C     =======
C     /  2/3(<0| [qj,H] |01R>  + <01L| [qj,H] |0>) \
C     |  1/6<j| H |01R>                            | * F1
C     |  2/3(<0| [qj+,H] |01R> + <01L| [qj+,H] |0>)|
C     \ -1/6<01L| H |j>                            /
C
C     Construct the density matrix <0L|..|0> + <0|..|0R>
C
      ILSYM  = IREFSY
      IRSYM  = MULD2H(IREFSY,ISYMV1)
      NCL    = MZCONF(1)
      NCR    = MZCONF(ISYMV1)
      KZVARL = MZCONF(1)
      KZVARR = MZYVAR(ISYMV1)
      LREF   = .TRUE.
      ISYMDN = MULD2H(ILSYM,IRSYM)
      CALL DZERO(DEN1,NASHT*NASHT)
      CALL RSPGDM(NSIM,ILSYM,IRSYM,NCL,NCR,KZVARL,KZVARR,
     *         WRK(KCREF),VEC1,OVLAP,DEN1,DEN2,ISPIN,ISPIN,TDM,
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
      NZYVEC = MZYVAR(ISYMV1)
      NZCVEC = MZCONF(ISYMV1)
      CALL DZERO(WRK(KRES),KZYVR)
      CALL DCOPY(NORBT*NORBT,H1,1,WRK(KFI),1)
      CALL RSP2CR(IGRSYM,ISPIN,WRK(KRES),VEC1,NZYVEC,NZCVEC,ISYMV1,
     *        ISYMDN,XINDX,OVLAP,DEN1,DEN2,WRK(KFI),WRK(KFREE),LFREE,
     *               KZYVR,LCON,LORB,LREF,IDAX0,INTSYM,
     *               ISPIN,ISPIN,ISPIN,IKLVL,DUMMY,IDUM,
     *               DUMMY,IDUM,DUMMY,IDUM,CMO,FC,MJWOP)
C
C     Multiply the configurational part by 1/6*F1
C     and the orbital part by 2/3*F1
C
      NZCONF = MZCONF(IGRSYM)
      NZWOPT = MZWOPT(IGRSYM)
      IOFFZC = 0
      IOFFZO = IOFFZC + NZCONF
      IOFFYC = IOFFZO + NZWOPT
      IOFFYO = IOFFYC + NZCONF
      CALL DAXPY(NZCONF,F1/D6,WRK(KRES+IOFFZC),1,E4TRS(1+IOFFZC),1)
      CALL DAXPY(NZCONF,F1/D6,WRK(KRES+IOFFYC),1,E4TRS(1+IOFFYC),1)
C
      CALL DAXPY(NZWOPT,D2*F1/D3,WRK(KRES+IOFFZO),1,E4TRS(1+IOFFZO),1)
      CALL DAXPY(NZWOPT,D2*F1/D3,WRK(KRES+IOFFYO),1,E4TRS(1+IOFFYO),1)
C
      CALL PRIRES(E4TRS,VECA,IGRSYM,'ECASE1b')
C
      CALL QEXIT('ECASE1')
      RETURN
      END
      SUBROUTINE ECASE2(VECA,VEC1,VEC2,VEC3,
     *                 E4TRS,XINDX,H1,ZYM1,ZYM2,ZYM3,DEN1,DEN2,
     *                 UDV,PV,WRK,LWRK,KZYVR,KZYV1,KZYV2,KZYV3,
     *                 IGRSYM,ISYMV1,ISYMV2,ISYMV3,CMO,FC,MJWOP)
C
C
#include <implicit.h>
C
      PARAMETER ( DH= 0.5D0, D0 = 0.0D0, D1 = 1.0D0 )
C
#include <maxorb.h>
#include <infdim.h>
#include <inforb.h>
#include <wrkrsp.h>
#include <infrsp.h>
#include <infpri.h>
#include <infvar.h>
#include <qrinf.h>
#include <infspi.h>
#include <infden.h>
C
      DIMENSION VECA(KZYVR),VEC1(KZYV1),VEC2(KZYV2),VEC3(KZYV3)
      DIMENSION E4TRS(KZYVR),MJWOP(2,MAXWOP,8),XINDX(*),H1(*)
      DIMENSION ZYM1(*),ZYM2(*),ZYM3(*)
      DIMENSION DEN1(*),DEN2(*),UDV(*),PV(*)
      DIMENSION WRK(*),CMO(*),FC(*)
C
      LOGICAL   LCON, LORB
      LOGICAL   TDM, NORHO2, LREF
C
      CALL QENTER('ECASE2')
C
C     Initialise variables and layout some workspace
C
      NSIM = 1
      IH1SYM = 1
      ISPIN  = 0
      TDM    = .TRUE.
      NORHO2 = .FALSE.
      IDAX0  = 0
      IKLVL = 1
C
      KFI     = 1
      KH1X    = KFI     + N2ORBX
      KCREF   = KH1X    + N2ORBX
      KRES    = KCREF   + MZCONF(1)
      KFREE   = KRES    + KZYVR
      LFREE   = LWRK    - KFREE
      IF (LFREE.LT.0) CALL ERRWRK('ECASE2',KFREE-1,LWRK)
C
C     Transform the integrals the first time
C
      IF ( MZWOPT(ISYMV3) .EQ. 0 ) GO TO 9999
C
      CALL GTZYMT(NSIM,VEC3,KZYV3,ISYMV3,ZYM3,MJWOP)
      IF (DIROIT) THEN
         IDAX = IDAX0
         INTSYM = 1
      ELSE
         JTRLVL= 3
         IPROIT  = MAX(IPRRSP/5,2)
         CALL RSPOIT(IDAX0,IDAX1,JTRLVL,ISYMV3,ZYM3,
     *            WRK,1,LFREE,IPROIT)
         IDAX = IDAX1
         INTSYM = ISYMV3
      END IF
      CALL RSPH1(H1,CMO,WRK,LWRK)
      CALL DZERO(WRK(KH1X),N2ORBX)
      CALL OITH1(ISYMV3,ZYM3,H1,WRK(KH1X),IH1SYM)
C
      IF (MZCONF(ISYMV1) .EQ. 0 .OR. MZCONF(ISYMV2) .EQ. 0) GO TO 9999
C
      CALL GETREF(WRK(KCREF),MZCONF(1))
C
C     Case 2a
C     =======
C     /  <01L| [qj,H(k3)] |02R>  + <02L| [qj,H(k3)] |01R>   \
C     |                       0                             |
C     |  <01L| [qj+,H(k3)] |02R> + <02L| [qj+,H(k3)] |01R>  |
C     \                       0                             /
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
      CALL DZERO(DEN1,NASHT*NASHT)
      CALL RSPGDM(NSIM,ILSYM,IRSYM,NCL,NCR,KZVARL,KZVARR,
     *         VEC1,VEC2,OVLAP,DEN1,DEN2,ISPIN,ISPIN,TDM,NORHO2,
     *         XINDX,WRK,KFREE,LFREE,LREF)
C
C     Make the gradient
C
      IF ( MZWOPT(IGRSYM) .GT. 0 ) THEN
         CALL DCOPY(NORBT*NORBT,WRK(KH1X),1,WRK(KFI),1)
         CALL RSP2CR(IGRSYM,ISPIN,E4TRS,DUMMY,IDUM,IDUM,IDUM,ISYMDN,
     *               XINDX,OVLAP,DEN1,DEN2,WRK(KFI),WRK(KFREE),LFREE,
     *               KZYVR,.FALSE.,.TRUE.,LREF,IDAX,INTSYM,
     *               ISPIN,ISPIN,ISPIN,IKLVL,ZYM3,ISYMV3,
     *               DUMMY,IDUM,DUMMY,IDUM,CMO,FC,MJWOP)
      END IF
C
      CALL PRIRES(E4TRS,VECA,IGRSYM,'ECASE2a')
C
C     Case 2b
C     =======
C     /               0                             \
C     | { 1/2<02L|H(k3)|0> + <0|H(k3)|02R> }*Sj(1)  |
C     |               0                             | + permutation
C     \ { <02L|H(k3)|0> + 1/2<0|H(k3)|02R> }*Sj(1)' /
C
C     F1R=<0|H(k3)|-01R>
C
      ILSYM  = IREFSY
      IRSYM  = MULD2H(IREFSY,ISYMV1)
      NCL    = MZCONF(1)
      NCR    = MZCONF(ISYMV1)
      IOFF   = 1
      ISYMDN = MULD2H(ILSYM,IRSYM)
      CALL DZERO(DEN1,NASHT*NASHT)
      CALL RSPDM(ILSYM,IRSYM,NCL,NCR,WRK(KCREF),VEC1(IOFF),
     *           DEN1,DEN2,ISPIN,ISPIN,TDM,NORHO2,XINDX,WRK,
     *           KFREE,LFREE)
C
      OVLAP = D0
      IF (ILSYM.EQ.IRSYM)
     *   OVLAP = DDOT(NCL,WRK(KCREF),1,VEC1(IOFF),1)
      CALL MELTWO(WRK(KH1X),DEN1,DEN2,OVLAP,ISYMDN,IDAX,ISYMV3,
     *            ISPIN,ISPIN,ISPIN,IKLVL,ZYM3,ISYMV3,
     *            DUMMY,IDUM,DUMMY,IDUM,
     *            F1R,CMO,WRK(KFREE),LFREE)
C
C     F2R=<0|H(k3)|-02R>
C
      IRSYM  = MULD2H(IREFSY,ISYMV2)
      NCR    = MZCONF(ISYMV2)
      ISYMDN = MULD2H(ILSYM,IRSYM)
      CALL DZERO(DEN1,NASHT*NASHT)
      CALL RSPDM(ILSYM,IRSYM,NCL,NCR,WRK(KCREF),VEC2(IOFF),
     *           DEN1,DEN2,ISPIN,ISPIN,TDM,NORHO2,XINDX,WRK,
     *           KFREE,LFREE)
C
      OVLAP = D0
      IF (ILSYM.EQ.IRSYM)
     *   OVLAP = DDOT(NCL,WRK(KCREF),1,VEC2(IOFF),1)
      CALL MELTWO(WRK(KH1X),DEN1,DEN2,OVLAP,ISYMDN,IDAX,ISYMV3,
     *            ISPIN,ISPIN,ISPIN,IKLVL,ZYM3,ISYMV3,
     *            DUMMY,IDUM,DUMMY,IDUM,
     *            F2R,CMO,WRK(KFREE),LFREE)
C
C     F1L=<01L|H(k3)|0>
C
      ILSYM  = MULD2H(IREFSY,ISYMV1)
      IRSYM  = IREFSY
      NCL    = MZCONF(ISYMV1)
      NCR    = MZCONF(1)
      IOFF   = MZVAR(ISYMV1) + 1
      ISYMDN = MULD2H(ILSYM,IRSYM)
      CALL DZERO(DEN1,NASHT*NASHT)
      CALL RSPDM(ILSYM,IRSYM,NCL,NCR,VEC1(IOFF),WRK(KCREF),
     *           DEN1,DEN2,ISPIN,ISPIN,TDM,NORHO2,XINDX,WRK,
     *           KFREE,LFREE)
C
      OVLAP = D0
      IF (ILSYM.EQ.IRSYM)
     *   OVLAP = DDOT(NCL,WRK(KCREF),1,VEC1(IOFF),1)
      CALL MELTWO(WRK(KH1X),DEN1,DEN2,OVLAP,ISYMDN,IDAX,ISYMV3,
     *            ISPIN,ISPIN,ISPIN,IKLVL,ZYM3,ISYMV3,
     *            DUMMY,IDUM,DUMMY,IDUM,
     *            F1L,CMO,WRK(KFREE),LFREE)
C
C     F2L=<02L|H(k3)|0>
C
      ILSYM  = MULD2H(IREFSY,ISYMV2)
      NCL    = MZCONF(ISYMV2)
      IOFF   = MZVAR(ISYMV2) + 1
      ISYMDN = MULD2H(ILSYM,IRSYM)
      CALL DZERO(DEN1,NASHT*NASHT)
      CALL RSPDM(ILSYM,IRSYM,NCL,NCR,VEC2(IOFF),WRK(KCREF),
     *           DEN1,DEN2,ISPIN,ISPIN,TDM,NORHO2,XINDX,WRK,
     *           KFREE,LFREE)
C
      OVLAP = D0
      IF (ILSYM.EQ.IRSYM)
     *   OVLAP = DDOT(NCL,WRK(KCREF),1,VEC2(IOFF),1)
      CALL MELTWO(WRK(KH1X),DEN1,DEN2,OVLAP,ISYMDN,IDAX,ISYMV3,
     *            ISPIN,ISPIN,ISPIN,IKLVL,ZYM3,ISYMV3,
     *            DUMMY,IDUM,DUMMY,IDUM,
     *            F2L,CMO,WRK(KFREE),LFREE)
C
      NZCONF = MZCONF(IGRSYM)
      NZVAR  = MZVAR(IGRSYM)
      IF (IGRSYM.EQ.ISYMV1) THEN
         FACT   = DH*F2L - F2R
         CALL DAXPY(NZCONF,FACT,VEC1,1,E4TRS,1)
         FACT   = -DH*F2R + F2L
         CALL DAXPY(NZCONF,FACT,VEC1(NZVAR+1),1,E4TRS(NZVAR+1),1)
      END IF
      IF (IGRSYM.EQ.ISYMV2) THEN
         FACT = DH*F1L - F1R
         CALL DAXPY(NZCONF,FACT,VEC2,1,E4TRS,1)
         FACT = -DH*F1R + F1L
         CALL DAXPY(NZCONF,FACT,VEC2(NZVAR+1),1,E4TRS(NZVAR+1),1)
      END IF
C
      CALL PRIRES(E4TRS,VECA,IGRSYM,'ECASE2b')
C
C     Case 2c
C     =======
C     /   <0| [qj+,H(k3)] |0> \
C     | 1/2<j| H(k3) |0>      | * ( S(1)S(2)' + S(1)'S(2) )
C     |   <0| [qj ,H(k3)] |0> |
C     \ -1/2<0| H(k3) |j>     /
C
C     FACT = S(1)S(2)' + S(1)'S(2)
C
      IF (ISYMV1.NE.ISYMV2) GO TO 9999
C
      NZCONF = MZCONF(ISYMV1)
      NZVAR  = MZVAR(ISYMV1)
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
      CALL DZERO(WRK(KRES),KZYVR)
      CALL DCOPY(NORBT*NORBT,WRK(KH1X),1,WRK(KFI),1)
      CALL RSP2CR(IGRSYM,ISPIN,WRK(KRES),DUMMY,IDUM,IDUM,IDUM,ISYMDN,
     *            XINDX,OVLAP,UDV,PV,WRK(KFI),WRK(KFREE),LFREE,
     *            KZYVR,LCON,LORB,LREF,IDAX,INTSYM,
     *            ISPIN,ISPIN,ISPIN,IKLVL,ZYM3,ISYMV3,
     *            DUMMY,IDUM,DUMMY,IDUM,CMO,FC,MJWOP)
C
      NZCONF = MZCONF(IGRSYM)
      NZWOPT = MZWOPT(IGRSYM)
      IOFFZC = 0
      IOFFZO = IOFFZC + NZCONF
      IOFFYC = IOFFZO + NZWOPT
      IOFFYO = IOFFYC + NZCONF
      CALL DAXPY(NZCONF,DH*FACT,WRK(KRES+IOFFZC),1,E4TRS(1+IOFFZC),1)
      CALL DAXPY(NZCONF,DH*FACT,WRK(KRES+IOFFYC),1,E4TRS(1+IOFFYC),1)
C
      CALL DAXPY(NZWOPT,FACT,WRK(KRES+IOFFZO),1,E4TRS(1+IOFFZO),1)
      CALL DAXPY(NZWOPT,FACT,WRK(KRES+IOFFYO),1,E4TRS(1+IOFFYO),1)
C
      CALL PRIRES(E4TRS,VECA,IGRSYM,'ECASE2c')
C
 9999 CALL QEXIT('ECASE2')
      RETURN
      END
      SUBROUTINE ECASE3(VECA,VEC1,VEC2,VEC3,
     *                 E4TRS,XINDX,H1,ZYM1,ZYM2,ZYM3,DEN1,DEN2,
     *                 UDV,PV,WRK,LWRK,KZYVR,KZYV1,KZYV2,KZYV3,
     *                 IGRSYM,ISYMV1,ISYMV2,ISYMV3,CMO,FC,MJWOP)
#include <implicit.h>
#include <dummy.h>
C
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0, DH = 0.5D0 )
C
#include <maxorb.h>
#include <infdim.h>
#include <inforb.h>
#include <wrkrsp.h>
#include <infrsp.h>
#include <infpri.h>
#include <infvar.h>
#include <qrinf.h>
#include <infspi.h>
#include <infden.h>
C
      DIMENSION VECA(KZYVR),VEC1(KZYV1),VEC2(KZYV2),VEC3(KZYV3)
      DIMENSION E4TRS(KZYVR),MJWOP(2,MAXWOP,8),XINDX(*),H1(*)
      DIMENSION ZYM1(*),ZYM2(*),ZYM3(*)
      DIMENSION DEN1(*),DEN2(*),UDV(*),PV(*)
      DIMENSION WRK(*),CMO(*),FC(*)
C
      LOGICAL   LCON, LORB
      LOGICAL   TDM, NORHO2, LREF
C
      CALL QENTER('ECASE3')
C
C     Initialise variables and layout some workspace
C
      NSIM = 1
      IH1SYM = 1
      ISPIN  = 0
      TDM    = .TRUE.
      NORHO2 = .FALSE.
      IDAX0  = 0
      IDAX1  = 1
      IKLVL = 2
C
      KFI     = 1
      KH1X    = KFI     + N2ORBX
      KH1XX   = KH1X    + N2ORBX
      KCREF   = KH1XX   + N2ORBX
      KRES    = KCREF   + MZCONF(1)
      KFREE   = KRES    + KZYVR
      LFREE   = LWRK    - KFREE
      IF (LFREE.LT.0) CALL ERRWRK('E4DRV',KFREE-1,LWRK)
C
      IF ( MZWOPT(ISYMV3).EQ.0 .OR. MZWOPT(ISYMV2).EQ.0 ) GO TO 9999
C
C     Transform the two-electron integrals a second time
C
      CALL GTZYMT(NSIM,VEC2,KZYV2,ISYMV2,ZYM2,MJWOP)
      CALL GTZYMT(NSIM,VEC3,KZYV3,ISYMV3,ZYM3,MJWOP)
      IF (DIROIT) THEN
         IDAX = IDAX0
         INTSYM =1
      ELSE
         JTRLVL= 2
         IPROIT  = MAX(IPRRSP/5,2)
         CALL RSPOIT(IDAX1,IDAX2,JTRLVL,ISYMV2,ZYM2,
     *            WRK,1,LFREE,IPROIT)
         IDAX = IDAX2
         INTSYM = MULD2H(ISYMV2,ISYMV3)
      END IF
C
C     Transform the one-electron integrals two times
C
      CALL RSPH1(H1,CMO,WRK,LWRK)
      CALL DZERO(WRK(KH1X),N2ORBX)
      CALL OITH1(ISYMV3,ZYM3,H1,WRK(KH1X),IH1SYM)
      CALL DZERO(WRK(KH1XX),N2ORBX)
      CALL OITH1(ISYMV2,ZYM2,WRK(KH1X),WRK(KH1XX),ISYMV3)
C
      IF (MZCONF(ISYMV1) .EQ. 0) GO TO 9999
C
      CALL GETREF(WRK(KCREF),MZCONF(1))
C
C     Case 3a
C     =======
C     /   <0| [qj,H(k2,k3)] |01R>  + <01L| [qj,H(k2,k3)] |0>  \
C     |   <j| H(k2,k3) |01R>                                  | * 1/2
C     |   <0| [qj+,H(k2,k3)] |01R> + <01L| [qj+,H(k2,k3)] |0> |
C     \  -<01L| H(k2,k3) |j>                                  /
C
C     Construct the density matrix <0L|..|0> + <0|..|0R>
C
      ILSYM  = IREFSY
      IRSYM  = MULD2H(IREFSY,ISYMV1)
      NCL    = MZCONF(1)
      NCR    = MZCONF(ISYMV1)
      KZVARL = MZCONF(1)
      KZVARR = MZYVAR(ISYMV1)
      LREF   = .TRUE.
      ISYMDN = MULD2H(ILSYM,IRSYM)
      CALL DZERO(DEN1,NASHT*NASHT)
      CALL RSPGDM(NSIM,ILSYM,IRSYM,NCL,NCR,KZVARL,KZVARR,
     *         WRK(KCREF),VEC1,OVLAP,DEN1,DEN2,ISPIN,ISPIN,TDM,
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
      NZYVEC = MZYVAR(ISYMV1)
      NZCVEC = MZCONF(ISYMV1)
      CALL DZERO(WRK(KRES),KZYVR)
      CALL DCOPY(NORBT*NORBT,WRK(KH1XX),1,WRK(KFI),1)
      CALL RSP2CR(IGRSYM,ISPIN,WRK(KRES),VEC1,NZYVEC,NZCVEC,ISYMV1,
     *     ISYMDN,XINDX,OVLAP,DEN1,DEN2,WRK(KFI),WRK(KFREE),LFREE,
     *            KZYVR,LCON,LORB,LREF,IDAX,INTSYM,
     *            ISPIN,ISPIN,ISPIN,IKLVL,ZYM3,ISYMV3,ZYM2,ISYMV2,
     *            DUMMY,IDUMMY,CMO,FC,MJWOP)
C
      CALL DAXPY(KZYVR,DH,WRK(KRES),1,E4TRS,1)
C
      CALL PRIRES(E4TRS,VECA,IGRSYM,'ECASE3a')
C
C     Case 3b
C     =======
C     /   0    \
C     | Sj(1)  | * 1/2<0| H(k2,k3) |0>
C     |   0    |
C     \ Sj(1)' /
C
      IF (IGRSYM.EQ.ISYMV1) THEN
         ISYMDN = 1
         OVLAP = D1
         CALL MELTWO(WRK(KH1XX),UDV,PV,OVLAP,ISYMDN,IDAX,INTSYM,
     *               ISPIN,ISPIN,ISPIN,IKLVL,ZYM3,ISYMV3,
     *               ZYM2,ISYMV2,DUMMY,IDUM,
     *               FACT,CMO,WRK(KFREE),LFREE)
         NZCONF = MZCONF(IGRSYM)
         NZVAR  = MZVAR(IGRSYM)
         CALL DAXPY(NZCONF,DH*FACT,VEC1,1,E4TRS,1)
         CALL DAXPY(NZCONF,DH*FACT,VEC1(NZVAR+1),1,E4TRS(NZVAR+1),1)
      END IF
C
      CALL PRIRES(E4TRS,VECA,IGRSYM,'ECASE3b')
C
 9999 CALL QEXIT('ECASE3')
      RETURN
      END
      SUBROUTINE ECASE4(VECA,VEC1,VEC2,VEC3,
     *                 E4TRS,XINDX,H1,ZYM1,ZYM2,ZYM3,DEN1,DEN2,
     *                 UDV,PV,WRK,LWRK,KZYVR,KZYV1,KZYV2,KZYV3,
     *                 IGRSYM,ISYMV1,ISYMV2,ISYMV3,CMO,FC,MJWOP)
#include <implicit.h>
C
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0, D6 = 6.0D0 )
C
#include <maxorb.h>
#include <infdim.h>
#include <inforb.h>
#include <wrkrsp.h>
#include <infrsp.h>
#include <infpri.h>
#include <infvar.h>
#include <qrinf.h>
#include <infspi.h>
#include <infden.h>
C
      DIMENSION VECA(KZYVR),VEC1(KZYV1),VEC2(KZYV2),VEC3(KZYV3)
      DIMENSION E4TRS(KZYVR),MJWOP(2,MAXWOP,8),XINDX(*),H1(*)
      DIMENSION ZYM1(*),ZYM2(*),ZYM3(*)
      DIMENSION DEN1(*),DEN2(*),UDV(*),PV(*)
      DIMENSION WRK(*),CMO(*),FC(*)
C
      LOGICAL   LCON, LORB
      LOGICAL   TDM, NORHO2, LREF
C
      CALL QENTER('ECASE4')
C
C     Initialise variables and layout some workspace
C
      NSIM = 1
      IH1SYM = 1
      ISPIN  = 0
      TDM    = .TRUE.
      NORHO2 = .FALSE.
      IDAX0  = 0
      IDAX1  = 1
      IDAX2  = 2
      IKLVL = 3
C
      KH1X    = 1
      KH1XX   = KH1X    + N2ORBX
      KH1XXX  = KH1XX   + N2ORBX
      KRES    = KH1XXX  + N2ORBX
      KFREE   = KRES    + KZYVR
      LFREE   = LWRK    - KFREE
      IF (LFREE.LT.0) CALL ERRWRK('E4DRV',KFREE-1,LWRK)
C
      IF ((MZWOPT(ISYMV1) .EQ. 0).OR.(MZWOPT(ISYMV2) .EQ. 0)
     *     .OR.(MZWOPT(ISYMV3) .EQ. 0)) THEN
         IF (.NOT.DIROIT) THEN
            CALL OITCLO(IDAX2,'DELETE')
         END IF
         GO TO 9999
      END IF
C
C     Transform the integrals a third time
C
      CALL GTZYMT(NSIM,VEC1,KZYV1,ISYMV1,ZYM1,MJWOP)
      CALL GTZYMT(NSIM,VEC2,KZYV2,ISYMV2,ZYM2,MJWOP)
      CALL GTZYMT(NSIM,VEC3,KZYV3,ISYMV3,ZYM3,MJWOP)
C
      IF (DIROIT) THEN
         IDAX = IDAX0
         INTSYM = 1
      ELSE
         JTRLVL= 1
         IPROIT  = MAX(IPRRSP/5,2)
         CALL RSPOIT(IDAX2,IDAX3,JTRLVL,ISYMV1,ZYM1,
     *            WRK,1,LFREE,IPROIT)
         CALL OITCLO(IDAX2,'DELETE')
         IDAX = IDAX3
         INTSYM = IGRSYM
      END IF
C
C     Transform the one-electron integrals three times
C
      CALL RSPH1(H1,CMO,WRK,LWRK)
      CALL DZERO(WRK(KH1X),N2ORBX)
      CALL OITH1(ISYMV3,ZYM3,H1,WRK(KH1X),IH1SYM)
      CALL DZERO(WRK(KH1XX),N2ORBX)
      CALL OITH1(ISYMV2,ZYM2,WRK(KH1X),WRK(KH1XX),ISYMV3)
      CALL DZERO(WRK(KH1XXX),N2ORBX)
      ISYM23 = MULD2H(ISYMV2,ISYMV3)
      CALL OITH1(ISYMV1,ZYM1,WRK(KH1XX),WRK(KH1XXX),ISYM23)
C
C     We have the density matrices over the reference state already
C
      OVLAP = D1
      ISYMDN = 1
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
      CALL DZERO(WRK(KRES),KZYVR)
      CALL RSP2CR(IGRSYM,ISPIN,WRK(KRES),DUMMY,IDUM,IDUM,IDUM,
     *     ISYMDN,XINDX,OVLAP,UDV,PV,WRK(KH1XXX),WRK(KFREE),LFREE,
     *            KZYVR,LCON,LORB,LREF,IDAX,INTSYM,ISPIN,ISPIN,ISPIN,
     *            IKLVL,ZYM3,ISYMV3,ZYM2,ISYMV2,ZYM1,ISYMV1,CMO,FC,
     *            MJWOP)
C
      CALL DAXPY(KZYVR,(D1/D6),WRK(KRES),1,E4TRS,1)
C
      CALL PRIRES(E4TRS,VECA,IGRSYM,'ECASE4')
C
      IF (.NOT.DIROIT) THEN
         CALL OITCLO(IDAX3,'DELETE')
      END IF
C
 9999 CALL QEXIT('ECASE4')
      RETURN
      END
      SUBROUTINE PRIRES(RESV,VECA,ISYMV,TEXT)
C
#include <implicit.h>
C
#include <priunit.h>
#include <infrsp.h>
#include <qrinf.h>
C
      CHARACTER*(*) TEXT
      DIMENSION RESV(*),VECA(*)
C
C
      IF (IPRRSP.GE.100) THEN
         NZVAR  = MZVAR(ISYMV)
         NZCONF = MZCONF(ISYMV)
         NZWOPT = MZWOPT(ISYMV)
C
         VALORB = DDOT(NZWOPT,VECA(NZCONF+1),1,RESV(NZCONF+1),1) +
     &      DDOT(NZWOPT,VECA(NZVAR+NZCONF+1),1,RESV(NZVAR+NZCONF+1),1)
         VALCNF = DDOT(NZCONF,VECA(1),1,RESV(1),1) +
     &      DDOT(NZCONF,VECA(NZVAR+1),1,RESV(NZVAR+1),1)
         WRITE(LUPRI,'(/2A,2(/A,F20.8))') 'Accumulated result in ',TEXT,
     &        'Orbital part:', VALORB, 'Config part :', VALCNF
      END IF
C
      RETURN
      END
