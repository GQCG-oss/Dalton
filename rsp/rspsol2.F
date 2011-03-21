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
#ifdef REVLOG
===========================================================================
Revision 1.2  2000/05/24 19:09:55  hjj
inserted Dalton header, some CREF changes for CSF w. trilet response
===========================================================================
#endif
C
      SUBROUTINE E4SOL(VECA,VEC1,VEC2,VEC3,ETRS,XINDX,ZYM1,ZYM2,ZYM3,
     *              DEN1,UDV,WRK,LWRK,KZYVR,KZYV1,KZYV2,KZYV3,
     *              IGRSYM,ISYMV1,ISYMV2,ISYMV3,CMO,MJWOP,ISYRLM)
C
C     Purpose:
C     Outer driver routine for solvent contribution 
C     to E[4] times three vectors.
C
#include <implicit.h>
C
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infdim.h>
#include <infinp.h>
#include <infvar.h>
C
      DIMENSION ETRS(KZYVR),XINDX(*)
      DIMENSION UDV(NASHDI,NASHDI),DEN1(NASHDI,NASHDI)
      DIMENSION ZYM1(*),ZYM2(*),ZYM3(*),WRK(*),CMO(*)
      DIMENSION VEC1(KZYV1),VEC2(KZYV2),VEC3(KZYV3),VECA(KZYVR)
      DIMENSION MJWOP(2,MAXWOP,8), ISYRLM(2*LSOLMX+1)
C
      NSIM = 1
      KFREE = 1
      LFREE = LWRK
      CALL MEMGET('REAL',KW   ,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KA   ,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KA1  ,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KA2  ,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KB1  ,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KA12 ,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KA123,N2ORBX,WRK,KFREE,LFREE)
C
      CALL GTZYMT(NSIM,VEC1,KZYV1,ISYMV1,ZYM1,MJWOP)
      CALL GTZYMT(NSIM,VEC2,KZYV2,ISYMV2,ZYM2,MJWOP)
      CALL GTZYMT(NSIM,VEC3,KZYV3,ISYMV3,ZYM3,MJWOP)
C
      VAL = DDOT(KZYVR,VECA,1,ETRS,1)
C
      CALL W4MAT(VEC1,VEC2,VEC3,ZYM1,ZYM2,ZYM3,ISYMV1,ISYMV2,ISYMV3,
     *           DEN1,UDV,WRK(KFREE),LFREE,CMO,XINDX,ISYRLM,
     *           KZYVR,KZYV1,KZYV2,KZYV3,
     *           WRK(KW),WRK(KA),WRK(KA1),WRK(KA2),
     *           WRK(KB1),WRK(KA12),WRK(KA123))
C
      CALL WCASE1(VECA,VEC1,VEC2,VEC3,ETRS,XINDX,ISYRLM,
     *            WRK(KW),WRK(KA),WRK(KA1),WRK(KA2),WRK(KA12),
     *            DEN1,UDV,WRK(KFREE),LFREE,ZYM1,ZYM2,
     *            KZYVR,KZYV1,KZYV2,KZYV3,IGRSYM,ISYMV1,ISYMV2,ISYMV3)
C
      CALL WCASE2(VECA,VEC1,VEC2,VEC3,ETRS,XINDX,ISYRLM,
     *            WRK(KW),WRK(KA),WRK(KA1),WRK(KA2),WRK(KA12),
     *            DEN1,UDV,WRK(KFREE),LFREE,CMO,MJWOP,ZYM1,ZYM2,
     *            KZYVR,KZYV1,KZYV2,KZYV3,IGRSYM,ISYMV1,ISYMV2,ISYMV3)
C
      CALL WCASE3(VECA,VEC1,VEC2,VEC3,ETRS,XINDX,ISYRLM,
     *            WRK(KW),WRK(KA1),
     *            DEN1,UDV,WRK(KFREE),LFREE,CMO,MJWOP,ZYM1,
     *            KZYVR,KZYV1,KZYV2,KZYV3,IGRSYM,ISYMV1,ISYMV2,ISYMV3)
C
      CALL WCASE4(VECA,VEC1,VEC2,VEC3,ETRS,XINDX,ISYRLM,
     *     WRK(KW),WRK(KA),WRK(KA1),WRK(KB1),WRK(KA12),WRK(KA123),
     *            DEN1,UDV,WRK(KFREE),LFREE,CMO,MJWOP,ZYM1,ZYM2,ZYM3,
     *            KZYVR,KZYV1,KZYV2,KZYV3,IGRSYM,ISYMV1,ISYMV2,ISYMV3)
C
      CALL W4MAT(VEC3,VEC1,VEC2,ZYM3,ZYM1,ZYM2,ISYMV3,ISYMV1,ISYMV2,
     *           DEN1,UDV,WRK(KFREE),LFREE,CMO,XINDX,ISYRLM,
     *           KZYVR,KZYV3,KZYV1,KZYV2,
     *           WRK(KW),WRK(KA),WRK(KA1),WRK(KA2),
     *           WRK(KB1),WRK(KA12),WRK(KA123))
C
      CALL WCASE1(VECA,VEC3,VEC1,VEC2,ETRS,XINDX,ISYRLM,
     *            WRK(KW),WRK(KA),WRK(KA1),WRK(KA2),WRK(KA12),
     *            DEN1,UDV,WRK(KFREE),LFREE,ZYM3,ZYM1,
     *            KZYVR,KZYV3,KZYV1,KZYV2,IGRSYM,ISYMV3,ISYMV1,ISYMV2)
C
      CALL WCASE2(VECA,VEC3,VEC1,VEC2,ETRS,XINDX,ISYRLM,
     *            WRK(KW),WRK(KA),WRK(KA1),WRK(KA2),WRK(KA12),
     *            DEN1,UDV,WRK(KFREE),LFREE,CMO,MJWOP,ZYM3,ZYM1,
     *            KZYVR,KZYV3,KZYV1,KZYV2,IGRSYM,ISYMV3,ISYMV1,ISYMV2)
C
      CALL WCASE3(VECA,VEC3,VEC1,VEC2,ETRS,XINDX,ISYRLM,
     *            WRK(KW),WRK(KA1),
     *            DEN1,UDV,WRK(KFREE),LFREE,CMO,MJWOP,ZYM3,
     *            KZYVR,KZYV3,KZYV1,KZYV2,IGRSYM,ISYMV3,ISYMV1,ISYMV2)
C
      CALL WCASE4(VECA,VEC3,VEC1,VEC2,ETRS,XINDX,ISYRLM,
     *     WRK(KW),WRK(KA),WRK(KA1),WRK(KB1),WRK(KA12),WRK(KA123),
     *            DEN1,UDV,WRK(KFREE),LFREE,CMO,MJWOP,ZYM3,ZYM1,ZYM2,
     *            KZYVR,KZYV3,KZYV1,KZYV2,IGRSYM,ISYMV3,ISYMV1,ISYMV2)
C
      CALL W4MAT(VEC2,VEC3,VEC1,ZYM2,ZYM3,ZYM1,ISYMV2,ISYMV3,ISYMV1,
     *           DEN1,UDV,WRK(KFREE),LFREE,CMO,XINDX,ISYRLM,
     *           KZYVR,KZYV2,KZYV3,KZYV1,
     *           WRK(KW),WRK(KA),WRK(KA1),WRK(KA2),
     *           WRK(KB1),WRK(KA12),WRK(KA123))
C
      CALL WCASE1(VECA,VEC2,VEC3,VEC1,ETRS,XINDX,ISYRLM,
     *            WRK(KW),WRK(KA),WRK(KA1),WRK(KA2),WRK(KA12),
     *            DEN1,UDV,WRK(KFREE),LFREE,ZYM2,ZYM3,
     *            KZYVR,KZYV2,KZYV3,KZYV1,IGRSYM,ISYMV2,ISYMV3,ISYMV1)
C
      CALL WCASE2(VECA,VEC2,VEC3,VEC1,ETRS,XINDX,ISYRLM,
     *            WRK(KW),WRK(KA),WRK(KA1),WRK(KA2),WRK(KA12),
     *            DEN1,UDV,WRK(KFREE),LFREE,CMO,MJWOP,ZYM2,ZYM3,
     *            KZYVR,KZYV2,KZYV3,KZYV1,IGRSYM,ISYMV2,ISYMV3,ISYMV1)
C
      CALL WCASE3(VECA,VEC2,VEC3,VEC1,ETRS,XINDX,ISYRLM,
     *            WRK(KW),WRK(KA1),
     *            DEN1,UDV,WRK(KFREE),LFREE,CMO,MJWOP,ZYM2,
     *            KZYVR,KZYV2,KZYV3,KZYV1,IGRSYM,ISYMV2,ISYMV3,ISYMV1)
C
      CALL WCASE4(VECA,VEC2,VEC3,VEC1,ETRS,XINDX,ISYRLM,
     *     WRK(KW),WRK(KA),WRK(KA1),WRK(KB1),WRK(KA12),WRK(KA123),
     *            DEN1,UDV,WRK(KFREE),LFREE,CMO,MJWOP,ZYM2,ZYM3,ZYM1,
     *            KZYVR,KZYV2,KZYV3,KZYV1,IGRSYM,ISYMV2,ISYMV3,ISYMV1)
C
      WRITE (LUPRI,'(A,F20.12)') 'Total contribution from E4SOL:',
     *                DDOT(KZYVR,VECA,1,ETRS,1) - VAL
      RETURN
      END
      SUBROUTINE W4MAT(VEC1,VEC2,VEC3,ZYM1,ZYM2,ZYM3,
     *                 ISYMV1,ISYMV2,ISYMV3,
     *                 DEN1,UDV,WRK,LFREE,CMO,XINDX,ISYRLM,
     *                 KZYVR,KZYV1,KZYV2,KZYV3,
     *                 W,A,A1,A2,B1,A12,A123)
C
C     Purpose:
C
#include <implicit.h>
C
      PARAMETER ( DH = 0.5D0, D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0 )
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
#include <infinp.h>
#include <dummy.h>
C
      DIMENSION UDV(NASHDI,NASHDI),DEN1(NASHDI,NASHDI)
      DIMENSION VEC1(*),VEC2(*),VEC3(*),ZYM1(*),ZYM2(*),ZYM3(*)
      DIMENSION WRK(*),CMO(*),XINDX(*),ISYRLM(*)
      DIMENSION W(N2ORBX),A(N2ORBX),A1(N2ORBX),A2(N2ORBX)
      DIMENSION B1(N2ORBX),A12(N2ORBX),A123(N2ORBX)
C
      LOGICAL   TDM, LREF, NORHO2
C
C     Initialise variables
C
      JSPIN  = 0
      TDM    = .TRUE.
      KFREE = 1
      NORHO2 = .TRUE.
      NSIM = 1
      IPRLVL = 10
C
      CALL DZERO(W,N2ORBX)
      CALL DZERO(A,N2ORBX)
      CALL DZERO(A1,N2ORBX)
      CALL DZERO(A2,N2ORBX)
      CALL DZERO(B1,N2ORBX)
      CALL DZERO(A12,N2ORBX)
      CALL DZERO(A123,N2ORBX)
C
      CALL MEMGET('REAL',KCREF,MZCONF(1),WRK,KFREE,LFREE)
      CALL GETREF(WRK(KCREF),MZCONF(1))
C
C
C A = 2*sum(l,m) f(l) <0| TE(l,m) |0> TE(l,m)
C
      CALL TOPGET(DUMMY,DUMMY,DUMMY,0,D1,1,
     *            IDUM,IDUM,IDUM,A,UDV,EPSOL,.FALSE.,
     *            ISYRLM,CMO,WRK(KFREE),LFREE,IPRLVL,'A',0)
c     call TOPGET(ZYM1,ZYM2,ZYM3,IKLVL,OVLAP,ISYMDN,
c    *            ISYMV1,ISYMV2,ISYMV3,TOP,DEN1,EPS,NUCFLG,
c    *            ISYRLM,CMO,WRK,LWRK,IPRLVL,PRSTR,ADDFLG)
      CALL DSCAL(N2ORBX,D2,A,1)
C
C W = 2*sum(l,m) f(l)*( <0|TE(l,m)|0> - Tn(l,m) )*TE(l,m) 
C
      IF (INERSI) THEN
         EPS = EPSTAT
      ELSE
         EPS = EPSOL
      END IF
      CALL TOPGET(DUMMY,DUMMY,DUMMY,0,DUMMY,1,
     *            IDUM,IDUM,IDUM,W,DUMMY,EPS,.TRUE.,
     *            ISYRLM,CMO,WRK(KFREE),LFREE,IPRLVL,'W',0)
      CALL DSCAL(N2ORBX,-D1,W,1)
      CALL TOPGET(DUMMY,DUMMY,DUMMY,0,D1,1,
     *            IDUM,IDUM,IDUM,W,UDV,EPS,.FALSE.,
     *            ISYRLM,CMO,WRK(KFREE),LFREE,IPRLVL,'W',0)
      CALL DSCAL(N2ORBX,D2,W,1)
C
C A1 = 2*sum(l,m) f(l) <0| TE(l,m)(k1) |0> TE(l,m) + ...
C
      IF (MZWOPT(ISYMV1).GT.0) THEN
         CALL TOPGET(ZYM1,DUMMY,DUMMY,1,D1,1,
     *        ISYMV1,IDUM,IDUM,A1,UDV,EPSOL,.FALSE.,
     *        ISYRLM,CMO,WRK(KFREE),LFREE,IPRLVL,'A1 cont1',0)
      END IF
C
C ... + 2*sum(l,m) f(l) ( <01L| TE(l,m) |0> + <0| TE(l,m) |01R> ) TE(l,m) 
C
      IF (MZCONF(ISYMV1).GT.0) THEN
C
C     Construct the density matrix <01L|..|0> + <0|..|01R>
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
     *        WRK(KCREF),VEC1,OVLAP,DEN1,DUMMY,JSPIN,JSPIN,TDM,
     *        NORHO2,XINDX,WRK,KFREE,LFREE,LREF)
C
         CALL TOPGET(DUMMY,DUMMY,DUMMY,0,OVLAP,ISYMDN,
     *        IDUM,IDUM,IDUM,B1,DEN1,EPSOL,.FALSE.,
     *        ISYRLM,CMO,WRK(KFREE),LFREE,IPRLVL,'A1 cont2',0)
      END IF
C
      CALL DAXPY(N2ORBX,D1,B1,1,A1,1)
      CALL DSCAL(N2ORBX,D2,A1,1)
C
C B1 = -2/3*sum(l,m) f(l) ( <01L| TE(l,m) |0> + <0| TE(l,m) |01R> ) TE(l,m) 
C
      CALL DSCAL(N2ORBX,-D2/D3,B1,1)
C
C
C A2 = 2*sum(l,m) f(l) <0| TE(l,m)(k2) |0> TE(l,m) + ...
C
      IF (MZWOPT(ISYMV2).GT.0) THEN
         CALL TOPGET(ZYM2,DUMMY,DUMMY,1,D1,1,
     *        ISYMV2,IDUM,IDUM,A2,UDV,EPSOL,.FALSE.,
     *        ISYRLM,CMO,WRK(KFREE),LFREE,IPRLVL,'A2 cont1',0)
      END IF
C
C ... + 2*sum(l,m) f(l) ( <02L| TE(l,m) |0> + <0| TE(l,m) |02R> ) TE(l,m) 
C
      IF (MZCONF(ISYMV2).GT.0) THEN
C
C     Construct the density matrix <02L|..|0> + <0|..|02R>
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
     *        WRK(KCREF),VEC2,OVLAP,DEN1,DUMMY,JSPIN,JSPIN,TDM,
     *        NORHO2,XINDX,WRK,KFREE,LFREE,LREF)
C
         CALL TOPGET(DUMMY,DUMMY,DUMMY,0,OVLAP,ISYMDN,
     *        IDUM,IDUM,IDUM,A2,DEN1,EPSOL,.FALSE.,
     *        ISYRLM,CMO,WRK(KFREE),LFREE,IPRLVL,'A2 cont2',0)
      END IF
C
      CALL DSCAL(N2ORBX,D2,A2,1)
C
C
C A12 = sum(l,m) f(l) <0| TE(l,m)(k1,k2) |0> TE(l,m)
C     + sum(l,m) f(l) <0| TE(l,m)(k2,k1) |0> TE(l,m) + ...
C
      IF (MZWOPT(ISYMV1).GT.0 .AND. MZWOPT(ISYMV2).GT.0) THEN
         CALL TOPGET(ZYM1,ZYM2,DUMMY,2,D1,1,
     *               ISYMV1,ISYMV2,IDUM,A12,UDV,EPSOL,.FALSE.,
     *               ISYRLM,CMO,WRK(KFREE),LFREE,IPRLVL,'A12 cont1a',0)
         CALL TOPGET(ZYM2,ZYM1,DUMMY,2,D1,1,
     *               ISYMV2,ISYMV1,IDUM,A12,UDV,EPSOL,.FALSE.,
     *               ISYRLM,CMO,WRK(KFREE),LFREE,IPRLVL,'A12 cont1b',0)
      END IF
C
C ... + 2*sum(l,m) f(l) ( <01L| TE(l,m)(k2) |0> + 
C                           <0| TE(l,m)(k2) |01R> ) TE(l,m) + ...
C
C     Put the factor two into one of the vectors.
C
      CALL DSCAL(KZYV1,D2,VEC1,1)
      CALL DSCAL(NORBT*NORBT,D2,ZYM1,1)
C
      IF (MZCONF(ISYMV1).GT.0 .AND. MZWOPT(ISYMV2).GT.0) THEN
C
C     Construct the density matrix <01L|..|0> + <0|..|01R>
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
     *         WRK(KCREF),VEC1,OVLAP,DEN1,DUMMY,JSPIN,JSPIN,TDM,
     *         NORHO2,XINDX,WRK,KFREE,LFREE,LREF)
C
         CALL TOPGET(ZYM2,DUMMY,DUMMY,1,OVLAP,ISYMDN,
     *               ISYMV2,IDUM,IDUM,A12,DEN1,EPSOL,.FALSE.,
     *               ISYRLM,CMO,WRK(KFREE),LFREE,IPRLVL,'A12 cont2a',0)
      END IF
C
C ... + 2*sum(l,m) f(l) ( <02L| TE(l,m)(k1) |0> + 
C                           <0| TE(l,m)(k1) |02R> ) TE(l,m) + ...
C
C     The factor two is already included in one of the vectors.
C
      IF (MZCONF(ISYMV2).GT.0 .AND. MZWOPT(ISYMV1).GT.0) THEN
C
C     Construct the density matrix <02L|..|0> + <0|..|02R>
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
     *         WRK(KCREF),VEC2,OVLAP,DEN1,DUMMY,JSPIN,JSPIN,TDM,
     *         NORHO2,XINDX,WRK,KFREE,LFREE,LREF)
C
         CALL TOPGET(ZYM1,DUMMY,DUMMY,1,OVLAP,ISYMDN,
     *               ISYMV1,IDUM,IDUM,A12,DEN1,EPSOL,.FALSE.,
     *               ISYRLM,CMO,WRK(KFREE),LFREE,IPRLVL,'A12 cont2b',0)
      END IF
C     
C ... + 2*sum(l,m) f(l) ( <01L| TE(l,m) |02R> + 
C                         <02L| TE(l,m) |01R> ) TE(l,m) + ...
C
C     The factor two is already included in one of the vectors.
C
      IF (MZCONF(ISYMV1) .GT. 0 .AND. MZCONF(ISYMV2) .GT. 0) THEN
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
     *         VEC1,VEC2,OVLAP,DEN1,DUMMY,JSPIN,JSPIN,TDM,NORHO2,
     *         XINDX,WRK,KFREE,LFREE,LREF)
C
         CALL TOPGET(DUMMY,DUMMY,DUMMY,0,OVLAP,ISYMDN,
     *               IDUM,IDUM,IDUM,A12,DEN1,EPSOL,.FALSE.,
     *               ISYRLM,CMO,WRK(KFREE),LFREE,IPRLVL,'A12 cont3',0)
      END IF
C
C     Restore the scaled vector
C
      CALL DSCAL(KZYV1,DH,VEC1,1)
      CALL DSCAL(NORBT*NORBT,DH,ZYM1,1)
C
C
C A123 = sum(l,m) f(l) <0| 1/3*TE(l,m)(k1,k2,k3) |0> TE(l,m)
C      + sum(l,m) f(l) <0| 1/3*TE(l,m)(k1,k3,k2) |0> TE(l,m) + ...
C
      IF (MZWOPT(ISYMV1).GT.0 .AND. MZWOPT(ISYMV2).GT.0
     *    .AND. MZWOPT(ISYMV3).GT.0) THEN
         CALL TOPGET(ZYM1,ZYM2,ZYM3,3,D1,1,
     *               ISYMV1,ISYMV2,ISYMV3,A123,UDV,EPSOL,.FALSE.,
     *               ISYRLM,CMO,WRK(KFREE),LFREE,IPRLVL,'A123 cont1a',0)
         CALL TOPGET(ZYM1,ZYM3,ZYM2,3,D1,1,
     *               ISYMV1,ISYMV3,ISYMV2,A123,UDV,EPSOL,.FALSE.,
     *               ISYRLM,CMO,WRK(KFREE),LFREE,IPRLVL,'A123 cont1b',0)
         CALL DSCAL(N2ORBX,D1/D3,A123,1)
      END IF
C
C ... + sum(l,m) f(l) ( <01L| TE(l,m)(k2,k3) |0> + 
C                         <0| TE(l,m)(k2,k3) |01R> ) TE(l,m) + ...
C
C
      IF (MZCONF(ISYMV1).GT.0 .AND. MZWOPT(ISYMV2).GT.0
     *    .AND. MZWOPT(ISYMV3).GT.0) THEN
C
C     Construct the density matrix <01L|..|0> + <0|..|01R>
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
     *         WRK(KCREF),VEC1,OVLAP,DEN1,DUMMY,JSPIN,JSPIN,TDM,
     *         NORHO2,XINDX,WRK,KFREE,LFREE,LREF)
C
         CALL TOPGET(ZYM2,ZYM3,DUMMY,2,OVLAP,ISYMDN,
     *               ISYMV2,ISYMV3,IDUM,A123,DEN1,EPSOL,.FALSE.,
     *               ISYRLM,CMO,WRK(KFREE),LFREE,IPRLVL,'A123 cont2a',0)
C
C ... + sum(l,m) f(l) ( <01L| TE(l,m)(k3,k2) |0> + 
C                         <0| TE(l,m)(k3,k2) |01R> ) TE(l,m) + ...
C
         CALL TOPGET(ZYM3,ZYM2,DUMMY,2,OVLAP,ISYMDN,
     *               ISYMV3,ISYMV2,IDUM,A123,DEN1,EPSOL,.FALSE.,
     *               ISYRLM,CMO,WRK(KFREE),LFREE,IPRLVL,'A123 cont2b',0)
      END IF
C     
C ... + 2*sum(l,m) f(l) ( <02L| TE(l,m)(k1) |03R> + 
C                         <03L| TE(l,m)(k1) |02R> ) TE(l,m) + ...
C
      IF (MZCONF(ISYMV2) .GT. 0 .AND. MZCONF(ISYMV3) .GT. 0
     *    .AND. MZWOPT(ISYMV1).GT.0) THEN
C
C     Construct <02L|..|03R> + <03L|..|02R> density
C
         ILSYM  = MULD2H(IREFSY,ISYMV2)
         IRSYM  = MULD2H(IREFSY,ISYMV3)
         NCL    = MZCONF(ISYMV2)
         NCR    = MZCONF(ISYMV3)
         KZVARL = MZYVAR(ISYMV2)
         KZVARR = MZYVAR(ISYMV3)
         LREF   = .FALSE.
         ISYMDN = MULD2H(ILSYM,IRSYM)
         CALL DZERO(DEN1,NASHT*NASHT)
         CALL RSPGDM(NSIM,ILSYM,IRSYM,NCL,NCR,KZVARL,KZVARR,
     *         VEC2,VEC3,OVLAP,DEN1,DUMMY,JSPIN,JSPIN,TDM,NORHO2,
     *         XINDX,WRK,KFREE,LFREE,LREF)
C
         CALL TOPGET(ZYM1,DUMMY,DUMMY,1,OVLAP,ISYMDN,
     *               ISYMV1,IDUM,IDUM,A123,DEN1,EPSOL,.FALSE.,
     *               ISYRLM,CMO,WRK(KFREE),LFREE,IPRLVL,'A123 cont3',0)
         CALL TOPGET(ZYM1,DUMMY,DUMMY,1,OVLAP,ISYMDN,
     *               ISYMV1,IDUM,IDUM,A123,DEN1,EPSOL,.FALSE.,
     *               ISYRLM,CMO,WRK(KFREE),LFREE,IPRLVL,'A123 cont3',0)
      END IF
C
      RETURN
      END
      SUBROUTINE WCASE1(VECA,VEC1,VEC2,VEC3,ETRS,XINDX,ISYRLM,
     *            W,A,A1,A2,
     *            A12,
     *            DEN1,UDV,WRK,LFREE,ZYM1,ZYM2,
     *            KZYVR,KZYV1,KZYV2,KZYV3,IGRSYM,ISYMV1,ISYMV2,ISYMV3)
C
#include <implicit.h>
C
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0,
     *            D4 = 4.0D0, DH = 0.5D0 )
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
#include <infinp.h>
#include <dummy.h>
C
      DIMENSION UDV(NASHDI,NASHDI),DEN1(NASHDI,NASHDI)
      DIMENSION ZYM1(*),ZYM2(*)
      DIMENSION VEC1(KZYV1),VEC2(KZYV2),VEC3(KZYV3),ETRS(KZYVR)
      DIMENSION XINDX(*),ISYRLM(*),WRK(*),VECA(KZYVR)
      DIMENSION W(N2ORBX),A(N2ORBX),A1(N2ORBX),A2(N2ORBX)
      DIMENSION A12(N2ORBX)
C
      LOGICAL   TDM, LREF, NORHO2
C
C     Initialise variables
C
      JSPIN  = 0
      TDM    = .TRUE.
      KFREE = 1
      IPRONE = 100
      NORHO2 = .TRUE.
      NSIM = 1
      F1 = D0
      F1L = D0
      F1R = D0
      F2L = D0
      F2R = D0
C
      CALL MEMGET('REAL',KCREF,MZCONF(1),WRK,KFREE,LFREE)
      CALL GETREF(WRK(KCREF),MZCONF(1))
C
      CALL MEMGET('REAL',KTMP,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KWT,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KRES,KZYVR ,WRK,KFREE,LFREE)
      CALL DZERO(WRK(KTMP),N2ORBX)
C
      IF (MZCONF(ISYMV3).LE.0 .OR. IGRSYM.NE.ISYMV3) RETURN
C
C     /   0              \
C     | (F1+F2) * Sj(3)  |
C     |   0              |
C     \ (F1+F3) * Sj(3)' /
C
      IF (MZCONF(ISYMV1) .GT. 0 .AND. MZCONF(ISYMV2) .GT. 0) THEN
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
     *         VEC1,VEC2,OVLAP,DEN1,DUMMY,JSPIN,JSPIN,TDM,NORHO2,
     *         XINDX,WRK,KFREE,LFREE,LREF)
C
C        F1 = 1/2*<01L|W|02R> + perm. + ...
C
         CALL MELONE(W,1,DEN1,OVLAP,F1,IPRONE,'F1 in WCASE1')
         F1 = F1*DH
C
C        ... + FACT*<0|2/3*W+A|0>
C
         CALL DCOPY(N2ORBX,W,1,WRK(KTMP),1)
         CALL DSCAL(N2ORBX,D2/D3,WRK(KTMP),1)
         CALL DAXPY(N2ORBX,D1,A,1,WRK(KTMP),1)
         FACT = DDOT(MZCONF(ISYMV1),VEC1,1,VEC2(MZVAR(ISYMV2)+1),1) +
     *          DDOT(MZCONF(ISYMV2),VEC2,1,VEC1(MZVAR(ISYMV1)+1),1)
         CALL DSCAL(N2ORBX,FACT,WRK(KTMP),1)
      END IF
C
C        ... + <0|1/2*WT12+A12+A1T2|0> + perm. 
C
      IF (MZWOPT(ISYMV1).GT.0 .AND. MZWOPT(ISYMV2).GT.0) THEN
         CALL DZERO(WRK(KWT),N2ORBX)
         CALL OITH1(ISYMV1,ZYM1,W,WRK(KWT),1)
         CALL DSCAL(N2ORBX,DH,WRK(KWT),1)
         CALL OITH1(ISYMV2,ZYM2,WRK(KWT),WRK(KTMP),ISYMV1)
         CALL DZERO(WRK(KWT),N2ORBX)
         CALL OITH1(ISYMV2,ZYM2,W,WRK(KWT),1)
         CALL DSCAL(N2ORBX,DH,WRK(KWT),1)
         CALL OITH1(ISYMV1,ZYM1,WRK(KWT),WRK(KTMP),ISYMV2)
      END IF
C
      IF (MZWOPT(ISYMV1).GT.0)
     *     CALL OITH1(ISYMV1,ZYM1,A2,WRK(KTMP),ISYMV2)
      IF (MZWOPT(ISYMV2).GT.0)
     *     CALL OITH1(ISYMV2,ZYM2,A1,WRK(KTMP),ISYMV1)
C     
      CALL DAXPY(N2ORBX,D1,A12,1,WRK(KTMP),1)
C     
      CALL MELONE(WRK(KTMP),MULD2H(ISYMV1,ISYMV2),UDV,D1,FACT,
     *            IPRONE,'FACT in WCASE1')
      F1 = F1 + FACT
C
      IF (MZCONF(ISYMV1).GT.0) THEN
         CALL DCOPY(N2ORBX,A2,1,WRK(KTMP),1)
         IF (MZWOPT(ISYMV2).GT.0)
     *      CALL OITH1(ISYMV2,ZYM2,W,WRK(KTMP),1)
C
C     F1R=<0|WT2 + A2|-01R>
C
         ILSYM  = IREFSY
         IRSYM  = MULD2H(IREFSY,ISYMV1)
         NCL    = MZCONF(1)
         NCR    = MZCONF(ISYMV1)
         IOFF   = 1
         ISYMDN = MULD2H(ILSYM,IRSYM)
         CALL DZERO(DEN1,NASHT*NASHT)
         CALL RSPDM(ILSYM,IRSYM,NCL,NCR,WRK(KCREF),VEC1(IOFF),
     *        DEN1,DUMMY,JSPIN,JSPIN,TDM,NORHO2,XINDX,WRK,
     *        KFREE,LFREE)
C
         OVLAP = D0
         IF (ILSYM.EQ.IRSYM)
     *        OVLAP = DDOT(NCL,WRK(KCREF),1,VEC1(IOFF),1)
         CALL MELONE(WRK(KTMP),ISYMV2,DEN1,OVLAP,F1R,IPRONE,
     *               'F1R in WCASE1')
C
C     F1L=<01L|WT2 + A2|0>
C
         ILSYM  = MULD2H(IREFSY,ISYMV1)
         IRSYM  = IREFSY
         NCL    = MZCONF(ISYMV1)
         NCR    = MZCONF(1)
         IOFF   = MZVAR(ISYMV1) + 1
         ISYMDN = MULD2H(ILSYM,IRSYM)
         CALL DZERO(DEN1,NASHT*NASHT)
         CALL RSPDM(ILSYM,IRSYM,NCL,NCR,VEC1(IOFF),WRK(KCREF),
     *        DEN1,DUMMY,JSPIN,JSPIN,TDM,NORHO2,XINDX,WRK,
     *        KFREE,LFREE)
C
         OVLAP = D0
         IF (ILSYM.EQ.IRSYM)
     *        OVLAP = DDOT(NCL,WRK(KCREF),1,VEC1(IOFF),1)
         CALL MELONE(WRK(KTMP),ISYMV2,DEN1,OVLAP,F1L,IPRONE,
     *               'F1L in WCASE1')
      END IF
C
      IF (MZCONF(ISYMV2).GT.0) THEN
         CALL DCOPY(N2ORBX,A1,1,WRK(KTMP),1)
         IF (MZWOPT(ISYMV1).GT.0)
     *      CALL OITH1(ISYMV1,ZYM1,W,WRK(KTMP),1)
C
C     F2R=<0|WT1 + A1|-02R>
C
         ILSYM  = IREFSY
         IRSYM  = MULD2H(IREFSY,ISYMV2)
         NCL    = MZCONF(1)
         NCR    = MZCONF(ISYMV2)
         IOFF   = 1
         ISYMDN = MULD2H(ILSYM,IRSYM)
         CALL DZERO(DEN1,NASHT*NASHT)
         CALL RSPDM(ILSYM,IRSYM,NCL,NCR,WRK(KCREF),VEC2(IOFF),
     *        DEN1,DUMMY,JSPIN,JSPIN,TDM,NORHO2,XINDX,WRK,
     *        KFREE,LFREE)
C
         OVLAP = D0
         IF (ILSYM.EQ.IRSYM)
     *        OVLAP = DDOT(NCL,WRK(KCREF),1,VEC2(IOFF),1)
         CALL MELONE(WRK(KTMP),ISYMV1,DEN1,OVLAP,F2R,IPRONE,
     *               'F2R in WCASE1')
C
C     F2L=<02L|WT1 + A1|0>
C
         ILSYM  = MULD2H(IREFSY,ISYMV2)
         IRSYM  = IREFSY
         NCL    = MZCONF(ISYMV2)
         NCR    = MZCONF(1)
         IOFF   = MZVAR(ISYMV2) + 1
         ISYMDN = MULD2H(ILSYM,IRSYM)
         CALL DZERO(DEN1,NASHT*NASHT)
         CALL RSPDM(ILSYM,IRSYM,NCL,NCR,VEC2(IOFF),WRK(KCREF),
     *        DEN1,DUMMY,JSPIN,JSPIN,TDM,NORHO2,XINDX,WRK,
     *        KFREE,LFREE)
C
         OVLAP = D0
         IF (ILSYM.EQ.IRSYM)
     *        OVLAP = DDOT(NCL,WRK(KCREF),1,VEC2(IOFF),1)
         CALL MELONE(WRK(KTMP),ISYMV1,DEN1,OVLAP,F2L,IPRONE,
     *               'F2L in WCASE1')
      END IF
C
      NZCONF = MZCONF(IGRSYM)
      NZVAR  = MZVAR(IGRSYM)
      FACT = F1 + DH*(F1L + F2L) - (F1R + F2R)
      CALL DAXPY(NZCONF,FACT,VEC3,1,ETRS,1)
      FACT = F1 + (F1L + F2L) - DH*(F1R + F2R)
      CALL DAXPY(NZCONF,FACT,VEC3(NZVAR+1),1,ETRS(NZVAR+1),1)
C
      CALL PRIRES(ETRS,VECA,IGRSYM,'WCASE1')
C     
      RETURN
      END
      SUBROUTINE WCASE2(VECA,VEC1,VEC2,VEC3,ETRS,XINDX,ISYRLM,
     *            W,A,A1,A2,A12,
     *            DEN1,UDV,WRK,LFREE,CMO,MJWOP,ZYM1,ZYM2,
     *            KZYVR,KZYV1,KZYV2,KZYV3,IGRSYM,ISYMV1,ISYMV2,ISYMV3)
C
#include <implicit.h>
#include <dummy.h>
C
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0,
     *            D4 = 4.0D0, DH = 0.5D0 )
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
#include <infinp.h>
C
      DIMENSION UDV(NASHDI,NASHDI),DEN1(NASHDI,NASHDI)
      DIMENSION ZYM1(*),ZYM2(*)
      DIMENSION VEC1(KZYV1),VEC2(KZYV2),VEC3(KZYV3),ETRS(KZYVR)
      DIMENSION XINDX(*),ISYRLM(*),WRK(*),CMO(*),MJWOP(*),VECA(KZYVR)
      DIMENSION W(N2ORBX),A(N2ORBX),A1(N2ORBX),A2(N2ORBX),A12(N2ORBX)
C
      LOGICAL LCON, LORB
      LOGICAL   TDM, LREF, NORHO2
C
C     Initialise variables
C
      JSPIN  = 0
      TDM    = .TRUE.
      KFREE = 1
      IPRONE = 100
      NORHO2 = .TRUE.
      NSIM = 1
C
      CALL MEMGET('REAL',KCREF,MZCONF(1),WRK,KFREE,LFREE)
      CALL GETREF(WRK(KCREF),MZCONF(1))
C
      CALL MEMGET('REAL',KTMP,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KWT,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KRES,KZYVR,WRK,KFREE,LFREE)
C
      IF (MZCONF(ISYMV3).LE.0) RETURN
C
C     Construct the density matrix <03L|..|0> + <0|..|03R>
C
      ILSYM  = IREFSY
      IRSYM  = MULD2H(IREFSY,ISYMV3)
      NCL    = MZCONF(1)
      NCR    = MZCONF(ISYMV3)
      KZVARL = MZCONF(1)
      KZVARR = MZYVAR(ISYMV3)
      LREF   = .TRUE.
      ISYMDN = MULD2H(ILSYM,IRSYM)
      CALL DZERO(DEN1,NASHT*NASHT)
      CALL RSPGDM(NSIM,ILSYM,IRSYM,NCL,NCR,KZVARL,KZVARR,
     *         WRK(KCREF),VEC3,OVLAP,DEN1,DUMMY,JSPIN,JSPIN,TDM,
     *         NORHO2,XINDX,WRK,KFREE,LFREE,LREF)
C
C     /   <0| [qj,T] |03R>  + <03L| [qj,T] |0>  \
C     |   <j| T |03R>                           |
C     |   <0| [qj+,T] |03R> + <03L| [qj+,T] |0> |
C     \  -<03L| T |j>                           /
C
C     T = 1/2*WT12 + A12 + A1T2
C
         CALL DZERO(WRK(KTMP),N2ORBX)
         CALL DZERO(WRK(KWT),N2ORBX)
         IF (MZWOPT(ISYMV1).GT.0)
     *      CALL OITH1(ISYMV1,ZYM1,W,WRK(KWT),1)
         CALL DSCAL(N2ORBX,DH,WRK(KWT),1)
         IF (MZWOPT(ISYMV2).GT.0)
     *      CALL OITH1(ISYMV2,ZYM2,WRK(KWT),WRK(KTMP),ISYMV1)
C
         CALL DZERO(WRK(KWT),N2ORBX)
         IF (MZWOPT(ISYMV2).GT.0)
     *      CALL OITH1(ISYMV2,ZYM2,W,WRK(KWT),1)
         CALL DSCAL(N2ORBX,DH,WRK(KWT),1)
         IF (MZWOPT(ISYMV1).GT.0)
     *      CALL OITH1(ISYMV1,ZYM1,WRK(KWT),WRK(KTMP),ISYMV2)
C
         CALL DAXPY(N2ORBX,D1,A12,1,WRK(KTMP),1)
         IF (MZWOPT(ISYMV1).GT.0)
     *      CALL OITH1(ISYMV1,ZYM1,A2,WRK(KTMP),ISYMV2)
         IF (MZWOPT(ISYMV2).GT.0)
     *      CALL OITH1(ISYMV2,ZYM2,A1,WRK(KTMP),ISYMV1)
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
      NZYVEC = MZYVAR(ISYMV3)
      NZCVEC = MZCONF(ISYMV3)
      CALL RSP1GR(NSIM,KZYVR,IDUM,JSPIN,IGRSYM,JSPIN,ISYMV3,ETRS,
     *            VEC3,NZYVEC,NZCVEC,OVLAP,ISYMDN,DEN1,WRK(KTMP),
     *            XINDX,MJWOP,WRK(KFREE),LFREE,LORB,LCON,LREF)
C
      CALL PRIRES(ETRS,VECA,IGRSYM,'WCASE2a')
C
      IF (MZCONF(ISYMV1) .LE. 0 .OR. MZCONF(ISYMV2) .LE. 0
     *    .OR. IGRSYM .NE. ISYMV3) RETURN
C
C
C     /   <0| [qj, 2/3*W+A ] |03R>  + <03L| [qj, 2/3*W+A ] |0>  \
C     |   1/4*<j|  2/3*W+A |03R>                                |*FACT
C     |   <0| [qj+, 2/3*W+A ] |03R> + <03L| [qj+, 2/3*W+A ] |0> |
C     \  -1/4*<03L|  2/3*W+A |j>                                /
C
C     FACT = S(1)S(2)' + S(1)'S(2)
C
      CALL DCOPY(N2ORBX,W,1,WRK(KTMP),1)
      CALL DSCAL(N2ORBX,D2/D3,WRK(KTMP),1)
      CALL DAXPY(N2ORBX,D1,A,1,WRK(KTMP),1)
C
C     Make the gradient (we already have the density and variables)
C
      CALL DZERO(WRK(KRES),KZYVR)
      CALL RSP1GR(NSIM,KZYVR,IDUM,JSPIN,IGRSYM,JSPIN,ISYMV3,WRK(KRES),
     *            VEC3,NZYVEC,NZCVEC,OVLAP,ISYMDN,DEN1,WRK(KTMP),
     *            XINDX,MJWOP,WRK(KFREE),LFREE,LORB,LCON,LREF)
C
      FACT = DDOT(MZCONF(ISYMV1),VEC1,1,VEC2(MZVAR(ISYMV2)+1),1) +
     *       DDOT(MZCONF(ISYMV2),VEC2,1,VEC1(MZVAR(ISYMV1)+1),1)
      NZCONF = MZCONF(IGRSYM)
      NZWOPT = MZWOPT(IGRSYM)
      IOFFZC = 0
      IOFFZO = IOFFZC + NZCONF
      IOFFYC = IOFFZO + NZWOPT
      IOFFYO = IOFFYC + NZCONF
      CALL DAXPY(NZCONF,FACT/D4,WRK(KRES+IOFFZC),1,ETRS(1+IOFFZC),1)
      CALL DAXPY(NZCONF,FACT/D4,WRK(KRES+IOFFYC),1,ETRS(1+IOFFYC),1)
C
      CALL DAXPY(NZWOPT,FACT,WRK(KRES+IOFFZO),1,ETRS(1+IOFFZO),1)
      CALL DAXPY(NZWOPT,FACT,WRK(KRES+IOFFYO),1,ETRS(1+IOFFYO),1)
C
      CALL PRIRES(ETRS,VECA,IGRSYM,'WCASE2b')
C
C     /       0       \
C     |  <j| A |03R>  | * ( S(1)S(2)' + S(1)'S(2) ) * 3/4
C     |       0       |
C     \ -<03L| A |j>  /
C
C     Make the gradient (we already have the density and variables)
C
      IF (LCON) THEN
         CALL DZERO(WRK(KRES),KZYVR)
         CALL CONSX(NSIM,KZYVR,IGRSYM,A,VEC3,NZYVEC,NZCVEC,
     *              ISYMV3,MZCONF(IGRSYM),ISYMST,LREF,WRK(KRES),
     *              XINDX,WRK(KFREE),LFREE)
         CALL DAXPY(KZYVR,FACT*D3/D4,WRK(KRES),1,ETRS,1)
      END IF
C
      CALL PRIRES(ETRS,VECA,IGRSYM,'WCASE2c')
C     
      RETURN
      END
      SUBROUTINE WCASE3(VECA,VEC1,VEC2,VEC3,ETRS,XINDX,ISYRLM,
     *            W,A1,
     *            DEN1,UDV,WRK,LFREE,CMO,MJWOP,ZYM1,
     *            KZYVR,KZYV1,KZYV2,KZYV3,IGRSYM,ISYMV1,ISYMV2,ISYMV3)
C
#include <implicit.h>
C
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0,
     *            D4 = 4.0D0, DH = 0.5D0 )
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
#include <infinp.h>
#include <dummy.h>
C
      DIMENSION UDV(NASHDI,NASHDI),DEN1(NASHDI,NASHDI)
      DIMENSION VEC1(KZYV1),VEC2(KZYV2),VEC3(KZYV3),ETRS(KZYVR)
      DIMENSION XINDX(*),ISYRLM(*),WRK(*),CMO(*),MJWOP(*),VECA(KZYVR)
      DIMENSION W(N2ORBX),A1(N2ORBX),ZYM1(*)
C
      LOGICAL   TDM, LREF, NORHO2
C
C     Initialise variables
C
      JSPIN  = 0
      TDM    = .TRUE.
      KFREE = 1
      IPRONE = 100
      NORHO2 = .TRUE.
      NSIM = 1
C
      CALL MEMGET('REAL',KCREF,MZCONF(1),WRK,KFREE,LFREE)
      CALL GETREF(WRK(KCREF),MZCONF(1))
C
      CALL MEMGET('REAL',KTMP,N2ORBX,WRK,KFREE,LFREE)
C
      CALL DCOPY(N2ORBX,A1,1,WRK(KTMP),1)
      CALL OITH1(ISYMV1,ZYM1,W,WRK(KTMP),1)
C
      IF (MZCONF(ISYMV2) .EQ. 0 .OR. MZCONF(ISYMV3) .EQ. 0) RETURN
C
C     /   <02L| [qj,WT1+A1] |03R>  + <03L| [qj,WT1+A1] |02R>  \
C     |                       0                               |
C     |   <02L| [qj+,WT1+A1] |03R> + <03L| [qj+,WT1+A1] |02R> |
C     \                       0                               /
C
C     Construct <02L|..|03R> + <03L|..|02R> density
C
      ILSYM  = MULD2H(IREFSY,ISYMV2)
      IRSYM  = MULD2H(IREFSY,ISYMV3)
      NCL    = MZCONF(ISYMV2)
      NCR    = MZCONF(ISYMV3)
      KZVARL = MZYVAR(ISYMV2)
      KZVARR = MZYVAR(ISYMV3)
      LREF   = .FALSE.
      ISYMDN = MULD2H(ILSYM,IRSYM)
      CALL DZERO(DEN1,NASHT*NASHT)
      CALL RSPGDM(NSIM,ILSYM,IRSYM,NCL,NCR,KZVARL,KZVARR,
     *         VEC2,VEC3,OVLAP,DEN1,DUMMY,JSPIN,JSPIN,TDM,NORHO2,
     *         XINDX,WRK,KFREE,LFREE,LREF)
C
C     Make the gradient
C
      IF ( MZWOPT(IGRSYM) .GT. 0 ) THEN
         CALL ORBSX(NSIM,IGRSYM,KZYVR,ETRS,WRK(KTMP),OVLAP,ISYMDN,
     *              DEN1,MJWOP,WRK(KFREE),LFREE) 
      END IF
C
      CALL PRIRES(ETRS,VECA,IGRSYM,'WCASE3')
C     
      RETURN
      END
      SUBROUTINE WCASE4(VECA,VEC1,VEC2,VEC3,ETRS,XINDX,ISYRLM,
     *            W,A,A1,B1,
     *            A12,A123,
     *            DEN1,UDV,WRK,LFREE,CMO,MJWOP,ZYM1,ZYM2,ZYM3,
     *            KZYVR,KZYV1,KZYV2,KZYV3,IGRSYM,ISYMV1,ISYMV2,ISYMV3)
C
#include <implicit.h>
C
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0,
     *            D4 = 4.0D0, D6 = 6.0D0, DH = 0.5D0 )
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
#include <infinp.h>
#include <dummy.h>
C
      DIMENSION UDV(NASHDI,NASHDI),DEN1(NASHDI,NASHDI)
      DIMENSION ZYM1(*),ZYM2(*),ZYM3(*)
      DIMENSION VEC1(KZYV1),VEC2(KZYV2),VEC3(KZYV3),ETRS(KZYVR)
      DIMENSION XINDX(*),ISYRLM(*),WRK(*),CMO(*),MJWOP(*),VECA(KZYVR)
      DIMENSION W(N2ORBX),A(N2ORBX),A1(N2ORBX)
      DIMENSION B1(N2ORBX)
      DIMENSION A12(N2ORBX),A123(N2ORBX)
C
      LOGICAL   TDM, LREF, LCON, LORB, NORHO2
C
C     Initialise variables
C
      JSPIN  = 0
      TDM    = .TRUE.
      KFREE = 1
      IPRONE = 100
      NORHO2 = .TRUE.
      NSIM = 1
      DATA F1/D0/,F1L/D0/,F1R/D0/,F2L/D0/,F2R/D0/
C
      CALL MEMGET('REAL',KCREF,MZCONF(1),WRK,KFREE,LFREE)
      CALL GETREF(WRK(KCREF),MZCONF(1))
C
      CALL MEMGET('REAL',KTMP,N2ORBX,WRK,KFREE,LFREE)
      CALL DZERO(WRK(KTMP),N2ORBX)
      CALL MEMGET('REAL',KRES,KZYVR ,WRK,KFREE,LFREE)
C
      IF (MZCONF(ISYMV2) .EQ. 0 .OR. MZCONF(ISYMV3) .EQ. 0
     *     .OR. ISYMV2.NE.ISYMV3) GOTO 100
C
C     /       0          \
C     |  <j| WT1+A1 |0>  | * ( S(2)S(3)' + S(2)'S(3) ) * -1/2
C     |       0          |
C     \ -<0| WT1+A1 |j>  /
C
C     FACT = S(2)S(3)' + S(2)'S(3)
C
      FACT = DDOT(MZCONF(ISYMV2),VEC2,1,VEC3(MZVAR(ISYMV3)+1),1) +
     *       DDOT(MZCONF(ISYMV3),VEC3,1,VEC2(MZVAR(ISYMV2)+1),1)
C
      CALL OITH1(ISYMV1,ZYM1,W,WRK(KTMP),1)
C
      CALL DAXPY(N2ORBX,D1,A1,1,WRK(KTMP),1)
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
         CALL CONSX(NSIM,KZYVR,IGRSYM,WRK(KTMP),WRK(KCREF),
     *              MZCONF(1),MZCONF(1),IREFSY,MZCONF(IGRSYM),ISYMST,
     *              LREF,WRK(KRES),XINDX,WRK(KFREE),LFREE)
         CALL DAXPY(KZYVR,-FACT*DH,WRK(KRES),1,ETRS,1)
      END IF
C
      CALL OITH1(ISYMV1,ZYM1,A,WRK(KTMP),1)
      CALL DAXPY(N2ORBX,D1,A1,1,WRK(KTMP),1)
      CALL DAXPY(N2ORBX,D1,B1,1,WRK(KTMP),1)
      CALL DSCAL(N2ORBX,FACT,WRK(KTMP),1)
C
 100  CONTINUE
C
      CALL DAXPY(N2ORBX,D1,A123,1,WRK(KTMP),1)
      CALL OITH1(ISYMV3,ZYM3,A12,WRK(KTMP),MULD2H(ISYMV1,ISYMV2))
      CALL DZERO(A123,N2ORBX)
      CALL DZERO(A,N2ORBX)
      CALL OITH1(ISYMV1,ZYM1,W,A123,1)
      CALL DSCAL(N2ORBX,D1/D6,A123,1)
      CALL DAXPY(N2ORBX,DH,A1,1,A123,1)
      CALL OITH1(ISYMV2,ZYM2,A123,A,ISYMV1)
      CALL OITH1(ISYMV3,ZYM3,A,WRK(KTMP),MULD2H(ISYMV1,ISYMV2))
      CALL DZERO(A,N2ORBX)
      CALL OITH1(ISYMV3,ZYM3,A123,A,ISYMV1)
      CALL OITH1(ISYMV2,ZYM2,A,WRK(KTMP),MULD2H(ISYMV1,ISYMV3))
C
C     / <0| [qj ,T] |0> \
C     | <j| T |0>       |
C     | <0| [qj+,T] |0> |
C     \ -<0| T |j>      /
C
      ISYMDN = 1
      OVLAP  = D1
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
      CALL RSP1GR(NSIM,KZYVR,IDUM,JSPIN,IGRSYM,JSPIN,ISYMV,ETRS,
     *            WRK(KCREF),NZYVEC,NZCVEC,OVLAP,ISYMDN,UDV,WRK(KTMP),
     *            XINDX,MJWOP,WRK(KFREE),LFREE,LORB,LCON,LREF)
C
      CALL PRIRES(ETRS,VECA,IGRSYM,'WCASE4')
C     
      RETURN
      END
      SUBROUTINE C4SOL(VECA,VEC1,VEC2,VEC3,ETRS,XINDX,ZYM1,ZYM2,ZYM3,
     *              UDV,WRK,LWRK,KZYVR,KZYV1,KZYV2,KZYV3,
     *              IGRSYM,ISYMV1,ISYMV2,ISYMV3,CMO,MJWOP,ISYRLM)
C
C     Purpose:
C     Memeory efficient routine for computing solvent contribution 
C     to E[4] times three vectors. Replaces E4SOL in SCF calculations.
C
#include <implicit.h>
C
      PARAMETER ( D0=0.0D0, D1=1.0D0 )
C
#include <maxorb.h>
#include <inforb.h>
#include <infdim.h>
#include <infinp.h>
#include <infvar.h>
#include <infrsp.h>
#include <infpri.h>
#include "rspprp.h"
#include <infcr.h>
#include <inftap.h>
#include <qrinf.h>
C
      DIMENSION ETRS(KZYVR),XINDX(*)
      DIMENSION UDV(NASHDI,NASHDI)
      DIMENSION ZYM1(*),ZYM2(*),ZYM3(*),WRK(*),CMO(*)
      DIMENSION VEC1(KZYV1),VEC2(KZYV2),VEC3(KZYV3),VECA(KZYVR)
      DIMENSION MJWOP(2,MAXWOP,8), ISYRLM(2*LSOLMX+1)
C
      LOGICAL LCON, LORB, LREF
C
      NSIM = 1
      KFREE = 1
      LFREE = LWRK
      CALL MEMGET('REAL',KTRES ,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KTELM ,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KTLMA ,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KTLMB ,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KFLST ,NLMSOL,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KFLOP ,NLMSOL,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KTNLM ,NLMSOL,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KRLMAO,NNBASX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KUCMO ,NORBT*NBAST,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KCREF ,NCREF,WRK,KFREE,LFREE)
C
C Get the reference state
C
      CALL GETREF(WRK(KCREF),MZCONF(1))
C
C Zero the final effective operator
C
      CALL DZERO(WRK(KTRES),N2ORBX)
C
C Unpack the response vectors
C
      CALL GTZYMT(NSIM,VEC1,KZYV1,ISYMV1,ZYM1,MJWOP)
      CALL GTZYMT(NSIM,VEC2,KZYV2,ISYMV2,ZYM2,MJWOP)
      CALL GTZYMT(NSIM,VEC3,KZYV3,ISYMV3,ZYM3,MJWOP)
C
C Unpack symmetry blocked CMO
C
      CALL UPKCMO(CMO,WRK(KUCMO))
C
C Calculate f(l) factors.
C
      CALL SOLFL(WRK(KFLST),EPSTAT,RSOL,LSOLMX)
      CALL SOLFL(WRK(KFLOP),EPSOL,RSOL,LSOLMX)
C
C Read nuclear contributions TN(l,m).
C
      IF (LUSOL .LE. 0) CALL GPOPEN(LUSOL,FNSOL,'UNKNOWN',' ',
     &     'UNFORMATTED',IDUMMY,.FALSE.)
      REWIND LUSOL
      CALL MOLLAB('SOLVRLM ',LUSOL,LUERR)
      READ (LUSOL)
      CALL READT(LUSOL,NLMSOL,WRK(KTNLM))
C
C Loop over l,m expansion.
C
      LM = 0
      DO 520 L = 0,LSOLMX
         READ (LUSOL) L1,(ISYRLM(M),M=1,2*L+1)
      DO 500 M = -L,L
         LM = LM + 1
         ISYMT = ISYRLM(L+M+1)
         IF (ISYMT.NE.1 .AND. ISYMT.NE.ISYMV1 .AND.
     &        ISYMT.NE.ISYMV2 .AND. ISYMT.NE.MULD2H(ISYMV1,ISYMV2)
     &        .AND. ISYMT.NE.MULD2H(ISYMV1,MULD2H(ISYMV2,ISYMV3))) THEN
            READ (LUSOL)
            GO TO 500
         END IF
C
C Read R(l,m) in ao basis, transform to mo basis and unpack.
C Electronic contribution TE(l,m).
C
         CALL READT(LUSOL,NNBASX,WRK(KRLMAO))
         CALL UTHU(WRK(KRLMAO),WRK(KTELM),WRK(KUCMO),WRK(KFREE),
     &             NBAST,NORBT)
         CALL DCOPY(NNORBX,WRK(KTELM),1,WRK(KTLMA),1)
         CALL DSPTSI(NORBT,WRK(KTLMA),WRK(KTELM))
C
C Create the effective operator:        
C 
C     TRES = sum(l,m)[ W(k1,k2,k3) + A1(k2,k3) + A12(k3) + A123 ]
C 
C      W(k1,k2,k3) = g(l)*( F1 + F2 )*TELM(k1,k2,k3)
C        A1(k2,k3) = g(l)*F3*TELM(k2,k3)
C          A12(k3) = g(l)*F4*TELM(k3)
C             A123 = g(l)*F5*TELM
C
C      F1 = -TNLM*1/3
C      F2 = 1/3*<0| TELM |0>
C      F3 = <0| TELM(k1) |0>
C      F4 = <0| TELM(k1,k2) |0>
C      F5 = 1/3*<0| TELM(k1,k2,k3) |0>
C
         F1=D0
         F2=D0
         F3=D0
         F4=D0
         F5=D0
C     
         IF (ISYMT.EQ.1) THEN
            F1 = -WRK((KTNLM-1)+LM) 
            CALL DCOPY(N2ORBX,WRK(KTELM),1,WRK(KTLMA),1)
            CALL MELONE(WRK(KTLMA),1,UDV,D1,F2,200,'C4SOL')
            F1 = F1/3
            F2 = F2/3
         END IF
         IF (ISYMT.EQ.ISYMV1) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL OITH1(ISYMV1,ZYM1,WRK(KTELM),WRK(KTLMA),ISYMT)
            CALL MELONE(WRK(KTLMA),1,UDV,D1,F3,200,'C4SOL')
         END IF
         IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL DZERO(WRK(KTLMB),N2ORBX)
            CALL OITH1(ISYMV1,ZYM1,WRK(KTELM),WRK(KTLMA),ISYMT)
            CALL OITH1(ISYMV2,ZYM2,WRK(KTLMA),WRK(KTLMB),ISYMV2)
            CALL MELONE(WRK(KTLMB),1,UDV,D1,F4,200,'C4SOL')
         END IF
         IF (ISYMT.EQ.MULD2H(ISYMV1,MULD2H(ISYMV2,ISYMV3))) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL DZERO(WRK(KTLMB),N2ORBX)
            CALL OITH1(ISYMV1,ZYM1,WRK(KTELM),WRK(KTLMA),ISYMT)
            CALL OITH1(ISYMV2,ZYM2,WRK(KTLMA),WRK(KTLMB),
     &                 MULD2H(ISYMT,ISYMV1))
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL OITH1(ISYMV3,ZYM3,WRK(KTLMB),WRK(KTLMA),ISYMV3)
            CALL MELONE(WRK(KTLMA),1,UDV,D1,F5,200,'C4SOL')
            F5 = F5/3
         END IF
C
         FLST = WRK((KFLST-1)+LM)
         FLOP = WRK((KFLOP-1)+LM)
         IF (ISYMT.EQ.MULD2H(ISYMV1,MULD2H(ISYMV2,ISYMV3))) THEN
            FACT = FLOP*F5
            CALL DAXPY(N2ORBX,FACT,WRK(KTELM),1,WRK(KTRES),1)
         END IF
         IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
            FACT = FLOP*F4
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL OITH1(ISYMV3,ZYM3,WRK(KTELM),WRK(KTLMA),ISYMT)
            CALL DAXPY(N2ORBX,FACT,WRK(KTLMA),1,WRK(KTRES),1)
         END IF
         IF (ISYMT.EQ.ISYMV1) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL DZERO(WRK(KTLMB),N2ORBX)
            CALL OITH1(ISYMV2,ZYM2,WRK(KTELM),WRK(KTLMA),ISYMT)
            CALL OITH1(ISYMV3,ZYM3,WRK(KTLMA),WRK(KTLMB),
     &                 MULD2H(ISYMT,ISYMV2))
            FACT = FLOP*F3
            CALL DAXPY(N2ORBX,FACT,WRK(KTLMB),1,WRK(KTRES),1)
         END IF
         IF (ISYMT.EQ.1) THEN
            IF (INERSI) THEN
               FACT = FLST*(F1+F2)
            ELSE
               FACT = FLOP*(F1+F2)
            END IF
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL DZERO(WRK(KTLMB),N2ORBX)
            CALL OITH1(ISYMV1,ZYM1,WRK(KTELM),WRK(KTLMA),ISYMT)
            CALL OITH1(ISYMV2,ZYM2,WRK(KTLMA),WRK(KTLMB),
     &                 MULD2H(ISYMT,ISYMV1))
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL OITH1(ISYMV3,ZYM3,WRK(KTLMB),WRK(KTLMA),
     &                 MULD2H(ISYMT,MULD2H(ISYMV1,ISYMV2)))
            CALL DAXPY(N2ORBX,FACT,WRK(KTLMA),1,WRK(KTRES),1)
         END IF
  500 CONTINUE
  520 CONTINUE
C
C       Make the gradient
C
C     / <0| [qj ,TRES] |0> \
C     |          0         |
C     | <0| [qj+,TRES] |0> |
C      \         0         /
C
      ISYMDN = 1
      OVLAP  = D1
      JSPIN = 0
      ISYMV  = IREFSY
      ISYMST = MULD2H(IGRSYM,IREFSY)
      IF ( ISYMST .EQ. IREFSY ) THEN
         LCON = ( MZCONF(IGRSYM) .GT. 1 )
      ELSE
         LCON = ( MZCONF(IGRSYM) .GT. 0 )
      END IF
      LORB   = ( MZWOPT(IGRSYM) .GT. 0 )
      LREF = .TRUE.
      NZYVEC = NCREF
      NZCVEC = NCREF
      CALL RSP1GR(NSIM,KZYVR,IDUM,JSPIN,IGRSYM,JSPIN,ISYMV,ETRS,
     *            WRK(KCREF),NZYVEC,NZCVEC,OVLAP,ISYMDN,UDV,WRK(KTRES),
     *            XINDX,MJWOP,WRK(KFREE),LFREE,LORB,LCON,LREF)
C
      CALL GPCLOSE(LUSOL,'KEEP')
      RETURN
      END
