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
Revision 1.2  2000/05/24 19:09:06  hjj
inserted Dalton header
some changes for triplet response with CSF
===========================================================================
#endif
C
      SUBROUTINE E3SOL(VECA, VEC1, VEC2,ETRS,XINDX,ZYM1,ZYM2,
     *              DEN1,UDV,WRK,LFREE,KZYVR,KZYV1,KZYV2,
     *              IGRSYM,ISYMV1,ISYMV2,CMO,MJWOP,ISYRLM)
C
C     Purpose:
C     Outer driver routine for solvent contribution 
C     to E[3] times two vectors.
C
#include <implicit.h>
C
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infdim.h>
#include <infinp.h>
#include <infvar.h>
#include <infrsp.h>
#include <infpri.h>
#include "rspprp.h"
#include <infcr.h>
C
      DIMENSION ETRS(KZYVR),XINDX(*)
      DIMENSION UDV(NASHDI,NASHDI),DEN1(NASHDI,NASHDI)
      DIMENSION ZYM1(*),ZYM2(*),WRK(*),CMO(*)
      DIMENSION VEC1(KZYV1),VEC2(KZYV2),VECA(KZYVR)
      DIMENSION MJWOP(2,MAXWOP,8)
      INTEGER   DBGFLG(10)
      
C
      NSIM = 1
      KFREE = 1
      CALL MEMGET('REAL',KTA ,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KTB ,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KTB1,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KTB2,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KTC1,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KTC2,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KTD1,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KTD2,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KTE ,N2ORBX,WRK,KFREE,LFREE)
C
      CALL DZERO(WRK(KTA),N2ORBX)
      CALL DZERO(WRK(KTB),N2ORBX)
      CALL DZERO(WRK(KTB1),N2ORBX)
      CALL DZERO(WRK(KTB2),N2ORBX)
      CALL DZERO(WRK(KTC1),N2ORBX)
      CALL DZERO(WRK(KTC2),N2ORBX)
      CALL DZERO(WRK(KTD1),N2ORBX)
      CALL DZERO(WRK(KTD2),N2ORBX)
      CALL DZERO(WRK(KTE),N2ORBX)
C
      CALL GTZYMT(NSIM,VEC1,KZYV1,ISYMV1,ZYM1,MJWOP)
      CALL GTZYMT(NSIM,VEC2,KZYV2,ISYMV2,ZYM2,MJWOP)
C
C DBGFLG initialization
C                  1  2  3  4  5  6  7  8  9  10
C     DATA DBGFLG/-1,-2,-3,-4,-5,-6,-7,-8,-9,-10/
      DATA DBGFLG/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10/
C
C     VECA is only available if E3TEST is set
C
      VAL = DDOT(KZYVR,VECA,1,ETRS,1)
C
      CALL TCASE1(VECA, VEC1, VEC2,WRK(KTA),WRK(KTB),
     *              ETRS,XINDX,ZYM1,ZYM2,ISYRLM,
     *              DEN1,UDV,WRK(KFREE),LFREE,
     *              KZYVR,KZYV1,KZYV2,IGRSYM,ISYMV1,ISYMV2,CMO,MJWOP,
     *              DBGFLG)
C
      CALL TCASE2(VECA,VEC1,VEC2,WRK(KTC1),
     *              ETRS,XINDX,ZYM1,ZYM2,ISYRLM,
     *              DEN1,UDV,WRK(KFREE),LFREE,
     *              KZYVR,KZYV1,KZYV2,IGRSYM,ISYMV1,ISYMV2,CMO,MJWOP,
     *              DBGFLG)
C
      CALL TCASE2(VECA,VEC2,VEC1,WRK(KTC2),
     *              ETRS,XINDX,ZYM2,ZYM1,ISYRLM,
     *              DEN1,UDV,WRK(KFREE),LFREE,
     *              KZYVR,KZYV2,KZYV1,IGRSYM,ISYMV2,ISYMV1,CMO,MJWOP,
     *              DBGFLG)
C
      CALL TCASE3(VECA,VEC1,VEC2,WRK(KTB),WRK(KTB1),WRK(KTC1),WRK(KTD1),
     *              ETRS,XINDX,ZYM1,ZYM2,ISYRLM,
     *              DEN1,UDV,WRK(KFREE),LFREE,
     *              KZYVR,KZYV1,KZYV2,IGRSYM,ISYMV1,ISYMV2,CMO,MJWOP,
     *              DBGFLG)
C
      CALL TCASE3(VECA,VEC2,VEC1,WRK(KTB),WRK(KTB2),WRK(KTC2),WRK(KTD2),
     *              ETRS,XINDX,ZYM2,ZYM1,ISYRLM,
     *              DEN1,UDV,WRK(KFREE),LFREE,
     *              KZYVR,KZYV2,KZYV1,IGRSYM,ISYMV2,ISYMV1,CMO,MJWOP,
     *              DBGFLG)
C
      CALL TCASE4(VECA, VEC1, VEC2,WRK(KTA),WRK(KTB1),
     *              WRK(KTB2),WRK(KTC1),WRK(KTC2),WRK(KTE),
     *              ETRS,XINDX,ZYM1,ZYM2,ISYRLM,
     *              DEN1,UDV,WRK(KFREE),LFREE,
     *              KZYVR,KZYV1,KZYV2,IGRSYM,ISYMV1,ISYMV2,CMO,MJWOP,
     *              DBGFLG)
C
      IF (CRCAL .OR. E3TEST) THEN
         WRITE (LUPRI,'(A,F20.12)') 'Total contribution from E3SOL:',
     *                  DDOT(KZYVR,VECA,1,ETRS,1) - VAL
      END IF
C
      RETURN
      END
      SUBROUTINE TOPGET(ZYM1,ZYM2,ZYM3,IKLVL,OVLAP,ISYMDN,
     *                  ISYMV1,ISYMV2,ISYMV3,TOP,DEN1,EPS,NUCFLG,
     *                  ISYRLM,CMO,WRK,LWRK,IPRLVL,PRSTR,ADDFLG)
C
C Output:
C
C TOP = sum(l,m) f(l) <L| T(l,m)(k1,k2,..) |R> TE(l,m)
C
C Input:
C
C OVLAP = <L|R>, DEN1 = <L|...|R> 
C NUCFLG indicates if T(l,m) = TN(l,m) or T(l,m) = TE(l,m)
C IKLVL is the number of times T(l,m) is to be one-index tranformed
C
#include <implicit.h>
#include <dummy.h>
C
#include <maxorb.h>
#include <priunit.h>
#include <infdim.h>
#include <inforb.h>
#include <infpri.h>
#include <inftap.h>
#include <infinp.h>
#include <infrsp.h>
C
      LOGICAL NUCFLG
      INTEGER ADDFLG
      CHARACTER*(*)PRSTR
C
      DIMENSION DEN1(NASHDI,NASHDI)
      DIMENSION WRK(*), CMO(*)
      DIMENSION ZYM1(NORBT,NORBT),ZYM2(NORBT,NORBT),ZYM3(NORBT,NORBT)
      DIMENSION TOP(N2ORBX)
      DIMENSION ISYRLM(*)
C
      NSIM = 1
C
      IF (IKLVL.EQ.0) ISYM = 1
      IF (IKLVL.EQ.1) ISYM = ISYMV1
      IF (IKLVL.EQ.2) ISYM = MULD2H(ISYMV1,ISYMV2)
      IF (IKLVL.EQ.3) ISYM = MULD2H(MULD2H(ISYMV1,ISYMV2),ISYMV3)
      ISYMT = MULD2H(ISYM,ISYMDN)
C
      KFREE = 1
      CALL MEMGET('REAL',KTNLM ,NLMSOL,WRK,KFREE,LWRK)
      CALL MEMGET('REAL',KRLMAO,NNBASX,WRK,KFREE,LWRK)
      CALL MEMGET('REAL',KUCMO ,NORBT*NBAST,WRK,KFREE,LWRK)
      CALL MEMGET('REAL',KTELM ,N2ORBX,WRK,KFREE,LWRK)
      CALL MEMGET('REAL',KTLMA ,N2ORBX,WRK,KFREE,LWRK)
      CALL MEMGET('REAL',KTLMB ,N2ORBX,WRK,KFREE,LWRK)
Clf
      CALL MEMGET('REAL',KTLMC ,N2ORBX,WRK,KFREE,LWRK)
Clf
      CALL MEMGET('REAL',KFLVEC,NLMSOL,WRK,KFREE,LWRK)
Clf      
      CALL DZERO(WRK(KTLMC),N2ORBX)

C
C     Unpack symmetry blocked CMO
C
      CALL UPKCMO(CMO,WRK(KUCMO))
C
C     Calculate f(l) factors.
C
      CALL SOLFL(WRK(KFLVEC),EPS,RSOL,LSOLMX)
C
C     Read nuclear contributions TN(l,m).
C
      IF (LUSOL.LE.0)
     &CALL GPOPEN(LUSOL,FNSOL,'OLD',' ','UNFORMATTED',IDUMMY,.FALSE.)
      REWIND LUSOL
      CALL MOLLAB('SOLVRLM ',LUSOL,LUERR)
      READ (LUSOL)
      CALL READT(LUSOL,NLMSOL,WRK(KTNLM))
C
C     Loop over l,m expansion.
C
      LM = 0
      DO 520 L = 0,LSOLMX
         READ (LUSOL) L1,(ISYRLM(M),M=1,2*L+1)
      DO 500 M = -L,L
         LM = LM + 1
         IF (ISYRLM(L+M+1) .NE. ISYMT) THEN
            READ (LUSOL)
            GO TO 500
         END IF
C
C     Read R(l,m) in ao basis, transform to mo basis and unpack.
C     Electronic contribution TE(l,m).
C
         CALL READT(LUSOL,NNBASX,WRK(KRLMAO))
         CALL UTHU(WRK(KRLMAO),WRK(KTELM),WRK(KUCMO),WRK(KFREE),
     *             NBAST,NORBT)
         CALL DCOPY(NNORBX,WRK(KTELM),1,WRK(KTLMA),1)
         CALL DSPTSI(NORBT,WRK(KTLMA),WRK(KTELM))
C
C     One-index transform TE(l,m) IKLVL times.
C     The result will be in WRK(KTLMA) and of symmetry ISYM.
C     (ISYM should equal ISYMDN.)
C
         CALL DCOPY(N2ORBX,WRK(KTELM),1,WRK(KTLMA),1)
         ISYM = ISYMT
         IF (IKLVL.GE.1) THEN 
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL OITH1(ISYMV1,ZYM1,WRK(KTELM),WRK(KTLMA),ISYM)
            ISYM = MULD2H(ISYM,ISYMV1)
         END IF
         IF (IKLVL.GE.2) THEN
            CALL DZERO(WRK(KTLMB),N2ORBX)
            CALL OITH1(ISYMV2,ZYM2,WRK(KTLMA),WRK(KTLMB),ISYM)
            CALL DCOPY(N2ORBX,WRK(KTLMB),1,WRK(KTLMA),1)
            ISYM = MULD2H(ISYM,ISYMV2)
         END IF 
         IF (IKLVL.GE.3) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL OITH1(ISYMV3,ZYM3,WRK(KTLMB),WRK(KTLMA),ISYM)
            ISYM = MULD2H(ISYM,ISYMV3)
         END IF
C
C     Add the contribution from TE(l,m) or TN(l,m) to the effective operator.
C
         IF (NUCFLG) THEN
            FACT = WRK((KTNLM-1)+LM) 
         ELSE
            CALL MELONE(WRK(KTLMA),ISYM,DEN1,OVLAP,FACT,200,'TOPGET')
         END IF
         FACT = WRK((KFLVEC-1)+LM)*FACT
         CALL DAXPY(N2ORBX,FACT,WRK(KTELM),1,WRK(KTLMC),1)
         IF (ADDFLG.GT.0) THEN
            CALL DAXPY(N2ORBX,FACT,WRK(KTELM),1,TOP,1)
         END IF
C
  500 CONTINUE
  520 CONTINUE
C
      CALL GPCLOSE(LUSOL,'KEEP')
C
      IF (IPRRSP.GE.IPRLVL) THEN
         WRITE(LUPRI,'(/3A,2D22.14)') 'Norm of TOPGET in ', PRSTR,
     *        ' : ',DNRM2(N2ORBX,wrk(ktlmc),1),DNRM2(N2ORBX,top,1)
      END IF
C
      RETURN
C
      END
      SUBROUTINE TCASE1(VECA, VEC1, VEC2,TA,TB,
     *              ETRS,XINDX,ZYM1,ZYM2,ISYRLM,
     *              DEN1,UDV,WRK,LFREE,KZYVR,KZYV1,KZYV2,
     *              IGRSYM,ISYMV1,ISYMV2,CMO,MJWOP,
     *     DBGFLG)
#include <implicit.h>
C
      PARAMETER ( D1 = 1.0D0, DM1 = -1.0D0, D2 = 2.0D0, DM2 = -2.0D0 )
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
      DIMENSION ETRS(KZYVR),XINDX(*),ISYRLM(2*LSOLMX+1)
      DIMENSION UDV(NASHDI,NASHDI),DEN1(NASHDI,NASHDI)
      DIMENSION ZYM1(*),ZYM2(*),WRK(*),CMO(*)
      DIMENSION TA(NORBT,NORBT),TB(NORBT,NORBT)
      DIMENSION VEC1(KZYV1),VEC2(KZYV2),VECA(KZYVR)
      DIMENSION MJWOP(2,MAXWOP,8)
      INTEGER DBGFLG(10)
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
C
C TA = 2*sum(l,m) f(l) <0| TE(l,m) |0> TE(l,m)
C
      CALL DZERO(TA,N2ORBX)
      CALL TOPGET(DUMMY,DUMMY,DUMMY,0,D1,1,
     *            IDUMMY,IDUMMY,IDUMMY,TA,UDV,EPSOL,.FALSE.,
     *            ISYRLM,CMO,WRK(KFREE),LFREE,100,'TA',
     $     DBGFLG(1))
      CALL DSCAL(N2ORBX,D2,TA,1)
C
C TB = 2*sum(l,m) f(l)*( <0|TE(l,m)|0> - Tn(l,m) )*TE(l,m) 
C
      IF (INERSI) THEN
         EPS = EPSTAT
      ELSE
         EPS = EPSOL
      END IF
      CALL DZERO(TB,N2ORBX)
      CALL TOPGET(DUMMY,DUMMY,DUMMY,0,DUMMY,1,
     *            IDUMMY,IDUMMY,IDUMMY,TB,DUMMY,EPS,.TRUE.,
     *            ISYRLM,CMO,WRK(KFREE),LFREE,100,'TB',
     $     DBGFLG(2))
clf
      CALL DSCAL(N2ORBX,DM1,TB,1)
      CALL TOPGET(DUMMY,DUMMY,DUMMY,0,D1,1,
     *            IDUMMY,IDUMMY,IDUMMY,TB,UDV,EPS,.FALSE.,
     *            ISYRLM,CMO,WRK(KFREE),LFREE,100,'TB',
     $     DBGFLG(3))
clf
      CALL DSCAL(N2ORBX,D2,TB,1)
C
      IF (MZCONF(ISYMV1) .EQ. 0 .OR. MZCONF(ISYMV2) .EQ. 0) RETURN
C
C     /   <01L| [qj,TB] |02R>  + <02L| [qj,TB] |01R>  \
C     |                       0                       |
C     |   <01L| [qj+,TB] |02R> + <02L| [qj+,TB] |01R> |
C     \                       0                       /
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
C     Make the gradient
C
      IF ( MZWOPT(IGRSYM) .GT. 0 ) THEN
         CALL ORBSX(NSIM,IGRSYM,KZYVR,ETRS,TB,OVLAP,ISYMDN,
     *              DEN1,MJWOP,WRK(KFREE),LFREE) 
      END IF
C
      CALL PRIRES(ETRS,VECA,IGRSYM,'TCASE1')
C
      RETURN
      END
      SUBROUTINE TCASE2(VECA, VEC1, VEC2,TC1,
     *              ETRS,XINDX,ZYM1,ZYM2,ISYRLM,
     *              DEN1,UDV,WRK,LFREE,KZYVR,KZYV1,KZYV2,
     *              IGRSYM,ISYMV1,ISYMV2,CMO,MJWOP,
     *     DBGFLG)
#include <implicit.h>
C
      PARAMETER ( D1 = 1.0D0, D2 = 2.0D0 )
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
      DIMENSION ETRS(KZYVR),XINDX(*),ISYRLM(2*LSOLMX+1)
      DIMENSION UDV(NASHDI,NASHDI),DEN1(NASHDI,NASHDI)
      DIMENSION ZYM1(*),ZYM2(*),WRK(*),CMO(*)
      DIMENSION TC1(NORBT,NORBT)
      DIMENSION VEC1(KZYV1),VEC2(KZYV2),VECA(KZYVR)
      DIMENSION MJWOP(2,MAXWOP,8)
      INTEGER DBGFLG(10)
C
      LOGICAL   TDM, LREF, NORHO2
C
C     Initialise variables
C
      JSPIN = 0
      TDM    = .TRUE.
      IPRONE = 200
      KFREE = 1
      NORHO2 = .TRUE.
      NSIM = 1
C
      CALL MEMGET('REAL',KCREF,MZCONF(1),WRK,KFREE,LFREE)
      CALL GETREF(WRK(KCREF),MZCONF(1))
C
C TC1 = 2*sum(l,m) f(l) <0| TE(l,m)(k1) |0> TE(l,m) + ...
C
      CALL DZERO(TC1,N2ORBX)
      IF (MZWOPT(ISYMV1).GT.0) THEN
         CALL TOPGET(ZYM1,DUMMY,DUMMY,1,D1,1,
     *        ISYMV1,IDUMMY,IDUMMY,TC1,UDV,EPSOL,.FALSE.,
     *        ISYRLM,CMO,WRK(KFREE),LFREE,100,'TC1 cont1',
     $        DBGFLG(4))
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
     *        IDUMMY,IDUMMY,IDUMMY,TC1,DEN1,EPSOL,.FALSE.,
     *        ISYRLM,CMO,WRK(KFREE),LFREE,100,'TC1 cont2',
     $        DBGFLG(5))
      END IF
C
      CALL DSCAL(N2ORBX,D2,TC1,1)
C
      IF (MZCONF(ISYMV2).LE.0) RETURN
C
C     /   0    \
C     | Sj(2)  | * <0| TC1 |0>
C     |   0    |
C     \ Sj(2)' /
C
      IF (IGRSYM.EQ.ISYMV2) THEN
         OVLAP = D1
         CALL MELONE(TC1,1,UDV,OVLAP,FACT,IPRONE,'FACT in TCASE2')
         NZCONF = MZCONF(IGRSYM)
         NZVAR  = MZVAR(IGRSYM)
         CALL DAXPY(NZCONF,FACT,VEC2,1,ETRS,1)
         CALL DAXPY(NZCONF,FACT,VEC2(NZVAR+1),1,ETRS(NZVAR+1),1)
      END IF
C
      CALL PRIRES(ETRS,VECA,IGRSYM,'TCASE2')
C
      RETURN
      END
      SUBROUTINE TCASE3(VECA, VEC1, VEC2,TB,TB1,TC1,TD1,
     *              ETRS,XINDX,ZYM1,ZYM2,ISYRLM,
     *              DEN1,UDV,WRK,LFREE,KZYVR,KZYV1,KZYV2,
     *              IGRSYM,ISYMV1,ISYMV2,CMO,MJWOP,
     *     DBGFLG)
#include <implicit.h>
C
      PARAMETER ( D1 = 1.0D0 )
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
      DIMENSION ETRS(KZYVR),XINDX(*),ISYRLM(2*LSOLMX+1)
      DIMENSION UDV(NASHDI,NASHDI),DEN1(NASHDI,NASHDI)
      DIMENSION ZYM1(*),ZYM2(*),WRK(*),CMO(*)
      DIMENSION TB(NORBT,NORBT),TB1(NORBT,NORBT)
      DIMENSION TC1(NORBT,NORBT),TD1(NORBT,NORBT)
      DIMENSION VEC1(KZYV1),VEC2(KZYV2),VECA(KZYVR)
      DIMENSION MJWOP(2,MAXWOP,8)
      INTEGER DBGFLG(10)
C
      LOGICAL   LCON, LORB
      LOGICAL   TDM, LREF, NORHO2
C
C     Initialise variables
C
      JSPIN  = 0
      TDM    = .TRUE.
      KFREE = 1
      IPRONE = -1
      NORHO2 = .TRUE.
      NSIM = 1
C
C TD1 = TB1 + TC1, TB1 = TB(k1)
C
      CALL DZERO(TB1,N2ORBX)
      CALL DZERO(TD1,N2ORBX)
      CALL OITH1(ISYMV1,ZYM1,TB,TB1,1)
      CALL DAXPY(N2ORBX,D1,TB1,1,TD1,1)
      CALL DAXPY(N2ORBX,D1,TC1,1,TD1,1)
C
      IF (MZCONF(ISYMV2).LE.0) RETURN
C
      CALL MEMGET('REAL',KCREF,MZCONF(1),WRK,KFREE,LFREE)
      CALL GETREF(WRK(KCREF),MZCONF(1))
C
C     /   <0| [qj,TD1] |02R>  + <02L| [qj,TD1] |0>  \
C     |   <j| TD1 |02R>                             |
C     |   <0| [qj+,TD1] |02R> + <02L| [qj+,TD1] |0> |
C     \  -<02L| TD1 |j>                             /
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
      CALL RSP1GR(NSIM,KZYVR,IDUMMY,JSPIN,IGRSYM,JSPIN,ISYMV2,ETRS,
     *            VEC2,NZYVEC,NZCVEC,OVLAP,ISYMDN,DEN1,TD1,
     *            XINDX,MJWOP,WRK(KFREE),LFREE,LORB,LCON,LREF)
C
      CALL PRIRES(ETRS,VECA,IGRSYM,'TCASE3')
C
      RETURN
      END
      SUBROUTINE TCASE4(VECA, VEC1, VEC2,TA,TB1,TB2,TC1,TC2,TE,
     *              ETRS,XINDX,ZYM1,ZYM2,ISYRLM,
     *              DEN1,UDV,WRK,LFREE,KZYVR,KZYV1,KZYV2,
     *              IGRSYM,ISYMV1,ISYMV2,CMO,MJWOP,
     *     DBGFLG)
#include <implicit.h>
C
      PARAMETER ( D1 = 1.0D0, D2 = 2.0D0, DH = 0.5D0 )
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
      DIMENSION ETRS(KZYVR),XINDX(*),ISYRLM(2*LSOLMX+1)
      DIMENSION UDV(NASHDI,NASHDI),DEN1(NASHDI,NASHDI)
      DIMENSION ZYM1(*),ZYM2(*),WRK(*),CMO(*)
      DIMENSION TA(NORBT,NORBT),TE(NORBT,NORBT)
      DIMENSION TB1(NORBT,NORBT),TB2(NORBT,NORBT)
      DIMENSION TC1(NORBT,NORBT),TC2(NORBT,NORBT)
      DIMENSION VEC1(KZYV1),VEC2(KZYV2),VECA(KZYVR)
      DIMENSION MJWOP(2,MAXWOP,8)
      INTEGER DBGFLG(10)
C
      LOGICAL   LCON, LORB
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
C TE = 1/2 * TB1(k2) + 1/2 * TB2(k1) + TC1(k2) + TC2(k1) + ...
C
      CALL DZERO(TE,N2ORBX)
      CALL OITH1(ISYMV2,ZYM2,TB1,TE,ISYMV1)
      CALL OITH1(ISYMV1,ZYM1,TB2,TE,ISYMV2)
      CALL DSCAL(N2ORBX,DH,TE,1)
      CALL OITH1(ISYMV2,ZYM2,TC1,TE,ISYMV1)
      CALL OITH1(ISYMV1,ZYM1,TC2,TE,ISYMV2)
C
C ... + ( S(1)S(2)' + S(2)S(1)' ) * TA + ...
C
      IF ((ISYMV1.EQ.ISYMV2) .AND. (MZCONF(ISYMV1).GT.0)) THEN
         NZCONF = MZCONF(ISYMV1)
         NZVAR = MZVAR(ISYMV1)
         FACT = DDOT(NZCONF,VEC1,1,VEC2(NZVAR+1),1) +
     *        DDOT(NZCONF,VEC2,1,VEC1(NZVAR+1),1)
         CALL DAXPY(N2ORBX,FACT,TA,1,TE,1)
      END IF
C
C ... + sum(l,m) f(l) <0| TE(l,m)(k1,k2) |0> TE(l,m)
C     + sum(l,m) f(l) <0| TE(l,m)(k2,k1) |0> TE(l,m) + ...
C
      IF (MZWOPT(ISYMV1).GT.0 .AND. MZWOPT(ISYMV2).GT.0) THEN
         CALL TOPGET(ZYM1,ZYM2,DUMMY,2,D1,1,
     *               ISYMV1,ISYMV2,IDUMMY,TE,UDV,EPSOL,.FALSE.,
     *               ISYRLM,CMO,WRK(KFREE),LFREE,100,'TE cont1a',
     $        DBGFLG(6))
         CALL TOPGET(ZYM2,ZYM1,DUMMY,2,D1,1,
     *               ISYMV2,ISYMV1,IDUMMY,TE,UDV,EPSOL,.FALSE.,
     *               ISYRLM,CMO,WRK(KFREE),LFREE,100,'TE cont1b',
     $        DBGFLG(7))
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
     *               ISYMV2,IDUMMY,IDUMMY,TE,DEN1,EPSOL,.FALSE.,
     *               ISYRLM,CMO,WRK(KFREE),LFREE,100,'TE cont2a',
     $        DBGFLG(8))
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
     *               ISYMV1,IDUMMY,IDUMMY,TE,DEN1,EPSOL,.FALSE.,
     *               ISYRLM,CMO,WRK(KFREE),LFREE,100,'TE cont2b',
     $        DBGFLG(9))
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
     *               IDUMMY,IDUMMY,IDUMMY,TE,DEN1,EPSOL,.FALSE.,
     *               ISYRLM,CMO,WRK(KFREE),LFREE,100,'TE cont3',
     $        DBGFLG(10))
      END IF
C
C     / <0| [qj ,TE] |0> \
C     | <j| TE |0>       |
C     | <0| [qj+,TE] |0> |
C     \ -<0| TE |j>      /
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
      CALL RSP1GR(NSIM,KZYVR,IDUMMY,JSPIN,IGRSYM,JSPIN,ISYMV,ETRS,
     *            WRK(KCREF),NZYVEC,NZCVEC,OVLAP,ISYMDN,UDV,TE,
     *            XINDX,MJWOP,WRK(KFREE),LFREE,LORB,LCON,LREF)
C
C     Restore the vector
C
      CALL DSCAL(KZYV1,DH,VEC1,1)
C
      CALL PRIRES(ETRS,VECA,IGRSYM,'TCASE4')
C
      RETURN
      END
      SUBROUTINE C3SOL(VECA,VEC1,VEC2,ETRS,XINDX,ZYM1,ZYM2,
     &                 UDV,WRK,LWRK,KZYVR,KZYV1,KZYV2,
     &                 IGRSYM,ISYMV1,ISYMV2,CMO,MJWOP,ISYRLM)
C
C     Purpose:
C     Memeory efficient routine for computing solvent contribution 
C     to E[3] times two vectors. Replaces E3SOL in SCF calculations.
C
#include <implicit.h>
#include <dummy.h>
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
      DIMENSION ZYM1(*),ZYM2(*),WRK(*),CMO(*)
      DIMENSION VEC1(KZYV1),VEC2(KZYV2),VECA(KZYVR)
      DIMENSION MJWOP(2,MAXWOP,8),ISYRLM(2*LSOLMX+1)
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
     &        ISYMT.NE.ISYMV2 .AND. ISYMT.NE.MULD2H(ISYMV1,ISYMV2)) THEN
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
C     TRES = sum(l,m)[ W(k1,k2) + A1(k2) + A12 ]
C 
C      W(k1,k2) = g(l)*( F1 + F2 )*TELM(k1,k2)
C        A1(k2) = g(l)*F3*TELM(k2)
C           A12 = g(l)*F4*TELM
C
C      F1 = -TNLM
C      F2 = <0| TELM |0>
C      F3 = 2*<0| TELM(k1) |0>
C      F4 = <0| TELM(k1,k2) |0>
C
         F1=D0
         F2=D0
         F3=D0
         F4=D0
C     
         IF (ISYMT.EQ.1) THEN
            F1 = -WRK((KTNLM-1)+LM) 
            CALL DCOPY(N2ORBX,WRK(KTELM),1,WRK(KTLMA),1)
            CALL MELONE(WRK(KTLMA),1,UDV,D1,F2,200,'C3SOL')
         END IF
         IF (ISYMT.EQ.ISYMV1) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL OITH1(ISYMV1,ZYM1,WRK(KTELM),WRK(KTLMA),ISYMT)
            CALL MELONE(WRK(KTLMA),1,UDV,D1,F3,200,'C3SOL')
            F3 = 2*F3
         END IF
         IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL DZERO(WRK(KTLMB),N2ORBX)
            CALL OITH1(ISYMV1,ZYM1,WRK(KTELM),WRK(KTLMA),ISYMT)
            CALL OITH1(ISYMV2,ZYM2,WRK(KTLMA),WRK(KTLMB),ISYMV2)
            CALL MELONE(WRK(KTLMB),1,UDV,D1,F4,200,'C3SOL')
         END IF
C
         FLST = WRK((KFLST-1)+LM)
         FLOP = WRK((KFLOP-1)+LM)
         IF (ISYMT.EQ.MULD2H(ISYMV1,ISYMV2)) THEN
            FACT = FLOP*F4
            CALL DAXPY(N2ORBX,FACT,WRK(KTELM),1,WRK(KTRES),1)
         END IF
         IF (ISYMT.EQ.ISYMV1) THEN
            CALL DZERO(WRK(KTLMA),N2ORBX)
            CALL OITH1(ISYMV2,ZYM2,WRK(KTELM),WRK(KTLMA),ISYMT)
            FACT = FLOP*F3
            CALL DAXPY(N2ORBX,FACT,WRK(KTLMA),1,WRK(KTRES),1)
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
            CALL DAXPY(N2ORBX,FACT,WRK(KTLMB),1,WRK(KTRES),1)
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
      CALL RSP1GR(NSIM,KZYVR,IDUMMY,JSPIN,IGRSYM,JSPIN,ISYMV,ETRS,
     *            WRK(KCREF),NZYVEC,NZCVEC,OVLAP,ISYMDN,UDV,WRK(KTRES),
     *            XINDX,MJWOP,WRK(KFREE),LFREE,LORB,LCON,LREF)

C     
      CALL GPCLOSE(LUSOL,'KEEP')
      RETURN
      END
