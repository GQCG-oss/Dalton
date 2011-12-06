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
C=============================================================================
C    /* Deck E3IEF */
C=============================================================================
      SUBROUTINE E3IEF(VECA, VEC1, VEC2,ETRS,XINDX,ZYM1,ZYM2,
     *              DEN1,UDV,WRK,LFREE,KZYVR,KZYV1,KZYV2,
     *              IGRSYM,ISYMV1,ISYMV2,CMO,MJWOP)
C
C     Purpose:
C     Outer driver routine for IEF-PCM solvent contribution 
C     to E[3] times two vectors. 
C     Derived from its "older brother", rspsol1.F
C
C     Main comments:
C     1) POPGET contains two cycles over tesserae because the terms
C        are no longer diagonal: the charge on each tessera 
C        does not depend only on the potential on that tessera
C     2) For the time beeing I have not implemented the C3IEF
C        for the memory efficent SCF calculations
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
      CALL DZERO(WRK(KTA), N2ORBX)
      CALL DZERO(WRK(KTB), N2ORBX)
      CALL DZERO(WRK(KTB1),N2ORBX)
      CALL DZERO(WRK(KTB2),N2ORBX)
      CALL DZERO(WRK(KTC1),N2ORBX)
      CALL DZERO(WRK(KTC2),N2ORBX)
      CALL DZERO(WRK(KTD1),N2ORBX)
      CALL DZERO(WRK(KTD2),N2ORBX)
      CALL DZERO(WRK(KTE), N2ORBX)
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
      CALL PCASE1(VECA, VEC1, VEC2,WRK(KTA),WRK(KTB),
     *              ETRS,XINDX,ZYM1,ZYM2,
     *              DEN1,UDV,WRK(KFREE),LFREE,
     *              KZYVR,KZYV1,KZYV2,IGRSYM,ISYMV1,ISYMV2,CMO,MJWOP,
     *              DBGFLG)
C
      CALL PCASE2(VECA,VEC1,VEC2,WRK(KTC1),
     *              ETRS,XINDX,ZYM1,ZYM2,
     *              DEN1,UDV,WRK(KFREE),LFREE,
     *              KZYVR,KZYV1,KZYV2,IGRSYM,ISYMV1,ISYMV2,CMO,MJWOP,
     *              DBGFLG)
C
      CALL PCASE2(VECA,VEC2,VEC1,WRK(KTC2),
     *              ETRS,XINDX,ZYM2,ZYM1,
     *              DEN1,UDV,WRK(KFREE),LFREE,
     *              KZYVR,KZYV2,KZYV1,IGRSYM,ISYMV2,ISYMV1,CMO,MJWOP,
     *              DBGFLG)
C
      CALL PCASE3(VECA,VEC1,VEC2,WRK(KTB),WRK(KTB1),WRK(KTC1),WRK(KTD1),
     *              ETRS,XINDX,ZYM1,ZYM2,
     *              DEN1,UDV,WRK(KFREE),LFREE,
     *              KZYVR,KZYV1,KZYV2,IGRSYM,ISYMV1,ISYMV2,CMO,MJWOP)
C
      CALL PCASE3(VECA,VEC2,VEC1,WRK(KTB),WRK(KTB2),WRK(KTC2),WRK(KTD2),
     *              ETRS,XINDX,ZYM2,ZYM1,
     *              DEN1,UDV,WRK(KFREE),LFREE,
     *              KZYVR,KZYV2,KZYV1,IGRSYM,ISYMV2,ISYMV1,CMO,MJWOP)
C
      CALL PCASE4(VECA, VEC1, VEC2,WRK(KTA),WRK(KTB1),
     *              WRK(KTB2),WRK(KTC1),WRK(KTC2),WRK(KTE),
     *              ETRS,XINDX,ZYM1,ZYM2,
     *              DEN1,UDV,WRK(KFREE),LFREE,
     *              KZYVR,KZYV1,KZYV2,IGRSYM,ISYMV1,ISYMV2,CMO,MJWOP,
     *              DBGFLG)
C
      IF (CRCAL .OR. E3TEST) THEN
         WRITE (LUPRI,'(A,F20.12)') 'Total contribution from E3IEF:',
     *        DDOT(KZYVR,VECA,1,ETRS,1) - VAL
      END IF
C
      RETURN
      END
C=============================================================================
C    /* Deck POPGET */
C=============================================================================
      SUBROUTINE POPGET(ZYM1,ZYM2,ZYM3,IKLVL,OVLAP,ISYMDN,
     *                  ISYMV1,ISYMV2,ISYMV3,TOP,DEN1,NUCFLG,
     *                  CMO,WRK,LWRK,IPRLVL,PRSTR,DFCTR,
     *                  ADDFLG)
C
C Output:
C
C TOP = sum(l,m) f(l) <L| T(l,m)(k1,k2,..) |R> TE(l,m)
C
C POP = sum_i <L|q_i(k1,k2,..)|R> V_i
C
C Input:
C
C OVLAP = <L|R>, DEN1 = <L|...|R> 
C NUCFLG indicates if T(l,m) = TN(l,m) or T(l,m) = TE(l,m)
C IKLVL is the number of times T(l,m) is to be one-index tranformed
C
C PCM Notes
C  1) The two T (TE and TN) are now replaced by charges q and
C     potentials V
C  2) Potentials are always the electronic potentials derived from
C     the MOs (J1 integrals)
C  3) Charges can be - qtot = qed + qei +qn (total charges)
C                    - qe   = qed + qei     (total electronic charges)
C                    - qed  (dynamic electronic charges)
C  4) All transformations are carried out on the charge part and
C     never on the potential part but since the final term is always Vq
C     this is not relevant
C  5) For the response calculation starting from a nonequilibrium
C     state some coding will be needed in order to get all the charges
C     properly: qei in particular and the total charges

C
C
#include <implicit.h>
#include <dummy.h>
C
#include <maxorb.h>
#include <mxcent.h>
#include <pcmdef.h>
#include <priunit.h>
#include <infdim.h>
#include <orgcom.h>
#include <inforb.h>
#include <infpri.h>
#include "infpar.h"
#include <inftap.h>
#include <infinp.h>
#include <infrsp.h>
#include <pcm.h>
#include <pcmlog.h>
C
Clf PCM NUCFLG
C       NUCFLG = 0 MULTIPLY BY ELECTRONIC GHARGES
C       NUCFLG = 1 MULTIPLY BY NUCLEAR GHARGES
C       NUCFLG = 2 MULTIPLY BY ELECTRONIC + NUCLEAR GHARGES
      CHARACTER*(*)PRSTR
      LOGICAL FNDLAB, CLCCHG, EXP1VL, TOFILE, TRIMAT
Clf
      INTEGER ADDFLG
C
      DIMENSION DEN1(NASHDI,NASHDI), INTREP(9*MXCENT), INTADR(9*MXCENT)
      DIMENSION WRK(*), CMO(*)
      DIMENSION ZYM1(NORBT,NORBT),ZYM2(NORBT,NORBT),ZYM3(NORBT,NORBT)
      DIMENSION TOP(N2ORBX)
      CHARACTER*8 LABINT(9*MXCENT)
C
      NSIM = 1

C   When IKLVL is passed as -1 we initialize CLCCHG as false
      IF(IKLVL.LT.0) THEN
         CLCCHG = .FALSE.
      ELSE
         CLCCHG = .TRUE.
      END IF
      IF (IKLVL.LE.0) ISYM = 1
      IF (IKLVL.EQ.1) ISYM = ISYMV1
      IF (IKLVL.EQ.2) ISYM = MULD2H(ISYMV1,ISYMV2)
      IF (IKLVL.EQ.3) ISYM = MULD2H(MULD2H(ISYMV1,ISYMV2),ISYMV3)
      ISYMT = MULD2H(ISYM,ISYMDN)
C
      KFREE = 1
      CALL MEMGET('REAL',KCHGAO,NSYM*NNBASX,WRK,KFREE,LWRK)
      CALL MEMGET('REAL',KUCMO ,NORBT*NBAST,WRK,KFREE,LWRK)
      CALL MEMGET('REAL',KTELM ,N2ORBX,WRK,KFREE,LWRK)
      CALL MEMGET('REAL',KTLMA ,N2ORBX,WRK,KFREE,LWRK)
      CALL MEMGET('REAL',KTLMB ,N2ORBX,WRK,KFREE,LWRK)
      CALL MEMGET('REAL',KTLMC ,N2ORBX,WRK,KFREE,LWRK)
      CALL MEMGET('REAL',KTCHG ,NTS,WRK,KFREE,LWRK)
      CALL MEMGET('REAL',KQCHG ,NTS,WRK,KFREE,LWRK)
      CALL MEMGET('REAL',KMULT ,NTS,WRK,KFREE,LWRK)
C
      CALL DZERO(WRK(KTLMC),N2ORBX)
Clf
C
C     Unpack symmetry blocked CMO
C
      CALL UPKCMO(CMO,WRK(KUCMO))
C
C     Loop over tesserae
C
ckr      REWIND (LUPROP)
C
C     One-index transform charges IKLVL times.
C     The result will be in WRK(KTLMA) and of symmetry ISYM.
C     (ISYM should equal ISYMDN.)
C     
C
C We calculate charges only if necessary
C Otherwise we use QSE or QSN
C
      IF (CLCCHG) THEN
         CALL DZERO(WRK(KTCHG),NTS)
         XI = DIPORG(1)
         YI = DIPORG(2)
         ZI = DIPORG(3)
#if defined (VAR_MPI)
         IF (NODTOT .GE. 1) THEN
            CALL J2XP(IKLVL,ISYMV1,ISYMV2,ISYMV3,ISYMDN,ZYM1,ZYM2,ZYM3,
     &                WRK(KTCHG),OVLAP,WRK(KUCMO),DEN1,WRK(KFREE),LWRK)
         ELSE
#endif
         DO ITS=1,NTSIRR
            CALL DZERO(WRK(KCHGAO),NNBASX)
            L = 1
            NCOMP = NSYM
            DIPORG(1) = XTSCOR(ITS)
            DIPORG(2) = YTSCOR(ITS)
            DIPORG(3) = ZTSCOR(ITS)
            EXP1VL    = .FALSE.
            TOFILE    = .FALSE.
            KPATOM    = 0
            TRIMAT    = .TRUE.
            CALL GET1IN(WRK(KCHGAO),'NPETES ',NCOMP,WRK(KFREE),LWRK,
     &               LABINT,INTREP,INTADR,L,TOFILE,KPATOM,TRIMAT,
     &               DUMMY,EXP1VL,DUMMY,IPRRSP)
            JCHGAO = KCHGAO
            DO ILOP = 1, NSYM
               ISYM = ILOP
               JTS = (ILOP - 1)*NTSIRR + ITS
               CALL UTHU(WRK(JCHGAO),WRK(KTLMA),WRK(KUCMO),WRK(KFREE),
     $                   NBAST,NORBT)
               CALL DSPTSI(NORBT,WRK(KTLMA),WRK(KTELM))
C with no transformation we copy the array from KTELM to KTLMA
               IF (IKLVL.EQ.0) THEN 
                  CALL DZERO(WRK(KTLMA),N2ORBX)
                  CALL DCOPY(N2ORBX,WRK(KTELM),1,WRK(KTLMA),1)
               END IF
C First transformation of charges: q^e_{ab} --> q^e_{ab}({}^1\kappa) 
               IF (IKLVL.GE.1) THEN
                  CALL DZERO(WRK(KTLMA),N2ORBX)
                  CALL OITH1(ISYMV1,ZYM1,WRK(KTELM),WRK(KTLMA),ISYM)
                  ISYM = MULD2H(ISYM,ISYMV1)
               END IF
C Second transformation of charges: q^e_{ab}({}^1\kappa) --> q^e_{ab}({}^1\kappa {}^2\kappa) 
               IF (IKLVL.GE.2) THEN
                  CALL DZERO(WRK(KTLMB),N2ORBX)
                  CALL OITH1(ISYMV2,ZYM2,WRK(KTLMA),WRK(KTLMB),ISYM)
                  CALL DCOPY(N2ORBX,WRK(KTLMB),1,WRK(KTLMA),1)
                  ISYM = MULD2H(ISYM,ISYMV2)
               END IF 
C Third transformation of charges: hope you can figure out the formula.....
               IF (IKLVL.GE.3) THEN
                  CALL DZERO(WRK(KTLMA),N2ORBX)
                  CALL OITH1(ISYMV3,ZYM3,WRK(KTLMB),WRK(KTLMA),ISYM)
                  ISYM = MULD2H(ISYM,ISYMV3)
               END IF
C Contract transformed charges with the density
               CALL MELONE(WRK(KTLMA),ISYM,DEN1,OVLAP,FACT,200,'POPGET')
               WRK(KTCHG + JTS - 1) = FACT
               JCHGAO = JCHGAO + NNBASX
            END DO
         END DO
#if defined (VAR_MPI)
         END IF
#endif
         CALL V2Q(WRK(KFREE),WRK(KTCHG),WRK(KQCHG),QTEXS,NEQRSP)
      END IF
      IF(LUPCMD .GT. 0) THEN
         CALL GPCLOSE(LUPCMD,'KEEP')
      END IF
C     
C     Make charges of the n-index transformed potentials
C
C     
C Construction of the operator (J, X, or whatever....)
C
      IF (CLCCHG) THEN
         CALL DCOPY(NTSIRR,WRK(KQCHG+(ISYMT - 1)*NTSIRR),1,WRK(KMULT),1)
      ELSE IF (NUCFLG .EQ. 0) THEN
         IF (NEQRSP.AND.(ABS(ADDFLG).EQ.1)) THEN
            CALL DCOPY(NTS,QSENEQ,1,WRK(KMULT),1)
         ELSE
            CALL DCOPY(NTS,QSE,1,WRK(KMULT),1)
         END IF
         CALL DSCAL(NTS,-1.0D0,WRK(KMULT),1)
      ELSE IF (NUCFLG .EQ. 1) THEN
         CALL DCOPY(NTS,QSN,1,WRK(KMULT),1)
         CALL DSCAL(NTS,-1.0D0,WRK(KMULT),1)
      ELSE IF (NUCFLG .EQ. 2) THEN
         CALL DCOPY(NTS,QSN,1,WRK(KMULT),1)
         CALL DAXPY(NTS,1.0D0,QSE,1,WRK(KMULT),1)
         CALL DSCAL(NTS,-1.0D0,WRK(KMULT),1)
      ELSE 
         CALL QUIT('CASE NOT DEFINED IN POPGET')
      END IF
      NOSIM = 1
      CALL DZERO(WRK(KCHGAO),NSYM*NNBASX)
      CALL J1INT(WRK(KMULT),.FALSE.,WRK(KCHGAO),NOSIM,.FALSE.,'NPETES ',
     &           ISYMT,WRK(KFREE),LWRK)
      CALL UTHU(WRK(KCHGAO),WRK(KTLMA),WRK(KUCMO),WRK(KFREE),
     &          NBAST,NORBT)
      CALL DSPTSI(NORBT,WRK(KTLMA),WRK(KTLMC))
      IF (ADDFLG .GT. 0) CALL DAXPY(N2ORBX,DFCTR,WRK(KTLMC),1,TOP,1)
C
      DIPORG(1) = XI
      DIPORG(2) = YI
      DIPORG(3) = ZI
C
      IF (IPRRSP.GE.IPRLVL) THEN
         WRITE(LUPRI,'(/3A,2D22.14)') 'Norm of TOPGET in ', PRSTR,
     *        ' : ',DNRM2(N2ORBX,wrk(ktlmc),1),DNRM2(N2ORBX,TOP,1)
      END IF
C
      RETURN
C
      END
C=============================================================================
C    /* Deck PCASE1 */
C=============================================================================
      SUBROUTINE PCASE1(VECA, VEC1, VEC2,TA,TB,
     *              ETRS,XINDX,ZYM1,ZYM2,
     *              DEN1,UDV,WRK,LFREE,KZYVR,KZYV1,KZYV2,
     *              IGRSYM,ISYMV1,ISYMV2,CMO,MJWOP,
     *              DBGFLG)
#include <implicit.h>
#include <dummy.h>
C
      PARAMETER ( D1 = 1.0D0, DM1 = -1.0D0, D2 = 2.0D0, DM2 = -2.0D0,
     $            DP5= 0.5D0)
C
#include <maxorb.h>
C#include <pcmdef.h>
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
      DIMENSION ETRS(KZYVR),XINDX(*)
      DIMENSION UDV(NASHDI,NASHDI),DEN1(NASHDI,NASHDI)
      DIMENSION ZYM1(*),ZYM2(*),WRK(*),CMO(*)
      DIMENSION TA(NORBT,NORBT),TB(NORBT,NORBT)
      DIMENSION VEC1(KZYV1),VEC2(KZYV2),VECA(KZYVR)
      DIMENSION MJWOP(2,MAXWOP,8)
      INTEGER   DBGFLG(10)
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
C PCM  sum_its QSE(its)*V(its)
C
Clf      IF ((ISYMV1.EQ.ISYMV2) .AND. (MZCONF(ISYMV1).GT.0)) THEN
         CALL DZERO(TA,N2ORBX)
         IKLVL = -1
Clf Call n.1
Clf
         CALL POPGET(DUMMY,DUMMY,DUMMY,IKLVL,D1,1,
     *        IDUMMY,IDUMMY,IDUMMY,TA,UDV,0,
     *        CMO,WRK(KFREE),LFREE,100,'TA',D1,
     *        DBGFLG(1))
C
C Scaling not needed for PCM
Clf         CALL DSCAL(N2ORBX,D2,TA,1)
Clf      END IF
C
C PCM     TB = 2*sum(its) [qse(its)+qsn(its)]*TE(its) (TE is the potential operator).
C
C Onsager:TB = 2*sum(l,m) f(l)*( <0|TE(l,m)|0> - Tn(l,m) )*TE(l,m) 
C
Clf      IF (NEQRSP) THEN
Clf         EPSTOP = EPSINF
Clf      ELSE
Clf         EPSTOP = EPS
Clf      END IF
      CALL DZERO(TB,N2ORBX)
Clf Call n.2
Clf
      IKLVL = -1
      CALL POPGET(DUMMY,DUMMY,DUMMY,IKLVL,DUMMY,1,
     *            IDUMMY,IDUMMY,IDUMMY,TB,DUMMY,1,
     *            CMO,WRK(KFREE),LFREE,100,'TB',D1,
     *            DBGFLG(2))
Clf      CALL DSCAL(N2ORBX,DM2,TB,1)
Clf      CALL DSCAL(N2ORBX,DM1,TB,1)
Clf Call n.3
Clf
      IKLVL = -1
      CALL POPGET(DUMMY,DUMMY,DUMMY,IKLVL,D1,1,
     *            IDUMMY,IDUMMY,IDUMMY,TB,UDV,0,
     *            CMO,WRK(KFREE),LFREE,100,'TB',D1,
     *            DBGFLG(3))
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
      CALL PRIRES(ETRS,VECA,IGRSYM,'PCASE1')
C
      RETURN
      END
C=============================================================================
C    /* Deck PCASE2 */
C=============================================================================
      SUBROUTINE PCASE2(VECA, VEC1, VEC2,TC1,
     *              ETRS,XINDX,ZYM1,ZYM2,
     *              DEN1,UDV,WRK,LFREE,KZYVR,KZYV1,KZYV2,
     *              IGRSYM,ISYMV1,ISYMV2,CMO,MJWOP,DBGFLG)
#include <implicit.h>
#include <dummy.h>
C
      PARAMETER ( D1 = 1.0D0, D2 = 2.0D0 , DM1 = -1.0D0 )
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
      DIMENSION ETRS(KZYVR),XINDX(*)
      DIMENSION UDV(NASHDI,NASHDI),DEN1(NASHDI,NASHDI)
      DIMENSION ZYM1(*),ZYM2(*),WRK(*),CMO(*)
      DIMENSION TC1(NORBT,NORBT)
      DIMENSION VEC1(KZYV1),VEC2(KZYV2),VECA(KZYVR)
      DIMENSION MJWOP(2,MAXWOP,8)
      INTEGER   DBGFLG(10)
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
C Onsager: TC1 = 2*sum(l,m) f(l) <0| TE(l,m) (k1) |0> TE(l,m) + ...
C PCM:     TC1 = 2*sum(its)      <0| q_e(its)(k1) |0> V(its)
C
      CALL DZERO(TC1,N2ORBX)
Clf      CALL OUTPUT(zym1,1,NORBT,1,NORBT,NORBT,NORBT,1,
Clf     &     2)
      IF (MZWOPT(ISYMV1).GT.0) THEN
Clf Call n.4
Clf
         IKLVL = 1
         CALL POPGET(ZYM1,DUMMY,DUMMY,IKLVL,D1,1,
     *        ISYMV1,IDUMMY,IDUMMY,TC1,UDV,IDUMMY,
     *        CMO,WRK(KFREE),LFREE,100,'TC1 cont1',DM1,
     *        DBGFLG(4))
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
Clf Call n.5
Clf
         IKLVL = 0
         CALL POPGET(DUMMY,DUMMY,DUMMY,IKLVL,OVLAP,ISYMDN,
     *        IDUMMY,IDUMMY,IDUMMY,TC1,DEN1,0,
     *        CMO,WRK(KFREE),LFREE,100,'TC1 cont2',DM1,
     *        DBGFLG(5))
      END IF
C
Clf      CALL DSCAL(N2ORBX,D2,TC1,1)
Clf      CALL DSCAL(N2ORBX,DM1,TC1,1)
C
      IF (MZCONF(ISYMV2).LE.0) RETURN
C
C
C     /   0    \
C     | Sj(2)  | * <0| TC1 |0>
C     |   0    |
C     \ Sj(2)' /
C
      IF (IGRSYM.EQ.ISYMV2) THEN
         OVLAP = D1
         CALL MELONE(TC1,1,UDV,OVLAP,FACT,IPRONE,'FACT in PCASE2')
         NZCONF = MZCONF(IGRSYM)
         NZVAR  = MZVAR(IGRSYM)
         CALL DAXPY(NZCONF,FACT,VEC2,1,ETRS,1)
         CALL DAXPY(NZCONF,FACT,VEC2(NZVAR+1),1,ETRS(NZVAR+1),1)
      END IF
C
      CALL PRIRES(ETRS,VECA,IGRSYM,'PCASE2')
C
      RETURN
      END
C=============================================================================
C    /* Deck PCASE3 */
C=============================================================================
      SUBROUTINE PCASE3(VECA, VEC1, VEC2,TB,TB1,TC1,TD1,
     *              ETRS,XINDX,ZYM1,ZYM2,
     *              DEN1,UDV,WRK,LFREE,KZYVR,KZYV1,KZYV2,
     *              IGRSYM,ISYMV1,ISYMV2,CMO,MJWOP)
#include <implicit.h>
#include <dummy.h>
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
      DIMENSION ETRS(KZYVR),XINDX(*)
      DIMENSION UDV(NASHDI,NASHDI),DEN1(NASHDI,NASHDI)
      DIMENSION ZYM1(*),ZYM2(*),WRK(*),CMO(*)
      DIMENSION TB(NORBT,NORBT),TB1(NORBT,NORBT)
      DIMENSION TC1(NORBT,NORBT),TD1(NORBT,NORBT)
      DIMENSION VEC1(KZYV1),VEC2(KZYV2),VECA(KZYVR)
      DIMENSION MJWOP(2,MAXWOP,8)
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
C
Clf Moved the check before the two DAXPY: TD1 is used only in this subroutine)
C
      IF (MZCONF(ISYMV2).LE.0) RETURN
C
      CALL DAXPY(N2ORBX,D1,TB1,1,TD1,1)
      CALL DAXPY(N2ORBX,D1,TC1,1,TD1,1)
C the check was here
      CALL MEMGET('REAL',KCREF,MZCONF(1),WRK,KFREE,LFREE)
      CALL GETREF(WRK(KCREF),MZCONF(1))
C
C MCSCF to be checked!!!
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
      CALL PRIRES(ETRS,VECA,IGRSYM,'PCASE3')
C
      RETURN
      END
C=============================================================================
C    /* Deck PCASE4 */
C=============================================================================
      SUBROUTINE PCASE4(VECA, VEC1, VEC2,TA,TB1,TB2,TC1,TC2,TE,
     *              ETRS,XINDX,ZYM1,ZYM2,
     *              DEN1,UDV,WRK,LFREE,KZYVR,KZYV1,KZYV2,
     *              IGRSYM,ISYMV1,ISYMV2,CMO,MJWOP,DBGFLG)
#include <implicit.h>
#include <dummy.h>
C
      PARAMETER ( D1 = 1.0D0, D2 = 2.0D0, DH = 0.5D0 , DMP5 = -0.5D0,
     $     DM1 = -1.0D0)
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
#include <priunit.h>
C
      DIMENSION ETRS(KZYVR),XINDX(*)
      DIMENSION UDV(NASHDI,NASHDI),DEN1(NASHDI,NASHDI)
      DIMENSION ZYM1(*),ZYM2(*),WRK(*),CMO(*)
      DIMENSION TA(NORBT,NORBT),TE(NORBT,NORBT)
      DIMENSION TB1(NORBT,NORBT),TB2(NORBT,NORBT)
      DIMENSION TC1(NORBT,NORBT),TC2(NORBT,NORBT)
      DIMENSION VEC1(KZYV1),VEC2(KZYV2),VECA(KZYVR)
      DIMENSION MJWOP(2,MAXWOP,8)
      INTEGER   DBGFLG(10)
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
Clf Call n.6
Clf
         IKLVL = 2
         CALL POPGET(ZYM1,ZYM2,DUMMY,IKLVL,D1,1,
     *               ISYMV1,ISYMV2,IDUMMY,TE,UDV,IDUMMY,
     *               CMO,WRK(KFREE),LFREE,100,'TE cont1a',DMP5,
     *               DBGFLG(6))
Clf Call n.7
Clf
         IKLVL = 2
         CALL POPGET(ZYM2,ZYM1,DUMMY,IKLVL,D1,1,
     *               ISYMV2,ISYMV1,IDUMMY,TE,UDV,IDUMMY,
     *               CMO,WRK(KFREE),LFREE,100,'TE cont1b',DMP5,
     *               DBGFLG(7))
      END IF
C
C ... + 2*sum(l,m) f(l) ( <01L| TE(l,m)(k2) |0> + 
C                           <0| TE(l,m)(k2) |01R> ) TE(l,m) + ...
C
CLF WE DO NOT NEED THIS FOR PCM
C     Put the factor two into one of the vectors.
C      CALL DSCAL(KZYV1,D2,VEC1,1)
C      CALL DSCAL(NORBT*NORBT,D2,ZYM1,1)
C
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
         IKLVL = 1
Clf Call n.8
Clf
         CALL POPGET(ZYM2,DUMMY,DUMMY,IKLVL,OVLAP,ISYMDN,
     *               ISYMV2,IDUMMY,IDUMMY,TE,DEN1,IDUMMY,
     *               CMO,WRK(KFREE),LFREE,100,'TE cont2a',DM1,
     *               DBGFLG(8))
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
Clf Call n.9
Clf
         IKLVL = 1
         CALL POPGET(ZYM1,DUMMY,DUMMY,IKLVL,OVLAP,ISYMDN,
     *               ISYMV1,IDUMMY,IDUMMY,TE,DEN1,IDUMMY,
     *               CMO,WRK(KFREE),LFREE,100,'TE cont2b',DM1,
     *               DBGFLG(9))
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
Clf Call n.10
Clf
         IKLVL = 0
         CALL POPGET(DUMMY,DUMMY,DUMMY,IKLVL,OVLAP,ISYMDN,
     *               IDUMMY,IDUMMY,IDUMMY,TE,DEN1,0,
     *               CMO,WRK(KFREE),LFREE,100,'TE cont3',DM1,
     *               DBGFLG(10))
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
Clf and consequently we do not need this.....
CC     Restore the vector
C
C      CALL DSCAL(KZYV1,DH,VEC1,1)
CC
Clf
C
      CALL PRIRES(ETRS,VECA,IGRSYM,'PCASE4')
C
      RETURN
      END

