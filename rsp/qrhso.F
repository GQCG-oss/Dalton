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
#if defined (VAR_DEBUG)
#define DEBUG .TRUE.
#else
#define DEBUG .FALSE.
#endif

C
#ifdef REV_LOG
C=======================================================================
C May 2000 hjj:
C Use MZCONF(1) instead NCREF for CREF (for triplet response w. CSF)
C Removed delete of LUMHSO which caused trouble because LUMHSO was neeeded later!
C=======================================================================
#endif
C  /* Deck x2hso */
      SUBROUTINE X2HSO(NSIM,KZYVR,KZYV,IGRSYM,IGRSPI,ISYMV,ISPINV,
     *                 OPLBL,IOPSYM,IOPSPI,
     *                 GRVEC,X2TRS,VEC,UDV,PV,OPMAT,DEN1,DEN2,XINDX,CMO,
     *                 MJWOP,WRK,LWRK)
C
C HSO[2]*N X2 style, i.e.
C
C                       (  <0|[q,HSO(N)]|0>  )
C         HSO[2]*N =  - (    <j|HSO(N)|0>    )    (Case 1)
C                       (  <0|[q+,HSO(N)]|0> )
C                       (   -<0|HSO(N)|j>    )
C
C
C                       (  <L|[q,HSO]|0> + <0|[q,HSO]|R>   )
C                     - (            <j|HSO|R>             )  (Case 2)
C                       (  <L|[q+,HSO]|0> + <0|[q+,HSO]|R> )
C                       (            -<L|HSO|j>            )
C
C
C                       (  0  )
C          - <0|HSO|0>  ( Sj' )  (Case 3)
C                       (  0  )
C                       ( Sj  )
C
C  Case 3 does not contribute for orbitally non-degenerate states (assumed)
C
#include <implicit.h>
#include <dummy.h>
C
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DM1 = -1.0D0 )
C
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <inforb.h>
#include <infind.h>
#include <infvar.h>
#include <infdim.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include "rspprp.h"
#include <infhyp.h>
#include <qrinf.h>
#include <infpri.h>
#include <inftap.h>
#include <infspi.h>
#include <infhso.h>
#include <trhso.h>
#include <codata.h>
C
C Used from common:
C
      CHARACTER*8 OPLBL
      DIMENSION WRK(*)
      DIMENSION X2TRS(KZYVR), MJWOP(2,MAXWOP,8)
      DIMENSION GRVEC(KZYVR)
      DIMENSION VEC(KZYV)
      DIMENSION DEN1(NASHDI,NASHDI), DEN2(*)
      DIMENSION OPMAT(NORBT,NORBT)
      DIMENSION XINDX(*), CMO(*)
      DIMENSION UDV(NASHDI,NASHDI)
C
      LOGICAL NORHO2, TDM, LORB, LCON, LREFST, LREF
      LOGICAL DO1, DO2
      CHARACTER*8 LABEL,HSOLBL(3)
      LOGICAL ACALL
      COMMON /ACALL/ACALL
      ACALL = .FALSE.
      DATA HSOLBL/'X SPNORB','Y SPNORB','Z SPNORB'/

      IPRHSO=IPRRSP
C
C If expicit one-electron part go normal way
C
      HSOFAC=ALPHAC**2/4
      IF (OPLBL(2:2).EQ.'1') THEN
         CALL X2DRV(NSIM,KZYVR,KZYV,IGRSYM,IGRSPI,ISYMV,ISPINV,
     *           OPLBL,IOPSYM,IOPSPI,GRVEC,X2TRS,
     *           VEC,UDV,OPMAT,DEN1,
     *           XINDX,CMO,MJWOP,WRK,LWRK)
         CALL DSCAL(KZYVR,HSOFAC,X2TRS,1)
         RETURN
      END IF
C
C Everything from AO integrals
C
      IF (TDHF) THEN
         CALL X2HSOAO(
     &      OPLBL,CMO,MJWOP,
     &      KZYV,ISYMV,ISPINV,VEC,
     &      KZYVR,IGRSYM,IGRSPI,X2TRS,
     &      WRK,LWRK)
         CALL DSCAL(KZYVR,-HSOFAC,X2TRS,1)
         RETURN
      END IF
      IF (OPLBL(2:2).EQ.'2') THEN
         DO1 = .FALSE.
         DO2 = .TRUE.
      END IF
      IF (OPLBL(2:2).EQ.' ') THEN
         DO1 = DOSO1
         DO2 = DOSO2
      END IF
      LABEL = OPLBL
      LABEL(2:2) = '1'
      CALL QENTER('X2HSO')
C
C     Layout some workspace
C
      KCREF  = 1
      KZYM   = KCREF + MZCONF(1)
      KX2X   = KZYM  + N2ORBX
      KFREE  = KX2X  + N2ORBX
      LFREE  = LWRK  - KFREE + 1
      IF (LFREE.LT.0) CALL ERRWRK('X2HSO',KFREE-1,LWRK)
C
C Transform if necessary: last transformed component: ILXYZ
C
      IF (OPLBL(1:1).NE.HSOLBL(ILXYZ)(1:1)) THEN
         IF (OPLBL(1:1) .EQ. 'X') ILXYZ=1
         IF (OPLBL(1:1) .EQ. 'Y') ILXYZ=2
         IF (OPLBL(1:1) .EQ. 'Z') ILXYZ=3
         KSYMSO = IOPSYM
         ITRLSO = 2
         IF (IPRHSO.GT.0) WRITE(LUPRI,'(/A,I3)')
     &      ' X2HSO: Transforming spin-orbit component', ILXYZ
         CALL GPOPEN(LUAHSO,'AO2SOINT','OLD',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         CALL TRAHSO(ITRLSO,CMO,WRK(KFREE),LFREE)
         CALL GPCLOSE(LUAHSO,'KEEP')
      END IF
      IF (X2TEST) THEN
         X2TMZC = D0
         X2TMZO = D0
         X2TMYC = D0
         X2TMYO = D0
         X2TM = D0
         KZX2O = MZWOPT(IGRSYM)
         KZX2C = MZCONF(IGRSYM)
         KZC = 1
         KZO = KZC + KZX2C
         KYC = KZO + KZX2O
         KYO = KYC + KZX2C
      END IF
C
      KSYMOP = IOPSYM
      TDM    = .TRUE.
      NORHO2 = .NOT.DOSO2
      IDAX   = LUMHSO
C
C     Get the operator matrix
C     Put the minus sign in the whole term at the end
C
      IF (DO1) THEN
         KSYMP = -1
         CALL PRPGET(LABEL,CMO,OPMAT,KSYMP,ANTSYM,WRK,LWRK,IPRHSO)
         IF (KSYMP .NE. KSYMOP) CALL QUIT(
     &      'ERROR: Unexpected symmetry of property matrix')
      ELSE
         CALL DZERO(OPMAT,N2ORBX)
      END IF
C
C     Case 1
C     ======
      IF (MZWOPT(ISYMV) .EQ. 0 ) GOTO 2000
C
      ISPIN = MULSP(ISPINV,IOPSPI)
      IF (ISPIN.NE.IGRSPI) CALL QUIT('X2HSO: SPIN ERROR')
C
C     Transform the operator
C
      CALL GTZYMT(NSIM,VEC,KZYV,ISYMV,WRK(KZYM),MJWOP)
      CALL DZERO(WRK(KX2X),N2ORBX)
      CALL OITH1(ISYMV,WRK(KZYM),OPMAT,WRK(KX2X),IOPSYM)
C
C     Make copies of the MC density matrix in DEN1
C
      CALL DCOPY(NASHT*NASHT,UDV,1,DEN1,1)
      OVLAP = D1
      ISYMDN = 1
      IKLVL = 1
C
C     Make the gradient
C
      ISYMR  = IREFSY
      INTSYM = KSYMOP
      ISYMST = MULD2H(IGRSYM,IREFSY)
      IF ( ISYMST .EQ. IREFSY ) THEN
         LCON = ( MZCONF(IGRSYM) .GT. 1 )
      ELSE
         LCON = ( MZCONF(IGRSYM) .GT. 0 )
      END IF
      LORB   = ( MZWOPT(IGRSYM) .GT. 0 )
      NZYVEC = MZCONF(1)
      NZCVEC = MZCONF(1)
      LREFST = .TRUE.
      CALL IPSET(0,0,1,1)
      CALL HSO2CR(IGRSYM,IGRSPI,X2TRS,VECDUM,NZYVEC,NZCVEC,ISYMV,ISYMDN,
     *              XINDX,OVLAP,UDV,PV,WRK(KX2X),WRK(KFREE),LFREE,KZYVR,
     *              LCON,LORB,LREFST,IDAX,INTSYM,IOPSPI,ISPINV,
     *              IKLVL,WRK(KZYM),ISYMV,ZYMC,ISYMC,MJWOP)
C     CALL HSO2CR(IGRSYM,IGRSPI,ETRS,VEC,NZYVEC,NZCVEC,ISYMV,ISYMDN,
C    *                  XINDX,OVLAP,DEN1,DEN2,FI,WRK,LWRK,KZYVR,
C    *                  LCON,LORB,LREFST,IDAX,INTSYM,ISPIN1,ISPIN2,
C    *                  IKLVL,ZYMAT1,ISYM1,ZYMAT2,ISYM2,MJWOP)
C     CALL RSP1GR(NSIM,KZYVR,INTSYM,ISPIN,IGRSYM,IGRSPI,ISYMR,X2TRS,
C    *            WRK(KCREF),NZYVEC,NZCVEC,OVLAP,ISYMDN,DEN1,WRK(KX2X),
C    *            XINDX,WRK(KFREE),LFREE,LORB,LCON,LREFST)
C
      IF (X2TEST) THEN
         X2ZO = -DDOT(KZX2O,GRVEC(KZO),1,X2TRS(KZO),1)
         X2ZC = -DDOT(KZX2C,GRVEC(KZC),1,X2TRS(KZC),1)
         X2YO = -DDOT(KZX2O,GRVEC(KYO),1,X2TRS(KYO),1)
         X2YC = -DDOT(KZX2C,GRVEC(KYC),1,X2TRS(KYC),1)
         X2TOT = X2ZO+X2ZC+X2YO+X2YC
         WRITE(LUPRI,'(A,F24.8)') ' X2TEST Case 1 ZO :',X2ZO - X2TMZO
         WRITE(LUPRI,'(A,F24.8)') ' X2TEST Case 1 ZC :',X2ZC - X2TMZC
         WRITE(LUPRI,'(A,F24.8)') ' X2TEST Case 1 YO :',X2YO - X2TMYO
         WRITE(LUPRI,'(A,F24.8)') ' X2TEST Case 1 YC :',X2YC - X2TMYC
         WRITE(LUPRI,'(A,F24.8)') ' X2TEST Case 1    :',X2TOT - X2TM
         X2TMZO = X2ZO
         X2TMZC = X2ZC
         X2TMYO = X2YO
         X2TMYC = X2YC
         X2TM = X2TOT
      END IF
      IF ( IPRRSP .GT. 100 ) THEN
         WRITE(LUPRI,'(A)') ' Case 1 for X2 term'
         CALL OUTPUT(X2TRS,1,KZYVR,1,NSIM,KZYVR,NSIM,1,LUPRI)
      END IF
C
C     Case 2
C     ======
 2000 IF (MZCONF(ISYMV) .EQ. 0 ) GOTO 3000
C     Construct the density matrix <02L|..|0> + <0|..|02R>
C
C
      CALL GETREF(WRK(KCREF),MZCONF(1))
C
      CALL DCOPY(N2ORBX,OPMAT,1,WRK(KX2X),1)
      CALL DZERO(DEN1,NASHT*NASHT)
      OVLAP = D0
      IKLVL = 0
C
      ILSYM  = IREFSY
      IRSYM  = MULD2H(IREFSY,ISYMV)
      NCL    = MZCONF(1)
      NCR    = MZCONF(ISYMV)
      KZVARL = MZCONF(1)
      KZVARR = MZYVAR(ISYMV)
      LREF   = .TRUE.
      ISYMDN = MULD2H(ILSYM,IRSYM)
C
C Densities for orbital gradient: <e(+,-)>,<e(-,+)> for IGRSPI = 0
C Densities for orbital gradient: <e(-,-)>,<e(+,+)> for IGRSPI = 1
C
      CALL RSPGDM(NSIM,ILSYM,IRSYM,NCL,NCR,KZVARL,KZVARR,
     *            WRK(KCREF),VEC,OVLAP,DEN1,DEN2,IGRSPI,1,
     *            TDM,NORHO2,XINDX,WRK,KFREE,LFREE,LREF)
      IF (IGRSPI.EQ.1) THEN
         IS1 = 0
         IS2 = 0
      ELSE
         IS1 = 1
         IS2 = IGRSPI
      END IF
      CALL RSPGDM(NSIM,ILSYM,IRSYM,NCL,NCR,KZVARL,KZVARR,
     *            WRK(KCREF),VEC,OVLAP,DEN1,DEN2(1+N2ASHX*N2ASHX),
     *            IS1,IS2,
     *            TDM,NORHO2,XINDX,
     *            WRK,KFREE,LFREE,LREF)
      CALL IPSET(IGRSPI,1,IS1,IS2)
C
C     Make the gradient
C
      ISYMR  = ISYMV
      INTSYM = KSYMOP
      ISYMST = MULD2H(IGRSYM,IREFSY)
      IF ( ISYMST .EQ. IREFSY ) THEN
         LCON = ( MZCONF(IGRSYM) .GT. 1 )
      ELSE
         LCON = ( MZCONF(IGRSYM) .GT. 0 )
      END IF
      LORB   = ( MZWOPT(IGRSYM) .GT. 0 )
      NZYVEC = MZYVAR(ISYMV)
      NZCVEC = MZCONF(ISYMV)
      LREFST = .FALSE.
      CALL HSO2CR(IGRSYM,IGRSPI,X2TRS,VEC,NZYVEC,NZCVEC,ISYMV,ISYMDN,
     *              XINDX,OVLAP,DEN1,DEN2,WRK(KX2X),
     *              WRK(KFREE),LFREE,KZYVR,
     *              LCON,LORB,LREFST,IDAX,INTSYM,IOPSPI,0,
     *              IKLVL,ZYMB,ISYMB,ZYMC,ISYMC,MJWOP)
C     CALL RSP1GR(NSIM,KZYVR,INTSYM,ISPIN,IGRSYM,IGRSPI,ISYMR,X2TRS,
C    *            VEC,NZYVEC,NZCVEC,OVLAP,ISYMDN,DEN1,OPMAT,
C    *            XINDX,WRK,LWRK,LORB,LCON,LREFST)
C
      IF (X2TEST) THEN
         X2ZO = -DDOT(KZX2O,GRVEC(KZO),1,X2TRS(KZO),1)
         X2ZC = -DDOT(KZX2C,GRVEC(KZC),1,X2TRS(KZC),1)
         X2YO = -DDOT(KZX2O,GRVEC(KYO),1,X2TRS(KYO),1)
         X2YC = -DDOT(KZX2C,GRVEC(KYC),1,X2TRS(KYC),1)
         X2TOT = X2ZO+X2ZC+X2YO+X2YC
         WRITE(LUPRI,'(A,F24.8)') ' X2TEST Case 2 ZO :',X2ZO - X2TMZO
         WRITE(LUPRI,'(A,F24.8)') ' X2TEST Case 2 ZC :',X2ZC - X2TMZC
         WRITE(LUPRI,'(A,F24.8)') ' X2TEST Case 2 YO :',X2YO - X2TMYO
         WRITE(LUPRI,'(A,F24.8)') ' X2TEST Case 2 YC :',X2YC - X2TMYC
         WRITE(LUPRI,'(A,F24.8)') ' X2TEST Case 2    :',X2TOT - X2TM
         X2TMZO = X2ZO
         X2TMZC = X2ZC
         X2TMYO = X2YO
         X2TMYC = X2YC
         X2TM = X2TOT
      END IF
      IF ( IPRRSP .GT. 100 ) THEN
         WRITE(LUPRI,'(A)') ' Case 2 for X2 term'
         CALL OUTPUT(X2TRS,1,KZYVR,1,NSIM,KZYVR,NSIM,1,LUPRI)
      END IF
C
C Skip expectation value of imaginary operator
C
C     Case 3
C     ======
C     Do Sj(2) <0 |op| 0>
C
 3000 CONTINUE
      IF (X2TEST) THEN
         WRITE(LUPRI,'(/A,F24.8)')' X2TEST Final result:  ZO', X2ZO
         WRITE(LUPRI,'(A,F24.8)') ' X2TEST Final result:  ZC', X2ZC
         WRITE(LUPRI,'(A,F24.8)') ' X2TEST Final result:  YO', X2YO
         WRITE(LUPRI,'(A,F24.8)') ' X2TEST Final result:  YC', X2YC
         WRITE(LUPRI,'(A,F24.8)') ' X2TEST Final result:    ', X2TOT
      END IF
C
C  Take care of the minus sign of the whole term
C  and transform to atomic units
C
      CALL DSCAL(KZYVR,-HSOFAC,X2TRS,1)
      CALL QEXIT('X2HSO')
      RETURN
      END
C  /* Deck a2hso */
      SUBROUTINE A2HSO(NSIM,KZYVR,KZYV,IGRSYM,IGRSPI,ISYMV,ISPINV,
     *                 OPLBL,IOPSYM,IOPSPI,
     *                 GRVEC,A2TRS,VEC,UDV,PV,OPMAT,XINDX,CMO,
     *                 MJWOP,WRK,LWRK)
C
C
C HSO[2]*N A2 style, i.e.
C
C                       (  1/2 <0|[q+,HSO(N)]|0> )
C         HSO[2]*N =  - (        <0|HSO(N)|j>    )    (Case 1)
C                       (  1/2 <0|[q,HSO(N)]|0>  )
C                       (        <j|HSO(N)|0>    )
C
C
C                       (        0      )
C                     - (  - <L|HSO|j>  )  (Case 2)
C                       (        0      )
C                       (    <j|HSO|R>  )
C
C
C                       (  0  )
C      - 1/2 <0|HSO|0>  ( Sj' )  (Case 3)
C                       (  0  )
C                       ( Sj  )
C
C  Case 3 does not contribute for orbitally non-degenerate states (assumed)
C
C     Drive the computation of A[2] times one vector
C
#include <implicit.h>
#include <dummy.h>
C
      PARAMETER( D0 = 0.0D0, D1 = 1.0D0, DH = 0.5D0, DM1 = -1.0D0)
      PARAMETER (D2 = 2.0D0)
C
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <infdim.h>
#include <inftap.h>
#include <inforb.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <infvar.h>
#include <infind.h>
#include <qrinf.h>
#include <infpri.h>
#include <infspi.h>
#include <trhso.h>
#include <infhso.h>
#include <codata.h>
C
      CHARACTER*8 OPLBL,HSOLBL(3)
      DATA HSOLBL/'X SPNORB','Y SPNORB','Z SPNORB'/
C
      DIMENSION WRK(*)
      DIMENSION A2TRS(KZYVR), MJWOP(2,MAXWOP,8)
      DIMENSION GRVEC(KZYVR)
      DIMENSION VEC(KZYV)
      DIMENSION OPMAT(NORBT,NORBT)
      DIMENSION XINDX(*)
      DIMENSION UDV(NASHDI,NASHDI)
      DIMENSION PV(*)
      DIMENSION CMO(*)
C
      CHARACTER*8 LABEL
      LOGICAL NORHO2, TDM, LORB, LCON, LREFST, LREF
      LOGICAL DO1, DO2
      LOGICAL ACALL
      COMMON /ACALL/ACALL
      ACALL = .TRUE.
C
C If expicit one-electron part go normal way
C
      HSOFAC=ALPHAC**2/4
      IF (OPLBL(2:2).EQ.'1') THEN
         CALL A2DRV(NSIM,KZYVR,KZYV,IGRSYM,IGRSPI,ISYMV,ISPINV,
     *                 OPLBL,IOPSYM,IOPSPI,
     *                 GRVEC,A2TRS,VEC,UDV,OPMAT,XINDX,CMO,MJWOP,
     *                 WRK,LWRK)
         CALL DSCAL(KZYVR,HSOFAC,A2TRS,1)
         RETURN
      END IF
C
C Everything from AO integrals
C
      IF (TDHF) THEN
         CALL X2HSOAO(
     &      OPLBL,CMO,MJWOP,
     &      KZYV,ISYMV,ISPINV,VEC,
     &      KZYVR,IGRSYM,IGRSPI,A2TRS,
     &      WRK,LWRK)
         CALL DSCAL(KZYVR,-HSOFAC/2,A2TRS,1)
         CALL DSWAP(MZVAR(IGRSYM),A2TRS,1,
     *           A2TRS(MZVAR(IGRSYM)+1),1)
         RETURN
      END IF
      IF (OPLBL(2:2).EQ.'2') THEN
         DO1 = .FALSE.
         DO2 = .TRUE.
      END IF
      IF (OPLBL(2:2).EQ.' ') THEN
         DO1 = DOSO1
         DO2 = DOSO2
      END IF
      LABEL = OPLBL
      LABEL(2:2) = '1'
      CALL QENTER('A2HSO')
C
C     Layout some workspace
C
      KZYM   = 1
      KA2X   = KZYM  + N2ORBX
      KFREE  = KA2X  + N2ORBX
      LFREE  = LWRK  - KFREE
      IF (LFREE.LT.0) CALL ERRWRK('A2HSO',KFREE-1,LWRK)
C
C Transform if necessary: last transformed component: ILXYZ
C
      IF (OPLBL(1:1).NE.HSOLBL(ILXYZ)(1:1)) THEN
         IF (OPLBL(1:1) .EQ. 'X') ILXYZ=1
         IF (OPLBL(1:1) .EQ. 'Y') ILXYZ=2
         IF (OPLBL(1:1) .EQ. 'Z') ILXYZ=3
         KSYMSO = IOPSYM
         ITRLSO = 2
         IF (IPRHSO.GT.0) WRITE(LUPRI,'(/A,I3)')
     &      ' A2HSO: Transforming spin-orbit component', ILXYZ
         CALL GPOPEN(LUAHSO,'AO2SOINT','OLD',' ','UNFORMATTED',
     &               IDUMMY,.FALSE.)
         CALL TRAHSO(ITRLSO,CMO,WRK(KFREE),LFREE)
         CALL GPCLOSE(LUAHSO,'KEEP')
      END IF
C
      TDM    = .TRUE.
      NORHO2 = .NOT.DO2
C
      IF (A2TEST) THEN
         A2TMZC = D0
         A2TMZO = D0
         A2TMYC = D0
         A2TMYO = D0
         A2TM = D0
         KZA2O = MZWOPT(IGRSYM)
         KZA2C = MZCONF(IGRSYM)
         KZC = 1
         KZO = KZC + KZA2C
         KYC = KZO + KZA2O
         KYO = KYC + KZA2C
      END IF
C
C     Get the operator matrix
C     910408-hjaaj: transpose the operator matrix to
C                   get right sign for imaginary operators
C
      KSYMOP = IOPSYM
      IDAX = LUMHSO
      IF (DO1) THEN
         KSYMP = -1
         CALL PRPGET(LABEL,CMO,OPMAT,KSYMP,ANTSYM,WRK,LWRK,IPRHSO)
         IF (KSYMP .NE. KSYMOP) CALL QUIT(
     &      'ERROR: Unexpected symmetry of property matrix')
         CALL DGETRN(OPMAT,NORBT,NORBT)
      ELSE
         CALL DZERO(OPMAT,N2ORBX)
      END IF
C
C     Case 1:
C     =======
C     / 1/2 <0| [qj+,A(K)] |0> \
C     |   - <0| A(K) |j>       |
C     | 1/2 <0| [qj ,A(K)] |0> |
C     \     <j| A(K) |0>       /
C
      IF ( MZWOPT(ISYMV) .EQ. 0 ) GOTO 2000
C
      ISPIN = MULSP(ISPINV,IOPSPI)
      IF (ISPIN.NE.IGRSPI) CALL QUIT('A2HSO: SPIN ERROR')
C
C     Transform the operator
C
      INTSYM = IOPSYM
      CALL GTZYMT(NSIM,VEC,KZYV,ISYMV,WRK(KZYM),MJWOP)
      CALL DZERO(WRK(KA2X),N2ORBX)
      CALL OITH1(ISYMV,WRK(KZYM),OPMAT,WRK(KA2X),IOPSYM)
C
C     Make the gradient
C
      OVLAP  = D1
      IKLVL = 1
      ISYMDN = 1
      ISYMR  = IREFSY
      ISYMST = MULD2H(IGRSYM,IREFSY)
      IF(ISYMST .EQ. IREFSY ) THEN
         LCON = ( MZCONF(IGRSYM) .GT. 1 )
      ELSE
         LCON = ( MZCONF(IGRSYM) .GT. 0 )
      END IF
      LORB = (MZWOPT(IGRSYM) .GT. 0 )
      NZYVEC = MZCONF(1)
      NZCVEC = MZCONF(1)
      LREFST = .TRUE.
      CALL IPSET(0,0,1,1)
      IF (LCON) THEN
         CALL HSO2CR(IGRSYM,IGRSPI,A2TRS,VECDUM,NZYVEC,NZCVEC,
     *              ISYMV,ISYMDN,
     *              XINDX,OVLAP,UDV,PV,WRK(KA2X),WRK(KFREE),LFREE,KZYVR,
     *              LCON,.FALSE.,LREFST,IDAX,INTSYM,IOPSPI,ISPINV,
     *              IKLVL,WRK(KZYM),ISYMV,ZYMC,ISYMC,MJWOP)
C        CALL RSP1GR(NSIM,KZYVR,INTSYM,ISPIN,IGRSYM,IGRSPI,ISYMR,A2TRS,
C    *               WRK(KCREF),NZYVEC,NZCVEC,OVLAP,ISYMDN,UDV,
C    *               WRK(KA2X),XINDX,WRK(KFREE),LFREE,
C    *               .FALSE.,LCON,LREFST)
      END IF
      IF (LORB) THEN
C        multiply A(K) with 1/2 for orbital part
C        multiply ZYMAT with 1/2 for orbital part
         CALL DSCAL(N2ORBX,DH,WRK(KZYM),1)
         CALL DZERO(WRK(KA2X),N2ORBX)
         CALL OITH1(ISYMV,WRK(KZYM),OPMAT,WRK(KA2X),IOPSYM)
         CALL HSO2CR(IGRSYM,IGRSPI,A2TRS,VECDUM,NZYVEC,NZCVEC,
     *           ISYMV,ISYMDN,
     *           XINDX,OVLAP,UDV,PV,WRK(KA2X),WRK(KFREE),LFREE,KZYVR,
     *           .FALSE.,LORB,LREFST,IDAX,INTSYM,IOPSPI,ISPINV,
     *           IKLVL,WRK(KZYM),ISYMV,ZYMC,ISYMC,MJWOP)
C        CALL RSP1GR(NSIM,KZYVR,INTSYM,ISPIN,IGRSYM,IGRSPI,ISYMR,A2TRS,
C    *               WRK(KCREF),NZYVEC,NZCVEC,OVLAP,ISYMDN,UDV,
C    *               WRK(KA2X),XINDX,WRK(KFREE),LFREE,
C    *               LORB,.FALSE.,LREFST)
      END IF
C
      IF (A2TEST) THEN
         A2ZO = DDOT(KZA2O,GRVEC(KZO),1,A2TRS(KYO),1)
         A2ZC = DDOT(KZA2C,GRVEC(KZC),1,A2TRS(KYC),1)
         A2YO = DDOT(KZA2O,GRVEC(KYO),1,A2TRS(KZO),1)
         A2YC = DDOT(KZA2C,GRVEC(KYC),1,A2TRS(KZC),1)
         A2TOT = A2ZO+A2ZC+A2YO+A2YC
         WRITE(LUPRI,'(A,F24.8)') ' A2TEST Case 1 ZO :',A2ZO - A2TMZO
         WRITE(LUPRI,'(A,F24.8)') ' A2TEST Case 1 ZC :',A2ZC - A2TMZC
         WRITE(LUPRI,'(A,F24.8)') ' A2TEST Case 1 YO :',A2YO - A2TMYO
         WRITE(LUPRI,'(A,F24.8)') ' A2TEST Case 1 YC :',A2YC - A2TMYC
         WRITE(LUPRI,'(A,F24.8)') ' A2TEST Case 1    :',A2TOT - A2TM
         A2TMZO = A2ZO
         A2TMZC = A2ZC
         A2TMYO = A2YO
         A2TMYC = A2YC
         A2TM   = A2TOT
      END IF
      IF ( IPRRSP .GT. 30 ) THEN
         WRITE(LUPRI,'(A)') ' Case 1 for A2 term'
         CALL OUTPUT(A2TRS,1,KZYVR,1,NSIM,KZYVR,NSIM,1,LUPRI)
      END IF
C
C     Case 2:
C     =======
2000  IF (MZCONF(ISYMV) .EQ. 0 ) GO TO 4000
C
C     Multiply the factor of one half into the operator matrix
C     for evaluation of Case 2 and 3.
C
C     I don't like it, but we scale the vector for cases 2 and 3
C
C     CALL DSCAL(NORBT*NORBT,DH,OPMAT,1)
      CALL DSCAL(MZCONF(ISYMV),DH,VEC,1)
      CALL DSCAL(MZCONF(ISYMV),DH,VEC(1+MZVAR(ISYMV)),1)
C
      ISPIN = IOPSPI
C
C     Make the gradient
C
      OVLAP  = D0
      ISYMDN = 1
      INTSYM = IOPSYM
      IKLVL = 0
      ISYMR  = ISYMV
      NZYVEC = MZYVAR(ISYMV)
      NZCVEC = MZCONF(ISYMV)
      LORB   = .FALSE.
      ISYMST = MULD2H(IGRSYM,IREFSY)
      IF(ISYMST .EQ. IREFSY ) THEN
         LCON = ( MZCONF(IGRSYM) .GT. 1 )
      ELSE
         LCON = ( MZCONF(IGRSYM) .GT. 0 )
      END IF
      LREFST = .FALSE.
      CALL HSO2CR(IGRSYM,IGRSPI,A2TRS,VEC,NZYVEC,NZCVEC,ISYMV,ISYMDN,
     *              XINDX,OVLAP,UDV,PV,OPMAT,WRK(KFREE),LFREE,KZYVR,
     *              LCON,LORB,LREFST,IDAX,INTSYM,IOPSPI,0,
     *              IKLVL,ZYMB,ISYMB,ZYMC,ISYMC,MJWOP)
C     CALL RSP1GR(NSIM,KZYVR,INTSYM,ISPIN,IGRSYM,IGRSPI,ISYMR,A2TRS,
C    *            VEC,NZYVEC,NZCVEC,OVLAP,ISYMDN,UDV,OPMAT,
C    *            XINDX,WRK,LWRK,LORB,LCON,LREFST)
C
      CALL DSCAL(MZCONF(ISYMV),D2,VEC,1)
      CALL DSCAL(MZCONF(ISYMV),D2,VEC(1+MZVAR(ISYMV)),1)
      IF (A2TEST) THEN
         A2ZO = DDOT(KZA2O,GRVEC(KZO),1,A2TRS(KYO),1)
         A2ZC = DDOT(KZA2C,GRVEC(KZC),1,A2TRS(KYC),1)
         A2YO = DDOT(KZA2O,GRVEC(KYO),1,A2TRS(KZO),1)
         A2YC = DDOT(KZA2C,GRVEC(KYC),1,A2TRS(KZC),1)
         A2TOT = A2ZO+A2ZC+A2YO+A2YC
         WRITE(LUPRI,'(A,F24.8)') ' A2TEST Case 2 ZO :',A2ZO - A2TMZO
         WRITE(LUPRI,'(A,F24.8)') ' A2TEST Case 2 ZC :',A2ZC - A2TMZC
         WRITE(LUPRI,'(A,F24.8)') ' A2TEST Case 2 YO :',A2YO - A2TMYO
         WRITE(LUPRI,'(A,F24.8)') ' A2TEST Case 2 YC :',A2YC - A2TMYC
         WRITE(LUPRI,'(A,F24.8)') ' A2TEST Case 2    :',A2TOT - A2TM
         A2TMZO = A2ZO
         A2TMZC = A2ZC
         A2TMYO = A2YO
         A2TMYC = A2YC
         A2TM   = A2TOT
      END IF
      IF ( IPRRSP .GT. 100 ) THEN
         WRITE(LUPRI,'(A)') ' Case 2 for A2 term'
         CALL OUTPUT(A2TRS,1,KZYVR,1,NSIM,KZYVR,NSIM,1,LUPRI)
      END IF
C
C Skip expectation value of imaginary operator
C
C     Case 3:
C     =======
C     Do Sj(2) <0 |op| 0>
C
4000  CONTINUE
      IF (A2TEST) THEN
         WRITE(LUPRI,'(/A,F24.8)')' A2TEST Final result:  ZO', A2ZO
         WRITE(LUPRI,'(A,F24.8)') ' A2TEST Final result:  ZC', A2ZC
         WRITE(LUPRI,'(A,F24.8)') ' A2TEST Final result:  YO', A2YO
         WRITE(LUPRI,'(A,F24.8)') ' A2TEST Final result:  YC', A2YC
         WRITE(LUPRI,'(A,F24.8)') ' A2TEST Final result:    ', A2TOT
      END IF
      IF ( IPRRSP .GT.  150 ) THEN
         WRITE(LUPRI,'(//A)') ' Gradient before swapping in A2DRV'
         CALL OUTPUT(A2TRS,1,KZYVR,1,NSIM,KZYVR,NSIM,1,LUPRI)
      END IF
C
C     Swap the final result to conform with the notation for A[2]
C
      CALL DSWAP(MZVAR(IGRSYM),A2TRS,1,
     *           A2TRS(MZVAR(IGRSYM)+1),1)
C
      IF ( IPRRSP .GT.  100 ) THEN
         WRITE(LUPRI,'(//A)') ' Gradient after swapping in A2DRV'
         CALL OUTPUT(A2TRS,1,KZYVR,1,NSIM,KZYVR,NSIM,1,LUPRI)
      END IF
      CALL DSCAL(KZYVR,HSOFAC,A2TRS,1)
      CALL QEXIT('A2HSO')
C
      RETURN
      END
C  /* Deck hso2cr */
      SUBROUTINE HSO2CR(IGRSYM,IGRSPI,ETRS,VEC,NZYVEC,NZCVEC,ISYMV,
     *                  ISYMDN,
     *                  XINDX,OVLAP,DEN1,DEN2,FI,WRK,LWRK,KZYVA,
     *                  LCON,LORB,LREFST,IDAX,INTSYM,IOPSPI,ISPINV,
     *                  IKLVL,ZYMAT1,ISYM1,ZYMAT2,ISYM2,MJWOP)
C
C     Layout the core for RSP2GR
C
#include <implicit.h>
C
#include <maxorb.h>
#include <infdim.h>
#include <inforb.h>
#include <infvar.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <qrinf.h>
#include <infpri.h>
#include <infspi.h>
C
      DIMENSION WRK(*)
      DIMENSION ETRS(KZYVA), MJWOP(2,MAXWOP,8)
      DIMENSION XINDX(*)
      DIMENSION DEN1(*), DEN2(N2ASHX*N2ASHX,*)
      DIMENSION FI(*)
      DIMENSION VEC(NZYVEC)
      DIMENSION ZYMAT1(*), ZYMAT2(*)
C
      LOGICAL LCON, LORB, LREFST, HSORUN
C
C     Layout of WRK:
C     We keep the H2AX  array at the beginning, in order to
C     release extra workspace in the gradient routine after processing
C     the orbital part of the linear transformation.
C
      KH2AX  = 1
      IF (LCON) THEN
         IF (DIROIT) THEN
            KFA    = KH2AX  + N2ASHX * N2ASHX * 2
         ELSE
            KFA    = KH2AX  + N2ASHX * N2ASHX
         END IF
      ELSE
         KFA    = KH2AX
      END IF
      KQA    = KFA    + NORBT  * NORBT
      KQB    = KQA    + NORBT  * NASHDI
      KPVD   = KQB    + NORBT  * NASHDI
      KH2    = KPVD   + N2ASHX * N2ASHX
      KH2X   = KH2    + N2ORBX
      KWRKO  = KH2X   + N2ORBX
C
C     Get free workspace for orbital part of linear transformation
C
      LWRKO  = LWRK   - KWRKO
      IF (LWRKO.LT.0) CALL ERRWRK('HSO2CR 1',KWRKO-1,LWRK)
      CALL DZERO(WRK,KWRKO)
C
C     Get space that can be used in configuration part
C     We need the arrays H2XAC and FI, so keep them intact
C
      KWRKC  = KFA
      LWRKC  = LWRK   - KWRKC
      IF (LWRKC.LT.0) CALL ERRWRK('HSO2CR 2',KWRKC-1,LWRK)
C
      HSORUN = .TRUE.
      CALL RSP2GR (IGRSYM,IGRSPI,ETRS,VEC,NZYVEC,NZCVEC,ISYMV,ISYMDN,
     *             FI,WRK(KFA),WRK(KQA),WRK(KQB),WRK(KH2),WRK(KH2X),
     *             OVLAP,DEN1,DEN2,WRK(KPVD),WRK(KH2AX),XINDX,
     *             WRK(KWRKO),LWRKO,WRK(KWRKC),LWRKC,KZYVA,
     *             LCON,LORB,LREFST,IDAX,INTSYM,IOPSPI,ISPINV,IDUM,
     *             IKLVL,ZYMAT1,ISYM1,ZYMAT2,ISYM2,DUM,IDUM,HSORUN,
     *             DUM,DUM,MJWOP)
C
      RETURN
      END
C  /* Deck hsofxd */
      SUBROUTINE HSOFXD (FI,FA,QA,QB,H2AX,
     *                  IDAX,INTSYM, ISYMDN,DEN1,DEN2,
     *                  PVDEN,H2,H2X,WRK,KFREE,LFREE,
     *                  LCON,LORB,IPRFX,IGRSPI,IOPSPI,ISPINV,
     *                  IKLVL,ZYMAT1,ISYM1,ZYMAT2,ISYM2)
C
C
C 920206 Olav Vahtras
C
C     This routine computes the FX matrices, that is, the Fock
C     type matrices FI, FA, QA, and QB with transformed ("X")
C     integrals.  In addition, If LCON true then the H2AX matrix
C     with transformed integrals is extracted for the CI routines.
C
C  The operator is of the form
C  HSO = h(r,s)T(r,s) + (rs|t^u)(e(+,-)  + 2e(-,+))(rstu)
C
C
#include <implicit.h>
C
      PARAMETER ( D0 = 0.0D0, DH = 0.50D0, D2 = 2.0D0 )
      PARAMETER ( IPLUS = 0, IMINUS = 1 )
C
C Used from common blocks:
C  ?
C
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <inforb.h>
#include <infind.h>
#include <infdim.h>
#include <infrsp.h>
#include <infpri.h>
#include <orbtypdef.h>
#include <inftap.h>
#include "rspprp.h"
#include <infhyp.h>
#include <infspi.h>
#include <trhso.h>
#include <cbgetdis.h>
#include <infhso.h>
#include <infdis.h>
C
C
      DIMENSION FI(NORBT,NORBT), FA(NORBT,NORBT)
      DIMENSION QA(NORBT,NASHDI),QB(NORBT,NASHDI)
      DIMENSION H2(NORBT,NORBT), H2X(NORBT,NORBT)
      DIMENSION H2AX(N2ASHX*N2ASHX,*)
      DIMENSION DEN1(NASHDI,NASHDI), DEN2(*), PVDEN(NASHDI,NASHDI)
      DIMENSION ZYMAT1(NORBT,NORBT),ZYMAT2(NORBT,NORBT)
      DIMENSION WRK(*)
C
      DIMENSION NEEDED(6)
      DATA NEEDED /1,1,1,1,1,0/
C
      LOGICAL LCON, LORB
      LOGICAL ACALL
      COMMON /ACALL/ACALL
C
C     Put up the structure for reading the (transformed) integrals:
C     IEND < 0 flags the completeness of the distributions
C
      CALL QENTER('HSOFXD')
      CALL DZERO(H2,N2ORBX)
      IF (IKLVL.EQ.1) NEEDED(6) = 1
      IF (IPRFX.GT.200) THEN
         WRITE(LUPRI,'(/A)') 'ENTER HSOFXD'
         WRITE(LUPRI,'(A,I5)') 'IDAX =',IDAX
         WRITE(LUPRI,'(A,I5)') 'INTSYM =',INTSYM
         WRITE(LUPRI,'(A,I5)') 'ISYMDN =',ISYMDN
         WRITE(LUPRI,'(A,L1)') 'LCON =',LCON
         WRITE(LUPRI,'(A,L2)') 'LORB =',LORB
         WRITE(LUPRI,'(A,I5)') 'IPRFX =',IPRFX
         WRITE(LUPRI,'(A,I5)') 'IGRSPI =',IGRSPI
         WRITE(LUPRI,'(A,I5)') 'IOPSPI =',IOPSPI
         WRITE(LUPRI,'(A,I5)') 'ISPINV =',ISPINV
         WRITE(LUPRI,'(A,I5)') 'IKLVL =',IKLVL
         WRITE(LUPRI,'(A,I5)') 'ISYM1 =',ISYM1
         WRITE(LUPRI,'(A,I5)') 'ISYM2 =',ISYM2
         WRITE(LUPRI,'(A)') 'One-electron FI'
         CALL OUTPUT(FI,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         IF (LORB) THEN
            WRITE(LUPRI,'(A)') 'One-electron density'
            CALL OUTPUT(DEN1,1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
            WRITE(LUPRI,'(A)') 'Two-electron density 1'
            CALL PRIAC2(DEN2,NASHT,LUPRI)
            WRITE(LUPRI,'(A)') 'Two-electron density 2'
            CALL PRIAC2(DEN2(1+N2ASHX*N2ASHX),NASHT,LUPRI)
         END IF
      END IF
      IF (IOPSPI.NE.1) THEN
         WRITE(LUPRI,'(/A)') 'HSOFXD: WRONG VALUE OF IOPSPI'
         WRITE(LUPRI,'(A,I5)') 'IOPSPI =', IOPSPI
         CALL QTRACE(0)
         CALL QUIT('HSOFXD: WRONG VALUE OF IOPSPI')
      END IF
      IDIST = 0
      IPERCD = -1
      KWRK = KFREE
      LWRK = LFREE
      IF (NASHT.GT.0) THEN
         CALL MEMGET('REAL',KH2XAC,N2ASHX,WRK,KWRK,LWRK)
         CALL DZERO(WRK(KH2XAC),N2ASHX)
         IF (LORB) THEN
            LSCR = MAX(NORBT,N2ASHX)
            CALL MEMGET('REAL',KSCR,LSCR,WRK,KWRK,LWRK)
         END IF
      END IF
      IF (.NOT.DOSO2) GOTO 2000
1000  CALL NXTHSO(IC,ID,H2,NEEDED,WRK,KWRK,LWRK,IDIST)
      IF ( IDIST .LT. 0) GO TO 2000
      IF (IC.EQ.ID) GO TO 1000
C
C     Transpose the operator for A[2] calcualtion
C
      IF (ACALL) THEN
         ITMP = IC
         IC = ID
         ID = ITMP
      END IF
C
C     Determine type of distribution
C
      ICTYP = IOBTYP(IC)
      IDTYP = IOBTYP(ID)
      ICDTYP = IDBTYP(ICTYP,IDTYP)
C
C     Determine symmetry of distribution
C
      ICSYM = ISMO( IC )
      IDSYM = ISMO( ID )
      ICDSYM = MULD2H(ICSYM,IDSYM)
C
C
C
      IF ( IPRFX .GT. 300) THEN
         WRITE(LUPRI,'(//A)') 'Characteristics of distribution:'
         WRITE(LUPRI,'(A)')   '================================'
         WRITE(LUPRI,'(A,2I5)')'IC ID   = ', IC,ID
         WRITE(LUPRI,'(3A)')   'TYP     = ', COBTYP(ICTYP),COBTYP(IDTYP)
         WRITE(LUPRI,'(A,2I5)')'Symmetry= ', ICSYM, IDSYM
         WRITE(LUPRI,'(A)')   '================================'
C
         CALL HEADER('Two-electron integrals from this distribution',3)
         CALL OUTPUT(H2,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
      IF (IDAX .EQ. LUMHSO) THEN
C
C Regular MO-integrals
C
         IF (IKLVL.EQ.0) THEN
C
C ( |^)e(+,-)
C
            IF (LORB) THEN
               IPDA = IPDENS(IGRSPI,1)
               IPDB = IPDENS(0,MULSP(IGRSPI,1))
            END IF
            CALL EXTRAC(H2,WRK(KH2XAC))
            CALL FQDIS1(INTSYM,IC,ID,ICDTYP,ICDSYM,ICSYM,IDSYM,
     *                  FI,FA,QA,QB,H2,WRK(KH2XAC),H2AX,
     *                  DEN1,DEN2(IPDA),DEN2(IPDB),ISYMDN,PVDEN,
     *                  LORB,LCON,
     *                  IGRSPI,IPLUS,IMINUS,IPERCD)
C
C 2( |^)e(-,+)
C
            IF (LORB) THEN
               IPDA = IPDENS(MULSP(IGRSPI,1),0)
               IPDB = IPDENS(1,IGRSPI)
            END IF
            CALL DSCAL(N2ORBX,D2,H2,1)
            CALL DSCAL(N2ASHX,D2,WRK(KH2XAC),1)
            CALL FQDIS1(INTSYM,IC,ID,ICDTYP,ICDSYM,ICSYM,IDSYM,
     *                  FI,FA,QA,QB,H2,WRK(KH2XAC),H2AX,
     *                  DEN1,DEN2(IPDA),DEN2(IPDB),ISYMDN,PVDEN,
     *                  LORB,.FALSE.,
     *                  IGRSPI,IMINUS,IPLUS,IPERCD)
         ELSE IF (IKLVL.EQ.1) THEN
            IH2SYM = MULD2H(ICDSYM,INTSYM)
C
C Transform left (~| )
C
            CALL DZERO(H2X,N2ORBX)
            CALL OITH1(ISYM1,ZYMAT1,H2,H2X,IH2SYM)
            INTSY1 = MULD2H(INTSYM,ISYM1)
            IF (IPRFX.GT.300) THEN
                WRITE(LUPRI,'(/A,I5)')'Left transform H2 of symmetry',
     *              IH2SYM
                CALL OUTPUT(H2,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
                WRITE(LUPRI,'(/A,I5)') 'with ZYMAT1 of symmetry',ISYM1
                CALL OUTPUT(ZYMAT1,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
                WRITE(LUPRI,'(/A)') 'to H2X '
                CALL OUTPUT(H2X,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            END IF
            IF (LORB) THEN
               IPDA = IPDENS(MULSP(IGRSPI,ISPINV),1)
               IPDB = IPDENS(ISPINV,MULSP(IGRSPI,1))
            END IF
            CALL EXTRAC(H2X,WRK(KH2XAC))
C
C  (~|^)e(~+,-)
C
            CALL FQDIS1(INTSY1,IC,ID,ICDTYP,ICDSYM,ICSYM,IDSYM,
     *                     FI,FA,QA,QB,H2X,WRK(KH2XAC),H2AX,
     *                     DEN1,DEN2(IPDA),DEN2(IPDB),ISYMDN,PVDEN,
     *                     LORB,LCON,
     *                     IGRSPI,ISPINV,IMINUS,IPERCD)
            IF (LORB) THEN
               IPDA = IPDENS(0,0)
               IPDB = IPDENS(MULSP(ISPINV,1),IGRSPI)
            END IF
            CALL DSCAL(N2ORBX,D2,H2X,1)
            CALL DSCAL(N2ASHX,D2,WRK(KH2XAC),1)
C
C 2(~|^)e(~-,+)
C
            CALL FQDIS1(INTSY1,IC,ID,ICDTYP,ICDSYM,ICSYM,IDSYM,
     *                     FI,FA,QA,QB,H2X,WRK(KH2XAC),H2AX,
     *                     DEN1,DEN2(IPDA),DEN2(IPDB),ISYMDN,PVDEN,
     *                     LORB,.FALSE.,
     *                     IGRSPI,MULSP(ISPINV,IMINUS),IPLUS,IPERCD)
C
C Transform right ( |~^)
C
            IF (LORB) THEN
               IPDA = IPDENS(IGRSPI,MULSP(ISPINV,1))
               IPDB = IPDENS(0,0)
            END IF
            CALL EXTRAC(H2,WRK(KH2XAC))
C
C ( |~^)e(+,~-)
C
            CALL FQDIS2(INTSYM,IC,ID,ICDTYP,ICDSYM,ICSYM,IDSYM,
     *                     FI,FA,QA,QB,H2,WRK(KH2XAC),H2AX(1,2),
     *                     ZYMAT1,ISYM1,
     *                     DEN1,DEN2(IPDA),DEN2(IPDB),ISYMDN,PVDEN,
     *                     LORB,LCON,
     *                     IGRSPI,IPLUS,MULSP(ISPINV,IMINUS),
     *                     IPERCD,WRK(KSCR))
            IF (LORB) THEN
               IPDA = IPDENS(MULSP(IGRSPI,1),ISPINV)
               IPDB = IPDENS(1,MULSP(IGRSPI,ISPINV))
            END IF
            CALL DSCAL(N2ORBX,D2,H2,1)
            CALL DSCAL(N2ASHX,D2,WRK(KH2XAC),1)
C
C 2 ( |~^)e(-,~+)
C
            CALL FQDIS2(INTSYM,IC,ID,ICDTYP,ICDSYM,ICSYM,IDSYM,
     *                     FI,FA,QA,QB,H2,WRK(KH2XAC),H2AX(1,2),
     *                     ZYMAT1,ISYM1,
     *                     DEN1,DEN2(IPDA),DEN2(IPDB),ISYMDN,PVDEN,
     *                     LORB,.FALSE.,
     *                     IGRSPI,IMINUS,ISPINV,IPERCD,WRK(KSCR))
         ELSE
            WRITE(LUPRI,'(/A)') ' WRONG VALUE OF IKLVL IN HSOFXD'
            WRITE(LUPRI,'(/A,I5)') ' IKLVL =',IKLVL
            CALL QUIT(' WRONG VALUE OF IKLVL IN HSOFXD')
         END IF
      ELSE
         WRITE(LUPRI,'(/A,I5)') 'HSOFXD: WRONG VALUE OF UNIT IDAX',IDAX
         WRITE(LUPRI,'(/A,2I5)') ' IDAX .NE. LUMHSO',IDAX,LUMHSO
         CALL QUIT('HSOFXD: WRONG VALUE OF UNIT IDAX')
      END IF
C
C
C
         IF( IPRFX .GT. 300 ) THEN
            WRITE(LUPRI,'(/A)') ' PARTIAL MATRICES IN  HSOFXD'
            WRITE(LUPRI,'(A)')  ' ==========================='
            WRITE(LUPRI,'(//A)') ' Inactive fock matrix'
            CALL OUTPUT(FI,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            WRITE(LUPRI,'(//A)') ' Active fock matrix'
            CALL OUTPUT(FA,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            WRITE(LUPRI,'(//A)') ' Qa matrix'
            CALL OUTPUT(QA,1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
            WRITE(LUPRI,'(//A)') ' Qb matrix'
            CALL OUTPUT(QB,1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
            WRITE(LUPRI,'(/A)') ' Active two-electron integrals'
            CALL OUTPUT(H2AX,1,N2ASHX,1,N2ASHX,N2ASHX,N2ASHX,1,LUPRI)
            IF (IKLVL.EQ.1) THEN
               WRITE(LUPRI,'(/A)') ' Active two-electron integrals'
               CALL OUTPUT(H2AX(1,2),1,N2ASHX,1,N2ASHX,N2ASHX,N2ASHX,1,
     &                     LUPRI)
            END IF
         END IF
C
C
C     We have now completed processing this load, so get the next
C
      GO TO 1000
2000  CONTINUE
C
C Account for 1/2 in definition of two-electron integrals
C
C and for the fact that we calculated 2 times FA
C
      CALL DSCAL(NORBT*NORBT,DH,FA,1)
C
C Set DISTYP for the cases we may have ci-gradients
C
      IF (IKLVL.EQ.0) THEN
         INFDIS(1) = 23
         INFDIS(2) = 24
      END IF
      IF (IKLVL.EQ.1) THEN
         INFDIS(1) = 19
         INFDIS(2) = 20
      END IF
      IF (IPRFX.GT.200) THEN
         WRITE(LUPRI,'(/A)') ' Completed matrices in HSOFXD'
         WRITE(LUPRI,'(/A)') ' Inactive Fock matrix'
         CALL OUTPUT(FI,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         IF (LORB) THEN
            WRITE(LUPRI,'(/A)') ' Active Fock matrix'
            CALL OUTPUT(FA,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            WRITE(LUPRI,'(/A)') ' QA matrix'
            CALL OUTPUT(QA,1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
            WRITE(LUPRI,'(/A)') ' QB matrix'
            CALL OUTPUT(QB,1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
         END IF
         IF (LCON)THEN
            WRITE(LUPRI,'(/A)') ' Active two-electron integrals'
            CALL PRIAC2(H2AX,NASHT,LUPRI)
            IF (IKLVL.EQ.1) THEN
               WRITE(LUPRI,'(/A)') ' Active two-electron integrals'
               CALL PRIAC2(H2AX(1,2),NASHT,LUPRI)
            END IF
         END IF
      END IF
C
      CALL QEXIT('HSOFXD')
      RETURN
      END
C  /* Deck ipdens */
      FUNCTION IPDENS(ISPIN1,ISPIN2)
C
C  Get a pointer to the density <e(ISPIN1,ISPIN2)>
C  The vector INFDEN(4) should be reset when a two-electron density
C  is constructed. In a triplet calculation we normally have two densities
C  stored after one another, e.g. ( <e(+,+)> <e(-,-)> )
C  Referring to these densities INFDEN should be set to (0,0,1,1)
C
#include <infden.h>
#include <inforb.h>
      CALL QENTER('IPDENS')
      IF (ISPIN1.EQ.INFDEN(1) .AND. ISPIN2.EQ.INFDEN(2)) THEN
         IPDENS = 1
      ELSE IF (ISPIN1.EQ.INFDEN(3) .AND. ISPIN2.EQ.INFDEN(4)) THEN
         IPDENS = N2ASHX*N2ASHX + 1
      ELSE
         CALL QTRACE(0)
         WRITE(0,'(/A,4I3)') 'INFDEN =',(INFDEN(I), I=1,4)
         CALL QUIT('IPDENS: REQUIRED DENSITIES ARE NOT STORED')
      END IF
      CALL QEXIT('IPDENS')
      RETURN
      END
C  /* Deck ipset */
      SUBROUTINE IPSET(I,J,K,L)
#include <infden.h>
      INFDEN(1) = I
      INFDEN(2) = J
      INFDEN(3) = K
      INFDEN(4) = L
      RETURN
      END
C
C /* Deck x2hsoao */
C
      SUBROUTINE X2HSOAO(
     &   OPLBL,CMO,MJWOP,
     &   LVEC,VECSYM,VECSPIN,
     &   VEC,
     &   LGR,GRSYM,GRSPIN,
     &   GR,
     &   WORK,LWORK
     &   )
C
C Built gradient of one-index transformed spin-orbit operator
C from AO integrals (single determinant)
C
      IMPLICIT NONE
C
C Input
C
      CHARACTER*8 OPLBL
      REAL*8  CMO(*)
      INTEGER MJWOP(*)
      INTEGER LVEC, VECSYM, VECSPIN
      REAL*8  VEC(*)
      INTEGER LGR, GRSYM, GRSPIN
C
C Output
C
       REAL*8 GR
C
C Work
C
      INTEGER LWORK
      REAL*8 WORK(LWORK)
#include <inforb.h>
#include <priunit.h>
#include <infrsp.h>
#include <maxorb.h>
#include <infinp.h>
#include <aovec.h>
#include <dftcom.h>
#include <dummy.h>
C
C Local
C
      INTEGER KFREE,LFREE,KAO,KMO,KS,KD,KD1,KD2,KF,KF1,KF2
      LOGICAL TRIPLET(2)
      REAL*8  D1
      PARAMETER (D1 = 1.0D0)
C
      INTEGER IPRINT, NDMAT
      INTEGER KJSTRS , KNPRIM , KNCONT , KIORBS , KJORBS , NPAO , KKORBS
      REAL*8       TIMEND,TIMSTR,SO1WAL,SO1CPU
      REAL*8       SO2WAL,SO2CPU,SOCPU,SOWAL, HFXSAV
      INTEGER I2TYP
      INTEGER IPRTWO,IPRNTA,IPRNTB,IPRNTC,IPRNTD
      INTEGER IRNTYP , NUMDIS, IFCTYP(2), ISYMDM(2)
      LOGICAL TKTIME,RTNTWO
C
      IPRINT=IPRRSP
      NDMAT=2
      IF (NASHT.GT.0) 
     &   CALL QUIT('X2HSOAO not implemented for open shells')

      CALL QENTER('X2HSOAO')
      KFREE=1
      LFREE=LWORK
      CALL MEMGET('REAL',KAO,N2BASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KMO,N2ORBX,WORK,KFREE,LFREE)
C
C Unpack vectors
C
      CALL DZERO(WORK(KMO),N2ORBX)
      CALL GTZYMT(1,VEC,LVEC,VECSYM,WORK(KMO),MJWOP)
C
C Ao-tranformed kappa
C
      CALL DZERO(WORK(KAO),N2BASX)
      CALL MOTOAO(WORK(KMO),WORK(KAO),CMO,VECSYM,WORK(KFREE),LFREE)
      CALL MEMREL('MOTOAO',WORK,KMO,KMO,KFREE,LFREE)
      IF (DEBUG) THEN
         CALL MAOPRI(WORK(KAO),'X2HSOAO:kappa(ao)')
      END IF
C
C Get overlap
C
      CALL MEMGET('REAL',KS,N2BASX,WORK,KFREE,LFREE)
      CALL GET_H1(WORK(KS),'OVERLAP',WORK(KFREE),LFREE)
C
C Allocate densities and Fock
C
      CALL MEMGET('REAL',KF,2*N2BASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KD,2*N2BASX,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KF),2*N2BASX)
      CALL DZERO(WORK(KD),2*N2BASX)
      KD1=KD
      KD2=KD+N2BASX
      KF1=KF
      KF2=KF+N2BASX
C
C Transform density 
C
      CALL GTDMSO(WORK(KFREE),CMO,WORK(KD1),WORK(KFREE),WORK(KFREE))
      CALL D_K(NBAST,WORK(KAO),WORK(KD1),WORK(KD2),WORK(KS),
     &         WORK(KFREE),LFREE)
      IF (DEBUG) THEN
         CALL MAOPRI(WORK(KD1),'X2HSOAO:Density 1')
         CALL MAOPRI(WORK(KD2),'X2HSOAO:Density 2')
      END IF
C
C Build fock
C
      IF (DIRFCK) THEN
         IF (OPLBL(1:1).EQ.'X') IFCTYP(1)=1
         IF (OPLBL(1:1).EQ.'Y') IFCTYP(1)=2
         IF (OPLBL(1:1).EQ.'Z') IFCTYP(1)=3
         IFCTYP(2) = IFCTYP(1)
         IF (GRSPIN.NE.VECSPIN) IFCTYP(1) = -IFCTYP(1)
         IF (GRSPIN.EQ.1)       IFCTYP(2) = -IFCTYP(2)
C
C        Transform density to AO basis 
C
         CALL DSOTAO(WORK(KD1),WORK(KF1),NBAST,0,IPRINT)
         CALL DSOTAO(WORK(KD2),WORK(KF2),NBAST,VECSYM-1,IPRINT)
         CALL DCOPY(NDMAT*N2BASX,WORK(KF),1,WORK(KD),1)
C
C Setup for TWOINT
C
         NPAO = MXSHEL*MXAOVC
         CALL MEMGET('INTE',KJSTRS,NPAO*2,WORK,KFREE,LFREE)     
         CALL MEMGET('INTE',KNPRIM,NPAO*2,WORK,KFREE,LFREE)     
         CALL MEMGET('INTE',KNCONT,NPAO*2,WORK,KFREE,LFREE)     
         CALL MEMGET('INTE',KIORBS,NPAO  ,WORK,KFREE,LFREE)     
         CALL MEMGET('INTE',KJORBS,NPAO  ,WORK,KFREE,LFREE)     
         CALL MEMGET('INTE',KKORBS,NPAO  ,WORK,KFREE,LFREE)     
         CALL PAOVEC(WORK(KJSTRS),WORK(KNPRIM),WORK(KNCONT),
     &               WORK(KIORBS),WORK(KJORBS),WORK(KKORBS),0,
     &               .FALSE.,IPRRSP)
         CALL MEMREL('HERINT.PAOVEC',WORK,1,KJORBS,KFREE,LFREE)
         CALL TIMER('START ',TIMSTR,TIMEND)
         CALL GETTIM(SO1CPU,SO1WAL)
         I2TYP = 0
         IRNTYP = -20
         IPRTWO = 0
         TKTIME = .FALSE.
         CALL DZERO(WORK(KF),NDMAT*N2BASX)
C        Always full exchange in spin-orbit integrals:
         HFXSAV=HFXFAC
         HFXFAC=D1
         CALL TWOINT(WORK(KFREE),LFREE,WORK(KFREE),
     &               WORK(KF),WORK(KD),NDMAT,
     &               IDUMMY, IFCTYP,
     &               DUMMY,IDUMMY,NUMDIS,1,IRNTYP,0,0,.TRUE.,.TRUE.,
     &               .FALSE.,TKTIME,IPRTWO,IPRNTA,IPRNTB,IPRNTC,
     &               IPRNTD,RTNTWO,IDUMMY,I2TYP,WORK(KJSTRS),
     &               WORK(KNPRIM),WORK(KNCONT),WORK(KIORBS),
     &               IDUMMY,IDUMMY,DUMMY,DUMMY,DUMMY,DUMMY,
     &               .FALSE.)
         CALL MEMREL('HSOFCKAO.TWOINT',WORK,KF,KJSTRS,KFREE,LFREE)
         HFXFAC=HFXSAV
      IF (DEBUG) THEN
         CALL MAOPRI(WORK(KF1),'X2HSOAO(TWOINT):Fock 1')
         CALL MAOPRI(WORK(KF2),'X2HSOAO(TWOINT):Fock 2')
      END IF
C
C  Transform to SO
C
         ISYMDM(1)=0
         ISYMDM(2)=VECSYM-1
         CALL SKLFCK(WORK(KF),DUMMY,WORK(KFREE),LFREE,
     &               IPRTWO,.FALSE.,
     &               .FALSE.,.FALSE.,.FALSE.,.TRUE.,IDUMMY,.TRUE.,NDMAT,
     &               ISYMDM,IFCTYP,IDUMMY,.TRUE.)
C
         CALL MEMCHK('HSOFCKAO.SKLFCK',WORK,KF)
         CALL GETTIM(SO2CPU,SO2WAL)
         SOCPU=SO2CPU-SO1CPU
         SOWAL=SO2WAL-SO1WAL
         CALL TIMER('TWOINT',TIMSTR,TIMEND)
         CALL FLSHFO(LUPRI)
         WRITE(LUPRI,'(2(/A),2(/A,F5.0,A))')
     &      '   Two-electron spin-orbit integrals',
     &   '   =================================',
     &   '   Spin-orbit 2-electron CPU  time ',SOCPU,' seconds',
     &   '   Spin-orbit 2-electron wall time ',SOWAL,' seconds'
      ELSE
         TRIPLET(1) = GRSPIN.NE.VECSPIN
         TRIPLET(2) = GRSPIN.EQ.1
         CALL GETFMATX(OPLBL,TRIPLET,WORK(KF),WORK(KD),2)
      END IF
C
C One-electron
C
      IF (OPLBL(2:2).EQ.' ') THEN
         CALL GET_P(OPLBL(1:1)//'1'//OPLBL(3:8),WORK(KD1))
         CALL DAXPY(N2BASX,D1,WORK(KD1),1,WORK(KF1),1)
      END IF
      IF (DEBUG) THEN
         CALL MAOPRI(WORK(KF1),'X2HSOAO:Fock 1')
         CALL MAOPRI(WORK(KF2),'X2HSOAO:Fock 2')
      END IF
C
C Final transform (adds [k,f1] to f2
C
      CALL F_K(NBAST,WORK(KAO),WORK(KF1),WORK(KF2),WORK(KS),
     &         WORK(KFREE),LFREE)
      IF (DEBUG) THEN
         CALL MAOPRI(WORK(KF2),'X2HSOAO:Final fock')
      END IF
C
C Mo basis (add to one-electron in input)
C
      CALL DZERO(WORK(KF1),N2ORBX)
      CALL AOTOMO(WORK(KF2),WORK(KF1),CMO,GRSYM,WORK(KFREE),LFREE)
      IF (DEBUG) THEN
         CALL MAOPRI(WORK(KF1),'X2HSOAO:Final fock (mo)')
      END IF
C
C Final gradient
C
      TRPLET=GRSPIN.EQ.VECSPIN
      CALL ORBEX(GRSYM,1,LGR,
     &   WORK(KF1),WORK(KFREE),WORK(KFREE),WORK(KFREE),
     &   GR,D1,WORK(KFREE),MJWOP)
      IF (DEBUG) THEN
         CALL OUTPUT(GR,1,LGR/2,1,2,LGR/2,2,1,LUPRI)
      END IF
      CALL MEMREL('X2HSOAO',WORK,1,1,KFREE,LFREE)
      CALL QEXIT('X2HSOAO')
      END
C
C /* Deck get_p */
C
      SUBROUTINE GET_P(LABEL,P)
C
C Simplified PRPGET
C
#include <implicit.h>
      CHARACTER*8 LABEL
      DIMENSION P(*)
C
C External
C
#include <inforb.h>
#include <inftap.h>
#include <dummy.h>
      LOGICAL FNDLB2, CLOSED
      EXTERNAL FNDLB2
C
C Local
C
      DIMENSION TMP(N2BASX)
      CHARACTER*8 RTNLBL(2)
C
      CALL QENTER('GET_P')
C
C Leave AOPROPER in the same state (open/closed)
C
      CLOSED=LUPROP.LT.0
      IF (CLOSED) THEN
         CALL GPOPEN(LUPROP,'AOPROPER','OLD',' ','UNFORMATTED',IDUMMY,
     &   .FALSE.)
      END IF
      REWIND LUPROP
      IF ( FNDLB2(LABEL,RTNLBL,LUPROP)) THEN
         IF (RTNLBL(2).EQ.'SQUARE  ') THEN
            CALL READT(LUPROP,N2BASX,P)
         ELSE IF (RTNLBL(2).EQ.'SYMMETRI') THEN
            CALL READT(LUPROP,NNBASX,TMP)
            CALL DSPTSI(NBAST,TMP,P)
         ELSE IF (RTNLBL(2).EQ.'ANTISYMM') THEN
            CALL READT(LUPROP,NNBASX,TMP)
            CALL DAPTGE(NBAST,TMP,P)
         ELSE
            CALL QUIT('Error: No antisymmetry label on LUPROP')
         END IF
      ELSE
         CALL QUIT('GET_P: '//LABEL//' not found on LUPROP')
      END IF
      IF (CLOSED) CALL GPCLOSE(LUPROP,'KEEP')
      CALL QEXIT('GET_P')
      END
