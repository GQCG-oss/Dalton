C
C...   Copyright (c) 1997 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of 
C...   "Dalton, an ab initio electronic structure program, Release 1.0
C...   (1997), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, T. Andersen, K. L. Bak, V. Bakken,
C...   O. Christiansen, P. Dahle, E. K. Dalskov, T. Enevoldsen,
C...   H. Heiberg, D. Jonsson, S. Kirpekar, R. Kobayashi, H. Koch,
C...   K. V. Mikkelsen, P. Norman, M. J. Packer, T.Saue,
C...   P. R. Taylor, and O. Vahtras"
C...
C...   This source code is provided under a written licence and may be 
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may 
C...   be distributed outside the research group of the licence holder. 
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence. 
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html 
C
C
C$Id: qfock.F,v 1.1.1.1 2001-02-08 13:33:14 hjj Exp $
C
#include <single.h>
#ifdef UNDEF
/* Comdeck log */
TO DO/941223-hjaaj:
(1)transpose matrices in QDENS[12]; new QDENS[12] as the new LFKDEN ?
(2)check what triplet means for IFCTYP and for QDENS[12].
(3)990505:Note that QDENS[12] have been replaced by CDENS[21] during
   implementation of CR.
990505-hjaaj: Use d.m. symm. in IFCTYP to save time in TWOINT
941221-hjaaj
 new parameter list (including spin values)
 implemented new SIRFCK
#endif
C  /* Deck qfock */
      SUBROUTINE QFOCK(NOSIM,MSYMB,MSYMC,IGRSPI,ISPIN1,ISPIN2,
     &                 ZYMATB,ZYMATC,FC,CMO,FDTI,WRK,KFREE,LFREE)
C
C 26-FEB 1992 : Ha  / last revision 23-Dec-1994 hjaaj
C
C PURPOSE:
C To calculate Fock matrix (FDTI) with double one-index transformed
C integrals to be used in double-direct quadratic RPA
C
C FLOW:
C 1 Call QDENS1 once to get D-1 and QDENS2 twice to get D-2,D-3
C 2 Call SIRFCK to get fock matrices F-1,F-2,F-3 in AO basis
C 3 Transform F-1 to MO basis
C 4 MO- transform, one-index transform fock matrix F-2 with kappa-2
C   and add to F-1
C 5 MO-transform F-3
C 6 One-index transform the inactive fock matrix FC with kappa-1 and
C   add result to F-3
C 7 one-index transform F-3 with kappa-1 and to add to F-1 =
C   = FDTI which is returned
C
C FDTI should be called with NOSIM=1 (but programmed here for .LT. 1)
C
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION FC(2),ZYMATB(NORBT,NORBT,2),ZYMATC(NORBT,NORBT,2)
      DIMENSION FDTI(NORBT,NORBT,2),WRK(2),CMO(2)
C
C  INFDIM : NASHDI
C  INFINP : DIRFCK
C  WRKRSP : KSYMOP
C
#include <thrzer.h>
#include <maxorb.h>
#include <maxash.h>
#include <inforb.h>
#include <infinp.h>
#include <wrkrsp.h>
#include <infpri.h>
C
C -- local constants
C
      PARAMETER ( DUMMY = 1.0D+20 )
      PARAMETER ( D1 = 1.0D0, D2 = 2.0D0 )
      PARAMETER ( MXOSIM = 20)
      DIMENSION ISYMDM(3*MXOSIM), IFCTYP(3*MXOSIM)
C
      IF (NOSIM .GT. MXOSIM) THEN
         WRITE (LUPRI,*) 'QFOCK ERROR: NOSIM .gt. '//
     &        'MXOSIM; NOSIM,MXOSIM =',NOSIM,MXOSIM
         CALL QUIT('error: NOSIM .gt. MXOSIM in QFOCK')
      END IF
C
C ALLOCATE WORK SPACE
C
      KFXCAO = 1
      KDXCAO = KFXCAO + 3 * N2BASX * NOSIM
      KWRKMU = KDXCAO + 3 * N2BASX * NOSIM
      LWRKMU = LFREE - KWRKMU
C
C     We will use WRK(KWRKMU) to keep a matrix of dimension 
      IF (LWRKMU.LT.N2ORBX) CALL ERRWRK('QFOCK',LWRKMU,LFREE)
C MO-Fock matrices
      KFMO1  = KDXCAO
      KFMO2  = KFMO1 + N2ORBX * NOSIM
C
C 3*NOSIM density matrices and 3*NOSIM Fock matrcies put consecutive
C in core.
C Construct density matrices
C
      IOFF2 = NOSIM * N2BASX
      IOFF3 = 2 * NOSIM * N2BASX
      KOFFAO = 0
      CALL DZERO(WRK(KFXCAO),(6*NOSIM*N2BASX))
      DO 100 IOSIM = 1,NOSIM
       CALL CDENS2(MSYMB,MSYMC,CMO,ZYMATB(1,1,IOSIM),ZYMATC(1,1,IOSIM),
     *             WRK(KDXCAO+KOFFAO),WRK(KFXCAO+KOFFAO),
     *	           WRK(KFXCAO+IOFF2+KOFFAO),WRK(KFXCAO+IOFF3+KOFFAO))
       CALL CDENS1(MSYMB,CMO,ZYMATB(1,1,IOSIM),WRK(KDXCAO+IOFF2+KOFFAO),
     *             WRK(KWRKMU),LWRKMU)
       CALL CDENS1(MSYMC,CMO,ZYMATC(1,1,IOSIM),WRK(KDXCAO+IOFF3+KOFFAO),
     *             WRK(KWRKMU),LWRKMU)
       KOFFAO = KOFFAO + N2BASX
 100  CONTINUE
C scale by 2
      CALL DSCAL((3*NOSIM*N2BASX),D2,WRK(KDXCAO),1)
C
C Construct Fock matrices
C
      CALL DZERO(FDTI,(NOSIM* N2ORBX))
C
C     Note: KSYMOP can be different for the different density
C     matrices. This must be corrected.
C941223-hjaaj: is this note by Poul Joergensen or Kenneth Ruud?
C
C     Symmetries of the three densities are b x c, b, and c
C
      NFMAT = 3 * NOSIM
      KSYMP = MULD2H(MSYMB,MSYMC)
      DO IOSIM = 1,NOSIM
         ISYMDM(IOSIM) = KSYMP
         ISYMDM(NOSIM+IOSIM) = MSYMB
         ISYMDM(2*NOSIM+IOSIM) = MSYMC
      END DO
C
C 941223-hjaaj: check what triplet means for IFCTYP and for QDENS[12].
C
      IF (ISPIN1 .NE. 0 .OR. ISPIN2 .NE. 0 .OR. IGRSPI .NE. 0) THEN
         CALL QUIT('QFOCK error:triplet not implemented for DIRFCK yet')
      END IF
C
C     IFCTYP = XY
C       X indicates symmetry about diagonal
C         X = 0 No symmetry
C         X = 1 Symmetric
C         X = 2 Anti-symmetric
C       Y indicates contributions
C         Y = 0 no contribution !
C         Y = 1 Coulomb
C         Y = 2 Exchange
C         Y = 3 Coulomb + Exchange
C
C     Check if density matrix is unsymmetric (IX=0),
C     symmetric (IX=10), antisymmetric (IX=20), or zero matrix (IX=30)
C     to threshold THRZER
C
      JDXCAO = KDXCAO
      DO I = 1,NFMAT
         IX = 10 * MATSYM(NBAST,NBAST,WRK(JDXCAO),THRZER)
C        INTEGER FUNCTION MATSYM(N,NDIM,AMAT,THRZER)
         IF (IX .EQ. 30) THEN
C           zero density matrix, do nothing !
            IFCTYP(I) = 0
CCCC     ELSE IF ("triplet Fock matrix") THEN
C           only exchange !triplet not implemented for DIRFCK yet
CCCC        IFCTYP(I) = IX + 2
         ELSE
C           Coulomb+exchange
            IFCTYP(I) = IX + 3
         END IF
         JDXCAO = JDXCAO + N2BASX
      END DO
      CALL DZERO(WRK(KFXCAO),NFMAT*N2BASX)
      CALL SIRFCK(WRK(KFXCAO),WRK(KDXCAO),NFMAT,ISYMDM,IFCTYP,DIRFCK,
     &            WRK(KWRKMU),LWRKMU)
C
C ***** Transform F-1 to MO basis
C
      CALL DZERO(WRK(KDXCAO),(3 * N2BASX * NOSIM))
      IF (NOSIM .GT. 0 .AND. NISHT .GT. 0) THEN
      KOFFAO = 0
      DO 400 IOSIM = 1,NOSIM
         DO 500 ISYM=1,NSYM
            KKSYM   = MULD2H(ISYM,MSYMB)
            JSYM   = MULD2H(KKSYM,MSYMC)
            NORBI  = NORB(ISYM)
            NORBJ  = NORB(JSYM)
            IF (NORBI.EQ.0 .OR. NORBJ.EQ.0) GO TO 500
C
            CALL AUTPV(ISYM,JSYM,CMO(ICMO(ISYM)+1),CMO(ICMO(JSYM)+1),
     &      WRK(KFXCAO+KOFFAO),NBAS,NBAST,FDTI(1,1,IOSIM),NORB,
     &      NORBT,WRK(KWRKMU),LWRKMU)
  500    CONTINUE
C
         KOFFAO = KOFFAO + N2BASX
  400 CONTINUE
      END IF
C
C ***** Transform F-2 to MO basis
C
      KOFFAO = 0
      KOFFMO = 0
      IF (NOSIM .GT. 0 .AND. NISHT .GT. 0) THEN
      DO 200 IOSIM = 1,NOSIM
         DO 300 ISYM=1,NSYM
            JSYM   = MULD2H(ISYM,MSYMB)
            NORBI  = NORB(ISYM)
            NORBJ  = NORB(JSYM)
            IF (NORBI.EQ.0 .OR. NORBJ.EQ.0) GO TO 300
            CALL AUTPV(ISYM,JSYM,CMO(ICMO(ISYM)+1),CMO(ICMO(JSYM)+1),
     &      WRK(KFXCAO+IOFF2+KOFFAO),NBAS,NBAST,WRK(KFMO1+KOFFMO),NORB,
     &      NORBT,WRK(KWRKMU),LWRKMU)
  300    CONTINUE
C
         KOFFMO = KOFFMO + N2ORBX
         KOFFAO = KOFFAO + N2BASX
  200 CONTINUE
      END IF
C
C One-index transform fock matrix F-2 with kappa-2
C and add to F-1
C
      CALL DZERO(WRK(KFMO2),NOSIM*N2ORBX)
      KOFFMO = 0
      DO 310 IOSIM = 1, NOSIM
        CALL OITH1(MSYMC,ZYMATC(1,1,IOSIM),WRK(KFMO1 + KOFFMO),
     &             WRK(KFMO2 + KOFFMO),MSYMB)
        KOFFMO = KOFFMO + N2ORBX
  310 CONTINUE
C
      CALL DAXPY((N2ORBX*NOSIM),D1,WRK(KFMO2),1,FDTI(1,1,1),1)
C
C ***** Transform F-3 to MO basis
C
      CALL DZERO(WRK(KFMO1),NOSIM*N2ORBX)
      KOFFMO = 0
      KOFFAO = 0
      IF (NOSIM .GT. 0 .AND. NISHT .GT. 0) THEN
      DO 600 IOSIM = 1,NOSIM
         DO 700 ISYM=1,NSYM
            JSYM   = MULD2H(ISYM,MSYMC)
            NORBI  = NORB(ISYM)
            NORBJ  = NORB(JSYM)
            IF (NORBI.EQ.0 .OR. NORBJ.EQ.0) GO TO 700
            CALL AUTPV(ISYM,JSYM,CMO(ICMO(ISYM)+1),CMO(ICMO(JSYM)+1),
     &      WRK(KFXCAO+IOFF3+KOFFAO),NBAS,NBAST,WRK(KFMO1+KOFFMO),NORB,
     &      NORBT,WRK(KWRKMU),LWRKMU)
  700    CONTINUE
C
         KOFFMO = KOFFMO + N2ORBX
         KOFFAO = KOFFAO + N2BASX
  600 CONTINUE
      END IF
C
C One-index transform usual inactive fock matrix with Kappa-2 and add to F-3
C
      KSAFE  = KSYMOP
      KSYMOP = MSYMC
      CALL FCKOIN(NOSIM,FC,DUMMY,ZYMATC,WRK(KFMO1),DUMMY)
      KSYMOP = KSAFE
C
C one-index transform F-3 with  Kappa-1
C
      CALL DZERO(WRK(KFMO2),(NOSIM*N2ORBX))
      KOFFMO = 0
      DO 710 IOSIM = 1, NOSIM
        CALL OITH1(MSYMB,ZYMATB(1,1,IOSIM),WRK(KFMO1 + KOFFMO),
     &             WRK(KFMO2 + KOFFMO),MSYMC)
        KOFFMO = KOFFMO + N2ORBX
  710 CONTINUE
C
C  add to F-1
C
      CALL DAXPY((N2ORBX*NOSIM),D1,WRK(KFMO2),1,FDTI,1)
C
C *** end of subroutine QFOCK
C
      RETURN
      END
