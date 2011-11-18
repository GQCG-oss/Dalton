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
#ifdef UNDEF
===========================================================================
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
===========================================================================
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
      DIMENSION FC(*),ZYMATB(NORBT,NORBT,*),ZYMATC(NORBT,NORBT,*)
      DIMENSION FDTI(NORBT,NORBT,*),WRK(*),CMO(*)
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
     *             WRK(KFXCAO+IOFF2+KOFFAO),WRK(KFXCAO+IOFF3+KOFFAO))
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
         ELSE IF (IX .EQ. 20) THEN
C           Only exchange if antisymmetric density matrix
            IFCTYP(I) = IX + 2
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
c/* Deck q3fock */
      SUBROUTINE Q3FOCK(
     &   VEC1, VEC2, LEN1 ,LEN2 , ISYM1,ISYM2, ISPIN1,ISPIN2, CMO,MJWOP,
     &   FCIN,FOIN,FC,FO,
     &   WORK,LWORK
     &   )
#include <implicit.h>
C
C Input: response vectors, length, symmetry and spin
C
      DIMENSION VEC1(*), VEC2(*), CMO(*)
      INTEGER LEN1, LEN2, ISYM1, ISYM2, ISPIN1, ISPIN2, MJWOP(*)
      DIMENSION FCIN(*), FOIN(*)
C

C Output: closed and open-shell MO high spin Fock matrices 
C
      DIMENSION FC(*), FO(*)
C
      INTEGER LWORK
      DIMENSION WORK(LWORK)
C
#ifdef VAR_DEBUG
#define DEBUG .TRUE.
#else
#define DEBUG .FALSE.
#endif
C
C External data
C
#include <inforb.h>
#include <maxorb.h>
#include <infinp.h>
#include <thrzer.h>
#include <priunit.h>
#include <infvar.h>
C
C External procedures
C
      INTEGER MATSYM
      EXTERNAL MATSYM
C   
C Local
C
      PARAMETER (DH = 0.5D0, D1=1.0D0, D2=2.0D0)
      INTEGER IFCTYP(6), ISYMDM(6), NDMAT
      INTEGER KD,KDC,KDC1,KDC2,KDC12,KDO,KDO1,KDO2,KDO12
      INTEGER KF,KFC,KFC1,KFC2,KFC12,KFO,KFO1,KFO2,KFO12,KTMP1,KTMP2
      INTEGER K1, K2, K1AO, K2AO, KUDV, KS, ID, IX, I, ISYM12
      INTEGER KFREE, LFREE
C
      CALL QENTER('Q3FOCK')
      IF (NASHT.GT.1 .AND. .NOT. HSROHF) 
     &   CALL QUIT('Q3FOCK called in a non-hispin context')
      IF (DEBUG) THEN
         CALL HEADER('Q3FOCK:Input VECB,VECC,CMO',-1)
         CALL OUTPUT(VEC1,1,LEN1/2,1,2,LEN1/2,2,1,LUPRI)
         CALL OUTPUT(VEC2,1,LEN2/2,1,2,LEN2/2,2,1,LUPRI)
         CALL PRORB(CMO,.FALSE.,LUPRI)
      END IF
      IF (NASHT.EQ.0) THEN
         NDMAT = 3
      ELSE
         NDMAT = 6
      END IF
      KFREE = 1
      LFREE = LWORK
      CALL MEMGET('REAL',K1AO,N2BASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',K2AO,N2BASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KS,N2BASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',K1,N2ORBX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',K2,N2ORBX,WORK,KFREE,LFREE)
C
C Unpack vectors to matrices to generate transformed densities
C
      CALL GTZYMT(1,VEC1,LEN1,ISYM1,WORK(K1),MJWOP)
      CALL GTZYMT(1,VEC2,LEN2,ISYM2,WORK(K2),MJWOP)
      IF (DEBUG) THEN
         CALL HEADER('Q3FOCK:Unpacked vectors K1,K2',-1)
         CALL OUTPUT(WORK(K1),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         CALL OUTPUT(WORK(K2),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
C Transform original Fock save in final
C
      CALL DZERO(FC,N2ORBX)
      CALL MEMGET('REAL',KTMP1,N2ORBX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTMP2,N2ORBX,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KTMP1),N2ORBX)
      CALL DZERO(WORK(KTMP2),N2ORBX)
      CALL PKSYM1(WORK(KTMP1),FCIN,NORB,NSYM,-1)
      CALL DSPTSI(NORBT,WORK(KTMP1),WORK(KTMP2))
      CALL DZERO(WORK(KTMP1),N2ORBX)
      CALL OITH1(ISYM1,WORK(K1),WORK(KTMP2),WORK(KTMP1),1)
      CALL OITH1(ISYM2,WORK(K2),WORK(KTMP1),FC,ISYM1)
      CALL DZERO(WORK(KTMP1),N2ORBX)
      CALL OITH1(ISYM2,WORK(K2),WORK(KTMP2),WORK(KTMP1),1)
      CALL OITH1(ISYM1,WORK(K1),WORK(KTMP1),FC,ISYM2)
      CALL DSCAL(N2ORBX,DH,FC,1)
      IF (DEBUG) THEN
         CALL MMOPRI(WORK(KTMP2),'FC (in)')
         CALL OUTPKB(FCIN,NORB,NSYM,1,LUPRI)
         CALL MMOPRI(FC,'12F')
      END IF

      IF (NASHT.GT.0) THEN
         CALL DZERO(FO,N2ORBX)
         CALL MEMGET('REAL',KTMP1,N2ORBX,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KTMP2,N2ORBX,WORK,KFREE,LFREE)
         CALL DZERO(WORK(KTMP1),N2ORBX)
         CALL DZERO(WORK(KTMP2),N2ORBX)
         CALL PKSYM1(WORK(KTMP1),FOIN,NORB,NSYM,-1)
         CALL DSPTSI(NORBT,WORK(KTMP1),WORK(KTMP2))
         CALL DZERO(WORK(KTMP1),N2ORBX)
         CALL OITH1(ISYM1,WORK(K1),WORK(KTMP2),WORK(KTMP1),1)
         CALL OITH1(ISYM2,WORK(K2),WORK(KTMP1),FO,ISYM1)
         CALL DZERO(WORK(KTMP1),N2ORBX)
         CALL OITH1(ISYM2,WORK(K2),WORK(KTMP2),WORK(KTMP1),1)
         CALL OITH1(ISYM1,WORK(K1),WORK(KTMP1),FO,ISYM2)
         CALL DSCAL(N2ORBX,DH,FO,1)
         IF (DEBUG) THEN
            CALL MMOPRI(WORK(KTMP2),'FO (in)')
            CALL OUTPKB(FOIN,NORB,NSYM,1,LUPRI)
            CALL MMOPRI(FO,'12F')
         END IF
      END IF
C
C Kappa matrices in AO
C
      CALL DZERO(WORK(K1AO),N2BASX)
      CALL DZERO(WORK(K2AO),N2BASX)
      CALL MOTOAO(WORK(K1),WORK(K1AO),CMO,ISYM1,WORK(KFREE),LFREE)
      CALL MOTOAO(WORK(K2),WORK(K2AO),CMO,ISYM2,WORK(KFREE),LFREE)
      IF (DEBUG) THEN
         CALL HEADER('Q3FOCK:Unpacked vectors K1,K2(AO)',-1)
         CALL OUTPUT(WORK(K1AO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL OUTPUT(WORK(K2AO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C
C AO overlap full square
C
      CALL GET_H1(WORK(KS),'OVERLAP',WORK(KFREE),LFREE)
      IF (DEBUG) THEN
         CALL HEADER('Q3FOCK:OVERLAP',-1)
         CALL OUTPUT(WORK(KS),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
      CALL MEMREL('MOTOAO',WORK,1,K1,KFREE,LFREE)
C
C Density matrix in AO basis (NOTE: HIGH SPIN)
C
      CALL MEMGET('REAL',KUDV,N2ASHX,WORK,KFREE,LFREE)
      IF (NASHT.GT.0) THEN
         CALL DUNIT(WORK(KUDV),NASHT)
      END IF
C
C Allocate AO density matrices 
C
C     D1(r,s) = <[k1,E(r,s)]>
C     D2(r,s) = <[k2,E(r,s)]>
C     D12(r,s) = 1/2 (<[k1,[k2,E(r,s)]]>  + <[k2,[k1,E(r,s)]]>)
C
C
C Sequential occupation of all density/fock matrices for TWOINT
C
      CALL MEMGET('REAL',KD,NDMAT*N2BASX,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KD),NDMAT*N2BASX)
      KDC1    = KD
      KDC2    = KDC1  + N2BASX
      KDC12   = KDC2  + N2BASX
      KDO1    = KDC12 + N2BASX
      KDO2    = KDO1  + N2BASX
      KDO12   = KDO2  + N2BASX
      CALL MEMGET('REAL',KDC,N2BASX,WORK,KFREE,LFREE)
      IF (NASHT.GT.0) THEN
         CALL MEMGET('REAL',KDO,N2BASX,WORK,KFREE,LFREE)
      ELSE
         CALL MEMGET('REAL',KDO,0,WORK,KFREE,LFREE)
      END IF
C
      CALL GTDMSO(WORK(KUDV),CMO,WORK(KDC),WORK(KDO),WORK(KFREE))
C
C For open shells have total density in 1 and negative 
C spin density (=active) in 2
C
      IF (NASHT.GT.0) THEN
         CALL DAXPY(N2BASX,D1,WORK(KDO),1,WORK(KDC),1)
         CALL DSCAL(N2BASX,-D1,WORK(KDO),1)
      END IF
      IF (DEBUG) THEN
         CALL HEADER('Q3FOCK:GTDMSO',-1)
         CALL OUTPUT(WORK(KDC),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         IF (NASHT.GT.0) THEN
            CALL OUTPUT(WORK(KDO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         END IF
      END IF
      CALL MEMCHK('FCKDEN',WORK,1)
C
C One-index transformed densities in AO
C
      CALL D_K(NBAST,WORK(K1AO),WORK(KDC), WORK(KDC1), WORK(KS),
     &         WORK(KFREE),LFREE)
      CALL D_K(NBAST,WORK(K2AO),WORK(KDC), WORK(KDC2), WORK(KS),
     &         WORK(KFREE),LFREE)
      CALL D_K(NBAST,WORK(K1AO),WORK(KDC2),WORK(KDC12),WORK(KS),
     &         WORK(KFREE),LFREE)
      CALL D_K(NBAST,WORK(K2AO),WORK(KDC1),WORK(KDC12),WORK(KS),
     &         WORK(KFREE),LFREE)
      CALL DSCAL(N2BASX,DH,WORK(KDC12),1)
      IF (NASHT.GT.0) THEN
         CALL D_K(NBAST,WORK(K1AO),WORK(KDO), WORK(KDO1), WORK(KS),
     &            WORK(KFREE),LFREE)
         CALL D_K(NBAST,WORK(K2AO),WORK(KDO), WORK(KDO2), WORK(KS),
     &     WORK(KFREE),LFREE)
         CALL D_K(NBAST,WORK(K1AO),WORK(KDO2),WORK(KDO12),WORK(KS),
     &            WORK(KFREE),LFREE)
         CALL D_K(NBAST,WORK(K2AO),WORK(KDO1),WORK(KDO12),WORK(KS),
     &            WORK(KFREE),LFREE)
         CALL DSCAL(N2BASX,DH,WORK(KDO12),1)
      END IF
      IF (DEBUG)THEN
         CALL HEADER('Q3FOCK:Inactive densities',-1)
         CALL OUTPUT(WORK(KDC1),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL OUTPUT(WORK(KDC2),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL OUTPUT(WORK(KDC12),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         IF (NASHT.GT.0) THEN
            CALL HEADER('Q3FOCK:Active densities',-1)
            CALL OUTPUT(WORK(KDO1),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
            CALL OUTPUT(WORK(KDO2),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
            CALL OUTPUT(WORK(KDO12),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         END IF
      END IF
      CALL MEMCHK('D_K',WORK,1)
C
C Sequential occupation of all density/fock matrices for TWOINT
C
      CALL MEMGET('REAL',KF,NDMAT*N2BASX,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KF),NDMAT*N2BASX)
      KFC1    = KF
      KFC2    = KFC1  + N2BASX
      KFC12   = KFC2  + N2BASX
      KFO1    = KFC12 + N2BASX
      KFO2    = KFO1  + N2BASX
      KFO12   = KFO2  + N2BASX
C
C Set fock type numbers 
C
      ID=KD
      DO I = 1,NDMAT
         IX = 10 * MATSYM(NBAST,NBAST,WORK(ID),THRZER)
         IF (IX .EQ. 30) THEN
C           zero density matrix, do nothing !
            IFCTYP(I) = 0
         ELSE IF (IX .EQ. 20) THEN
C           Only exchange if antisymmetric matrix
            IFCTYP(I) = IX + 2
         ELSE
C           Coulomb+exchange
            IFCTYP(I) = IX + 3
         END IF
         ID = ID + N2BASX
      END DO
      ISYM12=MULD2H(ISYM1,ISYM2)
      ISYMDM(1)=ISYM1
      ISYMDM(2)=ISYM2
      ISYMDM(3)=ISYM12
      IF (NASHT.GT.0) THEN
         ISYMDM(4)=ISYM1
         ISYMDM(5)=ISYM2
         ISYMDM(6)=ISYM12
      END IF
C
C Triplet densities - only exchange 
C
      IF (ISPIN1.EQ.1)               IFCTYP(1) = 10*(IFCTYP(1)/10)+2
      IF (ISPIN2.EQ.1)               IFCTYP(2) = 10*(IFCTYP(2)/10)+2
      IF (MOD(ISPIN1+ISPIN2,2).EQ.1) IFCTYP(3) = 10*(IFCTYP(3)/10)+2
      IF (ISPIN1.NE.1)               IFCTYP(4) = 10*(IFCTYP(4)/10)+2
      IF (ISPIN2.NE.1)               IFCTYP(5) = 10*(IFCTYP(5)/10)+2
      IF (MOD(ISPIN1+ISPIN2,2).NE.1) IFCTYP(6) = 10*(IFCTYP(6)/10)+2
      IF (DEBUG) THEN
         WRITE(LUPRI,'(A,8I3)')'IFCTYP ', (IFCTYP(I),I=1,NDMAT)
         WRITE(LUPRI,'(A,8I3)')'ISYMDM ', (ISYMDM(I),I=1,NDMAT)
      END IF
C
C Calculate Fock matrices
C
      CALL SIRFCK(
     &   WORK(KF),WORK(KD),NDMAT,ISYMDM,IFCTYP,DIRFCK,WORK(KFREE),LFREE
     &   )
C
C Modify fock to match high spin convention (fi+q, fa-q)
C
      IF (NASHT.GT.0) THEN
         CALL DAXPY(3*N2BASX,-D1,WORK(KFO1),1,WORK(KFC1),1)
      END IF
C
C Release densities
C
      IF (DEBUG) THEN
         CALL HEADER('Q3FOCK:Fock matrices',15)
         CALL OUTPUT(WORK(KFC1),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL OUTPUT(WORK(KFC2),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL OUTPUT(WORK(KFC12),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         IF (NASHT.GT.0) THEN
            CALL HEADER('Q3FOCK:Active Fock matrices',15)
            CALL OUTPUT(WORK(KFO1),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
            CALL OUTPUT(WORK(KFO2),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
            CALL OUTPUT(WORK(KFO12),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         END IF
      END IF
      CALL MEMREL('SIRFCK',WORK,1,KD,KFREE,LFREE)
C
C Remaining transformations add to F_12
C
C [k1,F(k2)]
C
      CALL F_K(NBAST,WORK(K1AO),WORK(KFC2),WORK(KFC12),WORK(KS),
     &         WORK(KFREE),LFREE)
      IF (NASHT.GT.0)
     &   CALL F_K(NBAST,WORK(K1AO),WORK(KFO2),WORK(KFO12),WORK(KS),
     &            WORK(KFREE),LFREE)
C
C [k2,F(k1)]
C
      CALL F_K(NBAST,WORK(K2AO),WORK(KFC1),WORK(KFC12),WORK(KS),
     &         WORK(KFREE),LFREE)
      IF (NASHT.GT.0)
     &   CALL F_K(NBAST,WORK(K2AO),WORK(KFO1),WORK(KFO12),WORK(KS),
     &            WORK(KFREE),LFREE)
C
C Add F(k1,k2)
C
      IF (DEBUG) THEN
         CALL MAOPRI(WORK(KFC12),'Q3FOCK:Final Fock C (AO)')
         IF (NASHT.GT.0)
     &      CALL MAOPRI(WORK(KFO12),'Q3FOCK:Final Fock O (AO)')
      END IF
      CALL AOTOMO(WORK(KFC12),FC,CMO,ISYM12,WORK(KFREE),LFREE)
      IF (NASHT.GT.0)
     &   CALL AOTOMO(WORK(KFO12),FO,CMO,ISYM12,WORK(KFREE),LFREE)
      IF (DEBUG) THEN
         CALL MMOPRI(FC,'Q3FOCK:Final Fock C (MO)')
         IF (NASHT.GT.0) 
     &      CALL MMOPRI(FO,'Q3FOCK:Final Fock O (MO)')
      END IF
      CALL MEMREL('Q3FOCK',WORK,1,1,KFREE,LFREE)
      CALL QEXIT('Q3FOCK')
      END
C
C /* Deck MOTOAO */
C
      SUBROUTINE MOTOAO(XMO,XAO,CMO,XSYM,WRK,LWRK)
#include <implicit.h>
#include <inforb.h>
      DIMENSION WRK(LWRK)
      CALL QENTER('MOTOAO')
      KFREE = 1
      LFREE = LWRK
      CALL MEMGET('REAL',KCX,NORBT*NBAST,WRK,KFREE,LFREE)
      CALL MOTOAO1(XMO,XAO,WRK(KCX),CMO,XSYM)
      CALL MEMREL('MOTOAO',WRK,1,1,KFREE,LFREE)
      CALL QEXIT('MOTOAO')
      END
      SUBROUTINE MOTOAO1(XMO,XAO,CX,CMO,XSYM)
C
C MO-AO Transforms XMO, add result to XAO (is not initialized to zero)
C
#include <implicit.h>
C
C Input XMO  
C       CMO
C       XSYM
C
#include <inforb.h>
      DIMENSION XMO(NORBT,NORBT), CMO(*)
      INTEGER XSYM
C
C Output XAO
C
      DIMENSION XAO(NBAST,NBAST)
C
C External: INFORB::NORBT,NBAST,ICMO,NBAS,IBAS,NORB,IORB
C
C
C Local
C
      DIMENSION CX(NBAST,NORBT)
      PARAMETER (D1=1.0D0, D0=0.0D0)
      INTEGER ISYM,JSYM,ICMOI,ICMOJ
      INTEGER NBASI,NBASJ,IBASI,IBASJ
      INTEGER NORBI,NORBJ,IORBI,IORBJ
C
      DO ISYM=1,NSYM
         JSYM=MULD2H(ISYM,XSYM)
         NBASI=NBAS(ISYM)
         NBASJ=NBAS(JSYM)
         IF (NBASI.GT.0 .AND. NBASJ.GT.0) THEN
            ICMOI=ICMO(ISYM)+1
            ICMOJ=ICMO(JSYM)+1
            IBASI=IBAS(ISYM)+1
            IBASJ=IBAS(JSYM)+1
            NORBI=NORB(ISYM)
            NORBJ=NORB(JSYM)
            IORBI=IORB(ISYM)+1
            IORBJ=IORB(JSYM)+1
            CALL DGEMM('N','N',NBASI,NORBJ,NORBI,
     &         D1,CMO(ICMOI),NBASI,
     &            XMO(IORBI,IORBJ),NORBT,
     &         D0,CX(IBASI,IORBJ),NBAST
     &         )
            CALL DGEMM('N','T',NBASI,NBASJ,NORBJ,
     &         D1,CX(IBASI,IORBJ),NBAST,
     &            CMO(ICMOJ),NBASJ,
     &         D1,XAO(IBASI,IBASJ),NBAST
     &         )
         END IF
      END DO
      END
C
C /* Deck AOTOMO */
C
      SUBROUTINE AOTOMO(XAO,XMO,CMO,XSYM,WRK,LWRK)
#include <implicit.h>
#include <inforb.h>
      DIMENSION WRK(LWRK)
      CALL QENTER('AOTOMO')
      KFREE = 1
      LFREE = LWRK
      CALL MEMGET('REAL',KCX,NORBT*NBAST,WRK,KFREE,LFREE)
      CALL AOTOMO1(XAO,XMO,WRK(KCX),CMO,XSYM)
      CALL MEMREL('AOTOMO',WRK,1,1,KFREE,LFREE)
      CALL QEXIT('AOTOMO')
      END
      SUBROUTINE AOTOMO1(XAO,XMO,CX,CMO,XSYM)
C
C MO-AO Transforms XAO, add result to XMO (is not initialized to zero)
C
#include <implicit.h>
C
C Input XMO  
C       CMO
C       XSYM
C
#include <inforb.h>
      DIMENSION XAO(NBAST,NBAST), CMO(*)
      INTEGER XSYM
C
C Output XAO
C
      DIMENSION XMO(NORBT,NORBT)
C
C External: INFORB::NORBT,NBAST,ICMO,NBAS,IBAS,NORB,IORB
C
C
C Local
C
      DIMENSION CX(NORBT,NBAST)
      PARAMETER (D1=1.0D0, D0=0.0D0)
      INTEGER ISYM,JSYM,ICMOI,ICMOJ
      INTEGER NBASI,NBASJ,IBASI,IBASJ
      INTEGER NORBI,NORBJ,IORBI,IORBJ
C
      DO ISYM=1,NSYM
         JSYM=MULD2H(ISYM,XSYM)
         NBASI=NBAS(ISYM)
         NBASJ=NBAS(JSYM)
         IF (NBASI.GT.0 .AND. NBASJ.GT.0) THEN
            ICMOI=ICMO(ISYM)+1
            ICMOJ=ICMO(JSYM)+1
            IBASI=IBAS(ISYM)+1
            IBASJ=IBAS(JSYM)+1
            NORBI=NORB(ISYM)
            NORBJ=NORB(JSYM)
            IORBI=IORB(ISYM)+1
            IORBJ=IORB(JSYM)+1
            CALL DGEMM('T','N',NORBI,NBASJ,NBASI,
     &         D1,CMO(ICMOI),NBASI,
     &            XAO(IBASI,IBASJ),NBAST,
     &         D0,CX(IORBI,IBASJ),NORBT
     &         )
            CALL DGEMM('N','N',NORBI,NORBJ,NBASJ,
     &          D1,CX(IORBI,IBASJ),NORBT,
     &             CMO(ICMOJ),NBASJ,
     &          D1,XMO(IORBI,IORBJ),NORBT
     &          )
         END IF
      END DO
      END
      SUBROUTINE D_K(N,E,D,DK,S,X,LWRK)
C
C  DK = D*S*E - E*S*D  ([k,D](T) in ao basis)
C
#include <implicit.h>
      INTEGER N
      DIMENSION E(N,N), D(N,N), DK(N,N), S(N,N)
      PARAMETER (D0=0.0D0, D1=1.0D0, DM1=-D1)
      if(LWRK.LT.N*N) CALL STOPIT('D_K','X',1,LWRK)
      CALL DGEMM('N','N',N,N,N,D1, D,N,S,N,D0,X, N)
      CALL DGEMM('N','N',N,N,N,D1, X,N,E,N,D1,DK,N)
      CALL DGEMM('N','N',N,N,N,D1, E,N,S,N,D0,X, N)
      CALL DGEMM('N','N',N,N,N,DM1,X,N,D,N,D1,DK,N)
      END
      SUBROUTINE F_K(N,E,F,FK,S,WRK,LFREE)
C
C  FK = S*E*F - F*E*S  ([k,F] in ao basis)
C
#include <implicit.h>
      INTEGER N
      DIMENSION E(N,N), F(N,N), FK(N,N), S(N,N)
      DIMENSION WRK(LFREE)
      CALL D_K(N,F,S,FK,E,WRK,LFREE) 
      END

      SUBROUTINE GET_H1(S,KEY,WRK,LWRK)
#include <implicit.h>
#include <inforb.h>
      DIMENSION S(NBAST,NBAST)
      CHARACTER*(*) KEY
      DIMENSION WRK(LWRK)
#include <priunit.h>
C
C Local (automatic)
C
      !DIMENSION SBT(NNBAST), SFT(NNBASX)
      LOGICAL FOUND
      INTEGER KFREE,LFREE
      KFREE = 1
      LFREE = LWRK
      CALL MEMGET('REAL',KSBT,NNBAST,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KSFT,NNBASX,WRK,KFREE,LFREE)
      CALL RDONEL(KEY,FOUND,WRK(KSBT),NNBAST)
      IF (.NOT.FOUND) CALL QUIT('RDONEL ERROR')
      CALL PKSYM1(WRK(KSFT),WRK(KSBT),NBAS,NSYM,-1)
      CALL DSPTSI(NBAST,WRK(KSFT),S)
      CALL MEMREL('GET_H1',WRK,1,1,KFREE,LFREE)
      END

      SUBROUTINE MAOPRI(A,TEXT)
#include <implicit.h>
#include <inforb.h>
      DIMENSION A(NBAST,NBAST)
      CHARACTER*(*) TEXT
      CALL MSQPRI(A,NBAST,TEXT)
      END

      SUBROUTINE MMOPRI(A,TEXT)
#include <implicit.h>
#include <inforb.h>
      DIMENSION A(NBAST,NBAST)
      CHARACTER*(*) TEXT
      CALL MSQPRI(A,NORBT,TEXT)
      END

      SUBROUTINE MSQPRI(A,N,TEXT)
#include <implicit.h>
      INTEGER N
      DIMENSION A(N,N)
      CHARACTER*(*) TEXT
      INTEGER INDENT
      PARAMETER (INDENT = 14)
#include <priunit.h>
      IF (TEXT .NE. ' ') CALL HEADER(TEXT,INDENT)
      CALL OUTPUT(A,1,N,1,N,N,N,1,LUPRI)
      END
