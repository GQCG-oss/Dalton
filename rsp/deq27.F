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
C  /* Deck deq27 */
      SUBROUTINE DEQ27(CMO,UBO,DV,DXCAO,DXVAO,WRK,LFRSAV)
C
C JAN-20 1992, H.A.
C
C Purpose:
C To compute eq.27 for the inactive density matrix and its counterpart
C for the active density matrix, following
C the direct RPA article of Jensen et al C CP,119, 297 (1988).
C That is to compute Inactive and Active  density matrices
C later to be used for calculating one-index transformed
C anti-symmetric Fock matrices over Atomic Orbitals.
C
C It follows TR1DEN in SIRIUS, which in turn
C follows appendix C in Chem. Phys. 104 (1986) 229.
C
C DXCAO and DXVAO must be reserved as square (unfolded) matrices
C
C General formula: DX(p,q) = sum(t) Bo(t,p) D(t,q)
C Thus: DXC(p,i) = 2*Bo(i,p)
C       DXV(p,u) = sum(v) DV(v,u) Bo(v,p) = sum(v) DV(u,v) Bo(v,p)
C
C Input:
C   CMO(*)  molecular orbital coefficients
C   UBO(*)  kappa matrix (orbital part of B), unpacked
C   DV(*)   active part of one-electron density matrix
C           (over MO's)
C
C Scratch:
C   WRK(LFRSAV)
C
#include <implicit.h>
      DIMENSION CMO(*),UBO(NORBT,*),DV(NASHT,*),DXCAO(*),DXVAO(*),
     *          WRK(*)
C
      PARAMETER (HALF = 0.5D0, D2 = 2.0D0, D4 = 4.0D0)
      PARAMETER (HALFM = -0.5D0, D1 = 1.D0, DM1 = -1.D0)
C
C Used from common blocks:
C  INFORB : NSYM,NASHT,...
C  INFIND : IROW(*)
C  INFVAR : ??????
C  WRKRSP : KSYMOP
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infind.h>
#include <infdim.h>
#include <infrsp.h>
#include <infvar.h>
#include <wrkrsp.h>
C
#if defined (VAR_DEBUG)
#include <idbg.h>
#endif
C
      CALL QENTER('DEQ27')
C
C     MAERKE 931116; mjp+hjaaj
      IF (SOPPA) THEN
         CALL QDUMP(6)
         CALL QUIT('DEQ27 fatal error: SOPPA not implemented yet!')
      END IF
C
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LFRSAV
      LDXAO1 = MAX(NASHT,NISHT)*NBAST
      CALL MEMGET('REAL',KDXAO1,LDXAO1,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KDXAO2,N2BASX,WRK,KFREE,LFREE)
C
C     First DXVAO, if nasht.gt.0:
C
      IF (NASHT .GT. 0) THEN
C
C     *************************************************
C     D-I and D-II:  D = D-I - D-II   Active matrices
C     *************************************************
C     DXVAO-I: loop over symmetries then put result in DXVAO
C
         CALL MEMGET('REAL',KDXV ,NASHT*NORBT,WRK,KFREE,LFREE)
         CALL DZERO(WRK(KDXAO2),N2BASX)
      DO 2000 ISYM = 1,NSYM
         NASHI = NASH(ISYM)
         NISHI = NISH(ISYM)
         JSYM  = MULD2H(ISYM,KSYMOP)
         ICMOI = ICMO(ISYM)
         IORBI = IORB(ISYM)
         NORBI = NORB(ISYM)
         NBASI = NBAS(ISYM)
         ICMOJ = ICMO(JSYM)
         IORBJ = IORB(JSYM)
         NORBJ = NORB(JSYM)
         NBASJ = NBAS(JSYM)
C **     step 1 of active dm:
C *****  first calculate one-index transformation of second index
C *****  of DV(uv), the active density matrix.
C        DXV(p,u) = sum(v) Bo(v,p) DV(v,u)
         IASHI = IASH(ISYM)
         IF (NASHI .EQ. 0 .OR. NORBJ .EQ. 0) GO TO 2000
         CALL DGEMM('T','N',NORBJ,NASHI,NASHI,1.D0,
     &              UBO(IORBI+NISHI+1,IORBJ+1),NORBT,
     &              DV(IASHI+1,IASHI+1),NASHT,0.D0,
     &              WRK(KDXV),NORBJ)
C **     step 2 of active dm:
         CALL DGEMM('N','N',NBASJ,NASHI,NORBJ,1.D0,
     &              CMO(ICMOJ+1),NBASJ,
     &              WRK(KDXV),NORBJ,0.D0,
     &              WRK(KDXAO1),NBASJ)
         IOFMOV = ICMOI + 1 + NISHI*NBASI
         IDXAO2 = KDXAO2 + IBAS(JSYM)*NBAST + IBAS(ISYM)
         CALL DGEMM('N','T',NBASI,NBASJ,NASHI,1.D0,
     &              CMO(IOFMOV),NBASI,
     &              WRK(KDXAO1),NBASJ,0.D0,
     &              WRK(IDXAO2),NBAST)
C **  this symmetry block finished
 2000 CONTINUE
C
         CALL DCOPY(N2BASX,WRK(KDXAO2),1,DXVAO,1)
C
C     DXVAO-II: loop over symmetries then add results to DXVAO
C               only chnage from previous loop is that first
C               MPATB is here MPAB
C
         CALL DZERO(WRK(KDXAO2),N2BASX)
      DO 3000 ISYM = 1,NSYM
         NASHI = NASH(ISYM)
         NISHI = NISH(ISYM)
         JSYM  = MULD2H(ISYM,KSYMOP)
         ICMOI = ICMO(ISYM)
         IORBI = IORB(ISYM)
         NORBI = NORB(ISYM)
         NBASI = NBAS(ISYM)
         NORBJ = NORB(JSYM)
         ICMOJ = ICMO(JSYM)
         IORBJ = IORB(JSYM)
         NORBJ = NORB(JSYM)
         NBASJ = NBAS(JSYM)
         NASHJ = NASH(JSYM)
         NISHJ = NISH(JSYM)
C **     step 1 of active dm:
C *****  first calculate one-index transformation of second index
C *****  of DV(uv), the active density matrix.
C        DXV(p,u) = sum(v) Bo(v,p) DV(v,u)
         IASHI = IASH(ISYM)
         IASHJ = IASH(JSYM)
         IF( NASHJ .EQ. 0 .OR. NORBI .EQ. 0) GO TO 3000
         CALL DGEMM('N','N',NORBI,NASHJ,NASHJ,1.D0,
     &              UBO(IORBI+1,IORBJ+NISHJ+1),NORBT,
     &              DV(IASHJ+1,IASHJ+1),NASHT,0.D0,
     &              WRK(KDXV),NORBI)
C **     step 2 of active dm:
         CALL DGEMM('N','N',NBASI,NASHJ,NORBI,1.D0,
     &              CMO(ICMOI+1),NBASI,
     &              WRK(KDXV),NORBI,0.D0,
     &              WRK(KDXAO1),NBASI)
         IOFMOV = ICMOJ + 1 + NISHJ*NBASJ
         IDXAO2 = KDXAO2 + IBAS(JSYM)*NBAST + IBAS(ISYM)
         CALL DGEMM('N','T',NBASI,NBASJ,NASHJ,1.D0,
     &              WRK(KDXAO1),NBASI,
     &              CMO(IOFMOV),NBASJ,0.D0,
     &              WRK(IDXAO2),NBAST)
C **  this symmetry block finished
 3000 CONTINUE
C
C subtract D_II from  D-I
         CALL DAXPY(N2BASX,DM1,WRK(KDXAO2),1,DXVAO,1)
      END IF
C
C     DXVAO finished, now DXCAO, if nisht.gt.0:
C
      IF (NISHT .GT. 0) THEN
C
C     *************************************************
C     D-I and D-II:  D = D-I - D-II   Inactive matrices
C     *************************************************
C     DXCAO-I: loop over symmetries then put result in DXCAO
C     see eq. 22
         CALL DZERO(WRK(KDXAO2),N2BASX)
      DO 4000 ISYM = 1,NSYM
         NISHI = NISH(ISYM)
      IF (NISHI .EQ. 0) GO TO 4000
         JSYM  = MULD2H(ISYM,KSYMOP)
         ICMOI = ICMO(ISYM)
         IORBI = IORB(ISYM)
         NORBI = NORB(ISYM)
         NBASI = NBAS(ISYM)
         ICMOJ = ICMO(JSYM)
         IORBJ = IORB(JSYM)
         NORBJ = NORB(JSYM)
         NBASJ = NBAS(JSYM)
C
C **     the inactive one-index transformed dm:
C        DXC(p,i) = 2*Bo(i,p)
         IF (NBASI*NBASJ*NORBJ .NE. 0) THEN
            CALL DGEMM('N','T',NBASJ,NISHI,NORBJ,1.D0,
     &                 CMO(ICMOJ+1),NBASJ,
     &                 UBO(IORBI+1,IORBJ+1),NORBT,0.D0,
     &                 WRK(KDXAO1),NBASJ)
            IDXAO2 = KDXAO2 + IBAS(JSYM)*NBAST + IBAS(ISYM)
            CALL DGEMM('N','T',NBASI,NBASJ,NISHI,1.D0,
     &                 CMO(ICMOI+1),NBASI,
     &                 WRK(KDXAO1),NBASJ,0.D0,
     &                 WRK(IDXAO2),NBAST)
         END IF
C **  this symmetry block finished
 4000 CONTINUE
         CALL DCOPY(N2BASX,WRK(KDXAO2),1,DXCAO,1)
C
C
C     DXCAO-II: loop over symmetries then subtract result from DXCAO
C     see eq. 25
C
         CALL DZERO(WRK(KDXAO2),N2BASX)
      DO 5000 ISYM = 1,NSYM
         NISHI = NISH(ISYM)
      IF (NISHI .EQ. 0) GO TO 5000
         JSYM  = MULD2H(ISYM,KSYMOP)
         ICMOI = ICMO(ISYM)
         IORBI = IORB(ISYM)
         NORBI = NORB(ISYM)
         NBASI = NBAS(ISYM)
         ICMOJ = ICMO(JSYM)
         IORBJ = IORB(JSYM)
         NORBJ = NORB(JSYM)
         NBASJ = NBAS(JSYM)
C
C **     the inactive one-index transformed dm:
C        UBO not transposed for D-II (see eq. 25)
C
C        DXC(p,i) = 2*Bo(i,p)
         IF (NBASI*NBASJ*NORBJ .NE. 0) THEN
            CALL DGEMM('N','N',NBASJ,NISHI,NORBJ,1.D0,
     &                 CMO(ICMOJ+1),NBASJ,
     &                 UBO(IORBJ+1,IORBI+1),NORBT,0.D0,
     &                 WRK(KDXAO1),NBASJ)
            IDXAO2 = KDXAO2 + IBAS(ISYM)*NBAST + IBAS(JSYM)
            CALL DGEMM('N','T',NBASJ,NBASI,NISHI,1.D0,
     &                 WRK(KDXAO1),NBASJ,
     &                 CMO(ICMOI+1),NBASI,0.D0,
     &                 WRK(IDXAO2),NBAST)
         END IF
C **  this symmetry block finished
 5000 CONTINUE
C
C do the subtraction in eq.27.
C
         CALL DAXPY(N2BASX,DM1,WRK(KDXAO2),1,DXCAO,1)
C
      END IF
C
      CALL MEMREL('DEQ27',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('DEQ27')
      RETURN
C
C *** end of subroutine DEQ27
C
      END
