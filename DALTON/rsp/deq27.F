!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2015 (2015), see http://daltonprogram.org"
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
!
! FILE: deq27.F
!
C
C  /* Deck deq27 */
      SUBROUTINE DEQ27(CMO,UBO,DV,DXCAO,DXVAO,WRK,LFRSAV)
C
C JAN-20 1992, H.A.
C
C Purpose:
C To compute eq.27 for the inactive density matrix and its counterpart
C for the active density matrix, following
C the direct RPA article of Jensen et al., Chem. Phys. 119, 297 (1988).
C That is to compute Inactive and Active density matrices
C to be used for calculating one-index transformed
C Fock matrices over Atomic Orbitals.
C
C It follows TR1DEN in SIRIUS, which in turn
C follows appendix C in Chem. Phys. 104 (1986) 229.
C
C DXCAO and DXVAO must be reserved as square matrices (NBAST,NBAST)
C
C General formula: DX(p,q) = sum(t) Bo(t,p) D(t,q) - D(p,t) Bo(q,t)
C Thus: DXC(p,i) =   2 Bo(i,p)
C       DXC(i,p) = - 2 Bo(p,i)
C       DXV(p,u) =   sum(v) DV(v,u) Bo(v,p)
C       DXV(u,p) = - sum(v) DV(u,v) Bo(p,v)
C
C Input:
C   CMO(*)  molecular orbital coefficients
C   UBO(*)  kappa matrix (orbital part of B), unpacked
C   DV(*)   active part of one-electron density matrix
C           (over MO's)
C
C Output:
C   if NISHT>0: DXCAO(*) = DXC in AO basis
C   if NASHT>0: DXVAO(*) = DXV in AO basis
C
C Scratch:
C   WRK(LFRSAV)
C
#include "implicit.h"
      DIMENSION CMO(*),UBO(NORBT,*),DV(NASHT,*),DXCAO(*),DXVAO(*),
     *          WRK(*)
C
      PARAMETER ( DM1 = -1.D0 )
C
C Used from common blocks:
C  INFORB : NSYM,NASHT,...
C  WRKRSP : KSYMOP
C
#include "maxash.h"
#include "maxorb.h"
#include "inforb.h"
#include "infrsp.h"
#include "wrkrsp.h"
C
#if defined (VAR_DEBUG)
#include "idbg.h"
#endif
C
      CALL QENTER('DEQ27')
C
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LFRSAV
      LDXAO1 = MAX(NASHT,NISHT)*NBAST
      CALL MEMGET2('REAL','DXAO1',KDXAO1,LDXAO1,WRK,KFREE,LFREE)
      CALL MEMGET2('REAL','DXAO2',KDXAO2,N2BASX,WRK,KFREE,LFREE)
C
C     First DXVAO, if nasht.gt.0:
C
      IF (NASHT .GT. 0) THEN
C
C     *************************************************
C     D_I and D_II:  D = D_I - D_II   Active matrices
C     *************************************************
C     DXVAO-I: loop over symmetries then put result in DXVAO
C
         CALL MEMGET2('REAL','DXV',KDXV ,NASHT*NORBT,WRK,KFREE,LFREE)
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
C        DXV_I(p.u) = sum(v) Bo(v,p) DV(v,u)
         IASHI = IASH(ISYM)
         IF (NASHI .EQ. 0 .OR. NORBJ .EQ. 0) GO TO 2000
         CALL DGEMM('T','N',NORBJ,NASHI,NASHI,1.D0,
     &              UBO(IORBI+NISHI+1,IORBJ+1),NORBT, ! (v,p)T
     &              DV(IASHI+1,IASHI+1),NASHT,0.D0,   ! (v,u)
     &              WRK(KDXV),NORBJ)                  ! (p,u)
C **     step 2 of active dm:
         CALL DGEMM('N','N',NBASJ,NASHI,NORBJ,1.D0,
     &              CMO(ICMOJ+1),NBASJ,               ! (pi,p)
     &              WRK(KDXV),NORBJ,0.D0,             ! (p,u)
     &              WRK(KDXAO1),NBASJ)                ! (pi,u)
         IOFMOV = ICMOI + 1 + NISHI*NBASI
         IDXAO2 = KDXAO2 + IBAS(JSYM)*NBAST + IBAS(ISYM)
         CALL DGEMM('N','T',NBASI,NBASJ,NASHI,1.D0,
     &              CMO(IOFMOV),NBASI,                ! (upsilon,u)
     &              WRK(KDXAO1),NBASJ,0.D0,           ! (pi,u)T
     &              WRK(IDXAO2),NBAST)                ! (pi,upsilon)
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
C        DXV_II(p,u) = sum(v) Bo(p,v) DV(v,u)
         IASHI = IASH(ISYM)
         IASHJ = IASH(JSYM)
         IF( NASHJ .EQ. 0 .OR. NORBI .EQ. 0) GO TO 3000
         CALL DGEMM('N','N',NORBI,NASHJ,NASHJ,1.D0,
     &              UBO(IORBI+1,IORBJ+NISHJ+1),NORBT,  ! (p,v)
     &              DV(IASHJ+1,IASHJ+1),NASHT,0.D0,    ! (v,u)
     &              WRK(KDXV),NORBI)                   ! (p,u)
C **     step 2 of active dm:
         CALL DGEMM('N','N',NBASI,NASHJ,NORBI,1.D0,
     &              CMO(ICMOI+1),NBASI,                ! (pi,p)
     &              WRK(KDXV),NORBI,0.D0,              ! (p,u)
     &              WRK(KDXAO1),NBASI)                 ! (pi,u)
         IOFMOV = ICMOJ + 1 + NISHJ*NBASJ
         IDXAO2 = KDXAO2 + IBAS(JSYM)*NBAST + IBAS(ISYM)
         CALL DGEMM('N','T',NBASI,NBASJ,NASHJ,1.D0,
     &              WRK(KDXAO1),NBASI,                 ! (pi,u)
     &              CMO(IOFMOV),NBASJ,0.D0,            ! (beta,u)T
     &              WRK(IDXAO2),NBAST)                 ! (pi,beta)
C **  this symmetry block finished
 3000 CONTINUE
C
C subtract D_II from  D_I
         CALL DAXPY(N2BASX,DM1,WRK(KDXAO2),1,DXVAO,1)
      END IF ! NASHT .gt. 0
C
C     DXVAO finished, now DXCAO, if nisht.gt.0:
C
      IF (NISHT .GT. 0) THEN
C
C     *************************************************
C     D_I and D_II:  D = D_I - D_II   Inactive matrices
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
C        DXC_I(i,p) = 2 Bo(i,p) -> DXCAO_I(iota,pi)

         IF (NORBJ .NE. 0) THEN
            CALL DGEMM('N','T',NBASJ,NISHI,NORBJ,2.D0,  ! 2
     &                 CMO(ICMOJ+1),NBASJ,              ! (pi,p)
     &                 UBO(IORBI+1,IORBJ+1),NORBT,0.D0, ! (i,p)T
     &                 WRK(KDXAO1),NBASJ)               ! (pi,i)
            IDXAO2 = KDXAO2 + IBAS(JSYM)*NBAST + IBAS(ISYM)
            CALL DGEMM('N','T',NBASI,NBASJ,NISHI,1.D0,
     &                 CMO(ICMOI+1),NBASI,              ! (iota,i)
     &                 WRK(KDXAO1),NBASJ,0.D0,          ! (pi,i)T
     &                 WRK(IDXAO2),NBAST)               ! (iota,pi)
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
C        DXC_II(p,i) = 2*Bo(p,i) -> DXCAO_II(pi,iota)

         IF (NORBJ .NE. 0) THEN
            CALL DGEMM('N','N',NBASJ,NISHI,NORBJ,2.D0,    ! 2
     &                 CMO(ICMOJ+1),NBASJ,                ! (pi,p)
     &                 UBO(IORBJ+1,IORBI+1),NORBT,0.D0,   ! (p,i)
     &                 WRK(KDXAO1),NBASJ)                 ! (pi,i)
            IDXAO2 = KDXAO2 + IBAS(ISYM)*NBAST + IBAS(JSYM)
            CALL DGEMM('N','T',NBASJ,NBASI,NISHI,1.D0,
     &                 WRK(KDXAO1),NBASJ,                 ! (pi,i)
     &                 CMO(ICMOI+1),NBASI,0.D0,           ! (iota,i)T
     &                 WRK(IDXAO2),NBAST)                 ! (pi,iota)
         END IF
C **  this symmetry block finished
 5000 CONTINUE
C
C do the subtraction in eq. 27.
C
         CALL DAXPY(N2BASX,DM1,WRK(KDXAO2),1,DXCAO,1)
C
      END IF ! NISHT .gt. 0
C
      CALL MEMREL('DEQ27',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('DEQ27')
      RETURN
C
C *** end of subroutine DEQ27
C
      END
C  /* Deck deq27mo */
      SUBROUTINE DEQ27MO(IBOSYM,UBO,UDV,DXC,DXV,WRK,LFRSAV)
C
C 2-Oct-2013 Hans Joergen Aa. Jensen.
C Based on DEQ27.
C The difference is that this routine returns the one-index
C transformed density matrices DXC and DXV in MO basis,
C while DEQ27 transforms them to AO-basis (DXCAO and DXVAO).
C Furthermore, we do not assume UDV is symmetric, thus DEQ27MO
C can also be used for non-symmetric transition density matrices.
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
C DXC and DXV are unpacked, i.e. (NORBT,NORBT).
C
C General formula: DX(p,q) = sum(t) Bo(t,p) D(t,q) - D(p,t) Bo(q,t)
C                          = DX-I(p,q) - DX-II(p,q)
C where DXC-I (p,i) = 2 Bo(i,p)
C       DXC-II(i,q) = 2 Bo(q,i)
C       DXV-I (p,u) = sum(v) UDV(v,u) Bo(v,p)
C       DXV-II(u,q) = sum(v) UDV(u,v) Bo(q,v)
C 
C
C Input:
C   UBO(NORBT,NORBT)  kappa matrix (orbital part of B), unpacked, symmetry IBOSYM
C   UDV(NASHT,NASHT)  active part of one-electron density matrix (over MOs), symmetry 1
C
C Output:
C   DXC(NORBT,NORBT)  one-index transformed inactive density matrix "DC"
C   DXV(NORBT,NORBT)  one-index transformed UDV
C
C Scratch:
C   WRK(LFRSAV)
C
#include "implicit.h"
      DIMENSION UBO(NORBT,*),UDV(NASHT,*),DXC(NORBT,*),DXV(NORBT,*),
     *          WRK(*)
C
C Used from common blocks:
C  INFORB : NSYM,NASHT,...
C
#include "maxash.h"
#include "maxorb.h"
#include "inforb.h"
#include "infrsp.h"
C
      CALL QENTER('DEQ27MO')
C
      IF (SOPPA) THEN
         CALL QDUMP(LUPRI)
         CALL QUIT('DEQ27MO fatal error: SOPPA not implemented yet!')
      END IF
C
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LFRSAV
C
C     First DXV, if nasht.gt.0:
C
      IF (NASHT .GT. 0) THEN
C
C     *************************************************
C     D-I and D-II:  D = D-I - D-II   Active matrices
C     *************************************************
C
      DO 2000 ISYM = 1,NSYM
         JSYM  = MULD2H(ISYM,IBOSYM)
         NISHI = NISH(ISYM)
         NISHJ = NISH(JSYM)
         IASHI = IASH(ISYM)
         NASHI = NASH(ISYM)
         IASHJ = IASH(JSYM)
         NASHJ = NASH(JSYM)
         IORBI = IORB(ISYM)
         NORBI = NORB(ISYM)
         IORBJ = IORB(JSYM)
         NORBJ = NORB(JSYM)
C *****  first calculate one-index transformation of second index
C *****  of UDV(uv), the active density matrix.
C        DXV-I(p,u) = sum(v) Bo(v,p) UDV(v,u)
         IF (NASHI .EQ. 0 .OR. NORBJ .EQ. 0) GO TO 2000
         CALL DGEMM('T','N',NORBJ,NASHI,NASHI,1.D0,
     &              UBO(IORBI+NISHI+1,IORBJ+1),NORBT,
     &              UDV(IASHI+1,IASHI+1),NASHT,1.D0,
     &              DXV(IORBJ+1,IORBI+NISHI+1),NORBT)
C
C
C        DXV-II: subtact results to DXV
C *****  first calculate one-index transformation of second index
C *****  of UDV(uv), the active density matrix.
C        DXV-II(p,u) = sum(v) Bo(v,p) UDV(v,u)
C        subtract D_II from  D-I
         IF( NASHJ .EQ. 0 .OR. NORBI .EQ. 0) GO TO 2000
         CALL DGEMM('N','N',NORBI,NASHJ,NASHJ,-1.D0,
     &              UBO(IORBI+1,IORBJ+NISHJ+1),NORBT,
     &              UDV(IASHJ+1,IASHJ+1),NASHT,1.D0,
     &              DXV(IORBJ+NISHJ+1,IORBI+1),NORBT)
 2000 CONTINUE
C
      END IF
C
C     DXV finished, now DXC, if nisht.gt.0:
C
      IF (NISHT .GT. 0) THEN
C
C     *************************************************
C     D-I and D-II:  D = D-I - D-II   Inactive matrices
C     *************************************************
C     Loop over symmetries then put result in DXC (see eq. 22)
C
      DO 4000 ISYM = 1,NSYM
         IORBI = IORB(ISYM)
         NISHI = NISH(ISYM)
         JSYM  = MULD2H(ISYM,IBOSYM)
         IORBJ = IORB(JSYM)
         NORBJ = NORB(JSYM)
C
C **     the inactive one-index transformed dm:
C
         DO I = IORBI+1,IORBI+NISHI
            DO J = IORBJ+1,IORBJ+NORBJ
C              DXC-I (p,i) = 2*Bo(i,p)
               DXC(J,I) = DXC(J,I) + 2.0D0 * UBO(I,J)
C              DXC-II(i,p) = -2*Bo(p,i)
               DXC(I,J) = DXC(I,J) - 2.0D0 * UBO(J,I)
            END DO
         END DO
C
 4000 CONTINUE
C
      END IF
C
      CALL QEXIT('DEQ27MO')
      RETURN
C
C *** end of subroutine DEQ27MO
C
      END
! -- end of file deq27.F --
