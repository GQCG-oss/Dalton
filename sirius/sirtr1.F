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
c900115: TODO goer noget ved PRINRM
c        CALL PRINRM(DTV(1,ICSIM),NNASHX,LUPRI)
C  /* Deck sirtr1 */
      SUBROUTINE SIRTR1(NCSIM1,NOSIM,CMO,
     &                  BCVECS,CREF,INDXCI,DTV,FTV,FTQ,
     &                  UBO,DV,PV,FC,FV,FXC,FXV,FXQ,EMYX,H2XAC,
     &                  WRK,KFRSAV,LFRSAV, ORBLIN,ADDTR1)
C
C 29-Nov-1989 Hans Joergen Aa. Jensen
C
C Master routine for LINTRN one-index transformations and
C transition Fock matrices.
C
C PARAMETERS:
C
C   input : ORBLIN,NCSIM1,NOSIM,CMO (always needed)
C           ADDTR1 true : Add FTV,FTQ,FXC,FXV,FXQ,H2XAC to previous
C                         contents
C                  false: Initialize matrices to zero
C   if (.NOT.ORBLIN) then
C     input : BCVECS,CREF,INDXCI
C     output: DTV,FTV,FTQ
C   endif
C   if (NOSIM .gt. 0) then
C     input : UBO,DV,PV,FC,FV
C     output: FXC,FXV,FXQ
C     if (.NOT.ORBLIN) then
C        output: EMYX,H2XAC
C     endif
C   endif
C   scratch: WRK(KFRSAV:LFRSAV)
C
C
#include <implicit.h>
      DIMENSION CMO(*)
      DIMENSION BCVECS(NCONDI,*), CREF(*), INDXCI(*),
     &          DTV(NNASHX,*), FTV(N2ORBX,*), FTQ(NASHT*NORBT,*)
      DIMENSION UBO(N2ORBX,*), DV(*), PV(*), FC(*), FV(*),
     &          FXC(N2ORBX,*), FXV(N2ORBX,*), FXQ(NASHT*NORBT,*),
     &          EMYX(*), H2XAC(NNASHX*NNASHX,*)
      DIMENSION WRK(LFRSAV)
      LOGICAL   ADDTR1, ORBLIN
C
C *** local constants
C
      PARAMETER ( DP5 = 0.5D0 )
      LOGICAL MOFCFV
      PARAMETER ( MOFCFV = .FALSE. )
C
C Used from common blocks:
C   INFINP : FLAG(*)
C   INFORB : NNASHX,NASHT,NORBT,...
C   INFVAR : NWOPT,JWOPSY,JWOP(2,*)
C   INFLIN : NCONRF,IPRLIN
C   INFDIM : NCONDI
C   INFTRA : USEDRC
C   INFPRI : P6FLAG(*), ?
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <infvar.h>
#include <inflin.h>
#include <infdim.h>
#include <inftra.h>
#include <infpri.h>
C
C
      CALL QENTER('SIRTR1')
C
C *** Check input
C
      IF (NWOPT .NE. NWOPDI .AND. NOSIM .GT. 1) THEN
C        SIRIUS program error, correct code.
         WRITE (LUPRI,'(//2A,2I8)') ' >>> FATAL ERROR (SIRTR1)',
     *      ' dimension error; NWOPT and NWOPDI =',NWOPT,NWOPDI
         CALL QTRACE(LUPRI)
         CALL QUIT('FATAL ERROR (SIRTR1) NWOPDI DIMENSION ERROR')
      END IF
C     ... MOFCFV always false in this version (920514-hjaaj)
      IF (MOFCFV .AND. .NOT.USEDRC) THEN
         WRITE (LUPRI,'(//2A)') ' >>> FATAL ERROR (SIRTR1)',
     &      ' MOFCFV true and USEDRC false'
         CALL QTRACE(LUPRI)
         CALL QUIT('FATAL ERROR (SIRTR1) MOFCFV true and USEDRC false')
      END IF
C
C *** Set parameters for ORBLIN or LINTRN (i.e. full Hessian)
C
C   ORBLIN flags only sigma vectors of orbital Hessian wanted
C     (frozen CI coefficients). Then EMYX and H2XAC
C     will not be needed and are not calculated.
C
      IF (ORBLIN) THEN
         IF (NCSIM1 .GT. 0) THEN
            WRITE (LUPRI,'(/A)')
     &         ' SIRTR1 error, ORBLIN true and NCSIM .gt. 0'
            WRITE (LUPRI,*) 'NCSIM ',NCSIM1
            CALL QTRACE(LUPRI)
            CALL QUIT('SIRTR1 error, ORBLIN true and NCSIM .gt. 0')
         END IF
         NCSIM = -1
         MCSIM = 0
         MOSIM = 0
      ELSE
         NCSIM = NCSIM1
         MCSIM = NCSIM1
         MOSIM = NOSIM
      END IF
C
C *** Allocate work area for MAKTDM,TR1FCK,TR1H2*
C
      KPTV   = KFRSAV
      KWRK1  = KPTV   + MCSIM*NNASHX*NNASHX
      LWRK1  = LFRSAV - KWRK1
      IF (KWRK1.GT.LFRSAV) CALL ERRWRK('SIRTR1',-KWRK1,LFRSAV)
C
C
C *** Calculate the NCSIM transition DV and PV density matrices
C
      IF (NCSIM.GT.0) THEN
         CALL DZERO( DTV, NCSIM*NNASHX )
         CALL DZERO( WRK(KPTV), NCSIM*NNASHX*NNASHX )
C
         CALL MAKTDM(NCSIM,BCVECS,CREF,DTV,WRK(KPTV),
     *               INDXCI,WRK(KWRK1),LWRK1)
C        CALL MAKTDM(NCSIM,BCVECS,CREF,DTV,PTV,INDXCI,WRK,LFREE)
C
         IF (P6FLAG(20) .OR. P6FLAG(27) .OR. IPRLIN .GE. 20) THEN
            DO 110 ICSIM = 1,NCSIM
               IF (P6FLAG(20) .OR. IPRLIN .GE. 20) THEN
                  WRITE (LUPRI,1010) ICSIM,NCSIM
CTODO             CALL PRINRM(DTV(1,ICSIM),NNASHX,LUPRI)
                  CALL OUTPAK(DTV(1,ICSIM),NASHT,1,LUPRI)
               END IF
               IF (P6FLAG(27) .OR. IPRLIN .GE. 20) THEN
                  JPTV = KPTV + (ICSIM-1)*NNASHX*NNASHX
                  WRITE (LUPRI,1020) ICSIM,NCSIM
CTODO             CALL PRINRM(WRK(JPTV),NNASHX*NNASHX,LUPRI)
                  IF (P6FLAG(27) .OR. IPRLIN .GE. 27)
     &               CALL OUTPUT(WRK(JPTV),
     &                  1,NNASHX,1,NNASHX,NNASHX,NNASHX,-1,LUPRI)
               END IF
  110       CONTINUE
         END IF
      END IF
 1010 FORMAT (/' DTV transformed one-el. d.m. (no.',I3,' of',I3,')')
 1020 FORMAT (/' PTV transformed two-el. d.m. (no.',I3,' of',I3,')')
C
C
C *** Print one-index transformation matrices, if desired
C
      IF (P6FLAG(19) .OR. IPRLIN .GE. 19) THEN
         DO 210 IOSIM = 1,NOSIM
            WRITE (LUPRI,2110) IOSIM,NOSIM
            CALL OUTPUT(UBO(1,IOSIM),1,NORBT,1,NORBT,NORBT,NORBT,
     &                  1,LUPRI)
  210    CONTINUE
      END IF
 2110 FORMAT (/,' SIRTR1: one-index transformation matrix',
     *  ' (no.',I3,' of',I3,')')
C
C     Initialization of output matrices
C
      IF (ADDTR1) THEN
         IF (NASHT .GT. 1 .AND. MOSIM .GT. 0) THEN
            CALL DSCAL((MOSIM*NNASHX*NNASHX),DP5,H2XAC,1)
C           ... divide by two because H2XAC and its transposed
C               are added in TR1H2M
         END IF
      ELSE
         CALL DZERO(FXC,NOSIM*N2ORBX)
         CALL DZERO(FXV,NOSIM*N2ORBX)
         CALL DZERO(FXQ,(NOSIM*NASHT*NORBT))
         CALL DZERO(H2XAC,(MOSIM*NNASHX*NNASHX))
         LFTV = NCSIM*N2ORBX
         IF (DOMCSRDFT) LFTV = 2*LFTV
         CALL DZERO(FTV,LFTV)
         CALL DZERO(FTQ,(NCSIM*NASHT*NORBT))
      END IF
C
C
C ***  Calculate contributions from MO two-electron integrals
C
C
      IF (MOFCFV .OR. NASHT .GT. 1 .AND. .NOT.HSROHF) THEN
         CALL TR1H2M(NCSIM,NOSIM,FTQ,H2XAC,FXQ,
     *               DTV,WRK(KPTV),DV,PV,UBO,
     *               FTV,FXC,FXV,MOFCFV,.NOT.USEDRC,WRK(KWRK1),LWRK1)
         IF (USEDRC) THEN
            CALL TR1H2D(NCSIM,NOSIM,DTV,DV,PV,UBO,
     &                  FTV,FXC,FXV,FXQ,MOFCFV,WRK(KWRK1),LWRK1)
         END IF
C
C        CALL TR1H2M(NCSIM,NOSIM,FTQ,H2XAC,FXQ,DTV,PTV,DV,PV,UBO,
C    *               FTV,FXC,FXV,DOFCFV,NODRC,WRK,LFRSAV)
      END IF
C
C
C *** Calculate various transformed Fock matrices from AO integrals
C     and untransformed Fock matrices FC and FV.
C
      IF (.NOT.MOFCFV) THEN
         CALL TR1FCK(NCSIM,NOSIM,FXC,FXV,FTV,UBO,
     *               CMO, DV,DTV, WRK(KWRK1), LWRK1)
      END IF
C
C     Print FT* and FX*  matrices if requested.
C
      DO 1100 ICSIM = 1,NCSIM
         IF (P6FLAG(21) .OR. IPRLIN.GE.21) THEN
            WRITE (LUPRI,4020) ICSIM,NCSIM
            CALL OUTPUT(FTV(1,ICSIM),1,NORBT,1,NORBT,
     &                  NORBT,NORBT,1,LUPRI)
            write(lupri,*) 'hjaaj: domcsrdft =',domcsrdft
            call flshfo(lupri)
            IF (DOMCSRDFT) THEN
               WRITE (LUPRI,4021) ICSIM,NCSIM
               CALL OUTPUT(FTV(1,NCSIM+ICSIM),1,NORBT,1,NORBT,
     &                  NORBT,NORBT,1,LUPRI)
            call flshfo(lupri)
            END IF
         END IF
         IF (P6FLAG(33) .OR. IPRLIN.GE.21) THEN
            WRITE (LUPRI,3010) ICSIM,NCSIM
            CALL OUTPUT(FTQ(1,ICSIM),1,NORBT,1,NASHT,NORBT,NASHT,
     &                     1,LUPRI)
         END IF
 1100 CONTINUE
 4020 FORMAT(/' FTV matrix (no.',I3,' of',I3,')')
 4021 FORMAT(/' FTC matrix for MC SRDFT (no.',I3,' of',I3,')')
 3010 FORMAT(/' FTQ contribution to transition Fock matrix (no.',
     *  I3,' of',I3,')')
C
      DO 1200 IOSIM = 1,NOSIM
         IF (P6FLAG(33) .OR. IPRLIN.GE.21 .AND. NASHT.GT.0) THEN
            WRITE (LUPRI,3011) IOSIM,NOSIM
            CALL OUTPUT(FXQ(1,IOSIM),1,NORBT,1,NASHT,NORBT,NASHT,
     &                  1,LUPRI)
         END IF
         IF (P6FLAG(21) .OR. IPRLIN.GE.21) THEN
            WRITE (LUPRI,5020) IOSIM,NOSIM
            CALL OUTPUT(FXC(1,IOSIM),1,NORBT,1,NORBT,
     &                  NORBT,NORBT,1,LUPRI)
            IF (NASHT.GT.0) THEN
               WRITE (LUPRI,5030) IOSIM,NOSIM
               CALL OUTPUT(FXV(1,IOSIM),1,NORBT,1,NORBT,
     &                     NORBT,NORBT,1,LUPRI)
            END IF
         END IF
 1200 CONTINUE
 3011 FORMAT(/' FXQ contribution to one-index transformed Fock matrix ',
     *  '(no.',I3,' of',I3,')')
 5020 FORMAT (/' FXC matrix without FC transf. (no.',I3,' of',I3,')')
 5030 FORMAT (/' FXV matrix without FV transf. (no.',I3,' of',I3,')')
C
C
      IF (NOSIM .GT. 0) THEN
         CALL TR1FD(NOSIM,FXC,FXV,UBO,FC,FV,WRK(KWRK1),LWRK1)
      END IF
C
C
C *** Calculate EMYX(IOSIM), the inactive energy based on the one-
C     index transformed Hamiltonian.
C
      IF (.NOT. ORBLIN .AND. NOSIM .GT. 0) CALL TR1EMY(NOSIM,FXC,EMYX)
C
C ***
C
      CALL QEXIT('SIRTR1')
      RETURN
      END
C  /* Deck tr1emy */
      SUBROUTINE TR1EMY(NOSIM,FXC,EMYX)
C
C Jan 90 hjaaj
C
C *** Calculate EMYX(IOSIM), the inactive energy based on the one-
C     index transformed Hamiltonian.
C     Multiply with EMYX(IOSIM) with two, the value of DC(ii)
C
#include <implicit.h>
      DIMENSION FXC(NORBT,NORBT,NOSIM), EMYX(NOSIM)
      PARAMETER ( D0 = 0.0D0, D2 = 2.0D0 )
C
C Used from common blocks:
C   INFORB : NORBT
C   INFIND : ISX(:)
C   INFVAR : JWOPSY
C   INFLIN : IPRLIN
C   INFPRI : P6FLAG(11)
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infind.h>
#include <infvar.h>
#include <inflin.h>
#include <infpri.h>
C
      DO 300 IOSIM = 1,NOSIM
         EMYX(IOSIM) = D0
#ifdef TODO
CHJAaJ TODO June 2001: the calculated EMYX is wrong;
C only half of the one-index transformed contribution is
C included in FXC (error found by Joern Thyssen).
C Luckily it doesn't matter because the CREF component
C is projected out of the sigma vector in the current
C implementation !!!
C
C The correct is EMYX = sum(K) (HX1(K,K) + FXC(K,K))
C        and not EMYX = sum(K) 2*FXC(K,K)
C
         IF (JWOPSY.EQ.1) THEN
            DO 200 KW = 1,NISHT
               K = ISX(KW)
               EMYX(IOSIM) = EMYX(IOSIM) + FXC(K,K,IOSIM)
  200       CONTINUE
            EMYX(IOSIM) = EMYX(IOSIM)*D2
            IF (P6FLAG(11) .OR. IPRLIN.GE.11)
     &         WRITE (LUPRI,5010) IOSIM,EMYX(IOSIM)
 5010 FORMAT(/' EMYX(',I2,') = ',F25.15)
         END IF
#endif
  300 CONTINUE
      RETURN
      END
C  /* Deck tr1den */
      SUBROUTINE TR1DEN(CMO,UBO,DV,DXCAO,DXVAO,WRK,LFRSAV)
C
C Written Nov-1983 by Hans Jorgen Aa. Jensen in Uppsala, Sweden.
C Completely rewritten Nov. 1989 hjaaj
C
C Revisions
C   2-May-1983 hjaaj.
C  30-Jan-1985 hjaaj (use UTBO instead of UBO)
C   3-Nov-1989 hjaaj (use UBO(NORBT,NORBT) instead of UTBO(N2ORBT)
C
C This subroutine is called by TR1SUP to calculate the one-index
C transformed inactive and active one-electron density matrices
C over AO's, contravariant transformation.
C Those d.m.'s are needed to calculate the contribution to
C the one-index transformed inactive and active Fock matrices
C from the atomic supermatrix integrals .
C The implementation follows appendix C in
C Chem. Phys. 104 (1986) 229.
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
      DIMENSION CMO(*),UBO(NORBT,*),DV(*)
      DIMENSION DXCAO(NBAST,*),DXVAO(NBAST,*), WRK(*)
C
      PARAMETER (HALF = 0.5D0, D1 = 1.0D0, D2 = 2.0D0)
C
C Used from common blocks:
C  INFORB : NSYM,NASHT,...
C  INFVAR : JWOPSY
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infvar.h>
C
#if defined (VAR_DEBUG)
#include <idbg.h>
#endif
C
      CALL QENTER('TR1DEN')
C
C
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LFRSAV
      CALL MEMGET('REAL',KDXAO1,N2BASX,WRK,KFREE,LFREE)
C
C     ************************************
C     DXVAO: loop over symmetries
C
      IF (NASHT .GT. 0) THEN
         CALL MEMGET('REAL',KUDV ,N2ASHX,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KDXV1,NASHT*NBAST,WRK,KFREE,LFREE)
         CALL DSPTSI(NASHT,DV,WRK(KUDV))
         CALL DZERO(DXVAO,N2BASX)
      DO 2000 ISYM = 1,NSYM
         NASHI = NASH(ISYM)
      IF (NASHI .EQ. 0) GO TO 2000
         NISHI = NISH(ISYM)
         JSYM  = MULD2H(ISYM,JWOPSY)
         ICMOI = ICMO(ISYM)
         IORBI = IORB(ISYM)
         IASHI = IASH(ISYM)
         NORBI = NORB(ISYM)
         NBASI = NBAS(ISYM)
         ICMOJ = ICMO(JSYM)
         IORBJ = IORB(JSYM)
         NORBJ = NORB(JSYM)
         NBASJ = NBAS(JSYM)
C *****  calculate one-index transformation of second index
C *****  of DV(uv), the active density matrix and backtransform
C *****  to AO basis
C        (950201-hjaaj: summations reordered to avoid NASHI
C         in inner loop in order to improve vectorization)
C        DXV(p,u) = sum(v) DV(v,u) Bo(v,p) = sum(v) DV(u,v) Bo(v,p)
C        DXVAO(a,b) =
C        sum(u) CMO(b,u) ( sum(v) (sum(p) CMO(a,p) Bo(v,p)) DV(v,u) )
         IF (NBASJ*NASHI*NORBJ .NE. 0) THEN
            CALL DGEMM('N','T',NBASJ,NASHI,NORBJ,1.D0,
     &                 CMO(ICMOJ+1),NBASJ,
     &                 UBO(IORBI+NISHI+1,IORBJ+1),NORBT,0.D0,
     &                 WRK(KDXAO1),NBASJ)
            CALL DGEMM('N','N',NBASJ,NASHI,NASHI,1.D0,
     &                 WRK(KDXAO1),NBASJ,
     &                 WRK(KUDV+(NASHT+1)*IASHI),NASHT,0.D0,
     &                 WRK(KDXV1),NBASJ)
         END IF
         IOFMOV = ICMOI + 1 + NISHI*NBASI
C
         IF (NBASJ*NBASI*NASHI .NE. 0)
     &        CALL DGEMM('N','T',NBASJ,NBASI,NASHI,1.D0,
     &                   WRK(KDXV1),NBASJ,
     &                   CMO(IOFMOV),NBASI,0.D0,
     &                   DXVAO(1+IBAS(JSYM),1+IBAS(ISYM)),NBAST)
C **  this symmetry block finished
 2000 CONTINUE
C **  add transposed matrix to get final result
         CALL MTRSP(NBAST,NBAST,DXVAO,NBAST,WRK(KDXAO1),NBAST)
         CALL DAXPY(N2BASX,D1,WRK(KDXAO1),1,DXVAO,1)
      END IF
C
C
C     ************************************
C     DXCAO: loop over symmetries
C
      IF (NISHT .GT. 0) THEN
         CALL DZERO(DXCAO,N2BASX)
      DO 4000 ISYM = 1,NSYM
         NISHI = NISH(ISYM)
      IF (NISHI .EQ. 0) GO TO 4000
         JSYM  = MULD2H(ISYM,JWOPSY)          
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
         IF (NBASJ*NORBJ*NBASI .NE. 0) THEN
            CALL DGEMM('N','T',NBASJ,NISHI,NORBJ,1.D0,
     &                 CMO(ICMOJ+1),NBASJ,
     &                 UBO(IORBI+1,IORBJ+1),NORBT,0.D0,
     &                 WRK(KDXAO1),NBASJ)
            CALL DGEMM('N','T',NBASJ,NBASI,NISHI,1.D0,
     &                 WRK(KDXAO1),NBASJ,
     &                 CMO(ICMOI+1),NBASI,0.D0,
     &                 DXCAO(1+IBAS(JSYM),1+IBAS(ISYM)),NBAST)
         END IF
C **  this symmetry block finished
 4000 CONTINUE
C     
C **  add transposed matrix to get final result
         CALL MTRSP(NBAST,NBAST,DXCAO,NBAST,WRK(KDXAO1),NBAST)
         CALL DAXPY(N2BASX,D1,WRK(KDXAO1),1,DXCAO,1)
C **     multiply with DC(i,i) = 2.0 to get final result
         CALL DSCAL(N2BASX,D2,DXCAO,1)
      END IF
C
C
      CALL MEMREL('TR1DEN',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('TR1DEN')
      RETURN
C
C *** end of subroutine TR1DEN
C
      END
C  /* Deck tr1fck */
      SUBROUTINE TR1FCK(NCSIM,NOSIM,FXC,FXV,FTV,
     *                  UBO,CMO,DV,DTV,WRK,LFRSAV)
C
C 890102 - hjaaj
C New interface routine for Tr1Fck
C 891103 - hjaaj : UBO(NORBT,NORBT) instead of UTBO(N2ORBT)
C
C MOTECC-90: The algorithms used in this module, TR1FCK, are
C            described in Chapter 8 Appendix 8B of MOTECC-90
C            "Calculation of Generalized Fock Matrices from
C            Regular and One-Index Transformed Integrals"
C
#include <implicit.h>
      DIMENSION FXC(N2ORBX,*),FXV(N2ORBX,*)
      DIMENSION FTV(N2ORBX,*),UBO(N2ORBX,*),CMO(*)
      DIMENSION DV(*),DTV(*), WRK(LFRSAV)
      PARAMETER (D2 = 2.0D0)
C
C Used from common blocks:
C   INFINP : DIRFCK
C   INFORB : NNORBX, N2ORBX, N2BASX, NISHT,NASHT, ...
C   INFVAR : JWOPSY
C   DFTCOM : DFT_SPINDNS
C
#include <maxorb.h>
#include <infinp.h>
#include <inforb.h>
#include <infvar.h>
#include <infpri.h>
#include <dftcom.h>
C
C
      CALL QENTER('TR1FCK')
C
C
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LFRSAV
      MCSIM  = MAX(0,NCSIM)
      MOSIM  = MAX(0,NOSIM)
      NDMAT  = 0
      IF (NISHT .GT. 0 .OR. ADDSRI) NDMAT = NDMAT + MOSIM
C     ... for FXCAO
      IF (NASHT .GT. 0) NDMAT = NDMAT + MOSIM
C     ... for FXVAO
      IF (DFT_SPINDNS)  NDMAT = NDMAT + MOSIM
C     ... for FXSAO
                        NDMAT = NDMAT + MCSIM
C     ... for FTVAO (note that NCSIM .eq. 0 if NASHT .eq. 0)
C     ... the IF test is necessary: NCSIM may be .lt. 0 to signal
C         only SOVECS wanted
      IF (DOMCSRDFT)    NDMAT = NDMAT + MCSIM
C     ... for MCSRDFT we also need "FTCAO" with sr-dft contribution
C         from transition density matrix DTV
      IF (DFT_SPINDNS)  NDMAT = NDMAT + MCSIM
C     ... for FTSAO
      LDFAO  = 2*NDMAT*N2BASX
      IF (DODFT .OR. DOHFSRDFT .OR. DOMCSRDFT) THEN
C     ... allocate for DTOTAO and DXTOTAO for DFTADD calls in SIRFCK
         LDAO  = 2*N2BASX
         LSRAO = 0
         IF (SRHYBR .AND. .NOT.DOHFSRDFT) THEN
            LDAO = LDAO + 2*N2BASX
            LSRAO = N2BASX
         END IF
C        ... if SRHYBR we also need space for DVAO and DXVAO in SIRFCK
C            and SRHYBR extra valence Fock type matrix
         LDFAO = LDFAO + LDAO + LSRAO
      ELSE
         LDAO = 0
         LSRAO = 0
      END IF
      CALL MEMGET('REAL',KDFAO,LDFAO,WRK,KFREE,LFREE)
C   Note: DXCAO, DXVAO, DTVAO must be contiguous.
C   Note: FXCAO, FXVAO, FTVAO must be contiguous.
      KDXCAO = KDFAO
      IF (NISHT .GT. 0 .OR. ADDSRI) THEN
         KDXVAO = KDXCAO + MOSIM*N2BASX
         IF (NISHT.EQ.0) CALL DZERO(WRK(KDXCAO),MOSIM*N2BASX)
C        ... because DXCAO not calculated in TR1DEN when NISHT.eq.0
C            but we need it for SRDFT to handle the DFT part
      ELSE
         KDXVAO = KDXCAO
      END IF
      IF (NASHT .GT. 0) THEN
         KDXSAO = KDXVAO + MOSIM*N2BASX
         IF (DFT_SPINDNS) THEN
            KDTCAO = KDXSAO + MOSIM*N2BASX
         ELSE
            KDTCAO = KDXSAO
         END IF
         IF (DOMCSRDFT) THEN
            KDTVAO = KDTCAO + MCSIM*N2BASX
            CALL DZERO(WRK(KDTCAO),MCSIM*N2BASX)
C           because DTCAO never calculated (we know always zero matrix)
C           but we need it for MCSRDFT to handle the DFT part
         ELSE
            KDTVAO = KDTCAO
         END IF
         KDTSAO = KDTVAO + MCSIM*N2BASX
         IF (DFT_SPINDNS) THEN
            KDAO   = KDTSAO + MCSIM*N2BASX
         ELSE
            KDAO   = KDTSAO
         END IF
         KFXCAO = KDAO + LDAO
      ELSE
         KDXSAO = KDXVAO
         KDTCAO = KDXVAO
         KDTVAO = KDTCAO
         KDTSAO = KDTCAO
         KDAO   = KDTVAO
         KFXCAO = KDAO + LDAO
      END IF
      IF (NISHT .GT. 0 .OR. ADDSRI) THEN
         KFXVAO = KFXCAO + MOSIM*N2BASX
      ELSE
         KFXVAO = KFXCAO
      END IF
      IF (NASHT .GT. 0) THEN
         KFXSAO = KFXVAO + MOSIM*N2BASX
         IF (DFT_SPINDNS) THEN
            KFTVAO = KFXSAO + MOSIM*N2BASX
         ELSE
            KFTVAO = KFXSAO
         END IF
         KFTCAO = KFTVAO + MCSIM*N2BASX
         IF (DOMCSRDFT) THEN
            KFTSAO = KFTCAO + MCSIM*N2BASX
         ELSE
            KFTSAO = KFTCAO
         END IF
         IF (DFT_SPINDNS) THEN
            KFSCR  = KFTSAO + MCSIM*N2BASX
         ELSE
            KFSCR  = KFTSAO
         END IF
      ELSE
         KFXSAO = KFXVAO
         KFTVAO = KFXVAO
         KFTCAO = KFTVAO
         KFTSAO = KFTVAO
         KFSCR  = KFTVAO
      END IF
      CALL MEMGET('INTE',KSYMDM,NDMAT,WRK,KFREE,LFREE)
      CALL MEMGET('INTE',KFCTYP,NDMAT,WRK,KFREE,LFREE)
C
C     If this is DFT (i.e. LDAO .gt. 0) then we must
C     get DTOTAO = DCAO + DVAO in WRK(KDAO) for SIRFCK
C     DFT call with DOATR .true. /Nov 2003
C
      IF (LDAO .GT. 0) THEN
         CALL FCKDEN(.TRUE.,(NASHT.GT.0),WRK(KDAO),
     *               WRK(KDAO+N2BASX),CMO,DV,WRK(KFREE),LFREE)
C        CALL FCKDEN(GETDC,GETDV,DCAO,DVAO,CMO,DV,WRK,LWRK)
         IF (NASHT .GT. 0) THEN
C           make DTOTAO in WRK(KDAO)
            CALL DAXPY(N2BASX,1.0D0,WRK(KDAO+N2BASX),1,WRK(KDAO),1)
         END IF
      END IF
C
      CALL TR1FCY(DIRFCK,NCSIM,NOSIM,FXC,FXV,FTV,
     *            UBO,CMO,DV,DTV,WRK(KSYMDM),WRK(KFCTYP),
     *            WRK(KDXCAO),WRK(KDXVAO),WRK(KDTVAO),
     *            WRK(KFXCAO),WRK(KFXVAO),WRK(KFTVAO),
     *            WRK,KFREE,LFREE)
C
      CALL MEMREL('TR1FCK',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('TR1FCK')
      RETURN
      END
C  /* Deck tr1fcy */
      SUBROUTINE TR1FCY(AODRCT,NCSIM,NOSIM,FXC,FXV,FTV,
     *                  UBO,CMO,DV,DTV,ISYMDM,IFCTYP,DXCAO,DXVAO,DTVAO,
     *                  FXCAO,FXVAO,FTVAO,WRK,KFREE,LFREE)
C
C Written 21-Nov-1983 Hans Jorgen Aa. Jensen.
C Revisions:
C  3-May-1984 -- hjaaj
C  16-Oct-1984 ha/hjaaj (I*2 + FPS for reading LUSUPM)
C  29-Dec-1984 hjaaj (NCONF=0 option for call by ORBLIN)
C  15-Jan-1985 hjaaj (moved calc. of UTBO outside TR1FCY)
C  30 Jan '85 hjaaj  (ver 2; separated conf. and orb. trial vectors)
C  31-Jan-1985 hjaaj (NCONF=0 option changed to NCSIM.lt.0)
C 891103 - hjaaj : UBO(NORBT,NORBT) instead of UTBO(N2ORBT)
C 891106 - hjaaj : only FXCAO if NISHT.gt.0
C 921201 - hjaaj : Fock matrices from AO integrals
C                  (partly following Hans Agren 911202)
C DFT modification tuh, fixed by pawsa
C
C Purpose:
C
C  Construction of NOSIM one-index transformed, inactive and active
C  Fock matrices (FXC and FXV) and construction of NCSIM transition
C  active Fock matrices FTV.  If nasht.eq.0 then FXV is a zero matrix.
C
C Algorithm:
C
C  Construction of inactive and active parts of the
C  Fock matrix in AO-basis (FXCAO, FXVAO, and FTVAO) using a
C  list of atomic super matrix elements. The matrices
C  are transformed to MO-basis (FXC, FXV, and FTV) in AUTPV.
C  This completes FTV, to complete the one-index transformation
C  we add the one-index transformed FC and FV matrices to
C  FXC and FXV, resp. (done outside).
C
C
C Scratch:
C  WRK();
C
C       Note: FXCAO, FXVAO, FTCAO, FTVAO are assumed to be contiguous
C             DXCAO, DXVAO, DTCAO, DTVAO are assumed to be contiguous
C       ( DTCAO is zero matrix, but space is needed for MC SR-DFT
C         effective one-electron matrix for transition DTV )
C
C
C
#include <implicit.h>
      DIMENSION FXC(N2ORBX,*),FXV(N2ORBX,*)
      DIMENSION FTV(N2ORBX,*),UBO(N2ORBX,*),CMO(*)
      DIMENSION DV(*),DTV(*),ISYMDM(*),IFCTYP(*)
      DIMENSION DXCAO(N2BASX,*),DXVAO(N2BASX,*),DTVAO(N2BASX,*)
      DIMENSION FXCAO(N2BASX,*),FXVAO(N2BASX,*),FTVAO(N2BASX,*)
      DIMENSION WRK(*)
      LOGICAL   AODRCT, P6SAVE, DFTADX
C
#include <dummy.h>
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0, D1=1.0D0, D2 = 2.0D0)
C
C Used from common blocks:
C   INFORB : N2ORBX,N2ORBT,N2BASX,NNASHX
C   INFIND : ?
C   INFVAR : JWOPSY
C   INFLIN : IPRLIN
C   INFDIM : NBASMA,?
C   INFTAP : LUSUPM
C   DFTCOM : DFT_SPINDNS, ?
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infind.h>
#include <infinp.h>
#include <infvar.h>
#include <inflin.h>
#include <infdim.h>
#include <inftap.h>
#include <infpri.h>
#include <dftcom.h>
#include <abainf.h>
Cholesky
#include <ccdeco.h>
Cholesky
C
#if defined (VAR_DEBUG)
#include <idbg.h>
#endif
#ifndef FOR_THE_FUTURE
c  this is a hack, to transfer KSYMOP to DEQ27 inside DFTATR
#include "wrkrsp.h"
#endif
C
      CALL QENTER('TR1FCY')
      KFRSAV = KFREE
C
C *** Prepare for reading LUSUPM by constructing ao density
C     matrices one B vector at a time.
C
      DO 750 IOSIM = 1,NOSIM
C
C ***** Construct modified density matrices needed to calculate
C       the Fock matrices with one-index transformed integrals.
C
        CALL TR1DEN (CMO,UBO(1,IOSIM),DV,DXCAO(1,IOSIM),DXVAO(1,IOSIM),
     &               WRK(KFREE),LFREE)
        IF (HSROHF) THEN
           CALL DAXPY(N2BASX,D1,DXVAO(1,NOSIM),1,DXCAO(1,NOSIM),1)
           CALL DSCAL(N2BASX,-D1,DXVAO(1,NOSIM),1)
        END IF
        IF (P6FLAG(24) .OR. IPRLIN .GE. 24) THEN
          IF (NISHT.GT.0 .OR. ADDSRI) THEN
            WRITE (LUPRI,2010) IOSIM,NOSIM
CTODO       CALL PRINRM(DXCAO(1,IOSIM),N2BASX,LUPRI)
            CALL OUTPUT(DXCAO(1,IOSIM),1,NBAST,1,NBAST,NBAST,NBAST,
     &         -1,LUPRI)
          END IF
          IF (NASHT.GT.0) THEN
            WRITE (LUPRI,2020) IOSIM,NOSIM
CTODO       CALL PRINRM(DXVAO(1,IOSIM),N2BASX,LUPRI)
            CALL OUTPUT(DXVAO(1,IOSIM),1,NBAST,1,NBAST,NBAST,NBAST,
     &         -1,LUPRI)
          END IF
        END IF
  750 CONTINUE
C
 2010 FORMAT(/,' DXCAO = one-index transformed inactive one-el.',
     *  ' d.m., AO basis',/5X,'(no.',I3,' of',I3,')')
 2020 FORMAT(/,' DXVAO = one-index transformed active one-el. d.m.'
     *  ,', AO basis',/5X,'(no.',I3,' of',I3,')')
C
C ***** Construct the transition AO density matrices for active
C       orbitals to DTVAO
C       (P6FLAG(18) prints inside FCKDEN, thus disabled here)
C
      IF (NCSIM .GT. 0) THEN
         P6SAVE = P6FLAG(18)
         P6FLAG(18) = .FALSE.
         DO 800 ICSIM = 1,NCSIM
            JDTV = 1 + (ICSIM-1)*NNASHX
            CALL FCKDEN(.FALSE.,.TRUE.,DUMMY,DTVAO(1,ICSIM),
     *                  CMO,DTV(JDTV),WRK(KFREE),LFREE)
C           CALL FCKDEN(GETDC,GETDV,DCAO,DVAO,CMO,DV,WRK,LWRK)
C
            IF (P6FLAG(24) .OR. IPRLIN.GE.24) THEN
               WRITE (LUPRI,2030) ICSIM,NCSIM
               CALL OUTPUT(DTVAO(1,ICSIM),1,NBAST,1,NBAST,
     &                     NBAST,NBAST,-1,LUPRI)
            END IF
C
  800    CONTINUE
         P6FLAG(18) = P6SAVE
      END IF
C
C
 2030 FORMAT(/,' DTVAO = the transition active one-el. d.m., AO basis',
     *  5X,'(no.',I3,' of',I3,')')
C
C *************************************
C
C       Note: FXCAO, FXVAO, FTCAO, FTVAO are assumed to be contiguous
C             DXCAO, DXVAO, DTCAO, DTVAO are assumed to be contiguous
C
C ***** count and initialize atomic Fock matrices
C
      MOSIM = MAX(0,NOSIM)
      MCSIM = MAX(0,NCSIM)
      NDMAT = 0
      IF (NISHT .GT. 0 .OR. ADDSRI) NDMAT = NDMAT + MOSIM
C     ... for FXCAO
      IF (NASHT .GT. 0) NDMAT = NDMAT + MOSIM
C     ... for FXVAO
      IF (DFT_SPINDNS ) NDMAT = NDMAT + MOSIM
C     ... for FXSAO
                        NDMAT = NDMAT + MCSIM
C     ... for FTVAO (note that MCSIM .eq. 0 if NASHT .eq. 0)
C     ... the MCSIM is necessary: NCSIM may be .lt. 0 to signal
C         only SOVECS wanted
      IF (DOMCSRDFT   ) NDMAT = NDMAT + MCSIM
C     ... for MCSRDFT we also need "FTCAO" with sr-dft contribution
C         from transition density matrix DTV
      IF (DFT_SPINDNS ) NDMAT = NDMAT + MCSIM
C     ... for FTSAO
C
C
C ***** TWO-ELECTRON PART (all one-electron contributions
C ***** are accounted for through the one-index transformation
C ***** of the inactive Fock matrix above).
C
      IF (NDMAT .GT. 0) THEN
         CALL DZERO(FXCAO,NDMAT*N2BASX)
      IF (LUSUPM .NE. -1 .AND.
     &   .NOT. (HSROHF .OR. DFTRUN .OR. DOCISRDFT .OR. DOHFSRDFT .OR.
     &    DOMCSRDFT)) THEN
         CALL RDSUPM(NDMAT,FXCAO,DXCAO,WRK(KFREE),LFREE)
C        CALL RDSUPM(NSIM,FMAT,DMAT,WORK,LFREE)
      ELSE
         DO 400 IDMAT = 1,NDMAT
            ISYMDM(IDMAT) = JWOPSY
            IFCTYP(IDMAT) = 13
            IF (CHOINT) IFCTYP(IDMAT) = 3
  400    CONTINUE
         IF (HSROHF) THEN
            DO IDMAT = NOSIM+1,NDMAT
               IFCTYP(IDMAT) = 12
               IF (CHOINT) IFCTYP(IDMAT) = 2
            END DO
         END IF
         DFTADX = DFTADD
         DFTADD = .FALSE.
C        ... Note : (DODFT.AND.NOT.DFTADD) now tricks 2'nd order 
C          F        DFT contributions in SIRFCK /mc-srdft work
         CALL SIRFCK(FXCAO,DXCAO,NDMAT,ISYMDM,IFCTYP,
     *               AODRCT,WRK(KFREE),LFREE)
         DFTADD = DFTADX
         IF (HSROHF) THEN
            CALL DAXPY(NOSIM*N2BASX,-D1,FXVAO,1,FXCAO,1)
         END IF
C        ... if (AODRCT) construct Fock matrix directly
C            else construct Fock matrix from LUINTA
C            ISYMDM : symmetry of density matrix (DUNP)
C            IFCTYP = 13 : singlet Fock matrix and DUNP symmetric
         DFTADD = DFTADX
      END IF
      END IF
#ifdef FOR_THE_FUTURE
C
C add DFT contribution in AO BASIS
C
      call quit(' sss ')
      IF (DFTADD) THEN
         IF (NASHT.GT.0) CALL QUIT('TR1FCY: open shell DFT not yet')
c           TODO         DFTIMG = .FALSE.
C           DXCAO contains effective
            CALL DFTATR_new_in_ao(NDMAT,FXCAO,DXCAO,CMO,.FALSE.,
     &                    JWOPSY,WRK(KFREE),LFREE)
c-rsp       CALL DFTATR(FCONE(1,1,IOSIM),CMO,ZYMAT(1,IOSIM),TRPLET,
c-rsp&                    KSYMOP,WRK(KWRK0),LFREE)
      END IF      
#endif
C
C
      IF (P6FLAG(24) .OR. IPRLIN.GE.24) THEN
        IF (DOMCSRDFT) THEN
           ICOFF = NCSIM
        ELSE
           ICOFF = 0
        END IF
        DO 210 ICSIM = 1,NCSIM
          IF (DOMCSRDFT) THEN
             WRITE (LUPRI,3031) ICSIM,NCSIM
             CALL OUTPUT(FTVAO(1,ICSIM),1,NBAST,1,NBAST,
     &                   NBAST,NBAST,-1,LUPRI)
          END IF
          WRITE (LUPRI,3030) ICSIM,NCSIM
          CALL OUTPUT(FTVAO(1,ICOFF+ICSIM),1,NBAST,1,NBAST,
     &                NBAST,NBAST,-1,LUPRI)
  210   CONTINUE
        DO 215 IOSIM = 1,NOSIM
          IF (NISHT.GT.0 .OR. ADDSRI) THEN
            WRITE (LUPRI,3010) IOSIM,NOSIM
            CALL OUTPUT(FXCAO(1,IOSIM),1,NBAST,1,NBAST,
     &                  NBAST,NBAST,-1,LUPRI)
          END IF
          IF (NASHT.GT.0) THEN
            WRITE (LUPRI,3020) IOSIM,NOSIM
            CALL OUTPUT(FXVAO(1,IOSIM),1,NBAST,1,NBAST,
     &                  NBAST,NBAST,-1,LUPRI)
          END IF
  215   CONTINUE
      END IF
C
 3010 FORMAT (/' FXCAO matrix (no.',I3,' of',I3,')')
 3020 FORMAT (/' FXVAO matrix (no.',I3,' of',I3,')')
 3030 FORMAT (/' FTVAO matrix (no.',I3,' of',I3,')')
 3031 FORMAT (/' FTCAO matrix for MC SR-DFT (no.',I3,' of',I3,')')
C
C
C *** Calculate the active Fock matrix FTV based on the transition
C     density matrix DTV :
C     Transform FTVAO to MO basis to form FTV
C
      IF (DOMCSRDFT) THEN
C        ... we have both "FTCAO" and FTVAO
         MCSIM = 2*NCSIM
      ELSE
C        ... we have only FTVAO
         MCSIM = NCSIM
      END IF
      DO 4900 ICSIM = 1,MCSIM
         IF (DOMCSRDFT) THEN
C           .. we swap order from "FTCAO", FTVAO
C              to FTV, "FTC" in order that general sirius
C              routines can work as before
            IF (ICSIM .GT. NCSIM) THEN
               JCSIM = ICSIM - NCSIM
            ELSE
               JCSIM = ICSIM + NCSIM
            END IF
         ELSE
            JCSIM = ICSIM
         END IF
         DO 4300 ISYM=1,NSYM
            JSYM = MULD2H(ISYM,JWOPSY)
            NORBI=NORB(ISYM)
            NORBJ=NORB(JSYM)
         IF (NORBI.EQ.0 .OR. NORBJ.EQ.0) GO TO 4300
C
            CALL AUTPV(ISYM,JSYM,CMO(ICMO(ISYM)+1),CMO(ICMO(JSYM)+1),
     &                 FTVAO(1,ICSIM),NBAS,NBAST,
     &                 FTV(1,JCSIM),NORB,NORBT, WRK(KFREE),LFREE)
C           CALL AUTPV(ISYM,JSYM,U,V,PRPAO,NBAS,NBAST,PRPMO,NORB,NORBT,
C    &                 WRK,LWRK)
C
 4300    CONTINUE
 4900 CONTINUE
C
C
C
C ***** Transform FXCAO and FXVAO to MO basis and add to FXC and
C ***** FXV, resp. (the one-index transformations of FC and FV are
C ***** added after this routine).
C           CALL AUTPV(ISYM,JSYM,U,V,
C    &                 PRPAO,NBAS,NBAST,PRPMO,NORB,NORBT,WRK,LWRK)
C
C
      DO 5900 IOSIM = 1,NOSIM
         DO 300 ISYM=1,NSYM
            JSYM   = MULD2H(ISYM,JWOPSY)
            NORBI  = NORB(ISYM)
            NORBJ  = NORB(JSYM)
        IF (NORBI.EQ.0 .OR. NORBJ.EQ.0) GO TO 300
            IF (NISHT .GT. 0 .OR. ADDSRI) THEN
               CALL AUTPV(ISYM,JSYM,CMO(ICMO(ISYM)+1),CMO(ICMO(JSYM)+1),
     &                    FXCAO(1,IOSIM),NBAS,NBAST,
     &                    FXC(1,IOSIM),NORB,NORBT,
     &                    WRK(KFREE),LFREE)
            END IF
            IF (NASHT .GT. 0) THEN
               CALL AUTPV(ISYM,JSYM,CMO(ICMO(ISYM)+1),CMO(ICMO(JSYM)+1),
     &                    FXVAO(1,IOSIM),NBAS,NBAST,
     &                    FXV(1,IOSIM),NORB,NORBT,
     &                    WRK(KFREE),LFREE)
            END IF
  300    CONTINUE
C
#if 0 && !defined(FOR_THE_FUTURE)
C
C     add DFT contribution in MO BASIS, this has been deprecated by
c        the call in sirlintrn.F
C
      IF (DODFT .AND. .NOT. MOLHES) THEN
         IF (NASHT.GT.0) CALL QUIT('TR1FCY: open shell DFT not yet')
c        TODO         DFTIMG = .FALSE.
         ksymop = jwopsy
         call dft_lin_respf(1,fxc(1,1),cmo,ubo(1,1),.FALSE.,
     &                      jwopsy,wrk(kfree),lfree)
c         CALL DFTATR(FXC(1,1),CMO,UBO(1,1),.false.,
c                     JWOPSY,WRK(KFREE),LFREE)
      END IF      
#endif
       IF (SRHYBR .AND. MCTYPE.NE.0) THEN
C hjaaj TODO
          WRITE(LUPRI,*) 'WARNING> SRHYBR_XV matrix missing!'
       END IF
#ifdef FOR_THE_PAST
C ... we have moved this to SIRFCK 26.nov.2003 jkp+hjaaj
C
C     add DFT contribution in MO BASIS
C
      IF (DFTADD) THEN
         DFTIMG = .FALSE.
         ksymop = jwopsy
         CALL DFTATR(FXC(1,IOSIM),CMO,DV,UBO(1,IOSIM),.false.,
     &               JWOPSY,WRK(KFREE),LFREE)
c-rsp    CALL DFTATR(FCONE(1,1,IOSIM),CMO,ZYMAT(1,IOSIM),TRPLET,
c-rsp&               KSYMOP,WRK(KWRK0),LFREE)
      END IF
#endif
 5900 CONTINUE
C
C *** End of subroutine TR1FCY
C
      CALL QEXIT('TR1FCY')
      RETURN
      END
C  /* Deck tr1fd */
      SUBROUTINE TR1FD(NOSIM,FXC,FXV, UBO,FC,FV, WRK,LFRSAV)
C
C Written Feb 1990 Hans Jorgen Aa. Jensen.
C
C Purpose:
C
C  To complete the construction of NOSIM one-index transformed,
C  inactive and active Fock matrices (FXC and FXV) by
C  adding the one-index transformed FC and FV matrices to
C  FXC and FXV, resp.
C
#include <implicit.h>
      DIMENSION FXC(N2ORBX,*),FXV(N2ORBX,*)
      DIMENSION UBO(N2ORBX,*), FC(*),FV(*)
      DIMENSION WRK(LFRSAV)
C
      PARAMETER ( DP5 = 0.5D0 )
C
C Used from common blocks:
C   INFORB : N2ORBX,N2ORBT,N2BASX,NNASHX
C   INFLIN : IPRLIN
C   INFPRI : P6FLAG(*)
C
#include <priunit.h>
#include <inforb.h>
#include <inflin.h>
#include <infpri.h>
C
      CALL QENTER('TR1FD ')
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LFRSAV
C
C
C ***** Calculate one-index transformation of FC and FV, the inactive
C ***** and active Fock matrices, resp., and add to FXC and FXV.
C ***** This completes FXC and FXV.
C
C
      CALL MEMGET('REAL',KFTMP,NNORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KFC  ,N2ORBX,WRK,KFREE,LFREE)
      CALL PKSYM1(WRK(KFTMP),FC,NORB,NSYM,-1)
      CALL DSPTSI(NORBT,WRK(KFTMP),WRK(KFC))
      IF (NASHT .GT. 0) THEN
         CALL MEMGET('REAL',KFV  ,N2ORBX,WRK,KFREE,LFREE)
         CALL PKSYM1(WRK(KFTMP),FV,NORB,NSYM,-1)
         CALL DSPTSI(NORBT,WRK(KFTMP),WRK(KFV))
      ELSE
         KFV = KFC
C        ... to avoid compiler messages
      END IF
      DO 5900 IOSIM = 1,NOSIM
         CALL DSCAL(N2ORBX,DP5,FXC(1,IOSIM),1)
         CALL TR1UH1(UBO(1,IOSIM),WRK(KFC),FXC(1,IOSIM),1)
         IF (NASHT .GT. 0) THEN
            CALL DSCAL(N2ORBX,DP5,FXV(1,IOSIM),1)
            CALL TR1UH1(UBO(1,IOSIM),WRK(KFV),FXV(1,IOSIM),1)
         END IF
C
         IF (P6FLAG(21) .OR. IPRLIN.GE.21) THEN
            WRITE (LUPRI,5020) IOSIM,NOSIM
            CALL OUTPUT(FXC(1,IOSIM),1,NORBT,1,NORBT,
     &                  NORBT,NORBT,1,LUPRI)
            IF (NASHT .GT. 0) THEN
               WRITE (LUPRI,5030) IOSIM,NOSIM
               CALL OUTPUT(FXV(1,IOSIM),1,NORBT,1,NORBT,
     &                     NORBT,NORBT,1,LUPRI)
            END IF
         END IF
C
 5020 FORMAT (/' FXC matrix (no.',I3,' of',I3,')')
 5030 FORMAT (/' FXV matrix (no.',I3,' of',I3,')')
C
 5900 CONTINUE
      CALL MEMREL('TR1FD.FXC/V',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
C
C *** End of subroutine TR1FD
C
      CALL QEXIT('TR1FD ')
      RETURN
      END
C  /* Deck tr1h1 */
      SUBROUTINE TR1H1(NOSIM,H1,H1X,UKAP,WRK)
C
C Written 31-Jan-1985 by Hans Joergen Aa. Jensen
C Revised 891006-hjaaj
C
C Purpose:
C  Add the one-index transformation of H1 to H1X,
C  transforming with second index of UKAP.
C
C  H1X(i,j) = H1X(i,j) +
C
C             sum(p) ( UKAP(i,p) * H1(p,j) + UKAP(j,p) * H1(i,p) )
C
C  H1 is symmetric and lower triangle packed, no assumptions
C  on UKAP.
C
#include <implicit.h>
      DIMENSION H1(*), H1X(NNORBT,*), UKAP(N2ORBX,*), WRK(*)
      PARAMETER ( D2 = 2.0D0 )
C
C Used from common blocks:
C   INFORB : NSYM,NISH,NORB,NBAS,IIORB,I2ORB,NNORBT,N2ORBT,...
C
C
#include <inforb.h>
C
      CALL QENTER('TR1H1 ')
C
C     Core allocation (using N2ORBT .ge. NORBMA*NORBMA)
C
      KH1   = 1
      KH1XA = KH1   + N2ORBT
      KH1XB = KH1XA + N2ORBT
C
      DO 900 ISYM=1,NSYM
         NORBI=NORB(ISYM)
      IF (NORBI.LE.0) GO TO 900
         JH1  = IIORB(ISYM) + 1
         JUTK = IORB(ISYM)*(NORBT+1) + 1
C
C --- unpack H1
C
         CALL DSPTSI(NORBI,H1(JH1),WRK(KH1))
C
         DO 800 IOSIM = 1,NOSIM
C           H1XA(i,j) = sum(p) UKAP(i,p) * H1(p,j)
            CALL DGEMM('N','N',NORBI,NORBI,NORBI,1.D0,
     &                 UKAP(JUTK,IOSIM),NORBT,
     &                 WRK(KH1),NORBI,0.D0,
     &                 WRK(KH1XA),NORBI)
C           H1XB(ij) = 0.5 * (H1XA(i,j) + H1XA(j,i))
            CALL DGETSP(NORBI,WRK(KH1XA),WRK(KH1XB))
C           H1X(ij) = H1X(ij) + 2 * H1XB(ij)
            CALL DAXPY(NNORB(ISYM),D2,WRK(KH1XB),1,H1X(JH1,IOSIM),1)
C
  800    CONTINUE
C 800    END DO IOSIM
  900 CONTINUE
C 900 END DO ISYM
C
C *** end of subroutine TR1H1
C
      CALL QEXIT('TR1H1 ')
      RETURN
      END
C  /* Deck tr1uh1 */
      SUBROUTINE TR1UH1(UBOVEC,H1,H1X,IH1SYM)
C
C 891012-hjaaj
C
C One-index transform unpacked H1(norbt,norbt) of symmetry IH1SYM
C to H1X(norbt,norbt) using UBOVEC(norbt,norbt) of symmetry
C JWOPSY.
C
C Input : H1(a,b) of symmetry IH1SYM
C Output: H1X(a,b) = H1X(a,b)
C                  + H1(a,b) one-index transformed with UBOVEC
C                  = H1X(a,b)
C                  + sum(c) [ UBOVEC(a,c) H1(c,b)
C                           + UBOVEC(b,c) H1(a,c) ]
C
#include <implicit.h>
      DIMENSION UBOVEC(NORBT,NORBT), H1(NORBT,NORBT), H1X(NORBT,NORBT)
C
C Used from common blocks:
C  INFORB : NSYM, MULD2H, NORBT, ...
C  INFVAR : JWOPSY
C
#include <maxorb.h>
#include <inforb.h>
#include <infvar.h>
C
C
C ONE INDEX TRANSFORM FIRST INTEGRAL INDEX
C
C  (P~B) =  SUM(A) UBOVEC(P,A)*(AB)
C
      DO 1100 IASYM = 1,NSYM
         IBSYM = MULD2H(IASYM,IH1SYM)
         IPSYM = MULD2H(IASYM,JWOPSY)
         IQSYM = MULD2H(IBSYM,JWOPSY)
         NORBA = NORB(IASYM)
         NORBB = NORB(IBSYM)
         NORBP = NORB(IPSYM)
         IF ((NORBP.NE.0) .AND. (NORBA.NE.0) .AND. (NORBB.NE.0)) THEN
            IAST  = IORB(IASYM) + 1
            IBST  = IORB(IBSYM) + 1
            IPST  = IORB(IPSYM) + 1
            CALL DGEMM('N','N',NORBP,NORBB,NORBA,1.D0,
     &                 UBOVEC(IPST,IAST),NORBT,
     &                 H1(IAST,IBST),NORBT,1.D0,
     &                 H1X(IPST,IBST),NORBT)
         END IF
 1100 CONTINUE
C
C  ADD ONE-INDEX TRANSFORM FROM SECOND INTEGRAL INDEX
C
C   H1X(P,B) = H1X(P,B) + H1X(B,P)  P.ge.B
C
      IH1XSY = MULD2H(IH1SYM,JWOPSY)
      DO 1200 IBSYM = 1,NSYM
         IPSYM = MULD2H(IBSYM,IH1XSY)
      IF (IPSYM .LT. IBSYM) GO TO 1200
         NORBB = NORB(IBSYM)
         NORBP = NORB(IPSYM)
      IF (NORBB.EQ.0 .OR. NORBP.EQ.0) GO TO 1200
         IORBB = IORB(IBSYM)
         IORBP = IORB(IPSYM)
         DO 1250 IP = IORBP+1,IORBP+NORBP
            IBEND = MIN(IP,IORBB+NORBB)
            DO 1260 IB = IORBB+1,IBEND
               H1X(IP,IB) = H1X(IP,IB) + H1X(IB,IP)
               H1X(IB,IP) = H1X(IP,IB)
 1260       CONTINUE
 1250    CONTINUE
 1200 CONTINUE
C
C     End of TR1UH1
      RETURN
      END
C  /* Deck adfxqm */
      SUBROUTINE ADFXQM(IC,ID,FXQ,H2XCD,PV,WRK,LWRK)
C
C Jan 90 hjaaj
C
C Add (p~ v~ | x y) contributions to FXQ from Mulliken integrals.
C
#include <implicit.h>
      DIMENSION FXQ(NORBT,NASHT), H2XCD(NORBT,NORBT),
     &          PV(NNASHX,NNASHX),WRK(LWRK)
C
C Used from common blocks:
C  INFORB : NORBT, NASHT, N2ASHX, NNASHX, MULD2H(8,8)
C  INFIND : IROW(:),ISMO(:),NSM(:),ICH(:),IOBTYP(:)
C  INFVAR : JWOPSY
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infind.h>
#include <infvar.h>
C
      KH2AC1 = 1
      KH2ACF = KH2AC1 + N2ASHX
      KNEED  = KH2ACF + NNASHX
      IF (KNEED .GT. LWRK) CALL ERRWRK('ADFXQM',KNEED,LWRK)
C
      CALL GETAC1(H2XCD,WRK(KH2AC1))
C     CALL GETAC1(HH, HHAC)
      CALL DGEFSP(NASHT,WRK(KH2AC1),WRK(KH2ACF))
C     CALL DGEFSP(N,AGE,ASP)
C
C We now have the folded one-index transformed integrals in H2ACF:
C   H2ACF(uv) = (2-DELTAuv) * (UX VX / C D)
C where
C   (AX BX / C D) = (AX B / C D) + (A BX / C D)
C   (A  BX / C D) = SUM(r)  BORB(b,r) * H2CD(a,r)
C
C
C Add PV(t,d;a,b)*(c d;ax bx) contributions to FXQ(c,t)
C
      IC1 = IC
      ID1 = ID
      IROUND = 0
  420    IROUND = IROUND + 1
         ITYPD  = IOBTYP(ID1)
      IF (ITYPD .EQ. JTACT) THEN
         ICSYM  = ISMO(IC1)
         NDW    = ICH(ID1)
         DO 450 NTW = 1,NASHT
            ITSYM = NSM(NTW)
         IF (MULD2H(ITSYM,ICSYM) .NE. JWOPSY) GO TO 450
            IF (NDW .GE. NTW) THEN
               NDTW = IROW(NDW) + NTW
            ELSE
               NDTW = IROW(NTW) + NDW
            END IF
            FXQ(IC1,NTW) = FXQ(IC1,NTW) +
     *         DDOT(NNASHX,WRK(KH2ACF),1,PV(1,NDTW),1)
  450    CONTINUE
      END IF
C
C        If (cd) active-active and c .ne. d we need to
C        consider (dc) distribution also.
C
      IF ( IROUND .EQ. 1 .AND. IC .NE. ID ) THEN
         IC1 = ID
         ID1 = IC
         GO TO 420
      END IF
      RETURN
      END
C  /* Deck autpv */
      SUBROUTINE AUTPV(ISYM,JSYM,U,V,PRPAO,NBAS,NBAST,PRPMO,NORB,NORBT,
     &                 WRK,LWRK)
C
C Jan. 90, Hans Joergen Aa. Jensen
C
C TRANSFORM (ISYM,JSYM) SYMMETRY BLOCK OF
C THE MATRIX PRPAO FROM AO BASIS TO MO BASIS
C
C Input
C  PRPAO ONE-ELECTRON MATRIX OVER SYMMETRY AO'S
C  U     MO COEFFICIENTS FOR SYMMETRY ISYM
C  V     MO COEFFICIENTS FOR SYMMETRY JSYM
C
C Output
C  PRPMO (ISYM,JSYM) BLOCK OF THE ONE-ELECTRON
C        MATRIX TRANSFORMED TO MO BASIS. PRPMO HAS
C        DIMENSION NORBT*NORBT
C
C Scratch
C  WRK  dimension: NBAS(ISYM)
C
#include <implicit.h>
      DIMENSION U(*), V(*), NBAS(8), NORB(8)
      DIMENSION PRPAO(NBAST,*),PRPMO(NORBT,*),WRK(*)
C
      PARAMETER (D0 = 0.0D0)
C
      IF (LWRK .LT. NBAST) CALL ERRWRK('AUTPV',NBAST,LWRK)
      IORBI = 0
      IBASI = 0
      DO 10 I = 1,ISYM-1
         IORBI = IORBI + NORB(I)
         IBASI = IBASI + NBAS(I)
   10 CONTINUE
      NORBI = NORB(ISYM)
      NBASI = NBAS(ISYM)
      IORBJ = 0
      IBASJ = 0
      DO 20 J = 1,JSYM-1
         IORBJ = IORBJ + NORB(J)
         IBASJ = IBASJ + NBAS(J)
   20 CONTINUE
      NORBJ = NORB(JSYM)
      NBASJ = NBAS(JSYM)
C
C  TRANSFORM TO MO BASIS
C
      DO 400 I=1,NORBI
         JU = 1 + (I-1)*NBASI
         DO 200 JB=1,NBASJ
            WRK(JB)=DDOT(NBASI,U(JU),1,PRPAO(IBASI+1,IBASJ+JB),1)
 200     CONTINUE
         DO 300 J=1,NORBJ
            JV = 1 + (J-1)*NBASJ
            PRPMO(IORBI+I,IORBJ+J)=PRPMO(IORBI+I,IORBJ+J)
     &         + DDOT(NBASJ,WRK,1,V(JV),1)
 300     CONTINUE
 400  CONTINUE
      RETURN
      END
C  /* Deck fxdfin */
      SUBROUTINE FXDFIN(NSIM,FCOCO,FVOCO,FCOEX,FVOEX)
C 1)
C SCALE FCOCO,FVOCO AND FVOEX TO HAVE RIGHT FACTORS
C 2)
C ADD CONTRIBUTIONS FROM OCCUPIED SECONDARY TO FCOCO AND FVOCO
C (USE THAT COULOMB CONTRIBUTION IS SYMMETRIC) AND
C 3)
C SYMMETRISE FCOEX AND FVOEX
C 4)
C SUM UP ONE INDEX TRANSFORMED DENSITY CONTRIBUTIONS TO FOCK
C MATRICES ( FCOEX  = FCOCO  + FCOEX )
C
#include <implicit.h>
C
      DIMENSION FCOCO(NORBT,NORBT,*),FVOCO(NORBT,NORBT,*)
      DIMENSION FCOEX(NORBT,NORBT,*),FVOEX(NORBT,NORBT,*)
C
      PARAMETER ( D4 = 4.0D0 , DM2 = -2.0D0 , DP5 = 0.5D0 )
C
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <inflin.h>
#include <inforb.h>
#include <infind.h>
#include <infvar.h>
C
      IF(NISHT.GT.0) CALL DSCAL(NORBT*NORBT*NSIM,D4,FCOCO,1)
      IF(NASHT.GT.0) CALL DSCAL(NORBT*NORBT*NSIM,DM2,FVOCO,1)
      DO 100 IPSYM = 1,NSYM
         IORBP = IORB(IPSYM) + 1
         NOCCP = NOCC(IPSYM)
         IQSYM = MULD2H(IPSYM,JWOPSY)
         IORBQ = IORB(IQSYM) + 1
         NOCCQ = NOCC(IQSYM)
         NORBQ = NORB(IQSYM)
         IF ((NOCCP.EQ.0).OR.((NORBQ-NOCCQ).EQ.0)) GO TO 100
         DO 400 ISIM = 1,NSIM
            DO 200 IP = IORBP,IORBP+NOCCP-1
               DO 300 IQ = IORBQ+NOCCQ,IORBQ+NORBQ-1
                  FCOCO(IP,IQ,ISIM) = FCOCO(IQ,IP,ISIM)
                  IF (NASHT.GT.0) FVOCO(IP,IQ,ISIM) = FVOCO(IQ,IP,ISIM)
 300           CONTINUE
 200        CONTINUE
 400     CONTINUE
 100  CONTINUE
      IF ( IPRLIN.GT.40 ) THEN
         DO 3010 ISIM = 1, NSIM
            WRITE(LUPRI,*)
     *      ISIM,' FCOCO matrix of',NSIM,'after occ-uocc block added'
            CALL OUTPUT(FCOCO(1,1,ISIM),1,NORBT,1,NORBT,
     *                  NORBT,NORBT,1,LUPRI)
            WRITE(LUPRI,*)
     *      ISIM,' FVOCO matrix of',NSIM,'after occ-uocc block added'
            CALL OUTPUT(FVOCO(1,1,ISIM),1,NORBT,1,NORBT,
     *                  NORBT,NORBT,1,LUPRI)
 3010    CONTINUE
      END IF
      DO 500 IPSYM = 1,NSYM
         IORBP = IORB(IPSYM) + 1
         NORBP = NORB(IPSYM)
         IQSYM = MULD2H(IPSYM,JWOPSY)
         IORBQ = IORB(IQSYM) + 1
         NORBQ = NORB(IQSYM)
         IF ((NORBQ.EQ.0).OR.(NORBP.EQ.0)) GO TO 500
         DO 800 ISIM = 1,NSIM
            DO 600 IP = IORBP,IORBP+NORBP-1
               DO 700 IQ = IORBQ,IORBQ+NORBQ-1
                  IF (IP.GE.IQ) THEN
                     FCOEX(IP,IQ,ISIM) = FCOEX(IP,IQ,ISIM)
     *                                 + FCOEX(IQ,IP,ISIM)
                     FCOEX(IQ,IP,ISIM) = FCOEX(IP,IQ,ISIM)
                     IF (NASHT.GT.0) THEN
                        FVOEX(IP,IQ,ISIM) = ( FVOEX(IP,IQ,ISIM)
     *                                    + FVOEX(IQ,IP,ISIM) ) * DP5
                        FVOEX(IQ,IP,ISIM) = FVOEX(IP,IQ,ISIM)
                     END IF
                  END IF
 700           CONTINUE
 600        CONTINUE
 800     CONTINUE
 500  CONTINUE
      DO 900 IPSYM = 1,NSYM
         IORBP = IORB(IPSYM) + 1
         NORBP = NORB(IPSYM)
         IQSYM = MULD2H(IPSYM,JWOPSY)
         IORBQ = IORB(IQSYM) + 1
         NORBQ = NORB(IQSYM)
         IF ((NORBQ.EQ.0).OR.(NORBP.EQ.0)) GO TO 900
         DO 910 ISIM = 1,NSIM
            DO 920 IP = IORBP,IORBP+NORBP-1
               DO 930 IQ = IORBQ,IORBQ+NORBQ-1
                  FCOEX(IP,IQ,ISIM) = FCOCO(IP,IQ,ISIM)
     *                              + FCOEX(IP,IQ,ISIM)
                  IF (NASHT.GT.0) FVOEX(IP,IQ,ISIM) = FVOCO(IP,IQ,ISIM)
     *                              + FVOEX(IP,IQ,ISIM)
 930           CONTINUE
 920        CONTINUE
 910     CONTINUE
 900  CONTINUE
      RETURN
      END
C  /* Deck adfxdm */
      SUBROUTINE ADFXDM(NOSIM,ICI1,IDI1,H2,
     *                  FCOEX,FVOEX,UBOVEC,UDXV,WRK,LWRK)
C
C CALCULATE CONTRIBUTINS TO INACTIVE AND ACTIVE FOCK MATRICES
C WHICH ORIGINATE FROM OCCUPIED-OCCUPIED COULOMN DISTRIBUTIONS
C (ONLY THE EXCHANGE PART OF THE ONE INDEX TRANSFORMED
C  DENSITY CONTRIBUTE)
C
C  FCOEX(P,Q) = FCOEX(P,Q) + SUM(R) (JQ/PR)*UBOVEC(R,J)  (Q OCC)
C
C  FVOEX(P,Q) = FVOEX(P,Q) + SUM(R) (XQ/PR)*UDXV(R,X)   (Q OCC)
C
#include <implicit.h>
C
      DIMENSION H2(NORBT,*)
      DIMENSION FCOEX(NORBT,NORBT,*),FVOEX(NORBT,NORBT,*)
      DIMENSION UBOVEC(NORBT,NORBT,*)
      DIMENSION UDXV(NORBT,NASHT,*)
C
C Used from common blocks:
C   INFORB :
C   INFIND : ISMO,IOBTYP
C
#include <maxorb.h>
#include <maxash.h>
#include <inflin.h>
#include <inforb.h>
#include <infind.h>
#include <infvar.h>
#include <orbtypdef.h>
C
C     Order (C,D) index such that C .ge. D
C     in inactive-active-secondary order (using ISW)
C
      IF (ISW(ICI1) .GE. ISW(IDI1)) THEN
         ICI = ICI1
         IDI = IDI1
      ELSE
         ICI = IDI1
         IDI = ICI1
      END IF
C
C     Find distribution type ITYPCD =
C     1:inactive-inactive  2:active-inactive  3:active-active
C     4:secondary-inactive 5:secondary-active 6:secondary-secondary
C     We only need occupied-occupied distributions, i.e. itypcd .le. 3
C
      ITYPC  = IOBTYP(ICI)
      ITYPD  = IOBTYP(IDI)
      ITYPCD = IDBTYP(ITYPC,ITYPD)
      IF (ITYPCD .GT. 3) RETURN
      IF (ITYPC .EQ. JTACT) NCIW = ICH(ICI)
      IF (ITYPD .EQ. JTACT) NDIW = ICH(IDI)
C
      ICSYM = ISMO(ICI)
      IDSYM = ISMO(IDI)
      IOFFC = IORB(ICSYM) + 1
      IOFFD = IORB(IDSYM) + 1
      ICDSYM = MULD2H(ICSYM,IDSYM)
      DO 100 ISIM = 1,NOSIM
         DO 200 IRSYM = 1,NSYM
            IPSYM = MULD2H(IRSYM,ICDSYM)
            IOFFR = IORB(IRSYM) + 1
            IOFFP = IORB(IPSYM) + 1
            NORBR = NORB(IRSYM)
            NORBP = NORB(IPSYM)
            IF ( (NORBR.EQ.0) .OR. (NORBP.EQ.0) ) GO TO 200
C
            ICRSYM = MULD2H(ICSYM,IRSYM)
            IDRSYM = MULD2H(IDSYM,IRSYM)
            IF ( JWOPSY.EQ.IDRSYM ) THEN
C
               IF ((ITYPCD.EQ.1) .OR. (ITYPCD.EQ.2)) THEN
C
C  FCOEX(P,Q) = FCOEX(P,Q) + SUM(R) (JQ/PR)*UBOVEC(R,J)  (Q OCC)
C                                    QJ/PR
                  CALL DGEMM('N','N',NORBP,1,NORBR,1.D0,
     &                       H2(IOFFP,IOFFR),NORBT,
     &                       UBOVEC(IOFFR,IDI,ISIM),NORBT,1.D0,
     &                       FCOEX(IOFFP,ICI,ISIM),NORBT)
C
               END IF
            END IF
            IF (ICRSYM.EQ.JWOPSY) THEN
               IF ((ITYPCD.EQ.1) .AND. (ICI.NE.IDI)) THEN
C
C  FCOEX(P,Q) = FCOEX(P,Q) + SUM(R) (JQ/PR)*UBOVEC(R,J)  (Q OCC)
C                                    QJ/PR
                  CALL DGEMM('N','N',NORBP,1,NORBR,1.D0,
     &                       H2(IOFFP,IOFFR),NORBT,
     &                       UBOVEC(IOFFR,ICI,ISIM),NORBT,1.D0,
     &                       FCOEX(IOFFP,IDI,ISIM),NORBT)
C
               END IF
            END IF
            IF ( ICRSYM.EQ.JWOPSY ) THEN
               IF ((ITYPCD.EQ.2).OR.(ITYPCD.EQ.3)) THEN
C
C  FVOEX(P,Q) = FVOEX(P,Q) + SUM(R) (XQ/PR)*UDXV(R,X)   (Q OCC)
C
                  CALL DGEMM('N','N',NORBP,1,NORBR,1.D0,
     &                       H2(IOFFP,IOFFR),NORBT,
     &                       UDXV(IOFFR,NCIW,ISIM),NORBT,1.D0,
     &                       FVOEX(IOFFP,IDI,ISIM),NORBT)
C
               END IF
            END IF
            IF (IDRSYM.EQ.JWOPSY) THEN
               IF ((ITYPCD.EQ.3).AND.(ICI.NE.IDI)) THEN
C
C  FVOEX(P,Q) = FVOEX(P,Q) + SUM(R) (XQ/PR)*UDXV(R,X)   (Q OCC)
C
                  CALL DGEMM('N','N',NORBP,1,NORBR,1.D0,
     &                       H2(IOFFP,IOFFR),NORBT,
     &                       UDXV(IOFFR,NDIW,ISIM),NORBT,1.D0,
     &                       FVOEX(IOFFP,ICI,ISIM),NORBT)
C
               END IF
            END IF
 200     CONTINUE
 100  CONTINUE
      RETURN
      END
C  /* Deck adfxdd */
      SUBROUTINE ADFXDD(NSIM,ICI1,IDI1,H2D,FCOCO,FVOCO,
     *                  FCOEX,FVOEX,UBOVEC,UDXV,WRK,LWRK)
C
C CALCULATE CONTRIBUTIONS TO INACTIVE AND ACTIVE FOCK MATRICES
C WHICH ORIGINATE FROM OCCUPIED-OCCUPIED DIRAC DISTRIBUTIONS
C
C F(CV)OEX CONTAIN EXCHANCE CONTRIBUTION
C F(CV)OMU CONTAIN COULOMN CONTRIBUTIONS
C
C  FCOEX(P,Q) = FCOEX(P,Q) + SUM(R) <JP/QR>*UBOVEC(R,J) P OCC Q SEC
C
C  FVOEX(P,Q) = FVOEX(P,Q) + SUM(R) <XP/QR>*UDXV(R,X)  P OCC Q SEC
C
C  FCOCO(P,Q) = FCOCO(P,Q) + SUM(R) <QJ/PR>*UBOVEC(J,R)
C
C  FVOCO(P,Q) = FVOCO(P,Q) + SUM(R) <QY/PR>*UDXV(R,Y)
C
#include <implicit.h>
C
      DIMENSION H2D(NORBT,*)
      DIMENSION FCOCO(NORBT,NORBT,*),FVOCO(NORBT,NORBT,*)
      DIMENSION FCOEX(NORBT,NORBT,*),FVOEX(NORBT,NORBT,*)
      DIMENSION UBOVEC(NORBT,NORBT,*)
      DIMENSION UDXV(NORBT,NASHT,*),WRK(*)
C
      PARAMETER ( D2 = 2.0D0)
C
C Used from common blocks:
C   INFORB :
C   INFIND : ISMO,IOBTYP
C
#include <maxorb.h>
#include <maxash.h>
#include <inflin.h>
#include <inforb.h>
#include <infind.h>
#include <infvar.h>
#include <orbtypdef.h>
C
C     Order (C,D) index such that C .ge. D
C     in inactive-active-secondary order (using ISW)
C
      NDITR = 1
      IF (ISW(ICI1) .GE. ISW(IDI1)) THEN
         ICI = ICI1
         IDI = IDI1
      ELSE
         ICI = IDI1
         IDI = ICI1
         CALL DGETRN(H2D,NORBT,NORBT)
         NDITR = -NDITR
      END IF
C
C     Find distribution type ITYPCD =
C     1:inactive-inactive  2:active-inactive  3:active-active
C     4:secondary-inactive 5:secondary-active 6:secondary-secondary
C     We only need occupied-occupied distributions, i.e. itypcd .le. 3
C
      ITYPC  = IOBTYP(ICI)
      ITYPD  = IOBTYP(IDI)
      ITYPCD = IDBTYP(ITYPC,ITYPD)
      IF (ITYPCD .GT. 3) RETURN
      IF (ITYPC .EQ. JTACT) NCIW = ICH(ICI)
      IF (ITYPD .EQ. JTACT) NDIW = ICH(IDI)
C
      ICSYM = ISMO(ICI)
      IDSYM = ISMO(IDI)
      IOFFC = IORB(ICSYM) + 1
      IOFFD = IORB(IDSYM) + 1
      ICDSYM = MULD2H(ICSYM,IDSYM)
      DO 100 ISIM = 1,NSIM
         DO 200 IRSYM = 1,NSYM
            IPSYM = MULD2H(IRSYM,ICDSYM)
            NORBR = NORB(IRSYM)
            NSSHP = NSSH(IPSYM)
            IF (NORBR.EQ.0) GO TO 200
            IOFFR = IORB(IRSYM) + 1
            IOFFP = IORB(IPSYM) + 1
            NORBP = NORB(IPSYM)
            IF (NORBP.EQ.0) GO TO 200
            NOCCP = NOCC(IPSYM)
            IDRSYM = MULD2H(IDSYM,IRSYM)
            ICRSYM = MULD2H(ICSYM,IRSYM)
            IF ( IDRSYM.EQ.JWOPSY) THEN
               IF (ITYPCD.LE.2) THEN
C
C  FCOCO(P,Q) = FCOCO(P,Q) + SUM(R) <QJ/PR>*UBOVEC(J,R)
C
                  DO 1100 IR=1,NORBR
                     WRK(IR) = UBOVEC(IDI,IOFFR-1+IR,ISIM)
 1100             CONTINUE
                  CALL DGEMM('N','N',NORBP,1,NORBR,1.D0,
     &                       H2D(IOFFP,IOFFR),NORBT,
     &                       WRK,NORBT,1.D0,
     &                       FCOCO(IOFFP,ICI,ISIM),NORBT)
               END IF
            END IF
            IF ( ICRSYM.EQ.JWOPSY) THEN
               IF (ITYPCD.GE.2) THEN
C
C  FVOCO(P,Q) = FVOCO(P,Q) + SUM(R) <QY/PR>*UDXV(R,Y)
C                                    YQ/RP       R,Y
C
                  DO 1200 IR=1,NORBR
                     WRK(IR) = UDXV(IOFFR-1+IR,NCIW,ISIM)
 1200             CONTINUE
                  CALL DGEMM('T','N',NORBP,1,NORBR,1.D0,
     &                       H2D(IOFFR,IOFFP),NORBT,
     &                       WRK,NORBT,1.D0,
     &                       FVOCO(IOFFP,IDI,ISIM),NORBT)
               END IF
            END IF
            IF (NSSHP.EQ.0) GO TO 200
            IF ( IDRSYM.EQ.JWOPSY) THEN
               IF ((ITYPCD.EQ.1) .OR. (ITYPCD.EQ.2)) THEN
C
C  FCOEX(P,Q) = FCOEX(P,Q) + SUM(R) <JP/QR>*UBOVEC(R,J) P OCC Q SEC
C                                    PJ RQ
C
                  CALL DGEMM('T','N',1,NSSHP,NORBR,1.D0,
     &                       UBOVEC(IOFFR,IDI,ISIM),NORBT,
     &                       H2D(IOFFR,IOFFP+NOCCP),NORBT,1.D0,
     &                       FCOEX(ICI,IOFFP+NOCCP,ISIM),NORBT)
               END IF
            END IF
            IF ( ICRSYM.EQ.JWOPSY) THEN
               IF ((ITYPCD.EQ.2).OR.(ITYPCD.EQ.3)) THEN
C
C  FVOEX(P,Q) = FVOEX(P,Q) + SUM(R) <XP/QR>*UDXV(R,X)   (P OCC)
C
                  CALL DGEMM('N','N',NSSHP,1,NORBR,1.D0,
     &                       H2D(IOFFP+NOCCP,IOFFR),NORBT,
     &                       UDXV(IOFFR,NCIW,ISIM),NORBT,0.D0,
     &                       WRK(1),NORBT)
                  DO 350 IQ = IOFFP+NOCCP,IOFFP+NORBP-1
                     IQ1 = IQ - IOFFP - NOCCP + 1
                     FVOEX(IDI,IQ,ISIM) = FVOEX(IDI,IQ,ISIM) + WRK(IQ1)
 350              CONTINUE
               END IF
            END IF
 200     CONTINUE
 100  CONTINUE
      IF (((ITYPCD.EQ.1).OR.(ITYPCD.EQ.3)) .AND. (ICI.NE.IDI) ) THEN
         CALL DGETRN(H2D,NORBT,NORBT)
         NDITR = -NDITR
         DO 110 ISIM = 1,NSIM
            DO 210 IRSYM = 1,NSYM
               IPSYM = MULD2H(IRSYM,ICDSYM)
               NORBR = NORB(IRSYM)
               NSSHP = NSSH(IPSYM)
               IF (NORBR.EQ.0) GO TO 210
               IOFFR = IORB(IRSYM) + 1
               IOFFP = IORB(IPSYM) + 1
               NORBP = NORB(IPSYM)
               IF (NORBP.EQ.0) GO TO 210
               NOCCP = NOCC(IPSYM)
C
               IDRSYM = MULD2H(IDSYM,IRSYM)
               ICRSYM = MULD2H(ICSYM,IRSYM)
               IF (ICRSYM.EQ.JWOPSY) THEN
                  IF (ITYPCD.EQ.1) THEN
C
C  FCOCO(P,Q) = FCOCO(P,Q) + SUM(R) <QJ/PR>*UBOVEC(J,R)
C
                     DO 1101 IR=1,NORBR
                        WRK(IR) = UBOVEC(ICI,IOFFR-1+IR,ISIM)
 1101                CONTINUE
                     CALL DGEMM('N','N',NORBP,1,NORBR,1.D0,
     &                          H2D(IOFFP,IOFFR),NORBT,
     &                          WRK,NORBT,1.D0,
     &                          FCOCO(IOFFP,IDI,ISIM),NORBT)
C
                  END IF
               END IF
               IF (IDRSYM.EQ.JWOPSY) THEN
                  IF (ITYPCD.EQ.3) THEN
C
C  FVOCO(P,Q) = FVOCO(P,Q) + SUM(R) <QY/PR>*UDXV(R,Y)
C                                    YQ/RP       R,Y
C
                     DO 1201 IR=1,NORBR
                        WRK(IR) = UDXV(IOFFR-1+IR,NDIW,ISIM)
 1201                CONTINUE
                     CALL DGEMM('T','N',NORBP,1,NORBR,1.D0,
     &                          H2D(IOFFR,IOFFP),NORBT,
     &                          WRK,NORBT,1.D0,
     &                          FVOCO(IOFFP,ICI,ISIM),NORBT)
C
                  END IF
               END IF
               IF (NSSHP.EQ.0) GO TO 210
               IF (ICRSYM.EQ.JWOPSY) THEN
                  IF (ITYPCD.EQ.1) THEN
C
C  FCOEX(P,Q) = FCOEX(P,Q) + SUM(R) <JP/QR>*UBOVEC(R,J) P OCC Q SEC
C                                    PJ RQ
C
                     CALL DGEMM('T','N',1,NSSHP,NORBR,1.D0,
     &                          UBOVEC(IOFFR,ICI,ISIM),NORBT,
     &                          H2D(IOFFR,IOFFP+NOCCP),NORBT,1.D0,
     &                          FCOEX(IDI,IOFFP+NOCCP,ISIM),NORBT)
C
                  END IF
               END IF
               IF (IDRSYM.EQ.JWOPSY) THEN
                  IF (ITYPCD.EQ.3) THEN
C
C  FVOEX(P,Q) = FVOEX(P,Q) + SUM(R) <XP/QR>*UDXV(R,X)   (P OCC)
C
                     CALL DGEMM('N','N',NSSHP,1,NORBR,1.D0,
     &                          H2D(IOFFP+NOCCP,IOFFR),NORBT,
     &                          UDXV(IOFFR,NDIW,ISIM),NORBT,0.D0,
     &                          WRK(1),NORBT)
                     DO 353 IQ = IOFFP+NOCCP,IOFFP+NORBP-1
                        IQ1 = IQ - IOFFP - NOCCP + 1
                        FVOEX(ICI,IQ,ISIM) = FVOEX(ICI,IQ,ISIM)
     *                                      + WRK(IQ1)
 353                 CONTINUE
                  END IF
               END IF
C
 210        CONTINUE
 110     CONTINUE
      END IF
      IF ( NDITR.LT.0) CALL DGETRN(H2D,NORBT,NORBT)
      RETURN
      END
C  /* Deck adftdd */
      SUBROUTINE ADFTDD(NCSIM,ICI1,IDI1,FVTD,UDVT,H2D)
C
C CALCULATE EXCHANGE CONTRIBUTINS TO ACTIVE FOCK MATRIX WITH A
C TRANSITION DENSITY MATRIX
C
C  FV(P,Q) = FV(P,Q) -0.5*SUM(X,Y)*<XY/QP>*DVT(X,Y)
C THE SUM X,Y ARE TAKEN BY READING IN ALL DIRAC DISTRIBUTIONS
C
#include <implicit.h>
C
      DIMENSION FVTD(NORBT,NORBT,*),UDVT(NASHT,NASHT,*)
      DIMENSION H2D(NORBT,*)
      PARAMETER ( DMP5 = -0.5D0 )
C
C Used from common blocks:
C   INFORB :
C   INFIND : ISMO,IOBTYP
C
#include <maxorb.h>
#include <maxash.h>
#include <inflin.h>
#include <inforb.h>
#include <infind.h>
#include <infvar.h>
#include <orbtypdef.h>
C
C     Order (C,D) index such that C .ge. D
C     in inactive-active-secondary order (using ISW)
C
      NDITR = 1
      IF (ISW(ICI1) .GE. ISW(IDI1)) THEN
         ICI = ICI1
         IDI = IDI1
      ELSE
         ICI = IDI1
         IDI = ICI1
         CALL DGETRN(H2D,NORBT,NORBT)
         NDITR = -NDITR
      END IF
C
C     Find distribution type ITYPCD =
C     1:inactive-inactive  2:active-inactive  3:active-active
C     4:secondary-inactive 5:secondary-active 6:secondary-secondary
C     We only need occupied-occupied distributions, i.e. itypcd .le. 3
C
      ITYPC  = IOBTYP(ICI)
      ITYPD  = IOBTYP(IDI)
      ITYPCD = IDBTYP(ITYPC,ITYPD)
      IF (ITYPCD .NE. 3) RETURN
      IF (ITYPC .EQ. JTACT) NCIW = ICH(ICI)
      IF (ITYPD .EQ. JTACT) NDIW = ICH(IDI)
C
      ICSYM = ISMO(ICI)
      IDSYM = ISMO(IDI)
      ICDSYM = MULD2H(ICSYM,IDSYM)
      IF ( ICDSYM.EQ.JWOPSY) THEN
         DO 100 ISIM = 1,NCSIM
            FAC = DMP5*UDVT(NDIW,NCIW,ISIM)
            CALL DAXPY(NORBT*NORBT,FAC,H2D,1,FVTD(1,1,ISIM),1)
 100     CONTINUE
         IF (NCIW.NE.NDIW) THEN
            CALL DGETRN(H2D,NORBT,NORBT)
            NDITR = -NDITR
            DO 200 ISIM = 1,NCSIM
               FAC = DMP5*UDVT(NCIW,NDIW,ISIM)
               CALL DAXPY(NORBT*NORBT,FAC,H2D,1,FVTD(1,1,ISIM),1)
 200        CONTINUE
         END IF
      END IF
      IF ( NDITR.LT.0) CALL DGETRN(H2D,NORBT,NORBT)
      RETURN
      END
C  /* Deck adftdm */
      SUBROUTINE ADFTDM(NCSIM,ICI1,IDI1,FVTD,DVT,H2)
C
C CALCULATE EXCHANGE CONTRIBUTINS TO ACTIVE FOCK MATRIX WITH A
C TRANSITION DENSITY MATRIX
C
C  FV(P,Q) = FV(P,Q) + SUM(X,Y)*(PQ/XY)*DVT(X,Y)
C THE SUM X,Y ARE TAKEN BY READING IN ALL MULLIKEN DISTRIBUTIONS
C
#include <implicit.h>
C
      DIMENSION FVTD(NORBT,NORBT,*),DVT(NASHT,NASHT,*)
      DIMENSION H2(NORBT,*)
C
C Used from common blocks:
C   INFORB :
C   INFIND : ISMO,IOBTYP
C
#include <maxorb.h>
#include <maxash.h>
#include <inflin.h>
#include <inforb.h>
#include <infind.h>
#include <infvar.h>
#include <orbtypdef.h>
C
C     Order (C,D) index such that C .ge. D
C     in inactive-active-secondary order (using ISW)
C
      NDITR = 1
      IF (ISW(ICI1) .GE. ISW(IDI1)) THEN
         ICI = ICI1
         IDI = IDI1
      ELSE
         ICI = IDI1
         IDI = ICI1
      END IF
C
C     Find distribution type ITYPCD =
C     1:inactive-inactive  2:active-inactive  3:active-active
C     4:secondary-inactive 5:secondary-active 6:secondary-secondary
C     We only need occupied-occupied distributions, i.e. itypcd .le. 3
C
      ITYPC  = IOBTYP(ICI)
      ITYPD  = IOBTYP(IDI)
      ITYPCD = IDBTYP(ITYPC,ITYPD)
      IF (ITYPCD .NE. 3) RETURN
      IF (ITYPC .EQ. JTACT) NCIW = ICH(ICI)
      IF (ITYPD .EQ. JTACT) NDIW = ICH(IDI)
C
      ICSYM = ISMO(ICI)
      IDSYM = ISMO(IDI)
      ICDSYM = MULD2H(ICSYM,IDSYM)
      DO 100 ISIM = 1,NCSIM
         IF ( ICDSYM.EQ.JWOPSY) THEN
            IF (NCIW.NE.NDIW) THEN
               FAC = DVT(NCIW,NDIW,ISIM)+DVT(NDIW,NCIW,ISIM)
            ELSE
               FAC = DVT(NCIW,NDIW,ISIM)
            ENDIF
            CALL DAXPY(NORBT*NORBT,FAC,H2,1,FVTD(1,1,ISIM),1)
         END IF
 100  CONTINUE
      RETURN
      END
C  /* Deck adfxqd */
      SUBROUTINE ADFXQD(NOSIM,NXW,NYW,FXQ,UBOVEC,H2,
     *                  H2X,PVXY,PVYX)
C
C 891006-hjaaj, based on ABAQDI and RESPON.QONEDI
C
C PURPOSE:
C  CALCULATE THE CONTRIBUTIONS TO FXQ  WHICH
C  ARE DETERMINED FROM ACTIVE-ACTIVE DIRAC INTEGRAL DISTRIBUTION
C  <AB,XY> X.ge.Y (DIRAC NOTATION).
C  CONTRIBUTIONS FROM MULLIKEN DISTRIBUTION ARE CALCULATED IN ABAQMU
C
#include <implicit.h>
C
      DIMENSION FXQ(NORBT,NASHT,*)
      DIMENSION UBOVEC(NORBT,NORBT,*),H2(NORBT,*),H2X(NORBT,*)
      DIMENSION PVXY(NASHT,*),PVYX(NASHT,*)
C
C Used from common blocks:
C   INFORB : NORBT, NASHT, ...
C   INFVAR : JWOPSY
C   INFIND : ISW, NSM
C
#include <maxorb.h>
#include <maxash.h>
#include <inflin.h>
#include <inforb.h>
#include <infvar.h>
#include <infind.h>
C
      IXYSYM = MULD2H ( NSM(NXW) , NSM(NYW) )
C
      DO 1000 IOSIM = 1,NOSIM
         CALL DZERO(H2X,NORBT*NORBT)
         DO 1100 IASYM = 1,NSYM
            IBSYM = MULD2H(IASYM,IXYSYM)
            IPSYM = MULD2H(IASYM,JWOPSY)
            NORBA = NORB(IASYM)
            NORBB = NORB(IBSYM)
            NORBP = NORB(IPSYM)
            IAST  = IORB(IASYM) + 1
            IBST  = IORB(IBSYM) + 1
            IPST  = IORB(IPSYM) + 1
C
C ONE INDEX TRANSFORM FIRST INTEGRAL INDEX
C
C  <P~B,XY> =  SUM(A) UBOVEC(P,A)*<AB,XY>
C
            IF ( (NORBA.GT.0) .AND. (NORBB.GT.0) .AND. (NORBP.GT.0) )
     *         CALL DGEMM('N','N',NORBP,NORBB,NORBA,1.D0,
     &                    UBOVEC(IPST,IAST,IOSIM),NORBT,
     &                    H2(IAST,IBST),NORBT,1.D0,
     &                    H2X(IPST,IBST),NORBT)
 1100    CONTINUE
C
C
         DO 1200 IUSYM = 1,NSYM
            IPSYM = MULD2H(IUSYM,JWOPSY)
            IVSYM = MULD2H(IUSYM,IXYSYM)
            NORBP = NORB(IPSYM)
            NASHU = NASH(IUSYM)
            NASHV = NASH(IVSYM)
            IF ((NORBP.GT.0).AND.(NASHU.GT.0).AND.(NASHV.GT.0)) THEN
               IPST = IORB(IPSYM) + 1
               IVST = IORB(IVSYM) + NISH(IVSYM) + 1
               NVWST  = ISW( IVST ) - NISHT
               NUWST  = ISW( IORB(IUSYM) + NISH(IUSYM) +1 ) - NISHT
C
C ADD SUM(V) <V~P,XY> * [UY,VX] TO FXQ(P,U)
C             (V~P)   *  (VU)
C
               CALL DGEMM('T','N',NORBP,NASHU,NASHV,1.D0,
     &                    H2X(IVST,IPST),NORBT,
     &                    PVXY(NVWST,NUWST),NASHT,1.D0,
     &                    FXQ(IPST,NUWST,IOSIM),NORBT)
            END IF
 1200    CONTINUE
         IF (NXW.NE.NYW) THEN
            CALL DZERO(H2X,NORBT*NORBT)
            DO 1400 IASYM = 1,NSYM
               IBSYM = MULD2H(IASYM,IXYSYM)
               IPSYM = MULD2H(IASYM,JWOPSY)
               NORBA = NORB(IASYM)
               NORBB = NORB(IBSYM)
               NORBP = NORB(IPSYM)
               IAST = IORB(IASYM) + 1
               IBST = IORB(IBSYM) + 1
               IPST = IORB(IPSYM) + 1
C
C ONE INDEX TRANSFORM FIRST INTEGRAL INDEX
C
C  <P~B,YX> =  SUM(A) UBOVEC(P,A)*<BA,XY>
C
               IF ( (NORBA.GT.0) .AND. (NORBB.GT.0) .AND. (NORBP.GT.0) )
     *             CALL DGEMM('N','T',NORBP,NORBB,NORBA,1.D0,
     &                        UBOVEC(IPST,IAST,IOSIM),NORBT,
     &                        H2(IBST,IAST),NORBT,1.D0,
     &                        H2X(IPST,IBST),NORBT)
 1400       CONTINUE
C
C
            DO 1500 IUSYM = 1,NSYM
               IPSYM = MULD2H(IUSYM,JWOPSY)
               IVSYM = MULD2H(IUSYM,IXYSYM)
               NORBP = NORB(IPSYM)
               NASHU = NASH(IUSYM)
               NASHV = NASH(IVSYM)
               IF ( (NORBP.GT.0) .AND. (NASHU.GT.0) .AND. (NASHV.GT.0) )
     *                                                             THEN
                  IPST = IORB(IPSYM) + 1
                  IVST = IORB(IVSYM) + NISH(IVSYM) + 1
                  NVWST  = ISW( IVST ) - NISHT
                  NUWST  = ISW( IORB(IUSYM) + NISH(IUSYM) +1 ) - NISHT
C
C ADD SUM(V) <V~P,YX> * [UX,VY] TO FXQ(P,U)
C             (V~P)   *  (VU)
C
                  CALL DGEMM('T','N',NORBP,NASHU,NASHV,1.D0,
     &                       H2X(IVST,IPST),NORBT,
     &                       PVYX(NVWST,NUWST),NASHT,1.D0,
     &                       FXQ(IPST,NUWST,IOSIM),NORBT)
               END IF
 1500       CONTINUE
         END IF
 1000 CONTINUE
C
C END OF ADFXQD
C
      RETURN
      END
C  /* Deck tr1h2m */
      SUBROUTINE TR1H2M(NCSIM,NOSIM,FTQ,H2XAC,FXQ,DTV,PTV,DV,PV,UBO,
     *                  FTV,FXC,FXV,DOFCFV,NODRC,WRK,LFRSAV)
C
C Written 2-3. Nov 1983 by Hans Jorgen Aa. Jensen.
C Rewritten Oct. 1989 hjaaj
C Revisions:
C   6 may '84 ha/hjaaj (for new guga)
C  12-Oct-1984 hjaaj   (act-act rotations)
C  29-Dec-1984 hjaaj   (NCONF.eq.0 option)
C  31-Jan-1985 hjaaj   (ver 2: separated conf. and orb. trial vectors;
C                       also changed NCONF.eq.0 option to NCSIM.lt.0)
C   7-Feb-1985 hjaaj   (corrected some errors)
C  22-Aug-1988 hjaaj   (removed FXCAC, FXIJ, FXJI)
C  19-Sep-1989 hjaaj   (squared PV,
C                       use matrix multiply for Fock Q matrices)
C
C Purpose:
C   Construct transformed 2-electron integrals (*u/vx) .
C   To construct FXQ and FTQ :
C     FXQ(p,t) = SUM(uvx): PV(tu:vx)*(pXuX:vXxX)
C     FTQ(p,t) = SUM(uvx): PTV(tu:vx)*(pu:vx).
C   If NCSIM .lt. 0, only calculate what needed for orbital
C     Hessian transformation (FXQ), that is, do not calculate H2XAC, FTQ
C   NOTE: FXQ, FTQ, and H2XAC are added to input
C
C Input:
C   Two-electron density matrices PV and PTV
C   Kappa as a full matrix in UBO
C   If (DOFCFV) one-electron density matrices DV and DTV
C
C Output:
C   "FXQ" matrix in FXQ(NORBT,NASHT,NOSIM)
C   "FTQ" matrix in FTQ(NORBT,NASHT,NCSIM)
C   H2XAC One-index transformed integrals over active orbitals
C   If (DOFCFV) FXC,FXV,FTV Fock matrices
C
C Scratch:
C   WRK(KFRSAV:LFRSAV)
C
C
C Externals:
C   DDOT, DZERO, MOLLAB, QTRACE
#include <implicit.h>
      DIMENSION PV(NNASHX,*),PTV(NNASHX,NNASHX,*),UBO(N2ORBX,*)
      DIMENSION H2XAC(NNASHX,NNASHX,*)
      DIMENSION FXQ(NASHT*NORBT,*),FTQ(NASHT*NORBT,*)
      DIMENSION DV(NNASHX), DTV(NNASHX,*)
      DIMENSION FTV(N2ORBX,*), FXC(N2ORBX,*), FXV(N2ORBX,*)
      DIMENSION WRK(LFRSAV)
      LOGICAL   DOFCFV, NODRC
C
C -- local constants
C
      PARAMETER ( D2 = 2.0D0 , D4 = 4.0D0 )
      PARAMETER ( THRH2X=1.D-10 )
C
C Used from common blocks
C   INFINP : FLAG(*)
C   INFORB : NSYM,NNASHX,NORBT,...
C   INFIND : IROW(*),ICH(*),IOBTYP(*),?
C   INFVAR : JWOPSY
C   INFLIN : IPRLIN
C   INFPRI : P6FLAG(*),?
C CBGETDIS : IADH2X
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <infind.h>
#include <infvar.h>
#include <inflin.h>
#include <infpri.h>
#include <cbgetdis.h>
C
C     Local variables
C
      DIMENSION NEEDMU(6)

#include <orbtypdef.h>

      CALL QENTER('TR1H2M')
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LFRSAV
C
C     Check input
C
      IF (.NOT.DOFCFV .AND. NASHT .LE. 1) GO TO 9999
C     (FXQ and FTQ zero matrix if NASHT .le. 1)
C
      IF (DOFCFV .AND. NODRC) THEN
         WRITE (LUPRI,'(/A)')
     &      ' TR1H2M error, DOFCFV and NODRC not compatible'
         CALL QTRACE(LUPRI)
         CALL QUIT('TR1H2M error, DOFCFV and NODRC both true')
      END IF
C
C     Set NEEDMU array
C
      DO 10 I = 1,6
         NEEDMU(I) = 0
   10 CONTINUE
      NEEDMU(3) = 1
C     ... active-active distributions always needed
      IF (NOSIM .GT. 0) THEN
         IF (NODRC) THEN
C           ... all distributions involving active orbitals needed
C               for FXQ
            NEEDMU(2) = 1
            NEEDMU(5) = 1
         END IF
         IF (DOFCFV) THEN
C           ... inactive-occupied needed for FXC, FXV
            NEEDMU(1) = 1
            NEEDMU(2) = 1
         END IF
      END IF
C
C     Allocate work memory
C
      CALL MEMGET('REAL',KPVCD ,N2ASHX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KH2CD ,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KH2XCD,N2ORBX,WRK,KFREE,LFREE)
      IF (DOFCFV) THEN
         IF (NCSIM .GT. 0) THEN
            CALL MEMGET('REAL',KUDTV,NCSIM*N2ASHX,WRK,KFREE,LFREE)
            DO 40 ICSIM = 1,NCSIM
               JUDTV = KUDTV + (ICSIM-1)*N2ASHX
               CALL DSPTSI(NASHT,DTV(1,ICSIM),WRK(JUDTV))
   40       CONTINUE
         END IF
         IF (NOSIM .GT. 0) THEN
            CALL MEMGET('REAL',KUDXV,NOSIM*NASHT*NORBT,WRK,KFREE,LFREE)
            IF (NASHT .GT. 0) THEN
               CALL MEMGET('REAL',KUDV,N2ASHX,WRK,KFREE,LFREE)
               CALL DSPTSI(NASHT,DV,WRK(KUDV))
               CALL TR1DV(NOSIM,WRK(KUDV),UBO,WRK(KUDXV))
               CALL MEMREL('TR1H2M.UDV',WRK,KFRSAV,KUDV,KFREE,LFREE)
            END IF
         END IF
      END IF
C
      JDIST = 0
      IF (P6FLAG(36) .OR. IPRLIN .GE. 40) THEN
         WRITE (LUPRI,'(///A)')
     &      ' Test output of Mulliken distributions from TR1H2M.'
      END IF
C
C ****************************************************************
C     Loop over Mulliken distributions allowed in NEEDMU(6)
C
      IDIST = 0
  100 CALL NXTH2M(IC,ID,WRK(KH2CD),NEEDMU,WRK,KFREE,LFREE,IDIST)
      IF (IDIST .LT. 0) GO TO 800
C     ... if idist .lt. 0 then no more distributions
      IF (P6FLAG(36) .OR. IPRLIN .GE. 40) THEN
         JDIST = JDIST + 1
         WRITE (LUPRI,'(/A,I5,A,2I5)')
     &      ' Mulliken distribution no.',JDIST,', IC and ID:',IC,ID
         WRITE (LUPRI,'(A,I5)') ' IDIST form NXTH2M :',IDIST
         CALL OUTPUT(WRK(KH2CD),1,NORBT,1,NORBT,NORBT,NORBT,-1,LUPRI)
      END IF
C
C     IF DOFCFV then call Fock matrix routines
C
         IF (DOFCFV) THEN
            IF (NOSIM .GT. 0) THEN
               CALL ADFXDM(NOSIM,IC,ID,WRK(KH2CD),FXC,FXV,
     &                     UBO,WRK(KUDXV),WRK(KFREE),LFREE)
            END IF
            IF (NCSIM .GT. 0) THEN
               CALL ADFTDM(NCSIM,IC,ID,FTV,WRK(KUDTV),WRK(KH2CD))
            END IF
C           CALL ADFXDM(NOSIM,ICI1,IDI1,H2,
C    *                  FCOEX,FVOEX,UBOVEC,UDXV,WRK,LWRK)
C           CALL ADFTDM(NCSIM,ICI1,IDI1,FVTD,DVT,H2)
         END IF
C
C     If NASHT.le.1 then no FTQ and FXQ is zero matrix, thus:
C
      IF (NASHT .LE. 1) GO TO 100
C
C     Find distribution type (only distributions with at least
C     one active index are needed):
C
         ITYPC  = IOBTYP(IC)
         ITYPD  = IOBTYP(ID)
         ITYPCD = IDBTYP(ITYPC,ITYPD)
      IF (ITYPC .NE. JTACT .AND. ITYPD .NE. JTACT) GO TO 100
      IF (NOSIM.LE.0 .AND. ITYPCD.NE.JTACAC) GO TO 100
C     ... only distributions with at least one active index needed
C     ... only active-active distributions are needed for FTQ
C         (previous two executable statements should never become
C          active because of NEEDMU)
         ISWAP = 0
         IF (ITYPCD .EQ. JTACAC) THEN
            IF (ISW(IC) .LT. ISW(ID)) ISWAP = 1
C           ... we want NCW .ge. NDW for active-active
         ELSE
            IF (IOBTYP(IC) .EQ. JTACT) ISWAP = 1
C           ... we want IC to be non-active index
         END IF
         IF (ISWAP .EQ. 1) THEN
            ISWAP = IC
            IC = ID
            ID = ISWAP
         END IF
         IF (ITYPCD .NE. JTACAC) THEN
            ICXSYM = MULD2H(JWOPSY,ISMO(IC))
            IF (NASH(ICXSYM).EQ.0) GO TO 100
C     ^-------------------------------------
C           if no active orbitals in ICXSYM (new symmetry after
C           one-index transformation of IC), we don't need this
C           distribution
         END IF
C
         ICSYM  = ISMO(IC)
         IDSYM  = ISMO(ID)
         ICDSYM = MULD2H(ICSYM,IDSYM)
         IF (ITYPCD .EQ. JTACAC) THEN
            NCW    = ICH(IC)
            NDW    = ICH(ID)
            NCDW   = IROW(NCW) + NDW
C
C
C           Add contributions from this active-active distribution
C           to transition Fock Q matrix (FTQ):
C
            DO 299 ICSIM = 1,NCSIM
               CALL SIRPVD(NCW,NDW,WRK(KPVCD),PTV(1,1,ICSIM),1)
               IF (NCW .NE. NDW) CALL DSCAL(N2ASHX,D2,WRK(KPVCD),1)
C              ... C,D and D,C contributions are identical
               CALL ADDFQ(NCW,NDW,FTQ(1,ICSIM),WRK(KH2CD),WRK(KPVCD),1)
C              CALL ADDFQ(NXW,NYW,FQ,H2XY,PVXY,JH2SYM)
C
  299       CONTINUE
C           end do 299 icsim = 1,ncsim
C
C
            IF (NOSIM .GT. 0) THEN
               CALL SIRPVD(NCW,NDW,WRK(KPVCD),PV,1)
               IF (NCW .NE. NDW) CALL DSCAL(N2ASHX,D2,WRK(KPVCD),1)
C              ... C,D and D,C contributions are identical
C
            DO 399 IOSIM = 1,NOSIM
               CALL DZERO(WRK(KH2XCD),N2ORBX)
               CALL TR1UH1(UBO(1,IOSIM),WRK(KH2CD),WRK(KH2XCD),ICDSYM)
C              CALL TR1UH1(UBOVEC,H1,H1X,IH1SYM)
C
C We now have the one-index transformed integrals in H2XCD(*,*):
C (AX BX / C D) = (AX B / C D) + (A BX / C D)
C (A  BX / C D) = SUM(r)  BORB(b,r) * H2CD(a,r)
C
C              Add contributions to FXQ
C
               CALL ADDFQ(NCW,NDW,FXQ(1,IOSIM),
     &                    WRK(KH2XCD),WRK(KPVCD),JWOPSY)
C              CALL ADDFQ(NXW,NYW,FQ,H2XY,PVXY,JH2SYM)
C
C              If (calculate h2xac) then
C              add active-active elements to H2XAC
C
               IF (NCSIM .GE. 0) THEN
                  IUVSYM = MULD2H(ICDSYM,JWOPSY)
                  CALL ADH2AC(H2XAC(1,NCDW,IOSIM),WRK(KH2XCD),IUVSYM)
               END IF
C
  399       CONTINUE
C           end do 399 iosim = 1,nosim
            END IF
C
         END IF
C        end of if ( (cd) active-active ) then ... end if
C
C        ... If NODRC then
C            now the (p~ v~ | x y) contribution to FXQ
C
C        ... (CD) inactive-active, active-active
C            or secondary-active distribution
C
         IF (NODRC .AND.
     &      (ITYPC.EQ.JTACT .OR. ITYPD.EQ.JTACT)) THEN
            DO 499 IOSIM = 1,NOSIM
               CALL DZERO(WRK(KH2XCD),N2ORBX)
               CALL TR1UH1(UBO(1,IOSIM),WRK(KH2CD),WRK(KH2XCD),ICDSYM)
C              CALL TR1UH1(UBOVEC,H1,H1X,IH1SYM)
               CALL ADFXQM(IC,ID,FXQ(1,IOSIM),WRK(KH2XCD),PV,
     &                     WRK(KFREE),LFREE)
C
  499       CONTINUE
C           end do 499 iosim = 1,nosim
C
         END IF
C
C        Go to 100 to get next needed Mulliken distribution
C
      GO TO 100
C
C     arrive at 800 when finished with all needed Mulliken distributions
C
  800 CONTINUE
      IF (P6FLAG(36) .OR. IPRLIN .GE. 40) THEN
         WRITE (LUPRI,'(//A/,2(/A,I5))')
     &     ' End of test output of Mulliken distributions from TR1H2M.',
     &     ' Total number of distributions treated  :',JDIST,
     &     ' Total number of distributions (NNORBX) :',NNORBX
      END IF
C
C
C ****************************************************************
C
C  Now H2XAC(ij,kl,iosim) = (it jt / k l),
C  finish H2XAC by symmetrizing
C  (remember (it jt / kt lt) = (it jt / k l) + (kt lt / i j) ).
C
C  If requested, print H2XAC
C
      IF (NCSIM.GE.0) THEN
      DO 1800 IOSIM = 1,NOSIM
C
         DO 1720 IJ = 1,NNASHX
            DO 1710 KL = 1,IJ
               X = H2XAC(IJ,KL,IOSIM) + H2XAC(KL,IJ,IOSIM)
               H2XAC(KL,IJ,IOSIM) = X
               H2XAC(IJ,KL,IOSIM) = X
 1710       CONTINUE
 1720    CONTINUE
         IF (P6FLAG(28) .OR. IPRLIN.GE.28) THEN
            WRITE (LUPRI,1020) IOSIM,NOSIM
            CALL OUTPUT(H2XAC(1,1,IOSIM),1,NNASHX,1,NNASHX,
     *                  NNASHX,NNASHX,-1,LUPRI)
         END IF
 1020 FORMAT (/' H2XAC matrix (no.',I3,' of',I3,')')
C
 1800 CONTINUE
      END IF
C
C *** IADH2X .lt. 0 means that H2XAC is in core.
C
      IADH2X = -1
C
C *** end of subroutine TR1H2M
C
 9999 IF (KFRSAV .NE. KFREE)
     &   CALL MEMREL('TR1H2M',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('TR1H2M')
      RETURN
      END
C  /* Deck tr1h2d */
      SUBROUTINE TR1H2D(NCSIM,NOSIM,DTV,DV,PV,UBO,
     &                  FTV,FXC,FXV,FXQ,DOFCFV,WRK,LFRSAV)
C
C Jan 90 Hans Joergen Aa. Jensen
C
C Purpose:
C   Calculate FXQ(:,1:NOSIM)
C   If (DOFCFV) calculate FTV(:,1:NCSIM), FXC(:,1:NOSIM),
C      and FXV(:,1:NOSIM)
C
#include <implicit.h>
#include <dummy.h>
      DIMENSION DTV(NNASHX,*), DV(*),PV(NNASHX,*)
      DIMENSION UBO(N2ORBX,*), FXQ(NASHT*NORBT,*)
      DIMENSION FTV(N2ORBX,*), FXC(N2ORBX,*), FXV(N2ORBX,*)
      DIMENSION WRK(LFRSAV)
      LOGICAL   DOFCFV
C
C Used from common blocks:
C  INFORB : NASHT,...
C  INFIND : IOBTYP(*),?
C  INFLIN : IPRLIN
C  INFPRI : P6FLAG(*)
C
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <inforb.h>
#include <infind.h>
#include <inflin.h>
#include <inftra.h>
#include <infpri.h>
C
C     Local variables
C
      LOGICAL   DOFXQ
      DIMENSION NEEDDI(6)
C
      CALL QENTER('TR1H2D')
C
      DOFXQ = (NASHT.GT.1)
C
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LFRSAV
      IF (.NOT.DOFCFV .AND. .NOT.DOFXQ) GO TO 9999
C
C     Set NEEDDI array
C
      IF (DOFCFV) THEN
C     ..OCCUPIED-OCCUPIED distributions needed
         NEEDDI(1) = 1
         NEEDDI(2) = 1
      ELSE
         NEEDDI(1) = 0
         NEEDDI(2) = 0
      END IF
      NEEDDI(3) = 1
      NEEDDI(4) = 0
      NEEDDI(5) = 0
      NEEDDI(6) = 0
C
C     Allocate work memory
C
      CALL MEMGET('REAL',KFCOCO,N2ORBX*NOSIM,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KFVOCO,N2ORBX*NOSIM,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KUDXV ,NORBT*NASHT*NOSIM,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KPVCD ,N2ASHX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KPVDC ,N2ASHX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KUDV  ,N2ASHX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KH2CD ,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KH2XCD,N2ORBX,WRK,KFREE,LFREE)
C
      IF ( DOFCFV ) THEN
         IF (NCSIM .GT. 0) THEN
            CALL MEMGET('REAL',KUDTV,NCSIM*N2ASHX,WRK,KFREE,LFREE)
            DO 40 ICSIM = 1,NCSIM
               JUDTV = KUDTV + (ICSIM-1)*N2ASHX
               CALL DSPTSI(NASHT,DTV(1,ICSIM),WRK(JUDTV))
   40       CONTINUE
         END IF
         IF ( NOSIM .GT. 0 ) THEN
            CALL DZERO(WRK(KFCOCO),N2ORBX*NOSIM)
            CALL DZERO(WRK(KFVOCO),N2ORBX*NOSIM)
         END IF
      END IF
      IF ( (NASHT.GT.0) .AND. (NOSIM.GT.0) ) THEN
         CALL DSPTSI(NASHT,DV,WRK(KUDV))
         CALL TR1DV(NOSIM,WRK(KUDV),UBO,WRK(KUDXV))
      END IF
C
      IF (P6FLAG(37) .OR. IPRLIN .GE. 40) THEN
         JDIST = 0
         WRITE (LUPRI,'(///A)')
     &      ' Test output of Dirac distributions from TR1H2D.'
      END IF
C
C ****************************************************************
C     Loop over DIRAC distributions allowed in NEEDDI(6)
C
      IDIST = 0
      LUINTD = -1
      IF (.NOT. NEWTRA) CALL GPOPEN(LUINTD,'DRCINT','OLD',' ',
     &                              'UNFORMATTED',IDUMMY,.FALSE.)
  100 CALL NXTH2D(IC,ID,WRK(KH2CD),NEEDDI,LUINTD,WRK,KFREE,LFREE,IDIST)
      IF (IDIST .LT. 0) GO TO 800
C     ... if idist .lt. 0 then no more distributions
C
      IF (P6FLAG(37) .OR. IPRLIN .GE. 40) THEN
         JDIST = JDIST + 1
         WRITE (LUPRI,'(/A,I5,A,2I5)')
     &      ' Dirac distribution no.',JDIST,', IC and ID:',IC,ID
         CALL OUTPUT(WRK(KH2CD),1,NORBT,1,NORBT,NORBT,NORBT,-1,LUPRI)
      END IF
C
C     Find distribution type (only distributions with at least
C     one active index are needed):
C
         ICSYM  = ISMO(IC)
         IDSYM  = ISMO(ID)
         ICDSYM = MULD2H(ICSYM,IDSYM)
         ITYPC  = IOBTYP(IC)
         ITYPD  = IOBTYP(ID)
         IF (NOSIM.GT.0) THEN
            IF ( DOFXQ .AND. (ITYPC.EQ.JTACT) .AND. (ITYPD.EQ.JTACT) )
     *         THEN
               NCW    = ICH(IC)
               NDW    = ICH(ID)
               CALL SIRPVD(NDW,NCW,WRK(KPVCD),PV,2)
               IF ( IC.NE.ID) THEN
                  CALL MTRSP(NASHT,NASHT,WRK(KPVCD),NASHT,
     &                       WRK(KPVDC),NASHT)
               END IF
               CALL ADFXQD(NOSIM,NCW,NDW,FXQ,UBO,WRK(KH2CD),
     *                     WRK(KH2XCD),WRK(KPVCD),WRK(KPVDC))
C              CALL ADFXQD(NOSIM,NXW,NYW,FXQ,UBOVEC,H2,
C    *                     H2X,PVXY,PVYX)
            END IF
            IF ( DOFCFV ) THEN
               CALL ADFXDD(NOSIM,IC,ID,WRK(KH2CD),WRK(KFCOCO),
     *                     WRK(KFVOCO),FXC,FXV,UBO,WRK(KUDXV),
     *                     WRK(KFREE),LFREE)
C              CALL ADFXDD(NSIM,ICI1,IDI1,H2D,FCOCO,FVOCO,
C    *                     FCOEX,FVOEX,UBOVEC,UDXV,WRK,LWRK)
            END IF
         END IF
         IF ( (NCSIM.GT.0) .AND. DOFCFV )
     *      CALL ADFTDD(NCSIM,IC,ID,FTV,WRK(KUDTV),WRK(KH2CD))
C           CALL ADFTDD(NSIM,ICI1,IDI1,FVTD,DVT,H2D)
      GO TO 100
 800  CONTINUE
      IF (LUINTD .GT. 0) CALL GPCLOSE(LUINTD,'KEEP')
C
      IF (P6FLAG(37) .OR. IPRLIN .GE. 40) THEN
         WRITE (LUPRI,'(//A/,2(/A,I5))')
     &     ' End of test output of Dirac distributions from TR1H2D.',
     &     ' Total number of distributions treated  :',JDIST,
     &     ' Total number of distributions (NNORBX) :',NNORBX
      END IF
C
      IF ( (NOSIM.GT.0) .AND. DOFCFV) THEN
         CALL FXDFIN(NOSIM,WRK(KFCOCO),WRK(KFVOCO),FXC,FXV)
      END IF
C
C *** end of subroutine TR1H2D
C
 9999 IF (KFRSAV .NE. KFREE)
     &   CALL MEMREL('TR1H2D',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('TR1H2D')
      RETURN
      END
C  /* Deck tr1dv */
      SUBROUTINE TR1DV(NOSIM,UDV,UBOVEC,UDXV)
C
C PURPOSE:
C   CALCULATE MODIFIED DENSITY MATRICES
C
C   UDXV(P,Y) = SUM(U) UBOVEC(P,U) * UDV(U,Y)
C
C   THE MATRIX UDXV IS STORED WITHOUT USING SYMMETRY
C
#include <implicit.h>
#include <maxash.h>
#include <maxorb.h>
C
      DIMENSION UDV(NASHT,*),UBOVEC(NORBT,NORBT,*)
      DIMENSION UDXV(NORBT,NASHT,*)
C
#include <priunit.h>
#include <inforb.h>
#include <infind.h>
#include <inflin.h>
#include <infvar.h>
C
      CALL DZERO(UDXV,NORBT*NASHT*NOSIM)
      DO 100 IOSIM=1,NOSIM
         DO 1200 IUSYM = 1,NSYM
            IPSYM = MULD2H(IUSYM,JWOPSY)
            NORBIP = NORB(IPSYM)
            NASHIU = NASH(IUSYM)
            IF ((NORBIP.GT.0).AND.(NASHIU.GT.0)) THEN
               IPST = IORB(IPSYM) + 1
               NISHIU = NISH(IUSYM)
               IUST   = IORB(IUSYM) + 1 + NISHIU
               NUWST  = ISW( IORB(IUSYM) + NISH(IUSYM) +1 ) - NISHT
C
C   UDXV(P,Y) = SUM(U) UBOVEC(P,U) * UDV(U,Y)
C
               CALL DGEMM('N','N',NORBIP,NASHIU,NASHIU,1.D0,
     &                    UBOVEC(IPST,IUST,IOSIM),NORBT,
     &                    UDV(NUWST,NUWST),NASHT,1.D0,
     &                    UDXV(IPST,NUWST,IOSIM),NORBT)
            END IF
 1200    CONTINUE
 100  CONTINUE
      IF ( IPRLIN.GT.70 ) THEN
         DO 3010 ISIM = 1, NOSIM
            WRITE(LUPRI,'(/I5,A,I5)') ISIM,' UDXV matrix out of',NOSIM
            CALL OUTPUT(UDXV(1,1,ISIM),1,NORBT,1,NASHT,
     *                  NORBT,NASHT,1,LUPRI)
 3010    CONTINUE
      END IF
      RETURN
      END
