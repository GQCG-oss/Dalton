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
#ifdef REVISION_LOG
===========================================================================
Revision 1.2  2001/01/17 10:41:33  vebjornb
Calls to *MPA*B* in arhpack.F have been replaced with DGEMM calls

Revision 1.6  2000/05/05 11:15:57  hjj
bugfix in updgrd, ci diagonal was not calculated because LUIT2 .lt. 0
has until today meant do not calculate CI diagonal. Now LUIT2 is temp.
set to zero to tell CIDIAG not to save diagonal on file.

Revision 1.4  2000/05/01 10:49:20  hjj
If FULODG then call MAKE_KLWOP to get KLWOP for full orbital diagonal.
(KLWOP is not in infvar.h any more)
===========================================================================
#endif
C  /* Deck cigrad */
      SUBROUTINE CIGRAD(CREF,FCAC,H2AC,INDXCI,GCI,EMCACT,WRK,LFREE)
C
C Written Nov-1983/May-1984 by Hans Agren/Hans Jorgen Aa. Jensen
C Rewritten 15-Jan-1988 hjaaj
C
C Purpose:
C   Compute the total energy and the CI gradient at the beginning
C   of each macroiteration. (The inactive energy, EMCMY, is not
C   included here as it would just be added and subtracted again.)
C
C Output:
C   Active MCSCF energy stored in EMCACT.
C   CI gradient is stored in GCI.
C
C Scratch:
C   WRK(LFREE)
C
C Externals:
C   CISIGC
C
C
#include <implicit.h>
      DIMENSION CREF(*),FCAC(*),H2AC(*), INDXCI(*),GCI(*),WRK(*)
      PARAMETER ( D2 = 2.0D0 )
C
C Used from common blocks:
C   INFVAR : NCONF
C CBGETDIS : DISTYP,IADINT,IADH2
C
#include <maxorb.h>
#include <infvar.h>
#include <cbgetdis.h>

#ifdef SIRGRAD_DEBUG
#include <priunit.h>
#endif
C
      CALL QENTER('CIGRAD')
C
C *** 1 *** GCI(I) = SUM(J) OF LCI(I,J)*CREF(J)
C           calculate CI matrix times reference vector
C
      DISTYP = 1
      IADINT = IADH2
      CALL CISIGC(1,CREF,GCI,NCONF,FCAC,H2AC,INDXCI,WRK,LFREE)
C     CALL CISIGC(NSIM,BCVECS,SCVECS,LSCVEC,FCAC,H2AC,INDXCI,WRK,LFREE)
C
C
C *** 2 *** Calculate total energy:
C
      EMCACT = DDOT(NCONF,CREF,1,GCI,1)
C
C *** 3 *** Finish the CI gradient calculation:
C           GCI(I) = 2*( GCI(I) - EMCACT * CREF(I) )
C
      DO 100 II = 1,NCONF
         GCI(II) = D2*( GCI(II) - EMCACT * CREF(II) )
  100 CONTINUE
C
C *** End of subroutine CIGRAD
C
      CALL QEXIT('CIGRAD')
      RETURN
      END
C  /* Deck fckden */
      SUBROUTINE FCKDEN(GETDC,GETDV,DCAO,DVAO,CMO,DV,WRK,LFRSAV)
C
C Jan. 1990 Hans Joergen Aa. Jensen
C l.r. 1997 hjaaj (d.m. now square matrices)
C
C Backtransform one-electron density matrices for calculation
C of Fock matrices in AO basis ( contravariant transformation ).
C
C Input:
C   GETDC   if true calculate DCAO
C   GETDV   if true calculate DVAO
C   CMO(*)  molecular orbital coefficients
C   DV(*)   active part of one-electron density matrix
C           (over MOs)
C Scratch:
C   WRK(LFRSAV)
C
#include <implicit.h>
      DIMENSION CMO(*),DV(NNASHX,*),DCAO(NBAST,*),DVAO(NBAST,NBAST,*), 
     &          WRK(*)
      LOGICAL   GETDC, GETDV
C
      PARAMETER (D2 = 2.0D0)
C
C Used from common blocks:
C  INFORB : NSYM,MULD2H(8,8),NASHT,NNBASX,N2BASX,...
C  INFVAR : JWOPSY
C  INFPRI : P6FLAG()
C
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infvar.h>
#include <infpri.h>
C dftcom.h : DFT_SPINDNS
#include <dftcom.h>
C
C
      CALL QENTER('FCKDEN')
C
C
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LFRSAV
C
C     ************************************
C     DVAO matrix
C
      IF (GETDV) THEN
         CALL DZERO(DVAO,N2BASX)
      IF (NASHT .GT. 0) THEN
         CALL MEMGET('REAL',KDAO1,NASHT*NBAST,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KUDV ,N2ASHX,WRK,KFREE,LFREE)
         CALL DSPTSI(NASHT,DV,WRK(KUDV))
         DO 2000 ISYM = 1,NSYM
            JSYM  = MULD2H(ISYM,JWOPSY)
            NASHI = NASH(ISYM)
            NASHJ = NASH(JSYM)
         IF (NASHI .EQ. 0 .OR. NASHJ .EQ. 0) GO TO 2000
            NORBI = NORB(ISYM)
            NBASI = NBAS(ISYM)
            NORBJ = NORB(JSYM)
            NBASJ = NBAS(JSYM)
            JCMO  = ICMO(JSYM) + 1 + NISH(JSYM)*NBASJ
            JUDV  = KUDV + IASH(ISYM)*NASHT + IASH(JSYM)
            CALL DGEMM('N','N',NBASJ,NASHI,NASHJ,1.D0,
     &                 CMO(JCMO),NBASJ,
     &                 WRK(JUDV),NASHT,0.D0,
     &                 WRK(KDAO1),NBASJ)
            JCMO  = ICMO(ISYM) + 1 + NISH(ISYM)*NBASI
            CALL DGEMM('N','T',NBASJ,NBASI,NASHI,1.D0,
     &                 WRK(KDAO1),NBASJ,
     &                 CMO(JCMO),NBASI,0.D0,
     &                 DVAO(IBAS(JSYM)+1,IBAS(ISYM)+1,1),NBAST) !JT
 2000    CONTINUE
      CALL MEMREL('FCKDEN DVAO',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      END IF
      END IF  ! IF (GETDV) THEN
C
C     ************************************
C     DVSAO matrix
C
C JT Aug 09 beg
C do spin density matrix
      IF (GETDV .AND. DFT_SPINDNS) THEN
         CALL DZERO(DVAO(1,1,2),N2BASX)
      IF (NASHT .GT. 0) THEN
         CALL MEMGET('REAL',KDAO1,NASHT*NBAST,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KUDV ,N2ASHX,WRK,KFREE,LFREE)
         CALL DSPTSI(NASHT,DV(1,2),WRK(KUDV))
         DO 2010 ISYM = 1,NSYM
            JSYM  = MULD2H(ISYM,JWOPSY)
            NASHI = NASH(ISYM)
            NASHJ = NASH(JSYM)
         IF (NASHI .EQ. 0 .OR. NASHJ .EQ. 0) GO TO 2010
            NORBI = NORB(ISYM)
            NBASI = NBAS(ISYM)
            NORBJ = NORB(JSYM)
            NBASJ = NBAS(JSYM)
            JCMO  = ICMO(JSYM) + 1 + NISH(JSYM)*NBASJ
            JUDV  = KUDV + IASH(ISYM)*NASHT + IASH(JSYM)
            CALL DGEMM('N','N',NBASJ,NASHI,NASHJ,1.D0,
     &                 CMO(JCMO),NBASJ,
     &                 WRK(JUDV),NASHT,0.D0,
     &                 WRK(KDAO1),NBASJ)
            JCMO  = ICMO(ISYM) + 1 + NISH(ISYM)*NBASI
            CALL DGEMM('N','T',NBASJ,NBASI,NASHI,1.D0,
     &                 WRK(KDAO1),NBASJ,
     &                 CMO(JCMO),NBASI,0.D0,
     &                 DVAO(IBAS(JSYM)+1,IBAS(ISYM)+1,2),NBAST)
 2010    CONTINUE
         CALL MEMREL('FCKDEN DSAO',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      END IF  ! IF (NASHT .GT. 0) THEN
      END IF  ! IF (GETDV .AND. DFT_SPINDNS) THEN
C JT Aug 09 end
C
C
C     ************************************
C     DCAO matrix
C
      IF (GETDC) THEN
         CALL DZERO(DCAO,N2BASX)
      IF (NISHT .GT. 0 .AND. JWOPSY .EQ. 1) THEN
         DO 4000 ISYM = 1,NSYM
            NISHI = NISH(ISYM)
         IF (NISHI .GT. 0) THEN
            NBASI = NBAS(ISYM)
            JCMO  = ICMO(ISYM) + 1
            CALL DGEMM('N','T',NBASI,NBASI,NISHI,1.D0,
     &                 CMO(JCMO),NBASI,
     &                 CMO(JCMO),NBASI,0.D0,
     &                 DCAO(IBAS(ISYM)+1,IBAS(ISYM)+1),NBAST)
         END IF
 4000    CONTINUE
C
C **     multiply with DC(i,i) = 2.0 factor to get final result
C        (if transition density matrix, then DCAO must be
C         multiplied with the <L | R> overlap outside this routine)
C
         CALL DSCAL(N2BASX,D2,DCAO,1)
      END IF
      END IF
C
C
C
      IF (P6FLAG(18)) THEN
        IF (GETDC) THEN
          WRITE(LUPRI,'(/A)')
     &    ' FCKDEN: DCAO = Density matrix, inactive part (AO-basis)'
          CALL OUTPUT(DCAO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
        END IF
        IF (GETDV) THEN
          WRITE(LUPRI,'(/A)')
     &    ' FCKDEN: DVAO = Density matrix,   active part (AO-basis)'
          CALL OUTPUT(DVAO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
        IF (DFT_SPINDNS) THEN
          WRITE(LUPRI,'(/A)')
     &    ' FCKDEN: DSAO = Spin density matrix, active part (AO-basis)'
          CALL OUTPUT(DVAO(1,1,2),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
        END IF
        END IF
      END IF
C
C
      CALL QEXIT('FCKDEN')
      RETURN
C
C *** end of subroutine FCKDEN
C
       END
C  /* Deck gtdmso */
      SUBROUTINE GTDMSO(UDV,CMO,DCAO,DVAO,WRK) 
C
C     Compute inactive (DCAO) and active (DVAO) density matrices in AO basis
C     from molecular orbitals given in CMO matrix and the occupation
C     data as contained in NISH array (the closed shells) and the open
C     shell occupation data in _quadratic_ UDV matrix.
C
#include <implicit.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0)
#include <inforb.h>
      DIMENSION UDV(NASHT,NASHT)
      DIMENSION CMO(*)
      DIMENSION DCAO(NBAST,NBAST)
      DIMENSION DVAO(NBAST,NBAST)
      DIMENSION WRK(*)
C
      CALL QENTER('GTDMSO')
C
C Input : MO coefficents in CMO and active density matrix in MO basis
C Output: inactive and active density matrix AO (actually SO) basis
C
C
      CALL DZERO(DCAO,N2BASX)
      IF (NASHT.GT.0) CALL DZERO(DVAO,N2BASX)
C
      DO ISYM = 1, NSYM
         ICMOI = ICMO(ISYM) + 1
         IORBI = IORB(ISYM) + 1
         IASHI = IASH(ISYM) + 1
         IBASI = IBAS(ISYM) + 1
         NORBI = NORB(ISYM)
         NASHI = NASH(ISYM)
         NISHI = NISH(ISYM)
         NBASI = NBAS(ISYM)
C
C Inactive density
C
         IF (NISHI.GT.0) THEN
            CALL DGEMM(
     &        'N', 'T', NBASI, NBASI, NISHI,
     &        D2, CMO(ICMOI), NBASI,
     &            CMO(ICMOI), NBASI,
     &        D1, DCAO(IBASI,IBASI),NBAST
     &        )
         END IF
C
C Active density
C
        ICMOA = ICMOI + NBASI*NISHI
         IF (NASHI.GT.0) THEN
            CALL DGEMM(
     &         'N', 'N', NBASI, NASHI, NASHI,
     &         D1, CMO(ICMOA), NBASI,
     &             UDV(IASHI,IASHI), NASHT,
     &         D0, WRK,NBASI
     &      )
            CALL DGEMM(
     &         'N', 'T', NBASI, NBASI, NASHI,
     &         D1, WRK,       NBASI,
     &             CMO(ICMOA),NBASI,
     &         D1, DVAO(IBASI,IBASI),NBAST
     &      )
         END IF
      END DO
      CALL QEXIT('GTDMSO')
      END
C  /* Deck fckmat */
      SUBROUTINE FCKMAT(ONLYFC,DV,CMO,EMCMY,FC,FV,WRK,LFRSAV)
C
C 20-Jan-1988 Hans Joergen Aa. Jensen
C l.r. 930518 hjaaj: inserted F?AO to F? transformation
C
C     Driver for Fock matrix constructions.
C
C MOTECC-90: The algorithms used in this module, FCKMAT, are
C            described in Chapter 8 Appendix 8B of MOTECC-90
C            "Calculation of Generalized Fock Matrices from
C            Regular and One-Index Transformed Integrals"
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION DV(*), CMO(*), FC(*), FV(NNORBT,*), WRK(*)
C     ... FV(*,1) is the real FV, FV(*,2) is FS for DFT_SPINDNS
C         or FV(*,2) is the Kohn-Sham matrix for LSRHYBR
C         / hjaaj Feb 2010
      LOGICAL   ONLYFC, NOFV, LSRHYBR
C
C Used from common blocks:
C   INFORB : NASHT, NNBASX, NSYM, ...
C   INFPRI : P6FLAG()
C
#include <maxorb.h>
#include <infinp.h>
#include <inforb.h>
#include <infpri.h>
      PARAMETER( D1 = 1.0D0 )
#include <dftcom.h>
#include <dfterg.h>
C
      CALL QENTER('FCKMAT')
      NOFV   = ( (NASHT .EQ. 0) .OR. ONLYFC )
      LSRHYBR = SRHYBR .AND. .NOT. MCTYPE.EQ.0
C
      IF (DFT_SPINDNS .AND. LSRHYBR) THEN
         CALL QUIT('SRHYBR and non-singlet not implemented')
      END IF
C
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LFRSAV
      IF (.NOT. NOFV) THEN
         IF (DFT_SPINDNS) THEN
C        FCAO, FVAO, and FSAO must be consecutive
C        DCAO, DVAO, and DSAO must be consecutive
            ND = 3
            NF = 3
         ELSE
C        FCAO and FVAO must be consecutive
C        DCAO and DVAO must be consecutive
            ND = 2
            NF = 2
         END IF
         IF (LSRHYBR) THEN
C           we need an extra FSRAO and DSRAO for hybrid SRDFT-DFT
C           in order to evaluate both SR-DFT term and DFT term
            NF = NF + 1
         END IF
         CALL MEMGET('REAL',KFCAO,NF*N2BASX,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KDCAO,ND*N2BASX,WRK,KFREE,LFREE)
         KFVAO = KFCAO + N2BASX
         KFSAO = KFVAO + N2BASX
         KSRAO = KFCAO + (NF-1)*N2BASX
C        ... FSRAO(1) = WRK(KSRAO) is only for LSRHYBR /hjaaj
         KDVAO = KDCAO + N2BASX
         KDSAO = KDVAO + N2BASX
      ELSE
         ND = 1
         NF = 1
         CALL MEMGET('REAL',KFCAO,N2BASX,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KFVAO,0,WRK,KFREE,LFREE)
         KFSAO = KFVAO
         CALL MEMGET('REAL',KDCAO,N2BASX,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KDVAO,0,WRK,KFREE,LFREE)
         KDSAO = KDVAO
      END IF
C
C
C ***** Construct folded inactive and active parts of the one-electron
C       density matrix in AO-basis (DCAO and DVAO). DV is the active
C       part of the one-electron density matrix in MO-basis.
C
      CALL FCKDEN((NISHT.GT.0),.NOT.NOFV,
     *            WRK(KDCAO),WRK(KDVAO),CMO,DV,WRK(KFREE),LFREE)
C     CALL FCKDEN(GETDC,GETDV,DCAO,DVAO,CMO,DV,WRK,LWRK)
C
      IF (HSROHF) THEN
         CALL DAXPY(N2BASX,D1,WRK(KDVAO),1,WRK(KDCAO),1)
         CALL DSCAL(N2BASX,-D1,WRK(KDVAO),1)
      END IF
      CALL FCKMAO(NOFV,EMCMY,WRK(KFCAO),WRK(KFVAO),
     *            WRK(KDCAO),WRK(KDVAO),DV,CMO,WRK(KFREE),LFREE)
C
C ***** Transform inactive and active Fock-matrices to MO basis
C ***** (FC and FV) using DCAO as scratch space
C
      CALL DCOPY(N2BASX,WRK(KFCAO),1,WRK(KDCAO),1)
      CALL DGETSP(NBAST,WRK(KDCAO),WRK(KFCAO))
      IF (NSYM .GT. 1) CALL PKSYM1(WRK(KFCAO),WRK(KFCAO),NBAS,NSYM,2)
      CALL UTHUB(WRK(KFCAO),FC,CMO,WRK(KDCAO),NSYM,NBAS,NORB)
C
      IF (.NOT.NOFV) THEN
         CALL DCOPY(N2BASX,WRK(KFVAO),1,WRK(KDCAO),1)
         CALL DGETSP(NBAST,WRK(KDCAO),WRK(KFVAO))
         IF (NSYM .GT. 1) CALL PKSYM1(WRK(KFVAO),WRK(KFVAO),NBAS,NSYM,2)
         CALL UTHUB(WRK(KFVAO),FV,CMO,WRK(KDCAO),NSYM,NBAS,NORB)
         IF (DFT_SPINDNS) THEN
            CALL DCOPY(N2BASX,WRK(KFSAO),1,WRK(KDCAO),1)
            CALL DGETSP(NBAST,WRK(KDCAO),WRK(KFSAO))
            IF (NSYM .GT. 1)
     &         CALL PKSYM1(WRK(KFSAO),WRK(KFSAO),NBAS,NSYM,2)
            CALL UTHUB(WRK(KFSAO),FV(1,2),CMO,WRK(KDCAO),NSYM,NBAS,NORB)
         END IF
      END IF
      IF (LSRHYBR) THEN
         CALL DCOPY(N2BASX,WRK(KSRAO),1,WRK(KDCAO),1)
         CALL DGETSP(NBAST,WRK(KDCAO),WRK(KSRAO))
         IF (NSYM .GT. 1) CALL PKSYM1(WRK(KSRAO),WRK(KSRAO),NBAS,NSYM,2)
         CALL UTHUB(WRK(KSRAO),FV(1,2),CMO,WRK(KDCAO),
     &      NSYM,NBAS,NORB)
      END IF
C
      IF (P6FLAG(14)) THEN
        WRITE(LUPRI,1200)
        CALL OUTPKB(FC,NORB,NSYM,1,LUPRI)
        IF (.NOT.NOFV) THEN
          WRITE(LUPRI,1300)
          CALL OUTPKB(FV,NORB,NSYM,1,LUPRI)
        END IF
        IF (DFT_SPINDNS) THEN
          WRITE(LUPRI,1350)
          CALL OUTPKB(FV(1,2),NORB,NSYM,1,LUPRI)
        END IF
        IF (LSRHYBR) THEN
          WRITE(LUPRI,1400)
          CALL OUTPKB(FV(1,2),NORB,NSYM,1,LUPRI)
        END IF
      END IF
 1200 FORMAT(/' FC = Fock matrix, inactive part (MO-basis)')
 1300 FORMAT(/' FV = Fock matrix, active part (MO-basis)')
 1350 FORMAT(/' FS = Fock matrix, spin-density part (MO-basis)')
 1400 FORMAT(/' SRHYBR extra matrix (MO-basis)')
C
      CALL MEMREL('FCKMAT',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
C
C     In our srDFT local spin model, the DS is a function of DV;
C     thus we can handle this by adding the appropriate matrix function
C     of FS and DV to FC, creating an effective FC for the local spin
C     model.  / hjaaj Feb 2010
C
      IF (DFT_LOCALSPIN .AND. DFT_SPINDNS) THEN
#ifdef MOD_SRDFT
         CALL FCK_MAKE_FCEFF_LOCALSPIN(FC,FV(1,2),DV,WRK(KFREE),LFREE)
#else
         call quit('Sorry, srDFT local spin not in this version')
#endif
      END IF
C
C     Add srdft correction to EMCMY:
C
      EMCMY = EMCMY + ESRDFTY
C
      CALL QEXIT('FCKMAT')
      RETURN
      END
C  /* Deck fckmao */
      SUBROUTINE FCKMAO(ONLYFC,EMCMY,
     *                  FCAO,FVAO,DCAO,DVAO,DV,CMO,WRK,LFRSAV)
C
C Written 6-Nov-1983 Hans Jorgen Aa. Jensen
C Last revision 18-May-1993 hjaaj (moved AO to MO transf. outside)
C March 1997 - tsaue Screening
C DFT modifications T. Helgaker
C
C Purpose:
C  Construction of inactive and active parts of the
C  Fock matrix in AO-basis (FCAO and FVAO).
C
C  The Fock matrices are constructed as follows:
C  IF (DIRFCK) directly from AO integrals calculated on the fly,
C  else if (LUSUPM .eq. -1) from AO integral file,
C  else using atomic P super matrix elements.
C
C  If (onlyfc) only inactive Fock matrix FCAO is calculated.
C
C
#include <implicit.h>
      DIMENSION FCAO(*),DCAO(*), FVAO(*),DVAO(*),DV(*),CMO(*)
      DIMENSION WRK(LFRSAV)
      LOGICAL   ONLYFC
C
      PARAMETER (D0=0.0D0, HALF=0.5D0, D1 = 1.0D0)
C
C Used from common blocks:
C   INFINP : DIRFCK
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <r12int.h>
#include <infinp.h>
#include <inforb.h>
#include <infind.h>
#include <infdim.h>
#include <infpri.h>
#include <dfterg.h>
C
      LOGICAL NOFV
C
      CALL QENTER('FCKMAO')
C
#ifdef SIRGRAD_DEBUG
      P6FLAG(11) = .true. !Cdbg DEBUG
#endif
C
C
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LFRSAV
C
C ***** TWO-ELECTRON PART calculated in FCAO, FVAO
C       (FCAO, FVAO zeroed in FCK2AO)
C
      ESRDFTY = D0
C     ... make sure ESRDFTY is reset for CISRDFT following HFSRDFT
C         (ESRDFTY is in common block and is calculated inside FCK2AO)
      CALL FCK2AO(ONLYFC,FCAO,DCAO,WRK,KFREE,LFREE)
      IF(DODFT) THEN
         IF(NASHT.NE.0) CALL QUIT('Gradient not implemented for RO-DFT')
         EXCTRO = HALF*DDOT(N2BASX,DCAO,1,FCAO,1)
         CALL DFTKSMb(DCAO,FCAO,EDFTY,WRK(KFREE),LFREE,IPRFCK)
         EDFTY = EDFTY + EXCTRO - HALF*DDOT(N2BASX,DCAO,1,FCAO,1)
         IF (P6FLAG(11)) WRITE (LUPRI,'(A,F25.12)')
     &   ' Kohn-Sham DFT xc energy      :',EDFTY
      ENDIF
      IF (ESRDFTY .ne. D0 .and. P6FLAG(11))
     &   WRITE (LUPRI,'(A,F25.12)')
     &   ' srDFT energy correction      :',ESRDFTY
C
C     If this is a direct run, DCAO is returned in SO basis (to save
C     memory). Thus, we need to recreate it. K.Ruud, March-98
C
! hjaaj obsolete for srDFT, because DCAO is now never destroyed in SIRFCK
! when srDFT (instead a copy is made which is destroyed) /hjaaj Feb 2010
! NB! DV cannot be removed, needed for HSROHF and for direct MCSCF
      IF (DIRFCK .AND. .NOT.(DOHFSRDFT .OR. DOMCSRDFT)) THEN
         NOFV  = ( (NASHT .EQ. 0) .OR. ONLYFC )
         CALL FCKDEN((NISHT.GT.0),.NOT.NOFV,
     &               DCAO,DVAO,CMO,DV,WRK(KFREE),LFREE)
      END IF
      IF (HSROHF) THEN
         CALL DAXPY(N2BASX,D1,DVAO,1,DCAO,1)
         CALL DSCAL(N2BASX,-D1,DVAO,1)
      END IF
C
C ***** one-electron part in WRK(KH1AO)
C
      CALL MEMGET('REAL',KH1AO,N2BASX,WRK,KFREE,LFREE)
      CALL SIRH1(WRK(KH1AO),WRK(KFREE),LFREE)
      CALL MEMGET('REAL',KTMP,NNBASX,WRK,KFREE,LFREE)
      CALL PKSYM1(WRK(KTMP),WRK(KH1AO),NBAS,NSYM,-1)
      CALL DSPTGE(NBAST,WRK(KTMP),WRK(KH1AO))
      CALL MEMREL('FCKMAO.H1AO',WRK,KFRSAV,KTMP,KFREE,LFREE)
C
C ***** INACTIVE ENERGY : EMCMY; ESRDFTY is added in FCKMAT
C
      IF (NISHT .GT. 0) THEN
         EMCMY = DDOT(N2BASX,DCAO,1,WRK(KH1AO),1)
         IF (P6FLAG(11)) WRITE (LUPRI,'(A,F25.12)')
     &   ' Inactive one-electron energy :',EMCMY
         EMCMY = EMCMY + HALF*DDOT(N2BASX,DCAO,1,FCAO,1)+EDFTY
      ELSE
         EMCMY = 0.0D0
      END IF
      IF (P6FLAG(11)) WRITE (LUPRI,'(A,F25.12)')
     &   ' Inactive total energy        :',EMCMY
C
C ***** add one-electron contribution to FCAO
C
C     One-electron contribution is not added for the construction of 
C     the exchange operator in the auxilliary basis (WK/UniKA/04-11-2002).
      IF (.NOT. (R12TRA .OR. R12ECO))
     &     CALL DAXPY(N2BASX,D1,WRK(KH1AO),1,FCAO,1)
C
C     ***** Test output of Fock matrices *****
C
      IF (P6FLAG(18)) THEN
        IF (.NOT. (R12TRA .OR. R12ECO)) THEN
          WRITE(LUPRI,1000)
        ELSE
          WRITE(LUPRI,1001)
        END IF
        CALL OUTPUT(FCAO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
        IF (.NOT.ONLYFC) THEN
          WRITE(LUPRI,1100)
          CALL OUTPUT(FVAO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
        END IF
      END IF
 1000 FORMAT(/' FCAO = Fock matrix, inactive part (AO-basis)')
 1001 FORMAT(/' FCAO = Fock matrix, inactive part (auxilliary basis)')
 1100 FORMAT(/' FVAO = Fock matrix, active part (AO-basis)')
C
C *** End of subroutine FCKMAO
C
      CALL QEXIT('FCKMAO')
      RETURN
      END
C  /* Deck fck2ao */
      SUBROUTINE FCK2AO(ONLYFC,FCAO,DCAO,WRK,KFREE,LFREE)
C
C 2-Apr-1997 Hans Joergen Aa. Jensen
C     : 2-electron part extracted from FCKMAO
C DFT modifications T. Helgaker
C
C Purpose:
C  Construction 2-electron part of inactive and active parts of the
C  Fock matrix in AO-basis (FCAO and FVAO).
C
C  The Fock matrices are constructed as follows:
C  IF (DIRFCK) directly from AO integrals calculated on the fly,
C  else if (LUSUPM .eq. -1) from AO integral file,
C  else using atomic P super matrix elements.
C
C  If (onlyfc) only inactive Fock matrix FCAO is calculated.
C
C
#include <implicit.h>
      DIMENSION FCAO(*),DCAO(*), WRK(*)
      LOGICAL   ONLYFC
C
C Used from common blocks:
C   INFINP : DIRFCK
C   INFORB : NISHT
C   INFTAP : LUSUPM
C
#include <maxash.h>
#include <maxorb.h>
#include <infinp.h>
#include <inforb.h>
#include <infind.h>
#include <inftap.h>
#include <infdim.h>
#include <infpri.h>
#include <dftcom.h>
#include <gnrinf.h>
#include <priunit.h>
C
      DIMENSION ISYMDM(3), IFCTYP(3)    !JT Aug 09
C
      CALL QENTER('FCK2AO')
      KFRSAV = KFREE
C
C
      ISYMDM(1) = 1
      ISYMDM(2) = 1
      IF (HFXFAC .eq. 0.0D0) THEN
C        ... only Coulomb term
         IFCTYP(1) = 11
         IFCTYP(2) = 11
      ELSE
C        ... Coulomb and Exchange
         IFCTYP(1) = 13
         IFCTYP(2) = 13
      END IF
C
C JT Aug 09 beg
      IF (DFT_SPINDNS) THEN
         ISYMDM(3) = 1
         IFCTYP(3) = -12
      END IF
C JT Aug 09 end
C
      IF (ONLYFC) THEN
         NDMAT = 1
         IF (NISHT .GT. 0) THEN
            NFMAT = 1
         ELSE
            NFMAT = 0
         END IF
         IOFFAO = 1
      ELSE IF (NISHT .GT. 0 .OR. ADDSRI) THEN
         IF (DFT_SPINDNS) THEN !JT Aug 09
            NDMAT = 3
            NFMAT = 3
         ELSE
            NDMAT = 2
            NFMAT = 2
         ENDIF
         IOFFAO = 1
         IF (NISHT.EQ.0) CALL DZERO(DCAO,N2BASX)
C        ... zero DCAO for ADDSRI with NISHT.eq.0
      ELSE IF (NASHT.EQ.1 .AND. NACTEL.EQ.1) THEN 
C     one-electron system
         NDMAT = 2
         IF (DODFT .OR. DOHFSRDFT .OR. DOMCSRDFT) THEN
            NFMAT = 2
         ELSE
            NFMAT = 0
         END IF
         IOFFAO = 1
      ELSE
C        DCAO zero matrix, only FVAO contribution
C        (except for RO-DFT and srDFT we need FCAO for VxcAO)
         NDMAT = 2
         IF (DODFT .OR. DOHFSRDFT .OR. DOMCSRDFT) THEN
            NFMAT  = 2
            IOFFAO = 1
         ELSE
            NFMAT  = 1
            IOFFAO = 1 + N2BASX
C           ... off-set to FVAO,DVAO in FCAO,DCAO
         END IF
      END IF
      IF (HSROHF .OR. (DODFT .AND. NASHT .EQ. 1)) THEN
C
C This gives the exchange part of the active Fock matrix;
C we need only the difference between Fa and Q which is
C the negative of the exchange part
C
         IFCTYP(2)=12
      END IF
      CALL DZERO(FCAO,NDMAT*N2BASX)
C     ... note: NDMAT = # of FCAO matrices; NFMAT = # of non-zero FCAO matrices
      IF (NFMAT .GT. 0) THEN
      IF (LUSUPM .NE. -1 .AND.
     &   .NOT. (HSROHF .OR. DODFT .OR.
     &          DOCISRDFT .OR. DOHFSRDFT .OR. DOMCSRDFT)) THEN
         CALL RDSUPM(NFMAT,FCAO(IOFFAO),DCAO(IOFFAO),WRK(KFREE),LFREE)
      ELSE
         CALL SIRFCK(FCAO(IOFFAO),DCAO(IOFFAO),NFMAT,ISYMDM,IFCTYP,
     *               DIRFCK,WRK(KFREE),LFREE)
C        ... if (DIRFCK) construct Fock matrix directly
C            else construct Fock matrix from LUINTA
C            ISYMDM(1:NFMAT) = 1 : symmetry of density matrix (DUNP)
C            IFCTYP(1:NFMAT) = 13 : singlet, symmetric Fock matrix
C            NFMAT is max = 2 here
C
      END IF
      END IF
C
C *** End of subroutine FCK2AO
C
      CALL QEXIT('FCK2AO')
      RETURN
      END
C  /* Deck geth2 */
      SUBROUTINE GETH2(DV,PV,H2AC,FQ,ORDIAG,FULODG,WRK,LFRSAV)
C
C Fall '83, Hans Joergen Aa. Jensen
C Rewritten Oct 89 hjaaj
C Modified  11-Nov-89 for abacus.dertra, hjaaj
C
C Purpose:
C
C    To get array of 2-electron integrals with active indices
C    for CI-use.
C    To calculate FQ Fock matrix contribution
C    (FQ(p,t) = SUM(uvx): PV(tu:vx)*(pu:vx) )
C    To add non-Fock matrix contributions to orbital diagonal.
C
C Input:  DV (active 1-dm)
C         PV (active 2-dm)
C         ORDIAG
C         FULODG true : full orbital diagonal to be calculated
C
C Output: H2AC (active 2-el integrals)
C         FQ   (Fock Q-matrix)
C         ORDIAG (only changed if FULODG true)
C
C Scratch: WRK
C
C ****************************************************************
#include <implicit.h>
      DIMENSION H2AC(NNASHX,NNASHX), FQ(NORBT,NASHT), ORDIAG(*)
      DIMENSION DV(NNASHX), PV(NNASHX,NNASHX), WRK(LFRSAV)
      LOGICAL   FULODG
C
      PARAMETER ( D2=2.0D0 )
C
C Used from common blocks:
C   INFORB : NASHT,NORBT,NNASHX,...
C   INFIND : IROW,ISMO,ICH,IOBTYP
C   INFVAR : JWOPSY
C   INFPRI : IPRSIR
C CBGETDIS : IADH2
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infind.h>
#include <infvar.h>
#include <infpri.h>
#include <cbgetdis.h>
C
C     Local variables
C
      DIMENSION NEEDMU(6)
C
C     Definition of orbital types
C
#include <orbtypdef.h>

      CALL QENTER('GETH2 ')
      IF (NASHT.LE.1) THEN
         CALL QUIT('PROGRAM ERROR, GETH2 called with NASHT .le. 1')
      END IF
      KFRSAV= 1
      KFREE = KFRSAV
      LFREE = LFRSAV
C     if full orbital diagonal then only secondary-secondary
C     distributions not needed (but we don't need all distribu-
C     tions with inactive orbitals -- code -1),
C     else only active-active distributions needed.
      IF (FULODG) THEN
         DO 10 I = 1,5
            NEEDMU(I) = -1
   10    CONTINUE
         NEEDMU(3) = 1
         NEEDMU(6) = 0
         NKLWOP = NOCCT*NORBT
      ELSE
         DO 20 I = 1,6
            NEEDMU(I) = 0
   20    CONTINUE
         NEEDMU(3) = 1
         NKLWOP = 0
      END IF
C
C *** Allocate work space
C
      CALL MEMGET('REAL',KH2CD,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KPVCD,N2ASHX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KLWOP,NKLWOP,WRK,KFREE,LFREE)
      IF (IPRSIR .GT. 60) THEN
         WRITE (LUPRI,*)
         WRITE (LUPRI,*) '>>> test output from GETH2 <<<'
         WRITE (LUPRI,*) 'FULODG =',FULODG
         WRITE (LUPRI,*) 'LFRSAV =',LFRSAV
         WRITE (LUPRI,*) 'NEEDMU =',NEEDMU
         WRITE (LUPRI,*) 'KH2CD  =',KH2CD
         WRITE (LUPRI,*) 'KPVCD  =',KPVCD
         WRITE (LUPRI,*) 'KLWOP  =',KLWOP
         WRITE (LUPRI,*) 'KFREE  =',KFREE
         WRITE (LUPRI,*) 'LFREE  =',LFREE
      END IF
C
C *** Initialize H2AC
C
      CALL DZERO(H2AC,NNASHX*NNASHX)
      CALL DZERO(FQ,NORBT*NASHT)
      IF (FULODG) CALL MAKE_KLWOP(WRK(KLWOP))
C
C ****************************************************************
C     Loop over active-active distributions
C
      IDIST = 0
  100 CALL NXTH2M(IC,ID,WRK(KH2CD),NEEDMU,WRK,KFREE,LFREE,IDIST)
      IF (IPRSIR .GT. 65) THEN
         WRITE (LUPRI,*)
         WRITE (LUPRI,*)
     &      '***GETH2 Next H2CD(a,b) = (ab|cd) distribution'
         WRITE (LUPRI,*) 'IDIST,IC,ID = ',IDIST,IC,ID
         WRITE (LUPRI,*) 'KFREE,LFREE = ',KFREE,LFREE
      END IF
      IF (IDIST .LT. 0) GO TO 800
C     ... if idist .lt. 0 then no more distributions
C
C        Find symmetry of integrals
C
         IABMAX = IDAMAX(N2ORBX,WRK(KH2CD),1)
         IA     = (IABMAX-1)/NORBT + 1
         IB     = IABMAX - (IA-1)*NORBT
         IABSYM = MULD2H(ISMO(IA),ISMO(IB))
C
         IF (IC.LT.ID) THEN
            ISWAP = IC
            IC = ID
            ID = ISWAP
         END IF
         ICSYM  = ISMO(IC)
         IDSYM  = ISMO(ID)
         ICDSYM = MULD2H(ICSYM,IDSYM)
         INTSYM = MULD2H(IABSYM,ICDSYM)
         ITYPC  = IOBTYP(IC)
         ITYPD  = IOBTYP(ID)
         ITYPCD = IDBTYP(ITYPC,ITYPD)
         IF (FULODG) THEN
            IF (ITYPCD.EQ.6) ITYPCD = 0
            IF (ITYPCD.EQ.1 .AND. IC.NE.ID) ITYPCD = 0
            IF (ITYPCD.EQ.4 .AND. ICDSYM.NE.JWOPSY) ITYPCD = 0
         ELSE
            IF (ITYPCD.NE.3) ITYPCD = 0
         END IF
         IF (IPRSIR .GT. 65) THEN
            WRITE (LUPRI,*)
            WRITE (LUPRI,*)'IABSYM,ICDSYM,INTSYM ',IABSYM,ICDSYM,INTSYM
            WRITE (LUPRI,*)'IDIST, IC,    ID     ',IDIST, IC,    ID
            WRITE (LUPRI,*)'ITYPC, ITYPD, ITYPCD ',ITYPC, ITYPD, ITYPCD
         END IF
      IF (ITYPCD.EQ.0) GO TO 100
         IF (IPRSIR .GT. 65) THEN
            WRITE (LUPRI,*)
            WRITE (LUPRI,*) 'H2CD matrix for this |cd) distribution :'
            CALL OUTPUT(WRK(KH2CD),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         END IF
         IF (ITYPCD .EQ. 3) THEN
C        ... |cd) active-active distribution
            NCW  = ICH(IC)
            NDW  = ICH(ID)
            NCDW = IROW(NCW) + NDW
            IF (IPRSIR .GT. 65) THEN
               WRITE (LUPRI,*)
               WRITE (LUPRI,*) 'active-active |cd) distribution'
               WRITE (LUPRI,*) 'NCW,   NDW,   NCDW   ',NCW,NDW,NCDW
            END IF
            CALL SIRPVD(NCW,NDW,WRK(KPVCD),PV,1)
            IF (NCW .NE. NDW) CALL DSCAL(N2ASHX,D2,WRK(KPVCD),1)
C           ... C,D and D,C contributions are identical
C
C           Add active-active elements to H2AC
C
            CALL ADH2AC(H2AC(1,NCDW),WRK(KH2CD),IABSYM)
C
C           Add the contribution to transposed FQ matrix from this
C           distribution
C
C           Note the active-active block of the FQ-matrix is also
C           calculated, it is used for
C           a) for gradient correction to orbital Hessian
C           b) if FLAG(23), that is if active-active orbital rotations
C           c) if FLAG(25), write converged Fock matrix to LUSIFC
C                           (information for geometry optimization)
C
            CALL ADDFQ(NCW,NDW,FQ,WRK(KH2CD),WRK(KPVCD),INTSYM)
C           CALL ADDFQ(NXW,NYW,FQ,H2XY,PVXY,JH2SYM)
C
         END IF
C
C        If exact (full) orbital diagonal requested,
C        then call ADDODG
C
         IF (FULODG) THEN
            CALL ADDODG(IC,ID,WRK(KLWOP),WRK(KH2CD),DV,PV,ORDIAG)
C           CALL ADDODG(IC0,ID0,KLWOP,H2CD,DV,PV,ORDIAG)
         END IF
C
C
C        Go get next active-active distributions
      GO TO 100
C
C     arrive at 800 when finished with all active-active distributions
C
  800 CONTINUE
C
C *** IADH2 .lt. 0 means that H2AC is in core.
C
      IADH2 = -1
C
C
C *** end of subroutine GETH2
C
      CALL MEMREL('GETH2 ',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('GETH2 ')
      RETURN
      END
C  /* Deck addodg */
      SUBROUTINE ADDODG(IC0,ID0,KLWOP,H2CD,DV,PV,ORDIAG)
C
C  Written by Hans Agren and Hans Joergen Aa. Jensen / March 1984
C  Last revision 26-Jul-1984 hjaaj (optimized ITYPCD=4 code)
C                 7-May-1984 ha/hjaaj
C                11-Oct-1984 ha/hjaaj (act-act orbdiag)
C  Revised Oct 89 hjaaj to only include ORDIAG contributions.
C
C Purpose:
C
C    To add 2-el. integral contributions to the diagonal
C    of the orbital hessian (except active-active rotations).
C
C Input:
C       IC0,ID0 : two fixed indices of (**|cd) integrals in H2CD
C       H2CD : (**|cd) integrals, where c is IC0 and d is ID0
C       DV : the 1-el. active density matrix.
C       PV : the 2-el. active density matrix.
C Output:
C       ORDIAG : Add non-Fock matrix contributions to
C                the hessian orbital-part diagonal
C                (except active-active rotations)
C Scratch:
C
C    The subroutine is divided into 5 main blocks,
C    branched according to the type of first pair of transformed orbital
C    index (CD) that is encountered:
C
C    Block   CD-distr. AB-distr.       Contributing to:
C
C     100      (ii/     /aa)             Ordiag(i,a)
C
C     200      (xi/     /vi)             Ordiag(i,u)
C              (ui/     /vi)             Ordiag(i,u)
C
C     300      (xy/     /ii)             Ordiag(i,u)
C                       /aa)             Ordiag(u,a)
C
C     400      (ai/     /ai)             Ordiag(i,a)
C
C     500      (av/     /ax)             Ordiag(u,a)
C
C   Symmetry of Ordiag is JWOPSY.
C
C   Reference: Appendix A in H.J.Aa.Jensen and H.Agren,
C              Chem.Phys. 104 (1986) 229-250
C
C
C ****************************************************************
C
#include <implicit.h>
      DIMENSION KLWOP(NOCCT,NORBT)
      DIMENSION DV(*),PV(NNASHX,*),H2CD(NORBT,*), ORDIAG(*)
C
      PARAMETER (D1=1.D0, D2 = 2.D0, D4 =4.D0, D12 = 12.D0)
C
C Used from common blocks:
C   INFORB : NNASHX,NORBT,...
C   INFIND : JTACT,JTSEC,ISW(*),IOBTYP(*),...
C   INFVAR : NWOPT,JWOP(2,*),JWOPSY
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infind.h>
#include <infvar.h>
C
#include <orbtypdef.h>
C
C
C ****************************************************************
C
C  ITYPCD-values:  1=i*i :  2=t*i : 3=t*t : 4=a*i : 5=a*t : 6=a*a
C
      IC     = IC0
      ID     = ID0
      ITYPC  = IOBTYP(IC)
      ITYPD  = IOBTYP(ID)
      ITYPCD = IDBTYP(ITYPC,ITYPD)
C
      GO TO (100,200,300,400,500,600) ITYPCD
C
C*******************************************************************
C
C   **SECTION 1** (cd/ inactive-inactive distribution
C
C    Process  (cd/ab)=(ii/aa) integrals.
C    ORDIAG(i,a) = ORDIAG(i,a) - 4 (ii/aa)     eq. (A.2)
C    (contributes to ORDIAG(i,a) when symmetry of (ia) is JWOPSY
C     and when IC = ID)
C
  100 CONTINUE
      IF (IC .NE. ID) GO TO 9999
      ICW   = ISW(IC)
      ISYMC = ISMO(IC)
      ISYMA = MULD2H(JWOPSY,ISYMC)
      NSSHA = NSSH(ISYMA)
      IF (NSSHA .GT. 0) THEN
         IOFFA = IORB(ISYMA) + NOCC(ISYMA)
         DO 105 IA = IOFFA + 1, IOFFA + NSSHA
            INDXO = KLWOP(ICW,ISW(IA)-NISHT)
            IF (INDXO.GT.0) THEN
               ORDIAG(INDXO) = ORDIAG(INDXO) - D4*H2CD(IA,IA)
            END IF
  105    CONTINUE
      END IF
      GO TO 9999
C
C*******************************************************************
C
C   **SECTION 2** (cd/ active-inactive distribution
C
C    Process  (cd/ab)=(ui/vi) and (cd/ab)=(xi/vi) integrals.
C    Ordiag(i,u) = Ordiag(i,u) + 4 sum(v) P(ux,uv)(vi/xi)
C                + 12 sum(v) (DELuv - D(u,v)) (vi/ui)
C    -- Eq. (A.1)
C    (contributes to ORDIAG only for sym(c) .eq. sym(a) and for
C     b = d)
C
  200 CONTINUE
C
C     Swap IC and ID ?
C
      IF (ITYPC.NE.JTACT) THEN
         ISWAP = IC
         IC    = ID
         ID    = ISWAP
      END IF
C
      ISYMC  = ISMO(IC)
      ISYMD  = ISMO(ID)
      ISYMCD = MULD2H(ISYMC,ISYMD)
      ICW    = ISW(IC)
      IDW    = ISW(ID)
      NCW    = ICW - NISHT
      NASHC  = NASH(ISYMC)
      IASHC  = IASH(ISYMC)
      NASHD  = NASH(ISYMD)
      IASHD  = IASH(ISYMD)
C     ... note that IBW = IDW, IDW used below
      DO 205 IAW = NISHT+1,ICW
         IA  = ISX(IAW)
         NAW = IAW - NISHT
C
C  First add (vi/xi) = (cd/ad) contribution to Ordiag(i,u)
C
         DO 207 NUW = IASHD+1,IASHD+NASHD
            INDXO = KLWOP(IDW,NUW)
         IF (INDXO.EQ.0) GO TO 207
            IF (NUW .GT. NAW) THEN
               INDPUX = IROW(NUW) + NAW
            ELSE
               INDPUX = IROW(NAW) + NUW
            END IF
            IF (IAW .EQ. ICW) THEN
               TERM = PV(INDPUX,INDPUX)
            ELSE
               IF (NUW .GT. NCW) THEN
                  INDPUV = IROW(NUW) + NCW
               ELSE
                  INDPUV = IROW(NCW) + NUW
               END IF
               TERM = D2 * PV(INDPUV,INDPUX)
            END IF
            ORDIAG(INDXO) = ORDIAG(INDXO) + TERM*D4*H2CD(IA,ID)
  207    CONTINUE
C
C     Then add (vi/ui) (and (ui/vi))
C
      IF (ISYMCD .NE. JWOPSY) GO TO 205
         INDXO = KLWOP(IDW,NCW)
         IF (IAW .EQ. ICW) THEN
            IF (INDXO.NE.0) ORDIAG(INDXO) = ORDIAG(INDXO) +
     *         D12*(D1 - DV(IROW(NCW+1))) * H2CD(IA,ID)
         ELSE
            TERM = D12*DV(IROW(NCW) + NAW) * H2CD(IA,ID)
            IF (INDXO.NE.0) ORDIAG(INDXO) = ORDIAG(INDXO) - TERM
            INDXO = KLWOP(IDW,NAW)
            IF (INDXO.NE.0) ORDIAG(INDXO) = ORDIAG(INDXO) - TERM
         END IF
C
C
  205 CONTINUE
      GO TO 9999
C
C*******************************************************************
C
C   **SECTION 3** (cd/ active-active distribution
C
C    Process  (cd/ab)=(vx/ii),(uv/ii), and (vx/aa) integrals.
C    Ordiag(i,u) = Ordiag(i,u) + 2 P(uu,vx)(vx/ii)
C                - 4 (DELuv - D(u,v)) (uv/ii)
C    -- Eq. (A.2)
C    Ordiag(u,a) = Ordiag(u,a) + 2 P(uu,vx)(vx/aa)
C    -- Eq. (A.4)
C    (contribute to ORDIAG only for sym(c) .eq. sym(d) and for
C     a = b, sym(ca) = jwopsy)
C
  300 CONTINUE
C
C     Swap IC and ID ?
C
      IF (IC.LT.ID) THEN
         ISWAP = IC
         IC    = ID
         ID    = ISWAP
      END IF
C
      ISYMC  = ISMO(IC)
      ISYMD  = ISMO(ID)
      ISYMCD = MULD2H(ISYMC,ISYMD)
      IF (ISYMCD .NE. 1) GO TO 9999
      ICW   = ISW(IC)
      IDW   = ISW(ID)
      NCW   = ICW - NISHT
      NDW   = IDW - NISHT
      NASHC = NASH(ISYMC)
      IASHC = IASH(ISYMC)
      NASHD = NASH(ISYMD)
      IASHD = IASH(ISYMD)
C
      INDPVX = IROW(NCW) + NDW
      DO 305 IAW = 1,NORBT
      IF (IAW .GT. NISHT .AND. IAW .LE. NOCCT) GO TO 305
C     ... only /ii) and /aa)
         IA    = ISX(IAW)
         ISYMA = ISMO(IA)
         NASHA = NASH(ISYMA)
         NSSHA = NSSH(ISYMA)
C
         FAC = D2*H2CD(IA,IA)
         IF (NCW .NE. NDW) FAC = D2 * FAC
         DO 307 NUW=IASH(ISYMA)+1,IASH(ISYMA)+NASHA
            IF (IAW .GT. NOCCT) THEN
               INDXO = KLWOP(NISHT+NUW,(IAW-NISHT))
            ELSE
               INDXO = KLWOP(IAW,NUW)
            END IF
            IF (INDXO.NE.0) THEN
               INDPUU = IROW(NUW+1)
               ORDIAG(INDXO) = ORDIAG(INDXO) + FAC*PV(INDPUU,INDPVX)
            END IF
  307    CONTINUE
C
C      Add (uv/ii) contributions.
C
       IF(IAW.GT.NISHT) GO TO 305
          ISYMAC = MULD2H(ISYMA,ISYMC)
       IF(ISYMAC .NE. JWOPSY) GO TO 305
          INDUV = IROW(NCW) + NDW
          INDXO = KLWOP(IAW,NCW)
          IF (NCW.GT.NDW) THEN
             TERM = D4*DV(INDUV) * H2CD(IA,IA)
             IF (INDXO.NE.0) ORDIAG(INDXO) = ORDIAG(INDXO) + TERM
             INDXO = KLWOP(IAW,NDW)
             IF (INDXO.NE.0) ORDIAG(INDXO) = ORDIAG(INDXO) + TERM
          ELSE IF (INDXO.NE.0) THEN
             ORDIAG(INDXO) = ORDIAG(INDXO) +
     *          D4*(DV(INDUV) - D1) * H2CD(IA,IA)
          END IF
C
  305 CONTINUE
      GO TO 9999
C
C*******************************************************************
C
C   **SECTION 4** (cd/ secondary-inactive distribution
C
C    Process  (cd/ab)=(ai/ai) integrals.
C    Ordiag(i,a) = Ordiag(i,a) + 12 (ai/ai)    eq. (A.2)
C    (contributes to ORDIAG(i,a) when symmetry of (ia) is JWOPSY)
C
  400 CONTINUE
C
C     Swap IC and ID ?
C
      IF (ITYPC.NE.JTSEC) THEN
         ISWAP = IC
         IC    = ID
         ID    = ISWAP
      END IF
C
      ISYMC  = ISMO(IC)
      ISYMD  = ISMO(ID)
      ISYMCD = MULD2H(ISYMC,ISYMD)
      IF (ISYMCD .NE. JWOPSY) GO TO 9999
C
      ICW   = ISW(IC)
      IDW   = ISW(ID)
      INDXO = KLWOP(IDW,ICW-NISHT)
      IF (INDXO.GT.0) THEN
         ORDIAG(INDXO) = ORDIAG(INDXO) + D12*H2CD(IC,ID)
      END IF
      GO TO 9999
C
C*******************************************************************
C
C   **SECTION 5** (cd/ active-secondary distribution
C
C    Process  (cd/ab)=(ax/av) integrals.
C    Ordiag(u,a) = Ordiag(u,a) + 4 sum(v) P(ux,uv)(ax/av)
C    -- Eq. (A.4)
C    (contributes to ORDIAG only for symmetry equality of x and v)
C
  500 CONTINUE
C
C     Swap IC and ID ?
C
      IF (ITYPC.NE.JTSEC) THEN
         ISWAP = IC
         IC    = ID
         ID    = ISWAP
      END IF
C
      ISYMC  = ISMO(IC)
      ISYMD  = ISMO(ID)
      ISYMCD = MULD2H(ISYMC,ISYMD)
      ICW   = ISW(IC)
      IDW   = ISW(ID)
      NCW   = ICW - NISHT
      NDW   = IDW - NISHT
      NASHC = NASH(ISYMC)
      IASHC = IASH(ISYMC)
      NASHD = NASH(ISYMD)
      IASHD = IASH(ISYMD)
C     ... note that IAW = ICW = 'a', ICW used below
      DO 505 IBW = NISHT+1,IDW
        IB   = ISX(IBW)
        NBW  = IBW - NISHT
        FAC  = D4*H2CD(IB,IC)
        DO 507 NUW=IASHC+1,IASHC+NASHC
          INDXO = KLWOP((NISHT+NUW),NCW)
        IF (INDXO.EQ.0) GO TO 507
          IF(NUW .GT. NBW) THEN
             INDPUX = IROW(NUW) + NBW
          ELSE
             INDPUX = IROW(NBW) + NUW
          END IF
          IF(NDW .EQ. NBW) THEN
             TERM = PV(INDPUX,INDPUX)
          ELSE
             IF(NUW .GT. NDW) THEN
                INDPUV = IROW(NUW) + NDW
             ELSE
                INDPUV = IROW(NDW) + NUW
             END IF
             TERM = D2*PV(INDPUV,INDPUX)
          END IF
          ORDIAG(INDXO) = ORDIAG(INDXO) + TERM*FAC
  507   CONTINUE
  505 CONTINUE
      GO TO 9999
C
C*******************************************************************
C
C   **SECTION 6** (cd/ secondary-secondary distribution
C
C   -- not needed --
C
  600 CONTINUE
      GO TO 9999
C
C*******************************************************************
C
C End of subroutine ADDODG
C
 9999 RETURN
      END
C  /* Deck grad */
      SUBROUTINE GRAD (CREF,CMO,INDXCI,DV,PV,FC,FV,FQ,FCAC,H2AC,
     *                 G,EMCMY,EMCACT,WRK,LFREE)
C
C Revisions:
C   15-May-1984 hjaaj
C    9-Nov-1984 hjaaj (average GORB and ORDIAG if requested)
C
C Purpose:
C   This subroutine computes the energy and the gradient.
C   The diagonal of the CI-matrix is calculated and stored on LUIT2
C   in CIDIAG.
C   The matrices DV, PV, H2AC, FC, FV, FQ, FCAC are also returned.
C
C MOTECC-90: The purpose of this module, GRAD, and the algorithms
C            used  are described in Chapter 8 Section C.1 of MOTECC-90
C            "Orbital Classes, Gradient and NEO Transformations"
C
C Input:
C   CREF, the CI vector for the current expansion point (CEP).
C   CMO, the orbitals for the CEP
C
C Output:
C   EMCACT, energy of the CEP.
C   G, gradient of the CEP.
C   EMCMY, inactive energy part of EMCSCF.
C   DV, one-electron charge (and spin if srDFT) d.m. for CREF:
C   PV, two-electron d.m. for CREF.
C   H2AC, the two-electron integrals over active orbitals.
C   FC, the inactive Fock matrix.
C   FCAC, part of FC with active orbitals indices.
C   FV, the active Fock matrix.
C   FQ, the "Q" Fock matrix
C
C Scratch:
C   WRK
C
#include <implicit.h>
      DIMENSION CREF(*),CMO(*),DV(NNASHX,*),PV(*),FC(*),FV(*),FQ(*)
      DIMENSION FCAC(*),H2AC(*)
      DIMENSION INDXCI(*), G(*), WRK(LFREE)
#include <iratdef.h>
C
C -- local constants
C
      PARAMETER (D0 = 0.0D0, D1=1.0D0)
      PARAMETER (THRSML = 1.D-8)
#include <thrzer.h>
C
C Used from common blocks:
C   INFINP : LSYM,SUPSYM,?
C   INFVAR : MAXOCC,NCONF,NWOPT,NVAR,JWOP(2,*)
C   INFDIM : ?
C   INFOPT : GNORM(3),?
C   -- GNORM, norm of the CI, the orbital, and the total gradient.
C   INFTAP : LUIT2,?
C   INFTIM : NCALLS,TIMCPU,TIMWAL;   IDTIM above is index for these
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <infvar.h>
#include <infind.h>
#include <inforb.h>
#include <infdim.h>
#include <infopt.h>
#include <inftap.h>
#include <infpri.h>
#include <dftcom.h>
C ---
      PARAMETER (IDTIM = 1)
#include <inftim.h>
C
C -- local variables
C
      LOGICAL FULODG, AVEODG
      CHARACTER*8 TABLE(3)
C
C -- DATA statements
C
      DATA TABLE/'********','ORBDIAG ','EODATA  '/
C
      CALL QENTER('GRAD  ')

C
C     If spin-dependent srDFT then both charge and spin density 1-el matrices
      IF (DFT_SPINDNS) THEN
         NDV = 2
      ELSE
         NDV = 1
      END IF
#ifdef SIRGRAD_DEBUG
CFRAN 18/02/2010
      WRITE(LUPRI,*) ' DFT_SPINDNS,NDV= ', DFT_SPINDNS,NDV
CFRAN 18/02/2010 END
#endif
C
C ******************************************************************
C Step 0: WRK space allocation:
C
C .. 1.0 .. (all subroutines except MAKDM, FCKMAT)
C .. 1.2 .. (FCKMAT)
      KWRKF  = 1
      LWRKF  = LFREE - KWRKF
C .. 1.3 .. (GETH2,ORDIAG and ORBGRD)
      KORDIA = 1
      KW3    = KORDIA + NWOPT
      LW3    = LFREE  - KW3
C .. 1.3.1 .. (ORBGRD)
      KUDV   = KW3
      KLAST  = MAX(KWRKF,KUDV + NDV*N2ASHX)
C .. 1.4 .. (CIGRAD)
      KWRK4  = 1
      LWRK4  = LFREE- KWRK4
C Test if WRK is big enough:
      IF (KLAST.GT.LFREE) CALL ERRWRK('GRAD',-KLAST,LFREE)
C
C ******************************************************************
C
      FULODG = FLAG(12)
      AVEODG = .NOT.FULODG
C
C ******************************************************************
C Step 1: calculate one- and two-electron density matrices for
C         active orbital indices.
C
C
#ifdef SIRGRAD_DEBUG
CFRAN 18/02/2010
      WRITE(LUPRI,*) ' NDV= ', NDV
CFRAN 18/02/2010 END
#endif
      CALL GETTIM (T0,W0)
      IF (NWOPT.EQ.0) FULODG = .FALSE.
      IF (NASHT.EQ.0) GO TO 1999
      IF (NASHT .EQ. 1 .AND. NACTEL .LT. 2) THEN
         DV(1,1) = NACTEL
         PV(1) = D0
         IF (NDV.GE.2) DV(1,2) = NACTEL
      ELSE IF (HSROHF) THEN
         CALL DZERO(DV,NDV*NNASHX)
         DO I=1,NASHT
            II=I*(I+1)/2
            DV(II,1)=D1
            IF (NDV.GE.2) DV(II,2)=D1
         END DO
         CALL DZERO(PV,NNASHX*NNASHX)
      ELSE
         CALL MAKDM(CREF,DV,PV,INDXCI,WRK,LFREE)
C        CALL MAKDM(CREF,DV,PV,INDXCI,WRK,LFREE)
         IF (NDV.GE.2) THEN
            IF (DFT_LOCALSPIN) THEN
#ifdef MOD_SRDFT
               CALL MAKE_DS_LOCALSPIN(DV(1,1),DV(1,2),WRK(KWRKF),LWRKF)
#else
               CALL QUIT('Sorry! Local spin not in this version')
#endif
            ELSE
               CALL MAKDSV(CREF,DV(1,2),INDXCI,WRK,LFREE)
            END IF
         END IF
      END IF
C
#ifdef SIRGRAD_DEBUG
CFRAN 17/02/2010
      WRITE(LUPRI,*) ' I SET P6FLAG(15) = .TRUE. !'
      P6FLAG(15) = .TRUE.
CFRAN 17/02/2010 END
#endif
      IF (P6FLAG(15)) THEN
         WRITE (LUPRI,'(/A)')
     &   ' GRAD: DV = One-el. density matrix, active part, MO-basis'
         CALL OUTPAK(DV,NASHT,1,LUPRI)
         IF (NDV.GE.2) THEN
            WRITE (LUPRI,'(/A)')
     &      ' GRAD: DSV = One-el. spin den. mat., active part, MO-basis'
            CALL OUTPAK(DV(1,2),NASHT,1,LUPRI)
         END IF
      END IF
      IF (P6FLAG(25)) THEN
         WRITE (LUPRI,'(/A)')
     &   ' PV = Two-el. density matrix, active part, MO-basis'
         CALL OUTPUT(PV,1,NNASHX,1,NNASHX,NNASHX,NNASHX,-1,LUPRI)
      END IF
C
 1999 CONTINUE
      CALL GETTIM(T1,W1)
C
C ******************************************************************
C Step 2: calculate inactive and active Fock matrices.
C
      CALL FCKMAT((NASHT.EQ.0),DV,CMO,EMCMY,FC,FV,WRK(KWRKF),LWRKF)
C     CALL FCKMAT(ONLYFC,DV,CMO,EMCMY,FC,FV,WRK,LFREE)
      IF (NASHT .EQ. 0) CALL DZERO(FV,NNORBT)
C     ... ORBGRD and ORBSIG always reference FV
C
      CALL GETTIM(T2,W2)
C
C ******************************************************************
C Step 3: calculate "FQ" matrix and two-electron integrals with
C         active orbital indices.
C         If fulodg also add the 2-electron non-fock matrix
C         contributions to the diagonal of orbital Hessian.
C
C
      IF (NWOPT.GT.0) CALL DZERO(WRK(KORDIA),NWOPT)
      IF (NASHT .GE. 1) THEN
         CALL GETAC(FC,FCAC)
         IF (P6FLAG(21)) THEN
            WRITE(LUPRI,'(/A)')
     &      ' Active block of inactive Fock matrix FCAC'
            CALL OUTPAK (FCAC,NASHT,1,LUPRI)
         END IF
      END IF
      IF (NASHT .GT. 1 .AND. .NOT. HSROHF) THEN
         CALL GETH2(DV,PV,H2AC,FQ,WRK(KORDIA),FULODG,
     &              WRK(KW3),LW3)
C        CALL GETH2(DV,PV,H2AC,FQ,ORDIAG,FULODG,WRK,LFRSAV)
         IF (P6FLAG(23) .AND. IPRI6.GT.60 .AND. FULODG) THEN
            WRITE (LUPRI,'(4(/A))')
     &         ' Non-Fock type contributions to inact-act and act-sec',
     &         ' ----------------------------------------------------',
     &         ' diagonal of orbital Hessian',
     &         '----------------------------'
            CALL PRKAP (NWOPT,WRK(KORDIA),-0.1D0,LUPRI)
         END IF
C
         IF (P6FLAG(26)) THEN
            WRITE(LUPRI,'(/A)') ' Two-electron integrals H2AC'
            CALL OUTPUT(H2AC,1,NNASHX,1,NNASHX,NNASHX,NNASHX,-1,LUPRI)
         END IF
C
C        *** calculate active energy, before we (maybe) modify FCAC
C
         EMCACT = ENRACT(DV,FCAC,FQ)
C
      ELSE IF (HSROHF) THEN
         IF (ADDSRI) THEN
            call quit('Something needs to be done for srdft and hsrohf')
         END IF
         CALL DZERO(FQ,NORBT*NASHT)
         CALL GETAC(FV,FCAC)
         EMCACT = -ENRACT(DV,FCAC,FQ)/2
      ELSE IF (NASHT .EQ. 1) THEN
         CALL DZERO(FQ,NORBT)
         EMCACT = ENRACT(DV,FCAC,FQ)
      ELSE
         EMCACT = D0
      END IF
      CALL GETTIM (T3,W3)
C
C
C ******************************************************************
C Step 3: calculate Fock matrix contributions to the diagonal of the
C         orbital Hessian.
C
      IF (NWOPT .GT. 0) THEN
         IF (HSROHF) CALL DAXPY(NNORBT,-D1,FV,1,FC,1)
         CALL ORDIAG(DV,PV,FC,FV,FQ,H2AC,WRK(KORDIA),
     &               AVEODG,IPRSIR)
C        CALL ORDIAG(DV,PV,FC,FV,FQ,H2AC,DIAG,AVEODG,IPRINT)
         IF (P6FLAG(23)) THEN
            IF (FULODG) THEN
               WRITE (LUPRI,'(/A/A)')
     &         ' Diagonal of orbital Hessian',
     &         ' ---------------------------'
            ELSE
               WRITE (LUPRI,'(/A/A)')
     &         ' Approximate diagonal of orbital Hessian',
     &         ' ---------------------------------------'
            END IF
            IF (SUPSYM) WRITE (LUPRI,'(/A)')
     &         ' (Diagonal is microcanonically averaged.)'
            CALL PRKAP (NWOPT,WRK(KORDIA),-0.1D0,LUPRI)
         END IF
      END IF
      CALL GETTIM(T4,W4)
C
C
C ******************************************************************
C Step 5: put the pieces together to the orbital gradient.
C
C
#if defined (VAR_OLDREDCHCK)
 5100 CONTINUE
#endif
      IF (NWOPT.EQ.0) GO TO 5999
C     Unpack DV into UDV ( SP to SI format)
      IF (NASHT .GT. 0) THEN
         CALL DSPTSI(NASHT,DV,WRK(KUDV))
         IF (NDV.GE.2) CALL DSPTSI(NASHT,DV(1,2),WRK(KUDV+N2ASHX))
      END IF
C     Initialize gradient
      CALL DZERO(G(NCONF+1),NWOPH)
      CALL ORBGRD(NWOPH,WRK(KUDV),FC,FV,FQ,G(NCONF+1))
C     CALL ORBGRD(NWOPX,UDV,FC,FV,FQ,GORB)
C
C
C ******************************************************************
C
#if defined (VAR_OLDREDCHCK)
CHJ-870208: THESE WHICH WILL BE CAUGHT WITH THIS, WE USUALLY DO NOT
C           WANT TO REMOVE (MAY E.G. BE INVOLVED IN AVERAGING).
C
C ******************************************************************
C Step 5b: Unless FLAG(??) then check and remove some redundant
C          orbital rotations.
C
C  A redundant orbital rotation (because spanned by configurations)
C  is identified through zero gradient value and zero value on
C  Hessian diagonal.  The test is not performed in first macro
C  iteration because orbital rotations may be temporarily redundant
C  because of Hartree-Fock or similar "pathological" reference
C  points. This check will in general not catch redundant
C  active-active rotations.
C
      IF (.NOT.FLAG(??) .AND. ITMAC.GT.1) THEN
         FLAG(??) = .TRUE.
         IREDUN = 0
         DO 5900 IG = 1,NWOPT
            IF (ABS(G(NCONF+IG)).LT.THRZER .AND.
     &          ABS(WRK(KORDIA-1+IG)).LT.THRZER) THEN
               IREDUN = IREDUN + 1
               IF (P4FLAG(8)) THEN
                  IF (IREDUN.EQ.1) WRITE (LUW4,5911)
                  WRITE (LUW4,5921) IG,JWOP(1,IG),JWOP(2,IG)
               END IF
            END IF
 5900    CONTINUE
         IF (IREDUN.EQ.0) THEN
C           WRITE (LUW4,'(/A)')
C    &         ' No redundant orbital variations found.'
         ELSE
            DO 135 J = 1,MXCORB
               DO 135 I = 1,MAXOCC
                  KLWOP(I,J) = 0
  135       CONTINUE
            IGNEW = 0
            DO 5910 IG = 1,NWOPT
               IF (ABS(G(NCONF+IG)).GT.THRZER .OR.
     *             ABS(WRK(KORDIA-1+IG)).GT.THRZER) THEN
                  IGNEW = IGNEW + 1
                  K = JWOP(1,IG)
                  L = JWOP(2,IG)
                  JWOP(1,IGNEW) = K
                  JWOP(2,IGNEW) = L
                  KLWOP(ISW(K),ISW(L)-NISHT) = IG
                  WRK(KORDIA-1+IGNEW) = WRK(KORDIA-1+IG)
               ELSE
                  ISK = ISMO(JWOP(1,IG))
               END IF
 5910       CONTINUE
            WRITE (LUW4,5941) IREDUN,NWOPT,IGNEW
            NWOPT = IGNEW
            NVAR = NCONF + IGNEW
            JWOP(1,IGNEW+1) = 0
            JWOP(2,IGNEW+1) = 0
            GO TO 5100
         END IF
      END IF
 5911 FORMAT(/,' Redundant orbital variations removed:',
     *       /,T10,'Old variable no.',T30,'  r    s',
     *       /,T10,'----------------',T30,'---- ----')
 5921 FORMAT(T15,I5,T30,I3,I5)
 5941 FORMAT(/,1X,I5,' redundant orbitals found,',/,
     *  5X,'Number of orbital variations reduced from',I5,
     *  ' to',I5,'.',/,5X,
     *  'Redo gradient calculation with remaining variables.')
#endif
C
C ******************************************************************
C Step 5c: Write diagonal of orbital Hessian to LUIT2
C
      REWIND LUIT2
      WRITE (LUIT2) TABLE(1),TABLE(1),TABLE(1),TABLE(2)
      NWOPX = MAX(4,NWOPT)
      CALL WRITT (LUIT2,NWOPX,WRK(KORDIA))
      WRITE (LUIT2) TABLE(1),TABLE(1),TABLE(1),TABLE(3)
C
C ******************************************************************
C
 5999 CONTINUE
      CALL GETTIM(T5,W5)
C
C ******************************************************************
C Step 6: calculate the CI part of the gradient.
C
      IF (NCONF.GT.1) THEN
         CALL CIDIAG(LSYM,.FALSE.,FCAC,H2AC,INDXCI,G,WRK(KWRK4),LWRK4)
C        CALL CIDIAG(ICSYM,NOH2,FCAC,H2AC,XNDXCI,DIAGC,WRK,LFREE)
         CALL CIGRAD(CREF,FCAC,H2AC,INDXCI,
     *               G,EMCACX,WRK(KWRK4),LWRK4)
C        CALL CIGRAD(CREF,FCAC,H2AC,INDXCI,
C    *               GCI,EMCACT,WRK,LFREE)
         DNCONF = NCONF
         DIFFOK = SQRT(DNCONF)*THRZER
         IF (ABS(EMCACT) .GT. THRZER) THEN
            DIFFOK = DIFFOK * ABS(EMCACT)
         END IF
         IF (NWOPT .GT. 0 .AND. ABS(EMCACT-EMCACX) .GT. DIFFOK) THEN
            IF (ABS(EMCACT-EMCACX) .GT. THRSML) THEN
               WRITE(LUPRI,6010) ' FATAL ERROR', EMCACT, EMCACX
               CALL QUIT(
     &           'FATAL ERROR (GRAD) EMCACT(DV,PV) .ne. EMCACT(CIGRAD)')
            ELSE
               NWARN = NWARN + 1
               WRITE(LUERR,6010) ' WARNING', EMCACT, EMCACX
               WRITE(LUPRI,6010) ' WARNING', EMCACT, EMCACX
            END IF
 6010       FORMAT(1P,//,A,T15,'EMCACT from DV, PV:',D25.15,
     *                 /,  T15,'EMCACT from CIGRAD:',D25.15,/)
         END IF
      ELSE
         G(1)   = D0
      END IF
      CALL GETTIM(T6,W6)
C
C save timing information
C
      NCALLS(  IDTIM) = NCALLS(  IDTIM) + 1
      TIMCPU(1,IDTIM) = TIMCPU(1,IDTIM) + T6 - T0
      TIMCPU(2,IDTIM) = TIMCPU(2,IDTIM) + T1 - T0
      TIMCPU(3,IDTIM) = TIMCPU(3,IDTIM) + T2 - T1
      TIMCPU(4,IDTIM) = TIMCPU(4,IDTIM) + T3 - T2
      TIMCPU(5,IDTIM) = TIMCPU(5,IDTIM) + T4 - T3
      TIMCPU(6,IDTIM) = TIMCPU(6,IDTIM) + T5 - T4
      TIMCPU(7,IDTIM) = TIMCPU(7,IDTIM) + T6 - T5
      TIMWAL(1,IDTIM) = TIMWAL(1,IDTIM) + W6 - W0
      TIMWAL(2,IDTIM) = TIMWAL(2,IDTIM) + W1 - W0
      TIMWAL(3,IDTIM) = TIMWAL(3,IDTIM) + W2 - W1
      TIMWAL(4,IDTIM) = TIMWAL(4,IDTIM) + W3 - W2
      TIMWAL(5,IDTIM) = TIMWAL(5,IDTIM) + W4 - W3
      TIMWAL(6,IDTIM) = TIMWAL(6,IDTIM) + W5 - W4
      TIMWAL(7,IDTIM) = TIMWAL(7,IDTIM) + W6 - W5
C
C     end of GRAD.
C
      CALL QEXIT('GRAD  ')
      RETURN
      END
C  /* Deck grdinf */
      SUBROUTINE GRDINF(GNORM,GTOT)
C
C 18-Jan-1988 hjaaj
C
C ******************************************************************
C Calculate norm of the gradient
C
#include <implicit.h>
#include <infvar.h>
      DIMENSION GNORM(3), GTOT(NVAR)
      PARAMETER ( D0 = 0.0D0 )
C
C Used from common blocks:
C   INFVAR : NCONF, NWOPT, NVAR
C   ITINFO : DINFO(*), IINFO(*)
C
#include <maxorb.h>
#include <itinfo.h>
C
      CALL QENTER('GRDINF')
      IF (NCONF .GT. 1) THEN
         GCINRM = DNORM2(NCONF,GTOT,1)
         INDGCM = IDAMAX(NCONF,GTOT,1)
         GCIMAX = GTOT(INDGCM)
      ELSE
         GCINRM = D0
         INDGCM = 0
         GCIMAX = D0
      END IF
      IF (NWOPT .GT. 0) THEN
         GOBNRM = DNORM2(NWOPT,GTOT(NCONF+1),1)
         INDGOM = IDAMAX(NWOPT,GTOT(NCONF+1),1)
         GORMAX = GTOT(NCONF+INDGOM)
      ELSE
         GOBNRM = D0
         INDGOM = 0
         GORMAX = D0
      END IF
      GNORM(1) = GCINRM
      GNORM(2) = GOBNRM
      GNORM(3) = SQRT(GCINRM**2+GOBNRM**2)
C
C     save information for final summary print-out
C
      IINFO(11) = INDGCM
      IINFO(12) = INDGOM
      DINFO(11) = GNORM(1)
      DINFO(12) = GNORM(2)
      DINFO(13) = GNORM(3)
      DINFO(14) = GCIMAX
      DINFO(15) = GORMAX
C
C     end of GRDINF.
C
      CALL QEXIT('GRDINF')
      RETURN
      END
C  /* Deck orbgrd */
      SUBROUTINE ORBGRD(NWOPX,UDV,FC,FV,FQ,GORB)
C
C Written 29-Nov-1983 by Hans Agren and Hans Jorgen Aa. Jensen
C Revisions
C   10-Oct-1984/13-Jan-1985 hjaaj (act-act rotations)
C    6-Jan-1990 hjaaj
C
C Purpose:
C   To construct the orbital gradient in GORB
C
C MOTECC-90: The algorithms used in this module, ORBGRD, are
C            described in Chapter 8 Appendix 8A of MOTECC-90
C            "Calculation of the Orbital Part of the Gradient
C            Vectors and Sigma Vectors"
C
C Input:
C   The inactive and active Fock matrices FC and FV
C   The active one-electron density matrix UDV
C   The Fock "Q" matrix in FQ
C
C Output:
C   Orbital gradient in GORB
C
C
#include <implicit.h>
      DIMENSION GORB(*), UDV(NASHDI,*)
      DIMENSION FC(*), FV(*), FQ(NORBT,*)
C
C Used from common blocks:
C   INFVAR : NWOPT,NWOPH,JWOP(2,*)
C   INFIND : IROW(*),?
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <infvar.h>
#include <inforb.h>
#include <infind.h>
#include <infdim.h>
#include <infpri.h>
C
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0, D4 = 4.0D0)
C
      CALL QENTER('ORBGRD')
C
      IF (NWOPX .NE. NWOPT .AND. NWOPX .NE. NWOPH) THEN
         WRITE (LUPRI,'(//A/A,2I10)')
     &      ' ORBGRD error, NWOPX .ne. NWOPT,NWOPH',
     &      ' NWOPX, NWOPT, NWOPH :',NWOPX,NWOPT,NWOPH
         CALL QUIT('ORBGRD error, NWOPX .ne. NWOPT,NWOPH')
      END IF
C
      IF (P6FLAG(16) .AND. NASHT .GT. 0) THEN
        WRITE (LUPRI,1000)
        CALL OUTPUT(FQ,1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
      END IF
 1000 FORMAT(/' FQ contribution to Fock matrix')
C
C ** Calculate gradient from the pieces we now have
C
      JSYM = 0
      DO 300 IG = 1,NWOPX
         K = JWOP(1,IG)
         L = JWOP(2,IG)
         ISYM = ISMO(K)
C
         IF (ISYM .NE. JSYM) THEN
            JSYM = ISYM
            NISHI = NISH(ISYM)
            NASHI = NASH(ISYM)
            IASHI = IASH(ISYM)
            IORBI = IORB(ISYM)
            IIORBI = IIORB(ISYM)
         END IF
C
         NK    = K - IORBI
         NL    = L - IORBI
         ITYPK = IOBTYP(K)
         ITYPL = IOBTYP(L)
         IF (ITYPK .EQ. JTINAC) THEN
C ** first index inactive:
            LK = IIORBI + IROW(NL) + NK
            GORB(IG) = GORB(IG) + D4*(FC(LK) + FV(LK))
         ELSE
C ** first index active
C    (only other possibility, inactive, is already eliminated):
            NKW = ISW(K) - NISHT
            TEMP = D0
            NX = NISHI
            DO 100 NXW = IASHI+1,IASHI+NASHI
               NX = NX + 1
               IF (NX.LE.NL) THEN
                  FCXL = FC(IIORBI+IROW(NL)+NX)
               ELSE
                  FCXL = FC(IIORBI+IROW(NX)+NL)
               END IF
               TEMP = TEMP + UDV(NXW,NKW)*FCXL
  100       CONTINUE
            GORB(IG) = GORB(IG) + D2*(TEMP + FQ(L,NKW))
         END IF
         IF (ITYPL .EQ. JTACT) THEN
C ** second index active:
            NLW = ISW(L) - NISHT
            TEMP = D0
            NX = NISHI
            DO 200 NXW = IASHI+1,IASHI+NASHI
               NX = NX + 1
               IF (NX.LE.NK) THEN
                  FCXK = FC(IIORBI+IROW(NK)+NX)
               ELSE
                  FCXK = FC(IIORBI+IROW(NX)+NK)
               END IF
               TEMP = TEMP + UDV(NXW,NLW)*FCXK
  200       CONTINUE
            GORB(IG) = GORB(IG) - D2*(TEMP + FQ(K,NLW))
         END IF
C *** next IG
  300 CONTINUE
C
C Make microcanonical average on GORB
C (AVERAG returns immediately if no averageing requested)
C
      CALL AVERAG(GORB,NWOPT,1)
      CALL QEXIT('ORBGRD')
      RETURN
      END
C  /* Deck ordiac */
      FUNCTION ORDIAC(NU,NV,PV,H2AC)
C
C  Written 11-Oct-1984 by Hans Joergen Aa. Jensen
C  Rewritten 18-Jan-1988 Hans Agren; H2AC is not modified any more.
C  Rerewritten 12-Sep-1989 hjaaj: PV not packed any more.
C    Restrictions NU < NV and sym(NU) = sym(NV) removed.
C
C MOTECC-90: The formulas used in this module, ORDIAC, are
C            described in Chapter 8 Appendix 8C of MOTECC-90
C            "The Diagonal of the Orbital Hessian Matrix"
C
C Purpose:
C  Calculate contributions from active 2-el. integrals
C  to diagonal of orbital Hessian for active-active
C  rotations Lorb,orb(uv,uv) (Eq. A.3 in Chem.Phys. 104(1986)229)
C  (k = min(u,v), l = max(u,v))
C    ORDIAC(u,v) = ORDIAC(k,l)
C    =  sum(x,y) [    4 PV(lx,ly) (kx|ky) +  4 PV(kx,ky) (lx|ly)
C                  -  8 PV(kx,ly) (kx|ly)
C                  +  2 PV(ll,xy) (kk|xy) +  2 PV(kk,xy) (ll|xy)
C                  -  4 PV(kl,xy) (kl|xy) ]
C    =  sum(x>y) [  8 ( PV(lx,ly) (kx|ky) + PV(kx,ky) (lx|ly)
C                     - PV(kx,ly) (kx|ly) - PV(lx,ky) (lx|ky)
C                     - PV(kl,xy) (kl|xy) )
C                +  4 ( PV(ll,xy) (kk|xy) + PV(kk,xy) (ll|xy) ) ]
C    +  sum(x)   [  4 ( PV(lx,lx) (kx|kx) + PV(kx,kx) (lx|lx)
C                     - PV(kl,xx) (kl|xx) )
C                -  8   PV(kx,lx) (kx|lx)
C                +  2 ( PV(ll,xx) (kk|xx) + PV(kk,xx) (ll|xx) ) ]
C
C  Note: PV(uv,xy) = 0.5 * <0| Euv,xy + Euv,yx |0>
C        The P in Appendix A in Chem.Phys. 104(1986)229 is by mistake
C        P(uv,xy) = 0.5 * <0| Euv,xy |0>, in disagreement with
C        the definition in Eq. (3.5).
C
C Input:
C  NU,NV;
C  PV; 2-el. active density matrix
C  H2AC; 2-el. active integrals.
C
#include <implicit.h>
      DIMENSION PV(NNASHX,*), H2AC(NNASHX,*)
C
C Used from common blocks:
C   INFORB : NNASHX,MULD2H,...
C   INFIND : IROW(*)
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infind.h>
C
      PARAMETER (D0=0.0D0, D2=2.0D0, D4 = 4.0D0, D8=8.0D0)
C
      CALL QENTER('ORDIAC')
C
      IF (NU .LT. NV) THEN
         NK = NU
         NL = NV
      ELSE IF (NU .GT. NV) THEN
         NK = NV
         NL = NU
      ELSE
         CALL QUIT('ERROR, ORDIAC called with NU = NV')
      END IF
C
      NKR   = IROW(NK)
      NKK   = NKR + NK
      NLR   = IROW(NL)
      NLK   = NLR + NK
      NLL   = NLR + NL
      KLSYM = MULD2H( NSM(NK) , NSM(NL) )
C
      TERM = D0
      DO 400 NX = 1,NASHT
        ISYMX = NSM(NX)
        ISYMY = MULD2H ( KLSYM , ISYMX )
      IF (ISYMY .GT. ISYMX) GO TO 400
C
        NXR   = IROW(NX)
        NXX   = NXR + NX
C
        IF (NX.LE.NK) THEN
C         NL > NK >= NX
          NKX  = NKR + NX
          NLX  = NLR + NX
        ELSE IF (NX.LT.NL) THEN
C         NL > NX > NK
          NKX  = NXR + NK
          NLX  = NLR + NX
        ELSE
C         NX >= NL > NK
          NKX  = NXR + NK
          NLX  = NXR + NL
        END IF
C
        NYS   = IASH(ISYMY) + 1
        NYE   = MIN( IASH(ISYMY)+NASH(ISYMY) , NX - 1 )
        DO 200 NY = NYS,NYE
          IF (NY.LE.NK) THEN
            NLY  = NLR + NY
            NKY  = NKR + NY
          ELSE IF (NY.LT.NL) THEN
            NLY  = NLR + NY
            NKY  = IROW(NY) + NK
          ELSE
            NLY  = IROW(NY) + NL
            NKY  = IROW(NY) + NK
          END IF
          NXY  = NXR + NY
C
C      +  sum(y<x) [ 8 ( PV(lx,ly) (kx|ky) + PV(kx,ky) (lx|ly)
C                      - PV(kx,ly) (kx|ly) - PV(lx,ky) (lx|ky)
C                      - PV(kl,xy) (kl|xy) )
C                  + 4 ( PV(ll,xy) (kk|xy) + PV(kk,xy) (ll|xy) ) ]
C
          TERM = TERM +  D8*(PV(NLY,NLX)*H2AC(NKY,NKX)
     &                     + PV(NKY,NKX)*H2AC(NLY,NLX)
     &                     - PV(NLY,NKX)*H2AC(NKY,NLX))
     &                     - PV(NKY,NLX)*H2AC(NLY,NKX)
     &                     - PV(NXY,NLK)*H2AC(NXY,NLK)
     &                +  D4*(PV(NXY,NLL)*H2AC(NXY,NKK)
     &                     + PV(NXY,NKK)*H2AC(NXY,NLL))
  200     CONTINUE
C
        IF (ISYMY .EQ. ISYMX) THEN
C
C      +  sum(x)   [  4 ( PV(lx,lx) (kx|kx) + PV(kx,kx) (lx|lx)
C                       - PV(kl,xx) (kl|xx) )
C                  -  8   PV(kx,lx) (kx|lx)
C                  +  2 ( PV(ll,xx) (kk|xx) + PV(kk,xx) (ll|xx) ) ]
C
          TERM = TERM + D4*(PV(NLX,NLX)*H2AC(NKX,NKX)
     &                    + PV(NKX,NKX)*H2AC(NLX,NLX))
     &                - D8* PV(NKX,NLX)*H2AC(NKX,NLX)
     &                + D2*(PV(NXX,NLL)*H2AC(NXX,NKK)
     &                    + PV(NXX,NKK)*H2AC(NXX,NLL))
     &                - D4* PV(NXX,NLK)*H2AC(NXX,NLK)
        END IF
C
  400 CONTINUE
C
      ORDIAC = TERM
C
      CALL QEXIT('ORDIAC')
      RETURN
      END
C  /* Deck ordiag */
      SUBROUTINE ORDIAG(DV,PV,FC,FV,FQ,H2AC,DIAG,AVEODG,IPRINT)
C
C Written by Hans Agren 84-03-27
C Last revision, including active-active part, 11-Oct-1984 Ha
C 13-Feb-1987 hjaaj (revised for RESPONSE, where JWOP(1,i) and
C                    JWOP(2,i) may have different symmetries)
C *** 13-FEB-1987 CHANGE IS ONLY FINISHED FOR HARTREE-FOCK ***
C 19-Sep-1989 hjaaj (removed jwop(1,i) .lt. jwop(2,i) restriction)
C Nov-1989 hjaaj : finished JWOPSY .ne. 1 changes
C                  prepared for AVEODG
C                  New parameter list, moved do 955 section to here
C                  from geth2.
C
C Purpose:
C   This subroutine completes the calculation of the diagonal part of
C   the orbital Hessian by adding the the Fock matrix contributions.
C
C MOTECC-90: The formulas used in this module, ORDIAG, are
C            described in Chapter 8 Appendix 8C of MOTECC-90
C            "The Diagonal of the Orbital Hessian Matrix"
C
C Input:
C   DV , 1-st order density matrix.
C   FC , the inactive Fock matrix.
C   FV , the active Fock matrix.
C
C Output:
C   DIAG , the diagonal of the orbital Hessian.
C
C   This version is directed by the gradient index.
C
#include <implicit.h>
      DIMENSION DV(*),PV(*),FC(*),FV(*),FQ(NORBT,*),H2AC(*),
     &          DIAG(*)
      LOGICAL   AVEODG
C
      PARAMETER ( D2 = 2.0D0, D4 = 4.0D0 )
C
C Used from common blocks:
C   INFORB : NISH(8),NASH(8),...
C   INFIND : ISW(*),ICH(*),?
C   INFVAR : NWOPT,JWOP(2,*),
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infind.h>
#include <infvar.h>
C
      CALL QENTER('ORDIAG')
C
      DO 100 ID = 1,NWOPT
C
         K     = JWOP(1,ID)
         L     = JWOP(2,ID)
         IF (ISW(K) .GT. ISW(L)) THEN
            I = K
            K = L
            L = I
         END IF
         ISYMK = ISMO(K)
         ISYML = ISMO(L)
C
         FCKK  = FC( IIORB(ISYMK) + IROW(K-IORB(ISYMK)+1) )
         FDKK  = FCKK
     *         + FV( IIORB(ISYMK) + IROW(K-IORB(ISYMK)+1) )
         FCLL  = FC( IIORB(ISYML) + IROW(L-IORB(ISYML)+1) )
         FDLL  = FCLL
     *         + FV( IIORB(ISYML) + IROW(L-IORB(ISYML)+1) )
C
         IF (IOBTYP(K).EQ.JTINAC) THEN
            DIAG(ID) = DIAG(ID) + D4*(FDLL - FDKK)
C i*u section
            IF (IOBTYP(L).EQ.JTACT) THEN
               NISHL  = NISH(ISYML)
               NASHL  = NASH(ISYML)
               IASHL  = IASH(ISYML)
               IIORBL = IIORB(ISYML)
               NL     = ICH(L)
C loop for F(l,l) = sum(u) DV(l,u)FC(u,l) + FQ(l,l)
               FOCKLL  = FQ(L,NL)
               DO 105 NU = IASHL+1,IASHL+NASHL
                  IF (NL.GE.NU) THEN
                     FOCKLL = FOCKLL + DV(IROW(NL)+NU) *
     *               FC(IIORBL+IROW(NISHL-IASHL+NL)+NISHL-IASHL+NU)
                  ELSE
                     FOCKLL = FOCKLL + DV(IROW(NU)+NL) *
     *               FC(IIORBL+IROW(NISHL-IASHL+NU)+NISHL-IASHL+NL)
                  END IF
  105          CONTINUE
C
               DVLL     = DV( IROW( NL + 1 ) )
               IF (AVEODG) THEN
                  DIAG(ID) = DIAG(ID) + D2*DVLL*FDKK - D2*FOCKLL
               ELSE
                  DIAG(ID) = DIAG(ID) + D2*DVLL*FCKK - D2*FOCKLL
               END IF
C
            END IF
C
C
         ELSE IF (IOBTYP(K).EQ.JTACT) THEN
            NISHK  = NISH(ISYMK)
            NASHK  = NASH(ISYMK)
            IASHK  = IASH(ISYMK)
            IIORBK = IIORB(ISYMK)
            NK     = ICH(K)
            IF (IOBTYP(L).EQ.JTSEC) THEN
C u*a section
C
C loop for F(k,k) = sum(u) DV(k,u)FC(u,k) + FQ(k,k)
               FOCKKK = FQ(K,NK)
               DO 205 NU = IASHK+1,IASHK+NASHK
                  IF (NK.GE.NU) THEN
                     FOCKKK = FOCKKK + DV(IROW(NK)+NU) *
     *               FC(IIORBK+IROW(NISHK-IASHK+NK)+NISHK-IASHK+NU)
                  ELSE
                     FOCKKK = FOCKKK + DV(IROW(NU)+NK) *
     *               FC(IIORBK+IROW(NISHK-IASHK+NU)+NISHK-IASHK+NK)
                  END IF
  205          CONTINUE
C
               DVKK = DV(IROW(NK+1))
               IF (AVEODG) THEN
                  DIAG(ID) = DIAG(ID) + D2*DVKK*FDLL - D2*FOCKKK
               ELSE
                  DIAG(ID) = DIAG(ID) + D2*DVKK*FCLL - D2*FOCKKK
               END IF
C
            ELSE IF (IOBTYP(L).EQ.JTACT) THEN
C u*v section
               NISHL  = NISH(ISYML)
               NASHL  = NASH(ISYML)
               IASHL  = IASH(ISYML)
               IIORBL = IIORB(ISYML)
               NL     = ICH(L)
               IORBK  = IORB(ISYMK)
C -- 2 loops for F(k,k) + F(l,l) elements
C where F(t,t) = sum(u) DV(t,u)FC(u,t) + FQ(t,t)
               FCKSUM = FQ(K,NK) + FQ(L,NL)
               DO 305 NU = IASHK+1,IASHK+NASHK
                  IF (NK.GE.NU) THEN
                     FCKSUM = FCKSUM + DV(IROW(NK)+NU) *
     *               FC(IIORBK+IROW(NISHK-IASHK+NK)+NISHK-IASHK+NU)
                  ELSE
                     FCKSUM = FCKSUM + DV(IROW(NU)+NK) *
     *               FC(IIORBK+IROW(NISHK-IASHK+NU)+NISHK-IASHK+NK)
                  END IF
  305          CONTINUE
               DO 315 NU = IASHL+1,IASHL+NASHL
                  IF (NL.GE.NU) THEN
                     FCKSUM = FCKSUM + DV(IROW(NL)+NU) *
     *               FC(IIORBL+IROW(NISHL-IASHL+NL)+NISHL-IASHL+NU)
                  ELSE
                     FCKSUM = FCKSUM + DV(IROW(NU)+NL) *
     *               FC(IIORBL+IROW(NISHL-IASHL+NU)+NISHL-IASHL+NL)
                  END IF
  315          CONTINUE
C -- no AVEODG here, u*v setion is always calculated exactly
C -- Non-Fock contributions from active 2-el. integrals to active-active
C  part of orbital Hessian diagonal are calculated and added using
C  FUNCTION ORDIAC.
               DIAG(ID) = DIAG(ID) + D2*( DV(IROW(NK+1)) * FCLL
     *                                  + DV(IROW(NL+1)) * FCKK )
     *                  - D2*FCKSUM + ORDIAC(NK,NL,PV,H2AC)
               IF (ISYMK .EQ. ISYML) THEN
                  FCKL   = FC( IIORBK + IROW(L-IORBK) + K - IORBK )
                  DIAG(ID) = DIAG(ID) - D4 * DV(IROW(NL)+NK) * FCKL
               END IF
C
            END IF
         END IF
C
  100 CONTINUE
C
C
C*******************************************************************
C
C Make microcanonical average on diagonal
C (AVERAG returns immediately if no averageing requested)
C
      CALL AVERAG(DIAG,NWOPT,1)
C
      CALL QEXIT('ORDIAG')
      RETURN
      END
C  /* Deck or1dia */
      SUBROUTINE OR1DIA(DV,FP,DIAG,IPRINT)
C
C 6-May-1990 Hans Joergen Aa. Jensen
C (based on ORDIAG)
C
C Purpose:
C   Orbital Hessian diagonal for a one-electron operator (e.g. solvent)
C   is added to DIAG.
C   Symmetry is determined implicitly through the JWOP operator.
C
C Input:
C   DV(NNASHX) , 1-st order density matrix.
C   FP(NNBASX) , the property matrix for the one-electron operator
C
C Output:
C   DIAG , the diagonal of the orbital Hessian.
C
#include <implicit.h>
      DIMENSION DV(*),FP(*), DIAG(*)
C
      PARAMETER ( D0 = 0.0D0, D2 = 2.0D0, D4 = 4.0D0 )
C
C Used from common blocks:
C   INFORB : NISH(8),NASH(8),...
C   INFIND : ISW(*),ICH(*),?
C   INFVAR : NWOPT,JWOP(2,*),
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infind.h>
#include <infvar.h>
C
      CALL QENTER('OR1DIA')
C
      DO 100 ID = 1,NWOPT
C
         K     = JWOP(1,ID)
         L     = JWOP(2,ID)
         IF (ISW(K) .GT. ISW(L)) THEN
            I = K
            K = L
            L = I
         END IF
         ISYMK = ISMO(K)
         ISYML = ISMO(L)
C
         FPKK  = FP( IROW(K+1) )
         FPLL  = FP( IROW(L+1) )
C
         IF (IOBTYP(K).EQ.JTINAC) THEN
            DIAG(ID) = DIAG(ID) + D4*(FPLL - FPKK)
C i*u section
            IF (IOBTYP(L).EQ.JTACT) THEN
               IASHL  = IASH(ISYML)
               NASHL  = NASH(ISYML)
               IOFFL  = IORB(ISYML) + NISH(ISYML) - IASHL
               NL     = ICH(L)
C loop for F(l,l) = sum(u) DV(l,u)FP(u,l)
               FOCKLL  = D0
               DO 105 NU = IASHL+1,IASHL+NASHL
                  IF (NL.GE.NU) THEN
                     FOCKLL = FOCKLL + DV(IROW(NL)+NU) *
     *                  FP(IROW(IOFFL+NL) + IOFFL+NU)
                  ELSE
                     FOCKLL = FOCKLL + DV(IROW(NU)+NL) *
     *                  FP(IROW(IOFFL+NU) + IOFFL+NL)
                  END IF
  105          CONTINUE
C
               DIAG(ID) = DIAG(ID) + D2*DV(IROW(NL+1))*FPKK - D2*FOCKLL
C
            END IF
C
C
         ELSE IF (IOBTYP(K).EQ.JTACT) THEN
            IASHK  = IASH(ISYMK)
            NASHK  = NASH(ISYMK)
            IOFFK  = IORB(ISYMK) + NISH(ISYMK) - IASHK
            NK     = ICH(K)
            IF (IOBTYP(L).EQ.JTSEC) THEN
C u*a section
C
C loop for F(k,k) = sum(u) DV(k,u)FP(u,k)
               FOCKKK = D0
               DO 205 NU = IASHK+1,IASHK+NASHK
                  IF (NK.GE.NU) THEN
                     FOCKKK = FOCKKK + DV(IROW(NK)+NU) *
     *                  FP(IROW(IOFFK+NK) + IOFFK+NU)
                  ELSE
                     FOCKKK = FOCKKK + DV(IROW(NU)+NK) *
     *                  FP(IROW(IOFFK+NU) + IOFFK+NK)
                  END IF
  205          CONTINUE
C
               DIAG(ID) = DIAG(ID) + D2*DV(IROW(NK+1))*FPLL - D2*FOCKKK
C
            ELSE IF (IOBTYP(L).EQ.JTACT) THEN
C u*v section
               IASHL  = IASH(ISYML)
               NASHL  = NASH(ISYML)
               IOFFL  = IORB(ISYML) + NISH(ISYML) - IASHL
               NL     = ICH(L)
C -- 2 loops for F(k,k) + F(l,l) elements
C where F(t,t) = sum(u) DV(t,u)FP(u,t)
               FCKSUM = D0
               DO 305 NU = IASHK+1,IASHK+NASHK
                  IF (NK.GE.NU) THEN
                     FCKSUM = FCKSUM + DV(IROW(NK)+NU) *
     *                  FP(IROW(IOFFK+NK) + IOFFK+NU)
                  ELSE
                     FCKSUM = FCKSUM + DV(IROW(NU)+NK) *
     *                  FP(IROW(IOFFK+NU) + IOFFK+NK)
                  END IF
  305          CONTINUE
               DO 315 NU = IASHL+1,IASHL+NASHL
                  IF (NL.GE.NU) THEN
                     FCKSUM = FCKSUM + DV(IROW(NL)+NU) *
     *                  FP(IROW(IOFFL+NL) + IOFFL+NU)
                  ELSE
                     FCKSUM = FCKSUM + DV(IROW(NU)+NL) *
     *                  FP(IROW(IOFFL+NU) + IOFFL+NL)
                  END IF
  315          CONTINUE
               DIAG(ID) = DIAG(ID) + D2*( DV(IROW(NK+1)) * FPLL
     *                                  + DV(IROW(NL+1)) * FPKK )
     *                  - D2*FCKSUM
C
               IF (K .LT. L) THEN
                  DIAG(ID) = DIAG(ID)
     &                     - D4 * DV(IROW(NL)+NK) * FP(IROW(L)+K)
               ELSE
                  DIAG(ID) = DIAG(ID)
     &                     - D4 * DV(IROW(NK)+NL) * FP(IROW(K)+L)
               END IF
C
            END IF
         END IF
C
  100 CONTINUE
C
C
C*******************************************************************
C
C Make microcanonical average on diagonal
C (AVERAG returns immediately if no averageing requested)
C
      CALL AVERAG(DIAG,NWOPT,1)
C
      CALL QEXIT('OR1DIA')
      RETURN
      END
C  /* Deck updgrd */
      SUBROUTINE UPDGRD (ICTL,CREF,CMO,DV,PV,G,EMCMY,EMCACT,
     *                   FC,FV,FQ,DIAOR,INDXCI,WRK,LFREE)
C
C 14-May-1985 Hans Jorgen Aa. Jensen
C (based on GRAD)
C Last revision: 7-Nov-1985 hjaaj (ver 2.0)
C               20-Jan-1988 hjaaj (nasht=1 case)
C                9-Oct-1989 hjaaj (indxci in param.list)
C
C Purpose:
C   This subroutine computes the energy and the gradient
C   (in C representation).
C   UPDGRD will be called when using update methods or
C   orbital absorption in the MCSCF optimization.
C
C Input:
C   ICTL = 0, new CI vector (calculate new DV and PV)
C        = 1, old CI vector (i.e. DV and PV are correct)
C        = 2, as = 1 but only calculate orbital gradient.
C        all above are for SIRUPD.
C   ICTL = 3, calculate orbital gradient for SIRORB.
C   CREF, the CI vector for the current expansion point (CEP).
C   CMO, the orbitals for the CEP
C
C Output:
C   G, gradient of the CEP.
C   EMCMY, inactive energy of the CEP.
C   EMCACT, active energy of the CEP.
C   DV, one-electron d.m. for CREF.
C   PV, two-electron d.m. for CREF.
C   FQ, Fock Q matrix
C   DIAOR, diagonal of orbital Hessian, if ICTL = 3
C
C Scratch:
C   WRK
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION CREF(*),CMO(*),DV(*),PV(*), G(*),
     *          FC(*),FV(*),FQ(*), DIAOR(*), INDXCI(*), WRK(LFREE)
#include <iratdef.h>
C
C -- local variables and constants
C
      LOGICAL DODVPV, DOCIGR, FULODG, AVEODG
      PARAMETER (D0 = 0.0D0 , D1 = 1.0D0)
      PARAMETER (THRSML = 1.D-8)
#include <thrzer.h>
C
C Used from common blocks:
C   INFINP : FLAG(*),LSYM,?
C   INFVAR : NCONF,NWOPT,NVAR,JWOP(2,*)
C   INFTIM : NCALLS,TIMCPU,TIMWAL;   IDTIM above is index for these
C   INFOPT : GNORM(3),?
C   -- GNORM, norm of the CI, the orbital, and the total gradient.
C   INFDIM : ?
C   INFTAP : LUIT2,?
C
#include <maxorb.h>
#include <infinp.h>
#include <infvar.h>
#include <inforb.h>
#include <infdim.h>
#include <infopt.h>
#include <inftap.h>
#include <infpri.h>
#include <itinfo.h>
C ---
      PARAMETER (IDTIM = 3)
C     ... GRAD has IDTIM = 1
#include <inftim.h>
C
C
      CALL QENTER('UPDGRD')
      CALL GETTIM (T0,W0)
C
      IF (ADDSRI) THEN
         call quit('UPDGRD not implemented yet for srDFT')
      END IF
C
      IF (ICTL.EQ.0) THEN
         DODVPV = .TRUE.
         DOCIGR = .TRUE.
      ELSE IF (ICTL.EQ.1) THEN
         DODVPV = .FALSE.
         DOCIGR = .TRUE.
      ELSE IF ( (ICTL.EQ.2 .OR. ICTL.EQ.3) .AND. NWOPT .GT. 0 ) THEN
         DODVPV = .FALSE.
         DOCIGR = .FALSE.
      ELSE
         WRITE (LUPRI,'(/A,I10)')
     *      ' *** SOFTWARE ERROR (UPDGRD), invalid ICTL =',ICTL
         IF ( (ICTL.EQ.2 .OR. ICTL.EQ.3) .AND. NWOPT .EQ. 0 )
     *      WRITE (LUPRI,'(/2A)')
     *      ' This ICTL requests the orbital gradient only,',
     *      ' but number of orbitals rotations is zero.'
         CALL QTRACE(LUPRI)
         CALL QUIT('*** SOFTWARE ERROR (UPDGRD), invalid ICTL.')
      END IF
C
C     Control evaluation of orbital diagonal.
C     IORDIA = 0 : do not calculate orbital diagonal
C              1 : calculate approximate orbital diagonal
C              2 : calculate exact orbital diagonal
      IF (ICTL.EQ.3) THEN
C        sirorb: calculate orbital diagonal
C        (else sirupd: do not calculate orbital diagonal)
         IF (NWOPT .GT. 0) CALL DZERO(DIAOR,NWOPT)
         IF (FLAG(12)) THEN
            FULODG = .TRUE.
            AVEODG = .FALSE.
         ELSE
            FULODG = .FALSE.
            AVEODG = .TRUE.
         END IF
      ELSE
         FULODG = .FALSE.
C        set FULODG for GETH2 for sirupd
      END IF
C
      IF (NASHT .LE. 1) THEN
         IF (NASHT .EQ. 1 .AND. DODVPV) THEN
C        ... one open shell special case
            DV(1) = D1
            PV(1) = D0
         END IF
         DODVPV = .FALSE.
         DOCIGR = .FALSE.
      END IF
C
C ******************************************************************
C Step 0: WRK space allocation:
C
C .. 1.1 .. (MAKDM, CIGRAD)
C .. 1.2 .. (FCKMAT)
      KFCAC  = 1
      KWRKF  = KFCAC + NNASHX
      LWRKF  = LFREE - KWRKF
C .. 1.3 .. (GETH2,  ORBGRD, CIGRAD)
      KH2AC  = KWRKF
      KWRK3  = KH2AC + NNASHX*NNASHX
      LWRK3  = LFREE - KWRK3
C .. 1.3.1 .. (ORBGRD)
      KUDV   = KWRK3
      KWRK   = MAX(KWRKF, KUDV + N2ASHX)
C .. 1.4 .. (CIGRAD)
      KWRK4  = KWRK3
      LWRK4  = LFREE- KWRK4
C Test if WRK is big enough:
      IF (KWRK.GT.(LFREE+1)) CALL ERRWRK('UPDGRD',-KWRK,LFREE)
C
C ******************************************************************
C Step 1: calculate one- and two-electron density matrices for
C         active orbital indices.
C
      IF (.NOT. DODVPV) GO TO 1999
      IF (NWOPT.EQ.0) THEN
         CALL DZERO(DV,NNASHX)
         CALL DZERO(PV,NNASHX*NNASHX)
         GO TO 1999
      END IF
      CALL MAKDM(CREF,DV,PV,INDXCI,WRK(KFCAC),(LFREE-KFCAC))
C     CALL MAKDM(CREF,DV,PV,INDXCI,WRK,LFREE)
C
      IF (P6FLAG(15)) THEN
         WRITE (LUPRI,1100)
         CALL OUTPAK(DV,NASHT,1,LUPRI)
      END IF
      IF (P6FLAG(25)) THEN
         WRITE (LUPRI,1200)
         CALL OUTPUT(PV,1,NNASHX,1,NNASHX,NNASHX,NNASHX,-1,LUPRI)
      END IF
 1100 FORMAT(/' DV = One-el. density matrix, active part, MO-basis')
 1200 FORMAT(/' PV = Two-el. density matrix, active part, MO-basis')
 1999 CONTINUE
      CALL GETTIM(T1,W1)
C
C ******************************************************************
C Step 2: calculate inactive and active Fock matrices.
C
      CALL FCKMAT((NASHT.EQ.0),DV,CMO,EMCMY,FC,FV,WRK(KWRKF),LWRKF)
C     CALL FCKMAT(ONLYFC,DV,CMO,EMCMY,FC,FV,WRK,LFREE)
      IF (NASHT .EQ. 0) CALL DZERO(FV,NNORBT)
C     ... currently ORBGRD and ORBSIG always references FV
C
      IF (NASHT .GE. 1) THEN
         CALL GETAC(FC,WRK(KFCAC))
         IF (P6FLAG(21)) THEN
            WRITE(LUPRI,'(/A)')
     *         ' Active block of inactive Fock matrix FCAC'
            CALL OUTPAK (WRK(KFCAC),NASHT,1,LUPRI)
         END IF
      END IF
      CALL GETTIM(T2,W2)
C
C ******************************************************************
C Step 3: calculate Fock "Q" matrix and two-electron integrals with
C         active orbital indices.
C
C
      IF (NASHT .GT. 1) THEN
         CALL GETH2(DV,PV,WRK(KH2AC),FQ,DIAOR,FULODG,
     &              WRK(KWRK3),LWRK3)
C        CALL GETH2(DV,PV,H2AC,FQ,ORDIAG,FULODG,WRK,LFRSAV)
         IF (P6FLAG(23) .AND. IPRI6.GT.60 .AND. FULODG) THEN
            WRITE (LUPRI,'(4(/A))')
     &         ' Non-Fock type contributions to inact-act and act-sec',
     &         ' ----------------------------------------------------',
     &         ' diagonal of orbital Hessian',
     &         ' ---------------------------'
            CALL PRKAP (NWOPT,DIAOR,-0.1D0,LUPRI)
         END IF
C
         IF (P6FLAG(26)) THEN
            WRITE(LUPRI,'(/A)') ' Two-electron integrals H2AC (UPDGRD)'
            CALL OUTPUT(WRK(KH2AC),1,NNASHX,1,NNASHX,
     *                  NNASHX,NNASHX,-1,LUPRI)
         END IF
C
C        *** calculate active energy, before we (maybe) modify FCAC
C
         EMCACT = ENRACT(DV,WRK(KFCAC),FQ)
      ELSE IF (NASHT .EQ. 1) THEN
         CALL DZERO(FQ,NORBT)
         EMCACT = ENRACT(DV,WRK(KFCAC),FQ)
      ELSE
         EMCACT = D0
      END IF
      CALL GETTIM(T3,W3)
C
C ******************************************************************
C Step 4: calculate Fock matrix contributions to the diagonal of the
C         orbital Hessian.
C
      IF (ICTL .EQ. 3) THEN
         CALL ORDIAG(DV,PV,FC,FV,FQ,WRK(KH2AC),DIAOR,
     *               AVEODG,IPRSIR)
C        CALL ORDIAG(DV,PV,FC,FV,FQ,H2AC,DIAOR,AVEODG,IPRSIR)
C
         IF (P6FLAG(23)) THEN
            IF (FLAG(12)) THEN
               WRITE (LUPRI,'(/A/A)')
     &         ' Diagonal of orbital Hessian',
     &         ' ---------------------------'
            ELSE
               WRITE (LUPRI,'(/A/A)')
     &         ' Approximate diagonal of orbital Hessian',
     &         ' ---------------------------------------'
            END IF
            IF (SUPSYM) WRITE (LUPRI,'(/A)')
     &         ' (Diagonal is microcanonically averaged.)'
            CALL PRKAP (NWOPT,DIAOR,-0.1D0,LUPRI)
         END IF
      END IF
      CALL GETTIM(T4,W4)
C
C
C ******************************************************************
C Step 5: put the pieces together for the orbital gradient.
C
C
      IF (ICTL .EQ. 3) THEN
         IGOFF = 0
      ELSE
         IGOFF = NCONF
      END IF
      IF (NWOPT .GT. 0) THEN
C        Unpack DV into UDV ( SP to SI format)
         IF (NASHT .GT. 0) CALL DSPTSI(NASHT,DV,WRK(KUDV))
C        Initialize gradient
         CALL DZERO(G(IGOFF+1),NWOPH)
         CALL ORBGRD(NWOPH,WRK(KUDV),FC,FV,FQ,G(IGOFF+1))
C        CALL ORBGRD(NWOPX,UDV,FC,FV,FQ,GORB)
      END IF
      CALL GETTIM(T5,W5)
C
C
C ******************************************************************
C Step 6: calculate the CI part of the gradient.
C         LUIT2 temp. set to zero, to tell CIDIAG not to write
C         CI diagonals to LUIT2.
C
      IF (DOCIGR) THEN
         LUIT2S = LUIT2
         LUIT2  = -1
         CALL CIDIAG(LSYM,.FALSE.,WRK(KFCAC),
     *               WRK(KH2AC),INDXCI,G,WRK(KWRK4),LWRK4)
C        CALL CIDIAG(ICSYM,NOH2,FCAC,H2AC,INDXCI,DIAGC,WRK,LFREE)
         LUIT2  = LUIT2S
C
         CALL CIGRAD(CREF,WRK(KFCAC),WRK(KH2AC),
     *               INDXCI,G,EMCACX,WRK(KWRK4),LWRK4)
C        CALL CIGRAD(CREF,FCAC,H2AC,INDXCI,
C    *               GCI,EMCACT,WRK,LFREE)
         DNCONF = NCONF
         DIFFOK = SQRT(DNCONF)*THRZER
         IF (ABS(EMCACT) .GT. THRZER) THEN
            DIFFOK = DIFFOK * ABS(EMCACT)
         END IF
         IF (NWOPT .GT. 0 .AND. ABS(EMCACT-EMCACX) .GT. DIFFOK) THEN
            IF (ABS(EMCACT-EMCACX) .GT. THRSML) THEN
               WRITE(LUPRI,6010) ' FATAL ERROR', EMCACT, EMCACX
               CALL QUIT(
     &            'FATAL ERROR, EMCACT(DV, PV) .ne. EMCACT(CIGRAD)')
            ELSE
               NWARN = NWARN + 1
               WRITE(LUERR,6010) ' WARNING', EMCACT, EMCACX
               WRITE(LUPRI,6010) ' WARNING', EMCACT, EMCACX
            END IF
 6010       FORMAT(1P,//,A,T15,'EMCACT from DV, PV:',D25.15,
     *                 /,  T15,'EMCACT from CIGRAD:',D25.15,/)
         END IF
      ELSE
         IF (IGOFF .EQ. 1) G(1) = D0
      END IF
      CALL GETTIM(T6,W6)
C
C ******************************************************************
C
C save timing information
C
      NCALLS(  IDTIM) = NCALLS(  IDTIM) + 1
      TIMCPU(1,IDTIM) = TIMCPU(1,IDTIM) + T6 - T0
      TIMCPU(2,IDTIM) = TIMCPU(2,IDTIM) + T1 - T0
      TIMCPU(3,IDTIM) = TIMCPU(3,IDTIM) + T2 - T1
      TIMCPU(4,IDTIM) = TIMCPU(4,IDTIM) + T3 - T2
      TIMCPU(5,IDTIM) = TIMCPU(5,IDTIM) + T4 - T3
      TIMCPU(6,IDTIM) = TIMCPU(6,IDTIM) + T5 - T4
      TIMCPU(7,IDTIM) = TIMCPU(7,IDTIM) + T6 - T5
      TIMWAL(1,IDTIM) = TIMWAL(1,IDTIM) + W6 - W0
      TIMWAL(2,IDTIM) = TIMWAL(2,IDTIM) + W1 - W0
      TIMWAL(3,IDTIM) = TIMWAL(3,IDTIM) + W2 - W1
      TIMWAL(4,IDTIM) = TIMWAL(4,IDTIM) + W3 - W2
      TIMWAL(5,IDTIM) = TIMWAL(5,IDTIM) + W4 - W3
      TIMWAL(6,IDTIM) = TIMWAL(6,IDTIM) + W5 - W4
      TIMWAL(7,IDTIM) = TIMWAL(7,IDTIM) + W6 - W5
C
C     end of UPDGRD.
C
      CALL QEXIT('UPDGRD')
      RETURN
      END
C  /* Deck enract */
      FUNCTION ENRACT(DV,FCAC,FQ)
C
C PURPOSE:
C  DETERMINE ACTIVE ELECTRONIC ENERGY OF MCSCF STATE
C
C  EACTIV = SUM(U,V) DV(U,V)*FI(U,V) + SUM(UVXY) P(UV,XY)*(UV/XY)
C
C         = SUM(U,V) DV(U,V)*FI(U,V) + 1/2 SUM(U) FQ(U,U)
C
#include <implicit.h>
C
      DIMENSION DV(*),FCAC(*),FQ(NORBT,*)
      PARAMETER ( DP5 = 0.5D0,  D2=2.0D0 )
C
C Used from common blocks:
C  INFORB : NORBT,NNASHX,NASHT,NISHT
C  INFIND : IROW(),ISX()
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infind.h>
C
      EACT = D2 * DDOT(NNASHX,DV,1,FCAC,1)
C     diagonal elements counted twice now.
      DO 100 NI = 1,NASHT
         NII = IROW(NI+1)
         I   = ISX(NISHT+NI)
         EACT = EACT + DP5*FQ(I,NI) - DV(NII)*FCAC(NII)
 100  CONTINUE
      ENRACT = EACT
      RETURN
      END
C  /* Deck adh2ac */
      SUBROUTINE ADH2AC(H2ACXY,H2XY,IUVSYM)
C
C 891013-hjaaj
C
C Add contributions from (xy) distribution to
C H2AC(uv,xy) from H2XY(uv); (uv) distribution has symmetry IUVSYM
C
C Subroutine can also be used for contributions to H2XAC.
C
#include <implicit.h>
      DIMENSION H2ACXY(*),H2XY(NORBT,NORBT)
C
C Used from common blocks:
C   INFORB : NSYM, MULD2H, NORBT, NASHT, IASH, NASH
C   INFIND : ISX, NSM, IROW
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infind.h>
C
      DO 200 NVW = 1,NASHT
         IVSYM = NSM(NVW)
         IUSYM = MULD2H(IVSYM, IUVSYM)
      IF (IUSYM .GT. IVSYM) GO TO 200
         NVWROW= IROW(NVW)
         IV    = ISX(NISHT+NVW)
         NUWST = IASH(IUSYM) + 1
         NUWEND= MIN(NVW,IASH(IUSYM)+NASH(IUSYM))
         DO 100 NUW = NUWST,NUWEND
            IU   = ISX(NISHT+NUW)
            NUVW = NVWROW + NUW
            H2ACXY(NUVW) = H2ACXY(NUVW) + H2XY(IU,IV)
  100    CONTINUE
  200 CONTINUE
      RETURN
      END
C  /* Deck addfq */
      SUBROUTINE ADDFQ(NXW,NYW,FQ,H2XY,PVXY,JH2SYM)
C
C 891006-hjaaj
C 920424-hjaaj (revised for general PVXY, i.e. not symmetric)
C
C PURPOSE:
C  Add the contributions to FQ which are determined
C  from integral distribution  (**,XY) (Mulliken notation)
C
C  FQ(p,u) = FQ(p,u) + PV(uv,xy) (pv,xy)
C
#include <implicit.h>
      DIMENSION FQ(NORBT,NASHT), PVXY(NASHT,NASHT), H2XY(NORBT,NORBT)
C
C Used from common blocks:
C   INFORB : NORBT, NASHT, ...
C   INFVAR : JWOPSY
C   INFIND : NSM
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infvar.h>
#include <infind.h>
C
C     Fock q matrix has symmetry jwopsy
C     Integrals from H2 : (pv|xy), symmetry jh2sym
C     Elements  from PV : [uv|xy], symmetry muld2h(jwopsy,jh2sym)
C
      IXYSYM = MULD2H ( NSM(NXW) , NSM(NYW) )
      IPVSYM = MULD2H ( IXYSYM   , JH2SYM   )
C
      DO 100 IUSYM = 1,NSYM
         IPSYM = MULD2H(IUSYM,JWOPSY)
         IVSYM = MULD2H(IPVSYM,IPSYM)
         NORBP = NORB(IPSYM)
         NASHU = NASH(IUSYM)
         NASHV = NASH(IVSYM)
         IF ( (NORBP.GT.0) .AND. (NASHU.GT.0) .AND. (NASHV.GT.0) ) THEN
            IPST  = IORB(IPSYM) + 1
            IVST  = IORB(IVSYM) + 1 + NISH(IVSYM)
            NUWST = IASH(IUSYM) + 1
            NVWST = IASH(IVSYM) + 1
C
C add SUM(V) (PV,XY) * [UV,XY]  to FQ(P,U)
C NOTE: PVXY may not be symmetric!
C
            CALL DGEMM('N','T',NORBP,NASHU,NASHV,1.D0,
     &                 H2XY(IPST,IVST),NORBT,
     &                 PVXY(NUWST,NVWST),NASHT,1.D0,
     &                 FQ(IPST,NUWST),NORBT)
C
         END IF
  100 CONTINUE
C
C END OF ADDFQ
C
      RETURN
      END
C  /* Deck sirpvd */
      SUBROUTINE SIRPVD(K,L,PVKL,PV,IPVDIS)
C
C GET 2-ELECTRON DENSITY DISTRIBUTIONS OF VARIOUS TYPE FROM PV
C
C     IPVDIS = 1  [**,KL]*2*(1+DELTA(K,L))-1 IN PVKL[*,*]
C                 TRIANGULAR PACKED DISTRIBUTIONS IN PV
C
C     IPVDIS = 2  [*K,*L] IN PVKL[*,*]
C                  P  Q            Q P
C                 NOTE ORDERING OF Q AND P IN PVKL
C                 TRIANGULAR PACKED DISTRIBUTIONS IN PV
C
#include <implicit.h>
      DIMENSION PVKL(NASHT,*),PV(NNASHX,NNASHX)
C
      PARAMETER ( D2 = 2.0D0 )
C
C Used from common blocks:
C   INFORB : NNASHX, NASHT, N2ASHX
C   INFPRI : IPRSIR
C
#include <priunit.h>
#include <inforb.h>
#include <infpri.h>
C
C
      IPRLIN = IPRSIR
      IF ( IPRLIN.GT.250 ) THEN
         WRITE(LUPRI,'(/A)') ' ********* SIRPVD **********'
         WRITE(LUPRI,'(/A)') ' TWO-BODY DENSITY MATRIX'
         CALL OUTPUT(PV(1,1),1,NNASHX,1,NNASHX,
     *               NNASHX,NNASHX,-1,LUPRI)
      END IF
C
      GO TO (1,2),IPVDIS
         WRITE(LUPRI,'(/A,I5)')
     *      ' SIRPVD: INCORRECT SPECIFICATION OF IPVDIS ,IPVDIS:',IPVDIS
         CALL QUIT('SIRPVD: INCORRECT SPECIFICATION OF IPVDIS')
 1    CONTINUE
C
C     IPVDIS = 1  [**,KL]*2*(1+DELTA(K,L))-1 IN PVKL[*,*]
C                 TRIANGULAR PACKED DISTRIBUTIONS IN PV
      IF (K.GT.L) THEN
         IKL = K*(K-1)/2 + L
      ELSE
         IKL = L*(L-1)/2 + K
      END IF
      CALL DSPTSI(NASHT,PV(1,IKL),PVKL)
      GO TO 100
 2    CONTINUE
C
C     IPVDIS = 2  [*K,*L] IN PVKL[*,*]
C                  P  Q            Q P
C                 NOTE ORDERING OF Q AND P IN PVKL
C                 TRIANGULAR PACKED DISTRIBUTIONS IN PV
C
      DO 2000 IP = 1,NASHT
         IF (IP.GT.K) THEN
            IPK = IP*(IP-1)/2 + K
         ELSE
            IPK = K*(K-1)/2   + IP
         END IF
         DO 2010 IQ = 1,NASHT
            IF (IQ.GT.L) THEN
               IQL = IQ*(IQ-1)/2 + L
            ELSE
               IQL = L*(L-1)/2   + IQ
            END IF
            PVKL (IQ,IP) = PV(IPK,IQL) * D2
 2010    CONTINUE
 2000 CONTINUE
      GO TO 100
 100  CONTINUE
      IF (IPRLIN.GT.90) THEN
         WRITE(LUPRI,'(/A,I5,A,2I5)')
     &      ' SIRPVD distribution type',IPVDIS,
     &      '     distribution: K,L',K,L
         CALL OUTPUT(PVKL,1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck fckden2 */
      SUBROUTINE FCKDEN2(GETDC,GETDV,DCAO,DVAO,CMO,DV,WRK,LFRSAV)
C
C Jan. 1990 Hans Joergen Aa. Jensen
C l.r. 1997 hjaaj (d.m. now square matrices)
C
C Backtransform one-electron density matrices for
C calculation of Fock matrices in AO basis
C ( contravariant transformation ).
C
C Input:
C   GETDC   if true calculate DCAO
C   GETDV   if true calculate DVAO
C   CMO(*)  molecular orbital coefficients
C   DV(*)   active part of one-electron density matrix
C           (over MOs) unpacked
C
C Scratch:
C   WRK(LFRSAV)
C
#include <implicit.h>
      DIMENSION CMO(*),DV(*),DCAO(NBAST,*),DVAO(NBAST,*), WRK(*)
      LOGICAL   GETDC, GETDV
C
      PARAMETER (D2 = 2.0D0)
C
C Used from common blocks:
C  INFORB : NSYM,MULD2H(8,8),NASHT,NNBASX,N2BASX,...
C  INFVAR : JWOPSY
C  INFPRI : P6FLAG()
C
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infvar.h>
#include <infpri.h>
C
C
      CALL QENTER('FCKDEN2')
C
C
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LFRSAV
C
C     ************************************
C     DVAO matrix
C
      IF (GETDV) THEN
         CALL DZERO(DVAO,N2BASX)
      IF (NASHT .GT. 0) THEN
         CALL MEMGET('REAL',KDAO1,NASHT*NBAST,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KUDV ,N2ASHX,WRK,KFREE,LFREE)
         CALL DCOPY(N2ASHX,DV,1,WRK(KUDV),1)
Clf         CALL DSPTSI(NASHT,DV,WRK(KUDV))
         DO 2000 ISYM = 1,NSYM
            JSYM  = MULD2H(ISYM,JWOPSY)
            NASHI = NASH(ISYM)
            NASHJ = NASH(JSYM)
         IF (NASHI .EQ. 0 .OR. NASHJ .EQ. 0) GO TO 2000
            NORBI = NORB(ISYM)
            NBASI = NBAS(ISYM)
            NORBJ = NORB(JSYM)
            NBASJ = NBAS(JSYM)
            JCMO  = ICMO(JSYM) + 1 + NISH(JSYM)*NBASJ
            JUDV  = KUDV + IASH(ISYM)*NASHT + IASH(JSYM)
            CALL DGEMM('N','N',NBASJ,NASHI,NASHJ,1.D0,
     &                 CMO(JCMO),NBASJ,
     &                 WRK(JUDV),NASHT,0.D0,
     &                 WRK(KDAO1),NBASJ)
            JCMO  = ICMO(ISYM) + 1 + NISH(ISYM)*NBASI
            CALL DGEMM('N','T',NBASJ,NBASI,NASHI,1.D0,
     &                 WRK(KDAO1),NBASJ,
     &                 CMO(JCMO),NBASI,0.D0,
     &                 DVAO(IBAS(JSYM)+1,IBAS(ISYM)+1),NBAST)
 2000    CONTINUE
      CALL MEMREL('FCKDEN2',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      END IF
      END IF
C
C
C     ************************************
C     DCAO matrix
C
      IF (GETDC) THEN
         CALL DZERO(DCAO,N2BASX)
      IF (NISHT .GT. 0 .AND. JWOPSY .EQ. 1) THEN
         DO 4000 ISYM = 1,NSYM
            NISHI = NISH(ISYM)
         IF (NISHI .GT. 0) THEN
            NBASI = NBAS(ISYM)
            JCMO  = ICMO(ISYM) + 1
            CALL DGEMM('N','T',NBASI,NBASI,NISHI,2.D0,
     &                 CMO(JCMO),NBASI,
     &                 CMO(JCMO),NBASI,0.D0,
     &                 DCAO(IBAS(ISYM)+1,IBAS(ISYM)+1),NBAST)
         END IF
 4000    CONTINUE
C
C **     in DGEMM above: multiply with DC(i,i) = 2.0 factor to get final result
C        (if transition density matrix, then DCAO must be
C         multiplied with the <L | R> overlap outside this routine)
C
      END IF
      END IF
C
C
C
C
      IF (P6FLAG(18)) THEN
        IF (GETDC) THEN
          WRITE(LUPRI,'(/A)')
     &    ' FCKDEN2: DCAO = Density matrix, inactive part (AO-basis)'
          CALL OUTPUT(DCAO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
        END IF
        IF (GETDV) THEN
          WRITE(LUPRI,'(/A)')
     &    ' FCKDEN2: DVAO = Density matrix,   active part (AO-basis)'
          CALL OUTPUT(DVAO,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
        END IF
      END IF
C
C
      CALL QEXIT('FCKDEN2')
      RETURN
C
C *** end of subroutine FCKDEN2
C
       END
C  /* Deck Make_ds_localspin */
      SUBROUTINE MAKE_DS_LOCALSPIN(DT,DS,WORK,LWORK)
C
C FRAN & HJJ - Feb 2010 / Oct 2010
C
C PURPOSE:
C  Do the Spin-Density as DS = f(DT) * DT
C     model 1: f(DT) = (2-DT)*DT
C
#include "implicit.h"
#include "priunit.h"
#include "inforb.h"
#include "infpri.h"
      DIMENSION DT(NASHT,NASHT),DS(NASHT,NASHT), WORK(LWORK)
C
!     IF (P6FLAG(15)) THEN
         WRITE(LUPRI,'(/A)') 'DV from MAKE_DS_LOCALSPIN'
         CALL OUTPAK(DT,NASHT,1,LUPRI)
!     END IF

!     IF (model 1) then   ! only model implemented so far
         KUDT = 1
         KUDS = KUDT + N2ASHX
         KUDI = KUDS + N2ASHX ! intermediate matrix
         CALL DSPTGE(NASHT,DT,WORK(KUDT))
!>       step 1 : D_S = D_I = D_T D_T
         CALL DGEMM('N','N',NASHT,NASHT,NASHT,1.0D0,
     &      WORK(KUDT),NASHT,
     &      WORK(KUDT),NASHT, 0.0D0,
     &      WORK(KUDI),NASHT)
         CALL DCOPY(N2ASHX,WORK(KUDI),1,WORK(KUDS),1)
!>       step 2 : D_S = - D_T D_I + 2 D_S = (2-D_T) D_T D_T
         CALL DGEMM('N','N',NASHT,NASHT,NASHT,-1.0D0,
     &      WORK(KUDT),NASHT,
     &      WORK(KUDI),NASHT, 2.0D0,
     &      WORK(KUDS),NASHT)
         CALL DGETSP(NASHT,WORK(KUDS),DS)
!     END IF ! models
C
C     IF (P6FLAG(15)) THEN
         WRITE (LUPRI,1000)
         CALL OUTPUT(WORK(KUDS),1,NASHT,1,nasht,nasht,nasht,1,LUPRI)
C     END IF
 1000 FORMAT(/' Spin. den. mat.: DS = (2-Dtot)Dtot Dtot , MO-basis')
C
C END OF MAKE_DS_LOCALSPIN
C
      RETURN
      END
