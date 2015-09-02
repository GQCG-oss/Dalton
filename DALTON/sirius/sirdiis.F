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
C
C
C FILE : sirius/sirdiis.F
C Principal author : Hans Jørgen Aa. Jensen
C First version    : 1993
C
#ifdef OLD_REV_LOG
!===========================================================================
!se efter HJTODO
!HJTODO 940816: if .FREEZE implemented in FCKEIG (requires new ORDRSS)
!  then check all setting of MXDIIS = 0 and MAXFCK = 0 (I think MAXFCK
!  also will work if .FREEZE implemented in FCKEIG).
!
!Revision 1.9  2001/02/01 17:30:12  hjj
!restart DIIS if DIIS is stalled ...
!
!Revision 1.6  2000/05/24 12:28:27  hjj
!1) fixed solvent calculations (DIFDEN must be false for solvent)
!2) fixed WR_SIRIFC call for open shell and solvent
!   (required change in SOLFCK also)
!=============================================================================
!970326-hjaaj -- adjusted abort DIIS algorithm (allow more NEWOCC)
!951130-hjaaj
!DIIS_CTL: moved 'reset to not converged' for solvent and writing SIRIFC to
!        DIIS_CTL from SIRCTL; it is a deficiency in the DIIS routines that
!        it cannot be done in them. Look for MAERKE.
!950824-hjaaj
!DIIS_CTL: output for AUTOCC
!FCKEIG: also AUTOCC for open shell plus smaller changes.
!9505-k.ruud implemented AUTOCC option
!940701-hjaaj
!SRDIIS: let MXERRV to default to MAX(2,MIN(NWOPT,10)) (was MXDIIS).
!  This fixes linear dependency problem for very small debug runs:
!  MXERRV=10 gave lin.dep. for 2 el. in 2 orbitals,
!  MXERRV=2 gave quadratic convergence !
!940511-hjaaj
!SRDIIS,DIIS_CTL: call flshfo so progress can be followed in output file;
!   moved heading printing from DIIS_CTL to SRDIIS
!940505-hjaaj
!DIIS_CTL: exit if energy increasing or convergence is too slow.
!940501-hjaaj
!FCKEIG: exit if ".FREEZE" not compatible with ORDRSS.
!    ORDRSS now always called (this also fixes error in prev. version).
!940408-hjaaj
!DIIS_CTL: s/EMCACT/EACTIV/ to define active energy in /INFOPT/
!940311-hjaaj
!DIIS_CTL: s/EVCNRM/GRDNRM/ (so gradient norm is transferred to SIROUT)
!        nice output to LUW4 if LUW4 .ne. LUPRI
!931126-hjaaj
!DIIS_CTL: call NEWORB after each iteration (for restart)
!930720-hjaaj
!SRDIIS: new parameter PROCC .FALSE. in CALL PRORB
!930715-hjaaj
!------: improved some print
!SRDIIS: exit MXDGSS .GT. 0
!930623-hjaaj
!FCKEIG: inserted PARAMETER(D0=0.0D0)
!===========================================================================
#endif
C  /* Deck srdiis */
      SUBROUTINE SRDIIS (ICONV,CMO,WRK,KFRSAV,LFRSAV)
C
C Copyright (c) 1993, 1997 Hans Joergen Aa. Jensen
C
C Driver for DIIS iterations in SIRIUS.
C

#include "implicit.h"
      DIMENSION CMO(*), WRK(*)
C
C
C Used from common blocks:
C  INFINP : DIRFCK, FLAG(*)
C  SCBRHF : MAXEVC
C  INFORB : N2BAST,NNBAST,N2BASX,NCMOT,NSYM, ...,MXDGSS
C  INFVAR : NWOPT
C  INFPRI : IPRRHF
C
#include "maxorb.h"
#include "priunit.h"
#include "infinp.h"
#include "scbrhf.h"
#include "inforb.h"
#include "infvar.h"
#include "infpri.h"
#include "gnrinf.h"
C
Cholesky
#include "ccdeco.h"
#include "choscf.h"
C
      PARAMETER (NBLOC = 1)
      INTEGER LUBLOC(NBLOC)
      DATA LUBLOC / 21 /
Cholesky
C
      PARAMETER (D4 = 4.0D0)
      LOGICAL C2DIIS
C
      CALL QENTER('SRDIIS')
      ICONV = 0
      IF (MXDGSS .GT. 1) THEN
         NINFO = NINFO + 1
         WRITE (LUPRI,'(//A/)')
     &   ' INFO: no DIIS iterations because degenerate'//
     &   'supersymmetries is not implemented for DIIS'
         GO TO 9999
      END IF
C
      KFREE  = KFRSAV
      LFREE  = LFRSAV

C
CHJMAERKE: define IPR_DIIS separately?
C
      MXERRV = MAXEVC
      IF (MXERRV .LE. 0) THEN
         MXERRV = MIN(NWOPT,8)
C        ... fixes linear dependency problem for very small debug runs
C            MXERRV=10 gave lin.dep. for 2 el. in 2 orbitals,
C            MXERRV=2 gave quadratic convergence ! (940701-hjaaj)
C        (130919-hjaaj: decreased default from 10 to 5, because
C         with 10 the span of BMAT may be >1.D15, giving problems
C         in the DIIS solver - an example has been seen !!!)
C        (140307-jmho: increased default from 5 to 8, because
C         with 5 it often resulted in 1-2 extra SCF iterations. At 8
C         we should still avoid the problem mentioned above.)
         MXERRV = MAX(2,MXERRV)
      END IF
      IPR_DIIS = IPRRHF
      C2DIIS = FLAG(22)
C
      LUWOUT = LUPRI
  10  WRITE (LUWOUT,'(//A/A/A)')
     &      ' ***********************************************',
     &      ' ***** DIIS acceleration of SCF iterations *****',
     &      ' ***********************************************'
Cholesky
      IF (CHOINT) THEN
         WRITE(LUWOUT,'(A)')
     &      '   -- using Cholesky decomposed integrals'
         IF (.NOT. CCMODSK) THEN
            WRITE(LUWOUT,'(A,1P,D9.2,A)')
     &      '   -- using Cholesky decomposed density, thr = ',THRDC
         END IF
      END IF
Cholesky
      IF (C2DIIS) THEN
         WRITE (LUWOUT,'(/A,I5)')
     &      ' C2-DIIS algorithm; max error vectors =',MXERRV
      ELSE
         WRITE (LUWOUT,'(/A,I5)')
     &      ' C1-DIIS algorithm; max error vectors =',MXERRV
      END IF
      CALL FLSHFO(LUWOUT)
      IF (LUWOUT .NE. LUW4) THEN
         LUWOUT = LUW4
         GO TO 10
      END IF
C
C     Get AO overlap matrix; transform left index to CMO1 basis.
C
      CALL MEMGET2('REAL','SMOAO',KSMOAO,N2BAST,WRK,KFREE,LFREE)
      KREL = KFREE
      CALL MEMGET2('REAL','SAO' ,KSAO ,NNBAST,WRK,KFREE,LFREE)
      CALL MEMGET2('REAL','STMP',KSTMP,N2BAST,WRK,KFREE,LFREE)
      CALL RDONEL('OVERLAP',.TRUE.,WRK(KSAO),NNBAST)
      IF (IPR_DIIS .GT. 25) THEN
         WRITE (LUPRI,'(/A)') ' SRDIIS test output of SAO:'
         CALL OUTPKB(WRK(KSAO),NBAS,NSYM,-1,LUPRI)
         WRITE (LUPRI,'(/A)') ' SRDIIS test output of CMO1:'
         CALL PRORB(CMO,.FALSE.,LUPRI)
C        CALL PRORB(CMO,PROCC,IOUT)
      END IF
      DO 100 I = 1,NSYM
         NORBI = NORB(I)
         IF (NORBI .GT. 0) THEN
            NBASI = NBAS(I)
            CALL DSPTSI(NBASI,WRK(KSAO+IIBAS(I)),WRK(KSTMP))
            CALL DGEMM('T','N',NORBI,NBASI,NBASI,1.D0,
     &                 CMO(1+ICMO(I)),NBASI,
     &                 WRK(KSTMP),NBASI,0.D0,
     &                 WRK(KSMOAO+I2BAS(I)),NORBI)
         END IF
  100 CONTINUE
C     scale by 4 to get correct orbital gradient from DIIS_EVC
      CALL DSCAL(N2BAST,D4,WRK(KSMOAO),1)
      CALL MEMREL('SRDIIS.SMOAO',WRK,KFRSAV,KREL,KFREE,LFREE)
C
      IF (NASHT .EQ. 0) THEN
         NDMAT = 1
      ELSE
         NDMAT = 2
      END IF
      CALL MEMGET2('REAL','CMO1',KCMO1 ,NCMOT, WRK,KFREE,LFREE)
      CALL MEMGET2('REAL','DCAO',KDCAO ,NDMAT*N2BASX,WRK,KFREE,LFREE)
      CALL MEMGET2('REAL','FCAO',KFCAO ,NDMAT*N2BASX,WRK,KFREE,LFREE)
      CALL MEMGET2('REAL','H1AO',KH1AO ,NNBAST,WRK,KFREE,LFREE)
      CALL MEMGET2('REAL','DSAV',KDSAV ,NDMAT*NNBAST,WRK,KFREE,LFREE)
      IF (FLAG(16)) THEN
         CALL MEMGET2('REAL','TNLM',KTNLM,2*NLMSOL,WRK,KFREE,LFREE)
      ELSE
         CALL MEMGET2('REAL','TNLM',KTNLM,0,WRK,KFREE,LFREE)
      END IF
C     estimate max value for MXERRV (951130-hjaaj) :
      IF (DIRFCK) THEN
C        -- direct SCF:
C          NDMAT*N2BASX used in HERFCK(hermit)
C                N2BASX used in SKLFC1(hermit)
C          1000000      (1 mw estimate for TWOINT(hermit) and the rest)
         NEVC  = (LFREE - (NDMAT+1)*N2BASX - 1000000) / NNBAST
      ELSE
C        -- conventional SCF:
C                N2BASX used in FCKDEN
C       somewhat arbitrary estimate of memory needed in RDSUPM etc.
C       estimate quite big to be on the safe side /960627-hjaaj
         NEVC  = (LFREE - N2BASX - 100000) / NNBAST
      END IF
      NEVC  = NEVC / (NDMAT + 1)
      IF (NEVC .LT. MXERRV) THEN
         IF (NEVC .LT. 3) THEN
            NWARN = NWARN + 1
            WRITE (LUPRI,'(/A/A,I5)')
     &      ' SRDIIS WARNING: insufficient memory for DIIS',
     &      '                 max # of error vectors only',NEVC
            GO TO 9990
         END IF
         NINFO = NINFO + 1
         WRITE (LUPRI,'(/A,I5,A,I5/A)')
     &   ' SRDIIS INFO: MXERRV reduced from',MXERRV,' to',NEVC,
     &   ' SRDIIS INFO: because of memory limitations.'
         MXERRV = NEVC
      END IF
C
      LBMAT = MXERRV*MXERRV
      ! was LBMAT = (MXERRV+1)*MXERRV / 2
      ! but in opt-solvers the square matrix is used
      ! and memory requirement is small so no case testing here ...
      CALL MEMGET2('REAL','BMAT',KBMAT ,LBMAT, WRK,KFREE,LFREE)
      CALL MEMGET2('REAL','CVEC',KCVEC ,MXERRV+1,WRK,KFREE,LFREE)
      CALL MEMGET2('REAL','GVEC',KGVEC ,MXERRV+1,WRK,KFREE,LFREE)
      CALL MEMGET2('REAL','FSAV',KFSAV ,NDMAT*MXERRV*NNBAST,
     &   WRK,KFREE,LFREE)
      CALL MEMGET2('REAL','EVSAV',KEVSAV,MXERRV*NNORBT,WRK,KFREE,LFREE)
      CALL MEMGET2('REAL','EGSAV',KEGSAV,MXERRV*2     ,WRK,KFREE,LFREE)
      CALL MEMGET2('REAL','DV',KDV,NNASHX,WRK,KFREE,LFREE)
      CALL DCOPY(NCMOT,CMO,1,WRK(KCMO1),1)
      CALL DIIS_CTL(CMO,WRK(KCMO1),WRK(KGVEC),WRK(KBMAT),WRK(KH1AO),
     &            WRK(KFSAV),WRK(KEVSAV),WRK(KEGSAV),MXERRV,C2DIIS,
     &            WRK(KSMOAO),WRK(KDCAO),WRK(KFCAO),WRK(KDSAV),WRK(KDV),
     &            WRK(KCVEC),WRK(KTNLM),
     &            WRK,KFREE,LFREE, IPR_DIIS,ICONV)
C
 9990 CALL MEMREL('SRDIIS',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
 9999 CALL FLSHFO(LUW4)
      IF (LUPRI .NE. LUW4) CALL FLSHFO(LUPRI)
      CALL QEXIT('SRDIIS')
      RETURN
      END
C  /* Deck DIIS_CTL */
      SUBROUTINE DIIS_CTL(CMO,CMO1,GVEC,BMAT,H1AO,FAOSAV,EVCSAV,EGSAV,
     *                    MXERRV,C2DIIS,SMOAO,DCAO,FCAO,DAOSAV,DV,
     *                    CVEC,TNLM,WRK,KFRSAV,LFRSAV,
     *                    IPR_DIIS,ICONV)

#ifdef ENABLE_EFS
      use efs_interface
#endif
C
C L.r. 5-May-1994 hjaaj / 1-Apr-1997
C
C DFT modifications T.  Helgaker
C
C
      use qmcmm, only: getdim_relmat, comp_relmat, write_relmat,
     &                 getdim_mmmat, comp_mmrelmat
#ifdef ENABLE_XCINT
      use dalton_xcint_interface, only: xcint_integrate_rks_scf,
     &                                  xcint_uses_own_grid,
     &                                  xcint_is_active
#endif
      use qmcmm_fock, only: qmnpmm_fock
#if defined (HAS_PCMSOLVER)
      use pcm_scf
      use pcmmod_cfg
#endif
      use pelib_interface, only: use_pelib, pelib_ifc_fock,
     &                           pelib_ifc_energy
#include "implicit.h"
#include "dummy.h"
      DIMENSION CMO(*), CMO1(*), GVEC(*), BMAT(*), H1AO(*)
      DIMENSION FAOSAV(NNBAST,MXERRV,2), EVCSAV(NNORBT,MXERRV)
      DIMENSION SMOAO(*), DCAO(*), FCAO(*), DAOSAV(NNBAST,*), CVEC(*)
      DIMENSION DV(*)
      DIMENSION EGSAV(2,MXERRV), TNLM(*), WRK(*)

#ifdef ENABLE_XCINT
      real(8)              :: xc_energy(1)
      real(8), allocatable :: dmat_ao(:)
      real(8), allocatable :: xc_mat(:)
      real(8), allocatable :: buffer1(:)
      real(8), allocatable :: buffer2(:)
      real(8), allocatable :: xc_mat_tri(:)
      logical              :: file_exists
#endif

#if defined (HAS_PCMSOLVER)
      real(8), allocatable :: density_ao(:), density_ao_symm(:)
      real(8), allocatable :: fock_pcm_ao(:), fock_pcm_ao_symm(:)
#endif

C
C Used from common blocks:
C  INFINP : FLAG(*), HSROHF, DODFT; ?
C  INFORB : NCMOT,NNBAST,NSYM
C  INFOPT : EPOT,EMCSCF,EMY,EACTIV,GRDNRM
C  SCBRHF : THRRHF,MXDIIS, ?
C  DFTCOM : HFXFAC, DODFTD
C  INFVAR : NWOPH
C  SCBRHF : IOPRHF,AUTOCC
C  CBIREA : LMULBS
C  R12INT : R12CAL
C
#include "maxorb.h"
#include "priunit.h"
#include "mxcent.h"

#include "molde.h"
#include "infinp.h"
#include "pcmlog.h"
#include "inforb.h"
#include "infvar.h"
#include "infopt.h"
#include "scbrhf.h"
#include "dftcom.h"
#include "cbihr2.h"
#include "infpri.h"
#include "gnrinf.h"
#include "dfterg.h"
#include "cbirea.h"
#include "r12int.h"
#include "qm3.h"
#include "incore.h"
#include "qmmm.h"
#include "infloc.h"
#include "mmtimes.h"
#include "infpar.h"
C
Cholesky
#include "ccdeco.h"
#include "sirftim.h"
#include "choscf.h"
C
#include "center.h"
C
Cholesky
C
C------------------------
#include "pcmdef.h"
#include "pcm.h"

! QM/NP/MM code
#include "qmnpmm.h"

      PARAMETER (THDDEF = 1.0D-5)
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0, D1 = 1.0D0, DM1 = -1.0D0)
      PARAMETER (THREVC = 0.01D0, THRINC = 1.0D-10, CNVFAC = 0.8D0)
C
      LOGICAL   ONLYFC, C2DIIS, DIFDEN, DODIFDEN, NEWTHR, RHFWOP, OLDWOP
      LOGICAL   LCONV, LNEWOCC, HSROHF_save
      DIMENSION NISH_save(8), NCAN(8)
      CHARACTER*8 CMOLBL

      REAL(8) :: PE_ENERGY
      REAL(8), DIMENSION(:), ALLOCATABLE :: PE_DENMAT, PE_FCKMAT

      real(8), allocatable :: relay_matrix(:, :, :)
      real(8), allocatable :: mm_relay_matrix(:)
      real(8), allocatable :: qmcmm_sol(:)

!#define EXPERIMENTAL_DMAT_RESTART
#ifdef EXPERIMENTAL_DMAT_RESTART
      ! TODO: do not allocate extra matrix, not needed
      !       but read directly into DCAO
      !       then DSCAL by 2.0 to match Dalton
      real(8), allocatable :: dmat_restart(:)
      logical              :: is_gcbasis
#endif
C
      CALL QENTER('DIIS_CTL')
      KFREE  = KFRSAV
      LFREE  = LFRSAV
      IFTHSV = IFTHRS
C
C     save current IOPRHF for check below
      IOPRHF_old = IOPRHF
      IF (NSYM .EQ. 1) AUTOCC = .FALSE.

!        Sep 2013 hjaaj: now HSROHF code here also for NASHT.eq.1 because we
!        have some evidence that the open shell Fock matrix from HSROHF
!        gives better convergence than the one from the older "single open shell".
!        At return we restore HSROHF to .false. for this case, because some
!        ABACUS and other modules require the old "single open shell"
!        to work correctly.

      HSROHF_save = HSROHF
      HSROHF = HSROHF .OR. (NASHT .EQ. 1 .AND. .NOT. AUTOCC)

C
C     Get one-electron Hamiltonian
C
      CALL SIRH1(H1AO,WRK(KFREE),LFREE)

CAMT  If requested calculate the DFT-D dispersion energy correction
      IF (DODFT .AND. DODFTD) THEN
         NDERIV=0
         CALL DFT_D_DAL_IFC(EDISP,NDERIV,WRK(KFREE),LFREE)
      ELSE
         EDISP = 0.0D0
      END IF

C
      IF (LUW4 .NE. LUPRI .OR. IPR_DIIS .LE. 1) THEN
#if defined (HAS_PCMSOLVER)
      IF (FLAG(16) .OR. PCM .OR. QM3 .OR. QMMM .OR. QMNPMM .OR.
     &    pcmmod_is_pcm_calculation) THEN
#else              
      IF (FLAG(16) .OR. PCM .OR. QM3 .OR. QMMM .OR. QMNPMM) THEN
#endif
C     this also includes the case of MMPCM   
      IF (.NOT.AUTOCC) THEN
         WRITE (LUW4,'(/A)')
     &      ' Iter      Total energy       Solvation energy    '//
     &      'Error norm    Delta(E)'
      ELSE IF (IOPRHF .EQ. 0) THEN
         WRITE (LUW4,'(/A,I4,A//A)')
     &      ' Automatic occupation of symmetries with',NRHFEL,
     &      ' electrons.',
     &      ' Iter     Total energy      Solvation energy  Error norm'//
     &      '  Delta(E)    SCF occupation'
      ELSE
         WRITE (LUW4,'(/A,I4,A/A//A)')
     &      ' Automatic occupation of symmetries with',NRHFEL,
     &      ' electrons.',
     &      ' "sym." is the symmetry of the one open shell orbital '//
     &      'in table below.',
     &      ' Iter     Total energy      Solvation energy  Error norm'//
     &      '  Delta(E)  Sym.  Closed shell occupation'
      END IF
      ELSE IF (USE_PELIB()) THEN
        WRITE(LUW4,'(/A)')
     &      ' Iter      Total energy     Embedding energy      '//
     &      'Error norm     Delta(E)'
      ELSE ! no solvation energy
      IF (.NOT.AUTOCC) THEN
         WRITE (LUW4,'(/A)')
     &      ' Iter      Total energy        Error norm    Delta(E)'//
     &      '  DIIS dim.'
      ELSE IF (IOPRHF .EQ. 0) THEN
         WRITE (LUW4,'(/A,I4,A//A)')
     &      ' Automatic occupation of symmetries with',NRHFEL,
     &      ' electrons.',
     &      ' Iter     Total energy    Error norm  Delta(E)  '//
     &      '  SCF occupation'
      ELSE
         WRITE (LUW4,'(/A,I4,A/A//A)')
     &      ' Automatic occupation of symmetries with',NRHFEL,
     &      ' electrons.',
     &      ' "sym." is the symmetry of the one open shell orbital '//
     &      'in table below.',
     &      ' Iter     Total energy    Error norm  Delta(E)  '//
     &      'Sym.  Closed shell occupation'
      END IF
      END IF
      END IF
C
      IF (NISHT.EQ.0) CALL DZERO(DAOSAV,NNBAST)
      ONLYFC = (NASHT .EQ. 0)
      IF (ONLYFC) THEN
C        DV(1) = D0
         NDMAT = 1
      ELSE IF (NASHT .GT. 0) THEN
         CALL DZERO(DV,NNASHX)
         DO I = 1, NASHT
            II = I*(I+1)/2
            DV(II) = D1
         END DO
         NDMAT = 2
      ELSE
         CALL QUIT('DIIS_CTL called with NASHT .gt. 1')
      END IF
      THDIIS = THRRHF
      IF (THDIIS .LE. D0) THDIIS = THDDEF
C
C     Level shift SHFTLVL ("restricted step" type step reduction)
C     - initial value: see next code block
C     - used to modify FD matrix: FDAO(shifted) = FDAO + SHFTLVL*DCAO
C       - this corresponds to increasing occ-unocc orb.en. gap by SHFTLVL
C         and approximately to add abs(SHFTLVL) to the approximate
C         Hessian diagonal implicit in the Roothaan diagonalization.
C     -- Feb. 2001, Hans Joergen Aa. Jensen.
C
c     DEFLVL is read from SCF INPUT/.SHIFT
      IF (DEFLVL .gt. 0.9D0) THEN ! .SHIFT not set in input
         SHFTLVL = -0.20D0
         DEFLVL  =  0.0D0 ! no default level shift in e.g. next geometry
      ELSE ! .SHIFT has ben set in input
         SHFTLVL = DEFLVL
      END IF
      SHFTFAC = 0.5D0

C ---  hjaaj: find safe screening factor depending on
C             convergence threshold:
C     SCRFAC = 100*NBAST
      SCRFAC = 100*N2BAST
      IFTMAX = INT(-LOG10(THDIIS/SCRFAC)) + 1
      IFTMIN = INT(-LOG10( 0.1D0/SCRFAC)) + 1
Chjaaj: screening factor for safe Fock matrix construction
C       is chosen to be 100 (we want accuracy to 0.1 times THDIIS)
C       F_ij ~ sum(kl) (ij|kl) * D_kl, where sum(kl) is over
C       max N2BAST elements (for totally symmetric "kl");
C       statistically we expect error to be ca. sqrt(n2bast) < nbast
C       The 100*NBAST is thus with an extra factor 10 factor for
C       safety, and because IFTHRS only can change screening
C       with a factor 10.   /01-Feb-2001 hjaaj
C
C       Yes, but for the gradient norm we in addition will have a sum
C       of ca. NOCCT*NVIRT/NSYM elements, thus to be safe wrt
C       gradient norm let us use 10*N2BAST in SCRFAC instead
C       of 100*NBAST.
C       Also minimum accuracy corresponding to error in
C       GRDNORM ~ 0.1 (IFTMIN). /31-Jul-2006 hjaaj
C
C       In iterations we want 0.001 times GRDNRM because we do
C       keep vectors until GRDNRM lowered with a factor 100
      DAMP   = D0
      EMCLOW = D0
      GRDNRM = -1.0D0
      ITNOCC = 0
      NEWOCC = 0
      LNEWOCC = AUTOCC
C     if (LNEWOCC) then write new BASINFO to SIRIUS.RST in NEWORB
      I_STALL = 0
C------------------
      DODIFDEN = .NOT.FLAG(49) .AND. .NOT.HSROHF
     &     .AND. .NOT.DODFT    .AND. .NOT.DOHFSRDFT
     &     .AND. .NOT.AOSAVE   .AND. .NOT.CHOINT
     &     .AND. .NOT.EMBEDDING
#ifdef ENABLE_EFS
     &     .AND. .NOT. use_efs
#endif
C
#if defined (HAS_PCMSOLVER)
      if (pcmmod_is_pcm_calculation) then
              dodifden = .false.
      end if
#endif
C------------------
Chjaaj: DIFDEN does not work with solvent/PCM/QM3
C       because FC_sol is added to the old FC matrix, but not the new!
Ctuh:   DIFDEN does not work with DFT either
Cov:    DIFDEN disabled for high spin HF
C-tbp: turn off DIFDEN for Cholesky decomposition (CHOINT true)
      DIFDEN   = DODIFDEN
C------------------
      IF (.NOT. USRSCR) THEN
C        AOSAVE or Cholesky (CHOINT): No screening.
         IF (AOSAVE .OR. CHOINT .OR. .NOT. DIRFCK) THEN
            IFTHRS = 20 ! code for no screening
         ELSE IF (FLAG(13)) THEN
Chj-aug99: MO's read from input, assume good guess
            IFTHRS = IFTMAX
         ELSE
Chj-aug99: MO's from huckel or H1DIAG, assume poor guess
            IFTHRS = IFTMIN
         END IF
      END IF
      IF (ICI0 .EQ. 6) THEN
Chj      ... if (GEOWLK) then
         EMCGEO = EMCOLD
         DEPGEO = DEPRED
      END IF

C
      ITDIIS = 0
      JTDIIS = 0
C
C Projection operators for removing non-rotated orbitals (.FREEZE option)
C A -> P_f*A*P_f(T), P_f = 1-S*D_f
C
      IF (LNOROT) THEN
         CALL MEMGET('REAL',KPFROZ,N2BASX,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KDFROZ,N2BASX,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KS,N2BASX,WRK,KFREE,LFREE)
         CALL GET_H1(WRK(KS),'OVERLAP',WRK(KFREE),LFREE)
         CALL DUNIT(WRK(KPFROZ),NBAST)
         CALL DZERO(WRK(KDFROZ),N2BASX)
         KDFI=KDFROZ
         DO ISYM=1,NSYM
            IORBI=IORB(ISYM)
            NORBI=NORB(ISYM)
            DO I=1,NORBI
               IF (NOROT(IORBI+I).NE.0) THEN
                  NBASI=NBAS(ISYM)
                  ICMOI=ICMO(ISYM)+NBASI*(I-1)+1
                  ! C*C.T()
                  CALL DGEMM(
     &               'N','T',NBASI,NBASI,1,
     &               D1,CMO(ICMOI),NBASI,
     &                  CMO(ICMOI),NBASI,
     &               D1,WRK(KDFI),NBAST
     &               )
               END IF
            END DO
            KDFI=KDFI+NBAST*NBAS(ISYM)+NBAS(ISYM)
         END DO
         CALL DGEMM('N','N',NBAST,NBAST,NBAST,
     &      DM1,WRK(KS),NBAST,
     &          WRK(KDFROZ),NBAST,
     &      D1, WRK(KPFROZ),NBAST
     &      )
         IF (IPR_DIIS.GT.20) THEN
            CALL HEADER('Start orbitals',-1)
            CALL PRORB(CMO,.FALSE.,LUPRI)
            CALL HEADER('Frozen orbital density',-1)
            CALL OUTPUT(
     &         WRK(KDFROZ),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI
     &         )
            CALL HEADER('Frozen orbital projector',-1)
            CALL OUTPUT(
     &         WRK(KPFROZ),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI
     &         )
         END IF
      END IF

      ITDIIS = 0
      JTDIIS = 0

  100 CONTINUE
         WRITE(LUW4,'(A)')
     &      ' -----------------------------------------------------'//
     &      '------------------------'
         ITDIIS = ITDIIS + 1
         JTDIIS = JTDIIS + 1
         ITNOCC = ITNOCC + 1
         INDEVC = MOD(JTDIIS-1,MXERRV) + 1
         NEVC   = MIN(MXERRV,JTDIIS)
C
C
         IF (IPR_DIIS .GT. 1) WRITE (LUPRI,'(/A,I5)')
     &         ' *** DIIS iteration number',ITDIIS
         CALL FLSHFO(LUPRI)
C
C  ***** Construct folded inactive and active parts of the one-electron
C        density matrix in AO-basis (DCAO and DVAO). DV is the active
C        part of the one-electron density matrix in MO-basis.
C
         CALL FCKDEN((NISHT.GT.0.OR.HSROHF),(.NOT.ONLYFC),
     &                DCAO,DCAO(1+N2BASX),CMO,DV,WRK(KFREE),LFREE)
C        CALL FCKDEN(GETDC,GETDV,DCAO,DVAO,CMO,DV,WRK,LWRK)

#ifdef EXPERIMENTAL_DMAT_RESTART
! radovan: experimental and untested feature, use at own risk
!          reads AO dmat from LSDalton's dens.restart

         if (itdiis == 1) then
            iunit = -1
            CALL GPOPEN(iunit,'dens.restart','UNKNOWN',' ',
     &                 'UNFORMATTED',IDUMMY,
     &                 .FALSE.)
            nrow = 0
            ncol = 0
            read(iunit) nrow, ncol

            ! check that dimensions match
            if (nrow*ncol /= n2basx) then
               call quit('dimensions do not match in restart')
            end if

            ! check that there is no symmetry
            if (nsym > 1) then
               call quit('nsym > 1 in restart')
            end if

            allocate(dmat_restart(nrow*ncol))
            read(iunit) (dmat_restart(I),I=1,nrow*ncol)
            read(iunit) is_gcbasis

            ! check that gcbasis is .false.
            if (is_gcbasis == .true.) then
               call quit('is_gcbasis is true in restart')
            end if

            CALL GPCLOSE(iunit, 'KEEP')
            dmat_restart = 2.0*dmat_restart
            CALL DCOPY(nrow*ncol,dmat_restart,1,DCAO,1)
            deallocate(dmat_restart)
         end if
#endif /* EXPERIMENTAL_DMAT_RESTART */

C
C     For a restricted open shell DFT calculation, or high spin
C     Let the matrices be Fc=Fi+Fa and Fc-Fo=Fa-Q=-Fa(exch)
C     generated by total and spin (=active) densities
C     This choice avoids the calculation of a third Q matrix
C
         IF (HSROHF) THEN
            CALL DAXPY(N2BASX,D1,DCAO(1+N2BASX),1,DCAO,1)
            CALL DSCAL(N2BASX,-D1,DCAO(1+N2BASX),1)
         END IF
C
C-tbp: For Cholesky, save a copy of the occupied part of CMO
C      on disk (for use in SIRFCK)
C
          IF (CHOINT .AND. CCMODSK) THEN
             CALL WOPEN2(LUCCMO,FCCMO,64,0)
             IADR2 = 1
             DO ISYM = 1,NSYM
                NDVCS(ISYM) = NISH(ISYM)
                IF (NISH(ISYM) .GT. 0) THEN
                   JCMO = ICMO(ISYM) + 1
                   LEN = NBAS(ISYM)*NISH(ISYM)
                   CALL PUTWA2(LUCCMO,FCCMO,CMO(JCMO),IADR2,LEN)
                   IADR2 = IADR2 + LEN
                END IF
             ENDDO
             CALL WCLOSE2(LUCCMO,FCCMO,'KEEP')
         END IF
C
         CALL MEMGET2('REAL','TMP1',KTMP1,NNBASX,WRK,KFREE,LFREE)
         CALL MEMGET2('REAL','TMP2',KTMP2,N2BASX,WRK,KFREE,LFREE)
         IF (NISHT.GT.0 .OR. HSROHF) THEN
            CALL DGEFSP(NBAST,DCAO,WRK(KTMP1))
            IF (.NOT.DIFDEN .OR. NEVC .EQ. 1) THEN
               CALL PKSYM1(WRK(KTMP1),DAOSAV,NBAS,NSYM,1)
C              ... save DCAO in DAOSAV
               DCOVLP = D0
            ELSE
               CALL PKSYM1(WRK(KTMP1),WRK(KTMP2),NBAS,NSYM,1)
               DCOVLP = DDOT(NNBAST,DAOSAV,1,DAOSAV,1)
               DCOVLP = DDOT(NNBAST,WRK(KTMP2),1,DAOSAV,1) / DCOVLP
               CALL PKSYM1(WRK(KTMP1),DAOSAV,NBAS,NSYM,-1)
               CALL DCOPY(NNBAST,WRK(KTMP2),1,DAOSAV,1)
C              ... save new DCAO in DAOSAV
               CALL DUNFLD(NBAST,WRK(KTMP1),WRK(KTMP2))
               CALL DAXPY(N2BASX,(-DCOVLP),WRK(KTMP2),1,DCAO,1)
C              ... difden: DCAO = DCAO - DCOVLP*DAOSAV
            END IF
         END IF
         IF (.NOT.ONLYFC) THEN
            CALL DGEFSP(NBAST,DCAO(1+N2BASX),WRK(KTMP1))
            IF (.NOT.DIFDEN .OR. NEVC .EQ. 1) THEN
               CALL PKSYM1(WRK(KTMP1),DAOSAV(1,2),NBAS,NSYM,1)
C              ... save "DVAO" in DAOSAV(,2)
               DVOVLP = D0
            ELSE
               CALL PKSYM1(WRK(KTMP1),WRK(KTMP2),NBAS,NSYM,1)
               DVOVLP = DDOT(NNBAST,DAOSAV(1,2),1,DAOSAV(1,2),1)
               DVOVLP = DDOT(NNBAST,WRK(KTMP2),1,DAOSAV(1,2),1)
     &                / DVOVLP
               CALL PKSYM1(WRK(KTMP1),DAOSAV(1,2),NBAS,NSYM,-1)
               CALL DCOPY(NNBAST,WRK(KTMP2),1,DAOSAV(1,2),1)
C              ... save new "DVAO" in DAOSAV(,2)
               CALL DUNFLD(NBAST,WRK(KTMP1),WRK(KTMP2))
               CALL DAXPY(N2BASX,(-DVOVLP),WRK(KTMP2),1,
     &            DCAO(1+N2BASX),1)
C              ... difden: "DVAO" = "DVAO" - DVOVLP*DVOSAV(,2)
            END IF
         END IF
         CALL MEMREL('DIIS_CTL.DIFDEN',WRK,1,KTMP1,KFREE,LFREE)
C
C  ***** Calculate 2-electron part of Fock matrix with differential density
C

         CALL GETTIM(t1_fck,w1_fck)

#if defined ENABLE_EFS
         if (use_efs) then
            call mpixbcast (const_fock_update, 1, 'INTEGER', 0)
            iprint = 0
            call mpixbcast (iprint, 1, 'INTEGER', 0)
            call efs_add_2efock (FCAO, DCAO, IFTHRS)
         else
#endif
         CALL FCK2AO(ONLYFC,FCAO,DCAO,WRK,KFREE,LFREE)
#if defined ENABLE_EFS
         end if
#endif

         IF (DIRFCK .AND. IFTHRS.LT.20) THEN
            IF (IPR_DIIS .GE. 0) THEN
               CALL GETTIM(t2_fck,w2_fck)
               WRITE (LUPRI,'(I4,A,2I5,L5,1P,2D10.2)')
     &  ITDIIS,'  Screening settings (-IFTHRS, JTDIIS, DIFDEN, times)',
     &  -IFTHRS, JTDIIS,DIFDEN,T2_FCK-T1_FCK,W2_FCK-W1_FCK
            END IF
         END IF

C
C        FAOSAV(1,INDEVC,1) = FD2AO = FC2AO + FV2AO
C        FAOSAV(1,INDEVC,2) = FV2AO
C
         IF (NSYM .GT. 1) THEN
            CALL MEMGET2('REAL','TMP1',KTMP1,NNBASX,WRK,KFREE,LFREE)
            CALL DGETSP(NBAST,FCAO,WRK(KTMP1))
            CALL PKSYM1(WRK(KTMP1),FAOSAV(1,INDEVC,1),NBAS,NSYM,1)
            IF (.NOT.ONLYFC) THEN
               CALL DGETSP(NBAST,FCAO(1+N2BASX),WRK(KTMP1))
               CALL PKSYM1(WRK(KTMP1),FAOSAV(1,INDEVC,2),NBAS,NSYM,1)
            END IF
            CALL MEMREL('DIIS_CTL.FAOSAV',WRK,1,KTMP1,KFREE,LFREE)
         ELSE
            CALL DGETSP(NBAST,FCAO,FAOSAV(1,INDEVC,1))
            IF (.NOT.ONLYFC) THEN
               CALL DGETSP(NBAST,FCAO(1+N2BASX),FAOSAV(1,INDEVC,2))
            END IF
         END IF

         IF (DIFDEN .AND. NEVC .GT. 1) THEN
            INDEVP = MOD(JTDIIS-2,MXERRV) + 1
            CALL DAXPY(NNBAST,DCOVLP,FAOSAV(1,INDEVP,1),1,
     &                               FAOSAV(1,INDEVC,1),1)
C           .. add DCOVLP*(FCAO_previous + FVAO_previous)
            IF (.NOT.ONLYFC) THEN
               IF (.NOT.HSROHF)
     &         CALL DAXPY(NNBAST,-DCOVLP,FAOSAV(1,INDEVP,2),1,
     &                                   FAOSAV(1,INDEVC,1),1)
C           .. subtract DCOVLP*FVAO_previous to get DCOVLP*FCAO_previous
               CALL DAXPY(NNBAST, DVOVLP,FAOSAV(1,INDEVP,2),1,
     &                                   FAOSAV(1,INDEVC,2),1)
            END IF
         END IF
         EMY = DDOT(NNBAST,DAOSAV,1,H1AO,1)
     &       + DDOT(NNBAST,DAOSAV,1,FAOSAV(1,INDEVC,1),1) * DP5
         IF (DOHFSRDFT) EMY = EMY + EDFTY + ESRDFTY
C
         IF (ONLYFC) THEN
            EACTIV = D0
         ELSE
            IF (HSROHF) THEN
C              FAOSAV(:,INDEVC,1) = FDAO = FCAO + FVAO
C              FAOSAV(:,INDEVC,2) = -FVAO_exchange_only
               EACTIV = DP5 *
     &            DDOT(NNBAST,DAOSAV(1,2),1,FAOSAV(1,INDEVC,2),1)
            ELSE ! NASHT.eq.1, i.e. PV(1) = 0, i.e. no H2AC contribution.
C              FAOSAV(:,INDEVC,1) = FCAO
C              FAOSAV(:,INDEVC,2) = FVAO
C              EACTIV = DV(u,u)FC(u,u) = FC(u,u)
C                     = DVAO(p,q)FC2AO(p,q) + DVAO(p,q)H1AO(p,q)
               EACTIV = DDOT(NNBAS(IOPRHF),DAOSAV(1+IIBAS(IOPRHF),2),1,
     &                           FAOSAV(1+IIBAS(IOPRHF),INDEVC,1),1)
     &                + DDOT(NNBAS(IOPRHF),DAOSAV(1+IIBAS(IOPRHF),2),1,
     &                           H1AO(1+IIBAS(IOPRHF)),1)

C              make FAOSAV(:,INDEVC,1) = FDAO = FCAO + FVAO
               CALL DAXPY(NNBAST,D1,FAOSAV(1,INDEVC,2),1,
     &                              FAOSAV(1,INDEVC,1),1)
            END IF
         END IF

         IF (DODFT) then
            IF (HSROHF) THEN
C
C              Reset densities from Dtot, -DV to DC, DV for dft
C
               CALL DSCAL(N2BASX,-D1,DCAO(1+N2BASX),1)
               CALL DAXPY(N2BASX,-D1,DCAO(1+N2BASX),1,DCAO,1)
            END IF

#ifdef ENABLE_XCINT
            if (xcint_is_active()) then
               if (.not. xcint_uses_own_grid()) then
                  ! ulgy hack to get dalton grid, we don't care about the integrated result
                  inquire(file='num_grid', exist=file_exists)
                  if (.not. file_exists) then
                     allocate(buffer1(nbast*nbast))
                     allocate(buffer2(nbast*nbast))
                     dummy_energy = 0.0d0
                     call dft_add_ks(onlyfc, dummy_energy, dcao, ndmat,
     &                               buffer1, buffer2,
     &                               wrk(kfree), lfree, ipr_diis)
                     deallocate(buffer1)
                     deallocate(buffer2)
                  end if
               end if

               allocate(dmat_ao(nbast*nbast))
               allocate(xc_mat(nbast*nbast))
               allocate(xc_mat_tri(nnbast))

               dmat_ao     = 0.0d0
               xc_mat     = 0.0d0
               xc_mat_tri = 0.0d0
               xc_energy   = 0.0d0
               xc_num_electrons = 0.0d0
               call daxpy(nbast*nbast, 0.5d0, dcao, 1,
     &                    dmat_ao, 1)

               call xcint_integrate_rks_scf(dmat_ao,
     &                                xc_energy,
     &                                xc_mat,
     &                                xc_num_electrons)
               emy = emy + xc_energy(1)

               call dgetsp(nbast, xc_mat, xc_mat_tri)
               call daxpy(nnbast, 1.0d0, xc_mat_tri, 1,
     &                    faosav(1, indevc, 1), 1)

               deallocate(dmat_ao)
               deallocate(xc_mat)
               deallocate(xc_mat_tri)
            else
#endif /* ENABLE_XCINT */
               CALL DFT_ADD_KS(ONLYFC,EMY,DCAO,NDMAT,
     &                    FAOSAV(1,INDEVC,1),FAOSAV(1,INDEVC,2),
     &                    WRK(KFREE),LFREE,IPR_DIIS)
#ifdef ENABLE_XCINT
            end if
#endif
         END IF

CAMT If calculated add the empirical dispersion energy correction
         EMY = EMY + EDISP

C
c        we need to divide the solvent calls into mmpcm and not mmpcm
c        since we in the latter iterate between the mm and pcm for a given QM
c        density

         ESOLT = D0
         IF (.NOT. MMPCM) THEN ! do as we usually do ...
           IF (PCM .AND. (FLAG(16) .OR. QMMM .OR. QMNPMM .OR. QM3)) THEN
             CALL QUIT("PCM not compatible with other solvent models")
           END IF
           IF (QM3 .AND. .NOT. DOCCSD) THEN
             CALL MEMGET2('REAL','FCQM3',KFCSOL,NNBAST,WRK,KFREE,LFREE)
             CALL DZERO(WRK(KFCSOL),NNBAST)
             CALL QM3FCK(DAOSAV(1,1),DAOSAV(1,2),
     &                   WRK(KFCSOL),ESOLT,ENSQM3,EPOQM3,EELEL,
     &                   WRK(KFREE),LFREE,IPR_DIIS)
             CALL DAXPY(NNBAST,D1,WRK(KFCSOL),1,FAOSAV(1,INDEVC,1),1)
             CALL MEMREL('DISCTL.QM3FCK',WRK,1,KFCSOL,KFREE,LFREE)
           ELSE
             EELEL = 0.0D0
           END IF

           IF (QMMM) THEN
              IF (MMTIME) DTIME = SECOND()
              CALL MEMGET2('REAL','FQMMM',KFCSOL,NNBAST,WRK,KFREE,LFREE)
              CALL DZERO(WRK(KFCSOL),NNBAST)
              CALL QMMMFCK(DAOSAV(1,1),DAOSAV(1,2),
     &                     WRK(KFCSOL),ESOLT,
     &                     WRK(KFREE),LFREE,IPQMMM)
              CALL DAXPY(NNBAST,D1,WRK(KFCSOL),1,FAOSAV(1,INDEVC,1),1)
              CALL MEMREL('DIIS_CTL.QMMMFCK',WRK,1,
     &                    KFCSOL,KFREE,LFREE)
              IF (MMTIME) THEN
                 DTIME = SECOND() - DTIME
                 TMMFCK = TMMFCK + DTIME
              END IF
           END IF

           IF (QMNPMM) THEN

             idim = GETDIM_RELMAT(.TRUE.)
             idim_side = idim**0.5
             allocate(relay_matrix(idim_side, idim_side, 1))
             CALL COMP_RELMAT(relay_matrix,WRK(KFREE),LFREE)
             CALL WRITE_RELMAT(relay_matrix)

             IF (DOMMSUB.AND.DOMMPOL) THEN
                CALL GETDIM_MMMAT(IDIM, .TRUE.)
                allocate(mm_relay_matrix(idim))
                CALL COMP_MMRELMAT(mm_relay_matrix,WRK(KFREE),LFREE)
             END IF

             allocate(qmcmm_sol(nnbast))
             qmcmm_sol = 0.0d0
             CALL QMNPMM_FOCK(DAOSAV(1,1),DAOSAV(1,2),relay_matrix,
     &                        mm_relay_matrix,
     &                        qmcmm_sol,ESOLT,WRK(KFREE),LFREE)
             CALL DAXPY(NNBAST,1.0d0,qmcmm_sol,1,FAOSAV(1,INDEVC,1),1)

             deallocate(relay_matrix)
             if (allocated(mm_relay_matrix)) deallocate(mm_relay_matrix)
             deallocate(qmcmm_sol)

           END IF

           IF (FLAG(16)) THEN
             CALL MEMGET2('REAL','FCSOL',KFCSOL,NNBAST,WRK,KFREE,LFREE)
             CALL SOLFCK(DAOSAV(1,1),DAOSAV(1,2),
     &                   WRK(KFCSOL),DUMMY,TNLM,.FALSE.,ESOLT,
     &                   WRK(KFREE),LFREE,IPR_DIIS)
             CALL DAXPY(NNBAST,D1,WRK(KFCSOL),1,FAOSAV(1,INDEVC,1),1)
             CALL MEMREL('DISCTL.SOLFCK',WRK,1,KFCSOL,KFREE,LFREE)

           ELSE IF (PCM) THEN
             CALL MEMGET2('REAL','FCPCM',KFCSOL,NNBASX,WRK,KFREE,LFREE)
             CALL PCMFCK(DAOSAV(1,1),DAOSAV(1,2),WRK(KFCSOL),
     &                   DUMMY,.FALSE.,ESOLT,WRK(KFREE),
     &                   LFREE,IPR_DIIS)
             CALL DAXPY(NNBAST,DM1,WRK(KFCSOL),1,FAOSAV(1,INDEVC,1),1)
           END IF
         ELSE  ! the mmpcm interface is only implemented for hf/dft so
C                no reference to doccsd
           IF (.NOT. QM3 .OR. .NOT. PCM) THEN
             CALL QUIT("MMPCM have to be used with both PCM and QM3")
           END IF
C          parameters on input: MXITMP and THRSMP
           LCONV = .FALSE.
           XTEST1 = D0
           CALL MEMGET2('REAL','FCSOL',KFCSOL,NNBAST,WRK,KFREE,LFREE)
           CALL MEMGET2('REAL','MMPCM',KMMPCM,NNBAST,WRK,KFREE,LFREE)
           DO 998 ITER = 1, MXITMP
C            First the MM part ...
             CALL DZERO(WRK(KFCSOL),NNBAST)
             CALL DZERO(WRK(KMMPCM),NNBAST)
             CALL QM3FCK(DAOSAV(1,1),DAOSAV(1,2),
     &                   WRK(KFCSOL),ESOLT,ENSQM3,EPOQM3,EELEL,
     &                   WRK(KFREE),LFREE,IPR_DIIS)
             CALL DCOPY(NNBAST,WRK(KFCSOL),1,WRK(KMMPCM),1)

             TEMP = ESOLT

C            Now the PCM part ...

             CALL DZERO(WRK(KFCSOL),NNBAST)
             CALL PCMFCK(DAOSAV(1,1),DAOSAV(1,2),WRK(KFCSOL),
     &                   DUMMY,.FALSE.,ESOLT,WRK(KFREE),
     &                   LFREE,IPR_DIIS)
             CALL DAXPY(NNBAST,DM1,WRK(KFCSOL),1,WRK(KMMPCM),1)

             ESOLT = ESOLT + TEMP
             XTEST2 = ESOLT
             IF ( ((ABS(XTEST1-XTEST2)) .LT. THRSMP) .OR. LOSPC) THEN
               LCONV = .TRUE.
               CALL DAXPY(NNBAST,D1,WRK(KMMPCM),1,FAOSAV(1,INDEVC,1),1)
               write(lupri,*) 'MMPCM converged for given density'
               write(lupri,*) 'No. of iterations: ',ITER
               GOTO 999
             ELSE
               XTEST1 = XTEST2
             END IF
 998       CONTINUE
 999       CONTINUE
           IF (.NOT. LCONV) write(lupri,*) 'MMPCM not converged !'
           CALL MEMREL('DISCTL.PCMQM3FCK',WRK,1,KFCSOL,KFREE,LFREE)
         END IF
C
         IF (USE_PELIB()) THEN
           IF (FLAG(16) .OR. PCM .OR. QMMM .OR. QMNPMM .OR. QM3) THEN
            CALL QUIT('.PEQM incompatible with other embedding models')
           END IF
           ALLOCATE(PE_DENMAT(NNBASX), PE_FCKMAT(NNBASX))
           IF (NASHT == 0 .OR. HSROHF) THEN
             PE_DENMAT = DAOSAV(:,1)
           ELSE
             PE_DENMAT = DAOSAV(:,1) + DAOSAV(:,2)
           END IF
           CALL PELIB_IFC_FOCK(PE_DENMAT, PE_FCKMAT, PE_ENERGY)
           ESOLT = PE_ENERGY
           FAOSAV(:,INDEVC,1) = FAOSAV(:,INDEVC,1) + PE_FCKMAT
           DEALLOCATE(PE_DENMAT, PE_FCKMAT)
         END IF
C
#if defined (HAS_PCMSOLVER)
         if (pcmmod_is_pcm_calculation) then
! Allocate scratch space for the PCM Fock matrix contribution in AO basis.
           allocate(fock_pcm_ao(nnbasx))
           fock_pcm_ao = 0.0d0
! Allocate scratch space for the density matrix to be passed to the PCM
! subroutines
           allocate(density_ao_symm(nnbasx))
           density_ao_symm = 0.0d0
! Allocate the PCM Fock matrix contribution in AO basis. 
! It is of dimension NNBASX and we take the memory from the WRK array. 
! We zero it out afterwards.
!           call memget('REAL', kfcsol, nnbasx, wrk, kfree, lfree)
!           call dzero(wrk(kfcsol), nnbasx)
! DAOSAV(:, 1) = DCAO (density matrix for closed shells)
! DAOSAV(:, 2) = DVAO (density matrix for open shells)
! PCM doesn't care about closed/open shells, so we decide at this
! higher level which one of them to pass, based on the wavefunction asked.          
           if (nasht == 0 .or. (nasht > 0 .and. dodft) .or. hsrohf) then
             density_ao_symm(:) = daosav(:, 1)
! We pass just DCAO, i.e. density_matrix = DCAO
           else
! We sum DCAO and DVAO, i.e. density_matrix = DCAO + DVAO
             density_ao_symm(:) = daosav(:, 1) + daosav(:, 2)
           end if
! Call the energy driver with the symmetry-unpacked density matrix 
           if (nsym > 1) then
             allocate(density_ao(nnbasx))
             density_ao = 0.0d0
! Symmetry un-pack density_ao_symm to density_ao             
             call pksym1(density_ao, density_ao_symm, nbas, nsym, -1)
             call pcm_energy_driver(density_ao, esolt, wrk(kfree), 
     &                               lfree)
           else
             call pcm_energy_driver(density_ao_symm, esolt, wrk(kfree), 
     &                              lfree)
           end if
! Call the Fock driver           
           call pcm_oper_ao_driver(fock_pcm_ao, 'TotASC'//char(0), 
     &                             wrk(kfree), lfree)
! Copy the PCM Fock matrix contribution on top of the AO basis Fock matrix
           if (nsym > 1) then
             allocate(fock_pcm_ao_symm(nnbasx))
             fock_pcm_ao_symm = 0.0d0
! Symmetry-pack fock_pcm_ao to fock_pcm_ao_symm 
             call pksym1(fock_pcm_ao, fock_pcm_ao_symm, nbas, nsym, +1)
             call daxpy(nnbast, -1.0d0, fock_pcm_ao_symm, 1, 
     &                                  faosav(1, indevc, 1), 1)
           else
             call daxpy(nnbast, -1.0d0, fock_pcm_ao, 1, 
     &                                  faosav(1, indevc, 1), 1)
           end if
! We can print it here, if needed.
!lf            write(lupri,*) 'extpcm fock matrix'
!lf            call outpak(wrk(kfcsol), nbast, 1, lupri)
! Clean-up
           deallocate(fock_pcm_ao)
           deallocate(density_ao_symm)
           if (nsym > 1) then
             deallocate(fock_pcm_ao_symm)
             deallocate(density_ao)
           end if
         end if
#endif

C
         EMCDIF = EMCSCF
         EMCSCF = EMY + EACTIV + EPOT + ESOLT

         IF (QM3 .AND. (.NOT. OLDTG)) EMCSCF = EMCSCF - EELEL
C        this also includes the case of MMPCM

         EGSAV(1,INDEVC) = EMCSCF

         EMCLOW = MIN(EMCLOW,EMCSCF)
         EMCDIF = EMCSCF - EMCDIF
         EMCDIFERR = THRINC*ABS(EMCLOW)
         GRDNRM_old = GRDNRM
C
         CALL DIIS_EVC(SMOAO,DAOSAV(1,1),DAOSAV(1,2),H1AO,
     &               FAOSAV(1,INDEVC,1),FAOSAV(1,INDEVC,2),CMO1,
     &               EVCSAV(1,INDEVC),WRK,KFREE,LFREE,IPR_DIIS)
C        CALL DIIS_EVC(SMOAO,DCAO,DVAO,H1AO,FDAO,FVAO,CMO1,
C    &               ERRVEC,WRK,KFRSAV,LFRSAV,IPR_DIIS)
         IF (IPR_DIIS .GE. 20) THEN
            WRITE (LUPRI,*) ' Test output of all of EVCSAV, newest is',
     &         INDEVC
            CALL OUTPUT(EVCSAV,1,NNORBT,1,NEVC,NNORBT,NEVC,-1,LUPRI)
         ELSE IF (IPR_DIIS .GE. 10) THEN
            WRITE (LUPRI,*) ' Error vector, DIIS iteration',ITDIIS
            CALL OUTPKB(EVCSAV(1,INDEVC),NORB,NSYM,-1,LUPRI)
         END IF
C
C Project out .FREEZE'ed orbitals from the gradient matrix
C
         IF (LNOROT) CALL PHPT(WRK(KPFROZ),WRK(KTMP1),NBAST)
         GRDNRM = DNRM2(NNORBT,EVCSAV(1,INDEVC),1)
         GRD_CHANGE = GRDNRM_old / GRDNRM
         EGSAV(2,INDEVC) = GRDNRM
         IF (LUW4 .NE. LUPRI .OR. IPR_DIIS .GT. 1) THEN
            WRITE (LUPRI,'(/A/A,I4,1P,G22.12,D15.5,D12.2,I5,L5)')
     &' ITDIIS, energy, error norm, change in energy, scr.thr., difden:'
     &,'@ DIIS iter.',ITDIIS,EMCSCF,GRDNRM,EMCDIF,-IFTHRS,DIFDEN
            IF (AUTOCC) THEN
               WRITE (LUPRI,'(/A,8I4)')
     &         '@ AUTO OCC: current SCF occupation',(NISH(I),I=1,NSYM)
               IF (IOPRHF .NE. 0) WRITE (LUPRI,'(A,I4)')
     &         '@ AUTO OCC: current open shell symmetry',IOPRHF
            END IF
            CALL FLSHFO(LUPRI)
         END IF
         IF (LUW4 .NE. LUPRI .OR. IPR_DIIS .LE. 1) THEN
#if defined (HAS_PCMSOLVER)
            IF (EMBEDDING .or. pcmmod_is_pcm_calculation) then
#else                
            IF (EMBEDDING) THEN
#endif                    
C           this also includes the case of MMPCM
               IF (.NOT.AUTOCC) THEN
                  WRITE (LUW4,'(A,I3,1P,2G20.12,D15.5,D12.2)')
     &                 '@',ITDIIS,EMCSCF,ESOLT,GRDNRM,EMCDIF
               ELSE IF (IOPRHF .EQ. 0) THEN
                  WRITE (LUW4,'(A,I3,1P,2G20.12,2D11.2,2X,8I4)')
     &                 '@',ITDIIS,EMCSCF,ESOLT,GRDNRM,EMCDIF,
     &                 (NISH(I),I=1,NSYM)
               ELSE
                  WRITE (LUW4,'(A,I3,1P,2G20.12,2D11.2,I4,2X,8I4)')
     &                 '@',ITDIIS,EMCSCF,ESOLT,GRDNRM,EMCDIF,
     &                 IOPRHF,(NISH(I),I=1,NSYM)
               END IF
            ELSE
               IF (.NOT.AUTOCC) THEN
                  WRITE (LUW4,'(A,I3,1P,G22.12,D15.5,D12.2,I5)')
     &                 '@',ITDIIS,EMCSCF,GRDNRM,EMCDIF,NEVC
               ELSE IF (IOPRHF .EQ. 0) THEN
                  WRITE (LUW4,'(A,I3,1P,G20.12,2D11.2,2X,8I4)')
     &                 '@',ITDIIS,EMCSCF,GRDNRM,EMCDIF,
     &                 (NISH(I),I=1,NSYM)
               ELSE
                  WRITE (LUW4,'(A,I3,1P,G20.12,2D11.2,I4,2X,8I4)')
     &                 '@',ITDIIS,EMCSCF,GRDNRM,EMCDIF,
     &                 IOPRHF,(NISH(I),I=1,NSYM)
               END IF
            END IF
            CALL FLSHFO(LUW4)
         END IF
         IF (MOLDEN) CALL MOLDEN_SCFCON(ITDIIS,EMCSCF,.FALSE.)

         IF (GRDNRM .LE. THDIIS) THEN  ! converged

            ICONV = 1

C           No level shift at convergence, otherwise orbital energies
C           for e.g. CC will not be correct
            SHFTLVL = 0.0D0

            IF (MOLDEN) CALL MOLDEN_SCFCON(ITDIIS,EMCSCF,.TRUE.)

            WRITE (LUPRI,4110) ITDIIS,EMCSCF,GRDNRM
            IF (LUW4 .NE. LUPRI) WRITE (LUW4,4110) ITDIIS,EMCSCF,GRDNRM
            IF (IPR_DIIS .GE. -1) THEN
              IF (CHOINT .AND. (IPR_DIIS .GE. 3)) THEN
              !  Cholesky
                 WRITE(LUPRI,4116) 'Coulomb part ',CSIRF
                 WRITE(LUPRI,4116) 'Exchange part',XSIRF
                 WRITE(LUPRI,4116) 'Cholesky I/O ',RSIRF
                 WRITE(LUPRI,4116) 'MO transform.',OSIRF
                 WRITE(LUPRI,4116) 'MO sorting   ',SSIRF
              END IF
              WRITE(LUPRI,4115) TSIRF
              IF (QMMM) CALL QMMMTIMES('SIRIUS')
              IF (LUW4 .NE. LUPRI) THEN
                 IF (CHOINT .AND. (IPR_DIIS .GE. 3)) THEN
                    WRITE(LUPRI,4116) 'Coulomb part ',CSIRF
                    WRITE(LUPRI,4116) 'Exchange part',XSIRF
                    WRITE(LUPRI,4116) 'Cholesky I/O ',RSIRF
                    WRITE(LUPRI,4116) 'MO transform.',OSIRF
                    WRITE(LUPRI,4116) 'MO sorting   ',SSIRF
                 END IF
                 WRITE(LUW4,4115) TSIRF
              END IF
            END IF
C           Save converged CMO
            WRITE (CMOLBL,'(A4,I4)') 'DIIS',ITDIIS
            CALL NEWORB(CMOLBL,CMO,LNEWOCC)
C           ... REWIT1 false: do not destroy any GEOWALK information
            GO TO 1000
         ELSE IF (ITDIIS .GE. MXDIIS) THEN
            ICONV = 0
            WRITE (LUPRI,4120)
            IF (LUW4 .NE. LUPRI) WRITE (LUW4,4120)
C
C           Too many DIIS iterations, save final CMO before exit
            WRITE (CMOLBL,'(A4,I4)') 'DIIS',ITDIIS
            CALL NEWORB(CMOLBL,CMO,LNEWOCC)
C           ... REWIT1 false: do not destroy any GEOWALK information
            GO TO 1100
         END IF
 4110 FORMAT(/'@ *** DIIS converged in',I4,' iterations !',
     &       /'@     Converged SCF energy, gradient:',F20.12,1P,D12.2)
 4115 FORMAT( '    - total time used in SIRFCK :' ,F18.2,' seconds')
 4116 FORMAT( '    - total time used for ',A13,':',F12.2,' seconds')
 4120 FORMAT(/
     & 'WARNING !!! DIIS aborted because max DIIS iterations reached !')
C
Chj 990819: use DIFDEN more often
Chj was: IF (IFTHRS .GE. 12 .AND. .NOT. FLAG(49)) DIFDEN = .TRUE.
         NEWTHR = .FALSE.
         IF (DIRFCK .AND. .NOT. USRSCR) THEN
            ITHRS = INT(-LOG10(GRDNRM/SCRFAC))+3
C           ... "+3" instead of "+1" because we must
C               have sufficient accuracy for lowering the
C               gradient norm with a factor 100 and for
C               a significant contribution in CVEC later ...
C               (comment inserted Feb 2001/Jul 2006/hjaaj)
            IF (DODIFDEN) THEN
               ILIM = 2
Chj-aug99: only NEWTHR with skip of 2 for DIFDEN
               IF (ITHRS .GE. IFTMAX-2) ITHRS = IFTMAX
Chj-aug99: for effective DIFDEN; note this test is OK w. ILIM.eq.2
            ELSE
               ILIM = 1
            END IF
            IF (.NOT. AOSAVE .AND. ABS(ITHRS-IFTHRS) .GE. ILIM) THEN
               NEWTHR = .TRUE.
               IFTHRS = MIN(IFTMAX,MAX(IFTMIN,ITHRS))
            END IF
         END IF
Chj 990819/000103 new:
         DIFDEN = .NOT. NEWTHR .AND. DODIFDEN
C        DIFDEN true if not new threshold and not "no difden"
C        -- if we have tightened screening threshold we must
C           calculate w/o DIFDEN or the error in the Fock matrix
C           would correspond to the previous screening threshold.
C
#ifdef USE_OPT_SOLVERS
         CALL  DIISWEIGHTS(METHOD,ENERG,H1AO,SAOUNP,
     &                     DAOSAV,FAOSAV,EVCSAV,
     &                     INDEVC,NEVC,GVEC,BMAT,CVEC,MXERRV,
     &                     WRK,KFREE,LFREE, WRK(KPFROZ))

#else    /* i.e. do not USE_OPT-SOLVERS */
C
C        Test if convergence is stopping because of numerical problems
C        (EMCDIF test should make sure we are in local region)
C
ckr         IF (ITNOCC .GT. 3 .AND. ABS(EMCDIF) .LT. THDDEF .AND.
ckr     &       GRDNRM .GT. CNVFAC*GRDNRM_OLD .AND. .NOT. NEWTHR) THEN
ckr            ICONV = 0
ckr            WRITE (LUPRI,4140)
ckr            IF (LUW4 .NE. LUPRI) WRITE (LUW4,4140)
ckr            GO TO 1100
ckr         END IF
ckr 4140 FORMAT(' DIIS aborted, convergence too slow !')
C
C hjaaj apr 2002
C   - check if a previous iteration had a lower gradient and a higher energy,
C     if yes: eliminate this entry and earlier entries, because they
C     are not relevant any more and may cause slower convergence
C     or - even worse - that the convergence is stalled.
C     Reason: DIIS tries to minimize gradient norm, not energy.
C   - hjaaj mar 2004: only remove entries with lower gradient and
C     higher energy, not "earlier entries"
C
         J_removed = 0
         IF (INDEVC.EQ.NEVC .AND. ITDIIS .GT. 3 .AND.
     &       GRDNRM .LE. 10.0D0) THEN
C        ... no wrap around yet - INDEVC .lt. NEVC requires more programming
C        ITDIIS/GRDNRM check: give initial oscillations a chance to settle down
C          before checking ... experience shows this will often be advantageous
C          (see below) /hjaaj
            IEND = NEVC - 1
            DO I = 1,IEND
C           IF (I .NE. INDEVC) THEN
               IF (EGSAV(1,I) .GT. EMCSCF + EMCDIFERR .AND.
     &             EGSAV(2,I) .LT. GRDNRM) THEN
C                  ... Oops! a previous iter with higher E and lower grd
                  J_removed = J_removed + 1
                  EGSAV(1,I) = D0
               END IF
C           END IF
            END DO
C
C           remove entries with higher E and lower grd, if any
C
            IF (J_removed .GT. 0) THEN
               WRITE(LUPRI,'(/A/A,I3,A)')
     &   'Info: SCF gradient has been lower than now,',
     &   '      therefore',J_removed,
     &   ' old iterations are removed from DIIS.'
               K = 0
               DO J = 1, NEVC
                  IF (EGSAV(1,J) .NE. D0) THEN
                     K = K + 1
                  IF (K .NE. J) THEN
                     EGSAV(1,K) = EGSAV(1,J)
                     EGSAV(2,K) = EGSAV(2,J)
                     CALL DCOPY(NNORBT,EVCSAV(1,J),1,
     &                                 EVCSAV(1,K),1)
                     CALL DCOPY(NNBAST,FAOSAV(1,J,1),1,
     &                                 FAOSAV(1,K,1),1)
                     IF (.NOT. ONLYFC)
     &               CALL DCOPY(NNBAST,FAOSAV(1,J,2),1,
     &                                 FAOSAV(1,K,2),1)
                     IF (J .LT. NEVC) THEN
C                       ... the new entry (i.e. J.EQ.NEVC)
C                           is added in DIIS_RED below
                        KROW = K*(K-1) / 2
                        DO L = 1,K
                           BMAT(KROW+L) =
     &                        DDOT(NNORBT,EVCSAV(1,L),1,EVCSAV(1,K),1)
                        END DO
                     END IF
                  END IF
                  END IF
               END DO
               JTDIIS  = K
               NEVC    = K
               INDEVC  = K
            END IF
         END IF
C
         CALL DIIS_RED(C2DIIS,BMAT,INDEVC,NEVC,DAMP,THREVC,
     &                 CVEC,XLMBDA,EVCSAV,NNORBT,
     &                 WRK,KFREE,LFREE,IPR_DIIS)
C        CALL DIIS_RED(C2DIIS,BMAT,INDEVC,NEVC,DAMP,THREVC,
C    &                 CVEC,XLMBDA,ERRVEC,L_ERRVEC,
C    &                 WRK,KFRSAV,LFRSAV,IPR_DIIS)
C
C
C        Check to see if DIIS iterations are stalled,
C          if "yes" then restart DIIS /Feb-2001 hjaaj
C        JTDIIS check: we need at least two Fock matrices to do this ...
C        ITDIIS/GRDNRM check: give initial oscillations a chance to settle down
C          before checking ... experience shows this will often be advantageous
C          (no. 3 will often be lower in energy than no. 1, even if no. 2
C           is higher in energy with lower gradient, especially if MOSTART H1DIAG)
C        Note that we only "stall" if energy goes up while gradient goes
C          down or vice versa because DIIS minimizes gradient and it is
C          OK wrt. DIIS if both goes up or both goes down /Feb-2004 hjaaj
!  Mar. 2011 hjaaj: if abs(EMCDIF) .lt. 1.0D-5 the error can be grid related in DFT
         IF (DODFT .AND. ABS(EMCDIF) .LT. 1.0D-5) THEN
             EMCDIF_test = 0.0D0
         ELSE
             EMCDIF_test = EMCDIF
         END IF
C        Do not try that with open-shell - it does not work./pawsa
!        IF (ONLYFC .AND. JTDIIS.GT.1 .AND. ITDIIS.GT.3 .AND.
         IF (JTDIIS.GT.1 .AND. ITDIIS.GT.3 .AND. J_removed.EQ.0 .AND.
     &       GRDNRM .LE. 10.0D0 .AND.
     &       ( ABS(CVEC(INDEVC)) .LT. 0.1D0 .OR.
     &         ABS(CVEC(INDEVC)) .GT. 5.0D0 .OR.
     &         EMCDIF_test*(GRDNRM-GRDNRM_OLD) .LT. -1.0D-10 )
     &        .AND. BCKSTP)  THEN
C              .. test EMCDIF*grddif against -1.0D-10 instead of 0.0D0
C                 to avoid round-off problems here /hjaaj Sep 2005

             IF (I_STALL .LT. 1) THEN ! do not stop first time energy goes up
                I_STALL = I_STALL + 1
                GO TO 300
             END IF
             I_STALL = 0

             WRITE(LUPRI,'(1P/A,2(/A,G11.3),/A,I3,A)')
     &       '!!! Info: DIIS restarted because it was stalled ...',
     &       ' - energy changed by  ',EMCDIF,
     &       ' - gradient changed by',GRDNRM - GRDNRM_OLD,
     &       ' - or strange C vector coeff. for current index (=',
     &       INDEVC,') :'
             WRITE (LUPRI,'(10F12.6)') (CVEC(I),I=1,NEVC)
             IF (EMCDIF .LE. ABS(EMCSCF)*1.0D-12) THEN
                JNDEVC = INDEVC
             ELSE
C               ... use previous Fock matrix if energy increased
                WRITE(LUPRI,'(/A)') '! Backstep to previous Fock'//
     &             ' matrix because energy increased.'
                JNDEVC = MOD(JTDIIS-2,MXERRV) + 1
C
C               restore old CMO from LUIT1 (for backstep);
C               no DIFDEN because DAOSAV is for current, not previous iteration
C
                CALL READMO(CMO,9,DUMMY,1)
                DIFDEN = .FALSE.
             END IF
             IF (JNDEVC .NE. 1) THEN
                CALL DCOPY(NNORBT,EVCSAV(1,JNDEVC),1,
     &                            EVCSAV(1,1     ),1)
                CALL DCOPY(NNBAST,FAOSAV(1,JNDEVC,1),1,
     &                            FAOSAV(1,1     ,1),1)
                IF (.NOT. ONLYFC)
     &          CALL DCOPY(NNBAST,FAOSAV(1,JNDEVC,2),1,
     &                            FAOSAV(1,1     ,2),1)
             END IF
             EGSAV(1,1) = EGSAV(1,JNDEVC)
             EGSAV(2,1) = EGSAV(2,JNDEVC)
             CVEC(1) = D1
             JTDIIS  = 1
             NEVC    = 1
             INDEVC  = 1
             BMAT(1) = BMAT( (JNDEVC*(JNDEVC+1))/2 )
C            also restrict orbital rotation with a level shift
c            the formula is empirical EMCDIF is positive here.
             SHFTLVL = MAX(-12D0,SHFTLVL - MIN(4.0D0,ABS(EMCDIF)*10.D0))
c            and slow down  the level shifting annealing...
             SHFTFAC = MIN(SHFTFAC + (1.D0-SHFTFAC)*0.25,0.85D0)
         ELSE IF (ITDIIS .GT. 1) THEN
             IF (GRD_CHANGE .GE. 2.0D0) SHFTLVL =  SHFTLVL*SHFTFAC
c            reduce level shift if convergence is good
             IF (GRD_CHANGE .GE. 4.0D0) SHFTLVL =  SHFTLVL*SHFTFAC
c            double reduce level shift if convergence is excellent
             IF (SHFTLVL .GE. -1.D-2) SHFTLVL = 0.0D0
         END IF
  300    CONTINUE
#endif    /* USE_OPT-SOLVER */

C        Level shift is only working for closed shell calculations ...
C          (open shell problems reported by users ...) / hjaaj 22-Jun-2005
!        IF ( NASHT .GT. 0 ) SHFTLVL = 0.0D0
!        enabled again 6-Nov-2014 after testing a problematic open-shell DFT
!        calculation which converged with level shift, but not without. /hjaaj

C        Iteration accepted; save current CMO for restart
         WRITE (CMOLBL,'(A4,I4)') 'DIIS',ITDIIS
         CALL NEWORB(CMOLBL,CMO,LNEWOCC)
C        ... REWIT1 false: do not destroy any GEOWALK information
         LNEWOCC = .FALSE.

C        print atomic populations in each iteration ? /hjaaj
         if (lim_poppri .gt. 0)
     &      call sirpop('DIIS ',DV,wrk(kfree),lfree)

C
         IF (HSROHF) THEN
C
C Update FCAO with the FV correction in AO basis
C
            CALL DAXPY(NNBAST,D1,H1AO,1,FAOSAV(1,INDEVC,1),1)
            CALL DCOPY(NNBAST,FAOSAV(1,INDEVC,1),1,FCAO,1)
            CALL FVCORR(1,D1,DAOSAV,FAOSAV(1,INDEVC,1),MXERRV,FCAO,
     &            WRK,KFREE,LFREE )
C
C Update the latest FAOSAV with the effective fock matrix
C
            CALL DCOPY(NNBAST,FCAO,1,FAOSAV(1,INDEVC,1),1)
            CALL DZERO(FCAO,NNBAST)
         ELSE
            CALL DCOPY(NNBAST,H1AO,1,FCAO,1)
         END IF
         CALL DGEMM('N','N',NNBAST,1,NEVC,1.D0,
     &              FAOSAV(1,1,1),NNBAST,
     &              CVEC,NEVC,1.D0,
     &              FCAO,NNBAST)
         IF (.NOT. ONLYFC .AND..NOT.HSROHF) THEN
         CALL DGEMM('N','N',NNBAST,1,NEVC,1.D0,
     &              FAOSAV(1,1,2),NNBAST,
     &              CVEC,NEVC,0.D0,
     &              FCAO(1+NNBAST),NNBAST)
         END IF
C
C        Transform FDAO to FD (FD saved in DCAO).
C        Find FV and make open-shell Fock matrix
C
         CALL UTHUB(FCAO,DCAO,CMO,WRK(KFREE),NSYM,NBAS,NORB)
         IF (HSROHF .OR. ONLYFC) THEN
            IF (IPR_DIIS .GE. 15) THEN
               IF (HSROHF) THEN
                  WRITE (LUPRI,'(/A)')
     &            ' Fock matrix after FV correction'
               ELSE
                  WRITE (LUPRI,'(/A)')
     &            ' Fock matrix before diagonalization'
               END IF
               CALL OUTPKB(DCAO,NORB,NSYM,-1,LUPRI)
            END IF
         ELSE
           DO 330 ISYM = 1,NSYM
            IF (NASH(ISYM) .EQ. 0) GO TO 330
            IF (IPR_DIIS .GT. 15) THEN
               WRITE (LUPRI,'(/A,I3)')
     &            ' Fock matrix before FV correction, symmetry',ISYM
               CALL OUTPAK(DCAO(1+IIORB(ISYM)),NORB(ISYM),-1,LUPRI)
            END IF
            CALL UTHU(FCAO(1+NNBAST+IIBAS(ISYM)),DCAO(1+NNBAST),
     &         CMO(1+ICMO(ISYM)),WRK(KFREE),NBAS(ISYM),NORB(ISYM))
            DO JACT = 1, NASH(ISYM)
            KACT = NISH(ISYM) + JACT
            KROW = KACT*(KACT-1)/2
            DO 310 J = 1, NISH(ISYM)
               JFKJ = IIORB(ISYM) + KROW + J
               DCAO(JFKJ) = DCAO(JFKJ) + DCAO(NNBAST+KROW+J)
  310       CONTINUE
            DO 320 J = NISH(ISYM)+NASH(ISYM) + 1,NORB(ISYM)
               JROW = J*(J-1)/2
               JFJK = IIORB(ISYM) + JROW + KACT
               DCAO(JFJK) = DCAO(JFJK) - DCAO(NNBAST+JROW+KACT)
  320       CONTINUE
            END DO
            IF (IPR_DIIS .GE. 15) THEN
               WRITE (LUPRI,'(/A,I3)')
     &            ' Fock matrix after FV correction, symmetry',ISYM
               CALL OUTPAK(DCAO(1+IIORB(ISYM)),NORB(ISYM),-1,LUPRI)
            END IF
  330      CONTINUE
         END IF
C
C        If requested (by SHFTLVL .ne. 0) :
C        Perform "restricted step" step reduction (i.e. level shift)
C        by adding 0.5*SHFTLVL*DENS(i,i) to all diagonal elements;
C        /hjaaj-Feb.2001

         IF (SHFTLVL .NE. D0) THEN
            WRITE(LUPRI,'(I4,A,1P,D10.2)') ITDIIS,'  Level shift: '//
     &      'doubly occupied orbital energies shifted by',SHFTLVL
            IF (NASHT .GT. 0) WRITE(LUPRI,'(T16,A,1P,D10.2)')'and '//
     &      'singly occupied orbital energies shifted by',0.5D0*SHFTLVL
         END IF
chjaaj: improve comment with something about "restricted step"
chjaaj: and "0.5*occ(i)* " /23may09 hjaaj
C
C        Diagonalize (level-shifted) FC (saved in DCAO)
C
         CALL ICOPY(8,NISH,1,NISH_save,1)
         JOPRHF = IOPRHF
         IF (GRDNRM_old > 0.0D0) THEN
            GRDNRM_factor = MAX(12.0D0, GRDNRM_old/GRDNRM)
         ELSE
            GRDNRM_factor = 12.0D0
         END IF
         IF ( GRDNRM .LE. MAX(1.D-5, THDIIS*GRDNRM_factor) ) THEN
            THR_FCKEIG = 0.8D0*THDIIS
            CALL ICOPY(8,NORB,1,NCAN,1) ! we need all orbitals canonical for CC
         ELSE
            THR_FCKEIG = 0.01D0*GRDNRM
            CALL ICOPY(8,NOCC,1,NCAN,1) ! we only need oocupied orbitals canonical in global region
         END IF

         call gettim(t1,w1)
         CALL FCKEIG(CMO,DCAO,SHFTLVL,THR_FCKEIG,NCAN,IPR_DIIS,
     &               WRK(KFREE),LFREE)
         IF (IPR_DIIS .GT. 2) THEN
            call gettim(t2,w2)
            write(lupri,'(A,1P,3D10.2)')
     &         ' CPU and wall times used in FCKEIG, THR_FCKEIG: ',
     &         t2-t1,w2-w1,thr_fckeig
         END IF

C        Test if AUTOCC has changed occupation :

         IF (AUTOCC) THEN
          NTEST = 0
          DO ISYM = 1,8
            NTEST = NTEST + ABS(NISH(ISYM)-NISH_save(ISYM))
          END DO
          NTEST = NTEST + ABS(IOPRHF-JOPRHF)
          IF (NTEST .GT. 0) THEN
            NEWOCC = NEWOCC + 1
            LNEWOCC = .TRUE.
            ITNOCC = 0
            IF (NEWOCC .GT. 6) THEN
               WRITE (LUPRI,4210) NEWOCC
               IF (LUW4 .NE. LUPRI) WRITE (LUW4,4210) NEWOCC
               CALL QUIT(
     &           'FATAL ERROR: problems with automatic SCF occupation')
            END IF
Chj-sep99: bugfix for open shell when open shell symmetry changed
C          (NASH(ISYM) used for test in DIIS_CTL !!!)
            IF (NASHT .EQ. 1) THEN
               CALL IZERO(NASH,8)
               NASH(IOPRHF) = 1
            END IF
Chj-sep99: and change index arrays to new occupation (for RHFENR and ?)
            CALL SETORB
          END IF
         END IF
 4210 FORMAT(/' DIIS has changed occupation numbers',I2,' times now.'
     &       /' Program aborts because this indicates problems with'//
     &        ' automatic occupation.')

C     --> Next DIIS iteration
      GO TO 100
C
C 1000: converged
C 1100: not converged
 1000 CONTINUE
C
C     QCHF currently needed for solvent and for writing SIRIFC:
C     reset to not converged in these cases. 951130-hjaaj.
C     We write to SIRIFC here instead, thus no need to reset ICONV,K.Ruud-May97
C     hjaaj aug99: revised code for SIRIFC to fix some potential problems,
C                  and call tractl if requested with ITRFIN
C
C     If converged (ICONV .eq. 1) and if this is final level of wave function:
C
C     Print orbital energy analysis for RHF if QCHF not called
C     because ICONV = 1 and if either IPR_DIIS .gt. 0 or
C     this is final level of wave function.
C
C
C
C           Check if we need new RHFWOP because of AUTOCC
C           (inserted by Kenneth Ruud May 95, open shell hjaaj Aug 95)
C           Moved from SIRCTL to here since we do not necessarily use
C           quadratic HF after the DIIS anymore....
C
      RHFWOP = .TRUE.
      IF (IOPRHF .NE. IOPRHF_old) RHFWOP = .FALSE.
      DO ISYM = 1, NSYM
         IF (NRHF(ISYM) .NE. NISH(ISYM)) THEN
            NRHF(ISYM) = NISH(ISYM)
            RHFWOP = .FALSE.
         END IF
      END DO
      IF (.NOT.RHFWOP) THEN
         OLDWOP = .FALSE.
C
         IF (.NOT.DOMC) WRINDX = .TRUE.
C        ... We need to correct orbital rotation information on file
C            for response and ABACUS
C            if occupation has changed and RHF determines orbitals
C
         JWOPSY = 1
         CALL SIRSET(WRK(KFREE),LFREE,OLDWOP)
         IAVERR = 0
         CALL AVECHK(WRK(KFREE),LFREE,IAVERR)
         IF (IAVERR .NE. 0) CALL QUIT(
     &      'SIRCTL DIIS error: inconsistency in sup.sym. after AUTOCC')
         RHFWOP = .TRUE.
      END IF
C
      IF ( ICONV .NE. 0. .AND.
     &     ( IPR_DIIS .GT. 0 .OR.
     &       CMOPRI .OR.
     &       .NOT.(DOMP2 .OR. DOCI .OR. DOMC)) ) THEN
         IPRENR = 1
      ELSE
         IPRENR = 0
      END IF

!     IF ( ( FLAG(25) .OR. INERSI ) .AND.
!    &     (ICONV .NE. 0. .AND. .NOT. DOMC) ) THEN
      IF (ICONV .NE. 0) THEN
      ! hjaaj nov 2010: now always write SIRIFC if converged - so we
      ! have it available afterwards for post-SCF programs. This is more
      ! modular than making dependent on explicit checking for a
      ! post-SCF program, and we can use it in a restart of Dalton.
         IWRIFC = 1
      ELSE
         IWRIFC = 0
      END IF
C
      IF (IPRENR .GT. 0 .OR. IWRIFC .GT. 0) THEN
C
C        Construct FDAO (FDAO saved in DCAO)
C        Transform FDAO to FD (FD saved in FCAO).
C
         CALL DCOPY(NNBAST,FAOSAV(1,INDEVC,1),1,DCAO,1)
         CALL DAXPY(NNBAST,D1,H1AO,1,DCAO,1)
         CALL UTHUB(DCAO,FCAO,CMO,WRK(KFREE),
     &              NSYM,NBAS,NORB)
      END IF
C
C     if qm3 .or. pcm then also mmpcm automatically
      IF (IPRENR .GT. 0) THEN
C
         WRITE (LUW4,'(//A)') ' *** SCF orbital energy analysis ***'
         IF (FLAG(16) .OR. PCM .OR. QM3 .OR. QMMM .OR. QMNPMM)
     &     WRITE (LUW4,'(A)') '    (incl. solvent contribution)'
         IF (USE_PELIB()) WRITE(LUW4,'(A)')
     &    '   (incl. contribution from polarizable embedding potential)'
         CALL RHFENR(IPRI4,LUW4,FCAO, WRK(KFREE),LFREE)
C-------------------
         IF (LUPRI .NE. LUW4) THEN
            WRITE (LUPRI,'(//A)') ' *** SCF orbital energy analysis ***'
            IF (FLAG(16) .OR. PCM .OR. QM3 .OR. QMMM .OR. QMNPMM)
     &        WRITE (LUPRI,'(A)') '    (incl. solvent contribution)'
            IF (USE_PELIB()) WRITE(LUPRI,'(A)')
     &    '   (incl. contribution from polarizable embedding potential)'
            CALL RHFENR(IPRI6,LUPRI,FCAO, WRK(KFREE),LFREE)
         END IF
C-------------------
C        CALL RHFENR(IPRINT,LUPRI,FC,FV,SCRA,LSCRA)
C
         IF (IPR_DIIS .GT. 10) THEN
            WRITE(LUPRI,'(/A)')
     &      'Final SCF orbitals after DIIS iterations:'
            CALL PRORB(CMO,(IPR_DIIS .LT. 20),LUPRI)
         END IF
C
      END IF
      IF (USE_PELIB()) THEN
        ALLOCATE(PE_DENMAT(NNBASX))
        IF (NASHT == 0 .OR. HSROHF) THEN
            PE_DENMAT = DAOSAV(:,1)
        ELSE
            PE_DENMAT = DAOSAV(:,1) + DAOSAV(:,2)
        END IF
        CALL PELIB_IFC_ENERGY(PE_DENMAT)
        DEALLOCATE(PE_DENMAT)
      END IF
C
C     If (FLAG(25) .OR. INERSI) write to SIRIFC
C       FLAG(25) is true for .ABACUS (geometry optimization and many properties)
C                         or .RESPONS (response calculation)
C                         or .INTERFACE (request for write of SIRIFC)
C                         or .WESTA
C                         or .CC
C       INERSI is true if this is initial state calculation
C              in a solvent calculation with inertial polarization
C
      IF ( IWRIFC .EQ. 1 ) THEN
         IF (IPRI6 .GE. 0) WRITE (LUPRI,'(/A)')
     &      ' >>> Writing SIRIFC interface file'
         IF (ICI0 .EQ. 6) THEN
Chj      ... if (GEOWLK) then
            IF (.NOT. OPTNEW) WRITE (LUPRI,'(//A/)')
     &            ' Check ratio for geometry walk with converged SCF :'
            JWSTEP = 1
            EMCOLD = EMCGEO
            DEPRED = DEPGEO
            CALL SIRSTP(JWSTEP,DUMMY,DUMMY,DUMMY,DUMMY,1)
Chj         ... WLKREJ in gnrinf.h will be true if geom. step rejected.
         END IF

         LPV = NNASHX*NNASHX
         CALL MEMGET2('REAL','PV',KPV,LPV,WRK,KFREE,LFREE)
         IF (NASHT .EQ. 1) THEN
Chj         PV(1) is 0.0D0 for a single open shell RHF
            WRK(KPV) = D0
         ELSE IF (NASHT .GT. 1) THEN ! HSROHF
            ! we fill PV with a very large number, to make sure it gives
            ! problems if used later. As it is now, HSROHF is done
            ! completely in AO basis, thus PV (and FQ below) are not
            ! used. This filling makes it very clear that there are
            ! problems if someone programs use of PV for HSROHF without
            ! checking code.
            WRK(KPV:KPV+LPV-1) = -9.87654321D123
         END IF

         CALL MEMGET2('REAL','GORB',KGORB,NWOPH,WRK,KFREE,LFREE)
         CALL DZERO(WRK(KGORB),NWOPH)
Chj      ... set orbital gradient to zero (i.e. completely converged!)

         LFQ = NASHT*NORBT
         CALL MEMGET2('REAL','FQ',KFQ,LFQ,WRK,KFREE,LFREE)
         IF (NASHT.GT.0) THEN
C           1. put FV in DCAO:
            IF (HSROHF .AND. .NOT.HSROHF_save) THEN
            !  calculate correct FV for interface SIRIFC
            !  (full FV has not been calculated before when HSROHF true, as
            !  FAOSAV(:,INDEVC,2) contains minus the exchange part of FVAO)
            !  Standard FV is necessary for NASHT.eq.1 because many places
            !  in ABACUS and RESPONS the test is on NASHT.EQ.1 instead of IOPRHF.gt.0 /141123-hjaaj

               ! retrieve final DVAO in DCAO
               IF (NSYM.GT.1) THEN
                  CALL PKSYM1(DCAO(1+N2BASX),DAOSAV(1,2),
     &                        NBAS,NSYM,-1)
                  CALL DUNFLD(NBAST,DCAO(1+N2BASX),DCAO)
               ELSE
                  CALL DUNFLD(NBAST,DAOSAV(1,2),DCAO)
               ENDIF
               DCAO(1:N2BASX) = -DCAO(1:N2BASX)
               HSROHF = .FALSE.
               CALL FCK2AO(.TRUE.,FCAO(1+N2BASX),DCAO,WRK,KFREE,LFREE) ! calculate FVAO in FCAO(1+N2BASX)
               ! ( first parameter .TRUE. tells FCK2AO only one density matrix )
               CALL DGETSP(NBAST,FCAO(1+N2BASX),WRK(KFREE))
               CALL PKSYM1(WRK(KFREE),FCAO(1+N2BASX),NBAS,NSYM,1)
               CALL UTHUB(FCAO(1+N2BASX),DCAO,CMO,WRK(KFREE),          ! and now FV is in DCAO
     &                    NSYM,NBAS,NORB)

            ELSE
               CALL UTHUB(FAOSAV(1,INDEVC,2),DCAO,CMO,WRK(KFREE),
     &                    NSYM,NBAS,NORB)
            END IF

C           2. form FC' = FD - FV' in FCAO, both for standard FV'=FV and for HSROHF FV'=-FV_exch
            CALL DAXPY(NNORBT,DM1,DCAO,1,FCAO,1)

C           3. FQ matrix is zero matrix for open shell RHF
            CALL DZERO(WRK(KFQ),NASHT*NORBT)
         END IF  ! (NASHT.GT.0)
C
C *** LOCALIZATION
C
         IF (BOYORB .OR. PIPORB) CALL SIRLOC(CMO,WRK(KFREE),LFREE)
         IF (BOYSEL) CALL LOCCTL(CMO,WRK(KFREE),LFREE)
C
         IF (LMULBS .OR. R12CAL) THEN
C           Construct orthonormal auxiliary basis set (WK/UniKA/04-11-2002).
            CALL TIMER('START ',TIMSTR,TIMEND)
            CALL R12AUX(WRK(KFREE),LFREE)
            CALL TIMER('R12AUX',TIMSTR,TIMEND)
         END IF
C
C *** and now we are ready to call WR_SIRIFC to write SIRIFC
         CALL MEMGET2('REAL','DUM',KDUM,0,WRK,KFREE,LFREE)
         CALL WR_SIRIFC(CMO,DV,WRK(KPV),FCAO,DCAO,
     &               WRK(KFQ),WRK(KDUM),WRK(KDUM),WRK(KDUM),
     &               WRK(KFREE),LFREE,WRK(KGORB),TNLM,.FALSE.,WRK(KDUM))
C        CALL WR_SIRIFC(CMO,DV,PV,FC,FV,FQ,CREF,FCAC,H2AC,WRK,LFREE,
C    &               GORB,ERLM,ORBHES,XINDX)
C
         CALL MEMREL('DIIS_CTL.WR_SIRIFC',WRK,1,KPV,KFREE,LFREE)
      END IF
C
C        Transform integrals for abacus or response,
C        if converged and if not followed by higher level,
C        and if not Cholesky decomposition
C
      IF (ICONV .NE. 0 .AND. ABS(ITRFIN) .LE. 10 .AND.
     &    .NOT.(DOMP2 .OR. DOCI .OR. DOMC) .AND.
     &    .NOT. CHOINT) THEN
         IF (IPRI6 .GE. 0) WRITE (LUPRI,'(/A,I3)')
     &      ' >>> Transforming 2-el. integrals acc. to'//
     &      ' .FINAL TRANSFORMATION =',ITRFIN
         JTRLVL = ITRFIN
         CALL TRACTL(JTRLVL,CMO,WRK(KFREE),LFREE)
      END IF
C
C ***********************************************************
C
 1100 CONTINUE
C     restore original values
      IFTHRS = IFTHSV
      HSROHF = HSROHF_save
      CALL QEXIT('DIIS_CTL')
      RETURN
      END
C  /* Deck DIIS_EVC */
      SUBROUTINE DIIS_EVC(SMOAO,DCAO,DVAO,H1AO,FDAO,FVAO,CMO1,
     &                    ERRVEC,WRK,KFRSAV,LFRSAV,IPR_DIIS)
C
C
#include "implicit.h"
      DIMENSION DCAO(NNBAST), DVAO(NNBAST), FDAO(NNBAST), FVAO(NNBAST)
      DIMENSION H1AO(NNBAST)
      DIMENSION SMOAO(N2BAST), CMO1(NCMOT), ERRVEC(NNORBT), WRK(*)
C
C Used from common blocks:
C  INFINP : LNOROT,NOROT(), HSROHF
C  INFORB : N2BAST, NNBAST, NCMOT, NSYM, ..., NFROT, ...
C  INFDIM : N2BASM
C
#include "gnrinf.h"
#include "maxorb.h"
#include "priunit.h"
#include "infinp.h"
#include "inforb.h"
#include "infdim.h"
#include "infpri.h"
#include "scbrhf.h"
C
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DM1 = -1.0D0)
C
      CALL QENTER('DIIS_EVC')
      KFREE = KFRSAV
      LFREE = LFRSAV
      CALL MEMGET2('REAL','TMP1',KTMP1,N2BASM,WRK,KFREE,LFREE)
      CALL MEMGET2('REAL','TMP2',KTMP2,N2BASM,WRK,KFREE,LFREE)
      CALL MEMGET2('REAL','TMP3',KTMP3,N2BASM,WRK,KFREE,LFREE)
      IF (HSROHF) THEN
         CALL DAXPY(NNBAST,D1,DVAO,1,DCAO,1)
         CALL DSCAL(NNBAST,DM1,DVAO,1)
      END IF
      DO 800 ISYM = 1,NSYM
         NOCCI = NOCC(ISYM)
         NORBI = NORB(ISYM)
      IF (NOCCI .EQ. 0) THEN
         IF (NORBI .GT. 0) THEN
            CALL DZERO(ERRVEC(1+IIORB(ISYM)),NNORB(ISYM))
         END IF
         GO TO 800
      END IF
         NBASI = NBAS(ISYM)
         IF (IPR_DIIS .GE. 35) THEN
            WRITE (LUPRI,*) 'Debug output from DIIS_EVC, ISYM=',ISYM
            WRITE (LUPRI,*) 'NORBI,NBASI:',NORBI,NBASI
            WRITE (LUPRI,*) 'DIIS_EVC CMO1:'
            CALL OUTPUT(CMO1(1+ICMO(ISYM)),1,NBASI,1,NORBI,
     &                  NBASI,NORBI,-1,LUPRI)
            WRITE (LUPRI,*) 'DIIS_EVC SMOAO = 4 Ct SAO:'
            CALL OUTPUT(SMOAO(1+I2BAS(ISYM)),1,NORBI,1,NBASI,
     &                  NORBI,NBASI,-1,LUPRI)
         IF (IPR_DIIS .GT. 50) THEN
            CALL DGEMM('N','N',NORBI,NORBI,NBASI,1.D0,
     &                 SMOAO(1+I2BAS(ISYM)),NORBI,
     &                 CMO1(1+ICMO(ISYM)),NBASI,0.D0,
     &                 WRK(KTMP3),NORBI)
            WRITE (LUPRI,*) 'DIIS_EVC SMO = 4 Ct SAO C:'
            CALL OUTPUT(WRK(KTMP3),1,NORBI,1,NORBI,
     &                  NORBI,NORBI,-1,LUPRI)
         END IF
         END IF
         CALL DSPTSI(NBASI,FDAO(1+IIBAS(ISYM)),WRK(KTMP3))
         CALL DSPTSI(NBASI,H1AO(1+IIBAS(ISYM)),WRK(KTMP2))
         CALL DAXPY(N2BAS(ISYM),D1,WRK(KTMP2),1,WRK(KTMP3),1)
         IF (NISH(ISYM) .GT. 0) THEN
            CALL DUNFLD(NBASI,DCAO(1+IIBAS(ISYM)),WRK(KTMP2))
            CALL DGEMM('N','N',NBASI,NBASI,NBASI,1.D0,
     &                 WRK(KTMP2),NBASI,
     &                 WRK(KTMP3),NBASI,0.D0,
     &                 WRK(KTMP1),NBASI)
            IF (IPR_DIIS .GE. 35) THEN
               WRITE (LUPRI,*) 'DIIS_EVC DCAO unfolded:'
               CALL OUTPUT(WRK(KTMP2),1,NBASI,1,NBASI,
     &                     NBASI,NBASI,-1,LUPRI)
               WRITE (LUPRI,*) 'DIIS_EVC FDAO:'
               CALL OUTPUT(WRK(KTMP3),1,NBASI,1,NBASI,
     &                     NBASI,NBASI,-1,LUPRI)
               WRITE (LUPRI,*) 'DIIS_EVC TMP1 = DCAO FDAO:'
               CALL OUTPUT(WRK(KTMP1),1,NBASI,1,NBASI,
     &                     NBASI,NBASI,-1,LUPRI)
            END IF
         ELSE
            CALL DZERO(WRK(KTMP1),NBASI*NBASI)
         END IF
         IF (NASH(ISYM) .GT. 0) THEN
            CALL DSPTSI(NBASI,FVAO(1+IIBAS(ISYM)),WRK(KTMP2))
            CALL DAXPY(NBASI*NBASI,DM1,WRK(KTMP2),1,WRK(KTMP3),1)
            CALL DUNFLD(NBASI,DVAO(1+IIBAS(ISYM)),WRK(KTMP2))
            CALL DGEMM('N','N',NBASI,NBASI,NBASI,1.D0,
     &                 WRK(KTMP2),NBASI,
     &                 WRK(KTMP3),NBASI,1.D0,
     &                 WRK(KTMP1),NBASI)
            IF (IPR_DIIS .GE. 35) THEN
               WRITE (LUPRI,*) 'DIIS_EVC DVAO unfolded:'
               CALL OUTPUT(WRK(KTMP2),1,NBASI,1,NBASI,
     &                     NBASI,NBASI,-1,LUPRI)
               WRITE (LUPRI,*) 'DIIS_EVC FCAO:'
               CALL OUTPUT(WRK(KTMP3),1,NBASI,1,NBASI,
     &                     NBASI,NBASI,-1,LUPRI)
               WRITE (LUPRI,*) 'DIIS_EVC TMP1 = DCAO FDAO + DVAO FCAO:'
               CALL OUTPUT(WRK(KTMP1),1,NBASI,1,NBASI,
     &                     NBASI,NBASI,-1,LUPRI)
            END IF
         END IF
         CALL DGEMM('N','N',NORBI,NBASI,NBASI,1.D0,
     &              SMOAO(1+I2BAS(ISYM)),NORBI,
     &              WRK(KTMP1),NBASI,0.D0,
     &              WRK(KTMP2),NORBI)
         CALL DGEMM('N','N',NORBI,NORBI,NBASI,1.D0,
     &              WRK(KTMP2),NORBI,
     &              CMO1(1+ICMO(ISYM)),NBASI,0.D0,
     &              WRK(KTMP1),NORBI)
C
C        Zero rows and columns corresponding to frozen orbitals
C
         DO 220 I = 1,NFRO(ISYM)
            JOFF = KTMP1 - 1 + (I-1)*NORBI
            DO 210 J = 1,NORBI
               WRK(KTMP1-1+(I-1)*NORBI+J) = D0
               WRK(KTMP1-1+(J-1)*NORBI+I) = D0
  210       CONTINUE
  220    CONTINUE
         IF (LNOROT) THEN
            IORBI = IORB(ISYM)
            DO 320 I = 1,NORBI
            IF (NOROT(IORBI+I) .NE. 0) THEN
               DO 310 J = 1,NORBI
                  WRK(KTMP1-1+(I-1)*NORBI+J) = D0
                  WRK(KTMP1-1+(J-1)*NORBI+I) = D0
  310          CONTINUE
            END IF
  320       CONTINUE
         END IF
C
C        note: now TMP1 = 4 Fmo = 4 Ct SAO (DCAO FDAO + DVAO FCAO) C
C              (SMOAO contains "4 Ct SAO")
C        note: DGETAP gives 0.5 (TMP1 - TMP1^t) = 2 ( Fmo - Fmo^t) = gradient
C
         CALL DGETAP(NORBI,WRK(KTMP1),ERRVEC(1+IIORB(ISYM)))
C
         IF (IPR_DIIS .GE. 15) THEN
            WRITE (LUPRI,'(/A/A,I3/A)')
     &      ' (DIIS_EVC) 4 Fmo = 4 (DCmo FDmo + DVmo FCmo)',
     &      '        = 4 Ct SAO (DCAO FDAO + DVAO FCAO) C for symmetry',
     &      ISYM,
     &      ' (error vector = gradient = 2 (Fmo - Fmo(transposed))'
            CALL OUTPUT(WRK(KTMP1),1,NORBI,1,NORBI,NORBI,NORBI,-1,LUPRI)
         END IF
  800 CONTINUE
      IF (HSROHF) THEN
         CALL DAXPY(NNBAST,D1,DVAO,1,DCAO,1)
         CALL DSCAL(NNBAST,DM1,DVAO,1)
      END IF
      CALL MEMREL('DIIS_EVC',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('DIIS_EVC')
      RETURN
      END
C  /* Deck DIIS_RED */
      SUBROUTINE DIIS_RED (C2DIIS,BMAT,INDEVC,NEVC,DAMP,THREVC,
     &                     CVEC,XLMBDA,ERRVEC,L_ERRVEC,
     &                     WRK,KFRSAV,LFRSAV,IPR_DIIS)
C
C 19-May-1993 HJAaJ+HA (based on KAPRED)
C
C Input:
C  C2DIIS, to select C1-DIIS or C2-DIIS:
C        FALSE  SOLVE LEVEL SHIFTED LINEAR SET OF EQUATIONS
C               (LEVEL SHIFT DAMP) TO FIND IMPROVED ORBITAL PARAMETERS
C        TRUE   SOLVE LEVEL SHIFTED LINEAR SET OF EQUATIONS
C               AS AN EIGENVALUE PROBLEM. ADJUST LEVEL SHIFT
C               TO OBTAIN STEP LENGTH RTRUST
C  BMAT,   the B matrix
C  ERRVEC, the NEVC error vectors
C  NEVC,   number of error vectors in ERRVEC
C  THREVC, threshold for acceptable solution in C2DIIS
C  DAMP,   damping factor IN LINEAR SET OF EQUATIONS
C
C Output:
C  BMAT,  the new, extended reduced orbital Hessian-matrix
C  CVEC  :.not.C2DIIS - SOLUTION TO LINEAR SET OF EQUATIONS
C        :     C2DIIS - LOWEST acceptable EIGENVECTOR
C
#include "implicit.h"
#include "priunit.h"
      DIMENSION BMAT(*), CVEC(*), ERRVEC(L_ERRVEC,*), WRK(*)
      LOGICAL   C2DIIS
C
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
C
      IROW(I) = (I*(I-1))/2
C
      CALL QENTER('DIIS_RED')
      KFREE = KFRSAV
      LFREE = LFRSAV
C
      IF (IPR_DIIS .GE. 5) WRITE (LUPRI,7555) C2DIIS,NEVC
 7555 FORMAT(/' (DIIS_RED) Control Parameters C2DIIS =',L10,
     &        ' NEVC =',I5)
C
C
C Section 1: extend BMAT with a new row corresponding to
C            the new error vector
C
      K = INDEVC
      KROW = IROW(K)
      DO, L = 1,K
         BMAT(KROW+L) = DDOT(L_ERRVEC,ERRVEC(1,L),1,ERRVEC(1,K),1)
      END DO
      DO, L = K+1,NEVC
         BMAT(IROW(L)+K) = DDOT(L_ERRVEC,ERRVEC(1,L),1,ERRVEC(1,K),1)
      END DO
      IF ( IPR_DIIS.GE.5 ) THEN
         WRITE (LUPRI,'(//A)') ' (DIIS_RED) B matrix:'
         CALL OUTPAK(BMAT,NEVC,-1,LUPRI)
      END IF
C
      LBMAT  = IROW(NEVC+1)
      IF (C2DIIS) THEN
         LBAUG  = LBMAT
      ELSE
         LBAUG  = IROW(NEVC+2)
      END IF
      CALL MEMGET2('REAL','BTMP',KBTMP,LBAUG,WRK,KFREE,LFREE)
      CALL DCOPY(LBMAT,BMAT,1,WRK(KBTMP),1)
      IF (.NOT. C2DIIS) GO TO 2000
C
C *******************************************************
C
C C2DIIS: SOLVE DIIS AS AN EIGENVALUE EQUATION
C
C
      CALL MEMGET2('REAL','EVEC',KEVEC,NEVC*NEVC,WRK,KFREE,LFREE)
      CALL MEMGET2('REAL','WJ'  ,KWJ  ,NEVC     ,WRK,KFREE,LFREE)
      CALL MEMGET2('INTE','IWJ' ,KIWJ ,NEVC     ,WRK,KFREE,LFREE)
C
      CALL DUNIT(WRK(KEVEC),NEVC)
      CALL JACO_THR(WRK(KBTMP),WRK(KEVEC),NEVC,NEVC,NEVC,1.0D-12)
      DO 150 I = 1,NEVC
         II = KBTMP-1+IROW(I+1)
         WRK(KBTMP-1+I) = WRK(II)
 150  CONTINUE
      CALL ORDER (WRK(KEVEC),WRK(KBTMP),NEVC,NEVC)
      IOK = 0
      DO 170 I = 1,NEVC
         EVCSUM = DSUM(NEVC,WRK(KEVEC+(I-1)*NEVC),1)
         IF (ABS(EVCSUM) .GT. THREVC) THEN
            IOK = I
            XLMBDA = WRK(KBTMP-1+I)
            CALL DCOPY(NEVC,WRK(KEVEC+(I-1)*NEVC),1,CVEC,1)
            CALL DSCAL(NEVC,(D1/EVCSUM),CVEC,1)
            GO TO 171
         END IF
 170  CONTINUE
 171  CONTINUE
      IF ( IPR_DIIS.GE.3 .OR. IOK.NE.1) THEN
         WRITE(LUPRI,'(//A,I5)')
     &      ' (DIIS_RED) C2-DIIS B matrix eigenvalues. IOK =',IOK
         CALL OUTPUT(WRK(KBTMP),1,1,1,NEVC,1,NEVC,-1,LUPRI)
      END IF
      IF ( IPR_DIIS.GE.4) THEN
         WRITE(LUPRI,'(/A)') ' - and B matrix eigenvectors:'
         CALL OUTPUT(WRK(KEVEC),1,NEVC,1,NEVC,NEVC,NEVC,-1,LUPRI)
      END IF
C
      GO TO 9999
C
C ********************************************
C
C DIIS: SOLVE LEVEL SHIFTED linear EQUATIONS
C
 2000 CONTINUE
C
      NLEQ  = NEVC + 1
      CALL MEMGET2('INTE','IPVT',KIPVT,NLEQ,WRK,KFREE,LFREE)
C
      IF (DAMP .NE. D0) THEN
         CALL QUIT('DIIS_RED: DAMP not implemented yet')
#if defined (VAR_KAPREDCODE)
  >>>>>>>>>>>>>>>> NOTE: Hamilton and Pulay uses '* (D1 + DAMP)'
                                        and not  ' + DAMP'
         DO 510 I = 1,NEVC
            II = KBTMP - 1 + IROW(I+1)
            WRK(II) = WRK(II) + DAMP
 510     CONTINUE
         IF ( IPR_DIIS.GE.3 ) THEN
            WRITE(LUPRI,'(A,1P,D15.8)')
     &      ' (DIIS_RED) B matrix damping factor =',DAMP
         END IF
#endif
      END IF
      DO 520 I = 1,NEVC
         WRK(KBTMP-1+LBMAT+I) = D1
  520 CONTINUE
      WRK(KBTMP-1+LBAUG) = D0
      CALL DZERO(CVEC,NEVC)
      CVEC(NLEQ) = D1
      CALL DSPSOL(NLEQ,1,WRK(KBTMP),CVEC,WRK(KIPVT),INFO)
      XLMBDA = CVEC(NLEQ)
      IF ( IPR_DIIS.GE.7 ) THEN
         WRITE(LUPRI,'(/A)')
     &      ' (DIIS_RED) Solutions to C1-DIIS set of linear equations'
         WRITE (LUPRI,'(10F12.6)') (CVEC(I),I=1,NLEQ)
      END IF
      IF (INFO.NE.0) THEN
         WRITE (LUPRI,8500) INFO
         WRITE (LUPRI,'(//A)') ' (DIIS_RED) B matrix:'
         CALL OUTPAK(BMAT,NEVC,-1,LUPRI)
         CALL QTRACE(LUPRI)
         CALL QUIT('DIIS_RED: no solution to linear equations')
      END IF
 8500 FORMAT(/' (DIIS_RED) Solution not obtained to linear equations'
     &       /T11,'Check if matrix is singular, LINPACK DSP code =',I3)
C
C *** End of subroutine DIIS_RED
C
 9999 CONTINUE
      IF (IPR_DIIS .GT. 1) THEN
         WRITE (LUPRI,'(/A,L2,A,1P,D10.2)')
     &      ' DIIS C vector; C2DIIS =',C2DIIS,', LAMBDA =',XLMBDA
         WRITE (LUPRI,'(10F12.6)') (CVEC(I),I=1,NEVC)
      END IF
      CALL MEMREL('DIIS_RED',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('DIIS_RED')
      RETURN
      END
c /* DFT_ADD_KS */
      SUBROUTINE DFT_ADD_KS(ONLYFC,EMY,DCAO,NDMAT,FCAO,FVAO,
     &                      WRK,LWRK,IPRFCK)
c     TMP1 is of NNBASX size, TMP2 is of N2BASX size.  They are passed
c     because they are available but they could as well be allocated
c     here.
#include "implicit.h"
c
#include "inforb.h"
#include "dfterg.h"
#include "priunit.h"
c
      PARAMETER (DP5 = 0.5D0, D1 = 1.0D0)
c
      DIMENSION DCAO(N2BASX*NDMAT),FCAO(N2BASX),FVAO(N2BASX)
      DIMENSION WRK(LWRK)
      LOGICAL ONLYFC
c
      CALL QENTER('DFT_ADD_KS')
c
      KTMP1 = 1
      KTMP2 = KTMP1 + NNBASX
      KLST  = KTMP2 + N2BASX
      LFREE = LWRK  - KLST +1
      IF(LFREE.LT.0) CALL STOPIT('ADDKS1','SIRDIIS',KLST,LFREE)
C     EXCTRO = DP5*DDOT(NNBASX,DCAO,1,FCAO,1)
      IF(ONLYFC) THEN
         CALL DZERO(WRK(KTMP2),N2BASX)
         CALL DFTKSMb(DCAO,WRK(KTMP2),EDFTY,WRK(KLST),LFREE,IPRFCK)
         CALL DGETSP(NBAST,WRK(KTMP2),WRK(KTMP1))
         IF (NSYM .GT. 1) THEN
            KKOHN = KLST
            KLST  = KKOHN + NNBAST
            IF (KLST.GT.LWRK) CALL STOPIT('ADDKS','SIRDIIS',KLST,LWRK)
            CALL PKSYM1(WRK(KTMP1),WRK(KKOHN),NBAS,NSYM,1)
            CALL DAXPY(NNBAST,1D0,WRK(KKOHN),1,FCAO,1)
         ELSE
            CALL DAXPY(NNBAST,1D0,WRK(KTMP1),1,FCAO,1)
         END IF
      ELSE
C     memory allocation
         KDENA = KLST
         KKSMA = KDENA +N2BASX*2
         KKSPA = KKSMA +N2BASX*2
         KKSPB = KKSPA +NNBAST
         KFREE = KKSPB +NNBAST
         IF (KFREE.GT.LFREE) CALL STOPIT('ADDKS1','SIRDIIS',KFREE,LFREE)
         CALL DZERO(WRK(KDENA),N2BASX*2)
         CALL DZERO(WRK(KKSMA),N2BASX*2)
         CALL DZERO(WRK(KKSPA),NNBAST)
         CALL DZERO(WRK(KKSPB),NNBAST)
C     Kohn-Sham contribution evaluation
         IF (NISHT .GT. 0) THEN
            CALL DAXPY(N2BASX,DP5,DCAO,1,WRK(KDENA),1)
            CALL DAXPY(N2BASX,DP5,DCAO,1,WRK(KDENA+N2BASX),1)
         END IF
         CALL DAXPY(N2BASX,1.D0,DCAO(1+N2BASX),1,WRK(KDENA),1)
C
C  For specific case of molecules bearing only alpha electron(s),
C  like hydrogen atom, we call "old" open-shell code, which
C  can handle zero beta density.
C
         IF (NISHT .EQ. 0) THEN
            CALL dft_kohn_shamab(WRK(KDENA),WRK(KKSMA),EDFTY,
     &                           WRK(KFREE),LFREE,IPRFCK)
         ELSE
            CALL dft_kohn_shamab_b(WRK(KDENA),WRK(KKSMA),EDFTY,
     &                             WRK(KFREE),LFREE,IPRFCK)
         END IF
C     packing
         CALL DGETSP(NBAST,WRK(KKSMA),WRK(KTMP1))
         IF (NSYM .GT. 1) THEN
            CALL PKSYM1(WRK(KTMP1),WRK(KKSPA),NBAS,NSYM,1)
         ELSE
            CALL DAXPY(NNBAST,D1,WRK(KTMP1),1,WRK(KKSPA),1)
         END IF
         CALL DGETSP(NBAST,WRK(KKSMA+N2BASX),WRK(KTMP1))
         IF (NSYM .GT. 1) THEN
            CALL PKSYM1(WRK(KTMP1),WRK(KKSPB),NBAS,NSYM,1)
         ELSE
            CALL DAXPY(NNBAST,D1,WRK(KTMP1),1,WRK(KKSPB),1)
         END IF
c     WRITE(LUPRI,'(/A)') 'ksma: '
c     CALL OUTPKB(WRK(KKSPA),NORB,NSYM,-1,LUPRI)
c     WRITE(LUPRI,'(/A)') 'ksmB: '
c     CALL OUTPKB(WRK(KKSPB),NORB,NSYM,-1,LUPRI)
c
         CALL DAXPY(NNBAST, DP5,WRK(KKSPA),1,FCAO,1)
         CALL DAXPY(NNBAST, DP5,WRK(KKSPB),1,FCAO,1)
c
         CALL DAXPY(NNBAST,-DP5,WRK(KKSPA),1,FVAO,1)
         CALL DAXPY(NNBAST, DP5,WRK(KKSPB),1,FVAO,1)
      END IF
C
      EMY = EMY + EDFTY
C     EDFTY = EDFTY +EXCTRO-DP5*DDOT(N2BASX,DCAO,1,FCAO,1)
      CALL QEXIT('DFT_ADD_KS')
      RETURN
      END
C  /* Deck fckeig */
      SUBROUTINE FCKEIG(CMO,FC,SHFTLVL,THR_FCKEIG,NCAN,IPRINT_in,
     &                  SCRA,LSCRA)
C
C Written 20-May-1993 by Hans Jorgen Aa. Jensen
C
C Purpose:
C  If requested (by SHFTLVL .ne. 0) :
C        Perform "restricted step" step reduction (i.e. level shift)
C        by adding SHFTLVL*DENS(i,i) to all diagonal elements;
C  Diagonalize Fock matrix
C
C Input:
C  CMO; initial molecular orbitals used to build Fock matrix,
C       assumed to be orthonormal.
C  FC;  the inactive Fock matrix
C  NCAN(8); number of orbitals in each symmetry to  make canonical
C
C Output:
C  CMO; molecular orbitals diagonalizing Fock matrix
C  FC;  the orbital energies
C
C Scratch:
C  SCRA; general scratch area
C
#include "implicit.h"
      DIMENSION CMO(*),FC(*),NCAN(8),SCRA(*)
C
C Used from common blocks:
C  INFINP : SUPSYM, LNOROT, NOROT()
C  SCBRHF : NFRRHF(*), AUTOCC
C  INFIND : IROW(*),...,ISSMO(*),?
C
#include "maxash.h"
#include "maxorb.h"
#include "priunit.h"
#include "infinp.h"
#include "inforb.h"
#include "scbrhf.h"
#include "infind.h"
#include "infpri.h"
#include "mxcent.h"
#include "dftcom.h"
#include "dftacb.h"
#include "dummy.h"
C
      LOGICAL LSAVE4,LSAVE6
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0)
C
      CALL QENTER('FCKEIG')
      call gettim(t1,w1)
C
      IPR_FCKEIG = IPRINT_in
      LSAVE4 = P4FLAG(9)
      LSAVE6 = P6FLAG(6)
      P4FLAG(9) = .FALSE.
      P6FLAG(6) = .FALSE.
C
C     Some memory allocation for SCF occupation determination
C
      KEIG  = 1
      KSYMS = KEIG  + NORBT
      KLAST = KSYMS + NORBT
C
C     Step 1: Diagonalize (level-shifted) Fock-matrix:
C
      K = 0
      DO, ISYM = 1,NSYM
         K = K + NOCC(ISYM)*(NORB(ISYM)-NOCC(ISYM))
      END DO
      FAC_THR_JACO = 16*K
      ! factor 16 because gradient ia = 2 F_ia = 4 FC_ia
      FAC_THR_JACO = 1.0D0 / SQRT(FAC_THR_JACO)
      ! factor to ensure gradient is less than THR_FCKEIG if all
      ! occ-vir elements are less than THR_FCKEIG/FAC_THR_JACO /hjaaj

      DO 200 ISYM = 1,NSYM
         NCANI = NCAN(ISYM)
      IF ( .NOT.AUTOCC .AND. NCANI.EQ.0 ) GO TO 200
         NORBI = NORB(ISYM)
      IF (NORBI.EQ.0) GO TO 200
         IORBI = IORB(ISYM)
         NBASI = NBAS(ISYM)
         NFRZI = NFRRHF(ISYM)
         IFSYM = IIORB(ISYM)
         ICSYM = ICMO(ISYM) + 1
C
C        If requested (by SHFTLVL .ne. 0) :
C        Perform "restricted step" step reduction (i.e. level shift)
C        by adding 0.5*SHFTLVL*DENS(i,i) to all diagonal elements:
C        /hjaaj March 2001
C
         IF (SHFTLVL .NE. D0) THEN
            DO K = 1,NISH(ISYM)
               FC(IFSYM+IROW(K+1)) = FC(IFSYM+IROW(K+1)) + SHFTLVL
            END DO
            DO K = NISH(ISYM)+1, NISH(ISYM)+NASH(ISYM)
               FC(IFSYM+IROW(K+1)) = FC(IFSYM+IROW(K+1)) + DP5*SHFTLVL
            END DO
         END IF
C        zero rows and columns corresponding to frozen orbitals
         JFRZI = 0
         IF (NFRZI .GT. 0 .OR. LNOROT) THEN
            JNFRZ = 0
            DO 130 K = 1,NORBI
               IK=IORBI+K
               IF (K .LE. NFRZI .OR. NOROT(IK) .NE. 0) THEN
                  IF (IOBTYP(IK).NE.JTACT) THEN
C              ... enable freezing an active orbital
                     IF (JNFRZ .GT. 0) GO TO 9000
C              ... exit if free orbital below this frozen orbital
C                  because then ORDRSS cannot be called
                     JFRZI = JFRZI + 1
                  END IF
                  KROW  = IROW(K)
                  DO 110 L = 1,K-1
                     FC(IFSYM+KROW+L) = D0
  110             CONTINUE
                  DO 120 L = K+1,NORBI
                     KL = (IFSYM+K)+IROW(L)
                     FC(KL) = D0
  120             CONTINUE
                  IF (NOROT(IK).NE.0 .AND. IOBTYP(IK).EQ.JTACT) THEN
                     KK = (IFSYM+K)+IROW(K)
                     FC(KK) = DUMMY
C
C This assignment puts it last in the sorting  - put back in place later
C
                  END IF
               ELSE
                  JNFRZ = JNFRZ + 1
               END IF
  130       CONTINUE
         END IF
C
         IF (NCANI .gt. 0) THEN
            THR_JACO = THR_FCKEIG * FAC_THR_JACO
            THR_JACO = MIN( THR_JACO, 1.D-3)
         ELSE IF (AUTOCC) THEN
            THR_JACO = 1.D-3
         ELSE
            GO TO 200
         END IF

         NMAX = MIN(NORBI, NCANI + 5)
         CALL JACO_THR(FC(IFSYM+1),CMO(ICSYM),NORBI,NMAX,NBASI,THR_JACO)
C        CALL JACO_THR(F,VEC,NB,NMAX,NROWV,THR_CONV)
C
         IF (IPR_FCKEIG .GT. 20) THEN
            write (lupri,*) 'FCKEIG debug: FC after JACO_THR',NMAX,NORBI
            call outpak(FC(IFSYM+1),NORBI,-1,LUPRI)
         END IF
         DO, I=1,NORBI
            FC(IORBI+I) = FC(IFSYM+IROW(I+1))
         END DO
C
         IF (AUTOCC) THEN
            CALL DCOPY(NORBI,FC(IORBI + 1),1,SCRA(KEIG + IORBI),1)
            DO, IK = 1, NORBI
               SCRA(KSYMS + IORBI + IK - 1) = ISYM
            END DO
         END IF
C
C
         JCSYM = ICSYM + JFRZI*NBASI
         JORBI = IORBI + JFRZI + 1
         NNOTFR= NORBI - JFRZI
         CALL ORDRSS(CMO(JCSYM),FC(JORBI),ISSMO(JORBI),NNOTFR,NBASI)
C
C         CALL HEADER('Sorted orbitals',-1)
C         CALL PRORB(CMO,.FALSE.,LUPRI)
C
C If there are frozen open shell orbitals they are last - put back in place
C
         IF (LNOROT) THEN
            NISHI=NISH(ISYM)
            NASHI=NASH(ISYM)
            JA=NISHI+1
            KA=IORBI+NISHI+1
            DO IA=1,NASHI
               JA=IA+NISHI
               KA=JA+IORBI
               IF (NOROT(KA).NE.0) THEN
                  CALL ROTCOL(CMO(ICSYM),NBASI,NORBI,JA,NORBI)
                  CALL ROTCOL(FC(IORBI+1),1,NORBI,JA,NORBI)
               END IF
            END DO
         END IF
C
C         CALL HEADER('Reordered orbitals',-1)
C         CALL PRORB(CMO,.FALSE.,LUPRI)
         IF (IPR_FCKEIG .GE. 4) WRITE (LUPRI,'(/A,I2/,(5F15.5))')
     &      ' DIIS Fock eigenvalues in symmetry',ISYM,
     &      (FC(IORBI+I), I = 1,NORBI)


  200 CONTINUE
      IF (SUPSYM) CALL AVEORD()
C     ... remake ISSORD() as ISSMO() may have changed in ORDRSS
C
C     Step 2: Reorthogonalize new mo's
C
      KSAO  = KLAST
      KSCR1 = KSAO + NNBAST
      LSCR1 = LSCRA - KSCR1

      IF (DFTASC) CALL GET_HOMO_ASC(FC,SCRA(KSCR1),LSCR1)

      CALL ORTHO(CMO,SCRA(KSAO),SCRA(KSCR1),LSCR1)
      IF (SUPSYM) THEN
         KFREE = 1
         LFREE = LSCRA
         CALL AVECPH(IPHCHA,CMO,SCRA(KLAST),KFREE,LFREE)
      END IF
C
C     Step 3: Reorder Hartree-Fock occupation just in case last
C             occupation suggestion was wrong
C
      IF (AUTOCC) THEN
         CALL ORDER(SCRA(KSYMS),SCRA(KEIG),NORBT,1)
         CALL IZERO(NISH,8)
         MOCC = NRHFEL/2
         DO 98 IK = 1, MOCC
            ISYM = NINT(SCRA(KSYMS + IK -1))
            NISH(ISYM) = NISH(ISYM) + 1
 98      CONTINUE
         CALL ICOPY(8,NISH,1,NOCC,1)
         IF (2*MOCC .NE. NRHFEL) THEN
            IOPRHF = NINT(SCRA(KSYMS + MOCC))
            LSYM   = IOPRHF
            CALL IZERO(NASH,8)
            CALL IZERO(IASH,8)
            NASH(IOPRHF) = 1
            NOCC(IOPRHF) = NOCC(IOPRHF) + 1
            DO 96 ISYM = IOPRHF + 1, 8
               IASH(ISYM) = 1
 96         CONTINUE
         END IF
      END IF
C
C *** end of subroutine FCKEIG
C
      IF (IPR_FCKEIG > 0) THEN
         call gettim(t2,w2)
         write(lupri,'(A,1P,3D10.2)')
     &      ' Time used in FCKEIG, thr_jaco: ',t2-t1,w2-w1,thr_jaco
      END IF
      CALL QEXIT('FCKEIG')
      RETURN
C
 9000 CONTINUE
      WRITE (LUPRI,'(//A/A)')
     &  ' ERROR in FCKEIG, FCKEIG cannot handle ".FREEZE" for orbitals',
     &  ' unless they also could have been frozen with ".FROZEN"'
      CALL QTRACE(LUPRI)
      CALL QUIT('FCKEIG error, cannot handle orbitals '//
     &          'frozen with ".FREEZE"')

      CONTAINS

      SUBROUTINE ROTCOL(A,NR,NC,I1,I2)
C
C Cyclic permutation of columns in A
C I2>I1 is moved to I1 and rest shifted right
C
      IMPLICIT NONE
      DOUBLE PRECISION A,B
      INTEGER NR,NC,I1,I2
      DIMENSION A(NR,NC),B(NR)

      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: TMP
      INTEGER I,MC
      IF (I2.LE.I1) CALL QEXIT('ROTCOL ERROR: I1<I2 not fulfilled')
      IF (I1.LE.0 .OR. I2.LE.0 .OR. I1.GT.NC.OR.I2.GT.NC)
     &   CALL QEXIT('ROTCOL ERROR: I1 or I2 out of bounds')
      MC=(I2-I1+1)
      ALLOCATE(TMP(NR,MC))
      TMP(:,1)=A(:,I2)
      TMP(:,2:)=A(:,I1:I2-1)
      A(:,I1:I2)=TMP(:,:)
      DEALLOCATE(TMP)
      RETURN
      END SUBROUTINE ROTCOL

      END   ! subroutine FCKEIG

      SUBROUTINE FVCORR(NVEC,W,DAO,FAO,MXERRV,F,WRK,KFREE,LFREE)
      IMPLICIT NONE
C
C Form projected (effective) fock matrix in AO basis for each accumulated
C iteration and sum up with diis weights.
C
C Input: NVEC accumulated Fock (FAO) and density (DAO)
C        matrices from previous DIIS iterations
C
C Output: sum of weighted effective fock matrices to be diagonalized
C
C
      INTEGER NVEC, MXERRV, KFREE, LFREE
      REAL*8 W(NVEC)
#include "inforb.h"
#include "priunit.h"
C     REAL*8 DAO(NNBAST,MXERRV,2), FAO(NNBAST,MXERRV,2), F(NNBAST)
      REAL*8 DAO(NNBAST,1,2),      FAO(NNBAST,MXERRV,2), F(NNBAST)
      REAL*8 WRK(*)
C
C Local
C
      REAL*8 D0, D1, D2
      PARAMETER (D0=0.0D0, D1=1.0D0, D2=2.0D0)
      INTEGER ISYM, KS, KTMP1, KTMP2, KDI, KFI, KSI
      INTEGER NBASI, NNBASI, IIBASI, N2BASI, IBLK, IVEC, I, II
      INTEGER IMXVEC
      EXTERNAL IMXVEC
      LOGICAL FOUND
C
      CALL QENTER('FVCORR')

      CALL MEMGET2('REAL','S',KS,NNBAST,WRK,KFREE,LFREE)
      N2BASI=IMXVEC(N2BAS,NSYM)
      CALL MEMGET2('REAL','DI',KDI,N2BASI,WRK,KFREE,LFREE)
      CALL MEMGET2('REAL','SI',KSI,N2BASI,WRK,KFREE,LFREE)
      CALL MEMGET2('REAL','FI',KFI,N2BASI,WRK,KFREE,LFREE)
      CALL MEMGET2('REAL','TMP1',KTMP1,N2BASI,WRK,KFREE,LFREE)
      CALL MEMGET2('REAL','TMP2',KTMP2,N2BASI,WRK,KFREE,LFREE)
      FOUND = .FALSE.
      CALL RDONEL('OVERLAP',FOUND,WRK(KS),NNBAST)
      IF (.NOT.FOUND) THEN
         CALL QUIT('FVCORR:Error reading overlap')
      END IF
      DO ISYM=1,NSYM
         NBASI=NBAS(ISYM)
         NNBASI=NNBAS(ISYM)
         N2BASI=N2BAS(ISYM)
         IIBASI=IIBAS(ISYM)
         IBLK=IIBAS(ISYM)+1
         IF (NBASI.EQ.0) GO TO 100
         CALL DSPTSI(NBASI,WRK(KS+IIBASI),WRK(KSI))
         DO IVEC=1,NVEC
            CALL DUNFLD(NBASI,DAO(IBLK,IVEC,2),WRK(KDI))
         ! -s*(-do)
            CALL DGEMM(
     &         'N','N',NBASI,NBASI,NBASI,
     &        -D1,WRK(KSI),NBASI,
     &            WRK(KDI),NBASI,
     &         D0,WRK(KTMP1),NBASI
     &         )
      ! (s*do)*(fc-fo)
            CALL DSPTSI(NBASI,FAO(IBLK,IVEC,2),WRK(KFI))
            CALL DGEMM(
     &         'N','N',NBASI,NBASI,NBASI,
     &         D1,WRK(KTMP1),NBASI,
     &            WRK(KFI),NBASI,
     &         D0,WRK(KTMP2),NBASI
     &         )
        ! (do+dc)*s
            !CALL DUNFLD(NBASI,DAO(IBLK,IVEC,1),WRK(KTMP1))
            !CALL DAXPY(N2BASI,D1,WRK(KTMP1),1,WRK(KDI),1)
            CALL DUNFLD(NBASI,DAO(IBLK,IVEC,1),WRK(KDI))
            CALL DGEMM(
     &         'N','N',NBASI,NBASI,NBASI,
     &         D1,WRK(KDI),NBASI,
     &            WRK(KSI),NBASI,
     &         D0,WRK(KTMP1),NBASI
     &         )
        ! (do+dc)*s-1
            DO I=1,NBASI
               II=KTMP1+(NBASI+1)*(I-1)
               WRK(II) = WRK(II) - D1
            END DO
        !final correction
        ! (s*do)*(fc-fo)*((do+dc)*s-1)
            CALL DGEMM(
     &         'N','N',NBASI,NBASI,NBASI,
     &         D1,WRK(KTMP2),NBASI,
     &            WRK(KTMP1),NBASI,
     &         D0,WRK(KDI),NBASI
     &         )
            CALL DGETSP(NBASI,WRK(KDI),WRK(KFI))
            CALL DAXPY(NNBASI,D2*W(IVEC),WRK(KFI),1,F(IBLK),1)
         END DO
 100     CONTINUE
      END DO
      CALL MEMREL('FVCORR',WRK,KS,KS,KFREE,LFREE)
      CALL QEXIT('FVCORR')
      END

      SUBROUTINE GET_HOMO_ASC(FD,SCRA,LSCRA)

#include "implicit.h"
#include "priunit.h"
      DIMENSION FD(*),SCRA(*)
      PARAMETER (DBIG = 1.D+12)
#include "maxash.h"
#include "maxorb.h"
#include "mxcent.h"
#include "inforb.h"
#include "infinp.h"
#include "infind.h"
#include "dftacb.h"
#include "dftcom.h"

      logical   exchanged
C
      KORBEN = 1
      LNEED  = KORBEN + NORBT
      IF (LNEED .GT. LSCRA) CALL ERRWRK('GET_HOMO_ASC',LNEED,LSCRA)
      CALL DZERO(SCRA,NORBT)

      call dcopy(NORBT,FD,1,SCRA,1)
100   continue
      exchanged=.false.
      do i=1,norbt-1
        if (scra(i+1) .lt. scra(i) ) then
          tmp=scra(i+1)
          scra(i+1)=scra(i)
          scra(i)=tmp
          exchanged=.true.
        end if
      end do
      if (exchanged) goto 100

      EhomoA =  scra(NOCCT)
      ehomo  =  ehomoa
      ehomob =  scra(NishT)

      ElumoA =  scra(NOCCT+1)
      ElumoB =  scra(NishT+1)

      SHIFA = DFTIPTA + EHOMOA - SHF

      write(LUPRI,10)' E_HOMO Alpha ',EhomoA,' E_HOMO Beta ',
     &                EhomoB, 'v_xc(inf)', SHIFA
      write(LUPRI,10)' E_LUMO Alpha ',ElumoA,' E_LUMO Beta ',
     &                ElumoB
10    format(5x,A14,F20.12,A14,F20.12,A14,F20.12)

      RETURN
      END

      SUBROUTINE SDFCOMM(F,D,SAOUNP,COMM,KFRSAV,LFRSAV,WRK)
C     In: F apk,D apk and S unp
C     Out: SDF-FDS

#include "implicit.h"
#include "inforb.h"
      DOUBLE PRECISION F(NNBAST),D(NNBAST),SAOUNP(NBAST,NBAST)
      DOUBLE PRECISION COMM(NBAST,NBAST),WRK(*)
      PARAMETER (D1=1.0D0, D0=0.0D0, DM1=-1.0D0)

      KFREE = KFRSAV
      LFREE = LFRSAV
      CALL MEMGET('REAL',KFTMP2,N2BASX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KDTMP2,N2BASX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KDForFD,N2BASX,WRK,KFREE,LFREE)

      IF (NSYM.GT.1) THEN
         CALL MEMGET('REAL',KFTMP,NNBASX,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KDTMP,NNBASX,WRK,KFREE,LFREE)
      ENDIF

      IF (NSYM.GT.1) THEN
         CALL PKSYM1(WRK(KFTMP),F,NBAS,NSYM,-1)
         CALL DSPTSI(NBAST,WRK(KFTMP),WRK(KFTMP2))
         CALL PKSYM1(WRK(KDTMP),D,NBAS,NSYM,-1)
         CALL DUNFLD(NBAST,WRK(KDTMP),WRK(KDTMP2))
      ELSE
         CALL DSPTSI(NBAST,F,WRK(KFTMP2))
         CALL DUNFLD(NBAST,D,WRK(KDTMP2))
      ENDIF

      CALL DGEMM('N','N',NBAST,NBAST,NBAST,D1,WRK(KFTMP2),NBAST,
     &     WRK(KDTMP2),NBAST,D0,WRK(KDForFD),NBAST)
      CALL DGEMM('N','N',NBAST,NBAST,NBAST,D1,WRK(KDForFD),NBAST,
     &     SAOUNP,NBAST,D0,COMM,NBAST)
           CALL DGEMM('N','N',NBAST,NBAST,NBAST,D1,WRK(KDTMP2),NBAST,
     &     WRK(KFTMP2),NBAST,D0,WRK(KDForFD),NBAST)
      CALL DGEMM('N','N',NBAST,NBAST,NBAST,D1,SAOUNP,NBAST,
     &     WRK(KDForFD),NBAST,DM1,COMM,NBAST)

      CALL MEMREL('SDFCOMM',WRK,KFRSAV,KFTMP2,KFREE,LFREE)
      END
! -- end of sirdiis.F --
