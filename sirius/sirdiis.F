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
#include "single.h"
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
!Revision 1.8  2001/01/31 13:49:00  hjj
!allocate FQ matrix for open shell before WR_SIRIFC
!
!Revision 1.6  2000/05/24 12:28:27  hjj
!1) fixed solvent calculations (DIFDEN must be false for solvent)
!2) fixed WR_SIRIFC call for open shell and solvent
!   (required change in SOLFCK also)
!=============================================================================
!970326-hjaaj -- adjusted abort DIIS algorithm (allow more NEWOCC and more
!  oscillations)
!951130-hjaaj
!DIIS_CTL: moved 'reset to not converged' for solvent and writing SIRIFC to
!        DIIS_CTL from SIRCTL; it is a deficiency in the DIIS routines that
!        it cannot be done in them. Look for MAERKE.
!950824-hjaaj
!DIIS_CTL: output for AUTOCC
!FCKEIG: also AUTOCC for open shell plus smaller changes.
!9505-k.ruud
!c: implemented AUTOCC option
!940916-hjaaj
!DIIS_CTL: changed print of energy from F18.12 to F22.12
!  (to avoid stars for e.g. Pb)
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
         WRITE (LUPRI,'(/A/A)')
     &   ' SRDIIS INFO: no DIIS iterations because degenerate',
     &   ' SRDIIS INFO: supersymmetries not implemented for DIIS'
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
         MXERRV = MIN(NWOPT,10)
C        ... fixes linear dependency problem for very small debug runs
C            MXERRV=10 gave lin.dep. for 2 el. in 2 orbitals,
C            MXERRV=2 gave quadratic convergence ! (940701-hjaaj)
         MXERRV = MAX(2,MXERRV)
      END IF
      IPR_DIIS = IPRRHF
      C2DIIS = FLAG(22)
C
      LUWOUT = LUPRI
  10  WRITE (LUWOUT,'(//A)')
     &      ' >>>>> DIIS optimization of Hartree-Fock <<<<<'
Cholesky
      IF (CHOINT) THEN
         WRITE(LUWOUT,'(A)')
     &      '     (using Cholesky decomposed integrals)'
         IF (.NOT. CCMODSK) THEN
            WRITE(LUWOUT,'(A,1P,D9.2,A)')
     &      '     (decomposing density, thr = ',THRDC,')'
         END IF
      ENDIF
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
      CALL MEMGET('REAL',KSMOAO,N2BAST,WRK,KFREE,LFREE)
      KREL = KFREE
      CALL MEMGET('REAL',KSAO ,NNBAST,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KSTMP,N2BAST,WRK,KFREE,LFREE)
      CALL RDONEL('OVERLAP',.TRUE.,WRK(KSAO),NNBAST)
      IF (IPR_DIIS .GT. 25) THEN
         WRITE (LUPRI,'(/A)') ' SRDIIS test output of SAO:'
         CALL OUTPKB(WRK(KSAO),NBAS,NSYM,1,LUPRI)
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
      CALL MEMGET('REAL',KCMO1 ,NCMOT, WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KDCAO ,NDMAT*N2BASX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KFCAO ,NDMAT*N2BASX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KH1AO ,NNBAST,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KDSAV ,NDMAT*NNBAST,WRK,KFREE,LFREE)
      IF (FLAG(16)) THEN
         CALL MEMGET('REAL',KTNLM,2*NLMSOL,WRK,KFREE,LFREE)
      ELSE
         CALL MEMGET('REAL',KTNLM,0,WRK,KFREE,LFREE)
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
      LBMAT = (MXERRV+1)*MXERRV / 2
      CALL MEMGET('REAL',KBMAT ,LBMAT, WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KCVEC ,MXERRV+1,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KFSAV ,NDMAT*MXERRV*NNBAST,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KEVSAV,MXERRV*NNORBT,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KEGSAV,MXERRV*2     ,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KDV,NNASHX,WRK,KFREE,LFREE)
      CALL DCOPY(NCMOT,CMO,1,WRK(KCMO1),1)
      CALL DIIS_CTL(CMO,WRK(KCMO1),WRK(KBMAT),WRK(KH1AO),
     &            WRK(KFSAV),WRK(KEVSAV),WRK(KEGSAV),MXERRV,C2DIIS,
     &            WRK(KSMOAO),WRK(KDCAO),WRK(KFCAO),WRK(KDSAV),WRK(KDV),
     &            WRK(KCVEC),WRK(KTNLM),WRK,KFREE,LFREE,IPR_DIIS,ICONV)
C
 9990 CALL MEMREL('SRDIIS',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
 9999 CALL FLSHFO(LUW4)
      IF (LUPRI .NE. LUW4) CALL FLSHFO(LUPRI)
      CALL QEXIT('SRDIIS')
      RETURN
      END
C  /* Deck DIIS_CTL */
      SUBROUTINE DIIS_CTL(CMO,CMO1,BMAT,H1AO,FAOSAV,EVCSAV,EGSAV,MXERRV,
     *                   C2DIIS,SMOAO,DCAO,FCAO,DAOSAV,DV,CVEC,
     *                   TNLM,WRK,KFRSAV,LFRSAV,IPR_DIIS,ICONV)
C
C L.r. 5-May-1994 hjaaj / 1-Apr-1997
C
C DFT modifications T.  Helgaker
C
C
      use polarizable_embedding, only: peqm, pe_mep, pe_master,
     &                                 pe_savden, pe_save_density
#include "implicit.h"
#include "dummy.h"
      DIMENSION CMO(*), CMO1(*), BMAT(*), H1AO(*)
      DIMENSION FAOSAV(NNBAST,MXERRV,2), EVCSAV(NNORBT,MXERRV)
      DIMENSION SMOAO(*), DCAO(*), FCAO(*), DAOSAV(NNBAST,*), CVEC(*)
      DIMENSION DV(*)
      DIMENSION EGSAV(2,MXERRV), TNLM(*), WRK(*)

      real*8, dimension(1) :: pe_energy
      real*8, dimension(:), allocatable :: pefock
C
C Used from common blocks:
C  INFINP : FLAG(*)
C  INFORB : NCMOT,NNBAST,NSYM
C  INFOPT : EPOT,EMCSCF,EMY,EACTIV,GRDNRM
C  SCBRHF : THRRHF,MXDIIS, ?
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
#include "nuclei.h"
C
Cholesky
C
C------------------------
C
      PARAMETER (THDDEF = 1.0D-5)
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0, D1 = 1.0D0, DM1 = -1.0D0)
      PARAMETER (THREVC = 0.01D0, THRINC = 1.0D-10, CNVFAC = 0.8D0)
C
      LOGICAL   ONLYFC, C2DIIS, DIFDEN, DODIFDEN, NEWTHR, RHFWOP, OLDWOP
      LOGICAL   HSFOCK, LNEWOCC
      DIMENSION MISH(8)
      CHARACTER*8 CMOLBL
C
      CALL QENTER('DIIS_CTL')
      KFREE  = KFRSAV
      LFREE  = LFRSAV
      IFTHSV = IFTHRS
C
C     save current IOPRHF for check below
      IOPRH1 = IOPRHF
      IF (NSYM .EQ. 1) AUTOCC = .FALSE.
C
C     Get one-electron Hamiltonian
C
      CALL SIRH1(H1AO,WRK(KFREE),LFREE)
C
      IF (LUW4 .NE. LUPRI .OR. IPR_DIIS .LE. 1) THEN
      IF (FLAG(16) .OR. PCM .OR. QM3 .OR. QMMM) THEN
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
      else if (peqm) then
        write(luw4,'(/a)')
     &      ' Iter      Total energy       PE energy           '//
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
      ELSE IF (NASHT .EQ. 1 .OR. HSROHF) THEN
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
C     Level shift ("restricted step" type step reduction)
C     - initial value: zero
C     - used to modify FD matrix: FDAO(shifted) = FDAO + SHFTLVL*DCAO
C       - this corresponds to increasing occ-unocc orb.en. gap by
C         SHFTLVL and approximately to add -SHFTLVL to the
C         approximate Hessian diagonal implicit in the Roothaan diag.
C     -- Feb. 2001, Hans Joergen Aa. Jensen.
C
c     DEFLVL is read from SCF INPUT/.SHIFT
      SHFTLVL = DEFLVL
      SHFTFAC = 0.7D0 
C ---  hjaaj: find safe screening factor depending on
C             convergence threshold:
C     SCRFAC = 100*NBAST
      SCRFAC = 100*N2BAST
      IFTMAX = INT(-LOG10(THDIIS/SCRFAC)) + 1
      IFTMIN = INT(-LOG10(1.0D-1/SCRFAC)) + 1
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
C       keep vectors until GRDNRM fallen a factor 100
      DAMP   = D0
      EMCLOW = D0
      GNRMSV = D0
      ITDIIS = 0
      JTDIIS = 0
      ITNOCC = 0
      NEWOCC = 0
      LNEWOCC = AUTOCC
C     if (LNEWOCC) then write new BASINFO to SIRIUS.RST in NEWORB
      NOSCIL = 0
C------------------
      DODIFDEN = .NOT.FLAG(49) .AND. .NOT.HSROHF
     &     .AND. .NOT.DODFT    .AND. .NOT.DOHFSRDFT
     &     .AND. .NOT.FLAG(16) .AND. .NOT.PCM
     &     .AND. .NOT.QM3      .AND. .NOT. QMMM
     &     .AND. .NOT.AOSAVE   .AND. .NOT. CHOINT
     &     .and. .not. peqm    .and. .not. pe_savden
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
         IF (AOSAVE .OR. CHOINT) THEN
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
  100 CONTINUE
         ITDIIS = ITDIIS + 1
         JTDIIS = JTDIIS + 1
         ITNOCC = ITNOCC + 1
         INDEVC = MOD(JTDIIS-1,MXERRV) + 1
         NEVC   = MIN(MXERRV,JTDIIS)
         WRITE(LUW4,'(A)')
     &      ' -----------------------------------------------------'//
     &      '------------------------'
C
C
         IF (DIRFCK .AND. IFTHRS.LT.20) THEN
            IF (IPR_DIIS .GT. 1) THEN
               WRITE (LUPRI,'(/A,2I5/A,I5,L5)')
     &         ' *** START DIIS iteration ',ITDIIS,JTDIIS,
     &         ' Screening settings (-IFTHRS, DIFDEN)',
     &         -IFTHRS, DIFDEN
            ELSE IF (IPR_DIIS .GE. 0) THEN
               WRITE (LUPRI,'(I4,A,I5,L5)') 
     &         ITDIIS,'  Screening settings (-IFTHRS, DIFDEN)',
     &         -IFTHRS, DIFDEN
            END IF
         END IF
         CALL FLSHFO(LUPRI)
         HSFOCK=((DODFT .AND. NASHT .GE. 1) .OR. HSROHF)
C
C  ***** Construct folded inactive and active parts of the one-electron
C        density matrix in AO-basis (DCAO and DVAO). DV is the active
C        part of the one-electron density matrix in MO-basis.
C
         CALL FCKDEN((NISHT.GT.0.OR.HSFOCK),(.NOT.ONLYFC),
     *                DCAO,DCAO(1+N2BASX),CMO,DV,WRK(KFREE),LFREE)
C        CALL FCKDEN(GETDC,GETDV,DCAO,DVAO,CMO,DV,WRK,LWRK)
C
C     For a restricted open shell DFT calculation, or high spin
C     Let the matrices be Fc=Fi+Fa and Fc-Fo=Fa-Q=-Fa(exch)
C     generated by total and spin (=active) densities
C     This choice avoids the calculation of a third Q matrix
C
         IF (HSFOCK) THEN
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
                ENDIF
             ENDDO
             CALL WCLOSE2(LUCCMO,FCCMO,'KEEP')
         ENDIF
C
         CALL MEMGET('REAL',KTMP1,NNBASX,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KTMP2,N2BASX,WRK,KFREE,LFREE)
         IF (NISHT.GT.0 .OR. HSFOCK) THEN
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
         CALL MEMREL('DIIS_CTL.DIFDEN',WRK,KFRSAV,KTMP1,KFREE,LFREE)
C
C  ***** Calculate 2-electron part of Fock matrix with
C        differential density
C
Cef begin
C         write(lupri,*) 'DCAO: '
C         call output(DCAO,1,nbast,1,nbast,nbast,nbast,1,lupri)
Cef end
         CALL FCK2AO(ONLYFC,FCAO,DCAO,WRK,KFREE,LFREE)
         IF (HSFOCK) THEN
C
C Reset densities for call to dft
C
            CALL DSCAL(N2BASX,-D1,DCAO(1+N2BASX),1)
            CALL DAXPY(N2BASX,-D1,DCAO(1+N2BASX),1,DCAO,1)
         END IF
C
C        FAOSAV(1,INDEVC,1) = FD2AO = FC2AO + FV2AO
C        FAOSAV(1,INDEVC,2) = FV2AO
C
         IF (NSYM .GT. 1) THEN
            CALL MEMGET('REAL',KTMP1,NNBASX,WRK,KFREE,LFREE)
            CALL DGETSP(NBAST,FCAO,WRK(KTMP1))
            CALL PKSYM1(WRK(KTMP1),FAOSAV(1,INDEVC,1),NBAS,NSYM,1)
            IF (.NOT.ONLYFC) THEN
               CALL DGETSP(NBAST,FCAO(1+N2BASX),WRK(KTMP1))
               CALL PKSYM1(WRK(KTMP1),FAOSAV(1,INDEVC,2),NBAS,NSYM,1)
            END IF
            CALL MEMREL('DIIS_CTL.FAOSAV',WRK,KFRSAV,KTMP1,KFREE,LFREE)
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
               IF (.NOT.HSFOCK)
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
C           EACTIV = DV(u,u)FC(u,u) = FC(u,u)
C                  = DVAO(p,q)FC2AO(p,q) + DVAO(p,q)H1AO(p,q)
            IF (HSFOCK) THEN
               EACTIV = D0
               DO ISYM = 1, NSYM
                  EACT1 = +DDOT(NNBAS(ISYM),
     &                       DAOSAV(1+IIBAS(ISYM),2),1,
     &                       FAOSAV(1+IIBAS(ISYM),INDEVC,2),1) * DP5
                  EACTIV = EACTIV + EACT1 
               END DO
            ELSE
            EACTIV = DDOT(NNBAS(IOPRHF),DAOSAV(1+IIBAS(IOPRHF),2),1,
     &                           FAOSAV(1+IIBAS(IOPRHF),INDEVC,1),1)
     &             + DDOT(NNBAS(IOPRHF),DAOSAV(1+IIBAS(IOPRHF),2),1,
     &                           H1AO(1+IIBAS(IOPRHF)),1)
            CALL DAXPY(NNBAST,D1,FAOSAV(1,INDEVC,2),1,
     &                 FAOSAV(1,INDEVC,1),1)
            END IF
         END IF
         IF (DODFT) CALL DFT_ADD_KS(ONLYFC,EMY,DCAO,NDMAT,
     &                  FAOSAV(1,INDEVC,1),FAOSAV(1,INDEVC,2),
     &                  WRK(KFREE),LFREE,IPR_DIIS)
C
         ESOLT = D0

         IF (QM3 .AND. .NOT. DOCCSD) THEN
            CALL MEMGET('REAL',KFCSOL,NNBAST,WRK,KFREE,LFREE)
            CALL DZERO(WRK(KFCSOL),NNBAST)
            CALL QM3FCK(DAOSAV(1,1),DAOSAV(1,2),
     &                  WRK(KFCSOL),ESOLT,ENSQM3,EPOQM3,EELEL,
     &                  WRK(KFREE),LFREE,IPR_DIIS)
            CALL DAXPY(NNBAST,D1,WRK(KFCSOL),1,FAOSAV(1,INDEVC,1),1)
            CALL MEMREL('DIIS_CTL.QM3FCK',WRK,KFRSAV,KFCSOL,KFREE,LFREE)
         ELSE
            EELEL = 0.0D0
         END IF

         IF (QMMM) THEN
            IF (MMTIME) DTIME = SECOND()
            CALL MEMGET('REAL',KFCSOL,NNBAST,WRK,KFREE,LFREE)
            CALL DZERO(WRK(KFCSOL),NNBAST)
            CALL QMMMFCK(DAOSAV(1,1),DAOSAV(1,2),
     &                  WRK(KFCSOL),ESOLT,
     &                  WRK(KFREE),LFREE,IPQMMM)
            CALL DAXPY(NNBAST,D1,WRK(KFCSOL),1,FAOSAV(1,INDEVC,1),1)
            CALL MEMREL('DIIS_CTL.QMMMFCK',WRK,KFRSAV,
     &                  KFCSOL,KFREE,LFREE)
            IF (MMTIME) THEN
              DTIME = SECOND() - DTIME
              TMMFCK = TMMFCK + DTIME
            ENDIF
         ENDIF

         if (peqm) then
           allocate(pefock(nnbasx))
           if (nasht .gt. 0 .and. .not. (dodft .or. hsrohf)) then
             call pe_master(runtype='fock',
     &                      denmats=(daosav(:,1) + daosav(:,2)),
     &                      fckmats=pefock, nmats=1, Epe=pe_energy,
     &                      dalwrk=wrk(kfree:lfree))
           else
             call pe_master(runtype='fock', denmats=daosav(:,1),
     &                      fckmats=pefock, nmats=1, Epe=pe_energy,
     &                      dalwrk=wrk(kfree:lfree))
           end if
           esolt = pe_energy(1)
           call daxpy(nnbasx, 1.0d0, pefock, 1, faosav(1,indevc,1), 1)
           deallocate(pefock)
         end if

         

         IF (FLAG(16)) THEN
            CALL MEMGET('REAL',KFCSOL,NNBAST,WRK,KFREE,LFREE)
            CALL SOLFCK(DAOSAV(1,1),DAOSAV(1,2),
     &                  WRK(KFCSOL),DUMMY,TNLM,.FALSE.,ESOLT,
     &                  WRK(KFREE),LFREE,IPR_DIIS)
            CALL DAXPY(NNBAST,D1,WRK(KFCSOL),1,FAOSAV(1,INDEVC,1),1)
            CALL MEMREL('DIIS_CTL.SOLFCK',WRK,KFRSAV,KFCSOL,KFREE,LFREE)
         ELSE IF (PCM) THEN
            CALL MEMGET('REAL',KFCSOL,NNBASX,WRK,KFREE,LFREE)
C
            CALL PCMFCK(DAOSAV(1,1),DAOSAV(1,2),WRK(KFCSOL),
     &                  DUMMY,.FALSE.,ESOLT,WRK(KFREE),
     &                  LFREE,IPR_DIIS)
C
            CALL DAXPY(NNBAST,DM1,WRK(KFCSOL),1,FAOSAV(1,INDEVC,1),1)
            CALL MEMREL('DIIS_CTL.PCMFCK',WRK,KFRSAV,KFCSOL,KFREE,LFREE)
         END IF
C
         EMCDIF = EMCSCF
         EMCSCF = EMY + EACTIV + EPOT + ESOLT

         IF (QM3 .AND. (.NOT. OLDTG)) EMCSCF = EMCSCF - EELEL

         EMCDIF = EMCSCF - EMCDIF
         EGSAV(1,INDEVC) = EMCSCF
C
         CALL DIIS_EVC(SMOAO,DAOSAV(1,1),DAOSAV(1,2),H1AO,
     &               FAOSAV(1,INDEVC,1),FAOSAV(1,INDEVC,2),CMO1,
     &               EVCSAV(1,INDEVC),WRK,KFREE,LFREE,IPR_DIIS)
C        CALL DIIS_EVC(SMOAO,DCAO,DVAO,H1AO,FDAO,FVAO,CMO1,
C    &               ERRVEC,WRK,KFRSAV,LFRSAV,IPR_DIIS)
         IF (IPR_DIIS .GT. 20) THEN
            WRITE (LUPRI,*) ' Test output of EVCSAV:'
            CALL OUTPUT(EVCSAV,1,NNORBT,1,NEVC,NNORBT,NEVC,1,LUPRI)
         ELSE IF (IPR_DIIS .GT. 10) THEN
            WRITE (LUPRI,*) ' Error vector, DIIS iteration',ITDIIS
            CALL OUTPKB(EVCSAV(1,INDEVC),NORB,NSYM,1,LUPRI)
         END IF
         GRDNRM = DNRM2(NNORBT,EVCSAV(1,INDEVC),1)
         EGSAV(2,INDEVC) = GRDNRM
         IF (LUW4 .NE. LUPRI .OR. IPR_DIIS .GT. 1) THEN
            WRITE (LUPRI,'(/A/A,I5,1P,G22.12,D15.5,D12.2,I5,L5)')
     &' ITDIIS, energy, error norm, change in energy, scr.thr., difden:'
     & ,' DIIS iter.',ITDIIS,EMCSCF,GRDNRM,EMCDIF,-IFTHRS,DIFDEN
            IF (AUTOCC) THEN
               WRITE (LUPRI,'(/A,8I4)')
     &         ' AUTOCC: current SCF occupation',(NISH(I),I=1,NSYM)
               IF (IOPRHF .NE. 0) WRITE (LUPRI,'(A,I4)')
     &         ' AUTOCC: current open shell symmetry',IOPRHF
            END IF
            CALL FLSHFO(LUPRI)
         END IF
         IF (LUW4 .NE. LUPRI .OR. IPR_DIIS .LE. 1) THEN
            IF (FLAG(16) .OR. PCM .OR. QM3 .OR. QMMM) THEN
               IF (.NOT.AUTOCC) THEN
                  WRITE (LUW4,'(I4,1P,G20.12,G20.12,D15.5,D12.2)')
     &               ITDIIS,EMCSCF,ESOLT,GRDNRM,EMCDIF
               ELSE IF (IOPRHF .EQ. 0) THEN
                 WRITE (LUW4,'(I4,1P,G20.12,G20.12,D11.2,D11.2,2X,8I3)')
     &                 ITDIIS,EMCSCF,ESOLT,GRDNRM,EMCDIF,
     &                 (NISH(I),I=1,NSYM)
               ELSE
                 WRITE (LUW4,'(I4,1P,G20.12,G20.12,D11.2,D11.2,
     &                 I4,2X,8I3)')
     &                 ITDIIS,EMCSCF,ESOLT,GRDNRM,EMCDIF,IOPRHF,
     &                 (NISH(I),I=1,NSYM)
               END IF
            else if (peqm) then
               write(luw4,'(i4,1p,g20.12,g20.12,d15.5,d12.2)')
     &          itdiis, emcscf, esolt, grdnrm, emcdif
            ELSE
               IF (.NOT.AUTOCC) THEN
                  WRITE (LUW4,'(I4,1P,G22.12,D15.5,D12.2,I5)')
     &                 ITDIIS,EMCSCF,GRDNRM,EMCDIF,NEVC
               ELSE IF (IOPRHF .EQ. 0) THEN
                  WRITE (LUW4,'(I4,1P,G20.12,D11.2,D11.2,2X,8I3)')
     &                 ITDIIS,EMCSCF,GRDNRM,EMCDIF,(NISH(I),I=1,NSYM)
               ELSE
                  WRITE (LUW4,'(I4,1P,G20.12,D11.2,D11.2,I4,2X,8I3)')
     &                 ITDIIS,EMCSCF,GRDNRM,EMCDIF,IOPRHF,
     &                 (NISH(I),I=1,NSYM)
               END IF
            END IF
            CALL FLSHFO(LUW4)
         END IF
         IF (MOLDEN) CALL MOSCFCON(ITDIIS,EMCSCF,.FALSE.)
C
C        Test and exit if energy is increasing
C        (from itdiis .ge. 7, second iteration is Roothaan step
C        and often increases in energy; fourth iteration has been
C        observed to increase in energy (perhaps in connection with
C        new occupation for AUTOCC)).
C        961113-hjaaj: use new ITNOCC>2 for test; ITNOCC counts
C        number of diis iterations since last change of occupation.
C        970326-hjaaj: changed to ITNOCC>4 test; have today observed
C        increase in energy in iter 4 with fixed HF OCCUP; let us
C        allow two oscillations before exiting.
C        970402-hjaaj: Always allow two oscillations (NOSCIL)
C        after a "NEWOCC"
C
         EMCLOW = MIN(EMCLOW,EMCSCF)
         EMCDIFERR = THRINC*ABS(EMCLOW)
C        ... if EMCDIF positive, but increase < THRINC*abs(EMCLOW) then
C            the increase are in digits which are not reliable
C        ... by comparing to EMCDIFERR, an estimate of the numerical round-off
C            error, this modified test should only be positive if the exact
C            energy has increased.
C
         IF (EMCDIF .GT. EMCDIFERR .AND. ITNOCC .GT. 4) THEN
            NOSCIL = NOSCIL + 1
         IF (NOSCIL .GT. 10) THEN
            IF (AUTOCC .AND. NSYM .GT. 1 .AND. ITNOCC .LE. 6) THEN
C              ... ITNOCC test because do not abort if same
C                  occupation has been used in 6 iterations/970402-hjaaj
               WRITE (LUPRI,4010)
               IF (LUW4 .NE. LUPRI) WRITE (LUW4,4010)
               CALL QUIT(
     &         'FATAL ERROR: DIIS diverging and automatic occupation')
            END IF
            ICONV = 0
            WRITE (LUPRI,4020)
            IF (LUW4 .NE. LUPRI) WRITE (LUW4,4020)
            CALL READMO(CMO,9,DUMMY,1)
C           CALL READMO(CMO,JRDMO,WRK,LFREE)
            GO TO 1100
         END IF
         END IF
C
 4010 FORMAT(' DIIS aborted because energy increased in this iteration.'
     &      /' Program aborted because this means automatic symmetry'
     &      /' designation of occupied orbitals is not possible!',
     &      /' Please send this output to hjj@dou.dk !!!',
     &      /' You may try again with explicit specfication of .DOUBLY')
 4020 FORMAT(' DIIS aborted because energy increased in this iteration!'
     &      /' Backstep to orbitals from previous iteration.')
C
         IF (GRDNRM .LE. THDIIS) THEN
            IF (MOLDEN) CALL MOSCFCON(ITDIIS,EMCSCF,.TRUE.)
            ICONV = 1
            WRITE (LUPRI,4110) ITDIIS,EMCSCF,GRDNRM
            IF (LUW4 .NE. LUPRI) WRITE (LUW4,4110) ITDIIS,EMCSCF,GRDNRM
C
Cholesky
            IF (IPR_DIIS .GE. -1) THEN
              IF (CHOINT .AND. (IPR_DIIS .GE. 3)) THEN
                 WRITE(LUPRI,4116) 'Coulomb part ',CSIRF
                 WRITE(LUPRI,4116) 'Exchange part',XSIRF
                 WRITE(LUPRI,4116) 'Cholesky I/O ',RSIRF
                 WRITE(LUPRI,4116) 'MO transform.',OSIRF
                 WRITE(LUPRI,4116) 'MO sorting   ',SSIRF
              ENDIF
              WRITE(LUPRI,4115) TSIRF
              IF (QMMM) CALL QMMMTIMES('SIRIUS')
              IF (LUW4 .NE. LUPRI) THEN
                 IF (CHOINT .AND. (IPR_DIIS .GE. 3)) THEN
                    WRITE(LUPRI,4116) 'Coulomb part ',CSIRF
                    WRITE(LUPRI,4116) 'Exchange part',XSIRF
                    WRITE(LUPRI,4116) 'Cholesky I/O ',RSIRF
                    WRITE(LUPRI,4116) 'MO transform.',OSIRF
                    WRITE(LUPRI,4116) 'MO sorting   ',SSIRF
                 ENDIF
                 WRITE(LUW4,4115) TSIRF
              ENDIF
            ENDIF
Cholesky
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
 4110 FORMAT(/' *** DIIS converged in',I4,' iterations !',
     &       /'@    Converged SCF energy, gradient:',F20.12,1P,D12.2)
 4115 FORMAT( '   - total time used in SIRFCK :' ,F18.2,' seconds')
 4116 FORMAT( '   - total time used for ',A13,':',F12.2,' seconds')
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
C
C        Test if convergence is stopping because of numerical problems
C        (EMCDIF test should make sure we are in local region)
C
ckr         IF (ITNOCC .GT. 3 .AND. ABS(EMCDIF) .LT. THDDEF .AND.
ckr     &       GRDNRM .GT. CNVFAC*GNRMSV .AND. .NOT. NEWTHR) THEN
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
         IF (INDEVC.EQ.NEVC .AND. ITDIIS .GT. 3 .AND.
     &       GRDNRM .LE. 10.0D0) THEN
C        ... no wrap around yet - INDEVC .lt. NEVC requires more programming
C        ITDIIS/GRDNRM check: give initial oscillations a chance to settle down
C          before checking ... experience shows this will often be advantageous
C          (see below) /hjaaj
            IEND = NEVC - 1
            J = 0
            DO I = 1,IEND
C           IF (I .NE. INDEVC) THEN
               IF (EGSAV(1,I) .GT. EMCSCF + EMCDIFERR .AND.
     &             EGSAV(2,I) .LT. GRDNRM) THEN
C                  ... Oops! a previous iter with higher E and lower grd
                  J = J + 1
                  EGSAV(1,I) = D0
               END IF
C           END IF
            END DO
C
C           remove entries with higher E and lower grd, if any
C
            IF (J .GT. 0) THEN
               WRITE(LUPRI,'(/A/A,I3,A)')
     &   'Info: SCF gradient has been lower than now,',
     &   '      therefore',J,' old iterations are removed from DIIS.'
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
     &               CVEC,XMLBDA,EVCSAV,WRK,KFREE,LFREE,IPR_DIIS)
C        CALL DIIS_RED(C2DIIS,BMAT,INDEVC,NEVC,DAMP,THREVC,
C    *               CVEC,XLMBDA,ERRVEC,WRK,KFRSAV,LFRSAV,IPR_DIIS)
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
C        Do not try that with open-shell - it does not work./pawsa
#if 0 /* I am not sure what the test for stalling detection should be */
         IF (JTDIIS.GT.1 .AND.
     &       EMCDIF.GT.5D-5) THEN
#else
!  Mar. 2011 hjaaj: if abs(EMCDIF) .lt. 1.0D-5 the error can be grid related in DFT 
          IF (ONLYFC .AND. JTDIIS.GT.1 .AND. ITDIIS.GT.3 .AND.
     &        GRDNRM .LE. 10.0D0 .AND. ABS(EMCDIF) .GT. 1.0D-5 .AND.
     &       ( ABS(CVEC(INDEVC)) .LT. 0.1D0 .OR.
     &         ABS(CVEC(INDEVC)) .GT. 5.0D0 .OR.
     &         EMCDIF*(GRDNRM-GNRMSV) .LT. -1.0D-10) .AND. BCKSTP)
     &       THEN
C              .. test EMCDIF*grddif against -1.0D-10 instead of 0.0D0
C                 to avoid round-off problems here /hjaaj Sep 2005
#endif
c        IF (ITDIIS.GT.1 .AND. 
c    &       ( ABS(CVEC(INDEVC)) .LT. 0.2D0 .OR.
c    &         ABS(CVEC(INDEVC)) .GT. 2.0D0 ) .AND.
c    &       EMCDIF .LT. D0 .AND. GRDNRM .GT. GNRMSV) THEN
             WRITE(LUPRI,'(1P/A,2(/A,G11.3),/A,I3,A)')
     &       '!!! Info: DIIS restarted because it was stalled ...',
     &       ' - energy changed by  ',EMCDIF,
     &       ' - gradient changed by',GRDNRM - GNRMSV,
     &       ' - or strange C vector coeff. for current index (=',
     &       INDEVC,') :'
             WRITE (LUPRI,'(6F12.6)') (CVEC(I),I=1,NEVC)
             IF (EMCDIF .LT. D0) THEN
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
C            also restrict orbital rotation with a level shift
c            the formula is empirical EMCDIF is positive here.
             SHFTLVL = MAX(-12D0,SHFTLVL - MIN(4.0D0,ABS(EMCDIF)*10.D0))
c            and slow down  the level shifting annealing...
             SHFTFAC = MIN(SHFTFAC + (1.D0-SHFTFAC)*0.25,0.85D0)
         ELSE
             SHFTLVL =  SHFTLVL*SHFTFAC
c            make sure we get right HOMO-LUMO gap in the end.
             IF(ITDIIS.GT.1.AND.ABS(EMCDIF).LE.5.D-5) SHFTLVL = 0.D0
         END IF
C hjaaj 22-Jun-20005-B
C        Level shift is only working for closed shell calculations ...
C          (open shell problems reported by users ...)
         IF (NASHT .GT. 0) SHFTLVL = 0.0D0
C hjaaj 22-Jun-20005-E
         GNRMSV = GRDNRM
C
C        Iteration accepted; save current CMO for restart
         WRITE (CMOLBL,'(A4,I4)') 'DIIS',ITDIIS
         CALL NEWORB(CMOLBL,CMO,LNEWOCC)
C        ... REWIT1 false: do not destroy any GEOWALK information
         LNEWOCC = .FALSE.

C        print atomic populations in each iteration ? /hjaaj
         if (lim_poppri .gt. 0)
     &      call sirpop('DIIS ',DV,wrk(kfree),lfree)

C
         IF (HSFOCK) THEN
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
         IF (.NOT. ONLYFC .AND..NOT.HSFOCK) THEN
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
         IF (HSFOCK .OR. ONLYFC) THEN
            IF (IPR_DIIS .GE. 15) THEN
               IF (HSFOCK) THEN
                  WRITE (LUPRI,'(/A)')
     &            ' Fock matrix after FV correction'
               ELSE
                  WRITE (LUPRI,'(/A)')
     &            ' Fock matrix before diagonalization'
               END IF
               CALL OUTPKB(DCAO,NORB,NSYM,1,LUPRI)
            END IF
         ELSE   
           DO 300 ISYM = 1,NSYM
            IF (NASH(ISYM) .EQ. 0) GO TO 300
            IF (IPR_DIIS .GT. 15) THEN
               WRITE (LUPRI,'(/A,I3)')
     &            ' Fock matrix before FV correction, symmetry',ISYM
               CALL OUTPAK(DCAO(1+IIORB(ISYM)),NORB(ISYM),1,LUPRI)
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
               CALL OUTPAK(DCAO(1+IIORB(ISYM)),NORB(ISYM),1,LUPRI)
            END IF
  300      CONTINUE
         END IF
C
C        If requested (by SHFTLVL .ne. 0) :
C        Perform "restricted step" step reduction (i.e. level shift)
C        by adding 0.5*SHFTLVL*DENS(i,i) to all diagonal elements;
C        /hjaaj-Feb.2001

         IF (SHFTLVL .NE. D0) WRITE(LUPRI,'(I4,A,1P,D10.2)') ITDIIS,
     &   '  Level shift: occupied orbital energies shifted by',SHFTLVL
chjaaj: improve comment with something about "restricted step"
chjaaj: and "0.5*occ(i)* " /23may09 hjaaj
C
C        Diagonalize (level-shifted) FC (saved in DCAO)
C
         CALL ICOPY(8,NISH,1,MISH,1)
         JOPRHF = IOPRHF
         CALL FCKEIG(CMO,DCAO,SHFTLVL,WRK(KFREE),LFREE,IPR_DIIS)
C        test if AUTOCC has changed occupation
         IF (AUTOCC) THEN
          NTEST = 0
          DO ISYM = 1,8
            NTEST = NTEST + ABS(NISH(ISYM)-MISH(ISYM))
          END DO
          NTEST = NTEST + ABS(IOPRHF-JOPRHF)
          IF (NTEST .GT. 0) THEN
            NEWOCC = NEWOCC + 1
            LNEWOCC = .TRUE.
            ITNOCC = 0
            NOSCIL = 0
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
C
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
      IF (IOPRHF .NE. IOPRH1) RHFWOP = .FALSE.
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
      IF (ICONV .NE. 0. .AND. (IPR_DIIS .GT. 0 .OR.
     &    .NOT.(DOMP2 .OR. DOCI .OR. DOMC)) ) THEN
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

       ! needed to get specific PE energy contributions for printout
       ! when converged
       if (peqm) then
         if (hsrohf .or. (dodft .and. nasht .ge. 1)) then
           call pe_master(runtype='energy',
     &                    denmats=(daosav(:,1) + daosav(:,2)),
     &                    nmats=1, dalwrk=wrk(kfree:lfree))
         else
           call pe_master(runtype='energy',
     &                    denmats=daosav(:,1), nmats=1,
     &                    dalwrk=wrk(kfree:lfree))
         end if
       end if
       if (pe_mep) then
        call pe_master(runtype='mep', denmats=daosav(:,1), nmats=1,
     &                 dalwrk=wrk(kfree:lfree))
       end if
       ! save density matrix 
       if (pe_savden) then
        if (nsym > 1) then
            call quit('PEQM cannot be combined with symmetry.')
        end if
        call pe_save_density(daosav(:,1), fcao(1:norbt*(norbt+1)/2),
     &                       cmo(1:nbast*norbt), nbast, nisht, norbt,
     &                       cord(:,1:natoms), charge(1:natoms),
     &                       wrk(kfree:lfree))
       end if
C
       IF (IPRENR .GT. 0) THEN
         WRITE (LUW4,'(//A)') ' *** SCF orbital energy analysis ***'
         IF (FLAG(16) .OR. PCM .OR. QM3 .OR. QMMM)
     &     WRITE (LUW4,'(A)') '    (incl. solvent contribution)'
         if (peqm) write(luw4,'(a)')
     &    '   (incl. contribution from Polarizable Embedding potential)'
         CALL RHFENR(IPRI4,LUW4,FCAO, WRK(KFREE),LFREE)
C-------------------
         IF (LUPRI .NE. LUW4) THEN
            WRITE (LUPRI,'(//A)') ' *** SCF orbital energy analysis ***'
            IF (FLAG(16) .OR. PCM .OR. QM3 .OR. QMMM)
     &        WRITE (LUPRI,'(A)') '    (incl. solvent contribution)'
            if (peqm) write(lupri,'(A)')
     &    '   (incl. contribution from Polarizable Embedding potential)'
            CALL RHFENR(IPRI6,LUPRI,FCAO, WRK(KFREE),LFREE)
         END IF
C-------------------
C        CALL RHFENR(IPRINT,LUPRI,FC,FV,SCRA,LSCRA)
      END IF
C
C     If (FLAG(25) .OR. INERSI) write to SIRIFC
C       FLAG(25) is true for .ABACUS (geometry optimization)
C                         or .RESPONS (response calculation)
C                         or .INTERFACE (request for write of SIRIFC)
C       INERSI is true if this is initial state calculation
C              in a solvent calculation with inertial polarization
C
      IF ( IWRIFC .EQ. 1 ) THEN
         IF (IPRI6 .GE. 0) WRITE (LUPRI,'(/A)')
     *      ' >>> Writing SIRIFC interface file <<<'
         IF (ICI0 .EQ. 6) THEN
Chj      ... if (GEOWLK) then
            IF (.NOT. OPTNEW)
     *         WRITE (LUW4,'(//A/)')
     *            ' Check ratio for geometry walk with converged SCF :'
            JWSTEP = 1
            EMCOLD = EMCGEO
            DEPRED = DEPGEO
            CALL SIRSTP(JWSTEP,DUMMY,DUMMY,DUMMY,DUMMY,1)
Chj         ... WLKREJ in gnrinf.h will be true if geom. step rejected.
         END IF

         LPV = NNASHX*NNASHX
         CALL MEMGET('REAL',KPV,LPV,WRK,KFREE,LFREE)
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

         CALL MEMGET('REAL',KGORB,NWOPH,WRK,KFREE,LFREE)
         CALL DZERO(WRK(KGORB),NWOPH)
Chj      ... set orbital gradient to zero (i.e. completely converged!)

         LFQ = NASHT*NORBT
         CALL MEMGET('REAL',KFQ,LFQ,WRK,KFREE,LFREE)
         IF (NASHT.GT.0) THEN
C           put FV in DCAO
            CALL UTHUB(FAOSAV(1,INDEVC,2),DCAO,CMO,WRK(KFREE),
     &                 NSYM,NBAS,NORB)
C           form FC = FD - FV in FCAO
            CALL DAXPY(NNORBT,DM1,DCAO,1,FCAO,1)
            IF (NASHT .EQ. 1) THEN
               CALL DZERO(WRK(KFQ),NASHT*NORBT)
C              FQ matrix is zero matrix for single open shell RHF
            ELSE
               WRK(KFQ:KFQ+LFQ-1) = -9.87654321D123
               ! see notes in connection with PV above
            END IF
         END IF
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
         CALL MEMGET('REAL',KDUM,0,WRK,KFREE,LFREE)
         CALL WR_SIRIFC(CMO,DV,WRK(KPV),FCAO,DCAO,
     &               WRK(KFQ),WRK(KDUM),WRK(KDUM),WRK(KDUM),
     &               WRK(KFREE),LFREE,WRK(KGORB),TNLM,.FALSE.,WRK(KDUM))
C        CALL WR_SIRIFC(CMO,DV,PV,FC,FV,FQ,CREF,FCAC,H2AC,WRK,LFREE,
C    &               GORB,ERLM,ORBHES,XINDX)
C
         CALL MEMREL('DIIS_CTL.WR_SIRIFC',WRK,KFRSAV,KPV,KFREE,LFREE)
      END IF
C
C        Transform integrals for abacus or response, 
C        if converged and if not followed by higher level,
C        and if not Cholesky decomposition
C
      IF (ICONV .NE. 0 .AND. ABS(ITRFIN) .LE. 10 .AND.
     &    .NOT.(DOMP2 .OR. DOCI .OR. DOMC) .AND.
     &    .NOT. CHOINT) THEN
         IF (IPRI6 .GE. 0) WRITE (LUPRI,'(/A,I3,A)')
     &      ' >>> Transforming 2-el. integrals acc. to'//
     &      ' .FINAL TRANSFORMATION =',ITRFIN,' <<<'
         JTRLVL = ITRFIN
         CALL TRACTL(JTRLVL,CMO,WRK(KFREE),LFREE)
      END IF
C
C ***********************************************************
C
 1100 CONTINUE
C     hjaaj-aug99: restore original IFTHRS value
      IFTHRS = IFTHSV
      CALL QEXIT('DIIS_CTL')
      RETURN
      END
C  /* Deck DIIS_EVC */
      SUBROUTINE DIIS_EVC(SMOAO,DCAO,DVAO,H1AO,FDAO,FVAO,CMO1,
     &                  ERRVEC,WRK,KFRSAV,LFRSAV,IPR_DIIS)
C
C
#include "implicit.h"
      DIMENSION DCAO(NNBAST), DVAO(NNBAST), FDAO(NNBAST), FVAO(NNBAST)
      DIMENSION H1AO(NNBAST)
      DIMENSION SMOAO(N2BAST), CMO1(NCMOT), ERRVEC(NNORBT), WRK(*)
C
C Used from common blocks:
C  INFINP : LNOROT,NOROT()
C  INFORB : N2BAST, NNBAST, NCMOT, NSYM, ..., NFROT, ...
C  INFDIM : N2BASM
C
#include "maxorb.h"
#include "priunit.h"
#include "infinp.h"
#include "inforb.h"
#include "infdim.h"
#include "infpri.h"
C
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DM1 = -1.0D0)
C
      CALL QENTER('DIIS_EVC')
      KFREE = KFRSAV
      LFREE = LFRSAV
      CALL MEMGET('REAL',KTMP1,N2BASM,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KTMP2,N2BASM,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KTMP3,N2BASM,WRK,KFREE,LFREE)
      IF (HSROHF.OR.(DODFT.AND.NASHT.GE.1)) THEN
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
     &                  NBASI,NORBI,1,LUPRI)
            WRITE (LUPRI,*) 'DIIS_EVC SMOAO = 4 Ct SAO:'
            CALL OUTPUT(SMOAO(1+I2BAS(ISYM)),1,NORBI,1,NBASI,
     &                  NORBI,NBASI,1,LUPRI)
         IF (IPR_DIIS .GT. 50) THEN
            CALL DGEMM('N','N',NORBI,NORBI,NBASI,1.D0,
     &                 SMOAO(1+I2BAS(ISYM)),NORBI,
     &                 CMO1(1+ICMO(ISYM)),NBASI,0.D0,
     &                 WRK(KTMP3),NORBI)
            WRITE (LUPRI,*) 'DIIS_EVC SMO = 4 Ct SAO C:'
            CALL OUTPUT(WRK(KTMP3),1,NORBI,1,NORBI,
     &                  NORBI,NORBI,1,LUPRI)
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
     &                     NBASI,NBASI,1,LUPRI)
               WRITE (LUPRI,*) 'DIIS_EVC FDAO:'
               CALL OUTPUT(WRK(KTMP3),1,NBASI,1,NBASI,
     &                     NBASI,NBASI,1,LUPRI)
               WRITE (LUPRI,*) 'DIIS_EVC TMP1 = DCAO FDAO:'
               CALL OUTPUT(WRK(KTMP1),1,NBASI,1,NBASI,
     &                     NBASI,NBASI,1,LUPRI)
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
     &                     NBASI,NBASI,1,LUPRI)
               WRITE (LUPRI,*) 'DIIS_EVC FCAO:'
               CALL OUTPUT(WRK(KTMP3),1,NBASI,1,NBASI,
     &                     NBASI,NBASI,1,LUPRI)
               WRITE (LUPRI,*) 'DIIS_EVC TMP1 = DCAO FDAO + DVAO FCAO:'
               CALL OUTPUT(WRK(KTMP1),1,NBASI,1,NBASI,
     &                     NBASI,NBASI,1,LUPRI)
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
         CALL DGETAP(NORBI,WRK(KTMP1),ERRVEC(1+IIORB(ISYM)))
         IF (IPR_DIIS .GE. 15) THEN
            WRITE (LUPRI,'(/A/A/A,I3)')
     &      ' (DIIS_EVC) error vector = orbital gradient vector =',
     &      '          4 (DCmo FDmo + DVmo FCmo) = ',
     &      '          4 Ct SAO (DCAO FDAO + DVAO FCAO) C for symmetry',
     &      ISYM
            CALL OUTPUT(WRK(KTMP1),1,NORBI,1,NORBI,NORBI,NORBI,1,LUPRI)
         END IF
  800 CONTINUE
      IF (HSROHF.OR.(DODFT.AND.NASHT.GE.1)) THEN
         CALL DAXPY(NNBAST,D1,DVAO,1,DCAO,1)
         CALL DSCAL(NNBAST,DM1,DVAO,1)
      END IF
      CALL MEMREL('DIIS_EVC',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('DIIS_EVC')
      RETURN
      END
C  /* Deck DIIS_RED */
      SUBROUTINE DIIS_RED (C2DIIS,BMAT,INDEVC,NEVC,DAMP,THREVC,
     *                   CVEC,XLMBDA,ERRVEC,WRK,KFRSAV,LFRSAV,IPR_DIIS)
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
      DIMENSION BMAT(*),CVEC(*), ERRVEC(NNORBT,*), WRK(*)
      LOGICAL C2DIIS
C
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
C
C Used from common blocks:
C   INFORB : NNORBT
C
#include "priunit.h"
#include "inforb.h"
#include "infpri.h"
C
      IROW(I) = I*(I-1)/2
C
      CALL QENTER('DIIS_RED')
      KFREE = KFRSAV
      LFREE = LFRSAV
C
      IF (IPR_DIIS .GE. 5) WRITE (LUPRI,7555) C2DIIS,NEVC
 7555 FORMAT(/' (DIIS_RED) CONTROL PARAMETERS C2DIIS =',L10,
     *        ' NEVC =',I5)
C
C
C Section 1: extend BMAT with a new row corresponding to
C            the new error vector
C
      K = INDEVC
      KROW = IROW(K)
      DO 40 L = 1,K
         BMAT(KROW+L) = DDOT(NNORBT,ERRVEC(1,L),1,ERRVEC(1,K),1)
   40 CONTINUE
      DO 50 L = K+1,NEVC
         BMAT(IROW(L)+K) = DDOT(NNORBT,ERRVEC(1,L),1,ERRVEC(1,K),1)
   50 CONTINUE
C
      LBMAT  = IROW(NEVC+1)
      IF (C2DIIS) THEN
         LBAUG  = LBMAT
      ELSE
         LBAUG  = IROW(NEVC+2)
      END IF
      CALL MEMGET('REAL',KBTMP,LBAUG,WRK,KFREE,LFREE)
      CALL DCOPY(LBMAT,BMAT,1,WRK(KBTMP),1)
      IF ( IPR_DIIS.GE.5 ) THEN
         WRITE (LUPRI,'(//A)') ' (DIIS_RED) B matrix:'
         CALL OUTPAK(WRK(KBTMP),NEVC,1,LUPRI)
      END IF
      IF (.NOT. C2DIIS) GO TO 2000
C
C *******************************************************
C
C SOLVE LINEAR SET OF EQUATIONS AS AN EIGENVALUE EQUATION
C
C
      CALL MEMGET('REAL',KEVEC,NEVC*NEVC,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KWJ  ,NEVC     ,WRK,KFREE,LFREE)
      CALL MEMGET('INTE',KIWJ ,NEVC     ,WRK,KFREE,LFREE)
C
      CALL DUNIT(WRK(KEVEC),NEVC)
      CALL JACO(WRK(KBTMP),WRK(KEVEC),NEVC,NEVC,NEVC,WRK(KWJ),WRK(KIWJ))
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
     *      ' (DIIS_RED) C2-DIIS B matrix eigenvalues. IOK =',IOK
         CALL OUTPUT(WRK(KBTMP),1,1,1,NEVC,1,NEVC,1,LUPRI)
      END IF
      IF ( IPR_DIIS.GE.4) THEN
         WRITE(LUPRI,'(/A)') ' - and B matrix eigenvectors:'
         CALL OUTPUT(WRK(KEVEC),1,NEVC,1,NEVC,NEVC,NEVC,1,LUPRI)
      END IF
C
      GO TO 9999
C
C ********************************************
C
C SOLVE LEVEL SHIFTED linear EQUATIONS
C
 2000 CONTINUE
C
      NLEQ  = NEVC + 1
      CALL MEMGET('REAL',KRHS ,NLEQ,WRK,KFREE,LFREE)
      CALL MEMGET('INTE',KIPVT,NLEQ,WRK,KFREE,LFREE)
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
     *      ' (DIIS_RED) B matrix damping factor =',DAMP
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
     *      ' (DIIS_RED) Solutions to C1-DIIS set of linear equations'
         CALL OUTPUT(CVEC,1,NLEQ,1,1,NLEQ,1,1,LUPRI)
      END IF
      IF (INFO.NE.0) THEN
         WRITE(LUPRI,8500) INFO
         CALL QTRACE(LUPRI)
         CALL QUIT('DIIS_RED: no solution to linear equations')
      END IF
 8500 FORMAT(/' (DIIS_RED) Solution not obtained to linear equations'
     *       /T11,'Check if matrix is singular, LINPACK DSP code =',I3)
C
C *** End of subroutine DIIS_RED
C
 9999 CONTINUE
      IF (IPR_DIIS .GT. 1) THEN
         WRITE (LUPRI,'(/A,L2,A,1P,D10.2)')
     &      ' DIIS C vector; C2DIIS =',C2DIIS,' LAMBDA =',XLMBDA
         WRITE (LUPRI,'(6F12.6)') (CVEC(I),I=1,NEVC)
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
         ENDIF         
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
c     CALL OUTPKB(WRK(KKSPA),NORB,NSYM,1,LUPRI)
c     WRITE(LUPRI,'(/A)') 'ksmB: '
c     CALL OUTPKB(WRK(KKSPB),NORB,NSYM,1,LUPRI)
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
      SUBROUTINE FCKEIG(CMO,FC,SHFTLVL,SCRA,LSCRA,IPR_DIIS)
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
C
C Output:
C  CMO; molecular orbitals diagonalizing Fock matrix
C  FC;  the orbital energies
C
C Scratch:
C  SCRA; general scratch area
C
#include "implicit.h"
      DIMENSION CMO(*),FC(*),SCRA(*)
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
#ifdef MOD_DFTASC
#include "dftacb.h"
#endif
C
      LOGICAL LSAVE4,LSAVE6
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0)
C
      CALL QENTER('FCKEIG')
C
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
      DO 200 ISYM = 1,NSYM
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
            IF (NASH(ISYM) .GT. 0) THEN
               K = NISH(ISYM) + 1
               FC(IFSYM+IROW(K+1)) = FC(IFSYM+IROW(K+1)) + DP5*SHFTLVL
            END IF
         END IF
C        zero rows and columns corresponding to frozen orbitals
         JFRZI = 0
         IF (NFRZI .GT. 0 .OR. LNOROT) THEN
            JNFRZ = 0
            DO 130 K = 1,NORBI
            IF (K .LE. NFRZI .OR. NOROT(K) .NE. 0) THEN
               IF (JNFRZ .GT. 0) GO TO 9000
C              ... exit if free orbital below this frozen orbital
C                  because then ORDRSS cannot be called
               JFRZI = JFRZI + 1
               KROW  = IROW(K)
               DO 110 L = 1,K-1
                  FC(IFSYM+KROW+L) = D0
  110          CONTINUE
               DO 120 L = K+1,NORBI
                  KL = (IFSYM+K)+IROW(L)
                  FC(KL) = D0
  120          CONTINUE
            ELSE
               JNFRZ = JNFRZ + 1
            END IF
  130       CONTINUE
         END IF
C
         CALL JACO(FC(IFSYM+1),CMO(ICSYM),NORBI,NORBI,NBASI,
     *             SCRA(KLAST),SCRA(KLAST+NORBI))
C        CALL JACO (F,VEC,NB,NMAX,NROWV,BIG,JBIG)
C
         DO 175 I=1,NORBI
            FC(IORBI+I)=FC(IFSYM+IROW(I+1))
  175    CONTINUE
C
         IF (AUTOCC) THEN
            CALL DCOPY(NORBI,FC(IORBI + 1),1,SCRA(KEIG + IORBI),1)
            DO 99 IK = 1, NORBI
               SCRA(KSYMS + IORBI + IK - 1) = ISYM
 99         CONTINUE
         END IF
C
C
         JCSYM = ICSYM + JFRZI*NBASI
         JORBI = IORBI + JFRZI + 1
         NNOTFR= NORBI - JFRZI
         CALL ORDRSS(CMO(JCSYM),FC(JORBI),ISSMO(JORBI),NNOTFR,NBASI)
C
         IF (IPR_DIIS .GE. 4) WRITE (LUPRI,'(/A,I2/,(5F15.5))')
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
#ifdef MOD_DFTASC
      IF (DFTASC) CALL DFTASC_GETHOMO(FC)
#endif
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
      CALL QEXIT('FCKEIG')
      RETURN
C
 9000 CONTINUE
      WRITE (LUW4,'(/A/A)')
     &  ' ERROR in FCKEIG, FCKEIG cannot handle ".FREEZE" for orbitals',
     &  ' unless they also could have been frozen with ".FROZEN"'
      CALL QTRACE(LUW4)
      CALL QUIT('FCKEIG error, cannot handle orbitals '//
     &          'frozen with ".FREEZE"')
      END
#ifdef MOD_DFTASC
 -- hjaaj: DFTASC is not active (June 09)
C  /* Deck gethom */
      SUBROUTINE DFTASC_GETHOMO(FD)
C find homo energy. tuh 16.05.03
#include "implicit.h"
      DIMENSION FD(*)
      PARAMETER (DBIG = 1.D+12)
#include "maxash.h"
#include "maxorb.h"
#include "mxcent.h"
#include "inforb.h"
#include "infind.h"
#include "dftacb.h"
C
      EHOMO  = -DBIG
      DO ISYM = 1,NSYM
         II = NOCC(ISYM)
         IF (II.GT.0) THEN
            EHOMOI = FD(IIORB(ISYM) + IROW(II) + II)
            EHOMO = MAX(EHOMO,EHOMOI)
         END IF
      END DO
      RETURN
      END
#endif
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

      CALL MEMGET('REAL',KS,NNBAST,WRK,KFREE,LFREE)
      N2BASI=IMXVEC(N2BAS,NSYM)
      CALL MEMGET('REAL',KDI,N2BASI,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KSI,N2BASI,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KFI,N2BASI,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KTMP1,N2BASI,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KTMP2,N2BASI,WRK,KFREE,LFREE)
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
