C
C...   Copyright (c) 2001 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 1.2
C...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
C...   O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
C...   T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
C...   H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
C...   H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
C...   T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
C...   B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
C...   and O. Vahtras"
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
C$Id: sirdiis.F,v 1.3 2001-10-01 13:22:20 vebjornb Exp $
C
#include <single.h>
C
#ifdef REV_LOG
se efter HJMAERKE
HJMAERKE 940816: if .FREEZE implemented in FCKEIG (requires new ORDRSS)
  then check all setting of MXDIIS = 0 and MAXFCK = 0 (I think MAXFCK
  also will work if .FREEZE implemented in FCKEIG).

$Log: not supported by cvs2svn $
Revision 1.2  2001/05/31 06:38:34  vebjornb
Updated copyright, authors and reference in all Fortran-files - VB

Revision 1.1.1.1  2001/02/08 13:33:32  hjj
*** Dalton Master Version *** cvs started 8-Feb-2001.

Revision 1.2  2001/01/17 10:41:33  vebjornb
Calls to *MPA*B* in arhpack.F have been replaced with DGEMM calls

Revision 1.1.1.1  2001/01/15 13:53:35  vebjornb
Dalton 1.2

Revision 1.6  2000/05/24 12:28:27  hjj
1) fixed solvent calculations (DIFDEN must be false for solvent)
2) fixed WRSIFC call for open shell and solvent
   (required change in SOLFCK also)

Revision 1.5  2000/05/16 08:31:17  hjj
bugfix, s/WORK,LWORK/WRK(KFREE),LFREE/ (bug found with ftnchek)

Revision 1.4  2000/05/01 13:22:04  hjj
removed obsolete KAVER test.

revision 1.2 + 1.3 April 2000 hjj
1)bug fix from bug made by hjaaj in aug 99:
  DIFDEN was not reset to false when threshold was changed.
2)corrected error in print statement
3)comments and polishing
4)restructured test of if new LUINDX file needs to be written
 (removed again IOPRH1 and FLAGSV in parameter list, this
 is more "clean").
=============================================================================
970326-hjaaj -- adjusted abort DIIS algorithm (allow more NEWOCC and more
  oscillations)
951130-hjaaj
DISCTL: moved 'reset to not converged' for solvent and writing SIRIFC to
        DISCTL from SIRCTL; it is a deficiency in the DIIS routines that
        it cannot be done in them. Look for MAERKE.
950824-hjaaj
DISCTL: output for AUTOCC
FCKEIG: also AUTOCC for open shell plus smaller changes.
9505-k.ruud
c: implemented AUTOCC option
950517-hjaaj
DISCTL: fixed MEMCHK error caused in RHFENR section; new call to READMO
940916-hjaaj
DISCTL: changed print of energy from F18.12 to F22.12
  (to avoid stars for e.g. Pb)
940701-hjaaj
SRDIIS: let MXERRV to default to MAX(2,MIN(NWOPT,10)) (was MXDIIS).
  This fixes linear dependency problem for very small debug runs:
  MXERRV=10 gave lin.dep. for 2 el. in 2 orbitals,
  MXERRV=2 gave quadratic convergence !
940511-hjaaj
SRDIIS,DISCTL: call flshfo so progress can be followed in output file;
   moved heading printing from DISCTL to SRDIIS
940505-hjaaj
DISCTL: exit if energy increasing or convergence is too slow.
940501-hjaaj
FCKEIG: exit if ".FREEZE" not compatible with ORDRSS.
    ORDRSS now always called (this also fixes error in prev. version).
940408-hjaaj
DISCTL: s/EMCACT/EACTIV/ to define active energy in /INFOPT/
940311-hjaaj
DISCTL: s/EVCNRM/GRDNRM/ (so gradient norm is transferred to SIROUT)
        nice output to LUW4 if LUW4 .ne. LUPRI
931126-hjaaj
DISCTL: call NEWORB after each iteration (for restart)
930720-hjaaj
SRDIIS: new parameter PROCC .FALSE. in CALL PRORB
930715-hjaaj
------: improved some print
SRDIIS: exit MXDGSS .GT. 0
930623-hjaaj
FCKEIG: inserted PARAMETER(D0=0.0D0)
#endif
C  /* Deck srdiis */
      SUBROUTINE SRDIIS (ICONV,CMO,WRK,KFRSAV,LFRSAV)
C
C Copyright (c) 1993 Hans Joergen Aa. Jensen 
C
C Last revision 1-Apr-1997 hjaaj
C
C
#include <implicit.h>
      DIMENSION CMO(*), WRK(*)
C
C
C Used from common blocks:
C  INFINP : MXDIIS,MAXEVC,DIRFCK
C  INFORB : N2BAST,NNBAST,N2BASX,NCMOT,NSYM, ...,MXDGSS
C  INFVAR : NWOPT
C  INFPRI : IPRRHF
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <infvar.h>
#include <infpri.h>
C
      PARAMETER (D4 = 4.0D0)
      LOGICAL C2DIIS
C
      CALL QENTER('SRDIIS')
      ICONV = 0
      IF (MXDGSS .GT. 1) THEN
         NINFO = NINFO + 1
         WRITE (LUPRI,'(/A/A)') ' SRDIIS INFO: no DIIS iterations',
     &'    because degenerate supersymmetries not implemented for DIIS'
         GO TO 9999
      END IF
C
      KFREE  = KFRSAV
      LFREE  = LFRSAV
C
CHJMAERKE: define IPRDIS separately?
C
      MXERRV = MAXEVC
      IF (MXERRV .LE. 0) THEN
         MXERRV = MIN(NWOPT,10)
C        ... fixes linear dependency problem for very small debug runs
C            MXERRV=10 gave lin.dep. for 2 el. in 2 orbitals,
C            MXERRV=2 gave quadratic convergence ! (940701-hjaaj)
         MXERRV = MAX(2,MXERRV)
      END IF
      IPRDIS = IPRRHF
      C2DIIS = FLAG(22)
C
      LUWOUT = LUPRI
  10  WRITE (LUWOUT,'(//A)')
     &      ' >>>>> DIIS optimization of Hartree-Fock <<<<<'
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
      IF (IPRDIS .GT. 25) THEN
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
C     scale by 4 to get correct orbital gradient from DISEVC
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
     &   '              because of memory limitations.'
         MXERRV = NEVC
      END IF
C
      LBMAT = (MXERRV+1)*MXERRV / 2
      CALL MEMGET('REAL',KBMAT ,LBMAT, WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KCVEC ,MXERRV+1,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KH1AO ,NNBAST,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KFSAV ,NDMAT*MXERRV*NNBAST,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KEVSAV,MXERRV*NNORBT,WRK,KFREE,LFREE)
      CALL DCOPY(NCMOT,CMO,1,WRK(KCMO1),1)
      CALL DISCTL(CMO,WRK(KCMO1),WRK(KBMAT),WRK(KH1AO),
     &            WRK(KFSAV),WRK(KEVSAV),MXERRV,C2DIIS,
     &            WRK(KSMOAO),WRK(KDCAO),WRK(KFCAO),WRK(KDSAV),
     &            WRK(KCVEC),WRK(KTNLM),WRK,KFREE,LFREE,IPRDIS,ICONV)
C
 9990 CALL MEMREL('SRDIIS',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
 9999 CALL FLSHFO(LUW4)
      IF (LUPRI .NE. LUW4) CALL FLSHFO(LUPRI)
      CALL QEXIT('SRDIIS')
      RETURN
      END
C  /* Deck disctl */
      SUBROUTINE DISCTL (CMO,CMO1,BMAT,H1AO,FAOSAV,EVCSAV,MXERRV,
     *                   C2DIIS,SMOAO,DCAO,FCAO,DAOSAV,CVEC,
     *                   TNLM,WRK,KFRSAV,LFRSAV,IPRDIS,ICONV)
C
C L.r. 5-May-1994 hjaaj / 1-Apr-1997
C
C DFT modifications T.  Helgaker
C
C
#include <implicit.h>
#include <dummy.h>
      DIMENSION CMO(*), CMO1(*), BMAT(*), H1AO(*)
      DIMENSION FAOSAV(NNBAST,MXERRV,2), EVCSAV(NNORBT,MXERRV)
      DIMENSION SMOAO(*), DCAO(*), FCAO(*), DAOSAV(NNBAST,*), CVEC(*)
      DIMENSION TNLM(*), WRK(*)
C
C Used from common blocks:
C  INFINP : THRRHF,MXDIIS
C  INFORB : NCMOT,NNBAST,NSYM
C  INFOPT : EPOT,EMCSCF,EMY,EACTIV,GRDNRM
C  INFVAR : NWOPH
C  SCBRHF : IOPRHF,AUTOCC
C
#include <maxorb.h>
#include <priunit.h>
#include <molde.h>
#include <infinp.h>
#include <inforb.h>
#include <infvar.h>
#include <infopt.h>
#include <scbrhf.h>
#include <cbihr2.h>
#include <infpri.h>
#include <gnrinf.h>
#include <dftcom.h>
#include <dfterg.h>
C
      PARAMETER (THDDEF = 1.0D-5)
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0, D1 = 1.0D0, DM1 = -1.0D0)
      PARAMETER (THREVC = 0.01D0, THRINC = 1.0D-10, CNVFAC = 0.8D0)
C
      LOGICAL   ONLYFC, C2DIIS, DIFDEN, DODIFDEN, NEWTHR, RHFWOP, OLDWOP
      DIMENSION MISH(8)
      CHARACTER*8 CMOLBL
C
      CALL QENTER('DISCTL')
      KFREE  = KFRSAV
      LFREE  = LFRSAV
      IFTHSV = IFTHRS
C
C     save current IOPRHF for check below
      IOPRH1 = IOPRHF
C
      IF (LUW4 .NE. LUPRI .OR. IPRDIS .LE. 1) THEN
      IF (FLAG(16)) THEN
      IF (.NOT.AUTOCC) THEN
         WRITE (LUW4,'(/A/A)')
     &      ' Iter      Total energy       Solvation energy    '//
     &      'Error norm    Delta(E)',
     &      ' -------------------------------------------------'//
     &      '----------------------'
      ELSE IF (IOPRHF .EQ. 0) THEN
         WRITE (LUW4,'(/A,I4,A//A/A)')
     &      ' Automatic occupation of symmetries with',NRHFEL,
     &      ' electrons.',
     &      ' Iter     Total energy      Solvation energy  Error norm'//
     &      '  Delta(E)    HF occupation',
     &      ' -----------------------------------------------------'//
     &      '-----------------------------'
      ELSE
         WRITE (LUW4,'(/A,I4,A/A//A/A)')
     &      ' Automatic occupation of symmetries with',NRHFEL,
     &      ' electrons.',
     &      ' "sym." is the symmetry of the one open shell orbital '//
     &      'in table below.',
     &      ' Iter     Total energy      Solvation energy  Error norm'//
     &      '  Delta(E)  sym.  closed shell occupation',
     &      ' -----------------------------------------------------'//
     &      '-------------------------------------------'
      END IF
      ELSE
      IF (.NOT.AUTOCC) THEN
         WRITE (LUW4,'(/A/A)')
     &      ' Iter      Total energy        Error norm    Delta(E)',
     &      ' ----------------------------------------------------'
      ELSE IF (IOPRHF .EQ. 0) THEN
         WRITE (LUW4,'(/A,I4,A//A/A)')
     &      ' Automatic occupation of symmetries with',NRHFEL,
     &      ' electrons.',
     &      ' Iter     Total energy    Error norm  Delta(E)  '//
     &      '  HF occupation',
     &      ' -----------------------------------------------------'//
     &      '------------------------'
      ELSE
         WRITE (LUW4,'(/A,I4,A/A//A/A)')
     &      ' Automatic occupation of symmetries with',NRHFEL,
     &      ' electrons.',
     &      ' "sym." is the symmetry of the one open shell orbital '//
     &      'in table below.',
     &      ' Iter     Total energy    Error norm  Delta(E)  '//
     &      'sym.  closed shell occupation',
     &      ' -----------------------------------------------------'//
     &      '------------------------'
      END IF
      END IF
      END IF
C
C     Get one-electron Hamiltonian
C
      CALL SIRH1(H1AO,WRK(KFREE),LFREE)
C
      ONLYFC = (NASHT .EQ. 0)
      IF (ONLYFC) THEN
         DV = D0
      ELSE IF (NASHT .EQ. 1) THEN
         DV = D1
      ELSE
         CALL QUIT('DISCTL called with NASHT .gt. 1')
      END IF
      THDIIS = THRRHF
      IF (THDIIS .LE. D0) THDIIS = THDDEF
      DAMP   = D0
      EMCSCF = D0
      EMCLOW = D0
      GNRMSV = D0
      ITDIIS = 0
      ITNOCC = 0
      NEWOCC = 0
      NOSCIL = 0
      DODIFDEN = .NOT. FLAG(49) .AND. .NOT. FLAG(16) .AND. .NOT.DFTRUN
Chjaaj: DIFDEN does not work with solvent
C       because FC_sol is added to the old FC matrix, but not the new!
Ctuh:   DIFDEN does not work with DFT either   
      DIFDEN   = DODIFDEN
      IF (.NOT. USRSCR) THEN
         IF (FLAG(13)) THEN
Chj-aug99: MO's read from input, assume good guess
            IFTHRS = 14
         ELSE
Chj-aug99: MO's from huckel or H1DIAG, assume poor guess
            IFTHRS = 7
         END IF
      END IF
      IF (ICI0 .EQ. 6) THEN
Chj      ... if (GEOWLK) then
         EMCGEO = EMCOLD
         DEPGEO = DEPRED
      END IF
  100 CONTINUE
         ITDIIS = ITDIIS + 1
         ITNOCC = ITNOCC + 1
         INDEVC = MOD(ITDIIS-1,MXERRV) + 1
         NEVC   = MIN(MXERRV,ITDIIS)
C
         IF (DIRFCK .AND. (IPRDIS .GE. 0 .OR. IPRDIS .LE. 1)) THEN
            WRITE (LUPRI,'(A,I5,L5)') 
     &          ' Screening settings in next iteration (IFTHRS/DIFDEN)',
     &           IFTHRS, DIFDEN
            CALL FLSHFO(LUPRI)
         END IF
C
C  ***** Construct folded inactive and active parts of the one-electron
C        density matrix in AO-basis (DCAO and DVAO). DV is the active
C        part of the one-electron density matrix in MO-basis.
C
         CALL FCKDEN((NISHT.GT.0),(.NOT.ONLYFC),
     *                DCAO,DCAO(1+N2BASX),CMO,DV,WRK(KFREE),LFREE)
C        CALL FCKDEN(GETDC,GETDV,DCAO,DVAO,CMO,DV,WRK,LWRK)
C
         CALL MEMGET('REAL',KTMP1,NNBASX,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KTMP2,N2BASX,WRK,KFREE,LFREE)
         IF (NISHT.GT.0) THEN
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
         CALL MEMREL('DISCTL.DIFDEN',WRK,KFRSAV,KTMP1,KFREE,LFREE)
C
C  ***** Calculate 2-electron part of Fock matrix with
C        differential density
C
         CALL FCK2AO(ONLYFC,FCAO,FCAO(1+N2BASX),
     *               DCAO,DCAO(1+N2BASX),WRK,KFREE,LFREE)
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
            CALL MEMREL('DISCTL.FAOSAV',WRK,KFRSAV,KTMP1,KFREE,LFREE)
         ELSE
            CALL DGETSP(NBAST,FCAO,FAOSAV(1,INDEVC,1))
            IF (.NOT.ONLYFC) THEN
               CALL DGETSP(NBAST,FCAO(1+N2BASX),FAOSAV(1,INDEVC,2))
            END IF
         END IF
         IF (DIFDEN .AND. NEVC .GT. 1) THEN
            INDEVP = MOD(ITDIIS-2,MXERRV) + 1
            CALL DAXPY(NNBAST,DCOVLP,FAOSAV(1,INDEVP,1),1,
     &                               FAOSAV(1,INDEVC,1),1)
C           .. add DCOVLP*(FCAO_previous + FVAO_previous)
            IF (.NOT.ONLYFC) THEN
               CALL DAXPY(NNBAST,-DCOVLP,FAOSAV(1,INDEVP,2),1,
     &                                   FAOSAV(1,INDEVC,1),1)
C           .. subtract DCOVLP*FVAO_previous to get DCOVLP*FCAO_previous
               CALL DAXPY(NNBAST, DVOVLP,FAOSAV(1,INDEVP,2),1,
     &                                   FAOSAV(1,INDEVC,2),1)
            END IF
         END IF
         EMY = DDOT(NNBAST,DAOSAV,1,H1AO,1)
     &       + DDOT(NNBAST,DAOSAV,1,FAOSAV(1,INDEVC,1),1) * DP5
     &       + EDFTY
C
         IF (ONLYFC) THEN
            EACTIV = D0
         ELSE
C           EACTIV = DV(u,u)FC(u,u) = FC(u,u)
C                  = DVAO(p,q)FC2AO(p,q) + DVAO(p,q)H1AO(p,q)
            EACTIV = DDOT(NNBAS(IOPRHF),DAOSAV(1+IIBAS(IOPRHF),2),1,
     &                           FAOSAV(1+IIBAS(IOPRHF),INDEVC,1),1)
     &             + DDOT(NNBAS(IOPRHF),DAOSAV(1+IIBAS(IOPRHF),2),1,
     &                           H1AO(1+IIBAS(IOPRHF)),1)
            CALL DAXPY(NNBAST,D1,FAOSAV(1,INDEVC,2),1,
     &                 FAOSAV(1,INDEVC,1),1)
         END IF
C
         ESOLT = D0
         IF (FLAG(16)) THEN
            CALL MEMGET('REAL',KFCSOL,NNBAST,WRK,KFREE,LFREE)
            CALL SOLFCK(DAOSAV(1,1),DAOSAV(1,2),
     &                  WRK(KFCSOL),DUMMY,TNLM,.FALSE.,ESOLT,
     &                  WRK(KFREE),LFREE,IPRDIS)
            CALL DAXPY(NNBAST,D1,WRK(KFCSOL),1,FAOSAV(1,INDEVC,1),1)
            CALL MEMREL('DISCTL.SOLFCK',WRK,KFRSAV,KFCSOL,KFREE,LFREE)
         END IF
C
         EMCDIF = EMCSCF
         EMCSCF = EMY + EACTIV + EPOT + ESOLT
         EMCDIF = EMCSCF - EMCDIF
C
         CALL DISEVC(SMOAO,DAOSAV(1,1),DAOSAV(1,2),H1AO,
     &               FAOSAV(1,INDEVC,1),FAOSAV(1,INDEVC,2),CMO1,
     &               EVCSAV(1,INDEVC),WRK,KFREE,LFREE,IPRDIS)
C        CALL DISEVC(SMOAO,DCAO,DVAO,H1AO,FDAO,FVAO,CMO1,
C    &               ERRVEC,WRK,KFRSAV,LFRSAV,IPRDIS)
         IF (IPRDIS .GT. 20) THEN
            WRITE (LUPRI,*) ' Test output of EVCSAV:'
            CALL OUTPUT(EVCSAV,1,NNORBT,1,NEVC,NNORBT,NEVC,1,LUPRI)
         ELSE IF (IPRDIS .GT. 10) THEN
            WRITE (LUPRI,*) ' Error vector, DIIS iteration',ITDIIS
            CALL OUTPKB(EVCSAV(1,INDEVC),NORB,NSYM,1,LUPRI)
         END IF
         GRDNRM = DNRM2(NNORBT,EVCSAV(1,INDEVC),1)
         IF (LUW4 .NE. LUPRI .OR. IPRDIS .GT. 1) THEN
            WRITE (LUPRI,'(/A/A,I5,F22.12,1P,D15.5,D12.2,I5,L5)')
     &' ITDIIS, energy, error norm, change in energy, scr.thr., difden:'
     & ,' DIIS iter.',ITDIIS,EMCSCF,GRDNRM,EMCDIF,-IFTHRS,DIFDEN
            IF (AUTOCC) THEN
               WRITE (LUPRI,'(/A,8I4)')
     &         ' AUTOCC: current HF occupation',(NISH(I),I=1,NSYM)
               IF (IOPRHF .NE. 0) WRITE (LUPRI,'(A,I4)')
     &         ' AUTOCC: current open shell symmetry',IOPRHF
            END IF
            CALL FLSHFO(LUPRI)
         END IF
         IF (LUW4 .NE. LUPRI .OR. IPRDIS .LE. 1) THEN
            IF (FLAG(16)) THEN
               IF (.NOT.AUTOCC) THEN
                  WRITE (LUW4,'(I5,F20.12,F20.12,1P,D15.5,D12.2)')
     &               ITDIIS,EMCSCF,ESOLT,GRDNRM,EMCDIF
               ELSE IF (IOPRHF .EQ. 0) THEN
                 WRITE (LUW4,'(I4,F20.12,F20.12,1P,D11.2,D11.2,2X,8I3)')
     &                 ITDIIS,EMCSCF,ESOLT,GRDNRM,EMCDIF,
     &                 (NISH(I),I=1,NSYM)
               ELSE
                 WRITE (LUW4,'(I4,F20.12,F20.12,1P,D11.2,D11.2,
     &                 I4,2X,8I3)')
     &                 ITDIIS,EMCSCF,ESOLT,GRDNRM,EMCDIF,IOPRHF,
     &                 (NISH(I),I=1,NSYM)
               END IF
            ELSE
               IF (.NOT.AUTOCC) THEN
                  WRITE (LUW4,'(I5,F22.12,1P,D15.5,D12.2)')
     &                 ITDIIS,EMCSCF,GRDNRM,EMCDIF
               ELSE IF (IOPRHF .EQ. 0) THEN
                  WRITE (LUW4,'(I4,F20.12,1P,D11.2,D11.2,2X,8I3)')
     &                 ITDIIS,EMCSCF,GRDNRM,EMCDIF,(NISH(I),I=1,NSYM)
               ELSE
                  WRITE (LUW4,'(I4,F20.12,1P,D11.2,D11.2,I4,2X,8I3)')
     &                 ITDIIS,EMCSCF,GRDNRM,EMCDIF,IOPRHF,
     &                 (NISH(I),I=1,NSYM)
               END IF
            END IF
            CALL FLSHFO(LUW4)
         END IF
         IF (MOLDEN) CALL MOSCFCON(ITDIIS,EMCSCF,.FALSE.)
C
C        Test and exit if energy is increasing
C        (from itdiis .eq. 7, second iteration is Roothaan step
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
         EMCDIF = EMCSCF - EMCLOW
         IF (EMCDIF .GT. THRINC*ABS(EMCLOW) .AND. ITNOCC .GT. 4) THEN
C            ... if increase < THRINC*abs(EMCLOW) then the increase
C                are in digits which are not reliable
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
         EMCLOW = EMCSCF
C
 4010 FORMAT(' DIIS aborted because energy increased in this iteration.'
     &      /' Program aborted because this means automatic symmetry'
     &      /' designation of occupied orbitals is not possible!',
     &      /' Please send this output to hjj@dou.dk !!!',
     &      /' You may try again with explicit specfication of .HF OCC')
 4020 FORMAT(' DIIS aborted because energy increased in this iteration!'
     &      /' Backstep to orbitals from previous iteration.')
C
C        Iteration accepted; save current CMO for restart
         WRITE (CMOLBL,'(A4,I4)') 'DIIS',ITDIIS
         CALL NEWORB(CMOLBL,CMO,.FALSE.)
C        ... REWIT1 false: do not destroy any GEOWALK information
C
         IF (GRDNRM .LE. THDIIS) THEN
            IF (MOLDEN) CALL MOSCFCON(ITDIIS,EMCSCF,.TRUE.)
            ICONV = 1
            WRITE (LUPRI,4110) ITDIIS
            IF (LUW4 .NE. LUPRI) WRITE (LUW4,4110) ITDIIS
            GO TO 1000
         ELSE IF (ITDIIS .GE. MXDIIS) THEN
            ICONV = 0
            WRITE (LUPRI,4120)
            IF (LUW4 .NE. LUPRI) WRITE (LUW4,4120)
            GO TO 1100
         END IF
 4110 FORMAT(' DIIS converged in',I4,' iterations !')
 4120 FORMAT(' DIIS aborted because max DIIS iterations reached !')
C
Chj 990819: use DIFDEN more often
Chj was: IF (IFTHRS .GE. 12 .AND. .NOT. FLAG(49)) DIFDEN = .TRUE.
         NEWTHR = .FALSE.
         IF (DIRFCK .AND. .NOT. USRSCR) THEN
            ITHRS = NINT(-LOG10(GRDNRM))+9
            IF (DODIFDEN) THEN
               ILIM = 2
Chj-aug99: only NEWTHR with skip of 2 for DIFDEN
               IF (ITHRS .GE. 12) ITHRS = 14
Chj-aug99: for effective DIFDEN; note this test is OK w. ILIM.eq.2
            ELSE
               ILIM = 1
            END IF
            IF (ABS(ITHRS-IFTHRS) .GE. ILIM) THEN 
               NEWTHR = .TRUE.
               IFTHRS = MIN(14,MAX(7,ITHRS))
            END IF
         END IF
Chj 990819/000103 new:
         DIFDEN = .NOT. NEWTHR .AND. DODIFDEN
C        DIFDEN true if not new threshold and not "no difden"
C        -- if we have tightened screening threshold we must
C           calculate w/o DIFDEN or the error in the Fock matrix
C           would correspond to the previous screening threshold.
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
 4140 FORMAT(' DIIS aborted, convergence too slow !')
C
         GNRMSV = GRDNRM
C
         CALL DISRED(C2DIIS,BMAT,INDEVC,NEVC,DAMP,THREVC,
     &               CVEC,XMLBDA,EVCSAV,WRK,KFREE,LFREE,IPRDIS)
C        CALL DISRED(C2DIIS,BMAT,INDEVC,NEVC,DAMP,THREVC,
C    *               CVEC,XLMBDA,ERRVEC,WRK,KFRSAV,LFRSAV,IPRDIS)
C
         CALL DCOPY(NNBAST,H1AO,1,FCAO,1)
         CALL DGEMM('N','N',NNBAST,1,NEVC,1.D0,
     &              FAOSAV(1,1,1),NNBAST,
     &              CVEC,NEVC,1.D0,
     &              FCAO,NNBAST)
         IF (.NOT. ONLYFC) THEN
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
         IF (.NOT. ONLYFC) THEN
         DO 300 ISYM = 1,NSYM
            IF (NASH(ISYM) .EQ. 0) GO TO 300
            IF (IPRDIS .GE. 15) THEN
               WRITE (LUPRI,'(/A,I3)')
     &            ' Fock matrix before FV correction, symmetry',ISYM
               CALL OUTPAK(DCAO(1+IIORB(ISYM)),NORB(ISYM),1,LUPRI)
            END IF
            CALL UTHU(FCAO(1+NNBAST+IIBAS(ISYM)),DCAO(1+NNBAST),
     &         CMO(1+ICMO(ISYM)),WRK(KFREE),NBAS(ISYM),NORB(ISYM))
            KACT = NISH(ISYM) + 1
            KROW = KACT*(KACT-1)/2
            DO 310 J = 1,KACT-1
               JFKJ = IIORB(ISYM) + KROW + J
               DCAO(JFKJ) = DCAO(JFKJ) + DCAO(NNBAST+KROW+J)
  310       CONTINUE
            DO 320 J = KACT+1,NORB(ISYM)
               JROW = J*(J-1)/2
               JFJK = IIORB(ISYM) + JROW + KACT
               DCAO(JFJK) = DCAO(JFJK) - DCAO(NNBAST+JROW+KACT)
  320       CONTINUE
            IF (IPRDIS .GE. 15) THEN
               WRITE (LUPRI,'(/A,I3)')
     &            ' Fock matrix after FV correction, symmetry',ISYM
               CALL OUTPAK(DCAO(1+IIORB(ISYM)),NORB(ISYM),1,LUPRI)
            END IF
  300    CONTINUE
         END IF
C        Diagonalize FC (saved in DCAO)
         CALL ICOPY(8,NISH,1,MISH,1)
         JOPRHF = IOPRHF
         CALL FCKEIG(CMO,DCAO,WRK(KFREE),LFREE)
C        test if AUTOCC has changed occupation
         NTEST = 0
         DO ISYM = 1,8
            NTEST = NTEST + ABS(NISH(ISYM)-MISH(ISYM))
         END DO
         NTEST = NTEST + ABS(IOPRHF-JOPRHF)
         IF (NTEST .GT. 0) THEN
            NEWOCC = NEWOCC + 1
            ITNOCC = 0
            NOSCIL = 0
            IF (NEWOCC .GT. 6) THEN
               WRITE (LUPRI,4210) NEWOCC
               IF (LUW4 .NE. LUPRI) WRITE (LUW4,4210) NEWOCC
               CALL QUIT(
     &           'FATAL ERROR: problems with automatic HF occupation')
            END IF
Chj-sep99: bugfix for open shell when open shell symmetry changed
C          (NASH(ISYM) used for test in DISCTL !!!)
            IF (NASHT .EQ. 1) THEN
               CALL IZERO(NASH,8)
               NASH(IOPRHF) = 1
            END IF
Chj-sep99: and change index arrays to new occupation (for RHFENR and ?)
            CALL SETORB
         END IF
         GO TO 100
 4210 FORMAT(/' DIIS has changed occupation numbers',I2,' times now.'
     &      /' Program aborts because this indicates problems with'//
     &       ' automatic occupation.')
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
C     because ICONV = 1 and if either IPRDIS .gt. 0 or
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
      DO 11 ISYM = 1, NSYM
         IF (NRHF(ISYM) .NE. NISH(ISYM)) THEN
            NRHF(ISYM) = NISH(ISYM)
            RHFWOP = .FALSE.
         END IF
 11   CONTINUE
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
      IPRENR = 0
      IF (ONLYFC .AND. ICONV .NE. 0. .AND. (IPRDIS .GT. 0 .OR.
     &    .NOT.(DOMP2 .OR. DOCI .OR. DOMC)) ) IPRENR = 1
      IWRIFC = 0
Chj aug99:
      IF ( ( FLAG(25) .OR. INERSI ) .AND.
     &     (ICONV .NE. 0. .AND. .NOT. DOMC) ) THEN
         IWRIFC = 1
Chj   solvnt (flag(16)) not implemented for flag(25) here yet !!! aug 99
Chj MAERKE
      END IF
C
      IF (IPRENR .GT. 0 .OR. IWRIFC .GT. 0) THEN
C
C        Transform FDAO to FD (FD saved in FCAO).
C        Zero FV (in DCAO)
C
         CALL DAXPY(NNBAST,D1,H1AO,1,FAOSAV(1,INDEVC,1),1)
         CALL UTHUB(FAOSAV(1,INDEVC,1),FCAO,CMO,WRK(KFREE),
     &              NSYM,NBAS,NORB)
         CALL DZERO(DCAO,NNORBT)
Chj from siropt.F, should be modified
       END IF
C
       IF (IPRENR .GT. 0) THEN
         WRITE (LUW4,'(//A)') ' *** RHF orbital energy analysis ***'
         IF (FLAG(16))
     &     WRITE (LUW4,'(A)') '    (incl. solvent contribution)'
         CALL RHFENR(IPRI4,LUW4,FCAO,DCAO, WRK(KFREE),LFREE)
Chj      CALL RHFENR(IPRI4,LUW4,WRK(KFCSOL),WRK(KFV),
Chj  &               WRK(KFREE),LFREE)
         IF (LUPRI .NE. LUW4) THEN
            WRITE (LUPRI,'(//A)') ' *** RHF orbital energy analysis ***'
            IF (FLAG(16))
     &        WRITE (LUPRI,'(A)') '    (incl. solvent contribution)'
            CALL RHFENR(IPRI6,LUPRI,FCAO,DCAO, WRK(KFREE),LFREE)
         END IF
C        CALL RHFENR(IPRINT,LUPRI,FC,FV,SCRA,LSCRA)
      END IF
C
C     If (FLAG(25) .OR. INERSI) write LUSIFC
C       FLAG(25) is true for .ABACUS (geometry optimization)
C                         or .RESPONS (response calculation)
C                         or .INTERFACE (request for LUSIFC)
C       INERSI is true if this is initial state calculation
C              in a solvent calculation with inertial polarization
C
      IF ( IWRIFC .EQ. 1 ) THEN
         IF (ICI0 .EQ. 6) THEN
Chj      ... if (GEOWLK) then
            IF (.NOT. OPTNEW)
     *         WRITE (LUW4,'(//A/)')
     *            ' Check ratio for geometry walk with converged RHF :'
            JWSTEP = 1
            EMCOLD = EMCGEO
            DEPRED = DEPGEO
            CALL SIRSTP(JWSTEP,DUMMY,DUMMY,DUMMY,DUMMY,1)
Chj         ... WLKREJ in gnrinf.h will be true if geom. step rejected.
         END IF
Chj      make PV(1) = D0 for open shell RHF
         CALL MEMGET('REAL',KPV,1,WRK,KFREE,LFREE)
         WRK(KPV) = D0
         CALL MEMGET('REAL',KGORB,NWOPH,WRK,KFREE,LFREE)
         CALL DZERO(WRK(KGORB),NWOPH)
Chj      ... set orbital gradient to zero (i.e. completely converged!)
C
         IF (NASHT.GT.0) THEN
C           put FV in DCAO
            CALL UTHUB(FAOSAV(1,INDEVC,2),DCAO,CMO,WRK(KFREE),
     &                 NSYM,NBAS,NORB)
C           form FC = FD - FV in FCAO
            CALL DAXPY(NNORBT,DM1,DCAO,1,FCAO,1)
            CALL MEMGET('REAL',KDUM,NASHT*NORBT,WRK,KFREE,LFREE)
            CALL DZERO(WRK(KDUM),NASHT*NORBT)
         ELSE
            CALL MEMGET('REAL',KDUM,0,WRK,KFREE,LFREE)
         END IF
C
         CALL WRSIFC(CMO,DV,WRK(KPV),FCAO,DCAO,
     &               WRK(KDUM),WRK(KDUM),WRK(KDUM),WRK(KDUM),
     &               WRK(KFREE),LFREE,WRK(KGORB),TNLM,.FALSE.,WRK(KDUM))
C        CALL WRSIFC(CMO,DV,PV,FC,FV,FQ,CREF,FCAC,H2AC,WRK,LFREE,
C    *               GORB,ERLM,ORBHES,XINDX)
C
         CALL MEMREL('DISCTL.WRSIFC',WRK,KFRSAV,KPV,KFREE,LFREE)
      END IF
C
C        Transform integrals for abacus or response, 
C        if converged and if not followed by higher level.
C
      IF (ICONV .NE. 0 .AND. ABS(ITRFIN) .LE. 10 .AND.
     &    .NOT.(DOMP2 .OR. DOCI .OR. DOMC) ) THEN
         IF (IPRI6 .GE. 5) WRITE (LUPRI,'(/2A,I3,A)')
     *      ' >>> Transforming integrals for Abacus/response with',
     *      ' ITRLVL =',ITRFIN,' <<<'
         JTRLVL = ITRFIN
         KCMO   = KFREE
         KWRK10 = KCMO  + NCMOT
         LWRK10 = LFREE + 1 - KWRK10
         CALL TRACTL(JTRLVL,WRK(KCMO),WRK(KWRK10),LWRK10,MWORK)
      END IF
C
C ***********************************************************
C
 1100 CONTINUE
C     hjaaj-aug99: restore original IFTHRS value
      IFTHRS = IFTHSV
      CALL QEXIT('DISCTL')
      RETURN
      END
C  /* Deck disevc */
      SUBROUTINE DISEVC(SMOAO,DCAO,DVAO,H1AO,FDAO,FVAO,CMO1,
     &                  ERRVEC,WRK,KFRSAV,LFRSAV,IPRDIS)
C
C
#include <implicit.h>
      DIMENSION DCAO(NNBAST), DVAO(NNBAST), FDAO(NNBAST), FVAO(NNBAST)
      DIMENSION H1AO(NNBAST)
      DIMENSION SMOAO(N2BAST), CMO1(NCMOT), ERRVEC(NNORBT), WRK(*)
C
C Used from common blocks:
C  INFINP : LNOROT,NOROT()
C  INFORB : N2BAST, NNBAST, NCMOT, NSYM, ..., NFROT, ...
C  INFDIM : N2BASM
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <infdim.h>
#include <infpri.h>
C
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DM1 = -1.0D0)
C
      CALL QENTER('DISEVC')
      KFREE = KFRSAV
      LFREE = LFRSAV
      CALL MEMGET('REAL',KTMP1,N2BASM,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KTMP2,N2BASM,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KTMP3,N2BASM,WRK,KFREE,LFREE)
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
         IF (IPRDIS .GE. 35) THEN
            WRITE (LUPRI,*) 'Debug output from DISEVC, ISYM=',ISYM
            WRITE (LUPRI,*) 'NORBI,NBASI:',NORBI,NBASI
            WRITE (LUPRI,*) 'DISEVC CMO1:'
            CALL OUTPUT(CMO1(1+ICMO(ISYM)),1,NBASI,1,NORBI,
     &                  NBASI,NORBI,1,LUPRI)
            WRITE (LUPRI,*) 'DISEVC SMOAO = 4 Ct SAO:'
            CALL OUTPUT(SMOAO(1+I2BAS(ISYM)),1,NORBI,1,NBASI,
     &                  NORBI,NBASI,1,LUPRI)
         IF (IPRDIS .GT. 50) THEN
            CALL DGEMM('N','N',NORBI,NORBI,NBASI,1.D0,
     &                 SMOAO(1+I2BAS(ISYM)),NORBI,
     &                 CMO1(1+ICMO(ISYM)),NBASI,0.D0,
     &                 WRK(KTMP3),NORBI)
            WRITE (LUPRI,*) 'DISEVC SMO = 4 Ct SAO C:'
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
            IF (IPRDIS .GE. 35) THEN
               WRITE (LUPRI,*) 'DISEVC DCAO unfolded:'
               CALL OUTPUT(WRK(KTMP2),1,NBASI,1,NBASI,
     &                     NBASI,NBASI,1,LUPRI)
               WRITE (LUPRI,*) 'DISEVC FDAO:'
               CALL OUTPUT(WRK(KTMP3),1,NBASI,1,NBASI,
     &                     NBASI,NBASI,1,LUPRI)
               WRITE (LUPRI,*) 'DISEVC TMP1 = DCAO FDAO:'
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
            IF (IPRDIS .GE. 35) THEN
               WRITE (LUPRI,*) 'DISEVC DVAO unfolded:'
               CALL OUTPUT(WRK(KTMP2),1,NBASI,1,NBASI,
     &                     NBASI,NBASI,1,LUPRI)
               WRITE (LUPRI,*) 'DISEVC FCAO:'
               CALL OUTPUT(WRK(KTMP3),1,NBASI,1,NBASI,
     &                     NBASI,NBASI,1,LUPRI)
               WRITE (LUPRI,*) 'DISEVC TMP1 = DCAO FDAO + DVAO FCAO:'
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
         IF (IPRDIS .GE. 15) THEN
            WRITE (LUPRI,'(/A/A/A,I3)')
     &      ' (DISEVC) error vector = orbital gradient vector =',
     &      '          4 (DCmo FDmo + DVmo FCmo) = ',
     &      '          4 Ct SAO (DCAO FDAO + DVAO FCAO) C for symmetry',
     &      ISYM
            CALL OUTPUT(WRK(KTMP1),1,NORBI,1,NORBI,NORBI,NORBI,1,LUPRI)
         END IF
  800 CONTINUE
      CALL MEMREL('DISEVC',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('DISEVC')
      RETURN
      END
C  /* Deck disred */
      SUBROUTINE DISRED (C2DIIS,BMAT,INDEVC,NEVC,DAMP,THREVC,
     *                   CVEC,XLMBDA,ERRVEC,WRK,KFRSAV,LFRSAV,IPRDIS)
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
C          Threshold for acceptable solution in C2-DIIS
C
C Output:
C  BMAT,  the new, extended reduced orbital Hessian-matrix
C  CVEC  :.not.C2DIIS SOLUTION TO LINEAR SET OF EQUATIONS
C        :     C2DIIS LOWEST acceptable EIGENVECTOR
C
#include <implicit.h>
      DIMENSION BMAT(*),CVEC(*), ERRVEC(NNORBT,*), WRK(*)
      LOGICAL C2DIIS
C
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
C
C Used from common blocks:
C   INFORB : NNORBT
C
#include <priunit.h>
#include <inforb.h>
#include <infpri.h>
C
      IROW(I) = I*(I-1)/2
C
      CALL QENTER('DISRED')
      KFREE = KFRSAV
      LFREE = LFRSAV
C
      IF (IPRDIS .GE. 9) WRITE (LUPRI,7555) C2DIIS,NEVC
 7555 FORMAT(/' (DISRED) CONTROL PARAMETERS C2DIIS =',L10,
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
      IF ( IPRDIS.GE.5 ) THEN
         WRITE (LUPRI,'(//A)') ' (DISRED) B matrix:'
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
      IF ( IPRDIS.GE.3 .OR. IOK.NE.1) THEN
         WRITE(LUPRI,'(//A,I5)')
     *      ' (DISRED) C2-DIIS B matrix eigenvalues. IOK =',IOK
         CALL OUTPUT(WRK(KBTMP),1,1,1,NEVC,1,NEVC,1,LUPRI)
      END IF
      IF ( IPRDIS.GE.4) THEN
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
         CALL QUIT('DISRED: DAMP not implemented yet')
#if defined (VAR_KAPREDCODE)
  >>>>>>>>>>>>>>>> NOTE: Hamilton and Pulay uses '* (D1 + DAMP)'
                                        and not  ' + DAMP'
         DO 510 I = 1,NEVC
            II = KBTMP - 1 + IROW(I+1)
            WRK(II) = WRK(II) + DAMP
 510     CONTINUE
         IF ( IPRDIS.GE.3 ) THEN
            WRITE(LUPRI,'(A,1P,D15.8)')
     *      ' (DISRED) B matrix damping factor =',DAMP
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
      IF ( IPRDIS.GE.3 ) THEN
         WRITE(LUPRI,'(/A)')
     *      ' (DISRED) Solutions to C1-DIIS set of linear equations'
         CALL OUTPUT(CVEC,1,NLEQ,1,1,NLEQ,1,1,LUPRI)
      END IF
      IF (INFO.NE.0) THEN
         WRITE(LUPRI,8500) INFO
         CALL QTRACE(LUPRI)
         CALL QUIT('DISRED: no solution to linear equations')
      END IF
 8500 FORMAT(/' (DISRED) Solution not obtained to linear equations'
     *       /T11,'Check if matrix is singular, LINPACK DSP INFO =',I3)
C
C *** End of subroutine DISRED
C
 9999 CONTINUE
      IF (IPRDIS .GT. 1) THEN
         WRITE (LUPRI,'(/A,L2,A,1P,D10.2)')
     &      ' DIIS C vector; C2DIIS =',C2DIIS,' LAMBDA =',XLMBDA
         WRITE (LUPRI,'(6F12.6)') (CVEC(I),I=1,NEVC)
      END IF
      CALL MEMREL('DISRED',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('DISRED')
      RETURN
      END
C  /* Deck fckeig */
      SUBROUTINE FCKEIG(CMO,FC,SCRA,LSCRA)
C
C Written 20-May-1993 by Hans Jorgen Aa. Jensen
C
C Purpose:
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
#include <implicit.h>
      DIMENSION CMO(*),FC(*),SCRA(*)
C
C Used from common blocks:
C  INFINP : SUPSYM, LNOROT, NOROT()
C  SCBRHF : NFRRHF(*), AUTOCC
C  INFIND : IROW(*),...,ISSMO(*),?
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <scbrhf.h>
#include <infind.h>
#include <infpri.h>
C
      LOGICAL LSAVE4,LSAVE6
      PARAMETER (D0 = 0.0D0)
C
      CALL QENTER('FCKEIG')
C
      LSAVE4 = P4FLAG(9)
      LSAVE6 = P6FLAG(6)
      P4FLAG(9) = .FALSE.
      P6FLAG(6) = .FALSE.
C
C     Some memory allocation for HF occupation determination
C
      KEIG  = 1
      KSYMS = KEIG  + NORBT
      KLAST = KSYMS + NORBT
C
C     Step 1: Diagonalize Fock-matrix:
C
      DO 200 ISYM = 1,NSYM
         NORBI = NORB(ISYM)
      IF (NORBI.EQ.0) GO TO 200
         IORBI = IORB(ISYM)
         NBASI = NBAS(ISYM)
         NFRZI = NFRRHF(ISYM)
         IFSYM = IIORB(ISYM) + 1
         ICSYM = ICMO(ISYM) + 1
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
                  FC(IFSYM-1+KROW+L) = D0
  110          CONTINUE
               DO 120 L = K+1,NORBI
                  KL = (IFSYM-1+K)+IROW(L)
                  FC(KL) = D0
  120          CONTINUE
            ELSE
               JNFRZ = JNFRZ + 1
            END IF
  130       CONTINUE
         END IF
C
         CALL JACO(FC(IFSYM),CMO(ICSYM),NORBI,NORBI,NBASI,
     *             SCRA(KLAST),SCRA(KLAST+NORBI))
C        CALL JACO (F,VEC,NB,NMAX,NROWV,BIG,JBIG)
C
         DO 175 I=1,NORBI
            FC(IORBI+I)=FC(IFSYM-1+IROW(I+1))
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
  200 CONTINUE
C
C     Step 2: Reorthogonalize new mo's
C
      KSAO  = KLAST
      KSCR1 = KSAO + NNBAST
      LSCR1 = LSCRA - KSCR1
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
