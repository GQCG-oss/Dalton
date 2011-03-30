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
C  /* Deck sirctl */
      SUBROUTINE SIRCTL(ICONV,WORK,LWORK)
C
C *** This is SIRIUS, a direct second order MCSCF program ***
C
C  (c) Copyright Hans Joergen Aa. Jensen and Hans Agren
C
C  Last revision 1-July-1994/Nov 1996 hjaaj
C
C  Original version
C  written in Uppsala late `83 and Aarhus early `84 by
C
C     Hans Joergen Aa. Jensen
C
C   Theoretical Division
C   Institute of Chemistry
C   University of Aarhus
C   DK-8000 Aarhus C
C   Denmark
C
C     Hans Agren
C
C   Institute of Quantum Chemistry
C   University of Uppsala
C   Box 518, S-75120 Uppsala
C   Sweden
C
C

      use lucita_orbital_spaces

#include <implicit.h>
C
      DIMENSION WORK(LWORK)
#include <iratdef.h>
#include <dummy.h>
#include <thrldp.h>
C
C Used from common blocks:
C   INFINP : DO*, MAXMAC
C   INFORB : NCMOT,
C   INFVAR : NCONF,JWOPSY,?
C   INFTRA : USEDRC
C   PRIUNIT : IPRSTAT, NWARN
C   SCBRHF : RHFCAN, MXHFMA, MXHFMI, NFRRHF(8), IOPRHF, MAXFCK, MXDIIS, THRRHF
C   GNRINF : WRINDX
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <infvar.h>
#include <inftap.h>
#include <inftra.h>
#include <infpri.h>
#include <scbrhf.h>
#include <gnrinf.h>
#include <mxcent.h>
#include <pcmdef.h>
#include <pcm.h>
#include <pcmlog.h>
! dftcom.h: DFT_SPINDNS, DFT_LOCALSPIN
#include <dftcom.h>
C
      LOGICAL   FLAGSV(NFLAG), OLDWOP, RHFWOP, PCMSV
      LOGICAL   DFT_SPINDNS_SAVE, DFT_LOCALSPIN_SAVE
      DIMENSION NASHSV(8), NISHSV(8), NFROSV(8)
      real(8), allocatable :: tmp_space(:)
 
#include <sirbkd.h>
C
C *** Call setup to define limits and I/O units
C     and maybe open some files (depending on host)
C
      CALL QENTER('SIRCTL')
      CALL GETTIM(TSTART,WSTART)
      TIMSIR = TSTART
      WALSIR = WSTART
C
C *** Input section
C     SIRINP processes input from luinp
C
      CALL SIRINP(WORK,LWORK)
      CALL GETTIM(TEND,WEND)
      TINP = TEND - TSTART
      WINP = WEND - WSTART
      IF (IPRSTAT .GT. 0) THEN
         WRITE (LUSTAT,'(/A,2F12.3)')
     &      ' >>>> CPU and wall time for SIRINP :',TINP,WINP
      END IF
      CALL FLSHFO(LUW4)
      CALL FLSHFO(LUPRI)
      CALL FLSHFO(LUSTAT)
C
C ***************************************************
C *** Optimize for RHF, MP2, CI, MCSCF wave functions
C ***************************************************
C
C     (ICONV = 1 if converged, otherwise ICONV = 0)
C
      ICONV = -1
      NUMRUN = 0
C
C To prepare combination runs: Save flags for MCSCF. Nullify flags
C for RHF, MP2, CI
C
C - save
      DO 2 IFLAG = 1, NFLAG
         FLAGSV(IFLAG) = FLAG(IFLAG)
    2 CONTINUE
      MCTPSV = MCTYPE
C - change
C     Now disable non-symmetric IWOPSY for properties
C     (only wanted in first call of SETWOP)
      IF (DOMC) FLAG(11) = .FALSE.
C     ... flag(11) kept if RHF determines orbitals,
C         because AUTOCC may change orbital occupation
C         in which case we need to update orb.rot. info /hjaaj.apr2000
      FLAGSV(11) = .FALSE.
C     if MCSCF then orb.rot. info written to file is OK /hjaaj,apr2000
C     No optimal orbital trial vectors nor orbital absorption
      FLAG(41) = .FALSE.
      FLAG(42) = .FALSE.
      FLAG(51) = .FALSE.
      FLAG(52) = .FALSE.
      FLAG(53) = .FALSE.
C     No active-active rotations
      FLAG(23) = .FALSE.
C
C     Do some solvent initializations
C
Cbm initialized PCMSV
      PCMSV   = PCM
c     WRITE(LUPRI,*)'INIZIALIZZO PCMSV',PCMSV
c     IF (PCM) WRITE(LUPRI,*) ' CHIAMA SOLVNT'
      IF (PCM) CALL PCMSOLVNT(WORK,LWORK)
C
C
C ***************************************************
C        *** set RHF parameters for RHF and MP2 ***
C
      IF (DORHF .OR. DOMP2) THEN
         MCTYPE   = 0
         IF (DOCINO .OR. DOCI .OR. DOMC) THEN
            FLAG(25)=.FALSE.
            IF (DOMP2) THEN
               FLAG(16)= .FALSE.
               PCMSV   = PCM
               PCM     = .FALSE.
            END IF
         END IF
C        ... do not write interface file if RHF or MP2 not final calc.
C        ... no solvent if RHF followed by MP2 and CI or MC
C            because solvent not implemented in MP2
C
C
C        Define closed shell RHF wave function
C        Reset below for open shell RHF
C
         ISTASV = ISTATE
         NROOSV = NROOTS
         LROOSV = LROOTS
         IROOSV = IROOT(1)
         LSYMSV = LSYM
         ISPISV = ISPIN
         NACTSV = NACTEL
C
         ISTATE = 1
         NROOTS = 1
         LROOTS = 1
         IROOT(1) = 1
         LSYM   = 1
         ISPIN  = 1
         NACTEL = 0
C
C        Correct for orbital input if MCSCF and/or CI
C
         RHFWOP = .TRUE.
         DO 10 ISYM = 1, NSYM
            NFROSV(ISYM) = NFRO(ISYM)
            NISHSV(ISYM) = NISH(ISYM)
            NASHSV(ISYM) = NASH(ISYM)
            IF (NFRO(ISYM) .NE. NFRRHF(ISYM)) THEN
               NFRO(ISYM) = NFRRHF(ISYM)
               RHFWOP = .FALSE.
            END IF
            IF (NISH(ISYM) .NE. NRHF(ISYM)) THEN
               NISH(ISYM) = NRHF(ISYM)
               RHFWOP = .FALSE.
            END IF
            IF (ISYM .EQ. IOPRHF) THEN
C           ... this is one open shell RHF; reset wave function def.
               NASH(IOPRHF) = 1
               LSYM   = IOPRHF
               ISPIN  = 2
               NACTEL = 1
            ELSE
               NASH(ISYM) = 0
            END IF
            IF (NASH(ISYM) .NE. NASHSV(ISYM)) RHFWOP = .FALSE.
   10    CONTINUE
         IF (HSROHF) THEN
            LSYM = 1
            DO ISYM = 2, NSYM
               IF (MOD(NROHF(ISYM),2).EQ.1) LSYM = MULD2H(LSYM,ISYM)
            END DO
            NASHT = ISUM(NSYM,NROHF,1)
            ISPIN = NASHT+1
            NACTEL = NASHT
            CALL ICOPY(8,NROHF,1,NASH,1)
         END IF
         IF (LSYMSV .NE. LSYM  ) RHFWOP = .FALSE.
         IF (ISPISV .NE. ISPIN ) RHFWOP = .FALSE.
         IF (NACTSV .NE. NACTEL) RHFWOP = .FALSE.
C
C        call sirset to set correct RHF orbital etc. information
C        if necessary
C
         IF (.NOT.RHFWOP) THEN
            OLDWOP = .FALSE.
            JWOPSY = 1
            CALL SIRSET(WORK,LWORK,OLDWOP)
            IAVERR = 0
            CALL AVECHK(WORK,LWORK,IAVERR)
            IF (IAVERR .NE. 0) CALL QUIT(
     &         'SIRCTL RHF error: inconsistency in sup.sym. averaging')
            RHFWOP = .TRUE.
         END IF
      END IF
C
C ***************************************************
C        *** RHF ***
C
      IF (DORHF) THEN
         NUMRUN = NUMRUN + 1
         IF (IPRRHF .GE. 0) THEN
            IPRSIR = IPRRHF
            IPRI6  = IPRRHF
            IPRI4  = IPRRHF
         ELSE
            IPRRHF = IPRSIR
C           ... use global print level in DIIS if IPRRHF not specified
         END IF
         FLAG(2)  = .FALSE.
         FLAG(4)  = .FALSE.
         FLAG(9)  = .FALSE.
C        Fock approximation to orbital Hessian diagonal
         FLAG(12) = .FALSE.
C        Integral transformations not needed
         FLAG(14) = .TRUE.
         FLAG(34) = .TRUE.
C        Use always NEO for RHF (920221-hjaaj)
C        unless .NR ALWAYS (FLAG(39)) is set, e.g. for core relax calc.
         IF (.NOT.FLAG(39)) FLAG(38) = .TRUE.
C
C
         CALL GETTIM(TSTART,WSTART)
C
C        Call FCMO to perform MAXFCK Roothaan Fock iterations
C
         IF (MAXFCK .GT. 0) THEN
            FLAG(5) = .FALSE.
C              we don't call SIROUT when no 2:nd order optimization
            JRDMO = 9
            KCMO  = 1
            KFC   = KCMO + NCMOT
            KWRK1 = KFC  + NNBAST
C           (FC needs NNBAST because it is also used for SAO)
            LWRK1 = LWORK - KWRK1
            CALL READMO(WORK(KCMO),JRDMO,WORK(KWRK1),LWRK1)
            CALL FCMO(MAXFCK,WORK(KCMO),WORK(KFC),
     &                WORK(KWRK1),LWRK1)
C           CALL FCMO(MAXFCK,CMO,FC,SCRA,LSCRA)
            CALL NEWORB('FCMO    ',WORK(KCMO),.TRUE.)
         END IF
C
C        Call SRDIIS to perform MXDIIS DIIS iterations
C
         IF (MXDIIS .GT. 0) THEN
            FLAG(5) = FLAGSV(5)
            JRDMO = 9
            KCMO  = 1
            KWRK1 = KCMO  + NCMOT
            KFREE = 1
            LFREE = LWORK - KWRK1
            CALL READMO(WORK(KCMO),JRDMO,WORK(KWRK1),LFREE)
            CALL SRDIIS(ICONV,WORK(KCMO),WORK(KWRK1),KFREE,LFREE)
C
C     We overwrite in case of symmetry breaking
C
            IF (AUTOCC) THEN
               CALL NEWORB('FOCKDIIS',WORK(KCMO),.TRUE.)
            ELSE
               CALL NEWORB('FOCKDIIS',WORK(KCMO),.FALSE.)
            END IF
C           ... REWIT1 false: do not destroy any GEOWALK information
C           QCSCF currently needed for solvent and for writing SIRIFC
C           QCSCF for writing SIRIFC no longer needed, except for SOPPA 
C           calculations, K.Ruud May-97
C           QCSCF no longer needed for solvent calculations, kr Feb-99
C           QCSCF no longer needed !!! hjaaj may 2000
C
         END IF
         FLAG(11) = .FALSE.
C     ... flag(11) definitely not needed after DIIS /hjaaj.apr2000
C
C        Call SIROPT to perform QC-SCF iterations
C
         ICONVI = ICONV 
         IF (FLAG(21) .AND. ICONV .LE. 0) THEN
            FLAG( 2) = .TRUE.
            FLAG( 5) = FLAGSV(5)
            FLAG(15) = RHFCAN
            MXMASV = MAXMAC
            MXJTSV = MAXJT
            MXMISV = MAXMIC
            MAXMAC = MXHFMA
            MAXJT  = MXHFMI
            MAXMIC = MIN(1000,MAXMAC*MAXJT)
            THRGRD = THRRHF
C
            thrldg = nwopt
            THRLDG = 1.0d2*thrldg*THRLDP
            THRLDG = SQRT(THRLDG)
            IF (THRGRD .LT. THRLDG) THEN
               WRITE(LUPRI,'(//A,2(/A,D12.2))')
     *         ' SIRCTL FATAL ERROR: SCF threshold below lin.dep. limit'
     *        ,' SCF convergence threshold:',THRGRD
     *        ,' linear dependency limit  :',THRLDG
               CALL QUIT(
     &            'SIRCTL error: SCF threshold below lin.dep. limit')
            END IF
            IF (SUPSYM) THEN
C              ... orbital super symmetries may have been reordered
               WRINDX = .TRUE.
               OLDWOP = .FALSE.
               CALL SIRSET(WORK,LWORK,OLDWOP)
               IAVERR = 0
               CALL AVECHK(WORK,LWORK,IAVERR)
               IF (IAVERR .NE. 0) CALL QUIT(
     &         'SIRCTL QCHF error: inconsistency in sup.sym. averaging')
         END IF
C
            CALL SIROPT(WORK,LWORK,ICONV)
            FLAG(15) = FLAGSV(15)
            MAXMAC = MXMASV
            MAXMIC = MXMISV
            MAXJT  = MXJTSV
         END IF
         CALL GETTIM(TEND,WEND)
         TUSED = TEND - TSTART
         WUSED = WEND - WSTART
         IF (IPRSTAT .GT. 0) WRITE (LUSTAT,'(/A,2F12.3)')
     &      ' >>>> CPU and wall time for SCF :',TUSED,WUSED
         IF (IPRSIR .GE. 0) WRITE (LUPRI,'(/A,2F12.3)')
     &      ' >>>> CPU and wall time for SCF :',TUSED,WUSED
         CALL FLSHFO(LUW4)
         CALL FLSHFO(LUPRI)
         CALL FLSHFO(LUSTAT)

         IPRSIR = MPRSIR
         IPRI6  = MPRI6
         IPRI4  = MPRI4
         FLAG(34)   = FLAGSV(34)
         FLAG(14)   = .FALSE.
         FLAGSV(14) = .FALSE.
      END IF
C
C ***************************************************
C        *** MP2 ***
C
      IF (DOMP2) THEN
         IF (DFT_SPINDNS .OR. DFT_LOCALSPIN) THEN
            WRITE(LUPRI,'(//A)') 'WARNING - DFT spin density is not '//
     &      'implemented in MP2 module and is ignored.'
         END IF
         DFT_SPINDNS_SAVE = DFT_SPINDNS
         DFT_SPINDNS = .FALSE.
         DFT_LOCALSPIN_SAVE = DFT_LOCALSPIN
         DFT_LOCALSPIN = .FALSE.
         NUMRUN = NUMRUN + 1
         IF (DORHF .AND. ICONV .LE. 0) THEN
            WRITE (LUPRI,7020)
            WRITE (LUERR ,7020)
            GO TO 8000
         END IF
 7020    FORMAT(/'ERROR: Sirius MP2 calculation aborted because',
     *      ' preceding Hartree-Fock not converged.')
C
         FLAG( 9) = .TRUE.
         FLAG( 2) = .FALSE.
         FLAG( 4) = .FALSE.
         FLAG(14) = .FALSE.
C
         FLAG(21) = .FALSE.
C
         CALL GETTIM(TSTART,WSTART)
         CALL MP2CTL(WORK,LWORK)
         CALL GETTIM(TEND,WEND)
         TUSED = TEND - TSTART
         WUSED = WEND - WSTART
         IF (IPRSTAT .GT. 0) WRITE (LUSTAT,'(/A,2F12.3)')
     &      ' >>>> CPU and wall time for MP2CTL :',TUSED,WUSED
         IF (IPRSIR .GE. 0) WRITE (LUPRI,'(/A,2F12.3)')
     &      ' >>>> CPU and wall time for MP2 :',TUSED,WUSED
         CALL FLSHFO(LUW4)
         CALL FLSHFO(LUPRI)
         CALL FLSHFO(LUSTAT)
         FLAG(14)   = .FALSE.
         FLAGSV(14) = .FALSE.
         DFT_SPINDNS = DFT_SPINDNS_SAVE
         DFT_LOCALSPIN = DFT_LOCALSPIN_SAVE
      END IF
C
C ***************************************************
C        *** reset to CI/MCSCF flags ***
C
      IF ( (DORHF .OR. DOMP2) .AND.
     &     (DOMC .OR. DOCI .OR. DOCINO .OR. FCVORB) ) THEN
         HSROHF = .FALSE.
C        ... otherwise error in some MCSCF routines
C            where HSROHF specifications erroneously will be used
         FLAG(25)   = FLAGSV(25)
         DO 25 ISYM = 1, NSYM
            NFRO(ISYM) = NFROSV(ISYM)
            NISH(ISYM) = NISHSV(ISYM)
            NASH(ISYM) = NASHSV(ISYM)
   25    CONTINUE
         ISTATE = ISTASV
         NROOTS = NROOSV
         LROOTS = LROOSV
         IROOT(1) = IROOSV
         NACTEL = NACTSV
         ISPIN  = ISPISV
         LSYM   = LSYMSV
C
C        Correct orbital information if RHF and/or MP2 was run before.
C
         IF (SUPSYM) THEN
C           ... orbital super symmetries may have been reordered
            WRINDX = .TRUE.
            OLDWOP = .FALSE.
         ELSE
            OLDWOP = .TRUE.
         END IF
         JWOPSY = 1
         MCTYPE = MCTPSV
         CALL SIRSET(WORK,LWORK,OLDWOP)
         IAVERR = 0
         CALL AVECHK(WORK,LWORK,IAVERR)
         IF (IAVERR .NE. 0) CALL QUIT(
     &      'SIRCTL CI/MC error: inconsistency in sup.sym. averaging')
C
      END IF
C
C ***************************************************
C        *** FCVORB module ***
C        (6-May-1994 hjaaj)
C        (moved from DORHF so that NISH is reset to MCSCF/CI values)
C
C        (2-Oct-1986)
C        With canonical Fock orbitals, it is very likely that
C        any diffuse orbitals will be in the active set, which
C        then must be rotated out afterwards in the MC optimization.
C        Therefore, transform virtual orbitals so they diagonalize
C        the one-electron Hamiltonian, this prevents these unwanted
C        diffuse orbitals from being among the active orbitals.
C
C        (28-Aug-1995 hjaaj)
C        Generalized to use modified FC
C        (following C.W. Bauschlicher, JCP 72 (1980) 880 )
C
      IF (FCVORB) THEN
         JRDMO = 9
         KCMO  = 1
         KWRK1 = KCMO  + NCMOT
         LWRK1 = LWORK - KWRK1
         CALL READMO(WORK(KCMO),JRDMO,WORK(KWRK1),LWRK1)
         CALL FCVIRT(WORK(KCMO),WORK(KWRK1),LWRK1)
         CALL NEWORB('FCVORB  ',WORK(KCMO),.TRUE.)
         FLAG(14) = .FALSE.
         FLAGSV(14) = .FALSE.
         IF (SUPSYM) THEN
C           ... orbital super symmetries may have been reordered
            WRINDX = .TRUE.
            OLDWOP = .FALSE.
            CALL SIRSET(WORK,LWORK,OLDWOP)
            IAVERR = 0
            CALL AVECHK(WORK,LWORK,IAVERR)
            IF (IAVERR .NE. 0) CALL QUIT(
     &      'SIRCTL FCVORB error: inconsistency in sup.sym. averaging')
         END IF
      END IF
C
C ***************************************************
C        *** CI ***
C
      IF (DOCI .OR. DOCINO) THEN
         IF (DFT_SPINDNS .OR. DFT_LOCALSPIN) THEN
            WRITE(LUPRI,'(//A)') 'WARNING - DFT spin density is not '//
     &      'implemented in CI module and is ignored.'
         END IF
         DFT_SPINDNS_SAVE = DFT_SPINDNS
         DFT_SPINDNS = .FALSE.
         DFT_LOCALSPIN_SAVE = DFT_LOCALSPIN
         DFT_LOCALSPIN = .FALSE.
         NUMRUN = NUMRUN + 1
         DOHFSRDFT = .FALSE.
C        ... no more HFSRDFT
         FLAG( 4) = .TRUE.
         FLAG( 2) = .FALSE.
         FLAG(21) = .FALSE.
         FLAG( 9) = .FALSE.
         FLAG(15) = (ICICNO .GT. 0)
C        ... if CI it is again OK to call sirout.
         FLAG( 5) = FLAGSV( 5)
         MCTYPE = MCTPSV
         IF (DOMC) FLAG(25) = .FALSE.
C        ... do not write interface file if followed by MCSCF
         thrldg = nconf
         THRLDG = 1.0d2*thrldg*THRLDP
         THRLDG = SQRT(THRLDG)
         IF (THRCI .LT. THRLDG) THEN
            WRITE(LUPRI,'(//A,2(/A,D12.2))')
     *      ' SIRCTL FATAL ERROR: CI threshold below lin.dep. limit',
     *      ' CI convergence threshold:',THRCI,
     *      ' linear dependency limit :',THRLDG
            CALL QUIT('SIRCTL error: CI threshold below lin.dep. limit')
         END IF
C
         CALL GETTIM(TSTART,WSTART)
         CALL CIMAIN(WORK,LWORK,ICONV)
         CALL GETTIM(TEND,WEND)
         TUSED = TEND - TSTART
         WUSED = WEND - WSTART
         IF (IPRSTAT .GT. 0) WRITE (LUSTAT,'(/A,2F12.3)')
     &      ' >>>> CPU and wall time for CIMAIN :',TUSED,WUSED
         IF (IPRSIR .GE. 0) WRITE (LUPRI,'(/A,2F12.3)')
     &      ' >>>> CPU and wall time for CI :',TUSED,WUSED
         CALL FLSHFO(LUW4)
         CALL FLSHFO(LUPRI)
         CALL FLSHFO(LUSTAT)

         IF (DOCINO) THEN
            FLAG(14) = .FALSE.
            FLAGSV(14) = .FALSE.
            IF (SUPSYM) THEN
C              ... orbital super symmetries may have been reordered
               WRINDX = .TRUE.
               OLDWOP = .FALSE.
               CALL SIRSET(WORK,LWORK,OLDWOP)
               IAVERR = 0
               CALL AVECHK(WORK,LWORK,IAVERR)
               IF (IAVERR .NE. 0) CALL QUIT(
     &         'SIRCTL CINO error: inconsistency in sup.sym. averaging')
            END IF
         END IF
         DFT_SPINDNS = DFT_SPINDNS_SAVE
         DFT_LOCALSPIN = DFT_LOCALSPIN_SAVE
      END IF
C
C ***************************************************
C        *** MCSCF ***
C
      IF (DOMC) THEN
         DOHFSRDFT = .FALSE.
         DOCISRDFT = .FALSE.
C        ... no more HFSRDFT nor CISRDFT
         NUMRUN = NUMRUN + 1
C
C Recover MCSCF flags and settings
         DO 3 IFLAG = 1, NFLAG
            FLAG(IFLAG) = FLAGSV(IFLAG)
    3    CONTINUE
         PCM    = PCMSV
Cbm
c        WRITE(LUPRI,*)'rimetto PCM = PCMSV',PCM
         MCTYPE = MCTPSV
C
         FLAG( 2) = .TRUE.
         FLAG(21) = .FALSE.
         FLAG( 4) = .FALSE.
         FLAG( 9) = .FALSE.
         IF (IMCCNO .GT. 0) THEN
            FLAG(15) = .TRUE.
         ELSE IF (IMCCNO .EQ. 0) THEN
            FLAG(15) = .TRUE.
C           IF (MCTYPE .EQ. 2) FLAG(48) = .TRUE.
C           HJMAERKE 921216: use FOCKONLY (FLAG(48))
C              for RAS ?????????
         END IF
C
         IF (DOCINO) THEN
            ICI0 = 1
         ELSE IF (DOCI) THEN
            ICI0 = 4
         END IF
         THRGRD = THRMC
C
C        thrldg = 100*nvar*thrldp
C        ^- gave negative number with both gfortran and ifort
C           when nvar = 37 000 000. Thus the new code. /aug08-hjaaj
         thrldg = nvar
         THRLDG = 1.0d2*thrldg*THRLDP
         THRLDG = SQRT(THRLDG)
         IF (THRGRD .LT. THRLDG) THEN
            WRITE(LUPRI,'(//A,2(/A,D12.2))')
     *      ' SIRCTL FATAL ERROR: MC threshold below lin.dep. limit',
     *      ' MC convergence threshold:',THRGRD,
     *      ' linear dependency limit :',THRLDG
            CALL QUIT('SIRCTL error: MC threshold below lin.dep. limit')
         END IF
C
         CALL GETTIM(TSTART,WSTART)
         CALL SIROPT(WORK,LWORK,ICONV)
         CALL GETTIM(TEND,WEND)
         TUSED = TEND - TSTART
         WUSED = WEND - WSTART
         IF (IPRSTAT .GT. 0) WRITE (LUSTAT,'(/A,2F12.3)')
     &      ' >>>> CPU and wall time for MCSCF :',TUSED,WUSED
         IF (IPRSIR .GE. 0) WRITE (LUPRI,'(/A,2F12.3)')
     &      ' >>>> CPU and wall time for MCSCF :',TUSED,WUSED
         CALL FLSHFO(LUW4)
         CALL FLSHFO(LUPRI)
         CALL FLSHFO(LUSTAT)
      END IF
C
      IF (ICONV .EQ. -1) THEN
C        *** set ICONV to converged ***
         ICONV = 1
      END IF
      IF (ICONV .EQ. 0) THEN
         NWARN = NWARN + 1
         WRITE (LUERR,7070)
         WRITE (LUPRI,7070)
         IF (LUW4.NE.LUPRI) WRITE (LUW4,7070)
      END IF
 7070 FORMAT(//' WARNING: wave function not converged')
C
C *** Output sections
C
C *** ORBITAL OUTPUT
      IF (LUPROP .GE. 0) CALL GPCLOSE(LUPROP,'KEEP')
      IF (FLAG(5) .AND. NUMRUN .GT. 0) THEN
         CALL SIROUT(ICONV,WORK,LWORK)
         CALL FLSHFO(LUW4)
         IF (LUPRI.NE.LUW4) CALL FLSHFO(LUPRI)
         IF (LUSTAT.NE.LUW4 .AND. LUSTAT.NE.LUPRI) CALL FLSHFO(LUSTAT)
      END IF
      IF (LUINF .GT.0) CALL GPCLOSE(LUINF,'DELETE')
C
C ***************************************************
C     *** NEVPT ***
C
      IF (DONEVPT) THEN
       IF (ICONV .EQ. 0) THEN
         WRITE (LUPRI,'(//A)')
     &      ' WARNING: NEVPT2 skipped because MCSCF not converged.'
       ELSE
         JRDMO = 9
         KCMO  = 1
         KWRK1 = KCMO  + NCMOT
         LWRK1 = LWORK - KWRK1
         CALL READMO(WORK(KCMO),JRDMO,WORK(KWRK1),LWRK1)
         USEDRC = .TRUE.
         IF (NEWTRA) THEN
C        ... koopro4 uses (ia/ (inactive-secondary) Mulliken distributions;
C            they are only calculated for level 5 in "NEWTRA"
C            transformation  /hjaaj June 09
            CALL TRACTL(5,WORK(KCMO),WORK(KWRK1),LWRK1)
         ELSE
            CALL TRACTL(4,WORK(KCMO),WORK(KWRK1),LWRK1)
         END IF

         CALL GETTIM(TSTART,WSTART)
         CALL KOOPRO4(WORK,LWORK)
         CALL GETTIM(TEND,WEND)
         TUSED = TEND - TSTART
         WUSED = WEND - WSTART

         IF (IPRSTAT .GT. 0) WRITE (LUSTAT,'(/A,2F12.3)')
     &      ' >>>> CPU and wall time for NEVPT2 :',TUSED,WUSED
         IF (IPRSIR .GE. 0) WRITE (LUPRI,'(/A,2F12.3)')
     &      ' >>>> CPU and wall time for NEVPT2 :',TUSED,WUSED
         CALL FLSHFO(LUW4)
         CALL FLSHFO(LUPRI)
         CALL FLSHFO(LUSTAT)

       END IF
      END IF
C
C ***************************************************
C     *** LUCITA ***
C
      IF (DOLUCITA) THEN
       IF (ICONV .EQ. 0) THEN
         WRITE (LUPRI,'(//A)')
     &      ' WARNING: LUCITA skipped because HF/MCSCF not converged.'
       ELSE
         kfirst_lucita = 1
         JRDMO  = 9
         KFRSAV = 1
         KFREE  = KFRSAV
         LFREE  = LWORK
         CALL MEMGET('REAL',KCMO,NCMOT,WORK,KFREE,LFREE)
         CALL READMO(WORK(KCMO),JRDMO,WORK(KFREE),LFREE)

         CALL SET_INFORB_LUCITA('LUCITA') ! define inforb.h with LUCITA
                                          ! orbital occupations for TRACTL
         JTRALVL = 0 ! all active orbitals for all four indices
         USEDRC  = .FALSE.

!        integral transformation
         CALL TRACTL(JTRALVL,WORK(KCMO),WORK(KFREE),LFREE)

!        prepare (aa|aa) 1e-/2e- integral file to be read within LUCITA
         CALL MAKE_LUCITA_INTEGRALS(WORK(KCMO),WORK(KFREE),LFREE)

         CALL MEMREL('lucita-tra.done',work,KFRSAV,KFRSAV,KFREE,LFREE)

!        perform large-scale (GAS)CI calculation with LUCITA
         CALL GETTIM(TSTART,WSTART)
         CALL LUCITA(WORK,LFREE)
         CALL GETTIM(TEND,WEND)
         TUSED = TEND - TSTART
         WUSED = WEND - WSTART

         IF (IPRSTAT .GT. 0) WRITE (LUSTAT,'(/A,2F12.3)')
     &      ' >>>> CPU and wall time for LUCITA CI :',TUSED,WUSED
         IF (IPRSIR .GE. 0) WRITE (LUPRI,'(/A,2F12.3)')
     &      ' >>>> CPU and wall time for LUCITA CI :',TUSED,WUSED
         CALL FLSHFO(LUW4)
         CALL FLSHFO(LUPRI)
         CALL FLSHFO(LUSTAT)

         CALL SET_INFORB_LUCITA('RESET ')
       END IF
      END IF
C
C *** POPULATION ANALYSIS AND PROPERTIES
C
      IF (FLAG(6) .AND. NUMRUN .GT. 0) THEN
         KDV   = 1
         KWRK1 = KDV + NNASHX
         LWRK1 = LWORK - KWRK1
         CALL SIRPOP('FINAL',WORK(KDV),WORK(KWRK1),LWRK1)
         CALL FLSHFO(LUW4)
         IF (LUPRI.NE.LUW4) CALL FLSHFO(LUPRI)
         IF (LUSTAT.NE.LUW4 .AND. LUSTAT.NE.LUPRI) CALL FLSHFO(LUSTAT)
      END IF
C
 8000 CALL GETTIM(TEND,WEND)
      TIMSIR = TEND - TIMSIR
      WALSIR = WEND - WALSIR
      WRITE (LUPRI,7100)  TIMSIR, WALSIR
      CALL TSTAMP(' ',LUPRI)
      IF (LUW4 .NE. LUPRI) THEN
         WRITE (LUW4,7100)  TIMSIR, WALSIR
         CALL TSTAMP(' ',LUW4)
      END IF
      IF (IPRSTAT .GT. 0)  THEN
         WRITE (LUSTAT,7100) TIMSIR, WALSIR
         CALL TSTAMP(' ',LUSTAT)
      END IF
 7100 FORMAT(//
     &   /' >>>> Total CPU  time used in SIRIUS :',F10.2,' seconds',
     &   /' >>>> Total wall time used in SIRIUS :',F10.2,' seconds')
C
      IF (NINFO .GT. 0) THEN
         WRITE (LUPRI,7210) NINFO
         IF (LUW4.NE.LUPRI) WRITE (LUW4,7210) NINFO
         WRITE (LUSTAT,7210) NINFO
         WRITE (LUERR,7210) NINFO
      END IF
      IF (NWARN .GT. 0) THEN
         WRITE (LUPRI,7200) NWARN
         IF (LUW4.NE.LUPRI) WRITE (LUW4,7200) NWARN
         WRITE (LUSTAT,7200) NWARN
         WRITE (LUERR,7200) NWARN
      END IF
 7200 FORMAT(/' NOTE:',I5,' warnings have been issued.',
     &       /' Check output, result, and error files for "WARNING".')
 7210 FORMAT(/' NOTE:',I5,' informational messages have been issued.',
     &       /' Check output, result, and error files for "INFO".')
C
      IF (FLAG(16)) CALL GPCLOSE(LUSOL,'KEEP')
      CALL GPCLOSE(LUIT1,'KEEP')
      IF (LUH2AC .GT. 0) CALL GPCLOSE(LUH2AC,'DELETE')
      IF (LUW4 .NE. LUPRI) CALL GPCLOSE(LUW4,'KEEP')
C
      IF (PCM) CALL GPCLOSE(LUPCMD,'KEEP')
      CALL QEXIT('SIRCTL')
      RETURN
C *** end of SIRCTL.
      END
C
C  /* Deck sbdorb */
      BLOCK DATA SBDORB
C     Initialize MULD2H in common block /INFORB/
#include <implicit.h>
#include <inforb.h>
C
C     MULTIPLICATION TABLE FOR SYMMETRIES (MULD2H is in /INFORB/)
C
      DATA MULD2H/1,2,3,4,5,6,7,8,
     *            2,1,4,3,6,5,8,7,
     *            3,4,1,2,7,8,5,6,
     *            4,3,2,1,8,7,6,5,
     *            5,6,7,8,1,2,3,4,
     *            6,5,8,7,2,1,4,3,
     *            7,8,5,6,3,4,1,2,
     *            8,7,6,5,4,3,2,1/
      END
C  /* Deck sbdgetd */
      BLOCK DATA SBDGETD
C 900302-hjaaj
C Define default values for CBGETD
C To be safe, include EXTERNAL SBDGETD in main program.
#include <implicit.h>
#include <cbgetdis.h>
      DATA DISTYP,IADINT,IADH2,IADH2X,IADH2D /1,4*-1/
      END
C  /* Deck sbdtra */
      BLOCK DATA SBDTRA
C 951130-hjaaj
C Initialize variables in /INFTRA/
C Initialization moved to here because NEWTRA may be set
C both in dalton and sirius input modules.
C To be safe, include EXTERNAL SBDTRA in sirctl
C
#include <implicit.h>
#include <inftra.h>
      DATA NEWTRA /.FALSE./
      END
C  /* Deck sbddim */
      BLOCK DATA SBDDIM
C 010321 hjaaj
C Initialize variables in /INFDIM/
C
#include <implicit.h>
#include <infdim.h>
      DATA MWORK /8 000 000/
C
C     Re MWORK:
C     default 8 mill. words for number of simultaneous AO distributions
C     in "old" integral transformation - defined for SORTA. The size
C     should not be bigger than resident memory, to avoid paging.
C     If MWORK too big, then SORTA may define more simult. dist. than
C     can be treated at later calls to TRACTL !!! /Mar 2001 hjaaj
C
      END
