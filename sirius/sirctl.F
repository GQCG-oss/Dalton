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
C   INFDIM : MWORK
C   INFTRA : USEDRC
C   INFPRI : IPRERR, NWARN
C   SCBRHF : RHFCAN, MXHFMA, MXHFMI, NFRRHF(8), IOPRHF, MAXFCK, MXDIIS, THRRHF
C   GNRINF : WRINDX
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <infvar.h>
#include <inftap.h>
#include <infdim.h>
#include <inftra.h>
#include <infpri.h>
#include <scbrhf.h>
#include <gnrinf.h>
C
      LOGICAL   FLAGSV(NFLAG), OLDWOP, RHFWOP
      DIMENSION NASHSV(8), NISHSV(8), NFROSV(8)
C
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
      IF (IPRERR .GT. 0) THEN
         WRITE (LUERR,'(//A,2F10.2)')
     *      ' >>>> CPU and wall time for SIRINP :',TINP,WINP
      END IF
      CALL FLSHFO(LUW4)
      IF (LUPRI.NE.LUW4) CALL FLSHFO(LUPRI)
      IF (LUERR.NE.LUW4 .AND. LUERR.NE.LUPRI) CALL FLSHFO(LUERR)
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
C ***************************************************
C        *** set RHF parameters for RHF and MP2 ***
C
      IF (DORHF .OR. DOMP2) THEN
         MCTYPE   = 0
         IF (DOCINO .OR. DOCI .OR. DOMC) THEN
            FLAG(25)=.FALSE.
            IF (DOMP2) FLAG(16)=.FALSE. 
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
         IF (FLAG(21) .AND. ICONV .LE. 0) THEN
            FLAG( 2) = .TRUE.
            FLAG( 5) = FLAGSV(5)
            FLAG(15) = RHFCAN
            MXMASV = MAXMAC
            MXJTSV = MAXJT
            MXMISV = MAXMIC
            MAXMAC = MXHFMA
            MAXJT  = MXHFMI
            MAXMIC = MIN(200,MAXMAC*MAXJT)
            THRGRD = THRRHF
C
            THRLDG = 100*NWOPT*THRLDP
            THRLDG = SQRT(THRLDG)
            IF (THRGRD .LT. THRLDG) THEN
               WRITE(LUPRI,'(//A,2(/A,D12.2))')
     *         ' SIRCTL error: SCF threshold below lin.dep. limit',
     *         ' SCF convergence threshold:',THRGRD,
     *         ' lin.dep. limit           :',THRLDG
               CALL QUIT(
     &            'SIRCTL error: SCF threshold below lin.dep. limit')
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
         IF (IPRERR .GT. 0) WRITE (LUERR,'(//A,2F10.2)')
     *      ' >>>> CPU and wall time for SCF :',TUSED,WUSED
         CALL FLSHFO(LUW4)
         IF (LUPRI.NE.LUW4) CALL FLSHFO(LUPRI)
         IF (LUERR.NE.LUW4 .AND. LUERR.NE.LUPRI) CALL FLSHFO(LUERR)
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
         NUMRUN = NUMRUN + 1
         IF (DORHF .AND. ICONV .LE. 0) THEN
            WRITE (LUPRI,7020)
            WRITE (LUERR ,7020)
            GO TO 8000
         END IF
 7020    FORMAT(//' Sirius MP2 calculation aborted because',
     *      ' preceeding Hartree-Fock not converged.')
C
         FLAG( 9) = .TRUE.
         FLAG( 2) = .FALSE.
         FLAG( 4) = .FALSE.
         FLAG(14) = .FALSE.
C
C        Transform to canonical orbitals, unless we know that
C        the orbitals already are canonical orbitals.
         IF (.NOT.DORHF .OR. .NOT.RHFCAN) THEN
            IPRCAN = 0
            KFRCAN = 1
            LFRCAN = LWORK
            CALL MEMGET('REAL',KCREF,NCONF,WORK,KFRCAN,LFRCAN)
            CALL MEMGET('REAL',KCMO ,NCMOT,WORK,KFRCAN,LFRCAN)
            CALL MEMGET('REAL',KOCC ,NORBT,WORK,KFRCAN,LFRCAN)
            CALL SIRCAN(IPRCAN,WORK(KCREF),WORK(KCMO),WORK(KOCC),
     *                  WORK,KFRCAN,LFRCAN)
            CALL MEMREL('SIRCTL.MP2.SIRCAN',WORK,1,1, KFRCAN,LFRCAN)
         END IF
         FLAG(21) = .FALSE.
C
         CALL GETTIM(TSTART,WSTART)
         CALL MP2CTL(WORK,LWORK)
         CALL GETTIM(TEND,WEND)
         TUSED = TEND - TSTART
         WUSED = WEND - WSTART
         IF (IPRERR .GT. 0) WRITE (LUERR,'(//A,2F10.2)')
     *      ' >>>> CPU and wall time for MP2CTL :',TUSED,WUSED
         CALL FLSHFO(LUW4)
         IF (LUPRI.NE.LUW4) CALL FLSHFO(LUPRI)
         IF (LUERR.NE.LUW4 .AND. LUERR.NE.LUPRI) CALL FLSHFO(LUERR)
         FLAG(14)   = .FALSE.
         FLAGSV(14) = .FALSE.
      END IF
C
C ***************************************************
C        *** reset to CI/MCSCF flags ***
C
      IF ( (DORHF .OR. DOMP2) .AND.
     &     (DOMC .OR. DOCI .OR. DOCINO .OR. FCVORB) ) THEN
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
         OLDWOP = .TRUE.
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
      END IF
C
C ***************************************************
C        *** CI ***
C
      IF (DOCI .OR. DOCINO) THEN
         NUMRUN = NUMRUN + 1
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
         THRLDG = 100*NCONF*THRLDP
         THRLDG = SQRT(THRLDG)
         IF (THRCI .LT. THRLDG) THEN
            WRITE(LUPRI,'(//A,2(/A,D12.2))')
     *      ' SIRCTL error: CI threshold below lin.dep. limit',
     *      ' CI convergence threshold:',THRCI,
     *      ' lin.dep. limit          :',THRLDG
            CALL QUIT('SIRCTL error: CI threshold below lin.dep. limit')
         END IF
C
         CALL GETTIM(TSTART,WSTART)
         CALL CIMAIN(WORK,LWORK,ICONV)
         CALL GETTIM(TEND,WEND)
         TUSED = TEND - TSTART
         WUSED = WEND - WSTART
         IF (IPRERR .GT. 0) WRITE (LUERR,'(//A,2F10.2)')
     *      ' >>>> CPU and wall time for CIMAIN :',TUSED,WUSED
         CALL FLSHFO(LUW4)
         IF (LUPRI.NE.LUW4) CALL FLSHFO(LUPRI)
         IF (LUERR.NE.LUW4 .AND. LUERR.NE.LUPRI) CALL FLSHFO(LUERR)
         IF (DOCINO) THEN
            FLAG(14) = .FALSE.
            FLAGSV(14) = .FALSE.
         END IF
      END IF
C
C ***************************************************
C        *** MCSCF ***
C
      IF (DOMC) THEN
         NUMRUN = NUMRUN + 1
C
C Recover MCSCF flags and settings
         DO 3 IFLAG = 1, NFLAG
            FLAG(IFLAG) = FLAGSV(IFLAG)
    3    CONTINUE
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
         THRLDG = 100*NVAR*THRLDP
         THRLDG = SQRT(THRLDG)
         IF (THRGRD .LT. THRLDG) THEN
            WRITE(LUPRI,'(//A,2(/A,D12.2))')
     *      ' SIRCTL error: MC threshold below lin.dep. limit',
     *      ' MC convergence threshold:',THRGRD,
     *      ' lin.dep. limit          :',THRLDG
            CALL QUIT('SIRCTL error: MC threshold below lin.dep. limit')
         END IF
C
         CALL GETTIM(TSTART,WSTART)
         CALL SIROPT(WORK,LWORK,ICONV)
         CALL GETTIM(TEND,WEND)
         TUSED = TEND - TSTART
         WUSED = WEND - WSTART
         IF (IPRERR .GT. 0) WRITE (LUERR,'(//A,2F10.2)')
     *      ' >>>> CPU and wall time for SIROPT :',TUSED,WUSED
         CALL FLSHFO(LUW4)
         IF (LUPRI.NE.LUW4) CALL FLSHFO(LUPRI)
         IF (LUERR.NE.LUW4 .AND. LUERR.NE.LUPRI) CALL FLSHFO(LUERR)
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
      IF (FLAG(5) .AND. NUMRUN .GT. 0) THEN
         CALL SIROUT(ICONV,WORK,LWORK)
         CALL FLSHFO(LUW4)
         IF (LUPRI.NE.LUW4) CALL FLSHFO(LUPRI)
         IF (LUERR.NE.LUW4 .AND. LUERR.NE.LUPRI) CALL FLSHFO(LUERR)
      END IF
      IF (LUINF .GT.0) CALL GPCLOSE(LUINF,'DELETE')
C
C ***************************************************
C        *** NEVPT ***
C
      IF (DONEVPT) THEN
         JRDMO = 9
         KCMO  = 1
         KWRK1 = KCMO  + NCMOT
         LWRK1 = LWORK - KWRK1
         CALL READMO(WORK(KCMO),JRDMO,WORK(KWRK1),LWRK1)
         USEDRC = .TRUE.
         CALL TRACTL(4,WORK(KCMO),WORK(KWRK1),LWRK1,MWORK)
         CALL KOOPRO4(WORK,LWORK)
      ENDIF
C
C *** POPULATION ANALYSIS AND PROPERTIES
C
      IF (FLAG(6) .AND. NUMRUN .GT. 0) THEN
         CALL SIRPOP(WORK,LWORK)
         CALL FLSHFO(LUW4)
         IF (LUPRI.NE.LUW4) CALL FLSHFO(LUPRI)
         IF (LUERR.NE.LUW4 .AND. LUERR.NE.LUPRI) CALL FLSHFO(LUERR)
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
      IF (IPRERR .GT. 0)  THEN
         WRITE (LUERR,7100) TIMSIR, WALSIR
         CALL TSTAMP(' ',LUERR)
      END IF
 7100 FORMAT(//
     *   /' >>>> Total CPU  time used in SIRIUS :',F10.2,' seconds',
     *   /' >>>> Total wall time used in SIRIUS :',F10.2,' seconds')
C
      IF (NINFO .GT. 0) THEN
         WRITE (LUPRI,7210) NINFO
         IF (LUW4.NE.LUPRI) WRITE (LUW4,7210) NINFO
         WRITE (LUERR,7210) NINFO
      END IF
      IF (NWARN .GT. 0) THEN
         WRITE (LUPRI,7200) NWARN
         IF (LUW4.NE.LUPRI) WRITE (LUW4,7200) NWARN
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
      CALL QEXIT('SIRCTL')
      RETURN
C *** end of SIRCTL.
      END
