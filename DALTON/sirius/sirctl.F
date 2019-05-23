!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!
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

      use lucita_cfg
      use gasci_input_cfg
      use lucita_mcscf_ci_cfg
      use mcscf_or_gasci_2_define_cfg
      use lucita_orbital_spaces, only: set_inforb_lucita
      use lucita_setup, only: setup_lucita_mcci_wrkspc_dimensions
      use lucita_mcscf_ci_interface_procedures
#ifdef ENABLE_EFS
      use efs_interface
#endif
#ifdef HAS_PCMSOLVER
      use pcm_config, only: pcm_configuration, pcm_cfg
#endif
#ifdef ENABLE_XCINT
      use dalton_xcint_interface, only: xcint_is_active, xcint_init
#endif

#include "implicit.h"
C
      DIMENSION WORK(LWORK)
#include "iratdef.h"
#include "dummy.h"
#include "thrldp.h"
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
#include "maxorb.h"
#include "priunit.h"
#include "infinp.h"
#include "inforb.h"
#include "infvar.h"
#include "inftap.h"
#include "inftra.h"
#include "infpri.h"
#include "scbrhf.h"
#include "gnrinf.h"
#include "mxcent.h"
#include "pcmdef.h"
#include "pcm.h"
#include "pcmlog.h"
! dftcom.h: DFT_SPINDNS, DFT_LOCALSPIN
#include "dftcom.h"
#include "exeinf.h"
#include "stex.h"
! for THRQ2
#include "cbisor.h"


C
      LOGICAL   FLAGSV(NFLAG), OLDWOP, RHFWOP, PCM_save
      LOGICAL   DFT_SPINDNS_save, DFT_LOCALSPIN_save, USEDRC_save
      LOGICAL   dftadd_s, addsri_s
      DIMENSION NASHSV(8), NISHSV(8), NFROSV(8)
      real(8), allocatable :: tmp_space(:)
      integer              :: xctype1, xctype2
 
#include "sirbkd.h"
#include "center.h"
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
     &      ' CPU and wall time for SIRINP :',TINP,WINP
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
      MCTYPE_save = MCTYPE
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
Cbm initialized PCM_save
      PCM_save   = PCM
      IF (PCM) THEN
         CALL PCMSOLVNT(WORK,LWORK)
      END IF


#ifdef ENABLE_XCINT
      if (xcint_is_active()) then
        ! initialize XCint interface
        call xcint_init()
      end if
#endif

C
C
C ***************************************************
C        *** set SCF parameters for RHF/DFT and MP2 ***
C
      IF (DOSCF .OR. DOMP2) THEN
         MCTYPE   = 0
         IF (HSROHF) MCTYPE = -1
         IF (DOMP2 ) THEN
            FLAG(25) = .TRUE. ! SIRIFC from RHF needed for MP2 if embedding
         END IF
         IF (DOCINO .OR. DOCI .OR. DOMC .OR. DOLUCITA) THEN
C        ... no continuum models if RHF followed by MP2 and CI or MC
            FLAG(16) = .FALSE.
            PCM = .FALSE.
         END IF
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
         NASHT = 0
         DO ISYM = 1, NSYM
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
C           ... this is one open shell RHF; reset wave function definition
               NASH(IOPRHF) = 1
               LSYM   = IOPRHF
               ISPIN  = 2
               NACTEL = 1
               NASHT  = 1
            ELSE
               NASH(ISYM) = 0
            END IF
            IF (NASH(ISYM) .NE. NASHSV(ISYM)) RHFWOP = .FALSE.
         END DO
         IF (HSROHF) THEN
            LSYM = 1
            DO ISYM = 2, NSYM
               IF (MOD(NROHF(ISYM),2).EQ.1) LSYM = MULD2H(LSYM,ISYM)
            END DO
            NASHT  = ISUM(NSYM,NROHF,1)
            ISPIN  = NASHT + 1
            NACTEL = NASHT
            NASH(1:8) = NROHF(1:8)
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
            CALL AVECHK(IAVERR)
            IF (IAVERR .NE. 0) CALL QUIT(
     &         'SIRCTL RHF error: inconsistency in sup.sym. averaging')
            RHFWOP = .TRUE.
         END IF
      END IF
C
C ***************************************************
C        *** RHF or DFT ***
C
      IF (DOSCF) THEN
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
C        Fock approximation to orbital Hessian diagonal:
         FLAG(12) = .FALSE.
C        Integral transformations not needed:
         FLAG(14) = .TRUE.
         FLAG(34) = .TRUE.
C        Use always NEO for RHF (920221-hjaaj)
C        unless .NR ALWAYS (FLAG(39)) is set, e.g. for core relax calc.
         IF (.NOT.FLAG(39)) FLAG(38) = .TRUE.
C
C        Also high-spin Fock matrices for all open-shell calcs with dft
         HSROHF = HSROHF .OR. (DODFT .AND. NASHT .GE. 1)
C
Ck
         CALL GETTIM(TSTART,WSTART)

#ifdef ENABLE_EFS
         if (use_efs) then
            call mpixbcast (const_generate_basis, 1, 'INTEGER', 0)
            iprint = 0
            call mpixbcast (iprint, 1, 'INTEGER', 0)
            call efs_generate_basis_product(THRQ2) ! gaussian product threshold

            call mpixbcast (const_init_2efock  , 1, 'INTEGER', 0)
            iprint = 0
            call mpixbcast (iprint, 1, 'INTEGER', 0)
            call efs_init_2efock(HFXFAC)
         end if
#endif

C
C        Modifications for open-shell HF-srDFT ...
C        (do not switch off DIIS if DFT_SPINDNS but NASHT.eq.0, i.e. closed-shell SCF)
C
         IF (DOHFSRDFT .and. DFT_SPINDNS .AND. NASHT.GT.0) THEN
         IF (MXDIIS .gt. 0) THEN ! DIIS not implemented for HFSRDFT and open shell SCF
            MXDIIS = 0
            NINFO  = NINFO + 1
            WRITE(LUPRI,'(/A)') '@ INFO: DIIS disabled '//
     &      'for open-shell HF-srDFT because not implemented.'
         END IF
         IF (.NOT.DIRCAL .OR. .NOT.DIRFCK) THEN
            NINFO  = NINFO + 1
            WRITE (LUPRI,'(/A)')
     &      '@ INFO: switched to AO-direct Fock matrices because'
     &      //' open-shell HF-srDFT will not work otherwise.'
            DIRCAL = .TRUE.
            DIRFCK = .TRUE.
         END IF
         END IF
C
C        Call FCMO to perform MAXFCK Roothaan Fock iterations
C
         IF (MAXFCK .GT. 0) THEN
            FLAG(5) = .FALSE.
C           we do not call SIROUT when no 2:nd order optimization
            JRDMO = 9
            KCMO  = 1
            KFC   = KCMO + NCMOT
            KWRK1 = KFC  + NNBAST
C           (FC needs NNBAST because it is also used for SAO)
            LWRK1 = LWORK - KWRK1
            CALL READMO(WORK(KCMO),JRDMO)
            CALL FCMO(MAXFCK,WORK(KCMO),WORK(KFC),
     &                WORK(KWRK1),LWRK1)
C           CALL FCMO(MAXFCK,CMO,FC,SCRA,LSCRA)
            CALL NEWORB('FCMO    ',WORK(KCMO),.TRUE.)
         END IF
C
C        Call SRDIIS to perform MXDIIS DIIS iterations
C
         IF (MXDIIS .gt. 0) THEN
            FLAG(5) = FLAGSV(5)
            JRDMO = 9
            KCMO  = 1
            KWRK1 = KCMO  + NCMOT
            KFREE = 1
            LFREE = LWORK - KWRK1
            CALL READMO(WORK(KCMO),JRDMO)
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
C        ... flag(11) definitely not needed after DIIS /hjaaj.apr2000
C
C        Call SIROPT to perform QC-SCF iterations
C
         ICONVI = ICONV
         IF (FLAG(21) .AND. ICONV .LE. 0) THEN
#ifdef HAS_PCMSOLVER
            if (pcm_cfg%do_pcm) call quit('ERROR: DIIS '//
     &      'failed to converge and QC-SCF is not implemented with '//
     &      'PCMSolver')
#endif
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
     &         ' SIRCTL FATAL ERROR: SCF threshold below lin.dep. limit'
     &        ,' SCF convergence threshold:',THRGRD
     &        ,' linear dependency limit  :',THRLDG
               CALL QUIT(
     &            'SIRCTL error: SCF threshold below lin.dep. limit')
            END IF
            IF (SUPSYM) THEN
C              ... orbital super symmetries may have been reordered
               WRINDX = .TRUE.
               OLDWOP = .FALSE.
               CALL SIRSET(WORK,LWORK,OLDWOP)
               IAVERR = 0
               CALL AVECHK(IAVERR)
               IF (IAVERR .NE. 0) CALL QUIT('SIRCTL error: '//
     &         'inconsistency in sup.sym. averaging before QC-SCF call')
            END IF

! Tell srdft routines to use MCSRDFT version, otherwise qc-hf-srdft will not work for open shells /May2019-hjaaj
            IF (DOHFSRDFT .AND. NASHT.GT.0) DOMCSRDFT = .TRUE.
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
     &      ' CPU and wall time for SCF :',TUSED,WUSED
         IF (IPRSIR .GE. 0) WRITE (LUPRI,'(/A,2F12.3)')
     &      ' CPU and wall time for SCF :',TUSED,WUSED
         CALL FLSHFO(LUW4)
         CALL FLSHFO(LUPRI)
         CALL FLSHFO(LUSTAT)

         IPRSIR = MPRSIR
         IPRI6  = MPRI6
         IPRI4  = MPRI4
         FLAG(25)   = FLAGSV(25)
         FLAG(34)   = FLAGSV(34)
         FLAG(14)   = .FALSE.
         FLAGSV(14) = .FALSE.

!#define AO_DENSITIES_ON_FILE
#ifdef AO_DENSITIES_ON_FILE
         KDV   = 1
         call write_aodens(work(kdv),nisht,nasht,nbast,nnashx,
     &                     ncmot,nsym,nbas,ibas,lupri,hsrohf,"HF")
#undef AO_DENSITIES_ON_FILE
#endif
      END IF
C
C ***************************************************
C        *** SELACT ***
C
      IF (ACTSEL .AND. (.NOT. DOCCSD)) THEN
         NUMRUN = NUMRUN + 1
         IF (DOSCF .AND. ICONV .LE. 0) THEN
            WRITE (LUPRI,7021)
            WRITE (LUERR ,7021)
            GO TO 8000
         END IF
 7021    FORMAT(//' Selection of active orbitals aborted because',
     *      ' preceeding Hartree-Fock not converged.')
C
         CALL CC_SELACT(WORK,LWORK)
      END IF
C
C ***************************************************
C        *** MP2 ***
C
      IF (DOMP2) THEN
         IF (DOCINO .OR. DOCI .OR. DOMC .OR. DOLUCITA) THEN
            FLAG(25) = .FALSE. ! do not write to SIRIFC if MP2 followed by MC or CI
         END IF
         IF (DFT_SPINDNS .OR. DFT_LOCALSPIN) THEN
            NWARN = NWARN + 1
            WRITE(LUPRI,'(//A)') 'WARNING - DFT spin density is not '//
     &      'implemented in MP2 module and is ignored.'
         END IF
         DFT_SPINDNS_save = DFT_SPINDNS
         DFT_SPINDNS = .FALSE.
         DFT_LOCALSPIN_save = DFT_LOCALSPIN
         DFT_LOCALSPIN = .FALSE.
         NUMRUN = NUMRUN + 1
         IF (DOSCF .AND. ICONV .LE. 0) THEN
            WRITE (LUPRI,7020)
            WRITE (LUERR,7020)
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
     &      ' CPU and wall time for MP2CTL :',TUSED,WUSED
         IF (IPRSIR .GE. 0) WRITE (LUPRI,'(/A,2F12.3)')
     &      ' CPU and wall time for MP2 :',TUSED,WUSED
         CALL FLSHFO(LUW4)
         CALL FLSHFO(LUPRI)
         CALL FLSHFO(LUSTAT)
         FLAG(14)   = .FALSE.
         FLAGSV(14) = .FALSE.
         DFT_SPINDNS = DFT_SPINDNS_save
         DFT_LOCALSPIN = DFT_LOCALSPIN_save
      END IF
C
C ***************************************************
C        *** reset to CI/MCSCF flags ***
C
      IF ( (DOSCF .OR. DOMP2) .AND.
     &     (DOMC .OR. DOCI .OR. DOCINO .OR. FCVORB) ) THEN
         HSROHF = .FALSE.
C        ... otherwise error in some MCSCF routines
C            where HSROHF specifications erroneously will be used
         IOPRHF  = -1
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
         MCTYPE = MCTYPE_save
         IF  (DOMC .OR. DOCI .OR. DOCINO .OR. FCVORB) THEN
            ! do not call SIRSET if next step is LUCITA:
            ! SETCI will return NCONF=0 and cause error exit.
            CALL SIRSET(WORK,LWORK,OLDWOP)
            IAVERR = 0
            CALL AVECHK(IAVERR)
            IF (IAVERR .NE. 0) CALL QUIT(
     &      'SIRCTL CI/MC error: inconsistency in sup.sym. averaging')
         END IF
C
      END IF
C
C ***************************************************
C        *** FCVORB module ***
C        (6-May-1994 hjaaj)
C        (moved from DOSCF so that NISH is reset to MCSCF/CI values)
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
         CALL READMO(WORK(KCMO),JRDMO)
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
            CALL AVECHK(IAVERR)
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
            NWARN = NWARN + 1
            WRITE(LUPRI,'(//A)') 'WARNING - DFT spin density is not '//
     &      'implemented in CI module and is ignored.'
         END IF
         DFT_SPINDNS_save = DFT_SPINDNS
         DFT_SPINDNS = .FALSE.
         DFT_LOCALSPIN_save = DFT_LOCALSPIN
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
         MCTYPE = MCTYPE_save
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
     &      ' CPU and wall time for CIMAIN :',TUSED,WUSED
         IF (IPRSIR .GE. 0) WRITE (LUPRI,'(/A,2F12.3)')
     &      ' CPU and wall time for CI :',TUSED,WUSED
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
               CALL AVECHK(IAVERR)
               IF (IAVERR .NE. 0) CALL QUIT(
     &         'SIRCTL CINO error: inconsistency in sup.sym. averaging')
            END IF
         END IF
         DFT_SPINDNS = DFT_SPINDNS_save
         DFT_LOCALSPIN = DFT_LOCALSPIN_save
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
         PCM    = PCM_save
         MCTYPE = MCTYPE_save
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
            ICI0 = 1 ! orbitals transformed to NO, thus CI coefficients cannot be used
         ELSE IF (DOCI) THEN
            ICI0 = 4 ! orbitals are not transformed, thus CI coefficients from CI module are OK
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
     &      ' CPU and wall time for MCSCF :',TUSED,WUSED
         IF (IPRSIR .GE. 0) WRITE (LUPRI,'(/A,2F12.3)')
     &      ' CPU and wall time for MCSCF :',TUSED,WUSED
         CALL FLSHFO(LUW4)
         CALL FLSHFO(LUPRI)
         CALL FLSHFO(LUSTAT)
      END IF
C
      IF (ICONV .EQ. 0) THEN
         NWARN = NWARN + 1
         WRITE (LUERR,7070)
         WRITE (LUPRI,7070)
         IF (LUW4.NE.LUPRI) WRITE (LUW4,7070)
      END IF
 7070 FORMAT(//'@ ******* WARNING: wave function not converged *******')
C
C     -- truncate secondary/virtual orbital space
C
      IF (DO_VIRTRUNC) THEN
         CALL SIR_VIRTRUNC(WORK,LWORK)
      END IF
C
C *** Output sections
C
C *** ORBITAL OUTPUT
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
     &      '@ WARNING: NEVPT2 skipped because MCSCF not converged.'
       ELSE
         JRDMO = 9
         KCMO  = 1
         KWRK1 = KCMO  + NCMOT
         LWRK1 = LWORK - KWRK1
         CALL READMO(WORK(KCMO),JRDMO)
         USEDRC_save = USEDRC
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
     &      ' CPU and wall time for NEVPT2 :',TUSED,WUSED
         IF (IPRSIR .GE. 0) WRITE (LUPRI,'(/A,2F12.3)')
     &      ' CPU and wall time for NEVPT2 :',TUSED,WUSED
         CALL FLSHFO(LUW4)
         CALL FLSHFO(LUPRI)
         CALL FLSHFO(LUSTAT)
         USEDRC = USEDRC_save

       END IF
      END IF
C
C ***************************************************
C     *** LUCITA ***
C
      IF (DOLUCITA) THEN
       IF (ICONV .EQ. 0) THEN
         WRITE (LUPRI,'(//A)') '@ WARNING: '//
     &      'LUCITA skipped because requested HF/MCSCF not converged.'
       ELSE
         addsri_s      = addsri
         dftadd_s      = dftadd
         DOHFSRDFT     = .FALSE.
         ADDSRI        = .FALSE.
C        ... Dont risk adding DFT contributions in later call of FCKMAT
C            The SR-DFT terms needed in CI-DFT is handled by SRFMAT.
C            Maybe this should be added to the generel code and not 
C            just the CIDFT code if one wants to run a CI with KS orbitals !
C            (if the CISRDFT was specified together with a presceding
C             HFSRDFT call, we assume the HFSRDFT is done by now!)
         DFTADD        = .FALSE.
         kfirst_lucita = 1
         JRDMO         = 9
         KFRSAV        = 1
         KFREE         = KFRSAV
         LFREE         = LWORK
         CALL MEMGET('REAL',KCMO,NCMOT,WORK,KFREE,LFREE)
         CALL READMO(WORK(KCMO),JRDMO)

!       define LUCITA control variables and orbital spaces provided via GASCI input (irun_type == '2')
        irun_type =  2 

        xctype1   = -1
        xctype2   = -1

!       step 1: dynamic variables
        call define_lucita_cfg_dynamic(gasci_input_ptg_symmetry,
     &                                 gasci_input_ptg_symmetry,
     &                                  0, ! no istaci; converge all roots being asked for
     &                                  0,
     &                                  0,
     &                                  gasci_input_nr_roots,
     &                                  gasci_input_is_spin_multiplett,
     &                                  gasci_input_max_nr_dav_ci_iter,
     &                                  gasci_input_density_calc_lvl,
     &                                  gasci_input_global_print_lvl,
     &                                  gasci_input_accepted_truncation,
     &                                  gasci_input_convergence_thr,
     &                                  gasci_input_aux_convergence_thr,
     &                                  gasci_input_natural_orb_occ_nr,
     &                                  gasci_input_analyze_cvec,
     &                                  .false.,
     &                                  .false.,
     &                                  .false.,
     &                                  .false.,
     &                                  xctype1, ! vector exchange type1
     &                                  xctype2, ! vector exchange type2
     &                                  .false., ! vector exchange active in parallel runs in mc2lu interface (both mc-->lu and lu-->mc)
     &                                  .false.) ! no vector exchange between lucita/mc using io2io (here: no vector exchange at all...)
!        step 2: static variables
         call define_lucita_cfg_static(irun_type)

         CALL SET_INFORB_LUCITA('LUCITA') ! define inforb.h with LUCITA
                                          ! orbital occupations for TRACTL
         JTRALVL = 0 ! only active orbitals for all four indices
         USEDRC_save = USEDRC
         USEDRC  = .FALSE. ! Dirac format 2-el. integrals not needed

         FTRCTL = .true.
         if (.not.gasci_input_skip_4index_trafo) then

!          integral transformation
           CALL TRACTL(JTRALVL,WORK(KCMO),WORK(KFREE),LFREE)

!          prepare (aa|aa) 1e-/2e- integral file to be read within LUCITA
           CALL MAKE_LUCITA_INTEGRALS(WORK(KCMO),WORK(KFREE),LFREE)
         end if

         CALL MEMREL('lucita-tra.done',work,KFRSAV,KFRSAV,KFREE,LFREE)

!        perform large-scale (GAS)CI calculation with LUCITA
!        ---------------------------------------------------

         CALL GETTIM(TSTART,WSTART)


!        if FCI dump active - nothing else to do...
!        ------------------------------------------
         if(.not.gasci_input_fci_dump)then
!          a1. calculate dimensions of vector and resolution matrices
           call mcscf_lucita_interface(vdummy,vdummy,vdummy,vdummy,
     &                               'return CIdim',work,lfree,iprsir)
!          a2. calculate dimensions of vector and resolution matrices
           call setup_lucita_mcci_wrkspc_dimensions('ijkl resort ',0)
           call mcscf_lucita_interface(vdummy,vdummy,vdummy,vdummy,
     &                                 'ijkl resort ',work,lfree,iprsir)
 
!          b. set vector, resolution matrices and integral array length in common MC/CI module
           call setup_lucita_mcci_wrkspc_dimensions('standard ci ',0)
!          c. allocate vectors (dimensions are provided by the module lucita_mcscf_ci_cfg)
           call memget('REAL',K_VEC1,len_cref_mc2lu,WORK,KFREE,LFREE)
           call memget('REAL',K_VEC2,  len_hc_mc2lu,WORK,KFREE,LFREE)

!          d. ------------------- start of LUCITA machinery ---------------------
           CALL LUCITA('standard ci ',work(k_vec1),work(k_vec2),vdummy,
     &                 vdummy,work(kfree),lfree)
!          ---------------------- end of LUCITA machinery -----------------------

           CALL MEMREL('lucita-ci.done',work,kfrsav,kfrsav,kfree,lfree)
         end if
         CALL GETTIM(TEND,WEND)

!        timing statistics
         TUSED = TEND - TSTART
         WUSED = WEND - WSTART

         IF (IPRSTAT .GT. 0) WRITE (LUSTAT,'(/A,2F12.3)')
     &      ' CPU and wall time for LUCITA CI :',TUSED,WUSED
         IF (IPRSIR .GE. 0) WRITE (LUPRI,'(/A,2F12.3)')
     &      ' CPU and wall time for LUCITA CI :',TUSED,WUSED
         CALL FLSHFO(LUW4)
         CALL FLSHFO(LUPRI)
         CALL FLSHFO(LUSTAT)
         USEDRC = USEDRC_save

         CALL SET_INFORB_LUCITA('RESET ')
         addsri        = addsri_s
         dftadd        = dftadd_s
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

      IF (DOSTEX) THEN
         CALL STXCTL(WORK,LWORK)
      END IF

#if defined(BUILD_GEN1INT)
C... added by Bin Gao, to generate cube files
      if (DO_CUBE) then
        call gen1int_host_get_cube(LFREE, WORK(KFREE), LUPRI, IPRSIR)
        call gen1int_host_cube_finalize()
      end if
#endif
C
 8000 CALL GETTIM(TEND,WEND)
      TIMSIR = TEND - TIMSIR
      WALSIR = WEND - WALSIR
      WRITE (LUPRI,'()')
      CALL TIMTXT(' Total CPU  time used in SIRIUS :',TIMSIR,LUPRI)
      CALL TIMTXT(' Total wall time used in SIRIUS :',WALSIR,LUPRI)
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
     &   /' Total CPU  time used in SIRIUS :',F10.2,' seconds',
     &   /' Total wall time used in SIRIUS :',F10.2,' seconds')
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
      IF (LUPROP .GT. 0) CALL GPCLOSE(LUPROP,'KEEP') ! AOPROPER   - property integrals
      IF (LUSOL  .GT. 0) CALL GPCLOSE(LUSOL ,'KEEP') ! AOSOLINT   - multipole integrals for spherical cavity
      IF (LUPCMD .GT. 0) CALL GPCLOSE(LUPCMD,'KEEP') ! PCMDATA    - PCM data
      IF (LUIT1  .GT. 0) CALL GPCLOSE(LUIT1,'KEEP')  ! SIRIUS.RST - restart info
      IF (LUW4 .NE. LUPRI) CALL GPCLOSE(LUW4,'KEEP')
      IF (LUH2AC .GT. 0) CALL GPCLOSE(LUH2AC,'DELETE') ! H2AC on disk for large CI
C
      CALL QEXIT('SIRCTL')
      RETURN
C *** end of SIRCTL.
      END
C
C  /* Deck sbdorb */
      BLOCK DATA SBDORB
C     Initialize MULD2H in common block /INFORB/
#include "implicit.h"
#include "inforb.h"
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
#include "implicit.h"
#include "cbgetdis.h"
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
#include "implicit.h"
#include "inftra.h"
      DATA FCKTRA_TYPE/-1/, NEWTRA /.FALSE./
      DATA THRP /1.D-15/
      END
C  /* Deck sbddim */
      BLOCK DATA SBDDIM
C 010321 hjaaj
C Initialize variables in /INFDIM/
C
#include "implicit.h"
#include "infdim.h"
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
