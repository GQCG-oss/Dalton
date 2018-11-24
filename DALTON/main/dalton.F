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
C
C=========================================================================
C old log:
C010131-vebjornb: Old dalton.F split into dalton.F (program) and
C                 dalgnr.F (subroutines).
C951009-kenneth : Changed name from PSISTAR to DALTON
C950621-vebjornb: First working optimization, GRAND changed to PSISTAR
C950516-vebjornb: First effort to link all programs: GRAND
C=========================================================================

! define DALTON_VERSION
#include "dalton_config.h"

#ifdef VAR_CHEMSHELL
      subroutine dalton
#else
      PROGRAM DALTON2017_DALTON
#endif
C
      use parallel_models_mpi
#ifdef ENABLE_XCINT
      use dalton_xcint_interface, only: xcint_is_active, xcint_finalize
#endif
C
#include "implicit.h"
#include "dummy.h"
#include "mtags.h"
#include "priunit.h"
#include "maxorb.h"
#include "infpar.h"
#include "infpri.h"
#include "gnrinf.h"
#include "inftap.h"
#include "prpc.h"
#include "incore.h"

#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
      LOGICAL FINISH
#endif
#ifdef VAR_OMP
      INTEGER omp_get_num_procs
#endif
      CHARACTER LABEL1*8, STHELP*10
C defined parallel calculation types  
#include "iprtyp.h"

!     this is to echo input files to the output
      character(100)     :: line
      logical            :: file_exists, file_open

!     memory-id tags for field 1 and 2+LMWORK in work (allocated in main driver routines)
      DATA WRKDLM/12345678D0/

      ! we use this to read env variables
      character*8 :: env_value

      ! number of OMP threads
      integer :: omp_num_threads
C
C
C Start MPI, if parallel, in MPIXINIT; define MYNUM, MASTER, NODTOT
C If not MPI, MPIXINIT defines: MYNUM  = 0; MASTER = 0; NODTOT = 0
C
      CALL MPIXINIT

      call parallel_models_initialize_mpi()

      CALL QENTER('DALTON main')
      PARIO = .FALSE.
C
C     Define SLAVE variable (only false for MASTER, MYNUM.eq.0)
C     Initialize file pointers; open error file for MASTER
C
      CALL GPIO_INI_DALTON
      IF (MYNUM .EQ. MASTER) THEN
         SLAVE = .FALSE.
         CALL GPOPEN(LUSTAT,'DALTON.STAT','UNKNOWN',' ','FORMATTED',
     &               IDUMMY,.FALSE.)
      ELSE
         SLAVE = .TRUE.
      END IF

#ifdef SYS_UNIX
      LUERR = 0 ! stderr
#else
      LUERR = 6 ! stdout (which is not DALTON.OUT !)
#endif

C
C     Read environment variables for LMWORK, LMWORK_NODES, and BASDIR
C     If LMWORK .le. 0 use default value (LWORK_DEF)
C
      CALL GETMMBAS(LMWORK,BASDIR,LBASDIR)
C
      NBYTES  = (LMWORK+2) * 8
      XMBYTES = (LMWORK+2) * 8.D0
      XMBYTES = XMBYTES / (1.024D3**2)
      GBYTES  = XMBYTES / 1.024D3
      
      IF (MYNUM .LE. 1) THEN
         IF (GBYTES .GE. 1.D0) THEN
            WRITE(*,'(/A,I12,A,F7.3,A,I4)')
     &        ' Work memory size (LMWORK+2):',
     &        LMWORK+2,' =',GBYTES,' gigabytes; node',MYNUM
         ELSE
            WRITE(*,'(/A,I12,A,F8.2,A,I4)')
     &        ' Work memory size (LMWORK+2):',
     &        LMWORK+2,' =',XMBYTES,' megabytes; node',MYNUM
         END IF
         IF (MYNUM .EQ. 1 .AND. NODTOT .GT. 1) THEN
            WRITE(*,'(A,I0,A)')
     &        ' Same work memory size for all ',NODTOT,' slaves.'
         END IF
      END IF
 
Cef begin
C If AOSAVE .eq. .TRUE. then run the incore-version of Dalton where integrals
C are stored in memory on each slave for reuse. Take a look at the header
C file DALTON_SOURCE_DIR/include/incore.h  and the variables: MMWORK and IEIR. 
C MMWORK (we used MMWORK=80000000 when testing ) determines the size of 
C the chunk of array WORK which is reserved for incore-calculations. 
C The variable IEIR in incore.h must also be set to a large enough positive 
C integer (we used IEIR=3000000).

      JSCORE    = 1
      ISCORE    = JSCORE
      LMCORE    = MMWORK 
      MMCORE    = LMCORE
      N_SHL     = 1
      AOSAVE    = (MMWORK.GT.1)
      INITX     = .FALSE.
      MSAVE     = .TRUE.
      LINTSV    = .FALSE.
      LINTMP    = .FALSE.
      INDX_SHL1 = 0
      INDX_SHL2 = 0
      INDX_SHL3 = 0
      INDX_SHL4 = 0
      CALL IZERO(INDX_SHL,MXTSK)
      IF (MMWORK .NE. MMCORE) THEN
         CALL QUIT('MMCORE must equal MMWORK.')
      END IF
Cef end

      LEN_BASDIR = MAX(1,LNBLNK(BASDIR))
      IF (MYNUM .LE. 1) THEN
         WRITE(*,'(/I4,A/5X,A)')
     &      MYNUM,': Directories for basis set searches:',
     &      BASDIR(1:LEN_BASDIR)
      END IF
C
C     Allowing now for parallelization over different molecular geometries,
C     we must first check if this calculation indeed is of that kind.
C     We do this by letting the Master node read DALTON.INP, searching for
C     the keyword .PARNMD. If found, read DALTON/MOLECULE.INP and send them
C     to the slaves, who writes them out in their own files
C
c      IF (NODTOT .GT. 0) THEN
c         IF (MYNUM .EQ. MASTER) THEN
c            CALL PARIOT
c         ELSE
c            CALL PARION
c         END IF
c      END IF
C
      CALL GPOPEN(LUPRI,'DALTON.OUT','NEW',' ','FORMATTED',
     &            IDUMMY,.FALSE.)
      IF (SLAVE) LUSTAT = LUPRI
C     ... GPOPEN will change name from DALTON.OUT to $WRKDIR/$OUTFIL
C         for Master if OUTFIL is an environment variable;
C         for slaves the name will be changed to
C         "DALTON.OUT.n<slave no>" in GPOPEN

      IF (MYNUM.EQ.MASTER) THEN
      CALL GPOPEN(LUNDPF,'DALTON.PROP','NEW',' ','FORMATTED',
     &            IDUMMY,.FALSE.)
      REWIND(LUNDPF)
      CALL GPOPEN(LUNMPF,'midasifc.prop','NEW',' ','FORMATTED',
     &            IDUMMY,.FALSE.)
      REWIND(LUNMPF)
      END IF
C
      CALL MEMINI(LUPRI,LUERR)
C
      IF (MYNUM .EQ. MASTER) THEN
C         
         CALL GPOPEN(LUCME,'DALTON.CM','NEW',' ','FORMATTED',
     &            IDUMMY,.FALSE.)
C
         CALL TITLER('Dalton - An Electronic Structure Program',
     &               '*',115)

         WRITE (LUPRI,'(T5,A)') 'This is output from DALTON release '//
     &DALTON_VERSION
         WRITE (LUPRI,'(T10,A)')'( Web site: http://daltonprogram.org )'
         WRITE (LUPRI,'(/3X,A/)') SEPARATOR(1:76)
         WRITE (LUPRI,'(T5,A)')
     &'NOTE:',
     &' ',
     &'Dalton is an experimental code for the evaluation of molecular',
     &'properties using (MC)SCF, DFT, CI, and CC wave functions.',
     &'The authors accept no responsibility for the performance of',
     &'the code or for the correctness of the results.',
     &' ',
     &'The code (in whole or part) is provided under a licence and',
     &'is not to be reproduced for further distribution without',
     &'the written permission of the authors or their representatives.',
     &' ',
     &'See the home page "http://daltonprogram.org"'
     &//' for further information.',
     &' ',
     &'If results obtained with this code are published,',
     &'the appropriate citations would be both of:',
     &' ',
     &'   K. Aidas, C. Angeli, K. L. Bak, V. Bakken, R. Bast,',
     &'   L. Boman, O. Christiansen, R. Cimiraglia, S. Coriani,',
     &'   P. Dahle, E. K. Dalskov, U. Ekstroem,',
     &'   T. Enevoldsen, J. J. Eriksen, P. Ettenhuber, B. Fernandez,',
     &'   L. Ferrighi, H. Fliegl, L. Frediani, K. Hald, A. Halkier,',
     &'   C. Haettig, H. Heiberg, T. Helgaker, A. C. Hennum,',
     &'   H. Hettema, E. Hjertenaes, S. Hoest, I.-M. Hoeyvik,',
     &'   M. F. Iozzi, B. Jansik, H. J. Aa. Jensen, D. Jonsson,',
     &'   P. Joergensen, J. Kauczor, S. Kirpekar,',
     &'   T. Kjaergaard, W. Klopper, S. Knecht, R. Kobayashi, H. Koch,',
     &'   J. Kongsted, A. Krapp, K. Kristensen, A. Ligabue,',
     &'   O. B. Lutnaes, J. I. Melo, K. V. Mikkelsen, R. H. Myhre,',
     &'   C. Neiss, C. B. Nielsen, P. Norman, J. Olsen,',
     &'   J. M. H. Olsen, A. Osted, M. J. Packer, F. Pawlowski,',
     &'   T. B. Pedersen, P. F. Provasi, S. Reine, Z. Rinkevicius,',
     &'   T. A. Ruden, K. Ruud, V. Rybkin, P. Salek, C. C. M. Samson,',
     &'   A. Sanchez de Meras, T. Saue, S. P. A. Sauer,',
     &'   B. Schimmelpfennig, K. Sneskov, A. H. Steindal,',
     &'   K. O. Sylvester-Hvid, P. R. Taylor, A. M. Teale,',
     &'   E. I. Tellgren, D. P. Tew, A. J. Thorvaldsen, L. Thoegersen,',
     &'   O. Vahtras, M. A. Watson, D. J. D. Wilson, M. Ziolkowski',
     &'   and H. Agren,',
     &'   "The Dalton quantum chemistry program system",',
     &'   WIREs Comput. Mol. Sci. 2014, 4:269–284 '//
     &    '(doi: 10.1002/wcms.1172)',
     &' ',
     &'and',
     &' ',
     &'   Dalton, a molecular electronic structure program,',
     &'   Release '//
     &DALTON_VERSION
     &//', see http://daltonprogram.org'
C
         WRITE (LUPRI,'(3X,A)') SEPARATOR(1:76)
         WRITE (LUPRI,'(/T5,A//,(T3,A,T29,A,T59,A,T71,A))')
     *'Authors in alphabetical order '//
     &  '(major contribution(s) in parenthesis):',
     *'Kestutis Aidas,','Vilnius University,','Lithuania',
     &  '(QM/MM)',
     *'Celestino Angeli,','University of Ferrara,','Italy',
     &  '(NEVPT2)',
     *'Keld L. Bak,','UNI-C,','Denmark',
     &  '(AOSOPPA, non-adiabatic coupling, magnetic properties)',
     *'Vebjoern Bakken,','University of Oslo,','Norway',
     &  '(DALTON; geometry optimizer, symmetry detection)',
     *'Radovan Bast,','UiT The Arctic U. of Norway,','Norway',
     &  '(DALTON installation and execution frameworks)',
     *'Pablo Baudin,','University of Valencia,','Spain',
     &  '(Cholesky excitation energies)',
     *'Linus Boman,','NTNU,', 'Norway',
     &  '(Cholesky decomposition and subsystems)',
     *'Ove Christiansen,','Aarhus University,','Denmark',
     &  '(CC module)',
     *'Renzo Cimiraglia,','University of Ferrara,','Italy',
     &  '(NEVPT2)',
     *'Sonia Coriani,','Technical Univ. of Denmark,','Denmark',
     &  '(CC module, MCD in RESPONS)',
     *'Janusz Cukras,','University of Trieste,','Italy',
     &  '(MChD in RESPONS)',
     *'Paal Dahle,','University of Oslo,','Norway',
     &  '(Parallelization)',
     *'Erik K. Dalskov,','UNI-C,','Denmark',
     &  '(SOPPA)',
     *'Thomas Enevoldsen,','Univ. of Southern Denmark,','Denmark',
     &  '(SOPPA)',
     *'Janus J. Eriksen,','Aarhus University,','Denmark',
     &  '(Polarizable embedding model, TDA)',
     *'Rasmus Faber,','University of Copenhagen,','Denmark',
     &  '(Vib.avg. NMR with SOPPA, parallel AO-SOPPA)',
     *'Tobias Fahleson,','KTH Stockholm,','Sweden',
     &  '(Damped cubic response)',
     *'Berta Fernandez,','U. of Santiago de Compostela,','Spain',
     &  '(doublet spin, ESR in RESPONS)',
     *'Lara Ferrighi,','Aarhus University,','Denmark',
     &  '(PCM Cubic response)',
     *'Heike Fliegl,','University of Oslo,','Norway',
     &  '(CCSD(R12))',
     *'Luca Frediani,','UiT The Arctic U. of Norway,','Norway',
     &  '(PCM)',
     *'Bin Gao,','UiT The Arctic U. of Norway,','Norway',
     &  '(Gen1Int library)',
     *'Christof Haettig,','Ruhr-University Bochum,','Germany',
     &  '(CC module)',
     *'Kasper Hald,', 'Aarhus University,','Denmark',
     &  '(CC module)',
     *'Asger Halkier,','Aarhus University,','Denmark',
     &  '(CC module)',
     *'Frederik Beyer Hansen,','University of Copenhagen,','Denmark',
     &  '(Parallel AO-SOPPA)',
     *'Erik D. Hedegaard,','Univ. of Southern Denmark,','Denmark',
     & '(Polarizable embedding model, QM/MM)',
     *'Hanne Heiberg,','University of Oslo,','Norway',
     &  '(geometry analysis, selected one-electron integrals)',
     *'Trygve Helgaker,','University of Oslo,','Norway',
     &  '(DALTON; ABACUS, ERI, DFT modules, London, and much more)',
     *'Alf Christian Hennum,','University of Oslo,','Norway',
     &  '(Parity violation)',
     *'Hinne Hettema,','University of Auckland,','New Zealand',
     &  '(quadratic response in RESPONS; SIRIUS supersymmetry)',
     *'Eirik Hjertenaes,', 'NTNU,','Norway',
     &  '(Cholesky decomposition)',
     *'Pi A. B. Haase,','University of Copenhagen,','Denmark',
     &  '(Triplet AO-SOPPA)',
     *'Maria Francesca Iozzi,','University of Oslo,','Norway',
     &  '(RPA)',
     *'Christoph Jacob','TU Braunschweig',
     &  'Germany','(Frozen density embedding model)',
     *'Brano Jansik','Technical Univ. of Ostrava','Czech Rep.',
     &  '(DFT cubic response)',
     *'Hans Joergen Aa. Jensen,','Univ. of Southern Denmark,','Denmark',
     &'(DALTON; SIRIUS, RESPONS, ABACUS modules, London, and much more)'
     *,'Dan Jonsson,','UiT The Arctic U. of Norway,','Norway',
     &  '(cubic response in RESPONS module)',
     *'Poul Joergensen,','Aarhus University,','Denmark',
     &  '(RESPONS, ABACUS, and CC modules)',
     *'Maciej Kaminski,','University of Warsaw,','Poland',
     &  '(CPPh in RESPONS)',
     *'Joanna Kauczor,','Linkoeping University,','Sweden',
     &  '(Complex polarization propagator (CPP) module)',
     *'Sheela Kirpekar,','Univ. of Southern Denmark,','Denmark',
     &  '(Mass-velocity & Darwin integrals)',
     *'Wim Klopper,','KIT Karlsruhe,','Germany',
     &  '(R12 code in CC, SIRIUS, and ABACUS modules)',
     *'Stefan Knecht,','ETH Zurich,','Switzerland',
     &  '(Parallel CI and MCSCF)',
     *'Rika Kobayashi,','Australian National Univ.,','Australia',
     &  '(DIIS in CC, London in MCSCF)',
     *'Henrik Koch,','NTNU,', 'Norway',
     &  '(CC module, Cholesky decomposition)',
     *'Jacob Kongsted,','Univ. of Southern Denmark,','Denmark',
     &  '(Polarizable embedding model, QM/MM)',
     *'Andrea Ligabue,','University of Modena,','Italy',
     &  '(CTOCD, AOSOPPA)',
     *'Nanna H. List','Univ. of Southern Denmark,','Denmark',
     &  '(Polarizable embedding model)',
     *'Ola B. Lutnaes,','University of Oslo,','Norway',
     &  '(DFT Hessian)',
     *'Juan I. Melo,','University of Buenos Aires,','Argentina',
     &'(LRESC, Relativistic Effects on NMR Shieldings)',
     *'Kurt V. Mikkelsen,','University of Copenhagen,','Denmark',
     &  '(MC-SCRF and QM/MM)',
     *'Rolf H. Myhre,', 'NTNU,','Norway',
     &  '(Subsystems and CC3)',
     *'Christian Neiss,','Univ. Erlangen-Nuernberg,','Germany',
     &  '(CCSD(R12))',
     *'Christian B. Nielsen,','University of Copenhagen,','Denmark',
     &  '(QM/MM)',
     *'Patrick Norman,','KTH Stockholm,','Sweden',
     &  '(Cubic response and complex frequency response in RESPONS)',
     *'Jeppe Olsen,','Aarhus University,','Denmark',
     &  '(SIRIUS CI/density modules)',
     *'Jogvan Magnus H. Olsen,','Univ. of Southern Denmark,', 'Denmark',
     &  '(Polarizable embedding model, QM/MM)',
     *'Anders Osted,','Copenhagen University,','Denmark',
     &  '(QM/MM)',
     *'Martin J. Packer,','University of Sheffield,','UK',
     &  '(SOPPA)',
     *'Filip Pawlowski,','Kazimierz Wielki University,','Poland',
     &  '(CC3)',
     *'Morten N. Pedersen,','Univ. of Southern Denmark,','Denmark',
     &  '(Polarizable embedding model)',
     *'Thomas B. Pedersen,','University of Oslo,','Norway',
     &  '(Cholesky decomposition)',
     *'Patricio F. Provasi,','University of Northeastern,','Argentina',
     &  '(Analysis of coupling constants in localized orbitals)',
     *'Zilvinas Rinkevicius,','KTH Stockholm,','Sweden',
     &  '(open-shell DFT, ESR)',
     *'Elias Rudberg,','KTH Stockholm,','Sweden',
     &  '(DFT grid and basis info)',
     *'Torgeir A. Ruden,','University of Oslo,','Norway',
     &  '(Numerical derivatives in ABACUS)',
     *'Kenneth Ruud,','UiT The Arctic U. of Norway,','Norway',
     &  '(DALTON; ABACUS magnetic properties and much more)',
     *'Pawel Salek,','KTH Stockholm,','Sweden',
     &  '(DALTON; DFT code)',
     *'Claire C. M. Samson','University of Karlsruhe','Germany',
     &  '(Boys localization, r12 integrals in ERI)',
     *'Alfredo Sanchez de Meras,','University of Valencia,','Spain',
     &  '(CC module, Cholesky decomposition)',
     *'Trond Saue,','Paul Sabatier University,','France',
     &  '(direct Fock matrix construction)',
     *'Stephan P. A. Sauer,','University of Copenhagen,','Denmark',
     &  '(SOPPA(CCSD), SOPPA prop., AOSOPPA, vibrational g-factors)',
     *'Andre S. P. Gomes,','CNRS/Universite de Lille,',
     &  'France','(Frozen density embedding model)',
     *'Bernd Schimmelpfennig,','Forschungszentrum Karlsruhe,',
     &  'Germany','(AMFI module)',
     *'Kristian Sneskov,','Aarhus University,','Denmark',
     &  '(Polarizable embedding model, QM/MM)',
     *'Arnfinn H. Steindal,','UiT The Arctic U. of Norway,',
     &  'Norway','(parallel QM/MM, Polarizable embedding model)',
     *'Casper Steinmann,','Univ. of Southern Denmark,','Denmark',
     &  '(QFIT, Polarizable embedding model)',
     *'K. O. Sylvester-Hvid,','University of Copenhagen,','Denmark',
     &  '(MC-SCRF)',
     *'Peter R. Taylor,','VLSCI/Univ. of Melbourne,','Australia',
     &  '(Symmetry handling ABACUS, integral transformation)',
     *'Andrew M. Teale,','University of Nottingham,','England',
     &  '(DFT-AC, DFT-D)',
     *'David P. Tew,','University of Bristol,','England',
     &  '(CCSD(R12))',
     *'Olav Vahtras,','KTH Stockholm,','Sweden',
     &  '(triplet response, spin-orbit, ESR, TDDFT, open-shell DFT)',
     *'Lucas Visscher,','Vrije Universiteit Amsterdam,',
     &  'Netherlands','(Frozen density embedding model)',
     *'David J. Wilson,','La Trobe University,','Australia',
     &  '(DFT Hessian and DFT magnetizabilities)',
     *'Hans Agren,','KTH Stockholm,','Sweden',
     &  '(SIRIUS module, RESPONS, MC-SCRF solvation model)'
C
         WRITE (LUPRI,'(1X,A)') SEPARATOR
C
C     Stamp date and time and hostname to output
C
         CALL TSTAMP('INIT',LUPRI)
         IF (GBYTES .GE. 1.D0) THEN
           WRITE(LUPRI,'(/A,I12,A,F7.3,A)')
     &        ' * Work memory size             :',
     &        LMWORK,' =',GBYTES,' gigabytes.'
         ELSE
           WRITE(LUPRI,'(/A,I12,A,F8.2,A)')
     &        ' * Work memory size             :',
     &        LMWORK,' =',XMBYTES,' megabytes.'
         END IF
         IF (AOSAVE) WRITE(LUPRI,'(A,I12)')
     &        ' + memory for in-core integrals :',MMCORE
C
         WRITE(LUPRI,'(/A)')
     &   ' * Directories for basis set searches:'
         J_END = LNBLNK(BASDIR)
         IF (J_END .LE. 0) THEN
            WRITE(LUPRI,'(5X,A)') 'NONE SPECIFIED!'
         ELSE
C           structure of BASDIR: 'dir1:dir2:dir3'
            I = 0
            I_END = -1
   20       CONTINUE
               I_ST = I_END + 2
               I_END = INDEX(BASDIR(I_ST:J_END),':') - 1
               IF (I_END .LE. 0) THEN
                  I_END = J_END
               ELSE
                  I_END = I_ST - 1 + I_END
               END IF
               IF (I_END .GT. I_ST) THEN
                  I = I + 1
                  WRITE(LUPRI,'(I4,2A)') I,') ',BASDIR(I_ST:I_END)
                  GO TO 20
               END IF
         END IF

         write(lupri, '(//a/a/)')
     &     'Compilation information',
     &     '-----------------------'
         call print_binary_info(lupri)

#if defined (VAR_MPI)
         IF (NODTOT .GT. 0) THEN
            WRITE(LUPRI,'(/A,I0,A)')
     &      " * MPI parallel run using ", NODTOT+1, " processes."
         ELSE
            WRITE(LUPRI,'(/A,I5,A)')
     &      " * Sequential calculation."
         END IF
#endif

#ifdef VAR_OMP
         call getenv('OMP_NUM_THREADS', env_value)
         read(env_value, '(i6)') omp_num_threads
         if (omp_num_threads > 1) then
            WRITE(LUPRI,'(/A,I0,A)')
     &      " * OpenMP run using ", omp_num_threads, " threads."
         end if
#endif

!        print content of input file to output
         inquire(file = 'DALTON.INP', exist = file_exists,
     &           opened = file_open)
         if (file_exists) then
            if (.not. file_open) then
               LUCMD = -1
               CALL GPOPEN(LUCMD,'DALTON.INP','OLD',' ','FORMATTED',
     &                     IDUMMY,.FALSE.)
            end if
            rewind(LUCMD)
         else
            call quit('DALTON.INP file not found')
         end if
         write(lupri, '(//A/A/)')
     &      '   Content of the .dal input file',
     &      ' ----------------------------------'
         do while (.true.)
            read(LUCMD, '(a)', end=31) line
            write(lupri, '(a)') trim(line)
         end do
 31      continue
         if (.not. file_open) CALL GPCLOSE(LUCMD,'KEEP')

!        print content of molecule file to output
         inquire(file = 'MOLECULE.INP', exist = file_exists)
         if (file_exists) then
            LUMOL = -1
            CALL GPOPEN(LUMOL,'MOLECULE.INP','OLD',' ','FORMATTED',
     &                  IDUMMY,.FALSE.)
            write(lupri, '(//a/a/)')
     &         '   Content of the .mol file',
     &         ' ----------------------------'
            rewind(LUMOL)
            do while (.true.)
               read(LUMOL, '(a)', end=32) line
               write(lupri, '(a)') trim(line)
            end do
 32         continue
            CALL GPCLOSE(LUMOL,'KEEP')
         else
!           call quit('MOLECUL.INP file not found')
! somebody forgot with this 'call quit' that it is OK
! to have the "MOLECULE.INP" input inside the "DALTON.INP" file,
! thus the *nevpt2* tests fail with this 'call quit' uncommented.
         end if

C        Read general input for Dalton
C
         CALL DALTON_GNRLINP

#if !defined (VAR_MPI)
      END IF ! mynum .eq. master
#else /* VAR_MPI */
         CALL RECVNAMES

         CALL MOLDEN_HEAD
!
!        call main driver routine
         call dalton_exedrv(lmwork,wrkdlm)

         CALL MOLDEN_TAIL

         CALL FLSHFO(LUPRI)

         LABEL1 = 'ALL_DONE'
         STHELP = 'THE_END   '
         ZERO   = 0.0D0
         CALL WRIPRO(0.0D0,STHELP,666,
     *               LABEL1,LABEL1,LABEL1,LABEL1,
     *               ZERO,ZERO,ZERO,1,0,0,0)
         IPRLOC = MAX(IPRUSR,IPRRHF)
C        IF (IPRLOC .GT. 2) CALL PRPRPC(LUNDPF,1,DUMMY,NPRMI)
         IF (IPRLOC .GT. 2) CALL PRPRPC(LUNMPF,1,DUMMY,NPRMI)
C
C
C     Send ending signal to slaves
C
         IPRTYP = ENDING_SIGNAL
         CALL MPIXBCAST(IPRTYP,1,'INTEGER',MYNUM)

         CALL GPCLOSE(LUCME,'KEEP')
         CALL GPCLOSE(LUNDPF,'KEEP')
         CALL GPCLOSE(LUNMPF,'KEEP')
         CALL GPCLOSE(LUPRI,'KEEP')
         CALL GPCLOSE(LUSTAT,'KEEP')
      ELSE
C     ... I am a "slave" ... (MYNUM .ne. MASTER)
C     (only print on first 10 slaves)
!        IF (MYNUM .LE. 10) THEN
           CALL TITLER('DALTON - An electronic structure program',
     &               '*',111)
           WRITE(LUPRI,'(/A)') ' * Release: '//
     &     DALTON_VERSION
           WRITE(LUPRI,'(/A,I6)') ' * Output from slave no.',MYNUM
           CALL TSTAMP('INIT',LUPRI)
           IF (GBYTES .GE. 1.D0) THEN
             WRITE(LUPRI,'(/A,I12,A,F7.3,A)')' * Work memory size :',
     &        LMWORK,' =',GBYTES,' gigabytes.'
           ELSE
             WRITE(LUPRI,'(/A,I12,A,F8.2,A)')' * Work memory size :',
     &        LMWORK,' =',XMBYTES,' megabytes.'
           END IF
           IF (AOSAVE) THEN
             WRITE(LUPRI,'(A,I12)')' + memory for in-core integrals :',
     &        MMCORE
           END IF
           CALL FLSHFO(LUPRI)
!        END IF
         IF (PARIO) THEN
            CALL DALTON_GNRLINP

            CALL MOLDEN_HEAD

!           call main driver routine (co-workers do it as well if PARIO == .true.)
            call dalton_exedrv(lmwork,wrkdlm)

            CALL MOLDEN_TAIL

            LABEL1 = 'ALL_DONE'
            STHELP = 'THE_END   '
            ZERO   = 0.0D0
            CALL WRIPRO(0.0D0,STHELP,666,
     *                  LABEL1,LABEL1,LABEL1,LABEL1,
     *                  ZERO,ZERO,ZERO,1,0,0,0)
            IPRLOC = MAX(IPRUSR,IPRRHF)
C           IF (IPRLOC .GT. 2) CALL PRPRPC(LUNDPF,1,DUMMY,NPRMI)
            IF (IPRLOC .GT. 2) CALL PRPRPC(LUNMPF,1,DUMMY,NPRMI)
C
C     Send ending signal to slaves
C
            FINISH = .TRUE.
            CALL MPI_BCAST(FINISH,1,MPI_LOGICAL,NODEID(MASTER),
     &           MPI_COMM_WORLD,IERR)
            CALL GPCLOSE(LUCME,'KEEP')
            CALL GPCLOSE(LUNDPF,'KEEP')
            CALL GPCLOSE(LUNMPF,'KEEP')
            CALL GPCLOSE(LUSTAT,'KEEP')
            LUSTAT = LUPRI ! LUSTAT may be used in GPCLOSE
            CALL GPCLOSE(LUPRI,'KEEP')
         ELSE
            CALL SENDNAME(MASTER,MYNUM)

!           call main co-workers driver routine
            call dalton_nodedriver(lmwork,wrkdlm)
         END IF
      END IF

#if defined(BUILD_GEN1INT)
C...  added by Bin Gao, May 13, 2012
C...  terminates Gen1Int interface
      call gen1int_host_finalize()
#endif

      call parallel_models_finalize_mpi(nodtot+1)

      CALL MPIXFINALIZE
C
C     All below is not relevant for MPI as all proc. has stopped
C     with CALL MPI_FINALIZE(IERR)
C
#endif /* VAR_MPI */
#if defined (VAR_PVM)
C
C     Spawn all PVM processes, as we now know how many there are
C
      IF (PARCAL) THEN
         CALL PVMXMYTID(MYTID)
         CALL PVMXSPAWN('node.x',NODTOT,NODEID(1))
         CALL RECVNAMES
         NODEID(0) = MYTID
      END IF
#endif
#if !defined (VAR_MPI)

      CALL MOLDEN_HEAD

!     call main co-workers driver routine
      call dalton_exedrv(lmwork,wrkdlm)

      CALL MOLDEN_TAIL

      LABEL1 = 'ALL_DONE'
      STHELP = 'THE_END   '
      ZERO   = 0.0D0
      CALL WRIPRO(0.0D0,STHELP,666,
     *            LABEL1,LABEL1,LABEL1,LABEL1,
     *            ZERO,ZERO,ZERO,1,0,0,0)
      IPRLOC = MAX(IPRUSR,IPRRHF)
C     IF (IPRLOC .GT. 2) CALL PRPRPC(LUNDPF,1,DUMMY,NPRMI)
      IF (IPRLOC .GT. 2) CALL PRPRPC(LUNMPF,1,DUMMY,NPRMI)
      CALL GPCLOSE(LUCME,'KEEP')
      CALL GPCLOSE(LUNDPF,'KEEP')
      CALL GPCLOSE(LUNMPF,'KEEP')
      CALL GPCLOSE(LUSTAT,'KEEP')
      LUSTAT = LUPRI ! LUSTAT may be used in GPCLOSE
      CALL GPCLOSE(LUPRI,'KEEP')
#if defined(BUILD_GEN1INT)
C...  added by Bin Gao, Oct. 2, 2011
C...  terminates Gen1Int interface
      call gen1int_host_finalize()
#endif
#endif
#if defined (VAR_PVM)
      IF (PARCAL) THEN
C
C        Tell slaves to quit.
C        --------------------
C
         IPRTYP = ENDING_SIGNAL
         CALL PVMXINITSEND(NCODE)
         CALL PVMXPACK(IPRTYP,1,'INTEGER')
         CALL PVMXMCAST(NODTOT,NODEID(1),MTAG1)
         CALL PVMXEXIT
      END IF
#endif


      ! radovan: close (hopefully) all open files
      !          otherwise we run out of units
      !          if we execute Dalton many times in sequence
      !          (ChemShell for example)
      do iunit = 1, 99
         inquire(iunit, opened=file_open)
         if (file_open) close(iunit, status='keep')
      end do
      ! clear the corresponding table and hope for the best
      call clear_iuntab()


#ifdef ENABLE_XCINT
      if (xcint_is_active()) then
        call xcint_finalize()
      end if
#endif


      CALL QEXIT('DALTON main')
      END

      SUBROUTINE GETMMBAS(LMWORK,BASDIR,LBASDIR)

!     module dependencies
      use parallel_communication_models_mpi

#include "implicit.h"
#include "priunit.h"
      CHARACTER*(*) BASDIR
      CHARACTER*20  WRKMEM, NODE_WRKMEM
#include "maxorb.h"
#include "infpar.h"
!     this is to scan the input file for asynchronous memory allocation (RMA model parallelization)
      character(100)     :: line
      logical            :: file_exists, file_open, rma_allocation
C     If LMWORK .le. 0 use default value (LWORK_DEF)
      PARAMETER (LWORK_DEF = INSTALL_WRKMEM)

      rma_allocation = .false.
C
C     Dynamic memory allocation; check for WRKMEM and NODE_WRKMEM in environment
C
      IF (.NOT. SLAVE) THEN
C     ... only master reads environment, transferred
C         to slave with MPI_BCAST.
C         On some MPI implementations it requires extra
C         options to transfer UNIX environment variables to
C         slave nodes. /hjaaj June 2005
        WRKMEM = ' '
        NODE_WRKMEM = ' '
#if defined (SYS_T3D) || defined (SYS_T90)
        CALL PXFGETENV('WRKMEM',6,WRKMEM,LWRKMM,IERR)
        CALL PXFGETENV('NODE_WRKMEM',11,NODE_WRKMEM,LWRKMM,IERR)
        CALL PXFGETENV('BASDIR',6,BASDIR,LBASDR,IERR)
#else
        CALL GETENV('WRKMEM',WRKMEM)
        CALL GETENV('NODE_WRKMEM',NODE_WRKMEM)
        CALL GETENV('BASDIR',BASDIR)
        !debug! write (0,*) 'GETENV WRKMEM ',WRKMEM
        !debug! write (0,*) 'GETENV BASDIR ',BASDIR
#endif

!       WRKMEM -> LMWORK

        READ(WRKMEM, '(I20)', ERR=10) LMWORK
        IF (LMWORK .GE. 0) GO TO 20
   10     WRITE(LUERR,'(/3A/A,I20)')
     &    ' DALTON: WRKMEM conversion error; WRKMEM = "',WRKMEM,'"',
     &    ' DALTON: read as LMWORK = ',LMWORK
          LMWORK = 0
   20   CONTINUE
C
C       Either default memory size,
C
        IF (LMWORK .LE. 0) THEN
          LMWORK = LWORK_DEF
          WRITE (*,'(/A,I20)')
     &      ' DALTON: default work memory size used.',LMWORK
        ELSE
C
C       ... or user specified memory size.
C
          LWMEM = MAX(1,LNBLNK(WRKMEM))
          WRITE(*,'(/A/3A)')
     &    ' DALTON: user specified work memory size used,',
     &    '         environment variable WRKMEM = "',WRKMEM(1:LWMEM),'"'
        END IF

!       NODE_WRKMEM -> LMWORK_NODES

      IF (NODTOT .GT. 0) THEN ! not relevant if no slaves
        READ(NODE_WRKMEM, '(I20)', ERR=60) LMWORK_NODES
        IF (LMWORK_NODES .GE. 0) GO TO 70
   60     WRITE(LUERR,'(/3A/A,I20)')
     &    ' DALTON: NODE_WRKMEM conversion error; NODE_WRKMEM = "',
     &      NODE_WRKMEM,'"',
     &    ' DALTON: read as LMWORK_NODES = ',LMWORK_NODES
          LMWORK_NODES = 0
   70   CONTINUE
C
C       Either default memory size,
C
        IF (LMWORK_NODES .LE. 0) THEN
          LMWORK_NODES = LMWORK
          WRITE (*,'(/A,I20)')
     &    ' DALTON: master work memory size also used for slaves.',
     &    LMWORK_NODES
        ELSE
C
C       ... or user specified memory size.
C
          LWMEM = MAX(1,LNBLNK(NODE_WRKMEM))
          WRITE(*,'(/A/3A)')
     &    ' DALTON: user specified work memory size used,',
     &    '         environment variable NODE_WRKMEM = "',
     &    NODE_WRKMEM(1:LWMEM),'"'
        END IF
      END IF ! NODTOT .gt. 0

      END IF ! .not. slave

      IF (LNBLNK(BASDIR) .LT. 1) BASDIR = './:'
#ifdef VAR_MPI

!     scan input for RMA (remote memory access model parallelization) keyword
!     reason: we need asynchronous memory allocation as this is our default RMA scheme in:
!             --> hermit
!             --> mcscf (soon)
!             --> lucita/gasci (soon)
      if(mynum == master)then
!       print content of input file to output
        inquire(file='DALTON.INP',exist=file_exists,opened=file_open)
        if(file_exists) then
          if(.not. file_open) then
            LUCMD = -1
            CALL GPOPEN(LUCMD,'DALTON.INP','OLD',' ','FORMATTED',
     &                  IDUMMY,.FALSE.)
          end if
            rewind(LUCMD)
         else
            call quit('GETMMBAS: DALTON.INP file not found')
        end if
!       scan for .RMA-MD
        do while (.true.)
           read(LUCMD, '(a50)', end=81) line
           if(line(1:7) == '.RMA-MD')then
             rma_allocation = .true.
             goto 81
           end if
        end do
 81     continue
        if (.not. file_open) CALL GPCLOSE(LUCMD,'KEEP')
      end if
!     end of keyword scan

!     Transfer LMWORK_NODES and BASDIR to slaves
!     ------------------------------------------

      CALL MPIXBCAST(LMWORK_NODES,1,'INTEGE',MASTER)
      CALL MPIXBCAST(LMWORK      ,1,'INTEGE',MASTER)
      CALL MPIXBCAST(rma_allocation,1,'LOGICA',MASTER)
      CALL MPIXBCAST(BASDIR,LBASDIR,'STRING',MASTER)

      IF(SLAVE)THEN
        if(rma_allocation)then
!         node/NUMA master get the same memory as the master only real slaves get less
          if(communication_info_mpi%my_shmem_node_id /= 0)then 
            LMWORK = LMWORK_NODES
          end if
        else
          LMWORK = LMWORK_NODES
        end if
      END IF
#endif /* defined(VAR_MPI) */
      END
C  /* Deck GPIO_INI_DALTON */
      SUBROUTINE GPIO_INI_DALTON
C
C 1-Mar-2000 K.Ruud
C 10-Apr-2000 HJAaJ: changed initialization to -9xxx from 0,
C                    because unit 0 is stderr on many systems.
C                    The number -9xxx is used for easier debugging.
C
#include "mxcent.h"
#include "eritap.h"
#include "priunit.h"
#include "inftap.h"
#include "ccinftap.h"
#include "r12int.h"
C
C     priunit.h:
C
      LUERR  = -8000
      LUSTAT = -8001
      LUW4   = -8004
      LUCMD  = -8005
      LUPRI  = -8006
      LUPOT  = -8009
      NINFO  = 0
      NWARN  = 0
      IPRSTAT= 0
C
C     inftap.h:
C
      LUCME  = -9003
      LUMOL  = -9004
      LUPROP = -9005
      LUSOL  = -9006
      LUINTA = -9007
      LUONEL = -9008
      LUSUPM = -9009
      LUTLM  = -9010
      LUDA1  = -9011
      LUITMP = -9012
      LU2DER = -9013
      LUDASP = -9014
      LURDR  = -9016
      LURDI  = -9017
      LUGDR  = -9018
      LUGDI  = -9019
      LUGDT  = -9020
      LURDT  = -9021
      LUDFCK = -9022
      LUSFDA = -9023
      LUFDC  = -9024
      LUWLK  = -9025
      LUPAO  = -9026
      LUPAS  = -9027
      LUNR1  = -9028
      LUNR3  = -9029
      LUNR5  = -9030
      LUINTR = -9032
      LUPMOM = -9034
      LUMOM  = -9035
      LUEIND = -9036
      LUENUC = -9037
      LUESITE= -9038
      LUEOBAR= -9039
      LUVDWSE= -9040
      LUENSA = -9041
      LUQM3E = -9042
      LUQM3P = -9043
      LUOSCR = -9044
      LUMMOV = -9045
      LUOVER = -9046
      LUNDPF = -9047
      LUSRINT =-9048
      LUNMPF = -9049
      LUMOLDEN    =-9050
      LUMOLDEN_MOS=-9051
C
C     /R12INT/ (WK/UniKA/26-11-2002).
C
      DO I = 1, 5
         LUR12(I) = -9100 - I
      ENDDO
C
C     /ERITAP/
C
      DO I = 0, MXCOOR
         LUAORC(I) = -10000 - I
      ENDDO
C
C     /RSPTAP/
C
      LUAHSO = -9201
      LUCRV1 = -9202
      LUCRV2 = -9203
      LUXYVE = -9204
      LUCRVE = -9205
      LURSP3 = -9206
      LURSP4 = -9207
      LURSP5 = -9208
      LUMHSO = -9209
      LURSP  = -9210
C
C     /SIRTAP/
C
      LUINTM = -9301
      LUIT1  = -9302
      LUIT2  = -9303
      LUIT3  = -9304
      LUIT5  = -9305
      LUINF  = -9306
      LUH2AC = -9307
      LUSIFC = -9308
C     next two for NEWTRA
      LUORDA = -9309
      LUMINT = -9310
C
      FNSOL  = 'AOSOLINT'
      ABARDR = 'ABACUS.RD'
      ABARDI = 'ABACUS.RDI'
      ABAGDR = 'ABACUS.GD'
      ABAGDI = 'ABACUS.GDI'
      ABAGDT = 'ABACUS.GDT'
      ABARDT = 'ABACUS.RDT'
      ABADFK = 'ABACUS.DFK'
      ABASF  = 'ABACUS.SF'
      ABATLM = 'ABACUS.TLM'
      ABAWLK = 'DALTON.WLK'
      ABAIRC = 'DALTON.IRC'
      ABATRJ = 'DALTON.TRJ'
      ABANR1 = 'ABAENR.RST'
      ABANR3 = 'ABAENR.BVC'
      ABANR5 = 'ABAENR.SVC'
      FNINTM = 'MOTWOINT'
      FNSUPM = 'AOSUPINT'
      FNONEL = 'AOONEINT'
      FNSIFC = 'SIRIFC'
      LBSIFC = 'SIR IPH '
C
      RETURN
      END
