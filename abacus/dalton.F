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
#include <single.h>
C
C=========================================================================
C old log:
C010131-vebjornb: Old dalton.F split into dalton.F (program) and
C                 dalgnr.F (subroutines).
C951009-kenneth : Changed name from PSISTAR to DALTON
C950621-vebjornb: First working optimization, GRAND changed to PSISTAR
C950516-vebjornb: First effort to link all programs: GRAND
C=========================================================================
      PROGRAM DALTON2011
C
#include <implicit.h>
#include <dummy.h>
#include <mtags.h>
#include <priunit.h>
#include <maxorb.h>
#include <infpar.h>
#include <infpri.h>
#include <gnrinf.h>
#include <inftap.h>
#include <prpc.h>
Cef begin
#include <incore.h>
Cef end
#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
      LOGICAL FINISH
#endif
      CHARACTER LABEL1*8, STHELP*10
      PARAMETER (LWORK_DEF = INSTALL_WRKMEM)

!     this is the "static" work array
      real(8), allocatable :: WORK(:)

!     memory-id tags for field 1 and 2+LMWORK
      DATA WRKDLM/12345678D0/
C
C
C Start MPI, if parallel, in MPIXINIT; define MYNUM, MASTER, NODTOT
C If not MPI, MPIXINIT defines: MYNUM  = 0; MASTER = 0; NODTOT = 0
C
      CALL MPIXINIT

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

#ifdef SYS_LINUX
      LUERR = 0 ! stderr
#else
      LUERR = 6 ! stdout (which is not DALTON.OUT !)
#endif

C
C     Read environment variables for LMWORK and BASDIR
C     If LMWORK .le. 0 use default value (LWORK_DEF)
C
      LMWORK = 0
      CALL GETMMBAS(LMWORK,BASDIR,LBASDIR)
      IF (LMWORK .LE. 0) LMWORK = LWORK_DEF
C
      NBYTES  = (LMWORK+2) * 8
      XMBYTES = (LMWORK+2) * 8.D0
      XMBYTES = XMBYTES / (1.024D3**2)
      GBYTES  = XMBYTES / 1.024D3
      
      IF (MYNUM .LE. 1) THEN
         IF (GBYTES .GE. 1.D0) THEN
            WRITE(LUERR,'(/A,I12,A,F7.3,A,I4)')
     &        ' Work memory size (LMWORK+2):',
     &        LMWORK+2,' =',GBYTES,' gigabytes; node',MYNUM
         ELSE
            WRITE(LUERR,'(/A,I12,A,F8.2,A,I4)')
     &        ' Work memory size (LMWORK+2):',
     &        LMWORK+2,' =',XMBYTES,' megabytes; node',MYNUM
         END IF
      END IF
 
!     allocate memory using f90 utilities
      allocate(work(lmwork+2),stat=i)

      IF(I .NE. 0)THEN
         WRITE (LUERR,*) MYNUM,': ALLOCATE command failed to allocate'//
     &        ' the requested memory'
         WRITE (LUERR,*) 'Reduce the memory demands and be welcome back'
         CALL QUIT('Failed to allocate memory')
      END IF

Cef begin
C If AOSAVE .eq. .TRUE. then run the incore-version of dalton where integrals
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
#if defined (VAR_MPI)
      CALL IZERO(INDX_SHL,MXTSK)
#endif
C      print *,'dalton.F: After IZERO: Printing NODTOT: ',NODTOT
      IF (MMWORK .NE. MMCORE) THEN
         CALL QUIT('MMCORE must equal MMWORK.')
      END IF
Cef end

      LEN_BASDIR = MAX(1,LNBLNK(BASDIR))
      IF (MYNUM .LE. 1) THEN
         WRITE(LUERR,'(/I4,A/5X,A)')
     &      MYNUM,': Directories for basis set searches:',
     &      BASDIR(1:LEN_BASDIR)
      END IF
C
C     Set memory traps and initialize unit numbers
C
      WORK(1)        = WRKDLM
      WORK(2+LMWORK) = WRKDLM
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
C         "DALTON.OUT#<slave no>" in GPOPEN

      CALL GPOPEN(LUNDPF,'DALTON.PROP','NEW',' ','FORMATTED',
     &            IDUMMY,.FALSE.)
      REWIND(LUNDPF)
      CALL GPOPEN(LUNMPF,'midasifc.prop','NEW',' ','FORMATTED',
     &            IDUMMY,.FALSE.)
      REWIND(LUNMPF)
C
      CALL MEMINI(LUPRI,LUERR)
C
      IF (MYNUM .EQ. MASTER) THEN
C         
         CALL GPOPEN(LUCME,'DALTON.CM','NEW',' ','FORMATTED',
     &            IDUMMY,.FALSE.)
C
         CALL TITLER('DALTON2011 - An electronic structure program',
     &               '*',115)

         WRITE (LUPRI,'(T5,A)')
     &'This is output from DALTON Release 2011 (Rev. 0, May 2011)'

         WRITE (LUPRI,'(3X,A)') SEPARATOR(1:76)
         WRITE (LUPRI,'(T5,A)')
     &'NOTE:',
     &' ',
     &'This is an experimental code for the evaluation of molecular',
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
     &'an appropriate citation would be:',
     &' ',
     &'"Dalton, a molecular electronic structure program,',
     &'Release DALTON2011 (2011), see http://daltonprogram.org"'
C
         WRITE (LUPRI,'(3X,A)') SEPARATOR(1:76)
         WRITE (LUPRI,'(/T5,A//,(T3,A,T29,A,T59,A,T71,A))')
     *'Authors in alphabetical order '//
     &  '(major contribution(s) in parenthesis):',
     *'Celestino Angeli,','University of Ferrara,','Italy',
     &  '(NEVPT2)',
     *'Keld L. Bak,','UNI-C,','Denmark',
     &  '(AOSOPPA, non-adiabatic coupling, magnetic properties)',
     *'Vebjoern Bakken,','University of Oslo,','Norway',
     &  '(DALTON; geometry optimizer, symmetry detection)',
     *'Ove Christiansen,','Aarhus University,','Denmark',
     &  '(CC module)',
     *'Renzo Cimiraglia,','University of Ferrara,','Italy','(NEVPT2)',
     *'Sonia Coriani,','University of Trieste,','Italy',
     &  '(CC module, MCD in RESPONS)',
     *'Paal Dahle,','University of Oslo,','Norway','(Parallelization)',
     *'Erik K. Dalskov,','UNI-C,','Denmark','(SOPPA)',
     *'Thomas Enevoldsen,','SDU - Odense University,','Denmark',
     &  '(SOPPA)',
     *'Berta Fernandez,','U. of Santiago de Compostela,','Spain',
     &  '(doublet spin, ESR in RESPONS)',
     *'Lara Ferrighi,','Aarhus University,','Denmark',
     &  '(PCM Cubic response)',
     *'Heike Fliegl,','University of Helsinki,','Finland',
     &  '(CCSD(R12))',
     *'Luca Frediani,','University of Tromsoe,','Norway','(PCM)',
     *'Christof Haettig,','Ruhr University Bochum,','Germany',
     &  '(CC module)',
     *'Kasper Hald,', 'Aarhus University,','Denmark','(CC module)',
     *'Asger Halkier,','Aarhus University,','Denmark','(CC module)',
     *'Hanne Heiberg,','University of Oslo,','Norway',
     &  '(geometry analysis, selected one-electron integrals)',
     *'Trygve Helgaker,','University of Oslo,','Norway',
     &  '(DALTON; ABACUS, ERI, DFT modules, London, and much more)',
     *'Hinne Hettema,','University of Auckland,','New Zealand',
     &  '(quadratic response in RESPONS; SIRIUS supersymmetry)',
     *'Brano Jansik','University of Aarhus','Denmark',
     &  '(DFT cubic response)',
     *'Hans Joergen Aa. Jensen,','Univ. of Southern Denmark,','Denmark',
     &'(DALTON; SIRIUS, RESPONS, ABACUS modules, London, and much more)'
     *,'Dan Jonsson,','University of Tromsoe,','Norway',
     &  '(cubic response in RESPONS module)',
     *'Poul Joergensen,','Aarhus University,','Denmark',
     &  '(RESPONS, ABACUS, and CC modules)',
     *'Sheela Kirpekar,','SDU - Odense University,','Denmark',
     &  '(Mass-velocity & Darwin integrals)',
     *'Wim Klopper,','University of Karlsruhe,','Germany',
     &  '(R12 code in CC, SIRIUS, and ABACUS modules)',
     *'Stefan Knecht,','Univ. of Southern Denmark,','Denmark',
     &  '(Parallel CI)',
     *'Rika Kobayashi,','ANU Supercomputer Facility,','Australia',
     &  '(DIIS in CC, London in MCSCF)',
     *'Jacob Kongsted,','Univ. of Southern Denmark,','Denmark',
     &  '(QM/MM code)',
     *'Henrik Koch,','University of Trondheim,','Norway',
     &  '(CC module, Cholesky decomposition)',
     *'Andrea Ligabue,','University of Modena,','Italy',
     &  '(CTOCD, AOSOPPA)',
     *'Ola B. Lutnaes,','University of Oslo,','Norway','(DFT Hessian)',
     *'Kurt V. Mikkelsen,','University of Copenhagen,','Denmark',
     &  '(MC-SCRF and QM/MM code)',
     *'Christian Neiss,','Univ. Erlangen-Nürnberg,','Germany',
     &  '(CCSD(R12))',
     *'Christian B. Nielsen,','University of Copenhagen,','Denmark',
     &  '(QM/MM code)',
     *'Patrick Norman,','University of Linkoeping,','Sweden',
     &  '(cubic response and complex response in RESPONS)',
     *'Jeppe Olsen,','Aarhus University,','Denmark',
     &  '(SIRIUS CI/density modules)',
     *'Anders Osted,','Copenhagen University,','Denmark','(QM/MM code)',
     *'Martin J. Packer,','University of Sheffield,','UK',
     &  '(SOPPA)',
     *'Thomas B. Pedersen,','University of Oslo,','Norway',
     &  '(Cholesky decomposition)',
     *'Zilvinas Rinkevicius,','KTH Stockholm,','Sweden',
     &  '(open-shell DFT, ESR)',
     *'Elias Rudberg,','KTH Stockholm,','Sweden',
     &  '(DFT grid and basis info)',
     *'Torgeir A. Ruden,','University of Oslo,','Norway',
     &  '(Numerical derivatives in ABACUS)',
     *'Kenneth Ruud,','University of Tromsoe,','Norway',
     &  '(DALTON; ABACUS magnetic properties and  much more)',
     *'Pawel Salek,','KTH Stockholm,','Sweden','(DALTON; DFT code)',
     *'Claire C.M. Samson','University of Karlsruhe','Germany',
     &  '(Boys localization, r12 integrals in ERI)',
     *'Alfredo Sanchez de Meras,','University of Valencia,','Spain',
     &  '(CC module, Cholesky decomposition)',
     *'Trond Saue,','CNRS/ULP Toulouse,','France',
     &  '(direct Fock matrix construction)',
     *'Stephan P. A. Sauer,','University of Copenhagen,','Denmark',
     &  '(SOPPA(CCSD), SOPPA prop., AOSOPPA, vibrational g-factors)',
     *'Bernd Schimmelpfennig,','Forschungszentrum Karlsruhe,',
     &  'Germany','(AMFI module)',
     *'Arnfinn H. Steindal,','University of Tromsoe,','Norway',
     &  '(parallel QM/MM)',
     *'K. O. Sylvester-Hvid,','University of Copenhagen,','Denmark',
     &  '(MC-SCRF)',
     *'Peter R. Taylor,','VLSCI/Univ. of Melbourne,','Australia',
     &  '(Symmetry handling ABACUS, integral transformation)',
     *'David P. Tew,','University of Bristol,','England',
     &  '(CCSD(R12))',
     *'Olav Vahtras,','KTH Stockholm,','Sweden',
     &  '(triplet response, spin-orbit, ESR, TDDFT, open-shell DFT)',
     *'David J. Wilson,','La Trobe University,','Australia',
     &  '(DFT Hessian and DFT magnetizabilities)',
     *'Hans Agren,','KTH Stockholm,','Sweden',
     &  '(SIRIUS module, MC-SCRF solvation model)'
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
C
C        Read general input for Dalton
C
         CALL DALTON_GNRLINP
C
#if !defined (VAR_MPI)
      END IF
#else /* VAR_MPI */
         CALL RECVNAMES
!
!        call main driver routine
         call exedrv(work,lmwork,wrkdlm)

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
C     Send ending signal to slaves
C
         IPRTYP = 0
         CALL MPIXBCAST(IPRTYP,1,'INTEGER',MYNUM)
         CALL GPCLOSE(LUCME,'KEEP')
         CALL GPCLOSE(LUSTAT,'KEEP')
         CALL GPCLOSE(LUMOLDEN,'KEEP')
         CALL GPCLOSE(LUNDPF,'KEEP')
         CALL GPCLOSE(LUNMPF,'KEEP')
         CALL GPCLOSE(LUPRI,'KEEP')
      ELSE
C     ... I am a "slave" ... (MYNUM .ne. MASTER)
         CALL TITLER('DALTON - An electronic structure program',
     &               '*',111)
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
         IF (PARIO) THEN
            CALL DALTON_GNRLINP

!           call main driver routine (co-workers do it as well if PARIO == .true.)
            call exedrv(work,lmwork,wrkdlm)

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
            CALL GPCLOSE(LUSTAT,'KEEP')
            CALL GPCLOSE(LUMOLDEN,'KEEP')
            CALL GPCLOSE(LUNDPF,'KEEP')
            CALL GPCLOSE(LUNMPF,'KEEP')
            CALL GPCLOSE(LUPRI,'KEEP')
         ELSE
            CALL SENDNAME(MASTER,MYNUM)

!           call main co-workers driver routine
            call nodedriver(work,lmwork,wrkdlm)
         END IF
      END IF

!     release dynamically allocated memory
      deallocate(work)

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
!     allocate memory using f90 utilities
      allocate(work(lmwork+2),stat=i)
      work(1)        = wrkdlm
      work(2+lmwork) = wrkdlm

!     call main co-workers driver routine
      call exedrv(work,lmwork,wrkdlm)

!     release dynamically allocated memory
      deallocate(work)

      LABEL1 = 'ALL_DONE'
      STHELP = 'THE_END   '
      ZERO   = 0.0D0
      CALL WRIPRO(0.0D0,STHELP,666,
     *            LABEL1,LABEL1,LABEL1,LABEL1,
     *            ZERO,ZERO,ZERO,1,0,0,0)
      IPRLOC = MAX(IPRUSR,IPRRHF)
C     IF (IPRLOC .GT. 2) CALL PRPRPC(LUNDPF,1,DUMMY,NPRMI)
      IF (IPRLOC .GT. 2) CALL PRPRPC(LUNMPF,1,DUMMY,NPRMI)
      CALL GPCLOSE(LUSTAT,'KEEP')
      CALL GPCLOSE(LUCME,'KEEP')
      CALL GPCLOSE(LUMOLDEN,'KEEP')
      CALL GPCLOSE(LUNDPF,'KEEP')
      CALL GPCLOSE(LUNMPF,'KEEP')
      CALL GPCLOSE(LUPRI,'KEEP')
#endif
#if defined (VAR_PVM)
      IF (PARCAL) THEN
C
C        Tell slaves to quit.
C        --------------------
C
         IPRTYP = 0
         CALL PVMXINITSEND(NCODE)
         CALL PVMXPACK(IPRTYP,1,'INTEGER')
         CALL PVMXMCAST(NODTOT,NODEID(1),MTAG1)
         CALL PVMXEXIT
      END IF
#endif

      CALL QEXIT('DALTON main')
      END

      SUBROUTINE GETMMBAS(LMWORK,BASDIR,LBASDIR)
#include <implicit.h>
#include <priunit.h>
      CHARACTER*(*) BASDIR
      CHARACTER*20  WRKMEM
#include <maxorb.h>
#include <infpar.h>
C
C     Dynamic memory allocation; check for WRKMEM in environment
C     The T3D thingie might work for DEC_ALPHA as well? Try it out, kr-dec96
C
      IF (.NOT. SLAVE) THEN
C     ... only master reads environment, transferred
C         to slave with MPI_BCAST.
C         On some MPI implementations it requires extra
C         options to transfer UNIX environment variables to
C         slave nodes. /hjaaj June 2005
        WRKMEM = '                    '
#if defined (SYS_T3D) || defined (SYS_T90)
        CALL PXFGETENV('WRKMEM',6,WRKMEM,LWRKMM,IERR)
        CALL PXFGETENV('BASDIR',6,BASDIR,LBASDR,IERR)
#else
        CALL GETENV('WRKMEM',WRKMEM)
        CALL GETENV('BASDIR',BASDIR)
        !debug! write (0,*) 'GETENV WRKMEM ',WRKMEM
        !debug! write (0,*) 'GETENV BASDIR ',BASDIR
#endif
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
          WRITE (LUERR,'(/A)') ' DALTON: default work memory size used.'
        ELSE
C
C       ... or user specified memory size.
C
          LWMEM = MAX(1,LNBLNK(WRKMEM))
          WRITE(LUERR,'(/A/3A)')
     &    ' DALTON: user specified work memory size used,',
     &    '         environment variable WRKMEM = "',WRKMEM(1:LWMEM),'"'
        END IF
      END IF
      IF (LNBLNK(BASDIR) .LT. 1) BASDIR = './:'
#if defined(VAR_MPI)
C
C     Transfer LMWORK and BASDIR to slaves
C
      CALL MPIXBCAST(LMWORK,1,'INTEGE',MASTER)
      CALL MPIXBCAST(BASDIR,LBASDIR,'STRING',MASTER)
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
#include <mxcent.h>
#include <eribuf.h>
#include <priunit.h>
#include <inftap.h>
#include <ccinftap.h>
#include <r12int.h>
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
      LUMOLDEN=-9033
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
C
C     /R12INT/ (WK/UniKA/26-11-2002).
C
      DO I = 1, 5
         LUR12(I) = -9100 - MXCOOR - I
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
