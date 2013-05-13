!> @file
!> Contains single subroutine for initializing lsitem and config structures 
!> from information in DALTON.INP and MOLECULE.INP

!> \brief Initilize lsitem and config structures from DALTON.INP and MOLECULE.INP
!> \author Kasper Kristensen (based on the old lsdalton subroutine)
!> \date December 2011
module init_lsdalton_mod
  use precision
  use configurationType, only: configitem
  use typedeftype, only: lsitem
  use matrix_module, only: matrix
  use configuration, only: config_set_default_config, config_read_input, &
       & set_final_config_and_print
  use files, only: lsopen,lsclose
  use linsca_debug, only: sparsetest
  use lstiming, only: lstimer
  use daltoninfo, only: ls_init
  use IIDFTINT, only: II_DFTDISP
#ifdef HAVE_BSM
  use SCFLOOP_module, only: ls_initbsm
#endif
  use matrix_operations, only: mat_no_of_matmuls, mat_pass_info, no_of_matmuls
#ifdef THIS_IS_CMAKE_BUILD
  use gitrevinfo, only: print_git_revision_info
 ! use compinfo, only: print_compilation_info
#endif
  use lsmpi_type, only: lsmpi_finalize, lsmpi_print
contains

!> \brief Initilize lsitem and config structures from DALTON.INP and MOLECULE.INP
!> \author Kasper Kristensen (based on the old lsdalton subroutine)
!> \date December 2011
SUBROUTINE init_lsdalton_and_get_lsitem(lupri,luerr,nbast,ls,config,mem_monitor)

  implicit none
  !> Logical unit number for LSDALTON.OUT
  integer,intent(in) :: lupri
  !> Logical unit number for LSDALTON.ERR
  integer,intent(in) :: luerr
  !> Number of basis functions
  integer,intent(inout) :: nbast
  !> lsitem structure initialized according to information in DALTON.INP and MOLECULE.INP
  TYPE(lsitem),target,intent(inout) :: ls
  !> configuration structure initialized according to information in DALTON.INP and MOLECULE.INP
  type(configItem),intent(inout)    :: config
  !> Monitor memory - mostly for debugging, set to false by default
  logical, intent(inout) :: mem_monitor
  logical             :: doDFT
  REAL(REALK) :: TIMSTR,TIMEND
  real(realk) :: DUMMY(1,1)

  ! Initializations 
#ifdef VAR_DEBUGINT
  mem_monitor = .true.  !Mostly for memory debugging
#else
  mem_monitor = .false. !Mostly for memory debugging
#endif
  CALL PRINT_INTRO(LUPRI)
  call lsmpi_print(lupri)
#ifdef THIS_IS_CMAKE_BUILD
  !call print_compilation_info(lupri)
  call print_git_revision_info(lupri)
#endif

  ! Timing of individual steps
  CALL LSTIMER('START ',TIMSTR,TIMEND,lupri)

  call mat_no_of_matmuls(no_of_matmuls)

  ! Set default configurations for config structure
  call config_set_default_config(config)

  ! Read input and change default configurations, if requested
  call config_read_input(config,lupri,luerr)
  ls%input%dalton = config%integral

  doDFT = config%opt%calctype.EQ.config%opt%dftcalc

  ! Initialize lsitem structure
  call ls_init(ls,lupri,luerr,nbast,config%integral,doDFT,config%doDEC,.true.)
  call set_final_config_and_print(lupri,config,ls,nbast)
  
  ! eventually empirical dispersion correction in case of dft
  CALL II_DFTDISP(LS%SETTING,DUMMY,1,1,0,LUPRI,1)

  call mat_pass_info(LUPRI,config%opt%info_matop,mem_monitor)

  CALL LSTIMER('*INPUT',TIMSTR,TIMEND,lupri)

#if defined(VAR_INT64)
  if (bit_size(nbast)==64) then
     write(lupri,*) 'Using 64-bit integers!'
  endif
#endif

 !Initialize grand canonical electronic configuration data table
 call EcData_init 

 ! Grand-canonical (GC) basis?
 if (config%decomp%cfg_gcbasis) then
    call trilevel_basis(config%opt,ls)
    CALL LSTIMER('*ATOM ',TIMSTR,TIMEND,lupri)
 endif
  
  if (config%opt%cfg_prefer_bsm) then
#ifdef HAVE_BSM
     CALL ls_initbsm(ls%input%BASIS%REGULAR,ls%input)
#else
     call lsquit('.BLOCK requested but BSM is not compiled',lupri)
#endif
  endif

  if (config%sparsetest) then
    call sparsetest(ls%setting, lupri)
  endif

end SUBROUTINE init_lsdalton_and_get_lsitem


!> \brief Initilize lsitem from input using init_lsdalton_and_get_lsitem.
!> Intended to be used for standalone tools.
!> \author Kasper Kristensen 
!> \date December 2011
SUBROUTINE get_lsitem_from_input(ls)

  implicit none
  !> lsitem structure initialized according to information in DALTON.INP and MOLECULE.INP
  TYPE(lsitem),target,intent(inout) :: ls
  type(configItem)    :: config
  integer :: lupri,luerr,nbast
  logical :: mem_monitor

  ! Open LSDALTON files
  call open_lsdalton_files(lupri,luerr)

  ! Init lsitem structure
  call init_lsdalton_and_get_lsitem(lupri,luerr,nbast,ls,config,mem_monitor)

  ! Close files LSDALTON.OUT and LSDALTON.ERR again
  call lsclose(lupri,'KEEP')
  call lsclose(luerr,'KEEP')

end SUBROUTINE get_lsitem_from_input


!> \brief Open LSDALTON.OUT and LSDALTON.ERR files on file units
!> lupri and luerr, respectively.
!> \author Kasper Kristensen 
!> \date January 2012
SUBROUTINE open_lsdalton_files(lupri,luerr)

  implicit none
  !> Logical unit number for LSDALTON.OUT
  integer,intent(inout) :: lupri
  !> Logical unit number for LSDALTON.ERR
  integer,intent(inout) :: luerr

  LUPRI=-1
  LUERR=-1
  CALL LSOPEN(LUPRI,'LSDALTON.OUT','NEW','FORMATTED')
  CALL LSOPEN(LUERR,'LSDALTON.ERR','UNKNOWN','FORMATTED')

end SUBROUTINE open_lsdalton_files


!> \brief Print author list for stand-alone f90 linear scaling SCF to LSDALTON.OUT
!> \author T. Kjaergaard, S. Reine
!> \date 2008-10-26
SUBROUTINE print_intro(lupri)
implicit none
!> Logical unit number for file DALTON.OUT
integer,intent(in)        :: lupri
integer                   :: I

WRITE(LUPRI,*)' '
WRITE(LUPRI,*)'    ******************************************************************    '
WRITE(LUPRI,*)'    **********  LSDALTON - An electronic structure program  **********    '
WRITE(LUPRI,*)'    ******************************************************************    '
WRITE(LUPRI,*)' '


         WRITE (LUPRI,'(5X,A)')&
     &  ' This is output from LSDALTON (Release Dalton2013)',&
     &  '                                           ',&
     &  ' Authors in alphabetical order (main contribution(s) in parenthesis)',&
     &  ' -------------------------------------------------------------------',&
     &  ' Vebjoern Bakken,         University of Oslo,     Norway    (geometry optimizer)',&
     &  ' Sonia Coriani,           University of Trieste,  Italy     (response)',&
     &  ' Patrick Ettenhuber,      Aarhus University,      Denmark   (CCSD)',&
     &  ' Stinne Hoest,            Aarhus University,      Denmark   (SCF optimization)',&
     &  ' Ida-Marie Hoeyvik,       Aarhus University,      Denmark   (orbital localization, SCF optimization)',&
     &  ' Branislav Jansik,        Aarhus University,      Denmark   (trilevel, orbital localization)',&
     &  ' Joanna Kauczor,          Aarhus University,      Denmark   (response solver)',&
     &  ' Thomas Kjaergaard,       Aarhus University,      Denmark   (response, integrals)',&
     &  ' Andreas Krapp,           University of Oslo,     Norway    (FMM, dispersion-corrected DFT)',&
     &  ' Kasper Kristensen,       Aarhus University,      Denmark   (response, DEC)',&
     &  ' Patrick Merlot,          University of Oslo,     Norway    (PARI)',&
     &  ' Simen Reine,             University of Oslo,     Norway    (integrals, geometry optimizer)',&
     &  ' Vladimir Rybkin,         University of Oslo,     Norway    (geometry optimizer, dynamics)',&
     &  ' Pawel Salek,             KTH Stockholm,          Sweden    (FMM, DFT functionals)',&
     &  ' Andreas J. Thorvaldsen,  University of Tromsoe,  Norway    (response)',&
     &  ' Lea Thoegersen,          Aarhus University,      Denmark   (SCF optimization)',&
     &  ' Mark Watson,             University of Oslo,     Norway    (FMM)',&
     &  ' Marcin Ziolkowski,       Aarhus University,      Denmark   (DEC)'

     WRITE(LUPRI,*)
         WRITE (LUPRI,'(5X,A)')&
     &  ' Mentors',&
     &  ' -------',&
     &  ' Trygve Helgaker,        University of Oslo,     Norway ',&
     &  ' Poul Joergensen,        Aarhus University,      Denmark',&
     &  ' Jeppe Olsen,            Aarhus University,      Denmark'

     WRITE(LUPRI,*)
     WRITE(LUPRI,*)
         WRITE (LUPRI,'(5X,A)')&
     &'NOTE:',&
     &' ',&
     &'This is an experimental code for the evaluation of molecular',&
     &'energies and properties using various electronic structure models.',&
     &'The authors accept no responsibility for the performance of the code or',&
     &'for the correctness of the results.',&
     &' ',&
     &'The code (in whole or part) is provided under a licence and',&
     &'is not to be reproduced for further distribution without',&
     &'the written permission of the authors or their representatives.',&
     &' ',&
     &'See the home page "http://daltonprogram.org"',&
     &'for further information.',&
     &' ',&
     &'If results obtained with this code are published,',&
     &'an appropriate citation would be:',&
     &' ',&
     &'"LSDalton, a molecular electronic structure program, Release 2013,',&
     &' written by <INSERT AUTHOR LIST>"'
        write(lupri,*)

END SUBROUTINE

end module init_lsdalton_mod

