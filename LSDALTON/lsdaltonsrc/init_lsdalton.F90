!> @file
!> Contains single subroutine for initializing lsitem and config structures 
!> from information in LSDALTON.INP and MOLECULE.INP

!> \brief Initilize lsitem and config structures from LSDALTON.INP and MOLECULE.INP
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
  use daltoninfo, only: ls_init,ls_free
  use IIDFTD, only: II_DFT_DISP
  use matrix_operations, only: mat_no_of_matmuls, mat_pass_info, no_of_matmuls
  use lsmpi_type, only: lsmpi_finalize, lsmpi_print
  use memory_handling, only: Print_Memory_info
  use tensor_interface_module, only:lspdm_free_global_buffer
  private
  public :: open_lsdalton_files,init_lsdalton_and_get_lsitem, &
       & get_lsitem_from_input,finalize_lsdalton_driver_and_free
contains

!> \brief Initilize lsitem and config structures from LSDALTON.INP and MOLECULE.INP
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
   !> lsitem structure initialized according to information in LSDALTON.INP and MOLECULE.INP
   TYPE(lsitem),target,intent(inout) :: ls
   !> configuration structure initialized according to information in LSDALTON.INP and MOLECULE.INP
   type(configItem),intent(inout)    :: config
   !> Monitor memory - mostly for debugging, set to false by default
   logical, intent(inout) :: mem_monitor
   logical             :: doDFT
   REAL(REALK) :: TIMSTR,TIMEND
   real(realk) :: DUMMY(1,1)

   ! Initializations 
!#ifdef VAR_LSDEBUG
   mem_monitor = .true.  !Mostly for memory debugging
!#else
!   mem_monitor = .false. !Mostly for memory debugging
!#endif
   CALL PRINT_INTRO(LUPRI)
   call lsmpi_print(lupri)
#ifdef BINARY_INFO_AVAILABLE
   call print_binary_info(lupri)
#endif

   ! Timing of individual steps
   CALL LSTIMER('START ',TIMSTR,TIMEND,lupri)

   call mat_no_of_matmuls(no_of_matmuls)

   ! Set default configurations for config structure
   call config_set_default_config(config)

   ! Read input and change default configurations, if requested
   call config_read_input(config,lupri,luerr)
   CALL Print_Memory_info(lupri,'at (almost) the Beginning of LSDALTON')
   ls%input%dalton = config%integral

   doDFT = config%opt%calctype.EQ.config%opt%dftcalc

   ! Initialize lsitem structure
   call ls_init(ls,lupri,luerr,nbast,config%integral,doDFT,config%doDEC,.true.)
   call set_final_config_and_print(lupri,config,ls,nbast)

   ! eventually empirical dispersion correction in case of dft
   CALL II_DFT_DISP(LS%SETTING,DUMMY,1,1,0,LUPRI)

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
      CALL Print_Memory_info(lupri,'before the GCBASIS calc')
      call trilevel_basis(config%opt,ls)
      CALL Print_Memory_info(lupri,'after the GCBASIS calc')
      CALL LSTIMER('*ATOM ',TIMSTR,TIMEND,lupri)
   endif

   if (config%sparsetest) then
      call sparsetest(ls%setting, lupri)
   endif

   call Test_if_64bit_integer_required(nbast,nbast)

end SUBROUTINE init_lsdalton_and_get_lsitem
!> \author Patrick Ettenhuber
!> \brief free all the structures that could only be initialized after the
!configuration was read. The other structures are initialized and freed outside
!of LSALTON_DRIVER
!> \date December 2015
SUBROUTINE finalize_lsdalton_driver_and_free(lupri,luerr,ls,config,meminfo_slaves)

   implicit none
   !> Logical unit number for LSDALTON.OUT
   integer,intent(in) :: lupri
   !> Logical unit number for LSDALTON.ERR
   integer,intent(in) :: luerr
   !> lsitem structure initialized according to information in LSDALTON.INP and MOLECULE.INP
   TYPE(lsitem),target,intent(inout) :: ls
   !> configuration structure initialized according to information in LSDALTON.INP and MOLECULE.INP
   type(configItem),intent(inout)    :: config
   logical, intent(inout) :: meminfo_slaves
   call ls_free(ls)

   if(config%opt%cfg_prefer_PDMM)then
      ! Free the background buffer used with PDMM
      call lspdm_free_global_buffer(.true.)
   endif

   meminfo_slaves = config%mpi_mem_monitor
   call Print_Memory_info(lupri,'End of LSDALTON_DRIVER')

end SUBROUTINE finalize_lsdalton_driver_and_free


!> \brief Initilize lsitem from input using init_lsdalton_and_get_lsitem.
!> Intended to be used for standalone tools.
!> \author Kasper Kristensen 
!> \date December 2011
SUBROUTINE get_lsitem_from_input(ls)

  implicit none
  !> lsitem structure initialized according to information in LSDALTON.INP and MOLECULE.INP
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
  !
  logical :: file_exists
  LUPRI=-1
  INQUIRE(FILE='LSDALTON.OUT',EXIST=file_exists)
  IF(file_exists)THEN
     CALL LSOPEN(LUPRI,'LSDALTON.OUT','OLD','FORMATTED')
     call lsclose(LUPRI,'DELETE')
     LUPRI=-1
  ENDIF
  CALL LSOPEN(LUPRI,'LSDALTON.OUT','NEW','FORMATTED')

  LUERR=-1
  INQUIRE(FILE='LSDALTON.ERR',EXIST=file_exists)
  IF(file_exists)THEN
     CALL LSOPEN(LUERR,'LSDALTON.ERR','OLD','FORMATTED')
     call lsclose(LUERR,'DELETE')
     LUERR=-1
  ENDIF
  CALL LSOPEN(LUERR,'LSDALTON.ERR','NEW','FORMATTED')

end SUBROUTINE open_lsdalton_files


!> \brief Print author list for stand-alone f90 linear scaling SCF to LSDALTON.OUT
!> \author T. Kjaergaard, S. Reine
!> \date 2008-10-26
SUBROUTINE print_intro(lupri)
  implicit none
  !> Logical unit number for file DALTON.OUT
  integer,intent(in)        :: lupri
  integer                   :: I

#ifdef BINARY_INFO_AVAILABLE
! sets DALTON_VERSION
#include "dalton_config.h"
#endif
  WRITE(LUPRI,*)' '
  WRITE(LUPRI,*)'    ******************************************************************    '
  WRITE(LUPRI,*)'    **********  LSDALTON - An electronic structure program  **********    '
  WRITE(LUPRI,*)'    ******************************************************************    '
  WRITE(LUPRI,*)' '
  write(LUPRI,*)' '

#ifdef BINARY_INFO_AVAILABLE
  WRITE (LUPRI,'(5X,A)')' This is output from LSDALTON '//DALTON_VERSION
#else
  WRITE (LUPRI,'(5X,A)')' This is output from LSDALTON '
#endif
  write(LUPRI,*)' '
  write(lupri,*)' '

  WRITE(LUPRI,'(5X,A)')&
       &' IF RESULTS OBTAINED WITH THIS CODE ARE PUBLISHED,',&
       &' THE FOLLOWING PAPER SHOULD BE CITED:',&
       &' ',&
       & ' K. Aidas, C. Angeli, K. L. Bak, V. Bakken, R. Bast,',&
       & ' L. Boman, O. Christiansen, R. Cimiraglia, S. Coriani,',&
       & ' P. Dahle, E. K. Dalskov, U. Ekstroem, T. Enevoldsen,',&
       & ' J. J. Eriksen, P. Ettenhuber, B. Fernandez,',&
       & ' L. Ferrighi, H. Fliegl, L. Frediani, K. Hald,',&
       & ' A. Halkier, C. Haettig, H. Heiberg,',&
       & ' T. Helgaker, A. C. Hennum, H. Hettema,',&
       & ' E. Hjertenaes, S. Hoest, I.-M. Hoeyvik,',&
       & ' M. F. Iozzi, B. Jansik, H. J. Aa. Jensen,',&
       & ' D. Jonsson, P. Joergensen, J. Kauczor,',&
       & ' S. Kirpekar, T. Kjaergaard, W. Klopper,',&
       & ' S. Knecht, R. Kobayashi, H. Koch, J. Kongsted,',&
       & ' A. Krapp, K. Kristensen, A. Ligabue,',&
       & ' O. B. Lutnaes, J. I. Melo, K. V. Mikkelsen, R. H. Myhre,',&
       & ' C. Neiss, C. B. Nielsen, P. Norman,',&
       & ' J. Olsen, J. M. H. Olsen, A. Osted,',&
       & ' M. J. Packer, F. Pawlowski, T. B. Pedersen,',&
       & ' P. F. Provasi, S. Reine, Z. Rinkevicius,',&
       & ' T. A. Ruden, K. Ruud, V. Rybkin,',&
       & ' P. Salek, C. C. M. Samson, A. Sanchez de Meras,',&
       & ' T. Saue, S. P. A. Sauer, B. Schimmelpfennig,',&
       & ' K. Sneskov, A. H. Steindal, K. O. Sylvester-Hvid,',&
       & ' P. R. Taylor, A. M. Teale, E. I. Tellgren,',&
       & ' D. P. Tew, A. J. Thorvaldsen, L. Thoegersen,',&
       & ' O. Vahtras, M. A. Watson, D. J. D. Wilson,',&
       & ' M. Ziolkowski, and H. AAgren,',&
       & ' "The Dalton quantum chemistry program system",',&
       & ' WIREs Comput. Mol. Sci. (doi: 10.1002/wcms.1172)'
  write(lupri,*)' '
  write(lupri,*)' '


  WRITE(LUPRI,'(5X,A)')&
  &  '                                           ',&
       &  ' LSDALTON authors in alphabetical order (main contribution(s) in parenthesis)',&
       &  ' ----------------------------------------------------------------------------',&
       &  ' Vebjoern Bakken,         University of Oslo,        Norway    (geometry optimizer)',&
       &  ' Radovan Bast,            KTH Stockholm,             Sweden    (CMake, Testing)',&
       &  ' Sonia Coriani,           University of Trieste,     Italy     (response)',&
       &  ' Patrick Ettenhuber,      Aarhus University,         Denmark   (CCSD)',&
       &  ' Trygve Helgaker,         University of Oslo,        Norway    (supervision)',&
       &  ' Stinne Hoest,            Aarhus University,         Denmark   (SCF optimization)',&
       &  ' Ida-Marie Hoeyvik,       Aarhus University,         Denmark   (orbital localization, SCF optimization)',&
       &  ' Robert Izsak,            University of Oslo,        Norway    (ADMM)',&
       &  ' Branislav Jansik,        Aarhus University,         Denmark   (trilevel, orbital localization)',&
       &  ' Poul Joergensen,         Aarhus University,         Denmark   (supervision)', &
       &  ' Joanna Kauczor,          Aarhus University,         Denmark   (response solver)',&
       &  ' Thomas Kjaergaard,       Aarhus University,         Denmark   (response, integrals)',&
       &  ' Andreas Krapp,           University of Oslo,        Norway    (FMM, dispersion-corrected DFT)',&
       &  ' Kasper Kristensen,       Aarhus University,         Denmark   (response, DEC)',&
       &  ' Patrick Merlot,          University of Oslo,        Norway    (ADMM)',&
       &  ' Cecilie Nygaard,         Aarhus University,         Denmark   (SOEO)',&
       &  ' Jeppe Olsen,             Aarhus University,         Denmark   (supervision)', &
       &  ' Simen Reine,             University of Oslo,        Norway    (integrals, geometry optimizer)',&
       &  ' Vladimir Rybkin,         University of Oslo,        Norway    (geometry optimizer, dynamics)',&
       &  ' Pawel Salek,             KTH Stockholm,             Sweden    (FMM, DFT functionals)',&
       &  ' Andrew M. Teale,         University of Nottingham   England   (E-coefficients)',&
       &  ' Erik Tellgren,           University of Oslo,        Norway    (density fitting, E-coefficients)',&
       &  ' Andreas J. Thorvaldsen,  University of Tromsoe,     Norway    (response)',&
       &  ' Lea Thoegersen,          Aarhus University,         Denmark   (SCF optimization)',&
       &  ' Mark Watson,             University of Oslo,        Norway    (FMM)',&
       &  ' Marcin Ziolkowski,       Aarhus University,         Denmark   (DEC)'

  WRITE(LUPRI,*)' '
  WRITE(LUPRI,*)' '
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
       &'for further information.'
  WRITE(LUPRI,*)' '
  WRITE(LUPRI,*)' '


END SUBROUTINE print_intro

end module init_lsdalton_mod

