!dalton_copyright_start
!
!
!dalton_copyright_end

module lucita_mcscf_ci_interface_procedures

! stefan: this module contains all interface procedures called before/after 
!         a particular MCSCF-CI run task

  use lucita_setup, only: setup_lucita_mcci_wrkspc_dimensions
  use lucita_mcscf_vector_exchange, only : vector_exchange_driver
  use lucita_mcscf_ci_cfg
  use file_type_module, only : file_info, file_type, file_type_active_fh_reset_lucipar

  implicit none

  public mcscf_lucita_interface

  private

  integer, public    :: restore_cref_vector_switch = -1
  logical, public    :: restore_cref               = .false.
  logical, public    :: accumulate_hc              = .false.
  integer, parameter, public :: print_lvl_limit = 1000

contains

  subroutine mcscf_lucita_interface(c_or_cr,          &
                                    hc_or_cl,         &
                                    int1_or_rho1,     &
                                    int2_or_rho2,     &
                                    run_type,         &
                                    work,             &
                                    lwork,            &
                                    print_lvl)
!*******************************************************************************
!
!    purpose: provide a generic interface to LUCITA for various MCSCF-CI tasks.
!              
!
!*******************************************************************************
#include "priunit.h"
      integer,   intent(in)    :: lwork
      integer,   intent(in)    :: print_lvl
      real(8),   intent(inout) :: work(lwork)
! hjaaj Apr 2016: intent(inout) does not work with paramter vdummy in parameter list
!     real(8),   intent(inout) :: c_or_cr(*)
!     real(8),   intent(inout) :: hc_or_cl(*)
!     real(8),   intent(inout) :: int1_or_rho1(*)
!     real(8),   intent(inout) :: int2_or_rho2(*)
      real(8)                  :: c_or_cr(*)
      real(8)                  :: hc_or_cl(*)
      real(8)                  :: int1_or_rho1(*)
      real(8)                  :: int2_or_rho2(*)
      character, intent(in)    :: run_type*(*)
!
!-------------------------------------------------------------------------------
      integer                  :: lfree
!-------------------------------------------------------------------------------

      lfree = lwork

!     a. pre-process
      call mcscf_pre_lucita_processing(c_or_cr,             &
                                       hc_or_cl,            &
                                       int1_or_rho1,        &
                                       int2_or_rho2,        &
                                       run_type,            &
                                       work,                &
                                       lfree,               &
                                       print_lvl)

!     b. run LUCITA
      call lucita(run_type,                                 &
                  c_or_cr,                                  &
                  hc_or_cl,                                 &
                  int1_or_rho1,                             &
                  int2_or_rho2,                             &
                  work,                                     &
                  lfree)

!     c. post-process
      call mcscf_post_lucita_processing(c_or_cr,            &
                                        hc_or_cl,           &
                                        int1_or_rho1,       &
                                        int2_or_rho2,       &
                                        run_type,           &
                                        work,               &
                                        lfree,              &
                                        print_lvl)

  end subroutine mcscf_lucita_interface
!*******************************************************************************

  subroutine mcscf_pre_lucita_processing(c_or_cr,          &
                                         hc_or_cl,         &
                                         int1_or_rho1,     &
                                         int2_or_rho2,     &
                                         run_type,         &
                                         work,             &
                                         lwork,            &
                                         print_lvl)
!*******************************************************************************
!
!    purpose: provide a generic interface to pre-LUCITA processing for various 
!             MCSCF-CI tasks.
!              
!
!*******************************************************************************
#include "priunit.h"
      integer,   intent(inout) :: lwork
      integer,   intent(in)    :: print_lvl
      real(8),   intent(inout) :: work(lwork)
! hjaaj Apr 2016: intent(inout) does not work with paramter vdummy in parameter list
!     real(8),   intent(inout) :: c_or_cr(*)
!     real(8),   intent(inout) :: hc_or_cl(*)
!     real(8),   intent(inout) :: int1_or_rho1(*)
!     real(8),   intent(inout) :: int2_or_rho2(*)
      real(8)                  :: c_or_cr(*)
      real(8)                  :: hc_or_cl(*)
      real(8)                  :: int1_or_rho1(*)
      real(8)                  :: int2_or_rho2(*)
      character, intent(in)    :: run_type*(*)
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

!     a. define allocations for a given CI runtype needed by slaves in parallel runs
      call setup_lucita_mcci_wrkspc_dimensions(run_type,         &
                                               print_lvl)
!     b. initialize file handles
      call file_type_active_fh_reset_lucipar(file_info)

!     c. select pre-process task
      select case(run_type)
        case('return CIdia', 'ijkl resort ', 'fci dump    ')
!         do nothing
        case('srdft   ci  ')
          call mcscf_pre_lucita_srdft(print_lvl)
        case('return CIdim')
          call mcscf_pre_lucita_setci(print_lvl)
        case('initial ci  ')
          call mcscf_pre_lucita_cistart(print_lvl)
        case('Xp-density m')
          call mcscf_pre_lucita_xpdens(c_or_cr,hc_or_cl,print_lvl)
        case('sigma vec   ')
          call mcscf_pre_lucita_return_e2b(c_or_cr,hc_or_cl,print_lvl)
        case('analyze Cvec')
          call mcscf_pre_lucita_bvec_analyze(c_or_cr,print_lvl)
        case('rotate  Cvec')
          call mcscf_pre_lucita_rotate_cref(c_or_cr,print_lvl)
        case default
          call quit('undefined pre LUCITA processing step in MCSCF process flow')
      end select

  end subroutine mcscf_pre_lucita_processing
!*******************************************************************************

  subroutine mcscf_post_lucita_processing(c_or_cr,          &
                                          hc_or_cl,         &
                                          int1_or_rho1,     &
                                          int2_or_rho2,     &
                                          run_type,         &
                                          work,             &
                                          lwork,            &
                                          print_lvl)
!*******************************************************************************
!
!    purpose: provide a generic interface to post-LUCITA processing for various 
!             MCSCF-CI tasks.
!              
!
!*******************************************************************************
#include "priunit.h"
      integer,   intent(inout) :: lwork
      integer,   intent(in)    :: print_lvl
      real(8),   intent(inout) :: work(lwork)
! hjaaj Apr 2016: intent(inout) does not work with paramter vdummy in parameter list
!     real(8),   intent(inout) :: c_or_cr(*)
!     real(8),   intent(inout) :: hc_or_cl(*)
!     real(8),   intent(inout) :: int1_or_rho1(*)
!     real(8),   intent(inout) :: int2_or_rho2(*)
      real(8)                  :: c_or_cr(*)
      real(8)                  :: hc_or_cl(*)
      real(8)                  :: int1_or_rho1(*)
      real(8)                  :: int2_or_rho2(*)
      character, intent(in)    :: run_type*(*)
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

!     a. initialize file handles
      call file_type_active_fh_reset_lucipar(file_info)

!     b. choose post-process task
      select case(run_type)
        case('analyze Cvec', 'ijkl resort ', 'fci dump    ')
!         do nothing
        case('srdft   ci  ')
          call mcscf_post_lucita_srdft(c_or_cr,print_lvl)
        case('return CIdim')
          call mcscf_post_lucita_setci(print_lvl)
        case('initial ci  ')
          call mcscf_post_lucita_cistart(c_or_cr,print_lvl)
        case('return CIdia')
          call mcscf_post_lucita_hdiag(c_or_cr,print_lvl)
        case('sigma vec   ')
          call mcscf_post_lucita_return_e2b(hc_or_cl,c_or_cr,print_lvl)
        case('rotate  Cvec')
          call mcscf_post_lucita_rotate_cref(c_or_cr,print_lvl)
        case('Xp-density m')
          call mcscf_post_lucita_xpdens(hc_or_cl,c_or_cr,print_lvl)
        case default
          call quit('undefined post LUCITA processing step in MCSCF process flow')
      end select

!     write(lupri,*) ' restore_cref, restore_cref_vector_switch', &
!                      restore_cref, restore_cref_vector_switch
!     c. restore cref
      select case(restore_cref)
        case(.true.)
          call restore_cref_lucita(restore_cref,restore_cref_vector_switch,c_or_cr,hc_or_cl)
        case(.false.)
!         nothing to do          
      end select

  end subroutine mcscf_post_lucita_processing
!*******************************************************************************

  subroutine restore_cref_lucita(do_restore,switch,c_or_cr,hc_or_cl)
!*******************************************************************************
!
!    purpose: restore CREF in parallel runs if necessary
!              
!
!*******************************************************************************
#include "priunit.h"
! hjaaj Apr 2016: intent(inout) does not work with paramter vdummy in parameter list
!     real(8),   intent(inout) :: c_or_cr(*)
!     real(8),   intent(inout) :: hc_or_cl(*)
      real(8)                  :: c_or_cr(*)
      real(8)                  :: hc_or_cl(*)
      integer,   intent(inout) :: switch
      logical,   intent(inout) :: do_restore
!-------------------------------------------------------------------------------

!     write(lupri,*) 'restore cref in action...',switch
      select case(switch)
        case(1)
          call rdcref(c_or_cr ,.false.)
        case(2)
          call rdcref(hc_or_cl,.false.)
      end select
      do_restore = .false.
      switch     = -1
      
  end subroutine restore_cref_lucita
!*******************************************************************************

  subroutine mcscf_pre_lucita_bvec_analyze(c_or_cr,print_lvl)
!*******************************************************************************
!
!    purpose: pre-LUCITA processing for CI task: analyze Cvec
!              
!
!*******************************************************************************
! lucita
  use lucita_cfg
! sirius
! nothing
#include "priunit.h"
! hjaaj Apr 2016: intent(inout) does not work with paramter vdummy in parameter list
!     real(8),   intent(inout) :: c_or_cr(*)
      real(8)                  :: c_or_cr(*)
      integer,   intent(in)    :: print_lvl
!-------------------------------------------------------------------------------
      integer                  :: push_pull = 2
!-------------------------------------------------------------------------------
      
#ifdef VAR_MPI
!     stefan: always enforce vector_update in order to write the full vector to the sequential file on the master
      vector_update_mc2lu_lu2mc((2-1)*vector_exchange_types+vector_exchange_type1) = .true.
#endif

!     rhs vector (the first argument '2' refers to the direction of transfer: mcscf ==> lucita)
      call vector_exchange_driver(push_pull,vector_exchange_type1,lucita_cfg_nr_roots,lucita_cfg_csym,                &
                                  io2io_vector_exchange_mc2lu_lu2mc,                                                  &
                                  vector_update_mc2lu_lu2mc((2-1)*vector_exchange_types+vector_exchange_type1),       &
                                  .false.,c_or_cr)

      if(print_lvl >= print_lvl_limit)then
        write(lupri,'(/a)') ' *** LUCITA-MCSCF interface reports: ***'
        write(lupri,*) ' rhs vector saved on lucita files'
        write(lupri,'(a/)') ' *** end of LUCITA-MCSCF interface   ***'
      end if

  end subroutine mcscf_pre_lucita_bvec_analyze
!*******************************************************************************

  subroutine mcscf_pre_lucita_srdft(print_lvl)
!*******************************************************************************
!
!    purpose: pre-LUCITA processing for CI task: srdft ci  
!              
!
!*******************************************************************************
! lucita
! use lucita_mcscf_ci_cfg ! globally included
! sirius
#include "cicb01.h"
#include "maxorb.h"
#include "infinp.h"
#include "priunit.h"
      integer,   intent(in)    :: print_lvl
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      
      einact_mc2lu           = emy_ci ! inactive energy
      docisrdft_mc2lu        = docisrdft                                     

      if(print_lvl >= print_lvl_limit)then
        write(lupri,'(/a)') ' *** LUCITA-MCSCF interface reports: ***'
        write(lupri,*) ' inactive energy obtained from MCSCF environment ==> ', einact_mc2lu
        write(lupri,'(a/)') ' *** end of LUCITA-MCSCF interface   ***'
      end if

  end subroutine mcscf_pre_lucita_srdft
!*******************************************************************************

  subroutine mcscf_pre_lucita_cistart(print_lvl)
!*******************************************************************************
!
!    purpose: pre-LUCITA processing for CI task: initial ci  
!              
!
!*******************************************************************************
! lucita
! use lucita_mcscf_ci_cfg ! globally included
! sirius
#include "cicb01.h"
#include "priunit.h"
      integer,   intent(in)    :: print_lvl
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      
      einact_mc2lu           = emy_ci ! inactive energy

      if(print_lvl >= print_lvl_limit)then
        write(lupri,'(/a)') ' *** LUCITA-MCSCF interface reports: ***'
        write(lupri,*) ' inactive energy obtained from MCSCF environment ==> ', einact_mc2lu
        write(lupri,'(a/)') ' *** end of LUCITA-MCSCF interface   ***'
      end if

  end subroutine mcscf_pre_lucita_cistart
!*******************************************************************************

  subroutine mcscf_pre_lucita_return_e2b(c_or_cr,hc_or_cl,print_lvl)
!*******************************************************************************
!
!    purpose: pre-LUCITA processing for CI task: sigma vec   
!              
!
!*******************************************************************************
! lucita
  use lucita_cfg
! sirius
#include "ciinfo.h"
! nothing
#include "priunit.h"
! hjaaj Apr 2016: intent(inout) does not work with paramter vdummy in parameter list
!     real(8),   intent(inout) :: c_or_cr(*)
!     real(8),   intent(inout) :: hc_or_cl(*)
      real(8)                  :: c_or_cr(*)
      real(8)                  :: hc_or_cl(*)
      integer,   intent(in)    :: print_lvl
!-------------------------------------------------------------------------------
      integer                  :: push_pull = 2
      real(8)                  :: norm_hc
      real(8), external        :: ddot
!-------------------------------------------------------------------------------
      
!     potentially restore cref after the ci task
!     write(lupri,*) ' check for restorage, vector_exchange_type1',vector_exchange_type1
      if(vector_exchange_type1 == 1)then
        restore_cref_vector_switch = 1
        restore_cref               = .true.
!       write(lupri,*) ' restore cref activated'
      end if

!     rhs vector (the first argument '2' refers to the direction of transfer: mcscf ==> lucita)
      call vector_exchange_driver(push_pull,vector_exchange_type1,lucita_cfg_nr_roots,lucita_cfg_csym,                &
                                  io2io_vector_exchange_mc2lu_lu2mc,                                                  &
                                  vector_update_mc2lu_lu2mc((2-1)*vector_exchange_types+vector_exchange_type1),       &
                                  .true.,c_or_cr)

!     check for e2b vector accumulation
      norm_hc = 0.0d0
      norm_hc = ddot(ndtasm(lucita_cfg_hcsym),hc_or_cl,1,hc_or_cl,1)
      if(norm_hc > 0.0d0)then
        accumulate_hc         = .true.
        vector_exchange_type1 = 4
        vector_update_mc2lu_lu2mc((push_pull-1)*vector_exchange_types+vector_exchange_type1) = .true.

!       (the first argument '2' refers to the direction of transfer: mcscf ==> lucita)
        call vector_exchange_driver(push_pull,vector_exchange_type1,lucita_cfg_nr_roots,lucita_cfg_hcsym,             &
                                    io2io_vector_exchange_mc2lu_lu2mc,                                                &
                                    vector_update_mc2lu_lu2mc((2-1)*vector_exchange_types+vector_exchange_type1),     &
                                    .true.,hc_or_cl)

      end if

      if(print_lvl >= print_lvl_limit)then
        write(lupri,'(/a)') ' *** LUCITA-MCSCF interface reports: ***'
        write(lupri,*) ' rhs vector saved on lucita files'
        write(lupri,'(a/)') ' *** end of LUCITA-MCSCF interface   ***'
      end if

  end subroutine mcscf_pre_lucita_return_e2b
!*******************************************************************************

  subroutine mcscf_pre_lucita_rotate_cref(c_or_cr,print_lvl)
!*******************************************************************************
!
!    purpose: pre-LUCITA processing for CI task: rotate  Cvec
!              
!
!*******************************************************************************
! lucita
  use lucita_cfg
! sirius
! nothing
#include "priunit.h"
! hjaaj Apr 2016: intent(inout) does not work with paramter vdummy in parameter list
!     real(8),   intent(inout) :: c_or_cr(*)
      real(8)                  :: c_or_cr(*)
      integer,   intent(in)    :: print_lvl
!-------------------------------------------------------------------------------
      integer                  :: push_pull = 2
!-------------------------------------------------------------------------------
      
!     rhs vector (the first argument '2' refers to the direction of transfer: mcscf ==> lucita)
      call vector_exchange_driver(push_pull,vector_exchange_type1,lucita_cfg_nr_roots,lucita_cfg_csym,                &
                                  io2io_vector_exchange_mc2lu_lu2mc,                                                  &
                                  vector_update_mc2lu_lu2mc((2-1)*vector_exchange_types+vector_exchange_type1),       &
                                  .true.,c_or_cr)

      if(print_lvl >= print_lvl_limit)then
        write(lupri,'(/a)') ' *** LUCITA-MCSCF interface reports: ***'
        write(lupri,*) ' rhs vector saved on lucita files'
        write(lupri,'(a/)') ' *** end of LUCITA-MCSCF interface   ***'
      end if

  end subroutine mcscf_pre_lucita_rotate_cref
!*******************************************************************************

  subroutine mcscf_pre_lucita_setci(print_lvl)
!*******************************************************************************
!
!    purpose: pre-LUCITA processing for CI task: return CIdim 
!              
!
!*******************************************************************************
! lucita
  use lucita_cfg
  use parallel_task_distribution_type_module
  use parallel_models_mpi, only: lucita_models_enabled
#include "priunit.h"
! sirius
      integer,   intent(in)    :: print_lvl
!
!-------------------------------------------------------------------------------

!     turn on common block (re-)initialization + co-worker synchronization of static lucita variables 
      lucita_cfg_initialize_cb = .true.
!     free parallel distribution object and thereby turn on (re-)calculation of static task list distribution
      call parallel_task_distribution_free_lucipar(ptask_distribution)
!     keep track of enabled parallel models in lucita
      lucita_models_enabled = .true.
      
  end subroutine mcscf_pre_lucita_setci
!*******************************************************************************

  subroutine mcscf_pre_lucita_xpdens(c_or_cr,hc_or_cl,print_lvl)
!*******************************************************************************
!
!    purpose: pre-LUCITA processing for CI task: Xp-density m
!              
!
!*******************************************************************************
! lucita
  use lucita_cfg
! sirius
! nothing
#include "priunit.h"
! hjaaj Apr 2016: intent(inout) does not work with paramter vdummy in parameter list
!     real(8),   intent(inout) :: c_or_cr(*)
!     real(8),   intent(inout) :: hc_or_cl(*)
      real(8)                  :: c_or_cr(*)
      real(8)                  :: hc_or_cl(*)
      integer,   intent(in)    :: print_lvl
!-------------------------------------------------------------------------------
      integer                  :: push_pull = 2
!-------------------------------------------------------------------------------

!     write(lupri,*) ' check for restorage, vector_exchange_type1/2',vector_exchange_type1,vector_exchange_type2
!     potentially restore cref after the ci task in parallel runs
!#ifdef MOD_SRDFT
!      if(srdft_ci_1pdens_cref_restore)then
!#endif
        if(vector_exchange_type1 == 1 .or. vector_exchange_type2 == 1)then
          restore_cref               = .true.
          restore_cref_vector_switch = 1
          if(vector_exchange_type2 == 1) restore_cref_vector_switch = 2
!         write(lupri,*) ' xpdens: restore_cref ==>', restore_cref
        end if
!#ifdef MOD_SRDFT
!      end if
!#endif
      
!     rhs vector (the first argument '2' refers to the direction of transfer: mcscf ==> lucita)
      call vector_exchange_driver(push_pull,vector_exchange_type1,lucita_cfg_nr_roots,lucita_cfg_csym,                &
                                  io2io_vector_exchange_mc2lu_lu2mc,                                                  &
                                  vector_update_mc2lu_lu2mc((2-1)*vector_exchange_types+vector_exchange_type1),       &
                                  .true.,c_or_cr)

!     lhs vector (if transition density matrix) - otherwise we take care of it internally in lucita
      call vector_exchange_driver(push_pull,vector_exchange_type2,lucita_cfg_nr_roots,lucita_cfg_hcsym,               &
                                  io2io_vector_exchange_mc2lu_lu2mc,                                                  &
                                  vector_update_mc2lu_lu2mc((2-1)*vector_exchange_types+vector_exchange_type2),       &
                                  .true.,hc_or_cl)


      if(print_lvl >= print_lvl_limit)then
        write(lupri,'(/a)') ' *** LUCITA-MCSCF interface reports: ***'
        write(lupri,*) ' rhs/lhs vectors saved on lucita files'
        write(lupri,'(a/)') ' *** end of LUCITA-MCSCF interface   ***'
      end if

  end subroutine mcscf_pre_lucita_xpdens
!*******************************************************************************

#ifdef MOD_SRDFT
  subroutine mcscf_post_lucita_srdft(c_or_cr,print_lvl)
!*******************************************************************************
!
!    purpose: post-LUCITA processing for CI task: srdft ci  
!              
!
!*******************************************************************************
! lucita
  use lucita_mcscf_srdftci_cfg
#include "cstate.inc"
! sirius
#include "cicb01.h"
#include "priunit.h"
#include "dfterg.h"
      real(8),   intent(out)   :: c_or_cr(*)
      integer,   intent(in)    :: print_lvl
!
!-------------------------------------------------------------------------------
      integer                  :: i
      integer                  :: push_pull = 1
!-------------------------------------------------------------------------------
      
      call dcopy(nroot,eroot,1,energy_root,1)
      call dcopy(nroot,root_residual,1,residual_croot,1)
      call icopy(nroot,root_converged,1,jcroot,1)
      
      jconv_c = 0

      emy_ci = einact_mc2lu + emydft_mc2lu ! inactive energy
      do i = 1, nroot
        if(jcroot(i) > 0) jconv_c = jconv_c + 1
        energy_root(i)            = energy_root(i) - emy_ci
      end do

      ncired = nfinal_vec


!     type 2 - put CI start vector/ vectors on file LUCITA_CVECS.x where x refers to the present symmetry id (a,b,c,...)
!              to file luit3 if io2io_vector_exchange_mc2lu_lu2mc == .true. 
      vector_exchange_type1                                                        = 2
      vector_update_mc2lu_lu2mc((push_pull-1)*vector_exchange_types+vector_exchange_type1) = .true.

!     (the first argument '1' refers to the direction of transfer: lucita ==> mcscf)
      call vector_exchange_driver(push_pull,vector_exchange_type1,nroot,irefsm_c,                                     &
                                  io2io_vector_exchange_mc2lu_lu2mc,                                                  &
                                  vector_update_mc2lu_lu2mc((1-1)*vector_exchange_types+vector_exchange_type1),       &
                                  .false.,c_or_cr)

      if(print_lvl >= print_lvl_limit)then
        write(lupri,'(/a)') ' *** LUCITA-MCSCF interface reports: ***'
        write(lupri,*) ' number of non-converged roots ==> ', jconv_c
        do i = 1, nroot
          write(lupri,*) ' energies of eigenstates(orig) ==> ', eroot(i)
          write(lupri,*) ' energies of eigenstates       ==> ', energy_root(i)
          write(lupri,*) ' energies of eigenstates+shift ==> ', energy_root(i) + einact_mc2lu
          write(lupri,*) ' residual of eigenstates       ==> ', residual_croot(i)
        end do
        write(lupri,'(a/)') ' *** end of LUCITA-MCSCF interface   ***'
      end if

  end subroutine mcscf_post_lucita_srdft
!*******************************************************************************
#endif

  subroutine mcscf_post_lucita_cistart(c_or_cr,print_lvl)
!*******************************************************************************
!
!    purpose: post-LUCITA processing for CI task: initial ci  
!              
!
!*******************************************************************************
! lucita
#include "cstate.inc"
! sirius
#include "cicb01.h"
#include "priunit.h"
      real(8),   intent(out)   :: c_or_cr(*)
      integer,   intent(in)    :: print_lvl
!
!-------------------------------------------------------------------------------
      integer                  :: i
      integer                  :: push_pull = 1
!-------------------------------------------------------------------------------
      
      call dcopy(nroot,eroot,1,energy_root,1)
      call dcopy(nroot,root_residual,1,residual_croot,1)
      call icopy(nroot,root_converged,1,jcroot,1)
      
      jconv_c = 0

      do i = 1, nroot
        if(jcroot(i) > 0) jconv_c = jconv_c + 1
        energy_root(i)            = energy_root(i) - einact_mc2lu
      end do

      ncired = nfinal_vec

!     type 2 - put CI start vector/ vectors on file LUCITA_CVECS.x where x refers to the present symmetry id (a,b,c,...)
!              to file luit3 if io2io_vector_exchange_mc2lu_lu2mc == .true. 
      vector_exchange_type1                                                        = 2
      vector_update_mc2lu_lu2mc((push_pull-1)*vector_exchange_types+vector_exchange_type1) = .true.

!     (the first argument '1' refers to the direction of transfer: lucita ==> mcscf)
      call vector_exchange_driver(push_pull,vector_exchange_type1,nroot,irefsm_c,                                     &
                                  io2io_vector_exchange_mc2lu_lu2mc,                                                  &
                                  vector_update_mc2lu_lu2mc((1-1)*vector_exchange_types+vector_exchange_type1),       &
                                  .false.,c_or_cr)

      if(print_lvl >= print_lvl_limit)then
        write(lupri,'(/a)') ' *** LUCITA-MCSCF interface reports: ***'
        write(lupri,*) ' number of non-converged roots ==> ', jconv_c
        do i = 1, nroot
          write(lupri,*) ' energies of eigenstates(orig) ==> ', eroot(i)
          write(lupri,*) ' energies of eigenstates       ==> ', energy_root(i)
          write(lupri,*) ' energies of eigenstates+shift ==> ', energy_root(i) + einact_mc2lu
          write(lupri,*) ' residual of eigenstates       ==> ', residual_croot(i)
        end do
        write(lupri,'(a/)') ' *** end of LUCITA-MCSCF interface   ***'
      end if

  end subroutine mcscf_post_lucita_cistart
!*******************************************************************************

  subroutine mcscf_post_lucita_hdiag(c_or_cr,print_lvl)
!*******************************************************************************
!
!    purpose: post-LUCITA processing for CI task: return CIdia
!              
!
!*******************************************************************************
! lucita
  use lucita_cfg
! sirius
! nothing
#include "priunit.h"
      real(8),   intent(out)   :: c_or_cr(*)
      integer,   intent(in)    :: print_lvl
!
!-------------------------------------------------------------------------------
      integer                  :: i
      integer                  :: push_pull = 1
!-------------------------------------------------------------------------------
      
!     type 4 - pull H diag from lucita file LUCITA_HDIAG.x where x refers to the present symmetry id (a,b,c,...)
!              to mcscf core memory 
      vector_exchange_type1 = 4
      vector_update_mc2lu_lu2mc((push_pull-1)*vector_exchange_types+vector_exchange_type1) = .true.

!     (the first argument '1' refers to the direction of transfer: lucita ==> mcscf)
      call vector_exchange_driver(push_pull,vector_exchange_type1,lucita_cfg_nr_roots,lucita_cfg_csym,                &
                                  io2io_vector_exchange_mc2lu_lu2mc,                                                  &
                                  vector_update_mc2lu_lu2mc((1-1)*vector_exchange_types+vector_exchange_type1),       &
                                  .true.,c_or_cr)

      if(print_lvl >= print_lvl_limit)then
        write(lupri,'(/a)') ' *** LUCITA-MCSCF interface reports: ***'
        write(lupri,*) ' H diag pushed into mcscf core memory '
        write(lupri,'(a/)') ' *** end of LUCITA-MCSCF interface   ***'
      end if

  end subroutine mcscf_post_lucita_hdiag
!*******************************************************************************

  subroutine mcscf_post_lucita_return_e2b(hc_or_cl,c_or_cr,print_lvl)
!*******************************************************************************
!
!    purpose: post-LUCITA processing for CI task: sigma vec   
!              
!
!*******************************************************************************
! lucita
  use lucita_cfg
! sirius
#include "ciinfo.h"
! nothing
#include "priunit.h"
! hjaaj Apr 2016: intent(inout) does not work with paramter vdummy in parameter list
!     real(8),   intent(inout) :: hc_or_cl(*)
!     real(8),   intent(inout) :: c_or_cr(*)
      real(8)                  :: hc_or_cl(*)
      real(8)                  :: c_or_cr(*)
      integer,   intent(in)    :: print_lvl
!-------------------------------------------------------------------------------
      integer                  :: push_pull = 1
!-------------------------------------------------------------------------------
      
      vector_exchange_type1                      = 3
      vector_update_mc2lu_lu2mc((push_pull-1)*vector_exchange_types+vector_exchange_type1) = .true. ! pull the outgoing lhs vector from LUCITA files to mc core-memory

!     lhs vector (the first argument '1' refers to the direction of transfer: lucita ==> mcscf)
      call vector_exchange_driver(push_pull,vector_exchange_type1,lucita_cfg_nr_roots,lucita_cfg_hcsym,                &
                                  io2io_vector_exchange_mc2lu_lu2mc,                                                   &
                                  vector_update_mc2lu_lu2mc((1-1)*vector_exchange_types+vector_exchange_type1),        &
                                  .true.,hc_or_cl)

      if(accumulate_hc)then
        vector_exchange_type1                    = 4
        vector_update_mc2lu_lu2mc((push_pull-1)*vector_exchange_types+vector_exchange_type1) = .true. ! pull the outgoing lhs vector from LUCITA files to mc core-memory

!       lhs vector (the first argument '1' refers to the direction of transfer: lucita ==> mcscf)
        call vector_exchange_driver(push_pull,vector_exchange_type1,lucita_cfg_nr_roots,lucita_cfg_hcsym,              &
                                    io2io_vector_exchange_mc2lu_lu2mc,                                                 &
                                    vector_update_mc2lu_lu2mc((1-1)*vector_exchange_types+vector_exchange_type1),      &
                                    .true.,c_or_cr)

!       write(lupri,*) 'ndtasm(lucita_cfg_hcsym) ',ndtasm(lucita_cfg_hcsym)
        call daxpy(ndtasm(lucita_cfg_hcsym),1.0d0,c_or_cr,1,hc_or_cl,1)
        accumulate_hc = .false.
!       write(lupri,*) 'accumulated hc....'
      end if

      if(.not.restore_cref)then ! bvec == ci trial vector: need to restore for reduced matrix calculation
!       write(lupri,*) ' restore lhs activated',vector_exchange_type1
        vector_exchange_type1                      = 2
        vector_update_mc2lu_lu2mc((push_pull-1)*vector_exchange_types+vector_exchange_type1) = .true. ! pull the outgoing lhs vector from LUCITA files to mc core-memory

!       rhs vector (the first argument '1' refers to the direction of transfer: lucita ==> mcscf)
        call vector_exchange_driver(push_pull,vector_exchange_type1,lucita_cfg_nr_roots,lucita_cfg_csym,                 &
                                    io2io_vector_exchange_mc2lu_lu2mc,                                                   &
                                    vector_update_mc2lu_lu2mc((1-1)*vector_exchange_types+vector_exchange_type1),        &
                                    .true.,c_or_cr)
      end if

      if(print_lvl >= print_lvl_limit)then
        write(lupri,'(/a)') ' *** LUCITA-MCSCF interface reports: ***'
        write(lupri,*) ' lhs vector pushed to mc core-memory'
        write(lupri,'(a/)') ' *** end of LUCITA-MCSCF interface   ***'
      end if

  end subroutine mcscf_post_lucita_return_e2b
!*******************************************************************************

  subroutine mcscf_post_lucita_rotate_cref(c_or_cr,print_lvl)
!*******************************************************************************
!
!    purpose: post-LUCITA processing for CI task: rotate  Cvec
!              
!
!*******************************************************************************
! lucita
  use lucita_cfg
! sirius
! nothing
#include "priunit.h"
! hjaaj Apr 2016: intent(inout) does not work with paramter vdummy in parameter list
!     real(8),   intent(inout) :: c_or_cr(*)
      real(8)                  :: c_or_cr(*)
      integer,   intent(in)    :: print_lvl
!-------------------------------------------------------------------------------
      integer                  :: push_pull = 1
!-------------------------------------------------------------------------------
#ifdef VAR_MPI      
      vector_exchange_type1 = 1
#else
      vector_exchange_type1 = 3
#endif
      vector_update_mc2lu_lu2mc((push_pull-1)*vector_exchange_types+vector_exchange_type1) = .true. ! pull the outgoing lhs vector from LUCITA files to mc core-memory

!     lhs vector (the first argument '1' refers to the direction of transfer: lucita ==> mcscf)
      call vector_exchange_driver(push_pull,vector_exchange_type1,lucita_cfg_nr_roots,lucita_cfg_hcsym,               &
                                  io2io_vector_exchange_mc2lu_lu2mc,                                                  &
                                  vector_update_mc2lu_lu2mc((1-1)*vector_exchange_types+vector_exchange_type1),       &
                                  .true.,c_or_cr)

      if(print_lvl >= print_lvl_limit)then
        write(lupri,'(/a)') ' *** LUCITA-MCSCF interface reports: ***'
        write(lupri,*) ' lhs vector pushed to mc core-memory'
        write(lupri,'(a/)') ' *** end of LUCITA-MCSCF interface   ***'
      end if

  end subroutine mcscf_post_lucita_rotate_cref
!*******************************************************************************

  subroutine mcscf_post_lucita_xpdens(hc_or_cl,c_or_cr,print_lvl)
!*******************************************************************************
!
!    purpose: post-LUCITA processing for CI task: Xpdens m
!              
!
!*******************************************************************************
! lucita
  use lucita_cfg
#ifdef MOD_SRDFT
  use lucita_mcscf_srdftci_cfg
#endif
! sirius
! nothing
#include "priunit.h"
! hjaaj Apr 2016: intent(inout) does not work with paramter vdummy in parameter list
!     real(8),   intent(inout) :: hc_or_cl(*)
!     real(8),   intent(inout) :: c_or_cr(*)
      real(8)                  :: hc_or_cl(*)
      real(8)                  :: c_or_cr(*)
      integer,   intent(in)    :: print_lvl
!-------------------------------------------------------------------------------
      integer                  :: push_pull = 1
!-------------------------------------------------------------------------------
      
!     write(lupri,*) ' post xpdens: restore_cref, switch ==>', restore_cref, restore_cref_vector_switch
 
      if(srdft_ci_1pdens_cref_restore)then
        restore_cref = .false.
      end if
        if(restore_cref_vector_switch == 1)then

          vector_exchange_type1                      = 2
          vector_update_mc2lu_lu2mc((push_pull-1)*vector_exchange_types+vector_exchange_type1) = .true. ! pull the outgoing lhs vector from LUCITA files to mc core-memory

!         write(lupri,*) ' restore lhs activated',vector_exchange_type1

!         lhs vector (the first argument '1' refers to the direction of transfer: lucita ==> mcscf)
          call vector_exchange_driver(push_pull,vector_exchange_type1,lucita_cfg_nr_roots,lucita_cfg_hcsym,              &
                                      io2io_vector_exchange_mc2lu_lu2mc,                                                 &
                                      vector_update_mc2lu_lu2mc((1-1)*vector_exchange_types+vector_exchange_type1),      &
                                      .true.,hc_or_cl)

        else if(restore_cref_vector_switch == 2)then

          vector_exchange_type1                      = 2
          vector_update_mc2lu_lu2mc((push_pull-1)*vector_exchange_types+vector_exchange_type1) = .true. ! pull the outgoing lhs vector from LUCITA files to mc core-memory

!         write(lupri,*) ' restore rhs activated',vector_exchange_type1

!         rhs vector (the first argument '1' refers to the direction of transfer: lucita ==> mcscf)
          call vector_exchange_driver(push_pull,vector_exchange_type1,lucita_cfg_nr_roots,lucita_cfg_csym,                 &
                                      io2io_vector_exchange_mc2lu_lu2mc,                                                   &
                                      vector_update_mc2lu_lu2mc((1-1)*vector_exchange_types+vector_exchange_type1),        &
                                      .true.,c_or_cr)
        end if

      if(print_lvl >= print_lvl_limit)then
        write(lupri,'(/a)') ' *** LUCITA-MCSCF interface reports: ***'
        write(lupri,*) ' lhs vector pushed to mc core-memory'
        write(lupri,'(a/)') ' *** end of LUCITA-MCSCF interface   ***'
      end if

  end subroutine mcscf_post_lucita_xpdens
!*******************************************************************************

  subroutine mcscf_post_lucita_setci(print_lvl)
!*******************************************************************************
!
!    purpose: post-LUCITA processing for CI task: return CIdim 
!              
!
!*******************************************************************************
! lucita
  use lucita_cfg
#include "priunit.h"
! sirius
#include "ciinfo.h"
! lucita
#include "mxpdim.inc"
#include "cicisp.inc"
      integer,   intent(in)    :: print_lvl
!
!-------------------------------------------------------------------------------
      integer                  :: i
!-------------------------------------------------------------------------------

!     turn off common block (re-)initialization + co-worker synchronization of static lucita variables 
      lucita_cfg_initialize_cb        = .false. 
!     lucita_cfg_initialize_cb        = .true. 
!debuglucita_cfg_initialize_cb        = .true. 
      
      mxndt = 0
      do i = 1, mxpirr
        ndtasm(i) = nint(xispsm(i,1))
        ncsasm(i) = 0 ! no CSFs in LUCITA yet
        mxndt     = max(mxndt, ndtasm(i))
      end do

      if(print_lvl >= print_lvl_limit)then
        write(lupri,'(/a)') ' *** LUCITA-MCSCF interface reports: ***'
        write(lupri,*) ' maximum number of determinants in all spaces: ', MXNDT
        write(lupri,'(a/)') ' *** end of LUCITA-MCSCF interface   ***'
      end if

  end subroutine mcscf_post_lucita_setci
!*******************************************************************************

end module
