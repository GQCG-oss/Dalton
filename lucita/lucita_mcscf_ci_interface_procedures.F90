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

  implicit none

  public mcscf_lucita_interface

  private

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
      real(8),   intent(inout) :: c_or_cr(*)
      real(8),   intent(inout) :: hc_or_cl(*)
      real(8),   intent(inout) :: int1_or_rho1(*)
      real(8),   intent(inout) :: int2_or_rho2(*)
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
      real(8),   intent(inout) :: c_or_cr(*)
      real(8),   intent(inout) :: hc_or_cl(*)
      real(8),   intent(inout) :: int1_or_rho1(*)
      real(8),   intent(inout) :: int2_or_rho2(*)
      character, intent(in)    :: run_type*(*)
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

!     a. define allocations for a given CI runtype needed by slaves in parallel runs
      call setup_lucita_mcci_wrkspc_dimensions(run_type,         &
                                               print_lvl)

!     b. select pre-process task
      select case(run_type)
        case('return CIdim', 'return CIdia')
!         do nothing
        case('initial ci  ')
          call mcscf_pre_lucita_cistart(print_lvl)
        case('Xp-density m')
          call mcscf_pre_lucita_xpdens(c_or_cr,hc_or_cl,print_lvl)
        case('sigma vec   ')
          call mcscf_pre_lucita_return_e2b(c_or_cr,print_lvl)
        case('analyze Cvec')
          call mcscf_pre_lucita_bvec_analyze(c_or_cr,print_lvl)
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
      real(8),   intent(inout) :: c_or_cr(*)
      real(8),   intent(inout) :: hc_or_cl(*)
      real(8),   intent(inout) :: int1_or_rho1(*)
      real(8),   intent(inout) :: int2_or_rho2(*)
      character, intent(in)    :: run_type*(*)
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

      select case(run_type)
        case('Xp-density m', 'analyze Cvec')
!         do nothing
        case('return CIdim')
          call mcscf_post_lucita_setci(print_lvl)
        case('initial ci  ')
          call mcscf_post_lucita_cistart(c_or_cr,print_lvl)
        case('return CIdia')
          call mcscf_post_lucita_hdiag(c_or_cr,print_lvl)
        case('sigma vec   ')
          call mcscf_post_lucita_return_e2b(hc_or_cl,print_lvl)
        case default
          call quit('undefined post LUCITA processing step in MCSCF process flow')
      end select

  end subroutine mcscf_post_lucita_processing
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
      real(8),   intent(inout) :: c_or_cr(*)
      integer,   intent(in)    :: print_lvl
!-------------------------------------------------------------------------------
      integer                  :: exchange_type1
!-------------------------------------------------------------------------------
      
      exchange_type1                      = 2
      vector_update_mc2lu(exchange_type1) = .true. ! save the incoming rhs vector to LUCITA files

      if(vector_update_mc2lu(exchange_type1))then
!       rhs vector (the first argument '2' refers to the direction of transfer: mcscf ==> lucita)
        call vector_exchange_driver(2,exchange_type1,lucita_cfg_nr_roots,lucita_cfg_csym,         &
                                    io2io_vector_exchange_mc2lu_lu2mc,c_or_cr)
      end if

      if(print_lvl >= -1)then
        write(lupri,'(/a)') ' *** LUCITA-MCSCF interface reports: ***'
        write(lupri,*) ' rhs vector saved on lucita files'
        write(lupri,'(a/)') ' *** end of LUCITA-MCSCF interface   ***'
      end if

!     reset
      vector_update_mc2lu(exchange_type1) = .false.

  end subroutine mcscf_pre_lucita_bvec_analyze
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
      vector_update_lu2mc(2) = .true. ! copy the final C-vectors from LUCITA files into the MCSCF input array  

      if(print_lvl >= -1)then
        write(lupri,'(/a)') ' *** LUCITA-MCSCF interface reports: ***'
        write(lupri,*) ' inactive energy obtained from MCSCF environment ==> ', einact_mc2lu
        write(lupri,'(a/)') ' *** end of LUCITA-MCSCF interface   ***'
      end if

  end subroutine mcscf_pre_lucita_cistart
!*******************************************************************************

  subroutine mcscf_pre_lucita_return_e2b(c_or_cr,print_lvl)
!*******************************************************************************
!
!    purpose: pre-LUCITA processing for CI task: sigma vec   
!              
!
!*******************************************************************************
! lucita
  use lucita_cfg
! sirius
! nothing
#include "priunit.h"
      real(8),   intent(inout) :: c_or_cr(*)
      integer,   intent(in)    :: print_lvl
!-------------------------------------------------------------------------------
      integer                  :: exchange_type1
!-------------------------------------------------------------------------------
      
      exchange_type1                      = 2
      vector_update_mc2lu(exchange_type1) = .true. ! save the incoming rhs vector to LUCITA files

      if(vector_update_mc2lu(exchange_type1))then
!       rhs vector (the first argument '2' refers to the direction of transfer: mcscf ==> lucita)
        call vector_exchange_driver(2,exchange_type1,lucita_cfg_nr_roots,lucita_cfg_csym,         &
                                    io2io_vector_exchange_mc2lu_lu2mc,c_or_cr)
      end if

      if(print_lvl >= -1)then
        write(lupri,'(/a)') ' *** LUCITA-MCSCF interface reports: ***'
        write(lupri,*) ' rhs vector saved on lucita files'
        write(lupri,'(a/)') ' *** end of LUCITA-MCSCF interface   ***'
      end if

!     reset
      vector_update_mc2lu(exchange_type1) = .false.

  end subroutine mcscf_pre_lucita_return_e2b
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
      real(8),   intent(inout) :: c_or_cr(*)
      real(8),   intent(inout) :: hc_or_cl(*)
      integer,   intent(in)    :: print_lvl
!-------------------------------------------------------------------------------
      integer                  :: exchange_type1
      integer                  :: exchange_type2
!-------------------------------------------------------------------------------
      
      exchange_type1                      = 2
      exchange_type2                      = 3
      vector_update_mc2lu(exchange_type1) = .true. ! save the incoming rhs vector to LUCITA files
      vector_update_mc2lu(exchange_type2) = .true. ! save the incoming lhs vector to LUCITA files

      if(vector_update_mc2lu(exchange_type1))then
!       rhs vector (the first argument '2' refers to the direction of transfer: mcscf ==> lucita)
        call vector_exchange_driver(2,exchange_type1,lucita_cfg_nr_roots,lucita_cfg_csym,         &
                                    io2io_vector_exchange_mc2lu_lu2mc,c_or_cr)
      end if

      if(vector_update_mc2lu(exchange_type2) .and. lucita_cfg_transition_densm)then
!       lhs vector (if transition density matrix) - otherwise we take care of it internally in lucita
        call vector_exchange_driver(2,exchange_type2,lucita_cfg_nr_roots,lucita_cfg_hcsym,        &
                                    io2io_vector_exchange_mc2lu_lu2mc,hc_or_cl)
      end if

      if(print_lvl >= -1)then
        write(lupri,'(/a)') ' *** LUCITA-MCSCF interface reports: ***'
        write(lupri,*) ' rhs/lhs vectors saved on lucita files'
        write(lupri,'(a/)') ' *** end of LUCITA-MCSCF interface   ***'
      end if

!     reset
      vector_update_mc2lu(exchange_type1) = .false.
      vector_update_mc2lu(exchange_type2) = .false.

  end subroutine mcscf_pre_lucita_xpdens
!*******************************************************************************

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
      integer                  :: exchange_type
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
      exchange_type = 2
      if(vector_update_lu2mc(exchange_type))then
!       (the first argument '1' refers to the direction of transfer: lucita ==> mcscf)
        call vector_exchange_driver(1,exchange_type,nroot,irefsm_c,io2io_vector_exchange_mc2lu_lu2mc,c_or_cr)
      end if

      if(print_lvl >= -1)then
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

!     reset
      vector_update_lu2mc(exchange_type) = .false.

  end subroutine mcscf_post_lucita_cistart
!*******************************************************************************

  subroutine mcscf_post_lucita_return_e2b(hc_or_cl,print_lvl)
!*******************************************************************************
!
!    purpose: pre-LUCITA processing for CI task: sigma vec   
!              
!
!*******************************************************************************
! lucita
  use lucita_cfg
! sirius
! nothing
#include "priunit.h"
      real(8),   intent(inout) :: hc_or_cl(*)
      integer,   intent(in)    :: print_lvl
!-------------------------------------------------------------------------------
      integer                  :: exchange_type1
!-------------------------------------------------------------------------------
      
      exchange_type1                      = 3
      vector_update_lu2mc(exchange_type1) = .true. ! pull the outgoing lhs vector from LUCITA files to mc core-memory

      if(vector_update_lu2mc(exchange_type1))then
!       lhs vector (the first argument '1' refers to the direction of transfer: lucita ==> mcscf)
        call vector_exchange_driver(1,exchange_type1,lucita_cfg_nr_roots,lucita_cfg_hcsym,        &
                                    io2io_vector_exchange_mc2lu_lu2mc,hc_or_cl)
      end if

      if(print_lvl >= -1)then
        write(lupri,'(/a)') ' *** LUCITA-MCSCF interface reports: ***'
        write(lupri,*) ' lhs vector pushed to mc core-memory'
        write(lupri,'(a/)') ' *** end of LUCITA-MCSCF interface   ***'
      end if

!     reset
      vector_update_lu2mc(exchange_type1) = .false.

  end subroutine mcscf_post_lucita_return_e2b
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
      integer                  :: exchange_type
!-------------------------------------------------------------------------------
      
!     type 4 - pull H diag from lucita file LUCITA_HDIAG.x where x refers to the present symmetry id (a,b,c,...)
!              to mcscf core memory 
      exchange_type = 4
      vector_update_lu2mc(exchange_type) = .true.
      if(vector_update_lu2mc(exchange_type))then
!       (the first argument '1' refers to the direction of transfer: lucita ==> mcscf)
        call vector_exchange_driver(1,exchange_type,lucita_cfg_nr_roots,lucita_cfg_csym,          &
                                    io2io_vector_exchange_mc2lu_lu2mc,c_or_cr)
      end if

      if(print_lvl >= -1)then
        write(lupri,'(/a)') ' *** LUCITA-MCSCF interface reports: ***'
        write(lupri,*) ' H diag pushed into mcscf core memory '
        write(lupri,'(a/)') ' *** end of LUCITA-MCSCF interface   ***'
      end if

!     reset
      vector_update_lu2mc(exchange_type) = .false.

  end subroutine mcscf_post_lucita_hdiag
!*******************************************************************************

  subroutine mcscf_post_lucita_setci(print_lvl)
!*******************************************************************************
!
!    purpose: post-LUCITA processing for CI task: return CIdim 
!              
!
!*******************************************************************************
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
      
      mxndt = 0
      do i = 1, mxpirr
        ndtasm(i) = nint(xispsm(i,1))
        ncsasm(i) = 0 ! no CSFs in LUCITA yet
        mxndt     = max(mxndt, ndtasm(i))
      end do

      if(print_lvl >= 1)then
        write(lupri,'(/a)') ' *** LUCITA-MCSCF interface reports: ***'
        write(lupri,*) ' maximum number of determinants in all spaces: ', MXNDT
        write(lupri,'(a/)') ' *** end of LUCITA-MCSCF interface   ***'
      end if

  end subroutine mcscf_post_lucita_setci
!*******************************************************************************

end module
