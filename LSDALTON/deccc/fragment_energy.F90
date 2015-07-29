!> @file
!> Module contains operations involving fragment energy calculations
!> \author Ida-Marie Hoeyvik, Kasper Kristensen and Pablo Baudin

!> Module contains operations involving fragment energy calculations
module fragment_energy_module

  use fundamental
  use precision
  use lstiming!, only: lstimer
  use typedeftype!, only: lsitem
  use memory_handling!, only: mem_alloc,mem_dealloc,collect_thread_memory,&
  !       & mem_TurnOffThread_Memory,mem_TurnONThread_Memory,init_threadmemvar
  use dec_typedef_module
  use rpa_module

  ! DEC DEPENDENCIES (within deccc directory)
  ! ****************************************
  use dec_fragment_utils
  use array2_simple_operations
  use array4_simple_operations
  use orbital_operations
  use mp2_module !,only: max_batch_dimension,get_vovo_integrals, &
  !       & mp2_integrals_and_amplitudes
  use rimp2_module
  use dec_ls_thc_rimp2_module
  use atomic_fragment_operations!  ,only: atomic_fragment_init_basis_part, &
  !       & get_fragmentt1_AOSAOS_from_full, extract_specific_fragmentt1, &
  !       & update_full_t1_from_atomic_frag,which_pairs, &
  !       & update_full_t1_from_atomic_frag, update_full_t1_from_pair_frag
  use mp2_gradient_module ,only: single_calculate_mp2gradient_driver,&
         & pair_calculate_mp2gradient_driver
  use ccdriver, only: mp2_solver,fragment_ccsolver
#ifdef MOD_UNRELEASED
  use f12_integrals_module
  use rif12_integrals_module
  use ccsd_gradient_module
  use ccsdpt_module, only:ccsdpt_driver,ccsdpt_energy_e5_frag,&
         & ccsdpt_energy_e5_pair, ccsdpt_decnp_e5_frag
#endif

public :: optimize_atomic_fragment, pair_driver_singles, atomic_driver, &
        & pair_driver, atomic_driver_advanced, plot_pair_energies 
private

contains


  !> \brief Wrapper for atomic_driver with the following special features:
  !> 1. It is assumed that the input fragment has been initialized but that
  !> the fragment basis information (expensive box in decfrag type) has not been set.
  !> 2. The fragment basis information is calculated here and then freed again.
  !> 3. This wrapper can also attach exisiting full molecular singles amplitudes
  !>    to the fragment structure and update new improved full molecular singles amplitudes
  !>    by the calculated fragment singles amplitudes.
  !> \author Kasper Kristensen
  !> \date March 2012
  subroutine atomic_driver_advanced(nocc,nvirt,&
       & OccOrbitals,virtOrbitals,MyLsitem,MyMolecule,MyFragment,grad,&
       & t1old,t1new)

    implicit none
    !> Atomic fragment 
    type(decfrag), intent(inout) :: myfragment
    !> MP2 gradient structure (only calculated if DECinfo%first_order is turned on)
    type(mp2grad),intent(inout) :: grad
    !> Number of occupied orbitals in full molecule
    integer,intent(in) :: nOcc
    !> Number of virtupied orbitals in full molecule
    integer,intent(in) :: nvirt
    !> Occupied orbitals for full molecule
    type(decorbital), intent(in) :: OccOrbitals(nocc)
    !> virtupied orbitals for full molecule
    type(decorbital), intent(in) :: virtOrbitals(nvirt)
    !> Full molecule info (will be unchanged at output)
    type(fullmolecule), intent(inout) :: MyMolecule
    !> LSDalton info
    type(lsitem), intent(inout) :: mylsitem
    !> Existing full molecular t1 amplitudes
    type(array2),intent(in),optional :: t1old
    !> New full molecular t1 amplitudes which will be updated with the contribution from "MyFragment"
    type(array2),intent(inout),optional :: t1new

    ! Init fragment basis information
    call atomic_fragment_init_basis_part(nvirt, nocc, OccOrbitals,&
         & virtOrbitals,MyMolecule,mylsitem,MyFragment)
     
    ! Attach fragments singles amplitudes to fragment structure
    ! if long-range singles polarization effects are requested.
    if(DECinfo%SinglesPolari) then
       call get_fragmentt1_AOSAOS_from_full(MyFragment,t1old)
    end if

    ! Call main driver to get energy (and possibly density or gradient)
    call atomic_driver(MyMolecule,mylsitem,OccOrbitals,virtOrbitals,&
         & MyFragment,grad=grad)
    
    ! Update full molecular singles amplitudes with (virt EOS,occ EOS) fragment contributions
    if(DECinfo%SinglesPolari) then
       ! Extract (virt EOS,occ EOS) indices from fragment
       call extract_specific_fragmentt1(MyFragment,.true.,.true.)
       ! Update full molecular singles amplitudes
       call update_full_t1_from_atomic_frag(MyFragment,t1new)
    end if

    ! Free basis info and t1 info again
    call atomic_fragment_free_basis_info(MyFragment)
    if(DECinfo%SinglesPolari) call free_fragment_t1(MyFragment)

  end subroutine atomic_driver_advanced


  !> \brief Driver for calculating atomic fragment energy for a given fragment using
  !> the Lagrangian approach.
  !> If requested, first order properties (MP2 density or gradient) are also calculated.
  !> \author Kasper Kristensen
  subroutine atomic_driver(MyMolecule,mylsitem,OccOrbitals,virtOrbitals,&
       & MyFragment,grad)

    implicit none
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> LS item info                                                                                    
    type(lsitem), intent(inout) :: mylsitem
    !> Information about DEC occupied orbitals                                                         
    type(decorbital), dimension(MyMolecule%nocc), intent(in) :: OccOrbitals
    !> Information about DEC virtupied orbitals
    type(decorbital), dimension(MyMolecule%nvirt), intent(in) :: virtOrbitals
    !> Atomic fragment 
    type(decfrag), intent(inout) :: myfragment
    !> MP2 gradient structure (only calculated if DECinfo%first_order is turned on)
    type(mp2grad),intent(inout),optional :: grad

    ! Sanity check
    if(DECinfo%first_order) then
       if(.not. present(grad)) then
          call lsquit('atomic_driver: Gradient argument not present!',-1)
       end if
       call fragment_energy_and_prop(MyFragment,grad=grad)
    else
       call fragment_energy_and_prop(MyFragment)
    end if
 
  end subroutine atomic_driver


  !> \brief Driver for calculating atomic fragment energy for a given fragment using the Lagrangian approach.
  !> If requested, first order properties (MP2 density or gradient) are also calculated and saved.
  !> \author Kasper Kristensen
  subroutine fragment_energy_and_prop(MyFragment,Fragment1,Fragment2,grad,add_MP2_opt)

    implicit none
    !> Atomic or pair fragment
    type(decfrag), intent(inout) :: myfragment
    !> Fragments 1 and 2 used to form pair fragment - only used if myfragment is pair fragment
    type(decfrag), intent(in),optional :: Fragment1, Fragment2
    !> MP2 gradient structure (only calculated if DECinfo%first_order is turned on)
    type(mp2grad),intent(inout),optional :: grad
    !> Get MP2 info on top of actual CC model:
    logical, intent(in), optional :: add_mp2_opt
    type(tensor) :: t1, ccsdpt_t1, m1, eta
    !> integrals
    type(tensor) :: VOVO,VOVOocc,VOVOvirt,VOOO,VOVV,VOVOvirtTMP
    !> doubles amplitudes
    type(tensor) :: t2,u,t2occ,t2virt,t2occp,t2virtp,t2MP2,t2MP2o,t2MP2v,t2_occEOS
    !> ccsdpt arrays
    type(tensor) :: ccsdpt_t2, ccsdpt_t2occ, ccsdpt_t2virt
    type(tensor) :: m2,m2occ,m2virt

    real(realk) :: tcpu, twall,debugenergy
    ! timings are allocated and deallocated behind the curtains
    real(realk),pointer :: times_ccsd(:), times_pt(:)
    logical :: print_frags,abc,pair,get_mp2, use_bg
    integer :: a,b,i,j,k,l, ccmodel, nvA,noA,nvE,noE

    get_mp2 = .false.
    if (present(add_mp2_opt)) get_mp2 = add_mp2_opt
    use_bg  = mem_is_background_buf_init()

    ! Pairfragment?
    pair = MyFragment%pairfrag
    if(pair) then
       if( (.not. present(Fragment1)) .or. (.not. present(Fragment2)) ) then
          call lsquit('fragment_energy_and_prop: Missing arguments for pair fragment!',-1)
       end if

       if (DECinfo%DECNP) then
          call lsquit("ERROR(fragment_energy_and_prop): DECNP inconsistent input",DECinfo%output)
       end if
    end if


    times_ccsd => null()
    times_pt   => null()
    call LSTIMER('START',tcpu,twall,DECinfo%output)

    nvA = MyFragment%nvirtAOS
    noA = MyFragment%noccAOS

    nvE = MyFragment%nvirtEOS
    noE = MyFragment%noccEOS

    ! Which model? MP2,CC2, CCSD etc.
    ! *******************************
    !MODIFY FOR NEW MODEL
    WhichCCmodel: select case(MyFragment%ccmodel)

    case(MODEL_NONE) ! SKip calculation

       return

    case(MODEL_MP2) ! MP2 calculation

          if(DECinfo%first_order .and. (.not. DECinfo%unrelaxed) ) then  
             ! calculate also MP2 density integrals
             ! ************************************
             call MP2_integrals_and_amplitudes(MyFragment,VOVOocc,t2occ,VOVOvirt,t2virt,VOOO,VOVV)
             print *,"1. Myfragment%energies(FRAGMODEL_OCCMP2) (fragment energy)",  Myfragment%energies(FRAGMODEL_OCCMP2)

          else if (DECinfo%DECNP) then
             ! Get MP2 amplitudes and integrals using ccsolver
             ! ***********************************************
             call mp2_solver(MyFragment,VOVO,t2,.false.)
             print *,"2. Myfragment%energies(FRAGMODEL_OCCMP2) (fragment energy)",  Myfragment%energies(FRAGMODEL_OCCMP2)

             ! Extract EOS indices for amplitudes and integrals
             ! ************************************************
             ! Note: EOS space is different for DECNP !!
             call tensor_extract_decnp_indices(t2,MyFragment,t2occ,t2virt)
             call tensor_extract_decnp_indices(VOVO,MyFragment,VOVOocc,VOVOvirt)
           
             ! free stuff
             ! **********
             call tensor_free(VOVO)
             call tensor_free(t2)
          else 
             ! calculate only MP2 energy integrals and MP2 amplitudes
             ! ******************************************************
             call MP2_integrals_and_amplitudes(MyFragment,VOVOocc,t2occ,VOVOvirt,t2virt)
             print *,"3. Myfragment%energies(FRAGMODEL_OCCMP2) (fragment energy)",  Myfragment%energies(FRAGMODEL_OCCMP2)
          end if
#ifdef MOD_UNRELEASED
          ! MP2-F12 Code
          MP2F12: if(DECinfo%F12) then    
             if(pair) then
                if(MyFragment%isopt) then
                   call get_f12_fragment_energy(MyFragment, t2occ%elm4, t1%elm2, MyFragment%ccmodel,&
                      &  Fragment1, Fragment2)
                end if
             else
                ! Calculate F12 only for optimized fragment or if it has been requested by input
                if(MyFragment%isopt .or. DECinfo%F12fragopt) then
                   call get_f12_fragment_energy(MyFragment, t2occ%elm4, t1%elm2, MyFragment%ccmodel)
                end if
             end if
             !> Free cabs after each calculation
             call free_cabs()
          end if MP2F12
#endif

    case(MODEL_RIMP2) ! RIMP2 calculation

       if(DECinfo%first_order .and. (.not. DECinfo%unrelaxed) ) then  
          ! calculate also RIMP2 density integrals
          call RIMP2_integrals_and_amplitudes(MyFragment,VOVOocc,t2occ,VOVOvirt,t2virt,VOOO,VOVV)
       else if (DECinfo%DECNP) then
          ! calculate also RIMP2 integrals (full AOS)
          call decnp_RIMP2_integrals_and_amplitudes(MyFragment,VOVO,t2)
           
          ! Extract EOS indices for amplitudes and integrals
          ! ************************************************
          ! Note: EOS space is different for DECNP !!
          call tensor_extract_decnp_indices(t2,MyFragment,t2occ,t2virt)
          call tensor_extract_decnp_indices(VOVO,MyFragment,VOVOocc,VOVOvirt)
          
          ! free stuff
          ! **********
          call tensor_free(VOVO)
          call tensor_free(t2)
       else
          ! calculate also RIMP2 integrals
          call RIMP2_integrals_and_amplitudes(MyFragment,VOVOocc,t2occ,VOVOvirt,t2virt)
       endif

#ifdef MOD_UNRELEASED
          ! MP2-F12 Code
          RIMP2F12: if(DECinfo%F12) then    
             if(pair) then
                if(MyFragment%isopt) then
                   call get_rif12_fragment_energy(MyFragment, t2occ%elm4, t1%elm2, MyFragment%ccmodel,&
                      &  Fragment1, Fragment2)
                end if
             else
                ! Calculate F12 only for optimized fragment or if it has been requested by input
                if(MyFragment%isopt .or. DECinfo%F12fragopt) then
                   call get_rif12_fragment_energy(MyFragment, t2occ%elm4, t1%elm2, MyFragment%ccmodel)
                end if
             end if
             !> Free cabs after each calculation
             call free_cabs()
          end if RIMP2F12
#endif




    case(MODEL_LSTHCRIMP2) ! LSTHCRIMP2 calculation

       if(DECinfo%first_order)call lsquit('no first order LSTHCRIMP2',-1)       
       call LSTHCRIMP2_integrals_and_amplitudes(MyFragment,VOVOocc,t2occ,VOVOvirt,t2virt)

    case(MODEL_CC2,MODEL_CCSD,MODEL_CCSDpT,MODEL_RPA,MODEL_SOSEX) ! higher order CC (-like)

       call dec_fragment_time_init(times_ccsd)

       ! Solve CC equation to calculate amplitudes and integrals 
       ! *******************************************************
       ! Here all output indices in t1,t2, and VOVO are AOS indices.
       ! calculate also MP2 density integrals
#ifdef MOD_UNRELEASED
       if(DECinfo%first_order.or.DECinfo%CCSDmultipliers) then  
          call fragment_ccsolver(MyFragment,t1,t2,VOVO,m1=m1,m2=m2)
          !extract the indices here for later use in CCSD first order properties
          !note: the t2occp and t2virtp are the t2 without further addmixture of t1!
          call tensor_extract_eos_indices(m2,MyFragment,tensor_occEOS=m2occ,tensor_virtEOS=m2virt)
          call tensor_extract_eos_indices(t2,MyFragment,tensor_occEOS=t2occp,tensor_virtEOS=t2virtp)
       else
#endif
          call fragment_ccsolver(MyFragment,t1,t2,VOVO)
          if (get_mp2) then
             call mp2_solver(MyFragment,VOVO,t2MP2,.true.)
          end if
#ifdef MOD_UNRELEASED
       endif

       if(DECinfo%first_order) then
          !call z_vec_rhs_ccsd(MyFragment,eta,m1_arr=m1,m2_arr=m2,t1_arr=t1,t2_arr=t2)
          !call tensor_free(m2)
          !call tensor_free(m1)
       endif
#endif

       call dec_fragment_time_get(times_ccsd)

#ifdef MOD_UNRELEASED
       ! calculate ccsd(t) fragment energies
       ! ***********************************
       if(MyFragment%ccmodel==MODEL_CCSDpT) then
          call dec_fragment_time_init(times_pt)

          if (pair) then
             call get_ccsdpt_fragment_energies(MyFragment,VOVO,t2,t1,pair,fragment1,fragment2)
          else
             call get_ccsdpt_fragment_energies(MyFragment,VOVO,t2,t1,pair)
          end if

          call dec_fragment_time_get(times_pt)
       end if

       ! CCSD-F12 Code
       CCSDF12: if(DECinfo%F12) then

          call tensor_extract_eos_indices(t2,MyFragment,tensor_occEOS=t2_occEOS)
          if(pair) then
             if(MyFragment%isopt) then
                call get_f12_fragment_energy(MyFragment, t2_occEOS%elm4, t1%elm2, MyFragment%ccmodel)  
             end if
          else
             ! Calculate F12 only for optimized fragment or if it has been requested by input
             if(MyFragment%isopt .or. DECinfo%F12fragopt) then
                call get_f12_fragment_energy(MyFragment, t2_occEOS%elm4, t1%elm2, MyFragment%ccmodel)  
             end if
          end if

          !> Free cabs after each calculation
          call tensor_free(t2_occEOS)
          call free_cabs()

       endif CCSDF12
#endif

       ! Calculate combined single+doubles amplitudes
       ! ********************************************
       ! u(a,i,b,j) = t2(a,i,b,j) + t1(a,i)*t1(b,j)          
       call get_combined_SingleDouble_amplitudes(t1,t2,u)


       if (DECinfo%DECNP) then
          ! Extract EOS indices for amplitudes and integrals
          ! ************************************************
          ! Note: EOS space is different for DECNP !!
          call tensor_extract_decnp_indices(u,MyFragment,t2occ,t2virt)
          call tensor_extract_decnp_indices(VOVO,MyFragment,VOVOocc,VOVOvirt)

       else
          ! Extract EOS indices for amplitudes and integrals
          ! ************************************************
          ! Note, t2occ and t2virt also contain singles contributions
          call tensor_extract_eos_indices(u,MyFragment,tensor_occEOS=t2occ,tensor_virtEOS=t2virt)
          call tensor_extract_eos_indices(VOVO,MyFragment,tensor_occEOS=VOVOocc, &
             & tensor_virtEOS=VOVOvirt)

          ! Extract also EOS indices from MP2 amplitudes:
          if (get_mp2) then
             call tensor_extract_eos_indices(t2MP2,MyFragment,tensor_occEOS=t2MP2o,tensor_virtEOS=t2MP2v)
             call tensor_free(t2MP2)
          end if
           
       end if

       ! free stuff
       ! **********
       call tensor_free(u)
       if(DECinfo%use_singles)then
          call tensor_free(t1)
       endif
       ! the order of deallocation is important!
       call tensor_free(t2)
       call tensor_free(VOVO)

    case default

       call lsquit("ERROR(fragment_energy_and_prop):MODEL not implemented",-1)

    end select WhichCCmodel


    ! Calcuate atomic or pair fragment energy
    ! ***************************************

    ! MODIFY FOR NEW MODEL!
    ! Two possible situations:
    ! (1) Your new model fits into the standard CC energy expression (this includes MP2). 
    !     Things should work out of the box by calling get_atomic_fragment_energy or
    !     get_pair_fragment_energy below.
    ! (2) Your model does NOT fit into the standard CC energy expression.
    !     In this case you need to make new atomic/pair fragment energy subroutines and call it from
    !     here instead of calling get_atomic_fragment_energy and get_pair_fragment_energy.
    !
    ! For frozen core and first order properties we need to remove core indices from VOVOvirt, 
    ! since they, "rather incoveniently", are required for the gradient but not for the energy
    if(DECinfo%frozencore .and. DECinfo%first_order .and. (.not. DECinfo%unrelaxed)) then
     
       call remove_core_orbitals_from_last_index(MyFragment,VOVOvirt,VOVOvirtTMP)
       if(pair) then
          ! Pair fragment
          call get_pair_fragment_energy(VOVOocc,VOVOvirtTMP,t2occ,t2virt,&
               & Fragment1, Fragment2, MyFragment)
             
          !For sosex contribution, as the residual
          !is the same for sosex and drpa
          !the energies can be calculated in one step
          if(MyFragment%ccmodel == MODEL_RPA) then
            call get_pair_fragment_energy(VOVOocc,VOVOvirtTMP,t2occ,t2virt,&
              & Fragment1, Fragment2, MyFragment,.true.)
          endif
       else
          ! Atomic fragment
          if (get_mp2) then
             ccmodel = MyFragment%ccmodel
             MyFragment%ccmodel = MODEL_MP2
             call get_atomic_fragment_energy(VOVOocc,VOVOvirt,t2MP2o,t2MP2v,MyFragment)
             call tensor_free(t2MP2o)
             call tensor_free(t2MP2v)
             MyFragment%ccmodel = ccmodel
          end if
          call get_atomic_fragment_energy(VOVOocc,VOVOvirtTMP,t2occ,t2virt,MyFragment)
          !For sosex contribution, as the residual
          !is the same for sosex and drpa
          !the energies can be calculated in one step
          if(MyFragment%ccmodel == MODEL_RPA) then
            call get_atomic_fragment_energy(VOVOocc,VOVOvirtTMP,t2occ,&
               & t2virt,MyFragment,doSOS=.true.)
          endif
       end if
       call tensor_free(VOVOvirtTMP)
     
    else if (DECinfo%DECNP) then

       ! Pair and fragment energy contributions are calculated together
       ! **************************************************************
       call  get_decnp_fragment_energy(MyFragment,VOVOocc,VOVOvirt,t2occ,t2virt)

    else ! use VOVOvirt as it is -- MOST COMMON CASE
     
       if(pair) then
          ! Pair fragment
          call get_pair_fragment_energy(VOVOocc,VOVOvirt,t2occ,t2virt,&
               & Fragment1, Fragment2, MyFragment)
          !For sosex contribution, as the residual
          !is the same for sosex and drpa
          !the energies can be calculated in one step
          if(MyFragment%ccmodel == MODEL_RPA) then
             call get_pair_fragment_energy(VOVOocc,VOVOvirt,t2occ,t2virt,&
               & Fragment1, Fragment2, MyFragment,.true.)
          endif
       else
          ! Atomic fragment
          if (get_mp2) then
             ccmodel = MyFragment%ccmodel
             MyFragment%ccmodel = MODEL_MP2
             call get_atomic_fragment_energy(VOVOocc,VOVOvirt,t2MP2o,t2MP2v,MyFragment)
             call tensor_free(t2MP2o)
             call tensor_free(t2MP2v)
             MyFragment%ccmodel = ccmodel
          end if
          call get_atomic_fragment_energy(VOVOocc,VOVOvirt,t2occ,t2virt,MyFragment)
          !For sosex contribution, as the residual
          !is the same for sosex and drpa
          !the energies can be calculated in one step
          if(MyFragment%ccmodel == MODEL_RPA) then
             call get_atomic_fragment_energy(VOVOocc,VOVOvirt,&
                & t2occ,t2virt,MyFragment,doSOS=.true.)
          endif
       end if
     
    end if

    !> Memory Stats
    if(DECinfo%F12debug) then
       WRITE(DECinfo%output,*) "Memstats before F12 fragment_energy calculation"  
       call stats_globalmem(DECinfo%output)
    end if


    ! MODIFY FOR NEW CORRECTION!
    ! If you implement a new correction (e.g. F12) please insert call to subroutine either here
    ! or above below the respective models.
    ! The energy correction should be stored in myfragment%energies(?),
    ! see dec_readme file and FRAGMODEL_* definitions in dec_typedef.F90.


    call LSTIMER('SINGLE L.ENERGY',tcpu,twall,DECinfo%output)

    ! First order properties
    ! **********************
    if(DECinfo%first_order) then
#ifdef MOD_UNRELEASED
       if(DECinfo%ccmodel == MODEL_CCSD)then
         if(pair) then
           ! Pair fragment
           call pair_calculate_CCSDgradient_driver&
                             &(Fragment1,Fragment2,MyFragment,t2occp,t2virtp,m2occ,m2virt,m1,&
                                       &VOOO,VOVV,VOVOocc,VOVOvirt,grad)
           call tensor_free(m2occ)
           call tensor_free(m2virt) 
           
         else
           ! Atomic fragment
           call single_calculate_CCSDgradient_driver&
                             &(MyFragment,t2occp,t2virtp,m2occ,m2virt,m1,&
                                       &VOOO,VOVV,VOVOocc,VOVOvirt,grad)
           call tensor_free(m2occ)
           call tensor_free(m2virt)
         endif
       endif
#endif
       if(DECinfo%ccmodel == MODEL_MP2.OR.DECinfo%ccmodel == MODEL_RIMP2)then
          if(pair) then
             ! Pair fragment
             call pair_calculate_mp2gradient_driver(Fragment1,Fragment2,MyFragment,&
                  & t2occ,t2virt,VOOO,VOVV,VOVOocc,VOVOvirt,grad)
          else
             ! Atomic fragment
             call single_calculate_mp2gradient_driver(MyFragment,t2occ,t2virt,VOOO,&
                  & VOVV,VOVOocc,VOVOvirt,grad)
          end if
!          if(.not. DECinfo%unrelaxed) then
!             call tensor_free(VOOO)
!             call tensor_free(VOVV)
!          end if
       end if
    end if
    !
    if(MyFragment%ccmodel /= MODEL_MP2.AND.MyFragment%ccmodel /= MODEL_RIMP2)then
       call dec_time_evaluate_efficiency_frag(MyFragment,times_ccsd,MODEL_CCSD,'CCSD part')
    endif
    if(MyFragment%ccmodel == MODEL_CCSDpT)then
       call dec_time_evaluate_efficiency_frag(MyFragment,times_pt,MODEL_CCSDpT,'(T)  part')
    endif

    ! Free remaining arrays
    call tensor_free(VOVOocc)
    call tensor_free(t2occ)
    call tensor_free(VOVOvirt)
    call tensor_free(t2virt)

    if(DECinfo%first_order) then
       if(DECinfo%ccmodel == MODEL_MP2.OR.DECinfo%ccmodel == MODEL_RIMP2)then
          if(.not. DECinfo%unrelaxed) then
             call tensor_free(VOOO)
             call tensor_free(VOVV)
          end if
       end if
    end if

    !print *,"s1",VOVOocc%initialized,associated(VOVOocc%elm1)
    !print *,"s2",VOVOvirt%initialized,associated(VOVOvirt%elm1)
    !print *,"s3",t2virt%initialized,associated(t2virt%elm1)
    !print *,"s4",t2occ%initialized,associated(t2occ%elm1)
    !print *,"s5",VOVO%initialized,associated(VOVO%elm1)
    !print *,"s6",t2%initialized,associated(t2%elm1)
    !print *,"s7",t1%initialized,associated(t1%elm1)
    !print *,"s8",ccsdpt_t1%initialized,associated(ccsdpt_t1%elm1)
    !print *,"s9",ccsdpt_t2%initialized,associated(ccsdpt_t2%elm1)
    !print *,"s10",u%initialized,associated(u%elm1)

  end subroutine fragment_energy_and_prop


  !> Purpose: Wrapper to ccsdpt driver and fragment energy routines
  !> Author:  Pablo Baudin
  !> Date:    May 2015
  subroutine get_ccsdpt_fragment_energies(frag,VOVO,t2,t1,pair,frag1,frag2)

     implicit none

     !> DEC fragment 
     type(decfrag), intent(inout) :: frag
     !> Two-electron integrals (a i | b j)
     type(tensor), intent(inout) :: VOVO
     !> CCSD amplitudes:
     type(tensor), intent(inout) :: t2, t1
     !> Is it a pair fragment?
     logical, intent(in) :: pair
     type(decfrag), intent(in), optional :: frag1, frag2

     type(tensor) :: ccsdpt_t1, ccsdpt_t2
     type(tensor) :: ccsdpt_t2occ, ccsdpt_t2virt, t2occ, t2virt
     logical :: abc,bg
#ifdef MOD_UNRELEASED
     bg  = mem_is_background_buf_init()

     ! Get (T) singles and doubles intermediates:
     ! ******************************************

     !FIXME: this should be decided based on the amount of memory available !!!
     abc = DECinfo%abc

     ! init ccsd(t) singles and ccsd(t) doubles (*T1 and *T2)
     if (abc) then
        call tensor_init(ccsdpt_t1,[frag%noccAOS,frag%nvirtAOS],2,bg=bg)
        call tensor_init(ccsdpt_t2,[frag%noccAOS,frag%noccAOS,&
             &frag%nvirtAOS,frag%nvirtAOS],4,bg=bg)
     else
        call tensor_init(ccsdpt_t1, [frag%nvirtAOS,frag%noccAOS],2,bg=bg)
        call tensor_init(ccsdpt_t2, [frag%nvirtAOS,frag%nvirtAOS,&
             &frag%noccAOS,frag%noccAOS],4,bg=bg)
     endif

     ! call ccsd(t) driver and single fragment evaluation
     call ccsdpt_driver(frag%noccAOS,frag%nvirtAOS,frag%nbasis,frag%ppfock,&
          & frag%qqfock,frag%Co,frag%Cv,frag%mylsitem,VOVO,t2,&
          & DECinfo%print_frags,abc,ccsdpt_singles=ccsdpt_t1,ccsdpt_doubles=ccsdpt_t2)

     if (abc) then
        call tensor_reorder(ccsdpt_t2,[3,1,4,2]) ! ccsdpt_doubles in the order (a,i,b,j)
        call tensor_reorder(ccsdpt_t1,[2,1])     ! ccsdpt_singles in the order (a,i)
     else
        call tensor_reorder(ccsdpt_t2,[1,3,2,4]) ! ccsdpt_doubles in the order (a,i,b,j)
     endif


     ! Calculate (T) fragment energy contributions:
     ! ********************************************
     if (DECinfo%DECNP) then
        ! Extract EOS indices, Note: EOS space is different for DECNP !!
        call tensor_extract_decnp_indices(t2,frag,t2occ,t2virt)
        call tensor_extract_decnp_indices(ccsdpt_t2,frag, &
           & ccsdpt_t2occ,ccsdpt_t2virt)
        ! release ccsd(t) doubles amplitudes
        call tensor_free(ccsdpt_t2)

        call get_decnp_fragment_energy(frag,ccsdpt_t2occ,ccsdpt_t2virt,t2occ,t2virt,.true.)
        call ccsdpt_decnp_e5_frag(frag,t1,ccsdpt_t1)
     else if(pair) then
        ! extract EOS indices from doubles (CCSD and (T)):
        call tensor_extract_eos_indices(ccsdpt_t2,frag,tensor_occEOS=ccsdpt_t2occ, &
           & tensor_virtEOS=ccsdpt_t2virt)
        call tensor_extract_eos_indices(t2,frag,tensor_occEOS=t2occ, &
           & tensor_virtEOS=t2virt)
        ! release ccsd(t) doubles amplitudes
        call tensor_free(ccsdpt_t2)

        call get_pair_fragment_pT4_energy(ccsdpt_t2occ,ccsdpt_t2virt,t2occ,t2virt,&
           & frag1,frag2,frag)
        call ccsdpt_energy_e5_pair(frag,t1,ccsdpt_t1)
     else
        ! extract EOS indices from doubles (CCSD and (T)):
        call tensor_extract_eos_indices(ccsdpt_t2,frag,tensor_occEOS=ccsdpt_t2occ, &
           & tensor_virtEOS=ccsdpt_t2virt)
        call tensor_extract_eos_indices(t2,frag,tensor_occEOS=t2occ, &
           & tensor_virtEOS=t2virt)
        ! release ccsd(t) doubles amplitudes
        call tensor_free(ccsdpt_t2)

        call get_single_fragment_pT4_energy(ccsdpt_t2occ,ccsdpt_t2virt,t2occ,t2virt,frag)
        call ccsdpt_energy_e5_frag(frag,t1,ccsdpt_t1)
     end if

     ! free stuff
     call tensor_free(ccsdpt_t1)
     call tensor_free(ccsdpt_t2occ)
     call tensor_free(ccsdpt_t2virt)
     call tensor_free(t2occ)
     call tensor_free(t2virt)
#else
     call lsquit("ERROR(get_ccsdpt_fragment_energies): not implemented in this version",-1)
#endif

  end subroutine get_ccsdpt_fragment_energies


  !> \brief Contract amplitudes, multipliers, and integrals to calculate atomic fragment Lagrangian energy.
  !> \author Kasper Kristensen
  !> \date August 2011
  subroutine get_atomic_fragment_energy(gocc,gvirt,t2occ,t2virt,MyFragment,doSOS)

    implicit none
    !> Two-electron integrals (a i | b j), only occ orbitals on central atom, virt AOS orbitals
    type(tensor), intent(in) :: gocc
    !> Two-electron integrals (a i | b j), only virt orbitals on central atom, occ AOS orbitals
    type(tensor), intent(in) :: gvirt
    !> amplitudes, only occ orbitals on central atom, virt AOS orbitals
    !  (combine singles and doubles if singles are present)
    type(tensor), intent(in) :: t2occ
    !> amplitudes, only virt orbitals on central atom, occ AOS orbitals
    !  (combine singles and doubles if singles are present)
    type(tensor), intent(in) :: t2virt
    !> Atomic fragment 
    type(decfrag), intent(inout) :: myfragment
    !>SOSEX cont
    logical,intent(in),optional :: doSOS

    real(realk) :: tcpu1, twall1, tcpu2,twall2, tcpu,twall
    real(realk) :: Eocc,lag_occ,Evirt,lag_virt,tmp,multaibj
    real(realk) :: prefac_coul,prefac_k
    real(realk),pointer :: occ_tmp(:),virt_tmp(:)
    real(realk),pointer :: t(:,:,:,:), g(:,:,:,:), f(:,:)
    logical ::  something_wrong,SOS,PerformLag! ,doOccPart, doVirtPart
    integer :: noccEOS,nvirtEOS,noccAOS,nvirtAOS
    integer :: i,j,k,a,b,c

    PerformLag = MyFragment%ccmodel==MODEL_MP2 .OR. MyFragment%ccmodel==MODEL_RIMP2
    ! Lagrangian energy can be split into four contributions:
    ! The first two (Eocc and lag_occ) use occupied EOS orbitals and virtual AOS orbitals.
    ! The last two (Evirt and lag_virt) use virtual EOS orbitals and occupied AOS orbitals.
    !
    ! With the restrictions above the Lagrangian energy is given by:
    ! energy  = Eocc + lag_occ + Evirt + lag_virt
    !
    ! Eocc     = sum_{ijab} t_{ij}^{ab} [ 2g_{aibj} - g_{biaj} ]
    ! lag_occ  = 1/2 sum_{ijabc} mult_{ij}^{ab} [ t_{ij}^{cb} F_{ac} + t_{ij}^{ac} F_{bc} ]
    ! Evirt    = 1/2 sum_{ijab} mult_{ij}^{ab} g_{aibj}
    ! lag_virt = - 1/2 sum_{ijkab} mult_{ij}^{ab} [ t_{kj}^{ab} F_{ki} + t_{ik}^{ab} F_{kj} ]
    !
    ! IF -and only if- the fragment is the full molecule, the following simple relations hold:
    ! energy = Eocc = Evirt = -lag_occ - lag_virt
    !
    ! Important: For MP2 the multipliers can be determined directly from the amplitudes
    !            using the simple relation:
    !
    ! mult_{ij}^{ab} = 4*t_{ij}^{ab} - 2*t_{ij}^{ba}
    !
    ! In the calculations below the factor 1/2 in the equations
    ! above is absorbed in the multipliers.


    ! MyFragment%OccContribs and MyFragment%VirtContribs contains the contributions from
    ! each individual occupied and virtual orbital -- e.g. MyFragment%OccContribs(i)
    ! is the estimated change in the energy if occupied orbital with index
    ! MyFragment%occAOSidx(i) is removed from the fragment.

    call LSTIMER('START',tcpu,twall,DECinfo%output)
    call LSTIMER('START',tcpu1,twall1,DECinfo%output)

    ! Init stuff
    noccEOS  = MyFragment%noccEOS
    nvirtEOS = MyFragment%nvirtEOS
    noccAOS  = MyFragment%noccAOS
    nvirtAOS = MyFragment%nvirtAOS
    Eocc     = 0E0_realk
    lag_occ  = 0E0_realk
    Evirt    = 0E0_realk
    lag_virt = 0E0_realk
    SOS = .false.
    if(present(doSOS)) SOS = doSOS

    if(MyFragment%ccmodel==MODEL_RPA) then
      prefac_coul=1._realk
      prefac_k=0.0_realk
      if(SOS) prefac_k=0.5_realk
    elseif(MyFragment%ccmodel==MODEL_SOSEX) then
      prefac_coul=1._realk
      prefac_k=0.5_realk
    else
       prefac_coul=2._realk
       prefac_k=1._realk
    endif


    ! Sanity checks
    ! *************
    something_wrong=.false.
    if(t2occ%dims(1) /= nvirtAOS) something_wrong=.true.
    if(t2occ%dims(2) /= noccEOS) something_wrong=.true.
    if(t2occ%dims(3) /= nvirtAOS) something_wrong=.true.
    if(t2occ%dims(4) /= noccEOS) something_wrong=.true.
    
    if(gocc%dims(1) /= nvirtAOS) something_wrong=.true.
    if(gocc%dims(2) /= noccEOS) something_wrong=.true.
    if(gocc%dims(3) /= nvirtAOS) something_wrong=.true.
    if(gocc%dims(4) /= noccEOS) something_wrong=.true.

    if(t2virt%dims(1) /= nvirtEOS) something_wrong=.true.
    if(t2virt%dims(2) /= noccAOS) something_wrong=.true.
    if(t2virt%dims(3) /= nvirtEOS) something_wrong=.true.
    if(t2virt%dims(4) /= noccAOS) something_wrong=.true.
    
    if(gvirt%dims(1) /= nvirtEOS) something_wrong=.true.
    if(gvirt%dims(2) /= noccAOS) something_wrong=.true.
    if(gvirt%dims(3) /= nvirtEOS) something_wrong=.true.
    if(gvirt%dims(4) /= noccAOS) something_wrong=.true.

    if(something_wrong) then
       print *, 't2occ%dims =', t2occ%dims
       print *, 'gocc%dims  =', gocc%dims
       print *, 't2virt%dims=', t2virt%dims
       print *, 'gvirt%dims =', gvirt%dims
       print *, 'noccEOS  = ', noccEOS
       print *, 'noccAOS  = ', noccAOS
       print *, 'nvirtEOS = ', nvirtEOS
       print *, 'nvirtAOS = ', nvirtAOS
       call lsquit('get_atomic_fragment_energy: &
            & Input dimensions do not match!',-1)
    end if


    ! Calculate Eocc and lag_occ
    ! **************************

    ! Point to tensor structure to avoid OMP problems:
    t => t2occ%elm4
    g => gocc%elm4
    f => MyFragment%qqfock

    call mem_TurnONThread_Memory()
    !$OMP PARALLEL DEFAULT(none) PRIVATE(multaibj,tmp,j,b,i,a,c,virt_tmp) &
    !$OMP SHARED(t,g,f,nvirtAOS,noccEOS,prefac_coul,prefac_k,myfragment, &
    !$OMP PerformLag) REDUCTION(+:Eocc,lag_occ)
    call init_threadmemvar()

    ! Contributions from each individual virtual orbital
    call mem_alloc(virt_tmp,nvirtAOS)
    virt_tmp = 0.0E0_realk       

    !$OMP DO SCHEDULE(dynamic,1) COLLAPSE(3)
    do j=1,noccEOS
       do b=1,nvirtAOS
          do i=1,noccEOS
             do a=1,nvirtAOS


                ! Contribution 1
                ! --------------

                ! Energy contribution for orbitals (j,b,i,a)
                tmp = t(a,i,b,j)*(prefac_coul*g(a,i,b,j) -prefac_k*g(b,i,a,j))

                ! Update total atomic fragment energy contribution 1
                Eocc = Eocc + tmp

                ! Update contribution from orbital a
                virt_tmp(a) = virt_tmp(a) + abs(tmp)
                ! Update contribution from orbital b (only if different from a to avoid double counting)
                if(a/=b) virt_tmp(b) = virt_tmp(b) + abs(tmp)


                ! Contribution 2
                ! --------------

                ! Skip contribution 2 for anything but MP2 and RIMP2
                if(PerformLag) then
                   ! Multiplier (multiplied by one half)
                   multaibj = prefac_coul*t(a,i,b,j) - prefac_k*t(b,i,a,j)

                   do c=1,nvirtAOS

                      ! Energy contribution for orbitals (j,b,i,a,c)
                      tmp = t(c,i,b,j)*f(c,a) + t(a,i,c,j)*f(c,b)
                      tmp = multaibj*tmp

                      ! Update total atomic fragment energy contribution 2
                      lag_occ = lag_occ + tmp

                      ! Update contribution from orbital a
                      !virt_tmp(a) = virt_tmp(a) + tmp
                      ! Update contribution from orbital b (only if different from a to avoid double counting)
                      !if(a/=b) virt_tmp(b) = virt_tmp(b) + tmp
                      ! Update contribution from orbital c (only if different from a and b)
                      !if( (a/=c) .and. (b/=c) ) virt_tmp(c) = virt_tmp(c) + tmp

                   end do

                end if


             end do
          end do
       end do
    end do
    !$OMP END DO NOWAIT

    ! Update total virtual contributions to fragment energy
    !$OMP CRITICAL
    do a=1,nvirtAOS
       MyFragment%VirtContribs(a) = MyFragment%VirtContribs(a) + virt_tmp(a)
    end do
    !$OMP END CRITICAL

    call mem_dealloc(virt_tmp)
    call collect_thread_memory()
    !$OMP END PARALLEL
    call mem_TurnOffThread_Memory()


    ! Calculate Evirt and lag_virt
    ! ****************************

    ! Point to tensor structure to avoid OMP problems:
    t => t2virt%elm4
    g => gvirt%elm4
    f => MyFragment%ppfock

    call mem_TurnONThread_Memory()
    !$OMP PARALLEL DEFAULT(none) PRIVATE(multaibj,tmp,j,b,i,a,k,occ_tmp) &
    !$OMP SHARED(t,g,f,noccAOS,nvirtEOS,prefac_coul,prefac_k,myfragment, &
    !$OMP Performlag) REDUCTION(+:Evirt,lag_virt)
    call init_threadmemvar()

    ! Contributions from each individual occupied orbital
    call mem_alloc(occ_tmp,noccAOS)
    occ_tmp = 0.0E0_realk

    !$OMP DO SCHEDULE(dynamic,1) COLLAPSE(3)
    do j=1,noccAOS
       do b=1,nvirtEOS
          do i=1,noccAOS
             do a=1,nvirtEOS


                ! Contribution 3
                ! --------------

                ! Multiplier (multiplied by one half)
                multaibj = prefac_coul*t(a,i,b,j) -prefac_k*t(b,i,a,j)


                ! Energy contribution for orbitals (j,b,i,a)
                tmp = multaibj*g(a,i,b,j)
                ! Update total atomic fragment energy contribution 3
                Evirt = Evirt + tmp

                ! Update contribution from orbital i
                occ_tmp(i) = occ_tmp(i) + abs(tmp)
                ! Update contribution from orbital j (only if different from i to avoid double counting)
                if(i/=j) occ_tmp(j) = occ_tmp(j) + abs(tmp)

                ! Contribution 4
                ! --------------

                ! Skip contribution 4 for anything but MP2 and RIMP2
                if(PerformLag) then

                   do k=1,noccAOS

                      tmp =  t(a,k,b,j)*f(k,i) + t(a,i,b,k)*f(k,j)
                      tmp = -multaibj*tmp

                      ! Update total atomic fragment energy contribution 4
                      lag_virt = lag_virt + tmp

                      ! Update contribution from orbital i
                      occ_tmp(i) = occ_tmp(i) + tmp
                      ! Update contribution from orbital j (only if different from i to avoid double counting)
                      !if(i/=j) occ_tmp(j) = occ_tmp(j) + tmp
                      ! Update contribution from orbital k (only if different from i and j)
                      !if( (i/=k) .and. (j/=k) ) occ_tmp(k) = occ_tmp(k) + tmp

                   end do

                end if

             end do
          end do
       end do
    end do
    !$OMP END DO NOWAIT

    !$OMP CRITICAL
    ! Update total occupied contributions to fragment energy
    do i=1,noccAOS
       MyFragment%OccContribs(i) = MyFragment%OccContribs(i) + occ_tmp(i)
    end do
    !$OMP END CRITICAL

    call mem_dealloc(occ_tmp)
    call collect_thread_memory()
    !$OMP END PARALLEL
    call mem_TurnOffThread_Memory()

    t => null()
    g => null()
    f => null()


    ! Total atomic fragment energy
    ! ****************************
    ! Lagrangian energy only implemented for MP2 and RIMP2 so it gets special treatment
    if(MyFragment%ccmodel==MODEL_MP2) then
       MyFragment%energies(FRAGMODEL_LAGMP2) = Eocc + lag_occ + Evirt + lag_virt
    elseif(MyFragment%ccmodel==MODEL_RIMP2) then
       MyFragment%energies(FRAGMODEL_LAGRIMP2) = Eocc + lag_occ + Evirt + lag_virt
    end if

    ! Put occupied (Eocc) and virtual (Evirt) scheme energies into fragment energies array
    call put_fragment_energy_contribs_wrapper(Eocc,Evirt,MyFragment,SOS)

    ! Print out contributions
    ! ***********************

    if( myfragment%isopt )then
       write(DECinfo%output,*)
       write(DECinfo%output,*)
       write(DECinfo%output,*) '**********************************************************************'
       write(DECinfo%output,'(1X,a,i7)') 'Energy summary for fragment: ', &
          & MyFragment%EOSatoms(1)
       write(DECinfo%output,*) '**********************************************************************'
       write(DECinfo%output,'(1X,a,g20.10)') 'Single occupied energy = ', Eocc
       write(DECinfo%output,'(1X,a,g20.10)') 'Single virtual  energy = ', Evirt
       write(DECinfo%output,'(1X,a,g20.10)') 'Single Lagrangian occ term  = ', lag_occ
       write(DECinfo%output,'(1X,a,g20.10)') 'Single Lagrangian virt term = ', lag_virt

       write(DECinfo%output,*)
       write(DECinfo%output,*)
    endif


    call LSTIMER('START',tcpu2,twall2,DECinfo%output)
    call LSTIMER('L.ENERGY CONTR',tcpu,twall,DECinfo%output)


  end subroutine get_atomic_fragment_energy


  !> \brief  Contract amplitudes, and (T) intermediates to get 4th order contributions
  !          to single fragment occupied and virtual energies.
  !> \author Pablo Baudin (Based on Kasper's routine)
  !> \date   May 2015
  subroutine get_single_fragment_pT4_energy(gocc,gvirt,t2occ,t2virt,MyFragment)

    implicit none
    !> (T) intermediates T_aibj, only occ orbitals on central atom, virt AOS orbitals
    type(tensor), intent(in) :: gocc
    !> (T) intermediates T_aibj, only virt orbitals on central atom, occ AOS orbitals
    type(tensor), intent(in) :: gvirt
    !> doubles amplitudes, only occ orbitals on central atom, virt AOS orbitals
    type(tensor), intent(in) :: t2occ
    !> doubles amplitudes, only virt orbitals on central atom, occ AOS orbitals
    type(tensor), intent(in) :: t2virt
    !> Atomic fragment 
    type(decfrag), intent(inout) :: myfragment

    real(realk) :: Eocc,Evirt,tmp,multaibj
    real(realk) :: prefac_coul,prefac_k
    real(realk),pointer :: occ_tmp(:),virt_tmp(:)
    real(realk),pointer :: t(:,:,:,:), g(:,:,:,:), f(:,:)
    logical ::  something_wrong
    integer :: noccEOS,nvirtEOS,noccAOS,nvirtAOS
    integer :: i,j,k,a,b,c

    ! MyFragment%OccContribs and MyFragment%VirtContribs contains the contributions from
    ! each individual occupied and virtual orbital -- e.g. MyFragment%OccContribs(i)
    ! is the estimated change in the energy if occupied orbital with index
    ! MyFragment%occAOSidx(i) is removed from the fragment.

    ! Init stuff
    noccEOS  = MyFragment%noccEOS
    nvirtEOS = MyFragment%nvirtEOS
    noccAOS  = MyFragment%noccAOS
    nvirtAOS = MyFragment%nvirtAOS
    Eocc     = 0E0_realk
    Evirt    = 0E0_realk
    prefac_coul=2._realk
    prefac_k=1._realk


    ! init energy reals to be on the safe side
    ! note: OccEnergyPT and VirtEnergyPT is also initialized from in here
    !       as this (e4) routine is called before the e5 routine
    MyFragment%energies(FRAGMODEL_OCCpT) = 0.0E0_realk
    MyFragment%energies(FRAGMODEL_VIRTpT) = 0.0E0_realk
    MyFragment%energies(FRAGMODEL_OCCpT4) = 0.0E0_realk
    MyFragment%energies(FRAGMODEL_VIRTpT4) = 0.0E0_realk

    ! Sanity checks
    ! *************
    something_wrong=.false.
    if(t2occ%dims(1) /= nvirtAOS) something_wrong=.true.
    if(t2occ%dims(2) /= noccEOS) something_wrong=.true.
    if(t2occ%dims(3) /= nvirtAOS) something_wrong=.true.
    if(t2occ%dims(4) /= noccEOS) something_wrong=.true.
    
    if(gocc%dims(1) /= nvirtAOS) something_wrong=.true.
    if(gocc%dims(2) /= noccEOS) something_wrong=.true.
    if(gocc%dims(3) /= nvirtAOS) something_wrong=.true.
    if(gocc%dims(4) /= noccEOS) something_wrong=.true.

    if(t2virt%dims(1) /= nvirtEOS) something_wrong=.true.
    if(t2virt%dims(2) /= noccAOS) something_wrong=.true.
    if(t2virt%dims(3) /= nvirtEOS) something_wrong=.true.
    if(t2virt%dims(4) /= noccAOS) something_wrong=.true.
    
    if(gvirt%dims(1) /= nvirtEOS) something_wrong=.true.
    if(gvirt%dims(2) /= noccAOS) something_wrong=.true.
    if(gvirt%dims(3) /= nvirtEOS) something_wrong=.true.
    if(gvirt%dims(4) /= noccAOS) something_wrong=.true.

    if(something_wrong) then
       print *, 't2occ%dims =', t2occ%dims
       print *, 'gocc%dims  =', gocc%dims
       print *, 't2virt%dims=', t2virt%dims
       print *, 'gvirt%dims =', gvirt%dims
       print *, 'noccEOS  = ', noccEOS
       print *, 'noccAOS  = ', noccAOS
       print *, 'nvirtEOS = ', nvirtEOS
       print *, 'nvirtAOS = ', nvirtAOS
       call lsquit('get_single_fragment_pT_energy: &
            & Input dimensions do not match!',-1)
    end if


    ! Calculate Eocc
    ! **************

    ! Point to tensor structure to avoid OMP problems:
    t => t2occ%elm4
    g => gocc%elm4
    f => MyFragment%qqfock

    call mem_TurnONThread_Memory()
    !$OMP PARALLEL DEFAULT(none) PRIVATE(tmp,j,b,i,a,c,virt_tmp) &
    !$OMP SHARED(t,g,f,nvirtAOS,noccEOS,prefac_coul,prefac_k,myfragment) &
    !$OMP REDUCTION(+:Eocc)
    call init_threadmemvar()

    ! Contributions from each individual virtual orbital
    call mem_alloc(virt_tmp,nvirtAOS)
    virt_tmp = 0.0E0_realk       

    !$OMP DO SCHEDULE(dynamic,1) COLLAPSE(3)
    do j=1,noccEOS
       do b=1,nvirtAOS
          do i=1,noccEOS
             do a=1,nvirtAOS

                ! Energy contribution for orbitals (j,b,i,a)
                tmp = 2.0E0_realk*t(a,i,b,j)*(prefac_coul*g(a,i,b,j) - prefac_k*g(b,i,a,j))

                ! Update total atomic fragment energy contribution 1
                Eocc = Eocc + tmp

                ! Update contribution from orbital a
                virt_tmp(a) = virt_tmp(a) + abs(tmp)
                ! Update contribution from orbital b (only if different from a to avoid double counting)
                if(a/=b) virt_tmp(b) = virt_tmp(b) + abs(tmp)

             end do
          end do
       end do
    end do
    !$OMP END DO NOWAIT

    ! Update total virtual contributions to fragment energy
    !$OMP CRITICAL
    do a=1,nvirtAOS
       MyFragment%VirtContribs(a) = MyFragment%VirtContribs(a) + virt_tmp(a)
    end do
    !$OMP END CRITICAL

    call mem_dealloc(virt_tmp)
    call collect_thread_memory()
    !$OMP END PARALLEL
    call mem_TurnOffThread_Memory()


    ! Calculate Evirt
    ! ***************

    ! Point to tensor structure to avoid OMP problems:
    t => t2virt%elm4
    g => gvirt%elm4
    f => MyFragment%ppfock

    call mem_TurnONThread_Memory()
    !$OMP PARALLEL DEFAULT(none) PRIVATE(multaibj,tmp,j,b,i,a,k,occ_tmp) &
    !$OMP SHARED(t,g,f,noccAOS,nvirtEOS,prefac_coul,prefac_k,myfragment) &
    !$OMP REDUCTION(+:Evirt)
    call init_threadmemvar()

    ! Contributions from each individual occupied orbital
    call mem_alloc(occ_tmp,noccAOS)
    occ_tmp = 0.0E0_realk

    !$OMP DO SCHEDULE(dynamic,1) COLLAPSE(3)
    do j=1,noccAOS
       do b=1,nvirtEOS
          do i=1,noccAOS
             do a=1,nvirtEOS

                ! Multiplier (multiplied by one half)
                multaibj = prefac_coul*t(a,i,b,j) -prefac_k*t(b,i,a,j)

                ! Energy contribution for orbitals (j,b,i,a)
                tmp = 2.0E0_realk*multaibj*g(a,i,b,j)
                ! Update total atomic fragment energy contribution 3
                Evirt = Evirt + tmp

                ! Update contribution from orbital i
                occ_tmp(i) = occ_tmp(i) + abs(tmp)
                ! Update contribution from orbital j (only if different from i to avoid double counting)
                if(i/=j) occ_tmp(j) = occ_tmp(j) + abs(tmp)

             end do
          end do
       end do
    end do
    !$OMP END DO NOWAIT

    !$OMP CRITICAL
    ! Update total occupied contributions to fragment energy
    do i=1,noccAOS
       MyFragment%OccContribs(i) = MyFragment%OccContribs(i) + occ_tmp(i)
    end do
    !$OMP END CRITICAL

    call mem_dealloc(occ_tmp)
    call collect_thread_memory()
    !$OMP END PARALLEL
    call mem_TurnOffThread_Memory()

    t => null()
    g => null()
    f => null()


    ! get total fourth--order energy contribution
    MyFragment%energies(FRAGMODEL_OCCpT4) = Eocc
    MyFragment%energies(FRAGMODEL_VIRTpT4) = Evirt
   
    MyFragment%energies(FRAGMODEL_OCCpT) = MyFragment%energies(FRAGMODEL_OCCpT) &
       & + MyFragment%energies(FRAGMODEL_OCCpT4)
    MyFragment%energies(FRAGMODEL_VIRTpT) = MyFragment%energies(FRAGMODEL_VIRTpT) &
       & + MyFragment%energies(FRAGMODEL_VIRTpT4)


  end subroutine get_single_fragment_pT4_energy


  !> Purpose: Contract amplitudes and integrals to calculate DECNP fragment 
  !           occupied and virtual energy:
  !
  !> Author:  Pablo Baudin
  !> Date:    Feb. 2015
  subroutine get_decnp_fragment_energy(MyFragment,gocc,gvirt,t2occ,t2virt,pT_contrib)

     implicit none

     !> Atomic fragment
     type(decfrag), intent(inout) :: MyFragment
     !> Energy integrals (ai|bj)
     type(tensor), intent(inout) :: gocc
     type(tensor), intent(inout) :: gvirt
     !> double amplitudes from CC calculation (MP2, CCSD...)
     !  singles might be included as:
     !  t2(a,i,b,j) := t2(a,i,b,j) + t1(a,i)*t1(b,j)
     type(tensor), intent(inout) :: t2occ
     type(tensor), intent(inout) :: t2virt
     logical, intent(in), optional :: pT_contrib

     real(realk), pointer :: occ_tmp(:), virt_tmp(:)
     real(realk), pointer :: t(:,:,:,:), g(:,:,:,:)

     real(realk) :: tcpu1, twall1, tcpu2,twall2, tcpu,twall
     real(realk) :: Eocc, Evirt
     real(realk) :: prefac_coul, prefac_k, tmp
     logical ::  something_wrong, dopT
     integer :: noccEOS, nvirtEOS, noccAOS, nvirtAOS
     integer :: i,j,k,a,b

     ! Get occupied and virtual DECNP energy contributions:
     ! ----------------------------------------------------
     !
     ! Eocc  = sum_{ijab} t_{ij}^{ab} [ 2g_{aibj} - g_{biaj} ]
     ! where index i is restricted to the occupied EOS space,
     ! index j is restricted to the occupied AOS space,
     ! and indices a and b are restricted to the virtual AOS space
     !
     ! Evirt = sum_{ijab} t_{ij}^{ab} [ 2g_{aibj} - g_{biaj} ]
     ! where index a is restricted to the virtual EOS space,
     ! indices i and j are restricted to the occupied AOS space,
     ! and index b is restricted to the virtual AOS space
     ! 
     ! MyFragment%OccContribs and MyFragment%VirtContribs contains the contributions from
     ! each individual occupied and virtual orbital -- e.g. MyFragment%OccContribs(i)
     ! is the estimated change in the energy if occupied orbital with index
     ! MyFragment%occAOSidx(i) is removed from the fragment.
      
     call LSTIMER('START',tcpu,twall,DECinfo%output)
     call LSTIMER('START',tcpu1,twall1,DECinfo%output)
      
     ! Init stuff
     noccEOS  = MyFragment%noccEOS
     nvirtEOS = MyFragment%nvirtEOS
     noccAOS  = MyFragment%noccAOS
     nvirtAOS = MyFragment%nvirtAOS
     Eocc     = 0E0_realk
     Evirt    = 0E0_realk

     ! Just in case, zero individual orbital contributions for fragment
     MyFragment%OccContribs=0E0_realk
     MyFragment%VirtContribs=0E0_realk
      
     ! are we calculating the (T) contributions?
     dopT = .false.
     if (present(pT_contrib)) dopT=pT_contrib

     prefac_coul=2._realk
     prefac_k=1._realk
      
      
     ! Sanity checks
     ! *************
     something_wrong=.false.
     if(t2occ%dims(1) /= nvirtAOS) something_wrong=.true.
     if(t2occ%dims(2) /= noccEOS) something_wrong=.true.
     if(t2occ%dims(3) /= nvirtAOS) something_wrong=.true.
     if(t2occ%dims(4) /= noccAOS) something_wrong=.true.
     
     if(gocc%dims(1) /= nvirtAOS) something_wrong=.true.
     if(gocc%dims(2) /= noccEOS) something_wrong=.true.
     if(gocc%dims(3) /= nvirtAOS) something_wrong=.true.
     if(gocc%dims(4) /= noccAOS) something_wrong=.true.
      
     if(t2virt%dims(1) /= nvirtEOS) something_wrong=.true.
     if(t2virt%dims(2) /= noccAOS) something_wrong=.true.
     if(t2virt%dims(3) /= nvirtAOS) something_wrong=.true.
     if(t2virt%dims(4) /= noccAOS) something_wrong=.true.
     
     if(gvirt%dims(1) /= nvirtEOS) something_wrong=.true.
     if(gvirt%dims(2) /= noccAOS) something_wrong=.true.
     if(gvirt%dims(3) /= nvirtAOS) something_wrong=.true.
     if(gvirt%dims(4) /= noccAOS) something_wrong=.true.
      
     if(something_wrong) then
        print *, 't2occ%dims =', t2occ%dims
        print *, 'gocc%dims  =', gocc%dims
        print *, 't2virt%dims=', t2virt%dims
        print *, 'gvirt%dims =', gvirt%dims
        print *, 'noccEOS  = ', noccEOS
        print *, 'noccAOS  = ', noccAOS
        print *, 'nvirtEOS = ', nvirtEOS
        print *, 'nvirtAOS = ', nvirtAOS
        call lsquit('get_decnp_fragment_energy: &
             & Input dimensions do not match!',-1)
     end if
      
     ! Point to tensor structure to avoid OMP problems:
     t => t2occ%elm4
     g => gocc%elm4

     call mem_TurnONThread_Memory()
     !$OMP PARALLEL DEFAULT(none) PRIVATE(tmp,j,b,i,a,virt_tmp) &
     !$OMP SHARED(t,g,noccAOS,nvirtAOS,noccEOS,prefac_coul,prefac_k,myfragment) &
     !$OMP REDUCTION(+:Eocc)
     call init_threadmemvar()
     ! Contributions from each individual virtual orbital
     call mem_alloc(virt_tmp,nvirtAOS)
     virt_tmp = 0.0E0_realk       
      
     ! Calculate Eocc
     ! **************
     !$OMP DO SCHEDULE(dynamic,1) COLLAPSE(3)
     do j=1,noccAOS
        do b=1,nvirtAOS
           do i=1,noccEOS ! EOS index
              do a=1,nvirtAOS
      
      
                 ! Contribution 1
                 ! --------------
      
                 ! Energy contribution for orbitals (j,b,i,a)
                 tmp = t(a,i,b,j)*(prefac_coul*g(a,i,b,j) -prefac_k*g(b,i,a,j))
      
                 ! Update total atomic fragment energy contribution 1
                 Eocc = Eocc + tmp
      
                 ! Update contribution from orbital a
                 virt_tmp(a) = virt_tmp(a) + tmp
                 ! Update contribution from orbital b (only if different from a to avoid double counting)
                 if(a/=b) virt_tmp(b) = virt_tmp(b) + tmp
      
      
              end do
           end do
        end do
     end do
     !$OMP END DO NOWAIT
      
     ! Update total virtual contributions to fragment energy
     !$OMP CRITICAL
     do a=1,nvirtAOS
        MyFragment%VirtContribs(a) =MyFragment%VirtContribs(a) + virt_tmp(a)
     end do
     !$OMP END CRITICAL
      
     call mem_dealloc(virt_tmp)
     call collect_thread_memory()
     !$OMP END PARALLEL
     call mem_TurnOffThread_Memory()


     ! Point to tensor structure to avoid OMP problems:
     t => t2virt%elm4
     g => gvirt%elm4

     call mem_TurnONThread_Memory()
     !$OMP PARALLEL DEFAULT(none) PRIVATE(tmp,j,b,i,a,occ_tmp) &
     !$OMP SHARED(t,g,noccAOS,nvirtAOS,nvirtEOS,prefac_coul,prefac_k,myfragment) &
     !$OMP REDUCTION(+:Evirt)
     call init_threadmemvar()
     ! Contributions from each individual occupied orbital
     call mem_alloc(occ_tmp,noccAOS)
     occ_tmp = 0.0E0_realk

     ! Calculate Evirt
     ! ***************
     !$OMP DO SCHEDULE(dynamic,1) COLLAPSE(3)
     do j=1,noccAOS
        do b=1,nvirtAOS
           do i=1,noccAOS
              do a=1,nvirtEOS ! EOS index


                 ! Contribution 2
                 ! --------------

                 ! Energy contribution for orbitals (j,b,i,a)
                 tmp = (prefac_coul*t(a,i,b,j) -prefac_k*t(a,j,b,i))*g(a,i,b,j)

                 ! Update total atomic fragment energy contribution 3
                 Evirt = Evirt + tmp

                 ! Update contribution from orbital i
                 occ_tmp(i) = occ_tmp(i) + tmp
                 ! Update contribution from orbital j (only if different from i to avoid double counting)
                 if(i/=j) occ_tmp(j) = occ_tmp(j) + tmp


              end do
           end do
        end do
     end do
     !$OMP END DO NOWAIT

     !$OMP CRITICAL
     ! Update total occupied contributions to fragment energy
     do i=1,noccAOS
        MyFragment%OccContribs(i) = MyFragment%OccContribs(i) + occ_tmp(i)
     end do
     !$OMP END CRITICAL

     call mem_dealloc(occ_tmp)
     call collect_thread_memory()
     !$OMP END PARALLEL
     call mem_TurnOffThread_Memory()

     t => null() 
     g => null() 
      
     ! Total atomic fragment energy
     ! ****************************
     if (dopT) then
        MyFragment%energies(FRAGMODEL_OCCpT4) = 2.0e0_realk*Eocc
        MyFragment%energies(FRAGMODEL_OCCpT) = MyFragment%energies(FRAGMODEL_OCCpT) &
           &+ MyFragment%energies(FRAGMODEL_OCCpT4)

        MyFragment%energies(FRAGMODEL_VIRTpT4) = 2.0e0_realk*Evirt
        MyFragment%energies(FRAGMODEL_VIRTpT) = MyFragment%energies(FRAGMODEL_VIRTpT) &
           &+ MyFragment%energies(FRAGMODEL_VIRTpT4)
     else
        call put_fragment_energy_contribs_wrapper(Eocc,Evirt,MyFragment)
     end if


     ! Print out contributions
     ! ***********************
     if( myfragment%isopt )then
        write(DECinfo%output,*)
        write(DECinfo%output,*)
        write(DECinfo%output,*) '**********************************************************************'
        write(DECinfo%output,'(1X,a,i7)') 'DECNP Energy summary for fragment: ', &
           & MyFragment%EOSatoms(1)
        write(DECinfo%output,*) '**********************************************************************'
        write(DECinfo%output,'(1X,a,g20.10)') 'Single occupied energy = ', Eocc
        write(DECinfo%output,'(1X,a,g20.10)') 'Single virtual  energy = ', Evirt
      
        write(DECinfo%output,*)
        write(DECinfo%output,*)
     endif
      
      
     call LSTIMER('START',tcpu2,twall2,DECinfo%output)
     call LSTIMER('DECNP.ENERGY CONTR',tcpu,twall,DECinfo%output)
     

  end subroutine get_decnp_fragment_energy


  !> \brief Wrapper for pair_driver where can attach existing full
  !> molecular singles amplitudes to the fragment structure and update new
  !> improved full molecular singles amplitudes by the calculated fragment singles amplitudes.
  !> \author Kasper Kristensen
  !> \date March 2012
  subroutine pair_driver_singles(natoms,nocc,nvirt,&
       & OccOrbitals,virtOrbitals,MyLsitem,MyMolecule,&
       & Fragment1,Fragment2,PairFragment,t1old,t1new)

    implicit none
    !> Number of atoms for full molecule
    integer, intent(in) :: natoms
    !> Number of occupied orbitals in full molecule
    integer,intent(in) :: nOcc
    !> Number of virtupied orbitals in full molecule
    integer,intent(in) :: nvirt
    !> Occupied orbitals for full molecule
    type(decorbital), intent(in) :: OccOrbitals(nocc)
    !> virtupied orbitals for full molecule
    type(decorbital), intent(in) :: virtOrbitals(nvirt)
    !> Full molecule info (will be unchanged at output)
    type(fullmolecule), intent(inout) :: MyMolecule
    !> LSDalton info
    type(lsitem), intent(inout) :: mylsitem
    !> Fragment 1 in the pair fragment
    type(decfrag),intent(inout) :: Fragment1
    !> Fragment 2 in the pair fragment
    type(decfrag),intent(inout) :: Fragment2
    !> Pair fragment
    type(decfrag), intent(inout) :: PairFragment
    !> Existing full molecular t1 amplitudes
    type(array2),intent(in) :: t1old
    !> New full molecular t1 amplitudes which will be updated with the contribution from pair fragment
    type(array2),intent(inout) :: t1new
    type(mp2grad) :: grad
    logical,pointer :: dopair(:,:)


    ! Attach singles amplitudes to pair fragment structure
    ! if long-range singles polarization effects are requested.
    call get_fragmentt1_AOSAOS_from_full(PairFragment,t1old)

    ! Call main driver to get energy
    call pair_driver(Fragment1,Fragment2,PairFragment,grad)

    ! Update full molecular singles amplitudes with (virt EOS,occ EOS) fragment contributions
    ! ***************************************************************************************

    ! Extract (virt EOS,occ EOS) indices from fragment
    call extract_specific_fragmentt1(PairFragment,.true.,.true.)

    ! Determine which orbital pairs to consider (avoid double counting)
    call mem_alloc(dopair,natoms,natoms)
    call which_pairs(Fragment1, Fragment2, natoms,dopair)

    ! Update full molecular singles amplitudes
    ! with pair (i,j) interaction contribution.
    call update_full_t1_from_pair_frag(PairFragment,nocc,nvirt,&
         & natoms,dopair,OccOrbitals,virtOrbitals,t1new)
    call mem_dealloc(dopair)

    ! Free t1 info again
    call free_fragment_t1(PairFragment)

  end subroutine pair_driver_singles


  !> \brief Driver for calculating pair interaction energy for a given fragment using
  !> the Lagrangian approach.
  !> If requested, first order properties (MP2 density or gradient) are also calculated and saved.
  !> \author Kasper Kristensen
  subroutine pair_driver(Fragment1,Fragment2,PairFragment,grad)

    implicit none
    !> Fragment 1 in the pair fragment
    type(decfrag),intent(in) :: Fragment1
    !> Fragment 2 in the pair fragment
    type(decfrag),intent(in) :: Fragment2
    !> Atomic fragment 
    type(decfrag), intent(inout) :: PairFragment
    !> MP2 gradient structure (only calculated if DECinfo%first_order is turned on)
    type(mp2grad),intent(inout) :: grad


    ! Run calculation using input fragment
    call fragment_energy_and_prop(PairFragment,Fragment1=Fragment1,Fragment2=Fragment2,&
         & grad=grad)

  end subroutine pair_driver


  !> \brief Contract amplitudes, multipliers, and integrals to calculate pair interaction
  !> Lagrangian energy. Heavily inspired by get_atomic_fragment_energy, but it is necessary
  !> to keep it separate for clairity.
  !> \author Kasper Kristensen
  !> \date August 2011
  subroutine get_pair_fragment_energy(gocc,gvirt,t2occ,t2virt,&
      & Fragment1, Fragment2, PairFragment,doSOS)


     implicit none
     !> Two-electron integrals (a i | b j), only occ orbitals on central atom, virt AOS orbitals
     type(tensor), intent(in) :: gocc
     !> Two-electron integrals (a i | b j), only virt orbitals on central atom, occ AOS orbitals
     type(tensor), intent(in) :: gvirt
     !> MP2 amplitudes, only occ orbitals on central atom, virt AOS orbitals
     type(tensor), intent(in) :: t2occ
     !> MP2 amplitudes, only virt orbitals on central atom, occ AOS orbitals
     type(tensor), intent(in) :: t2virt
     !> Fragment 1 in the pair fragment
     type(decfrag),intent(in) :: Fragment1
     !> Fragment 2 in the pair fragment
     type(decfrag),intent(in) :: Fragment2
     !> Pair fragment formed from fragment 1 and 2
     type(decfrag), intent(inout) :: PairFragment
    !>SOSEX cont
     logical,intent(in),optional :: doSOS

     real(realk) :: tcpu,twall,tcpu1,tcpu2,twall1,twall2
     real(realk) :: pairdist,tmp,multaibj,prefac_coul,prefac_k
     real(realk) :: Eocc,lag_occ,Evirt,lag_virt
     real(realk),pointer :: t(:,:,:,:), g(:,:,:,:), f(:,:)
     logical,pointer :: dopair_occ(:,:), dopair_virt(:,:)
     logical :: something_wrong, do_non_pdm,SOS,PerformLag
     integer :: noccEOS,nvirtEOS,noccAOS,nvirtAOS
     integer :: i,j,k,a,b,c

     PerformLag = pairfragment%ccmodel==MODEL_MP2 .OR. pairfragment%ccmodel==MODEL_RIMP2
     ! Pair interaction Lagrangian energy can be split into four contributions:
     ! The first two (Eocc and lag_occ) use occupied EOS orbitals and virtual AOS orbitals.
     ! The last two (Evirt and lag_virt) use virtual EOS orbitals and occupied AOS orbitals.
     !
     ! With the restrictions above the Lagrangian energy is given by:
     ! energy  = Eocc + lag_occ + Evirt + lag_virt
     !
     ! Eocc     = sum_{ijab} t_{ij}^{ab} [ 2g_{aibj} - g_{biaj} ]
     ! lag_occ  = 1/2 sum_{ijabc} mult_{ij}^{ab} [ t_{ij}^{cb} F_{ac} + t_{ij}^{ac} F_{bc} ]
     ! Evirt    = 1/2 sum_{ijab} mult_{ij}^{ab} g_{aibj}
     ! lag_virt = - 1/2 sum_{ijkab} mult_{ij}^{ab} [ t_{kj}^{ab} F_{ki} + t_{ik}^{ab} F_{kj} ]
     !
     ! Additional restriction to avoid double counting:
     ! Eocc  and lag_occ:  i and j must belong to two different atoms,
     !                     e.g., i belongs to atom1 and j belongs to atom2 - or vice versa
     ! Evirt and lag_virt: a and b must belong to two different atoms,
     !                     e.g., a belongs to atom1 and b belongs to atom2 - or vice versa
     !
     ! Important: For MP2 the multipliers can be determined directly from the amplitudes
     !            using the simple relation:
     !
     ! mult_{ij}^{ab} = 4*t_{ij}^{ab} - 2*t_{ij}^{ba}
     !
     ! In the calculations below the factor 1/2 in the equations
     ! above is absorbed in the multipliers.

     call LSTIMER('START',tcpu,twall,DECinfo%output)
     call LSTIMER('START',tcpu1,twall1,DECinfo%output)

     ! Init stuff
     noccEOS  = PairFragment%noccEOS
     nvirtEOS = PairFragment%nvirtEOS
     noccAOS  = PairFragment%noccAOS
     nvirtAOS = PairFragment%nvirtAOS
     Eocc     = 0E0_realk
     lag_occ  = 0E0_realk
     Evirt    = 0E0_realk
     lag_virt = 0E0_realk
     SOS      = .false.

     if(present(doSOS)) SOS = doSOS


     ! Distance between fragments in Angstrom
     pairdist = bohr_to_angstrom*PairFragment%pairdist
     if(PairFragment%ccmodel==MODEL_RPA) then
       prefac_coul = 1._realk
       prefac_k=0.0_realk
       if(SOS) prefac_k=0.5_realk
     elseif(PairFragment%ccmodel==MODEL_SOSEX) then
       prefac_coul = 1._realk
       prefac_k=0.5_realk
     else
        prefac_coul =2._realk
        prefac_k = 1._realk
     endif

     ! Which "interaction pairs" to include for occ and virt space (avoid double counting)
     call mem_alloc(dopair_occ,noccEOS,noccEOS)
     call mem_alloc(dopair_virt,nvirtEOS,nvirtEOS)
     call which_pairs_occ(Fragment1,Fragment2,PairFragment,dopair_occ)
     call which_pairs_virt(Fragment1,Fragment2,PairFragment,dopair_virt)

     ! Sanity checks
     ! *************
     something_wrong=.false.
     IF(.NOT.DECinfo%OnlyVirtPart)THEN
        if(t2occ%dims(1) /= nvirtAOS) something_wrong = .true.
        if(t2occ%dims(2) /= noccEOS)  something_wrong = .true.
        if(t2occ%dims(3) /= nvirtAOS) something_wrong = .true.
        if(t2occ%dims(4) /= noccEOS)  something_wrong = .true.

        if(gocc%dims(1) /= nvirtAOS) something_wrong = .true.
        if(gocc%dims(2) /= noccEOS)  something_wrong = .true.
        if(gocc%dims(3) /= nvirtAOS) something_wrong = .true.
        if(gocc%dims(4) /= noccEOS)  something_wrong = .true.
     ENDIF

     IF(.NOT.DECinfo%OnlyOccPart)THEN
        if(t2virt%dims(1) /= nvirtEOS) something_wrong = .true.
        if(t2virt%dims(2) /= noccAOS)  something_wrong = .true.
        if(t2virt%dims(3) /= nvirtEOS) something_wrong = .true.
        if(t2virt%dims(4) /= noccAOS)  something_wrong = .true.

        if(gvirt%dims(1) /= nvirtEOS) something_wrong = .true.
        if(gvirt%dims(2) /= noccAOS)  something_wrong = .true.
        if(gvirt%dims(3) /= nvirtEOS) something_wrong = .true.
        if(gvirt%dims(4) /= noccAOS)  something_wrong = .true.
     ENDIF
     if(something_wrong) then
        print *, 't2occ%dims =', t2occ%dims
        print *, 'gocc%dims  =', gocc%dims
        print *, 't2virt%dims=', t2virt%dims
        print *, 'gvirt%dims =', gvirt%dims
        print *, 'noccEOS  = ', noccEOS
        print *, 'noccAOS  = ', noccAOS
        print *, 'nvirtEOS = ', nvirtEOS
        print *, 'nvirtAOS = ', nvirtAOS
        call lsquit('get_pair_fragment_energy: Input dimensions do not match!',-1)
     end if

     do_non_pdm = .false.
     if( .not.DECinfo%OnlyVIRTPart )then
        do_non_pdm = do_non_pdm .or. (t2occ%itype == TT_DENSE .and. gocc%itype == TT_DENSE )
     endif
     if(.not.DECinfo%OnlyoccPart)then
        do_non_pdm = do_non_pdm .or. (t2virt%itype ==  TT_DENSE .and. gvirt%itype == TT_DENSE)
     endif

     if( do_non_pdm )then

        if(.not. DECinfo%onlyVirtpart) then

           ! Calculate Eocc and lag_occ
           ! **************************

           ! Point to tensor structure to avoid OMP problems:
           t => t2occ%elm4
           g => gocc%elm4
           f => PairFragment%qqfock

           call mem_TurnONThread_Memory()
           !$OMP PARALLEL DEFAULT(none) PRIVATE(multaibj,tmp,j,b,i,a,c) &
           !$OMP SHARED(t,g,f,nvirtAOS,noccEOS,prefac_coul,prefac_k,pairfragment, &
           !$OMP PerformLag,dopair_occ) REDUCTION(+:Eocc,lag_occ)
           call init_threadmemvar()

           !$OMP DO SCHEDULE(dynamic,1) COLLAPSE(3)
           do j=1,noccEOS
              do b=1,nvirtAOS
                 do i=1,noccEOS

                    ! Only update for "interaction orbital pairs" - see which_pairs_occ
                    if( dopair_occ(i,j) ) then  !DoPair1and2

                       do a=1,nvirtAOS

                          ! Update pair interaction energy contribution 1
                          Eocc = Eocc + t(a,i,b,j)*(prefac_coul*g(a,i,b,j) -prefac_k*g(b,i,a,j))


                          ! Skip contribution 2 for anything but MP2 and RIMP2
                          if(PerformLag) then

                             ! Multiplier (multiplied by one half)
                             multaibj = prefac_coul*t(a,i,b,j) - prefac_k*t(b,i,a,j)

                             tmp = 0E0_realk
                             do c=1,nvirtAOS
                                tmp = tmp + t(c,i,b,j)*f(c,a) + t(a,i,c,j)*f(c,b)
                             end do

                             ! Update pair interaction energy contribution 2
                             lag_occ = lag_occ + multaibj*tmp

                          end if


                       end do

                    end if

                 end do
              end do
           end do
           !$OMP END DO NOWAIT
           call collect_thread_memory()
           !$OMP END PARALLEL
           call mem_TurnOffThread_Memory()
        ELSE
           Eocc    = 0.0E0_realk
           lag_occ = 0.0E0_realk
        ENDIF

        if(.not. DECinfo%onlyoccpart) then

           ! Calculate Evirt and lag_virt
           ! ****************************

           ! Point to tensor structure to avoid OMP problems:
           t => t2virt%elm4
           g => gvirt%elm4
           f => PairFragment%ppfock

           call mem_TurnONThread_Memory()
           !$OMP PARALLEL DEFAULT(none) PRIVATE(multaibj,tmp,j,b,i,a,k) &
           !$OMP SHARED(t,g,f,noccAOS,nvirtEOS,prefac_coul,prefac_k,pairfragment, &
           !$OMP Performlag,dopair_virt) REDUCTION(+:Evirt,lag_virt)
           call init_threadmemvar()

           !$OMP DO SCHEDULE(dynamic,1) COLLAPSE(3)
           do j=1,noccAOS
              do b=1,nvirtEOS
                 do i=1,noccAOS
                    do a=1,nvirtEOS

                       ! Only update for "interaction orbital pairs" - see which_pairs_virt
                       if( dopair_virt(a,b) ) then !Dopair3and4

                          ! Multiplier (multiplied by one half)
                          multaibj = prefac_coul*t(a,i,b,j) - prefac_k*t(b,i,a,j)

                          ! Update total atomic fragment energy contribution 3
                          Evirt = Evirt + multaibj*g(a,i,b,j)


                          ! Skip contribution 4 for anything but MP2 and RIMP2
                          if(PerformLag) then

                             tmp=0E0_realk
                             do k=1,noccAOS
                                tmp = tmp + t(a,k,b,j)*f(k,i) + t(a,i,b,k)*f(k,j)
                             end do

                             ! Update pair interaction energy contribution 4
                             lag_virt = lag_virt - multaibj*tmp

                          end if

                       end if

                    end do
                 end do
              end do
           end do
           !$OMP END DO NOWAIT
           call collect_thread_memory()
           !$OMP END PARALLEL
           call mem_TurnOffThread_Memory()
        ELSE
           Evirt    = 0.0E0_realk
           lag_virt = 0.0E0_realk      
        ENDIF

     else
        call lsquit("ERROR(get_pair_fragment_energy) PDM version not yetimplemented",-1)
     endif

     ! Total pair interaction energy
     ! *****************************
     ! Lagrangian energy only implemented for MP2 so it gets special treatment
     if(PairFragment%ccmodel==MODEL_MP2) then
        PairFragment%energies(FRAGMODEL_LAGMP2) = Eocc + lag_occ + Evirt + lag_virt
     elseif(PairFragment%ccmodel==MODEL_RIMP2) then
        PairFragment%energies(FRAGMODEL_LAGRIMP2) = Eocc + lag_occ + Evirt + lag_virt
     end if

     ! Put occupied (Eocc) and virtual (Evirt) scheme energies into fragment energies array
     call put_fragment_energy_contribs_wrapper(Eocc,Evirt,PairFragment,SOS)

     call mem_dealloc(dopair_occ)
     call mem_dealloc(dopair_virt)

     ! Print out contributions
     ! ***********************

     write(DECinfo%output,*)
     write(DECinfo%output,*)
     write(DECinfo%output,*) '*****************************************************************************'
     write(DECinfo%output,'(1X,a,2i7)') 'Energy summary for pair fragment: ', &
        & Fragment1%EOSatoms(1), Fragment2%EOSatoms(1)
     write(DECinfo%output,*) '*****************************************************************************'

     if(.not. DECinfo%onlyvirtpart) then
        write(DECinfo%output,'(1X,a,g16.5,g20.10)') 'Distance(Ang), pair occ energy  = ', pairdist,Eocc
     endif
     if(.not. DECinfo%onlyoccpart) then
        write(DECinfo%output,'(1X,a,g16.5,g20.10)') 'Distance(Ang), pair virt energy = ', pairdist,Evirt
     endif
     if(.not. (DECinfo%onlyoccpart.or. DECinfo%onlyvirtpart)) then
        write(DECinfo%output,'(1X,a,g16.5,g20.10)') 'Distance(Ang), pair lagr. occ term  = ', pairdist,lag_occ
        write(DECinfo%output,'(1X,a,g16.5,g20.10)') 'Distance(Ang), pair lagr. virt term = ', pairdist,lag_virt
     end if
     write(DECinfo%output,*)
     write(DECinfo%output,*)

     call LSTIMER('L.ENERGY CONTR',tcpu,twall,DECinfo%output)
     call LSTIMER('START',tcpu2,twall2,DECinfo%output)

  end subroutine get_pair_fragment_energy


  !> \brief  Contract amplitudes, and (T) intermediates to get 4th order contributions
  !          to pair fragment occupied and virtual energies.
  !> \author Pablo Baudin (Based on Kasper's routine)
  !> \date   May 2015
  subroutine get_pair_fragment_pT4_energy(gocc,gvirt,t2occ,t2virt,&
      & Fragment1, Fragment2, PairFragment)


     implicit none
     !> (T) doubles intermediates, only occ orbitals on central atom, virt AOS orbitals
     type(tensor), intent(in) :: gocc
     !> (T) doubles intermediates, only virt orbitals on central atom, occ AOS orbitals
     type(tensor), intent(in) :: gvirt
     !> doubles amplitudes, only occ orbitals on central atom, virt AOS orbitals
     type(tensor), intent(in) :: t2occ
     !> doubles amplitudes, only virt orbitals on central atom, occ AOS orbitals
     type(tensor), intent(in) :: t2virt
     !> Fragment 1 in the pair fragment
     type(decfrag),intent(in) :: Fragment1
     !> Fragment 2 in the pair fragment
     type(decfrag),intent(in) :: Fragment2
     !> Pair fragment formed from fragment 1 and 2
     type(decfrag), intent(inout) :: PairFragment

     real(realk) :: multaibj,prefac_coul,prefac_k
     real(realk) :: Eocc,Evirt
     real(realk),pointer :: t(:,:,:,:), g(:,:,:,:), f(:,:)
     logical,pointer :: dopair_occ(:,:), dopair_virt(:,:)
     logical :: something_wrong, do_non_pdm
     integer :: noccEOS,nvirtEOS,noccAOS,nvirtAOS
     integer :: i,j,k,a,b,c


     ! Init stuff
     noccEOS  = PairFragment%noccEOS
     nvirtEOS = PairFragment%nvirtEOS
     noccAOS  = PairFragment%noccAOS
     nvirtAOS = PairFragment%nvirtAOS
     Eocc     = 0E0_realk
     Evirt    = 0E0_realk
     prefac_coul =2._realk
     prefac_k = 1._realk

     ! Which "interaction pairs" to include for occ and virt space (avoid double counting)
     call mem_alloc(dopair_occ,noccEOS,noccEOS)
     call mem_alloc(dopair_virt,nvirtEOS,nvirtEOS)
     call which_pairs_occ(Fragment1,Fragment2,PairFragment,dopair_occ)
     call which_pairs_virt(Fragment1,Fragment2,PairFragment,dopair_virt)

     ! Sanity checks
     ! *************
     something_wrong=.false.
     IF(.NOT.DECinfo%OnlyVirtPart)THEN
        if(t2occ%dims(1) /= nvirtAOS) something_wrong = .true.
        if(t2occ%dims(2) /= noccEOS)  something_wrong = .true.
        if(t2occ%dims(3) /= nvirtAOS) something_wrong = .true.
        if(t2occ%dims(4) /= noccEOS)  something_wrong = .true.

        if(gocc%dims(1) /= nvirtAOS) something_wrong = .true.
        if(gocc%dims(2) /= noccEOS)  something_wrong = .true.
        if(gocc%dims(3) /= nvirtAOS) something_wrong = .true.
        if(gocc%dims(4) /= noccEOS)  something_wrong = .true.
     ENDIF

     IF(.NOT.DECinfo%OnlyOccPart)THEN
        if(t2virt%dims(1) /= nvirtEOS) something_wrong = .true.
        if(t2virt%dims(2) /= noccAOS)  something_wrong = .true.
        if(t2virt%dims(3) /= nvirtEOS) something_wrong = .true.
        if(t2virt%dims(4) /= noccAOS)  something_wrong = .true.

        if(gvirt%dims(1) /= nvirtEOS) something_wrong = .true.
        if(gvirt%dims(2) /= noccAOS)  something_wrong = .true.
        if(gvirt%dims(3) /= nvirtEOS) something_wrong = .true.
        if(gvirt%dims(4) /= noccAOS)  something_wrong = .true.
     ENDIF
     if(something_wrong) then
        print *, 't2occ%dims =', t2occ%dims
        print *, 'gocc%dims  =', gocc%dims
        print *, 't2virt%dims=', t2virt%dims
        print *, 'gvirt%dims =', gvirt%dims
        print *, 'noccEOS  = ', noccEOS
        print *, 'noccAOS  = ', noccAOS
        print *, 'nvirtEOS = ', nvirtEOS
        print *, 'nvirtAOS = ', nvirtAOS
        call lsquit('get_pair_fragment_pT4_energy: Input dimensions do not match!',-1)
     end if

     do_non_pdm = .false.
     if( .not.DECinfo%OnlyVIRTPart )then
        do_non_pdm = do_non_pdm .or. (t2occ%itype == TT_DENSE .and. gocc%itype == TT_DENSE )
     endif
     if(.not.DECinfo%OnlyoccPart)then
        do_non_pdm = do_non_pdm .or. (t2virt%itype ==  TT_DENSE .and. gvirt%itype == TT_DENSE)
     endif

     if( do_non_pdm )then

        if(.not. DECinfo%onlyVirtpart) then

           ! Calculate Eocc
           ! **************

           ! Point to tensor structure to avoid OMP problems:
           t => t2occ%elm4
           g => gocc%elm4
           f => PairFragment%qqfock

           call mem_TurnONThread_Memory()
           !$OMP PARALLEL DEFAULT(none) PRIVATE(j,b,i,a,c) &
           !$OMP SHARED(t,g,f,nvirtAOS,noccEOS,prefac_coul,prefac_k,pairfragment, &
           !$OMP dopair_occ) REDUCTION(+:Eocc)
           call init_threadmemvar()

           !$OMP DO SCHEDULE(dynamic,1) COLLAPSE(3)
           do j=1,noccEOS
              do b=1,nvirtAOS
                 do i=1,noccEOS

                    ! Only update for "interaction orbital pairs" - see which_pairs_occ
                    if( dopair_occ(i,j) ) then  !DoPair1and2

                       do a=1,nvirtAOS

                          ! Update pair interaction energy contribution 1
                          Eocc = Eocc + 2.0e0_realk*t(a,i,b,j)*(prefac_coul*g(a,i,b,j) -prefac_k*g(b,i,a,j))

                       end do

                    end if

                 end do
              end do
           end do
           !$OMP END DO NOWAIT
           call collect_thread_memory()
           !$OMP END PARALLEL
           call mem_TurnOffThread_Memory()
        ELSE
           Eocc    = 0.0E0_realk
        ENDIF

        if(.not. DECinfo%onlyoccpart) then

           ! Calculate Evirt
           ! ***************

           ! Point to tensor structure to avoid OMP problems:
           t => t2virt%elm4
           g => gvirt%elm4
           f => PairFragment%ppfock

           call mem_TurnONThread_Memory()
           !$OMP PARALLEL DEFAULT(none) PRIVATE(multaibj,j,b,i,a,k) &
           !$OMP SHARED(t,g,f,noccAOS,nvirtEOS,prefac_coul,prefac_k,pairfragment, &
           !$OMP dopair_virt) REDUCTION(+:Evirt)
           call init_threadmemvar()

           !$OMP DO SCHEDULE(dynamic,1) COLLAPSE(3)
           do j=1,noccAOS
              do b=1,nvirtEOS
                 do i=1,noccAOS
                    do a=1,nvirtEOS

                       ! Only update for "interaction orbital pairs" - see which_pairs_virt
                       if( dopair_virt(a,b) ) then !Dopair3and4

                          ! Multiplier (multiplied by one half)
                          multaibj = prefac_coul*t(a,i,b,j) - prefac_k*t(b,i,a,j)

                          ! Update total atomic fragment energy contribution 3
                          Evirt = Evirt + 2.0e0_realk*multaibj*g(a,i,b,j)

                       end if

                    end do
                 end do
              end do
           end do
           !$OMP END DO NOWAIT
           call collect_thread_memory()
           !$OMP END PARALLEL
           call mem_TurnOffThread_Memory()
        ELSE
           Evirt    = 0.0E0_realk
        ENDIF

     else
        call lsquit("ERROR(get_pair_fragment_pT4_energy) PDM version not yetimplemented",-1)
     endif

     ! Put occupied (Eocc) and virtual (Evirt) scheme energies into fragment energies array
     PairFragment%energies(FRAGMODEL_OCCpT4) = Eocc
     PairFragment%energies(FRAGMODEL_OCCpT) = PairFragment%energies(FRAGMODEL_OCCpT) &
        &+ PairFragment%energies(FRAGMODEL_OCCpT4)

     PairFragment%energies(FRAGMODEL_VIRTpT4) = Evirt
     PairFragment%energies(FRAGMODEL_VIRTpT) = PairFragment%energies(FRAGMODEL_VIRTpT) &
        &+ PairFragment%energies(FRAGMODEL_VIRTpT4)


     call mem_dealloc(dopair_occ)
     call mem_dealloc(dopair_virt)

  end subroutine get_pair_fragment_pT4_energy


  !> \brief Print a simple "ascii-art" plot of the largest pair energies.
  !> \author Kasper Kristensen
  !> \date October 2012
  subroutine plot_pair_energies(natoms,paircut,FragEnergies,MyMolecule,dofrag)

    implicit none
    !> Number of atoms in molecule
    integer,intent(in) :: natoms
    !> Pair cut off distance used (in a.u.)
    real(realk),intent(in) :: paircut
    !> Fragment energies 
    real(realk),dimension(natoms,natoms),intent(in) :: FragEnergies
    !> Full molecule structure
    type(fullmolecule),intent(in) :: MyMolecule
    !> dofrag(P) is true if P has orbitals assigned
    logical,intent(in) :: dofrag(natoms)
    real(realk),pointer :: DistAng(:,:), xpoints(:), ypoints(:), xpoints2(:), ypoints2(:)
    integer :: mindist, maxdist, tmp, distInt,idx,npoints2
    real(realk) :: cutAng,endpoint,ln10
    integer :: i,j,npoints,k,interval,npairfrags
    character(len=16) :: xlabel, ylabel
    logical,pointer :: anypoints(:)

    ! Get number of pair fragments in relevant plotting interval
    npairfrags = get_num_of_pair_fragments(MyMolecule,dofrag,DECinfo%PairMinDist,&
         & DECinfo%pair_distance_threshold)
    ! Only relevant to plot at least two points
    if(npairfrags < 2) return


    ! Get distance table and pair cut off in Angstrom
    call mem_alloc(DistAng,natoms,natoms)
    DistAng = bohr_to_angstrom*MyMolecule%DistanceTable
    cutAng = bohr_to_angstrom*paircut

    ! Minimum and maximum pair distances
    mindist=1000
    maxdist=0
    iloop: do i=1,natoms
       if(.not. dofrag(i)) cycle iloop
       jloop: do j=i+1,natoms

          if(DistAng(i,j) > cutAng) cycle jloop
          if(.not. dofrag(j)) cycle jloop
          if(MyMolecule%ccmodel(i,j)==MODEL_NONE) cycle jloop

          

          ! Nearest integer smaller than actual pair distance +0.5Angstrom
          ! (this is most appropriate when we consider intervals of 1 Angstrom below,
          !  because we then ensure that we don't have any intervals with very few points)
          tmp = floor(DistAng(i,j) + 0.4999E0_realk)

          ! Minimum
          ! Safety precaution: Never consider pairs separated by less than 1 Ang
          if(mindist > tmp .and. tmp>0) mindist = tmp

          ! Max
          if(maxdist < tmp) maxdist = tmp

       end do jloop
    end do iloop


    ! We now consider the following intervals:
    !
    ! [mindist,mindist+1]               ! interval 1
    ! ]mindist+1,mindist+2]             ! interval 2
    ! ]mindist+2,mindist+3]             ! interval 3
    ! ...
    ! ]maxdist-1,maxdist]               ! interval maxdist-mindist

    ! Number of points
    npoints = maxdist - mindist

    ! Sanity check
    if(npoints==0) then
       write(DECinfo%output,*) 'WARNING: No pair points to plot!'
       call mem_dealloc(DistAng)
       return
    end if

    ! x points to plot (distances in Angstrom)
    call mem_alloc(xpoints,npoints) 
    ! xpoints(k) is set to be the beginning point of interval k
    ! (they will be reset to the middle of the interval before plotting)
    do k=1,npoints
       interval = mindist+k-1
       xpoints(k) = real(interval) 
    end do

    ! y points to plot (absolute pair interaction energies in a.u.)
    call mem_alloc(ypoints,npoints) 
    ypoints = 0.0E0_realk

    ! Keep track of whether there are in fact data points in each interval
    call mem_alloc(anypoints,npoints) 
    anypoints=.false.


    ! Determine ypoints(k) as maximum (absolute) pair energy in interval k
    do i=1,natoms
       if(.not. dofrag(i)) cycle
       do j=i+1,natoms

          if(DistAng(i,j) > real(maxdist)) cycle
          if(.not. dofrag(j)) cycle

          ! Determine which interval pair energy (i,j) belongs to
          interval=0
          IntervalLoop: do k=1,npoints
             endpoint = xpoints(k)+1.0E0_realk ! end point for interval k

             ! Is pair (i,j) in interval k?
             if(DistAng(i,j) > xpoints(k) .and. DistAng(i,j)<=endpoint) then
                interval=k
                anypoints(interval) = .true.
                exit IntervalLoop
             end if

          end do IntervalLoop

          ! Sanity check
          if(interval/=0) then

             ! Replace current energy in interval if pair energy (i,j) is larger
             if( abs(FragEnergies(i,j)) > ypoints(interval)) then
                ypoints(interval) = abs(FragEnergies(i,j))
             end if

          end if

       end do
    end do


    ! Remove empty intervals
    ! **********************
    npoints2 = count(anypoints) ! Number of points with actual data points

    ! Sanity check
    if(npoints2<2) then
       write(DECinfo%output,*) 'Skipping pair plot - less than 2 points to plot'
       call mem_dealloc(xpoints)
       call mem_dealloc(ypoints)
       call mem_dealloc(anypoints)
       call mem_dealloc(DistAng)
       return
    end if

    call mem_alloc(xpoints2,npoints2) 
    call mem_alloc(ypoints2,npoints2)
    idx=0
    do i=1,npoints  ! loop over number of intervals (including empty ones)
       if(anypoints(i)) then ! Put values into new x and y vectors of reduced size
          idx=idx+1
          xpoints2(idx) = xpoints(i)
          ypoints2(idx) = ypoints(i)
       end if
    end do
    ! From now on we work with xpoints2 and ypoints2


    ! Change x points to be in the middel of the interval
    do k=1,npoints2
       xpoints2(k) = xpoints2(k) + 0.5E0_realk
    end do

    ! Take 10-logarithm of x and y points
    ! (Note: Fortran-log is natural logarithm, ln, and
    !        log_10(x) = ln(a) / ln(10)
    ln10 = log(10.0E0_realk)
    do k=1,npoints2
       xpoints2(k) = log(xpoints2(k))/ln10
       ypoints2(k) = log(ypoints2(k))/ln10
    end do


    ! Make "ascii-art" plot
    ! *********************
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a)') '-------------------------------------------------------------'
    write(DECinfo%output,'(1X,a)') '                 PAIR INTERACTION ENERGY PLOT                '
    write(DECinfo%output,'(1X,a)') '-------------------------------------------------------------'
    write(DECinfo%output,'(1X,a)') 'Plot contains maximum pair energies in intervals of 1 Angstrom.'
    write(DECinfo%output,*)

    ! Set x and y labels
    xlabel = 'log[R/Ang]      '
    ylabel = 'log[|E|/Hartree]'
    call simple_ascii_plot(npoints2,xpoints2,ypoints2,xlabel,ylabel)

    call mem_dealloc(xpoints)
    call mem_dealloc(ypoints)
    call mem_dealloc(xpoints2)
    call mem_dealloc(ypoints2)
    call mem_dealloc(anypoints)
    call mem_dealloc(DistAng)


  end subroutine plot_pair_energies


  ! ===================================================================!
  !                      FRAGMENT OPTIMIZATION                         !
  ! ===================================================================!
  
  !> Purpose: Routine that optimizes an atomic fragment by (1) expanding to include 
  !           a given number of orbital based on a priority list. Checking that energy 
  !           change is smaller than FOT. (2) For converged expanded fragment we 
  !           reduce fragment again by removing individual orbitals based on (more accurate)
  !           priority list. 
  !           In (2) it is ensured that the energy error introduced by reducing the 
  !           fragment (compared to the fragment obtained in (1)) is below the FOT.
  !
  !> Author:  Pablo Baudin, (based on previous work by Ida-Marie Hoeyvik,
  !           Kasper Kristensen & Thomas Kjaergaard)
  !> Date:    July 2014
  subroutine optimize_atomic_fragment(MyAtom,AtomicFragment,natoms, &
           & OccOrbitals,no,VirOrbitals,nv,MyMolecule,mylsitem,freebasisinfo)

     implicit none

     !> Number of occupied orbitals in molecule
     integer, intent(in) :: no
     !> Number of virtual orbitals in molecule
     integer, intent(in) :: nv
     !> Number of atoms in molecule
     integer, intent(in) :: natoms
     !> Central atom in molecule
     integer, intent(in) :: MyAtom
     !> Atomic fragment to be optimized
     type(decfrag), intent(inout) :: AtomicFragment
     !> All occupied orbitals
     type(decorbital), dimension(no), intent(in) :: OccOrbitals
     !> All virtual orbitals
     type(decorbital), dimension(nv), intent(in) :: VirOrbitals
     !> Full molecule information
     type(fullmolecule), intent(inout) :: MyMolecule
     !> Integral information
     type(lsitem), intent(inout) :: mylsitem
     !> Delete fragment basis information ("expensive box in decfrag type") at exit?
     logical, intent(in) :: freebasisinfo

     !> number of occ/vir orbitals in each atom:
     integer, dimension(natoms) :: nocc_per_atom, nvir_per_atom
     !> Logical vector telling which orbital is include in the fragment
     logical, pointer :: Occ_AOS(:), Vir_AOS(:)
     !> Expansion priority lists:
     integer, pointer :: exp_list_occ(:), exp_list_vir(:)
     !> # occ/vir orbitals added in each step of the expansion:
     integer :: nexp_occ, nexp_vir
     !> # occ/vir orbitals added to EOS space to define the initial fragment:
     integer :: ninit_occ, ninit_vir
     !> Reduction priority list
     integer, pointer :: red_list_occ(:)
     integer, pointer :: red_list_vir(:)
     !> minimum gap in number of orbital allowed between the 
     !  last two steps of the binary search.
     integer :: nred_occ, nred_vir
     !> true if/when the fragment include the full molecule:
     logical :: full_mol, dored
     !> Timings for frag opt
     real(realk), pointer :: times_fragopt(:)  
     !> number of core orbitals
     integer :: nc
     integer :: idx,i,a,expit,redit
     real(realk), dimension(ndecenergies) :: Eexp
     logical :: add_MP2_opt


     !==================================================================================!
     !                       Initialization of various things...                        !
     !==================================================================================!
     !
     ! Start timings:
     times_fragopt => null()
     call dec_fragment_time_init(times_fragopt)
     add_MP2_opt = .false.
     Eexp = 0E0_realk

     expit=0
     redit=0
     if( DECinfo%print_small_calc )then
        write(DECinfo%output,'(a)')    ' FOP'
        write(DECinfo%output,'(a)')    ' FOP ==============================================='
        write(DECinfo%output,'(a,i4)') ' FOP  Site fragment generator for fragment,',MyAtom
        write(DECinfo%output,'(a)')    ' FOP ==============================================='
        write(DECinfo%output,'(a)')    ' FOP'
     endif


     ! Get number of orbital in EOS spaces:
     nocc_per_atom=get_number_of_orbitals_per_atom(OccOrbitals,no,natoms,.true.)
     nvir_per_atom=get_number_of_orbitals_per_atom(VirOrbitals,nv,natoms,.true.)

     ! Only do fragment optimization if there are orbitals assigned to central atom.
     if( (nocc_per_atom(MyAtom) == 0) .and. (nvir_per_atom(MyAtom) == 0) ) then
        write(DECinfo%output,*) 'FOP Skipping optimization of fragment ', MyAtom
        AtomicFragment%energies = 0E0_realk
        AtomicFragment%ccmodel  = DECinfo%ccmodel

        ! Get final timings before quiting
        call dec_fragment_time_get(times_fragopt)
        call dec_time_evaluate_efficiency_frag(AtomicFragment,times_fragopt,&
           & AtomicFragment%ccmodel,'Fragment optmization')
        return
     end if

     ! Debug case: Full molecule included in fragment
     ! --> we then skip energy calculation here and just init fragment
     if(DECinfo%InclFullMolecule .or. DECinfo%simulate_full) then
        call fragopt_include_fullmolecule(MyAtom,AtomicFragment, &
           & OccOrbitals,no,VirOrbitals,nv,MyMolecule,mylsitem,freebasisinfo)

        ! Get final timings before quiting
        call dec_fragment_time_get(times_fragopt)
        call dec_time_evaluate_efficiency_frag(AtomicFragment,times_fragopt,&
           & AtomicFragment%ccmodel,'Fragment optmization')
        return
     end if

     ! Overwrite ccmodel for myatom with the expansion required ccmodel
     MyMolecule%ccmodel(MyAtom,Myatom) = DECinfo%fragopt_exp_model
     select case(DECinfo%fragopt_exp_model)
     case(MODEL_CC2,MODEL_CCSD,MODEL_CCSDpT)
        add_MP2_opt = .true.
     end select


     ! Get information on how the expansion should be performed:
     call mem_alloc(exp_list_occ,no)
     call mem_alloc(exp_list_vir,nv)


     !Define the EOS spaces in the fragment to be able to calculate priority lists
     !============================================================================

     if(DECinfo%frozencore)then
        nc = MyMolecule%ncore
     else
        nc = 0
     endif

     ! Occ orbitals not in the core
     idx = 0
     do i=nc+1,MyMolecule%nocc
        if(OccOrbitals(i)%centralatom == MyAtom)then
           idx = idx + 1
        endif
     enddo

     AtomicFragment%noccEOS   = idx
     AtomicFragment%nvirtEOS = nvir_per_atom(MyAtom)
     call mem_alloc(AtomicFragment%occEOSidx,AtomicFragment%noccEOS)
     call mem_alloc(AtomicFragment%virtEOSidx,AtomicFragment%nvirtEOS)

     idx = 1
     do i=nc+1,MyMolecule%nocc
        if(OccOrbitals(i)%centralatom == MyAtom)then
           AtomicFragment%occEOSidx(idx) = i
           idx = idx + 1
        endif
     enddo
     if(idx-1/=AtomicFragment%noccEOS)then
        call lsquit("ERROR finding the right amount of occupied EOS orbitals",-1)
     endif

     idx = 1
     do a=1,MyMolecule%nvirt
        if(VirOrbitals(a)%centralatom == MyAtom)then
           AtomicFragment%virtEOSidx(idx) = a
           idx = idx + 1
        endif
     enddo
     if(idx-1/=AtomicFragment%nvirtEOS)then
        call lsquit("ERROR finding the right amount of occupied EOS orbitals",-1)
     endif


     call define_frag_expansion(no,nv,natoms,MyAtom,MyMolecule,AtomicFragment, &
        & exp_list_occ,exp_list_vir,nexp_occ,nexp_vir,ninit_occ,ninit_vir)

     call mem_dealloc(AtomicFragment%occEOSidx)
     call mem_dealloc(AtomicFragment%virtEOSidx)

     !==================================================================================!
     !                               Initialize Fragment                                !
     !==================================================================================!
     !
     ! Get Initial fragment for Myatom, we include the EOS space plus Frag_Init_Size 
     ! occ/vir orbitals from the priority lists:
     call mem_alloc(Occ_AOS,no)
     call mem_alloc(Vir_AOS,nv)

     ! Get logical list Occ_AOS/Vir_AOS to know which orbitals to include:
     call expand_fragment(no,nv,exp_list_occ,exp_list_vir,ninit_occ,ninit_vir,MyAtom, &
        & MyMolecule,OccOrbitals,VirOrbitals,Occ_AOS,Vir_AOS,full_mol,frag_init=.true.)
     ! Initialize fragment base on the two orbital lists Occ_AOS/Vir_AOS:
     call atomic_fragment_init_orbital_specific(MyAtom,nv,no,Vir_AOS,Occ_AOS,OccOrbitals, &
        & VirOrbitals,MyMolecule,mylsitem,AtomicFragment,.true.,.false.) 
     ! Get Energy for the initialized fragment
     call fragment_energy_and_prop(AtomicFragment,add_mp2_opt=add_mp2_opt) 
     ! Print initial fragment information
     if (full_mol) write(DECinfo%output,*) 'FOP Expansion Include Full Molecule !!!'
     call fragopt_print_info(AtomicFragment,AtomicFragment%ccmodel,0)




     !==================================================================================!
     !                             Enter Fragment Expansion                             !
     !==================================================================================!
     !
     ! Do expansion only if the initial fragment does not include the full molecule:
     if (.not.full_mol) then
        call fragment_expansion_procedure(Occ_AOS,Vir_AOS,AtomicFragment,no,nv, &
           & exp_list_occ,exp_list_vir,nexp_occ,nexp_vir,MyAtom,MyMolecule, &
           & OccOrbitals,VirOrbitals,expit,add_MP2_opt,mylsitem)
     end if



     !==================================================================================!
     !                  Transition from expansion to reduction loop                     !
     !==================================================================================!

     if( DECinfo%print_small_calc )then
        write(DECinfo%output,*) 'FOP'
        write(DECinfo%output,*) 'FOP ============================================='
        write(DECinfo%output,*) 'FOP  Expansion has converged. We start reduction '
        write(DECinfo%output,*) 'FOP ============================================='
        write(DECinfo%output,*) 'FOP'
     endif

     ! Deallocate exp list and allocate red ones:
     call mem_dealloc(exp_list_occ)
     call mem_dealloc(exp_list_vir)
     call mem_alloc(red_list_occ,no)
     call mem_alloc(red_list_vir,nv)

     ! Overwrite ccmodel for myatom with the reduction required ccmodel
     MyMolecule%ccmodel(MyAtom,Myatom) = DECinfo%fragopt_red_model
     AtomicFragment%ccmodel = DECinfo%fragopt_red_model
     select case(DECinfo%fragopt_red_model)
     case(MODEL_CC2,MODEL_CCSD,MODEL_CCSDpT)
        add_MP2_opt = .true.
     end select

     ! if expansion and reduction ccmodels are different, we need to calculate
     ! the energy of the expanded fragment using the new model:
     if(DECinfo%fragopt_exp_model /= DECinfo%fragopt_red_model) then
        call fragment_energy_and_prop(AtomicFragment,add_mp2_opt=add_mp2_opt)
        if( DECinfo%print_small_calc )then
           write(DECinfo%output,'(2a)') ' FOP Calculated ref atomic fragment energy for relevant CC model: ', &
              & DECinfo%cc_models(MyMolecule%ccmodel(MyAtom,Myatom))
           call fragopt_print_info(AtomicFragment,AtomicFragment%ccmodel,0)
        endif
     end if

     !Store the energies of the expanded fragment
     Eexp = AtomicFragment%energies

     ! Get information on how the reduction should be performed:
     call define_frag_reduction(no,nv,natoms,MyAtom,MyMolecule,AtomicFragment, &
        & red_list_occ,red_list_vir,nred_occ,nred_vir)


     ! For DECNP no pairs are calculated so we do reduction only if:
     !    CC-model /= Exp-model .and. CC-model /= Red-model
     dored = .true.
     if (DECinfo%DECNP) dored = ( (DECinfo%ccmodel/=DECinfo%fragopt_exp_model) .and. &
        & (DECinfo%ccmodel/=DECinfo%fragopt_red_model) )


     ! Perform reduction:
     if (dored) then
        call fragment_reduction_procedure_wrapper(AtomicFragment,no,nv,red_list_occ, &
              & red_list_vir,Occ_AOS,Vir_AOS,MyAtom,MyMolecule,OccOrbitals, &
              & VirOrbitals,mylsitem,nred_occ,nred_vir,DECinfo%FOT,redit,add_MP2_opt)
     end if

     !==================================================================================!
     !                              Finalize subroutine                                 !
     !==================================================================================!

     !Store the error estimates for the fragments
     call get_fragopt_energy_error(AtomicFragment,AtomicFragment%ccmodel,Eexp)

     ! Deallocation:
     call mem_dealloc(red_list_occ)
     call mem_dealloc(red_list_vir)
     call mem_dealloc(Occ_AOS)
     call mem_dealloc(Vir_AOS)

     if(freebasisinfo) then
        call atomic_fragment_free_basis_info(AtomicFragment)
     end if

     ! Get final Timings:
     call dec_fragment_time_get(times_fragopt)
     call dec_time_evaluate_efficiency_frag(AtomicFragment,times_fragopt, &
        & AtomicFragment%ccmodel,'Fragment optmization')

     ! Restore the original CC model 
     ! (only relevant if expansion and/or reduction was done using the MP2 model, but it doesn't hurt)
     MyMolecule%ccmodel(MyAtom,Myatom) = DECinfo%ccmodel

     ! Fragment has been optimized
     AtomicFragment%isopt = .true.

     ! Final print:
     write(DECinfo%output,'(/A,I4,A,I4,A,I4,A/)') ' FOP SUMMARY: Fragment,',MyAtom, &
        & '  converged after ',expit,' expansion and ',redit,' reduction steps'

  end subroutine optimize_atomic_fragment


  !> Purpose: Perform expansion procedure based on occ/vir priority list. In each step
  !           the fragment AOS is increased using the nexp most important orbitals from
  !           occupied and virtual priority lists.
  !           The expansion is converged when the difference in energie(s) between two
  !           step is smaller than the FOT.
  !
  !> Author:  Pablo Baudin (based on previous work by Kasper Kristensen & Thomas Kjaergaard)
  !> Date:    July 2014
  subroutine fragment_expansion_procedure(Occ_AOS,Vir_AOS,AtomicFragment,no,nv, &
           & occ_priority_list,vir_priority_list,nexp_occ,nexp_vir,MyAtom,MyMolecule, &
           & OccOrbitals,VirOrbitals,iter,add_MP2_opt,mylsitem)

     implicit none

     !> Number of occupied orbitals in molecule
     integer, intent(in) :: no
     !> Number of virtual orbitals in molecule
     integer, intent(in) :: nv
     !> Atomic fragment to be optimized
     type(decfrag),intent(inout) :: AtomicFragment
     !> Priority list of orbitals:
     integer, intent(in) :: occ_priority_list(no)
     integer, intent(in) :: vir_priority_list(nv)
     !> Number of orbital to add in each expansion step
     integer, intent(in) :: nexp_occ, nexp_vir
     !> Central Atom of the current fragment
     integer, intent(in) :: MyAtom
     !> Full molecule information
     type(fullmolecule), intent(inout) :: MyMolecule
     !> All occupied orbitals
     type(decorbital), dimension(no), intent(in) :: OccOrbitals
     !> All virtupied orbitals
     type(decorbital), dimension(nv), intent(in) :: VirOrbitals
     !> Logical vector telling which orbital is include in the fragment
     logical, intent(inout) :: Occ_AOS(no), Vir_AOS(nv)
     !> Out: number of iteration used in the expansion
     integer, intent(out) :: iter
     !> Get MP2 information on top of actual ccmodel:
     logical, intent(in) :: add_MP2_opt
     !> Integral information
     type(lsitem), intent(inout)       :: mylsitem
     
     logical :: full_mol, expansion_converged, MP2_converged
     real(realk), dimension(ndecenergies) :: Eold



     EXPANSION_LOOP: do iter=1,DECinfo%maxiter

        ! Save current fragment energies
        Eold = AtomicFragment%energies

        ! Expand fragment:
        call expand_fragment(no,nv,occ_priority_list,vir_priority_list,nexp_occ,nexp_vir, &
           & MyAtom,MyMolecule,OccOrbitals,VirOrbitals,Occ_AOS,Vir_AOS,full_mol)

        ! Free old fragment and initialize the expanded fragment:
        call atomic_fragment_free(AtomicFragment)
        call atomic_fragment_init_orbital_specific(MyAtom,nv,no,Vir_AOS,Occ_AOS,OccOrbitals, &
           & VirOrbitals,MyMolecule,mylsitem,AtomicFragment,.true.,.false.)

        ! Get new fragment energy:
        call fragment_energy_and_prop(AtomicFragment,add_mp2_opt=add_mp2_opt)

        ! Energy differences
        call get_fragopt_energy_error(AtomicFragment,AtomicFragment%ccmodel,Eold)

        ! print expanded fragment information:
        if (full_mol) write(DECinfo%output,*) 'FOP Expansion Include Full Molecule !!!'
        call fragopt_print_info(AtomicFragment,AtomicFragment%ccmodel,iter)

        ! Test if fragment energy (or energies) are converged to FOT precision
        call fragopt_check_convergence(AtomicFragment%ccmodel,AtomicFragment, &
           & DECinfo%FOT,expansion_converged)


        ! Test if fragment MP2 energies are converged
        if (add_mp2_opt) then

           call get_fragopt_energy_error(AtomicFragment,MODEL_MP2,Eold)
           call fragopt_check_convergence(MODEL_MP2,AtomicFragment,DECinfo%FOT,MP2_converged)
           call fragopt_print_info(AtomicFragment,MODEL_MP2,iter,add_mp2_opt)

           ! require both MP2 and CC energy to be converged
           if (.not.MP2_converged) expansion_converged = .false.
        end if

        ! Set the expansion to be converged if the current fragment include the full molecule:
        if (full_mol) expansion_converged = .true.


        ! Exit loop if we are converged
        ExpansionConvergence: if(expansion_converged) then
           if( DECinfo%print_small_calc )then
              write(DECinfo%output,*) 'FOP Fragment expansion converged in iteration ', iter
           endif
           exit EXPANSION_LOOP
        end if ExpansionConvergence

     end do EXPANSION_LOOP

     ! Check that expansion loop is converged
     if(.not. expansion_converged) then
        write(DECinfo%output,*) 'Number of expansion steps = ', DECinfo%MaxIter
        call lsquit('Fragment expansion did not converge! Try to increase the number of &
           & expansion steps using the .MaxIter keyword',DECinfo%output)
     end if

  end subroutine fragment_expansion_procedure


  !> \brief   Wrapper for reducing fragment by removing orbitals until error
  !>          is of size FOT. This is also done for reduced fragment, i.e. larger FOTs.
  !> \author Kasper Kristensen
  !> \date November 2014
  subroutine fragment_reduction_procedure_wrapper(AtomicFragment,no,nv,occ_priority_list, &
       & vir_priority_list,Occ_AOS,Vir_AOS,MyAtom,MyMolecule,OccOrbitals, &
       & VirOrbitals,mylsitem,no_gap,nv_gap,FOT,iter,add_MP2_opt)

    implicit none

    !> Atomic fragment to be optimized
    type(decfrag),intent(inout)        :: AtomicFragment
    !> Number of occupied orbitals in molecule
    integer, intent(in) :: no
    !> Number of virtual orbitals in molecule
    integer, intent(in) :: nv
    !> Priority list of orbitals:
    integer, intent(in) :: occ_priority_list(no)
    integer, intent(in) :: vir_priority_list(nv)
    !> Logical vector telling which orbital is include in the fragment
    logical, intent(inout) :: Occ_AOS(no), Vir_AOS(nv)
    !> Central Atom of the current fragment
    integer, intent(in) :: MyAtom
    !> Full molecule information
    type(fullmolecule), intent(in) :: MyMolecule
    !> All occupied orbitals
    type(decorbital), dimension(no), intent(in) :: OccOrbitals
    !> All virtupied orbitals
    type(decorbital), dimension(nv), intent(in) :: VirOrbitals
    !> Integral information
    type(lsitem), intent(inout)       :: mylsitem
    !> minimum gap in number of orbital allowed between the 
    !  last two steps of the binary search.
    integer,intent(in) :: no_gap, nv_gap
    !> Fragment optimization threshold to use in reduction
    real(realk),intent(in) :: FOT
    !> Out: number of iteration used in the reduction
    integer, intent(out) :: iter
    !> Get MP2 info on top of actual cc model:
    logical, intent(in) :: add_MP2_opt

    !> MP2 amplitudes and integrals for fragment adapted reduction:
    type(tensor) :: g, t2
    real(realk), dimension(ndecenergies) :: Eexp
    real(realk) :: FOTincreased
    type(decfrag) :: ReducedFragment
    integer :: i,noccAOSprev, nvirtAOSprev
    integer,pointer :: occAOSidxprev(:),virtAOSidxprev(:)

    ! Sanity check
    if(DECinfo%FOTscaling<1.0_realk) then
       print *, 'Scaling factor ',DECinfo%FOTscaling
       call lsquit('fragment_reduction_procedure_wrapper: Requires scaling factor &
          & larger than one! ',-1)
    end if


    ! Reduce fragment space using fragment adapted orbital:
    ! =====================================================
    FragAdapt: if(DECinfo%fragadapt) then

       ! Get MP2 amplitudes and energy for fragment
       call mp2_solver(AtomicFragment,g,t2,.false.)

       ! Get correlation density matrix for atomic fragment
       call calculate_MP2corrdens_frag(t2,AtomicFragment)
       call tensor_free(t2)
       call tensor_free(g)

       ! Reduce using fragment-adapted orbitals
       call fragopt_reduce_FOs(MyAtom,AtomicFragment,OccOrbitals,no,VirOrbitals,nv, &
          & MyMolecule,mylsitem,t2,iter)

       return
    end if FragAdapt
    

    ! At this point AtomicFragment correspondings to the expanded fragment
    ! We store the atomic fragment energies of the expanded fragment
    Eexp = AtomicFragment%energies

    iter=0
    ! Determine AtomicFragment according to main FOT
    call fragment_reduction_procedure(AtomicFragment,no,nv,occ_priority_list, &
         & vir_priority_list,Occ_AOS,Vir_AOS,MyAtom,MyMolecule,OccOrbitals, &
         & VirOrbitals,mylsitem,no_gap,nv_gap,FOT,iter,add_MP2_opt)

    if(DECinfo%print_small_calc)then
       write(DECinfo%output,'(1X,a,i7,g14.3,3i7)') 'FOP reduction: Atom,FOT,O,V,B',MyAtom,FOT,&
          & AtomicFragment%noccAOS,AtomicFragment%nvirtAOS,AtomicFragment%nbasis
    endif

    ! Store AOS space for converged AtomicFragment 
    noccAOSprev=AtomicFragment%noccAOS
    nvirtAOSprev=AtomicFragment%nvirtAOS
    call mem_alloc(occAOSidxprev,noccAOSprev)
    call mem_alloc(virtAOSidxprev,nvirtAOSprev)
    occAOSidxprev = AtomicFragment%occAOSidx
    virtAOSidxprev = AtomicFragment%virtAOSidx

    
    ! Loop over different increased FOTs to determine reduced spaces
    ! *************************************************************
    FOTincreased=FOT
    do i=1,DECinfo%nFRAGSred

       ! Initialize ReducedFragment identical fragment for previous FOT
       ! (converged fragment for i=1).
       call atomic_fragment_init_integer_list(MyAtom,nv, no, nvirtAOSprev,&
            & noccAOSprev,virtAOSidxprev,occAOSidxprev,&
            & OccOrbitals,VirOrbitals,MyMolecule,mylsitem,ReducedFragment,.true.,.false.)

       ! Set initial energies for ReducedFragment equal to the ones for the expanded fragment
       AtomicFragment%energies = Eexp

       ! Increase FOT by scaling factor
       FOTincreased = FOTincreased*DECinfo%FOTscaling

       ! Determine ReducedFragment according to the scaled FOT
       call fragment_reduction_procedure(ReducedFragment,no,nv,occ_priority_list, &
            & vir_priority_list,Occ_AOS,Vir_AOS,MyAtom,MyMolecule,OccOrbitals, &
            & VirOrbitals,mylsitem,no_gap,nv_gap,FOTincreased,iter,add_MP2_opt)
       
       ! Save information about ReducedFragment AOS in AtomicFragment structure
       AtomicFragment%REDfrags(i)%noccAOS   = ReducedFragment%noccAOS
       AtomicFragment%REDfrags(i)%nvirtAOS = ReducedFragment%nvirtAOS
       call mem_alloc(AtomicFragment%REDfrags(i)%occAOSidx,AtomicFragment%REDfrags(i)%noccAOS)
       AtomicFragment%REDfrags(i)%occAOSidx = ReducedFragment%occAOSidx
       call mem_alloc(AtomicFragment%REDfrags(i)%virtAOSidx,AtomicFragment%REDfrags(i)%nvirtAOS)
       AtomicFragment%REDfrags(i)%virtAOSidx = ReducedFragment%virtAOSidx
       AtomicFragment%REDfrags(i)%FOT = FOTincreased

       ! Print summary (delete this at some point but nice to have for analysis now)
       if(DECinfo%print_small_calc)then
          write(DECinfo%output,'(1X,a,i7,g14.3,3i7)') 'FOP reduction: Atom,FOT,O,V,B',MyAtom,&
             & FOTincreased, ReducedFragment%noccAOS,ReducedFragment%nvirtAOS,&
             & ReducedFragment%nbasis
       endif

       ! Store AOS information to use as starting point in fragment for next FOT
       call mem_dealloc(occAOSidxprev)
       call mem_dealloc(virtAOSidxprev)
       noccAOSprev=ReducedFragment%noccAOS
       nvirtAOSprev=ReducedFragment%nvirtAOS
       call mem_alloc(occAOSidxprev,noccAOSprev)
       call mem_alloc(virtAOSidxprev,nvirtAOSprev)
       occAOSidxprev = ReducedFragment%occAOSidx
       virtAOSidxprev = ReducedFragment%virtAOSidx

       ! Done with reduced fragment
       call atomic_fragment_free(ReducedFragment)

    end do

    call mem_dealloc(occAOSidxprev)
    call mem_dealloc(virtAOSidxprev)


  end subroutine fragment_reduction_procedure_wrapper


  !> Purpose: Perform reduction procedure based on occ/vir priority list. We perform
  !           a binary search on the priority list and accept a step based on energy
  !           criterions (dE_occ, dE_vir). The binary search stops for one space when
  !           the difference in number of orbitals between two steps is lower than
  !           no_gap/nv_gap. By default the occupied spaced is reduced 1st.
  !
  !> Author:  Pablo Baudin
  !> Date:    July 2014
  subroutine fragment_reduction_procedure(AtomicFragment,no,nv,occ_priority_list, &
           & vir_priority_list,Occ_AOS,Vir_AOS,MyAtom,MyMolecule,OccOrbitals, &
           & VirOrbitals,mylsitem,no_gap,nv_gap,FOT,totit,add_MP2_opt)

     implicit none

     !> Atomic fragment to be optimized
     type(decfrag),intent(inout)        :: AtomicFragment
     !> Number of occupied orbitals in molecule
     integer, intent(in) :: no
     !> Number of virtual orbitals in molecule
     integer, intent(in) :: nv
     !> Priority list of orbitals:
     integer, intent(in) :: occ_priority_list(no)
     integer, intent(in) :: vir_priority_list(nv)
     !> Logical vector telling which orbital is include in the fragment
     logical, intent(inout) :: Occ_AOS(no), Vir_AOS(nv)
     !> Central Atom of the current fragment
     integer, intent(in) :: MyAtom
     !> Full molecule information
     type(fullmolecule), intent(in) :: MyMolecule
     !> All occupied orbitals
     type(decorbital), dimension(no), intent(in) :: OccOrbitals
     !> All virtupied orbitals
     type(decorbital), dimension(nv), intent(in) :: VirOrbitals
     !> Integral information
     type(lsitem), intent(inout)       :: mylsitem
     !> minimum gap in number of orbital allowed between the 
     !  last two steps of the binary search.
     integer,intent(in) :: no_gap, nv_gap
     !> Fragment optimization threshold to use in reduction
     real(realk),intent(in) :: FOT
     !> Total number of iteration used in the reduction
     integer, intent(inout) :: totit
     !> Get MP2 info on top of actual cc model:
     logical, intent(in) :: add_MP2_opt

     !> energy error acceptance in reduction steps:
     real(realk) :: dE_occ, dE_vir
     integer :: no_exp, nv_exp, no_min, nv_min, no_max, nv_max
     integer :: no_old, nv_old, no_new, nv_new, iter
     logical :: redocc, redvir, step_accepted, MP2_accepted
     logical :: reduction_converged, occ_red_conv, vir_red_conv
     logical, pointer :: OccAOS_old(:), VirAOS_old(:)
     real(realk), dimension(ndecenergies) :: Eexp, Eold


     ! INITIALIZATION:
     ! ***************
     call mem_alloc(OccAOS_old,no)
     call mem_alloc(VirAOS_old,nv)

     no_exp = AtomicFragment%noccAOS
     nv_exp = AtomicFragment%nvirtAOS
     Eexp = AtomicFragment%energies

     ! Set boundaries of the list (min = EOS, max = Expanded fragment AOS)
     ! note: if frozen core then core orbitals are already excluded from those numbers
     no_min = AtomicFragment%noccEOS
     nv_min = AtomicFragment%nvirtEOS
     no_max = AtomicFragment%noccAOS
     nv_max = AtomicFragment%nvirtAOS

     ! The initial (or old) condition correspond to the expanded fragment
     no_old = no_exp
     nv_old = nv_exp
     OccAOS_old = Occ_AOS
     VirAOS_old = Vir_AOS
     Eold = Eexp

     ! Define specific reduction parameters:
     ! -------------------------------------
     if (DECinfo%Frag_red_occ.and.(.not.DECinfo%Frag_red_virt)) then
        ! Start reducing occupied space:
        if(DECinfo%print_small_calc)then
           write(DECinfo%output,'(1X,a,/)') 'FOP: User chose to reduce occupied space first'
        endif
        redocc = .true.
        redvir = .false.
        dE_occ = DECinfo%frag_red1_thr*FOT
        dE_vir = DECinfo%frag_red2_thr*FOT
     else if (DECinfo%Frag_red_virt.and.(.not.DECinfo%Frag_red_occ)) then
        ! Start reducing virtual space:
        if(DECinfo%print_small_calc)then
           write(DECinfo%output,'(1X,a,/)') 'FOP: User chose to reduce virtual space first'
        endif
        redvir = .true.
        redocc = .false.
        dE_vir = DECinfo%frag_red1_thr*FOT
        dE_occ = DECinfo%frag_red2_thr*FOT
     else
        ! Default, reduce both spaces from the begining:
        redocc = .true.
        redvir = .true.
        dE_occ = FOT
        dE_vir = FOT
     end if


     step_accepted = .false.
     occ_red_conv  = .false.
     vir_red_conv  = .false.
     reduction_converged = .false.

     ! Check That fragment can be reduced:
     if (nv_min==nv_max) then ! virtual space cannot be reduced:
        vir_red_conv = .true.
        redvir       = .false.
     end if
     if (no_min==no_max) then ! occupied space cannot be reduced
        occ_red_conv = .true.
        redocc       = .false.
     end if
     if (vir_red_conv.and.occ_red_conv) reduction_converged = .true.


     REDUCTION_LOOP: do iter=1,DECinfo%maxiter

        ! Exit before doing anything if the fragment cannot be reduced:
        if (reduction_converged) exit REDUCTION_LOOP


        ! REDUCE FRAGMENT:
        ! ****************
        if (redocc.and.(.not.redvir)) then
           ! keep old virtual space and change occupied one
           nv_new = nv_old
           ! The new condition is set by removing half of the occ 
           ! orbital between max and min:
           no_new = no_max - (no_max - no_min)/2
        else if (redvir.and.(.not.redocc)) then
           ! keep old occupied space and change virtual one
           no_new = no_old
           ! The new condition is set by removing half of the vir
           ! orbital between max and min:
           nv_new = nv_max - (nv_max - nv_min)/2
        else if (redocc.and.redvir) then
           ! The new condition is set by removing half of the occ
           ! and vir orbital between max and min:
           no_new = no_max - (no_max - no_min)/2
           nv_new = nv_max - (nv_max - nv_min)/2
        else
           call lsquit("ERROR(fragment_reduction_procedure): space to reduce not defined", &
              & DECinfo%output)
        end if

        call reduce_fragment(AtomicFragment,MyMolecule,no,no_new,Occ_AOS, &
           & occ_priority_list,.true.)
        call reduce_fragment(AtomicFragment,MyMolecule,nv,nv_new,Vir_AOS, &
           & vir_priority_list,.false.)

        ! Free old fragment and initialize the reduced fragment:
        call atomic_fragment_free(AtomicFragment)
        call atomic_fragment_init_orbital_specific(MyAtom,nv,no,Vir_AOS,Occ_AOS,OccOrbitals, &
           & VirOrbitals,MyMolecule,mylsitem,AtomicFragment,.true.,.false.)

        ! Get new fragment energy:
        call fragment_energy_and_prop(AtomicFragment,add_mp2_opt=add_mp2_opt)

        ! Energy differences
        call get_fragopt_energy_error(AtomicFragment,AtomicFragment%ccmodel,Eexp)

        ! print reduced fragment information:
        call fragopt_print_info(AtomicFragment,AtomicFragment%ccmodel,iter)


        ! CHECK CONVERGENCES:
        ! *******************
        ! Check if reduction step is accepted (Energy criterion):
        if (redocc.and.(.not.redvir)) then
           call fragopt_check_convergence(AtomicFragment%ccmodel,AtomicFragment,dE_occ,step_accepted)
        else if (redvir.and.(.not.redocc)) then
           call fragopt_check_convergence(AtomicFragment%ccmodel,AtomicFragment,dE_vir,step_accepted)
        else if (redocc.and.redvir) then
           call fragopt_check_convergence(AtomicFragment%ccmodel,AtomicFragment,FOT,step_accepted)
        end if


        ! Test if fragment MP2 energies are converged
        if (add_mp2_opt) then

           call get_fragopt_energy_error(AtomicFragment,MODEL_MP2,Eexp)

           if (redocc.and.(.not.redvir)) then
              call fragopt_check_convergence(MODEL_MP2,AtomicFragment,dE_occ,MP2_accepted)
           else if (redvir.and.(.not.redocc)) then
              call fragopt_check_convergence(MODEL_MP2,AtomicFragment,dE_vir,MP2_accepted)
           else if (redocc.and.redvir) then
              call fragopt_check_convergence(MODEL_MP2,AtomicFragment,FOT,MP2_accepted)
           end if

           call fragopt_print_info(AtomicFragment,MODEL_MP2,iter,add_mp2_opt)

           ! require both MP2 and CC energy to be converged
           if (.not.MP2_accepted) then
              step_accepted = .false.
           end if
        end if


        ! If the step is accepted we need to save the frag info:
        if (step_accepted) then
           Eold = AtomicFragment%energies
           OccAOS_old = Occ_AOS
           VirAOS_old = Vir_AOS
        end if


        ! Check convergence of spaces:
        SpaceConvergence: if (redocc.and.(.not.redvir)) then
           call check_red_space_convergence(no_old,no_new,no_min,no_max,no_gap,occ_red_conv,vir_red_conv, &
              & reduction_converged,redocc,step_accepted)
        else if (redvir.and.(.not.redocc)) then
           call check_red_space_convergence(nv_old,nv_new,nv_min,nv_max,nv_gap,vir_red_conv,occ_red_conv, &
              & reduction_converged,redocc,step_accepted)
        else if (redvir.and.redocc) then
           call check_red_space_convergence(no_old,no_new,no_min,no_max,no_gap,occ_red_conv,vir_red_conv, &
              & reduction_converged,redocc,step_accepted)
           call check_red_space_convergence(nv_old,nv_new,nv_min,nv_max,nv_gap,vir_red_conv,occ_red_conv, &
              & reduction_converged,.not.redocc,step_accepted)
        end if SpaceConvergence

        ! If everything has converged then we keep the last valid
        ! information and quit the loop:
        FullConvergence: if (reduction_converged) then
           if (DECinfo%PL > 1 .and. DECinfo%print_small_calc) write(DECinfo%output,*) 'BIN SEARCH: REDUCTION CONVERGED'

           if (.not.step_accepted) then
              ! The last binary search step did not succeed so  
              ! we init everything to the last successful step
              call atomic_fragment_free(AtomicFragment)
              call atomic_fragment_init_orbital_specific(MyAtom,nv,no,VirAOS_old,OccAOS_old, &
                   & OccOrbitals,VirOrbitals,MyMolecule,mylsitem,AtomicFragment,.true.,.false.)
              AtomicFragment%energies = Eold
              Occ_AOS = OccAOS_old
              Vir_AOS = VirAOS_old
           end if
           exit REDUCTION_LOOP
        end if FullConvergence


        ! SETTING PARAMETERS FOR NEXT ITERATION:
        ! **************************************
        no_old = no_new
        nv_old = nv_new
        if (step_accepted) then
           if (DECinfo%PL > 1.and. DECinfo%print_small_calc) write(DECinfo%output,*) 'BIN SEARCH: Step accepted'
           ! The current number of occ and virt is good or too high:
           ! Setting new maximum:
           if (redocc.and.(.not.redvir)) then
              no_max = no_new
           else if (redvir.and.(.not.redocc)) then
              nv_max = nv_new
           else if (redvir.and.redocc) then
              no_max = no_new
              nv_max = nv_new
           end if
        else
           if (DECinfo%PL > 1.and. DECinfo%print_small_calc) write(DECinfo%output,*) 'BIN SEARCH: Step NOT accepted'
           ! The current number of occ and virt is NOT good enough:
           ! Setting new minimum:
           if (redocc.and.(.not.redvir)) then
              no_min = no_new
           else if (redvir.and.(.not.redocc)) then
              nv_min = nv_new
           else if (redvir.and.redocc) then
              no_min = no_new
              nv_min = nv_new
           end if
        end if

        ! swicht spaces to be reduce only when the initial one is fully converged
        if (vir_red_conv) then
           redocc = .true.
           redvir = .false.
        else if (occ_red_conv) then
           redocc = .false.
           redvir = .true.
        end if

        ! Sanity checks:
        if (no_min > no_max) call lsquit('FOP ERROR: no_min > no_max', DECinfo%output)
        if (nv_min > nv_max) call lsquit('FOP ERROR: nv_min > nv_max', DECinfo%output)

     end do REDUCTION_LOOP
     totit = totit + iter

     ! CHECK THAT REDUCTION CONVERGED:
     ! *******************************
     if (.not.reduction_converged) then
       write(DECinfo%output,*) "FATAL ERROR: Reduction not converged for fragment",MyAtom
       write(DECinfo%output,*) "          Exit loop after",iter,"iterations"
       call lsquit('Fragment reduction not converged',DECinfo%output)
     end if

     call mem_dealloc(OccAOS_old)
     call mem_dealloc(VirAOS_old)

     ! PRINT INFO FOR FINAL (REDUCED) FRAGMENT:
     ! ****************************************
     if(DECinfo%print_small_calc)then
        write(DECinfo%output,'(1X,a,/)') 'FOP'
        write(DECinfo%output,'(1X,a)') 'FOP========================================================='
        write(DECinfo%output,'(1X,a,i4)') 'FOP    LOCAL REDUCTION HAS CONVERGED FOR SITE',MyAtom
        write(DECinfo%output,'(1X,a)') 'FOP---------------------------------------------------------'
        write(DECinfo%output,'(1X,a,i4)')    'FOP Done: Fragment number                  :', MyAtom
        write(DECinfo%output,'(1X,a,i4)')    'FOP Done: Number of orbitals in virt total :', &
           & AtomicFragment%nvirtAOS
        write(DECinfo%output,'(1X,a,i4)')    'FOP Done: Number of orbitals in occ total  :', &
           & AtomicFragment%noccAOS
        write(DECinfo%output,'(1X,a,i4)')     'FOP Done: Number of basis functions        :', &
           & AtomicFragment%nbasis
        if(.not. DECinfo%onlyvirtpart) then
           write(DECinfo%output,'(1X,a,f16.10)') 'FOP Done: Occupied Fragment energy         :',&
              & get_occ_energy_fragopt(AtomicFragment%energies,AtomicFragment%ccmodel)
        endif
        if(.not. DECinfo%onlyoccpart) then
           write(DECinfo%output,'(1X,a,f16.10)') 'FOP Done: Virtual Fragment energy          :',&
              & get_virt_energy_fragopt(AtomicFragment%energies,AtomicFragment%ccmodel)
        endif
        if(.NOT.(DECinfo%onlyoccpart.OR. DECinfo%onlyvirtpart))then
           write(DECinfo%output,'(1X,a,f16.10)') 'FOP Done: Lagrangian Fragment energy       :',&
              & get_lag_energy_fragopt(AtomicFragment%energies,AtomicFragment%ccmodel)
        end if
        write(DECinfo%output,'(1X,a,i6,a,i6,a,f5.2,a)')  'FOP Done: Occupied reduction removed ', &
           & no_exp-AtomicFragment%noccAOS, ' of ', no_exp, ' orbitals ( ', &
           & (no_exp-AtomicFragment%noccAOS)*100.0_realk/no_exp, ' %)'
        write(DECinfo%output,'(1X,a,i6,a,i6,a,f5.2,a)')  'FOP Done: Virtual  reduction removed ', &
           & nv_exp-AtomicFragment%nvirtAOS, ' of ', nv_exp, ' orbitals ( ', &
           & (nv_exp-AtomicFragment%nvirtAOS)*100.0_realk/nv_exp, ' %)'
        write(DECinfo%output,'(1X,a)') 'FOP========================================================='
        write(DECinfo%output,'(1X,a,/)') 'FOP'
     endif


  end subroutine fragment_reduction_procedure


  !> Given a converged atomic fragment using local orbitals, determine fragment
  !> of reduced size, where the energy error compared to the original converged fragment is
  !> below the FOT.
  !> \date September 2013
  !> \author Kasper Kristensen
  subroutine fragopt_reduce_FOs(MyAtom,AtomicFragment,OccOrbitals,nOcc, &
           & virtOrbitals,nvirt,MyMolecule,mylsitem,t2tens,iter)
    implicit none
    !> Number of occupied orbitals in molecule
    integer, intent(in) :: nOcc
    !> Number of virtupied orbitals in molecule
    integer, intent(in) :: nvirt
    !> Central atom in molecule
    integer, intent(in) :: MyAtom
    !> Atomic fragment to be optimized
    type(decfrag),intent(inout)        :: AtomicFragment
    !> All occupied orbitals
    type(decorbital), dimension(nOcc), intent(in)      :: OccOrbitals
    !> All virtupied orbitals
    type(decorbital), dimension(nvirt), intent(in)    :: virtOrbitals
    !> Full molecule information
    type(fullmolecule), intent(in) :: MyMolecule
    !> Integral information
    type(lsitem), intent(inout)       :: mylsitem
    !> Doubles amplitudes for converged fragment (local orbital basis, index order: a,i,b,j)
    !> At output, the virtual indices will have been transformed to the (smaller) FO orbital space
    type(tensor),intent(inout) :: t2tens
    !> Maximum number of reduction steps
    integer,intent(out) :: iter
    integer :: ov,Nold,Nnew,nocc_exp,nvirt_exp
    logical :: reduction_converged,ReductionPossible(2)
    real(realk), dimension(ndecenergies) :: Eold
    real(realk) :: FOT
    type(array4)  :: t2
    character(4) :: stens_atype
    integer :: os,vs,nnod
    nnod=1
#ifdef VAR_MPI
    nnod=infpar%lg_nodtot
#endif

    ! Init stuff
    ! **********
    FOT = DECinfo%FOT
    ! Reference energies for converged fragment
    Eold = AtomicFragment%energies
    ! Dimensions for converged fragment
    nocc_exp = AtomicFragment%noccAOS
    nvirt_exp = AtomicFragment%nvirtAOS


    
    ! ======================================================================
    !                             Reduction loop
    ! ======================================================================
    ! Here we carry out unitary transformation of occupied and virtual AOSs 
    ! to adapt the AOS to the EOS orbitals. The orbitals assigned to the central atom
    ! of course requires special treatment, see calculate_corrdens for details.
    ! We term the new set of orbitals "fragment-adapted orbitals".
    ReductionPossible=.false.


    ! First reduce virtual space (ov=2) then occupied space (ov=1)
    OCC_OR_VIRT: do ov=2,1,-1   

       reduction_converged=.false.

       ! When we have selected a set of virtual FOs, we want to adapt the
       ! occupied correlation density matrix to this set of virtual FOs.
       ! Diagonalization of the occupied correlation density matrix defines a 
       ! set of occupied FOs which describe amplitudes t(a,i,b,j) 
       ! (where a,b are virtual FOs, NOT local orbitals) as compactly as possible.
       ! Note: Only do this for occupied partitioning scheme, otherwise
       !       we are not allowed to mix the virtual orbitals assigned to the central atom
       !       with the other atoms.
       if(ov==1 .and. DECinfo%onlyoccpart) then

          !FIXME: dirty hack
          t2 = array4_init(t2tens%dims)
          call tensor_convert(t2tens, t2%val)
          stens_atype = t2tens%atype
          call tensor_free(t2tens)

          call transform_virt_amp_to_FOs(t2,AtomicFragment)
          call calculate_corrdens_AOS_occocc(t2,AtomicFragment)

          !FIXME: dirty hack-restore input tensor with the new dimensions and data
          call get_symm_tensor_segmenting_simple(nnod,t2%dims(2),t2%dims(1),os,vs)
          call tensor_minit(t2tens,t2%dims,4, atype = stens_atype, tdims=[vs,os,vs,os])
          call tensor_convert(t2%val,t2tens)
          call array4_free(t2)

       end if
       if(DECinfo%onlyVirtpart) then
          WRITE(DECinfo%output,*)'WARNING FOs usign onlyVirtpart not tested'
          call lsquit('Error FOs usign onlyVirtpart not tested',-1)
       endif
       REDUCTION_LOOP: do iter=1,DECinfo%maxiter

          if(ov==1) then
             write(DECinfo%output,*) 'FOP Starting occ reduction step ', iter
          else
             write(DECinfo%output,*) 'FOP Starting virt reduction step ', iter
          end if


          ! Set rejection threshold and some bookkeeping of fragment size
          ! *************************************************************
          if(iter == 1) then
             ! Number of occ or virt orbitals in converged fragment
             if(ov==1) then
                Nold = AtomicFragment%noccAOS
             else
                Nold = AtomicFragment%nvirtAOS
             end if
             ! Threshold for throwing away fragment-adapted orbitals (start with FOT value)
             AtomicFragment%RejectThr(ov) = FOT
          else  
             ! Number of occ or virt fragment-adapted orbitals from previous step
             if(ov==1) then
                Nold = AtomicFragment%noccAOS
             else
                Nold = AtomicFragment%nvirtAOS
             end if
             ! decrease rejection threshold by a factor 10 in each step
             AtomicFragment%RejectThr(ov) = AtomicFragment%RejectThr(ov)/10.0_realk
          end if
          if(ov==1) then
             write(DECinfo%output,'(a,ES13.5)') ' FOP occ rejection threshold  : ',&
                  & AtomicFragment%RejectThr(ov)
          else
             write(DECinfo%output,'(a,ES13.5)') ' FOP virt rejection threshold  : ',&
                  & AtomicFragment%RejectThr(ov)
          end if


          ! Make fragment-adapted fragment with smaller AOS according to rejection threshold
          ! *********************************************************************************
          call fragment_adapted_transformation_matrices(AtomicFragment)
          if(ov==1) then
             Nnew = AtomicFragment%noccAOS
          else
             Nnew = AtomicFragment%nvirtAOS
          end if


          ! Cycle if the AOS was not decreased in the rejection procedure
          ! *************************************************************
          if( Nold==Nnew ) then
             if(iter == DECinfo%maxiter) then
                if(ov==2) then  ! not possible to reduce virtual space, try occupied space
                   cycle OCC_OR_VIRT

                else ! not possible to decrease occ space
                   if (.not. any(ReductionPossible)) then
                      write(DECinfo%output,*) "FOP No reduction possible. Use original &
                         & converged fragment"
                      AtomicFragment%RejectThr=0.0_realk
                   end if
                end if
             else
                write(DECinfo%output,*) 'FOP No reduction occurred - try again...'
                cycle REDUCTION_LOOP
             end if
          end if


          ! Get fragment energy for fragment using FOs
          ! ******************************************
          call fragment_energy_and_prop(AtomicFragment)


          ! Test if fragment energy (or energies) are converged to FOT precision
          ! ********************************************************************
          call get_fragopt_energy_error(AtomicFragment,AtomicFragment%ccmodel,Eold)
          call fragopt_print_info(AtomicFragment,AtomicFragment%ccmodel,iter)
          call fragopt_check_convergence(AtomicFragment%ccmodel,AtomicFragment,FOT,reduction_converged)


          ! Quit if fragment reduction is converged
          ! ***************************************          
          if (reduction_converged) then
             ReductionPossible(ov)=.true.
             if(ov==2) then  ! virtual reduction converged 
                write(DECinfo%output,*) 'FOP Virt reduction of fragment converged in step', iter
                cycle OCC_OR_VIRT  ! try to reduce occ space

             else ! both occ and virtual reductions have converged
                write(DECinfo%output,*) 'FOP Occ  reduction of fragment converged in step', iter
                exit
             end if
          end if


          ! Special case: No reduction possible for neither occ nor virt space
          ! ******************************************************************
          if ( (iter == DECinfo%maxiter) .and. (ov==1) .and. &
               & (.not. any(ReductionPossible) ) ) then
             write(DECinfo%output,*) "FOP No reduction possible. Use original converged fragment"
             AtomicFragment%RejectThr=0.0_realk
          end if

       end do REDUCTION_LOOP

    end do OCC_OR_VIRT
    

    ! Print out info
    ! **************
    write(DECinfo%output,*)'FOP'
    write(DECinfo%output,'(1X,a)') 'FOP========================================================='
    write(DECinfo%output,'(1X,a,i4)') 'FOP    FO REDUCTION HAS CONVERGED FOR SITE',MyAtom
    write(DECinfo%output,'(1X,a)') 'FOP---------------------------------------------------------'
    write(DECinfo%output,'(1X,a,i4)')    'FOP Done: Fragment number                  :', MyAtom
    write(DECinfo%output,'(1X,a,i4)')    'FOP Done: Number of orbitals in virt total :', &
         & AtomicFragment%nvirtFA
    write(DECinfo%output,'(1X,a,i4)')    'FOP Done: Number of orbitals in occ total  :', &
         & AtomicFragment%noccFA
    write(DECinfo%output,'(1X,a,i4)')     'FOP Done: Number of basis functions        :', &
         & AtomicFragment%nbasis
    if(.not. DECinfo%onlyvirtpart) then
       write(DECinfo%output,'(1X,a,f16.10)') 'FOP Done: Occupied Fragment energy         :',&
         & get_occ_energy_fragopt(AtomicFragment%energies,AtomicFragment%ccmodel)
    endif
    if(.not. DECinfo%onlyoccpart) then
       write(DECinfo%output,'(1X,a,f16.10)') 'FOP Done: Virtual Fragment energy          :',&
         & get_virt_energy_fragopt(AtomicFragment%energies,AtomicFragment%ccmodel)
    endif
    if(.NOT.(DECinfo%onlyoccpart.OR. DECinfo%onlyvirtpart))then
       write(DECinfo%output,'(1X,a,f16.10)') 'FOP Done: Lagrangian Fragment energy       :',&
         & get_lag_energy_fragopt(AtomicFragment%energies,AtomicFragment%ccmodel)
    end if
    write(DECinfo%output,'(1X,a,g14.2)')  'FOP Done: Occupied reduction threshold     :', &
         & AtomicFragment%RejectThr(1)
    write(DECinfo%output,'(1X,a,g14.2)')  'FOP Done: Virtual  reduction threshold     :', &
         & AtomicFragment%RejectThr(2)
    write(DECinfo%output,'(1X,a,i6,a,i6,a,f5.2,a)')  'FOP Done: Occupied reduction removed ', &
         & nocc_exp-AtomicFragment%noccFA, ' of ', nocc_exp, ' orbitals ( ', &
         & (nocc_exp-AtomicFragment%noccFA)*100.0_realk/nocc_exp, ' %)'
    write(DECinfo%output,'(1X,a,i6,a,i6,a,f5.2,a)')  'FOP Done: Virtual  reduction removed ', &
         & nvirt_exp-AtomicFragment%nvirtFA, ' of ', nvirt_exp, ' orbitals ( ', &
         & (nvirt_exp-AtomicFragment%nvirtFA)*100.0_realk/nvirt_exp, ' %)'
    write(DECinfo%output,'(1X,a,/)') 'FOP========================================================='
    write(DECinfo%output,*) 'FOP'

  end subroutine fragopt_reduce_FOs


  !> \brief Check if fragment energy (or energies) is converged to FOT precision
  !> \author Kasper Kristensen
  !> \date September 2013
  subroutine fragopt_check_convergence(ccmodel,frag,FOT,converged)
    implicit none
    !> CC model
    integer :: ccmodel
    !> Atomic fragment with energy errors for a given model 
    !  (usually frag%ccmodel but not always !!!)
    type(decfrag), intent(in) :: frag
    !> Fragment optimization threshold
    real(realk),intent(in) :: FOT
    !> Is fragment energy (or energies) converged?
    logical,intent(inout) :: converged
    logical :: lag_converged, occ_converged, virt_converged

    ! Lagrangian
    TEST_CONVERGENCE_LAG: if(DECinfo%OnlyOccPart) then 
       ! do not consider Lagrangian error if we are only interested in occupied partitioning scheme
       ! --> just set Lagrangian to be converged always
       lag_converged=.true.
    elseif(DECinfo%OnlyVirtPart) then 
       ! do not consider Lagrangian error if we are only interested in virtual partitioning scheme
       ! --> just set Lagrangian to be converged always
       lag_converged=.true.
    else
       if  (abs(frag%Elag_err) < FOT) then
          if( DECinfo%print_small_calc )then
             write(DECinfo%output,'(1X,a,F14.9)') 'FOP: Lagrangian energy converged, energydiff =', &
                & abs(frag%Elag_err)
          endif
          lag_converged=.true.
       else
          if( DECinfo%print_small_calc )then
             write(DECinfo%output,*) 'FOP: Lagrangian energy NOT converged'
          endif
          lag_converged=.false.
       end if
    end if TEST_CONVERGENCE_LAG

    ! Occupied 
    TEST_CONVERGENCE_OCC: if(DECinfo%OnlyVirtPart) then
       occ_converged=.true.
    else
       if  (abs(frag%Eocc_err) < FOT) then
          if( DECinfo%print_small_calc )then
             write(DECinfo%output,'(1X,a,F14.9)') 'FOP: Occupied energy converged, energydiff   =',&
                & abs(frag%Eocc_err)
          endif
          occ_converged=.true.
       else
          if( DECinfo%print_small_calc )then
             write(DECinfo%output,*) 'FOP: Occupied energy NOT converged'
          endif
          occ_converged=.false.
       end if
    endif TEST_CONVERGENCE_OCC

    ! Virtual 
    TEST_CONVERGENCE_VIRT: if(DECinfo%OnlyOccPart) then
       virt_converged=.true.
    else
       if  (abs(frag%Evir_err) < FOT) then
          if( DECinfo%print_small_calc )then
             write(DECinfo%output,'(1X,a,F14.9)') 'FOP: Virtual energy converged, energydiff    =',&
                & abs(frag%Evir_err)
          endif
          virt_converged=.true.
       else
          if( DECinfo%print_small_calc )then
             write(DECinfo%output,*) 'FOP: Virtual energy NOT converged'
          endif
          virt_converged=.false.
       end if
    end if TEST_CONVERGENCE_VIRT

    ! We are converged only if ALL three energies are converged
    ExpansionConvergence: if(lag_converged .and. occ_converged .and. virt_converged) then
       converged=.true.
    else
       converged=.false.
    end if ExpansionConvergence


  end subroutine fragopt_check_convergence


  !> Purpose: Check if the last two steps of the binary search are close enough
  !           for the reduction to be stoped.
  !
  !> Author:  Pablo Baudin
  !> Date:    July 2014
  subroutine check_red_space_convergence(Nold,Nnew,Nmin,Nmax,gap,Space1_conv,Space2_conv, &
           & converged,reduce_occ,step_accepted)

     implicit none

     integer, intent(inout) :: Nold, Nnew, Nmin, Nmax
     integer, intent(in)    :: gap
     logical, intent(inout) :: Space1_conv, Space2_conv, converged
     logical, intent(in)    :: reduce_occ, step_accepted

     ! if we have reach a small enough gap (# of orbs) between the 
     ! last two step it means we have converged:
     if (abs(Nold-Nnew) <= gap) then
        Space1_conv = .true.

        ! If the last step was accepted, we keep it:
        if (step_accepted) then
           Nold = Nnew
           Nmin = Nnew
           Nmax = Nnew
        ! Else, we keep the last converged step:
        else
           Nnew = Nmax
           Nold = Nmax
           Nmin = Nmax
        end if

        if( DECinfo%print_small_calc )then
           if (reduce_occ) then
              write(DECinfo%output,*) &
                 & 'BIN SEARCH: OCCUPIED REDUCTION CONVERGED', gap
           else
              write(DECinfo%output,*) &
                 & 'BIN SEARCH: VIRTUAL REDUCTION CONVERGED', gap
           end if
        endif

        if (Space2_conv) converged = .true.
     end if

  end subroutine check_red_space_convergence


  !> Fragment optimization, special case where we want to include full molecular system
  ! in fragment orbital space (only for debugging purposes).
  !> \date September 2013
  !> \author Kasper Kristensen
  subroutine fragopt_include_fullmolecule(MyAtom,AtomicFragment, &
       &OccOrbitals,nOcc,virtOrbitals,nvirt, &
       &MyMolecule,mylsitem,freebasisinfo,t1full)
    implicit none
    !> Number of occupied orbitals in molecule
    integer, intent(in) :: nOcc
    !> Number of virtupied orbitals in molecule
    integer, intent(in) :: nvirt
    !> Central atom in molecule
    integer, intent(in) :: MyAtom
    !> Atomic fragment to be optimized
    type(decfrag),intent(inout)        :: AtomicFragment
    !> All occupied orbitals
    type(decorbital), dimension(nOcc), intent(in)      :: OccOrbitals
    !> All virtupied orbitals
    type(decorbital), dimension(nvirt), intent(in)    :: virtOrbitals
    !> Full molecule information
    type(fullmolecule), intent(inout) :: MyMolecule
    !> Integral information
    type(lsitem), intent(inout)       :: mylsitem
    !> Delete fragment basis information ("expensive box in decfrag type") at exit?
    logical,intent(in) :: freebasisinfo
    !> t1 amplitudes for full molecule to be updated (only used when DECinfo%SinglesPolari is set)
    type(array2),intent(inout),optional :: t1full
    type(tensor) :: g, t2


    write(DECinfo%output,*) 'FOP Fragment includes all orbitals and fragment optimization is skipped'

    ! Init fragment will orbital space for full molecule
    call fragment_init_simulate_full(MyAtom,nvirt, nocc, OccOrbitals,virtOrbitals,&
         & MyMolecule,mylsitem,AtomicFragment,.true.)

    ! Set information for fragment-adapted orbitals
    FullFragAdapt: if(DECinfo%fragadapt) then

       ! Set model to be MP2
       MyMolecule%ccmodel(MyAtom,MyAtom) = MODEL_MP2

       ! Integrals (ai|bj)
       !call get_VOVO_integrals(AtomicFragment%mylsitem,AtomicFragment%nbasis,&
       !     & AtomicFragment%noccAOS,AtomicFragment%nvirtAOS,&
       !     & AtomicFragment%Cv, AtomicFragment%Co, g)
       ! MP2 amplitudes
       !call mp2_solver(AtomicFragment%noccAOS,AtomicFragment%nvirtAOS,&
       !     & AtomicFragment%ppfock,AtomicFragment%qqfock,g,t2)
       call mp2_solver(AtomicFragment,g,t2,.false.)
       call tensor_free(g)
       ! Get correlation density matrix
       call calculate_MP2corrdens_frag(t2,AtomicFragment)
       call tensor_free(t2)

       ! Set MO coefficient matrices corresponding to fragment-adapted orbitals
       call fragment_adapted_transformation_matrices(AtomicFragment)

       ! Restore the original CC model 
       MyMolecule%ccmodel(MyAtom,MyAtom) = DECinfo%ccmodel

    end if FullFragAdapt

    if(freebasisinfo) then
       call atomic_fragment_free_basis_info(AtomicFragment)
    end if
    AtomicFragment%isopt=.true.

  end subroutine fragopt_include_fullmolecule


  !> \brief prints info for atomic fragment in a given iteration of the fragment optimization
  !> using the Lagrangian partitioning scheme.
  !> \date: august-2011
  !> \author: Ida-Marie Hoeyvik
  subroutine fragopt_print_info(frag,ccmodel,iter,add_mp2_opt)
    implicit none
    !> Atomic fragment to be optimized
    type(decfrag),intent(in) :: frag
    !> CC model
    integer, intent(in) :: ccmodel
    !> extra mp2?
    logical, intent(in), optional :: add_mp2_opt
    !> Iteration number
    integer,intent(in) :: iter
    logical :: mp2

    mp2 = .false.
    if (present(add_mp2_opt)) mp2 = add_mp2_opt

    if(DECinfo%print_small_calc)then
       write(DECinfo%output,*)'FOP'
       if (mp2) write(DECinfo%output,*) 'FOP              MP2 INFORMATION (also requested):'
       write(DECinfo%output,*)'FOP'
       write(DECinfo%output,'(1X,a)') 'FOP========================================================='
       write(DECinfo%output,'(1X,a,i4)') 'FOP              Fragment information, loop', iter
       write(DECinfo%output,'(1X,a)') 'FOP---------------------------------------------------------'
       write(DECinfo%output,'(1X,a,i4)')    'FOP Loop: Fragment number                  :', &
          & frag%EOSatoms(1)
       write(DECinfo%output,'(1X,a,i4)')    'FOP Loop: Number of orbitals in virt total :', &
          & Frag%nvirtAOS
       write(DECinfo%output,'(1X,a,i4)')    'FOP Loop: Number of orbitals in occ total  :', &
          & Frag%noccAOS
       write(DECinfo%output,'(1X,a,i4)')    'FOP Loop: Number of basis functions        :', &
          & Frag%nbasis
       if(.not. DECinfo%OnlyVirtPart) then
          write(DECinfo%output,'(1X,a,f16.10)') 'FOP Loop: Occupied Fragment energy         :',&
             & get_occ_energy_fragopt(frag%energies,ccmodel)
       endif
       if(.not. DECinfo%OnlyOccPart) then
          write(DECinfo%output,'(1X,a,f16.10)') 'FOP Loop: Virtual Fragment energy          :',&
             & get_virt_energy_fragopt(frag%energies,ccmodel)
       endif
       if(.not. (DECinfo%OnlyOccPart.or.DECinfo%OnlyVirtPart)) then
          write(DECinfo%output,'(1X,a,f16.10)') 'FOP Loop: Lagrangian Fragment energy       :',&
             & get_lag_energy_fragopt(frag%energies,ccmodel)
       end if

       PrintDiff: if(iter/=0) then ! no point in printing energy differences for initial energy
          if(.not. DECinfo%OnlyVirtPart) then
             write(DECinfo%output,'(1X,a,f16.10)') 'FOP Loop: Occupied energy diff             :',&
               & abs(frag%Eocc_err)
          endif
          if(.not. DECinfo%OnlyOccPart) then
             write(DECinfo%output,'(1X,a,f16.10)') 'FOP Loop: Virtual energy diff              :',&
               & abs(frag%Evir_err)
          endif
          if(.not. (DECinfo%OnlyOccPart.or.DECinfo%OnlyVirtPart)) then
             write(DECinfo%output,'(1X,a,f16.10)') 'FOP Loop: Lagrangian energy diff           :',&
               & abs(frag%Elag_err)
          end if
       end if PrintDiff

       write(DECinfo%output,'(1X,a)') 'FOP========================================================='
       write(DECinfo%output,*) 'FOP'
    endif


  end subroutine fragopt_print_info


  !> Purpose: For a given model, get the corresponding occupied fragment energy
  !> Author:  Pablo Baudin
  !> Date:    April 2015
  function get_occ_energy_fragopt(energies,ccmodel) result(Eocc)
    implicit none
    real(realk) :: Eocc
    real(realk), dimension(ndecenergies), intent(in) :: energies
    integer, intent(in) :: ccmodel
    ! MODIFY FOR NEW MODEL
    ! If you implement a new model, please copy fragment%energies(?) for your model,
    ! see decfrag type def to determine the "?".

    Eocc = 0.0_realk

    select case(ccmodel)
    case(MODEL_MP2)
       ! MP2
       Eocc = energies(FRAGMODEL_OCCMP2)
    case(MODEL_CC2)
       ! CC2
       Eocc = energies(FRAGMODEL_OCCCC2)
    case(MODEL_RPA)
       ! RPA
       Eocc = energies(FRAGMODEL_OCCRPA)
    case(MODEL_SOSEX)
       ! RPA
       Eocc = energies(FRAGMODEL_OCCSOS)
    case(MODEL_CCSD)
       ! CCSD
       Eocc = energies(FRAGMODEL_OCCCCSD)
#ifdef MOD_UNRELEASED
    case(MODEL_CCSDpT)
       ! CCSD(T): CCSD contribution + (T) contribution
       Eocc = energies(FRAGMODEL_OCCCCSD) + energies(FRAGMODEL_OCCpT)
#endif
    case(MODEL_RIMP2)
       ! RI-MP2
       Eocc = energies(FRAGMODEL_OCCRIMP2)
    case(MODEL_LSTHCRIMP2)
       ! LS-THC-RI-MP2
       Eocc = energies(FRAGMODEL_OCCLSTHCRIMP2)
    case default
       write(DECinfo%output,*) 'WARNING: get_occ_energy_fragopt needs implementation &
            & for model:', ccmodel
    end select

  end function get_occ_energy_fragopt


  !> Purpose: For a given model, get the corresponding virtual fragment energy
  !> Author:  Pablo Baudin
  !> Date:    April 2015
  function get_virt_energy_fragopt(energies,ccmodel) result(Evirt)
    implicit none
    real(realk) :: Evirt
    real(realk), dimension(ndecenergies), intent(in) :: energies
    integer, intent(in) :: ccmodel
    ! MODIFY FOR NEW MODEL
    ! If you implement a new model, please copy fragment%energies(?) for your model,
    ! see decfrag type def to determine the "?".

    Evirt = 0.0_realk

    select case(ccmodel)
    case(MODEL_MP2)
       ! MP2
       Evirt = energies(FRAGMODEL_VIRTMP2)
    case(MODEL_CC2)
       ! CC2
       Evirt = energies(FRAGMODEL_VIRTCC2)
    case(MODEL_RPA)
       ! RPA
       Evirt = energies(FRAGMODEL_VIRTRPA)
    case(MODEL_SOSEX)
       ! RPA
       Evirt = energies(FRAGMODEL_VIRTSOS)
    case(MODEL_CCSD)
       ! CCSD
       Evirt = energies(FRAGMODEL_VIRTCCSD)
#ifdef MOD_UNRELEASED
    case(MODEL_CCSDpT)
       ! CCSD(T): CCSD contribution + (T) contribution
       Evirt = energies(FRAGMODEL_VIRTCCSD) + energies(FRAGMODEL_VIRTpT)
#endif
    case(MODEL_RIMP2)
       ! RI-MP2
       Evirt = energies(FRAGMODEL_VIRTRIMP2)
    case(MODEL_LSTHCRIMP2)
       ! LS-THC-RI-MP2
       Evirt = energies(FRAGMODEL_VIRTLSTHCRIMP2)
    case default
       write(DECinfo%output,*) 'WARNING: get_virt_energy_fragopt needs implementation &
            & for model:', ccmodel
    end select

  end function get_virt_energy_fragopt


  !> Purpose: For a given model, get the corresponding lagrangian fragment energy
  !> Author:  Pablo Baudin
  !> Date:    April 2015
  function get_lag_energy_fragopt(energies,ccmodel) result(Elag)
    implicit none
    real(realk) :: Elag
    real(realk), dimension(ndecenergies), intent(in) :: energies
    integer, intent(in) :: ccmodel
    ! MODIFY FOR NEW MODEL
    ! If you implement a new model, please copy fragment%energies(?) for your model,
    ! see decfrag type def to determine the "?".

    Elag = 0.0_realk

    select case(ccmodel)
    case(MODEL_MP2)
       ! MP2
       Elag = energies(FRAGMODEL_OCCMP2)
    case(MODEL_CC2)
       ! CC2
       Elag = energies(FRAGMODEL_OCCCC2)
    case(MODEL_RPA)
       ! RPA
       Elag = energies(FRAGMODEL_OCCRPA)
    case(MODEL_SOSEX)
       ! RPA
       Elag = energies(FRAGMODEL_OCCSOS)
    case(MODEL_CCSD)
       ! CCSD
       Elag = energies(FRAGMODEL_OCCCCSD)
#ifdef MOD_UNRELEASED
    case(MODEL_CCSDpT)
       ! CCSD(T): CCSD contribution + (T) contribution
       Elag = energies(FRAGMODEL_OCCCCSD) + energies(FRAGMODEL_OCCpT)
#endif
    case(MODEL_RIMP2)
       ! RI-MP2
       Elag = energies(FRAGMODEL_OCCRIMP2)
    case(MODEL_LSTHCRIMP2)
       ! LS-THC-RI-MP2
       Elag = energies(FRAGMODEL_OCCLSTHCRIMP2)
    case default
       write(DECinfo%output,*) 'WARNING: get_lag_energy_fragopt needs implementation &
            & for model:', ccmodel
    end select

  end function get_lag_energy_fragopt

  
  !> Purpose: Calculate fragment energy error for a given model
  !
  !> Author:  Pablo Baudin
  !> Date:    April 2015
  subroutine get_fragopt_energy_error(frag,ccmodel,Eref)
     implicit none
     !> Atomic fragment to optimized
     type(decfrag), intent(inout) :: frag
     !> CC model
     integer, intent(in) :: ccmodel
     !> Reference energies
     real(realk), dimension(ndecenergies), intent(in) :: Eref

      frag%Elag_err = (get_lag_energy_fragopt(Eref,ccmodel) - &
         & get_lag_energy_fragopt(frag%energies,ccmodel))
      frag%Eocc_err = (get_occ_energy_fragopt(Eref,ccmodel) - &
         & get_occ_energy_fragopt(frag%energies,ccmodel))
      frag%Evir_err = (get_virt_energy_fragopt(Eref,ccmodel) - &
         & get_virt_energy_fragopt(frag%energies,ccmodel))

   end subroutine get_fragopt_energy_error


  !> \brief When fragment energies have been calculated, put them
  !> into the energies array in fragment structure according to the given model
  !> (see FRAGMODEL_* in dec_typedef.F90).
  !> \author Kasper Kristensen
  !> \date October 2013
  subroutine put_fragment_energy_contribs_wrapper(Eocc,Evirt,MyFragment,doSOS)
    implicit none
    !> Occupied and virtual partitioning scheme energies
    real(realk),intent(in) :: Eocc, Evirt
    !> Atomic or pair fragment
    type(decfrag),intent(inout) :: MyFragment
    !> SOS cont
    logical,intent(in),optional :: doSOS
    logical :: SOS

    !initalize SOS
    SOS = .false.
    if(present(doSOS)) SOS = doSOS

    ! Put energies into their proper place in the MyFragment%energies array
    ! according to the CC model used for the fragment
    call put_fragment_energy_contribs(MyFragment%ccmodel,Eocc,Evirt,MyFragment%energies,SOS)


    ! Special case: When some pairs are treated at the MP2 level, while the
    ! actual CC model is higher in the hierarchy (e.g. CCSD), we need to 
    ! put the MP2 fragment energies into the energy array for the more accurate model.
    if(MyFragment%ccmodel==MODEL_MP2 .and. DECinfo%ccmodel/=MODEL_MP2) then
       call put_fragment_energy_contribs(DECinfo%ccmodel,Eocc,Evirt,MyFragment%energies,SOS)
    end if

  end subroutine put_fragment_energy_contribs_wrapper


  ! MODIFY FOR NEW MODEL
  !> \brief When fragment energies have been calculated, put them
  !> into the energies array according to the given model
  !> (see FRAGMODEL_* in dec_typedef.F90).
  !> \author Kasper Kristensen
  !> \date October 2013
  subroutine put_fragment_energy_contribs(ccmodel,Eocc,Evirt,energies,SOS)
    implicit none
    !> Which CC model
    integer,intent(in) :: ccmodel
    !> Occupied and virtual partitioning scheme energies
    real(realk),intent(in) :: Eocc, Evirt
    !> Energies array
    real(realk),intent(inout) :: energies(ndecenergies)
    !> SOS cont
    logical,intent(in) :: SOS

    ! Put energies into their proper place in the energies array
    select case(ccmodel)
    case(MODEL_MP2)
       ! MP2
       energies(FRAGMODEL_OCCMP2) = Eocc   ! occupied
       energies(FRAGMODEL_VIRTMP2) = Evirt   ! virtual
    case(MODEL_CC2)
       ! CC2
       energies(FRAGMODEL_OCCCC2) = Eocc   ! occupied
       energies(FRAGMODEL_VIRTCC2) = Evirt   ! virtual
    case(MODEL_RPA)
       ! RPA
       if(SOS) then
         energies(FRAGMODEL_OCCSOS) = Eocc   ! occupied
         energies(FRAGMODEL_VIRTSOS) = Evirt   ! virtual
       else
         energies(FRAGMODEL_OCCRPA) = Eocc   ! occupied
         energies(FRAGMODEL_VIRTRPA) = Evirt   ! virtual
       endif
    case(MODEL_SOSEX)
       !SOSEX
       energies(FRAGMODEL_OCCSOS) = Eocc   ! occupied
       energies(FRAGMODEL_VIRTSOS) = Evirt   ! virtual
    case(MODEL_CCSD)
       ! CCSD
       energies(FRAGMODEL_OCCCCSD) = Eocc   ! occupied
       energies(FRAGMODEL_VIRTCCSD) = Evirt   ! virtual
#ifdef MOD_UNRELEASED
    case(MODEL_CCSDpT)
       ! Save also CCSD contribution for CCSD(T)
       energies(FRAGMODEL_OCCCCSD) = Eocc   ! occupied
       energies(FRAGMODEL_VIRTCCSD) = Evirt   ! virtual
       !endif mod_unreleased
#endif
    case(MODEL_RIMP2)
       ! RI-MP2
       energies(FRAGMODEL_OCCRIMP2) = Eocc     ! occupied
       energies(FRAGMODEL_VIRTRIMP2) = Evirt   ! virtual
    case(MODEL_LSTHCRIMP2)
       ! LS-THC-RI-MP2
       energies(FRAGMODEL_OCCLSTHCRIMP2) = Eocc     ! occupied
       energies(FRAGMODEL_VIRTLSTHCRIMP2) = Evirt   ! virtual
    case default
       call lsquit('case unknown in put_fragment_energy_contribs',-1)
    end select

  end subroutine put_fragment_energy_contribs


end module fragment_energy_module
