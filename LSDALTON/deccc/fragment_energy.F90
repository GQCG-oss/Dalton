!> @file
!> Module contains operations involving fragment energy calculations
!> \author Marcin Ziolkowski and Kasper Kristensen

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
#ifdef MOD_UNRELEASED
  use f12_integrals_module
#endif
  use dec_fragment_utils
  use array2_simple_operations
  use array4_simple_operations
  use orbital_operations
  use mp2_module !,only: max_batch_dimension,get_vovo_integrals, &
!       & mp2_integrals_and_amplitudes
  use atomic_fragment_operations!  ,only: atomic_fragment_init_basis_part, &
  !       & get_fragmentt1_AOSAOS_from_full, extract_specific_fragmentt1, &
  !       & update_full_t1_from_atomic_frag,which_pairs, &
  !       & update_full_t1_from_atomic_frag, update_full_t1_from_pair_frag
#ifdef MOD_UNRELEASED
  use ccsdpt_module, only:ccsdpt_driver,ccsdpt_energy_e4_frag,ccsdpt_energy_e5_frag,&
       ccsdpt_energy_e4_pair,ccsdpt_energy_e5_pair
!endif mod_unreleased
#endif
  use mp2_gradient_module ,only: single_calculate_mp2gradient_driver,&
       & pair_calculate_mp2gradient_driver
  use ccdriver, only: mp2_solver,fragment_ccsolver

public :: optimize_atomic_fragment, pair_driver_singles, atomic_driver, &
     & pair_driver,atomic_driver_advanced,plot_pair_energies,&
     & Full_DECMP2_calculation,estimate_energy_error
private

contains

  !> \brief Construct new atomic fragment based on info in occ_atoms and Unocc_atoms,
  !> and calculate fragment energy. Energy contributions from each individual orbital
  !> is also calculated and stored in MyFragment%OccContribs and MyFragment%VirtContribs for
  !> the occupied and virtual orbitals, respectively.
  !> \author Kasper Kristensen
  !> \date January 2011
  subroutine get_fragment_and_Energy(MyAtom,natoms,occ_atoms,Unocc_atoms,&
       & MyMolecule,MyLsitem,nocc_tot,nunocc_tot,OccOrbitals,UnoccOrbitals,&
       & MyFragment)

    implicit none
    !> Central atom in fragment
    integer, intent(in) :: MyAtom
    !> Number of atoms in the molecule
    integer,intent(in) :: natoms
    !> Which atoms are in the occupied fragment space
    logical, dimension(natoms),intent(in) :: occ_atoms
    !> Which atoms are in the unoccupied fragment space
    logical, dimension(natoms),intent(in) :: Unocc_atoms
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> Full molecule lsitem
    type(lsitem), intent(inout) :: mylsitem
    !> Number of occupied orbitals in full molecule
    integer, intent(in) :: nocc_tot
    !> Number of unoccupied orbitals in full molecule
    integer, intent(in) :: nunocc_tot
    !> Occupied orbitals for full molecule
    type(decorbital), dimension(nOcc_tot), intent(in) :: OccOrbitals
    !> Unoccupied orbitals for full molecule
    type(decorbital), dimension(nUnocc_tot), intent(in) :: UnoccOrbitals
    !> Atomic fragment to be determined  (NOT pair fragment)
    type(decfrag), intent(inout) :: MyFragment
    real(realk) :: tcpu, twall
    logical :: DoBasis
    logical :: ForcePrint
    ForcePrint=.TRUE.

    ! **************************************************************
    ! *                       RUN CALCULATION                      *
    ! **************************************************************
    call LSTIMER('START',tcpu,twall,DECinfo%output)

    ! Initialize fragment
    ! *******************
    DoBasis=.true.
    call atomic_fragment_init_atom_specific(MyAtom,natoms,Unocc_atoms, &
         & occ_atoms,nocc_tot,nunocc_tot,OccOrbitals,UnoccOrbitals, &
         & MyMolecule,mylsitem,MyFragment,DoBasis,.false.)

    ! Calculate fragment energies
    ! ***************************
    call fragment_Energy_and_prop(MyFragment)
    call LSTIMER('FRAG: L.ENERGY',tcpu,twall,DECinfo%output,ForcePrint)

  end subroutine get_fragment_and_Energy

  !> \brief Construct new fragment based on list of orbitals in OccAOS and UnoccAOS,
  !> and calculate fragment energy. 
  !> \author Kasper Kristensen
  !> \date March 2013
  subroutine get_fragment_and_Energy_orb_specific(MyAtom,MyMolecule,mylsitem,&
       & OccOrbitals,UnoccOrbitals,OccAOS,UnoccAOS, MyFragment)

    implicit none
    !> Central atom in fragment
    integer, intent(in) :: MyAtom
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> Full molecule lsitem
    type(lsitem), intent(inout) :: mylsitem
    !> Occupied orbitals for full molecule
    type(decorbital), dimension(MyMolecule%nocc), intent(in) :: OccOrbitals
    !> Unoccupied orbitals for full molecule
    type(decorbital), dimension(MyMolecule%nunocc), intent(in) :: UnoccOrbitals
    !> Which atoms are in the occupied fragment space
    logical, dimension(MyMolecule%nocc),intent(inout) :: OccAOS
    !> Which atoms are in the unoccupied fragment space
    logical, dimension(MyMolecule%nunocc),intent(in) :: UnoccAOS
    !> Atomi Fragment to be determined (NOT pair fragment)
    type(decfrag), intent(inout) :: MyFragment
    real(realk) :: tcpu, twall
    logical :: DoBasis


    ! **************************************************************
    ! *                       RUN CALCULATION                      *
    ! **************************************************************
    call LSTIMER('START',tcpu,twall,DECinfo%output)


    ! Init fragment type based on logical vectors
    DoBasis=.true.
    call atomic_fragment_init_orbital_specific(MyAtom,MyMolecule%nunocc, &
         & MyMolecule%nocc, UnoccAOS, &
         & occAOS,OccOrbitals,UnoccOrbitals,MyMolecule,mylsitem,MyFragment,DoBasis,.false.)

    ! Calculate fragment energies
    call fragment_energy_and_prop(MyFragment)

    call LSTIMER('FRAG: L.ENERGY',tcpu,twall,DECinfo%output)

  end subroutine get_fragment_and_Energy_orb_specific


  !> \brief Wrapper for atomic_driver with the following special features:
  !> 1. It is assumed that the input fragment has been initialized but that
  !> the fragment basis information (expensive box in decfrag type) has not been set.
  !> 2. The fragment basis information is calculated here and then freed again.
  !> 3. This wrapper can also attach exisiting full molecular singles amplitudes
  !>    to the fragment structure and update new improved full molecular singles amplitudes
  !>    by the calculated fragment singles amplitudes.
  !> \author Kasper Kristensen
  !> \date March 2012
  subroutine atomic_driver_advanced(nocc,nunocc,&
       & OccOrbitals,UnoccOrbitals,MyLsitem,MyMolecule,MyFragment,grad,&
       & t1old,t1new)

    implicit none
    !> Atomic fragment 
    type(decfrag), intent(inout) :: myfragment
    !> MP2 gradient structure (only calculated if DECinfo%first_order is turned on)
    type(mp2grad),intent(inout) :: grad
    !> Number of occupied orbitals in full molecule
    integer,intent(in) :: nOcc
    !> Number of unoccupied orbitals in full molecule
    integer,intent(in) :: nUnocc
    !> Occupied orbitals for full molecule
    type(decorbital), intent(in) :: OccOrbitals(nocc)
    !> Unoccupied orbitals for full molecule
    type(decorbital), intent(in) :: UnoccOrbitals(nunocc)
    !> Full molecule info (will be unchanged at output)
    type(fullmolecule), intent(inout) :: MyMolecule
    !> LSDalton info
    type(lsitem), intent(inout) :: mylsitem
    !> Existing full molecular t1 amplitudes
    type(array2),intent(in),optional :: t1old
    !> New full molecular t1 amplitudes which will be updated with the contribution from "MyFragment"
    type(array2),intent(inout),optional :: t1new

    ! Init fragment basis information
    call atomic_fragment_init_basis_part(nunocc, nocc, OccOrbitals,&
         & UnoccOrbitals,MyMolecule,mylsitem,MyFragment)
     
    ! Attach fragments singles amplitudes to fragment structure
    ! if long-range singles polarization effects are requested.
    if(DECinfo%SinglesPolari) then
       call get_fragmentt1_AOSAOS_from_full(MyFragment,t1old)
    end if

    ! Call main driver to get energy (and possibly density or gradient)
    call atomic_driver(MyMolecule,mylsitem,OccOrbitals,UnoccOrbitals,&
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
  subroutine atomic_driver(MyMolecule,mylsitem,OccOrbitals,UnoccOrbitals,&
       & MyFragment,grad)

    implicit none
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> LS item info                                                                                    
    type(lsitem), intent(inout) :: mylsitem
    !> Information about DEC occupied orbitals                                                         
    type(decorbital), dimension(MyMolecule%nocc), intent(in) :: OccOrbitals
    !> Information about DEC unoccupied orbitals
    type(decorbital), dimension(MyMolecule%nunocc), intent(in) :: UnoccOrbitals
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
  subroutine fragment_energy_and_prop(MyFragment,Fragment1,Fragment2,grad)

    implicit none
    !> Atomic or pair fragment
    type(decfrag), intent(inout) :: myfragment
    !> Fragments 1 and 2 used to form pair fragment - only used if myfragment is pair fragment
    type(decfrag), intent(in),optional :: Fragment1, Fragment2
    !> MP2 gradient structure (only calculated if DECinfo%first_order is turned on)
    type(mp2grad),intent(inout),optional :: grad
    type(tensor) :: t1, ccsdpt_t1, m1
    type(tensor) :: VOVO,VOVOocc,VOVOvirt,t2occ,t2virt,VOOO,VOVV,t2,u,VOVOvirtTMP,ccsdpt_t2,m2
    type(tensor) :: t2_occEOS
    real(realk) :: tcpu, twall,debugenergy
    ! timings are allocated and deallocated behind the curtains
    real(realk),pointer :: times_ccsd(:), times_pt(:)
    logical :: print_frags,abc,pair
    type(tensor) :: t2f_local, VOVO_local
    integer :: a,b,i,j,k,l

    ! Pairfragment?
    pair = MyFragment%pairfrag
    if(pair) then
       if( (.not. present(Fragment1)) .or. (.not. present(Fragment2)) ) then
          call lsquit('fragment_energy_and_prop: Missing arguments for pair fragment!',-1)
       end if
    end if

    times_ccsd => null()
    times_pt   => null()
    call LSTIMER('START',tcpu,twall,DECinfo%output)


    ! Which model? MP2,CC2, CCSD etc.
    ! *******************************
    !MODIFY FOR NEW MODEL
    WhichCCmodel: select case(MyFragment%ccmodel)

    case(MODEL_NONE) ! SKip calculation

       return

    case(MODEL_MP2) ! MP2 calculation

       if(DECinfo%first_order) then  ! calculate also MP2 density integrals
          call MP2_integrals_and_amplitudes(MyFragment,VOVOocc,t2occ,VOVOvirt,t2virt,VOOO,VOVV)
       else ! calculate only MP2 energy integrals and MP2 amplitudes
          call MP2_integrals_and_amplitudes(MyFragment,VOVOocc,t2occ,VOVOvirt,t2virt)
       end if

#ifdef MOD_UNRELEASED
       ! MP2-F12 Code
       if(DECinfo%F12) then    
          if(pair) then
             call get_f12_fragment_energy(MyFragment, t2occ%elm4, t1%elm2, MyFragment%ccmodel,&
                  &  Fragment1, Fragment2)
          else
             call get_f12_fragment_energy(MyFragment, t2occ%elm4, t1%elm2, MyFragment%ccmodel)   ! WANGY t2_justdoublesEOS
          end if
             !> Free cabs after each calculation
             call free_cabs()
       endif
#endif

    case(MODEL_RIMP2) ! RIMP2 calculation

       if(DECinfo%first_order)call lsquit('no first order RIMP2',-1)
       call RIMP2_integrals_and_amplitudes(MyFragment,VOVOocc,t2occ,VOVOvirt,t2virt)

    case(MODEL_CC2,MODEL_CCSD,MODEL_CCSDpT,MODEL_RPA) ! higher order CC (-like)

       call dec_fragment_time_init(times_ccsd)

       ! Solve CC equation to calculate amplitudes and integrals 
       ! *******************************************************
       ! Here all output indices in t1,t2, and VOVO are AOS indices.
       ! calculate also MP2 density integrals
#ifdef MOD_UNRELEASED
       if(DECinfo%first_order) then  
          call fragment_ccsolver(MyFragment,t1,t2,VOVO,m1=m1,m2=m2)
       else
#endif
          call fragment_ccsolver(MyFragment,t1,t2,VOVO)
#ifdef MOD_UNRELEASED
       endif
#endif

       ! Extract EOS indices for integrals
       ! *********************************
       call tensor_extract_eos_indices(VOVO,MyFragment,tensor_occEOS=VOVOocc,tensor_virtEOS=VOVOvirt)

#ifdef MOD_UNRELEASED
       if(DECinfo%first_order) then
          !          call z_vec_rhs_ccsd(MyFragment,eta,m1_arr=m1,m2_arr=m2,t1_arr=t1,t2_arr=t2)
          call tensor_free(m2)
          call tensor_free(m1)
       endif
#endif

       ! Calculate combined single+doubles amplitudes
       ! ********************************************
       ! u(a,i,b,j) = t2(a,i,b,j) + t1(a,i)*t1(b,j)          
       call get_combined_SingleDouble_amplitudes(t1,t2,u)

       ! Extract EOS indices for amplitudes
       ! **********************************
       call tensor_extract_eos_indices(u,MyFragment,tensor_occEOS=t2occ,tensor_virtEOS=t2virt)
       ! Note, t2occ and t2virt also contain singles contributions
       call tensor_free(u)

       call dec_fragment_time_get(times_ccsd)

#ifdef MOD_UNRELEASED
       ! calculate ccsd(t) fragment energies
       ! ***********************************

       ! we calculate (T) contribution to single  fragment energy 
       ! and store in MyFragment%energies(FRAGMODEL_OCCpT) and MyFragment%energies(FRAGMODEL_VIRTpT)
       if(MyFragment%ccmodel==MODEL_CCSDpT) then

          !FIXME: (T) DEPENDING ON MANY V^2O^2 ALLOCATIONS ON A LOCAL NODE
          !THIS IS A WORKAROUND:
          call tensor_init(VOVO_local,VOVO%dims,4)
          call tensor_init(t2f_local,t2%dims,4)

          call tensor_add(VOVO_local,1.0E0_realk, VOVO, a = 0.0E0_realk ) 
          call tensor_add(t2f_local, 1.0E0_realk, t2,   a = 0.0E0_realk ) 

          call dec_fragment_time_init(times_pt)

          print_frags = DECinfo%print_frags
!!! this should be decided based on the amount of memory available !!!
          abc = DECinfo%abc

          ! init ccsd(t) singles and ccsd(t) doubles (*T1 and *T2)
          if (abc) then

             call tensor_reorder(VOVO_local,[2,4,1,3]) ! vovo integrals in the order (i,j,a,b)
             call tensor_reorder(t2f_local,[2,4,1,3]) ! ccsd_doubles in the order (i,j,a,b)

             call tensor_init(ccsdpt_t1,[MyFragment%noccAOS,MyFragment%nunoccAOS],2)
             call tensor_init(ccsdpt_t2,[MyFragment%noccAOS,MyFragment%noccAOS,&
                  &MyFragment%nunoccAOS,MyFragment%nunoccAOS],4)

          else

             call tensor_reorder(VOVO_local,[1,3,2,4]) ! vovo integrals in the order (a,b,i,j)
             call tensor_reorder(t2f_local,[1,3,2,4]) ! ccsd_doubles in the order (a,b,i,j)

             call tensor_init(ccsdpt_t1, [MyFragment%nunoccAOS,MyFragment%noccAOS],2)
             call tensor_init(ccsdpt_t2, [MyFragment%nunoccAOS,MyFragment%nunoccAOS,&
                  &MyFragment%noccAOS,MyFragment%noccAOS],4)

          endif

          ! call ccsd(t) driver and single fragment evaluation
          call ccsdpt_driver(MyFragment%noccAOS,MyFragment%nunoccAOS,&
               & MyFragment%nbasis,MyFragment%ppfock,&
               & MyFragment%qqfock,MyFragment%Co,&
               & MyFragment%Cv,MyFragment%mylsitem,&
               & VOVO_local,t2f_local,ccsdpt_t1,print_frags,abc,ccsdpt_doubles=ccsdpt_t2)
          if (abc) then

             call tensor_reorder(t2f_local,[3,4,1,2]) ! ccsd_doubles in the order (a,b,i,j)
             call tensor_reorder(ccsdpt_t2,[3,4,1,2]) ! ccsdpt_doubles in the order (a,b,i,j)
             call tensor_reorder(ccsdpt_t1,[2,1]) ! ccsdpt_singles in the order (a,i)

          endif

          if(pair) then
             call ccsdpt_energy_e4_pair(Fragment1,Fragment2,MyFragment,t2f_local,ccsdpt_t2)
             call ccsdpt_energy_e5_pair(MyFragment,t1,ccsdpt_t1)
          else
             call ccsdpt_energy_e4_frag(MyFragment,t2f_local,ccsdpt_t2,&
                  & MyFragment%OccContribs,MyFragment%VirtContribs)
             call ccsdpt_energy_e5_frag(MyFragment,t1,ccsdpt_t1)
          end if

          ! release ccsd(t) singles and doubles amplitudes
          call tensor_free(ccsdpt_t1)
          call tensor_free(ccsdpt_t2)

          call dec_fragment_time_get(times_pt)

          call tensor_free(VOVO_local)
          call tensor_free(t2f_local)
       end if
#endif 

#ifdef MOD_UNRELEASED
       ! CCSD-F12 Code
       if(DECinfo%F12) then
          call tensor_extract_eos_indices(t2,MyFragment,tensor_occEOS=t2_occEOS)
          call get_f12_fragment_energy(MyFragment, t2_occEOS%elm4, t1%elm2, MyFragment%ccmodel)   ! WANGY t2_justdoublesEOS

          !> Free cabs after each calculation
          call tensor_free(t2_occEOS)
          call free_cabs()
       endif
#endif
       ! free vovo integrals
       call tensor_free(VOVO)

       if(DECinfo%use_singles)then
          call tensor_free(t1)
       endif
       call tensor_free(t2)


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

    ! For frozen core and first order properties we need to remove core indices from VOVOvirt, 
    ! since they, "rather incoveniently", are required for the gradient but not for the energy
    if(DECinfo%frozencore .and. DECinfo%first_order) then

       call remove_core_orbitals_from_last_index(MyFragment,VOVOvirt,VOVOvirtTMP)
       if(pair) then
          ! Pair fragment
          call get_pair_fragment_energy(VOVOocc,VOVOvirtTMP,t2occ,t2virt,&
               & Fragment1, Fragment2, MyFragment)
       else
          ! Atomic fragment
          call get_atomic_fragment_energy(VOVOocc,VOVOvirtTMP,t2occ,t2virt,MyFragment)
       end if
       call tensor_free(VOVOvirtTMP)

    else ! use VOVOvirt as it is -- MOST COMMON CASE

       if(pair) then
          ! Pair fragment
          call get_pair_fragment_energy(VOVOocc,VOVOvirt,t2occ,t2virt,&
               & Fragment1, Fragment2, MyFragment)
       else
          ! Atomic fragment
          call get_atomic_fragment_energy(VOVOocc,VOVOvirt,t2occ,t2virt,MyFragment)
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
          !first order properties
       endif
#endif
       if(DECinfo%ccmodel == MODEL_MP2)then
          if(pair) then
             ! Pair fragment
             call pair_calculate_mp2gradient_driver(Fragment1,Fragment2,MyFragment,&
                  & t2occ,t2virt,VOOO,VOVV,VOVOocc,VOVOvirt,grad)
          else
             ! Atomic fragment
             call single_calculate_mp2gradient_driver(MyFragment,t2occ,t2virt,VOOO,&
                  & VOVV,VOVOocc,VOVOvirt,grad)
          end if
          call tensor_free(VOOO)
          call tensor_free(VOVV)
       end if
    end if
    !
    if(MyFragment%ccmodel /= MODEL_MP2.AND.MyFragment%ccmodel.NE.MODEL_RIMP2)then
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


  !> \brief Contract amplitudes, multipliers, and integrals to calculate atomic fragment Lagrangian energy.
  !> \author Kasper Kristensen
  !> \date August 2011
  subroutine get_atomic_fragment_energy(gocc,gvirt,t2occ,t2virt,MyFragment)

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
    integer :: noccEOS,nvirtEOS,noccAOS,nvirtAOS
    integer :: i,j,k,a,b,c
    real(realk) :: tcpu1, twall1, tcpu2,twall2, tcpu,twall
    real(realk) :: e1, e2, e3, e4,tmp,multaibj
    logical ::  something_wrong! ,doOccPart, doVirtPart
    real(realk) :: Eocc, lag_occ,Evirt,lag_virt
    real(realk),pointer :: occ_tmp(:),virt_tmp(:)
    real(realk) :: prefac_coul,prefac_k
    
    ! Lagrangian energy can be split into four contributions:
    ! The first two (e1 and e2) use occupied EOS orbitals and virtual AOS orbitals.
    ! The last two (e3 and e4) use virtual EOS orbitals and occupied AOS orbitals.
    !
    ! With the restrictions above the Lagrangian energy is given by:
    ! energy = e1 + e2 + e3 + e4
    ! e1 = sum_{ijab} t_{ij}^{ab} [ 2g_{aibj} - g_{biaj} ]
    ! e2 = 1/2 sum_{ijabc} mult_{ij}^{ab} [ t_{ij}^{cb} F_{ac} + t_{ij}^{ac} F_{bc} ]
    ! e3 = 1/2 sum_{ijab} mult_{ij}^{ab} g_{aibj}
    ! e4 = - 1/2 sum_{ijkab} mult_{ij}^{ab} [ t_{kj}^{ab} F_{ki} + t_{ik}^{ab} F_{kj} ]
    !
    ! IF -and only if- the fragment is the full molecule, the following simple relations hold:
    ! energy = e1 = e2 = -e3-e4
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
    nvirtEOS = MyFragment%nunoccEOS
    noccAOS  = MyFragment%noccAOS
    nvirtAOS = MyFragment%nunoccAOS
    Eocc     = 0E0_realk
    lag_occ  = 0E0_realk
    Evirt    = 0E0_realk
    lag_virt = 0E0_realk
    ! Just in case, zero individual orbital contributions for fragment
    MyFragment%OccContribs=0E0_realk
    MyFragment%VirtContribs=0E0_realk
    if(MyFragment%ccmodel==MODEL_RPA) then
      prefac_coul=1._realk
      prefac_k=0.0_realk
      if(Decinfo%SOS) prefac_k=0.5_realk
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

    call mem_TurnONThread_Memory()
    !$OMP PARALLEL DEFAULT(shared) PRIVATE(multaibj,tmp,e1,e2,j,b,i,a,c,virt_tmp)
    call init_threadmemvar()
    e1=0E0_realk
    e2=0E0_realk
    ! Contributions from each individual virtual orbital
    call mem_alloc(virt_tmp,nvirtAOS)
    virt_tmp = 0.0E0_realk       

    !$OMP DO SCHEDULE(dynamic,1)

    ! Calculate e1 and e2
    ! *******************
    do j=1,noccEOS
       do b=1,nvirtAOS
          do i=1,noccEOS
             do a=1,nvirtAOS


                ! Contribution 1
                ! --------------

                ! Energy contribution for orbitals (j,b,i,a)
                tmp = t2occ%elm4(a,i,b,j)*(prefac_coul*gocc%elm4(a,i,b,j) -prefac_k*gocc%elm4(b,i,a,j))

                ! Update total atomic fragment energy contribution 1
                e1 = e1 + tmp

                ! Update contribution from orbital a
                virt_tmp(a) = virt_tmp(a) + tmp
                ! Update contribution from orbital b (only if different from a to avoid double counting)
                if(a/=b) virt_tmp(b) = virt_tmp(b) + tmp


                ! Contribution 2
                ! --------------

                ! Skip contribution 2 for anything but MP2
                if(MyFragment%ccmodel==MODEL_MP2) then
                   ! Multiplier (multiplied by one half)
                   multaibj = prefac_coul*t2occ%elm4(a,i,b,j) - prefac_k*t2occ%elm4(b,i,a,j)

                   do c=1,nvirtAOS

                      ! Energy contribution for orbitals (j,b,i,a,c)
                      tmp = t2occ%elm4(c,i,b,j)*MyFragment%qqfock(c,a) + t2occ%elm4(a,i,c,j)*MyFragment%qqfock(c,b)
                      tmp = multaibj*tmp

                      ! Update total atomic fragment energy contribution 2
                      e2 = e2 + tmp

                      ! Update contribution from orbital a
                      virt_tmp(a) = virt_tmp(a) + tmp
                      ! Update contribution from orbital b (only if different from a to avoid double counting)
                      if(a/=b) virt_tmp(b) = virt_tmp(b) + tmp
                      ! Update contribution from orbital c (only if different from a and b)
                      if( (a/=c) .and. (b/=c) ) virt_tmp(c) = virt_tmp(c) + tmp

                   end do

                end if


             end do
          end do
       end do
    end do
    !$OMP END DO NOWAIT

    ! Total e1, e2, and individual virt atomic contributions are found by summing all thread contributions
    !$OMP CRITICAL
    Eocc = Eocc + e1
    lag_occ = lag_occ + e2

    ! Update total virtual contributions to fragment energy
    do a=1,nvirtAOS
       MyFragment%VirtContribs(a) =MyFragment%VirtContribs(a) + virt_tmp(a)
    end do
    !$OMP END CRITICAL

    call mem_dealloc(virt_tmp)
    call collect_thread_memory()
    !$OMP END PARALLEL
    call mem_TurnOffThread_Memory()


    ! Calculate e3 and e4
    ! *******************

    call mem_TurnONThread_Memory()
    !$OMP PARALLEL DEFAULT(shared) PRIVATE(multaibj,tmp,e3,e4,j,b,i,a,k,occ_tmp)
    call init_threadmemvar()
    e3=0E0_realk
    e4=0E0_realk

    ! Contributions from each individual occupied orbital
    call mem_alloc(occ_tmp,noccAOS)
    occ_tmp = 0.0E0_realk


    !$OMP DO SCHEDULE(dynamic,1)
    do j=1,noccAOS
       do b=1,nvirtEOS
          do i=1,noccAOS
             do a=1,nvirtEOS


                ! Contribution 3
                ! --------------

                ! Multiplier (multiplied by one half)
                multaibj = prefac_coul*t2virt%elm4(a,i,b,j) -prefac_k*t2virt%elm4(b,i,a,j)


                ! Energy contribution for orbitals (j,b,i,a)
                tmp = multaibj*gvirt%elm4(a,i,b,j)
                ! Update total atomic fragment energy contribution 3
                e3 = e3 + tmp

                ! Update contribution from orbital i
                occ_tmp(i) = occ_tmp(i) + tmp
                ! Update contribution from orbital j (only if different from i to avoid double counting)
                if(i/=j) occ_tmp(j) = occ_tmp(j) + tmp

                ! Contribution 4
                ! --------------

                ! Skip contribution 4 for anything but MP2
                if(MyFragment%ccmodel==MODEL_MP2) then

                   do k=1,noccAOS

                      tmp =  t2virt%elm4(a,k,b,j)*MyFragment%ppfock(k,i) &
                           & + t2virt%elm4(a,i,b,k)*MyFragment%ppfock(k,j)
                      tmp = -multaibj*tmp

                      ! Update total atomic fragment energy contribution 4
                      e4 = e4 + tmp

                      ! Update contribution from orbital i
                      occ_tmp(i) = occ_tmp(i) + tmp
                      ! Update contribution from orbital j (only if different from i to avoid double counting)
                      if(i/=j) occ_tmp(j) = occ_tmp(j) + tmp
                      ! Update contribution from orbital k (only if different from i and j)
                      if( (i/=k) .and. (j/=k) ) occ_tmp(k) = occ_tmp(k) + tmp

                   end do

                end if

             end do
          end do
       end do
    end do
    !$OMP END DO NOWAIT

    ! Total e3, e4, and individual occ atomic contributions are found by summing all thread contributions
    !$OMP CRITICAL
    Evirt = Evirt + e3
    lag_virt = lag_virt + e4

    ! Update total occupied contributions to fragment energy
    do i=1,noccAOS
       MyFragment%OccContribs(i) = MyFragment%OccContribs(i) + occ_tmp(i)
    end do
    !$OMP END CRITICAL

    call mem_dealloc(occ_tmp)
    call collect_thread_memory()
    !$OMP END PARALLEL
    call mem_TurnOffThread_Memory()


    ! Total atomic fragment energy
    ! ****************************
    if(MyFragment%ccmodel==MODEL_MP2) then
       ! Lagrangian energy only implemented for MP2 so it gets special treatment
       MyFragment%energies(FRAGMODEL_LAGMP2) = Eocc + lag_occ + Evirt + lag_virt
    end if
    ! Put occupied (Eocc) and virtual (Evirt) scheme energies into fragment energies array
    call put_fragment_energy_contribs_main(Eocc,Evirt,MyFragment)

    ! Set energies used by fragment optimization
    call get_occ_virt_lag_energies_fragopt(MyFragment)

    ! Print out contributions
    ! ***********************

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

    call LSTIMER('START',tcpu2,twall2,DECinfo%output)
    call LSTIMER('L.ENERGY CONTR',tcpu,twall,DECinfo%output)


  end subroutine get_atomic_fragment_energy



  !> \brief Wrapper for pair_driver where can attach existing full
  !> molecular singles amplitudes to the fragment structure and update new
  !> improved full molecular singles amplitudes by the calculated fragment singles amplitudes.
  !> \author Kasper Kristensen
  !> \date March 2012
  subroutine pair_driver_singles(natoms,nocc,nunocc,&
       & OccOrbitals,UnoccOrbitals,MyLsitem,MyMolecule,&
       & Fragment1,Fragment2,PairFragment,t1old,t1new)

    implicit none
    !> Number of atoms for full molecule
    integer, intent(in) :: natoms
    !> Number of occupied orbitals in full molecule
    integer,intent(in) :: nOcc
    !> Number of unoccupied orbitals in full molecule
    integer,intent(in) :: nUnocc
    !> Occupied orbitals for full molecule
    type(decorbital), intent(in) :: OccOrbitals(nocc)
    !> Unoccupied orbitals for full molecule
    type(decorbital), intent(in) :: UnoccOrbitals(nunocc)
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
    call update_full_t1_from_pair_frag(PairFragment,nocc,nunocc,&
         & natoms,dopair,OccOrbitals,UnoccOrbitals,t1new)
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
        & Fragment1, Fragment2, PairFragment)


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
     integer :: noccEOS,nvirtEOS,noccAOS,nvirtAOS
     integer :: i,j,k,a,b,c
     real(realk) :: tcpu, twall,pairdist
     real(realk) :: e1, e2, e3, e4,tmp,multaibj
     real(realk) :: tcpu1,tcpu2,twall1,twall2
     logical,pointer :: dopair_occ(:,:), dopair_virt(:,:)
     real(realk) :: Eocc, lag_occ,Evirt,lag_virt
     logical :: something_wrong, do_non_pdm
     real(realk) :: prefac_coul,prefac_k


     ! Pair interaction Lagrangian energy can be split into four contributions:
     ! The first two (e1 and e2) use occupied EOS orbitals and virtual AOS orbitals.
     ! The last two (e3 and e4) use virtual EOS orbitals and occupied AOS orbitals.
     !
     ! With the restrictions above the Lagrangian energy is given by:
     ! energy = e1 + e2 + e3 + e4
     ! e1 = sum_{ijab} t_{ij}^{ab} [ 2g_{aibj} - g_{biaj} ]
     ! e2 = 1/2 sum_{ijabc} mult_{ij}^{ab} [ t_{ij}^{cb} F_{ac} + t_{ij}^{ac} F_{bc} ]
     ! e3 = 1/2 sum_{ijab} mult_{ij}^{ab} g_{aibj}
     ! e4 = - 1/2 sum_{ijkab} mult_{ij}^{ab} [ t_{kj}^{ab} F_{ki} + t_{ik}^{ab} F_{kj} ]
     !
     ! Additional restriction to avoid double counting:
     ! e1 and e2: i and j must belong to two different atoms,
     !            e.g., i belongs to atom1 and j belongs to atom2 - or vice versa
     ! e1 and e2: a and b must belong to two different atoms,
     !            e.g., a belongs to atom1 and b belongs to atom2 - or vice versa
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
     noccEOS = PairFragment%noccEOS
     nvirtEOS = PairFragment%nunoccEOS
     noccAOS = PairFragment%noccAOS
     nvirtAOS = PairFragment%nunoccAOS
     Eocc=0E0_realk
     lag_occ=0E0_realk
     Evirt=0E0_realk
     lag_virt=0E0_realk
     ! Distance between fragments in Angstrom
     pairdist = bohr_to_angstrom*PairFragment%pairdist
     if(PairFragment%ccmodel==MODEL_RPA) then
        prefac_coul = 1._realk
        prefac_k=0.5_realk
     else
        prefac_coul =2._realk
        prefac_k = 1._realk
     endif

     ! Which "interaction pairs" to include for occ and unocc space (avoid double counting)
     call mem_alloc(dopair_occ,noccEOS,noccEOS)
     call mem_alloc(dopair_virt,nvirtEOS,nvirtEOS)
     call which_pairs_occ(Fragment1,Fragment2,PairFragment,dopair_occ)
     call which_pairs_unocc(Fragment1,Fragment2,PairFragment,dopair_virt)

     ! Sanity checks
     ! *************
     something_wrong=.false.
     IF(.NOT.DECinfo%OnlyVirtPart)THEN
        if(t2occ%dims(1) /= nvirtAOS) something_wrong=.true.
        if(t2occ%dims(2) /= noccEOS) something_wrong=.true.
        if(t2occ%dims(3) /= nvirtAOS) something_wrong=.true.
        if(t2occ%dims(4) /= noccEOS) something_wrong=.true.

        if(gocc%dims(1) /= nvirtAOS) something_wrong=.true.
        if(gocc%dims(2) /= noccEOS) something_wrong=.true.
        if(gocc%dims(3) /= nvirtAOS) something_wrong=.true.
        if(gocc%dims(4) /= noccEOS) something_wrong=.true.
     ENDIF

     IF(.NOT.DECinfo%OnlyOccPart)THEN
        if(t2virt%dims(1) /= nvirtEOS) something_wrong=.true.
        if(t2virt%dims(2) /= noccAOS) something_wrong=.true.
        if(t2virt%dims(3) /= nvirtEOS) something_wrong=.true.
        if(t2virt%dims(4) /= noccAOS) something_wrong=.true.

        if(gvirt%dims(1) /= nvirtEOS) something_wrong=.true.
        if(gvirt%dims(2) /= noccAOS) something_wrong=.true.
        if(gvirt%dims(3) /= nvirtEOS) something_wrong=.true.
        if(gvirt%dims(4) /= noccAOS) something_wrong=.true.
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
        call lsquit('get_pair_fragment_energy: &
           & Input dimensions do not match!',-1)
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

           ! Calculate e1 and e2
           ! *******************

           call mem_TurnONThread_Memory()
           !$OMP PARALLEL DEFAULT(shared) PRIVATE(multaibj,tmp,e1,e2,j,b,i,a,c)
           call init_threadmemvar()
           e1=0E0_realk
           e2=0E0_realk

           !$OMP DO SCHEDULE(dynamic,1)

           do j=1,noccEOS
              do b=1,nvirtAOS
                 do i=1,noccEOS

                    ! Only update for "interaction orbital pairs" - see which_pairs_occ
                    if( dopair_occ(i,j) ) then  !DoPair1and2

                       do a=1,nvirtAOS

                          ! Update pair interaction energy contribution 1
                          e1 = e1 + t2occ%elm4(a,i,b,j)*(prefac_coul*gocc%elm4(a,i,b,j) -prefac_k*gocc%elm4(b,i,a,j))


                          ! Skip contribution 2 for anything but MP2
                          if(pairfragment%ccmodel==MODEL_MP2) then

                             ! Multiplier (multiplied by one half)
                             multaibj = prefac_coul*t2occ%elm4(a,i,b,j) - prefac_k*t2occ%elm4(b,i,a,j)

                             tmp = 0E0_realk
                             do c=1,nvirtAOS
                                tmp = tmp + t2occ%elm4(c,i,b,j)*PairFragment%qqfock(c,a) &
                                   & + t2occ%elm4(a,i,c,j)*PairFragment%qqfock(c,b)
                             end do

                             ! Update pair interaction energy contribution 2
                             e2 = e2 + multaibj*tmp

                          end if


                       end do

                    end if

                 end do
              end do
           end do

           !$OMP END DO NOWAIT

           ! Total e1 and e2 contributions are found by summing all thread contributions
           !$OMP CRITICAL
           Eocc = Eocc + e1
           lag_occ = lag_occ + e2
           !$OMP END CRITICAL

           call collect_thread_memory()
           !$OMP END PARALLEL
           call mem_TurnOffThread_Memory()
        ELSE
           Eocc = 0.0E0_realk
           lag_occ = 0.0E0_realk
        ENDIF

        if(.not. DECinfo%onlyoccpart) then

           ! Calculate e3 and e4
           ! *******************

           call mem_TurnONThread_Memory()
           !$OMP PARALLEL DEFAULT(shared) PRIVATE(multaibj,tmp,e3,e4,j,b,i,a,k)
           call init_threadmemvar()
           e3=0e0_realk
           e4=0e0_realk

           !$OMP DO SCHEDULE(dynamic,1)

           do j=1,noccAOS
              do b=1,nvirtEOS
                 do i=1,noccAOS
                    do a=1,nvirtEOS

                       ! Only update for "interaction orbital pairs" - see which_pairs_unocc
                       if( dopair_virt(a,b) ) then !Dopair3and4

                          ! Multiplier (multiplied by one half)
                          multaibj = prefac_coul*t2virt%elm4(a,i,b,j) - prefac_k*t2virt%elm4(b,i,a,j)

                          ! Update total atomic fragment energy contribution 3
                          e3 = e3 + multaibj*gvirt%elm4(a,i,b,j)


                          ! Skip contribution 4 for anything but MP2
                          if(pairfragment%ccmodel==MODEL_MP2) then

                             tmp=0E0_realk
                             do k=1,noccAOS
                                tmp = tmp + t2virt%elm4(a,k,b,j)*PairFragment%ppfock(k,i) &
                                   & + t2virt%elm4(a,i,b,k)*PairFragment%ppfock(k,j)
                             end do

                             ! Update pair interaction energy contribution 4
                             e4 = e4 - multaibj*tmp

                          end if

                       end if

                    end do
                 end do
              end do
           end do

           !$OMP END DO NOWAIT

           ! Total e3 and e4 contributions are found by summing all thread contributions
           !$OMP CRITICAL
           Evirt = Evirt + e3
           lag_virt = lag_virt + e4
           !$OMP END CRITICAL

           call collect_thread_memory()
           !$OMP END PARALLEL
           call mem_TurnOffThread_Memory()
        ELSE
           Evirt = 0.0E0_realk
           lag_virt = 0.0E0_realk      
        ENDIF

     else
        call lsquit("ERROR(get_pair_fragment_energy) PDM version not yetimplemented",-1)
     endif

     ! Total pair interaction energy
     ! *****************************
     if(PairFragment%ccmodel==MODEL_MP2) then
        ! Lagrangian energy only implemented for MP2 so it gets special treatment
        PairFragment%energies(FRAGMODEL_LAGMP2) = Eocc + lag_occ + Evirt + lag_virt
     end if
     ! Put occupied (Eocc) and virtual (Evirt) scheme energies into fragment energies array
     call put_fragment_energy_contribs_main(Eocc,Evirt,PairFragment)

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

  !> \brief Full DEC calculation where exact single and pair energies are calculated using Lagrangian approach.
  !> Only implemented for MP2.
  !> \author Kasper Kristensen
  !> \date April 2011
  subroutine Full_DECMP2_calculation(MyMolecule,mylsitem,Ecorr)

    implicit none
    !> Molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> LS Dalton info
    type(lsitem), intent(inout) :: mylsitem
    !> Total correlation energy
    real(realk),intent(inout) :: Ecorr
    logical,dimension(MyMolecule%natoms) :: orbitals_assigned
    real(realk),pointer :: Cocc(:,:), Cvirt(:,:)
    type(tensor) :: t2_t, g_t
    type(array4) :: t2, g
    real(realk) :: energy_matrix(MyMolecule%natoms,MyMolecule%natoms), multaibj, multbiaj
    integer :: nthreads, idx, nbatchINT, intstep, ncore,offset
    integer :: i,j,k,a,b,c,atomI,atomJ,atomA,atomB
    real(realk) :: intMEM, solMEM,OO,VV,AA,BB,mem_required
    real(realk) :: singleenergy, pairenergy, tmp, tmp2, InteractionEcorr
    real(realk),dimension(MyMolecule%natoms,MyMolecule%natoms) :: e1,e2,e3,e4,e1_tmp,e2_tmp,e3_tmp,e4_tmp
    integer, dimension(4) :: dims
    real(realk), pointer :: gval(:,:,:),t2val(:,:,:),ppfock(:,:)
#ifdef VAR_OMP
    integer, external :: OMP_GET_MAX_THREADS
#endif
    type(decorbital), pointer :: OccOrbitals(:)
    type(decorbital), pointer :: UnoccOrbitals(:)
    integer :: nocc,nunocc,nbasis,natoms,lupri

    write(DECinfo%output,*) 'Using DEC-MP2 debug routine for full molecular system...'

    ! Only for MP2
    if(DECinfo%ccModel/=MODEL_MP2) then
       call lsquit('Full_DECMP2_calculation: Only implemented for MP2!', DECinfo%output)
    end if


    ncore = MyMolecule%ncore
    nOcc = MyMolecule%nocc
    nUnocc = MyMolecule%nunocc
    nBasis = MyMolecule%nbasis
    nAtoms = MyMolecule%natoms

    ! -- Analyze basis and create orbitals 
    call mem_alloc(OccOrbitals,nOcc)
    call mem_alloc(UnoccOrbitals,nUnocc)
    call GenerateOrbitals_driver(MyMolecule,mylsitem,nocc,nunocc,natoms, &
         & OccOrbitals, UnoccOrbitals)


    if(DECinfo%frozencore) then
       ! Frozen core: Only valence orbitals
       nocc = MyMolecule%nval
    else
       nocc = MyMolecule%nocc
    end if


    ! Initialize stuff
    ! ****************
    energy_matrix(:,:) = 0E0_realk
    dims = [nunocc, nocc, nunocc, nocc]
    call mem_alloc(ppfock,nocc,nocc)
    if(DECinfo%frozencore) then
       ! Only copy valence orbitals into array2 structure
       call mem_alloc(Cocc,nbasis,nocc)
       do i=1,nocc
          Cocc(:,i) = MyMolecule%Co(:,i+Ncore)
       end do

       ! Fock valence
       do j=1,nocc
          do i=1,nocc
             ppfock(i,j) = MyMolecule%ppfock(i+Ncore,j+Ncore)
          end do
       end do
       offset = ncore
    else
       ! No frozen core, simply copy elements for all occupied orbitals
       call mem_alloc(Cocc,nbasis,nocc)
       Cocc=MyMolecule%Co
       ppfock = MyMolecule%ppfock
       offset=0
    end if
    call mem_alloc(Cvirt,nbasis,nunocc)
    Cvirt = MyMolecule%Cv
    e1=0E0_realk
    e2=0E0_realk
    e3=0E0_realk
    e4=0E0_realk

    ! Number of OMP threads
#ifdef VAR_OMP
    ! LSDALTON compiled with OMP
    nthreads=OMP_GET_MAX_THREADS()
#else
    ! No OMP, set number of threads to one
    nthreads=1
#endif

    ! Orbital assignment
    orbitals_assigned=.false.
    do i=1,nocc
       idx = OccOrbitals(i)%centralatom
       orbitals_assigned(idx) = .true.
    end do
    do i=1,nunocc
       idx = UnoccOrbitals(i)%centralatom
       orbitals_assigned(idx) = .true.
    end do


    ! Estimate memory consumption
    ! ***************************
    nbatchINT = max_batch_dimension(mylsitem,nbasis)
    OO=nocc ! Number of occupied orbitals (as real)
    VV=nunocc ! Number of virtual orbitals (as real)
    BB=nbatchINT ! Maximum batch dimension (as real)
    AA=nbasis ! Number of atomic orbitals (as real)

    call estimate_memory_for_mp2_energy(nthreads,OO,VV,AA,BB,intMEM,intStep,solMEM)
    mem_required = max(intMEM,solMEM)
    mem_required = mem_required + DECinfo%fullmolecule_memory
    if(mem_required > DECinfo%memory) then
       DECinfo%array4OnFile =.true.
    else
       DECinfo%array4OnFile =.false.
    end if
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) '========================================================'
    write(DECinfo%output,*) '           FULL DEC: MEMORY ESTIMATE SUMMARY'
    write(DECinfo%output,*) '--------------------------------------------------------'
    write(DECinfo%output,'(1X,a,i9)') 'ME: Number of OMP threads         = ', &
         & nthreads
    write(DECinfo%output,'(1X,a,i9)') 'ME: Number of occupied orbitals   = ', &
         & nocc
    write(DECinfo%output,'(1X,a,i9)') 'ME: Number of virtual orbitals    = ', &
         & nunocc
    write(DECinfo%output,'(1X,a,i9)') 'ME: Number of atomic orbitals     = ', &
         & nbasis
    write(DECinfo%output,'(1X,a,i9)') 'ME: Maximum number of batches     = ', nbatchINT
    write(DECinfo%output,'(1X,a,g16.6,a,i3)') &
         & 'ME: Integral memory required (GB) = ', intMEM, ' in step ', intStep
    write(DECinfo%output,'(1X,a,g16.6)') 'ME: Memory in use, full mol  (GB) = ', &
         & DECinfo%fullmolecule_memory
    write(DECinfo%output,'(1X,a,g16.6)') 'ME: Solver memory required   (GB) = ', solMEM
    write(DECinfo%output,'(1X,a,g16.6)') 'ME: Maximum memory required  (GB) = ', mem_required
    write(DECinfo%output,'(1X,a,g16.6)') 'ME: Maximum memory available (GB) = ', DECinfo%memory
    write(DECinfo%output,*)
    if(intMEM > solMEM) then
       if(DECinfo%array4OnFile) then
          write(DECinfo%output,*) 'ME: INTEGRAL ROUTINE DEFINES MEMORY CONSUMPTION &
               & --> STORE ON FILE'
       else
          write(DECinfo%output,*) 'ME: INTEGRAL ROUTINE DEFINES MEMORY CONSUMPTION &
               & --> STORE IN MEMORY'
       end if

    else
       if(DECinfo%array4OnFile) then
          write(DECinfo%output,*) 'ME: SOLVER ROUTINE DEFINES MEMORY CONSUMPTION &
               & --> STORE ON FILE'
       else
          write(DECinfo%output,*) 'ME: SOLVER ROUTINE DEFINES MEMORY CONSUMPTION &
               & --> STORE IN MEMORY'
       end if
    end if
    write(DECinfo%output,*) '--------------------------------------------------------'
    write(DECinfo%output,*)
    write(DECinfo%output,*)


    ! Get (C K | D L) integrals stored in the order (C,K,D,L)
    ! *******************************************************
    !call get_VOVO_integrals(mylsitem,nbasis,nocc,nunocc,Cvirt,Cocc,g)
    call mem_dealloc(Cocc)
    call mem_dealloc(Cvirt)


    ! Get t2 amplitudes
    ! *****************
    call mp2_solver(MyMolecule,mylsitem,g_t,t2_t,.false.)


    !FIXME:THIS IS A DIRTY WORKAROUND NOT TO CHANGE TOO MUCH IN THE SUBROUTINE
    !-------------------------------------------------------------------------
    t2 = array4_init(t2%dims)                                                !
    g  = array4_init(g%dims)                                                 !
                                                                             !
    call tensor_convert(t2_t,t2%val)                                         !
    call tensor_convert(g_t,g%val)                                           !
                                                                             !
    call tensor_free(t2_t)                                                   !
    call tensor_free(g_t)                                                    !
    !-------------------------------------------------------------------------


    ! STATUS: Now integrals (g) and amplitudes (t) have been determined
    ! for the full molecular system. The individual fragment contributions - solved in the
    ! total orbital space - can now be determined by extracting EOS indices from the full indicies.



    ! Energy contributions (see fragment_energy_and_prop subroutine)
    ! ***************************************************************************

    if(DECinfo%array4OnFile) then ! Integrals and amplitudes on file
       call array4_open_file(g)
       call array4_open_file(t2)
    end if
    call mem_alloc(gval,nunocc,nocc,nunocc)
    call mem_alloc(t2val,nunocc,nocc,nunocc)


    ! Calculate e1 and e2
    ! -------------------
    do j=1,nocc
       ! Central atom for orbital j
       atomJ = OccOrbitals(j+offset)%CentralAtom

       if(DECinfo%array4OnFile) then ! Read values from file
          ! g
          call array4_read_file_type1(g,j,&
               & gval(1:nunocc,1:nocc,1:nunocc),nunocc,nocc,nunocc)
          ! t2
          call array4_read_file_type1(t2,j,&
               & t2val(1:nunocc,1:nocc,1:nunocc),nunocc,nocc,nunocc)

       else ! Values exist in memory

          gval(1:nunocc,1:nocc,1:nunocc) = g%val(1:nunocc,1:nocc,1:nunocc,j)
          t2val(1:nunocc,1:nocc,1:nunocc) = t2%val(1:nunocc,1:nocc,1:nunocc,j)

       end if

       call mem_TurnONThread_Memory()
       !$OMP PARALLEL DEFAULT(shared) PRIVATE(b,atomB,i,atomI,a,atomA,multaibj,multbiaj,&
       !$OMP tmp,tmp2,e1_tmp,e2_tmp,e3_tmp,e4_tmp,k,c)
       call init_threadmemvar()
       e1_tmp=0E0_realk
       e2_tmp=0E0_realk
       e3_tmp=0E0_realk
       e4_tmp=0E0_realk

       !$OMP DO SCHEDULE(dynamic,1)

       do b=1,nunocc
          atomB = UnoccOrbitals(b)%CentralAtom
          do i=1,nocc
             atomI = OccOrbitals(i+offset)%CentralAtom
             do a=1,nunocc
                atomA = UnoccOrbitals(a)%CentralAtom

                ! Multiplier (multiplied by one half)
                multaibj = 2.0_realk*t2val(a,i,b) - t2val(b,i,a)
                multbiaj = 2.0_realk*t2val(b,i,a) - t2val(a,i,b)


                ! Contribution 1
                ! ''''''''''''''
                tmp = t2val(a,i,b)*(2.0_realk*gval(a,i,b) - gval(b,i,a))
                ! Update total atomic fragment energy contribution 1
                e1_tmp(atomI,atomJ) = e1_tmp(atomI,atomJ) + tmp


                ! Contribution 2
                ! ''''''''''''''
                do c=1,nunocc
                   tmp = t2val(c,i,b)*MyMolecule%qqfock(a,c) + t2val(a,i,c)*MyMolecule%qqfock(b,c)
                   tmp = multaibj*tmp
                   ! Update total atomic fragment energy contribution 2
                   e2_tmp(atomI,atomJ) = e2_tmp(atomI,atomJ) + tmp
                end do


                ! Contribution 3
                ! ''''''''''''''
                tmp = multaibj*gval(a,i,b)
                ! Update total atomic fragment energy contribution 3
                e3_tmp(atomA,atomB) = e3_tmp(atomA,atomB) + tmp


                ! Contribution 4
                ! --------------

                ! Note: Here we have renamed i <--> j in the second term compared to
                ! the second term of e4 in get_atomic_fragment_energy to be able
                ! to store the arrays on file in a simple manner.

                do k=1,nocc

                   tmp =  t2val(a,k,b)*ppfock(k,i)
                   tmp = multaibj*tmp

                   tmp2 = t2val(b,k,a)*ppfock(k,i)
                   tmp2 = multbiaj*tmp2

                   ! Update total atomic fragment energy contribution 4
                   e4_tmp(atomA,atomB) = e4_tmp(atomA,atomB) - tmp - tmp2

                end do


             end do
          end do
       end do

       !$OMP END DO NOWAIT

       ! Sum up contributions from each thread
       !$OMP CRITICAL
       e1=e1+e1_tmp
       e2=e2+e2_tmp
       e3=e3+e3_tmp
       e4=e4+e4_tmp
       !$OMP END CRITICAL

       call collect_thread_memory()
       !$OMP END PARALLEL
       call mem_TurnOffThread_Memory()

    end do


    if(DECinfo%array4OnFile) then ! Integrals and amplitudes on file
       call array4_close_file(g,'DELETE')
       call array4_close_file(t2,'DELETE')
    end if

    call array4_free(t2)
    call array4_free(g)
    call mem_dealloc(t2val)
    call mem_dealloc(gval)


    ! Total fragment correlation density matrix
    Ecorr =0.0_realk
    do atomI=1,natoms
       do atomJ=1,natoms
          energy_matrix(atomI,atomJ) = e1(atomI,atomJ) &
               & + e2(atomI,atomJ) + e3(atomI,atomJ) + e4(atomI,atomJ)
          ! Calculate correlation energy using occ scheme
          ! (of course we get the same for the other schemes).
          Ecorr = Ecorr + e1(atomI,atomJ)
       end do
    end do
    
    ! Only consider pairs IJ where J>I; thus, move contributions
    do atomI=1,natoms
       do atomJ=atomI+1,natoms
          energy_matrix(atomI,atomJ) =energy_matrix(atomI,atomJ) + energy_matrix(atomJ,atomI)
          e1(atomI,atomJ) =e1(atomI,atomJ) + e1(atomJ,atomI)
          e2(atomI,atomJ) =e2(atomI,atomJ) + e2(atomJ,atomI)
          e3(atomI,atomJ) =e3(atomI,atomJ) + e3(atomJ,atomI)
          e4(atomI,atomJ) =e4(atomI,atomJ) + e4(atomJ,atomI)
       end do
    end do


    ! Print stuff
    if(.not. DECinfo%onlyvirtpart) then
       call print_atomic_fragment_energies(natoms,e1,orbitals_assigned,&
            & 'MP2 occupied single energies','AF_MP2_OCC')
    endif
    if(.not.DECinfo%onlyoccpart) then
       call print_atomic_fragment_energies(natoms,e3,orbitals_assigned,&
            & 'MP2 virtual single energies','AF_MP2_VIR')
    endif
    if(.not.(DECinfo%onlyoccpart.or. DECinfo%onlyvirtpart)) then
       call print_atomic_fragment_energies(natoms,energy_matrix,orbitals_assigned,&
            & 'MP2 Lagrangian single energies','AF_MP2_LAG')
    end if
    if(.not. DECinfo%onlyvirtpart) then
       call print_pair_fragment_energies(natoms,e1,orbitals_assigned,&
            & MyMolecule%DistanceTable, 'MP2 occupied pair energies','PF_MP2_OCC')
    endif
    if(.not.DECinfo%onlyoccpart) then
       call print_pair_fragment_energies(natoms,e3,orbitals_assigned,&
            & MyMolecule%DistanceTable, 'MP2 virtual pair energies','PF_MP2_VIR')
    endif
    if(.not.(DECinfo%onlyoccpart.or. DECinfo%onlyvirtpart)) then
       call print_pair_fragment_energies(natoms,energy_matrix,orbitals_assigned,&
            & MyMolecule%DistanceTable, 'MP2 Lagrangian pair energies','PF_MP2_LAG')
    end if
    call mem_dealloc(ppfock)

    do i=1,nOcc
       call orbital_free(OccOrbitals(i))
    end do
    do i=1,nUnocc
       call orbital_free(UnoccOrbitals(i))
    end do
    call mem_dealloc(OccOrbitals)
    call mem_dealloc(UnoccOrbitals)

  end subroutine Full_DECMP2_calculation




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


  !> Routine that optimizes an atomic fragment by (1) expanding to include neighbour atoms and
  !> checking that energy change is smaller than FOT, (2) for converged fragment we 
  !> reduce fragment again by removing individual orbitals - either (2a) removing
  !> local orbitals according to their energy contributions OR (2b) generating
  !> fragment-adapted orbitals and removing those according to their eigenvalues.
  !> In (2) it is ensured that the energy error introduced by reducing the fragment (compared to 
  !> the fragment obtained in (1)) is below the FOT.
  !> \date February 2013
  !> \author Ida-Marie Hoeyvik & Kasper Kristensen & Thomas Kjaergaard
  subroutine optimize_atomic_fragment(MyAtom,AtomicFragment,nAtoms, &
        &OccOrbitals,nOcc,UnoccOrbitals,nUnocc,&
        &MyMolecule,mylsitem,freebasisinfo,t1full)
     implicit none
     !> Number of occupied orbitals in molecule
     integer, intent(in) :: nOcc
     !> Number of unoccupied orbitals in molecule
     integer, intent(in) :: nunocc
     !> Number of atoms in molecule
     integer, intent(in) :: natoms
     !> Central atom in molecule
     integer, intent(inout) :: MyAtom !> Wangy hack intent inout
     !> Atomic fragment to be optimized
     type(decfrag),intent(inout)        :: AtomicFragment
     !> All occupied orbitals
     type(decorbital), dimension(nOcc), intent(in)      :: OccOrbitals
     !> All unoccupied orbitals
     type(decorbital), dimension(nUnocc), intent(in)    :: UnoccOrbitals
     !> Full molecule information
     type(fullmolecule), intent(inout) :: MyMolecule
     !> Integral information
     type(lsitem), intent(inout)       :: mylsitem
     !> Delete fragment basis information ("expensive box in decfrag type") at exit?
     logical,intent(in) :: freebasisinfo
     !> t1 amplitudes for full molecule to be updated (only used when DECinfo%SinglesPolari is set)
     type(array2),intent(inout),optional :: t1full
     real(realk)                    :: LagEnergyDiff, OccEnergyDiff,VirtEnergyDiff
     real(realk)                    :: LagEnergyOld, OccEnergyOld, VirtEnergyOld, FOT, FOT2
     real(realk)                    :: init_Occradius, init_Virtradius
     logical, dimension(natoms)     :: Occ_atoms,Virt_atoms,OccOld,VirtOld
     real(realk),dimension(natoms)  :: DistMyAtom,SortedDistMyAtom
     integer,dimension(natoms)      :: DistTrackMyAtom, nocc_per_atom,nunocc_per_atom
     integer      :: iter,i,idx,StepsizeLoop3,StepsizeLoop2,NF,IT,StepsizeLoop4
     integer      :: max_iter_red,nAtomsWithOccOrb,nAtomsWithVirtOrb,ncore
     logical :: expansion_converged,ExpandFragmentConverged,DistanceRemoval
     logical :: OccUnchanged,VirtUnchanged,ExpandVirt,ExpandOcc
     logical :: BinarySearch,SeperateExpansion,OrbDistanceSpec,TestOcc
     type(tensor) :: t2,g
     real(realk)  :: FmaxOcc(nocc),FmaxVirt(nunocc)
     real(realk),pointer :: OccContribs(:),VirtContribs(:)    
     real(realk),pointer :: times_fragopt(:)
     real(realk),pointer :: SortedDistanceTableOrbAtomOcc(:)
     real(realk),pointer :: SortedDistanceTableOrbAtomVirt(:)
     integer,pointer :: OrbOccDistTrackMyAtom(:),OrbOccFockTrackMyAtom(:)
     integer,pointer :: OrbVirtDistTrackMyAtom(:),OrbVirtFockTrackMyAtom(:)
     logical,pointer :: OccAOS(:),VirtAOS(:),OldOccAOS(:),OldVirtAOS(:)
     logical :: BruteForce,FockMatrixOrdering

!!$     !! HACK for testing purposes, F12 code, Do not remove 
!!$     !! ****
!!$     !! All virtual, change occupied
!!$
!!$     MyAtom = 1
!!$
!!$     do i=1,natoms
!!$        Occ_Atoms(i) = .False.
!!$        Virt_Atoms(i) = .True.
!!$     enddo
!!$
!!$     do i=1,natoms
!!$        Occ_Atoms(1:i) = .True.
!!$        print *, "-------------------------------------------------"
!!$        print *, "     All virtual", "Number of Occupied Orbitals", i
!!$        print *, "-------------------------------------------------"
!!$
!!$        call get_fragment_and_Energy(MyAtom,natoms,Occ_Atoms,Virt_Atoms,&
!!$             & MyMolecule,MyLsitem,nocc,nunocc,OccOrbitals,UnoccOrbitals,&
!!$             & AtomicFragment)
!!$        call atomic_fragment_free(AtomicFragment)
!!$     end do
!!$
!!$     !! All occupied, change virtual
!!$     do i=1,natoms
!!$        Occ_Atoms(i) = .True.
!!$        Virt_Atoms(i) = .False.
!!$     enddo
!!$
!!$     do i=1,natoms
!!$        Virt_Atoms(1:i) = .True.
!!$        print *, "-------------------------------------------------"
!!$        print *, "     All occ", "Number of Virtual Orbitals", i
!!$        print *, "-------------------------------------------------"
!!$
!!$        call get_fragment_and_Energy(MyAtom,natoms,Occ_Atoms,Virt_Atoms,&
!!$             & MyMolecule,MyLsitem,nocc,nunocc,OccOrbitals,UnoccOrbitals,&
!!$             & AtomicFragment)
!!$        call atomic_fragment_free(AtomicFragment)
!!$     end do
!!$
!!$     do i=1,natoms
!!$        Occ_Atoms(i) = .False.
!!$        Virt_Atoms(i) = .False.
!!$     enddo
!!$
!!$     do i=1,natoms
!!$        Virt_Atoms(1:i) = .True.
!!$        Occ_Atoms(1:i) = .True.
!!$        print *, "-------------------------------------------------"
!!$        print *, "     Occ and Virt", i
!!$        print *, "-------------------------------------------------"
!!$        call get_fragment_and_Energy(MyAtom,natoms,Occ_Atoms,Virt_Atoms,&
!!$             & MyMolecule,MyLsitem,nocc,nunocc,OccOrbitals,UnoccOrbitals,&
!!$             & AtomicFragment)
!!$        call atomic_fragment_free(AtomicFragment)
!!$     end do
!!$
!!$     stop 'KK/Wangy HACK'

     if (.not.DECinfo%no_orb_based_fragopt) then
       call optimize_atomic_fragment_clean(MyAtom,AtomicFragment,nAtoms, &
            & OccOrbitals,nOcc,UnoccOrbitals,nUnocc,&
            & MyMolecule,mylsitem,freebasisinfo)
       return
     end if

     !Number of core occupied orbitals in molecule
     ncore = MyMolecule%ncore
     times_fragopt => null()
     call dec_fragment_time_init(times_fragopt)
     
     write(DECinfo%output,'(a)')    ' FOP'
     write(DECinfo%output,'(a)')    ' FOP ==============================================='
     write(DECinfo%output,'(a,i4)') ' FOP  Site fragment generator for fragment,',MyAtom
     write(DECinfo%output,'(a)')    ' FOP ==============================================='
     write(DECinfo%output,'(a)')    ' FOP'

     BinarySearch = .FALSE.      !Use Binary Search in Reduction Scheme
     SeperateExpansion = .FALSE. !Expansion for Occ and Virt is done seperately 
     OrbDistanceSpec = .FALSE.   !Orbital Specific expansion
     DistanceRemoval = .FALSE.   !Use Distance to remove orbitals when no Energy Contribs 
     TestOcc = .FALSE.           !converge Occ first 
     BruteForce = .FALSE.        !BruteForce 
     FockMatrixOrdering = .FALSE.!Fock BruteForce 
     write(DECinfo%output,'(a)') ' FOP  Fragment optimization scheme '
     IF(DECinfo%Frag_Exp_Scheme.EQ.1)THEN
        write(DECinfo%output,'(a)') ' FOP  Standard Fragment optimization scheme is used'
     ELSEIF(DECinfo%Frag_Exp_Scheme.EQ.2)THEN
        write(DECinfo%output,'(a)') ' FOP Binary search is used in Fragment reduction scheme '        
        BinarySearch = .TRUE.
     ELSEIF(DECinfo%Frag_Exp_Scheme.EQ.3)THEN
        write(DECinfo%output,'(a)') ' FOP Binary search is used in Fragment reduction scheme '        
        write(DECinfo%output,'(a)') ' FOP Expansion for Occ and Virt is done seperately'        
        BinarySearch = .TRUE.
        SeperateExpansion = .TRUE.        
     ELSEIF(DECinfo%Frag_Exp_Scheme.EQ.4)THEN
        write(DECinfo%output,'(a)') ' FOP Binary search is used in Fragment reduction scheme '        
        write(DECinfo%output,'(a)') ' FOP Orbital distance specific expansion is used in Fragment expansion scheme '        
        BinarySearch = .TRUE.
        OrbDistanceSpec = .TRUE. 
     ELSEIF(DECinfo%Frag_Exp_Scheme.EQ.5)THEN
        write(DECinfo%output,'(a)') ' FOP Binary search is used in Fragment reduction scheme '
        write(DECinfo%output,'(a)') ' FOP Orbital distance specific expansion is used in Fragment expansion scheme '        
        IF(DECInfo%OnlyOccPart)THEN
           write(DECinfo%output,'(a)') ' Energy Contribution is used to remove Virtual Orbitals in Reduction Scheme'
           write(DECinfo%output,'(a)') ' Distance is used to remove Occupied Orbitals in  Reduction Scheme'        
           DistanceRemoval = .TRUE. 
        ELSEIF(DECInfo%OnlyVirtPart)THEN
           write(DECinfo%output,'(a)') ' Energy Contribution is used to remove Occupied Orbitals in Reduction Scheme'
           write(DECinfo%output,'(a)') ' Distance is used to remove Virtual Orbitals in  Reduction Scheme'        
           DistanceRemoval = .TRUE. 
        ELSE
           write(DECinfo%output,'(a)') ' This schem is identical to Scheme = 4'
           write(DECinfo%output,'(a)') ' unless ONLYOCCPART or ONLYVIRTPART is specified which you have not'
        ENDIF
        BinarySearch = .TRUE.
        OrbDistanceSpec = .TRUE. 
     ELSEIF(DECinfo%Frag_Exp_Scheme.EQ.6)THEN
        write(DECinfo%output,'(a)') ' FOP Binary search is used in Fragment reduction scheme '
        write(DECinfo%output,'(a)') ' FOP Orbital distance specific expansion is used in Fragment expansion scheme '        
        IF(DECInfo%OnlyOccPart)THEN
           write(DECinfo%output,'(a)') ' Energy Contribution is used to remove Virtual Orbitals in Reduction Scheme'
           write(DECinfo%output,'(a)') ' Distance is used to remove Occupied Orbitals in  Reduction Scheme'        
           DistanceRemoval = .TRUE. 
        ELSEIF(DECInfo%OnlyVirtPart)THEN
           write(DECinfo%output,'(a)') ' Energy Contribution is used to remove Occupied Orbitals in Reduction Scheme'
           write(DECinfo%output,'(a)') ' Distance is used to remove Virtual Orbitals in  Reduction Scheme'        
           DistanceRemoval = .TRUE. 
        ELSE
           write(DECinfo%output,'(a)') ' This schem is identical to Scheme = 4'
           write(DECinfo%output,'(a)') ' unless ONLYOCCPART or ONLYVIRTPART is specified which you have not'
        ENDIF
        write(DECinfo%output,'(a)') ' This scheme is identical to Scheme = 5 but first converges the Occ - then Virt'
        BinarySearch = .TRUE.
        OrbDistanceSpec = .TRUE. 
        TestOcc = .TRUE.
     ELSEIF(DECinfo%Frag_Exp_Scheme.EQ.7)THEN
        write(DECinfo%output,'(a)') ' FOP Orbital distance specific expansion is used in Fragment expansion scheme '
        IF(DECInfo%OnlyOccPart)THEN
           write(DECinfo%output,'(a)') ' FOP Brute Force Energy Contribution is used to Remove Occupied Orbitals'
        ELSEIF(DECInfo%OnlyVirtPart)THEN
           write(DECinfo%output,'(a)') ' FOP Brute Force Energy Contribution is used to Remove Virtual Orbitals'
        ELSE
!           write(DECinfo%output,'(a)') ' This schem is identical to Scheme = 4'
!           write(DECinfo%output,'(a)') ' unless ONLYOCCPART or ONLYVIRTPART is specified which you have not'
        ENDIF
        write(DECinfo%output,'(a)') ' FOP Binary search is used in Fragment reduction scheme '
        BinarySearch = .TRUE.
        OrbDistanceSpec = .TRUE. 
        BruteForce = .TRUE.
     ELSEIF(DECinfo%Frag_Exp_Scheme.EQ.8)THEN
        write(DECinfo%output,'(a)') ' FOP Orbital Fock Matrix Ordering specific expansion is used in Fragment expansion scheme '
        IF(DECInfo%OnlyOccPart)THEN
           write(DECinfo%output,'(a)') ' FOP Fock Matrix Ordering is used to Remove Occupied Orbitals'
        ELSEIF(DECInfo%OnlyVirtPart)THEN
           write(DECinfo%output,'(a)') ' FOP Fock Matrix Ordering is used to Remove Virtual Orbitals'
        ELSE
 !          write(DECinfo%output,'(a)') ' This schem is identical to Scheme = 4'
 !          write(DECinfo%output,'(a)') ' unless ONLYOCCPART or ONLYVIRTPART is specified which you have not'
        ENDIF
        write(DECinfo%output,'(a)') ' FOP Binary search is used in Fragment reduction scheme '
        BinarySearch = .TRUE.
        OrbDistanceSpec = .TRUE. 
        FockMatrixOrdering = .TRUE.
        call mem_alloc(OrbOccFockTrackMyAtom,nocc)
        call mem_alloc(OrbVirtFockTrackMyAtom,nunocc)
        call GetSortedFockMaxElements(FmaxOcc,FmaxVirt,nocc,nunocc,MyMolecule,ncore,OccOrbitals,&
             & UnOccOrbitals,MyAtom,OrbOccFockTrackMyAtom,OrbVirtFockTrackMyAtom)
     ELSE
        call lsquit('Frag_Opt_Scheme unknown',-1)
     ENDIF

     ! Sanity check for singles polarization
     if(DECinfo%SinglesPolari) then
        if(.not. present(t1full)) then
           call lsquit('optimize_atomic_fragment: Full singles polarization is requrested &
              & but t1 argument is not present!',DECinfo%output)
        end if
     end if


     ! ======================================================================
     !                    Initialization of various things...
     ! ======================================================================

     iter=0
     LagEnergyDiff=0.0_realk
     OccEnergyDiff=0.0_realk
     VirtEnergyDiff=0.0_realk
     expansion_converged=.false.
     max_iter_red=15   ! set to 100 for binary search where more steps might be needed
     if (BinarySearch) max_iter_red=100
     FOT = DECinfo%FOT
     DistMyAtom= mymolecule%DistanceTable(:,MyAtom)   ! distance vector for central atom
     ! Sort atoms according to distance from central atom
     call GetSortedList(SortedDistMyAtom,DistTrackMyAtom,mymolecule%DistanceTable,natoms,natoms,MyAtom)
     nocc_per_atom=get_number_of_orbitals_per_atom(OccOrbitals,nocc,natoms,.true.)
     nunocc_per_atom=get_number_of_orbitals_per_atom(UnoccOrbitals,nunocc,natoms,.true.)
     ! Only do fragment optimization if there are orbitals assigned to central atom.
     if( (nocc_per_atom(MyAtom) == 0) .and. (nunocc_per_atom(MyAtom) == 0) ) then
        write(DECinfo%output,*) 'FOP Skipping optimization of fragment ', MyAtom
        AtomicFragment%LagFOP=0E0_realk
        AtomicFragment%EoccFOP=0E0_realk
        AtomicFragment%EvirtFOP=0E0_realk
        AtomicFragment%energies=0E0_realk
        AtomicFragment%ccmodel = DECinfo%ccmodel
        call dec_fragment_time_get(times_fragopt)
        call dec_time_evaluate_efficiency_frag(AtomicFragment,times_fragopt,&
           &AtomicFragment%ccmodel,'Fragment optmization')
        return
     end if

     ! Debug case: Full molecule included in fragment
     ! --> we then skip energy calculation here and just init fragment
     if(DECinfo%InclFullMolecule .or. DECinfo%simulate_full) then
        call fragopt_include_fullmolecule(MyAtom,AtomicFragment, &
           &OccOrbitals,nOcc,UnoccOrbitals,nUnocc, &
           &MyMolecule,mylsitem,freebasisinfo,t1full)
        call dec_fragment_time_get(times_fragopt)
        call dec_time_evaluate_efficiency_frag(AtomicFragment,times_fragopt,&
           &AtomicFragment%ccmodel,'Fragment optmization')
        return
     end if
     IF(OrbDistanceSpec)THEN
        call mem_alloc(SortedDistanceTableOrbAtomOcc,nocc)
        call mem_alloc(OrbOccDistTrackMyAtom,nocc)
        call GetSortedList(SortedDistanceTableOrbAtomOcc,OrbOccDistTrackMyAtom,&
             & mymolecule%DistanceTableOrbAtomOcc,nocc,natoms,MyAtom)
        call mem_alloc(SortedDistanceTableOrbAtomVirt,nunocc)
        call mem_alloc(OrbVirtDistTrackMyAtom,nunocc)
        call GetSortedList(SortedDistanceTableOrbAtomVirt,OrbVirtDistTrackMyAtom,&
             & mymolecule%DistanceTableOrbAtomVirt,nunocc,natoms,MyAtom)
     ENDIF

     ! Do fragment expansion at different level than target model?
     if(DECinfo%fragopt_exp_model /= DECinfo%ccmodel) then
        MyMolecule%ccmodel(MyAtom,Myatom) = DECinfo%fragopt_exp_model
     end if

     ! ======================================================================
     !                            Initial fragment
     ! ======================================================================
     ! Start fragment optimization by calculating initial fragment 
     IF(DECinfo%onlyoccpart) then
        !All Occupied orbitals assigned to atoms within 1.0 Angstrom of central atom are included
        init_Occradius = 1.0_realk/bohr_to_angstrom
     ELSE
        !All Occupied orbitals assigned to atoms within 3.0 Angstrom of central atom are included
        IF(FOT.GT.2.0E-5_realk)THEN 
           init_Occradius = 3.0_realk/bohr_to_angstrom
        ELSE
           init_Occradius = 3.0_realk/bohr_to_angstrom
        ENDIF
     ENDIF
     IF(DECinfo%onlyvirtpart) then
        !All Virtual orbitals assigned to atoms within 1.0 Angstrom of central atom are included
        init_Virtradius = 1.0_realk/bohr_to_angstrom
     ELSE
        !All Virtual orbitals assigned to atoms within 3.0 Angstrom of central atom are included
        IF(FOT.GT.2.0E-5_realk)THEN 
           init_Virtradius = 3.0_realk/bohr_to_angstrom
        ELSE
           init_Virtradius = 3.0_realk/bohr_to_angstrom
        ENDIF
     ENDIF

     IF(FockMatrixOrdering)THEN
        call mem_alloc(OccAOS,nocc)
        call mem_alloc(VirtAOS,nunocc)
        call mem_alloc(OldOccAOS,nocc)
        call mem_alloc(OldVirtAOS,nunocc)
        call InitialFragmentFockSpec(natoms,nocc,nunocc,&
             & FmaxOcc,OrbOccFockTrackMyAtom,&
             & FmaxVirt,OrbVirtFockTrackMyAtom,&
             & OccAOS,VirtAOS,OccOrbitals,UnOccOrbitals,MyAtom)
        call atomic_fragment_init_orbital_specific(MyAtom,nunocc, nocc, VirtAOS, &
             & OccAOS,OccOrbitals,UnoccOrbitals,MyMolecule,mylsitem,AtomicFragment,.true.,.false.)
        call fragment_energy_and_prop(AtomicFragment)
     ELSEIF(OrbDistanceSpec)THEN
        call mem_alloc(OccAOS,nocc)
        call mem_alloc(VirtAOS,nunocc)
        call mem_alloc(OldOccAOS,nocc)
        call mem_alloc(OldVirtAOS,nunocc)
        call InitialFragmentOrbitalSpec(natoms,nocc,nunocc,&
             & SortedDistanceTableOrbAtomOcc,OrbOccDistTrackMyAtom,&
             & SortedDistanceTableOrbAtomVirt,OrbVirtDistTrackMyAtom,&
             & init_Occradius,init_Virtradius,OccAOS,VirtAOS,&
             & OccOrbitals,UnOccOrbitals,MyAtom)
        call atomic_fragment_init_orbital_specific(MyAtom,nunocc, nocc, VirtAOS, &
             & OccAOS,OccOrbitals,UnoccOrbitals,MyMolecule,mylsitem,AtomicFragment,.true.,.false.)
        call fragment_energy_and_prop(AtomicFragment)
     ELSE
        call InitialFragment(natoms,nocc_per_atom,nunocc_per_atom,DistMyatom,&
             & init_Occradius, init_Virtradius, Occ_atoms,Virt_atoms)
        call get_fragment_and_Energy(MyAtom,natoms,Occ_Atoms,Virt_Atoms,&
             & MyMolecule,MyLsitem,nocc,nunocc,OccOrbitals,UnoccOrbitals,&
             & AtomicFragment)
     ENDIF
     ! Print initial fragment information
     call fragopt_print_info(AtomicFragment,LagEnergyDiff,OccEnergyDiff,VirtEnergyDiff,iter)

     call GetnAtomsWithOrb(natoms,nocc_per_atom,&
          & nunocc_per_atom,nAtomsWithOccOrb,nAtomsWithVirtOrb)

     StepsizeLoop2 = DECinfo%Frag_Exp_Size
     ExpandVirt = .TRUE.    !Expand both Occupied and Virtual Space untill convergence
     ExpandOcc  = .TRUE. 
     call FragmentExpansionProcedure(MyAtom,AtomicFragment,nAtoms, &
          & OccOrbitals,nOcc,UnoccOrbitals,nUnocc,&
          & MyMolecule,mylsitem,freebasisinfo,t1full,ExpandOcc,ExpandVirt,&
          & Occ_atoms,Virt_atoms,FOT,DistMyAtom,SortedDistMyAtom,&
          & DistTrackMyAtom, nocc_per_atom,nunocc_per_atom,&
          & StepsizeLoop2,LagEnergyOld, OccEnergyOld, VirtEnergyOld,&
          & nAtomsWithOccOrb,nAtomsWithVirtOrb,&
          & OrbDistanceSpec,OccAOS,VirtAOS,OldOccAOS,OldVirtAOS,&
          & SortedDistanceTableOrbAtomOcc,OrbOccDistTrackMyAtom,&
          & SortedDistanceTableOrbAtomVirt,OrbVirtDistTrackMyAtom,&
          & OrbOccFockTrackMyAtom,OrbVirtFockTrackMyAtom,&
          & FockMatrixOrdering)

     ! ======================================================================
     !             Transition from expansion to reduction loop
     ! ======================================================================

     ! Save contributions from individual local orbitals for current fragment
     ! **********************************************************************
     ! Note 1: The current fragment is larger than the converged fragment
     ! Note 2: This information is only used for local reduction procedure below (not fragment-adapted)
     call mem_alloc(OccContribs,nocc)
     call mem_alloc(VirtContribs,nunocc)
     OccContribs=0.0E0_realk
     VirtContribs=0.0E0_realk

     ! Contributions from local occupied orbitals 
     do i=1,AtomicFragment%noccAOS
        ! index of occupied AOS orbital "i" in list of ALL occupied orbitals in the molecule
        idx=AtomicFragment%occAOSidx(i)
        OccContribs(idx) = AtomicFragment%OccContribs(i)
     end do

     ! Contributions from local virtual orbitals
     do i=1,AtomicFragment%nunoccAOS
        ! index of virtual AOS orbital "i" in list of ALL virtual orbitals in the molecule
        idx=AtomicFragment%unoccAOSidx(i)
        VirtContribs(idx) = AtomicFragment%VirtContribs(i)
     end do


     ! Set AtomicFragment to be the converged fragment                                        
     ! ***********************************************
     ! Delete current fragment (which was too large)
     IF(OrbDistanceSpec)THEN
        IF(DECinfo%fragopt_exp_model .NE. DECinfo%fragopt_red_model)THEN
!           LagEnergyOldFull = AtomicFragment%LagFOP
!           OccEnergyOldFull = AtomicFragment%EoccFOP
!           VirtEnergyOldFull = AtomicFragment%EvirtFOP
           call atomic_fragment_free(AtomicFragment)
           call atomic_fragment_init_orbital_specific(MyAtom,nunocc, nocc,VirtAOS,&
                & OccAOS,OccOrbitals,UnoccOrbitals,MyMolecule,mylsitem,AtomicFragment,.true.,.false.)
        ELSE
           LagEnergyOld = AtomicFragment%LagFOP
           OccEnergyOld = AtomicFragment%EoccFOP
           VirtEnergyOld = AtomicFragment%EvirtFOP
        ENDIF
     ELSE
        call atomic_fragment_free(AtomicFragment)
        ! Init fragment with converged size
        call atomic_fragment_init_atom_specific(MyAtom,natoms,Virt_Atoms, &
             & Occ_Atoms,nocc,nunocc,OccOrbitals,UnoccOrbitals, &
             & MyMolecule,mylsitem,AtomicFragment,.true.,.false.)
     ENDIF
     ! Information for fragment-adapted orbitals
     ! *****************************************
     ! For practical reasons we now simply repeat the MP2 calculation to get all AOS amplitudes
     ! When properly tested, this might be fixed such that we do not need to repeat calcs.
     FragAdapt: if(DECinfo%fragadapt) then

        ! Get MP2 amplitudes for fragment
        ! *******************************
        ! Integrals (ai|bj)
        !call get_VOVO_integrals(AtomicFragment%mylsitem,AtomicFragment%nbasis,&
        !   & AtomicFragment%noccAOS,AtomicFragment%nunoccAOS,&
        !   & AtomicFragment%Cv, AtomicFragment%Co, g)
        !! Amplitudes
        call mp2_solver(AtomicFragment,g,t2,.false.)
        !call array4_free(g)

        call tensor_free(g)
        ! Get correlation density matrix for atomic fragment
        call calculate_MP2corrdens_frag(t2,AtomicFragment)
     end if FragAdapt


     ! Which model for reduction loop?
     ! *******************************
     MyMolecule%ccmodel(MyAtom,Myatom) = DECinfo%fragopt_red_model


     ! Save energies in converged space of local orbitals
     ! **************************************************
     if(DECinfo%fragopt_exp_model .eq. DECinfo%fragopt_red_model) then
        ! If the same model is used for expansion and reduction
        ! (e.g. using MP2 for expansion AND reduction  
        !   - or - using CCSD for expansion and reduction)
        ! then we can simply copy reference energy from expansion calculation
        AtomicFragment%LagFOP = LagEnergyOld
        AtomicFragment%EoccFOP = OccEnergyOld
        AtomicFragment%EvirtFOP = VirtEnergyOld

     else
        ! Different model in expansion and reduction steps - Calculate new reference energy
        ! for converged fragment from expansion loop.
        AtomicFragment%ccmodel = MyMolecule%ccmodel(MyAtom,Myatom)
        call fragment_energy_and_prop(AtomicFragment)
        LagEnergyDiff=0.0_realk
        OccEnergyDiff=0.0_realk
        VirtEnergyDiff=0.0_realk
        iter=0
        write(DECinfo%output,'(2a)') 'FOP Calculated ref atomic fragment energy for relevant CC model: ', &
           & DECinfo%cc_models(MyMolecule%ccmodel(MyAtom,Myatom))
        call fragopt_print_info(AtomicFragment,LagEnergyDiff,OccEnergyDiff,VirtEnergyDiff,iter)
        LagEnergyOld = AtomicFragment%LagFOP
        OccEnergyOld= AtomicFragment%EoccFOP
        VirtEnergyOld = AtomicFragment%EvirtFOP
     end if

     ! ======================================================================
     !                             Reduction loop
     ! ======================================================================


     write(DECinfo%output,*) ' FOP'
     write(DECinfo%output,*) ' FOP *************************************************'
     write(DECinfo%output,*) ' FOP ** Expansion has converged. We start reduction **'
     write(DECinfo%output,*) ' FOP *************************************************'
     write(DECinfo%output,*) ' FOP'


     WhichReductionScheme: if(DECinfo%fragadapt) then
        ! Reduce using fragment-adapted orbitals
        print *,"before",t2%initialized
        call fragopt_reduce_FOs(MyAtom,AtomicFragment, &
           &OccOrbitals,nOcc,UnoccOrbitals,nUnocc, &
           &MyMolecule,mylsitem,freebasisinfo,t2,max_iter_red)
        print *,"freeing tesn",t2%initialized
        call tensor_free(t2)
     else
        ! Reduce using local orbitals
        if(present(t1full)) then
           call fragopt_reduce_local_orbitals(MyAtom,AtomicFragment, &
                & OccOrbitals,nOcc,UnoccOrbitals,nUnocc,ncore,BinarySearch, &
                & MyMolecule,mylsitem,freebasisinfo,max_iter_red,OccContribs,&
                & VirtContribs,&
                & SortedDistanceTableOrbAtomOcc,OrbOccDistTrackMyAtom,&
                & SortedDistanceTableOrbAtomVirt,OrbVirtDistTrackMyAtom,&
                & OccAOS,VirtAOS,OldOccAOS,OldVirtAOS,OrbDistanceSpec,DistanceRemoval,&
                & TestOcc,BruteForce,FockMatrixOrdering,&
                & FmaxOcc,OrbOccFockTrackMyAtom,FmaxVirt,OrbVirtFockTrackMyAtom,&
                & t1full=t1full)
        else
           call fragopt_reduce_local_orbitals(MyAtom,AtomicFragment, &
              & OccOrbitals,nOcc,UnoccOrbitals,nUnocc,ncore,BinarySearch, &
              & MyMolecule,mylsitem,freebasisinfo,max_iter_red,OccContribs,VirtContribs,&
              & SortedDistanceTableOrbAtomOcc,OrbOccDistTrackMyAtom,&
              & SortedDistanceTableOrbAtomVirt,OrbVirtDistTrackMyAtom,&
              & OccAOS,VirtAOS,OldOccAOS,OldVirtAOS,OrbDistanceSpec,DistanceRemoval,&
              & TestOcc,BruteForce,FockMatrixOrdering,&
              & FmaxOcc,OrbOccFockTrackMyAtom,FmaxVirt,OrbVirtFockTrackMyAtom)
        end if

     end if WhichReductionScheme

     if(freebasisinfo) then
        call atomic_fragment_free_basis_info(AtomicFragment)
     end if
     call mem_dealloc(OccContribs)
     call mem_dealloc(VirtContribs)
     IF(OrbDistanceSpec)THEN
        call mem_dealloc(SortedDistanceTableOrbAtomOcc)
        call mem_dealloc(SortedDistanceTableOrbAtomVirt)
        call mem_dealloc(OrbOccDistTrackMyAtom)
        call mem_dealloc(OrbVirtDistTrackMyAtom)
        call mem_dealloc(OccAOS)
        call mem_dealloc(VirtAOS)
        call mem_dealloc(OldOccAOS)
        call mem_dealloc(OldVirtAOS)
     ENDIF
     IF(FockMatrixOrdering)THEN
        call mem_dealloc(OrbOccFockTrackMyAtom)
        call mem_dealloc(OrbVirtFockTrackMyAtom)
     ENDIF
     !IF((DECinfo%fragopt_exp_model.eq.DECinfo%fragopt_red_model).AND.(.NOT.DECinfo%fragadapt))THEN
     !   !the most correct energies - calculated for the biggest fragments used in expansion
     !   AtomicFragment%LagFOP = LagEnergyOldFull
     !   AtomicFragment%EoccFOP = OccEnergyOldFull 
     !   AtomicFragment%EvirtFOP = VirtEnergyOldFull
     !ENDIF
     ! Ensure that energies in fragment are set consistently
     call set_energies_decfrag_structure_fragopt(AtomicFragment)

     call dec_fragment_time_get(times_fragopt)
     call dec_time_evaluate_efficiency_frag(AtomicFragment,times_fragopt,AtomicFragment%ccmodel,'Fragment optmization')


     ! Restore the original CC model 
     ! (only relevant if expansion and/or reduction was done using the MP2 model, but it doesn't hurt)
     MyMolecule%ccmodel(MyAtom,Myatom) = DECinfo%ccmodel
     ! call lsquit('TEST DONE',-1)

  end subroutine optimize_atomic_fragment

  subroutine GetSortedFockMaxElements(FmaxOcc,FmaxVirt,nocc,nunocc,&
       & MyMolecule,ncore,OccOrbitals,UnOccOrbitals,MyAtom,&
       & OrbOccFockTrackMyAtom,OrbVirtFockTrackMyAtom)
    implicit none
    integer,intent(in) :: nocc,nunocc,ncore,MyAtom
    !> All occupied orbitals
    type(decorbital), dimension(nOcc), intent(in)      :: OccOrbitals
    !> All unoccupied orbitals
    type(decorbital), dimension(nUnocc), intent(in)    :: UnoccOrbitals
    !> Full molecule information          
    type(fullmolecule), intent(in) :: MyMolecule
    real(realk),intent(inout) :: FmaxOcc(nocc),FmaxVirt(nunocc)
    integer,intent(inout) :: OrbOccFockTrackMyAtom(nocc),OrbVirtFockTrackMyAtom(nunocc)
    !
    integer :: noccEOS,j,centralatom,nvirtEOS
    integer,pointer :: occEOSidx(:),virtEOSidx(:)

    noccEOS = 0
    do j=1,nOcc !note this include core orbitals
       CentralAtom=OccOrbitals(j)%centralatom
       if( CentralAtom .EQ. MyAtom)THEN
          noccEOS         = noccEOS + 1
       endif
    enddo
    call mem_alloc(occEOSidx,noccEOS)
    occEOSidx = 0
    noccEOS = 0
    do j=1,nOcc !note this include core orbitals
       CentralAtom=OccOrbitals(j)%centralatom
       if( CentralAtom .EQ. MyAtom)THEN
          noccEOS         = noccEOS + 1
          occEOSidx(noccEOS)=j
       endif
       OrbOccFockTrackMyAtom(j) = j
    enddo
    call GetFockMaxElement(FmaxOcc,nocc,MyMolecule%ppFock,noccEOS,occEOSidx,ncore)
    call mem_dealloc(occEOSidx)
    call real_inv_sort_with_tracking(FmaxOcc,OrbOccFockTrackMyAtom,nocc)

    nvirtEOS = 0
    do j=1,nunocc !note this include core orbitals
       CentralAtom=UnOccOrbitals(j)%centralatom
       if( CentralAtom .EQ. MyAtom)THEN
          nvirtEOS         = nvirtEOS + 1
       endif
    enddo
    call mem_alloc(virtEOSidx,nvirtEOS)
    virtEOSidx = 0
    nvirtEOS = 0
    do j=1,nunocc !note this include core orbitals
       CentralAtom=UnOccOrbitals(j)%centralatom
       if( CentralAtom .EQ. MyAtom)THEN
          nvirtEOS         = nvirtEOS + 1
          virtEOSidx(nvirtEOS)=j
       endif
       OrbVirtFockTrackMyAtom(j) = j
    enddo
    call GetFockMaxElement(FmaxVirt,nunocc,MyMolecule%qqFock,nvirtEOS,virtEOSidx,0)
    call mem_dealloc(virtEOSidx)
    call real_inv_sort_with_tracking(FmaxVirt,OrbVirtFockTrackMyAtom,nunocc)
  end subroutine GetSortedFockMaxElements

  subroutine InitialFragmentOrbitalSpec(natoms,nocc,nvirt,&
       & SortedDistanceTableOrbAtomOcc,OrbOccDistTrackMyAtom,&
       & SortedDistanceTableOrbAtomVirt,OrbVirtDistTrackMyAtom,&
       & init_Occradius,init_Virtradius,OccAOS,VirtAOS,&
       & OccOrbitals,UnoccOrbitals,MyAtom)
    implicit none
    integer,intent(in) :: natoms,nocc,nvirt
    real(realk),intent(in) :: SortedDistanceTableOrbAtomOcc(nocc)
    real(realk),intent(in) :: SortedDistanceTableOrbAtomVirt(nvirt)
    integer,intent(in) :: OrbOccDistTrackMyAtom(nocc)
    integer,intent(in) :: OrbVirtDistTrackMyAtom(nvirt)
    logical,intent(inout) :: OccAOS(nocc), VirtAOS(nvirt)
    real(realk),intent(in) :: init_Occradius,init_Virtradius
    type(decorbital), dimension(nOcc), intent(in) :: OccOrbitals
    type(decorbital), dimension(nvirt), intent(in) :: UnoccOrbitals
    integer,intent(in) :: MyAtom
    !
    integer :: i,ii
    OccAOS=.false.
    do i=1,nocc
       IF(init_Occradius.LT.SortedDistanceTableOrbAtomOcc(i))EXIT
       ii = OrbOccDistTrackMyAtom(i)
       OccAOS(ii)=.true.
    enddo
    VirtAOS=.false.
    do i=1,nvirt
       IF(init_Virtradius.LT.SortedDistanceTableOrbAtomVirt(i))EXIT
       ii = OrbVirtDistTrackMyAtom(i)
       VirtAOS(ii)=.true.
    enddo
    !include EOS space
    do i=1,nocc
       if(OccOrbitals(i)%centralatom.EQ.MyAtom) then
          OccAOS(i)=.true.
       end if
    end do
    do i=1,nvirt
       if(UnOccOrbitals(i)%centralatom.EQ.MyAtom) then
          VirtAOS(i)=.true.
       end if
    end do
  end subroutine InitialFragmentOrbitalSpec

  subroutine InitialFragmentFockSpec(natoms,nocc,nvirt,&
       & SortedFmaxOcc,OrbOccFockTrackMyAtom,&
       & SortedFmaxVirt,OrbVirtFockTrackMyAtom,&
       & OccAOS,VirtAOS,OccOrbitals,UnoccOrbitals,MyAtom)
    implicit none
    integer,intent(in) :: natoms,nocc,nvirt
    real(realk),intent(in) :: SortedFmaxOcc(nocc)
    real(realk),intent(in) :: SortedFmaxVirt(nvirt)
    integer,intent(in) :: OrbOccFockTrackMyAtom(nocc)
    integer,intent(in) :: OrbVirtFockTrackMyAtom(nvirt)
    logical,intent(inout) :: OccAOS(nocc), VirtAOS(nvirt)
    type(decorbital), dimension(nOcc), intent(in) :: OccOrbitals
    type(decorbital), dimension(nvirt), intent(in) :: UnoccOrbitals
    integer,intent(in) :: MyAtom
    !
    integer :: i,ii
    !For FOT=3 include all orbitals that has a Max Fock matrix element
    !F(i,eosorbital) > 10**-4
    OccAOS=.false.
    do i=1,nocc
       IF(SortedFmaxOcc(i).GT.DECinfo%FOT*1.0E-1_realk)EXIT
       ii = OrbOccFockTrackMyAtom(i)
       OccAOS(ii)=.true.
    enddo
    VirtAOS=.false.
    do i=1,nvirt
       IF(SortedFmaxVirt(i).GT.DECinfo%FOT*1.0E-1_realk)EXIT
       ii = OrbVirtFockTrackMyAtom(i)
       VirtAOS(ii)=.true.
    enddo
    !include EOS space
    do i=1,nocc
       if(OccOrbitals(i)%centralatom.EQ.MyAtom) then
          OccAOS(i)=.true.
       end if
    end do
    do i=1,nvirt
       if(UnOccOrbitals(i)%centralatom.EQ.MyAtom) then
          VirtAOS(i)=.true.
       end if
    end do
  end subroutine InitialFragmentFockSpec

  subroutine ExpandFragmentOrbitalSpec(natoms,nocc,nvirt,&
       & SortedDistanceTableOrbAtomOcc,OrbOccDistTrackMyAtom,&
       & SortedDistanceTableOrbAtomVirt,OrbVirtDistTrackMyAtom,&
       & StepsizeLoop,OccAOS,VirtAOS,ExpandFragmentConverged)
    implicit none
    integer,intent(in) :: natoms,nocc,nvirt,StepsizeLoop
    real(realk),intent(in) :: SortedDistanceTableOrbAtomOcc(nocc)
    real(realk),intent(in) :: SortedDistanceTableOrbAtomVirt(nvirt)
    integer,intent(in) :: OrbOccDistTrackMyAtom(nocc)
    integer,intent(in) :: OrbVirtDistTrackMyAtom(nvirt)
    logical,intent(inout) :: OccAOS(nocc), VirtAOS(nvirt),ExpandFragmentConverged
    !
    integer :: i,ii,StepsizeLoop6,nOrb1,nOrb2,Xocc,Yvirt,jj
    real(realk) :: nov, maxvirtdist    
    integer :: AverageOrbitalPerAtoms
    AverageOrbitalPerAtoms = CEILING((nocc+nvirt*1.0E0_realk)/natoms)
    StepsizeLoop6 = StepsizeLoop*AverageOrbitalPerAtoms 

    IF(DECinfo%onlyOccPart)THEN
       !increase with StepsizeLoop6 virtual orbitals
       nOrb2 = COUNT(VirtAOS)
       do i=nOrb2+1,MIN(nOrb2+StepsizeLoop6,nvirt)
          ii = OrbVirtDistTrackMyAtom(i)
          VirtAOS(ii)=.true.
       enddo
       !increase Occ with Radius determined by the biggest Virtual radius
       jj = MIN(nOrb2+StepsizeLoop6,nvirt)
       maxvirtdist = SortedDistanceTableOrbAtomVirt(jj)
       nOrb1 = COUNT(OccAOS)
       do i=1,nocc
          if(SortedDistanceTableOrbAtomOcc(i) .le. maxvirtdist) then
             ii = OrbOccDistTrackMyAtom(i)
             OccAOS(ii)=.true.
          end if
       enddo
    ELSEIF(DECinfo%onlyVirtPart)THEN
       !increase with StepsizeLoop6 Occupied orbitals
       nOrb1 = COUNT(OccAOS)
       do i=nOrb1+1,MIN(nOrb1+StepsizeLoop6,nocc)
          ii = OrbOccDistTrackMyAtom(i)
          OccAOS(ii)=.true.
       enddo
       !increase Virt with Radius determined by the biggest Occupied radius
       jj = MIN(nOrb1+StepsizeLoop6,nocc)
       nOrb2 = COUNT(VirtAOS)
       do i=nOrb2+1,MIN(nOrb2+StepsizeLoop6,nvirt)
          IF(SortedDistanceTableOrbAtomOcc(jj).LT.SortedDistanceTableOrbAtomVirt(i))EXIT
          ii = OrbVirtDistTrackMyAtom(i)
          VirtAOS(ii)=.true.
       enddo
    ELSE
       nov = (nocc*1.0E0_realk/(nvirt*1.0E0_realk))
       Xocc = NINT(2*StepsizeLoop6*nov)
       Yvirt = 2*StepsizeLoop6-Xocc
       !increase with Yvirt Virtual Orbitals
       nOrb2 = COUNT(VirtAOS)
       do i=nOrb2+1,MIN(nOrb2+Yvirt,nvirt)
          ii = OrbVirtDistTrackMyAtom(i)
          VirtAOS(ii)=.true.
       enddo
       !increase with Xocc Occupied Orbitals
       nOrb1 = COUNT(OccAOS)
       do i=nOrb1+1,MIN(nOrb1+Xocc,nocc)
          ii = OrbOccDistTrackMyAtom(i)
          OccAOS(ii)=.true.
       enddo
    ENDIF
    ExpandFragmentConverged=COUNT(OccAOS).EQ.nOcc.AND.COUNT(VirtAOS).EQ.nvirt
  end subroutine ExpandFragmentOrbitalSpec

  subroutine ExpandFragmentOrbitalSpecFock(natoms,nocc,nvirt,&
       & OrbOccFockTrackMyAtom,OrbVirtFockTrackMyAtom,&
       & StepsizeLoop,OccAOS,VirtAOS,ExpandFragmentConverged)
    implicit none
    integer,intent(in) :: natoms,nocc,nvirt,StepsizeLoop
    integer,intent(in) :: OrbOccFockTrackMyAtom(nocc)
    integer,intent(in) :: OrbVirtFockTrackMyAtom(nvirt)
    logical,intent(inout) :: OccAOS(nocc), VirtAOS(nvirt),ExpandFragmentConverged
    !
    integer :: i,ii,StepsizeLoop6,nOrb1,nOrb2,Xocc,Yvirt,jj
    real(realk) :: nov    
    integer :: AverageOrbitalPerAtoms
    AverageOrbitalPerAtoms = CEILING((nocc+nvirt*1.0E0_realk)/natoms)
    StepsizeLoop6 = StepsizeLoop*AverageOrbitalPerAtoms 
    nov = (nocc*1.0E0_realk/(nvirt*1.0E0_realk))
    Xocc = NINT(2*StepsizeLoop6*nov)
    Yvirt = 2*StepsizeLoop6-Xocc
    !increase with Yvirt Virtual Orbitals
    nOrb2 = COUNT(VirtAOS)
    do i=nOrb2+1,MIN(nOrb2+Yvirt,nvirt)
       ii = OrbVirtFockTrackMyAtom(i)
       VirtAOS(ii)=.true.
    enddo
    !increase with Xocc Occupied Orbitals
    nOrb1 = COUNT(OccAOS)
    do i=nOrb1+1,MIN(nOrb1+Xocc,nocc)
       ii = OrbOccFockTrackMyAtom(i)
       OccAOS(ii)=.true.
    enddo
    ExpandFragmentConverged=COUNT(OccAOS).EQ.nOcc.AND.COUNT(VirtAOS).EQ.nvirt
  end subroutine ExpandFragmentOrbitalSpecFock
  
  subroutine FragmentExpansionProcedure(MyAtom,AtomicFragment,nAtoms, &
          & OccOrbitals,nOcc,UnoccOrbitals,nUnocc,&
          & MyMolecule,mylsitem,freebasisinfo,t1full,ExpandOcc,ExpandVirt,&
          & Occ_atoms,Virt_atoms,FOT,DistMyAtom,SortedDistMyAtom,&
          & DistTrackMyAtom, nocc_per_atom,nunocc_per_atom,&
          & StepsizeLoop2,LagEnergyOld, OccEnergyOld, VirtEnergyOld,&
          & nAtomsWithOccOrb,nAtomsWithVirtOrb,&
          & OrbDistanceSpec,OccAOS,VirtAOS,OldOccAOS,OldVirtAOS,&
          & SortedDistanceTableOrbAtomOcc,OrbOccDistTrackMyAtom,&
          & SortedDistanceTableOrbAtomVirt,OrbVirtDistTrackMyAtom,&
          & OrbOccFockTrackMyAtom,OrbVirtFockTrackMyAtom,&
          & FockMatrixOrdering)
       implicit none
     !> Number of occupied orbitals in molecule
     integer, intent(in) :: nOcc
     !> Number of unoccupied orbitals in molecule
     integer, intent(in) :: nunocc
     !> Number of atoms in molecule
     integer, intent(in) :: natoms
     !> Central atom in molecule
     integer, intent(in) :: MyAtom
     !> Atomic fragment to be optimized
     type(decfrag),intent(inout)        :: AtomicFragment
     !> All occupied orbitals
     type(decorbital), dimension(nOcc), intent(in)      :: OccOrbitals
     !> All unoccupied orbitals
     type(decorbital), dimension(nUnocc), intent(in)    :: UnoccOrbitals
     !> Full molecule information
     type(fullmolecule), intent(inout) :: MyMolecule
     !> Integral information
     type(lsitem), intent(inout)       :: mylsitem
     !> Delete fragment basis information ("expensive box in decfrag type") at exit?
     logical,intent(in) :: freebasisinfo
     !> t1 amplitudes for full molecule to be updated (only used when DECinfo%SinglesPolari is set)
     type(array2),intent(inout),optional :: t1full
     !> Should we expand Occupied space
     logical,intent(in) :: ExpandOcc
     !> Should we expand Virtual space
     logical,intent(in) :: ExpandVirt,OrbDistanceSpec
     logical, dimension(natoms)     :: Occ_atoms,Virt_atoms
     logical, dimension(natoms)     :: OccOld,VirtOld       !previous Fragment
     real(realk),intent(in)         :: FOT
     integer,intent(in)             :: StepsizeLoop2,nAtomsWithOccOrb,nAtomsWithVirtOrb
     real(realk),dimension(natoms),intent(in)  :: DistMyAtom,SortedDistMyAtom
     integer,dimension(natoms),intent(in)      :: DistTrackMyAtom, nocc_per_atom,nunocc_per_atom
     real(realk),intent(inout)      :: LagEnergyOld, OccEnergyOld, VirtEnergyOld
     logical :: FockMatrixOrdering
     logical,pointer :: OccAOS(:),VirtAOS(:),OldOccAOS(:),OldVirtAOS(:) !not always allocated only if OrbDistanceSpec 
     real(realk),pointer :: SortedDistanceTableOrbAtomOcc(:)
     real(realk),pointer :: SortedDistanceTableOrbAtomVirt(:)
     integer,pointer :: OrbOccDistTrackMyAtom(:)
     integer,pointer :: OrbVirtDistTrackMyAtom(:)
     integer,pointer :: OrbOccFockTrackMyAtom(:)
     integer,pointer :: OrbVirtFockTrackMyAtom(:)
     !Local variables
     real(realk)  :: LagEnergyDiff, OccEnergyDiff,VirtEnergyDiff,EnergyDiff
     integer      :: iter,i,idx,StepsizeLoop3
     integer      :: max_iter_red
     logical :: expansion_converged,ExpandFragmentConverged

     StepsizeLoop3 = StepsizeLoop2
     IF(OrbDistanceSpec)THEN
        OldOccAOS = OccAOS
        OldVirtAOS = VirtAOS
     ENDIF
     if(DECinfo%frozencore)THEN
      IF(FockMatrixOrdering.OR.OrbDistanceSpec)THEN
         !exclude core orbitals in OccAOS
         !they are removed in atomic_fragment_init_orbital_specific
         !anyway
         do iter=1,MyMolecule%ncore
            OccAOS(iter) = .FALSE.
         enddo
      ENDIF
     ENDIF
 ! ======================================================================
 !                   Expansion loop
 ! ======================================================================
     EXPANSION_LOOP: do iter = 1,DECinfo%maxiter

        ! Save information for current fragment (in case current fragment is the final one)
        OccOld=Occ_atoms;VirtOld=Virt_atoms
        LagEnergyOld = AtomicFragment%LagFOP
        OccEnergyOld = AtomicFragment%EoccFOP
        VirtEnergyOld = AtomicFragment%EvirtFOP

        IF(FockMatrixOrdering)THEN
           call ExpandFragmentOrbitalSpecFock(natoms,nocc,nunocc,&
                & OrbOccFockTrackMyAtom,OrbVirtFockTrackMyAtom,&
                & StepsizeLoop3,OccAOS,VirtAOS,ExpandFragmentConverged)
           call atomic_fragment_free(AtomicFragment)
           call atomic_fragment_init_orbital_specific(MyAtom,nunocc,nocc,VirtAOS, &
                & OccAOS,OccOrbitals,UnoccOrbitals,MyMolecule,mylsitem,AtomicFragment,.true.,.false.)
           call fragment_energy_and_prop(AtomicFragment)
        ELSEIF(OrbDistanceSpec)THEN
           call ExpandFragmentOrbitalSpec(natoms,nocc,nunocc,&
                & SortedDistanceTableOrbAtomOcc,OrbOccDistTrackMyAtom,&
                & SortedDistanceTableOrbAtomVirt,OrbVirtDistTrackMyAtom,&
                & StepsizeLoop3,OccAOS,VirtAOS,ExpandFragmentConverged)
           call atomic_fragment_free(AtomicFragment)
           call atomic_fragment_init_orbital_specific(MyAtom,nunocc,nocc,VirtAOS, &
                & OccAOS,OccOrbitals,UnoccOrbitals,MyMolecule,mylsitem,AtomicFragment,.true.,.false.)
           call fragment_energy_and_prop(AtomicFragment)
        ELSE
           ! Expand fragment and get new energy
           call Expandfragment(Occ_atoms,Virt_atoms,DistTrackMyAtom,natoms,&
                & nocc_per_atom,nunocc_per_atom,ExpandFragmentConverged,&
                & nAtomsWithOccOrb,nAtomsWithVirtOrb,StepsizeLoop3,&
                & ExpandOcc,ExpandVirt)           
           call atomic_fragment_free(AtomicFragment)
           call get_fragment_and_Energy(MyAtom,natoms,Occ_Atoms,Virt_Atoms,&
                & MyMolecule,MyLsitem,nocc,nunocc,OccOrbitals,UnoccOrbitals,&
                & AtomicFragment)
        ENDIF
        ! Energy differences
        LagEnergyDiff=abs(LagEnergyOld-AtomicFragment%LagFOP)
        OccEnergyDiff=abs(OccEnergyOld-AtomicFragment%EoccFOP)
        VirtEnergyDiff=abs(VirtEnergyOld-AtomicFragment%EvirtFOP)

        call fragopt_print_info(AtomicFragment,LagEnergyDiff,OccEnergyDiff,VirtEnergyDiff,iter)

        ! Test if fragment energy (or energies) are converged to FOT precision
        call fragopt_check_convergence(LagEnergyDiff,OccEnergyDiff,VirtEnergyDiff,&
           &FOT,expansion_converged)
        IF(ExpandFragmentConverged)THEN
           IF(.NOT.expansion_converged)THEN
              WRITE(DECinfo%output,*)'Expansion Include Full Molecule'
              OccOld = Occ_atoms;VirtOld = Virt_atoms
              LagEnergyOld = AtomicFragment%LagFOP
              OccEnergyOld = AtomicFragment%EoccFOP
              VirtEnergyOld = AtomicFragment%EvirtFOP
              expansion_converged = .TRUE.
           ENDIF
        ENDIF

        ! Exit loop if we are converged
        ExpansionConvergence: if(expansion_converged) then
           Occ_atoms = OccOld;Virt_atoms = VirtOld
           IF(ExpandOcc.AND.ExpandVirt)THEN
              write(DECinfo%output,*) 'FOP Fragment expansion converged in iteration ', iter
           ELSEIF(ExpandOcc)THEN
              write(DECinfo%output,*) 'FOP Occupied Fragment expansion converged in iteration ', iter
           ELSEIF(ExpandVirt)THEN
              write(DECinfo%output,*) 'FOP Virtual Fragment expansion converged in iteration ', iter
           ENDIF
           exit EXPANSION_LOOP
        else
           IF(OrbDistanceSpec)THEN
              OldOccAOS = OccAOS
              OldVirtAOS = VirtAOS
           ENDIF
        end if ExpansionConvergence
     end do EXPANSION_LOOP


     ! Check that expansion loop is converged
     if(.not. expansion_converged) then
        write(DECinfo%output,*) 'Number of expansion steps = ', DECinfo%MaxIter
        call lsquit('Fragment expansion did not converge! &
        & Try to increase the number of expansion steps using the .MaxIter keyword',DECinfo%output)
     end if

   end subroutine FragmentExpansionProcedure

  !> Given a converged atomic fragment using local orbitals, determine fragment
  !> of reduced size, where the energy error compared to the original converged fragment is
  !> below the FOT.
  !> \date September 2013
  !> \author Kasper Kristensen
  subroutine fragopt_reduce_FOs(MyAtom,AtomicFragment, &
       &OccOrbitals,nOcc,UnoccOrbitals,nUnocc, &
       &MyMolecule,mylsitem,freebasisinfo,t2tens,max_iter_red)
    implicit none
    !> Number of occupied orbitals in molecule
    integer, intent(in) :: nOcc
    !> Number of unoccupied orbitals in molecule
    integer, intent(in) :: nunocc
    !> Central atom in molecule
    integer, intent(in) :: MyAtom
    !> Atomic fragment to be optimized
    type(decfrag),intent(inout)        :: AtomicFragment
    !> All occupied orbitals
    type(decorbital), dimension(nOcc), intent(in)      :: OccOrbitals
    !> All unoccupied orbitals
    type(decorbital), dimension(nUnocc), intent(in)    :: UnoccOrbitals
    !> Full molecule information
    type(fullmolecule), intent(in) :: MyMolecule
    !> Integral information
    type(lsitem), intent(inout)       :: mylsitem
    !> Delete fragment basis information ("expensive box in decfrag type") at exit?
    logical,intent(in) :: freebasisinfo
    !> Doubles amplitudes for converged fragment (local orbital basis, index order: a,i,b,j)
    !> At output, the virtual indices will have been transformed to the (smaller) FO orbital space
    type(tensor),intent(inout) :: t2tens
    !> Maximum number of reduction steps
    integer,intent(in) :: max_iter_red
    integer :: ov,iter,Nold,Nnew,nocc_exp,nvirt_exp
    logical :: reduction_converged,ReductionPossible(2)
    real(realk) :: FOT,LagEnergyOld,OccEnergyOld,VirtEnergyOld
    real(realk)                    :: LagEnergyDiff, OccEnergyDiff,VirtEnergyDiff
    type(array4)  :: t2
    character(4) :: stens_atype
    integer :: os,vs

    ! Init stuff
    ! **********
    FOT = DECinfo%FOT
    ! Reference energies for converged fragment
    LagEnergyOld = AtomicFragment%LagFOP
    OccEnergyOld = AtomicFragment%EoccFOP
    VirtEnergyOld = AtomicFragment%EvirtFOP
    ! Dimensions for converged fragment
    nocc_exp = AtomicFragment%noccAOS
    nvirt_exp = AtomicFragment%nunoccAOS


    
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
          call get_symm_tensor_segmenting_simple(t2%dims(2),t2%dims(1),os,vs)
          call tensor_minit(t2tens,t2%dims,4, atype = stens_atype, tdims=[vs,os,vs,os])
          call tensor_convert(t2%val,t2tens)
          call array4_free(t2)

       end if
       if(DECinfo%onlyVirtpart) then
          WRITE(DECinfo%output,*)'WARNING FOs usign onlyVirtpart not tested'
          call lsquit('Error FOs usign onlyVirtpart not tested',-1)
       endif
       REDUCTION_LOOP: do iter=1,max_iter_red

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
                Nold = AtomicFragment%nunoccAOS
             end if
             ! Threshold for throwing away fragment-adapted orbitals (start with FOT value)
             AtomicFragment%RejectThr(ov) = FOT
          else  
             ! Number of occ or virt fragment-adapted orbitals from previous step
             if(ov==1) then
                Nold = AtomicFragment%noccAOS
             else
                Nold = AtomicFragment%nunoccAOS
             end if
             ! decrease rejection threshold by a factor 10 in each step
             AtomicFragment%RejectThr(ov) = AtomicFragment%RejectThr(ov)/10.0_realk
          end if
          if(ov==1) then
             write(DECinfo%output,'(a,ES13.5)') 'FOP occ rejection threshold  : ',&
                  & AtomicFragment%RejectThr(ov)
          else
             write(DECinfo%output,'(a,ES13.5)') 'FOP virt rejection threshold  : ',&
                  & AtomicFragment%RejectThr(ov)
          end if


          ! Make fragment-adapted fragment with smaller AOS according to rejection threshold
          ! *********************************************************************************
          call fragment_adapted_transformation_matrices(AtomicFragment)
          if(ov==1) then
             Nnew = AtomicFragment%noccAOS
          else
             Nnew = AtomicFragment%nunoccAOS
          end if


          ! Cycle if the AOS was not decreased in the rejection procedure
          ! *************************************************************
          if( Nold==Nnew ) then
             if(iter == max_iter_red) then
                if(ov==2) then  ! not possible to reduce virtual space, try occupied space
                   cycle OCC_OR_VIRT

                else ! not possible to decrease occ space
                   goto 1000
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
          LagEnergyDiff=abs(AtomicFragment%LagFOP-LagEnergyOld)
          OccEnergyDiff=abs(AtomicFragment%EoccFOP-OccEnergyOld)
          VirtEnergyDiff=abs(AtomicFragment%EvirtFOP-VirtEnergyOld)
          call fragopt_print_info(AtomicFragment,LagEnergyDiff,OccEnergyDiff,VirtEnergyDiff,iter)
          call fragopt_check_convergence(LagEnergyDiff,OccEnergyDiff,VirtEnergyDiff,&
               &FOT,reduction_converged)


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

1000      continue


          ! Special case: No reduction possible for neither occ nor virt space
          ! ******************************************************************
          if ( (iter == max_iter_red) .and. (ov==1) .and. &
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
         & AtomicFragment%nunoccFA
    write(DECinfo%output,'(1X,a,i4)')    'FOP Done: Number of orbitals in occ total  :', &
         & AtomicFragment%noccFA
    write(DECinfo%output,'(1X,a,i4)')     'FOP Done: Number of basis functions        :', &
         & AtomicFragment%nbasis
    if(.not. DECinfo%onlyvirtpart) then
       write(DECinfo%output,'(1X,a,f16.10)') 'FOP Done: Occupied Fragment energy         :', &
            & AtomicFragment%EoccFOP
    endif
    if(.not. DECinfo%onlyoccpart) then
       write(DECinfo%output,'(1X,a,f16.10)') 'FOP Done: Virtual Fragment energy          :', &
            & AtomicFragment%EvirtFOP
    endif
    if(.NOT.(DECinfo%onlyoccpart.OR. DECinfo%onlyvirtpart))then
       write(DECinfo%output,'(1X,a,f16.10)') 'FOP Done: Lagrangian Fragment energy       :', &
            & AtomicFragment%LagFOP
    end if
    write(DECinfo%output,'(1X,a,g14.2)')  'FOP Done: Occupied reduction threshold     :', &
         & AtomicFragment%RejectThr(1)
    write(DECinfo%output,'(1X,a,g14.2)')  'FOP Done: Virtual  reduction threshold     :', &
         & AtomicFragment%RejectThr(2)
    write(DECinfo%output,'(1X,a,i6,a,i6,a,f5.2,a)')  'FOP Done: Occupied reduction removed ', &
         & nocc_exp-AtomicFragment%noccFA, ' of ', nocc_exp, ' orbitals ( ', &
         & (nocc_exp-AtomicFragment%noccFA)*100.0_realk/nocc_exp, ' %)'
    write(DECinfo%output,'(1X,a,i6,a,i6,a,f5.2,a)')  'FOP Done: Virtual  reduction removed ', &
         & nvirt_exp-AtomicFragment%nunoccFA, ' of ', nvirt_exp, ' orbitals ( ', &
         & (nvirt_exp-AtomicFragment%nunoccFA)*100.0_realk/nvirt_exp, ' %)'
    write(DECinfo%output,'(1X,a,/)') 'FOP========================================================='
    write(DECinfo%output,*) 'FOP'

  end subroutine fragopt_reduce_FOs



  !> Given a converged atomic fragment, remove individual orbitals with the smallest energy
  !> contributions, ensuring that the error introduced compared to the converged fragment
  !> is smaller than the FOT.
  !> \date September 2013
  !> \author Kasper Kristensen
  subroutine fragopt_reduce_local_orbitals(MyAtom,AtomicFragment, &
       & OccOrbitals,nOcc,UnoccOrbitals,nUnocc,ncore,BinarySearch, &
       & MyMolecule,mylsitem,freebasisinfo,max_iter_red,OccContribs,VirtContribs,&
       & SortedDistanceTableOrbAtomOcc,OrbOccDistTrackMyAtom,&
       & SortedDistanceTableOrbAtomVirt,OrbVirtDistTrackMyAtom,&
       & ExpOccAOS,ExpVirtAOS,ExpOldOccAOS,ExpOldVirtAOS,OrbDistanceSpec,&
       & DistanceRemoval,TestOcc,BruteForce,FockMatrixOrdering,&
       & FmaxOcc,OrbOccFockTrackMyAtom,FmaxVirt,OrbVirtFockTrackMyAtom,t1full)
    implicit none
    !> Number of occupied orbitals in molecule
    integer, intent(in) :: nOcc
    !> Number of unoccupied orbitals in molecule
    integer, intent(in) :: nunocc
    !> Number of core occupied orbitals in molecule
    integer, intent(in) :: ncore
    !> Central atom in molecule
    integer, intent(in) :: MyAtom
    !> Atomic fragment to be optimized
    type(decfrag),intent(inout)        :: AtomicFragment
    !> All occupied orbitals
    type(decorbital), dimension(nOcc), intent(in)      :: OccOrbitals
    !> All unoccupied orbitals
    type(decorbital), dimension(nUnocc), intent(in)    :: UnoccOrbitals
    !> Perform binary search? 
    logical,intent(in) :: BinarySearch
    !> Full molecule information
    type(fullmolecule), intent(in) :: MyMolecule
    !> Integral information
    type(lsitem), intent(inout)       :: mylsitem
    !> Delete fragment basis information ("expensive box in decfrag type") at exit?
    logical,intent(in) :: freebasisinfo
    !> Maximum number of reduction steps
    integer,intent(in) :: max_iter_red
    !> Contributions from individual occupied and virtual orbitals
    real(realk),intent(in) :: OccContribs(nocc), VirtContribs(nunocc)
    !> t1 amplitudes for full molecule to be updated (only used when DECinfo%SinglesPolari is set)
    type(array2),intent(inout),optional :: t1full
    real(realk)  :: RejectThresh,FOT
    logical,intent(in) :: OrbDistanceSpec,DistanceRemoval,TestOcc,BruteForce,FockMatrixOrdering
    logical,pointer :: ExpOccAOS(:),ExpVirtAOS(:),ExpOldOccAOS(:),ExpOldVirtAOS(:) !only alloc if OrbDistanceSpec
    logical,pointer :: OccAOS_old(:), VirtAOS_old(:), OccAOS_new(:), VirtAOS_new(:), &
         & OccAOS_orig(:),VirtAOS_orig(:)
    real(realk),pointer :: SortedDistanceTableOrbAtomOcc(:)
    real(realk),pointer :: SortedDistanceTableOrbAtomVirt(:)
    integer,pointer :: OrbOccDistTrackMyAtom(:)
    integer,pointer :: OrbVirtDistTrackMyAtom(:)
    logical :: reduction_converged
    integer :: i,iter,nocc_old,nvirt_old,nocc_new,nvirt_new,nocc_orig,nvirt_orig,ii
    integer :: nHigherOcc,nHigherVirt,nLowerOcc,nLowerVirt,loopI,nDimOcc,nDimVirt
    logical :: bin_reduction_converged,ModVirt,bin_virt_conv,bin_occ_conv,O2V2choice
    real(realk)  :: LagEnergyDiff, OccEnergyDiff,VirtEnergyDiff
    real(realk)  :: LagEnergyOld, OccEnergyOld,VirtEnergyOld
    real(realk),pointer :: SortedOccContribs(:),SortedVirtContribs(:)
    integer,pointer :: TrackSortedOccContribs(:),TrackSortedVirtContribs(:)
    real(realk),pointer :: BruteForceOccContribs(:),BruteForceVirtContribs(:)
    integer,pointer :: BruteForceTrackListOcc(:),BruteForceTrackListVirt(:)
    real(realk) :: FmaxOcc(nocc),FmaxVirt(nunocc)
    integer,pointer :: OrbOccFockTrackMyAtom(:),OrbVirtFockTrackMyAtom(:)
    integer :: k,j,noccEOS,nunoccEOS,nO2V2_ModOcc,nO2V2_ModVirt
    real(realk)  :: LagEnergy,OccEnergy,VirtEnergy

    ! Initialize logical vectors controlling occupied and virtual AOS during reduction scheme
    call mem_alloc(OccAOS_old,nocc)
    call mem_alloc(VirtAOS_old,nunocc)
    call mem_alloc(OccAOS_new,nocc)
    call mem_alloc(VirtAOS_new,nunocc)
    call mem_alloc(OccAOS_orig,nocc)
    call mem_alloc(VirtAOS_orig,nunocc)
    OccAOS_old=.false.
    VirtAOS_old=.false.
    OccAOS_new=.false.
    VirtAOS_new=.false.
    OccAOS_orig=.false.
    VirtAOS_orig=.false.
    do i=1,AtomicFragment%noccAOS
       OccAOS_orig(AtomicFragment%occAOSidx(i)) = .true.
    end do
    do i=1,AtomicFragment%nunoccAOS
       VirtAOS_orig(AtomicFragment%unoccAOSidx(i)) = .true.
    end do
    OccAOS_old = OccAOS_orig
    VirtAOS_old = VirtAOS_orig
    nocc_orig = AtomicFragment%noccAOS
    nvirt_orig = AtomicFragment%nunoccAOS
    LagEnergyOld = AtomicFragment%LagFOP
    OccEnergyOld = AtomicFragment%EoccFOP
    VirtEnergyOld = AtomicFragment%EvirtFOP

    noccEOS = AtomicFragment%noccEOS
    nunoccEOS = AtomicFragment%nunoccEOS

    O2V2choice = .FALSE.

    IF(BruteForce)THEN
       call mem_alloc(BruteForceOccContribs,nocc)
       BruteForceOccContribs = 0.0E0_realk
       call mem_alloc(BruteForceTrackListOcc,nocc)
       do I=1,nocc
          BruteForceTrackListOcc(I) = I
       enddo
!       IF(DECinfo%OnlyOccPart)THEN
          !determine Brute Force Occupied Energy contribution
          WRITE(DECinfo%output,*)'COUNT(ExpOccAOS)',COUNT(ExpOccAOS)
          do I=1,nocc
             if(DECinfo%frozencore.AND.I.LE.ncore)THEN
              IF(ExpOccAOS(I))THEN                      
                WRITE(DECinfo%output,'(A,I5)')'Orbital i - core orbital (set to 400000 + i)  i=',I
                BruteForceOccContribs(I) = 400000.0E0_realk + I !core orbital                      
               ELSE
                WRITE(DECinfo%output,'(A,I5)')'Orbital i - core orbital and not included in AOS space i=',I
                BruteForceOccContribs(I) = 0.0E0_realk
              ENDIF
             ELSEIF(ExpOccAOS(I))THEN
                VirtAOS_new = ExpVirtAOS
                OccAOS_new = ExpOccAOS
                WRITE(DECinfo%output,*)'COUNT(OccAOS_new)',COUNT(OccAOS_new)
                WRITE(DECinfo%output,*)'REMOVE number ',I
                OccAOS_new(I) = .FALSE.
                WRITE(DECinfo%output,*)'After Removal COUNT(OccAOS_new)',COUNT(OccAOS_new)
                call SanityCheckOrbAOS(AtomicFragment,nocc,'O',OccAOS_new)
                WRITE(DECinfo%output,*)'After SanityChack COUNT(OccAOS_new)',COUNT(OccAOS_new)
                IF(COUNT(OccAOS_new).NE.COUNT(ExpOccAOS))THEN
                   bin_reduction_converged = .FALSE.
                   call atomic_fragment_free(AtomicFragment)
                   call atomic_fragment_init_orbital_specific(MyAtom,nunocc, nocc, VirtAOS_new, &
                        & OccAOS_new,OccOrbitals,UnoccOrbitals,MyMolecule,mylsitem,AtomicFragment,.true.,.false.)
                   call fragment_energy_and_prop(AtomicFragment)
                   WRITE(DECinfo%output,'(A,I5)')'Orbital i included in AOS space (CALC)  i=',I
                   BruteForceOccContribs(I) = ABS(AtomicFragment%EoccFOP-OccEnergyOld)
                ELSE
                   WRITE(DECinfo%output,'(A,I5)')'Orbital i included in EOS space (set to 300000 + i)  i=',I
                   BruteForceOccContribs(I) = 300000.0E0_realk + I !EOS orbital
                ENDIF
             ELSE
                WRITE(DECinfo%output,'(A,I5)')'Orbital i not included in AOS space from Expansion (set to 0) i=',I
                BruteForceOccContribs(I) = 0.0E0_realk
             ENDIF
          ENDDO
          do I=1,nocc
             WRITE(DECinfo%output,'(A,I5,A,ES16.8)')'Non Sorted BruteForceOccContribs(',I,')',BruteForceOccContribs(I)
          enddo
          call real_inv_sort_with_tracking(BruteForceOccContribs,BruteForceTrackListOcc,nocc)
          do I=1,nocc
             II = BruteForceTrackListOcc(I)
             WRITE(DECinfo%output,'(A,I5,A,ES16.8)')'Sorted BruteForceOccContribs(',II,')',BruteForceOccContribs(I)
          enddo
          do I=1,nocc
             II = OrbOccDistTrackMyAtom(I)
             WRITE(DECinfo%output,'(A,I5,A,ES16.8)')'Sorted Distance(',II,')',SortedDistanceTableOrbAtomOcc(I)
          enddo
!          do I=1,nocc
!             II = OrbOccFockTrackMyAtom(I)
!             WRITE(DECinfo%output,'(A,I5,A,ES16.8)')'Sorted FockMatrix(',II,')',FmaxOcc(I)
!          enddo

!       ENDIF

       call mem_alloc(BruteForceVirtContribs,nunocc)
       BruteForceVirtContribs = 0.0E0_realk
       call mem_alloc(BruteForceTrackListVirt,nunocc)       
       do I=1,nunocc
          BruteForceTrackListVirt(I) = I
       enddo
       IF(DECinfo%OnlyVirtPart)THEN
          call lsquit('Brute force and virtonly not implemented',-1)
       ENDIF
!!$!       IF(DECinfo%OnlyVirtPart)THEN
!!$          !determine Brute Force Occupied Energy contribution
!!$          WRITE(DECinfo%output,*)'COUNT(ExpVirtAOS)',COUNT(ExpVirtAOS)
!!$          do I=1,COUNT(ExpVirtAOS)
!!$             VirtAOS_new = ExpVirtAOS
!!$             OccAOS_new = ExpOccAOS
!!$             WRITE(DECinfo%output,*)'COUNT(VirtAOS_new)',COUNT(VirtAOS_new)
!!$             WRITE(DECinfo%output,*)'REMOVE number ',I
!!$             K=0
!!$             DO J=1,SIZE(VirtAOS_new)
!!$                IF(VirtAOS_new(J))K = K + 1
!!$                IF(K.EQ.I)THEN
!!$                   VirtAOS_new(J) = .FALSE.
!!$                   EXIT
!!$                ENDIF
!!$             ENDDO
!!$             WRITE(DECinfo%output,*)'After Removal COUNT(VirtAOS_new)',COUNT(VirtAOS_new)
!!$
!!$             bin_reduction_converged = .FALSE.
!!$             call atomic_fragment_free(AtomicFragment)
!!$             call atomic_fragment_init_orbital_specific(MyAtom,nunocc, nocc, VirtAOS_new, &
!!$                  & OccAOS_new,OccOrbitals,UnoccOrbitals,MyMolecule,mylsitem,AtomicFragment,.true.,.false.)
!!$             call fragment_energy_and_prop(AtomicFragment)
!!$             BruteForceOccContribs(I) = ABS(AtomicFragment%EoccFOP)            
!!$          ENDDO
!!$          call real_inv_sort_with_tracking(BruteForceOccContribs,BruteForceTrackListOcc,nocc)
!!$!       ENDIF
    ENDIF
    ! Set FOT and convergence control
    FOT = DECinfo%FOT
    reduction_converged=.false.

    IF(BinarySearch)THEN
       call mem_alloc(SortedOccContribs,nocc)
       call mem_alloc(TrackSortedOccContribs,nocc)
       do i=1,nocc
          TrackSortedOccContribs(i) = i          
       enddo
       SortedOccContribs = ABS(OccContribs)
       call real_inv_sort_with_tracking(SortedOccContribs,TrackSortedOccContribs,nocc)

!       IF(DECinfo%PL.GT.1)THEN
          IF(BruteForce.OR.FockMatrixOrdering)THEN
             IF(FockMatrixOrdering)THEN
                do I=1,nocc
                   II = OrbOccFockTrackMyAtom(I)
                   WRITE(DECinfo%output,'(A,I5,A,ES16.8)')'Sorted OccFockMatrix(',II,')',FmaxOcc(I)
                enddo
             ENDIF
             do I=1,nocc
                II = TrackSortedOccContribs(I)
                WRITE(DECinfo%output,'(A,I3,A,ES16.8)')'Sorted OccContribs(',II,')',SortedOccContribs(I)
             enddo
             do I=1,nocc
                II = OrbOccDistTrackMyAtom(I)
                WRITE(DECinfo%output,'(A,I5,A,ES16.8)')'Sorted Distance(',II,')',SortedDistanceTableOrbAtomOcc(I)
             enddo
          ENDIF
!       ENDIF

       call mem_alloc(SortedVirtContribs,nunocc)
       call mem_alloc(TrackSortedVirtContribs,nunocc)
       do i=1,nunocc
          TrackSortedVirtContribs(i) = i          
       enddo
       SortedVirtContribs = ABS(VirtContribs)
       call real_inv_sort_with_tracking(SortedVirtContribs,TrackSortedVirtContribs,nunocc)

!       call BinarySearch(PriorityListOcc,PriorityListVirt,...)

       IF(OrbDistanceSpec)THEN
          !the highest number of Occupied that has been verified to be acceptable
          nHigherOcc  = COUNT(ExpOccAOS)  ! (not nocc_orig - it does not include frozencore orbs)
       ELSE
          nHigherOcc  = nocc_orig         
       ENDIF
       !the highest number of Virtuals that has been verified to be acceptable
       nHigherVirt = nvirt_orig 

       IF(OrbDistanceSpec)THEN
          !In the LowerLimit we include the first N orbitals (sorted by Contribs)
          !which was included in the second last fragment 
          nLowerOcc = 0 
          IF(.NOT.BruteForce)THEN
             nLowerOcc = noccEOS
             !DO I=1,nocc
             !   II = TrackSortedOccContribs(I)
             !   IF(ExpOldOccAOS(II))THEN
             !      nLowerOcc = nLowerOcc + 1
             !   ELSE
             !      EXIT
             !   ENDIF
             !ENDDO
          ELSE
             nLowerOcc = 1 
          ENDIF
          !nLowerVirt = 0
          !DO I=1,nunocc
          !   II = TrackSortedVirtContribs(I)
          !   IF(ExpOldVirtAOS(II))THEN
          !      nLowerVirt = nLowerVirt + 1
          !   ELSE
          !      EXIT
          !   ENDIF
          !ENDDO
          nLowerVirt = nunoccEOS
       ELSE
          !RejectThresh a little bigger than FOT to ensure that this step would be rejected 
          !and provide a good lower guess
          RejectThresh = FOT*1.01E0_realk 
          OccAOS_new = OccAOS_orig
          VirtAOS_new = VirtAOS_orig
          !Reduce occupied/virtual AOS according to rejection threshold
          write(DECinfo%output,*) ' FOP OCC: '
          call ReduceSpace_orbitalspecific(AtomicFragment,nocc,OccContribs,'O',&
               & RejectThresh,OccAOS_new,nocc_new)
          write(DECinfo%output,*) ' FOP VIRT: '
          call ReduceSpace_orbitalspecific(AtomicFragment,nunocc,VirtContribs,'V',&
               & RejectThresh,VirtAOS_new,nvirt_new)
          nLowerOcc  = nocc_new  !highest number of Occupied that has been verified to be NOT acceptable
          nLowerVirt = nvirt_new !highest number of Virtuals that has been verified to be NOT acceptable  
       ENDIF
       !nLowerOcc = noccEOS
       !nLowerVirt = nunoccEOS
       call mem_dealloc(SortedOccContribs)
       call mem_dealloc(SortedVirtContribs)

       reduction_converged = .FALSE.
       bin_reduction_converged = .FALSE.

       nDimOcc  = nHigherOcc-nLowerOcc
       nDimVirt = nHigherVirt-nLowerVirt
       
       !IF(O2V2choice)THEN
       !   nvirt_new = nHigherVirt - (nHigherVirt-nLowerVirt)/2
       !   nocc_new  = nHigherOcc -  (nHigherOcc-nLowerOcc)/2
       !   nO2V2_ModOcc = (nocc_new*nocc_new*i8)*(nHigherVirt*nHigherVirt*i8)
       !   nO2V2_ModVirt = (nvirt_new*nvirt_new*i8)*(nHigherVirt*nHigherVirt*i8)
       !   ModVirt = nO2V2_ModOcc.GT.nO2V2_ModVirt
       !ELSE
       !   ModVirt = nDimVirt.GE.nDimOcc 
       !ENDIF
       ModVirt = .true.

       IF(ModVirt)THEN
          nvirt_old = nLowerVirt
          nocc_old  = nHigherOcc
       ELSE
          nvirt_old = nHigherVirt
          nocc_old  = nLowerOcc
       ENDIF
       bin_reduction_converged = .FALSE.
       bin_virt_conv = .FALSE.
       bin_occ_conv = .FALSE.
       BINARY_REDUCTION_LOOP: do iter=1,max_iter_red
          nDimOcc  = nHigherOcc-nLowerOcc     !diff between converged space and nonconverged space
          nDimVirt = nHigherVirt-nLowerVirt
          !determine if we should do a step in Virtual or Occupied space
          !IF(O2V2choice)THEN
          !   !take the step that have biggest potential to reduce the dim O**2*V**2
          !   IF(bin_reduction_converged)THEN
          !      nvirt_new = nLowerVirt + (nHigherVirt-nLowerVirt)/2
          !      nocc_new  = nLowerOcc  + (nHigherOcc-nLowerOcc)/2
          !   ELSE
          !      nvirt_new = nHigherVirt - (nHigherVirt-nLowerVirt)/2
          !      nocc_new  = nHigherOcc -  (nHigherOcc-nLowerOcc)/2
          !   ENDIF
          !   nO2V2_ModOcc = (nocc_new*nocc_new*i8)*(nHigherVirt*nHigherVirt*i8)
          !   nO2V2_ModVirt = (nvirt_new*nvirt_new*i8)*(nHigherVirt*nHigherVirt*i8)
          !   ModVirt = nO2V2_ModOcc.GT.nO2V2_ModVirt
          !ELSEIF(DistanceRemoval.AND.DECinfo%onlyOccPart)THEN
          !   !First converge Virtual space - then remove Occupied 
          !   ModVirt = .TRUE.
          !   IF(bin_virt_conv)ModVirt = .FALSE.
          !ELSEIF(DistanceRemoval.AND.DECinfo%onlyVirtPart)THEN
          !   !First converge Occupied space - then remove Vitual 
          !   ModVirt = .FALSE.
          !   IF(bin_occ_conv)ModVirt = .TRUE.
          !ELSE
          !   !take the step that have biggest potential to 
          !   !reduce the number of orbitals Occ + Virt
          !   ModVirt = nDimVirt.GE.nDimOcc 
          !ENDIF

          ! Start reducing occ space only when virtual space is fully converged
          ModVirt = .true.
          if (bin_virt_conv) ModVirt = .false.

          IF(TestOcc)THEN
             ModVirt = .FALSE.
             IF(bin_occ_conv)ModVirt = .TRUE.
          ENDIF          
          IF(DECinfo%PL.GT.1)THEN
             IF(ModVirt)THEN
                WRITE(DECinfo%output,*)'BIN SEARCH Modify Virtual Orbital Space'
             ELSE
                WRITE(DECinfo%output,*)'BIN SEARCH Modify Occupied Orbital Space'
             ENDIF
          ENDIF
          nocc_new = nHigherOcc
          nvirt_new = nHigherVirt
          IF(ModVirt)THEN
             nvirt_new = nLowerVirt + (nHigherVirt-nLowerVirt)/2
          ELSE
             nocc_new  = nLowerOcc  + (nHigherOcc-nLowerOcc)/2
          ENDIF
          IF(DECinfo%PL.GT.1)THEN
             WRITE(DECinfo%output,*)'BIN SEARCH nHigherOcc,nHigherVirt',nHigherOcc,nHigherVirt
             WRITE(DECinfo%output,*)'BIN SEARCH nLowerOcc,nLowerVirt  ',nLowerOcc,nLowerVirt
             WRITE(DECinfo%output,*)'BIN SEARCH nocc_new,nvirt_new    ',nocc_new,nvirt_new
          ENDIF
          IF(nLowerOcc.GT.nHigherOcc)CALL LSQUIT('nLowerOcc.GT.nHigherOcc',-1)
          IF(nLowerVirt.GT.nHigherVirt)CALL LSQUIT('nLowerVirt.GT.nHigherVirt',-1)
          IF(BruteForce)THEN
             IF(DECinfo%onlyOccPart)THEN
                call ReduceSpace_binary(AtomicFragment,nocc,TrackSortedOccContribs,BruteForce,&
                     &BruteForceTrackListOcc,'O',OccAOS_new,nocc_new)
             ELSE
                !use OccContribs/distance to do it
                call ReduceSpace_binary(AtomicFragment,nocc,TrackSortedOccContribs,DistanceRemoval,&
                     &OrbOccDistTrackMyAtom,'O',OccAOS_new,nocc_new)
             ENDIF
!             IF(DECinfo%onlyVirtPart)THEN
!                call ReduceSpace_binary(AtomicFragment,nunocc,TrackSortedVirtContribs,BruteForce,&
!                     &BruteForceTrackListVirt,'V',VirtAOS_new,nvirt_new)
!             ELSE
                call ReduceSpace_binary(AtomicFragment,nunocc,TrackSortedVirtContribs,DistanceRemoval,&
                     &OrbVirtDistTrackMyAtom,'V',VirtAOS_new,nvirt_new)
!             ENDIF
          ELSEIF(FockMatrixOrdering)THEN
             IF(DECinfo%onlyOccPart)THEN
                !use Fock matrix ordering to remove Occupied
                call ReduceSpace_binary(AtomicFragment,nocc,TrackSortedOccContribs,FockMatrixOrdering,&
                     & OrbOccFockTrackMyAtom,'O',OccAOS_new,nocc_new)
             ELSE
                !use OccContribs/distance to do it
                call ReduceSpace_binary(AtomicFragment,nocc,TrackSortedOccContribs,DistanceRemoval,&
                     &OrbOccDistTrackMyAtom,'O',OccAOS_new,nocc_new)
             ENDIF
             IF(DECinfo%onlyVirtPart)THEN
                !use Fock matrix ordering to remove Virtual 
                call ReduceSpace_binary(AtomicFragment,nunocc,TrackSortedVirtContribs,FockMatrixOrdering,&
                     &OrbVirtFockTrackMyAtom,'V',VirtAOS_new,nvirt_new)
             ELSE
                call ReduceSpace_binary(AtomicFragment,nunocc,TrackSortedVirtContribs,DistanceRemoval,&
                     &OrbVirtDistTrackMyAtom,'V',VirtAOS_new,nvirt_new)
             ENDIF
          ELSE
             call ReduceSpace_binary(AtomicFragment,nocc,TrackSortedOccContribs,DistanceRemoval,&
                  &OrbOccDistTrackMyAtom,'O',OccAOS_new,nocc_new)
             call ReduceSpace_binary(AtomicFragment,nunocc,TrackSortedVirtContribs,DistanceRemoval,&
                  &OrbVirtDistTrackMyAtom,'V',VirtAOS_new,nvirt_new)
          ENDIF

          IF(DECinfo%PL.GT.1)THEN
             WRITE(DECinfo%output,*)'BIN SEARCH COUNT(OccAOS_new)',COUNT(OccAOS_new)
             WRITE(DECinfo%output,*)'BIN SEARCH COUNT(VirtAOS_new)',COUNT(VirtAOS_new)
          ENDIF

          ! Get new atomic fragment
          ! ***************************
          bin_reduction_converged = .FALSE.
          call atomic_fragment_free(AtomicFragment)
          call atomic_fragment_init_orbital_specific(MyAtom,nunocc, nocc, VirtAOS_new, &
               & OccAOS_new,OccOrbitals,UnoccOrbitals,MyMolecule,mylsitem,AtomicFragment,.true.,.false.)
          call fragment_energy_and_prop(AtomicFragment)
          
          ! Check if reduced fragment energy is converged to FOT precision
          ! **************************************************************
          LagEnergyDiff=abs(AtomicFragment%LagFOP-LagEnergyOld)
          OccEnergyDiff=abs(AtomicFragment%EoccFOP-OccEnergyOld)
          VirtEnergyDiff=abs(AtomicFragment%EvirtFOP-VirtEnergyOld)
          call fragopt_print_info(AtomicFragment,LagEnergyDiff,OccEnergyDiff,VirtEnergyDiff,iter)
          call fragopt_check_convergence(LagEnergyDiff,OccEnergyDiff,VirtEnergyDiff,&
               &FOT,bin_reduction_converged)
          IF(bin_reduction_converged)THEN
             !Save information for later
             LagEnergy = AtomicFragment%LagFOP
             OccEnergy = AtomicFragment%EoccFOP
             VirtEnergy = AtomicFragment%EvirtFOP
             OccAOS_old = OccAOS_new
             VirtAOS_old = VirtAOS_new
          ENDIF

          if(ModVirt)THEN
             IF(ABS(nvirt_old-nvirt_new).LE.1)THEN
                bin_virt_conv = .TRUE.
                IF(bin_reduction_converged)THEN
                   nvirt_old  = nvirt_new
                   nLowerVirt = nvirt_new
                   nHigherVirt= nvirt_new
                ELSE
                   nvirt_old = nHigherVirt
                   nvirt_new = nHigherVirt
                   nLowerVirt= nHigherVirt   
                ENDIF
                IF(DECinfo%PL.GT.1)WRITE(DECinfo%output,*)'VIRTUAL REDUCTION CONVERGED'
                !IF(bin_occ_conv)reduction_converged = .TRUE.
             ENDIF
             !IF((ABS(nocc_old-nocc_new).LE.1).AND.((nHigherOcc-nLowerOcc).LE.1))THEN
             !   bin_occ_conv = .TRUE. !somehow the occupied also converged
             !   IF(bin_reduction_converged)THEN
             !      nOcc_old  = nOcc_new
             !      nLowerOcc = nocc_new
             !      nHigherOcc= nocc_new
             !   ELSE
             !      nocc_old = nHigherOcc
             !      nocc_new = nHigherOcc
             !      nLowerOcc= nHigherOcc
             !   ENDIF
             !   reduction_converged = .TRUE.
             !ENDIF                
          else
             IF(ABS(nocc_old-nocc_new).LE.1)THEN
                bin_occ_conv = .TRUE.
                IF(bin_reduction_converged)THEN
                   nOcc_old  = nOcc_new
                   nLowerOcc = nocc_new
                   nHigherOcc= nocc_new
                ELSE
                   nocc_old = nHigherOcc
                   nocc_new = nHigherOcc
                   nLowerOcc= nHigherOcc
                ENDIF
                IF(DECinfo%PL.GT.1)WRITE(DECinfo%output,*)'OCCUPIED REDUCTION CONVERGED'
                reduction_converged = .TRUE.
             ENDIF
             !IF((ABS(nvirt_old-nvirt_new).LE.1).AND.((nHigherVirt-nLowerVirt).LE.1))THEN
             !   bin_virt_conv = .TRUE. !somehow the virtual also converged
             !   IF(bin_reduction_converged)THEN
             !      nvirt_old  = nvirt_new
             !      nLowerVirt = nvirt_new
             !      nHigherVirt= nvirt_new
             !   ELSE
             !      nvirt_old = nHigherVirt
             !      nvirt_new = nHigherVirt
             !      nLowerVirt= nHigherVirt   
             !   ENDIF
             !   reduction_converged = .TRUE.
             !ENDIF
          endif

          IF(reduction_converged)THEN
             IF(DECinfo%PL.GT.1)WRITE(DECinfo%output,*)'reduction_converged'
             IF(.not.bin_reduction_converged)THEN
                !the last binary search step did not succeed so  
                !we init everything to the last successful step
                call atomic_fragment_free(AtomicFragment)
                call atomic_fragment_init_orbital_specific(MyAtom,nunocc, nocc, VirtAOS_old, &
                     & OccAOS_old,OccOrbitals,UnoccOrbitals,MyMolecule,mylsitem,AtomicFragment,.true.,.false.)
                AtomicFragment%LagFOP = LagEnergy
                AtomicFragment%EoccFOP = OccEnergy
                AtomicFragment%EvirtFOP = VirtEnergy
             ENDIF
             EXIT BINARY_REDUCTION_LOOP
          ENDIF

          nocc_old = nocc_new
          nvirt_old = nvirt_new
          IF(bin_reduction_converged)THEN
             IF(DECinfo%PL.GT.1)WRITE(DECinfo%output,*)'Step accepted'
             !the current number of occ and virt is good or too high
             IF(ModVirt)THEN
                nHigherVirt = nvirt_new    !New Upper limit 
             ELSE
                nHigherOcc = nocc_new      !New Upper limit 
             ENDIF
          ELSE
             IF(DECinfo%PL.GT.1)WRITE(DECinfo%output,*)'Step not accepted'
             !the current number of occ and virt is NOT good enough 
             IF(ModVirt)THEN
                nLowerVirt = nvirt_new     !New Lower limit 
             ELSE
                nLowerOcc = nocc_new       !New Lower limit 
             ENDIF
          ENDIF
       END DO BINARY_REDUCTION_LOOP
       call mem_dealloc(TrackSortedOccContribs)
       call mem_dealloc(TrackSortedVirtContribs)

       IF(BruteForce)THEN
          call mem_dealloc(BruteForceOccContribs)
          call mem_dealloc(BruteForceTrackListOcc)
          call mem_dealloc(BruteForceVirtContribs)
          call mem_dealloc(BruteForceTrackListVirt)
       ENDIF
    ELSE
       REDUCTION_LOOP: do iter=1,max_iter_red
          
          write(DECinfo%output,*) 'FOP Starting reduction step ', iter
          if(iter == 1) then
             RejectThresh = FOT
          else  ! decrease rejection threshold by a factor 5 in each step
             RejectThresh = RejectThresh/5.0E0_realk
          end if
          write(DECinfo%output,'(a,ES13.5)') 'FOP Rejection threshold  : ',RejectThresh
          
          ! Number of occupied and virtual AOS orbitals in current fragment
          nocc_old = AtomicFragment%noccAOS
          nvirt_old = AtomicFragment%nunoccAOS
          
          
          ! Set new logical vectors controlling which occupied and virtual AOS orbitals are included
          ! ****************************************************************************************
          
          ! 1. Copy information from optimal fragment from expansion loop
          OccAOS_new = OccAOS_old
          VirtAOS_new = VirtAOS_old
          
          ! 2. Reduce occupied/virtual AOS according to rejection threshold
          write(DECinfo%output,*) ' FOP OCC: '
          call ReduceSpace_orbitalspecific(AtomicFragment,nocc,OccContribs,'O',&
               & RejectThresh,OccAOS_new,nocc_new)
          write(DECinfo%output,*) ' FOP VIRT: '
          call ReduceSpace_orbitalspecific(AtomicFragment,nunocc,VirtContribs,'V',&
               & RejectThresh,VirtAOS_new,nvirt_new)
          
          
          ! Check if the orbital space was reduced
          ! **************************************
          
          ! Special case: No reduction has occurred in any step
          if(iter == max_iter_red .and. (nocc_new==nocc_orig) .and. (nvirt_new==nvirt_orig) ) then
             
             write(DECinfo%output,*) "FOP No reduction possible. Use original converged fragment"
             
             ! Go back to old fragment
             ! ***********************
             call atomic_fragment_free(AtomicFragment)
             call atomic_fragment_init_orbital_specific(MyAtom,nunocc, nocc, VirtAOS_orig, &
                  & OccAOS_orig,OccOrbitals,UnoccOrbitals,MyMolecule,mylsitem,&
                  & AtomicFragment,.true.,.false.)
             
             call fragment_energy_and_prop(AtomicFragment)

             ! we set the reduction to be converged to avoid quiting:
             reduction_converged = .true.
             exit REDUCTION_LOOP
          end if
          
          ! Cycle if no atoms were rejected by the rejection procedure
          if ((nocc_old == nocc_new) .and. (nvirt_old == nvirt_new)) then
             write(DECinfo%output,*) 'FOP No reduction occurred - we try to reduce again'
             cycle
          end if
          
          
          ! Get reduced atomic fragment
          ! ***************************
          call atomic_fragment_free(AtomicFragment)
          call atomic_fragment_init_orbital_specific(MyAtom,nunocc, nocc, VirtAOS_new, &
               & OccAOS_new,OccOrbitals,UnoccOrbitals,MyMolecule,mylsitem,AtomicFragment,.true.,.false.)
          call fragment_energy_and_prop(AtomicFragment)
          
          
          ! Check if reduced fragment energy is converged to FOT precision
          ! **************************************************************
          LagEnergyDiff=abs(AtomicFragment%LagFOP-LagEnergyOld)
          OccEnergyDiff=abs(AtomicFragment%EoccFOP-OccEnergyOld)
          VirtEnergyDiff=abs(AtomicFragment%EvirtFOP-VirtEnergyOld)
          call fragopt_print_info(AtomicFragment,LagEnergyDiff,OccEnergyDiff,VirtEnergyDiff,iter)
          call fragopt_check_convergence(LagEnergyDiff,OccEnergyDiff,VirtEnergyDiff,&
               &FOT,reduction_converged)
          
          ! Quit if we are converged
          if(reduction_converged) then
             write(DECinfo%output,*) 'FOP Reduction of fragment converged in step', iter
             exit REDUCTION_LOOP
          end if
       end do REDUCTION_LOOP
    ENDIF

    ! CHECK THAT REDUCTION CONVERGED:
    if (.not.reduction_converged) then
      write(DECinfo%output,*) 'FOP'
      write(DECinfo%output,*) "WARNING:: Reduction not converged for fragment",MyAtom
      write(DECinfo%output,*) "          Exit loop after",iter,"iterations"
      call lsquit('Fragment reduction not converged',DECinfo%output)
    end if

    ! Update t1 amplitudes for full molecule
    ! **************************************
    ! KK: These are leftovers for treating long-range t1 interactions.
    !     Some it may be reused but it general the lone-range t1 implementation requires 
    !     rethinking the basic strategy!

    if(DECinfo%SinglesPolari) then

       ! Extract (virt AOS,occ EOS) indices from fragment
       ! (This is necessary to avoid double counting)
       call extract_specific_fragmentt1(AtomicFragment,.true.,.false.)

       ! Update
       call update_full_t1_from_atomic_frag(AtomicFragment,t1full)

       ! Delete t1 information in fragment
       call free_fragment_t1(AtomicFragment)

    end if



    ! Print out info
    ! **************
    write(DECinfo%output,*)'FOP'
    write(DECinfo%output,'(1X,a)') 'FOP========================================================='
    write(DECinfo%output,'(1X,a,i4)') 'FOP    LOCAL REDUCTION HAS CONVERGED FOR SITE',MyAtom
    write(DECinfo%output,'(1X,a)') 'FOP---------------------------------------------------------'
    write(DECinfo%output,'(1X,a,i4)')    'FOP Done: Fragment number                  :', MyAtom
    write(DECinfo%output,'(1X,a,i4)')    'FOP Done: Number of orbitals in virt total :', &
         & AtomicFragment%nunoccAOS
    write(DECinfo%output,'(1X,a,i4)')    'FOP Done: Number of orbitals in occ total  :', &
         & AtomicFragment%noccAOS
    write(DECinfo%output,'(1X,a,i4)')     'FOP Done: Number of basis functions        :', &
         & AtomicFragment%nbasis
    if(.not. DECinfo%onlyvirtpart) then
       write(DECinfo%output,'(1X,a,f16.10)') 'FOP Done: Occupied Fragment energy         :', &
            & AtomicFragment%EoccFOP
    endif
    if(.not. DECinfo%onlyoccpart) then
       write(DECinfo%output,'(1X,a,f16.10)') 'FOP Done: Virtual Fragment energy          :', &
            & AtomicFragment%EvirtFOP
    endif
    if(.NOT.(DECinfo%onlyoccpart.OR. DECinfo%onlyvirtpart))then
       write(DECinfo%output,'(1X,a,f16.10)') 'FOP Done: Lagrangian Fragment energy       :', &
            & AtomicFragment%LagFOP
    end if
    IF(.NOT.BinarySearch)THEN
       write(DECinfo%output,'(1X,a,g14.2)')  'FOP Done: Reduction threshold              :', RejectThresh
    ENDIF
    write(DECinfo%output,'(1X,a,i6,a,i6,a,f5.2,a)')  'FOP Done: Occupied reduction removed ', &
         & nocc_orig-AtomicFragment%noccAOS, ' of ', nocc_orig, ' orbitals ( ', &
         & (nocc_orig-AtomicFragment%noccAOS)*100.0_realk/nocc_orig, ' %)'
    write(DECinfo%output,'(1X,a,i6,a,i6,a,f5.2,a)')  'FOP Done: Virtual  reduction removed ', &
         & nvirt_orig-AtomicFragment%nunoccAOS, ' of ', nvirt_orig, ' orbitals ( ', &
         & (nvirt_orig-AtomicFragment%nunoccAOS)*100.0_realk/nvirt_orig, ' %)'
    write(DECinfo%output,'(1X,a,/)') 'FOP========================================================='
    write(DECinfo%output,*) 'FOP'


    call mem_dealloc(OccAOS_orig)
    call mem_dealloc(VirtAOS_orig)
    call mem_dealloc(OccAOS_old)
    call mem_dealloc(VirtAOS_old)
    call mem_dealloc(OccAOS_new)
    call mem_dealloc(VirtAOS_new)

  end subroutine fragopt_reduce_local_orbitals

  subroutine GetFockMaxElement(Fmax,nocc,OccFock,noccEOS,occEOSidx,ncore)
    implicit none
    integer,intent(in)  :: nocc,noccEOS,ncore
    integer,intent(in)  :: occEOSidx(noccEOS)
    real(realk),intent(in) :: OccFock(nocc,nocc) 
    real(realk),intent(inout) :: Fmax(nocc)
    !
    integer :: i,j,ii
    do j=1,nocc
       Fmax(j) = 0.0E0_realk
       !loop over EOS orbitals 
       do i=1,noccEOS
          ii = occEOSidx(i)
          Fmax(j) = MAX(Fmax(j),ABS(OccFock(ii,j)))          
       end do
    enddo
    !set the EOS orbitals to an artifical high number 
    !so that they are always included 
    do i=1,noccEOS
       ii = occEOSidx(i)
       Fmax(ii) = 300000.0E0_realk + ii !EOS orbital                      
    enddo
    !set the core orbitals to an artifical high number 
    !so that they are always included (removed when builing fragment) 
    do i=1,ncore
       Fmax(i) = 400000.0E0_realk + I !core orbital                      
    enddo
  end subroutine GetFockMaxElement

  !> \brief Check if fragment energy (or energies) is converged to FOT precision
  !> \author Kasper Kristensen
  !> \date September 2013
  subroutine fragopt_check_convergence(LagEnergyDiff,OccEnergyDiff,VirtEnergyDiff,FOT,converged)
    implicit none
    !> ABSOLUTE Energy differences for Lagrangian, occupied, and virtual partitioning schemes
    !> (difference between current and previous iteration)
    real(realk),intent(in) :: LagEnergyDiff,OccEnergyDiff,VirtEnergyDiff
    !> Fragment optimization threshold
    real(realk),intent(in) :: FOT
    !> Is fragment energy (or energies) converged?
    logical,intent(inout) :: converged
    logical :: lag_converged, occ_converged, virt_converged

    ! Default: Test convergence for both Lagrangian, occupied, and virtual energies
    ! If DECinfo%onlyoccpart: Only test for energy using occupied partitioning scheme
    ! If DECinfo%onlyvirtpart: Only test for energy using virtual partitioning scheme


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
       if  (LagEnergyDiff < FOT) then
          write(DECinfo%output,'(1X,a,F14.9)') 'FOP: Lagrangian energy converged, energydiff =', &
               & LagEnergyDiff
          lag_converged=.true.
       else
          write(DECinfo%output,*) 'FOP: Lagrangian energy NOT converged'
          lag_converged=.false.
       end if
    end if TEST_CONVERGENCE_LAG

    ! Occupied 
    TEST_CONVERGENCE_OCC: if(DECinfo%OnlyVirtPart) then
       occ_converged=.true.
    else
       if  (OccEnergyDiff < FOT) then
          write(DECinfo%output,'(1X,a,F14.9)') 'FOP: Occupied energy converged, energydiff   =', OccEnergyDiff
          occ_converged=.true.
       else
          write(DECinfo%output,*) 'FOP: Occupied energy NOT converged'
          occ_converged=.false.
       end if
    endif TEST_CONVERGENCE_OCC

    ! Virtual 
    TEST_CONVERGENCE_VIRT: if(DECinfo%OnlyOccPart) then
       virt_converged=.true.
    else
       if  (VirtEnergyDiff < FOT) then
          write(DECinfo%output,'(1X,a,F14.9)') 'FOP: Virtual energy converged, energydiff    =', VirtEnergyDiff
          virt_converged=.true.
       else
          write(DECinfo%output,*) 'FOP: Virtual energy NOT converged'
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



  !> Fragment optimization, special case where we want to include full molecular system
  ! in fragment orbital space (only for debugging purposes).
  !> \date September 2013
  !> \author Kasper Kristensen
  subroutine fragopt_include_fullmolecule(MyAtom,AtomicFragment, &
       &OccOrbitals,nOcc,UnoccOrbitals,nUnocc, &
       &MyMolecule,mylsitem,freebasisinfo,t1full)
    implicit none
    !> Number of occupied orbitals in molecule
    integer, intent(in) :: nOcc
    !> Number of unoccupied orbitals in molecule
    integer, intent(in) :: nunocc
    !> Central atom in molecule
    integer, intent(in) :: MyAtom
    !> Atomic fragment to be optimized
    type(decfrag),intent(inout)        :: AtomicFragment
    !> All occupied orbitals
    type(decorbital), dimension(nOcc), intent(in)      :: OccOrbitals
    !> All unoccupied orbitals
    type(decorbital), dimension(nUnocc), intent(in)    :: UnoccOrbitals
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
    call fragment_init_simulate_full(MyAtom,nunocc, nocc, OccOrbitals,UnoccOrbitals,&
         & MyMolecule,mylsitem,AtomicFragment,.true.)

    ! Set information for fragment-adapted orbitals
    FullFragAdapt: if(DECinfo%fragadapt) then

       ! Set model to be MP2
       MyMolecule%ccmodel(MyAtom,MyAtom) = MODEL_MP2

       ! Integrals (ai|bj)
       !call get_VOVO_integrals(AtomicFragment%mylsitem,AtomicFragment%nbasis,&
       !     & AtomicFragment%noccAOS,AtomicFragment%nunoccAOS,&
       !     & AtomicFragment%Cv, AtomicFragment%Co, g)
       ! MP2 amplitudes
       !call mp2_solver(AtomicFragment%noccAOS,AtomicFragment%nunoccAOS,&
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

  end subroutine fragopt_include_fullmolecule


  !> \brief prints info for atomic fragment in a given iteration of the fragment optimization
  !> using the Lagrangian partitioning scheme.
  !> \date: august-2011
  !> \author: Ida-Marie Hoeyvik
  subroutine fragopt_print_info(Fragment,LagEnergyDiff,OccEnergyDiff,VirtEnergyDiff,iter)
    implicit none
    !> Atomic fragment
    type(decfrag),intent(inout) :: Fragment
    !> Lagrangian energy difference since last loop in fragment optimization
    real(realk),intent(in) :: LagEnergyDiff
    !> Occupied energy difference since last loop in fragment optimization
    real(realk),intent(in) :: VirtEnergyDiff
    !> Virtual energy difference since last loop in fragment optimization
    real(realk),intent(in) :: OccEnergyDiff
    !> Iteration number
    integer,intent(in) :: iter


    write(DECinfo%output,*)'FOP'
    write(DECinfo%output,'(1X,a)') 'FOP========================================================='
    write(DECinfo%output,'(1X,a,i4)') 'FOP              Fragment information, loop', iter
    write(DECinfo%output,'(1X,a)') 'FOP---------------------------------------------------------'
    write(DECinfo%output,'(1X,a,i4)')    'FOP Loop: Fragment number                  :', &
         & fragment%EOSatoms(1)
    write(DECinfo%output,'(1X,a,i4)')    'FOP Loop: Number of orbitals in virt total :', &
         & Fragment%nunoccAOS
    write(DECinfo%output,'(1X,a,i4)')    'FOP Loop: Number of orbitals in occ total  :', &
         & Fragment%noccAOS
    write(DECinfo%output,'(1X,a,i4)')    'FOP Loop: Number of basis functions        :', &
         & Fragment%nbasis
    if(.not. DECinfo%OnlyVirtPart) then
       write(DECinfo%output,'(1X,a,f16.10)') 'FOP Loop: Occupied Fragment energy         :', &
            & Fragment%EoccFOP
    endif
    if(.not. DECinfo%OnlyOccPart) then
       write(DECinfo%output,'(1X,a,f16.10)') 'FOP Loop: Virtual Fragment energy          :', &
            & Fragment%EvirtFOP
    endif
    if(.not. (DECinfo%OnlyOccPart.or.DECinfo%OnlyVirtPart)) then
       write(DECinfo%output,'(1X,a,f16.10)') 'FOP Loop: Lagrangian Fragment energy       :', &
            & Fragment%LagFOP
    end if

    PrintDiff: if(iter/=0) then ! no point in printing energy differences for initial energy
       if(.not. DECinfo%OnlyVirtPart) then
          write(DECinfo%output,'(1X,a,f16.10)') 'FOP Loop: Occupied energy diff             :', &
               & OccEnergyDiff
       endif
       if(.not. DECinfo%OnlyOccPart) then
          write(DECinfo%output,'(1X,a,f16.10)') 'FOP Loop: Virtual energy diff              :', &
               & VirtEnergyDiff
       endif
       if(.not. (DECinfo%OnlyOccPart.or.DECinfo%OnlyVirtPart)) then
          write(DECinfo%output,'(1X,a,f16.10)') 'FOP Loop: Lagrangian energy diff           :', &
               & LagEnergyDiff
       end if
    end if PrintDiff

    write(DECinfo%output,'(1X,a)') 'FOP========================================================='
    write(DECinfo%output,*) 'FOP'


  end subroutine fragopt_print_info



  !> \brief Given a prioritized list of atoms obtained from an "old" fragment
  !> set logical atom vector defining fragment AOS for "new" fragment, which is smaller than
  !> (or the same size) as the old one.
  !> \author Kasper Kristensen
  !> \date April 2013
  subroutine set_reduced_fragmentAOS(MyAtom,natomsAOS_old,natomsfull,natomsAOS_new,atomlist,&
       & occ_atoms,unocc_atoms)

    implicit none
    !> Central atom in fragment (will of course always be included in AOS)
    integer,intent(in) :: MyAtom
    !> Number of AOS atoms in old fragment
    integer,intent(in) :: natomsAOS_old
    !> Number of atoms in full molecule
    integer,intent(in) :: natomsfull
    !> Number of AOS atoms in new fragment
    integer,intent(in) :: natomsAOS_new
    !> Prioritized list of atoms in old fragment
    integer,intent(in),dimension(natomsAOS_old) :: atomlist
    !> Which AOS atoms to include in new fragment (occ space)
    logical,intent(inout),dimension(natomsfull) :: occ_atoms
    !> Which AOS atoms to include in new fragment (unocc space)
    !> (occ_atoms = unocc_atoms but it is convenient to have both outputs in fragment optimization)
    logical,intent(inout),dimension(natomsfull) :: unocc_atoms
    integer :: counter,idx,i

    ! Always include central atom in fragment
    occ_atoms=.false.
    occ_atoms(MyAtom) = .true.

    ! Include atoms until there are natomsAOS_new atoms in AOS
    ! --------------------------------------------------------
    counter = 1  ! Central atom has already been include

    IncludeLoop: do i=1,natomsAOS_old
       if(counter==natomsAOS_new) exit IncludeLoop
       idx = atomlist(i)
       if(.not. occ_atoms(idx)) then
          occ_atoms(idx) = .true.
          counter = counter +1
       end if
    end do IncludeLoop

    ! Sanity check
    if(counter/=natomsAOS_new) then
       print *, 'counter / natomsAOS_new', counter,natomsAOS_new
       call lsquit('set_reduced_fragmentAOS: Counter is wrong',-1)
    end if
    unocc_atoms = occ_atoms

  end subroutine set_reduced_fragmentAOS


  !> Expand both occ and virt site fragment
  !> date: august-2011
  !> author: Ida-Marie Hoyvik
  subroutine Expandfragment(Occ,Virt,Track,natoms,&
       & nocc_per_atom,nunocc_per_atom,Converged,&
       & nAtomsWithOccOrb,nAtomsWithVirtOrb,Inputstep,&
       & ExpandOcc,ExpandVirt)
    implicit none
    integer,intent(in)         :: natoms
    logical,dimension(natoms)  :: Occ,Virt,Temp
    integer,intent(in)         :: nocc_per_atom(natoms)
    integer,intent(in)         :: nunocc_per_atom(natoms)
    integer,intent(in)         :: Track(natoms)
    logical,intent(inout)      :: Converged
    logical,intent(in)         :: ExpandOcc,ExpandVirt
    integer,intent(in)         :: nAtomsWithOccOrb,nAtomsWithVirtOrb
    integer,optional,intent(in):: InputStep
    !
    integer :: i,indx,counter,step
    logical :: VirtConverged,OccConverged
    OccConverged = .FALSE.
    VirtConverged = .FALSE.
    IF(PRESENT(InputStep))THEN
       step = InputStep
    ELSE
       step = DECinfo%Frag_Exp_Size
    ENDIF

    IF(ExpandOcc)THEN
       !Exand occ space first
       counter = 0
       Temp = Occ
       do i=1,natoms
          indx=Track(i)
          if ((.not. Occ(indx)) .and. &
               & (nocc_per_atom(indx) > 0)) then
             Temp(indx)=.true.
             counter = counter+1
             if (counter == step) exit
          end if
       end do
       Occ = Temp
       IF(COUNT(Occ).EQ.nAtomsWithOccOrb)THEN
          OccConverged = .TRUE.
       ENDIF
    ELSE
       OccConverged = .TRUE.       
    ENDIF

    !Expand virt space
    IF(ExpandVirt)THEN
       counter =0
       Temp = Virt
       do i=1,natoms
          indx= Track(i)
          if ((.not. Virt(indx)) .and. &
               & (nunocc_per_atom(indx)> 0)) then
             Temp(indx)=.true.
             counter=counter+1
             if (counter == step) exit
          end if
       end do
       Virt = Temp
       IF(COUNT(Virt).EQ.nAtomsWithVirtOrb)THEN
          VirtConverged = .TRUE.
       ENDIF
    ELSE
       VirtConverged = .TRUE.       
    ENDIF

    IF(VirtConverged.AND.OccConverged)THEN
       Converged = .TRUE.
    ELSE
       Converged = .FALSE.       
    ENDIF

  end subroutine Expandfragment

  !> Expand both occ and virt site fragment
  !> date: august-2011
  !> author: Ida-Marie Hoyvik
  subroutine GetnAtomsWithOrb(natoms,nocc_per_atom,&
       & nunocc_per_atom,nAtomsWithOccOrb,nAtomsWithVirtOrb)
    implicit none
    integer,intent(in)        :: natoms
    integer,intent(in)        :: nocc_per_atom(natoms)
    integer,intent(in)        :: nunocc_per_atom(natoms)
    integer,intent(inout)     :: nAtomsWithOccOrb,nAtomsWithVirtOrb
    !
    integer                   :: indx
    nAtomsWithOccOrb = 0
    nAtomsWithVirtOrb = 0
    do indx=1,natoms
       if (nocc_per_atom(indx) > 0) then
          nAtomsWithOccOrb = nAtomsWithOccOrb + 1
       end if
       if (nunocc_per_atom(indx) > 0) then
          nAtomsWithVirtOrb = nAtomsWithVirtOrb + 1
       end if
    end do
  end subroutine GetnAtomsWithOrb

  !> \brief Set logical vectors defining fragment to contain a fixed number of AOS atoms
  !> based on a prioritized track list.
  !> \date April 2013
  !> \author Kasper Kristensen
  subroutine Set_fragment_fixed_AOSatoms(natomsAOS,natoms,&
       & nocc_per_atom,nunocc_per_atom,Track,Occ,Unocc)
    implicit none
    !> Number of requested AOS atoms
    integer,intent(in)        :: natomsAOS
    !> Number of atoms in full molecule
    integer,intent(in)        :: natoms
    !> Number of occ orbitals assigned to each atom
    integer,intent(in)        :: nocc_per_atom(natoms)
    !> Number of unocc orbitals assigned to each atom
    integer,intent(in)        :: nunocc_per_atom(natoms)
    !> Prioritized track list of atoms (the first natomsAOS from this list are included)
    integer,intent(in)        :: Track(natoms)
    !> Which AOS atoms to include for occ space
    logical,dimension(natoms),intent(inout) :: Occ
    !> Which AOS atoms to include for unocc space 
    !> (will in general be identical to occ but keep it general)
    logical,dimension(natoms),intent(inout) :: Unocc
    integer                   :: i,indx,counter

    ! Sanity check
    if(natomsAOS>natoms) then
       print *, 'natoms,natomsAOS',natoms,natomsAOS
       call lsquit('Set_fragment_fixed_AOSatoms: #AOS atoms larger than # atoms in molecule!',-1)
    end if

    ! Include natomsAOS atoms for occ space
    counter = 0
    Occ=.false.
    IncludeOcc: do i=1,natoms
       indx=Track(i)
       if (nocc_per_atom(indx) > 0) then
          Occ(indx)=.true.
          counter = counter+1
          if (counter == natomsAOS) exit IncludeOcc
       end if
    end do IncludeOcc

    ! Include natomsAOS atoms for unocc space
    counter = 0
    Unocc=.false.
    IncludeUnocc: do i=1,natoms
       indx=Track(i)
       if (nunocc_per_atom(indx) > 0) then
          Unocc(indx)=.true.
          counter = counter+1
          if (counter == natomsAOS) exit IncludeUnocc
       end if
    end do IncludeUnocc


  end subroutine Set_fragment_fixed_AOSatoms

  subroutine ReduceSpace_binary(MyFragment,nocc_full,TrackSortedOccContribs,&
       & NonOccContribsRemoval,SecondaryTrackList,OccOrVirt,OccOrbAOS,nocc)
    implicit none
    !> Fragment info
    type(decfrag),intent(inout) :: MyFragment
    !> Number of orbitals (occ/virt) in full molecule
    integer,intent(in)        :: nocc_full
    !> Tracking info of Sorted Contributions to the Virtual fragment 
    !> energy from each individual Occ/virt orbital
    integer,intent(in),dimension(nocc_full) :: TrackSortedOcccontribs
    !> Tracking info of Sorted Distances from Orbital to MyAtom   - OrbOccDistTrackMyAtom
    integer,intent(in),dimension(nocc_full) :: SecondaryTrackList 
    !> Occupied ('O') or virtual orbitals ('V') under consideration
    character(len=1),intent(in) :: OccOrVirt
    !> Logical vector telling which Occ orbitals are included in AOS 
    !> (true) and not included (false)
    logical,intent(inout),dimension(nocc_full) :: occOrbAOS
    !> Number of desired occ orbitals AFTER the reduction has been carried out
    integer,intent(in) :: nOcc
    !> Use Distance compared to Energy contribution to reduce space
    logical,intent(in) :: NonOccContribsRemoval
    !local orbitals
    integer :: i,ii    
    !loop over the first nocc (with highest OccContribs(i)) 
    OccOrbAOS = .false.
    call SanityCheckOrbAOS(MyFragment,nocc_full,OccOrVirt,OccOrbAOS)
    IF(NonOccContribsRemoval.AND.(DECinfo%onlyOccPart.AND.(OccOrVirt.EQ.'O')))THEN
       !add/remove Occupied Orbitals based on Distance
       do i=1,nocc
          ii = SecondaryTrackList(i)
          OccOrbAOS(ii) = .true.
       enddo
    ELSEIF(NonOccContribsRemoval.AND.(DECinfo%onlyVirtPart.AND.(OccOrVirt.EQ.'V')))THEN
       !add/remove Virtual Orbitals based on Distance
       do i=1,nocc
          ii = SecondaryTrackList(i)
          OccOrbAOS(ii) = .true.
       enddo
    ELSE
       !add/remove Occupied/Virtual Orbitals based on energy contribution
       do i=1,nocc
          ii = TrackSortedOccContribs(i)
          OccOrbAOS(ii) = .true.
       enddo
    ENDIF
  end subroutine ReduceSpace_binary

  !> \brief Reduce occupied or virtual AOS for fragment using the individual orbital
  !> contributions according to the input threshold.
  !> \author Kasper Kristensen
  !> \date December 2011
  subroutine ReduceSpace_orbitalspecific(MyFragment,norb_full,contributions,OccOrVirt,&
       &Thresh,OrbAOS,Nafter)

    implicit none
    !> Fragment info
    type(decfrag),intent(inout) :: MyFragment
    !> Number of orbitals (occ OR virt) in full molecule
    integer,intent(in)        :: norb_full
    !> Contributions to the fragment energy from each individual orbital
    real(realk),intent(in),dimension(norb_full) :: contributions
    !> Threshold - all orbitals with contributions smaller than this are excluded
    real(realk),intent(in) :: Thresh
    !> Occupied ('O') or virtual orbitals ('V') under consideration
    character(len=1),intent(in) :: OccOrVirt
    !> Logical vector telling which orbitals are included in AOS (true) and not included (false)
    logical,intent(inout),dimension(norb_full) :: OrbAOS
    !> Number of orbitals AFTER the reduction has been carried out
    integer,intent(inout) :: Nafter
    integer :: i,Nbefore,Nexcl
    logical :: doOccPart,doVirtPart

    ! Number of orbitals at input
    Nbefore = count(OrbAOS)

    !For MP2 we always do both occ and virt
    doOccPart = (.NOT.DECinfo%OnlyVirtPart) .OR.(MyFragment%ccmodel==MODEL_MP2)
    doVirtPart = (.NOT.DECinfo%OnlyOccPart) .OR.(MyFragment%ccmodel==MODEL_MP2)


    ! Exclude orbitals with contributions smaller than threshold
    if(OccOrVirt=='O') then ! checking occupied orbitals
       !contributions = OccContribs 
       IF(doVirtPart)THEN
          do i=1,norb_full
             if (abs(contributions(i)) < Thresh) OrbAOS(i)=.false.
          end do
       ELSE
          !Occupied partitioning only:
          !OccContribs have not been calculated and we cannot remove any
          !Occupied Orbitals based on this 
       ENDIF
    elseif(OccOrVirt=='V') then
       !contributions = VirtContribs
       IF(doOccPart)THEN
          do i=1,norb_full
             if (abs(contributions(i)) < Thresh) OrbAOS(i)=.false.
          end do
       ELSE
          !Virtual partitioning only:
          !VirtContribs have not been calculated and we cannot remove any
          !Virtual Orbitals based on this 
       ENDIF
    endif

    ! Sanity check: The orbitals assigned to the central atom in the fragment should ALWAYS be included
    call SanityCheckOrbAOS(MyFragment,norb_full,OccOrVirt,OrbAOS)

    Nafter = count(OrbAOS)
    Nexcl = Nbefore-Nafter

    write(DECinfo%output,'(a,i4)') ' FOP Number of orbitals excluded: ', Nexcl

  end subroutine ReduceSpace_orbitalspecific


  !> \brief For a given model, get the occupied, virtual and Lagragian fragment energies
  !> to use for fragment optimization, i.e. simply copy the relevant energies from
  !> "fragment%energies" to fragment%EoccFOP, fragment%EvirtFOP, and fragment%LagFOP.
  !> (See description of EoccFOP, EvirtFOP, and LagFOP in decfrag type definition).
  !> \author Kasper Kristensen
  !> \date March 2013
  subroutine get_occ_virt_lag_energies_fragopt(Fragment)
    implicit none
    type(decfrag),intent(inout) :: fragment
    ! MODIFY FOR NEW MODEL
    ! If you implement a new model, please copy fragment%energies(?) for your model,
    ! see decfrag type def to determine the "?".

    fragment%EoccFOP = 0.0_realk
    fragment%EoccFOP_Corr = 0.0_realk
    fragment%EvirtFOP = 0.0_realk
    fragment%LagFOP = 0.0_realk

    select case(fragment%ccmodel)
    case(MODEL_MP2)
       ! MP2
       fragment%LagFOP = fragment%energies(FRAGMODEL_LAGMP2)
       fragment%EoccFOP = fragment%energies(FRAGMODEL_OCCMP2)
       fragment%EvirtFOP = fragment%energies(FRAGMODEL_VIRTMP2)

#ifdef MOD_UNRELEASED         
       if(Decinfo%F12) then
          ! MP2-F12: F12-correction
          ! fragment%EoccFOP_Corr = fragment%energies(FRAGMODEL_MP2f12)      
       endif
#endif

    case(MODEL_CC2)
       ! CC2
       fragment%EoccFOP = fragment%energies(FRAGMODEL_OCCCC2)
       fragment%EvirtFOP = fragment%energies(FRAGMODEL_VIRTCC2)
       ! simply use average of occ and virt energies since Lagrangian is not yet implemented
       fragment%LagFOP =  0.5_realk*(fragment%EoccFOP+fragment%EvirtFOP)   
    case(MODEL_RPA)
       ! RPA
       fragment%EoccFOP = fragment%energies(FRAGMODEL_OCCRPA)
       fragment%EvirtFOP = fragment%energies(FRAGMODEL_VIRTRPA)
       !print *,"JOHANNES: CURRENTLY LAGRANGIAN ENERGY IS NOT CONSIDERED; PLEASE IMPLEMENT"
       ! simply use average of occ and virt energies since Lagrangian is not yet implemented
       fragment%LagFOP =  0.5_realk*(fragment%EoccFOP+fragment%EvirtFOP)   
       !fragment%LagFOP = fragment%energies(FRAGMODEL_LAGRPA)
    case(MODEL_CCSD)
       ! CCSD
       fragment%EoccFOP = fragment%energies(FRAGMODEL_OCCCCSD)
       fragment%EvirtFOP = fragment%energies(FRAGMODEL_VIRTCCSD)
       ! simply use average of occ and virt energies since Lagrangian is not yet implemented
       fragment%LagFOP =  0.5_realk*(fragment%EoccFOP+fragment%EvirtFOP)

#ifdef MOD_UNRELEASED         
       if(Decinfo%F12) then
          ! CCSD-F12: F12-correction
          fragment%EoccFOP_Corr = fragment%energies(FRAGMODEL_CCSDf12)      
       endif
#endif

#ifdef MOD_UNRELEASED
    case(MODEL_CCSDpT)
       ! CCSD(T): CCSD contribution + (T) contribution
       fragment%EoccFOP = fragment%energies(FRAGMODEL_OCCCCSD) + fragment%energies(FRAGMODEL_OCCpT)
       fragment%EvirtFOP = fragment%energies(FRAGMODEL_VIRTCCSD) + fragment%energies(FRAGMODEL_VIRTpT)
       ! simply use average of occ and virt energies since Lagrangian is not yet implemented
       fragment%LagFOP =  0.5_realk*(fragment%EoccFOP+fragment%EvirtFOP)
       !endif mod_unreleased
#endif
    case(MODEL_RIMP2)
       ! RI-MP2
       fragment%EoccFOP = fragment%energies(FRAGMODEL_OCCRIMP2)
       fragment%EvirtFOP = fragment%energies(FRAGMODEL_VIRTRIMP2)

    case default
       write(DECinfo%output,*) 'WARNING: get_occ_virt_lag_energies_fragopt needs implementation &
            & for model:', fragment%ccmodel
    end select

  end subroutine get_occ_virt_lag_energies_fragopt


  !> \brief After the fragment optimization, make sure that the energies stored
  !> in fragment%LagFOP, fragment%EoccFOP, and fragment%EvirtFOP are copied correctly 
  !> to the general fragment%energies arrays.
  !> This is effectively the inverse routine of get_occ_virt_lag_energies_fragopt.
  !> (See description of energies, EoccFOP, EvirtFOP, and LagFOP in decfrag type definition).
  !> \author Kasper Kristensen
  !> \date March 2013
  subroutine set_energies_decfrag_structure_fragopt(Fragment)
    implicit none
    type(decfrag),intent(inout) :: fragment
    ! MODIFY FOR NEW MODEL 
    ! If you implement a new model, please set fragment%energies(?) for your model,
    ! see FRAGMODEL_* definitions in dec_typedef.F90 to determine the "?".
    
    select case(fragment%ccmodel)
    case(MODEL_MP2)
       ! MP2
       fragment%energies(FRAGMODEL_LAGMP2) = fragment%LagFOP 
       fragment%energies(FRAGMODEL_OCCMP2) = fragment%EoccFOP
       fragment%energies(FRAGMODEL_VIRTMP2) = fragment%EvirtFOP 

#ifdef MOD_UNRELEASED 
       if(DECinfo%F12) then
          ! F12: F12-correction or the Single Fragments
          fragment%energies(FRAGMODEL_MP2f12) = fragment%EoccFOP_Corr
       endif
#endif

    case(MODEL_CC2)
       ! CC2
       fragment%energies(FRAGMODEL_OCCCC2) = fragment%EoccFOP
       fragment%energies(FRAGMODEL_VIRTCC2) = fragment%EvirtFOP
    case(MODEL_RPA)
       ! RPA
       fragment%energies(FRAGMODEL_OCCRPA) = fragment%EoccFOP
       fragment%energies(FRAGMODEL_VIRTRPA) = fragment%EvirtFOP
    case(MODEL_CCSD)
       ! CCSD
       fragment%energies(FRAGMODEL_OCCCCSD) = fragment%EoccFOP 
       fragment%energies(FRAGMODEL_VIRTCCSD) = fragment%EvirtFOP

#ifdef MOD_UNRELEASED 
       if(DECinfo%F12) then
          ! F12: F12-correction or the Single Fragments
          fragment%energies(FRAGMODEL_CCSDf12) = fragment%EoccFOP_Corr
       endif
#endif

#ifdef MOD_UNRELEASED 
    case(MODEL_CCSDpT)
       ! (T) contribution =  CCSD(T) contribution minus CCSD contribution
       fragment%energies(FRAGMODEL_OCCpT) = fragment%EoccFOP - fragment%energies(FRAGMODEL_OCCCCSD)
       fragment%energies(FRAGMODEL_VIRTpT) = fragment%EvirtFOP - fragment%energies(FRAGMODEL_VIRTCCSD)
!endif mod_unreleased
#endif
    case(MODEL_RIMP2)
       ! RI-MP2
       fragment%energies(FRAGMODEL_OCCRIMP2) = fragment%EoccFOP
       fragment%energies(FRAGMODEL_VIRTRIMP2) = fragment%EvirtFOP 

    case default
       write(DECinfo%output,*) 'WARNING: get_occ_virt_lag_energies_fragopt needs implementation &
            & for model:', fragment%ccmodel
    end select

  end subroutine set_energies_decfrag_structure_fragopt


  !> \brief When fragment energies have been calculated, put them
  !> into the energies array in fragment structure according to the given model
  !> (see FRAGMODEL_* in dec_typedef.F90).
  !> \author Kasper Kristensen
  !> \date October 2013
  subroutine put_fragment_energy_contribs_main(Eocc,Evirt,MyFragment)
    implicit none
    !> Occupied and virtual partitioning scheme energies
    real(realk),intent(in) :: Eocc, Evirt
    !> Atomic or pair fragment
    type(decfrag),intent(inout) :: MyFragment


    ! Put energies into their proper place in the MyFragment%energies array
    ! according to the CC model used for the fragment
    call put_fragment_energy_contribs(MyFragment%ccmodel,Eocc,Evirt,MyFragment%energies)


    ! Special case: When some pairs are treated at the MP2 level, while the
    ! actual CC model is higher in the hierarchy (e.g. CCSD), we need to 
    ! put the MP2 fragment energies into the energy array for the more accurate model.
    if(MyFragment%ccmodel==MODEL_MP2 .and. DECinfo%ccmodel/=MODEL_MP2) then
       call put_fragment_energy_contribs(DECinfo%ccmodel,Eocc,Evirt,MyFragment%energies)
    end if

  end subroutine put_fragment_energy_contribs_main



  ! MODIFY FOR NEW MODEL
  !> \brief When fragment energies have been calculated, put them
  !> into the energies array according to the given model
  !> (see FRAGMODEL_* in dec_typedef.F90).
  !> \author Kasper Kristensen
  !> \date October 2013
  subroutine put_fragment_energy_contribs(ccmodel,Eocc,Evirt,energies)
    implicit none
    !> Which CC model
    integer,intent(in) :: ccmodel
    !> Occupied and virtual partitioning scheme energies
    real(realk),intent(in) :: Eocc, Evirt
    !> Energies array
    real(realk),intent(inout) :: energies(ndecenergies)


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
       energies(FRAGMODEL_OCCRPA) = Eocc   ! occupied
       energies(FRAGMODEL_VIRTRPA) = Evirt   ! virtual
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
    case default
       call lsquit('case unknown in put_fragment_energy_contribs',-1)
    end select

  end subroutine put_fragment_energy_contribs

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
   subroutine optimize_atomic_fragment_CLEAN(MyAtom,AtomicFragment,natoms, &
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
      logical :: full_mol
      !> Timings for frag opt
      real(realk), pointer :: times_fragopt(:)  


      !==================================================================================!
      !                       Initialization of various things...                        !
      !==================================================================================!
      !
      ! Start timings:
      times_fragopt => null()
      call dec_fragment_time_init(times_fragopt)

      write(DECinfo%output,'(a)')    ' FOP'
      write(DECinfo%output,'(a)')    ' FOP ==============================================='
      write(DECinfo%output,'(a,i4)') ' FOP  Site fragment generator for fragment,',MyAtom
      write(DECinfo%output,'(a)')    ' FOP ==============================================='
      write(DECinfo%output,'(a)')    ' FOP'


      ! Get number of orbital in EOS spaces:
      nocc_per_atom=get_number_of_orbitals_per_atom(OccOrbitals,no,natoms,.true.)
      nvir_per_atom=get_number_of_orbitals_per_atom(VirOrbitals,nv,natoms,.true.)

      ! Only do fragment optimization if there are orbitals assigned to central atom.
      if( (nocc_per_atom(MyAtom) == 0) .and. (nvir_per_atom(MyAtom) == 0) ) then
         write(DECinfo%output,*) 'FOP Skipping optimization of fragment ', MyAtom
         AtomicFragment%LagFOP   = 0E0_realk
         AtomicFragment%EoccFOP  = 0E0_realk
         AtomicFragment%EvirtFOP = 0E0_realk
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

      ! Get information on how the expansion should be performed:
      call mem_alloc(exp_list_occ,no)
      call mem_alloc(exp_list_vir,nv)

      call define_frag_expansion(no,nv,natoms,MyAtom,MyMolecule,AtomicFragment, &
         & exp_list_occ,exp_list_vir,nexp_occ,nexp_vir,ninit_occ,ninit_vir)



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
      call fragment_energy_and_prop(AtomicFragment) 
      ! Print initial fragment information
      if (full_mol) write(DECinfo%output,*) 'FOP Expansion Include Full Molecule !!!'
      call fragopt_print_info(AtomicFragment,0.0E0_realk,0.0E0_realk,0.0E0_realk,0)



      !==================================================================================!
      !                             Enter Fragment Expansion                             !
      !==================================================================================!
      !
      ! Do expansion only if the initial fragment does not include the full molecule:
      if (.not.full_mol) then
         call fragment_expansion_procedure(Occ_AOS,Vir_AOS,AtomicFragment,no,nv, &
            & exp_list_occ,exp_list_vir,nexp_occ,nexp_vir,MyAtom,MyMolecule, &
            & OccOrbitals,VirOrbitals,mylsitem)
      end if



      !==================================================================================!
      !                  Transition from expansion to reduction loop                     !
      !==================================================================================!

      write(DECinfo%output,*) 'FOP'
      write(DECinfo%output,*) 'FOP ============================================='
      write(DECinfo%output,*) 'FOP  Expansion has converged. We start reduction '
      write(DECinfo%output,*) 'FOP ============================================='
      write(DECinfo%output,*) 'FOP'

      ! Deallocate exp list and allocate red ones:
      call mem_dealloc(exp_list_occ)
      call mem_dealloc(exp_list_vir)
      call mem_alloc(red_list_occ,no)
      call mem_alloc(red_list_vir,nv)

      ! Overwrite ccmodel for myatom with the reduction required ccmodel
      MyMolecule%ccmodel(MyAtom,Myatom) = DECinfo%fragopt_red_model
      AtomicFragment%ccmodel = DECinfo%fragopt_red_model

      ! if expansion and reduction ccmodels are different, we need to calculate
      ! the energy of the expanded fragment using the new model:
      if(DECinfo%fragopt_exp_model /= DECinfo%fragopt_red_model) then
         call fragment_energy_and_prop(AtomicFragment)
         write(DECinfo%output,'(2a)') ' FOP Calculated ref atomic fragment energy for relevant CC model: ', &
           & DECinfo%cc_models(MyMolecule%ccmodel(MyAtom,Myatom))
         call fragopt_print_info(AtomicFragment,0.0E0_realk,0.0E0_realk,0.0E0_realk,0)
      end if

      ! Get information on how the reduction should be performed:
      call define_frag_reduction(no,nv,natoms,MyAtom,MyMolecule,AtomicFragment, &
         & red_list_occ,red_list_vir,nred_occ,nred_vir)

      ! Perform reduction:
      call fragment_reduction_procedure_wrapper(AtomicFragment,no,nv,red_list_occ, &
            & red_list_vir,Occ_AOS,Vir_AOS,MyAtom,MyMolecule,OccOrbitals, &
            & VirOrbitals,mylsitem,nred_occ,nred_vir,DECinfo%FOT)

      !==================================================================================!
      !                              Finalize subroutine                                 !
      !==================================================================================!


      ! Deallocation:
      call mem_dealloc(red_list_occ)
      call mem_dealloc(red_list_vir)
      call mem_dealloc(Occ_AOS)
      call mem_dealloc(Vir_AOS)

      if(freebasisinfo) then
         call atomic_fragment_free_basis_info(AtomicFragment)
      end if

      ! Ensure that energies in fragment are set consistently
      call set_energies_decfrag_structure_fragopt(AtomicFragment)

      ! Get final Timings:
      call dec_fragment_time_get(times_fragopt)
      call dec_time_evaluate_efficiency_frag(AtomicFragment,times_fragopt, &
         & AtomicFragment%ccmodel,'Fragment optmization')
 
      ! Restore the original CC model 
      ! (only relevant if expansion and/or reduction was done using the MP2 model, but it doesn't hurt)
      MyMolecule%ccmodel(MyAtom,Myatom) = DECinfo%ccmodel

   end subroutine optimize_atomic_fragment_CLEAN


   ! Purpose: Perform expansion procedure based on occ/vir priority list. In each step
   !          the fragment AOS is increased using the nexp most important orbitals from
   !          occupied and virtual priority lists.
   !          The expansion is converged when the difference in energie(s) between two
   !          is smaller than the FOT.
   !
   ! Author:  Pablo Baudin (based on previous work by Kasper Kristensen & Thomas Kjaergaard)
   ! Date:    July 2014
   subroutine fragment_expansion_procedure(Occ_AOS,Vir_AOS,AtomicFragment,no,nv, &
            & occ_priority_list,vir_priority_list,nexp_occ,nexp_vir,MyAtom,MyMolecule, &
            & OccOrbitals,VirOrbitals,mylsitem)

      implicit none

      !> Number of occupied orbitals in molecule
      integer, intent(in) :: no
      !> Number of virtual orbitals in molecule
      integer, intent(in) :: nv
      !> Atomic fragment to be optimized
      type(decfrag),intent(inout)        :: AtomicFragment
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
      !> All unoccupied orbitals
      type(decorbital), dimension(nv), intent(in) :: VirOrbitals
      !> Logical vector telling which orbital is include in the fragment
      logical, intent(inout) :: Occ_AOS(no), Vir_AOS(nv)
      !> Integral information
      type(lsitem), intent(inout)       :: mylsitem
      
      integer :: iter
      logical :: full_mol, expansion_converged
      real(realk) :: LagEnergy_old, OccEnergy_old, VirEnergy_old
      real(realk) :: LagEnergy_dif, OccEnergy_dif, VirEnergy_dif


      EXPANSION_LOOP: do iter=1,DECinfo%maxiter

         ! Save current fragment energies
         LagEnergy_old = AtomicFragment%LagFOP
         OccEnergy_old = AtomicFragment%EoccFOP
         VirEnergy_old = AtomicFragment%EvirtFOP

         ! Expand fragment:
         call expand_fragment(no,nv,occ_priority_list,vir_priority_list,nexp_occ,nexp_vir, &
            & MyAtom,MyMolecule,OccOrbitals,VirOrbitals,Occ_AOS,Vir_AOS,full_mol)

         ! Free old fragment and initialize the expanded fragment:
         call atomic_fragment_free(AtomicFragment)
         call atomic_fragment_init_orbital_specific(MyAtom,nv,no,Vir_AOS,Occ_AOS,OccOrbitals, &
            & VirOrbitals,MyMolecule,mylsitem,AtomicFragment,.true.,.false.)

         ! Get new fragment energy:
         call fragment_energy_and_prop(AtomicFragment)

         ! Energy differences
         LagEnergy_dif = abs(LagEnergy_old - AtomicFragment%LagFOP)
         OccEnergy_dif = abs(OccEnergy_old - AtomicFragment%EoccFOP)
         VirEnergy_dif = abs(VirEnergy_old - AtomicFragment%EvirtFOP)
                                                                                                                 
         ! print expanded fragment information:
         if (full_mol) write(DECinfo%output,*) 'FOP Expansion Include Full Molecule !!!'
         call fragopt_print_info(AtomicFragment,LagEnergy_dif,OccEnergy_dif,VirEnergy_dif,iter)    

         ! Test if fragment energy (or energies) are converged to FOT precision
         call fragopt_check_convergence(LagEnergy_dif,OccEnergy_dif,VirEnergy_dif, &
            & DECinfo%FOT,expansion_converged)

         ! Set the expansion to be converged if the current fragment include the full molecule:
         if (full_mol) expansion_converged = .true.


         ! Exit loop if we are converged
         ExpansionConvergence: if(expansion_converged) then
            write(DECinfo%output,*) 'FOP Fragment expansion converged in iteration ', iter
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
        & VirOrbitals,mylsitem,no_gap,nv_gap,FOT)

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
     !> All unoccupied orbitals
     type(decorbital), dimension(nv), intent(in) :: VirOrbitals
     !> Integral information
     type(lsitem), intent(inout)       :: mylsitem
     !> minimum gap in number of orbital allowed between the 
     !  last two steps of the binary search.
     integer,intent(in) :: no_gap, nv_gap
     !> Fragment optimization threshold to use in reduction
     real(realk),intent(in) :: FOT
     real(realk) :: LagEnergy_exp,OccEnergy_exp,VirEnergy_exp,FOTincreased
     type(decfrag) :: ReducedFragment
     integer :: i

     ! At this point AtomicFragment correspondings to the expanded fragment
     ! We store the atomic fragment energies of the expanded fragment
     LagEnergy_exp = AtomicFragment%LagFOP
     OccEnergy_exp = AtomicFragment%EoccFOP
     VirEnergy_exp = AtomicFragment%EvirtFOP

     ! Determine AtomicFragment according to main FOT
     call fragment_reduction_procedure(AtomicFragment,no,nv,occ_priority_list, &
          & vir_priority_list,Occ_AOS,Vir_AOS,MyAtom,MyMolecule,OccOrbitals, &
          & VirOrbitals,mylsitem,no_gap,nv_gap,FOT)
     write(DECinfo%output,'(1X,a,g14.3,4i8)') 'FOP reduction: Atom,FOT,O,V,B',MyAtom,FOT,&
          & AtomicFragment%noccAOS,AtomicFragment%nunoccAOS,AtomicFragment%nbasis
     
     ! Loop over different increased FOTs to determine reduced spaces
     ! *************************************************************
     FOTincreased=FOT
     do i=1,DECinfo%nFRAGSred

        ! Initialize ReducedFragment identical to AtomicFragment
        call atomic_fragment_init_integer_list(MyAtom,nv, no, AtomicFragment%nunoccAOS,&
             & AtomicFragment%noccAOS,AtomicFragment%unoccAOSidx,AtomicFragment%occAOSidx,&
             & OccOrbitals,VirOrbitals,MyMolecule,mylsitem,ReducedFragment,.true.,.false.)

        ! Set initial energies for ReducedFragment equal to the ones for the expanded fragment
        ReducedFragment%LagFOP   = LagEnergy_exp
        ReducedFragment%EoccFOP  = OccEnergy_exp
        ReducedFragment%EvirtFOP = VirEnergy_exp

        ! Increase FOT by scaling factor
        FOTincreased = FOTincreased*DECinfo%FOTscaling

        ! Determine ReducedFragment according to the scaled FOT
        call fragment_reduction_procedure(ReducedFragment,no,nv,occ_priority_list, &
             & vir_priority_list,Occ_AOS,Vir_AOS,MyAtom,MyMolecule,OccOrbitals, &
             & VirOrbitals,mylsitem,no_gap,nv_gap,FOTincreased)
        
        ! Save information about ReducedFragment AOS in AtomicFragment structure
        AtomicFragment%REDfrags(i)%noccAOS   = ReducedFragment%noccAOS
        AtomicFragment%REDfrags(i)%nunoccAOS = ReducedFragment%nunoccAOS
        call mem_alloc(AtomicFragment%REDfrags(i)%occAOSidx,AtomicFragment%REDfrags(i)%noccAOS)
        AtomicFragment%REDfrags(i)%occAOSidx = ReducedFragment%occAOSidx
        call mem_alloc(AtomicFragment%REDfrags(i)%unoccAOSidx,AtomicFragment%REDfrags(i)%nunoccAOS)
        AtomicFragment%REDfrags(i)%unoccAOSidx = ReducedFragment%unoccAOSidx
        AtomicFragment%REDfrags(i)%FOT = FOTincreased

        ! Print summary (delete this at some point but nice to have for analysis now)
        write(DECinfo%output,'(1X,a,g14.3,4i8)') 'FOP reduction: Atom,FOT,O,V,B',MyAtom,&
             & FOTincreased, ReducedFragment%noccAOS,ReducedFragment%nunoccAOS,&
             & ReducedFragment%nbasis

        ! Done with reduced fragment
         call atomic_fragment_free(ReducedFragment)

     end do


   end subroutine fragment_reduction_procedure_wrapper




   ! Purpose: Perform reduction procedure based on occ/vir priority list. We perform
   !          a binary search on the priority list and accept a step based on energy
   !          criterions (dE_occ, dE_vir). The binary search stops for one space when
   !          the difference in number of orbitals between two steps is lower than
   !          no_gap/nv_gap. By default the occupied spaced is reduced 1st.
   !
   ! Author:  Pablo Baudin
   ! Date:    July 2014
   subroutine fragment_reduction_procedure(AtomicFragment,no,nv,occ_priority_list, &
            & vir_priority_list,Occ_AOS,Vir_AOS,MyAtom,MyMolecule,OccOrbitals, &
            & VirOrbitals,mylsitem,no_gap,nv_gap,FOT)

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
      !> All unoccupied orbitals
      type(decorbital), dimension(nv), intent(in) :: VirOrbitals
      !> Integral information
      type(lsitem), intent(inout)       :: mylsitem
      !> minimum gap in number of orbital allowed between the 
      !  last two steps of the binary search.
      integer,intent(in) :: no_gap, nv_gap
      !> Fragment optimization threshold to use in reduction
      real(realk),intent(in) :: FOT

      !> energy error acceptance in reduction steps:
      real(realk) :: dE_occ, dE_vir
      integer :: no_exp, nv_exp, no_min, nv_min, no_max, nv_max
      integer :: no_old, nv_old, no_new, nv_new, iter
      logical :: redocc, redvir, step_accepted, reduction_converged, occ_red_conv, vir_red_conv
      logical, pointer :: OccAOS_old(:), VirAOS_old(:)
      real(realk) :: LagEnergy_exp, OccEnergy_exp, VirEnergy_exp
      real(realk) :: LagEnergy_old, OccEnergy_old, VirEnergy_old
      real(realk) :: LagEnergy_dif, OccEnergy_dif, VirEnergy_dif


      ! INITIALIZATION:
      ! ***************
      call mem_alloc(OccAOS_old,no)
      call mem_alloc(VirAOS_old,nv)

      no_exp = AtomicFragment%noccAOS
      nv_exp = AtomicFragment%nunoccAOS
      LagEnergy_exp = AtomicFragment%LagFOP
      OccEnergy_exp = AtomicFragment%EoccFOP
      VirEnergy_exp = AtomicFragment%EvirtFOP

      ! Set boundaries of the list (min = EOS, max = Expanded fragment AOS)
      ! note: if frozen core then core orbitals are already excluded from those numbers
      no_min = AtomicFragment%noccEOS
      nv_min = AtomicFragment%nunoccEOS
      no_max = AtomicFragment%noccAOS
      nv_max = AtomicFragment%nunoccAOS

      ! The initial (or old) condition correspond to the expanded fragment
      no_old = no_exp
      nv_old = nv_exp
      OccAOS_old = Occ_AOS
      VirAOS_old = Vir_AOS
      LagEnergy_old = LagEnergy_exp
      OccEnergy_old = OccEnergy_exp
      VirEnergy_old = VirEnergy_exp


      ! Define specific reduction parameters:
      ! -------------------------------------
      if (DECinfo%Frag_red_occ.and.(.not.DECinfo%Frag_red_virt)) then
         ! Start reducing occupied space:
         write(DECinfo%output,'(1X,a,/)') 'FOP: User chose to reduce occupied space first'
         redocc = .true.
         redvir = .false.
         dE_occ = DECinfo%frag_red1_thr*FOT
         dE_vir = DECinfo%frag_red2_thr*FOT
      else if (DECinfo%Frag_red_virt.and.(.not.DECinfo%Frag_red_occ)) then
         ! Start reducing virtual space:
         write(DECinfo%output,'(1X,a,/)') 'FOP: User chose to reduce virtual space first'
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
         call fragment_energy_and_prop(AtomicFragment)

         ! Energy differences
         LagEnergy_dif = abs(LagEnergy_exp - AtomicFragment%LagFOP)
         OccEnergy_dif = abs(OccEnergy_exp - AtomicFragment%EoccFOP)
         VirEnergy_dif = abs(VirEnergy_exp - AtomicFragment%EvirtFOP)

         ! print reduced fragment information:
         call fragopt_print_info(AtomicFragment,LagEnergy_dif,OccEnergy_dif,VirEnergy_dif,iter)


         ! CHECK CONVERGENCES:
         ! *******************
         ! Check if reduction step is accepted (Energy criterion):
         if (redocc.and.(.not.redvir)) then
            call fragopt_check_convergence(LagEnergy_dif,OccEnergy_dif,VirEnergy_dif, &
               & dE_occ,step_accepted)
         else if (redvir.and.(.not.redocc)) then
            call fragopt_check_convergence(LagEnergy_dif,OccEnergy_dif,VirEnergy_dif, &
               & dE_vir,step_accepted)
         else if (redocc.and.redvir) then
            call fragopt_check_convergence(LagEnergy_dif,OccEnergy_dif,VirEnergy_dif, &
               & FOT,step_accepted)
         end if

         ! If the step is accepted we need to save the frag info:
         if (step_accepted) then
            LagEnergy_old = AtomicFragment%LagFOP
            OccEnergy_old = AtomicFragment%EoccFOP
            VirEnergy_old = AtomicFragment%EvirtFOP
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
            if (DECinfo%PL > 1) write(DECinfo%output,*) 'BIN SEARCH: REDUCTION CONVERGED'

            if (.not.step_accepted) then
               ! The last binary search step did not succeed so  
               ! we init everything to the last successful step
               call atomic_fragment_free(AtomicFragment)
               call atomic_fragment_init_orbital_specific(MyAtom,nv,no,VirAOS_old,OccAOS_old, &
                    & OccOrbitals,VirOrbitals,MyMolecule,mylsitem,AtomicFragment,.true.,.false.)
               AtomicFragment%LagFOP   = LagEnergy_old
               AtomicFragment%EoccFOP  = OccEnergy_old
               AtomicFragment%EvirtFOP = VirEnergy_old
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
            if (DECinfo%PL > 1) write(DECinfo%output,*) 'BIN SEARCH: Step accepted'
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
            if (DECinfo%PL > 1) write(DECinfo%output,*) 'BIN SEARCH: Step NOT accepted'
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
      write(DECinfo%output,'(1X,a,/)') 'FOP'
      write(DECinfo%output,'(1X,a)') 'FOP========================================================='
      write(DECinfo%output,'(1X,a,i4)') 'FOP    LOCAL REDUCTION HAS CONVERGED FOR SITE',MyAtom
      write(DECinfo%output,'(1X,a)') 'FOP---------------------------------------------------------'
      write(DECinfo%output,'(1X,a,i4)')    'FOP Done: Fragment number                  :', MyAtom
      write(DECinfo%output,'(1X,a,i4)')    'FOP Done: Number of orbitals in virt total :', &
           & AtomicFragment%nunoccAOS
      write(DECinfo%output,'(1X,a,i4)')    'FOP Done: Number of orbitals in occ total  :', &
           & AtomicFragment%noccAOS
      write(DECinfo%output,'(1X,a,i4)')     'FOP Done: Number of basis functions        :', &
           & AtomicFragment%nbasis
      if(.not. DECinfo%onlyvirtpart) then
         write(DECinfo%output,'(1X,a,f16.10)') 'FOP Done: Occupied Fragment energy         :', &
              & AtomicFragment%EoccFOP
      endif
      if(.not. DECinfo%onlyoccpart) then
         write(DECinfo%output,'(1X,a,f16.10)') 'FOP Done: Virtual Fragment energy          :', &
              & AtomicFragment%EvirtFOP
      endif
      if(.NOT.(DECinfo%onlyoccpart.OR. DECinfo%onlyvirtpart))then
         write(DECinfo%output,'(1X,a,f16.10)') 'FOP Done: Lagrangian Fragment energy       :', &
              & AtomicFragment%LagFOP
      end if
      write(DECinfo%output,'(1X,a,i6,a,i6,a,f5.2,a)')  'FOP Done: Occupied reduction removed ', &
           & no_exp-AtomicFragment%noccAOS, ' of ', no_exp, ' orbitals ( ', &
           & (no_exp-AtomicFragment%noccAOS)*100.0_realk/no_exp, ' %)'
      write(DECinfo%output,'(1X,a,i6,a,i6,a,f5.2,a)')  'FOP Done: Virtual  reduction removed ', &
           & nv_exp-AtomicFragment%nunoccAOS, ' of ', nv_exp, ' orbitals ( ', &
           & (nv_exp-AtomicFragment%nunoccAOS)*100.0_realk/nv_exp, ' %)'
      write(DECinfo%output,'(1X,a)') 'FOP========================================================='
      write(DECinfo%output,'(1X,a,/)') 'FOP'


   end subroutine fragment_reduction_procedure


   ! Purpose: Check if the last two steps of the binary search are close enough
   !          for the reduction to be stoped.
   !
   ! Author:  Pablo Baudin
   ! Date:    July 2014
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

         if (reduce_occ) then
            write(DECinfo%output,*) &
            & 'BIN SEARCH: OCCUPIED REDUCTION CONVERGED', gap
         else
            write(DECinfo%output,*) &
            & 'BIN SEARCH: VIRTUAL REDUCTION CONVERGED', gap
         end if

         if (Space2_conv) converged = .true.
      end if

   end subroutine check_red_space_convergence

end module fragment_energy_module
