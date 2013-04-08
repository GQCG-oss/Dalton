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


  ! DEC DEPENDENCIES (within deccc directory)                                                         
  ! ****************************************
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
  use ccsdpt_module , only:ccsdpt_driver,ccsdpt_energy_e4_frag,ccsdpt_energy_e5_frag,&
       ccsdpt_energy_e4_pair,ccsdpt_energy_e5_pair
  use mp2_gradient_module ,only: single_calculate_mp2gradient_driver,&
       & pair_calculate_mp2gradient_driver
  use ccdriver, only: mp2_solver,fragment_ccsolver

public :: optimize_atomic_fragment, pair_driver_singles, atomic_driver, &
     & pair_driver,atomic_driver_advanced,dec_energy_control_center,&
     & Full_DEC_calculation_Lagrangian,estimate_energy_error
private

contains



  !> \brief Construct new atomic fragment based on info in OccAtoms and UnoccAtoms,
  !> and calculate fragment energy. Energy contributions from each individual orbital
  !> is also calculated and stored in MyFragment%OccContribs and MyFragment%VirtContribs for
  !> the occupied and virtual orbitals, respectively.
  !> \author Kasper Kristensen
  !> \date January 2011
  subroutine get_fragment_and_Energy(MyAtom,natoms,OccAtoms,UnoccAtoms,&
       & MyMolecule,MyLsitem,nocc_tot,nunocc_tot,OccOrbitals,UnoccOrbitals,&
       & MyFragment)

    implicit none
    !> Central atom in fragment
    integer, intent(in) :: MyAtom
    !> Number of atoms in the molecule
    integer,intent(in) :: natoms
    !> Which atoms are in the occupied fragment space
    logical, dimension(natoms),intent(in) :: OccAtoms
    !> Which atoms are in the unoccupied fragment space
    logical, dimension(natoms),intent(in) :: UnoccAtoms
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> Full molecule lsitem
    type(lsitem), intent(inout) :: mylsitem
    !> Number of occupied orbitals in full molecule
    integer, intent(in) :: nocc_tot
    !> Number of unoccupied orbitals in full molecule
    integer, intent(in) :: nunocc_tot
    !> Occupied orbitals for full molecule
    type(ccorbital), dimension(nOcc_tot), intent(in) :: OccOrbitals
    !> Unoccupied orbitals for full molecule
    type(ccorbital), dimension(nUnocc_tot), intent(in) :: UnoccOrbitals
    !> Atomic fragment to be determined  (NOT pair fragment)
    type(ccatom), intent(inout) :: MyFragment
    real(realk) :: tcpu, twall
    logical :: DoBasis


    ! **************************************************************
    ! *                       RUN CALCULATION                      *
    ! **************************************************************
    call LSTIMER('START',tcpu,twall,DECinfo%output)



    ! Initialize fragment
    ! *******************
    DoBasis=.true.
    call atomic_fragment_init_atom_specific(MyAtom,natoms,UnoccAtoms, &
         & OccAtoms,nocc_tot,nunocc_tot,OccOrbitals,UnoccOrbitals, &
         & MyMolecule,mylsitem,MyFragment,DoBasis,.false.)

    ! Calculate fragment energies
    ! ***************************
    call single_lagrangian_energy_and_prop(MyFragment)

    call LSTIMER('FRAG: L.ENERGY',tcpu,twall,DECinfo%output)


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
    type(ccorbital), dimension(MyMolecule%numocc), intent(in) :: OccOrbitals
    !> Unoccupied orbitals for full molecule
    type(ccorbital), dimension(MyMolecule%numvirt), intent(in) :: UnoccOrbitals
    !> Which atoms are in the occupied fragment space
    logical, dimension(MyMolecule%numocc),intent(in) :: OccAOS
    !> Which atoms are in the unoccupied fragment space
    logical, dimension(MyMolecule%numvirt),intent(in) :: UnoccAOS
    !> Atomic Fragment to be determined (NOT pair fragment)
    type(ccatom), intent(inout) :: MyFragment
    real(realk) :: tcpu, twall
    logical :: DoBasis


    ! **************************************************************
    ! *                       RUN CALCULATION                      *
    ! **************************************************************
    call LSTIMER('START',tcpu,twall,DECinfo%output)


    ! Init fragment type based on logical vectors
    DoBasis=.true.
    call atomic_fragment_init_orbital_specific(MyAtom,MyMolecule%numvirt, &
         & MyMolecule%numocc, UnoccAOS, &
         & occAOS,OccOrbitals,UnoccOrbitals,MyMolecule,mylsitem,MyFragment,DoBasis,.false.)


    ! Calculate fragment energies
    call single_lagrangian_energy_and_prop(MyFragment)

    call LSTIMER('FRAG: L.ENERGY',tcpu,twall,DECinfo%output)


  end subroutine get_fragment_and_Energy_orb_specific




  !> \brief Wrapper for atomic_driver with the following special features:
  !> 1. It is assumed that the input fragment has been initialized but that
  !> the fragment basis information (expensive box in ccatom type) has not been set.
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
    type(ccatom), intent(inout) :: myfragment
    !> MP2 gradient structure (only calculated if DECinfo%first_order is turned on)
    type(mp2grad),intent(inout) :: grad
    !> Number of occupied orbitals in full molecule
    integer,intent(in) :: nOcc
    !> Number of unoccupied orbitals in full molecule
    integer,intent(in) :: nUnocc
    !> Occupied orbitals for full molecule
    type(ccorbital), intent(in) :: OccOrbitals(nocc)
    !> Unoccupied orbitals for full molecule
    type(ccorbital), intent(in) :: UnoccOrbitals(nunocc)
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
    type(ccorbital), dimension(MyMolecule%numocc), intent(in) :: OccOrbitals
    !> Information about DEC unoccupied orbitals
    type(ccorbital), dimension(MyMolecule%numvirt), intent(in) :: UnoccOrbitals
    !> Atomic fragment 
    type(ccatom), intent(inout) :: myfragment
    !> MP2 gradient structure (only calculated if DECinfo%first_order is turned on)
    type(mp2grad),intent(inout),optional :: grad
    type(ccatom) :: FOfragment

    ! Sanity check
    if(DECinfo%first_order) then
       if(.not. present(grad)) then
          call lsquit('atomic_driver: Gradient argument not present!',-1)
       end if
    end if


    if(DECinfo%fragadapt) then  ! Use fragment-adapted orbitals

       ! Init new fragment with fragment-adapted orbitals
       call init_fragment_adapted(MyMolecule,mylsitem,OccOrbitals,UnoccOrbitals,&
            & MyFragment,FOfragment)

       ! Run calculation using fragment with fragment-adapted orbitals
       if(DECinfo%first_order) then
          call single_lagrangian_energy_and_prop(FOfragment,grad=grad)
       else
          call single_lagrangian_energy_and_prop(FOfragment)
       end if

       ! Copy stuff from FA fragment to original fragment
       call copy_mpi_main_info_from_FOfragment(FOfragment,MyFragment)

       ! Ensure that energies within fragment structure are set consistently
       call get_occ_virt_lag_energies_fragopt(MyFragment)

       call atomic_fragment_free(FOfragment)

    else

       ! Run calculation with input fragment
       if(DECinfo%first_order) then
          call single_lagrangian_energy_and_prop(MyFragment,grad=grad)
       else
          call single_lagrangian_energy_and_prop(MyFragment)
       end if

    end if

  end subroutine atomic_driver



  !> \brief Driver for calculating atomic fragment energy for a given fragment using the Lagrangian approach.
  !> If requested, first order properties (MP2 density or gradient) are also calculated and saved.
  !> \author Kasper Kristensen
  subroutine single_lagrangian_energy_and_prop(MyFragment,grad)

    implicit none
    !> Atomic fragment
    type(ccatom), intent(inout) :: myfragment
    !> MP2 gradient structure (only calculated if DECinfo%first_order is turned on)
    type(mp2grad),intent(inout),optional :: grad
    type(array2) :: t1, ccsdpt_t1
    type(array4) :: VOVO,VOVOocc,VOVOvirt,t2occ,t2virt,VOOO,VOVV,t2,u,VOVOvirtTMP,ccsdpt_t2
    real(realk) :: tcpu, twall

    call LSTIMER('START',tcpu,twall,DECinfo%output)

    ! Which model? MP2,CC2, CCSD etc.
    WhichCCmodel: if(DECinfo%ccModel==1) then ! MP2 calculation

       if(DECinfo%first_order) then  ! calculate also MP2 density integrals
          call MP2_integrals_and_amplitudes(MyFragment,VOVOocc,t2occ,VOVOvirt,t2virt,VOOO,VOVV)
       else ! calculate only MP2 energy integrals and MP2 amplitudes
          call MP2_integrals_and_amplitudes(MyFragment,VOVOocc,t2occ,VOVOvirt,t2virt)
       end if


    else ! higher order CC (currently CC2 or CCSD)


       ! Solve CC equation to calculate amplitudes and integrals 
       ! *******************************************************
       ! Here all output indices in t1,t2, and VOVO are AOS indices.
       call fragment_ccsolver(MyFragment,t1,t2,VOVO)


       ! Extract EOS indices for integrals
       ! *********************************
       call array4_extract_eos_indices_both_schemes(VOVO, &
            & VOVOocc, VOVOvirt, MyFragment)
       call array4_free(VOVO)
       

       ! Calculate combined single+doubles amplitudes
       ! ********************************************
       ! u(a,i,b,j) = t2(a,i,b,j) + t1(a,i)*t1(b,j)          
       call get_combined_SingleDouble_amplitudes(t1,t2,u)


       ! Extract EOS indices for amplitudes
       ! **********************************
       call array4_extract_eos_indices_both_schemes(u, &
            & t2occ, t2virt, MyFragment)
       ! Note, t2occ and t2virt also contain singles contributions
       call array4_free(u)


       ! calculate ccsd(t) fragment energies
       ! ***********************************

       ! we calculate (T) contribution to single  fragment energy 
       ! and store in MyFragment%energies(8) and MyFragment%energies(9)

       if(DECinfo%ccModel==4) then

          ! init ccsd(t) singles and ccsd(t) doubles (*T1 and *T2)
          ccsdpt_t1 = array2_init([MyFragment%nunoccAOS,MyFragment%noccAOS])
          ccsdpt_t2 = array4_init([MyFragment%nunoccAOS,MyFragment%nunoccAOS,&
               &MyFragment%noccAOS,MyFragment%noccAOS])

          ! call ccsd(t) driver and single fragment evaluation
          call ccsdpt_driver(MyFragment%noccAOS,MyFragment%nunoccAOS,&
                             & MyFragment%number_basis,MyFragment%ppfock,&
                             & MyFragment%qqfock,MyFragment%ypo,&
                             & MyFragment%ypv,MyFragment%mylsitem,&
                             & t2,ccsdpt_t1,ccsdpt_t2)
          call ccsdpt_energy_e4_frag(MyFragment,t2,ccsdpt_t2,&
                             & MyFragment%OccContribs,MyFragment%VirtContribs)
          call ccsdpt_energy_e5_frag(MyFragment,t1,ccsdpt_t1)

          ! release ccsd(t) singles and doubles amplitudes
          call array2_free(ccsdpt_t1)
          call array4_free(ccsdpt_t2)

       end if

       call array2_free(t1)
       call array4_free(t2)

    end if WhichCCmodel


    ! Calcuate atomic fragment energy using Lagrangian scheme
    ! *******************************************************

    ! MODIFY FOR NEW MODEL!
    ! If you implement a new model, which does fit into the standard CC energy
    ! expression implemented in get_atomic_fragment_energy, then please make a new atomic
    ! fragment energy subroutine and call it from here instead of 
    ! calling get_atomic_fragment_energy.
    ! If the energy in your new model does fit into the standard CC energy expression,
    ! then grep for
    ! "MODIFY FOR NEW MODEL THAT FITS INTO STANDARD CC ENERGY EXPRESSION"
    ! The calculated pair interaction energy should be saved in myfragment%energies(?),
    ! see "energies" in ccatom type definition to determine the "?".

    ! For frozen core and first order properties we need to remove core indices from VOVOvirt, 
    ! since they, "rather incoveniently", are required for the gradient but not for the energy
    if(DECinfo%frozencore .and. DECinfo%first_order) then
       call remove_core_orbitals_from_last_index(MyFragment,VOVOvirt,VOVOvirtTMP)
       call get_atomic_fragment_energy(VOVOocc,VOVOvirtTMP,t2occ,t2virt,MyFragment)
       call array4_free(VOVOvirtTMP)
    else ! use VOVOvirt as it is
       call get_atomic_fragment_energy(VOVOocc,VOVOvirt,t2occ,t2virt,MyFragment)
    end if

    ! MODIFY FOR NEW CORRECTION!
    ! If you implement a new correction (e.g. F12) please insert call to subroutine here
    ! which calculates atomic fragment contribution and saves it in myfragment%energies(?),
    ! see dec_readme file.


    call LSTIMER('SINGLE L.ENERGY',tcpu,twall,DECinfo%output)


    ! First order properties
    ! **********************
    if(DECinfo%first_order) then
       call single_calculate_mp2gradient_driver(MyFragment,t2occ,t2virt,VOOO,VOVV,VOVOocc,VOVOvirt,grad)
       call array4_free(VOOO)
       call array4_free(VOVV)
    end if

    ! Free remaining arrays
    call array4_free(VOVOocc)
    call array4_free(VOVOvirt)
    call array4_free(t2occ)
    call array4_free(t2virt)


  end subroutine single_lagrangian_energy_and_prop



  !> \brief Contract amplitudes, multipliers, and integrals to calculate atomic fragment Lagrangian energy.
  !> \author Kasper Kristensen
  !> \date August 2011
  subroutine get_atomic_fragment_energy(gocc,gvirt,t2occ,t2virt,MyFragment)

    implicit none
    !> Two-electron integrals (a i | b j), only occ orbitals on central atom, virt AOS orbitals
    type(array4), intent(in) :: gocc
    !> Two-electron integrals (a i | b j), only virt orbitals on central atom, occ AOS orbitals
    type(array4), intent(in) :: gvirt
    !> MP2 amplitudes, only occ orbitals on central atom, virt AOS orbitals
    type(array4), intent(in) :: t2occ
    !> MP2 amplitudes, only virt orbitals on central atom, occ AOS orbitals
    type(array4), intent(in) :: t2virt
    !> Atomic fragment 
    type(ccatom), intent(inout) :: myfragment
    integer :: noccEOS,nvirtEOS,noccAOS,nvirtAOS
    integer :: i,j,k,a,b,c
    real(realk) :: tcpu1, twall1, tcpu2,twall2, tcpu,twall
    real(realk) :: e1, e2, e3, e4,tmp,multaibj
    logical ::  something_wrong
    real(realk) :: e1_final, e2_final,e3_final,e4_final
    real(realk),pointer :: occ_tmp(:),virt_tmp(:),VirtMat_tmp(:,:),OccMat_tmp(:,:)

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
    noccEOS = MyFragment%noccEOS
    nvirtEOS = MyFragment%nunoccEOS
    noccAOS = MyFragment%noccAOS
    nvirtAOS = MyFragment%nunoccAOS
    e1_final=0E0_realk
    e2_final=0E0_realk
    e3_final=0E0_realk
    e4_final=0E0_realk
    ! Just in case, zero individual orbital contributions for fragment
    MyFragment%OccContribs=0E0_realk
    MyFragment%VirtContribs=0E0_realk


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

    ! Init correlation density matrices
    if(.not. MyFragment%CDset) then
       call mem_alloc(MyFragment%OccMat,MyFragment%noccAOS,MyFragment%noccAOS)
       MyFragment%OccMat=0.0_realk
       call mem_alloc(MyFragment%VirtMat,MyFragment%nunoccAOS,MyFragment%nunoccAOS)
       MyFragment%VirtMat=0.0_realk
       MyFragment%CDSet=.true. ! correlation density matrices have been set
    end if



call mem_TurnONThread_Memory()
!$OMP PARALLEL DEFAULT(shared) PRIVATE(multaibj,tmp,e1,e2,j,b,i,a,c,virt_tmp,VirtMat_tmp)
call init_threadmemvar()
e1=0E0_realk
e2=0E0_realk
! Contributions from each individual virtual orbital
call mem_alloc(virt_tmp,nvirtAOS)
virt_tmp = 0.0E0_realk

! Virtual correlation density matrix
call mem_alloc(VirtMat_tmp,nvirtAOS,nvirtAOS)
VirtMat_tmp = 0.0_realk

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
                tmp = t2occ%val(a,i,b,j)*(2.0_realk*gocc%val(a,i,b,j) - gocc%val(b,i,a,j))

                ! Update total atomic fragment energy contribution 1
                e1 = e1 + tmp

                ! Update contribution from orbital a
                virt_tmp(a) = virt_tmp(a) + tmp
                ! Update contribution from orbital b (only if different from a to avoid double counting)
                if(a/=b) virt_tmp(b) = virt_tmp(b) + tmp


                ! Contribution 2
                ! --------------

                ! Skip contribution 2 for hybrid scheme
                if(.not. DECinfo%HybridScheme) then
                   ! Multiplier (multiplied by one half)
                   multaibj = 2.0_realk*t2occ%val(a,i,b,j) - t2occ%val(b,i,a,j)

                   do c=1,nvirtAOS

                      ! Energy contribution for orbitals (j,b,i,a,c)
                      tmp = t2occ%val(c,i,b,j)*MyFragment%qqfock(c,a) + t2occ%val(a,i,c,j)*MyFragment%qqfock(c,b)
                      tmp = multaibj*tmp

                      ! Update total atomic fragment energy contribution 2
                      e2 = e2 + tmp

                      ! Update contribution from orbital a
                      virt_tmp(a) = virt_tmp(a) + tmp
                      ! Update contribution from orbital b (only if different from a to avoid double counting)
                      if(a/=b) virt_tmp(b) = virt_tmp(b) + tmp
                      ! Update contribution from orbital c (only if different from a and b)
                      if( (a/=c) .and. (b/=c) ) virt_tmp(c) = virt_tmp(c) + tmp

                      ! Virtual correlation density matrix
                      VirtMat_tmp(a,b) = VirtMat_tmp(a,b) + t2occ%val(a,i,c,j)&
                           &*(4.0_realk*t2occ%val(b,i,c,j) - 2.0_realk*t2occ%val(b,j,c,i))

                   end do

                end if


             end do
          end do
       end do
    end do

!$OMP END DO NOWAIT

! Total e1, e2, and individual virt atomic contributions are found by summing all thread contributions
!$OMP CRITICAL
e1_final = e1_final + e1
e2_final = e2_final + e2

! Update total virtual contributions to fragment energy
do a=1,nvirtAOS
   MyFragment%VirtContribs(a) =MyFragment%VirtContribs(a) + virt_tmp(a)

   ! Virtual correlation density matrix
   do b=1,nvirtAOS
      MyFragment%VirtMat(b,a) = MyFragment%VirtMat(b,a) + VirtMat_tmp(b,a)
   end do
end do

!$OMP END CRITICAL

call mem_dealloc(virt_tmp)
call mem_dealloc(VirtMat_tmp)
call collect_thread_memory()
!$OMP END PARALLEL
call mem_TurnOffThread_Memory()




    ! Calculate e3 and e4
    ! *******************

call mem_TurnONThread_Memory()
!$OMP PARALLEL DEFAULT(shared) PRIVATE(multaibj,tmp,e3,e4,j,b,i,a,k,occ_tmp,OccMat_tmp)
call init_threadmemvar()
e3=0E0_realk
e4=0E0_realk

! Contributions from each individual occupied orbital
call mem_alloc(occ_tmp,noccAOS)
occ_tmp = 0.0E0_realk

! Occupied correlation density matrix
call mem_alloc(OccMat_tmp,noccAOS,noccAOS)
OccMat_tmp = 0.0_realk

!$OMP DO SCHEDULE(dynamic,1)

    do j=1,noccAOS
       do b=1,nvirtEOS
          do i=1,noccAOS
             do a=1,nvirtEOS


                ! Contribution 3
                ! --------------

                ! Multiplier (multiplied by one half)
                multaibj = 2.0_realk*t2virt%val(a,i,b,j) - t2virt%val(b,i,a,j)


                ! Energy contribution for orbitals (j,b,i,a)
                tmp = multaibj*gvirt%val(a,i,b,j)

                ! Update total atomic fragment energy contribution 3
                e3 = e3 + tmp
                
                ! Update contribution from orbital i
                occ_tmp(i) = occ_tmp(i) + tmp
                ! Update contribution from orbital j (only if different from i to avoid double counting)
                if(i/=j) occ_tmp(j) = occ_tmp(j) + tmp


                ! Contribution 4
                ! --------------

                ! Skip contribution 4 for hybrid scheme
                if(.not. DECinfo%HybridScheme) then

                   do k=1,noccAOS

                      tmp =  t2virt%val(a,k,b,j)*MyFragment%ppfock(k,i) &
                           & + t2virt%val(a,i,b,k)*MyFragment%ppfock(k,j)
                      tmp = -multaibj*tmp

                      ! Update total atomic fragment energy contribution 4
                      e4 = e4 + tmp

                      ! Update contribution from orbital i
                      occ_tmp(i) = occ_tmp(i) + tmp
                      ! Update contribution from orbital j (only if different from i to avoid double counting)
                      if(i/=j) occ_tmp(j) = occ_tmp(j) + tmp
                      ! Update contribution from orbital k (only if different from i and j)
                      if( (i/=k) .and. (j/=k) ) occ_tmp(k) = occ_tmp(k) + tmp

                      ! Occupied correlation density matrix
                      OccMat_tmp(i,j) = OccMat_tmp(i,j) + t2virt%val(a,i,b,k)&
                           &*(4.0_realk*t2virt%val(a,j,b,k) - 2.0_realk*t2virt%val(a,k,b,j))

                   end do

                end if

             end do
          end do
       end do
    end do

!$OMP END DO NOWAIT

! Total e3, e4, and individual occ atomic contributions are found by summing all thread contributions
!$OMP CRITICAL
e3_final = e3_final + e3
e4_final = e4_final + e4

! Update total occupied contributions to fragment energy
do i=1,noccAOS
MyFragment%OccContribs(i) = MyFragment%OccContribs(i) + occ_tmp(i)

! Occupied correlation density matrix
do j=1,noccAOS
   MyFragment%OccMat(j,i) = MyFragment%OccMat(j,i) + OccMat_tmp(j,i)
end do

end do
!$OMP END CRITICAL

call mem_dealloc(occmat_tmp)
call mem_dealloc(occ_tmp)
call collect_thread_memory()
!$OMP END PARALLEL
call mem_TurnOffThread_Memory()



    ! Total atomic fragment energy
    ! ****************************
    ! MODIFY FOR NEW MODEL THAT FITS INTO STANDARD CC ENERGY EXPRESSION
    select case(DECinfo%ccmodel)
    case(1)
       ! MP2
       MyFragment%energies(1) = e1_final + e2_final + e3_final + e4_final  ! Lagrangian
       MyFragment%energies(2) = e1_final   ! occupied
       MyFragment%energies(3) = e3_final   ! virtual
    case(2)
       ! CC2
       MyFragment%energies(4) = e1_final   ! occupied
       MyFragment%energies(5) = e3_final   ! virtual
    case(3)
       ! CCSD
       MyFragment%energies(6) = e1_final   ! occupied
       MyFragment%energies(7) = e3_final   ! virtual
    case(4)
       ! Save also CCSD contribution for CCSD(T)
       MyFragment%energies(6) = e1_final   ! occupied
       MyFragment%energies(7) = e3_final   ! virtual
    end select
    ! Energy contributions other than MP2,CC2, and CCSD are calculated elsewhere

    ! Set energies used by fragment optimization
    call get_occ_virt_lag_energies_fragopt(MyFragment)

    ! Print out contributions
    ! ***********************

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) '**********************************************************************'
    write(DECinfo%output,'(1X,a,i7)') 'Energy summary for fragment: ', &
         & MyFragment%atomic_number
    write(DECinfo%output,*) '**********************************************************************'
    write(DECinfo%output,'(1X,a,g20.10)') 'Single occupied energy = ', e1_final
    write(DECinfo%output,'(1X,a,g20.10)') 'Single virtual  energy = ', e3_final
    write(DECinfo%output,'(1X,a,g20.10)') 'Single Lagrangian occ term  = ', e2_final
    write(DECinfo%output,'(1X,a,g20.10)') 'Single Lagrangian virt term = ', e4_final

    write(DECinfo%output,*)
    write(DECinfo%output,*)

    call LSTIMER('START',tcpu2,twall2,DECinfo%output)
    DECinfo%energy_time_wall = DECinfo%energy_time_wall + (twall2-twall1)
    DECinfo%energy_time_cpu = DECinfo%energy_time_cpu + (tcpu2-tcpu1)
    call LSTIMER('L.ENERGY CONTR',tcpu,twall,DECinfo%output)



  end subroutine get_atomic_fragment_energy



  !> \brief Wrapper for pair_driver where can attach existing full
  !> molecular singles amplitudes to the fragment structure and update new
  !> improved full molecular singles amplitudes by the calculated fragment singles amplitudes.
  !> \author Kasper Kristensen
  !> \date March 2012
  subroutine pair_driver_singles(natoms,nocc,nunocc,distancetable,&
       & OccOrbitals,UnoccOrbitals,MyLsitem,MyMolecule,&
       & Fragment1,Fragment2,PairFragment,t1old,t1new)

    implicit none
    !> Number of atoms for full molecule
    integer, intent(in) :: natoms
    real(realk), intent(in) :: DistanceTable(natoms,natoms)
    !> Number of occupied orbitals in full molecule
    integer,intent(in) :: nOcc
    !> Number of unoccupied orbitals in full molecule
    integer,intent(in) :: nUnocc
    !> Occupied orbitals for full molecule
    type(ccorbital), intent(in) :: OccOrbitals(nocc)
    !> Unoccupied orbitals for full molecule
    type(ccorbital), intent(in) :: UnoccOrbitals(nunocc)
    !> Full molecule info (will be unchanged at output)
    type(fullmolecule), intent(inout) :: MyMolecule
    !> LSDalton info
    type(lsitem), intent(inout) :: mylsitem
    !> Fragment 1 in the pair fragment
    type(ccatom),intent(inout) :: Fragment1
    !> Fragment 2 in the pair fragment
    type(ccatom),intent(inout) :: Fragment2
    !> Pair fragment
    type(ccatom), intent(inout) :: PairFragment
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
    call pair_driver(MyMolecule,mylsitem,OccOrbitals,UnoccOrbitals, &
         & Fragment1,Fragment2,natoms,distancetable, PairFragment,grad)

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
  subroutine pair_driver(MyMolecule,mylsitem,OccOrbitals,UnoccOrbitals,&
       & Fragment1,Fragment2,natoms,distancetable,PairFragment,grad)

    implicit none
    !> Full molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> LS item info
    type(lsitem), intent(inout) :: mylsitem
    !> Information about DEC occupied orbitals
    type(ccorbital), dimension(MyMolecule%numocc), intent(in) :: OccOrbitals
    !> Information about DEC unoccupied orbitals
    type(ccorbital), dimension(MyMolecule%numvirt), intent(in) :: UnoccOrbitals
    !> Fragment 1 in the pair fragment
    type(ccatom),intent(inout) :: Fragment1
    !> Fragment 2 in the pair fragment
    type(ccatom),intent(inout) :: Fragment2
    !> Number of atoms for full molecule
    integer, intent(in) :: natoms
    real(realk), intent(in) :: DistanceTable(natoms,natoms)
    !> Atomic fragment 
    type(ccatom), intent(inout) :: PairFragment
    !> MP2 gradient structure (only calculated if DECinfo%first_order is turned on)
    type(mp2grad),intent(inout) :: grad
    type(ccatom) :: FOfragment


    if(DECinfo%fragadapt) then
       ! Init new fragment with fragment-adapted orbitals'
       call init_fragment_adapted(MyMolecule,mylsitem,OccOrbitals,UnoccOrbitals,&
            & PairFragment,FOfragment)

       ! Run calculation using fragment with fragment-adapted orbitals
       call pair_lagrangian_energy_and_prop(Fragment1,Fragment2, &
            & natoms, DistanceTable, FOfragment,grad)

       ! Copy stuff from FO fragment to original fragment
       call copy_mpi_main_info_from_FOfragment(FOfragment,PairFragment)

       call atomic_fragment_free(FOfragment)
    else

       ! Run calculation using input fragment
       call pair_lagrangian_energy_and_prop(Fragment1,Fragment2, &
            & natoms, DistanceTable, PairFragment,grad)       
    end if

  end subroutine pair_driver



  !> \brief Driver for calculating pair interaction energy for a given
  !> pair fragment using the Lagrangian approach.
  !> \author Kasper Kristensen
  subroutine pair_lagrangian_energy_and_prop(Fragment1,Fragment2, &
       & natoms, DistanceTable, PairFragment,grad)

    implicit none
    !> Fragment 1 in the pair fragment
    type(ccatom),intent(inout) :: Fragment1
    !> Fragment 2 in the pair fragment
    type(ccatom),intent(inout) :: Fragment2
    !> Number of atoms for full molecule
    integer, intent(in) :: natoms
    !> Distance table for all atoms in the molecule
    real(realk), intent(in) :: DistanceTable(natoms,natoms)
    !> Pair fragment formed from fragment 1 and 2 
    type(ccatom), intent(inout) :: PairFragment
    !> MP2 gradient structure (only calculated if DECinfo%first_order is turned on)
    type(mp2grad),intent(inout) :: grad
    type(array2) :: t1, ccsdpt_t1
    type(array4) :: g,VOVOocc,VOVOvirt,t2occ,t2virt,VOOO,VOVV,t2,u,VOVO,ccsdpt_t2,VOVOvirtTMP
    real(realk) :: tcpu, twall
    real(realk) :: tmp_energy

    call LSTIMER('START',tcpu,twall,DECinfo%output)


    WhichCCModel: if(DECinfo%ccModel==1) then ! MP2

       if(DECinfo%first_order) then  ! calculate also MP2 density integrals
          call MP2_integrals_and_amplitudes(PairFragment,VOVOocc,t2occ,VOVOvirt,t2virt,VOOO,VOVV)
       else ! calculate only MP2 energy integrals and MP2 amplitudes
          call MP2_integrals_and_amplitudes(PairFragment,VOVOocc,t2occ,VOVOvirt,t2virt)
       end if

    else ! higher order CC (currently CC2 or CCSD)


       ! Solve CC equation to calculate amplitudes and integrals 
       ! *******************************************************
       ! Here all output indices in t1,t2, and VOVO are AOS indices. 
       call fragment_ccsolver(PairFragment,t1,t2,VOVO)


       ! Extract EOS indices for integrals 
       ! *********************************
       call array4_extract_eos_indices_both_schemes(VOVO, &
            & VOVOocc, VOVOvirt, PairFragment)
       call array4_free(VOVO)


       ! Calculate combined single+doubles amplitudes
       ! ********************************************
       ! u(a,i,b,j) = t2(a,i,b,j) + t1(a,i)*t1(b,j) 
       call get_combined_SingleDouble_amplitudes(t1,t2,u)


       ! Extract EOS indices for amplitudes
       ! **********************************
       call array4_extract_eos_indices_both_schemes(u, &
            & t2occ, t2virt, PairFragment)
       ! Note, t2occ and t2virt also contain singles contributions
       call array4_free(u)

    end if WhichCCmodel


    ! Calculate pair interaction energy using Lagrangian scheme
    ! *********************************************************


    ! MODIFY FOR NEW MODEL!
    ! If you implement a new model, which does fit into the standard CC energy
    ! expression implemented in get_pair_fragment_energy, then please make a new pair
    ! fragment energy subroutine and call it from here instead of 
    ! calling get_pair_fragment_energy.
    ! If the energy in your new model does fit into the standard CC energy expression,
    ! then grep for
    ! "MODIFY FOR NEW MODEL THAT FITS INTO STANDARD CC ENERGY EXPRESSION"
    ! The calculated pair interaction energy should be saved in pairfragment%energies(?),
    ! see "energies" in ccatom type definition to determine the "?".

    ! For frozen core and first order properties we need to remove core indices from VOVOvirt, 
    ! since they, "rather incoveniently", are required for the gradient but not for the energy
    if(DECinfo%frozencore .and. DECinfo%first_order) then
       call remove_core_orbitals_from_last_index(PairFragment,VOVOvirt,VOVOvirtTMP)
       call get_pair_fragment_energy(VOVOocc,VOVOvirtTMP,t2occ,t2virt,&
            & Fragment1, Fragment2, PairFragment,natoms,DistanceTable)
       call array4_free(VOVOvirtTMP)
    else ! use VOVOvirt as it is
       call get_pair_fragment_energy(VOVOocc,VOVOvirt,t2occ,t2virt,&
            & Fragment1, Fragment2, PairFragment,natoms,DistanceTable)
    end if

    ! MODIFY FOR NEW CORRECTION!
    ! If you implement a new correction (e.g. F12) please insert call to subroutine here
    ! which calculates pair fragment contribution and saves it in pairfragment%energies(?),
    ! see dec_readme file.

    call LSTIMER('PAIR L.ENERGY',tcpu,twall,DECinfo%output)


    ! First order properties
    ! **********************
    if(DECinfo%first_order) then
       call pair_calculate_mp2gradient_driver(Fragment1,Fragment2,PairFragment,&
            & t2occ,t2virt,VOOO,VOVV,VOVOocc,VOVOvirt,grad)
       call array4_free(VOOO)
       call array4_free(VOVV)
    end if
 

    ! calculate ccsd(t) pair interaction energies
    ! *******************************************

    ! (T) energy contributions are stored 
    ! in PairFragment%energies(8) and PairFragment%energies(9) 

    if (DECinfo%CCModel .eq. 4) then

       ! init ccsd(t) singles and ccsd(t) doubles
       ccsdpt_t1 = array2_init([PairFragment%nunoccAOS,PairFragment%noccAOS])
       ccsdpt_t2 = array4_init([PairFragment%nunoccAOS,PairFragment%nunoccAOS,&
            & PairFragment%noccAOS,PairFragment%noccAOS])

       ! call ccsd(t) driver and pair fragment evaluation
       call ccsdpt_driver(PairFragment%noccAOS,PairFragment%nunoccAOS,&
                          & PairFragment%number_basis,PairFragment%ppfock,&
                          & PairFragment%qqfock,PairFragment%ypo,&
                          & PairFragment%ypv,PairFragment%mylsitem,&
                          & t2,ccsdpt_t1,ccsdpt_t2)
       call ccsdpt_energy_e4_pair(Fragment1,Fragment2,&
                          &PairFragment,t2,ccsdpt_t2)
       call ccsdpt_energy_e5_pair(Fragment1,Fragment2,&
                          &PairFragment,t1,ccsdpt_t1)

       ! release ccsd(t) singles and doubles amplitudes
       call array2_free(ccsdpt_t1)
       call array4_free(ccsdpt_t2)

    end if


    if(DECinfo%ccmodel/=1) then
       call array2_free(t1)
       call array4_free(t2)
    end if

    ! Free remaining arrays
    call array4_free(VOVOocc)
    call array4_free(VOVOvirt)
    call array4_free(t2occ)
    call array4_free(t2virt)

  end subroutine pair_lagrangian_energy_and_prop



  !> \brief Contract amplitudes, multipliers, and integrals to calculate pair interaction
  !> Lagrangian energy. Heavily inspired by get_atomic_fragment_energy, but it is necessary
  !> to keep it separate for clairity.
  !> \author Kasper Kristensen
  !> \date August 2011
  subroutine get_pair_fragment_energy(gocc,gvirt,t2occ,t2virt,&
       & Fragment1, Fragment2, PairFragment,natoms,DistanceTable)


    implicit none
    !> Two-electron integrals (a i | b j), only occ orbitals on central atom, virt AOS orbitals
    type(array4), intent(in) :: gocc
    !> Two-electron integrals (a i | b j), only virt orbitals on central atom, occ AOS orbitals
    type(array4), intent(in) :: gvirt
    !> MP2 amplitudes, only occ orbitals on central atom, virt AOS orbitals
    type(array4), intent(in) :: t2occ
    !> MP2 amplitudes, only virt orbitals on central atom, occ AOS orbitals
    type(array4), intent(in) :: t2virt
    !> Fragment 1 in the pair fragment
    type(ccatom),intent(inout) :: Fragment1
    !> Fragment 2 in the pair fragment
    type(ccatom),intent(inout) :: Fragment2
    !> Pair fragment formed from fragment 1 and 2
    type(ccatom), intent(inout) :: PairFragment
    !> Number of atoms for full molecule
    integer, intent(in) :: natoms
    !> Distance table for all atoms in the molecule
    real(realk), intent(in) :: DistanceTable(natoms,natoms)
    integer :: noccEOS,nvirtEOS,noccAOS,nvirtAOS
    integer :: i,j,k,a,b,c
    real(realk) :: tcpu, twall,pairdist,au_to_angstrom
    real(realk) :: e1, e2, e3, e4,tmp,multaibj
    real(realk) :: tcpu1,tcpu2,twall1,twall2
    logical,pointer :: dopair_occ(:,:), dopair_virt(:,:)
    real(realk) :: e1_final, e2_final,e3_final,e4_final
    logical :: something_wrong


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
    e1_final=0E0_realk
    e2_final=0E0_realk
    e3_final=0E0_realk
    e4_final=0E0_realk
    ! Distance between fragments in Angstrom
    pairdist = get_distance_between_fragments(Fragment1,Fragment2,natoms,DistanceTable)
    au_to_angstrom = bohr_to_angstrom
    pairdist = au_to_angstrom*pairdist

    ! Which "interaction pairs" to include for occ and unocc space (avoid double counting)
    call mem_alloc(dopair_occ,noccEOS,noccEOS)
    call mem_alloc(dopair_virt,nvirtEOS,nvirtEOS)
    call which_pairs_occ(Fragment1,Fragment2,PairFragment,dopair_occ)
    call which_pairs_unocc(Fragment1,Fragment2,PairFragment,dopair_virt)


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
       call lsquit('get_pair_fragment_energy: &
            & Input dimensions do not match!',-1)
    end if




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
                   e1 = e1 + t2occ%val(a,i,b,j)*(2.0_realk*gocc%val(a,i,b,j) - gocc%val(b,i,a,j))


                ! Skip contribution 2 for hybrid scheme
                if(.not. DECinfo%HybridScheme) then

                   ! Multiplier (multiplied by one half)
                   multaibj = 2.0_realk*t2occ%val(a,i,b,j) - t2occ%val(b,i,a,j)

                   tmp = 0E0_realk
                   do c=1,nvirtAOS
                      tmp = tmp + t2occ%val(c,i,b,j)*PairFragment%qqfock(c,a) &
                           & + t2occ%val(a,i,c,j)*PairFragment%qqfock(c,b)
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
e1_final = e1_final + e1
e2_final = e2_final + e2
!$OMP END CRITICAL

call collect_thread_memory()
!$OMP END PARALLEL
call mem_TurnOffThread_Memory()





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
                   multaibj = 2.0_realk*t2virt%val(a,i,b,j) - t2virt%val(b,i,a,j)

                   ! Update total atomic fragment energy contribution 3
                   e3 = e3 + multaibj*gvirt%val(a,i,b,j)


                ! Skip contribution 4 for hybrid scheme
                if(.not. DECinfo%HybridScheme) then

                   tmp=0E0_realk
                   do k=1,noccAOS
                      tmp = tmp + t2virt%val(a,k,b,j)*PairFragment%ppfock(k,i) &
                           & + t2virt%val(a,i,b,k)*PairFragment%ppfock(k,j)
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
e3_final = e3_final + e3
e4_final = e4_final + e4
!$OMP END CRITICAL

call collect_thread_memory()
!$OMP END PARALLEL
call mem_TurnOffThread_Memory()



    ! Total pair interaction energy
    ! *****************************
    ! MODIFY FOR NEW MODEL THAT FITS INTO STANDARD CC ENERGY EXPRESSION
    select case(DECinfo%ccmodel)
    case(1)
       ! MP2
       PairFragment%energies(1) = e1_final + e2_final + e3_final + e4_final  ! Lagrangian
       PairFragment%energies(2) = e1_final   ! occupied
       PairFragment%energies(3) = e3_final   ! virtual
    case(2)
       ! CC2
       PairFragment%energies(4) = e1_final   ! occupied
       PairFragment%energies(5) = e3_final   ! virtual
    case(3)
       ! CCSD
       PairFragment%energies(6) = e1_final   ! occupied
       PairFragment%energies(7) = e3_final   ! virtual
    case(4)
       ! save CCSD contribution for CCSD(T)
       PairFragment%energies(6) = e1_final   ! occupied
       PairFragment%energies(7) = e3_final   ! virtual
    end select
    ! Energy contributions other than MP2,CC2, and CCSD are calculated elsewhere
    call mem_dealloc(dopair_occ)
    call mem_dealloc(dopair_virt)


    ! Print out contributions
    ! ***********************

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) '*****************************************************************************'
    write(DECinfo%output,'(1X,a,2i7)') 'Energy summary for pair fragment: ', &
         & Fragment1%atomic_number, Fragment2%atomic_number
    write(DECinfo%output,*) '*****************************************************************************'

    write(DECinfo%output,'(1X,a,g20.10)') 'Pair occupied energy = ', e1_final
    write(DECinfo%output,'(1X,a,g20.10)') 'Pair virtual  energy = ', e3_final
    write(DECinfo%output,'(1X,a,g20.10)') 'Pair Lagrangian occ term  = ', e2_final
    write(DECinfo%output,'(1X,a,g20.10)') 'Pair Lagrangian virt term = ', e4_final


    write(DECinfo%output,'(1X,a,g16.5,g20.10)') 'Distance(Ang), pair occ energy  = ', pairdist,e1_final
    write(DECinfo%output,'(1X,a,g16.5,g20.10)') 'Distance(Ang), pair virt energy = ', pairdist,e3_final
    write(DECinfo%output,'(1X,a,g16.5,g20.10)') 'Distance(Ang), pair lagr. occ term  = ', pairdist,e2_final
    write(DECinfo%output,'(1X,a,g16.5,g20.10)') 'Distance(Ang), pair lagr. virt term = ', pairdist,e4_final
    write(DECinfo%output,*)
    write(DECinfo%output,*)


    call LSTIMER('L.ENERGY CONTR',tcpu,twall,DECinfo%output)
    call LSTIMER('START',tcpu2,twall2,DECinfo%output)
    DECinfo%energy_time_wall = DECinfo%energy_time_wall + (twall2-twall1)
    DECinfo%energy_time_cpu = DECinfo%energy_time_cpu + (tcpu2-tcpu1)



  end subroutine get_pair_fragment_energy



  !> \brief Full DEC calculation where exact single and pair energies are calculated using Lagrangian approach.
  !> Only implemented for MP2.
  !> \author Kasper Kristensen
  !> \date April 2011
  subroutine Full_DEC_calculation_Lagrangian(MyMolecule,mylsitem,OccOrbitals,UnoccOrbitals, &
       & natoms,nocctot,nunocc, DistanceTable,Ecorr)

    implicit none
    !> Number of atoms in full molecule
    integer,intent(in) :: natoms
    !> Number of occupied orbitals in molecule
    integer,intent(in) :: nocctot
    !> Number of unoccupied orbitals in molecule
    integer,intent(in) :: nunocc
    !> Molecule info
    type(fullmolecule), intent(in) :: MyMolecule
    !> LS Dalton info
    type(lsitem), intent(inout) :: mylsitem
    !> Occupied MOs
    type(ccorbital), intent(in) :: OccOrbitals(nocctot)
    !> UnOccupied MOs
    type(ccorbital), intent(in) :: UnoccOrbitals(nunocc)
    !> Distance table with interatomic distances for atoms in molecule
    real(realk),intent(in) :: DistanceTable(natoms,natoms)
    !> Total correlation energy
    real(realk),intent(inout) :: Ecorr
    logical,dimension(natoms) :: orbitals_assigned
    type(array2) :: Cocc, Cvirt
    type(array4) :: t2, g
    real(realk) :: energy_matrix(natoms,natoms), multaibj, multbiaj
    integer :: nthreads, idx, nbatchINT, intstep, nbasis,ncore,offset
    integer :: i,j,k,a,b,c,atomI,atomJ,atomA,atomB,nocc
    real(realk) :: intMEM, solMEM,OO,VV,AA,BB,mem_required
    real(realk) :: singleenergy, pairenergy, tmp, tmp2
    real(realk),dimension(natoms,natoms) :: e1,e2,e3,e4,e1_tmp,e2_tmp,e3_tmp,e4_tmp
    integer, dimension(4) :: dims
    integer, dimension(2) :: occ_dims,virt_dims
    real(realk), pointer :: gval(:,:,:),t2val(:,:,:),ppfock(:,:)
#ifdef VAR_OMP
    integer, external :: OMP_GET_MAX_THREADS
#endif

    ! Only for MP2
    if(DECinfo%ccModel/=1) then
       call lsquit('Full_DEC_calculation: Only implemented for MP2!', DECinfo%output)
    end if

    if(DECinfo%frozencore) then
       ! Frozen core: Only valence orbitals
       nocc = MyMolecule%nval
    else
       nocc = MyMolecule%numocc
    end if


    ! Initialize stuff
    ! ****************
    ncore = MyMolecule%ncore
    nbasis=MyMolecule%nbasis
    energy_matrix(:,:) = 0E0_realk
    dims = [nunocc, nocc, nunocc, nocc]
    occ_dims = [nbasis,nocc]
    virt_dims = [nbasis,nunocc]
    call mem_alloc(ppfock,nocc,nocc)
    if(DECinfo%frozencore) then
       ! Only copy valence orbitals into array2 structure
       Cocc=array2_init(occ_dims)
       do i=1,nocc
          Cocc%val(:,i) = MyMolecule%ypo(:,i+Ncore)
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
       Cocc=array2_init(occ_dims,MyMolecule%ypo)
       ppfock = MyMolecule%ppfock
       offset=0
    end if
    Cvirt=array2_init(virt_dims,MyMolecule%ypv)
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
    call get_VOVO_integrals(mylsitem,nbasis,nocc,nunocc,Cvirt,Cocc,g)
    call array2_free(Cocc)
    call array2_free(Cvirt)

    ! Get t2 amplitudes
    ! *****************
    call mp2_solver(nocc,nunocc,ppfock,MyMolecule%qqfock,g,t2)



    ! STATUS: Now integrals (g) and amplitudes (t) have been determined
    ! for the full molecular system. The individual fragment contributions - solved in the
    ! total orbital space - can now be determined by extracting EOS indices from the full indicies.



    ! Energy contributions (see get_atomic_fragment_energy_and_prop subroutine)
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
    do atomI=1,natoms
       do atomJ=1,natoms
          energy_matrix(atomI,atomJ) = e1(atomI,atomJ) &
               & + e2(atomI,atomJ) + e3(atomI,atomJ) + e4(atomI,atomJ)
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


    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a)') '*****************************************************************************'
    write(DECinfo%output,'(1X,a)') '*               Full DEC Lagrangian calculation is done!                    *'
    write(DECinfo%output,'(1X,a)') '*****************************************************************************'
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(8X,a)') '-- Atomic fragments'
    write(DECinfo%output,'(8X,a)') '------    --------------------'
    write(DECinfo%output,'(8X,a)') ' Atom            Energy '
    write(DECinfo%output,'(8X,a)') '------    --------------------'

    singleenergy=0E0_realk
    do i=1,natoms
       if(orbitals_assigned(i)) then
          write(DECinfo%output,'(1X,a,i6,4X,g20.10)') '#SING#', i, energy_matrix(i,i)
          singleenergy = singleenergy + energy_matrix(i,i)
       end if
    end do

    write(DECinfo%output,*)
    write(DECinfo%output,*)

    write(DECinfo%output,'(8X,a)') '-- Pair fragments'
    write(DECinfo%output,'(8X,a)') '------    ------    ----------    --------------------'
    write(DECinfo%output,'(8X,a)') '   P         Q        R(Ang)           deltaE(PQ) '
    write(DECinfo%output,'(8X,a)') '------    ------    ----------    --------------------'
    pairenergy=0E0_realk
    do i=1,natoms
       do j=i+1,natoms

          ! print increments only if pair interaction energy is nonzero
          if( orbitals_assigned(i) .and. orbitals_assigned(j) ) then
             write(DECinfo%output,'(1X,a,i6,4X,i6,4X,g10.4,4X,g20.10)') '#PAIR#',i,j,&
                  &bohr_to_angstrom*DistanceTable(i,j), energy_matrix(i,j)
             pairenergy = pairenergy + energy_matrix(i,j)
          end if

       end do
    end do
    Ecorr = singleenergy + pairenergy

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a)')   '------------------------------------------------------'
    write(DECinfo%output,'(1X,a,g20.10)') 'Atomic fragment energy   = ', singleenergy
    write(DECinfo%output,'(1X,a,g20.10)') 'Pair atomic fragment en. = ', pairenergy
    write(DECinfo%output,'(1X,a,g20.10)') 'Total correlation energy = ', Ecorr
    write(DECinfo%output,'(1X,a)')    '-----------------------------------------------------'

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a)') '-- Atomic fragments - the four contributions'
    write(DECinfo%output,'(1X,a)') '********************************************'
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(7X,a)') '-- Contribution 1'
    write(DECinfo%output,'(7X,a)') '------    --------------------'
    write(DECinfo%output,'(7X,a)') ' Atom            Energy '
    write(DECinfo%output,'(7X,a)') '------    --------------------'

    do i=1,natoms
       if(orbitals_assigned(i)) then
          write(DECinfo%output,'(a,i6,4X,g20.10)') '#Se1#',i, e1(i,i)
       end if
    end do


    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(7X,a)') '-- Contribution 2'
    write(DECinfo%output,'(7X,a)') '------    --------------------'
    write(DECinfo%output,'(7X,a)') ' Atom            Energy '
    write(DECinfo%output,'(7X,a)') '------    --------------------'

    do i=1,natoms
       if(orbitals_assigned(i)) then
          write(DECinfo%output,'(a,i6,4X,g20.10)') '#Se2#',i, e2(i,i)
       end if
    end do

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(7X,a)') '-- Contribution 3'
    write(DECinfo%output,'(7X,a)') '------    --------------------'
    write(DECinfo%output,'(7X,a)') ' Atom            Energy '
    write(DECinfo%output,'(7X,a)') '------    --------------------'

    do i=1,natoms
       if(orbitals_assigned(i)) then
          write(DECinfo%output,'(a,i6,4X,g20.10)') '#Se3#',i, e3(i,i)
       end if
    end do


    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(7X,a)') '-- Contribution 4'
    write(DECinfo%output,'(7X,a)') '------    --------------------'
    write(DECinfo%output,'(7X,a)') ' Atom            Energy '
    write(DECinfo%output,'(7X,a)') '------    --------------------'

    do i=1,natoms
       if(orbitals_assigned(i)) then
          write(DECinfo%output,'(a,i6,4X,g20.10)') '#Se4#',i, e4(i,i)
       end if
    end do


    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a)') '-- Pair fragments - the four contributions'
    write(DECinfo%output,'(1X,a)') '******************************************'
    write(DECinfo%output,*)
    write(DECinfo%output,*)

    write(DECinfo%output,'(7X,a)') '-- Contribution 1'
    write(DECinfo%output,'(7X,a)') '------    ------    ----------    --------------------'
    write(DECinfo%output,'(7X,a)') '   P         Q        R(Ang)           deltaE(PQ) '
    write(DECinfo%output,'(7X,a)') '------    ------    ----------    --------------------'
    pairenergy=0E0_realk
    do i=1,natoms
       do j=i+1,natoms

          ! print increments only if pair interaction energy is nonzero
          if( orbitals_assigned(i) .and. orbitals_assigned(j) ) then
             write(DECinfo%output,'(1X,a,i6,4X,i6,4X,g10.4,4X,g20.10)') '#Pe1#', i,j,&
                  &bohr_to_angstrom*DistanceTable(i,j), e1(i,j)
          end if

       end do
    end do
    write(DECinfo%output,*)
    write(DECinfo%output,*)

    write(DECinfo%output,'(7X,a)') '-- Contribution 2'
    write(DECinfo%output,'(7X,a)') '------    ------    ----------    --------------------'
    write(DECinfo%output,'(7X,a)') '   P         Q        R(Ang)           deltaE(PQ) '
    write(DECinfo%output,'(7X,a)') '------    ------    ----------    --------------------'
    pairenergy=0E0_realk
    do i=1,natoms
       do j=i+1,natoms

          ! print increments only if pair interaction energy is nonzero
          if( orbitals_assigned(i) .and. orbitals_assigned(j) ) then
             write(DECinfo%output,'(1X,a,i6,4X,i6,4X,g10.4,4X,g20.10)') '#Pe2#', i,j,&
                  &bohr_to_angstrom*DistanceTable(i,j), e2(i,j)
          end if

       end do
    end do
    write(DECinfo%output,*)
    write(DECinfo%output,*)


    write(DECinfo%output,'(7X,a)') '-- Contribution 3'
    write(DECinfo%output,'(7X,a)') '------    ------    ----------    --------------------'
    write(DECinfo%output,'(7X,a)') '   P         Q        R(Ang)           deltaE(PQ) '
    write(DECinfo%output,'(7X,a)') '------    ------    ----------    --------------------'
    pairenergy=0E0_realk
    do i=1,natoms
       do j=i+1,natoms

          ! print increments only if pair interaction energy is nonzero
          if( orbitals_assigned(i) .and. orbitals_assigned(j) ) then
             write(DECinfo%output,'(1X,a,i6,4X,i6,4X,g10.4,4X,g20.10)') '#Pe3#', i,j,&
                  &bohr_to_angstrom*DistanceTable(i,j), e3(i,j)
          end if

       end do
    end do
    write(DECinfo%output,*)
    write(DECinfo%output,*)


    write(DECinfo%output,'(7X,a)') '-- Contribution 4'
    write(DECinfo%output,'(7X,a)') '------    ------    ----------    --------------------'
    write(DECinfo%output,'(7X,a)') '   P         Q        R(Ang)           deltaE(PQ) '
    write(DECinfo%output,'(7X,a)') '------    ------    ----------    --------------------'
    pairenergy=0E0_realk
    do i=1,natoms
       do j=i+1,natoms

          ! print increments only if pair interaction energy is nonzero
          if( orbitals_assigned(i) .and. orbitals_assigned(j) ) then
             write(DECinfo%output,'(1X,a,i6,4X,i6,4X,g10.4,4X,g20.10)') '#Pe4#', i,j,&
                  &bohr_to_angstrom*DistanceTable(i,j), e4(i,j)
          end if

       end do
    end do
    write(DECinfo%output,*)
    write(DECinfo%output,*)

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    call mem_dealloc(ppfock)

  end subroutine Full_DEC_calculation_Lagrangian



  !> \brief Energy control and printout center.
  !> 1. Estimate error of correlation energy Eerr as difference between energies
  !>    for the three partitioning schemes (see PCCP, 14, 15706 (2012))
  !> 2. Estimate contibutions from neglected pairs.
  !> 3. Accept current pair cut off if estimated contributions from neglected pairs
  !>    are less than half of Err for all three partitioning schemes.
  !> \author Kasper Kristensen
  !> \date October 2012
  subroutine dec_energy_control_center(natoms,paircut,dofrag,DistanceTable,FragEnergiesAll,&
       & Fragments, morepairsneeded, newpaircut,Eerr)

    implicit none
    !> Number of atoms in molecule
    integer,intent(in) :: natoms
    !> Pair cut off distance used
    !> (i.e. FragEnergies must contain pair energies with pair distances from 0 to paircut, while
    !> pairs separated by more than paircut have not been treated and need to be extrapolated)
    real(realk),intent(in) :: paircut
    !> dofrag(P) is true if P has orbitals assigned
    logical,intent(in) :: dofrag(natoms)
    !> Distances between atoms
    real(realk),dimension(natoms,natoms),intent(in) :: DistanceTable
    !> Fragment energies as determined in main_fragment_driver.
    real(realk),dimension(natoms,natoms,ndecenergies),intent(in) :: FragEnergiesAll
    !> Single  fragments
    type(ccatom),intent(inout) :: Fragments(natoms)
    !> Do we need to include more pairs...? (Answer to point 3 above)
    logical, intent(inout) :: morepairsneeded
    !> New pair cutoff distance (identical to paircut if morepairsneeded=.false.)
    real(realk),intent(inout) :: newpaircut
    !> Estimated intrinsic energy error of DEC calculation
    real(realk),intent(inout) :: Eerr
    real(realk) :: Lmiss, EOCCmiss, EVIRTmiss,maxmiss
    ! Number of terms in fitting function, currently just nterms=1: f(x) = a(1)*x^{-6}
    ! (of nterms=2, then f(x) = a(1)*x^{-6} + a(2)*x^{-7} etc.)
    integer,parameter :: nterms=1
    ! Fitting parameters
    real(realk),dimension(nterms) :: alag,aocc,avirt
    real(realk) :: L,Eocc,Evirt
    integer :: i,j,npairsdone,npairstot, nnewpairsL, nnewpairsO, nnewpairsV, nnewpairs
    logical :: allpairs
    real(realk) :: Lnewpaircut, Onewpaircut, Vnewpaircut, Etarget,dist
    real(realk),pointer :: FragEnergiesModel(:,:,:)

    ! Extract fragment energies for model under consideration
    call mem_alloc(FragEnergiesModel,natoms,natoms,3)
    call extract_fragenergies_for_model(natoms,FragEnergiesAll,FragEnergiesModel)

    ! Total number of pairs: N*(N-1)/2
    ! where N is the number of atoms with orbitals assigned
    npairstot = count(dofrag)*(count(dofrag)-1)/2


    ! Calculate current energies using the three part. schemes
    ! ********************************************************
    L=0.0E0_realk
    Eocc=0.0E0_realk
    Evirt=0.0E0_realk
    npairsdone=0
    do i=1,natoms
       if(.not. dofrag(i)) cycle  ! no orbitals assigned to atom "i"
       L = L + FragEnergiesModel(i,i,1)
       Eocc = Eocc + FragEnergiesModel(i,i,2)
       Evirt = Evirt + FragEnergiesModel(i,i,3)
       do j=i+1,natoms
          if(.not. dofrag(j)) cycle  ! no orbitals assigned to atom "j"
          if(DistanceTable(i,j) > paircut) cycle   ! only pairs below cutoff have been determined
          npairsdone = npairsdone +1
          L = L + FragEnergiesModel(i,j,1)
          Eocc = Eocc + FragEnergiesModel(i,j,2)
          Evirt = Evirt + FragEnergiesModel(i,j,3)
       end do
    end do



    ! SANITY CHECKS
    ! *************

    ! Sanity check 1
    if(DECinfo%simulate_full) then ! no pairs in calculation, and therefore no pairs are needed
       write(DECinfo%output,*) 'Simulated full calculation - skipping pair analysis'
       morepairsneeded=.false.
       call mem_dealloc(FragEnergiesModel)
       return
    end if

    ! Sanity check 2
    if(npairsdone==npairstot) then ! all pairs have already been calculated
       allpairs=.true.
       morepairsneeded=.false.  ! no more pairs needed
    else
       allpairs=.false.
    end if

    ! Sanity check 3
    if(npairsdone==0 .and. (count(dofrag)>1) ) then 
       ! No pairs have been calculated, but there is more than one fragment,
       ! and thus at least one pair.
       ! The remaining part of the routine will be meaningless.
       ! This is presumably a debug calculation.
       write(DECinfo%output,*) 'WARNING: Pair energy estimate could not be carried out, there are no pairs!'
       morepairsneeded=.false.
       call mem_dealloc(FragEnergiesModel)
       return
    end if

    ! Sanity check 4
    if(npairstot==1) then 
       ! Regression plot is meaningless if there is only one pair
       write(DECinfo%output,*) 'Pair energy estimate could not be carried out, there is only one pair!'
       if(npairsdone==1) then  ! The single pair has already been calculated
          morepairsneeded=.false.
       else ! require the single pair to be calculated
          morepairsneeded=.true.
       end if
       call mem_dealloc(FragEnergiesModel)
       return
    end if

    ! Sanity check 5
    ! We only perform regression for pair distances beyond DECinfo%PairMinDist
    ! in estimate_remaining_pair_energies, and therefore paircutoff beyond this distance is meaningless
    if(paircut < DECinfo%PairMinDist) then
       write(DECinfo%output,*) 'Pair cut off smaller than 3 Angstrom, skipping pair regression.'
       morepairsneeded=.true.
       call mem_dealloc(FragEnergiesModel)
       return
    end if


    ! Estimate error as largest different between these energies
    ! **********************************************************
    Eerr = max(L,Eocc,Evirt) - min(L,Eocc,Evirt)


    ! Estimate contributions from neglected pairs
    ! *******************************************
    ! Lagrangian
    call estimate_remaining_pair_energies(natoms,paircut,DistanceTable,FragEnergiesModel(:,:,1),dofrag,&
         & Lmiss, nterms, alag)
    ! Occupied
    call estimate_remaining_pair_energies(natoms,paircut,DistanceTable,FragEnergiesModel(:,:,2),dofrag,&
         & EOCCmiss, nterms, aocc)
    ! Virtual
    call estimate_remaining_pair_energies(natoms,paircut,DistanceTable,FragEnergiesModel(:,:,3),dofrag,&
         & EVIRTmiss, nterms, avirt)

    ! Maximum missing contribution (all these should be negative)
    maxmiss = max( abs(Lmiss), abs(EOCCmiss), abs(EVIRTmiss) )



    ! Ensure that missing contribution is well below estimated error
    ! **************************************************************
    ! Current criteria: Calculate more pairs if maxmiss is larger than half of Eerr
    Etarget = 0.5E0_realk*Eerr

    if(.not. allpairs) then  ! only relevant if some pairs have not been considered
       if( maxmiss > Etarget) then
          morepairsneeded=.true.
       else
          morepairsneeded=.false.
       end if
    end if


    ! If more pairs are needed, determine new pair cutoff
    ! ***************************************************

    NeedMorePairs: if(morepairsneeded) then
       ! New pair cutoff estimated by Lagrangian scheme
       call determine_new_paircutoff(natoms,paircut,DistanceTable,FragEnergiesModel(:,:,1),dofrag,Etarget,&
            & Lmiss,nterms, alag,Lnewpaircut,nnewpairsL)

       ! New pair cutoff estimated by occupied scheme
       call determine_new_paircutoff(natoms,paircut,DistanceTable,FragEnergiesModel(:,:,2),dofrag,Etarget,&
            & EOCCmiss,nterms, aocc,Onewpaircut,nnewpairsO)

       ! New pair cutoff estimated by virtual scheme
       call determine_new_paircutoff(natoms,paircut,DistanceTable,FragEnergiesModel(:,:,3),dofrag,Etarget,&
            & EVIRTmiss,nterms, avirt,Vnewpaircut,nnewpairsV)

       ! Be on the safe side and use maximum value of new pair cut offs from the 3 schemes
       newpaircut = max(Lnewpaircut,Onewpaircut,Vnewpaircut)
       ! Number of new pairs
       nnewpairs = max(nnewpairsL,nnewpairsO,nnewpairsV)

    else
       ! Everything is fine with current pair cutoff - set new pair cutoff equal to existing one
       newpaircut = paircut
       nnewpairs=0
    end if NeedMorePairs

    ! Print out energy statistics
    call dec_energy_pairanalysis_printout(natoms,paircut,DistanceTable,FragEnergiesModel,&
         & morepairsneeded,newpaircut, Lmiss,EOCCmiss,EVIRTmiss,Eerr,L,Eocc,Evirt,alag,aocc,avirt,&
         & allpairs,dofrag,nnewpairs)

    if(DECinfo%NoExtraPairs .and. morepairsneeded) then
       write(DECinfo%output,*) 'WARNING! You have turned off pair control!'
       write(DECinfo%output,*) 'However, the pair analysis indicates that more pairs are needed!'
       write(DECinfo%output,*) 'Your final results may be inaccurate!'
       write(DECinfo%output,*) 'I recommend that you restart the calculation WITHOUT'
       write(DECinfo%output,*) 'the .NoExtraPairs keyword and WITH the .restart option.'
       Eerr=0.0_realk
       newpaircut = paircut
       morepairsneeded = .false.
    end if

    call mem_dealloc(FragEnergiesModel)

  end subroutine dec_energy_control_center




  !> \brief Print out information from DEC energy control center (see dec_energy_control_center).
  !> \author Kasper Kristensen
  !> \date October 2012
  subroutine dec_energy_pairanalysis_printout(natoms,paircut,DistanceTable,FragEnergies,&
       & morepairsneeded,newpaircut,Lmiss,EOCCmiss,EVIRTmiss,Eerr,L,Eocc,Evirt,&
       & alag,aocc,avirt,allpairs,dofrag,nnewpairs)

    implicit none
    !> Number of atoms in molecule
    integer,intent(in) :: natoms
    !> Pair cut off distance used
    !> (i.e. FragEnergies must contain pair energies with pair distances from 0 to paircut, while
    !> pairs separated by more than paircut have not been treated and need to be extrapolated)
    real(realk),intent(in) :: paircut
    !> Distances between atoms
    real(realk),dimension(natoms,natoms),intent(in) :: DistanceTable
    !> Fragment energies for Lagrangian (:,:,1), occupied (:,:,2), and virtual (:,:,3) schemes
    real(realk),dimension(natoms,natoms,3),intent(in) :: FragEnergies
    !> Do we need to include more pairs...? (see point 3 in dec_energy_control_center)
    logical, intent(in) :: morepairsneeded
    !> New pair cut off (only printed if morepairsneeded is true)
    real(realk),intent(in) :: newpaircut
    !> Energy from missing pairs for Lagrangian partitioning scheme
    real(realk),intent(in) :: Lmiss
    !> Energy from missing pairs for occupied partitioning scheme
    real(realk),intent(in) :: EOCCmiss
    !> Energy from missing pairs for virtual partitioning scheme
    real(realk),intent(in) :: EVIRTmiss
    !> Estimated intrinsic energy error of DEC calculation
    real(realk),intent(in) :: Eerr
    !> Correlation energy for Lagrangian scheme (including pairs up to paircut)
    real(realk),intent(in) :: L
    !> Correlation energy for occupied scheme (including pairs up to paircut)
    real(realk),intent(in) :: Eocc
    !> Correlation energy for virtual scheme (including pairs up to paircut)
    real(realk),intent(in) :: Evirt
    !> Fitting parameters for Lagrangian scheme
    real(realk),intent(in) :: alag(1)
    !> Fitting parameters for occupied scheme
    real(realk),intent(in) :: aocc(1)
    !> Fitting parameters for virtual scheme
    real(realk),intent(in) :: avirt(1)
    !> Have all pairs been included or not?
    logical,intent(in) :: allpairs
    !> dofrag(P) is true if P has orbitals assigned
    logical,intent(in) :: dofrag(natoms)
    !> Number of new pairs to calculate
    integer,intent(in) :: nnewpairs

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a)') '=================================================================='
    write(DECinfo%output,'(1X,a)') '|                 DEC PAIR ENERGY CONTROL CENTER                 |'
    write(DECinfo%output,'(1X,a)') '=================================================================='
    write(DECinfo%output,*)
    if(allpairs) then
       write(DECinfo%output,'(1X,a)') 'All pairs have been included!'
    else
       write(DECinfo%output,'(1X,a,f10.2,a,f10.2,a)') 'Pair cutoff distance = ', paircut, ' a.u. = ', &
            & bohr_to_angstrom*paircut, ' Angstrom'
    end if
    write(DECinfo%output,*)
    if(.not. allpairs) then
       write(DECinfo%output,'(1X,a)') '---------------------------------------------------------------'
       write(DECinfo%output,'(1X,a)') ' Calculated energies with pair contributions up to pair cutoff '
       write(DECinfo%output,'(1X,a)') '---------------------------------------------------------------'
    end if
    write(DECinfo%output,'(1X,a,g20.9)') 'Lagrangian scheme energy      : ', L
    write(DECinfo%output,'(1X,a,g20.9)') 'Occupied scheme energy        : ', Eocc
    write(DECinfo%output,'(1X,a,g20.9)') 'Virtual scheme energy         : ', Evirt
    write(DECinfo%output,'(1X,a,g20.9)') '*** Estimated intrinsic error : ', Eerr

    if(.not. allpairs) then
       write(DECinfo%output,*)
       write(DECinfo%output,'(1X,a)') '----------------------------------------------------------------'
       write(DECinfo%output,'(1X,a)') '  Regression-estimated energy contributions from missing pairs  '
       write(DECinfo%output,'(1X,a)') '----------------------------------------------------------------'
       write(DECinfo%output,'(1X,a,g20.9)') 'Lagrangian scheme missing  : ', Lmiss
       write(DECinfo%output,'(1X,a,g20.9)') 'Occupied scheme missing    : ', EOCCmiss
       write(DECinfo%output,'(1X,a,g20.9)') 'Virtual scheme missing     : ', EVIRTmiss
       write(DECinfo%output,'(1X,a,g20.9)') '*** Max (absolute) missing : ', &
            & max( abs(Lmiss), abs(EOCCmiss), abs(EVIRTmiss) )

       write(DECinfo%output,*)
       write(DECinfo%output,'(1X,a)') '----------------------'
       write(DECinfo%output,'(1X,a)') 'PAIR CUTOFF CONCLUSION'
       write(DECinfo%output,'(1X,a)') '----------------------'
       if(morepairsneeded) then
          write(DECinfo%output,'(1X,a)') 'The maxium estimated missing pair energy contribution'
          write(DECinfo%output,'(1X,a)') 'is larger than half of the estimated intrinsic energy error.' 
          write(DECinfo%output,'(1X,a)') 'PairConclusion: Include more pairs to adapt to precision!'
          write(DECinfo%output,'(1X,a,f10.2,a,f10.2,a)') 'NEW pair cutoff distance = ', newpaircut, &
               & ' a.u. = ', bohr_to_angstrom*newpaircut, ' Angstrom'
          write(DECinfo%output,'(1X,a,i8)') 'Number of new pairs to calculate : ', nnewpairs
       else
          write(DECinfo%output,'(1X,a)') 'The maxium estimated missing pair energy contribution'
          write(DECinfo%output,'(1X,a)') 'is smaller than half of the estimated intrinsic energy error.' 
          write(DECinfo%output,'(1X,a)') 'PairConclusion: No need to include more pairs!'
       end if


       write(DECinfo%output,*)
       write(DECinfo%output,'(1X,a)') '-------------------------------------------------------------'
       write(DECinfo%output,'(1X,a)') '             Estimated total energies (all pairs)            '
       write(DECinfo%output,'(1X,a)') '-------------------------------------------------------------'
       write(DECinfo%output,'(1X,a,g20.9)') 'Lagrangian scheme estimat. : ', L+Lmiss
       write(DECinfo%output,'(1X,a,g20.9)') 'Occupied scheme estimat.   : ', Eocc+EOCCmiss
       write(DECinfo%output,'(1X,a,g20.9)') 'Virtual scheme estimat.    : ', Evirt+EVIRTmiss

    end if

    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a)') '-------------------------------------------------------------'
    write(DECinfo%output,'(1X,a)') '           Pair energy regression parameters                 '
    write(DECinfo%output,'(1X,a)') '-------------------------------------------------------------'
    write(DECinfo%output,'(1X,a)') 'The calculated pair energies were fitted to the function:'
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a)') 'f(x) = a*x^{-6} '
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a)') 'where x is the pair distance.'
    write(DECinfo%output,'(1X,a)') 'Only pairs in the dispersion range'
    write(DECinfo%output,'(1X,a)') '(3 Angstrom to pair cutoff distance) were included in the fit'
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a)') 'The x^{-6} term is the leading order term as shown in:'
    write(DECinfo%output,'(1X,a)') 'Hoyvik et al., J. Chem. Phys. 136, 014105 (2012), app. C     '
    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a,g18.8,a,g18.8,a)') 'Lag. scheme : a = ',alag(1), &
         & ' Hartree*bohr^6 = ', alag(1)*(bohr_to_angstrom**6), ' Hartree*Angstrom^6'
    write(DECinfo%output,'(1X,a,g18.8,a,g18.8,a)') 'Occ. scheme : a = ',aocc(1), &
         & ' Hartree*bohr^6 = ', aocc(1)*(bohr_to_angstrom**6), ' Hartree*Angstrom^6'
    write(DECinfo%output,'(1X,a,g18.8,a,g18.8,a)') 'Vir. scheme : a = ',avirt(1), &
         & ' Hartree*bohr^6 = ', avirt(1)*(bohr_to_angstrom**6), ' Hartree*Angstrom^6'
    write(DECinfo%output,*)


    write(DECinfo%output,*)
    write(DECinfo%output,'(1X,a)') '-------------------------------------------------------------'
    write(DECinfo%output,'(1X,a)') '                 PAIR INTERACTION ENERGY PLOT                '
    write(DECinfo%output,'(1X,a)') '-------------------------------------------------------------'
    if(DECinfo%ccmodel==1) then ! MP2: Lagrangian
       write(DECinfo%output,'(1X,a)') 'Plot contains maximum pair energies for Lagrangian scheme'
       write(DECinfo%output,'(1X,a)') 'in intervals of 1 Angstrom.'
    else  ! higher order CC: currently averaged occupied-virtual scheme
       write(DECinfo%output,'(1X,a)') 'Plot contains maximum pair energies for averaged occ-virt scheme'
       write(DECinfo%output,'(1X,a)') 'in intervals of 1 Angstrom.'
    end if
    write(DECinfo%output,*)
    call plot_pair_energies(natoms,paircut,FragEnergies(:,:,1),DistanceTable,dofrag)


    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*)

  end subroutine dec_energy_pairanalysis_printout



  !> \brief Print a simple "ascii-art" plot of the largest pair energies.
  !> \author Kasper Kristensen
  !> \date October 2012
  subroutine plot_pair_energies(natoms,paircut,FragEnergies,DistanceTable,dofrag)

    implicit none
    !> Number of atoms in molecule
    integer,intent(in) :: natoms
    !> Pair cut off distance used (in a.u.)
    real(realk),intent(in) :: paircut
    !> Fragment energies 
    real(realk),dimension(natoms,natoms),intent(in) :: FragEnergies
    !> Distances between atoms
    real(realk),dimension(natoms,natoms),intent(in) :: DistanceTable
    !> dofrag(P) is true if P has orbitals assigned
    logical,intent(in) :: dofrag(natoms)
    real(realk),pointer :: DistAng(:,:), xpoints(:), ypoints(:), xpoints2(:), ypoints2(:)
    integer :: mindist, maxdist, tmp, distInt,idx,npoints2
    real(realk) :: cutAng,endpoint,ln10
    integer :: i,j,npoints,k,interval
    character(len=16) :: xlabel, ylabel
    logical,pointer :: anypoints(:)

    ! Get distance table and pair cut off in Angstrom
    call mem_alloc(DistAng,natoms,natoms)
    DistAng = bohr_to_angstrom*DistanceTable
    cutAng = bohr_to_angstrom*paircut

    ! Minimum and maximum pair distances
    mindist=1000
    maxdist=0
    do i=1,natoms
       if(.not. dofrag(i)) cycle
       do j=i+1,natoms

          if(DistAng(i,j) > cutAng) cycle
          if(.not. dofrag(j)) cycle

          ! Nearest integer smaller than actual pair distance +0.5Angstrom
          ! (this is most appropriate when we consider intervals of 1 Angstrom below,
          !  because we then ensure that we don't have any intervals with very few points)
          tmp = floor(DistAng(i,j) + 0.4999E0_realk)

          ! Minimum
          ! Safety precaution: Never consider pairs separated by less than 1 Ang
          if(mindist > tmp .and. tmp>0) mindist = tmp

          ! Max
          if(maxdist < tmp) maxdist = tmp

       end do
    end do


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



  !> \brief Assuming that pair energies for pair distances from 0 Å to some pair cutoff distance
  !> have been calculated, use regression to estimate contributions from remaining pairs.
  !> This is done by fitting exiting pair energies to function: 
  !> f(x) = a(1)*x^{-6} + a(2)*x^{-7} + a(3)*x^{-8} + ...  (x is pair distance)
  !> The theoretical foundation for this fit is given in Appendix C of JCP 136,014105 (2012).
  !> \author Kasper Kristensen
  !> \date October 2012
  subroutine estimate_remaining_pair_energies(natoms,paircut,DistanceTable,FragEnergies,dofrag,Emiss,&
       & nterms, a)

    implicit none
    !> Number of atoms in molecule
    integer,intent(in) :: natoms
    !> Pair cut off distance used
    !> (i.e. FragEnergies must contain pair energies with pair distances from 0 to paircut, while
    !> pairs separated by more than paircut have not been treated and need to be extrapolated)
    real(realk),intent(in) :: paircut
    !> Distances between atoms
    real(realk),dimension(natoms,natoms),intent(in) :: DistanceTable
    !> Fragment energies 
    real(realk),dimension(natoms,natoms),intent(in) :: FragEnergies
    !> dofrag(P) is true if P has orbitals assigned
    logical,intent(in) :: dofrag(natoms)
    !> Estimated energy contribution from missing pairs (separated by more than paircut)
    real(realk),intent(inout) :: Emiss
    !> Number of terms in fitting function, e.g. for nterms=2: f(x) = a(1)*x^{-6} + a(2)*x^{-7}
    integer,intent(in) :: nterms
    !> Fitting parameters
    real(realk),intent(inout) :: a(nterms)
    real(realk),dimension(nterms) :: powers
    integer :: i,j,k,npairsdone,idx
    real(realk),pointer :: R(:),E(:)
    real(realk) :: mindist,errornorm


    ! Only consider pairs separated by more than 3 Angstrom
    mindist = DECinfo%PairMinDist


    ! Count how many pairs have already been calculated (number of pair distances below paircut)
    npairsdone=0
    do i=1,natoms
       if(.not. dofrag(i)) cycle  ! only consider fragment with orbital assigned
       do j=i+1,natoms
          if(.not. dofrag(j)) cycle
          if(DistanceTable(i,j) < mindist) cycle
          if(DistanceTable(i,j) < paircut) npairsdone = npairsdone +1
       end do
    end do

    ! Sanity check - if no pairs are done, then return
    if(npairsdone == 0) then
       Emiss = 0.0E0_realk
       a = 0.0E0_realk
       return
    end if

    ! Collect distances (R) and energies (E) for calculated pairs
    call mem_alloc(R,npairsdone)
    call mem_alloc(E,npairsdone)
    idx=0
    do i=1,natoms
       if(.not. dofrag(i)) cycle  ! only consider fragment with orbital assigned
       do j=i+1,natoms
          if(.not. dofrag(j)) cycle
          if(DistanceTable(i,j) < mindist) cycle
          if(DistanceTable(i,j) < paircut) then
             idx=idx+1
             R(idx) = DistanceTable(i,j)
             E(idx) = FragEnergies(i,j)
          end if
       end do
    end do


    ! *************************************************************
    ! For calculated (R,E) points, make regression fit to function:
    ! f(x) = a(1)*x^{-6} + a(2)*x^{-7} + ...    (for nterms)
    ! *************************************************************

    ! Set powers: -6,-7,-8,...
    do i=1,nterms
       powers(i) = -real(6+i-1)
    end do

    ! Perform regression
    call dec_regression(npairsdone,R,E,nterms,powers,a,errornorm)
    if(DECinfo%PL>1) write(DECinfo%output,'(1X,a,g20.8)') 'Pair regression, error norm:', errornorm



    ! Estimate contributions from missing pairs based on fit
    ! ******************************************************
    Emiss = 0.0E0_realk
    do i=1,natoms
       if(.not. dofrag(i)) cycle  ! only consider fragment with orbital assigned
       do j=i+1,natoms
          if(.not. dofrag(j)) cycle

          if(DistanceTable(i,j) < mindist) cycle

          AddPairEstimate: if(DistanceTable(i,j) > paircut) then 

             ! Estimate pair (i,j) contribution using fitted function 
             ! f(x) = a(1)*x^{-6} + a(2)*x^{-7} + ...    (for nterms)
             ! and add contribution to missing energy.
             do k=1,nterms
                Emiss = Emiss + a(k)*DistanceTable(i,j)**(powers(k))
             end do

          end if AddPairEstimate

       end do
    end do



    call mem_dealloc(R)
    call mem_dealloc(E)

  end subroutine estimate_remaining_pair_energies




  !> \brief For a given set of fragment energies and a given pair cutoff distance, 
  !> use data from pair regression fit (see estimate_remaining_pair_energies) to
  !> determine new pair cutoff distance beyond which the error associated with the 
  !> neglected pairs is below some target energy error.
  !> (To be on the safe side we add 1 a.u. to the actual estimated distance).
  !> \author Kasper Kristensen
  !> \date October 2012
  subroutine determine_new_paircutoff(natoms,paircut,DistanceTable,FragEnergies,dofrag,Etarget,Emiss,&
       & nterms, a,newpaircut,nnewpairs)

    implicit none
    !> Number of atoms in molecule
    integer,intent(in) :: natoms
    !> Pair cut off distance currently used
    real(realk),intent(in) :: paircut
    !> Distances between atoms
    real(realk),dimension(natoms,natoms),intent(in) :: DistanceTable
    !> Fragment energies (atomic units)
    real(realk),dimension(natoms,natoms),intent(in) :: FragEnergies
    !> dofrag(P) is true if P has orbitals assigned
    logical,intent(in) :: dofrag(natoms)
    !> Target energy error
    real(realk),intent(in) :: Etarget
    !> Estimated energy contribution from ALL missing pairs
    real(realk),intent(in) :: Emiss
    !> Number of terms in fitting function, e.g. for nterms=2: f(x) = a(1)*x^{-6} + a(2)*x^{-7}
    integer,intent(in) :: nterms
    !> Fitting parameters (in atomic units)
    real(realk),intent(in) :: a(nterms)
    !> New pair cut off distance
    real(realk),intent(inout) :: newpaircut
    !> Number of additional pairs to include
    integer,intent(inout) :: nnewpairs
    real(realk),dimension(nterms) :: powers
    integer :: i,j,k,maxexpansion,l
    real(realk) :: Enewpairs, pairerror
    logical :: conv


    ! Init stuff
    ! **********

    ! Sanity check
    if(Etarget <=0.0E0_realk) then
       call lsquit('determine_new_paircutoff: target energy error must be positive!',-1)
    end if
    
    ! quit if pair distance increment of 40 bohr is not enough (in which case something is very wrong)
    maxexpansion = 40 
    newpaircut = paircut

    ! Fitting function: f(x) = a(1)*x^{-6} + a(2)*x^{-7} + ...    (see estimate_remaining_pair_energies)
    ! Set powers: -6,-7,-8,...
    do i=1,nterms
       powers(i) = -real(6+i-1)
    end do



    ! Increase pair cut off until estimated pair error is below target energy error
    ! *****************************************************************************
    conv=.false.
    nnewpairs=0
    IncreasePairCutOff: do k=1,maxexpansion

       ! Zero accumulated energy for new pairs to include
       Enewpairs = 0.0E0_realk

       ! Increase new pair cutoff by 1 bohr
       newpaircut = newpaircut + real(k)

       do i=1,natoms
          if(.not. dofrag(i)) cycle  ! only consider fragment with orbital assigned
          do j=i+1,natoms
             if(.not. dofrag(j)) cycle

             ! Add additional pair energy contributions obtained with new pair cutoff,
             ! i.e. pair energy contributions for distances between current and new pair cutoffs.
             AddPairEstimate: if(DistanceTable(i,j) > paircut &
                  & .and. DistanceTable(i,j) < newpaircut ) then 

                ! Estimate pair (i,j) contribution using fitted function 
                ! f(x) = a(1)*x^{-6} + a(2)*x^{-7} + ...    (for nterms)
                ! and add contribution to missing energy.
                do l=1,nterms
                   Enewpairs = Enewpairs + a(l)*DistanceTable(i,j)**(powers(l))
                end do

                ! Update number of new pairs
                nnewpairs= nnewpairs+1

             end if AddPairEstimate

          end do
       end do

       ! Estimated pair cutoff error using new pair cutoff, i.e.:
       ! E(all missing pairs) - E(pairs between current and new cutoff)
       pairerror = abs(Emiss - Enewpairs)
       
       ! Exit if new estimate is smaller than requested target energy
       if(pairerror < Etarget) then
          conv=.true.
          exit IncreasePairCutOff
       end if

    end do IncreasePairCutOff

    ! Sanity check
    if(.not. conv .or. DECinfo%PL>1) then
       write(DECinfo%output,*) 'Missing pair energy: ', Emiss
       write(DECinfo%output,*) 'Pair error ', Emiss
       write(DECinfo%output,*) 'Old pair distance ', paircut
       write(DECinfo%output,*) 'New pair distance ', newpaircut
    end if
    if(.not. conv) then
       call lsquit('determine_new_paircutoff: Pair cut off expansion did not converge!',-1)
    end if


    ! Now newpaircutoff contains the estimated pair cutoff distance to ensure
    ! that the pair cutoff error is smaller than the target error.
    ! Since this is just a rough estimate we add 1 bohr to be on the safe side
    newpaircut = newpaircut + 1.0E0_realk    


  end subroutine determine_new_paircutoff




  ! ===================================================================!
  !                      FRAGMENT OPTIMIZATION                         !
  ! ===================================================================!
  

  !> Routine that optimizes an atomic fragment using the Lagrangian partitioning scheme
  !> and using a reduction scheme based on the individual orbitals.
  !> \date November 2011
  !> \author Ida-Marie Hoeyvik
  subroutine optimize_atomic_fragment(MyAtom,AtomicFragment,nAtoms, &
       &OccOrbitals,nOcc,UnoccOrbitals,nUnocc,DistanceTable, &
       &MyMolecule,mylsitem,freebasisinfo,t1full)
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
    type(ccatom),intent(inout)        :: AtomicFragment
    !> All occupied orbitals
    type(ccorbital), dimension(nOcc), intent(in)      :: OccOrbitals
    !> All unoccupied orbitals
    type(ccorbital), dimension(nUnocc), intent(in)    :: UnoccOrbitals
    !> Distance table for all atoms in molecule
    real(realk), dimension(nAtoms,nAtoms), intent(in) :: DistanceTable
    !> Full molecule information
    type(fullmolecule), intent(in) :: MyMolecule
    !> Integral information
    type(lsitem), intent(inout)       :: mylsitem
    !> Delete fragment basis information ("expensive box in ccatom type") at exit?
    logical,intent(in) :: freebasisinfo
    !> t1 amplitudes for full molecule to be updated (only used when DECinfo%SinglesPolari is set)
    type(array2),intent(inout),optional :: t1full
    integer :: savemodel

    ! Save existing model and do fragment optimization with MP2
    savemodel = DECinfo%ccmodel
    DECinfo%ccmodel = 1

    if(DECinfo%SinglesPolari) then 
       ! currently we store full amplitudes but not AOS amplitudes,
       ! to be modified!
       call optimize_atomic_fragment_main(MyAtom,AtomicFragment,nAtoms, &
            &OccOrbitals,nOcc,UnoccOrbitals,nUnocc,DistanceTable, &
            &MyMolecule,mylsitem,freebasisinfo,t1full=t1full)
    else   


       if(DECinfo%fragadapt) then
          call optimize_atomic_fragment_FA(MyAtom,AtomicFragment,nAtoms, &
               &OccOrbitals,nOcc,UnoccOrbitals,nUnocc,DistanceTable, &
               &MyMolecule,mylsitem,freebasisinfo)
       else
          call optimize_atomic_fragment_main(MyAtom,AtomicFragment,nAtoms, &
               &OccOrbitals,nOcc,UnoccOrbitals,nUnocc,DistanceTable, &
               &MyMolecule,mylsitem,freebasisinfo)
       end if


    end if

    DECinfo%ccmodel = savemodel

  end subroutine optimize_atomic_fragment


  !> Routine that optimizes an atomic fragment using the Lagrangian partitioning scheme
  !> and using a reduction scheme based on the individual orbitals.
  !> \date November 2011
  !> \author Ida-Marie Hoeyvik
  subroutine optimize_atomic_fragment_main(MyAtom,AtomicFragment,nAtoms, &
       &OccOrbitals,nOcc,UnoccOrbitals,nUnocc,DistanceTable, &
       &MyMolecule,mylsitem,freebasisinfo,t1full)
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
    type(ccatom),intent(inout)        :: AtomicFragment
    !> All occupied orbitals
    type(ccorbital), dimension(nOcc), intent(in)      :: OccOrbitals
    !> All unoccupied orbitals
    type(ccorbital), dimension(nUnocc), intent(in)    :: UnoccOrbitals
    !> Distance table for all atoms in molecule
    real(realk), dimension(nAtoms,nAtoms), intent(in) :: DistanceTable
    !> Full molecule information
    type(fullmolecule), intent(in) :: MyMolecule
    !> Integral information
    type(lsitem), intent(inout)       :: mylsitem
    !> Delete fragment basis information ("expensive box in ccatom type") at exit?
    logical,intent(in) :: freebasisinfo
    !> t1 amplitudes for full molecule to be updated (only used when DECinfo%SinglesPolari is set)
    type(array2),intent(inout),optional :: t1full
    real(realk)                    :: LagEnergyDiff, OccEnergyDiff,VirtEnergyDiff
    real(realk)                    :: LagEnergyOld, OccEnergyOld, VirtEnergyOld
    real(realk)                    :: FOT,REDFOT
    logical, dimension(natoms)     :: Occ_atoms,Virt_atoms,OccOld,VirtOld
    real(realk),dimension(natoms)  :: DistMyAtom,SortedDistMyAtom
    integer,dimension(natoms)      :: DistTrackMyAtom
    integer,dimension(natoms)      :: nocc_per_atom,nunocc_per_atom
    integer      :: iter,i,idx
    integer      :: nocc_new, nvirt_new, nocc_old, nvirt_old, max_iter_red
    real(realk)  :: RejectThresh
    logical      :: converged,REDconverged
    logical :: expansion_converged, lag_converged, occ_converged, virt_converged
    real(realk),pointer :: OccContribs(:),VirtContribs(:)
    logical,pointer :: OccAOS_old(:), VirtAOS_old(:), OccAOS_new(:), VirtAOS_new(:)
    integer :: REDocc, REDvirt
    real(realk) :: slavetime, flops_slaves
    integer,pointer :: REDoccIDX(:), REDvirtIDX(:)

    write(DECinfo%output,'(a)')    ' FOP'
    write(DECinfo%output,'(a)')    ' FOP ==============================================='
    write(DECinfo%output,'(a,i4)') ' FOP  Site fragment generator for fragment,',MyAtom
    write(DECinfo%output,'(a)')    ' FOP ==============================================='
    write(DECinfo%output,'(a)')    ' FOP'


    ! Sanity check for singles polarization
    if(DECinfo%SinglesPolari) then
       if(.not. present(t1full)) then
          call lsquit('optimize_atomic_fragment: Full singles polarization is requrested &
               & but t1 argument is not present!',DECinfo%output)
       end if
    end if


    ! ******************************************
    ! **   Initialization of stuff needed..   **
    ! ******************************************
    converged=.false.
    REDconverged=.false.
    max_iter_red=15   ! allow 15 reduction steps (should be more than enough)
    FOT = DECinfo%FOT
    REDFOT = 10E0_realk*FOT
    DistMyAtom= DistanceTable(:,MyAtom)
    call GetSortedList(SortedDistMyAtom,DistTrackMyAtom,DistanceTable,natoms,MyAtom)
    nocc_per_atom=get_number_of_orbitals_per_atom(OccOrbitals,nocc,natoms)
    nunocc_per_atom=get_number_of_orbitals_per_atom(UnoccOrbitals,nunocc,natoms)

    ! MPI fragment statistics
    flops_slaves = 0.0E0_realk
    slavetime = 0.0E0_realk

    if( (nocc_per_atom(MyAtom) == 0) .and. (nunocc_per_atom(MyAtom) == 0) ) then
       write(DECinfo%output,*) 'FOP Skipping optimization of fragment ', MyAtom
       AtomicFragment%energies=0E0_realk
       return
    end if

    ! ******************************************
    ! **  Starting computation of fragment    **
    ! ******************************************
    if(DECinfo%simulate_full .or. DECinfo%InclFullMolecule) then
       ! skip calculation of energy and just initiate fragment including ALL orbitals
       write(DECinfo%output,*) 'FOP Fragment includes all orbitals and fragment optimization is skipped'
       call fragment_init_simulate_full(MyAtom,nunocc, nocc, OccOrbitals,UnoccOrbitals,&
            & MyMolecule,mylsitem,AtomicFragment,.true.)
       call PrintInfo_Lagrangian(AtomicFragment,0.0E0_realk,0.0E0_realk,0.0E0_realk,0)
       if(freebasisinfo) then
          call atomic_fragment_free_basis_info(AtomicFragment)
       end if
       return
    else
       call InitialFragment_L(Occ_atoms,Virt_atoms,DistMyAtom,nAtoms)
       call get_fragment_and_Energy(MyAtom,natoms,Occ_Atoms,Virt_Atoms,&
            & MyMolecule,MyLsitem,nocc,nunocc,OccOrbitals,UnoccOrbitals,AtomicFragment)
       ! MPI fragment statistics
       slavetime = slavetime +  AtomicFragment%slavetime
       flops_slaves = flops_slaves + AtomicFragment%flops_slaves
    end if


    ! Print initial fragment information
    write(DECinfo%output,*)'FOP'
    write(DECinfo%output,'(1X,a)') 'FOP========================================================='
    write(DECinfo%output,'(1X,a)') 'FOP               Initial fragment information'
    write(DECinfo%output,'(1X,a)') 'FOP---------------------------------------------------------'
    write(DECinfo%output,'(1X,a,i4)')    'FOP Init: Fragment number                  :', MyAtom
    write(DECinfo%output,'(1X,a,i4)')    'FOP Init: Number of orbitals in virt total :', &
         & AtomicFragment%nunoccAOS
    write(DECinfo%output,'(1X,a,i4)')    'FOP Init: Number of orbitals in occ total  :', &
         & AtomicFragment%noccAOS
    write(DECinfo%output,'(1X,a,f16.10)') 'FOP Init: Lagrangian Fragment energy       :', &
         & AtomicFragment%LagFOP
    write(DECinfo%output,'(1X,a,f16.10)') 'FOP Init: Occupied Fragment energy         :', &
         & AtomicFragment%EoccFOP
    write(DECinfo%output,'(1X,a,f16.10)') 'FOP Init: Virtual Fragment energy          :', &
         & AtomicFragment%EvirtFOP
    write(DECinfo%output,'(1X,a,i4)')    'FOP Init: Number of basis functions        :', &
         & AtomicFragment%number_basis
    write(DECinfo%output,'(1X,a,/)') 'FOP========================================================='
    write(DECinfo%output,*) 'FOP'



    ! ======================================================================
    !                             Expansion loop
    ! ======================================================================

    expansion_converged=.false.
    lag_converged=.false.
    occ_converged=.false.
    virt_converged=.false.
    EXPANSION_LOOP: do iter = 1,DECinfo%MaxIter
       OccOld=Occ_atoms;VirtOld=Virt_atoms
       LagEnergyOld = AtomicFragment%LagFOP
       OccEnergyOld = AtomicFragment%EoccFOP
       VirtEnergyOld = AtomicFragment%EvirtFOP

       call Expandfragment(Occ_atoms,Virt_atoms,DistTrackMyAtom,natoms,&
            & nocc_per_atom,nunocc_per_atom)
       call atomic_fragment_free(AtomicFragment)
       call get_fragment_and_Energy(MyAtom,natoms,Occ_Atoms,Virt_Atoms,&
            & MyMolecule,MyLsitem,nocc,nunocc,OccOrbitals,UnoccOrbitals,AtomicFragment)
       ! MPI fragment statistics
       slavetime = slavetime +  AtomicFragment%slavetime
       flops_slaves = flops_slaves + AtomicFragment%flops_slaves
       LagEnergyDiff=abs(LagEnergyOld-AtomicFragment%LagFOP)
       OccEnergyDiff=abs(OccEnergyOld-AtomicFragment%EoccFOP)
       VirtEnergyDiff=abs(VirtEnergyOld-AtomicFragment%EvirtFOP)
       call PrintInfo_Lagrangian(AtomicFragment,LagEnergyDiff,OccEnergyDiff,VirtEnergyDiff,iter)

       ! Test convergence for both Lagrangian, occupied, and virtual energies
       ! ********************************************************************

       ! Lagrangian
       TEST_CONVERGENCE_LAG: if  (LagEnergyDiff < FOT) then
          write(DECinfo%output,'(1X,a,F14.9)') 'FOP: Lagrangian energy converged, energydiff =', LagEnergyDiff
          lag_converged=.true.
       else
          write(DECinfo%output,*) 'FOP: Lagrangian energy NOT converged'
          lag_converged=.false.
       end if TEST_CONVERGENCE_LAG

       ! Occupied
       TEST_CONVERGENCE_OCC: if  (OccEnergyDiff < FOT) then
          write(DECinfo%output,'(1X,a,F14.9)') 'FOP: Occupied energy converged, energydiff   =', OccEnergyDiff
          occ_converged=.true.
       else
          write(DECinfo%output,*) 'FOP: Occupied energy NOT converged'
          occ_converged=.false.
       end if TEST_CONVERGENCE_OCC

       ! Virtual
       TEST_CONVERGENCE_VIRT: if  (VirtEnergyDiff < FOT) then
          write(DECinfo%output,'(1X,a,F14.9)') 'FOP: Virtual energy converged, energydiff    =', VirtEnergyDiff
          virt_converged=.true.
       else
          write(DECinfo%output,*) 'FOP: Virtual energy NOT converged'
          virt_converged=.false.
       end if TEST_CONVERGENCE_VIRT

       ! We are converged only if ALL three energies are converged
       ExpansionConvergence: if(lag_converged .and. occ_converged .and. virt_converged) then
          expansion_converged=.true.
          Occ_atoms = OccOld;Virt_atoms = VirtOld
          write(DECinfo%output,*) 'FOP Fragment expansion converged in iteration ', iter
          exit
       end if ExpansionConvergence

    end do EXPANSION_LOOP


    ! Check that expansion loop is converged
    if(.not. expansion_converged) then
       write(DECinfo%output,*) 'Number of expansion steps = ', DECinfo%MaxIter
       call lsquit('Lagrangian fragment expansion did not converge! &
            & Try to increase the number of expansion steps using the .MaxIter keyword',DECinfo%output)
    end if


    ! IMPORTANT - TO AVOID CONFUSION
    ! ******************************
    ! The size of the optimal fragment form the expansion loop is defined
    ! by the atoms in occ_atoms and virtatoms, and the energies for this fragment are stored
    ! in LagEnergyOld, OccEnergyOld, and VirtEnergyOld.
    ! However, the fragment information currently stored in "AtomicFragment" describes a fragment
    ! which is LARGER than the fragment defined by occ_atoms and virtatoms - in fact, the fragment
    ! expansion procedure concluded that this fragment was larger than required by the predefined
    ! accuracy (FOT).
    ! When we reduce fragments below we consider contributions from the individual orbitals
    ! to the Lagrangian fragment energy. The information stored in AtomicFragment
    ! contains more (and slightly more accurate) contributions from the individual orbitals
    ! than the current optimal fragment (defined by occ_atoms and virtatoms).
    ! Therefore, we might as well use the more accurate estimates contained in AtomicFragment
    ! when trying to reduce the fragment size below.
    call mem_alloc(OccContribs,nocc)
    call mem_alloc(VirtContribs,nunocc)
    OccContribs=0.0E0_realk
    VirtContribs=0.0E0_realk

    ! Contributions from occupied orbitals
    do i=1,AtomicFragment%noccAOS
       ! index of occupied AOS orbital "i" in list of ALL occupied orbitals in the molecule
       idx=AtomicFragment%occAOSidx(i)
       OccContribs(idx) = AtomicFragment%OccContribs(i)
    end do

    ! Contributions from virtual orbitals
    do i=1,AtomicFragment%nunoccAOS
       ! index of virtual AOS orbital "i" in list of ALL virtual orbitals in the molecule
       idx=AtomicFragment%unoccAOSidx(i)
       VirtContribs(idx) = AtomicFragment%VirtContribs(i)
    end do


    ! Set AtomicFragment to be the current optimal fragment
    call atomic_fragment_free(AtomicFragment)
    call atomic_fragment_init_atom_specific(MyAtom,natoms,Virt_Atoms, &
         & Occ_Atoms,nocc,nunocc,OccOrbitals,UnoccOrbitals, &
         & MyMolecule,mylsitem,AtomicFragment,.true.,.false.)
    ! Set energies correctly without having to do a new calculation
    AtomicFragment%LagFOP = LagEnergyOld
    AtomicFragment%EoccFOP = OccEnergyOld
    AtomicFragment%EvirtFOP = VirtEnergyOld


    ! Initialize logical vectors controlling occupied and virtual AOS during reduction scheme
    call mem_alloc(OccAOS_old,nocc)
    call mem_alloc(VirtAOS_old,nunocc)
    call mem_alloc(OccAOS_new,nocc)
    call mem_alloc(VirtAOS_new,nunocc)
    OccAOS_old=.false.
    VirtAOS_old=.false.
    OccAOS_new=.false.
    VirtAOS_new=.false.
    do i=1,AtomicFragment%noccAOS
       OccAOS_old(AtomicFragment%occAOSidx(i)) = .true.
    end do
    do i=1,AtomicFragment%nunoccAOS
       VirtAOS_old(AtomicFragment%unoccAOSidx(i)) = .true.
    end do


    write(DECinfo%output,*) ' FOP'
    write(DECinfo%output,*) ' FOP *************************************************'
    write(DECinfo%output,*) ' FOP ** Expansion has converged. We start reduction **'
    write(DECinfo%output,*) ' FOP *************************************************'
    write(DECinfo%output,*) ' FOP'


    ! ======================================================================
    !                             Reduction loop
    ! ======================================================================


    lag_converged=.false.
    occ_converged=.false.
    virt_converged=.false.
    converged=.false.
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

       ! Cycle if no atoms were rejected by the rejection procedure
       if ((nocc_old == nocc_new) .and. (nvirt_old == nvirt_new)) then
          if(iter == max_iter_red) then

             write(DECinfo%output,*) "FOP No reduction possible. Use original converged fragment"

             ! Go back to old fragment
             ! ***********************
             Occ_atoms=OccOld
             Virt_atoms=VirtOld
             call atomic_fragment_free(AtomicFragment)
             call get_fragment_and_Energy(MyAtom,natoms,Occ_Atoms,Virt_Atoms,&
                     & MyMolecule,MyLsitem,nocc,nunocc,OccOrbitals,UnoccOrbitals,&
                     & AtomicFragment)
             ! MPI fragment statistics
             slavetime = slavetime +  AtomicFragment%slavetime
             flops_slaves = flops_slaves + AtomicFragment%flops_slaves
             converged = .true.

          else
             write(DECinfo%output,*) 'FOP No reduction occurred - we try to reduce again'
             cycle
          end if
       end if

       ! Get reduced atomic fragment
       call atomic_fragment_free(AtomicFragment)
       call atomic_fragment_init_orbital_specific(MyAtom,nunocc, nocc, VirtAOS_new, &
            & OccAOS_new,OccOrbitals,UnoccOrbitals,MyMolecule,mylsitem,AtomicFragment,.true.,.false.)

       call single_lagrangian_energy_and_prop(AtomicFragment)

       ! MPI fragment statistics
       slavetime = slavetime +  AtomicFragment%slavetime
       flops_slaves = flops_slaves + AtomicFragment%flops_slaves

       LagEnergyDiff=abs(AtomicFragment%LagFOP-LagEnergyOld)
       OccEnergyDiff=abs(AtomicFragment%EoccFOP-OccEnergyOld)
       VirtEnergyDiff=abs(AtomicFragment%EvirtFOP-VirtEnergyOld)
       call PrintInfo_Lagrangian(AtomicFragment,LagEnergyDiff,OccEnergyDiff,VirtEnergyDiff,iter)


       ! Test convergence for both Lagrangian, occupied, and virtual energies
       ! ********************************************************************

       ! Lagrangian
       TEST_REDUCTION_LAG: if  (LagEnergyDiff < FOT) then
          write(DECinfo%output,'(1X,a,F14.9)') 'FOP: Lagrangian reduction converged, energydiff =', LagEnergyDiff
          lag_converged=.true.
       else
          write(DECinfo%output,*) 'FOP: Lagrangian reduction NOT converged'
          lag_converged=.false.
       end if TEST_REDUCTION_LAG

       ! Occupied
       TEST_REDUCTION_OCC: if  (OccEnergyDiff < FOT) then
          write(DECinfo%output,'(1X,a,F14.9)') 'FOP: Occupied reduction converged, energydiff   =', OccEnergyDiff
          occ_converged=.true.
       else
          write(DECinfo%output,*) 'FOP: Occupied reduction NOT converged'
          occ_converged=.false.
       end if TEST_REDUCTION_OCC

       ! Virtual
       TEST_REDUCTION_VIRT: if  (VirtEnergyDiff < FOT) then
          write(DECinfo%output,'(1X,a,F14.9)') 'FOP: Virtual reduction converged, energydiff    =', VirtEnergyDiff
          virt_converged=.true.
       else
          write(DECinfo%output,*) 'FOP: Virtual reduction NOT converged'
          virt_converged=.false.
       end if TEST_REDUCTION_VIRT
       if (lag_converged .and. occ_converged .and. virt_converged) converged=.true.


       ! Test for fragment of lower accuracy REDFOT=10*FOT
       ! *************************************************
       if(.not. REDconverged) then
          if( (LagEnergyDiff<REDFOT) .and. (OccEnergyDiff<REDFOT) .and. (VirtEnergyDiff<REDFOT) ) then
             REDconverged=.true.
             write(DECinfo%output,*) 'FOP Lower-accuracy fragment converged in step', iter

             ! Save information for reduced fragment
             REDocc = AtomicFragment%noccAOS
             REDvirt = AtomicFragment%nunoccAOS
             call mem_alloc(REDoccIDX,REDocc)
             call mem_alloc(REDvirtIDX,REDvirt)
             do i=1,REDocc
                REDoccIDX(i) = AtomicFragment%occAOSidx(i)
             end do
             do i=1,REDvirt
                REDvirtIDX(i) = AtomicFragment%unoccAOSidx(i)
             end do
          end if
       end if


       ! Quit if fragment reduction is converged
       ! ***************************************
       if(converged) then
          write(DECinfo%output,*) 'FOP Reduction of fragment converged in step', iter
          exit
       else if ((iter == max_iter_red)) then
          write(DECinfo%output,*) "FOP No reduction possible. Use original converged fragment"

          ! Go back to old fragment
          ! ***********************
          Occ_atoms=OccOld
          Virt_atoms=VirtOld
          call atomic_fragment_free(AtomicFragment)
          call get_fragment_and_Energy(MyAtom,natoms,Occ_Atoms,Virt_Atoms,&
               & MyMolecule,MyLsitem,nocc,nunocc,OccOrbitals,UnoccOrbitals,&
               & AtomicFragment)
          ! MPI fragment statistics
          slavetime = slavetime +  AtomicFragment%slavetime
          flops_slaves = flops_slaves + AtomicFragment%flops_slaves
          converged = .true.
       end if


    end do REDUCTION_LOOP


    if (converged) then



       ! Set info for reduced fragment of lower accuracy
       ! ***********************************************
       call mem_dealloc(AtomicFragment%REDoccAOSidx)
       call mem_dealloc(AtomicFragment%REDunoccAOSidx)
       AtomicFragment%REDnoccAOS = REDocc
       AtomicFragment%REDnunoccAOS = REDvirt
       call mem_alloc(AtomicFragment%REDoccAOSidx,REDocc)
       call mem_alloc(AtomicFragment%REDunoccAOSidx,REDvirt)
       AtomicFragment%REDoccAOSidx = REDoccIDX
       AtomicFragment%REDunoccAOSidx = REDvirtIDX


       ! Update t1 amplitudes for full molecule
       ! **************************************

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
       write(DECinfo%output,'(1X,a,i4)') 'FOP      SITE OPTIMIZATION HAS CONVERGED FOR SITE',MyAtom
       write(DECinfo%output,'(1X,a)') 'FOP---------------------------------------------------------'
       write(DECinfo%output,'(1X,a,i4)')    'FOP Done: Fragment number                  :', MyAtom
       write(DECinfo%output,'(1X,a,i4)')    'FOP Done: Number of orbitals in virt total :', &
            & AtomicFragment%nunoccAOS
       write(DECinfo%output,'(1X,a,i4)')    'FOP Done: Number of orbitals in occ total  :', &
            & AtomicFragment%noccAOS
       write(DECinfo%output,'(1X,a,f16.10)') 'FOP Done: Lagrangian Fragment energy       :', &
            & AtomicFragment%LagFOP
       write(DECinfo%output,'(1X,a,f16.10)') 'FOP Done: Occupied Fragment energy         :', &
            & AtomicFragment%EoccFOP
       write(DECinfo%output,'(1X,a,f16.10)') 'FOP Done: Virtual Fragment energy          :', &
            & AtomicFragment%EvirtFOP
       write(DECinfo%output,'(1X,a,i4)')    'FOP Done: Number of basis functions        :', &
            & AtomicFragment%number_basis
       write(DECinfo%output,'(1X,a,/)') 'FOP========================================================='
       write(DECinfo%output,*) 'FOP'
    else
       call lsquit('optimize_atomic_fragment_main: Optimization NOT converged',DECinfo%output)
    end if

    call mem_dealloc(OccContribs)
    call mem_dealloc(VirtContribs)
    call mem_dealloc(OccAOS_old)
    call mem_dealloc(VirtAOS_old)
    call mem_dealloc(OccAOS_new)
    call mem_dealloc(VirtAOS_new)
    call mem_dealloc(REDoccIDX)
    call mem_dealloc(REDvirtIDX)

    if(freebasisinfo) then
       call atomic_fragment_free_basis_info(AtomicFragment)
    end if

    ! MPI fragment statistics
    AtomicFragment%slavetime  = slavetime
    AtomicFragment%flops_slaves = flops_slaves

    ! Ensure that energies in fragment are set consistently
    call set_energies_ccatom_structure_fragopt(AtomicFragment)

  end subroutine optimize_atomic_fragment_main




  !> Routine that optimizes an atomic fragment by (1) expanding based on orbital interaction matrix,
  !> (2) checking that energy change is smaller than FOT, (3) for converged fragment we 
  !> make unitary transformations of occupied and virtual AOS (NOT including EOS) to 
  !> generate fragment-adapted orbitals which describe AOS with as few orbitals as possible.
  !> \date February 2013
  !> \author Ida-Marie Hoeyvik & Kasper Kristensen
!!$  subroutine optimize_atomic_fragment_FA_experimental(MyAtom,AtomicFragment,nAtoms, &
!!$       &OccOrbitals,nOcc,UnoccOrbitals,nUnocc,DistanceTable, &
!!$       &MyMolecule,mylsitem,freebasisinfo,t1,t2,t1full)
!!$    implicit none
!!$    !> Number of occupied orbitals in molecule
!!$    integer, intent(in) :: nOcc
!!$    !> Number of unoccupied orbitals in molecule
!!$    integer, intent(in) :: nunocc
!!$    !> Number of atoms in molecule
!!$    integer, intent(in) :: natoms
!!$    !> Central atom in molecule
!!$    integer, intent(in) :: MyAtom
!!$    !> Atomic fragment to be optimized
!!$    type(ccatom),intent(inout)        :: AtomicFragment
!!$    !> All occupied orbitals
!!$    type(ccorbital), dimension(nOcc), intent(in)      :: OccOrbitals
!!$    !> All unoccupied orbitals
!!$    type(ccorbital), dimension(nUnocc), intent(in)    :: UnoccOrbitals
!!$    !> Distance table for all atoms in molecule
!!$    real(realk), dimension(nAtoms,nAtoms), intent(in) :: DistanceTable
!!$    !> Full molecule information
!!$    type(fullmolecule), intent(in) :: MyMolecule
!!$    !> Integral information
!!$    type(lsitem), intent(inout)       :: mylsitem
!!$    !> Delete fragment basis information ("expensive box in ccatom type") at exit?
!!$    logical,intent(in) :: freebasisinfo
!!$    !> t1 amplitudes for AOS (not full molecule)
!!$    type(array2),intent(inout),optional :: t1
!!$    !> t2 amplitudes for AOS (not full molecule)
!!$    type(array4),intent(inout),optional :: t2
!!$    !> t1 amplitudes for full molecule to be updated (only used when DECinfo%SinglesPolari is set)
!!$    type(array2),intent(inout),optional :: t1full
!!$    real(realk)                    :: LagEnergyDiff, OccEnergyDiff,VirtEnergyDiff
!!$    real(realk)                    :: LagEnergyOld, OccEnergyOld, VirtEnergyOld
!!$    real(realk)                    :: FOT
!!$    integer,dimension(natoms)      :: nocc_per_atom,nunocc_per_atom
!!$    type(ccatom) :: FOfragment
!!$    integer      :: iter,ov,occsize,unoccsize
!!$    integer      :: Nnew,Nold, max_iter_red,nocc_exp,nvirt_exp
!!$    logical      :: converged,storeamp,ampset,ReductionPossible(2)
!!$    logical :: expansion_converged, lag_converged, occ_converged, virt_converged,allorbs
!!$    real(realk) :: slavetime, flops_slaves
!!$    real(realk),pointer :: OccMat(:,:), VirtMat(:,:),occval(:),unoccval(:)
!!$    logical,pointer :: OccAOS(:),UnoccAOS(:),OccOld(:),UnoccOld(:)
!!$    integer,pointer :: occidx(:),unoccidx(:)
!!$
!!$
!!$    write(DECinfo%output,'(a)')    ' FOP'
!!$    write(DECinfo%output,'(a)')    ' FOP ==============================================='
!!$    write(DECinfo%output,'(a,i4)') ' FOP  Site fragment generator for fragment,',MyAtom
!!$    write(DECinfo%output,'(a)')    ' FOP ==============================================='
!!$    write(DECinfo%output,'(a)')    ' FOP'
!!$
!!$
!!$    ! Sanity check for singles polarization
!!$    if(DECinfo%SinglesPolari) then
!!$       if(.not. present(t1full)) then
!!$          call lsquit('optimize_atomic_fragment: Full singles polarization is requrested &
!!$               & but t1 argument is not present!',DECinfo%output)
!!$       end if
!!$    end if
!!$
!!$
!!$    ! ******************************************
!!$    ! **   Initialization of stuff needed..   **
!!$    ! ******************************************
!!$    allorbs=.false.
!!$    ampset=.false.
!!$    storeamp=.false.
!!$    if(DECinfo%ccmodel==4 .and. present(t1) .and. present(t2) ) then 
!!$       storeamp=.true. ! store t1 and t2 amplitudes for (T) model
!!$    end if
!!$    converged=.false.
!!$    max_iter_red=15   ! allow 15 reduction steps (should be more than enough)
!!$    FOT = DECinfo%FOT
!!$    nocc_per_atom=get_number_of_orbitals_per_atom(OccOrbitals,nocc,natoms)
!!$    nunocc_per_atom=get_number_of_orbitals_per_atom(UnoccOrbitals,nunocc,natoms)
!!$
!!$    ! MPI fragment statistics
!!$    flops_slaves = 0.0E0_realk
!!$    slavetime = 0.0E0_realk
!!$
!!$    ! Sanity quit
!!$    if( (nocc_per_atom(MyAtom) == 0) .and. (nunocc_per_atom(MyAtom) == 0) ) then
!!$       write(DECinfo%output,*) 'FOP Skipping optimization of fragment ', MyAtom
!!$       AtomicFragment%LagFOP=0E0_realk
!!$       AtomicFragment%EoccFOP=0E0_realk
!!$       AtomicFragment%EvirtFOP=0E0_realk
!!$       return
!!$    end if
!!$
!!$    ! ******************************************
!!$    ! **  Starting computation of fragment    **
!!$    ! ******************************************
!!$    if(DECinfo%simulate_full .or. DECinfo%InclFullMolecule) then
!!$       ! skip calculation of energy and just initiate fragment including ALL orbitals
!!$       write(DECinfo%output,*) 'FOP Fragment includes all orbitals and fragment optimization is skipped'
!!$       call fragment_init_simulate_full(MyAtom,nunocc, nocc, OccOrbitals,UnoccOrbitals,&
!!$            & MyMolecule,mylsitem,AtomicFragment,.true.)
!!$
!!$       if(DECinfo%fragadapt) then
!!$          write(DECinfo%output,*) 'WARNING! For fragment-adapted orbitals we cannot skip energy calculation! '
!!$          call single_lagrangian_energy_and_prop(AtomicFragment)
!!$       end if
!!$
!!$       call PrintInfo_Lagrangian(AtomicFragment,0.0E0_realk,0.0E0_realk,0.0E0_realk,0)
!!$       if(freebasisinfo) then
!!$          call atomic_fragment_free_basis_info(AtomicFragment)
!!$       end if
!!$       return
!!$
!!$    else  ! regular fragment optimization
!!$
!!$       call mem_alloc(OccOld,nocc)
!!$       call mem_alloc(UnoccOld,nunocc)
!!$       call mem_alloc(OccAOS,nocc)
!!$       call mem_alloc(UnoccAOS,nunocc)
!!$       call mem_alloc(occidx,nocc)
!!$       call mem_alloc(unoccidx,nunocc)
!!$       call mem_alloc(occval,nocc)
!!$       call mem_alloc(unoccval,nunocc)
!!$
!!$       ! Prioritize all orbitals based on interactions with orbitals assigned to MyAtom
!!$       call prioritize_orbitals_based_on_atom(MyAtom,MyMolecule,OccOrbitals,UnoccOrbitals,&
!!$            & occidx,unoccidx,occval,unoccval)
!!$
!!$       ! Number of orbitals to include in expansion loop (under investigation)
!!$       call get_expansion_parameters_for_fragopt(MyMolecule,nocc_per_atom(MyAtom),&
!!$            & nunocc_per_atom(MyAtom),occsize,unoccsize)
!!$
!!$       ! Get initial fragment logical vectors
!!$       call init_fragment_logical_vectors_from_orbital_interactions(MyAtom,&
!!$            & nocc,nunocc,OccOrbitals,UnoccOrbitals,occidx,unoccidx,MyMolecule,&
!!$            & occsize,unoccsize,occAOS,unoccAOS)
!!$
!!$       ! Init fragment type based on logical vectors and calculate associated fragment energies
!!$       call get_fragment_and_Energy_orb_specific(MyAtom,MyMolecule,mylsitem,&
!!$            & OccOrbitals,UnoccOrbitals,OccAOS,UnoccAOS, AtomicFragment)
!!$
!!$       ! MPI fragment statistics
!!$       slavetime = slavetime +  AtomicFragment%slavetime
!!$       flops_slaves = flops_slaves + AtomicFragment%flops_slaves
!!$
!!$    end if
!!$
!!$    ! Print initial fragment information
!!$    write(DECinfo%output,*)'FOP'
!!$    write(DECinfo%output,'(1X,a)') 'FOP========================================================='
!!$    write(DECinfo%output,'(1X,a)') 'FOP               Initial fragment information'
!!$    write(DECinfo%output,'(1X,a)') 'FOP---------------------------------------------------------'
!!$    write(DECinfo%output,'(1X,a,i4)')    'FOP Init: Fragment number                  :', MyAtom
!!$    write(DECinfo%output,'(1X,a,i4)')    'FOP Init: Number of orbitals in virt total :', &
!!$         & AtomicFragment%nunoccAOS
!!$    write(DECinfo%output,'(1X,a,i4)')    'FOP Init: Number of orbitals in occ total  :', &
!!$         & AtomicFragment%noccAOS
!!$    write(DECinfo%output,'(1X,a,f16.10)') 'FOP Init: Lagrangian Fragment energy       :', &
!!$         & AtomicFragment%LagFOP
!!$    write(DECinfo%output,'(1X,a,f16.10)') 'FOP Init: Occupied Fragment energy         :', &
!!$         & AtomicFragment%EoccFOP
!!$    write(DECinfo%output,'(1X,a,f16.10)') 'FOP Init: Virtual Fragment energy          :', &
!!$         & AtomicFragment%EvirtFOP
!!$    write(DECinfo%output,'(1X,a,i4)')    'FOP Init: Number of basis functions        :', &
!!$         & AtomicFragment%number_basis
!!$    write(DECinfo%output,'(1X,a,/)') 'FOP========================================================='
!!$    write(DECinfo%output,*) 'FOP'
!!$
!!$
!!$
!!$    ! ======================================================================
!!$    !                             Expansion loop
!!$    ! ======================================================================
!!$
!!$    expansion_converged=.false.
!!$    lag_converged=.false.
!!$    occ_converged=.false.
!!$    virt_converged=.false.
!!$    EXPANSION_LOOP: do iter = 1,DECinfo%MaxIter
!!$
!!$       ! Save current AOS vectors
!!$       OccOld=OccAOS
!!$       UnoccOld=UnoccAOS
!!$       LagEnergyOld = AtomicFragment%LagFOP
!!$       OccEnergyOld = AtomicFragment%EoccFOP
!!$       VirtEnergyOld = AtomicFragment%EvirtFOP
!!$
!!$       ! Save correlation density matrices
!!$       call mem_alloc(OccMat,AtomicFragment%noccAOS,AtomicFragment%noccAOS)
!!$       call mem_alloc(VirtMat,AtomicFragment%nunoccAOS,AtomicFragment%nunoccAOS)
!!$       OccMat = AtomicFragment%OccMat
!!$       VirtMat = AtomicFragment%VirtMat
!!$
!!$
!!$       ! Expand fragment!
!!$       ! ----------------
!!$       ! Expand logical vectors
!!$       call expand_fragment_logical_vectors_from_orbital_interactions(nocc,nunocc,&
!!$            & occidx,unoccidx,MyMolecule,occsize,unoccsize,occAOS,unoccAOS)
!!$       ! Free existing fragment
!!$       call atomic_fragment_free(AtomicFragment)
!!$       ! Init fragment structure with new AOS and calculate associated fragment energies
!!$       call get_fragment_and_Energy_orb_specific(MyAtom,MyMolecule,mylsitem,&
!!$            & OccOrbitals,UnoccOrbitals,OccAOS,UnoccAOS, AtomicFragment)
!!$       ! MPI fragment statistics
!!$       slavetime = slavetime +  AtomicFragment%slavetime
!!$       flops_slaves = flops_slaves + AtomicFragment%flops_slaves
!!$       LagEnergyDiff=abs(LagEnergyOld-AtomicFragment%LagFOP)
!!$       OccEnergyDiff=abs(OccEnergyOld-AtomicFragment%EoccFOP)
!!$       VirtEnergyDiff=abs(VirtEnergyOld-AtomicFragment%EvirtFOP)
!!$       call PrintInfo_Lagrangian(AtomicFragment,LagEnergyDiff,OccEnergyDiff,VirtEnergyDiff,iter)
!!$
!!$       ! Test convergence for both Lagrangian, occupied, and virtual energies
!!$       ! ********************************************************************
!!$
!!$       ! Lagrangian
!!$       TEST_CONVERGENCE_LAG: if  (LagEnergyDiff < FOT) then
!!$          write(DECinfo%output,'(1X,a,F14.9)') 'FOP: Lagrangian energy converged, energydiff =', LagEnergyDiff
!!$          lag_converged=.true.
!!$       else
!!$          write(DECinfo%output,*) 'FOP: Lagrangian energy NOT converged'
!!$          lag_converged=.false.
!!$       end if TEST_CONVERGENCE_LAG
!!$
!!$       ! Occupied
!!$       TEST_CONVERGENCE_OCC: if  (OccEnergyDiff < FOT) then
!!$          write(DECinfo%output,'(1X,a,F14.9)') 'FOP: Occupied energy converged, energydiff   =', OccEnergyDiff
!!$          occ_converged=.true.
!!$       else
!!$          write(DECinfo%output,*) 'FOP: Occupied energy NOT converged'
!!$          occ_converged=.false.
!!$       end if TEST_CONVERGENCE_OCC
!!$
!!$       ! Virtual
!!$       TEST_CONVERGENCE_VIRT: if  (VirtEnergyDiff < FOT) then
!!$          write(DECinfo%output,'(1X,a,F14.9)') 'FOP: Virtual energy converged, energydiff    =', VirtEnergyDiff
!!$          virt_converged=.true.
!!$       else
!!$          write(DECinfo%output,*) 'FOP: Virtual energy NOT converged'
!!$          virt_converged=.false.
!!$       end if TEST_CONVERGENCE_VIRT
!!$
!!$       ! We are converged only if ALL three energies are converged
!!$       ExpansionConvergence: if(lag_converged .and. occ_converged .and. virt_converged) then
!!$          expansion_converged=.true.
!!$          OccAOS = OccOld
!!$          UnoccAOS = UnoccOld
!!$          write(DECinfo%output,'(a,i6)') 'FOP Fragment expansion converged in iteration ', iter
!!$          exit
!!$       end if ExpansionConvergence
!!$
!!$       ! Special case for small molecules - is the whole molecule included?
!!$       if(AtomicFragment%noccAOS==nocc .and. AtomicFragment%nunoccAOS==nunocc) then
!!$          expansion_converged=.true.
!!$          allorbs=.true.   ! all orbitals included in fragment!
!!$          write(DECinfo%output,'(a,i6)') 'FOP The full molecule has been included in iteration ', iter
!!$          exit
!!$       end if
!!$
!!$       call mem_dealloc(OccMat)
!!$       call mem_dealloc(VirtMat)
!!$
!!$    end do EXPANSION_LOOP
!!$
!!$
!!$    ! Check that expansion loop is converged
!!$    if(.not. expansion_converged) then
!!$       write(DECinfo%output,*) 'Number of expansion steps = ', DECinfo%MaxIter
!!$       call lsquit('Lagrangian fragment expansion did not converge! &
!!$            & Try to increase the number of expansion steps using the .MaxIter keyword',DECinfo%output)
!!$    end if
!!$
!!$
!!$    ! Set AtomicFragment to be the converged fragment
!!$    ! --> but do not do this for special case where current 
!!$    !     fragment is already the whole molecule
!!$    if(.not. allorbs) then
!!$       call atomic_fragment_free(AtomicFragment)
!!$       call atomic_fragment_init_orbital_specific(MyAtom,nunocc,nocc, UnoccAOS, &
!!$            & occAOS,OccOrbitals,UnoccOrbitals,MyMolecule,mylsitem,AtomicFragment,.true.,.false.)
!!$
!!$       ! Set energies and correlation densities correctly without having to do a new calculation
!!$       AtomicFragment%LagFOP = LagEnergyOld
!!$       AtomicFragment%EoccFOP = OccEnergyOld
!!$       AtomicFragment%EvirtFOP = VirtEnergyOld
!!$
!!$       ! correlation density matrices
!!$       call mem_alloc(AtomicFragment%OccMat,AtomicFragment%noccAOS,AtomicFragment%noccAOS)
!!$       call mem_alloc(AtomicFragment%VirtMat,AtomicFragment%nunoccAOS,AtomicFragment%nunoccAOS)
!!$       AtomicFragment%CDSet=.true. ! correlation density matrices have been set
!!$       AtomicFragment%OccMat = OccMat
!!$       AtomicFragment%VirtMat = VirtMat
!!$    end if
!!$       
!!$    call mem_dealloc(OccMat)
!!$    call mem_dealloc(VirtMat)
!!$
!!$    ! Save dimensions for statistics
!!$    nocc_exp = AtomicFragment%noccAOS
!!$    nvirt_exp = AtomicFragment%nunoccAOS
!!$
!!$
!!$    write(DECinfo%output,*) ' FOP'
!!$    write(DECinfo%output,*) ' FOP *************************************************'
!!$    write(DECinfo%output,*) ' FOP ** Expansion has converged. We start reduction **'
!!$    write(DECinfo%output,*) ' FOP *************************************************'
!!$    write(DECinfo%output,*) ' FOP'
!!$
!!$
!!$
!!$    ! ======================================================================
!!$    !                             Reduction loop
!!$    ! ======================================================================
!!$    ! Here we carry out unitary transformation of occupied and virtual AOSs 
!!$    ! (NOT including EOS orbitals) to adapt the AOS to the EOS orbitals.
!!$    ! We term the new set of orbitals "fragment-adapted orbitals".
!!$    ReductionPossible=.false.
!!$
!!$
!!$    ! First reduce virtual space (ov=2) then occupied space (ov=1)
!!$    OCC_OR_VIRT: do ov=2,1,-1   
!!$
!!$       lag_converged=.false.
!!$       occ_converged=.false.
!!$       virt_converged=.false.
!!$       converged=.false.
!!$
!!$       REDUCTION_LOOP: do iter=1,max_iter_red
!!$
!!$          if(ov==1) then
!!$             write(DECinfo%output,*) 'FOP Starting occ reduction step ', iter
!!$          else
!!$             write(DECinfo%output,*) 'FOP Starting virt reduction step ', iter
!!$          end if
!!$
!!$          if(iter == 1) then
!!$             ! Number of occ or virt orbitals in converged fragment
!!$             if(ov==1) then
!!$                Nold = AtomicFragment%noccAOS
!!$             else
!!$                Nold = AtomicFragment%nunoccAOS
!!$             end if
!!$             ! Threshold for throwing away fragment-adapted orbitals (start with FOT value)
!!$             AtomicFragment%RejectThr(ov) = FOT
!!$             if(ov==1) call atomic_fragment_free(FOfragment)
!!$          else  
!!$             ! Number of occ or virt fragment-adapted orbitals from previous step
!!$             if(ov==1) then
!!$                Nold = FOfragment%noccAOS
!!$             else
!!$                Nold = FOfragment%nunoccAOS
!!$             end if
!!$             ! decrease rejection threshold by a factor 10 in each step
!!$             AtomicFragment%RejectThr(ov) = AtomicFragment%RejectThr(ov)/10.0_realk
!!$             call atomic_fragment_free(FOfragment)
!!$          end if
!!$          if(ov==1) then
!!$             write(DECinfo%output,'(a,ES13.5)') 'FOP occ rejection threshold  : ',&
!!$                  & AtomicFragment%RejectThr(ov)
!!$          else
!!$             write(DECinfo%output,'(a,ES13.5)') 'FOP virt rejection threshold  : ',&
!!$                  & AtomicFragment%RejectThr(ov)
!!$          end if
!!$
!!$
!!$          ! Make fragment-adapted fragment with smaller AOS according to rejection threshold
!!$          call fragment_adapted_driver(MyMolecule,mylsitem,OccOrbitals,UnoccOrbitals,&
!!$               & AtomicFragment,FOfragment)
!!$          if(ov==1) then
!!$             Nnew = FOfragment%noccAOS
!!$          else
!!$             Nnew = FOfragment%nunoccAOS
!!$          end if
!!$
!!$          ! Cycle if the AOS was not decreased in the rejection procedure
!!$          if( Nold==Nnew ) then
!!$
!!$             ! special case if max number of iterations was reached
!!$             if(iter == max_iter_red) then   
!!$                if(ov==2) then  ! not possible to reduce virtual space, try occupied space
!!$                   cycle OCC_OR_VIRT
!!$
!!$                else ! not possible to decrease occ space
!!$                   goto 1000
!!$                end if
!!$             else
!!$                write(DECinfo%output,*) 'FOP No reduction occurred - try again...'
!!$                cycle REDUCTION_LOOP
!!$             end if
!!$
!!$          end if
!!$          if(storeamp) then
!!$             if(ampset) then  ! free existing amplitudes before calculating new ones
!!$                call array2_free(t1)
!!$                call array4_free(t2)
!!$             end if
!!$             call single_lagrangian_energy_and_prop(FOfragment,t1=t1,t2=t2)
!!$             ampset=.true.
!!$          else
!!$             call single_lagrangian_energy_and_prop(FOfragment)
!!$          end if
!!$
!!$
!!$          ! MPI fragment statistics
!!$          slavetime = slavetime +  FOfragment%slavetime
!!$          flops_slaves = flops_slaves + FOfragment%flops_slaves
!!$
!!$          LagEnergyDiff=abs(FOfragment%LagFOP-LagEnergyOld)
!!$          OccEnergyDiff=abs(FOfragment%EoccFOP-OccEnergyOld)
!!$          VirtEnergyDiff=abs(FOfragment%EvirtFOP-VirtEnergyOld)
!!$          call PrintInfo_Lagrangian(FOfragment,LagEnergyDiff,OccEnergyDiff,VirtEnergyDiff,iter)
!!$
!!$
!!$          ! Test convergence for both Lagrangian, occupied, and virtual energies
!!$          ! ********************************************************************
!!$
!!$          ! Lagrangian
!!$          TEST_REDUCTION_LAG: if  (LagEnergyDiff < FOT) then
!!$             if(ov==1) then
!!$                write(DECinfo%output,'(1X,a,F14.9)') &
!!$                     & 'FOP: Lagrangian energy converged (occ reduction), diff =', LagEnergyDiff
!!$             else
!!$                write(DECinfo%output,'(1X,a,F14.9)') &
!!$                     & 'FOP: Lagrangian energy converged (virt reduction), diff =', LagEnergyDiff
!!$             end if
!!$             lag_converged=.true.
!!$          else
!!$             if(ov==1) then
!!$                write(DECinfo%output,*) 'FOP: Lagrangian energy (occ reduction) NOT converged'
!!$             else
!!$                write(DECinfo%output,*) 'FOP: Lagrangian energy (virt reduction) NOT converged'
!!$             end if
!!$             lag_converged=.false.
!!$          end if TEST_REDUCTION_LAG
!!$
!!$          ! Occupied
!!$          TEST_REDUCTION_OCC: if  (OccEnergyDiff < FOT) then
!!$             if(ov==1) then
!!$                write(DECinfo%output,'(1X,a,F14.9)') &
!!$                     & 'FOP: Occupied energy converged (occ reduction), diff =', OccEnergyDiff
!!$             else
!!$                write(DECinfo%output,'(1X,a,F14.9)') &
!!$                     & 'FOP: Occupied energy converged (virt reduction), diff =', OccEnergyDiff
!!$             end if
!!$             occ_converged=.true.
!!$          else
!!$             if(ov==1) then
!!$                write(DECinfo%output,*) 'FOP: Occupied energy (occ reduction) NOT converged'
!!$             else
!!$                write(DECinfo%output,*) 'FOP: Occupied energy (virt reduction) NOT converged'
!!$             end if
!!$             occ_converged=.false.
!!$          end if TEST_REDUCTION_OCC
!!$
!!$          ! Virtual
!!$          TEST_REDUCTION_VIRT: if  (VirtEnergyDiff < FOT) then
!!$             if(ov==1) then
!!$                write(DECinfo%output,'(1X,a,F14.9)') &
!!$                     & 'FOP: Virtual energy converged (occ reduction), diff =', VirtEnergyDiff
!!$             else
!!$                write(DECinfo%output,'(1X,a,F14.9)') &
!!$                     & 'FOP: Virtual energy converged (virt reduction), diff =', VirtEnergyDiff
!!$             end if
!!$             virt_converged=.true.
!!$          else
!!$             if(ov==1) then
!!$                write(DECinfo%output,*) 'FOP: Virtual energy (occ reduction) NOT converged'
!!$             else
!!$                write(DECinfo%output,*) 'FOP: Virtual energy (virt reduction) NOT converged'
!!$             end if
!!$             virt_converged=.false.
!!$          end if TEST_REDUCTION_VIRT
!!$          if (lag_converged .and. occ_converged .and. virt_converged) then
!!$             converged=.true.
!!$             ReductionPossible(ov)=.true.
!!$          end if
!!$
!!$
!!$
!!$          ! Quit if fragment reduction is converged
!!$          ! ***************************************
!!$          if(converged) then
!!$
!!$             if(ov==2) then  ! virtual reduction converged 
!!$                write(DECinfo%output,*) 'FOP Virt reduction of fragment converged in step', iter
!!$                AtomicFragment%LagFOP = FOfragment%LagFOP
!!$                AtomicFragment%EoccFOP = FOfragment%EoccFOP
!!$                AtomicFragment%EvirtFOP = FOfragment%EvirtFOP
!!$                cycle OCC_OR_VIRT  ! try to reduce occ space
!!$
!!$             else ! both occ and virtual reductions have converged
!!$                write(DECinfo%output,*) 'FOP Occ  reduction of fragment converged in step', iter
!!$                AtomicFragment%LagFOP = FOfragment%LagFOP
!!$                AtomicFragment%EoccFOP = FOfragment%EoccFOP
!!$                AtomicFragment%EvirtFOP = FOfragment%EvirtFOP
!!$                exit
!!$             end if
!!$
!!$          end if
!!$
!!$1000      continue
!!$
!!$          if ( (iter == max_iter_red) .and. (ov==1) .and. &
!!$               & (.not. any(ReductionPossible) ) ) then
!!$             write(DECinfo%output,*) "FOP No reduction possible. Use original converged fragment"
!!$
!!$             ! Go back to old fragment
!!$             ! ***********************
!!$             ! Recalculate! This is extra work but it will never be done in practice, only
!!$             ! relevant for small debug systems.
!!$             occAOS=.true.
!!$             unoccAOS=.true.
!!$             call atomic_fragment_free(AtomicFragment)
!!$             if(storeamp) then
!!$                if(ampset) then  ! free existing amplitudes before calculating new ones
!!$                   call array2_free(t1)
!!$                   call array4_free(t2)
!!$                end if
!!$                call get_fragment_and_Energy_orb_specific(MyAtom,MyMolecule,mylsitem,&
!!$                     & OccOrbitals,UnoccOrbitals,OccAOS,UnoccAOS, AtomicFragment,t1=t1,t2=t2)
!!$                ampset=.true.
!!$             else
!!$                call get_fragment_and_Energy_orb_specific(MyAtom,MyMolecule,mylsitem,&
!!$                     & OccOrbitals,UnoccOrbitals,OccAOS,UnoccAOS, AtomicFragment)
!!$             end if
!!$             call atomic_fragment_free(FOfragment)
!!$             call fragment_adapted_driver(MyMolecule,mylsitem,OccOrbitals,UnoccOrbitals,&
!!$                  & AtomicFragment,FOfragment)
!!$             ! MPI fragment statistics
!!$             slavetime = slavetime +  AtomicFragment%slavetime
!!$             flops_slaves = flops_slaves + AtomicFragment%flops_slaves
!!$          end if
!!$
!!$       end do REDUCTION_LOOP
!!$
!!$    end do OCC_OR_VIRT
!!$
!!$
!!$    call mem_dealloc(OccOld)
!!$    call mem_dealloc(UnoccOld)
!!$    call mem_dealloc(OccAOS)
!!$    call mem_dealloc(UnoccAOS)
!!$    call mem_dealloc(occidx)
!!$    call mem_dealloc(unoccidx)
!!$    call mem_dealloc(occval)
!!$    call mem_dealloc(unoccval)
!!$
!!$
!!$    ! Print out info
!!$    ! **************
!!$    write(DECinfo%output,*)'FOP'
!!$    write(DECinfo%output,'(1X,a)') 'FOP========================================================='
!!$    write(DECinfo%output,'(1X,a,i4)') 'FOP      SITE OPTIMIZATION HAS CONVERGED FOR SITE',MyAtom
!!$    write(DECinfo%output,'(1X,a)') 'FOP---------------------------------------------------------'
!!$    write(DECinfo%output,'(1X,a,i4)')    'FOP Done: Fragment number                  :', MyAtom
!!$    write(DECinfo%output,'(1X,a,i4)')    'FOP Done: Number of orbitals in virt total :', &
!!$         & AtomicFragment%nunoccFA
!!$    write(DECinfo%output,'(1X,a,i4)')    'FOP Done: Number of orbitals in occ total  :', &
!!$         & AtomicFragment%noccFA
!!$    write(DECinfo%output,'(1X,a,f16.10)') 'FOP Done: Lagrangian Fragment energy       :', &
!!$         & AtomicFragment%LagFOP
!!$    write(DECinfo%output,'(1X,a,f16.10)') 'FOP Done: Occupied Fragment energy         :', &
!!$         & AtomicFragment%EoccFOP
!!$    write(DECinfo%output,'(1X,a,f16.10)') 'FOP Done: Virtual Fragment energy          :', &
!!$         & AtomicFragment%EvirtFOP
!!$    write(DECinfo%output,'(1X,a,i4)')     'FOP Done: Number of basis functions        :', &
!!$         & AtomicFragment%number_basis
!!$    write(DECinfo%output,'(1X,a,g14.2)')  'FOP Done: Occupied reduction threshold     :', &
!!$         & AtomicFragment%RejectThr(1)
!!$    write(DECinfo%output,'(1X,a,g14.2)')  'FOP Done: Virtual  reduction threshold     :', &
!!$         & AtomicFragment%RejectThr(2)
!!$    write(DECinfo%output,'(1X,a,i6,a,i6,a,f5.2,a)')  'FOP Done: Occupied reduction removed ', &
!!$         & nocc_exp-FOfragment%noccAOS, ' of ', nocc_exp, ' orbitals ( ', &
!!$         & (nocc_exp-FOfragment%noccAOS)*100.0_realk/nocc_exp, ' %)'
!!$    write(DECinfo%output,'(1X,a,i6,a,i6,a,f5.2,a)')  'FOP Done: Virtual  reduction removed ', &
!!$         & nvirt_exp-FOfragment%nunoccAOS, ' of ', nvirt_exp, ' orbitals ( ', &
!!$         & (nvirt_exp-FOfragment%nunoccAOS)*100.0_realk/nvirt_exp, ' %)'
!!$    write(DECinfo%output,'(1X,a,/)') 'FOP========================================================='
!!$    write(DECinfo%output,*) 'FOP'
!!$    call atomic_fragment_free(FOfragment)
!!$
!!$    if(freebasisinfo) then
!!$       call atomic_fragment_free_basis_info(AtomicFragment)
!!$    end if
!!$
!!$    ! MPI fragment statistics
!!$    AtomicFragment%slavetime  = slavetime
!!$    AtomicFragment%flops_slaves = flops_slaves
!!$
!!$    ! Ensure that energies within fragment structure are set consistently
!!$    call set_energies_ccatom_structure_fragopt(AtomicFragment)
!!$
!!$
!!$  end subroutine optimize_atomic_fragment_FA_experimental



  !> Routine that optimizes an atomic fragment by (1) expanding to include neighbour atoms,
  !> (2) checking that energy change is smaller than FOT, (3) for converged fragment we 
  !> make unitary transformations of occupied and virtual AOS (NOT including EOS) to 
  !> generate fragment-adapted orbitals which describe AOS with as few orbitals as possible.
  !> \date February 2013
  !> \author Ida-Marie Hoeyvik & Kasper Kristensen
  subroutine optimize_atomic_fragment_FA(MyAtom,AtomicFragment,nAtoms, &
       &OccOrbitals,nOcc,UnoccOrbitals,nUnocc,DistanceTable, &
       &MyMolecule,mylsitem,freebasisinfo,t1,t2,t1full)
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
    type(ccatom),intent(inout)        :: AtomicFragment
    !> All occupied orbitals
    type(ccorbital), dimension(nOcc), intent(in)      :: OccOrbitals
    !> All unoccupied orbitals
    type(ccorbital), dimension(nUnocc), intent(in)    :: UnoccOrbitals
    !> Distance table for all atoms in molecule
    real(realk), dimension(nAtoms,nAtoms), intent(in) :: DistanceTable
    !> Full molecule information
    type(fullmolecule), intent(in) :: MyMolecule
    !> Integral information
    type(lsitem), intent(inout)       :: mylsitem
    !> Delete fragment basis information ("expensive box in ccatom type") at exit?
    logical,intent(in) :: freebasisinfo
    !> t1 amplitudes for AOS (not full molecule)
    type(array2),intent(inout),optional :: t1
    !> t2 amplitudes for AOS (not full molecule)
    type(array4),intent(inout),optional :: t2
    !> t1 amplitudes for full molecule to be updated (only used when DECinfo%SinglesPolari is set)
    type(array2),intent(inout),optional :: t1full
    real(realk)                    :: LagEnergyDiff, OccEnergyDiff,VirtEnergyDiff
    real(realk)                    :: LagEnergyOld, OccEnergyOld, VirtEnergyOld
    real(realk)                    :: FOT
    logical, dimension(natoms)     :: Occ_atoms,Virt_atoms,OccOld,VirtOld
    real(realk),dimension(natoms)  :: DistMyAtom,SortedDistMyAtom
    integer,dimension(natoms)      :: DistTrackMyAtom
    integer,dimension(natoms)      :: nocc_per_atom,nunocc_per_atom
    type(ccatom) :: FOfragment
    integer      :: iter,i,ov
    integer      :: Nnew,Nold, max_iter_red,nocc_exp,nvirt_exp
    logical      :: converged,ReductionPossible(2)
    logical :: expansion_converged, lag_converged, occ_converged, virt_converged
    real(realk) :: slavetime, flops_slaves
    real(realk),pointer :: OccMat(:,:), VirtMat(:,:)

    write(DECinfo%output,'(a)')    ' FOP'
    write(DECinfo%output,'(a)')    ' FOP ==============================================='
    write(DECinfo%output,'(a,i4)') ' FOP  Site fragment generator for fragment,',MyAtom
    write(DECinfo%output,'(a)')    ' FOP ==============================================='
    write(DECinfo%output,'(a)')    ' FOP'


    ! Sanity check for singles polarization
    if(DECinfo%SinglesPolari) then
       if(.not. present(t1full)) then
          call lsquit('optimize_atomic_fragment: Full singles polarization is requrested &
               & but t1 argument is not present!',DECinfo%output)
       end if
    end if


    ! ******************************************
    ! **   Initialization of stuff needed..   **
    ! ******************************************
    converged=.false.
    max_iter_red=15   ! allow 15 reduction steps (should be more than enough)
    FOT = DECinfo%FOT
    DistMyAtom= DistanceTable(:,MyAtom)
    call GetSortedList(SortedDistMyAtom,DistTrackMyAtom,DistanceTable,natoms,MyAtom)
    nocc_per_atom=get_number_of_orbitals_per_atom(OccOrbitals,nocc,natoms)
    nunocc_per_atom=get_number_of_orbitals_per_atom(UnoccOrbitals,nunocc,natoms)

    ! MPI fragment statistics
    flops_slaves = 0.0E0_realk
    slavetime = 0.0E0_realk

    if( (nocc_per_atom(MyAtom) == 0) .and. (nunocc_per_atom(MyAtom) == 0) ) then
       write(DECinfo%output,*) 'FOP Skipping optimization of fragment ', MyAtom
       AtomicFragment%LagFOP=0E0_realk
       AtomicFragment%EoccFOP=0E0_realk
       AtomicFragment%EvirtFOP=0E0_realk
       return
    end if

    ! ******************************************
    ! **  Starting computation of fragment    **
    ! ******************************************
    if(DECinfo%simulate_full .or. DECinfo%InclFullMolecule) then
       ! skip calculation of energy and just initiate fragment including ALL orbitals
       write(DECinfo%output,*) 'FOP Fragment includes all orbitals and fragment optimization is skipped'
       call fragment_init_simulate_full(MyAtom,nunocc, nocc, OccOrbitals,UnoccOrbitals,&
            & MyMolecule,mylsitem,AtomicFragment,.true.)

       if(DECinfo%fragadapt) then
          write(DECinfo%output,*) 'WARNING! For fragment-adapted orbitals we cannot skip energy calculation! '
          call single_lagrangian_energy_and_prop(AtomicFragment)
       end if
       call fragment_adapted_transformation_matrices(AtomicFragment)

       call PrintInfo_Lagrangian(AtomicFragment,0.0E0_realk,0.0E0_realk,0.0E0_realk,0)
       if(freebasisinfo) then
          call atomic_fragment_free_basis_info(AtomicFragment)
       end if
       return
    else
       call InitialFragment_L(Occ_atoms,Virt_atoms,DistMyAtom,nAtoms)
       call get_fragment_and_Energy(MyAtom,natoms,Occ_Atoms,Virt_Atoms,&
            & MyMolecule,MyLsitem,nocc,nunocc,OccOrbitals,UnoccOrbitals,&
            & AtomicFragment)
       ! MPI fragment statistics
       slavetime = slavetime +  AtomicFragment%slavetime
       flops_slaves = flops_slaves + AtomicFragment%flops_slaves
    end if

    ! Print initial fragment information
    write(DECinfo%output,*)'FOP'
    write(DECinfo%output,'(1X,a)') 'FOP========================================================='
    write(DECinfo%output,'(1X,a)') 'FOP               Initial fragment information'
    write(DECinfo%output,'(1X,a)') 'FOP---------------------------------------------------------'
    write(DECinfo%output,'(1X,a,i4)')    'FOP Init: Fragment number                  :', MyAtom
    write(DECinfo%output,'(1X,a,i4)')    'FOP Init: Number of orbitals in virt total :', &
         & AtomicFragment%nunoccAOS
    write(DECinfo%output,'(1X,a,i4)')    'FOP Init: Number of orbitals in occ total  :', &
         & AtomicFragment%noccAOS
    write(DECinfo%output,'(1X,a,f16.10)') 'FOP Init: Lagrangian Fragment energy       :', &
         & AtomicFragment%LagFOP
    write(DECinfo%output,'(1X,a,f16.10)') 'FOP Init: Occupied Fragment energy         :', &
         & AtomicFragment%EoccFOP
    write(DECinfo%output,'(1X,a,f16.10)') 'FOP Init: Virtual Fragment energy          :', &
         & AtomicFragment%EvirtFOP
    write(DECinfo%output,'(1X,a,i4)')    'FOP Init: Number of basis functions        :', &
         & AtomicFragment%number_basis
    write(DECinfo%output,'(1X,a,/)') 'FOP========================================================='
    write(DECinfo%output,*) 'FOP'



    ! ======================================================================
    !                             Expansion loop
    ! ======================================================================

    expansion_converged=.false.
    lag_converged=.false.
    occ_converged=.false.
    virt_converged=.false.
    EXPANSION_LOOP: do iter = 1,DECinfo%MaxIter
       OccOld=Occ_atoms;VirtOld=Virt_atoms
       LagEnergyOld = AtomicFragment%LagFOP
       OccEnergyOld = AtomicFragment%EoccFOP
       VirtEnergyOld = AtomicFragment%EvirtFOP
       ! correlation density matrices
       call mem_alloc(OccMat,AtomicFragment%noccAOS,AtomicFragment%noccAOS)
       call mem_alloc(VirtMat,AtomicFragment%nunoccAOS,AtomicFragment%nunoccAOS)
       OccMat = AtomicFragment%OccMat
       VirtMat = AtomicFragment%VirtMat

       call Expandfragment(Occ_atoms,Virt_atoms,DistTrackMyAtom,natoms,&
            & nocc_per_atom,nunocc_per_atom)
       call atomic_fragment_free(AtomicFragment)
       call get_fragment_and_Energy(MyAtom,natoms,Occ_Atoms,Virt_Atoms,&
            & MyMolecule,MyLsitem,nocc,nunocc,OccOrbitals,UnoccOrbitals,&
            & AtomicFragment)
       ! MPI fragment statistics
       slavetime = slavetime +  AtomicFragment%slavetime
       flops_slaves = flops_slaves + AtomicFragment%flops_slaves
       LagEnergyDiff=abs(LagEnergyOld-AtomicFragment%LagFOP)
       OccEnergyDiff=abs(OccEnergyOld-AtomicFragment%EoccFOP)
       VirtEnergyDiff=abs(VirtEnergyOld-AtomicFragment%EvirtFOP)
       call PrintInfo_Lagrangian(AtomicFragment,LagEnergyDiff,OccEnergyDiff,VirtEnergyDiff,iter)

       ! Test convergence for both Lagrangian, occupied, and virtual energies
       ! ********************************************************************

       ! Lagrangian
       TEST_CONVERGENCE_LAG: if  (LagEnergyDiff < FOT) then
          write(DECinfo%output,'(1X,a,F14.9)') 'FOP: Lagrangian energy converged, energydiff =', LagEnergyDiff
          lag_converged=.true.
       else
          write(DECinfo%output,*) 'FOP: Lagrangian energy NOT converged'
          lag_converged=.false.
       end if TEST_CONVERGENCE_LAG

       ! Occupied
       TEST_CONVERGENCE_OCC: if  (OccEnergyDiff < FOT) then
          write(DECinfo%output,'(1X,a,F14.9)') 'FOP: Occupied energy converged, energydiff   =', OccEnergyDiff
          occ_converged=.true.
       else
          write(DECinfo%output,*) 'FOP: Occupied energy NOT converged'
          occ_converged=.false.
       end if TEST_CONVERGENCE_OCC

       ! Virtual
       TEST_CONVERGENCE_VIRT: if  (VirtEnergyDiff < FOT) then
          write(DECinfo%output,'(1X,a,F14.9)') 'FOP: Virtual energy converged, energydiff    =', VirtEnergyDiff
          virt_converged=.true.
       else
          write(DECinfo%output,*) 'FOP: Virtual energy NOT converged'
          virt_converged=.false.
       end if TEST_CONVERGENCE_VIRT

       ! We are converged only if ALL three energies are converged
       ExpansionConvergence: if(lag_converged .and. occ_converged .and. virt_converged) then
          expansion_converged=.true.
          Occ_atoms = OccOld;Virt_atoms = VirtOld
          write(DECinfo%output,*) 'FOP Fragment expansion converged in iteration ', iter
          exit
       end if ExpansionConvergence

       call mem_dealloc(OccMat)
       call mem_dealloc(VirtMat)

    end do EXPANSION_LOOP


    ! Check that expansion loop is converged
    if(.not. expansion_converged) then
       write(DECinfo%output,*) 'Number of expansion steps = ', DECinfo%MaxIter
       call lsquit('Lagrangian fragment expansion did not converge! &
            & Try to increase the number of expansion steps using the .MaxIter keyword',DECinfo%output)
    end if


    ! Set AtomicFragment to be the converged fragment
    call atomic_fragment_free(AtomicFragment)
    call atomic_fragment_init_atom_specific(MyAtom,natoms,Virt_Atoms, &
         & Occ_Atoms,nocc,nunocc,OccOrbitals,UnoccOrbitals, &
         & MyMolecule,mylsitem,AtomicFragment,.true.,.false.)
    ! Set energies correctly without having to do a new calculation
    AtomicFragment%LagFOP = LagEnergyOld
    AtomicFragment%EoccFOP = OccEnergyOld
    AtomicFragment%EvirtFOP = VirtEnergyOld

    ! correlation density matrices
    call mem_alloc(AtomicFragment%OccMat,AtomicFragment%noccAOS,AtomicFragment%noccAOS)
    call mem_alloc(AtomicFragment%VirtMat,AtomicFragment%nunoccAOS,AtomicFragment%nunoccAOS)
    AtomicFragment%CDSet=.true. ! correlation density matrices have been set
    AtomicFragment%OccMat = OccMat
    AtomicFragment%VirtMat = VirtMat

    call mem_dealloc(OccMat)
    call mem_dealloc(VirtMat)

    ! Save dimensions for statistics
    nocc_exp = AtomicFragment%noccAOS
    nvirt_exp = AtomicFragment%nunoccAOS


    write(DECinfo%output,*) ' FOP'
    write(DECinfo%output,*) ' FOP *************************************************'
    write(DECinfo%output,*) ' FOP ** Expansion has converged. We start reduction **'
    write(DECinfo%output,*) ' FOP *************************************************'
    write(DECinfo%output,*) ' FOP'



    ! ======================================================================
    !                             Reduction loop
    ! ======================================================================
    ! Here we carry out unitary transformation of occupied and virtual AOSs 
    ! (NOT including EOS orbitals) to adapt the AOS to the EOS orbitals.
    ! We term the new set of orbitals "fragment-adapted orbitals".
    ReductionPossible=.false.


    ! First reduce virtual space (ov=2) then occupied space (ov=1)
    OCC_OR_VIRT: do ov=2,1,-1   

       lag_converged=.false.
       occ_converged=.false.
       virt_converged=.false.
       converged=.false.

       REDUCTION_LOOP: do iter=1,max_iter_red

          if(ov==1) then
             write(DECinfo%output,*) 'FOP Starting occ reduction step ', iter
          else
             write(DECinfo%output,*) 'FOP Starting virt reduction step ', iter
          end if

          if(iter == 1) then
             ! Number of occ or virt orbitals in converged fragment
             if(ov==1) then
                Nold = AtomicFragment%noccAOS
             else
                Nold = AtomicFragment%nunoccAOS
             end if
             ! Threshold for throwing away fragment-adapted orbitals (start with FOT value)
             AtomicFragment%RejectThr(ov) = FOT
             if(ov==1) call atomic_fragment_free(FOfragment)
          else  
             ! Number of occ or virt fragment-adapted orbitals from previous step
             if(ov==1) then
                Nold = FOfragment%noccAOS
             else
                Nold = FOfragment%nunoccAOS
             end if
             ! decrease rejection threshold by a factor 10 in each step
             AtomicFragment%RejectThr(ov) = AtomicFragment%RejectThr(ov)/10.0_realk
             call atomic_fragment_free(FOfragment)
          end if
          if(ov==1) then
             write(DECinfo%output,'(a,ES13.5)') 'FOP occ rejection threshold  : ',&
                  & AtomicFragment%RejectThr(ov)
          else
             write(DECinfo%output,'(a,ES13.5)') 'FOP virt rejection threshold  : ',&
                  & AtomicFragment%RejectThr(ov)
          end if


          ! Make fragment-adapted fragment with smaller AOS according to rejection threshold
          call fragment_adapted_driver(MyMolecule,mylsitem,OccOrbitals,UnoccOrbitals,&
               & AtomicFragment,FOfragment)
          if(ov==1) then
             Nnew = FOfragment%noccAOS
          else
             Nnew = FOfragment%nunoccAOS
          end if

          ! Cycle if the AOS was not decreased in the rejection procedure
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
          call single_lagrangian_energy_and_prop(FOfragment)


          ! MPI fragment statistics
          slavetime = slavetime +  FOfragment%slavetime
          flops_slaves = flops_slaves + FOfragment%flops_slaves

          LagEnergyDiff=abs(FOfragment%LagFOP-LagEnergyOld)
          OccEnergyDiff=abs(FOfragment%EoccFOP-OccEnergyOld)
          VirtEnergyDiff=abs(FOfragment%EvirtFOP-VirtEnergyOld)
          call PrintInfo_Lagrangian(FOfragment,LagEnergyDiff,OccEnergyDiff,VirtEnergyDiff,iter)


          ! Test convergence for both Lagrangian, occupied, and virtual energies
          ! ********************************************************************

          ! Lagrangian
          TEST_REDUCTION_LAG: if  (LagEnergyDiff < FOT) then
             if(ov==1) then
                write(DECinfo%output,'(1X,a,F14.9)') &
                     & 'FOP: Lagrangian energy converged (occ reduction), diff =', LagEnergyDiff
             else
                write(DECinfo%output,'(1X,a,F14.9)') &
                     & 'FOP: Lagrangian energy converged (virt reduction), diff =', LagEnergyDiff
             end if
             lag_converged=.true.
          else
             if(ov==1) then
                write(DECinfo%output,*) 'FOP: Lagrangian energy (occ reduction) NOT converged'
             else
                write(DECinfo%output,*) 'FOP: Lagrangian energy (virt reduction) NOT converged'
             end if
             lag_converged=.false.
          end if TEST_REDUCTION_LAG

          ! Occupied
          TEST_REDUCTION_OCC: if  (OccEnergyDiff < FOT) then
             if(ov==1) then
                write(DECinfo%output,'(1X,a,F14.9)') &
                     & 'FOP: Occupied energy converged (occ reduction), diff =', OccEnergyDiff
             else
                write(DECinfo%output,'(1X,a,F14.9)') &
                     & 'FOP: Occupied energy converged (virt reduction), diff =', OccEnergyDiff
             end if
             occ_converged=.true.
          else
             if(ov==1) then
                write(DECinfo%output,*) 'FOP: Occupied energy (occ reduction) NOT converged'
             else
                write(DECinfo%output,*) 'FOP: Occupied energy (virt reduction) NOT converged'
             end if
             occ_converged=.false.
          end if TEST_REDUCTION_OCC

          ! Virtual
          TEST_REDUCTION_VIRT: if  (VirtEnergyDiff < FOT) then
             if(ov==1) then
                write(DECinfo%output,'(1X,a,F14.9)') &
                     & 'FOP: Virtual energy converged (occ reduction), diff =', VirtEnergyDiff
             else
                write(DECinfo%output,'(1X,a,F14.9)') &
                     & 'FOP: Virtual energy converged (virt reduction), diff =', VirtEnergyDiff
             end if
             virt_converged=.true.
          else
             if(ov==1) then
                write(DECinfo%output,*) 'FOP: Virtual energy (occ reduction) NOT converged'
             else
                write(DECinfo%output,*) 'FOP: Virtual energy (virt reduction) NOT converged'
             end if
             virt_converged=.false.
          end if TEST_REDUCTION_VIRT
          if (lag_converged .and. occ_converged .and. virt_converged) then
             converged=.true.
             ReductionPossible(ov)=.true.
          end if



          ! Quit if fragment reduction is converged
          ! ***************************************
          if(converged) then
             if(ov==2) then  ! virtual reduction converged 
                write(DECinfo%output,*) 'FOP Virt reduction of fragment converged in step', iter
                AtomicFragment%LagFOP = FOfragment%LagFOP
                AtomicFragment%EoccFOP = FOfragment%EoccFOP
                AtomicFragment%EvirtFOP = FOfragment%EvirtFOP
                cycle OCC_OR_VIRT  ! try to reduce occ space

             else ! both occ and virtual reductions have converged
                write(DECinfo%output,*) 'FOP Occ  reduction of fragment converged in step', iter
                AtomicFragment%LagFOP = FOfragment%LagFOP
                AtomicFragment%EoccFOP = FOfragment%EoccFOP
                AtomicFragment%EvirtFOP = FOfragment%EvirtFOP
                exit
             end if
          end if

1000      continue

          if ( (iter == max_iter_red) .and. (ov==1) .and. &
               & (.not. any(ReductionPossible) ) ) then
             write(DECinfo%output,*) "FOP No reduction possible. Use original converged fragment"

             ! Go back to old fragment
             ! ***********************
             ! This is extra work but it will never be done in practice, only
             ! relevant for small debug systems.
             call atomic_fragment_free(AtomicFragment)
             call get_fragment_and_Energy(MyAtom,natoms,Occ_atoms,Virt_atoms,&
                  & MyMolecule,MyLsitem,nocc,nunocc,OccOrbitals,UnoccOrbitals,&
                  & AtomicFragment)
             call atomic_fragment_free(FOfragment)
             call fragment_adapted_driver(MyMolecule,mylsitem,OccOrbitals,UnoccOrbitals,&
                  & AtomicFragment,FOfragment)
             ! MPI fragment statistics
             slavetime = slavetime +  AtomicFragment%slavetime
             flops_slaves = flops_slaves + AtomicFragment%flops_slaves
          end if

       end do REDUCTION_LOOP

    end do OCC_OR_VIRT

    ! Print out info
    ! **************
    write(DECinfo%output,*)'FOP'
    write(DECinfo%output,'(1X,a)') 'FOP========================================================='
    write(DECinfo%output,'(1X,a,i4)') 'FOP      SITE OPTIMIZATION HAS CONVERGED FOR SITE',MyAtom
    write(DECinfo%output,'(1X,a)') 'FOP---------------------------------------------------------'
    write(DECinfo%output,'(1X,a,i4)')    'FOP Done: Fragment number                  :', MyAtom
    write(DECinfo%output,'(1X,a,i4)')    'FOP Done: Number of orbitals in virt total :', &
         & AtomicFragment%nunoccFA
    write(DECinfo%output,'(1X,a,i4)')    'FOP Done: Number of orbitals in occ total  :', &
         & AtomicFragment%noccFA
    write(DECinfo%output,'(1X,a,f16.10)') 'FOP Done: Lagrangian Fragment energy       :', &
         & AtomicFragment%LagFOP
    write(DECinfo%output,'(1X,a,f16.10)') 'FOP Done: Occupied Fragment energy         :', &
         & AtomicFragment%EoccFOP
    write(DECinfo%output,'(1X,a,f16.10)') 'FOP Done: Virtual Fragment energy          :', &
         & AtomicFragment%EvirtFOP
    write(DECinfo%output,'(1X,a,i4)')     'FOP Done: Number of basis functions        :', &
         & AtomicFragment%number_basis
    write(DECinfo%output,'(1X,a,g14.2)')  'FOP Done: Occupied reduction threshold     :', &
         & AtomicFragment%RejectThr(1)
    write(DECinfo%output,'(1X,a,g14.2)')  'FOP Done: Virtual  reduction threshold     :', &
         & AtomicFragment%RejectThr(2)
    write(DECinfo%output,'(1X,a,i6,a,i6,a,f5.2,a)')  'FOP Done: Occupied reduction removed ', &
         & nocc_exp-FOfragment%noccAOS, ' of ', nocc_exp, ' orbitals ( ', &
         & (nocc_exp-FOfragment%noccAOS)*100.0_realk/nocc_exp, ' %)'
    write(DECinfo%output,'(1X,a,i6,a,i6,a,f5.2,a)')  'FOP Done: Virtual  reduction removed ', &
         & nvirt_exp-FOfragment%nunoccAOS, ' of ', nvirt_exp, ' orbitals ( ', &
         & (nvirt_exp-FOfragment%nunoccAOS)*100.0_realk/nvirt_exp, ' %)'
    write(DECinfo%output,'(1X,a,/)') 'FOP========================================================='
    write(DECinfo%output,*) 'FOP'
    call atomic_fragment_free(FOfragment)

    if(freebasisinfo) then
       call atomic_fragment_free_basis_info(AtomicFragment)
    end if

    ! MPI fragment statistics
    AtomicFragment%slavetime  = slavetime
    AtomicFragment%flops_slaves = flops_slaves

    ! Ensure that energies in fragment are set consistently
    call set_energies_ccatom_structure_fragopt(AtomicFragment)

  end subroutine optimize_atomic_fragment_FA



  !> \brief prints info for atomic fragment in a given iteration of the fragment optimization
  !> using the Lagrangian partitioning scheme.
  !> \date: august-2011
  !> \author: Ida-Marie Hoeyvik
  subroutine PrintInfo_Lagrangian(Fragment,LagEnergyDiff,OccEnergyDiff,VirtEnergyDiff,iter)
    implicit none
    !> Atomic fragment
    type(ccatom),intent(inout) :: Fragment
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
   write(DECinfo%output,'(1X,a,i4)')    'FOP Loop: Fragment number                  :', fragment%atomic_number
   write(DECinfo%output,'(1X,a,i4)')    'FOP Loop: Number of orbitals in virt total :', Fragment%nunoccAOS
   write(DECinfo%output,'(1X,a,i4)')    'FOP Loop: Number of orbitals in occ total  :', Fragment%noccAOS
   write(DECinfo%output,'(1X,a,f16.10)') 'FOP Loop: Lagrangian Fragment energy       :', Fragment%LagFOP
   write(DECinfo%output,'(1X,a,f16.10)') 'FOP Loop: Occupied Fragment energy         :', Fragment%EoccFOP
   write(DECinfo%output,'(1X,a,f16.10)') 'FOP Loop: Virtual Fragment energy          :', Fragment%EvirtFOP
   write(DECinfo%output,'(1X,a,f16.10)') 'FOP Loop: Lagrangian energy diff           :', LagEnergyDiff
   write(DECinfo%output,'(1X,a,f16.10)') 'FOP Loop: Occupied energy diff             :', OccEnergyDiff
   write(DECinfo%output,'(1X,a,f16.10)') 'FOP Loop: Virtual energy diff              :', VirtEnergyDiff
   write(DECinfo%output,'(1X,a,i4)')    'FOP Loop: Number of basis functions        :', Fragment%number_basis
   write(DECinfo%output,'(1X,a,/)') 'FOP========================================================='
   write(DECinfo%output,*) 'FOP'


  end subroutine PrintInfo_Lagrangian



  !> Expand both occ and virt site fragment
  !> date: august-2011
  !> author: Ida-Marie Hoyvik
  subroutine  Expandfragment(Occ,Virt,Track,natoms,&
       & nocc_per_atom,nunocc_per_atom)
    implicit none
    integer,intent(in)        :: natoms
    logical,dimension(natoms) :: Occ,Virt,Temp
    integer,intent(in)        :: nocc_per_atom(natoms)
    integer,intent(in)        :: nunocc_per_atom(natoms)
    integer,intent(in)        :: Track(natoms)
    integer                   :: i,indx,counter,step

    step = DECinfo%LagStepSize

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


    !Expand virt space
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

  end subroutine Expandfragment


  !> \brief Reduce occupied or virtual AOS for fragment using the individual orbital
  !> contributions according to the input threshold.
  !> \author Kasper Kristensen
  !> \date December 2011
  subroutine ReduceSpace_orbitalspecific(MyFragment,norb_full,contributions,OccOrVirt,&
                            &Thresh,OrbAOS,Nafter)

    implicit none
    !> Fragment info
    type(ccatom),intent(inout) :: MyFragment
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

    ! Number of orbitals at input
    Nbefore = count(OrbAOS)


    ! Exclude orbitals with contributions smaller than threshold
    do i=1,norb_full
       if (abs(contributions(i)) < Thresh) OrbAOS(i)=.false.
    end do

    ! Sanity check: The orbitals assigned to the central atom in the fragment should ALWAYS be included
    if(OccOrVirt=='O') then ! checking occupied orbitals
       do i=1,MyFragment%noccEOS
          OrbAOS(MyFragment%occEOSidx(i)) = .true.
       end do
    elseif(OccOrVirt=='V') then
       do i=1,MyFragment%nunoccEOS
          OrbAOS(MyFragment%unoccEOSidx(i)) = .true.
       end do
    else
       call lsquit('ReduceSpace_orbitalspecific: OccOrVirt input must be O or V',DECinfo%output)
    end if

    Nafter = count(OrbAOS)
    Nexcl = Nbefore-Nafter

    write(DECinfo%output,'(a,i4)') ' FOP Number of orbitals excluded: ', Nexcl

  end subroutine ReduceSpace_orbitalspecific



  !> \brief For a given model, get the occupied, virtual and Lagragian fragment energies
  !> to use for fragment optimization, i.e. simply copy the relevant energies from
  !> "fragment%energies" to fragment%EoccFOP, fragment%EvirtFOP, and fragment%LagFOP.
  !> (See description of EoccFOP, EvirtFOP, and LagFOP in ccatom type definition).
  !> \author Kasper Kristensen
  !> \date March 2013
  subroutine get_occ_virt_lag_energies_fragopt(Fragment)
    implicit none
    type(ccatom),intent(inout) :: fragment
    ! MODIFY FOR NEW MODEL
    ! If you implement a new model, please copy fragment%energies(?) for your model,
    ! see ccatom type def to determine the "?".

    fragment%EoccFOP = 0.0_realk
    fragment%EvirtFOP = 0.0_realk
    fragment%LagFOP = 0.0_realk

    select case(DECinfo%ccmodel)
    case(1)
       ! MP2
       fragment%LagFOP = fragment%energies(1)
       fragment%EoccFOP = fragment%energies(2)
       fragment%EvirtFOP = fragment%energies(3)
    case(2)
       ! CC2
       fragment%EoccFOP = fragment%energies(4)
       fragment%EvirtFOP = fragment%energies(5)
       ! simply use average of occ and virt energies since Lagrangian is not yet implemented
       fragment%LagFOP =  0.5_realk*(fragment%EoccFOP+fragment%EvirtFOP)   
    case(3)
       ! CCSD
       fragment%EoccFOP = fragment%energies(6)
       fragment%EvirtFOP = fragment%energies(7)
       ! simply use average of occ and virt energies since Lagrangian is not yet implemented
       fragment%LagFOP =  0.5_realk*(fragment%EoccFOP+fragment%EvirtFOP)
    case(4)
       ! CCSD(T)
       fragment%EoccFOP = fragment%energies(8)
       fragment%EvirtFOP = fragment%energies(9)
       ! simply use average of occ and virt energies since Lagrangian is not yet implemented
       fragment%LagFOP =  0.5_realk*(fragment%EoccFOP+fragment%EvirtFOP)
    case default
       write(DECinfo%output,*) 'WARNING: get_occ_virt_lag_energies_fragopt needs implementation &
            & for model:', DECinfo%ccmodel
    end select

  end subroutine get_occ_virt_lag_energies_fragopt




  !> \brief After the fragment optimization, make sure that the energies stored
  !> in fragment%LagFOP, fragment%EoccFOP, and fragment%EvirtFOP are copied correctly 
  !> to the general fragment%energies arrays.
  !> This is effectively the inverse routine of get_occ_virt_lag_energies_fragopt.
  !> (See description of energies, EoccFOP, EvirtFOP, and LagFOP in ccatom type definition).
  !> \author Kasper Kristensen
  !> \date March 2013
  subroutine set_energies_ccatom_structure_fragopt(Fragment)
    implicit none
    type(ccatom),intent(inout) :: fragment
    ! MODIFY FOR NEW MODEL 
    ! If you implement a new model, please set fragment%energies(?) for your model,
    ! see ccatom type def to determine the "?".

    select case(DECinfo%ccmodel)
    case(1)
       ! MP2
       fragment%energies(1) = fragment%LagFOP 
       fragment%energies(2) = fragment%EoccFOP
       fragment%energies(3) = fragment%EvirtFOP 
    case(2)
       ! CC2
       fragment%energies(4) = fragment%EoccFOP
       fragment%energies(5) = fragment%EvirtFOP
    case(3)
       ! CCSD
       fragment%energies(6) = fragment%EoccFOP 
       fragment%energies(7) = fragment%EvirtFOP 
    case(4)
       ! CCSD(T)
       fragment%energies(8) = fragment%EoccFOP 
       fragment%energies(9) = fragment%EvirtFOP
    case default
       write(DECinfo%output,*) 'WARNING: get_occ_virt_lag_energies_fragopt needs implementation &
            & for model:', DECinfo%ccmodel
    end select

  end subroutine set_energies_ccatom_structure_fragopt



  !> \brief Given the list of all fragment energies as defined in main_fragment_driver,
  !> extract the fragment energies for Lagrangian, occupied, and virtual part. schemes
  !> for the given CC model (see "energies" in ccatom typedef).
  !> \author Kasper Kristensen
  !> \date March 2013
  subroutine extract_fragenergies_for_model(natoms,FragEnergiesAll,FragEnergiesModel)
    
    implicit none

    !> Number of atoms in molecule
    integer,intent(in) :: natoms
    !> Fragment energies for all models as determined in main_fragment_driver.
    real(realk),dimension(natoms,natoms,ndecenergies),intent(in) :: FragEnergiesAll
    !> Fragment energies for model under consideration (DECinfo%ccmodel) in the order:
    !> Lagrangian (:,:,1)   Occupied(:,:,2)    Virtual(:,:,3)
    real(realk),dimension(natoms,natoms,3),intent(inout) :: FragEnergiesModel
    integer :: i,j
    ! MODIFY FOR NEW MODEL
    ! If you implement a new model, please introduce your model below,
    ! see "energies" in ccatom type def to determine how to copy energies for your model.


    FragEnergiesModel=0.0_realk

    do i=1,natoms
       do j=1,natoms
          select case(DECinfo%ccmodel)
          case(1)
             ! MP2

             ! Lagrangian MP2 energy stored in entry 1 (see "energies" in ccatom type def)
             FragEnergiesModel(i,j,1) = FragEnergiesAll(i,j,1)

             ! Occupied MP2 energies stored in entry 2
             FragEnergiesModel(i,j,2) = FragEnergiesAll(i,j,2)

             ! Virtual MP2 energies stored in entry 3
             FragEnergiesModel(i,j,3) = FragEnergiesAll(i,j,3)
          case(2)
             ! CC2

             ! Occupied CC2 energies stored in entry 4
             FragEnergiesModel(i,j,2) = FragEnergiesAll(i,j,4)

             ! Virtual CC2 energies stored in entry 5
             FragEnergiesModel(i,j,3) = FragEnergiesAll(i,j,5)

             ! Lagrangian CC2 energy not implemented, simply use average of occ and virt energies
             FragEnergiesModel(i,j,1) = 0.5_realk*(FragEnergiesModel(i,j,2) + &
                  & FragEnergiesModel(i,j,3) )
          case(3)
             ! CCSD

             ! Occupied CCSD energies stored in entry 6
             FragEnergiesModel(i,j,2) = FragEnergiesAll(i,j,6)

             ! Virtual CCSD energies stored in entry 7
             FragEnergiesModel(i,j,3) = FragEnergiesAll(i,j,7)

             ! Lagrangian CCSD energy not implemented, simply use average of occ and virt energies
             FragEnergiesModel(i,j,1) = 0.5_realk*(FragEnergiesModel(i,j,2) + &
                  & FragEnergiesModel(i,j,3) )

          case(4)
             ! CCSD(T)

             ! Occupied CCSD energies stored in entry 6 + occupied (T) energies stored in entry 8 
             FragEnergiesModel(i,j,2) = FragEnergiesAll(i,j,6) + FragEnergiesAll(i,j,8)

             ! Virtual CCSD energies stored in entry 7 + virtual (T) energies stored in entry 9
             FragEnergiesModel(i,j,2) = FragEnergiesAll(i,j,7) + FragEnergiesAll(i,j,9)

             ! Lagrangian CCSD(T) energy not implemented, simply use average of occ and virt energies
             FragEnergiesModel(i,j,1) = 0.5_realk*(FragEnergiesModel(i,j,2) + &
                  & FragEnergiesModel(i,j,3) )

          case default
             write(DECinfo%output,*) 'WARNING: extract_fragenergies_for_model: Needs implementation &
                  & for model:', DECinfo%ccmodel
          end select

       end do
    end do


  end subroutine extract_fragenergies_for_model

end module fragment_energy_module
