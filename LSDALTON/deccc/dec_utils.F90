!> @file
!> Utils for DEC subroutines
!> \author Marcin Ziolkowski (modified by Kasper Kristensen)
module dec_fragment_utils
  use,intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
  use fundamental
  use precision
  use lstiming
  use ls_util!,only: dgemm_ts
  use typedeftype!,only: lsitem
  use molecule_module!, only: get_geometry
  use files!,only:lsopen,lsclose
  use DALTONINFO!, only: ls_free
  use dec_typedef_module
  use dec_workarounds_module
  use background_buffer_module
  use memory_handling!, only: mem_alloc, mem_dealloc, mem_allocated_global,&
  !       & stats_mem, get_avaiLable_memory
  use matrix_module!, only:matrix
  use matrix_operations
  use tensor_interface_module
  use array4_simple_operations
  use IntegralInterfaceMOD!, only: ii_get_h1, ii_get_nucpot
  use BUILDAOBATCH
  use II_XC_interfaceModule,only:II_get_xc_energy
#ifdef VAR_MPI
  use infpar_module
  use lsmpi_op
#endif
  use papi_module
  use gpu_interfaces  

  ! F12 DEPENDENCIES 
  ! *****************************************
  !use f12_routines_module
  
  ! DEC DEPENDENCIES (within deccc directory)
  ! *****************************************

  !> Read 64 bit integer(s) from file and convert to 32 bit
  interface read_64bit_to_int
     module procedure read_64bit_to_32bit_singleinteger
     module procedure read_64bit_to_32bit_vectorinteger
     module procedure read_64bit_to_32bit_singlelogical
     module procedure read_64bit_to_32bit_vectorlogical
  end interface read_64bit_to_int
  interface read_64bit_to_32bit
     module procedure read_64bit_to_32bit_singleinteger
     module procedure read_64bit_to_32bit_vectorinteger
     module procedure read_64bit_to_32bit_singlelogical
     module procedure read_64bit_to_32bit_vectorlogical
  end interface read_64bit_to_32bit
  interface read_32bit_to_64bit
     module procedure read_32bit_to_64bit_singleinteger
     module procedure read_32bit_to_64bit_vectorinteger
     module procedure read_32bit_to_64bit_singlelogical
     module procedure read_32bit_to_64bit_vectorlogical
  end interface read_32bit_to_64bit
  interface get_combined_SingleDouble_amplitudes
     module procedure get_combined_SingleDouble_amplitudes_newarr
     module procedure get_combined_SingleDouble_amplitudes_oldarr
  end interface get_combined_SingleDouble_amplitudes
  interface remove_core_orbitals_from_last_index
     module procedure remove_core_orbitals_from_last_index_oldarr
     module procedure remove_core_orbitals_from_last_index_newarr
  end interface remove_core_orbitals_from_last_index

private

public :: dec_time_evaluate_efficiency_frag, &
     & dec_fragment_time_init, dec_fragment_time_get, &
     & remove_repeted_entries, ExpandBufferKraken, &
     & TempExpBuffer, RejectAtoms, &
     & StepwiseInclusionOfAtoms, &
     & AtomsToIncludeStepwise, &
     & ListOcc, CountNonZeroElements, &
     & ExpandTarget, count_atoms, &
     & sort_track, sort_track_vector, &
     & adjust_basis_matrix, adjust_basis_matrix2, &
     & adjust_square_matrix, adjust_square_matrix2, &
     & adjust_square_matrix_mo, GetSubSystemIndex, &
     & GetDistances, solve_linear_equations, &
     & invert_matrix, get_power_of_symmetric_matrix, &
     & ExcludeIfNoOrbs, FindMaxDistance, &
     & ReduceBuffer, absorb_logical_vector, &
     & get_logical_pair_vector, InitialFragment, &
     & dec_read_mat_from_file, dec_simple_dgemm, &
     & dec_simple_dgemm_update, calculate_fragment_memory, &
     & mypointer_init, start_flop_counter, &
     & end_flop_counter, fragment_print, &
     & atomic_fragment_free, atomic_fragment_free_simple, &
     & atomic_fragment_free_basis_info, fragmentAOS_type_free, &
     & atomic_fragment_free_f12, free_fragment_t1, orbital_free, &
     & get_memory_for_dec_calculation, get_currently_available_memory, &
     & dec_regression, dec_regression_get_powers, &
     & simple_ascii_plot, init_SPgridbox, &
     & free_SPgridbox, get_density_from_occ_orbitals, &
     & read_64bit_to_32bit,read_32bit_to_64bit, &
     & read_64bit_to_int,&
     & save_fragment_t1_AOSAOSamplitudes, &
     & which_pairs, which_pairs_occ, which_pairs_virt, &
     & which_pairs_occ_virt, write_fragment_joblist_to_file, &
     & read_fragment_joblist_from_file, init_joblist, &
     & free_joblist, put_job_into_joblist, &
     & estimate_memory_for_mp2_energy, &
     & dec_simple_basis_transform1, &
     & dec_simple_basis_transform2, &
     & dec_diff_basis_transform1, &
     & dec_diff_basis_transform2, &
     & add_dec_energies, estimate_energy_of_skipped_pairs, &
     & project_onto_MO_space, project_out_MO_space, &
     & orthogonalize_MOs, print_total_energy_summary, &
     & print_total_energy_summary_lupri, print_all_fragment_energies, &
     & print_fragment_energies_full, print_atomic_fragment_energies, &
     & print_pair_fragment_energies, print_spec_pair_fragment_energies, &
     & get_occfragenergies, get_virtfragenergies, &
     & get_estimated_energy_error, fragment_basis_point_to_FOs, &
     & fragment_basis_point_to_LOs, fragment_init_dimension_pointers, &
     & get_combined_SingleDouble_amplitudes, &
     & remove_core_orbitals_from_last_index, &
     & secondary_assigning, general_distance_table, &
     & GetOrbAtomDistances, basic_write_jobs_and_fragment_energies_for_restart,&
     & basic_read_jobs_and_fragment_energies_for_restart, &
     & restart_sanity_check, unique_entries, get_distance_between_two_points,&
     & get_matrix_position, count_number_of_nonhydrogen_atoms,&
     & get_distance_between_fragments, get_minimum_distance_between_groups_of_atoms,&
     & get_HF_energy, get_HF_energy_fullmolecule, get_dft_energy_fullmolecule,&
     & fragment_restart_file_exist, get_total_number_of_fragments,&
     & get_fragenergy_restart_filename, get_fragenergy_restart_filename_backup,&
     & get_num_of_pair_fragments, SD_dotproduct, max_batch_dimension

contains
!> \brief Get maximum batch dimension encountered in integral program.
!> \author Kasper Kristensen
!> \date February 2011
function max_batch_dimension(mylsitem,nbasis) result(maxdim)  
  implicit none 
  !> LS item info
  type(lsitem), intent(inout) :: mylsitem
  !> Number of basis function
  integer, intent(in) :: nbasis
  integer :: maxdim
  integer, pointer :: orb2batch(:), batchdim(:)
  integer :: i, nbatches
  
  ! Initialize stuff
  nullify(orb2batch)
  nullify(batchdim)
  call mem_alloc(orb2batch,nbasis)
  
  ! Get batch info
  call II_getBatchOrbitalInfo(mylsitem%setting,nbasis,&
       & orb2Batch,nbatches,DECinfo%output,DECinfo%output)
  
  ! Vector containing dimensions for each batch
  call mem_alloc(batchdim,nbatches)
  batchdim = 0
  do i=1,nbasis
     batchdim(orb2batch(i)) = batchdim(orb2batch(i))+1
  end do
  
  ! Find maximum batch dimension
  maxdim=0
  do i=1,nbatches
     if( batchdim(i) > maxdim ) maxdim=batchdim(i)
  end do
  
  ! Clean up
  call mem_dealloc(batchdim)
  call mem_dealloc(orb2batch)
  
end function max_batch_dimension

  subroutine dec_time_evaluate_efficiency_frag(frag,t,ccmodel,label)
     implicit none
     type(decfrag), intent(inout) :: frag
     real(realk), pointer, intent(inout) :: t(:)
     integer, intent(in) :: ccmodel
     character*(*), intent(in) :: label
     integer(kind=ls_mpik) :: nnod,nod
     real(realk) :: tottime_work, tottime_idle, tottime_comm
     real(realk) :: max_work, max_comm, max_idle, time_tts, time_tot
     real(realk),pointer :: time_tot_node(:)
     nnod = 1
#ifdef VAR_MPI
     nnod = infpar%lg_nodtot
#endif
     call mem_alloc(time_tot_node,nnod)

     tottime_work = 0.0E0_realk
     tottime_comm = 0.0E0_realk
     tottime_idle = 0.0E0_realk
     time_tts     = t(PHASE_WORK_IDX) + t(PHASE_COMM_IDX) + t(PHASE_IDLE_IDX)

     !skip local master time, this time will be added outside
     do nod = 1, nnod
        tottime_work = tottime_work + t((nod-1)*nphases + PHASE_WORK_IDX)
        tottime_comm = tottime_comm + t((nod-1)*nphases + PHASE_COMM_IDX)
        tottime_idle = tottime_idle + t((nod-1)*nphases + PHASE_IDLE_IDX)
        time_tot_node(nod) = &
           &t((nod-1)*nphases + PHASE_WORK_IDX) + &
           &t((nod-1)*nphases + PHASE_COMM_IDX) + &
           &t((nod-1)*nphases + PHASE_IDLE_IDX)
     enddo

     time_tot = tottime_work + tottime_comm + tottime_idle

     !timings are not precise of course, test if they are correct within 1%
     if( abs(time_tot - time_tts*nnod) > time_tot*1.0E-2_realk )then
        print *,"WARNING(dec_time_evaluate_efficiency_frag)&
           &MAYBE SOMETHING WRONG WITH TIMERS",time_tot,time_tts,nnod
        print *,"TOTAL TIMES ON NODES",time_tot_node
     endif

     frag%slavetime_work(ccmodel) = tottime_work
     frag%slavetime_comm(ccmodel) = tottime_comm
     frag%slavetime_idle(ccmodel) = tottime_idle

     if(DECinfo%PL>0.and.DECinfo%print_small_calc)then
        if(time_tot>0.1E-3)then
           write(DECinfo%output,'("Portion time spent working       in ",a," is: ",g10.3,"%")')label,tottime_work/time_tot*100
           write(DECinfo%output,'("Portion time spent communicating in ",a," is: ",g10.3,"%")')label,tottime_comm/time_tot*100
           write(DECinfo%output,'("Portion time spent idle          in ",a," is: ",g10.3,"%")')label,tottime_idle/time_tot*100
        endif
     endif

     call mem_dealloc(time_tot_node)
     call mem_dealloc(t)
  end subroutine dec_time_evaluate_efficiency_frag


  subroutine dec_fragment_time_init(t)
     implicit none
     real(realk), pointer, intent(inout) :: t(:)

#ifdef VAR_MPI
     call init_slave_timers(t,infpar%lg_comm)
#else
     call mem_alloc(t,nphases)
     call time_start_phase(PHASE_WORK, &
        &swinit = t(PHASE_INIT_IDX) ,&
        &swwork = t(PHASE_WORK_IDX) ,&
        &swcomm = t(PHASE_COMM_IDX) ,&
        &swidle = t(PHASE_IDLE_IDX) )
#endif

  end subroutine dec_fragment_time_init

  subroutine dec_fragment_time_get(t)
     implicit none
     real(realk), pointer, intent(inout) :: t(:)
#ifdef VAR_MPI
     call get_slave_timers(t,infpar%lg_comm)
#else
     call time_start_phase(PHASE_WORK, &
        &dwinit = t(PHASE_INIT_IDX) ,&
        &dwwork = t(PHASE_WORK_IDX) ,&
        &dwcomm = t(PHASE_COMM_IDX) ,&
        &dwidle = t(PHASE_IDLE_IDX) )
#endif
  end subroutine dec_fragment_time_get

  !> \brief Returns number of unique elements in a vector
  function unique_entries(vector,size_of_vector) result(num)

    implicit none
    integer, intent(in) :: size_of_vector
    integer, dimension(size_of_vector), intent(in) :: vector
    integer :: num,i

    num=size_of_vector
    do i=2,size_of_vector
       if(vector(i)==vector(i-1)) num=num-1
    end do

    return
  end function unique_entries

  !> \brief Remove repeted entries from a vector
  !> Remove repeted entries from a vector input(input_size) and returns
  !> output(output_size) where output_size is number of unique elements
  subroutine remove_repeted_entries(input_size,input,output_size,output)

    implicit none
    integer, intent(in) :: input_size,output_size
    integer, dimension(input_size), intent(in) :: input
    integer, dimension(output_size), intent(inout) :: output
    integer :: i,unique

    if(input_size <= 1) then
       return
    else
       unique=1
       output(unique)=input(unique)
       do i=2,input_size
          if(input(i)/=input(i-1)) then
             unique=unique+1
             output(unique)=input(i)
          end if
       end do
    end if

    if(unique/=output_size) print *,'something went wrong'

    return
  end subroutine remove_repeted_entries


  subroutine ExpandBufferKraken(BufferVec,TrackMatrix,&
       & MyAtom,natoms,zetacount)
    implicit none
    integer, intent(in) :: natoms,MyAtom,zetacount
    !> logical vector. T if atom is included in buffer
    logical,intent(inout) :: BufferVec(natoms)
    logical :: TempBuffer(natoms)
    !> Vector with sorted list of atoms
    integer,intent(in) :: TrackMatrix(zetacount)
    !> STEP determines how many atoms we include each expansion
    integer :: i,counter,STEP,atoms_buffer,atoms_buffer_new


    call count_atoms(BufferVec,atoms_buffer,natoms)

    counter = 0
    STEP = 4

    !> Go from natoms--> since zeta matrix was arranged from low-->high
    do i=zetacount,1,-1
       if (.not. BufferVec(TrackMatrix(i))) then
          BufferVec(TrackMatrix(i)) = .true.
          counter = counter+1
       end if
       if (counter == STEP) exit
    end do


  end subroutine ExpandBufferKraken


  subroutine TempExpBuffer(BufferVec,Track,&
       & MyAtom,natoms)
    implicit none
    integer, intent(in) :: natoms,MyAtom
    !> logical vector. T if atom is included in buffer
    logical,intent(inout) :: BufferVec(natoms)
    logical :: TempBuffer(natoms)
    !> Vector with sorted list of atoms
    integer,intent(in) :: Track(natoms)
    !> STEP determines how many atoms we include each expansion
    integer :: i,counter,STEP,atoms_buffer,atoms_buffer_new


    call count_atoms(BufferVec,atoms_buffer,natoms)

    counter = 0
    STEP = 4

    !> Go from natoms--> since zeta matrix was arranged from low-->high
    do i=natoms,1,-1
       if (.not. BufferVec(Track(i))) then
          BufferVec(Track(i)) = .true.
          counter = counter+1
       end if
       if (counter == STEP) exit
    end do


  end subroutine TempExpBuffer

  !>Subroutine that takes a vector of energy contributions from all atoms
  !> set to "true" in UnoccEOS_atoms, and sets value to "false" if
  !> energy contribution  for that atom is less than certain thresh
  subroutine RejectAtoms(EOSvector, EnergyContributions, ListOfAtoms,&
       &RejectThresh, EOSToStore,DimOfVect,NafterExcl,Nexcl, natoms)
    implicit none
    integer :: DimOfVect, NAfterExcl, Nexcl, natoms,NEOS
    logical :: EOSvector(natoms), EOSToStore(natoms)
    real(realk), intent(in) :: EnergyContributions(DimOfVect)
    real(realk), intent(in) :: RejectThresh
    integer, intent(in) :: ListOfAtoms(DimOfVect)
    integer :: i

    EOSToStore = EOSvector

    do i=1, DimOfVect
       if (abs(EnergyContributions(i))<RejectThresh) then
          EOSToStore(ListOfAtoms(i))= .false.
       end if
    end do

    call count_atoms(EOSvector,NEOS,natoms)
    call count_atoms(EOSToStore,NafterExcl,natoms)
    NExcl = NEOS-NAfterExcl



  end subroutine RejectAtoms


  !> Subroutine that minimizes size of EOS after convergence have been reached
  !> author: Ida-Marie Hoeyvik
  subroutine StepwiseInclusionOfAtoms(SortedList, SortedContributions,Nexcl,&
       &l,EOS_excl,i,natoms)
    implicit none
    integer, intent(in)        :: Nexcl,natoms
    real(realk), intent(inout) :: SortedContributions(Nexcl)
    integer, intent(inout)     :: SortedList(Nexcl)
    integer, intent(inout)     :: l
    logical, intent(inout)     :: EOS_excl(natoms)
    integer, intent(in)        :: i

    if (l < (Nexcl-1)) then
       !Include two new atoms in the fragment
       l = l+1
       EOS_excl(SortedList(l)) = .true.
       l = l+1
       EOS_excl(SortedList(l)) = .true.
    elseif (l==(Nexcl-1)) then
       l= l+1
       EOS_excl(SortedList(l))= .true.
    else
       write(DECinfo%output,*)'All atoms excluded have been included again.&
            & Since fragment is not converged, something is wrong either with the inclusion&
            & or the energy of the fragment referenced!'
    end if

  end subroutine StepwiseInclusionOfAtoms



  !>Subroutine find atoms that should be included stepwise in StepwiseInclusion
  !> author: Ida-Marie Hoeyvik
  subroutine AtomsToIncludeStepwise(UnoccEOS_atoms, UnoccEOS_excl,list_of_atoms,&
       &energy_contributions,ListOfExclAtoms,ContrForExclAtoms,NEOS, Nexcl, NafterExcl,natoms)
    implicit none
    integer, intent(in)     :: NAfterExcl, NEOS, Nexcl,natoms
    logical, intent(in)     :: UnoccEOS_atoms(natoms), UnoccEOS_excl(natoms)
    real(realk), intent(in) :: energy_contributions(NEOS)
    integer, intent(in)     :: list_of_atoms(NEOS)
    real(realk),intent(out) :: ContrForExclAtoms(Nexcl)
    integer, intent(out)    :: ListOfExclAtoms(Nexcl)
    integer                 :: l,i,k

    l=0
    do i=1,natoms
       if (UnoccEOS_atoms(i) .and. (.not. UnoccEOS_excl(i)) ) then
          l=l+1
          ListOfExclAtoms(l) = i
          do k=1,NEOS
             if (list_of_atoms(k)==i) then
                ContrForExclAtoms(l)=energy_contributions(k)
             end if
          end do
       end if
    end do

  end subroutine AtomsToIncludeStepwise

  !Make list over atoms that does not have occupied orbitals assigned
  !Logical value in List is .true. if atom has occ orbitals and .false. if it has not
  !Based on values of Fik/Fkk, which are zero when no occ orb exist on atom
  subroutine ListOcc(MyAtom,natoms,List,vec)
    implicit none
    integer, intent(in) :: natoms,MyAtom
    real (realk), intent(in) :: vec(natoms)
    logical :: List(natoms)
    real(realk) :: thresh
    integer :: i

    thresh = 1E-20_realk
    do i =1, natoms
       if (vec(i)< thresh) then
          List(i) = .false.
       else
          List(i) = .true.
       end if
    end do



  end subroutine ListOcc

  !>Subroutine that counts the number of elements of a vector, vec,  greater than zero.
  !>Author:Ida-Marie Hoeyvik
  subroutine CountNonZeroElements(vec, NonZeroNumb,MyAtom, natoms)
    implicit none
    real(realk), intent(in) :: vec(natoms)
    integer, intent(out)    :: NonZeroNumb
    integer, intent(in)     :: MyAtom, natoms
    real(realk), parameter  :: thresh=1E-20_realk
    integer                 :: i

    NonZeroNumb = 0
    do i=1,natoms
       if (vec(i) > thresh ) then
          NonZeroNumb = NonZeroNumb + 1
       end if
    end do



  end subroutine CountNonZeroElements




  !> Subroutine used to expand target space in fragment optimization
  !> \author Ida-Marie Hoeyvik
  subroutine ExpandTarget(EOSvector, TrackMat, natoms)
    implicit none
    logical, intent(inout) :: EOSvector(natoms)
    logical		:: TempEOSvector(natoms)
    integer, intent(in) 	:: natoms
    integer, intent(in)	:: TrackMat(natoms)
    integer               :: counter
    integer 		:: i,increase,AtomsOut,AtomsIn

    counter = 1
    call count_atoms(EOSvector,AtomsIn, natoms)
    AtomsOut = AtomsIn
    TempEOSvector= EOSvector
    do i=1,natoms
       if ((AtomsOut > (AtomsIn+8)) .or. (AtomsOut == natoms)) exit
       if (counter .le. (natoms-1)) then
          counter = counter + 1
          TempEOSvector(TrackMat(i)) = .true.
       end if
       call count_atoms(TempEOSvector, AtomsOut, natoms)
    end do

    EOSvector = TempEOSvector

    if (AtomsOut == (natoms-1)) EOSvector = .true.


  end subroutine ExpandTarget



  !>Subroutine to count the number of atoms in given space (Target, OccBuffer or UnoccBuffer).
  !>/author: Ida-Marie Hoeyvik
  subroutine count_atoms(atoms_vector, number_of_atoms,natoms)
    implicit none
    integer, intent(in)	:: natoms
    logical, intent(in) :: atoms_vector(natoms)
    integer, intent(out):: number_of_atoms
    integer		:: i

    number_of_atoms = 0

    do i = 1, natoms
       if ( atoms_vector(i) ) then
          number_of_atoms = number_of_atoms + 1
       end if
    end do

  end subroutine count_atoms



  !> \brief Sort with keeping track of the origial indices
  !> sort  elements of matrix, DistVecToSort,  from low to high.
  !> index matrix, TrackVec, to keep track on the position of the elements
  !> author: Ida-Marie Hoeyvik
  subroutine sort_track(DistVecToSort,TrackVec,natoms)

    implicit none
    integer, intent(in) :: natoms
    real(realk), dimension(natoms,natoms), intent(inout) :: DistVecToSort
    integer, dimension(natoms,natoms), intent(out) :: TrackVec
    real(realk) :: temp
    integer :: temp_track,i,j
    logical :: swp

    do j=1,natoms
       !create trackvec
       do i=1, natoms
          TrackVec(i,j)=i
       end do


       swp=.true.
       do while (swp)
          swp=.false.
          do i=1,natoms-1
             if(DistVecToSort(i,j) > DistVecToSort(i+1,j)) then

                temp = DistVecToSort(i+1,j)
                DistVecToSort(i+1,j) = DistVecToSort(i,j)
                DistVecToSort(i,j) = temp

                temp_track = TrackVec(i+1,j)
                TrackVec(i+1,j) = TrackVec(i,j)
                TrackVec(i,j) = temp_track

                swp=.true.
             end if
          end do
       end do
    end do

    return
  end subroutine sort_track


  !> \brief Sort with keeping track of the origial indices
  !> sort  elements of matrix, DistVecToSort,  from low to high.
  !> index matrix, TrackVec, to keep track on the position of the elements
  !> author: Ida-Marie Hoeyvik
  subroutine sort_track_vector(DistVecToSort,TrackVec,Vec_dim)

    implicit none
    integer, intent(in) :: Vec_dim
    real(realk) :: DistVecToSort(Vec_dim)
    integer :: TrackVec(Vec_dim)
    real(realk) :: temp
    integer :: temp_track,i
    logical :: swp


    swp=.true.
    do while (swp)
       swp=.false.
       do i=1,Vec_dim-1
          if(DistVecToSort(i) > DistVecToSort(i+1)) then

             temp = DistVecToSort(i+1)
             DistVecToSort(i+1) = DistVecToSort(i)
             DistVecToSort(i) = temp

             temp_track = TrackVec(i+1)
             TrackVec(i+1) = TrackVec(i)
             TrackVec(i) = temp_track

             swp=.true.
          end if
       end do
    end do


  end subroutine sort_track_vector


  !> \brief Adjust full molecular basis to a list of atoms and orbitals
  !> FullMatrix(nbasis,norbitals) -> SmallMatrix(nbasis_small,norbitals_small)
  !> for atomic indices in list atoms(natoms_small) using atoms_size(natoms),
  !> atoms_start(natoms), atoms_end(natoms)
  subroutine adjust_basis_matrix(FullMatrix,SmallMatrix,orbitals,atoms, &
       atom_size,atom_start,atom_end,nbasis,norbitals,natoms, &
       nbasis_small,norbitals_small,natoms_small)

    implicit none
    real(realk), dimension(nbasis,norbitals), intent(in) :: FullMatrix
    real(realk), dimension(nbasis_small,norbitals_small), &
         intent(inout) :: SmallMatrix
    integer, intent(in) :: nbasis,norbitals,nbasis_small,norbitals_small,natoms,natoms_small
    integer, dimension(norbitals_small), intent(in) :: orbitals
    integer, dimension(natoms_small), intent(in) :: atoms
    integer, dimension(natoms), intent(in) :: atom_size, atom_start, atom_end
    integer :: offset_start, offset_end, i,j

    do i=1,norbitals_small
       offset_start=1
       do j=1,natoms_small
          offset_end=offset_start+atom_size(atoms(j))-1
          SmallMatrix(offset_start:offset_end,i)= &
               FullMatrix(atom_start(atoms(j)):atom_end(atoms(j)),orbitals(i))
          offset_start=offset_end+1
       end do
    end do

  end subroutine adjust_basis_matrix

  !> \brief Adjust full molecular basis to a list of atoms and orbitals
  !> FullMatrix(nbasis,norbitals) -> SmallMatrix(nbasis_small,norbitals_small)
  !> for atomic indices in list atoms(natoms_small) using atoms_size(natoms),
  !> atoms_start(natoms), atoms_end(natoms)
  subroutine adjust_basis_matrix2(FullMatrix,SmallMatrix,orbitals,nbasis,&
       & norbitals,nbasis_small,norbitals_small,basis_idx)

    implicit none
    integer, intent(in) :: nbasis,norbitals,nbasis_small,norbitals_small
    real(realk), dimension(nbasis,norbitals), intent(in) :: FullMatrix
    real(realk), dimension(nbasis_small,norbitals_small), intent(inout) :: SmallMatrix
    integer, dimension(norbitals_small), intent(in) :: orbitals
    integer, dimension(nbasis_small), intent(in) :: basis_idx
    integer :: i,j

    do i=1,norbitals_small
       do j=1,nbasis_small
          SmallMatrix(j,i)= FullMatrix(basis_idx(j),orbitals(i))
       end do
    end do

  end subroutine adjust_basis_matrix2

  !> \brief Adjust full molecular suqare matrix to a list of atoms (AO matrix)
  subroutine adjust_square_matrix(FullMatrix,SmallMatrix,atoms,atom_size, &
       atom_start,atom_end,nbasis,natoms,nbasis_small,natoms_small)

    implicit none
    real(realk), dimension(nbasis,nbasis), intent(in) :: FullMatrix
    real(realk), dimension(nbasis_small,nbasis_small), intent(inout) :: SmallMatrix
    integer, intent(in) :: nbasis,nbasis_small,natoms,natoms_small
    integer, dimension(natoms_small), intent(in) :: atoms
    integer, dimension(natoms), intent(in) :: atom_size,atom_start,atom_end
    integer :: i_offset_start,i_offset_end, j_offset_start,j_offset_end, i,j

    i_offset_start=1
    do i=1,natoms_small
       i_offset_end=i_offset_start+atom_size(atoms(i))-1
       j_offset_start=1
       do j=1,natoms_small
          j_offset_end=j_offset_start+atom_size(atoms(j))-1

          SmallMatrix(i_offset_start:i_offset_end,j_offset_start:j_offset_end) = &
               FullMatrix(atom_start(atoms(i)):atom_end(atoms(i)),atom_start(atoms(j)):atom_end(atoms(j)))

          j_offset_start=j_offset_end+1
       end do
       i_offset_start=i_offset_end+1
    end do

  end subroutine adjust_square_matrix

  !> \brief Adjust full molecular suqare matrix to a list of AOs  (AO matrix)
  subroutine adjust_square_matrix2(FullMatrix,SmallMatrix,aos,nbasis,nbasis_small)

    implicit none
    real(realk), dimension(nbasis,nbasis), intent(in) :: FullMatrix
    real(realk), dimension(nbasis_small,nbasis_small), intent(inout) :: SmallMatrix
    integer, intent(in) :: nbasis,nbasis_small
    integer, dimension(nbasis_small), intent(in) :: aos
    integer :: i,j,BigJ
!    !$OMP PARALLEL DO DEFAULT(none) PRIVATE(J,BigJ,I) SHARED(nbasis_small,FullMatrix,SmallMatrix,aos)
    do j=1,nbasis_small
       BigJ = aos(j)
       do i=1,nbasis_small
          SmallMatrix(i,j) = FullMatrix(aos(i),BigJ)
       end do
    end do
!    !$OMP END PARALLEL DO
  end subroutine adjust_square_matrix2

  !> \brief Adjust full molecular basis to a list of orbitals (MO matrix)
  subroutine adjust_square_matrix_mo(FullMatrix,SmallMatrix,idx,norb,norb_small)
    implicit none
    real(realk), dimension(norb,norb), intent(in) :: FullMatrix
    real(realk), dimension(norb_small,norb_small), intent(inout) :: SmallMatrix
    integer, intent(in) :: norb, norb_small
    integer, dimension(norb_small), intent(in) :: idx
    integer :: i,j,idx_i,idx_j

    do j=1,norb_small
       idx_j = idx(j)
       do i=1,norb_small
          idx_i = idx(i)
          SmallMatrix(i,j) = FullMatrix(idx_i,idx_j)
       end do
    end do

  end subroutine adjust_square_matrix_mo

  !> \brief Get SubSystem indexes 
  subroutine GetSubSystemIndex(SubSystemIndex,nAtoms,mylsitem,int_output)
    implicit none
    type(lsitem), intent(inout) :: mylsitem
    integer, intent(in) :: nAtoms,int_output
    integer, intent(inout) :: SubSystemIndex(nAtoms)
    ! local variables
    integer :: i

    IF(mylsitem%input%molecule%nSubSystems>1)THEN
       do i=1,nAtoms
          SubSystemIndex(i)=mylsitem%input%molecule%ATOM(I)%SubSystemIndex
       end do
    ELSE
       SubSystemIndex=1
    ENDIF

  end subroutine GetSubSystemIndex

  !> \brief Get a table with interatomic distances (or interorbital for DECCO)
  subroutine GetDistances(MyMolecule,mylsitem,int_output)

    implicit none
    type(fullmolecule),intent(inout) :: MyMolecule
    type(lsitem), intent(inout) :: mylsitem
    integer, intent(in) :: int_output
    real(realk), pointer :: geometry(:,:)
    real(realk) :: dist
    integer :: i,j,k

    MyMolecule%DistanceTable=0.0E0_realk

    ! get geometry
    call mem_alloc(geometry,MyMolecule%nfrags,3)
    geometry=0.0E0_realk
    if(DECinfo%DECCO) then
       ! Geometry corresponds to center of charge for occupied MOs.
       do j=1,MyMolecule%nfrags
          do i=1,3
             geometry(j,i) = MyMolecule%carmomocc(i,j)
          end do
       end do
    else
       ! Atombased DEC (nfrags=natoms)
       call get_geometry(int_output,0,mylsitem%input%molecule,MyMolecule%nfrags,geometry(:,1), &
            geometry(:,2),geometry(:,3))
    end if

    do i=1,MyMolecule%nfrags
       do j=1,i
          dist=0.0E0_realk
          do k=1,3
             dist=dist+(geometry(i,k)-geometry(j,k))**2
          end do
          dist=sqrt(dist)
          MyMolecule%DistanceTable(i,j)=dist
          MyMolecule%DistanceTable(j,i)=dist
       end do
    end do

    call mem_dealloc(geometry)
  end subroutine GetDistances

  !> \brief distance between two points r1=(x1,y1,z1) and r2=(x2,y2,z2)
  !> \author Kasper Kristensen
  !> \date August 2011
  function get_distance_between_two_points(r1,r2) result(distance)
    implicit none
    !> Distance between points r1 and r2
    real(realk) :: distance
    !> Point 1
    real(realk), dimension(3),intent(in) :: r1
    !> Point 2
    real(realk), dimension(3),intent(in) :: r2
    integer :: i

    distance=0E0_realk
    do i=1,3
       distance = distance + ( r1(i)-r2(i) )**2
    end do
    distance=sqrt(distance)

  end function get_distance_between_two_points

  !> \brief Solve system of linear equations using DGESV lapack routine
  subroutine solve_linear_equations(A,x,B,n)

    implicit none
    real(realk), dimension(n,n), intent(in) :: A
    real(realk), dimension(n), intent(in) :: B
    real(realk), dimension(n), intent(inout) :: x
    integer, intent(in) :: n

    real(realk), pointer :: tmpA(:,:)
    real(realk), dimension(n) :: tmp
    integer, dimension(n) :: ipiv
    integer :: infoLAPACK
    external dgesv
    infoLAPACK = 0

    call mem_alloc(tmpA,n,n)  ! dgesv destroyes original A
    tmpA = A

    tmp = B
    x = 0.0E0_realk
    call dgesv(n,1,tmpA,n,ipiv,tmp,n,infoLAPACK)
    if(infoLAPACK /= 0) then
       print *, 'Error in DGESV: ', infoLAPACK
       call lsquit('solve_linear_equations: error in LAPACK DGESV routine',-1)
    end if
    x = tmp

    call mem_dealloc(tmpA)

    return
  end subroutine solve_linear_equations

  !> \brief Invert matrix
  subroutine invert_matrix(input,output,n)

    implicit none
    integer, intent(in) :: n
    real(realk), dimension(n,n), intent(in) :: input
    real(realk), dimension(n,n), intent(out) :: output
    real(realk), pointer :: output_full(:,:)

    real(realk), pointer :: work(:)
    integer, pointer :: ipiv(:)
    integer :: nrow,ncol
    integer :: info
    external dgetrf,dgetri

    call mem_alloc(output_full,n,n)
    call mem_alloc(work,n)
    call mem_alloc(ipiv,n)

    output_full = 0.0E0_realk
    output_full = input
    info=0
    call dgetrf(n,n,output_full,n,ipiv,info)
    if(info /= 0) then
       print *, 'info=', info
       call lsquit('error1 :: invert_matrix',-1)
    end if
    info=0
    call dgetri(n,output_full,n,ipiv,work,n,info)
    if(info /= 0) then
       print *, 'info=', info
       call lsquit('error2 :: invert_matrix',-1)
    end if

    output = 0.0E0_realk
    output = output_full

    call mem_dealloc(output_full)
    call mem_dealloc(work)
    call mem_dealloc(ipiv)

    return
  end subroutine invert_matrix



  !> \brief Simple function for getting the position in memory
  !> of the (a,b)th element of a matrix -
  !> i.e. idx = nrow*(b-1) + a
  !> \author Kasper Kristensen
  !> \date November 2010
  function get_matrix_position(a,b,nrow,ncol) result(idx)
    implicit none
    !> Position of element in memory
    integer :: idx
    !> Row index
    integer, intent(in) :: a
    !> Column index
    integer, intent(in) :: b
    !> Number of rows in matrix
    integer, intent(in) :: nrow
    !> Number of columns in matrix (only needed for sanity check)
    integer, intent(in) :: ncol

    idx = (b-1)*nrow + a

    ! Sanity check
    if(idx > nrow*ncol) then
       write(DECinfo%output,*) 'get_matrix_position: Index position exceeds matrix dimensions!'
       write(DECinfo%output,*) 'nrow*ncol :', nrow*ncol
       write(DECinfo%output,*) 'idx       :', idx
       write(DECinfo%output,*) 'a         :', a
       write(DECinfo%output,*) 'b         :', b
       write(DECinfo%output,*) 'nrow      :', nrow
       write(DECinfo%output,*) 'ncol      :', ncol
       call lsquit('get_matrix_position: Index number exceeds matrix dimensions!',-1)
    end if

  end function get_matrix_position


  !> \brief Get symmetric matrix A to power m.
  !> \author Kasper Kristensen
  !> \date February 2011
  subroutine get_power_of_symmetric_matrix(n,m,A,Am)
    implicit none

    !> Dimension of matrix A
    integer,intent(in) :: n
    !> Power to which A should be computed
    real(realk) :: m
    !> Matrix A
    real(realk),intent(in) :: A(n,n)
    !> Output matrix A^m
    real(realk),intent(inout) :: Am(n,n)
    real(realk),pointer :: S(:,:), C(:,:), eivalM(:,:), CeivalM(:,:)
    real(realk),pointer :: eival(:)
    integer :: i

    call mem_alloc(S,n,n)
    call mem_alloc(C,n,n)
    call mem_alloc(eivalM,n,n)
    call mem_alloc(CeivalM,n,n)
    call mem_alloc(eival,n)


    ! Solve eigenvalue problem for A with unit overlap
    ! ************************************************
    S=0.0E0_realk
    do i=1,n
       S(i,i) = 1.0E0_realk
    end do
    call solve_eigenvalue_problem(n,A,S,eival,C)



    ! Diagonal matrix with eigenvalues to power m on the diagonal
    ! ***********************************************************
    eivalM=0.0E0_realk
    do i=1,n
       eivalM(i,i) = eival(i)**m
    end do


    ! Calculate A^m:  A^m = C eival^m C^T
    ! ***********************************
    ! C eival^m
    call dgemm('n','n',n,n,n,1.0E0_realk,C,n, &
         & eivalM,n,0.0E0_realk,CeivalM,n)

    ! (C eival^m) C^T
    call dgemm('n','t',n,n,n,1.0E0_realk,CeivalM,n, &
         & C,n,0.0E0_realk,Am,n)


    call mem_dealloc(S)
    call mem_dealloc(C)
    call mem_dealloc(eivalM)
    call mem_dealloc(CeivalM)
    call mem_dealloc(eival)


  end subroutine get_power_of_symmetric_matrix





  subroutine ExcludeIfNoOrbs(Track,NewTrack,natoms,norb_per_atom)
    implicit none
    integer, intent(in)  :: natoms
    integer              :: Track(natoms),NewTrack(natoms),Recollect(natoms)
    integer              :: norb_per_atom(natoms)
    integer              :: i,counter,counter2


    NewTrack = 0
    Recollect = 0

    counter2 = 0
    counter  = 0
    do i=1,natoms
       if ((norb_per_atom(Track(i)) .ne. 0) .and.&
            & (counter .ne.natoms)) then
          counter = counter + 1
          NewTrack(counter) = Track(i)
       else
          counter2 = counter2 + 1
          Recollect(counter2) = Track(i)
       end if
    end do

    if ((counter2 .ne. 0).and. (counter2 .ne. natoms)) then
       do i = 1,natoms
          if (Recollect(i) == 0) exit
          counter = counter + 1
          NewTrack(counter) = Recollect(i)
       end do
    end if

  end subroutine ExcludeIfNoOrbs

  !Find maximum distance between MyAtom and an atom in EOS.
  ! EOSVec: UnoccEOS if occupied partitioning of corr.energy
  ! EOSvec: OccEOS if virtual partitioning of corr.energy
  subroutine FindMaxDistance(EOSVec,MyAtom,DistanceTable,natoms,MaxDist)
    implicit none
    integer, intent(in) :: natoms,MyAtom
    logical,intent(in)  :: EOSVec(natoms)
    real(realk)         :: MaxDist
    real(realk)         :: DistanceTable(natoms,natoms)
    integer             :: i
    real(realk)         :: tmp

    MaxDist=0.0E0_realk

    do i=1,natoms
       if (EOSVec(i)) then
          tmp = DistanceTable(i,MyAtom)
          if (tmp > MaxDist) MaxDist = tmp
       end if
    end do

    !Convert to Angstrom
    MaxDist = MaxDist/(1.88973)

  end subroutine FindMaxDistance

  subroutine ReduceBuffer(BufferVec,NewBuffer,DistTrackMatrix,MyAtom,natoms)
    implicit none
    integer, intent(in) :: natoms,MyAtom
    logical, intent(in) :: BufferVec(natoms)
    logical,intent(inout) :: NewBuffer(natoms)
    integer, intent(in) :: DistTrackMatrix(natoms,natoms)
    integer :: i, counter,indx

    NewBuffer = BufferVec
    counter = 0

    do i=natoms,1,-1
       indx = DistTrackMatrix(i,MyAtom)
       if (BufferVec(indx)) then
          NewBuffer(indx) = .false.
          counter = counter + 1
          if (counter == 3) exit
       end if
    end do

  end subroutine ReduceBuffer


  !> \brief Absorb logical vector 2 into logical vector 1
  !> such that at output vec1(i) is true if the original vec1(i) OR vec2(i) is true,
  !> while all vec2(i) entries are false at output.
  !> \author Kasper Kristensen
  !> \date August 2011
  subroutine absorb_logical_vector(n,vec1,vec2)
    implicit none

    !> Dimension of logical vectors
    integer,intent(in) :: n
    !> Vector 1
    logical,dimension(n), intent(inout) :: vec1
    !> Vector 2
    logical,dimension(n), intent(inout) :: vec2
    logical,dimension(n) :: tmp1,tmp2
    integer :: i

    tmp1(1:n)=vec1(1:n)
    tmp2(1:n)=vec2(1:n)
    vec1=.false.
    vec2=.false.
    do i=1,n
       if(tmp1(i) .or. tmp2(i)) then
          vec1(i) = .true.
       end if
    end do


  end subroutine absorb_logical_vector


  !> \brief Get logical pair vector where vec12(i) is true
  !> if the input vec1(i) OR vec2(i) is true.
  !> \author Kasper Kristensen
  !> \date August 2011
  subroutine get_logical_pair_vector(n,vec1,vec2,vec12)
    implicit none

    !> Dimension of logical vectors
    integer,intent(in) :: n
    !> Vector 1
    logical,dimension(n), intent(in) :: vec1
    !> Vector 2
    logical,dimension(n), intent(in) :: vec2
    !> Pair vector
    logical,dimension(n), intent(inout) :: vec12
    integer :: i

    vec12=.false.
    do i=1,n
       if(vec1(i) .or. vec2(i)) vec12(i) = .true.
    end do

  end subroutine get_logical_pair_vector



  !> \brief Subroutine that creates initial fragment
  !> \date august 2011
  !> \author Ida-Marie Hoyvik
  subroutine InitialFragment(natoms,nocc_per_atom,nvirt_per_atom,DistMyatom,&
       & init_Occradius,init_Virtradius,Occ,Virt)
    implicit none
    !> number of atoms in MOLECULE
    Integer, intent(in)    :: natoms
    !> Number of occupied / virtupied orbitals per atom
    integer,intent(in), dimension(natoms) :: nocc_per_atom, nvirt_per_atom
    !> Distances from central atom to other atoms
    real(realk),intent(in) :: DistMyAtom(natoms)
    !> Include Occ orbitals assigned to atoms within this distance of central atom
    real(realk),intent(in) :: init_Occradius
    !> Include Virt orbitals assigned to atoms within this distance of central atom
    real(realk),intent(in) :: init_Virtradius
    !> Which AOS atoms to include for occ and virt spaces
    !> (entry i is T if orbitals on atom i is included)
    !> (In practice occ and virt will be identical but we keep it general)
    logical,intent(inout)  :: Occ(natoms),Virt(natoms)
    integer                :: i
    real(realk)            :: FOT

    FOT=DECinfo%FOT

    write(DECinfo%output,'(a,f5.2)') " FOP Occ Radius for initial fragment: ",init_Occradius
    write(DECinfo%output,'(a,f5.2)') " FOP Virt Radius for initial fragment: ",init_Virtradius

    ! Include atoms within init_radius
    Occ=.false.
    Virt=.false.
    do i=1,natoms

       ! Skip if no orbitals are assigned
       if (DistMyAtom(i) .le. init_Occradius .and. (nocc_per_atom(i)/=0)) then
          Occ(i) = .true.
       end if
       if (DistMyAtom(i) .le. init_Virtradius .and. (nvirt_per_atom(i)/=0)) then
          Virt(i) = .true.
       end if

    end do

  end subroutine InitialFragment



  !> \brief Read matrix from file using simple Fortran array.
  !> Thus, if SYS_REAL is set (realk=4), then the matrix on file is
  !> in double precision, while the output matrix is in single precision.
  !> \author Kasper Kristensen
  !> \date September 2011
  subroutine dec_read_mat_from_file(filename,dim1,dim2,mat)

    implicit none
    !> Name of file
    character(*), intent(in) :: filename
    !> Dimension 1 of matrix
    integer, intent(in) :: dim1
    !> Dimension 2 of matrix
    integer, intent(in) :: dim2
    !> Output matrix
    real(realk), intent(inout), dimension(dim1,dim2) :: mat
    !real(realk), dimension(dim1,dim2) :: tmp
    real(realk),pointer :: tmp(:,:)
    logical :: file_exist
    integer :: funit, i,j
    integer(kind=long) :: i64,j64


    file_exist=.false.
    inquire(file=filename,exist=file_exist)
    if(file_exist) then
      call mem_alloc(tmp,dim1,dim2)

       funit=-1
       call lsopen(funit,filename,'OLD','UNFORMATTED')

       ! Read dimensions stored on file
       ! files always written using 64 bit integers
       read (funit) i64,j64
       i=int(i64,4)
       j=int(j64,4)

       ! Sanity check
       if( (i /= dim1) .or. (j/=dim2) ) then
          write(DECinfo%output,*) 'Filename            : ', filename
          write(DECinfo%output,*) 'Input dims          : ', dim1,dim2
          write(DECinfo%output,*) 'Dims read from file : ', i,j
          call lsquit('Dec_read_mat_from_file: Something wrong with dimensions',DECinfo%output)
       end if

       ! Read values from file
       read(funit) tmp
       call lsclose(funit,'KEEP')

       ! Convert values to realk format and store in output matrix
       mat = real(tmp,realk)

    else
       write(DECinfo%output,*) 'File does not exist: ', filename
       call lsquit('dec_read_mat_from_file: File does not exist',DECinfo%output)
    end if
    call mem_dealloc(tmp)

  end subroutine dec_read_mat_from_file




  !> \brief Calculate simple matrix product
  !> (without any update or multiplication factors) by calling dgemm:
  !> C = op(A) op(B)
  !> where A and B are fortran arrays and
  !> op(A) is either A og A^T (transpose), depending on the input
  !> \author Kasper Kristensen
  !> \date November 2011
  subroutine dec_simple_dgemm(m,k,n,A,B,C,TransA,TransB,use_thread_safe)
    implicit none
    !> Number of rows in op(A) = Number of rows in C
    integer,intent(in) :: m
    !> Number of columns in op(A) = Number of rows in op (B)
    integer,intent(in) :: k
    !> Number of columns in op(B) = Number of columns in C
    integer,intent(in) :: n
    !> Input matrix A in matrix product C = op(A) op(B)
    real(realk),intent(in),dimension(*) :: A
    !> Input matrix B in matrix product C = op(A) op(B)
    real(realk),intent(in),dimension(*) :: B
    !> Output matrix C in matrix product C = op(A) op(B)
    real(realk),intent(inout),dimension(*) :: C
    !> Transpose A [TransA='t' and op(A)=A^T] or not [TransA='n' and op(A)=A]
    character(len=1), intent(in) ::TransA
    !> Transpose B [TransB='t' and op(B)=B^T] or not [TransB='n' and op(B)=B]
    character(len=1), intent(in) ::TransB
    !> Enforce thread safe version of dgemm?
    !> (encouraged if dgemm is called from inside OMP loop)
    logical,intent(in),optional :: use_thread_safe
    integer :: lda,ldb,ldc
    logical :: ts


    ! Init stuff
    ! **********

    ! Use thread safe dgemm (careful with optional arguments)
    ts=.false.
    if(present(use_thread_safe)) then
       if(use_thread_safe) then
          ts=.true.
       end if
    end if

    ! Leading dimension of A
    if( TransA=='n' .or. TransA=='N' ) then
       lda = m
    elseif( TransA=='t' .or. TransA=='T' ) then
       lda = k
    else
       call lsquit('dec_simple_dgemm: TransA must be n,N,t, or T!', DECinfo%output)
    end if

    ! Leading dimension of B
    if( TransB=='n' .or. TransB=='N' ) then
       ldb = k
    elseif( TransB=='t' .or. TransB=='T' ) then
       ldb = n
    else
       call lsquit('dec_simple_dgemm: TransB must be n,N,t, or T!', DECinfo%output)
    end if

    ! Leading dimension of C
    ldc = m

    ! Call dgemm to calculate C = A B
    ! *******************************
    if(ts) then ! use threadsafe version
       call DGEMM_TS(TransA,TransB,m,n,k,1.0E0_realk,A,lda,B,ldb,0.0E0_realk,C,ldc)
    else ! use standard dgemm
       call DGEMM(TransA,TransB,m,n,k,1.0E0_realk,A,lda,B,ldb,0.0E0_realk,C,ldc)
    end if


  end subroutine dec_simple_dgemm


  !> \brief Add matrix product to existing matrix:
  !> C = C + op(A) op(B)
  !> where A and B are fortran arrays and
  !> op(A) is either A og A^T (transpose), depending on the input
  !> \author Kasper Kristensen
  !> \date November 2011
  subroutine dec_simple_dgemm_update(m,k,n,A,B,C,TransA,TransB)
    implicit none
    !> Number of rows in op(A) = Number of rows in C
    integer,intent(in) :: m
    !> Number of columns in op(A) = Number of rows in op (B)
    integer,intent(in) :: k
    !> Number of columns in op(B) = Number of columns in C
    integer,intent(in) :: n
    !> Input matrix A in matrix product C = op(A) op(B)
    real(realk),intent(in),dimension(*) :: A
    !> Input matrix B in matrix product C = op(A) op(B)
    real(realk),intent(in),dimension(*) :: B
    !> Output matrix C in matrix product C = op(A) op(B)
    real(realk),intent(inout),dimension(*) :: C
    !> Transpose A [TransA='t' and op(A)=A^T] or not [TransA='n' and op(A)=A]
    character(len=1), intent(in) ::TransA
    !> Transpose B [TransB='t' and op(B)=B^T] or not [TransB='n' and op(B)=B]
    character(len=1), intent(in) ::TransB
    integer :: lda,ldb,ldc



    ! Init stuff
    ! **********

    ! Leading dimension of A
    if( TransA=='n' .or. TransA=='N' ) then
       lda = m
    elseif( TransA=='t' .or. TransA=='T' ) then
       lda = k
    else
       call lsquit('dec_simple_dgemm: TransA must be n,N,t, or T!', DECinfo%output)
    end if

    ! Leading dimension of B
    if( TransB=='n' .or. TransB=='N' ) then
       ldb = k
    elseif( TransB=='t' .or. TransB=='T' ) then
       ldb = n
    else
       call lsquit('dec_simple_dgemm: TransB must be n,N,t, or T!', DECinfo%output)
    end if

    ! Leading dimension of C
    ldc = m

    ! Call dgemm to calculate C = A B
    ! *******************************
    call DGEMM(TransA,TransB,m,n,k,1.0E0_realk,A,lda,B,ldb,1.0E0_realk,C,ldc)


  end subroutine dec_simple_dgemm_update


  !> \brief Calculate how much memory is used for the fragment.
  !> Only two-dimensional arrays are considered and memory for
  !> vectors and single elements are ignored.
  !> \author Kasper Kristensen
  !> \date February 2011
  subroutine calculate_fragment_memory(MyFragment,fragmem)

    implicit none
    !> Fragment information
    type(decfrag),intent(inout) :: MyFragment
    real(realk), intent(inout) :: fragmem
    real(realk) :: O,V,A,tmp,GB


    ! Number of occupied (O), Virtual (V), atomic basis functions (A)
    ! ***************************************************************
    O = MyFragment%noccAOS
    V = MyFragment%nvirtAOS
    A = MyFragment%nbasis
    GB = 1.000E9_realk ! 1 GB

    ! Use type decfrag to calculate memory use
    ! ***************************************

    ! Y matrices
    fragmem = 2E0_realk*V*A + 2E0_realk*O*A

    ! Fock matrix in AO basis
    tmp = A*A
    fragmem = fragmem + tmp

    ! ppfock and qqfock
    tmp = O*O + V*V
    fragmem = fragmem + tmp

    ! Convert to GB
    fragmem = realk*fragmem/GB

  end subroutine calculate_fragment_memory



  !> \brief For a logical vector of length natoms (the number of atoms in the molecule),
  !> count the number "true" entries which are NOT hydrogen atoms.
  !> \author Kasper Kristensen
  !> \date November 2011
  function count_number_of_nonhydrogen_atoms(mylsitem,natoms,vec) result(n)

    implicit none
    !> Number of true entries in "vec" which are NOT hydrogen atoms
    integer :: n
    !> LS item info for full molecule
    type(lsitem), intent(inout) :: mylsitem
    !> Number of atoms in the molecule
    integer,intent(in) :: natoms
    !> Logical vector in question
    logical,dimension(natoms) :: vec
    integer :: i, atomnumber

    n=0
    do i=1,natoms

       ! Atomic number for atomic index "i"
       atomnumber = MyLsitem%input%molecule%atom(i)%atomic_number
       if(vec(i) .and. (atomnumber/=1) ) then
          ! Increase counter if vector index "i" is true and atom is NOT hydrogen
          n=n+1
       end if

    end do

  end function count_number_of_nonhydrogen_atoms

  !> \brief Initialize simple pointer structure
  !> \author Kasper Kristensen
  !> \date May 2012
  subroutine mypointer_init(arrdim,arr,start,N,thepointer)

    implicit none
    !> Dimension for reference array
    integer(kind=long),intent(in) :: arrdim
    !> Larger array where the pointer should point to elements "start to start+N-1"
    real(realk), dimension(arrdim),target :: arr
    !> Start index in larger array that the pointer points to
    integer(kind=long),intent(inout) :: start
    !> Size of my pointer
    integer(kind=long),intent(in) :: N
    !> Simple pointer structure
    type(mypointer),intent(inout) :: thepointer
    integer(kind=long),parameter :: IL512=512,IL1=1
#ifndef VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN
    if(start < 1) then
       print *, 'start = ', start
       call lsquit('mypointer_init: start value is smaller than 1!',-1)
    end if

    ! Ensure that we start at optimal place in memory: 512+integer
    do while(mod(start-IL1,IL512) /= 0)
       start=start+1
    end do

    ! Set simple pointer info
    thepointer%start = start
    thepointer%N = N
    thepointer%end = start + N -1

    ! Sanity check
    if(thepointer%end > arrdim) then
       print *, 'end   = ', thepointer%end
       print *, 'Array size = ', arrdim
       call lsquit('Pointer dimensions exceed actual array!',-1)
    end if

    ! Make pointer point to requested values
    nullify(thepointer%p)
    thepointer%p => arr(thepointer%start:thepointer%end)
#else
    call mem_alloc(thepointer%p,N)
    ! Set simple pointer info
    thepointer%start = 1
    thepointer%N = N
    thepointer%end = start + N -1
#endif

  end subroutine mypointer_init



  !> \brief Start up PAPI flop counting
  !> Assumes that mypapi_init has been called with global "eventset" parameter.
  !> \author Kasper Kristensen
  subroutine start_flop_counter()

    implicit none
    integer :: retval

    retval=0
#ifdef VAR_PAPI
    call PAPIf_start(eventset, retval)
#endif
    call init_FLOPonGPUaccouting()

  end subroutine start_flop_counter


  !> \brief Stop and read PAPI FLOP counter.
  !> If LSDALTON is not linked to PAPI, we just return 0.
  !> \author Kasper Kristensen
  subroutine end_flop_counter(flops,FLOPonGPU)

    implicit none
    !> Flops (simplest to save as real)
    real(realk),intent(inout) :: flops
    !> Flop executed on the GPU
    real(realk),intent(inout) :: FLOPonGPU
    !> "FLOPS" used by PAPI must be hardcoded 64 bit integer
    integer(kind=8) :: flops_int
    integer :: retval

    flops_int=0
    retval=0
#ifdef VAR_PAPI
    call PAPIf_stop(eventset,flops_int,retval)
#endif
    flops = real(flops_int)
    call extract_FLOPonGPUaccouting(FLOPonGPU)

  end subroutine end_flop_counter



  !> \brief Print fragment information to file unit
  !> referenced by DECinfo%output.
  !> \author Kasper Kristensen
  subroutine fragment_print(MyFragment,printlevel)

    implicit none
    !> Atomic (or pair) fragment
    type(decfrag),intent(inout) :: MyFragment
    !> How much to print?
    !> 1. Basic fragment info (size and energies, i.e. no pointers are printed)
    !> 2. Basic fragment info AND EOS/AOS indices
    !> 3. Print almost all information including MO coeff and Fock matrix
    integer,intent(in) :: printlevel
    logical :: pair
    integer :: i

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) '*************************************************************'
    write(DECinfo%output,*) '*                    FRAGMENT INFORMATION                   *'
    write(DECinfo%output,*) '*************************************************************'
    write(DECinfo%output,*)
    if(MyFragment%nEOSatoms==2) then
       pair=.false.
       write(DECinfo%output,'(1X,a,i6)') 'SINGLE FRAGMENT, atomic number:', MyFragment%EOSatoms(1)
    else
       pair=.true.
       write(DECinfo%output,'(1X,a,2i6)') 'PAIR FRAGMENT, atomic numbers:', &
            & MyFragment%EOSatoms(1),MyFragment%EOSatoms(2)
       write(DECinfo%output,'(1X,a,f12.2)') 'Pair distance = ', MyFragment%pairdist
    end if
    write(DECinfo%output,'(1X,a,i8)') 'Size occ EOS   =', MyFragment%noccEOS
    write(DECinfo%output,'(1X,a,i8)') 'Size virt EOS  =', MyFragment%nvirtEOS
    write(DECinfo%output,'(1X,a,i8)') 'Size occ AOS   =', MyFragment%noccAOS
    write(DECinfo%output,'(1X,a,i8)') 'Size virt AOS  =', MyFragment%nvirtAOS
    write(DECinfo%output,'(1X,a,i8)') '#core orbitals =', MyFragment%ncore
    write(DECinfo%output,'(1X,a,i8)') 'Red Frag: Size occ AOS  =', MyFragment%noccAOS
    write(DECinfo%output,'(1X,a,i8)') 'Red Frag: Size virt AOS =', MyFragment%nvirtAOS
    write(DECinfo%output,'(1X,a,i8)') 'Atoms in atomic extent =', MyFragment%natoms
    write(DECinfo%output,'(1X,a,i8)') 'Number of basis functions =', MyFragment%nbasis
    do i=1,ndecenergies
       write(DECinfo%output,'(1X,a,i5,f16.10)') 'Idx, Fragenergy ', i,MyFragment%energies(i)
    end do

    write(DECinfo%output,*)
    write(DECinfo%output,*) 'Fragment-adapted orbital (FO) information:'
    write(DECinfo%output,*) 'Using FOs: ', MyFragment%fragmentadapted
    write(DECinfo%output,'(1X,a,i8)') 'Size occ FOs   = ', MyFragment%noccFA
    write(DECinfo%output,'(1X,a,i8)') 'Size virt FOs = ', MyFragment%nvirtFA
    write(DECinfo%output,'(1X,a,g12.2)') 'Occ rejection thr   ', MyFragment%RejectThr(2)
    write(DECinfo%output,'(1X,a,g12.2)') 'Unocc rejection thr ', MyFragment%RejectThr(1)
    write(DECinfo%output,*)


    if(MyFragment%pairfrag) then
       write(DECinfo%output,*) 'This is a pair fragment!'
       write(DECinfo%output,*) 'Number of EOS atoms =', MyFragment%nEOSatoms
       write(DECinfo%output,*) 'List of EOS atoms = ', MyFragment%EOSatoms
    else
       write(DECinfo%output,*) 'This is a standard fragment!'
    end if

    if(printlevel>1) then

       write(DECinfo%output,*) 'Occ EOS indices (frag,full)'
       do i=1,MyFragment%noccEOS
          write(DECinfo%output,*) i, MyFragment%occEOSidx(i)
       end do
       write(DECinfo%output,*)

       write(DECinfo%output,*) 'Virt EOS indices (frag,full)'
       do i=1,MyFragment%nvirtEOS
          write(DECinfo%output,*) i, MyFragment%virtEOSidx(i)
       end do
       write(DECinfo%output,*)

       write(DECinfo%output,*) 'Occ AOS indices (frag,full)'
       do i=1,MyFragment%noccAOS
          write(DECinfo%output,*) i, MyFragment%occAOSidx(i)
       end do
       write(DECinfo%output,*)

       write(DECinfo%output,*) 'Virt AOS indices (frag,full)'
       do i=1,MyFragment%nvirtAOS
          write(DECinfo%output,*) i, MyFragment%virtAOSidx(i)
       end do
       write(DECinfo%output,*)

       write(DECinfo%output,*) 'Atomic extent indices (frag,full in terms of atoms)'
       do i=1,MyFragment%natoms
          write(DECinfo%output,*) i, MyFragment%atoms_idx(i)
       end do
       write(DECinfo%output,*)

    end if

    if(printlevel > 2 .and. MyFragment%BasisInfoIsSet) then

       write(DECinfo%output,*) 'Occ MO coefficients (column, elements in column)'
       do i=1,MyFragment%noccAOS
          write(DECinfo%output,*) i, MyFragment%Co(:,i)
       end do
       write(DECinfo%output,*)

       write(DECinfo%output,*) 'Virt MO coefficients (column, elements in column)'
       do i=1,MyFragment%nvirtAOS
          write(DECinfo%output,*) i, MyFragment%Cv(:,i)
       end do
       write(DECinfo%output,*)

       write(DECinfo%output,*) 'AO fock matrix (column, elements in column)'
       do i=1,MyFragment%nbasis
          write(DECinfo%output,*) i, MyFragment%fock(:,i)
       end do
       write(DECinfo%output,*)

       write(DECinfo%output,*) 'Occ-occ fock matrix (column, elements in column)'
       do i=1,MyFragment%noccAOS
          write(DECinfo%output,*) i, MyFragment%ppfock(:,i)
       end do
       write(DECinfo%output,*)

       write(DECinfo%output,*) 'Virt-virt fock matrix (column, elements in column)'
       do i=1,MyFragment%nvirtAOS
          write(DECinfo%output,*) i, MyFragment%qqfock(:,i)
       end do
       write(DECinfo%output,*)

       if(MyFragment%ncore>0) then
          write(DECinfo%output,*) 'Core-core fock matrix (column, elements in column)'
          do i=1,MyFragment%ncore
             write(DECinfo%output,*) i, MyFragment%ccfock(:,i)
          end do
          write(DECinfo%output,*)
       end if

       if(MyFragment%CDset) then
          write(DECinfo%output,*) 'Occ correlation density (column, elements in column)'
          do i=1,MyFragment%noccAOS
             write(DECinfo%output,*) i, MyFragment%OccMat(:,i)
          end do
          write(DECinfo%output,*)

          write(DECinfo%output,*) 'Virt correlation density (column, elements in column)'
          do i=1,MyFragment%nvirtAOS
             write(DECinfo%output,*) i, MyFragment%VirtMat(:,i)
          end do
          write(DECinfo%output,*)
       end if

       if(MyFragment%FAset) then
          write(DECinfo%output,*) 'Occupied FO coefficients (column, elements in column)'
          do i=1,MyFragment%noccFA
             write(DECinfo%output,*) i, MyFragment%CoFA(:,i)
          end do
          write(DECinfo%output,*)

          write(DECinfo%output,*) 'Virtual FO coefficients (column, elements in column)'
          do i=1,MyFragment%nvirtFA
             write(DECinfo%output,*) i, MyFragment%CvFA(:,i)
          end do
          write(DECinfo%output,*)

          if(.not. MyFragment%pairfrag) then
             write(DECinfo%output,*) 'Occupied corrdens eigenvalues'
             do i=1,MyFragment%noccFA
                write(DECinfo%output,*) i, MyFragment%CDocceival(i)
             end do
             write(DECinfo%output,*)

             write(DECinfo%output,*) 'Unoccupied corrdens eigenvalues'
             do i=1,MyFragment%nvirtFA
                write(DECinfo%output,*) i, MyFragment%CDvirteival(i)
             end do
             write(DECinfo%output,*)

          end if

       end if

    end if

    write(DECinfo%output,*) '*************************************************************'
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*)

  end subroutine fragment_print


  !> \brief Delete atomic fragment
  !> \author Marcin Ziolkowski
  subroutine atomic_fragment_free(fragment)

    implicit none
    !> Atomic fragment to be freed
    type(decfrag),intent(inout) :: fragment

    ! Free everything in fragment - including basis info (fock matrix, MO coefficients, lsitem etc.)
    call atomic_fragment_free_simple(fragment)
    call atomic_fragment_free_basis_info(fragment)
    call free_fragment_t1(fragment)
   
  end subroutine atomic_fragment_free

  !> \brief Delete the "simple" part of the atomic fragment structure, i.e.
  !> everything not deleted in atomic_fragment_free_basis_info.
  !> \author Kasper Kristensen
  subroutine atomic_fragment_free_simple(fragment)

    implicit none
    !> Atomic fragment to be freed
    type(decfrag),intent(inout) :: fragment
    integer :: i

    deallocate(fragment%noccLOC)
    deallocate(fragment%nvirtLOC)
    deallocate(fragment%noccFA)
    deallocate(fragment%nvirtFA)
    nullify(fragment%noccLOC)
    nullify(fragment%nvirtLOC)
    nullify(fragment%noccFA)
    nullify(fragment%nvirtFA)
    nullify(fragment%noccAOS)
    nullify(fragment%nvirtAOS)

    if(associated(fragment%occEOSidx)) then
       call mem_dealloc(fragment%occEOSidx)
       fragment%occEOSidx => null()
    end if

    if(associated(fragment%occAOSidx)) then
       call mem_dealloc(fragment%occAOSidx)
       fragment%occAOSidx => null()
    end if

    if(associated(fragment%virtEOSidx)) then
       call mem_dealloc(fragment%virtEOSidx)
       fragment%virtEOSidx => null()
    end if

    if(associated(fragment%virtAOSidx)) then
       call mem_dealloc(fragment%virtAOSidx)
       fragment%virtAOSidx => null()
    end if


    if(associated(fragment%coreidx)) then
       call mem_dealloc(fragment%coreidx)
       fragment%coreidx => null()
    end if

    ! indices
    if(associated(fragment%idxo)) then
       call mem_dealloc(fragment%idxo)
       fragment%idxo => null()
    end if

    if(associated(fragment%idxu)) then
       call mem_dealloc(fragment%idxu)
       fragment%idxu => null()
    end if


    if(associated(fragment%EOSatoms)) then
       call mem_dealloc(fragment%EOSatoms)
       fragment%EOSatoms => null()
    end if

    if(associated(fragment%OccContribs)) then
       call mem_dealloc(fragment%OccContribs)
       fragment%OccContribs => null()
    end if

    if(associated(fragment%VirtContribs)) then
       call mem_dealloc(fragment%VirtContribs)
       fragment%VirtContribs => null()
    end if

    if(fragment%CDset) then
       call mem_dealloc(fragment%OccMat)
       call mem_dealloc(fragment%VirtMat)
       fragment%CDset=.false.
    end if

    if(fragment%FAset) then
       call mem_dealloc(fragment%CoFA)
       call mem_dealloc(fragment%CvFA)

       !Check if associated because that might be different on master and slaves
       !due to the "Expensive Box"
       if(associated(fragment%ppfockFA))then
         call mem_dealloc(fragment%ppfockFA)
       endif

       if(associated(fragment%qqfockFA))then
         call mem_dealloc(fragment%qqfockFA)
       endif

       if(.not. fragment%pairfrag) then
          call mem_dealloc(fragment%CDocceival)
          call mem_dealloc(fragment%CDvirteival)
       end if
       fragment%FAset=.false.
    end if

    if(associated(fragment%atoms_idx)) then
       call mem_dealloc(fragment%atoms_idx)
    end if

    if(associated(fragment%basis_idx)) then
       call mem_dealloc(fragment%basis_idx)
    end if

    if(associated(fragment%cabsbasis_idx)) then
       call mem_dealloc(fragment%cabsbasis_idx)
    end if

    ! DEC orbitals
    if(associated(fragment%occAOSorb)) then
       do i=1,size(fragment%occAOSorb)
          call orbital_free(fragment%occAOSorb(i))
       end do
       call mem_dealloc(fragment%occAOSorb)
    end if

    if(associated(fragment%virtAOSorb)) then
       do i=1,size(fragment%virtAOSorb)
          call orbital_free(fragment%virtAOSorb(i))
       end do
       call mem_dealloc(fragment%virtAOSorb)
    end if

    ! Reduced fragment info
    do i=1,DECinfo%nFRAGSred
       call fragmentAOS_type_free(fragment%REDfrags(i))
    end do
    call mem_dealloc(fragment%REDfrags)

  end subroutine atomic_fragment_free_simple

  
  !> \brief Delete basis set infomation for atomic fragment (fock matrix, MO coefficients, lsitem etc.)
  !> \author Kasper Kristensen
  subroutine atomic_fragment_free_basis_info(fragment)

    implicit none
    !> Atomic fragment to be freed
    type(decfrag),intent(inout) :: fragment
    integer :: i

    if(associated(fragment%S)) then
       call mem_dealloc(fragment%S)
    end if

    ! Transformation matrices
    nullify(fragment%Co,fragment%Cv)
    if(associated(fragment%CoLOC)) then
       call mem_dealloc(fragment%CoLOC)
    end if

    if(associated(fragment%CvLOC)) then
       call mem_dealloc(fragment%CvLOC)
    end if

    ! Free CABS MOs !
    if(associated(fragment%Ccabs)) then
       call mem_dealloc(fragment%Ccabs)
    end if

    ! Free CABS RI MOs !
    if(associated(fragment%Cri)) then
       call mem_dealloc(fragment%Cri)
    end if

    if(associated(fragment%coreMO)) then
       call mem_dealloc(fragment%coreMO)
    end if

    if(associated(fragment%fock)) then
       call mem_dealloc(fragment%fock)
    end if

    ! Fock matrices
    nullify(fragment%ppfock,fragment%qqfock)
    if(associated(fragment%qqfockLOC)) then
       call mem_dealloc(fragment%qqfockLOC)
    end if

    if(associated(fragment%ppfockLOC)) then
       call mem_dealloc(fragment%ppfockLOC)
    end if

    if(associated(fragment%ccfock)) then
       call mem_dealloc(fragment%ccfock)
    end if

    ! delete dalton input
#ifdef VAR_MPI
    ! Quick fix such that lsitem is never handled for global master                                
    ! as this will destroy the overall MPI framework.
    if(infpar%mynum/=infpar%master) call ls_free(fragment%mylsitem)
#else
    call ls_free(fragment%mylsitem)
#endif

    ! Internal control of whether basis info is set or not
    fragment%BasisInfoIsSet=.false.

   if(DECinfo%F12) then
       call atomic_fragment_free_f12(fragment)
    end if
    
  end subroutine atomic_fragment_free_basis_info


  !> Free pointers in fragmentAOS type
  subroutine fragmentAOS_type_free(MyFragmentAOS)
    implicit none
    !> FragmentAOS information
    type(fragmentAOS),intent(inout) :: MyFragmentAOS

    if(associated(MyFragmentAOS%occAOSidx)) then
       call mem_dealloc(MyFragmentAOS%occAOSidx)
    end if
    if(associated(MyFragmentAOS%virtAOSidx)) then
       call mem_dealloc(MyFragmentAOS%virtAOSidx)
    end if

  end subroutine fragmentAOS_type_free



  !> \brief Free the f12 fragment free matrices
  !> \author Yang Min Wang
  subroutine atomic_fragment_free_f12(fragment)
    
    implicit none
    !> Atomic fragment to be freed
    type(decfrag),intent(inout) :: fragment
    
    ! Free CABS MOs !
    if(associated(fragment%Ccabs)) then
       call mem_dealloc(fragment%Ccabs)
    end if
    
    ! Free CABS RI MOs !
    if(associated(fragment%Cri)) then
       call mem_dealloc(fragment%Cri)
    end if

    if(associated(fragment%hJir)) then
       call mem_dealloc(fragment%hJir)
       fragment%hJir => null()
    end if
    
    if(associated(fragment%Krs)) then
       call mem_dealloc(fragment%Krs)
       fragment%Krs => null()
    end if
    
    if(associated(fragment%Frs)) then
       call mem_dealloc(fragment%Frs)
       fragment%Frs => null()
    end if

    if(associated(fragment%Fac)) then
       call mem_dealloc(fragment%Fac)
       fragment%Fac => null()
    end if

    if(associated(fragment%Frm)) then
       call mem_dealloc(fragment%Frm)
       fragment%Frm => null()
    end if

    if(associated(fragment%Fcp)) then
       call mem_dealloc(fragment%Fcp)
       fragment%Fcp => null()
    end if
    
    if(associated(fragment%Fij)) then
       call mem_dealloc(fragment%Fij)
       fragment%Fij => null()
    end if
    
  end subroutine atomic_fragment_free_f12


  !> \brief Free and nullify fragment information related to t1 amplitudes.
  !> \author Kasper Kristensen
  !> \date January 2012
  subroutine free_fragment_t1(Fragment)
    implicit none
    !> Fragment info (only t1-related information will be modified here)
    type(decfrag), intent(inout) :: Fragment

    if(associated(fragment%t1)) then
       call mem_dealloc(fragment%t1)
       fragment%t1 => null()
    end if

    if(associated(fragment%t1_occidx)) then
       call mem_dealloc(fragment%t1_occidx)
       fragment%t1_occidx => null()
    end if

    if(associated(fragment%t1_virtidx)) then
       call mem_dealloc(fragment%t1_virtidx)
       fragment%t1_virtidx => null()
    end if

  end subroutine free_fragment_t1


  !> \brief Destroy orbital
  subroutine orbital_free(myorbital)

    implicit none
    type(decorbital), intent(inout) :: myorbital

    if(associated(myorbital%aos)) then
       call mem_dealloc(myorbital%aos)
       myorbital%aos => null()
    end if

  end subroutine orbital_free


  !> \brief Determine memory for DEC calculation and store in DECinfo%memory.
  !> If memory was set manually in input, nothing is done here.
  !> Otherwise a system call is used to determine memory.
  !> \author Kasper Kristensen, modified by Pablo Baudin
  !> \date August 2012
  subroutine get_memory_for_dec_calculation()
    implicit none
    real(realk) :: mem, MemInUse
    logical :: memfound

    memfound=.false.
    if(DECinfo%memory_defined) then ! do nothing 
       write(DECinfo%output,'(1X,a,g12.4,a)') 'Memory set in input to be: ', DECinfo%memory, ' GB'

       ! sanity check
       MemInUse = 1.0E-9_realk*mem_allocated_global
       if (DECinfo%memory<MemInUse) then
          call get_available_memory(DECinfo%output,Mem,memfound)
          DECinfo%memory = Mem
          write(DECinfo%output,*) 'WARNING! Specified memory for DEC too small!'
          write(DECinfo%output,'(1X,a,g12.4,a)') 'Memory set by default to be: ', DECinfo%memory, ' GB'
       end if

    else ! using system call

       call get_available_memory(DECinfo%output,Mem,memfound)

       if(memfound) then
          DECinfo%memory = Mem
          write(DECinfo%output,'(1X,a,g12.4,a)') 'Memory found by system call to be: ', DECinfo%memory, ' GB'
       else
          write(DECinfo%output,*) 'WARNING! Memory could not be found using system call!'
          write(DECinfo%output,*) 'For optimal performace specify memory (in GB) manually using &
               & .Memory keyword in *DEC section.'
          write(DECinfo%output,'(1X,a,g12.4,a)') 'Memory set by default to be: ', DECinfo%memory, ' GB'
       end if

    end if

  end subroutine get_memory_for_dec_calculation

  !> \brief Get currently available memory for DEC calculation
  !> \author Kasper Kristensen
  !> \date August 2012
  subroutine get_currently_available_memory(mem)
    implicit none
    !> Currently available memory in GB
    real(realk),intent(inout) :: mem
    real(realk) :: MemInUse
    logical :: memfound

    if(DECinfo%use_system_memory_info)then
       ! Use the system available memory information accessible via
       ! /proc/meminfo or the top command on mac os X. Especially the latter is
       ! error prone and not recommended. The /proc/meminfo can be preferred
       ! over the input setting on some systems

       call get_available_memory(DECinfo%output,mem,memfound,.true.)

       if(.not.memfound)then
          call lsquit("ERROR(get_currently_available_memory):system call failed,&
          & use .MEMORY keyword to specify memory in LSDALTON.INP",-1)
       endif

    else
       ! Total memory was determined at the beginning of DEC calculaton (DECinfo%memory)
       ! Memory currently in use is mem_allocated_global (see memory.f90)
       ! Memory available: Total memory - memory in use
       ! Mem in use in GB
       MemInUse = 1.0E-9_realk*mem_allocated_global
       mem = DECinfo%memory - MemInUse

    endif

    ! Sanity check
    if(mem < 0.0E0_realk) then
       write(DECinfo%output,*) 'Out of memory or something wrong with memory bookkeping!'
       write(DECinfo%output,'(1X,a,g12.4)') 'Total memory:  ', DECinfo%memory
       write(DECinfo%output,'(1X,a,g12.4)') 'Memory in use: ', MemInUse
       write(DECinfo%output,*)
       write(DECinfo%output,*) 'I print memory statistics overview before quitting...'
       call stats_mem(DECinfo%output)
       call lsquit('get_currently_available_memory: &
            & Out of memory or something wrong with memory book keeping!',-1)
    end if


  end subroutine get_currently_available_memory



  !> Fit set of {x,y} data to function y ~ f(x) = a1*x^{-p1} + a2*x^{-p2} + ...
  !> where p1,p2,... are inputs, and a1,a2,... are to be determined.
  !> Intended to be used for pair energy regression estimates, for example:
  !> f(x) = a2*x^{-6} + a3*x^{-7} + a4*x^{-8}
  !> \author Kasper Kristensen
  !> \date October 2012
  subroutine dec_regression(m,xval,yval,nterms,powers,a,errornorm)
    implicit none
    !> Number of (x,y) points to fit
    integer,intent(in) :: m
    !> x values (will be distances for pair energy fits)
    real(realk),dimension(m),intent(in) :: xval
    !> y values (will be pair interaction energies for pair energy fits)
    real(realk),dimension(m),intent(in) :: yval
    !> Number of terms to include in function (will be 3 for f(x) = a1*x^{-6} + a2*x^{-7} + a3*x^{-8})
    integer,intent(in) :: nterms
    !> Which powers to use in expansion (will be -6.0,-7.0, and -8.0 for f(x) example above)
    real(realk),dimension(nterms),intent(in) :: powers
    !> Fitted parameters (a1,a2,a3 for f(x) above)
    real(realk),dimension(nterms),intent(inout) :: a
    !> Norm of difference between y data values and fitted function f(x)
    real(realk),intent(inout) :: errornorm
    real(realk),pointer :: X(:,:), XT(:,:), XTy(:),XTX(:,:)
    integer :: i,j
    real(realk) :: fx


    ! Init stuff
    ! **********
    call mem_alloc(X,m,nterms)
    call mem_alloc(XT,nterms,m)
    call mem_alloc(XTy,nterms)
    call mem_alloc(XTX,nterms,nterms)

    ! Generate set of vectors with powers of x values: X(i,j) = xval(i)**powers(j) 
    call dec_regression_get_powers(m,xval,nterms,powers,X,XT)



    !=============================================================
    !                   Linear regression scheme                 !
    !=============================================================

    ! Having generated X which contains powers of the x values, 
    ! the problem is reduced to a linear least squares fit, which corresponds to solving
    ! the following equation for a:
    !
    ! (X^T X) a = (X^T yval)
    !
    ! See e.g. http://en.wikipedia.org/wiki/Linear_least_squares_(mathematics).


    ! Calculate X^T yval,  dimension: (nterms,1)
    call dec_simple_dgemm(nterms,m,1,XT,yval,XTy,'N','N')

    ! Calculate X^T X,  dimension: (nterms,nterms)
    call dec_simple_dgemm(nterms,m,nterms,XT,X,XTX,'N','N')

    ! Solve equation to get a: (X^T X) a = (X^T yval)
    call solve_linear_equations(XTX,a,XTy,nterms)


    ! Calculate norm of difference between real values y and fitted function f(x)
    errornorm=0.0E0_realk
    do i=1,m

       ! f(x) = a1*x^{-p1} + a2*x^{-p2}
       fx = 0.0E0_realk
       do j=1,nterms
          ! Recall: X(i,j) = xval(i)**powers(j)
          fx = fx + a(j)*X(i,j)
       end do

       ! Update error norm: (yval-f(x))**2
       errornorm = errornorm + (yval(i)-fx)**2
    end do
    errornorm = sqrt(errornorm)


    call mem_dealloc(X)
    call mem_dealloc(XT)
    call mem_dealloc(XTy)
    call mem_dealloc(XTX)


  end subroutine dec_regression


  !> Get various powers of data points xval, generate both matrix with x power values itself and its transpose.
  !> \author Kasper Kristensen
  !> \date October 2012
  subroutine dec_regression_get_powers(m,xval,nterms,powers,X,XT)
    implicit none
    !> Number of (x,y) points to fit
    integer,intent(in) :: m
    !> x values 
    real(realk),dimension(m),intent(in) :: xval
    !> Number of different powers that we are interested in
    integer,intent(in) :: nterms
    !> Which powers to use in expansion 
    real(realk),dimension(nterms),intent(in) :: powers
    !> Powers of x value data points: X(i,j) = xval(i)**powers(j)
    real(realk),dimension(m,nterms),intent(inout) :: X
    !> Transpose of X (needed explicitly in dec_regression)
    real(realk),dimension(nterms,m),intent(inout) :: XT
    integer :: i,j

    ! X(i,j) = xval(i)**powers(j)
    do j=1,nterms
       do i=1,m
          X(i,j) = xval(i)**powers(j)
          XT(j,i) = X(i,j)
       end do
    end do


  end subroutine dec_regression_get_powers



  !> Make simple "ascii-art" plot of a set of points (x,y)
  !> \author Kasper Kristensen
  !> \date October 2012
  subroutine simple_ascii_plot(npoints,x,y,xlabel,ylabel)

    implicit none
    !> Number of points to plot
    integer,intent(in) :: npoints
    !> x data points
    real(realk),intent(in),dimension(npoints) :: x
    !> y data points
    real(realk),intent(in),dimension(npoints) :: y
    !> Label for x
    character(len=16), intent(in) :: xlabel
    !> Label for y
    character(len=16), intent(in) :: ylabel
    integer,parameter :: length = 82
    integer,parameter :: lengthx = 100
    integer,parameter :: height = 32
    character(len=length) :: string,tmp
    character(len=lengthx) :: xaxis
    integer :: i,j,xminplot,xmaxplot,yminplot,ymaxplot,xmidplot,ymidplot
    real(realk) :: xmax,xmin,ymax,ymin,xscaling,yscaling,xmid,ymid
    logical,dimension(length,height) :: plot
    integer :: xplot(npoints),yplot(npoints),xmaxidx(1),ymaxidx(1),xminidx(1),yminidx(1)

    ! Note: 
    ! To make a "pretty" plot several values have been hardcoded.
    ! If modifications are to be made, then make sure that 
    ! hardcoded values are changed consistently.

    ! Sanity check
    if(npoints<2) then  ! Skip plot if only 1 point
       write(DECinfo%output,*) 'Ascii-plot will be skipped because there is less than 2 points!'
       return
    end if



    ! Init string to print
    do i=1,length

       ! First elements will be part of yaxis
       if(i==14) then  ! y axis hardcoded to be in column 14
          string(i:i) = '*'
       else
          string(i:i) = ' '
       end if

    end do

    ! Special treatment for string for xaxis
    do i=1,lengthx
       if(i<=length .and. i>13) then  ! x axis
          xaxis(i:i)='*'
       else
          xaxis(i:i)=' '  ! make room for xlabel
       end if
    end do
    ! Set xlabel
    xaxis(lengthx-len(xlabel)+1:lengthx) = xlabel


    ! Minimum and maximum values of original points
    xmin = minval(x)
    xmax = maxval(x)
    ymin = minval(y)
    ymax = maxval(y)
    ! Index positions for min and max values
    xminidx=minloc(x)
    yminidx=minloc(y)
    xmaxidx=maxloc(x)
    ymaxidx=maxloc(y)

    ! Also find mid points
    xmid = 0.5E0_realk*(xmax-xmin) + xmin
    ymid = 0.5E0_realk*(ymax-ymin) + ymin

    ! Max and min values in actual plot 
    xminplot = 17   ! Start two spaces next two yaxis
    yminplot = 2    ! start 1 space above xaxis
    xmaxplot = length -3  ! A little empty space at the end of x axis
    ymaxplot = height -1  ! A little empty space at the end of y axis
    ! Midpoints
    xmidplot = (xmaxplot-xminplot)/2 + xminplot
    ymidplot = (ymaxplot-yminplot)/2 + yminplot

    ! Scaling factor to generate points to plot from original data
    xscaling = real(xmaxplot-xminplot)/(xmax-xmin)
    yscaling = real(ymaxplot-yminplot)/(ymax-ymin)


    ! Set logical vector where plot(k,l) is true if "pixel" (k,l) represent a data point
    plot=.false.
    do i=1,npoints
       ! x/y data values converted to actual plot using scaling factor
       xplot(i) = int(xscaling*(x(i)-xmin) + real(xminplot))
       yplot(i) = int(yscaling*(y(i)-ymin) + real(yminplot))

       ! Ensure that min and max values are represented 
       ! (i.e. purefy possible round-offs issues)
       if(i==xminidx(1)) xplot(xminidx(1)) = xminplot
       if(i==yminidx(1)) yplot(yminidx(1)) = yminplot
       if(i==xmaxidx(1)) xplot(xmaxidx(1)) = xmaxplot
       if(i==ymaxidx(1)) yplot(ymaxidx(1)) = ymaxplot

       plot(xplot(i),yplot(i)) = .true.
    end do


    ! Start plotting
    ! **************

    write(DECinfo%output,*) 
    write(DECinfo%output,*) 
    write(DECinfo%output,'(1X,a)') ylabel
    write(DECinfo%output,*) 

    ! Loop over all y values in plot 
    do j=height,1,-1

       ! Save default string
       tmp = string

       ! Insert plot point 'X' if plot points (i,j) represents data value
       do i=1,length
          if(plot(i,j)) string(i:i) = 'X'
       end do

       ! Write maximum y data value next to largest y value point
       if(j==ymaxplot) then
          write(string(1:12),'(g12.2)') ymax
       end if

       ! Write midpoint y value at the midpoint of plotted y axis
       if(j==ymidplot) then
          write(string(1:12),'(g12.2)') ymid
       end if

       ! Write minimum y data value at the bottom of y axis
       if(j==yminplot) then
          write(string(1:12),'(g12.2)') ymin
       end if

       ! Write string to output file
       write(DECinfo%output,'(1X,a)') string

       ! Restore default string
       string = tmp

    end do

    ! Write x axis
    write(DECinfo%output,'(1X,a)') xaxis


    ! Write min, mid, and max x data values
    do i=1,lengthx
       xaxis(i:i) = ' '
    end do
    write(xaxis(xminplot-6:xminplot+5),'(g12.2)') xmin
    write(xaxis(xmidplot-6:xmidplot+5),'(g12.2)') xmid
    write(xaxis(xmaxplot-6:xmaxplot+5),'(g12.2)') xmax
    write(DECinfo%output,'(1X,a)') xaxis


    write(DECinfo%output,*) 
    write(DECinfo%output,*) 


  end subroutine simple_ascii_plot



  !> \brief Calculate distance between two fragments defined as the
  !> smallest distance between an atom in fragment 1 and an atom in fragment 2.
  !> (Also works for standard fragments).
  !> \author Kasper Kristensen
  !> \date November 2011
  function get_distance_between_fragments(Fragment1,Fragment2,natoms,DistanceTable) result(fragdist)


    implicit none
    !> Distance between fragments
    real(realk) :: fragdist
    !>  fragment 1
    type(decfrag),intent(inout) :: Fragment1
    !>  fragment 2
    type(decfrag),intent(inout) :: Fragment2
    !> Number of atoms for full molecule
    integer, intent(in) :: natoms
    !> Distance table for all atoms in the molecule
    real(realk), intent(in) :: DistanceTable(natoms,natoms)
    integer :: i,j
    real(realk) :: dist

    ! In principle this can also give the "distance" between two pair fragments
    ! as the maximu distance between an atom in pair1 and an atom in pair2.
    fragdist = get_minimum_distance_between_groups_of_atoms(Fragment1%nEOSatoms,&
         & Fragment2%nEOSatoms, Fragment1%EOSatoms,  Fragment2%EOSatoms, natoms,DistanceTable)

  end function get_distance_between_fragments



  !> \brief Calculate the minimum distance between an atom in atomlist1 and
  !> an atom in atomlist2.
  !> \author Kasper Kristensen
  !> \date November 2011
  function get_minimum_distance_between_groups_of_atoms(dim1,dim2,atomlist1,atomlist2,&
       & natoms,DistanceTable) result(minimum_dist)


    implicit none
    !> Number of atoms in atomlist 1
    integer,intent(in) :: dim1
    !> Number of atoms in atomlist 2
    integer,intent(in) :: dim2
    !> List of atomic indices 1
    integer,intent(in) :: atomlist1(dim1)
    !> List of atomic indices 2
    integer,intent(in) :: atomlist2(dim2)
    !> Minimum distance between atoms in the two list above
    real(realk) :: minimum_dist
    !> Number of atoms for full molecule
    integer, intent(in) :: natoms
    !> Distance table for all atoms in the molecule
    real(realk), intent(in) :: DistanceTable(natoms,natoms)
    integer :: i,j
    real(realk) :: dist


    minimum_dist=huge(1.0)

    do i=1,dim1
       do j=1,dim2

          ! Distance between atom "i" in list 1 and atom "j" in list 2
          dist = DistanceTable(atomlist1(i),atomlist2(j))

          ! Find smallest distance
          if(dist < minimum_dist) minimum_dist = dist

       end do
    end do


  end function get_minimum_distance_between_groups_of_atoms


  !> \brief Initialize single precision gridbox (see type SPgridbox)
  !> \author Kasper Kristensen
  !> \date November 2012
  subroutine init_SPgridbox(center,delta,n,mygrid)
    implicit none
    !> Center of grid box
    real(4),intent(in) :: center(3)
    !> Distance between neighbouring points in grid box
    real(4),intent(in) :: delta
    !> Number of grid points in each x,y,z direction measured from center
    integer,intent(in) :: n
    !> Grid box
    type(SPgridbox),intent(inout) :: mygrid

    ! Basic info
    mygrid%center = center
    mygrid%delta = delta
    mygrid%n = n

    ! Number of points in each x,y, or z direction: 2n+1 (see type SPgridbox)
    mygrid%nd = 2*n+1
    call mem_alloc(mygrid%val,mygrid%nd,mygrid%nd,mygrid%nd)
    mygrid%val = 0.0_4

  end subroutine init_SPgridbox



  !> \brief Free single precision gridbox
  !> \author Kasper Kristensen
  !> \date November 2012
  subroutine free_SPgridbox(mygrid)
    implicit none
    !> Grid box
    type(SPgridbox),intent(inout) :: mygrid

    call mem_dealloc(mygrid%val)

  end subroutine free_SPgridbox

  !> Get density D = Cocc Cocc^T from occupied orbitals
  !> \author Kasper Kristensen, mod by PE
  !> \date November 2012
  subroutine get_density_from_occ_orbitals(nbasis,nocc,Cocc,dens,Cocc2)
    implicit none
    !> Number of basis functions
    integer,intent(in) :: nbasis
    !> Number of occupied orbitals (can be only valence orbitals for frozen core)
    integer,intent(in) :: nocc
    !> Occupied MO coefficients (can be only valence orbitals for frozen core)
    real(realk),intent(in),dimension(nbasis,nocc) :: Cocc
    real(realk),intent(in),dimension(nbasis,nocc),optional :: Cocc2
    !> Density
    real(realk),intent(inout),dimension(nbasis,nbasis) :: dens
    real(realk),pointer :: Cocc_copy(:,:)
    integer :: i,j

    if(present(Cocc2))then
       ! density = Cocc Cocc^T 
       call dec_simple_dgemm(nbasis,nocc,nbasis,Cocc,Cocc2,dens,'n','t')
    else
       ! Cocc copy (avoid passing the same element into dgemm twice)
       call mem_alloc(Cocc_copy,nbasis,nocc)
       Cocc_copy = Cocc

       ! density = Cocc Cocc^T 
       call dec_simple_dgemm(nbasis,nocc,nbasis,Cocc,Cocc_copy,dens,'n','t')
       call mem_dealloc(Cocc_copy)

    endif

  end subroutine get_density_from_occ_orbitals



  !> \brief Simple routine for getting HF energy from density and Fock matrix.
  !> Assumes closed-shell system
  !> \author Kasper Kristensen
  !> \date July 2011
  function get_HF_energy(D,F,Mylsitem) result(Ehf)

    implicit none
    !> HF energy
    real(realk) :: Ehf
    !> Density matrix (full molecule)
    type(matrix),intent(in) :: D
    !> Fock matrix (full molecule)
    type(matrix),intent(in) :: F
    !> LSDALTON INFO
    type(lsitem), intent(inout) :: MyLsitem
    type(matrix) :: h
    real(realk) :: Enuc,hD, FD, fac
    integer :: nbasis


    ! One electron matrix
    nbasis = D%nrow
    call mat_init(h,nbasis,nbasis)
    call mat_zero(h)
    call ii_get_h1(DECinfo%output,DECinfo%output,mylsitem%setting,h)

    ! Nuclear repulsion energy
    CALL II_get_nucpot(DECinfo%output, DECinfo%output,mylsitem%setting,Enuc)

    !Tr(hD)
    hD = mat_dotproduct(D,h)
    call mat_free(h)
    !Tr(FD)
    FD = mat_dotproduct(D,F)

    ! Since F=h+G(D) one "h" is already included and the HF energy is found from:
    ! Ehf = Tr(hD) + Tr(FD) + Enuc
    Ehf = hd + FD + Enuc

  end function get_HF_energy



  !> \brief Get HF energy from full molecule structure.
  !> Assumes closed-shell system
  !> \author Kasper Kristensen
  !> \date December 2012
  function get_HF_energy_fullmolecule(MyMolecule,Mylsitem,D) result(Ehf)

    implicit none
    !> HF energy
    real(realk) :: Ehf
    !> Full molecule structure
    type(fullmolecule),intent(in) :: MyMolecule
    !> LSDALTON INFO
    type(lsitem), intent(inout) :: MyLsitem
    !> HF Density matrix
    type(matrix),intent(in) :: D
    !> get E_DFT
    type(matrix) :: F,h
    real(realk)  :: exchangeFactor,enuc,edft(1)
    real(realk), pointer :: fock(:,:)

    if(DECinfo%noaofock) then
       write(DECinfo%output,*) 'Warning: NOFOCKAO keyword is set and HF energy is &
            & set to zero'
       Ehf = 0.0_realk
       return
    end if

    ! Init Fock matrix in matrix form
    call mat_init(F,MyMolecule%nbasis,MyMolecule%nbasis)



    if(DECinfo%DFTreference) then
      !Needs the fock matrix from the KS density
      exchangeFactor = mylsitem%SETTING%SCHEME%exchangeFactor
      mylsitem%SETTING%SCHEME%exchangeFactor=1.0_realk
      !This was from the beginning set to zero, has
      !to be 1.0 for evaluation of exact exchange.

      call mat_init(h,MyMolecule%nbasis,MyMolecule%nbasis)
      call mat_zero(h)
      call mat_zero(F)

      call II_get_h1(DECinfo%output,DECinfo%output,mylsitem%setting,h)
      call II_get_coulomb_and_exchange_mat(DECinfo%output, &
        & DECinfo%output,mylsitem%setting,D,F,1)
      !adds one-particle part
      call mat_daxpy(1.0_realk,h,F)
      !reset exchangeFactor for future use
      mylsitem%SETTING%SCHEME%exchangeFactor=exchangeFactor 

      call mat_free(h)
    else

       if( MyMolecule%mem_distributed )then
          call mem_alloc(fock,MyMolecule%nbasis,MyMolecule%nbasis)
          call tensor_gather(1.0E0_realk,MyMolecule%fock,0.0E0_realk,fock,MyMolecule%nbasis**2*i8)
       else
          fock => MyMolecule%fock%elm2
       endif

       call mat_set_from_full(fock, 1E0_realk, F)

       if( MyMolecule%mem_distributed )then
          call mem_dealloc( fock )
       endif

    endif

    ! Get HF energy
    Ehf = get_HF_energy(D,F,Mylsitem) 

    call mat_free(F)

  end function get_HF_energy_fullmolecule

  !> \brief Get KS energy from full molecule structure.
  !> Assumes closed-shell system
  !> \author J. Rekkedal
  !> \date December 2012
  function get_dft_energy_fullmolecule(MyMolecule,Mylsitem,D) result(Ehf)

    implicit none
    !> HF energy
    real(realk) :: Ehf
    !> Full molecule structure
    type(fullmolecule),intent(in) :: MyMolecule
    !> LSDALTON INFO
    type(lsitem), intent(inout) :: MyLsitem
    !> HF Density matrix
    type(matrix),intent(in) :: D
    !> get E_DFT
    type(matrix) :: F,K,h
    type(matrix) :: Dtmp(1)
    integer :: igrid
    real(realk)  :: DFTELS,enuc,edft(1)

    ! Init Fock matrix in matrix form
    call mat_init(F,MyMolecule%nbasis,MyMolecule%nbasis)
    call mat_init(K,MyMolecule%nbasis,MyMolecule%nbasis)
    call mat_init(h,MyMolecule%nbasis,MyMolecule%nbasis)
    call mat_init(Dtmp(1),MyMolecule%nbasis,MyMolecule%nbasis)
    call mat_copy(1.0_realk,D,Dtmp(1))
    call mat_zero(h)
    call mat_zero(F)
    call mat_zero(K)


    MyLsitem%setting%scheme%DFT%griddone = 0
    igrid = MyLsitem%setting%scheme%DFT%igrid
    MyLsitem%setting%scheme%DFT%GridObject(igrid)%GRIDDONE = 0
    DFTELS = MyLsItem%setting%scheme%DFT%DFTELS
    MyLsItem%setting%scheme%DFT%DFTELS = 100.0E0_realk

    call II_get_xc_energy(DECinfo%output,DECinfo%output,&
      & mylsitem%setting,MyMolecule%nbasis,Dtmp,Edft,1)

    MyLsItem%setting%scheme%DFT%DFTELS = DFTELS

    call mat_free(Dtmp(1))
    !Get one-electron part of the Fock Matrix
    call II_get_h1(DECinfo%output,DECinfo%output,mylsitem%setting,h)
    !Get the Coulomb matrix
    call II_get_coulomb_mat(DECinfo%output, &
      & DECinfo%output,mylsitem%setting,D,F,1)
    !Get the long-range exchange matrix for the range-separated case
    call II_get_exchange_mat(DECinfo%output, &
      & DECinfo%output,mylsitem%setting,D,1,.true.,K)
    
    call mat_daxpy(2.0_realk,h,F)
    CALL II_get_nucpot(DECinfo%output, DECinfo%output,mylsitem%setting,Enuc)
    Ehf =Edft(1) + mat_dotproduct(D,F) + Enuc + mat_dotproduct(D,K)
    
    call mat_free(h)
    call mat_free(K)
    call mat_free(F)


  end function get_dft_energy_fullmolecule




  !> Check whether fragment restart file exist.
  !> \author Kasper Kristensen
  !> \date December 2012
  function fragment_restart_file_exist(first_order,esti) result(file_exist)

    implicit none
    !> First order calculation?
    logical,intent(in) :: first_order
    !> Use estimated fragment energies
    logical,intent(in) :: esti
    logical :: file_exist
    character(len=40) :: FileName

    if(first_order) then  ! first order calculation
       filename = 'mp2grad.info'
    else ! energy calculation
       filename = get_fragenergy_restart_filename(esti) 
    end if

    inquire(file=FileName,exist=file_exist)


  end function fragment_restart_file_exist

  !> Read 64 bit integer from file and convert to 32 bit integer
  !> \author Kasper Kristensen
  !> \date December 2012
  subroutine read_64bit_to_32bit_singleinteger(funit,myint)
    implicit none
    !> File unit to read from (assumes file is already open)
    integer,intent(in) :: funit
    !> Output 32 bit integer
    integer,intent(inout) :: myint
    integer(kind=8) :: myint_long

    read(funit) myint_long
    myint = int(myint_long)

  end subroutine read_64bit_to_32bit_singleinteger

  !> Read 32 bit integer from file and convert to 64 bit integer
  !> \author Patrick Ettenhuber (adapted from Kasper)
  !> \date December 2012
  subroutine read_32bit_to_64bit_singleinteger(funit,myint)
    implicit none
    !> File unit to read from (assumes file is already open)
    integer,intent(in) :: funit
    !> Output 32 bit integer
    integer ,intent(inout) :: myint
    integer(kind=4) :: myint_long

    read(funit) myint_long
    myint = int(myint_long)

  end subroutine read_32bit_to_64bit_singleinteger


  !> Read 64 bit logical from file and convert to 32 bit logical
  !> \author Kasper Kristensen
  !> \date December 2012
  subroutine read_64bit_to_32bit_singlelogical(funit,mylog)
    implicit none
    !> File unit to read from (assumes file is already open)
    integer,intent(in) :: funit
    !> Output 32 bit integer
    logical,intent(inout) :: mylog
    logical(8) :: mylog_long

    read(funit) mylog_long
    mylog = mylog_long

  end subroutine read_64bit_to_32bit_singlelogical

  !> Read 32 bit logical from file and convert to 64 bit logical
  !> \author Patrick Ettenhuber (adapted from Kasper)
  !> \date December 2012
  subroutine read_32bit_to_64bit_singlelogical(funit,mylog)
    implicit none
    !> File unit to read from (assumes file is already open)
    integer,intent(in) :: funit
    !> Output 32 bit integer
    logical,intent(inout) :: mylog
    logical(kind=4) :: mylog_long

    read(funit) mylog_long
    mylog = mylog_long

  end subroutine read_32bit_to_64bit_singlelogical



  !> Read 64 bit integer vector from file and convert to 32 bit integer vector
  !> \author Kasper Kristensen
  !> \date December 2012
  subroutine read_64bit_to_32bit_vectorinteger(funit,n,myint)
    implicit none
    !> File unit to read from (assumes file is already open)
    integer,intent(in) :: funit
    !> Length of integer vector
    integer,intent(in) :: n
    !> Output 32 bit integer vector
    integer,intent(inout),dimension(n) :: myint
    integer(8),pointer :: myint_long(:)
    integer :: i

    call mem_alloc(myint_long,n)

    read(funit) myint_long
    do i=1,n
       myint(i) = int(myint_long(i))
    end do

    call mem_dealloc(myint_long)

  end subroutine read_64bit_to_32bit_vectorinteger

  !> Read 32 bit integer vector from file and convert to 64 bit integer vector
  !> \author Patrick Ettenhuber
  !> \date December 2012
  subroutine read_32bit_to_64bit_vectorinteger(funit,n,myint)
    implicit none
    !> File unit to read from (assumes file is already open)
    integer,intent(in) :: funit
    !> Length of integer vector
    integer,intent(in) :: n
    !> Output 32 bit integer vector
    integer,intent(inout),dimension(n) :: myint
    integer(kind=4),pointer :: myint_long(:)
    integer :: i

    call mem_alloc(myint_long,n)

    read(funit) myint_long
    do i=1,n
       myint(i) = int(myint_long(i))
    end do

    call mem_dealloc(myint_long)

  end subroutine read_32bit_to_64bit_vectorinteger


  !> Read 64 bit logical vector from file and convert to 32 bit logical vector
  !> \author Kasper Kristensen
  !> \date December 2012
  subroutine read_64bit_to_32bit_vectorlogical(funit,n,mylog)
    implicit none
    !> File unit to read from (assumes file is already open)
    integer,intent(in) :: funit
    !> Length of logical vector
    integer,intent(in) :: n
    !> Output 32 bit logical vector
    logical,intent(inout),dimension(n) :: mylog
    logical(8),pointer :: mylog_long(:)
    integer :: i

    allocate(mylog_long(n))

    read(funit) mylog_long
    do i=1,n
       mylog(i) = mylog_long(i)
    end do

    deallocate(mylog_long)

  end subroutine read_64bit_to_32bit_vectorlogical

  !> Read 32 bit logical vector from file and convert to 64 bit logical vector
  !> \author Kasper Kristensen
  !> \date December 2012
  subroutine read_32bit_to_64bit_vectorlogical(funit,n,mylog)
    implicit none
    !> File unit to read from (assumes file is already open)
    integer,intent(in) :: funit
    !> Length of logical vector
    integer,intent(in) :: n
    !> Output 32 bit logical vector
    logical,intent(inout),dimension(n) :: mylog
    logical(kind=4),pointer :: mylog_long(:)
    integer :: i

    allocate(mylog_long(n))

    read(funit) mylog_long
    do i=1,n
       mylog(i) = mylog_long(i)
    end do

    deallocate(mylog_long)

  end subroutine read_32bit_to_64bit_vectorlogical



  !> \brief Initialize t1 information in fragment structure and set
  !> t1 amplitude elements in MyFragment%t1 equal to input t1.
  !> It is assumed that amplitudes have dimension (virtual AOS, occupied AOS).
  !> \author Kasper Kristensen
  !> \date January 2012
  subroutine save_fragment_t1_AOSAOSamplitudes(MyFragment,t1)
    implicit none
    !> Fragment info (only t1-related information will be modified here)
    type(decfrag), intent(inout) :: MyFragment
    !> Singles amplitudes to be stored (stored as virtual,occupied)
    real(realk),intent(in) :: t1(:,:)
    integer :: nocc,nvirt,i,a,ix,ax

    ! Init dimensions
    nocc = MyFragment%noccAOS   ! occupied AOS dimension
    nvirt = MyFragment%nvirtAOS   ! virtual AOS dimension

    ! Sanity check
    !if( (nvirt/=t1%dims(1)) .or. (nocc/=t1%dims(2)) ) then
    !   write(DECinfo%output,*) 'Fragment virt,occ', nvirt,nocc
    !   write(DECinfo%output,*) 't1 input virt,occ', t1%dims
    !   call lsquit('save_fragment_t1_AOSAOSamplitudes &
    !        & AOS dimension mismatch!',DECinfo%output)
    !end if

    ! Free t1 stuff (in case old ampltiudes are already stored)
    call free_fragment_t1(MyFragment)

    ! Init t1 stuff according to fragment information
    MyFragment%t1_stored = .true.   ! t1 amplitudes will be stored
    MyFragment%t1dims(1) = nvirt
    MyFragment%t1dims(2) = nocc
    call mem_alloc(MyFragment%t1_occidx,nocc)
    call mem_alloc(MyFragment%t1_virtidx,nvirt)
    call mem_alloc(MyFragment%t1,nvirt,nocc)
    MyFragment%t1_occidx = MyFragment%occAOSidx ! occupied AOS indices
    MyFragment%t1_virtidx = MyFragment%virtAOSidx ! virtual AOS indices


    ! Save amplitudes and indices
    do i=1,nocc
       do a=1,nvirt
          MyFragment%t1(a,i) = t1(a,i)
       end do
    end do


  end subroutine save_fragment_t1_AOSAOSamplitudes



  !> \brief Construct logical array telling which atomic pair combinations should be
  !> included when calculating pair interaction energies and other pair interaction
  !> properties - to avoid double counting.
  !> It also works for standard fragments but in that case the output is trivial.
  !> \author Kasper Kristensen
  !> \date September 2011
  subroutine which_pairs(Fragment1,Fragment2,natoms,dopair)

    implicit none
    ! Fragment 1 in pair
    type(decfrag),intent(inout) :: Fragment1
    ! Fragment 2 in pair
    type(decfrag),intent(inout) :: Fragment2
    integer, intent(in) :: natoms
    logical,dimension(natoms,natoms),intent(inout) :: dopair
    integer :: a,b,ax,bx

    ! dopair tells which (atom1,atom2) combinations to include when calculating
    ! pair interaction properties.
    !
    ! Example:
    ! ********
    ! fragment 1 has been constructed from standard atoms 3 and 5
    ! fragment 2 has been constructed from standard atoms 7 and 8.
    ! In this case we should include all pair combinaitons of 3 and 5
    ! with 7 and 8 - and of course not the diagonal elements, since that would
    ! introduce double countings.
    ! Thus, in this example we would have:
    ! dopair(3,7) = dopair(7,3) = .true.
    ! dopair(3,8) = dopair(8,3) = .true.
    ! dopair(5,7) = dopair(7,5) = .true.
    ! dopair(5,8) = dopair(8,5) = .true.
    ! dopair(i,j) = .false. for all other pairs, including the diagonal elements.


    ! Set which atoms to consider for pair
    dopair=.false.
    do a=1,fragment1%nEOSatoms  ! Loop over atoms in fragment 1 
       do b=1,fragment2%nEOSatoms ! Loop over atoms in fragment 2

          ax=fragment1%EOSatoms(a) ! Atom index for atom in fragment 1
          bx=fragment2%EOSatoms(b) ! Atom index for atom in fragment 2
          dopair(ax,bx)=.true.
          dopair(bx,ax)=.true.

          ! Sanity check
          if(ax==bx) then
             call lsquit('which_pairs: &
                  & Something wrong with atomic indices in pair',DECinfo%output)
          end if

       end do
    end do



  end subroutine which_pairs



  !> \brief Construct logical array telling which occupied orbital pair combinations should be
  !> included when calculating pair interaction energies and other pair interaction
  !> properties - to avoid double counting.
  !> \author Kasper Kristensen
  !> \date September 2011
  subroutine which_pairs_occ(Fragment1,Fragment2,PairFragment,dopair)

    implicit none
    ! Fragment 1 in pair
    type(decfrag),intent(in) :: Fragment1
    ! Fragment 2 in pair
    type(decfrag),intent(in) :: Fragment2
    !> Pair fragment
    type(decfrag),intent(in) :: PairFragment
    !> Do pair or not - dimension: (noccEOS,noccEOS) for PAIR
    logical,dimension(PairFragment%noccEOS,PairFragment%noccEOS),&
         & intent(inout) :: dopair
    integer :: a,b,ax,bx,p1,p2,i,j

    ! dopair tells which (orb1,orb2) combinations to include when calculating
    ! pair interaction properties. This is best explained by an example.
    !
    ! Example:
    ! ********
    ! Pair fragment contains 5 occupied EOS orbitals with the following
    ! indices in the full molecular list of orbitals:
    ! 6,8,10,12,15
    ! These indices are saved in PairFragment%occEOSidx where
    ! they are simply ordered according to size:
    !
    ! Index in PairFragment%occEOSidx           Full index     (*)
    !               1                                   6
    !               2                                   8
    !               3                                   10
    !               4                                   12
    !               5                                   15
    !
    ! Fragment 1 contains occupied EOS indices: 6,10
    ! These indices are saved in Fragment1%occEOSidx
    !
    ! Fragment 2 contains occupied EOS indices: 8,12,15
    ! These indices are saved in Fragment2%occEOSidx
    !
    ! The "pair interaction indices" are then the combinations of
    ! indices from fragment 1 and fragment 2:
    ! (6,8)  (6,12)  (6,15)  (10,8)  (10,12)  (10,15)
    !
    ! Seen from the perspective of PairFragment%occEOSidx (*), we
    ! have following "pair interaction indices":
    ! (1,2)  (1,4)  (1,5)  (3,2)  (3,4)  (3,5)
    !
    ! Thus, the dopair output (dimension: 5,5) will be:
    ! dopair(1,2) = dopair(2,1) = .true.
    ! dopair(1,4) = dopair(4,1) = .true.
    ! dopair(1,5) = dopair(5,1) = .true.
    ! dopair(3,2) = dopair(2,3) = .true.
    ! dopair(3,4) = dopair(4,3) = .true.
    ! dopair(3,5) = dopair(5,3) = .true.
    ! dopair(i,j) = .false. for all other pairs, including the diagonal elements.


    ! Set which atoms to consider for pair
    dopair=.false.

    ! Skip this when only virtual part. is requested
    DoCheckOcc: if(.not. DECinfo%onlyvirtpart) then
       do a=1,fragment1%noccEOS   ! Occupied EOS for fragment 1
          do b=1,fragment2%noccEOS ! Occupied EOS for fragment 2
             
             ax=fragment1%occEOSidx(a)  ! index in full list orbitals
             bx=fragment2%occEOSidx(b)  ! index in full list orbitals
             p1=0
             p2=0
             
             ! Index for fragment 1 in pair fragment list
             do i=1,PairFragment%noccEOS
                
                ! "ax" index in PairFragment%occEOSidx list
                if(PairFragment%occEOSidx(i) == ax) p1 = i
                
                ! "bx" index in PairFragment%occEOSidx list
                if(PairFragment%occEOSidx(i) == bx) p2 = i
                
             end do
             
             ! Sanity check
             if(p1==p2 .or. p1==0 .or. p2==0 ) then
                write(DECinfo%output,'(1X,a,4i6)') 'ax,bx,p1,p2', ax,bx,p1,p2
                call lsquit('which_pairs_occ: &
                     & Something wrong with indices in pair',DECinfo%output)
             end if
             
             
             ! Pair interaction for (p1,p2) index pair
             dopair(p1,p2)=.true.
             dopair(p2,p1)=.true.
             
          end do
       end do
    endif DoCheckOcc

  end subroutine which_pairs_occ



  !> \brief Construct logical array telling which virtupied orbital pair combinations should be
  !> included when calculating pair interaction energies and other pair interaction
  !> properties - to avoid double counting.
  !> \author Kasper Kristensen
  !> \date September 2011
  subroutine which_pairs_virt(Fragment1,Fragment2,PairFragment,dopair)

    implicit none
    ! Fragment 1 in pair
    type(decfrag),intent(in) :: Fragment1
    ! Fragment 2 in pair
    type(decfrag),intent(in) :: Fragment2
    !> Pair fragment
    type(decfrag),intent(in) :: PairFragment
    !> Do pair or not - dimension: (nvirtEOS,nvirtEOS) for PAIR
    logical,dimension(PairFragment%nvirtEOS,PairFragment%nvirtEOS),&
         & intent(inout) :: dopair
    integer :: a,b,ax,bx,p1,p2,i

    ! This is the same as which_pairs_occ, but for the virtupied space.
    ! See example in which_pairs_occ with "occ" replaced by "virt".

    ! Set which atoms to consider for pair
    dopair=.false.

    ! Skip this when only occupied part. is requested
    DoCheck: if(.not. DECinfo%onlyoccpart) then
       do a=1,fragment1%nvirtEOS   ! Unoccupied EOS for fragment 1
          do b=1,fragment2%nvirtEOS ! Unoccupied EOS for fragment 2

             ax=fragment1%virtEOSidx(a)  ! index in full list orbitals
             bx=fragment2%virtEOSidx(b)  ! index in full list orbitals
             p1=0
             p2=0

             ! Index for fragment 1 in pair fragment list
             do i=1,PairFragment%nvirtEOS

                ! "ax" index in PairFragment%virtEOSidx list
                if(PairFragment%virtEOSidx(i) == ax) p1 = i

                ! "bx" index in PairFragment%virtEOSidx list
                if(PairFragment%virtEOSidx(i) == bx) p2 = i

             end do

             ! Sanity check
             if(p1==p2 .or. p1==0 .or. p2==0 ) then
                write(DECinfo%output,'(1X,a,4i6)') 'ax,bx,p1,p2', ax,bx,p1,p2
                call lsquit('which_pairs_virt: &
                     & Something wrong with indices in pair',DECinfo%output)
             end if


             ! Pair interaction for (p1,p2) index pair
             dopair(p1,p2)=.true.
             dopair(p2,p1)=.true.

          end do
       end do

    end if DoCheck

  end subroutine which_pairs_virt

  !> \brief Construct logical array telling which virtupied/occupied combinations should be
  !> included when calculating pair interaction energies and other pair interaction
  !> properties - to avoid double counting.
  !> \author Dmytro Bykov
  !> \date October 2014
  subroutine which_pairs_occ_virt(Fragment1,Fragment2,PairFragment,dopair)

    implicit none
    ! Fragment 1 in pair
    type(decfrag),intent(in) :: Fragment1
    ! Fragment 2 in pair
    type(decfrag),intent(in) :: Fragment2
    !> Pair fragment
    type(decfrag),intent(inout) :: PairFragment
    !> Do pair or not - dimension: (nvirtEOS,nvirtEOS) for PAIR
    logical,dimension(PairFragment%nvirtEOS,PairFragment%noccEOS),&
         & intent(inout) :: dopair
    integer :: a,b,ax,bx,p1,p2,i,ix,pox,pux

    ! This is the same as which_pairs_virt, but for the case where one needs 
    ! virtupied/occupied combination.

    ! Set which atoms to consider for pair
    dopair=.false.

    ! occ fragment1 - virt fragment2 part:
    do a=1,fragment1%nvirtEOS   ! Unoccupied EOS for fragment 1
      do i=1,fragment2%noccEOS   ! Occupied   EOS for fragment 2

        ax=fragment1%virtEOSidx(a)  ! index in full list orbitals
        ix=fragment2%occEOSidx(i)    ! index in full list orbitals
        p1=0                     ! pair virt full list index
        p2=0                     ! pair occ   full list index

        ! loop to find pair indices  
        do pux=1,PairFragment%nvirtEOS
          do pox=1,PairFragment%noccEOS

            if(PairFragment%virtEOSidx(pux) == ax) p1 = pux
            if(PairFragment%occEOSidx(pox)   == ix) p2 = pox

          end do
        end do

        ! Sanity check
!        if(p1==p2 .or. p1==0 .or. p2==0 ) then
!          write(DECinfo%output,'(1X,a,4i6)') 'ax,ix,p1,p2', ax,ix,p1,p2
!          call lsquit('which_pairs_occ_virt: &
!                     & Something wrong with indices in pair',DECinfo%output)
!        end if

        ! Pair interaction for (p1,p2) index pair
        dopair(p1,p2)=.true.

      end do
    end do

    ! occ fragment2 - virt fragment1 part:
    do a=1,fragment2%nvirtEOS   ! Unoccupied EOS for fragment 1
      do i=1,fragment1%noccEOS   ! Occupied   EOS for fragment 2

        ax=fragment2%virtEOSidx(a)  ! index in full list orbitals
        ix=fragment1%occEOSidx(i)    ! index in full list orbitals
        p1=0                     ! pair virt full list index
        p2=0                     ! pair occ   full list index

        ! loop to find pair indices  
        do pux=1,PairFragment%nvirtEOS
          do pox=1,PairFragment%noccEOS
            
            if(PairFragment%virtEOSidx(pux) == ax) p1 = pux
            if(PairFragment%occEOSidx(pox)   == ix) p2 = pox
            
          end do
        end do
        
        ! Sanity check
!        if(p1==p2 .or. p1==0 .or. p2==0 ) then
!          write(DECinfo%output,'(1X,a,4i6)') 'ax,ix,p1,p2', ax,ix,p1,p2
!          call lsquit('which_pairs_occ_virt: &
!                     & Something wrong with indices in pair',DECinfo%output)
!        end if
        
        ! Pair interaction for (p1,p2) index pair
        dopair(p1,p2)=.true.
        
      end do
    end do

  end subroutine which_pairs_occ_virt


  !> Write fragment job list to file.
  !> Also write pair cutoff distance in case that was changed during the calculation.
  !> \author Kasper Kristensen
  !> \date November 2012
  subroutine write_fragment_joblist_to_file(jobs,funit)
    implicit none
    !> Job list of fragments
    type(joblist),intent(in) :: jobs
    !> File unit number to write to (of course assumes that file is open)
    integer,intent(in) :: funit
    logical(8) :: jobsdone64(jobs%njobs),dofragopt64(jobs%njobs),esti64(jobs%njobs)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! IMPORTANT: ALWAYS WRITE AND READ INTEGERS AND LOGICALS WITH 64BIT!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    jobsdone64  = jobs%jobsdone
    dofragopt64 = jobs%dofragopt
    esti64      = jobs%esti

    write(funit) int(jobs%njobs,kind=8)
    write(funit) int(jobs%atom1,kind=8)
    write(funit) int(jobs%atom2,kind=8)
    write(funit) int(jobs%jobsize,kind=8)
    write(funit) jobsdone64
    write(funit) dofragopt64
    write(funit) esti64

    ! MPI fragment statistics
    write(funit) int(jobs%nslaves,kind=8)
    write(funit) int(jobs%nocc,kind=8)
    write(funit) int(jobs%nvirt,kind=8)
    write(funit) int(jobs%nbasis,kind=8)
    write(funit) int(jobs%ntasks,kind=8)
    write(funit) jobs%flops
    write(funit) jobs%LMtime
    write(funit) jobs%workt
    write(funit) jobs%commt
    write(funit) jobs%idlet
    write(funit) jobs%comm_gl_master_time
    write(funit) jobs%gpu_flops

  end subroutine write_fragment_joblist_to_file


  !> Read fragment job list from file assuming that joblist has already been initialized 
  !> with the proper dimensions.
  !> Also read pair cutoff distance.
  !> \author Kasper Kristensen
  !> \date November 2012
  subroutine read_fragment_joblist_from_file(jobs,funit)
    implicit none
    !> Job list of fragments
    type(joblist),intent(inout) :: jobs
    !> File unit number to read from (of course assumes that file is open)
    integer,intent(in) :: funit
    integer :: njobs


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! IMPORTANT: ALWAYS WRITE AND READ INTEGERS AND LOGICALS WITH 64BIT!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    call read_64bit_to_int(funit,njobs)
    if(njobs/=jobs%njobs) then
       print *, 'Number of jobs in job list   : ', jobs%njobs
       print *, 'Number of jobs read from file: ', njobs
       call lsquit('read_fragment_joblist_from_file1: Error in number of jobs!',-1)
    end if
    call read_64bit_to_int(funit,njobs,jobs%atom1)
    call read_64bit_to_int(funit,njobs,jobs%atom2)
    call read_64bit_to_int(funit,njobs,jobs%jobsize)
    call read_64bit_to_int(funit,njobs,jobs%jobsdone)
    call read_64bit_to_int(funit,njobs,jobs%dofragopt)
    call read_64bit_to_int(funit,njobs,jobs%esti)
    call read_64bit_to_int(funit,njobs,jobs%nslaves)
    call read_64bit_to_int(funit,njobs,jobs%nocc)
    call read_64bit_to_int(funit,njobs,jobs%nvirt)
    call read_64bit_to_int(funit,njobs,jobs%nbasis)
    call read_64bit_to_int(funit,njobs,jobs%ntasks)

    read(funit) jobs%flops
    read(funit) jobs%LMtime
    read(funit) jobs%workt
    read(funit) jobs%commt
    read(funit) jobs%idlet
    read(funit) jobs%comm_gl_master_time
    read(funit) jobs%gpu_flops

    write(DECinfo%output,*)
    write(DECinfo%output,*) 'JOB LIST RESTART'
    write(DECinfo%output,'(1X,a,i10)') '-- total number of jobs : ', jobs%njobs
    write(DECinfo%output,'(1X,a,i10)') '-- number of jobs done  : ', count(jobs%jobsdone)
    write(DECinfo%output,'(1X,a,i10)') '-- number of jobs to do : ', jobs%njobs-count(jobs%jobsdone)
    write(DECinfo%output,*)


  end subroutine read_fragment_joblist_from_file



  !> \brief Initialize fragment job list
  !> \author Kasper Kristensen
  !> \date May 2012
  subroutine init_joblist(njobs,jobs)
    implicit none
    !> Number of jobs
    integer,intent(in) :: njobs
    !> Job list
    type(joblist),intent(inout) ::  jobs

    ! Number of jobs
    jobs%njobs = njobs

    ! Set all pointers to be of size njobs and equal to 0
    if (njobs>0) then
       call mem_alloc(jobs%atom1,njobs)
       call mem_alloc(jobs%atom2,njobs)
       call mem_alloc(jobs%jobsize,njobs)
       call mem_alloc(jobs%jobsdone,njobs)
       call mem_alloc(jobs%dofragopt,njobs)
       call mem_alloc(jobs%esti,njobs)
       jobs%atom1     = 0
       jobs%atom2     = 0
       jobs%jobsize   = 0
       jobs%jobsdone  = .false. ! no jobs are done
       jobs%dofragopt = .false. 
       jobs%esti      = .false.
        
       ! MPI fragment statistics
       call mem_alloc(jobs%nslaves,njobs)
       call mem_alloc(jobs%nocc,njobs)
       call mem_alloc(jobs%nvirt,njobs)
       call mem_alloc(jobs%nbasis,njobs)
       call mem_alloc(jobs%ntasks,njobs)
       call mem_alloc(jobs%flops,njobs)
       call mem_alloc(jobs%LMtime,njobs)
       call mem_alloc(jobs%commt,njobs)
       call mem_alloc(jobs%workt,njobs)
       call mem_alloc(jobs%idlet,njobs)
       call mem_alloc(jobs%gpu_flops,njobs)
       call mem_alloc(jobs%comm_gl_master_time,njobs)
       jobs%nslaves = 0
       jobs%nocc    = 0
       jobs%nvirt   = 0
       jobs%nbasis  = 0
       jobs%ntasks  = 0
       jobs%flops     = 0.0E0_realk
       jobs%LMtime    = 0.0E0_realk
       jobs%commt     = 0.0E0_realk
       jobs%workt     = 0.0E0_realk
       jobs%idlet     = 0.0E0_realk
       jobs%gpu_flops = 0.0E0_realk
       jobs%comm_gl_master_time = 0.0E0_realk
    end if

  end subroutine init_joblist



  !> \brief Free fragment job list
  !> \author Kasper Kristensen
  !> \date May 2012
  subroutine free_joblist(jobs)
    implicit none
    !> Job list
    type(joblist),intent(inout) ::  jobs

    if (jobs%njobs>0) then
    ! Deallocate pointers and nullify
    if(associated(jobs%atom1)) then
       call mem_dealloc(jobs%atom1)
       nullify(jobs%atom1)
    end if

    if(associated(jobs%atom2)) then
       call mem_dealloc(jobs%atom2)
       nullify(jobs%atom2)
    end if

    if(associated(jobs%jobsize)) then
       call mem_dealloc(jobs%jobsize)
       nullify(jobs%jobsize)
    end if

    if(associated(jobs%jobsdone)) then
       call mem_dealloc(jobs%jobsdone)
       nullify(jobs%jobsdone)
    end if

    if(associated(jobs%dofragopt)) then
       call mem_dealloc(jobs%dofragopt)
       nullify(jobs%dofragopt)
    end if

    if(associated(jobs%esti)) then
       call mem_dealloc(jobs%esti)
       nullify(jobs%esti)
    end if

    if(associated(jobs%nslaves)) then
       call mem_dealloc(jobs%nslaves)
       nullify(jobs%nslaves)
    end if

    if(associated(jobs%nocc)) then
       call mem_dealloc(jobs%nocc)
       nullify(jobs%nocc)
    end if

    if(associated(jobs%nvirt)) then
       call mem_dealloc(jobs%nvirt)
       nullify(jobs%nvirt)
    end if

    if(associated(jobs%nbasis)) then
       call mem_dealloc(jobs%nbasis)
       nullify(jobs%nbasis)
    end if

    if(associated(jobs%ntasks)) then
       call mem_dealloc(jobs%ntasks)
       nullify(jobs%ntasks)
    end if

    if(associated(jobs%flops)) then
       call mem_dealloc(jobs%flops)
       nullify(jobs%flops)
    end if

    if(associated(jobs%gpu_flops)) then
       call mem_dealloc(jobs%gpu_flops)
       nullify(jobs%gpu_flops)
    end if

    if(associated(jobs%LMtime)) then
       call mem_dealloc(jobs%LMtime)
       nullify(jobs%LMtime)
    end if

    if(associated(jobs%workt)) then
       call mem_dealloc(jobs%workt)
       nullify(jobs%workt)
    end if

    if(associated(jobs%commt)) then
       call mem_dealloc(jobs%commt)
       nullify(jobs%commt)
    end if

    if(associated(jobs%idlet)) then
       call mem_dealloc(jobs%idlet)
       nullify(jobs%idlet)
    end if
    if(associated(jobs%comm_gl_master_time)) then
       call mem_dealloc(jobs%comm_gl_master_time)
       nullify(jobs%comm_gl_master_time)
    end if
    end if

  end subroutine free_joblist


  !> Put job list info for a given job into bigger job list
  !> \author Kasper Kristensen
  !> \date November 2012
  subroutine put_job_into_joblist(singlejob,position,jobs)
    implicit none
    !> Job "list" containing a single job to be put into bigger job list
    type(joblist),intent(in) :: singlejob
    !> Position in big job list to overwrite with new job info for single job
    integer,intent(in) :: position
    !> Big job list
    type(joblist),intent(inout) :: jobs

    ! Sanity check 1: singlejob should really contain just a single job
    if(singlejob%njobs /= 1) then
       write(DECinfo%output,*) 'Number of jobs in single job list: ', singlejob%njobs
       call lsquit('put_job_into_joblist: singlejob does not contain ONE job!',-1)
    end if

    ! Sanity check 2: Position must not exceed job size of big job list
    if(position > jobs%njobs) then
       write(DECinfo%output,*) 'Position / # jobs ', position,jobs%njobs
       call lsquit('put_job_into_joblist: Input position exceed size of job list!',-1)
    end if

    ! Copy info from single job into big job list
    jobs%atom1(position)     = singlejob%atom1(1)
    jobs%atom2(position)     = singlejob%atom2(1)
    jobs%jobsize(position)   = singlejob%jobsize(1)
    jobs%jobsdone(position)  = singlejob%jobsdone(1)
    jobs%dofragopt(position) = singlejob%dofragopt(1)
    jobs%esti(position)      = singlejob%esti(1)
    jobs%nslaves(position)   = singlejob%nslaves(1)
    jobs%nocc(position)      = singlejob%nocc(1)
    jobs%nvirt(position)     = singlejob%nvirt(1)
    jobs%nbasis(position)    = singlejob%nbasis(1)
    jobs%ntasks(position)    = singlejob%ntasks(1)
    jobs%flops(position)     = singlejob%flops(1)
    jobs%LMtime(position)    = singlejob%LMtime(1)
    jobs%workt(position)     = singlejob%workt(1)
    jobs%commt(position)     = singlejob%commt(1)
    jobs%idlet(position)     = singlejob%idlet(1)
    jobs%gpu_flops(position) = singlejob%gpu_flops(1)
    jobs%comm_gl_master_time(position) = singlejob%comm_gl_master_time(1)

  end subroutine put_job_into_joblist


  !> \brief Estimate memory used in MP2 energy calculation
  !> \author Kasper Kristensen
  !> \date February 2011
  subroutine estimate_memory_for_mp2_energy(nthreads,O,V,A,B,intMEM,intStep,solMEM)

    implicit none
    !> Number of OMP threads
    integer, intent(in) :: nthreads
    !> Number of occupied orbitals (as a real)
    real(realk), intent(in) :: O
    !> Number of virtual orbitals (as a real)
    real(realk), intent(in) :: V
    !> Number of atomic orbitals (as a real)
    real(realk), intent(in) :: A
    !> Maximum batch dimension (as a real)
    real(realk), intent(in) :: B
    !> Maximum memory for integrals (in GB)
    real(realk),intent(inout) :: intMEM
    !> Step in integral routine where memory use is greatest
    integer, intent(inout) :: intStep
    !> Maximum memory for solver (in GB)
    real(realk),intent(inout) :: solMEM
    real(realk) :: tmp,GB

    ! Initialize stuff
    ! ****************
    GB = 1.000E9_realk ! 1 GB


    ! Memory for integrals (in GB)
    ! ****************************

    ! In different places of the get_VOVO_integrals_mem routine, different amounts of memory
    ! are allocated. Here we find the maximum.
    ! Roughly speaking, we can talk about 5 different memory comsumptions.
    intMEM = O*V*O*A + V*O*A*B + (A*A*B*B + V*A*B*B + V*O*B*B)*nthreads ! 1
    intStep=1

    tmp = O*V*O*A + 2E0_realk*V*O*A*B ! 2
    if(tmp > intMEM) then
       intMEM=tmp
       intStep=2
    end if

    tmp = O*V*O*A + V*O*A*B + O*V*O*B ! 3
    if(tmp > intMEM) then
       intMEM=tmp
       intStep=3
    end if

    tmp = 2E0_realk*O*V*O*A ! 4
    if(tmp > intMEM) then
       intMEM=tmp
       intStep=4
    end if

    tmp = O*V*O*A + V*O*V*O !5
    if(tmp > intMEM) then
       intMEM=tmp
       intStep=5
    end if

    ! Multiply intMEM by realk (8) and divide by GB to get memory in GB
    intMEM = realk*intMEM/GB


    ! Memory for solver (in GB)
    ! *************************
    ! Maximum memory allocated in solver is three arrays of dimensions (V,O,V,O)
    solMEM = 3E0_realk*(realk*V*O*V*O)/GB


  end subroutine estimate_memory_for_mp2_energy



  !> \brief Simple basis transform: matB = C^T matA C
  !> using simple Fortran arrays.
  !> \author Kasper Kristensen
  !> \date March 2013
  subroutine dec_simple_basis_transform1(nA,nB,C,matA,matB)
    implicit none
    !> A dimension
    integer,intent(in) :: nA
    !> B dimension
    integer,intent(in) :: nB
    !> B coefficient matrix
    real(realk),dimension(nA,nB),intent(in) :: C
    !> Matrix in A basis
    real(realk),dimension(nA,nA),intent(in) :: matA
    !> Matrix in B basis
    real(realk),dimension(nB,nB),intent(inout) :: matB
    real(realk),pointer :: tmp(:,:)

    ! tmp = C^T matA
    call mem_alloc(tmp,nB,nA)
    call dec_simple_dgemm(nB,nA,nA,C,matA,tmp,'t','n')

    ! matB = tmp C
    call dec_simple_dgemm(nB,nA,nB,tmp,C,matB,'n','n')
    call mem_dealloc(tmp)

  end subroutine dec_simple_basis_transform1


  !> \brief Simple transform: matB = C matA C^T
  !> using simple Fortran arrays.
  !> \author Kasper Kristensen
  !> \date March 2013
  subroutine dec_simple_basis_transform2(nB,nA,C,matA,matB)
    implicit none
    !> B dimension
    integer,intent(in) :: nB
    !> A dimension
    integer,intent(in) :: nA
    !> A coefficient matrix
    real(realk),dimension(nB,nA),intent(in) :: C
    !> Matrix in A basis
    real(realk),dimension(nA,nA),intent(in) :: matA
    !> Matrix in B basis
    real(realk),dimension(nB,nB),intent(inout) :: matB
    real(realk),pointer :: tmp(:,:)

    ! tmp = C matA
    call mem_alloc(tmp,nB,nA)
    call dec_simple_dgemm(nB,nA,nA,C,matA,tmp,'n','n')

    ! matB = tmp C^T
    call dec_simple_dgemm(nB,nA,nB,tmp,C,matB,'n','t')
    call mem_dealloc(tmp)

  end subroutine dec_simple_basis_transform2



  !> \brief Transformation using different coefficient matrices
  !> (e.g. C1 could be occupied MOs and C2 could be virtual MOs):
  !> matB = C1^T matA C2
  !> using simple Fortran arrays.
  !> \author Kasper Kristensen
  !> \date March 2013
  subroutine dec_diff_basis_transform1(nA,nB1,nB2,C1,C2,matA,matB)
    implicit none
    !> A dimension
    integer,intent(in) :: nA
    !> B dimension 1
    integer,intent(in) :: nB1
    !> B dimension 2
    integer,intent(in) :: nB2
    !> A coefficient matrix 1
    real(realk),dimension(nA,nB1),intent(in) :: C1
    !> A coefficient matrix 2
    real(realk),dimension(nA,nB2),intent(in) :: C2
    !> Matrix in A basis
    real(realk),dimension(nA,nA),intent(in) :: matA
    !> Matrix in B basis
    real(realk),dimension(nB1,nB2),intent(inout) :: matB
    real(realk),pointer :: tmp(:,:)

    ! tmp = C1^T matA
    call mem_alloc(tmp,nB1,nA)
    call dec_simple_dgemm(nB1,nA,nA,C1,matA,tmp,'t','n')

    ! matB = tmp C2
    call dec_simple_dgemm(nB1,nA,nB2,tmp,C2,matB,'n','n')
    call mem_dealloc(tmp)

  end subroutine dec_diff_basis_transform1



  !> \brief Transformation using different coefficient matrices
  !> (e.g. C1 could be occupied MOs and C2 could be virtual MOs):
  !> matB = C1 matA C2^T
  !> using simple Fortran arrays.
  !> \author Kasper Kristensen
  !> \date March 2013
  subroutine dec_diff_basis_transform2(nB,nA1,nA2,C1,C2,matA,matB)
    implicit none
    !> B dimension
    integer,intent(in) :: nB
    !> A dimension 1
    integer,intent(in) :: nA1
    !> A dimension 2
    integer,intent(in) :: nA2
    !> A coefficient matrix 1
    real(realk),dimension(nB,nA1),intent(in) :: C1
    !> A coefficient matrix 2
    real(realk),dimension(nB,nA1),intent(in) :: C2
    !> Matrix in A basis
    real(realk),dimension(nA1,nA2),intent(in) :: matA
    !> Matrix in B basis
    real(realk),dimension(nB,nB),intent(inout) :: matB
    real(realk),pointer :: tmp(:,:)

    ! tmp = C1 matA
    call mem_alloc(tmp,nB,nA2)
    call dec_simple_dgemm(nB,nA1,nA2,C1,matA,tmp,'n','n')

    ! matB = tmp C2^T
    call dec_simple_dgemm(nB,nA2,nB,tmp,C2,matB,'n','t')
    call mem_dealloc(tmp)

  end subroutine dec_diff_basis_transform2


  !> Add DEC energies: E = sum_P E_P  +  sum_{P>Q} dE_PQ
  !> taking into account that not all atoms have orbitals assigned.
  !> \author Kasper Kristensen
  !> \date March 2013
  subroutine add_dec_energies(natoms,FragEnergies,orbitals_assigned,E)

    implicit none
    !> Number of atoms in molecule
    integer,intent(in) :: natoms
    !> Fragment energies (E_P on diagonal, dE_PQ on off-diagonal)
    real(realk),dimension(natoms,natoms),intent(in) :: FragEnergies
    !> Which atoms have orbitals assigned?
    logical,dimension(natoms) :: orbitals_assigned
    !> Total energy E = sum_P E_P  +  sum_{P>Q} dE_PQ 
    real(realk),intent(inout) :: E
    integer :: P,Q

    E = 0.0_realk
    do P=1,natoms
       if(orbitals_assigned(P)) then
          do Q=1,P
             if(orbitals_assigned(Q)) then
                E = E + FragEnergies(P,Q)
             end if
          end do
       end if
    end do

  end subroutine add_dec_energies

  !> Add estimated DEC energies for pairs which are skipped from the calculation
  !> to estimate the associated error.
  !> \author Kasper Kristensen
  !> \date October 2013
  subroutine estimate_energy_of_skipped_pairs(natoms,FragEnergies,orbitals_assigned,MyMolecule,E)

    implicit none
    !> Number of atoms in molecule
    integer,intent(in) :: natoms
    !> Estimated fragment energies (E_P on diagonal, dE_PQ on off-diagonal)
    real(realk),dimension(natoms,natoms),intent(in) :: FragEnergies
    !> Which atoms have orbitals assigned?
    logical,dimension(natoms) :: orbitals_assigned
    !> Full molecule structure, where MyMolecule%ccmodel defines the model to be used for
    !> each pair. If the model is MODEL_NONE (see dec_typedef.F90), then the pair
    !> will be skipped in the DEC calculation.
    type(fullmolecule),intent(in) :: MyMolecule
    !> Estimated energy contribution from skipped pairs
    real(realk),intent(inout) :: E
    integer :: P,Q

    E = 0.0_realk
    do P=1,natoms
       if(orbitals_assigned(P)) then
          do Q=1,P-1
             if(orbitals_assigned(Q) .and. MyMolecule%ccmodel(P,Q)==MODEL_NONE ) then
                E = E + FragEnergies(P,Q)
             end if
          end do
       end if
    end do

  end subroutine estimate_energy_of_skipped_pairs


  !> \brief Project orbitals onto MO space defined by input (see details inside subroutine).
  !> \author Kasper Kristensen
  !> \date April 2013
  subroutine project_onto_MO_space(nMOC,nMOZ,nAO,C,S,Z)
    implicit none

    !> MO dimension for C coefficients (see below)
    integer,intent(in) :: nMOC
    !> MO dimension for Z coefficients (see below)
    integer,intent(in) :: nMOZ
    !> AO dimension 
    integer,intent(in) :: nAO
    !> MO coefficients for orbitals {phi} to use in projector (see details below)
    real(realk),intent(in),dimension(nAO,nMOC) :: C
    !> AO overlap matrix
    real(realk),intent(in),dimension(nAO,nAO) :: S
    !> MO coefficients for orbitals {psi} to be projected (see details below)
    real(realk),intent(inout),dimension(nAO,nMOZ) :: Z
    real(realk),pointer :: tmp(:,:),tmp2(:,:),M(:,:),Minv(:,:)


    ! The orbitals to be projected |psi_r> are written in terms of AOs {chi}:
    !
    ! |psi_r> = sum_{alpha} Z_{alpha r} |chi_alpha>
    ! 
    ! The projector is: 
    ! 
    ! P = sum_{pq} |phi_p> (M^-1)_pq <phi_q|
    ! 
    ! where the MOs to project against {phi} are given in terms of AOs as:
    !
    ! |phi_p> = sum_{mu} C_{mu p} |chi_mu> 
    !
    ! and M is the MO overlap matrix:
    !
    ! M_pq = <phi_p | phi_q>.
    ! 
    ! Effectively we do the projection:
    ! 
    ! |psi_r> --> P |psi_r> = ( C M^-1 C^T S Z )_{mu r} |chi_mu>
    !
    ! where S is the AO overlap matrix: S_{mu nu} = <chi_mu | chi_nu>
    !
    ! Thus, the task of this subroutine is to change the input Z to (C M^-1 C^T S Z).


    ! Get inverse overlap for phi orbitals: M^-1 = (C^T S C)^-1
    ! *********************************************************
    call mem_alloc(M,nMOC,nMOC)
    call dec_simple_basis_transform1(nAO,nMOC,C,S,M)
    ! Minv = M^-1
    call mem_alloc(Minv,nMOC,nMOC)

    call invert_matrix(M,Minv,nMOC)
    call mem_dealloc(M)

    ! tmp = S Z
    call mem_alloc(tmp,nAO,nMOZ)
    call dec_simple_dgemm(nAO,nAO,nMOZ,S,Z,tmp,'n','n')

    ! tmp2 = C^T S Z
    call mem_alloc(tmp2,nMOC,nMOZ)
    call dec_simple_dgemm(nMOC,nAO,nMOZ,C,tmp,tmp2,'t','n')
    call mem_dealloc(tmp)

    ! tmp = M^-1 C^T S Z
    call mem_alloc(tmp,nMOC,nMOZ)
    call dec_simple_dgemm(nMOC,nMOC,nMOZ,Minv,tmp2,tmp,'n','n')
    call mem_dealloc(tmp2)

    ! Z --> C M^-1 C^T S Z
    call dec_simple_dgemm(nAO,nMOC,nMOZ,C,tmp,Z,'n','n')
    call mem_dealloc(tmp)
    call mem_dealloc(Minv)

  end subroutine project_onto_MO_space



  !> \brief Project {phi} MO space defined by input out of {psi} MO space (details inside subroutine).
  !> \author Kasper Kristensen
  !> \date April 2013
  subroutine project_out_MO_space(nMOC,nMOZ,nAO,C,S,Z)
    implicit none

    !> MO dimension for C coefficients (see below)
    integer,intent(in) :: nMOC
    !> MO dimension for Z coefficients (see below)
    integer,intent(in) :: nMOZ
    !> AO dimension 
    integer,intent(in) :: nAO
    !> MO coefficients for orbitals {phi} to project out (see details below)
    real(realk),intent(in),dimension(nAO,nMOC) :: C
    !> AO overlap matrix
    real(realk),intent(in),dimension(nAO,nAO) :: S
    !> MO coefficients for orbitals {psi} to be projected (see details below)
    real(realk),intent(inout),dimension(nAO,nMOZ) :: Z
    real(realk),pointer :: PZ(:,:)

    integer :: i,j

    ! The orbitals to be projected |psi_r> are written in terms of AOs {chi}:
    !
    ! |psi_r> = sum_{alpha} Z_{alpha r} |chi_alpha>
    ! 
    ! The projector is: 
    ! 
    ! P = 1 - sum_{pq} |phi_p> (M^-1)_pq <phi_q|
    ! 
    ! where the MOs to projected out {phi} are given in terms of AOs as:
    !
    ! |phi_p> = sum_{mu} C_{mu p} |chi_mu> 
    !
    ! and M is the MO overlap matrix:
    !
    ! M_pq = <phi_p | phi_q>.
    ! 
    ! Effectively we do the projection:
    ! 
    ! |psi_r> --> (1 - P) |psi_r> = (Z  -  C M^-1 C^T S Z )_{mu r} |chi_mu>
    !
    ! where S is the AO overlap matrix: S_{mu nu} = <chi_mu | chi_nu>
    !
    ! Thus, the task of this subroutine is to change the input Z to:
    ! Z - PZ = (Z  -  C M^-1 C^T S Z).

    ! Copy Z and calculate projection on Z: PZ = (C M^-1 C^T S) Z
    call mem_alloc(PZ,nAO,nMOZ)
    PZ = Z 
    call project_onto_MO_space(nMOC,nMOZ,nAO,C,S,PZ)

    ! Set output Z as: Z - PZ
    do i=1,nAO
       do j=1,nMOZ
          Z(i,j) = Z(i,j) - PZ(i,j)
       end do
    end do
    call mem_dealloc(PZ)

  end subroutine project_out_MO_space


  !> \brief Orthogonalize MOs by (i) setting up MO overlap matrix,
  !> (ii) diagonalizing MO overlap matrix, (iii) writing new orthogonalized MOs
  !> in terms of eigenvalues and eigenvalues of MO overlap matrix.
  !> \author Kasper Kristensen
  !> \date April 2013
  subroutine orthogonalize_MOs(nMO,nAO,S,C)
    implicit none

    !> MO dimension 
    integer,intent(in) :: nMO
    !> AO dimension 
    integer,intent(in) :: nAO
    !> AO overlap matrix
    real(realk),intent(in),dimension(nAO,nAO) :: S
    !> MO coefficients to be orthogonalized
    real(realk),intent(inout),dimension(nAO,nMO) :: C
    real(realk),pointer :: M(:,:), lambda(:), T(:,:), CT(:,:)
    integer :: mu,p
    real(realk) :: lambdascale

    ! The MOs are orthogonalized as follows:
    !
    ! (i) Setup MO overlap matrix:  M = C^T S C
    !
    ! (ii) Diagonalize M:  lambda = T^T M T    (lambda is diagonal matrix, T is unitary)
    !
    ! (iii) Labelling the input MOs as {phi} and the AOs as {chi}, we may write {phi} as
    !
    ! phi_p = sum_{mu} C_{mu p} |chi_mu>
    ! 
    ! and the orthogonalized MOs {psi} are given as:
    !
    ! psi_p = sum_mu lambda_p^{-1/2} (C T)_{mu p} |chi_mu>
    !
    ! Thus, the task of this subroutine is to change the input C_{mu p} 
    ! to lambda_p^{-1/2} (C T)_{mu p}.


    ! (i) M = C^T S C
    ! ***************
    call mem_alloc(M,nMO,nMO)
    call dec_simple_basis_transform1(nAO,nMO,C,S,M)


    ! (ii) Diagonalize MO overlap matrix: lambda = T^T M T
    ! ****************************************************
    call mem_alloc(lambda,nMO)
    call mem_alloc(T,nMO,nMO)
    call solve_eigenvalue_problem_unitoverlap(nMO,M,lambda,T)

    ! All lambda's should be positive. However, to be complete sure we do not do something
    ! dirty with square roots and negative numbers, we take the absolute value...
    lambda=abs(lambda)


    ! (iii) Get final orthogonalized MOs
    ! **********************************

    ! C T
    call mem_alloc(CT,nAO,nMO)
    call dec_simple_dgemm(nAO,nMO,nMO,C,T,CT,'n','n')

    ! C_{mu p} --> lambda_p^{-1/2} (C T)_{mu p}
    do p=1,nMO
       lambdascale = 1.0_realk/sqrt(lambda(p)) ! Lambda scaling factor: lambda_p^{-1/2}
       do mu=1,nAO
          C(mu,p) = CT(mu,p)*lambdascale
       end do
    end do


    call mem_dealloc(M)
    call mem_dealloc(lambda)
    call mem_dealloc(T)
    call mem_dealloc(CT)

  end subroutine orthogonalize_MOs

  !> \brief Print energy summary for CC calculation to both standard output and LSDALTON.OUT.
  subroutine print_total_energy_summary(EHF,Edft,Ecorr,EF12sing,dE_est1,dE_est2,dE_est3,doSOS)
    implicit none
    !> HF energy
    real(realk),intent(in) :: EHF,Edft
    !> Correlation energy
    real(realk),intent(in) :: Ecorr
    !> F12 singles correction
    real(realk),intent(in) :: EF12sing
    !> Estimated intrinsic DEC energy error
    real(realk),intent(in) :: dE_est1,dE_est2,dE_est3
    logical,intent(in),optional :: doSOS
    integer :: lupri

    lupri=6
    call print_total_energy_summary_lupri(EHF,Edft,Ecorr,EF12sing,dE_est1,dE_est2,&
         & dE_est3,lupri,doSOS=doSOS)

    lupri=DECinfo%output
    call print_total_energy_summary_lupri(EHF,Edft,Ecorr,EF12sing,dE_est1,dE_est2,&
         & dE_est3,lupri,doSOS=doSOS)


  end subroutine print_total_energy_summary

  !> \brief Print short energy summary (both HF and correlation) to specific logical unit number.
  !> (Necessary to place here because it is used both for DEC and for full calculation).
  !> \author Kasper Kristensen
  !> \date April 2013
  subroutine print_total_energy_summary_lupri(EHF,Edft,Ecorr,EF12sing,&
       & dE_est1,dE_est2,dE_est3,lupri,doSOS)
    implicit none
    !> HF energy
    real(realk),intent(in) :: EHF,Edft
    !> Correlation energy
    real(realk),intent(in) :: Ecorr
    !> F12 singles correction
    real(realk),intent(in) :: EF12sing
    !> Estimated intrinsic DEC energy error
    real(realk),intent(in) :: dE_est1,dE_est2,dE_est3
    !> Logical unit number to print to
    integer,intent(in) :: lupri
    !> SOS cont
    logical,intent(in),optional :: doSOS
    logical :: SOS
    real(realk) :: EF12singles

    if(DECinfo%F12singles) then
       EF12singles=EF12sing
    else
       EF12singles=0.0_realk
    end if

    SOS = .false.
    if(present(doSOS)) SOS = doSOS

    ! MODIFY FOR NEW MODEL

    ! Print summary
    write(lupri,*)
    write(lupri,*)
    write(lupri,*)
    write(lupri,*)
    write(lupri,'(13X,a)') '**********************************************************'
    if(DECinfo%full_molecular_cc) then
       write(lupri,'(13X,a,19X,a,20X,a)') '*', 'CC ENERGY SUMMARY', '*'
    else
       write(lupri,'(13X,a,19X,a,19X,a)') '*', 'DEC ENERGY SUMMARY', '*'
    end if
    write(lupri,'(13X,a)') '**********************************************************'
    write(lupri,*)
    if(DECinfo%first_order) then
       IF(.NOT.DECinfo%DFTreference)THEN
          write(lupri,'(15X,a,f20.10)') 'G: Hartree-Fock energy      :', Ehf
       ENDIF
       IF(DECinfo%DFTreference)THEN
          write(lupri,'(15X,a,f20.10)') 'G: HF energy (KS orb)       :', Ehf
          write(lupri,'(15X,a,f20.10)') 'G: DFT energy               :', Edft
       ENDIF
       write(lupri,'(15X,a,f20.10)')    'G: Correlation energy       :', Ecorr

       if(DECinfo%F12 .and. DECinfo%F12singles) then
          write(lupri,'(15X,a,f20.10)') 'G: F12 singles              :', EF12singles
       end if

       if(DECinfo%ccmodel==MODEL_MP2) then
          if (DECinfo%F12) then
             write(lupri,'(15X,a,f20.10)') 'G: Total MP2-F12 energy     :', Ehf+Ecorr+EF12singles
          else          
             write(lupri,'(15X,a,f20.10)') 'G: Total MP2 energy         :', Ehf+Ecorr      
          endif
       elseif(DECinfo%ccmodel==MODEL_RIMP2) then
          if (DECinfo%F12) then
             write(lupri,'(15X,a,f20.10)') 'G: Total RIMP2-F12 energy   :', Ehf+Ecorr+EF12singles
          else
             write(lupri,'(15X,a,f20.10)') 'G: Total RIMP2 energy       :', Ehf+Ecorr
          endif   
       elseif(DECinfo%ccmodel==MODEL_LSTHCRIMP2) then
          write(lupri,'(15X,a,f20.10)') 'G: Total LS-THC-RIMP2 energy:', Ehf+Ecorr
       elseif(DECinfo%ccmodel==MODEL_CC2) then
          write(lupri,'(15X,a,f20.10)') 'G: Total CC2 energy         :', Ehf+Ecorr
       elseif(DECinfo%ccmodel==MODEL_CCSD) then
          if (DECinfo%F12) then
             write(lupri,'(15X,a,f20.10)') 'G: Total CCSD-F12 energy     :', Ehf+Ecorr+EF12singles
          else    
             write(lupri,'(15X,a,f20.10)') 'G: Total CCSD energy         :', Ehf+Ecorr
          endif
       elseif(DECinfo%ccmodel==MODEL_CCSDpT) then
          write(lupri,'(15X,a,f20.10)') 'G: Total CCSD(T) energy     :', Ehf+Ecorr
       elseif(DECinfo%ccmodel==MODEL_SOSEX) then
         write(lupri,'(15X,a,f20.10)')  'G: HF + SOSEX energy        :', Ehf+Ecorr
         IF(DECinfo%DFTreference) then
           write(lupri,'(15X,a,f20.10)')'G: KS + SOSEX energy        :',&
             & Edft+Ecorr
         endif
       elseif(DECinfo%ccmodel==MODEL_RPA) then
         if(.not. SOS) then
           write(lupri,'(15X,a,f20.10)')  'G: HF + dRPA energy         :', Ehf+Ecorr
           IF(DECinfo%DFTreference) then
             write(lupri,'(15X,a,f20.10)')'G: KS + dRPA energy         :',&
              & Edft+Ecorr
           endif
         else
           write(lupri,'(15X,a,f20.10)')  'G: HF + SOSEX energy        :', Ehf+Ecorr
           IF(DECinfo%DFTreference) then
             write(lupri,'(15X,a,f20.10)')'G: KS + SOSEX energy        :',&
              & Edft+Ecorr
           endif
         endif
       else
          write(lupri,'(15X,A,I4,A,I4)') 'G: Unknown Energy DECinfo%ccmodel',DECinfo%ccmodel
       end if

       ! skip error print for full calculation (0 by definition)
       if(.not.DECinfo%full_molecular_cc)then  
          if(.not.(DECinfo%onlyoccpart.or.DECinfo%onlyvirtpart))then
             write(lupri,'(15X,a,f20.10)') 'G: Estimated DEC err 1 :', dE_est1
          endif
          write(lupri,'(15X,a,f20.10)') 'G: Estimated DEC err 2 :', dE_est2
          write(lupri,'(15X,a,f20.10)') 'G: Estimated DEC err 3 :', dE_est3
       end if

    else
       IF(.NOT.DECinfo%DFTreference)THEN
          write(lupri,'(15X,a,f20.10)')    'E: Hartree-Fock energy :', Ehf
       ENDIF
       IF(DECinfo%DFTreference)THEN
          write(lupri,'(15X,a,f20.10)')    'E: HF energy (KS orb)  :', Ehf
          write(lupri,'(15X,a,f20.10)')    'E: DFT energy          :', Edft
       ENDIF
       if(SOS) then
         write(lupri,'(15X,a,f20.10)')     'E: SOSEX energy        :', Ecorr
       else
         write(lupri,'(15X,a,f20.10)')     'E: Correlation energy  :', Ecorr
       endif

       if(DECinfo%F12 .and. DECinfo%F12singles) then
          write(lupri,'(15X,a,f20.10)')    'E: F12 singles         :', EF12singles
       end if

       if(DECinfo%ccmodel==MODEL_MP2) then
          if (DECinfo%F12) then
             write(lupri,'(15X,a,f20.10)') 'E: Total MP2-F12 energy:', Ehf+Ecorr+EF12singles
          else          
             write(lupri,'(15X,a,f20.10)') 'E: Total MP2 energy    :', Ehf+Ecorr      
          endif
       elseif(DECinfo%ccmodel==MODEL_RIMP2) then
          if (DECinfo%F12) then
             write(lupri,'(15X,a,f20.10)') 'E: Total RI-MP2F12 energy:', Ehf+Ecorr+EF12singles
          else          
             write(lupri,'(15X,a,f20.10)') 'E: Total RIMP2 energy     :', Ehf+Ecorr
          endif
       elseif(DECinfo%ccmodel==MODEL_LSTHCRIMP2) then
          write(lupri,'(15X,a,f20.10)') 'E: Total LS-THC-RIMP2 energy:', Ehf+Ecorr
       elseif(DECinfo%ccmodel==MODEL_MP3) then
          write(lupri,'(15X,a,f20.10)')    'E: Total MP3 energy    :', Ehf+Ecorr
       elseif(DECinfo%ccmodel==MODEL_CC2) then
          write(lupri,'(15X,a,f20.10)') 'E: Total CC2 energy         :', Ehf+Ecorr
       elseif(DECinfo%ccmodel==MODEL_CCSD) then
          if (DECinfo%F12) then
             write(lupri,'(15X,a,f20.10)') 'E: Tot CCSD-F12 energy :', Ehf+Ecorr+EF12singles
          else          
             write(lupri,'(15X,a,f20.10)') 'E: Total CCSD energy   :', Ehf+Ecorr
          endif
       elseif(DECinfo%ccmodel==MODEL_CCSDpT) then
          write(lupri,'(15X,a,f20.10)')    'E: Total CCSD(T) energy:', Ehf+Ecorr
       elseif(DECinfo%ccmodel==MODEL_SOSEX) then
         write(lupri,'(15X,a,f20.10)') 'E: HF + SOSEX energy   :', Ehf+Ecorr
         IF(DECinfo%DFTreference) then
           write(lupri,'(15X,a,f20.10)') 'E: KS + SOSEX energy   :',&
             & Edft+Ecorr
         endif
       elseif(DECinfo%ccmodel==MODEL_RPA) then
          if(.not. SOS) then
             write(lupri,'(15X,a,f20.10)') 'E: HF + dRPA energy    :', Ehf+Ecorr
           IF(DECinfo%DFTreference) then
             write(lupri,'(15X,a,f20.10)') 'E: KS + dRPA energy    :',&
              & Edft+Ecorr
           endif
          else
             write(lupri,'(15X,a,f20.10)') 'E: HF + SOSEX energy   :', Ehf+Ecorr
           IF(DECinfo%DFTreference) then
             write(lupri,'(15X,a,f20.10)') 'E: KS + SOSEX energy   :',&
              & Edft+Ecorr
           endif
          endif
       else
          write(lupri,'(15X,A,I4)') 'E: Unknown Energy DECinfo%ccmodel',DECinfo%ccmodel
       end if

       ! skip error print for full calculation (0 by definition)
       if(.not.DECinfo%full_molecular_cc)then  
          if(.not.(DECinfo%onlyoccpart.or.DECinfo%onlyvirtpart))then
             write(lupri,'(15X,a,f20.10)') 'E: Estimated DEC err 1 :', dE_est1
          endif
          write(lupri,'(15X,a,f20.10)')    'E: Estimated DEC err 2 :', dE_est2
          write(lupri,'(15X,a,f20.10)')    'E: Estimated DEC err 3 :', dE_est3
       end if

    end if
    write(lupri,*)
    write(lupri,*)


  end subroutine print_total_energy_summary_lupri

  !> \brief Print all fragment energies for given CC model.
  !> \author Kasper Kristensen
  !> \date March 2013
  subroutine print_all_fragment_energies(natoms,FragEnergies,dofrag,&
       & DistanceTable,energies)
    implicit none
    !> Number of atoms in full molecule
    integer,intent(in) :: nAtoms
    ! Fragment energies as listed in decfrag type def "energies"
    real(realk),intent(in) :: FragEnergies(natoms,natoms,ndecenergies)
    !> Which atoms are associated with a fragment?
    logical,intent(in) :: dofrag(natoms)
    !> Distances between all atoms (not changed at output, is intent(inout) for MPI purposes)
    real(realk),intent(inout) :: DistanceTable(natoms,natoms)
    !> Total DEC energies (sum of frag energies)
    real(realk),intent(in) :: energies(ndecenergies)
    !local variables 
    character(len=30) :: CorrEnergyString
    integer :: iCorrLen
    logical :: print_pair

    ! Print Header:
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) '             ==============================================='
    write(DECinfo%output,*) '             |   Print single and pair fragment energies   |'
    write(DECinfo%output,*) '             ==============================================='
    write(DECinfo%output,*)

    CorrEnergyString = 'correlation energy            '
    iCorrLen = 18
    print_pair = count(dofrag)>1 .and. (.not. DECinfo%no_pairs)
    
    select case(DECinfo%ccmodel)
    case(MODEL_MP2)
       if(.not.DECinfo%onlyvirtpart) then  
          call print_atomic_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_OCCMP2),dofrag,&
               & 'MP2 occupied single energies','AF_MP2_OCC')
       endif
       if(.not. DECinfo%onlyoccpart) then  
          call print_atomic_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_VIRTMP2),dofrag,&
               & 'MP2 virtual single energies','AF_MP2_VIR')
       endif
       if(.not.(DECinfo%onlyoccpart.or.DECinfo%onlyvirtpart)) then  
          call print_atomic_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_LAGMP2),dofrag,&
               & 'MP2 Lagrangian single energies','AF_MP2_LAG')
       end if

       if((.not.DECinfo%onlyvirtpart).and.print_pair) then  
          call print_pair_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_OCCMP2),dofrag,&
               & DistanceTable, 'MP2 occupied pair energies','PF_MP2_OCC')
       endif
       if((.not. DECinfo%onlyoccpart).and.print_pair) then  
          call print_pair_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_VIRTMP2),dofrag,&
               & DistanceTable, 'MP2 virtual pair energies','PF_MP2_VIR')          
       endif
       if((.not.(DECinfo%onlyoccpart.or.DECinfo%onlyvirtpart)).and.print_pair) then  
          call print_pair_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_LAGMP2),dofrag,&
               & DistanceTable, 'MP2 Lagrangian pair energies','PF_MP2_LAG')
       end if

       write(DECinfo%output,*)
       if(.not.DECinfo%onlyvirtpart) then  
          write(DECinfo%output,'(1X,A,A,A,g20.10)') &
               & 'MP2 occupied   ',CorrEnergyString(1:iCorrLen),' : ',energies(FRAGMODEL_OCCMP2)
       endif
       if(.not.DECinfo%onlyoccpart) then
          write(DECinfo%output,'(1X,a,a,a,g20.10)') 'MP2 virtual    ',CorrEnergyString(1:iCorrLen),' : ', &
               & energies(FRAGMODEL_VIRTMP2)
       endif
       if(.not.(DECinfo%onlyoccpart.or.DECinfo%onlyvirtpart)) then  
          write(DECinfo%output,'(1X,a,a,a,g20.10)') 'MP2 Lagrangian ',CorrEnergyString(1:iCorrLen),' : ', &
               & energies(FRAGMODEL_LAGMP2)
       end if
       write(DECinfo%output,*)
    case(MODEL_CC2)
       if(.not.DECinfo%onlyvirtpart) then  
          call print_atomic_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_OCCCC2),dofrag,&
               & 'CC2 occupied single energies','AF_CC2_OCC')
       endif
       if(.not.DECinfo%onlyoccpart) then
          call print_atomic_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_VIRTCC2),dofrag,&
               & 'CC2 virtual single energies','AF_CC2_VIR')
       end if

       if((.not.DECinfo%onlyvirtpart).and.print_pair) then  
          call print_pair_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_OCCCC2),dofrag,&
               & DistanceTable, 'CC2 occupied pair energies','PF_CC2_OCC')
       endif
       if((.not.DECinfo%onlyoccpart).and.print_pair) then
          call print_pair_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_VIRTCC2),dofrag,&
               & DistanceTable, 'CC2 virtual pair energies','PF_CC2_VIR')
       end if

       write(DECinfo%output,*)
       write(DECinfo%output,'(1X,a,a,a,g20.10)') 'CC2 occupied   ',CorrEnergyString(1:iCorrLen),' : ', &
            & energies(FRAGMODEL_OCCCC2)
       if(.not.DECinfo%onlyoccpart) then
          write(DECinfo%output,'(1X,a,a,a,g20.10)') 'CC2 virtual    ',CorrEnergyString(1:iCorrLen),' : ', &
               & energies(FRAGMODEL_VIRTCC2)
       end if
       write(DECinfo%output,*)

    case(MODEL_SOSEX)

       if(.not.DECinfo%onlyvirtpart) then  
          call print_atomic_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_OCCSOS),dofrag,&
               & 'SOSEX occupied single energies','AF_SOS_OCC')
       endif
       if(.not.DECinfo%onlyoccpart) then
          call print_atomic_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_VIRTSOS),dofrag,&
               & 'SOSEX virtual single energies','AF_SOS_VIR')
       endif

       if((.not.DECinfo%onlyvirtpart).and.print_pair) then  
          call print_pair_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_OCCSOS),dofrag,&
               & DistanceTable, 'SOSEX occupied pair energies','PF_SOS_OCC')
       endif
       if((.not.DECinfo%onlyoccpart).and.print_pair) then
          call print_pair_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_VIRTSOS),dofrag,&
               & DistanceTable, 'SOSEX virtual pair energies','PF_SOS_VIR')
       endif

       write(DECinfo%output,*)
       write(DECinfo%output,'(1X,a,a,a,g20.10)') 'SOSEX occupied   ',CorrEnergyString(1:iCorrLen),' : ', &
            & energies(FRAGMODEL_OCCSOS)
       if(.not.DECinfo%onlyoccpart) then
          write(DECinfo%output,'(1X,a,a,a,g20.10)') 'SOSEX virtual    ',CorrEnergyString(1:iCorrLen),' : ', &
               & energies(FRAGMODEL_VIRTSOS)
       end if
       write(DECinfo%output,*)

    case(MODEL_RPA)

       if(.not.DECinfo%onlyvirtpart) then  
          call print_atomic_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_OCCRPA),dofrag,&
               & 'dRPA occupied single energies','AF_RPA_OCC')
       endif
       if(.not.DECinfo%onlyoccpart) then
          call print_atomic_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_VIRTRPA),dofrag,&
               & 'dRPA virtual single energies','AF_RPA_VIR')
       endif

       if((.not.DECinfo%onlyvirtpart).and.print_pair) then  
          call print_pair_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_OCCRPA),dofrag,&
               & DistanceTable, 'dRPA occupied pair energies','PF_RPA_OCC')
       endif
       if((.not.DECinfo%onlyoccpart).and.print_pair) then
          call print_pair_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_VIRTRPA),dofrag,&
               & DistanceTable, 'dRPA virtual pair energies','PF_RPA_VIR')          
       endif

       write(DECinfo%output,*)
       write(DECinfo%output,'(1X,a,a,a,g20.10)') 'dRPA occupied   ',CorrEnergyString(1:iCorrLen),' : ', &
            & energies(FRAGMODEL_OCCRPA)
       if(.not.DECinfo%onlyoccpart) then
          write(DECinfo%output,'(1X,a,a,a,g20.10)') 'dRPA virtual    ',CorrEnergyString(1:iCorrLen),' : ', &
               & energies(FRAGMODEL_VIRTRPA)
       end if
       write(DECinfo%output,*)
       write(DECinfo%output,*)

       if(.not.DECinfo%onlyvirtpart) then  
          call print_atomic_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_OCCSOS),dofrag,&
               & 'SOSEX occupied single energies','AF_SOS_OCC')
       endif
       if(.not.DECinfo%onlyoccpart) then
          call print_atomic_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_VIRTSOS),dofrag,&
               & 'SOSEX virtual single energies','AF_SOS_VIR')
       endif

       if((.not.DECinfo%onlyvirtpart).and.print_pair) then  
          call print_pair_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_OCCSOS),dofrag,&
               & DistanceTable, 'SOSEX occupied pair energies','PF_SOS_OCC')
       endif
       if((.not.DECinfo%onlyoccpart).and.print_pair) then
          call print_pair_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_VIRTSOS),dofrag,&
               & DistanceTable, 'SOSEX virtual pair energies','PF_SOS_VIR') 
       endif

       write(DECinfo%output,*)
       write(DECinfo%output,'(1X,a,a,a,g20.10)') 'SOSEX occupied   ',CorrEnergyString(1:iCorrLen),' : ', &
            & energies(FRAGMODEL_OCCSOS)
       if(.not.DECinfo%onlyoccpart) then
          write(DECinfo%output,'(1X,a,a,a,g20.10)') 'SOSEX virtual    ',CorrEnergyString(1:iCorrLen),' : ', &
               & energies(FRAGMODEL_VIRTSOS)
       end if
       write(DECinfo%output,*)

    case(MODEL_CCSD)
       if(.not.DECinfo%CCDhack)then
          if(.not.DECinfo%onlyvirtpart) then  
             call print_atomic_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_OCCCCSD),dofrag,&
                  & 'CCSD occupied single energies','AF_CCSD_OCC')
          endif
          if(.not.DECinfo%onlyoccpart) then
             call print_atomic_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_VIRTCCSD),dofrag,&
                  & 'CCSD virtual single energies','AF_CCSD_VIR')
          end if

          if((.not.DECinfo%onlyvirtpart).and.print_pair) then  
             call print_pair_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_OCCCCSD),dofrag,&
                  & DistanceTable, 'CCSD occupied pair energies','PF_CCSD_OCC')
          endif
          if((.not.DECinfo%onlyoccpart).and.print_pair) then
             call print_pair_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_VIRTCCSD),dofrag,&
                  & DistanceTable, 'CCSD virtual pair energies','PF_CCSD_VIR')
          end if

          write(DECinfo%output,*)
          if(.not.DECinfo%onlyvirtpart) then  
             write(DECinfo%output,'(1X,a,a,a,g20.10)') 'CCSD occupied   ',CorrEnergyString(1:iCorrLen),' : ', &
                  & energies(FRAGMODEL_OCCCCSD)
          endif
          if(.not.DECinfo%onlyoccpart) then
             write(DECinfo%output,'(1X,a,a,a,g20.10)') 'CCSD virtual    ',CorrEnergyString(1:iCorrLen),' : ', &
                  & energies(FRAGMODEL_VIRTCCSD)
          end if
          write(DECinfo%output,*)
       else
          if(.not.DECinfo%onlyvirtpart) then  
             call print_atomic_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_OCCCCSD),dofrag,&
                  & 'CCD occupied single energies','AF_CCD_OCC')
          endif
          if(.not.DECinfo%onlyoccpart) then
             call print_atomic_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_VIRTCCSD),dofrag,&
                  & 'CCD virtual single energies','AF_CCD_VIR')
          end if

          if((.not.DECinfo%onlyvirtpart).and.print_pair) then  
             call print_pair_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_OCCCCSD),dofrag,&
                  & DistanceTable, 'CCD occupied pair energies','PF_CCD_OCC')
          endif
          if((.not.DECinfo%onlyoccpart).and.print_pair) then
             call print_pair_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_VIRTCCSD),dofrag,&
                  & DistanceTable, 'CCD virtual pair energies','PF_CCD_VIR')
          end if

          write(DECinfo%output,*)
          write(DECinfo%output,'(1X,a,a,a,g20.10)') 'CCD occupied   ',CorrEnergyString(1:iCorrLen),' : ', &
               & energies(FRAGMODEL_OCCCCSD)
          if(.not.DECinfo%onlyoccpart) then
             write(DECinfo%output,'(1X,a,a,a,g20.10)') 'CCD virtual    ',CorrEnergyString(1:iCorrLen),' : ', &
                  & energies(FRAGMODEL_VIRTCCSD)
          end if
          write(DECinfo%output,*)
       endif

    case(MODEL_CCSDpT)

       if(.not.DECinfo%onlyvirtpart) then  
          call print_atomic_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_OCCCCSD),dofrag,&
               & 'CCSD occupied single energies','AF_CCSD_OCC')
       endif
       if(.not.DECinfo%onlyoccpart) then
          call print_atomic_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_VIRTCCSD),dofrag,&
               & 'CCSD virtual single energies','AF_CCSD_VIR')
       end if

       if((.not.DECinfo%onlyvirtpart).and.print_pair) then  
          call print_pair_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_OCCCCSD),dofrag,&
               & DistanceTable, 'CCSD occupied pair energies','PF_CCSD_OCC')
       endif
       if((.not.DECinfo%onlyoccpart).and.print_pair) then
          call print_pair_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_VIRTCCSD),dofrag,&
               & DistanceTable, 'CCSD virtual pair energies','PF_CCSD_VIR')
       end if

       if(.not.DECinfo%onlyvirtpart) then  
          call print_atomic_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_OCCpT),dofrag,&
               & '(T) occupied single energies','AF_ParT_OCC_BOTH')
       endif
       if(.not.DECinfo%onlyoccpart) then
          call print_atomic_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_VIRTpT),dofrag,&
               & '(T) virtual single energies','AF_ParT_VIR_BOTH')
       end if

       if((.not.DECinfo%onlyvirtpart)) then  
          call print_atomic_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_OCCpT4),dofrag,&
               & '(T) occupied single energies (fourth order)','AF_ParT_OCC4')
       endif
       if((.not.DECinfo%onlyoccpart)) then
          call print_atomic_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_VIRTpT4),dofrag,&
               & '(T) virtual single energies (fourth order)','AF_ParT_VIR4')
       end if

       if((.not.DECinfo%onlyvirtpart)) then  
          call print_atomic_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_OCCpT5),dofrag,&
               & '(T) occupied single energies (fifth order)','AF_ParT_OCC5')
       endif
       if((.not.DECinfo%onlyoccpart)) then
          call print_atomic_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_VIRTpT5),dofrag,&
               & '(T) virtual single energies (fifth order)','AF_ParT_VIR5')
       end if

       if((.not.DECinfo%onlyvirtpart).and.print_pair) then  
          call print_pair_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_OCCpT),dofrag,&
               & DistanceTable, '(T) occupied pair energies','PF_ParT_OCC_BOTH')
       endif
       if((.not.DECinfo%onlyoccpart).and.print_pair) then
          call print_pair_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_VIRTpT),dofrag,&
               & DistanceTable, '(T) virtual pair energies','PF_ParT_VIR_BOTH')
       end if

       if((.not.DECinfo%onlyvirtpart).and.print_pair) then  
          call print_pair_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_OCCpT4),dofrag,&
               & DistanceTable, '(T) occupied pair energies (fourth order)','PF_ParT_OCC4')
       endif
       if((.not.DECinfo%onlyoccpart).and.print_pair) then
          call print_pair_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_VIRTpT4),dofrag,&
               & DistanceTable, '(T) virtual pair energies (fourth order)','PF_ParT_VIR4')
       end if

       if((.not.DECinfo%onlyvirtpart).and.print_pair) then  
          call print_pair_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_OCCpT5),dofrag,&
               & DistanceTable, '(T) occupied pair energies (fifth order)','PF_ParT_OCC5')
       endif
       if((.not.DECinfo%onlyoccpart).and.print_pair) then
          call print_pair_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_VIRTpT5),dofrag,&
               & DistanceTable, '(T) virtual pair energies (fifth order)','PF_ParT_VIR5')
       end if

       write(DECinfo%output,*)
       write(DECinfo%output,*)
       write(DECinfo%output,*)
       if(.not.DECinfo%onlyvirtpart) then  
          write(DECinfo%output,'(1X,a,a,a,g20.10)') 'CCSD occupied ',CorrEnergyString(1:iCorrLen),' : ', &
               & energies(FRAGMODEL_OCCCCSD)
          write(DECinfo%output,'(1X,a,g20.10)') '(T) occupied correlation energy  : ', energies(FRAGMODEL_OCCpT)
          write(DECinfo%output,'(1X,a,g20.10)') '(T) occupied 4th order energy    : ', energies(FRAGMODEL_OCCpT4)
          write(DECinfo%output,'(1X,a,g20.10)') '(T) occupied 5th order energy    : ', energies(FRAGMODEL_OCCpT5)
       endif
       if(.not.DECinfo%onlyoccpart) then
          write(DECinfo%output,'(1X,a,a,a,g20.10)') 'CCSD virtual  ',CorrEnergyString(1:iCorrLen),' : ', &
               & energies(FRAGMODEL_VIRTCCSD)
          write(DECinfo%output,'(1X,a,g20.10)') '(T) virtual  correlation energy  : ', energies(FRAGMODEL_VIRTpT)
          write(DECinfo%output,'(1X,a,g20.10)') '(T) virtual  4th order energy    : ', energies(FRAGMODEL_VIRTpT4)
          write(DECinfo%output,'(1X,a,g20.10)') '(T) virtual  5th order energy    : ', energies(FRAGMODEL_VIRTpT5)
       end if
       write(DECinfo%output,*)
       if(.not.DECinfo%onlyvirtpart) then  
          write(DECinfo%output,'(1X,a,a,a,g20.10)') 'Total CCSD(T) occupied ',CorrEnergyString(1:iCorrLen),' : ', &
               & energies(FRAGMODEL_OCCCCSD)+energies(FRAGMODEL_OCCpT)
       endif
       if(.not.DECinfo%onlyoccpart) then
          write(DECinfo%output,'(1X,a,a,a,g20.10)') 'Total CCSD(T) virtual  ',CorrEnergyString(1:iCorrLen),' : ', &
               & energies(FRAGMODEL_VIRTCCSD)+energies(FRAGMODEL_VIRTpT)
       end if
       write(DECinfo%output,*)

    case(MODEL_RIMP2)
       if(.not.DECinfo%onlyvirtpart) then  
          call print_atomic_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_OCCRIMP2),dofrag,&
               & 'RI-MP2 occupied single energies','AF_RI_MP2_OCC')
       endif
       if(.not. DECinfo%onlyoccpart) then  
          call print_atomic_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_VIRTRIMP2),dofrag,&
               & 'RI-MP2 virtual single energies','AF_RI_MP2_VIR')
       endif
       if(.not.(DECinfo%onlyoccpart.or.DECinfo%onlyvirtpart)) then  
          call print_atomic_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_LAGRIMP2),dofrag,&
               & 'RI-MP2 Lagrangian single energies','AF_RI_MP2_LAG')
       end if

       if((.not.DECinfo%onlyvirtpart).and.print_pair) then  
          call print_pair_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_OCCRIMP2),dofrag,&
               & DistanceTable, 'RI-MP2 occupied pair energies','PF_RI_MP2_OCC')
       endif
       if((.not. DECinfo%onlyoccpart).and.print_pair) then  
          call print_pair_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_VIRTRIMP2),dofrag,&
               & DistanceTable, 'RI-MP2 virtual pair energies','PF_RI_MP2_VIR')          
       endif
       if((.not.(DECinfo%onlyoccpart.or.DECinfo%onlyvirtpart)).and.print_pair) then  
          call print_pair_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_LAGRIMP2),dofrag,&
               & DistanceTable, 'RI-MP2 Lagrangian pair energies','PF_RI_MP2_LAG')
       end if

       write(DECinfo%output,*)
       if(.not.DECinfo%onlyvirtpart) then  
          write(DECinfo%output,'(1X,A,A,A,g20.10)') &
               & 'RI-MP2 occupied   ',CorrEnergyString(1:iCorrLen),' : ',energies(FRAGMODEL_OCCRIMP2)
       endif
       if(.not.DECinfo%onlyoccpart) then
          write(DECinfo%output,'(1X,a,a,a,g20.10)') &
               & 'RI-MP2 virtual    ',CorrEnergyString(1:iCorrLen),' : ', energies(FRAGMODEL_VIRTRIMP2)
       endif
       if(.not.(DECinfo%onlyoccpart.or.DECinfo%onlyvirtpart)) then  
          write(DECinfo%output,'(1X,a,a,a,g20.10)') &
               &'RI-MP2 Lagrangian ',CorrEnergyString(1:iCorrLen),' : ', energies(FRAGMODEL_LAGRIMP2)
       end if
       write(DECinfo%output,*)

    case(MODEL_LSTHCRIMP2)
       if(.not.DECinfo%onlyvirtpart) then  
          call print_atomic_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_OCCLSTHCRIMP2),dofrag,&
               & 'LS-THC-RI-MP2 occupied single energies','AF_RI_MP2_OCC')
       endif
       if(.not. DECinfo%onlyoccpart) then  
          call print_atomic_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_VIRTLSTHCRIMP2),dofrag,&
               & 'LS-THC-RI-MP2 virtual single energies','AF_RI_MP2_VIR')
       endif
       if(.not.(DECinfo%onlyoccpart.or.DECinfo%onlyvirtpart)) then  
          call print_atomic_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_LAGLSTHCRIMP2),dofrag,&
               & 'LS-THC-RI-MP2 Lagrangian single energies','AF_RI_MP2_LAG')
       end if

       if((.not.DECinfo%onlyvirtpart).and.print_pair) then  
          call print_pair_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_OCCLSTHCRIMP2),dofrag,&
               & DistanceTable, 'LS-THC-RI-MP2 occupied pair energies','PF_RI_MP2_OCC')
       endif
       if((.not. DECinfo%onlyoccpart).and.print_pair) then  
          call print_pair_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_VIRTLSTHCRIMP2),dofrag,&
               & DistanceTable, 'LS-THC-RI-MP2 virtual pair energies','PF_RI_MP2_VIR')          
       endif
       if((.not.(DECinfo%onlyoccpart.or.DECinfo%onlyvirtpart)).and.print_pair) then  
          call print_pair_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_LAGLSTHCRIMP2),dofrag,&
               & DistanceTable, 'LS-THC-RI-MP2 Lagrangian pair energies','PF_RI_MP2_LAG')
       end if

       write(DECinfo%output,*)
       if(.not.DECinfo%onlyvirtpart) then  
          write(DECinfo%output,'(1X,A,A,A,g20.10)') &
               & 'LS-THC-RI-MP2 occupied   ',CorrEnergyString(1:iCorrLen),' : ',energies(FRAGMODEL_OCCLSTHCRIMP2)
       endif
       if(.not.DECinfo%onlyoccpart) then
          write(DECinfo%output,'(1X,a,a,a,g20.10)') &
               & 'LS-THC-RI-MP2 virtual    ',CorrEnergyString(1:iCorrLen),' : ', energies(FRAGMODEL_VIRTLSTHCRIMP2)
       endif
       if(.not.(DECinfo%onlyoccpart.or.DECinfo%onlyvirtpart)) then  
          write(DECinfo%output,'(1X,a,a,a,g20.10)') &
               &'LS-THC-RI-MP2 Lagrangian ',CorrEnergyString(1:iCorrLen),' : ', energies(FRAGMODEL_LAGLSTHCRIMP2)
       end if
       write(DECinfo%output,*)
    case default
       ! MODIFY FOR NEW MODEL
       ! If you implement new model, please print the fragment energies here,
       ! see FRAGMODEL_* definitions in dec_typedef.F90.
       write(DECinfo%output,*) 'WARNING: print_all_fragment_energies needs implementation &
            & for model: ', DECinfo%ccmodel
    end select

#ifdef MOD_UNRELEASED
    ! MODIFY FOR NEW CORRECTION
    if(DECInfo%F12) then

       select case(DECinfo%ccmodel)

       case(MODEL_RIMP2)
          call print_atomic_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_RIMP2f12),dofrag,&
             & 'RIMP2F12 occupied single energies','AF_MP2f12_OCC')

          if (print_pair) call print_pair_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_RIMP2f12),dofrag,&
             & DistanceTable, 'RIMP2f12 occupied pair energies','PF_MP2F12f12_OCC')
          
          write(DECinfo%output,*)   
          write(DECinfo%output,'(1X,a,f20.10)') 'RIMP2 CORRECTION TO ENERGY:    ', energies(FRAGMODEL_OCCRIMP2)  
          write(DECinfo%output,'(1X,a,f20.10)') 'F12 CORRECTION TO MP2 ENERGY:  ', energies(FRAGMODEL_RIMP2f12)
          write(DECinfo%output,'(1X,a,f20.10)') 'RIMP2-F12 CORRELATION ENERGY:  ', &
             & energies(FRAGMODEL_OCCRIMP2) + energies(FRAGMODEL_RIMP2f12)
          write(DECinfo%output,*)       


       case(MODEL_MP2)
          call print_atomic_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_MP2f12),dofrag,&
             & 'MP2F12 occupied single energies','AF_MP2f12_OCC')
          if (print_pair) call print_pair_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_MP2f12),dofrag,&
             & DistanceTable, 'MP2f12 occupied pair energies','PF_MP2F12f12_OCC')
         
          write(DECinfo%output,*)   
          write(DECinfo%output,'(1X,a,f20.10)') 'MP2 CORRECTION TO ENERGY:      ', energies(FRAGMODEL_OCCMP2)  
          write(DECinfo%output,'(1X,a,f20.10)') 'F12 CORRECTION TO MP2 ENERGY:  ', energies(FRAGMODEL_MP2f12)
          write(DECinfo%output,'(1X,a,f20.10)') 'MP2-F12 CORRELATION ENERGY:    ', &
             & energies(FRAGMODEL_OCCMP2) + energies(FRAGMODEL_MP2f12)
          write(DECinfo%output,*)       

       case(MODEL_CCSD)
          call print_atomic_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_CCSDf12),dofrag,&
             & 'CCSDF12 occupied single energies','AF_CCSDf12_OCC')
          if (print_pair) call print_pair_fragment_energies(natoms,FragEnergies(:,:,FRAGMODEL_CCSDf12),dofrag,&
             & DistanceTable, 'CCSDf12 occupied pair energies','PF_CCSDf12_OCC')

          write(DECinfo%output,*)   
          write(DECinfo%output,'(1X,a,f20.10)') 'CCSD CORRECTION TO ENERGY:     ', energies(FRAGMODEL_OCCCCSD)
          write(DECinfo%output,'(1X,a,f20.10)') 'F12 CORRECTION TO CCSD ENERGY: ', energies(FRAGMODEL_CCSDf12)
          write(DECinfo%output,'(1X,a,f20.10)') 'CCSD-F12 CORRELATION ENERGY:   ', &
             & energies(FRAGMODEL_OCCCCSD) + energies(FRAGMODEL_CCSDf12)

       end select

       endif
#endif

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) '============================================================================='

  end subroutine print_all_fragment_energies


  !> \brief: print out solver fragment and pair interaction energies for full molecule 
  !          calculation. This routine should print the information in the same 
  !          way kasper's routine, print_all_fragment_energies does for DEC.
  !
  !> \author: Janus Juul Eriksen and Pablo Baudin
  !> \date: February 2013
  subroutine print_fragment_energies_full(nfrags,FragEnergies,ccenergies,dofrag,distancetable)

    implicit none

    !> number of atoms in molecule
    integer, intent(in) :: nfrags
    !> matrices containing Frag. energies and interatomic distances
    real(realk), intent(in) :: FragEnergies(nfrags,nfrags,8), distancetable(nfrags,nfrags)
    !> Total cc energies:
    real(realk), intent(in) :: ccenergies(8)
    !> vector handling how the orbitals are assigned?
    logical, intent(inout) :: dofrag(nfrags)

    !> local variables 
    character(len=30) :: CorrEnergyString
    integer :: iCorrLen, cc_sol_o, pT_4_o, pT_5_o, pT_full_o
    integer :: cc_sol_v, pT_4_v, pT_5_v, pT_full_v
    logical :: print_pair

    print_pair = count(dofrag)>1 .and. (.not. DECinfo%no_pairs)
    CorrEnergyString = 'correlation energy            '
    iCorrLen = 18
    cc_sol_o  = 1
    cc_sol_v  = 2
    pT_full_o = 3
    pT_full_v = 4
    pT_4_o    = 5
    pT_4_v    = 6
    pT_5_o    = 7
    pT_5_v    = 8

    ! Print Header:
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) '             ==============================================='
    write(DECinfo%output,*) '             |   Print single and pair fragment energies   |'
    write(DECinfo%output,*) '             ==============================================='
    write(DECinfo%output,*)

    if(.not.DECinfo%CCDhack)then
       if( DECinfo%ccmodel == MODEL_RPA)then
          if (.not.DECinfo%OnlyVirtPart) then
             call print_atomic_fragment_energies(nfrags,FragEnergies(:,:,cc_sol_o),dofrag,&
                & 'RPA occupied single energies','AF_RPA_OCC')
          end if
          if (.not.DECinfo%OnlyOccPart) then
             call print_atomic_fragment_energies(nfrags,FragEnergies(:,:,cc_sol_v),dofrag,&
                & 'RPA virtual single energies','AF_RPA_VIR')
          end if
          if (.not.DECinfo%OnlyVirtPart .and. print_pair) then
             call print_pair_fragment_energies(nfrags,FragEnergies(:,:,cc_sol_o),&
                & dofrag,Distancetable, 'RPA occupied pair energies','PF_RPA_OCC')
          end if
          if (.not.DECinfo%OnlyOccPart .and. print_pair) then
             call print_pair_fragment_energies(nfrags,FragEnergies(:,:,cc_sol_v),&
                & dofrag,Distancetable, 'RPA virtual pair energies','PF_RPA_VIR')
          end if

          write(DECinfo%output,*)
          write(DECinfo%output,*)
          write(DECinfo%output,*)
          write(DECinfo%output,'(1X,A,A,A,g20.10)') 'dRPA ', &
             & CorrEnergyString(1:iCorrLen),' : ',ccenergies(cc_sol_o)
          write(DECinfo%output,*)

       else if( DECinfo%ccmodel == MODEL_SOSEX)then
          if (.not.DECinfo%OnlyVirtPart) then
             call print_atomic_fragment_energies(nfrags,FragEnergies(:,:,cc_sol_o),dofrag,&
                & 'SOSEX occupied single energies','AF_SOS_OCC')
          end if
          if (.not.DECinfo%OnlyOccPart) then
             call print_atomic_fragment_energies(nfrags,FragEnergies(:,:,cc_sol_v),dofrag,&
                & 'SOSEX virtual single energies','AF_SOS_VIR')
          end if
          if (.not.DECinfo%OnlyVirtPart .and. print_pair) then
             call print_pair_fragment_energies(nfrags,FragEnergies(:,:,cc_sol_o),&
                & dofrag,Distancetable, 'SOSEX occupied pair energies','PF_SOS_OCC')
          end if
          if (.not.DECinfo%OnlyOccPart .and. print_pair) then
             call print_pair_fragment_energies(nfrags,FragEnergies(:,:,cc_sol_v),&
                & dofrag,Distancetable, 'SOSEX virtual pair energies','PF_SOS_VIR')
          end if

          write(DECinfo%output,*)
          write(DECinfo%output,*)
          write(DECinfo%output,*)
          write(DECinfo%output,'(1X,A,A,A,g20.10)') 'SOSEX ', &
             & CorrEnergyString(1:iCorrLen),' : ',ccenergies(cc_sol_o)
          write(DECinfo%output,*)

       else if( DECinfo%ccmodel == MODEL_MP2)then
          if (.not.DECinfo%OnlyVirtPart) then
             call print_atomic_fragment_energies(nfrags,FragEnergies(:,:,cc_sol_o),dofrag,&
                & 'MP2 occupied single energies','AF_MP2_OCC')
          end if
          if (.not.DECinfo%OnlyOccPart) then
             call print_atomic_fragment_energies(nfrags,FragEnergies(:,:,cc_sol_v),dofrag,&
                & 'MP2 virtual single energies','AF_MP2_VIR')
          end if
          if (.not.DECinfo%OnlyVirtPart .and. print_pair) then
             call print_pair_fragment_energies(nfrags,FragEnergies(:,:,cc_sol_o),&
                & dofrag,Distancetable, 'MP2 occupied pair energies','PF_MP2_OCC')
          end if
          if (.not.DECinfo%OnlyOccPart .and. print_pair) then
             call print_pair_fragment_energies(nfrags,FragEnergies(:,:,cc_sol_v),&
                & dofrag,Distancetable, 'MP2 virtual pair energies','PF_MP2_VIR')
          end if

          write(DECinfo%output,*)
          write(DECinfo%output,*)
          write(DECinfo%output,*)
          write(DECinfo%output,'(1X,A,A,A,g20.10)') 'MP2 ', &
             & CorrEnergyString(1:iCorrLen),' : ',ccenergies(cc_sol_o)
          write(DECinfo%output,*)

       else if( DECinfo%ccmodel == MODEL_CC2 )then
          if (.not.DECinfo%OnlyVirtPart) then
             call print_atomic_fragment_energies(nfrags,FragEnergies(:,:,cc_sol_o),dofrag,&
                & 'CC2 occupied single energies','AF_CC2_OCC')
          end if
          if (.not.DECinfo%OnlyOccPart) then
             call print_atomic_fragment_energies(nfrags,FragEnergies(:,:,cc_sol_v),dofrag,&
                & 'CC2 virtual single energies','AF_CC2_VIR')
          end if
          if (.not.DECinfo%OnlyVirtPart .and. print_pair) then
             call print_pair_fragment_energies(nfrags,FragEnergies(:,:,cc_sol_o),&
                & dofrag,Distancetable, 'CC2 occupied pair energies','PF_CC2_OCC')
          end if
          if (.not.DECinfo%OnlyOccPart .and. print_pair) then
             call print_pair_fragment_energies(nfrags,FragEnergies(:,:,cc_sol_v),&
                & dofrag,Distancetable, 'CC2 virtual pair energies','PF_CC2_VIR')
          end if

          write(DECinfo%output,*)
          write(DECinfo%output,*)
          write(DECinfo%output,*)
          write(DECinfo%output,'(1X,A,A,A,g20.10)') 'CC2 ', &
             & CorrEnergyString(1:iCorrLen),' : ',ccenergies(cc_sol_o)
          write(DECinfo%output,*)

       else if( DECinfo%ccmodel == MODEL_CCSD )then 
          if (.not.DECinfo%OnlyVirtPart) then
             call print_atomic_fragment_energies(nfrags,FragEnergies(:,:,cc_sol_o),dofrag,&
                & 'CCSD occupied single energies','AF_CCSD_OCC')
          end if
          if (.not.DECinfo%OnlyOccPart) then
             call print_atomic_fragment_energies(nfrags,FragEnergies(:,:,cc_sol_v),dofrag,&
                & 'CCSD virtual single energies','AF_CCSD_VIR')
          end if
          if (.not.DECinfo%OnlyVirtPart .and. print_pair) then
             call print_pair_fragment_energies(nfrags,FragEnergies(:,:,cc_sol_o),&
                & dofrag,Distancetable, 'CCSD occupied pair energies','PF_CCSD_OCC')
          end if
          if (.not.DECinfo%OnlyOccPart .and. print_pair) then
             call print_pair_fragment_energies(nfrags,FragEnergies(:,:,cc_sol_v),&
                & dofrag,Distancetable, 'CCSD virtual pair energies','PF_CCSD_VIR')
          end if

          write(DECinfo%output,*)
          write(DECinfo%output,*)
          write(DECinfo%output,*)
          write(DECinfo%output,'(1X,A,A,A,g20.10)') 'CCSD ', &
             & CorrEnergyString(1:iCorrLen),' : ',ccenergies(cc_sol_o)
          write(DECinfo%output,*)

       else if( DECinfo%ccmodel == MODEL_CCSDpT )then
          ! CCSD part single fragment
          if (.not.DECinfo%OnlyVirtPart) then
             call print_atomic_fragment_energies(nfrags,FragEnergies(:,:,cc_sol_o),dofrag,&
                & 'CCSD occupied single energies','AF_CCSD_OCC')
          end if
          if (.not.DECinfo%OnlyOccPart) then
             call print_atomic_fragment_energies(nfrags,FragEnergies(:,:,cc_sol_v),dofrag,&
                & 'CCSD virtual single energies','AF_CCSD_VIR')
          end if
          ! CCSD part pair fragment
          if (.not.DECinfo%OnlyVirtPart .and. print_pair) then
             call print_pair_fragment_energies(nfrags,FragEnergies(:,:,cc_sol_o),&
                & dofrag,Distancetable, 'CCSD occupied pair energies','PF_CCSD_OCC')
          end if
          if (.not.DECinfo%OnlyOccPart .and. print_pair) then
             call print_pair_fragment_energies(nfrags,FragEnergies(:,:,cc_sol_v),&
                & dofrag,Distancetable, 'CCSD virtual pair energies','PF_CCSD_VIR')
          end if
          ! (T) full single fragment
          if (.not.DECinfo%OnlyVirtPart) then
             call print_atomic_fragment_energies(nfrags,FragEnergies(:,:,pT_full_o),dofrag,&
                & '(T) occupied single energies','AF_ParT_OCC_BOTH')
          end if
          if (.not.DECinfo%OnlyOccPart) then
             call print_atomic_fragment_energies(nfrags,FragEnergies(:,:,pT_full_v),dofrag,&
                & '(T) virtual single energies','AF_ParT_VIR_BOTH')
          end if
          ! [4] single fragment
          if (.not.DECinfo%OnlyVirtPart) then
             call print_atomic_fragment_energies(nfrags,FragEnergies(:,:,pT_4_o),dofrag,&
                & '(T) occupied single energies (fourth order)','AF_ParT_OCC4')
          end if
          if (.not.DECinfo%OnlyOccPart) then
             call print_atomic_fragment_energies(nfrags,FragEnergies(:,:,pT_4_v),dofrag,&
                & '(T) virtual single energies (fourth order)','AF_ParT_VIR4')
          end if
          ! [5] single fragment
          if (.not.DECinfo%OnlyVirtPart) then
             call print_atomic_fragment_energies(nfrags,FragEnergies(:,:,pT_5_o),dofrag,&
                & '(T) occupied single energies (fifth order)','AF_ParT_OCC5')
          end if
          if (.not.DECinfo%OnlyOccPart) then
             call print_atomic_fragment_energies(nfrags,FragEnergies(:,:,pT_5_v),dofrag,&
                & '(T) virtual single energies (fifth order)','AF_ParT_VIR5')
          end if

          ! (T) full pair fragment
          if (.not.DECinfo%OnlyVirtPart .and. print_pair) then
             call print_pair_fragment_energies(nfrags,FragEnergies(:,:,pT_full_o),&
                & dofrag,Distancetable, '(T) occupied pair energies','PF_ParT_OCC_BOTH')
          end if
          if (.not.DECinfo%OnlyOccPart .and. print_pair) then
             call print_pair_fragment_energies(nfrags,FragEnergies(:,:,pT_full_v),&
                & dofrag,Distancetable, '(T) virtual pair energies','PF_ParT_VIR_BOTH')
          end if
          ! [4] pair fragment
          if (.not.DECinfo%OnlyVirtPart .and. print_pair) then
             call print_pair_fragment_energies(nfrags,FragEnergies(:,:,pT_4_o),&
                & dofrag,Distancetable, '(T) occupied pair energies (fourth order)','PF_ParT_OCC4')
          end if
          if (.not.DECinfo%OnlyOccPart .and. print_pair) then
             call print_pair_fragment_energies(nfrags,FragEnergies(:,:,pT_4_v),&
                & dofrag,Distancetable, '(T) virtual pair energies (fourth order)','PF_ParT_VIR4')
          end if
          ! [5] pair fragment
          if (.not.DECinfo%OnlyVirtPart .and. print_pair) then
             call print_pair_fragment_energies(nfrags,FragEnergies(:,:,pT_5_o),&
                & dofrag,Distancetable, '(T) occupied pair energies (fifth order)','PF_ParT_OCC5')
          end if
          if (.not.DECinfo%OnlyOccPart .and. print_pair) then
             call print_pair_fragment_energies(nfrags,FragEnergies(:,:,pT_5_v),&
                & dofrag,Distancetable, '(T) virtual pair energies (fifth order)','PF_ParT_VIR5')
          end if

          write(DECinfo%output,*)
          write(DECinfo%output,*)
          write(DECinfo%output,*)
          write(DECinfo%output,'(1X,a,a,a,g20.10)') 'CCSD ', &
             & CorrEnergyString(1:iCorrLen),' : ', ccenergies(cc_sol_o)
          write(DECinfo%output,'(1X,a,g20.10)') '(T) correlation energy  : ', &
             & ccenergies(pT_full_o)
          write(DECinfo%output,'(1X,a,g20.10)') '(T) 4th order energy    : ', &
             & ccenergies(pT_4_o)
          write(DECinfo%output,'(1X,a,g20.10)') '(T) 5th order energy    : ', &
             & ccenergies(pT_5_o)
          write(DECinfo%output,*)
          write(DECinfo%output,'(1X,a,a,a,g20.10)') 'Total CCSD(T) ', &
             & CorrEnergyString(1:iCorrLen),' : ', ccenergies(cc_sol_o)+ccenergies(pT_full_o)
          write(DECinfo%output,*)

       else
          call lsquit("ERROR(print_fragment_energies_full) model not implemented",-1)
       endif
    else
       if (.not.DECinfo%OnlyVirtPart) then
          call print_atomic_fragment_energies(nfrags,FragEnergies(:,:,cc_sol_o),dofrag,&
             & 'CCD occupied single energies','AF_CCD_OCC')
       end if
       if (.not.DECinfo%OnlyOccPart) then
          call print_atomic_fragment_energies(nfrags,FragEnergies(:,:,cc_sol_v),dofrag,&
             & 'CCD virtual single energies','AF_CCD_VIR')
       end if
       if (.not.DECinfo%OnlyVirtPart .and. print_pair) then
          call print_pair_fragment_energies(nfrags,FragEnergies(:,:,cc_sol_o),&
             & dofrag,Distancetable, 'CCD occupied pair energies','PF_CCD_OCC')
       end if
       if (.not.DECinfo%OnlyOccPart .and. print_pair) then
          call print_pair_fragment_energies(nfrags,FragEnergies(:,:,cc_sol_v),&
             & dofrag,Distancetable, 'CCD virtual pair energies','PF_CCD_VIR')
       end if

       write(DECinfo%output,*)
       write(DECinfo%output,*)
       write(DECinfo%output,*)
       write(DECinfo%output,'(1X,A,A,A,g20.10)') 'CCD ', &
          & CorrEnergyString(1:iCorrLen),' : ',ccenergies(cc_sol_o)
       write(DECinfo%output,*)

    endif

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*)'============================================================================='

  end subroutine print_fragment_energies_full


  !> \brief Print atomic fragment energies
  !> \author Kasper Kristensen
  !> \date September 2012
  subroutine print_atomic_fragment_energies(natoms,FragEnergies,dofrag,headline,greplabel)

    implicit none
    !> Number of atoms in full molecule
    integer,intent(in) :: nAtoms
    ! Fragment energies 
    real(realk),intent(in) :: FragEnergies(natoms,natoms)
    !> Which atoms are associated with a fragment?
    logical,intent(in) :: dofrag(natoms)
    !> Character string to print as headline
    character(*),intent(in) :: headline
    !> Label to print after each energy for easy grepping
    character(*),intent(in) :: greplabel
    integer :: i

    IF(DECinfo%RepeatAF)THEN 
       write(DECinfo%output,*)
       write(DECinfo%output,*)
       write(DECinfo%output,*) '================================================================='
       write(DECinfo%output,*) trim(headline)
       write(DECinfo%output,*) '================================================================='
       write(DECinfo%output,*)
       write(DECinfo%output,*) 'Fragment       Energy'
       do i=1,natoms
          if(.not. dofrag(i)) cycle
          write(DECinfo%output,'(I6,3X,g20.10,2a)') i,FragEnergies(i,i), "    ",greplabel
       end do
    ELSE
       write(DECinfo%output,*)
       write(DECinfo%output,*)
       write(DECinfo%output,*) '================================================================='
       write(DECinfo%output,*) trim(headline)
       write(DECinfo%output,*) '================================================================='
       write(DECinfo%output,*)
       write(DECinfo%output,*) 'Fragment       Energy'
       do i=1,natoms
          if(.not. dofrag(i)) cycle
          write(DECinfo%output,'(I6,3X,g20.10,2a)') i,FragEnergies(i,i), "    ",greplabel
       end do
    ENDIF

  end subroutine print_atomic_fragment_energies


  !> \brief Print pair fragment energies
  !> \author Kasper Kristensen
  !> \date September 2012
  subroutine print_pair_fragment_energies(natoms,FragEnergies,dofrag,&
       & DistanceTable, headline, greplabel)

    implicit none

    !> Number of atoms in full molecule
    integer,intent(in) :: nAtoms
    ! Fragment energies 
    real(realk),intent(in) :: FragEnergies(natoms,natoms)
    !> Which atoms are associated with a fragment?
    logical,intent(in) :: dofrag(natoms)
    !> Distances between all atoms (a.u.)
    real(realk),intent(in) :: DistanceTable(natoms,natoms)
    !> Character string to print as headline
    character(*),intent(in) :: headline
    !> Label to print after each energy for easy grepping
    character(*),intent(in) :: greplabel
    integer :: i,j
    real(realk) :: pairdist,thr


    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) '================================================================='
    write(DECinfo%output,*) trim(headline)
    write(DECinfo%output,*) '================================================================='
    write(DECinfo%output,*)
    write(DECinfo%output,'(2X,a)') 'Frag1  Frag2     Dist(Ang)        Energy'
    thr=1.0E-15_realk
    do i=1,natoms
       do j=i+1,natoms

          ! Skip if no fragment
          if(.not. dofrag(i)) cycle
          if(.not. dofrag(j)) cycle

          pairdist = DistanceTable(i,j)

          ! Only print if pair distance is below threshold and nonzero
          DistanceCheck: if(pairdist < DECinfo%pair_distance_threshold &
               & .and. abs(FragEnergies(i,j)) > thr  ) then
             write(DECinfo%output,'(I6,2X,I6,2X,g14.5,2X,g20.10,2a)') &
                  & i,j,pairdist*bohr_to_angstrom, FragEnergies(i,j),"    ",greplabel
          end if DistanceCheck

       end do
    end do

  end subroutine print_pair_fragment_energies


  !> \brief Print specific set  of pair fragment energies (modified version 
  !         of kasper's print_pair_fragment_energies routine)
  !> \author Pablo Baudin
  !> \date June 2014
  subroutine print_spec_pair_fragment_energies(natoms,npairs,pair_set,FragEnergies,&
       & dofrag, DistanceTable, headline, greplabel)

    implicit none

    !> Number of pairs to print:
    integer,intent(in) :: npairs
    !> Number of atoms in the molecule:
    integer,intent(in) :: natoms
    !> Atomic indices of the pairs
    integer,intent(in) :: pair_set(npairs,2)
    ! Fragment energies 
    real(realk),intent(in) :: FragEnergies(natoms,natoms)
    !> Logical vector describing which atoms have orbitals assigned 
    !> (i.e., which atoms to consider in atomic fragment calculations)
    logical,intent(in) :: dofrag(natoms)
    !> Distances between all atoms (a.u.)
    real(realk),intent(in) :: DistanceTable(natoms,natoms)
    !> Character string to print as headline
    character(*),intent(in) :: headline
    !> Label to print after each energy for easy grepping
    character(*),intent(in) :: greplabel
    integer :: i,P,Q
    real(realk) :: pairdist, thr


    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) '================================================================='
    write(DECinfo%output,*) trim(headline)
    write(DECinfo%output,*) '================================================================='
    write(DECinfo%output,*)
    write(DECinfo%output,'(2X,a)') 'Atom1  Atom2     Dist(Ang)        Energy'
    thr=1.0E-15_realk
    do i=1,npairs
      P=pair_set(i,1)
      Q=pair_set(i,2)

      ! Skip if no fragment
      if(.not. dofrag(P)) cycle
      if(.not. dofrag(Q)) cycle

      pairdist = DistanceTable(P,Q)

      ! Only print if pair distance is below threshold and nonzero
      DistanceCheck: if(pairdist < DECinfo%pair_distance_threshold &
           & .and. abs(FragEnergies(P,Q)) > thr  ) then

        write(DECinfo%output,'(I6,2X,I6,2X,g14.5,2X,g20.10,2a)') &
             & P,Q,pairdist*bohr_to_angstrom, FragEnergies(P,Q),"    ",greplabel
      end if DistanceCheck

    end do

  end subroutine print_spec_pair_fragment_energies


  !> \brief Get total number of atomic fragments + pair fragments
  !> \author Kasper Kristensen
  !> \date October 2013
  function get_total_number_of_fragments(natoms,dofrag,DistanceTable) result(njobs)
    implicit none
    !> Number of atoms in molecule
    integer,intent(in) :: natoms
    !> Logical vector describing which atoms have orbitals assigned 
    !> (i.e., which atoms to consider in atomic fragment calculations)
    logical,intent(in) :: dofrag(natoms)
    !> Distance table with interatomic distances
    real(realk),intent(in) :: DistanceTable(natoms,natoms)
    integer :: njobs
    integer :: naf,npf,i,j


    ! Number of atomic fragments
    naf = count(dofrag)

    ! Number of pair fragments (do not include pairs with interatomic distance beyond pair threshold)
    npf=0
    iloop: do i=1,natoms
       if(.not. dofrag(i)) cycle iloop
       jloop: do j=i+1,natoms
          if(.not. dofrag(j)) cycle iloop

          ! Pair distance below threshold?
          if(DistanceTable(i,j)<DECinfo%pair_distance_threshold) then
             npf=npf+1
          end if

       end do jloop
    end do iloop

    ! Number of jobs: Atomic frags + pair frags
    njobs = naf+npf

  end function get_total_number_of_fragments


  !> \brief Extract fragment energies for occupied partitioning for given CC model
  !> from set of all fragment energies.
  !> \author Kasper Kristensen
  !> \date October 2013
  subroutine get_occfragenergies(natoms,ccmodel,FragEnergiesAll,FragEnergies)
    implicit none
    !> Number of atoms in molecule
    integer,intent(in) :: natoms
    !> CC model according to MODEL_* conventions in dec_typedef.F90
    integer,intent(in) :: ccmodel
    !> Fragment energies for all models (see FRAGMODEL_* in dec_typedef.F90)
    real(realk),intent(in) :: FragEnergiesAll(natoms,natoms,ndecenergies)
    !> Fragment energies for occupied partitioning for given CC model
    real(realk),intent(inout) :: FragEnergies(natoms,natoms)

    ! MODIFY FOR NEW MODEL
    select case(ccmodel)
    case(MODEL_MP2)
       FragEnergies=FragEnergiesAll(:,:,FRAGMODEL_OCCMP2)

    case(MODEL_CC2)
       FragEnergies=FragEnergiesAll(:,:,FRAGMODEL_OCCCC2)

    case(MODEL_RPA)

       FragEnergies=FragEnergiesAll(:,:,FRAGMODEL_OCCRPA)

    case(MODEL_SOSEX)

       FragEnergies=FragEnergiesAll(:,:,FRAGMODEL_OCCSOS)

    case(MODEL_CCSD)
       FragEnergies=FragEnergiesAll(:,:,FRAGMODEL_OCCCCSD)

    case(MODEL_CCSDpT)
       ! CCSD(T): Add CCSD and (T) contributions using occupied partitioning
       FragEnergies=FragEnergiesAll(:,:,FRAGMODEL_OCCCCSD) &
            & + FragEnergiesAll(:,:,FRAGMODEL_OCCpT)

    case(MODEL_RIMP2)
       FragEnergies=FragEnergiesAll(:,:,FRAGMODEL_OCCRIMP2)
    case(MODEL_LSTHCRIMP2)
       FragEnergies=FragEnergiesAll(:,:,FRAGMODEL_OCCLSTHCRIMP2)

    case default
       print *, 'Model is: ', ccmodel
       call lsquit('get_occfragenergies: Model needs implementation!',-1)
    end select

  end subroutine get_occfragenergies

  !> \brief Extract fragment energies for occupied partitioning for given CC model
  !> from set of all fragment energies.
  !> \author Kasper Kristensen - TK
  !> \date October 2013
  subroutine get_virtfragenergies(natoms,ccmodel,FragEnergiesAll,FragEnergies)
    implicit none
    !> Number of atoms in molecule
    integer,intent(in) :: natoms
    !> CC model according to MODEL_* conventions in dec_typedef.F90
    integer,intent(in) :: ccmodel
    !> Fragment energies for all models (see FRAGMODEL_* in dec_typedef.F90)
    real(realk),intent(in) :: FragEnergiesAll(natoms,natoms,ndecenergies)
    !> Fragment energies for occupied partitioning for given CC model
    real(realk),intent(inout) :: FragEnergies(natoms,natoms)

    ! MODIFY FOR NEW MODEL
    select case(ccmodel)
    case(MODEL_MP2)
       FragEnergies=FragEnergiesAll(:,:,FRAGMODEL_VIRTMP2)

    case(MODEL_CC2)
       FragEnergies=FragEnergiesAll(:,:,FRAGMODEL_VIRTCC2)

    case(MODEL_RPA)

       FragEnergies=FragEnergiesAll(:,:,FRAGMODEL_VIRTRPA)

    case(MODEL_SOSEX)

       FragEnergies=FragEnergiesAll(:,:,FRAGMODEL_VIRTSOS)

    case(MODEL_CCSD)
       FragEnergies=FragEnergiesAll(:,:,FRAGMODEL_VIRTCCSD)

    case(MODEL_CCSDpT)
       ! CCSD(T): Add CCSD and (T) contributions using occupied partitioning
       FragEnergies=FragEnergiesAll(:,:,FRAGMODEL_VIRTCCSD) &
            & + FragEnergiesAll(:,:,FRAGMODEL_VIRTpT)

    case(MODEL_RIMP2)
       FragEnergies=FragEnergiesAll(:,:,FRAGMODEL_VIRTRIMP2)

    case(MODEL_LSTHCRIMP2)
       FragEnergies=FragEnergiesAll(:,:,FRAGMODEL_VIRTLSTHCRIMP2)

    case default
       print *, 'Model is: ', ccmodel
       call lsquit('get_virtfragenergies: Model needs implementation!',-1)
    end select

  end subroutine get_virtfragenergies

  !> \brief Estimate (absolute) energy error in DEC calculation.
  !> \author Kasper Kristensen
  !> \date October 2013
  subroutine get_estimated_energy_error(nfrags,energies,dE_est1,dE_est2,doSOS)
    implicit none
    !> Number of fragments in molecule
    integer,intent(in) :: nfrags
    !> SUM of fragment energies for all models (see FRAGMODEL_* in dec_typedef.F90)
    real(realk),intent(in) :: energies(ndecenergies)
    !> Estimated (absolute) energy error
    real(realk),intent(inout) :: dE_est1,dE_est2
    logical,intent(in),optional :: doSOS
    real(realk) :: Eocc,Evirt
    logical :: SOS

    SOS = .false.
    if(present(doSOS)) SOS=doSOS
    ! MODIFY FOR NEW MODEL
    select case(DECinfo%ccmodel)
    case(MODEL_MP2)
       ! Energy error = max difference between occ,virt, and Lag energies
       dE_est1 = max(energies(FRAGMODEL_LAGMP2),energies(FRAGMODEL_OCCMP2),energies(FRAGMODEL_VIRTMP2)) &
            & - min(energies(FRAGMODEL_LAGMP2),energies(FRAGMODEL_OCCMP2),energies(FRAGMODEL_VIRTMP2))

    case(MODEL_CC2)
       ! Energy error = difference between occ and virt energies
       dE_est1 = abs(energies(FRAGMODEL_OCCCC2) - energies(FRAGMODEL_VIRTCC2))

    case(MODEL_RPA)
       ! Energy error = difference between occ and virt energies
       if(SOS) then
         dE_est1 = abs(energies(FRAGMODEL_OCCSOS) - energies(FRAGMODEL_VIRTSOS))
       else
         dE_est1 = abs(energies(FRAGMODEL_OCCRPA) - energies(FRAGMODEL_VIRTRPA))
       endif

    case(MODEL_SOSEX)
         ! Energy error = difference between occ and virt energies
         dE_est1 = abs(energies(FRAGMODEL_OCCSOS) - energies(FRAGMODEL_VIRTSOS))

    case(MODEL_CCSD)
       dE_est1 = abs(energies(FRAGMODEL_OCCCCSD) - energies(FRAGMODEL_VIRTCCSD))

    case(MODEL_CCSDpT)
       ! CCSD(T): Add CCSD and (T) contributions and find diff
       Eocc = energies(FRAGMODEL_OCCCCSD) + energies(FRAGMODEL_OCCpT)
       Evirt = energies(FRAGMODEL_VIRTCCSD) + energies(FRAGMODEL_VIRTpT)
       dE_est1 = abs(Eocc-Evirt)

    case(MODEL_RIMP2)
       ! Energy error = difference between occ and virt energies
       dE_est1 = max(energies(FRAGMODEL_LAGRIMP2),energies(FRAGMODEL_OCCRIMP2),energies(FRAGMODEL_VIRTRIMP2)) &
            & - min(energies(FRAGMODEL_LAGRIMP2),energies(FRAGMODEL_OCCRIMP2),energies(FRAGMODEL_VIRTRIMP2))

    case(MODEL_LSTHCRIMP2)
       ! Energy error = difference between occ and virt energies
       dE_est1 = max(energies(FRAGMODEL_LAGLSTHCRIMP2),energies(FRAGMODEL_OCCLSTHCRIMP2),energies(FRAGMODEL_VIRTLSTHCRIMP2)) &
            & - min(energies(FRAGMODEL_LAGLSTHCRIMP2),energies(FRAGMODEL_OCCLSTHCRIMP2),energies(FRAGMODEL_VIRTLSTHCRIMP2))

    case default
       print *, 'Model is: ', DECinfo%ccmodel
       call lsquit('get_estimated_energy_error: Model needs implementation!',-1)
    end select

    dE_est2 = 2.0E0_realk * nfrags * DECinfo%FOT

  end subroutine get_estimated_energy_error

  !> \brief Get filename for for fragment energy restart file 
  !> \author Kasper Kristensen
  !> \date October 2013
  function get_fragenergy_restart_filename(esti) result(filename)
    implicit none
    !> Is this estimated fragment energies?
    logical,intent(in) :: esti
    character(len=40) :: FileName

    if(esti) then
       FileName='estimated_fragenergies.info'
    else
       FileName='fragenergies.info'
    end if

  end function get_fragenergy_restart_filename

  !> \brief Get filename for for fragment energy restart file used for backing up existing file
  !> \author Kasper Kristensen
  !> \date October 2013
  function get_fragenergy_restart_filename_backup(esti) result(filename)
    implicit none
    !> Is this estimated fragment energies?
    logical,intent(in) :: esti
    character(len=40) :: FileName

    if(esti) then
       FileName='estimated_fragenergies.backup'
    else
       FileName='fragenergies.backup'
    end if

  end function get_fragenergy_restart_filename_backup


  !> \brief Get number of pair fragments with interatomic distances in the range [r1,r2].
  !> This is done by counting the number of entries in MyMolecule%ccmodel(i,j), where j>i 
  !> AND ccmodel(i,j) is not an empty model (MODEL_NONE) AND atoms "i" and "j"
  !> both have orbitals assigned AND the distance between "i" and "j" is in the range [r1,r2].
  function get_num_of_pair_fragments(MyMolecule,dofrag,r1,r2) result(npairfrags)
    implicit none
    !> Full molecule structure
    type(fullmolecule),intent(in) :: MyMolecule
    !> List of which atoms have orbitals assigned
    logical,dimension(MyMolecule%natoms),intent(in) :: dofrag
    !> Minimum distance to consider
    real(realk),intent(in) :: r1
    !> Maximum distance to consider
    real(realk),intent(in) :: r2
    integer :: npairfrags
    integer :: i,j

    npairfrags=0
    iloop: do i=1,MyMolecule%natoms
       if(.not. dofrag(i)) cycle iloop
       jloop: do j=i+1,MyMolecule%natoms
          if(.not. dofrag(j)) cycle jloop
          if(MyMolecule%ccmodel(i,j) /= MODEL_NONE) then
             if( MyMolecule%DistanceTable(i,j) .ge. r1 &
                  & .and. MyMolecule%DistanceTable(i,j) .le. r2 ) then
                ! Pair fragment (i,j) needs to be calculated
                npairfrags = npairfrags+1
             end if
          end if
       end do jloop
    end do iloop

  end function get_num_of_pair_fragments


  !> \brief Make pointers for MO coefficients and MO Fock matrices in decfrag type
  !> point to the fragment-adapted orbitals
  !> \author Kasper Kristensen
  !> \date November 2013
  subroutine fragment_basis_point_to_FOs(MyFragment,skipfock)
    implicit none
    !> Atomic or pair fragment
    type(decfrag),intent(inout) :: MyFragment
    !> Do not change pointers for Fock matrix
    logical,intent(in),optional :: skipfock
    logical :: skipf

    skipf=.false.
    if(present(skipfock)) then
       if(skipfock) skipf=.true.
    end if

    ! Sanity check: Fragment-adapted MO coefficients have been set
    if(.not. MyFragment%FAset) then
       call lsquit('fragment_basis_point_to_FOs: Fragment-adapted MO coefficients &
            & have not been set!',-1)
    end if

    ! Dimensions for fragment-adapted orbitals
    MyFragment%noccAOS => MyFragment%noccFA
    MyFragment%nvirtAOS => MyFragment%nvirtFA

    ! Total number of occupied orbitals
    if(DECinfo%frozencore) then
       ! core + valence (AOS)
       Myfragment%nocctot = Myfragment%ncore + Myfragment%noccAOS
    else
       ! AOS already contains core orbitals
       Myfragment%nocctot = Myfragment%noccAOS
    end if

    ! MO coefficients
    MyFragment%Co => MyFragment%CoFA
    MyFragment%Cv => MyFragment%CvFA

    ! MO Fock
    if(.not. skipf) then
       MyFragment%ppfock => MyFragment%ppfockFA
       MyFragment%qqfock => MyFragment%qqfockFA    
    end if

    ! Pointers point to FO data
    MyFragment%fragmentadapted=.true.

  end subroutine fragment_basis_point_to_FOs



  !> \brief Make pointers for MO coefficients and MO Fock matrices in decfrag type
  !> point to the local orbitals
  !> \author Kasper Kristensen
  !> \date November 2013
  subroutine fragment_basis_point_to_LOs(MyFragment)
    implicit none
    !> Atomic or pair fragment
    type(decfrag),intent(inout) :: MyFragment

    ! Dimensions for fragment-adapted orbitals
    MyFragment%noccAOS => MyFragment%noccLOC
    MyFragment%nvirtAOS => MyFragment%nvirtLOC

    ! Total number of occupied orbitals
    if(DECinfo%frozencore) then
       ! core + valence (AOS)
       Myfragment%nocctot = Myfragment%ncore + Myfragment%noccAOS
    else
       ! AOS already contains core orbitals
       Myfragment%nocctot = Myfragment%noccAOS
    end if

    ! MO coefficients
    MyFragment%Co => MyFragment%CoLOC
    MyFragment%Cv => MyFragment%CvLOC

    ! MO Fock
    MyFragment%ppfock => MyFragment%ppfockLOC
    MyFragment%qqfock => MyFragment%qqfockLOC   

    ! Pointers do not point to FO data
    MyFragment%fragmentadapted=.false.

  end subroutine fragment_basis_point_to_LOs


  !> \brief Initialize pointers in fragment structure handling some dimensions.
  !> These pointers are the only pointers that are not arrays.
  !> \author Kasper Kristensen
  !> \date November 2013
  subroutine fragment_init_dimension_pointers(fragment)
    implicit none
    !> Atomic or pair fragment
    type(decfrag),intent(inout) :: fragment

    nullify(fragment%noccAOS)
    nullify(fragment%nvirtAOS)
    nullify(fragment%noccLOC)
    nullify(fragment%nvirtLOC)
    nullify(fragment%noccFA)
    nullify(fragment%nvirtFA)
    allocate(fragment%noccLOC)
    allocate(fragment%nvirtLOC)
    allocate(fragment%noccFA)
    allocate(fragment%nvirtFA)

  end subroutine fragment_init_dimension_pointers

  !> \brief Calculate combined single+double amplitudes:
  !> u(a,i,b,j) = t2(a,i,b,j) + t1(a,i)*t1(b,j)
  !> \author Kasper Kristensen
  !> \date January 2012
  subroutine get_combined_SingleDouble_amplitudes_oldarr(t1,t2,u)
     implicit none
     !> Singles amplitudes t1(a,i)
     type(array2),intent(in) :: t1
     !> Doubles amplitudes t2(a,i,b,j)
     type(array4),intent(in) :: t2
     !> Combined single+double amplitudes
     type(array4),intent(inout) :: u
     integer :: i,j,a,b,nocc,nvirt

     ! Number of occupied/virtual orbitals assuming index ordering given above
     nocc=t1%dims(2)
     nvirt=t1%dims(1)

     ! Init combined amplitudes
     u = array4_init(t2%dims)

     if(DECinfo%use_singles)then
        do j=1,nocc
           do b=1,nvirt
              do i=1,nocc
                 do a=1,nvirt
                    u%val(a,i,b,j) = t2%val(a,i,b,j) + t1%val(a,i)*t1%val(b,j)
                 end do
              end do
           end do
        end do
     else
#ifdef VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN
        call assign_in_subblocks(u%val,'=',t2%val,t2%nelements)
#else
        !$OMP WORKSHARE
        u%val = t2%val
        !$OMP END WORKSHARE
#endif
     endif


  end subroutine get_combined_SingleDouble_amplitudes_oldarr

  !> \brief Calculate combined single+double amplitudes:
  !> u(a,i,b,j) = t2(a,i,b,j) + t1(a,i)*t1(b,j)
  !> \author Patrick Ettenhuber adapted from Kasper Kristensen
  !> \date January 2014
  subroutine get_combined_SingleDouble_amplitudes_newarr(t1,t2,u)
     implicit none
     !> Singles amplitudes t1(a,i)
     type(tensor),intent(in) :: t1
     !> Doubles amplitudes t2(a,i,b,j)
     type(tensor),intent(in) :: t2
     !> Combined single+double amplitudes
     type(tensor),intent(inout) :: u
     integer :: i,j,a,b,nocc,nvirt

     logical :: bg

     bg = mem_is_background_buf_init()

     ! Number of occupied/virtual orbitals assuming index ordering given above
     nocc  = t2%dims(2)
     nvirt = t2%dims(1)

     select case(t2%itype)
     case(TT_DENSE,TT_REPLICATED)

        ! Init combined amplitudes
        call tensor_init(u,t2%dims,4,bg=bg)

        if(DECinfo%use_singles)then

           do j=1,nocc
              do b=1,nvirt
                 do i=1,nocc
                    do a=1,nvirt
                       u%elm4(a,i,b,j) = t2%elm4(a,i,b,j) + t1%elm2(a,i)*t1%elm2(b,j)
                    end do
                 end do
              end do
           end do

        else

#ifdef VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN
           call assign_in_subblocks(u%elm1,'=',t2%elm1,t2%nelms)
#else
           !$OMP WORKSHARE
           u%elm1 = t2%elm1
           !$OMP END WORKSHARE
#endif

        endif
     case( TT_TILED_DIST )


        call tensor_init(u, t2%dims, int(t2%mode), tensor_type = t2%itype,&
           &pdm = t2%access_type, tdims = int(t2%tdim), fo = int(t2%offset), bg=bg )

        if(DECinfo%use_singles)then

          !This was put outside of this test, it did not make sense
          !as t1 may not have been initialized.
          if(t1%itype /= TT_DENSE .and. t1%itype /= TT_REPLICATED)then
            call lsquit("ERROR(get_combined_SingleDouble_amplitudes_newarr): only dense and replicated t1 implemented",-1)
          endif

          call lspdm_get_combined_SingleDouble_amplitudes( t1, t2, u )
        else
           !this is just copying t2 to u
           call tensor_add( u, 1.0E0_realk, t2, a = 0.0E0_realk)
        endif

     case default
        call lsquit("ERROR(get_combined_SingleDouble_amplitudes_newarr) no PDM version implemented yet",-1)
     end select


  end subroutine get_combined_SingleDouble_amplitudes_newarr

  !> Assuming that the last index in the array A contains core+valence indices,
  !> remove the core indices. Only to be used for frozen core approx.
  !> Assumes core indices are placed BEFORE valence indices!
  !> \author Kasper Kristensen
  !> \date December 2012
  subroutine remove_core_orbitals_from_last_index_oldarr(MyFragment,A,B)
    implicit none
    !> Atomic fragment
    type(decfrag),intent(inout) :: MyFragment
    !> Original array
    type(array4),intent(in) :: A
    !> New array where core indices for the last index are removed
    type(array4),intent(inout) :: B
    integer :: dims(4), i,j,k,l

    ! Sanity check 1: Frozen core.
    if(.not. DECinfo%frozencore) then
       call lsquit('remove_core_orbitals_from_last_index: Only works for frozen core approx!',-1)
    end if

    ! Sanity check 2: Correct dimensions
    if(A%dims(4) /= MyFragment%nocctot) then
       print *, 'Array dim, #occ orbitals', A%dims(4), MyFragment%nocctot
       call lsquit('remove_core_orbitals_from_last_index: Dimension mismatch!',-1)
    end if

    ! Init with new dimensions - same as before, except for last index which is only valence indices
    dims(1) = A%dims(1)
    dims(2) = A%dims(2)
    dims(3) = A%dims(3)
    dims(4) = MyFragment%noccAOS
    B = array4_init_standard(dims) 


    ! Copy elements from A to B, but only valence for last index
    do l=1,B%dims(4)
       do k=1,B%dims(3)
          do j=1,B%dims(2)
             do i=1,B%dims(1)
                B%val(i,j,k,l) = A%val(i,j,k,l+MyFragment%ncore)
             end do
          end do
       end do
    end do


  end subroutine remove_core_orbitals_from_last_index_oldarr
  !> Assuming that the last index in the array A contains core+valence indices,
  !> remove the core indices. Only to be used for frozen core approx.
  !> Assumes core indices are placed BEFORE valence indices!
  !> \author Kasper Kristensen
  !> \date December 2012
  subroutine remove_core_orbitals_from_last_index_newarr(MyFragment,A,B)
    implicit none
    !> Atomic fragment
    type(decfrag),intent(inout) :: MyFragment
    !> Original array
    type(tensor),intent(in) :: A
    !> New array where core indices for the last index are removed
    type(tensor),intent(inout) :: B
    integer :: dims(4), i,j,k,l

    ! Sanity check 1: Frozen core.
    if(.not. DECinfo%frozencore) then
       call lsquit('remove_core_orbitals_from_last_index: Only works for frozen core approx!',-1)
    end if

    ! Sanity check 2: Correct dimensions
    if(A%dims(4) /= MyFragment%nocctot) then
       print *, 'Array dim, #occ orbitals', A%dims(4), MyFragment%nocctot
       call lsquit('remove_core_orbitals_from_last_index: Dimension mismatch!',-1)
    end if

    ! Init with new dimensions - same as before, except for last index which is only valence indices
    dims(1) = A%dims(1)
    dims(2) = A%dims(2)
    dims(3) = A%dims(3)
    dims(4) = MyFragment%noccAOS

    if( A%itype == TT_DENSE )then
       call tensor_init(B, dims,4) 

       ! Copy elements from A to B, but only valence for last index
       do l=1,B%dims(4)
           B%elm4(:,:,:,l) = A%elm4(:,:,:,l+MyFragment%ncore)
       end do
    else
       call lsquit("ERROR(remove_core_orbitals_from_last_index_newarr)NO PDM VERSION",-1)
    endif


 end subroutine remove_core_orbitals_from_last_index_newarr


 
 !> Set logical arrays according to secondary assignment
 !>
 !> SEC_occ(i) = .true.   
 !>
 !> if the secondary atom of occ orbital "i" is identical to the central atom
 !> defining the atomic fragments (or to one of the central atoms for pair fragments).
 !> Otherwise SEC_occ=.false.  (and similarly for SEC_virt).
 subroutine secondary_assigning(MyFragment,SEC_occ,SEC_virt)

   implicit none

   !> fragment info
   type(decfrag), intent(inout) :: MyFragment
   !> Logical arrays defined as described above
   logical,intent(inout) :: SEC_occ(MyFragment%noccAOS), SEC_virt(MyFragment%nvirtAOS)
   integer :: noccAOS,nvirtAOS,i,P


   noccAOS = MyFragment%noccAOS
   nvirtAOS = MyFragment%nvirtAOS
   SEC_occ=.false.
   SEC_virt=.false.
   do i=1,noccAOS
      ! Set SEC_occ(i) to true if secondary atom equals (one of the) atom(s) 
      ! defining atomic (pair) fragment
      do P=1,MyFragment%nEOSatoms
         if(MyFragment%occAOSorb(i)%secondaryatom == MyFragment%EOSatoms(P)) then
            SEC_occ(i)=.true.
         end if
      end do
   end do

   ! Same for virt orbitals
   do i=1,nvirtAOS
      do P=1,MyFragment%nEOSatoms
         if(MyFragment%virtAOSorb(i)%secondaryatom == MyFragment%EOSatoms(P)) then
            SEC_virt(i)=.true.
         end if
      end do
   end do


 end subroutine secondary_assigning



 !> For two sets of points in space, make table with distances between the points of the two sets
 subroutine general_distance_table(n1,n2,list1,list2,DistanceTable)
   implicit none
   !> Dimensions of the two lists
   integer,intent(in) :: n1,n2
   !> The two lists, e.g. list1(i,j) is the x (i=1), y (i=2), or z (i=3) coordinate
   !> of the jth point in list1.
   real(realk),intent(in) :: list1(3,n1), list2(3,n2)
   !> Distance table described above
   real(realk),intent(inout) :: DistanceTable(n1,n2)
   integer :: i,j

   do j=1,n2
      do i=1,n1
           DistanceTable(i,j) = get_distance_between_two_points( list1(1:3,i) , list2(1:3,j) )
      end do
   end do
   
 end subroutine general_distance_table


 !> \brief Get distance between atoms and orbitals.
 subroutine GetOrbAtomDistances(norb,natoms,&
      & Carmom,AtomCenters,DistanceTableOrbAtom)
   implicit none
   !> Number of orbitals and atoms
   integer,intent(in) :: norb,natoms
   !> Position of atoms
   real(realk),intent(in) :: AtomCenters(3,nAtoms)
   !> Positions of orbitals
   real(realk),intent(in) :: Carmom(3,norb)
   !> Distances between orbitals and atoms
   real(realk),intent(inout) :: DistanceTableOrbAtom(norb,natoms)
   integer :: iatom,i
   real(realk) :: Xa,Ya,Za

   do iatom=1,nAtoms
      Xa = -AtomCenters(1,iatom)
      Ya = -AtomCenters(2,iatom)
      Za = -AtomCenters(3,iatom)
      do i=1,norb
         DistanceTableOrbAtom(i,iatom)= sqrt( (Xa+Carmom(1,i))*(Xa+Carmom(1,i))  &
              & + (Ya+Carmom(2,i))*(Ya+Carmom(2,i)) + (Za+Carmom(3,i))*(Za+Carmom(3,i)) )
      end do
   end do

 end subroutine GetOrbAtomDistances


  !> \brief The actual writing of job list and fragment energies
  !> \author Kasper Kristensen
  !> \date May 2012
 subroutine basic_write_jobs_and_fragment_energies_for_restart(nfrags,FragEnergies,jobs,funit,filename)

   implicit none
   !> Number of fragments
   integer,intent(in) :: nfrags
   !> Fragment energies (see decfrag type def)
   real(realk),dimension(nfrags,nfrags,ndecenergies),intent(in) :: FragEnergies
   !> Job list of fragment jobs
   type(joblist),intent(in) :: jobs
   !> File unit number
   integer,intent(in) :: funit
   !> File name
   character(len=40) :: FileName
   integer :: i,j
   logical :: file_exist
   integer(8) :: ndecenergies_file

   call write_fragment_joblist_to_file(jobs,funit)
   ndecenergies_file = ndecenergies
   write(funit) ndecenergies_file

   do j=1,ndecenergies
      do i=1,nfrags
         write(funit) FragEnergies(:,i,j)
         flush(funit)
      end do
   end do

 end subroutine basic_write_jobs_and_fragment_energies_for_restart


  !> \brief The actual reading of fragment energies for fragment from file 
  !> (assumes file is already opened)
  !> \author Kasper Kristensen
  !> \date May 2012
  subroutine basic_read_jobs_and_fragment_energies_for_restart(nfrags,FragEnergies,jobs,funit,FileName)

    implicit none
    !> Number of fragments
    integer,intent(in) :: nfrags
    !> Fragment energies (see decfrag type def)
    real(realk),dimension(nfrags,nfrags,ndecenergies),intent(inout) :: FragEnergies
    !> Job list of fragments
    type(joblist),intent(inout) :: jobs
    !> File unit number
    integer,intent(in) :: funit
    !> File name
    character(len=40),intent(in) :: FileName
    integer :: i,j
    integer(8) :: ndecenergies_file

    ! Read job list and fragment energies from file
    call read_fragment_joblist_from_file(jobs,funit)
    read(funit) ndecenergies_file

    ! Sanity check
    if(ndecenergies_file /= ndecenergies) then

       call restart_sanity_check(ndecenergies_file)

       ! Zero all fragenergies, just in case
       do j=1,ndecenergies
          do i=1,nfrags
             FragEnergies(:,i,j) = 0.0_realk
          end do
       end do

    end if

    do j=1,ndecenergies_file
       do i=1,nfrags
          read(funit) FragEnergies(:,i,j)
       end do
    end do

  end subroutine basic_read_jobs_and_fragment_energies_for_restart



  !> \brief Sanity check for reading fragment energies
  !> \author Kasper Kristensen
  !> \date May 2012
  subroutine restart_sanity_check(ndecenergies_file)

    implicit none
    !> Number of DEC energies read from file
    integer(8),intent(in) :: ndecenergies_file

    write(DECinfo%output,*) 'WARNING! Different ndecenergies: file/current:', &
         & ndecenergies_file,ndecenergies

    ! Sanity check - do we really want to restart?
    if( (ndecenergies_file > ndecenergies) .or. (ndecenergies_file<1) ) then
       call lsquit('restart_sanity_check: bad ndecenergies read from file!',-1)
    end if
    if(.not. DECinfo%EnforceRestart) then
       write(DECinfo%output,*) 'Quitting due to wrong ndecenergies read from file!'
       write(DECinfo%output,*) 'You can enforce restart using .ENFORCERESTART keyword'
       write(DECinfo%output,*) '--> USE AT OWN RISK!'
       call lsquit('Wrong ndecenergies read from file!',-1)
    end if

  end subroutine restart_sanity_check


  !> For two tensors b1 and b2, calculate sum of the following dot products:
  !> SD_ddot = dotproduct(b1,b1) + dotproduct(b2,b2)
  !> Intended to be used when b1/b2 decsribes singles/doubles quantities, but
  !> can also be used for other purposes.
  !> \author Kasper Kristensen
  !> \date June 2015
  function SD_dotproduct(b1,b2) result(SD_ddot)
    implicit none
    !> Result as described above
    real(realk) :: SD_ddot
    !> Input tensors may be of different dimensions
    type(tensor),intent(in) :: b1,b2
    real(realk) :: ddot1, ddot2 

    call print_norm(b1,nrm=ddot1,returnsquared=.true.)
    call print_norm(b2,nrm=ddot2,returnsquared=.true.)
    SD_ddot = ddot1+ddot2

  end function SD_dotproduct


end module dec_fragment_utils
