!> @file
!> Utils for DEC subroutines
!> \author Marcin Ziolkowski (modified by Kasper Kristensen)
module dec_fragment_utils

  use precision
  use ls_util!,only: dgemm_ts
  use typedeftype!,only: lsitem
  use molecule_module!, only: get_geometry
  use files!,only:lsopen,lsclose
  use DALTONINFO!, only: ls_free
  use dec_typedef_module
  use memory_handling!, only: mem_alloc, mem_dealloc, mem_allocated_global,&
!       & stats_mem, get_avaiLable_memory
  use,intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
  use matrix_module!, only:matrix
  use matrix_operations
  use IntegralInterfaceMOD!, only: ii_get_h1, ii_get_nucpot
  use BUILDAOBATCH
#ifdef VAR_LSMPI
  use infpar_module
#endif

#ifdef VAR_PAPI
    use papi_module
#endif

    ! DEC DEPENDENCIES (within deccc directory)
    ! *****************************************

  !> Maximum number of files to be opened at the same time
  ! NOTE: If you change this number, you have to change
  ! max_file accordingly in crayio.c.
  integer, parameter :: max_number_files=250
  !> Number of opened files using C file handling
  integer,save :: files_opened=0
  !> Keeping control of which file units are available
  logical,save :: available_file_units(max_number_files)=.true.

  !> Read 64 bit integer(s) from file and convert to 32 bit
  interface read_64bit_to_32bit
     module procedure read_64bit_to_32bit_singleinteger
     module procedure read_64bit_to_32bit_vectorinteger
     module procedure read_64bit_to_32bit_singlelogical
     module procedure read_64bit_to_32bit_vectorlogical
  end interface
  interface read_32bit_to_64bit
     module procedure read_32bit_to_64bit_singleinteger
     module procedure read_32bit_to_64bit_vectorinteger
     module procedure read_32bit_to_64bit_singlelogical
     module procedure read_32bit_to_64bit_vectorlogical
  end interface

contains

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

  !> \brief Sort first 'n' elements of vector a with 'm' elements
  subroutine int_sort(a,n,m)

    implicit none
    integer, dimension(m), intent(inout) :: a
    integer, intent(in) :: m,n
    integer :: i
    integer :: tmp
    logical :: swp

    swp=.true.
    do while (swp)
       swp=.false.
       do i=1,n-1
          if(a(i)>a(i+1)) then

             tmp=a(i+1)
             a(i+1)=a(i)
             a(i)=tmp

             swp=.true.
          endif
       end do
    end do

    return
  end subroutine int_sort


  !> \brief Sort real vector, keeping track of the original indices.
  !> Note: Largest elements first!
  subroutine real_inv_sort_with_tracking(to_sort,to_track,n)

    implicit none
    integer, intent(in) :: n
    real(realk), dimension(n), intent(inout) :: to_sort
    integer, dimension(n), intent(inout) :: to_track
    real(realk) :: tmp
    integer :: tmp1,i
    logical :: swp

    ! Set original track order
    do i=1,n
       to_track(i)=i
    end do

    swp=.true.
    do while (swp)
       swp=.false.
       do i=1,n-1
          if(to_sort(i) < to_sort(i+1)) then ! reverse order

             tmp = to_sort(i+1)
             to_sort(i+1) = to_sort(i)
             to_sort(i) = tmp

             tmp1 = to_track(i+1)
             to_track(i+1) = to_track(i)
             to_track(i) = tmp1

             swp=.true.
          end if
       end do
    end do
    return
  end subroutine real_inv_sort_with_tracking


  !> \brief Sort integer vector, keeping track of the original indices.
  !> Note: Largest elements first!
  !> \author Kasper Kristensen (based on real_inv_sort_with_tracking)
  !> \date February 2011
  subroutine integer_inv_sort_with_tracking(to_sort,to_track,n)

    implicit none
    !> Dimension of vector to sort
    integer, intent(in) :: n
    !> Vector to sort
    integer, dimension(n), intent(inout) :: to_sort
    !> List of sorted original indices
    integer, dimension(n), intent(inout) :: to_track
    integer :: tmp,tmp1,i
    logical :: swp

    ! Set original track order
    do i=1,n
       to_track(i)=i
    end do

    swp=.true.
    do while (swp)
       swp=.false.
       do i=1,n-1
          if(to_sort(i) < to_sort(i+1)) then ! reverse order

             tmp = to_sort(i+1)
             to_sort(i+1) = to_sort(i)
             to_sort(i) = tmp

             tmp1 = to_track(i+1)
             to_track(i+1) = to_track(i)
             to_track(i) = tmp1

             swp=.true.
          end if
       end do
    end do

  end subroutine integer_inv_sort_with_tracking


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

  !> Subroutine to find initial fragment based on a radius around MyAtom
  !> author: Ida-Marie Hoeyvik
  !> EOSvector:logical vector that controls EOS space
  !> BufferVector: logical vector that controls buffer space
  subroutine initial_fragment(MyAtom, SortedDistTable,TrackVec,&
       &EOSvector,BufferVector, natoms,counter)
    implicit none
    integer, intent(in) 	:: MyAtom, natoms
    logical, intent(inout)	:: EOSvector(natoms)
    logical                     :: BufferVector(natoms)
    real(realk)			:: SortedDistTable(natoms,natoms)
    integer			:: i
    integer			:: counter(natoms)
    real(realk)         	:: init_radius
    integer, intent(in)		:: TrackVec(natoms,natoms)

    if (DECinfo%FOT > 5E-3_realk) then
       init_radius = 5.0
    else
       init_radius = 7.0
    end if

    !>initialize counter
    do i=1,natoms
       counter(i)=1
    end do

    !>initialize logical vectors
    EOSvector(:)=.false.
    BufferVector(:)=.false.

    BufferVector(MyAtom)=.true.
    EOSvector(MyAtom) = .true.

    counter(MyAtom) = 0
    do i=1, natoms
       if (SortedDistTable(i,MyAtom) .le. init_radius) then
          EOSvector(TrackVec(i,MyAtom)) = .true.
          counter(MyAtom) = counter(MyAtom) + 1
       end if
    end do

  end subroutine initial_fragment

  !>Subroutine that takes the logicalvector containg information on what atoms
  !>are included and make a compressedlist of these
  !>Author: Ida-Marie Hoeyvik
  subroutine atoms_included(UnoccEOS_atoms,atoms_in_vec, list_of_atoms, natoms)
    implicit none
    integer, intent(in) :: natoms, atoms_in_vec
    logical, intent(in) :: UnoccEOS_atoms(natoms)
    integer, intent(out):: list_of_atoms(atoms_in_vec)
    integer		::i,j

    j=0
    do i=1,natoms
       If (UnoccEOS_atoms(i)) then
          j=j+1
          list_of_atoms(j)=i
       end if
    end do

  end subroutine atoms_included



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

  !> \brief Adjust full molecular basis to a list of orbitals (MO matrix)
  subroutine adjust_square_matrix_mo(FullMatrix,SmallMatrix,idx,norb,norb_small)
    implicit none
    real(realk), dimension(norb,norb), intent(in) :: FullMatrix
    real(realk), dimension(norb_small,norb_small), intent(inout) :: SmallMatrix
    integer, intent(in) :: norb, norb_small
    integer, dimension(norb_small), intent(in) :: idx
    integer :: i,j,idx_i,idx_j

    do i=1,norb_small
       do j=1,norb_small
          idx_i = idx(i)
          idx_j = idx(j)
          SmallMatrix(i,j) = FullMatrix(idx_i,idx_j)
       end do
    end do

    return
  end subroutine adjust_square_matrix_mo

  !> \brief Get a table with interatomic distances
  subroutine GetDistances(DistanceTable,nAtoms,mylsitem,int_output)

    implicit none
    type(lsitem), intent(inout) :: mylsitem
    integer, intent(in) :: nAtoms,int_output
    real(realk), dimension(nAtoms,nAtoms), intent(inout) :: DistanceTable
    real(realk), pointer :: geometry(:,:)
    real(realk) :: dist
    integer :: i,j,k

    DistanceTable=0.0E0_realk

    ! get geometry
    call mem_alloc(geometry,nAtoms,3)
    geometry=0.0E0_realk
    call get_geometry(int_output,0,mylsitem%input%molecule,nAtoms,geometry(:,1), &
         geometry(:,2),geometry(:,3))

    do i=1,nAtoms
       do j=1,i
          dist=0.0E0_realk
          do k=1,3
             dist=dist+(geometry(i,k)-geometry(j,k))**2
          end do
          dist=sqrt(dist)
          DistanceTable(i,j)=dist
          DistanceTable(j,i)=dist
       end do
    end do

    call mem_dealloc(geometry)
    return
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



  !> \brief Open file using C routine to be able to access arbitrary address in file.
  !> \author Kasper Kristensen
  !> \date September 2010
  subroutine openfile(funit,filename)
    implicit none
    !> File unit number
    integer, intent(in) :: funit
    !> Name of file
    character(*), intent(in) :: filename
    integer :: length,io_err,status
    logical :: is_open

    status=0
    io_err=0

    ! Length of file (avoid blank spaces)
    length = len(trim(filename))

    ! If file is opened in fortran, then close it
    inquire(file=filename,opened=is_open)
    if(is_open) close(funit,status='KEEP')

    ! Open file using C routine
    call wopen(funit,filename(1:length),length,status,io_err)

    ! Check that everything went fine
    if(io_err /= 0) then
       write(DECinfo%output,*) 'DEC openfile: Something wrong when opening file'
       write(DECinfo%output,*) 'Filename:', filename(1:length)
       write(DECinfo%output,*) 'File unit:', funit
       write(DECinfo%output,*) 'Total number of opened files:', files_opened
       CALL lsQUIT('DEC openfile: Something went wrong when opening file!',DECinfo%output)
    end if

    files_opened = files_opened+1

  end subroutine openfile


  !> \brief Find available file unit for opening file
  !> using C file handling routine.
  !> \author Kasper Kristensen
  !> \date October 2010
  subroutine get_available_file_unit(funit)
    implicit none
    !> File unit
    integer, intent(inout) :: funit
    integer :: i

    funit=0

    do i=1,max_number_files

       ! Avoid numbers 5 and 6 and output file unit.
       ! (I am quite sure this is redundant as these file units
       !  would only be problematic if the array4 file handling
       !  was done using Fortran. However, just in case,
       ! we omit these numbers here...)
       if(i==5 .or. i==6 .or. i==DECinfo%output) cycle

       if(available_file_units(i)) then ! File unit i is available
          funit=i
          exit
       end if


    end do

    ! Check that an available file unit was found
    if(funit == 0) then
       call lsquit('get_available_file_unit: &
            & No available file unit was found for opening file using &
            & C file handling. Try increasing max_number_files &
            & in dec_utils.f90.',DECinfo%output)
    end if

    ! Now funit is no longer available to open another file
    available_file_units(funit) = .false.

  end subroutine get_available_file_unit




  !> \brief Close file using C routine.
  !> \author Kasper Kristensen
  !> \date September 2010
  subroutine closefile(funit,keep_or_delete)
    implicit none
    !> File unit number
    integer, intent(in) :: funit
    !> Status for file? 'KEEP' or 'DELETE'
    character(*), intent(in) :: keep_or_delete
    integer ::io_err, status_for_file

    io_err=0

    ! The integer status_for_file is set to:
    ! 0 if the file is to be deleted
    ! 1 if the file is to be kept
    ! In this way the communication with C routine wclose is consistent
    if(keep_or_delete=='delete' .or. keep_or_delete=='DELETE') then
       status_for_file=0
    elseif(keep_or_delete=='keep' .or. keep_or_delete=='KEEP') then
       status_for_file=1
    else
       CALL lsQUIT('DEC closefile: keep_or_delete can only be KEEP or DELETE!',DECinfo%output)
    end if

    ! Close file using C routine
    call wclose(funit,io_err,status_for_file)

    ! Check that everything went fine
    if(io_err /= 0) then
       write(DECinfo%output,*) 'DEC closefile: Something wrong when closing file unit:', funit
       CALL lsQUIT('DEC closefile: Something went wrong when closing file!',DECinfo%output)
    end if

    files_opened = files_opened-1
    ! Now funit is again available to open another file
    available_file_units(funit) = .true.


  end subroutine closefile

  !> \brief Writes vector to file at given address using C routine
  !> \author Kasper Kristensen
  !> \date October 2010
  !> Note: The input is a VECTOR and elements are written in "Fortran order".
  !> I.e. if one inputs an two-dimensional real array as the "vector"
  !> then it is written to file column by column.
  subroutine writevector(funit,begin_add,nelements,vector)
    implicit none
    !> Logical unit number for file
    integer, intent(in) :: funit
    !> Begin address (where we start writing the vector elements)
    integer(kind=long), intent(in) :: begin_add
    !> Number of elements to write
    integer(kind=long), intent(in) :: nelements
    !> Real elements to write to file
    real(realk), dimension(nelements) :: vector
    integer :: io_err

    io_err=0

    ! Call C routine to write to file at begin_add
    call putwa(funit,vector,begin_add,nelements,io_err)

    ! Check that everything went fine
    if(io_err /= 0) then
       write(DECinfo%output,*) 'DEC writevector Something wrong &
            &when writing to file unit:', funit
       CALL lsQUIT('DEC writevector: Something went wrong when writing to file!',DECinfo%output)
    end if

  end subroutine writevector


  !> \brief Reads vector from file at given address using C routine
  !> \author Kasper Kristensen
  !> \date October 2010
  !> Note: The output is a VECTOR and elements are read in Fortran order.
  !> I.e. if one inputs an two-dimensional real array as the "vector" ,
  !> then the elements in the file is read into the vector column by column.
  subroutine readvector(funit,begin_add,nelements,vector)
    implicit none
    !> Logical unit number for file
    integer, intent(in) :: funit
    !> Begin address (where we start reading the vector elements)
    integer(kind=long), intent(in) :: begin_add
    !> Number of elements to read
    integer(kind=long), intent(in) :: nelements
    !> Real elements to read from file
    real(realk), dimension(nelements) :: vector
    integer :: io_err

    io_err=0

    ! Call C routine to write to file at begin_add
    call getwa(funit,vector,begin_add,nelements,io_err)

    ! Check that everything went fine
    if(io_err /= 0) then
       write(DECinfo%output,*) 'DEC readvector Something wrong &
            &when reading from file unit:', funit
       CALL lsQUIT('DEC readvector: Something went wrong &
            & when reading from file! &
            & Perhaps the file has not been opened before reading...',DECinfo%output)
    end if

  end subroutine readvector



  !> \brief Copy file using C routine.
  !> \author Kasper Kristensen
  !> \date October 2010
  subroutine copyfile(source_file, destination_file)
    implicit none
    !> Name of (old) source file to copy
    character(*), intent(in) :: source_file
    !> Name of (new) destination file
    character(*), intent(in) :: destination_file
    integer :: source_length, destination_length

    ! Length of source file (avoid blank spaces)
    source_length = len(trim(source_file))

    ! Length of destination file (avoid blank spaces)
    destination_length = len(trim(destination_file))

    ! Copy file using C routine
    call filecopy_c(source_file(1:source_length),source_length, &
         & destination_file(1:destination_length),destination_length )

  end subroutine copyfile



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



  !> \brief Solve eigenvalue problem: F*C = S*C*eival
  !> \author Kasper Kristensen
  !> \date February 2011
  subroutine solve_eigenvalue_problem(n,F,S,eival,C)
    implicit none

    !> Dimension of matrices (Number of eigenvalues)
    integer,intent(in) :: n
    !> F matrix in eigenvalue problem (typically Fock matrix)
    real(realk),intent(in) :: F(n,n)
    !> Overlap matrix
    real(realk),intent(in) :: S(n,n)
    !> Eigenvectors
    real(realk),intent(inout) :: C(n,n)
    real(realk),intent(inout) :: eival(n)
    real(realk), pointer :: tmp(:,:)

    ! We must use temporary matrix as overlap matrix input because it is overwritten
    call mem_alloc(tmp,n,n)
    tmp(1:n,1:n) = S(1:n,1:n)

    ! Copy F elements to C, then use C as "F input".
    ! At the end, C is then overwritten by the eigenvectors.
    C(1:n,1:n) = F(1:n,1:n)

    ! Solve eigenvalue problem
    call my_DSYGV(N,C,tmp,eival,"DEC_SOLVE_EIGENVALUE")

    ! Free stuff
    call mem_dealloc(tmp)

  end subroutine solve_eigenvalue_problem



  !> \brief Solve eigenvalue problem: F*C = C*eival   (overlap matrix is the unit matrix)
  !> \author Kasper Kristensen
  !> \date February 2011
  subroutine solve_eigenvalue_problem_unitoverlap(n,F,eival,C)
    implicit none

    !> Dimension of matrices (Number of eigenvalues)
    integer,intent(in) :: n
    !> F matrix in eigenvalue problem (typically Fock matrix)
    real(realk),intent(in) :: F(n,n)
    !> Eigenvectors
    real(realk),intent(inout) :: C(n,n)
    real(realk),intent(inout) :: eival(n)
    real(realk), pointer :: tmp(:,:)
    integer :: i

    ! Overlap matrix is unit matrix
    call mem_alloc(tmp,n,n)
    tmp = 0.0_realk
    do i=1,n
       tmp(i,i) = 1.0_realk
    end do

    ! Copy F elements to C, then use C as "F input".
    ! At the end, C is then overwritten by the eigenvectors.
    C(1:n,1:n) = F(1:n,1:n)

    ! Solve eigenvalue problem
    call my_DSYGV(n,C,tmp,eival,"DEC_SOLVE_EIGENVALU2")

    ! Free stuff
    call mem_dealloc(tmp)

  end subroutine solve_eigenvalue_problem_unitoverlap


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
  subroutine InitialFragment(natoms,nocc_per_atom,nunocc_per_atom,DistMyatom,Occ,Virt)
    implicit none
    !> number of atoms in MOLECULE
    Integer, intent(in)    :: natoms
    !> Number of occupied / unoccupied orbitals per atom
    integer,intent(in), dimension(natoms) :: nocc_per_atom, nunocc_per_atom
    !> Distances from central atom to other atoms
    real(realk),intent(in) :: DistMyAtom(natoms)
    !> Which AOS atoms to include for occ and virt spaces
    !> (entry i is T if orbitals on atom i is included)
    !> (In practice occ and virt will be identical but we keep it general)
    logical,intent(inout)  :: Occ(natoms),Virt(natoms)
    real(realk)            :: init_radius
    integer                :: i
    real(realk)            :: FOT

    FOT=DECinfo%FOT

    ! Larger init radius for tighter FOT
    if (FOT > 5.0e-4) then
       init_radius = 6.0 
    else
       init_radius = 8.0 
    end if

    write(DECinfo%output,'(a,f5.2)') " FOP Radius for initial fragment: ",init_radius

    ! Include atoms within init_radius
    Occ=.false.
    Virt=.false.
    do i=1,natoms

       ! Skip if no orbitals are assigned - do NOT modify this line.
       if(nocc_per_atom(i)==0 .or. nunocc_per_atom(i)==0) cycle

       if (DistMyAtom(i) .le. init_radius) then
          Occ(i) = .true.
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
    real(realk), dimension(dim1,dim2) :: tmp
    logical :: file_exist
    integer :: funit, i,j
    integer(kind=8) :: i64,j64
    integer(kind=4) :: i32,j32


    file_exist=.false.
    inquire(file=filename,exist=file_exist)
    if(file_exist) then

       funit=-1
       call lsopen(funit,filename,'OLD','UNFORMATTED')

       ! Read dimensions stored on file
       if(DECinfo%convert64to32) then
          ! file uses 64 bit integers but current run uses 32 bit integers
          read (funit) i64,j64
          i=int(i64,4)
          j=int(j64,4)
       elseif(DECinfo%convert32to64) then
          ! file uses 64 bit integers but current run uses 32 bit integers
          read (funit) i32,j32
          i=i32
          j=j32
       else
          read (funit) i,j
       end if

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

  end subroutine dec_read_mat_from_file



  !> \brief Print all elements of four-dimensional array to LSDALTON.OUT.
  !> Only to be used for testing purposes!
  !> \author Kasper Kristensen
  !> \date November 2011
  subroutine print_4_dimensional_array(dims,A,label)

    implicit none
    !> Dimensions of 4-dimensional array
    integer,dimension(4),intent(in) :: dims
    !> 4-dimensional array to be printed
    real(realk),intent(in) :: A(dims(1),dims(2),dims(3),dims(4))
    !> Label for array
    character(*), intent(in) :: label
    integer :: i,j,k,l

    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*) '***********************************************************'
    write(DECinfo%output,*) '             ARRAY LABEL: ', label
    write(DECinfo%output,*) '***********************************************************'
    write(DECinfo%output,*)
    write(DECinfo%output,'(8X,a,8X,a,8X,a,8X,a,12X,a)') 'i','j','k','l', 'value'

    do i=1,dims(1)
       do j=1,dims(2)
          do k=1,dims(3)
             do l=1,dims(4)
                write(DECinfo%output,'(4i9,5X,g18.10)') i,j,k,l,A(i,j,k,l)
             end do
          end do
       end do
    end do

    write(DECinfo%output,*)
    write(DECinfo%output,*)


  end subroutine print_4_dimensional_array



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
    type(ccatom),intent(inout) :: MyFragment
    real(realk), intent(inout) :: fragmem
    real(realk) :: O,V,A,tmp,GB


    ! Number of occupied (O), Virtual (V), atomic basis functions (A)
    ! ***************************************************************
    O = MyFragment%noccAOS
    V = MyFragment%nunoccAOS
    A = MyFragment%number_basis
    GB = 1.000E9_realk ! 1 GB

    ! Use type ccatom to calculate memory use
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
  !> NOTE!! Only start,end, and N are set here. Due to weird PGI problems
  !> the pointer association must be executed where the pointer is needed!
  !> E.g. using:
  !> CALL c_f_pointer(c_loc(arr(thepointer%start)),thepointer%p,[thepointer%N])
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
!    nullify(thepointer%p)
!    CALL c_f_pointer(c_loc(arr(thepointer%start)),thepointer%p,[thepointer%N])


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

  end subroutine start_flop_counter


  !> \brief Stop and read PAPI FLOP counter.
  !> If LSDALTON is not linked to PAPI, we just return 0.
  !> \author Kasper Kristensen
  subroutine end_flop_counter(flops)

    implicit none
    !> Flops (simplest to save as real)
    real(realk),intent(inout) :: flops
    !> "FLOPS" used by PAPI must be hardcoded 64 bit integer
    integer(kind=8) :: flops_int
    integer :: retval

flops_int=0
retval=0
#ifdef VAR_PAPI
    call PAPIf_stop(eventset,flops_int,retval)
#endif
    flops = real(flops_int)

  end subroutine end_flop_counter



  !> \brief Print fragment information to file unit
  !> referenced by DECinfo%output.
  !> \author Kasper Kristensen
  subroutine fragment_print(MyFragment,printlevel)

    implicit none
    !> Atomic (or pair) fragment
    type(ccatom),intent(inout) :: MyFragment
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
       write(DECinfo%output,'(1X,a,i6)') 'SINGLE FRAGMENT, atomic number:', MyFragment%atomic_number
    else
       pair=.true.
       write(DECinfo%output,'(1X,a,2i6)') 'PAIR FRAGMENT, atomic numbers:', &
            & MyFragment%EOSatoms(1),MyFragment%EOSatoms(2)
       write(DECinfo%output,'(1X,a,f12.2)') 'Pair distance = ', MyFragment%pairdist
    end if
    write(DECinfo%output,'(1X,a,i8)') 'Size occ EOS   =', MyFragment%noccEOS
    write(DECinfo%output,'(1X,a,i8)') 'Size virt EOS  =', MyFragment%nunoccEOS
    write(DECinfo%output,'(1X,a,i8)') 'Size occ AOS   =', MyFragment%noccAOS
    write(DECinfo%output,'(1X,a,i8)') 'Size virt AOS  =', MyFragment%nunoccAOS
    write(DECinfo%output,'(1X,a,i8)') '#core orbitals =', MyFragment%ncore
    write(DECinfo%output,'(1X,a,i8)') 'Red Frag: Size occ AOS  =', MyFragment%noccAOS
    write(DECinfo%output,'(1X,a,i8)') 'Red Frag: Size virt AOS =', MyFragment%nunoccAOS
    write(DECinfo%output,'(1X,a,i8)') 'Atoms in atomic extent =', MyFragment%number_atoms
    write(DECinfo%output,'(1X,a,i8)') 'Number of basis functions =', MyFragment%number_basis
    do i=1,ndecenergies
       write(DECinfo%output,'(1X,a,i5,f16.10)') 'Idx, Fragenergy ', i,MyFragment%energies(i)
    end do

    write(DECinfo%output,*)
    write(DECinfo%output,*) 'Fragment-adapted orbital (FO) information:'
    write(DECinfo%output,*) 'Using FOs: ', MyFragment%fragmentadapted
    write(DECinfo%output,'(1X,a,i8)') 'Size occ FOs   = ', MyFragment%noccFA
    write(DECinfo%output,'(1X,a,i8)') 'Size unocc FOs = ', MyFragment%nunoccFA
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
       do i=1,MyFragment%nunoccEOS
          write(DECinfo%output,*) i, MyFragment%unoccEOSidx(i)
       end do
       write(DECinfo%output,*)

       write(DECinfo%output,*) 'Occ AOS indices (frag,full)'
       do i=1,MyFragment%noccAOS
          write(DECinfo%output,*) i, MyFragment%occAOSidx(i)
       end do
       write(DECinfo%output,*)

       write(DECinfo%output,*) 'Virt AOS indices (frag,full)'
       do i=1,MyFragment%nunoccAOS
          write(DECinfo%output,*) i, MyFragment%unoccAOSidx(i)
       end do
       write(DECinfo%output,*)

       write(DECinfo%output,*) 'Reduced Frag: Occ AOS indices (frag,full)'
       do i=1,MyFragment%REDnoccAOS
          write(DECinfo%output,*) i, MyFragment%REDoccAOSidx(i)
       end do
       write(DECinfo%output,*)

       write(DECinfo%output,*) 'Reduced Frag: Virt AOS indices (frag,full)'
       do i=1,MyFragment%REDnunoccAOS
          write(DECinfo%output,*) i, MyFragment%REDunoccAOSidx(i)
       end do
       write(DECinfo%output,*)

       write(DECinfo%output,*) 'Atomic extent indices (frag,full in terms of atoms)'
       do i=1,MyFragment%number_atoms
          write(DECinfo%output,*) i, MyFragment%atoms_idx(i)
       end do
       write(DECinfo%output,*)

    end if

    if(printlevel > 2 .and. MyFragment%BasisInfoIsSet) then
       
       write(DECinfo%output,*) 'Occ MO coefficients (column, elements in column)'
       do i=1,MyFragment%noccAOS
          write(DECinfo%output,*) i, MyFragment%ypo(:,i)
       end do
       write(DECinfo%output,*)

       write(DECinfo%output,*) 'Virt MO coefficients (column, elements in column)'
       do i=1,MyFragment%nunoccAOS
          write(DECinfo%output,*) i, MyFragment%ypv(:,i)
       end do
       write(DECinfo%output,*)

       write(DECinfo%output,*) 'AO fock matrix (column, elements in column)'
       do i=1,MyFragment%number_basis
          write(DECinfo%output,*) i, MyFragment%fock(:,i)
       end do
       write(DECinfo%output,*)

       write(DECinfo%output,*) 'Occ-occ fock matrix (column, elements in column)'
       do i=1,MyFragment%noccAOS
          write(DECinfo%output,*) i, MyFragment%ppfock(:,i)
       end do
       write(DECinfo%output,*)

       write(DECinfo%output,*) 'Virt-virt fock matrix (column, elements in column)'
       do i=1,MyFragment%nunoccAOS
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
          do i=1,MyFragment%nunoccAOS
             write(DECinfo%output,*) i, MyFragment%VirtMat(:,i)
          end do
          write(DECinfo%output,*)
       end if

       if(MyFragment%FAset) then
          write(DECinfo%output,*) 'Occupied FO coefficients (column, elements in column)'
          do i=1,MyFragment%noccFA
             write(DECinfo%output,*) i, MyFragment%CoccFA(:,i)
          end do
          write(DECinfo%output,*)

          write(DECinfo%output,*) 'Virtual FO coefficients (column, elements in column)'
          do i=1,MyFragment%nunoccFA
             write(DECinfo%output,*) i, MyFragment%CunoccFA(:,i)
          end do
          write(DECinfo%output,*)
       end if

    end if

    write(DECinfo%output,*) '*************************************************************'
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*)



  End subroutine fragment_print



  !> \brief Delete atomic fragment
  !> \author Marcin Ziolkowski
  subroutine atomic_fragment_free(fragment)

    implicit none
    !> Atomic fragment to be freed
    type(ccatom),intent(inout) :: fragment

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
    type(ccatom),intent(inout) :: fragment
    integer :: i

    if(associated(fragment%occEOSidx)) then
       call mem_dealloc(fragment%occEOSidx)
       fragment%occEOSidx => null()
    end if

    if(associated(fragment%occAOSidx)) then
       call mem_dealloc(fragment%occAOSidx)
       fragment%occAOSidx => null()
    end if

    if(associated(fragment%unoccEOSidx)) then
       call mem_dealloc(fragment%unoccEOSidx)
       fragment%unoccEOSidx => null()
    end if

    if(associated(fragment%unoccAOSidx)) then
       call mem_dealloc(fragment%unoccAOSidx)
       fragment%unoccAOSidx => null()
    end if


    if(associated(fragment%REDoccAOSidx)) then
       call mem_dealloc(fragment%REDoccAOSidx)
       fragment%REDoccAOSidx => null()
    end if

    if(associated(fragment%REDunoccAOSidx)) then
       call mem_dealloc(fragment%REDunoccAOSidx)
       fragment%REDunoccAOSidx => null()
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
       call mem_dealloc(fragment%CoccFA)
       call mem_dealloc(fragment%CunoccFA)
       fragment%FAset=.false.
    end if

    if(associated(fragment%atoms_idx)) then
       call mem_dealloc(fragment%atoms_idx)
    end if

    if(associated(fragment%basis_idx)) then
       call mem_dealloc(fragment%basis_idx)
    end if

    ! DEC orbitals
    if(associated(fragment%occAOSorb)) then
       do i=1,size(fragment%occAOSorb)
          call orbital_free(fragment%occAOSorb(i))
       end do
       call mem_dealloc(fragment%occAOSorb)
    end if

    if(associated(fragment%unoccAOSorb)) then
       do i=1,size(fragment%unoccAOSorb)
          call orbital_free(fragment%unoccAOSorb(i))
       end do
       call mem_dealloc(fragment%unoccAOSorb)
    end if

  end subroutine atomic_fragment_free_simple



  !> \brief Delete basis set infomation for atomic fragment (fock matrix, MO coefficients, lsitem etc.)
  !> \author Kasper Kristensen
  subroutine atomic_fragment_free_basis_info(fragment)

    implicit none
    !> Atomic fragment to be freed
    type(ccatom),intent(inout) :: fragment
    integer :: i

    if(associated(fragment%S)) then
       call mem_dealloc(fragment%S)
    end if

    ! Transformation matrices
    if(associated(fragment%ypo)) then
       call mem_dealloc(fragment%ypo)
    end if

    if(associated(fragment%ypv)) then
       call mem_dealloc(fragment%ypv)
    end if

    if(associated(fragment%coreMO)) then
       call mem_dealloc(fragment%coreMO)
    end if

    if(associated(fragment%fock)) then
       call mem_dealloc(fragment%fock)
    end if

    if(associated(fragment%qqfock)) then
       call mem_dealloc(fragment%qqfock)
    end if

    if(associated(fragment%ppfock)) then
       call mem_dealloc(fragment%ppfock)
    end if

    if(associated(fragment%ccfock)) then
       call mem_dealloc(fragment%ccfock)
    end if

    ! delete dalton input
#ifdef VAR_LSMPI
    ! Quick fix such that lsitem is never handled for global master                                
    ! as this will destroy the overall MPI framework.
    if(infpar%mynum/=infpar%master) call ls_free(fragment%mylsitem)
#else
    call ls_free(fragment%mylsitem)
#endif

    ! Internal control of whether basis info is set or not
    fragment%BasisInfoIsSet=.false.

  end subroutine atomic_fragment_free_basis_info


  !> \brief Free and nullify fragment information related to t1 amplitudes.
  !> \author Kasper Kristensen
  !> \date January 2012
  subroutine free_fragment_t1(Fragment)
    implicit none
    !> Fragment info (only t1-related information will be modified here)
    type(ccatom), intent(inout) :: Fragment

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
    type(ccorbital), intent(inout) :: myorbital

    if(associated(myorbital%atoms)) then
       call mem_dealloc(myorbital%atoms)
       myorbital%atoms => null()
    end if

  end subroutine orbital_free


  !> \brief Determine memory for DEC calculation and store in DECinfo%memory.
  !> If memory was set manually in input, nothing is done here.
  !> Otherwise a system call is used to determine memory.
  !> \author Kasper Kristensen
  !> \date August 2012
  subroutine get_memory_for_dec_calculation()
    implicit none
    real(realk) :: mem
    logical :: memfound

    memfound=.false.
    if(DECinfo%memory_defined) then ! do nothing 
       write(DECinfo%output,'(1X,a,g12.4,a)') 'Memory set in input to be: ', DECinfo%memory, ' GB'

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
    
    ! Total memory was determined at the beginning of DEC calculaton (DECinfo%memory)
    ! Memory currently in use is mem_allocated_global (see memory.f90)
    ! Memory available: Total memory - memory in use
    ! Mem in use in GB
    MemInUse = 1.0E-9_realk*mem_allocated_global
    mem = DECinfo%memory - MemInUse

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



  !> \brief Get list of atoms sorted according to distance from "MyAtom".
  !> \author Ida-Marie Hoeyvik
  subroutine GetSortedList(ListMyAtom,ListTrack,ToSort,&
       & natoms,MyAtom)
    implicit none
    integer,intent(in)     :: natoms,MyAtom
    real(realk),intent(in) :: ToSort(natoms,natoms)
    real(realk)            :: ListMyAtom(natoms), TempList(natoms)
    integer                :: ListTrack(natoms), TempTrack(natoms)
    integer                :: counter,i

    ListMyAtom(:)=ToSort(:,MyAtom)
    ! Sort large--> small
    call real_inv_sort_with_tracking(ListMyAtom,ListTrack,natoms)

    TempList = 0.0E0_realk
    TempTrack = 0
    counter = 1

    ! change to small-->large
    do i=natoms,1,-1
       TempList(counter) = ListMyAtom(i)
       TempTrack(counter)= ListTrack(i)
       counter = counter + 1
    end do

    ListMyAtom=TempList
    ListTrack=TempTrack

  end subroutine GetSortedList



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
    if(npoints==1) then  ! Skip plot if only 1 point
       write(DECinfo%output,*) 'Ascii-plot will be skipped because there is only one point!'
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
    type(ccatom),intent(inout) :: Fragment1
    !>  fragment 2
    type(ccatom),intent(inout) :: Fragment2
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
  !> \author Kasper Kristensen
  !> \date November 2012
  subroutine get_density_from_occ_orbitals(nbasis,nocc,Cocc,dens)
    implicit none
    !> Number of basis functions
    integer,intent(in) :: nbasis
    !> Number of occupied orbitals (can be only valence orbitals for frozen core)
    integer,intent(in) :: nocc
    !> Occupied MO coefficients (can be only valence orbitals for frozen core)
    real(realk),intent(in),dimension(nbasis,nocc) :: Cocc
    !> Density
    real(realk),intent(inout),dimension(nbasis,nbasis) :: dens
    real(realk),pointer :: Cocc_copy(:,:)
    integer :: i,j

    ! Cocc copy (avoid passing the same element into dgemm twice)
    call mem_alloc(Cocc_copy,nbasis,nocc)
    Cocc_copy = Cocc
    
    ! density = Cocc Cocc^T 
    call dec_simple_dgemm(nbasis,nocc,nbasis,Cocc,Cocc_copy,dens,'n','t')
    call mem_dealloc(Cocc_copy)


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
    type(matrix) :: F
    
    ! Init Fock matrix in matrix form
    call mat_init(F,MyMolecule%nbasis,MyMolecule%nbasis)
    call mat_set_from_full(MyMolecule%fock, 1E0_realk, F)
    
    ! Get HF energy
    Ehf = get_HF_energy(D,F,Mylsitem) 

    call mat_free(F)

  end function get_HF_energy_fullmolecule


  !> Check whether fragment restart file exist.
  !> \author Kasper Kristensen
  !> \date December 2012
  function fragment_restart_file_exist(first_order) result(file_exist)

    implicit none
    !> First order calculation?
    logical,intent(in) :: first_order
    logical :: file_exist
    character(len=40) :: FileName

    if(first_order) then  ! first order calculation
       filename = 'mp2grad.info'

    else ! energy calculation
       filename = 'fragenergies.info'
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
    myint = myint_long

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
    type(ccatom), intent(inout) :: MyFragment
    !> Singles amplitudes to be stored (stored as virtual,occupied)
    type(array2),intent(in) :: t1
    integer :: nocc,nvirt,i,a,ix,ax

    ! Init dimensions
    nocc = MyFragment%noccAOS   ! occupied AOS dimension
    nvirt = MyFragment%nunoccAOS   ! virtual AOS dimension

    ! Sanity check
    if( (nvirt/=t1%dims(1)) .or. (nocc/=t1%dims(2)) ) then
       write(DECinfo%output,*) 'Fragment virt,occ', nvirt,nocc
       write(DECinfo%output,*) 't1 input virt,occ', t1%dims
       call lsquit('save_fragment_t1_AOSAOSamplitudes &
            & AOS dimension mismatch!',DECinfo%output)
    end if

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
    MyFragment%t1_virtidx = MyFragment%unoccAOSidx ! virtual AOS indices


    ! Save amplitudes and indices
    do i=1,nocc
       do a=1,nvirt
          MyFragment%t1(a,i) = t1%val(a,i)
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
    type(ccatom),intent(inout) :: Fragment1
    ! Fragment 2 in pair
    type(ccatom),intent(inout) :: Fragment2
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
    type(ccatom),intent(inout) :: Fragment1
    ! Fragment 2 in pair
    type(ccatom),intent(inout) :: Fragment2
    !> Pair fragment
    type(ccatom),intent(inout) :: PairFragment
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


  end subroutine which_pairs_occ



  !> \brief Construct logical array telling which unoccupied orbital pair combinations should be
  !> included when calculating pair interaction energies and other pair interaction
  !> properties - to avoid double counting.
  !> \author Kasper Kristensen
  !> \date September 2011
  subroutine which_pairs_unocc(Fragment1,Fragment2,PairFragment,dopair)

    implicit none
    ! Fragment 1 in pair
    type(ccatom),intent(inout) :: Fragment1
    ! Fragment 2 in pair
    type(ccatom),intent(inout) :: Fragment2
    !> Pair fragment
    type(ccatom),intent(inout) :: PairFragment
    !> Do pair or not - dimension: (nunoccEOS,nunoccEOS) for PAIR
    logical,dimension(PairFragment%nunoccEOS,PairFragment%nunoccEOS),&
         & intent(inout) :: dopair
    integer :: a,b,ax,bx,p1,p2,i

    ! This is the same as which_pairs_occ, but for the unoccupied space.
    ! See example in which_pairs_occ with "occ" replaced by "unocc".

    ! Set which atoms to consider for pair
    dopair=.false.
    do a=1,fragment1%nunoccEOS   ! Unoccupied EOS for fragment 1
       do b=1,fragment2%nunoccEOS ! Unoccupied EOS for fragment 2

          ax=fragment1%unoccEOSidx(a)  ! index in full list orbitals
          bx=fragment2%unoccEOSidx(b)  ! index in full list orbitals
          p1=0
          p2=0

          ! Index for fragment 1 in pair fragment list
          do i=1,PairFragment%nunoccEOS

             ! "ax" index in PairFragment%unoccEOSidx list
             if(PairFragment%unoccEOSidx(i) == ax) p1 = i

             ! "bx" index in PairFragment%unoccEOSidx list
             if(PairFragment%unoccEOSidx(i) == bx) p2 = i

          end do

          ! Sanity check
          if(p1==p2 .or. p1==0 .or. p2==0 ) then
             write(DECinfo%output,'(1X,a,4i6)') 'ax,bx,p1,p2', ax,bx,p1,p2
             call lsquit('which_pairs_unocc: &
                  & Something wrong with indices in pair',DECinfo%output)
          end if


          ! Pair interaction for (p1,p2) index pair
          dopair(p1,p2)=.true.
          dopair(p2,p1)=.true.

       end do
    end do


  end subroutine which_pairs_unocc



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

    ! Pair cutoff
    write(funit) DECinfo%pair_distance_threshold

    write(funit) jobs%njobs
    write(funit) jobs%atom1
    write(funit) jobs%atom2
    write(funit) jobs%jobsize
    write(funit) jobs%jobsdone

    ! MPI fragment statistics
    write(funit) jobs%nslaves
    write(funit) jobs%nocc
    write(funit) jobs%nvirt
    write(funit) jobs%nbasis
    write(funit) jobs%ntasks
    write(funit) jobs%flops
    write(funit) jobs%LMtime
    write(funit) jobs%load

  end subroutine write_fragment_joblist_to_file


  !> Read fragment job list from file (includes initialization of job list).
  !> Also read pair cutoff distance in case that was changed during the original calculation.
  !> \author Kasper Kristensen
  !> \date November 2012
  subroutine read_fragment_joblist_from_file(jobs,funit)
    implicit none
    !> Job list of fragments
    type(joblist),intent(inout) :: jobs
    !> File unit number to read from (of course assumes that file is open)
    integer,intent(in) :: funit
    integer :: njobs

    ! Pair cutoff
    read(funit) DECinfo%pair_distance_threshold
    write(DECinfo%output,'(1X,a,g20.8)') 'Pair cutoff read from file:', DECinfo%pair_distance_threshold

    if(DECinfo%convert64to32) then
       call read_64bit_to_32bit(funit,njobs)
       call init_joblist(njobs,jobs)
       call read_64bit_to_32bit(funit,njobs,jobs%atom1)
       call read_64bit_to_32bit(funit,njobs,jobs%atom2)
       call read_64bit_to_32bit(funit,njobs,jobs%jobsize)
       call read_64bit_to_32bit(funit,njobs,jobs%jobsdone)
       call read_64bit_to_32bit(funit,njobs,jobs%nslaves)
       call read_64bit_to_32bit(funit,njobs,jobs%nocc)
       call read_64bit_to_32bit(funit,njobs,jobs%nvirt)
       call read_64bit_to_32bit(funit,njobs,jobs%nbasis)
       call read_64bit_to_32bit(funit,njobs,jobs%ntasks)
    elseif(DECinfo%convert32to64) then
       call read_32bit_to_64bit(funit,njobs)
       call init_joblist(njobs,jobs)
       call read_32bit_to_64bit(funit,njobs,jobs%atom1)
       call read_32bit_to_64bit(funit,njobs,jobs%atom2)
       call read_32bit_to_64bit(funit,njobs,jobs%jobsize)
       call read_32bit_to_64bit(funit,njobs,jobs%jobsdone)
       call read_32bit_to_64bit(funit,njobs,jobs%nslaves)
       call read_32bit_to_64bit(funit,njobs,jobs%nocc)
       call read_32bit_to_64bit(funit,njobs,jobs%nvirt)
       call read_32bit_to_64bit(funit,njobs,jobs%nbasis)
       call read_32bit_to_64bit(funit,njobs,jobs%ntasks)
    else
       read(funit) njobs
       call init_joblist(njobs,jobs)
       read(funit) jobs%atom1
       read(funit) jobs%atom2
       read(funit) jobs%jobsize
       read(funit) jobs%jobsdone

       read(funit) jobs%nslaves
       read(funit) jobs%nocc
       read(funit) jobs%nvirt
       read(funit) jobs%nbasis
       read(funit) jobs%ntasks
    end if

    read(funit) jobs%flops
    read(funit) jobs%LMtime
    read(funit) jobs%load

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
    nullify(jobs%atom1,jobs%atom2,jobs%jobsize)
    call mem_alloc(jobs%atom1,njobs)
    call mem_alloc(jobs%atom2,njobs)
    call mem_alloc(jobs%jobsize,njobs)
    call mem_alloc(jobs%jobsdone,njobs)
    jobs%atom1=0
    jobs%atom2=0
    jobs%jobsize=0
    jobs%jobsdone=.false. ! no jobs are done

    ! MPI fragment statistics
    call mem_alloc(jobs%nslaves,njobs)
    call mem_alloc(jobs%nocc,njobs)
    call mem_alloc(jobs%nvirt,njobs)
    call mem_alloc(jobs%nbasis,njobs)
    call mem_alloc(jobs%ntasks,njobs)
    call mem_alloc(jobs%flops,njobs)
    call mem_alloc(jobs%LMtime,njobs)
    call mem_alloc(jobs%load,njobs)
    jobs%nslaves=0
    jobs%nocc=0
    jobs%nvirt=0
    jobs%nbasis=0
    jobs%ntasks=0
    jobs%flops=0.0E0_realk
    jobs%LMtime=0.0E0_realk
    jobs%load=0.0E0_realk

  end subroutine init_joblist



  !> \brief Free fragment job list
  !> \author Kasper Kristensen
  !> \date May 2012
  subroutine free_joblist(jobs)
    implicit none
    !> Job list
    type(joblist),intent(inout) ::  jobs

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

    if(associated(jobs%LMtime)) then
       call mem_dealloc(jobs%LMtime)
       nullify(jobs%LMtime)
    end if

    if(associated(jobs%load)) then
       call mem_dealloc(jobs%load)
       nullify(jobs%load)
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
    jobs%atom1(position) = singlejob%atom1(1)
    jobs%atom2(position) = singlejob%atom2(1)
    jobs%jobsize(position) = singlejob%jobsize(1)
    jobs%jobsdone(position) = singlejob%jobsdone(1)
    jobs%nslaves(position) = singlejob%nslaves(1)
    jobs%nocc(position) = singlejob%nocc(1)
    jobs%nvirt(position) = singlejob%nvirt(1)
    jobs%nbasis(position) = singlejob%nbasis(1)
    jobs%ntasks(position) = singlejob%ntasks(1)
    jobs%flops(position) = singlejob%flops(1)
    jobs%LMtime(position) = singlejob%LMtime(1)
    jobs%load(position) = singlejob%load(1)

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


  !> \brief Print short energy summary (both HF and correlation)
  !> (Necessary to place here because it is used both for DEC and for full calculation).
  !> \author Kasper Kristensen
  !> \date April 2013
  subroutine print_total_energy_summary(EHF,Ecorr,Eerr)
    implicit none
    !> HF energy
    real(realk),intent(in) :: EHF
    !> Correlation energy
    real(realk),intent(in) :: Ecorr
    !> Estimated intrinsic DEC energy error
    real(realk),intent(in) :: Eerr

    ! Print summary
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,*)
    write(DECinfo%output,'(13X,a)') '**********************************************************'
    write(DECinfo%output,'(13X,a,19X,a,19X,a)') '*', 'DEC ENERGY SUMMARY', '*'
    write(DECinfo%output,'(13X,a)') '**********************************************************'
    write(DECinfo%output,*)
    if(DECinfo%first_order) then
       write(DECinfo%output,'(15X,a,f20.10)') 'G: Hartree-Fock energy :', Ehf
       write(DECinfo%output,'(15X,a,f20.10)') 'G: Correlation energy  :', Ecorr
       ! skip error print for full calculation (0 by definition)
       if(.not. DECinfo%full_molecular_cc) then  
          write(DECinfo%output,'(15X,a,f20.10)') 'G: Estimated DEC error :', Eerr
       end if
       if(DECinfo%ccmodel==1) then
          write(DECinfo%output,'(15X,a,f20.10)') 'G: Total MP2 energy    :', Ehf+Ecorr
       elseif(DECinfo%ccmodel==2) then
          write(DECinfo%output,'(15X,a,f20.10)') 'G: Total CC2 energy    :', Ehf+Ecorr
       elseif(DECinfo%ccmodel==3) then
          write(DECinfo%output,'(15X,a,f20.10)') 'G: Total CCSD energy   :', Ehf+Ecorr
       elseif(DECinfo%ccmodel==4) then
          write(DECinfo%output,'(15X,a,f20.10)') 'G: Total CCSD(T) energy:', Ehf+Ecorr
       end if
    else
       write(DECinfo%output,'(15X,a,f20.10)') 'E: Hartree-Fock energy :', Ehf
       write(DECinfo%output,'(15X,a,f20.10)') 'E: Correlation energy  :', Ecorr
       ! skip error print for full calculation (0 by definition)
       if(.not. DECinfo%full_molecular_cc) then  
          write(DECinfo%output,'(15X,a,f20.10)') 'E: Estimated DEC error :', Eerr
       end if
       if(DECinfo%ccmodel==1) then
          write(DECinfo%output,'(15X,a,f20.10)') 'E: Total MP2 energy    :', Ehf+Ecorr
       elseif(DECinfo%ccmodel==2) then
          write(DECinfo%output,'(15X,a,f20.10)') 'E: Total CC2 energy    :', Ehf+Ecorr
       elseif(DECinfo%ccmodel==3) then
          write(DECinfo%output,'(15X,a,f20.10)') 'E: Total CCSD energy   :', Ehf+Ecorr
       elseif(DECinfo%ccmodel==4) then
          write(DECinfo%output,'(15X,a,f20.10)') 'E: Total CCSD(T) energy:', Ehf+Ecorr
       end if
    end if
    write(DECinfo%output,*)
    write(DECinfo%output,*)


  end subroutine print_total_energy_summary


end module dec_fragment_utils
