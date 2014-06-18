!this file contains the lower-order tools as opposed to the dec_utils which
!conains higher routines
module dec_tools_module
  use precision
  use dec_typedef_module
  use fundamental
  use memory_handling
  use files
  use BUILDAOBATCH


  !> Maximum number of files to be opened at the same time
  ! NOTE: If you change this number, you have to change
  ! max_file accordingly in crayio.c.
  integer, parameter :: max_number_files=250
  !> Number of opened files using C file handling
  integer,save :: files_opened=0
  !> Keeping control of which file units are available
  logical,save :: available_file_units(max_number_files)=.true.


  contains

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

  subroutine init_batch_info(mylsitem,batch,max_allowed,nb)
     implicit none
     type(lsitem), intent(inout) :: mylsitem
     type(int_batch),intent(inout) :: batch
     integer, intent(in)  :: max_allowed,nb
     integer :: i,iorb,k
     logical :: master
     master = .true.
#ifdef VAR_MPI
     master = (infpar%lg_mynum == infpar%master)
#endif

     ! Orbital to batch information
     ! ----------------------------
     call mem_alloc(batch%orb2batch,nb)
     call build_batchesofAOS(DECinfo%output,mylsitem%setting,max_allowed,&
        & nb,batch%max_dim,batch%batchsize,batch%batchdim,batch%batchindex,&
        &batch%nbatches,batch%orb2batch,'R')

     ! Translate batchindex to orbital index
     ! -------------------------------------
     call mem_alloc(batch%batch2orb,batch%nbatches)
     do i=1,batch%nbatches
        call mem_alloc(batch%batch2orb(i)%orbindex,batch%batchdim(i))
        batch%batch2orb(i)%orbindex = 0
        batch%batch2orb(i)%norbindex = 0
     end do
     do iorb=1,nb
        i = batch%orb2batch(iorb)
        batch%batch2orb(i)%norbindex = batch%batch2orb(i)%norbindex+1
        K = batch%batch2orb(i)%norbindex
        batch%batch2orb(i)%orbindex(K) = iorb
     end do
  end subroutine init_batch_info

  subroutine free_batch_info(batch)
     type(int_batch),intent(inout) :: batch
     integer :: i
     call mem_dealloc(batch%orb2batch)
     call mem_dealloc(batch%batchdim)
     call mem_dealloc(batch%batchsize)
     call mem_dealloc(batch%batchindex)
     do i=1,batch%nbatches
        call mem_dealloc(batch%batch2orb(i)%orbindex)
        batch%batch2orb(i)%orbindex => null()
     end do
     call mem_dealloc(batch%batch2orb)
  end subroutine free_batch_info
end module dec_tools_module
