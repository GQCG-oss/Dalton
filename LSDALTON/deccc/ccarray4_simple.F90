!> @file
!> Simple structure to store four-dimensional arrays
!> \author Marcin Ziolkowski

!> Four-dimensional array operations
module array4_simple_operations


  ! Outside DEC directory
  use memory_handling
  use precision
  use dec_typedef_module
  use files!,only: lsopen,lsclose
  use LSTIMING!,only:lstimer
  use reorder_frontend_module

  ! DEC DEPENDENCIES (within deccc directory)   
  ! *****************************************
  use array4_memory_manager
  use dec_tools_module
  use crayio_tools_module

  !> Number of array
  integer(kind=long) :: ArrayNumber=0
  !> Number of created arrays
  integer(kind=long) :: CreatedArrays=0
  !> Number of destroyed arrays
  integer(kind=long) :: DestroyedArrays=0


  !> Overloaded operator for adding arrays
  interface operator(+)
     module procedure array4_add_simple
  end interface

  !> Overloaded operator for dot-product
  interface operator(*)
     module procedure array4_dotproduct
  end interface

  !> Write array4 elements to file
  interface array4_write_file
     module procedure array4_write_file_type1
     module procedure array4_write_file_type2
     module procedure array4_write_file_type3
  end interface

  !> Read array4 elements from file
  interface array4_read_file
     module procedure array4_read_file_type1
     module procedure array4_read_file_type2
     module procedure array4_read_file_type3
  end interface


  !> Init arrays
  interface array4_init
     module procedure array4_init_standard
     module procedure array4_init_file
  end interface

  !> transposition routines
  !interface transp
  !   module procedure partial_transpose_ip
  !   module procedure partial_transpose
  !   module procedure mat_transpose
  !end interface

  ! for timing
  real(realk) :: time_array4_init = 0.0E0_realk
  real(realk) :: time_array4_free = 0.0E0_realk
  real(realk) :: time_array4_duplicate = 0.0E0_realk
  real(realk) :: time_array4_alloc = 0.0E0_realk
  real(realk) :: time_array4_dealloc = 0.0E0_realk
  real(realk) :: time_array4_write_data = 0.0E0_realk
  real(realk) :: time_array4_read_data = 0.0E0_realk
  real(realk) :: time_array4_write = 0.0E0_realk
  real(realk) :: time_array4_read = 0.0E0_realk
  real(realk) :: time_array4_add = 0.0E0_realk
  real(realk) :: time_array4_add_to = 0.0E0_realk
  real(realk) :: time_array4_scale = 0.0E0_realk
  real(realk) :: time_array4_dotproduct = 0.0E0_realk
  real(realk) :: time_array4_norm = 0.0E0_realk
  real(realk) :: time_array4_contract1 = 0.0E0_realk
  real(realk) :: time_array4_contract2 = 0.0E0_realk
  real(realk) :: time_array4_contract2_middle = 0.0E0_realk
  real(realk) :: time_array4_contract_array2 = 0.0E0_realk
  real(realk) :: time_array4_contract3 = 0.0E0_realk
  real(realk) :: time_array4_mat_multiply = 0.0E0_realk
  real(realk) :: time_array4_reorder = 0.0E0_realk
  real(realk) :: time_array4_mat_transpose = 0.0E0_realk

  real(realk) :: time_array4_contract_with_ao = 0.0E0_realk



contains






  !> \brief Initialize four-dimensional array
  !> If we save array4 on file, we do NOT allocate memory here.
  !> \author Marcin Ziolkowski
  !> \param dims Dimensions
  !> \return Four-dimensional array
  function array4_init_standard(dims,distribution) result(array)
    implicit none
    type(array4) :: array
    integer, dimension(2), intent(in), optional :: distribution
    integer, dimension(4), intent(in) :: dims
    logical :: file_exist
    integer :: ierr
    real(realk) :: t0,t1
    
    call cpu_time(t0)
    CreatedArrays = CreatedArrays+1
    array%dims=dims
    array%order=[1,2,3,4]
 
    call memory_allocate_4d(array%val,dims)
    call ls_dzero(array%val,size(array%val))
 
    ! By default: Do not create files when arrays are stored in memory.
    array%funit=0
    array%filename = 'NoFilename'
 
 
    ! Nullify information needed for storing on file
    ! Note: If you want to store the array4 on file, then 
    ! array4_init_file should be used intead!
 
    ! Address counter set to zero
    array%address_counter=0
 
    ! storing_type=0 --> No storing on file
    array%storing_type=0
 
    ! nelements in each chunk on file = 0 --> No storing on file
    array%nelements=0
 
    ! Nullify address pointer for file
    nullify(array%address)
    call cpu_time(t1)
    time_array4_init = time_array4_init + (t1-t0)

    return
  end function array4_init_standard


  !> \brief Initialize four-dimensional array and set information
  !> used when storing array on file.
  !> Should ONLY be used when DECinfo%array4OnFile=.true.
  !> Note: No memory for the array values is allocated here!
  !> \author Kasper Kristensen
  !> \date October 2010
  function array4_init_file(dims,storing_type,zero_elements) result(array)

    implicit none
    !> Array to initialize
    type(array4) :: array
    !> Dimensions for array
    integer, dimension(4), intent(in) :: dims
    !> Storing type when writing/reading to/from file (see type array4 for details)
    integer, intent(in) :: storing_type
    !> Set all array elements in file to zero
    !> (If this is false, the array file is created but it is empty)
    logical, intent(in) :: zero_elements
    logical :: file_exist
    integer :: ierr
    real(realk) :: t0,t1

    call cpu_time(t0)

    ! Sanity check
    if(.not. DECinfo%array4OnFile) then
       write(DECinfo%output,*) 'array4_init_file should not be called &
            & when array4OnFile=.false.'
       call lsquit('array4_init_file should not be called when &
            &array4OnFile=.false.', DECinfo%output)
    end if

    array%dims=dims
    array%order=[1,2,3,4]

    ArrayNumber = ArrayNumber+1
    CreatedArrays = CreatedArrays+1

    ! file
    write(array%FileName,'("tensor_",i6.6,".data")') ArrayNumber


    ! Get file unit for array
    call get_available_file_unit(array%funit)

    ! in case the file already exists
    inquire(file=array%FileName,exist=file_exist)
    if(file_exist) then
       open(array%FUnit,file=array%FileName,form='unformatted',&
            access='sequential',status='old')
       close(array%FUnit,status='delete')
    end if

    ! Nullify values in array (just in case)
    nullify(array%val)


    ! Information needed for storing on file
    ! **************************************

    ! Address counter set to 1
    array%address_counter=1

    ! Storing type (see type array4 for details)
    array%storing_type=storing_type

    ! Number of elements in each chunk on file
    ! This depends on the storing type (see type array4 for details)

    StoringType: select case(storing_type)

    case(1)
       ! For storing_type=1 we write to file chunks of the form array%val(:,:,:,n)
       ! I.e. in each chunk there are dims(1)*dims(2)*dims(3) elements.
       array%nelements = array%dims(1)*array%dims(2)*array%dims(3)

       ! The addresses identifiers are stored in array%address(:,:,:,:)
       ! array%address is by definition a four dimensional array but here we only use the
       ! first index to define a given address.
       ! We therefore set the following dimensions:
       ! address(1), address(2), address(3) : 1
       ! address(4) : dims(4)
       call mem_alloc(array%address,1,1,1,dims(4) )

    case(2)
       ! For storing_type=2 we write to file chunks of the form array%val(:,:,n1,n2)
       ! I.e. in each chunk there are array%dims(1)*array%dims(2) elements.
       array%nelements = array%dims(1)*array%dims(2)

       ! The addresses identifiers are stored in array%address(:,:,:,:)
       ! array%address is by definition a four dimensional array but we only use the
       ! last two indices to define a given address.
       ! We therefore set the following dimensions:
       ! address(1) and address(2) : 1    (this is necessary to keep the structure general)
       ! address(3) and address(4) : dims(3) and dims(4)
       call mem_alloc(array%address,1,1,dims(3),dims(4) )

    case(3)
       ! For storing_type=3 we write to file chunks of the form array%val(:,n1,n2,n3)
       ! I.e. in each chunk there are array%dims(1) elements.
       array%nelements = array%dims(1)

       ! The addresses identifiers are stored in array%address(:,:,:,:)
       ! array%address is by definition a four dimensional array but we only use the
       ! last three indices to define a given address.
       ! We therefore set the following dimensions:
       ! address(1) : 1    (this is necessary to keep the structure general)
       ! address(2), address(3), and address(4) : dims(2), dims(3) and dims(4)
       call mem_alloc(array%address,1,dims(2),dims(3),dims(4) )


    case default
       write(DECinfo%output,*) 'array4_init_file: Requested storing type is not &
            & implemented:', array%storing_type
       call lsquit('array4_init_file: Requested storing type is not &
            & implemented.', DECinfo%output)

    end select StoringType


    ! To keep the general framework we may want to write zero elements to the array file
    ! Admittedly, this is not very pretty - but it's working.
    ! Sometime this should be replaced by a more elegant solution.
    ! In any case we set the addresses (array%address) are set to default addresses,
    ! which are the "natural choices" of addresses for Fortran -
    ! e.g. for storing type 2 we store in the order:
    ! (1,1) (2,1) (3,1) .... (n1,1)  -->
    ! (1,2) (2,2) (3,2) .... (n1,2)  -->
    ! ....                          -->
    ! (1,n2) (2,n2) (3,n2) .... (n1,n2)
    ! To catch up possible errors we first set the addresses to 0.
    array%address=0
    if(zero_elements) then
       call array4_zero_file(array)
    else
       call array4_set_standard_address(array)
    end if

    call cpu_time(t1)
    time_array4_init = time_array4_init + (t1-t0)

    return
  end function array4_init_file


  !> \bried Destroy four-dimensional array
  !> \author Marcin Ziolkowski
  !> \param array Four-dimensional array to be deleted
  !> \param keep Keep file
  subroutine array4_free(array,keep)

    implicit none
    type(array4), intent(inout) :: array
    logical, optional :: keep
    real(realk) :: t0,t1
    logical :: file_exists

    call cpu_time(t0)

    DestroyedArrays = DestroyedArrays + 1

    call memory_deallocate_4d(array%val)

    inquire(file=array%FileName,exist=file_exists)


    FileExists: if(file_exists) then


       if(DECinfo%array4OnFile) then ! always delete file using C routine

          if(present(keep)) then
             if(.not. keep) then
                call array4_delete_file(array)
             end if
          else
             call array4_delete_file(array)
          end if

          ! Also deallocate address elements
          if(associated(array%address)) then
             call mem_dealloc(array%address)
             nullify(array%address)
          endif

       else
          call lsopen(array%FUnit,array%FileName,'OLD','UNFORMATTED')
          if(present(keep)) then
             if(keep) then
                call lsclose(array%FUnit,'KEEP')
             else
                call lsclose(array%FUnit,'DELETE')
             end if
          else
             call lsclose(array%FUnit,'DELETE')
          end if
       end if

       ! File associated with array is no longer used so that file unit is again available.
       !available_file_units(array%funit)=.true.
       array%funit=0

    end if FileExists

    call cpu_time(t1)
    time_array4_free  = time_array4_free + (t1-t0)

    return
  end subroutine array4_free




  !> \brief Create a copy of an array, but let
  !> the old and the new arrays point to the SAME FILE.
  !> WARNING!!! Do not use this routine unless you are completely sure
  !> that the arrays "this" and "res" are never used at the same time!!!
  !> \author Kasper Kristensen
  !> \param this Four-dimensional array to be copied
  !> \return New array that is a duplicate of the requested array
  function array4_duplicate_same_file(this) result(res)
    implicit none
    !> Array to be copied
    type(array4), intent(in) :: this
    !> New array, identical to this, and referencing the same file.
    type(array4) :: res
    integer(kind=long) :: vec_size64
    integer :: vec_size
    logical :: file_exists


    ! Only use when array values are stored on file,
    ! otherwise it is meaningless.
    ArrayOnFile: if(.not. DECinfo%array4OnFile) then
       call lsquit('array4_duplicate_same_file: &
            & This routine should only be called when the .array4OnFile is used!',-1)
    end if ArrayOnFile
    res = array4_init_file(this%dims, this%storing_type, .false.)

    ! Delete newly contructed file for array res - because we want
    ! associate array res with the same file as array this.
    call array4_delete_file(res)

    ! Set file addresses
    res%address = this%address

    ! Associate res with the SAME FILE as this
    inquire(file=this%Filename,exist=file_exists)
    if(.not. file_exists) then
       call lsquit('array4_duplicate: File for array does not exist!', DECinfo%output)
    end if
    res%FileName = this%FileName

    ! Same ordering for the two arrays
    res%order = this%order

    ! copy data
    DataInMemory: if(associated(this%val)) then

       vec_size64 = int(this%dims(1)*this%dims(2)* &
            this%dims(3)*this%dims(4),kind=8)
       if(vec_size64>MAXINT)then
          call lsquit('ERROR(array4_duplicate_same_file): size of array cannot be &
                       &described by current integer type, please try another &
                       &compilation or fix this routine', DECinfo%output)
       endif
       vec_size = this%dims(1)*this%dims(2)* &
            & this%dims(3)*this%dims(4)
       call dcopy(vec_size,this%val,1,res%val,1)

    end if DataInMemory

  end function array4_duplicate_same_file




  !> \brief Create a copy of an array
  !> \author Marcin Ziolkowski (modified by Kasper Kristensen)
  !> \param this Four-dimensional array to be copied
  !> \return New array that is a duplicate of the requested array
  function array4_duplicate(this) result(res)
    implicit none
    type(array4), intent(in) :: this
    type(array4) :: res
    real(realk) :: t0,t1
    integer :: vec_size
    integer(kind=long) :: vec_size64
    logical :: file_exists

    call cpu_time(t0)

    ! create and allocate
    ArrayOnFile: if(DECinfo%array4OnFile) then
       res = array4_init_file(this%dims, this%storing_type, .false.)

       ! Set file addresses
       res%address = this%address

       ! Copy file
       inquire(file=this%Filename,exist=file_exists)
       if(.not. file_exists) then
          call lsquit('array4_duplicate: File for array does not exist!', DECinfo%output)
       end if
       call copyfile(this%FileName, res%FileName)

    else
       res = array4_init(this%dims)
    end if ArrayOnFile

    ! Same ordering for the two arrays
    res%order = this%order

    ! copy data
    DataInMemory: if(associated(this%val)) then

       vec_size64 = int(this%dims(1)*this%dims(2)*this%dims(3)*this%dims(4),kind=8)
       if(vec_size64>MAXINT)then
          call lsquit('ERROR(array4_duplicate): size of array cannot be &
                       &described by current integer type, please try another &
                       &compilation or fix this routine', DECinfo%output)
       endif
       vec_size = this%dims(1)*this%dims(2)*this%dims(3)*this%dims(4)
       call dcopy(vec_size,this%val,1,res%val,1)

    end if DataInMemory


    call cpu_time(t1)
    time_array4_duplicate = time_array4_duplicate + (t1-t0)

  end function array4_duplicate

  !> \brief Allocate memory for data assiciated with the array
  !> \author Marcin Ziolkowski
  !> \param array Array with deallocated memory we will allocate
  subroutine array4_alloc(array)
    implicit none
    type(array4), intent(inout) :: array
    real(realk) :: t0,t1
    integer(kind=long) :: vec_size64
    integer :: vec_size

    call cpu_time(t0)
    call memory_allocate_4d(array%val,array%dims)
    vec_size64 = int(array%dims(1)*array%dims(2)* &
         array%dims(3)*array%dims(4),kind=8)
    if(vec_size64>MAXINT)then
       call lsquit('ERROR(array4_alloc): size of array cannot be &
                    &described by current integer type, please try another &
                    &compilation or fix this routine', DECinfo%output)
    endif
    vec_size = array%dims(1)*array%dims(2)* &
         array%dims(3)*array%dims(4)
    call ls_dzero(array%val,vec_size)
    call cpu_time(t1)
    time_array4_alloc = time_array4_alloc + (t1-t0)
    return
  end subroutine array4_alloc

  !> \breif Deallocate memory (release data associated with this array)
  subroutine array4_dealloc(array)
    implicit none
    type(array4), intent(inout) :: array
    real(realk) :: t0,t1

    call cpu_time(t0)
    call memory_deallocate_4d(array%val)
    call cpu_time(t1)
    time_array4_dealloc = time_array4_dealloc + (t1-t0)
    return
  end subroutine array4_dealloc

  !> \brief Write data to the file assiciated with this array
  subroutine array4_write_data(array)

    implicit none
    type(array4), intent(inout) :: array
    real(realk) :: t0,t1
    logical :: is_open, file_exist

    call cpu_time(t0)
    inquire(file=array%FileName,exist=file_exist)
    inquire(file=array%FileName,opened=is_open)
    if(.not.is_open) then
       if(file_exist) then
          call lsopen(array%FUnit,array%FileName,'OLD','UNFORMATTED')
       else
          call lsopen(array%FUnit,array%FileName,'REPLACE','UNFORMATTED')
       end if
    end if

    rewind(array%FUnit)
    write(array%FUnit) array%val
    call lsclose(array%FUnit,'KEEP')
    call cpu_time(t1)
    time_array4_write_data = time_array4_write_data + (t1-t0)
    return
  end subroutine array4_write_data

  !> \brief Read data from file
  subroutine array4_read_data(array)

    implicit none
    type(array4), intent(inout) :: array
    real(realk) :: t0,t1
    logical :: is_open

    call cpu_time(t0)
    inquire(file=array%FileName,opened=is_open)
    if(.not.is_open) then
       call lsopen(array%FUnit,array%FileName,'OLD','UNFORMATTED')
    end if
    rewind(array%FUnit)
    read(array%FUnit) array%val
    call lsclose(array%FUnit,'KEEP')
    call cpu_time(t1)
    time_array4_read_data = time_array4_read_data + (t1-t0)
    return
  end subroutine array4_read_data

  !> \brief Write data and deallocate memory (does not destroy an array)
  subroutine array4_write(array)
    implicit none
    type(array4), intent(inout) :: array
    real(realk) :: t0,t1
    call cpu_time(t0)
    call array4_write_data(array)
    call array4_dealloc(array)
    call cpu_time(t1)
    time_array4_write = time_array4_write + (t1-t0)
    return
  end subroutine array4_write

  !> \brief Allocate memory and read data (does not create an array)
  subroutine array4_read(array)
    implicit none
    type(array4), intent(inout) :: array
    real(realk) :: t0,t1
    call cpu_time(t0)
    call array4_alloc(array)
    call array4_read_data(array)
    call cpu_time(t1)
    time_array4_read = time_array4_read + (t1-t0)
    return
  end subroutine array4_read

  !> \brief Simple addition
  function array4_add_simple(A,B) result(C)
    implicit none
    type(array4), intent(in) :: A,B
    type(array4) :: C
    C = array4_add(1.0E0_realk,A,1.0E0_realk,B)
    return
  end function array4_add_simple


  !> \brief Add scaled arrays: array = alpha*A + beta*B
  function array4_add(alpha,A,beta,B) result(array)
    implicit none
    type(array4) :: array
    type(array4), intent(in) :: A,B
    real(realk), intent(in) :: alpha,beta
    real(realk) :: t0,t1
    integer :: i

    call cpu_time(t0)


    ! Sanity checks: Dimensions and orders match
    do i=1,4
       if(A%dims(i) /= B%dims(i) .or. A%order(i) /= B%order(i)) then
          call lsquit('array4_add: Dimensions or orders of &
               & A and B arrays do not match',DECinfo%output)
       end if
    end do

    if(DECinfo%array4OnFile) then ! array elements stored on file
       array = array4_add_file(alpha,A,beta,B)
    else ! array elements kept in memory
       array = array4_add_memory(alpha,A,beta,B)
    end if

    call cpu_time(t1)
    time_array4_add = time_array4_add + (t1-t0)

  end function array4_add


  !> \brief Add scaled arrays when array elements are kept in memory
  function array4_add_memory(alpha,A,beta,B) result(array)
    implicit none
    type(array4) :: array
    type(array4), intent(in) :: A,B
    real(realk), intent(in) :: alpha,beta
    integer(kind=long) :: nelements64
    integer :: nelements

    array=array4_init(A%dims)
    !#ifdef USE_BLAS
    nelements64=array%dims(1)*array%dims(2)*array%dims(3)*array%dims(4)
    if(nelements64>MAXINT)then
       call lsquit('ERROR(array4_add_memory): size of array cannot be &
                    &described by current integer type, please try another &
                    &compilation or fix this routine', DECinfo%output)
    endif
    nelements=array%dims(1)*array%dims(2)*array%dims(3)*array%dims(4)
    call daxpy(nelements,alpha,A%val,1,array%val,1)
    call daxpy(nelements,beta,B%val,1,array%val,1)
    !#else
    !    array%val=alpha*A%val+beta*B%val
    !#endif

  end function array4_add_memory


  !> \brief Add scaled arrays when array elements are stored on file
  function array4_add_file(alpha,A,beta,B) result(array)
    implicit none
    type(array4) :: array
    type(array4), intent(in) :: A,B
    real(realk), intent(in) :: alpha,beta


    ! Sanity check
    ! ************

    ! 1. Array A and B use the same storing type
    if(A%storing_type /= B%storing_type) then
       call lsquit('array4_add_to_file: Only implemented when &
            & array A and B use the same storing type!',DECinfo%output)
    end if

    ! Choose storing type
    StoringType: select case(A%storing_type)

    case(2)
       array = array4_add_file_type2(alpha,A,beta,B)

    case default
       call lsquit('array4_add_file: Requested &
            & storing type is not implemented!',DECinfo%output)

    end select StoringType

  end function array4_add_file



  !> \brief Add scaled arrays when elements are stored on file
  !>  using storing type 2.
  !> \author Kasper Kristensen
  !> \date October 2010
  function array4_add_file_type2(alpha,A,beta,B) result(array)
    implicit none
    type(array4) :: array
    type(array4), intent(in) :: A,B
    real(realk), intent(in) :: alpha,beta
    integer :: k,l
    integer :: nelements,dim1,dim2,dim3,dim4
    real(realk),pointer :: Aval(:,:), Bval(:,:), arrayval(:,:)


    ! Initialize stuff
    ! ****************
    dim1 = A%dims(1)
    dim2 = A%dims(2)
    dim3 = A%dims(3)
    dim4 = A%dims(4)
    call mem_alloc(Aval,dim1,dim2)
    call mem_alloc(Bval,dim1,dim2)
    call mem_alloc(arrayval,dim1,dim2)
    nelements=dim1*dim2
    ! Initialize array using storing type 2
    array = array4_init(A%dims,2,.false.)

    ! Open files
    ! ----------
    call array4_open_file(A)
    call array4_open_file(B)
    call array4_open_file(array)



    ! Carry out summation: array = alpha*A + beta*B
    ! *********************************************


    l_loop: do l=1,dim4
       k_loop: do k=1,dim3

          arrayval = 0E0_realk

          ! Read in A(:,:,k,l) and B(:,:,k,l)
          call array4_read_file(A,k,l,Aval,dim1,dim2)
          call array4_read_file(B,k,l,Bval,dim1,dim2)

          ! arrayval = alpha*A
          call daxpy(nelements,alpha,Aval,1,arrayval,1)
          ! arrayval = arrayval + beta*B = alpha*A + beta*B
          call daxpy(nelements,beta,Bval,1,arrayval,1)

          ! Now arrayval contains elements (:,:,k,l) of the
          ! the 4-dimensional array=alpha*A+beta*B.
          ! Write these elements to the file referenced by array
          call array4_write_file(array,k,l,arrayval,dim1,dim2)


       end do k_loop
    end do l_loop


    ! Free stuff
    call mem_dealloc(Aval)
    call mem_dealloc(Bval)
    call mem_dealloc(arrayval)
    call array4_close_file(A,'keep')
    call array4_close_file(B,'keep')
    call array4_close_file(array,'keep')


  end function array4_add_file_type2




  !> \brief Scale B and add to A
  subroutine array4_add_to(A,beta,B)
    implicit none
    type(array4), intent(inout) :: A
    type(array4), intent(in) :: B
    real(realk), intent(in) :: beta
    real(realk) :: t0,t1

    call cpu_time(t0)

    if(DECinfo%array4OnFile) then
       call array4_add_to_file(A,beta,B)
    else
       call array4_add_to_memory(A,beta,B)
    end if

    call cpu_time(t1)
    time_array4_add_to = time_array4_add_to + (t1-t0)
    return
  end subroutine array4_add_to


  !> \brief Scale B and add to A when array elements are stored in memory
  subroutine array4_add_to_memory(A,beta,B)
    implicit none
    type(array4), intent(inout) :: A
    type(array4), intent(in) :: B
    real(realk), intent(in) :: beta
    integer(kind=long) :: nelements64
    integer :: nelements

    nelements64=A%dims(1)*A%dims(2)*A%dims(3)*A%dims(4)
    if(nelements64>MAXINT)then
       call lsquit('ERROR(array4_alloc): size of array cannot be &
                    &described by current integer type, please try another &
                    &compilation or fix this routine', DECinfo%output)
    endif
    nelements=A%dims(1)*A%dims(2)*A%dims(3)*A%dims(4)
    call daxpy(nelements,beta,B%val,1,A%val,1)

  end subroutine array4_add_to_memory



  !> \brief Scale B and add to A when array elements are stored only on file
  subroutine array4_add_to_file(A,beta,B)
    implicit none
    type(array4), intent(inout) :: A
    type(array4), intent(in) :: B
    real(realk), intent(in) :: beta

    ! Sanity checks
    ! *************

    ! 1. Array A and B use the same storing type
    if(A%storing_type /= B%storing_type) then
       call lsquit('array4_add_to_file: Only implemented when &
            & array A and B use the same storing type!',DECinfo%output)
    end if

    ! 2. Array A and B have the same dimensions
    if( A%dims(1) /= B%dims(1) .or. A%dims(2) /= B%dims(2) .or. &
         & A%dims(3) /= B%dims(3) .or. A%dims(4) /= B%dims(4) ) then
       call lsquit('array4_add_to_file: &
            & Array A and B have different dimensions!',DECinfo%output)
    end if

    ! Choose storing type
    StoringType: select case(A%storing_type)

    case(2)
       call array4_add_to_file_type2(A,beta,B)

    case default
       call lsquit('array4_add_to_file: Requested &
            & storing type is not implemented!',DECinfo%output)

    end select StoringType

  end subroutine array4_add_to_file


  !> \brief Scale B and add to A when array elements are stored on file
  !> using storing type 2.
  subroutine array4_add_to_file_type2(A,beta,B)
    implicit none
    type(array4), intent(inout) :: A
    type(array4), intent(in) :: B
    real(realk), intent(in) :: beta
    integer :: i,j,k,l,nelements
    integer :: dim1,dim2,dim3,dim4
    real(realk),pointer :: Aval(:,:), Bval(:,:)
    real(realk) :: tcpu, twall


    ! Initialize stuff
    ! ****************
    dim1 = A%dims(1)
    dim2 = A%dims(2)
    dim3 = A%dims(3)
    dim4 = A%dims(4)
    call mem_alloc(Aval,dim1,dim2)
    call mem_alloc(Bval,dim1,dim2)
    nelements=dim1*dim2


    ! Open files
    ! **********
    call array4_open_file(A)
    call array4_open_file(B)



    ! Carry out summation: A --> A + beta*B
    ! *************************************

    l_loop: do l=1,dim4
       k_loop: do k=1,dim3

          ! Read in A(:,:,k,l) and B(:,:,k,l)
          call array4_read_file(A,k,l,Aval,dim1,dim2)
          call array4_read_file(B,k,l,Bval,dim1,dim2)

          ! Aval = Aval + beta*Bval
          call daxpy(nelements,beta,Bval,1,Aval,1)

          ! Now Aval contains elements (:,:,k,l) of the new A, i.e.
          ! Aval = Aold(:,:,k,l) + beta*B(:,:,k,l)
          ! Write these elements to the file referenced by A.
          ! NOTE!! It is necessary to set the address counter here, otherwise
          ! it will be changed by arary4_write_file.
          A%address_counter=A%address(1,1,k,l)
          call array4_write_file(A,k,l,Aval,dim1,dim2)

       end do k_loop
    end do l_loop


    ! Free stuff
    call mem_dealloc(Aval)
    call mem_dealloc(Bval)
    call array4_close_file(A,'keep')
    call array4_close_file(B,'keep')

  end subroutine array4_add_to_file_type2




  !> \brief Scale array values by factor
  subroutine array4_scale(this,factor)
    implicit none
    type(array4), intent(inout) :: this
    real(realk), intent(in) :: factor
    real(realk) :: t0,t1
    call cpu_time(t0)

    if(DECinfo%array4OnFile) then ! arrays are stored on file
       call array4_scale_file(this,factor)

    else ! arrays are kept in memory, just scale them by factor
       this%val = factor*this%val
    end if

    call cpu_time(t1)
    time_array4_scale = time_array4_scale + (t1-t0)

  end subroutine array4_scale

  !> \brief Scale array values by factor when array values are stored on file.
  !> \author Kasper Kristensen
  !> \date October 2010
  subroutine array4_scale_file(this,factor)
    implicit none
    type(array4), intent(inout) :: this
    real(realk), intent(in) :: factor

    StoringType: select case(this%storing_type)
    case(2)
       call array4_scale_file_type2(this,factor)
    case default
       write(DECinfo%output,*) 'array4_scale_file: Requested storing type is not &
            & implemented:', this%storing_type
       call lsquit('array4_dotproduct_file: Requested storing type is not &
            & implemented.', DECinfo%output)
    end select StoringType

  end subroutine array4_scale_file


  !> \brief Scale array values by factor when array values are stored on file
  !> using storing type 2.
  !> \author Kasper Kristensen
  !> \date October 2010
  subroutine array4_scale_file_type2(this,factor)
    implicit none
    !> Array to scale
    type(array4), intent(inout) :: this
    !> Factor to scale array by
    real(realk), intent(in) :: factor
    real(realk),pointer :: values(:,:)
    integer :: idx3,idx4,dim1,dim2


    ! Initialize stuff
    ! ****************

    ! Temporary 2-dimensional array for storing values
    dim1 = this%dims(1)
    dim2 = this%dims(2)
    call mem_alloc(values,dim1,dim2)

    ! Open file for original array
    call array4_open_file(this)

    ! Set address counter for new scaled array to 1
    this%address_counter=1

    ! Loop over index 3 and 4 and scale values
    idx4_loop: do idx4=1,this%dims(4)
       idx3_loop: do idx3=1,this%dims(3)

          ! Read original file
          call array4_read_file(this,idx3,idx4,values,dim1,dim2)

          ! Scale values
          values=factor*values

          ! Write values to file (overwrite old non-scaled values)
          this%address_counter=this%address(1,1,idx3,idx4)
          call array4_write_file(this,idx3,idx4,values,dim1,dim2)

       end do idx3_loop
    end do idx4_loop

    ! Close and keep file for new scaled array
    call array4_close_file(this,'keep')

    ! Deallocate memory
    call mem_dealloc(values)


  end subroutine array4_scale_file_type2


  !> \brief Dot-product
  function array4_dotproduct(A,B) result(res)
    implicit none
    type(array4), intent(in) :: A,B
    real(realk) :: res
    integer :: i,j,k
    real(realk) :: t0,t1

    res=0E0_realk
    call cpu_time(t0)

    if(DECinfo%array4OnFile) then ! arrays stored on file
       res=array4_dotproduct_file(A,B)
    else ! arrays stored in memory
       res=array4_dotproduct_memory(A,B)
    end if

    call cpu_time(t1)
    time_array4_dotproduct = time_array4_dotproduct + (t1-t0)
    return
  end function array4_dotproduct


  !> \brief Dot-product
  function array4_dotproduct_memory(A,B) result(res)
    implicit none
    type(array4), intent(in) :: A,B
    real(realk) :: res
    integer :: i,j,k
    integer(kind=long) :: nelements64
    integer :: nelements
    real(realk), external :: ddot


    res=0.0E0_realk
    !#ifdef USE_BLAS
    nelements64=A%dims(1)*A%dims(2)*A%dims(3)*A%dims(4)
    if(nelements64>MAXINT)then
       call lsquit('ERROR(array4_dotproduct_memory): size of array cannot be &
                    &described by current integer type, please try another &
                    &compilation or fix this routine', DECinfo%output)
    endif
    nelements=A%dims(1)*A%dims(2)*A%dims(3)*A%dims(4)
    res=ddot(nelements,A%val,1,B%val,1)
!!$#else
!!$    do k=1,A%dims(4)
!!$       do j=1,A%dims(3)
!!$          do i=1,A%dims(2)
!!$             res=res+dot_product(A%val(:,i,j,k),B%val(:,i,j,k))
!!$          end do
!!$       end do
!!$    end do
!!$#endif


    return
  end function array4_dotproduct_memory



  !> \brief Get dot product of arrays A and B when elements are stored on file
  !> \author Kasper Kristensen
  !> \date October 2010
  function array4_dotproduct_file(A,B) result(res)
    implicit none
    !> Array for which norm is requested
    type(array4), intent(in) :: A,B
    !> Dot product
    real(realk) :: res

    ! Sanity checks
    ! *************

    ! 1. Currently only implemented when A and B use the same storing type
    if(A%storing_type /= B%storing_type) then
       call lsquit('array4_dotproduct_file: Dot products only implemented &
            & when array A and B use the same storing type!', DECinfo%output)
    end if

    ! 2. A and B must have the same dimensions
    if( A%dims(1) /= B%dims(1) .or. A%dims(2) /= B%dims(2) &
         & .or. A%dims(3) /= B%dims(3) .or. A%dims(4) /= B%dims(4) ) then
       call lsquit('array4_dotproduct_file: &
            & Dimensions of A and B do not match!', DECinfo%output)
    end if

    ! If A and B refer to the same file (i.e. they are the same array)
    ! then we should call array4_norm instead of array4_dotproduct
    AeqB: if(A%filename == B%filename) then

       call array4_norm_file(A,res)

    else ! A/=B and we call dotproduct routines

       StoringType: select case(A%storing_type)

       case(2)
          call array4_dotproduct_file_type2(A,B,res)

       case(3)
          call array4_dotproduct_file_type3(A,B,res)

       case default
          write(DECinfo%output,*) 'array4_dotproduct_file: Requested storing type is not &
               & implemented:', A%storing_type
          call lsquit('array4_dotproduct_file: Requested storing type is not &
               & implemented.', DECinfo%output)

       end select StoringType

    end if AeqB


  end function array4_dotproduct_file


  !> \brief Get dot product A*B when the elements are stored on file
  !> using storing type 2 (see type array4_init_file).
  !> \author Kasper Kristensen
  !> \date October 2010
  subroutine array4_dotproduct_file_type2(A,B,res)
    implicit none
    !> Array for which norm is requested
    type(array4), intent(in) :: A,B
    !> Dotproduct
    real(realk), intent(inout) :: res
    real(realk),pointer :: Avalues(:,:), Bvalues(:,:)
    integer :: idx3,idx4,dim1,dim2,nelements,i
    real(realk), external :: ddot

    ! Initialize stuff
    res = 0E0_realk
    dim1 = A%dims(1)
    dim2 = A%dims(2)
    nelements=dim1*dim2
    call mem_alloc(Avalues,dim1,dim2)
    call mem_alloc(Bvalues,dim1,dim2)

    ! Open file for arrays
    call array4_open_file(A)
    call array4_open_file(B)

    ! Loop over index 3 and 4 and take the dot products
    idx4_loop: do idx4=1,A%dims(4)
       idx3_loop: do idx3=1,A%dims(3)

          ! Read in A and B
          call array4_read_file_type2(A,idx3,idx4,Avalues,dim1,dim2)
          call array4_read_file_type2(B,idx3,idx4,Bvalues,dim1,dim2)

          !#ifdef USE_BLAS
          res = res + ddot(nelements,Avalues,1,Bvalues,1)
          !#else
          !          do i=1,dim2
          !             res=res+dot_product(Avalues(:,i),Bvalues(:,i))
          !          end do
          !#endif

       end do idx3_loop
    end do idx4_loop



    ! Close files for arrays
    call array4_close_file(A,'keep')
    call array4_close_file(B,'keep')

    ! Deallocate memory
    call mem_dealloc(Avalues)
    call mem_dealloc(Bvalues)

  end subroutine array4_dotproduct_file_type2


  !> \brief Get dotproduct A*B when the elements are stored on file
  !> using storing type 3 (see type array4_init_file).
  !> \author Kasper Kristensen
  !> \date October 2010
  subroutine array4_dotproduct_file_type3(A,B,res)
    implicit none
    !> Array for which norm is requested
    type(array4), intent(in) :: A,B
    !> Dot product
    real(realk), intent(inout) :: res
    real(realk),pointer :: Avalues(:), Bvalues(:)
    integer :: idx2, idx3,idx4,nelements,i

    ! Initialize stuff
    res = 0E0_realk
    nelements = A%dims(1)
    call mem_alloc(Avalues,nelements)
    call mem_alloc(Bvalues,nelements)

    ! Open files for arrays
    call array4_open_file(A)
    call array4_open_file(B)


    ! Loop over index 2,3 and 4 and take the dot products
    idx4_loop: do idx4=1,A%dims(4)
       idx3_loop: do idx3=1,A%dims(3)
          idx2_loop: do idx2=1,A%dims(2)

             call array4_read_file_type3(A,idx2,idx3,idx4,Avalues,nelements)
             call array4_read_file_type3(B,idx2,idx3,idx4,Bvalues,nelements)
             res=res+dot_product(Avalues,Bvalues)

          end do idx2_loop
       end do idx3_loop
    end do idx4_loop



    ! Close files for arrays
    call array4_close_file(A,'keep')
    call array4_close_file(B,'keep')

    ! Deallocate memory
    call mem_dealloc(Avalues)
    call mem_dealloc(Bvalues)

  end subroutine array4_dotproduct_file_type3


  !> \brief Square of a 2norm
  function array4_norm(A) result(res)
    implicit none
    type(array4), intent(in) :: A
    real(realk) :: res
    real(realk) :: t0,t1
    call cpu_time(t0)

    res=0.0E0_realk

    ArrayOnFile: if(DECinfo%array4OnFile) then
       call array4_norm_file(A,res)
    else
       res=array4_dotproduct(A,A)
       call cpu_time(t1)
       time_array4_norm = time_array4_norm + (t1-t0)

    end if ArrayOnFile


    return
  end function array4_norm

  !> \brief Contraction of array4 and array2 summing the first index
  subroutine array4_contract1(A,B,C,zero)

    implicit none
    type(array4), intent(in) :: A
    type(array4), intent(inout) :: C
    type(array2), intent(in) :: B
    logical, intent(in) :: zero
    real(realk) :: MaxElement1,MaxElement2
    real(realk) :: beta
    real(realk) :: tcpu1,twall1,tcpu2,twall2
    integer :: k,l
    integer(kind=long) :: vec_size64
    integer :: vec_size

    call LSTIMER('START',tcpu1,twall1,DECinfo%output)

    if(A%dims(1) /= B%dims(1)) then
       print *,'index 1 -> ',A%dims(1),B%dims(1)
       stop 'Error :: Contraction array4_contract1 :: Dimensions do not match'
    end if
    vec_size64 = C%dims(1)*C%dims(2)*C%dims(3)*C%dims(4)
    if(vec_size64>MAXINT)then
       call lsquit('ERROR(array4_contract1): size of array cannot be &
                    &described by current integer type, please try another &
                    &compilation or fix this routine', DECinfo%output)
    endif
    vec_size = C%dims(1)*C%dims(2)*C%dims(3)*C%dims(4)
    if(zero) call ls_dzero(C%val,vec_size)

    if( (A%dims(2) /= C%dims(2)) .or. (A%dims(3) /= C%dims(3)) &
         & .or. (A%dims(4) /= C%dims(4)) ) then
       write(DECinfo%output,*) 'A%dims', A%dims
       write(DECinfo%output,*) 'C%dims', C%dims
       call lsquit('array4_contract1: A and C dimensions do not match',-1)
    end if

    do k=1,A%dims(4)
       do l=1,A%dims(3)
          call dgemm('t','n',B%dims(2),A%dims(2),B%dims(1), &
               1.0E0_realk,B%val,B%dims(1),A%val(:,:,l,k),A%dims(1), &
               1.0E0_realk,C%val(:,:,l,k),C%dims(1))
       end do
    end do

    call LSTIMER('START',tcpu2,twall2,DECinfo%output)

  end subroutine array4_contract1



  !> \brief Contraction of array4 and array2 summing the first index.
  !> A elements are read from file,
  !> whereas B is assumed to be available in memory.
  !> The result is written to the file belonging to C (but not to memory).
  !> \author Kasper Kristensen
  !> \date October 2010
  subroutine array4_contract1_file(A,B,C)

    implicit none
    type(array4), intent(in) :: A
    type(array4), intent(inout) :: C
    type(array2), intent(in) :: B
    real(realk) :: tcpu1,twall1,tcpu2,twall2

    call LSTIMER('START',tcpu1,twall1,DECinfo%output)


    ! Sanity checks
    ! *************

    if(A%dims(1) /= B%dims(1)) then
       print *,'index 1 -> ',A%dims(1),B%dims(1)
       stop 'Error :: Contraction array4_contract1_file :: Dimensions do not match'
    end if

    if(A%storing_type /= C%storing_type ) then
       write(DECinfo%output,*) 'Error in array4_contract1_file: &
            & A%storing_type /= C%storing_type'
       write(DECinfo%output,*) 'A%storing_type:', A%storing_type
       write(DECinfo%output,*) 'C%storing_type:', C%storing_type
       call lsquit('Error in array4_contract1_file: &
            & A%storing_type /= C%storing_type', DECinfo%output)
    end if


    StoringType: select case(A%storing_type)

    case(2)
       call array4_contract1_file_type2(A,B,C)
    case default
       call lsquit('array4_contract1_file: &
            & Requested storing type not implemented', DECinfo%output)

    end select StoringType

    call LSTIMER('START',tcpu2,twall2,DECinfo%output)

  end subroutine array4_contract1_file



  !> \brief Contraction of array4 and array2 summing the first index using
  !> storing type 2 (see type array4).
  !> A elements are read from file,
  !> whereas B is assumed to be available in memory.
  !> The result is written to the file belonging to C (but not to memory).
  !> \author Kasper Kristensen
  !> \date October 2010
  subroutine array4_contract1_file_type2(A,B,C)

    implicit none
    type(array4), intent(in) :: A
    type(array4), intent(inout) :: C
    type(array2), intent(in) :: B
    real(realk), pointer :: Avalues(:,:), Cvalues(:,:)
    integer :: k,l
    integer :: nelements
    integer(kind=long) :: nelements64

    nelements64 = C%dims(1)*C%dims(2)
    if(nelements64>MAXINT)then
       call lsquit('array4_contract1_file_type2: size of array cannot be &
                    &described by current integer type, please try another &
                    &compilation or fix this routine', DECinfo%output)
    endif
    nelements = C%dims(1)*C%dims(2)

    ! Sanity checks
    ! *************

    if(A%dims(1) /= B%dims(1)) then
       print *,'index 1 -> ',A%dims(1),B%dims(1)
       stop 'Error :: Contraction array4_contract1_file_type2 :: Dimensions do not match'
    end if

    if(A%storing_type /= C%storing_type ) then
       write(DECinfo%output,*) 'Error in array4_contract1_file: &
            & A%storing_type /= C%storing_type'
       write(DECinfo%output,*) 'A%storing_type:', A%storing_type
       write(DECinfo%output,*) 'C%storing_type:', C%storing_type
       call lsquit('Error in array4_contract1_file: &
            & A%storing_type /= C%storing_type', DECinfo%output)
    end if


    ! Initialize stuff
    ! ****************

    ! Temporary two-dimensional arrays to keep values for
    ! each set of (index 3, index 4).
    call mem_alloc(Avalues,A%dims(1),A%dims(2) )
    call mem_alloc(Cvalues,C%dims(1),C%dims(2) )

    ! Open file for arrays A
    call array4_open_file(A)

    ! Delete file for array C if present
    call array4_delete_file(C)
    ! Open new file for C
    call array4_open_file(C)



    ! Carry out transformations
    ! *************************

    !#ifdef USE_BLAS

    do k=1,A%dims(4)
       do l=1,A%dims(3)

          call ls_dzero(Cvalues,nelements)

          ! Read in the A(:,:,l,k) elements
          call array4_read_file_type2(A,l,k,Avalues,A%dims(1),A%dims(2))

          ! Calculate C(:,:,l,k) and store in Cvalues
          call dgemm('t','n',B%dims(2),A%dims(2),B%dims(1), &
               1.0E0_realk,B%val,B%dims(1),Avalues,A%dims(1), &
               1.0E0_realk,Cvalues,C%dims(1))

          ! Write to file the C(:,:,l,k) elements
          call array4_write_file_type2( C,l,k,Cvalues,C%dims(1),C%dims(2) )

       end do
    end do


    ! Free memory
    call mem_dealloc(Avalues)
    call mem_dealloc(Cvalues)

    ! Close file for array2
    call array4_close_file(A,'keep')
    call array4_close_file(C,'keep')


  end subroutine array4_contract1_file_type2




  !> \brief Contract with two indices (not first, but explicitly middle index)
  subroutine array4_contract2_middle(A,B,C)

    implicit none
    type(array4), intent(in) :: A,B
    type(array4), intent(inout) :: C
    integer :: i,j,k,l,m,n
    integer :: dim1,dim2,dim3
    real(realk) :: tcpu1,twall1,tcpu2,twall2

    call LSTIMER('START',tcpu1,twall1,DECinfo%output)

    ! contractions
    if(A%dims(3) /= B%dims(1) .or. A%dims(4) /= B%dims(2)) then
       do i=1,2
          print *,'index ',i+2,i,' -> ',A%dims(i+2),B%dims(i)
       end do
       stop 'error :: contraction 2 :: dimensions do not match'
    end if

    ! check the dimensions of the final array
    if(C%dims(1) /= A%dims(1) .or. C%dims(2) /= A%dims(2) .or. &
         C%dims(3) /= B%dims(3) .or. C%dims(4) /= B%dims(4)) then
       print *,'-- Final array --'
       print *,' 1 C ',C%dims(1), ' /= A(3) ',A%dims(1)
       print *,' 2 C ',C%dims(2), ' /= A(4) ',A%dims(2)
       print *,' 3 C ',C%dims(3), ' /= B(3) ',B%dims(3)
       print *,' 4 C ',C%dims(4), ' /= B(4) ',B%dims(4)
       stop 'error :: contraction 2 :: final arrays dimensions do not match'
    end if

    dim1 = A%dims(1)*A%dims(2)
    dim2 = A%dims(3)*A%dims(4)
    dim3 = B%dims(3)*B%dims(4)
    call dgemm('n','n',dim1,dim3,dim2, &
         1.0E0_realk,A%val,dim2,B%val,dim2,0.0E0_realk,C%val,dim1)

    call LSTIMER('START',tcpu2,twall2,DECinfo%output)

  end subroutine array4_contract2_middle

  !> \brief Contract with two indices (two first)
  subroutine array4_contract2(A,B,C)

    implicit none
    type(array4), intent(in) :: A,B
    type(array4), intent(inout) :: C
    integer :: i,j,k,l,m,n
    integer :: dim1,dim2,dim3
    real(realk) :: tcpu1,twall1,tcpu2,twall2

    call LSTIMER('START',tcpu1,twall1,DECinfo%output)

    ! contractions
    if(A%dims(1) /= B%dims(1) .or. A%dims(2) /= B%dims(2)) then
       do i=1,2
          print *,'index ',i,' -> ',A%dims(i),B%dims(i)
       end do
       stop 'error :: contraction 2 :: dimensions do not match'
    end if

    ! check the dimensions of the final array
    if(C%dims(1) /= A%dims(3) .or. C%dims(2) /= A%dims(4) .or. &
         C%dims(3) /= B%dims(3) .or. C%dims(4) /= B%dims(4)) then
       print *,'-- Final array --'
       print *,' 1 C ',C%dims(1), ' /= A(3) ',A%dims(3)
       print *,' 2 C ',C%dims(2), ' /= A(4) ',A%dims(4)
       print *,' 3 C ',C%dims(3), ' /= B(3) ',B%dims(3)
       print *,' 4 C ',C%dims(4), ' /= B(4) ',B%dims(4)
       stop 'error :: contraction 2 :: final arrays dimensions do not match'
    end if

    !#ifdef USE_BLAS
    dim1=A%dims(3)*A%dims(4)
    dim2=A%dims(1)*A%dims(2)
    dim3=B%dims(3)*B%dims(4)
    call dgemm('t','n',dim1,dim3,dim2, &
         1.0E0_realk,A%val,dim2,B%val,dim2,0.0E0_realk,C%val,dim1)
!!$#else
!!$    do n=1,B%dims(4)
!!$       do m=1,B%dims(3)
!!$          do l=1,A%dims(4)
!!$             do k=1,A%dims(3)
!!$
!!$                ! Contraction
!!$                do j=1,A%dims(2)
!!$                   do i=1,A%dims(1)
!!$                      C%val(k,l,m,n) = C%val(k,l,m,n) + &
!!$                           A%val(i,j,k,l)*B%val(i,j,m,n)
!!$                   end do
!!$                end do
!!$
!!$             end do
!!$          end do
!!$       end do
!!$    end do
!!$#endif

    call LSTIMER('START',tcpu2,twall2,DECinfo%output)

  end subroutine array4_contract2

  !> \brief Contract over both indices from array2
  subroutine array4_contract_array2(A,B,C)


    implicit none
    type(array4), intent(in) :: A
    type(array2), intent(in) :: B
    type(array2), intent(inout) :: C
    integer :: dim_con, dim_c1, dim_c2
    integer :: i
    real(realk) :: tcpu1,twall1,tcpu2,twall2

    call LSTIMER('START',tcpu1,twall1,DECinfo%output)

    ! sanity check
    if( A%dims(1) /= B%dims(1) .or. A%dims(2) /= B%dims(2) ) then
       do i=1,2
          print *,'index ',i,' -> ',A%dims(i),B%dims(i)
       end do
       stop 'error :: total contraction array4-array2 -> dimensions dont match'
    end if

    ! multiply
    dim_con = A%dims(1)*A%dims(2)
    dim_c1 = A%dims(3)*A%dims(4)
    dim_c2 = 1
    call mat_multiply(A%val,B%val,C%val,dim_con,dim_c1,dim_c2)

    call LSTIMER('START',tcpu2,twall2,DECinfo%output)

  end subroutine array4_contract_array2

  !> \brief Contract over 3 common indices
  subroutine array4_contract3(A,B,C)

    implicit none
    type(array4), intent(in) :: A,B
    type(array2), intent(inout) :: C
    integer :: dim_con, dim_c1, dim_c2
    integer :: i
    real(realk) :: tcpu1,twall1,tcpu2,twall2

    call LSTIMER('START',tcpu1,twall1,DECinfo%output)

    ! sanity check
    if( A%dims(1) /= B%dims(1) .or. A%dims(2) /= B%dims(2) &
         .or. A%dims(3) /= B%dims(3) ) then
       do i=1,3
          print *,'index ',i,' -> ',A%dims(i),B%dims(i)
       end do
       stop 'error :: contraction3 -> dimensions dont match'
    end if

    ! multiply
    dim_con = A%dims(1)*A%dims(2)*A%dims(3)
    dim_c1 = A%dims(4)
    dim_c2 = B%dims(4)
    call mat_multiply(A%val,B%val,C%val,dim_con,dim_c1,dim_c2)

    call LSTIMER('START',tcpu2,twall2,DECinfo%output)

  end subroutine array4_contract3

  !> \brief contract two matrices a(dim_con,dim_c1)*b(dim_con,dim_c2) over dim_con (the first index)
  subroutine mat_multiply(a,b,c,dim_con,dim_c1,dim_c2)
    implicit none
    integer, intent(in) :: dim_con,dim_c1,dim_c2
    real(realk), dimension(dim_con,dim_c1), intent(in) :: a
    real(realk), dimension(dim_con,dim_c2), intent(in) :: b
    real(realk), dimension(dim_c1,dim_c2), intent(inout) :: c
    real(realk) :: t0,t1
    call cpu_time(t0)
    !#ifdef USE_BLAS
    call dgemm('t','n',dim_c1,dim_c2,dim_con,1.0E0_realk,a,dim_con,&
         b,dim_con,0.0E0_realk,c,dim_c1)
!!$#else
!!$    c = matmul(transpose(a),b)
!!$#endif
    call cpu_time(t1)
    time_array4_mat_multiply = time_array4_mat_multiply + (t1-t0)
    return
  end subroutine mat_multiply

  !> \brief Print info about an array
  subroutine array4_info(array,output)

    implicit none
    type(array4), intent(in) :: array
    integer, intent(in) :: output
    real(realk) :: memory
    integer :: i

    write(DECinfo%output,'(/,a)')   '    Array4 Info     '
    write(DECinfo%output,'(a)')     '--------------------'
    write(DECinfo%output,'(a,i4)')  '    Dim1 : ',array%dims(1)
    write(DECinfo%output,'(a,i4)')  '    Dim2 : ',array%dims(2)
    write(DECinfo%output,'(a,i4)')  '    Dim3 : ',array%dims(3)
    write(DECinfo%output,'(a,i4)')  '    Dim4 : ',array%dims(4)
    write(DECinfo%output,'(a,4i4)') '   Order : ',(array%order(i),i=1,4)
    write(DECinfo%output,'(a,i4)')  'FileUnit : ',array%FUnit
    write(DECinfo%output,'(a,a)')   'FileName : ',array%FileName

    memory=( dble(array%dims(1))*dble(array%dims(2))*dble(array%dims(3))* &
         dble(array%dims(4))*8 )/dble(2**20)
    write(DECinfo%output,'(a,f10.2,a)') '  Memory : ',memory,' MB'

    return
  end subroutine array4_info

  !> \brief Reorder indices with additional memory allocation
  !> \author Marcin Ziolkowski
  subroutine array4_reorder(array,order)

    implicit none
    type(array4), intent(inout) :: array
    integer, dimension(4), intent(in) :: order
    integer, dimension(4) :: new_dims,new_order,order1,order2
    real(realk), pointer :: new_data(:,:,:,:)
    integer :: a,b,c,d
    integer :: dim1,dim2,dim3,dim4
    integer :: i,j
    integer :: aa,bb,cc,dd
    integer :: order_type,m,n
    real(realk) :: tcpu1,twall1,tcpu2,twall2

    call LSTIMER('START',tcpu1,twall1,DECinfo%output)

    ! get new order
    do i=1,4
!!$#ifdef VERBOSE
!!$       print *, array%order(i),'(',array%dims(i), &
!!$            ') => ',array%order(order(i)),' (',array%dims(order(i)),')'
!!$#endif
       new_order(i) = array%order(order(i))
       new_dims(i) = array%dims(order(i))
    end do

    ! Allocate space for reordered data

    ! Note: To make things work when array4_reorder is called in parallel
    !       we need to make a manual allocation of new_data (not using memory_allocate_4d).
    call mem_alloc(new_data,new_dims(1),new_dims(2),new_dims(3),new_dims(4) )
    call ls_dzero(new_data,size(new_data))

    call array_reorder_4d(1.0E0_realk,array%val,array%dims(1),array%dims(2),&
    &array%dims(3),array%dims(4),order,0.0E0_realk,new_data)

    call memory_deallocate_4d(array%val)
    array%dims=new_dims
    array%order(1)=1
    array%order(2)=2
    array%order(3)=3
    array%order(4)=4
    call memory_allocate_4d(array%val,array%dims)
    array%val = new_data
    call mem_dealloc(new_data)


    ! Set file address information according to new order (see array4_init_file)
    array%nelements = array%dims(1)*array%dims(2)*array%dims(3)
    if(DECinfo%array4OnFile) then
       if(associated(array%address)) then
          call mem_dealloc(array%address)
          nullify(array%address)
          call mem_alloc(array%address,1,1,1,array%dims(4) )
          array%address=0
       end if
    end if

    call LSTIMER('START',tcpu2,twall2,DECinfo%output)

  end subroutine array4_reorder



  subroutine array_reorder_4d_debug(tensor_in,dims,max_dim,order,tensor_out)
    implicit none
    integer :: max_dim
    real(realk), dimension(max_dim):: tensor_in, tensor_out
    integer, dimension(4), intent(in) :: order
    integer, dimension(4) :: new_order,order1,order2,dims
    integer :: a,b,c,d
    integer :: dim1,dim2,dim3,dim4,dim1b,dim2b,dim3b,vdim
    integer :: i,j,k,l
    integer :: aa,bb,cc,dd,block_size
    integer :: order_type,m,n
    real(realk) :: tcpu1,twall1,tcpu2,twall2

    call LSTIMER('START',tcpu1,twall1,DECinfo%output)


    dim1 = dims(1)
    dim2 = dims(2)*dims(1)
    dim3 = dims(3)*dims(2)*dims(1)


    do d = 1, dims(4),1
      do c = 1, dims(3),1
        do b = 1, dims(2),1
          do a = 1, dims(1),1
            ! old order
            order1 = [a,b,c,d]
            !print *, order1
            ! new order
            i = order1(order(1))
            j = order1(order(2))
            k = order1(order(3))
            l = order1(order(4))
            !print *,order2
            ! reorder
            dim1b=dims(order(1))
            dim2b=dims(order(1))*dims(order(2))
            dim3b=dims(order(1))*dims(order(2))*dims(order(3))
            tensor_out(i+(j-1)*dim1b+(k-1)*dim2b+(l-1)*dim3b) = tensor_in(a+(b-1)*dim1+(c-1)*dim2+(d-1)*dim3)
           enddo
         enddo
      enddo
    enddo

    call LSTIMER('START',tcpu2,twall2,DECinfo%output)

  end subroutine array_reorder_4d_debug

  !subroutine print_norm(p,x,d,e)
  !  integer :: d, e, i,j
  !  real(realk),intent(in) :: x(d,e),p
  !  real(realk) :: s
  !  i=0
  !  j=0
  !  s = 0.d0
  !  do i = 1, d,1
  !    do j = 1, e,1
 !      s = s +(p*x(i,j))**2
 !     enddo
 !   enddo
 !   print *, sqrt(s)
 ! end subroutine print_norm


  !> \brief Transpose data
!  subroutine mat_transp(mat_in,dim1,dim2,mat_out)
!    implicit none
!    integer,intent(in) :: dim1,dim2
!    real(realk), dimension(dim1,dim2), intent(in) :: mat_in
!    real(realk), dimension(dim2,dim1), intent(inout) :: mat_out
!    real(realk) :: t0,t1
!    integer :: i
!
!    call cpu_time(t0)
!    !#ifdef USE_BLAS
!    do i=1,dim2
!       call dcopy(dim1,mat_in(1,i),1,mat_out(i,1),dim2)
!    end do
!    !#else
!    !    mat_out = transpose(mat_in)
!    !#endif
!    call cpu_time(t1)
!    time_array4_mat_transpose = time_array4_mat_transpose + (t1-t0)
!
!    return
!  end subroutine mat_transp


  !> \brief partial transposition of multi index quantity
  subroutine partial_transpose_ip(x,m,n,o)
    integer,intent(in) :: m, n, o 
    real(realk), dimension(m*n*o),intent(inout) :: x
    integer :: i
    real(realk), dimension(m*n) :: y
    do i=1, o, 1
      call mat_transpose(m,n,1.0E0_realk,x((o-1)*m*n),0.0E0_realk,y)
      do j=1, m*n, 1
        x(j+(o-1)*m*n) = y(j)
      enddo
    enddo
  end subroutine partial_transpose_ip

!  !> \brief partial transposition of multi index quantity
!  subroutine partial_transpose(x,m,n,o,y)
!    integer,intent(in) :: m, n, o
!    real(realk), dimension(m*n*o),intent(inout) :: x
!    integer :: i
!    real(realk), dimension(m*n*o) :: y
!    do i=1, o, 1
!      call mat_transpose(x((o-1)*m*n),m,n,y((o-1)*m*n))
!    enddo
!  end subroutine partial_transpose
!
!  !> \author Patrick Ettenhuber
!  !> \brief Transpose data fast with adding to the destination matrix
!  subroutine mat_transpose_p(x,r,c,y)
!    integer,intent(in) :: r, c
!    real(realk), dimension(r,c), intent(in) :: x
!    real(realk), dimension(c,r), intent(inout) :: y
!    integer :: i,j, s, p, k, a, b,e,f
!    integer :: nr, nc, rr, rc, IOK, kb,d, iwrk
!    real(realk) :: val
!    real(realk), dimension(:)  , pointer :: sm
!    integer,     dimension(:)  , pointer :: h
!    kb = 512 !cache size of cpu --> /proc/cpuinfo
!    s = kb * 1024 / 8
!    k = int(0.9*sqrt(float(s)))
!    nr=r/k
!    nc=c/k
!    rr=mod(r,k)
!    rc=mod(c,k)
!    d = max(rr,rc)
!    iwrk = (d+k) / 4
!    call mem_alloc(sm,k*k)
!    call mem_alloc(h,iwrk)
!
!    do a = 0, nr-1, 1
!      do b = 0, nc-1, 1
!        do j=1,k,1
!          do i=1,k,1
!            sm((j-1)*k+i) = x(i+a*k,j+b*k)
!          enddo
!        enddo
!        call alg513(sm,k,k,k*k,h,iwrk,IOK)
!        do j=1,k,1
!          do i=1,k,1
!            y(i+b*k,j+a*k) = y(i+b*k,j+a*k) + sm((j-1)*k+i)
!          enddo
!        enddo
!      enddo
!    enddo
!
!    if (rc>0) then
!      iwrk = (rc+k)/4
!      do a = 0, nr-1, 1
!        do j=1,rc,1
!          do i=1,k,1
!            sm((j-1)*k+i) = x(i+a*k,j+nc*k)
!          enddo
!        enddo
!        call alg513(sm,k,rc,k*rc,h,iwrk,IOK)
!        do i=1,rc,1
!          do j=1,k,1
!            y(i+nc*k,j+a*k) = y(i+b*k,j+a*k) + sm((j-1)*rc+i)
!          enddo
!        enddo
!      enddo
!    endif
!
!    if (rr>0) then
!      iwrk = (rr+k)/4
!      do b = 0, nc-1, 1
!        do j=1,k,1
!          do i=1,rr,1
!            sm((j-1)*rr+i) = x(i+nr*k,j+b*k)
!          enddo
!        enddo
!        call alg513(sm,rr,k,k*rr,h,iwrk,IOK)
!        do j=1,rr,1
!          do i=1,k,1
!            y(i+b*k,j+nr*k) = y(i+b*k,j+nr*k) + sm((j-1)*k+i)
!          enddo
!        enddo
!      enddo
!    endif
!
!    if (rr>0 .and. rc>0) then
!      iwrk = (rc+rr)/4
!      do j=1,rc,1
!        do i=1,rr,1
!          sm((j-1)*rr+i) = x(i+nr*k,j+nc*k)
!        enddo
!      enddo
!      call alg513(sm,rr,rc,rc*rr,h,iwrk,IOK)
!      do j=1,rr,1
!        do i=1,rc,1
!          y(i+nc*k,j+nr*k) = y(i+nc*k,j+nr*k) + sm((j-1)*rc+i)
!        enddo
!      enddo
!    endif
!    call mem_dealloc(h)
!    call mem_dealloc(sm)
!  end subroutine mat_transpose_p
!
!
!  !> \author Patrick Ettenhuber
!  !> \brief Transpose data fast
!  subroutine mat_transpose(x,r,c,y)
!    integer,intent(in) :: r, c
!    real(realk), dimension(r,c), intent(in) :: x
!    real(realk), dimension(c,r), intent(inout) :: y
!    integer :: i,j, s, p, k, a, b,e,f
!    integer :: nr, nc, rr, rc, IOK, kb,d, iwrk
!    real(realk) :: val
!    real(realk), dimension(:)  , pointer :: sm
!    integer,     dimension(:)  , pointer :: h
!    kb = 512 !cache size of cpu --> /proc/cpuinfo
!    s = kb * 1024 / 8
!    k = int(0.9*sqrt(float(s)))
!    nr=r/k
!    nc=c/k
!    rr=mod(r,k)
!    rc=mod(c,k)
!    d = max(rr,rc)
!    iwrk = (d+k) / 4
!    call mem_alloc(sm,k*k)
!    call mem_alloc(h,iwrk)
!
!    do a = 0, nr-1, 1
!      do b = 0, nc-1, 1
!        do j=1,k,1
!          do i=1,k,1
!            sm((j-1)*k+i) = x(i+a*k,j+b*k)
!          enddo
!        enddo
!        call alg513(sm,k,k,k*k,h,iwrk,IOK)
!        do j=1,k,1
!          do i=1,k,1
!            y(i+b*k,j+a*k) = sm((j-1)*k+i)
!          enddo
!        enddo
!      enddo
!    enddo
!
!    if (rc>0) then
!      iwrk = (rc+k)/4
!      do a = 0, nr-1, 1
!        do j=1,rc,1
!          do i=1,k,1
!            sm((j-1)*k+i) = x(i+a*k,j+nc*k)
!          enddo
!        enddo
!        call alg513(sm,k,rc,k*rc,h,iwrk,IOK)
!        do i=1,rc,1
!          do j=1,k,1
!            y(i+nc*k,j+a*k) = sm((j-1)*rc+i)
!          enddo
!        enddo
!      enddo
!    endif
!
!    if (rr>0) then
!      iwrk = (rr+k)/4
!      do b = 0, nc-1, 1
!        do j=1,k,1
!          do i=1,rr,1
!            sm((j-1)*rr+i) = x(i+nr*k,j+b*k)
!          enddo
!        enddo
!        call alg513(sm,rr,k,k*rr,h,iwrk,IOK)
!        do j=1,rr,1
!          do i=1,k,1
!            y(i+b*k,j+nr*k) = sm((j-1)*k+i)
!          enddo
!        enddo
!      enddo
!    endif
!
!    if (rr>0 .and. rc>0) then
!      iwrk = (rc+rr)/4
!      do j=1,rc,1
!        do i=1,rr,1
!          sm((j-1)*rr+i) = x(i+nr*k,j+nc*k)
!        enddo
!      enddo
!      call alg513(sm,rr,rc,rc*rr,h,iwrk,IOK)
!      do j=1,rr,1
!        do i=1,rc,1
!          y(i+nc*k,j+nr*k) = sm((j-1)*rc+i)
!        enddo
!      enddo
!    endif
!    call mem_dealloc(h)
!    call mem_dealloc(sm)
!  end subroutine mat_transpose



  !> \brief Print statistics of array4 objects
  subroutine array4_print_statistics(output)

    implicit none
    integer, intent(in) :: output

    write(DECinfo%output,'(/,a)')    '  Array4 statistics       '
    write(DECinfo%output,'(a)')      ' =================================================='
    write(DECinfo%output,'(a,i8)')   ' Number of created arrays   : ',CreatedArrays
    write(DECinfo%output,'(a,i8)')   ' Number of destroyed arrays : ',DestroyedArrays
    write(DECinfo%output,'(a,i8)')   ' Orphaned arrays            : ',CreatedArrays-DestroyedArrays

    write(DECinfo%output,'(/,a)')    '  Array4 time statistics     '
    write(DECinfo%output,'(a)')      ' =================================================='
    write(DECinfo%output,'(a,f22.2,a)') ' array4_init              :',time_array4_init            , 's '
    write(DECinfo%output,'(a,f22.2,a)') ' array4_free              :',time_array4_free            , 's '
    write(DECinfo%output,'(a,f22.2,a)') ' array4_duplicate         :',time_array4_duplicate       , 's '
    write(DECinfo%output,'(a,f22.2,a)') ' array4_alloc          :',time_array4_alloc        , 's '
    write(DECinfo%output,'(a,f22.2,a)') ' array4_dealloc        :',time_array4_dealloc      , 's '
    write(DECinfo%output,'(a,f22.2,a)') ' array4_write_data        :',time_array4_write_data      , 's '
    write(DECinfo%output,'(a,f22.2,a)') ' array4_read_data         :',time_array4_read_data       , 's '
    write(DECinfo%output,'(a,f22.2,a)') ' array4_write             :',time_array4_write           , 's '
    write(DECinfo%output,'(a,f22.2,a)') ' array4_read              :',time_array4_read            , 's '
    write(DECinfo%output,'(a,f22.2,a)') ' array4_add               :',time_array4_add             , 's '
    write(DECinfo%output,'(a,f22.2,a)') ' array4_add_to            :',time_array4_add_to          , 's '
    write(DECinfo%output,'(a,f22.2,a)') ' array4_scale             :',time_array4_scale           , 's '
    write(DECinfo%output,'(a,f22.2,a)') ' array4_dotproduct        :',time_array4_dotproduct      , 's '
    write(DECinfo%output,'(a,f22.2,a)') ' array4_norm              :',time_array4_norm            , 's '
    write(DECinfo%output,'(a,f22.2,a)') ' array4_contract1         :',time_array4_contract1       , 's '
    write(DECinfo%output,'(a,f22.2,a)') ' array4_contract2         :',time_array4_contract2       , 's '
    write(DECinfo%output,'(a,f22.2,a)') ' array4_contract2_middle  :',time_array4_contract2_middle, 's '
    write(DECinfo%output,'(a,f22.2,a)') ' array4_contract3         :',time_array4_contract3       , 's '
    write(DECinfo%output,'(a,f22.2,a)') ' array4_contract_array2   :',time_array4_contract_array2 , 's '
    write(DECinfo%output,'(a,f22.2,a)') ' array4_mat_multiply      :',time_array4_mat_multiply    , 's '
    write(DECinfo%output,'(a,f22.2,a)') ' array4_reorder           :',time_array4_reorder         , 's '

    write(DECinfo%output,'(a,f22.2,a)') ' array4_contract_with_ao  :',time_array4_contract_with_ao, 's '

    return
  end subroutine array4_print_statistics


  !> \brief Open file for array A using file handling in C.
  !> Important: A must be initialized before this call!
  !> \author Kasper Kristensen
  !> \date October 2010
  subroutine array4_open_file(A)
    implicit none
    !> Array for which file is opened
    type(array4),intent(in) :: A
    integer :: funit

    call openfile(A%funit,A%Filename)

  end subroutine array4_open_file


  !> \brief Close file for array A using file handling in C.
  !> Important: A must be initialized before this call!
  !> \author Kasper Kristensen
  !> \date October 2010
  subroutine array4_close_file(A,keep_or_delete)
    implicit none
    !> Array for which file is closed
    type(array4) :: A
    !> Status for file? 'KEEP' or 'DELETE'
    character(*), intent(in) :: keep_or_delete

    call closefile(A%Funit,keep_or_delete)

  end subroutine array4_close_file


  !> \brief Delete file for array A using file handling in C.
  !> \author Kasper Kristensen
  !> \date October 2010
  subroutine array4_delete_file(A)
    implicit none
    !> Array for which file is deleted
    type(array4) :: A
    logical :: file_exists

    ! Check that file for array A exists at all
    inquire(file=A%Filename,exist=file_exists)

    if(file_exists) then
       call array4_open_file(A)
       call array4_close_file(A,'delete')
    end if

  end subroutine array4_delete_file




  !> \brief Write array4 values to file using storing type 1 (see type array4_init_file).
  !> I.e. elements A(:,:,:,i) for a given "i" are written to file.
  !> \author Kasper Kristensen
  !> \date February 2011
  subroutine array4_write_file_type1(A,idx4,values,dim1,dim2,dim3)
    implicit none
    !> Full array information (but no memory allocated for values)
    type(array4), intent(inout) :: A
    !> Fourth index in full array that we want to write to file
    integer, intent(in) :: idx4
    !> Dimension 1 for array
    integer,intent(in) :: dim1
    !> Dimension 2 for array
    integer,intent(in) :: dim2
    !> Dimension 3 for array
    integer,intent(in) :: dim3
    !> Values to write to file
    real(realk) :: values(dim1,dim2,dim3)


    ! Sanity checks
    ! *************

    ! 1. idx4 is smaller than or equal to fourth dimension of A
    if(idx4 > A%dims(4)) then
       call lsquit('array4_write_file_type1: idx4 > A%dims(4)', DECinfo%output)
    end if

    ! 2. Number of elements to write to disk must equal dims(1)*dims(2)*dims(3)
    if( A%nelements /= dim1*dim2*dim3 ) then
       call lsquit('array4_write_file_type1: &
            & A%nelements /= dim1*dim2*dim3', DECinfo%output)
    end if

    ! 3. Consistent input values
    if( A%dims(1) /= dim1 .or. A%dims(2) /= dim2 .or. A%dims(3) /= dim3 ) then
       call lsquit('array4_write_file_type1: &
            & Array dimensions inconsistent with input values', DECinfo%output)
    end if

    ! 4. Address counter value is meaningful
    if(A%address_counter < 1) then
       call lsquit('array4_write_file_type1: A%address_counter < 1', DECinfo%output)
    end if


    ! Write values to file at given address
    ! *************************************
    call writevector( A%Funit,A%address_counter,A%nelements,values )

    ! Bookkeeping of which array values are stored where
    A%address(1,1,1,idx4) = A%address_counter

    ! Increase address counter
    A%address_counter = A%address_counter + A%nelements



  end subroutine array4_write_file_type1



  !> \brief Write array4 values to file using storing type 2 (see type array4_init_file).
  !> I.e. for given values of the third and fourth array indices
  !> the input values are written to file.
  !> \author Kasper Kristensen
  !> \date October 2010
  subroutine array4_write_file_type2(A,idx3,idx4,values,dim1,dim2)
    implicit none
    !> Full array information (but no memory allocated for values)
    type(array4), intent(inout) :: A
    !> Third index in full array that we want to write to file
    integer, intent(in) :: idx3
    !> Fourth index in full array that we want to write to file
    integer, intent(in) :: idx4
    !> Number of rows for values array
    integer,intent(in) :: dim1
    !> Number of colums for values array
    integer,intent(in) :: dim2
    !> Values to write to file
    real(realk) :: values(dim1,dim2)


    ! Sanity checks
    ! *************

    ! 1. idx3 and idx4 are smaller than or equal to dimensions 3 and 4.
    if(idx3 > A%dims(3)) then
       call lsquit('array4_write_file_type2: dx3 > A%dims(3)', DECinfo%output)
    end if
    if(idx4 > A%dims(4)) then
       call lsquit('array4_write_file_type2: idx4 > A%dims(4)', DECinfo%output)
    end if

    ! 2. Number of elements to write to disk must equal dims(1)*dims(2)
    if( A%nelements /= dim1*dim2 ) then
       call lsquit('array4_write_file_type2: &
            & A%nelements /= dim1*dim2', DECinfo%output)
    end if

    ! 3. Consistent input values
    if( A%dims(1) /= dim1 .or. A%dims(2) /= dim2 ) then
       call lsquit('array4_write_file_type2: &
            & Array dimensions inconsistent with input values', DECinfo%output)
    end if

    ! 4. Address counter value is meaningful
    if(A%address_counter < 1) then
       call lsquit('array4_write_file_type2: A%address_counter < 1', DECinfo%output)
    end if


    ! Write values to file at given address
    ! *************************************
    call writevector( A%Funit,A%address_counter,A%nelements,values )

    ! Bookkeeping of which array values are stored where
    A%address(1,1,idx3,idx4) = A%address_counter

    ! Increase address counter
    A%address_counter = A%address_counter + A%nelements



  end subroutine array4_write_file_type2



  !> \brief Write array4 values to file using storing type 3 (see type array4_init_file).
  !> I.e. for given values of second, third and fourth array indices
  !> the full set of values for the first index is written to file.
  !> \author Kasper Kristensen
  !> \date October 2010
  subroutine array4_write_file_type3(A,idx2,idx3,idx4,values,dim1)
    implicit none
    !> Full array information (but no memory allocated for values)
    type(array4), intent(inout) :: A
    !> Second index in full array that we want to write to file
    integer, intent(in) :: idx2
    !> Third index in full array that we want to write to file
    integer, intent(in) :: idx3
    !> Fourth index in full array that we want to write to file
    integer, intent(in) :: idx4
    !> Number of elements in values array
    integer,intent(in) :: dim1
    !> Values to write to file
    real(realk) :: values(dim1)


    ! Sanity checks
    ! *************

    ! 1. idx2, idx3, and idx4 are smaller than or equal to dimensions 3 and 4.
    if(idx2 > A%dims(2)) then
       call lsquit('array4_write_file_type3: dx2 > A%dims(2)', DECinfo%output)
    end if
    if(idx3 > A%dims(3)) then
       call lsquit('array4_write_file_type3: dx3 > A%dims(3)', DECinfo%output)
    end if
    if(idx4 > A%dims(4)) then
       call lsquit('array4_write_file_type3: idx4 > A%dims(4)', DECinfo%output)
    end if

    ! 2. Number of elements to write to disk must equal dim1
    if( A%nelements /= dim1 ) then
       call lsquit('array4_write_file_type3: &
            & A%nelements /= dim1', DECinfo%output)
    end if

    ! 3. Consistent input values
    if( A%dims(1) /= dim1 ) then
       call lsquit('array4_write_file_type3: &
            & Array dimensions inconsistent with input values', DECinfo%output)
    end if

    ! 4. Address counter value is meaningful
    if(A%address_counter < 1) then
       call lsquit('array4_write_file_type3: A%address_counter < 1', DECinfo%output)
    end if


    ! Write values to file at given address
    ! *************************************
    call writevector( A%Funit,A%address_counter,A%nelements,values )

    ! Bookkeeping of which array values are stored where
    A%address(1,idx2,idx3,idx4) = A%address_counter

    ! Increase address counter
    A%address_counter = A%address_counter + A%nelements



  end subroutine array4_write_file_type3



  !> \brief Read array4 values from file using storing type 1 (see type array4).
  !> I.e. elements A(:,:,:,i) for a given "i" are read from file.
  !> \author Kasper Kristensen
  !> \date February 2011
  subroutine array4_read_file_type1(A,idx4,values,dim1,dim2,dim3)
    implicit none
    !> Full array information (but no memory allocated for values)
    type(array4), intent(in) :: A
    !> Fourth index in full array that we want to read from file
    integer, intent(in) :: idx4
    !> Dimension 1 for array values
    integer,intent(in) :: dim1
    !> Dimension 2 for array values
    integer,intent(in) :: dim2
    !> Dimension 3 for array values
    integer,intent(in) :: dim3
    !> Values to read from file
    real(realk) :: values(dim1,dim2,dim3)
    integer(kind=long) :: address

    ! Initialize
    values=0E0_realk

    ! Set address on file based on input indices
    address = A%address(1,1,1,idx4)


    ! Sanity checks
    ! *************

    ! 1. idx4 is smaller than or equal to fourth dimension of A.
    if(idx4 > A%dims(4)) then
       call lsquit('array4_read_file_type1: idx4 > A%dims(4)', DECinfo%output)
    end if

    ! 2. Number of elements to read from disk must equal dims(1)*dims(2)dims(3)
    if( A%nelements /= dim1*dim2*dim3 ) then
       call lsquit('array4_read_file_type1: &
            & A%nelements /= dim1*dim2*dim3', DECinfo%output)
    end if

    ! 3. Consistent input values
    if( A%dims(1) /= dim1 .or. A%dims(2) /= dim2 .or. A%dims(3) /= dim3 ) then
       call lsquit('array4_read_file_type1: &
            & Array dimensions inconsistent with input values', DECinfo%output)
    end if

    ! 4. Address value is meaningful
    if(address < 1) then
       call lsquit('array4_read_file_type1: A%address_counter < 1', DECinfo%output)
    end if

    ! Read values (:,:,idx3,idx4) to file at given address
    ! *****************************************************
    call readvector( A%Funit,address,A%nelements,values )


  end subroutine array4_read_file_type1




  !> \brief Read array4 values from file using storing type 2 (see type array4).
  !> I.e. for given values of the third and fourth array indices
  !> the corresponding values are read from file into the values array.
  !> \author Kasper Kristensen
  !> \date October 2010
  subroutine array4_read_file_type2(A,idx3,idx4,values,dim1,dim2)
    implicit none
    !> Full array information (but no memory allocated for values)
    type(array4), intent(in) :: A
    !> Third index in full array that we want to read from file
    integer, intent(in) :: idx3
    !> Fourth index in full array that we want to read from file
    integer, intent(in) :: idx4
    !> Number of rows for values array
    integer,intent(in) :: dim1
    !> Number of colums for values array
    integer,intent(in) :: dim2
    !> Values to read from file
    real(realk) :: values(dim1,dim2)
    integer(kind=long) :: address

    ! Initialize
    values=0E0_realk


    ! Set address on file based on input indices
    address = A%address(1,1,idx3,idx4)


    ! Sanity checks
    ! *************

    ! 1. idx3 and idx4 are smaller than or equal to dimensions 3 and 4.
    if(idx3 > A%dims(3)) then
       call lsquit('array4_read_file_type2: idx3 > A%dims(3)', DECinfo%output)
    end if
    if(idx4 > A%dims(4)) then
       call lsquit('array4_read_file_type2: idx4 > A%dims(4)', DECinfo%output)
    end if

    ! 2. Number of elements to read from disk must equal dims(1)*dims(2)
    if( A%nelements /= dim1*dim2 ) then
       call lsquit('array4_read_file_type2: &
            & A%nelements /= dim1*dim2', DECinfo%output)
    end if

    ! 3. Consistent input values
    if( A%dims(1) /= dim1 .or. A%dims(2) /= dim2 ) then
       call lsquit('array4_read_file_type2: &
            & Array dimensions inconsistent with input values', DECinfo%output)
    end if

    ! 4. Address value is meaningful
    if(address < 1) then
       call lsquit('array4_read_file_type2: A%address_counter < 1', DECinfo%output)
    end if

    ! Read values (:,:,idx3,idx4) to file at given address
    ! *****************************************************
    call readvector( A%Funit,address,A%nelements,values )


  end subroutine array4_read_file_type2




  !> \brief Read array4 values from file using storing type 3 (see type array4).
  !> I.e. for given values of the second, third, and fourth array indices
  !> the corresponding values for the first index are read from file.
  !> \author Kasper Kristensen
  !> \date October 2010
  subroutine array4_read_file_type3(A,idx2,idx3,idx4,values,dim1)
    implicit none
    !> Full array information (but no memory allocated for values)
    type(array4), intent(in) :: A
    !> Second index in full array that we want to read from file
    integer, intent(in) :: idx2
    !> Third index in full array that we want to read from file
    integer, intent(in) :: idx3
    !> Fourth index in full array that we want to read from file
    integer, intent(in) :: idx4
    !> Number of elements in values vector
    integer,intent(in) :: dim1
    !> Values to read from file
    real(realk) :: values(dim1)
    integer(kind=long) :: address

    ! Initialize
    values=0E0_realk

    ! Set address on file based on input indices
    address = A%address(1,idx2,idx3,idx4)


    ! Sanity checks
    ! *************

    ! 1. idx2, idx3 and idx4 are smaller than or equal to dimensions 3 and 4.
    if(idx2 > A%dims(2)) then
       call lsquit('array4_read_file_type3: idx3 > A%dims(3)', DECinfo%output)
    end if
    if(idx3 > A%dims(3)) then
       call lsquit('array4_read_file_type3: idx3 > A%dims(3)', DECinfo%output)
    end if
    if(idx4 > A%dims(4)) then
       call lsquit('array4_read_file_type3: idx4 > A%dims(4)', DECinfo%output)
    end if

    ! 2. Number of elements to read from disk must equal dim1
    if( A%nelements /= dim1 ) then
       call lsquit('array4_read_file_type3: &
            & A%nelements /= dim1*dim2', DECinfo%output)
    end if

    ! 3. Consistent input values
    if( A%dims(1) /= dim1 ) then
       call lsquit('array4_read_file_type3: &
            & Array dimensions inconsistent with input values', DECinfo%output)
    end if

    ! 4. Address value is meaningful
    if(address < 1) then
       call lsquit('array4_read_file_type3: A%address_counter < 1', DECinfo%output)
    end if

    ! Read values (:,:,idx3,idx4) to file at given address
    ! *****************************************************
    call readvector( A%Funit,address,A%nelements,values )


  end subroutine array4_read_file_type3





  !> \brief Get norm of array4 when the elements are stored on file and not in memory.
  !> Norm = dot(A,A)   (We do NOT take the square root!)
  !> Thus, here we assume that all values have been written to file!
  subroutine array4_norm_file(A,res)
    implicit none
    !> Array for which norm is requested
    type(array4), intent(in) :: A
    !> Norm of array
    real(realk), intent(inout) :: res

    StoringType: select case(A%storing_type)

    case(1)
       call array4_norm_file_type1(A,res)

    case(2)
       call array4_norm_file_type2(A,res)

    case(3)
       call array4_norm_file_type3(A,res)

    case default
       write(DECinfo%output,*) 'array4_norm_file: Requested storing type is not &
            & implemented:', A%storing_type
       call lsquit('array4_norm_file: Requested storing type is not &
            & implemented.', DECinfo%output)

    end select StoringType


  end subroutine array4_norm_file



  !> \brief Get norm of array4 when the elements are stored on file
  !> using storing type 1 (see type array4).
  !> Norm = dot(A,A)  (We do NOT take the square root!)
  !> Thus, here we assume that all values have been written to file!
  subroutine array4_norm_file_type1(A,res)
    implicit none
    !> Array for which norm is requested
    type(array4), intent(in) :: A
    !> Norm of array
    real(realk), intent(inout) :: res
    real(realk),pointer :: values(:,:,:)
    integer(kind=long) :: nelements64
    integer :: nelements
    integer :: idx4,dim1,dim2,dim3,dim4,i
    real(realk), external :: ddot

    ! Initialize stuff
    res = 0E0_realk
    dim1 = A%dims(1)
    dim2 = A%dims(2)
    dim3 = A%dims(3)
    dim4 = A%dims(4)
    nelements64=dim1*dim2*dim3
    if(nelements64>MAXINT)then
       call lsquit('ERROR(array4_norm_file_type1): size of array cannot be &
                    &described by current integer type, please try another &
                    &compilation or fix this routine', DECinfo%output)
    endif
    nelements=dim1*dim2*dim3
    call mem_alloc(values,dim1,dim2,dim3)

    ! Open file for array
    call array4_open_file(A)

    ! Loop over index 4 and take the dot products
    idx4_loop: do idx4=1,dim4

       call array4_read_file_type1(A,idx4,values,dim1,dim2,dim3)

       !#ifdef USE_BLAS
       res = res + ddot(nelements,values,1,values,1)
       !#else
       !          do i=1,dim2
       !             res=res+dot_product(values(:,i),values(:,i))
       !          end do
       !#endif

    end do idx4_loop

    ! Close file for array
    call array4_close_file(A,'keep')

    ! Deallocate memory
    call mem_dealloc(values)

  end subroutine array4_norm_file_type1



  !> \brief Get norm of array4 when the elements are stored on file
  !> using storing type 2 (see type array4).
  !> Norm = dot(A,A)  (We do NOT take the square root!)
  !> Thus, here we assume that all values have been written to file!
  subroutine array4_norm_file_type2(A,res)
    implicit none
    !> Array for which norm is requested
    type(array4), intent(in) :: A
    !> Norm of array
    real(realk), intent(inout) :: res
    real(realk),pointer :: values(:,:)
    integer :: idx3,idx4,dim1,dim2,nelements,i
    real(realk), external :: ddot

    ! Initialize stuff
    res = 0E0_realk
    dim1 = A%dims(1)
    dim2 = A%dims(2)
    nelements=dim1*dim2
    call mem_alloc(values,A%dims(1),A%dims(2) )

    ! Open file for array
    call array4_open_file(A)

    ! Loop over index 3 and 4 and take the dot products
    idx4_loop: do idx4=1,A%dims(4)
       idx3_loop: do idx3=1,A%dims(3)


          call array4_read_file_type2(A,idx3,idx4,values,dim1,dim2)

          !#ifdef USE_BLAS
          res = res + ddot(nelements,values,1,values,1)
          !#else
          !          do i=1,dim2
          !             res=res+dot_product(values(:,i),values(:,i))
          !          end do
          !#endif

       end do idx3_loop
    end do idx4_loop



    ! Close file for array
    call array4_close_file(A,'keep')

    ! Deallocate memory
    call mem_dealloc(values)

  end subroutine array4_norm_file_type2


  !> \brief Get norm of array4 when the elements are stored on file
  !> using storing type 3 (see type array4_init_file).
  !> Thus, here we assume that all values have been written to file!
  !> Norm = dot(A,A)  (We do NOT take the square root!)
  subroutine array4_norm_file_type3(A,res)
    implicit none
    !> Array for which norm is requested
    type(array4), intent(in) :: A
    !> Norm of array
    real(realk), intent(inout) :: res
    real(realk),pointer :: values(:)
    integer :: idx2, idx3,idx4,nelements,i

    ! Initialize stuff
    res = 0E0_realk
    nelements = A%dims(1)
    call mem_alloc(values,nelements)
    ! Open file for array
    call array4_open_file(A)


    ! Loop over index 2,3 and 4 and take the dot products
    idx4_loop: do idx4=1,A%dims(4)
       idx3_loop: do idx3=1,A%dims(3)
          idx2_loop: do idx2=1,A%dims(2)

             call array4_read_file_type3(A,idx2,idx3,idx4,values,nelements)
             res=res+dot_product(values,values)

          end do idx2_loop
       end do idx3_loop
    end do idx4_loop



    ! Close file for array
    call array4_close_file(A,'keep')

    ! Deallocate memory
    call mem_dealloc(values)

  end subroutine array4_norm_file_type3


  !> \brief Set all elements in A to zero,
  !> assuming that the A elements are stored on file.
  !> \author Kasper Kristensen
  !> \date October 2010
  subroutine array4_zero_file(A)
    implicit none
    type(array4), intent(inout) :: A

    StoringType: select case(A%storing_type)

    case(2)
       call array4_zero_file_type2(A)

    case(3)
       call array4_zero_file_type3(A)

    case default

       call lsquit('array4_zero_file: Requested storing type not implemented!', DECinfo%output)

    end select StoringType

  end subroutine array4_zero_file


  !> \brief Set all elements of array A to zero and write
  !> to file using storing type 3 (see type array4_init_file).
  !> \author Kasper Kristensen
  !> \date October 2010
  subroutine array4_zero_file_type2(A)
    implicit none
    !> Full array information (but no memory allocated for values)
    type(array4), intent(inout) :: A
    !> Values to write to file
    real(realk),pointer :: values(:,:)
    integer :: idx3,idx4

    ! Sanity check
    ! ************

    ! Number of elements to write to disk must equal dims(1)*dims(2)
    if( A%nelements /= A%dims(1)*A%dims(2) ) then
       call lsquit('array4_zero_file_type2: &
            &  A%nelements /= dim1*dim2', DECinfo%output)
    end if



    ! Initialize stuff
    ! ****************

    ! Set address counter to 1
    A%address_counter = 1

    ! Allocate values and set to zero
    call mem_alloc(values,A%dims(1),A%dims(2) )
    values=0.0E0_realk


    ! Open file
    call array4_open_file(A)


    ! Write values (all equal zero) to file
    ! *************************************

    do idx4=1,A%dims(4)
       do idx3=1,A%dims(3)
          call array4_write_file_type2(A,idx3,idx4,values,A%dims(1),A%dims(2))
       end do
    end do

    ! Set address counter back to 1
    A%address_counter = 1

    ! Free values and close file
    call mem_dealloc(values)
    call array4_close_file(A,'keep')


  end subroutine array4_zero_file_type2



  !> \brief Set all elements of array A to zero and write
  !> to file using storing type 2 (see type array4_init_file).
  !> \author Kasper Kristensen
  !> \date October 2010
  subroutine array4_zero_file_type3(A)
    implicit none
    !> Full array information (but no memory allocated for values)
    type(array4), intent(inout) :: A
    !> Values to write to file
    real(realk),pointer :: values(:)
    integer :: idx2,idx3,idx4

    ! Sanity check
    ! ************

    ! Number of elements to write to disk must equal dims(1)
    if( A%nelements /= A%dims(1) ) then
       call lsquit('array4_zero_file_type2: &
            & A%nelements /= A%dims(1)', DECinfo%output)
    end if



    ! Initialize stuff
    ! ****************

    ! Set address counter to 1
    A%address_counter = 1

    ! Allocate values and set to zero
    call mem_alloc(values,A%dims(1) )
    values=0.0E0_realk


    ! Open file
    call array4_open_file(A)


    ! Write values (all equal zero) to file
    ! *************************************

    do idx4=1,A%dims(4)
       do idx3=1,A%dims(3)
          do idx2=1,A%dims(2)
             call array4_write_file_type3(A,idx2,idx3,idx4,values,A%dims(1))
          end do
       end do
    end do

    ! Free values and close file
    call mem_dealloc(values)
    call array4_close_file(A,'keep')

    ! Set address counter back to 1
    A%address_counter = 1


  end subroutine array4_zero_file_type3





  !> \brief Set all elements of array A to one and write
  !> to file using storing type 2 (see type array4).
  !> ONLY FOR DEBUGGING PURPOSES!!!
  !> \author Kasper Kristensen
  !> \date October 2010
  subroutine array4_zeroone_file_type2(A)
    implicit none
    !> Full array information (but no memory allocated for values)
    type(array4), intent(inout) :: A
    !> Values to write to file
    real(realk),pointer :: values(:,:)
    integer :: idx3,idx4

    ! Sanity check
    ! ************

    ! Number of elements to write to disk must equal dims(1)*dims(2)
    if( A%nelements /= A%dims(1)*A%dims(2) ) then
       call lsquit('array4_zero_file_type2: &
            & A%nelements /= dim1*dim2', DECinfo%output)
    end if



    ! Initialize stuff
    ! ****************

    ! Set address counter to 1
    A%address_counter = 1

    ! Allocate values and set to zero
    call mem_alloc(values,A%dims(1),A%dims(2) )
    values=1.0E0_realk


    ! Open file
    call array4_open_file(A)


    ! Write values (all equal zero) to file
    ! *************************************

    do idx4=1,A%dims(4)
       do idx3=1,A%dims(3)
          call array4_write_file_type2(A,idx3,idx4,values,A%dims(1),A%dims(2))
       end do
    end do

    ! Free values and close file
    call mem_dealloc(values)
    call array4_close_file(A,'keep')


  end subroutine array4_zeroone_file_type2



  !> \brief Set all elements of array A to one and write
  !> to file using storing type 3 (see type array4_init_file).
  !> ONLY FOR DEBUGGING PURPOSES!!!
  !> \author Kasper Kristensen
  !> \date October 2010
  subroutine array4_zeroone_file_type3(A)
    implicit none
    !> Full array information (but no memory allocated for values)
    type(array4), intent(inout) :: A
    !> Values to write to file
    real(realk),pointer :: values(:)
    integer :: idx2,idx3,idx4

    ! Sanity check
    ! ************

    ! Number of elements to write to disk must equal dims(1)
    if( A%nelements /= A%dims(1) ) then
       call lsquit('array4_zero_file_type2: &
            &  A%nelements /= A%dims(1)', DECinfo%output)
    end if



    ! Initialize stuff
    ! ****************

    ! Set address counter to 1
    A%address_counter = 1

    ! Allocate values and set to zero
    call mem_alloc(values,A%dims(1) )
    values=1.0E0_realk


    ! Open file
    call array4_open_file(A)


    ! Write values (all equal zero) to file
    ! *************************************

    do idx4=1,A%dims(4)
       do idx3=1,A%dims(3)
          do idx2=1,A%dims(2)
             call array4_write_file_type3(A,idx2,idx3,idx4,values,A%dims(1))
          end do
       end do
    end do

    ! Free values and close file
    call mem_dealloc(values)
    call array4_close_file(A,'keep')

    ! Set address counter back to 1
    A%address_counter = 1


  end subroutine array4_zeroone_file_type3


  !> \brief Set all address in array A to the "standard addresses"
  !> corresponding to storing arrays in "fortran order".
  !> See function array4_init_file.
  !> \author Kasper Kristensen
  !> \date December 2010
  subroutine array4_set_standard_address(A)
    implicit none
    !> Array
    type(array4), intent(inout) :: A

    StoringType: select case(A%storing_type)

    case(1)
       call array4_set_standard_address_type1(A)

    case(2)
       call array4_set_standard_address_type2(A)

    case(3)
       call array4_set_standard_address_type3(A)

    case default

       call lsquit('array4_set_standard_address: &
            & Requested storing type not implemented!', DECinfo%output)

    end select StoringType

  end subroutine array4_set_standard_address


  !> \brief Set all address in array A to the "standard addresses"
  !> corresponding to storing arrays in "fortran order" using storing type 1.
  !> See function array4_init_file.
  !> \author Kasper Kristensen
  !> \date December 2010
  subroutine array4_set_standard_address_type1(A)
    implicit none
    !> Full array information (but no memory allocated for values)
    type(array4), intent(inout) :: A
    integer :: i
    integer(kind=long) :: counter


    ! Sanity check
    ! ************
    ! Number of elements to write to disk must equal dims(1)*dims(2)*dims(3)
    if( A%nelements /= A%dims(1)*A%dims(2)*A%dims(3) ) then
       call lsquit('array4_set_standard_address_type1: &
            & A%nelements /= dim1*dim2*dim3', DECinfo%output)
    end if


    ! Initialize stuff
    ! ****************
    A%address_counter = 1
    counter=1


    ! Set address values
    ! ******************

    do i=1,A%dims(4)
       A%address(1,1,1,i) = counter
       counter = counter + A%nelements
    end do


  end subroutine array4_set_standard_address_type1



  !> \brief Set all address in array A to the "standard addresses"
  !> corresponding to storing arrays in "fortran order" using storing type 2.
  !> See function array4_init_file.
  !> \author Kasper Kristensen
  !> \date December 2010
  subroutine array4_set_standard_address_type2(A)
    implicit none
    !> Full array information (but no memory allocated for values)
    type(array4), intent(inout) :: A
    integer :: idx3,idx4
    integer(kind=long) :: counter

    ! Sanity check
    ! ************

    ! Number of elements to write to disk must equal dims(1)*dims(2)
    if( A%nelements /= A%dims(1)*A%dims(2) ) then
       call lsquit('array4_set_standard_address_type2: &
            & A%nelements /= dim1*dim2', DECinfo%output)
    end if


    ! Initialize stuff
    ! ****************
    A%address_counter = 1
    counter=1


    ! Set address values
    ! ******************

    do idx4=1,A%dims(4)
       do idx3=1,A%dims(3)
          A%address(1,1,idx3,idx4) = counter
          counter = counter + A%nelements
       end do
    end do


  end subroutine array4_set_standard_address_type2




  !> \brief Set all address in array A to the "standard addresses"
  !> corresponding to storing arrays in "fortran order" using storing type 3.
  !> See function array4_init_file.
  !> \author Kasper Kristensen
  !> \date December 2010
  subroutine array4_set_standard_address_type3(A)
    implicit none
    !> Full array information (but no memory allocated for values)
    type(array4), intent(inout) :: A
    integer :: idx2,idx3,idx4
    integer(kind=long) :: counter

    ! Sanity check
    ! ************

    ! Number of elements to write to disk must equal dims(1)
    if( A%nelements /= A%dims(1) ) then
       call lsquit('array4_set_standard_address_type2: &
            & A%nelements /= A%dims(1)', DECinfo%output)
    end if



    ! Initialize stuff
    ! ****************
    A%address_counter = 1
    counter = 1



    ! Set address values
    ! ******************

    do idx4=1,A%dims(4)
       do idx3=1,A%dims(3)
          do idx2=1,A%dims(2)
             A%address(1,idx2,idx3,idx4) = counter
             counter = counter + A%nelements
          end do
       end do
    end do


  end subroutine array4_set_standard_address_type3



  !> \brief Extract EOS indices from array4 and deletes buffer indices.
  !> NOTE: Assumes that the incoming array Arr has dimensions
  !> (nvirt,nocc,nvirt,nocc) and then extracts occupied EOS indices if
  !> DECinfo%VirtualScheme=.false. and
  !> extracts the virtual EOS indices otherwise.
  !> \author Kasper Kristensen
  !> \date December 2010
  subroutine array4_extract_eos_indices(Arr,MyFragment,occupied)


    implicit none
    !> Array where EOS indices are extracted
    type(array4),intent(inout) :: Arr
    !> Atomic fragment
    type(decfrag),intent(inout) :: MyFragment
    !> Occupied (true) or virtual (false) scheme
    logical,intent(in) :: occupied
    type(array4) :: tensor_copy
    integer :: nEOS,i
    integer, pointer :: EOS_idx(:)

    ! Copy of original array
    ! **********************

    if(DECinfo%array4OnFile) then
       ! Copy array and use same file as for original array
       tensor_copy = array4_duplicate_same_file(Arr)
    else
       ! Copy array elements stored in memory
       tensor_copy = array4_duplicate(Arr)
    end if

    ! Free original array - but save file if arrays are stored on file
    if(DECinfo%array4OnFile) then
       call array4_free(Arr,keep=.true.)
    else
       call array4_free(Arr)
    end if


    if(occupied) then ! occupied space partitioning

       ! Number of occupied EOS orbitals
       nEOS = MyFragment%noccEOS

       ! Occupied EOS orbital indices in the AOS (EOS+buffer) list of orbitals
       call mem_alloc(EOS_idx,nEOS)
       do i=1,nEOS
          EOS_idx(i) = MyFragment%idxo(i)
       end do

       if(DECinfo%array4OnFile) then ! array values stored on file
          call array4_extract_eos_indices_occ_file(Arr,tensor_copy,nEOS,EOS_idx)
       else ! array values are kept in memory
          call array4_extract_eos_indices_occ_memory(Arr,tensor_copy,nEOS,EOS_idx)
       end if

    else ! Virtual space partitioing

       ! Number of virtual EOS orbitals
       nEOS = MyFragment%nvirtEOS

       ! Virtual EOS orbital indices in the AOS (EOS+buffer) list of orbitals
       call mem_alloc(EOS_idx,nEOS)
       do i=1,nEOS
          EOS_idx(i) = MyFragment%idxu(i)
       end do

       if(DECinfo%array4OnFile) then ! array values stored on file
          call array4_extract_eos_indices_virt_file(Arr,tensor_copy,nEOS,EOS_idx)
       else ! array values are kept in memory
          call array4_extract_eos_indices_virt_memory(Arr,tensor_copy,nEOS,EOS_idx)
       end if

    end if


    ! Free stuff
    call mem_dealloc(EOS_idx)
    call array4_free(tensor_copy)

  end subroutine array4_extract_eos_indices




  !> \brief Extract EOS indices from array4 for both occupied and virtual partitioning schemes:
  !> 1. tensor_occEOS: The occupied orbitals not assigned to the central atom are removed while
  !>                the virtual indices are unchanged.
  !> 2. tensor_virtEOS: The virtual orbitals not assigned to the central atom are removed while
  !>                 the occupied indices are unchanged.
  !> \author Kasper Kristensen
  !> \date August 2011
  subroutine array4_extract_eos_indices_both_schemes(tensor_orig,tensor_occEOS,tensor_virtEOS,MyFragment)


    implicit none
    !> Array where occupied EOS indices are extracted
    type(array4),intent(inout) :: tensor_occEOS
    !> Array where virtual EOS indices are extracted
    type(array4),intent(inout) :: tensor_virtEOS
    !> Original array with AOS fragment indices for both occ and virt spaces
    type(array4),intent(in) :: tensor_orig
    !> Atomic fragment
    type(decfrag),intent(inout) :: MyFragment
    integer :: nocc, nvirt

    ! Number of occ and virt orbitals on central atom in fragment
    nocc = MyFragment%noccEOS
    nvirt = MyFragment%nvirtEOS


    ! Extract virtual EOS indices and leave occupied indices untouched
    ! ****************************************************************


    IF(.NOT.DECinfo%OnlyOccPart)THEN
       if(DECinfo%array4OnFile) then ! array values stored on file
          call array4_extract_eos_indices_virt_file(tensor_virtEOS,tensor_orig,&
               & nvirt,MyFragment%idxu(1:nvirt))
       else ! array values are kept in memory
          call array4_extract_eos_indices_virt_memory(tensor_virtEOS,tensor_orig,&
               & nvirt,MyFragment%idxu(1:nvirt))
       end if
    ENDIF

    ! Extract occupied EOS indices and leave virtual indices untouched
    ! ****************************************************************

    IF(.NOT.DECinfo%OnlyVIRTPart)THEN
       if(DECinfo%array4OnFile) then ! array values stored on file
          call array4_extract_eos_indices_occ_file(tensor_occEOS,tensor_orig,&
               & nocc, MyFragment%idxo(1:nocc))
       else ! array values are kept in memory
          call array4_extract_eos_indices_occ_memory(tensor_occEOS,tensor_orig,&
               & nocc, MyFragment%idxo(1:nocc))
       end if
    ENDIF

  end subroutine array4_extract_eos_indices_both_schemes



  !> \brief Extract virtual EOS indices from array4 and deletes buffer indices.
  !> NOTE: Assumes that the incoming array Arr has dimensions (nvirt,nocc,nvirt,nocc)
  !> and that the array values are kept in memory.
  !> \author Kasper Kristensen
  !> \date December 2010
  subroutine array4_extract_eos_indices_virt_memory(Arr,tensor_copy,nEOS,EOS_idx)


    implicit none
    !> Array output where EOS indices are extracted
    type(array4),intent(inout) :: Arr
    !> Original array
    type(array4),intent(in) :: tensor_copy
    !> Number of EOS indices
    integer,intent(in) :: nEOS
    !> List of EOS indices in the total (EOS+buffer) list of orbitals
    integer, dimension(nEOS),intent(in) :: EOS_idx
    integer :: nocc,nvirt,i,a,b,j,ax,bx
    integer, dimension(4) :: new_dims



    ! Initialize stuff
    ! ****************
    nocc = tensor_copy%dims(2)  ! Total number of occupied orbitals
    nvirt = tensor_copy%dims(1)  ! Total number of virtual orbitals
    new_dims=[nEOS,nocc,nEOS,nocc] ! nEOS=Number of virtual EOS orbitals



    ! Sanity checks
    ! *************

    ! 1. Positive number of orbitals
    if( (nocc<1) .or. (nvirt<1) ) then
       write(DECinfo%output,*) 'nocc = ', nocc
       write(DECinfo%output,*) 'nvirt = ', nvirt
       call lsquit('array4_extract_eos_indices_virt_memory: &
            & Negative or zero number of orbitals!',DECinfo%output)
    end if

    ! 2. Array structure is (virt,occ,virt,occ)
    if( (nvirt/=tensor_copy%dims(3)) .or. (nocc/=tensor_copy%dims(4)) ) then
       write(DECinfo%output,*) 'tensor_copy%dims(1) = ', tensor_copy%dims(1)
       write(DECinfo%output,*) 'tensor_copy%dims(2) = ', tensor_copy%dims(2)
       write(DECinfo%output,*) 'tensor_copy%dims(3) = ', tensor_copy%dims(3)
       write(DECinfo%output,*) 'tensor_copy%dims(4) = ', tensor_copy%dims(4)
       call lsquit('array4_extract_eos_indices_virt_memory: &
            & Arr dimensions does not match (virt,occ,virt,occ) structure!',DECinfo%output)
    end if

    ! 3. EOS dimension must be smaller than (or equal to) total number of virt orbitals
    if(nEOS > nvirt) then
       write(DECinfo%output,*) 'nvirt = ', nvirt
       write(DECinfo%output,*) 'nEOS  = ', nEOS
       call lsquit('array4_extract_eos_indices_virt_memory: &
            & Number of EOS orbitals must be smaller than (or equal to) total number of &
            & virtual orbitals!',DECinfo%output)
    end if

    ! 4. EOS indices must not exceed total number of virtual orbitals
    do i=1,nEOS
       if(EOS_idx(i) > nvirt) then
          write(DECinfo%output,'(a,i6,a)') 'EOS index number ', i, ' is larger than nvirt!'
          write(DECinfo%output,*) 'nvirt   = ', nvirt
          write(DECinfo%output,*) 'EOS_idx = ', EOS_idx(i)
          call lsquit('array4_extract_eos_indices_virt_memory: &
               & EOS index value larger than nvirt!',DECinfo%output)
       end if
    end do


    ! Extract virtual EOS indices and store in Arr
    ! ********************************************

    ! Initiate Arr with new dimensions (nvirt_EOS,nocc,nvirt_EOS,nocc)
    Arr=array4_init(new_dims)

    ! Set Arr equal to the EOS indices of the original Arr array (tensor_copy)
    do j=1,nocc
       do b=1,nEOS
          bx=EOS_idx(b)
          do i=1,nocc
             do a=1,nEOS
                ax=EOS_idx(a)
                Arr%val(a,i,b,j) = tensor_copy%val(ax,i,bx,j)
             end do
          end do
       end do
    end do


  end subroutine array4_extract_eos_indices_virt_memory


  !> \brief Extract virtual EOS indices from array4 and deletes buffer indices.
  !> NOTE: Assumes that the incoming array Arr has dimensions (nvirt,nocc,nvirt,nocc)
  !> and that the array values are stored on file.
  !> \author Kasper Kristensen
  !> \date December 2010
  subroutine array4_extract_eos_indices_virt_file(Arr,tensor_copy,nEOS,EOS_idx)


    implicit none
    !> Array output where EOS indices are extracted
    type(array4),intent(inout) :: Arr
    !> Original array
    type(array4),intent(in) :: tensor_copy
    !> Number of EOS indices
    integer,intent(in) :: nEOS
    !> List of EOS indices in the total (EOS+buffer) list of orbitals
    integer, dimension(nEOS),intent(in) :: EOS_idx
    integer :: nocc,nvirt,i,j,ax,bx,a,b
    integer, dimension(4) :: new_dims
    real(realk),pointer :: values(:,:,:)


    ! Initialize stuff
    ! ****************
    nocc = tensor_copy%dims(2)  ! Total number of occupied orbitals
    nvirt = tensor_copy%dims(1)  ! Total number of virtual orbitals
    new_dims=[nEOS,nocc,nEOS,nocc] ! nEOS=Number of virtual EOS orbitals


    ! Sanity checks
    ! *************

    ! 1. Positive number of orbitals
    if( (nocc<1) .or. (nvirt<1) ) then
       write(DECinfo%output,*) 'nocc = ', nocc
       write(DECinfo%output,*) 'nvirt = ', nvirt
       call lsquit('array4_extract_eos_indices_virt_file: &
            & Negative or zero number of orbitals!',DECinfo%output)
    end if

    ! 2. Array structure is (virt,occ,virt,occ)
    if( (nvirt/=tensor_copy%dims(3)) .or. (nocc/=tensor_copy%dims(4)) ) then
       write(DECinfo%output,*) 'tensor_copy%dims(1) = ', tensor_copy%dims(1)
       write(DECinfo%output,*) 'tensor_copy%dims(2) = ', tensor_copy%dims(2)
       write(DECinfo%output,*) 'tensor_copy%dims(3) = ', tensor_copy%dims(3)
       write(DECinfo%output,*) 'tensor_copy%dims(4) = ', tensor_copy%dims(4)
       call lsquit('array4_extract_eos_indices_virt_file: &
            & Arr dimensions does not match (virt,occ,virt,occ) structure!',DECinfo%output)
    end if

    ! 3. EOS dimension must be smaller than (or equal to) total number of virt orbitals
    if(nEOS > nvirt) then
       write(DECinfo%output,*) 'nvirt = ', nvirt
       write(DECinfo%output,*) 'nEOS = ', nEOS
       call lsquit('array4_extract_eos_indices_virt_file: &
            & Number of EOS orbitals must be smaller than (or equal to) total number of &
            & virtual orbitals!',DECinfo%output)
    end if

    ! 4. EOS indices must not exceed total number of virtual orbitals
    do i=1,nEOS
       if(EOS_idx(i) > nvirt) then
          write(DECinfo%output,'(a,i6,a)') 'EOS index number ', i, ' is larger than nvirt!'
          write(DECinfo%output,*) 'nvirt = ', nvirt
          write(DECinfo%output,*) 'EOS_idx = ', EOS_idx(i)
          call lsquit('array4_extract_eos_indices_virt_file: &
               & EOS index value larger than nocc!',DECinfo%output)
       end if
    end do

    ! 5. Currently only implemented for storing type 1
    if(tensor_copy%storing_type /= 1) then
       write(DECinfo%output,*) 'Array storing type: ', tensor_copy%storing_type
       call lsquit('array4_extract_eos_indices_virt_file: &
            & Only implemented for storing type 1!',DECinfo%output)
    end if


    ! Extract virtual EOS indices and store in Arr
    ! ********************************************

    ! tensor_copy contains all necessary information, and we can therefore free Arr
    ! and re-initialize. This is done manually not to screw up array bookkeeping.


    ! Re-initialize Arr manually (as in array4_init_standard).
    ! (Do not use array4_init_standard, that would destroy the array bookkeeping).
    Arr%dims=new_dims
    Arr%order=[1,2,3,4]
    call memory_allocate_4d(Arr%val,Arr%dims)
    Arr%funit=0
    Arr%filename = 'NoFilename'
    Arr%address_counter=0
    Arr%storing_type=0
    Arr%nelements=0

    ! Temporary array
    call mem_alloc(values,nvirt,nocc,nvirt)

    ! Open array file
    call array4_open_file(tensor_copy)


    do j=1,nocc

       ! Read in (a,i,b,j) values for all AIB for given j
       call array4_read_file_type1(tensor_copy,j,values,nvirt,nocc,nvirt)

       do b=1,nEOS
          bx=EOS_idx(b) ! Only consider occupied EOS indices
          do i=1,nocc
             do a=1,nEOS
                ax=EOS_idx(a)
                Arr%val(a,i,b,j) = values(ax,i,bx)
             end do
          end do
       end do
    end do


    call array4_close_file(tensor_copy,'KEEP')
    call mem_dealloc(values)



  end subroutine array4_extract_eos_indices_virt_file





  !> \brief Extract occupied EOS indices from array4 and deletes buffer indices.
  !> NOTE: Assumes that the incoming array Arr has dimensions (nvirt,nocc,nvirt,nocc)
  !> and that the array values are kept in memory.
  !> \author Kasper Kristensen
  !> \date December 2010
  subroutine array4_extract_eos_indices_occ_memory(Arr,tensor_copy,nEOS,EOS_idx)


    implicit none
    !> Array where EOS indices where are extracted
    type(array4),intent(inout) :: Arr
    !> Original array
    type(array4),intent(in) :: tensor_copy
    !> Number of EOS indices
    integer,intent(in) :: nEOS
    !> List of EOS indices in the total (EOS+buffer) list of orbitals
    integer, dimension(nEOS),intent(in) :: EOS_idx
    integer :: nocc,nvirt,i,a,b,j,ix,jx
    integer, dimension(4) :: new_dims


    ! Initialize stuff
    ! ****************
    nocc = tensor_copy%dims(2)  ! Total number of occupied orbitals
    nvirt = tensor_copy%dims(1)  ! Total number of virtual orbitals
    new_dims=[nvirt,nEOS,nvirt,nEOS] ! nEOS=Number of occupied EOS orbitals



    ! Sanity checks
    ! *************

    ! 1. Positive number of orbitals
    if( (nocc<1) .or. (nvirt<1) ) then
       write(DECinfo%output,*) 'nocc = ', nocc
       write(DECinfo%output,*) 'nvirt = ', nvirt
       call lsquit('array4_extract_eos_indices_occ_memory: &
            & Negative or zero number of orbitals!',DECinfo%output)
    end if

    ! 2. Array structure is (virt,occ,virt,occ)
    if( (nvirt/=tensor_copy%dims(3)) .or. (nocc/=tensor_copy%dims(4)) ) then
       write(DECinfo%output,*) 'tensor_copy%dims(1) = ', tensor_copy%dims(1)
       write(DECinfo%output,*) 'tensor_copy%dims(2) = ', tensor_copy%dims(2)
       write(DECinfo%output,*) 'tensor_copy%dims(3) = ', tensor_copy%dims(3)
       write(DECinfo%output,*) 'tensor_copy%dims(4) = ', tensor_copy%dims(4)
       call lsquit('array4_extract_eos_indices_occ_memory: &
            & Arr dimensions does not match (virt,occ,virt,occ) structure!',DECinfo%output)
    end if

    ! 3. EOS dimension must be smaller than (or equal to) total number of occ orbitals
    if(nEOS > nocc) then
       write(DECinfo%output,*) 'nocc = ', nocc
       write(DECinfo%output,*) 'nEOS = ', nEOS
       call lsquit('array4_extract_eos_indices_occ_memory: &
            & Number of EOS orbitals must be smaller than (or equal to) total number of &
            & occupied orbitals!',DECinfo%output)
    end if

    ! 4. EOS indices must not exceed total number of occupied orbitals
    do i=1,nEOS
       if(EOS_idx(i) > nocc) then
          write(DECinfo%output,'(a,i6,a)') 'EOS index number ', i, ' is larger than nocc!'
          write(DECinfo%output,*) 'nocc = ', nocc
          write(DECinfo%output,*) 'EOS_idx = ', EOS_idx(i)
          call lsquit('array4_extract_eos_indices_occ_memory: &
               & EOS index value larger than nocc!',DECinfo%output)
       end if
    end do


    ! Extract occupied EOS indices and store in Arr
    ! *********************************************


    ! Initiate Arr with new dimensions (nvirt,nocc_EOS,nvirt,nocc_EOS)
    Arr=array4_init(new_dims)

    ! Set Arr equal to the EOS indices of the original Arr array (tensor_copy)
    do j=1,nEOS
       jx=EOS_idx(j)
       do b=1,nvirt
          do i=1,nEOS
             ix=EOS_idx(i)
             do a=1,nvirt
                Arr%val(a,i,b,j) = tensor_copy%val(a,ix,b,jx)
             end do
          end do
       end do
    end do


  end subroutine array4_extract_eos_indices_occ_memory


  !> \brief Extract occupied EOS indices from array4 and deletes buffer indices.
  !> NOTE: Assumes that the incoming array Arr has dimensions (nvirt,nocc,nvirt,nocc)
  !> and that the array values are stored on file in the INPUT array.
  !> The OUTPUT array values are stored in memory.
  !> \author Kasper Kristensen
  !> \date December 2010
  subroutine array4_extract_eos_indices_occ_file(Arr,tensor_copy,nEOS,EOS_idx)


    implicit none
    !> Array where EOS indices are extracted
    type(array4),intent(inout) :: Arr
    !> Original array
    type(array4),intent(in) :: tensor_copy
    !> Number of EOS indices
    integer,intent(in) :: nEOS
    !> List of EOS indices in the total (EOS+buffer) list of orbitals
    integer, dimension(nEOS),intent(in) :: EOS_idx
    integer :: nocc,nvirt,i,a,b,j,ix,jx
    integer, dimension(4) :: new_dims
    real(realk),pointer :: values(:,:,:)


    ! Initialize stuff
    ! ****************
    nocc = tensor_copy%dims(2)  ! Total number of occupied orbitals
    nvirt = tensor_copy%dims(1)  ! Total number of virtual orbitals
    new_dims=[nvirt,nEOS,nvirt,nEOS] ! nEOS=Number of occupied EOS orbitals


    ! Sanity checks
    ! *************

    ! 1. Positive number of orbitals
    if( (nocc<1) .or. (nvirt<1) ) then
       write(DECinfo%output,*) 'nocc = ', nocc
       write(DECinfo%output,*) 'nvirt = ', nvirt
       call lsquit('array4_extract_eos_indices_occ_file: &
            & Negative or zero number of orbitals!',DECinfo%output)
    end if

    ! 2. Array structure is (virt,occ,virt,occ)
    if( (nvirt/=tensor_copy%dims(3)) .or. (nocc/=tensor_copy%dims(4)) ) then
       write(DECinfo%output,*) 'tensor_copy%dims(1) = ', tensor_copy%dims(1)
       write(DECinfo%output,*) 'tensor_copy%dims(2) = ', tensor_copy%dims(2)
       write(DECinfo%output,*) 'tensor_copy%dims(3) = ', tensor_copy%dims(3)
       write(DECinfo%output,*) 'tensor_copy%dims(4) = ', tensor_copy%dims(4)
       call lsquit('array4_extract_eos_indices_occ_file: &
            & Arr dimensions does not match (virt,occ,virt,occ) structure!',DECinfo%output)
    end if

    ! 3. EOS dimension must be smaller than (or equal to) total number of occ orbitals
    if(nEOS > nocc) then
       write(DECinfo%output,*) 'nocc = ', nocc
       write(DECinfo%output,*) 'nEOS = ', nEOS
       call lsquit('array4_extract_eos_indices_occ_file: &
            & Number of EOS orbitals must be smaller than (or equal to) total number of &
            & occupied orbitals!',DECinfo%output)
    end if

    ! 4. EOS indices must not exceed total number of occupied orbitals
    do i=1,nEOS
       if(EOS_idx(i) > nocc) then
          write(DECinfo%output,'(a,i6,a)') 'EOS index number ', i, ' is larger than nocc!'
          write(DECinfo%output,*) 'nocc = ', nocc
          write(DECinfo%output,*) 'EOS_idx = ', EOS_idx(i)
          call lsquit('array4_extract_eos_indices_occ_file: &
               & EOS index value larger than nocc!',DECinfo%output)
       end if
    end do

    ! 5. Currently only implemented for storing type 1
    if(tensor_copy%storing_type /= 1) then
       write(DECinfo%output,*) 'Array storing type: ', tensor_copy%storing_type
       call lsquit('array4_extract_eos_indices_occ_file: &
            & Only implemented for storing type 1!',DECinfo%output)
    end if


    ! Extract occupied EOS indices and store in Arr
    ! *********************************************


    ! tensor_copy now contains all necessary information, and we can therefore free Arr
    ! and re-initialize. This is done manually not to screw up array bookkeeping.


    ! Re-initialize Arr manually (as in array4_init_standard).
    ! (Do not use array4_init_standard, that would destroy the array bookkeeping).
    Arr%dims=new_dims
    Arr%order=[1,2,3,4]
    call memory_allocate_4d(Arr%val,Arr%dims)
    Arr%funit=0
    Arr%filename = 'NoFilename'
    Arr%address_counter=0
    Arr%storing_type=0
    Arr%nelements=0

    ! Temporary array
    call mem_alloc(values,nvirt,nocc,nvirt)

    ! Open array file
    call array4_open_file(tensor_copy)


    do j=1,nEOS
       jx=EOS_idx(j)

       ! Read in (A,I,B,JX) values for all AIB for given EOS index JX
       call array4_read_file_type1(tensor_copy,jx,values,nvirt,nocc,nvirt)

       do b=1,nvirt
          do i=1,nEOS
             ix=EOS_idx(i)  ! Only consider occupied EOS indices
             do a=1,nvirt
                Arr%val(a,i,b,j) = values(a,ix,b)
             end do
          end do
       end do
    end do


    call array4_close_file(tensor_copy,'KEEP')
    call mem_dealloc(values)


  end subroutine array4_extract_eos_indices_occ_file



  !> \brief Dump array values to file and free array memory.
  !> Note 1: Filename in array structure is set equal to the filename input.
  !> Note 2: Sets the array storing type to 1!
  !> \author Kasper Kristensen
  !> \date September 2011
  subroutine dump_array4_to_file(Arr,filename)
    implicit none
    !> Array to write to file
    type(array4),intent(inout) :: Arr
    !> Filename where array is written from
    character(*), intent(in) :: filename
    integer :: i

    ! Set array info to write to file "filename" using storing type 1 (see array4_init_file)
    arr%storing_type=1
    if(associated(arr%address)) then
       call mem_dealloc(arr%address)
       nullify(arr%address)
    endif
    call mem_alloc(arr%address,1,1,1,arr%dims(4) )
    arr%nelements = arr%dims(1)*arr%dims(2)*arr%dims(3)
    arr%filename=filename
    call array4_set_standard_address(arr)
    ! Just in case, delete file "filename", if is already exists
    call array4_delete_file(arr)

    ! Sanity check
    if(.not. associated(Arr%val)) then
       call lsquit('dump_array4_to_file: Array values are not associated!',DECinfo%output)
    end if

    ! Open array file
    call array4_open_file(Arr)

    ! Write array to file in chuncks of Arr%val(:,:,:,i) for a fixed "i".
    do i=1,Arr%dims(4)
       call array4_write_file_type1(Arr,i, Arr%val(1:Arr%dims(1), 1:Arr%dims(2), 1:Arr%dims(3),i), &
            & Arr%dims(1),Arr%dims(2),Arr%dims(3) )
    end do

    ! Close array file
    call array4_close_file(Arr,'KEEP')

    ! Deallocate array4 pointer
    call memory_deallocate_4d(Arr%val)


  end subroutine dump_array4_to_file



  !> \brief Retrive array values from file and reallocate array memory.
  !> Note 1: Reads array values from the filename associated with the array4 structure.
  !> Note 2: Assumes that array storing type is 1!
  !> \author Kasper Kristensen
  !> \date September 2011
  subroutine retrieve_array4_from_file(Arr,keep_or_delete)
    implicit none
    !> Array to write to file
    type(array4),intent(inout) :: Arr
    !> Status for array file after reading - 'KEEP' og 'DELETE'
    character(*), intent(in) :: keep_or_delete
    integer :: i

    ! Sanity checks
    if(Arr%storing_type /= 1) then
       write(DECinfo%output,*) 'Array storing type = ', Arr%storing_type
       call lsquit('retrieve_array4_from_file: Only implemented for storing type 1!',DECinfo%output)
    end if

    if(associated(Arr%val)) then
       call lsquit('retrive_array4_from_file: Array values are should not associated!',DECinfo%output)
    end if


    ! Allocate array memory
    call memory_allocate_4d(Arr%val,Arr%dims)


    ! Open array file
    call array4_open_file(Arr)

    ! Read array values from file in chuncks of Arr%val(:,:,:,i) for a fixed "i".
    do i=1,Arr%dims(4)
       call array4_read_file_type1(Arr,i, Arr%val(1:Arr%dims(1), 1:Arr%dims(2), 1:Arr%dims(3),i), &
            & Arr%dims(1),Arr%dims(2),Arr%dims(3) )
    end do


    ! Close array file
    call array4_close_file(Arr,keep_or_delete)


  end subroutine retrieve_array4_from_file


  !> \brief Print all elements of four-dimensional array to LSDALTON.OUT.
  !> Only to be used for testing purposes!
  !> \author Kasper Kristensen
  !> \date November 2011
  subroutine array4_print(A,label)
    implicit none

    ! Array to be printed
    type(array4) :: A
    !> Label for array
    character(*), intent(in) :: label
    integer,dimension(4) :: d

    d=A%dims
    call print_4_dimensional_array( d, A%val(1:d(1),1:d(2),1:d(3),1:d(4)), label )

  end subroutine array4_print


  !> \brief Transform virtual orbital indices in local basis (or whatever the input basis is) 
  !> for doubles amplitudes to fragment-adapted orbital basis, leaving occupied indices untouched.
  !> \author Kasper Kristensen
  !> \date September 2013
  subroutine transform_virt_amp_to_FOs(t2,MyFragment)
    implicit none
    !> Doubles amplitudes stored as (a,i,b,j)
    type(array4),intent(inout) :: t2
    !> Atomic (or pair) fragment; the virtual correlation density which implicitly defines 
    !> the fragment-adapted orbital basis is stored in MyFragment%VirtMat
    type(decfrag),intent(in) :: MyFragment
    real(realk),pointer :: U(:,:),eival(:),Ured(:,:),tmp(:,:,:,:),tmp2(:,:,:,:)
    integer :: nvirt,nvirtFO,i,idx,nocc,dims(4)
    real(realk) :: thr

    ! Sanity check, this will mix all virtual orbitals, so it can only be used 
    ! for the occupied partitioning scheme.
    if(.not. DECinfo%OnlyOccPart) then
       call lsquit('transform_virt_amp_to_FOs: Will mix all virtual orbitals, so it can &
            & only be used for the occupied partitioning scheme! Suggestion: Turn on .ONLYOCCPART.',-1)
    end if
    nvirt = t2%dims(1)
    nocc = t2%dims(2)
    if(nvirt/=MyFragment%nvirtAOS .or. nocc/=MyFragment%noccAOS) then
       print *, 't2%dims', t2%dims
       print *, 'nvirt,nocc', MyFragment%nvirtAOS,MyFragment%noccAOS
       call lsquit('transform_virt_amp_to_FOs: Dimension mismatch for amplitudes and&
            & associated fragment',-1)
    end if


    ! Determine transformation matrix from local basis to FO basis
    ! ************************************************************
    call mem_alloc(U,nvirt,nvirt)
    call mem_alloc(eival,nvirt)

    ! Diagonalize virtual correlation density matrix
    call solve_eigenvalue_problem_unitoverlap(nvirt,MyFragment%VirtMat,eival,U)

    ! Throw away eigenvalues smaller than rejection threshold
    thr = MyFragment%RejectThr(2)
    
    ! Number of eigenvalues larger than threshold
    nvirtFO=0
    do i=1,nvirt
       if(eival(i)>thr) then
          nvirtFO = nvirtFO+1
       end if
    end do

    ! Sanity check
    if(nvirtFO==0) then
       write(DECinfo%output,*) 'WARNING! Number of virtual FOs is zero for thr= ', thr
       write(DECinfo%output,*) '--> to avoid artificial behaviour I set nvirtFO=nvirt!'
       nvirtFO = nvirt
    end if
    
    ! Transformation matrix for (local,FO) indices:
    call mem_alloc(Ured,nvirt,nvirtFO)
    idx=0
    do i=1,nvirt
       if(eival(i)>thr) then
          idx=idx+1
          Ured(:,idx) = U(:,i)
       end if
    end do
    call mem_dealloc(eival)
    call mem_dealloc(U)


    ! Transform the two virtual indices
    ! *********************************
    ! We denote the original (probably local) indices by a,b,i,j and the FO indices by A,B.
    !
    ! We want to make the transformation: t2(a,i,b,j) --> t2(A,i,B,j)

    ! tmp(A,i,b,j) = sum_a Ured^T(A,a) t2(a,i,b,j)
    call mem_alloc(tmp,nvirtFO,nocc,nvirt,nocc)
    call dgemm('t','n',nvirtFO,nocc*nvirt*nocc,nvirt,1.0E0_realk,Ured,nvirt,t2%val,nvirt,&
         & 0.0E0_realk,tmp,nvirtFO)
    call array4_free(t2)

    ! Reorder: tmp(A,i,b,j) --> t2(b,j,A,i)
    call mem_alloc(tmp2,nvirt,nocc,nvirtFO,nocc)
    call mat_transpose(nvirtFO*nocc,nvirt*nocc,1.0E0_realk,tmp,0.0E0_realk,tmp2)
    call mem_dealloc(tmp)

    ! t2(B,j,A,i) = sum_b Ured^T(B,b) tmp2(b,j,A,i)
    dims(1) = nvirtFO
    dims(2) = nocc
    dims(3) = nvirtFO
    dims(4) = nocc
    t2 = array4_init(dims)
    call dgemm('t','n',nvirtFO,nocc*nvirtFO*nocc,nvirt,1.0E0_realk,Ured,nvirt,tmp2,nvirt,&
         & 0.0E0_realk,t2%val,nvirtFO)
    call mem_dealloc(tmp2)

    ! Since t2(B,j,A,i) = t2(A,i,B,j), the t2 output array now contains the desired amplitudes.

    call mem_dealloc(Ured) 

  end subroutine transform_virt_amp_to_FOs

end module array4_simple_operations
