!dalton_copyright_start
!
!
!dalton_copyright_end
!
! stefan oct 2011: this module contains the type definition for the
!                  MPI-file-list/handle used in the MCSCF-LUCITA interface.
!
module file_type_module

  implicit none

  public file_type
  public file_init_lucipar
  public file_free_lucipar
  public file_set_list_lucipar
  public file_type_active_fh_reset_lucipar


! type definition
  type file_type

    integer ::                  &
      active_nr_f_lucipar,      &             ! current active number of MPI files in lucita
      current_file_fh_seqf(2),  &             ! present non-parallel file handle for in/outgoing vector(s) #1 and #2
      current_file_nr_active1,  &             ! 
      current_file_nr_active2,  &             ! 
      current_file_nr_bvec,     &      
      current_file_nr_diag,     &      
      max_list_length_fac,      &             ! 
      max_list_length,          &             ! 
      max_list_length_bvec                    !

    logical ::                  &
      file_type_init       = .false.          ! status of type file_type
    logical ::                  &
      file_handle_seq_init = .false.          ! status of initialization of sequential file handle declaration
    logical ::                  &
      file_type_mc         = .false.          ! set file types for MCSCF

    integer, allocatable ::     &
      iluxlist(:,:),            &             ! file entry list for all active parallel MPI files
      ilublist(:,:),            &             ! file entry list for all active parallel MPI files
      facofffl(:),              &             ! file offset factor list for all active parallel MPI files
      fh_lu(:)                                ! file handle list

    integer(8), allocatable ::  &
      file_offsets(:)

  end type file_type

! ttss block type object
  type(file_type),    public, save  :: file_info
  integer, parameter, public  :: mx_nr_files_lucipar    = 10 ! max number of MPI files in lucita (general upper limit)
  integer, parameter, private :: mx_nr_files_lucipar_ci = 10 ! max number of MPI files in lucita/ci
  integer, parameter, private :: mx_nr_files_lucipar_mc =  5 ! max number of MPI files in lucita/mcscf

contains 

  subroutine file_init_lucipar(A, mx_vector, nr_eigenstates)

!   ----------------------------------------------------------------------------
    type(file_type)      :: A
    integer, intent(in)  :: mx_vector, nr_eigenstates
!   ----------------------------------------------------------------------------

!   reset old type information
    call file_free_lucipar(A)

    A%file_type_init          = .true.
    A%active_nr_f_lucipar     = -1
    A%current_file_nr_active1 = -1
    A%current_file_nr_active2 = -1
    A%current_file_nr_bvec    = -1
    A%current_file_nr_diag    = -1
    A%max_list_length_fac     = -1
    A%max_list_length         = -1
    A%max_list_length_bvec    = -1

    if(A%file_type_mc)then
      A%active_nr_f_lucipar = mx_nr_files_lucipar_mc
    else
      A%active_nr_f_lucipar = mx_nr_files_lucipar_ci
    end if

    allocate(A%facofffl(A%active_nr_f_lucipar))
    allocate(A%fh_lu(A%active_nr_f_lucipar))
    allocate(A%file_offsets(A%active_nr_f_lucipar))

    call file_set_factors_lucipar(A, mx_vector, nr_eigenstates)

  end subroutine file_init_lucipar

  subroutine file_set_list_lucipar(A, number_of_active_ttss_blocks_per_process, number_of_ttss_blocks)

!   ----------------------------------------------------------------------------
    type(file_type)      :: A
    integer, intent(in)  :: number_of_active_ttss_blocks_per_process
    integer, intent(in)  :: number_of_ttss_blocks
!   ----------------------------------------------------------------------------

    A%max_list_length      = number_of_active_ttss_blocks_per_process * A%max_list_length_fac
    A%max_list_length_bvec = number_of_ttss_blocks

    allocate(A%iluxlist(A%max_list_length,A%active_nr_f_lucipar-1))
    allocate(A%ilublist(A%max_list_length_bvec,1))

  end subroutine file_set_list_lucipar

  subroutine file_free_lucipar(A)

!   ----------------------------------------------------------------------------
    type(file_type) :: A
!   ----------------------------------------------------------------------------

    if(.not. A%file_type_init) return

    A%file_type_init          = .false.
    A%file_type_mc            = .false.
    A%active_nr_f_lucipar     = -1

    if(allocated(A%iluxlist))     deallocate(A%iluxlist)
    if(allocated(A%ilublist))     deallocate(A%ilublist)
    if(allocated(A%facofffl))     deallocate(A%facofffl)
    if(allocated(A%fh_lu))        deallocate(A%fh_lu)
    if(allocated(A%file_offsets)) deallocate(A%file_offsets)

  end subroutine file_free_lucipar

  subroutine file_set_factors_lucipar(A, mx_vector, nr_eigenstates)

!   ----------------------------------------------------------------------------
    type(file_type)      :: A
    integer, intent(in)  :: mx_vector
    integer, intent(in)  :: nr_eigenstates
!   ----------------------------------------------------------------------------

    if(A%file_type_mc)then
      A%facofffl( 1)        = 1
      A%facofffl( 2)        = 1
      A%facofffl( 3)        = 1
      A%facofffl( 4)        = 1
      A%facofffl( 5)        = 0

      A%max_list_length_fac = 1
    else
      A%facofffl( 1)        = 1
      A%facofffl( 2)        = nr_eigenstates * 4
      A%facofffl( 3)        = nr_eigenstates + mx_vector
      A%facofffl( 4)        = 1
      A%facofffl( 5)        = 0
      A%facofffl( 6)        = nr_eigenstates + mx_vector
      A%facofffl( 7)        = nr_eigenstates + mx_vector
      A%facofffl( 8)        = nr_eigenstates + mx_vector
      A%facofffl( 9)        = 1
      A%facofffl(10)        = mx_vector
      
      A%max_list_length_fac = nr_eigenstates + mx_vector
    end if

    A%current_file_nr_diag  = 4
    A%current_file_nr_bvec  = 5

  end subroutine file_set_factors_lucipar

  subroutine file_type_active_fh_reset_lucipar(A)

!   ----------------------------------------------------------------------------
    type(file_type) :: A
!   ----------------------------------------------------------------------------

    if(A%file_handle_seq_init) return
    A%current_file_fh_seqf(1) = -1
    A%current_file_fh_seqf(2) = -1
    A%current_file_nr_active1 = -1
    A%current_file_nr_active2 = -1

  end subroutine file_type_active_fh_reset_lucipar
 
end module file_type_module
