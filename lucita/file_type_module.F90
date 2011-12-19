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


! type definition
  type file_type_lucipar

    integer ::                  &
      active_nr_f_lucipar!      &             ! current active number of MPI files in lucita

    logical ::                  &
      file_type_init = .false.                ! status of type file_type
    logical ::                  &
      file_type_mc   = .false.                ! set file types for MCSCF

    integer, allocatable ::     &
      iluxlist(:,:),            &             ! file entry list for all active parallel MPI files
      facofffl(:),              &             ! file offset factor list for all active parallel MPI files
      fh_lu(:)                                ! file handle list

    integer(8), allocatable ::  &
      file_offsets(:)

  end type file_type_lucipar

! ttss block type object
  type(file_type), public     :: file_info
  integer, parameter, private :: mx_nr_files_lucipar    = 10 ! max number of MPI files in lucita (general upper limit)
  integer, parameter, private :: mx_nr_files_lucipar_ci = 10 ! max number of MPI files in lucita/ci
  integer, parameter, private :: mx_nr_files_lucipar_mc =  5 ! max number of MPI files in lucita/mcscf

contains 

  subroutine file_init_lucipar(A, mx_ttss, nirreps, sum_of_ci_spaces)

!   ----------------------------------------------------------------------------
    type(file_type)      :: A
    integer, intent(in)  :: mx_files, nirreps, sum_of_ci_spaces
!   ----------------------------------------------------------------------------

!   reset old type information
    call file_free_lucipar(A)

    nullify(A%iluxlist)
    nullify(A%facofffl)
    nullify(A%fh_lu)

    A%file_type_init      = .true.
    A%mx_nr_files_lucipar = mx_ttss
    A%active_nr_f_lucipar = nirreps
    A%nr_ci_spaces        = sum_of_ci_spaces
    A%present_sym_irrep   = -1
    A%present_ci_space    = -1
    A%total_present_ttss  = -1 
    A%total_present_vec   = -1 

    allocate(A%ttss_block_length(A%mx_nr_ttss,A%nr_ptg_irreps,A%nr_ci_spaces))
    allocate(A%ttss_block_nr(A%nr_ptg_irreps,A%nr_ci_spaces))
    allocate(A%ttss_vec_length(A%nr_ptg_irreps,A%nr_ci_spaces))

  end subroutine ttss_init

  subroutine ttss_free(A)

!   ----------------------------------------------------------------------------
    type(file_type) :: A
!   ----------------------------------------------------------------------------

    if(.not. A%ttss_info_init) return

    A%ttss_info_init = .false.
    deallocate(A%ttss_vec_length)
    deallocate(A%ttss_block_length)
    deallocate(A%ttss_block_nr)
    nullify(A%ttss_vec_length)
    nullify(A%ttss_block_length)
    nullify(A%ttss_block_nr)

  end subroutine ttss_free

end module file_type_module
