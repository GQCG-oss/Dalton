!dalton_copyright_start
!
!
!dalton_copyright_end
!
! stefan oct 2011: this module contains the type definition for the ttss_block_structure 
!                  used in the MCSCF-LUCITA interface.
!
module ttss_block_module

  implicit none

  public ttss_block_structure
  public ttss_init
  public ttss_free
  public ttss_switch


! type definition
  type ttss_block_structure

    integer ::                  &
      mx_nr_ttss,               &             ! max number of ttss blocks in all symmetry irreps
      nr_ptg_irreps,            &             ! number of point group irreps
      nr_ci_spaces,             &             ! number of CI spaces (currently always 1)
      present_sym_irrep,        &             ! current active symmetry irrep in use
      present_ci_space,         &             ! current active CI space in use (currently always 1, see above)
      total_present_vec,        &             ! vector length for the present_sym_irrep and active CI space in use (the latter currently always 1)
      total_present_ttss                      ! current total number of ttss blocks for the present_sym_irrep

    logical ::                  &
      ttss_info_init = .false.                ! status of ttss block type

!   real(8), pointer ::  &
!     xxx(:),            & ! 
!     xxx(:,:,:)           !
    integer, allocatable ::     &
      ttss_block_length(:,:,:), &             ! length of each ttss block per symmetry irrep and CI space (the latter currently always 1)
      ttss_vec_split(:,:),      &             ! vector split into ttss blocks for a given symmetry and CI space
      ttss_block_type(:),       &             ! ttss block type in a given symmetry
      ttss_block_nr(:,:),       &             ! total number of ttss blocks per symmetry irrep and CI space (the latter currently always 1)
      ttss_vec_length(:,:)                    ! total vector length per symmetry irrep and CI space (the latter currently always 1)

  end type ttss_block_structure

! ttss block type object
  type(ttss_block_structure), public , save:: ttss_info

contains 

  subroutine ttss_init(A, mx_ttss, nirreps, sum_of_ci_spaces)

!   ----------------------------------------------------------------------------
    type(ttss_block_structure) :: A
    integer, intent(in)        :: mx_ttss, nirreps, sum_of_ci_spaces
!   ----------------------------------------------------------------------------

!   reset old type information
    call ttss_free(A)

    A%ttss_info_init     = .true.
    A%mx_nr_ttss         = mx_ttss
    A%nr_ptg_irreps      = nirreps
    A%nr_ci_spaces       = sum_of_ci_spaces
    A%present_sym_irrep  = -1
    A%present_ci_space   = -1
    A%total_present_ttss = -1 
    A%total_present_vec  = -1 

    allocate(A%ttss_block_length(A%mx_nr_ttss,A%nr_ptg_irreps,A%nr_ci_spaces))
    allocate(A%ttss_block_nr(A%nr_ptg_irreps,A%nr_ci_spaces))
    allocate(A%ttss_vec_length(A%nr_ptg_irreps,A%nr_ci_spaces))
    allocate(A%ttss_vec_split(8,A%mx_nr_ttss))
    allocate(A%ttss_block_type(A%nr_ptg_irreps))

  end subroutine ttss_init

  subroutine ttss_free(A)

!   ----------------------------------------------------------------------------
    type(ttss_block_structure) :: A
!   ----------------------------------------------------------------------------

    if(.not. A%ttss_info_init) return

    A%ttss_info_init = .false.
    deallocate(A%ttss_block_type)
    deallocate(A%ttss_vec_split)
    deallocate(A%ttss_vec_length)
    deallocate(A%ttss_block_length)
    deallocate(A%ttss_block_nr)

  end subroutine ttss_free

  subroutine ttss_switch(A, current_irrep, current_ci_space)

!   ----------------------------------------------------------------------------
    type(ttss_block_structure) :: A
    integer, intent(in)        :: current_irrep
    integer, intent(in)        :: current_ci_space
!   ----------------------------------------------------------------------------

    A%present_sym_irrep  = current_irrep
    A%present_ci_space   = current_ci_space
    A%total_present_ttss = A%ttss_block_nr(A%present_sym_irrep,A%present_ci_space)
    A%total_present_vec  = A%ttss_vec_length(A%present_sym_irrep,A%present_ci_space)

  end subroutine ttss_switch

end module ttss_block_module
