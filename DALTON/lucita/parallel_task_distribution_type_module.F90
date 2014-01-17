!dalton_copyright_start
!
!
!dalton_copyright_end
!
! stefan dec 2011: this module contains the type definition for the
!                  MPI communcator handle used in the MCSCF-LUCITA interface + LUCITA program.
!
module parallel_task_distribution_type_module

  implicit none

  public parallel_task_distribution
  public parallel_task_distribution_init_lucipar
  public parallel_task_distribution_free_lucipar
  public parallel_task_distribution_switch_lucipar


! type definition
  type parallel_task_distribution


    integer            ::       &
      active_csym                       = -1,&    ! active csym
      active_csym_max                   = -1      ! max active csym

    logical ::                  &
      parallel_task_distribution_init   = .false. ! status of the communicator type
    logical ::                  &
      parallel_task_distribution_set    = .false. ! status of the parallel_task_list

    integer, allocatable ::     &
      c2s_connections(:,:),     &                 ! 
      parallel_task_list(:,:)                     ! 

  end type parallel_task_distribution

! ttss block type object
  type(parallel_task_distribution), public, save :: ptask_distribution
  integer, parameter, private              :: max_symmetry_distributions_active =  8 ! place holder for # symmetry irrep 
!                                                                                    - if we activate LUCITA for respone this needs to be adapted...

contains 

  subroutine parallel_task_distribution_init_lucipar(A, number_of_blocks, number_of_irreps)

!   ----------------------------------------------------------------------------
    type(parallel_task_distribution) :: A
    integer, intent(in)              :: number_of_blocks
    integer, intent(in)              :: number_of_irreps
!   ----------------------------------------------------------------------------

!   reset old type information
    call parallel_task_distribution_free_lucipar(A)

!   if(number_of_irreps > max_symmetry_distributions_active) call quit('no multi-sym distribution active yet')

    A%parallel_task_distribution_init = .true.
    A%parallel_task_distribution_set  = .false.
    A%active_csym                     = number_of_irreps
    A%active_csym_max                 = number_of_irreps

    if(number_of_irreps > max_symmetry_distributions_active) call quit('# of sym-irrep exceeds max #')

    allocate(A%c2s_connections(number_of_blocks,number_of_irreps))
    allocate(A%parallel_task_list(number_of_blocks,number_of_irreps))
    A%c2s_connections      =  0
    A%parallel_task_list   = -2

  end subroutine parallel_task_distribution_init_lucipar

  subroutine parallel_task_distribution_free_lucipar(A)

!   ----------------------------------------------------------------------------
    type(parallel_task_distribution) :: A
!   ----------------------------------------------------------------------------

    if(.not. A%parallel_task_distribution_init) return

    A%parallel_task_distribution_init = .false.
    A%parallel_task_distribution_set  = .false.
    deallocate(A%c2s_connections)
    deallocate(A%parallel_task_list)

  end subroutine parallel_task_distribution_free_lucipar

  subroutine parallel_task_distribution_switch_lucipar(A,active_irrep)

!   ----------------------------------------------------------------------------
    type(parallel_task_distribution)  :: A
    integer, intent(in)               :: active_irrep
!   ----------------------------------------------------------------------------

!   CALL THIS ROUTINE HERE IF LUCITA IS ACTIVATED FOR RESPONSE...
    A%active_csym = active_irrep
    if(A%active_csym > max_symmetry_distributions_active) call quit('# of sym-irrep exceeds max #')

  end subroutine parallel_task_distribution_switch_lucipar

end module parallel_task_distribution_type_module
