!dalton_copyright_start
!
!
!dalton_copyright_end

module par_cfg

! stefan: this module contains all types needed for the i/o model of the
!         parMCmod library.
!         

  implicit none

  save

! parameters

  integer, parameter, public :: max_number_of_parallel_files = 20

! file type
  type parallel_files

    integer            ::   &
      nr_parallel_files       ! number of parallel files

    integer,   pointer ::   &
      file_handle(:),       & ! 
      file_offset_factor(:),& !
      nr_file_entries(:)      !

    integer(8), pointer ::  &
      file_offsets(:)         !

  end type parallel_files

end module
