!dalton_copyright_start
!
!
!dalton_copyright_end

module vector_xc_file_type

   implicit none

   public exchange_files

!  type definition
   type exchange_files

     integer            ::               &
       present_sym_irrep,                &  ! current active symmetry irrep in use
       present_vector_type,              &  ! current active vector_type
       present_fh_lu,                    &  ! file handle for LUCITA file
       present_fh_par,                   &  ! MPI file handle for LUCITA file
       present_fh_mc,                    &  ! file handle for MCSCF  file
       push_pull_switch,                 &  ! push or pull data to/from lucita i/o files (1 or 2)
       total_nr_vectors,                 &  ! total # of vectors of a given type
       my_process_id                        ! process ID

     logical ::                          &
       exchange_file_init     = .false., &  ! status of exchange_files_info
       exchange_file_io2io    = .false., &  ! file exchange using i/o 2 i/o [.false. ==> i/o 2 core-mem (ex-type 1) / core-mem 2 i/o (ex-type 2)]
       exchange_file_open(2)                ! leave status 'open' of file(s) after i/o operations

     character(len= 1) ::                &
       exchange_files_f_extension(8)        ! extension for generic exchange file names (LUCITA_CVECS.x, LUCITA_HCVEC.x, etc) 8 <--> D2h symmetry
     character(len=14) ::                &  ! generic exchange file names
       exchange_files_generic(8)

   end type exchange_files

!  exchange_files object
   type(exchange_files), public, save :: exchange_f_info

end module
