!#define LUCI_DEBUG
!dalton_copyright_start
!
!
!dalton_copyright_end

module lucita_mcscf_vector_exchange

   use ttss_block_module
   use file_type_module, only : file_info, file_type
   use vector_xc_file_type
#ifdef VAR_MPI
#ifdef USE_MPI_MOD_F90
   use mpi
#endif
   use parallel_task_distribution_type_module
   use sync_coworkers
   use lucita_cfg
#endif

   implicit none

#ifdef VAR_MPI
#ifndef USE_MPI_MOD_F90
#include "mpif.h"
#endif
   public vector_exchange_interface_cw
#endif
   public vector_exchange_driver

   private

#ifdef VAR_MPI
   integer(kind=MPI_INTEGER_KIND)         :: my_MPI_REAL8 = MPI_REAL8
   integer(kind=MPI_INTEGER_KIND)         :: ierr_mpi, fh_mpi, len_mpi, my_mpi_sum, root_mpi = 0
#endif

   integer, parameter, private :: mc_offset = 4 ! offset to mc types in exchange_f... ==> must be equal to max #/2 (1/2 <= lucita; 1/2 >  mcscf)

!  print unit
#include "priunit.h"

contains

!*******************************************************************************

   subroutine vector_exchange_driver(exchange_type,              &
                                     vector_type,                &
                                     nr_vectors,                 &
                                     vector_symmetry,            &
                                     io2io_exchange,             &
                                     do_vector_exchange,         &
                                     parallel_vector_xc,         &
                                     xmat,                       &
                                     ymat)
                                     
!-------------------------------------------------------------------------------
!
!  purpose: driver for vector exchange (I/O) based between the MCSCF and LUCITA 
!           programs.
!
!           exchange_type == 1: copy vector from LUCITA to MCSCF
!           exchange_type == 2: copy vector from  MCSCF to LUCITA
!
!           vector type   == 1: reference vector / right-hand side vector
!           vector type   == 2: CI start vector(s) / right-hand side vector
!           vector type   == 3: sigma vector / left-hand side vector
!           vector type   == 4: H diagonal
!
!-------------------------------------------------------------------------------
      real(8), intent(inout)              :: xmat(*)
      real(8), intent(inout), optional    :: ymat(*)            ! prepare for core-mem 2 core-mem + core-mem 2 i/o exchange in one shot...
      integer, intent(inout)              :: exchange_type      ! inout because of MPI calling structure
      integer, intent(inout)              :: vector_type        ! inout because of MPI calling structure 
      integer, intent(inout)              :: nr_vectors         ! inout because of MPI calling structure
      integer, intent(inout)              :: vector_symmetry    ! inout because of MPI calling structure
      logical, intent(inout)              :: io2io_exchange     ! inout because of MPI calling structure
      logical, intent(inout)              :: do_vector_exchange ! inout because of MPI calling structure
      logical, intent(in)                 :: parallel_vector_xc
!-------------------------------------------------------------------------------

!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
      write(lupri,*) ' do_vector_exchange ==> ',do_vector_exchange
      write(lupri,*) ' vector_type        ==> ',vector_type
      write(lupri,*) ' nr_vectors         ==> ',nr_vectors
      write(lupri,*) ' vector_symmetry    ==> ',vector_symmetry
#undef LUCI_DEBUG
#endif

!     switch to exchange_files info for symmetry vector_symmetry
      call exchange_f_switch(exchange_f_info,                    &
                             file_info,                          &
                             exchange_type,                      &
                             vector_symmetry,                    &
                             vector_type,                        &
                             nr_vectors,                         &
                             io2io_exchange)

!     switch to ttss-block info for symmetry vector_symmetry in CI space 1
      call ttss_switch(ttss_info,                                &
                       vector_symmetry,                          &
                       1)

      if(do_vector_exchange)then
#ifdef VAR_MPI
!       get hold of co-workers and provide them with the mandatory details
        if(parallel_vector_xc) call vector_exchange_wake_up_cw()
#endif

        call push_pull_vector_exchange_driver(xmat,              &
                                              nr_vectors,        &
                                              parallel_vector_xc,&
                                              exchange_f_info,   &
                                              ttss_info)

      end if

!     reset vector exchange control variables
      do_vector_exchange = .false.
      vector_type        = -1

   end subroutine vector_exchange_driver
!*******************************************************************************

#ifdef VAR_MPI
   subroutine vector_exchange_wake_up_cw()
!-------------------------------------------------------------------------------
!
!       purpose: get hold of co-workers and provide them with the mandatory details
!-------------------------------------------------------------------------------

      lucita_ci_run_id   = 'xc vector   '
      sync_ctrl_array(2) = .true.
      sync_ctrl_array(5) = .true.
      sync_ctrl_array(6) = .true.

      call lucita_start_cw(1)
      call sync_coworkers_cfg()

   end subroutine vector_exchange_wake_up_cw
!*******************************************************************************

   subroutine vector_exchange_interface_cw(active_xc_vector_type,      &
                                           do_vector_exchange,         &
                                           xmat,                       &
                                           ymat)
                                     
!-------------------------------------------------------------------------------
!
!  purpose: co-workers driver for vector exchange (I/O) based between the MCSCF and LUCITA 
!           programs.
!
!           active_xc_vector_type == 1: reference vector / right-hand side vector
!           active_xc_vector_type == 2: CI start vector(s) / right-hand side vector
!           active_xc_vector_type == 3: sigma vector / left-hand side vector
!           active_xc_vector_type == 4: H diagonal
!
!-------------------------------------------------------------------------------
      real(8), intent(inout)              :: xmat(*)
      real(8), intent(inout), optional    :: ymat(*)         ! prepare for core-mem 2 core-mem + core-mem 2 i/o exchange in one shot...
      integer, intent(inout)              :: active_xc_vector_type
      logical, intent(inout)              :: do_vector_exchange
!-------------------------------------------------------------------------------

!     switch to exchange_files info for symmetry vector_symmetry
      call exchange_f_switch(exchange_f_info,                                          &
                             file_info,                                                &
                             exchange_f_info%push_pull_switch,                         &
                             exchange_f_info%present_sym_irrep,                        &
                             active_xc_vector_type,                                    &
                             exchange_f_info%total_nr_vectors,                         &
                             exchange_f_info%exchange_file_io2io)

!     switch to ttss-block info for symmetry vector_symmetry in CI space 1
      call ttss_switch(ttss_info,                                                      &
                       exchange_f_info%present_sym_irrep,                              &
                       1)

#ifdef LUCI_DEBUG
      write(lupri,*) ' parallel node xc-driver call? ==> ',do_vector_exchange
      write(lupri,*) ' parallel node xc-driver active_xc_vector_type? ==> ',active_xc_vector_type
#endif

      if(do_vector_exchange)then
        call push_pull_vector_exchange_driver(xmat,                                    &
                                              exchange_f_info%total_nr_vectors,        &
                                              .true.,                                  &
                                              exchange_f_info,                         &
                                              ttss_info)
      end if

!     reset vector exchange control variables
      do_vector_exchange    = .false.
      active_xc_vector_type = -1

   end subroutine vector_exchange_interface_cw
!*******************************************************************************
#endif

   subroutine push_pull_vector_exchange_driver(xmat,              &
                                               nr_vectors,        &
                                               parallel_vector_xc,&
                                               A,                 &
                                               B)
                                     
!-------------------------------------------------------------------------------
!
!  purpose: vector exchange driver [(I/O) - memory or (I/O) - (I/O) based] between 
!           the MCSCF and LUCITA programs - push-pull version. 
!
!           pull: transfer vector(s) from LUCITA i/o files to MCSCF core-memory / i/o
!           push: transfer vector(s) from MCSCF core-memory / i/o to LUCITA i/o files 
!
!-------------------------------------------------------------------------------
      real(8), intent(inout)     :: xmat(*)
      integer, intent(in)        :: nr_vectors
      logical, intent(in)        :: parallel_vector_xc
      type(exchange_files)       :: A     
      type(ttss_block_structure) :: B     
!-------------------------------------------------------------------------------
      integer                    :: ioff
      integer                    :: current_vector
      integer                    :: current_offset
!------------------------------------------------------------------------------

      current_vector = 0
      current_offset = 0
      ioff           = 1
!#define LUCI_DEBUG

      do ! loop over vectors

        current_vector = current_vector + 1

        if(current_vector > nr_vectors) exit

        select case(A%exchange_file_io2io)
          case(.true.)
#ifdef LUCI_DEBUG
            write(lupri,*)'io2io exchange'
#endif
#ifdef VAR_MPI
            if(parallel_vector_xc)then
!             call push_pull_vector_exchange_io2io_parallel()
            else
#endif
              call push_pull_vector_exchange_io2io(xmat,                      &
                                                   current_vector,            &
                                                   A,                         &
                                                   B)
#ifdef VAR_MPI
            end if
#endif
          case(.false.)
            current_offset = current_offset + (current_vector - 1) * B%total_present_vec

#ifdef LUCI_DEBUG
            write(lupri,*)'parallel_vector_xc',parallel_vector_xc
            write(lupri,*)'current_offset    ',current_offset
#endif

#ifdef VAR_MPI
            if(parallel_vector_xc)then
              call push_pull_vector_exchange_memory_parallel(xmat,               &
                                                             current_vector,     &
                                                             current_offset,     &
                                                             A,                  &
                                                             B,                  &
                                                             ptask_distribution, &
                                                             file_info)
            else
#endif
              call push_pull_vector_exchange_memory(xmat,                     &
                                                    current_vector,           &
                                                    current_offset,           &
                                                    A,                        &
                                                    B)
#ifdef VAR_MPI
            end if
#endif
        end select
!
      end do

#undef LUCI_DEBUG
      if(A%my_process_id == 0)then
!       close file(s) if they were opened only for the present task
        if(.not. A%exchange_file_open(1)) call gpclose(A%present_fh_lu,'KEEP')
        if(.not. A%exchange_file_open(2) .and. A%exchange_file_io2io) call gpclose(A%present_fh_mc,'KEEP')
      end if

   end subroutine push_pull_vector_exchange_driver
!*******************************************************************************

   subroutine push_pull_vector_exchange_memory(xmat,             &
                                               current_vector,   &
                                               current_offset,   &
                                               A,                &
                                               B)
                                     
!-------------------------------------------------------------------------------
!
!  purpose: vector exchange (I/O) - memory based between the MCSCF and LUCITA 
!           programs - push-pull version. 
!
!           pull: transfer vector(s) from LUCITA i/o files to MCSCF core-memory
!           push: transfer vector(s) from MCSCF core-memory to LUCITA i/o files
!
!-------------------------------------------------------------------------------
      real(8), intent(inout)     :: xmat(*)
      integer, intent(in)        :: current_vector
      integer, intent(in)        :: current_offset
      type(exchange_files)       :: A     
      type(ttss_block_structure) :: B     
!-------------------------------------------------------------------------------
      integer                    :: ioff
      integer                    :: current_block
      integer                    :: is0
      integer                    :: packing
      integer                    :: block_length_rw
      integer, parameter         :: no_zeroing  = 0
!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
#include "priunit.h"
#endif
!------------------------------------------------------------------------------

      ioff          = 1 + current_offset
      current_block = 0 

      do ! loop over blocks

        current_block  = current_block + 1 
        if(current_block > B%total_present_ttss) exit

        select case(A%push_pull_switch)
          case(1) 
!           pull block from file to core memory
            call ifrmds(block_length_rw,1,-1,A%present_fh_lu)
            call frmdsc2(xmat(ioff),block_length_rw,-1,A%present_fh_lu,is0,packing,no_zeroing)
          case(2) 
!           push block from core memory to file
            call itods(B%ttss_block_length(current_block,B%present_sym_irrep,B%present_ci_space),                 &
                       1,-1,A%present_fh_lu)
            call todsc_luci(xmat(ioff),B%ttss_block_length(current_block,B%present_sym_irrep,B%present_ci_space), &
                            -1,A%present_fh_lu)
        end select
#ifdef LUCI_DEBUG
            write(lupri,*) 'seq  pull/push (1 or 2) block from/to file to/from offset ==> ',A%push_pull_switch,   &
                             current_block,ioff
!           call wrtmatmn(xmat(ioff),1,B%ttss_block_length(current_block,B%present_sym_irrep,B%present_ci_space),1, &
!                         B%ttss_block_length(current_block,B%present_sym_irrep,B%present_ci_space),lupri)
#undef LUCI_DEBUG
#endif

!       keep track of core-memory offset
        ioff = ioff + B%ttss_block_length(current_block,B%present_sym_irrep,B%present_ci_space)

      end do

      select case(A%push_pull_switch)
        case(1) 
!         skip end-of-vector marker on file A%present_fh_lu
          call ifrmds(block_length_rw,1,-1,A%present_fh_lu)
        case(2) 
!         set end-of-vector marker on file A%present_fh_lu
          call itods(-1,1,-1,A%present_fh_lu)
      end select

!
   end subroutine push_pull_vector_exchange_memory
!*******************************************************************************

#ifdef VAR_MPI
   subroutine push_pull_vector_exchange_memory_parallel(xmat,             &
                                                        current_vector,   &
                                                        current_offset,   &
                                                        A,                &
                                                        B,                &
                                                        C,                &
                                                        D)
                                     
!-------------------------------------------------------------------------------
!
!  purpose: vector exchange (I/O) - memory based between the MCSCF and LUCITA 
!           programs - push-pull version. 
!
!           pull: transfer vector(s) from LUCITA i/o files to MCSCF core-memory
!           push: transfer vector(s) from MCSCF core-memory to LUCITA i/o files
!
!-------------------------------------------------------------------------------
      real(8), intent(inout)           :: xmat(*)
      integer, intent(in)              :: current_vector
      integer, intent(in)              :: current_offset
      type(exchange_files)             :: A
      type(ttss_block_structure)       :: B
      type(parallel_task_distribution) :: C
      type(file_type)                  :: D
!-------------------------------------------------------------------------------
      real(8), external                :: ddot
      real(8)                          :: checkdot
      integer                          :: ioff
      integer                          :: my_STATUS(mpi_status_size)
      integer                          :: current_block
      integer                          :: block_length_rw
      integer(kind=mpi_offset_kind)    :: ioffset
      integer(kind=mpi_offset_kind)    :: ioffset_scratch
      integer                          :: ioffset_int
#include "parluci.h"
!------------------------------------------------------------------------------

      select case(A%push_pull_switch)
        case(1) 
          call dzero(xmat(1+current_offset),B%total_present_vec)
        case(2) 
          len_mpi = B%total_present_vec
          call mpi_bcast(xmat(1+current_offset),len_mpi,my_MPI_REAL8,         &
                         root_mpi,mpi_comm_world,ierr_mpi)
      end select

!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
      write(lupri,*) ' present fh handle ==> ',A%present_fh_par, A%push_pull_switch
#undef LUCI_DEBUG
#endif
!
!     initialize
      ioff             = 1 + current_offset
      current_block    = 0 
      ioffset          = 0
      ioffset_int      = 0

!     calculate offset
      ioffset          = D%file_offsets(A%present_vector_type) + (current_vector - 1) * my_vec2_ioff
      ioffset_int      = 1                                     + (current_vector - 1) * my_act_blk2

      do ! loop over blocks

        current_block   = current_block + 1 
        if(current_block > B%total_present_ttss) exit

        block_length_rw = B%ttss_block_length(current_block,B%present_sym_irrep,B%present_ci_space)
!       keep track of core-memory offset
        ioff            = ioff + block_length_rw

        if(C%parallel_task_list(current_block,A%present_sym_irrep) /= A%my_process_id) cycle

        ioffset_scratch = block_length_rw

        select case(A%push_pull_switch)
          case(1) ! pull block from file to core memory
            if(D%iluxlist(ioffset_int,A%present_vector_type) > 0)then
              fh_mpi  = A%present_fh_par
              len_mpi = block_length_rw
              call mpi_file_read_at(fh_mpi,ioffset,xmat(ioff-block_length_rw),          &
                                    len_mpi,my_mpi_real8,my_STATUS,ierr_mpi)
            end if
          case(2) ! push block from core memory to file
            D%iluxlist(ioffset_int,A%present_vector_type) = 0
            checkdot                                      = ddot(block_length_rw,xmat(ioff-block_length_rw),1,        &
                                                                                 xmat(ioff-block_length_rw),1)
            if(checkdot > 0.0d0)then
              fh_mpi  = A%present_fh_par
              len_mpi = block_length_rw
              call mpi_file_write_at(fh_mpi,ioffset,xmat(ioff-block_length_rw),         &
                                     len_mpi,my_mpi_real8,my_STATUS,ierr_mpi)
              D%iluxlist(ioffset_int,A%present_vector_type) = 1
            end if
        end select
!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
            write(lupri,*) ' pull/push (1 or 2) block from/to file to/from fh handle offset ==> ',A%push_pull_switch,         &
                             current_block,A%present_fh_par,ioff-block_length_rw
!           call wrtmatmn(xmat(ioff-block_length_rw),1,block_length_rw,1,block_length_rw,lupri)
#undef LUCI_DEBUG
#endif

!       keep track of file / file-list offsets
        ioffset     = ioffset     + ioffset_scratch
        ioffset_int = ioffset_int + 1

      end do

      select case(A%push_pull_switch)
        case(1) 
          if(A%my_process_id > 0)then
!           write(lupri,*) 'A%my_process_id ',A%my_process_id, B%total_present_vec, current_offset
            len_mpi = B%total_present_vec
            my_mpi_sum = mpi_sum
            call mpi_reduce(xmat(1+current_offset),mpi_in_place,len_mpi,my_MPI_REAL8,      &
                            my_mpi_sum,root_mpi,mpi_comm_world,ierr_mpi)
          else
!           write(lupri,*) 'A%my_process_id ',A%my_process_id, B%total_present_vec, current_offset
            len_mpi = B%total_present_vec
            my_mpi_sum = mpi_sum
            call mpi_reduce(mpi_in_place,xmat(1+current_offset),len_mpi,my_MPI_REAL8,      &
                            my_mpi_sum,root_mpi,mpi_comm_world,ierr_mpi)
          end if
        case(2) 
          ! nothing to do
      end select
!
   end subroutine push_pull_vector_exchange_memory_parallel
!*******************************************************************************
#endif /* VAR_MPI */

   subroutine push_pull_vector_exchange_io2io(xmat,              &
                                              current_vector,    &
                                              A,                 &
                                              B)
                                     
!-------------------------------------------------------------------------------
!
!  purpose: vector exchange (I/O) - (I/O) based between the MCSCF and LUCITA 
!           programs - push-pull version. 
!
!           pull: transfer vector(s) from LUCITA i/o files to MCSCF i/o files   
!           push: transfer vector(s) from MCSCF i/o files to LUCITA i/o files
!
!-------------------------------------------------------------------------------
      real(8), intent(inout)     :: xmat(*)
      integer, intent(in)        :: current_vector
      type(exchange_files)       :: A     
      type(ttss_block_structure) :: B     
!-------------------------------------------------------------------------------
      integer                    :: ioff
      integer                    :: current_block
      integer                    :: is0
      integer                    :: packing
      integer                    :: block_length_rw
      integer, parameter         :: no_zeroing  = 0
!------------------------------------------------------------------------------

      ioff          = 1
      current_block = 0 

      select case(A%push_pull_switch)
        case(1) 
!         nothing to do... 
        case(2)
!         read current vector from mc file
          call readt(A%present_fh_mc,B%total_present_vec,xmat)
      end select

      do ! loop over blocks

        current_block  = current_block + 1 
        if(current_block > B%total_present_ttss) exit

        select case(A%push_pull_switch)
          case(1) 
!           pull block from lucita file to mc core memory
            call ifrmds(block_length_rw,1,-1,A%present_fh_lu)
            call frmdsc2(xmat(ioff),block_length_rw,-1,A%present_fh_lu,is0,packing,no_zeroing)
          case(2) 
!           push block from mc core memory to lucita file
            call itods(B%ttss_block_length(current_block,B%present_sym_irrep,B%present_ci_space),                 &
                       1,-1,A%present_fh_lu)
            call todsc_luci(xmat(ioff),B%ttss_block_length(current_block,B%present_sym_irrep,B%present_ci_space), &
                            -1,A%present_fh_lu)
        end select

!       keep track of core-memory offset
        ioff = ioff + B%ttss_block_length(current_block,B%present_sym_irrep,B%present_ci_space)

      end do

      select case(A%push_pull_switch)
        case(1) 
!         write current vector to mc file
          call writt(A%present_fh_mc,max(4,B%total_present_vec),xmat)
!         skip end-of-vector marker on file A%present_fh_lu
          call ifrmds(block_length_rw,1,-1,A%present_fh_lu)
        case(2) 
!         set end-of-vector marker on file A%present_fh_lu
          call itods(-1,1,-1,A%present_fh_lu)
      end select
!
   end subroutine push_pull_vector_exchange_io2io
!*******************************************************************************

  subroutine exchange_f_init(A)

#include "maxorb.h"
#include "infpar.h"
!   ----------------------------------------------------------------------------
    type(exchange_files) :: A
!   ----------------------------------------------------------------------------

    A%exchange_file_init            = .true.
    A%my_process_id                 = mytid
    A%exchange_file_open(1)         = .false.
    A%exchange_file_open(2)         = .false.
    A%present_fh_lu                 = -1
    A%present_fh_par                = -1
    A%present_fh_mc                 = -1

    if(A%my_process_id == 0)then
      A%exchange_file_io2io           = .false.
      A%push_pull_switch              = -1
      A%total_nr_vectors              = -1
      A%present_sym_irrep             = -1
      A%present_vector_type           = -1
    end if

    A%exchange_files_generic(1)     = 'LUCITA_CEPVC.x'
    A%exchange_files_generic(2)     = 'LUCITA_CVECS.x'
    A%exchange_files_generic(3)     = 'LUCITA_HCVEC.x'
    A%exchange_files_generic(4)     = 'LUCITA_HDIAG.x'
    A%exchange_files_generic(5)     = 'SIRIUS.RST    '
    A%exchange_files_generic(6)     = 'SIRIUS.BVECS  '
    A%exchange_files_generic(7)     = 'SIRIUS.E2BVECS'
    A%exchange_files_generic(8)     = 'SIRIUS.E2DIAG '

    A%exchange_files_f_extension(1) = 'a'
    A%exchange_files_f_extension(2) = 'b'
    A%exchange_files_f_extension(3) = 'c'
    A%exchange_files_f_extension(4) = 'd'
    A%exchange_files_f_extension(5) = 'e'
    A%exchange_files_f_extension(6) = 'f'
    A%exchange_files_f_extension(7) = 'g'
    A%exchange_files_f_extension(8) = 'h'

  end subroutine exchange_f_init
!*******************************************************************************

  subroutine exchange_f_switch(A,                       &
                               B,                       &
                               exchange_type,           &
                               vector_symmetry,         &
                               vector_type,             &
                               vector_quantity,         &
                               io2io_exchange)

!   ----------------------------------------------------------------------------
    type(exchange_files)   :: A
    type(file_type)        :: B
    integer, intent(inout) :: exchange_type    ! inout because co-workers are calling with the type argument in the calling list
    integer, intent(inout) :: vector_symmetry  ! inout because co-workers are calling with ...
    integer, intent(inout) :: vector_type      ! inout because co-workers are calling with ...
    integer, intent(inout) :: vector_quantity  ! inout because co-workers are calling with ...
    logical, intent(inout) :: io2io_exchange   ! inout because co-workers are calling with ...
!   ----------------------------------------------------------------------------
    integer                :: dummy
!   ----------------------------------------------------------------------------

    ! check for initialization of type exchange_files
    if(.not. A%exchange_file_init) call exchange_f_init(A)

    A%present_sym_irrep       = vector_symmetry
    A%present_vector_type     = vector_type
    A%exchange_file_io2io     = io2io_exchange
    A%push_pull_switch        = exchange_type
    A%total_nr_vectors        = vector_quantity
    A%exchange_file_open(1)   = .false.
    A%exchange_file_open(2)   = .false.
    A%present_fh_lu           = -1
    A%present_fh_mc           = -1
    A%present_fh_par          = -1
#ifdef VAR_MPI
    if(allocated(B%fh_lu))then
    A%present_fh_par          = B%fh_lu(A%present_vector_type)
    else
    call quit('LUCITA: access attempt to non-open MPI file')
    end if
#endif


    if(A%my_process_id == 0)then

!     set parallel file handle in type file_type to be potentially used inside LUCITA
      if(B%current_file_nr_active1 < 0)then
         B%current_file_nr_active1 = A%present_vector_type
      else
         B%current_file_nr_active2 = A%present_vector_type
      end if

!     step 1: LUCITA file
      write(A%exchange_files_generic(vector_type),'(a13,a1)') A%exchange_files_generic(vector_type),            &
                                                              A%exchange_files_f_extension(A%present_sym_irrep)

!     write(lupri,*) ' file name set for sym ==> ',A%present_sym_irrep,A%exchange_files_generic(vector_type)
  
      inquire(opened=A%exchange_file_open(1),file=A%exchange_files_generic(vector_type),number=A%present_fh_lu)
  
!     write(lupri,*) ' file found and status ==> ',A%present_fh_lu,A%exchange_file_open(1)

      if(.not. A%exchange_file_open(1))then
        call gpopen(A%present_fh_lu,A%exchange_files_generic(vector_type),'OLD',' ','UNFORMATTED',dummy,.FALSE.)
      end if
      rewind A%present_fh_lu

!     set sequential file handle in type file_type to be potentially used inside LUCITA
      if(B%current_file_fh_seqf(1) < 0)then
         B%current_file_fh_seqf(1) = A%present_fh_lu
      else
         B%current_file_fh_seqf(2) = A%present_fh_lu
      end if
    
!     step 2: MCSCF file
      if(A%exchange_file_io2io)then

        inquire(opened=A%exchange_file_open(2),file=A%exchange_files_generic(vector_type+mc_offset), & 
                number=A%present_fh_mc)

!       write(lupri,*) ' file found and status ==> ',A%present_fh_mc,A%exchange_file_open(2)

        if(.not. A%exchange_file_open(2))then
          call gpopen(A%present_fh_mc,A%exchange_files_generic(vector_type+mc_offset),'OLD',' ','UNFORMATTED',    &
                    dummy,.FALSE.)
        end if
        rewind A%present_fh_mc
      end if
    end if ! global master switch

  end subroutine exchange_f_switch
!*******************************************************************************

end module
!#undef LUCI_DEBUG
