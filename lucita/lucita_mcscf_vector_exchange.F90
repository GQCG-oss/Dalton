!dalton_copyright_start
!
!
!dalton_copyright_end

module lucita_mcscf_vector_exchange

   use ttss_block_module

   implicit none

   public vector_exchange_driver

!  type definition
   type exchange_files

     integer            ::               &
       present_sym_irrep,                &  ! current active symmetry irrep in use
       present_fh_lu,                    &  ! file handle for LUCITA file
       present_fh_mc,                    &  ! file handle for MCSCF  file
       mc_offset,                        &  ! offset to mc types in exchange_files_generic ==> must be equal to max #/2 (1/2 <= lucita; 1/2 >  mcscf)
       push_pull_switch                     ! push or pull data to/from lucita i/o files (1 or 2)

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
   type(exchange_files) :: exchange_f_info


   private

!  print unit
#include "priunit.h"

contains

!*******************************************************************************

   subroutine vector_exchange_driver(exchange_type,              &
                                     vector_type,                &
                                     nr_vectors,                 &
                                     vector_symmetry,            &
                                     io2io_exchange,             &
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
      real(8), intent(inout), optional    :: ymat(*)         ! prepare for core-mem 2 core-mem + core-mem 2 i/o exchange in one shot...
      integer, intent(in)                 :: exchange_type
      integer, intent(in)                 :: vector_type
      integer, intent(in)                 :: nr_vectors
      integer, intent(in)                 :: vector_symmetry
      logical, intent(in)                 :: io2io_exchange
!-------------------------------------------------------------------------------

!     check for initialization of exchange_files_info
      if(.not. exchange_f_info%exchange_file_init) call exchange_f_init(exchange_f_info)

!     switch to ttss-block info for symmetry vector_symmetry in CI space 1
      call ttss_switch(ttss_info,                                &
                       vector_symmetry,                          &
                       1)
!     switch to exchange_files info for symmetry vector_symmetry
      call exchange_f_switch(exchange_f_info,                    &
                             exchange_type,                      &
                             vector_symmetry,                    &
                             vector_type,                        &
                             io2io_exchange)

       call push_pull_vector_exchange_driver(xmat,               &
                                             nr_vectors,         &
                                             exchange_f_info,    &
                                             ttss_info)

   end subroutine vector_exchange_driver
!*******************************************************************************

   subroutine push_pull_vector_exchange_driver(xmat,              &
                                               nr_vectors,        &
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

      do ! loop over vectors

        current_vector = current_vector + 1

        if(current_vector > nr_vectors) exit

        select case(A%exchange_file_io2io)
          case(.true.)
            call push_pull_vector_exchange_io2io(xmat,                      &
                                                 current_vector,            &
                                                 A,                         &
                                                 B)
          case(.false.)
            current_offset = current_offset + (current_vector - 1) * B%total_present_vec
            call push_pull_vector_exchange_memory(xmat,                     &
                                                  current_vector,           &
                                                  current_offset,           &
                                                  A,                        &
                                                  B)
        end select
!
      end do

!     close file(s) if they were opened only for the present task
      if(.not. A%exchange_file_open(1)) call gpclose(A%present_fh_lu,'KEEP')
      if(.not. A%exchange_file_open(2) .and. A%exchange_file_io2io) call gpclose(A%present_fh_mc,'KEEP')

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
            write(lupri,*) ' pull/push (1 or 2) block from/to file to/from offset ==> ',A%push_pull_switch,         &
                             current_block,ioff
            call wrtmatmn(xmat(ioff),1,B%ttss_block_length(current_block,B%present_sym_irrep,B%present_ci_space),1, &
                          B%ttss_block_length(current_block,B%present_sym_irrep,B%present_ci_space),lupri)
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

!   ----------------------------------------------------------------------------
    type(exchange_files) :: A
!   ----------------------------------------------------------------------------

    A%exchange_file_init            = .true.
    A%exchange_file_io2io           = .false.
    A%exchange_file_open(1)         = .false.
    A%exchange_file_open(2)         = .false.
    A%present_sym_irrep             = -1
    A%present_fh_lu                 = -1
    A%present_fh_mc                 = -1
    A%push_pull_switch              = -1
    A%mc_offset                     =  4 ! must be equal to the max number of generic files / 2 (1/2 == lucita; 1/2 == mcscf)

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
                               exchange_type,           &
                               vector_symmetry,         &
                               vector_type,             &
                               io2io_exchange)

!   ----------------------------------------------------------------------------
    type(exchange_files) :: A
    integer, intent(in)  :: exchange_type
    integer, intent(in)  :: vector_symmetry
    integer, intent(in)  :: vector_type
    logical, intent(in)  :: io2io_exchange
!   ----------------------------------------------------------------------------
    integer              :: dummy
!   ----------------------------------------------------------------------------

    A%present_sym_irrep     = vector_symmetry
    A%exchange_file_io2io   = io2io_exchange
    A%push_pull_switch      = exchange_type
    A%exchange_file_open(1) = .false.
    A%exchange_file_open(2) = .false.
    A%present_fh_lu         = -1
    A%present_fh_mc         = -1

!   step 1: LUCITA file
    write(A%exchange_files_generic(vector_type),'(a13,a1)') A%exchange_files_generic(vector_type),            &
                                                            A%exchange_files_f_extension(A%present_sym_irrep)

    write(lupri,*) ' file name set for sym ==> ',A%present_sym_irrep,A%exchange_files_generic(vector_type)

    inquire(opened=A%exchange_file_open(1),file=A%exchange_files_generic(vector_type),number=A%present_fh_lu)

    write(lupri,*) ' file found and status ==> ',A%present_fh_lu,A%exchange_file_open(1)

    if(.not. A%exchange_file_open(1))then
      call gpopen(A%present_fh_lu,A%exchange_files_generic(vector_type),'OLD',' ','UNFORMATTED',dummy,.FALSE.)
    end if
    rewind A%present_fh_lu

!   step 2: MCSCF file
    if(A%exchange_file_io2io)then

      inquire(opened=A%exchange_file_open(2),file=A%exchange_files_generic(vector_type+A%mc_offset), & 
              number=A%present_fh_mc)

      write(lupri,*) ' file found and status ==> ',A%present_fh_mc,A%exchange_file_open(2)

      if(.not. A%exchange_file_open(2))then
        call gpopen(A%present_fh_mc,A%exchange_files_generic(vector_type+A%mc_offset),'OLD',' ','UNFORMATTED',    &
                    dummy,.FALSE.)
      end if
      rewind A%present_fh_mc
    end if

  end subroutine exchange_f_switch
!*******************************************************************************

end module
