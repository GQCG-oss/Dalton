!dalton_copyright_start
!
!
!dalton_copyright_end

#ifdef VAR_MPI
module file_io_model

! stefan: - this module provides all necessary functionality
!           to setup a file i/o model in parallel mcscf/ci 
!           calculations.
!
!           written by sknecht, may 2007 for DIRAC MCSCF/KR-CI/LUCITA
!           adapted for DALTON by sknecht, november 2010.
  use mpi

#if defined (VAR_INT64)
#define my_MPI_INTEGER MPI_INTEGER8
#else
#define my_MPI_INTEGER MPI_INTEGER4
#endif

  implicit none

  public setup_file_io_model
  public close_file_io_model

  private

  save

  integer, private                       :: istat(MPI_STATUS_SIZE)
  integer, private                       :: ierr
  integer, parameter                     :: flabel_length = 14

contains 
 
  subroutine setup_file_io_model(communicator_io_group,          &
                                 nr_files,                       &
                                 fh_array,                       &
                                 group_id,                       &
                                 group_io_size,                  &
                                 file_identification,            &
                                 print_unit)
!******************************************************************************
!
!    purpose:  return valid fh for parallel MPI-I/O in CI/MCSCF runs:
!                - provides: 
!                           a. file handles
!                           b. opened files ready for reading and writing
!                - requires: 
!                           a. global (group) communication handle
!                           a. file identification string
!
!    for a detailed description of the I/O model see references:
!    S. Knecht, H. J. Aa. Jensen, and T. Fleig
!       JCP, 128, 014108 (2008)
!       JCP, 132, 014108 (2008)
!
!*******************************************************************************
     integer,           intent(in)    :: communicator_io_group
     integer,           intent(in)    :: nr_files
     integer,           intent(inout) :: fh_array(nr_files)
     integer,           intent(in)    :: group_id
     integer,           intent(in)    :: group_io_size
     integer,           intent(in)    :: print_unit
     character (len=5), intent(in)    :: file_identification
!-------------------------------------------------------------------------------
     integer                          :: i  
     integer                          :: j  
     integer                          :: file_info_obj
     integer(kind=MPI_OFFSET_KIND)    :: displacement
     character (len=  4)              :: file_info_groupsz
     character (len=  4)              :: fstring
     character (len=  4)              :: gstring
     character (len= flabel_length)   :: flabel
!-------------------------------------------------------------------------------

!     initial displacement in newly created file
      displacement = 0

!     group id appended to each file
      call int2char_converter(group_id,gstring) 

!     file info object - provide hints to the MPI implementation
      call mpi_info_create(file_info_obj,ierr)

!     1. number of processes sharing the following MPI-I/O files
      write(file_info_groupsz,'(i4)') group_io_size
      call mpi_info_set(file_info_obj,"nb_proc",file_info_groupsz,ierr)
!
#ifdef VAR_PFS
!     2. extra information on IBMs GPFS to enhance I/O performance
      call mpi_info_set(file_info_obj,"IBM_largeblock_io","true",ierr)
#endif
 
      do i = 1, nr_files

!       step a. setup individual file identifier
        call int2char_converter(i,fstring)

!       step b. determine full file name
        write(flabel,'(a5,a4,a1,a4)') file_identification,fstring,'.',gstring
      
!       step c. open the file
        call mpi_file_open(communicator_io_group,flabel(1:flabel_length),              &
                           MPI_MODE_CREATE + MPI_MODE_RDWR + MPI_MODE_DELETE_ON_CLOSE, &
                           file_info_obj,fh_array(i),ierr)
!       step d. set fileview
        call mpi_file_set_view(fh_array(i),displacement,MPI_REAL8,MPI_REAL8,           &
                               "native",file_info_obj,ierr)
      end do

!     free info object
      call mpi_info_free(file_info_obj,ierr)

  end subroutine setup_file_io_model
!*******************************************************************************

  subroutine close_file_io_model(nr_files,                &
                                 fh_offset,               &
                                 fh_array)                 
!*******************************************************************************
!
!    purpose: close MPI-I/O files and "nullify" file handles.
!
!*******************************************************************************
     integer, intent(in )   :: nr_files
     integer, intent(in )   :: fh_offset
     integer, intent(inout) :: fh_array(nr_files+fh_offset)
!-------------------------------------------------------------------------------
     integer                :: i
!-------------------------------------------------------------------------------

      do i = 1, nr_files
        call mpi_file_close(fh_array(i+fh_offset),ierr)
      end do
     
  end subroutine close_file_io_model
!*******************************************************************************
  
  subroutine int2char_converter(int_number,string_rep)
!*******************************************************************************
!     
!  purpose: convert positive integer number int_number (< 10 000)
!           into the 4-byte string string_rep.
!     
!           based on the routine num2str originally written 
!           by C.V. Larsen in Dirac.
!
!*******************************************************************************
     integer,             intent(in )   :: int_number
     character (len = 4), intent(inout) :: string_rep
!-------------------------------------------------------------------------------
     character (len = 1)                :: tmp_str(4)
     integer                            :: num1
     integer                            :: num2
     integer                            :: num3
     integer                            :: num4
!-------------------------------------------------------------------------------

      num3 = int_number
      num2 = 1
      num1 = 1000

      do
        num4          = num3/num1
        tmp_str(num2) = char(num4 + 48)
        num3          = mod(num3,num1)
        num2          = num2 + 1
        num1          = num1/10

        if(num1 < 1) exit
      end do
 
      string_rep=tmp_str(1)//tmp_str(2)//tmp_str(3)//tmp_str(4)
 
  end subroutine int2char_converter

#else 
module dummy_fhio_model
#endif
end module
