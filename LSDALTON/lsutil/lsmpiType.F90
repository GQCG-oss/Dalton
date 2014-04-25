module lsmpi_type
  use precision
  use ptr_assoc_module !, only: ass_44I2to1,ass_48I2to1,ass_84I2to1,ass_88I2to1,&
!      & ass_44I3to1,ass_48I3to1,ass_84I3to1,ass_88I3to1,&
!      & ass_44I4to1,ass_48I4to1,ass_84I4to1,ass_88I4to1,&
!      & ass_4D2to1,ass_8D2to1,ass_4D3to1,ass_8D3to1,&
!      & ass_4D4to1,ass_8D4to1,ass_D2to1,ass_D3to1,ass_D4to1
  use,intrinsic :: iso_c_binding,only:c_ptr,c_f_pointer,c_associated,c_null_ptr
  use Integralparameters
  use memory_handling, only: mem_alloc,mem_dealloc, max_mem_used_global,&
       & longintbuffersize, print_maxmem, stats_mem, copy_from_mem_stats,&
       & init_globalmemvar, stats_mpi_mem, copy_to_mem_stats, &
       & MemModParamPrintMemory, Print_Memory_info, &
       & MemModParamPrintMemorylupri
#ifdef VAR_MPI
  use infpar_module
  use lsmpi_module
#endif

  INTERFACE ls_mpibcast_chunks
    MODULE PROCEDURE ls_mpibcast_realkV_parts44,ls_mpibcast_realkV_parts48,&
                   & ls_mpibcast_realkV_parts84,ls_mpibcast_realkV_parts88
  END INTERFACE ls_mpibcast_chunks

  INTERFACE lsmpi_send
    MODULE PROCEDURE lsmpi_send_realkV_4,lsmpi_send_realkV_8
  END INTERFACE lsmpi_send

  INTERFACE lsmpi_recv
    MODULE PROCEDURE lsmpi_recv_realkV_4,lsmpi_recv_realkV_8
  END INTERFACE lsmpi_recv

  INTERFACE ls_mpibcast
     MODULE PROCEDURE ls_mpibcast_integer,& 
          &           ls_mpibcast_integerV,ls_mpibcast_integerV_wrapper8,&
          &           ls_mpibcast_longV, ls_mpibcast_longV_wrapper8, &
          &           ls_mpibcast_integerM,ls_mpibcast_integerM_wrapper8,&
          &           ls_mpibcast_longM, ls_mpibcast_longM_wrapper8, &
          &           ls_mpibcast_realk, ls_mpibcast_realkV,&
          &           ls_mpibcast_realkV_wrapper8,&
          &           ls_mpibcast_realkM, ls_mpibcast_realkT,&
          &           ls_mpibcast_realkQ,&
          &           ls_mpibcast_logical4,ls_mpibcast_logical8,&
          &           ls_mpibcast_logical4V,ls_mpibcast_logical4V_wrapper8, &
          &           ls_mpibcast_logical8V,ls_mpibcast_logical8V_wrapper8, &
          &           ls_mpibcast_logical4M4,ls_mpibcast_logical4M8,&
          &           ls_mpibcast_logical8M4,ls_mpibcast_logical8M8,&
          &           ls_mpibcast_short, &
          &           ls_mpibcast_charac, &
          &           ls_mpibcast_characV,ls_mpibcast_characV_wrapper8, &
          &           ls_mpibcast_characV2,ls_mpibcast_characV2_wrapper8, &
          &           ls_mpibcast_shortV,ls_mpibcast_shortV_wrapper8,&
          &           ls_mpibcast_long
  END INTERFACE

  INTERFACE lsmpi_reduction
     MODULE PROCEDURE lsmpi_reduction_realk_single, &
          &           lsmpi_reduction_realk, &
          &           lsmpi_reduction_realk_wrapper8, &
          &           lsmpi_reduction_realkM4,lsmpi_reduction_realkM8, &
          &           lsmpi_reduction_realkT4,lsmpi_reduction_realkT8, &
          &           lsmpi_reduction_integer4,lsmpi_reduction_integer4_wrapper8, &
          &           lsmpi_reduction_integer8,lsmpi_reduction_integer8_wrapper8, &
          &           lsmpi_reduction_integer4M4,lsmpi_reduction_integer4M8,&
          &           lsmpi_reduction_integer8M4,lsmpi_reduction_integer8M8
  END INTERFACE lsmpi_reduction

  INTERFACE ls_mpisendrecv
     MODULE PROCEDURE ls_mpisendrecv_integer,& 
          &           ls_mpisendrecv_integerV,ls_mpisendrecv_integerV_wrapper8,&
          &           ls_mpisendrecv_longV,ls_mpisendrecv_longV_wrapper8, &
          &           ls_mpisendrecv_realk, &
          &           ls_mpisendrecv_realkV,ls_mpisendrecv_realkV_wrapper8,&
          &           ls_mpisendrecv_realkM, ls_mpisendrecv_realkT,&
          &           ls_mpisendrecv_realkQ,&
          &           ls_mpisendrecv_logical4,ls_mpisendrecv_logical8, &
          &           ls_mpisendrecv_logical4V,ls_mpisendrecv_logical4V_wrapper8,&
          &           ls_mpisendrecv_logical8V,ls_mpisendrecv_logical8V_wrapper8,&
          &           ls_mpisendrecv_logical4M,ls_mpisendrecv_logical8M,&
          &           ls_mpisendrecv_short, &
          &           ls_mpisendrecv_charac, &
          &           ls_mpisendrecv_characV, ls_mpisendrecv_characV_wrapper8,&
          &           ls_mpisendrecv_characV2,ls_mpisendrecv_characV2_wrapper8, &
          &           ls_mpisendrecv_shortV, ls_mpisendrecv_shortV_wrapper8,&
          &           ls_mpisendrecv_long
  END INTERFACE ls_mpisendrecv

  INTERFACE ls_mpi_buffer
     MODULE PROCEDURE ls_mpi_buffer4_integer,ls_mpi_buffer8_integer,&
          &           ls_mpi_buffer_integer8V,ls_mpi_buffer_integer8V_buf4_wrapper,&
          &           ls_mpi_buffer_integer4V,ls_mpi_buffer_integer4V_buf4_wrapper,&
          &           ls_mpi_buffer_integer4M_wrapper4, ls_mpi_buffer_integer4M_wrapper8,&
          &           ls_mpi_buffer_integer8M_wrapper4, ls_mpi_buffer_integer8M_wrapper8,&
          &           ls_mpi_buffer_integer4T_wrapper4, ls_mpi_buffer_integer4T_wrapper8,&
          &           ls_mpi_buffer_integer8T_wrapper4, ls_mpi_buffer_integer8T_wrapper8,&
          &           ls_mpi_buffer_integer4Q_wrapper4, ls_mpi_buffer_integer4Q_wrapper8,&
          &           ls_mpi_buffer_integer8Q_wrapper4, ls_mpi_buffer_integer8Q_wrapper8,&
          &           ls_mpi_buffer_realk, &
          &           ls_mpi_buffer_realkV, ls_mpi_buffer_realkM, &
          &           ls_mpi_buffer_realkT,&
          &           ls_mpi_buffer_logical, ls_mpi_buffer_logicalV,&
          &           ls_mpi_buffer_logicalM,ls_mpi_buffer_shortinteger, &
          &           ls_mpi_buffer_charac, ls_mpi_buffer_characV, &
          &           ls_mpi_buffer_characV2,&
          &           ls_mpi_buffer_shortintegerV, ls_mpi_buffer_shortintegerM
  END INTERFACE ls_mpi_buffer

  INTERFACE lsmpi_int_reduction
    MODULE PROCEDURE lsmpi_int8_reduction,lsmpi_int4_reduction
  END INTERFACE lsmpi_int_reduction

  INTERFACE lsmpi_sho_reduction
    MODULE PROCEDURE lsmpi_shoV_reduction,lsmpi_shoV_reduction_wrapper8
  END INTERFACE lsmpi_sho_reduction

  INTERFACE lsmpi_local_reduction
     MODULE PROCEDURE lsmpi_local_reduction_realkVN4,&
          & lsmpi_local_reduction_realkVN8,&
          & lsmpi_local_reduction_realkM,lsmpi_local_reduction_realkT,&
          & lsmpi_local_reduction_realkQ,lsmpi_local_reduction_intV,&
          & lsmpi_local_reduction_realk,lsmpi_local_reduction_int, &
          & lsmpi_local_reduction_realkVN8_parts, &
          & lsmpi_local_reduction_realkVN4_parts
  END INTERFACE lsmpi_local_reduction


  INTERFACE lsmpi_allreduce
     MODULE PROCEDURE lsmpi_allreduce_D,&
          & lsmpi_allreduce_D1N4,lsmpi_allreduce_D1N8,&
          & lsmpi_allreduce_D2,&
          & lsmpi_allreduce_D3,lsmpi_allreduce_D4,&
          & lsmpi_allreduce_int4V,lsmpi_allreduce_int4V_wrapper8, &
          & lsmpi_allreduce_int8V,lsmpi_allreduce_int8V_wrapper8, &
          & lsmpi_allreduce_D1N8_parts, lsmpi_allreduce_D1N4_parts
  END INTERFACE lsmpi_allreduce


  interface lsmpi_local_allgatherv
    module procedure lsmpi_localallgatherv_realk4,lsmpi_localallgatherv_realk8
  end interface lsmpi_local_allgatherv

  interface lsmpi_win_create
    module procedure lsmpi_win_create_int4,lsmpi_win_create_int8,&
                   & lsmpi_win_create_realk8,lsmpi_win_create_realk4
  end interface lsmpi_win_create

  interface lsmpi_win_fence
    module procedure lsmpi_win_fence_simple,lsmpi_win_fence_special
  end interface lsmpi_win_fence

  interface lsmpi_put
    module procedure lsmpi_put_realk,&
                 &   lsmpi_put_realkV,lsmpi_put_realkV_wrapper8,lsmpi_put_realkV_parts,&
                 &   lsmpi_put_realkV_parts_wrapper8
  end interface lsmpi_put

  interface lsmpi_get
    module procedure lsmpi_get_realk,&
                 &   lsmpi_get_realkV,lsmpi_get_realkV_wrapper8,lsmpi_get_realkV_parts,&
                 &   lsmpi_get_realkV_parts_wrapper8,lsmpi_get_int4,lsmpi_get_int8
  end interface lsmpi_get

  interface lsmpi_acc
    module procedure lsmpi_acc_realk,&
                 &   lsmpi_acc_realkV,lsmpi_acc_realkV_wrapper8,lsmpi_acc_realkV_parts, &
                 &   lsmpi_acc_realkV_parts_wrapper8,lsmpi_acc_int4,lsmpi_acc_int8
  end interface lsmpi_acc

  interface lsmpi_get_acc
    module procedure lsmpi_get_acc_int444,lsmpi_get_acc_int888
  end interface lsmpi_get_acc
  !save

  


!!!!!!!!!!!!!!!!!!!!!!!!!
!Constants for MPIBUFFER!
!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef VAR_MPI
  integer,parameter     :: LSMPIBROADCAST       = 1
  integer,parameter     :: LSMPIREDUCTION       = 2
  integer,parameter     :: LSMPIREDUCTIONmaster = 3
  integer,parameter     :: LSMPISENDRECV        = 4


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!General MPI vars, aka junkbox!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=ls_mpik) :: MPI_COMM_LSDALTON
  logical               :: LSMPIASYNCP
  logical               :: lsmpi_enabled_comm_procs 

  !split mpi messages in case of 32bit mpi library to subparts, which are
  !describable by a 32bit integer and dividable by 8
  !integer,parameter     :: SPLIT_MPI_MSG = 2147483640
  integer,parameter     :: SPLIT_MPI_MSG      = 1000000000
  !The recommended size of message chunks
  integer,parameter     :: SPLIT_MSG_REC      =  100000000
  !split mpi one sided communication into 1GB msg, with CRAY workaround in 100MB
  !chunks
#ifndef VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN
  integer,parameter     :: MAX_SIZE_ONE_SIDED = 125000000
#else
  integer,parameter     :: MAX_SIZE_ONE_SIDED =  12500000
#endif

  !mpistatus
  integer(kind=ls_mpik) :: status(MPI_STATUS_SIZE) 

  type mpigroup
     integer(kind=ls_mpik)         :: groupsize
     integer(kind=ls_mpik),pointer :: ranks(:)
  end type mpigroup

#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!integer conversion factor!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef VAR_INT64
#ifdef VAR_MPI
#ifdef VAR_MPI_32BIT_INT
  integer,parameter :: int_to_short = 4 !int64,mpi & mpi32
#else
  integer,parameter :: int_to_short = 8 !int64,mpi nompi32
#endif
#else
  integer,parameter :: int_to_short = 8 !int64 nompi
#endif
#else
  integer,parameter :: int_to_short = 4 !no int64
#endif
  integer,parameter :: MaxIncreaseSize = 50000000 !0.4 GB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Checking and measuring variables!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  logical                     :: AddToBuffer
  integer(kind=long)          :: iLog,iDP,iInt4,iInt8,iSho,iCha
  integer(kind=long)          :: nLog,nDP,nShort,nInteger4,nInteger8,nCha
  real(realk),pointer         :: lsmpibufferDP(:)
  integer(kind=4),pointer     :: lsmpibufferInt4(:)
  integer(kind=8),pointer     :: lsmpibufferInt8(:)
  integer(kind=short),pointer :: lsmpibufferSho(:)
  logical,pointer             :: lsmpibufferLog(:)
  character,pointer           :: lsmpibufferCha(:)
  integer,parameter           :: incremLog=169,incremDP=100,incremInteger=626
  integer,parameter           :: incremCha=1510,incremShort=incremInteger*int_to_short
  real(realk)                 :: poketime=0.0E0_realk
  integer(kind=long)          :: poketimes = 0
  real(realk)                 :: time_win_unlock = 0.0E0_realk


!$OMP THREADPRIVATE(AddToBuffer,iLog,iDP,iInt4,iInt8,iSho,iCha,&
!$OMP nLog,nDP,nInteger4,nInteger8,nShort,nCha,lsmpibufferDP,lsmpibufferInt4,&
!$OMP lsmpibufferInt8,lsmpibufferSho,lsmpibufferLog,lsmpibufferCha)

contains
#ifdef VAR_MPI
    SUBROUTINE GET_MPI_COMM_SELF(outputcomm)
      implicit none
      integer(kind=ls_mpik),intent(inout) :: OutputComm   
      OutputComm = MPI_COMM_SELF !Predefined Communicator - the local processor
    END SUBROUTINE GET_MPI_COMM_SELF
#endif
    subroutine ls_mpibcast_integer(buffer,master,comm)
      implicit none
      integer(kind=4) :: buffer
      integer(kind=ls_mpik),intent(in) :: master
      integer(kind=ls_mpik),intent(in) :: comm   ! communicator
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: DATATYPE,n,IERR
      IERR=0
      DATATYPE = MPI_INTEGER4
      n = 1
      CALL MPI_BCAST(BUFFER,n,DATATYPE,master,comm,IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpibcast_integer

    subroutine ls_mpibcast_long(buffer,master,comm)
      implicit none
      integer(kind=8) :: buffer
      integer(kind=ls_mpik),intent(in) :: master
      integer(kind=ls_mpik),intent(in) :: comm   ! communicator
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: DATATYPE,n,IERR
      IERR=0
      DATATYPE = MPI_INTEGER8
      n = 1
      CALL MPI_BCAST(BUFFER,n,DATATYPE,master,comm,IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpibcast_long

    subroutine ls_mpibcast_short(buffer,master,comm)
      implicit none
      integer(kind=short) :: buffer
      integer(kind=ls_mpik) :: master
      integer(kind=ls_mpik) :: comm   ! communicator
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: ierr,n,datatype
      integer :: integerbuffer

      DATATYPE = MPI_INTEGER
      n = 1
      IERR=0

!     Convert from short integer to regular integer
      integerbuffer = buffer
!     Broadcast regular integer
      CALL MPI_BCAST(INTEGERBUFFER,n,DATATYPE,master,comm,IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
!     Convert back
      buffer        = integerbuffer
#endif
    end subroutine ls_mpibcast_short

    subroutine ls_mpibcast_shortV_wrapper8(buffer,n,master,comm)
      implicit none
      integer(kind=ls_mpik),intent(in) :: comm
      integer(kind=ls_mpik),intent(in) :: master
      integer(kind=8) :: n
      integer(kind=short) :: buffer(:)
#ifdef VAR_MPI
      integer(kind=4) :: n4
      integer(kind=8) :: i,k
      integer(kind=ls_mpik):: ierr,count,datatype
      integer(kind=4),pointer :: intbuffer(:)
      integer(kind=long) :: n1
      IERR=0

      !in case of a 32 bit mpi library, split the messages, else just send
      if(ls_mpik==4)then
        !loop over batches, which contain a number of elements,
        !describable by 32 bit integers, here 2E9
        k=SPLIT_MPI_MSG
        do i=1,n,k
          n4=k
          !if((n-i)<k)n4=mod(n,k)
          if(((n-i)<k).and.(mod(n-i+1,k)/=0))n4=mod(n,k)
          call ls_mpibcast_shortV(buffer(i:i+n4-1),n4,master,comm)
        enddo
      else

        DATATYPE = MPI_INTEGER
        COUNT = SIZE(buffer,kind=ls_mpik)
        n1 = SIZE(buffer,kind=long)
        IF(n1.NE.COUNT)call lsquit('lsmpi error1 in ls_mpibcast_ShortV',-1)
        if(COUNT.NE.n)call lsquit('lsmpi error2 in ls_mpibcast_ShortV',-1)

        IF (mod(n,int_to_short).NE.0) THEN
          write(*,*) 'Error in ls_mpibcast_shortV',n,int_to_short,mod(n,int_to_short)
          call lsquit('lsmpi nubf modular error in ls_mpibcast_ShortV',-1)
        ENDIF
        count = n/int_to_short
        CALL MPI_BCAST(BUFFER(1:n),count,DATATYPE,master,comm,IERR)
        IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
      endif
#endif
    end subroutine ls_mpibcast_shortV_wrapper8
    subroutine ls_mpibcast_shortV(buffer,nbuf,master,comm)
      implicit none
      integer(kind=short)  :: buffer(:)
      integer(kind=ls_mpik):: master
      integer(kind=4)      :: nbuf
      integer(kind=ls_mpik):: comm   ! communicator
#ifdef VAR_MPI
      integer(kind=ls_mpik):: ierr,count,datatype,n
      integer(kind=4),pointer :: intbuffer(:)
      integer(kind=long) :: n1
      IERR=0

      DATATYPE = MPI_INTEGER
      COUNT = SIZE(buffer,kind=ls_mpik)
      n1 = SIZE(buffer,kind=long)
      if(COUNT.NE.n1)call lsquit('lsmpi error1 in ls_mpibcast_ShortV',-1)
      if(COUNT.NE.nbuf)call lsquit('lsmpi error2 in ls_mpibcast_ShortV',-1)

      IF (mod(nbuf,int_to_short).NE.0) THEN
        write(*,*) 'Error in ls_mpibcast_shortV',nbuf,int_to_short,mod(nbuf,int_to_short)
        call lsquit('lsmpi nubf modular error in ls_mpibcast_ShortV',-1)
      ENDIF
      n = nbuf/int_to_short
      CALL MPI_BCAST(BUFFER(1:nbuf),n,DATATYPE,master,comm,IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpibcast_shortV


    subroutine ls_mpibcast_longV_wrapper8(buffer,n,master,comm)
      implicit none
      integer(kind=ls_mpik),intent(in) :: comm
      integer(kind=ls_mpik),intent(in) :: master
      integer(kind=8) :: n
      integer(kind=8) :: buffer(:)
#ifdef VAR_MPI
      integer(kind=4) :: n4
      integer(kind=8) :: i,k
      integer(kind=ls_mpik) :: ierr,count,datatype
      integer(kind=long) :: n1
      !loop over batches, which contain a number of elements,
      !describable by 32 bit integers, here 2E9
      IERR=0
      if(ls_mpik==4)then
        k=SPLIT_MPI_MSG
        do i=1,n,k
          n4=k
          !if((n-i)<k)n4=mod(n,k)
          if(((n-i)<k).and.(mod(n-i+1,k)/=0))n4=mod(n,k)
          call ls_mpibcast_longV(buffer(i:i+n4-1),n4,master,comm)
        enddo
      else
        DATATYPE = MPI_INTEGER8
        COUNT = SIZE(buffer,kind=ls_mpik)
        n1 = SIZE(buffer,kind=long)
        if(COUNT.NE.n1)call lsquit('lsmpi error1 in ls_mpibcast_longV_wrapper8',-1)

        if(COUNT.NE.n)call lsquit('lsmpi error in ls_mpibcast_longV_wrapper8',-1)
        !      COUNT = nbuf
        CALL MPI_BCAST(BUFFER,COUNT,DATATYPE,master,comm,IERR)
        IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
      endif
#endif
    end subroutine ls_mpibcast_longV_wrapper8
    subroutine ls_mpibcast_longV(buffer,nbuf,master,comm)
      implicit none
      integer(kind=8) :: buffer(:)
      integer(kind=4) :: nbuf
      integer(kind=ls_mpik) :: master
      integer(kind=ls_mpik) :: comm   ! communicator
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: ierr,datatype,n
      integer(kind=long) :: n1
      IERR=0
      DATATYPE = MPI_INTEGER8
      n = SIZE(buffer)
      n1 = SIZE(buffer,kind=long)
      if(n.NE.n1)call lsquit('lsmpi error1 in ls_mpibcast_longV',-1)

      if(n.NE.nbuf)call lsquit('lsmpi error in ls_mpibcast_longV',-1)
      !      COUNT = nbuf
      CALL MPI_BCAST(BUFFER,n,DATATYPE,master,comm,IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpibcast_longV


    subroutine ls_mpibcast_integerV_wrapper8(buffer,n,master,comm)
      implicit none
      integer(kind=ls_mpik),intent(in) :: comm
      integer(kind=ls_mpik),intent(in) :: master
      integer(kind=8) :: n
      integer(kind=4) :: buffer(:)
#ifdef VAR_MPI
      integer(kind=4) :: n4
      integer(kind=8) :: i,k
      !loop over batches, which contain a number of elements,
      !describable by 32 bit integers, here 2E9
      integer(kind=ls_mpik) :: ierr,count,datatype
      integer(kind=long) :: n1
      IERR=0
      if(ls_mpik==4)then
        k=SPLIT_MPI_MSG
        do i=1,n,k
          n4=k
          !if((n-i)<k)n4=mod(n,k)
          if(((n-i)<k).and.(mod(n-i+1,k)/=0))n4=mod(n,k)
          call ls_mpibcast_integerV(buffer(i:i+n4-1),n4,master,comm)
        enddo
      else
        DATATYPE = MPI_INTEGER4
        COUNT = SIZE(buffer,kind=ls_mpik)
        n1 = SIZE(buffer,kind=long)
        if(COUNT.NE.n1)call lsquit('lsmpi error1 in ls_mpibcast_longV',-1)
        if(COUNT.NE.n)call lsquit('lsmpi error in ls_mpibcast_integerV',-1)
        CALL MPI_BCAST(BUFFER,COUNT,DATATYPE,master,comm,IERR)
        IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
      endif
#endif
    end subroutine ls_mpibcast_integerV_wrapper8
    subroutine ls_mpibcast_integerV(buffer,nbuf,master,comm)
      implicit none
      integer(kind=4)       :: nbuf
      integer(kind=ls_mpik) :: master
      integer(kind=4)       :: buffer(:)
      integer(kind=ls_mpik) :: comm   ! communicator
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: ierr,count,datatype
      integer(kind=long) :: n1
      IERR=0
      DATATYPE = MPI_INTEGER4
      COUNT = SIZE(buffer,kind=ls_mpik)
      n1 = SIZE(buffer,kind=long)
      if(COUNT.NE.n1)call lsquit('lsmpi error1 in ls_mpibcast_integerV',-1)

      if(COUNT.NE.nbuf)call lsquit('lsmpi error in ls_mpibcast_integerV',-1)
      CALL MPI_BCAST(BUFFER,COUNT,DATATYPE,master,comm,IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpibcast_integerV



    subroutine ls_mpibcast_longM_wrapper8(buffer,n1,n2,master,comm)
      implicit none
      integer(kind=ls_mpik),intent(in) :: comm
      integer(kind=ls_mpik),intent(in) :: master
      integer(kind=8) :: n1,n2
      integer(kind=8) :: buffer(:,:)
#ifdef VAR_MPI
      integer(kind=8),pointer :: buffertmp(:)
      integer(kind=4) :: n4
      integer(kind=8) :: i,k,n
      integer(kind=ls_mpik) :: ierr,count,datatype
      !loop over batches, which contain a number of elements,
      !describable by 32 bit integers, here 2E9
      IERR=0
      n=n1*n2
      if(ls_mpik==4)then
         call ass_88I2to1(buffer,buffertmp,[n1,n2])
        k=SPLIT_MPI_MSG
        do i=1,n,k
          n4=k
          !if((n-i)<k)n4=mod(n,k)
          if(((n-i)<k).and.(mod(n-i+1,k)/=0))n4=mod(n,k)
          call ls_mpibcast_longV(buffertmp(i:i+n4-1),n4,master,comm)
        enddo
        nullify(buffertmp)
      else
        DATATYPE = MPI_INTEGER8
        COUNT = SIZE(buffer,kind=ls_mpik)
        n1 = SIZE(buffer,kind=long)
        if(COUNT.NE.n1)call lsquit('lsmpi error1 in ls_mpibcast_longM_wrapper8',-1)
        if(COUNT.NE.n)call lsquit('lsmpi error in ls_mpibcast_longM_wrapper8',-1)
        !      COUNT = nbuf
        CALL MPI_BCAST(BUFFER,COUNT,DATATYPE,master,comm,IERR)
        IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
      endif
#endif
    end subroutine ls_mpibcast_longM_wrapper8
    subroutine ls_mpibcast_longM(buffer,nbuf1,nbuf2,master,comm)
      implicit none
      integer(kind=8) :: buffer(:,:)
      integer(kind=4) :: nbuf1,nbuf2
      integer(kind=ls_mpik) :: master
      integer(kind=ls_mpik) :: comm   ! communicator
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: ierr,datatype,n
      integer(kind=long)    :: n1
      IERR=0
      DATATYPE = MPI_INTEGER8
      n = SIZE(buffer,kind=ls_mpik)
      n1 = SIZE(buffer,kind=long)
      if(n.NE.n1)call lsquit('lsmpi error1 in ls_mpibcast_longM',-1)

      if(n.NE.nbuf1*nbuf2)call lsquit('lsmpi error in ls_mpibcast_longM',-1)
      !      COUNT = nbuf
      CALL MPI_BCAST(BUFFER,n,DATATYPE,master,comm,IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpibcast_longM


    subroutine ls_mpibcast_integerM_wrapper8(buffer,n1,n2,master,comm)
      implicit none
      integer(kind=ls_mpik),intent(in) :: comm
      integer(kind=ls_mpik),intent(in) :: master
      integer(kind=8) :: n1,n2
      integer(kind=4) :: buffer(:,:)
#ifdef VAR_MPI
      integer(kind=4),pointer :: buffertmp(:)
      integer(kind=8) :: i,k,n
      integer(kind=4) :: n4
      integer(kind=long)    :: n1T
      !loop over batches, which contain a number of elements,
      !describable by 32 bit integers, here 2E9
      integer(kind=ls_mpik) :: ierr,count,datatype
      IERR=0
      n=n1*n2
      if(ls_mpik==4)then
         call ass_48I2to1(buffer,buffertmp,[n1,n2])
        k=SPLIT_MPI_MSG
        do i=1,n,k
          n4=k
          if(((n-i)<k).and.(mod(n-i+1,k)/=0))n4=mod(n,k)
          call ls_mpibcast_integerV(buffertmp(i:i+n4-1),n4,master,comm)
        enddo
        nullify(buffertmp)
      else
        DATATYPE = MPI_INTEGER4
        COUNT = SIZE(buffer,kind=ls_mpik)
        n1T = SIZE(buffer,kind=long)
        if(COUNT.NE.n1T)call lsquit('lsmpi error1 in ls_mpibcast_integerM_wrapper8',-1)
        if(COUNT.NE.n)call lsquit('lsmpi error in ls_mpibcast_integerM_wrapper8',-1)
        CALL MPI_BCAST(BUFFER,COUNT,DATATYPE,master,comm,IERR)
        IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
      endif
#endif
    end subroutine ls_mpibcast_integerM_wrapper8

    subroutine ls_mpibcast_integerM(buffer,nbuf1,nbuf2,master,comm)
      implicit none
      integer(kind=4)       :: nbuf1,nbuf2
      integer(kind=ls_mpik) :: master
      integer(kind=4)       :: buffer(:,:)
      integer(kind=ls_mpik) :: comm   ! communicator
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: ierr,count,datatype
      integer(kind=long)    :: n1
      IERR=0
      DATATYPE = MPI_INTEGER4
      COUNT = SIZE(buffer,kind=ls_mpik)
      n1 = SIZE(buffer,kind=long)
      if(COUNT.NE.n1)call lsquit('lsmpi error1 in ls_mpibcast_integerM',-1)

      if(COUNT.NE.nbuf1*nbuf2)call lsquit('lsmpi error in ls_mpibcast_integerM',-1)
      CALL MPI_BCAST(BUFFER,COUNT,DATATYPE,master,comm,IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpibcast_integerM



    subroutine ls_mpibcast_realk(buffer,master,comm)
      implicit none
      real(realk) :: buffer
      integer(kind=ls_mpik) :: master
      integer(kind=ls_mpik) :: comm   ! communicator
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: ierr,count,datatype
      IERR=0
      DATATYPE = MPI_DOUBLE_PRECISION
      COUNT = 1
      CALL MPI_BCAST(BUFFER,COUNT,DATATYPE,master,comm,IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpibcast_realk

    subroutine ls_mpibcast_realkV_wrapper8(buffer,n,master,comm)
      implicit none
      integer(kind=ls_mpik) :: master
      integer(kind=ls_mpik) :: comm   ! communicator
      integer(kind=8) :: n
      real(realk) :: buffer(:)
#ifdef VAR_MPI
      integer(kind=4) :: n4
      integer(kind=8) :: i,k
      integer(kind=long) :: n1
      integer(kind=ls_mpik) :: ierr,count,datatype      
      !loop over batches, which contain a number of elements,
      !describable by 32 bit integers, here 2E9
      IERR=0
      if(ls_mpik==4)then
        k=SPLIT_MPI_MSG
        do i=1,n,k
          n4=k
          if(((n-i)<k).and.(mod(n-i+1,k)/=0))n4=mod(n,k)
          call ls_mpibcast_realkV(buffer(i:i+n4-1),n4,master,comm)
        enddo
      else
        DATATYPE = MPI_DOUBLE_PRECISION
        COUNT = SIZE(buffer,kind=ls_mpik)
        n1 = SIZE(buffer,kind=long)
        if(COUNT.NE.n1)call lsquit('lsmpi error1 in ls_mpibcast_realkV_wrapper8',-1)
        if(COUNT.NE.n)THEN
           call lsquit('lsmpi error in ls_mpibcast_realkV_wrapper8',-1)
        endif
        CALL MPI_BCAST(BUFFER,COUNT,DATATYPE,master,comm,IERR)
        IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
      endif
#endif
    end subroutine ls_mpibcast_realkV_wrapper8
    subroutine ls_mpibcast_realkV(buffer,nbuf,master,comm)
      implicit none
      integer(kind=ls_mpik) :: master
      integer(kind=ls_mpik) :: comm   ! communicator
      integer(kind=4) :: nbuf
      real(realk) :: buffer(:)
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: ierr,count,datatype
      integer(kind=long)    :: n1
      IERR=0
      DATATYPE = MPI_DOUBLE_PRECISION
      COUNT = SIZE(buffer,kind=ls_mpik)
      n1 = SIZE(buffer,kind=long)
      if(COUNT.NE.n1)call lsquit('lsmpi error1 in ls_mpibcast_realkV',-1)
      if(COUNT.NE.nbuf)THEN
         call lsquit('lsmpi error in ls_mpibcast_realkV',-1)
      endif
      COUNT = nbuf
      CALL MPI_BCAST(BUFFER,COUNT,DATATYPE,master,comm,IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpibcast_realkV

    subroutine ls_mpibcast_realkM(buffer,nbuf1,nbuf2,master,comm)
      implicit none
      integer :: nbuf1,nbuf2
      integer(kind=ls_mpik) :: master
      real(realk) :: buffer(:,:)
      integer(kind=ls_mpik) :: comm   ! communicator
#ifdef VAR_MPI
      real(realk),pointer :: buf(:)
      call ass_8D2to1(buffer,buf,[i8*nbuf1,i8*nbuf2])
      call ls_mpibcast_realkV_wrapper8(buf,(i8*nbuf1)*nbuf2,master,comm)
      buf => null()
#endif
    end subroutine ls_mpibcast_realkM

    subroutine ls_mpibcast_realkT(buffer,nbuf1,nbuf2,nbuf3,master,comm)
      implicit none
      integer(kind=ls_mpik) :: master
      integer :: nbuf1,nbuf2,nbuf3
      real(realk) :: buffer(:,:,:)
      integer(kind=ls_mpik) :: comm   ! communicator
      real(realk),pointer :: buf(:)
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: ierr,count,datatype
      call ass_8D3to1(buffer,buf,[i8*nbuf1,i8*nbuf2,i8*nbuf3])
      call ls_mpibcast_realkV_wrapper8(buf,((i8*nbuf1)*nbuf2)*nbuf3,master,comm)
      buf => null()
#endif
    end subroutine ls_mpibcast_realkT

    subroutine ls_mpibcast_realkQ(buffer,nbuf1,nbuf2,nbuf3,nbuf4,master,comm)
      implicit none
      integer(kind=ls_mpik) :: master
      integer :: nbuf1,nbuf2,nbuf3,nbuf4
      real(realk) :: buffer(:,:,:,:)
      integer(kind=ls_mpik) :: comm   ! communicator
      real(realk),pointer :: buf(:)
#ifdef VAR_MPI
      call ass_8D4to1(buffer,buf,[i8*nbuf1,i8*nbuf2,i8*nbuf3,i8*nbuf4])
      call ls_mpibcast_realkV_wrapper8(buf,(((i8*nbuf1)*nbuf2)*nbuf3)*nbuf4,master,comm)
      buf => null()
#endif
    end subroutine ls_mpibcast_realkQ

    subroutine ls_mpibcast_logical8(buffer,master,comm)
      implicit none
      logical(kind=8),intent(inout) :: buffer
      integer(kind=ls_mpik),intent(in) :: master
      integer(kind=ls_mpik),intent(in) :: comm   ! communicator
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: ierr,count,datatype
      integer(kind=MPI_ADDRESS_KIND) :: mpi_logical_extent,lb
      logical(kind=4) :: buffer4
      DATATYPE = MPI_LOGICAL
      COUNT = 1
      IERR=0
      call MPI_TYPE_GET_EXTENT(MPI_LOGICAL,lb,mpi_logical_extent,ierr)
      IF(mpi_logical_extent.EQ.4)THEN
         !32 bit mpi logical
         BUFFER4 = BUFFER
         CALL MPI_BCAST(BUFFER4,COUNT,DATATYPE,master,comm,IERR)         
         BUFFER = BUFFER4
      ELSE
         !64 bit mpi logical
         CALL MPI_BCAST(BUFFER,COUNT,DATATYPE,master,comm,IERR)         
      ENDIF
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpibcast_logical8

    subroutine ls_mpibcast_logical4(buffer,master,comm)
      implicit none
      logical(kind=4),intent(inout) :: buffer
      integer(kind=ls_mpik),intent(in) :: master
      integer(kind=ls_mpik),intent(in) :: comm   ! communicator
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: ierr,count,datatype
      logical(kind=8) :: buffer8
      integer(kind=MPI_ADDRESS_KIND) :: mpi_logical_extent,lb
      IERR=0
      DATATYPE = MPI_LOGICAL
      COUNT = 1
      call MPI_TYPE_GET_EXTENT(MPI_LOGICAL,lb,mpi_logical_extent,ierr)
      IF(mpi_logical_extent.EQ.4)THEN
         !32 bit mpi logical
         CALL MPI_BCAST(BUFFER,COUNT,DATATYPE,master,comm,IERR)         
      ELSE
         !64 bit mpi logical
         BUFFER8 = BUFFER
         CALL MPI_BCAST(BUFFER8,COUNT,DATATYPE,master,comm,IERR)         
         BUFFER = BUFFER8
      ENDIF
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpibcast_logical4

    subroutine ls_mpibcast_logical8V_wrapper8(buffer,n,master,comm)
      implicit none
      integer(kind=ls_mpik),intent(in) :: master
      integer(kind=ls_mpik),intent(in) :: comm   ! communicator
      logical(kind=8),intent(inout) :: buffer(:)
      integer(kind=8),intent(in) :: n
#ifdef VAR_MPI
      integer(kind=4) :: n4
      integer(kind=8) :: i,k
      integer(kind=ls_mpik) :: ierr,count,datatype
      logical(kind=4),pointer :: buffer4(:)
      integer(kind=MPI_ADDRESS_KIND) :: mpi_logical_extent,lb
      integer(kind=long) :: n1
      !loop over batches, which contain a number of elements,
      !escribable by 32 bit integers, here 2E9
      IERR=0
      if(ls_mpik==4)then
        k=SPLIT_MPI_MSG
        do i=1,n,k
          n4=k
          if(((n-i)<k).and.(mod(n-i+1,k)/=0))n4=mod(n,k)
          call ls_mpibcast_logical8V(buffer(i:i+n4-1),n4,master,comm)
        enddo
      else
        DATATYPE = MPI_LOGICAL
        COUNT = SIZE(buffer,kind=ls_mpik)
        n1 = SIZE(buffer,kind=long)
        if(COUNT.NE.n1)call lsquit('lsmpi error1 in ls_mpibcast_logical8V_wrapper8',-1)
        if(COUNT.NE.n)THEN
           call lsquit('lsmpi error in ls_mpibcast_logical8V_wrapper8',-1)
        endif
        call MPI_TYPE_GET_EXTENT(MPI_LOGICAL,lb,mpi_logical_extent,ierr)
        IF(mpi_logical_extent.EQ.4)THEN
           !32 bit mpi logical
!           n4=n
           call mem_alloc(BUFFER4,n)
           DO I=1,n
              BUFFER4(I) = BUFFER(I)
           ENDDO
           CALL MPI_BCAST(BUFFER4,COUNT,DATATYPE,master,comm,IERR)         
           DO I=1,n
              BUFFER(I) = BUFFER4(I)
           ENDDO
           call mem_dealloc(BUFFER4)
        ELSE
           !64 bit mpi logical
           CALL MPI_BCAST(BUFFER,COUNT,DATATYPE,master,comm,IERR)         
        ENDIF
        IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
      endif
#endif
    end subroutine ls_mpibcast_logical8V_wrapper8
    subroutine ls_mpibcast_logical8V(buffer,nbuf,master,comm)
      implicit none
      integer(kind=4),intent(in) :: nbuf
      integer(kind=ls_mpik),intent(in) :: master
      integer(kind=ls_mpik),intent(in) :: comm   ! communicator
      logical(kind=8),intent(inout) :: buffer(:)
#ifdef VAR_MPI
      integer :: I
      logical(kind=4),pointer :: buffer4(:)
      integer(kind=ls_mpik) :: ierr,count,datatype
      integer(kind=MPI_ADDRESS_KIND) :: mpi_logical_extent,lb
      integer(kind=long) :: n1
      IERR=0
      DATATYPE = MPI_LOGICAL
      COUNT = SIZE(buffer,kind=ls_mpik)
      n1 = SIZE(buffer,kind=long)
      if(COUNT.NE.n1)call lsquit('lsmpi error1 in ls_mpibcast_logical8V',-1)
      if(COUNT.NE.nbuf)THEN
         call lsquit('lsmpi error in ls_mpibcast_logical8V',-1)
      endif
      call MPI_TYPE_GET_EXTENT(MPI_LOGICAL,lb,mpi_logical_extent,ierr)
      IF(mpi_logical_extent.EQ.4)THEN
         !32 bit mpi logical
         call mem_alloc(BUFFER4,nbuf)
         DO I=1,nbuf
            BUFFER4(I) = BUFFER(I)
         ENDDO
         CALL MPI_BCAST(BUFFER4,COUNT,DATATYPE,master,comm,IERR)         
         DO I=1,nbuf
            BUFFER(I) = BUFFER4(I)
         ENDDO
         call mem_dealloc(BUFFER4)
      ELSE
         !64 bit mpi logical
         CALL MPI_BCAST(BUFFER,COUNT,DATATYPE,master,comm,IERR)         
      ENDIF
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpibcast_logical8V

    subroutine ls_mpibcast_logical4V_wrapper8(buffer,n,master,comm)
      implicit none
      integer(kind=ls_mpik),intent(in) :: master
      integer(kind=ls_mpik),intent(in) :: comm   ! communicator
      logical(kind=4),intent(inout) :: buffer(:)
      integer(kind=8),intent(in) :: n
#ifdef VAR_MPI
      logical(kind=8),pointer :: buffer8(:)
      integer(kind=4) :: n4
      integer(kind=8) :: i,k
      integer(kind=ls_mpik) :: ierr,count,datatype
      integer(kind=MPI_ADDRESS_KIND) :: mpi_logical_extent,lb
      integer(kind=long) :: n1
      !loop over batches, which contain a number of elements,
      !escribable by 32 bit integers, here 2E9
      IERR=0
      if(ls_mpik==4)then
        k=SPLIT_MPI_MSG
        do i=1,n,k
          n4=k
          if(((n-i)<k).and.(mod(n-i+1,k)/=0))n4=mod(n,k)
          call ls_mpibcast_logical4V(buffer(i:i+n4-1),n4,master,comm)
        enddo
      else
        DATATYPE = MPI_LOGICAL
        COUNT = SIZE(buffer,kind=ls_mpik)
        n1 = SIZE(buffer,kind=long)
        if(COUNT.NE.n1)call lsquit('lsmpi error1 in ls_mpibcast_logicalV',-1)
        if(COUNT.NE.n)THEN
           call lsquit('lsmpi error in ls_mpibcast_logicalV',-1)
        endif
        call MPI_TYPE_GET_EXTENT(MPI_LOGICAL,lb,mpi_logical_extent,ierr)
        IF(mpi_logical_extent.EQ.4)THEN
           !32 bit mpi logical
           CALL MPI_BCAST(BUFFER,COUNT,DATATYPE,master,comm,IERR)         
        ELSE
           !64 bit mpi logical
           call mem_alloc(BUFFER8,n)
           DO I=1,n
              BUFFER8(I) = BUFFER(I)
           ENDDO
           CALL MPI_BCAST(BUFFER8,COUNT,DATATYPE,master,comm,IERR)         
           DO I=1,n
              BUFFER(I) = BUFFER8(I)
           ENDDO
           call mem_dealloc(BUFFER8)
        ENDIF
        IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
      endif
#endif
    end subroutine ls_mpibcast_logical4V_wrapper8
    subroutine ls_mpibcast_logical4V(buffer,nbuf,master,comm)
      implicit none
      integer(kind=4),intent(in) :: nbuf
      integer(kind=ls_mpik),intent(in) :: master
      integer(kind=ls_mpik),intent(in) :: comm   ! communicator
      logical(kind=4),intent(inout) :: buffer(:)
#ifdef VAR_MPI
      integer(kind=8) :: nbuf8
      logical(kind=8),pointer :: buffer8(:)
      integer(kind=ls_mpik) :: ierr,count,datatype
      integer(kind=MPI_ADDRESS_KIND) :: mpi_logical_extent,lb
      integer(kind=long) :: n1
      integer :: I
      IERR=0
      DATATYPE = MPI_LOGICAL
      COUNT = SIZE(buffer,kind=ls_mpik)
      n1 = SIZE(buffer,kind=long)
      if(COUNT.NE.n1)call lsquit('lsmpi error1 in ls_mpibcast_logicalV',-1)

      if(COUNT.NE.nbuf)THEN
         call lsquit('lsmpi error in ls_mpibcast_logicalV',-1)
      endif
      call MPI_TYPE_GET_EXTENT(MPI_LOGICAL,lb,mpi_logical_extent,ierr)
      IF(mpi_logical_extent.EQ.4)THEN
         !32 bit mpi logical
         CALL MPI_BCAST(BUFFER,COUNT,DATATYPE,master,comm,IERR)         
      ELSE
         !64 bit mpi logical
         nbuf8 = nbuf
         call mem_alloc(BUFFER8,nbuf8)
         DO I=1,nbuf8
            BUFFER8(I) = BUFFER(I)
         ENDDO
         CALL MPI_BCAST(BUFFER8,COUNT,DATATYPE,master,comm,IERR)         
         DO I=1,nbuf8
            BUFFER(I) = BUFFER8(I)
         ENDDO
         call mem_dealloc(BUFFER8)
      ENDIF
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpibcast_logical4V

    subroutine ls_mpibcast_logical4M4(buffer,nbuf1,nbuf2,master,comm)
      implicit none
      integer(kind=4),intent(in) :: nbuf1,nbuf2
      integer(kind=ls_mpik),intent(in) :: master
      logical(kind=4),intent(inout) :: buffer(:,:)
      integer(kind=ls_mpik),intent(in) :: comm   ! communicator
#ifdef VAR_MPI
      logical(kind=4),pointer :: buf(:)
      integer(kind=8) :: n
      n=nbuf1*nbuf2
      call ass_44L2to1(buffer,buf,[nbuf1,nbuf2])
      call ls_mpibcast_logical4V_wrapper8(buf,n,master,comm)
      nullify(buf)
#endif
    end subroutine ls_mpibcast_logical4M4

    subroutine ls_mpibcast_logical4M8(buffer,nbuf1,nbuf2,master,comm)
      implicit none
      integer(kind=8),intent(in) :: nbuf1,nbuf2
      integer(kind=ls_mpik),intent(in) :: master
      logical(kind=4),intent(inout) :: buffer(:,:)
      integer(kind=ls_mpik),intent(in) :: comm   ! communicator
#ifdef VAR_MPI
      logical(kind=4),pointer :: buf(:)
      integer(kind=8) :: n
      n=nbuf1*nbuf2
      call ass_48L2to1(buffer,buf,[nbuf1,nbuf2])
      call ls_mpibcast_logical4V_wrapper8(buf,n,master,comm)
      nullify(buf)
#endif
    end subroutine ls_mpibcast_logical4M8
    subroutine ls_mpibcast_logical8M4(buffer,nbuf1,nbuf2,master,comm)
      implicit none
      integer(kind=4),intent(in) :: nbuf1,nbuf2
      integer(kind=ls_mpik),intent(in) :: master
      logical(kind=8),intent(inout) :: buffer(:,:)
      integer(kind=ls_mpik),intent(in) :: comm   ! communicator
#ifdef VAR_MPI
      logical(kind=8),pointer :: buf(:)
      integer(kind=8) :: n
      n=nbuf1*nbuf2
      call ass_84L2to1(buffer,buf,[nbuf1,nbuf2])
      call ls_mpibcast_logical8V_wrapper8(buf,n,master,comm)
      nullify(buf)
#endif
    end subroutine ls_mpibcast_logical8M4
    subroutine ls_mpibcast_logical8M8(buffer,nbuf1,nbuf2,master,comm)
      implicit none
      integer(kind=8),intent(in) :: nbuf1,nbuf2
      integer(kind=ls_mpik),intent(in) :: master
      logical(kind=8),intent(inout) :: buffer(:,:)
      integer(kind=ls_mpik),intent(in) :: comm   ! communicator
#ifdef VAR_MPI
      logical(kind=8),pointer :: buf(:)
      integer(kind=8) :: n
      n=nbuf1*nbuf2
      call ass_88L2to1(buffer,buf,[nbuf1,nbuf2])
      call ls_mpibcast_logical8V_wrapper8(buf,n,master,comm)
      nullify(buf)
#endif
    end subroutine ls_mpibcast_logical8M8

    subroutine ls_mpibcast_charac(buffer,master,comm)
      implicit none
      character :: buffer
      integer(kind=ls_mpik) :: master
      integer(kind=ls_mpik) :: comm   ! communicator
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: ierr,count,datatype
      IERR=0
      DATATYPE = MPI_CHARACTER
      COUNT = 1
      CALL MPI_BCAST(BUFFER,COUNT,DATATYPE,master,comm,IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpibcast_charac

  subroutine ls_mpibcast_characV_wrapper8(buffer,n,master,comm)
    implicit none
    integer(kind=ls_mpik),intent(in) :: comm
    integer(kind=ls_mpik),intent(in) :: master
    integer(kind=8) :: n
    character*(*) :: buffer
#ifdef VAR_MPI
    integer(kind=4) :: n4
    integer(kind=8) :: i,k
    integer(kind=ls_mpik) :: ierr,count,datatype
    !loop over batches, which contain a number of elements,
    !describable by 32 bit integers, here 2E9
    IERR=0
    if(ls_mpik==4)then
      k=SPLIT_MPI_MSG
      do i=1,n,k
        n4=k
        if(((n-i)<k).and.(mod(n-i+1,k)/=0))n4=mod(n,k)
        call ls_mpibcast_characV(buffer(i:i+n4-1),n4,master,comm)
      enddo
    else
      DATATYPE = MPI_CHARACTER
      COUNT = n
      CALL MPI_BCAST(BUFFER,COUNT,DATATYPE,master,comm,IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
    endif
#endif
    end subroutine ls_mpibcast_characV_wrapper8
    subroutine ls_mpibcast_characV(buffer,nbuf,master,comm)
      implicit none
      integer(kind=4) :: nbuf
      integer(kind=ls_mpik),intent(in) :: master
      character*(*) :: buffer
      integer(kind=ls_mpik),intent(in) :: comm   ! communicator
#ifdef VAR_MPI
      character :: tmpbuffer
      integer(kind=ls_mpik) :: ierr,count,datatype
      IERR=0
      DATATYPE = MPI_CHARACTER
      COUNT = nbuf
      CALL MPI_BCAST(BUFFER,COUNT,DATATYPE,master,comm,IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpibcast_characV

  subroutine ls_mpibcast_characV2_wrapper8(buffer,n,master,comm)
    implicit none
    integer(kind=ls_mpik),intent(in) :: comm
    integer(kind=ls_mpik),intent(in) :: master
    integer(kind=8) :: n
    character :: buffer(n)
#ifdef VAR_MPI
    integer(kind=4) :: n4
    integer(kind=8) :: i,k
    integer(kind=ls_mpik) :: ierr,count,datatype
    !loop over batches, which contain a number of elements,
    !describable by 32 bit integers, here 2E9
      IERR=0
    if(ls_mpik==4)then
      k=SPLIT_MPI_MSG
      do i=1,n,k
        n4=k
        if(((n-i)<k).and.(mod(n-i+1,k)/=0))n4=mod(n,k)
        call ls_mpibcast_characV2(buffer(i:i+n4-1),n4,master,comm)
      enddo
    else
      DATATYPE = MPI_CHARACTER
      COUNT = n
      CALL MPI_BCAST(BUFFER,COUNT,DATATYPE,master,comm,IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
    endif
#endif
    end subroutine ls_mpibcast_characV2_wrapper8

    subroutine ls_mpibcast_characV2(buffer,nbuf,master,comm)
      implicit none
      integer(kind=4),intent(in) :: nbuf
      character,intent(inout) :: buffer(:)
      integer(kind=ls_mpik),intent(in) :: master
      integer(kind=ls_mpik),intent(in) :: comm   ! communicator
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: ierr,count,datatype,i
      IERR=0
      DATATYPE = MPI_CHARACTER
      COUNT = nbuf
      CALL MPI_BCAST(BUFFER,COUNT,DATATYPE,master,comm,IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpibcast_characV2

    subroutine ls_mpibcast_realkV_parts44(buffer,nelms,master,comm,batchsze)
      implicit none
      integer(kind=ls_mpik) :: master
      integer(kind=ls_mpik) :: comm   ! communicator
      integer(kind=4),intent(in) :: batchsze
      integer(kind=4),intent(in) :: nelms
      real(realk) :: buffer(:)
#ifdef VAR_MPI
      integer(kind=8) :: i,n
      do i=1,nelms,batchsze
        n=batchsze
        if(((nelms-i)<batchsze).and.&
          &(mod(nelms-i+1,batchsze)/=0))n=mod(nelms,batchsze)
        call ls_mpibcast_realkV_wrapper8(buffer(i:i+n-1),n,master,comm)
      enddo
#endif
    end subroutine ls_mpibcast_realkV_parts44
    subroutine ls_mpibcast_realkV_parts48(buffer,nelms,master,comm,batchsze)
      implicit none
      integer(kind=ls_mpik) :: master
      integer(kind=ls_mpik) :: comm   ! communicator
      integer(kind=4),intent(in) :: batchsze
      integer(kind=8),intent(in) :: nelms
      real(realk) :: buffer(:)
#ifdef VAR_MPI
      integer(kind=8) :: i,n
      do i=1,nelms,batchsze
        n=batchsze
        if(((nelms-i)<batchsze).and.&
          &(mod(nelms-i+1,batchsze)/=0))n=mod(nelms,batchsze)
        call ls_mpibcast_realkV_wrapper8(buffer(i:i+n-1),n,master,comm)
      enddo
#endif
    end subroutine ls_mpibcast_realkV_parts48
    subroutine ls_mpibcast_realkV_parts84(buffer,nelms,master,comm,batchsze)
      implicit none
      integer(kind=ls_mpik) :: master
      integer(kind=ls_mpik) :: comm   ! communicator
      integer(kind=8),intent(in) :: batchsze
      integer(kind=4),intent(in) :: nelms
      real(realk) :: buffer(:)
#ifdef VAR_MPI
      integer(kind=8) :: i,n
      do i=1,nelms,batchsze
        n=batchsze
        if(((nelms-i)<batchsze).and.&
          &(mod(nelms-i+1,batchsze)/=0))n=mod(nelms,batchsze)
        call ls_mpibcast_realkV_wrapper8(buffer(i:i+n-1),n,master,comm)
      enddo
#endif
    end subroutine ls_mpibcast_realkV_parts84
    subroutine ls_mpibcast_realkV_parts88(buffer,nelms,master,comm,batchsze)
      implicit none
      integer(kind=ls_mpik) :: master
      integer(kind=ls_mpik) :: comm   ! communicator
      integer(kind=8),intent(in) :: batchsze
      integer(kind=8),intent(in) :: nelms
      real(realk) :: buffer(:)
#ifdef VAR_MPI
      integer(kind=8) :: i,n
      do i=1,nelms,batchsze
        n=batchsze
        if(((nelms-i)<batchsze).and.&
          &(mod(nelms-i+1,batchsze)/=0))n=mod(nelms,batchsze)
        call ls_mpibcast_realkV_wrapper8(buffer(i:i+n-1),n,master,comm)
      enddo
#endif
    end subroutine ls_mpibcast_realkV_parts88
! ########################################################################
!                          MPI SEND/RECEIVE
! ########################################################################
    subroutine lsmpi_send_realkV_8(buffer,nbuf,comm,receiver)
      implicit none
      real(realk) :: buffer(:)
      integer(kind=8) :: nbuf
      integer(kind=ls_mpik) :: comm
      integer(kind=ls_mpik) :: receiver
      integer(kind=ls_mpik) :: tag
      integer(kind=ls_mpik) :: nel,dtype,ierr
      integer(kind=4) :: n4
      integer(kind=8) :: i,k
#ifdef VAR_MPI
      tag   = 124_ls_mpik
      ierr  = 0_ls_mpik
      nel   = int(nbuf,kind=ls_mpik)
      dtype = MPI_DOUBLE_PRECISION
      if(ls_mpik==4)then
        k=SPLIT_MPI_MSG
        do i=1,nbuf,k
          n4=k
          if(((nbuf-i)<k).and.(mod(nbuf-i+1,k)/=0))n4=mod(nbuf,k)
          call lsmpi_send_realkV_4(buffer(i:i+n4-1),n4,comm,receiver)
        enddo
      else
        call MPI_SEND(buffer,nbuf,dtype,receiver,tag,comm,ierr)
        IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
      endif
#endif 
    end subroutine lsmpi_send_realkV_8
    subroutine lsmpi_send_realkV_4(buffer,nbuf,comm,receiver)
      implicit none
      real(realk) :: buffer(:)
      integer(kind=4) :: nbuf
      integer(kind=ls_mpik) :: comm
      integer(kind=ls_mpik) :: receiver
      integer(kind=ls_mpik) :: tag
      integer(kind=ls_mpik) :: nel,dtype,ierr
#ifdef VAR_MPI
      tag   = 124_ls_mpik
      ierr  = 0_ls_mpik
      nel   = int(nbuf,kind=ls_mpik)
      dtype = MPI_DOUBLE_PRECISION
      call MPI_SEND(buffer,nbuf,dtype,receiver,tag,comm,ierr)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif 
    end subroutine lsmpi_send_realkV_4
    subroutine lsmpi_recv_realkV_8(buffer,nbuf,comm,sender)
      implicit none
      real(realk) :: buffer(:)
      integer(kind=8) :: nbuf
      integer(kind=ls_mpik) :: comm
      integer(kind=ls_mpik) :: sender
      integer(kind=ls_mpik) :: tag
      integer(kind=ls_mpik) :: nel,dtype,ierr
      integer(kind=4) :: n4
      integer(kind=8) :: i,k
#ifdef VAR_MPI
      tag   = 124_ls_mpik
      ierr  = 0_ls_mpik
      nel   = int(nbuf,kind=ls_mpik)
      dtype = MPI_DOUBLE_PRECISION
      if(ls_mpik==4)then
        k=SPLIT_MPI_MSG
        do i=1,nbuf,k
          n4=k
          if(((nbuf-i)<k).and.(mod(nbuf-i+1,k)/=0))n4=mod(nbuf,k)
          call lsmpi_recv_realkV_4(buffer(i:i+n4-1),n4,comm,sender)
        enddo
      else
        call MPI_RECV(buffer,nbuf,dtype,sender,tag,comm,status,ierr)
        IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
      endif
#endif 
    end subroutine lsmpi_recv_realkV_8
    subroutine lsmpi_recv_realkV_4(buffer,nbuf,comm,sender)
      implicit none
      real(realk) :: buffer(:)
      integer(kind=4) :: nbuf
      integer(kind=ls_mpik) :: comm
      integer(kind=ls_mpik) :: sender
      integer(kind=ls_mpik) :: tag
      integer(kind=ls_mpik) :: nel,dtype,ierr
#ifdef VAR_MPI
      tag   = 124_ls_mpik
      ierr  = 0_ls_mpik
      nel   = int(nbuf,kind=ls_mpik)
      dtype = MPI_DOUBLE_PRECISION
      call MPI_RECV(buffer,nbuf,dtype,sender,tag,comm,status,ierr)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif 
    end subroutine lsmpi_recv_realkV_4
    

    subroutine ls_mpisendrecv_integer(buffer,comm,sender,receiver)
      implicit none
      integer(kind=4) :: buffer
      integer(kind=ls_mpik) :: comm   ! communicator
      integer(kind=ls_mpik) :: sender,receiver
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: ierr,thesize,datatype
      integer(kind=ls_mpik) :: mynum,tag

      IERR=0
      tag=0
      ! Get rank within specific communicator
      call get_rank_for_comm(comm,mynum)
      DATATYPE = MPI_INTEGER4
      THESIZE = 1

      if(mynum.EQ.sender) then ! send stuff to receiver
         call MPI_SEND(buffer,thesize,DATATYPE,receiver,tag,comm,ierr)
      else if(mynum.EQ.receiver) then  ! receive stuff from sender
         call MPI_RECV(buffer,thesize,DATATYPE,sender,tag,comm,status,ierr)
      else ! Error: Node should be either sender or receiver
         print '(a,3i6)', 'Rank,sender,receiver',mynum,sender,receiver
         call lsquit('ls_mpisendrecv_integer: &
              & Rank is neither sender nor receiver',-1)
      end if
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpisendrecv_integer

    subroutine ls_mpisendrecv_long(buffer,comm,sender,receiver)
      implicit none
      integer(kind=8) :: buffer
      integer(kind=ls_mpik) :: comm   ! communicator
      integer(kind=ls_mpik) :: sender,receiver
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: ierr,thesize,datatype
      integer(kind=ls_mpik) :: mynum,tag,dummystat
      IERR=0
      tag=1
      dummystat=0
      ! Get rank within specific communicator
      call get_rank_for_comm(comm,mynum)
      DATATYPE = MPI_INTEGER8
      THESIZE = 1

      if(mynum.EQ.sender) then ! send stuff to receiver
         call MPI_SEND(buffer,thesize,DATATYPE,receiver,tag,comm,ierr)
      else if(mynum.EQ.receiver) then  ! receive stuff from sender
         call MPI_RECV(buffer,thesize,DATATYPE,sender,tag,comm,status,ierr)
      else ! Error: Node should be either sender or receiver
         print '(a,3i6)', 'Rank,sender,receiver',mynum,sender,receiver
         call lsquit('ls_mpisendrecv_long: &
              & Rank is neither sender nor receiver',-1)
      end if
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpisendrecv_long

    subroutine ls_mpisendrecv_short(buffer,comm,sender,receiver)
      implicit none
      integer(kind=short) :: buffer
      integer(kind=ls_mpik) :: comm   ! communicator
      integer(kind=ls_mpik) :: sender,receiver
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: ierr,thesize,datatype
      integer :: integerbuffer
      integer(kind=ls_mpik) :: mynum,tag,dummystat
      IERR=0
      tag=2
      dummystat=0
      ! Get rank within specific communicator
      call get_rank_for_comm(comm,mynum)
      DATATYPE = MPI_INTEGER
      THESIZE = 1

      !     Convert from short integer to regular integer
      integerbuffer = buffer

      !     Send/receive regular integer
      if(mynum.EQ.sender) then ! send stuff to receiver
         call MPI_SEND(integerbuffer,thesize,DATATYPE,receiver,tag,comm,ierr)
      else if(mynum.EQ.receiver) then  ! receive stuff from sender
         call MPI_RECV(integerbuffer,thesize,DATATYPE,sender,tag,comm,status,ierr)
      else ! Error: Node should be either sender or receiver
         print '(a,3i6)', 'Rank,sender,receiver',mynum,sender,receiver
         call lsquit('ls_mpisendrecv_short: &
              & Rank is neither sender nor receiver',-1)
      end if
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)

      !     Convert back
      buffer = integerbuffer
#endif
    end subroutine ls_mpisendrecv_short

    subroutine ls_mpisendrecv_shortV_wrapper8(buffer,nbuf,comm,sender,receiver)
      implicit none
      integer(kind=short)   :: buffer(:)
      integer(kind=8)       :: nbuf
      integer(kind=ls_mpik) :: comm   ! communicator
      integer(kind=ls_mpik) :: sender,receiver
      integer(kind=ls_mpik) :: ierr,thesize,datatype
      integer(kind=ls_mpik) :: mynum,tag
      integer(kind=8) :: i,k,n
      integer(kind=4) :: n4
      integer(kind=long) :: n1
#ifdef VAR_MPI
      IERR=0
      if(ls_mpik==4)then
        n=nbuf
        k=SPLIT_MPI_MSG
        do i=1,n,k
          n4=k
          !if((n-i)<k)n4=mod(n,k)
          if(((n-i)<k).and.(mod(n-i+1,k)/=0))n4=mod(n,k)
          call ls_mpisendrecv_shortV(buffer(i:i+n4-1),n4,comm,sender,receiver)
        enddo
      else
        tag=3
        ! Get rank within specific communicator
        call get_rank_for_comm(comm,mynum)
       
        DATATYPE = MPI_INTEGER
        THESIZE = SIZE(buffer,kind=ls_mpik)
        n1 = SIZE(buffer,kind=long)
        if(THESIZE.NE.n1)call lsquit('lsmpi error1 in ls_mpisendrecv_shortV',-1)
        if(THESIZE.NE.nbuf)call lsquit('lsmpi error in ls_mpisendrecv_ShortV',-1)
       
        IF (mod(nbuf,int_to_short).NE.0) THEN
          write(*,*) 'Error in ls_mpisendrecv_shortV',nbuf,int_to_short,mod(nbuf,int_to_short)
          call lsquit('lsmpi nubf modular error in ls_mpisendrecv_ShortV',-1)
        ENDIF
        n = nbuf/int_to_short
       
        !     Send/receive regular integer
        if(mynum.EQ.sender) then ! send stuff to receiver
           call MPI_SEND(buffer(1:nbuf),n,DATATYPE,receiver,tag,comm,ierr)
        else if(mynum.EQ.receiver) then  ! receive stuff from sender
           call MPI_RECV(buffer(1:nbuf),n,DATATYPE,sender,tag,comm,status,ierr)
        else ! Error: Node should be either sender or receiver
           print '(a,3i6)', 'Rank,sender,receiver',mynum,sender,receiver
           call lsquit('ls_mpisendrecv_shortV: &
                & Rank is neither sender nor receiver',-1)
        end if
        IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
      endif
#endif
    end subroutine ls_mpisendrecv_shortV_wrapper8
    subroutine ls_mpisendrecv_shortV(buffer,nbuf,comm,sender,receiver)
      implicit none
      integer(kind=short)   :: buffer(:)
      integer(kind=4)       :: nbuf
      integer(kind=ls_mpik) :: comm   ! communicator
      integer(kind=ls_mpik) :: sender,receiver
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: ierr,thesize,datatype,n
      integer(kind=ls_mpik) :: mynum,tag
      integer(kind=long)    :: n1
      IERR=0
      tag=3
      ! Get rank within specific communicator
      call get_rank_for_comm(comm,mynum)

      DATATYPE = MPI_INTEGER
      THESIZE = SIZE(buffer,kind=ls_mpik)
      n1 = SIZE(buffer,kind=long)
      if(THESIZE.NE.n1)call lsquit('lsmpi error1 in ls_mpisendrecv_shortV',-1)
      if(THESIZE.NE.nbuf)call lsquit('lsmpi error in ls_mpisendrecv_ShortV',-1)

      IF (mod(nbuf,int_to_short).NE.0) THEN
        write(*,*) 'Error in ls_mpisendrecv_shortV',nbuf,int_to_short,mod(nbuf,int_to_short)
        call lsquit('lsmpi nubf modular error in ls_mpisendrecv_ShortV',-1)
      ENDIF
      n = nbuf/int_to_short

      !     Send/receive regular integer
      if(mynum.EQ.sender) then ! send stuff to receiver
         call MPI_SEND(buffer(1:nbuf),n,DATATYPE,receiver,tag,comm,ierr)
      else if(mynum.EQ.receiver) then  ! receive stuff from sender
         call MPI_RECV(buffer(1:nbuf),n,DATATYPE,sender,tag,comm,status,ierr)
      else ! Error: Node should be either sender or receiver
         print '(a,3i6)', 'Rank,sender,receiver',mynum,sender,receiver
         call lsquit('ls_mpisendrecv_shortV: &
              & Rank is neither sender nor receiver',-1)
      end if
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpisendrecv_shortV


    subroutine ls_mpisendrecv_longV_wrapper8(buffer,nbuf,comm,sender,receiver)
      implicit none
      integer(kind=8) :: nbuf
      integer(kind=8) :: buffer(:)
      integer(kind=ls_mpik) :: comm   ! communicator
      integer(kind=ls_mpik) :: sender,receiver
      integer(kind=ls_mpik) :: ierr,thesize,datatype
      integer(kind=ls_mpik) :: mynum,tag
      integer(kind=8) :: i,k
      integer(kind=4) :: n4
      integer(kind=long)    :: n1
#ifdef VAR_MPI
      IERR=0
      if(ls_mpik==4)then
        k=SPLIT_MPI_MSG
        do i=1,nbuf,k
          n4=k
          !if((nbuf-i)<k)n4=mod(nbuf,k)
          if(((nbuf-i)<k).and.(mod(nbuf-i+1,k)/=0))n4=mod(nbuf,k)
          call ls_mpisendrecv_longV(buffer(i:i+n4-1),n4,comm,sender,receiver)
        enddo
      else
        tag=4
        ! Get rank within specific communicator
        call get_rank_for_comm(comm,mynum)
       
        DATATYPE = MPI_INTEGER8
        THESIZE = SIZE(buffer,kind=ls_mpik)
        n1 = SIZE(buffer,kind=long)
        if(THESIZE.NE.n1)call lsquit('lsmpi error1 in ls_mpisendrecv_longV',-1)
        if(THESIZE.NE.nbuf)call lsquit('lsmpi error in ls_mpisendrecv_longV',-1)
       
        !     Send/receive 
        if(mynum.EQ.sender) then ! send stuff to receiver
           call MPI_SEND(buffer,thesize,DATATYPE,receiver,tag,comm,ierr)
        else if(mynum.EQ.receiver) then  ! receive stuff from sender
           call MPI_RECV(buffer,thesize,DATATYPE,sender,tag,comm,status,ierr)
        else ! Error: Node should be either sender or receiver
           print '(a,3i6)', 'Rank,sender,receiver',mynum,sender,receiver
           call lsquit('ls_mpisendrecv_longV: &
                & Rank is neither sender nor receiver',-1)
        end if
        IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
      endif
#endif
    end subroutine ls_mpisendrecv_longV_wrapper8
    subroutine ls_mpisendrecv_longV(buffer,nbuf,comm,sender,receiver)
      implicit none
      integer(kind=4) :: nbuf
      integer(kind=8) :: buffer(:)
      integer(kind=ls_mpik) :: comm   ! communicator
      integer(kind=ls_mpik) :: sender,receiver
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: ierr,thesize,datatype
      integer(kind=ls_mpik) :: mynum,tag
      integer(kind=long)    :: n1
      IERR=0
      tag=4
      ! Get rank within specific communicator
      call get_rank_for_comm(comm,mynum)

      DATATYPE = MPI_INTEGER8
      THESIZE = SIZE(buffer,kind=ls_mpik)
      n1 = SIZE(buffer,kind=long)
      if(THESIZE.NE.n1)call lsquit('lsmpi error1 in ls_mpisendrecv_longV',-1)
      if(THESIZE.NE.nbuf)call lsquit('lsmpi error in ls_mpisendrecv_longV',-1)

      !     Send/receive 
      if(mynum.EQ.sender) then ! send stuff to receiver
         call MPI_SEND(buffer,thesize,DATATYPE,receiver,tag,comm,ierr)
      else if(mynum.EQ.receiver) then  ! receive stuff from sender
         call MPI_RECV(buffer,thesize,DATATYPE,sender,tag,comm,status,ierr)
      else ! Error: Node should be either sender or receiver
         print '(a,3i6)', 'Rank,sender,receiver',mynum,sender,receiver
         call lsquit('ls_mpisendrecv_longV: &
              & Rank is neither sender nor receiver',-1)
      end if
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpisendrecv_longV


    subroutine ls_mpisendrecv_integerV_wrapper8(buffer,nbuf,comm,sender,receiver)
      implicit none
      integer(kind=8) :: nbuf
      integer(kind=4) :: buffer(:)
      integer(kind=ls_mpik) :: comm   ! communicator
      integer(kind=ls_mpik) :: sender,receiver
      integer(kind=ls_mpik) :: ierr,thesize,datatype
      integer(kind=ls_mpik) :: mynum,tag
      integer(kind=8) :: i,k
      integer(kind=4) :: n4
      integer(kind=long)    :: n1
#ifdef VAR_MPI
      IERR=0
      if(ls_mpik==4)then
        k=SPLIT_MPI_MSG
        do i=1,nbuf,k
          n4=k
          !if((nbuf-i)<k)n4=mod(nbuf,k)
          if(((nbuf-i)<k).and.(mod(nbuf-i+1,k)/=0))n4=mod(nbuf,k)
          call ls_mpisendrecv_integerV(buffer(i:i+n4-1),n4,comm,sender,receiver)
        enddo
      else
        tag=5
        ! Get rank within specific communicator
        call get_rank_for_comm(comm,mynum)
        DATATYPE = MPI_INTEGER4
        THESIZE = SIZE(buffer,kind=ls_mpik)
        n1 = SIZE(buffer,kind=long)
        if(THESIZE.NE.n1)call lsquit('lsmpi error1 in ls_mpisendrecv_integerV',-1)
        if(THESIZE.NE.nbuf)call lsquit('lsmpi error in ls_mpisendrecv_integerV',-1)
       
        !     Send/receive 
        if(mynum.EQ.sender) then ! send stuff to receiver
           call MPI_SEND(buffer,thesize,DATATYPE,receiver,tag,comm,ierr)
        else if(mynum.EQ.receiver) then  ! receive stuff from sender
           call MPI_RECV(buffer,thesize,DATATYPE,sender,tag,comm,status,ierr)
        else ! Error: Node should be either sender or receiver
         print '(a,3i6)', 'Rank,sender,receiver',mynum,sender,receiver
         call lsquit('ls_mpisendrecv_integerV: Rank not sender nor receiver',-1)
        end if
        IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
      endif
#endif
    end subroutine ls_mpisendrecv_integerV_wrapper8
    subroutine ls_mpisendrecv_integerV(buffer,nbuf,comm,sender,receiver)
      implicit none
      integer(kind=4) :: nbuf
      integer(kind=4) :: buffer(:)
      integer(kind=ls_mpik) :: comm   ! communicator
      integer(kind=ls_mpik) :: sender,receiver
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: ierr,thesize,datatype
      integer(kind=ls_mpik) :: mynum,tag
      integer(kind=long)    :: n1
      IERR=0
      tag=5
      ! Get rank within specific communicator
      call get_rank_for_comm(comm,mynum)
      DATATYPE = MPI_INTEGER4
      THESIZE = SIZE(buffer,kind=ls_mpik)
      n1 = SIZE(buffer,kind=long)
      if(THESIZE.NE.n1)call lsquit('lsmpi error1 in ls_mpisendrecv_integerV',-1)
      if(THESIZE.NE.nbuf)call lsquit('lsmpi error in ls_mpisendrecv_integerV',-1)

      !     Send/receive 
      if(mynum.EQ.sender) then ! send stuff to receiver
         call MPI_SEND(buffer,thesize,DATATYPE,receiver,tag,comm,ierr)
      else if(mynum.EQ.receiver) then  ! receive stuff from sender
         call MPI_RECV(buffer,thesize,DATATYPE,sender,tag,comm,status,ierr)
      else ! Error: Node should be either sender or receiver
       print '(a,3i6)', 'Rank,sender,receiver',mynum,sender,receiver
       call lsquit('ls_mpisendrecv_integerV: Rank not sender nor receiver',-1)
      end if
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpisendrecv_integerV

    subroutine ls_mpisendrecv_realk(buffer,comm,sender,receiver)
      implicit none
      real(realk) :: buffer
      integer(kind=ls_mpik) :: comm   ! communicator
      integer(kind=ls_mpik) :: sender,receiver
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: ierr,thesize,datatype
      integer(kind=ls_mpik) :: mynum,tag
      IERR=0
      tag=6
      ! Get rank within specific communicator
      call get_rank_for_comm(comm,mynum)

      DATATYPE = MPI_DOUBLE_PRECISION
      THESIZE = 1

      if(mynum.EQ.sender) then ! send stuff to receiver
         call MPI_SEND(buffer,thesize,DATATYPE,receiver,tag,comm,ierr)
      else if(mynum.EQ.receiver) then  ! receive stuff from sender
         call MPI_RECV(buffer,thesize,DATATYPE,sender,tag,comm,status,ierr)
      else ! Error: Node should be either sender or receiver
         print '(a,3i6)', 'Rank,sender,receiver',mynum,sender,receiver
         call lsquit('ls_mpisendrecv_realk: &
              & Rank is neither sender nor receiver',-1)
      end if
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpisendrecv_realk

    subroutine ls_mpisendrecv_realkV_wrapper8(buffer,nbuf,comm,sender,receiver)
      implicit none
      integer(kind=8) :: nbuf
      real(realk) :: buffer(:)
      integer(kind=ls_mpik) :: comm   ! communicator
      integer(kind=ls_mpik) :: sender,receiver
      integer(kind=ls_mpik) :: ierr,thesize,datatype
      integer(kind=ls_mpik) :: mynum,tag
      integer(kind=8) :: i,k
      integer(kind=4) :: n4
      integer(kind=long)    :: n1
#ifdef VAR_MPI
      IERR=0
      if(ls_mpik==4)then
        k=SPLIT_MPI_MSG
        do i=1,nbuf,k
          n4=k
          if(((nbuf-i)<k).and.(mod(nbuf-i+1,k)/=0))n4=mod(nbuf,k)
          call ls_mpisendrecv_realkV(buffer(i:i+n4-1),n4,comm,sender,receiver)
        enddo
      else
        tag=7
        ! Get rank within specific communicator
        call get_rank_for_comm(comm,mynum)
        DATATYPE = MPI_DOUBLE_PRECISION
        THESIZE = SIZE(buffer,kind=ls_mpik)
        n1 = SIZE(buffer,kind=long)
        if(THESIZE.NE.n1)call lsquit('lsmpi error1 in ls_mpisendrecv_realkV_wrapper8',-1)
        if(THESIZE.NE.nbuf)THEN
           call lsquit('lsmpi error in ls_mpisendrecv_realkV_wrapper8',-1)
        endif
       
        !     Send/receive 
        if(mynum.EQ.sender) then ! send stuff to receiver
           call MPI_SEND(buffer,thesize,DATATYPE,receiver,tag,comm,ierr)
        else if(mynum.EQ.receiver) then  ! receive stuff from sender
           call MPI_RECV(buffer,thesize,DATATYPE,sender,tag,comm,status,ierr)
        else ! Error: Node should be either sender or receiver
           print '(a,3i6)', 'Rank,sender,receiver',mynum,sender,receiver
           call lsquit('ls_mpisendrecv_realkV: &
                & Rank is neither sender nor receiver',-1)
        end if
        IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
      endif
#endif
    end subroutine ls_mpisendrecv_realkV_wrapper8
    subroutine ls_mpisendrecv_realkV(buffer,nbuf,comm,sender,receiver)
      implicit none
      integer(kind=4) :: nbuf
      real(realk) :: buffer(:)
      integer(kind=ls_mpik) :: comm   ! communicator
      integer(kind=ls_mpik) :: sender,receiver
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: ierr,thesize,datatype
      integer(kind=ls_mpik) :: mynum,tag
      integer(kind=long)    :: n1
      IERR=0
      tag=7
      ! Get rank within specific communicator
      call get_rank_for_comm(comm,mynum)
      DATATYPE = MPI_DOUBLE_PRECISION
      THESIZE = SIZE(buffer,kind=ls_mpik)
      n1 = SIZE(buffer,kind=long)
      if(THESIZE.NE.n1)call lsquit('lsmpi error1 in ls_mpisendrecv_realkV',-1)
      if(THESIZE.NE.nbuf)THEN
         call lsquit('lsmpi error in ls_mpisendrecv_realkV',-1)
      endif

      !     Send/receive 
      if(mynum.EQ.sender) then ! send stuff to receiver
         call MPI_SEND(buffer,thesize,DATATYPE,receiver,tag,comm,ierr)
      else if(mynum.EQ.receiver) then  ! receive stuff from sender
         call MPI_RECV(buffer,thesize,DATATYPE,sender,tag,comm,status,ierr)
      else ! Error: Node should be either sender or receiver
         print '(a,3i6)', 'Rank,sender,receiver',mynum,sender,receiver
         call lsquit('ls_mpisendrecv_realkV: &
              & Rank is neither sender nor receiver',-1)
      end if
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpisendrecv_realkV

    subroutine ls_mpisendrecv_realkM(buffer,nbuf1,nbuf2,comm,sender,receiver)
      implicit none
      integer :: nbuf1,nbuf2
      real(realk) :: buffer(:,:)
      integer(kind=ls_mpik) :: comm   ! communicator
      integer(kind=ls_mpik) :: sender,receiver
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: ierr,thesize,datatype
      integer(kind=ls_mpik) :: mynum,tag
      integer ::  I,J,offset
      integer(kind=long)    :: n1
      IERR=0
      tag=8
      ! Get rank within specific communicator
      call get_rank_for_comm(comm,mynum)

      DATATYPE = MPI_DOUBLE_PRECISION
      THESIZE = SIZE(buffer(:,:),kind=ls_mpik)
      n1 = SIZE(buffer,kind=long)
      if(THESIZE.NE.n1)call lsquit('lsmpi error1 in ls_mpisendrecv_realkM',-1)
      if(THESIZE.NE.nbuf1*nbuf2) call lsquit('lsmpi error in ls_mpisendrecv_realkM',-1)

      !     Send/receive 
      if(mynum.EQ.sender) then ! send stuff to receiver
         call MPI_SEND(buffer,thesize,DATATYPE,receiver,tag,comm,ierr)
      else if(mynum.EQ.receiver) then  ! receive stuff from sender
         call MPI_RECV(buffer,thesize,DATATYPE,sender,tag,comm,status,ierr)
      else ! Error: Node should be either sender or receiver
         print '(a,3i6)', 'Rank,sender,receiver',mynum,sender,receiver
         call lsquit('ls_mpisendrecv_realkM: &
              & Rank is neither sender nor receiver',-1)
      end if
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpisendrecv_realkM

    subroutine ls_mpisendrecv_realkT(buffer,nbuf1,nbuf2,nbuf3,comm,sender,receiver)
      implicit none
      integer :: nbuf1,nbuf2,nbuf3
      real(realk) :: buffer(:,:,:)
      integer(kind=ls_mpik) :: comm   ! communicator
      integer(kind=ls_mpik) :: sender,receiver
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: ierr,thesize,datatype
      integer(kind=ls_mpik) :: mynum,tag
      integer :: I,J,K,offset
      integer(kind=long)    :: n1
      IERR=0
      tag=9
      ! Get rank within specific communicator
      call get_rank_for_comm(comm,mynum)

      DATATYPE = MPI_DOUBLE_PRECISION
      THESIZE = SIZE(buffer,kind=ls_mpik)
      n1 = SIZE(buffer,kind=long)
      if(THESIZE.NE.n1)call lsquit('lsmpi error1 in ls_mpisendrecv_realT',-1)
      if(THESIZE.NE.nbuf1*nbuf2*nbuf3)THEN
         call lsquit('lsmpi error in ls_mpisendrecv_realkT',-1)
      endif

      !     Send/receive 
      if(mynum.EQ.sender) then ! send stuff to receiver
         call MPI_SEND(buffer,thesize,DATATYPE,receiver,tag,comm,ierr)
      else if(mynum.EQ.receiver) then  ! receive stuff from sender
         call MPI_RECV(buffer,thesize,DATATYPE,sender,tag,comm,status,ierr)
      else ! Error: Node should be either sender or receiver
         print '(a,3i6)', 'Rank,sender,receiver',mynum,sender,receiver
         call lsquit('ls_mpisendrecv_realkT: &
              & Rank is neither sender nor receiver',-1)
      end if
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)

#endif
    end subroutine ls_mpisendrecv_realkT

    subroutine ls_mpisendrecv_realkQ(buffer,nbuf1,nbuf2,nbuf3,nbuf4,comm,sender,receiver)
      implicit none
      integer :: nbuf1,nbuf2,nbuf3,nbuf4
      real(realk) :: buffer(:,:,:,:)
      integer(kind=ls_mpik) :: comm   ! communicator
      integer(kind=ls_mpik) :: sender,receiver
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: ierr,thesize,datatype
      integer(kind=ls_mpik) :: mynum,tag
      integer :: I,J,K,offset
      integer(kind=long)    :: n1
      IERR=0
      tag=10
      ! Get rank within specific communicator
      call get_rank_for_comm(comm,mynum)

      DATATYPE = MPI_DOUBLE_PRECISION
      THESIZE = SIZE(buffer,kind=ls_mpik)
      n1 = SIZE(buffer,kind=long)
      if(THESIZE.NE.n1)call lsquit('lsmpi error1 in ls_mpisendrecv_realkQ',-1)
      if(THESIZE.NE.nbuf1*nbuf2*nbuf3*nbuf4)THEN
         call lsquit('lsmpi error in ls_mpisendrecv_realkQ',-1)
      endif

      if(mynum.EQ.sender) then ! send stuff to receiver
         call MPI_SEND(buffer,thesize,DATATYPE,receiver,tag,comm,ierr)
      else if(mynum.EQ.receiver) then  ! receive stuff from sender
         call MPI_RECV(buffer,thesize,DATATYPE,sender,tag,comm,status,ierr)
      else ! Error: Node should be either sender or receiver
         print '(a,3i6)', 'Rank,sender,receiver',mynum,sender,receiver
         call lsquit('ls_mpisendrecv_realkQ: &
              & Rank is neither sender nor receiver',-1)
      end if
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpisendrecv_realkQ

    subroutine ls_mpisendrecv_logical8(buffer,comm,sender,receiver)
      implicit none
      logical(kind=8) :: buffer
      integer(kind=ls_mpik) :: comm   ! communicator
      integer(kind=ls_mpik) :: sender,receiver
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: ierr,thesize,datatype
      integer(kind=ls_mpik) :: mynum,tag
      integer(kind=MPI_ADDRESS_KIND) :: mpi_logical_extent,lb
      logical(kind=4) :: buffer4
      IERR=0
      tag=11
      ! Get rank within specific communicator
      call get_rank_for_comm(comm,mynum)

      DATATYPE = MPI_LOGICAL
      THESIZE = 1
      !     Send/receive 
      call MPI_TYPE_GET_EXTENT(MPI_LOGICAL,lb,mpi_logical_extent,ierr)
      if(mynum.EQ.sender) then ! send stuff to receiver
         IF(mpi_logical_extent.EQ.4)THEN
            !32 bit mpi logical
            BUFFER4 = BUFFER
            call MPI_SEND(buffer4,thesize,DATATYPE,receiver,tag,comm,ierr)
         ELSE
            !64 bit mpi logical
            call MPI_SEND(buffer,thesize,DATATYPE,receiver,tag,comm,ierr)
         ENDIF
      else if(mynum.EQ.receiver) then  ! receive stuff from sender
         IF(mpi_logical_extent.EQ.4)THEN
            !32 bit mpi logical
            call MPI_RECV(buffer4,thesize,DATATYPE,sender,tag,comm,status,ierr)
            BUFFER = BUFFER4
         ELSE
            !64 bit mpi logical
            call MPI_RECV(buffer,thesize,DATATYPE,sender,tag,comm,status,ierr)
         ENDIF
      else ! Error: Node should be either sender or receiver
         print '(a,3i6)', 'Rank,sender,receiver',mynum,sender,receiver
         call lsquit('ls_mpisendrecv_logical: &
              & Rank is neither sender nor receiver',-1)
      end if
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpisendrecv_logical8
    subroutine ls_mpisendrecv_logical4(buffer,comm,sender,receiver)
      implicit none
      logical(kind=4) :: buffer
      integer(kind=ls_mpik) :: comm   ! communicator
      integer(kind=ls_mpik) :: sender,receiver
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: ierr,thesize,datatype
      integer(kind=ls_mpik) :: mynum,tag
      logical(kind=8) :: buffer8
      integer(kind=MPI_ADDRESS_KIND) :: mpi_logical_extent,lb
      IERR=0
      tag=11
      ! Get rank within specific communicator
      call get_rank_for_comm(comm,mynum)

      DATATYPE = MPI_LOGICAL
      THESIZE = 1

      !     Send/receive 
      call MPI_TYPE_GET_EXTENT(MPI_LOGICAL,lb,mpi_logical_extent,ierr)
      if(mynum.EQ.sender) then ! send stuff to receiver
         IF(mpi_logical_extent.EQ.4)THEN
            !32 bit mpi logical
            call MPI_SEND(buffer,thesize,DATATYPE,receiver,tag,comm,ierr)
         ELSE
            BUFFER8 = BUFFER
            call MPI_SEND(buffer8,thesize,DATATYPE,receiver,tag,comm,ierr)
         ENDIF
      else if(mynum.EQ.receiver) then  ! receive stuff from sender
         IF(mpi_logical_extent.EQ.4)THEN
            !32 bit mpi logical
            call MPI_RECV(buffer,thesize,DATATYPE,sender,tag,comm,status,ierr)
         ELSE
            call MPI_RECV(buffer8,thesize,DATATYPE,sender,tag,comm,status,ierr)
            BUFFER = BUFFER8
         ENDIF
      else ! Error: Node should be either sender or receiver
         print '(a,3i6)', 'Rank,sender,receiver',mynum,sender,receiver
         call lsquit('ls_mpisendrecv_logical: &
              & Rank is neither sender nor receiver',-1)
      end if
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpisendrecv_logical4
    subroutine ls_mpisendrecv_logical8V_wrapper8(buffer,nbuf,comm,sender,receiver)
      implicit none
      integer(kind=8) :: nbuf
      logical(kind=8) :: buffer(:)
      integer(kind=ls_mpik) :: comm   ! communicator
      integer(kind=ls_mpik) :: sender,receiver
      integer(kind=ls_mpik) :: ierr,thesize,datatype
      integer(kind=ls_mpik) :: mynum,tag
      integer(kind=8) :: i,k
      integer(kind=4) :: n4
#ifdef VAR_MPI
      logical(kind=4),pointer :: buffer4(:)
      integer(kind=MPI_ADDRESS_KIND) :: mpi_logical_extent,lb
      integer(kind=long)    :: n1
      IERR=0
      if(ls_mpik==4)then
        k=SPLIT_MPI_MSG
        do i=1,nbuf,k
          n4=k
          if(((nbuf-i)<k).and.(mod(nbuf-i+1,k)/=0))n4=mod(nbuf,k)
          call ls_mpisendrecv_logical8V(buffer(i:i+n4-1),n4,comm,sender,receiver)
        enddo
      else
        tag=12
        ! Get rank within specific communicator
        call get_rank_for_comm(comm,mynum)
       
        DATATYPE = MPI_LOGICAL
        THESIZE = SIZE(buffer,kind=ls_mpik)
        n1 = SIZE(buffer,kind=long)
        if(THESIZE.NE.n1)call lsquit('lsmpi error1 in ls_mpisendrecv_logical8V_wrapper',-1)
        if(THESIZE.NE.nbuf)THEN
           call lsquit('lsmpi error in ls_mpisendrecv_logical8V_wrapper8',-1)
        endif
        THESIZE = nbuf
       
        !     Send/receive 
        call MPI_TYPE_GET_EXTENT(MPI_LOGICAL,lb,mpi_logical_extent,ierr)
        if(mynum.EQ.sender) then ! send stuff to receiver
           IF(mpi_logical_extent.EQ.4)THEN
              !32 bit mpi logical
              call mem_alloc(buffer4,nbuf)
              do I = 1,thesize
                 buffer4(I) = buffer(I)
              enddo
              call MPI_SEND(buffer4,thesize,DATATYPE,receiver,tag,comm,ierr)
              call mem_dealloc(buffer4)
           ELSE
              !64 bit mpi logical
              call MPI_SEND(buffer,thesize,DATATYPE,receiver,tag,comm,ierr)
           ENDIF
        else if(mynum.EQ.receiver) then  ! receive stuff from sender
           IF(mpi_logical_extent.EQ.4)THEN
              !32 bit mpi logical
              call mem_alloc(buffer4,nbuf)
              call MPI_RECV(buffer4,thesize,DATATYPE,sender,tag,comm,status,ierr)
              do I = 1,thesize
                 buffer(I) = buffer4(I)
              enddo
              call mem_dealloc(buffer4)
           ELSE
              !64 bit mpi logical
              call MPI_RECV(buffer,thesize,DATATYPE,sender,tag,comm,status,ierr)
           ENDIF
        else ! Error: Node should be either sender or receiver
           print '(a,3i6)', 'Rank,sender,receiver',mynum,sender,receiver
           call lsquit('ls_mpisendrecv_logical8V_wrapper8: &
                & Rank is neither sender nor receiver',-1)
        end if
        IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
      endif
#endif
    end subroutine ls_mpisendrecv_logical8V_wrapper8

    subroutine ls_mpisendrecv_logical8V(buffer,nbuf,comm,sender,receiver)
      implicit none
      integer(kind=4) :: nbuf
      logical(kind=8) :: buffer(:)
      integer(kind=ls_mpik) :: comm   ! communicator
      integer(kind=ls_mpik) :: sender,receiver
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: ierr,thesize,datatype
      integer(kind=ls_mpik) :: mynum,tag
      logical(kind=4),pointer :: buffer4(:)
      integer :: I
      integer(kind=MPI_ADDRESS_KIND) :: mpi_logical_extent,lb
      integer(kind=long)    :: n1
      IERR=0
      tag=12
      ! Get rank within specific communicator
      call get_rank_for_comm(comm,mynum)

      DATATYPE = MPI_LOGICAL
      THESIZE = SIZE(buffer,kind=ls_mpik)
      n1 = SIZE(buffer,kind=long)
      if(THESIZE.NE.n1)call lsquit('lsmpi error1 in ls_mpisendrecv_logical8V',-1)
      if(THESIZE.NE.nbuf)THEN
         call lsquit('lsmpi error in ls_mpisendrecv_logical8V',-1)
      endif
      THESIZE = nbuf

      !     Send/receive 
      call MPI_TYPE_GET_EXTENT(MPI_LOGICAL,lb,mpi_logical_extent,ierr)
      if(mynum.EQ.sender) then ! send stuff to receiver
           IF(mpi_logical_extent.EQ.4)THEN
              !32 bit mpi logical
              call mem_alloc(buffer4,nbuf)
              do I = 1,nbuf
                 buffer4(I) = buffer(I)
              enddo
              call MPI_SEND(buffer4,thesize,DATATYPE,receiver,tag,comm,ierr)
              call mem_dealloc(buffer4)
           ELSE
              !64 bit mpi logical
              call MPI_SEND(buffer,thesize,DATATYPE,receiver,tag,comm,ierr)
           ENDIF
      else if(mynum.EQ.receiver) then  ! receive stuff from sender
           IF(mpi_logical_extent.EQ.4)THEN
              !32 bit mpi logical
              call mem_alloc(buffer4,nbuf)
              call MPI_RECV(buffer4,thesize,DATATYPE,sender,tag,comm,status,ierr)
              do I = 1,nbuf
                 buffer(I) = buffer4(I)
              enddo
              call mem_dealloc(buffer4)
           ELSE
              !64 bit mpi logical
              call MPI_RECV(buffer,thesize,DATATYPE,sender,tag,comm,status,ierr)
           ENDIF
      else ! Error: Node should be either sender or receiver
         print '(a,3i6)', 'Rank,sender,receiver',mynum,sender,receiver
         call lsquit('ls_mpisendrecv_logical8V: &
              & Rank is neither sender nor receiver',-1)
      end if
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpisendrecv_logical8V

    subroutine ls_mpisendrecv_logical4V_wrapper8(buffer,nbuf,comm,sender,receiver)
      implicit none
      integer(kind=8) :: nbuf
      logical(kind=4) :: buffer(:)
      integer(kind=ls_mpik) :: comm   ! communicator
      integer(kind=ls_mpik) :: sender,receiver
      integer(kind=ls_mpik) :: ierr,thesize,datatype
      integer(kind=ls_mpik) :: mynum,tag
      integer(kind=8) :: i,k
      integer(kind=4) :: n4
#ifdef VAR_MPI
      logical(kind=8),pointer :: buffer8(:)
      integer(kind=MPI_ADDRESS_KIND) :: mpi_logical_extent,lb 
      integer(kind=long)    :: n1
      IERR=0
      if(ls_mpik==4)then
        k=SPLIT_MPI_MSG
        do i=1,nbuf,k
          n4=k
          if(((nbuf-i)<k).and.(mod(nbuf-i+1,k)/=0))n4=mod(nbuf,k)
          call ls_mpisendrecv_logical4V(buffer(i:i+n4-1),n4,comm,sender,receiver)
        enddo
      else
        tag=12
        ! Get rank within specific communicator
        call get_rank_for_comm(comm,mynum)
       
        DATATYPE = MPI_LOGICAL
        THESIZE = SIZE(buffer,kind=ls_mpik)
        n1 = SIZE(buffer,kind=long)
        if(THESIZE.NE.n1)call lsquit('lsmpi error1 in ls_mpisendrecv_logicalV',-1)
        if(THESIZE.NE.nbuf)THEN
           call lsquit('lsmpi error in ls_mpisendrecv_logicalV',-1)
        endif
        THESIZE = nbuf
       
        call MPI_TYPE_GET_EXTENT(MPI_LOGICAL,lb,mpi_logical_extent,ierr)
        !     Send/receive 
        if(mynum.EQ.sender) then ! send stuff to receiver
           IF(mpi_logical_extent.EQ.4)THEN
              !32 bit mpi logical
              call MPI_SEND(buffer,thesize,DATATYPE,receiver,tag,comm,ierr)
           ELSE
              !64 bit mpi logical
              call mem_alloc(buffer8,nbuf)
              do I = 1,nbuf
                 buffer8(I) = buffer(I)
              enddo
              call MPI_SEND(buffer8,thesize,DATATYPE,receiver,tag,comm,ierr)
              call mem_dealloc(buffer8)
           ENDIF
        else if(mynum.EQ.receiver) then  ! receive stuff from sender
           IF(mpi_logical_extent.EQ.4)THEN
              !32 bit mpi logical
              call MPI_RECV(buffer,thesize,DATATYPE,sender,tag,comm,status,ierr)
           ELSE
              !64 bit mpi logical
              call mem_alloc(buffer8,nbuf)
              call MPI_RECV(buffer8,thesize,DATATYPE,sender,tag,comm,status,ierr)
              do I = 1,nbuf
                 buffer(I) = buffer8(I)
              enddo
              call mem_dealloc(buffer8)
           ENDIF
        else ! Error: Node should be either sender or receiver
           print '(a,3i6)', 'Rank,sender,receiver',mynum,sender,receiver
           call lsquit('ls_mpisendrecv_logicalV: &
                & Rank is neither sender nor receiver',-1)
        end if
        IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
      endif
#endif
    end subroutine ls_mpisendrecv_logical4V_wrapper8

    subroutine ls_mpisendrecv_logical4V(buffer,nbuf,comm,sender,receiver)
      implicit none
      integer(kind=4) :: nbuf
      logical(kind=4) :: buffer(:)
      integer(kind=ls_mpik) :: comm   ! communicator
      integer(kind=ls_mpik) :: sender,receiver
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: ierr,thesize,datatype
      integer(kind=ls_mpik) :: mynum,tag
      logical(kind=8),pointer :: buffer8(:)
      integer(kind=8) :: nbuf8
      integer :: I
      integer(kind=MPI_ADDRESS_KIND) :: mpi_logical_extent,lb
      integer(kind=long)    :: n1
      IERR=0
      tag=12
      ! Get rank within specific communicator
      call get_rank_for_comm(comm,mynum)

      DATATYPE = MPI_LOGICAL
      THESIZE = SIZE(buffer,kind=ls_mpik)
      n1 = SIZE(buffer,kind=long)
      if(THESIZE.NE.n1)call lsquit('lsmpi error1 in ls_mpisendrecv_logicalV',-1)
      if(THESIZE.NE.nbuf)THEN
         call lsquit('lsmpi error in ls_mpisendrecv_logicalV',-1)
      endif
      THESIZE = nbuf

      call MPI_TYPE_GET_EXTENT(MPI_LOGICAL,lb,mpi_logical_extent,ierr)

      !     Send/receive 
      if(mynum.EQ.sender) then ! send stuff to receiver
         IF(mpi_logical_extent.EQ.4)THEN
            !32 bit mpi logical
            call MPI_SEND(buffer,thesize,DATATYPE,receiver,tag,comm,ierr)
         ELSE
            !64 bit mpi logical
            nbuf8 = nbuf
            call mem_alloc(buffer8,nbuf8)
            do I = 1,nbuf8
               buffer8(I) = buffer(I)
            enddo
            call MPI_SEND(buffer8,thesize,DATATYPE,receiver,tag,comm,ierr)
            call mem_dealloc(buffer8)
         ENDIF
      else if(mynum.EQ.receiver) then  ! receive stuff from sender
         IF(mpi_logical_extent.EQ.4)THEN
            !32 bit mpi logical
            call MPI_RECV(buffer,thesize,DATATYPE,sender,tag,comm,status,ierr)
         ELSE
            !64 bit mpi logical
            nbuf8 = nbuf
            call mem_alloc(buffer8,nbuf8)
            call MPI_RECV(buffer8,thesize,DATATYPE,sender,tag,comm,status,ierr)
            do I = 1,nbuf8
               buffer(I) = buffer8(I)
            enddo
            call mem_dealloc(buffer8)
         ENDIF
      else ! Error: Node should be either sender or receiver
         print '(a,3i6)', 'Rank,sender,receiver',mynum,sender,receiver
         call lsquit('ls_mpisendrecv_logicalV: &
              & Rank is neither sender nor receiver',-1)
      end if
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpisendrecv_logical4V

    subroutine ls_mpisendrecv_logical4M(buffer,nbuf1,nbuf2,comm,sender,receiver)
      implicit none
      integer :: nbuf1,nbuf2
      logical(kind=4) :: buffer(:,:)
      integer(kind=ls_mpik) :: comm   ! communicator
      integer(kind=ls_mpik) :: sender,receiver
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: ierr,thesize,datatype
      integer(kind=ls_mpik) :: mynum,tag
      logical(kind=8),pointer :: buffer8(:,:)
      integer :: I,J
      integer(kind=MPI_ADDRESS_KIND) :: mpi_logical_extent,lb
      integer(kind=long)    :: n1
      IERR=0
      tag=13
      ! Get rank within specific communicator
      call get_rank_for_comm(comm,mynum)

      DATATYPE = MPI_LOGICAL
      THESIZE = SIZE(buffer,kind=ls_mpik)
      n1 = SIZE(buffer,kind=long)
      if(THESIZE.NE.n1)call lsquit('lsmpi error1 in ls_mpisendrecv_logical4M',-1)
      if(THESIZE.NE.nbuf1*nbuf2)THEN
         call lsquit('lsmpi error in ls_mpisendrecv_logical4M',-1)
      endif
      THESIZE = nbuf1*nbuf2

      call MPI_TYPE_GET_EXTENT(MPI_LOGICAL,lb,mpi_logical_extent,ierr)
      !     Send/receive 
      if(mynum.EQ.sender) then ! send stuff to receiver
         IF(mpi_logical_extent.EQ.4)THEN
            !32 bit mpi logical
            call MPI_SEND(buffer,thesize,DATATYPE,receiver,tag,comm,ierr)
         ELSE
            !64 bit mpi logical
            call mem_alloc(buffer8,nbuf1,nbuf2)
            do J = 1,nbuf2
               do I = 1,nbuf1
                  buffer8(I,J) = buffer(I,J)
               enddo
            enddo
            call MPI_SEND(buffer8,thesize,DATATYPE,receiver,tag,comm,ierr)
            call mem_dealloc(buffer8)
         ENDIF
      else if(mynum.EQ.receiver) then  ! receive stuff from sender
         IF(mpi_logical_extent.EQ.4)THEN
            !32 bit mpi logical
            call MPI_RECV(buffer,thesize,DATATYPE,sender,tag,comm,status,ierr)
         ELSE
            !64 bit mpi logical
            call mem_alloc(buffer8,nbuf1,nbuf2)
            call MPI_RECV(buffer8,thesize,DATATYPE,sender,tag,comm,status,ierr)
            do J = 1,nbuf2
               do I = 1,nbuf1
                  buffer(I,J) = buffer8(I,J)
               enddo
            enddo
            call mem_dealloc(buffer8)
         ENDIF
      else ! Error: Node should be either sender or receiver
         print '(a,3i6)', 'Rank,sender,receiver',mynum,sender,receiver
         call lsquit('ls_mpisendrecv_logical4M: &
              & Rank is neither sender nor receiver',-1)
      end if
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpisendrecv_logical4M

    subroutine ls_mpisendrecv_logical8M(buffer,nbuf1,nbuf2,comm,sender,receiver)
      implicit none
      integer :: nbuf1,nbuf2
      logical(kind=8) :: buffer(:,:)
      integer(kind=ls_mpik) :: comm   ! communicator
      integer(kind=ls_mpik) :: sender,receiver
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: ierr,thesize,datatype
      integer(kind=ls_mpik) :: mynum,tag
      logical(kind=4),pointer :: buffer4(:,:)
      integer :: I,J 
      integer(kind=MPI_ADDRESS_KIND) :: mpi_logical_extent,lb
      integer(kind=long)    :: n1
      IERR=0
      tag=13
      ! Get rank within specific communicator
      call get_rank_for_comm(comm,mynum)

      DATATYPE = MPI_LOGICAL
      THESIZE = SIZE(buffer,kind=ls_mpik)
      n1 = SIZE(buffer,kind=long)
      if(THESIZE.NE.n1)call lsquit('lsmpi error1 in ls_mpisendrecv_logical8M',-1)
      if(THESIZE.NE.nbuf1*nbuf2)THEN
         call lsquit('lsmpi error in ls_mpisendrecv_logical8M',-1)
      endif
      THESIZE = nbuf1*nbuf2

      call MPI_TYPE_GET_EXTENT(MPI_LOGICAL,lb,mpi_logical_extent,ierr)
      !     Send/receive 
      if(mynum.EQ.sender) then ! send stuff to receiver
         IF(mpi_logical_extent.EQ.4)THEN
            !32 bit mpi logical
            call mem_alloc(buffer4,nbuf1,nbuf2)
            do J = 1,nbuf2
               do I = 1,nbuf1
                  buffer4(I,J) = buffer(I,J)
               enddo
            enddo
            call MPI_SEND(buffer4,thesize,DATATYPE,receiver,tag,comm,ierr)
            call mem_dealloc(buffer4)
         ELSE
            !64 bit mpi logical
            call MPI_SEND(buffer,thesize,DATATYPE,receiver,tag,comm,ierr)
         ENDIF
      else if(mynum.EQ.receiver) then  ! receive stuff from sender
         IF(mpi_logical_extent.EQ.4)THEN
            !32 bit mpi logical
            call mem_alloc(buffer4,nbuf1,nbuf2)
            call MPI_RECV(buffer4,thesize,DATATYPE,sender,tag,comm,status,ierr)
            do J = 1,nbuf2
               do I = 1,nbuf1
                  buffer(I,J) = buffer4(I,J)
               enddo
            enddo
            call mem_dealloc(buffer4)
         ELSE
            !64 bit mpi logical
            call MPI_RECV(buffer,thesize,DATATYPE,sender,tag,comm,status,ierr)
         ENDIF
      else ! Error: Node should be either sender or receiver
         print '(a,3i6)', 'Rank,sender,receiver',mynum,sender,receiver
         call lsquit('ls_mpisendrecv_logical8M: &
              & Rank is neither sender nor receiver',-1)
      end if
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpisendrecv_logical8M

    subroutine ls_mpisendrecv_charac(buffer,comm,sender,receiver)
      implicit none
      character :: buffer
      integer(kind=ls_mpik) :: comm   ! communicator
      integer(kind=ls_mpik) :: sender,receiver
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: ierr,thesize,datatype
      integer(kind=ls_mpik) :: mynum,tag
      IERR=0
      tag=14
      ! Get rank within specific communicator
      call get_rank_for_comm(comm,mynum)

      DATATYPE = MPI_CHARACTER
      THESIZE = 1

      !     Send/receive 
      if(mynum.EQ.sender) then ! send stuff to receiver
         call MPI_SEND(buffer,thesize,DATATYPE,receiver,tag,comm,ierr)
      else if(mynum.EQ.receiver) then  ! receive stuff from sender
         call MPI_RECV(buffer,thesize,DATATYPE,sender,tag,comm,status,ierr)
      else ! Error: Node should be either sender or receiver
         print '(a,3i6)', 'Rank,sender,receiver',mynum,sender,receiver
         call lsquit('ls_mpisendrecv_charac: &
              & Rank is neither sender nor receiver',-1)
      end if
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpisendrecv_charac

    subroutine ls_mpisendrecv_characV_wrapper8(buffer,nbuf,comm,sender,receiver)
      implicit none
      integer(kind=8) :: nbuf
      character*(*) :: buffer
      integer(kind=ls_mpik) :: comm   ! communicator
      integer(kind=ls_mpik) :: sender,receiver
      character :: tmpbuffer
      integer(kind=ls_mpik) :: ierr,thesize,datatype
      integer(kind=ls_mpik) :: mynum,tag
      integer(kind=8) :: i,k
      integer(kind=4) :: n4
#ifdef VAR_MPI
      IERR=0
      if(ls_mpik==4)then
        k=SPLIT_MPI_MSG
        do i=1,nbuf,k
          n4=k
          if(((nbuf-i)<k).and.(mod(nbuf-i+1,k)/=0))n4=mod(nbuf,k)
          call ls_mpisendrecv_characV(buffer(i:i+n4-1),n4,comm,sender,receiver)
        enddo
      else
        tag=15
        ! Get rank within specific communicator
        call get_rank_for_comm(comm,mynum)
       
        DATATYPE = MPI_CHARACTER
        THESIZE = nbuf
       
        !     Send/receive 
        if(mynum.EQ.sender) then ! send stuff to receiver
           call MPI_SEND(buffer,thesize,DATATYPE,receiver,tag,comm,ierr)
        else if(mynum.EQ.receiver) then  ! receive stuff from sender
           call MPI_RECV(buffer,thesize,DATATYPE,sender,tag,comm,status,ierr)
        else ! Error: Node should be either sender or receiver
           print '(a,3i6)', 'Rank,sender,receiver',mynum,sender,receiver
           call lsquit('ls_mpisendrecv_characV: &
                & Rank is neither sender nor receiver',-1)
        end if
        IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
      endif
#endif
    end subroutine ls_mpisendrecv_characV_wrapper8
    subroutine ls_mpisendrecv_characV(buffer,nbuf,comm,sender,receiver)
      implicit none
      integer(kind=4) :: nbuf
      character*(*) :: buffer
      integer(kind=ls_mpik) :: comm   ! communicator
      integer(kind=ls_mpik) :: sender,receiver
      character :: tmpbuffer
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: ierr,thesize,datatype
      integer(kind=ls_mpik) :: mynum,tag
      IERR=0
      tag=15
      ! Get rank within specific communicator
      call get_rank_for_comm(comm,mynum)

      DATATYPE = MPI_CHARACTER
      THESIZE = nbuf

      !     Send/receive 
      if(mynum.EQ.sender) then ! send stuff to receiver
         call MPI_SEND(buffer,thesize,DATATYPE,receiver,tag,comm,ierr)
      else if(mynum.EQ.receiver) then  ! receive stuff from sender
         call MPI_RECV(buffer,thesize,DATATYPE,sender,tag,comm,status,ierr)
      else ! Error: Node should be either sender or receiver
         print '(a,3i6)', 'Rank,sender,receiver',mynum,sender,receiver
         call lsquit('ls_mpisendrecv_characV: &
              & Rank is neither sender nor receiver',-1)
      end if
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpisendrecv_characV

    subroutine ls_mpisendrecv_characV2_wrapper8(buffer,nbuf,comm,sender,receiver)
      implicit none
      integer(kind=8) :: nbuf
      character :: buffer(nbuf)
      integer(kind=ls_mpik) :: comm   ! communicator
      integer(kind=ls_mpik) :: sender,receiver
      character :: tmpbuffer
      integer(kind=ls_mpik) :: ierr,thesize,datatype
      integer(kind=ls_mpik) :: mynum,tag
      integer(kind=8) :: i,k
      integer(kind=4) :: n4
#ifdef VAR_MPI
      IERR=0
      if(ls_mpik==4)then
        k=SPLIT_MPI_MSG
        do i=1,nbuf,k
          n4=k
          if(((nbuf-i)<k).and.(mod(nbuf-i+1,k)/=0))n4=mod(nbuf,k)
          call ls_mpisendrecv_characV2(buffer(i:i+n4-1),n4,comm,sender,receiver)
        enddo
      else
        tag=16
        ! Get rank within specific communicator
        call get_rank_for_comm(comm,mynum)
       
        DATATYPE = MPI_CHARACTER
        THESIZE = nbuf
       
        !     Send/receive 
        if(mynum.EQ.sender) then ! send stuff to receiver
           call MPI_SEND(buffer,thesize,DATATYPE,receiver,tag,comm,ierr)
        else if(mynum.EQ.receiver) then  ! receive stuff from sender
           call MPI_RECV(buffer,thesize,DATATYPE,sender,tag,comm,status,ierr)
        else ! Error: Node should be either sender or receiver
           print '(a,3i6)', 'Rank,sender,receiver',mynum,sender,receiver
           call lsquit('ls_mpisendrecv_characV2: &
                & Rank is neither sender nor receiver',-1)
        end if
        IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
      endif 
#endif
    end subroutine ls_mpisendrecv_characV2_wrapper8
    subroutine ls_mpisendrecv_characV2(buffer,nbuf,comm,sender,receiver)
      implicit none
      integer(kind=4) :: nbuf
      character :: buffer(nbuf)
      integer(kind=ls_mpik) :: comm   ! communicator
      integer(kind=ls_mpik) :: sender,receiver
      character :: tmpbuffer
      integer(kind=ls_mpik) :: ierr,thesize,datatype
      integer(kind=ls_mpik) :: mynum,tag
#ifdef VAR_MPI
      IERR=0
      tag=16
      ! Get rank within specific communicator
      call get_rank_for_comm(comm,mynum)

      DATATYPE = MPI_CHARACTER
      THESIZE = nbuf

      !     Send/receive 
      if(mynum.EQ.sender) then ! send stuff to receiver
         call MPI_SEND(buffer,thesize,DATATYPE,receiver,tag,comm,ierr)
      else if(mynum.EQ.receiver) then  ! receive stuff from sender
         call MPI_RECV(buffer,thesize,DATATYPE,sender,tag,comm,status,ierr)
      else ! Error: Node should be either sender or receiver
         print '(a,3i6)', 'Rank,sender,receiver',mynum,sender,receiver
         call lsquit('ls_mpisendrecv_characV2: &
              & Rank is neither sender nor receiver',-1)
      end if
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
    end subroutine ls_mpisendrecv_characV2

!##########################################################
!#
!# SUBROUTINES THAT FILL OR EXTRACT FROM A BUFFER 
!#
!##########################################################
    subroutine ls_mpi_buffer4_integer(buffer,master)
      implicit none
      integer(kind=4) :: buffer
      integer(kind=ls_mpik) :: master
#ifdef VAR_MPI
      IF(AddToBuffer)THEN
         IF(iInt4+1 .GT. nInteger4)call increaselsmpibufferInt4(1_long)
         IF(iInt4+1.GT.size(lsmpibufferInt4,kind=long))call lsquit('errorTK',-1)
         lsmpibufferInt4(iInt4+1) = buffer
         iInt4 = iInt4 + 1
      ELSE
         IF(iInt4+1 .GT. nInteger4)call lsquit('ls_mpi_buffer_integer: error using buffer',-1)
         buffer = lsmpibufferInt4(iInt4+1)
         iInt4 = iInt4 + 1
      ENDIF
#endif
    end subroutine ls_mpi_buffer4_integer

    subroutine ls_mpi_buffer8_integer(buffer,master)
      implicit none
      integer(kind=8) :: buffer
      integer(kind=ls_mpik) :: master
#ifdef VAR_MPI
      IF(AddToBuffer)THEN
         IF(iInt8+1 .GT. nInteger8)call increaselsmpibufferInt8(1_long)
         lsmpibufferInt8(iInt8+1) = buffer
         iInt8 = iInt8 + 1
      ELSE
         IF(iInt8+1 .GT. nInteger8)call lsquit('ls_mpi_buffer_integer: error using buffer',-1)
         buffer = lsmpibufferInt8(iInt8+1)
         iInt8 = iInt8 + 1
      ENDIF
#endif
    end subroutine ls_mpi_buffer8_integer
    
    subroutine ls_mpi_buffer_integer8V_buf4_wrapper(buffer,nbuf,master)
      implicit none
      integer(kind=ls_mpik):: master
      integer(kind=4) :: nbuf
      integer(kind=8) :: buffer(:)
#ifdef VAR_MPI
      integer(kind=8) :: nbuf8
      nbuf8=nbuf
      call ls_mpi_buffer_integer8V(buffer,nbuf8,master)
#endif
    end subroutine ls_mpi_buffer_integer8V_buf4_wrapper

    subroutine ls_mpi_buffer_integer8V(buffer,nbuf,master)
      implicit none
      integer(kind=ls_mpik) :: master
      integer(kind=8) :: nbuf
      integer(kind=8) :: buffer(:)
#ifdef VAR_MPI
      integer :: I
      IF(AddToBuffer)THEN
         IF(iInt8+nbuf .GT. nInteger8) THEN
           call increaselsmpibufferInt8(nbuf)
         ENDIF
         DO I=1,nbuf
            lsmpibufferInt8(iInt8+I) = buffer(I)
         ENDDO
         iInt8 = iInt8 + nbuf
      ELSE
         IF(iInt8+nbuf .GT. nInteger8) THEN
           write(*,*) 'ls_mpi_buffer_integer8V:',infpar%mynum,iInt8,nbuf,nInteger8
           call lsquit('ls_mpi_buffer_integerV: error using buffer',-1)
         ENDIF
         DO I=1,nbuf
            buffer(I) = lsmpibufferInt8(iInt8+I)
         ENDDO
         iInt8 = iInt8 + nbuf
      ENDIF
#endif
    end subroutine ls_mpi_buffer_integer8V

    subroutine ls_mpi_buffer_integer4V_buf4_wrapper(buffer,nbuf,master)
      implicit none
      integer(kind=ls_mpik) :: master
      integer(kind=4) :: nbuf
      integer(kind=4) :: buffer(:)
#ifdef VAR_MPI
      integer(kind=8) :: nbuf8
      nbuf8=nbuf
      call ls_mpi_buffer_integer4V(buffer,nbuf8,master)
#endif
    end subroutine ls_mpi_buffer_integer4V_buf4_wrapper
    subroutine ls_mpi_buffer_integer4V(buffer,nbuf,master)
      implicit none
      integer(kind=ls_mpik) :: master
      integer(kind=8) :: nbuf
      integer(kind=4) :: buffer(:)
#ifdef VAR_MPI
      integer :: I
      IF(AddToBuffer)THEN
         IF(iInt4+nbuf .GT. nInteger4) THEN
           call increaselsmpibufferInt4(nbuf)
         ENDIF
         DO I=1,nbuf
            IF(iInt4+1.GT.size(lsmpibufferInt4,kind=long))call lsquit('errorTK',-1)
            lsmpibufferInt4(iInt4+I) = buffer(I)
         ENDDO
         iInt4 = iInt4 + nbuf
      ELSE
         IF(iInt4+nbuf .GT. nInteger4) THEN
           write(*,*) 'ls_mpi_buffer_integer4V:',infpar%mynum,iInt4,nbuf,nInteger4
           call lsquit('ls_mpi_buffer_integer4V: error using buffer',-1)
         ENDIF
         DO I=1,nbuf
            buffer(I) = lsmpibufferInt4(iInt4+I)
         ENDDO
         iInt4 = iInt4 + nbuf
      ENDIF
#endif
    end subroutine ls_mpi_buffer_integer4V

    subroutine ls_mpi_buffer_integer4M_wrapper4(buffer,nbuf1,nbuf2,master)
      implicit none
      integer(kind=ls_mpik) :: master
      integer(kind=4) :: nbuf1,nbuf2
      integer(kind=4) :: buffer(nbuf1,nbuf2)
#ifdef VAR_MPI
      integer(kind=4),pointer :: buf(:)
      integer(kind=8) :: nbuf
      nbuf=nbuf1*nbuf2
      call ass_44I2to1(buffer,buf,[nbuf1,nbuf2])
      call ls_mpi_buffer_integer4V(buf,nbuf,master)
      nullify(buf)
#endif
    end subroutine ls_mpi_buffer_integer4M_wrapper4
    subroutine ls_mpi_buffer_integer4M_wrapper8(buffer,nbuf1,nbuf2,master)
      implicit none
      integer(kind=ls_mpik) :: master
      integer(kind=8) :: nbuf1,nbuf2
      integer(kind=4) :: buffer(nbuf1*nbuf2)
#ifdef VAR_MPI
      integer(kind=4),pointer :: buf(:)
      integer(kind=8) :: nbuf
      nbuf=nbuf1*nbuf2
      call ass_48I2to1(buffer,buf,[nbuf1,nbuf2])
      call ls_mpi_buffer_integer4V(buf,nbuf,master)
      nullify(buf)
#endif
    end subroutine ls_mpi_buffer_integer4M_wrapper8
    subroutine ls_mpi_buffer_integer8M_wrapper4(buffer,nbuf1,nbuf2,master)
      implicit none
      integer(kind=ls_mpik) :: master
      integer(kind=4) :: nbuf1,nbuf2
      integer(kind=8) :: buffer(nbuf1,nbuf2)
#ifdef VAR_MPI
      integer(kind=8),pointer :: buf(:)
      integer(kind=8) :: nbuf
      nbuf=nbuf1*nbuf2
      call ass_84I2to1(buffer,buf,[nbuf1,nbuf2])
      call ls_mpi_buffer_integer8V(buf,nbuf,master)
      nullify(buf)
#endif
    end subroutine ls_mpi_buffer_integer8M_wrapper4
    subroutine ls_mpi_buffer_integer8M_wrapper8(buffer,nbuf1,nbuf2,master)
      implicit none
      integer(kind=ls_mpik) :: master
      integer(kind=8) :: nbuf1,nbuf2
      integer(kind=8) :: buffer(nbuf1,nbuf2)
#ifdef VAR_MPI
      integer(kind=8),pointer :: buf(:)
      integer(kind=8) :: nbuf
      nbuf=nbuf1*nbuf2
      call ass_88I2to1(buffer,buf,[nbuf1,nbuf2])
      call ls_mpi_buffer_integer8V(buf,nbuf,master)
      nullify(buf)
#endif
    end subroutine ls_mpi_buffer_integer8M_wrapper8


    subroutine ls_mpi_buffer_integer4T_wrapper4(buffer,nbuf1,nbuf2,nbuf3,master)
      implicit none
      integer(kind=ls_mpik) :: master
      integer(kind=4) :: nbuf1,nbuf2,nbuf3
      integer(kind=4) :: buffer(nbuf1,nbuf2,nbuf3)
#ifdef VAR_MPI
      integer(kind=4),pointer :: buf(:)
      integer(kind=8) :: nbuf
      nbuf=nbuf1*nbuf2*nbuf3
      call ass_44I3to1(buffer,buf,[nbuf1,nbuf2,nbuf3])
      call ls_mpi_buffer_integer4V(buf,nbuf,master)
      nullify(buf)
#endif
    end subroutine ls_mpi_buffer_integer4T_wrapper4
    subroutine ls_mpi_buffer_integer4T_wrapper8(buffer,nbuf1,nbuf2,nbuf3,master)
      implicit none
      integer(kind=ls_mpik) :: master
      integer(kind=8) :: nbuf1,nbuf2,nbuf3
      integer(kind=4) :: buffer(nbuf1,nbuf2,nbuf3)
#ifdef VAR_MPI
      integer(kind=4),pointer :: buf(:)
      integer(kind=8) :: nbuf
      nbuf=nbuf1*nbuf2*nbuf3
      call ass_48I3to1(buffer,buf,[nbuf1,nbuf2,nbuf3])
      call ls_mpi_buffer_integer4V(buf,nbuf,master)
      nullify(buf)
#endif
    end subroutine ls_mpi_buffer_integer4T_wrapper8
    subroutine ls_mpi_buffer_integer8T_wrapper4(buffer,nbuf1,nbuf2,nbuf3,master)
      implicit none
      integer(kind=ls_mpik) :: master
      integer(kind=4) :: nbuf1,nbuf2,nbuf3
      integer(kind=8) :: buffer(nbuf1,nbuf2,nbuf3)
#ifdef VAR_MPI
      integer(kind=8),pointer :: buf(:)
      integer(kind=8) :: nbuf
      nbuf=nbuf1*nbuf2*nbuf3
      call ass_84I3to1(buffer,buf,[nbuf1,nbuf2,nbuf3])
      call ls_mpi_buffer_integer8V(buf,nbuf,master)
      nullify(buf)
#endif
    end subroutine ls_mpi_buffer_integer8T_wrapper4
    subroutine ls_mpi_buffer_integer8T_wrapper8(buffer,nbuf1,nbuf2,nbuf3,master)
      implicit none
      integer(kind=ls_mpik) :: master
      integer(kind=8) :: nbuf1,nbuf2,nbuf3
      integer(kind=8) :: buffer(nbuf1,nbuf2,nbuf3)
#ifdef VAR_MPI
      integer(kind=8),pointer :: buf(:)
      integer(kind=8) :: nbuf
      nbuf=nbuf1*nbuf2*nbuf3
      call ass_88I3to1(buffer,buf,[nbuf1,nbuf2,nbuf3])
      call ls_mpi_buffer_integer8V(buf,nbuf,master)
      nullify(buf)
#endif
    end subroutine ls_mpi_buffer_integer8T_wrapper8


    subroutine ls_mpi_buffer_integer4Q_wrapper4(buffer,nbuf1,nbuf2,nbuf3,nbuf4,master)
      implicit none
      integer(kind=ls_mpik) :: master
      integer(kind=4) :: nbuf1,nbuf2,nbuf3,nbuf4
      integer(kind=4) :: buffer(nbuf1,nbuf2,nbuf3,nbuf4)
#ifdef VAR_MPI
      integer(kind=4),pointer :: buf(:)
      integer(kind=8) :: nbuf
      nbuf=nbuf1*nbuf2*nbuf3*nbuf4
      call ass_44I4to1(buffer,buf,[nbuf1,nbuf2,nbuf3,nbuf4])
      call ls_mpi_buffer_integer4V(buf,nbuf,master)
      nullify(buf)
#endif
    end subroutine ls_mpi_buffer_integer4Q_wrapper4
    subroutine ls_mpi_buffer_integer4Q_wrapper8(buffer,nbuf1,nbuf2,nbuf3,nbuf4,master)
      implicit none
      integer(kind=ls_mpik) :: master
      integer(kind=8) :: nbuf1,nbuf2,nbuf3,nbuf4
      integer(kind=4) :: buffer(nbuf1,nbuf2,nbuf3,nbuf4)
#ifdef VAR_MPI
      integer(kind=4),pointer :: buf(:)
      integer(kind=8) :: nbuf
      nbuf=nbuf1*nbuf2*nbuf3*nbuf4
      call ass_48I4to1(buffer,buf,[nbuf1,nbuf2,nbuf3,nbuf4])
      call ls_mpi_buffer_integer4V(buf,nbuf,master)
      nullify(buf)
#endif
    end subroutine ls_mpi_buffer_integer4Q_wrapper8
    subroutine ls_mpi_buffer_integer8Q_wrapper4(buffer,nbuf1,nbuf2,nbuf3,nbuf4,master)
      implicit none
      integer(kind=ls_mpik) :: master
      integer(kind=4) :: nbuf1,nbuf2,nbuf3,nbuf4
      integer(kind=8) :: buffer(nbuf1,nbuf2,nbuf3,nbuf4)
#ifdef VAR_MPI
      integer(kind=8),pointer :: buf(:)
      integer(kind=8) :: nbuf
      nbuf=nbuf1*nbuf2*nbuf3*nbuf4
      call ass_84I4to1(buffer,buf,[nbuf1,nbuf2,nbuf3,nbuf4])
      call ls_mpi_buffer_integer8V(buf,nbuf,master)
      nullify(buf)
#endif
    end subroutine ls_mpi_buffer_integer8Q_wrapper4
    subroutine ls_mpi_buffer_integer8Q_wrapper8(buffer,nbuf1,nbuf2,nbuf3,nbuf4,master)
      implicit none
      integer(kind=ls_mpik) :: master
      integer(kind=8) :: nbuf1,nbuf2,nbuf3,nbuf4
      integer(kind=8) :: buffer(nbuf1,nbuf2,nbuf3,nbuf4)
#ifdef VAR_MPI
      integer(kind=8),pointer :: buf(:)
      integer(kind=8) :: nbuf
      nbuf=nbuf1*nbuf2*nbuf3*nbuf4
      call ass_88I4to1(buffer,buf,[nbuf1,nbuf2,nbuf3,nbuf4])
      call ls_mpi_buffer_integer8V(buf,nbuf,master)
      nullify(buf)
#endif
    end subroutine ls_mpi_buffer_integer8Q_wrapper8


    subroutine ls_mpi_buffer_realk(buffer,master)
      implicit none
      real(realk) :: buffer
      integer(kind=ls_mpik) :: master
#ifdef VAR_MPI
      IF(AddToBuffer)THEN
         IF(iDP+1 .GT. nDP)call increaselsmpibufferDP(1)
         lsmpibufferDP(iDP+1) = buffer
         iDP = iDP + 1
      ELSE
         IF(iDP+1 .GT. nDP)call lsquit('ls_mpi_buffer_realk: error using buffer',-1)
         buffer = lsmpibufferDP(iDP+1)
         iDP = iDP + 1
      ENDIF
#endif
    end subroutine ls_mpi_buffer_realk

    subroutine ls_mpi_buffer_realkV(buffer,nbuf,master)
      implicit none
      integer :: nbuf
      integer(kind=ls_mpik) :: master
      real(realk) :: buffer(:)
#ifdef VAR_MPI
      integer :: I
      IF(AddToBuffer)THEN
         IF(iDP + nbuf.GT. nDP)call increaselsmpibufferDP(nbuf)
         DO I=1,nbuf
            lsmpibufferDP(iDP+I) = buffer(I)
         ENDDO
         iDP = iDP + nbuf
      ELSE
         IF(iDP+nbuf .GT. nDP)call lsquit('ls_mpi_buffer_realkV: error using buffer',-1)
         DO I=1,nbuf
            buffer(I) = lsmpibufferDP(iDP+I)
         ENDDO
         iDP = iDP + nbuf
      ENDIF
#endif
    end subroutine ls_mpi_buffer_realkV

    subroutine ls_mpi_buffer_realkM(buffer,nbuf1,nbuf2,master)
      implicit none
      integer :: nbuf1,nbuf2
      integer(kind=ls_mpik) :: master
      real(realk) :: buffer(:,:)
#ifdef VAR_MPI
      integer :: ierr,count,datatype,I,J,offset
      IERR=0
      IF(AddToBuffer)THEN
         IF(iDP + nbuf1*nbuf2.GT. nDP)call increaselsmpibufferDP(nbuf1*nbuf2)
         DO J=1,nbuf2
            offset = iDP+(J-1)*nbuf1
            DO I=1,nbuf1
               lsmpibufferDP(offset+I) = buffer(I,J)
            ENDDO
         ENDDO
         iDP = iDP + nbuf1*nbuf2
      ELSE
         IF(iDP+nbuf1*nbuf2 .GT. nDP)call lsquit('ls_mpi_buffer_realkM: error using buffer',-1)
         DO J=1,nbuf2
            offset = iDP+(J-1)*nbuf1
            DO I=1,nbuf1
               buffer(I,J) = lsmpibufferDP(offset+I)
            ENDDO
         ENDDO
         iDP = iDP + nbuf1*nbuf2
      ENDIF
#endif
    end subroutine ls_mpi_buffer_realkM

    subroutine ls_mpi_buffer_realkT(buffer,nbuf1,nbuf2,nbuf3,master)
      implicit none
      integer(kind=ls_mpik) :: master
      integer :: nbuf1,nbuf2,nbuf3
      real(realk) :: buffer(:,:,:)
#ifdef VAR_MPI
      integer :: ierr,count,datatype,I,J,K,offset
      IERR=0
      IF(AddToBuffer)THEN
         IF(iDP + nbuf1*nbuf2*nbuf3.GT. nDP)call increaselsmpibufferDP(nbuf1*nbuf2*nbuf3)
         DO K=1,nbuf3
            DO J=1,nbuf2
               offset = iDP+(J-1)*nbuf1+(K-1)*nbuf1*nbuf2
               DO I=1,nbuf1
                  lsmpibufferDP(offset+I) = buffer(I,J,K)
               ENDDO
            ENDDO
         ENDDO
         iDP = iDP + nbuf1*nbuf2*nbuf3
      ELSE
         IF(iDP+nbuf1*nbuf2*nbuf3 .GT. nDP)call lsquit('ls_mpi_buffer_realkT: error using buffer',-1)
         DO K=1,nbuf3
            DO J=1,nbuf2
               offset = iDP+(J-1)*nbuf1+(K-1)*nbuf1*nbuf2
               DO I=1,nbuf1
                  buffer(I,J,K) = lsmpibufferDP(offset+I)
               ENDDO
            ENDDO
         ENDDO
         iDP = iDP + nbuf1*nbuf2*nbuf3
      ENDIF
#endif
    end subroutine ls_mpi_buffer_realkT

    subroutine ls_mpi_buffer_logical(buffer,master)
      implicit none
      logical :: buffer
      integer(kind=ls_mpik) :: master
#ifdef VAR_MPI
      integer :: ierr,count,datatype
      IERR=0
      IF(AddToBuffer)THEN
         IF(iLog+1 .GT. nLog)call increaselsmpibufferLog(1)
         lsmpibufferLog(iLog+1) = buffer
         iLog = iLog + 1
      ELSE
         IF(iLog+1 .GT. nLog)call lsquit('ls_mpi_buffer_logical: error using buffer',-1)
         buffer = lsmpibufferLog(iLog+1)
         iLog = iLog + 1
      ENDIF
#endif
    end subroutine ls_mpi_buffer_logical

    subroutine ls_mpi_buffer_logicalV(buffer,nbuf,master)
      implicit none
      integer :: nbuf
      integer(kind=ls_mpik) :: master
      logical :: buffer(:)
#ifdef VAR_MPI
      integer :: ierr,count,datatype,I
      IERR=0
      IF(AddToBuffer)THEN
         IF(iLog + nbuf.GT. nLog)call increaselsmpibufferLog(nbuf)
         DO I=1,nbuf
            lsmpibufferLog(iLog+I) = buffer(I)
         ENDDO
         iLog = iLog + nbuf
      ELSE
         IF(iLog+nbuf .GT. nLog)call lsquit('ls_mpi_buffer_logicalV: error using buffer',-1)
         DO I=1,nbuf
            buffer(I) = lsmpibufferLog(iLog+I)
         ENDDO
         iLog = iLog + nbuf
      ENDIF
#endif
    end subroutine ls_mpi_buffer_logicalV

    subroutine ls_mpi_buffer_logicalM(buffer,nbuf1,nbuf2,master)
      implicit none
      integer :: nbuf1,nbuf2
      integer(kind=ls_mpik) :: master
      logical :: buffer(:,:)
#ifdef VAR_MPI
      integer :: ierr,count,datatype,I,J,offset
      IERR=0
      IF(AddToBuffer)THEN
         IF(iLog + nbuf1*nbuf2.GT. nLog)call increaselsmpibufferLog(nbuf1*nbuf2)
         DO J=1,nbuf2
            offset = iLog+(J-1)*nbuf1
            DO I=1,nbuf1
               lsmpibufferLog(offset+I) = buffer(I,J)
            ENDDO
         ENDDO
         iLog = iLog + nbuf1*nbuf2
      ELSE
         IF(iLog+nbuf1*nbuf2 .GT. nLog)call lsquit('ls_mpi_buffer_logicalM: error using buffer',-1)
         DO J=1,nbuf2
            offset = iLog+(J-1)*nbuf1
            DO I=1,nbuf1
               buffer(I,J) = lsmpibufferLog(offset+I)
            ENDDO
         ENDDO
         iLog = iLog + nbuf1*nbuf2
      ENDIF
#endif
    end subroutine ls_mpi_buffer_logicalM

    subroutine ls_mpi_buffer_charac(buffer,master)
      implicit none
      character :: buffer
      integer(kind=ls_mpik) :: master
#ifdef VAR_MPI
      integer :: ierr,count,datatype
      IERR=0
      IF(AddToBuffer)THEN
         IF(iCha+1 .GT. nCha)call increaselsmpibufferCha(1)
         lsmpibufferCha(iCha+1) = buffer
         iCha = iCha + 1
      ELSE
         IF(iCha+1 .GT. nCha)call lsquit('ls_mpi_buffer_charac: error using buffer',-1)
         buffer = lsmpibufferCha(iCha+1)
         iCha = iCha + 1
      ENDIF
#endif
    end subroutine ls_mpi_buffer_charac

    subroutine ls_mpi_buffer_characV(buffer,nbuf,master)
      implicit none
      integer :: nbuf
      integer(kind=ls_mpik) :: master
      character*(*) :: buffer
#ifdef VAR_MPI
      character :: tmpbuffer
      integer :: ierr,count,datatype,I
      IERR=0
      IF(AddToBuffer)THEN
         IF(iCha + nbuf.GT. nCha)call increaselsmpibufferCha(nbuf)
         DO I=1,nbuf
            lsmpibufferCha(iCha+I) = buffer(I:I)
         ENDDO
         iCha = iCha + nbuf
      ELSE
         IF(iCha+nbuf .GT. nCha)call lsquit('ls_mpi_buffer_characV: error using buffer',-1)
         DO I=1,nbuf
            buffer(I:I) = lsmpibufferCha(iCha+I)
         ENDDO
         iCha = iCha + nbuf
      ENDIF
#endif
    end subroutine ls_mpi_buffer_characV

    subroutine ls_mpi_buffer_characV2(buffer,nbuf,master)
      implicit none
      integer :: nbuf
      integer(kind=ls_mpik) :: master
      character :: buffer(nbuf)
#ifdef VAR_MPI
      character :: tmpbuffer
      integer :: ierr,count,datatype,I
      IERR=0
      IF(AddToBuffer)THEN
         IF(iCha + nbuf.GT. nCha)call increaselsmpibufferCha(nbuf)
         DO I=1,nbuf
            lsmpibufferCha(iCha+I) = buffer(I)
         ENDDO
         iCha = iCha + nbuf
      ELSE
         IF(iCha+nbuf .GT. nCha)call lsquit('ls_mpi_buffer_characV2: error using buffer',-1)
         DO I=1,nbuf
            buffer(I) = lsmpibufferCha(iCha+I)
         ENDDO
         iCha = iCha + nbuf
      ENDIF
#endif
    end subroutine ls_mpi_buffer_characV2

    subroutine ls_mpi_buffer_shortinteger(buffer,master)

      implicit none
      integer(kind=ls_mpik) :: master
      integer(kind=short) :: buffer
#ifdef VAR_MPI
      IF(AddToBuffer)THEN
         IF(iSho+1 .GT. nShort)call increaselsmpibufferSho(1)
         lsmpibufferSho(iSho+1) = buffer
         iSho = iSho + 1
      ELSE
         IF(iSho+1 .GT. nShort) then
           write(*,*) 'ls_mpi_buffer_shortinteger',iSho,nShort
           call lsquit('ls_mpi_buffer_shortinteger: error using buffer',-1)
         ENDIF
         buffer = lsmpibufferSho(iSho+1)
         iSho = iSho + 1
      ENDIF
#endif
    end subroutine ls_mpi_buffer_shortinteger

    subroutine ls_mpi_buffer_shortintegerV(buffer,nbuf,master)
      implicit none
      integer :: nbuf
      integer(kind=ls_mpik) :: master
      integer(kind=short) :: buffer(:)
#ifdef VAR_MPI
      integer :: I

      IF(AddToBuffer)THEN
         IF(iSho+nbuf .GT. nShort)call increaselsmpibufferSho(nbuf)
         DO I=1,nbuf
            lsmpibufferSho(iSho+I) = buffer(I)
         ENDDO
         iSho = iSho + nbuf
      ELSE
         IF(iSho+nbuf .GT. nShort) then
           write(*,*) 'ls_mpi_buffer_shortintegerV',iSho,nbuf,nShort
           call lsquit('ls_mpi_buffer_shortintegerV: error using buffer',-1)
         ENDIF
         DO I=1,nbuf
            buffer(I) = lsmpibufferSho(iSho+I)
         ENDDO
         iSho = iSho + nbuf
      ENDIF
#endif
    end subroutine ls_mpi_buffer_shortintegerV

    subroutine ls_mpi_buffer_shortintegerM(buffer,nbuf1,nbuf2,master)
      implicit none
      integer :: nbuf1,nbuf2
      integer(kind=ls_mpik) :: master
      integer(kind=short) :: buffer(:,:)
#ifdef VAR_MPI
      integer :: I,J,offset

      IF(AddToBuffer)THEN
         IF(iSho+nbuf1*nbuf2 .GT. nShort)call increaselsmpibufferSho(nbuf1*nbuf2)
         DO J=1,nbuf2
            offset = iSho+(J-1)*nbuf1
            DO I=1,nbuf1
               lsmpibufferSho(offset+I) = buffer(I,J)
            ENDDO
         ENDDO
         iSho = iSho + nbuf1*nbuf2
      ELSE
         IF(iSho+nbuf1*nbuf2 .GT. nShort) then
           write(*,*) 'ls_mpi_buffer_shortintegerM',iSho,nbuf1,nbuf2,nShort
           call lsquit('ls_mpi_buffer_shortintegerM: error using buffer',-1)
         ENDIF
         DO J=1,nbuf2
            offset = iSho+(J-1)*nbuf1
            DO I=1,nbuf1
               buffer(I,J) = lsmpibufferSho(offset+I)
            ENDDO
         ENDDO
         iSho = iSho + nbuf1*nbuf2
      ENDIF
#endif
    end subroutine ls_mpi_buffer_shortintegerM

!######################################################
!#
!#
!#
!######################################################
    subroutine ls_mpiInitBuffer(master,Job,comm,sender,receiver)
      implicit none
      integer(kind=ls_mpik) :: master
      integer :: Job
      integer(kind=ls_mpik) :: comm  ! communicator
      !> Only for Job=LSMPISENDRECV: rank for sender within comm group
      integer(kind=ls_mpik),intent(in),optional :: sender
      !> Only for Job=LSMPISENDRECV: rank for receiver within comm group
      integer(kind=ls_mpik),intent(in),optional :: receiver
#ifdef VAR_MPI
      integer(kind=8) :: ndim(6)
      integer(kind=ls_mpik) :: mynum  ! Number of node WITHIN group specified by communicator

      ! Get rank within specific communicator
      call get_rank_for_comm(comm,mynum)
      AddToBuffer = (mynum.EQ.master) .OR. (Job.EQ.LSMPIREDUCTION)

      ! Special settings for send/receive job 
      if(job.eq.LSMPISENDRECV) then

         ! Sanity check
         if(.not. present(sender)) call lsquit('ls_mpiInitBuffer: Sender not present!',-1)
         if(.not. present(receiver)) call lsquit('ls_mpiInitBuffer: Receiver not present!',-1)

         ! Set addtobuffer and check that current rank is meaningful
         IF(mynum.EQ.sender)THEN
            addtobuffer=.true.
         else if(mynum.EQ.receiver) then
            addtobuffer=.false.
         else
            print '(a,3i6)', 'Rank,sender,receiver',mynum,sender,receiver
            call lsquit('ls_mpiInitBuffer: Rank number is neither sender nor receiver!',-1)
         end if
      end if

      ! Nullify buffers
      call nullify_mpibuffer


      IF( (Job.EQ.LSMPIBROADCAST) .or. (Job.EQ.LSMPISENDRECV)) THEN
         
         ! job .eq. lsmpibroadcast
         ! -----------------------
         ! Bcast buffer from master to slaves
         ! If master: Put stuff into buffer
         ! If slave: Read stuff from buffer

         ! job .eq. lsmpisendrecv
         ! If sender: Put stuff into buffer
         ! If receiver: Receive stuff from buffer

         if(addtobuffer) then  ! put stuff info buffer
            IF(nLog.EQ.0) nLog = incremLog
            IF(nDP.EQ.0) nDP = incremDP
            IF(nInteger4.EQ.0) nInteger4 = incremInteger
            IF(nInteger8.EQ.0) nInteger8 = incremInteger
            IF(nShort.EQ.0) nShort = incremShort
            IF(nCha.EQ.0) nCha = incremCha
            IF(MemModParamPrintMemory)THEN
               print*,'MemModParamPrintMemory   mynum',mynum
               CALL Print_Memory_info(MemModParamPrintMemorylupri,&
                    & 'InitBuffer: Before allocation in AddtoBuffer')
               Write(MemModParamPrintMemorylupri,*)'# DoublePrecison elements',nDP
               Write(MemModParamPrintMemorylupri,*)'# Integer 4 elements     ',nInteger4
               Write(MemModParamPrintMemorylupri,*)'# Integer 8 elements     ',nInteger8
               Write(MemModParamPrintMemorylupri,*)'# Short integer elements ',nShort
               Write(MemModParamPrintMemorylupri,*)'# Logical elements       ',nLog
               Write(MemModParamPrintMemorylupri,*)'# Character elements     ',nCha
            ENDIF
            call mem_alloc(lsmpibufferDP,nDP)
            call mem_alloc(lsmpibufferInt4,nInteger4)
            call mem_alloc(lsmpibufferInt8,nInteger8)
            call mem_alloc(lsmpibufferSho,nShort)
            call mem_alloc(lsmpibufferLog,nLog)
            call mem_alloc(lsmpibufferCha,nCha)
             IF(MemModParamPrintMemory)THEN
               CALL Print_Memory_info(MemModParamPrintMemorylupri,&
                    &'InitBuffer: After allocation in AddtoBuffer')
            ENDIF
           iDP   = 0
            iInt4 = 0
            iInt8 = 0
            iSho  = 0
            iLog  = 0
            iCha  = 0
         ELSE ! read stuff from buffer
            if(job.eq.LSMPIBROADCAST) then  ! bcast
               call ls_mpibcast(ndim,6,master,COMM)
            else ! specific receiver
               call ls_mpisendrecv(ndim,6,comm,sender,receiver)
            end if
            nDP = ndim(1)
            nInteger4 = ndim(2)
            nInteger8 = ndim(3)
            nShort = ndim(4)
            nLog = ndim(5)
            nCha = ndim(6)
            IF(MemModParamPrintMemory)THEN
               print*,'MemModParamPrintMemory   mynum',mynum
               CALL Print_Memory_info(MemModParamPrintMemorylupri,'InitBuffer: Before allocation in ReadFromBuffer')
               Write(MemModParamPrintMemorylupri,*)'# DoublePrecison elements',nDP
               Write(MemModParamPrintMemorylupri,*)'# Integer 4 elements     ',nInteger4
               Write(MemModParamPrintMemorylupri,*)'# Integer 8 elements     ',nInteger8
               Write(MemModParamPrintMemorylupri,*)'# Short integer elements ',nShort
               Write(MemModParamPrintMemorylupri,*)'# Logical elements       ',nLog
               Write(MemModParamPrintMemorylupri,*)'# Character elements     ',nCha
            ENDIF
            if(ndp.gt.0) call mem_alloc(lsmpibufferDP,nDP)
            if(ninteger4.gt.0) call mem_alloc(lsmpibufferInt4,nInteger4)
            if(ninteger8.gt.0) call mem_alloc(lsmpibufferInt8,nInteger8)
            if(nshort .gt. 0) call mem_alloc(lsmpibufferSho,nShort)
            if(nlog.gt.0) call mem_alloc(lsmpibufferLog,nLog)
            if(ncha .gt. 0) call mem_alloc(lsmpibufferCha,nCha)
            IF(MemModParamPrintMemory)THEN
               CALL Print_Memory_info(MemModParamPrintMemorylupri,'InitBuffer: After allocation in ReadFromBuffer')
            ENDIF
            if(job.eq.lsmpibroadcast) then  ! BCAST
               IF(ndim(1).GT.0)THEN
                  call ls_mpibcast(lsmpibufferDP,nDP,master,COMM)
               ENDIF
               IF(ndim(2).GT.0)THEN
                  call ls_mpibcast(lsmpibufferInt4,nInteger4,master,COMM)
               ENDIF
               IF(ndim(3).GT.0)THEN
                  call ls_mpibcast(lsmpibufferInt8,nInteger8,master,COMM)
               ENDIF
               IF(ndim(4).GT.0)THEN
                  call ls_mpibcast(lsmpibufferSho,nShort,master,COMM)
               ENDIF
               IF(ndim(5).GT.0)THEN
                  call ls_mpibcast(lsmpibufferLog,nLog,master,COMM)
               ENDIF
               IF(ndim(6).GT.0)THEN
                  call ls_mpibcast(lsmpibufferCha,nCha,master,COMM)
               ENDIF

            else  ! Receiver: Receive job from sender (sent in ls_MpiFinalizeBuffer)
               IF(ndim(1).GT.0)THEN
                  call ls_mpisendrecv(lsmpibufferDP,nDP,comm,sender,receiver)
               ENDIF
               IF(ndim(2).GT.0)THEN
                  call ls_mpisendrecv(lsmpibufferInt4,nInteger4,comm,sender,receiver)
               ENDIF
               IF(ndim(3).GT.0)THEN
                  call ls_mpisendrecv(lsmpibufferInt8,nInteger8,comm,sender,receiver)
               ENDIF
               IF(ndim(4).GT.0)THEN
                  call ls_mpisendrecv(lsmpibufferSho,nShort,comm,sender,receiver)
               ENDIF
               IF(ndim(5).GT.0)THEN
                  call ls_mpisendrecv(lsmpibufferLog,nLog,comm,sender,receiver)
               ENDIF
               IF(ndim(6).GT.0)THEN
                  call ls_mpisendrecv(lsmpibufferCha,nCha,comm,sender,receiver)
               ENDIF
            end if
               
            iDP = 0
            iInt4 = 0
            iInt8 = 0
            iSho = 0
            iLog = 0
            iCha = 0
         ENDIF
      else if(Job.EQ.LSMPIREDUCTION)THEN
         IF(nLog.EQ.0) nLog = incremLog+1
         IF(nDP.EQ.0) nDP = incremDP+1
         IF(nInteger4.EQ.0) nInteger4 = incremInteger+1
         IF(nInteger8.EQ.0) nInteger8 = incremInteger+1
         IF(nShort.EQ.0) nShort = incremShort+int_to_short
         IF(nCha.EQ.0) nCha = incremCha+1
         IF(MemModParamPrintMemory)THEN
            print*,'MemModParamPrintMemory   mynum',mynum
            CALL Print_Memory_info(MemModParamPrintMemorylupri,'InitBuffer: Before allocation in LSMPIREDUCTION')
            Write(MemModParamPrintMemorylupri,*)'# DoublePrecison elements',nDP
            Write(MemModParamPrintMemorylupri,*)'# Integer 4 elements     ',nInteger4
            Write(MemModParamPrintMemorylupri,*)'# Integer 8 elements     ',nInteger8
            Write(MemModParamPrintMemorylupri,*)'# Short integer elements ',nShort
            Write(MemModParamPrintMemorylupri,*)'# Logical elements       ',nLog
            Write(MemModParamPrintMemorylupri,*)'# Character elements     ',nCha
         ENDIF
         call mem_alloc(lsmpibufferDP,nDP)
         call mem_alloc(lsmpibufferInt4,nInteger4)
         call mem_alloc(lsmpibufferInt8,nInteger8)
         call mem_alloc(lsmpibufferSho,nShort)
         call mem_alloc(lsmpibufferLog,nLog)
         call mem_alloc(lsmpibufferCha,nCha)
         IF(MemModParamPrintMemory)THEN
            CALL Print_Memory_info(MemModParamPrintMemorylupri,'InitBuffer: After allocation in ReadFromBuffer')
         ENDIF
         iDP = 0
         iInt4 = 0
         iInt8 = 0
         iSho = 0
         iLog = 0
         iCha = 0
      ELSE
         call lsquit('Job specifier invalid in ls_mpiInitBuffer',-1)
      ENDIF
#endif

    end subroutine ls_mpiInitBuffer

    subroutine ls_MpiFinalizeBuffer(master,Job,comm,sender,receiver)
      implicit none
      integer(kind=ls_mpik) :: master
      integer :: Job
      integer(kind=ls_mpik) :: comm  ! communicator
      !> Only for Job=LSMPISENDRECV: rank for sender
      integer(kind=ls_mpik),intent(in),optional :: sender
      !> Only for Job=LSMPISENDRECV: rank for receiver
      integer(kind=ls_mpik),intent(in),optional :: receiver
#ifdef VAR_MPI
      integer(kind=8) :: ndim(6)
      integer :: II,additional,modula
      integer(kind=ls_mpik) :: mynum  ! Number of node WITHIN group specified by communicator

      ! Get rank number within specific communicator
      call get_rank_for_comm(comm,mynum)

      ! Special checks for send/receive job 
      if(job.eq.LSMPISENDRECV) then

         ! Sanity check
         if(.not. present(sender)) call lsquit('ls_mpiInitBuffer: Sender not present!',-1)
         if(.not. present(receiver)) call lsquit('ls_mpiInitBuffer: Receiver not present!',-1)

         ! Set addtobuffer and check that current rank is meaningful
         IF(mynum/=sender .and. mynum/=receiver) THEN
            print '(a,3i6)', 'Rank,sender,receiver',mynum,sender,receiver
            call lsquit('ls_mpiInitBuffer: Rank number is neither sender nor receiver!',-1)
         end if
      end if


      IF( (Job.EQ.LSMPIBROADCAST) .or. (Job.EQ.LSMPISENDRECV)) THEN
         ! If addtobuffer was set to TRUE in ls_mpiInitBuffer, then
         ! the current rank should send buffer info here.
         ! If addtobuffer was set to FALSE in ls_mpiInitBuffer, then
         ! the current rank received in ls_mpiInitBuffer, and 
         ! the buffer should just be deallocd again here.

         IF(AddToBuffer) THEN
            ndim(1) = iDP
            ndim(2) = iInt4
            ndim(3) = iInt8
! Add modula of iSho divided by int_to_short in order to use integer mpibcast for short integer
            modula = mod(iSho,int_to_short)
            additional = 0
            IF (modula.GT.0) additional = int_to_short - modula
            ndim(4) = iSho + additional
            ndim(5) = iLog
            ndim(6) = iCha

            nDP = iDP
            nInteger4 = iInt4
            nInteger8 = iInt8
            nShort = iSho + additional
            nLoG = iLog
            nCha = iCha
            if(job .eq. LSMPIBROADCAST) then  ! communication via bcast
               call ls_mpibcast(ndim,6,master,COMM)
               IF(ndim(1).GT.0)THEN
                  call ls_mpibcast(lsmpibufferDP(1:ndim(1)),ndim(1),master,COMM)
               ENDIF
               IF(ndim(2).GT.0)THEN
                  call ls_mpibcast(lsmpibufferInt4(1:ndim(2)),ndim(2),master,COMM)
               ENDIF
               IF(ndim(3).GT.0)THEN
                  call ls_mpibcast(lsmpibufferInt8(1:ndim(3)),ndim(3),master,COMM)
               ENDIF
               IF(ndim(4).GT.0)THEN
                  call ls_mpibcast(lsmpibufferSho(1:ndim(4)),ndim(4),master,COMM)
               ENDIF
               IF(ndim(5).GT.0)THEN
                  call ls_mpibcast(lsmpibufferLog(1:ndim(5)),ndim(5),master,COMM)
               ENDIF
               IF(ndim(6).GT.0)THEN
                  call ls_mpibcast(lsmpibufferCha(1:ndim(6)),ndim(6),master,COMM)
               ENDIF

               call nullify_mpibuffer

            else  ! Sender: Send buffer to receiver
               call ls_mpisendrecv(ndim,6,comm,sender,receiver)

               IF(ndim(1).GT.0)THEN
                  call ls_mpisendrecv(lsmpibufferDP(1:ndim(1)),ndim(1),comm,sender,receiver)
               ENDIF
               IF(ndim(2).GT.0)THEN
                  call ls_mpisendrecv(lsmpibufferInt4(1:ndim(2)),ndim(2),comm,sender,receiver)
               ENDIF
               IF(ndim(3).GT.0)THEN
                  call ls_mpisendrecv(lsmpibufferInt8(1:ndim(3)),ndim(3),comm,sender,receiver)
               ENDIF
               IF(ndim(4).GT.0)THEN
                  call ls_mpisendrecv(lsmpibufferSho(1:ndim(4)),ndim(4),comm,sender,receiver)
               ENDIF
               IF(ndim(5).GT.0)THEN
                  call ls_mpisendrecv(lsmpibufferLog(1:ndim(5)),ndim(5),comm,sender,receiver)
               ENDIF
               IF(ndim(6).GT.0)THEN
                  call ls_mpisendrecv(lsmpibufferCha(1:ndim(6)),ndim(6),comm,sender,receiver)
               ENDIF
            end if

            call nullify_mpibuffer

         ELSE ! Receiver - dealloc buffer 
            modula = mod(iSho,int_to_short)
            additional = 0
            IF (modula.GT.0) additional = int_to_short - modula
            IF ((iDP.NE.nDP).OR.(iInt4.NE.nInteger4).OR.(iInt8.NE.nInteger8)&
                &.OR.(iSho.NE.(nShort-additional)).OR.&
     &          (iLog.NE.nLog).OR.(iCha.NE.nCha)) THEN
              IF(iDP.NE.nDP)                  write(*,*) 'The full Buffer has not been used DP',iDP,nDP
              IF(iInt4.NE.nInteger4)          write(*,*) 'The full Buffer has not been used Int4',iInt4,nInteger4
              IF(iInt8.NE.nInteger8)          write(*,*) 'The full Buffer has not been used Int8',iInt8,nInteger8
              IF(iSho.NE.(nShort-additional)) write(*,*) 'The full Buffer has not been used Short',iSho,nShort-additional
              IF(iLog.NE.nLog)                write(*,*) 'The full Buffer has not been used Log',iLog,nLog
              IF(iCha.NE.nCha)                write(*,*) 'The full Buffer has not been used Cha',iCha,nCha
            ENDIF

            IF(iDP.NE.nDP)call lsquit('The full Buffer has not been used DP',-1)
            IF(iInt4.NE.nInteger4) call lsquit('The full Buffer has not been used Int4',-1)
            IF(iInt8.NE.nInteger8) call lsquit('The full Buffer has not been used Int8',-1)
            IF(iSho.NE.(nShort-additional)) call lsquit('The full Buffer has not been used Short',-1)
            IF(iLog.NE.nLog) call lsquit('The full Buffer has not been used Log',-1)
            IF(iCha.NE.nCha) call lsquit('The full Buffer has not been used Cha',-1)
            call nullify_mpibuffer
         ENDIF
      else if(Job.EQ.LSMPIREDUCTION)THEN
         IF(iDP.GT.0)THEN
            call lsmpi_barrier(comm)
            call lsmpi_reduction(lsmpibufferDP(1:iDP),iDP,master,comm)
         ENDIF
         IF(iInt4.GT.0)THEN
            call lsmpi_barrier(comm)
            call lsmpi_int_reduction(lsmpibufferInt4(1:iInt4),iInt4,master,comm)
         ENDIF
         IF(iInt8.GT.0)THEN
            call lsmpi_barrier(comm)
            call lsmpi_int_reduction(lsmpibufferInt8(1:iInt8),iInt8,master,comm)
         ENDIF
         IF(iSho.GT.0)THEN
            call lsmpi_barrier(comm)
            call lsmpi_sho_reduction(lsmpibufferSho(1:iSho),iSho,master,comm)
         ENDIF
         IF(iLog.GT.0)THEN
            call lsquit('implement a logical reduction',-1)
            !call ls_log_reduction(lsmpibufferLog,iLog,master,COMM)
         ENDIF
         IF(iCha.GT.0)THEN
            call lsquit('implement a character reduction',-1)
            !call ls_mpibcast(lsmpibufferCha,iCha,master,COMM)
         ENDIF
         IF(mynum.NE.master)THEN
            IF(MemModParamPrintMemory)THEN
               print*,'MemModParamPrintMemory   mynum',mynum
               CALL Print_Memory_info(MemModParamPrintMemorylupri,'FinalBuffer: Before deallocation')
               Write(MemModParamPrintMemorylupri,*)'# DoublePrecison elements',size(lsmpibufferDP,kind=long)
               Write(MemModParamPrintMemorylupri,*)'# Integer 4 elements     ',size(lsmpibufferInt4,kind=long)
               Write(MemModParamPrintMemorylupri,*)'# Integer 8 elements     ',size(lsmpibufferInt8,kind=long)
               Write(MemModParamPrintMemorylupri,*)'# Short integer elements ',size(lsmpibufferSho,kind=long)
               Write(MemModParamPrintMemorylupri,*)'# Logical elements       ',size(lsmpibufferLog,kind=long)
               Write(MemModParamPrintMemorylupri,*)'# Character elements     ',size(lsmpibufferCha,kind=long)
            ENDIF
            call mem_dealloc(lsmpibufferDP)
            call mem_dealloc(lsmpibufferInt4)
            call mem_dealloc(lsmpibufferInt8)
            call mem_dealloc(lsmpibufferSho)
            call mem_dealloc(lsmpibufferLog)
            call mem_dealloc(lsmpibufferCha)
            IF(MemModParamPrintMemory)THEN
               print*,'MemModParamPrintMemory   mynum',mynum
               CALL Print_Memory_info(MemModParamPrintMemorylupri,'FinalBuffer: After deallocation')
            ENDIF
         ELSE
            AddToBuffer = .FALSE.
         ENDIF
         iDP = 0
         iInt4 = 0
         iInt8 = 0
         iSho = 0
         iLog = 0
         iCha = 0
      else if(Job.EQ.LSMPIREDUCTIONmaster)THEN
         IF(mynum.NE.master)THEN
            call lsquit('Programming error in ls_mpiFinalizeBuffer',-1)
         ENDIF
         iDP = 0
         iInt4 = 0
         iInt8 = 0
         iSho = 0
         iLog = 0
         iCha = 0
         IF(MemModParamPrintMemory)THEN
            print*,'MemModParamPrintMemory   mynum',mynum
            CALL Print_Memory_info(MemModParamPrintMemorylupri,'FinalBuffer: Before deallocation LSMPIREDUCTIONmaster')
            Write(MemModParamPrintMemorylupri,*)'# DoublePrecison elements',size(lsmpibufferDP,kind=long)
            Write(MemModParamPrintMemorylupri,*)'# Integer 4 elements     ',size(lsmpibufferInt4,kind=long)
            Write(MemModParamPrintMemorylupri,*)'# Integer 8 elements     ',size(lsmpibufferInt8,kind=long)
            Write(MemModParamPrintMemorylupri,*)'# Short integer elements ',size(lsmpibufferSho,kind=long)
            Write(MemModParamPrintMemorylupri,*)'# Logical elements       ',size(lsmpibufferLog,kind=long)
            Write(MemModParamPrintMemorylupri,*)'# Character elements     ',size(lsmpibufferCha,kind=long)
         ENDIF
         call mem_dealloc(lsmpibufferDP)
         call mem_dealloc(lsmpibufferInt4)
         call mem_dealloc(lsmpibufferInt8)
         call mem_dealloc(lsmpibufferSho)
         call mem_dealloc(lsmpibufferLog)
         call mem_dealloc(lsmpibufferCha)
         IF(MemModParamPrintMemory)THEN
            print*,'MemModParamPrintMemory   mynum',mynum
            CALL Print_Memory_info(MemModParamPrintMemorylupri,'FinalBuffer: After deallocation LSMPIREDUCTIONmaster')
         ENDIF
      ELSE
         call lsquit('Job specifier invalid in ls_mpiFinalizeBuffer',-1)
      ENDIF
#endif
    end subroutine ls_mpiFinalizeBuffer

#ifdef VAR_MPI

    !> \brief Dealloc (if associated) and nullify MPI buffers
    !> \author Kasper Kristensen
    !> \date March 2012
    subroutine nullify_mpibuffer

      use memory_handling

      implicit none

      if(associated(lsmpibufferDP)) then
         call mem_dealloc(lsmpibufferDP)
      end if

      if(associated(lsmpibufferInt4)) then
         call mem_dealloc(lsmpibufferInt4)
      end if

      if(associated(lsmpibufferInt8)) then
         call mem_dealloc(lsmpibufferInt8)
      end if

      if(associated(lsmpibufferSho)) then
         call mem_dealloc(lsmpibufferSho)
      end if

      if(associated(lsmpibufferLog)) then
         call mem_dealloc(lsmpibufferLog)
      end if

      if(associated(lsmpibufferCha)) then
         call mem_dealloc(lsmpibufferCha)
      end if

      nullify(lsmpibufferDP)
      nullify(lsmpibufferInt4)
      nullify(lsmpibufferInt8)
      nullify(lsmpibufferSho)
      nullify(lsmpibufferLog)
      nullify(lsmpibufferCha)


    end subroutine nullify_mpibuffer

#endif

#ifdef VAR_MPI
    subroutine increaselsmpibufferInt4(add)
      use memory_handling
      implicit none
      integer(kind=8) :: add
      integer(kind=4),pointer :: tmpbuffer(:)
!
      integer(kind=long) :: i,n
      n = size(lsmpibufferInt4,kind=long)
      call mem_alloc(tmpbuffer,n)
      DO i=1,n
         tmpbuffer(i) = lsmpibufferInt4(i)
      ENDDO
      call mem_dealloc(lsmpibufferInt4)
      nInteger4 = nInteger4 + MIN(MAX(incremInteger,add,nInteger4),MaxIncreaseSize)
      call mem_alloc(lsmpibufferInt4,nInteger4)
      DO i=1,n
         IF(i.GT.size(lsmpibufferInt4,kind=long))call lsquit('errorTK',-1)
         lsmpibufferInt4(i) = tmpbuffer(i)
      ENDDO
      call mem_dealloc(tmpbuffer)

    end subroutine increaselsmpibufferInt4

    subroutine increaselsmpibufferInt8(add)
      use memory_handling
      implicit none
      integer(kind=8) :: add
      integer(kind=8),pointer :: tmpbuffer(:)
      integer(kind=long) :: i,n
      n = size(lsmpibufferInt8,kind=long)
      call mem_alloc(tmpbuffer,n)
      DO i=1,n
         tmpbuffer(i) = lsmpibufferInt8(i)
      ENDDO
      call mem_dealloc(lsmpibufferInt8)
      nInteger8 = nInteger8 + MIN(MAX(incremInteger,add,nInteger8),MaxIncreaseSize)
      call mem_alloc(lsmpibufferInt8,nInteger8)
      DO i=1,n
         lsmpibufferInt8(i) = tmpbuffer(i)
      ENDDO
      call mem_dealloc(tmpbuffer)

    end subroutine increaselsmpibufferInt8

    subroutine increaselsmpibufferSho(add)
      use memory_handling
      implicit none
      integer :: add
      integer(kind=short),pointer :: tmpbuffer(:)
!
      integer(kind=long) :: i,n
      n = size(lsmpibufferSho,kind=long)
      call mem_alloc(tmpbuffer,n)
      DO i=1,n
         tmpbuffer(i) = lsmpibufferSho(i)
      ENDDO
      call mem_dealloc(lsmpibufferSho)
      nShort = nShort + MIN(MAX(incremShort,add*int_to_short,nShort),MaxIncreaseSize)
      call mem_alloc(lsmpibufferSho,nShort)
      DO i=1,n
         lsmpibufferSho(i) = tmpbuffer(i)
      ENDDO
      call mem_dealloc(tmpbuffer)

    end subroutine increaselsmpibufferSho

    subroutine increaselsmpibufferLog(add)
      use memory_handling
      implicit none
      integer :: add
      logical,pointer :: tmpbuffer(:)
!
      integer(kind=long) :: i,n
      n = size(lsmpibufferLog,kind=long)
      call mem_alloc(tmpbuffer,n)
      DO i=1,n
         tmpbuffer(i) = lsmpibufferLog(i)
      ENDDO
      call mem_dealloc(lsmpibufferLog)
      nLog = nLog + MIN(MAX(incremLog,add,nLog),MaxIncreaseSize)
      call mem_alloc(lsmpibufferLog,nLog)
      DO i=1,n
         lsmpibufferLog(i) = tmpbuffer(i)
      ENDDO
      call mem_dealloc(tmpbuffer)

    end subroutine increaselsmpibufferLog

    subroutine increaselsmpibufferCha(add)
      use memory_handling
      implicit none
      integer :: add
      character,pointer :: tmpbuffer(:)
!
      integer(kind=long) :: i,n
      n = size(lsmpibufferCha,kind=long)
      call mem_alloc(tmpbuffer,n)
      DO i=1,n
         tmpbuffer(i) = lsmpibufferCha(i)
      ENDDO
      call mem_dealloc(lsmpibufferCha)
      nCha = nCha + MIN(MAX(incremCha,add,nCha),MaxIncreaseSize)
      call mem_alloc(lsmpibufferCha,nCha)
      DO i=1,n
         lsmpibufferCha(i) = tmpbuffer(i)
      ENDDO
      call mem_dealloc(tmpbuffer)

    end subroutine increaselsmpibufferCha

    subroutine increaselsmpibufferDP(add)
      use memory_handling
      implicit none
      integer :: add
      real(realk),pointer :: tmpbuffer(:)
!
      integer(kind=long) :: i,n
      n = size(lsmpibufferDP,kind=long)
      call mem_alloc(tmpbuffer,n)
      DO i=1,n
         tmpbuffer(i) = lsmpibufferDP(i)
      ENDDO
      call mem_dealloc(lsmpibufferDP)
      nDP = n + MIN(MAX(incremDP,add,n),MaxIncreaseSize)
      call mem_alloc(lsmpibufferDP,nDP)
      DO i=1,n
         lsmpibufferDP(i) = tmpbuffer(i)
      ENDDO
      call mem_dealloc(tmpbuffer)

    end subroutine increaselsmpibufferDP
#endif



    subroutine LSMPI_MYFAIL(IERR)
      implicit none
      integer(kind=ls_mpik) :: IERR
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: IERR2,IERRCL
      CHARACTER(len=40) :: ERRBUF
      CALL MPI_ERROR_CLASS(IERR,IERRCL,IERR2)
      IF (IERRCL.EQ.MPI_SUCCESS) THEN
         ERRBUF = 'No error'
      ELSE IF (IERRCL.EQ.MPI_ERR_BUFFER) THEN
         ERRBUF = 'Invalid buffer pointer'
      ELSE IF (IERRCL.EQ.MPI_ERR_COUNT) THEN
         ERRBUF = 'Invalid count argument'
      ELSE IF (IERRCL.EQ.MPI_ERR_TYPE) THEN
         ERRBUF = 'Invalid datatype argument'
      ELSE IF (IERRCL.EQ.MPI_ERR_TAG) THEN
         ERRBUF = 'Invalid tag argument'
      ELSE IF (IERRCL.EQ.MPI_ERR_COMM) THEN
         ERRBUF = 'Invalid communicator'
      ELSE IF (IERRCL.EQ.MPI_ERR_RANK) THEN
         ERRBUF = 'Invalid rank'
      ELSE IF (IERRCL.EQ.MPI_ERR_REQUEST) THEN
         ERRBUF = 'Invalid request (handle)'
      ELSE IF (IERRCL.EQ.MPI_ERR_ROOT) THEN
         ERRBUF = 'Invalid root'
      ELSE IF (IERRCL.EQ.MPI_ERR_GROUP) THEN
         ERRBUF = 'Invalid group'
      ELSE IF (IERRCL.EQ.MPI_ERR_OP) THEN
         ERRBUF = 'Invalid operation'
      ELSE IF (IERRCL.EQ.MPI_ERR_TOPOLOGY) THEN
         ERRBUF = 'Invalid topology'
      ELSE IF (IERRCL.EQ.MPI_ERR_DIMS) THEN
         ERRBUF = 'Invalid dimension argument'
      ELSE IF (IERRCL.EQ.MPI_ERR_ARG) THEN
         ERRBUF = 'Invalid argument of some other kind'
      ELSE IF (IERRCL.EQ.MPI_ERR_UNKNOWN) THEN
         ERRBUF = 'Unknown error'
      ELSE IF (IERRCL.EQ.MPI_ERR_TRUNCATE) THEN
         ERRBUF = 'Message truncated on receive'
      ELSE IF (IERRCL.EQ.MPI_ERR_OTHER) THEN
         ERRBUF = 'Known error not in this list'
      ELSE IF (IERRCL.EQ.MPI_ERR_INTERN) THEN
         ERRBUF = 'Internal MPI (implementation) error'
      ELSE IF (IERRCL.EQ.MPI_ERR_IN_STATUS) THEN
         ERRBUF = 'Error code is in status'
      ELSE IF (IERRCL.EQ.MPI_ERR_PENDING) THEN
         ERRBUF = 'Pending request'
      ELSE IF (IERRCL.EQ.MPI_ERR_LASTCODE) THEN
         ERRBUF = 'Last error code'
      ELSE
!        something we didn't know ...
         WRITE(6,'(A,I4)') 'MPI error class',IERRCL
      END IF
!
      WRITE(6,'(/A)') ' ERROR in MPI : '//ERRBUF
!
      CALL LSQUIT('Error detected in MPI. Please consult dalton output!',-1)
!
#endif
    end subroutine LSMPI_MYFAIL

    subroutine lsmpi_reduction_integer4_wrapper8(buffer,n,master,comm)
    implicit none
    integer(kind=ls_mpik),intent(in) :: comm
    integer(kind=ls_mpik),intent(in) :: master
    integer(kind=8) :: n
    integer(kind=4) :: buffer(:)
#ifdef VAR_MPI
    integer(kind=4) :: n4
    integer(kind=8) :: i,k
    integer(kind=ls_mpik) :: ierr,mynum,thesize
    real(realk) :: null 
      IERR=0
    !loop over batches, which contain a number of elements,
    !describable by 32 bit integers, here 2E9
    if(ls_mpik==4)then
      k=SPLIT_MPI_MSG
      do i=1,n,k
        n4=k
        !if((n-i)<k)n4=mod(n,k)
        if(((n-i)<k).and.(mod(n-i+1,k)/=0))n4=mod(n,k)
        call lsmpi_reduction_integer4(buffer(i:i+n4-1),n4,master,comm)
      enddo
    else
      call get_rank_for_comm(comm,mynum)
      THESIZE = n
      !if(THESIZE.NE.n)THEN
      !   call lsquit('lsmpi error in lmpi_reduction_integer_wrapper8',-1)
      !endif
      IF(mynum.EQ.master) THEN
        CALL MPI_REDUCE(MPI_IN_PLACE,BUFFER,thesize,MPI_INTEGER4,MPI_SUM,&
             &master,comm,IERR)
        IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
      ELSE
        CALL MPI_REDUCE(BUFFER,NULL,thesize,MPI_INTEGER4,MPI_SUM,&
             &master,comm,IERR)
        IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
      ENDIF
    endif
#endif
    end subroutine lsmpi_reduction_integer4_wrapper8
    
    subroutine lsmpi_reduction_integer4(buffer,n1,master,comm)
    implicit none
    integer(kind=ls_mpik),intent(in) :: comm   ! communicator
    integer(kind=ls_mpik) :: master
    integer(kind=4) :: n1
    integer(kind=4) :: buffer(:)
#ifdef VAR_MPI
    integer(kind=ls_mpik) :: ierr,mynum,thesize
    real(realk) :: null 
      IERR=0
    call get_rank_for_comm(comm,mynum)
    THESIZE = n1
    IF(mynum.EQ.master) THEN
      CALL MPI_REDUCE(MPI_IN_PLACE,BUFFER,thesize,MPI_INTEGER4,MPI_SUM,&
           &master,comm,IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
    ELSE
      CALL MPI_REDUCE(BUFFER,NULL,thesize,MPI_INTEGER4,MPI_SUM,&
           &master,comm,IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
    ENDIF
#endif
  end subroutine lsmpi_reduction_integer4


  subroutine lsmpi_reduction_integer8_wrapper8(buffer,n,master,comm)
    implicit none
    integer(kind=ls_mpik),intent(in) :: comm
    integer(kind=ls_mpik),intent(in) :: master
    integer(kind=8) :: n
    integer(kind=8) :: buffer(:)
#ifdef VAR_MPI
    integer(kind=4) :: n4
    integer(kind=8) :: i,k
    integer(kind=ls_mpik) :: mynum
    integer(kind=ls_mpik) :: ierr
    real(realk) :: null 
      IERR=0
    !loop over batches, which contain a number of elements,
    !describable by 32 bit integers, here 2E9
    if(ls_mpik==4)then
      k=SPLIT_MPI_MSG
      do i=1,n,k
        n4=k
        !if((n-i)<k)n4=mod(n,k)
        if(((n-i)<k).and.(mod(n-i+1,k)/=0))n4=mod(n,k)
        call lsmpi_reduction_integer8(buffer(i:i+n4-1),n4,master,comm)
      enddo
    else
      call get_rank_for_comm(comm,mynum)
      IF(mynum.EQ.master) THEN
        CALL MPI_REDUCE(MPI_IN_PLACE,BUFFER,n,MPI_INTEGER8,MPI_SUM,&
             &master,comm,IERR)
        IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
      ELSE
        CALL MPI_REDUCE(BUFFER,NULL,n,MPI_INTEGER8,MPI_SUM,&
             &master,comm,IERR)
        IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
      ENDIF
    endif
#endif
    end subroutine lsmpi_reduction_integer8_wrapper8
    subroutine lsmpi_reduction_integer8(buffer,n1,master,comm)
    implicit none
    integer(kind=ls_mpik),intent(in) :: comm   ! communicator
    integer(kind=ls_mpik),intent(in) :: master
    integer(kind=4) :: n1
    integer(kind=8) :: buffer(:)
#ifdef VAR_MPI
    integer(kind=ls_mpik) :: mynum
    integer(kind=ls_mpik) :: ierr
    real(realk) :: null 
      IERR=0
    call get_rank_for_comm(comm,mynum)
    IF(mynum.EQ.master) THEN
      CALL MPI_REDUCE(MPI_IN_PLACE,BUFFER,n1,MPI_INTEGER8,MPI_SUM,&
           &master,comm,IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
    ELSE
      CALL MPI_REDUCE(BUFFER,NULL,n1,MPI_INTEGER8,MPI_SUM,&
           &master,comm,IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
    ENDIF
#endif
  end subroutine lsmpi_reduction_integer8

  subroutine lsmpi_reduction_integer4M4(buffer,n1,n2,master,comm)
    implicit none
    integer(kind=ls_mpik),intent(in) :: comm   ! communicator
    integer(kind=ls_mpik),intent(in) :: master
    integer(kind=4) :: n1,n2
    integer(kind=4) :: buffer(n1,n2)
#ifdef VAR_MPI
    integer(kind=8) :: n,i,k
    integer(kind=4) :: n4
    integer(kind=4), pointer :: buf(:)
    n=n1*n2
    call ass_44I2to1(buffer,buf,[n1,n2])
    if(ls_mpik==4)then
      k=SPLIT_MPI_MSG
      do i=1,n,k
        n4=k
        !if((n-i)<k)n4=mod(n,k)
        if(((n-i)<k).and.(mod(n-i+1,k)/=0))n4=mod(n,k)
        call lsmpi_reduction_integer4(buf(i:i+n4-1),n4,master,comm)
      enddo
    else
      call lsmpi_reduction_integer4_wrapper8(buf,n,master,comm)
    endif
    nullify(buf)
#endif
  end subroutine lsmpi_reduction_integer4M4
  subroutine lsmpi_reduction_integer4M8(buffer,n1,n2,master,comm)
    implicit none
    integer(kind=ls_mpik),intent(in) :: comm   ! communicator
    integer(kind=ls_mpik),intent(in) :: master
    integer(kind=8) :: n1,n2
    integer(kind=4) :: buffer(n1,n2)
#ifdef VAR_MPI
    integer(kind=8) :: n,i,k
    integer(kind=4) :: n4
    integer(kind=4), pointer :: buf(:)
    n=n1*n2
    call ass_48I2to1(buffer,buf,[n1,n2])
    if(ls_mpik==4)then
      k=SPLIT_MPI_MSG
      do i=1,n,k
        n4=k
        !if((n-i)<k)n4=mod(n,k)
        if(((n-i)<k).and.(mod(n-i+1,k)/=0))n4=mod(n,k)
        call lsmpi_reduction_integer4(buf(i:i+n4-1),n4,master,comm)
      enddo
    else
      call lsmpi_reduction_integer4_wrapper8(buf,n,master,comm)
    endif
    nullify(buf)
#endif
  end subroutine lsmpi_reduction_integer4M8
  subroutine lsmpi_reduction_integer8M4(buffer,n1,n2,master,comm)
    implicit none
    integer(kind=ls_mpik),intent(in) :: comm   ! communicator
    integer(kind=ls_mpik),intent(in) :: master
    integer(kind=4) :: n1,n2
    integer(kind=8) :: buffer(n1,n2)
#ifdef VAR_MPI
    integer(kind=8) :: n,i,k
    integer(kind=4) :: n4
    integer(kind=8), pointer :: buf(:)
    n=n1*n2
    call ass_84I2to1(buffer,buf,[n1,n2])
    if(ls_mpik==4)then
      k=SPLIT_MPI_MSG
      do i=1,n,k
        n4=k
        !if((n-i)<k)n4=mod(n,k)
        if(((n-i)<k).and.(mod(n-i+1,k)/=0))n4=mod(n,k)
        call lsmpi_reduction_integer8(buf(i:i+n4-1),n4,master,comm)
      enddo
    else
      call lsmpi_reduction_integer8_wrapper8(buf,n,master,comm)
    endif
    nullify(buf)
#endif
  end subroutine lsmpi_reduction_integer8M4
  subroutine lsmpi_reduction_integer8M8(buffer,n1,n2,master,comm)
    implicit none
    integer(kind=ls_mpik),intent(in) :: comm   ! communicator
    integer(kind=ls_mpik),intent(in) :: master
    integer(kind=8) :: n1,n2
    integer(kind=8) :: buffer(n1,n2)
#ifdef VAR_MPI
    integer(kind=8) :: n,i,k
    integer(kind=4) :: n4
    integer(kind=8), pointer :: buf(:)
    call ass_88I2to1(buffer,buf,[n1,n2])
    n=n1*n2
    if(ls_mpik==4)then
      k=SPLIT_MPI_MSG
      do i=1,n,k
        n4=k
        !if((n-i)<k)n4=mod(n,k)
        if(((n-i)<k).and.(mod(n-i+1,k)/=0))n4=mod(n,k)
        call lsmpi_reduction_integer8(buf(i:i+n4-1),n4,master,comm)
      enddo
    else
      call lsmpi_reduction_integer8_wrapper8(buf,n,master,comm)
    endif
    nullify(buf)
#endif
  end subroutine lsmpi_reduction_integer8M8


  subroutine lsmpi_reduction_realk_single(buffer,master,comm)
    implicit none
    integer(kind=ls_mpik),intent(in) :: comm   ! communicator
    integer(kind=ls_mpik),intent(in) :: master
    real(realk) :: buffer
#ifdef VAR_MPI
    integer(kind=ls_mpik) :: ierr,datatype,n
    real(realk) :: null
    integer(kind=ls_mpik) :: mynum,count
      IERR=0
    call get_rank_for_comm(comm,mynum)    
    datatype=MPI_DOUBLE_PRECISION  
    n=1
    IF(mynum.EQ.master) THEN
      CALL MPI_REDUCE(MPI_IN_PLACE,BUFFER,n,datatype,MPI_SUM,&
           &master,comm,IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
    ELSE
      CALL MPI_REDUCE(BUFFER,NULL,n,datatype,MPI_SUM,&
           &master,comm,IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
    ENDIF
#endif
  end subroutine lsmpi_reduction_realk_single


  subroutine lsmpi_reduction_realk_wrapper8(buffer,n,master,comm)
    implicit none
    integer(kind=ls_mpik),intent(in) :: comm   ! communicator
    integer(kind=ls_mpik),intent(in) :: master
    integer(kind=8),intent(in) :: n
    real(realk)                :: buffer(n)
#ifdef VAR_MPI
    integer(kind=4)    :: n4
    integer(kind=8)    :: i,k
    integer(kind=ls_mpik) :: ierr,mynum,nelms
    real(realk) :: null
    !loop over batches, which contain a number of elements,
    !describable by 32 bit integers, here 2E9
    if(ls_mpik==4)then
      k=SPLIT_MPI_MSG
      do i=1,n,k
        n4=k
        if(((n-i)<k).and.(mod(n-i+1,k)/=0))n4=mod(n,k)
        call lsmpi_reduction_realk(buffer(i:i+n4-1),n4,master,comm)
      enddo
    else
      IERR=0
      call get_rank_for_comm(comm,mynum)
      nelms=n
      IF(mynum.EQ.master) THEN
        CALL MPI_REDUCE(MPI_IN_PLACE,BUFFER,nelms,MPI_DOUBLE_PRECISION,MPI_SUM,&
             &master,comm,IERR)
        IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
      ELSE
        CALL MPI_REDUCE(BUFFER,NULL,nelms,MPI_DOUBLE_PRECISION,MPI_SUM,&
             &master,comm,IERR)
        IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
      ENDIF
    endif
#endif
  end subroutine lsmpi_reduction_realk_wrapper8

  subroutine lsmpi_reduction_realk(buffer,n1,master,comm)
    implicit none
    integer(kind=ls_mpik),intent(in) :: comm   ! communicator
    integer(kind=ls_mpik),intent(in) :: master
    integer(kind=4) :: n1
    real(realk) :: buffer(:)
#ifdef VAR_MPI
    integer(kind=ls_mpik) :: ierr,mynum,nelms
    real(realk) :: null
    IERR=0
    call get_rank_for_comm(comm,mynum)
    nelms=n1
    IF(mynum.EQ.master) THEN
      CALL MPI_REDUCE(MPI_IN_PLACE,BUFFER,nelms,MPI_DOUBLE_PRECISION,MPI_SUM,&
           &master,comm,IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
    ELSE
      CALL MPI_REDUCE(BUFFER,NULL,nelms,MPI_DOUBLE_PRECISION,MPI_SUM,&
           &master,comm,IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
    ENDIF
#endif
  end subroutine lsmpi_reduction_realk

  subroutine lsmpi_reduction_realkM4(buffer,n1,n2,master,comm)
    implicit none
    integer(kind=ls_mpik),intent(in) :: comm   ! communicator
    integer(kind=ls_mpik),intent(in) :: master
    integer(kind=4) ::n1,n2
    real(realk) :: buffer(n1,n2)
#ifdef VAR_MPI
    integer(kind=8) :: n
    real(realk),pointer :: buf(:)
    n=n1*n2
    call ass_4D2to1(buffer,buf,[n1,n2])
    call lsmpi_reduction_realk_wrapper8(buf,n,master,comm)
    nullify(buf)
#endif
  end subroutine lsmpi_reduction_realkM4
  subroutine lsmpi_reduction_realkM8(buffer,n1,n2,master,comm)
    implicit none
    integer(kind=ls_mpik),intent(in) :: comm   ! communicator
    integer(kind=ls_mpik),intent(in) :: master
    integer(kind=8) ::n1,n2
    real(realk) :: buffer(n1,n2)
#ifdef VAR_MPI
    integer(kind=8) :: n
    real(realk),pointer :: buf(:)
    n=n1*n2
    call ass_8D2to1(buffer,buf,[n1,n2])
    call lsmpi_reduction_realk_wrapper8(buf,n,master,comm)
    nullify(buf)
#endif
  end subroutine lsmpi_reduction_realkM8

  subroutine lsmpi_reduction_realkT4(buffer,n1,n2,n3,master,comm)
    implicit none
    integer(kind=ls_mpik),intent(in) :: comm   ! communicator
    integer(kind=ls_mpik),intent(in) :: master
    integer(kind=4) ::n1,n2,n3
    real(realk) :: buffer(n1,n2,n3)
#ifdef VAR_MPI
    integer(kind=8) :: n
    real(realk),pointer :: buf(:)
    n=n1*n2*n3
    call ass_4D3to1(buffer,buf,[n1,n2,n3])
    call lsmpi_reduction_realk_wrapper8(buf,n,master,comm)
    nullify(buf)
#endif
  end subroutine lsmpi_reduction_realkT4
  subroutine lsmpi_reduction_realkT8(buffer,n1,n2,n3,master,comm)
    implicit none
    integer(kind=ls_mpik),intent(in) :: comm   ! communicator
    integer(kind=ls_mpik),intent(in) :: master
    integer(kind=8) ::n1,n2,n3
    real(realk) :: buffer(n1,n2,n3)
#ifdef VAR_MPI
    integer(kind=8) :: n
    real(realk),pointer :: buf(:)
    n=n1*n2*n3
    call ass_8D3to1(buffer,buf,[n1,n2,n3])
    call lsmpi_reduction_realk_wrapper8(buf,n,master,comm)
    nullify(buf)
#endif
  end subroutine lsmpi_reduction_realkT8

  subroutine lsmpi_reduce_realk_min(buffer,dest,comm)
    implicit none
    real(realk), intent(inout):: buffer
    integer(kind=ls_mpik),intent(in) :: comm,dest
    integer(kind=ls_mpik) :: one_el, IERR
    real(realk) :: sendbuffer
    IERR = 0
    one_el = int(1,kind=ls_mpik)
#ifdef VAR_MPI

    sendbuffer = buffer
    CALL MPI_REDUCE(sendbuffer,BUFFER,one_el,MPI_DOUBLE_PRECISION,MPI_MIN,&
           & dest,comm, IERR)
#endif
  end subroutine lsmpi_reduce_realk_min
  subroutine lsmpi_reduce_realk_max(buffer,dest,comm)
    implicit none
    real(realk), intent(inout):: buffer
    integer(kind=ls_mpik),intent(in) :: comm,dest
    integer(kind=ls_mpik) :: one_el, IERR
    real(realk) :: sendbuffer
    IERR = 0
    one_el = int(1,kind=ls_mpik)
#ifdef VAR_MPI

    sendbuffer = buffer
    CALL MPI_REDUCE(sendbuffer,BUFFER,one_el,MPI_DOUBLE_PRECISION,MPI_MAX,&
           & dest,comm, IERR)
#endif
  end subroutine lsmpi_reduce_realk_max

  ! MPI all reduce within local group (infpar%lg_comm communicator)
  subroutine lsmpi_allreduce_int8V_wrapper8(rbuffer,n1,comm)
    implicit none
    integer(kind=8)                  :: n1
    integer(kind=8)                  :: rbuffer(n1)
    integer(kind=ls_mpik),intent(in) :: comm
#ifdef VAR_MPI
    integer(kind=4) :: n4,k,i
    integer(kind=ls_mpik) :: ierr,nelms
    if(ls_mpik==4)then
      k=SPLIT_MPI_MSG
      do i=1,n1,k
        n4=k
        if(((n1-i)<k).and.(mod(n1-i+1,k)/=0))n4=mod(n1,k)
        call lsmpi_allreduce_int8V(rbuffer(i:i+n4-1),n4,comm)
      enddo
    else
      IERR=0
      nelms=n1
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,RBUFFER,nelms,MPI_INTEGER8,MPI_SUM,&
           & comm, IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
    endif
#endif
  end subroutine lsmpi_allreduce_int8V_wrapper8
  subroutine lsmpi_allreduce_int8V(rbuffer,n1,comm)
    implicit none
    integer(kind=4)                  :: n1
    integer(kind=8)                  :: rbuffer(n1)
    integer(kind=ls_mpik),intent(in) :: comm
#ifdef VAR_MPI
    integer(kind=ls_mpik) :: ierr,nelms
    IERR=0
    nelms=n1
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,RBUFFER,nelms,MPI_INTEGER8,MPI_SUM,&
         &comm,IERR)
    IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
  end subroutine lsmpi_allreduce_int8V
  subroutine lsmpi_allreduce_int4V_wrapper8(rbuffer,n1,comm)
    implicit none
    integer(kind=8)                  :: n1
    integer(kind=4)                  :: rbuffer(n1)
    integer(kind=ls_mpik),intent(in) :: comm
#ifdef VAR_MPI
    integer(kind=4) :: n4,k,i
    integer(kind=ls_mpik) :: ierr,nelms
    if(ls_mpik==4)then
      k=SPLIT_MPI_MSG
      do i=1,n1,k
        n4=k
        if(((n1-i)<k).and.(mod(n1-i+1,k)/=0))n4=mod(n1,k)
        call lsmpi_allreduce_int4V(rbuffer(i:i+n4-1),n4,comm)
      enddo
    else
      IERR=0
      nelms=n1
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,RBUFFER,nelms,MPI_INTEGER4,MPI_SUM,&
         & comm,IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
    endif
#endif
  end subroutine lsmpi_allreduce_int4V_wrapper8
  subroutine lsmpi_allreduce_int4V(rbuffer,n1,comm)
    implicit none
    integer(kind=4)                  :: n1
    integer(kind=4)                  :: rbuffer(n1)
    integer(kind=ls_mpik),intent(in) :: comm
#ifdef VAR_MPI
    integer(kind=ls_mpik) :: ierr,nelms
    IERR=0
    nelms=n1
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,RBUFFER,nelms,MPI_INTEGER4,MPI_SUM,&
         &comm,IERR)
    IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
  end subroutine lsmpi_allreduce_int4V

  ! MPI all reduce within local group (infpar%lg_comm communicator)
  subroutine lsmpi_allreduce_D(rbuffer,comm)
    implicit none
    real(realk)                      :: sbuffer,rbuffer
    integer(kind=ls_mpik),intent(in) :: comm
#ifdef VAR_MPI
    integer(kind=ls_mpik) :: ierr,nel
    IERR=0
    nel=1
    sbuffer=rbuffer
    CALL MPI_ALLREDUCE(SBUFFER,RBUFFER,nel,MPI_DOUBLE_PRECISION&
    &,MPI_SUM,comm,IERR)
    IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
  end subroutine lsmpi_allreduce_D

  subroutine lsmpi_allreduce_D1N4(rbuffer,n1,comm)
    implicit none
    integer(kind=4)                  :: n1
    real(realk)                      :: rbuffer(n1)
    integer(kind=ls_mpik),intent(in) :: comm
#ifdef VAR_MPI
    integer(ls_mpik) :: ierr,nel
    IERR=0
    nel=n1
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,RBUFFER,nel,MPI_DOUBLE_PRECISION,&
         &MPI_SUM,comm,IERR)
    IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
#endif
  end subroutine lsmpi_allreduce_D1N4
  subroutine lsmpi_allreduce_D1N8_parts(rbuffer,nelms,comm,split)
    implicit none
    integer(kind=8),intent(in)       :: nelms
    real(realk)                      :: rbuffer(nelms)
    integer(kind=ls_mpik),intent(in) :: comm
    integer,intent(in) :: split
#ifdef VAR_MPI
    integer(kind=8) :: n,i
    integer(kind=8) :: n1,n2
    n2 = int(split,kind=8)
    do i=1,nelms,split
      n=split
      n1 = nelms-i+1
      if(((nelms-i)<split).and.(mod(n1,n2)/=0))n=mod(nelms,split)
      call lsmpi_allreduce_D1N8(rbuffer(i:i+n-1),n,comm)
    enddo
#endif
  end subroutine lsmpi_allreduce_D1N8_parts
  subroutine lsmpi_allreduce_D1N4_parts(rbuffer,nelms,comm,split)
    implicit none
    integer(kind=4),intent(in)       :: nelms
    real(realk)                      :: rbuffer(nelms)
    integer(kind=ls_mpik),intent(in) :: comm
    integer,intent(in)               :: split
#ifdef VAR_MPI
    integer(kind=8) :: n,i
    integer(kind=8) :: n1,n2
    n2 = int(split,kind=8)
    do i=1,nelms,split
      n=split
      n1 = nelms-i+1
      if(((nelms-i)<split).and.(mod(n1,n2)/=0))n=mod(nelms,split)
      call lsmpi_allreduce_D1N8(rbuffer(i:i+n-1),n,comm)
    enddo
#endif
  end subroutine lsmpi_allreduce_D1N4_parts
  subroutine lsmpi_allreduce_D1N8(rbuffer,n,comm)
    implicit none
    integer(kind=8)                  :: n
    real(realk)                      :: rbuffer(n)
    integer(kind=ls_mpik),intent(in) :: comm
#ifdef VAR_MPI
    integer(ls_mpik) :: ierr,nel
    integer(kind=8)  :: i,k
    integer(kind=4)  :: n4
    if(ls_mpik==4)then
      k=SPLIT_MPI_MSG
      do i=1,n,k
        n4=k
        !if((n-i)<k)n4=mod(n,k)
        if(((n-i)<k).and.(mod(n-i+1,k)/=0))n4=mod(n,k)
        call lsmpi_allreduce_D1N4(rbuffer(i:i+n4-1),n4,comm)
      enddo
    else
      IERR=0
      nel=n
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,RBUFFER,nel,MPI_DOUBLE_PRECISION,&
           &MPI_SUM,comm,IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
    endif
#endif
  end subroutine lsmpi_allreduce_D1N8

  subroutine lsmpi_allreduce_D2(rbuffer,n1,n2,comm)
    implicit none
    integer                          :: n1,n2
    real(realk), dimension(n1,n2)    :: rbuffer
    integer(kind=ls_mpik),intent(in) :: comm
    integer(kind=8) :: n
#ifdef VAR_MPI
    n=n1*n2
    call lsmpi_allreduce_D1N8(rbuffer,n,comm)
#endif
  end subroutine lsmpi_allreduce_D2

  subroutine lsmpi_allreduce_D3(rbuffer,n1,n2,n3,comm)
    implicit none
    integer                          :: n1,n2,n3
    real(realk), dimension(n1,n2,n3) :: rbuffer
    integer(kind=ls_mpik),intent(in) :: comm
    integer(kind=8) :: n
#ifdef VAR_MPI
    n=(i8*n1)*n2*n3
    call lsmpi_allreduce_D1N8(rbuffer,n,comm)
#endif
  end subroutine lsmpi_allreduce_D3

  subroutine lsmpi_allreduce_D4(rbuffer,n1,n2,n3,n4,comm)
    implicit none
    integer                             :: n1,n2,n3,n4
    real(realk), dimension(n1,n2,n3,n4) :: rbuffer
    integer(kind=ls_mpik),intent(in)    :: comm
    integer(kind=8) :: n
#ifdef VAR_MPI
    n=(i8*n1)*n2*n3*n4
    call lsmpi_allreduce_D1N8(rbuffer,n,comm)
#endif
  end subroutine lsmpi_allreduce_D4

  subroutine lsmpi_local_reduction_int(buffer,master)
    implicit none
    integer(kind=ls_mpik) :: master
    integer :: buffer
#ifdef VAR_MPI
    integer(kind=ls_mpik) :: ierr,n1
    real(realk) :: null
    IERR=0
    n1=1
    IF(infpar%lg_mynum.EQ.master) THEN
      CALL MPI_REDUCE(MPI_IN_PLACE,BUFFER,n1,MPI_INTEGER,MPI_SUM,&
           &master,infpar%lg_comm,IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
    ELSE
      CALL MPI_REDUCE(BUFFER,NULL,n1,MPI_INTEGER,MPI_SUM,&
           &master,infpar%lg_comm,IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
    ENDIF
#endif
  end subroutine lsmpi_local_reduction_int


  subroutine lsmpi_local_reduction_intV(buffer,n1,master)
    implicit none
    integer(kind=ls_mpik) :: master
    integer :: n1
    integer :: buffer(:)
#ifdef VAR_MPI
    integer(kind=ls_mpik) :: ierr,nel
    real(realk) :: null
    IERR=0
    nel=n1
    IF(infpar%lg_mynum.EQ.master) THEN
      CALL MPI_REDUCE(MPI_IN_PLACE,BUFFER,nel,MPI_INTEGER,MPI_SUM,&
           &master,infpar%lg_comm,IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
    ELSE
      CALL MPI_REDUCE(BUFFER,NULL,nel,MPI_INTEGER,MPI_SUM,&
           &master,infpar%lg_comm,IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
    ENDIF
#endif
  end subroutine lsmpi_local_reduction_intV

  ! MPI reduction within local group (infpar%lg_comm communicator)
  subroutine lsmpi_local_reduction_realk(buffer,master)
    implicit none
    integer(kind=ls_mpik) :: master
    real(realk) :: buffer
#ifdef VAR_MPI
    integer(kind=ls_mpik) :: ierr,nel
    real(realk) :: null
    IERR=0
    nel=1
    IF(infpar%lg_mynum.EQ.master) THEN
      CALL MPI_REDUCE(MPI_IN_PLACE,BUFFER,nel,MPI_DOUBLE_PRECISION,MPI_SUM,&
           &master,infpar%lg_comm,IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
    ELSE
      CALL MPI_REDUCE(BUFFER,NULL,nel,MPI_DOUBLE_PRECISION,MPI_SUM,&
           &master,infpar%lg_comm,IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
    ENDIF
#endif
  end subroutine lsmpi_local_reduction_realk

  subroutine lsmpi_local_reduction_realkVN8_parts(rbuffer,nelms,master,split)
    implicit none
    integer(kind=8),intent(in) :: nelms
    real(realk) :: rbuffer(nelms)
    integer(kind=ls_mpik) :: master
    integer,intent(in) :: split
#ifdef VAR_MPI
    integer(kind=8) :: n,i,n1,n2
    n2 = split
    do i=1,nelms,split
      n=split
      n1 = nelms-i+1
      if(((nelms-i)<split).and.&
        &(mod(n1,n2)/=0))n=mod(nelms,split)
      call lsmpi_local_reduction_realkVN8(rbuffer(i:i+n-1),n,master)
    enddo
#endif
  end subroutine lsmpi_local_reduction_realkVN8_parts

  subroutine lsmpi_local_reduction_realkVN4_parts(rbuffer,nelms,master,split)
    implicit none
    integer(kind=4),intent(in) :: nelms
    real(realk) :: rbuffer(nelms)
    integer(kind=ls_mpik) :: master
    integer,intent(in) :: split
#ifdef VAR_MPI
    integer(kind=8) :: n,i,n1,n2
    n2 = split
    do i=1,nelms,split
      n=split
      n1 = nelms-i+1
      if(((nelms-i)<split).and.&
        &(mod(n1,n2)/=0))n=mod(nelms,split)
      call lsmpi_local_reduction_realkVN8(rbuffer(i:i+n-1),n,master)
    enddo
#endif
  end subroutine lsmpi_local_reduction_realkVN4_parts
  ! MPI reduction within local group (infpar%lg_comm communicator)
  subroutine lsmpi_local_reduction_realkVN8(buffer,n,master)
    implicit none
    integer(kind=ls_mpik) :: master
    integer(kind=8) :: n
    real(realk) :: buffer(n)
#ifdef VAR_MPI
    integer(kind=ls_mpik) :: ierr,nel
    integer(kind=8)  :: i,k
    integer(kind=4)  :: n4
    real(realk) :: null
    if(ls_mpik==4)then
      k=SPLIT_MPI_MSG
      do i=1,n,k
        n4=k
        !if((n-i)<k)n4=mod(n,k)
        if(((n-i)<k).and.(mod(n-i+1,k)/=0))n4=mod(n,k)
        call lsmpi_local_reduction_realkVN4(buffer(i:i+n4-1),n4,master)
      enddo
    else
      IERR=0
      nel=n
      IF(infpar%lg_mynum.EQ.master) THEN
        CALL MPI_REDUCE(MPI_IN_PLACE,BUFFER,nel,MPI_DOUBLE_PRECISION,MPI_SUM,&
             &master,infpar%lg_comm,IERR)
        IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
      ELSE
        CALL MPI_REDUCE(BUFFER,NULL,nel,MPI_DOUBLE_PRECISION,MPI_SUM,&
             &master,infpar%lg_comm,IERR)
        IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
      ENDIF
    endif
#endif
  end subroutine lsmpi_local_reduction_realkVN8

  subroutine lsmpi_local_reduction_realkVN4(buffer,n1,master)
    implicit none
    integer(kind=ls_mpik) :: master
    integer(kind=4) :: n1
    real(realk) :: buffer(n1)
#ifdef VAR_MPI
    integer(kind=ls_mpik) :: ierr,nel
    real(realk) :: null
    IERR=0
    nel=n1
    IF(infpar%lg_mynum.EQ.master) THEN
      CALL MPI_REDUCE(MPI_IN_PLACE,BUFFER,nel,MPI_DOUBLE_PRECISION,MPI_SUM,&
           &master,infpar%lg_comm,IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
    ELSE
      CALL MPI_REDUCE(BUFFER,NULL,nel,MPI_DOUBLE_PRECISION,MPI_SUM,&
           &master,infpar%lg_comm,IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
    ENDIF
#endif
  end subroutine lsmpi_local_reduction_realkVN4

  ! MPI reduction within local group (infpar%lg_comm communicator)
  subroutine lsmpi_local_reduction_realkM(buffer,n1,n2,master)
    implicit none
    integer(kind=ls_mpik) :: master
    integer               :: n1,n2
    real(realk) :: buffer(n1,n2)
#ifdef VAR_MPI
    real(realk),pointer :: buf(:)
    integer(kind=8) :: n
    n=(i8*n1)*n2
    call ass_8D2to1(buffer,buf,[i8*n1,i8*n2])
    call lsmpi_local_reduction_realkVN8(buf,n,master)
    nullify(buf)
#endif
  end subroutine lsmpi_local_reduction_realkM

  ! MPI reduction within local group (infpar%lg_comm communicator)
    subroutine lsmpi_local_reduction_realkT(buffer,n1,n2,n3,master)
    implicit none
    integer(kind=ls_mpik) :: master
    integer :: n1,n2,n3
    real(realk) :: buffer(n1,n2,n3)
#ifdef VAR_MPI
    real(realk),pointer :: buf(:)
    integer(kind=8) :: n
    n=(i8*n1)*n2*n3
    call ass_8D3to1(buffer,buf,[i8*n1,i8*n2,i8*n3])
    call lsmpi_local_reduction_realkVN8(buf,n,master)
    nullify(buf)
#endif
  end subroutine lsmpi_local_reduction_realkT

  ! MPI reduction within local group (infpar%lg_comm communicator)
    subroutine lsmpi_local_reduction_realkQ(buffer,n1,n2,n3,n4,master)
    implicit none
    integer(kind=ls_mpik) :: master
    integer :: n1,n2,n3,n4
    real(realk) :: buffer(n1,n2,n3,n4)
#ifdef VAR_MPI
    real(realk),pointer :: buf(:)
    integer(kind=8) :: n
    n=((i8*n1)*n2)*n3*n4
    call ass_8D4to1(buffer,buf,[i8*n1,i8*n2,i8*n3,i8*n4])
    call lsmpi_local_reduction_realkVN8(buf,n,master)
    nullify(buf)
#endif
  end subroutine lsmpi_local_reduction_realkQ



  subroutine lsmpi_int4_reduction(buffer,n1,master,comm)
    implicit none
    integer(kind=ls_mpik) :: master,comm
    integer(kind=8) :: n1
    integer(kind=8) :: i,k
    integer(kind=ls_mpik) :: n4
    integer(kind=4) :: buffer(:)
#ifdef VAR_MPI
    integer(kind=4),pointer :: buffer2(:)
    integer(kind=ls_mpik) :: mynum
    integer(kind=ls_mpik) :: count,ierr
      IERR=0

    k=SPLIT_MPI_MSG

    call mem_alloc(buffer2,n1)
    call get_rank_for_comm(comm,mynum)

    IF(mynum.EQ.master)THEN
      buffer2=0
    ENDIF    
    COUNT = n1
    if(ls_mpik==4)then
      do i=1,n1,k
        n4=k
        if(((n1-i)<k).and.(mod(n1-i+1,k)/=0))n4=mod(n1,k)
        CALL MPI_REDUCE(BUFFER(i:i+n4-1),BUFFER2(i:i+n4-1),n4,MPI_INTEGER4,MPI_SUM,&
           &master,COMM,IERR)
        IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
      enddo
    else
      CALL MPI_REDUCE(BUFFER,BUFFER2,count,MPI_INTEGER4,MPI_SUM,&
          &master,COMM,IERR)
    endif

    IF(mynum.EQ.master)THEN
       DO I=1,n1
          buffer(I)=buffer2(I)
       ENDDO
    ENDIF    
    call mem_dealloc(buffer2)
#endif
  end subroutine lsmpi_int4_reduction

  subroutine lsmpi_int8_reduction(buffer,n1,master,comm)
    implicit none
    integer(kind=ls_mpik) :: master,comm
    integer(kind=8) :: n1
    integer(kind=8) :: i,k
    integer(kind=ls_mpik) :: n4
    integer(kind=8) :: buffer(:)
#ifdef VAR_MPI
    integer(kind=8),pointer :: buffer2(:)
    integer(kind=ls_mpik) :: mynum
    integer(kind=ls_mpik) :: count,ierr
      IERR=0

    k=SPLIT_MPI_MSG

    call mem_alloc(buffer2,n1)
    call get_rank_for_comm(comm,mynum)

    IF(mynum.EQ.master)THEN
      buffer2=0
    ENDIF    
    COUNT = n1
    if(ls_mpik==4)then
      do i=1,n1,k
        n4=k
        if(((n1-i)<k).and.(mod(n1-i+1,k)/=0))n4=mod(n1,k)
        CALL MPI_REDUCE(BUFFER(i:i+n4-1),BUFFER2(i:i+n4-1),n4,MPI_INTEGER8,MPI_SUM,&
           &master,COMM,IERR)
        IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
      enddo
    else
      CALL MPI_REDUCE(BUFFER,BUFFER2,count,MPI_INTEGER8,MPI_SUM,&
           &master,COMM,IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
    endif

    IF(mynum.EQ.master)THEN
       DO I=1,n1
          buffer(I)=buffer2(I)
       ENDDO
    ENDIF    
    call mem_dealloc(buffer2)
#endif
  end subroutine lsmpi_int8_reduction

  subroutine lsmpi_shoV_reduction_wrapper8(buffer,n1,master,comm)
    implicit none
    integer(kind=ls_mpik):: master,comm
    integer(kind=4)     :: n1
    integer(kind=8)     :: i,k
    integer(kind=short) :: buffer(*)
    integer(kind=8)     :: n
#ifdef VAR_MPI
    n=n1
    call lsmpi_shoV_reduction(buffer,n,master,comm)
#endif
  end subroutine lsmpi_shoV_reduction_wrapper8
  subroutine lsmpi_shoV_reduction(buffer,n,master,comm)
    implicit none
    integer(kind=ls_mpik):: master,comm
    integer(kind=8)     :: n
    integer(kind=8)     :: i,k
    integer(kind=short) :: buffer(*)
#ifdef VAR_MPI
#ifdef VAR_INT64
#ifdef VAR_MPI_32BIT_INT
    integer(kind=4),pointer :: buffer2(:),buffer1(:)
#else
    integer(kind=8),pointer :: buffer2(:),buffer1(:)
#endif
#else
    integer(kind=4),pointer :: buffer2(:),buffer1(:)
#endif
    integer(kind=ls_mpik) :: mynum,n4
    integer(kind=ls_mpik) :: count,ierr
      IERR=0
    k=SPLIT_MPI_MSG
    call get_rank_for_comm(comm,mynum)
    call mem_alloc(buffer1,n)
    call mem_alloc(buffer2,n)
    count=n
    DO I=1,n
       buffer1(I)=buffer(I)
    ENDDO
    IF(mynum.EQ.master)THEN
       DO I=1,n
          buffer2(I)=0
       ENDDO
    ENDIF  
    if(ls_mpik==4)then
      do i=1,n,k
        n4=k
        !if((n-i)<k)n4=mod(n,k)
        if(((n-i)<k).and.(mod(n-i+1,k)/=0))n4=mod(n,k)
        CALL MPI_REDUCE(BUFFER1(i:i+n4-1),BUFFER2(i:i+n4-1),n4,MPI_INTEGER,MPI_MAX,&
           &master,COMM,IERR)
        IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
      enddo
    else
      CALL MPI_REDUCE(BUFFER1,BUFFER2,count,MPI_INTEGER,MPI_MAX,&
           &master,COMM,IERR)
      IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
    endif
    IF(mynum.EQ.master)THEN
       DO I=1,n
          buffer(I)=buffer2(I)
       ENDDO
    ENDIF    
    call mem_dealloc(buffer1)
    call mem_dealloc(buffer2)
#endif
  end subroutine lsmpi_shoV_reduction

  subroutine lsmpi_max_int_reduction(buffer,master,comm)
    implicit none
    integer(kind=ls_mpik) :: master,comm
    integer :: buffer
#ifdef VAR_MPI
    integer :: buffer2
    integer(kind=ls_mpik) :: ierr,datatype,count,mynum
      IERR=0
    call get_rank_for_comm(comm,mynum)
    buffer2=0
    DATATYPE = MPI_INTEGER
    COUNT = 1
    CALL MPI_REDUCE(BUFFER,BUFFER2,count,MPI_INTEGER,MPI_MAX,&
         &master,COMM,IERR)
    IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
    IF(mynum.EQ.master)buffer = buffer2
#endif
  end subroutine lsmpi_max_int_reduction

  subroutine lsmpi_add_real_local_reduction(buffer,master)
    implicit none
    integer(kind=ls_mpik) :: master
    real(realk) :: buffer
#ifdef VAR_MPI
    real(realk) :: buffer2
    integer(kind=ls_mpik) :: ierr,datatype,count!,myop
      IERR=0
    buffer2=0.0E0_realk
    count=1
    DATATYPE = MPI_DOUBLE_PRECISION
    CALL MPI_REDUCE(BUFFER,BUFFER2,count,DATATYPE,MPI_SUM,&
         &master,infpar%lg_comm,IERR)
    IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
    IF(infpar%mynum.EQ.infpar%master)buffer = buffer2
#endif
  end subroutine lsmpi_add_real_local_reduction

  subroutine lsmpi_max_realk_reduction(buffer,master,comm)
    implicit none
    integer(kind=ls_mpik) :: master,comm
    real(realk) :: buffer
#ifdef VAR_MPI
    real(realk) :: buffer2
    integer(kind=ls_mpik) :: ierr,datatype,count,mynum!,myop
      IERR=0
    call get_rank_for_comm(comm,mynum)
    buffer2=0
    DATATYPE = MPI_INTEGER
    COUNT = 1
    CALL MPI_REDUCE(BUFFER,BUFFER2,count,MPI_DOUBLE_PRECISION,MPI_MAX,&
         &master,COMM,IERR)
    IF (IERR.GT. 0) CALL LSMPI_MYFAIL(IERR)
    IF(mynum.EQ.master)buffer = buffer2
#endif
  end subroutine lsmpi_max_realk_reduction

#ifdef VAR_MPI
  !> \brief Initialize MPI groups by setting information in the global
  !> type infpar. In particular, 
  !> infpar%lg_mynum, infpar%lg_nodtot, and infpar%lg_comm are set (see par_mod.f90).
  !> \author Kasper Kristensen
  !> \date March 2012
  subroutine init_mpi_groups(groupsize,lupri)
    implicit none
    !> Size of each group
    integer(kind=ls_mpik),intent(in) :: groupsize
    !> File unit number for LSDALTON.OUT
    integer,intent(in) :: lupri
    integer :: i
    integer(kind=ls_mpik) :: ierr,ngroups,worldgroup,localgroup
    integer(kind=ls_mpik) :: mygroup,hstatus
    type(mpigroup),allocatable :: lg(:)
    CHARACTER*(MPI_MAX_PROCESSOR_NAME) ::  hname
    integer(kind=ls_mpik) :: gmsize
    integer(kind=ls_mpik),pointer :: gmranks(:)
    IERR=0


    ! CONVENTION FOR CREATING GROUPS
    ! ******************************

    ! For now, the groups have equal size if the input groupsize is a multiplum of
    ! the total number of nodes disregarding the master node (infpar%nodtot-1).
    ! Otherwise the last local group will be correspondingly larger than the local groups.
    ! For simplicity, the master rank (0) gets special treatment and is not assigned to any
    ! local group.

    ! Example 1:  7 nodes in total (ranks 0,1,2,3,4,5,6); groupsize=2
    ! ---------------------------------------------------------------
    ! 
    ! There is one master (rank 0) and 3 local groups of size 2 (infpar%lg_nodtot=2) with ranks:
    ! First  local group: (1,2)    - rank numbers inside local group (infpar%lg_mynum): (0,1)
    ! Second local group: (3,4)    - rank numbers inside local group (infpar%lg_mynum): (0,1)
    ! Third  local group: (5,6)    - rank numbers inside local group (infpar%lg_mynum): (0,1)
    ! The local groups talk together via infpar%lg_comm.
    ! 
    ! Example 2:  8 nodes in total (ranks 0,1,2,3,4,5,6,7); groupsize=2
    ! -----------------------------------------------------------------
    !
    ! Same principle as above, except that now there will be three local groups 
    ! of size 2,2, and 3:
    ! First  local group: (1,2)    - rank numbers inside local group (infpar%lg_mynum): (0,1)
    ! Second local group: (3,4)    - rank numbers inside local group (infpar%lg_mynum): (0,1)
    ! Third  local group: (5,6,7)  - rank numbers inside local group (infpar%lg_mynum): (0,1,2)

    ! Sanity check 1: Size of group is smaller than total number of nodes
    ! (Maximum allowed group size is infpar%nodtot-1 because 1 node is reserved for master)
    ! Groupsize of course also have to be positive.
    if( (groupsize > infpar%nodtot-1) .or. (groupsize<1) ) then
       print *, 'Rank = ', infpar%mynum
       print *, 'Requested groupsize = ', groupsize
       print *, 'Number of nodes (excluding master)', infpar%nodtot-1
       call lsquit('init_mpi_groups: Requested group size is unacceptable!',-1)
    end if

    ! Sanity check 2: Groups only meaningful if there are least two nodes
    if(infpar%nodtot < 2) then
       print *, 'Rank = ', infpar%mynum
       print *, 'Number of nodes = ', infpar%nodtot
       call lsquit('init_mpi_groups: Cannot create groups - less than two nodes!',-1)
    end if


    ! For MASTER rank
    ! ===============
    if(infpar%mynum.EQ.infpar%master) then
       ! Wake up slaves to call this routine
       call ls_mpibcast(GROUPINIT,infpar%master,MPI_COMM_LSDALTON)

       ! Bcast groupsize to slave
       call ls_mpibcast(groupsize,infpar%master,MPI_COMM_LSDALTON)
    end if

    ! Extract group handle for LSDALTON group (world group)
    call MPI_COMM_GROUP(MPI_COMM_LSDALTON, worldgroup, ierr)


    ! **************************************************
    ! *                   Local groups                 *
    ! **************************************************

    ! Number of local groups
    ngroups = floor(real(infpar%nodtot-1)/real(groupsize))

    ! Local group structure
    allocate(lg(ngroups))
    call create_mpi_local_group_structure(groupsize,ngroups,mygroup,lg)

    ! For global master, only include master itself
    gmsize = 1
    call mem_alloc(gmranks,gmsize)
    gmranks=0

    ! Create local MPI group for this rank
    if(infpar%mynum /= infpar%master) then
       call MPI_GROUP_INCL(worldgroup, lg(mygroup)%groupsize, lg(mygroup)%ranks, localgroup, ierr)
    else
       ! Master is assigned its own group.
       call MPI_GROUP_INCL(worldgroup, gmsize, gmranks, localgroup, ierr)
    end if
    call mem_dealloc(gmranks)

    ! Create local MPI communicator
    call MPI_COMM_CREATE(MPI_COMM_LSDALTON, localgroup, infpar%lg_comm, ierr)

    ! Done with group information for local group (information is now contained in communicator)
    call MPI_GROUP_FREE(localgroup,ierr)

    if(infpar%mynum /= infpar%master) then

       ! Set global parameters for local group
       ! -------------------------------------

       ! number of nodes in group
       infpar%lg_nodtot = lg(mygroup)%groupsize   

       ! Rank within group
       call get_rank_for_comm(infpar%lg_comm,infpar%lg_mynum)

       ! Print out
       CALL MPI_GET_PROCESSOR_NAME(HNAME,HSTATUS,IERR)
       print '(a,4i6,1X,a)', 'World rank, local rank, group, groupsize, host:', infpar%mynum, &
            & infpar%lg_mynum, mygroup, lg(mygroup)%groupsize, hname

    else  ! Global master
       infpar%lg_nodtot = gmsize
       call get_rank_for_comm (infpar%lg_comm,infpar%lg_mynum)
    end if

    call MPI_GROUP_FREE(worldgroup,ierr)
    do i=1,ngroups
       deallocate(lg(i)%ranks)
       nullify(lg(i)%ranks)
    end do
    deallocate(lg)

  end subroutine init_mpi_groups


  !> \brief Create simple mpigroup structure for local groups (see init_mpi_groups for details).
  !> \author Kasper Kristensen
  !> \date March 2012
  subroutine create_mpi_local_group_structure(groupsize,ngroups,mygroup,lg)
    implicit none

    !> Requested group size (the last group may be larger than this, see init_mpi_groups)
    integer(kind=ls_mpik),intent(in) :: groupsize
    !> Number of groups
    integer(kind=ls_mpik),intent(in) :: ngroups
    !> Group number for rank under consideration
    integer(kind=ls_mpik),intent(inout) :: mygroup
    !> Local group information (size and list of ranks in each group)
    type(mpigroup),intent(inout) :: lg(ngroups)
    integer(kind=ls_mpik) :: i,j
    integer(kind=ls_mpik) :: count


    ! Loop over groups
    count=0
    mygroup=-1
    do i=1,ngroups
       if(i.EQ.ngroups) then
          ! Assign all the remaining nodes (except master) to the same local group
          lg(i)%groupsize = infpar%nodtot-1 - (ngroups-1)*groupsize
       else
          ! Input group size
          lg(i)%groupsize = groupsize
       end if
       nullify(lg(i)%ranks)
       allocate(lg(i)%ranks(lg(i)%groupsize))

       ! Define ranks in local group (see init_mpi_groups)
       do j=1,lg(i)%groupsize
          count = count+1
          lg(i)%ranks(j) = count
          if(count.EQ.infpar%mynum) mygroup=i
       end do

    end do

    ! Master (count=0) is not included in any local group
    if(infpar%mynum.EQ.infpar%master) mygroup=0

    ! Check
    if(count /= infpar%nodtot-1) then
       call lsquit('create_mpi_group_structure: &
            & Something wrong with accounting in group creation ',-1)
    end if
    if(mygroup.EQ.-1) then
       print *, 'rank = ', infpar%mynum
       call lsquit('create_mpi_group_structure: &
            & Group number for rank not identified!',-1)
    end if


  end subroutine create_mpi_local_group_structure



  ! \brief Print information about MPI master group and local group.
  !> Only intended to be called from the main master rank.
  ! \author Kasper Kristensen
  ! \date March 2012
  subroutine print_mpi_group_info(ngroups,lg,lupri)

    implicit none

    !> Number of groups
    integer,intent(in) :: ngroups
    !> Local group information 
    type(mpigroup),intent(in) :: lg(ngroups)
    !> File unit number for LSDALTON.OUT
    integer,intent(in) :: lupri
    integer :: j

    write(lupri,*) 
    write(lupri,*) 
    write(lupri,*) '********************************************************'
    write(lupri,*) '*                  MPI group information               *'
    write(lupri,*) '********************************************************'
    write(lupri,'(1X,a,i6,a)') 'There are ', ngroups, ' local groups'
    write(lupri,*) 
    write(lupri,*) 
    write(lupri,*) 'LOCAL GROUPS'
    write(lupri,*) '============'
    do j=1,ngroups
       write(lupri,*) 
       write(lupri,'(1X,a,i6)') 'Group:', j
       write(lupri,*) '--------------------'
       write(lupri,'(1X,a,i6)') 'Groupsize:', lg(j)%groupsize
       write(lupri,'(1X,a,1000i6)') 'Ranks:', lg(j)%ranks
    end do

  end subroutine print_mpi_group_info

  !> Get rank number within a specific communicator
  !> \author Kasper Kristensen
  !> \date March 2011
  subroutine get_rank_for_comm(comm,rank)
    implicit none
    !> Communicator 
    integer(kind=ls_mpik),intent(in) :: comm
    !> Rank number in communicator group
    integer(kind=ls_mpik),intent(inout) :: rank
    integer(kind=ls_mpik) :: ierr
    ierr=0
    call MPI_COMM_RANK(comm,rank,ierr)
    if(ierr/=0) then
       call lsquit('get_rank_for_comm: Something wrong!',-1)
    end if
  end subroutine get_rank_for_comm

  !> Get number of processors for a specific communicator
  !> \author Thomas Kjaergaard
  !> \date juni 2012
  subroutine get_size_for_comm(comm,nodtot)
    implicit none
    !> Communicator
    integer(kind=ls_mpik),intent(in) :: comm
    !> Rank number in communicator group
    integer(kind=ls_mpik),intent(inout) :: nodtot
    integer(kind=ls_mpik) :: ierr
    ierr=0
    call MPI_COMM_SIZE(comm,nodtot,ierr)
    if(ierr/=0) then
       call lsquit('get_size_for_comm: Something wrong!',-1)
    end if
  end subroutine get_size_for_comm



  !> \brief Divide local MPI group into ngroups smaller of dimensions defined in input.
  !> To do this the global parameters infpar%lg_mynum, infpar%lg_nodtot, and infpar%lg_comm 
  !> are redefined here.
  !> It is important to ensure that ALL MEMBERS of the local group call
  !> this routine at the same time - if not, some ranks will be waiting for a signal
  !> which never comes...
  !> In other words, synchronization of all nodes in the local groups must be done OUTSIDE this
  !> routine because their communication channel is redefined in here!
  !> \author Kasper Kristensen
  !> \date May 2012
  subroutine divide_local_mpi_group(ngroups,groupdims)
    implicit none

    !> Number of new local groups to create from existing local group
    integer(kind=ls_mpik),intent(in) :: ngroups
    !> Dimensions for these groups
    !> NOTE: The sum of these dimensions must equal the current group size:
    !> sum(groupsize) = infpar%lg_nodtot
    integer(kind=ls_mpik),dimension(ngroups),intent(in) :: groupdims
    integer(kind=ls_mpik) :: old_nodtot, old_mynum, oldgroup, ierr, newgroup,I
    integer(kind=ls_mpik), allocatable :: mygroup(:)
    integer(kind=ls_mpik) :: mysize,offset,counter,mygroupnumber,newcomm
      IERR=0

    ! EXAMPLE
    ! *******
    ! Assume that the current local group size is 7 and that we want
    ! to create two new groups of sizes 4 and 3:
    !
    ! ngroups = 2
    ! groupdims(1) = 4
    ! groupdims(2) = 3
    !
    ! The ranks are then reassigned to local groups as follows:
    !    Current local rank number               New local rank number 
    !    0 (groupsize=7)                         0  (new subgroup 1, new groupsize=4)
    !    1 (groupsize=7)                         1  (new subgroup 1, new groupsize=4)
    !    2 (groupsize=7)                         2  (new subgroup 1, new groupsize=4)
    !    3 (groupsize=7)                         3  (new subgroup 1, new groupsize=4)
    !    4 (groupsize=7)                         0  (new subgroup 2, new groupsize=3)
    !    5 (groupsize=7)                         1  (new subgroup 2, new groupsize=3)
    !    6 (groupsize=7)                         2  (new subgroup 2, new groupsize=3)
    !
    ! Local rank number is stored in infpar%lg_mynum,
    ! and the local group size is stored in infpar%lg_nodtot.
    ! The local groups communicate vta infpar%lg_comm
    ! These three global parameters are all modified here to describe the new local group.


    ! Save existing local group information before overwriting
    old_mynum = infpar%lg_mynum
    old_nodtot = infpar%lg_nodtot

    ! Sanity check: Sum of new group sizes equals old group size
    if(sum(groupdims) /= old_nodtot) then
       print *, 'Old size = ', old_nodtot
       print *, 'New sizes= ', groupdims
       call lsquit('divide_local_mpi_group: Sum of new group sizes &
            & does not equal old group size!',-1)
    end if


    ! Determine new local group for this rank
    ! =======================================
    counter=0
    mysize=0
    do i=1,ngroups
       counter = counter + groupdims(i)
       if(old_mynum < counter) then  

          ! this rank belongs to new group i (see example above)
          mygroupnumber  = i   

          ! size of new group
          mysize = groupdims(i)  

          ! number of ranks in other new groups having smaller group number
          offset = counter - mysize
          exit

       end if
    end do

    if(mysize.EQ.0) call lsquit('Something wrong in divide local groups',-1)
    allocate(mygroup(mysize))

    ! Define new local group
    do i=1,mysize
       ! Ranks in new group defined as [offset;offset+mysize]
       ! and then we subtract 1 to start counting from rank 0 in old local group.
       mygroup(i) = offset + i - 1
    end do



    ! Create new group
    ! ================

    ! Extract group handle (oldgroup) for existing local group
    call MPI_COMM_GROUP(infpar%lg_comm, oldgroup, ierr)

    ! Create group handle (newgroup) for new local group
    call MPI_GROUP_INCL(oldgroup, mysize, mygroup, newgroup, ierr)

    call MPI_COMM_CREATE(infpar%lg_comm, newgroup, newcomm, ierr)


    ! Done with group information for local group (information is now in communicator, newcomm)
    call MPI_GROUP_FREE(newgroup,ierr)
    call MPI_GROUP_FREE(oldgroup,ierr)
    ! Free old communcator
    call MPI_COMM_FREE(infpar%lg_comm, IERR)


    ! Set global parameters for new local group (overwrite old ones)
    ! ==============================================================

    ! number of nodes in group
    infpar%lg_nodtot = mysize

    ! Rank within group
    call get_rank_for_comm(newcomm,infpar%lg_mynum)

    ! Communicator
    infpar%lg_comm = newcomm

    if(infpar%lg_mynum==0) then
       print *, 'Redefining local group'
       print *, '======================'
       print '(1X,a,4i6)', 'Old rank, new rank, old size, new size ', &
            & old_mynum, infpar%lg_mynum, old_nodtot, infpar%lg_nodtot
       print *, 'New group = ', mygroup
    end if

    deallocate(mygroup)

  end subroutine divide_local_mpi_group
#endif


    subroutine lsmpi_print(lupri)
    implicit none
    integer :: lupri    
#ifdef VAR_MPI
    write(lupri,'(3X,A,I7,A)')'This is an MPI run using ',infpar%nodtot,' processes.'
#endif
  end subroutine lsmpi_print

    !> \brief Set rank, sizes and communcators for local groups
    !> to be identical to those for world group.
    !> \author Kasper Kristensen
    !> \date March 2012
    subroutine lsmpi_default_mpi_group

#ifdef VAR_MPI
      implicit none

      infpar%lg_comm = MPI_COMM_LSDALTON
      infpar%lg_mynum = infpar%mynum
      infpar%lg_nodtot = infpar%nodtot      
#endif

    end subroutine lsmpi_default_mpi_group


    !> \brief calling this routine from the master process will start up
    !communication processes on all the slaves
    !> \author Patrick Ettenhuber
    !> \date 2013
    subroutine local_group_start_comm_processes
      implicit none
#ifdef VAR_MPI
      print *,"STARTING UP THE COMMUNICATION PROCESSES"
      !impregnate the slaves
      call ls_mpibcast(GIVE_BIRTH,infpar%master,infpar%lg_comm)
      !just to find, that the master is also fertile
      call give_birth_to_child_process
#else
      call lsquit("ERROR(local_group_start_comm_processes): not available without mpi",-1)
#endif
    end subroutine local_group_start_comm_processes

    !> \brief calling this routine from the master process will kill all
    !communication processes on all the slaves
    !> \author Patrick Ettenhuber
    !> \date 2013
    subroutine local_group_kill_comm_processes
      implicit none
#ifdef VAR_MPI
      print *,"SHUTTING DOWN THE COMMUNICATION PROCESSES"
      !kill the babies of the slaves
      call ls_mpibcast(SLAVES_SHUT_DOWN_CHILD,infpar%master,infpar%lg_comm)
      !kill own baby
      call shut_down_child_process
#else
      call lsquit("ERROR(local_group_kill_comm_processes): not available without mpi",-1)
#endif
    end subroutine local_group_kill_comm_processes

    !> \brief make a new communicator for the child and parent process(es)
    !> \author Patrick Ettenhuber
    !> \date 2013
    subroutine get_parent_child_relation
      implicit none
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: ierr
      logical(kind=ls_mpik) :: last

      if( infpar%parent_comm == MPI_COMM_NULL ) then
        ! if i am a parent 
        last = .false.
        call MPI_INTERCOMM_MERGE( infpar%child_comm, last, infpar%pc_comm, ierr )
      else
        ! if i am a child
        last = .true.
        call MPI_INTERCOMM_MERGE( infpar%parent_comm, last, infpar%pc_comm, ierr )
      endif

      call get_rank_for_comm( infpar%pc_comm, infpar%pc_mynum  )
      call get_size_for_comm( infpar%pc_comm, infpar%pc_nodtot )

      lsmpi_enabled_comm_procs = .true.

      if( infpar%parent_comm == MPI_COMM_NULL .and.  infpar%pc_mynum/=infpar%master )&
      & call lsquit("ERROR(get_parent_child_relation)&
      & parent needs to have rank 0 in the pc_comm",-1)

      if( infpar%parent_comm /= MPI_COMM_NULL .and.  infpar%pc_mynum==infpar%master )&
      & call lsquit("ERROR(get_parent_child_relation)&
      & child cannot have rank 0 in the pc_comm",-1)
#endif
    end subroutine get_parent_child_relation


    subroutine give_birth_to_child_process
      implicit none
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: procs_to_spawn,root,errorcode(1),ierr
      logical               :: localdalton

          procs_to_spawn = int(1,kind=ls_mpik)
          root           = int(0,kind=ls_mpik)

          !check if the program is in the workdir, if not, spawning is not
          !possible
          inquire(file='lsdalton.x', exist=localdalton)

          !if lsdalton.x is in workdir, spawn another process
          if(localdalton)then
            call MPI_COMM_SPAWN('./lsdalton.x',MPI_ARGV_NULL,procs_to_spawn,MPI_INFO_NULL,&
             &root,MPI_COMM_SELF,infpar%child_comm,errorcode,ierr)
            call get_parent_child_relation
          else
            call lsquit("ERROR(give_birth_to_child_process):lsdalton.x was not&
             &found in the working directory, move it there and restart",-1)
          endif

#endif
    end subroutine give_birth_to_child_process


    !> \brief free everything related to child processes
    !> \author Patrick Ettenhuber
    !> \date 2013
    subroutine shut_down_child_process
      implicit none
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: ierr
      logical(kind=ls_mpik) :: have_priority

      if( infpar%parent_comm == MPI_COMM_NULL ) then
        call MPI_COMM_FREE(infpar%pc_comm,ierr)
        call MPI_COMM_FREE(infpar%child_comm,ierr)
      else
        ! if i am a child
        call MPI_COMM_FREE(infpar%pc_comm,ierr)
        call MPI_COMM_FREE(infpar%parent_comm,ierr)
      endif

      lsmpi_enabled_comm_procs = .false.

#endif
    end subroutine shut_down_child_process

    subroutine lsmpi_print_mem_info(lupri,mastercall)
      implicit none
      integer,intent(in) :: lupri
      logical,intent(in) :: mastercall
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: o,n,ierr,I,t(1),tag_meminfo,count,dest,tag,from,root
      integer(kind=long) :: recvbuffer
      real(realk) :: recvbuffer_real
      integer(kind=long),pointer :: longintbufferInt(:) 
#ifdef VAR_INT64
      integer, parameter :: i2l = 1
#else
      integer, parameter :: i2l = 2
#endif
      IERR=0

      tag=2001
      dest=0
      t=0
      o=1;n=0

      if ((infpar%mynum.eq.infpar%master).and.mastercall) THEN
         !wake up slaves
         call ls_mpibcast(LSMPIPRINTINFO,infpar%master,MPI_COMM_LSDALTON)
      ENDIF

      count = 1
      root = 0
      recvbuffer = 0
      !Total max_mem_used_global across all nodes
      CALL MPI_REDUCE(max_mem_used_global,recvbuffer,&
           & count,MPI_INTEGER8,MPI_SUM,root,MPI_COMM_LSDALTON,IERR)
      !IF direct communication is used the unlock times are interesting

      call lsmpi_reduction(time_win_unlock,infpar%master,MPI_COMM_LSDALTON)
  
      IF(infpar%mynum.eq.infpar%master) THEN
         WRITE(lupri,'(A)')'  The total memory used across all MPI nodes'
         call print_maxmem(lupri,recvbuffer,'TOTAL')
         Write(lupri,'("time spent in unlock: ",f19.10)')time_win_unlock
      ENDIF
      !Largest max_mem_used_global including Master
      CALL MPI_REDUCE(max_mem_used_global,recvbuffer,&
           & count,MPI_INTEGER8,MPI_MAX,root,MPI_COMM_LSDALTON,IERR)
      IF(infpar%mynum.eq.infpar%master) THEN
         WRITE(lupri,'(A)')'  Largest memory allocated on a single MPI slave node (including Master)'
         call print_maxmem(lupri,recvbuffer,'TOTAL')
      ENDIF
      !Largest Slave max_mem_used_global (not including Master)      
      IF(infpar%mynum.EQ.infpar%master)max_mem_used_global=0
      CALL MPI_REDUCE(max_mem_used_global,recvbuffer,&
           & count,MPI_INTEGER8,MPI_MAX,root,MPI_COMM_LSDALTON,IERR)
      IF(infpar%mynum.eq.infpar%master) THEN
         WRITE(lupri,'(A)')'  Largest memory allocated on a single MPI slave node (exclusing Master)'
         call print_maxmem(lupri,recvbuffer,'TOTAL')
      ENDIF


      IF(infpar%mynum.NE.infpar%master)then
         call stats_mem(6)
         call mem_alloc(longintbufferInt,longintbuffersize)
         call copy_from_mem_stats(longintbufferInt)
         count=longintbuffersize
         CALL MPI_SEND(longintbufferInt,count,MPI_INTEGER8,dest,tag,MPI_COMM_LSDALTON,IERR)
         call mem_dealloc(longintbufferInt)

      ELSE          
         DO I=1,infpar%nodtot-1
            call init_globalmemvar !WARNING removes all mem info on master 
            call mem_alloc(longintbufferInt,longintbuffersize)
            count = longintbuffersize
            CALL MPI_RECV(longintbufferInt,count,MPI_INTEGER8,I,tag,&
                 & MPI_COMM_LSDALTON,status,IERR)
            call copy_to_mem_stats(longintbufferInt)
            call mem_dealloc(longintbufferInt)
            WRITE(LUPRI,'("  Memory statistics for MPI node number ",i9," ")') I
            call stats_mpi_mem(lupri)
         ENDDO
      ENDIF

!!$!#ifdef VAR_LSDEBUG
!!$!     count=1
!!$!     root = infpar%master
!!$!     recvbuffer = 0
!!$!     recvbuffer_real = 0.0E0_realk
!!$!     CALL MPI_REDUCE(poketime,recvbuffer_real,&
!!$!            & count,MPI_DOUBLE_PRECISION,MPI_SUM,root,MPI_COMM_LSDALTON,IERR)
!!$!     CALL MPI_REDUCE(poketimes,recvbuffer,&
!!$!            & count,MPI_INTEGER8,MPI_SUM,root,MPI_COMM_LSDALTON,IERR)
!!$!     if(infpar%mynum==infpar%master)then
!!$!       print *,"CUMULATIVE MPI POKETIME",recvbuffer_real,recvbuffer,recvbuffer_real/(recvbuffer*1.0E0_realk)
!!$!     endif
!!$!#endif     
#endif
    end subroutine lsmpi_print_mem_info

    subroutine lsmpi_finalize(lupri,mastercall)
      implicit none
      integer,intent(in)    :: lupri
      logical,intent(in)     :: mastercall
#ifdef VAR_MPI
      integer(kind=ls_mpik) :: ierr
      ierr = 0

      if ((infpar%mynum.eq.infpar%master).and.mastercall.and.(infpar%parent_comm==MPI_COMM_NULL))&
       &call ls_mpibcast(LSMPIQUIT,infpar%master,MPI_COMM_LSDALTON)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !CHECK IF EVERYTHING IS OKAY INSIDE THE MODULE!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if(lsmpi_enabled_comm_procs)call lsquit("ERROR(lsmpi_finalize):&
      &comm processes were enabled on shutdown, this should never occur",-1)

#ifdef VAR_CHEMSHELL
      ! jump out of LSDALTON if a slave (instead of STOP)
      if (infpar%mynum.ne.infpar%master) call lsdaltonjumpout(99)
#else

      call MPI_FINALIZE(ierr)

      if(ierr/=0)then
        write (*,*) "mpi_finalize returned",ierr
        call LSMPI_MYFAIL(ierr)
        call lsquit("ERROR(MPI_FINALIZE):non zero exit)",-1)
      endif

#endif 
#endif 
    end subroutine lsmpi_finalize


    subroutine lsmpi_barrier(comm)
#ifdef VAR_MPI
    implicit none
    integer(kind=ls_mpik) :: ierr,comm
    IERR=0
    call MPI_BARRIER(comm,ierr)
#endif 
  end subroutine lsmpi_barrier
 

  !> \biref Slave routine for initializing MPI groups.
  !> Get groupsize from master and calls main group initialization routine.
  !> \author Kasper Kristensen
  !> \date March 2012
  subroutine init_mpi_groups_slave
#ifdef VAR_MPI
    implicit none
    integer(kind=ls_mpik) :: groupsize
    integer :: lupri

    ! Set output to standard output (currently not used for slave anyway)
    lupri=6

    ! Get groupsize from master
    call ls_mpibcast(groupsize,infpar%master,MPI_COMM_LSDALTON)
    ! Main group initialization routine
    call init_mpi_groups(groupsize,lupri)

#endif

  end subroutine init_mpi_groups_slave



!############################################################
!#
!# SPECIAL ROUTINES FOR ONE-SIDED COMMUNICATION 
!#  
!############################################################


  !> \brief simple mpi_win creation
  !> \author Patrick Ettenhuber
  !> \date September 2012
  subroutine lsmpi_win_create_realk8(darr,win,nel,comm)
    implicit none
    integer(kind=8), intent(in) :: nel
    real(realk),intent(in) :: darr(nel)
    integer(kind=ls_mpik),intent(inout) :: win,comm
#ifdef VAR_MPI
    integer(kind=ls_mpik) :: ierr,info,rk_len
    integer(kind=MPI_ADDRESS_KIND) :: mpi_realk,lb,bytes
    IERR=0
    info = MPI_INFO_NULL
    call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,mpi_realk,ierr)
    if(ierr.ne.0)then
      call lsquit("Error(lsmpi_localwin_create_realk)",ierr)
    endif
    bytes  = nel*mpi_realk
    rk_len = int(mpi_realk,kind=ls_mpik)
    call MPI_WIN_CREATE(darr,bytes,rk_len,info,comm,win,ierr)
    if(ierr.ne.0)then
      call lsquit("Error(lsmpi_localwin_create_realk)",ierr)
    endif
#endif
  end subroutine lsmpi_win_create_realk8
  !> \brief simple mpi_win creation
  !> \author Patrick Ettenhuber
  !> \date September 2012
  subroutine lsmpi_win_create_realk4(darr,win,nel,comm)
    implicit none
    integer(kind=4), intent(in) :: nel
    real(realk),intent(in) :: darr(nel)
    integer(kind=ls_mpik),intent(inout) :: win,comm
#ifdef VAR_MPI
    integer(kind=ls_mpik) :: ierr,info,rk_len
    integer(kind=MPI_ADDRESS_KIND) :: mpi_realk,lb,bytes
    IERR=0
    info = MPI_INFO_NULL
    call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,mpi_realk,ierr)
    if(ierr.ne.0)then
      call lsquit("Error(lsmpi_localwin_create_realk)",ierr)
    endif
    bytes  = nel*mpi_realk
    rk_len = int(mpi_realk,kind=ls_mpik)
    call MPI_WIN_CREATE(darr,bytes,rk_len,info,comm,win,ierr)
    if(ierr.ne.0)then
      call lsquit("Error(lsmpi_localwin_create_realk)",ierr)
    endif
#endif
  end subroutine lsmpi_win_create_realk4
  subroutine lsmpi_win_create_int8(iarr,win,nel,comm)
    implicit none
    integer(kind=8),intent(in) :: iarr(nel)
    integer(kind=ls_mpik),intent(inout) :: win
    integer, intent(in) :: nel
    integer(kind=ls_mpik) :: comm
#ifdef VAR_MPI
    integer(kind=ls_mpik) :: ierr,info,intlen
    integer(kind=MPI_ADDRESS_KIND) :: mpi_intlen,lb,bytes
    IERR=0

    info = MPI_INFO_NULL
    call MPI_TYPE_GET_EXTENT(MPI_INTEGER8,lb,mpi_intlen,ierr)
    if(ierr.ne.0)then
      call lsquit("Error(lsmpi_localwin_create_int8)",ierr)
    endif
    bytes  = nel*mpi_intlen
    intlen = int(mpi_intlen,kind=ls_mpik)
    call MPI_WIN_CREATE(iarr,bytes,intlen,info,comm,win,ierr)
    if(ierr.ne.0)then
      call lsquit("Error(lsmpi_localwin_create_int8)",ierr)
    endif
#endif
  end subroutine lsmpi_win_create_int8
  subroutine lsmpi_win_create_int4(iarr,win,nel,comm)
    implicit none
    integer(kind=4),intent(in) :: iarr(nel)
    integer(kind=ls_mpik),intent(inout) :: win
    integer, intent(in) :: nel
    integer(kind=ls_mpik) :: comm
#ifdef VAR_MPI
    integer(kind=ls_mpik) :: ierr,info,intlen
    integer(kind=MPI_ADDRESS_KIND) :: mpi_intlen,lb,bytes
    IERR=0

    info = MPI_INFO_NULL
    call MPI_TYPE_GET_EXTENT(MPI_INTEGER4,lb,mpi_intlen,ierr)
    if(ierr.ne.0)then
      call lsquit("Error(lsmpi_localwin_create_int4)",ierr)
    endif
    bytes  = nel*mpi_intlen
    intlen = int(mpi_intlen,kind=ls_mpik)
    call MPI_WIN_CREATE(iarr,bytes,intlen,info,comm,win,ierr)
    if(ierr.ne.0)then
      call lsquit("Error(lsmpi_localwin_create_int4)",ierr)
    endif
#endif
  end subroutine lsmpi_win_create_int4

  !> \brief simple mpi_win_free
  !> \author Patrick Ettenhuber
  !> \date September 2012
  subroutine lsmpi_win_free(win)
    implicit none
    integer(kind=ls_mpik),intent(inout) :: win
#ifdef VAR_MPI
    integer(kind=ls_mpik) :: ierr
    IERR=0
    call MPI_WIN_FREE(win,ierr)
    if(ierr.ne.0)then
      call lsquit("Error(lsmpi_win_free)",ierr)
    endif
#endif
  end subroutine lsmpi_win_free

  !> \brief simple mpi_win_fence call for dma/rma
  !> \author Patrick Ettenhuber
  !> \date September 2012
  subroutine lsmpi_win_fence_simple(win)
    implicit none
    integer(kind=ls_mpik),intent(in) :: win
#ifdef VAR_MPI
    integer(kind=ls_mpik) :: ierr,nr
    IERR=0
    nr=0
    call MPI_WIN_FENCE(nr,win,ierr)
    if(ierr.ne.0)then
      call lsquit("Error(lsmpi_win_fence)",ierr)
    endif
#endif
  end subroutine lsmpi_win_fence_simple

  !> \brief simple call to  mpi_win_fence for intitialization
  !> \author Patrick Ettenhuber
  !> \date Januar 2013
  subroutine lsmpi_win_fence_special(win,openwin)
    implicit none
    integer(kind=ls_mpik),intent(in) :: win
    logical,intent(in) :: openwin
#ifdef VAR_MPI
    integer(kind=ls_mpik) :: ierr
    IERR=0
    if(openwin)then
      call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE,win,ierr)
      if(ierr.ne.0)then
        call lsquit("Error(lsmpi_first_fence)",ierr)
      endif
    else
      IERR=0
      call MPI_WIN_FENCE(MPI_MODE_NOSUCCEED,win,ierr)
      if(ierr.ne.0)then
        call lsquit("Error(lsmpi_last_fence)",ierr)
      endif
    endif
#endif
  end subroutine lsmpi_win_fence_special

  subroutine lsmpi_win_lock(dest,win,typeoflock,ass)
    implicit none
    integer(kind=ls_mpik),intent(in) :: win
    integer(kind=ls_mpik),intent(in) :: dest
    integer(kind=ls_mpik),intent(in),optional :: ass
    character, intent(in) :: typeoflock
#ifdef VAR_MPI   
    integer(kind=ls_mpik) :: ierr, assert
    assert = 0
    if(present(ass))assert=ass
    ierr = 0
    if(typeoflock=='e')then
      CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,dest,assert,win,ierr)
    else if(typeoflock=='s')then
      CALL MPI_WIN_LOCK(MPI_LOCK_SHARED,dest,assert,win,ierr)
    else
      call lsquit("ERROR(lsmpi_win_lock): no valid lock type selected",-1)
    endif

    if(ierr.ne.0)then
      call lsquit("Error(lsmpi_win_lock): error in mpi",ierr)
    endif
#endif
  end subroutine lsmpi_win_lock
  subroutine lsmpi_win_unlock(dest,win)
    implicit none
    integer(kind=ls_mpik),intent(in) :: dest
    integer(kind=ls_mpik),intent(in) :: win
#ifdef VAR_MPI   
    integer(kind=ls_mpik) :: ierr
    real(realk) :: ta,te
    ierr=0
    ta = MPI_WTIME()
    call MPI_WIN_UNLOCK(dest,win,ierr)     
    te = MPI_WTIME()
    if(ierr.ne.0)then
      call lsquit("Error(lsmpi_win_unlock): error in mpi",ierr)
    endif
    time_win_unlock = time_win_unlock + te - ta
#endif
  end subroutine lsmpi_win_unlock

  subroutine lsmpi_win_flush(win,rank,local)
     implicit none
     integer(kind=ls_mpik) :: win
     integer(kind=ls_mpik), optional :: rank
     logical, optional :: local
#ifdef VAR_MPI   
#ifdef VAR_HAVE_MPI3
     integer(kind=ls_mpik) :: ierr
     logical :: loc
     loc  = .false.
     ierr = 0

     if(present(local))loc = local

     if(loc)then
        if(present(rank))then
           call MPI_WIN_FLUSH_LOCAL(rank,win,ierr)
        else
           call MPI_WIN_FLUSH_LOCAL_ALL(win,ierr)
        endif
     else
        if(present(rank))then
           call MPI_WIN_FLUSH(rank,win,ierr)
        else
           call MPI_WIN_FLUSH_ALL(win,ierr)
        endif
     endif

     if(ierr /= 0_ls_mpik)then
        call lsquit("ERROR(lsmpi_win_flush): non zero exit, error in mpi",-1)
     endif
#else
     print *,"WARNING(lsmpi_win_flush)called without MPI3, unlock should force&
     & the sync"
#endif
#endif
  end subroutine lsmpi_win_flush

  !=========================================================!
  !                   MPI PUT ROUTINES                      !
  !=========================================================!
  subroutine lsmpi_put_realk(buf,pos,dest,win)
    implicit none
    real(realk),intent(in) :: buf
    integer, intent(in) :: pos
    integer(kind=ls_mpik),intent(in) :: dest
    integer(kind=ls_mpik),intent(in) :: win
#ifdef VAR_MPI
    integer(kind=ls_mpik) :: n,ierr
    integer(kind=MPI_ADDRESS_KIND) :: offset
    ierr = 0
    n  = 1
    offset = int(pos-1,kind=MPI_ADDRESS_KIND)
    call MPI_PUT(buf,n,MPI_DOUBLE_PRECISION,dest, &
     & offset,n,MPI_DOUBLE_PRECISION,win,ierr)
    if(ierr.ne.0)then
      call lsquit("Error(lsmpi_put_realk)",ierr)
    endif
#endif
  end subroutine lsmpi_put_realk
  subroutine lsmpi_put_realkV_wrapper8(buf,nelms,pos,dest,win)
    implicit none
    real(realk),intent(in) :: buf(*)
    integer, intent(in) :: pos
    integer(kind=8) :: nelms
    integer(kind=ls_mpik),intent(in) :: dest
    integer(kind=ls_mpik),intent(in) :: win
#ifdef VAR_MPI
    integer(kind=ls_mpik) :: n,ierr
    integer(kind=4) :: n4,k,i
    integer(kind=MPI_ADDRESS_KIND) :: offset
    if(ls_mpik==4)then
      k=SPLIT_MPI_MSG
      do i=1,nelms,k
        n4=k
        if(((nelms-i)<k).and.(mod(nelms-i+1,k)/=0))n4=mod(nelms,k)
        call lsmpi_put_realkV(buf(i:i+n4-1),n4,pos+i-1,dest,win)
      enddo
    else
      ierr = 0
      n  = nelms
      offset = int(pos-1,kind=MPI_ADDRESS_KIND)
      call MPI_PUT(buf,n,MPI_DOUBLE_PRECISION,dest, &
       & offset,n,MPI_DOUBLE_PRECISION,win,ierr)
      if(ierr.ne.0)then
        call lsquit("Error(lsmpi_put_realk)",ierr)
      endif
    endif
#endif
  end subroutine lsmpi_put_realkV_wrapper8
  subroutine lsmpi_put_realkV(buf,nelms,pos,dest,win)
    implicit none
    real(realk),intent(in) :: buf(*)
    integer, intent(in) :: pos
    integer(kind=4) :: nelms
    integer(kind=ls_mpik),intent(in) :: dest
    integer(kind=ls_mpik),intent(in) :: win
#ifdef VAR_MPI
    integer(kind=ls_mpik) :: n,ierr
    integer(kind=MPI_ADDRESS_KIND) :: offset
    ierr = 0
    n  = nelms
    offset = int(pos-1,kind=MPI_ADDRESS_KIND)
    call MPI_PUT(buf,n,MPI_DOUBLE_PRECISION,dest, &
     & offset,n,MPI_DOUBLE_PRECISION,win,ierr)
    if(ierr.ne.0)then
      call lsquit("Error(lsmpi_put_realk)",ierr)
    endif
#endif
  end subroutine lsmpi_put_realkV
  subroutine lsmpi_put_realkV_parts_wrapper8(buf,nelms,pos,dest,win,batchsze,flush_it)
    implicit none
    real(realk),intent(in) :: buf(*)
    integer, intent(in) :: pos
    integer(kind=8) :: nelms
    integer(kind=ls_mpik),intent(in) :: dest
    integer(kind=ls_mpik),intent(in) :: win
    integer, intent(in) :: batchsze
    logical, intent(in), optional :: flush_it
#ifdef VAR_MPI
    integer :: newpos
    integer(kind=4) :: n4,k,i
    integer(kind=8) :: n,j
    logical :: fi

    fi = .false.
    if(present(flush_it))fi = flush_it

    if(ls_mpik==4)then

      k=SPLIT_MPI_MSG

      do i=1,nelms,k

        n4=k
        if(((nelms-i)<k).and.(mod(nelms-i+1,k)/=0))n4=mod(nelms,k)
        call lsmpi_put_realkV_parts(buf(i:i+n4-1),n4,pos+i-1,dest,win,batchsze,flush_it = flush_it)

      enddo

    else

      do j=1,nelms,batchsze

        n=batchsze

        if(((nelms-j)<batchsze).and.&
          &(mod(nelms-j+1,batchsze)/=0))n=mod(nelms,batchsze)

        newpos = pos+j-1

        call lsmpi_put_realkV_wrapper8(buf(j:j+n-1),n,newpos,dest,win)

#ifdef VAR_HAVE_MPI3
        if(fi)call lsmpi_win_flush(win,rank=dest,local=.true.)
#endif
      enddo

    endif
#endif
  end subroutine lsmpi_put_realkV_parts_wrapper8
  subroutine lsmpi_put_realkV_parts(buf,nelms,pos,dest,win,batchsze,flush_it)
    implicit none
    real(realk),intent(in) :: buf(*)
    integer, intent(in) :: pos
    integer(kind=4),intent(in) :: nelms
    integer(kind=ls_mpik),intent(in) :: dest
    integer(kind=ls_mpik),intent(in) :: win
    integer, intent(in) :: batchsze
    logical, intent(in), optional :: flush_it
#ifdef VAR_MPI
    integer :: newpos
    integer(kind=8) :: n,i
    logical :: fi

    fi = .false.
    if(present(flush_it))fi = flush_it

    do i=1,nelms,batchsze

      n=batchsze

      if(((nelms-i)<batchsze).and.&
        &(mod(nelms-i+1,batchsze)/=0))n=mod(nelms,batchsze)

      newpos = pos+i-1

      call lsmpi_put_realkV_wrapper8(buf(i:i+n-1),n,newpos,dest,win)

#ifdef VAR_HAVE_MPI3
      if(fi)call lsmpi_win_flush(win,rank=dest,local=.true.)
#endif

    enddo
#endif
  end subroutine lsmpi_put_realkV_parts


  !=========================================================!
  !                   MPI GET ROUTINES                      !
  !=========================================================!
  subroutine lsmpi_get_int8(buf,pos,dest,win)
    implicit none
    integer(kind=8),intent(in) :: buf
    integer, intent(in) :: pos
    integer(kind=ls_mpik),intent(in) :: dest
    integer(kind=ls_mpik),intent(in) :: win
#ifdef VAR_MPI
    integer(kind=ls_mpik) :: n,ierr
    integer(kind=MPI_ADDRESS_KIND) :: offset
    ierr = 0
    n  = 1
    offset = int(pos-1,kind=MPI_ADDRESS_KIND)
    call MPI_GET(buf,n,MPI_INTEGER8,dest, &
     & offset,n,MPI_INTEGER8,win,ierr)
    if(ierr.ne.0)then
      call lsquit("Error(lsmpi_get_int8)",ierr)
    endif
#endif
  end subroutine lsmpi_get_int8
  subroutine lsmpi_get_int4(buf,pos,dest,win)
    implicit none
    integer(kind=4),intent(in) :: buf
    integer, intent(in) :: pos
    integer(kind=ls_mpik),intent(in) :: dest
    integer(kind=ls_mpik),intent(in) :: win
#ifdef VAR_MPI
    integer(kind=ls_mpik) :: n,ierr
    integer(kind=MPI_ADDRESS_KIND) :: offset
    ierr = 0
    n  = 1
    offset = int(pos-1,kind=MPI_ADDRESS_KIND)
    call MPI_GET(buf,n,MPI_INTEGER4,dest, &
     & offset,n,MPI_INTEGER4,win,ierr)
    if(ierr.ne.0)then
      call lsquit("Error(lsmpi_get_int4)",ierr)
    endif
#endif
  end subroutine lsmpi_get_int4
  subroutine lsmpi_get_realk(buf,pos,dest,win)
    implicit none
    real(realk),intent(in) :: buf
    integer, intent(in) :: pos
    integer(kind=ls_mpik),intent(in) :: dest
    integer(kind=ls_mpik),intent(in) :: win
#ifdef VAR_MPI
    integer(kind=ls_mpik) :: n,ierr
    integer(kind=MPI_ADDRESS_KIND) :: offset
    ierr = 0
    n  = 1
    offset = int(pos-1,kind=MPI_ADDRESS_KIND)
    call MPI_GET(buf,n,MPI_DOUBLE_PRECISION,dest, &
     & offset,n,MPI_DOUBLE_PRECISION,win,ierr)
    if(ierr.ne.0)then
      call lsquit("Error(lsmpi_get_realk)",ierr)
    endif
#endif
  end subroutine lsmpi_get_realk
  subroutine lsmpi_get_realkV_wrapper8(buf,nelms,pos,dest,win)
    implicit none
    real(realk),intent(in) :: buf(*)
    integer, intent(in) :: pos
    integer(kind=8) :: nelms
    integer(kind=ls_mpik),intent(in) :: dest
    integer(kind=ls_mpik),intent(in) :: win
#ifdef VAR_MPI
    integer(kind=ls_mpik) :: n,ierr
    integer(kind=4) :: n4,k,i
    integer(kind=MPI_ADDRESS_KIND) :: offset
    if(ls_mpik==4)then
      k=SPLIT_MPI_MSG
      do i=1,nelms,k
        n4=k
        if(((nelms-i)<k).and.(mod(nelms-i+1,k)/=0))n4=mod(nelms,k)
        call lsmpi_get_realkV(buf(i:i+n4-1),n4,pos+i-1,dest,win)
      enddo
    else
      ierr = 0
      n  = nelms
      offset = int(pos-1,kind=MPI_ADDRESS_KIND)
      call MPI_GET(buf,n,MPI_DOUBLE_PRECISION,dest, &
       & offset,n,MPI_DOUBLE_PRECISION,win,ierr)
      if(ierr.ne.0)then
        call lsquit("Error(lsmpi_get_realk)",ierr)
      endif
    endif
#endif
  end subroutine lsmpi_get_realkV_wrapper8
  subroutine lsmpi_get_realkV(buf,nelms,pos,dest,win)
    implicit none
    real(realk),intent(in) :: buf(*)
    integer, intent(in) :: pos
    integer(kind=4) :: nelms
    integer(kind=ls_mpik),intent(in) :: dest
    integer(kind=ls_mpik),intent(in) :: win
#ifdef VAR_MPI
    integer(kind=ls_mpik) :: n,ierr
    integer(kind=MPI_ADDRESS_KIND) :: offset
    ierr = 0
    n  = nelms
    offset = int(pos-1,kind=MPI_ADDRESS_KIND)
    call MPI_GET(buf,n,MPI_DOUBLE_PRECISION,dest, &
     & offset,n,MPI_DOUBLE_PRECISION,win,ierr)
    if(ierr.ne.0)then
      call lsquit("Error(lsmpi_get_realk)",ierr)
    endif
#endif
  end subroutine lsmpi_get_realkV
  subroutine lsmpi_get_realkV_parts_wrapper8(buf,nelms,pos,dest,win,batchsze,flush_it)
    implicit none
    real(realk),intent(in) :: buf(*)
    integer, intent(in) :: pos
    integer(kind=8) :: nelms
    integer(kind=ls_mpik),intent(in) :: dest
    integer(kind=ls_mpik),intent(in) :: win
    integer, intent(in) :: batchsze
    logical, intent(in), optional :: flush_it
#ifdef VAR_MPI
    integer :: newpos
    integer(kind=4) :: n4,k,i
    integer(kind=8) :: n,j
    logical :: fi

    fi = .false.
    if(present(flush_it))fi = flush_it

    if(ls_mpik==4)then

      k=SPLIT_MPI_MSG

      do i=1,nelms,k
        n4=k
        if(((nelms-i)<k).and.(mod(nelms-i+1,k)/=0))n4=mod(nelms,k)

        call lsmpi_get_realkV_parts(buf(i:i+n4-1),n4,pos+i-1,dest,win,batchsze,flush_it=flush_it)
      enddo

    else

      do j=1,nelms,batchsze
        n=batchsze
        if(((nelms-j)<batchsze).and.&
          &(mod(nelms-j+1,batchsze)/=0))n=mod(nelms,batchsze)
        newpos = pos+j-1
        call lsmpi_get_realkV_wrapper8(buf(j:j+n-1),n,newpos,dest,win)
#ifdef VAR_HAVE_MPI3
      if(fi)call lsmpi_win_flush(win,rank=dest,local=.true.)
#endif
      enddo
    endif
#endif
  end subroutine lsmpi_get_realkV_parts_wrapper8
  subroutine lsmpi_get_realkV_parts(buf,nelms,pos,dest,win,batchsze,flush_it)
    implicit none
    real(realk),intent(in) :: buf(*)
    integer, intent(in) :: pos
    integer(kind=4),intent(in) :: nelms
    integer(kind=ls_mpik),intent(in) :: dest
    integer(kind=ls_mpik),intent(in) :: win
    integer, intent(in) :: batchsze
    logical, intent(in), optional :: flush_it
#ifdef VAR_MPI
    integer :: newpos
    integer(kind=8) :: n,i
    logical :: fi

    fi = .false.
    if(present(flush_it))fi = flush_it

    do i=1,nelms,batchsze

      n=batchsze

      if(((nelms-i)<batchsze).and.&
        &(mod(nelms-i+1,batchsze)/=0))n=mod(nelms,batchsze)

      newpos = pos+i-1

      call lsmpi_get_realkV_wrapper8(buf(i:i+n-1),n,newpos,dest,win)

#ifdef VAR_HAVE_MPI3
      if(fi)call lsmpi_win_flush(win,rank=dest,local=.true.)
#endif

    enddo
#endif
  end subroutine lsmpi_get_realkV_parts


  !=========================================================!
  !                   MPI GET ACC ROUTINES                  !
  !=========================================================!
  subroutine lsmpi_get_acc_int888(ibuf,obuf,dest,pos,win)
    implicit none
    integer(kind=8),intent(inout)    :: ibuf
    integer(kind=8),intent(inout)    :: obuf
    integer(kind=ls_mpik),intent(in) :: dest,win
    integer,intent(in)               :: pos
#ifdef VAR_MPI
    integer(kind=ls_mpik)            :: n,ierr
    integer(kind=MPI_ADDRESS_KIND)   :: offset
    
    n=1
    ierr = 0
    offset = int(pos-1,kind=MPI_ADDRESS_KIND)
#ifdef VAR_HAVE_MPI3
    call MPI_GET_ACCUMULATE(ibuf,n,MPI_INTEGER8,obuf,n,&
    &MPI_INTEGER8,dest,offset,n,MPI_INTEGER8,MPI_SUM,win,ierr)
#else
    call lsquit("ERROR(lsmpi_get_acc):you did not comile with an MPI3 enabled&
          &MPI library. Recompile.",-1)
#endif
#endif
  end subroutine lsmpi_get_acc_int888
  subroutine lsmpi_get_acc_int444(ibuf,obuf,dest,pos,win)
    implicit none
    integer(kind=4),intent(inout)    :: ibuf
    integer(kind=4),intent(inout)    :: obuf
    integer(kind=ls_mpik),intent(in) :: dest,win
    integer,intent(in)               :: pos
#ifdef VAR_MPI
    integer(kind=ls_mpik)            :: n,ierr
    integer(kind=MPI_ADDRESS_KIND)   :: offset
    
    n=1
    ierr = 0
    offset = int(pos-1,kind=MPI_ADDRESS_KIND)
#ifdef VAR_HAVE_MPI3
    call MPI_GET_ACCUMULATE(ibuf,n,MPI_INTEGER4,obuf,n,&
    &MPI_INTEGER4,dest,offset,n,MPI_INTEGER4,MPI_SUM,win,ierr)
#else
    call lsquit("ERROR(lsmpi_get_acc):you did not comile with an MPI3 enabled&
          &MPI library. Recompile.",-1)
#endif
#endif
  end subroutine lsmpi_get_acc_int444


  !=========================================================!
  !                   MPI ACC ROUTINES                      !
  !=========================================================!
  subroutine lsmpi_acc_int8(buf,pos,dest,win)
    implicit none
    integer(kind=8),intent(in) :: buf
    integer, intent(in) :: pos
    integer(kind=ls_mpik),intent(in) :: dest
    integer(kind=ls_mpik),intent(in) :: win
#ifdef VAR_MPI
    integer(kind=ls_mpik) :: n,ierr
    integer(kind=MPI_ADDRESS_KIND) :: offset
    ierr = 0
    n  = 1
    offset = int(pos-1,kind=MPI_ADDRESS_KIND)
    call MPI_ACCUMULATE(buf,n,MPI_INTEGER8,dest, &
     & offset,n,MPI_INTEGER8,MPI_SUM,win,ierr)
    if(ierr.ne.0)then
      call lsquit("Error(lsmpi_acc_realk)",ierr)
    endif
#endif
  end subroutine lsmpi_acc_int8
  subroutine lsmpi_acc_int4(buf,pos,dest,win)
    implicit none
    integer(kind=4),intent(in) :: buf
    integer, intent(in) :: pos
    integer(kind=ls_mpik),intent(in) :: dest
    integer(kind=ls_mpik),intent(in) :: win
#ifdef VAR_MPI
    integer(kind=ls_mpik) :: n,ierr
    integer(kind=MPI_ADDRESS_KIND) :: offset
    ierr = 0
    n  = 1
    offset = int(pos-1,kind=MPI_ADDRESS_KIND)
    call MPI_ACCUMULATE(buf,n,MPI_INTEGER4,dest, &
     & offset,n,MPI_INTEGER4,MPI_SUM,win,ierr)
    if(ierr.ne.0)then
      call lsquit("Error(lsmpi_acc_realk)",ierr)
    endif
#endif
  end subroutine lsmpi_acc_int4
  subroutine lsmpi_acc_realk(buf,pos,dest,win)
    implicit none
    real(realk),intent(in) :: buf
    integer, intent(in) :: pos
    integer(kind=ls_mpik),intent(in) :: dest
    integer(kind=ls_mpik),intent(in) :: win
#ifdef VAR_MPI
    integer(kind=ls_mpik) :: n,ierr
    integer(kind=MPI_ADDRESS_KIND) :: offset
    ierr = 0
    n  = 1
    offset = int(pos-1,kind=MPI_ADDRESS_KIND)
    call MPI_ACCUMULATE(buf,n,MPI_DOUBLE_PRECISION,dest, &
     & offset,n,MPI_DOUBLE_PRECISION,MPI_SUM,win,ierr)
    if(ierr.ne.0)then
      call lsquit("Error(lsmpi_acc_realk)",ierr)
    endif
#endif
  end subroutine lsmpi_acc_realk
  subroutine lsmpi_acc_realkV_wrapper8(buf,nelms,pos,dest,win)
    implicit none
    real(realk),intent(in) :: buf(*)
    integer, intent(in) :: pos
    integer(kind=8) :: nelms
    integer(kind=ls_mpik),intent(in) :: dest
    integer(kind=ls_mpik),intent(in) :: win
#ifdef VAR_MPI
    integer(kind=ls_mpik) :: n,ierr
    integer(kind=4) :: n4,k,i
    integer(kind=MPI_ADDRESS_KIND) :: offset
    if(ls_mpik==4)then
      k=SPLIT_MPI_MSG
      do i=1,nelms,k
        n4=k
        if(((nelms-i)<k).and.(mod(nelms-i+1,k)/=0))n4=mod(nelms,k)
        call lsmpi_acc_realkV(buf(i:i+n4-1),n4,pos+i-1,dest,win)
      enddo
    else
      ierr = 0
      n  = nelms
      offset = int(pos-1,kind=MPI_ADDRESS_KIND)
      call MPI_ACCUMULATE(buf,n,MPI_DOUBLE_PRECISION,dest, &
       & offset,n,MPI_DOUBLE_PRECISION,MPI_SUM,win,ierr)
      if(ierr.ne.0)then
        call lsquit("Error(lsmpi_acc_realk)",ierr)
      endif
    endif
#endif
  end subroutine lsmpi_acc_realkV_wrapper8
  subroutine lsmpi_acc_realkV(buf,nelms,pos,dest,win)
    implicit none
    real(realk),intent(in) :: buf(*)
    integer, intent(in) :: pos
    integer(kind=4) :: nelms
    integer(kind=ls_mpik),intent(in) :: dest
    integer(kind=ls_mpik),intent(in) :: win
#ifdef VAR_MPI
    integer(kind=ls_mpik) :: n,ierr
    integer(kind=MPI_ADDRESS_KIND) :: offset
    ierr = 0
    n  = nelms
    offset = int(pos-1,kind=MPI_ADDRESS_KIND)
    call MPI_ACCUMULATE(buf,n,MPI_DOUBLE_PRECISION,dest, &
     & offset,n,MPI_DOUBLE_PRECISION,MPI_SUM,win,ierr)
    if(ierr.ne.0)then
      call lsquit("Error(lsmpi_acc_realk)",ierr)
    endif
#endif
  end subroutine lsmpi_acc_realkV
  subroutine lsmpi_acc_realkV_parts_wrapper8(buf,nelms,pos,dest,win,batchsze,flush_it)
    implicit none
    real(realk),intent(in) :: buf(*)
    integer, intent(in) :: pos
    integer(kind=8) :: nelms
    integer(kind=ls_mpik),intent(in) :: dest
    integer(kind=ls_mpik),intent(in) :: win
    integer, intent(in) :: batchsze
    logical,optional, intent(in) :: flush_it
#ifdef VAR_MPI
    logical :: fi
    integer :: newpos
    integer(kind=4) :: n4,k,i
    integer(kind=8) :: n,j
    fi = .false.
    if(present(flush_it))fi = flush_it
    if(ls_mpik==4)then
      k=SPLIT_MPI_MSG
      do i=1,nelms,k
        n4=k
        if(((nelms-i)<k).and.(mod(nelms-i+1,k)/=0))n4=mod(nelms,k)
        call lsmpi_acc_realkV_parts(buf(i:i+n4-1),n4,pos+i-1,dest,&
        &win,batchsze,flush_it=flush_it)
      enddo
    else
      do j=1,nelms,batchsze
        n=batchsze
        if(((nelms-j)<batchsze).and.&
          &(mod(nelms-j+1,batchsze)/=0))n=mod(nelms,batchsze)
        newpos = pos+j-1
        call lsmpi_acc_realkV_wrapper8(buf(j:j+n-1),n,newpos,dest,win)
#ifdef VAR_HAVE_MPI3
        if(fi)call lsmpi_win_flush(win,rank=dest,local=.true.)
#endif
      enddo
    endif
#endif
  end subroutine lsmpi_acc_realkV_parts_wrapper8
  subroutine lsmpi_acc_realkV_parts(buf,nelms,pos,dest,win,batchsze,flush_it)
    implicit none
    real(realk),intent(in) :: buf(*)
    integer, intent(in) :: pos
    integer(kind=4),intent(in) :: nelms
    integer(kind=ls_mpik),intent(in) :: dest
    integer(kind=ls_mpik),intent(in) :: win
    integer, intent(in) :: batchsze
    logical,optional, intent(in) :: flush_it
#ifdef VAR_MPI
    logical :: fi
    integer :: newpos
    integer(kind=8) :: n,i

    fi = .false.
    if(present(flush_it))fi = flush_it

    do i=1,nelms,batchsze

      n=batchsze

      if(((nelms-i)<batchsze).and.&
        &(mod(nelms-i+1,batchsze)/=0))n=mod(nelms,batchsze)

      newpos = pos+i-1

      call lsmpi_acc_realkV_wrapper8(buf(i:i+n-1),n,newpos,dest,win)

#ifdef VAR_HAVE_MPI3
      if(fi)call lsmpi_win_flush(win,rank=dest,local=.true.)
#endif

    enddo
#endif
  end subroutine lsmpi_acc_realkV_parts
  

  subroutine lsmpi_localallgatherv_realk8(sendbuf,recbuf,reccounts,disps)
    implicit none
    real(realk), intent(in) :: sendbuf(:)
    real(realk), intent(inout) :: recbuf(:)
    integer(kind=8),intent(in) :: reccounts(:)
    integer(kind=8),intent(in) :: disps(:)
#ifdef VAR_MPI
    integer(kind=ls_mpik) :: ierr, dtype,n
    integer(kind=4) :: rc(infpar%lg_nodtot),dp(infpar%lg_nodtot)
    integer(kind=4) :: oldrc(infpar%lg_nodtot),i,k,n4,j
    integer :: node,nelms
    ierr = 0
    if(ls_mpik==4)then
      
      !get the maximum number of elements to loop over
      nelms=0
      do node=1,infpar%lg_nodtot
        nelms = max(nelms,reccounts(node))
      enddo


      k=SPLIT_MPI_MSG
      dp = disps
      oldrc = 0

      !loop over elements
      do i=1,nelms,k
    
        !get the displacements and number of elements to transfer 
        do node=0,infpar%lg_nodtot-1
          n4=k
          if(((reccounts(node+1)-i)<k).and.(mod(reccounts(node+1)-i+1,k)/=0))&
              &n4=mod(reccounts(node+1),k)
          if((reccounts(node+1)-i)<0)n4=0
          rc(node+1) = n4
          dp(node+1) = dp(node+1) + oldrc(node+1)
          oldrc(node+1) = rc(node+1)
        enddo
        
        !get the first element and number of elements on the current node
        if(i<reccounts(infpar%lg_mynum+1))then
          j=i
          n=rc(infpar%lg_mynum+1)
        else
          j=1
          n=0
        endif
        
        !do all gather
        call lsmpi_localallgatherv_realk4(sendbuf(j:j+n-1),recbuf,rc,dp)
      enddo
    else
      dtype = MPI_DOUBLE_PRECISION
      n = reccounts(infpar%lg_mynum+1)
     
#ifdef VAR_INT64
#ifdef VAR_MPI_32BIT_INT
  call lsquit('Error in lsmpi_localallgatherv_realk8 VAR_INT64 and VAR_MPI_32BIT_INT',-1)
#else
      call MPI_ALLGATHERV(sendbuf,n,dtype,recbuf,reccounts,&
                          &disps,dtype,infpar%lg_comm,ierr)
#endif
#else
  call lsquit('Error in lsmpi_localallgatherv_realk8 no VAR_INT64',-1)
#endif
     
      if(ierr/=0)then
        call lsquit("ERROR(lsmpi_localallgatherv_realk4):mpi is wrong",-1)
      endif
    endif
#endif
  end subroutine lsmpi_localallgatherv_realk8

  subroutine lsmpi_localallgatherv_realk4(sendbuf,recbuf,reccounts,disps)
    implicit none
    real(realk), intent(in) :: sendbuf(:)
    real(realk), intent(inout) :: recbuf(:)
    integer(kind=4),intent(in) :: reccounts(:)
    integer(kind=4),intent(in) :: disps(:)
#ifdef VAR_MPI
    integer(kind=ls_mpik) :: ierr, dtype,n
    ierr = 0
    dtype = MPI_DOUBLE_PRECISION
    n = reccounts(infpar%lg_mynum+1)

#ifdef VAR_INT64
#ifdef VAR_MPI_32BIT_INT
    call MPI_ALLGATHERV(sendbuf,n,dtype,recbuf,reccounts,&
                        &disps,dtype,infpar%lg_comm,ierr)
#else
    call lsquit('Error in lsmpi_localallgatherv_realk4 VAR_INT64 and not VAR_MPI_32BIT_INT',-1)
#endif
#else
    call MPI_ALLGATHERV(sendbuf,n,dtype,recbuf,reccounts,&
                        &disps,dtype,infpar%lg_comm,ierr)
#endif

    if(ierr/=0)then
      call lsquit("ERROR(lsmpi_localallgatherv_realk4):mpi is wrong",-1)
    endif
#endif
  end subroutine lsmpi_localallgatherv_realk4

  subroutine lsmpi_poke()
    implicit none
#ifdef VAR_MPI
    logical(kind=ls_mpik) :: flag
    integer(kind=ls_mpik) :: ierr
    real(realk) :: sta, sto
    ierr = 0
    if(.not.LSMPIASYNCP)then
      sta=MPI_WTIME()
      call mpi_iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,infpar%lg_comm,flag,status,ierr)
      sto=MPI_WTIME()
      poketime=poketime+sto-sta
      poketimes = poketimes + 1
    endif
#endif
  end subroutine lsmpi_poke

end module lsmpi_type


