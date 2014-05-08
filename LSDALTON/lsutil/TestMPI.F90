module lsmpi_test
  use precision
  use Integralparameters
  use memory_handling
  use lsmpi_type
#ifdef VAR_MPI
  use infpar_module
#endif

CONTAINS
  SUBROUTINE test_mpi(comm)
    implicit none
    integer(kind=ls_mpik) :: comm
#ifdef VAR_MPI_32BIT_INT
    integer :: single
    integer,pointer :: intbuffer(:)
    integer(kind=4) :: single4
    integer(kind=4),pointer :: intbuffer4(:)
    integer(kind=8) :: single8
    integer(kind=8),pointer :: intbuffer8(:)
    integer(kind=short) :: singleS
    integer(kind=short),pointer :: intbufferS(:)
    real(realk),pointer  :: DPbuffer(:)
    real(realk),parameter :: DM=101234567891234E0_realk
    real(realk),parameter :: DS=1.12345698765412E-16_realk
    real(realk),parameter :: D0=0.0E0_realk
    integer(kind=8),parameter :: S8R=1
    integer(kind=8),parameter :: S8S=9223372036854775806_8
    integer(kind=8),parameter :: S8M=-9223372036854775806_8
    integer(kind=8),parameter :: S80=0
    integer :: N, N3
    integer(kind=ls_mpik) :: i,master,mynum
    logical :: slave 
    integer(kind=4) ,parameter :: S4R=1
    integer(kind=4) ,parameter :: S4S=2147483647
    integer(kind=4) ,parameter :: S4M=-2147483647
    integer(kind=4) ,parameter :: S40=0
    integer(kind=short) ,parameter :: SSR=1
    integer(kind=short) ,parameter :: SSS=127
    integer(kind=short) ,parameter :: SSM=-128
    integer(kind=short) ,parameter :: SS0=0
    logical,pointer :: Lbuffer(:),Lbuffer2(:)
    character,pointer :: Cbuffer(:),Cbuffer2(:)
    character(len=19) :: Cbuffer3
    master = 0
    call get_rank_for_comm(comm,mynum)
    !wake up slaves    
    slave = mynum.NE.master
    IF(.NOT.slave)THEN
       !bcast integer
       call ls_mpibcast(LSMPITEST,master,comm)
    ENDIF

    !Test integer kind=8 bcast
    print*,'Test integer kind=8 bcast'
    N = 9
    N3 = 3
    CALL MEM_ALLOC(intbuffer8,N)
    IF(.NOT.slave)THEN
       do i=1,N3
          intbuffer8(i) = S8M+i
       enddo
       do i=N3+1,2*N3
          intbuffer8(i) = i
       enddo
       do i=2*N3+1,N
          intbuffer8(i) = S8S-i
       enddo
    ELSE
       do i=1,N
          intbuffer8(i) = S80
       enddo
    ENDIF
    call ls_mpibcast(intbuffer8,N,master,comm)
    IF(Slave)THEN
       do i=1,N3
          IF( intbuffer8(i) .NE. S8M+i)CALL LSQUIT('intBuffer error',-1)
       enddo
       do i=N3+1,2*N3
          IF( intbuffer8(i) .NE. i)CALL LSQUIT('intBuffer error',-1)
       enddo
       do i=2*N3+1,N
          IF( intbuffer8(i) .NE. S8S-i)CALL LSQUIT('intBuffer error',-1)
       enddo
    ENDIF
    CALL MEM_DEALLOC(intbuffer8)

    !Test integer kind=4 bcast
    print*,'Test integer kind=4 bcast'
    N = 9
    N3 = 3
    CALL MEM_ALLOC(intbuffer4,N)
    IF(.NOT.slave)THEN
       do i=1,N3
          intbuffer4(i) = S4M+i
       enddo
       do i=N3+1,2*N3
          intbuffer4(i) = i
       enddo
       do i=2*N3+1,N
          intbuffer4(i) = S4S-i
       enddo
    ELSE
       do i=1,N
          intbuffer4(i) = S40
       enddo
    ENDIF
    call ls_mpibcast(intbuffer4,N,master,comm)
    IF(Slave)THEN
       do i=1,N3
          IF( intbuffer4(i) .NE. S4M+i)CALL LSQUIT('intBuffer error',-1)
       enddo
       do i=N3+1,2*N3
          IF( intbuffer4(i) .NE. i)CALL LSQUIT('intBuffer error',-1)
       enddo
       do i=2*N3+1,N
          IF( intbuffer4(i) .NE. S4S-i)CALL LSQUIT('intBuffer error',-1)
       enddo
    ENDIF
    CALL MEM_DEALLOC(intbuffer4)

    !Test integer kind=short bcast
    print*,'Test integer kind=short bcast'
    N = 8
    N3 = 3
    CALL MEM_ALLOC(intbufferS,N)
    IF(.NOT.slave)THEN
       do i=1,N3
          intbufferS(i) = SSM+i
       enddo
       do i=N3+1,2*N3
          intbufferS(i) = i
       enddo
       do i=2*N3+1,N
          intbufferS(i) = SSS-i
       enddo
    ELSE
       do i=1,N
          intbufferS(i) = SS0
       enddo
    ENDIF
    call ls_mpibcast(intbufferS,N,master,comm)
    IF(Slave)THEN
       do i=1,N3
          IF( intbufferS(i) .NE. SSM+i)CALL LSQUIT('intBuffer error',-1)
       enddo
       do i=N3+1,2*N3
          IF( intbufferS(i) .NE. i)CALL LSQUIT('intBuffer error',-1)
       enddo
       do i=2*N3+1,N
          IF( intbufferS(i) .NE. SSS-i)CALL LSQUIT('intBuffer error',-1)
       enddo
    ENDIF
    CALL MEM_DEALLOC(intbufferS)

    !Test realk bcast
    print*,'Test realk bcast'
    N = 9
    N3 = 3
    CALL MEM_ALLOC(DPbuffer,N)
    IF(.NOT.slave)THEN
       do i=1,N3
          DPbuffer(i) = DM+i
       enddo
       do i=N3+1,2*N3
          DPbuffer(i) = i
       enddo
       do i=2*N3+1,N
          DPbuffer(i) = DS-i*DS
       enddo
    ELSE
       do i=1,N
          DPbuffer(i) = D0
       enddo
    ENDIF
    call ls_mpibcast(DPbuffer,N,master,comm)
    IF(Slave)THEN
       do i=1,N3
          IF( DPbuffer(i) .NE. DM+i)CALL LSQUIT('DPBuffer error',-1)
       enddo
       do i=N3+1,2*N3
          IF( DPbuffer(i) .NE. i)CALL LSQUIT('DPBuffer error',-1)
       enddo
       do i=2*N3+1,N
          IF( DPbuffer(i) .NE. DS-i*DS)CALL LSQUIT('DPBuffer error',-1)
       enddo
    ENDIF
    CALL MEM_DEALLOC(DPbuffer)

    !Test logical bcast
    print*,'Test logical bcast'
    N = 200
    CALL MEM_ALLOC(Lbuffer,N)
    IF(.NOT.slave)THEN
       do i=1,N
          Lbuffer(i) = .TRUE.
       enddo
       do i=1,N,17
          Lbuffer(i) = .FALSE.
       enddo
       do i=1,N,13
          Lbuffer(i) = .FALSE.
       enddo
       do i=1,N,11
          Lbuffer(i) = .FALSE.
       enddo
       do i=1,N,12
          Lbuffer(i) = .FALSE.
       enddo
    ELSE
       CALL MEM_ALLOC(Lbuffer2,N)
       do i=1,N
          Lbuffer2(i) = .TRUE.
       enddo
       do i=1,N,17
          Lbuffer2(i) = .FALSE.
       enddo
       do i=1,N,13
          Lbuffer2(i) = .FALSE.
       enddo
       do i=1,N,11
          Lbuffer2(i) = .FALSE.
       enddo
       do i=1,N,12
          Lbuffer2(i) = .FALSE.
       enddo
       do i=1,N
          Lbuffer(i) = .FALSE.
       enddo
    ENDIF
    call ls_mpibcast(Lbuffer,N,master,comm)
    IF(Slave)THEN
       do i=1,N
          IF( Lbuffer(i) .NEQV. Lbuffer2(i))CALL LSQUIT('LBuffer error',-1)
       enddo
       CALL MEM_DEALLOC(Lbuffer2)
    ENDIF
    CALL MEM_DEALLOC(Lbuffer)

    !Test charater bcast
    print*,'Test charater bcast'
    N = 20
    CALL MEM_ALLOC(Cbuffer,N)
    IF(.NOT.slave)THEN
       do i=1,N
          Cbuffer(i) = 'T'
       enddo
       do i=1,N,17
          Cbuffer(i) = 'A'
       enddo
       do i=1,N,13
          Cbuffer(i) = 'B'
       enddo
       do i=1,N,11
          Cbuffer(i) = 'C'
       enddo
       do i=1,N,12
          Cbuffer(i) = 'D'
       enddo
    ELSE
       CALL MEM_ALLOC(Cbuffer2,N)
       do i=1,N
          Cbuffer2(i) = 'T'
       enddo
       do i=1,N,17
          Cbuffer2(i) = 'A'
       enddo
       do i=1,N,13
          Cbuffer2(i) = 'B'
       enddo
       do i=1,N,11
          Cbuffer2(i) = 'C'
       enddo
       do i=1,N,12
          Cbuffer2(i) = 'D'
       enddo
       do i=1,N
          Cbuffer(i) = 'S'
       enddo
    ENDIF
    call ls_mpibcast(Cbuffer,N,master,comm)
    IF(Slave)THEN
       do i=1,N
          IF( Cbuffer(i) .NE. Cbuffer2(i))CALL LSQUIT('CBuffer error',-1)
       enddo
       CALL MEM_DEALLOC(Cbuffer2)
    ENDIF
    CALL MEM_DEALLOC(Cbuffer)
    N=19
    IF(slave)Cbuffer3 = 'Kasper is a MPI God'
    IF(.NOT.slave)Cbuffer3 = 'THOMAS is a MPI God'
    call ls_mpibcast(Cbuffer3,N,master,comm)
    IF(Slave)Then
       IF(Cbuffer3.NE.'THOMAS is a MPI God')CALL LSQUIT('CBuffer error2',-1)
    ENDIF
    call lsmpi_barrier(comm)
    print*,' START SEND RECV PART   mynum',mynum
    call lsmpi_barrier(comm)

    print*,'Test kind=8 send/recv  mynum',mynum
    !Test integer kind=8 bcast
    N = 9
    N3 = 3
    CALL MEM_ALLOC(intbuffer8,N)
    IF(.NOT.slave)THEN
       do i=1,N3
          intbuffer8(i) = S8M+i
       enddo
       do i=N3+1,2*N3
          intbuffer8(i) = i
       enddo
       do i=2*N3+1,N
          intbuffer8(i) = S8S-i
       enddo
    ELSE
       do i=1,N
          intbuffer8(i) = S80
       enddo
    ENDIF
    IF(Slave)THEN
       call ls_mpisendrecv(intbuffer8,N,comm,master,mynum)
    ELSE
       do i=1,infpar%nodtot-1
          call ls_mpisendrecv(intbuffer8,N,comm,mynum,i)
       enddo
    ENDIF
    IF(Slave)THEN
       do i=1,N3
          IF( intbuffer8(i) .NE. S8M+i)CALL LSQUIT('intBuffer error',-1)
       enddo
       do i=N3+1,2*N3
          IF( intbuffer8(i) .NE. i)CALL LSQUIT('intBuffer error',-1)
       enddo
       do i=2*N3+1,N
          IF( intbuffer8(i) .NE. S8S-i)CALL LSQUIT('intBuffer error',-1)
       enddo
    ENDIF
    CALL MEM_DEALLOC(intbuffer8)
    call lsmpi_barrier(comm)
    !Test integer kind=4 bcast
    print*,'Test kind=4 send/recv  mynum',mynum
    N = 9
    N3 = 3
    CALL MEM_ALLOC(intbuffer4,N)
    IF(.NOT.slave)THEN
       do i=1,N3
          intbuffer4(i) = S4M+i
       enddo
       do i=N3+1,2*N3
          intbuffer4(i) = i
       enddo
       do i=2*N3+1,N
          intbuffer4(i) = S4S-i
       enddo
    ELSE
       do i=1,N
          intbuffer4(i) = S40
       enddo
    ENDIF
    IF(Slave)THEN
       call ls_mpisendrecv(intbuffer4,N,comm,master,mynum)
    ELSE
       do i=1,infpar%nodtot-1
          call ls_mpisendrecv(intbuffer4,N,comm,mynum,i)
       enddo
    ENDIF
    IF(Slave)THEN
       do i=1,N3
          IF( intbuffer4(i) .NE. S4M+i)CALL LSQUIT('intBuffer error',-1)
       enddo
       do i=N3+1,2*N3
          IF( intbuffer4(i) .NE. i)CALL LSQUIT('intBuffer error',-1)
       enddo
       do i=2*N3+1,N
          IF( intbuffer4(i) .NE. S4S-i)CALL LSQUIT('intBuffer error',-1)
       enddo
    ENDIF
    CALL MEM_DEALLOC(intbuffer4)
    !Test integer kind=short bcast
    print*,'Test kind=short send/recv mynum',mynum
    N = 8
    N3 = 3
    CALL MEM_ALLOC(intbufferS,N)
    IF(.NOT.slave)THEN
       do i=1,N3
          intbufferS(i) = SSM+i
       enddo
       do i=N3+1,2*N3
          intbufferS(i) = i
       enddo
       do i=2*N3+1,N
          intbufferS(i) = SSS-i
       enddo
    ELSE
       do i=1,N
          intbufferS(i) = SS0
       enddo
    ENDIF
    IF(Slave)THEN
       call ls_mpisendrecv(intbufferS,N,comm,master,mynum)
    ELSE
       do i=1,infpar%nodtot-1
          call ls_mpisendrecv(intbufferS,N,comm,mynum,i)          
       enddo
    ENDIF
!    call ls_mpisendrecv(intbufferS,N,comm,master,mynum)
    IF(Slave)THEN
       do i=1,N3
          IF( intbufferS(i) .NE. SSM+i)CALL LSQUIT('intBuffer error',-1)
       enddo
       do i=N3+1,2*N3
          IF( intbufferS(i) .NE. i)CALL LSQUIT('intBuffer error',-1)
       enddo
       do i=2*N3+1,N
          IF( intbufferS(i) .NE. SSS-i)CALL LSQUIT('intBuffer error',-1)
       enddo
    ENDIF
    CALL MEM_DEALLOC(intbufferS)

    !Test realk bcast
    print*,'Test realk send/recv  mynum',mynum
    N = 9
    N3 = 3
    CALL MEM_ALLOC(DPbuffer,N)
    IF(.NOT.slave)THEN
       do i=1,N3
          DPbuffer(i) = DM+i
       enddo
       do i=N3+1,2*N3
          DPbuffer(i) = i
       enddo
       do i=2*N3+1,N
          DPbuffer(i) = DS-i*DS
       enddo
    ELSE
       do i=1,N
          DPbuffer(i) = D0
       enddo
    ENDIF
    IF(Slave)THEN
       call ls_mpisendrecv(DPbuffer,N,comm,master,mynum)
    ELSE
       do i=1,infpar%nodtot-1
          call ls_mpisendrecv(DPbuffer,N,comm,mynum,i)          
       enddo
    ENDIF
!    call ls_mpisendrecv(DPbuffer,N,comm,master,mynum)
!    call ls_mpibcast(DPbuffer,N,master,comm)
    IF(Slave)THEN
       do i=1,N3
          IF( DPbuffer(i) .NE. DM+i)CALL LSQUIT('DPBuffer error',-1)
       enddo
       do i=N3+1,2*N3
          IF( DPbuffer(i) .NE. i)CALL LSQUIT('DPBuffer error',-1)
       enddo
       do i=2*N3+1,N
          IF( DPbuffer(i) .NE. DS-i*DS)CALL LSQUIT('DPBuffer error',-1)
       enddo
    ENDIF
    CALL MEM_DEALLOC(DPbuffer)

    call lsmpi_barrier(comm)
    print*,'Test logical send/recv  mynum',mynum
    !Test logical bcast
    N = 200
    CALL MEM_ALLOC(Lbuffer,N)
    IF(.NOT.slave)THEN
       do i=1,N
          Lbuffer(i) = .TRUE.
       enddo
       do i=1,N,17
          Lbuffer(i) = .FALSE.
       enddo
       do i=1,N,13
          Lbuffer(i) = .FALSE.
       enddo
       do i=1,N,11
          Lbuffer(i) = .FALSE.
       enddo
       do i=1,N,12
          Lbuffer(i) = .FALSE.
       enddo
    ELSE
       CALL MEM_ALLOC(Lbuffer2,N)
       do i=1,N
          Lbuffer2(i) = .TRUE.
       enddo
       do i=1,N,17
          Lbuffer2(i) = .FALSE.
       enddo
       do i=1,N,13
          Lbuffer2(i) = .FALSE.
       enddo
       do i=1,N,11
          Lbuffer2(i) = .FALSE.
       enddo
       do i=1,N,12
          Lbuffer2(i) = .FALSE.
       enddo
       do i=1,N
          Lbuffer(i) = .FALSE.
       enddo
    ENDIF
    IF(Slave)THEN
       call ls_mpisendrecv(Lbuffer,N,comm,master,mynum)
    ELSE
       do i=1,infpar%nodtot-1
          call ls_mpisendrecv(Lbuffer,N,comm,mynum,i)          
       enddo
    ENDIF
    IF(Slave)THEN
       do i=1,N
          IF( Lbuffer(i) .NEQV. Lbuffer2(i))CALL LSQUIT('LBuffer error',-1)
       enddo
       CALL MEM_DEALLOC(Lbuffer2)
    ENDIF
    CALL MEM_DEALLOC(Lbuffer)

    call lsmpi_barrier(comm)
    !Test charater bcast
    print*,'Test character send/recv   mynum',mynum
    N = 20
    CALL MEM_ALLOC(Cbuffer,N)
    IF(.NOT.slave)THEN
       do i=1,N
          Cbuffer(i) = 'T'
       enddo
       do i=1,N,17
          Cbuffer(i) = 'A'
       enddo
       do i=1,N,13
          Cbuffer(i) = 'B'
       enddo
       do i=1,N,11
          Cbuffer(i) = 'C'
       enddo
       do i=1,N,12
          Cbuffer(i) = 'D'
       enddo
    ELSE
       CALL MEM_ALLOC(Cbuffer2,N)
       do i=1,N
          Cbuffer2(i) = 'T'
       enddo
       do i=1,N,17
          Cbuffer2(i) = 'A'
       enddo
       do i=1,N,13
          Cbuffer2(i) = 'B'
       enddo
       do i=1,N,11
          Cbuffer2(i) = 'C'
       enddo
       do i=1,N,12
          Cbuffer2(i) = 'D'
       enddo
       do i=1,N
          Cbuffer(i) = 'S'
       enddo
    ENDIF
    IF(Slave)THEN
       call ls_mpisendrecv(Cbuffer,N,comm,master,mynum)
    ELSE
       do i=1,infpar%nodtot-1
          call ls_mpisendrecv(Cbuffer,N,comm,mynum,i)          
       enddo
    ENDIF
!    call ls_mpisendrecv(Cbuffer,N,comm,master,mynum)
    IF(Slave)THEN
       do i=1,N
          IF( Cbuffer(i) .NE. Cbuffer2(i))CALL LSQUIT('CBuffer error',-1)
       enddo
       CALL MEM_DEALLOC(Cbuffer2)
    ENDIF
    CALL MEM_DEALLOC(Cbuffer)
    N=19
    IF(slave)Cbuffer3 = 'Kasper is a MPI God'
    IF(.NOT.slave)Cbuffer3 = 'THOMAS is a MPI God'
!    call ls_mpibcast(Cbuffer3,N,master,comm)
!    call ls_mpisendrecv(Cbuffer3,N,comm,master,mynum)
    IF(Slave)THEN
       call ls_mpisendrecv(Cbuffer3,N,comm,master,mynum)
    ELSE
       do i=1,infpar%nodtot-1
          call ls_mpisendrecv(Cbuffer3,N,comm,mynum,i)          
       enddo
    ENDIF
    IF(Slave)Then
       IF(Cbuffer3.NE.'THOMAS is a MPI God')CALL LSQUIT('CBuffer error2',-1)
    ENDIF
!!$    !Test integer(kind=8) buffer
!!$    N = 1243
!!$    CALL MEM_ALLOC(intbuffer8,N)
!!$    IF(.NOT.slave)THEN
!!$       do i=1,N
!!$          intbuffer8(i) = i*DR+(i-1)*DS
!!$       enddo
!!$    ELSE
!!$       intbuffer8(i) = D0
!!$    ENDIF
!!$    call ls_mpibcast(intbuffer8,N,master,comm)
!!$    IF(Slave)THEN
!!$       do i=1,N
!!$          IF( intbuffer8(i) .NE. i*DR+(i-1)*DS )CALL LSQUIT('intBuffer8 error',-1)
!!$       enddo
!!$    ENDIF

    !Test integer(kind=short) buffer


    call lsmpi_barrier(comm)

#endif

  end SUBROUTINE test_mpi

end module lsmpi_test
