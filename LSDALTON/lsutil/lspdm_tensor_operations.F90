!> @file
!> This is the file where all higher order pdm operations should go, especially,
!everything that has to do with arrays, allocating, feeing and so on, the one
!sided wrappers should go to lspdm_base_module
!> \author Patrick Ettenhuber
!> \date April 2013
module lspdm_tensor_operations_module


  ! Outside DEC directory
  use precision
  use ptr_assoc_module, only: ass_D1to4
  use dec_typedef_module
  use memory_handling
#ifdef VAR_MPI
  use infpar_module
  use lsmpi_type
#endif

  use tensor_basic_module
  use lspdm_basic_module


  !INTERFACES
  !**********
  interface array_get_tile
    module procedure array_gettile_combidx4,&
                    &array_gettile_combidx8,&
                    &array_gettile_modeidx
  end interface array_get_tile


  interface array_put_tile
    module procedure array_puttile_combidx4,&
                    &array_puttile_combidx8,&
                    &array_puttile_modeidx
  end interface array_put_tile

  interface array_accumulate_tile
    module procedure array_accumulate_tile_combidx4,&
                    &array_accumulate_tile_combidx8,&
                    &array_accumulate_tile_modeidx
  end interface array_accumulate_tile 

#ifdef COMPILER_UNDERSTANDS_FORTRAN_2003
  abstract interface
    subroutine put_acc_tile(arr,globtilenr,fort,nelms,lock_set)
      use precision
      import
      implicit none
      type(array),intent(in) :: arr
      integer,intent(in) :: globtilenr
#ifdef VAR_INT64
      integer(kind=8),intent(in) :: nelms
#else
      integer(kind=4),intent(in) :: nelms
#endif
      real(realk),intent(inout) :: fort(*)
      logical, optional, intent(in) :: lock_set
    end subroutine put_acc_tile
    subroutine put_acc_el(buf,pos,dest,win)
      use precision
      implicit none
      real(realk),intent(in) :: buf
      integer, intent(in) :: pos
      integer(kind=ls_mpik),intent(in) :: dest
      integer(kind=ls_mpik),intent(in) :: win
    end subroutine put_acc_el
    subroutine put_acc_vec(buf,nelms,pos,dest,win)
      use precision
      implicit none
      real(realk),intent(in) :: buf(*)
      integer, intent(in) :: pos
      integer(kind=8) :: nelms
      integer(kind=ls_mpik),intent(in) :: dest
      integer(kind=ls_mpik),intent(in) :: win
    end subroutine put_acc_vec
  end interface
#endif

  !interface array_accumulate_tile_nobuff
  !  module procedure array_accumulate_tile_combidx4_nobuff,&
  !                  &array_accumulate_tile_combidx8_nobuff,&
  !                  &array_accumulate_tile_modeidx_nobuff
  !end interface array_accumulate_tile_nobuff


  !Persistent array type definition
  !--------------------------------
  !> \brief the persistent array is a collection of n=500 arrays on each node
  !with some additional information. Here the storage fort tiled distributed and
  !replicated arrays is allocated, if 500  is not enough, please change here
  !> amount of arrays which are storable in the persistent array
  integer, parameter :: n_arrays = 500
  !> persistent array type-def
  type persistent_array
    !> collection of arrays
    type(array),pointer :: a(:)
    !> current address on node
    integer :: curr_addr_on_node  = 1
    !> counter for how many arrays were allocated in the persisten array
    integer :: arrays_allocated   = 0
    !> counter for how many arrays were deallocated
    integer :: arrays_deallocated = 0
    !> conter for the arrays currently in use
    integer :: arrays_in_use      = 0
    !> offset for the first tile allocation to get a better load distribution
    integer :: new_offset         = 0
    !> list of n logicals as indicator wheter an adress is free to allocate a
    !new array
    logical,pointer :: free_addr_on_node(:) => null()
  endtype persistent_array

  save
  
  ! job parameters for pdm jobs
  integer,parameter :: JOB_PC_ALLOC_DENSE      =  1
  integer,parameter :: JOB_PC_DEALLOC_DENSE    =  2
  integer,parameter :: JOB_FREE_ARR_STD        =  3
  integer,parameter :: JOB_INIT_ARR_TILED      =  4
  integer,parameter :: JOB_FREE_ARR_PDM        =  5
  integer,parameter :: JOB_INIT_ARR_REPLICATED =  6
  integer,parameter :: JOB_PRINT_MEM_INFO1     =  7
  integer,parameter :: JOB_PRINT_MEM_INFO2     =  8
  integer,parameter :: JOB_GET_NRM2_TILED      =  9
  integer,parameter :: JOB_DATA2TILED_DIST     = 10
  integer,parameter :: JOB_GET_TILE_SEND       = 11
  integer,parameter :: JOB_PRINT_TI_NRM        = 12
  integer,parameter :: JOB_SYNC_REPLICATED     = 13
  integer,parameter :: JOB_GET_NORM_REPLICATED = 14
  integer,parameter :: JOB_PREC_DOUBLES_PAR    = 15
  integer,parameter :: JOB_DDOT_PAR            = 16
  integer,parameter :: JOB_ADD_PAR             = 17
  integer,parameter :: JOB_CP_ARR              = 18
  integer,parameter :: JOB_ARRAY_ZERO          = 19
  integer,parameter :: JOB_GET_CC_ENERGY       = 20
  integer,parameter :: JOB_GET_FRAG_CC_ENERGY  = 21
  integer,parameter :: JOB_CHANGE_ACCESS_TYPE  = 22
  integer,parameter :: JOB_ARRAY_SCALE         = 23
  integer,parameter :: JOB_INIT_ARR_PC         = 24
  integer,parameter :: JOB_TEST_ARRAY          = 25

  !> definition of the persistent array 
  type(persistent_array) :: p_arr

  !> timing and measuring variables
  real(realk) :: time_pdm_acc          = 0.0E0_realk
  integer(kind=long) :: bytes_transferred_acc = 0
  integer(kind=long) :: nmsg_acc = 0
  real(realk) :: time_pdm_put          = 0.0E0_realk
  integer(kind=long) :: bytes_transferred_put = 0
  integer(kind=long) :: nmsg_put = 0
  real(realk) :: time_pdm_get          = 0.0E0_realk
  integer(kind=long) :: bytes_transferred_get = 0
  integer(kind=long) :: nmsg_get = 0

#ifdef VAR_MPI
#ifdef COMPILER_UNDERSTANDS_FORTRAN_2003
  procedure(array_acct4),pointer :: acc_ti4 
  procedure(array_acct8),pointer :: acc_ti8 
  procedure(array_gett4),pointer :: get_ti4 
  procedure(array_gett8),pointer :: get_ti8 
  procedure(array_putt4),pointer :: put_ti4 
  procedure(array_putt8),pointer :: put_ti8 
#endif
#endif

  !procedure(lsmpi_put_realkV_w8),pointer :: put_rk8 
  !procedure(lsmpi_get_realkV_w8),pointer :: get_rk8 
  !procedure(lsmpi_acc_realkV_w8),pointer :: acc_rk8 
  contains

  !>  \brief intitialize storage room for the tiled distributed arrays
  !> \author Patrick Ettenhuber
  !> \date May 2013
  subroutine init_persistent_array()
    implicit none
    call mem_alloc(p_arr%a,n_arrays)
    call mem_alloc(p_arr%free_addr_on_node,n_arrays)
    p_arr%free_addr_on_node=.true.
    !if( lspdm_use_comm_proc ) call lsquit("ERROR(init_persistent_array)&
    !& lspdm_use_comm_proc cannot be true at startup",-1)
    lspdm_use_comm_proc = .false.
  end subroutine init_persistent_array
  !>  \brief free storage room for the tiled distributed arrays
  !> \author Patrick Ettenhuber
  !> \date May 2013
  subroutine free_persistent_array()
    implicit none
    if(associated(p_arr%a))then
      call mem_dealloc(p_arr%a)
    endif
    if(associated(p_arr%free_addr_on_node))then
      call mem_dealloc(p_arr%free_addr_on_node)
    endif
    if( lspdm_use_comm_proc ) call lsquit("ERROR(free_persistent_array) &
    & lspdm_use_comm_proc has to be disabled at shutdown, otherwise there &
    & still might be processes running",-1)
  end subroutine free_persistent_array

  subroutine new_group_reset_persistent_array
    implicit none
    p_arr%new_offset = 0
  end subroutine new_group_reset_persistent_array
  


  subroutine lspdm_start_up_comm_procs
    implicit none
#ifdef VAR_MPI
    if(.not. lspdm_use_comm_proc)then

      if(infpar%lg_mynum == infpar%master .and. infpar%parent_comm == MPI_COMM_NULL)then
        write (*,'(55A)',advance='no')" STARTING UP THE COMMUNICATION PROCESSES (LSPDM) ..."
        !impregnate the slaves
        call ls_mpibcast(LSPDM_GIVE_BIRTH,infpar%master,infpar%lg_comm)
      endif

      lspdm_use_comm_proc = .true.
  
      if(infpar%parent_comm == MPI_COMM_NULL)then
        !all slaves and master get a baby where all the communicators are set
        call give_birth_to_child_process
        !call slaves to set lspdm_use_comm_proc to .true.
        call ls_mpibcast(LSPDM_GIVE_BIRTH,infpar%pc_mynum,infpar%pc_comm)
      endif

#ifdef VAR_LSDEBUG
      call lsmpi_barrier(infpar%pc_comm)
#endif

      if(infpar%parent_comm == MPI_COMM_NULL)then
#ifdef VAR_LSDEBUG
        call lsmpi_barrier(infpar%lg_comm)
#endif
        if(infpar%lg_mynum == infpar%master)then
          write(*,*) " SUCCESS"
        endif
      endif

    else
      print *,"WARNING(lspdm_startup_comm_procs): comm procs are already running"
    endif
#else
    call lsquit("ERROR(lspdm_start_up_comm_procs): not available without mpi",-1)
#endif
  end subroutine lspdm_start_up_comm_procs



  subroutine lspdm_shut_down_comm_procs
    implicit none
#ifdef VAR_MPI
    if( lspdm_use_comm_proc ) then

      if(infpar%lg_mynum == infpar%master .and. infpar%parent_comm == MPI_COMM_NULL )then

        print *,"SHUTTING DOWN THE COMMUNICATION PROCESSES (LSPDM)"
        !kill the babies of the slaves, i.e. get the slaves here
        call ls_mpibcast(LSPDM_SLAVES_SHUT_DOWN_CHILD,infpar%master,infpar%lg_comm)

      endif

      if( infpar%parent_comm == MPI_COMM_NULL )then
        ! slaves and master get their childs here
        call ls_mpibcast(LSPDM_SLAVES_SHUT_DOWN_CHILD,infpar%pc_mynum,infpar%pc_comm)
      endif

      !All master, slaves and children need to call the shut_down_child_process
      !routine to free communicators
      call shut_down_child_process
      lspdm_use_comm_proc = .false.
     

    else
      print *,"WARNING(lspdm_shutdown_comm_procs): no comm procs found running"
    endif
#else
    call lsquit("ERROR(lspdm_shut_down_comm_procs): not available without mpi",-1)
#endif
  end subroutine lspdm_shut_down_comm_procs


  
  !> \brief main subroutine for the communication of nodes on grid handling 
  !arr structures. this routine assumes, that always the process with rank 0 in
  !the given communicator is the "manager"
  !> \author Patrick Ettenhuber
  !> \date May 2012
  subroutine pdm_array_sync(comm,job,a,b,c,d,loc_addr)
    implicit none
    !> job is input for master and output for slaves, the arguments have to be
    !in the job paramenters list in top of this file
    integer                          :: job
    !> the communicator on which this routine should work
    integer(kind=ls_mpik),intent(in) :: comm
    !the array(s) to be passed to the slaves for which the operation is
    !performed
    type(array),optional             :: a,b,c,d
    logical,optional                 :: loc_addr

    !> comm arrays
    integer,pointer                  :: TMPI(:), dims(:)
    !character :: TMPC(12)
    integer :: i, j, context,modes(3),counter, stat,ierr,basic
    integer(kind=ls_mpik)            :: sendctr,root,me,nn
    logical                          :: loc
    modes=0
#ifdef VAR_MPI

    call get_rank_for_comm(comm,me)
    call get_size_for_comm(comm,nn)
    

    root  = infpar%master
    basic = 12

    loc = .false.
    if(present(loc_addr))loc = loc_addr

    IF( me == root) then
      !**************************************************************************************
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!code for MASTER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !**************************************************************************************
      !Wake up slaves
      call ls_mpibcast(PDMA4SLV, me, comm)
      !1     = JOB
      !2-5   = address in slot a-c
      !5-8   = modes a-c
      !9-13  = zero -> bool passed as integer for a-c
      !rest specifies dimensions
      counter = basic
      IF (PRESENT(A)) THEN
         counter     = counter+2*A%mode
      ENDIF
      IF (PRESENT(B)) THEN
         counter     = counter+2*B%mode
      ENDIF
      IF (PRESENT(C)) THEN
         counter     = counter+2*C%mode
      ENDIF
      IF (PRESENT(D)) THEN
         counter     = counter+2*D%mode
      ENDIF
      call ls_mpibcast(counter,root,comm)
      call mem_alloc(TMPI,counter)

      !change counter and basic for checking in the end
      TMPI(1) = counter
      counter = basic
      basic   = TMPI(1) 

      !get comm vector done
      TMPI    = 0
      TMPI(1) = job
      IF (PRESENT(A)) THEN
         TMPI(6)                        = A%mode
         if(A%zeros) TMPI(10)           = 1
         TMPI(counter+1:counter+A%mode) = A%dims
         counter = counter + A%mode
         TMPI(counter+1:counter+A%mode) = A%tdim
         counter = counter + A%mode
      ENDIF
      IF (PRESENT(B)) THEN
         TMPI(7)                        = B%mode
         if(B%zeros) TMPI(11)            = 1
         TMPI(counter+1:counter+B%mode) = B%dims
         counter = counter + B%mode
         TMPI(counter+1:counter+B%mode) = B%tdim
         counter = counter + B%mode
      ENDIF
      IF (PRESENT(C)) THEN
         TMPI(8)                        = C%mode
         if(C%zeros) TMPI(12)           = 1
         TMPI(counter+1:counter+C%mode) = C%dims
         counter = counter + C%mode
         TMPI(counter+1:counter+C%mode) = C%tdim
         counter = counter + C%mode
      ENDIF
      IF (PRESENT(D)) THEN
         TMPI(10)                       = D%mode
         if(D%zeros) TMPI(13)           = 1
         TMPI(counter+1:counter+D%mode) = D%dims
         counter = counter + D%mode
         TMPI(counter+1:counter+D%mode) = D%tdim
         counter = counter + D%mode
      ENDIF

      if(counter/=basic)call lsquit("ERROR(pdm_arr_sync):different number of&
      & elements for MASTER",DECinfo%output)
      if(loc)then
        if(nn>1.and.present(A).and..not.associated(A%addr_loc))&
        &call lsquit("ERROR(pdm_arr_sync):addr_loc for array A not associated",DECinfo%output)
        if(nn>1.and.present(B).and..not.associated(B%addr_loc))&
        &call lsquit("ERROR(pdm_arr_sync):addr_loc for array B not associated",DECinfo%output)
        if(nn>1.and.present(C).and..not.associated(C%addr_loc))&
        &call lsquit("ERROR(pdm_arr_sync):addr_loc for array C not associated",DECinfo%output)
        if(nn>1.and.present(D).and..not.associated(D%addr_loc))&
        &call lsquit("ERROR(pdm_arr_sync):addr_loc for array D not associated",DECinfo%output)
      else
        if(nn>1.and.present(A).and..not.associated(A%addr_p_arr))&
        &call lsquit("ERROR(pdm_arr_sync):addr_p_arr for array A not associated",DECinfo%output)
        if(nn>1.and.present(B).and..not.associated(B%addr_p_arr))&
        &call lsquit("ERROR(pdm_arr_sync):addr_p_arr for array B not associated",DECinfo%output)
        if(nn>1.and.present(C).and..not.associated(C%addr_p_arr))&
        &call lsquit("ERROR(pdm_arr_sync):addr_p_arr for array C not associated",DECinfo%output)
        if(nn>1.and.present(D).and..not.associated(D%addr_p_arr))&
        &call lsquit("ERROR(pdm_arr_sync):addr_p_arr for array D not associated",DECinfo%output)
      endif
      

      do sendctr=1,nn-1
        if(loc)then
          IF (PRESENT(A)) TMPI(2)  = A%addr_loc(sendctr+1)
          IF (PRESENT(B)) TMPI(3)  = B%addr_loc(sendctr+1)
          IF (PRESENT(C)) TMPI(4)  = C%addr_loc(sendctr+1)
          IF (PRESENT(D)) TMPI(5)  = D%addr_loc(sendctr+1)
        else
          IF (PRESENT(A)) TMPI(2)  = A%addr_p_arr(sendctr+1)
          IF (PRESENT(B)) TMPI(3)  = B%addr_p_arr(sendctr+1)
          IF (PRESENT(C)) TMPI(4)  = C%addr_p_arr(sendctr+1)
          IF (PRESENT(D)) TMPI(5)  = D%addr_p_arr(sendctr+1)
        endif
        call ls_mpisendrecv( TMPI, counter, comm, root, sendctr)
      enddo
      call mem_dealloc(TMPI)


    else  

      !**************************************************************************************
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!code for SLAVES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !**************************************************************************************
      call ls_mpibcast( counter, root, comm )
      call mem_alloc( TMPI, counter )
      call ls_mpisendrecv( TMPI, counter, comm, root, me)

      !get data from info vector; THIS COUNTER CONSTRUCTION HAS TO BE REWRITTEN
      !IF NEEDED FOR NOW IT IS CONVENIENT, BECAUSE IT IS SIMPLE
      counter = basic
      job = TMPI(1) !slaves needs to know what to do
      !1     = JOB
      !2-5   = address in slot a-c
      !6-9   = modes a-c
      !10-13 = zero -> bool passed as integer for a-c
      !rest specifies dimensions
      IF (TMPI(2).gt.0) THEN
         A = p_arr%a(TMPI(2))
      ELSE
         IF(TMPI(6).gt.0)THEN
           A%mode                = TMPI(6)
           if(TMPI(10)==1) A%zeros= .true.
           call arr_set_dims(A,TMPI(counter+1:counter+A%mode),A%mode)
           counter = counter + A%mode
           call arr_set_tdims(A,TMPI(counter+1:counter+A%mode),A%mode)
           counter = counter + A%mode
         ENDIF
      ENDIF
      IF (TMPI(3).gt.0) THEN
         B = p_arr%a(TMPI(3))
      ELSE
         IF(TMPI(7).gt.0)THEN
           B%mode                = TMPI(7)
           if(TMPI(11)==1)B%zeros= .true.
           call arr_set_dims(B,TMPI(counter+1:counter+B%mode),B%mode)
           counter = counter + B%mode
           call arr_set_tdims(B,TMPI(counter+1:counter+B%mode),B%mode)
           counter = counter + B%mode
         ENDIF
      ENDIF
      IF (TMPI(4).gt. 0) THEN
         C = p_arr%a(TMPI(4))
      ELSE
         IF(TMPI(8).gt.0)THEN
           C%mode                = TMPI(8)
           if(TMPI(12)==1)C%zeros= .true.
           call arr_set_dims(C,TMPI(counter+1:counter+C%mode),C%mode)
           counter = counter + C%mode
           call arr_set_tdims(C,TMPI(counter+1:counter+C%mode),C%mode)
           counter = counter + C%mode
         ENDIF
      ENDIF
      IF (TMPI(5).gt. 0) THEN
         !C = associate_to_p_arr(TMPI(4))
         D = p_arr%a(TMPI(5))
      ELSE
         IF(TMPI(9).gt.0)THEN
           D%mode                = TMPI(9)
           if(TMPI(13)==1)D%zeros= .true.
           call arr_set_dims(D,TMPI(counter+1:counter+D%mode),D%mode)
           counter = counter + D%mode
           call arr_set_tdims(D,TMPI(counter+1:counter+D%mode),D%mode)
           counter = counter + D%mode
         ENDIF
      ENDIF
      call mem_dealloc(TMPI)
    endif
#endif
  end subroutine pdm_array_sync

  !> \author Patrick Ettenhuber
  !> \date December 2012
  !> \brief get an array from the persistent array by specifying its address
  function get_arr_from_parr(addr) result(arr)
    implicit none
    !> the address of the array to extract
    integer,intent(in) :: addr
    !> array extracted from persisten array 
    type(array) :: arr
    arr=p_arr%a(addr)
  end function get_arr_from_parr

  function get_residence_of_tile(globaltilenumber,arr) result(rankofnode)
    implicit none
    type(array), intent(in) :: arr
    integer,intent(in) :: globaltilenumber
    integer :: rankofnode,nnod
    nnod=1
#ifdef VAR_MPI
    nnod=infpar%lg_nodtot
#endif      
    rankofnode=mod(globaltilenumber-1+arr%offset,nnod)
  end function get_residence_of_tile


  subroutine test_array(arr)
    type(array), intent(inout) :: arr
#ifdef VAR_MPI
    if(infpar%pc_mynum==0)then
      call pdm_array_sync(infpar%pc_comm,JOB_TEST_ARRAY,arr,loc_addr=.true.)
    endif
    print *,infpar%pc_mynum,"has da test",associated(arr%elm1),associated(arr%elm4),arr%addr_loc
#endif      
  end subroutine test_array


  !> \author Patrick Ettenhuber
  !> \date December 2012
  !> \brief calculate fragment eos cc energy in parallel (PDM)
  function get_fragment_cc_energy_parallel(t1,t2,gmo,occ_num,virt_num,occ_idx,virt_idx) result(fEc)
    implicit none
    !> singles amplitudes
    type(array), intent(inout) :: t1
    !> two electron integrals in the mo-basis
    type(array), intent(inout) :: gmo
    !> doubles amplitudes
    type(array), intent(in) :: t2
    !> number of occupied indices
    integer, intent(in) :: occ_num
    !> number of virtual indices
    integer, intent(in) :: virt_num
    !> referencing the occupied indices of the fragment to the full basis
    integer, intent(in) :: occ_idx(occ_num)
    !> referencing the virtueal indices of the fragment to the full basis
    integer, intent(in) :: virt_idx(virt_num)
    !> return-calue fEc contains the fragment correlation energy
    real(realk) :: Evirt,Eocc,fEc
    real(realk),pointer :: t(:,:,:,:)
    integer :: lt,i,j,a,b,o(t2%mode),fr_i,fr_j,fr_a,fr_b
    integer :: i_high,j_high,a_high,b_high

#ifdef VAR_MPI
    !Get the slaves to this routine
    if(infpar%lg_mynum==infpar%master)then
      call time_start_phase(PHASE_COMM)

      call pdm_array_sync(infpar%lg_comm,JOB_GET_FRAG_CC_ENERGY,t1,t2,gmo)
      call ls_mpiinitbuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)
      call ls_mpi_buffer(occ_num,infpar%master)
      call ls_mpi_buffer(occ_idx,occ_num,infpar%master)
      call ls_mpi_buffer(virt_num,infpar%master)
      call ls_mpi_buffer(virt_idx,virt_num,infpar%master)
      call ls_mpifinalizebuffer(infpar%master,LSMPIBROADCAST,infpar%lg_comm)

      call time_start_phase(PHASE_WORK)
    endif
    call memory_allocate_array_dense(gmo)

    call time_start_phase(PHASE_COMM)
    call cp_tileddata2fort(gmo,gmo%elm1,gmo%nelms,.true.)
    call time_start_phase(PHASE_WORK)

    Eocc  = 0.0E0_realk
    Evirt = 0.0E0_realk
    fEc   = 0.0E0_realk

    do lt=1,t2%nlti
      call ass_D1to4(t2%ti(lt)%t,t,t2%ti(lt)%d)
      !get offset for global indices
      call get_midx(t2%ti(lt)%gt,o,t2%ntpm,t2%mode)
      
      do j=1,t2%mode
        o(j)=(o(j)-1)*t2%tdim(j)
      enddo
    
      !find the limits of the current tile
      a_high= t2%ti(lt)%d(1)
      b_high= t2%ti(lt)%d(2)
      i_high= t2%ti(lt)%d(3)
      j_high= t2%ti(lt)%d(4)
      

      ! Energy using occupied scheme
      ! ****************************
      do j=1,occ_num
        fr_j = occ_idx(j)-o(4)   ! occupied EOS index in occupied AOS list
        if(j_high>=fr_j.and.fr_j>0)then 
          do i=1,occ_num
            fr_i = occ_idx(i)-o(3)
            if(i_high>=fr_i.and.fr_i>0)then
              do b=1,t2%ti(lt)%d(2)
                do a=1,t2%ti(lt)%d(1)
         
                  Eocc = Eocc + &
                   & ( t(a,b,fr_i,fr_j) + t1%elm2(a+o(1),fr_i+o(3))*t1%elm2(b+o(2),fr_j+o(4)) ) * &
                   & ( 2.0E0_realk*gmo%elm4(fr_i+o(3),a+o(1),fr_j+o(4),b+o(2)) &
                   & - gmo%elm4(fr_i+o(3),b+o(2),fr_j+o(4),a+o(1)) )
         
                end do
              end do
            endif
          end do
        endif
      end do
 
      ! Energy using virtual scheme
      do a=1,virt_num
        fr_a = virt_idx(a)-o(1)
        if(a_high>=fr_a.and.fr_a>0)then
          do j=1,t2%ti(lt)%d(4)
            do b=1,virt_num
              fr_b = virt_idx(b)-o(2)  ! virtual EOS index in occupied AOS list
              if(b_high>=fr_b.and.fr_b>0)then
                do i=1,t2%ti(lt)%d(3)

                  Evirt = Evirt + &
                    & ( t(fr_a,fr_b,i,j) + t1%elm2(fr_a+o(1),i+o(3))*t1%elm2(fr_b+o(2),j+o(4)) ) * &
                    & ( 2.0E0_realk*gmo%elm4(i+o(3),fr_a+o(1),j+o(4),fr_b+o(2)) &
                    &- gmo%elm4(i+o(3),fr_b+o(2),j+o(4),fr_a+o(1)) )

                end do
              endif
            end do
          end do
        endif
      end do

      nullify(t)

    enddo

      ! Hybrid scheme: Mixture of occupied and virtual partitioning schemes
      ! Ehybrid = 1/2* (Eocc + Evirt)

    call arr_deallocate_dense(gmo)
    
    call time_start_phase(PHASE_COMM)
    call lsmpi_local_reduction(Eocc,infpar%master)
    call lsmpi_local_reduction(Evirt,infpar%master)
    call time_start_phase(PHASE_WORK)

    fEc = 0.50E0_realk*(Eocc + Evirt)
#else
    fec = 0.0E0_realk
#endif
  end function get_fragment_cc_energy_parallel

  !> \author Patrick Ettenhuber
  !> \date December 2012
  !> \brief calculate aos cc energy in parallel (PDM)
  function get_cc_energy_parallel(t1,t2,gmo) result(Ec)
    implicit none
    !> singles amplitudes
    type(array), intent(inout) :: t1
    !> two electron integrals in the mo-basis
    type(array), intent(inout) :: gmo
    !> doubles amplitudes
    type(array), intent(in) :: t2
    !> on return Ec contains the correlation energy
    real(realk) :: E1,E2,Ec
    real(realk),pointer :: t(:,:,:,:)
    integer :: lt,i,j,a,b,o(t2%mode),da,db,di,dj

#ifdef VAR_MPI
    !Get the slaves to this routine
    if(infpar%lg_mynum==infpar%master)then
      call pdm_array_sync(infpar%lg_comm,JOB_GET_CC_ENERGY,t1,t2,gmo)
    endif
    call memory_allocate_array_dense(gmo)
    call cp_tileddata2fort(gmo,gmo%elm1,gmo%nelms,.true.)

    E1=0.0E0_realk
    E2=0.0E0_realk
    Ec=0.0E0_realk
    do lt=1,t2%nlti
      call ass_D1to4(t2%ti(lt)%t,t,t2%ti(lt)%d)
      !get offset for global indices
      call get_midx(t2%ti(lt)%gt,o,t2%ntpm,t2%mode)
      do j=1,t2%mode
        o(j)=(o(j)-1)*t2%tdim(j)
      enddo

      da = t2%ti(lt)%d(1)
      db = t2%ti(lt)%d(2)
      di = t2%ti(lt)%d(3)
      dj = t2%ti(lt)%d(4)
      !count over local indices
      !$OMP  PARALLEL DO DEFAULT(NONE) SHARED(gmo,o,t1,t,&
      !$OMP  da,db,di,dj) PRIVATE(i,j,a,b) REDUCTION(+:E1,E2) COLLAPSE(3)
      do j=1,dj
        do i=1,di
          do b=1,db
            do a=1,da
     
              E2 = E2 + t(a,b,i,j)*&
              & (2.0E0_realk*  gmo%elm4(i+o(3),a+o(1),j+o(4),b+o(2))-gmo%elm4(i+o(3),b+o(2),j+o(4),a+o(1)))
              E1 = E1 + ( t1%elm2(a+o(1),i+o(3))*t1%elm2(b+o(2),j+o(4)) ) * &
                   (2.0E0_realk*gmo%elm4(i+o(3),a+o(1),j+o(4),b+o(2))-gmo%elm4(i+o(3),b+o(2),j+o(4),a+o(1)))
   
            enddo 
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
      nullify(t)
    enddo

    call arr_deallocate_dense(gmo)
    
    call lsmpi_local_reduction(E1,infpar%master)
    call lsmpi_local_reduction(E2,infpar%master)

    Ec=E1+E2
#else
    Ec = 0.0E0_realk
#endif
  end function get_cc_energy_parallel

  !> \brief doubles preconditionning routine for pdm distributed doubles
  !amplitudes
  !> \author Patrick Ettenhuber
  !> \date december 2012
  subroutine precondition_doubles_parallel(omega2,ppfock,qqfock,prec)
    implicit none
    !> doubles residual, occupied and virtual blocks of the fock matrix
    type(array), intent(in) :: omega2,ppfock,qqfock
    !> output is the preconditioned doubles residual
    type(array), intent(inout) :: prec
    integer :: lt,a, b, i, j, dims(4)
    real(realk),pointer :: om(:,:,:,:),pp(:,:),qq(:,:),p(:,:,:,:)
    real(realk) :: nrm
    integer :: t(4),da,db,di,dj

#ifdef VAR_MPI

    !CHECK if the distributions are the same, if it becomes necessary, that they
    !are not, then this routine has to be rewritten
    if(omega2%tdim(1)/=prec%tdim(1).or.omega2%tdim(2)/=prec%tdim(2).or.&
    &omega2%tdim(3)/=prec%tdim(3).or.omega2%tdim(4)/=prec%tdim(4))then
      call lsquit("ERROR(precondition_doubles_parallel):omega2 and prec have&
      &different distributions", DECinfo%output)
    endif
    !Get the slaves to this routine
    if(infpar%lg_mynum==infpar%master)then
      call time_start_phase(PHASE_COMM)

      call pdm_array_sync(infpar%lg_comm,JOB_PREC_DOUBLES_PAR,omega2,ppfock,qqfock,prec)
      call time_start_phase(PHASE_WORK)
    endif

    dims=prec%dims

    
    !do a loop over the local tiles of the preconditioned matrix and get the
    !corresponding tiles of the residual to form the preconditioned residual
    do lt=1,prec%nlti

      call time_start_phase(PHASE_COMM)

      call array_get_tile(omega2,prec%ti(lt)%gt,prec%ti(lt)%t,prec%ti(lt)%e)

      call time_start_phase(PHASE_WORK)


      call ass_D1to4(prec%ti(lt)%t,om,prec%ti(lt)%d)
      
      !get offset for global indices
      call get_midx(prec%ti(lt)%gt,dims,prec%ntpm,prec%mode)
      do j=1,prec%mode
        dims(j)=(dims(j)-1)*prec%tdim(j)
      enddo

      !workaround for intel OMP error
      da = prec%ti(lt)%d(1)
      db = prec%ti(lt)%d(2)
      di = prec%ti(lt)%d(3)
      dj = prec%ti(lt)%d(4)

      !count over local indices
      !$OMP  PARALLEL DO DEFAULT(NONE)  SHARED(om,dims,&
      !$OMP  ppfock,qqfock,da,db,di,dj) PRIVATE(i,j,a,b) COLLAPSE(3)
      do j=1,dj
        do i=1,di
          do b=1,db
            do a=1,da
     
              om(a,b,i,j) = om(a,b,i,j) / &
                 ( ppfock%elm2(i+dims(3),i+dims(3)) - qqfock%elm2(a+dims(1),a+dims(1)) + &
                   ppfock%elm2(j+dims(4),j+dims(4)) - qqfock%elm2(b+dims(2),b+dims(2)) )
   
            enddo 
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
      nullify(om)
    enddo
    
    !crucial barrier, wait for all slaves to finish their jobs
    call time_start_phase(PHASE_IDLE)

    call lsmpi_barrier(infpar%lg_comm)

    call time_start_phase(PHASE_WORK)
#endif
  end subroutine precondition_doubles_parallel

  !> \brief calculate the dot product of two parallel distributed arrays. the
  !arrays must have the same tiling parameters, otherwise it is not implemented
  !> \author Patrick Ettenhuber
  !> \date december 2012
  function array_ddot_par(arr1,arr2,dest) result(res)
    implicit none
    !> the two arrays to calculate the dot-product from
    type(array),intent(in) :: arr1, arr2
    !> rank of the node to collect the result, -1 means all
    integer, intent(in) :: dest
    !> result
    real(realk) :: res
    real(realk),pointer :: buffer(:)
    integer :: lt,rem_els
    real(realk), external :: ddot
    integer(kind=ls_mpik) :: dest_mpi

#ifdef VAR_MPI
    !check if the init-types are the same
    if(arr1%access_type/=arr2%access_type)then
      call lsquit("ERROR(array_ddot_par):different init types of the&
      & arrays is not possible",DECinfo%output)
    endif

    !check if the destination to collet the resut makes sense in connection with
    !the access_type
    if(arr1%access_type==MASTER_ACCESS.and.dest/=0)then
      call lsquit("ERROR(array_ddot_par): the choice of destnation is&
      & useless",DECinfo%output)
    endif

    !get the slaves to this routine
    if(arr1%access_type==MASTER_ACCESS.and.infpar%lg_mynum==infpar%master)then
      call pdm_array_sync(infpar%lg_comm,JOB_DDOT_PAR,arr1,arr2)
    endif
    
    !zeroing the result
    res=0.0E0_realk
    
    !check for the same distribution of the arrays
    if(arr1%tdim(1)==arr2%tdim(1).and.arr1%tdim(2)==arr2%tdim(2).and.&
      &arr1%tdim(3)==arr2%tdim(3).and.arr1%tdim(4)==arr2%tdim(4))then
  
      !allocate buffer for the tiles
      call mem_alloc(buffer,arr1%tsize)
      buffer=0.0E0_realk
 
      !loop over local tiles of array2  and get the corresponding tiles of
      !array1
      do lt=1,arr2%nlti
        call array_get_tile(arr1,arr2%ti(lt)%gt,buffer,arr2%ti(lt)%e)
        res = res + ddot(arr2%ti(lt)%e,arr2%ti(lt)%t,1,buffer,1)
      enddo
      call mem_dealloc(buffer)
    else
      call lsquit("ERROR(array_ddot_par):NOT YET IMPLEMENTED, if the arrays have&
      & different distributions",DECinfo%output)
    endif

    !get result on the specified node/s
    if(dest==-1)then
      call lsmpi_allreduce(res,infpar%lg_comm)
    else
      dest_mpi=dest
      call lsmpi_local_reduction(res,dest_mpi)
    endif
#else
    res = 0.0E0_realk
#endif
  end function array_ddot_par

  !> \brief array addition routine for TILED_DIST arrays
  !> \author Patrick Ettenhuber
  !> \date January 2013
  subroutine array_add_par(x,b,y)
    implicit none
    !> array to collect the result in
    type(array), intent(inout) :: x
    !> array to add to x
    type(array), intent(in) :: y
    !> scale factor without intent, because it might be overwiritten for the slaves
    real(realk) :: b
    real(realk),pointer :: buffer(:)
    integer :: lt
#ifdef VAR_MPI

    !check if the access_types are the same
    if(x%access_type/=y%access_type)then
      call lsquit("ERROR(array_add_par):different init types&
      & impossible",DECinfo%output)
    endif

    !IF NOT MASTER_ACCESS all processes should know b on call-time, else b is
    !broadcasted here
    if(x%access_type==MASTER_ACCESS.and.infpar%lg_mynum==infpar%master)then
      call pdm_array_sync(infpar%lg_comm,JOB_ADD_PAR,x,y)
      call ls_mpibcast(b,infpar%master,infpar%lg_comm)
    else if(x%access_type==MASTER_ACCESS.and.infpar%lg_mynum/=infpar%master)then
      call ls_mpibcast(b,infpar%master,infpar%lg_comm)
    endif

    !check for the same distribution of the arrays
    if(x%tdim(1)==y%tdim(1).and.x%tdim(2)==y%tdim(2).and.&
      &x%tdim(3)==y%tdim(3).and.y%tdim(4)==y%tdim(4))then
      
      !allocate buffer for the tiles
      call mem_alloc(buffer,x%tsize)
  
      !lsoop over local tiles of array x
      do lt=1,x%nlti
        call array_get_tile(y,x%ti(lt)%gt,buffer,x%ti(lt)%e)
        call daxpy(x%ti(lt)%e,b,buffer,1,x%ti(lt)%t,1)
      enddo

      call mem_dealloc(buffer)
    else
      call lsquit("ERROR(array_add_par):NOT YET IMPLEMENTED, if the arrays have&
      & different distributions",DECinfo%output)
    endif

    !crucial barrier, because direct memory access is used
    call lsmpi_barrier(infpar%lg_comm)
#endif
  end subroutine array_add_par


  !> \brief array copying routine for TILED_DIST arrays
  !> \author Patrick Ettenhuber
  !> \date January 2013
  subroutine array_cp_tiled(from,to_ar)
    implicit none
    !> source, array to copy
    type(array), intent(in) :: from
    !> drain, the copied array
    type(array), intent(inout) :: to_ar
    real(realk),pointer :: buffer(:)
    integer :: lt
#ifdef VAR_MPI

    !check for the same access_types
    if(from%access_type/=to_ar%access_type)then
      call lsquit("ERROR(array_cp_tiled):different init types&
      & impossible",DECinfo%output)
    endif

    !get the slaves
    if(from%access_type==MASTER_ACCESS.and.infpar%lg_mynum==infpar%master)then
      call pdm_array_sync(infpar%lg_comm,JOB_CP_ARR,from,to_ar)
    endif

    !check for the same distributions
    if(from%tdim(1)==to_ar%tdim(1).and.from%tdim(2)==to_ar%tdim(2).and.&
      &from%tdim(3)==to_ar%tdim(3).and.to_ar%tdim(4)==to_ar%tdim(4))then
      do lt=1,to_ar%nlti
        call array_get_tile(from,to_ar%ti(lt)%gt,to_ar%ti(lt)%t,to_ar%ti(lt)%e)
      enddo
    else
      call lsquit("ERROR(array_cp_tiled):NOT YET IMPLEMENTED, if the arrato_ars have&
      & different distributions",DECinfo%output)
    endif

    !crucial barrier as remote direct memory access is used
    call lsmpi_barrier(infpar%lg_comm)
#endif
  end subroutine array_cp_tiled


  !> \brief zeroing routine for tiled distributed arrays
  !> \author Patrick Ettenhuber
  !> \date late 2012
  subroutine array_zero_tiled_dist(a)
    implicit none
    !> array to zero
    type(array),intent(inout) :: a
    integer :: lt
#ifdef VAR_MPI
    !get the slaves here
    if(a%access_type==MASTER_ACCESS.and.infpar%lg_mynum==infpar%master)then
      call pdm_array_sync(infpar%lg_comm,JOB_ARRAY_ZERO,a)
    endif

    !loop over local tiles and zero them individually
    do lt=1,a%nlti
      a%ti(lt)%t=0.0E0_realk
    enddo
#endif
  end subroutine array_zero_tiled_dist


  !> \author Patrick Ettenhuber
  !> \date January 2013
  !> \brief initialized a replicated matrix on each node
  function array_init_replicated(dims,nmodes,pdm)result(arr)
    implicit none
    !> array to be initialilzed
    type(array) :: arr
    !> number of modes and the dimensions of the array
    integer,intent(in) :: nmodes,dims(nmodes)
    !> integer specifying the access_type of the array
    integer,intent(in) :: pdm
    integer(kind=long) :: i,j
    integer :: addr,tdimdummy(nmodes)
    integer :: nelms
    integer(kind=ls_mpik) :: lg_nnodes,pc_nnodes
    integer(kind=ls_mpik) :: pc_me, lg_me
    integer, pointer :: lg_buf(:),pc_buf(:)
    logical :: master,pc_master,lg_master,child,parent

    !set the initial values and overwrite them later
    pc_nnodes               = 1
    pc_master               = .true.
    pc_me                   = 0
    lg_nnodes               = 1
    lg_master               = .true.
    child                   = .false.
    parent                  = .not.child

#ifdef VAR_MPI
    child     = (infpar%parent_comm /= MPI_COMM_NULL)
    parent    = .not.child

    !assign if master and the number of nodes in the local group
    if( lspdm_use_comm_proc ) then
      pc_me        = infpar%pc_mynum
      pc_nnodes    = infpar%pc_nodtot
      pc_master    = (infpar%parent_comm == MPI_COMM_NULL)
    endif

    if( parent )then
      lg_master = (infpar%lg_mynum==infpar%master)
      lg_nnodes = infpar%lg_nodtot
    endif
#endif

    master = (pc_master.and.lg_master)


    !allocate all pdm in p_arr therefore get free address and associate it with
    !the array, and increment the array counter
    p_arr%curr_addr_on_node = get_free_address(.true.)
    addr                    = p_arr%curr_addr_on_node
    p_arr%arrays_in_use     = p_arr%arrays_in_use + 1

    p_arr%a(addr)%access_type = pdm

    !SET MODE
    p_arr%a(addr)%mode = nmodes

    !SET DIMS
    call arr_set_dims(p_arr%a(addr),dims,nmodes)

    !SET ARRAY TYPE
    p_arr%a(addr)%itype=REPLICATED

    !SET NELMS
    nelms=1
    do i=1,nmodes
      nelms=nelms*dims(i)
    enddo
    p_arr%a(addr)%nelms=nelms

    !put 0 in tdim, since for the replicated array it is not important
    tdimdummy=0
    call arr_set_tdims(p_arr%a(addr),tdimdummy,nmodes)
    !SET NELMS

    !In the initialization the addess has to be set, since pdm_array_sync
    !depends on the  adresses, but setting them correctly is done later
    if( parent )then
      call mem_alloc(lg_buf,lg_nnodes)
      lg_buf = 0
    endif
    if( lspdm_use_comm_proc )then
      call mem_alloc(pc_buf,pc_nnodes)
      pc_buf = 0
    endif
    
    !if master init only master has to init the addresses addresses before
    !pdm syncronization
    if(lg_master .and. p_arr%a(addr)%access_type==MASTER_ACCESS.and.parent)then
      call arr_set_addr(p_arr%a(addr),lg_buf,lg_nnodes)
#ifdef VAR_MPI
      call pdm_array_sync(infpar%lg_comm,JOB_INIT_ARR_REPLICATED,p_arr%a(addr),loc_addr=.false.)
#endif
    endif

    if(pc_master .and.  p_arr%a(addr)%access_type==MASTER_ACCESS.and.lspdm_use_comm_proc)then
      call arr_set_addr(p_arr%a(addr),pc_buf,pc_nnodes,.true.)
#ifdef VAR_MPI
      call pdm_array_sync(infpar%pc_comm,JOB_INIT_ARR_REPLICATED,p_arr%a(addr),loc_addr=.true.)
#endif
    endif

    !if ALL_ACCESS all have to have the addresses allocated
    if(p_arr%a(addr)%access_type==ALL_ACCESS)then
      if(parent)call arr_set_addr(p_arr%a(addr),lg_buf,lg_nnodes)
      if(lspdm_use_comm_proc)call arr_set_addr(p_arr%a(addr),pc_buf,pc_nnodes,.true.)
    endif

#ifdef VAR_MPI
    !SET THE ADDRESSES ON ALL NODES     
    if( parent )then
      lg_buf(infpar%lg_mynum+1)=addr 
      call lsmpi_allreduce(lg_buf,lg_nnodes,infpar%lg_comm)
      call arr_set_addr(p_arr%a(addr),lg_buf,lg_nnodes,.false.)
    endif
    if( lspdm_use_comm_proc )then
      pc_buf(infpar%pc_mynum+1)=addr 
      call lsmpi_allreduce(pc_buf,pc_nnodes,infpar%pc_comm)
      call arr_set_addr(p_arr%a(addr),pc_buf,pc_nnodes,.true.)
    endif
#endif
  
    !ALLOCATE STORAGE SPACE FOR THE ARRAY
    call memory_allocate_array_dense(p_arr%a(addr))

    !RETURN THE CURRENLY ALLOCATE ARRAY
    arr=p_arr%a(addr)

    if(parent)call mem_dealloc(lg_buf)
    if(lspdm_use_comm_proc)call mem_dealloc(pc_buf)
  end function array_init_replicated


  !> \brief print the norm of a replicated array from each node, just a
  !debugging routine
  !> \author Patrick Ettenhuber
  !> \date January 2012
  function array_print_norm_repl(arr) result(nrm)
    implicit none
    !> replicated array to print the norm from
    type(array), intent(in) :: arr
    !return-value is the norm
    real(realk) :: nrm
    integer :: i
#ifdef VAR_MPI

    !get the slaves
    if(infpar%lg_mynum==infpar%master.and.arr%access_type==MASTER_ACCESS)then
      call pdm_array_sync(infpar%lg_comm,JOB_GET_NORM_REPLICATED,arr)
    endif

    !zero the norm an calculate it
    nrm =0.0E0_realk
    do i=1,arr%nelms
      nrm=nrm+arr%elm1(i)*arr%elm1(i)
    enddo
    print *,"on nodes",infpar%lg_mynum,sqrt(nrm)
#else
    nrm = 0.0E0_realk
#endif
  end function array_print_norm_repl


  !> \brief synchronize a replicated array from a source
  !> \author Patrick Ettenhuber
  !> \date cannot remember, 2012
  subroutine array_sync_replicated(arr,fromnode)
    implicit none
    !> array to synchronize
    type(array), intent(inout) :: arr
    !> specify the node which holds the original data that should be
    !synchronized to all nodes
    integer,optional, intent(in) :: fromnode
    integer(kind=ls_mpik) :: source
#ifdef VAR_MPI

    !give meaningful quit statement for useless input
    if(present(fromnode).and.arr%access_type==MASTER_ACCESS)then
      call lsquit("ERROR(array_sync_replicated): This combintion of input&
      &elements does not give sense",DECinfo%output)
      ! why would you want to collect the data on a node you cannot direcly
      ! access, or if you can access the data in the calling subroutine on the
      ! specified node, why is the init_tyep MASTER_ACCESS?
    endif

    ! get slaves
    if(infpar%lg_mynum==infpar%master.and.arr%access_type==MASTER_ACCESS)then
      call pdm_array_sync(infpar%lg_comm,JOB_SYNC_REPLICATED,arr)
    endif


    !specify the source of the data, by default master
    source = infpar%master
    if(present(fromnode))source=fromnode

    !do the synchronization
    call ls_mpibcast(arr%elm1,arr%nelms,source,infpar%lg_comm)
#endif    
  end subroutine array_sync_replicated


  !> \brief calculate the default tile-dimensions for the tiled dirtributed
  !array. default tile dimensions are currently a compromize between data
  !distribution and communication cost, it is cheaper to communicate large
  !chunks than many of them
  !> \author Patrick Ettenhuber
  !> \date march 2013
  subroutine array_default_batches(dims,nmodes,tdim,div)
    implicit none
    !> mode of the array
    integer :: nmodes
    !> dimensions in the modes
    integer :: dims(nmodes)
    !> divisor the last dimension whic is slict
    integer,intent(out) :: div
    !> tdim output 
    integer :: tdim(nmodes)
    integer :: i,j
    integer :: nlocalnodes
    integer :: cdims

    nlocalnodes=1
#ifdef VAR_MPI
    nlocalnodes=infpar%lg_nodtot
#endif    


    !calculate how many of the last modes have to be combined to get at least
    !the number of nodes tiles
    cdims=1
    do i=nmodes,1,-1
      if(cdims*dims(i)>nlocalnodes) exit
      cdims = cdims * dims(i)
    enddo

    !assing tiling dimensions
    do j=1,nmodes
      if(j<i)  tdim(j)=dims(j)
      if(j==i)then
        do div=1,dims(j)
          if(cdims*div>=nlocalnodes)exit
        enddo
        tdim(j)=(dims(j))/div
      endif
      if(j>i)  tdim(j)=1
    enddo

  end subroutine array_default_batches
  
  !> \brief calculate the number of tiles per mode
  !> \author Patrick Ettenhuber
  !> \date march 2013
  subroutine array_get_ntpm(dims,tdim,mode,ntpm,ntiles)
    implicit none
    !> number of modes and number of tiles
    integer :: mode,ntiles
    !> full dimensions, tile dimensinos, number of tiles per mode
    integer :: dims(mode),tdim(mode),ntpm(mode)
    integer :: i

    ntiles = 1

    do i=1,mode
      ntpm(i)= dims(i)/tdim(i)
      if(mod(dims(i),tdim(i))>0)then
        ntpm(i)=ntpm(i)+1
      endif
      ntiles = ntiles * ntpm(i)
    enddo

  end subroutine array_get_ntpm

  !> \author Patrick Ettenhuber
  !> \date September 2012
  !> \brief initialized a distributed tiled array
  function array_init_tiled(dims,nmodes,at,it,pdm,tdims,zeros_in_tiles,ps_d)result(arr)
    implicit none
    type(array) :: arr
    integer,intent(in) :: nmodes,dims(nmodes)
    character(4) :: at
    integer :: it, pdm
    integer,optional :: tdims(nmodes)
    logical, optional :: zeros_in_tiles
    logical, optional :: ps_d
    integer(kind=long) :: i,j
    integer ::addr,pdmt,k,div
    integer :: dflt(nmodes),cdims
    integer, pointer :: lg_buf(:),pc_buf(:)
    integer(kind=ls_mpik) :: lg_nnodes,pc_nnodes
    integer(kind=ls_mpik) :: pc_me, lg_me
    logical :: master,defdims, pseudo_dense
    logical :: pc_master,lg_master,child,parent
    integer :: infobuf(2)
   
    !set the initial values and overwrite them later
    pc_nnodes               = 1
    pc_master               = .true.
    pc_me                   = 0
    lg_nnodes               = 1
    lg_master               = .true.
    child                   = .false.
    parent                  = .not.child
    lg_me                   = 0
#ifdef VAR_MPI
    child     = (infpar%parent_comm /= MPI_COMM_NULL)
    parent    = .not.child

    !assign if master and the number of nodes in the local group
    if( lspdm_use_comm_proc ) then
      pc_me        = infpar%pc_mynum
      pc_nnodes    = infpar%pc_nodtot
      pc_master    = (infpar%parent_comm == MPI_COMM_NULL)
    endif

    if( parent )then
      lg_master = (infpar%lg_mynum==infpar%master)
      lg_nnodes = infpar%lg_nodtot
      lg_me     = infpar%lg_mynum
    endif
#endif
    pseudo_dense = .false.
    if(present(ps_d))pseudo_dense=ps_d

    master = (pc_master.and.lg_master)

    !allocate all tiled arrays in p_arr, get free
    p_arr%curr_addr_on_node = get_free_address(.true.)
    addr                    = p_arr%curr_addr_on_node
    p_arr%arrays_in_use     = p_arr%arrays_in_use + 1

    p_arr%a(addr)%access_type=pdm

    !INITIALIZE TILE STRUCTURE, if master from basics, if slave most is already
    !there
    defdims = .false.
    if(present(zeros_in_tiles)) p_arr%a(addr)%zeros = zeros_in_tiles

    !SET MODE
    p_arr%a(addr)%mode = nmodes

    !SET DIMS
    call arr_set_dims(p_arr%a(addr),dims,nmodes)

    if(present(tdims))then
      dflt=tdims
    else
      defdims=.true.
      !insert a routine for estimation of ts according to mem here
    endif

    !check if invalid numbers occur and fall back to default if so
    if(.not.defdims.and.present(tdims))then
      do i=1,nmodes
        if(dflt(i)<=0)then
          print *,"WARNING:INVALID NUMBER --> GET DEFAULT"
          defdims=.true.
          exit
        endif
      enddo
    endif

    !if needed, get default batch sizes, which are chosen such, that the best
    !distribution in terms of transfer speed and even distribution occur 
    !-> lots of consecutive !elements, big tiles, enough tiles
    if(defdims)then
      call array_default_batches(dims,nmodes,dflt,div)
    endif
    call arr_set_tdims(p_arr%a(addr),dflt,p_arr%a(addr)%mode)
    if (lg_master) then 
      p_arr%a(addr)%itype=it
      p_arr%a(addr)%atype=at
    end if    

    !divide A into tiles, according to dimensions
    !begin with counting the number of tiles needed in each mode
    dflt=0
    p_arr%a(addr)%nelms=1
    do i=1,p_arr%a(addr)%mode
      p_arr%a(addr)%nelms = p_arr%a(addr)%nelms * &
      &p_arr%a(addr)%dims(i)
      dflt(i)=p_arr%a(addr)%dims(i)/p_arr%a(addr)%tdim(i)
      if(mod(p_arr%a(addr)%dims(i),p_arr%a(addr)%tdim(i))>0)then
        dflt(i)=dflt(i)+1
      endif
    enddo
    call arr_set_ntpm(p_arr%a(addr),dflt,p_arr%a(addr)%mode)
    !print *,infpar%mynum,"ntpm:",arr%ntpm,arr%nelms
    !count the total number of tiles for the array and allocate in structure
    !calculate tilesize

    p_arr%a(addr)%ntiles = 1
    p_arr%a(addr)%tsize  = 1
    do i=1,p_arr%a(addr)%mode
      p_arr%a(addr)%ntiles = p_arr%a(addr)%ntiles * p_arr%a(addr)%ntpm(i)
      p_arr%a(addr)%tsize  = p_arr%a(addr)%tsize  * min(p_arr%a(addr)%tdim(i),p_arr%a(addr)%dims(i))
    enddo


    !In the initialization the addess has to be set, since pdm_array_sync
    !depends on the  adresses, but setting them correctly is done later
    if( parent )then
      call mem_alloc(lg_buf,2*lg_nnodes)
      lg_buf = 0
    endif
    if( lspdm_use_comm_proc )then
      call mem_alloc(pc_buf,pc_nnodes)
      pc_buf = 0
    endif
    
    !if master init only master has to get addresses
    if(lg_master .and. p_arr%a(addr)%access_type==MASTER_ACCESS.and.parent)then
      call arr_set_addr(p_arr%a(addr),lg_buf,lg_nnodes)
#ifdef VAR_MPI
      call pdm_array_sync(infpar%lg_comm,JOB_INIT_ARR_TILED,p_arr%a(addr))
#endif
    endif
    ! get child processes
    if(pc_master .and.  p_arr%a(addr)%access_type==MASTER_ACCESS.and.lspdm_use_comm_proc)then
      call arr_set_addr(p_arr%a(addr),pc_buf,pc_nnodes,.true.)
#ifdef VAR_MPI
      call pdm_array_sync(infpar%pc_comm,JOB_INIT_ARR_TILED,p_arr%a(addr),loc_addr=.true.)
#endif
    endif
#ifdef VAR_MPI
    if( lspdm_use_comm_proc ) then
       infobuf(1) = lg_me; infobuf(2) = 0; if(pseudo_dense) infobuf(2) = 1
       call ls_mpibcast(infobuf,2,infpar%master,infpar%pc_comm)
       lg_me = infobuf(1); pseudo_dense = (infobuf(2) == 1)
    endif
    call ls_mpibcast(p_arr%a(addr)%itype,infpar%master,infpar%lg_comm)
    call ls_mpibcast(p_arr%a(addr)%atype,4,infpar%master,infpar%lg_comm)
#endif

    !if ALL_ACCESS only all have to know the addresses
    if(p_arr%a(addr)%access_type==ALL_ACCESS)call arr_set_addr(p_arr%a(addr),lg_buf,lg_nnodes)

    call get_distribution_info(p_arr%a(addr))
#ifdef VAR_MPI
    if( parent )then
      lg_buf(infpar%lg_mynum+1)=addr 
      lg_buf(lg_nnodes+infpar%lg_mynum+1)=p_arr%a(addr)%offset
      call lsmpi_allreduce(lg_buf,2*lg_nnodes,infpar%lg_comm)
      call arr_set_addr(p_arr%a(addr),lg_buf,lg_nnodes)
      do i=1,lg_nnodes
        if(lg_buf(lg_nnodes+i)/=p_arr%a(addr)%offset)then
          print * ,infpar%lg_mynum,"found",lg_buf(lg_nnodes+i),p_arr%a(addr)%offset,i
          call lsquit("ERROR(array_init_tiled):offset &
          &is not the same on all nodes",DECinfo%output)
        endif
      enddo
    endif
    if( lspdm_use_comm_proc )then
      pc_buf(infpar%pc_mynum+1)=addr 
      call lsmpi_allreduce(pc_buf,pc_nnodes,infpar%pc_comm)
      call arr_set_addr(p_arr%a(addr),pc_buf,pc_nnodes,.true.)
    endif
#endif
 
    call arr_init_lock_set(p_arr%a(addr))
    call memory_allocate_tiles(p_arr%a(addr))

    if(pseudo_dense .and. lg_master)then
      call memory_allocate_array_dense(p_arr%a(addr))
    endif

    arr = p_arr%a(addr)
    !print *,infpar%lg_mynum,associated(arr%wi),"peristent",associated(p_arr%a(addr)%wi)

    if(parent)call mem_dealloc(lg_buf)
    if(lspdm_use_comm_proc)call mem_dealloc(pc_buf)
  end function array_init_tiled
  
  !> \brief add tiled distributed data to a basic fortran type array
  !> \author Patrick Ettenhuber
  !> date march 2013
  subroutine add_tileddata2fort(arr,b,fort,nelms,pdm,order)
    implicit none
    !> array to add to the input
    type(array),intent(in) :: arr
    !> basic fotran type array to which arr is added
    real(realk),intent(inout) :: fort(*)
    !> scaling factor for arr
    real(realk),intent(in) :: b
    !> nuber of elements in the array
    integer(kind=8), intent(in) :: nelms
    !> logical specifying whether the tiles are in pdm
    logical, intent(in) :: pdm
    !> reorder if the array is reorder with respect to the fortran array
    integer, intent(in), optional :: order(arr%mode)
    integer :: i,j,k,tmdidx(arr%mode), o(arr%mode)
    integer :: l,nelintile,tdim(arr%mode),fullfortdim(arr%mode)
    real(realk), pointer :: tmp(:)

    !check nelms
    if(nelms/=arr%nelms)call lsquit("ERROR(cp_tileddate2fort):array&
        &dimensions are not the same",DECinfo%output)

    !allocate space 
    if(pdm)then
      call mem_alloc(tmp,arr%tsize)
    endif

    do i=1,arr%mode
      o(i)=i
    enddo

    if(present(order))o=order

    do i = 1, arr%mode
      fullfortdim(i) = arr%dims(o(i))
    enddo

    do i=1,arr%ntiles
      call get_midx(i,tmdidx,arr%ntpm,arr%mode)
      call get_tile_dim(nelintile,i,arr%dims,arr%tdim,arr%mode)
      if(pdm)then
        call array_get_tile(arr,i,tmp,nelintile)
      else
        tmp => arr%ti(i)%t
      endif
      call tile_in_fort(b,tmp,i,arr%tdim,1.0E0_realk,fort,fullfortdim,arr%mode,o)
    enddo

    if(pdm)then
      call mem_dealloc(tmp)
    else
      nullify(tmp)
    endif
  end subroutine add_tileddata2fort

  subroutine cp_tileddata2fort(arr,fort,nelms,pdm,order)
    implicit none
    type(array),intent(in) :: arr
    real(realk),intent(inout) :: fort(*)
    integer(kind=8), intent(in) :: nelms
    logical, intent(in) :: pdm
    integer, intent(in), optional :: order(arr%mode)
    integer :: i,j,k,tmdidx(arr%mode),minimode(arr%mode),o(arr%mode)
    integer :: glbmodeidx(arr%mode),glbidx,l,nelintile,tdim(arr%mode),fullfortdim(arr%mode)
    real(realk), pointer :: tmp(:)


    if(nelms/=arr%nelms)call lsquit("ERROR(cp_tileddate2fort):array&
        &dimensions are not the same",DECinfo%output)
    if(pdm)then
      call mem_alloc(tmp,arr%tsize)
    endif
    do i=1,arr%mode
      o(i)=i
    enddo
    if(present(order))o=order

    do i = 1, arr%mode
      fullfortdim(i) = arr%dims(o(i))
    enddo

    do i=1,arr%ntiles
      call get_midx(i,tmdidx,arr%ntpm,arr%mode)
      call get_tile_dim(l,i,arr%dims,arr%tdim,arr%mode,2)
      call get_tile_dim(nelintile,i,arr%dims,arr%tdim,arr%mode)
      if(pdm)then
        call array_get_tile(arr,i,tmp,nelintile)
      else
        tmp => arr%ti(i)%t
      endif
      call tile_in_fort(1.0E0_realk,tmp,i,arr%tdim,0.0E0_realk,fort,fullfortdim,arr%mode,o)
    enddo

    if(pdm)then
      call mem_dealloc(tmp)
    else
      nullify(tmp)
    endif
  end subroutine cp_tileddata2fort



  subroutine array_scatteradd_densetotiled(arr,sc,A,nelms,nod,optorder)
    implicit none
    type(array),intent(inout) :: arr
    real(realk),intent(in) :: A(*)
    real(realk),intent(in) :: sc
    integer(kind=long),intent(in) :: nelms
    integer(kind=ls_mpik),intent(in) :: nod
    integer,intent(in),optional :: optorder(arr%mode)
    real(realk),pointer :: buf(:)
    integer :: nelmsit,i, order(arr%mode)
    integer :: ltidx,fullfortdims(arr%mode)
    integer(kind=ls_mpik) :: nnod,dest,me
#ifdef VAR_MPI

    !TRY TO INCLUDE MPI_PUT FOR THAT OPERATION,SO MAYBE MASTER_SLAVE DEPENDENCE
    !IN THIS ROUTINE IS GONE
    do i=1,arr%mode
      order(i)=i
    enddo
    if(present(optorder))order=optorder
    do i=1,arr%mode
      fullfortdims(order(i)) = arr%dims(i)
    enddo
   

    me = 0
    nnod=1
    me=infpar%lg_mynum
    nnod=infpar%lg_nodtot
    !begin with sanity checks
    ! corresponding elements
    call mem_alloc(buf,arr%tsize)

    
    do i=1,arr%ntiles
      dest=get_residence_of_tile(i,arr)
      call get_tile_dim(nelmsit,arr,i)
      if(dest==me.and.nod==me)then
        ltidx = (i - 1) /nnod + 1
        call tile_from_fort(1.0E0_realk,A,fullfortdims,arr%mode,0.0E0_realk,buf,i,arr%tdim,order)
        call daxpy(nelmsit,sc,buf,1,arr%ti(ltidx)%t,1)
      else if(nod==me)then
        call tile_from_fort(1.0E0_realk,A,fullfortdims,arr%mode,0.0E0_realk,buf,i,arr%tdim,order)
        call lsmpi_send(buf,nelmsit,infpar%lg_comm,dest)
      else if(dest==me)then
        ltidx = (i - 1) /nnod + 1
        call lsmpi_recv(buf,nelmsit,infpar%lg_comm,nod)
        call daxpy(nelmsit,sc,buf,1,arr%ti(ltidx)%t,1)
      endif
    enddo
    call mem_dealloc(buf)
#else
    call lsquit("ERROR(array_scatteradd_densetotiled):this routine is MPI only",-1)
#endif
  end subroutine array_scatteradd_densetotiled



  ! arr = pre1 * fort + pre2 * arr
  subroutine array_scatter(pre1,fort,pre2,arr,nelms,oo,wrk,iwrk)
    implicit none
    real(realk),intent(in)             :: pre1,pre2
    type(array),intent(in)             :: arr
    real(realk),intent(inout)          :: fort(*)
    integer(kind=long), intent(in)     :: nelms
    integer(kind=ls_mpik)              :: nod
    integer, intent(in), optional             :: oo(arr%mode)
    real(realk),intent(inout),target,optional :: wrk(*)
    integer(kind=8),intent(in),optional,target:: iwrk
    integer(kind=ls_mpik) :: src,me,nnod
    integer               :: i,ltidx,o(arr%mode)
    integer               :: nelintile,fullfortdim(arr%mode)
    real(realk), pointer  :: tmp(:)
    integer               :: tmps 
    logical               :: internal_alloc,lock_outside
    integer               :: maxintmp,b,e,minstart
#ifdef VAR_MPI
#ifdef COMPILER_UNDERSTANDS_FORTRAN_2003
    procedure(put_acc_tile), pointer :: put_acc => null()

    acc_ti8 => array_acct8
    acc_ti4 => array_acct4
    put_ti8 => array_putt8
    put_ti4 => array_putt4
  
    do i=1,arr%mode
      o(i)=i
    enddo
    if(present(oo))o=oo

#ifdef VAR_INT64
    if(pre2==0.0E0_realk) put_acc => put_ti8
    if(pre2/=0.0E0_realk) put_acc => acc_ti8
#else
    if(pre2==0.0E0_realk) put_acc => put_ti4
    if(pre2/=0.0E0_realk) put_acc => acc_ti4
#endif
   
    if(pre2/=0.0E0_realk.and.pre2/=1.0E0_realk)then
      call array_scale_td(arr,pre2)
    endif

#ifdef VAR_LSDEBUG
    if((present(wrk).and..not.present(iwrk)).or.(.not.present(wrk).and.present(iwrk)))then
      call lsquit("ERROR(array_scatter):both or neither wrk and iwrk have to &
                  &be given",-1)
    endif
#endif

    internal_alloc = .true.
    if(present(wrk).and.present(iwrk))then
      if(iwrk>arr%tsize)then
        internal_alloc=.false.
#ifdef VAR_LSDEBUG
      else
        print *,"WARNING(array_scatter):allocating internally, given buffer not large enough"
#endif
      endif
    endif

    me=0
    nnod=1
    me=infpar%lg_mynum
    nnod=infpar%lg_nodtot

#ifdef VAR_LSDEBUG
    if(nelms/=arr%nelms)call lsquit("ERROR(array_scatter):array&
        &dimensions are not the same",DECinfo%output)
#endif

    do i = 1, arr%mode
      fullfortdim(i) = arr%dims(o(i))
    enddo

    if(internal_alloc)then
      tmps = arr%tsize
      call mem_alloc(tmp,tmps)
    else
      tmps =  iwrk
      tmp  => wrk(1:tmps)
    endif
  
    maxintmp = tmps / arr%tsize

    do i=1,arr%ntiles
      if(i>maxintmp)then
       if(arr%lock_set(i-maxintmp)) call arr_unlock_win(arr,i-maxintmp)
      endif
      b = 1 + mod(i - 1, maxintmp) * arr%tsize
      e = b + arr%tsize -1
      call tile_from_fort(pre1,fort,fullfortdim,arr%mode,&
               &0.0E0_realk,tmp(b:e),i,arr%tdim,o)
      call get_tile_dim(nelintile,i,arr%dims,arr%tdim,arr%mode)
      call put_acc(arr,i,tmp(b:e),nelintile,arr%lock_set(i))
    enddo

    if(arr%ntiles - maxintmp >= 0)then
      minstart = arr%ntiles - maxintmp + 1
    else
      minstart = 1
    endif

    do i=minstart, arr%ntiles
      if(arr%lock_set(i))call arr_unlock_win(arr,i)
    enddo

    if(internal_alloc)then
      call mem_dealloc(tmp)
    else
      tmp  => null()
    endif
#else
    call lsquit("ERROR(array_scatter):this routine is FORTRAN 2003 only",-1)
#endif
#else
    call lsquit("ERROR(array_scatter):this routine is MPI only",-1)
#endif
  end subroutine array_scatter


  ! fort = pre1 * arr + pre2 *fort
  subroutine array_gather(pre1,arr,pre2,fort,nelms,oo,wrk,iwrk)
    implicit none
    real(realk),intent(in)             :: pre1,pre2
    type(array),intent(in)             :: arr
    real(realk),intent(inout)          :: fort(*)
    integer(kind=long), intent(in)     :: nelms
    integer(kind=ls_mpik)              :: nod
    integer, intent(in), optional             :: oo(arr%mode)
    real(realk),intent(inout),target,optional :: wrk(*)
    integer(kind=8),intent(in),optional,target:: iwrk
    integer(kind=ls_mpik) :: src,me,nnod
    integer               :: i,ltidx,o(arr%mode)
    integer               :: nelintile,fullfortdim(arr%mode)
    real(realk), pointer  :: tmp(:)
    integer               :: tmps 
    logical               :: internal_alloc,lock_outside,so
    integer               :: maxintmp,b,e,minstart
#ifdef VAR_MPI
  
    do i=1,arr%mode
      o(i)=i
    enddo
    if(present(oo))o=oo
 
    so = .true.
    do i=1,arr%mode
      if(o(i)/=i)so = .false.
    enddo

#ifdef VAR_LSDEBUG
    if((present(wrk).and..not.present(iwrk)).or.(.not.present(wrk).and.present(iwrk)))then
      call lsquit('ERROR(array_gather):both or neither wrk and iwrk have to &
                  &be given',-1)
    endif
#endif


    !CHECK IF INTERNAL MEMORY ALLOCATION IS NEEDED
    internal_alloc = .true.
    if(present(wrk).and.present(iwrk))then
      if(iwrk>arr%tsize)then
        internal_alloc=.false.
#ifdef VAR_LSDEBUG
      else
        print *,'WARNING(array_gather):allocating internally, given buffer not large enough'
#endif
      endif
    endif

    me=0
    nnod=1
    me=infpar%lg_mynum
    nnod=infpar%lg_nodtot

#ifdef VAR_LSDEBUG
    if(nelms/=arr%nelms)call lsquit('ERROR(array_gather):array&
        &dimensions are not the same',DECinfo%output)
#endif

    do i = 1, arr%mode
      fullfortdim(i) = arr%dims(o(i))
    enddo

  

    if(so.and.pre1==1.0E0_realk.and.pre2==0.0E0_realk)then
      b=1
      do i=1,arr%ntiles
        call get_tile_dim(nelintile,i,arr%dims,arr%tdim,arr%mode)
        e = b + nelintile - 1
        call array_get_tile(arr,i,fort(b:e),nelintile,arr%lock_set(i))
        b = e + 1
      enddo
    else

      if(internal_alloc)then
#ifdef VAR_LSDEBUG
        print *,'WARINING(array_gather):Allocating internally'
#endif
        tmps = arr%tsize
        call mem_alloc(tmp,tmps)
      else
        tmps =  iwrk
        tmp  => wrk(1:tmps)
      endif

      maxintmp = tmps / arr%tsize

      do i=1,arr%ntiles
        if(i>maxintmp)then
          b = 1 + mod(i - maxintmp - 1, maxintmp) * arr%tsize
          e = b + arr%tsize -1
          if(arr%lock_set(i-maxintmp))call arr_unlock_win(arr,i-maxintmp)
          call tile_in_fort(pre1,tmp(b:e),i-maxintmp,arr%tdim,&
                 &pre2,fort,fullfortdim,arr%mode,o)
        endif
        b = 1 + mod(i - 1, maxintmp) * arr%tsize
        e = b + arr%tsize -1
        call get_tile_dim(nelintile,i,arr%dims,arr%tdim,arr%mode)
        call array_get_tile(arr,i,tmp(b:e),nelintile,arr%lock_set(i))
      enddo
     
      if(arr%ntiles - maxintmp >= 0)then
        minstart = arr%ntiles - maxintmp + 1
      else
        minstart = 1
      endif
     
      do i=minstart, arr%ntiles
        b = 1 + mod(i - 1, maxintmp) * arr%tsize
        e = b + arr%tsize -1
        if(arr%lock_set(i))call arr_unlock_win(arr,i)
        call tile_in_fort(pre1,tmp(b:e),i,arr%tdim,&
               &pre2,fort,fullfortdim,arr%mode,o)
      enddo

      if(internal_alloc)then
        call mem_dealloc(tmp)
      else
        tmp  => null()
      endif
    endif

#else
    call lsquit('ERROR(array_gather):this routine is MPI only',-1)
#endif
  end subroutine array_gather

  subroutine array_two_dim_1batch(arr,o,op,fort,n2comb,fel,tl,lock_outside,debug)
    implicit none
    type(array),intent(in) :: arr
    real(realk),intent(inout) :: fort(*)
    character, intent(in) :: op
    integer, intent(in),target :: o(arr%mode)
    integer, intent(in) :: fel,tl,n2comb
    logical, intent(in) :: lock_outside
    logical, intent(in),optional :: debug
    integer :: fordims(arr%mode)
    integer,target :: fx(arr%mode)
    integer,target :: flx(arr%mode)
    integer :: oldidx(arr%mode)
    integer :: i,lel
    integer,target :: ro(arr%mode)
    integer :: comb1,comb2,c1,c2
    integer :: tidx(arr%mode),idxt(arr%mode)
    integer :: tlidx(arr%mode),lidxt(arr%mode)
    integer :: ctidx, cidxt, cidxf, st_tiling
    integer,pointer :: u_o(:),u_ro(:),tinfo(:,:)
    integer,pointer ::for3,for4
    real(realk), pointer :: p_fort3(:,:,:),p_fort2(:,:)
    integer :: tsze(arr%mode),mult1,mult2,dummy
    integer(kind=8) :: cons_el_in_t,cons_els,tl_max,tl_mod
    integer(kind=8) :: cons_el_rd
    integer(kind=8) :: part1,part2,split_in, diff_ord,modp1,modp2
    logical :: deb
#ifdef VAR_MPI
#ifdef COMPILER_UNDERSTANDS_FORTRAN_2003
    procedure(put_acc_el), pointer :: pga => null()
    procedure(put_acc_vec), pointer :: pgav => null()
    integer(kind=ls_mpik) :: source

    deb = .false.
    if(present(debug))deb = debug

    if( op/='a'.and.op/='p'.and.op/='g')then
      call lsquit("ERROR(array_two_dim_1batch):unknown choice of operator",-1)
    endif 
    if(op=='p')then
      pga  => lsmpi_put_realk
      pgav => lsmpi_put_realkV_w8
    else if(op=='g')then
      pga  => lsmpi_get_realk
      pgav => lsmpi_get_realkV_w8
    else if(op=='a')then
      pga  => lsmpi_acc_realk
      pgav => lsmpi_acc_realkV_w8
    endif

    do i = 1,arr%mode
      ro(o(i))      = i
    enddo

    if(op=='g')then
      u_o  => o(1:arr%mode)
      u_ro => ro(1:arr%mode)
    else
      u_o  => ro(1:arr%mode)
      u_ro => o(1:arr%mode)
    endif

  
    lel = fel + tl -1

    do i = 1,arr%mode
      fordims(i) = arr%dims(u_o(i))
    enddo

    comb1 = 1
    do i = 1, n2comb
      comb1   = comb1 * fordims(i)
    enddo

    comb2 = 1
    do i = n2comb + 1, arr%mode
      comb2          = comb2 * fordims(i)
    enddo

    if(arr%mode==4.and.n2comb==3.and.o(1)==1.and.o(2)==2.and.o(3)==3.and..not.deb)then
      !ATTENTION ONLY WORKS IF TL <= cons_el_in_t --> always given if order = 1,2,3,4
      !if modification needed for other types, compare the else if statement
      !where n2comb==2, this has been implemented generally

      cons_el_in_t = 1_long
      do i = 1, 4
        if(o(i)==i)then
          cons_el_in_t = cons_el_in_t * arr%tdim(i)
        else
          exit
        endif
      enddo
      cons_els = min(cons_el_in_t,int(tl,kind=8))

      st_tiling = 4
      do i = 1, 4
        if(arr%ntpm(i)/=1)then
          st_tiling = i
          exit
        endif
      enddo
      for4 => fx(4) 
      tidx = 1
      mult1 = arr%ntpm(1) * arr%ntpm(2)
      mult2 = arr%ntpm(1) * arr%ntpm(2) * arr%ntpm(3)
      !print *,"NEWOPTION",cons_els,st_tiling
       


      !precalculate tile dimensions and positions
      call mem_alloc(tinfo,arr%ntiles,7)
      do ctidx = 1, arr%ntiles
        tinfo(ctidx,1) = get_residence_of_tile(ctidx,arr)
        call get_tile_dim(tinfo(ctidx,2:5),arr,ctidx)
        tinfo(ctidx,6) = tinfo(ctidx,2) * tinfo(ctidx,3)
        tinfo(ctidx,7) = tinfo(ctidx,2) * tinfo(ctidx,3) * tinfo(ctidx,4)
      enddo
      call ass_D1to2(fort,p_fort2,[tl,fordims(4)])

      if(lock_outside) then
        do c1 = 1, tl,cons_els
          call get_midx(c1+fel-1,fx(1:n2comb),fordims(1:n2comb),n2comb)
          do for4 = 1, fordims(4)
       
            do i = st_tiling, 4
              tidx(i)   = (fx(u_ro(i))-1) / arr%tdim(i) + 1
            enddo
            
            do i = 1, 4
              idxt(i)   = mod((fx(u_ro(i))-1), arr%tdim(i)) + 1
            enddo
       
            !ctidx  = get_cidx(tidx,arr%ntpm,arr%mode)
            ctidx  = tidx(1) + (tidx(2)-1) * arr%ntpm(1) + (tidx(3)-1) *&
                     & mult1 + (tidx(4)-1) * mult2
            !cidxt  = get_cidx(idxt,tinfo(ctidx,2:5),arr%mode)
            cidxt  = idxt(1) + (idxt(2)-1) * tinfo(ctidx,2) + (idxt(3)-1) *&
                     &tinfo(ctidx,6) + (idxt(4)-1) * tinfo(ctidx,7)
       
            call pgav(p_fort2(c1:c1+cons_els-1,for4),cons_els,&
            &cidxt,int(tinfo(ctidx,1),kind=ls_mpik),arr%wi(ctidx))

          enddo
        enddo
      endif

      call mem_dealloc(tinfo)
      for3 => null()
      for4 => null()


    else if(arr%mode==4.and.n2comb==2.and..not.deb)then

      !CODE FOR 2 DIMENSIONS TO COMBINE IF A 4 MODE TENSOR IS GIVEN

      ! find the index in the tiles of the first tiled dimension
      st_tiling = 4
      do i = 1, 4
        if(arr%ntpm(i)/=1)then
          st_tiling = i
          exit
        endif
      enddo


      !find the number of consecutive elements in a tile and the index where the
      !order of the two arrays differ
      cons_el_in_t = 1_long
      diff_ord = arr%mode
      do i = 1, st_tiling
        if(u_o(i)==i)then
          cons_el_in_t = cons_el_in_t * arr%tdim(i)
        else
          diff_ord = i - 1
          exit
        endif
      enddo

      !find the consecutive elements in the reduced dimension of the unfolded
      !tensor, i.e. 1 and 2 as long as the order fits to the original tensor
      cons_el_rd = 1_long
      do i = 1, 2
        if(u_o(i)==i)then
          cons_el_rd = cons_el_rd * arr%tdim(i)
        else
          exit
        endif
      enddo


      !precalculate tile dimensions and positions such, that they may be read
      !afterwards
      call mem_alloc(tinfo,arr%ntiles,7)
      do ctidx = 1, arr%ntiles
        tinfo(ctidx,1) = get_residence_of_tile(ctidx,arr)
        call get_tile_dim(tinfo(ctidx,2:5),arr,ctidx)
        tinfo(ctidx,6) = tinfo(ctidx,2) * tinfo(ctidx,3)
        tinfo(ctidx,7) = tinfo(ctidx,2) * tinfo(ctidx,3) * tinfo(ctidx,4)
      enddo
      for3 => fx(3)
      for4 => fx(4) 
      tidx = 1
      tlidx = 1
      mult1 = arr%ntpm(1) * arr%ntpm(2)
      mult2 = arr%ntpm(1) * arr%ntpm(2) * arr%ntpm(3)
   
      !find the index of the first element (= first element of the combined two
      !first dimensions after the unfolding and splitting in stripes across the
      !nodes)
      fx=1
      call get_midx(fel,fx(1:n2comb),fordims(1:n2comb),n2comb)
      do i = 1, 4
        idxt(i)   = mod((fx(u_ro(i))-1), arr%tdim(i)) + 1
      enddo

      !DETERMINE THE PARTS: If more than one element can be transferred at a
      !time the chunks can be of two different sizes, depending on the overlap
      !of the different consecutive numbers of elements. the lengths of these
      !different parts is determined in the following
      part1 = 1
      part2 = 1
      !print *,infpar%lg_mynum,"RETARDO1:",cons_el_rd,tl
      if(cons_el_rd<tl)then
        cons_els = cons_el_rd
        do i = 1,min(diff_ord,2)
          split_in = i
          if(i==min(diff_ord,2))then
            part1 = part1 * (arr%tdim(i) - idxt(i) + 1)
            part2 = part2 * (idxt(i) - 1)
            
            exit
          else
            part1 = part1 * arr%tdim(i)
            part2 = part2 * arr%tdim(i)
          endif
        enddo
        !print *,infpar%lg_mynum,"first:",tl,arr%tdim(1), idxt(1),part1,part2
      else
        cons_els = tl
        
        part1 = min(arr%tdim(1) - idxt(1) + 1,tl)
        part2 = tl - part1
        !print *,infpar%lg_mynum,"second:",tl,arr%tdim(1), idxt(1),part1,part2
      endif

      tl_max = (tl / cons_els) * cons_els
      tl_mod = mod(tl ,cons_els)

      !do i=0,infpar%lg_nodtot-1
      !call lsmpi_barrier(infpar%lg_comm)
      !if(i==infpar%lg_mynum)then
      !  print *,i,"YAYYAYYAYYYAAAAAAYYYYYY:"
      !  print *,fel,st_tiling,cons_el_in_t,diff_ord,cons_el_rd
      !  print *,"tl",tl,"p1",part1,"p2",part2,"tlmax",tl_max,"tlmod",tl_mod
      !  print *,"fx     :",fx
      !  print *,"idxt   :",idxt
      !  print *,"a tdim :",arr%tdim
      !  print *,"a dims :",arr%dims
      !  print *,"fordims:",fordims
      !  print *,"u_o    :",u_o
      !  !stop 0 
      !endif
      !call lsmpi_barrier(infpar%lg_comm)
      !enddo


      call ass_D1to3(fort,p_fort3,[tl,fordims(3),fordims(4)])

      if(lock_outside) then
        
        !IF ONLY ONE ELEMENT CAN BE TRANSFERRED AT A TIME
        if(cons_els==1)then

          do c1 = 1, tl_max
            call get_midx(c1+fel-1,fx(1:n2comb),fordims(1:n2comb),n2comb)
            do for4 = 1, fordims(4)
              do for3 = 1, fordims(3)
        
                do i = st_tiling, 4
                  tidx(i)   = (fx(u_ro(i))-1) / arr%tdim(i) + 1
                enddo
               
                do i = 1, 4
                  idxt(i)   = mod((fx(u_ro(i))-1), arr%tdim(i)) + 1
                enddo
        
                ctidx  = tidx(1) + (tidx(2)-1) * arr%ntpm(1) + (tidx(3)-1) *&
                         & mult1 + (tidx(4)-1) * mult2

                cidxt  = idxt(1) + (idxt(2)-1) * tinfo(ctidx,2) + (idxt(3)-1) *&
                         &tinfo(ctidx,6) + (idxt(4)-1) * tinfo(ctidx,7)
        
                call pga(p_fort3(c1,for3,for4),cidxt,int(tinfo(ctidx,1),kind=ls_mpik),arr%wi(ctidx))
              enddo
            enddo
          enddo

        !IF MORE THAN ONE ELEMENT CAN BE TRANSFERRED AT A TIME
        else

          do c1 = 1, tl_max, cons_els
            call get_midx(c1+fel-1,fx(1:n2comb),fordims(1:n2comb),n2comb)
            do for4 = 1, fordims(4)
              flx(4) = fx(4) 
              do for3 = 1, fordims(3)
                do i = st_tiling, 4
                  tidx(i)   = (fx(u_ro(i))-1)  / arr%tdim(i) + 1
                enddo
                do i = 1, 4
                  idxt(i)   = mod((fx(u_ro(i))-1), arr%tdim(i)) + 1
                enddo
        
                ctidx  = tidx(1) + (tidx(2)-1) * arr%ntpm(1) + (tidx(3)-1) *&
                         & mult1 + (tidx(4)-1) * mult2
                cidxt  = idxt(1) + (idxt(2)-1) * tinfo(ctidx,2) + (idxt(3)-1) *&
                         &tinfo(ctidx,6) + (idxt(4)-1) * tinfo(ctidx,7)
        
                !FOR PART 1
                call pgav(p_fort3(c1:c1+part1-1,for3,for4),part1,&
                &cidxt,int(tinfo(ctidx,1),kind=ls_mpik),arr%wi(ctidx))
              enddo
            enddo
          enddo

          if(part2/=0)then
            do c1 = part1 + 1, tl_max, cons_els
              call get_midx(c1+fel-1,fx(1:n2comb),fordims(1:n2comb),n2comb)
              do for4 = 1, fordims(4)
                flx(4) = fx(4) 
                do for3 = 1, fordims(3)
                  do i = st_tiling, 4
                    tidx(i)   = (fx(u_ro(i))-1)  / arr%tdim(i) + 1
                  enddo
                  do i = 1, 4
                    idxt(i)   = mod((fx(u_ro(i))-1), arr%tdim(i)) + 1
                  enddo
         
                  ctidx  = tidx(1) + (tidx(2)-1) * arr%ntpm(1) + (tidx(3)-1) *&
                           & mult1 + (tidx(4)-1) * mult2
                  cidxt  = idxt(1) + (idxt(2)-1) * tinfo(ctidx,2) + (idxt(3)-1) *&
                           &tinfo(ctidx,6) + (idxt(4)-1) * tinfo(ctidx,7)
         
                  !FOR PART 2
                  call pgav(p_fort3(c1:c1+part2-1,for3,for4),part2,&
                  &cidxt,int(tinfo(ctidx,1),kind=ls_mpik),arr%wi(ctidx))
                enddo
              enddo
            enddo
          endif


          if(tl_mod/=0)then
            modp1=min(tl_mod,part1)
            modp2=tl_mod - modp1
            call get_midx(tl_max+fel,fx(1:n2comb),fordims(1:n2comb),n2comb)
            do for4 = 1, fordims(4)
              do for3 = 1, fordims(3)
          
                do i = st_tiling, 4
                  tidx(i)   = (fx(u_ro(i))-1) / arr%tdim(i) + 1
                enddo
               
                do i = 1, 4
                  idxt(i)   = mod((fx(u_ro(i))-1), arr%tdim(i)) + 1
                enddo
          
                ctidx  = tidx(1) + (tidx(2)-1) * arr%ntpm(1) + (tidx(3)-1) *&
                         & mult1 + (tidx(4)-1) * mult2
                cidxt  = idxt(1) + (idxt(2)-1) * tinfo(ctidx,2) + (idxt(3)-1) *&
                         &tinfo(ctidx,6) + (idxt(4)-1) * tinfo(ctidx,7)

                call pgav(p_fort3(tl_max+1:tl_max+modp1,for3,for4),modp1,&
                &cidxt,int(tinfo(ctidx,1),kind=ls_mpik),arr%wi(ctidx))
              enddo
            enddo

            if(tl_mod>part1)then
              call get_midx(tl_max+fel+modp1,fx(1:n2comb),fordims(1:n2comb),n2comb)
              do for4 = 1, fordims(4)
                do for3 = 1, fordims(3)
           
                  do i = st_tiling, 4
                    tidx(i)   = (fx(u_ro(i))-1) / arr%tdim(i) + 1
                  enddo
                 
                  do i = 1, 4
                    idxt(i)   = mod((fx(u_ro(i))-1), arr%tdim(i)) + 1
                  enddo
           
                  ctidx  = tidx(1) + (tidx(2)-1) * arr%ntpm(1) + (tidx(3)-1) *&
                           & mult1 + (tidx(4)-1) * mult2
                  cidxt  = idxt(1) + (idxt(2)-1) * tinfo(ctidx,2) + (idxt(3)-1) *&
                           &tinfo(ctidx,6) + (idxt(4)-1) * tinfo(ctidx,7)
             
                  call pgav(p_fort3(tl_max+modp1+1:tl_max+modp2,for3,for4),modp2,&
                  &cidxt,int(tinfo(ctidx,1),kind=ls_mpik),arr%wi(ctidx))
                enddo
              enddo
            endif
          endif

        endif

      else
      !CASE arr%mode==4,n2comb==2 and lock_outside = .false.
        if(cons_els==1)then

          do c1 = 1, tl
            call get_midx(c1+fel-1,fx(1:n2comb),fordims(1:n2comb),n2comb)
            do for4 = 1, fordims(4)
              do for3 = 1, fordims(3)
        
                do i = st_tiling, 4
                  tidx(i)   = (fx(u_ro(i))-1) / arr%tdim(i) + 1
                enddo
               
                do i = 1, 4
                  idxt(i)   = mod((fx(u_ro(i))-1), arr%tdim(i)) + 1
                enddo
        
                ctidx  = tidx(1) + (tidx(2)-1) * arr%ntpm(1) + (tidx(3)-1) *&
                         & mult1 + (tidx(4)-1) * mult2
                cidxt  = idxt(1) + (idxt(2)-1) * tinfo(ctidx,2) + (idxt(3)-1) *&
                         &tinfo(ctidx,6) + (idxt(4)-1) * tinfo(ctidx,7)
        
                call lsmpi_win_lock(int(tinfo(ctidx,1),kind=ls_mpik),arr%wi(ctidx),'s')
                call pga(p_fort3(c1,for3,for4),cidxt,int(tinfo(ctidx,1),kind=ls_mpik),arr%wi(ctidx))
                call lsmpi_win_unlock(int(tinfo(ctidx,1),kind=ls_mpik),arr%wi(ctidx))
              enddo
            enddo
          enddo

        else

          do c1 = 1, tl_max, cons_els
            call get_midx(c1+fel-1,fx(1:n2comb),fordims(1:n2comb),n2comb)
            do for4 = 1, fordims(4)
              flx(4) = fx(4) 
              do for3 = 1, fordims(3)
                do i = st_tiling, 4
                  tidx(i)   = (fx(u_ro(i))-1)  / arr%tdim(i) + 1
                enddo
                do i = 1, 4
                  idxt(i)   = mod((fx(u_ro(i))-1), arr%tdim(i)) + 1
                enddo
        
                ctidx  = tidx(1) + (tidx(2)-1) * arr%ntpm(1) + (tidx(3)-1) *&
                         & mult1 + (tidx(4)-1) * mult2
                cidxt  = idxt(1) + (idxt(2)-1) * tinfo(ctidx,2) + (idxt(3)-1) *&
                         &tinfo(ctidx,6) + (idxt(4)-1) * tinfo(ctidx,7)
        
                !FOR PART 1
                call lsmpi_win_lock(int(tinfo(ctidx,1),kind=ls_mpik),arr%wi(ctidx),'s')
                call pgav(p_fort3(c1:c1+part1-1,for3,for4),part1,&
                &cidxt,int(tinfo(ctidx,1),kind=ls_mpik),arr%wi(ctidx))
                call lsmpi_win_unlock(int(tinfo(ctidx,1),kind=ls_mpik),arr%wi(ctidx))
              enddo
            enddo
          enddo

          if(part2/=0)then
            do c1 = part1 + 1, tl_max, cons_els
              call get_midx(c1+fel-1,fx(1:n2comb),fordims(1:n2comb),n2comb)
              do for4 = 1, fordims(4)
                flx(4) = fx(4) 
                do for3 = 1, fordims(3)
                  do i = st_tiling, 4
                    tidx(i)   = (fx(u_ro(i))-1)  / arr%tdim(i) + 1
                  enddo
                  do i = 1, 4
                    idxt(i)   = mod((fx(u_ro(i))-1), arr%tdim(i)) + 1
                  enddo
         
                  ctidx  = tidx(1) + (tidx(2)-1) * arr%ntpm(1) + (tidx(3)-1) *&
                           & mult1 + (tidx(4)-1) * mult2
                  cidxt  = idxt(1) + (idxt(2)-1) * tinfo(ctidx,2) + (idxt(3)-1) *&
                           &tinfo(ctidx,6) + (idxt(4)-1) * tinfo(ctidx,7)
         
                  !FOR PART 2
                  call lsmpi_win_lock(int(tinfo(ctidx,1),kind=ls_mpik),arr%wi(ctidx),'s')
                  call pgav(p_fort3(c1:c1+part2-1,for3,for4),part2,&
                  &cidxt,int(tinfo(ctidx,1),kind=ls_mpik),arr%wi(ctidx))
                  call lsmpi_win_unlock(int(tinfo(ctidx,1),kind=ls_mpik),arr%wi(ctidx))
                enddo
              enddo
            enddo
          endif


          if(tl_mod/=0)then
            modp1=min(tl_mod,part1)
            modp2=tl_mod - modp1
            call get_midx(tl_max+fel,fx(1:n2comb),fordims(1:n2comb),n2comb)
            do for4 = 1, fordims(4)
              do for3 = 1, fordims(3)
          
                do i = st_tiling, 4
                  tidx(i)   = (fx(u_ro(i))-1) / arr%tdim(i) + 1
                enddo
               
                do i = 1, 4
                  idxt(i)   = mod((fx(u_ro(i))-1), arr%tdim(i)) + 1
                enddo
          
                ctidx  = tidx(1) + (tidx(2)-1) * arr%ntpm(1) + (tidx(3)-1) *&
                         & mult1 + (tidx(4)-1) * mult2
                cidxt  = idxt(1) + (idxt(2)-1) * tinfo(ctidx,2) + (idxt(3)-1) *&
                         &tinfo(ctidx,6) + (idxt(4)-1) * tinfo(ctidx,7)

                call lsmpi_win_lock(int(tinfo(ctidx,1),kind=ls_mpik),arr%wi(ctidx),'s')
                call pgav(p_fort3(tl_max+1:tl_max+part1,for3,for4),modp1,&
                &cidxt,int(tinfo(ctidx,1),kind=ls_mpik),arr%wi(ctidx))
                call lsmpi_win_unlock(int(tinfo(ctidx,1),kind=ls_mpik),arr%wi(ctidx))
              enddo
            enddo

            if(tl_mod>part1)then
              call get_midx(tl_max+fel+modp1,fx(1:n2comb),fordims(1:n2comb),n2comb)
              do for4 = 1, fordims(4)
                do for3 = 1, fordims(3)
           
                  do i = st_tiling, 4
                    tidx(i)   = (fx(u_ro(i))-1) / arr%tdim(i) + 1
                  enddo
                 
                  do i = 1, 4
                    idxt(i)   = mod((fx(u_ro(i))-1), arr%tdim(i)) + 1
                  enddo
           
                  ctidx  = tidx(1) + (tidx(2)-1) * arr%ntpm(1) + (tidx(3)-1) *&
                           & mult1 + (tidx(4)-1) * mult2
                  cidxt  = idxt(1) + (idxt(2)-1) * tinfo(ctidx,2) + (idxt(3)-1) *&
                           &tinfo(ctidx,6) + (idxt(4)-1) * tinfo(ctidx,7)
             
                  call lsmpi_win_lock(int(tinfo(ctidx,1),kind=ls_mpik),arr%wi(ctidx),'s')
                  call pgav(p_fort3(tl_max+modp1+1:tl_max+tl_mod,for3,for4),modp2,&
                  &cidxt,int(tinfo(ctidx,1),kind=ls_mpik),arr%wi(ctidx))
                  call lsmpi_win_unlock(int(tinfo(ctidx,1),kind=ls_mpik),arr%wi(ctidx))
                enddo
              enddo
            endif
          endif

        endif
      
      endif

      call mem_dealloc(tinfo)
      for3 => null()
      for4 => null()

    else
 
      !ONLY PRINT IF DEBUG IS NOT GIVEN, ELSE THE USER IS ASSUMED TO KNOW THAT
      !IT IS SLOWER
      if(.not.deb)print *,"WARINING(array_two_dim_1batch):this is a slow fallback option"

      do c2 = 1, comb2
        fx = 0
        call get_midx(c2,fx(n2comb+1:arr%mode),fordims(n2comb+1:arr%mode),arr%mode-n2comb)
 
        do c1 = fel, lel
          call get_midx(c1,fx(1:n2comb),fordims(1:n2comb),n2comb)
 
          !get the information about the required index in the context of the
          !array, i.e. which tile and which index in the tile
          !also get the position of the index in the batched matrix

          do i = 1, arr%mode
            tidx(i)   = (fx(u_ro(i))-1) / arr%tdim(i) + 1
            idxt(i)   = mod((fx(u_ro(i))-1), arr%tdim(i)) + 1
          enddo

          !get the combined indices, find all positions
          ctidx  = get_cidx(tidx,arr%ntpm,arr%mode)
          call get_tile_dim(tsze,arr,ctidx)
          cidxt  = get_cidx(idxt,tsze,arr%mode)
          source = get_residence_of_tile(ctidx,arr)
          cidxf  = get_cidx([c1-fel+1,c2],[tl,comb2],2)
     
          !get the element from the correct place to the correct place
          if(.not.lock_outside)call arr_lock_win(arr,ctidx,'s')
          call pga(fort(cidxf),cidxt,source,arr%wi(ctidx))
          if(.not.lock_outside)call arr_unlock_win(arr,ctidx)

        enddo
      enddo
    endif

    u_o  => null()
    u_ro => null()
    pga  => null()
#else
    call lsquit("ERROR(array_gather_to_two_dim_1batch):this routine is F2003 only",-1)
#endif
#else
    call lsquit("ERROR(array_gather_to_two_dim_1batch):this routine is MPI only",-1)
#endif
  end subroutine array_two_dim_1batch

  subroutine array_two_dim_2batch(arr,o,op,fort,n2comb,fel,tl,lock_outside)
    implicit none
    type(array),intent(in) :: arr
    real(realk),intent(inout) :: fort(*)
    character, intent(in) :: op
    integer, intent(in),target :: o(arr%mode)
    integer, intent(in) :: fel,tl,n2comb
    logical, intent(in) :: lock_outside
    integer :: fordims(arr%mode)
    integer,target :: fx(arr%mode)
    integer,target :: flx(arr%mode)
    integer :: oldidx(arr%mode)
    integer :: i,lel
    integer,target :: ro(arr%mode)
    integer :: comb1,comb2,c1,c2
    integer :: tidx(arr%mode),idxt(arr%mode)
    integer :: tlidx(arr%mode),lidxt(arr%mode)
    integer :: ctidx, cidxt, cidxf, st_tiling
    integer,pointer :: u_o(:),u_ro(:),tinfo(:,:)
    integer,pointer ::for3,for4
    real(realk), pointer :: p_fort3(:,:,:),p_fort2(:,:)
    integer :: tsze(arr%mode),mult1,mult2,dummy
    integer(kind=8) :: cons_el_in_t,cons_els,tl_max,tl_mod
    integer(kind=8) :: cons_el_rd
    integer(kind=8) :: part1,part2,split_in, diff_ord,modp1,modp2
    logical :: goto_default
#ifdef VAR_MPI
#ifdef COMPILER_UNDERSTANDS_FORTRAN_2003
    procedure(put_acc_el), pointer :: pga => null()
    procedure(put_acc_vec), pointer :: pgav => null()
    integer(kind=ls_mpik) :: source

    goto_default = .false.

    if( op/='a'.and.op/='p'.and.op/='g')then
      call lsquit("ERROR(array_two_dim_2batch):unknown choice of operator",-1)
    endif 
    if(op=='p')then
      pga  => lsmpi_put_realk
      pgav => lsmpi_put_realkV_w8
    else if(op=='g')then
      pga  => lsmpi_get_realk
      pgav => lsmpi_get_realkV_w8
    else if(op=='a')then
      pga  => lsmpi_acc_realk
      pgav => lsmpi_acc_realkV_w8
    endif

    do i = 1,arr%mode
      ro(o(i))      = i
    enddo

    if(op=='g')then
      u_o  => o(1:arr%mode)
      u_ro => ro(1:arr%mode)
    else
      u_o  => ro(1:arr%mode)
      u_ro => o(1:arr%mode)
    endif

  
    lel = fel + tl -1

    do i = 1,arr%mode
      fordims(i) = arr%dims(u_o(i))
    enddo

    comb1 = 1
    do i = 1, arr%mode-n2comb
      comb1   = comb1 * fordims(i)
    enddo

    comb2 = 1
    do i = arr%mode-n2comb+1 , arr%mode
      comb2          = comb2 * fordims(i)
    enddo
    
    if(arr%ntpm(1)>1)then
      goto_default = .true.
      print *,"WARNING(array_two_dim_2batch):going to slow default option"
    endif

    if(arr%mode==4.and.n2comb==3.and.o(1)==1.and..not.goto_default)then

      cons_el_in_t = 1_long
      do i = 1, 4
        if(o(i)==i)then
          cons_el_in_t = cons_el_in_t * arr%tdim(i)
        else
          exit
        endif
      enddo

      cons_els = int(arr%tdim(1),kind=8)

      st_tiling = 4
      do i = 1, 4
        if(arr%ntpm(i)/=1)then
          st_tiling = i
          exit
        endif
      enddo
      tidx = 1
      mult1 = arr%ntpm(1) * arr%ntpm(2)
      mult2 = arr%ntpm(1) * arr%ntpm(2) * arr%ntpm(3)


      !precalculate tile dimensions and positions
      call mem_alloc(tinfo,arr%ntiles,7)
      do ctidx = 1, arr%ntiles
        tinfo(ctidx,1) = get_residence_of_tile(ctidx,arr)
        call get_tile_dim(tinfo(ctidx,2:5),arr,ctidx)
        tinfo(ctidx,6) = tinfo(ctidx,2) * tinfo(ctidx,3)
        tinfo(ctidx,7) = tinfo(ctidx,2) * tinfo(ctidx,3) * tinfo(ctidx,4)
      enddo
      call ass_D1to2(fort,p_fort2,[tl,fordims(4)])

      if(lock_outside) then
        do c1 = 1, tl
          fx(1)=1
          call get_midx(c1+fel-1,fx(2:4),fordims(2:4),n2comb)

          do i = 1, 4
            tidx(i)   = (fx(u_ro(i))-1) / arr%tdim(i) + 1
          enddo
          
          do i = 1, 4
            idxt(i)   = mod((fx(u_ro(i))-1), arr%tdim(i)) + 1
          enddo

          ctidx  = tidx(1) + (tidx(2)-1) * arr%ntpm(1) + (tidx(3)-1) *&
                   & mult1 + (tidx(4)-1) * mult2
          cidxt  = idxt(1) + (idxt(2)-1) * tinfo(ctidx,2) + (idxt(3)-1) *&
                   &tinfo(ctidx,6) + (idxt(4)-1) * tinfo(ctidx,7)
       
          call pgav(fort(1+(c1-1)*cons_els:c1*cons_els),cons_els,&
          &cidxt,int(tinfo(ctidx,1),kind=ls_mpik),arr%wi(ctidx))

        enddo
      else
        do c1 = 1, tl
          fx(1)=1
          call get_midx(c1+fel-1,fx(2:4),fordims(2:4),n2comb)

          do i = st_tiling, 4
            tidx(i)   = (fx(u_ro(i))-1) / arr%tdim(i) + 1
          enddo
          
          do i = 1, 4
            idxt(i)   = mod((fx(u_ro(i))-1), arr%tdim(i)) + 1
          enddo

          ctidx  = tidx(1) + (tidx(2)-1) * arr%ntpm(1) + (tidx(3)-1) *&
                   & mult1 + (tidx(4)-1) * mult2
          cidxt  = idxt(1) + (idxt(2)-1) * tinfo(ctidx,2) + (idxt(3)-1) *&
                   &tinfo(ctidx,6) + (idxt(4)-1) * tinfo(ctidx,7)
       
          call arr_lock_win(arr,ctidx,'s')
          call pgav(fort(1+(c1-1)*cons_els:c1*cons_els),cons_els,&
          &cidxt,int(tinfo(ctidx,1),kind=ls_mpik),arr%wi(ctidx))
          call arr_unlock_win(arr,ctidx)

        enddo
      endif

      call mem_dealloc(tinfo)
      for3 => null()
      for4 => null()

    else
 
      print *,"WARINING(array_two_dim_2batch):this is a slow fallback option"

      do c1 = 1, comb1
        fx = 0
        call get_midx(c1,fx(1:arr%mode-n2comb),fordims(1:arr%mode-n2comb),arr%mode-n2comb)
        print *,c1,fx
        do c2 = fel, lel
          call get_midx(c2,fx(arr%mode-n2comb+1:arr%mode),fordims(arr%mode-n2comb+1:arr%mode),n2comb)
 
          !get the information about the required index in the context of the
          !array, i.e. which tile and which index in the tile
          !also get the position of the index in the batched matrix

          do i = 1, arr%mode
            tidx(i)   = (fx(u_ro(i))-1) / arr%tdim(i) + 1
            idxt(i)   = mod((fx(u_ro(i))-1), arr%tdim(i)) + 1
          enddo

          !get the combined indices, find all positions
          ctidx  = get_cidx(tidx,arr%ntpm,arr%mode)
          call get_tile_dim(tsze,arr,ctidx)
          cidxt  = get_cidx(idxt,tsze,arr%mode)
          source = get_residence_of_tile(ctidx,arr)
          cidxf  = get_cidx([c1,c2-fel+1],[comb1,tl],2)
     
          !get the element from the correct place to the correct place
          if(.not.lock_outside)call arr_lock_win(arr,ctidx,'s')
          call pga(fort(cidxf),cidxt,source,arr%wi(ctidx))
          if(.not.lock_outside)call arr_unlock_win(arr,ctidx)

        enddo
      enddo
    endif

    u_o  => null()
    u_ro => null()
    pga  => null()
#else
    call lsquit("ERROR(array_gather_to_two_dim_2batch):this routine is F2003 only",-1)
#endif
#else
    call lsquit("ERROR(array_gather_to_two_dim_2batch):this routine is MPI only",-1)
#endif
  end subroutine array_two_dim_2batch



  subroutine add_data2tiled_lowmem(arr,mult,A,dims,mode)
    implicit none
    type(array),intent(inout) :: arr
    real(realk),intent(in) :: A(*),mult
    integer,intent(in) :: mode, dims(mode)
    integer :: fib,lt,ce,j,step,mod_step,iter,nccblocks,st
    integer(kind=ls_mpik) :: nnod, me, dest, assert,ierr, act_step
    integer :: loc_ti,comp_ti,comp_el,i,nelms,fe_in_block
    integer, pointer :: elm_in_tile(:),in_tile_mode(:),orig_addr(:),remote_td(:)
    logical :: pdm
    assert = 0
    pdm=(arr%itype==TILED_DIST)

    !TRY TO INCLUDE MPI_PUT FOR THAT OPERATION,SO MAYBE MASTER_SLAVE DEPENDENCE
    !IN THIS ROUTINE IS GONE

    me = 0
    nnod=1
#ifdef VAR_MPI
    me=infpar%lg_mynum
    nnod=infpar%lg_nodtot
#endif
    !begin with sanity checks
    if(arr%mode/=mode)then
      print *,"ERROR(add_data2tiled_lowmem):mode of array does not match mode of tiled_array"
      stop 1
    endif
    do i=1,mode
      if(arr%dims(i)/=dims(i))then
        print *,"ERROR(add_data2tiled_lowmem):dims in input do not match dims of tiled_array"
        stop 1
      endif
    enddo
    ! corresponding elements
    call mem_alloc(elm_in_tile,arr%mode)
    call mem_alloc(in_tile_mode,arr%mode)
    call mem_alloc(orig_addr,arr%mode)
    call mem_alloc(remote_td,arr%mode)

    !find consecutive elements in both full and tiled matrix according to tile
    !dimensions --> determines largest memory blocks that can be transferred in
    !one go --> step
    step=arr%tdim(1)
    do i=1,arr%mode
      if(arr%tdim(i)==arr%dims(i).and.i<arr%mode)then
        step=step*arr%tdim(i+1)
      else
        exit
      endif
    enddo
    !determine the truncated block size --> mod_step is first the number of
    !elements in dimensions 1..i, this is used to determine the number of
    !iterations needed with the determined step size, and the mod gives the
    !remainder 
    mod_step=1
    do j=1,i
      mod_step = mod_step * arr%dims(j)
    enddo
    ! determine how many blocks, including truncated blocks fit into the
    ! (eventually comnbined) dimensions of the full matrix for if none
    ! of the following modes are considered
    iter = mod_step/step 
    mod_step = mod(mod_step,step)
    if(mod_step>0)iter=iter+1
    if(mod_step==0)mod_step=step
    ! determine how many consecutive blocks are in the full matrix, i.e.
    ! multiplying iter=number of blocks in the (combined) dimension by the
    ! number of elements in the non-consecutive modes of the array
    nccblocks = iter
    do j=i+1,arr%mode
      nccblocks = nccblocks * arr%dims(j)
    enddo
    !print *,step,mod_step,iter,nccblocks,arr%nelms
    
    if(mult/=1.0E0_realk)call dscal(nelms,mult,A,1)
    fe_in_block=1
    do i=1,nccblocks
      !print *,me,"in round", i
      if(mod(i,iter)>0) act_step=step
      if(mod(i,iter)==0)act_step=mod_step
      ! get the position for the tile of the first element in the block, the previous
      ! treatment ensures that all the following act_step elements are
      ! consecutive in the same tile
      call get_midx(fe_in_block,orig_addr,arr%dims,arr%mode)

      do j=1,arr%mode
        in_tile_mode(j)=(orig_addr(j)-1)/arr%tdim(j) + 1
      enddo
      comp_ti=get_cidx(in_tile_mode,arr%ntpm,arr%mode)
      !check where the current tile resides and jump the following steps if not
      !master where the full matrix resides or the destination slave
      !dest = mod(comp_ti-1+arr%offset,nnod) 
      if(pdm)dest = get_residence_of_tile(comp_ti,arr) 
    
      !get the dimensions of the remote tile
      call get_tile_dim(remote_td,comp_ti,arr%dims,arr%tdim,arr%mode)

      !now get position of the first element of the batch in the current tile
      do j=1,arr%mode
        elm_in_tile(j) = mod(orig_addr(j)-1,arr%tdim(j)) + 1
      enddo

      !get the one index element number for the remote tile
      comp_el=get_cidx(elm_in_tile,remote_td,arr%mode)

      !copy data to the identified places
      if(pdm)then
#ifdef VAR_MPI
        call lsmpi_win_lock(dest,arr%wi(comp_ti),'e')
        call lsmpi_acc(A(fe_in_block:fe_in_block+act_step-1),act_step,comp_el,dest,arr%wi(comp_ti))
        call lsmpi_win_unlock(dest,arr%wi(comp_ti))
#endif
      else
        call dcopy(act_step,A(fe_in_block),1,arr%ti(comp_ti)%t(comp_el),1)
      endif
      fe_in_block=fe_in_block + act_step
    enddo
    if(mult/=1.0E0_realk)call dscal(nelms,1.0E0_realk/mult,A,1)
    call mem_dealloc(remote_td)
    call mem_dealloc(elm_in_tile)
    call mem_dealloc(in_tile_mode)
    call mem_dealloc(orig_addr)
  end subroutine add_data2tiled_lowmem

  subroutine print_mem_per_node(output,allaccs,infoonmaster)
    implicit none
    integer, intent(in) :: output
    logical,intent(in)  :: allaccs
    real(realk),pointer,optional  :: infoonmaster(:)
    real(realk) :: get_mem(8)
    integer(kind=ls_mpik) :: i
    integer :: allallocd
    logical :: master
    real(realk) :: mb_acc,mb_put,mb_get,total(9),speed_acc,speed_get,speed_put
#ifdef VAR_MPI
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !NODE SPECIFIC ONE-SIDED INFO!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    mb_acc = (bytes_transferred_acc*1.0E0_realk)/(1024.0**2)
    mb_put = (bytes_transferred_put*1.0E0_realk)/(1024.0**2) 
    mb_get = (bytes_transferred_get*1.0E0_realk)/(1024.0**2)
    speed_acc=0.0E0_realk
    speed_put=0.0E0_realk
    speed_get=0.0E0_realk
    if(time_pdm_acc/=0.0E0_realk)speed_acc=mb_acc/time_pdm_acc
    if(time_pdm_put/=0.0E0_realk)speed_put=mb_put/time_pdm_put
    if(time_pdm_get/=0.0E0_realk)speed_get=mb_get/time_pdm_get
    
    master = .true.
    if(infpar%lg_mynum/=infpar%master)master=.false.

    if(.not.present(infoonmaster))then
      if(master.and..not.allaccs)then
        call pdm_array_sync(infpar%lg_comm,JOB_PRINT_MEM_INFO1)
      endif
      do i=1,infpar%lg_nodtot
        if(infpar%lg_mynum+1==i)then
          write(*,'("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")')
          write(*,'("Printing memory information for rank",I3)') infpar%lg_mynum
          write(*,'("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")')
          call array_print_memory_currents(output)
          write(*,'("")')
          write(*,'(" Printing one-sided transfer information for rank",I3)') infpar%lg_mynum
          write(*,'(" ***************************************************")')
          write(*,'(I6," acc: ",f15.4," MB in ",f15.4," s, bandwidth ",f15.4," MB/s")') &
          &nmsg_acc,mb_acc,time_pdm_acc,speed_acc
          write(*,'(I6," put: ",f15.4," MB in ",f15.4," s, bandwidth ",f15.4," MB/s")') &
          &nmsg_put,mb_put,time_pdm_put,speed_put
          write(*,'(I6," get: ",f15.4," MB in ",f15.4," s, bandwidth ",f15.4," MB/s")') &
          &nmsg_get,mb_get,time_pdm_get,speed_get
          write(*,'(" currently",I4," arrays allocated")')p_arr%arrays_in_use
          write(*,'("")')
        endif
        call lsmpi_barrier(infpar%lg_comm)
      enddo
    else
      if(master.and..not.allaccs)then
        call pdm_array_sync(infpar%lg_comm,JOB_PRINT_MEM_INFO2)
      endif
      call array_print_memory_currents(output,get_mem)

      if(master)then
        infoonmaster(1:8)=get_mem(1:8)
      endif
      do i=1,infpar%lg_nodtot-1
        if(infpar%lg_mynum==i)then
          call ls_mpisendrecv(get_mem,8,infpar%lg_comm,i,infpar%master)
        endif
        if(master)call ls_mpisendrecv(infoonmaster(i*8+1:i*8+8),8,infpar%lg_comm,i,infpar%master)
      enddo
      allallocd=p_arr%arrays_in_use
      call lsmpi_local_reduction(allallocd,infpar%master)
      if(master)then
        infoonmaster(1)=0.0E0_realk
        do i=1,infpar%lg_nodtot
          infoonmaster(1)=infoonmaster(1)+infoonmaster(7+(i-1)*8)
        enddo
        ! if no memory leaks are present infooonmaster is zero
        infoonmaster(1) =infoonmaster(1) + 1.0E0_realk*allallocd
        !write (output,'("SUM OF CURRENTLY ALLOCATED ARRAYS:",I5)'),allallocd
      endif
    endif


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !TOTAL ONE-SIDED INFO!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    total(1) = mb_acc
    total(2) = mb_put
    total(3) = mb_get
    total(4) = time_pdm_acc
    total(5) = time_pdm_put
    total(6) = time_pdm_get
    total(7) = nmsg_acc*1.0E0_realk
    total(8) = nmsg_put*1.0E0_realk
    total(9) = nmsg_get*1.0E0_realk
    call lsmpi_local_reduction(total,9,infpar%master)
    if(master)then
      speed_acc=0.0E0_realk
      speed_put=0.0E0_realk
      speed_get=0.0E0_realk
      if(total(4)/=0.0E0_realk)speed_acc=total(1)/total(4)
      if(total(5)/=0.0E0_realk)speed_put=total(2)/total(5)
      if(total(6)/=0.0E0_realk)speed_get=total(3)/total(6)
      if(.not.present(infoonmaster))then
        write(*,'("")')
        write(*,'("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")')
        write(*,'("Printing one-sided transfer information summary")')
        write(*,'("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")')
        write(*,'(I9," Total acc: ",f15.4," MB in ",f15.4," s, bandwidth ",f15.8," MB/s")') &
        &int(total(7)),total(1),total(4),speed_acc
        write(*,'(I9," Total put: ",f15.4," MB in ",f15.4," s, bandwidth ",f15.8," MB/s")') &
        &int(total(8)),total(2),total(5),speed_put
        write(*,'(I9," Total get: ",f15.4," MB in ",f15.4," s, bandwidth ",f15.8," MB/s")') &
        &int(total(9)),total(3),total(6),speed_get
      endif
    endif
#endif
  end subroutine print_mem_per_node

  subroutine add_data2tiled_intiles_nobuffer(arr,A,dims,mode,o)
    implicit none
    type(array),intent(inout) :: arr
    real(realk),intent(in) :: A(*)
    integer,intent(in) :: mode, dims(mode)
    integer,intent(in) :: o(mode)
    real(realk),pointer :: buf(:)
    integer ::nnod,me
    integer :: nelmsit,i
    integer :: fullfortdims(arr%mode)
    do i=1,arr%mode
      fullfortdims(o(i)) = arr%dims(i)
    enddo
    me = 0
    nnod=1
#ifdef VAR_MPI
    me=infpar%lg_mynum
    nnod=infpar%lg_nodtot
#endif
    !begin with sanity checks
    if(arr%mode/=mode)then
      print *,"ERROR(add_data2tiled_intiles):mode of array does not match mode of tiled_array"
      stop 1
    endif
    do i=1,mode
      if(arr%dims(i)/=dims(i))then
        print *,"ERROR(add_data2tiled_intiles):dims in input do not match dims of tiled_array"
        stop 1
      endif
    enddo
    ! corresponding elements
    call mem_alloc(buf,arr%tsize)

    
    do i=1,arr%ntiles
      call get_tile_dim(nelmsit,arr,i)
#ifdef VAR_MPI
      call array_accumulate_tile_combidx_nobuff(&
      &A,fullfortdims,arr,i,o,arr%lock_set(i))
#endif
    enddo
    call mem_dealloc(buf)
  end subroutine add_data2tiled_intiles_nobuffer

  subroutine add_data2tiled_intiles_stackbuffer(arr,mult,A,dims,mode,o)
    implicit none
    type(array),intent(inout) :: arr
    real(realk),intent(in) :: A(*),mult
    integer,intent(in) :: mode, dims(mode)
    integer,intent(in) :: o(mode)
    real(realk),pointer :: buf(:)
    integer ::nnod,me
    integer :: nelmsit,i
    integer :: fullfortdims(arr%mode)
    do i=1,arr%mode
      fullfortdims(o(i)) = arr%dims(i)
    enddo
    me = 0
    nnod=1
#ifdef VAR_MPI
    me=infpar%lg_mynum
    nnod=infpar%lg_nodtot
#endif
    !begin with sanity checks
    if(arr%mode/=mode)then
      print *,"ERROR(add_data2tiled_intiles):mode of array does not match mode of tiled_array"
      stop 1
    endif
    do i=1,mode
      if(arr%dims(i)/=dims(i))then
        print *,"ERROR(add_data2tiled_intiles):dims in input do not match dims of tiled_array"
        stop 1
      endif
    enddo
    ! corresponding elements
    call mem_alloc(buf,arr%tsize)

    
    do i=1,arr%ntiles
      call tile_from_fort(mult,A,fullfortdims,arr%mode,0.0E0_realk,buf,i,arr%tdim,o)
      call get_tile_dim(nelmsit,arr,i)
#ifdef VAR_MPI
      call array_accumulate_tile(arr,i,buf,nelmsit,lock_set=arr%lock_set(i))
#endif
    enddo
    call mem_dealloc(buf)
  end subroutine add_data2tiled_intiles_stackbuffer

  subroutine add_data2tiled_intiles_explicitbuffer(arr,mult,A,dims,mode,o,wrk,iwrk)
    implicit none
    type(array),intent(inout) :: arr
    real(realk),intent(in) :: A(*),mult
    integer,intent(in) :: mode, dims(mode)
    integer,intent(in) :: o(mode)
    integer(kind=8),intent(in) :: iwrk
    real(realk), intent(inout) :: wrk(*)
    integer ::nnod,me
    integer :: nelmsit
    integer(kind=8) ::i,b,e,maxntiinwrk,mod_el
    integer :: fullfortdims(arr%mode)

    do i=1,arr%mode
      fullfortdims(o(i)) = arr%dims(i)
    enddo
#ifdef VAR_MPI

    me   = infpar%lg_mynum
    nnod = infpar%lg_nodtot

    !begin with sanity checks
    if(arr%mode/=mode)then
      print *,"ERROR(add_data2tiled_intiles_explicitbuffer):&
      &mode of array does not match mode of tiled_array"
      stop 1
    endif
    do i=1,mode
      if(arr%dims(i)/=dims(i))then
        print *,"ERROR(add_data2tiled_intiles_explicitbuffer):&
        &dims in input do not match dims of tiled_array"
        stop 1
      endif
    enddo

    !compute the maximum number of tiles to be stored in the workspace
    maxntiinwrk = int(iwrk/arr%tsize,kind=8)

    if(maxntiinwrk == 0)then

      if(mult==1.0E0_realk)then
        print *,"WARNING(add_data2tiled_intiles_explicitbuffer)&
        &:not enough space in wrk, try more nodes (in a slot) -> redirecting to _nobuffer"
        call add_data2tiled_intiles_nobuffer(arr,A,dims,mode,o)
      else
        print *,"WARNING(add_data2tiled_intiles_explicitbuffer)&
        &:not enough space in wrk, try more nodes (in a slot) -> redirecting to _stackbuffer"
        call add_data2tiled_intiles_stackbuffer(arr,mult,A,dims,mode,o)
      endif

    else
      
      do i=1,arr%ntiles

        if(arr%lock_set(i))then
          if(i>maxntiinwrk) then
            call arr_unlock_win(arr,int(i-maxntiinwrk))
          endif
        endif

        call get_tile_dim(nelmsit,arr,i)
        b = 1       + mod(i-1,maxntiinwrk) * arr%tsize
        e = nelmsit + mod(i-1,maxntiinwrk) * arr%tsize
        call tile_from_fort(mult,A,fullfortdims,arr%mode,0.0E0_realk,wrk(b),int(i),arr%tdim,o)
        call array_accumulate_tile(arr,int(i),wrk(b:e),nelmsit,lock_set=arr%lock_set(i))
      enddo
    endif

#endif
  end subroutine add_data2tiled_intiles_explicitbuffer

  subroutine cp_data2tiled_lowmem(arr,A,dims,mode)
    implicit none
    type(array),intent(inout) :: arr
    real(realk),intent(in) :: A(*)
    integer,intent(in) :: mode, dims(mode)
    integer :: fib,lt,ce,j,step,mod_step,iter,nccblocks,st
    integer(kind=ls_mpik) :: nnod, me, dest, assert,ierr, act_step
    integer :: loc_ti,comp_ti,comp_el,i,nelms,fe_in_block
    integer, pointer :: elm_in_tile(:),in_tile_mode(:),orig_addr(:),remote_td(:)
    logical :: pdm

    pdm=(arr%itype==TILED_DIST)

    !TRY TO INCLUDE MPI_PUT FOR THAT OPERATION,SO MAYBE MASTER_SLAVE DEPENDENCE
    !IN THIS ROUTINE IS GONE

    me = 0
    nnod=1
#ifdef VAR_MPI
    me=infpar%lg_mynum
    nnod=infpar%lg_nodtot
#endif
    !begin with sanity checks
    if(arr%mode/=mode)then
      print *,"ERROR(cp_data2tiled_lowmem):mode of array does not match mode of tiled_array"
      stop 1
    endif
    do i=1,mode
      if(arr%dims(i)/=dims(i))then
        print *,"ERROR(cp_data2tiled_lowmem):dims in input do not match dims of tiled_array"
        stop 1
      endif
    enddo
    ! corresponding elements
    call mem_alloc(elm_in_tile,arr%mode)
    call mem_alloc(in_tile_mode,arr%mode)
    call mem_alloc(orig_addr,arr%mode)
    call mem_alloc(remote_td,arr%mode)

    !find consecutive elements in both full and tiled matrix according to tile
    !dimensions --> determines largest memory blocks that can be transferred in
    !one go --> step
    step=arr%tdim(1)
    do i=1,arr%mode
      if(arr%tdim(i)==arr%dims(i).and.i<arr%mode)then
        step=step*arr%tdim(i+1)
      else
        exit
      endif
    enddo
    !determine the truncated block size --> mod_step is first the number of
    !elements in dimensions 1..i, this is used to determine the number of
    !iterations needed with the determined step size, and the mod gives the
    !remainder 
    mod_step=1
    do j=1,i
      mod_step = mod_step * arr%dims(j)
    enddo
    ! determine how many blocks, including truncated blocks fit into the
    ! (eventually comnbined) dimensions of the full matrix for if none
    ! of the following modes are considered
    iter = mod_step/step 
    mod_step = mod(mod_step,step)
    if(mod_step>0)iter=iter+1
    if(mod_step==0)mod_step=step
    ! determine how many consecutive blocks are in the full matrix, i.e.
    ! multiplying iter=number of blocks in the (combined) dimension by the
    ! number of elements in the non-consecutive modes of the array
    nccblocks = iter
    do j=i+1,arr%mode
      nccblocks = nccblocks * arr%dims(j)
    enddo
    !print *,step,mod_step,iter,nccblocks,arr%nelms
    
    fe_in_block=1
    do i=1,nccblocks
      if(mod(i,iter)>0) act_step=step
      if(mod(i,iter)==0)act_step=mod_step
      ! get the position for the tile of the first element in the block, the previous
      ! treatment ensures that all the following act_step elements are
      ! consecutive in the same tile
      call get_midx(fe_in_block,orig_addr,arr%dims,arr%mode)

      do j=1,arr%mode
        in_tile_mode(j)=(orig_addr(j)-1)/arr%tdim(j) + 1
      enddo
      comp_ti=get_cidx(in_tile_mode,arr%ntpm,arr%mode)

      !check where the current tile resides and jump the following steps if not
      !master where the full matrix resides or the destination slave
      !dest = mod(comp_ti-1+arr%offset,nnod) 
      if(pdm)dest = get_residence_of_tile(comp_ti,arr) 

      !get the dimensions of the remote tile
      call get_tile_dim(remote_td,comp_ti,arr%dims,arr%tdim,arr%mode)

      !now get position of the first element of the batch in the current tile
      do j=1,arr%mode
        elm_in_tile(j) = mod(orig_addr(j)-1,arr%tdim(j)) + 1
      enddo

      !get the one index element number for the remote tile
      comp_el=get_cidx(elm_in_tile,remote_td,arr%mode)

      !copy data to the identified places
      if(pdm)then
#ifdef VAR_MPI
        call lsmpi_win_lock(dest,arr%wi(comp_ti),'e')
        call lsmpi_put(A(fe_in_block:fe_in_block+act_step-1),act_step,comp_el,dest,arr%wi(comp_ti))
        call lsmpi_win_unlock(dest,arr%wi(comp_ti))
#endif
      else
        call dcopy(act_step,A(fe_in_block),1,arr%ti(comp_ti)%t(comp_el),1)
      endif
      fe_in_block=fe_in_block + act_step
    enddo
    call mem_dealloc(remote_td)
    call mem_dealloc(elm_in_tile)
    call mem_dealloc(in_tile_mode)
    call mem_dealloc(orig_addr)
  end subroutine cp_data2tiled_lowmem

  subroutine cp_data2tiled_intiles(arr,A,dims,mode,optorder)
    implicit none
    type(array),intent(inout) :: arr
    real(realk),intent(in) :: A(*)
    integer,intent(in) :: mode, dims(mode)
    integer,intent(in),optional :: optorder(mode)
    real(realk),pointer :: buf(:)
    integer ::nnod,fib,lt,ce,j,me,dest,step,act_step,mod_step,iter,nccblocks,ierr,st
    integer :: nelmsit,loc_ti,comp_ti,comp_el,i,nelms,fe_in_block, order(mode)
    integer, pointer :: elm_in_tile(:),in_tile_mode(:),orig_addr(:),remote_td(:)
    integer :: fullfortdims(arr%mode)

    !TRY TO INCLUDE MPI_PUT FOR THAT OPERATION,SO MAYBE MASTER_SLAVE DEPENDENCE
    !IN THIS ROUTINE IS GONE
    do i=1,mode
      order(i)=i
    enddo
    if(present(optorder))order=optorder
   
    do i=1,arr%mode
      fullfortdims(order(i)) = arr%dims(i)
    enddo

    me = 0
    nnod=1
#ifdef VAR_MPI
    me=infpar%lg_mynum
    nnod=infpar%lg_nodtot
#endif
    !begin with sanity checks
    if(arr%mode/=mode)then
      print *,"ERROR(cp_data2tiled_intiles):mode of array does not match mode of tiled_array"
      stop 1
    endif
    do i=1,mode
      if(arr%dims(i)/=dims(i))then
        print *,"ERROR(cp_data2tiled_intiles):dims in input do not match dims of tiled_array"
        stop 1
      endif
    enddo
    ! corresponding elements
    call mem_alloc(buf,arr%tsize)

    
    do i=1,arr%ntiles
      call tile_from_fort(1.0E0_realk,A,fullfortdims,arr%mode,0.0E0_realk,buf,i,arr%tdim,order)
      call get_tile_dim(nelmsit,arr,i)
      !copy data to the identified places
#ifdef VAR_MPI
      call array_put_tile(arr,i,buf,nelmsit)
#endif
    enddo
    call mem_dealloc(buf)
  end subroutine cp_data2tiled_intiles




  subroutine array_gatheradd_tilestofort(arr,sc,fort,nelms,nod,optorder)
    implicit none
    type(array),intent(in) :: arr
    real(realk),intent(in) :: sc
    real(realk),intent(inout) :: fort(*)
    integer(kind=long), intent(in) :: nelms
    integer(kind=ls_mpik) :: nod
    integer, intent(in), optional :: optorder(arr%mode)
    integer(kind=ls_mpik) :: src,me,nnod
    integer :: i,ltidx,order(arr%mode)
    integer :: nelintile,fullfortdim(arr%mode)
    real(realk), pointer :: tmp(:)
#ifdef VAR_MPI

    do i=1,arr%mode
      order(i)=i
    enddo
    if(present(optorder))order=optorder
   
    me=0
    nnod=1
    me=infpar%lg_mynum
    nnod=infpar%lg_nodtot
    if(nelms/=arr%nelms)call lsquit("ERROR(cp_tileddate2fort):array&
        &dimensions are not the same",DECinfo%output)

    do i = 1, arr%mode
      fullfortdim(i) = arr%dims(order(i))
    enddo

    call mem_alloc(tmp,arr%tsize)

    do i=1,arr%ntiles
      src=get_residence_of_tile(i,arr)
      if(src==me.or.nod==me)then
        call get_tile_dim(nelintile,i,arr%dims,arr%tdim,arr%mode)
        if(src==me.and.nod==me)then
          ltidx = (i - 1) /nnod + 1
          call tile_in_fort(sc,arr%ti(ltidx)%t,i,arr%tdim,&
               &1.0E0_realk,fort,fullfortdim,arr%mode,order)
        else if(src==me)then
          ltidx = (i - 1) /nnod + 1
          call lsmpi_send(arr%ti(ltidx)%t,nelintile,infpar%lg_comm,nod)
        else if(nod==me)then
          call lsmpi_recv(tmp,nelintile,infpar%lg_comm,src)
          call tile_in_fort(sc,tmp,i,arr%tdim,&
               &1.0E0_realk,fort,fullfortdim,arr%mode,order)
        endif
      endif
    enddo

    call mem_dealloc(tmp)
#else
    call lsquit("ERROR(array_gatheradd_tilestofort):this routine is MPI only",-1)
#endif
  end subroutine array_gatheradd_tilestofort



  subroutine array_gather_tilesinfort(arr,fort,nelms,nod,optorder)
    implicit none
    type(array),intent(in) :: arr
    real(realk),intent(inout) :: fort(*)
    integer(kind=long), intent(in) :: nelms
    integer(kind=ls_mpik) :: nod
    integer, intent(in), optional :: optorder(arr%mode)
    integer(kind=ls_mpik) :: src,me,nnod
    integer :: i,j,k,ltidx
    integer :: nelintile,order(arr%mode)
    integer :: fullfortdim(arr%mode)
    real(realk), pointer :: tmp(:)
#ifdef VAR_MPI
    do i = 1, arr%mode
      order(i) = i
    enddo
    if(present(optorder))order=optorder
    do i = 1, arr%mode
      fullfortdim(i) = arr%dims(order(i))
    enddo
    me=0
    nnod=1
    me=infpar%lg_mynum
    nnod=infpar%lg_nodtot
    if(nelms/=arr%nelms)call lsquit("ERROR(cp_tileddate2fort):array&
        &dimensions are not the same",DECinfo%output)

    call mem_alloc(tmp,arr%tsize)

    do i=1,arr%ntiles
      src=get_residence_of_tile(i,arr)
      if(src==me.or.nod==me)then
        call get_tile_dim(nelintile,i,arr%dims,arr%tdim,arr%mode)
        if(src==me.and.nod==me)then
          ltidx = (i - 1) /nnod + 1
          call tile_in_fort(1.0E0_realk,arr%ti(ltidx)%t,i,arr%tdim,&
                           &0.0E0_realk,fort,fullfortdim,arr%mode,order)
        else if(src==me)then
          ltidx = (i - 1) /nnod + 1
          call lsmpi_send(arr%ti(ltidx)%t,nelintile,infpar%lg_comm,nod)
        else if(nod==me)then
          call lsmpi_recv(tmp,nelintile,infpar%lg_comm,src)
          call tile_in_fort(1.0E0_realk,tmp,i,arr%tdim,&
                           &0.0E0_realk,fort,fullfortdim,arr%mode,order)
        endif
      endif
    enddo

    call mem_dealloc(tmp)
#else
    call lsquit("ERROR(array_gather_tilesinfort):this routine is MPI only",-1)
#endif
  end subroutine array_gather_tilesinfort


  subroutine array_scatter_densetotiled(arr,A,nelms,nod,optorder)
    implicit none
    type(array),intent(inout) :: arr
    real(realk),intent(in) :: A(*)
    integer(kind=long),intent(in) :: nelms
    integer(kind=ls_mpik),intent(in) :: nod
    integer,intent(in),optional :: optorder(arr%mode)
    real(realk),pointer :: buf(:)
    integer :: nelmsit,i, order(arr%mode)
    integer :: ltidx
    integer(kind=ls_mpik) :: nnod,dest,me
    integer :: fullfortdims(arr%mode)
#ifdef VAR_MPI

    !TRY TO INCLUDE MPI_PUT FOR THAT OPERATION,SO MAYBE MASTER_SLAVE DEPENDENCE
    !IN THIS ROUTINE IS GONE
    do i=1,arr%mode
      order(i)=i
    enddo
    if(present(optorder))order=optorder

    do i=1,arr%mode
      fullfortdims(order(i)) = arr%dims(i)
    enddo
   

    me = 0
    nnod=1
    me=infpar%lg_mynum
    nnod=infpar%lg_nodtot
    !begin with sanity checks
    ! corresponding elements
    call mem_alloc(buf,arr%tsize)

    
    do i=1,arr%ntiles
      dest=get_residence_of_tile(i,arr)
      call get_tile_dim(nelmsit,arr,i)
      if(dest==me.and.nod==me)then
        ltidx = (i - 1) /nnod + 1
        call tile_from_fort(1.0E0_realk,A,fullfortdims,arr%mode,&
                           &0.0E0_realk,arr%ti(ltidx)%t,i,arr%tdim,order)
      else if(nod==me)then
        call tile_from_fort(1.0E0_realk,A,fullfortdims,arr%mode,0.0E0_realk,buf,i,arr%tdim,order)
        call lsmpi_send(buf,nelmsit,infpar%lg_comm,dest)
      else if(dest==me)then
        ltidx = (i - 1) /nnod + 1
        call lsmpi_recv(arr%ti(ltidx)%t,nelmsit,infpar%lg_comm,nod)
      endif
    enddo
    call mem_dealloc(buf)
#else
    call lsquit("ERROR(array_scatter_densetotiled):this routine is MPI only",-1)
#endif

  end subroutine array_scatter_densetotiled



  !\> \brief lock all windows of a tensor from the current node with the
  !specified lock and assertion
  !\> \author Patrick Ettenhuber
  !\> \date July 2013
#ifdef VAR_MPI
  subroutine arr_lock_win(arr,ti_idx,locktype,assert)
    implicit none
    type(array) :: arr
    integer,intent(in) :: ti_idx
    character, intent(in) :: locktype
    integer(kind=ls_mpik), optional,intent(in) :: assert
    integer(kind=ls_mpik) ::node

    node=get_residence_of_tile(ti_idx,arr)
    call lsmpi_win_lock(node,arr%wi(ti_idx),locktype,ass=assert)
    arr%lock_set(ti_idx)=.true.

  end subroutine arr_lock_win

  subroutine arr_unlock_win(arr,ti_idx)
    implicit none
    type(array) :: arr
    integer,intent(in) :: ti_idx
    integer(kind=ls_mpik) :: node

    node                 = get_residence_of_tile(ti_idx,arr)
    call lsmpi_win_unlock(node,arr%wi(ti_idx))
    arr%lock_set(ti_idx) = .false.

  end subroutine arr_unlock_win

  subroutine arr_lock_wins(arr,locktype,assert)
    implicit none
    type(array) :: arr
    character, intent(in) :: locktype
    integer(kind=ls_mpik), optional,intent(in) :: assert
    integer(kind=ls_mpik) :: node
    integer :: i

    do i=1,arr%ntiles
      node            = get_residence_of_tile(i,arr)
      call lsmpi_win_lock(node,arr%wi(i),locktype,ass=assert)
      arr%lock_set(i) = .true.
    enddo

  end subroutine arr_lock_wins

  !\> \brief unlock all windows of a tensor 
  !\> \author Patrick Ettenhuber
  !\> \date July 2013
  subroutine arr_unlock_wins(arr,check)
    implicit none
    type(array) :: arr
    logical, intent(in),optional:: check
    integer(kind=ls_mpik) :: node
    integer :: i
    logical :: ch

    ch = .false.
    if(present(check))ch=check

    if(ch)then

      do i=1,arr%ntiles
        if(arr%lock_set(i))then
          node=get_residence_of_tile(i,arr)
          call lsmpi_win_unlock(node,arr%wi(i))
          arr%lock_set(i)=.false.
        endif
      enddo

    else

      do i=1,arr%ntiles
        node=get_residence_of_tile(i,arr)
        call lsmpi_win_unlock(node,arr%wi(i))
        arr%lock_set(i)=.false.
      enddo

    endif

  end subroutine arr_unlock_wins
#endif


  subroutine pn(a,n)
    implicit none
    real(realk), intent(in) :: a(*)
    integer,intent(in) :: n
    integer :: i
    real(realk) :: nrm
    nrm = 0.0E0_realk
    do i=1,n
      nrm=nrm+a(i)*a(i)
    enddo
    nrm = sqrt(nrm)
    print *,"NORM:",nrm
  end subroutine

  subroutine arr_deallocate_dense(arr)
    implicit none
    type(array) :: arr
    integer :: a
    call memory_deallocate_array_dense(arr)
    a = 1
#ifdef VAR_MPI
    a = infpar%lg_mynum + 1
#endif
    if(associated(p_arr%a(arr%addr_p_arr(a))%elm1))then
      p_arr%a(arr%addr_p_arr(a))%elm1 => null()
    endif
  end subroutine arr_deallocate_dense


  subroutine array_free_pdm(arr)
    implicit none
    type(array) :: arr
    logical     :: parent
#ifdef VAR_MPI
    parent = (infpar%parent_comm == MPI_COMM_NULL)

    if( arr%access_type==MASTER_ACCESS &
    &   .and.infpar%lg_mynum==infpar%master &
    &   .and. parent                          )then

      call pdm_array_sync(infpar%lg_comm,JOB_FREE_ARR_PDM,arr)

    endif

    if( parent .and. lspdm_use_comm_proc ) then
      call pdm_array_sync(infpar%pc_comm,JOB_FREE_ARR_PDM,arr,loc_addr=.true.)
    endif

    if( parent )then
      p_arr%free_addr_on_node(arr%addr_p_arr(infpar%lg_mynum+1))=.true.
      p_arr%arrays_in_use = p_arr%arrays_in_use - 1 
      call array_free_basic(p_arr%a(arr%addr_p_arr(infpar%lg_mynum+1))) 
    else
      p_arr%free_addr_on_node(arr%addr_loc(infpar%lg_mynum+1))=.true.
      p_arr%arrays_in_use = p_arr%arrays_in_use - 1 
      call array_free_basic(p_arr%a(arr%addr_loc(infpar%lg_mynum+1))) 
    endif
    call array_nullify_pointers(arr)
#endif
  end subroutine array_free_pdm

  subroutine get_distribution_info(arr)
    implicit none
    type(array),intent(inout) :: arr
    integer :: i,ntiles2dis
    logical :: parent
    integer(kind=ls_mpik) :: lg_me,lg_nnod,pc_me,pc_nnod,buf(2)
#ifdef VAR_MPI
    lg_me   = infpar%lg_mynum
    lg_nnod = infpar%lg_nodtot
    if( lspdm_use_comm_proc ) then
      pc_me   = infpar%pc_mynum
      pc_nnod = infpar%pc_nodtot
      buf(1)  = lg_me 
      buf(2)  = lg_nnod
      call ls_mpibcast(buf,2,infpar%master,infpar%pc_comm)
      lg_me   = buf(1)
      lg_nnod = buf(2)
    endif
 
    if(arr%access_type==NO_PDM_ACCESS.or.arr%itype==TILED)then
      arr%offset       = 0
      p_arr%new_offset = 0
      arr%nlti         = arr%ntiles
    else
      arr%offset       = p_arr%new_offset
      p_arr%new_offset = mod(p_arr%new_offset+arr%ntiles,lg_nnod)
      arr%nlti         = arr%ntiles/lg_nnod
      if(mod(arr%ntiles,lg_nnod)>mod(lg_me+lg_nnod-arr%offset,lg_nnod))arr%nlti=arr%nlti+1
    endif
#endif
  end subroutine get_distribution_info

  !> \brief routine to get a free address in the persisten array
  !> \autor Patrick Ettenhuber
  function get_free_address(occ_addr) result(addr)
    implicit none
    !> retrurn value with the address
    integer :: addr
    !> logical which tells the routine to set the value of the found address to occupied
    logical, intent(in) :: occ_addr
    if(p_arr%arrays_in_use==n_arrays)then
      call lsquit("ERROR(get_free_address):max number of arrays in p_arr allocated, change&
      & the parameter n_arrays in lsutil/lspdm_tensor_operations.F90 and recompile",-1)
    endif

    do addr=1,p_arr%arrays_in_use+1
      if(p_arr%free_addr_on_node(addr))then
        if(occ_addr)p_arr%free_addr_on_node(addr) = .false.
        return
      endif
    enddo

  end function get_free_address

  !> \brief debugging routine to check the norms of individual tiles
  !> \author Patrick Ettenhuber
  subroutine array_tiled_pdm_print_ti_nrm(arr,globtinr,whichnode,nrm) 
    implicit none
    !> input array for which to check the tile
    type(array), intent(in) :: arr
    !> global index number of the tile
    integer, intent(in) :: globtinr
    !> optional input, return value for the destination of the tile
    integer, intent(out), optional :: whichnode
    !> optional input, return value for the norm
    real(realk), intent(out), optional :: nrm
    real(realk) :: norm
    integer :: i,j,loctinr,gtnr
    integer(kind=ls_mpik) :: dest
#ifdef VAR_MPI
    gtnr=globtinr
    if(arr%access_type==MASTER_ACCESS.and.infpar%lg_mynum==infpar%master)then
      call pdm_array_sync(infpar%lg_comm,JOB_PRINT_TI_NRM,arr)
    endif
    call ls_mpibcast(gtnr,infpar%master,infpar%lg_comm)

    dest=get_residence_of_tile(gtnr,arr)
    if(present(whichnode))whichnode=dest

    if(dest==infpar%lg_mynum)then
      loctinr=(gtnr-1)/infpar%lg_nodtot + 1
      norm=0.0E0_realk
      do j=1,arr%ti(loctinr)%e
        norm = norm + arr%ti(loctinr)%t(j) * arr%ti(loctinr)%t(j)
      enddo
      call ls_mpisendrecv(norm,infpar%lg_comm,infpar%lg_mynum,infpar%master)
    endif

    if(infpar%lg_mynum==0.and.infpar%lg_mynum/=dest)then
      call ls_mpisendrecv(norm,infpar%lg_comm,dest,infpar%master)
    endif

    !if nrm is present return the squared norm, else print the norm
    if(infpar%lg_mynum==0.and.present(nrm))then
      nrm = norm
    else if(infpar%lg_mynum==0)then
      write(DECinfo%output,'("LOCAL TILE NORM ON",I3,f20.15)') dest,sqrt(norm)
    endif
#endif
  end subroutine array_tiled_pdm_print_ti_nrm

  function array_tiled_pdm_get_nrm2(arr) result(nrm)
    implicit none
    type(array), intent(in) :: arr
    real(realk) :: nrm
    integer :: i,j,should
#ifdef VAR_MPI
    if(infpar%lg_mynum==infpar%master.and.arr%access_type==MASTER_ACCESS) then
      call pdm_array_sync(infpar%lg_comm,JOB_GET_NRM2_TILED,arr)
    endif
    nrm=0.0E0_realk
    do i=1,arr%nlti
      do j=1,arr%ti(i)%e
        nrm = nrm +(arr%ti(i)%t(j) * arr%ti(i)%t(j))
      enddo
    enddo
    if(arr%access_type==MASTER_ACCESS)call lsmpi_local_reduction(nrm,infpar%master)
    if(arr%access_type==ALL_ACCESS)call lsmpi_allreduce(nrm,infpar%lg_comm)
#else
    nrm = 0.0E0_realk
#endif
  end function array_tiled_pdm_get_nrm2

   subroutine change_access_type_td(arr,totype)
     implicit none
     type(array),intent(inout) :: arr
     integer,intent(in) :: totype
#ifdef VAR_MPI
     if(totype/=REPLICATED.and.totype/=DENSE.and.totype/=TILED_DIST.and.totype/=TILED)then
       call lsquit("ERROR(change_access_type_td): wrong type given",-1)
     endif
     if(infpar%lg_mynum==infpar%master.and.arr%access_type==MASTER_ACCESS) then
       call pdm_array_sync(infpar%lg_comm,JOB_CHANGE_access_type,arr)
     endif
     arr%access_type=totype
#endif
   end subroutine change_access_type_td



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!                 ACCUMULATE TILES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> \brief direct communication routine for the accumulation of arrays,
  !> interface to the combined index routine
  !> \author Patrick Ettenhuber
  subroutine array_accumulate_tile_modeidx(arr,modidx,fort,nelms,lock_set,flush_it)
    implicit none
    !> input array for which a tile should be accumulated
    type(array),intent(in) ::arr
    !> input, the index of the tile in modular form and the number of elements
    integer,intent(in) :: modidx(arr%mode),nelms
    !> input the fortan array which should be transferred to the tile
    real(realk),intent(inout) :: fort(*)
    logical, optional, intent(in) :: lock_set,flush_it
    integer :: cidx
    cidx=get_cidx(modidx,arr%ntpm,arr%mode)
    call array_accumulate_tile(arr,cidx,fort,nelms,lock_set=lock_set,flush_it=flush_it)
  end subroutine array_accumulate_tile_modeidx
  subroutine array_acct4(arr,globtilenr,fort,nelms,lock_set,flush_it)
    implicit none
    type(array),intent(in) :: arr
    integer,intent(in) :: globtilenr
    integer(kind=4),intent(in) :: nelms
    real(realk),intent(inout) :: fort(*)
    logical, optional, intent(in) :: lock_set,flush_it
    call array_accumulate_tile_combidx4(arr,globtilenr,fort,&
    &nelms,lock_set=lock_set,flush_it=flush_it)
  end subroutine array_acct4
  subroutine array_accumulate_tile_combidx4(arr,globtilenr,fort,nelms,lock_set,flush_it)
    implicit none
    type(array),intent(in) :: arr
    integer,intent(in) :: globtilenr
    integer(kind=4),intent(in) :: nelms
    !> input the fortan array which should be transferred to the tile
    real(realk),intent(inout) :: fort(*)
    logical, optional, intent(in) :: lock_set, flush_it
    integer(kind=ls_mpik) :: dest
    logical :: ls
    real(realk) :: sta,sto
#ifdef VAR_MPI
    integer :: maxsze
    maxsze = MAX_SIZE_ONE_SIDED

    ls = .false.
    if(present(lock_set))ls=lock_set

    dest = get_residence_of_tile(globtilenr,arr)
    sta  = MPI_WTIME()

    if(.not.ls)call lsmpi_win_lock(dest,arr%wi(globtilenr),'s')
    call lsmpi_acc(fort,nelms,1,dest,arr%wi(globtilenr),maxsze,flush_it=flush_it)
    if(.not.ls)CALL lsmpi_win_unlock(dest, arr%wi(globtilenr))

    sto          = MPI_WTIME()
    time_pdm_acc = time_pdm_acc + sto - sta
    bytes_transferred_acc = bytes_transferred_acc + nelms * 8_long
    nmsg_acc     = nmsg_acc + 1
#endif
  end subroutine array_accumulate_tile_combidx4
  subroutine array_acct8(arr,globtilenr,fort,nelms,lock_set,flush_it)
    implicit none
    type(array),intent(in) :: arr
    integer,intent(in) :: globtilenr
    integer(kind=8),intent(in) :: nelms
    real(realk),intent(inout) :: fort(*)
    logical, optional, intent(in) :: lock_set,flush_it
    call array_accumulate_tile_combidx8(arr,globtilenr,fort,nelms,lock_set=lock_set,flush_it=flush_it)
  end subroutine array_acct8
  subroutine array_accumulate_tile_combidx8(arr,globtilenr,fort,nelms,lock_set,flush_it)
    implicit none
    type(array),intent(in) :: arr
    integer,intent(in) :: globtilenr
    logical, optional, intent(in) :: lock_set,flush_it
    integer(kind=8),intent(in) :: nelms
    !> input the fortan array which should be transferred to the tile
    real(realk),intent(inout) :: fort(*)
    integer(kind=ls_mpik) :: dest
    logical :: ls
    real(realk) :: sta,sto
#ifdef VAR_MPI
    integer :: maxsze
    maxsze = MAX_SIZE_ONE_SIDED

    ls = .false.
    if(present(lock_set))ls=lock_set

    dest = get_residence_of_tile(globtilenr,arr)
    sta  = MPI_WTIME()

    if(.not.ls)call lsmpi_win_lock(dest,arr%wi(globtilenr),'s')
    call lsmpi_acc(fort,nelms,1,dest,arr%wi(globtilenr),maxsze,flush_it=flush_it)
    if(.not.ls)call lsmpi_win_unlock(dest,arr%wi(globtilenr))

    sto          = MPI_WTIME()
    time_pdm_acc = time_pdm_acc + sto - sta
    bytes_transferred_acc = bytes_transferred_acc + nelms * 8_long
    nmsg_acc     = nmsg_acc + 1
#endif
  end subroutine array_accumulate_tile_combidx8



  subroutine array_accumulate_tile_combidx_nobuff(A,dimsA,arr,globtilenr,o,lock_set)
    implicit none
    real(realk),intent(in) :: A(*)
    type(array),intent(in) :: arr
    integer,intent(in) :: dimsA(arr%mode)
    integer,intent(in) :: globtilenr
    integer :: o(arr%mode)
    !> input the fortan array which should be transferred to the tile
    integer(kind=ls_mpik) :: dest
    logical, intent(in) :: lock_set
    integer :: order_type
    integer :: ro(arr%mode),rtd(arr%mode),nel
    integer :: acttdim(arr%mode),tmodeidx(arr%mode)
    integer :: idxintile(arr%mode),fels(arr%mode)
    integer :: glbmodeidx(arr%mode),pos1,i,k,nelms,ntimes,ccels
    real(realk) :: sta,sto,bs

#ifdef VAR_MPI

    sta=MPI_WTIME()

    order_type = -1
    bs=int(((8000.0*1000.0)/(8.0*2.0))**(1.0/float(arr%mode)))
    !bs=5
    order_type=0
    do i=1,arr%mode
      if(o(i)/=i)order_type=-1
    enddo


    do i=1,arr%mode
      !get the reverse order information
      ro(o(i))=i
    enddo

    dest=get_residence_of_tile(globtilenr,arr)

    if(.not.lock_set)call lsmpi_win_lock(dest,arr%wi(globtilenr),'s')

    if(arr%mode==4)then
      if(o(1)==1.and.o(2)==2.and.o(3)==3.and.o(4)==4)order_type = 0
      if(o(1)==3.and.o(2)==4.and.o(3)==1.and.o(4)==2)order_type = 1
      if(o(1)==4.and.o(2)==1.and.o(3)==2.and.o(4)==3)order_type = 2
      if(o(1)==2.and.o(2)==3.and.o(3)==4.and.o(4)==1)order_type = 3
      if(o(1)==1.and.o(2)==2.and.o(3)==4.and.o(4)==3)order_type = 4
      if(o(1)==1.and.o(2)==4.and.o(3)==2.and.o(4)==3)order_type = 5
      if(o(1)==1.and.o(2)==3.and.o(3)==4.and.o(4)==2)order_type = 6
      if(o(1)==3.and.o(2)==1.and.o(3)==2.and.o(4)==4)order_type = 7
      if(o(1)==2.and.o(2)==3.and.o(3)==1.and.o(4)==4)order_type = 8
      if(o(1)==2.and.o(2)==1.and.o(3)==3.and.o(4)==4)order_type = 9
      if(o(1)==4.and.o(2)==3.and.o(3)==1.and.o(4)==2)order_type = 10
      if(o(1)==4.and.o(2)==2.and.o(3)==3.and.o(4)==1)order_type = 11
      if(o(1)==3.and.o(2)==4.and.o(3)==2.and.o(4)==1)order_type = 12
      if(o(1)==2.and.o(2)==4.and.o(3)==1.and.o(4)==3)order_type = 13
      if(o(1)==3.and.o(2)==2.and.o(3)==1.and.o(4)==4)order_type = 14
      if(o(1)==1.and.o(2)==3.and.o(3)==2.and.o(4)==4)order_type = 15
      if(o(1)==4.and.o(2)==1.and.o(3)==3.and.o(4)==2)order_type = 16
      if(o(1)==2.and.o(2)==1.and.o(3)==4.and.o(4)==3)order_type = 17
      if(o(1)==4.and.o(2)==3.and.o(3)==2.and.o(4)==1)order_type = 18
      if(o(1)==2.and.o(2)==4.and.o(3)==3.and.o(4)==1)order_type = 19
      if(o(1)==1.and.o(2)==4.and.o(3)==3.and.o(4)==2)order_type = 20
      if(o(1)==3.and.o(2)==1.and.o(3)==4.and.o(4)==2)order_type = 21
      if(o(1)==3.and.o(2)==2.and.o(3)==4.and.o(4)==1)order_type = 22
      if(o(1)==4.and.o(2)==2.and.o(3)==1.and.o(4)==3)order_type = 23
    endif

    call get_midx(globtilenr,tmodeidx,arr%ntpm,arr%mode)
    ntimes=1
    nel = 1
    do i=1,arr%mode
      fels(o(i)) = (tmodeidx(i)-1) * arr%tdim(i) + 1
      if(tmodeidx(i)*arr%tdim(i)>arr%dims(i))then
        acttdim(i)=mod(arr%dims(i),arr%tdim(i))
      else
        acttdim(i)=arr%tdim(i)
      endif
      if(i>1)ntimes=ntimes*acttdim(i)
      nel = acttdim(i) * nel
      rtd(o(i))     = acttdim(i)
    enddo

    select case(order_type)
      case(0)
        ccels=acttdim(1)
        !loop over the remaining not-consecutive dimensions
        do i=1,ntimes
          !get the mode-index in the remaining dimensions
          call get_midx(i,idxintile(2:arr%mode),acttdim(2:arr%mode),arr%mode-1)
          !get the position of the first element in the consecutive stretch
          idxintile(1)=1
          do k=1,arr%mode
            glbmodeidx(k)=idxintile(k) +(tmodeidx(k)-1) * arr%tdim(k)
          enddo
          pos1=get_cidx(glbmodeidx,arr%dims,arr%mode)
          !  call dcopy(ccels,fort(pos1),1,tileout(1+(i-1)*ccels),1)
          !  call daxpy(ccels,pre1,fort(pos1),1,tileout(1+(i-1)*ccels),1)
          call lsmpi_acc(A(pos1:pos1+ccels-1),ccels,1+(i-1)*ccels,dest,arr%wi(globtilenr))
        enddo
      !case(1)
      !  call manual_3412_reordering_f2t(bs,rtd,dimsA,fels,pre1,fort,pre2,tileout)
      !case(2)
      !  call manual_4123_reordering_f2t(bs,rtd,dimsA,fels,pre1,fort,pre2,tileout)
      !case(3)
      !  call manual_2341_reordering_f2t(bs,rtd,dimsA,fels,pre1,fort,pre2,tileout)
      !case(4)
      !  call manual_1243_reordering_f2t(bs,rtd,dimsA,fels,pre1,fort,pre2,tileout)
      !case(5)
      !  call manual_1423_reordering_f2t(bs,rtd,dimsA,fels,pre1,fort,pre2,tileout)
      !case(6)
      !  call manual_1342_reordering_f2t(bs,rtd,dimsA,fels,pre1,fort,pre2,tileout)
      !case(7)
      !  call manual_3124_reordering_f2t(bs,rtd,dimsA,fels,pre1,fort,pre2,tileout)
      !case(8)
      !  call manual_2314_reordering_f2t(bs,rtd,dimsA,fels,pre1,fort,pre2,tileout)
      !case(9)
      !  call manual_2134_reordering_f2t(bs,rtd,dimsA,fels,pre1,fort,pre2,tileout)
      !case(10)
      !  call manual_4312_reordering_f2t(bs,rtd,dimsA,fels,pre1,fort,pre2,tileout)
      !case(11)
      !  call manual_4231_reordering_f2t(bs,rtd,dimsA,fels,pre1,fort,pre2,tileout)
      !case(12)
      !  call manual_3421_reordering_f2t(bs,rtd,dimsA,fels,pre1,fort,pre2,tileout)
      !case(13)
      !  call manual_2413_reordering_f2t(bs,rtd,dimsA,fels,pre1,fort,pre2,tileout)
      !case(14)
      !  call manual_3214_reordering_f2t(bs,rtd,dimsA,fels,pre1,fort,pre2,tileout)
      !case(15)
      !  call manual_1324_reordering_f2t(bs,rtd,dimsA,fels,pre1,fort,pre2,tileout)
      !case(16)
      !  call manual_4132_reordering_f2t(bs,rtd,dimsA,fels,pre1,fort,pre2,tileout)
      !case(17)
      !  call manual_2143_reordering_f2t(bs,rtd,dimsA,fels,pre1,fort,pre2,tileout)
      !case(18)
      !  call manual_4321_reordering_f2t(bs,rtd,dimsA,fels,pre1,fort,pre2,tileout)
      !case(19)
      !  call manual_2431_reordering_f2t(bs,rtd,dimsA,fels,pre1,fort,pre2,tileout)
      !case(20)
      !  call manual_1432_reordering_f2t(bs,rtd,dimsA,fels,pre1,fort,pre2,tileout)
      !case(21)
      !  call manual_3142_reordering_f2t(bs,rtd,dimsA,fels,pre1,fort,pre2,tileout)
      !case(22)
      !  call manual_3241_reordering_f2t(bs,rtd,dimsA,fels,pre1,fort,pre2,tileout)
      !case(23)
      !  call manual_4213_reordering_f2t(bs,rtd,dimsA,fels,pre1,fort,pre2,tileout)
      case default
        print *,"expensive default tile_from_fort",o
        !print *,"order  :",o
        !print *,"rorder :",ro
        !print *,"atd    :",acttdim
        !count elements in the current tile for loop over elements
        !identify their original position and put them in tile
        nelms=1
        do i=1,arr%mode
          nelms = nelms * acttdim(i)
        enddo
        do i = 1,nelms
          !get mode index of element in tile
          call get_midx(i,idxintile,acttdim,arr%mode)
          !get global index of element, example order = 2 3 1 4 of new array with
          !respect to old --> element 54 3 27 8 of old goes to 3 27 54 8 of new -
          ! old with respect to new 3 1 2 4 
          do k=1,arr%mode
            glbmodeidx(o(k))=idxintile(k) + (tmodeidx(k)-1)*arr%tdim(k)
          enddo
          pos1=get_cidx(glbmodeidx,dimsA,arr%mode)
          !DECinfo%ccModel>2tileout(i)=pre2*tileout(i)+pre1*A(pos1)
          call lsmpi_acc(A(pos1:pos1),1,i,dest,arr%wi(globtilenr))
        enddo
    end select
    if(.not.lock_set) call lsmpi_win_unlock(dest,arr%wi(globtilenr))
    sto = MPI_WTIME()
    time_pdm_acc = time_pdm_acc + sto - sta
    bytes_transferred_acc = bytes_transferred_acc + nel * 8_long
    nmsg_acc = nmsg_acc + 1
#endif
  end subroutine array_accumulate_tile_combidx_nobuff


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!                   PUT TILES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine array_puttile_modeidx(arr,modidx,fort,nelms,lock_set,flush_it)
    implicit none
    type(array),intent(in) ::arr
    integer,intent(in) :: modidx(arr%mode),nelms
    real(realk),intent(inout) :: fort(*)
    logical, optional, intent(in) :: lock_set,flush_it
    logical :: ls
    integer :: cidx
    ls = .false.
    if(present(lock_set))ls=lock_set
    cidx=get_cidx(modidx,arr%ntpm,arr%mode)
    call array_put_tile(arr,cidx,fort,nelms,lock_set=lock_set,flush_it=flush_it)
  end subroutine array_puttile_modeidx

  subroutine array_putt8(arr,globtilenr,fort,nelms,lock_set,flush_it)
    implicit none
    type(array),intent(in) :: arr
    integer,intent(in) :: globtilenr
    integer(kind=8),intent(in) :: nelms
    real(realk),intent(inout) :: fort(*)
    logical, optional, intent(in) :: lock_set,flush_it
    call array_puttile_combidx8(arr,globtilenr,fort,nelms,lock_set=lock_set,flush_it=flush_it)
  end subroutine array_putt8
  subroutine array_puttile_combidx8(arr,globtilenr,fort,nelms,lock_set,flush_it)
    implicit none
    type(array),intent(in) :: arr
    integer,intent(in) :: globtilenr
    integer(kind=8),intent(in) :: nelms
    real(realk),intent(inout) :: fort(*)
    logical, optional, intent(in) :: lock_set,flush_it
    logical :: ls
    integer(kind=ls_mpik) :: dest
    real(realk) :: sta,sto
#ifdef VAR_MPI
    integer :: maxsze

    maxsze = MAX_SIZE_ONE_SIDED
    ls = .false.
    if(present(lock_set))ls=lock_set

    dest = get_residence_of_tile(globtilenr,arr)

    sta  = MPI_WTIME()


    if(.not.ls)call lsmpi_win_lock(dest,arr%wi(globtilenr),'s')
    call lsmpi_put(fort,nelms,1,dest,arr%wi(globtilenr),maxsze,flush_it=flush_it)
    if(.not.ls)call lsmpi_win_unlock(dest,arr%wi(globtilenr))

    sto = MPI_WTIME()

    time_pdm_put          = time_pdm_put + sto - sta
    bytes_transferred_put = bytes_transferred_put + nelms * 8_long
    nmsg_put              = nmsg_put + 1
#endif
  end subroutine array_puttile_combidx8
  subroutine array_putt4(arr,globtilenr,fort,nelms,lock_set,flush_it)
    implicit none
    type(array),intent(in) :: arr
    integer,intent(in) :: globtilenr
    integer(kind=4),intent(in) :: nelms
    real(realk),intent(inout) :: fort(*)
    logical, optional, intent(in) :: lock_set,flush_it
    call array_puttile_combidx4(arr,globtilenr,fort,nelms,lock_set=lock_set,flush_it=flush_it)
  end subroutine array_putt4
  subroutine array_puttile_combidx4(arr,globtilenr,fort,nelms,lock_set,flush_it)
    implicit none
    type(array),intent(in) :: arr
    integer,intent(in) :: globtilenr
    integer(kind=4),intent(in) :: nelms
    real(realk),intent(inout) :: fort(*)
    logical, optional, intent(in) :: lock_set,flush_it
    logical :: ls
    integer(kind=ls_mpik) :: dest
    real(realk) :: sta,sto
#ifdef VAR_MPI
    integer :: maxsze
    maxsze = MAX_SIZE_ONE_SIDED

    ls = .false.
    if(present(lock_set))ls=lock_set

    dest = get_residence_of_tile(globtilenr,arr)

    sta  = MPI_WTIME()

    if(.not.ls)call lsmpi_win_lock(dest,arr%wi(globtilenr),'s')
    call lsmpi_put(fort,nelms,1,dest,arr%wi(globtilenr),maxsze,flush_it = flush_it)
    if(.not.ls)call lsmpi_win_unlock(dest,arr%wi(globtilenr))

    sto = MPI_WTIME()

    time_pdm_put          = time_pdm_put + sto - sta
    bytes_transferred_put = bytes_transferred_put + nelms * 8_long
    nmsg_put              = nmsg_put + 1
#endif
  end subroutine array_puttile_combidx4



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!                   GET TILES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !interface to the array_gettile_combidx
  subroutine array_gettile_modeidx(arr,modidx,fort,nelms,lock_set,flush_it)
    implicit none
    type(array),intent(in) ::arr
    integer,intent(in) :: modidx(arr%mode),nelms
    real(realk),intent(inout) :: fort(*)
    logical, optional, intent(in) :: lock_set,flush_it
    logical :: ls
    integer :: cidx
    ls = .false.
    if(present(lock_set))ls=lock_set
    cidx=get_cidx(modidx,arr%ntpm,arr%mode)
    call array_get_tile(arr,cidx,fort,nelms,lock_set=ls,flush_it=flush_it)
  end subroutine array_gettile_modeidx
  subroutine array_gett8(arr,globtilenr,fort,nelms,lock_set,flush_it)
    implicit none
    type(array),intent(in) :: arr
    integer,intent(in) :: globtilenr
    integer(kind=8),intent(in) :: nelms
    real(realk),intent(inout) :: fort(*)
    logical, optional, intent(in) :: lock_set,flush_it
    call array_gettile_combidx8(arr,globtilenr,fort,nelms,lock_set=lock_set,flush_it=flush_it)
  end subroutine array_gett8
  subroutine array_gettile_combidx8(arr,globtilenr,fort,nelms,lock_set,flush_it)
    implicit none
    type(array),intent(in) :: arr
    integer,intent(in) :: globtilenr
    integer(kind=8),intent(in) :: nelms
    real(realk),intent(inout) :: fort(*)
    logical, optional, intent(in) :: lock_set,flush_it
    integer(kind=ls_mpik) :: source
    real(realk) :: sta,sto
    logical :: ls
#ifdef VAR_MPI
    integer :: maxsze
    maxsze = MAX_SIZE_ONE_SIDED

    ls = .false.
    if(present(lock_set))ls=lock_set

    source = get_residence_of_tile(globtilenr,arr)

    sta    = MPI_WTIME()

    if(.not.ls)call lsmpi_win_lock(source,arr%wi(globtilenr),'s')
    call lsmpi_get(fort,nelms,1,source,arr%wi(globtilenr),maxsze,flush_it=flush_it)
    if(.not.ls)call lsmpi_win_unlock(source,arr%wi(globtilenr))

    sto = MPI_WTIME()

    time_pdm_get          = time_pdm_get + sto - sta
    bytes_transferred_get = bytes_transferred_get + nelms * 8_long
    nmsg_get              = nmsg_get + 1
#endif
  end subroutine array_gettile_combidx8
  subroutine array_gett4(arr,globtilenr,fort,nelms,lock_set,flush_it)
    implicit none
    type(array),intent(in) :: arr
    integer,intent(in) :: globtilenr
    integer(kind=4),intent(in) :: nelms
    real(realk),intent(inout) :: fort(*)
    logical, optional, intent(in) :: lock_set,flush_it
    call array_gettile_combidx4(arr,globtilenr,fort,nelms,lock_set=lock_set,flush_it=flush_it)
  end subroutine array_gett4
  subroutine array_gettile_combidx4(arr,globtilenr,fort,nelms,lock_set,flush_it)
    implicit none
    type(array),intent(in) :: arr
    integer,intent(in) :: globtilenr
    integer(kind=4),intent(in) :: nelms
    real(realk),intent(inout) :: fort(*)
    logical, optional, intent(in) :: lock_set,flush_it
    integer(kind=ls_mpik) :: source
    real(realk) :: sta,sto
    logical :: ls
#ifdef VAR_MPI
    integer :: maxsze
    maxsze = MAX_SIZE_ONE_SIDED

    ls = .false.
    if(present(lock_set))ls=lock_set

    source = get_residence_of_tile(globtilenr,arr)
    sta    = MPI_WTIME()

    if(.not.ls)call lsmpi_win_lock(source,arr%wi(globtilenr),'s')
    call lsmpi_get(fort,nelms,1,source,arr%wi(globtilenr),maxsze,flush_it=flush_it)
    if(.not.ls)call lsmpi_win_unlock(source,arr%wi(globtilenr))

    sto = MPI_WTIME()

    time_pdm_get          = time_pdm_get + sto - sta
    bytes_transferred_get = bytes_transferred_get + nelms * 8_long
    nmsg_get              = nmsg_get + 1
#endif
  end subroutine array_gettile_combidx4

  subroutine get_int_dist_info(o2v2,firstintel,nintel,remoterank)
    implicit none
    integer(kind=long), intent(in) :: o2v2
    integer, intent(out) :: firstintel,nintel
    integer(kind=ls_mpik), intent(in), optional :: remoterank
    integer(kind=ls_mpik) :: nnod, me
    nnod = 1
    me   = 0
#ifdef VAR_MPI
    nnod = infpar%lg_nodtot
    if(.not.present(remoterank))then
      me = infpar%lg_mynum
    else
      me = remoterank
    endif
#endif
    nintel = o2v2/nnod
    firstintel = me*nintel + 1
    if(me<int(mod(o2v2,int(nnod,kind=long)),kind=ls_mpik))then
      nintel = nintel + 1
      firstintel = firstintel + int(me) 
    else if(me>=int(mod(o2v2,int(nnod,kind=long)),kind=ls_mpik))then
      firstintel = firstintel + int(mod(o2v2,int(nnod,kind=long))) 
    endif
  end subroutine get_int_dist_info

  subroutine dist_int_contributions(g,o2v2,win,lock_outside)
    implicit none
    integer(kind=long),intent(in) :: o2v2
    real(realk),intent(in) :: g(o2v2)
    logical :: lock_outside
    integer(kind=ls_mpik),intent(in) :: win
    integer(kind=ls_mpik) :: nnod,node,me
    integer :: fe,ne,msg_len_mpi
    real(realk) :: sta,sto
    fe=1
    ne=0
    nnod = 1

#ifdef VAR_MPI
#ifdef VAR_LSDEBUG
    msg_len_mpi=24
#else
    msg_len_mpi=SPLIT_MPI_MSG
#endif
    nnod = infpar%lg_nodtot
    me   = infpar%lg_mynum
    do node=0,nnod-1
      call get_int_dist_info(o2v2,fe,ne,node)
      sta=MPI_WTIME()
      !print *,infpar%lg_mynum,"distributing",fe,fe+ne-1,ne,o2v2,node
      if(.not.lock_outside)call lsmpi_win_lock(node,win,'s')
      call lsmpi_acc(g(fe:fe+ne-1),ne,1,node,win,msg_len_mpi,.true.)
      if(.not.lock_outside)call lsmpi_win_unlock(node,win)
      sto = MPI_WTIME()
      time_pdm_acc = time_pdm_acc + sto - sta
      bytes_transferred_acc = bytes_transferred_acc + ne * 8_long
      nmsg_acc = nmsg_acc + 1
    enddo
#endif
  end subroutine dist_int_contributions

  subroutine collect_int_contributions(g,o2v2,win)
    implicit none
    integer(kind=long),intent(in) :: o2v2
    real(realk),intent(in) :: g(o2v2)
    integer(kind=ls_mpik),intent(in) :: win
    integer(kind=ls_mpik) :: nnod,node,me
    integer :: fe,ne,msg_len_mpi
    real(realk) :: sta,sto
    fe=1
    ne=0
    nnod = 1
#ifdef VAR_LSDEBUG
    msg_len_mpi=24
#else
    msg_len_mpi=170000000
#endif
#ifdef VAR_MPI
    nnod = infpar%lg_nodtot
    me   = infpar%lg_mynum
    do node=0,nnod-1
      !print *,infpar%lg_mynum,"collecting",fe,fe+ne-1,ne,o2v2,node
      call get_int_dist_info(o2v2,fe,ne,node)
      sta=MPI_WTIME()
      call lsmpi_win_lock(node,win,'s')
      call lsmpi_get(g(fe:fe+ne-1),ne,1,node,win,msg_len_mpi)
      call lsmpi_win_unlock(node,win)
      sto = MPI_WTIME()
      time_pdm_get = time_pdm_get + sto - sta
      bytes_transferred_get = bytes_transferred_get + ne * 8_long
      nmsg_get = nmsg_get + 1
    enddo
#endif
  end subroutine collect_int_contributions

  !fenced routine does not need locks, but rewuires fences around call
  subroutine collect_int_contributions_f(g,o2v2,win)
    implicit none
    integer(kind=long),intent(in) :: o2v2
    real(realk),intent(in) :: g(o2v2)
    integer(kind=ls_mpik),intent(in) :: win
    integer(kind=ls_mpik) :: nnod,node,me
    integer :: fe,ne,msg_len_mpi
    real(realk) :: sta,sto
    fe=1
    ne=0
    nnod = 1
    !msg_len_mpi=17
#ifdef VAR_LSDEBUG
    msg_len_mpi=24
#else
    msg_len_mpi=170000000
#endif
#ifdef VAR_MPI
    nnod = infpar%lg_nodtot
    me   = infpar%lg_mynum
    do node=0,nnod-1
      !print *,infpar%lg_mynum,"collecting f",fe,fe+ne-1,ne,o2v2,node
      sta=MPI_WTIME()
      call get_int_dist_info(o2v2,fe,ne,node)
      call lsmpi_get(g(fe:fe+ne-1),ne,1,node,win,msg_len_mpi)
      sto = MPI_WTIME()
      time_pdm_get = time_pdm_get + sto - sta
      bytes_transferred_get = bytes_transferred_get + ne * 8_long
      nmsg_get = nmsg_get + 1
    enddo
#endif
  end subroutine collect_int_contributions_f



  subroutine array_scale_td(arr,sc)
    implicit none
    type(array) :: arr
    real(realk) :: sc
#ifdef VAR_MPI
    integer     :: i

    if(arr%access_type==MASTER_ACCESS)then
      call PDM_ARRAY_SYNC(infpar%lg_comm,JOB_ARRAY_SCALE,arr)
      call ls_mpibcast(sc,infpar%master,infpar%lg_comm)
    endif

    do i=1,arr%nlti
      call dscal(int(arr%ti(i)%e),sc,arr%ti(i)%t,1)
    enddo
#endif
  end subroutine array_scale_td

  
  subroutine memory_allocate_array_dense_pc(arr)
      implicit none
      type(array), intent(inout) :: arr
      logical :: parent
#ifdef VAR_MPI
      parent = (infpar%parent_comm == MPI_COMM_NULL)
      if(lspdm_use_comm_proc.and.parent.and.arr%access_type==MASTER_ACCESS)then
        call pdm_array_sync(infpar%pc_comm,JOB_PC_ALLOC_DENSE,arr,loc_addr=.true.)
      endif
#endif
      call memory_allocate_array_dense(arr)
  end subroutine memory_allocate_array_dense_pc

  subroutine memory_deallocate_array_dense_pc(arr)
      implicit none
      type(array), intent(inout) :: arr
      logical :: parent
#ifdef VAR_MPI
      parent = (infpar%parent_comm == MPI_COMM_NULL)
      if(lspdm_use_comm_proc.and.parent.and.arr%access_type==MASTER_ACCESS)then
        call pdm_array_sync(infpar%pc_comm,JOB_PC_DEALLOC_DENSE,arr,loc_addr=.true.)
      endif
#endif
      call arr_deallocate_dense(arr)
  end subroutine memory_deallocate_array_dense_pc

  subroutine lsmpi_put_realkV_w8(buf,nelms,pos,dest,win)
    implicit none
    real(realk),intent(in) :: buf(*)
    integer, intent(in) :: pos
    integer(kind=8) :: nelms
    integer(kind=ls_mpik),intent(in) :: dest
    integer(kind=ls_mpik),intent(in) :: win
#ifdef VAR_MPI
    call lsmpi_put_realkV_wrapper8(buf,nelms,pos,dest,win)
#endif
  end subroutine lsmpi_put_realkV_w8
  subroutine lsmpi_get_realkV_w8(buf,nelms,pos,dest,win)
    implicit none
    real(realk),intent(in) :: buf(*)
    integer, intent(in) :: pos
    integer(kind=8) :: nelms
    integer(kind=ls_mpik),intent(in) :: dest
    integer(kind=ls_mpik),intent(in) :: win
#ifdef VAR_MPI
    call lsmpi_get_realkV_wrapper8(buf,nelms,pos,dest,win)
#endif
  end subroutine lsmpi_get_realkV_w8
  subroutine lsmpi_acc_realkV_w8(buf,nelms,pos,dest,win)
    implicit none
    real(realk),intent(in) :: buf(*)
    integer, intent(in) :: pos
    integer(kind=8) :: nelms
    integer(kind=ls_mpik),intent(in) :: dest
    integer(kind=ls_mpik),intent(in) :: win
#ifdef VAR_MPI
    call lsmpi_acc_realkV_wrapper8(buf,nelms,pos,dest,win)
#endif
  end subroutine lsmpi_acc_realkV_w8


end module lspdm_tensor_operations_module




